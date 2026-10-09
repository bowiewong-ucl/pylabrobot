"""Agilent BioTek Cytation 7 plate-reader integration."""

from __future__ import annotations

import asyncio
import logging
import math
import re
import time
from dataclasses import dataclass
from datetime import datetime
from typing import AsyncIterator, Dict, List, Literal, Optional, Sequence, Tuple, Union

from pylabrobot.agilent.biotek.cytation.base import _CytationBase
from pylabrobot.agilent.biotek.plate_reader_base import AbsorbanceResult, FluorescenceResult
from pylabrobot.resources import Plate, Well

logger = logging.getLogger(__name__)

FluorescenceOptics = Literal["top", "bottom"]
LampEnergy = Literal["high", "low"]
ReadSpeed = Literal["normal", "sweep"]

_OPTICS_CODES: Dict[str, str] = {"bottom": "4", "top": "5"}
_LAMP_ENERGY_CODES: Dict[str, str] = {"low": "0", "high": "2"}
_READ_SPEED_CODES: Dict[str, str] = {"normal": "0", "sweep": "1"}

# Read timing block after the well range: 100 ms settling delay for "normal" reads, none for
# "sweep" reads.
_TIMING_NORMAL = "000120010000"
_TIMING_SWEEP = "000120000000"

# Schedule block for a single endpoint pass: one read, no interval, no shaking.
_SCHEDULE_ENDPOINT = "1100100000"

# Shake mode code in the kinetic schedule block: 03 is double orbital.
_SHAKE_DOUBLE_ORBITAL = "03"

# One result block per read step: `\r<hhh:mm:ss.s>,<temperature x 10>,<row,col,value...>\r\n`
# followed by the temperature again, `\x1a`, and a 3-digit code.
_RESULT_BLOCK = re.compile(
  rb"\r(\d{3}):(\d{2}):(\d{2}\.\d),(\d+),(.*?)\r\n\d+\x1a(?:\d{3})?", re.DOTALL
)

# A complete result block at the end of a kinetic stream buffer. A reader that held its results
# while disconnected sends them without the block's leading `\r`, so the block may also start
# the buffer.
_KINETIC_BLOCK = re.compile(
  rb"(?:\r|^)(\d{3}):(\d{2}):(\d{2}\.\d),(\d+),(.*?)\r\n\d+\x1a\d{3}$", re.DOTALL
)

# Upper bound on one result block: a full 96-well block is about 1.5 kB.
_MAX_BLOCK_BYTES = 64 * 1024

# Serial settings the Cytation 7 uses: 38461 baud, 8 data bits, 2 stop bits, RTS/CTS.
_BAUDRATE = 38_461
_SIO_RTS_CTS_HS = 0x1 << 8

Grid = List[List[Optional[float]]]


@dataclass(frozen=True)
class AbsorbanceStep:
  """Absorbance read step of a kinetic run: Normal speed, 100 ms delay, 8 measurements.

  Attributes:
    wavelengths: Wavelengths in nm; each is its own read step on the reader.
  """

  wavelengths: Sequence[int]


@dataclass(frozen=True)
class FluorescenceStep:
  """Fluorescence read step (one filter set) of a kinetic run, at Normal speed.

  Attributes:
    excitation_wavelength: Excitation center wavelength, in nm.
    emission_wavelength: Emission center wavelength, in nm.
    optics: Read from the "top" or the "bottom" of the plate.
    gain: Detector gain.
    excitation_bandwidth: Excitation bandwidth, in nm (0.1 nm resolution).
    emission_bandwidth: Emission bandwidth, in nm (0.1 nm resolution).
    lamp_energy: Xenon flash lamp energy.
    measurements_per_data_point: Measurements averaged per well.
  """

  excitation_wavelength: int
  emission_wavelength: int
  optics: FluorescenceOptics = "top"
  gain: int = 100
  excitation_bandwidth: float = 20
  emission_bandwidth: float = 20
  lamp_energy: LampEnergy = "high"
  measurements_per_data_point: int = 10


@dataclass(frozen=True)
class KineticRead:
  """One read of one step during a kinetic run.

  Attributes:
    cycle: 0-based kinetic cycle.
    elapsed: Time of the read since the run started, in seconds, as reported by the reader.
    step: The `AbsorbanceStep` or `FluorescenceStep` this read belongs to.
    result: Measured values, with the chamber temperature the reader reported for this read.
  """

  cycle: int
  elapsed: float
  step: Union[AbsorbanceStep, FluorescenceStep]
  result: Union[AbsorbanceResult, FluorescenceResult]


KineticSteps = Sequence[Union[AbsorbanceStep, FluorescenceStep]]
# One entry per read step on the reader: the user's step and, for absorbance, its wavelength.
ExpandedSteps = List[Tuple[Union[AbsorbanceStep, FluorescenceStep], Optional[int]]]


def _parse_records(records: bytes, plate: Plate) -> Grid:
  """Parse `row,col,value` records into a `[row][column]` grid.

  Asterisk values (overflow) become NaN; wells outside the read are `None`.

  Args:
    records: Comma-separated records of one result block.
    plate: Plate that was read.
  """
  fields = [f for f in records.replace(b"\r\n", b"").split(b",") if f]
  if len(fields) % 3 != 0:
    raise ValueError(f"Malformed result records: {records!r}")
  grid: Grid = [[None for _ in range(plate.num_items_x)] for _ in range(plate.num_items_y)]
  for i in range(0, len(fields), 3):
    row, column, raw = int(fields[i]), int(fields[i + 1]), fields[i + 2].decode()
    try:
      value = float("nan") if "*" in raw else float(raw)
    except ValueError:
      logger.warning("Unreadable value %r for row %d, column %d", raw, row, column)
      value = float("nan")
    grid[row - 1][column - 1] = value
  return grid


def _parse_result_blocks(body: bytes, plate: Plate) -> List[Grid]:
  """Parse every result block in an acquisition body into `[row][column]` grids.

  Args:
    body: Acquisition body, through the terminating `\\x03`.
    plate: Plate that was read.
  """
  return [_parse_records(match.group(5), plate) for match in _RESULT_BLOCK.finditer(body)]


def _merge_into(target: Grid, source: Grid) -> None:
  """Copy every measured (non-`None`) value of `source` into `target`."""
  for r, row in enumerate(source):
    for c, value in enumerate(row):
      if value is not None:
        target[r][c] = value


def _absorbance_steps(wavelengths: Sequence[int]) -> List[str]:
  """Encode one absorbance read step per wavelength."""
  return [f"0600008{wl:04}" for wl in wavelengths]


def _fluorescence_step(step: FluorescenceStep, read_speed: ReadSpeed) -> str:
  """Encode a fluorescence read step (one filter set) in the Cytation 7 layout.

  Args:
    step: Fluorescence settings.
    read_speed: "sweep" encodes a single measurement per well.
  """
  measurements = 1 if read_speed == "sweep" else step.measurements_per_data_point
  return (
    f"3{_OPTICS_CODES[step.optics]}{_READ_SPEED_CODES[read_speed]}{measurements:04}0"
    f"{round(step.excitation_bandwidth * 10):03}{round(step.emission_bandwidth * 10):03}"
    f"{step.excitation_wavelength:04}000{step.emission_wavelength:04}{'0' * 13}"
    f"{_LAMP_ENERGY_CODES[step.lamp_energy]}{step.gain:03}0"
  )


class Cytation7(_CytationBase):
  """Agilent BioTek Cytation 7 plate reader.

  Absorbance, fluorescence, temperature, and kinetic runs use Cytation 7 command layouts
  (firmware 1.02); the remaining operations use the shared BioTek commands.
  Microscopy is not initialized. The incubator heats but does not actively cool. Command
  checksums use the firmware rule, which a Cytation 7 enforces.
  """

  _model_name = "Agilent BioTek Cytation 7"

  def _checksum(self, payload: str, hundreds_digit: str, offset: int) -> str:
    """Return the firmware checksum: the byte sum after the 4-digit length, modulo 256.

    The Cytation 7 firmware (1.02) checks this checksum. It agrees with the shared reader's
    per-command checksum except where the sum modulo 256 falls outside that
    command's fixed hundreds digit, in which case the reader rejects the shared checksum.

    Args:
      payload: Command text from the 4-digit length through the last parameter digit.
      hundreds_digit: Unused; the firmware checksum does not depend on the command type.
      offset: Unused; the firmware checksum does not depend on the command type.
    """
    return f"{sum(payload[4:].encode()) % 256:03d}"

  def _optics_height(self, plate: Plate, focal_height: float, mode_digit: str) -> str:
    """Return the `t` parameter: the mode digit, then plate height plus focal height, in µm.

    For example, an 8 mm top fluorescence read on a 15.00 mm plate is `623000`. A lid is not
    included.

    Args:
      plate: Plate being read.
      focal_height: Focal height above the plate, in mm.
      mode_digit: Read mode digit (`6` for fluorescence, `3` for luminescence).
    """
    return f"{mode_digit}{round(1000 * (plate.get_size_z() + focal_height))}"

  def _read_command(
    self,
    rectangle: Sequence[int],
    timing: str,
    steps: Sequence[str],
    kinetic: bool = False,
    schedule: str = _SCHEDULE_ENDPOINT,
  ) -> str:
    """Return a complete `D` parameter: length, wells, schedule, steps, and checksum.

    Args:
      rectangle: 0-based min row, min column, max row, max column of the wells to read.
      timing: 12-digit read timing block (`_TIMING_NORMAL` or `_TIMING_SWEEP`).
      steps: Encoded read steps, one per absorbance wavelength or fluorescence filter set.
      kinetic: Whether the schedule repeats the steps as a kinetic run.
      schedule: Schedule block; a single endpoint pass unless given.
    """
    if not 1 <= len(steps) <= 9:
      raise ValueError(f"{self._model_name}: a read takes 1 to 9 steps, got {len(steps)}")
    min_row, min_col, max_row, max_col = rectangle
    body = (
      f"{'1' if kinetic else '0'}1"
      f"{min_row + 1:02}{min_col + 1:02}{max_row + 1:02}{max_col + 1:02}"
      f"{timing}{schedule}{len(steps)}{''.join(steps)}"
    )
    payload = f"{len(body) + 3:04}{body}"
    return payload + self._checksum(payload, "", 0) + "\x03"

  def _validate_fluorescence(
    self, step: FluorescenceStep, focal_height: Optional[float] = None
  ) -> None:
    """Raise `ValueError` for fluorescence settings outside the reader's ranges.

    Args:
      step: Fluorescence settings to check.
      focal_height: Read height for top reads, in mm, if one is used.
    """
    if focal_height is not None:
      min_fh, max_fh = self.focal_height_range
      if not min_fh <= focal_height <= max_fh:
        raise ValueError(f"{self._model_name}: focal height must be within {min_fh}-{max_fh}")
    min_ex, max_ex = self.excitation_range
    if not min_ex <= step.excitation_wavelength <= max_ex:
      raise ValueError(f"{self._model_name}: excitation wavelength must be {min_ex}-{max_ex}")
    min_em, max_em = self.emission_range
    if not min_em <= step.emission_wavelength <= max_em:
      raise ValueError(f"{self._model_name}: emission wavelength must be {min_em}-{max_em}")
    if not 0 <= step.gain <= 999:
      raise ValueError(f"{self._model_name}: gain must be within 0-999, got {step.gain}")
    for name, bandwidth in (
      ("excitation", step.excitation_bandwidth),
      ("emission", step.emission_bandwidth),
    ):
      if not 0 < round(bandwidth * 10) <= 999:
        raise ValueError(f"{self._model_name}: {name} bandwidth must be within 0.1-99.9 nm")
    if not 1 <= step.measurements_per_data_point <= 9999:
      raise ValueError(f"{self._model_name}: measurements per data point must be 1-9999")

  def _validate_absorbance(self, wavelengths: Sequence[int]) -> None:
    """Raise `ValueError` for absorbance wavelengths outside the reader's range."""
    if not wavelengths:
      raise ValueError(f"{self._model_name}: give at least one wavelength")
    min_abs, max_abs = self.abs_wavelength_range
    for wl in wavelengths:
      if not min_abs <= wl <= max_abs:
        raise ValueError(f"{self._model_name}: wavelength must be within {min_abs}-{max_abs}")

  @property
  def supports_heating(self) -> bool:
    """The Cytation 7 incubator can heat the read chamber."""
    return True

  @property
  def supports_cooling(self) -> bool:
    """The Cytation 7 has no active cooling."""
    return False

  async def set_temperature(self, temperature: float, gradient: float = 0) -> Optional[bytes]:
    """Set the incubator, optionally with a temperature gradient, without waiting.

    The parameter is the setpoint, then the gradient, each in 0.1 °C: `32505` is 32.5 °C with
    a 0.5 °C gradient.

    Args:
      temperature: Setpoint, in °C (0.1 °C resolution).
      gradient: Temperature gradient, in °C (0.1 °C resolution, at most 9.9 °C).
    """
    if temperature < 0:
      raise ValueError(f"{self._model_name}: temperature must not be negative, got {temperature}")
    if not 0 <= round(gradient * 10) <= 99:
      raise ValueError(f"{self._model_name}: gradient must be within 0-9.9 °C, got {gradient}")
    tmin, tmax = self.temperature_range
    if (tmin is not None and temperature < tmin) or (tmax is not None and temperature > tmax):
      raise ValueError(
        f"{self._model_name}: requested temperature {temperature} °C is outside {tmin}-{tmax} °C"
      )
    current_temperature = await self.request_current_temperature()
    if temperature < current_temperature and not self.supports_cooling:
      raise ValueError(f"{self.__class__.__name__}: Cooling is not supported.")
    return await self.send_command("g", f"{round(temperature * 10):03}{round(gradient * 10):02}")

  async def read_absorbance(
    self,
    plate: Plate,
    wells: List[Well],
    wavelength: Union[int, Sequence[int]],
  ) -> List[AbsorbanceResult]:
    """Read absorbance at one or more wavelengths in a single pass.

    Each wavelength is a read step of the same command. Settings are Normal speed, a 100 ms
    delay, and 8 measurements per data point.

    Args:
      plate: Plate to read.
      wells: Wells to read.
      wavelength: Wavelength in nm, or a sequence of up to 9 wavelengths.

    Returns:
      One result per wavelength, in the order given.
    """
    wavelengths = [wavelength] if isinstance(wavelength, int) else list(wavelength)
    self._validate_absorbance(wavelengths)

    logger.info(
      "[BioTek %s] read_absorbance: plate=%s, wavelengths=%s, wells=%d",
      self.io.device_id,
      plate.name,
      wavelengths,
      len(wells),
    )
    await self.set_plate(plate)

    all_data: List[Grid] = [
      [[None for _ in range(plate.num_items_x)] for _ in range(plate.num_items_y)]
      for _ in wavelengths
    ]
    steps = _absorbance_steps(wavelengths)
    for rectangle in self._get_min_max_row_col_tuples(wells, plate):
      command = self._read_command(rectangle, _TIMING_NORMAL, steps)
      body = await self._acquire(command, timeout=60 * 3 * len(wavelengths))
      grids = _parse_result_blocks(body, plate)
      if len(grids) != len(wavelengths):
        raise ValueError(f"Expected {len(wavelengths)} result blocks, got {len(grids)}")
      for target, grid in zip(all_data, grids):
        _merge_into(target, grid)

    try:
      temperature: Optional[float] = await self.request_current_temperature()
    except TimeoutError:
      temperature = None

    return [
      AbsorbanceResult(wavelength=wl, data=data, temperature=temperature, timestamp=datetime.now())
      for wl, data in zip(wavelengths, all_data)
    ]

  async def read_fluorescence(
    self,
    plate: Plate,
    wells: List[Well],
    excitation_wavelength: int,
    emission_wavelength: int,
    focal_height: float,
    optics: FluorescenceOptics = "top",
    gain: int = 100,
    excitation_bandwidth: float = 20,
    emission_bandwidth: float = 20,
    lamp_energy: LampEnergy = "high",
    read_speed: ReadSpeed = "normal",
    measurements_per_data_point: Optional[int] = None,
  ) -> List[FluorescenceResult]:
    """Read fluorescence intensity with one filter set.

    Reads use the Xenon flash at standard dynamic range. The optics are positioned only for
    top reads, so `focal_height` is sent and checked only when `optics` is "top".

    Args:
      plate: Plate to read.
      wells: Wells to read.
      excitation_wavelength: Excitation center wavelength, in nm.
      emission_wavelength: Emission center wavelength, in nm.
      focal_height: Read height above the plate for top reads, in mm.
      optics: Read from the "top" or the "bottom" of the plate.
      gain: Detector gain.
      excitation_bandwidth: Excitation bandwidth, in nm (0.1 nm resolution).
      emission_bandwidth: Emission bandwidth, in nm (0.1 nm resolution).
      lamp_energy: Xenon flash lamp energy.
      read_speed: "normal" waits 100 ms per well and averages several measurements; "sweep"
        takes a single measurement per well without waiting.
      measurements_per_data_point: Measurements averaged per well. Defaults to 10 for
        "normal"; "sweep" always takes 1.
    """
    if read_speed == "sweep" and measurements_per_data_point not in (None, 1):
      raise ValueError(f"{self._model_name}: a sweep read takes 1 measurement per data point")
    step = FluorescenceStep(
      excitation_wavelength=excitation_wavelength,
      emission_wavelength=emission_wavelength,
      optics=optics,
      gain=gain,
      excitation_bandwidth=excitation_bandwidth,
      emission_bandwidth=emission_bandwidth,
      lamp_energy=lamp_energy,
      measurements_per_data_point=(
        10 if measurements_per_data_point is None else measurements_per_data_point
      ),
    )
    self._validate_fluorescence(step, focal_height if optics == "top" else None)
    timing = _TIMING_SWEEP if read_speed == "sweep" else _TIMING_NORMAL

    logger.info(
      "[BioTek %s] read_fluorescence: plate=%s, excitation=%dnm, emission=%dnm, optics=%s, "
      "gain=%d, wells=%d",
      self.io.device_id,
      plate.name,
      excitation_wavelength,
      emission_wavelength,
      optics,
      gain,
      len(wells),
    )
    await self.set_plate(plate)
    if optics == "top":
      await self.send_command("t", self._optics_height(plate, focal_height, "6") + "\x03")

    encoded = _fluorescence_step(step, read_speed)
    all_data: Grid = [[None for _ in range(plate.num_items_x)] for _ in range(plate.num_items_y)]
    for rectangle in self._get_min_max_row_col_tuples(wells, plate):
      command = self._read_command(rectangle, timing, [encoded])
      body = await self._acquire(command, timeout=60 * 3)
      grids = _parse_result_blocks(body, plate)
      if len(grids) != 1:
        raise ValueError(f"Expected 1 result block, got {len(grids)}")
      _merge_into(all_data, grids[0])

    try:
      temperature: Optional[float] = await self.request_current_temperature()
    except TimeoutError:
      temperature = None

    return [
      FluorescenceResult(
        excitation_wavelength=excitation_wavelength,
        emission_wavelength=emission_wavelength,
        data=all_data,
        temperature=temperature,
        timestamp=datetime.now(),
      )
    ]

  async def _read_byte(self, deadline: float) -> bytes:
    """Return the next byte from the reader, raising `TimeoutError` after `deadline`."""
    while True:
      byte = await self.io.read(1)
      if byte:
        return bytes(byte)
      if time.time() > deadline:
        raise TimeoutError(f"{self._model_name}: no kinetic data before the deadline")
      await asyncio.sleep(0.01)

  def _expand_steps(
    self, steps: KineticSteps, focal_height: Optional[float]
  ) -> Tuple[List[str], ExpandedSteps]:
    """Validate and encode kinetic steps, one reader step per wavelength or filter set.

    Args:
      steps: Read steps per cycle, in order.
      focal_height: Read height for top fluorescence steps, in mm; `None` skips its check.

    Returns:
      The encoded reader steps, and for each the user's step and absorbance wavelength.
    """
    encoded: List[str] = []
    expanded: ExpandedSteps = []
    for step in steps:
      if isinstance(step, AbsorbanceStep):
        self._validate_absorbance(step.wavelengths)
        encoded.extend(_absorbance_steps(step.wavelengths))
        expanded.extend((step, wl) for wl in step.wavelengths)
      else:
        self._validate_fluorescence(step, focal_height if step.optics == "top" else None)
        encoded.append(_fluorescence_step(step, "normal"))
        expanded.append((step, None))
    return encoded, expanded

  async def _stream_kinetic(
    self, plate: Plate, expanded: ExpandedSteps, interval: float, first_read: int
  ) -> AsyncIterator[KineticRead]:
    """Yield kinetic reads from the reader's result stream until the run ends.

    Args:
      plate: Plate being read.
      expanded: Reader steps per cycle, from `_expand_steps`.
      interval: Seconds between cycle starts.
      first_read: Index (over all cycles and steps) of the next read the stream delivers.
    """
    # A cycle shakes for up to a whole interval before reading.
    block_timeout = interval + 30 * 60
    buffer = b""
    count = first_read
    deadline = time.time() + block_timeout
    while True:
      buffer += await self._read_byte(deadline)
      if time.time() > deadline or len(buffer) > _MAX_BLOCK_BYTES:
        raise TimeoutError(f"{self._model_name}: no complete kinetic result before the deadline")
      if buffer.endswith(b"\x03"):
        status = buffer[max(buffer.rfind(b"\x1a"), buffer.rfind(b"\x10")) + 1 :]
        if status != b"0000\x03":
          logger.warning("%s kinetic run ended with %r", self._model_name, buffer[-8:])
        return
      match = _KINETIC_BLOCK.search(buffer)
      if match is None:
        continue
      hours, minutes, seconds, temperature, records = match.groups()
      buffer = b""
      deadline = time.time() + block_timeout
      elapsed = int(hours) * 3600 + int(minutes) * 60 + float(seconds)
      cycle = count // len(expanded)
      if int(elapsed // interval) != cycle:
        logger.warning(
          "%s read at %.1f s falls in cycle %d, expected cycle %d",
          self._model_name,
          elapsed,
          int(elapsed // interval),
          cycle,
        )
      step, wavelength = expanded[count % len(expanded)]
      data = _parse_records(records, plate)
      result: Union[AbsorbanceResult, FluorescenceResult]
      if isinstance(step, AbsorbanceStep):
        assert wavelength is not None
        result = AbsorbanceResult(
          wavelength=wavelength,
          data=data,
          temperature=int(temperature) / 10,
          timestamp=datetime.now(),
        )
      else:
        result = FluorescenceResult(
          excitation_wavelength=step.excitation_wavelength,
          emission_wavelength=step.emission_wavelength,
          data=data,
          temperature=int(temperature) / 10,
          timestamp=datetime.now(),
        )
      yield KineticRead(cycle=cycle, elapsed=elapsed, step=step, result=result)
      count += 1

  async def run_kinetic(
    self,
    plate: Plate,
    wells: List[Well],
    steps: KineticSteps,
    reads: int,
    interval: float,
    shake_amplitude: int = 2,
    shake_duration: Optional[float] = None,
    focal_height: float = 7,
  ) -> AsyncIterator[KineticRead]:
    """Run a kinetic protocol on the reader and yield each read as the reader reports it.

    The whole run (schedule, double-orbital shaking, and read steps) is downloaded to the
    reader, which times, shakes, and reads on its own. Each cycle shakes, then reads every
    step; cycles start `interval` seconds apart. Set the incubator with `set_temperature`
    beforehand. Settings are checked when iteration starts. Leaving the loop stops listening;
    the run continues on the reader until `stop_kinetic` aborts it.

    Args:
      plate: Plate to read.
      wells: Wells to read; they must form one rectangle.
      steps: Read steps per cycle, in order. Each absorbance wavelength and each fluorescence
        step counts toward the reader's limit of 9.
      reads: Number of cycles, 2-999.
      interval: Time between the starts of consecutive cycles, in whole seconds.
      shake_amplitude: Double-orbital shake amplitude, in mm.
      shake_duration: Shake time before each cycle's reads, in seconds (a multiple of 10, at
        most 990). `None` shakes continuously between reads.
      focal_height: Read height above the plate for top fluorescence steps, in mm.
    """
    if not 2 <= reads <= 999:
      raise ValueError(f"{self._model_name}: reads must be within 2-999, got {reads}")
    if not math.isfinite(interval) or interval != int(interval) or not 1 <= interval <= 99999:
      raise ValueError(f"{self._model_name}: interval must be whole seconds within 1-99999")
    if not 1 <= shake_amplitude <= 6:
      raise ValueError(f"{self._model_name}: shake amplitude must be within 1-6 mm")
    if shake_duration is None:
      shake_tens = 0
    elif shake_duration % 10 != 0 or not 10 <= shake_duration <= 990:
      raise ValueError(f"{self._model_name}: shake duration must be a multiple of 10 s, 10-990")
    else:
      shake_tens = int(shake_duration // 10)
    rectangles = self._get_min_max_row_col_tuples(wells, plate)
    if len(rectangles) != 1:
      raise ValueError(f"{self._model_name}: kinetic wells must form one rectangle")

    encoded, expanded = self._expand_steps(steps, focal_height)
    top = any(isinstance(s, FluorescenceStep) and s.optics == "top" for s in steps)
    schedule = (
      f"23{shake_tens:02}{_SHAKE_DOUBLE_ORBITAL}{shake_amplitude}11{reads:03}{int(interval):05}"
    )
    command = self._read_command(rectangles[0], _TIMING_NORMAL, encoded, True, schedule)

    logger.info(
      "[BioTek %s] run_kinetic: plate=%s, wells=%d, steps=%d, reads=%d, interval=%ds",
      self.io.device_id,
      plate.name,
      len(wells),
      len(encoded),
      reads,
      interval,
    )
    await self.set_plate(plate)
    if top:
      await self.send_command("t", self._optics_height(plate, focal_height, "6") + "\x03")
    reply = await self.send_command("D", command)
    logger.info("%s kinetic download reply: %r", self._model_name, reply)
    started = await self.send_command("O")
    assert started == b"\x060000\x03"

    async for read in self._stream_kinetic(plate, expanded, interval, 0):
      yield read

  async def _attach(self) -> None:
    """Open the port to a reader that may be running a kinetic protocol, without disturbing it.

    The reader holds its results while no host is connected and sends them once the port
    asserts RTS again. A USB reset, a buffer purge, or a command would discard or garble
    them, and the reader answers no commands while a kinetic run is in progress, so this sets
    the serial line up and sends nothing.
    """
    await self.io.setup()
    await self.io.set_baudrate(_BAUDRATE)
    await self.io.set_line_property(8, 2, 0)
    await self.io.set_latency_timer(16)
    await self.io.set_flowctrl(_SIO_RTS_CTS_HS)
    await self.io.set_rts(True)
    self._shaking = False
    self._shaking_task = None

  async def resume_kinetic(
    self,
    plate: Plate,
    steps: KineticSteps,
    interval: float,
    reads_received: int = 0,
  ) -> AsyncIterator[KineticRead]:
    """Reconnect to a kinetic run in progress and yield its remaining reads.

    Use this on a new, unconnected instance instead of `setup` after the connection to a run
    started with `run_kinetic` was lost (process restart, sleep, unplugged cable). The reads
    the reader took while disconnected arrive first, then the rest of the run. Once the run
    has ended, the reader answers commands again.

    Args:
      plate: Plate passed to `run_kinetic`.
      steps: Steps passed to `run_kinetic`.
      interval: Interval passed to `run_kinetic`, in seconds.
      reads_received: Number of `KineticRead`s already received for this run, so cycles and
        steps continue from there.
    """
    _, expanded = self._expand_steps(steps, None)
    logger.info(
      "[BioTek %s] resume_kinetic: plate=%s, reads already received=%d",
      self.io.device_id,
      plate.name,
      reads_received,
    )
    await self._attach()
    async for read in self._stream_kinetic(plate, expanded, interval, reads_received):
      yield read

  async def _abort(self) -> None:
    """Abort the running assay and wait until the reader reports that it has stopped.

    The Cytation 7 answers `x` with `\\x10<status>\\x03` once it has stopped, about 14 s after
    aborting a shake, and answers no other command until then. A non-zero status, such as
    after aborting a kinetic run, is logged.
    """
    await self.send_command("x", wait_for_response=False)
    reply = await self._read_until(b"\x03", timeout=60)
    status = reply[reply.rfind(b"\x10") :]
    if status != b"\x100000\x03":
      logger.warning("%s abort reply: %r", self._model_name, status)

  async def stop_kinetic(self) -> None:
    """Abort a kinetic run and wait until the reader has stopped.

    Leave the `run_kinetic` or `resume_kinetic` loop first: the abort reply arrives on the same
    stream. Without a run in progress, the reader does not reply and this raises `TimeoutError`.
    """
    await self._abort()

  async def setup(self) -> None:
    """Warn about unverified C7 support and initialize the shared reader transport."""
    logger.warning(
      "%s plate-reader support has not been verified on hardware. "
      "Please contribute operation-specific results and firmware details after testing.",
      self._model_name,
    )
    await super().setup()
