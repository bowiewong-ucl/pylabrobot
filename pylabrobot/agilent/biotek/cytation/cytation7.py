"""Agilent BioTek Cytation 7 plate-reader integration."""

from __future__ import annotations

import logging

from pylabrobot.agilent.biotek.cytation.base import _CytationBase

logger = logging.getLogger(__name__)


class Cytation7(_CytationBase):
  """Agilent BioTek Cytation 7 using the shared BioTek plate-reader commands.

  Inherited reader operations have not been verified on Cytation 7 hardware.
  Microscopy is not initialized, and temperature setting remains disabled by the
  inherited heating/cooling flags. Wavelengths, focal heights, and shaking settings
  retain the shared driver's defaults and limits pending hardware verification.
  """

  _model_name = "Agilent BioTek Cytation 7"

  async def setup(self) -> None:
    """Warn about unverified C7 support and initialize the shared reader transport."""
    logger.warning(
      "%s plate-reader support has not been verified on hardware. "
      "Please contribute operation-specific results and firmware details after testing.",
      self._model_name,
    )
    await super().setup()
