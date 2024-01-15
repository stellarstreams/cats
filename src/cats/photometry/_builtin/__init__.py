"""CMD functions."""

from __future__ import annotations

from . import desi, gaia, ps1
from .desi import *
from .gaia import *
from .ps1 import *

__all__: list[str] = []
__all__ += desi.__all__
__all__ += gaia.__all__
__all__ += ps1.__all__
