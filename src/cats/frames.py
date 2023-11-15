"""Default Coordinate Frames."""

from __future__ import annotations

__all__ = ["galactocentric"]

import astropy.coordinates as coord
import astropy.units as u

galactocentric = coord.Galactocentric(
    galcen_distance=8.275 * u.kpc,
    galcen_v_sun=[8.4, 251.8, 8.4] * u.km / u.s,
)
