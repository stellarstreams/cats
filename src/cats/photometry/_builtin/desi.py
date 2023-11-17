"""DES year 6 Photometric Survey."""

from __future__ import annotations

__all__ = ["DESY6Phot"]

from typing import TYPE_CHECKING, ClassVar, TypedDict

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import QTable

from cats.photometry._base import AbstractPhotometricSurvey

if TYPE_CHECKING:
    from dustmaps.map_base import DustMap
    from numpy import bool_
    from numpy.typing import NDArray


class DESY6BandNames(TypedDict):
    """DESY6 band names."""

    WAVG_MAG_PSF_G: str
    WAVG_MAG_PSF_R: str


class DESY6ExtinctionCoeffs(TypedDict):
    """DESY6 extinction coefficients."""

    g: float
    r: float


class DESY6Phot(AbstractPhotometricSurvey):
    """DESY6 Photometric Survey."""

    band_names: ClassVar[DESY6BandNames] = {
        "WAVG_MAG_PSF_G": "g",
        "WAVG_MAG_PSF_R": "r",
    }
    # Schlafly+2011, Rv=3.1
    extinction_coeffs: ClassVar[DESY6ExtinctionCoeffs] = {"g": 3.237, "r": 2.176}
    custom_extinction: ClassVar[bool] = True

    def get_skycoord(self) -> SkyCoord:
        return SkyCoord(self.data["RA"] * u.deg, self.data["DEC"] * u.deg)

    def get_star_mask(self) -> NDArray[bool_]:
        return (self.data["EXT_FITVD"] >= 0) & (self.data["EXT_FITVD"] < 2)

    def get_ext_corrected_phot(
        self, dustmaps_cls: tuple[DustMap] | None = None
    ) -> QTable:
        if dustmaps_cls is None:
            dustmaps_cls = self.dustmaps_cls

        c = self.get_skycoord()
        ebv = dustmaps_cls().query(c)

        tbl = QTable()
        for short_name in self.band_names.values():
            # Compute extinction correction
            Ax = self.extinction_coeffs[short_name] * ebv

            # Apply correction
            tbl[f"A_{short_name}"] = Ax
            tbl[f"{short_name}0"] = self.data[f"BDF_MAG_{short_name.upper()}_CORRECTED"]

        # Metadata
        tbl.meta["band_names"] = [f"{k}0" for k in self.band_names.values()]
        tbl.meta["dustmap"] = dustmaps_cls.__class__.__name__

        return tbl
