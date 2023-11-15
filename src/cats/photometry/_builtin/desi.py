from __future__ import annotations

__all__ = ["DESY6Phot"]

from typing import ClassVar, TypedDict

import astropy.coordinates as coord
import astropy.units as u
from astropy.table import QTable

from cats.photometry._base import AbstractPhotometricSurvey


class DESY6BandNames(TypedDict):
    WAVG_MAG_PSF_G: str
    WAVG_MAG_PSF_R: str


class DESY6ExtinctionCoeffs(TypedDict):
    g: float
    r: float


class DESY6Phot(AbstractPhotometricSurvey):
    band_names: ClassVar[DESY6BandNames] = {
        "WAVG_MAG_PSF_G": "g",
        "WAVG_MAG_PSF_R": "r",
    }
    # Schlafly+2011, Rv=3.1
    extinction_coeffs: ClassVar[DESY6ExtinctionCoeffs] = {"g": 3.237, "r": 2.176}
    custom_extinction: ClassVar[bool] = True

    def get_skycoord(self):
        return coord.SkyCoord(self.data["RA"] * u.deg, self.data["DEC"] * u.deg)

    def get_star_mask(self):
        return (self.data["EXT_FITVD"] >= 0) & (self.data["EXT_FITVD"] < 2)

    def get_ext_corrected_phot(self, dustmaps_cls=None):
        if dustmaps_cls is None:
            dustmaps_cls = self.dustmaps_cls

        c = self.get_skycoord()
        ebv = dustmaps_cls().query(c)

        tbl = QTable()
        new_band_names = []
        for short_name in self.band_names.values():
            Ax = self.extinction_coeffs[short_name] * ebv
            tbl[f"A_{short_name}"] = Ax
            tbl[f"{short_name}0"] = self.data[f"BDF_MAG_{short_name.upper()}_CORRECTED"]

            #
            new_band_names.append(f"{short_name}0")

        # Metadata
        tbl.meta["band_names"] = new_band_names
        tbl.meta["dustmap"] = dustmaps_cls.__class__.__name__

        return tbl
