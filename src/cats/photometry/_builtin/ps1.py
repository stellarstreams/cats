from __future__ import annotations

__all__ = ["PS1Phot"]

from typing import ClassVar, TypedDict

import astropy.units as u
import numpy as np
import numpy.typing as npt
from astropy.coordinates import SkyCoord

from cats.photometry._base import AbstractPhotometricSurvey


class PS1BandNames(TypedDict):
    gMeanPSFMag: str
    rMeanPSFMag: str
    iMeanPSFMag: str
    zMeanPSFMag: str
    yMeanPSFMag: str


class PS1ExtinctionCoeffs(TypedDict):
    g: str
    r: str
    i: str
    z: str
    y: str


class PS1Phot(AbstractPhotometricSurvey):
    band_names: ClassVar[PS1BandNames] = {
        "gMeanPSFMag": "g",
        "rMeanPSFMag": "r",
        "iMeanPSFMag": "i",
        "zMeanPSFMag": "z",
        "yMeanPSFMag": "y",
    }

    # Schlafly+2011, Rv=3.1
    # TODO: load from a config file
    extinction_coeffs: ClassVar[PS1ExtinctionCoeffs] = {
        "g": 3.172,
        "r": 2.271,
        "i": 1.682,
        "z": 1.322,
        "y": 1.087,
    }

    def get_skycoord(self) -> SkyCoord:
        return SkyCoord(
            self.data["raMean"] << u.deg,
            self.data["decMean"] << u.deg,
            frame="icrs",
        )

    def get_star_mask(self) -> npt.NDArray[np.bool_]:
        """Star/galaxy separation for PS1.

        See:
        https://outerspace.stsci.edu/display/PANSTARRS/How+to+separate+stars+and+galaxies

        Returns
        -------
        star_mask : `numpy.ndarray`
            True where the stars are.
        """
        return (self.data["iMeanPSFMag"] - self.data["iMeanKronMag"]) < 0.05