import abc

import astropy.coordinates as coord
import astropy.table as at
import astropy.units as u
import numpy as np
from dustmaps.sfd import SFDQuery
from pyia import GaiaData

__all__ = ["PS1Phot", "GaiaDR3Phot", "DESY6Phot"]


class PhotometricSurvey(abc.ABC):
    band_names = {}
    extinction_coeffs = {}
    custom_extinction = False
    dustmaps_cls = SFDQuery

    def __init_subclass__(cls) -> None:

        if len(cls.band_names) == 0:
            raise ValueError(
                "You must define some photometric band names in band_names for any "
                "survey-specific subclass"
            )

        for short_name in cls.band_names.values():
            if not cls.custom_extinction and short_name not in cls.extinction_coeffs:
                raise ValueError(
                    "You must specify extinction coefficients for all photometric "
                    "bands in any survey-specific subclass"
                )

    def __init__(self, data) -> None:
        """

        Parameters
        ----------
        data : table-like, str
            Anything that can be passed into `astropy.table.Table` to construct an
            astropy table instance, or a string filename (that can be read into an
            astropy table instance).
        """
        if isinstance(data, str):
            data = at.Table.read(data)
        self.data = at.Table(data)

    @abc.abstractmethod
    def get_skycoord(self):
        """
        Return a SkyCoord object from the data table.
        """
        pass

    @abc.abstractmethod
    def get_star_mask(self):
        """
        Star-galaxy separation
        """
        pass

    def get_ext_corrected_phot(self, dustmaps_cls=None):
        if self.custom_extinction:
            raise RuntimeError("TODO")

        if dustmaps_cls is None:
            dustmaps_cls = self.dustmaps_cls

        c = self.get_skycoord()
        ebv = dustmaps_cls().query(c)

        tbl = at.Table()
        new_band_names = []
        for band, short_name in self.band_names.items():
            Ax = self.extinction_coeffs[short_name] * ebv
            tbl[f"A_{short_name}"] = Ax
            tbl[f"{short_name}0"] = self.data[band] - Ax
            new_band_names.append(f"{short_name}0")
        tbl.meta["band_names"] = new_band_names
        tbl.meta["dustmap"] = dustmaps_cls.__class__.__name__

        return tbl


class PS1Phot(PhotometricSurvey):
    band_names = {
        "gMeanPSFMag": "g",
        "rMeanPSFMag": "r",
        "iMeanPSFMag": "i",
        "zMeanPSFMag": "z",
        "yMeanPSFMag": "y",
    }

    # Schlafly+2011, Rv=3.1
    extinction_coeffs = {
        "g": 3.172,
        "r": 2.271,
        "i": 1.682,
        "z": 1.322,
        "y": 1.087,
    }

    def get_skycoord(self):
        return coord.SkyCoord(self.data["raMean"] * u.deg, self.data["decMean"] * u.deg)

    def get_star_mask(self):
        """
        Star/galaxy separation for PS1

        See:
        https://outerspace.stsci.edu/display/PANSTARRS/How+to+separate+stars+and+galaxies

        Returns
        -------
        star_mask : `numpy.ndarray`
            True where the stars are.
        """
        d_mag_mask = self.data["iMeanPSFMag"] - self.data["iMeanKronMag"] < 0.05
        return d_mag_mask


class GaiaDR3Phot(PhotometricSurvey):
    band_names = {
        "phot_g_mean_mag": "G",
        "phot_bp_mean_mag": "BP",
        "phot_rp_mean_mag": "RP",
    }
    custom_extinction = True

    def get_skycoord(self):
        return GaiaData(self.data).get_skycoord(distance=False)

    def get_star_mask(self):
        """
        Star-galaxy separation:
        """
        return np.ones(len(self.data), dtype=bool)

    def get_ext_corrected_phot(self, dustmaps_cls=None):
        if dustmaps_cls is None:
            dustmaps_cls = self.dustmaps_cls
        g = GaiaData(self.data)
        As = g.get_ext(dustmaps_cls=self.dustmaps_cls)
        As = {"G": As[0], "BP": As[1], "RP": As[2]}  # NOTE: assumption!

        tbl = at.Table()
        new_band_names = []
        for band, short_name in self.band_names.items():
            Ax = As[short_name]
            if hasattr(Ax, "value"):
                Ax = Ax.value
            tbl[f"A_{short_name}"] = Ax
            tbl[f"{short_name}0"] = self.data[band] - Ax
            new_band_names.append(f"{short_name}0")
        tbl.meta["band_names"] = new_band_names
        tbl.meta["dustmap"] = dustmaps_cls.__class__.__name__

        return tbl


class DESY6Phot(PhotometricSurvey):
    band_names = {
        "WAVG_MAG_PSF_G": "g",
        "WAVG_MAG_PSF_R": "r",
    }
    # Schlafly+2011, Rv=3.1
    extinction_coeffs = {
        "g": 3.237,
        "r": 2.176,
    }
    custom_extinction = True

    def get_skycoord(self):
        return coord.SkyCoord(self.data["RA"] * u.deg, self.data["DEC"] * u.deg)

    def get_star_mask(self):
        """
        Star-galaxy separation:
        """
        return (self.data["EXT_FITVD"] >= 0) & (self.data["EXT_FITVD"] < 2)

    def get_ext_corrected_phot(self, dustmaps_cls=None):
        if dustmaps_cls is None:
            dustmaps_cls = self.dustmaps_cls

        c = self.get_skycoord()
        ebv = dustmaps_cls().query(c)

        tbl = at.Table()
        new_band_names = []
        for short_name in self.band_names.values():
            Ax = self.extinction_coeffs[short_name] * ebv
            tbl[f"A_{short_name}"] = Ax
            tbl[f"{short_name}0"] = self.data[f"BDF_MAG_{short_name.upper()}_CORRECTED"]
            new_band_names.append(f"{short_name}0")
        tbl.meta["band_names"] = new_band_names
        tbl.meta["dustmap"] = dustmaps_cls.__class__.__name__

        return tbl
