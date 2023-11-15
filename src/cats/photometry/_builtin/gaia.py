from __future__ import annotations

__all__ = ["GaiaDR3Phot"]

from typing import ClassVar, TypedDict

import numpy as np
from astropy.table import QTable
from pyia import GaiaData

from cats.photometry._base import AbstractPhotometricSurvey


class GaiaDR3BandNames(TypedDict):
    phot_g_mean_mag: str
    phot_bp_mean_mag: str
    phot_rp_mean_mag: str


class GaiaDR3Phot(AbstractPhotometricSurvey):
    band_names: ClassVar[GaiaDR3BandNames] = {
        "phot_g_mean_mag": "G",
        "phot_bp_mean_mag": "BP",
        "phot_rp_mean_mag": "RP",
    }
    custom_extinction: ClassVar[bool] = True

    def get_skycoord(self):
        return GaiaData(self.data).get_skycoord(distance=False)

    def get_star_mask(self):
        return np.ones(len(self.data), dtype=bool)

    def get_ext_corrected_phot(self, dustmaps_cls=None):
        if dustmaps_cls is None:
            dustmaps_cls = self.dustmaps_cls

        g = GaiaData(self.data)
        As = g.get_ext(dustmaps_cls=self.dustmaps_cls)
        As = {"G": As[0], "BP": As[1], "RP": As[2]}  # NOTE: assumption!

        tbl = QTable()
        short_band_names: list[str] = []
        for band, short_name in self.band_names.items():
            # Get extinction coefficient
            Ax = As[short_name]
            if hasattr(Ax, "value"):
                Ax = Ax.value

            # Add to table
            tbl[f"A_{short_name}"] = Ax
            tbl[f"{short_name}0"] = self.data[band] - Ax

            # Record new band name
            short_band_names.append(f"{short_name}0")

        # Metadata
        tbl.meta["band_names"] = short_band_names
        tbl.meta["dustmap"] = dustmaps_cls.__class__.__name__

        return tbl