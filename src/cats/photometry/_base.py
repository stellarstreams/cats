from __future__ import annotations

__all__ = ["AbstractPhotometricSurvey"]

import abc
from dataclasses import dataclass
from typing import ClassVar

import numpy as np
import numpy.typing as npt
from astropy.coordinates import SkyCoord
from astropy.table import QTable, Table
from dustmaps.map_base import DustMap
from dustmaps.sfd import SFDQuery
from typing_extensions import Self


@dataclass(frozen=True)
class AbstractPhotometricSurvey(metaclass=abc.ABCMeta):
    """Photoemtric survey base class.

    Parameters
    ----------
    data : :class:`~astropy.table.QTable`
        Anything that can be passed into `astropy.table.Table` to construct an
        astropy table instance, or a string filename (that can be read into an
        astropy table instance).
    """

    band_names: ClassVar[dict[str, str]] = {}
    extinction_coeffs: ClassVar[dict[str, float]] = {}
    custom_extinction: ClassVar[bool] = False
    dustmaps_cls: ClassVar[type[DustMap]] = SFDQuery

    data: QTable

    def __init_subclass__(cls) -> None:
        if not cls.band_names:  # empty dict
            msg = (
                "You must define some photometric band names in band_names for any "
                "survey-specific subclass"
            )
            raise ValueError(msg)

        for short_name in cls.band_names.values():
            if not cls.custom_extinction and short_name not in cls.extinction_coeffs:
                msg = (
                    "You must specify extinction coefficients for all photometric "
                    "bands in any survey-specific subclass"
                )
                raise ValueError(msg)

    @classmethod
    def from_tablelike(cls: type[Self], data: str | Table) -> Self:
        if isinstance(data, str):
            return cls(QTable.read(data))
        return cls(data)

    @abc.abstractmethod
    def get_skycoord(self) -> SkyCoord:
        """Return a SkyCoord object from the data table."""

    @abc.abstractmethod
    def get_star_mask(self) -> npt.NDArray[np.bool_]:
        """Star-galaxy separation."""

    def get_ext_corrected_phot(
        self, dustmaps_cls: type[DustMap] | None = None
    ) -> QTable:
        """Get extinction-corrected photometry.

        Parameters
        ----------
        dustmaps_cls : type[:class:`~dustmaps.map_base.DustMap`], optional
            Dustmap class to use for extinction correction. Default is
            :class:`~dustmaps.sfd.SFDQuery`.

        Notes
        -----
        This is a default implementation. Most subclasses will need to override
        this method.
        """
        if self.custom_extinction:
            msg = "TODO"
            raise RuntimeError(msg)

        if dustmaps_cls is None:
            dustmaps_cls = self.dustmaps_cls

        c = self.get_skycoord()
        ebv = dustmaps_cls().query(c)

        tbl = QTable()
        short_band_names: list[str] = []
        for band, short_name in self.band_names.items():
            # Get extinction coefficient
            Ax = self.extinction_coeffs[short_name] * ebv

            # Add to table
            tbl[f"A_{short_name}"] = Ax
            tbl[f"{short_name}0"] = self.data[band] - Ax

            # Record new band name
            short_band_names.append(f"{short_name}0")

        # Metadata
        tbl.meta["band_names"] = short_band_names
        tbl.meta["dustmap"] = dustmaps_cls.__class__.__name__

        return tbl
