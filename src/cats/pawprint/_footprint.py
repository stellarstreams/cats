"""Footprint class."""

from __future__ import annotations

__all__ = ["Footprint2D"]

from typing import TYPE_CHECKING, Any

import astropy.table as apt
import numpy as np
from astropy.coordinates import SkyCoord
from matplotlib.path import Path as mpl_path

if TYPE_CHECKING:
    from astropy.coordinates import BaseCoordinateFrame
    from numpy import bool_
    from numpy.typing import NDArray
    from typing_extensions import Self


class Footprint2D(dict):
    """A 2D footprint."""

    def __init__(
        self,
        vertex_coordinates: Any,
        footprint_type: Any,
        stream_frame: BaseCoordinateFrame | None = None,
    ) -> None:
        if footprint_type == "sky":
            if isinstance(vertex_coordinates, SkyCoord):
                vc = vertex_coordinates
            else:
                vc = SkyCoord(vertex_coordinates)
            self.edges = vc
            self.vertices = np.array(
                [vc.transform_to(stream_frame).phi1, vc.transform_to(stream_frame).phi2]
            ).T

        elif footprint_type == "cartesian":
            self.edges = vertex_coordinates
            self.vertices = vertex_coordinates

        self.stream_frame = stream_frame
        self.footprint_type = footprint_type
        self.footprint = mpl_path(self.vertices)

    # ===============================================================

    @classmethod
    def from_vertices(
        cls: type[Self], vertex_coordinates: Any, footprint_type: Any
    ) -> Self:
        """Initialize from vertices."""
        return cls(vertex_coordinates, footprint_type)

    @classmethod
    def from_box(
        cls: type[Self],
        min1: float,
        max1: float,
        min2: float,
        max2: float,
        footprint_type: str,
    ) -> Self:
        """Initialize from a box."""
        vertices = get_vertices_from_box(min1, max1, min2, max2)
        return cls(vertices, footprint_type)

    @classmethod
    def from_file(cls: type[Self], fname: str) -> Self:
        """Initialize from a file."""
        with apt.Table.read(fname) as t:
            vertices = t["vertices"]
            footprint_type = t["footprint_type"]
        return cls(vertices, footprint_type)

    # ===============================================================

    def inside_footprint(self, data: SkyCoord | Any) -> NDArray[bool_] | None:
        """Check if a point is inside the footprint."""
        if isinstance(data, SkyCoord):
            if self.stream_frame is None:
                print("can't!")
                return None

            pts = np.array(
                [
                    data.transform_to(self.stream_frame).phi1.value,
                    data.transform_to(self.stream_frame).phi2.value,
                ]
            ).T
            return self.footprint.contains_points(pts)

        return self.footprint.contains_points(data)

    def export(self) -> dict[str, Any]:
        """Export the footprint to a dictionary."""
        data = {}
        data["stream_frame"] = self.stream_frame
        data["vertices"] = self.vertices
        data["footprint_type"] = self.footprint_type
        return data


def get_vertices_from_box(
    min1: float, max1: float, min2: float, max2: float
) -> list[list[float]]:
    """Get vertices from a box."""
    return [[min1, min2], [min1, max2], [max1, min2], [max1, max2]]
