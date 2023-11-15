from __future__ import annotations

__all__ = ["Footprint2D"]

import astropy.table as apt
import numpy as np
from astropy.coordinates import SkyCoord
from matplotlib.path import Path as mpl_path


class Footprint2D(dict):
    def __init__(self, vertex_coordinates, footprint_type, stream_frame=None):
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

    @classmethod
    def from_vertices(cls, vertex_coordinates, footprint_type):
        return cls(vertex_coordinates, footprint_type)

    @classmethod
    def from_box(cls, min1, max1, min2, max2, footprint_type):
        vertices = cls.get_vertices_from_box(min1, max1, min2, max2)
        return cls(vertices, footprint_type)

    @classmethod
    def from_file(cls, fname):
        with apt.Table.read(fname) as t:
            vertices = t["vertices"]
            footprint_type = t["footprint_type"]
        return cls(vertices, footprint_type)

    def get_vertices_from_box(self, min1, max1, min2, max2):
        return [[min1, min2], [min1, max2], [max1, min2], [max1, max2]]

    def inside_footprint(self, data):
        if isinstance(data, SkyCoord):
            if self.stream_frame is None:
                print("can't!")
                return None
            else:
                pts = np.array(
                    [
                        data.transform_to(self.stream_frame).phi1.value,
                        data.transform_to(self.stream_frame).phi2.value,
                    ]
                ).T
                return self.footprint.contains_points(pts)
        else:
            return self.footprint.contains_points(data)

    def export(self):
        data = {}
        data["stream_frame"] = self.stream_frame
        data["vertices"] = self.vertices
        data["footprint_type"] = self.footprint_type
        return data
