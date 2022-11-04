import os

import asdf
import astropy.table as apt
import astropy.units as u
import galstreams as gst
import numpy as np
from astropy.coordinates import SkyCoord
from gala.coordinates import GreatCircleICRSFrame
from matplotlib.path import Path as mpl_path


class Footprint2D:
    def __init__(self, vertex_coordinates, footprint_type, stream_frame=None):
        if footprint_type == "sky":
            if isinstance(vertex_coordinates, SkyCoord):
                vc = vertex_coordinates
            else:
                # check if the vertices come with units
                has_units = True
                for v in vertex_coordinates:
                    has_units *= isinstance(v[0], u.Quantity)
                    has_units *= isinstance(v[1], u.Quantity)
                if has_units:
                    # assume coordinates are in stream frame
                    vc = SkyCoord(vertex_coordinates, frame=stream_frame)
                else:
                    # assume units are degrees and frame is phi1/phi2
                    vc = SkyCoord(vertex_coordinates, unit="deg", frame=stream_frame)
            self.edges = vc
            self.vertices = np.array(
                [vc.transform_to(stream_frame).phi1, vc.transform_to(stream_frame).phi2]
            ).T

        # HELP - how do we do PMs properly?
        elif footprint_type == "pm":
            if isinstance(vertex_coordinates, SkyCoord):
                vc = vertex_coordinates
            else:
                # check if the vertices come with units
                has_units = True
                for v in vertex_coordinates:
                    has_units *= isinstance(v[0], u.Quantity)
                    has_units *= isinstance(v[1], u.Quantity)
                if has_units:
                    # assume coordinates are in stream frame
                    vc = SkyCoord(vertex_coordinates, frame=stream_frame)
                else:
                    # assume units are mas/yr and frame is phi1/phi2
                    vc = SkyCoord(
                        stream_frame,
                        lat=vertex_coordinates[0],
                        lon=vertex_coordinates[1],
                        unit="mas/yr",
                    )
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
    def from_vertices(cls, vertex_coordinates, footprint_type, stream_frame=None):
        return cls(vertices, footprint_type, stream_frame)

    @classmethod
    def from_box(cls, min1, max1, min2, max2, footprint_type, stream_frame=None):
        def get_vertices_from_box(min1, max1, min2, max2):
            return [[min1, min2], [min1, max2], [max1, min2], [max1, max2]]

        vertices = get_vertices_from_box(min1, max1, min2, max2)
        return cls(vertices, footprint_type, stream_frame)

    @classmethod
    def from_file(cls, fname):
        with Table.read(fname) as t:
            vertices = t["vertices"]
            footprint_type = t["footprint_type"]
        return cls(vertices, footprint_type)

    def inside_footprint(self, data):
        if isinstance(data, SkyCoord):
            if self.stream_frame is None:
                print("can't!")  # yeah this is my error catching right now
                return
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


class Pawprint(dict):
    """Dictionary class to store a "pawprint":
    polygons in multiple observational spaces that define the initial selection
    used for stream track modeling,
    membership calculation / density modeling, and background modeling.

    New convention: everything is in phi1 phi2 (don't cross the streams)

    """

    def __init__(self, data):

        self.stream_name = data["stream_name"]
        self.pawprint_ID = data["pawprint_ID"]
        self.stream_frame = data["stream_frame"]
        self.width = data["width"]
        self.skyprint = {
            "stream": Footprint2D(
                data["stream_vertices"],
                footprint_type="sky",
                stream_frame=self.stream_frame,
            ),
            "background": Footprint2D(
                data["background_vertices"],
                footprint_type="sky",
                stream_frame=self.stream_frame,
            ),
        }
        # print(self.skyprint)
        # WG3: how to implement distance dependence in isochrone selections?
        self.cmd_filters = data["cmd_filters"]

        if self.cmd_filters is not None:
            self.cmdprint = {}
            for k in data.cmd_filters.keys():
                self.cmdprint[k] = Footprint2D(
                    data["cmd_vertices"][k], footprint_type="cartesian"
                )
        else:
            self.cmdprint = None
        if data["pm_vertices"] is not None:
            self.pmprint = {}
            for k in data["pm_vertices"].keys():
                self.pmprint[k] = Footprint2D(
                    data["pm_vertices"][k],
                    footprint_type="sky",
                    stream_frame=self.stream_frame,
                )  # polygon(s) in proper-motion space mu_phi1, mu_phi2
        else:
            self.pmprint = None

        self.track = data["track"]

    @classmethod
    def from_file(cls, fname):
        def _make_track_file_name(stream_name, pawprint_ID):
            return (
                galstreams_tracks
                + "track.st."
                + stream_name
                + "."
                + pawprint_ID
                + ".ecsv"
            )

        def _make_summary_file_name(stream_name, pawprint_ID):
            return (
                galstreams_tracks
                + "track.st."
                + stream_name
                + "."
                + pawprint_ID
                + ".summary.ecsv"
            )

        import asdf

        data = {}
        with asdf.open(fname, copy_arrays=True) as a:
            # first transfer the stuff that goes directly
            data["stream_name"] = a["stream_name"]
            data["pawprint_ID"] = a["pawprint_ID"]
            data["stream_frame"] = a["stream_frame"]
            data["width"] = a["width"]
            data["cmd_filters"] = a["cmd_filters"]

            # now create footprints from vertices
            data["stream_vertices"] = a["on_stream"]["sky"]["vertices"][:]
            data["background_vertices"] = a["off_stream"]["vertices"][:]

            if a["on_stream"]["cmd"] is not None:
                data["cmd_vertices"] = dict(
                    [
                        (k, a["on_stream"]["cmd"][k]["vertices"])
                        for k in a["on_stream"]["cmd"].keys()
                    ]
                )

            if a["on_stream"]["pm"] is not None:
                data["pm_vertices"] = dict(
                    [
                        (k, a["on_stream"]["pm"][k]["vertices"])
                        for k in a["on_stream"]["pm"].keys()
                    ]
                )

            # right now getting track from galstreams since I can't save it yet
            galstreams_dir = os.path.dirname(gst.__file__)
            galstreams_tracks = os.path.join(galstreams_dir, "tracks/")
            track_file = _make_track_file_name(data["stream_name"], data["pawprint_ID"])
            summary_file = _make_summary_file_name(
                data["stream_name"], data["pawprint_ID"]
            )
            data["track"] = gst.Track6D(
                stream_name=data["stream_name"],
                track_name=data["pawprint_ID"],
                track_file=track_file,
                summary_file=summary_file,
            )

        return cls(data)

    @classmethod
    def from_galstreams(cls, stream_name, pawprint_ID):

        galstreams_dir = os.path.dirname(gst.__file__)
        galstreams_tracks = os.path.join(galstreams_dir, "tracks/")

        def _make_track_file_name(stream_name, pawprint_ID):
            return (
                galstreams_tracks
                + "track.st."
                + stream_name
                + "."
                + pawprint_ID
                + ".ecsv"
            )

        def _make_summary_file_name(stream_name, pawprint_ID):
            return (
                galstreams_tracks
                + "track.st."
                + stream_name
                + "."
                + pawprint_ID
                + ".summary.ecsv"
            )

        def _get_stream_frame_from_file(summary_file):
            t = apt.QTable.read(summary_file)

            x = dict()
            atts = [x.replace("mid.", "") for x in t.keys() if "mid" in x]
            for (
                att
            ) in (
                atts
            ):  # we're effectively looping over skycoords defined for mid here (ra, dec, ...)
                x[att] = t[f"mid.{att}"][
                    0
                ]  # <- make sure to set it up as a scalar. if not, frame conversions get into trouble
            mid_point = SkyCoord(**x)

            x = dict()
            atts = [x.replace("pole.", "") for x in t.keys() if "pole" in x]
            for (
                att
            ) in (
                atts
            ):  # we're effectively looping over skycoords defined for pole here (ra, dec, ...)
                x[att] = t[f"pole.{att}"][0]
            # Make sure to set the pole's distance attribute to 1 (zero causes problems, when transforming to stream frame coords)
            x["distance"] = (
                1.0 * u.kpc
            )  # it shouldn't matter, but if it's zero it does crazy things
            mid_pole = SkyCoord(**x)

            return GreatCircleICRSFrame(pole=mid_pole, ra0=mid_point.icrs.ra)

        def get_recommended_stream_width(stream_name):
            """as recommended by Cecilia.
            eventually pulled from galstreams as an attribute of the track"""
            if "Jhelum" in stream_name:
                if "Jhelum-a" in stream_name:
                    return 0.4 * u.deg
                else:
                    return 0.94 * u.deg
            elif "Fjorm" in stream_name:
                return 0.9 * u.deg
            elif "Pal-5" in stream_name:
                return 0.5 * u.deg
            elif "GD-1" in stream_name:
                return 0.53 * u.deg
            elif "PS1-A" in stream_name:
                return 0.45 * u.deg
            else:  # default
                return 1.0 * u.deg

        data = {}
        data["stream_name"] = stream_name
        data["pawprint_ID"] = pawprint_ID

        track_file = _make_track_file_name(stream_name, pawprint_ID)
        summary_file = _make_summary_file_name(stream_name, pawprint_ID)
        data["stream_frame"] = _get_stream_frame_from_file(summary_file)

        data["track"] = gst.Track6D(
            stream_name=data["stream_name"],
            track_name=data["pawprint_ID"],
            track_file=track_file,
            summary_file=summary_file,
        )
        data["width"] = 2.0 * get_recommended_stream_width(stream_name)
        data["stream_vertices"] = data["track"].create_sky_polygon_footprint_from_track(
            width=data["width"], phi2_offset=0.0 * u.deg
        )
        data["background_vertices"] = data[
            "track"
        ].create_sky_polygon_footprint_from_track(
            width=data["width"], phi2_offset=3.0 * u.deg
        )
        data["cmd_filters"] = None
        data["cmd_vertices"] = None
        data["pm_vertices"] = None

        return cls(data)

    def add_cmd_footprint(self, new_footprint, color, mag, name):
        if self.cmd_filters is None:
            self.cmd_filters = {}
            self.cmdprint = {}

        self.cmd_filters[name] = [color, mag]
        self.cmdprint[name] = new_footprint

    def add_pm_footprint(self, new_footprint, name):
        if self.pmprint is None:
            self.pmprint = {}
        self.pmprint[name] = new_footprint

    def save_pawprint(self):
        # WARNING this doesn't save the track yet - need schema
        # WARNING the stream frame doesn't save right either
        fname = self.stream_name + "." + self.pawprint_ID + ".asdf"
        tree = {
            "stream_name": self.stream_name,
            "pawprint_ID": self.pawprint_ID,
            "stream_frame": self.stream_frame,  # needs a schema to save properly
            "cmd_filters": self.cmd_filters,
            "width": self.width,
            "on_stream": {"sky": self.skyprint["stream"].export()},
            "off_stream": self.skyprint["background"].export(),
            #    'track':self.track   #TODO
        }
        if self.cmdprint is not None:
            tree["on_stream"]["cmd"] = dict(
                [(k, self.cmdprint[k].export()) for k in self.cmd_filters.keys()]
            )
        else:
            tree["on_stream"]["cmd"] = None

        if self.pmprint is not None:
            tree["on_stream"]["pm"] = dict(
                [(k, self.pmprint[k].export()) for k in self.pmprint.keys()]
            )
        else:
            tree["on_stream"]["pm"] = None

        out = asdf.AsdfFile(tree)
        out.write_to(fname)
