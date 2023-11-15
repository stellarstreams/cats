from __future__ import annotations

__all__ = ["Pawprint"]

import pathlib

import asdf
import astropy.table as apt
import astropy.units as u
import galstreams as gst
from astropy.coordinates import SkyCoord
from gala.coordinates import GreatCircleICRSFrame


class Pawprint(dict):
    """Dictionary class to store a "pawprint".

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
            self.pmprint = Footprint2D(
                data["pm_vertices"], footprint_type="cartesian"
            )  # polygon(s) in proper-motion space mu_phi1, mu_phi2
        else:
            self.pmprint = None
        if data["pm1_vertices"] is not None:
            self.pm1print = Footprint2D(
                data["pm1_vertices"], footprint_type="cartesian"
            )  # polygon(s) in proper-motion space mu_phi1, mu_phi2
        else:
            self.pm1print = None

        # Need to code this into CMD and proper motion stuff
        if data["pm2_vertices"] is not None:
            self.pm2print = Footprint2D(
                data["pm2_vertices"], footprint_type="cartesian"
            )  # polygon(s) in proper-motion space mu_phi1, mu_phi2
        else:
            self.pm2print = None

        self.track = data["track"]

    @classmethod
    def from_file(cls, fname):
        import asdf

        data = {}
        with asdf.open("fname") as a:
            # first transfer the stuff that goes directly
            data["stream_name"] = a["stream_name"]
            data["pawprint_ID"] = a["pawprint_ID"]
            data["stream_frame"] = a["stream_frame"]
            data["width"] = a["width"]
            data["cmd_filters"] = a["cmd_filters"]

            # now create footprints from vertices
            data["sky"] = {}

            if a["on_stream"]["cmd"] is not None:
                data["cmd_vertices"] = {
                    k: Footprint2D(
                        a["on_stream"]["cmd"][k]["vertices"],
                        a["on_stream"]["cmd"][k],
                    )["footprint_type"]
                    for k in a["on_stream"]["cmd"].keys()
                }

            if a["on_stream"]["pm"] is not None:
                data["cmd_vertices"] = {
                    k: Footprint2D(
                        a["on_stream"]["cmd"][k]["vertices"],
                        a["on_stream"]["cmd"][k],
                    )["footprint_type"]
                    for k in a["on_stream"]["cmd"].keys()
                }

        return cls(data)

    @classmethod
    def pawprint_from_galstreams(cls, stream_name, pawprint_ID, width):
        def _get_stream_frame_from_file(summary_file):
            t = apt.QTable.read(summary_file)

            x = {}
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

            x = {}
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

        data = {}
        data["stream_name"] = stream_name
        data["pawprint_ID"] = pawprint_ID

        galstreams_tracks = pathlib.Path(gst.__file__).resolve().parent / "tracks"
        track_file = galstreams_tracks / f"track.st.{stream_name}.{pawprint_ID}.ecsv"
        summary_file = (
            galstreams_tracks / f"track.st.{stream_name}.{pawprint_ID}.summary.ecsv"
        )
        data["stream_frame"] = _get_stream_frame_from_file(summary_file)

        data["track"] = gst.Track6D(
            stream_name=data["stream_name"],
            track_name=data["pawprint_ID"],
            track_file=track_file,
            summary_file=summary_file,
        )
        try:
            data["width"] = (
                2 * data["track"].track_width["width_phi2"]
            )  # one standard deviation on each side (is this wide enough?)
        except Exception:
            data["width"] = width
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
        data["pm1_vertices"] = None
        data["pm2_vertices"] = None

        return cls(data)

    def add_cmd_footprint(self, new_footprint, color, mag, name):
        if self.cmd_filters is None:
            self.cmd_filters = dict((name, [color, mag]))
            self.cmdprint = dict((name, new_footprint))
        else:
            self.cmd_filters[name] = [color, mag]
            self.cmdprint[name] = new_footprint

    def add_pm_footprint(self, new_footprint, name):
        if self.pmprint is None:
            self.pmprint = dict((name, new_footprint))
        else:
            self.pmprint[name] = new_footprint

    def save_pawprint(self):
        # WARNING this doesn't save the track yet - need schema
        # WARNING the stream frame doesn't save right either
        fname = self.stream_name + self.pawprint_ID + ".asdf"
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
            tree["on_stream"]["cmd"] = {
                k: self.cmdprint[k].export() for k in self.cmd_filters
            }
        else:
            tree["on_stream"]["cmd"] = None

        if self.pmprint is not None:
            tree["on_stream"]["pm"] = {
                k: self.pmprint[k].export() for k in self.pmprint
            }
        else:
            tree["on_stream"]["pm"] = None

        out = asdf.AsdfFile(tree)
        out.write_to(fname)
