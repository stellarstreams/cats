"""Combine PM and CMD cuts."""

from __future__ import annotations

from typing import Any

import astropy.table as at
import matplotlib.pyplot as plt
import pandas as pd

from cats.cmd import Isochrone
from cats.pawprint import Footprint2D, Pawprint

plt.rc("xtick", top=True, direction="in", labelsize=15)
plt.rc("ytick", right=True, direction="in", labelsize=15)
plt.rc("font", family="Arial")


def generate_isochrone_vertices(
    cat: Any, sky_poly: Any, pm_poly: Any, config: Any
) -> Any:
    """Generate Isochrone Vertices.

    Use the generated class to make a new polygon for the given catalog in CMD
    space given a sky and PM polygon.
    """
    o = Isochrone(
        cat,
        age=config.isochrone.age,  # Gyr
        feh=config.isochrone.feh,
        distance=config.isochrone.distance,  # kpc
        alpha=config.isochrone.alpha,
        sky_poly=sky_poly,
        pm_poly=pm_poly,
    )

    o.generate_isochrone()
    o.sel_sky()
    o.sel_pm()
    o.data_cmd(0.03, 0.2)
    o.correct_isochrone()

    fig = o.plot_CMD(tolerance=0.1)
    fig.savefig(
        "./jhelum_testcmd.png",
        dpi=300,
        bbox_inches="tight",
    )

    return o.simpleSln(0.1, 15, mass_thresh=0.83)[0]


def generate_pm_vertices() -> list[list[float]]:
    """Generate Proper Motion Vertices.

    Use the generated class to make a new polygon for the given catalog in PM
    space given a sky and CMD polygon.
    """
    return [[-7.0, 0.0], [-5.0, 0.0], [-5.0, 1.6], [-7.0, -1.6]]


def load_sky_region(fn: Any) -> tuple[list[float], list[float]]:
    """Load Sky Region."""
    sky_print = [
        [-5, -2],
        [+5, -2],
        [+5, +2],
        [-5, +2],
    ]
    bg_print = []
    return sky_print, bg_print


def main() -> int:
    """Run Script."""
    # load in config file, catalog from filename
    config = pd.read_json("config.json")
    cat = at.Table.read(config.streaminfo.cat_fn)

    # load in file with the sky footprint.
    sky_poly, _ = load_sky_region(config.streaminfo.sky_print)

    # have an initial selection for the PM region that is very wide
    # this could also be stored in a footprint
    pm_poly = config.proper_motion.init_region

    # Initialize a new pawprint for the stream
    p = Pawprint.from_galstreams(
        config.streaminfo.gst_name,
        config.streaminfo.gst_track,
    )
    ##### probably add the initial cuts #####

    # generate first iso mask, add to pawprint
    iso_vertices_1 = generate_isochrone_vertices(
        cat,
        sky_poly,
        pm_poly,
        config,
    )
    cmd_poly = Footprint2D.from_vertices(
        iso_vertices_1,
        footprint_type="cartesian",
    )
    p.add_cmd_footprint(
        cmd_poly,
        "g_r",
        "g",
        "iso_mask1",
    )

    # use this cut for a new PM cut
    pm_vertices_1 = generate_pm_vertices(
        cat,
        sky_poly,
        iso_vertices_1,
        config,
    )
    pm_poly = Footprint2D.from_vertices(
        pm_vertices_1,
        footprint_type="cartesian",
    )
    p.add_pm_footprint(
        pm_poly,
        "pm_mask1",
    )

    # generate refined iso mask, add to pawprint
    iso_vertices_2 = generate_isochrone_vertices(
        cat,
        sky_poly,
        pm_vertices_2,
        config,
    )
    cmd_poly = Footprint2D.from_vertices(
        iso_vertices_2,
        footprint_type="cartesian",
    )
    p.add_cmd_footprint(
        cmd_poly,
        "g_r",
        "g",
        "iso_mask2",
    )

    p.save_pawprint()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
