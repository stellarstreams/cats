from __future__ import annotations

import astropy.table as at
from CMD import Isochrone


def main() -> int:
    fn = "./joined-Jhelum.fits"
    cat = at.Table.read(fn)

    sky_poly = [
        [-5, -2],
        [+5, -2],
        [+5, +2],
        [-5, +2],
    ]
    pm_poly = [
        [6, -6],
        [7.5, -6],
        [7.5, -4.5],
        [6, -4.5],
    ]

    pm_poly = [
        [-7.0, 0.0],
        [-5.0, 0.0],
        [-5.0, 1.6],
        [-7.0, -1.6],
    ]

    o = Isochrone(
        cat,
        age=10.00,  # Gyr
        feh=-1.5,
        distance=13.2,  # kpc
        alpha=0.0,
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

    iso_patch, iso_mask, iso_model, iso_low, iso_high = o.simpleSln(
        0.1, 15, mass_thresh=0.83
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
