import astropy.table as at
import matplotlib.pyplot as plt
from CMD import Isochrone


def main() -> int:
    fn = "./joined-GD-1.fits"
    cat = at.Table.read(fn)

    sky_poly = [[-20, -2], [-40, -2], [-40, 2], [-20, 2]]

    pm_poly = [[-9, -2], [-9, 0.5], [-4, 1.5], [-4, -1]]

    o = Isochrone(
        cat,
        age=10.00,  # Gyr
        feh=-1.5,
        distance=8.3,  # kpc
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
        "./gd1_testcmd.png",
        dpi=300,
        bbox_inches="tight",
    )

    iso_patch, iso_mask, iso_model, iso_low, iso_high = o.simpleSln(
        0.1, 15, mass_thresh=0.83
    )
    print(o.x_shift)
    print(o.y_shift)
    # plt.figure()
    # plt.plot(iso_patch[:,0], iso_patch[:,1], '--k')
    # plt.gca().invert_yaxis()
    # plt.show()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
