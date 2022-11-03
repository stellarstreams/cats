import astropy.table as at
import astropy.units as u
import gala.coordinates as gc
import numpy as np
import scipy.interpolate as sci


def make_astro_photo_joined_data(gaia_data, phot_data, track6d):
    """
    Parameters
    ----------
    gaia_data : `pyia.GaiaData`
    phot_data : `cats.photometry.PhotometricSurvey`
    track6d : `galstreams.Track6D`

    """
    stream_fr = track6d.stream_frame

    track = track6d.track.transform_to(stream_fr)
    if np.all(track.distance.value == 0):
        raise ValueError(
            "A distance track is required: this stream has no distance information in "
            "the galstreams track."
        )

    # interpolator to get predicted distance from phi1
    dist_interp = sci.InterpolatedUnivariateSpline(
        track.phi1.degree, track.distance.value, k=1
    )

    # get stream coordinates for all stars, and reflex correct with predicted distance
    _c_tmp = gaia_data.get_skycoord(distance=False).transform_to(stream_fr)
    c = gaia_data.get_skycoord(
        distance=dist_interp(_c_tmp.phi1.degree) * track.distance.unit,
        radial_velocity=0 * u.km / u.s,
    )
    c_stream = c.transform_to(stream_fr)
    c_stream_refl = gc.reflex_correct(c_stream)

    # get extinction-corrected photometry and star/galaxy mask
    ext = phot_data.get_ext_corrected_phot()
    ext["star_mask"] = phot_data.get_star_mask()

    # start building the final joined table
    joined = gaia_data.data.copy()
    for name in ["phi1", "phi2", "pm_phi1_cosphi2", "pm_phi2"]:
        joined[name] = getattr(c_stream_refl, name)
        if name not in ["phi1", "phi2"]:
            joined[f"{name}_unrefl"] = getattr(c_stream, name)

    phot_full = at.hstack([phot_data.data, ext])
    cols = ["source_id", "star_mask"] + [b for b in ext.colnames if b.endswith("0")]
    phot_min = phot_full[cols]

    joined = at.join(joined, phot_min, keys="source_id")
    joined = at.unique(joined, keys="source_id")

    return joined
