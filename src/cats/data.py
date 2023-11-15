"""Helper functions for working with data."""

from __future__ import annotations

__all__ = ["make_astro_photo_joined_data"]

from typing import TYPE_CHECKING

import astropy.units as u
import gala.coordinates as gc
import numpy as np
from astropy.table import hstack, join, unique
from scipy.interpolate import InterpolatedUnivariateSpline

if TYPE_CHECKING:
    from astropy.table import QTable
    from galstreams import Track6D
    from pyia import GaiaData

    from cats.photometry._base import PhotometricSurvey


def make_astro_photo_joined_data(
    gaia_data: GaiaData, phot_data: PhotometricSurvey, track6d: Track6D
) -> QTable:
    """Join Gaia and photometry data, and transform to stream coordinates.

    Parameters
    ----------
    gaia_data : `pyia.GaiaData`
    phot_data : `cats.photometry.PhotometricSurvey`
    track6d : `galstreams.Track6D`

    Returns
    -------
    joined : `astropy.table.Table`
        Joined table of Gaia and photometry data, transformed to stream coordinates.
    """
    stream_frame = track6d.stream_frameame

    # --------------------------------------------
    # Process the track

    track = track6d.track.transform_to(stream_frame)
    if np.all(track.distance.value == 0):
        msg = "A distance track is required: this stream has no distance information."
        raise ValueError(msg)

    # Interpolator to get predicted distance from phi1
    dist_interp = InterpolatedUnivariateSpline(
        track.phi1.degree, track.distance.value, k=1
    )

    # Get stream coordinates for all stars, and reflex correct with predicted distance
    _c_tmp = gaia_data.get_skycoord(distance=False).transform_to(stream_frame)
    c = gaia_data.get_skycoord(
        distance=dist_interp(_c_tmp.phi1.degree) * track.distance.unit,
        radial_velocity=0 * u.km / u.s,
    )
    c_stream = c.transform_to(stream_frame)
    c_stream_refl = gc.reflex_correct(c_stream)

    # --------------------------------------------

    # Get extinction-corrected photometry and star/galaxy mask
    ext = phot_data.get_ext_corrected_phot()
    ext["star_mask"] = phot_data.get_star_mask()

    # Start building the final joined table
    joined = gaia_data.data.copy()
    for name in ("phi1", "phi2"):
        joined[name] = getattr(c_stream_refl, name)
    for name in ("pm_phi1_cosphi2", "pm_phi2"):
        joined[name] = getattr(c_stream_refl, name)
        joined[f"{name}_unrefl"] = getattr(c_stream, name)

    phot_full = hstack([phot_data.data, ext])
    cols = ["source_id", "star_mask"] + [b for b in ext.colnames if b.endswith("0")]
    phot_min = phot_full[cols]

    return unique(join(joined, phot_min, keys="source_id"), keys="source_id")
