from __future__ import annotations

from typing import Any

import astropy.units as u
from astropy.coordinates import Distance

stream_inputs: dict[str, dict[str, Any]] = {}
stream_inputs["GD-1"] = {
    # galstream stuff
    "short_name": "GD-1",
    "pawprint_id": "pricewhelan2018",
    # stream stuff
    "width": 2.0,  # full width in degrees (add units in pawprint)
    # data stuff
    "phot_survey": "PS1",
    "band1": "g",
    "band2": "r",
    "mag": "g0",
    "color1": "g0",
    "color2": "r0",
    "minmag": 16.0 * u.mag,
    "maxmag": 24.0 * u.mag,
    # isochrone stuff
    "age": 11.8 * u.Gyr,
    "feh": -1.5,
    "distance": Distance(8.3, u.kpc),
    "turnoff": 17.8 * u.mag,  # mag of MS turnoff
    "alpha": 0,  # don't think we actually use this
    "scale_err": 2,
    "base_err": 0.075,
    "bin_sizes": [0.03, 0.2],  # xbin and ybin width for CMD
}
stream_inputs["Pal5"] = {
    # galstream stuff
    "short_name": "Pal5",
    "pawprint_id": "pricewhelan2019",
    # stream stuff
    "width": 1.0,  # degrees (add units in pawprint)
    # data stuff
    "phot_survey": "PS1",
    "band1": "g",
    "band2": "r",
    "mag": "g0",
    "color1": "g0",
    "color2": "r0",
    "minmag": 16.0 * u.mag,
    "maxmag": 24.0 * u.mag,
    # isochrone stuff
    "age": 12 * u.Gyr,
    "feh": -1.4,
    "distance": Distance(20.9, u.kpc),
    "turnoff": 15 * u.mag,  # mag of MS turnoff
    "alpha": 0,  # don't think we actually use this
    "scale_err": 2,
    "base_err": 0.075,
    "bin_sizes": [0.03, 0.2],
}
stream_inputs["Jhelum"] = {
    # galstream stuff
    "short_name": "Jhelum-b",
    "pawprint_id": "bonaca2019",
    # stream stuff
    "width": 2.0,  # degrees (add units in pawprint)
    # data stuff
    "phot_survey": "des",
    "band1": "g",
    "band2": "r",
    "mag": "g0",
    "color1": "g0",
    "color2": "r0",
    "minmag": 16.0 * u.mag,
    "maxmag": 24.0 * u.mag,
    # isochrone stuff
    "age": 12 * u.Gyr,
    "feh": -1.7,
    "distance": Distance(13.2, u.kpc),
    "turnoff": 18.7 * u.mag,  # mag of MS turnoff
    "alpha": 0,  # don't think we actually use this
    "scale_err": 2,
    "base_err": 0.075,
    "bin_sizes": [0.03, 0.2],
}
stream_inputs["Fjorm-M68"] = {
    # galstream stuff
    "short_name": "M68-Fjorm",
    # "pawprint_id": 'ibata2021',
    "pawprint_id": "palau2019",
    # stream stuff
    "width": 1,  # TOTAL width degrees, recommend 2sigma if known
    # data stuff
    "phot_survey": "Gaia",
    "band1": "BP",
    "band2": "RP",
    "mag": "G0",
    "color1": "BP0",
    "color2": "RP0",
    "minmag": 16.0 * u.mag,
    "maxmag": 24.0 * u.mag,
    # isochrone stuff
    "age": 11.2 * u.Gyr,
    "feh": -2.2,
    "distance": Distance(6, u.kpc),
    "turnoff": 17.0 * u.mag,  # mag of MS turnoff
    "alpha": 0,  # don't think we actually use this
    "scale_err": 2,
    "base_err": 0.075,
    "bin_sizes": [0.03, 0.2],
}
