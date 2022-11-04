"""
Assumptions and helpers for working with coordinates
"""
import astropy.coordinates as coord
import astropy.units as u

galcen_frame = coord.Galactocentric(
    galcen_distance=8.275 * u.kpc, galcen_v_sun=[8.4, 251.8, 8.4] * u.km / u.s
)
