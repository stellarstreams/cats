pawprint = pawprintClass("Pal5", "PW19")

# Generate random points on the sphere in an area around the track
rao, raf = np.min(pawprint.track.ra.deg), np.max(pawprint.track.ra.deg)
deco, decf = np.min(pawprint.track.dec.deg), np.max(pawprint.track.dec.deg)
field_ra, field_dec = galstreams.get_random_spherical_angles(
    5000, az=[rao, raf], lat=[deco, decf], degree=True
)
field = coord.SkyCoord(ra=field_ra * u.deg, dec=field_dec * u.deg, frame="icrs")

stars = starClass(field)


plt.figure(1, figsize=(10, 8))
ax = plt.subplot(111)


# Plot the track
ax.plot(pawprint.track.ra, pawprint.track.dec, ".", ms=2.0, color="C0")
# plot the default polygon track stored in the library
ax.plot(
    pawprint.skyprint["stream"].icrs.ra,
    pawprint.skyprint["stream"].icrs.dec,
    ls="--",
    color="C0",
)
# Random background "field" points
ax.plot(stars.ra, stars.dec, "k.", ms=0.5)

# Select the field points inside the polygon footprint
on = stars.makeMask(pawprint, what="sky.stream")  # is a function of starlist
ax.plot(stars.ra[on], stars.dec[on], ".", ms=2.5, color="C0")

# Create a new polygon footprint off-stream, with a given offset and width, and select field points inside it
# <this will now be stored in the pawprint>
# off_poly = mwsts[st].create_sky_polygon_footprint_from_track(width=1.*u.deg, phi2_offset=3.5*u.deg)
off = stars.makeMask(pawprint, what="sky.background")
# Plot the off-stream polygon footprint and points selected inside it
ax.plot(
    pawprint.skyprint["background"].icrs.ra,
    pawprint.skyprint["background"].icrs.dec,
    ls=":",
    color="C1",
)
ax.plot(stars.ra[off], stars.dec[off], ".", ms=2.5, color="C1")

ax.set_xlim(rao, raf)
ax.set_ylim(deco, decf)
ax.set_xlabel("RA (deg)")
ax.set_ylabel("DEC (deg)")
ax.invert_xaxis()
