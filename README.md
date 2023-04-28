# CATS 🐈

Community Atlas of Tidal Streams

## Python dependencies

We have included a number of useful Python packages in the [requirements.txt](https://github.com/stellarstreams/cats/blob/main/requirements.txt) file — to install these, you can use `pip`:

    pip install -r requirements.txt

Python 3.11 was *just* released, but not all packages have built wheels for this, so we
recommend using Python 3.9 or 3.10.

## Contributor guide

- Fork the repository
- Clone the fork
- Set up pre-commit locally: `pre-commit install

## Kiyan Updates:

General Note: This pipeline is specifically designed to be applied to streams for which we have all (or at least
most) of the input information. It cannot be used to detect new streams and I would anticipate some issues
if it is used to derive parameters for which we have very broad priors. Luckily, galstreams is fantastic and we have the 
requisite phase space track information for nearly all the streams. And we generally have isochrone information from 
the discovery papers, although not always. Because of this, the process currently relies heavily on galstreams and 
assumes that it is accurate enough to be used as a good starting point for fine-tuning.

I believe there are two goals to using this process to make polygons and cuts:
1. Create a slightly cleaned version of the data that can be used for the membership probability model that Adrian and I have built for GD-1. Note: For GD-1, this used a very rough proper motion cut and a rough ischrone cut so that the overdensity in the sky position and proper motions became visible enough that the algorithm would work.
2. Create polygons that, when used together, give a relatively pure sample of the stream. It should be possible to get an initial guess at the stream members without going through the whole membership probability process.

#### Possible Next Improvements:
1. Use a background subtracted hess diagram when fitting for best isochrone instead of just the data. Just the data works fine once a proper motion cut has been applied in most cases but not beforehand.
2. Find a way to improve the extrapolation of galstreams tracks so that we can use the purest proper motion cuts (#4 below) in regions
3. How can I update one of the galstreams tracks without affecting the other tracks. (E.g. if I want to use the Starkman stream track in ($\phi1, \phi2$) but the Ibata track in the other coordinates.)
4. Using the functions that WG1 had originally written for adding footprints (I am currently just saying `pmprint= ...` instead of `add_new_pmprint(...)` which should improve things.
5. In general, cleaning up the code and making it legible/understandable for the community.
6. Testing on many other streams to see where improvements still need to be made.

#### More specifics on current code:

Making the polygons and stream member cuts is now a fully automated process that uses the pawprint class.
The process takes in the following inputs for each stream:
- Galstreams (uses `pawprint_from_galstreams` function to create initial sky polygon)
    - short_name
    - pawprint_id 
    - stream width (in degrees)
- Isochrone parameters:
    - Age (in Gyrs)
    - Metallicity
    - Distance (in kpc; don't think I actually use this, uses galstreams distance track instead)
    - alpha (not currently used either)
    - base error (half-width of CMD polygon if magnitudes had no errors)
    - scale error (scaling the width of the CMD polygon given errors on the magnitudes from data)
    - bin sizes (CMD color and magnitude bin sizes)
- Data parameters:
    - survey (name of survey)
    - band1 (generic name of band1 (e.g. 'g', 'r', 'i', etc))
    - band2
    - mag1 (specific name of band1 in survey (e.g. 'g0', 'r0', 'i0' for Panstarrs)
    - mag2
    - minmag
    - maxmag
    

The process will automatically make the following polygons to achieve these two objectives:
1. Rough PM cut polygon (see `proper_motions.py` for documentation)
    - Takes the minima and maxima of the two proper motion galstreams tracks
    - Creates a rectangle in proper motion space out of those four numbers
    - Adds a "buffer" to each side of the rectangle (buffer defaults to 2 mas/yr but is tunable)
2. Sky polygon
    - Created in `pawprint` class using galstreams track and the input width
2. Isochrone polygon (see `CMD.py` for documentation)
    - Take any existing masks in sky position and proper motion
        - In the first pass, these should come directly from the galstreams track and the rough proper motion cut desscribed above. In later passes, it can be more specific.
    - Plot the masked data in CMD space, accounting for distance changes along the stream (using the distance track from galstreams)
    - Compare the theoretical isochrone to the data and make any small shifts in color or magnitude necessary to create a better match
        - In the first pass (before pm cuts), I do not make any shifts because the stream population is not distinguishable from the background. So I take the isochrone parameters to be the truth. In later passes, I allow shifts.
    - Fit a spline to this corrected theoretical isochrone and create the polygon by making upper and lower color bounds from that spline fit
3. Proper motion polygon in pm space
    - Take any existing masks in sky position and CMD (including an "offstream mask")
    - Make an initial guess at the mean proper motion (either an input or taken to be proper motion in middle of stream)
    - Use this guess to fit for the mean and standard deviation in each proper motion coordinate. This uses the residual pm histogram between the on-stream and off-stream regions.
    - Create polygon by fitting the best ellipse around the mean proper motion using a 2D Gaussian
    - Then, take galstreams proper motion tracks and "move" the ellipse along the stream to create a mask (includes extrapolation). This does not create a polygon but it can make a useful mask. See the next polygon for a way to make polygons that encodes this information
4. Proper motion polygon in each proper motion coordinate
    - Take any existing masks in sky position and CMD (including an "offstream mask")
    - Make an initial guess at the mean proper motion (either an input or taken to be proper motion in middle of stream)
    - Use this guess to fit for the mean and standard deviation in each proper motion coordinate. This uses the residual pm histogram between the on-stream and off-stream regions.
    - Take the galstreams tracks in each pm coordinate and create a polygon with half-width equal to the fitted standard deviation.
    - Limitation for now: does not extrapolate to outside the galstreams track so would not capture any potential stream stars outside this region.
    
Notes: 
1. There is currently an issue when trying to use g vs g-i rather than g vs g-r for Pal 5. Not sure where this is coming from, I'll try to look into it but would appreciate any ideas.
2. In general, I think it's probably possible to do better with Pal 5 but I was trying to do this in as non-specific a way as possible, so that it could be generalized.
3. I've not added footprints in the way that it was originally designed. I couldn't figure out the add_new_footprint functions especially for the CMD footprints. So I am just updating the pawprint by saying `pawprint.cmdprint = new_cmdprint` (as opposed to `pawprint.add_new_cmdprint(new_cmd_print,...)` as originally designed. This can probably be easily updated by someone who originally helped write the pawprint functions.
4. I have only updated parts of the package that are related to the output that I was attempting to generate.
5. I can definitely improve the legibility of the code to make it more understandable. I will work on that next but could use help here also from those who have written packages before.