#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 11:22:24 2022

@author: Ani, Kiyan, Richard
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.interpolate import interp1d
from scipy.signal import correlate2d
from ugali.analysis.isochrone import factory as isochrone_factory
from astropy.convolution import Gaussian1DKernel, convolve

plt.rc(
    "xtick",
    top=True,
    direction="in",
    labelsize=15,
)
plt.rc(
    "ytick",
    right=True,
    direction="in",
    labelsize=15,
)
plt.rc(
    "font",
    family="Arial",
)


class Isochrone:
    def __init__(
        self, cat, age, feh, distance, alpha, sky_poly, pm_poly, back_sky_poly
    ):

        """
        Defining variables loaded into class.
        Note- the assumption here is that the on-stream and off-stream sky
              polygon are the same area

        ------------------------------------------------------------------

        Parameters:
        cat = Input catalogue.
        age = Input age of stream from galstreams and/or literature.
        feh = Input metallicity from galstreams and/or literature.
        distance = Input distance from galstreams and/or literature.
        alpha = alpha/Fe
        sky_poly = Polygon mask from sky selection in (RA, Dec).
        pm_poly = Polygon mask from proper motions in (pm1, pm2).
        back_sky_poly = Polygon mask for background
        """
        # Maybe include distance track?
        # Pull survey from catalog?
        self.cat = cat
        self.age = age
        self.feh = feh
        self.distance = distance
        self.alpha = alpha
        self.sky_poly = sky_poly
        self.pm_poly = pm_poly
        self.x_shift = 0
        self.y_shift = 0
        self.back_sky_poly = back_sky_poly
        self.survey = "ps1"
        self.cmd_mask = 0
        self.cmd_poly = 0

    def sel_sky(self):

        """
        Initialising the on-sky and off-sky polygon masks to return only contained/bkg sources.
        """
        cat = self.cat
        on_poly_patch = mpl.patches.Polygon(
            self.sky_poly, facecolor="none", edgecolor="k", linewidth=2
        )
        on_points = np.vstack((cat["phi1"], cat["phi2"])).T
        on_mask = on_poly_patch.get_path().contains_points(on_points)

        self.on_skymask = on_mask

        off_poly_patch = mpl.patches.Polygon(
            self.back_sky_poly, facecolor="none", edgecolor="k", linewidth=2
        )
        off_points = np.vstack((cat["phi1"], cat["phi2"])).T
        off_mask = off_poly_patch.get_path().contains_points(off_points)

        self.off_skymask = off_mask

    def sel_pm(self):

        """
        Initialising the proper motions polygon mask to return only contained sources.
        """

        cat = self.cat
        on_poly_patch = mpl.patches.Polygon(
            self.pm_poly, facecolor="none", edgecolor="k", linewidth=2
        )
        on_points = np.vstack((cat["pm_phi1_cosphi2"], cat["pm_phi2"])).T
        on_mask = on_poly_patch.get_path().contains_points(on_points)

        self.on_pmmask = on_mask

    def generate_isochrone(self):

        """
        load an isochrone, LF model for a given metallicity, age, distance
        """

        # Convert feh to z
        Y_p = 0.245  # Primordial He abundance (WMAP, 2003)
        c = 1.54  # He enrichment ratio
        ZX_solar = 0.0229
        z = (1 - Y_p) / ((1 + c) + (1 / ZX_solar) * 10 ** (-self.feh))

        iso = isochrone_factory(
            "Dotter",
            survey=self.survey,
            age=self.age,
            distance_modulus=5 * np.log10(self.distance * 1000) - 5,
            z=z,
            band_1="g",
            band_2="r",
        )

        iso.afe = self.alpha

        initial_mass, mass_pdf, actual_mass, mag_1, mag_2 = iso.sample(mass_steps=4e2)
        mag_1 = mag_1 + iso.distance_modulus
        mag_2 = mag_2 + iso.distance_modulus

        # Excise the horizontal branch
        turn_idx = scipy.signal.argrelextrema(mag_1, np.less)[0][0]
        initial_mass = initial_mass[0:turn_idx]
        mass_pdf = mass_pdf[0:turn_idx]
        actual_mass = actual_mass[0:turn_idx]
        mag_1 = mag_1[0:turn_idx]
        mag_2 = mag_2[0:turn_idx]

        mmag_1 = interp1d(initial_mass, mag_1, fill_value="extrapolate")
        mmag_2 = interp1d(initial_mass, mag_2, fill_value="extrapolate")
        mmass_pdf = interp1d(initial_mass, mass_pdf, fill_value="extrapolate")

        self.iso = iso
        self.mag1 = mag_1
        self.mag2 = mag_2
        self.masses = actual_mass
        self.mass_pdf = mass_pdf

        self.mmag_1 = mmag_1
        self.mmag_2 = mmag_2
        self.mmass_pdf = mmass_pdf

        # return iso, initial_mass, mass_pdf, actual_mass, mag_1, mag_2, mmag_1, mmag_2, \
        #            mmass_pdf

    def data_cmd(self, xbin, ybin, xrange=[-0.5, 1.0], yrange=[15, 22]):

        """
        Empirical CMD generated from the input catalogue.

        ------------------------------------------------------------------

        Parameters:
        xbin: Bin size for color.
        ybin: Bin size for magnitudes.
        xrange: Set the range of color values. Default is [-0.5, 1.0].
        yrange: Set the range of magnitude values. Default is [15, 22].
        """
        tab = self.cat
        x_bins = np.arange(xrange[0], xrange[1], xbin)  # Used 0.03 for Jhelum
        y_bins = np.arange(yrange[0], yrange[1], ybin)  # Used 0.2 for Jhelum

        data, xedges, yedges = np.histogram2d(
            (tab["g0"] - tab["r0"])[self.on_pmmask & self.on_skymask],
            tab["g0"][self.on_pmmask & self.on_skymask],
            bins=[x_bins, y_bins],
            density=True,
        )

        self.x_edges = xedges
        self.y_edges = yedges
        self.x_bins = x_bins
        self.y_bins = y_bins
        self.CMD_data = data.T

        background_data, _, _ = np.histogram2d(
            (tab["g0"] - tab["r0"])[self.on_pmmask & self.off_skymask],
            tab["g0"][self.on_pmmask & self.off_skymask],
            bins=[x_bins, y_bins],
            density=True,
        )

        self.background_CMD_data = background_data.T

        # Is this used anywhere currently? Can add as a plot further down... this should also be divided by the relative area
        self.subtracted_CMD_data = (self.CMD_data) - (self.background_CMD_data)

    def make_poly(self, iso_model, tolerance, maxmag=26, minmag=14):
        """
        Generate the CMD polygon that depends on the errors in the photometric data
        In other words, the width of the selection should increase for fainter g-magnitude.

        --------------------------------------------------------------------

        Parameters:
        iso_model: spline function describing best fitting theoretical isochrone model
        phot_errors: an 2D array or a spline that describes the size of the photometric errors as a function of g-magnitude
        maxmag: faint limit of theoretical isochrone, should be deeper than all data
        minmag: bright limit of theoretical isochrone, either include just MS and subgiant branch or whole isochrone

        Returns:
        cmd_poly: Polygon vertices in CMD space.
        cmd_mask: Boolean mask in CMD sapce.
        """

        magband1, magband2 = "g0", "r0"
        # mag_vals = np.arange(minmag, maxmag, 0.01)

        ##Use the photometric error to determine the width of the isochrone (tolerance)
        self.errFn_default()  # to get self.gerrs
        gerrs = np.array(self.gerrs(self.mag1))

        tol = np.sqrt(
            2 * gerrs**2 + tolerance**2
        )  # 1-sigma tolerance on the error added in quadrature
        self.tol = tol

        col_low_vals = iso_model(self.mag1) - self.tol
        col_high_vals = iso_model(self.mag1) + self.tol

        # from here it is the same as make_poly
        cmd_poly = np.concatenate(
            [
                np.array([col_low_vals, self.mag1]).T,
                np.flip(np.array([col_high_vals, self.mag1]).T, axis=0),
            ]
        )

        # to create the boolean mask, don't know if we want this
        cmd_poly_patch = mpl.patches.Polygon(
            cmd_poly, facecolor="none", edgecolor="k", linewidth=2
        )
        cmd_points = np.vstack(
            (self.cat[magband1] - self.cat[magband2], self.cat[magband1])
        ).T
        cmd_mask = cmd_poly_patch.get_path().contains_points(cmd_points)

        return cmd_poly, cmd_mask

    def correct_isochrone(self):
        """
        Correlate the 2D histograms from the data and the
        theoretical isochrone to find the shift in color
        and magnitude necessary for the best match
        """

        signal, xedges, yedges = np.histogram2d(
            self.mag1 - self.mag2,
            self.mag1,
            bins=[self.x_edges, self.y_edges],
            weights=np.ones(len(self.mag1)),
        )

        signal_counts, xedges, yedges = np.histogram2d(
            self.mag1 - self.mag2, self.mag1, bins=[self.x_edges, self.y_edges]
        )
        signal = signal / signal_counts
        signal[np.isnan(signal)] = 0.0
        signal = signal.T

        ccor2d = correlate2d(self.CMD_data, signal)
        y, x = np.unravel_index(np.argmax(ccor2d), ccor2d.shape)
        self.x_shift = (x - len(ccor2d[0]) / 2.0) * (self.x_edges[1] - self.x_edges[0])
        self.y_shift = (y - len(ccor2d) / 2.0) * (self.y_edges[1] - self.y_edges[0])
        self.signal = signal

    def simpleSln(
        self,
        tolerance,
        maxmag,
        mass_thresh=0.80,
    ):
        """
        Select the stars that are within the CMD polygon cut
        --------------------------------
        Parameters:
        - mag1: g-magnitudes of all stars in catalog (or other band for different surveys)
        - mag2: i-magnitudes of all stars in catalog (or other band for different surveys)
        - masses: theoretical mass from isochrone for each star in the catalog
        - tolerance: half the "width" of the created CMD polygon
        - maxmag: faint limit of created CMD polygon, should be deeper than all data
        - mass_thresh: upper limit for the theoretical mass that dictates the bright limit of the
                       theoretical isochrone used for polygon
        - coloff: shift in color from theoretical isochrone to data
        - magoff: shift in magnitude from theoretical isochrone to data

        Returns:
        - cmd_poly: vertices of the CMD polygon cut
        - cmd_mask: bitmask of stars that pass the polygon cut
        - iso_model: the theoretical isochrone after shifts
        - iso_low: the "left" bound of the CMD polygon cut made from theoretical isochrone
        - iso_high: the "right" bound of the CMD polygon cut made from theoretical isochrone
        """

        coloff = self.x_shift
        magoff = self.y_shift
        ind = self.masses < mass_thresh
        mag1 = self.mag1[ind]
        mag2 = self.mag2[ind]

        iso_low = interp1d(
            mag1 + magoff, mag1 - mag2 + coloff - tolerance, fill_value="extrapolate"
        )
        iso_high = interp1d(
            mag1 + magoff, mag1 - mag2 + coloff + tolerance, fill_value="extrapolate"
        )
        iso_model = interp1d(
            mag1 + magoff, mag1 - mag2 + coloff, fill_value="extrapolate"
        )

        cmd_poly, cmd_mask = self.make_poly(
            iso_model,
            tolerance,
            maxmag=24,
            minmag=14,
        )

        self.cmd_poly = cmd_poly
        self.cmd_mask = cmd_mask

        return cmd_poly, cmd_mask, iso_model, iso_low, iso_high

    def plot_CMD(self, tolerance):
        """
        Plot the shifted isochrone over a 2D histogram of the polygon-selected
        data.

        Returns matplotlib Figure.
        """
        cat = self.cat[self.on_pmmask & self.on_skymask]
        mag1 = self.mag1
        color = mag1 - self.mag2 + self.x_shift
        mag1 = mag1 + self.y_shift
        mass_pdf = self.masses
        band1 = "g"
        band2 = "r"
        bins = (np.linspace(-0.5, 1.5, 128), np.linspace(10, 22.5, 128))

        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111)

        ax.hist2d(
            cat["g0"] - cat["r0"],
            cat["g0"],
            bins=bins,
            norm=mpl.colors.LogNorm(),
            zorder=5,
        )
        ax.plot(
            color,
            mag1,
            color="k",
            ls="--",
            zorder=10,
        )

        ax.plot(
            color - tolerance,
            mag1,
            color="b",
            ls="--",
            zorder=10,
        )
        ax.plot(
            color + tolerance,
            mag1,
            color="b",
            ls="--",
            zorder=10,
        )

        ax.set_xlabel(
            f"{band1}-{band2}",
            fontsize=20,
        )
        ax.set_ylabel(
            f"{band1}",
            fontsize=20,
        )

        ax.set_ylim(21.7, 15)
        ax.set_xlim(-0.5, 1.0)

        return fig

    def convolve_1d(self, line, mag_err):

        """
        1D Gaussian convolution.

        ------------------------------------------------------------------

        Parameters:
        probabilities:
        mag_err: Uncertainty in the magnitudes.

        """

        sigma = mag_err / (self.y_bins[1] - self.y_bins[0])  # error in pixel units
        kernel = Gaussian1DKernel(sigma)
        convolved = convolve(line, kernel)

        return convolved

    def convolve_errors(self, intr_err=0.1):

        """

        1D Gaussian convolution of the data with uncertainities.

        ------------------------------------------------------------------

        Parameters:
        intr_err: Free to set. Default is 0.1.

        """

        self.smoothed_signal = np.zeros((len(self.signal), len(self.signal[0])))
        for i in range(len(self.smoothed_signal)):
            self.smoothed_signal[i] = self.convolve_1d(
                self.signal[i],
                np.sqrt(
                    self.gerrs(self.x_bins[i]) ** 2
                    + self.rerrs(self.y_bins[i]) ** 2
                    + intr_err**2
                ),
            )

    def errFn_default(self):

        """
        Generate the errors given input magnitude, based on the input survey
        """

        if self.survey == "ps1":
            self.gerrs = lambda x: 0.00363355415 + np.exp((x - 23.9127145) / 1.09685211)
            self.rerrs = lambda x: 0.00363355415 + np.exp((x - 23.9127145) / 1.09685211)
            print("using panstarrs")
        elif self.survey == "des":
            self.gerrs = lambda x: 0.0010908679647672335 + np.exp(
                (x - 27.091072029215375) / 1.0904624484538419
            )
            self.rerrs = lambda x: 0.0010908679647672335 + np.exp(
                (x - 27.091072029215375) / 1.0904624484538419
            )
        else:
            self.gerrs = lambda x: 0.0
            self.rerrs = lambda x: 0.0

    def errFn_data(self):
        gerrs = np.zeros(len(self.y_bins))
        rerrs = np.zeros(len(self.x_bins))

        for i in range(len(self.y_bins)):
            gerrs[i] = np.nanmedian(
                self.cat["g0"][abs(self.cat["g0"] - self.y_bins[i]) < self.ybin / 2]
            )
            rerrs[i] = np.nanmedian(
                self.cat["r0"][abs(self.cat["g0"] - self.x_bins[i]) < self.xbin / 2]
            )

        gerrs = interp1d(self.y_bins, gerrs, fill_value="extrapolate")
        rerrs = interp1d(self.x_bins, rerrs, fill_value="extrapolate")

        self.gerrs = gerrs
        self.rerrs = rerrs

    def hb_mask(self, width=0.5):
        """
        Returns a polygon selection along the BHB.
        Takes a BHB branch from E. Vasiliev and S. Koposov in BP/RP and
        transforms to the photometric system

        Parameters
        ----------
        width : 'float'
            Required width of the horizontal branch filter

        Returns:
        -------
        polygon : 'matplotlib.patches.Polygon'
            Polygon for PatchCollection (e.g.)
        hb_2d_array : 'array_like, shape (N, 2)'
           Horizontal branch widen by added width
        original_track : 'array_like, shape (N, 2)'
            Horizontal branch path from E. Vasiliev
            originaly from S. Koposov.

        """

        mag1, mag2 = "g0", "r0"
        bp_rp = [
            -0.15,
            -0.075,
            0.0,
            0.075,
            0.15,
            0.225,
            0.3,
            0.375,
            0.45,
            0.525,
            0.6,
            0.675,
            0.75,
            0.825,
            0.9,
        ]

        abs_g_mag = np.array(
            [
                1.56922028,
                1.30633187,
                0.98436866,
                0.8272674,
                0.73556216,
                0.66450369,
                0.64182914,
                0.5948267,
                0.62143181,
                0.61441487,
                0.62462876,
                0.60931419,
                0.60433558,
                0.57833912,
                0.56579838,
            ]
        )

        original_track = np.vstack([bp_rp, abs_g_mag]).T

        abs_g_mag_bright, abs_g_mag_dimmer = (abs_g_mag - width).tolist(), (
            abs_g_mag + width
        ).tolist()

        bp_rp_poly = np.array(bp_rp + bp_rp[::-1])
        abs_g_mag_poly = (
            np.array(abs_g_mag_bright + abs_g_mag_dimmer[::-1])
            + 5 * np.log10(self.distance * 1000)
            - 5
        )

        hb_2d_array = np.vstack([bp_rp_poly, abs_g_mag_poly]).T

        hb_polygon = mpl.patches.Polygon(
            hb_2d_array, facecolor="none", edgecolor="k", linewidth=2
        )

        hb_points = np.vstack((self.cat[mag1] - self.cat[mag2], self.cat[mag1])).T
        hb_mask = hb_polygon.get_path().contains_points(hb_points)

        return hb_polygon, hb_mask
