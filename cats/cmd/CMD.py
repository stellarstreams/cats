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
    def __init__(self, cat, age, feh, distance, alpha, sky_poly, pm_poly):
        """
        Defining variables loaded into class.

        ------------------------------------------------------------------

        Parameters:
        cat = Input catalogue.
        age = Input age of stream from galstreams and/or literature.
        feh = Input metallicity from galstreams and/or literature.
        distance = Input distance from galstreams and/or literature.
        alpha = alpha/Fe
        sky_poly = Polygon mask from sky selection in (RA, Dec).
        pm_poly = Polygon mask from proper motions in (pm1, pm2).
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
        self.survey = "ps1"

    def sel_sky(self):
        """
        Initialising the on-sky polygon mask to return only contained sources.
        """
        cat = self.cat
        on_poly_patch = mpl.patches.Polygon(
            self.sky_poly, facecolor="none", edgecolor="k", linewidth=2
        )
        on_points = np.vstack((cat["phi1"], cat["phi2"])).T
        on_mask = on_poly_patch.get_path().contains_points(on_points)

        self.on_skymask = on_mask

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
        self.CMD_data = data.T

    def make_poly(self, iso_low, iso_high, maxmag=26, minmag=14):
        """
        Generate the CMD polygon mask.

        ------------------------------------------------------------------

        Parameters:
        iso_low: spline function describing the "left" bound of the theorietical isochrone
        iso_high: spline function describing the "right" bound of the theoretical isochrone
        maxmag: faint limit of theoretical isochrone, should be deeper than all data
        minmag: bright limit of theoretical isochrone, either include just MS and subgiant branch or whole isochrone

        Returns:
        cmd_poly: Polygon vertices in CMD space.
        cmd_mask: Boolean mask in CMD sapce.

        """

        mag1, mag2 = "g0", "r0"
        mag_vals = np.arange(minmag, maxmag, 0.01)
        col_low_vals = iso_low(mag_vals)
        col_high_vals = iso_high(mag_vals)

        cmd_poly = np.concatenate(
            [
                np.array([col_low_vals, mag_vals]).T,
                np.flip(np.array([col_high_vals, mag_vals]).T, axis=0),
            ]
        )

        # to create the boolean mask, don't know if we want this
        cmd_poly_patch = mpl.patches.Polygon(
            cmd_poly, facecolor="none", edgecolor="k", linewidth=2
        )
        cmd_points = np.vstack((self.cat[mag1] - self.cat[mag2], self.cat[mag1])).T
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

        cmd_poly, cmd_mask = self.make_poly(iso_low, iso_high, maxmag=24, minmag=14)

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

    def convolve_1d(probabilities, mag_err):
        """
        1D Gaussian convolution.

        ------------------------------------------------------------------

        Parameters:
        probabilities:
        mag_err: Uncertainty in the magnitudes.

        """
        self.probabilities = probabilities
        self.mag_err = mag_err

        sigma = mag_err / self.ybin  # error in pixel units
        kernel = Gaussian1DKernel(sigma)
        convolved = convolve(probabilities, kernel)

        convolved.self = convolved

    def convolve_errors(self, g_errors, r_errors, intr_err=0.1):
        """

        1D Gaussian convolution of the data with uncertainities.

        ------------------------------------------------------------------

        Parameters:
        g_errors: g magnitude uncertainties.
        r_errors: r magnitude uncertainties.
        intr_err: Free to set. Default is 0.1.

        """

        for i in range(len(probabilities)):
            probabilities[i] = convolve_1d(
                probabilities[i],
                np.sqrt(
                    g_errors(self.x_bins[i]) ** 2
                    + r_errors(self.y_bins[i]) ** 2
                    + intr_err**2
                ),
                sel.fx_bins[1] - self.x_bins[0],
            )

        self.probabilities = probabilities

    def errFn(self):
        """
        Generate the errors for the magnitudes?
        """

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
