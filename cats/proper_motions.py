#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Sophia, Nora, Nondh, Lina, Bruno, Kiyan
"""
import sys

import gala.coordinates as gc
import galstreams
import matplotlib as mpl
import matplotlib.pyplot as plt

# %%
# import packages
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import InterpolatedUnivariateSpline as IUS
from scipy.interpolate import UnivariateSpline as US
from scipy.spatial import ConvexHull

sys.path.append("../")
from cats.inputs import stream_inputs as inputs
from cats.pawprint.pawprint import Footprint2D, Pawprint


def rough_pm_poly(pawprint, data, buffer=2):
    """
    Will return a polygon with a rough cut in proper motion space.
    This aims to be ~100% complete with no thoughts about purity.
    The goal is to use this cut in conjunction with the cmd cut in order
    to see the stream as a clear overdensity in (phi_1, phi_2), which will
    allow membership probability modeling
    """
    stream_fr = pawprint.track.stream_frame
    track = pawprint.track.track.transform_to(stream_fr)
    # track_refl = gc.reflex_correct(track)

    # use the galstream proper motion track
    track_pm1_min = np.min(track.pm_phi1_cosphi2.value)
    track_pm1_max = np.max(track.pm_phi1_cosphi2.value)
    track_pm2_min = np.min(track.pm_phi2.value)
    track_pm2_max = np.max(track.pm_phi2.value)

    # make rectangular box around this region with an extra 2 mas/yr on each side
    rough_pm_poly = np.array(
        [
            [track_pm1_min - buffer, track_pm2_min - buffer],
            [track_pm1_min - buffer, track_pm2_max + buffer],
            [track_pm1_max + buffer, track_pm2_max + buffer],
            [track_pm1_max + buffer, track_pm2_min - buffer],
        ]
    )

    pawprint.pmprint = Footprint2D(rough_pm_poly, footprint_type="cartesian")

    pm_points = np.vstack((data["pm_phi1_cosphi2_unrefl"], data["pm_phi2_unrefl"])).T
    rough_pm_mask = pawprint.pmprint.inside_footprint(pm_points)

    return pawprint.pmprint, rough_pm_mask


class ProperMotionSelection:
    def __init__(
        self,
        stream,
        data,
        pawprint,
        #                  CMD_mask=True,
        #                  spatial_mask_on=True,
        #                  spatial_mask_off=True,
        #                  pm_phi1_grad = None, # think we should take this from pawprint, or at least make that the default
        #                  pm_phi2_grad = None,
        best_pm_phi1_mean=None,
        best_pm_phi2_mean=None,
        best_pm_phi1_std=None,
        best_pm_phi2_std=None,
        cutoff=0.95,
        n_dispersion_phi1=1,
        n_dispersion_phi2=1,
        refine_factor=100,
    ):
        """
        stream_obj: galstream object that contains stream's proper motion tracks
        data:
        :param: stream_obj: from galstreams so far #TODO: generalize
        :param: CMD_mask: Used before the PM
        :param: spatial_mask_on:
        :param: spatial_mask_off:
        :param: best_pm_phi1_mean: best initial guess for mean of pm_phi1
        :param: best_pm_phi2_mean: best initial guess for mean of pm_phi2
        :param: best_pm_phi1_std: best initial guess for pm_phi1 standard deviation
        :param: best_pm_phi2_std: best initial guess for pm_phi2 standard deviation
        :param: n_dispersion_phi1: float, default set to 1 standard deviation around phi_1
        :param: n_dispersion_phi2: float, default set to 1 standard deviation around phi_2
        :param: refine_factor: int, default set to 100, how smooth are the edges of the polygons
        :param: cutoff: float, in [0,1], cutoff on the height of the pdf to keep the stars that have a probability to belong to the 2D gaussian above the cutoff value
        """

        # stream_obj starting as galstream but then should be replaced by best values that we find
        self.stream = stream
        self.stream_obj = pawprint.track
        self.data = data
        self.pawprint = pawprint
        self.distance = inputs[stream]["distance"]
        self.dist_mod = 5 * np.log10(self.distance * 1000) - 5

        self.cutoff = cutoff

        assert (
            self.cutoff <= 1 and self.cutoff >= 0
        ), "the value of self.cutoff put in does not make sense! It has to be between 0 and 1"

        # Get tracks from galstreams with splines
        spline_phi2, spline_pm1, spline_pm2, spline_dist = self.from_galstreams()

        ####################################
        ## Make a rough proper motion cut ##
        ####################################
        self.rough_pm_poly, self.rough_pm_mask = self.rough_pm(buffer=2)

        ##################################
        ## Different proper motion cuts ##
        ##################################
        if self.stream == "GD-1":
            # distance spline will instead be a distance modulus spline
            self.dist_mod_correct = spline_dist(self.data["phi1"]) - self.dist_mod
        else:
            self.dist_mod_correct = (
                5 * np.log10(spline_dist(self.data["phi1"]) * 1000) - 5
            ) - self.dist_mod

        ## for GD-1 specifically:
        # distmod_spl = np.poly1d([2.41e-4, 2.421e-2, 15.001])
        # self.dist_mod_correct = distmod_spl(self.cat["phi1"]) - self.dist_mod

        # SHOULD THE CMD CUT ALSO MAKE AN OFFSTREAM MASK? MAY BE USEFUL TO MAKE CUTS FOR SOME STREAMS
        self.initial_masks()
        self.pm_phi1_cosphi2 = self.data["pm_phi1_cosphi2_unrefl"][self.mask]
        self.pm_phi2 = self.data["pm_phi2_unrefl"][self.mask]

        ####################################
        ## Choosing an initial best guess ##
        ####################################

        mid_phi1 = np.median(self.track.phi1.value)
        print(mid_phi1)

        if best_pm_phi1_mean == None:
            # TODO: generalize this later to percentile_values = [16, 50, 84]

            #             if self.stream == 'Fjorm-M68':
            #                 self.best_pm_phi1_mean = 1
            #                 self.best_pm_phi2_mean = 4
            #             else:
            self.best_pm_phi1_mean = spline_pm1(mid_phi1)
            self.best_pm_phi2_mean = spline_pm2(mid_phi1)

            self.best_pm_phi1_std = 1  # np.std(pm_phi1_cosphi2)
            self.best_pm_phi2_std = 1  # np.std(pm_phi2)

        else:
            self.best_pm_phi1_mean = best_pm_phi1_mean
            self.best_pm_phi2_mean = best_pm_phi2_mean

            self.best_pm_phi1_std = best_pm_phi1_std
            self.best_pm_phi2_std = best_pm_phi2_std

        ##############################################
        ## Fitting for the best central pm location ##
        ##############################################

        if stream == "Fjorm-M68":
            print("Skipping peak pm fitting")
        else:
            print("Fitting for peak pm location")
            # Uses default best_pm_phi_mean to generate improved ones
            peak_locations = self.find_peak_location(
                self.data, x_width=3.0, y_width=3.0, draw_histograms=True
            )
            print(
                "Post-fitting (pm1_mean, pm2_mean, pm1_std, pm2_std): {} \n".format(
                    peak_locations
                )
            )
        #         except:
        #             print('Skipping peak pm fitting')

        ################################################
        ## Ellipse-like proper motion cut in PM space ##
        ################################################
        print("Producing the polygon and mask")
        self.pm_poly, self.pm_mask = self.build_poly_and_mask(
            n_dispersion_phi1=n_dispersion_phi1,
            n_dispersion_phi2=n_dispersion_phi2,
            refine_factor=refine_factor,
        )

        self.mask = self.pm_mask & self.spatial_mask_on & self.CMD_mask

        self.pm_phi1_cosphi2 = data["pm_phi1_cosphi2_unrefl"][self.mask]
        self.pm_phi2 = data["pm_phi2_unrefl"][self.mask]

        # Plot the ellipse-like cut
        #         self.plot_pms_scatter(self.data, mask=True,
        #                               n_dispersion_phi1=n_dispersion_phi1,
        #                               n_dispersion_phi2=n_dispersion_phi2)
        #         self.plot_pm_hist(self.data, pms=[self.best_pm_phi1_mean, self.best_pm_phi2_mean])

        ######################################################
        ## PM cut in PM space using PM gradient information ##
        ######################################################

        # Instead of one ellipse for the whole stream,
        #  change the center of the ellipse depending on location on stream
        #  Disadvantage: Can't do this with one polygon --> not compatible with pawprint

        self.pm12_mask = self.build_mask(
            self.data, spline_pm1, spline_pm2, self.pm_poly
        )
        self.mask = self.pm12_mask & self.spatial_mask_on & self.CMD_mask

        # Plot the cut that accounts for PM variations
        self.plot_pms_scatter(
            self.data,
            mask=True,
            n_dispersion_phi1=n_dispersion_phi1,
            n_dispersion_phi2=n_dispersion_phi2,
        )
        self.plot_pm_hist(
            self.data, pms=[self.best_pm_phi1_mean, self.best_pm_phi2_mean]
        )

        #################################################
        ## PM cut in (phi1, pm1) and (phi1, pm2) space ##
        #################################################
        (
            self.pm1_poly,
            self.pm2_poly,
            self.pm1_mask,
            self.pm2_mask,
        ) = self.build_pm12_polys_and_masks()
        self.mask = self.pm1_mask & self.pm2_mask & self.spatial_mask_on & self.CMD_mask

        # Plot the cut in (phi1, pm1) and (phi1, pm2) space
        #         self.plot_pms_scatter(self.data, mask=True,
        #                               n_dispersion_phi1=n_dispersion_phi1,
        #                               n_dispersion_phi2=n_dispersion_phi2)
        # self.plot_pm_hist(self.data, pms=[self.best_pm_phi1_mean, self.best_pm_phi2_mean])

        return None

    def from_galstreams(self):
        stream_fr = self.stream_obj.stream_frame
        self.track = self.stream_obj.track.transform_to(stream_fr)
        # track_refl = gc.reflex_correct(track)
        # self.track_refl = track_refl

        self.galstream_phi1 = self.track.phi1.value
        self.galstream_phi2 = self.track.phi2.value
        self.galstream_pm_phi1_cosphi2 = self.track.pm_phi1_cosphi2.value
        self.galstream_pm_phi2 = self.track.pm_phi2.value
        self.galstream_dist = self.track.distance.value

        # Making splines to tracks: should add weights from errors???
        if self.stream == "Pal5":
            spline_phi2 = IUS(
                self.galstream_phi1[20:-100], self.galstream_phi2[20:-100]
            )
            spline_pm1 = IUS(
                self.galstream_phi1[20:-100], self.galstream_pm_phi1_cosphi2[20:-100]
            )
            spline_pm2 = IUS(
                self.galstream_phi1[20:-100], self.galstream_pm_phi2[20:-100]
            )
        else:
            spline_phi2 = IUS(self.galstream_phi1, self.galstream_phi2)
            spline_pm1 = IUS(self.galstream_phi1, self.galstream_pm_phi1_cosphi2)
            spline_pm2 = IUS(self.galstream_phi1, self.galstream_pm_phi2)

        #         spline_phi2 = US(self.galstream_phi1, self.galstream_phi2, k=3, s=len(self.galstream_phi1)/1000)
        #         spline_pm1 = US(self.galstream_phi1, self.galstream_pm_phi1_cosphi2, k=3, s=len(self.galstream_phi1)/1000)
        #         spline_pm2 = US(self.galstream_phi1, self.galstream_pm_phi2, k=3, s=len(self.galstream_phi1)/1000)

        if self.stream == "GD-1":
            spline_dist = np.poly1d(
                [2.41e-4, 2.421e-2, 15.001]
            )  # from Adrian and Kiyan's separate work
        else:
            spline_dist = IUS(self.galstream_phi1, self.galstream_dist)

        return spline_phi2, spline_pm1, spline_pm2, spline_dist

    def sel_sky(self):
        """
        Initialising the on-sky polygon mask to return only contained sources.
        """
        on_poly_patch = mpl.patches.Polygon(
            self.pawprint.skyprint["stream"].vertices[::100],
            facecolor="none",
            edgecolor="k",
            linewidth=2,
        )
        on_points = np.vstack((self.data["phi1"], self.data["phi2"])).T
        on_mask = on_poly_patch.get_path().contains_points(on_points)

        off_poly_patch = mpl.patches.Polygon(
            self.pawprint.skyprint["background"].vertices[::100],
            facecolor="none",
            edgecolor="k",
            linewidth=2,
        )
        off_points = np.vstack((self.data["phi1"], self.data["phi2"])).T
        off_mask = off_poly_patch.get_path().contains_points(off_points)

        return on_mask, off_mask

    def sel_cmd(self):
        """
        Initialising the proper motions polygon mask to return only contained sources.
        """

        mag = inputs[self.stream]["mag"]
        color1 = inputs[self.stream]["color1"]
        color2 = inputs[self.stream]["color2"]

        cmd_points = np.vstack(
            (
                self.data[color1] - self.data[color2],
                self.data[mag] - self.dist_mod_correct,
            )
        ).T
        cmd_mask = self.pawprint.cmdprint.inside_footprint(cmd_points)

        return cmd_mask

    def initial_masks(self):
        """
        Generate the initial spatial, and CMD masks based on the input
        """
        self.spatial_mask_on, self.spatial_mask_off = self.sel_sky()
        self.CMD_mask = self.sel_cmd()
        self.mask = self.spatial_mask_on & self.CMD_mask
        self.off_mask = self.spatial_mask_off & self.CMD_mask

    def rough_pm(self, buffer=2):
        """
        Will return a polygon with a rough cut in proper motion space.
        This aims to be ~100% complete with no thoughts about purity.
        The goal is to use this cut in conjunction with the cmd cut in order
        to see the stream as a clear overdensity in (phi_1, phi_2), which will
        allow membership probability modeling
        """
        # use the galstream proper motion track
        track_pm1_min = np.min(self.galstream_pm_phi1_cosphi2)
        track_pm1_max = np.max(self.galstream_pm_phi1_cosphi2)
        track_pm2_min = np.min(self.galstream_pm_phi2)
        track_pm2_max = np.max(self.galstream_pm_phi2)

        # make rectangular box around this region with an extra 2 mas/yr on each side
        self.rough_pm_poly = np.array(
            [
                [track_pm1_min - buffer, track_pm2_min - buffer],
                [track_pm1_min - buffer, track_pm2_max + buffer],
                [track_pm1_max + buffer, track_pm2_max + buffer],
                [track_pm1_max + buffer, track_pm2_min - buffer],
            ]
        )

        pm_points = np.vstack(
            (self.data["pm_phi1_cosphi2_unrefl"], self.data["pm_phi2_unrefl"])
        ).T
        self.pawprint.pmprint = Footprint2D(
            self.rough_pm_poly, footprint_type="cartesian"
        )
        self.rough_pm_mask = self.pawprint.pmprint.inside_footprint(pm_points)

        return self.rough_pm_poly, self.rough_pm_mask

    @staticmethod
    def two_dimensional_gaussian(x, y, x0, y0, sigma_x, sigma_y):
        """
        Evaluates a two dimensional gaussian distribution in x, y, with means x0, y0, and dispersions sigma_x and sigma_y
        """

        return np.exp(
            -((x - x0) ** 2 / (2 * sigma_x**2) + (y - y0) ** 2 / (2 * sigma_y**2))
        )

    def build_poly_and_mask(
        self, n_dispersion_phi1=3, n_dispersion_phi2=3, refine_factor=100
    ):
        """
        Builds the mask of the proper motion with n_dispersion around the mean
        :param: n_dispersion_phi1: float, default set to 1 standard deviation around phi_1
        :param: n_dispersion_phi2: float, default set to 1 standard deviation around phi_2
        :param: refine_factor: int, default set to 100, how smooth are the edges of the polygons
        :param: cutoff: float, in [0,1], cutoff on the height of the pdf to keep the stars that have a probability to belong to the 2D gaussian above the cutoff value

        :output: is a list of points that are the vertices of a polygon
        """

        # First generate the 2D histograms
        pm_phi1_min, pm_phi1_max = (
            self.best_pm_phi1_mean - n_dispersion_phi1 * self.best_pm_phi1_std,
            self.best_pm_phi1_mean + n_dispersion_phi1 * self.best_pm_phi1_std,
        )

        pm_phi2_min, pm_phi2_max = (
            self.best_pm_phi2_mean - n_dispersion_phi2 * self.best_pm_phi2_std,
            self.best_pm_phi2_mean + n_dispersion_phi2 * self.best_pm_phi2_std,
        )

        pm_phi1_array = np.linspace(pm_phi1_min, pm_phi1_max, refine_factor)
        pm_phi2_array = np.linspace(pm_phi2_min, pm_phi2_max, refine_factor)

        pm_pdf = np.zeros((len(pm_phi1_array), refine_factor))
        points_x = np.zeros((len(pm_phi1_array), refine_factor))
        points_y = np.zeros((len(pm_phi1_array), refine_factor))

        for n, s in enumerate(pm_phi1_array):
            pm_pdf[n, :] = self.two_dimensional_gaussian(
                s,
                pm_phi2_array,
                self.best_pm_phi1_mean,
                self.best_pm_phi2_mean,
                self.best_pm_phi1_std,
                self.best_pm_phi2_std,
            )
            points_x[n, :] = np.array([s for _ in range(len(pm_phi2_array))])
            points_y[n, :] = pm_phi2_array

        cut = np.where(pm_pdf.flatten() > self.cutoff)[0]
        x_cut = points_x.flatten()[cut]
        y_cut = points_y.flatten()[cut]

        xy = np.transpose([x_cut, y_cut])
        hull = ConvexHull(xy)

        # self.vertices = xy[hull.vertices]
        self.pm_poly = xy[hull.vertices]
        pm_points = np.vstack(
            (self.data["pm_phi1_cosphi2_unrefl"], self.data["pm_phi2_unrefl"])
        ).T
        self.pawprint.pmprint = Footprint2D(self.pm_poly, footprint_type="cartesian")
        self.pm_mask = self.pawprint.pmprint.inside_footprint(pm_points)

        return self.pm_poly, self.pm_mask

    def build_pm12_polys_and_masks(self):
        """
        This assumes that galstreams is correct, which is maybe not a great assumption but will work for now.
        """
        self.pm1_poly = np.concatenate(
            [
                np.array(
                    [
                        self.galstream_phi1[::50],
                        self.galstream_pm_phi1_cosphi2[::50] - self.best_pm_phi1_std,
                    ]
                ).T,
                np.flip(
                    np.array(
                        [
                            self.galstream_phi1[::50],
                            self.galstream_pm_phi1_cosphi2[::50]
                            + self.best_pm_phi1_std,
                        ]
                    ).T,
                    axis=0,
                ),
            ]
        )

        self.pm2_poly = np.concatenate(
            [
                np.array(
                    [
                        self.galstream_phi1[::50],
                        self.galstream_pm_phi2[::50] - self.best_pm_phi2_std,
                    ]
                ).T,
                np.flip(
                    np.array(
                        [
                            self.galstream_phi1[::50],
                            self.galstream_pm_phi2[::50] + self.best_pm_phi2_std,
                        ]
                    ).T,
                    axis=0,
                ),
            ]
        )

        pm1_points = np.vstack(
            (self.data["phi1"], self.data["pm_phi1_cosphi2_unrefl"])
        ).T
        pm2_points = np.vstack((self.data["phi1"], self.data["pm_phi2_unrefl"])).T
        self.pawprint.pm1print = Footprint2D(self.pm1_poly, footprint_type="cartesian")
        self.pm1_mask = self.pawprint.pm1print.inside_footprint(pm1_points)
        self.pawprint.pm2print = Footprint2D(self.pm2_poly, footprint_type="cartesian")
        self.pm2_mask = self.pawprint.pm2print.inside_footprint(pm2_points)

        return self.pm1_poly, self.pm2_poly, self.pm1_mask, self.pm2_mask

    def build_mask(self, data, spline_pm1, spline_pm2, pm_poly):
        """
        This builds a mask (i.e. finds the data points satisfying pm constraints)
        that does not use the peak fitting used elsewhere.
        It relies on splines for pm_phi1_cosphi2 and pm_phi2 vs phi1 which must be given as inputs
        Most of the time, these will naturally come from galstreams
        """

        pm1_data_corrected = data["pm_phi1_cosphi2_unrefl"] - spline_pm1(data["phi1"])
        pm2_data_corrected = data["pm_phi2_unrefl"] - spline_pm2(data["phi1"])

        pm_vert_corrected = pm_poly - [self.best_pm_phi1_mean, self.best_pm_phi2_mean]

        pm_corrected_poly_patch = mpl.patches.Polygon(
            pm_vert_corrected, facecolor="none", edgecolor="k", linewidth=2
        )

        pm_points = np.vstack((pm1_data_corrected, pm2_data_corrected)).T
        pm_mask = pm_corrected_poly_patch.get_path().contains_points(pm_points)
        # self.pawprint.pmprint = Footprint2D(pm_vert_corrected, footprint_type='cartesian')
        # self.pm_mask = self.pawprint.pmprint.inside_footprint(pm_points)

        return pm_mask

    def plot_pms_scatter(
        self,
        data,
        save=True,
        mask=False,
        n_dispersion_phi1=1,
        n_dispersion_phi2=1,
        refine_factor=100,
        **kwargs,
    ):
        """
        Plot proper motions on stream and off stream scatter or hist2d plots
        :param: save: boolean, whether or not to save the figure
        :param: mask: boolean, if true, calls in the mask
        """
        data_on = data[self.mask]
        data_off = data[self.off_mask]

        fig, ax = plt.subplots(1, 2, figsize=(10, 5))

        # resize and fix column name
        scatter_size = 1.0 / data_on["pmra_error"]

        ax[0].scatter(
            data_on["pm_phi1_cosphi2_unrefl"],
            data_on["pm_phi2_unrefl"],
            c="k",
            s=scatter_size,
            alpha=0.2,
            **kwargs,
        )

        if mask:
            vertices, _ = self.build_poly_and_mask(
                n_dispersion_phi1=n_dispersion_phi1,
                n_dispersion_phi2=n_dispersion_phi2,
                refine_factor=refine_factor,
            )

            x, y = vertices.T
            x = np.append(x, x[0])
            y = np.append(y, y[0])
            vert = np.transpose([x, y])

            ax[0].plot(vert.T[0], vert.T[1], "k-")

        ax[0].set_xlim(-20, 20)
        ax[0].set_ylim(-20, 20)
        ax[0].set_xlabel("$\mu_{\phi_1}$ [mas yr$^{-1}$]")
        ax[0].set_ylabel("$\mu_{\phi_2}$ [mas yr$^{-1}$]")
        ax[0].set_title("Stream", fontsize="medium")

        # resize and fix column name
        scatter_size = 1.0 / data_off["pmra_error"]
        ax[1].scatter(
            data_off["pm_phi1_cosphi2_unrefl"],
            data_off["pm_phi2_unrefl"],
            c="k",
            s=scatter_size,
            alpha=0.2,
            **kwargs,
        )

        if mask:
            ax[1].plot(vert.T[0], vert.T[1], "k-")

        ax[1].set_xlim(-20, 20)
        ax[1].set_ylim(-20, 20)
        ax[1].set_xlabel("$\mu_{\phi_1}$ [mas yr$^{-1}$]")
        ax[1].set_ylabel("$\mu_{\phi_2}$ [mas yr$^{-1}$]")
        ax[1].set_title("Off stream", fontsize="medium")

        fig.tight_layout()

        if save:
            fig.savefig(
                f"proper_motions_{self.stream_obj.stream_name}_scatter.png",
                bbox_inches="tight",
            )
            fig.savefig(
                f"proper_motions_{self.stream_obj.stream_name}_scatter.pdf",
                bbox_inches="tight",
            )

        return fig

    def plot_pm_hist(
        self,
        data,
        dx=0.5,
        norm=1,
        save=0,
        pms=[None, None],
        match_norm=False,
        stream_coords=True,
        reflex_corr=True,
        zero_line=True,
        pm_lims=[-20, 20],
        **kwargs,
    ):
        # Code from Nora

        data_on = data[self.mask]
        data_off = data[self.off_mask]

        bins = np.arange(pm_lims[0], pm_lims[1] + dx / 2.0, dx)

        # data access depends on structure, needs to be adapted
        if stream_coords:
            h1 = np.histogram2d(
                data_on["pm_phi1_cosphi2_unrefl"], data_on["pm_phi2_unrefl"], bins
            )[0]
            h2 = np.histogram2d(
                data_off["pm_phi1_cosphi2_unrefl"], data_off["pm_phi2_unrefl"], bins
            )[0]
        else:
            if reflex_corr:
                h1 = np.histogram2d(data_on["pm_ra"], data_on["pm_dec"], bins)[0]
                h2 = np.histogram2d(data_off["pm_ra"], data_off["pm_dec"], bins)[0]
            else:
                h1 = np.histogram2d(data_on["PMRA0"], data_on["PMDEC0"], bins)[0]
                h2 = np.histogram2d(data_off["PMRA0"], data_off["PMDEC0"], bins)[0]

        # might need to normalise histogram for different areas of off stream mask for subtraction histogram

        h2 *= norm
        # print h1.sum(), h2.sum()
        diff = h1 - h2

        if match_norm:
            vmin, vmax = np.min([np.min(h1), np.min(h2), np.min(diff)]), np.max(
                [np.max(h1), np.max(h2), np.max(diff)]
            )
        else:
            vmin, vmax = None, None
        vmin = 0.0

        offset = -dx / 2.0

        fig, ((ax1, ax2, ax3)) = plt.subplots(
            1, 3, sharex="col", sharey="row", figsize=(15, 6)
        )
        im1 = ax1.imshow(
            h1.T,
            extent=[
                bins.min() - offset,
                bins.max() - offset,
                bins.min() - offset,
                bins.max() - offset,
            ],
            origin="lower",
            vmin=vmin,
            vmax=vmax,
            interpolation="none",
            **kwargs,
        )
        im2 = ax2.imshow(
            h2.T,
            extent=[
                bins.min() - offset,
                bins.max() - offset,
                bins.min() - offset,
                bins.max() - offset,
            ],
            origin="lower",
            vmin=vmin,
            vmax=vmax,
            interpolation="none",
            **kwargs,
        )
        im3 = ax3.imshow(
            diff.T,
            extent=[
                bins.min() - offset,
                bins.max() - offset,
                bins.min() - offset,
                bins.max() - offset,
            ],
            origin="lower",
            vmin=vmin,
            vmax=vmax,
            interpolation="none",
            **kwargs,
        )

        colorbar(im1)
        colorbar(im2)
        colorbar(im3)

        if (pms[0] is None) or (pms[1] is None):
            ax1.axvline(self.best_pm[0], ls="--", c="k", lw=1)
            ax2.axvline(self.best_pm[0], ls="--", c="k", lw=1)
            ax3.axvline(self.best_pm[0], ls="--", c="k", lw=1)
            ax1.axhline(self.best_pm[1], ls="--", c="k", lw=1)
            ax2.axhline(self.best_pm[1], ls="--", c="k", lw=1)
            ax3.axhline(self.best_pm[1], ls="--", c="k", lw=1)
        else:
            ax1.axvline(pms[0], ls="--", c="k", lw=1)
            ax2.axvline(pms[0], ls="--", c="k", lw=1)
            ax3.axvline(pms[0], ls="--", c="k", lw=1)
            ax1.axhline(pms[1], ls="--", c="k", lw=1)
            ax2.axhline(pms[1], ls="--", c="k", lw=1)
            ax3.axhline(pms[1], ls="--", c="k", lw=1)
        # f.suptitle(r'$\mathrm{%s}$' % stream['name'], fontsize=20)

        ax1.set_title(r"$\mathrm{on-stream}$")
        ax2.set_title(r"$\mathrm{off-stream}$")
        ax3.set_title(r"$\mathrm{residual}$")

        if stream_coords:
            ax1.set_xlabel(r"$\mu_1$")
            ax2.set_xlabel(r"$\mu_1$")
            ax3.set_xlabel(r"$\mu_1$")

            ax1.set_ylabel(r"$\mu_2$")
        else:
            if reflex_corr:
                ax1.set_xlabel(r"$\mu_{ra}$")
                ax2.set_xlabel(r"$\mu_{ra}$")
                ax3.set_xlabel(r"$\mu_{ra}$")

                ax1.set_ylabel(r"$\mu_{dec}$")

        if zero_line:
            ax1.axhline(0, ls="--", c="k", lw=0.5)
            ax2.axhline(0, ls="--", c="k", lw=0.5)
            ax3.axhline(0, ls="--", c="k", lw=0.5)

        ax1.set_xlim(bins.min(), bins.max())
        ax2.set_xlim(bins.min(), bins.max())
        ax3.set_xlim(bins.min(), bins.max())
        ax1.set_ylim(bins.min(), bins.max())
        ax2.set_ylim(bins.min(), bins.max())
        ax3.set_ylim(bins.min(), bins.max())

        fig.tight_layout()

        if save:
            fig.savefig(
                f"proper_motions_{self.stream_obj.stream_name}_hist.png",
                bbox_inches="tight",
            )
            fig.savefig(
                f"proper_motions_{self.stream_obj.stream_name}_hist.pdf",
                bbox_inches="tight",
            )
            # plt.close()
        return fig

        # ========================= added Nov 3 (need checking) ========================

    def find_peak_location(self, data, x_width=3.0, y_width=3.0, draw_histograms=True):
        """
        find peak location in the proper motion space
        :param: data: list of the stellar parameters to get the peak pm.
        :param: x_width: float, half x-size of zoomed region box, default set to 3.
        :param: y_width: float, half y-size of zoomed region box, default set to 3.
        :param: draw_histograms: print histograms, default set to True

        output: [pm_x_cen, pm_y_cen, x_std, y_std]: array
            pm_x_cen: peak proper motion in phi1
            pm_y_cen: peak proper motion in phi2
            x_std: standard deviation proper motion in phi1
            y_std: standard deviation proper motion in phi2
        """
        from astropy.modeling import fitting, models
        from matplotlib.colors import LogNorm

        x_center, y_center = self.best_pm_phi1_mean, self.best_pm_phi2_mean
        print("Pre-fitting mean PM values: {}, {}".format(x_center, y_center))
        xmin, xmax, ymin, ymax = (
            x_center - x_width,
            x_center + x_width,
            y_center - y_width,
            y_center + y_width,
        )

        # on-stream data (after masking) into 2D histogram
        H1, x_edges, y_edges = np.histogram2d(
            self.pm_phi1_cosphi2,
            self.pm_phi2,
            bins=(
                np.arange(xmin, xmax, x_width / 30),
                np.arange(ymin, ymax, y_width / 30),
            ),
        )  # explore differtent bin sizes

        stream_off = self.off_mask
        # off-stream data
        H2, x_edges, y_edges = np.histogram2d(
            data["pm_phi1_cosphi2_unrefl"][stream_off],
            data["pm_phi2_unrefl"][stream_off],
            bins=(
                np.arange(xmin, xmax, x_width / 30),
                np.arange(ymin, ymax, y_width / 30),
            ),
        )

        hist = (
            H1 - self.mask.sum() / self.off_mask.sum() * H2
        )  # check this scale factor --> self.mask.sum()/self.off_mask.sum()
        # Do we want to do based on counts or based on area, since we do expect more counts on stream (but maybe negligible)

        # fitting 2D gaussian (Code from Ani)
        # Find overdensity
        ind = np.unravel_index(hist.argmax(), hist.shape)
        # Fit 2d gaussian to a zoomed in histogram

        # hist_zoom = np.transpose(hist[(ind[0]-6):(ind[0]+7), (ind[1]-6):(ind[1]+7)])
        hist_zoom = np.transpose(
            hist[(ind[0] - 10) : (ind[0] + 11), (ind[1] - 10) : (ind[1] + 11)]
        )
        g_init = models.Gaussian2D(
            amplitude=1.0,
            x_mean=x_edges[ind[0]],
            y_mean=y_edges[ind[1]],
            x_stddev=0.5,
            y_stddev=0.5,
        )
        fit_g = fitting.LevMarLSQFitter()
        # x,y = np.meshgrid(x_edges[(ind[0]-6):(ind[0]+7)], y_edges[(ind[1]-6):(ind[1]+7)])
        x, y = np.meshgrid(x_edges[:-1], y_edges[:-1])

        # g = fit_g(g_init, x, y, hist_zoom)
        g = fit_g(g_init, x, y, hist)

        # the mean position of the peak
        pm_x_cen = g.x_mean + 0.5 * (x_edges[1] - x_edges[0])  # half width of the bin
        pm_y_cen = g.y_mean + 0.5 * (y_edges[1] - y_edges[0])
        # the std's of the peak
        x_std, y_std = g.x_stddev.value, g.y_stddev.value

        # draw proper motions of the on-stream, off-stream and residual histogram
        if draw_histograms:
            fig, axes = plt.subplots(
                1, 3, figsize=(15, 5), sharex=True, sharey=True, constrained_layout=True
            )

            bins = np.linspace(-20, 20, 128)
            H1, xe, ye, _ = axes[0].hist2d(
                self.pm_phi1_cosphi2, self.pm_phi2, bins=bins, norm=LogNorm()
            )
            H2, *_ = axes[1].hist2d(
                data["pm_phi1_cosphi2_unrefl"][stream_off],
                data["pm_phi2_unrefl"][stream_off],
                bins=bins,
                norm=LogNorm(),
            )

            axes[2].pcolormesh(
                xe,
                ye,
                (H1 - self.mask.sum() / self.off_mask.sum() * H2).T,
                cmap="RdBu",
                vmin=-10,
                vmax=10,
            )
            t = np.linspace(0, 2 * np.pi, 100)
            axes[2].plot(
                pm_x_cen + x_std * np.cos(t), pm_y_cen + y_std * np.sin(t), c="green"
            )
            # axes[2].set_xlim(xmin,xmax)
            # axes[2].set_ylim(ymin,ymax)

            for ax in axes[:2]:
                ax.plot(
                    self.galstream_pm_phi1_cosphi2, self.galstream_pm_phi2, color="cyan"
                )

                ax.set_xlabel(r"$\mu_{\phi_1}$", fontsize=20)
            axes[0].set_ylabel(r"$\mu_{\phi_2}$", fontsize=20)

        self.best_pm_phi1_mean = pm_x_cen
        self.best_pm_phi2_mean = pm_y_cen
        self.best_pm_phi1_std = x_std
        self.best_pm_phi2_std = y_std

        return [pm_x_cen, pm_y_cen, x_std, y_std]


def colorbar(mappable):
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    return fig.colorbar(mappable, cax=cax)
