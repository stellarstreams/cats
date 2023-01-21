#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Sophia, Nora, Nondh, Lina, Bruno
"""
#%%
# import packages
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import galstreams
import gala.coordinates as gc
from scipy.spatial import ConvexHull


# %%
# first CMD selection
## default CMD mask fixed but for distance
# take Gaia data crossmatched with photometric data
# apply CMD and spatial masks
# plot on stream and off stream proper motion 
# automated to find highest density location
# define our mask


class ProperMotionSelection():
    def __init__(self, stream_obj, 
                 CMD_mask=None, 
                 spatial_mask_on=None, 
                 spatial_mask_off=None, 
                 best_pm_phi1_mean = None,
                 best_pm_phi2_mean = None, 
                 best_pm_phi1_std = None,
                 best_pm_phi2_std = None,
                 cutoff = 0.95,
                 n_dispersion_phi1=1, 
                 n_dispersion_phi2=1, 
                 refine_factor = 100):
        '''
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
        '''
        # stream_obj starting as galstream but then should be replaced by best values that we find
        
        self.stream_obj = stream_obj
        
        self.spatial_mask_on = spatial_mask_on
        self.spatial_mask_off = spatial_mask_off
        self.mask = spatial_mask_on & CMD_mask
        self.off_mask = spatial_mask_off & CMD_mask
        self.cutoff = cutoff

        assert(self.cutoff <= 1 and self.cutoff >= 0), "the value of self.cutoff put in does not make sense! It has to be between 0 and 1"

        stream_fr = stream_obj.stream_frame
        track = stream_obj.track.transform_to(stream_fr)
        track_refl = gc.reflex_correct(track)
        
        pm_phi1_cosphi2 = track_refl.pm_phi1_cosphi2.value
        pm_phi2 = track_refl.pm_phi2.value
        
        if best_pm_phi1_mean == None:
            # TODO: generalize this later to percentile_values = [16, 50, 84]

            self.best_pm_phi1_mean = np.mean(pm_phi1_cosphi2)
            self.best_pm_phi2_mean = np.mean(pm_phi2)

            self.pm_phi1_std = np.std(pm_phi1_cosphi2)
            self.pm_phi2_std = np.std(pm_phi2)

        else:
            self.best_pm_phi1_mean = best_pm_phi1_mean
            self.best_pm_phi2_mean = best_pm_phi2_mean

            self.best_pm_phi1_std = best_pm_phi1_std
            self.best_pm_phi2_std = best_pm_phi2_std

        print("Producing the initial mask")
        self.pm_mask = self.build_mask(n_dispersion_phi1=n_dispersion_phi1, n_dispersion_phi2=n_dispersion_phi2, refine_factor = refine_factor)
        
        self.mask = self.mask & self.pm_mask

        # Nondh and Bruno added
        # from galstream
        self.galstream_pm_phi1_cosphi2 = track_refl.pm_phi1_cosphi2.value
        self.galstream_pm_phi2 = track_refl.pm_phi2.value
        self.initial_pm_phi1_cosphi2_center = np.median(self.galstream_pm_phi1_cosphi2)
        self.initial_pm_phi2_center = np.median(self.galstream_pm_phi2)
        
        self.pm_phi1_cosphi2 = data['pm_phi1_cosphi2'][self.mask]
        self.pm_phi2 = data['pm_phi2'][self.mask]

        return None
    
    @staticmethod
    def two_dimensional_gaussian(x, y, x0, y0, sigma_x, sigma_y):
        """
        Evaluates a two dimensional gaussian distribution in x, y, with means x0, y0, and dispersions sigma_x and sigma_y
        """

        return np.exp(- ( (x-x0)**2/(2*sigma_x**2) + (y-y0)**2/(2*sigma_y**2) ) )

    def build_mask(self, data, n_dispersion_phi1=1, n_dispersion_phi2=1, refine_factor = 100):
        """
        Builds the mask of the proper motion with n_dispersion around the mean
        :param: data, the initial data file to get the peak of the pm
        :param: n_dispersion_phi1: float, default set to 1 standard deviation around phi_1
        :param: n_dispersion_phi2: float, default set to 1 standard deviation around phi_2
        :param: refine_factor: int, default set to 100, how smooth are the edges of the polygons
        :param: cutoff: float, in [0,1], cutoff on the height of the pdf to keep the stars that have a probability to belong to the 2D gaussian above the cutoff value

        :output: is a list of points that are the vertices of a polygon
        """

        # First generate the 2D histograms
        pm_phi1_min, pm_phi1_max = (self.best_pm_phi1_mean - \
            n_dispersion_phi1*self.pm_phi1_std, \
            self.best_pm_phi1_mean + n_dispersion_phi1*self.pm_phi1_std)

        pm_phi2_min, pm_phi2_max = (self.best_pm_phi2_mean - \
            n_dispersion_phi2*self.pm_phi2_std, \
            self.best_pm_phi2_mean + n_dispersion_phi2*self.pm_phi2_std)

        pm_phi1_array = np.linspace(pm_phi1_min, pm_phi1_max, refine_factor)
        pm_phi2_array = np.linspace(pm_phi2_min, pm_phi2_max, refine_factor)

        pm_pdf = np.zeros((len(pm_phi1_array), refine_factor))
        points_x = np.zeros((len(pm_phi1_array), refine_factor))
        points_y = np.zeros((len(pm_phi1_array), refine_factor))

        for n, s in enumerate(pm_phi1_array):
            pm_pdf[n,:] = self.two_dimensional_gaussian(s, pm_phi2_array, self.best_pm_phi1_mean, self.best_pm_phi2_mean, self.pm_phi1_std, self.pm_phi2_std)
            points_x[n,:] = np.array([s for _ in range(len(pm_phi2_array))])
            points_y[n,:] = pm_phi2_array


        cut = np.where(pm_pdf.flatten()> self.cutoff)[0]
        x_cut = points_x.flatten()[cut]
        y_cut = points_y.flatten()[cut]

        xy = np.transpose([x_cut, y_cut])
        hull = ConvexHull(xy)

        # self.vertices = xy[hull.vertices]
        return xy[hull.vertices]

    def plot_pms_scatter(self, data, save=True, mask=False, n_dispersion_phi1=1, n_dispersion_phi2=1, refine_factor = 100, **kwargs):
        '''
        Plot proper motions on stream and off stream scatter or hist2d plots
        :param: save: boolean, whether or not to save the figure
        :param: mask: boolean, if true, calls in the mask 
        '''
        data_on = data[self.mask]
        data_off = data[self.off_mask]
                
        fig, ax = plt.subplots(1,2, figsize=(10,5))

        # resize and fix column name
        scatter_size = 1./data_on['pmra_error']
        
        ax[0].scatter(data_on['PMPHI1'], data_on['PMPHI2'], c='k', s=scatter_size, alpha=0.2, **kwargs)

        if mask:
            vertices = self.build_mask(n_dispersion_phi1=n_dispersion_phi1, n_dispersion_phi2=n_dispersion_phi2, refine_factor = refine_factor)

            x,y= vertices.T
            x=np.append(x, x[0])
            y=np.append(y, y[0])
            vert = np.transpose([x,y])

            ax[0].plot(vert.T[0],vert.T[1], 'k-')

        ax[0].set_xlim(-15,15)
        ax[0].set_ylim(-15,15)
        ax[0].set_xlabel('$\mu_{\phi_1}$ [mas yr$^{-1}$]')
        ax[0].set_ylabel('$\mu_{\phi_2}$ [mas yr$^{-1}$]')
        ax[0].set_title('Stream', fontsize='medium')

        # resize and fix column name
        scatter_size = 1./data_off['pmra_error']
        ax[1].scatter(data_off['PMPHI1'], data_off['PMPHI2'], c='k', s=scatter_size, alpha=0.2, **kwargs)

        if mask:
            ax[1].plot(vert.T[0],vert.T[1], 'k-')

        ax[1].set_xlim(-15,15)
        ax[1].set_ylim(-15,15)
        ax[1].set_xlabel('$\mu_{\phi_1}$ [mas yr$^{-1}$]')
        ax[1].set_ylabel('$\mu_{\phi_2}$ [mas yr$^{-1}$]')
        ax[1].set_title('Off stream', fontsize='medium')

        fig.tight_layout()

        if save:
            fig.savefig(f'proper_motions_{self.stream_obj.stream_name}_scatter.png',
                        bbox_inches='tight')
            fig.savefig(f'proper_motions_{self.stream_obj.stream_name}_scatter.pdf',
                        bbox_inches='tight')            
            
        return fig

    def plot_pm_hist(self, data, dx=0.5, norm=1, save=0, 
                pms=[None, None], match_norm=False, 
                stream_coords=True, reflex_corr=True, 
                zero_line=True, pm_lims=[-10, 10], **kwargs):
        # Code from Nora 
                    
        data_on = data[self.mask]
        data_off = data[self.off_mask]
        
        bins = np.arange(pm_lims[0], pm_lims[1] + dx / 2., dx)

        # data access depends on structure, needs to be adapted
        if stream_coords:
            h1 = np.histogram2d(data_on['PMPHI1'], data_on['PMPHI2'], bins)[0]
            h2 = np.histogram2d(data_off['PMPHI1'], data_off['PMPHI2'], bins)[0]
        else:
            if reflex_corr:
                h1 = np.histogram2d(data_on['PMRA'], data_on['PMDEC'], bins)[0]
                h2 = np.histogram2d(data_off['PMRA'], data_off['PMDEC'], bins)[0]
            else:
                h1 = np.histogram2d(data_on['PMRA0'], data_on['PMDEC0'], bins)[0]
                h2 = np.histogram2d(data_off['PMRA0'], data_off['PMDEC0'], bins)[0]

        
        # might need to normalise histogram for different areas of off stream mask for subtraction histogram
        
        h2 *= norm
        # print h1.sum(), h2.sum()
        diff = h1 - h2

        if match_norm:
            vmin, vmax = np.min([np.min(h1), np.min(h2), np.min(diff)]), np.max([np.max(h1), np.max(h2), np.max(diff)])
        else:
            vmin, vmax = None, None
        vmin = 0.

        offset = -dx / 2.

        fig, ((ax1, ax2, ax3)) = plt.subplots(1, 3, sharex='col', sharey='row', figsize=(15, 6))
        im1 = ax1.imshow(h1.T, extent=[bins.min() - offset, bins.max() - offset, bins.min() - offset, bins.max() - offset], origin='lower', vmin=vmin, vmax=vmax, interpolation='none', **kwargs)
        im2 = ax2.imshow(h2.T, extent=[bins.min() - offset, bins.max() - offset, bins.min() - offset, bins.max() - offset], origin='lower', vmin=vmin, vmax=vmax, interpolation='none', **kwargs)
        im3 = ax3.imshow(diff.T, extent=[bins.min() - offset, bins.max() - offset, bins.min() - offset, bins.max() - offset], origin='lower', vmin=vmin, vmax=vmax, interpolation='none', **kwargs)

        colorbar(im1)
        colorbar(im2)
        colorbar(im3)

        if (pms[0] is None) or (pms[1] is None):
            ax1.axvline(self.best_pm[0], ls='--', c='k', lw=1)
            ax2.axvline(self.best_pm[0], ls='--', c='k', lw=1)
            ax3.axvline(self.best_pm[0], ls='--', c='k', lw=1)
            ax1.axhline(self.best_pm[1], ls='--', c='k', lw=1)
            ax2.axhline(self.best_pm[1], ls='--', c='k', lw=1)
            ax3.axhline(self.best_pm[1], ls='--', c='k', lw=1)
        else:
            ax1.axvline(pms[0], ls='--', c='k', lw=1)
            ax2.axvline(pms[0], ls='--', c='k', lw=1)
            ax3.axvline(pms[0], ls='--', c='k', lw=1)
            ax1.axhline(pms[1], ls='--', c='k', lw=1)
            ax2.axhline(pms[1], ls='--', c='k', lw=1)
            ax3.axhline(pms[1], ls='--', c='k', lw=1)           
        # f.suptitle(r'$\mathrm{%s}$' % stream['name'], fontsize=20)

        ax1.set_title(r'$\mathrm{on-stream}$')
        ax2.set_title(r'$\mathrm{off-stream}$')
        ax3.set_title(r'$\mathrm{residual}$')

        if stream_coords:
            ax1.set_xlabel(r'$\mu_1$')
            ax2.set_xlabel(r'$\mu_1$')
            ax3.set_xlabel(r'$\mu_1$')

            ax1.set_ylabel(r'$\mu_2$')
        else:
            if reflex_corr:
                ax1.set_xlabel(r'$\mu_{ra}$')
                ax2.set_xlabel(r'$\mu_{ra}$')
                ax3.set_xlabel(r'$\mu_{ra}$')

                ax1.set_ylabel(r'$\mu_{dec}$')

        if zero_line:
            ax1.axhline(0, ls='--', c='k', lw=0.5)
            ax2.axhline(0, ls='--', c='k', lw=0.5)
            ax3.axhline(0, ls='--', c='k', lw=0.5)

        ax1.set_xlim(bins.min(), bins.max())
        ax2.set_xlim(bins.min(), bins.max())
        ax3.set_xlim(bins.min(), bins.max())
        ax1.set_ylim(bins.min(), bins.max())
        ax2.set_ylim(bins.min(), bins.max())
        ax3.set_ylim(bins.min(), bins.max())

        fig.tight_layout()

        if save:
            fig.savefig(f'proper_motions_{self.stream_obj.stream_name}_hist.png',
                        bbox_inches='tight')
            fig.savefig(f'proper_motions_{self.stream_obj.stream_name}_hist.pdf',
                        bbox_inches='tight')            
            # plt.close()
        return fig
    
        # ========================= added Nov 3 (need checking) ========================

    def find_peak_location(self, data, x_width=3., y_width=3., draw_histograms=True):
        '''
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
        '''
        from matplotlib.colors import LogNorm
        from astropy.modeling import models, fitting
        	
        x_center, y_center = self.best_pm_phi1_mean, self.best_pm_phi2_mean
        xmin, xmax, ymin, ymax = x_center-x_width, x_center+x_width, y_center-y_width, y_center+y_width
		
        # on-stream data into 2D histogram
        H1, x_edges, y_edges = np.histogram2d(self.pm_phi1_cosphi2, self.pm_phi2, bins = (np.arange(xmin, xmax, x_width/20),np.arange(ymin, ymax, y_width/20))) #explore differtent bin sizes
        
        stream_off = self.off_mask
        # off-stream data
        H2, x_edges, y_edges = np.histogram2d(data['pm_phi1_cosphi2'][stream_off], 
                                        data['pm_phi2'][stream_off], 
                                             bins = (np.arange(xmin, xmax, x_width/20), np.arange(ymin, ymax, y_width/20)))
        
        
        hist = H1-self.mask.sum()/self.off_mask.sum()*H2 #check this scale factor --> self.mask.sum()/self.off_mask.sum()

        # fitting 2D gaussian (Code from Ani)
        #Find overdensity
        ind = np.unravel_index(hist.argmax(), hist.shape)
        #Fit 2d gaussian to a zoomed in histogram
        hist_zoom = np.transpose(hist[(ind[0]-6):(ind[0]+7), (ind[1]-6):(ind[1]+7)])
        g_init = models.Gaussian2D(amplitude=1., x_mean=x_edges[ind[0]], y_mean=y_edges[ind[1]], x_stddev = 0.5, y_stddev = 0.5)
        fit_g = fitting.LevMarLSQFitter()
        x,y = np.meshgrid(x_edges[(ind[0]-6):(ind[0]+7)], y_edges[(ind[1]-6):(ind[1]+7)])
        g = fit_g(g_init, x, y, hist_zoom)

        # the mean position of the peak
        pm_x_cen = g.x_mean + 0.5*(x_edges[1]-x_edges[0])        # half with of the bin
        pm_y_cen = g.y_mean + 0.5*(y_edges[1]-y_edges[0]) 
        # the std's of the peak
        x_std, y_std = g.x_stddev.value, g.y_stddev.value
        
        # draw proper motions of the on-stream, off-stream and residual histogram
        if draw_histograms:
            fig, axes = plt.subplots(1, 3,
                figsize=(15, 5), sharex=True, sharey=True, 
                constrained_layout=True
            )

            bins = np.linspace(-15, 15, 128)
            H1, xe, ye, _ = axes[0].hist2d(
                self.pm_phi1_cosphi2,
                self.pm_phi2,
                bins=bins,
                norm=LogNorm()
            )
            H2, *_ = axes[1].hist2d(
                data['pm_phi1_cosphi2'][stream_off],
                data['pm_phi2'][stream_off],
                bins=bins,
                norm=LogNorm()
            )

            axes[2].pcolormesh(xe, ye, (H1-self.mask.sum()/self.off_mask.sum()*H2).T, cmap='RdBu', vmin=-10, vmax=10)
            t = np.linspace(0, 2*np.pi, 100)
            axes[2].plot(pm_x_cen + x_std*np.cos(t), pm_y_cen + y_std*np.sin(t), c='green')
            #axes[2].set_xlim(xmin,xmax)
            #axes[2].set_ylim(ymin,ymax)

            for ax in axes[:2]:
                ax.plot(
                    self.galstream_pm_phi1_cosphi2,
                    self.galstream_pm_phi2,
                    color='cyan'
                )


                ax.set_xlabel(r'$\mu_{\phi_1}$', fontsize=20)
            axes[0].set_ylabel(r'$\mu_{\phi_2}$', fontsize=20)
		
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