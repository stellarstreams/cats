import numpy as np
from matplotlib.path import Path
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.vizier import Vizier
from shapely import geometry
from shapely.ops import transform
from scipy.interpolate import InterpolatedUnivariateSpline
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import galstreams

class sky_polygon(geometry.Polygon):
    '''
    Wrapper for shapely's plotting
    '''
    def display(self, **kwargs):
        '''
        Add matplotlib path as patch to current axis
        '''
        plt.gca().add_patch(patches.PathPatch(self.path(), **kwargs))
    def spherical_area(self):
        '''
        compute area of on-sky polygon (uses Lambert equal area projection)
        '''
        return transform(lambda x, y: (np.sqrt(2./(1-np.sin(np.deg2rad(y))))*np.cos(np.deg2rad(y))*np.cos(np.deg2rad(x)), 
                                       np.sqrt(2./(1-np.sin(np.deg2rad(y))))*np.cos(np.deg2rad(y))*np.sin(np.deg2rad(x))), self).area
    def path(self):
        '''
        generate matplotlib path from shapely polygon
        '''
        return Path.make_compound_path(
                Path(np.asarray(self.exterior.coords)),
                *[Path(np.asarray(ring.coords)) for ring in self.interiors])
    def remove_hole(self, xy, r):
        '''
        remove a circular region from the polygon
        
        Parameters
        ----------
        xy : array
            centre of the hold 
        r : float
            radius
        '''
        p = geometry.point.Point(*xy)
        circle = p.buffer(r)
        return sky_polygon(self-circle)
    def __add__(self, a):
        '''
        Combine two overlapping polygons
        
        Parameters
        ----------
        a: sky_polygon
            polygon to combine
        '''
        return sky_polygon(unary_union([self,a]))
    def __subtract__(self, a):
        '''
        subtract polygon
        
        Parameters
        ----------
        a: sky_polygon
            polygon to combine
        '''
        return sky_polygon(self-a)

class sky_polygon_multi(object):
    def __init__(self,x):
        '''
        Multiply non-overlapping polygons
        
        Parameters
        ----------
        x : list of shapely.Polygons
        '''
        self.poly_list = x
        for i in range(len(self.poly_list)):
            flgg=[]
            for j in range(len(self.poly_list)):
                if i!=j and self.poly_list[i].intersects(self.poly_list[j]):
                    self.poly_list[i] = sky_polygon(unary_union([self.poly_list[i], self.poly_list[j]]))
                    flgg+=[j]
            self.poly_list = [self.poly_list[k] for k in range(len(self.poly_list)) if k not in flgg]
                    
    def display(self, **kwargs):
        '''
        Add matplotlib path as patch to current axis
        '''
        plt.gca().add_patch(patches.PathPatch(self.path(), **kwargs))
    def spherical_area(self):
        '''
        compute area of on-sky polygon (uses Lambert equal area projection)
        '''
        return np.sum([x.spherical_area() for x in self.poly_list])
    def path(self):
        '''
        generate matplotlib path from shapely multipolygon
        '''
        pth = self.poly_list[0].path()
        for i in range(1,len(self.poly_list)):
            pth=Path.make_compound_path(pth, self.poly_list[i].path())
        return pth
    
class extrapolated_track(object):
    
    def __init__(self, track, smoothing_scale = 0.2*u.deg):
        """
        Extrapolation routine for a track -- to stop wiggles at the ends dominating
        
        Parameters
        ----------
        track: galstreams track
        smoothing_scale: Quantity
            smoo
        """
        dt = int(smoothing_scale / np.nanmedian(np.diff(track.phi1)))
        
        diff10 = np.sign(track.phi2.degree[dt:]-track.phi2.degree[:-dt])
        trim_low = np.argwhere(np.diff(track.phi2.degree)*diff10[0]>=0.)[0][0]
        trim_high = np.argwhere(np.diff(track.phi2.degree)*diff10[-1]>=0.)[-1][0]
        
        self.extrap_track = InterpolatedUnivariateSpline(
            track.phi1.degree[trim_low:trim_high],
            track.phi2.degree[trim_low:trim_high],
            k=1
        )
        
    def __call__(self, x):
        '''
        Parameters
        ----------
        x : SkyCoord
            location of track evaluation
        '''
        return self.extrap_track(x.to(u.deg))*u.deg


class globular_cluster_table(object):
    
    def __init__(self):
        """
        Loading globular cluster data
        """
        self.data = Vizier(columns=['_RAJ2000','_DEJ2000','Rh']).get_catalogs(catalog='VII/202')[0]
        self.skycoord = SkyCoord(ra=self.data['_RAJ2000'],
                                 dec=self.data['_DEJ2000'],frame='icrs')
        self.rh = self.data['Rh']
        
        ## Fill in missing half-light radii
        self.rh[self.rh.mask]=np.nanmedian(self.rh[~self.rh.mask])
        
        self.nclusters = len(self.data)
        
        
class on_sky_galstream_selections(object):
    
    def __init__(self, stream_name, stream_ref=None, gs_table=None):
        """
        A class for constructing an on-sky polygon for on-/off-stream selections
        based on the galstreams tracks
        
        Parameters
        ----------
        stream_name : str
            The name of the stream to consider
        stream_ref : str
            The specific galstreams version/reference of the stream (e.g. PW19 for Pal5)
        gs_table : galstreams.MWStreams instance or None
            if None, initialize from galstreams
        """
        if gs_table is None:
            gs_table = galstreams.MWStreams()
    
        self.stream_name = stream_name
        if stream_ref is None:
            internal_name = gs_table.get_track_names_for_stream(self.stream_name, 
                                                                On_only=True)[0]
            self.stream_ref = internal_name.split('-')[-1]
        else:
            internal_name = self.stream_name + '-' + stream_ref
            self.stream_ref = stream_ref
            
        self.stream_frame = gs_table[internal_name].stream_frame
        
        self.track = extrapolated_track(gs_table[internal_name].track.transform_to(
                                        self.stream_frame))
        
        self.length =gs_table.summary.length[internal_name]*u.deg
    
        self.gc = globular_cluster_table()
    
    def add_hole(self, path, coord, size):
        """
        Add a circular cutout in a path
        
        Parameters
        ----------
        path: matplotlib.Path
            the overall on-sky polygon from which holes are removed
        coord: SkyCoord instance
            the centre of the hole
        size: Quantity
            the size of the hole
        """
        ## Let's just say that holes are circles in (phi1, phi2) -- seems OK unless stream is fat
        
        p12 = coord.transform_to(self.stream_frame)
        return path.remove_hole([p12.phi1.deg, p12.phi2.deg], 
                                size.to(u.deg).value)
        pth=Path.circle([p12.phi1.deg, p12.phi2.deg], 
                        size.to(u.deg).value)
        if path.contains_path(pth):
            path = Path.make_compound_path(path, pth)
        return path
    
    def add_gc_holes(self, path, rh_multiple=3.):
        """
        Loop over gc and add holes
        
        Parameters
        ----------
        path: matplotlib.Path
            the overall on-sky polygon from which holes are removed
        rh_multiple: 
            the multiple of the half-light radius to cutout
        """
        
        for i in range(self.gc.nclusters):
            path = self.add_hole(path, self.gc.skycoord[i], rh_multiple*self.gc.rh.quantity[i])
            
        return path
    
    def make_path_poly(self, length_range, width, offset):
        """
        Build an on-sky polygon
        
        Parameters
        ----------
        length_range: array
            array of lengths
        rh_multiple: 
            the multiple of the half-light radius to cutout
        """
        path = np.concatenate([np.vstack([length_range,self.track(length_range)+width*.5+offset]),
                               np.vstack([length_range,self.track(length_range)-width*.5+offset])[:,::-1]],axis=1)
        poly = sky_polygon(path.value.T)
        
        return poly
    
    def make_onsky_footprint(self, priority, offset):
        """
        Auto-make an on-sky footprint
        
        Parameters
        ----------
        priority : str
            either 'completeness' or 'purity' -- completeness on-sky footprints are more generous than purity
        offset : Quantity
            the offset of the polygon in phi2
        """
        
        width_options = {'completeness':1.*u.deg, 'purity':0.4*u.deg}
        length_options = {'completeness':1.2, 'purity':0.5}
        
        length_range = self.length * length_options[priority] * np.linspace(-0.5, 0.5, 300)
        
        return self.make_path_poly(length_range, width_options[priority], offset)
    
    def on_stream(self, priority='purity'):
        """
        Auto-make an on-stream on-sky footprint
        
        Parameters
        ----------
        priority : str
            either 'completeness' or 'purity' -- completeness on-sky footprints are more generous than purity
        """
        return self.add_gc_holes(self.make_onsky_footprint(priority, 0.*u.deg))
        
    def off_stream(self, shift=1.*u.deg, priority='purity'):
        """
        Auto-make an off-stream on-sky footprint
        
        Parameters
        ----------
        shift : Quantity
            offset of the off-stream polygon in phi2
        priority : str
            either 'completeness' or 'purity' -- completeness on-sky footprints are more generous than purity
        """
        return sky_polygon_multi([self.add_gc_holes(self.make_onsky_footprint(priority, -shift)),
                                  self.add_gc_holes(self.make_onsky_footprint(priority, shift))])
        
    def display(self,equal=False):
        """
        Display the on-/off-stream polygons for completeness and purity priority
        
        Parameters
        ----------
        equal : bool
            make axes equal
        """
        off_stream_complete = self.off_stream(priority='completeness')
        
        f,ax=plt.subplots()
        
        self.on_stream().display(facecolor='C1', lw=1)
        self.off_stream().display(facecolor='C0', lw=1)
        self.on_stream(priority='completeness').display(facecolor='C1', lw=1, alpha=0.3)
        self.off_stream(priority='completeness').display(facecolor='C0', lw=1, alpha=0.3)
        
        ax.set_xlim(np.min(off_stream_complete.path().vertices[:,0]),
                    np.max(off_stream_complete.path().vertices[:,0]))
        ax.set_ylim(np.min(off_stream_complete.path().vertices[:,1]),
                    np.max(off_stream_complete.path().vertices[:,1]))
        
        if equal:
            plt.gca().set_aspect('equal')
        plt.xlabel(r'$\phi_1$ [deg]')
        plt.ylabel(r'$\phi_2$ [deg]')