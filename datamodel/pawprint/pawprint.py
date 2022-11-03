import galstreams as gst
from matplotlib.path import Path as mpl_path
import astropy.table as apt
import asdf
import astropy.units as u
from astropy.coordinates import SkyCoord
from gala.coordinates import GreatCircleICRSFrame
import numpy as np


#class densityClass: #TODO: how to represent densities?

class Footprint2D(dict):
    def __init__(self, vertex_coordinates, footprint_type, stream_frame=None):
        if footprint_type=='sky':
            if isinstance(vertex_coordinates, SkyCoord):
                vc = vertex_coordinates
            else:
                vc = SkyCoord(vertex_coordinates)
            self.edges = vc
            self.vertices = np.array([vc.transform_to(stream_frame).phi1, vc.transform_to(stream_frame).phi2]).T

        elif footprint_type=='cartesian':
            self.edges = vertex_coordinates
            self.vertices = vertex_coordinates

        self.stream_frame = stream_frame
        self.footprint = mpl_path(self.vertices)
        
    @classmethod
    def from_vertices(cls, vertex_coordinates, footprint_type):
        return cls(vertices,footprint_type)

    @classmethod
    def from_box(cls, min1, max1, min2, max2, footprint_type):
        vertices = get_vertices_from_box(min1, max1, min2, max2)
        return cls(vertices,footprint_type)

    @classmethod
    def from_file(cls,fname):
        with Table.read(fname) as t:
            vertices = t['vertices']
            footprint_type = t['footprint_type']
        return cls(vertices,footprint_type)

    def get_vertices_from_box(min1, max1, min2, max2):
        return [[min1,min2],[min1,max2],[max1,min2],[max1,max2]]
    
    def inside_footprint(self,data):
        if isinstance(data, SkyCoord):
            if self.stream_frame is None:
                print("can't!")
                return
            else:
                pts = np.array([data.transform_to(self.stream_frame).phi1.value,data.transform_to(self.stream_frame).phi2.value]).T
                return self.footprint.contains_points(pts)
        else:
            return self.footprint.contains_points(data)

    def export(self):
        data = {}
        data['stream_frame']=self.stream_frame
        data['vertices']=self.vertices
        data['footprint_type']=self.footprint_type
        return data






class pawprintClass(dict):
    '''Dictionary class to store a "pawprint": 
        polygons in multiple observational spaces that define the initial selection 
        used for stream track modeling, 
        membership calculation / density modeling, and background modeling.
        
        New convention: everything is in phi1 phi2 (don't cross the streams)

        '''

    
    def __init__(self, data):
        
        #we need a way to specify and load the vertices for the footprints. 
        #How do we want to do it? 
        #I sketched here passing the name of the stream 
        #but we could also pass a list of vertices provided by WG2/3

        self.stream_name = data['stream_name']
        self.pawprint_ID = data['pawprint_ID']
        self.stream_frame = data['stream_frame']
        self.width = data['width']
        self.skyprint = {'stream':Footprint2D(data['stream_vertices'],footprint_type='sky',stream_frame=self.stream_frame), 
                        'background':Footprint2D(data['background_vertices'],footprint_type='sky',stream_frame=self.stream_frame)} 
        #WG3: how to implement distance dependence in isochrone selections?
        self.cmd_filters = data['cmd_filters']
        self.cmdprint = {}
        if self.cmd_filters is not None:
            for k in data.cmd_filters.keys():
                self.cmdprint[k] = Footprint2D(data['cmd_vertices'][k], footprint_type='cartesian')
        if data['pm_vertices'] is not None:
            self.pmprint = Footprint2D(data['pm_vertices'],footprint_type='sky') #polygon(s) in proper-motion space mu_phi1, mu_phi2
        
        self.track = data['track']

    @classmethod
    def load_pawprint(cls,fname):
        #TODO: load all the stuff from the file into a dictionary called data and call the main init function
        ...
        return cls(data)
        

    
    @classmethod
    def pawprint_from_galstreams(cls,stream_name,pawprint_ID):
        
        
        galstreams_tracks = '../../galstreams/galstreams/tracks/'
    
        def _make_track_file_name(stream_name,pawprint_ID):
            return galstreams_tracks+'track.st.'+stream_name+'.'+pawprint_ID+".ecsv"
    
        def _make_summary_file_name(stream_name,pawprint_ID):
            return galstreams_tracks+'track.st.'+stream_name+'.'+pawprint_ID+".summary.ecsv"
        
        def _get_stream_frame_from_file(summary_file):
            t = apt.QTable.read(summary_file)
        
            x = dict()
            atts = [x.replace('mid.','') for x in t.keys() if 'mid' in x ]
            for att in atts:  #we're effectively looping over skycoords defined for mid here (ra, dec, ...)
                x[att] = t[f'mid.{att}'][0]   #<- make sure to set it up as a scalar. if not, frame conversions get into trouble
            mid_point = SkyCoord(**x) 

            x = dict()
            atts = [x.replace('pole.','') for x in t.keys() if 'pole' in x ]
            for att in atts:  #we're effectively looping over skycoords defined for pole here (ra, dec, ...)
                x[att] = t[f'pole.{att}'][0]
            #Make sure to set the pole's distance attribute to 1 (zero causes problems, when transforming to stream frame coords)
            x["distance"] = 1.*u.kpc   #it shouldn't matter, but if it's zero it does crazy things
            mid_pole = SkyCoord(**x)

            return GreatCircleICRSFrame(pole=mid_pole, ra0=mid_point.icrs.ra)

        data = {}
        data['stream_name'] = stream_name
        data['pawprint_ID'] = pawprint_ID

        track_file = _make_track_file_name(stream_name,pawprint_ID) 
        summary_file = _make_summary_file_name(stream_name,pawprint_ID) 
        data['stream_frame'] = _get_stream_frame_from_file(summary_file)

        data['track'] = gst.Track6D(stream_name=data['stream_name'], track_name=data['pawprint_ID'], track_file=track_file, summary_file=summary_file)
        data['width'] = 1.0*u.deg
        data['stream_vertices'] = data['track'].create_sky_polygon_footprint_from_track(width=data['width'], phi2_offset=0.*u.deg)
        data['background_vertices'] = data['track'].create_sky_polygon_footprint_from_track(width=data['width'], phi2_offset=3.*u.deg)
        data['cmd_filters'] = None
        data['cmd_vertices'] = None
        data['pm_vertices'] = None


        return cls(data)


    def add_cmd_footprint(self, new_footprint, color, mag, name):
        if self.cmd_filters is None:
            self.cmd_filters = {name:[color, mag]}
            self.cmdprint = {name: new_footprint}
        else:
            self.cmd_filters[name] = [color,mag]
            self.cmdprint[name] = new_footprint

    def add_pm_footprint(self, new_footprint, name):
        if self.pmprint is None:
            self.pmprint = {name: new_footprint}
        else:
            self.pmprint[name] = new_footprint


    def save_pawprint(self):
        fname = stream_name+pawprint_ID+'.asdf'
        tree = {
            'stream_name':stream_name,
            'pawprint_ID':pawprint_ID,
            'stream_frame':stream_frame,
            'cmd_filters': cmd_filters,
            'width':width,
            'on_stream':{
                        'sky':skyprint['stream'].export(),
                        'cmd':cmdprint.export(),
                        'pm':pmprint.export()
                        },
            'off_stream':skyprint['background'].export(),
            'track':track
        }
        out = asdf.AsdfFile(tree)
        out.write_to(fname)
    


