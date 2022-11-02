import galstreams as gst
from matplotlib.path import Path as mpl_path
from astropy.coordinates import CoordFrame #shorthand
from astropy.table import Table

#use ecsv for polygon specification?
#asdf packages these together 


class densityClass: #TODO: how to represent densities?

class Footprint2D(dict):
    def __init__(self, vertex_coordinates, footprint_type):
        if footprint_type='sky':
            #steal Cecilia's implementation from galstreams and return mpl_path and vertices in skycoords
            self.vertices = SkyCoord(vertex_coordinates)

        if footprint_type='cartesian':
            #should return a mpl_path and vertices in regular old cartesian coords
            self.vertices = vertex_coordinates

        self.poly = mpl_path(vertices)

    @classmethod
    def from_vertices(cls, vertex_coordinates, footprint_type):
        return cls(vertices,footprint_type)

    @classmethod
    def from_box(cls, min1, max1, min2, max2, footprint_type):
        vertices = get_vertices_from_box(min1, max1, min2, max2)
        return cls(vertices,footprint_type)

    @classmethod
    def from_file(cls,fname):
        with Table.read(fname, format='ascii.ecsv') as t:
            vertices = t['vertices']
            footprint_type = t['footprint_type']
        return cls(vertices,footprint_type)

    def get_vertices_from_box(min1, max1, min2, max2):
        return [[min1,min2],[min1,max2],[max1,min2],[max1,max2]]

    def save(fname):
        #TODO - save as .ecsv




class pawprintClass(dict):
    '''Dictionary class to store a "pawprint": 
        polygons in multiple observational spaces that define the initial selection 
        used for stream track modeling, 
        membership calculation / density modeling, and background modeling.
        
        New convention: everything is in phi1 phi2 (don't cross the streams)

        '''

    def __init__(self, stream_name, pawprint_ID, from_galstreams=True):
        
        #we need a way to specify and load the vertices for the footprints. 
        #How do we want to do it? 
        #I sketched here passing the name of the stream 
        #but we could also pass a list of vertices provided by WG2/3

        self.stream_name = stream_name
        self.pawprint_ID = pawprint_ID
        self.stream_frame = CoordFrame(self,input_data_specifier?)

                #right now loading tracks via galstreams, later include updates to track
        if from_galstreams:
            track_file = _make_track_file_name(stream_name,pawprint_ID) #todo
            summary_file = _make_summary_file_name(stream_name,pawprint_ID) #todo
            self.track = gst.Track6D(stream_name, track_file, summary_file)

        self.width = 0.0#(lambda phi1: return width(phi1))



        stream_vertices,background_vertices = load_sky_vertices(self)

        self.skyprint = {'stream':Footprint2D(stream_vertices,type='sky'), 
                        'background':Footprint2D(background_vertices,type='sky')} 
                        #write new thing to rept footprint a la galstreams, but for any 2 

        

        #WG3: how to implement distance dependence in isochrone selections?
        self.cmd_filters = {'poly1':['g', 'bprp'], 'poly2':['h','jk']} # list of which filters are in which footprint in the list. need to match keys in starlist
        self.cmdprint = {'poly1':Footprint2D(load_cmd_vertices(self,'poly1'),type='cartesian'),...} #polygon(s) in cmd space
 

        self.pmprint = {} = Footprint2D(load_pm_vertices(self),type='sky') #polygon(s) in proper-motion space mu_phi1, mu_phi2


        

    
    def save_pawprint(self):
        '''TODO: save as YAML (Eduardo)'''

