import numpy as np
from astropy.coordinates import SkyCoord, CoordFrame

#notes on data types
dtypes_dict = {
	#identifiers
	'sourceID': 'u8', #Gaia DR3 source id, unsigned integer
	'sourceID_version': 'u2', #integer representing version of gaia catalog that the source ID refers to
	'streamID': 's10', #10-char string? check galstreams
	'crossmatches':dict(), #dictionary to store IDs in other surveys

	#phase space
	'w': {}#dict of phase space coordinates
	'w_uncert': {}#dictionary of uncertainties per coordinate
	#todo: add function to convert w to SkyCoord
	#todo: bring in pyia to account for transformations to uncertainties

	#precompute and store phi1 and phi2 for a desired coordinate frame
	'phi1': 'f8'
	'phi2': 'f8'
	'rotation': CoordFrame() #placeholder

	#flag for variable stars
	'variability': 'u2', #0 if not variable, 1 if variable

	#photometry from Gaia used for selections in pawprint
	'phot_g_mean_mag':  'f4',
	'phot_rp_mean_mag': 'f4',
	'phot_g_mean_mag_error': 'f4',
	'phot_rp_mean_mag_error':  'f4',

	#TODO WG3: how to store extinction?

	#consensus chemistry
	'feh_logeps': 'f4', #sun-independent value
	'feh': 'f4', #iron abundance as [Fe/H]
	'feh_solar': 'f4', #solar value of [Fe/H] for this star

	'alpha_logeps': 'f4', #sun-independent value
	'alpha_fe': 'f4', #alpha abundance as [alpha/Fe]
	'alpha_solar': 'f4', #solar value of [alpha/H] for this star


	#references: assumes that sky position, photometry, and PM are from Gaia
	'refs': {
		'distance': ['s19'] #ADS bibcode (or pointer to docs) for distance measurement; can be a list
		'rv': ['s19']
		'feh': ['s19']
		'alpha': ['s19']
		'variability': ['s19']
		'extinction': ['s19']
	}

	#TODO: membership likelihoods

}

def get_phase_space():
	'''queries gaia to initialize SkyCoords for phase-space position and uncertainty'''

def get_gaia_photometry():
	'''queries gaia to get photometry and uncertainties, returned as 32 bit floats'''

def get_abundances():
	'''queries ancillary data tables for spectroscopy, probably given some options that perhaps user can set'''


class starClass(dict):
	def __init__(self, streamID):
		#starclass 
		self.data = starDataClass(streamID)


class starDataClass(dict):
	'''dictionary class to store __measured__ attributes for one star in the catalog'''
    def __init__(self, streamID): 
    	
    	self.sourceID_version = np.uint(3) #default for now is DR3
    	self.streamID = np.str_(streamID)

    	self.sourceID = load_stream(self.streamID) #read from adrian's initial files

    	self.nstars = len(self.sourceID)

    	self.crossmatches = {}

    	self.w, self.w_uncert = phasespace_to_skycoords() #TODO:function to return sky coordinates and uncertainties in two skycoord objects by querying Gaia
    	
    	#flexible magnitudes - pin down standardised naming convention
    	#my proposal: [survey]_[filter]
    	#uncertainties and extinctions are specified by the same tags
    	self.mags = {'gaia_g':np.array(nstars,dtype='f4'),'gaia_rp':np.array(nstars,dtype='f4'), }
    	self.mag_uncert = {'gaia_g':np.array(nstars,dtype='f4'),'gaia_rp':np.array(nstars,dtype='f4'), }
    	self.ext = {'gaia_g':np.array(nstars,dtype='f4'),'gaia_rp':np.array(nstars,dtype='f4'), }


    	self.variability = np.array(nstars,dtype='u2')

    	self.feh = np.masked_array(nstars,dtype='f4') 
    	self.feh_logeps = np.masked_array(nstars,dtype='f4')
    	self.feh_solar = np.masked_array(nstars,dtype='f4')
    	self.alpha_logeps = np.masked_array(nstars,dtype='f4')
    	self.alpha_fe = np.masked_array(nstars,dtype='f4')
    	self.alpha_solar = np.masked_array(nstars,dtype='f4') 

    	self.refs = {
			'distance': np.array(nstars,dtype='s19') #ADS bibcode (or pointer to doi) for distance measurement; can be a list
			'rv': np.array(nstars,dtype='s19')
			'feh': np.array(nstars,dtype='s19')
			'alpha': np.array(nstars,dtype='s19')
			'variability': np.array(nstars,dtype='s19')
    	}
    	
    	#some stuff could be read directly and stored
    	if get_gaia:
	    	get_gaia_photometry(self) #load gaia photometry in from catalog
    		get_abundances(self) #load abundances from detailed spectroscopoc tables

class starDerivedClass(dict):
	stuff



def _inside_poly(data, vertices):
        '''This function takes a list of points (data) and returns a boolean mask that is True for all points inside the polygon defined by vertices'''
        return mpl_path(vertices).contains_points(data)

def makeMask(self, pawprint, what):
        '''take in some data and return masks for stuff in the pawprint (basically by successively applying _inside_poly)'''
        #returns mask with same dimension as data



