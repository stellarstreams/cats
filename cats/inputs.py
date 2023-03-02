from collections import OrderedDict as odict
# config file with all the inputs

stream_inputs = odict([
    ('GD-1',
     dict(
         #galstream stuff
         short_name='GD-1',
         pawprint_id='pricewhelan2018',
         #stream stuff
         width=2.0, # degrees (add units in pawprint)
         #data stuff
         survey='PS1',
         band1='g',
         band2='r',
         mag1='g0',
         mag2='r0',
         minmag=16.,
         maxmag=24.0,
         #isochrone stuff
         age=11.8, # Gyr
         feh=-1.5, 
         distance=8.3, #kpc
         alpha=0, #don't think we actually use this
         scale_err=2,
         base_err=0.075,
         bin_sizes=[0.03,0.2], # xbin and ybin width for CMD
     )),
     ('Pal5',
     dict(
         #galstream stuff
         short_name='Pal5',
         pawprint_id='pricewhelan2019',
         #stream stuff
         width=1.0, # degrees (add units in pawprint)
         #data stuff
         survey='PS1',
         band1='g',
         band2='i',
         mag1='g0',
         mag2='i0',
         minmag=16.,
         maxmag=24.0,
         #isochrone stuff
         age=12, # Gyr
         feh=-1.4, 
         distance=20.9, #kpc
         alpha=0, #don't think we actually use this
         scale_err=2,
         base_err=0.075,
         bin_sizes=[0.03,0.2])),
    ('Jhelum',
     dict(
         #galstream stuff
         short_name='Jhelum-b',
         pawprint_id='bonaca2019',
         #stream stuff
         width=2.0, # degrees (add units in pawprint)
         #data stuff
         survey='PS1',
         band1='g',
         band2='r',
         mag1='g0',
         mag2='r0',
         minmag=16.,
         maxmag=24.0,
         #isochrone stuff
         age=11, # Gyr
         feh=-1.5, 
         distance=13.2, #kpc
         alpha=0, #don't think we actually use this
         scale_err=2,
         base_err=0.075,
         bin_sizes=[0.03,0.2])),])
