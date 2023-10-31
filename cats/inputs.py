from collections import OrderedDict as odict
# config file with all the inputs

stream_inputs = odict([
    ('GD-1',
     dict(
         #galstream stuff
         short_name='GD-1',
         pawprint_id='pricewhelan2018',
         #stream stuff
         width=2.0, # full width in degrees (add units in pawprint)
         #data stuff
         phot_survey='PS1',
         band1='g',
         band2='r',
         mag='g0',
         color1='g0',
         color2='r0',
         minmag=16.,
         maxmag=24.0,
         #isochrone stuff
         age=11.8, # Gyr
         feh=-1.5, 
         distance=8.3, #kpc
         turnoff=17.8, # mag of MS turnoff
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
         phot_survey='PS1',
         band1='g',
         band2='r',
         mag='g0',
         color1='g0',
         color2='r0',
         minmag=16.,
         maxmag=24.0,
         #isochrone stuff
         age=12, # Gyr
         feh=-1.4, 
         distance=20.9, #kpc
         turnoff=15, # mag of MS turnoff
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
         phot_survey='des',
         band1='g',
         band2='r',
         mag='g0',
         color1='g0',
         color2='r0',
         minmag=16.,
         maxmag=24.0,
         #isochrone stuff
         age=12, # Gyr
         feh=-1.7, 
         distance=13.2, #kpc
         turnoff=18.7, # mag of MS turnoff
         alpha=0, #don't think we actually use this
         scale_err=2,
         base_err=0.075,
         bin_sizes=[0.03,0.2])),
    ('Fjorm-M68',
     dict(
         #galstream stuff
         short_name='M68-Fjorm',
         #pawprint_id='ibata2021',
         pawprint_id='palau2019',
         #stream stuff
         width=1, # TOTAL width degrees, recommend 2sigma if known
         #data stuff
         phot_survey='Gaia',
         band1='BP',
         band2='RP',
         mag='G0',
         color1='BP0',
         color2='RP0',
         minmag=16.,
         maxmag=24.0,
         #isochrone stuff
         age=11.2, # Gyr
         feh=-2.2, 
         distance=6, #kpc
         turnoff=17.0, # mag of MS turnoff
         alpha=0, #don't think we actually use this
         scale_err=2,
         base_err=0.075,
         bin_sizes=[0.03,0.2])),
    ('Gunnthra',
     dict(
         #galstream stuff
         short_name='Gunnthra',
         #pawprint_id='ibata2021',
         pawprint_id='ibata2021',
         #stream stuff
         width=5, # TOTAL width degrees, recommend 2sigma if known
         #data stuff
         phot_survey='Gaia',
         band1='BP',
         band2='RP',
         mag='G0',
         color1='BP0',
         color2='RP0',
         minmag=16.,
         maxmag=24.0,
         #isochrone stuff
         age=10, # Gyr
         feh=-0.7, 
         distance=3, #kpc
         turnoff=17.0, # mag of MS turnoff
         alpha=0, #don't think we actually use this
         scale_err=2,
         base_err=0.075,
         bin_sizes=[0.03,0.2]))])
