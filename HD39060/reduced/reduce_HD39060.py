"""

This script was written to run on CASA 6.4.1.12 to reduce the following datasets belonging to project 2022.1.00338.L 

Target: HD39069 (Beta Pic)

- ACA observed with 2 executions in Oct 2013
- 12m SB observed with 1 execution blocks and 2 fields in Dec 2013
- 12m LB observed with 3 execution blocks in Aug 2015

reducer: S. Marino

"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

### import functions for data reduction
#sys.path.append('../../reduction_scripts/')
#from functions_data_reduction import *
execfile('../../functions_data_reduction.py') # it needs to be run twice for correct import

### system parameters
system={'target'      : 'HD39060', # used to name files
        'field_name'  : 'beta_pictoris', # for ACA, field name in ms
        'Radius'      : 50., # arcsec, maximum radius of interest used to calculate maximum freq and time averaging. Set to 500 au or at least 2x disc outer edge or FWHM of Antenna primary beam
        'rv'       : 16.84, # km/s, radial velocity in barycentric reference frame from Gaia DR3
        'dvel'        : 50.0,      # km/s, line width. It is used to flag channels at +-dvel and to image a cube with a width of +-3dvel
        }


##############################
########## ACA ##############
#############################

cfg='ACA'
path_calibrated='../calibrated/{}/'.format(cfg)

### for ACA 
robust=2.0
Npix=256
dpix=0.4  # arcsec
beam=4.0  # arcsec

# ##################
# #### TASKS #######
# ##################

do_mstransform = False  # remove calibrators and change to BARY reference frame

# continuum
do_select_cont = False   # select continuum by flagging
do_fav         = False # ferquency average for continuum
do_tav         = False  # time average
do_extractvis  = False
do_clean       = False

# lines
do_split_tav    = False
do_subcont      = False
do_clean_lines  = False


# Nms=select_data(system,
#                 path_calibrated,
#                 cfg,
#                 do_mstransform=do_mstransform)

# continuum(system,
#           Nms=Nms,
#           cfg=cfg,
#           beam=beam,
#           flines=[f12CO21],
#           do_select_cont=do_select_cont,
#           do_fav=do_fav,
#           do_tav=do_tav,
#           do_extractvis=do_extractvis,
#           do_clean=do_clean,
#           nu=2.3e11,
#           robust=robust,
#           dpix=dpix,
#           Npix=Npix,
#           tol=0.5 # 0.5 by default for 3% loss
#           )


# lines(system,
#       Nms=Nms,
#       cfg=cfg,
#       beam=beam,
#       flines=[f12CO21],
#       tags=['12CO'],
#       do_split_tav=do_split_tav,
#       do_subcont=do_subcont,
#       do_clean_lines=do_clean_lines,
#       robust=robust,
#       dpix=dpix,
#       Npix=Npix
#       )

##############################
########## 12mSB ##############
#############################

cfg='12mSB'
path_calibrated='../calibrated/{}/'.format(cfg)

### for 12m 
robust=2.0
Npix=512
dpix=0.1  # arcsec
beam=0.9  # arcsec
uvtaper=''#0.9arcsec'#1.5arcsec'

##################
#### TASKS #######
##################

do_mstransform = False  # remove calibrators and change to BARY reference frame

# continuum
do_select_cont = False  # select continuum by flagging
do_fav         = False  # ferquency average for continuum
do_tav         = False  # time average
do_extractvis  = False
do_clean       = True

# lines
do_split_tav    = False
do_subcont      = False
do_clean_lines  = False


Nms=select_data(system,
                path_calibrated,
                cfg=cfg,
                do_mstransform=do_mstransform)


continuum(system,
          Nms=Nms,
          cfg=cfg,
          beam=beam,
          flines=[f12CO21],
          do_select_cont=do_select_cont,
          do_fav=do_fav,
          do_tav=do_tav,
          do_extractvis=do_extractvis,
          do_clean=do_clean,
          nu=2.3e11,
          robust=robust,
          dpix=dpix,
          Npix=Npix,
          tol=0.5, # 0.5 by default for 3% loss
          uvtaper=uvtaper,
          field='3',
          #gridder='mosaic'
          )

lines(system,
      Nms=Nms,
      cfg=cfg,
      beam=beam,
      do_split_tav=do_split_tav,
      do_subcont=do_subcont,
      do_clean_lines=do_clean_lines,
      robust=robust,
      dpix=dpix,
      Npix=Npix,
      flines=[f12CO21],
      tags=['12CO'],
      )





##############################
# ########## 12mLB ##############
# #############################

cfg='12mLB'
path_calibrated='../calibrated/{}/'.format(cfg)

### for 12m 
robust=2.0
Npix=1024
dpix=0.03  # arcsec
beam=0.2  # arcsec
uvtaper=''#0.9arcsec'#1.5arcsec'

##################
#### TASKS #######
##################

do_mstransform = False  # remove calibrators and change to BARY reference frame

# continuum
do_select_cont = False  # select continuum by flagging
do_fav         = False  # ferquency average for continuum
do_tav         = False  # time average
do_extractvis  = False
do_clean       = False

# lines
do_split_tav    = False
do_subcont      = False
do_clean_lines  = False


Nms=select_data(system,
                path_calibrated,
                cfg=cfg,
                do_mstransform=do_mstransform)

continuum(system,
          Nms=Nms,
          cfg=cfg,
          beam=beam,
          do_select_cont=do_select_cont,
          do_fav=do_fav,
          do_tav=do_tav,
          do_extractvis=do_extractvis,
          do_clean=do_clean,
          flines=[f12CO21],
          nu=2.3e11,
          robust=robust,
          dpix=dpix,
          Npix=Npix,
          tol=0.5, # 0.5 by default for 3% loss
          uvtaper=uvtaper
          )

lines(system,
      Nms=Nms,
      cfg=cfg,
      beam=beam,
      do_split_tav=do_split_tav,
      do_subcont=do_subcont,
      do_clean_lines=do_clean_lines,
      robust=robust,
      dpix=dpix,
      Npix=Npix,
      flines=[f12CO21],
      tags=['12CO'],
      )

"""

### image both data sets (without any corrections)

do_clean_combine = False

### for 12m 
robust=2.0
Npix=512
dpix=0.04  # arcsec
beam=0.4  # arcsec

if do_clean_combine:
    image_name=system['target']+'.12mACA.briggs.%1.1f.%i.%1.2f'%(robust, Npix, dpix)
    os.system('rm -r '+image_name+'*')

    tclean(
            vis            = ['TYC9340.ACA.cal.bary.continuum.fav.tav.ms', 'TYC9340.12m.cal.bary.continuum.fav.tav.ms'],
            imagename      = image_name,
            cell           = '{}arcsec'.format(dpix),
            imsize         = [Npix,Npix],
            specmode           = 'mfs',
            gridder='mosaic',
            weighting      = 'briggs',
            robust         = robust,
            uvtaper='',
            niter          = 10000,
            interactive    = True,
            pblimit=0.1,
            pbcor=True)


    # export fits file 
    fits_image=image_name+'.fits'
    fits_pbcor=image_name+'.pbcor.fits'
    fits_pb=image_name+'.pb.fits'
            
    exportfits( 
        imagename = image_name+'.image',
        fitsimage = fits_image,
        overwrite = True
            )
    exportfits( 
        imagename = image_name+'.image.pbcor',
        fitsimage = fits_pbcor,
        overwrite = True
        )
    exportfits(
        imagename = image_name+'.pb',
        fitsimage = fits_pb,
        overwrite = True
        )
    
    imview(image_name+'.image')
    
"""

