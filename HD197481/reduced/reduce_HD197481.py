"""

This script was written to run on CASA 6.4.1.12 to reduce the following datasets belonging to project 2022.1.00338.L 

Target: HD197481 (AU Mic) / archival

Archival bang 6 observations

- 12m SB observed with 1 execution block in March 2014
- 12m LB observed with 2 execution blocks in Aug 2014-Jun 2015

reducer: S. Marino

###### THIS SYSTEM REQUIRES FLAGGING FLARES ####

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
system={'target'      : 'HD197481', # used to name files
        'field_name'  : 'AU_Mic', # field name in ms
        #'Radius'      : 16.4, # arcsec, maximum radius of interest used to calculate maximum freq and time averaging. Set to FWHM of Antenna primary beam
        'Radius'      : 16.5, # arcsec, maximum radius of interest used to calculate maximum freq and time averaging. Set to 500 au or at least 2x disc outer edge or FWHM of Antenna primary beam
        'rv'       : -6.9, # km/s, radial velocity in barycentric reference frame from Gaia DR3
        'dvel'        : 50.0,      # km/s, line width. It is used to flag channels at +-dvel and to image a cube with a width of +-3dvel
        }


# #############################
# ########## 12mSB ##############
# #############################

cfg='12mSB'
path_calibrated='../calibrated/{}/'.format(cfg)

### for 12m SB
robust=2.0
Npix=512
dpix=0.08  # arcsec
beam=0.8 # arcsec
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
          nu=3.4e11,
          robust=robust,
          dpix=dpix,
          Npix=Npix,
          tol=0.5, # 0.5 by default for 3% loss
          uvtaper=uvtaper,
          flines=[f12CO21, f13CO21]
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
      flines=[f12CO21, f13CO21]
      )

# #############################
# ########## 12mLB ##############
# #############################

cfg='12mLB'
path_calibrated='../calibrated/{}/'.format(cfg)

### for 12m SB
robust=2.0
Npix=1024
dpix=0.02  # arcsec
beam=0.25  # arcsec
uvtaper=''#0.9arcsec'#1.5arcsec'

##################
#### TASKS #######
##################

do_mstransform = False  # remove calibrators and change to BARY reference frame
do_flag        = False
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

if do_flag: # flag scan 27 in obs 2 ro remove flare.
    msf='{}.{}.obs2.cal.bary.ms'.format(system['target'],cfg)
    flagdata(msf, scan='27', flagbackup=False)

continuum(system,
          Nms=Nms,
          cfg=cfg,
          beam=beam,
          do_select_cont=do_select_cont,
          do_fav=do_fav,
          do_tav=do_tav,
          do_extractvis=do_extractvis,
          do_clean=do_clean,
          nu=3.4e11,
          robust=robust,
          dpix=dpix,
          Npix=Npix,
          tol=0.5, # 0.5 by default for 3% loss
          uvtaper=uvtaper,
          flines=[f12CO21, f13CO21]

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
      flines=[f12CO21, f13CO21]
      )


#"""

### image both data sets (without any corrections)

do_clean_combine = False

### for 12m 
robust=0.5
Npix=1024
dpix=0.02  # arcsec
beam=0.4  # arcsec

if do_clean_combine:
    image_name=system['target']+'.12m.briggs.%1.1f.%i.%1.2f'%(robust, Npix, dpix)
    os.system('rm -r '+image_name+'*')

    tclean(
            vis            = [system['target']+'.12mSB.obs1.cal.bary.continuum.fav.tav.ms',
                              system['target']+'.12mLB.obs1.cal.bary.continuum.fav.tav.ms',
                              system['target']+'.12mLB.obs2.cal.bary.continuum.fav.tav.ms',
                              system['target']+'.12mLB.obs3.cal.bary.continuum.fav.tav.ms'],
            imagename      = image_name,
            cell           = '{}arcsec'.format(dpix),
            imsize         = [Npix,Npix],
            specmode           = 'mfs',
            #gridder='mosaic',
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
    
#"""

