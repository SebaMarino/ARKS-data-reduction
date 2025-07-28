"""

This script was written to run on CASA 6.4.1.12 to reduce the following datasets belonging to project 2022.1.00338.L 

Target: HD107146

band 7:
- ACA observed with 4 executions in 2016 and 2017.
- 12m c40-3 observed with 1 execution blocks in Dec 2016
- 12m c40-2 observed with 2 execution blocks in Dec 2019
- 12m c40-5 observed with 1 execution blocks in May 2021

band6:

- 12m observed with 5 execution in Apr 2017

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
system_1={'target'      : 'HD107146', # used to name files
        'field_name'  : 'hd107146', # for ACA, field name in ms
        'Radius'      : 16.4, # arcsec, maximum radius of interest used to calculate maximum freq and time averaging. Set to 500 au or at least 2x disc outer edge or FWHM of Antenna primary beam
        'rv'       : 1.76, # km/s, radial velocity in barycentric reference frame from Gaia DR3
        'dvel'        : 50.0,      # km/s, line width. It is used to flag channels at +-dvel and to image a cube with a width of +-3dvel
        }

system_2={'target'      : 'HD107146', # used to name files
        'field_name'  : 'HD_107146', # for ACA, field name in ms
        'Radius'      : 16.4, # arcsec, maximum radius of interest used to calculate maximum freq and time averaging. Set to 500 au or at least 2x disc outer edge or FWHM of Antenna primary beam
        'rv'       : 1.76, # km/s, radial velocity in barycentric reference frame from Gaia DR3
        'dvel'        : 50.0,      # km/s, line width. It is used to flag channels at +-dvel and to image a cube with a width of +-3dvel
        }

##############################
# ########## ACA ##############
# #############################

band='band7'
cfg='b7.ACA'
path_calibrated='../calibrated/{}/{}/'.format(band,cfg)

### for ACA 
robust=2.0
Npix=256
dpix=0.3  # arcsec
beam=3.0  # arcsec

# ##################
# #### TASKS #######
# ##################

do_mstransform = False  # remove calibrators and change to BARY reference frame

# continuum
do_select_cont = False   # select continuum by flagging
do_fav         = False   # ferquency average for continuum
do_tav         = False   # time average
do_extractvis  = False
do_clean       = False

# lines
do_split_tav    = False
do_subcont      = False
do_clean_lines  = False


Nms=select_data(system_1,
                path_calibrated,
                cfg,
                do_mstransform=do_mstransform)

continuum(system_1,
          Nms=Nms,
          cfg=cfg,
          beam=beam,
          flines=[f12CO32, f13CO32],
          do_select_cont=do_select_cont,
          do_fav=do_fav,
          do_tav=do_tav,
          do_extractvis=do_extractvis,
          do_clean=do_clean,
          nu=3.4e11,
          robust=robust,
          dpix=dpix,
          Npix=Npix,
          tol=0.5 # 0.5 by default for 3% loss
          )


lines(system_1,
      Nms=Nms,
      cfg=cfg,
      beam=beam,
      flines=[f12CO32, f13CO32],
      tags=['12CO','13CO'],
      do_split_tav=do_split_tav,
      do_subcont=do_subcont,
      do_clean_lines=do_clean_lines,
      robust=robust,
      dpix=dpix,
      Npix=Npix
      )

##############################
# ########## 12mSB ##############
# #############################

band='band7'
cfg='b7.12mSB'
path_calibrated='../calibrated/{}/{}/'.format(band,cfg)


### for 12m 
robust=2.0
Npix=1024
dpix=0.08  # arcsec
beam=0.6  # arcsec
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
do_clean       = False

# lines
do_split_tav    = False
do_subcont      = False
do_clean_lines  = False


Nms=select_data(system_2,
                path_calibrated,
                cfg=cfg,
                do_mstransform=do_mstransform)
if do_flag: # flag short baselines in obs 2 ro remove fringes.
    msf='{}.{}.obs1.cal.bary.ms'.format(system_2['target'],cfg)
    flagdata(msf, uvrange='<20klambda', flagbackup=False)
    msf='{}.{}.obs2.cal.bary.ms'.format(system_2['target'],cfg)
    flagdata(msf, uvrange='<20klambda', flagbackup=False)


continuum(system_2,
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
          uvtaper=uvtaper
          )

lines(system_2,
      Nms=Nms,
      cfg=cfg,
      beam=beam,
      do_split_tav=do_split_tav,
      do_subcont=do_subcont,
      do_clean_lines=do_clean_lines,
      robust=robust,
      dpix=dpix,
      Npix=Npix
      )



##############################
########## 12mMB ##############
#############################

band='band7'
cfg='b7.12mMB'
path_calibrated='../calibrated/{}/{}/'.format(band,cfg)


### for 12m 
robust=2.0
Npix=1024
dpix=0.04  # arcsec
beam=0.4  # arcsec
uvtaper=''#0.9arcsec'#1.5arcsec'

##################
#### TASKS #######
##################

do_mstransform = False # remove calibrators and change to BARY reference frame

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


Nms=select_data(system_1,
                path_calibrated,
                cfg=cfg,
                do_mstransform=do_mstransform)

continuum(system_1,
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
          uvtaper=uvtaper
          )

lines(system_1,
      Nms=Nms,
      cfg=cfg,
      beam=beam,
      do_split_tav=do_split_tav,
      do_subcont=do_subcont,
      do_clean_lines=do_clean_lines,
      robust=robust,
      dpix=dpix,
      Npix=Npix
      )



##############################
# ########## 12mLB ##############
# #############################

band='band7'
cfg='b7.12mLB'
path_calibrated='../calibrated/{}/{}/'.format(band,cfg)


### for 12m 
robust=2.0
Npix=2048
dpix=0.02  # arcsec
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


Nms=select_data(system_2,
                path_calibrated,
                cfg=cfg,
                do_mstransform=do_mstransform)

continuum(system_2,
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
          uvtaper=uvtaper
          )

lines(system_2,
      Nms=Nms,
      cfg=cfg,
      beam=beam,
      do_split_tav=do_split_tav,
      do_subcont=do_subcont,
      do_clean_lines=do_clean_lines,
      robust=robust,
      dpix=dpix,
      Npix=Npix
      )

# ##############################
# ########## band 6 12m ##############
# #############################

band='band6'
cfg='b6.12m'
path_calibrated='../calibrated/{}/{}/'.format(band,cfg)


### for 12m 
robust=2.0
Npix=1024
dpix=0.05  # arcsec
beam=0.7  # arcsec
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


Nms=select_data(system_2,
                path_calibrated,
                cfg=cfg,
                do_mstransform=do_mstransform)

continuum(system_2,
          Nms=Nms,
          cfg=cfg,
          beam=beam,
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
          uvtaper=uvtaper
          )

lines(system_2,
      Nms=Nms,
      cfg=cfg,
      beam=beam,
      do_split_tav=do_split_tav,
      do_subcont=do_subcont,
      do_clean_lines=do_clean_lines,
      robust=robust,
      dpix=dpix,
      Npix=Npix
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

