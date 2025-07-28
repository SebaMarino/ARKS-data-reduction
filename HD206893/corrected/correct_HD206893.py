"""
This script was written to run on CASA 6.4.1.12 to apply corrections to the following datasets belonging to project 2022.1.00338.L 

Target: HD197481 (AU Mic) / archival

Archival bang 6 observations

- 12m SB observed with 1 execution block in March 2014
- 12m LB observed with 2 execution blocks in Aug 2014-Jun 2015

reducer: S. Marino

"""
import os, sys
import numpy as np
import json
from astropy import units as u
from astropy.coordinates import Angle

execfile('../../functions_data_reduction.py') # it needs to be run twice for correct import


do_correct=True
do_concat=True
do_extractvis=True
do_clean=False

target='HD206893'
model='model_SMG_best'
cfgs=['12mLB', '12mSB']
obref='12mLB.obs1'
Nms=[4, 2]
obtags=[]


for i in range(len(cfgs)):
    for j in range(Nms[i]):
            obtags.append('{}.obs{}'.format(cfgs[i], j+1))

f=open('MCMC_results.json')
results=json.load(f)
f.close()
## add flux correction of 1 to reference EB
results['fflux-{}'.format(obref)]={'median':1.0}

##############################
# ##### APPLY CORRECTIONS #####
# #############################

if do_correct:
    for i, obi in enumerate(obtags):
        print(i, obi)

#         # with SMG subtraction
        apply_corrections('../reduced/{}.{}.cal.bary.continuum.fav.tav.ms'.format(target, obi),
                          dRA=results['deltaRA-{}'.format(obi)]['median']-results['deltaRA-{}'.format(obref)]['median'],
                          dDec=results['deltaDec-{}'.format(obi)]['median']-results['deltaDec-{}'.format(obref)]['median'],
                          fsigma=results['fsigma-{}'.format(obi)]['median'],
                          flux_offset=results['fflux-{}'.format(obi)]['median'],
                          model_vis='{}.{}.continuum.fav.tav.{}.txt'.format(target, obi, model),
                          extractvis=False,
                          cfg=obi
                         )
        
        
        # without SMG subtraction
        apply_corrections('../reduced/{}.{}.cal.bary.continuum.fav.tav.ms'.format(target, obi),
                          dRA=results['deltaRA-{}'.format(obi)]['median']-results['deltaRA-{}'.format(obref)]['median'],
                          dDec=results['deltaDec-{}'.format(obi)]['median']-results['deltaDec-{}'.format(obref)]['median'],
                          fsigma=results['fsigma-{}'.format(obi)]['median'],
                          flux_offset=results['fflux-{}'.format(obi)]['median'],
                          model_vis='',
                          extractvis=False,
                          cfg=obi)
        
        # no time averaging and no SMG subtraction
        apply_corrections('../reduced/{}.{}.cal.bary.continuum.fav.ms'.format(target, obi),
                          dRA=results['deltaRA-{}'.format(obi)]['median']-results['deltaRA-{}'.format(obref)]['median'],
                          dDec=results['deltaDec-{}'.format(obi)]['median']-results['deltaDec-{}'.format(obref)]['median'],
                          fsigma=results['fsigma-{}'.format(obi)]['median'],
                          flux_offset=results['fflux-{}'.format(obi)]['median'],
                          model_vis='',
                          extractvis=False,
                          cfg=obi)
        

        # apply corrections to lines
        apply_corrections('../reduced/{}.{}.cal.bary.12CO.tav.contsub.ms'.format(target, obi),
                          dRA=results['deltaRA-{}'.format(obi)]['median']-results['deltaRA-{}'.format(obref)]['median'],
                          dDec=results['deltaDec-{}'.format(obi)]['median']-results['deltaDec-{}'.format(obref)]['median'],
                          fsigma=results['fsigma-{}'.format(obi)]['median'],
                          flux_offset=results['fflux-{}'.format(obi)]['median'],
                          model_vis='',
                          extractvis=False,
                          cfg=obi)
        
#         apply_corrections('../reduced/{}.{}.cal.bary.13CO.tav.contsub.ms'.format(target, obi),
#                           dRA=results['deltaRA-{}'.format(obi)]['median']-results['deltaRA-{}'.format(obref)]['median'],
#                           dDec=results['deltaDec-{}'.format(obi)]['median']-results['deltaDec-{}'.format(obref)]['median'],
#                           fsigma=results['fsigma-{}'.format(obi)]['median'],
#                           flux_offset=results['fflux-{}'.format(obi)]['median'],
#                           model_vis='',
#                           extractvis=False,
#                           cfg=obi)

###############################################################################
# ##### concatenate and fix pointing #####
# ########################################
if do_concat:
    sufixes=['continuum.fav.tav.SMGsub.corrected',
             'continuum.fav.tav.corrected',
             'continuum.fav.corrected',
             '12CO.tav.contsub.corrected',
             #'13CO.tav.contsub.corrected'
            ]
    ## get ref direction
    msmd.open('{}.{}.cal.bary.{}.ms'.format(target, obref, sufixes[0]))
    refdir=msmd.phasecenter()
    msmd.close()
    
    ra=Angle(refdir['m0']['value'], refdir['m0']['unit'] ).wrap_at(360 * u.deg)
    dec=Angle(refdir['m1']['value'], refdir['m1']['unit'] )   
    direction='J2000 '+ra.to_string(unit=u.hour, precision=6)+' '+dec.to_string(unit=u.degree, precision=6)
    print(direction, refdir['refer'])
    
    for i in range(len(cfgs)):
        for sufix in sufixes:
            ms_files=['{}.{}.obs{}.cal.bary.{}.ms'.format(target, cfgs[i], j+1, sufix) for j in range(Nms[i])]
            concatenate(ms_files, concatvis='{}.{}.cal.bary.{}.ms'.format(target, cfgs[i], sufix), refdir=direction, reframe=refdir['refer'])
        
if do_extractvis:
    
    for i in range(len(cfgs)):
        
        uvtable_filename='{}.{}.continuum.fav.tav.SMGsub.corrected.txt'.format(target, cfgs[i])
        os.system('rm '+uvtable_filename)
        extract_visibilities('{}.{}.cal.bary.continuum.fav.tav.SMGsub.corrected.ms'.format(target, cfgs[i]), uvtable_filename)

        uvtable_filename='{}.{}.continuum.fav.tav.corrected.txt'.format(target, cfgs[i])
        os.system('rm '+uvtable_filename)
        extract_visibilities('{}.{}.cal.bary.continuum.fav.tav.corrected.ms'.format(target, cfgs[i]), uvtable_filename)

#################################
# ###### DO FINAL CLEANS #########
# ################################


if do_clean:
    
    #"""
    ### clean combined configurations 
    robust=0.5
    Npix=1024
    dpix=0.02  # arcsec
    tag='corrected'#.corrected'
    ms_image=[]
    
    for i in range(len(cfgs)):
        ms_image.append('{}.{}.cal.bary.continuum.fav.tav.{}.ms'.format(target, cfgs[i], tag))
        
    image_name='{}.combined.{}.briggs.{:.1f}.{}.{}'.format(target, tag, robust, Npix, dpix)

    os.system('rm -r '+image_name+'*')

    tclean(vis            = ms_image,
           imagename      = image_name,
           cell           = '{}arcsec'.format(dpix),
           imsize         = [Npix,Npix],
           specmode           = 'mfs',
           #gridder='mosaic',
           #deconvolver='multiscale',
           #scales=[0, 5, 15, 25],
           weighting      = 'briggs',
           robust         = robust,
           uvtaper='',
           niter          = 10000,
           interactive    = True,
           # datacolumn='corrected',
           pbcor=True)
    
    # export fits file 
    fits_image=image_name+'.fits'
    fits_pbcor=image_name+'.pbcor.fits'
    fits_pb=image_name+'.pb.fits'
    fits_model=image_name+'.model.fits'

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
    exportfits(
                imagename = image_name+'.model',
                fitsimage = fits_model,
                overwrite = True
            )
    imview(image_name+'.image')

    ### JvM correction
    do_JvM_correction_and_get_epsilon(image_name)

    fits_image=image_name+'.JvMcorr.fits'
    exportfits( 
                imagename = image_name+'.JvMcorr.image',
                fitsimage = fits_image,
                overwrite = True
                )
    
    
# if do_change_to_ICRS:
    
    
    
#     robusts=[0.5,2.0]
#     Npix=1024
#     dpix=0.04  # arcsec
#     tags=['SMGsub.corrected','corrected']#.corrected'
#     outputs=['.image', '.image.pbcor', '.pb', '.model', '.JvMcorr.image']
#     fitsoutputs=['.fits', '.pbcor.fits', 'pb.fits', '.model.fits', 'JvMcorr.fits']
#     for robust in robusts:
#         for tag in tags:
#             image_name0='{}.combined.{}.briggs.{:.1f}.{}.{}'.format(target, tag, robust, Npix, dpix)
            
#             for i, output in enumerate(outputs):
#                 image_name=image_name0+output
#                 # change to ICRS
#                 change_to_ICRS_image(image_name)
#                 # export to fits again
#                 exportfits(  
#                     imagename = image_name,
#                     fitsimage = image_name0+fitsoutputs[i],
#                     overwrite = True
#                 )
    
    
    
    
    
    
    ## do clean 12m only 

    """
    
    ### clean single configurations 
    #robust=2.0
    #Npix=256
    #dpix=0.3  # arcsec
    #tag='SMGsub.corrected'
    #cfg='ACA'
    
    robust=2.0
    Npix=1024
    dpix=0.04  # arcsec
    tag='SMGsub.corrected'
    cfg='12m'
        
    ms_image='{}.{}.cal.bary.continuum.fav.tav.{}.ms'.format(target, cfg, tag)
    
    image_name='{}.{}.{}.briggs.{:.1f}.{}.{}'.format(target, cfg, tag,robust, Npix, dpix)
    
    os.system('rm -r '+image_name+'*')


    tclean(vis            = ms_image,
           imagename      = image_name,
           cell           = '{}arcsec'.format(dpix),
           imsize         = [Npix,Npix],
           specmode           = 'mfs',
           #gridder='mosaic',
           deconvolver='multiscale',
           scales=[0, 5, 15, 25],
           weighting      = 'briggs',
           robust         = robust,
           uvtaper='',
           niter          = 10000,
           interactive    = True,
           # datacolumn='corrected',
           pbcor=True)


    # export fits file 
    fits_image=image_name+'.fits'
    fits_pbcor=image_name+'.pbcor.fits'
    fits_pb=image_name+'.pb.fits'
    fits_model=image_name+'.model.fits'

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
    exportfits(
                imagename = image_name+'.model',
                fitsimage = fits_model,
                overwrite = True
            )
    
    imview(image_name+'.image')

    ### JvM correction
    do_JvM_correction_and_get_epsilon(image_name)

    fits_image=image_name+'.JvMcorr.fits'
    exportfits( 
                imagename = image_name+'.JvMcorr.image',
                fitsimage = fits_image,
                overwrite = True
                )
    
    """
