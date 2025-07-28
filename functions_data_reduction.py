"""
Functions used for the data reduction of ARKS data.

Contributors: S. Marino


"""

import os, sys
import numpy as np
import shutil

clight_kms= 2.99792458e+5 # [km/s] Speed of light
clight_ms = 2.99792458e+8 # [m/s] Speed of light


f12CO32=3.4579599e11    # Hz
f13CO32=3.30587965e11 # Hz

f12CO21=2.3053800000e11 # Hz
f13CO21=2.2039868420e11    # Hz


#################################################################
############### MISCELLANEOUS FUNCTIONS  #########################
##################################################################

def get_spwids(msfile, field): 
    """
    Function to obtain list of valid spectral windows in an ms for the specified field

    Parameters
    ==========
    msfile, string, path to ms file
    field: int or strong, field id or name 
    Returns
    =======
    string with spw IDs separated by commas
    """
    msmd.open(msfile)
    spwIDs=msmd.spwsforfield(field)    
    Nspws=len(spwIDs)

    ret_spw=''
    
    if Nspws<=4:
         for ispw in spwIDs:
            ret_spw+=', {}'.format(ispw)
        
    else: ## if number>4 in calibrated ms of single EB, then select only those with more than 1 channel or those with FULL_RES in their name. This is necessary for ACA ms' of HR8799 since calibrated ms contained spw's repeated that spwsforfield returned. 
        for ispw in spwIDs:
            chan_freqs = msmd.chanfreqs(ispw)
            nchan = len(chan_freqs)
            if nchan>1:
                ret_spw+=', {}'.format(ispw)  
               
    msmd.close()

    return ret_spw[2:]

def get_linefreq(flines=[], dvel=50.0, rv=0.0):
    """
    Function to obtain frquencies where line emission could be present and that can be passed to spw argument.

    Parameters
    ==========
    flines: array of floats, line frequencies in Hz
    dvel: float, 1/2 of line width in km/s
    rv:   float, system radial velocity in same frame as the ms file
    
    Returns
    =======
    string with the range of frequencies of a line in format that can be passed to CASA functions as spw 
    """

    channels=''
    for fline in flines:    
        
        flinei=fline*(1 - rv/clight_kms - dvel/clight_kms) # Hz
        flinef=fline*(1 - rv/clight_kms + dvel/clight_kms) # Hz

        channels+='*:%1.3f~%1.3fMHz,'%(flinei/1.0e6, flinef/1.0e6)

    return channels[:-1]

def get_linechannels(msfile, flines=[], dvel=50.0, rv=0.0):
    """
    Function to obtain spws and channels where line emission could be present 

    Parameters
    ==========
    msfile, string, path to ms file
    flines: array of floats, line frequencies in Hz
    dvel: float, 1/2 of line width in km/s
    rv:   float, system radial velocity in same frame as the ms file
    
    Returns
    =======
    string with the spw's and channels where line emission may be present in format that can be passed to CASA functions as spw 

    """
    
    msmd.open(msfile)
    Nspws=msmd.nspw()

    channels=''

    
    for ispw in range(Nspws):
        chan_freqs = msmd.chanfreqs(ispw)
        nchan = len(chan_freqs)
        dnus  = msmd.chanwidths(ispw)

        fi=chan_freqs[0]
        ff=chan_freqs[-1]

        ## check if spw contain lines
        for fline in flines:

            ## correct line frequency due to RV
            fline=fline*(1 - rv/clight_kms) 
    
            if fi<fline<ff or ff<fline<fi:

                ## search for channels that contain line

                flinei=fline*(1 - dvel/clight_kms)
                flinef=fline*(1 + dvel/clight_kms)
                
                chani = int((flinei-fi)/dnus[0])
                chanf = int((flinef-fi)/dnus[0])

                if chani<chanf:
                    channels+='%i:%i~%i,'%(ispw, chani, chanf)
                else:
                    channels+='%i:%i~%i,'%(ispw, chanf, chani)
    msmd.close()

    return channels[:-1]


def get_spws_line(msfile, fline, rv=0.0):
    """
    Function to find spws where a single emission line could be present 

    Parameters
    ==========
    msfile, string, path to ms file
    flines: array of floats, line frequencies in Hz
    rv:   float, system radial velocity in same frame as the ms file
    
    Returns
    =======
    string with list of spw ID's that overlay with an emission line
    """
    
    msmd.open(msfile)
    Nspws=msmd.nspw()


    fline=fline*(1 - rv/clight_kms) 
    spws=''

    for ispw in range(Nspws):
        chan_freqs = msmd.chanfreqs(ispw)
        nchan = len(chan_freqs)
        dnus  = msmd.chanwidths(ispw)

        fi=chan_freqs[0]
        ff=chan_freqs[-1]

        ## check if spw contain lines
        if fi<fline<ff or ff<fline<fi:
            spws+='%i,'%ispw
    msmd.close()

    
    return spws[:-1]
                    

                    
def get_av_widths(msfile, max_freq_av=1.0e9):
    """
    function to find the widths that one needs to pass to split to do frequency averaging, such that all resulting spws have the same number of channels (important for exporting visibilities)
    
    Parameters
    ==========
    msfile, string, path to ms file
    max_freq_av: float, maximum amount of frequency averaging in Hz

    Returns
    =======
    ndarray with widths that should be passed to split to obtain a final frequency-averaged ms with spw's with the same number of channels.  
    """
    
    msmd.open(msfile)
    Nspws=msmd.nspw()
  
    bandwidths=[]
    nchans=[]
    for ispw in range(Nspws):

        chan_freqs = msmd.chanfreqs(ispw)
        nchan = len(chan_freqs)
        dnus  = msmd.chanwidths(ispw)

        bandwidths.append(np.abs(np.sum(dnus)))
        nchans.append(nchan)

    msmd.close()

    bandwidths=np.array(bandwidths)
    nchans=np.array(nchans)
    
    # find the largest bandwidth
    imax=np.argmax(bandwidths)
    print('maximum bandwidth of spws is {} GHz'.format(bandwidths[imax]/1.0e9))
    # calculate the number of channels (multiples of 2) after averaging
    for i in range(10):
        if bandwidths[imax]/2**(i)<1.1*max_freq_av:
            nchan_final=2**(i)
            break
    print('Will average channels to leave {} per spw'.format(nchan_final))
    ## calculate final widths
    widths=nchans/nchan_final
    
    return widths.astype(int)

def get_channel_widths(msfile):
    """
    function to find the channel frequency widths 
    
 
    Parameters
    ==========
    msfile, string, path to ms file

 
    Returns
    =======
    ndarray with the widths in Hz.
    """
    
    msmd.open(msfile)
    Nspws=msmd.nspw()
  
    fwidths=[]
    nchans=[]
    for ispw in range(Nspws):

        chan_freqs = msmd.chanfreqs(ispw)
        nchan = len(chan_freqs)
        dnus  = msmd.chanwidths(ispw)

        fwidths.append(np.abs(dnus[0]))

    msmd.close()
    
    return np.array(fwidths)

def get_fields(msfile, field_name):
    """
    function to find fields with data in spw. Useful to split mosaics into different pointings
 
    Parameters
    ==========
    msfile, string, path to ms file

 
    Returns
    =======
    ndarray with field IDs.
    """
    msmd.open(msfile)
    
    fieldids=msmd.fieldsforname(field_name)
    
    fields=[]
    for i in fieldids:
        
        if len(msmd.intentsforfield(i))>0:
            fields.append(i)
       
    msmd.close()
    
    return(fields)

def extract_visibilities(msfile, uvtable_filename, datacolumn='data', verbose=True, dualpol=True, fmt='%10.6e'):
    """
    function to extract visibilities and save them in a table.
    Adapted by S. Marino from uvplot package by M. Tazzari to put uv points in lambda, rather than m.

    Parameters
    ==========
    msfile, string, path to ms file
    uvtable_filename, string, name given to the output table
    datacolumn, string, column to extract (e.g. data or corrected)
    verbose, bool, whether to print some information
    dualpol, bool, whether to extract dual pollarization data
    fmt, fmt, format in which to extract the data.

    """

    tb.open(msfile)

    ## get coordinates (m)
    uvw = tb.getcol("UVW")
    u, v, w = [uvw[i, :] for i in range(3)]

    # get weights
    weights_orig = tb.getcol("WEIGHT")

    # get visibilities
    tb_columns = tb.colnames()
    #print(tb_columns)
    if datacolumn.upper() in tb_columns:
        data = tb.getcol(datacolumn.upper())
    else:
        raise KeyError("datacolumn {} is not available.".format(datacolumn))
        

    spw = tb.getcol("DATA_DESC_ID")
    nspw = len(np.unique(spw))


    #print(spw)
    #print(spw.shape)
    #print(u.shape)
    #print(data.shape)
        
    nchan = data.shape[1]

    # repeat u and v nchan times
    
    u = np.tile(u, nchan)
    v = np.tile(v, nchan)

    
    if verbose:
        print("Exporting {} channels for {} spws.".format(nchan, nspw))


    if dualpol:
        # dual polarisation: extract the polarised visibilities and weights
        V_XX = data[0, :, :].reshape(-1)
        V_YY = data[1, :, :].reshape(-1)
        weights_XX = weights_orig[0, :]
        weights_YY = weights_orig[1, :]
        if nchan > 1:
            weights_XX = np.tile(weights_XX, nchan)
            weights_YY = np.tile(weights_YY, nchan)

        # compute weighted average of the visibilities and weights
        V = (V_XX * weights_XX + V_YY * weights_YY) / (weights_XX + weights_YY)
        weights = weights_XX + weights_YY

    else:
        # single polarisation
        V = data[0, :, :].reshape(-1)
        weights = weights_orig
        if nchan > 1:
            weights = np.tile(weights, nchan)

    spw_path = tb.getkeyword('SPECTRAL_WINDOW').split()[-1]
    tb.close()

    # get the mean observing frequency
    tb.open(spw_path)
    freqs = tb.getcol('CHAN_FREQ')  # [Hz]
    tb.close()
    wle = clight_ms / freqs.mean()  # [m]

    ### convert u,v in m to lambda
    lams=np.zeros_like(data[0,:,:], dtype=float) # channel x spw
    for ichan in range(nchan):
        for irow in range(len(spw)):
            lams[ichan, irow] = clight_ms/(freqs[ichan, spw[irow]]) # m
    lams=lams.reshape(-1)
    u=u/lams
    v=v/lams


    # export to file as ascii
    if verbose:
        print("Exporting visibilities to {}...".format(uvtable_filename), end='')

    np.savetxt(uvtable_filename,
               np.column_stack([u, v, V.real, V.imag, weights, lams]),
               fmt=fmt, delimiter='\t',
               header='Extracted from {}.\nwavelength[m] = {}\nColumns:\tu[lambda]\tv[lambda]\tRe(V)[Jy]\tIm(V)[Jy]\tweight\tlambda[m]'.format(
                   msfile , wle))

    if verbose:
        print("done.")

### it might be useful to image residuals to have a function to subtract or put back visibility table into ms


################################################################
############### FUNCTIONS FOR DATA REDUCTION ###################
################################################################
        
def select_data(system, path_calibrated='', cfg='',  do_concat=True, do_mstransform=True):
    """
    function to select and mstransform individual calibrated observations (EB's) produced by the ALMA pipeline
    
    Parameters
    ==========
    system: dictionary, it contains system information to name files
    path_calibrated: string, path to directory with calibrated observations 
    cfg: str, tag to identify these observations with an alma antenna configuration (ACA, 12m, 12mSB, 12mLB, etc)
    do_mstransform: boolm whether to remove calibrators from observations and transform ms to the Barycentric reference frame (skip if already done)
    """
 

    # find files
    path_files=[path_calibrated+f for f in os.listdir(path_calibrated) if ( (f.endswith('split.cal') | f.endswith('.ms'))  & (f[:3]=='uid') ) ]
    Nms=len(path_files)
    
    for i in range(Nms):
        
        ms0=path_files[i]
        ms1='{}.{}.obs{}.cal.bary.ms'.format(system['target'],cfg,i+1)
        print(ms1)
        
        if do_mstransform:
            os.system('rm -rf '+ms1)
            
            spws=get_spwids(ms0, system['field_name'])
            print(ms0)
            print(spws)
            mstransform(vis=ms0,
                        outputvis=ms1,
                        field=system['field_name'],
                        datacolumn='data',
                        spw=spws,
                        regridms=True,
                        mode='channel',
                        outframe='BARY',
                        keepflags=False)
        
            listobs(ms1, listfile=ms1[:-3]+'.list', overwrite=True)
            
    return Nms

        
    
def continuum(system, Nms=1, cfg='', beam=1.0, flines=[f12CO32, f13CO32], do_select_cont=True, do_fav=True, do_tav=True, do_extractvis=False, do_clean=False, clean_individual=False, tol=0.5, nu=3.4e11, robust=0.5, dpix=0.1, Npix=512, uvtaper='', mosaic=False, **tclean_kwargs, ):

    """
    function to produce continuum ms with lines flagged, frequency averaged, time averaged, extract visibilities, and produce clean image
     

    Parameters
    ==========
    system: dictionary, it contains system information
    cfg: str, tag to identify these observations with an alma antenna configuration (ACA, 12m, 12mSB, 12mLB, etc)
    beam: float, approximate beam size in arcsec. This is used to calculate maximum amount of frequency and time averaging
    flines: array, frequency of lines in Hz to be flagged
    do_select_cont: bool, 
    do_fav: bool, whether to frequency average data
    do_tav: bool, whether to time average the data
    do_extractvis: bool, whether to extract the visibilities
    do_clean: bool, whether to clean continuum
    tol: float, it represents the maximum amount of averaging allowed. A value of 0.2 results in a loss of about than 1 per cent whereas a value of 0.5 results in a loss of about 5%. 
    """
    
    ms_image=[]
    
    
    
    for i in range(Nms):

        ms1  ='{}.{}.obs{}.cal.bary.ms'.format(system['target'],cfg, i+1)
        mscon='{}.{}.obs{}.cal.bary.continuum.ms'.format(system['target'], cfg, i+1)
        msconfav='{}.{}.obs{}.cal.bary.continuum.fav.ms'.format(system['target'], cfg, i+1)
        msconfavtav='{}.{}.obs{}.cal.bary.continuum.fav.tav.ms'.format(system['target'], cfg, i+1)
        uvtable_filename='{}.{}.obs{}.continuum.fav.tav.txt'.format(system['target'], cfg, i+1)

        ms_image.append(msconfavtav)

        if do_fav or do_tav or do_extractvis:

            fields=get_fields(ms1, system['field_name'])
            Nfields=len(fields)
        
            ## if multiple fields, split fav tav ms for modelling
            if Nfields>1 and mosaic:
                msconfavs=[]
                msconfavtavs=[]
                uvtable_filenames=[]
                for fieldi in fields:
                    msconfavs.append('{}.{}.obs{}.field{}.cal.bary.continuum.fav.ms'.format(system['target'], cfg, i+1, fieldi))
                    msconfavtavs.append('{}.{}.obs{}.field{}.cal.bary.continuum.fav.tav.ms'.format(system['target'], cfg, i+1, fieldi))
                    uvtable_filenames.append('{}.{}.obs{}.field{}.continuum.fav.tav.txt'.format(system['target'], cfg, i+1, fieldi))        
        
                        
        if do_select_cont:

            # identify spws that coincide with the 12CO or 13CO 3-2 lines. 
            spw_line=get_linefreq(flines=flines, dvel=system['dvel'], rv=system['rv'])
            # remove previous version
            os.system('rm -rf '+mscon)
            # copy ms before flagging
            os.system('cp -r '+ms1 +' '+mscon)
            ### flagging below if needed
            flagdata(mscon, spw=spw_line, flagbackup=False)
    
        if do_fav:

            ### calculate amount of averaging to avoid bandwidth and time-averaging smearing
            dnu=beam/system['Radius']*nu * tol # Hz
            max_freq_av=dnu # Hz

            print('Max frequency averaging: %1.2f GHz'%(max_freq_av/1.0e9))
            width_cont=get_av_widths(mscon, max_freq_av=max_freq_av)
            print('widths used for averaging channels ', width_cont)
            os.system('rm -rf '+msconfav)
            split(vis=mscon, outputvis=msconfav, width=width_cont, datacolumn='data', keepflags=False)
            listobs(msconfav, listfile=msconfav[:-3]+'.list', overwrite=True)
            
            if Nfields>1 and mosaic:
                for j in range(Nfields):
                    os.system('rm -rf '+msconfavs[j])
                    split(vis=msconfav, outputvis=msconfavs[j], field=fields[j], datacolumn='data')
                    listobs(msconfavs[j], listfile=msconfavs[j][:-3]+'.list', overwrite=True) 

        if do_tav:

            dt= tol * beam/system['Radius'] * 86400.0 /(2*np.pi) # s
            timebin='%1.0fs'%(min(dt,60))
            print('Max time averaging: '+timebin)
 
            os.system('rm -rf '+msconfavtav)
            split(vis=msconfav, outputvis=msconfavtav, timebin=timebin, datacolumn='data', field='')
            listobs(msconfavtav, listfile=msconfavtav[:-3]+'.list', overwrite=True)

            if Nfields>1 and mosaic:
                for j in range(Nfields):
                    os.system('rm -rf '+msconfavtavs[j])
                    split(vis=msconfavtav, outputvis=msconfavtavs[j], field=fields[j], datacolumn='data')
                    listobs(msconfavtavs[j], listfile=msconfavtavs[j][:-3]+'.list', overwrite=True) 

        if do_extractvis:

            if Nfields>1 and mosaic:
                for j in range(Nfields):
                    os.system('rm '+uvtable_filenames[j])
                    extract_visibilities(msconfavtavs[j], uvtable_filenames[j])
                    
            else:
                os.system('rm '+uvtable_filename)
                extract_visibilities(msconfavtav, uvtable_filename)
            
    
        if clean_individual:
            
            if uvtaper=='':
                image_name='{}.{}.obs{}.cal.bary.continuum.fav.tav.briggs.{:.1f}.{:d}.{:.3f}'.format(system['target'],cfg, i+1, robust, Npix, dpix)
            else:
                image_name='{}.{}.obs{}.cal.bary.continuum.fav.tav.briggs.{:.1f}.{}.{:d}.{:.3f}'.format(system['target'],cfg, i+1, robust, uvtaper, Npix, dpix)
            os.system('rm -r '+image_name+'*')

            tclean(
                vis            = msconfavtav,
                imagename      = image_name,
                cell           = '{}arcsec'.format(dpix),
                imsize         = [Npix,Npix],
                specmode           = 'mfs',
                # gridder=gridder,
                weighting      = 'briggs',
                robust         = robust,
                uvtaper=uvtaper,
                niter          = 10000,
                interactive    = True,
                pbcor=True,
                **tclean_kwargs
            )


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
            
    if do_clean:
        
        if uvtaper=='':
            image_name='{}.{}.cal.bary.continuum.fav.tav.briggs.{:.1f}.{:d}.{:.3f}'.format(system['target'],cfg, robust, Npix, dpix)
        else:
            image_name='{}.{}.cal.bary.continuum.fav.tav.briggs.{:.1f}.{}.{:d}.{:.3f}'.format(system['target'],cfg, robust, uvtaper, Npix, dpix)
        os.system('rm -r '+image_name+'*')

        print(tclean_kwargs)
        tclean(
            vis            = ms_image,
            imagename      = image_name,
            cell           = '{}arcsec'.format(dpix),
            imsize         = [Npix,Npix],
            specmode           = 'mfs',
            # gridder=gridder,
            weighting      = 'briggs',
            robust         = robust,
            uvtaper=uvtaper,
            niter          = 10000,
            interactive    = True,
            pblimit        = 0.1,
            pbcor=True,
            **tclean_kwargs
        )


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


def lines(system, Nms=1, cfg='', flines=[f12CO32, f13CO32], tags=['12CO', '13CO'], beam=1.0, do_split_tav=True, do_subcont=True, do_clean_lines=True, tol=0.5, robust=0.5, Npix=512, dpix=0.1, width_kms=0., mosaic=False):
    """
    function to produce line ms with continuum subtracted, time averaged and produce clean cubes


    Parameters
    ==========
    system: dictionary, it contains system information
    cfg: str, tag to identify these observations with an alma antenna configuration (ACA, 12m, 12mSB, 12mLB, etc)
    flines: array, frequency of lines in Hz to be split and imaged
    tags: array of strings, tags of the lines to name the files
    do_split_tav: bool, whether to split and time average the spw's with the individual lines in separate files
    do_subcont: bool, whether to subtract the continuum
    do_clean_lines: bool, whether whether to clean and produce line cubes.
    width_kms: float, channel width in km/s for line imaging. If =0 it will use the channel width of the spw.
    """

    ms_image=[]
    
    for i in range(Nms):
    
        ms1  ='{}.{}.obs{}.cal.bary.ms'.format(system['target'],cfg, i+1)

        if do_subcont:
            fields=get_fields(ms1, system['field_name'])
            Nfields=len(fields)
        
         
        if do_split_tav:


            dt= tol * beam/system['Radius'] * 8.6400e5 /(2*np.pi) # s
            timebin='%1.0fs'%(min(dt,60))
            print('Max time averaging: '+timebin)
        
            for i,linei in enumerate(flines):
                # identify spw's
                spws=get_spws_line(ms1, linei, rv=system['rv'])
                if spws!='':
                    mslinei=ms1[:-3]+'.'+tags[i]+'.tav.ms'
                    os.system('rm -r '+mslinei)
                    split(vis=ms1, outputvis=mslinei, spw=spws, datacolumn='data', timebin=timebin)


        if do_subcont:
    
            for i,linei in enumerate(flines):

                # find spw and channels containing line
                spw_line=get_linefreq(flines=[linei], dvel=system['dvel'], rv=system['rv'])
                print('Excluding the region '+spw_line+' for fitting '+tags[i]+' line at %1.3f GHz'%(linei/1.0e9))
        
                mslinei=ms1[:-3]+'.'+tags[i]+'.tav.ms'
            
                ## if multiple fields, split fav tav ms for modelling
                if Nfields>1 and mosaic:
                    mslinecontsubs=[]
                    for fieldi in fields:
                        mslinecontsubs.append('{}.{}.obs{}.field{}.cal.bary.'.format(system['target'], cfg, i+1, fieldi)+tags[i]+'.tav.contsub.ms')
       
                if os.path.isdir(mslinei):
                    os.system('rm -r '+mslinei+'.contsub')
                    os.system('rm -r '+mslinei[:-3]+'.contsub.ms')

                    uvcontsub(vis=mslinei, fitspw=spw_line, excludechans=True, fitorder=1)
                    # change name
                    os.system('mv {}.contsub {}.contsub.ms'.format(mslinei, mslinei[:-3]))

                    if Nfields>1 and mosaic:
                        for j in range(Nfields):
                            os.system('rm -rf '+mslinecontsubs[j])
                            split(vis=mslinei, outputvis=mslinecontsubs[j], field=fields[j], datacolumn='data')
                            listobs(mslinecontsubs[j], listfile=mslinecontsubs[j][:-3]+'.list', overwrite=True) 
                
    if do_clean_lines:
        
        
        for i,linei in enumerate(flines):
            
            ms_image=[]
            for ims in range(Nms):
                mslinei  ='{}.{}.obs{}.cal.bary.{}.tav.contsub.ms'.format(system['target'],cfg, ims+1, tags[i] )
                if os.path.isdir(mslinei):
                    ms_image.append(mslinei)
            
            if ms_image==[]:
                continue
            image_name='{}.{}.cal.bary.{}'.format(system['target'],cfg, tags[i])
            image_name+='.briggs.%1.1f.%i.%1.3f'%(robust, Npix, dpix)
        
            os.system('rm -r '+image_name+'*')

            ## determine number of channels necessary
            if width_kms==0.0:
                # use channel width
                width_inp=''
                width0=get_channel_widths(ms_image[0])[0] # channel width of first spw in Hz. It assumes all spw have the same channel widths
                width_kms=width0*clight_kms/linei
                
            width_inp='{}km/s'.format(width_kms)
            vi=system['rv']-system['dvel']
            dv=2*system['dvel']
            nchan=int(dv/width_kms)
            print('nchan: {}, vi: {} km/s'.format(nchan, vi))
        
            tclean( vis            = ms_image,
                    imagename      = image_name,
                    cell           = '{}arcsec'.format(dpix),
                    imsize         = [Npix,Npix],
                    specmode       = 'cube',
                    restfreq='{}GHz'.format(linei/1.0e9),
                    start          = '{}km/s'.format(vi),
                    width          = width_inp, 
                    nchan          = nchan,
                    weighting      = 'briggs',
                    robust         = robust, # 2= standard que maximiza S/N, 0.5 = balance standard para ganar resolucion y sacrificar un poco de S/N
                    niter          = 100000,
                    interactive    = True,
                    pbcor=True,
                    outframe='BARY')


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


def subtract_model(msfile, new_msfile, model_table, datacolumn='DATA', verbose=True):
    """
    function to subtract model from visibilities and save it as a new ms file

    Parameters
    ==========
    msfile: string, path to ms original file
    new_msfile: string, name of new ms file  
    uvtable_filename: string, path to visibility table with model
    datacolumn: string, column to extract (e.g. data or corrected)
    """

    print('subtracting model from visibilities..')

    # load model visibilities
    um, vm, Vrealm, Vimagm, wm, lamsm = np.require(np.loadtxt(model_table, unpack=True), requirements='C')
    Vmodel=Vrealm+Vimagm*1j
    
    # copy observations before modifying
    os.system('rm -r {}'.format(new_msfile))
    os.system('cp -r {} {}'.format(msfile, new_msfile))
    
    
    # open observations
    tb.open(new_msfile, nomodify=False)
    
    # get visibilities
    tb_columns = tb.colnames()
    #print(tb_columns)
    if datacolumn.upper() in tb_columns:
        data = tb.getcol(datacolumn.upper())
    else:
        raise KeyError("datacolumn {} is not available.".format(datacolumn))
    
    # reshape model visibilities
    shape=np.shape(data)
    nrows, ncol=shape[1], shape[2]
    Vmodel_reshaped=np.reshape(Vmodel, (nrows, ncol))
        
    vis_sub=data - Vmodel_reshaped 
    # save visibilities
    tb.putcol(datacolumn, vis_sub) # save modified data
    tb.close()
    print('done')

    
def simulate_model(msfile, new_msfile, model_table, datacolumn='DATA', verbose=True):
    """
    function to simulate model by replacing observed visibilities with model visibilities and save it as a new ms file

    Parameters
    ==========
    msfile: string, path to ms original file
    new_msfile: string, name of new ms file  
    uvtable_filename: string, path to visibility table with model
    datacolumn: string, column to extract (e.g. data or corrected)
    """

    print('replacing observed visibilities with model visibilities..')

    # load model visibilities
    um, vm, Vrealm, Vimagm, wm, lamsm = np.require(np.loadtxt(model_table, unpack=True), requirements='C')
    Vmodel=Vrealm+Vimagm*1j
    
    # copy observations before modifying
    os.system('rm -r {}'.format(new_msfile))
    os.system('cp -r {} {}'.format(msfile, new_msfile))
    
    
    # open observations
    tb.open(new_msfile, nomodify=False)
    
    # get visibilities
    tb_columns = tb.colnames()
    #print(tb_columns)
    if datacolumn.upper() in tb_columns:
        data = tb.getcol(datacolumn.upper())
    else:
        raise KeyError("datacolumn {} is not available.".format(datacolumn))
    
    # reshape model visibilities
    shape=np.shape(data)
    nrows, ncol=shape[1], shape[2]
    Vmodel_reshaped=np.reshape(Vmodel, (nrows, ncol))
            
    vis_model= data*0+Vmodel_reshaped 
    # save visibilities
    tb.putcol(datacolumn, vis_model) # save modified data
    tb.close()
    print('done')
    
    
def correct_vis(msfile, new_msfile, dRA=0., dDec=0., fsigma=1.0, flux_offset=1.0, datacolumn='DATA',):
    """
    function to apply phase center, amplitude and weight scale corrections to ms and save it as a new ms file. 

    Parameters
    ==========
    msfile: string, path to ms original file
    new_msfile: string, name of new ms file  
    dRa: float, source offset in arcsec from phase center. This function will apply the opposite shift to center observations.
    dDec: float, source offset in arcsec from phase center. This function will apply the opposite shift to center observations.
    fsigma: float, factor by which the uncertainties need to be multiplied to match dispersion. The weights will be divided by fsigma**2.
    flux_offset: float, amplitude offset to be applied to visibilities. The weights will be rescaled by 1/flux_offset**2.
    """
    print('applying corrections to visibilities..')

    # copy observations before modifying
    os.system('rm -r {}'.format(new_msfile))
    os.system('cp -r {} {}'.format(msfile, new_msfile))
    
    # open observations
    tb.open(new_msfile, nomodify=False)

    # get visibilities
    tb_columns = tb.colnames()

    if datacolumn.upper() in tb_columns:
        data = tb.getcol(datacolumn.upper())
    else:
        raise KeyError("datacolumn {} is not available.".format(datacolumn))
        
    spw = tb.getcol("DATA_DESC_ID")
    nspw = len(np.unique(spw))
    
    nrows, ncol = data.shape[1], data.shape[2]
    nchan=nrows
    
    # get weights
    weights_orig = tb.getcol("WEIGHT")
    
    # get uv points
    uvw = tb.getcol("UVW")
    u, v, w = [uvw[i, :] for i in range(3)]
    # repeat u and v nchan times
    u = np.tile(u, nchan).reshape((nrows,ncol))
    v = np.tile(v, nchan).reshape((nrows,ncol))

    ### convert u,v in m to lambda
    spw_path = tb.getkeyword('SPECTRAL_WINDOW').split()[-1]
    tb.close()
    tb.open(spw_path)
    freqs = tb.getcol('CHAN_FREQ')  # [Hz]
    tb.close()
    
    lams=np.zeros_like(data[0,:,:], dtype=float) # channel x spw
    for ichan in range(nchan):
        for irow in range(len(spw)):
            lams[ichan, irow] = clight_ms/(freqs[ichan, spw[irow]]) # m
    lams=lams.reshape((nrows,ncol))
    u=u/lams
    v=v/lams
    
    # apply phase center and amplitude corrections
    offra_rad=dRA*np.pi/180.0/3600.0 # RA offset in rad
    offdec_rad=dDec*np.pi/180.0/3600.0 # Dec offset in rad
    
    shift_correction=np.exp(-2.0*np.pi*1j*( u*offra_rad + v*offdec_rad ) ) # - sign is to be consistent with the definition of visibility = int I exp(+2pi ux+vy)dxdy. It has been tested and works well. This means that the source will be shifted by -offra and -offdec in the image space (well tested).
    
    data=data*shift_correction*flux_offset
    
    tb.open(new_msfile, nomodify=False)
    tb.putcol(datacolumn, data) # save modified data

    weights=weights_orig / flux_offset**2. / fsigma**2.
    tb.putcol("WEIGHT", weights)
    
    tb.close()
    print('done')


    
    
def apply_corrections(msfile, dRA=0., dDec=0., fsigma=1., flux_offset=1., model_vis='', extractvis=False, cfg='', clean=False, robust=0.5, dpix=0.1, Npix=512, uvtaper='' ):
    """
    function to subtract model and apply corrections to ms file and save it as a new ms file. 

    Parameters
    ==========
    msfile: string, path to ms original files
    dRa: floats, source offset in arcsec from phase center. This function will apply the opposite shift to center observations.
    dDec: float, source offset in arcsec from phase center. This function will apply the opposite shift to center observations.
    fsigma: float, factor by which the uncertainties need to be multiplied to match dispersion. The weights will be divided by fsigma**2.
    flux_offset: float, amplitude offset to be applied to visibilities. The weights will be rescaled by 1/flux_offset**2.
    model_vis: string, path to visibility table with model
    """
    
    
    tag=''
    print('correcting {}'.format(msfile))
    
    ### subtract model (background object or similar). The model should have the same offsets as the data before corrections.
    if model_vis!='':
        tag+='.SMGsub'
        new_msfile=msfile.split('/')[-1][:-3]+'.SMGsub.ms'
        subtract_model(msfile, new_msfile, model_vis)

        # redefine current ms file
        msfile=new_msfile

    ### correct visibilities if needed
    if dRA!=0. or dDec!=0. or fsigma!=1.0 or flux_offset!=1.0:
        tag+='.corrected'
        new_msfile=msfile.split('/')[-1][:-3]+'.corrected.ms'

        correct_vis(msfile, new_msfile, dRA=dRA, dDec=dDec, fsigma=fsigma, flux_offset=flux_offset)

    if extractvis:
        target=new_msfile.split('.')[0]
        uvtable_filename='{}.{}.continuum.fav.tav{}.txt'.format(target, cfg, tag) # assumes it is extracting fav.tav visibilities
        os.system('rm '+uvtable_filename)
        extract_visibilities(new_msfile, uvtable_filename)
    
    if clean:
        
        if uvtaper=='':
            image_name=new_msfile[:-3]+'.briggs.{:.1f}.{:d}.{:.3f}'.format(robust, Npix, dpix)
        else:
            image_name=new_msfile[:-3]+'.briggs.{:.1f}.{}.{:d}.{:.3f}'.format(robust, uvtaper, Npix, dpix)
        
        os.system('rm -r '+image_name+'*')

        tclean(
            vis            = new_msfile,
            imagename      = image_name,
            cell           = '{}arcsec'.format(dpix),
            imsize         = [Npix,Npix],
            specmode           = 'mfs',
            # gridder=gridder,
            weighting      = 'briggs',
            robust         = robust,
            uvtaper=uvtaper,
            niter          = 10000,
            interactive    = True,
            pbcor=False)
        
        imview(image_name+'.image')

#             # export fits file
#             fits_image=image_name+'.fits'
#             fits_pbcor=image_name+'.pbcor.fits'
#             fits_pb=image_name+'.pb.fits'
#             fits_model=image_name+'.model.fits'

#             exportfits(
#                 imagename = image_name+'.image',
#                 fitsimage = fits_image,
#                 overwrite = True
#             )
#             exportfits(
#                 imagename = image_name+'.image.pbcor',
#                 fitsimage = fits_pbcor,
#                 overwrite = True
#             )
#             exportfits(
#                 imagename = image_name+'.pb',
#                 fitsimage = fits_pb,
#                 overwrite = True
#             )
#             exportfits(
#                 imagename = image_name+'.model',
#                 fitsimage = fits_model,
#                 overwrite = True
#             )
            
            
    
def concatenate(msfiles, concatvis, refdir='', remove_files=False, reframe='J2000'): 
    
    for i, msi in enumerate(msfiles):
        ## fixplanet to change coordinates of phase center (no phase shift applied).
        fields=get_fields(msi, field_name='')
        fixplanets(vis=msi, direction=refdir, field=fields)#, field=field)
        if reframe=='ICRS':
            change_to_ICRS(msi)
    
#         listobs(msi, listfile=msi.replace('.ms', '.list'), overwrite=True)
        
    os.system('rm -r '+concatvis)
    concat(vis=msfiles, concatvis=concatvis, copypointing=False)
    if remove_files:
        for i, msi in enumerate(msfiles):
            os.system('rm -r '+msi )

def change_to_ICRS(msfile):

    tb.open(msfile+'/FIELD', nomodify=False)
    for colname in ['PhaseDir_Ref', 'DelayDir_Ref', 'RefDir_Ref']:
        a = tb.getcol(colname)
        for i in range(len(a)):
            if a[i]==0: # J2000
                a[i]=21 # ICRS
                print('Found J2000 '+colname+' for field '+str(i)+' and changed it to ICRS.')
        tb.putcol(colname, a)
    tb.close()

    tb.open(msfile+'/SOURCE', nomodify=False)
    x = tb.getcolkeywords('DIRECTION')   
    if x['MEASINFO']['Ref'] == 'J2000':
        x['MEASINFO']['Ref'] = 'ICRS'
        tb.putcolkeywords('DIRECTION', x)
        print('Found J2000 DIRECTION in SOURCE table and changed it to ICRS.')

    tb.close()

def change_to_ICRS_image(imagename):
    
    ia.open(imagename)
    mycs = ia.coordsys().torecord()
    if mycs['direction0']['conversionSystem'] == 'J2000':
        mycs['direction0']['conversionSystem'] = 'ICRS'
        print("Found J2000 conversion system and changed it to ICRS")
    if mycs['direction0']['system'] == 'J2000':
        mycs['direction0']['system'] = 'ICRS'
        print("Found J2000 direction system and changed it to ICRS")
    ia.setcoordsys(mycs)
    ia.close()


###########################
####### JvM correction ####
####### provided by Nicolas Kurtovic based on MAPS codes
###########################

def gaussian_eval(params, data, center):
    '''
    Returns a gaussian with the given parameters
    '''
    width_x, width_y, rotation = params
    rotation = 90-rotation
    width_x = float(width_x)
    width_y = float(width_y)

    rotation = np.deg2rad(rotation)
    x, y = np.indices(data.shape) - center

    xp = x * np.cos(rotation) - y * np.sin(rotation)
    yp = x * np.sin(rotation) + y * np.cos(rotation)
    g = 1. * np.exp( -((xp / width_x)**2 + \
                       (yp / width_y)**2) / 2.)
    return g


def gaussian2D(params, nrow):
    '''
    Returns a gaussian with the given parameters
    '''
    width_x, width_y, rotation = params
    rotation = 90-rotation

    rotation = np.deg2rad(rotation)
    x, y = np.indices((int(nrow*2+1),int(nrow*2+1))) - nrow

    xp = x * np.cos(rotation) - y * np.sin(rotation)
    yp = x * np.sin(rotation) + y * np.cos(rotation)
    g = 1. * np.exp( -((xp / width_x)**2 + \
                       (yp / width_y)**2) / 2.)
    return g


def beam_chi2(params, psf, nrow):
    '''
    Calculates the difference between a Gaussian and the PSF
    '''
    psf_ravel = psf[~np.isnan(psf)]
    gaussian = gaussian2D(params, nrow)[~np.isnan(psf)]
    chi2 = np.sum( (gaussian - psf_ravel)**2 )
    # Return
    return chi2


def do_JvM_correction_and_get_epsilon(root, taper_match=None):
    '''
    Applies the JvM correction and returns the epsilon value. This function is an adaptation of a function by Nicolas Jurtovic.
    '''
    # Get the psf file to fit
    psf_file = root + '.psf'
    model_file = root + '.model'
    residual_file = root + '.residual'
    npix_window = 301  # Bigger than the 201 default from MAPS

    # Open psf and read off the metadata
    ia.open(psf_file)
    psf_data_raw = ia.getregion()
    hdr = ia.summary(list=False)
    ia.close()

    delta = np.abs(hdr['incr'][0]*206265)
    try:
        rb = hdr['restoringbeam']
        major = rb['major']['value']
        minor = rb['minor']['value']
        phi = rb['positionangle']['value']
    except:
        major = hdr['perplanebeams']['beams']['*0']['*0']['major']['value']
        minor = hdr['perplanebeams']['beams']['*0']['*0']['minor']['value']
        phi = hdr['perplanebeams']['beams']['*0']['*0']['positionangle']['value']

    print('The CASA fitted beam is ' + str(major) + 'x' + str(minor) + ' at ' + str(phi) + 'deg')

    npix = psf_data_raw.shape[0]         # Assume image is square

    # Check if image cube, or just single psf; this example doesn't handle the full
    # polarization case - implicitly assumes we can drop Stokes
    # If single psf, add an axis so we can use a single loop
    psf_data = np.squeeze(psf_data_raw)
    if len(psf_data.shape) == 2:
        psf_data = np.expand_dims(psf_data, axis=2)

    # Roll the axes to make looping more straightforward
    psf_rolled = np.rollaxis(psf_data,2)

    # Window out the region we want to consider
    i_min = int(npix/2-(npix_window-1)/2)
    i_max = int(npix/2+(npix_window-1)/2 + 1)
#    for j in range(10): I don't know why it used to be 10
    for j in range(len(psf_rolled)):
# I don't know whythis used to be 12
#        psf_windowed = psf_rolled[12*j][i_min:i_max,i_min:i_max]
        psf_windowed = psf_rolled[j][i_min:i_max,i_min:i_max]
        if np.sum(psf_windowed) > 0:
            break

    # Mask out anything beyond the first null
    psf_windowed[psf_windowed<0.] = -1.
    psf_windowed = np.fft.fftshift(psf_windowed)

    for i in range(psf_windowed.shape[0]):
        left_edge = np.argmax(psf_windowed[i] < 0.)
        right_edge = npix_window-np.argmax(psf_windowed[i][::-1] < 0.)
        psf_windowed[i][left_edge:right_edge] = 0.

    psf_windowed = np.fft.fftshift(psf_windowed)

    # Create a clean beam to evaluate against
    clean_beam = gaussian_eval([major/2.355/delta, minor/2.355/delta, phi], psf_windowed, 
                               (npix_window-1)/2)

    # Calculate epsilon
    epsilon = np.sum(clean_beam)/np.sum(psf_windowed)
    print('Epsilon = ' + str(epsilon))

    # Check if taper match is requested
    if taper_match:
        # create the convolved model
        convolved_temp_image = '{:s}_convolved_model_temp.image'.format(root)
        imsmooth(imagename=model_file, \
                 major=str(major)+'arcsec', \
                 minor=str(minor)+'arcsec', \
                 pa=str(phi)+'deg', \
                 targetres=True, \
                 outfile=convolved_temp_image)

        # doing the correction
        try:
            shutil.rmtree(root+'.JvMcorr.temp.image')
        except:
            pass
        immath(imagename=[convolved_temp_image, residual_file], \
               expr='IM0 + ' + str(epsilon) + '*IM1', \
               outfile=root+'.JvMcorr.temp.image')

        # now smooth to taper_match
        try:
            shutil.rmtree(root+'.JvMcorr.image')
        except:
            pass

        imsmooth(imagename=root+'.JvMcorr.temp.image', \
                 major=str(taper_match)+'arcsec', \
                 minor=str(taper_match)+'arcsec', \
                 pa=str(phi)+'deg', \
                 targetres=True, \
                 outfile=root+'.JvMcorr.image')

        if not os.path.exists(root+'.JvMcorr.image'):
            imsmooth(imagename=root+'.JvMcorr.temp.image', \
                     major=str(taper_match*1.05)+'arcsec', \
                     minor=str(taper_match*1.05)+'arcsec', \
                     pa=str(phi)+'deg', \
                     targetres=True, \
                     outfile=root+'.JvMcorr.image')
            if not os.path.exists(root+'.JvMcorr.image'):
                imsmooth(imagename=root+'.JvMcorr.temp.image', \
                         major=str(taper_match*1.1)+'arcsec', \
                         minor=str(taper_match*1.1)+'arcsec', \
                         pa=str(phi)+'deg', \
                         targetres=True, \
                         outfile=root+'.JvMcorr.image')
        print('Wrote ' + root + '.JvMcorr.image')

        # clean up
        shutil.rmtree(convolved_temp_image)
        shutil.rmtree(root+'.JvMcorr.temp.image')

    else:
        # Regular JvM correction
        # create the convolved model
        convolved_temp_image = '{:s}_convolved_model_temp.image'.format(root)
        imsmooth(imagename=model_file, \
                 major=str(major)+'arcsec', minor=str(minor)+'arcsec',
                 pa=str(phi)+'deg', \
                 targetres=True, \
                 outfile=convolved_temp_image)
        # doing the correction
        try:
            shutil.rmtree(root+'.JvMcorr.image')
        except:
            pass
        immath(imagename=[convolved_temp_image, residual_file], \
               expr='IM0 + ' + str(epsilon) + '*IM1', \
               outfile=root+'.JvMcorr.image')
        print('Wrote ' + root + '.JvMcorr.image')

        # clean up
        shutil.rmtree(convolved_temp_image)
        # 'lowres' JvM correction
        # create the convolved model
        imsmooth(imagename=model_file, \
                 major=str(np.sqrt(1./epsilon)*major)+'arcsec', \
                 minor=str(np.sqrt(1./epsilon)*minor)+'arcsec', \
                 pa=str(phi)+'deg', \
                 targetres=True, \
                 outfile=convolved_temp_image)
        # doing the correction
        try:
            shutil.rmtree(root+'.JvMcorr_lowres.image')
        except:
            pass
        immath(imagename=[convolved_temp_image, residual_file], \
               expr='IM0 + IM1',
               outfile=root+'.JvMcorr_lowres.image')
        print('Wrote ' + root + '.JvMcorr_lowres.image')
        # clean up
        shutil.rmtree(convolved_temp_image)
    # Return
    return epsilon








