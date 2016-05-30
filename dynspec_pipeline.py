# run this pipeline from the directory for the particular observation you want to process
# (e.g., /data/jrv/15A-416/YZCMi/1)
# The only files in this directory should be directories named for the bands observed (e.g., P/, L/, S/).
#
# This pipeline should be run after the VLA pipeline.  It will look for ms files in locations such
#   as L/YZCMi_1L.ms, which should have the calibrated data column from the VLA pipeline.

overwrite = False # set this to True if you want existing .src.ms files, tar files, etc produced by pipeline
# to be overwritten (otherwise pipeline skips steps for which files have already been produced)

#bands = os.listdir(mydir)
bands = ['L','S'] # for testing

import os
from glob import glob
from dynspec.tbavg import tbavg,dyn_spec
from dynspec.extract_dynspec import saveTxt

mydir = os.getcwd()
tmp = mydir[1:].split('/')
proj = tmp[2]
src = tmp[3]
obs = tmp[4]

def get_phasecen(obs_str):
    # load phasecen dict from file phasecen.txt (in casa_utils) and find entry for this observation
    # format of one line of phasecen dict: 15A-416_YZCMi_1---J2000 7h44m39.810 3d33m1.91
    
    fname = '/data/jrv/casa_utils/dynspec/phasecen.txt'
    f = open(fname)
    lines = f.readlines()
    f.close()
    pcen_dict = {}
    for l in lines:
        [obs_name,coords]=l.rstrip().split('---')
        pcen_dict[obs_name] = coords
    return pcen_dict[obs_str]

def get_clean_spws(band):
    # relatively RFI-free spw's that will be used for imaging and subtracting bg srcs
    clean_spw_dict = {'P':'','L':'6,7,11~13','S':'16,19~26','C':''}
    return clean_spw_dict[band]

def get_clean_scans(band):
    # read clean (flare-free) scans from file clean_scans.txt in working directory
    # format of a line in this file should be:
    # L   0~13,15~17
    fname = mydir + '/clean_scans.txt'
    try:
        f = open(fname)
        lines = f.readlines()
        f.close()
    except:
        print 'no file exists at',fname,'- assuming all scans are flare-free'
        return ''
    clean_scan_dict = {}
    for l in lines:
        [band_name,scans]=l.rstrip().split()
        clean_scan_dict[band_name] = scans
    try:
        return clean_scan_dict[band]
    except KeyError:
        print 'no line for band', band, 'in file', fname, '- assuming all scans are flare-free'

def im_params(vis):
    # return a good cell size and imsize (in pixels) for the ms
    im.open(vis)
    returnval,pixels,cell,facets,phasecen=im.advise(takeadvice=False)
    im.done()
    if facets > 1:
        print 'Warning: recommended # of facets > 1.  Imaging will proceed with only one facet.'
    pixelsize = str(cell['value'] * 2/3) + cell['unit'] # so we will have at least 3 pixels per synthesized beam (perhaps more)
    npixels = int(pixels * 3/2 * 4) # imsize will be 4 * primary beam at center freq - useful at low freq b/c we get bright srcs outside beam
    return pixelsize,npixels

# read phase center from file /data/jrv/casa_utils/dynspec/phasecen.txt
obs_str = proj + '_' + src + '_' + obs
phasecen = get_phasecen(obs_str)

# create plots directory if it does not exist yet
plotdir = mydir + '/plots/'
if not os.path.exists(plotdir):
    os.system('mkdir '+plotdir)

for band in bands:
    os.chdir(band)
    msroot = src + '_' + obs + band
    msfile = msroot + '.ms'
    srcms = msroot + '.src.ms'
    srctbavg = msroot + '.tbavg.ms'
    smallms = msroot + '.small.ms'
    clean_spws = get_clean_spws(band)
    smalltbavg = msroot + '.small_tbavg.ms'
    smallim = msroot + '.smallim'

    print '\n----------------'+msroot+'--------------------'

    # split src-only ms from calibrated data
    if not os.path.exists(srcms) or overwrite:  # skip this step if srcms already exists
        if os.path.exists(srcms):
            rmtables(srcms)
        print 'creating src-only ms', srcms, 'from calibrated data column of', msfile
        split(vis=msfile,outputvis=srcms,field=src,keepflags=False)
        print 'shifting phase center of', srcms, 'to', phasecen
        fixvis(vis=srcms,outputvis=srcms,phasecenter=phasecen)
    else:
        print srcms, 'already exists, now deleting MODEL_DATA and CORRECTED_DATA and assuming phase center already shifted'
        # erase corrected data column to remove source subtraction from previous pipeline runs
        clearcal(srcms)

    # use backup as tar file
    tarfile = srcms + '.tar'
    if not os.path.exists(tarfile) or overwrite: # skip this step if tarfile already exists
        print 'saving',tarfile,'with fixvis applied'
        os.system('tar -cf ' + tarfile + ' ' + srcms)
    else:
        print tarfile, 'already exists, skipping the creation of this file'
    # to retrieve pointing center (in radians - final index 0 for RA, 1 for dec):
    #     ptc = vishead(vis=srcms,mode='get',hdkey='ptcs')[0]['r1'][0]

    # make small ms - clean spw's only, avg'd in time (30s) and freq (64 channels/whole spw)
    if not os.path.exists(smallms) or overwrite:
        if os.path.exists(smallms):
            rmtables(smallms)
        print 'avging clean spws ('+clean_spws+') in', srcms,'over 30s and 64 channels to create', smallms
        split(vis=srcms,outputvis=smallms,spw=clean_spws,datacolumn='data',width=64,timebin='30s',keepflags=False)
    else:
        print smallms, 'already exists, skipping the creation of this file but erasing MODEL_DATA and CORRECTED_DATA columns'
        # erase corrected data column to remove source subtraction from previous pipeline runs
        clearcal(smallms)

    # plot Re(Vis) and Im(Vis) time series and note any times that have obvious flares
    # save plots to a directory then ask user to review them
    # work from smallms b/c for some reason the avgspw=True option doesn't work after tbavg
    plotfile = plotdir + msroot + '_dirty_tseries.png'
    plotms(vis=smallms,xaxis='scan',yaxis='real',avgspw=True,avgbaseline=True,correlation='RR,LL',coloraxis='corr',plotfile=plotfile,overwrite=True,showgui=False)
    raw_input('Check plots/'+msroot+'_dirty_tseries.png to identify flare-free scans, add them to clean_scans.txt, then hit Return.')
    clean_scans = get_clean_scans(band)
    print 'Scans used as flare-free for',band,'band (blank is all):', clean_scans
    
    ## CREATE MODEL OF BG SOURCES ##
    
    # get parameters for imaging
    pixel_size,imsize = im_params(smallms)    

    # measure RMS of Stokes Q&U dirty image during non-flaring scans
    #  (when subtracting BG srcs later we assume unpolarized, so this will prevent overcleaning polarized BG src)
    dirtyim='dirty'
    if len(glob(dirtyim+'.*'))>0:
        rmtables(dirtyim + '.*') # have to delete old image files (esp. model) so clean will start from scratch
    print 'Creating dirty image during non-flaring scans - will use Stokes Q&U RMS to set threshold (3*RMS) for cleaning'
    clean(vis=smallms,imagename=dirtyim,scan=clean_scans,imsize=imsize,cell=pixel_size,niter=0,stokes='IQUV') # create dirty image
    rmsQ = imstat(imagename=dirtyim+'.image',stokes='Q')['sigma']
    rmsU = imstat(imagename=dirtyim+'.image',stokes='U')['sigma']
    rmsQU = max(rmsQ,rmsU)[0]*1000
    print 'RMS of Stokes Q or U dirty image (max of the two):',rmsQU, 'mJy'
    threshold = str(rmsQU * 3) + 'mJy'

    # clean smallms using only the non-flaring scans with nterms=2 since we have large fractional bandwidth (may need nterms=3?)
    if len(glob(smallim+'.*'))>0:
        rmtables(smallim + '.*') # have to delete old image files (esp. model) so clean will start from scratch
    print 'Imaging non-flaring scans in',smallms,'with nterms=2 - image file', smallim+'.image.tt0'
    clean(vis=smallms,imagename=smallim,scan=clean_scans,nterms=2,imsize=imsize,cell=pixel_size,niter=10000,threshold=threshold,cyclefactor=5)
    # high cyclefactor is critical to successful clean - this implies "bad PSF" --> poorly known beam map? (due to few antennas? wide bandwidth?)
    # default gain of 0.1 is fine - I tested 0.01 and 0.3 as well - not a huge difference in performance (0.1 had lowest RMS by small amount)

    # mask center of clean model - use immath to subtract center of image from whole image?
    model0 = smallim+'.model.tt0'
    model1 = smallim+'.model.tt1'
    cen_mask = 'center_mask'
    bgmodel0 = smallim + '.bgmodel0'
    bgmodel1 = smallim + '.bgmodel1'
    ncen = str(int(imsize/2)) # center pixel
    circle_cen = 'circle[['+ncen+'pix,'+ncen+'pix],5pix]' # a 5-pixel diameter circle in the image center
    print 'Excluding image center region', circle_cen, 'from model of bg srcs'
    makemask(mode='copy',inpimage=model0,output=cen_mask,inpmask=circle_cen,overwrite=True) # create mask with ones in center
    immath(imagename=[model0,cen_mask],outfile=bgmodel0,expr='IM0*(1-IM1)')
    immath(imagename=[model1,cen_mask],outfile=bgmodel1,expr='IM0*(1-IM1)')
    # bgmodel0,1 are the same as model0,1 but with the center pixels removed (radius 5 pixels)
    

    ## SMALLMS - REMOVE BG SRCS AND PLOT TSERIES, SAVE DYNSPEC ##
    
    # ft to populate model column with clean model - smallms
    ft(vis=smallms,nterms=2,model=[bgmodel0,bgmodel1],usescratch=True) # inspected Re(vis) for model vs data in plotms - looks good!
    
    # uvsub to subtract model column from data column and put it in corrected column - smallms
    uvsub(vis=smallms)
    
    # create dirty image of smallms w/o BG srcs and save png to plots directory for user to inspect
    dirty_nobg=msroot+'.dirty_nobg'
    imfile = plotdir+dirty_nobg+'.png'
    if len(glob(dirty_nobg+'.*'))>0:
        rmtables(dirty_nobg+'.*') # delete old image files
    print 'Creating dirty image during non-flaring scans with bg srcs subtracted:',imfile
    clean(vis=smallms,imagename=dirty_nobg,scan=clean_scans,imsize=imsize,cell=pixel_size,niter=0,stokes='IQUV')
    raster={'file':dirty_nobg+'.image','colormap':'Smooth 2'}
    imview(raster=raster,out=imfile)
    
    # plot tseries for smallms
    plotfile = plotdir + msroot + '_tseries.png'
    imfile = plotdir + msroot + '_im_tseries.png' # make sure there is no variation in Im(vis)
    print 'Plotting time series for', band,'band with bg srcs removed:',plotfile, '(+im in filename for Im(vis))'
    plotms(vis=smallms,xaxis='time',yaxis='real',ydatacolumn='corrected',avgspw=True,avgbaseline=True,correlation='RR,LL',coloraxis='corr',plotfile=plotfile,overwrite=True,showgui=False)
    plotms(vis=smallms,xaxis='time',yaxis='imag',ydatacolumn='corrected',avgspw=True,avgbaseline=True,correlation='RR,LL',coloraxis='corr',plotfile=imfile,overwrite=True,showgui=False)

    # avg small ms over all baselines
    if os.path.exists(smalltbavg):
        rmtables(smalltbavg)
    print 'running tbavg on', smallms, '(bg subtracted) to create', smalltbavg
    try:
        tbavg(split,smallms,smalltbavg,speed='fast',weight_mode='flat',datacolumn='corrected')
    except:
        print 'tbavg failed b/c table',smalltbavg, 'is already open in the CASA cache - restart CASA to fix this problem'
    
    # create dynspec for smallms
    spec = dyn_spec(smalltbavg)
    if os.path.exists(smalltbavg+'.dynspec'):
        os.system('rm -rf ' + smalltbavg + '.dynspec')
    saveTxt(spec,smalltbavg)
    
    # looking at this dynamic spectrum, it looks like it could use some RFI excision...
    # can I automate RFI flagging?  in previous experience, flagging by hand did not
    #  greatly improve the dynamic spectrum compared to the pipeline, I think
    # Let's go ahead for now and can come back to this if needed

    ## FULL SRC MS - SUBTRACT BG SRCS AND PLOT TSERIES, SAVE DYNSPEC ##

    # ft to populate model column with clean model - srcms
    print 'subtracting bg src model from', srcms, '(ft + uvsub)'
    ft(vis=srcms,nterms=2,model=[bgmodel0,bgmodel1],usescratch=True)
    
    # uvsub to subtract model column from data column and put it in corrected column - srcms
    uvsub(vis=srcms)
    
    # avg src ms over all baselines
    if os.path.exists(srctbavg):
        rmtables(srctbavg)
    print 'running tbavg on', srcms, '(bg subtracted) to create', srctbavg
    try:
        tbavg(split,srcms,srctbavg,speed='fast',weight_mode='flat',datacolumn='corrected')
    except:
        print 'tbavg failed b/c table',srctbavg, 'is already open in the CASA cache - restart CASA to fix this problem'
    
    # create dynspec for srcms
    spec = dyn_spec(srctbavg)
    if os.path.exists(srctbavg+'.dynspec'):
        os.system('rm -rf ' + srctbavg + '.dynspec')
    saveTxt(spec,srctbavg)

    os.chdir('..')
