# run this pipeline from the directory for the particular observation you want to process
# (e.g., /data/jrv/15A-416/YZCMi/1/L)
#
# This pipeline should be run after the VLA pipeline.  It will look for ms files in locations such
#   as YZCMi_1L.ms (where the name can be inferred from the directory location), which should have
#   the calibrated data column from the VLA pipeline.
#
# If you want to change default parameters for pipeline, create dictionary pipe_params before running this
#   script with execfile.  Default for pipe_params is: {'overwrite':False,'makebig':False,'interactive':False}

# Steps:
# 1) Definitions - parameters, field names, filenames, phasecen
# 2) ms file creation
#    a) Use mstransform to make srcms (calibrated science target only w/ phase center on star, single spw)
#    b) Run automated flagging (rflag) 2x on single-spw srcms
#    c) Split smallms (time- and freq-avg'd srcms)
# 3) Create BG src model from smallms
#    a) Plot raw time series (will show sidelobes of BG srcs) and read non-flaring scans/spws from file
#    b) Image and clean smallms in non-flare scans/spws w/ nterms>1 (2? 3?)
#    c) Mask central star out of model image to create BG model image
# 4) Subtract BG srcs and avg over baselines - run on smallms, and srcms if makebig=True
#    a) FT BG model image then uvsub
#    b) Average over all baselines
#    c) Extract dynamic spectrum from baseline-avg'd ms
#    d) Diagnostic plots of BG-subtracted smallms - time series & dirty image

import os
from glob import glob
from dynspec.tbavg import tbavg,dyn_spec
from dynspec.extract_dynspec import saveTxt
from dynspec.pipeline_utils import *



### 1) DEFINITIONS ###
# run parameters, field names, filenames, phasecen

try:
    pipe_params
except:
    pipe_params = {}
overwrite = get_params('overwrite',pipe_params)
makebig = get_params('makebig',pipe_params)
interactive = get_params('interactive',pipe_params)

names = get_names()
sb = names['sb']
band = names['band']
print '\n----------------'+sb+'--------------------'

plotdir = get_plot_dir() # file names
msfile = sb + '.ms'
rawsrc = sb + '.src.ms'
srcms = sb + '.star.ms'
srctbavg = sb + '.tbavg.ms'
smallms = sb + '.star.small.ms'
smalltbavg = sb + '.star.small_tbavg.ms'
smallim = sb + '.smallim'
dirtyim='dirty'
bgmodel0 = smallim + '.bgmodel0'
bgmodel1 = smallim + '.bgmodel1'

if os.path.exists(srcms):
    fields = get_fields(srcms,vishead) # field names
elif os.path.exists(msfile):
    fields = get_fields(msfile,vishead)
else:
    fields = get_fields(rawsrc,vishead)
srcname = fields.get('src')
bpcal = fields.get('bpcal')
gcal = fields.get('gcal')
phasecen=get_phasecen()

clean_scans,clean_spws = get_clean_scans()
threshold='0.2mJy'
imsize_d=1024

pixel_size,imsize = im_params(srcms)
if imsize > 8192:
    print 'Warning: requested imsize of',imsize,'pix is large, instead setting to imsize=8192 (may cut off full FOV)'
    imsize = 8192


### 2) MS FILE CREATION ###

## a) Use mstransform to make srcms (calibrated science target only w/ phase center on star, single spw) ##

# use mstransform to create 1-spw src-only ms from calibrated data
# in some directories, we have pipeline output (msfile) and in
#    others, we only have target data from pipeline output (rawsrc)
if band=='S':
    spw='16~31'
else:
    spw=''
if not os.path.exists(srcms):  # skip this step if srcms already exists
    if os.path.exists(srcms):
        rmtables(srcms)
    if os.path.exists(srcms+'.flagversions'):
        os.system('rm -rf '+srcms+'.flagversions')
    if os.path.exists(msfile): # msfile exists so use it to create srcms
        print 'creating 1-spw src-only ms', srcms, 'from calibrated data column of', msfile
        mstransform(vis=msfile,outputvis=srcms,field=srcname,combinespws=True,keepflags=False,spw=spw)
    else: # msfile does not exist so use rawsrc to create srcms
        print 'creating 1-spw src-only ms', srcms, 'from data column of', rawsrc, '(which is calibrated)'
        mstransform(vis=rawsrc,outputvis=srcms,combinespws=True,keepflags=False,datacolumn='data',spw=spw)
    print 'shifting phase center of', srcms, 'to', phasecen
    fixvis(vis=srcms,outputvis=srcms,phasecenter=phasecen)
    ## b) Run automated flagging (rflag) 2x on single-spw srcms ##
    #         This is so RFI doesn't impact BG model and also for dynspec quality
    summary_1 = field_flag_summary(srcms,flagdata)
    flagdata(vis=srcms,mode='rflag',freqdevscale=3.0) # freqdevscale=3.0: flag channels w/ local (3-chan) rms 3x greater than global rms
    summary_2 = field_flag_summary(srcms,flagdata)
    flagdata(vis=srcms,mode='rflag',freqdevscale=3.0)
    summary_3 = field_flag_summary(srcms,flagdata)
elif makebig:
    clearcal(srcms) # erase corrected data column to remove source subtraction from previous pipeline runs
    #fixvis(vis=srcms,outputvis=srcms,phasecenter=phasecen)
    if overwrite:
        print srcms, 'already exists, now deleting MODEL_DATA and CORRECTED_DATA and shifting phase center to', phasecen
        fixvis(vis=srcms,outputvis=srcms,phasecenter=phasecen)
    else:
        print srcms, 'already exists, now deleting MODEL_DATA and CORRECTED_DATA and assuming phase center already shifted'

# backup src-only ms as tar file (in case calibrated pipeline ms and srcms get corrupted)
tarfile = srcms + '.tar'
if not os.path.exists(tarfile): # skip this step if tarfile already exists
    print 'saving',tarfile,'with fixvis applied'
    os.system('tar -cf ' + tarfile + ' ' + srcms)
else:
    print tarfile, 'already exists, skipping the creation of this file - warning: phase center in tarfile may not be most recent and flagging may not have been applied'

## c) Split smallms - srcms avg'd in time (30s) and freq (16 channels) ##

if not os.path.exists(smallms) or overwrite:
    if os.path.exists(smallms):
        rmtables(smallms)
    print 'avging', srcms,'over 30s and 16 channels to create', smallms
    split(vis=srcms,outputvis=smallms,datacolumn='data',width=16,timebin='30s',keepflags=False)
else:
    print smallms, 'already exists, skipping the creation of this file but erasing MODEL_DATA and CORRECTED_DATA columns'
    # erase corrected data column to remove source subtraction from previous pipeline runs
    clearcal(smallms)

'''

### 3) CREATE BACKGROUND SOURCE MODEL FROM SMALLMS ###

## a) Plot raw time series (will show sidelobes of BG srcs) and read non-flaring scans/spws from file ##

plotfile = plotdir + sb + '_dirty_tseries.png'
im_plotfile = plotdir + sb + '_dirty_tseries_im.png'
plotms(vis=smallms,xaxis='scan',yaxis='real',avgchannel='64',avgbaseline=True,correlation='RR,LL',coloraxis='corr',plotfile=plotfile,overwrite=True,showgui=False)
plotms(vis=smallms,xaxis='scan',yaxis='imag',avgchannel='64',avgbaseline=True,correlation='RR,LL',coloraxis='corr',plotfile=im_plotfile,overwrite=True,showgui=False)
if interactive:
    raw_input('Check plots/'+sb+'_dirty_tseries[_im].png to identify flare-free scans, add them to dynspec/clean_scans.txt, then hit Return.')
clean_scans,clean_spws = get_clean_scans()
print 'Scans and spws used as flare-free for',sb,'(blank is all): scan=\''+clean_scans+'\', spw=\''+clean_spws+'\''


## b) Image and clean smallms in non-flare scans/spws w/ nterms=2 ##
# get parameters for imaging
pixel_size,imsize = im_params(srcms)
if imsize > 8192:
    print 'Warning: requested imsize of',imsize,'pix is large, instead setting to imsize=8192 (may cut off full FOV)'
    imsize = 8192

# measure RMS of Stokes Q&U dirty image during non-flaring scans
#  (when subtracting BG srcs later we assume unpolarized, so this will prevent overcleaning polarized BG src)
if len(glob(dirtyim+'.*'))>0:
    os.system('rm -rf '+dirtyim + '.*') # have to delete old image files (esp. model) so clean will start from scratch
print 'Creating dirty image during non-flaring scans - will use Stokes Q&U RMS to set threshold (3*RMS) for cleaning'
imsize_d = min(imsize,1024)
clean(vis=smallms,imagename=dirtyim,scan=clean_scans,spw=clean_spws,imsize=imsize_d,cell=pixel_size,niter=0,stokes='IQUV') # create dirty image

rmsQ = imstat(imagename=dirtyim+'.image',stokes='Q')['sigma']
rmsU = imstat(imagename=dirtyim+'.image',stokes='U')['sigma']
rmsV = imstat(imagename=dirtyim+'.image',stokes='V')['sigma']
rmsQUV = min([rmsQ,rmsU,rmsV])[0]*1000
print 'RMS of Stokes Q,U,V dirty image (min of the three):',rmsQUV, 'mJy'
threshold = str(rmsQUV * 3) + 'mJy'

clearcal(smallms)
# clean smallms using only the non-flaring scans with nterms=2 since we have large fractional bandwidth
#  - may need nterms=3 to help remove bright srcs from ADLeo, EVLac, but nterms=3 gave worse residuals than nterms=2 for YZCMi_1S
if len(glob(smallim+'.*'))>0:
    rmtables(smallim + '.*') # have to delete old image files (esp. model) so clean will start from scratch
print 'Imaging non-flaring scans in',smallms,'with nterms=2 - image file', smallim+'.image.tt0'
clean(vis=smallms,imagename=smallim,scan=clean_scans,spw=clean_spws,nterms=2,imsize=imsize,cell=pixel_size,niter=3000,threshold=threshold,cyclefactor=5,interactive=interactive)

# high cyclefactor is critical to successful clean - this implies "bad PSF" --> poorly known beam map? (due to few antennas? wide bandwidth?)
# default gain of 0.1 is fine - I tested 0.01 and 0.3 as well - not a huge difference in performance (0.1 had lowest RMS by small amount)
# YZCMi_1S: selfcal made RMS/residuals worse - not bright enough? overcleaned before selfcal? for now, don't do selfcal.

#os.system('rm -rf aproj.*')
#clean(vis=smallms,imagename='aproj',scan=clean_scans,spw=clean_spws,nterms=2,imsize=imsize,cell=pixel_size,niter=3000,threshold=threshold,cyclefactor=5,interactive=interactive,gridmode='aprojection',wprojplanes=-1)
# uses small circle for FOV - this is a problem...
# next step: try tclean w/ gridder='awproj' to see if this works better


## c) Mask central star out of model image to create BG model image ##

# mask center of clean model
model0 = smallim+'.model.tt0'
model1 = smallim+'.model.tt1'
cen_mask = 'center_mask'
if os.path.exists(cen_mask):
    rmtables(cen_mask)
bgmodel0 = smallim + '.bgmodel0'
bgmodel1 = smallim + '.bgmodel1'
ncen = str(int(imsize/2)) # center pixel
circle_cen = 'circle[['+ncen+'pix,'+ncen+'pix],5pix]' # a 5-pixel diameter circle in the image center
print 'Excluding image center region', circle_cen, 'from model of bg srcs'
makemask(mode='copy',inpimage=model0,output=cen_mask,inpmask=circle_cen,overwrite=True) # create mask with ones in center
immath(imagename=[model0,cen_mask],outfile=bgmodel0,expr='IM0*(1-IM1)')
immath(imagename=[model1,cen_mask],outfile=bgmodel1,expr='IM0*(1-IM1)')
# bgmodel0,1 are the same as model0,1 but with the center pixels removed (radius 5 pixels)

'''
### 4) SUBTRACT BG AND CREATE DYNSPEC ###

mslist = [smallms]
tblist = [smalltbavg]
if makebig:
    mslist.append(srcms)
    tblist.append(srctbavg)

for (ms,tb) in zip(mslist,tblist):
    
    ## a) FT background model image then uvsub ##
    
    # ft to populate MODEL_DATA column with clean model
    print 'subtracting bg src model from', ms, '(ft + uvsub)'
    ft(vis=ms,nterms=2,model=[bgmodel0,bgmodel1],usescratch=True)
    
    # uvsub to subtract model column from data column and put it in corrected column
    uvsub(vis=ms)
    
    ## b) Avg over all baselines ##

    # tbavg creates new "single-baseline" ms
    if os.path.exists(tb):
        os.system('rm -rf '+tb)
    print 'running tbavg on', ms, '(bg subtracted) to create', tb
    try:
        tbavg(split,ms,tb,speed='fast',weight_mode='flat',datacolumn='corrected')
    except:
        print 'tbavg failed, probably b/c table',tb, 'is already open in the CASA cache - restart CASA to fix this problem'

    ## c) Extract dynspec ##
    
    print 'Saving dynspec to files in directory',tb+'.dynspec'
    spec = dyn_spec(tb)
    if os.path.exists(tb+'.dynspec'):
        os.system('rm -rf ' + tb + '.dynspec')
    saveTxt(spec,tb)


## d) Diagnostic plots of BG-subtracted smallms - time series & dirty image ##
'''
# create dirty images of smallms w/o BG
dirty_nobg=sb+'.dirty_nobg'
imfile = plotdir+dirty_nobg+'.jpg'
imfile2 = plotdir+dirty_nobg+'.small.jpg'
if len(glob(dirty_nobg+'.*'))>0:
    rmtables(dirty_nobg+'.*') # delete old image files
print 'Creating image during non-flaring scans with bg srcs subtracted:',imfile
clean(vis=smallms,imagename=dirty_nobg,scan=clean_scans,spw=clean_spws,imsize=imsize_d,cell=pixel_size,niter=500,stokes='IQUV')
blc = int(imsize_d/2)-32
trc = int(imsize_d/2)+32
zoom={'blc':[blc,blc],'trc':[trc,trc]}
raster={'file':dirty_nobg+'.image','colormap':'Greyscale 1','range':[0,rmsQUV*0.005]}
imview(raster=raster,out=imfile)
imview(raster=raster,out=imfile2,zoom=zoom)
'''
# plot tseries for smallms w/o BG
plotfile = plotdir + sb + '_tseries.png'
plotfile_im = plotdir + sb + '_im_tseries.png' # make sure there is no variation in Im(vis)
print 'Plotting time series for', sb,'with bg srcs removed:',plotfile, '(+im in filename for Im(vis))'
plotms(vis=smallms,xaxis='time',yaxis='real',ydatacolumn='corrected',avgchannel='128',avgbaseline=True,correlation='RR,LL',coloraxis='corr',plotfile=plotfile,overwrite=True,showgui=False)
plotms(vis=smallms,xaxis='time',yaxis='imag',ydatacolumn='corrected',avgchannel='128',avgbaseline=True,correlation='RR,LL',coloraxis='corr',plotfile=plotfile_im,overwrite=True,showgui=False)

