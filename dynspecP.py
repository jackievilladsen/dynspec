# run this pipeline from the directory for the particular observation you want to process
# (e.g., /data/jrv/15A-416/YZCMi/1/P)
#
# This pipeline should be run after flagPband.py and calPband.py.  It will look for ms1spw (e.g., YZCMi_1P.1spw.ms),
#  where the calibrated data column should include the calibrated target field.
#
# If you want to change default parameters for pipeline, create dictionary pipe_params before running this
#   script with execfile.  Default for pipe_params is: {'overwrite':False,'makebig':False,'interactive':False}

# Steps:
# 1) Definitions - parameters, field names, filenames, phasecen
# 2) ms file creation
#    a) Use mstransform to make srcms (calibrated science target only w/ phase center on star, single spw)
# 3) Create BG src model from srcms
#    a) Plot raw time series (will show sidelobes of BG srcs) and read non-flaring scans/spws from file
#    b) Image and clean non-flare scans/spws w/ nterms>1 (2? 3?)
#    c) Mask central star out of model image to create BG model image
# 4) Subtract BG srcs and avg over baselines
#    a) FT BG model image then uvsub
#    b) Average over all baselines
#    c) Extract dynamic spectrum from baseline-avg'd ms
#    d) Diagnostic plots of BG-subtracted srcms - time series & dirty image

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
ms1spw = sb + '.1spw.ms'
srcms = sb + '.star.ms'
srctbavg = sb + '.tbavg.ms'
smallim = sb + '.smallim'

fields = get_fields(ms1spw,vishead)
srcname = fields.get('src')
bpcal = fields.get('bpcal')
gcal = fields.get('gcal')
phasecen=get_phasecen()


### 2) MS FILE CREATION ###

## a) Use mstransform to make srcms (calibrated science target only w/ phase center on star, single spw) ##

# use mstransform to create src-only ms from calibrated data column of ms1spw
if not os.path.exists(srcms) or overwrite:  # skip this step if srcms already exists
    if os.path.exists(srcms):
        rmtables(srcms)
    if os.path.exists(srcms+'.flagversions'):
        os.system('rm -rf '+srcms+'.flagversions')
    print 'creating src-only ms', srcms, 'from calibrated data column of', ms1spw
    mstransform(vis=ms1spw,outputvis=srcms,field=srcname,keepflags=False)
    print 'shifting phase center of', srcms, 'to', phasecen
    fixvis(vis=srcms,outputvis=srcms,phasecenter=phasecen)
else:
    print srcms, 'already exists, now deleting MODEL_DATA and CORRECTED_DATA and assuming phase center already shifted'
    # erase corrected data column to remove source subtraction from previous pipeline runs
    clearcal(srcms)

# backup src-only ms as tar file (in case calibrated pipeline ms and srcms get corrupted)
tarfile = srcms + '.tar'
if not os.path.exists(tarfile) or overwrite: # skip this step if tarfile already exists
    print 'saving',tarfile,'with fixvis applied'
    os.system('tar -cf ' + tarfile + ' ' + srcms)
else:
    print tarfile, 'already exists, skipping the creation of this file'



### 3) CREATE BACKGROUND SOURCE MODEL FROM SRCMS ###

## a) Plot raw time series (will show sidelobes of BG srcs) and read non-flaring scans/spws from file ##
print 'Plotting time series before source subtraction... a circularly polarized stellar burst should be positive in real part of XX and YY (orange and purple), and opposite signs in imaginary part of XY and YX (pink and black).'
plotfile = plotdir + sb + '_dirty_tseries.png'
implotfile = plotdir + sb + '_dirty_tseries_im.png'
plotms(vis=srcms,xaxis='scan',yaxis='real',avgchannel='2048',avgtime='1e8',avgbaseline=True,coloraxis='corr',plotfile=plotfile,overwrite=True,showgui=False)
plotms(vis=srcms,xaxis='scan',yaxis='imag',avgchannel='2048',avgtime='1e8',avgbaseline=True,coloraxis='corr',plotfile=implotfile,overwrite=True,showgui=False)
if interactive:
    raw_input('Check plots/'+sb+'_dirty_tseries[_im].png to identify flare-free scans, add them to dynspec/clean_scans.txt, then hit Return.')
clean_scans,clean_spws = get_clean_scans()
print 'Scans and spws used as flare-free for',sb,'(blank is all): scan=\''+clean_scans+'\', spw=\''+clean_spws+'\''


## b) Image and clean srcms in non-flare scans/spws w/ nterms=2 ##

# get parameters for imaging
pixel_size,imsize = im_params(srcms)

if imsize > 10000:
    print 'Warning: imsize>10000, setting to imsize=10000 (may cut off full FOV)'
    imsize = 10000

# measure RMS of Stokes Q&U dirty image during non-flaring scans
#  (when subtracting BG srcs later we assume unpolarized, so this will prevent overcleaning polarized BG src)
dirtyim='dirty'
if len(glob(dirtyim+'.*'))>0:
    rmtables(dirtyim + '.*') # have to delete old image files (esp. model) so clean will start from scratch
print 'Creating dirty image during non-flaring scans - will use Stokes Q&U RMS to set threshold (3*RMS) for cleaning'
imsize_d = min(imsize,1024)
clean(vis=srcms,imagename=dirtyim,scan=clean_scans,spw=clean_spws,imsize=imsize_d,cell=pixel_size,niter=0,stokes='IQUV') # create dirty image
rmsQ = imstat(imagename=dirtyim+'.image',stokes='Q')['sigma']
rmsU = imstat(imagename=dirtyim+'.image',stokes='U')['sigma']
rmsV = imstat(imagename=dirtyim+'.image',stokes='V')['sigma']
rmsQUV = min([rmsQ,rmsU,rmsV])[0]*1000
print 'RMS of Stokes Q,U,V dirty image (min of the three):',rmsQUV, 'mJy'
threshold = str(rmsQUV * 3) + 'mJy'

# clean srcms using only the non-flaring scans with nterms=2 since we have large fractional bandwidth (may need nterms=3?)
if len(glob(smallim+'.*'))>0:
    rmtables(smallim + '.*') # have to delete old image files (esp. model) so clean will start from scratch
print 'Imaging non-flaring scans in',srcms,'with nterms=2 - image file', smallim+'.image.tt0'
interactive=False
clean(vis=srcms,imagename=smallim,scan=clean_scans,spw=clean_spws,nterms=2,imsize=imsize,cell=pixel_size,niter=3000,threshold=threshold,cyclefactor=5,interactive=interactive)
# high cyclefactor is CRITICAL to successful clean if not using clean boxes, otherwise sidelobes end up in clean model - this means that the image plane
#   PSF used by clean is not very accurate (b/c of nterms=2? b/c of subarray --> fewer ants?)
# default gain of 0.1 is fine - higher gain (0.3) results in more spurious clean components


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
circle_cen = 'circle[['+ncen+'pix,'+ncen+'pix],10pix]' # a 10-pixel radius(?) circle in the image center
print 'Excluding image center region', circle_cen, 'from model of bg srcs'
makemask(mode='copy',inpimage=model0,output=cen_mask,inpmask=circle_cen,overwrite=True) # create mask with ones in center
immath(imagename=[model0,cen_mask],outfile=bgmodel0,expr='IM0*(1-IM1)')
immath(imagename=[model1,cen_mask],outfile=bgmodel1,expr='IM0*(1-IM1)')
# bgmodel0,1 are the same as model0,1 but with the center pixels removed (radius 10 pixels)

'''
### 4) SUBTRACT BG AND CREATE DYNSPEC ###

mslist = [srcms]
tblist = [srctbavg]

for (ms,tb) in zip(mslist,tblist):
    
    ## a) FT background model image then uvsub ##
    clearcal(ms)
    # ft to populate MODEL_DATA column with clean model
    print 'subtracting bg src model from', ms, '(ft + uvsub)'
    ft(vis=ms,nterms=2,model=[bgmodel0,bgmodel1],usescratch=True)
    ft(vis=ms,nterms=2,model=[bgmodel0,bgmodel1],usescratch=True)
    
    # uvsub to subtract model column from data column and put it in corrected column
    uvsub(vis=ms)
    #uvsub(vis=ms)
    
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
    


## d) Diagnostic plots of BG-subtracted srcms - time series & dirty image ##

# create dirty images of srcms w/o BG
dirty_nobg=sb+'.dirty_nobg'
imfile = plotdir+dirty_nobg+'.jpg'
imfile2 = plotdir+dirty_nobg+'.small.jpg'
if len(glob(dirty_nobg+'.*'))>0:
    rmtables(dirty_nobg+'.*') # delete old image files
print 'Creating dirty image during non-flaring scans with bg srcs subtracted:',imfile
clean(vis=srcms,imagename=dirty_nobg,scan=clean_scans,spw=clean_spws,imsize=imsize_d,cell=pixel_size,niter=0,stokes='IQUV')
blc = int(imsize_d/2)-64
trc = int(imsize_d/2)+64
zoom={'blc':[blc,blc],'trc':[trc,trc]}
raster={'file':dirty_nobg+'.image','colormap':'Greyscale 1','range':[0,rmsQUV*0.005]}
imview(raster=raster,out=imfile)
imview(raster=raster,out=imfile2,zoom=zoom)

# plot tseries for srcms w/o BG
plotfile = plotdir + sb + '_tseries.png'
plotfile_im = plotdir + sb + '_im_tseries.png' # make sure there is no variation in Im(vis)
print 'Plotting time series for', sb,'with bg srcs removed:',plotfile, '(+im in filename for Im(vis))'
plotms(vis=srcms,xaxis='scan',yaxis='real',ydatacolumn='corrected',avgchannel='2048',avgtime='1e8',avgbaseline=True,coloraxis='corr',plotfile=plotfile,overwrite=True,showgui=False)
plotms(vis=srcms,xaxis='scan',yaxis='imag',ydatacolumn='corrected',avgchannel='2048',avgtime='1e8',avgbaseline=True,coloraxis='corr',plotfile=plotfile_im,overwrite=True,showgui=False)

'''
