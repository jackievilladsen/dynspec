'''

calPband.py

Run this from the directory with the ms file in it: e.g., /data/jrv/15A-416/YZCMi/1/P

Run this after running flagPband.py.

General approach:
Work with low-res ms ([SBname].small.ms) to solve for calibration tables and inspect, then apply them to
full-res 1-spw ms ([SBname].1spw.ms) at the end.

Calibration steps (field used):
- load flux model (BP cal)
- delay (BP cal): K0
- phase-only gaincal (BP cal): G0
- bandpass (BP cal): B1  (B0 was initial BP cal in flagPband.py)
- polarization (BP cal)
---- cross-pol: Kc0
---- leakage: Df0
- redo delay and bandpass: K1, B2
- amp+phase gaincal on BP cal (G1)
- apply gains to target field and image (scans near BP cal only? or all?)
- selfcal on target field
- phase-only gaincal (gcal) then image pcal: G1
- image pcal, amp+phase selfcal (gcal): G2
- image science target, phase-only selfcal (target): G3
- apply all calibrations - K1,B2,Kc0,Df0,G3

- to do: add TEC corr and requantizer gains to initial processing in flagPband.py, also apply them now to YZCMi_1P
- see if TEC corr fixes the time-variable phases of BP calibrator

Print flag summaries throughout after applying each calibration step.
'''

from dynspec.pipeline_utils import *
import numpy as np

### DEFINE VARIABLES ###
# field names, file names, phasecen, refant

names = get_names() # file names
sb = names['sb']
name1spw = sb + '.1spw'
ms1spw = name1spw + '.ms'  # full-size ms produced by flagPband.py
smallms = sb + '.small.ms' # low-res ms avg'd in time and freq produced by flagPband.py
smallname = sb + '.small'
smallsrcname = sb + '.smallsrc'
smallsrcms = smallsrcname+'.ms'
gaintableGKB = [smallname+'.Gs0',smallname+'.K0',smallname+'.B1']
gaintableGKBKc = [smallname+'.Gs0',smallname+'.K0',smallname+'.B1',smallname+'.Kc0']
gaintableGKBKcDf = [smallname+'.Gs0',smallname+'.K0',smallname+'.B1',smallname+'.Kc0',smallname+'.Df0']
gaintableGBKcDf = [smallname+'.Gs0',smallname+'.B1',smallname+'.Kc0',smallname+'.Df0']
gaintableGKKcDf = [smallname+'.Gs0',smallname+'.K1',smallname+'.Kc0',smallname+'.Df0']
gaintableGKBKcDf2 = [smallname+'.Gs0',smallname+'.K1',smallname+'.B2',smallname+'.Kc0',smallname+'.Df0']
gaintableKBKcDf = [smallname+'.K1',smallname+'.B2',smallname+'.Kc0',smallname+'.Df0']
gaintableKBKcDfG = [smallname+'.K1',smallname+'.B2',smallname+'.Kc0',smallname+'.Df0',smallname+'.G1']
gaintableKBKcDfG2 = [smallname+'.K1',smallname+'.B2',smallname+'.Kc0',smallname+'.Df0',smallname+'.G2']

fields = get_fields(ms1spw,vishead) # field names
srcname = fields.get('src')
bpcal = fields.get('bpcal')
gcal = fields.get('gcal')

refanten = 'ea06' # reference antenna - ea06 is central ant on north arm - check plots after to make sure this is okay

phasecen=get_phasecen() # read phase center from file

plotdir = get_plot_dir()
'''

### DOCUMENT FLAGS BEFORE CALIBRATIONS ###

print '% flagged in',smallms,'before calibrations:'
summary_6 = field_flag_summary(smallms,flagdata)
np.save('summary6.npy',summary_6)


### LOAD FLUXCAL MODEL ###

# apply setjy (to account for source structure and flux/spectral index for BP cal)
# for now we are using L band model images since there are no images for P band (but SF12 flux standards do cover P band)
print 'Loading virtual flux model from',bpcal+'_L.im'
setjy(vis=smallms,field=bpcal,standard='Scaife-Heald 2012',model=bpcal+'_L.im')


### DELAY CALIBRATION: BPCAL ###

print 'Solving for delays...'
gaincal(vis=smallms,gaintype='K',field=bpcal,caltable=smallname+'.K0',refant=refanten,minsnr=3.0,parang=True)

# plot delays
plotcal(caltable=smallname+'.K0',xaxis='antenna',yaxis='delay',figfile=plotdir+'delay.png',showgui=False)
print 'Check plots dir for delay.png - all delays should be < few ns.'


### INITIAL BANDPASS CALIBRATION: BPCAL ###

# we will re-do bandpass calibration after a gaincal and polcal

print 'Solving for bandpass for',smallms,'- will redo bandpass solutions after gaincal and polcal'
bandpass(vis=smallms,caltable=smallname+'.B1',field=bpcal,refant=refanten,gaintable=[smallname+'.K0'])

plotcal(caltable=smallname+'.B1',xaxis='freq',yaxis='amp',subplot=441,iteration='antenna',showgui=False,figfile=plotdir+'BPamp.png')
plotcal(caltable=smallname+'.B1',xaxis='freq',yaxis='phase',subplot=441,iteration='antenna',showgui=False,figfile=plotdir+'BPphase.png')
print 'Check plots/BPamp.png and BPphase.png for bandpass solutions.'

# apply bandpass calibration to BP cal field
applycal(vis=smallms,field=bpcal,applymode='calflag',gaintable=[smallname+'.K0',smallname+'.B1'])

print '% flagged in',smallms,'after bandpass and delays applied to',bpcal+':'
summary_8 = field_flag_summary(smallms,flagdata)
np.save('summary8.npy',summary_8)

plotms(vis=smallms,field=bpcal,ydatacolumn='corrected',xaxis='freq',yaxis='amp',coloraxis='corr',correlation='XX,YY',plotrange=[0.2,0.5,0,150],plotfile=plotdir+'post_B1.png',showgui=False,overwrite=True,avgtime='1e8')
print 'Check plots/post_B1.png to confirm that bandpass-calibrated BPcal looks good.'


### GAINCAL: BPCAL ###

# phases of BP cal vary smoothly by tens of degrees over a single ~5-min BP cal obs
# I will apply these gains directly to the target field (skip the phase calibrator)

plotms(vis=smallms,field=bpcal,xaxis='time',yaxis='phase',ydatacolumn='corrected',coloraxis='baseline',antenna=refanten,correlation='XX',avgchannel='10000',plotfile=plotdir+'bpcal_phase_precal.png',showgui=False,overwrite=True)
plotms(vis=smallms,field=bpcal,xaxis='time',yaxis='amp',ydatacolumn='corrected',coloraxis='baseline',antenna=refanten,correlation='XX',avgchannel='10000',plotfile=plotdir+'bpcal_amp_precal.png',showgui=False,overwrite=True)
print 'Check plots/bpcal_phase[/amp]_precal.png to see how BP calibrator phases and amps vary with time (baselines with '+refanten+').'

gaincal(vis=smallms,gaintype='G',calmode='ap',field=bpcal,caltable=smallname+'.G0',refant=refanten,solint='int',minsnr=5,gaintable=[smallname+'.K0',smallname+'.B1'])

#plotcal(caltable=smallname+'.G0',xaxis='time',yaxis='phase',subplot=441,iteration='antenna',showgui=False,figfile=plotdir+'G0phase.png')
#plotcal(caltable=smallname+'.G0',xaxis='time',yaxis='amp',subplot=441,iteration='antenna',showgui=False,figfile=plotdir+'G0amp.png')
#print 'Check plots/G0amp.png and G0phase.png for BP calibrator gaincal solutions.'

# smooth BPcal gains so that gain solutions can be applied to target field even if nearest gain is flagged
# smoothtime = 60 min --> will have ~ same gain for whole time
smoothcal(vis=smallms,tablein=smallname+'.G0',caltable=smallname+'.Gs0', smoothtype='median', smoothtime = 60.*60.)
plotcal(caltable=smallname+'.Gs0', xaxis='time', yaxis='amp', iteration='antenna',subplot=441,showgui=False,figfile=plotdir+'Gs0amp.png')
plotcal(caltable=smallname+'.Gs0', xaxis='time', yaxis='phase', iteration='antenna',subplot=441,showgui=False,figfile=plotdir+'Gs0phase.png')
print 'Check plots/Gs0amp[/phase].png for smoothed BP calibration gains - should be one constant value per antenna per poln.'

print 'Applying delays, bandpass, and smoothed gains to',bpcal
applycal(vis=smallms,field=bpcal,applymode='calflag',gaintable=gaintableGKB)

plotms(vis=smallms,field=bpcal,xaxis='time',yaxis='phase',ydatacolumn='corrected',coloraxis='baseline',antenna=refanten,correlation='XX',avgchannel='10000',plotfile=plotdir+'bpcal_phase_postcal.png',showgui=False,overwrite=True)
plotms(vis=smallms,field=bpcal,xaxis='time',yaxis='amp',ydatacolumn='corrected',coloraxis='baseline',antenna=refanten,correlation='XX',avgchannel='10000',plotfile=plotdir+'bpcal_amp_postcal.png',showgui=False,overwrite=True)
print 'Check plots/bpcal_phase[/amp]_postcal.png to see how BP calibrator gains vary with time after smoothed gaincal.'
# should be < few degrees

# show % flagged after gaincal solns applied to BP cal
print '% flagged in',smallms,'after gaincal on BP cal (should only affect BPcal field if any):'
summary_9 = field_flag_summary(smallms,flagdata)
np.save('summary9.npy',summary_9)


### POLARIZATION CALIBRATION: 3C147 ###

# Perley et al. 2013: % (linear?) polarization @ 1.05 GHz, rising w/ frequency --> P band should be lower
pol_dict = {'3C48':'0.3%','3C138':'5.6%','3C147':'<0.05%','3C286':'8.6%'}
p = pol_dict.get(bpcal,'[unknown source name]')
print 'Perley+13:',bpcal,'1.05-GHz linear polarization is',p+', error on P-band pol cal should be less than this.'

## Cross-delay calibration (Kc0) ##

# cross-delay calibration is the delay between the two polarizations of the refant
gaincal(vis=smallms,field=bpcal,caltable=smallname+'.Kc0',gaintype='KCROSS',refant=refanten,parang=True,gaintable=gaintableGKB)

plotcal(caltable=smallname+'.Kc0',xaxis='antenna',yaxis='delay',showgui=False,figfile=plotdir+'crossdelay.png')
print 'Check plots/crossdelay.png to see delay between 2 polns of refant', refanten

## Leakage calibration ##

# leakage calibration removes freq-dependent leakage from X to Y and vice versa for each antenna ("D terms")

# poltype='Df': assume calibrator is unpolarized; solve for each channel
# preavg=1. - time interval for applying parallactic angle corrections
polcal(vis=smallms,field=bpcal,poltype='Df',caltable=smallname+'.Df0',refant=refanten,preavg=1.0,gaintable=gaintableGKBKc)

plotcal(caltable=smallname+'.Df0',xaxis='freq',yaxis='amp',iteration='antenna',subplot=441,plotrange=[200,500,0,1],showgui=False,figfile=plotdir+'leakage_preflag.png')
print 'Check plots/leakage_preflag.png for D terms (freq-dpdt leakage between polns for each antenna).  Should be <<1.'

# flag D terms (> 0.4? or sudden spikes?)
flagdata(vis=smallname+'.Df0',mode='clip',clipminmax=[-0.1,0.4],datacolumn='CPARAM')

plotcal(caltable=smallname+'.Df0',xaxis='freq',yaxis='amp',iteration='antenna',subplot=441,plotrange=[200,500,0,1],showgui=False,figfile=plotdir+'leakage_clip.png')
plotcal(caltable=smallname+'.Df0',xaxis='freq',yaxis='amp',iteration='antenna',subplot=441,showgui=False,figfile=plotdir+'leakage_clip_zoom.png')
print 'Check plots/leakage_clip[_zoom].png for D terms after autoflag (clip values > 0.4).  Should be <<1.'

print 'Applying poln caltables to', bpcal
applycal(vis=smallms,field=bpcal,applymode='calflag',gaintable=gaintableGKBKcDf)

print '% flagged in',smallms,'after applying poln calibration tables:'
summary_10 = field_flag_summary(smallms,flagdata)
np.save('summary10.npy',summary_10)

# inspect amplitude of cross-pols (XY,YX) for BP cal
plotms(vis=smallms,field=bpcal,xaxis='freq',yaxis='amp',ydatacolumn='corrected',coloraxis='corr',avgtime='1e8',avgscan=True,plotrange=[0.2,0.5,0,70],plotfile=plotdir+'bpcal_cross_corr.png',showgui=False,overwrite=True)
print 'Check plots/bpcal_cross_corr.png to verify no power in cross polns after polcal.  XX and YY are purple and orange, XY and YX are black and pink (respectively).'


### REDO BANDPASS AND DELAY CALIBRATION ###

print 'Solving for new delays w/ poln cal applied...'
# delay calibration
gaincal(vis=smallms,gaintype='K',field=bpcal,caltable=smallname+'.K1',refant=refanten,minsnr=3.0,parang=True,gaintable=gaintableGBKcDf)
plotcal(caltable=smallname+'.K1',xaxis='antenna',yaxis='delay',figfile=plotdir+'delay1.png',showgui=False)
print 'Check plots/delay1.png for updated delay solutions after polcal.'

print 'Solving for new bandpass w/ poln cal applied...'
# bandpass calibration
bandpass(vis=smallms,caltable=smallname+'.B2',field=bpcal,refant=refanten,gaintable=gaintableGKKcDf)
plotcal(caltable=smallname+'.B2',xaxis='freq',yaxis='amp',subplot=441,iteration='antenna',showgui=False,figfile=plotdir+'BPamp2.png')
plotcal(caltable=smallname+'.B2',xaxis='freq',yaxis='phase',subplot=441,iteration='antenna',showgui=False,figfile=plotdir+'BPphase2.png')
print 'Check plots/BPamp2.png and BPphase2.png for updated bandpass solutions after polcal.'

# flag bandpass caltable based on phases - for some reason a few channels have negative values for Real(bandpass)
#  - can catch this by flagging Real(bandpass) or phase(bandpass) ('ARG_ALL'=phase of all polns)
flagdata(vis=smallname+'.B2',mode='rflag',winsize=3,correlation='ARG_ALL',datacolumn='CPARAM',display='report')
plotcal(caltable=smallname+'.B2',xaxis='freq',yaxis='amp',subplot=441,iteration='antenna',showgui=False,figfile=plotdir+'BPamp3.png')
plotcal(caltable=smallname+'.B2',xaxis='freq',yaxis='phase',subplot=441,iteration='antenna',showgui=False,figfile=plotdir+'BPphase3.png')
print 'Check plots/BPamp3.png and BPphase3.png for updated B2 after rflag on B2 phase.'

# apply updated calibration to all fields
applycal(vis=smallms,applymode='calflag',gaintable=gaintableGKBKcDf2)

print '% flagged in',smallms,'after updated caltables (based on',bpcal+') applied to all fields:'
summary_11 = field_flag_summary(smallms,flagdata)
np.save('summary11.npy',summary_11)

plotms(vis=smallms,field=bpcal,ydatacolumn='corrected',xaxis='freq',yaxis='amp',coloraxis='corr',correlation='XX,YY',plotrange=[0.2,0.5,0,150],plotfile=plotdir+'post_K1B2.png',showgui=False,overwrite=True,avgtime='1e8')
print 'Check plots/post_K1B2.png to confirm that BPcal looks good after updated bandpass and delay after polcal.'
# tutorial re-solves for polcal, not delay - should I change this?


### FLAG w/ final BP calibration ###

# run autoflagger one more time with final calibrations based on BPcal
flagdata(vis=smallms,mode='rflag',datacolumn='corrected',winsize=5,display='report')

print '% flagged in',smallms,'after rflag with K1 and B2 applied:'
summary_12 = field_flag_summary(smallms,flagdata)
np.save('summary12.npy',summary_12)

plotms(vis=smallms,field=bpcal,ydatacolumn='corrected',xaxis='freq',yaxis='amp',coloraxis='corr',correlation='XX,YY',plotrange=[0.2,0.5,0,150],plotfile=plotdir+'post_K1B2_rflag.png',showgui=False,overwrite=True,avgtime='1e8')
print 'Check plots directory for post_K1B2_rflag.png to confirm that polarization- and bandpass-calibrated BPcal looks good after rflag.'
'''

### IMAGE TARGET FIELD NEAR BPCAL SCAN ###

# create new ms w/ calibrated target field data
split(vis=smallms,field=srcname,outputvis=smallsrcms)

im0='pix8192_3arcsec_n1'

#clean(vis=smallsrcms,imagename=im0,niter=10000,interactive=True,imsize=8192,cell='3arcsec',usescratch=True)
#ft(vis=smallsrcms,model=im0+'.model',usescratch=True) # b/c I forgot to do usescratch=True for clean

### PHASE-ONLY SELFCAL ON TARGET FIELD ###

# try single solution for whole scan
gaincal(vis=smallsrcms,refant=refanten,calmode='p',caltable=smallsrcname+'.ScG0',solint='60s')
plotcal(caltable=smallsrcname+'.ScG0',xaxis='time',yaxis='phase',iteration='antenna',subplot=441,showgui=False,figfile=plotdir+'ScG0phase.png')
applycal(vis=smallsrcms,gaintable=smallsrcname+'.ScG0')

### IMAGE AGAIN - QUALITY IMPROVED?? ###

im1 = 'pix8192_3arcsec_n1_v2'
clean(vis=smallsrcms,imagename=im1,mask=im0+'.mask',niter=3000,interactive=True,imsize=8192,cell='3arcsec',usescratch=True)

im2 = 'pix8192_3arcsec_n1_auto'
clean(vis=smallsrcms,imagename=im2,mask=im1+'.mask',niter=3000,cyclefactor=5,imsize=8192,cell='3arcsec')

im3 = 'pix8192_3arcsec_n2_auto'
clean(vis=smallsrcms,imagename=im3,mask=im1+'.mask',niter=3000,cyclefactor=5,nterms=2,imsize=8192,cell='3arcsec')


# are the coordinates accurate? try imaging the pcal w/ the BPcal gains applied
# I already applied the gains to smallms so now just need to image
clean(vis=smallms,imagename='pcal',field=gcal,niter=0,imsize=1024,cell='3arcsec') # simple params - just want to see if in center
clean(vis=smallms,imagename='pcal21',field=gcal,scan='21',niter=0,imsize=1024,cell='3arcsec') # simple params - just want to see if in center
clean(vis=smallms,imagename='pcal27',field=gcal,scan='27',niter=0,imsize=1024,cell='3arcsec') # simple params - just want to see if in center


'''
### GAIN CALIBRATION ###
# approach: initial phase-only gaincal w/ no model, followed by imaging and selfcal on gaincal field

## Initial phase-only gaincal ##
gaincal(vis=smallms,caltable=smallname+'.G1',field=bpcal+','+gcal,gaintable=gaintableKBKcDf,solint='inf',gaintype='G',calmode='p',refant=refanten)
# solint='inf' --> one solution per scan

# plot gain solutions - phase only (since that's all we solved for)
plotcal(caltable=smallname+'.G1',xaxis='time',yaxis='phase',subplot=441,iteration='antenna',showgui=False,figfile=plotdir+'G2phase.png')

# apply to gcal field so we can plot and image
applycal(vis=smallms,field=gcal,applymode='calflag',gaintable=gaintableKBKcDfG)

# should I break it back up into multiple spw's for gaincal? come back to this
# inspecting spectrum per scan in plotms, it looks mostly good - one scan has one bad spw, see if flagging gets this,
#    otherwise we may need to split it up

## Image gain calibrator field ##

# SKIPPED THIS SECTION #-----------------------------------

# Get imaging parameters
pixel_size,imsize = im_params(smallms)

# measure RMS of Stokes Q,U,V dirty images during non-flaring scans
dirtyim='dirty'
if len(glob(dirtyim+'.*'))>0:
    rmtables(dirtyim + '.*') # have to delete old image files (esp. model) so clean will start from scratch
print 'Creating dirty image of',gcal,'- will use Stokes QUV RMS to set threshold (3*RMS) for cleaning'
clean(vis=smallms,imagename=dirtyim,field=gcal,imsize=imsize,cell=pixel_size,niter=0,stokes='IQUV') # create dirty image
rmsQ = imstat(imagename=dirtyim+'.image',stokes='Q')['sigma']
rmsU = imstat(imagename=dirtyim+'.image',stokes='U')['sigma']
rmsV = imstat(imagename=dirtyim+'.image',stokes='V')['sigma']
rmsQUV = min([rmsQ,rmsU,rmsV])[0]*1000
print 'RMS of Stokes Q,U,V dirty image (min of the three):',rmsQUV, 'mJy'
threshold = str(rmsQUV * 3) + 'mJy'

# clean smallms with nterms=2 since we have large fractional bandwidth (may need nterms=3?)
gcalim = gcal+'.v0'
if len(glob(gcalim+'.*'))>0:
    rmtables(gcalim + '.*') # have to delete old image files (esp. model) so clean will start from scratch
print 'Imaging',gcal,'field with nterms=2 - image file', gcalim+'.image.tt0'
clean(vis=smallms,imagename=gcalim,field=gcal,nterms=2,imsize=imsize,cell=pixel_size,niter=3000,threshold=threshold,cyclefactor=5,interactive=True,usescratch=True)
# high cyclefactor forces major cycles (going back to visibilities) more often - important if PSF poorly known
# can also try changing cycle gain but in other bands the default gain of 0.1 has been fine

# add command to save image of field as .png

# use ft to populate MODEL_DATA column of ms from model images (.tt0 and .tt1) - don't need to do this if had usescratch=True for clean
#ft(vis=smallms,nterms=2,model=[gcal+'.v0.model.tt0',gcal+'.v0.model.tt1'],usescratch=True,field=gcal)

## Amp+phase gaincal using clean model of pcal field ##

gaincal(vis=smallms,caltable=smallname+'.G2',field=gcal,gaintable=gaintableKBKcDf,solint='120s',gaintype='G',calmode='ap',refant=refanten)
# I think it will use MODEL_DATA column by default?

# plot gains and compare to previous phase-only model-less gaincal
plotcal(caltable=smallname+'.G1',xaxis='time',yaxis='phase',subplot=441,iteration='antenna',showgui=False,figfile=plotdir+'G1phase.png')
plotcal(caltable=smallname+'.G2',xaxis='time',yaxis='phase',subplot=441,iteration='antenna',showgui=False,figfile=plotdir+'G2phase.png')
plotcal(caltable=smallname+'.G2',xaxis='time',yaxis='amp',subplot=441,iteration='antenna',showgui=False,figfile=plotdir+'G2amp.png')
# phases are a little different; amps vary from ~0.8-1.2 --> good thing we did amp gaincal

# END OF SKIPPED PORTION # --------------------------


### APPLY CALIBRATIONS TO TARGET ###

# then image target field

# then selfcal on target 

# if this works then we can just toss the pcal scans
# shame i spent time on it but it will save time for data reduction and future observations

gaintableKBKcDfG2 = [smallname+'.K0',smallname+'.B1',smallname+'.Kc0',smallname+'.Df0',smallname+'.G2']
applycal(vis=smallms,field=srcname,applymode='calflagstrict',gaintable=gaintableKBKcDfG2)
# do I need to define interp=nearest for some of these cal tables?  linear is default.
# backup flag file: before_applycal_10 (don't think it flagged any additional data though)

# inspect time series in plotms
# it looks like there is a RCP flare - peak flux of ~50 mJy
# Stokes V = i(XY-YX)
# need to do a bit more RFI flagging and source subtraction to improve SNR, remove outliers from time series



### CLEAN UV CETI FIELD ###

# use same params as when imaging gaincal field
clean(vis=smallms,imagename=srcname+'.v0',field=srcname,niter=1000,nterms=2,cell=['1.5arcsec'],imsize=[8192,8192],usescratch=True)
# clean model contains a lot of sidelobes...  use threshold and/or interactive clean to avoid this

# let's do an interactive clean
clean(vis=smallms,imagename=srcname+'.v1',field=srcname,niter=1000,nterms=2,cell=['1.5arcsec'],imsize=[8192,8192],usescratch=True,interactive=True)
# based on this, sources at the edge of the field move more than others; do I need to use peel?
#  the procedure for peel is: shift phase center to brightest source, selfcal, clean brightest source to level
#      of sidelobes of other sources, subtract brightest source, shift phase center to new brightest source, repeat
#  Since in the end we want a good phase calibration for UV Cet (which is near the center of the field, we should
#      probably leave a bright source near UV Cet for last and use selfcal on that source as the phase calibration for
#      UV Cet
#  I won't try peel yet - for now just selfcal based on the clean model of the whole field, then see what happens


## Selfcal on target field ##

# start with phase-only gaincal
gaincal(vis=smallms,caltable=smallname+'.G3',field=srcname,gaintable=gaintableKBKcDfG2,solint='inf',gaintype='G',calmode='p',refant=refanten)

# apply calibration (no flagging)
gaintableKBKcDfG23 = [smallname+'.K0',smallname+'.B1',smallname+'.Kc0',smallname+'.Df0',smallname+'.G2',smallname+'.G3']
applycal(vis=smallms,field=srcname,applymode='calonly',gaintable=gaintableKBKcDfG23)

# re-image to see if this is an improvement (just a dirty image)
clean(vis=smallms,imagename=srcname+'.v2',field=srcname,niter=0,nterms=2,cell=['1.5arcsec'],imsize=[8192,8192])
# perhaps a bit improved, but not fixed - some sources are stretched out (by phase errors?) and some are not
# peeling seems like a lot of work so now just clean it as well as I can

# interactive clean
clean(vis=smallms,imagename=srcname+'.v3',field=srcname,niter=1000,nterms=2,cell=['1.5arcsec'],imsize=[8192,8192],usescratch=True,interactive=True)

'''
