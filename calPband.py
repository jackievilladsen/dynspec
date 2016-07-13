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
gaintableGKB = [smallname+'.G0',smallname+'.K0',smallname+'.B1']
gaintableGKBKc = [smallname+'.G0',smallname+'.K0',smallname+'.B1',smallname+'.Kc0']
gaintableGKBKcDf = [smallname+'.G0',smallname+'.K0',smallname+'.B1',smallname+'.Kc0',smallname+'.Df0']
gaintableGBKcDf = [smallname+'.G0',smallname+'.B1',smallname+'.Kc0',smallname+'.Df0']
gaintableGKKcDf = [smallname+'.G0',smallname+'.K1',smallname+'.Kc0',smallname+'.Df0']
gaintableKBKcDf = [smallname+'.K1',smallname+'.B2',smallname+'.Kc0',smallname+'.Df0']
gaintableKBKcDfG = [smallname+'.K1',smallname+'.B2',smallname+'.Kc0',smallname+'.Df0',smallname+'.G1']
gaintableKBKcDfG2 = [smallname+'.K1',smallname+'.B2',smallname+'.Kc0',smallname+'.Df0',smallname+'.G2']

fields = get_fields(msfile,vishead) # field names
srcname = fields.get('src')
bpcal = fields.get('bpcal')
gcal = fields.get('gcal')

refanten = 'ea06' # reference antenna - ea06 is central ant on north arm - check plots after to make sure this is okay

phasecen=get_phasecen() # read phase center from file

plotdir = get_plot_dir()


### DOCUMENT FLAGS BEFORE CALIBRATIONS ###

print '% flagged in',smallms,'before calibrations:'
summary_6 = field_flag_summary(smallms,flagdata)
np.save('summary6.npy',summary_6)
# 5-10% lower than before avg'ing (hi-res ms) --> most flags covered multiple times/freqs


### LOAD FLUXCAL MODEL ###

# apply setjy (to account for source structure and flux/spectral index for BP cal)
# for now we are using L band model images since there are no images for P band (but SF12 flux standards do cover P band)
setjy(vis=smallms,field=bpcal,standard='Scaife-Heald 2012',usescratch=True,model=bpcal+'_L.im')


### DELAY CALIBRATION: BPCAL ###

gaincal(vis=smallms,gaintype='K',field=bpcal,caltable=smallname+'.K0',refant=refanten,minsnr=3.0,parang=True)
# delay soln was successful

# plot delays
plotcal(caltable=smallname+'.K0',xaxis='antenna',yaxis='delay',figfile=plotdir+'delay.png',showgui=False)
print 'Check plots dir for delay.png - all delays should be < few ns.'

# apply delay table - doing this so we can track how it affects flags
applycal(vis=smallms,applymode='calflag',gaintable=[smallname+'.K0'],interp=['linear'])

# show % flagged after delay solns applied
print '% flagged in',smallms,'after delay cal (no flags at this step if all delay solns were successful):'
summary_7 = field_flag_summary(smallms,flagdata)
np.save('summary7.npy',summary_7)


### PHASE-ONLY GAINCAL: BPCAL ###

# phases of BP cal vary smoothly by tens of degrees over a single ~5-min BP cal obs
# this could be due to ionospheric variation or bright BG srcs in BP cal field
# assuming it is ionospheric, then we should selfcal on BP cal before bandpass solutions

plotms(vis=smallms,field=bpcal,xaxis='time',yaxis='phase',ydatacolumn='corrected',coloraxis='baseline',antenna=refanten,correlation='XX',avgchannel='10000',plotfile=plotdir+'bpcal_phase_precal.png',showgui=False,overwrite=True)
print 'Check plots dir bpcal_phase_precal.png to see how BP calibrator phases vary with time.'

gaincal(vis=smallms,gaintype='G',calmode='p',field=bpcal,caltable=smallname+'.G0',refant=refanten,solint='int',minsnr=5,gaintable=[smallname+'.K0'])

plotcal(caltable=smallname+'.G0',xaxis='time',yaxis='phase',subplot=441,iteration='antenna',showgui=False,figfile=plotdir+'G0.png')
print 'Check plots dir G0.png for BP calibrator phase-only gaincal solutions.'

applycal(vis=smallms,field=bpcal,applymode='calflag',gaintable=[smallname+'.K0',smallname+'.G0'],interp=['linear'])

plotms(vis=smallms,field=bpcal,xaxis='time',yaxis='phase',ydatacolumn='corrected',coloraxis='baseline',antenna=refanten,correlation='XX',avgchannel='10000',plotfile=plotdir+'bpcal_phase_postcal.png',showgui=False,overwrite=True)
print 'Check plots dir bpcal_phase_postcal.png to see how BP calibrator phases vary with time after phase-only gaincal.'
# should be < few degrees

# show % flagged after gaincal solns applied to BP cal
print '% flagged in',smallms,'after phase-only gaincal on BP cal (should only affect BPcal field if any):'
summary_8 = field_flag_summary(smallms,flagdata)
np.save('summary8.npy',summary_8)


### BANDPASS CALIBRATION: BPCAL ###

bandpass(vis=smallms,caltable=smallname+'.B1',field=bpcal,refant=refanten,gaintable=[smallname+'.G0',smallname+'.K0'])

plotcal(caltable=smallname+'.B1',xaxis='freq',yaxis='amp',subplot=441,iteration='antenna',showgui=False,figfile=plotdir+'BPamp.png')
plotcal(caltable=smallname+'.B1',xaxis='freq',yaxis='phase',subplot=441,iteration='antenna',showgui=False,figfile=plotdir+'BPphase.png')
print 'Check plots dir BPamp.png and BPphase.png for bandpass solutions.'

# apply bandpass calibration to all fields
applycal(vis=smallms,applymode='calflag',gaintable=[smallname+'.G0',smallname+'.K0',smallname+'.B1'],interp=['nearest','nearest','nearest'])

print '% flagged in',smallms,'after bandpass applied to all fields:'
summary_9 = field_flag_summary(smallms,flagdata)
np.save('summary9.npy',summary_9)

plotms(vis=smallms,field=bpcal,ydatacolumn='corrected',xaxis='freq',yaxis='amp',coloraxis='corr',correlation='XX,YY',plotrange=[0.2,0.5,0,150],plotfile=plotdir+'post_B1.png',showgui=False,overwrite=True,avgtime='1e8')
print 'Check plots directory for post_B1.png to confirm that bandpass-calibrated BPcal looks good.'
# There are a few channels where phase or amp of BP solns looks bad, or where BPcal amp looks bad, but
#  another round of RFI auto-flagging should catch those hopefully.


### FLAG BASED ON BP CALIBRATOR ###

# run autoflagger one more time with better-calibrated data
flagdata(vis=smallms,mode='rflag',datacolumn='corrected',winsize=5,display='report')
# backup: flagdata_3

print '% flagged in',smallms,'after rflag with bandpass applied:'
summary_10 = field_flag_summary(smallms,flagdata)
np.save('summary10.npy',summary_10)

plotms(vis=smallms,field=bpcal,ydatacolumn='corrected',xaxis='freq',yaxis='amp',coloraxis='corr',correlation='XX,YY',plotrange=[0.2,0.5,0,150],plotfile=plotdir+'post_rflag3.png',showgui=False,overwrite=True,avgtime='1e8')
print 'Check plots directory for post_rflag3.png to confirm that bandpass-calibrated BPcal looks good after rflag.'


### POLARIZATION CALIBRATION: 3C147 ###

# Perley et al. 2013: % (linear?) polarization @ 1.05 GHz, rising w/ frequency --> P band should be lower
# 3C48: 0.3%
# 3C138: 5.6%
# 3C147: <0.05%
# 3C286: 8.6%
# We can assume few percent or less lin pol in these sources --> cross-delay calibration will be accurate
#   to within a few percent (few hundredths of a percent for 3C147, the calibrator I use most often).

## Cross-delay calibration ##

# cross-delay calibration is the delay between the two polarizations of the refant

# name of cross-delay table: smallname + '.Kc0'
gaincal(vis=smallms,field=bpcal,caltable=smallname+'.Kc0',gaintype='KCROSS',refant=refanten,parang=True,gaintable=gaintableGKB)

plotcal(caltable=smallname+'.Kc0',xaxis='antenna',yaxis='delay',showgui=False,figfile=plotdir+'crossdelay.png')
print 'Check plots/crossdelay.png to see delay between 2 polns of refant', refanten

## Leakage calibration ##

# leakage calibration removes freq-dependent leakage from X to Y and vice versa for each antenna ("D terms")

# inspect amplitude of cross-pols (XY,YX) for BP cal
plotms(vis=smallms,field=bpcal,xaxis='freq',yaxis='amp',ydatacolumn='corrected',coloraxis='corr',iteraxis='baseline',avgtime='1e8',avgscan=True,plotrange=[0.2,0.5,0,70],gridrows=3,gridcols=3,plotfile=plotdir+'bpcal_cross.png',showgui=False,exprange='all',overwrite=True)
print 'Check plots/bpcal_cross.png to see strength of cross polns.  XX and YY are purple and orange, XY and YX are black and pink (respectively).'

# poltype='Df': assume calibrator is unpolarized; solve for each channel
# preavg=1. - time interval for applying parallactic angle corrections
#polcal(vis=smallms,field=bpcal,poltype='Df',caltable=smallname+'.Df0',refant=refanten,preavg=1.0,gaintable=gaintableGKBKc)

plotcal(caltable=smallname+'.Df0',xaxis='freq',yaxis='amp',iteration='antenna',subplot=441,plotrange=[200,500,0,1],showgui=False,figfile=plotdir+'leakage_preflag.png')
print 'Check plots/leakage_preflag.png for D terms (freq-dpdt leakage between polns for each antenna).  Should be <<1.'

# flag D terms (> 0.4? or sudden spikes?)
flagdata(vis=smallname+'.Df0',mode='clip',clipminmax=[-0.1,0.4],datacolumn='CPARAM')

plotcal(caltable=smallname+'.Df0',xaxis='freq',yaxis='amp',iteration='antenna',subplot=441,plotrange=[200,500,0,1],showgui=False,figfile=plotdir+'leakage_clip.png')
plotcal(caltable=smallname+'.Df0',xaxis='freq',yaxis='amp',iteration='antenna',subplot=441,showgui=False,figfile=plotdir+'leakage_clip_zoom.png')
print 'Check plots/leakage_clip[_zoom].png for D terms after autoflag (clip values > 0.4).  Should be <<1.'

applycal(vis=smallms,field=bpcal,applymode='calflagstrict',gaintable=gaintableGKBKcDf)

print '% flagged in',smallms,'after applying poln calibration tables:'
summary_11 = field_flag_summary(smallms,flagdata)
np.save('summary11.npy',summary_11)

# inspect amplitude of cross-pols (XY,YX) for BP cal
plotms(vis=smallms,field=bpcal,xaxis='freq',yaxis='amp',ydatacolumn='corrected',coloraxis='corr',iteraxis='baseline',avgtime='1e8',avgscan=True,plotrange=[0.2,0.5,0,70],gridrows=3,gridcols=3,plotfile=plotdir+'bpcal_cross_corr.png',showgui=False,exprange='all',overwrite=True)
print 'Check plots/bpcal_cross_corr.png to see strength of cross polns after polcal.  XX and YY are purple and orange, XY and YX are black and pink (respectively).'

'''

# broke off here!

### REDO BANDPASS AND DELAY CALIBRATION ###

# delay calibration
gaincal(vis=smallms,gaintype='K',field=bpcal,caltable=smallname+'.K1',refant=refanten,minsnr=3.0,parang=True,gaintable=gaintableGBKcDf)
plotcal(caltable=smallname+'.K0',xaxis='antenna',yaxis='delay',figfile=plotdir+'delay.png',showgui=False)
print 'Check plots dir for delay.png - all delays should be < few ns.'


### GAIN CALIBRATION: J0204-1701 ###
# useful gaincal identifying info#
# field name: J0204-1701

## Initial phase-only gaincal ##

# apply all cals to gain calibrator (don't flag yet)
applycal(vis=smallms,field=gcal,applymode='calonly',gaintable=gaintableGKBKcDf)
# inspect in plotms: real(vis) has a wide spread (not positive), but >> for XX,YY than for XY,YX --> perhaps this is due
#   to having a number of roughly equally bright srcs in the field?

# let's do a phase-only gaincal (assuming pt src) then image gaincal field
gaincal(vis=smallms,caltable=smallname+'.G1',field=bpcal+','+gcal,gaintable=gaintableKBKcDf,solint='120s',gaintype='G',calmode='p',refant=refanten)

# apply this new calibration then inspect in plotms
applycal(vis=smallms,field=gcal,applymode='calflagstrict',gaintable=gaintableKBKcDfG)
# backup: before_applycal_9

# inspect in plotms
# frequency dependence of gain appears to change with time - perhaps better to break it back up into multiple spw's?
# don't do this for now but may come back to it


## Image gain calibrator field ##

# max uvwave ~ 60,000 --> 3.3'' res --> 1'' pixels
# imsize: 40'/0.25 GHz --> 160' --> 10,000 pixels (more if we want to image 2-3 x FOV)
clean(vis=smallms,imagename=gcal+'.v0',field=gcal,niter=1000,nterms=2,cell=['1.5arcsec'],imsize=[8192,8192])
# whoops I probably should have set "threshold"... oh well
# also I should have used usescratch=True to populate the MODEL_DATA column


## Use model image to populate model column and inspect ##

# use ft to populate MODEL_DATA column of ms from model images (.tt0 and .tt1)
ft(vis=smallms,nterms=2,model=[gcal+'.v0.model.tt0',gcal+'.v0.model.tt1'],usescratch=True,field=gcal)


## Amp+phase gaincal using clean model of pcal field ##

gaincal(vis=smallms,caltable=smallname+'.G2',field=gcal,gaintable=gaintableKBKcDf,solint='120s',gaintype='G',calmode='ap',refant=refanten)
# I think it will use MODEL_DATA column by default?

# plot gains and compare to previous phase-only model-less gaincal
plotcal(caltable=smallname+'.G1',xaxis='time',yaxis='phase',subplot=441,iteration='antenna',showgui=False,figfile=plotdir+'G1phase.png')
plotcal(caltable=smallname+'.G2',xaxis='time',yaxis='phase',subplot=441,iteration='antenna',showgui=False,figfile=plotdir+'G2phase.png')
plotcal(caltable=smallname+'.G2',xaxis='time',yaxis='amp',subplot=441,iteration='antenna',showgui=False,figfile=plotdir+'G2amp.png')
# phases are a little different; amps vary from ~0.8-1.2 --> good thing we did amp gaincal


### APPLY CALIBRATIONS TO TARGET ###

# we are skipping the step from the tutorial to interpolate between gain solutions because I don't think
#   there are any failed gain solutions  (failed gain solutions would result in flagged data at those times)

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
