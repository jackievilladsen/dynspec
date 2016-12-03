'''
Now I need to improve calibration of J0132 and by extension of UV Cet.  First steps:
- try concatenating the two data sets - can we force CASA to split them into scans?
- otherwise return to AIPS and try to use it to mash J0132 and UV Cet together after shifting pcen
- image J0132 and quiescent UV Cet
- try running an amplitude+phase gaincal on J0132 with no model, one soln per scan
- apply ap gains to J0132 & UV Cet and re-image, compare quality

'''

### Import J0132 and UVCet avg'd data sets ###

importuvfits(fitsfile='J0132.uvfits', vis='J0132.ms')
importuvfits(fitsfile='UVCet.uvfits', vis='UVCet.ms')


### Write my own scan numbers ###

from dynspec.tbavg import scan_reindex

scan_reindex(vis='J0132.ms',gap=50.)
scan_reindex(vis='UVCet.ms',gap=50.)


### Image J0132 and UVCet and assess sensitivity ###

# Expected sensitivity:
# time on source: J0132 - 74 scans x 50 sec = 62 min on source
#                 UVCet - 72 scans x 110 sec = 132 min on source
# EVN calculator, for 256 MHz dual pol:
#   62 min  --> 40 uJy (-1 ant), 36 uJy (all ants)
#   132 min --> 27 uJy (-1 ant), 24 uJy (all ants)
#   2 min   --> 200 uJy (all ants)
#   1 min   --> 280 uJy (all ants)

## Dirty images (to look at Stokes V RMS) ##

# beam: 2.8x1.2 mas --> 0.28 pixel size is fine
clean(vis='J0132.ms',imagename='J0132/J0132_precal_dirty',cell='0.00028arcsec',niter=0,stokes='IV')
clean(vis='UVCet.ms',imagename='UVCet/UVCet_precal_dirty',cell='0.00028arcsec',niter=0,stokes='IV')
# Stokes V RMS: J0132 - 1.3 mJy (dominated by source leakage from Stokes I)
#               UVCet - 62 uJy (40 in region w/o src) --> already pretty good!

# UV Cet looks resolved but is this due to phase calibration errors?
# The antennas where J0132 and J0140 gain phases are significantly different are: HN, MK
# Dirty image of UVCet w/o HN, MK:
clean(vis='UVCet.ms',imagename='UVCet/noHNMK_precal_dirty',cell='0.00028arcsec',niter=0,stokes='IV',antenna='!HN,MK')
# beam is 3.9x1.8 mas
# dirty image now looks less clearly resolved but if we clean and break it down by time it could still be

# create ms w/ pcen shifted to UVCet peak w/o HN, MK
# peak: 01:39:05.116, -17:56:51.78
fixvis(vis='UVCet.ms',outputvis='UVCet_pcen.ms',phasecenter='J2000 1h39m05.1157272 -17d56m51.7775896')


## Clean images (to look at Stokes I residual RMS) ##

# identify quiescent timerange for UVCet
t_quiet = '12:00:00~12:34:00,14:41:00~16:00:00'

# clean each field once interactively to generate mask
clean(vis='J0132.ms',imagename='J0132/J0132',cell='0.00028arcsec',interactive=True,niter=5000,cyclefactor=5,npercycle=10) # 90 iterations
clean(vis='UVCet.ms',imagename='UVCet/UVCet',timerange=t_quiet,antenna='!MK,HN',cell='0.00028arcsec',interactive=True,niter=5000,cyclefactor=5,npercycle=10)
pcal_mask = 'J0132/J0132.mask'
star_mask = 'UVCet/UVCet.mask'

# clean each field automatically using mask
clean(vis='J0132.ms',imagename='J0132/auto_precal',cell='0.00028arcsec',niter=90,cyclefactor=5,mask=pcal_mask,usescratch=True)
clean(vis='UVCet.ms',imagename='UVCet/auto_precal_v2',cell='0.00028arcsec',antenna='!MK,HN',timerange=t_quiet,threshold='0.15mJy',cyclefactor=5,mask=star_mask,usescratch=True)
# residual RMS: J0132 - 3.9 mJy (dynamic range: 1 Jy integrated, 700 mJy peak flux/3.9 = 180 --> UV Cet not dyn range limited if equally good)
#               UVCet - 70 uJy (--> ~2x expected) - not bad! much worse when used full time instead of quiet


### Amp+phase gaincal w/ and w/o model ###

# split data column only for gaincal w/o model
split(vis='J0132.ms',outputvis='J0132_v2.ms',datacolumn='data')

# try phase-only gaincal first
gaincal(vis='J0132.ms',refant='LA',calmode='p',caltable='J0132.ScP0') # ScP0: phase-only self-cal #0 (with model)
gaincal(vis='J0132_v2.ms',refant='LA',calmode='p',caltable='J0132.P0') # P0: phase-only gaincal #0 (no model)
# for no-model gaincal, a few solns were flagged due to low SNR; not for selfcal w/ model - all for ant 6 (based on # pts plotted by plotcal)

plotcal(caltable='J0132.ScP0', xaxis='time',yaxis='phase',iteration='antenna',subplot=431)
plotcal(caltable='J0132.P0', xaxis='time',yaxis='phase',iteration='antenna',subplot=431)
# MK6 pretty structured over time, significantly different shape w/ and w/o model

# try amp+phase gaincal w/o phase-only cal (can we get away w/ just doing one amp+phase gaincal?)
gaincal(vis='J0132.ms',refant='LA',calmode='ap',caltable='J0132.ScAP0') # ScAP0: amp+phase self-cal #0 (with model)
gaincal(vis='J0132_v2.ms',refant='LA',calmode='ap',caltable='J0132.AP0') # AP0: amp+phase gaincal #0 (no model)
# no-model gaincal reported 1-2 solns (probably 1 or 2 pols of ant 6) flagged due to low SNR for many scans

plotcal(caltable='J0132.ScAP0', xaxis='time',yaxis='phase',iteration='antenna',subplot=431)
plotcal(caltable='J0132.AP0', xaxis='time',yaxis='phase',iteration='antenna',subplot=431)
# now model-less phases agree w/ model phases, which agree w/ phase only selfcal phases - but model-less ap phases have no solns for 1st half of obs
# model-less amp gains very low, whereas model amp gains all centered around 1

# apply caltable to J0132 and image
for exten in ['ScP0','P0','ScAP0','AP0']:
   applycal(vis='J0132.ms',gaintable='J0132.'+exten,applymode='calonly')
   clean(vis='J0132.ms',imagename='J0132/dirty_'+exten,cell='0.00028arcsec',niter=0,stokes='IV')
   clean(vis='J0132.ms',imagename='J0132/clean_'+exten,cell='0.00028arcsec',mask=pcal_mask,niter=90,cyclefactor=5)
# dirty Stokes V RMS:
#  - precal, ScP0, P0: 1.28 mJy
#  - AP0: 661 uJy
#  - ScAP0: 75 uJy (2x expected --> pretty good!)
# clean Stokes I residual RMS:
#  - precal, ScP0: 3.9 mJy;  P0: 4.0 mJy; AP0: 5.2; ScAP0: 2.6 mJy (improving mask may improve rms, but there are clear cal errors)

# ScAP0 is winner on J0132 - apply ScAP0 to UVCet and image
applycal(vis='UVCet.ms',gaintable='J0132.ScAP0',applymode='calonly')
clean(vis='UVCet.ms',imagename='UVCet/v2_dirty_ScAP0',cell='0.00028arcsec',niter=0,stokes='IV')
clean(vis='UVCet.ms',imagename='UVCet/v2_clean_ScAP0_noMKHN',cell='0.00028arcsec',antenna='!MK,HN',timerange=t_quiet,threshold='0.15mJy',cyclefactor=5,mask=star_mask,usescratch=True)
clean(vis='UVCet.ms',imagename='UVCet/v2_clean_ScAP0',cell='0.00028arcsec',timerange=t_quiet,threshold='0.15mJy',cyclefactor=5,mask=star_mask,usescratch=True)
# v2 is b/c first time I imaged I forgot to run applycal
# dirty Stokes V RMS: precal: 62 uJy; ScAP0: 60 uJy
# clean Stokes I residual RMS: precal noHNMK: 69 uJy; precal all ants: 61 (image created below)
#                               ScAP0 noHNMK: 66;      ScAP0 all ants: 57
# --> amp+phase selfcal doesn't make huge difference to quality of UV Cet image - slight improvement, doesn't seem to hurt

# split off precal UVCet data and image
split(vis='UVCet.ms',outputvis='UVCet_nocal.ms',datacolumn='data')
clean(vis='UVCet_nocal.ms',imagename='UVCet/auto_precal_wHNMK',cell='0.00028arcsec',timerange=t_quiet,threshold='0.15mJy',cyclefactor=5,mask=star_mask)
# images: measured peak flux for UV Cet is ~20% lower w/ HN and MK than without, but integrated flux is higher

# based on this, I may want to use a uvrange cutoff or antenna cutoff when making time series

# peak flux location is same w/ and w/o HN/MK, and before and after amp+phase gaincal - we can use this as phase center for time series
# gaussfit params based on UVCet v2_clean_ScAP0 (w/ HN MK):
# --- ra:    01:39:05.1157390 +/- 0.0000039 s (0.0000561 arcsec along great circle)
# --- dec:  -17.56.51.7777031 +/- 0.0001366 arcsec
# --- major axis FWHM:     2.58 +/- 0.50 marcsec
# --- minor axis FWHM:     1.71 +/- 0.21 marcsec
# --- position angle: 9.3 +/- 28.1 deg
# --- Integrated:   1.87 +/- 0.18 mJy                 
# --- Peak:         873 +/- 60 uJy/beam
#
# gaussfit params based on UVCet v2_clean_ScAP0_noHNMK:
# --- ra:    01:39:05.1157421 +/- 0.0000034 s (0.0000491 arcsec along great circle)             
# --- dec: -017.56.51.7776495 +/- 0.0001383 arcsec
# --- major axis FWHM:     2.75 +/- 0.67 marcsec
# --- minor axis FWHM:     1.51 +/- 0.24 marcsec
# --- position angle: 9.3 +/- 20.2 deg
# --- Integrated:   1.78 +/- 0.14 mJy
# --- Peak:         1.150 +/- 0.058 mJy/beam
# Everything is consistent w/in <1sigma except peak flux (since that depends on clean beam size)

# now compare quiescent burst location to flare peak - burst peak is scan 
peak_scan = '37'
clean(vis='UVCet_nocal.ms',imagename='UVCet/burst_nocal_allants_dirty',cell='0.00028arcsec',niter=0,stokes='IV',scan=peak_scan)
clean(vis='UVCet.ms',imagename='UVCet/burst_APcal_allants_dirty',cell='0.00028arcsec',niter=0,stokes='IV',scan=peak_scan)
clean(vis='UVCet_nocal.ms',imagename='UVCet/burst_nocal_noHNMK_dirty',cell='0.00028arcsec',niter=0,stokes='IV',scan=peak_scan,antenna='!MK,HN')
clean(vis='UVCet.ms',imagename='UVCet/burst_APcal_noMKHN_dirty',cell='0.00028arcsec',niter=0,stokes='IV',scan=peak_scan,antenna='!MK,HN')
# MK & HN make it look like 2 lobes instead of one - I think the reason they did not do this for the quiescent emission
#  is that the quiescent emission is resolved, so the long baselines have lower amplitude and are downweighted

# conclusion: drop HN and MK for the time series
# uvrange to avoid resolved:
#   largest dimension of quiescent emission is 2.75 mas
#    this is 75 Mlambda, 2.7 Mm --> cutoff at 60 Mlambda, 2.2 Mm (this should include about half the data and will already cut out almost all MK,HN BLs)
# AD Leo I used 2e6m - we could just use this for both
# try imaging with this cutoff


### Try gaincal restricted to short baselines, see if it changes solutions ###

## iteration 0: < 2.2e6m ##
uvrange = '<2.2e6m'

gaincal(vis='J0132_v2.ms',refant='LA',calmode='p',caltable='J0132.shortP0',uvrange=uvrange) # shortP0: short-BL, phase-only gaincal #0 (no model)
plotcal(caltable='J0132.shortP0', xaxis='time',yaxis='phase',iteration='antenna',subplot=431)
# few to no solns for HN3, MK6, SC10
# phases look basically the same as J0132.ScP0 for other ants --> long BL antennas not messing up gain solns for short BL antennas

# image w/o HN3, MK6, SC10
anten='!HN,MK,SC'
applycal(vis='J0132.ms',gaintable='J0132.shortP0',applymode='calonly')
clean(vis='J0132.ms',imagename='J0132/shortP0_v2',cell='0.00028arcsec',niter=5000,interactive=True,npercycle=10,cyclefactor=5,usescratch=True,antenna=anten)
# I may have overcleaned... come back to this.  May be beneficial to apply single amplitude cal based on J0140 first to get better dynamic range

## iteration 1: <3e6m ##
uvrange='<3e6m'

# populate model column (I'm concerned that clean may have excluded antennas)
ft(vis='J0132.ms',model='J0132/shortP0_v2.model',usescratch=True)

gaincal(vis='J0132.ms',refant='LA',calmode='p',caltable='J0132.shortP1',uvrange=uvrange) # shortP1: short-BL, phase-only selfcal #1 (w/ model)
plotcal(caltable='J0132.shortP1', xaxis='time',yaxis='phase',iteration='antenna',subplot=431)
# no solutions for SC even though there's data - probably b/c of minblperant (default is 4)

# try lower values for minblperant
gaincal(vis='J0132.ms',refant='LA',calmode='p',caltable='J0132.shortP1bl',uvrange=uvrange,minblperant=1) # shortP1bl: short-BL, phase-only selfcal #1 (w/ model, minblperant=1)
plotcal(caltable='J0132.shortP1bl', xaxis='time',yaxis='phase',iteration='antenna',subplot=431)
# solns for HN3 and SC10 don't cover full time range: need to run applycal w/ flagging before imaging

# copy data set to a version I can flag
split(vis='J0132_v2.ms',outputvis='J0132_shortP1.ms',datacolumn='data')

# applycal w/ flagging
applycal(vis='J0132_shortP1.ms',gaintable='J0132.shortP1bl',applymode='calflagstrict')
clean(vis='J0132_shortP1.ms',imagename='J0132/shortP1bl',cell='0.00028arcsec',niter=5000,interactive=True,npercycle=10,cyclefactor=5,usescratch=True)

# ft to put model in full (unflagged) data set then inspect amp vs. uvdist in plotms for model and data
ft(vis='J0132.ms',model='J0132/shortP1bl.model',usescratch=True)
# data/model amp is good out to 5e6 m, phase to 3.5e6 m (but we are trying to fix data phase)
# MK is the only antenna that strongly diverges from model at this point - b/c we did not use it to make model!

## Iteration 2: <5e6m ##
uvrange='<5e6m'

gaincal(vis='J0132.ms',refant='LA',calmode='p',caltable='J0132.shortP2',uvrange=uvrange,minblperant=1) # shortP2: short-BL, phase-only selfcal #2 (w/ model)
plotcal(caltable='J0132.shortP2', xaxis='time',yaxis='phase',iteration='antenna',subplot=431)

gaincal(vis='J0132.ms',refant='LA',calmode='p',caltable='J0132.shortP2_4e6',uvrange='<4e6m',minblperant=1) # shortP2: short-BL, phase-only selfcal #2 (w/ model)
plotcal(caltable='J0132.shortP2_4e6', xaxis='time',yaxis='phase',iteration='antenna',subplot=431)
# result is quite different for 4e6 and 5e6 --> let's go with 4e6?

split(vis='J0132_v2.ms',outputvis='J0132_shortP2.ms',datacolumn='data')
applycal(vis='J0132_shortP2.ms',gaintable='J0132.shortP2_4e6',applymode='calflagstrict')
clean(vis='J0132_shortP2.ms',imagename='J0132/shortP2_4e6',cell='0.00028arcsec',niter=5000,interactive=True,npercycle=10,cyclefactor=5,usescratch=True)

# ft to put model in full (unflagged) data set then inspect amp vs. uvdist in plotms for model and data
ft(vis='J0132.ms',model='J0132/shortP2_4e6.model',usescratch=True)

## Iteration 3: full uvdist? ##

# compare <5e6m and full #

gaincal(vis='J0132.ms',refant='LA',calmode='p',caltable='J0132.shortP3_5e6',uvrange='<5e6m',minblperant=1) # phase-only #3, <5e6m
gaincal(vis='J0132.ms',refant='LA',calmode='p',caltable='J0132.shortP3') # phase-only #3, all baselines
plotcal(caltable='J0132.shortP3_5e6', xaxis='time',yaxis='phase',iteration='antenna',subplot=431)
plotcal(caltable='J0132.shortP3', xaxis='time',yaxis='phase',iteration='antenna',subplot=431)
# looks basically the same

# apply and image to supply model for amp+phase selfcal
split(vis='J0132_v2.ms',outputvis='J0132_shortP3.ms',datacolumn='data') # shouldn't flag anything but we'll use calflagstrict to check
applycal(vis='J0132_shortP3.ms',gaintable='J0132.shortP3',applymode='calflagstrict') # no additional flagging, yay!
clean(vis='J0132_shortP3.ms',imagename='J0132/shortP3',cell='0.00028arcsec',niter=5000,interactive=True,npercycle=10,cyclefactor=5,usescratch=True)

# ft and applycal and inspect
ft(vis='J0132.ms',model='J0132/shortP3.model',usescratch=True)
applycal(vis='J0132.ms',gaintable='J0132.shortP3')

### Amplitude+phase calibration ###

# do normalized amplitude calibration (so we don't reduce source flux)
gaincal(vis='J0132.ms',refant='LA',calmode='ap',caltable='J0132.shortAP0_v2',solnorm=True) # shortAP0: all baselines
gaincal(vis='J0132.ms',refant='LA',calmode='ap',caltable='J0132.shortAP0_22e6_v2',uvrange='<2.2e6m',solnorm=True) # shortAP0_22e6: start with baselines < 2.2e6m
#v2 - for first version I forgot to use solnorm=True

plotcal(caltable='J0132.shortAP0_v2', xaxis='time',yaxis='amp',iteration='antenna',subplot=431)
plotcal(caltable='J0132.shortAP0_22e6_v2', xaxis='time',yaxis='amp',iteration='antenna',subplot=431)
# including long baselines does affect relative amps of short BL ants --> probably better to iterate but I'm tired of this!

# apply and image
split(vis='J0132_v2.ms',outputvis='J0132_shortAP0.ms',datacolumn='data') # shouldn't flag anything but we'll use calflagstrict to check
applycal(vis='J0132_shortAP0.ms',gaintable='J0132.shortAP0_v2',applymode='calflagstrict') # only 0.7% flagged
clean(vis='J0132_shortAP0.ms',imagename='J0132/shortAP0',cell='0.00028arcsec',niter=5000,interactive=True,npercycle=10,cyclefactor=5,usescratch=True)
# Stokes I residual RMS: 3.9 mJy (vs 4.0 mJy for shortP3.residual)
# --> our amplitude calibration is not improving things much

# other thing to try: clean only center component then do ap cal - maybe we are dynamic range limited because im introducing fake structure?


### Compare UVCet w/ phase-only vs. amp+phase cal ###

split(vis='UVCet.ms',outputvis='UVCet_shortP3.ms',datacolumn='data')
split(vis='UVCet.ms',outputvis='UVCet_shortAP0.ms',datacolumn='data')
applycal(vis='UVCet_shortP3.ms',gaintable='J0132.shortP3',applymode='calflagstrict')
applycal(vis='UVCet_shortAP0.ms',gaintable='J0132.shortAP0_v2',applymode='calflagstrict')

peak_scan='37'

# dirty images of burst peak
# dirty image w/ no cal: UVCet/burst_nocal_allants_dirty
# dirty image w/ AP selfcal using all baselines: UVCet/burst_APcal_allants_dirty
clean(vis='UVCet_shortP3.ms',imagename='UVCet/shortP3_burst_dirty',cell='0.00028arcsec',niter=0,stokes='IV',scan=peak_scan)
clean(vis='UVCet_shortAP0.ms',imagename='UVCet/shortAP0_burst_dirty',cell='0.00028arcsec',niter=0,stokes='IV',scan=peak_scan)
# burst peak location much closer to quiescent emission peak --> believable; but it looks like the short BLs are still slightly offset from the long BLs (unless there is resolved quiescent emission with slight offset)


### Try selfcal on J0132 w/o southwest lobe ###

clean(vis='J0132_shortP3.ms',imagename='J0132/shortP3_nolobe',cell='0.00028arcsec',niter=5000,interactive=True,npercycle=10,cyclefactor=5,usescratch=True)
gaincal(vis='J0132_shortP3.ms',refant='LA',calmode='ap',caltable='J0132.shortAP0_nolobe',solnorm=True) # shortAP0: all baselines

plotcal(caltable='J0132.shortAP0_nolobe', xaxis='time',yaxis='amp',iteration='antenna',subplot=431)

split(vis='UVCet.ms',outputvis='UVCet_shortAP0_nolobe.ms',datacolumn='data')
applycal(vis='UVCet_shortAP0_nolobe.ms',gaintable='J0132.shortAP0_nolobe',applymode='calflagstrict')
clean(vis='UVCet_shortAP0_nolobe.ms',imagename='UVCet/shortAP0_nolobe_burst_dirty',cell='0.00028arcsec',niter=0,stokes='IV',scan=peak_scan)
# absolutely identical to shortAP0 w/ lobe

# clean image of burst time
clean(vis='UVCet_shortAP0.ms',imagename='UVCet/shortAP0_burst_clean',cell='0.00028arcsec',niter=5000,interactive=True,npercycle=10,cyclefactor=5,stokes='IV',scan=peak_scan)

# clean image of quiescent emission
t_quiet = '12:00:00~12:34:00,14:41:00~16:00:00'
clean(vis='UVCet_shortAP0.ms',imagename='UVCet/shortAP0_quiet_clean',timerange=t_quiet,cell='0.00028arcsec',interactive=True,niter=5000,cyclefactor=5,npercycle=10,stokes='IV')

# shift phase center to gaussfit location of burst
# --- ra:    01:39:05.11573485 +/- 0.00000079 s (0.00001122 arcsec along great circle)          
# --- dec: -017.56.51.77789101 +/- 0.00005123 arcsec
fixvis(vis='UVCet_shortAP0.ms',outputvis='UVCet_shortAP0.ms',phasecenter='J2000 01h39m05.11573485 -17d56m51.77789101')

mask = 'UVCet/shortAP0_quiet_clean.mask'
threshold='0.2mJy' # I'm guessing at what will be good
# break quiet times down into ~30 min intervals
# gave blocks of 8 scans letters (a: 1-8, b:9-16, etc)
# a-b,g-i: quiescent
# c: 17~24 - RR starting to increase
# d: 25~32 - RR almost flat at ~3 mJy
# e: 33~40 - burst peak, longer total duration b/c of bright cal obs after scan 36
# f: 41~48 - RR mostly low again
for nf in range(8,73,8):
    n0 = nf-7
    scan = str(n0)+'~'+str(nf)
    name = 'scan'+str(n0)+'_'+str(nf)
    clean(vis='UVCet_shortAP0.ms',imagename='UVCet/'+name,cell='0.00028arcsec',threshold=threshold,cyclefactor=5,mask=mask,usescratch=True,stokes='IV',scan=scan)
# uv cet appears to move by 1.1 mas during obs... check proper motion+parallax

# uv cet has huge proper motion

# load data set w/ proper motion + parallax corrected
importuvfits(fitsfile='../UVCet_shift.uvfits',vis='UVCet_shift.ms')
from dynspec.tbavg import scan_reindex
scan_reindex('UVCet_shift.ms')
applycal(vis='UVCet_shift.ms',gaintable='J0132.shortAP0_v2',applymode='calflagstrict')

for nf in range(8,73,8):
    n0 = nf-7
    scan = str(n0)+'~'+str(nf)
    name = 'scan'+str(n0)+'_'+str(nf)
    clean(vis='UVCet_shift.ms',imagename='UVCet/shift_dirty_'+name,cell='0.00028arcsec',niter=0,stokes='IV',scan=scan)

# create uncalibrated UVCet_shift to compare how src moves
split(vis='UVCet_shift.ms',outputvis='UVCet_shift_nocal.ms',datacolumn='data')
for nf in range(8,73,8):
    n0 = nf-7
    scan = str(n0)+'~'+str(nf)
    name = 'scan'+str(n0)+'_'+str(nf)
    clean(vis='UVCet_shift_nocal.ms',imagename='UVCet/shift_nocal_dirty_'+name,cell='0.00028arcsec',niter=0,stokes='IV',scan=scan)
# conclusion: motion of source during obs is not due to J0132.shortAP0 (however it could be due to problem with fringe fit that
#    J0132.shortAP0 has failed to correct)

# clean once interactively to create mask
clean(vis='UVCet_shift.ms',imagename='UVCet/shift_quiet',timerange=t_quiet,interactive=True,npercycle=10,cyclefactor=5,stokes='IV',cell='0.00028arcsec')
# wow!  Stokes V is positive on left side, negative on right side; both are offset slightly southeast of stokes I peak (which is direction towards burst)

clean(vis='UVCet_shift_nocal.ms',imagename='UVCet/shift_quiet_nocal',timerange=t_quiet,interactive=True,npercycle=10,cyclefactor=5,stokes='IV',cell='0.00028arcsec')
# Still there w/o cal - not some weird consequence of doing an amplitude cal

# now let's see if it's there over shorter time intervals or if it's a consequence of combining different times
mask = 'UVCet/shift_quiet.mask'
threshold='0.25mJy'
for nf in range(8,73,8):
    n0 = nf-7
    scan = str(n0)+'~'+str(nf)
    name = 'scan'+str(n0)+'_'+str(nf)
    clean(vis='UVCet_shift.ms',imagename='UVCet/v2_shift_clean_'+name,cell='0.00028arcsec',threshold=threshold,cyclefactor=5,mask=mask,stokes='IV',scan=scan)
# first 5 blocks have thing to south south east of UV Cet, in first 4 blocks peak appears to be moving in straight line - however it is at the level
#   of other artifacts/noise in the background so it could be a sidelobe that moves with time; perhaps due to a bad antenna that goes away after scan 5?

# there's a lot of mucky stuff - bad data? or bad cal?

# try a phase-only selfcal (no model) on UV Cet RR during burst peak
scan = '32~39'
gaincal(vis='UVCet_shift.ms',refant='LA',calmode='p',gaintable='J0132.shortAP0_v2',caltable='UVCet.P0',scan=scan)
plotcal(caltable='UVCet.P0', xaxis='time',yaxis='phase',iteration='antenna',subplot=431)
split(vis='UVCet_shift.ms',outputvis='UVCet_burstRR.ms',datacolumn='corrected',correlation='RR',scan='32~39')
gaincal(vis='UVCet_burstRR.ms',refant='LA',calmode='p',caltable='UVCet.P0_v2')
plotcal(caltable='UVCet.P0_v2', xaxis='time',yaxis='phase',iteration='antenna',subplot=431)
# need to shift phase center to star first

# Gaussfit on scan 33-40:
# --- ra:    01:39:05.11534620 +/- 0.00000100 s (0.00001431 arcsec along great circle)          
# --- dec: -017.56.51.77725822 +/- 0.00007058 arcsec
fixvis(vis='UVCet_burstRR.ms',outputvis='UVCet_burstRR.ms',phasecenter='J2000 1h39m5.11534620 -17d56m51.77725822')
gaincal(vis='UVCet_burstRR.ms',refant='LA',calmode='p',caltable='UVCet.P0_v3')
plotcal(caltable='UVCet.P0_v3', xaxis='time',yaxis='phase',iteration='antenna',subplot=431)
gaincal(vis='UVCet_burstRR.ms',refant='LA',calmode='p',caltable='UVCet.P0_v4',combine='scan')
plotcal(caltable='UVCet.P0_v4', xaxis='time',yaxis='phase',iteration='antenna',subplot=431)
gaincal(vis='UVCet_burstRR.ms',refant='LA',calmode='p',caltable='UVCet.P0_v5',combine='spw')
plotcal(caltable='UVCet.P0_v5', xaxis='time',yaxis='phase',iteration='antenna',subplot=431)
applycal(vis='UVCet_burstRR.ms',gaintable='UVCet.P0_v5',applymode='calonly',spwmap=[0,0,0,0,0,0,0,0])
clean(vis='UVCet_burstRR.ms',imagename='UVCet/burstRR',cell='0.00028arcsec',interactive=True,npercycle=10,cyclefactor=5,stokes='IV')
# the thing in the SSW goes away

# try doing selfcal on UV Cet throughout by combining spws
fixvis(vis='UVCet_shift.ms',outputvis='UVCet_shift.ms',phasecenter='J2000 1h39m5.11534620 -17d56m51.77725822')
split(vis='UVCet_shift.ms',outputvis='UVCet_shift_cal.ms',datacolumn='corrected')
gaincal(vis='UVCet_shift_cal.ms',refant='LA',calmode='p',caltable='UVCet.P0_v6',combine='spw')
plotcal(caltable='UVCet.P0_v6', xaxis='time',yaxis='phase',iteration='antenna',subplot=431)
gaincal(vis='UVCet_shift_cal.ms',refant='LA',calmode='p',caltable='UVCet.P0_v7',combine='spw,scan',solint='1440s')
plotcal(caltable='UVCet.P0_v7', xaxis='time',yaxis='phase',iteration='antenna',subplot=431)
smoothcal(vis='UVCet_shift.ms',tablein='UVCet.P0_v7',caltable='UVCet.P0_v7smooth',smoothtime=2880)
plotcal(caltable='UVCet.P0_v7smooth', xaxis='time',yaxis='phase',iteration='antenna',subplot=431)
