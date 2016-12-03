# fits files exported from AIPS on 11/19/16:
# - J0140.uvfits: full-resolution (2 sec, 0.5 MHz), dual pol
# - J0140_small.uvfits: avg'd (30 sec, 8 MHz), dual pol
# ditto for J0132 and UVCet - UVCet shifted to pcen before avg'ing, exporting

### Import J0140 data, both res ###

importuvfits(fitsfile='J0140.uvfits', vis='J0140.ms')
importuvfits(fitsfile='J0140_small.uvfits', vis='J0140_small.ms')


### Compare whether averaging affects sensitivity ###

# using J0140 for this because it's a small data set
# 2 measures of image sensitivity: dirty Stokes V RMS, clean Stokes I residuals RMS

# generate dirty images (to look at dirty Stokes V RMS)
clean(vis='J0140.ms',imagename='J0140_dirty',cell='0.00028arcsec',niter=0,stokes='IV')
clean(vis='J0140_small.ms',imagename='J0140_small_dirty',cell='0.00028arcsec',niter=0,stokes='IV')

# clean once w/ interactive to generate mask
clean(vis='J0140.ms',imagename='J0140',cell='0.00028arcsec',interactive=True,niter=5000,cyclefactor=5,npercycle=10)
# used 60 iterations to get to point where source flux low relative to sidelobes

# clean both data sets non-interactive w/ mask, 60 iterations
clean(vis='J0140_small.ms',imagename='J0140_small_auto',cell='0.00028arcsec',mask='J0140.mask',niter=60,cyclefactor=5,usescratch=True)
clean(vis='J0140.ms',imagename='J0140_auto',cell='0.00028arcsec',mask='J0140.mask',niter=60,cyclefactor=5)

# Results: (same for both --> ok to use small data set)
#   Peak flux:              320 mJy
#   Clean I residual RMS:   2.6 mJy --> dynamic range 120:1
#   Dirty V RMS:            0.91 mJy
#   Expected thermal noise: 0.18 mJy  (EVN calculator: 3 min, 256 MHz, dual pol, 1 ant missing)
# UV Cet peak flux is 7 mJy --> 7/120 = 0.06 mJy --> dynamic range limit may not be the issue; currently
#   the error on our AD Leo time series is 0.1 mJy in 30 seconds - this doesn't seem so bad - but that is
#   based on std(imag(I)), but the real(I) tseries seems to vary more than that? maybe it is real?


### Compare reducing soln interval vs. using model image ###

# reducing the solution interval should improve sensitivity on J0140 but may not help that much with 
# quality of phase referencing; using model image might help more with quality of phase referencing
# but I'm not sure if it will improve the sensitivity on J0140 (especially if I overcleaned, for example)

# Create 2 copies of J0140_small.ms - one w/o model image
split(vis='J0140_small.ms',outputvis='J0140_small2.ms',datacolumn='data')

## Approach 1: phase-only selfcal w/ model image, 1 soln per scan ##

gaincal(vis='J0140_small.ms',refant='LA',calmode='p',caltable='J0140.ScP0') # ScP0: phase-only self-cal #0
plotcal(caltable='J0140.ScP0', xaxis='time',yaxis='phase',iteration='antenna',subplot=431)
# few deg variation between scans, few tenths of deg between pols and spws

applycal(vis='J0140_small.ms',gaintable='J0140.ScP0',applymode='calonly')

# dirty image to look at Stokes V (.v2 b/c first version I accidentally let applycal flag data)
clean(vis='J0140_small.ms',imagename='J0140_dirty_ScP0.v2',cell='0.00028arcsec',niter=0,stokes='IV')
# Stokes V RMS is basically the same (901.8 uJy instead of 902.6)

# clean image to look at Stokes I residuals - try niter=60
clean(vis='J0140_small.ms',imagename='J0140_auto_ScP0.v2',cell='0.00028arcsec',mask='J0140.mask',niter=60,cyclefactor=5)
# residuals have decreased a few percent: to 2.50 mJy from 2.58

# try nterms=2
clean(vis='J0140_small.ms',imagename='J0140_nterms2_ScP0',cell='0.00028arcsec',mask='J0140.mask',niter=60,cyclefactor=5,nterms=2)
# basically identical --> nterms=2 not needed for this source


# try imaging single spw's
for i in range(8):
   spw = str(i)
   clean(vis='J0140_small.ms',imagename='J0140_spw'+spw+'_ScP0',cell='0.00028arcsec',mask='J0140.mask',niter=60,cyclefactor=5,spw=spw)
# RMS: 2.46e-3 --> using more spw's does not improve our sensitivity...
# lowest RMS: spw 3, 2.15 mJy; highest: spw 7, 3.9 mJy
# maybe an amplitude self-cal would help with this?  maybe some spw's have weird shapes? or maybe phase wrap is going
#  one direction in low spw's and one direction in high spw's and happens to be nearly flat in the middle?
# I used plotms to plot real(vis) avg'd over baselines vs. spw - min is at spw=3

# plotted amp vs. phase for J0140 - one baseline has lowest amps, non-zero phase, intermediate uvdist:
# BR & SC (1&10)

# clean image again with nterms=1, all spw's, to populate model column
clean(vis='J0140_small.ms',imagename='J0140_auto_ScP0.v3',cell='0.00028arcsec',mask='J0140.mask',niter=60,cyclefactor=5,usescratch=True)
# this is the only baseline where phase strongly departs from the model - let's flag it?

# merge into 1 spw since plotms can't avg over spw's for data imported from aips
mstransform(vis='J0140_small.ms',outputvis='J0140_1spw.ms',combinespws=True,datacolumn='all')

# look at phase vs. amp for J0132 - same problem? no. J0132 strongly resolved, a lot of baselines look similar to
#   BR&SC --> maybe BR&SC is real structure, I'm just not imaging it correctly?

# try flagging BR&SC then imaging
flagdata(vis='J0140_1spw.ms',antenna='BR&SC')  # flag backup is flagdata_1
clean(vis='J0140_1spw.ms',imagename='J0140_noBRSC_ScP0',cell='0.00028arcsec',mask='J0140.mask',niter=60,cyclefactor=5,usescratch=True) # rms slightly worse (2.52 mJy instead of 2.50)
flagmanager(vis='J0140_1spw.ms',mode='restore',versionname='flagdata_1') # restore flags since it didn't help

# try imaging without mask
clean(vis='J0140_small.ms',imagename='J0140_nomask_ScP0',cell='0.00028arcsec',niter=60,cyclefactor=5,usescratch=True)
# RMS is 2.3 and negative sidelobes are less deep --> my clean mask is probably missing real flux


## APPROACH 2: phase-only gaincal w/ no model, multiple solns per scan ##

gaincal(vis='J0140_small2.ms',refant='LA',calmode='p',caltable='J0140.p90s',solint='inf')
gaincal(vis='J0140_small2.ms',refant='LA',calmode='p',caltable='J0140.p30s',solint='30s') # p30s: phase-only gaincal 30s
gaincal(vis='J0140_small2.ms',refant='LA',calmode='p',caltable='J0140.p20s',solint='20s')
gaincal(vis='J0140_small2.ms',refant='LA',calmode='p',caltable='J0140.p10s',solint='10s')
gaincal(vis='J0140_small2.ms',refant='LA',calmode='p',caltable='J0140.p2s',solint='2s')


plotcal(caltable='J0140.p90s', xaxis='time',yaxis='phase',iteration='antenna',subplot=431)
plotcal(caltable='J0140.p30s', xaxis='time',yaxis='phase',iteration='antenna',subplot=431)
plotcal(caltable='J0140.p20s', xaxis='time',yaxis='phase',iteration='antenna',subplot=431)
plotcal(caltable='J0140.p10s', xaxis='time',yaxis='phase',iteration='antenna',subplot=431)
plotcal(caltable='J0140.p2s', xaxis='time',yaxis='phase',iteration='antenna',subplot=431)
# time variation of up to 20 deg w/in scans, few deg between spws and pols
# this makes me wish I had requested another phase center on a bright target w/in the FOV so I could selfcal w/in FOV
# could probably combine spws to improve SNR


for exten in ['p90s','p30s','p20s','p10s','p2s']:
   applycal(vis='J0140_small2.ms',gaintable='J0140.'+exten,applymode='calonly')
   clean(vis='J0140_small2.ms',imagename='dirty_'+exten,cell='0.00028arcsec',niter=0,stokes='IV')
   clean(vis='J0140_small2.ms',imagename='clean_'+exten,cell='0.00028arcsec',mask='J0140.mask',niter=60,cyclefactor=5)
# dirty V RMS: .908-.911 --> a bit more than w/o calibration (.904), but comparable
# clean I residual RMS: 90s-2.5, 30s-2.45, 20s-2.44, 10s,2s-2.43

# Again, we get only marginal improvement from this


### TRY ADDING AMPLITUDE CALIBRATION ###

## no model ##
gaincal(vis='J0140_small2.ms',refant='LA',calmode='ap',caltable='J0140.ap90s',solint='inf')
gaincal(vis='J0140_small2.ms',refant='LA',calmode='ap',caltable='J0140.ap30s',solint='30s') # ap30s: complex gaincal 30s
gaincal(vis='J0140_small2.ms',refant='LA',calmode='ap',caltable='J0140.ap20s',solint='20s')
gaincal(vis='J0140_small2.ms',refant='LA',calmode='ap',caltable='J0140.ap10s',solint='10s')
gaincal(vis='J0140_small2.ms',refant='LA',calmode='ap',caltable='J0140.ap2s',solint='2s')

plotcal(caltable='J0140.ap90s', xaxis='time',yaxis='amp',iteration='antenna',subplot=431)
plotcal(caltable='J0140.ap30s', xaxis='time',yaxis='amp',iteration='antenna',subplot=431)
plotcal(caltable='J0140.ap20s', xaxis='time',yaxis='amp',iteration='antenna',subplot=431)
plotcal(caltable='J0140.ap10s', xaxis='time',yaxis='amp',iteration='antenna',subplot=431)
plotcal(caltable='J0140.ap2s', xaxis='time',yaxis='amp',iteration='antenna',subplot=431)
# varies up to 15% between antennas, spws, scans for solint='inf'; upt to 35% over time with shorter solint

for exten in ['ap90s','ap30s','ap20s','ap10s','ap2s']:
   applycal(vis='J0140_small2.ms',gaintable='J0140.'+exten,applymode='calonly')
   clean(vis='J0140_small2.ms',imagename='dirty_'+exten,cell='0.00028arcsec',niter=0,stokes='IV')
   clean(vis='J0140_small2.ms',imagename='clean_'+exten,cell='0.00028arcsec',mask='J0140.mask',niter=60,cyclefactor=5)
# dirty V RMS: 90s-210 uJy, 30s to 2s-185 uJy --> we get improvement from 30s instead of scan but not beyond that
#    185 uJy is pretty much expected thermal noise, woohoo!
#    This means that Stokes V was also dynamic range limited due to differences between amplitudes of RR and LL which
#     came from using only an a priori amplitude calibration
# clean I residual RMS: 550-560 uJy, but clearly need to increase niter now

# interactive clean to decide how many iterations
clean(vis='J0140_small2.ms',imagename='interactive_'+exten,cell='0.00028arcsec',mask='J0140.mask',niter=600,cyclefactor=5,interactive=True,npercycle=10)
# niter~75 brings it down to 450 uJy
# could use threshold now instead: threshold=200 uJy?

clean(vis='J0140_small2.ms',imagename='thresh600_'+exten,cell='0.00028arcsec',mask='J0140.mask',niter=600,cyclefactor=5,threshold='0.6mJy')
clean(vis='J0140_small2.ms',imagename='thresh300_'+exten,cell='0.00028arcsec',mask='J0140.mask',niter=600,cyclefactor=5,threshold='0.3mJy')
clean(vis='J0140_small2.ms',imagename='thresh200_'+exten,cell='0.00028arcsec',mask='J0140.mask',niter=600,cyclefactor=5,threshold='0.2mJy')
# 600: 394 uJy RMS, 300: 335, 200: 332 (but maxed out at 600 iterations)
clean(vis='J0140_small2.ms',imagename='thresh200_'+exten+'.v3',cell='0.00028arcsec',mask='J0140.mask',niter=3000,cyclefactor=5,threshold='0.2mJy',usescratch=True)
# now reaches 319 uJy RMS --> pretty good!

# Now use our optimized clean parameters to compare different ap gain solution intervals
for exten in ['ap90s','ap30s','ap20s','ap10s','ap2s']:
   applycal(vis='J0140_small2.ms',gaintable='J0140.'+exten,applymode='calonly')
   clean(vis='J0140_small2.ms',imagename='v3clean_'+exten,cell='0.00028arcsec',mask='J0140.mask',cyclefactor=5,niter=3000,threshold='0.3mJy')
# 330, 430, 330, 410, 320 - the alternation is probably due to wrong clean threshold...
# nope same thing happens with thresh of 200 or 300... weird.

# anyways I am happy to go ahead with 1 complex gain solution per scan, it seems to fix most problems


### TRY AMP+PHASE SELFCAL NOW THAT WE HAVE NICE MODEL IMAGE ###

# apply per-scan ap gaincal, split corrected data, then image to populate model column
applycal(vis='J0140_small2.ms',gaintable='J0140.ap90s',applymode='calonly')
split(vis='J0140_small2.ms',outputvis='J0140_ap.ms')
clean(vis='J0140_ap.ms',imagename='J0140_ap_preselfcal',cell='0.00028arcsec',mask='J0140.mask',cyclefactor=5,niter=3000,threshold='0.3mJy',usescratch=True)

# ap gaincal
gaincal(vis='J0140_ap.ms',refant='LA',calmode='ap',caltable='J0140_ap.ap90s',solint='inf')
gaincal(vis='J0140_ap.ms',refant='LA',calmode='ap',caltable='J0140_ap.ap20s',solint='20s')
gaincal(vis='J0140_ap.ms',refant='LA',calmode='ap',caltable='J0140_ap.ap2s',solint='2s')

plotcal(caltable='J0140_ap.ap90s', xaxis='time',yaxis='amp',iteration='antenna',subplot=431)
plotcal(caltable='J0140_ap.ap20s', xaxis='time',yaxis='amp',iteration='antenna',subplot=431)
plotcal(caltable='J0140_ap.ap2s', xaxis='time',yaxis='amp',iteration='antenna',subplot=431)
# 90s: phases are all pretty much zero, amps are few percent for FD and SC, few tenths of percent for others

for exten in ['ap90s','ap20s','ap2s']:
   applycal(vis='J0140_ap.ms',gaintable='J0140_ap.'+exten,applymode='calonly')
   clean(vis='J0140_ap.ms',imagename='J0140_ap_'+exten,cell='0.00028arcsec',mask='J0140.mask',cyclefactor=5,niter=3000,threshold='0.3mJy')
# 325, 334, 315 (compare to 345 before) --> selfcal w/ 90s yields almost 10% improvement, could be helpful but not
#    as critical as doing an amplitude cal (however may become more important for J0132 since it's strongly resolved)


### Path ahead ###

'''
Now I need to improve calibration of J0132 and by extension of UV Cet.  First steps:
- try concatenating the two data sets - can we force CASA to split them into scans?
- otherwise return to AIPS and try to use it to mash J0132 and UV Cet together after shifting pcen
- image J0132 and quiescent UV Cet
- try running an amplitude+phase gaincal on J0132 with no model, one soln per scan
- apply ap gains to J0132 & UV Cet and re-image, compare quality

'''
