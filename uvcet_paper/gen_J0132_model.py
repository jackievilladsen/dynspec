'''
gen_J0132_model.py

Strategy: combine data sets of J0140 and adjacent scans of J0132 which are phase referenced to J0140.
Do additional amp+phase gaincal based on J0140 then use all 6 J0132 scans to image J0132.
Combine all J0132 observations (calibrated using self). Assess relative uv coverage of J0140-adjacent J0132 scans
and all J0132 data.  Use J0132 model image to do phase-only gaincal
(or amp+phase? maybe apply amp gains from J0140). Image J0132 using full data set and use this as model for
final amp+phase gaincal.

'''

from dynspec.tbavg import scan_reindex

obsdir = '/data/jrv/BV071/UVCet/'
refms = 'J0132_J0140.ms'
ap0 = obsdir+'J0140.AP0_v2'
refms_nocal = 'J0132_J0140.nocal.ms'

### Load and combine J0132 and J0140 data ###

mslist = []
obslist = ['3','4','5']
# scans 37, 74 of J0132_ref.uvfits are adjacent to J0140 scans
for obs in obslist:
    J0132_uvfits = obsdir+obs+'/J0132_ref.uvfits'
    J0140_uvfits = '/data/jrv/15A-416/UVCet/'+obs+'/X/J0140.uvfits'
    J0132ms = obsdir + obs+'/J0132_ref.ms'
    J0140ms = obsdir + obs+'/J0140.ms'
    importuvfits(fitsfile=J0132_uvfits,vis=J0132ms)
    importuvfits(fitsfile=J0140_uvfits,vis=J0140ms)
    scan_reindex(J0132ms)
    scan_reindex(J0140ms)
    J0132_split_ms = obsdir+obs+'/J0132_split.ms'
    split(vis=J0132ms,outputvis=J0132_split_ms,datacolumn='data',scan='37,74')
    mslist.append(J0132_split_ms)
    mslist.append(J0140ms)

mslist = ['/data/jrv/BV071/UVCet/3/J0132_split.ms',
 '/data/jrv/BV071/UVCet/3/J0140.ms',
 '/data/jrv/BV071/UVCet/4/J0132_split.ms',
 '/data/jrv/BV071/UVCet/4/J0140.ms',
 '/data/jrv/BV071/UVCet/5/J0132_split.ms',
 '/data/jrv/BV071/UVCet/5/J0140.ms']

# combine all J0140 and adjacent J0132 scans into one ms and renumber scans
tmpms = 'temp.ms'
concat(vis=mslist,concatvis=tmpms,timesort=True)
refms = 'J0132_J0140.ms'
split(tmpms,outputvis=refms,datacolumn='data')
# for some reason this split prevents a segmentation fault when I run applycal but then scan_reindex doesn't work right
scan_reindex(refms,gap=5.)

# flag first integration in each scan
flagdata(vis=refms,mode='quack',quackinterval=1.0)

# inspect amp vs phase for both calibrators in plotms
# antenna 5 has low amps for both calibrators for scan 7 - may need to flag it?
# for now just keep an eye out to see if it looks bad

### Selfcal on J0140 and apply to nearest J0132 scans ###

# inspected in plotms: looks like antenna 5 is bad for scan 6
flagdata(vis=refms,scan='6',antenna='LA')

# make copy of refms w/o calibration
refms_nocal = 'J0132_J0140.nocal.ms'
split(vis=refms,outputvis=refms_nocal,datacolumn='data')

# using FD2 as refant b/c of low LA5 gains for epoch 4

#p0 = obsdir + 'J0140.P0'
#gaincal(vis=refms,refant='FD',calmode='p',field='J0140-15',caltable=p0,minblperant=1,antenna='!BR&SC')
ap0 = obsdir + 'J0140.AP0_v2'
gaincal(vis=refms,refant='FD',calmode='ap',field='J0140-15',caltable=ap0,minblperant=1,antenna='!BR&SC',solnorm=True)
plotcal(caltable=ap0, xaxis='scan',yaxis='phase',iteration='antenna',subplot=431)
plotcal(caltable=ap0, xaxis='scan',yaxis='amp',iteration='antenna',subplot=431)

# apply calibration to J0132 and J0140
applycal(vis=refms,gaintable=ap0,interp='nearest',flagbackup=False)

# dirty image J0132
clean(vis=refms,field='J0132-16',imagename='images/J0132_dirty',niter=0,stokes='IV',cell='.00028arcsec')
clean(vis=refms_nocal,field='J0132-16',imagename='images/J0132_dirty_nocal',niter=0,stokes='IV',cell='.00028arcsec')

# clean J0132
clean(vis=refms_nocal,field='J0132-16',imagename='images/J0132_clean_nocal',niter=5000,stokes='IV',cell='.00028arcsec',npercycle=10,cyclefactor=5,interactive=True)
clean(vis=refms,field='J0132-16',imagename='images/J0132_clean',niter=5000,stokes='IV',cell='.00028arcsec',npercycle=10,cyclefactor=5,interactive=True,mask='images/J0132_clean_nocal.mask')

# clean automatically
mask = 'images/J0132_clean_nocal.mask'
clean(vis=refms_nocal,field='J0132-16',imagename='images/J0132_clean_nocal_auto',threshold='60mJy',cell='.00028arcsec',cyclefactor=5,mask=mask,usescratch=True)
clean(vis=refms,field='J0132-16',imagename='images/J0132_clean_auto',threshold='60mJy',cell='.00028arcsec',cyclefactor=5,mask=mask,usescratch=True)

# solve for phase-only gains for J0132 based on model
p1 = obsdir+'J0132.P1_v3'
gaincal(vis=refms,field='J0132-16',calmode='p',refant='FD',caltable=p1,gaintable=ap0,interp='nearest')
gaincal(vis=refms_nocal,field='J0132-16',calmode='p',refant='FD',caltable=p1_nocal)
plotcal(caltable=p1, xaxis='scan',yaxis='phase',iteration='antenna',subplot=431)

# apply phases
applycal(vis=refms,field='J0132-16',gaintable=[p1,ap0],interp='nearest',flagbackup=False)
refms2 = obsdir+'P1_v2.ms'
split(vis=refms,outputvis=refms2,datacolumn='corrected')

# image
clean(vis=refms2,field='J0132-16',imagename='images/P1_dirty_v2',niter=0,stokes='IV',cell='.00028arcsec')
mask='images/P1_clean.mask'
clean(vis=refms2,field='J0132-16',imagename='images/P1_clean_v2',threshold='60mJy',stokes='IV',cell='.00028arcsec',cyclefactor=5,mask=mask,usescratch=True)

p2 = obsdir+'J0132.P2_v2'
gaincal(vis=refms2,field='J0132-16',calmode='p',refant='FD',caltable=p2)
plotcal(caltable=p2, xaxis='scan',yaxis='phase',iteration='antenna',subplot=431)
applycal(vis=refms2,field='J0132-16',gaintable=p2,interp='nearest',flagbackup=False)

p2ms = obsdir + 'P2.ms'
split(vis=refms2,outputvis=p2ms,datacolumn='corrected')
clean(vis=p2ms,field='J0132-16',imagename='images/P2_dirty',niter=0,stokes='IV',cell='.00028arcsec')
clean(vis=p2ms,field='J0132-16',imagename='images/P2_clean_v2',threshold='60mJy',cell='.00028arcsec',cyclefactor=5,stokes='IV',mask=mask,usescratch=True)
# looks way better than P1

clean(vis=p2ms,field='J0132-16',imagename='images/P2_clean_manual',interactive=True,cell='.00028arcsec',cyclefactor=5,stokes='I',usescratch=True)
p3 = obsdir+'J0132.P3'
gaincal(vis=p2ms,field='J0132-16',calmode='p',refant='FD',caltable=p3)
plotcal(caltable=p3, xaxis='scan',yaxis='phase',iteration='antenna',subplot=431)
applycal(vis=p2ms,field='J0132-16',gaintable=p3,interp='nearest',flagbackup=False)

p3ms = obsdir+'P3.ms'
split(vis=p2ms,outputvis=p3ms,datacolumn='corrected')
clean(vis=p3ms,field='J0132-16',imagename='images/P3_dirty',niter=0,stokes='IV',cell='.00028arcsec')
clean(vis=p3ms,field='J0132-16',imagename='images/P3_clean_manual',interactive=True,npercycle=10,cell='.00028arcsec',cyclefactor=5,stokes='I',usescratch=True)
clean(vis=p3ms,field='J0132-16',imagename='images/P3_clean_manual_w_lobe',interactive=True,npercycle=10,cell='.00028arcsec',cyclefactor=5,stokes='I',usescratch=True)

p4 = obsdir+'J0132.P4'
gaincal(vis=p3ms,field='J0132-16',calmode='p',refant='FD',caltable=p4)
plotcal(caltable=p4, xaxis='scan',yaxis='phase',iteration='antenna',subplot=431)
applycal(vis=p3ms,field='J0132-16',gaintable=p4,interp='nearest',flagbackup=False)

p4ms = obsdir+'P4.ms'
split(vis=p3ms,outputvis=p4ms,datacolumn='corrected')
clean(vis=p4ms,field='J0132-16',imagename='images/P4_dirty',niter=0,stokes='IV',cell='.00028arcsec')
clean(vis=p4ms,field='J0132-16',imagename='images/P4_clean_manual',interactive=True,npercycle=10,cell='.00028arcsec',cyclefactor=5,stokes='I',usescratch=True)

# time for amp+phase gaincal
ap5 = obsdir+'J0132.AP5'
gaincal(vis=p4ms,field='J0132-16',calmode='ap',refant='FD',caltable=ap5,solnorm=True)
plotcal(caltable=ap5, xaxis='scan',yaxis='phase',iteration='antenna',subplot=431)
applycal(vis=p4ms,field='J0132-16',gaintable=ap5,interp='nearest',flagbackup=False)


p5ms = obsdir+'P5.ms' # w/ amp cal applied
split(vis=p4ms,outputvis=p5ms,datacolumn='corrected')
clean(vis=p5ms,field='J0132-16',imagename='images/P5_dirty',niter=0,stokes='IV',cell='.00028arcsec')
clean(vis=p5ms,field='J0132-16',imagename='images/P5_clean_manual',interactive=True,npercycle=10,cell='.00028arcsec',cyclefactor=5,stokes='I',usescratch=True)

ap6 = obsdir+'J0132.AP6'
gaincal(vis=p5ms,field='J0132-16',calmode='ap',refant='FD',caltable=ap6,solnorm=True)
plotcal(caltable=ap6, xaxis='scan',yaxis='phase',iteration='antenna',subplot=431)
applycal(vis=p5ms,field='J0132-16',gaintable=ap6,interp='nearest',flagbackup=False)
# woohoo it's starting to converge!

# Let's use images/P5_clean_manual.model as a model for gaincal of full J0132 data set

### CONCATENATE ALL J0132 DATA ###

J0132mslist = ['/data/jrv/15A-416/UVCet/3/X/J0132.ms',
               '/data/jrv/15A-416/UVCet/4/X/J0132.ms',
               '/data/jrv/15A-416/UVCet/5/X/J0132.ms']

for jms in J0132mslist:
    uvfits = jms[:-3]+'.uvfits'
    importuvfits(fitsfile=uvfits,vis=jms)
    scan_reindex(jms)

msJ0132 = obsdir+'full_J0132.ms'
concat(vis=J0132mslist,concatvis=msJ0132)

ms_nocal = obsdir+'full_J0132_nocal.ms'
split(vis=msJ0132,outputvis=ms_nocal,datacolumn='data')


### Collect model images and use them to populate model column for all J0132 data ###

import os
os.system('cp -r images/P5_clean_manual.model imA.model')
os.system('cp -r /data/jrv/15A-416/UVCet/3/X/old/J0132/shortAP0.model imB.model')
# currently imB residuals have better RMS (1.6 mJy instead of 2.4) but the image uses a lot more data
# imA: model image based on bootstrapping gains from J0140 to adjacent J0132 scans in all 3 epochs
# imB: model image based on imaging J0132 by itself, imaging and calibrating short BLs then working outwards, first epoch only

msA = obsdir+'A.ms'
msB = obsdir+'B.ms'
modelA = 'imA.model'
modelB = 'imB.model'

split(vis=msJ0132,outputvis=msA,datacolumn='data')
split(vis=msJ0132,outputvis=msB,datacolumn='data')

ft(vis=msA,model=modelA,usescratch=True)
ft(vis=msB,model=modelB,usescratch=True)

# antenna 5 has weird amps throughout obs2
# it looks like antenna 5 is bad for whole obs2 - unless some amplitude calibration step in AIPS failed that can be fixed in AIPS?
#flagdata
# try calibrating and if it is still weird then we will flag ant 5 during obs2

gaincal(vis=msA,refant='FD',calmode='ap',solnorm=True,caltable='A.ap0')
gaincal(vis=msB,refant='FD',calmode='ap',solnorm=True,caltable='B.ap0')

plotcal(caltable='A.ap0', xaxis='scan',yaxis='phase',iteration='antenna',subplot=431,plotrange=[0,230,-180,180])
plotcal(caltable='B.ap0', xaxis='scan',yaxis='phase',iteration='antenna',subplot=431,plotrange=[0,230,-180,180])
# amplitude solns are smoother for ant6 for B, which is probably more realistic

applycal(vis=msA,gaintable='A.ap0')
applycal(vis=msB,gaintable='B.ap0')
# scan 135 and 136 have some bad data, remember to flag adjacent UVCet scans as well if applycal doesn't do it
# ant 5 is corrected to proper amps but more scatter --> questionable data quality

# now that we know that ant5 can be corrected, let's do phase-only gaincal first

gaincal(vis=msA,refant='FD',calmode='p',solnorm=True,caltable='A.p0')
gaincal(vis=msB,refant='FD',calmode='p',solnorm=True,caltable='B.p0')

plotcal(caltable='A.p0', xaxis='scan',yaxis='phase',iteration='antenna',subplot=431,plotrange=[0,230,-180,180])
plotcal(caltable='B.p0', xaxis='scan',yaxis='phase',iteration='antenna',subplot=431,plotrange=[0,230,-180,180])

# flag scan 136
flagdata(vis=msA,scan='136')
flagdata(vis=msB,scan='136')

applycal(vis=msA,gaintable='A.p0')
applycal(vis=msB,gaintable='B.p0')

# split out calibrated data; average in time and frequency to speed up future imaging and calibration
msAp0 = 'A.p0.ms'
msBp0 = 'B.p0.ms'
split(vis=msA,outputvis=msAp0,datacolumn='corrected',width=4,timebin='1000s')
split(vis=msB,outputvis=msBp0,datacolumn='corrected',width=4,timebin='1000s')

## Iteration 1 ##

clean(vis=msAp0,imagename='images2/Ap0_dirty',niter=0,stokes='IV',cell='.00028arcsec')
clean(vis=msBp0,imagename='images2/Bp0_dirty',niter=0,stokes='IV',cell='.00028arcsec')
# one is offset by 1.6 mas compared to other...
# Stokes V RMS is 997 uJy for both, images look identical

clean(vis=msAp0,imagename='images2/Ap0_clean',niter=5000,stokes='I',cell='.00028arcsec',npercycle=10,cyclefactor=5,interactive=True,usescratch=True)
clean(vis=msBp0,imagename='images2/Bp0_clean',niter=5000,stokes='I',cell='.00028arcsec',npercycle=10,cyclefactor=5,interactive=True,usescratch=True)
# residual RMS is 4.6-4.5 mJy for both

gaincal(vis=msAp0,refant='FD',calmode='p',solnorm=True,caltable='A.p1')
gaincal(vis=msBp0,refant='FD',calmode='p',solnorm=True,caltable='B.p1')

plotcal(caltable='A.p1', xaxis='scan',yaxis='phase',iteration='antenna',subplot=431,plotrange=[0,230,-15,15])
plotcal(caltable='B.p1', xaxis='scan',yaxis='phase',iteration='antenna',subplot=431,plotrange=[0,230,-15,15])
# biggest different: MK6 still moves by ~6 deg for A.p1, centered around 0 for B.p1
# since A.p1 is still moving in phase, let's iterate one more time

applycal(vis=msAp0,gaintable='A.p1')
applycal(vis=msBp0,gaintable='B.p1')

msAp1 = 'A.p1.ms'
msBp1 = 'B.p1.ms'
split(vis=msAp0,outputvis=msAp1,datacolumn='corrected')
split(vis=msBp0,outputvis=msBp1,datacolumn='corrected')


## Iteration 2 ##

clean(vis=msAp1,imagename='images2/Ap1_clean',niter=5000,stokes='I',cell='.00028arcsec',npercycle=10,cyclefactor=5,interactive=True,usescratch=True)
clean(vis=msBp1,imagename='images2/Bp1_clean',niter=5000,stokes='I',cell='.00028arcsec',npercycle=10,cyclefactor=5,interactive=True,usescratch=True)
# residual RMS is 4.4 mJy for both

gaincal(vis=msAp1,refant='FD',calmode='p',solnorm=True,caltable='A.p2')
gaincal(vis=msBp1,refant='FD',calmode='p',solnorm=True,caltable='B.p2')

plotcal(caltable='A.p2', xaxis='scan',yaxis='phase',iteration='antenna',subplot=431,plotrange=[0,230,-7,5])
plotcal(caltable='B.p2', xaxis='scan',yaxis='phase',iteration='antenna',subplot=431,plotrange=[0,230,-7,5])
# A has MK6 at about -4 degrees, B has all ants centered on zero
# -4 degrees isn't so bad, let's make next iteration amp+phase

applycal(vis=msAp1,gaintable='A.p2')
applycal(vis=msBp1,gaintable='B.p2')

msAp2 = 'A.p2.ms'
msBp2 = 'B.p2.ms'
split(vis=msAp1,outputvis=msAp2,datacolumn='corrected')
split(vis=msBp1,outputvis=msBp2,datacolumn='corrected')


## Iteration 3 ##

clean(vis=msAp2,imagename='images2/Ap2_clean',niter=5000,stokes='I',cell='.00028arcsec',npercycle=10,cyclefactor=5,interactive=True,usescratch=True)
clean(vis=msBp2,imagename='images2/Bp2_clean',niter=5000,stokes='I',cell='.00028arcsec',npercycle=10,cyclefactor=5,interactive=True,usescratch=True)
# residual RMS is 4.4 mJy for both

gaincal(vis=msAp2,refant='FD',calmode='ap',solnorm=True,caltable='A.ap3')
gaincal(vis=msBp2,refant='FD',calmode='ap',solnorm=True,caltable='B.ap3')

plotcal(caltable='A.ap3', xaxis='scan',yaxis='phase',iteration='antenna',subplot=431,plotrange=[0,230,-3,3])
plotcal(caltable='B.ap3', xaxis='scan',yaxis='phase',iteration='antenna',subplot=431,plotrange=[0,230,-3,3])
plotcal(caltable='A.ap3', xaxis='scan',yaxis='amp',iteration='antenna',subplot=431,plotrange=[0,230,0.4,1.3])
plotcal(caltable='B.ap3', xaxis='scan',yaxis='amp',iteration='antenna',subplot=431,plotrange=[0,230,0.4,1.3])
# A MK6 avg ~-1 deg, not bad - less scatter vs. spw than B for MK6 phases?
# A is now slightly flatter than B for MK6 amps --> phase-only selfcal helped!

applycal(vis=msAp2,gaintable='A.ap3')
applycal(vis=msBp2,gaintable='B.ap3')

msAap3 = 'A.ap3.ms'
msBap3 = 'B.ap3.ms'
split(vis=msAp2,outputvis=msAap3,datacolumn='corrected')
split(vis=msBp2,outputvis=msBap3,datacolumn='corrected')
# inspected uvdist vs. amp, phase in plotms; amp looks basically the same, phase totally different b/c of variation
#  in source position between model images


## Iteration 4 ##

clean(vis=msAap3,imagename='images2/Aap3_dirty',niter=0,stokes='IV',cell='.00028arcsec')
clean(vis=msBap3,imagename='images2/Bap3_dirty',niter=0,stokes='IV',cell='.00028arcsec')
# dirty Stokes V: both are 67-68 uJy, Aap3 is slightly lower peak flux but higher RMS (68 instead of 67)

clean(vis=msAap3,imagename='images2/Aap3_clean',niter=5000,stokes='I',cell='.00028arcsec',npercycle=10,cyclefactor=5,interactive=True,usescratch=True)
clean(vis=msBap3,imagename='images2/Bap3_clean',niter=5000,stokes='I',cell='.00028arcsec',npercycle=10,cyclefactor=5,interactive=True,usescratch=True)
# residual RMS is 1.2 mJy for A, 1.1 for B
# however I think I cleaned B a little deeper

# try multi-scale clean and nterms=2 to improve quality now
# recommendations for multiscale clean:
# multiscale = [0,a,3*a,9*a,... no larger than shortest baseline] where a = (clean beam size in pixels) = (3x1.2 mas)/(0.28 mas)
clean(vis=msBap3,imagename='images2/Bap3_msclean',multiscale=[0,4,12,36],niter=5000,stokes='I',cell='.00028arcsec',npercycle=10,cyclefactor=5,interactive=True,usescratch=True)
# that turned out really poorly
# try again with smaller initial mask? or different scales?
clean(vis=msBap3,imagename='images2/Bap3_msclean_v2',multiscale=[0,4,12],niter=5000,stokes='I',cell='.00028arcsec',npercycle=10,cyclefactor=5,interactive=True,usescratch=True) # went a lot better but residuals are worse than normal clean still
clean(vis=msBap3,imagename='images2/Bap3_nterms2',niter=5000,stokes='I',cell='.00028arcsec',npercycle=10,cyclefactor=5,interactive=True,usescratch=True,nterms=2) # same residual RMS as normal clean

# put normal clean model back in model column
ft(vis=msBap3,model='images2/Bap3_clean.model',usescratch=True)

# let's see if gain calibration can still improve
gaincal(vis=msAap3,refant='FD',calmode='ap',solnorm=True,caltable='A.ap4')
gaincal(vis=msBap3,refant='FD',calmode='ap',solnorm=True,caltable='B.ap4')

plotcal(caltable='A.ap4', xaxis='scan',yaxis='phase',iteration='antenna',subplot=431)
plotcal(caltable='B.ap4', xaxis='scan',yaxis='phase',iteration='antenna',subplot=431)
plotcal(caltable='A.ap4', xaxis='scan',yaxis='amp',iteration='antenna',subplot=431)
plotcal(caltable='B.ap4', xaxis='scan',yaxis='amp',iteration='antenna',subplot=431)
# only get 1-2 deg/% variation in phase/amp --> no need to continue iterating

# We now have model images for J0132:
# images2/Aap3_clean.model, Bap3_clean.model
os.system('cp -r images2/Aap3_clean.model J0132_refJ0140.model')
os.system('cp -r images2/Bap3_clean.model J0132_self.model')
# copied these to BV071/UVCet/J0132_refJ0140.model and J0132_self.model
