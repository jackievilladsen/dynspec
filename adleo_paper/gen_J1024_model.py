from dynspec.tbavg import scan_reindex
import os

obsdir = '/data/jrv/BV071/ADLeo/'
pcal_dir = obsdir+'gen_pcal_model/'
refms = 'J0132_J0140.ms'
ap0 = obsdir+'J0140.AP0_v2'
refms_nocal = 'J0132_J0140.nocal.ms'

os.chdir(pcal_dir)

mslist = []
obslist = ['3','4','5']
for obs in obslist:
    J1024_uvfits = obsdir+obs+'/J1024.uvfits'
    J1024ms = obsdir+obs+'/J1024.ms'
    importuvfits(fitsfile=J1024_uvfits,vis=J1024ms)
    scan_reindex(J1024ms)
    mslist.append(J1024ms)

msJ1024 = pcal_dir+'full_J1024.ms'
concat(vis=mslist,concatvis=msJ1024)

#split to average to shorter time
pcalms = pcal_dir + 'full_J1024_avg.ms'
split(vis=msJ1024,outputvis=pcalms,timebin='1000s',width=4,datacolumn='data')
# inspecting in plotms: obs 5 has a lot of bad data (low amps) ... was src rising and setting? what is going on? grrr
# NL, HN, SC all have really low amps some of the time
# LA has low amps for obs 4 (scans 78~154)
flagdata(vis=pcalms,scan='155~229',antenna='HN,NL,SC')
flagdata(vis=pcalms,scan='78~154',antenna='LA')

# image J1024
clean(vis=pcalms,imagename='images/dirty0',cell='0.00028arcsec',niter=0,interactive=True,stokes='IV')
clean(vis=pcalms,imagename='images/image0_v2',cell='0.00028arcsec',niter=5000,interactive=True,cyclefactor=5,npercycle=10)

p0 = 'J1024.p0'
gaincal(vis=pcalms,refant='FD',solnorm=True,calmode='p',caltable=p0)
plotcal(caltable=p0, xaxis='scan',yaxis='phase',iteration='antenna',subplot=431)
applycal(vis=pcalms,gaintable=p0)

p0ms = 'p0.ms'
split(vis=pcalms,outputvis=p0ms,datacolumn='corrected')

## Iteration 1 ##

clean(vis=p0ms,imagename='images/p0',cell='0.00028arcsec',niter=5000,interactive=True,cyclefactor=5,npercycle=10)

p1 = 'J1024.p1'
gaincal(vis=p0ms,refant='FD',solnorm=True,calmode='p',caltable=p1)
plotcal(caltable=p1, xaxis='scan',yaxis='phase',iteration='antenna',subplot=431)
applycal(vis=p0ms,gaintable=p1)

p1ms = 'p1.ms'
split(vis=p0ms,outputvis=p1ms,datacolumn='corrected')

## Iteration 2 ##

clean(vis=p1ms,imagename='images/p1',cell='0.00028arcsec',niter=5000,interactive=True,cyclefactor=5,npercycle=10)

ap2 = 'J1024.ap2'
gaincal(vis=p1ms,refant='FD',solnorm=True,calmode='ap',caltable=ap2)
plotcal(caltable=ap2, xaxis='scan',yaxis='phase',iteration='antenna',subplot=431)
plotcal(caltable=ap2, xaxis='scan',yaxis='amp',iteration='antenna',subplot=431)
applycal(vis=p1ms,gaintable=ap2)

ap2ms = 'ap2.ms'
split(vis=p1ms,outputvis=ap2ms,datacolumn='corrected')

## Iteration 3 ##

clean(vis=ap2ms,imagename='images/ap2',cell='0.00028arcsec',niter=5000,interactive=True,cyclefactor=5,npercycle=10)

ap3 = 'J1024.ap3'
gaincal(vis=ap2ms,refant='FD',solnorm=True,calmode='ap',caltable=ap3)
plotcal(caltable=ap3, xaxis='scan',yaxis='phase',iteration='antenna',subplot=431)
plotcal(caltable=ap3, xaxis='scan',yaxis='amp',iteration='antenna',subplot=431)
applycal(vis=ap2ms,gaintable=ap3)

ap3ms = 'ap3.ms'
split(vis=ap2ms,outputvis=ap3ms,datacolumn='corrected')

## Iteration 4 ##

clean(vis=ap3ms,imagename='images/ap3',cell='0.00028arcsec',niter=5000,interactive=True,cyclefactor=5,npercycle=10)

ap4 = 'J1024.ap4'
gaincal(vis=ap3ms,refant='FD',solnorm=True,calmode='ap',caltable=ap4)
plotcal(caltable=ap4, xaxis='scan',yaxis='phase',iteration='antenna',subplot=431)
plotcal(caltable=ap4, xaxis='scan',yaxis='amp',iteration='antenna',subplot=431)
applycal(vis=ap3ms,gaintable=ap4)

ap4ms = 'ap4.ms'
split(vis=ap3ms,outputvis=ap4ms,datacolumn='corrected')

## Iteration 5 ##

clean(vis=ap4ms,imagename='images/ap4',cell='0.00028arcsec',niter=5000,interactive=True,cyclefactor=5,npercycle=10)
# ok this is pretty good

import os
os.system('cp -r images/ap4.model /data/jrv/BV071/ADLeo/J1024.model')
