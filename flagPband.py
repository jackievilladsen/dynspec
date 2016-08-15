'''

flagPband.py

Run this from the directory with the ms file in it: e.g., /data/jrv/15A-416/YZCMi/1/P

Run this before calPband.py.

General approach:
Apply antpos+requantizer corrections and do Hanning smoothing.
Perform tfcrop autoflag twice, then merge to a single
spw and do initial BP calibration followed by rflag autoflag twice.
Average the flagged ms in time and frequency to create smaller ms that will be used by calPband.py for calibration.

Steps:
1 - diagnostic info - check for dead ants/swapped pols & verify good refant
2 - initial processing - flag zeros, a priori corrections (antpos, requantizer gain), Hanning smoothing
3 - auto flagging - important to do on full resolution data to maximize sensitivity

This version has a lot of plotting removed from it to save time.  See flagPband_plots.py for code to produce plots.
'''

# split command used to create test data set (4 scans: pcal, src, pcal, bpcal) - YZCMi_3P (actually copied from YZCMi_1)
# split(vis='YZCMi_1P.ms',outputvis='test.ms',scan='21~22,27~28',datacolumn='data')

from dynspec.pipeline_utils import *
import numpy as np
import os

### DEFINE VARIABLES ###
# field names, file names, phasecen, refant

names = get_names() # file names
sb = names['sb']
msfile = sb + '.ms'         # raw, full-size ms file
ms_hs = sb + '.hs.ms'       # hanning-smoothed, full-size ms file (w/ antpos table applied if it exists)
name1spw = sb + '.1spw'
ms1spw = name1spw + '.ms'
smallms = sb + '.small.ms'
small2ms = sb + '.small2.ms'
smallname = sb + '.small'

if not os.path.exists(msfile):
    extract_ms(msfile)

fields = get_fields(msfile,vishead) # field names
srcname = fields.get('src')
bpcal = fields.get('bpcal')
gcal = fields.get('gcal')

refanten = 'ea06' # reference antenna - ea06 is central ant on north arm - check plots after to make sure this is okay
phasecen=get_phasecen() # read phase center from file
plotdir = get_plot_dir()

## Diagnostic info ##

listobs(vis=msfile,listfile=msfile+'.listobs')

plotants(vis=msfile,figfile=plotdir+'plotants.png') # plot ant locations to choose refant - should be near array center
print 'Check plots/plotants.png to verify that default refant',refanten,'is near array center.'

## Initial processing ##

print 'Flagging zeros (correlator errors)...'
flagdata(vis=msfile,mode='clip',clipzeros=True) # zeros indicate correlator malfunction - should be ~1% of data or less

gaintable = init_cal(msfile,gencal,plotcal)
#print 'Applying gaintables',gaintable, 'to',msfile
#applycal(vis=msfile,gaintable=gaintable)

# inspect BP cal in plotms to look for dead or swapped ants
amax = amp_plot_max(msfile,visstat,field=bpcal)
bmax = amax * 2.5
plotms(vis=msfile,field=bpcal,xaxis='freq',yaxis='amp',coloraxis='corr',iteraxis='baseline',avgtime='1e8',avgscan=True,plotrange=[0.2,0.5,0,amax],gridrows=3,gridcols=3,plotfile=plotdir+'bpcal_raw.png',exprange='all',showgui=False,overwrite=True)
print 'Check plots/bpcal_raw.png to look for dead ants or swapped polarizations and confirm whether', refanten,'is a good refant (low RFI, no weird behavior).  XX and YY are purple and orange, XY and YX are black and pink (respectively).'

# Hanning smooth on corrected data
print 'Hanning smoothing',msfile,'to create',ms_hs
hanningsmooth(vis=msfile,outputvis=ms_hs,field=bpcal+','+srcname)

# remove msfile (to save room) - it has CORRECTED_DATA column so is 2x the size of ms_hs
# only do this if the tarfile exists
if os.path.exists(msfile+'.tar'):
    rmtables(msfile)


## AUTOFLAG: TFCROP on bandpass calibrator ##

# print summary of flagged data
print 'Before auto RFI flagging (after flagging zeros):'
summary_1 = field_flag_summary(ms_hs,flagdata)
np.save('summary1.npy',summary_1)

# run auto-flagger w/ tfcrop (5-sigma threshold - fairly high)
# maxnpieces=5 is lower than default but we have fewer channels per spw than at higher freq
print 'First tfcrop...'
flagdata(vis=ms_hs,mode='tfcrop',timecutoff=5.0,freqcutoff=5.0,maxnpieces=5)

# run tfcrop again with same params to get more RFI - typically adds 1-2% more flagged
print 'Second tfcrop...'
flagdata(vis=ms_hs,mode='tfcrop',timecutoff=5.0,freqcutoff=5.0,maxnpieces=5)

print 'Flagged data summary after two tfcrop autoflags:'
summary_2 = field_flag_summary(ms_hs,flagdata)
np.save('summary2.npy',summary_2)


## COMBINE SPWS ##
# remove spw boundaries and apply initial bpcal to help with auto RFI flagging

print 'Combining all spws from',ms_hs,'to create',ms1spw
mstransform(vis=ms_hs,outputvis=ms1spw,combinespws=True,datacolumn='data') # remove spw boundaries - takes ~ 1 hr for 50 GB SB
rmtables(ms_hs) # save disk space

fluxmodel=bpcal+'_L.im'
print 'Loading virtual flux model',fluxmodel
setjy(vis=ms1spw,field=bpcal,standard='Scaife-Heald 2012',model=fluxmodel) # load bpcal model to acct for structure, spec index

print 'Solving for initial bandpass w/ minsnr=0.1.'
bandpass(vis=ms1spw,field=bpcal,caltable=name1spw+'.B0',refant=refanten,minsnr=0.1) # run BP cal on 1-spw BP ms
# minsnr=0.1 --> get solns even for channels w/ lots of RFI

print 'Check plots/'+name1spw+'.B0[amp/phase].png to inspect initial BP cal.'
plotcal(caltable=name1spw+'.B0',xaxis='freq',yaxis='amp',iteration='antenna',subplot=441,showgui=False,figfile=plotdir+name1spw+'.B0amp.png')
plotcal(caltable=name1spw+'.B0',xaxis='freq',yaxis='phase',iteration='antenna',subplot=441,showgui=False,figfile=plotdir+name1spw+'.B0phase.png')

# apply BP cal table
print 'Applying initial bandpass solution so spectrum is smooth for rflag.'
applycal(vis=ms1spw,gaintable=[name1spw+'.B0'])
# backup flag file: before_applycal_1

cmax = amp_plot_max(ms1spw,visstat,'corrected',field=bpcal)
plotms(vis=ms1spw,field=bpcal,ydatacolumn='corrected',xaxis='freq',yaxis='amp',coloraxis='corr',correlation='XX,YY',plotrange=[0.2,0.5,0,cmax],plotfile=plotdir+'post_B0.png',showgui=False,overwrite=True,avgtime='1e8')
print 'Check plots/post_B0.png to confirm that initial BP cal has worked.'

### AUTOFLAG: RFLAG ###

# mode='rflag': use sliding window in time & freq to flag (whereas tfcrop just uses avg & std over whole time/all freqs)
# winsize = 5: # of timesteps in sliding window
# run 2x to get some more RFI (1st round gets a lot more but 2nd round's cleanup seems to be important)
print '1st rflag...'
flagdata(vis=ms1spw,mode='rflag',datacolumn='corrected',winsize=5)
print '2nd rflag...'
flagdata(vis=ms1spw,mode='rflag',datacolumn='corrected',winsize=5,flagbackup=False)

plotms(vis=ms1spw,field=bpcal,xaxis='freq',yaxis='amp',ydatacolumn='corrected',coloraxis='corr',correlation='XX,YY',plotrange=[0.2,0.5,0,cmax],plotfile=plotdir+'post_rflag2.png',showgui=False,overwrite=True,avgtime='1e8')
print 'Check plots/post_rflag2.png to see RFI in BP cal after second rflag.'

print 'Flagged data summary after second rflag:'
summary_3 = field_flag_summary(ms1spw,flagdata)
np.save('summary3.npy',summary_3)

# back up final flag state - flags will also be saved by applycal in calPband.py but
#  this makes it easier to identify flagversion associated with end of flagPband.py
flagmanager(vis=ms1spw,mode='save',versionname='post_flagPband')
