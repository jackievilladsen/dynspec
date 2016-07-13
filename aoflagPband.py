'''

aoflagPband.py

Run this from the directory with the ms file in it: e.g., /data/jrv/15A-416/YZCMi/1/P

Run this before calPband.py.

General approach:
Apply antpos and requantizer corrections, do Hanning smoothing then run aoflagger once.
Then merge to single spw and apply initial bandpass - run aoflagger again after this?
Then average flagged ms in time and frequency to create smaller
ms that will be used by calPband.py for calibration.
'''

from dynspec.pipeline_utils import *
import numpy as np
 
### DEFINE VARIABLES ###
# field names, file names, phasecen, refant

names = get_names() # file names
sb = names['sb']
msfile = sb + '.ms'         # raw, full-size ms file
ms_hs = sb + '.hs.ms'       # hanning-smoothed, full-size ms file (w/ antpos table applied if it exists)
name1spw = sb + '.1spw'
ms1spw = name1spw + '.ms'
smallms = sb + '.small.ms'
smallname = sb + '.small'

fields = get_fields(msfile,vishead) # field names
srcname = fields.get('src')
bpcal = fields.get('bpcal')
gcal = fields.get('gcal')

refanten = 'ea06' # reference antenna - ea06 is central ant on north arm - check plots after to make sure this is okay

phasecen=get_phasecen() # read phase center from file

plotdir = get_plot_dir()

### FULL-SIZE MS ###

# steps:
# 1 - diagnostic info (check for dead ants/swapped pols & identify refant)
# 2 - initial processing (flag zeros, antpos corrections, Hanning smooth)
# 3 - auto flagging (important to do on full resolution data to maximize sensitivity)


# create BPcal-only ms for testing purposes
# split(vis='YZCMi_1P.ms',outputvis='bpcal.ms',field='3C147',datacolumn='data')
# zeros have been flagged but otherwise this is raw data set

## Diagnostic info ##

listobs(vis=msfile,listfile=msfile+'.listobs')

plotants(vis=msfile,figfile=plotdir+'plotants.png') # plot ant locations to choose refant - should be near array center

# inspect BP cal in plotms - is flux detected in XX and YY for all BLs and not in XY or YX?
amax = amp_plot_max(msfile,visstat)
plotms(vis=msfile,field=bpcal,xaxis='freq',yaxis='amp',coloraxis='corr',iteraxis='baseline',avgtime='1e8',avgscan=True,plotrange=[0.2,0.5,0,amax],gridrows=3,gridcols=3,plotfile=plotdir+'bpcal_raw.png',showgui=False,exprange='all',overwrite=True)
print 'Check plots/bpcal_raw.png to look for dead ants or swapped polarizations, and confirm whether',refanten,'is a good refant (low RFI, no weird behavior.  XX and YY are purple and orange, XY and YX are black and pink (respectively).'


## Initial processing ##

flagdata(vis=msfile,mode='clip',clipzeros=True) # zeros indicate correlator malfunction - should be ~1% of data or less

gaintable = init_cal(msfile,gencal,plotcal)
print 'Applying gaintables',gaintable, 'to',msfile
applycal(vis=msfile,gaintable=gaintable)

# inspect BP cal in plotms - is flux detected in XX and YY for all BLs and not in XY or YX?
amax = amp_plot_max(msfile,visstat,datacolumn='corrected')
bmax = amax * 2.5
plotms(vis=msfile,field=bpcal,xaxis='freq',yaxis='amp',ydatacolumn='corrected',coloraxis='corr',iteraxis='baseline',avgtime='1e8',avgscan=True,plotrange=[0.2,0.5,0,amax],gridrows=3,gridcols=3,plotfile=plotdir+'bpcal_corr.png',exprange='all',showgui=False,overwrite=True)
print 'Check plots/bpcal_corr.png to see how BPcal looks after initial cals (antpos, requantizer - no TEC yet).  XX and YY are purple and orange, XY and YX are black and pink (respectively).'

# Hanning smooth (using datacolumn='corrected' means antpos table will be applied to new ms, if it exists - otherwise will use raw data column)
hanningsmooth(vis=msfile,outputvis=ms_hs,datacolumn='corrected')


## AUTOFLAG: AOFLAGGER ##

# plotrange is no longer good - probably b/c of requantizer gains
plotms(vis=ms_hs,field=bpcal,xaxis='freq',yaxis='amp',coloraxis='spw',correlation='XX,YY',plotrange=[0.2,0.5,0,bmax],plotfile=plotdir+'preflag.png',showgui=False,overwrite=True) # plot amp vs. frequency of BP cal
print 'Check plots dir for preflag.png to see RFI in BP cal before auto flagging.'

# print summary of flagged data
print 'Before auto RFI flagging (after flagging zeros):'
summary_1 = field_flag_summary(ms_hs,flagdata)
np.save('summary1.npy',summary_1)

# transfer ms_hs to manwe, run aoflagger, transfer back

plotms(vis=ms_hs,field=bpcal,xaxis='freq',yaxis='amp',coloraxis='spw',correlation='XX,YY',plotrange=[0.2,0.5,0,bmax],plotfile=plotdir+'flagdata1.png',showgui=False,overwrite=True)
print 'Check plots dir for flagdata1.png to see RFI in BP cal after one round auto-flagging.'

# summarize flagged data
print 'Flagged data summary after one aoflagger run (w/o bandpass correction):'
summary_2 = field_flag_summary(ms_hs,flagdata)
np.save('summary2.npy',summary_2)


## COMBINE SPWS ##
# remove spw boundaries and apply initial bpcal to help with auto RFI flagging

mstransform(vis=ms_hs,outputvis=ms1spw,combinespws=True,datacolumn='data') # remove spw boundaries - takes ~ 1 hr for 50 GB SB

setjy(vis=ms1spw,field=bpcal,standard='Scaife-Heald 2012',usescratch=True,model=bpcal+'_L.im') # load bpcal model to acct for structure, spec index

bandpass(vis=ms1spw,field=bpcal,caltable=name1spw+'.B0',refant=refanten,minsnr=0.1) # run BP cal on 1-spw BP ms
# minsnr=0.1 --> get solns even for channels w/ lots of RFI

print 'Check plots/'+name1spw+'.B0amp.png to inspect initial BP cal.'
plotcal(caltable=name1spw+'.B0',xaxis='freq',yaxis='amp',iteration='antenna',subplot=441,showgui=False,figfile=plotdir+name1spw+'.B0amp.png')

# apply BP cal table
applycal(vis=ms1spw,gaintable=[name1spw+'.B0'])
# backup flag file: before_applycal_1

cmax = amp_plot_max(ms1spw,visstat,'corrected')
plotms(vis=ms1spw,field=bpcal,ydatacolumn='corrected',xaxis='freq',yaxis='amp',coloraxis='corr',correlation='XX,YY',plotrange=[0.2,0.5,0,cmax],plotfile=plotdir+'post_B0.png',showgui=False,overwrite=True,avgtime='1e8')
print 'Check plots/post_B0.png to confirm that initial BP cal has worked.'


# add a second round of aoflagger here??


### CREATE SMALLER MS BY AVG'ING OVER TIME AND FREQ ###

# avg over 8 channels --> output channels are 1 MHz
# avg over 15s intervals
split(vis=ms1spw,outputvis=smallms,datacolumn='data',width=8,timebin='15s')

# backup small ms as tar file
os.system('tar -cvf ' + smallms + '.tar ' + smallms)

