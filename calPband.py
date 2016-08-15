'''

calPband.py

Run this from the directory with the ms file in it: e.g., /data/jrv/15A-416/YZCMi/1/P

Run this after running flagPband.py.

General approach:
Solve for calibration tables based on BPcal field, apply them to target field, and create
a calibrated target-only ms.

Calibration steps:
- setjy was already run by flagPband.py
- delay: K0
- gaincal: G0
- bandpass: B1  (B0 was initial bandpass in flagPband.py)
- polarization
---- cross-pol: Kc0
---- leakage: Df0
- apply caltables to all fields
- small image of gcal field to confirm gcal at expected position
- split calibrated target field

'''

from dynspec.pipeline_utils import *
import numpy as np

print 'Running calPband.py, estimated duration of ~... min.'

### DEFINE VARIABLES ###
# field names, file names, phasecen, refant

names = get_names() # file names
sb = names['sb']
name1spw = sb + '.1spw'
ms1spw = name1spw + '.ms'  # full-size ms produced by flagPband.py
srcms = sb + '.src.ms'

fields = get_fields(ms1spw,vishead) # field names
srcname = fields.get('src')
bpcal = fields.get('bpcal')
gcal = fields.get('gcal')

refanten = 'ea06' # reference antenna - ea06 is central ant on north arm - check plots after to make sure this is okay
phasecen=get_phasecen() # read phase center from file
plotdir = get_plot_dir()

gaintableKB = [name1spw+'.K0',name1spw+'.B1']
gaintableKBKc = [name1spw+'.K0',name1spw+'.B1',name1spw+'.Kc0']
gaintableKBKcDf = [name1spw+'.K0',name1spw+'.B1',name1spw+'.Kc0',name1spw+'.Df0']
gaintableKBKcDfB = [name1spw+'.K0',name1spw+'.B1',name1spw+'.Kc0',name1spw+'.Df0',name1spw+'.B2']

### DELAY CALIBRATION ### - duration: ~1 min

print 'Solving for delays...'
gaincal(vis=ms1spw,gaintype='K',field=bpcal,caltable=name1spw+'.K0',refant=refanten,minsnr=3.0,parang=True)
plotcal(caltable=name1spw+'.K0',xaxis='antenna',yaxis='delay',figfile=plotdir+'delay.png',showgui=False)
print 'Check plots/delay.png - all delays should be < few ns.'


### BANDPASS CALIBRATION ### - duration: ~3.5 min

print 'Solving for bandpass table B1 for',ms1spw
bandpass(vis=ms1spw,caltable=name1spw+'.B1',field=bpcal,refant=refanten,gaintable=[name1spw+'.K0'])
flagdata(vis=name1spw+'.B1',mode='rflag',winsize=3,correlation='ARG_ALL',datacolumn='CPARAM',freqdevscale=4.0) # ARG_ALL = phase
plotcal(caltable=name1spw+'.B1',xaxis='freq',yaxis='amp',subplot=441,iteration='antenna',showgui=False,figfile=plotdir+'B1amp.png')
plotcal(caltable=name1spw+'.B1',xaxis='freq',yaxis='phase',subplot=441,iteration='antenna',showgui=False,figfile=plotdir+'B1phase.png')
print 'Check plots/B1[amp/phase].png for bandpass solutions B1 after rflag on phase.'


### POLARIZATION CALIBRATION: 3C147 ### - total duration: ~3 min

pol_dict = {'3C48':'0.3%','3C138':'5.6%','3C147':'<0.05%','3C286':'8.6%'}
p = pol_dict.get(bpcal,'[unknown source name]')
print 'Perley+13:',bpcal,'1.05-GHz linear polarization is',p+'; assuming unpolarized for P-band polarization calibration.', \
      'The main source of error in pol cal may be time variability.'

## Cross-delay (Kc0) ## - duration: ~ 1 min
gaincal(vis=ms1spw,field=bpcal,caltable=name1spw+'.Kc0',gaintype='KCROSS',refant=refanten,parang=True,gaintable=gaintableKB)
plotcal(caltable=name1spw+'.Kc0',xaxis='antenna',yaxis='delay',showgui=False,figfile=plotdir+'crossdelay.png')
print 'Check plots/crossdelay.png to see delay between 2 polns of refant', refanten

## D-terms: frequency dependent cross-pol leakage ## - duration: ~2 min
polcal(vis=ms1spw,poltype='Df',field=bpcal,caltable=name1spw+'.Df0',refant=refanten,preavg=1.e6,gaintable=gaintableKBKc)
flagdata(vis=name1spw+'.Df0',mode='clip',clipminmax=[-0.1,0.4],datacolumn='CPARAM') # flag high-amplitude D-terms
plotcal(caltable=name1spw+'.Df0',xaxis='freq',yaxis='amp',iteration='antenna',subplot=441,showgui=False,figfile=plotdir+'leakage.png')
print 'Check plots/leakage.png for D terms after autoflag (clip amplitude > 0.4).  Should be <<1.'


### APPLY CALIBRATIONS TO ALL FIELDS ### - duration: ~1+0.5*Nscans min

print '% flagged in',ms1spw,'before applying calibration:'
x = field_flag_summary(ms1spw,flagdata)

# duration: ~0.5 min/scan
applycal(vis=ms1spw,applymode='calflagstrict',gaintable=gaintableKBKcDf)

print '% flagged in',ms1spw,'after updated caltables (based on',bpcal+') applied to all fields:'
summary_4 = field_flag_summary(ms1spw,flagdata)
np.save('summary4.npy',summary_4)

# duration: ~1 min
plotms(vis=ms1spw,field=bpcal,ydatacolumn='corrected',xaxis='freq',yaxis='amp',coloraxis='corr',plotrange=[0.2,0.5,0,150],plotfile=plotdir+'BPcal_postcal.png',showgui=False,overwrite=True,avgtime='1e8')
print 'Check plots/BPcal_postcal.png to confirm that fully-calibrated BPcal looks good. XX and YY are purple and orange, XY and YX are black and pink (respectively).'


### FLAG w/ final BP calibration ### - duration: ~1 min + 3.5 min * Nscans

# run autoflagger one more time with final calibrations based on BPcal - duration: 3.5 min/target scan, 1 min/pcal scan
flagdata(vis=ms1spw,mode='rflag',datacolumn='corrected',winsize=5,display='report')
flagdata(vis=ms1spw,mode='extend',growaround=True) # growaround=True: flag pts w/ >4 flagged pts adjacent in dynspec

print '% flagged in',ms1spw,'after rflag+growaround on calibrated data:'
summary_5 = field_flag_summary(ms1spw,flagdata)
np.save('summary5.npy',summary_5)

plotms(vis=ms1spw,field=bpcal,ydatacolumn='corrected',xaxis='freq',yaxis='amp',coloraxis='corr',correlation='XX,YY',plotrange=[0.2,0.5,0,150],plotfile=plotdir+'BPcal_postcal_rflag.png',showgui=False,overwrite=True,avgtime='1e8')
print 'Check plots/BPcal_postcal_rflag.png to confirm that fully-calibrated BPcal looks good after rflag+growaround.'


### IMAGE GAIN CALIBRATOR TO CONFIRM LOCATION ###

clean(vis=ms1spw,field=gcal,imagename='gcal',cell='3arcsec',niter=0) # duration:
imview(raster={'file':'gcal.image'},out=plotdir+'gcal.png')
print 'Check plots/gcal.png to confirm that gain calibrator is roughly in center of image and not smeared, implying that calibrations',\
      'derived from BP calibrator are fairly good throughout entire SB.'

