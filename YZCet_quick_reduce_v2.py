split(vis='19B-222.sb37554451.eb37565636.58817.005467372685.ms',datacolumn='corrected',field='YZCet',outputvis='YZCet_small.ms',timebin='30s',width=4)

split(vis='19B-222.sb37555306.eb37565660.58818.00239415509.ms/',datacolumn='corrected',field='YZCet',outputvis='YZCet_small.ms',timebin='30s',width=4)

split(vis='19B-222.sb37555388.eb37565973.58818.99868625.ms/',datacolumn='corrected',field='YZCet',outputvis='YZCet_small.ms',timebin='30s',width=4)

concat(vis=['19B-222_2019_11_30_T06_49_39.524/YZCet_small.ms','19B-222_2019_12_01_T06_49_46.983/YZCet_small.ms','19B-222_2019_12_02_T06_49_27.743/YZCet_small.ms'],concatvis='YZCet_small_all.ms')

myvis='YZCet_small_all.ms'

### DIRTY IMAGE AND LOCATE SOURCES ###
tclean(vis=myvis,imagename='images/dirty_IQUV_v0/dirty_IQUV_v0',cell='3arcsec',imsize=1024,stokes='IQUV',niter=0,pblimit=-0.001)
# beam size: 37" x 25"
# fairly strong Stokes V (leakage? true pol?) for PMN: 4.3 mJy V vs. 148 mJy I --> ~3%
#  I'm guessing this is leakage because it's not in the center of the FOV - it's
#    more than I was hoping for :/

# source coords:
# PMN J0112-1658: 01 12 01.330 -16 58 09.90
# YZCet (epoch J2019.92, for Dec 1 2019): 01 12 32.31034 -16 59 43.6567
#     - RA: +237.935", dec: -48.9088" from current phase center
#     Nothing visible in dirty I, possible barest hint of something in V but could just be noise
#     (Going to shorter timescales might turn it into something significant)
pcen_PMN = 'J2000 01h12m01.330 -16d58m09.90'
pcen_YZCet = 'J2000 01h12m32.31034 -16d59m43.6567'
# in dirty image, PMN source position is consistent with these coords
# ran gaussfit on dirty image --> source size is marginally resolved (~1/2 of size of beam)
#   - this may just be due to it being a dirty image and having sidelobes, but even if it is
#     real, this is small enough that we probably won't gain that much from multi-scale clean
#     ---> don't use multi-scale clean at first

# use plotms to check frequency dependence of src - need nterms=3?
# choose only one scan (scan 4 is first scan)
plotms(vis=myvis,xaxis='freq',yaxis='real',ydatacolumn='data',avgbaseline=True,avgtime='1e8',correlation='RR,LL',coloraxis='corr',plotfile='plots/spectra/PMN.dirty.png',shift=[-208.683,45.268],scan='4')
# in the long run we will want nterms=3 I think, but nterms=2 should give us a decent first pass



### SHIFT PHASE CENTER TO TARGET ###

fixvis(vis=myvis,phasecenter=pcen_YZCet)


### INITIAL CLEAN IMAGE ###

myvis='YZCet_small_all.ms'
clearcal(myvis) # erase old selfcal from previous reduction

# use cleaning approach from YZCet_followup2

# clean full data set, non-interactively
# but use mask from previous reduction of this data set
tclean(vis=myvis,imagename='images/clean_I_n3_ms/clean_I_n3_ms',cell='1arcsec',imsize=3000,stokes='I',niter=10000,pblimit=-0.001,deconvolver='mtmfs',nterms=3,interactive=False,scales=[0,10,30],mask='old_images/clean_I_n3_scp0/clean_I_n3_scp0.mask',threshold='1mJy')
# made image too large - use cell='3arcsec' and imsize=1024 after this


### POPULATE MODEL COLUMN WITH MODEL IMAGE ###

tclean(vis=myvis,imagename='images/clean_I_n3_ms/clean_I_n3_ms',cell='1arcsec',imsize=3000,stokes='I',niter=0,pblimit=-0.001,deconvolver='mtmfs',nterms=3,interactive=False,scales=[0,10,30],threshold='1mJy',restart=True,savemodel='modelcolumn',calcpsf=False,calcres=False)



### SELF-CAL ###

# flag ea22 on day 2 (scan 92-265) because in previous reduction it had wild phases then
flagdata(vis=myvis,antenna='ea22',scan='92~265')

# combine spw's for selfcal
# (come back and apply solns to non-combined ms once done with all selfcal)
# (there was an error when I tried using combine='spw' in gaincal so this is the way to go for now)
vis1spw='YZCet_small_all.1spw.ms'
mstransform(vis=myvis,outputvis=vis1spw,combinespws=True,datacolumn='all')

scp0 = vis1spw + '.solint_scans.ScP0'
gaincal(vis=vis1spw,caltable=scp0,solint='inf',refant='ea24',calmode='p')
# initially tried ea01 as a refant but the solutions had more scatter and it was missing sometimes
# also tried 30-s integrations but in my opinion it had too much scatter

plotcal(caltable=scp0,xaxis='scan',yaxis='phase',iteration='antenna',subplot=551,plotrange=[0,0,-20,20])
# looks pretty good - could even go to a bit longer interval but this time interval takes out occasional
#  more rapid variations

applycal(vis=vis1spw,gaintable=scp0,applymode='calflag')
# 34.2 -> 34.5% flagged


### IMAGE POST-SELFCAL ###

# I initially ran it on the wrong data set! - IV.1mJy.mask is the mask from that run,
#  which includes a region in Stokes V over the bright background source
tclean(vis=vis1spw,imagename='images/clean_IV_n3_ms_ScP0/clean_IV_n3_ms_ScP0',cell='3arcsec',imsize=1024,stokes='IV',niter=10000,pblimit=-0.001,deconvolver='mtmfs',nterms=3,interactive=False,scales=[0,10,30],mask='images/IV.1mJy.mask',threshold='1mJy')

# resume clean, interactive to go deeper
tclean(vis=vis1spw,imagename='images/clean_IV_n3_ms_ScP0/clean_IV_n3_ms_ScP0',cell='3arcsec',imsize=1024,stokes='IV',niter=10000,pblimit=-0.001,deconvolver='mtmfs',nterms=3,interactive=True,scales=[0,10,30],restart=True)
# cleaned down to 0.28 mJy then stopped

# ft to populate model column
tclean(vis=vis1spw,imagename='images/clean_IV_n3_ms_ScP0/clean_IV_n3_ms_ScP0',cell='3arcsec',imsize=1024,stokes='IV',niter=0,pblimit=-0.001,deconvolver='mtmfs',nterms=3,interactive=False,scales=[0,10,30],restart=True,calcpsf=False,calcres=False,savemodel='modelcolumn')


### AMP + PHASE SELF-CAL ###

scAp0 = vis1spw + '.solint_scan.ScAP0'
gaincal(vis=vis1spw,caltable=scAp0,solint='inf',refant='ea24',calmode='ap',solnorm=True)
# no failed solutions

plotcal(caltable=scAp0,xaxis='scan',yaxis='phase',iteration='antenna',subplot=551,plotrange=[0,0,-20,20])
plotcal(caltable=scAp0,xaxis='scan',yaxis='amp',iteration='antenna',subplot=231) #,plotrange=[0,0,-50,50])
# some amps looked weird but majority okay - flagged weird ones
# In general the amps are much grosser than the phases - not a great data set based on these amps
# flagged:
#   - low amps on first scan (some antennas)
#   - weird ants on ea16 & ea17
# did not flag, but might want to: ea26 has low amps throughout (~70%)
# After applycal, run extend over pols since sometimes I maybe only flagged R or L gain, not both

applycal(vis=vis1spw,gaintable=scAp0,applymode='calflag')
# flag backup: applycal_2
# 34.5 -> 35.2% flagged

flagdata(vis=vis1spw,mode='extend',extendpols=True,growtime=100,growfreq=100)
# 35.2 -> 36.5% flagged
# remember to do this for myvis

# now apply to myvis
applycal(vis=myvis,gaintable=scAp0,applymode='calflag',spwmap=16*[0]) # applycal_1 - 34.0 -> 35.1% flagged
flagdata(vis=myvis,mode='extend',extendpols=True,growtime=100,growfreq=100) # flagdata_2


### IMAGE MULTI-SPW DATASET POST-SELFCAL ###

tclean(vis=myvis,imagename='images/clean_IV_n3_ms_ScAp0_myvis/clean_IV_n3_ms_ScAp0_myvis',cell='3arcsec',imsize=1500,stokes='IV',niter=10000,pblimit=-0.001,deconvolver='mtmfs',nterms=3,interactive=True,scales=[0,10,30],mask='images/clean_IV_n3_ms_ScP0/clean_IV_n3_ms_ScP0.mask',threshold='0.3mJy')

# restart and go deeper - interactive
tclean(vis=myvis,imagename='images/clean_IV_n3_ms_ScAp0_myvis/clean_IV_n3_ms_ScAp0_myvis',cell='3arcsec',imsize=1500,stokes='IV',niter=10000,pblimit=-0.001,deconvolver='mtmfs',nterms=3,interactive=True,scales=[0,10,30],restart=True)
# stopped at ~50 uJy

### FT, UVSUB, TBAVG ###

tclean(vis=myvis,imagename='images/clean_IV_n3_ms_ScAp0_myvis/clean_IV_n3_ms_ScAp0_myvis',cell='3arcsec',imsize=1500,stokes='IV',niter=0,pblimit=-0.001,deconvolver='mtmfs',nterms=3,interactive=False,scales=[0,10,30],restart=True,calcpsf=False,calcres=False,savemodel='modelcolumn')

uvsub(vis=myvis)


### TBAVG AND MAKE TSERIES ###

from dynspec.tbavg import tbavg

tbvis = 'YZCet_small_all.ScAp0.tbavg.ms'
tbavg(myvis,tbvis,datacolumn='corrected')

plotms(vis=tbvis,xaxis='scan',yaxis='real',ydatacolumn='data',avgtime='1e8',avgchannel='256',avgspw=True,correlation='RR,LL',coloraxis='corr',plotfile='plots/tseries/YZCet.ScAp0.png')































### IMAGE TIMES WITH EXCESS FLUX IN TIME SERIES ###

# 1) scan 332 (and 333-335 to lesser extent) - looks like noise but good to check - L pol only

tclean(vis=vis1spw,imagename='images/nobg_IQUV_scp0_scan332/nobg_IQUV_scp0_scan332',datacolumn='corrected',cell='3arcsec',imsize=1024,stokes='IQUV',niter=0,pblimit=-0.001,interactive=False,scan='332')

# 2) scan 170~190: R pol excess

tclean(vis=vis1spw,imagename='images/nobg_IQUV_scp0_scan170_190/nobg_IQUV_scp0_scan170_190',datacolumn='corrected',cell='3arcsec',imsize=1024,stokes='IQUV',niter=0,pblimit=-0.001,interactive=False,scan='170~190')
# I = 410 uJy, V = 320 uJy

# 3) scan 246~264: excess in both pols
tclean(vis=vis1spw,imagename='images/nobg_IQUV_scp0_scan246_264/nobg_IQUV_scp0_scan246_264',datacolumn='corrected',cell='3arcsec',imsize=1024,stokes='IQUV',niter=0,pblimit=-0.001,interactive=False,scan='246~264')

# also all of day 3 (~266-360) has excess in both pols --> brighter quiescent emission on that day

# 4) scan 155-159: possible RR excess at low freqs, kinda sketchy looking
# based on looking at the spectrum it just looks like RFI in some channels

