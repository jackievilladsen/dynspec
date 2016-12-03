'''
UVCet_compareAB.py

Purpose: use modelA and modelB to calibrate UV Cet and compare
structure/motion of UVCet over time with the two calibrations
'''

UVCet_3_uvfits = '/data/jrv/15A-416/UVCet/3/X/UVCet_shift.uvfits'
J0132_3_uvfits = '/data/jrv/15A-416/UVCet/3/X/J0132.uvfits'

modelA = '/data/jrv/BV071/UVCet/J0132_refJ0140.model'
modelB = '/data/jrv/BV071/UVCet/J0132_self.model'

working_dir = '/data/jrv/BV071/UVCet/testAB/'
import os
os.chdir(working_dir)

starms = 'UVCet_shift.ms'
calms = 'J0132.ms'
importuvfits(fitsfile=UVCet_3_uvfits,vis=starms)
importuvfits(fitsfile=J0132_3_uvfits,vis=calms)
from dynspec.tbavg import scan_reindex
scan_reindex(starms)
scan_reindex(calms)

ft(vis=calms,model=modelA)
gaincal(vis=calms,refant='FD',calmode='ap',solnorm=True,caltable='A.ap0')
ft(vis=calms,model=modelB)
gaincal(vis=calms,refant='FD',calmode='ap',solnorm=True,caltable='B.ap0')

plotcal(caltable='A.ap0', xaxis='time',yaxis='phase',iteration='antenna',subplot=431)
plotcal(caltable='B.ap0', xaxis='time',yaxis='phase',iteration='antenna',subplot=431)
# nothing scary looking, amps look same, phases different b/c of different location of J0132 in the 2 models

applycal(vis=starms,gaintable='A.ap0')
for nf in range(8,73,8):
    n0 = nf-7
    scan = str(n0)+'~'+str(nf)
    name = 'scan'+str(n0)+'_'+str(nf)
    clean(vis=starms,imagename='images/A_dirty_'+name,cell='0.00028arcsec',niter=0,stokes='IV',scan=scan)

applycal(vis=starms,gaintable='B.ap0')
for nf in range(8,73,8):
    n0 = nf-7
    scan = str(n0)+'~'+str(nf)
    name = 'scan'+str(n0)+'_'+str(nf)
    clean(vis=starms,imagename='images/B_dirty_'+name,cell='0.00028arcsec',niter=0,stokes='IV',scan=scan)

# motion of star between scans seems to be the same for both calibrations

# now do clean image of scan 37 (burst peak)
clean(vis=starms,imagename='images/B_clean_37'+name,cell='0.00028arcsec',niter=5000,stokes='IV',scan='37',interactive=True,npercycle=10,cyclefactor=5)
# yuck there is a 2-mJy copy of the burst 0.9 mas to the southwest...
# something in our calibration must be hideously wrong
# maybe bad interpolation of rates?

applycal(vis=starms,gaintable='A.ap0')
clean(vis=starms,imagename='images/A_clean_37'+name,cell='0.00028arcsec',niter=5000,stokes='IV',scan='37',interactive=True,npercycle=10,cyclefactor=5)
clean(vis=starms,imagename='images/A_clean_scan1_8'+name,cell='0.00028arcsec',niter=5000,stokes='IV',scan='1~8',interactive=True,npercycle=10,cyclefactor=5)
# equally meh

# inspect J0132 data
applycal(vis=calms,gaintable='B.ap0')

# try uvmodelfit versus time
# if we want SNR of 10, and sensitivity of one scan is ~200 uJy/beam according to calculator --> 4 scans should reach 100 uJy/beam

statwt(starms)
uvmodelfit(starms,scan='37',outfile='test.cl')
uvmodelfit(starms,scan='1~8',outfile='test1_8P.cl')

# try imaging w/o certain ants and see if image improves
for anten in ['BR','FD','HN','KP','LA','MK','NL','OV','PT','SC']:
    imagename = 'images/scan1_8/no'+anten+'_dirty_A'
    clean(vis=starms,imagename=imagename,cell='0.00028arcsec',niter=0,stokes='IV',antenna='!'+anten,scan='1~8')

clean(vis=starms,imagename='images/scan1_8/noSC_clean_A'+name,cell='0.00028arcsec',niter=5000,stokes='IV',scan='1~8',interactive=True,npercycle=10,cyclefactor=5,antenna='!SC')

starms_nocal = 'uvcet_nocal.ms'
split(vis=starms,outputvis=starms_nocal,datacolumn='data')
clean(vis=starms_nocal,imagename='images/scan1_8/noSC_clean_nocal',cell='0.00028arcsec',niter=5000,stokes='IV',scan='1~8',interactive=True,npercycle=10,cyclefactor=5,antenna='!SC')

# try flagging first integration of each scan
flagdata(vis=starms,mode='quack',quackinterval=1.0)

clean(vis=starms,imagename='images/scan1_8/spw0_3_dirty',cell='0.00028arcsec',niter=0,stokes='IV',scan='1~8',spw='0~3')
clean(vis=starms,imagename='images/scan1_8/spw4_7_dirty',cell='0.00028arcsec',niter=0,stokes='IV',scan='1~8',spw='4~7')
# burst and quiescent separation: 1.3 mas --> 5.6e6m

# nothing gets rid of the big uncleanable sidelobes!  worst dynamic range ever
