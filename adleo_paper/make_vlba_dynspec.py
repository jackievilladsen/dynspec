'''
make_vlba_dynspec.py

Starting with an uvfits file of VLBA data:
  - Load ms files for pcal and star and split into scans
  - Run amp+phase gaincal on pcal and apply to star and pcal
  - Split out calibrated star, filtering by uvrange and/or antenna as specified
  - Dirty image of star
  - Shift phase center to location identified by user (in pcen_dict)
  - 
'''

# user defines what obs to run this on
# maybe swap this out for a list of obs
obsnum = '5'
star = 'UVCet'
do_image = True

anten_dict = {'UVCet_4':'!LA',
	      'ADLeo_4':'!LA'}
uvrange_dict = {'UVCet_3':'<2e6m',
		'UVCet_4':'<2e6m',
		'UVCet_5':'<2e6m',
		'ADLeo_3':'<2e6m',
		'ADLeo_4':'<4e6m',
		'ADLeo_5':'<2e6m'}

stardir = '/data/jrv/BV071/'+star+'/'
obsdir = stardir+obsnum+'/'
obsname = star + '_' + obsnum

pcal_dict = {'UVCet':'J0132','ADLeo':'J1024'}
pcal = pcal_dict[star]
pcal_model = stardir+pcal+'.model'
uvfits_pcal = obsdir + pcal + '.uvfits'
ms_pcal = obsdir + pcal + '.ms'
caltab = obsdir + pcal + '.AP0'

msroot = obsdir + star
uvfits_star = msroot + '_shift.uvfits' # _shift.ms: version that has been shifted to account for parallax+proper+orbital motion
ms_star = msroot + '.ms'  
ms_postcal = msroot + '.postcal.ms'
ms_uvrange = msroot + '.uvrange.ms'
star_image = msroot + '.image'
ms_tbavg = msroot +'_'+obsnum+ 'X.tbavg.ms'

uvrange=uvrange_dict.get(obsname,'<2e6m')
anten=anten_dict.get(obsname,'')

from dynspec.tbavg import tbavg,dyn_spec,scan_reindex
from dynspec.extract_dynspec import saveTxt
import os
'''
### Load uvfits and split into scans ###

if os.path.exists(ms_pcal):
	os.system('rm -rf '+ms_pcal)
importuvfits(fitsfile=uvfits_pcal,vis=ms_pcal)
scan_reindex(ms_pcal)

if os.path.exists(ms_star):
	os.system('rm -rf '+ms_star)
importuvfits(fitsfile=uvfits_star,vis=ms_star)
scan_reindex(ms_star)


### Amp+phase gaincal on pcal ###

# load pcal model
ft(vis=ms_pcal,model=pcal_model)

if os.path.exists(caltab):
	os.system('rm -rf '+caltab)
gaincal(vis=ms_pcal,refant='FD',calmode='ap',caltable=caltab,solnorm=True)

plotcal(caltable=caltab, xaxis='time',yaxis='phase',iteration='antenna',subplot=431,showgui=False,figfile=obsdir+'AP0_phase.png')
plotcal(caltable=caltab, xaxis='time',yaxis='amp',iteration='antenna',subplot=431,showgui=False,figfile=obsdir+'AP0_amp.png')

if obsname == 'ADLeo_5':
	flagdata(vis=caltab,mode='clip',antenna='HN,NL,SC',clipminmax=[0.75,2],datacolumn='CPARAM') # clip caltable gains below 2
	plotcal(caltable=caltab, xaxis='time',yaxis='phase',iteration='antenna',subplot=431,showgui=False,figfile=obsdir+'AP0_phase_flagged.png')
	plotcal(caltable=caltab, xaxis='time',yaxis='amp',iteration='antenna',subplot=431,showgui=False,figfile=obsdir+'AP0_amp_flagged.png')

# apply to star and split out corrected data (w/ and w/o uvrange)
applycal(vis=ms_star,gaintable=caltab)
if os.path.exists(ms_postcal):
	os.system('rm -rf '+ms_postcal)
split(vis=ms_star,outputvis=ms_postcal,datacolumn='corrected')

if do_image:
	allBL_image = msroot+'.allBL'
	os.system('rm -rf '+allBL_image+'.*')
	clean(vis=ms_postcal,imagename=allBL_image,cell='0.00028arcsec',niter=5000,cyclefactor=5,interactive=True,npercycle=10,stokes='IV',antenna=anten)

if os.path.exists(ms_uvrange):
	os.system('rm -rf '+ms_uvrange)
split(vis=ms_star,outputvis=ms_uvrange,datacolumn='corrected',uvrange=uvrange,antenna=anten)


### Create image of star so user can identify phase center then shift phase center ###

if do_image:
	os.system('rm -rf '+ star_image + '.*')
	clean(vis=ms_uvrange,imagename=star_image,cell='0.00028arcsec',niter=5000,cyclefactor=5,interactive=True,npercycle=10,stokes='IV')
# AD Leo 5: there may be a bad antenna in here, it looks awful! scans 1-40 have best calibration

if obsname=='ADLeo_5':
	clean(vis=ms_uvrange,imagename='5/scan1_40',cell='0.00028arcsec',niter=5000,cyclefactor=5,interactive=True,npercycle=10,stokes='IV',scan='1~40')
# still looks bad, but better

pcen_dict = {'UVCet_3':'J2000 01h39m05.1153408 -17d56m51.7770071',
	     'UVCet_4':'J2000 01h39m05.1246492 -17d56m51.7916119',
	     'UVCet_5':'J2000 01h39m05.1472030 -17d56m51.8803255',
	     'ADLeo_3':'J2000 10h19m35.7223057 19d52m11.3708328',
	     'ADLeo_4':'J2000 10h19m35.7232319 19d52m11.3517484',
	     'ADLeo_5':'J2000 10h19m35.7288197 19d52m11.2845273'}
pcen = pcen_dict[obsname]

fixvis(vis=ms_postcal,phasecenter=pcen)
fixvis(vis=ms_uvrange,phasecenter=pcen)

star_im_dirty = msroot+'.dirty'
os.system('rm -rf '+ star_im_dirty + '.*')
clean(vis=ms_uvrange,imagename=star_im_dirty,cell='0.00028arcsec',niter=0,stokes='IV')
print 'Check',star_im_dirty,'to verify star at center.'

### Average over all baselines and export dynamic spectrum ###

if os.path.exists(ms_tbavg):
	os.system('rm -rf '+ms_tbavg)

tbavg(split, ms_uvrange, ms_tbavg)

if os.path.exists(ms_tbavg+'.dynspec'):
	os.system('rm -rf '+ms_tbavg+'.dynspec')
spec = dyn_spec(ms_tbavg)
saveTxt(spec,ms_tbavg)
'''
savedir = '/data/jrv/15A-416/'+star+'/'+obsnum+'/X'
if not os.path.exists(savedir):
	os.system('mkdir '+savedir)
os.system('cp -r '+ms_tbavg+'.dynspec '+savedir)
