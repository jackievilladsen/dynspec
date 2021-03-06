'''
image_vlba_intervals.py

Goal: Image location of star in few-scan intervals throughout whole dataset, comparing different calibration and uvrange filters
  - use ms files created by make_vlba_dynspec.py
'''

# user defines what obs to run this on
# maybe swap this out for a list of obs
obsnum = '3'
star = 'ADLeo'

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

msroot = obsdir + star
ms_star = msroot + '.ms'  
ms_postcal = msroot + '.postcal.ms' # arleady exists
ms_postcal_uvrange = msroot + '.uvrange.ms' # already exists
ms_precal = msroot + '.precal.ms' # need to create - same as ms_star data column, no uvrange
ms_precal_uvrange = msroot + '.precal.uvrange.ms' # need to create but limited to antennas specified by anten
ms_precal_uvrange = msroot + 'precal.uvrange.ms' # need to create but limited to antennas specified by anten

uvrange=uvrange_dict.get(obsname,'<2e6m')
anten=anten_dict.get(obsname,'')

from dynspec.tbavg import tbavg,dyn_spec,scan_reindex
from dynspec.extract_dynspec import saveTxt
import os
from numpy import arange
'''
pcen_dict = {'UVCet_3':'J2000 01h39m05.1153408 -17d56m51.7770071',
	     'UVCet_4':'J2000 01h39m05.1246492 -17d56m51.7916119',
	     'UVCet_5':'J2000 01h39m05.1472030 -17d56m51.8803255',
	     'ADLeo_3':'J2000 10h19m35.7223057 19d52m11.3708328',
	     'ADLeo_4':'J2000 10h19m35.7232319 19d52m11.3517484',
	     'ADLeo_5':'J2000 10h19m35.7288197 19d52m11.2845273'}
pcen = pcen_dict[obsname]

# ms_star: data column is uncalibrated, corrected column is calibrated, no anten or uvrange limits applied - phasecen has not been shifted
# ms_postcal: calibrated star data with no anten or uvrange limits applied - phasecen has been shifted
# ms_postcal_uvrange: calibrated star data with anten and uvrange limits applied - phasecen has been shifted

# create uncalibrated versions of ms w/ fixvis applied and w/ and w/o uvrange/anten applied
fixvis(vis=ms_star,phasecenter=pcen)
if os.path.exists(ms_precal):
	os.system('rm -rf '+ms_precal)
if os.path.exists(ms_precal_uvrange):
	os.system('rm -rf '+ms_precal_uvrange)
split(vis=ms_star,datacolumn='data',outputvis=ms_precal) # uncalibrated w/o uvrange or anten applied
split(vis=ms_star,datacolumn='data',outputvis=ms_precal_uvrange,antenna=anten,uvrange=uvrange)
'''
# make directories for images - one directory per ms version
imdir = obsdir+'images/'
if not os.path.exists(imdir):
	os.system('mkdir '+imdir)
dir_postcal = imdir + 'postcal/'
dir_postcal_uvrange = imdir + 'postcal_uvrange/'
dir_precal = imdir + 'precal/'
dir_precal_uvrange = imdir + 'precal_uvrange/'
for dirname in [dir_postcal,dir_postcal_uvrange,dir_precal,dir_precal_uvrange]:
	if not os.path.exists(dirname):
		os.system('mkdir '+dirname)
dir_dict = {ms_postcal:dir_postcal,ms_postcal_uvrange:dir_postcal_uvrange,ms_precal:dir_precal,ms_precal_uvrange:dir_precal_uvrange}

# plot dirty images for every 10 scans --> ~30 minute intervals
# total number of scans: ADLeo 3 and 4: 75 (38+37); ADLeo 5: 73 (37+36)
Nscan=10
'''
for msfile in dir_dict:
	dirname = dir_dict[msfile]
	#for nf in arange(Nscan,75+Nscan,Nscan):
        for nf in [80]:
		n0 = nf-Nscan+1
		scan = str(n0)+'~'+str(nf)
		name_dirty = dirname+'dirty_scan'+str(n0)+'_'+str(nf)
		os.system('rm -rf '+name_dirty+'.*')
		clean(vis=msfile,imagename=name_dirty,cell='0.00028arcsec',niter=0,stokes='IV',scan=scan)
# source clearly moves around in both ms_precal and ms_postcal and two peaks during flare, motion is identical
# source also clearly moves around in ms_postcal_uvrange and motion similar but not identical


# clean one 10-scan image to determine appropriate clean threshold and generate mask (use CASA viewer to make sure mask works on all dirty images above)
clean(vis=ms_postcal_uvrange,imagename=imdir+star+'_interactive',interactive=True,npercycle=10,cyclefactor=5,cell='0.00028arcsec',stokes='IV')
mask = imdir + star + '_interactive.mask'
clean(vis=ms_postcal_uvrange,imagename=imdir+star+'_interactive',interactive=True,mask=mask,npercycle=10,cyclefactor=5,cell='0.00028arcsec',stokes='IV',scan='51~60') # clean 10 scans to verify mask and see what kind of threshold to aim for

mask = imdir + star + '_interactive.mask'
threshold='0.2mJy'

# remove ms_precal_uvrange from dir_dict since I'm fairly satisfied that I'm going to go with postcal
dir_dict = {ms_postcal:dir_postcal,ms_postcal_uvrange:dir_postcal_uvrange,ms_precal:dir_precal}
for msfile in dir_dict:
	dirname = dir_dict[msfile]
	for nf in arange(Nscan,75+Nscan,Nscan):
		n0 = nf-Nscan+1
		scan = str(n0)+'~'+str(nf)
		name_clean = dirname+'clean_scan'+str(n0)+'_'+str(nf)
		os.system('rm -rf '+name_clean+'.*')
		clean(vis=msfile,imagename=name_clean,cell='0.00028arcsec',threshold=threshold,cyclefactor=5,mask=mask,stokes='IV',scan=scan)	
# conclusions: I'm happy with ms_postcal, I think it is real
# also there is a flare in scans 11-20 (AD Leo 3)
'''
# image full time duration
imname = dir_postcal + 'clean_full'
clean(vis=ms_postcal,imagename=imname,cell='0.00028arcsec',interactive=True,npercycle=10,cyclefactor=5,mask=mask,stokes='IV')


