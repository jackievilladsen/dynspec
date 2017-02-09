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

stardir = '/data/jrv/BV071/'+star+'/'
obsdir = stardir+obsnum+'/'
imdir = obsdir+'images/'
obsname = star + '_' + obsnum

msroot = obsdir + star
ms_star = msroot + '.ms'  
ms_postcal = msroot + '.postcal.ms' # already exists

anten=anten_dict.get(obsname,'')

from dynspec.tbavg import tbavg,dyn_spec,scan_reindex
from dynspec.extract_dynspec import saveTxt
import os
from numpy import arange

dirname = imdir + 'postcal/imflare/'
if not os.path.exists(dirname):
	os.system('mkdir '+dirname)

# clean interactively once to create mask
clean(vis=ms_postcal,imagename=imdir+star+'_Iinteractive',cell='0.00028arcsec',interactive=True,npercycle=10,cyclefactor=5,stokes='I')

mask = imdir + star + '_Iinteractive.mask'
threshold='0.3mJy'
# ADLeo 3: made flare_compare.pdf by plotting scans 41-50 (contour) and 51-60 (color)

# Scan ranges for flares, before, after
preflare1 = '0~25'
preflare2 = '26~49'
flare1 = '50~51'
between = '52~56'
flare2 = '57~58'
postflare = '59~75'

scanlist = [preflare1,preflare2,flare1,between,flare2,postflare]

for scan in scanlist:
	name_I = (dirname+'cleanI_scan'+scan).replace('~','_')
	name_V = (dirname+'cleanV_scan'+scan).replace('~','_')
	os.system('rm -rf '+name_I+'.*')
	os.system('rm -rf '+name_V+'.*')
	clean(vis=ms_postcal,imagename=name_I,cell='0.00028arcsec',threshold=threshold,cyclefactor=5,mask=mask,stokes='I',scan=scan)
	clean(vis=ms_postcal,imagename=name_V,cell='0.00028arcsec',threshold=threshold,cyclefactor=5,mask=mask,stokes='V',scan=scan)	
