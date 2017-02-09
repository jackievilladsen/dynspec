'''
compile_ds_VLAonly.py - Load all bands for a given epoch and calculate and save a Stokes V dynamic spectrum combining all the bands.
    Adjusted for VLA-only observations (doesn't try to load VLBA data).
'''

from copy import deepcopy
from dynspec.plot import Dynspec
from pylab import *
import os
import pickle

src = 'UVCet'
obs_dict = {'2013':['1'],'2015':['2']}
obs_dict = {'2015':['2']}
proj_dict = {'2013':'13A-423','2015':'15A-416'}
bandlist = ['L','S','C','P']

params = {'legend.fontsize': 'small',
          'axes.titlesize': 'medium',
          'axes.labelsize': 'small',
          'xtick.labelsize': 'x-small',
          'ytick.labelsize': 'x-small',
          'image.interpolation': 'hanning'}
rcParams.update(params)

savedir = '/data/jrv/'+src.lower()+'_paper/'
if not os.path.exists(savedir):
    os.system('mkdir '+savedir)

# Cycle through years
for year in obs_dict:
    epochlist = obs_dict[year]
    proj = proj_dict[year]
    # Cycle through epochs
    for epoch in epochlist:
        # Cycle through bands and load dynspec for each, adding to single large dynspec
        ds = None
        for band in bandlist:
            ds_file = '/data/jrv/'+proj+'/'+src+'/'+epoch+'/'+band+'/'+src+'_'+epoch+band+'.tbavg.ms.dynspec'
            print ds_file
            if not os.path.exists(ds_file):
                print 'File does not exist, skipping this band:',ds_file
                continue
            params={'filename':ds_file,'uniform':True}

            # load dynspec and calculate stokes V dynspec, delete other pols
            ds_band = Dynspec(params)
            if band == 'P':
                specV = (ds_band.spec['xy']-ds_band.spec['yx'])/(2.j)
                specI = (ds_band.spec['xx']+ds_band.spec['yy'])/2
            else:
                specV = (ds_band.spec['rr']-ds_band.spec['ll'])/2
                specI = (ds_band.spec['rr']+ds_band.spec['ll'])/2
            del ds_band.spec
            ds_band.spec = {'v':specV,'i':specI}
            del specV,specI

            # bin bands to have matching time-freq resolution - otherwise add_dynspec inserts blanks which makes file large
            if band=='P': # bin P band to 2-MHz resolution to match S band
                ds_band.mask_RFI(rmsfac=1.5)
                ds_band = ds_band.bin_dynspec(nt=1,nf=16,mask_partial=0.5)
            elif band=='S': # bin S band to 2-MHz resolution, 6-sec integrations to match current version of P band data reduction
                ds_band = ds_band.bin_dynspec(nt=6,nf=2,mask_partial=0.75) # 50% masked already in S band b/c every other 1-MHz channel blank
            elif band=='L': # bin L band to 2-MHz resolution, 6-sec integrations to match current version of P band data reduction
                #ds_band.mask_RFI(rmsfac=3.)
                ds_band = ds_band.bin_dynspec(nt=6,nf=2,mask_partial=0.5)
            elif band=='C': # bin C band to 6-sec integations (should already have 2-MHz resolution)
                ds_band = ds_band.bin_dynspec(nt=6,nf=1,mask_partial=0.5)

            # add bands together
            if ds is None:
                ds = deepcopy(ds_band)
            else:
                ds.add_dynspec(ds_band)
            del ds_band
            
        # save multi-band dynspec
        savefile = savedir + year + '_' + src + '_' + epoch + '.dynspec.pickle'
        print 'Saving multi-band dynspec to',savefile
        pickle.dump(ds, open(savefile,"wb"))
