'''
compile_adleo_ds.py - Load all bands for a given epoch and calculate and save a Stokes V dynamic spectrum combining all the bands.
    Cycle through and do this for all 3 epochs of AD Leo VLA+VLBA observations.
'''

from copy import deepcopy
from dynspec.plot import Dynspec
from pylab import *
import os
import pickle

src = 'ADLeo'
epochlist = ['3','4','5']
epochlist = ['3']
bandlist = ['L','S','P','X']

params = {'legend.fontsize': 'small',
          'axes.titlesize': 'medium',
          'axes.labelsize': 'small',
          'xtick.labelsize': 'x-small',
          'ytick.labelsize': 'x-small',
          'image.interpolation': 'hanning'}
rcParams.update(params)

savedir = '/data/jrv/adleo_paper/'
if not os.path.exists(savedir):
    os.system('mkdir '+savedir)

# Cycle through epochs
for epoch in epochlist:

    # Cycle through bands and load dynspec for each, adding to single large dynspec
    ds = None
    for band in bandlist:
        if band=='X':
            ds_file = '/data/jrv/15A-416/'+src+'/'+epoch+'/X/'+src+'.uv2e6.tbavg.ms.dynspec'
        else:
            ds_file = '/data/jrv/15A-416/'+src+'/'+epoch+'/'+band+'/'+src+'_'+epoch+band+'.tbavg.ms.dynspec'
        print ds_file    
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
            ds_band = ds_band.bin_dynspec(nt=1,nf=16,mask_partial=0.5)
        elif band=='X': # bin VLBA to 6-sec resolution to match P band (has 8-MHz resolution)
            ds_band = ds_band.bin_dynspec(nt=3,nf=1,mask_partial=0.5)
        elif band=='S': # bin S band to 2-MHz resolution, 6-sec integrations to match current version of P band data reduction
            ds_band = ds_band.bin_dynspec(nt=6,nf=2,mask_partial=0.75) # 50% masked already in S band b/c every other 1-MHz channel blank
        elif band=='L': # bin L band to 2-MHz resolution, 6-sec integrations to match current version of P band data reduction
            ds_band = ds_band.bin_dynspec(nt=6,nf=2,mask_partial=0.5)
        # add bands together
        if ds is None:
            ds = deepcopy(ds_band)
        else:
            ds.add_dynspec(ds_band)
        del ds_band
    
    # break apart into VLA and VLBA dynspec (merged them together to get them on the same time axis)
    print 'Separating VLA and VLBA dynspecs (merged them first so they would be on the same time axis)'
    dsVLA = ds.clip(fmax=4.0e9,trim_mask=False)
    dsVLBA = ds.clip(fmin=8.28e9,trim_mask=False)
    dsVLBA = dsVLBA.bin_dynspec(nt=1,nf=4) # return to 8-MHz resolution, getting rid of blank channels
    ds_dict = {'VLA':dsVLA,'VLBA':dsVLBA}
    
    # save multi-band dynspec
    savefile = savedir + src + '_' + epoch + '.dynspec.pickle'
    print 'Saving VLA & VLBA dynspecs to',savefile
    pickle.dump(ds_dict, open(savefile,"wb"))
