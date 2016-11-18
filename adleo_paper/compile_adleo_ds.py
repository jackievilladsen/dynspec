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
bandlist = ['L','S','P']

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
        ds_file = '/data/jrv/15A-416/'+src+'/'+epoch+'/'+band+'/'+src+'_'+epoch+band+'.tbavg.ms.dynspec'
        print ds_file    
        params={'filename':ds_file,'uniform':True}

        # load dynspec and calculate stokes V dynspec, delete other pols
        ds_band = Dynspec(params)
        if band == 'P':
            specV = (ds_band.spec['xy']-ds_band.spec['yx'])/(2.j)
        else:
            specV = (ds_band.spec['rr']-ds_band.spec['ll'])/2
        del ds_band.spec
        ds_band.spec = {'v':specV}
        del specV

        # bin bands to have matching time-freq resolution - otherwise add_dynspec inserts blanks which makes file large
        if band=='P': # bin P band to 2-MHz resolution to match S band
            ds_band = ds_band.bin_dynspec(nt=1,nf=16)
        else: # bin L, S band to 2-MHz resolution, 6-sec integrations to match current version of P band data reduction
            ds_band = ds_band.bin_dynspec(nt=6,nf=2)

        # add bands together
        if ds is None:
            ds = deepcopy(ds_band)
        else:
            ds.add_dynspec(ds_band)
        del ds_band
    
    # save multi-band dynspec
    savefile = savedir + src + '_' + epoch + '.dynspec.pickle'
    pickle.dump(ds, open(savefile,"wb"))
