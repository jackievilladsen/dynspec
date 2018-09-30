"""
compile_ds_with_Pband.py: Script to load ds files for each epoch (P, L, S band) and combine them, and flag RFI and average to a desired resolution before saving.

RUN IN IPYTHON NOT CASA, OTHERWISE THE DICTIONARY WON'T SAVE PROPERLY AND YOU WON'T BE ABLE TO LOAD IT IN CASA OR IPYTHON :(
"""

import dynspec.plot
reload(dynspec.plot)

from dynspec.plot import *
from dynspec import load_dict
from dynspec.pipeline_utils import load_burst_filelist
import numpy as np
import pickle

def get_band(filename):
    # get band from dynspec filename
    return filename.split('/')[6]

savefile = '/data/jrv/burst_paper/all_burst_epoch_dynspec_LSband.npy'
savefileP = '/data/jrv/burst_paper/all_burst_epoch_dynspec_Pband.npy'

filelist = load_burst_filelist()

try:
    ds_list = load_dict(savefile)
except:
    ds_list = {}
try:
    dsP_list = load_dict(savefileP)
except:
    dsP_list = {}

reload_obslist = filelist.keys()
#reload_obslist = ['/data/jrv/15A-416/EQPeg/2'] # change this to only reload one observation
for obs in reload_obslist:
    ds_files = filelist[obs]
    ds_obs = None
    band_list = [get_band(f) for f in ds_files]
    for f in ds_files:
        band = get_band(f)
        if band == 'X':
            continue
        
        # load dynamic spectrum
        params={'filename':f,'uniform':True,'convert_stokes':True}
        ds_band = Dynspec(params)
        
        # mask pixels and channels with strong RFI
        ds_band.mask_RFI_pixels(rmsfac=10.,func=imag)
        ds_band.mask_RFI(rmsfac=7.)
        
        # bin bands to have matching time-freq resolution
        #   - otherwise add_dynspec inserts blank channels/timestamps, making file large
        df_bin_MHz = 2                       # desired frequency resolution (this is the S band channel width)
        dt_bin = 6                           # desired time resolution (some of my P band ds's have this res)
        nt = int(round(dt_bin/ds_band.dt())) # number of integrations to bin together
        nf = int(round(df_bin_MHz/(ds_band.df()/1e6))) # number of channels to bin together
        if band=='S' and ds_band.df()==1.e6:
            mask_partial = 0.75   # if S band has "1-MHz" channels, half are blank/flagged
        else:
            mask_partial = 0.5
        if band == 'P':
            # don't bin in frequency for P band (so I can play with RFI flagging later)
            ds_band = ds_band.bin_dynspec(nt=nt,nf=1,mask_partial=mask_partial)
        else:
            ds_band = ds_band.bin_dynspec(nt=nt,nf=nf,mask_partial=mask_partial)
        
        # add bands together
        if ds_obs is None:
            ds_obs = deepcopy(ds_band)
        else:
            ds_obs.add_dynspec(ds_band)
        del ds_band
    
    # split apart L/S and P dynspecs if needed (they were first merged together to line up time axis)
    #   and save them to the dictionaries of bursts
    if 'P' in band_list:
        ds = ds_obs.clip(fmin=0.999e9,fmax=4.01e9,trim_mask=False)
        ds = ds.bin_dynspec(nt=1,nf=16)
        del ds.spec['q']
        del ds.spec['u']
        ds_list[obs] = ds
        dsP_list[obs] = ds_obs.clip(fmin=2.2e8,fmax=5.e8,trim_mask=False)
    else:
        ds_list[obs] = ds_obs

np.save(savefile,ds_list)
print '\nSaved dictionary of all L,S band dynspec to', savefile

np.save(savefileP,dsP_list)
print '\nSaved dictionary of all P band dynspec to', savefileP
