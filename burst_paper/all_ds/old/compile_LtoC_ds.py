"""
compile_ds.py: Script to load ds files for each epoch (L, S, C band only) and combine them, and flag RFI and average to a desired resolution before saving.

RUN IN IPYTHON NOT CASA, OTHERWISE THE DICTIONARY WON'T SAVE PROPERLY AND YOU WON'T BE ABLE TO LOAD IT IN CASA OR IPYTHON :(
"""

from dynspec.plot import *
from dynspec.pipeline_utils import load_burst_filelist
import numpy as np
import pickle

def get_band(filename):
    # get band from dynspec filename
    return filename.split('/')[6]

nt = 150
nf = 64

savefile = '/data/jrv/burst_paper/all_burst_dynspec.npy'

filelist = load_burst_filelist()

ds_list = {}
for obs in filelist:
    ds_files = filelist[obs]
    ds_obs = None
    for f in ds_files:
        band = get_band(f)
        if band == 'X' or band == 'P':
            continue
        params={'filename':f,'uniform':True}
        ds = Dynspec(params)
        ds.spec['i'] = (ds.spec['rr']+ds.spec['ll'])/2
        ds.spec['v'] = (ds.spec['rr']-ds.spec['ll'])/2
        del ds.spec['rr']
        del ds.spec['ll']
        if ds_obs is None:
            ds_obs = deepcopy(ds)
        else:
            ds_obs.add_dynspec(ds)
        del ds
    ds_obs.mask_RFI(rmsfac=3.) # so far have only plotted YZ CMi dynspec w/ this - only makes a minor diff - maybe bigger when there is less binning?
    ds_obs = ds_obs.bin_dynspec(nt=nt,nf=nf)
    ds_list[obs] = ds_obs

np.save(savefile,ds_list)
print '\nSaved dictionary of all dynspec to', savefile

'''
output = open(savefile[:-3] + 'pickle','wb')
# b is for binary
# w opens and overwrites existing file
pickle.dump(ds_list,output)
output.close()
'''
