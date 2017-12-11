"""
compile_Pds.py: compile Dynspec objects for all P band dynamic spectra, binned to low resolution

"""

import dynspec.pipeline_utils
reload(dynspec.pipeline_utils)

import subprocess
from dynspec.plot import *
from dynspec.pipeline_utils import load_burst_filelist

def get_band(filename):
    # get band from dynspec filename
    return filename.split('/')[6]

nt = 50
nf = 32

savefile = '/data/jrv/burst_paper/ds/all_burst_Pdynspec.npy'

filelist = load_burst_filelist(band='P')

ds_list = {}
for obs in filelist:
    ds_files = filelist[obs]
    for f in ds_files:
        band = get_band(f)
        if band != 'P':
            continue
        params={'filename':f,'uniform':True}
        ds = Dynspec(params)
        ds.spec['i'] = (ds.spec['xx']+ds.spec['yy'])/2
        ds.spec['v'] = (ds.spec['xy']-ds.spec['yx'])/(2.j)
        del ds.spec['xx']
        del ds.spec['yy']
        del ds.spec['xy']
        del ds.spec['yx']
    ds.mask_RFI(rmsfac=1.5)
    ds = ds.bin_dynspec(nt=nt,nf=nf)
    ds_list[obs] = ds

save(savefile,ds_list)
print '\nSaved dictionary of all dynspec to', savefile

