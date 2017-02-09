"""
test_plot.py: Script to test functionality of the Dynspec class as I'm
              modifying it to use complex visibilities.
"""

import subprocess
from dynspec.plot import *

# tested spline interpolation to smooth tseries
# From time series z (a masked array):
#    t = s.get_tlist()
#    ind = find(~s.mask)
#    ti = t[ind]
#    zi = z[ind]
#    us = scipy.interpolate.UnivariateSpline(ti,zi,s=0.5) #0.48 also works
#    z_smooth=us(ti)
# This gives us a pretty good smoothed version of the time series

def load_filelist():
    # returns a dictionary whose keys are file directories for a single obs and whose
    #  entries are lists of dynspec files (e.g. one for L band, one for S band)
    filelist = {}
    cmd = 'ls -d /data/jrv/*/*/*/*/*.tbavg.ms.dynspec | cut -d / -f 1-6 | uniq'
    obslist = subprocess.check_output(cmd, shell=True).rstrip().split('\n') # get list of every obs that has dynspec files
    for obs in obslist:
        cmd = 'ls -d ' + obs + '/*/*.tbavg.ms.dynspec'
        flist = subprocess.check_output(cmd,shell=True).rstrip().split('\n')
        filelist[obs] = flist
    return filelist

def get_srcname(filename):
    # get source name from dynspec filename
    return filename.split('/')[4]

def get_band(filename):
    # get band from dynspec filename
    return filename.split('/')[6]

nt = 150
nf = 64

savefile = '/data/jrv/burst_paper/all_dynspec.npy'

filelist = load_filelist()

# for testing:
#filelist={'/data/jrv/15A-416/YZCMi/1': ['/data/jrv/15A-416/YZCMi/1/L/YZCMi_1L.tbavg.ms.dynspec', '/data/jrv/15A-416/YZCMi/1/S/YZCMi_1S.tbavg.ms.dynspec']}
#keys = filelist.keys()[2:4]
#filelist = {k:filelist[k] for k in keys}
#print filelist
#filelist = {'/data/jrv/15A-416/UVCet/5':['/data/jrv/15A-416/UVCet/5/L/UVCet_5L.tbavg.ms.dynspec','/data/jrv/15A-416/UVCet/5/S/UVCet_5S.tbavg.ms.dynspec']}

ds_list = {}
for obs in filelist:
    src = get_srcname(obs)
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
    ds_obs = ds_obs.bin_dynspec(nt=nt,nf=nf)
    ds_list[obs] = ds_obs

save(savefile,ds_list)
print '\nSaved dictionary of all dynspec to', savefile

