"""
compile_Pds.py: compile Dynspec objects for all P band dynamic spectra, binned to low resolution

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

def load_P_filelist():
    # returns a dictionary whose keys are file directories for a single P band obs and whose
    #  entries are lists of dynspec files
    filelist = {}
    cmd = 'ls -d /data/jrv/15A-416/*/*/P/*.tbavg.ms.dynspec | cut -d / -f 1-6 | uniq'
    obslist = subprocess.check_output(cmd, shell=True).rstrip().split('\n') # get list of every obs that has dynspec files
    for obs in obslist:
        cmd = 'ls -d ' + obs + '/P/*.tbavg.ms.dynspec'
        flist = subprocess.check_output(cmd,shell=True).rstrip().split('\n')
        filelist[obs] = flist
    return filelist

def get_srcname(filename):
    # get source name from dynspec filename
    return filename.split('/')[4]

def get_band(filename):
    # get band from dynspec filename
    return filename.split('/')[6]

nt = 50
nf = 32

savefile = '/data/jrv/burst_paper/ds/all_Pdynspec.npy'

filelist = load_P_filelist()

# for testing:
#filelist={'/data/jrv/15A-416/YZCMi/1': ['/data/jrv/15A-416/YZCMi/1/L/YZCMi_1L.tbavg.ms.dynspec', '/data/jrv/15A-416/YZCMi/1/S/YZCMi_1S.tbavg.ms.dynspec']}
#keys = filelist.keys()[2:4]
#filelist = {k:filelist[k] for k in keys}
#print filelist
#filelist = {'/data/jrv/15A-416/UVCet/4': ['/data/jrv/15A-416/UVCet/4/P/UVCet_4P.tbavg.ms.dynspec'], '/data/jrv/15A-416/UVCet/5': ['/data/jrv/15A-416/UVCet/5/P/UVCet_5P.tbavg.ms.dynspec']}
#filelist = {'/data/jrv/15A-416/UVCet/5':['/data/jrv/15A-416/UVCet/5/L/UVCet_5L.tbavg.ms.dynspec','/data/jrv/15A-416/UVCet/5/S/UVCet_5S.tbavg.ms.dynspec']}

ds_list = {}
for obs in filelist:
    src = get_srcname(obs)
    ds_files = filelist[obs]
    for f in ds_files:
        band = get_band(f)
        params={'filename':f,'uniform':True}
        ds = Dynspec(params)
        ds.spec['i'] = (ds.spec['xx']+ds.spec['yy'])/2
        ds.spec['v'] = (ds.spec['xy']-ds.spec['yx'])/(2.j)
        del ds.spec['xx']
        del ds.spec['yy']
        del ds.spec['xy']
        del ds.spec['yx']
    ds = ds.bin_dynspec(nt=nt,nf=nf)
    ds_list[obs] = ds

save(savefile,ds_list)
print '\nSaved dictionary of all dynspec to', savefile

