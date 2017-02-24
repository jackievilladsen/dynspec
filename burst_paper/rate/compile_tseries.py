"""

compile_tseries.py: Script to compile Stokes I tseries for all P,L,S band data from VLA survey.

"""

import os
from dynspec.plot import *
integration_time = 600
survey = 'Plo'

# tested spline interpolation to smooth tseries
# From time series z (a masked array):
#    t = s.get_tlist()
#    ind = find(~s.mask)
#    ti = t[ind]
#    zi = z[ind]
#    us = scipy.interpolate.UnivariateSpline(ti,zi,s=0.5) #0.48 also works
#    z_smooth=us(ti)
# This gives us a pretty good smoothed version of the time series

def load_filelist(band):
    # makes a file dynspec_filelist.txt containing the locations of all dynspec files for this band
    # then returns a list of strings of all the filenames
    cmd = 'ls -d /data/jrv/*/*/*/'+band+'/*.tbavg.ms.dynspec > dynspec_filelist.txt'
    os.system(cmd)
    f = open('dynspec_filelist.txt')
    lines = f.readlines()
    f.close()
    filelist = [l.strip('\n') for l in lines]
    return filelist

def get_srcname(filename):
    # get source name from dynspec filename
    return filename.split('/')[4]

def get_dist(src):
    # return distance in pc to src
    # (will use this to convert flux to flux at 1 pc)
    dstar = {'ADLeo':4.69,
             'UVCet':2.62,
             'EQPeg':6.18,
             'EVLac':5.12,
             'YZCMi':5.96}
    return dstar[src]

def get_flims(survey):
    # return band,fmin,fmax for survey
    flims_dict = {'VLASS':      {'band':'S', 'fmin':2.e9,   'fmax':4.e9  },   # VLASS is all of VLA S-band
                  'VAST':       {'band':'L', 'fmin':1.13e9, 'fmax':1.43e9},
                  'VASTtest':   {'band':'L', 'fmin':1.13e9, 'fmax':1.43e9},
                  'ThunderKAT': {'band':'L', 'fmin':1.e9,   'fmax':1.75e9},
                  'Pband':      {'band':'P', 'fmin':1.e8,   'fmax':1.e9  },
                  'Lband':      {'band':'L', 'fmin':0.9e9,  'fmax':2.1e9 },
                  'Plo':        {'band':'P', 'fmin':2.4e8,  'fmax':3.4e8 },
                  'Phi':        {'band':'P', 'fmin':3.4e8,  'fmax':4.8e8 },
                  'Llo':        {'band':'L', 'fmin':1.0e9,  'fmax':1.4e9 },
                  'Lhi':        {'band':'L', 'fmin':1.4e9,  'fmax':2.0e9 },
                  'Slo':        {'band':'S', 'fmin':2.0e9,  'fmax':2.8e9 },
                  'Shi':        {'band':'S', 'fmin':2.8e9,  'fmax':4.0e9 }}
    flims = flims_dict[survey]
    return flims['band'], flims['fmin'],flims['fmax']

band,fmin,fmax = get_flims(survey)
print 'Generating tseries for', survey,'using', band,'band data from', fmin/1e9,'to', fmax/1e9, 'GHz\n'
savedir = '/data/jrv/burst_paper/rate/'
savefile = savedir + survey + '_tseries.npy'

filelist = load_filelist(band)
times = []
flux = []
flux_err = []
d = []
src_list = []
#filelist=['/data/jrv/15A-416/YZCMi/1/L/YZCMi_1L.tbavg.ms.dynspec'] # for testing
for f in filelist:
    src = get_srcname(f)
    dist = get_dist(src)
    params={'filename':f}
    ds = Dynspec(params)
    if ds.f[-1]<1.e9: # P band
        ds.spec['i'] = (ds.spec['xx']+ds.spec['yy'])/2
        #ds.mask_RFI(rmsfac=1.5) - conclusion: better w/o RFI masking
        #ds.spec['v'] = (ds.spec['xy']-ds.spec['yx'])/(2.j)
        del ds.spec['xx']
        del ds.spec['yy']
    else:
        ds.spec['i'] = (ds.spec['rr']+ds.spec['ll'])/2
        del ds.spec['rr']
        del ds.spec['ll']
    nt = int(integration_time/ds.dt())
    ds = ds.bin_dynspec(nt=nt,nf=1)
    tseries = ds.tseries(fmin,fmax)  # generates tseries weighted by 1/variance of im(vis) (downweights strong RFI to reduce scatter in tseries)
    t = tseries.time.mjds()
    s = tseries.spec['i']
    ind = find(~s.mask) # in case some times are masked, select only unmasked times
    ti = t[ind]
    si = s[ind].data
    s_1pc = si * dist**2
    rms = std(imag(s_1pc)) * ones(len(s_1pc))
    d += list(dist * ones(len(s_1pc)))
    src_list += [src] * len(s_1pc)
    times += list(t)
    flux += list(s_1pc)
    flux_err += list(rms)

times = array(times)
flux = array(flux)
flux_err = array(flux_err)
d = array(d)
save_dict  = {'times':times,'flux':flux,'flux_err':flux_err,'dist':d,'src':src_list,'t_int':integration_time}
save(savefile,save_dict)
print '\nSaved tseries for', survey, 'survey to', savefile

