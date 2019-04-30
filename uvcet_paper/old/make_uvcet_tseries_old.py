# -*- coding: utf-8 -*-
"""
display_dynspec.py: Script to load and plot dynamic spectra, using
the class dynspec.plot.Dynspec.
"""

from pylab import *
from numpy import *
#import matplotlib.pyplot as plt
import os
#from astropy.time import Time
#from copy import deepcopy
from dynspec.plot import *
import scipy.signal as signal
#import pickle

def mjds_to_hjds(mjds,sb):
    dt_min_dict = {'UVCet_2':    -5.753424060170966,
               'UVCet_3':    -0.48225985388936543,
               'UVCet_4':    1.2457881442377725,
               'UVCet_5':    6.288634966106827} # calculated for obs dates using http://www.physics.sfasu.edu/astro/javascript/hjd.html
    dt = dt_min_dict[sb]*60 # HJD-MJD (sec)
    return mjds + dt

try:
    overwrite
except:
    overwrite=False

close('all')

# font size for plots
rcParams.update({'font.size': 10})

#mydir = os.getcwd()
#tmp = mydir[1:].split('/')
#proj = tmp[2]
#src = tmp[3]
#obs = tmp[4]
#sb = proj + '_' + src + '_' + obs

mydir = '/data/jrv/15A-416/UVCet/'
src = 'UVCet'
obs_list = ['2','3','4','5']
bands = ['S','L']

filelist=[]
times = {'L':[],'S':[]}
flux = {'L':[],'S':[]}
i=0
for band in bands:
    leg = []
    i+=1
    subplot(2,1,i)
    for obs in obs_list:
        sb = src+'_'+obs
        filedir = mydir + obs + '/' + band + '/'
        savefile = filedir + band + '_tseries.npy'
        timefile = filedir + band + '_times.npy'
        if not os.path.exists(savefile) or not os.path.exists(timefile) or overwrite:
            dsfile = filedir + sb + band + '.tbavg.ms.dynspec'
            print 'Loading dynamic spectrum from', dsfile
            params={'filename':dsfile,'uniform':True}
            if band=='S' or band=='C':
                params['df']=2.e6
            ds_big = Dynspec(params)
            ds = ds_big.bin_dynspec(30,2) # rebin dynspec to 30 sec, 2 MHz resolution
            rms = std(ds.spec['ll'],0) # rms for each frequency
            tseries = average(ds.spec['rr'],1,1/rms**2) # compute variance-weighted RR tseries (so RFI doesn't dominate)
            t = ds.time.mjds()
            print 't0:',ds.t0()
            print 'Saving tseries to',savefile
            save(savefile,tseries)
            print 'Saving times to',timefile,'\n'
            save(timefile,t)
        else:
            t = load(timefile)
            tseries = load(savefile)
        t = mjds_to_hjds(t,sb)
        ind = find( tseries != 0 )
        times[band] += t[ind].tolist()
        flux[band] += tseries[ind].tolist()
        plot(t[ind]-t[0],tseries[ind])
        leg.append(obs+band)
    legend(leg)
    filename='tseries/'+band+'band_tseries.dat'
    savetxt(filename,column_stack((times[band],flux[band])))
savefig('tseries/tseries_band.pdf',bbox_inches='tight')

figure()
for band in bands:
    plot(times[band],flux[band])
legend(bands)
savefig('tseries/long.pdf',bbox_inches='tight')

# frequency range to search w/ LS periodogram
Pmin = 1. # hrs
Pmax = 25.
wmin = 2*pi/(Pmax * 3600)
dw = wmin/1000
wmax = 2*pi/(Pmin * 3600)
w = arange(wmin,wmax,dw) # range of omegas for which to compute periodogram
Prange = 2 * pi / w / 3600

figure()
# test Lomb-Scargle periodogram - generate sine wave w/ 5.5-hr period at times
#  specified by times list (just use one band), then run LS, see if we recover period
band = bands[0]
t = array(times[band])
P = 3600 * 5.5
f = 1.0 * (sin(2*pi/P * t) > 0.6)
pgram = signal.lombscargle(t, f, w)
subplot(311)
plot(t,f,'-.')
subplot(312)
plot(Prange,pgram,'-*')
subplot(313)
plot(Prange,pgram,'-*')
v=axis()
axis([5,6,v[2],v[3]])
savefig('tseries/test.pdf',bbox_inches='tight')

close('all')
figure()
# compute Lomb-Scargle periodogram for each band
i = 0
for band in bands:
    i += 1
    t = array(times[band])
    f = array(flux[band]) - mean(flux[band]) # traditional L-S periodogram assumes mean is zero
    pgram = signal.lombscargle(t,f,w)
    subplot(3,2,i)
    plot(t,f,'.')
    subplot(3,2,2+i)
    plot(Prange,pgram,'-*')
    subplot(3,2,4+i)
    plot(Prange,pgram,'-*')
    v=axis()
    axis([5,6,v[2],v[3]])
subplot(321)
title(bands[0])
subplot(322)
title(bands[1])
savefig('tseries/uvcet_ls.pdf',bbox_inches='tight')

