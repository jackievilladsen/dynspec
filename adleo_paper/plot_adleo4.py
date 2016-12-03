# -*- coding: utf-8 -*-
"""
display_dynspec.py: Script to load and plot dynamic spectra, using
the class dynspec.plot.Dynspec.
"""

from pylab import *
from numpy import *
import matplotlib.pyplot as plt
import os
from astropy.time import Time
from copy import deepcopy
from dynspec.plot import *

close('all')

# mpl.rcParams['image.interpolation'] = 'bilinear' # default setting
# other options: 'nearest', 'bilinear', 'bicubic', 'spline16', 'spline36', 'hanning', 'hamming', 'hermite', 'kaiser',
#        'quadric', 'catrom', 'gaussian', 'bessel', 'mitchell', 'sinc', 'lanczos'
mpl.rcParams['image.interpolation'] = 'hanning'

# add tailored parameters (such as time range) here
plot_params_dict = {
#'15A-416_ADLeo_4': {'flims':array([1.e9,1.9e9]),'dx':30,'smin':0.02,'smax':0.2,'nt':8}, # L band
'15A-416_ADLeo_4': {'dx':1.0,'dy':0.1,'smin':-0.03,'smax':0.03,'nt':4,'nf':32,'scale':'linear'}, # short burst
}
flims = array([1.6e9,4.0e9])
tlims = array([133.,145.])

##### USER-DEFINED PARAMETERS #####

# scaling factors for color plots
plot_scale_factor = 5

# font size for plots
rcParams.update({'font.size': 10})

# properties for binned dynspec plots
ar = 1.0 # aspect ratio

##### END USER-DEFINED PARAMETERS #####

mydir = os.getcwd()
tmp = mydir[1:].split('/')
proj = tmp[2]
src = tmp[3]
obs = tmp[4]
sb = proj + '_' + src + '_' + obs

pp = plot_params_dict.get(sb,{})
plot_params = pp

vpp = {'scale':'linear'}
ipp = {}
try:
    smax2 = pp['smax']/2
    vpp['smax'] = smax2
    vpp['smin'] = -smax2
    ipp['smax'] = smax2
except:
    pass
try:
    smin2 = pp['smin']/2
    ipp['smin'] = smin2
except:
    pass
V_plot_params = pp.copy()
V_plot_params.update(vpp)

I_plot_params = pp.copy()
I_plot_params.update(ipp)

rcpp = {'scale':'linear','smin':-1,'smax':1}
rc_plot_params = pp.copy()
rc_plot_params.update(rcpp)

bands = plot_params.get('bands',['S','L'])

# number of time and frequency channels to bin together
nt = plot_params.get('nt',4)
nf = plot_params.get('nf',8)

# parameter filelist is not currently used - will eventually be used to cycle through files to load
filelist={}
for band in bands:
    filename = band + '/' + src + '_' + obs + band + '.tbavg.ms.dynspec'
    filelist[band] = filename
print 'file list:', filelist

try:
    ds
except:
    ds = None
    for band in bands:
        params={'filename':filelist[band],'uniform':True}
        if band=='S' or band=='C':
            params['df']=2.e6
        ds_band = Dynspec(params)
        if band=='L':
            ds_band = ds_band.bin_dynspec(nt=1,nf=2) # bin L-band dynamic spectrum to 2-MHz channels
            ds_band.f -= 0.5e6   # this is a quick fix; figure out how to address this better later
        if ds is None:
            ds = deepcopy(ds_band)
        else:
            ds.add_dynspec(ds_band)
        ds_band = None

ds_bin = ds.bin_dynspec(nt,nf)
ds_bin = ds_bin.clip(tmin=tlims[0],tmax=tlims[1],fmin=flims[0],fmax=flims[1])
ds_bin.spec['v'] = (ds_bin.spec['rr']-ds_bin.spec['ll'])/2
ds_bin.spec['i'] = (ds_bin.spec['rr']+ds_bin.spec['ll'])/2
ds_bin.spec['rc'] = real(ds_bin.spec['v'])/real(ds_bin.spec['i'])
#ds_bin.mask_RFI(5.)
'''
pol_dict = {'ll':'LCP','rr':'RCP','v':'StokesV','i':'StokesI','rc':'rc'}
for pol in ['ll','rr','i','v','rc']:
    close('all')
    if pol is 'v':
        pp = V_plot_params
    elif pol is 'i':
        pp = I_plot_params
    elif pol is 'rc':
        pp = rc_plot_params
    else:
        pp = plot_params
    pp['pol'] = pol
    ds_bin.plot_dynspec(plot_params=pp)
    p = pol_dict[pol]
    title(p+' Dynamic Spectrum')
    plotfile = 'plots/'+p+'dynspec.png'
    if os.path.exists(plotfile):
        os.system('rm -f ' + plotfile)
    savefig(plotfile,bbox_inches='tight')
    
    clf()
    pp['func']=imag
    ds_bin.plot_dynspec(plot_params=pp)
    p = pol_dict[pol]
    title('Imag('+p+') Dynamic Spectrum')
    plotfile = 'plots/'+p+'dynspec_imag.png'
    if os.path.exists(plotfile):
        os.system('rm -f ' + plotfile)
    savefig(plotfile,bbox_inches='tight')
    pp['func']=real
'''

# clip out RCP burst dynspec
dsR = ds_bin.clip(tmin=2.5,tmax=4.5,fmin=1.6e9,fmax=2.3e9)
clf()
dsR.plot_dynspec(plot_params={'pol':'rr','smin':0.01,'df':0.1,'dx':1.0})
savefig('plots/dsR.pdf',bbox_inches='tight')

# calculate flux-weighted mean time
t = dsR.get_tlist() * dsR.dt() # time in seconds
f = dsR.f/(1.e6)  # frequency in MHz

def t_avg(spec,t):
   return sum(t[:,newaxis] * spec,0) / sum(spec,0)

spec_real = abs(real(dsR.spec['rr']))
spec_imag = abs(imag(dsR.spec['rr']))
treal = t_avg(spec_real,t)
timag = t_avg(spec_imag,t)

from scipy.stats import linregress
slope, intercept, r_value, p_value, std_err = linregress(f,treal)
slope_imag, intercept_imag, r_imag, p_imag, err_imag = linregress(f,timag)

yfit = intercept + slope * f
yfit_imag = intercept_imag + slope_imag * f

clf()
#dsR.plot_dynspec(plot_params={'pol':'rr','smin':0.01})
#plot(f/1.e3,treal,'.b')
#plot(f/1.e3,yfit,'-b')

drift_rate = 1/slope # MHz/sec
print 'Drift rate:', drift_rate, 'MHz/sec'

#show()


tseries1 = dsR.tseries(fmin=1.6e9,fmax=1.8e9)
tseries2 = dsR.tseries(fmin=1.8e9,fmax=2.0e9)
t = tseries1.get_tlist() * tseries1.dt()
plot(t,tseries1.spec['rr'])
plot(t,tseries2.spec['rr'])
xlabel('Time in seconds')
ylabel('Flux (mJy) in RCP')
legend(('1.6-1.8 GHz','1.8-2 GHz'))
savefig('plots/drift.pdf',bbox_inches='tight')

### LCP portion of burst ###

dsL = ds_bin.clip(tmin=2.5,tmax=4.5,fmin=2.5e9,fmax=4.0e9)
