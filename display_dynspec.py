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
#'13A-423_UVCet_1': {'bands':['S','L','C'],'tlims':array([75.,125.]),'dx':10,'smin':0.01,'smax':0.1}, # long
#'13A-423_UVCet_1': {'bands':['S','L','C'],'tlims':array([75.,125.]),'dx':10,'smin':0.0025,'smax':0.025,'nt':16,'nf':32}, # long bright
#'13A-423_UVCet_1': {'bands':['S','L','C'],'tlims':array([15.,45.]),'dx':5,'smin':0.0025,'smax':0.025,'nt':16,'nf':32}, # short
'13A-423_UVCet_1': {'bands':['S','L','C'],'smin':0.005,'smax':0.05,'nt':30,'nf':32}, # full
'13A-423_UVCet_2': {'bands':['S','L','C'],'smin':0.003,'smax':0.03,'nt':30,'nf':32},
'13A-423_UVCet_3': {'bands':['S','L','C'],'smin':0.001,'smax':0.01,'nt':150,'nf':64},
'13A-423_ADLeo_1': {'bands':['S','C'],'smin':0,'smax':0.005,'nt':150,'nf':64,'scale':'linear'},
'13A-423_ADLeo_2': {'bands':['S','L','C'],'smin':0,'smax':0.005,'nt':150,'nf':64,'scale':'linear'},
'13A-423_ADLeo_3': {'bands':['S','L','C'],'smin':0,'smax':0.005,'nt':150,'nf':64,'scale':'linear'},
'13A-423_ADLeo_4': {'bands':['S','L','C'],'smin':0,'smax':0.005,'nt':150,'nf':64,'scale':'linear'},
'15A-416_ADLeo_1': {'smin':0,'smax':0.005,'nt':150,'nf':64,'scale':'linear'},
'15A-416_ADLeo_2': {'smin':0,'smax':0.005,'nt':150,'nf':64,'scale':'linear'},
#'15A-416_ADLeo_3': {'flims':array([1.e9,1.9e9]),'dx':30,'smin':0.01,'smax':0.1,'nt':15,'nf':8},
#'15A-416_ADLeo_3': {'dx':30,'smin':0.01,'smax':0.1,'nt':15,'nf':8},
'15A-416_ADLeo_3': {'dx':30,'smin':0,'smax':0.005,'nt':300,'nf':64,'scale':'linear'},
#'15A-416_ADLeo_4': {'flims':array([1.e9,1.9e9]),'dx':30,'smin':0.02,'smax':0.2,'nt':8}, # L band
'15A-416_ADLeo_4': {'flims':array([1.6e9,4.0e9]),'tlims':array([133.,145.]),'dx':5,'smin':-0.03,'smax':0.03,'nt':5,'nf':32,'scale':'linear'}, # short burst
'15A-416_ADLeo_5': {'dx':30,'smin':0.005,'smax':0.05,'nt':60,'nf':32},
'15A-416_EQPeg_1': {'smin':0,'smax':0.005,'nt':150,'nf':64,'scale':'linear'},
'15A-416_EQPeg_2': {'smin':0.002,'smax':0.02,'nt':60,'nf':32},
#'15A-416_EVLac_1': {'tlims':array([75.,120.]),'smin':0,'smax':0.01,'nt':15,'nf':32,'scale':'linear'}, # inspect burst of RFI(?) seen in Stokes V
'15A-416_EVLac_1': {'smin':0,'smax':0.005,'nt':300,'nf':64,'scale':'linear'},
'15A-416_EVLac_2': {'smin':0,'smax':0.005,'nt':150,'nf':64,'scale':'linear'},
'15A-416_UVCet_1': {'tlims':array([40.,100.]),'dx':10,'smin':0.005,'smax':.5,'nt':8,'nf':8,'bands':['L']},
'15A-416_UVCet_2': {'tlims':array([60.,110.]),'dx':10,'smin':0.03,'smax':0.3,'nt':4,'nf':8},
'15A-416_UVCet_3': {'tlims':array([20.,120.]),'dx':30,'smin':0.02,'smax':0.2},
#'15A-416_UVCet_4': {'tlims':array([180.,240.]),'dx':10,'smin':0.02,'smax':0.2,'nt':16,'nf':16}, # main burst
#'15A-416_UVCet_4': {'tlims':array([0.,40.]),'dx':10,'smin':0.01,'smax':0.1,'nt':16,'nf':16}, # small
'15A-416_UVCet_4': {'smin':0.001,'smax':0.1,'nt':60,'nf':32}, # full
'15A-416_UVCet_5': {'tlims':array([60.,165.]),'dx':10,'smin':0.02,'smax':0.2},
'15A-416_YZCMi_1': {'tlims':array([61.,64.]),'flims':array([1.7e9,3.8e9]),'dx':1,'smin':.018,'smax':.18,'nt':1,'nf':16}, # leading edge
#'15A-416_YZCMi_1': {'tlims':array([58.,73.]),'dx':5,'smin':.030,'smax':.180,'nt':4,'nf':8}, # burst
#'15A-416_YZCMi_1': {'tlims':array([58.,73.]),'dx':5,'smin':.0030,'smax':.03,'nt':15,'nf':16}, # burst faint
#'15A-416_YZCMi_1': {'dx':5,'smin':.0030,'smax':.03,'nt':30,'nf':32}, # full obs
'15A-416_YZCMi_2': {'smin':0.003,'smax':0.03,'nt':30,'nf':32}
}

V_params_dict = {
#'13A-423_UVCet_1': {'tlims':array([75.,125.]),'dx':10,'smin':-0.05,'smax':0.05,'scale':'linear'}, # long
#'13A-423_UVCet_1': {'tlims':array([75.,125.]),'dx':10,'smin':-0.01,'smax':0.01,'scale':'linear'}, # long bright
#'13A-423_UVCet_1': {'tlims':array([15.,45.]),'dx':5,'smin':-0.01,'smax':0.01,'scale':'linear'}, # short
'13A-423_UVCet_1': {'smin':-0.025,'smax':0.025,'scale':'linear'}, # full
'13A-423_UVCet_2': {'smin':-0.015,'smax':0.015,'scale':'linear'},
'13A-423_UVCet_3': {'smin':-0.005,'smax':0.005,'scale':'linear'},
'13A-423_ADLeo_1': {'smin':-0.0025,'smax':0.0025,'scale':'linear'},
'13A-423_ADLeo_2': {'smin':-0.0025,'smax':0.0025,'scale':'linear'},
'13A-423_ADLeo_3': {'smin':-0.0025,'smax':0.0025,'scale':'linear'},
'13A-423_ADLeo_4': {'smin':-0.0025,'smax':0.0025,'scale':'linear'},
#'15A-416_ADLeo_3': {'flims':array([1.e9,1.9e9]),'dx':30,'smin':-0.025,'smax':0.025,'scale':'linear'},
#'15A-416_ADLeo_4': {'flims':array([1.e9,1.9e9]),'dx':30,'smin':-0.07,'smax':0.07,'scale':'linear'},
'15A-416_ADLeo_5': {'dx':30,'smin':-0.01,'smax':0.01,'scale':'linear'},
'15A-416_EQPeg_2': {'smin':-0.01,'smax':0.01,'scale':'linear'},
'15A-416_UVCet_1': {'pol':'v','tlims':array([40.,120.]),'dx':10,'smin':-.05,'smax':.05,'scale':'linear'},
'15A-416_UVCet_2': {'pol':'v','tlims':array([60.,110.]),'dx':10,'smin':-0.1,'smax':0.1,'scale':'linear'},
'15A-416_UVCet_3': {'tlims':array([20.,120.]),'dx':30,'smin':-0.1,'smax':0.1,'scale':'linear'},
#'15A-416_UVCet_4': {'tlims':array([180.,240.]),'dx':10,'smin':-0.1,'smax':0.1,'scale':'linear'}, # main burst
#'15A-416_UVCet_4': {'tlims':array([0.,40.]),'dx':10,'smin':-0.05,'smax':0.05,'scale':'linear'}, # small
'15A-416_UVCet_4': {'smin':-0.05,'smax':0.05,'scale':'linear'}, # full
'15A-416_UVCet_5': {'tlims':array([60.,165.]),'dx':10,'smin':-0.1,'smax':0.1,'scale':'linear'},
'15A-416_YZCMi_2': {'smin':-0.015,'smax':0.015,'scale':'linear'}
}

I_params_dict = {
#'13A-423_UVCet_1': {'tlims':array([75.,125.]),'dx':10,'smin':0.005,'smax':0.05}, # long
#'13A-423_UVCet_1': {'tlims':array([75.,125.]),'dx':10,'smin':0.0025,'smax':0.025}, # long bright
#'13A-423_UVCet_1': {'tlims':array([15.,45.]),'dx':5,'smin':0.001,'smax':0.01}, # short
'13A-423_UVCet_1': {'smin':0.0025,'smax':0.025}, # full
'13A-423_UVCet_2': {'smin':0.0015,'smax':0.015},
'13A-423_UVCet_3': {'smin':0.001,'smax':0.01},
'13A-423_ADLeo_1': {'smin':0,'smax':0.0025,'scale':'linear'},
'13A-423_ADLeo_2': {'smin':0,'smax':0.0025,'scale':'linear'},
'13A-423_ADLeo_3': {'smin':0,'smax':0.0025,'scale':'linear'},
'13A-423_ADLeo_4': {'smin':0,'smax':0.0025,'scale':'linear'},
#'15A-416_ADLeo_3': {'flims':array([1.e9,1.9e9]),'dx':30,'smin':0.0025,'smax':0.025},
#'15A-416_ADLeo_4': {'flims':array([1.e9,1.9e9]),'dx':30,'smin':0.01,'smax':0.1},
'15A-416_ADLeo_5': {'dx':30,'smin':0.0025,'smax':0.025},
'15A-416_EQPeg_2': {'smin':0.001,'smax':0.01},
'15A-416_UVCet_1': {'pol':'i','tlims':array([40.,120.]),'dx':10,'smin':.01,'smax':.05},
'15A-416_UVCet_2': {'pol':'i','tlims':array([60.,110.]),'dx':10,'smin':0.015,'smax':0.15},
'15A-416_UVCet_3': {'tlims':array([20.,120.]),'dx':30,'smin':0.01,'smax':0.1},
#'15A-416_UVCet_4': {'tlims':array([180.,240.]),'dx':10,'smin':0.01,'smax':0.1}, # main burst
#'15A-416_UVCet_4': {'tlims':array([0.,40.]),'dx':10,'smin':0.005,'smax':0.5}, # small
'15A-416_UVCet_4': {'smin':0.005,'smax':0.5}, # full
'15A-416_UVCet_5': {'tlims':array([60.,165.]),'dx':10,'smin':0.01,'smax':0.1},
'15A-416_YZCMi_2': {'smin':0.0015,'smax':0.015}
}

rc_params_dict = {
#'13A-423_UVCet_1': {'tlims':array([75.,125.]),'dx':10,'smin':-1,'smax':1,'scale':'linear'}, # long, long bright
#'13A-423_UVCet_1': {'tlims':array([15.,45.]),'dx':5,'smin':-1,'smax':1,'scale':'linear'}, # short
'13A-423_UVCet_1': {'smin':-1,'smax':1,'scale':'linear'}, # full
'13A-423_UVCet_2': {'smin':-1,'smax':1,'scale':'linear'},
'13A-423_UVCet_3': {'smin':-1,'smax':1,'scale':'linear'},
'13A-423_ADLeo_1': {'smin':-1,'smax':1,'scale':'linear'},
'13A-423_ADLeo_2': {'smin':-1,'smax':1,'scale':'linear'},
'13A-423_ADLeo_3': {'smin':-1,'smax':1,'scale':'linear'},
'13A-423_ADLeo_4': {'smin':-1,'smax':1,'scale':'linear'},
#'15A-416_ADLeo_3': {'flims':array([1.e9,1.9e9]),'dx':30,'smin':-1,'smax':1,'scale':'linear'},
#'15A-416_ADLeo_4': {'flims':array([1.e9,1.9e9]),'dx':30,'smin':-1,'smax':1,'scale':'linear'},
'15A-416_ADLeo_5': {'dx':30,'smin':-1,'smax':1,'scale':'linear'},
'15A-416_EQPeg_2': {'smin':-1,'smax':1,'scale':'linear'},
'15A-416_UVCet_1': {'tlims':array([40.,120.]),'dx':10,'smin':-1,'smax':1,'scale':'linear'},
'15A-416_UVCet_2': {'tlims':array([60.,110.]),'dx':10,'smin':-1,'smax':1,'scale':'linear'},
'15A-416_UVCet_3': {'tlims':array([20.,120.]),'dx':30,'smin':-1,'smax':1,'scale':'linear'},
#'15A-416_UVCet_4': {'tlims':array([180.,240.]),'dx':10,'smin':-1,'smax':1,'scale':'linear'}, # main burst
#'15A-416_UVCet_4': {'tlims':array([0.,40.]),'dx':10,'smin':-1,'smax':1,'scale':'linear'}, # small
'15A-416_UVCet_4': {'smin':-1,'smax':1,'scale':'linear'}, # full
'15A-416_UVCet_5': {'tlims':array([60.,165.]),'dx':10,'smin':-1,'smax':1,'scale':'linear'},
'15A-416_YZCMi_2': {'smin':-1,'smax':1,'scale':'linear'}
}

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
V_plot_params.update(V_params_dict.get(sb,{}))

I_plot_params = pp.copy()
I_plot_params.update(ipp)
I_plot_params.update(I_params_dict.get(sb,{}))

rcpp = {'scale':'linear','smin':-1,'smax':1}
rc_plot_params = pp.copy()
rc_plot_params.update(rcpp)
rc_plot_params.update(rc_params_dict.get(sb,{}))

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
ds_bin.spec['v'] = (ds_bin.spec['rr']-ds_bin.spec['ll'])/2
ds_bin.spec['i'] = (ds_bin.spec['rr']+ds_bin.spec['ll'])/2
ds_bin.spec['rc'] = real(ds_bin.spec['v'])/real(ds_bin.spec['i'])
#ds_bin.mask_RFI(5.)

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
# calculate time series for ADLeo 3, 3-4 GHz
rms = std(imag(ds_bin.spec['v']),0) # rms for each frequency
#tseries = average(ds_bin.spec['v'],1,1/rms**2) # compute variance-weighted V tseries (so RFI doesn't dominate)
#tseries = average(ds_bin.spec['v'],1)
ds2 = ds_bin.bin_dynspec(1,15)
tseriesV = ds2.spec['v'][:,2]
tseriesI = ds2.spec['i'][:,2]
tseriesRR = ds2.spec['rr'][:,2]
tseriesLL = ds2.spec['ll'][:,2]
t = ds2.time.mjds()                                                                                 
print 't0:',ds2.t0()
close('all')
subplot(211)
plot(t/3600.,tseriesV*1000,'.')
plot(t/3600.,tseriesI*1000,'.')
legend(('V','I'))
ylabel('Flux (mJy)')
subplot(212)
plot(t/3600.,tseriesRR*1000,'.')
plot(t/3600.,tseriesLL*1000,'.')
legend(('RR','LL'))
xlabel('Time in hours since '+ds2.t0())
ylabel('Flux (mJy)')
savefig('plots/tseries_3to4GHz.png')

'''


