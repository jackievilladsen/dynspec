'''
plot_adleo3_Lband_tseries.py - Plot tseries of coherent emission in L band (plot Stokes V to avoid bg contamination?)
'''

from dynspec.plot import Dynspec
from pylab import *
import pickle
#import os
'''
n_sec_VLA = 1
n_MHz = 16
mask_partial=0.75
weight_mode='rms'

src = 'ADLeo'
epoch='3'

params = {'legend.fontsize': 'small',
          'axes.titlesize': 'medium',
          'axes.labelsize': 'small',
          'xtick.labelsize': 'x-small',
          'ytick.labelsize': 'x-small',
          'image.interpolation': 'hanning'}
rcParams.update(params)

savedir = '/data/jrv/burst_paper/'
dsfile = '/data/jrv/15A-416/ADLeo/3/L/test_clean/ds_ap0_big_RR_n2_ms/tbavg.ms.dynspec/'

nt = n_sec_VLA  # number of integrations to bin together (current resolution is 1 sec)
nf = n_MHz  # number of channels to bin together (current resolution is 1 MHz)

close('all')
figure(figsize=(6,6))

# load multi-band dynspec
params={'filename':dsfile,'uniform':True}
ds = Dynspec(params)
ds.spec['i'] = (ds.spec['rr']+ds.spec['ll'])/2
ds.spec['v'] = (ds.spec['rr']-ds.spec['ll'])/2

# Bin dynspec to desired time resolution
#ds_bin = ds.bin_dynspec(nt=nt,nf=1,mask_partial=mask_partial)
ds_bin = ds

# make time series
ds12 = ds_bin.tseries(fmin=1.1e9,fmax=1.5e9,weight_mode=weight_mode)
t12 = ds12.get_tlist() * ds12.dt()/60.
i12 = ds12.spec['i']*1e3
v12 = ds12.spec['v']*1e3
i_err12 = std(imag(i12))
v_err12 = std(imag(v12))
'''

rr12 = ds12.spec['rr']*1e3
ll12 = ds12.spec['ll']*1e3
rr_err12 = std(imag(rr12))
ll_err12 = std(imag(ll12))

clf()
plot(t12,real(rr12),'b')
plot(t12,real(ll12),'g')
plot(t12,imag(rr12),'b--')
plot(t12,imag(ll12),'g--')
'''
plot(t12,real(i12),'b')
plot(t12,real(-v12),'g')
plot(t12,imag(i12),'b--')
plot(t12,imag(v12),'g--')
'''
xlabel('Time in minutes')
ylabel('Flux Density (mJy)')
legend(('I','V','Imag(I)','Imag(V)'))
gca().axhline(0,color='k')
savefig(savedir+'ADLeo3_tseries.pdf',bbox_inches='tight')
close('all')

# the burst does not appear to be composed of individual spikes, even at 1 second level
