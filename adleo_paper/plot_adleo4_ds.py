'''
plot_adleo_ds.py - Load ADLeo multi-band dynamic spectrum for a given epoch, bin to specified resolution, and plot to file
'''

#from dynspec import load_dict
#from dynspec.plot import Dynspec
from pylab import *
import pickle
#import os

n_sec = 120 # must be multiple of 6
n_MHz = 16  # must be multiple of 2
rmsfac = 10
smax = 0.02

src = 'ADLeo'
epoch='4'

params = {'legend.fontsize': 'small',
          'axes.titlesize': 'medium',
          'axes.labelsize': 'small',
          'xtick.labelsize': 'x-small',
          'ytick.labelsize': 'x-small',
          'image.interpolation': 'hanning'}
rcParams.update(params)

savedir = '/data/jrv/adleo_paper/'

nt = n_sec/6  # number of integrations to bin together (current resolution is 6 sec)
nf = n_MHz/2  # number of channels to bin together (current resolution is 2 MHz)


close('all')
figure(figsize=(6,6))

# load multi-band dynspec
savefile = savedir + src + '_' + epoch + '.dynspec.pickle'
ds_dict = pickle.load( open( savefile, "rb" ) )
ds = ds_dict['VLA']
#dsVLBA = ds_dict['VLBA']

# bin dynspec to improve signal-to-noise ratio
ds_bin = ds.bin_dynspec(nt,nf)
# nt=15, nf=8 comes from high-resolution plots of bright burst
# maybe use lower res for all bands then do second plot with high res of just bright burst

smax = 0.01
smin = -smax
pp = {'pol':'v','smin':smin,'smax':smax}

# plot tseries above dynspec (time axis aligned) and save figure
# to do: set x-axis limits for time series
clf()
ax=axes([0,0,1,1])
ds_bin.plot_dynspec(plot_params=pp)
ax=axes([0,0.96,0.915,0.2]) #2: 0.96
axhline(0,color='k')
errorbar(t,iflux,flux_err)
plot(t,vflux)
xlim([min(t),max(t)])
ylabel('VLBA Flux (mJy)')
title(src + ' ' + epoch)
gca().xaxis.set_visible(False)
figfile = savedir + src + '_' + epoch + '_dynspec.pdf'
savefig(figfile,bbox_inches='tight')


'''
# subplots to show im(vis) also

subplot(121)
ds_bin.plot_dynspec(plot_params=pp)
cb = gca().images[-1].colorbar
cb.remove()
gca().xaxis.set_label_coords(1.2,-0.15)
title('Stokes V')
subplot(122)
pp['func']=imag
ds_bin.plot_dynspec(plot_params=pp)
title('Imag(Stokes V)')
suptitle(src + ' ' + epoch)
'''
