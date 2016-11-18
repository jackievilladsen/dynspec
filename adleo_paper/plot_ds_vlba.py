'''
plot_ds_vlba.py - Load ADLeo VLA multi-band dynamic spectrum for a given epoch, bin to specified resolution, and plot with VLBA time series
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
epochlist = ['3','4','5']

params = {'legend.fontsize': 'small',
          'axes.titlesize': 'medium',
          'axes.labelsize': 'small',
          'xtick.labelsize': 'x-small',
          'ytick.labelsize': 'x-small',
          'image.interpolation': 'hanning'}
rcParams.update(params)

savedir = '/data/jrv/adleo_paper/'

# Cycle through epochs
for epoch in epochlist:

    # load multi-band dynspec
    savefile = savedir + src + '_' + epoch + '.dynspec.pickle'
    ds = pickle.load( open( savefile, "rb" ) )
    
    ### Plot dynspec ###

    # bin dynspec to improve signal-to-noise ratio
    nt = n_sec/6  # number of integrations to bin together (current resolution is 6 sec)
    nf = n_MHz/2  # number of channels to bin together (current resolution is 2 MHz)
    ds_bin = ds.bin_dynspec(nt,nf)
    # nt=15, nf=8 comes from high-resolution plots of bright burst
    # maybe use lower res for all bands then do second plot with high res of just bright burst

    #smax = min(ds_bin.rms_spec('v'))*rmsfac
    smin = -smax
    pp = {'pol':'v','smin':smin,'smax':smax}
    
    # create figure and save
    close('all')
    figure(figsize=(6,6))
    ds_bin.plot_dynspec(plot_params=pp)
    title(src + ' ' + epoch)
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
