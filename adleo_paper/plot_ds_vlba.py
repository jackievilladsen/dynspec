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
epochlist = ['3','4','5']
epochlist = ['3']

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

# Cycle through epochs
for epoch in epochlist:

    close('all')
    figure(figsize=(6,6))

    # load multi-band dynspec
    savefile = savedir + src + '_' + epoch + '.dynspec.pickle'
    ds_dict = pickle.load( open( savefile, "rb" ) )
    ds = ds_dict['VLA']
    dsVLBA = ds_dict['VLBA']
    
    # Calculate VLBA time series
    dsVLBA = dsVLBA.bin_dynspec(nt=nt,nf=1)
    dsVLBA = dsVLBA.tseries(weight_mode='flat')
    t = dsVLBA.get_tlist() * dsVLBA.dt()/60.
    iflux = dsVLBA.spec['i']*1e3
    vflux = dsVLBA.spec['v']*1e3
    flux_err = std(imag(iflux))
    
    # bin dynspec to improve signal-to-noise ratio
    #ds_bin = ds.bin_dynspec(nt,nf)
    ds_bin = ds.bin_dynspec(nt,nf,mask_partial=0.75) # will flag if 3/4 of contributing pixels flagged --> 50% sensitivity
    # nt=15, nf=8 comes from high-resolution plots of bright burst
    # maybe use lower res for all bands then do second plot with high res of just bright burst
    
    #smax = min(ds_bin.rms_spec('v'))*rmsfac
    smin = -smax
    pp = {'pol':'v','smin':smin,'smax':smax,'func':real}

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
