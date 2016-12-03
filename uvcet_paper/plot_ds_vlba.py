'''
plot_ds_vlba.py - Load ADLeo or UVCet multi-band dynamic spectrum for a given epoch, bin to specified resolution, and plot to file
'''

#from dynspec import load_dict
#from dynspec.plot import Dynspec
from pylab import *
import pickle
#import os

n_sec_P = 120 # must be multiple of 6
n_sec_VLA = 30
n_sec_VLBA = 60
n_MHz = 16  # must be multiple of 2
rmsfac = 10
smax = 0.08
smin = 0.002
scale = 'log'

src = 'UVCet'
epochlist = ['3','4','5']

params = {'legend.fontsize': 'small',
          'axes.titlesize': 'medium',
          'axes.labelsize': 'small',
          'xtick.labelsize': 'x-small',
          'ytick.labelsize': 'x-small',
          'image.interpolation': 'hanning'}
rcParams.update(params)

savedir = '/data/jrv/'+src.lower()+'_paper/'

# Cycle through epochs
for epoch in epochlist:

    close('all')
    figure(figsize=(6,6))
    
    # load multi-band dynspec
    savefile = savedir + src + '_' + epoch + '.dynspec.pickle'
    ds_dict = pickle.load( open( savefile, "rb" ) )
    ds = ds_dict['VLA']
    dsP = ds_dict['P']
    dsVLBA = ds_dict['VLBA']
    
    # Calculate number of dynspec pixels to bin together
    nt_P = int(round(n_sec_P/dsP.dt()))  # number of integrations to bin together (current resolution is 6 sec)
    nt_VLA = int(round(n_sec_VLA/ds.dt()))
    nt_VLBA = int(round(n_sec_VLBA/dsVLBA.dt()))
    nf = n_MHz/int(round(ds.df()/1e6))  # number of channels to bin together (current resolution is 2 MHz)
    nf_VLBA = n_MHz/int(round(dsVLBA.df()/1e6))
    
    # Calculate VLBA time series
    dsVLBA_bin = dsVLBA.bin_dynspec(nt=nt_VLBA,nf=nf_VLBA)
    tseriesVLBA = dsVLBA_bin.tseries(weight_mode='flat')
    t = tseriesVLBA.get_tlist() * tseriesVLBA.dt()/60.
    iflux = tseriesVLBA.spec['i']*1e3
    vflux = tseriesVLBA.spec['v']*1e3
    flux_err = std(imag(iflux))

    # bin dynspec to improve signal-to-noise ratio
    #ds_bin = ds.bin_dynspec(nt,nf)
    ds_bin = ds.bin_dynspec(nt_VLA,nf,mask_partial=0.75) # will flag if 3/4 of contributing pixels flagged --> 50% sensitivity
    dsP_bin = dsP.bin_dynspec(nt_P,nf,mask_partial=0.75)
    # nt=15, nf=8 comes from high-resolution plots of bright burst
    # maybe use lower res for all bands then do second plot with high res of just bright burst
    
    #smax = min(ds_bin.rms_spec('v'))*rmsfac
    pp = {'pol':'v','smin':smin,'smax':smax,'trim_mask':False,'ar0':0.8,'scale':scale,'axis_labels':['ylabel','cbar','cbar_label'],'dy':0.5}
    ppP = {'pol':'v','smin':smin,'smax':smax,'ar0':0.2,'trim_mask':False,'scale':scale,'axis_labels':['xlabel','cbar'],'dy':0.1}
    ppX = {'pol':'v','smin':smin,'smax':smax,'ar0':0.2,'trim_mask':False,'scale':scale,'axis_labels':['cbar'],'dy':0.1}
    # plot tseries above dynspec (time axis aligned) and save figure
    # to do: set x-axis limits for time series
    clf()
    ax=axes([0,0.18,1,0.81])
    ds_bin.plot_dynspec(plot_params=pp)
    ax.xaxis.set_visible(False)
    ax = axes([0,0.97,1,0.2])
    dsVLBA_bin.plot_dynspec(plot_params=ppX)
    cb = gca().images[-1].colorbar
    cb.remove()
    ax.xaxis.set_visible(False)
    ax=axes([0,1.17,0.914,0.2]) #2: 0.96
    axhline(0,color='k')
    errorbar(t,iflux,flux_err)
    plot(t,vflux)
    xlim([min(t),max(t)])
    ylabel('VLBA Flux (mJy)')
    title(src + ' ' + epoch)
    ax.xaxis.set_visible(False)
    ax = axes([0,0,1,0.2])
    dsP_bin.plot_dynspec(plot_params=ppP)
    cb = gca().images[-1].colorbar
    cb.remove()
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
