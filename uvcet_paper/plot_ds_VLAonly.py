'''
plot_ds_VLAonly.py - Load ADLeo or UVCet multi-band dynamic spectrum for a given epoch, bin to specified resolution, and plot to file
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

obs_dict = {'2013':['1'],'2015':['2']}
obs_dict = {'2015':['2']}
bandlist = ['L','S','C','P']

params = {'legend.fontsize': 'small',
          'axes.titlesize': 'medium',
          'axes.labelsize': 'small',
          'xtick.labelsize': 'x-small',
          'ytick.labelsize': 'x-small',
          'image.interpolation': 'hanning'}
rcParams.update(params)

savedir = '/data/jrv/'+src.lower()+'_paper/'

Pband = False
Cband = False
for year in obs_dict:
    epochlist = obs_dict[year]
    # Cycle through epochs
    for epoch in epochlist:

        close('all')
        figure(figsize=(6,6))

        # load multi-band dynspec
        savefile = savedir + year + '_' + src + '_' + epoch + '.dynspec.pickle'
        ds_full = pickle.load( open( savefile, "rb" ) )
        if min(ds_full.f) < 0.8e9: # if there is P band data
            dsP = ds_full.clip(fmax=5.5e8,trim_mask=False)
            ds = ds_full.clip(fmin=0.999e9,trim_mask=False)
            Pband = True
        else:
            ds = ds_full

        # bin dynspec to improve signal-to-noise ratio
        nt_VLA = int(round(n_sec_VLA/ds.dt()))
        nf = n_MHz/int(round(ds.df()/1e6))  # number of channels to bin together (current resolution is 2 MHz)
        ds_bin = ds.bin_dynspec(nt_VLA,nf,mask_partial=0.75) # will flag if 3/4 of contributing pixels flagged --> 50% sensitivity

        if Pband:
            nt_P = int(round(n_sec_P/dsP.dt()))  # number of integrations to bin together (current resolution is 6 sec)
            dsP_bin = dsP.bin_dynspec(nt_P,nf,mask_partial=0.75)    

        if max(ds.f)>5.e9:
            Cband = True
            ar_fac = 1.33
        else:
            ar_fac = 1
        pp = {'pol':'v','smin':smin,'smax':smax,'trim_mask':False,'ar0':0.8*ar_fac*2.25,'scale':scale,'axis_labels':['ylabel','cbar','cbar_label'],'dy':0.5}
        if not Pband:
            del pp['axis_labels']
        ppP = {'pol':'v','smin':smin,'smax':smax,'ar0':0.2*2,'trim_mask':False,'scale':scale,'axis_labels':['xlabel','cbar'],'dy':0.1}
        # plot tseries above dynspec (time axis aligned) and save figure
        # to do: set x-axis limits for time series
        clf()
        ax=axes([0,0.20,1,0.81*ar_fac])
        ds_bin.plot_dynspec(plot_params=pp)
        #ax = axes([0,0.97,1,0.2])
        title(year + ' ' + src + ' ' + epoch)
        if Pband:
            ax.xaxis.set_visible(False)
            ax = axes([0,0,1,0.18])
            dsP_bin.plot_dynspec(plot_params=ppP)
            cb = gca().images[-1].colorbar
            cb.remove()
        figfile = savedir + year + '_'+ src + '_' + epoch + '_dynspec.pdf'
        savefig(figfile,bbox_inches='tight')
