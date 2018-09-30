'''
plot_allband_ds.py - Load P,L,S band dynamic spectrum for a given epoch, bin to specified resolution, and plot to file
'''

import dynspec.plot
reload(dynspec.plot)

from dynspec import load_dict
from dynspec.plot import *
from pylab import *
import os, subprocess
import matplotlib.gridspec as gridspec

def get_obsname(obsfile):
    # take a file directory such as '/data/jrv/15A-416/YZCMi/1' and
    # convert to obs name such as '15A-416_YZCMi_1' and srcname 'YZCMi'
    names = obsfile.split('/')
    srcname = names[4]
    obsname = names[3]+'_'+names[4]+'_'+names[5]
    return obsname,srcname

params = {'legend.fontsize': 'small',
          'axes.titlesize': 'small',
          'axes.labelsize': 'small',
          'xtick.labelsize': 'x-small',
          'ytick.labelsize': 'x-small',
          'image.interpolation': 'nearest'}
rcParams.update(params)

loadfile = '/data/jrv/burst_paper/all_burst_epoch_dynspec_LSband.npy'
ds_list = load_dict(loadfile)
loadfileP = '/data/jrv/burst_paper/all_burst_epoch_dynspec_Pband.npy'
dsP_list = load_dict(loadfileP)

ds_dir = '/data/jrv/burst_paper/ds/' # where to save ds plots
if not os.path.exists(ds_dir):
    os.system('mkdir '+ds_dir)

close('all')

# default binning settings
tint_P = 300
tint_LS = 60
df_MHz_P = 16
df_MHz_LS = 16

### PLOT INDIVIDUAL OBSERVATIONS ###

obs_list = dsP_list.keys()  # only using this code to plot epochs with P band data
#obs_list = ['/data/jrv/15A-416/EQPeg/2'] # so I can work on just this event

for obs in obs_list:
    
    # load dynamic spectra for this observation
    print '\n', obs
    obsname,srcname = get_obsname(obs)
    ds = ds_list[obs]
    dsP = dsP_list[obs]
    
    #print 'Duration:',ds.get_tlist()[-1]*ds.dt(),'sec'
    # 2015_UVCet_5: 13164 sec
    # 2015_UVCet_2: 5940 sec
    # difference: 7224 sec
    # Pad 2015_UVCet_2 with zeros so that it's on the same time scale as the others
    #  (couldn't figure out how to do it with the axis ratios of the plots and get
    #   consistency)
    if obs == '/data/jrv/15A-416/UVCet/2':
        ds = ds.expand_tlims(t_add_left=0.,t_add_right=7224.)
        dsP = dsP.expand_tlims(t_add_left=0.,t_add_right=7224.)
    
    # bin LS band dynamic spectrum to desired resolution
    nt = int(round(tint_LS/ds.dt()))    # number of integrations to bin together
    nf = int(round(df_MHz_LS/(ds.df()/1e6))) # number of channels to bin together
    ds = ds.bin_dynspec(nt=nt,nf=nf,mask_partial=0.5)
    ds.mask_RFI(rmsfac=5.)

    # bin P band dynamic spectrum to desired resolution
    nt = int(round(tint_P/dsP.dt()))    # number of integrations to bin together
    nf = int(round(df_MHz_P/(dsP.df()/1e6))) # number of channels to bin together
    dsP = dsP.bin_dynspec(nt=nt,nf=nf,mask_partial=0.5)
    
    # create figure
    figname = ds_dir+obsname+'.pdf'
    fig=figure(figsize=(9,6.5))
    n_rows = 3
    n_cols = 3
    heights = [2.1,0.7,0.7]
    widths = [1,1,0.15]
    gs = gridspec.GridSpec(n_rows,n_cols,height_ratios=heights,width_ratios=widths) #,wspace=0.1,hspace=0.1)
    clf()
    
    # set flux limits for LS band
    smax = max(percentile(real(ds.spec['i']),99)*1.1,median(real(ds.spec['i']))*2)
    smin = -smax  # make colorbar symmetric about zero
    
    # set axis ratio for LS band plots
    ar0 = 0.9
    
    # plot Stokes I real, LS band
    fig.add_subplot(gs[0,0])
    pp = {'pol':'i','smin':smin,'smax':smax,'trim_mask':False,'axis_labels':['cbar','cbar_label','ylabel'],'ar0':ar0,'dy':0.5}
    plt,cbar_ticks,cbar_ticklbls = ds.plot_dynspec(plot_params=pp)
    cbLS = gca().images[-1].colorbar
    cbLS.remove()
    #gca().xaxis.set_visible(False)
    gca().yaxis.set_label_coords(-0.12,0)
    title('Stokes I, 1-4 GHz')
    
    # plot Stokes V real, LS band
    fig.add_subplot(gs[0,1])
    pp = {'pol':'v','smin':smin,'smax':smax,'trim_mask':False,'axis_labels':['cbar','cbar_label'],'ar0':ar0,'dy':0.5}
    ds.plot_dynspec(plot_params=pp)
    gca().yaxis.tick_right()
    cbLS = gca().images[-1].colorbar
    cbLS.remove()
    #gca().xaxis.set_visible(False)
    title('Stokes V, 1-4 GHz')

    # plot LS band colorbar
    gs_LSband_cbar = gridspec.GridSpecFromSubplotSpec(1,2, subplot_spec = gs[0,2], width_ratios=[0.8,0.3])
    ax = subplot(gs_LSband_cbar[1])
    cbar=colorbar(plt,cax=ax)
    cbar.set_ticks(cbar_ticks)
    cbar.set_ticklabels(cbar_ticklbls)
    #cbar.set_label('1-4 GHz Flux Density (mJy)')
    ax = cbar.ax
    ax.text(3.7,0.67,'1-4 GHz Flux Density (mJy)',rotation=90,fontsize='small')

    # set flux limits for P band
    smaxP = dsP.get_rms('i')*3.5
    #smaxP = 0.09
    sminP = -smaxP

    # set axis ratio for P band plots
    ar0 = 0.3
    
    # plot Stokes I real, P band
    fig.add_subplot(gs[1,0])
    pp = {'pol':'i','smin':sminP,'smax':smaxP,'trim_mask':False,'axis_labels':['cbar','cbar_label','xlabel'],'dy':0.05,'ar0':ar0}
    dsP.plot_dynspec(plot_params=pp)
    cbP = gca().images[-1].colorbar
    cbP.remove()
    title('Stokes I, 0.2-0.5 GHz')
    
    # plot Stokes V real, P band
    fig.add_subplot(gs[1,1])
    pp = {'pol':'v','smin':sminP,'smax':smaxP,'trim_mask':False,'axis_labels':['cbar','cbar_label'],'dy':0.05,'ar0':ar0}
    plt,cbar_ticks,cbar_ticklbls=dsP.plot_dynspec(plot_params=pp)
    gca().yaxis.tick_right()
    #gca().xaxis.set_visible(False)
    cbP = gca().images[-1].colorbar
    #gca().set_axis_bgcolor('#00ff00') #'0.5')
    cbP.remove()
    title('Apparent Stokes V, 0.2-0.5 GHz')
    
    # plot Stokes U real, P band
    fig.add_subplot(gs[2,1])
    pp = {'pol':'u','smin':sminP,'smax':smaxP,'trim_mask':False,'axis_labels':['cbar','cbar_label'],'dy':0.05,'ar0':ar0}
    gca().yaxis.tick_right()
    dsP.plot_dynspec(plot_params=pp)
    cbP = gca().images[-1].colorbar
    cbP.remove()
    title('Apparent Stokes U, 0.2-0.5 GHz')
    
    # plot P band colorbar
    gs_Pband_cbar = gridspec.GridSpecFromSubplotSpec(1,2, subplot_spec = gs[1:3,2], width_ratios=[0.8,0.3])
    ax = subplot(gs_Pband_cbar[1])
    cbar=colorbar(plt,cax=ax)
    cbar.set_ticks(cbar_ticks)
    cbar.set_ticklabels(cbar_ticklbls)
    ax = cbar.ax
    ax.text(3.7,0.8,'0.2-0.5 GHz Flux Density (mJy)',rotation=90,fontsize='small')

    date = dsP.t0().split()[0]
    suptitle(srcname[0:2]+' '+srcname[2:5]+' - '+date,y=1.05,fontsize='medium')
    tight_layout(pad=0.05)
    savefig(figname,bbox_inches='tight')
    clf()


'''
    #smax = min(ds_bin.rms_spec('v'))*rmsfac
#    pp = {'pol':'v','smin':smin,'smax':smax,'trim_mask':False,'ar0':0.8,'scale':scale,'axis_labels':['ylabel','cbar','cbar_label'],'dy':0.5}
#    ppP = {'pol':'v','smin':smin,'smax':smax,'ar0':0.2,'trim_mask':False,'scale':scale,'axis_labels':['xlabel','cbar'],'dy':0.1}
#    ppX = {'pol':'v','smin':smin,'smax':smax,'ar0':0.2,'trim_mask':False,'scale':scale,'axis_labels':['cbar'],'dy':0.1}
    
    # plot tseries above dynspec (time axis aligned) and save figure
    # to do: set x-axis limits for time series
    clf()
    ax=axes([0,0.18,1,0.81])                  # VLA LS-band dynspec
    ds_bin.plot_dynspec(plot_params=pp)
    ax.xaxis.set_visible(False)
    
    ax = axes([0,0.97,1,0.2])                 # VLBA dynspec
    dsVLBA_bin.plot_dynspec(plot_params=ppX)
    cb = gca().images[-1].colorbar
    cb.remove()
    ax.xaxis.set_visible(False)
    
    ax=axes([0,1.17,0.914,0.2]) #2: 0.96      # VLBA time series
    axhline(0,color='k')
    errorbar(t,iflux,flux_err)
    plot(t,vflux)
    xlim([min(t),max(t)])
    ylabel('VLBA Flux Density (mJy)')
    title(year + ' ' + src + ' ' + epoch)
    ax.xaxis.set_visible(False)
    
    ax = axes([0,0,1,0.2])                    # VLA P-band dynspec
    dsP_bin.plot_dynspec(plot_params=ppP)
    cb = gca().images[-1].colorbar
    cb.remove()

    figfile = savedir + src + '_' + epoch + '_dynspec.pdf'
    savefig(figfile,bbox_inches='tight')



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
