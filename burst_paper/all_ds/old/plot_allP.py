'''
plot_all_ds.py - Load all Stokes I,V dynamic spectra from file produced by compile_ds.py, then produce
                 a series of figures showing all of the dynamic spectra - Re(vis) and Im(vis).
'''

import dynspec.plot
reload(dynspec.plot)

from dynspec import load_dict
from dynspec.plot import *
from pylab import *
import os, subprocess
import matplotlib.gridspec as gridspec

rmsfac=3

def get_obsname(obsfile):
    # take a file directory such as '/data/jrv/15A-416/YZCMi/1' and
    # convert to obs name such as '15A-416_YZCMi_1' and srcname 'YZCMi'
    names = obsfile.split('/')
    srcname = names[4]
    obsname = names[3]+'_'+names[4]+'_'+names[5]
    return obsname,srcname

params = {'legend.fontsize': 'small',
          'axes.titlesize': 'small',
          'axes.labelsize': 'x-small',
          'xtick.labelsize': 'xx-small',
          'ytick.labelsize': 'xx-small',
          'image.interpolation': 'hanning'}
rcParams.update(params)

savefile = '/data/jrv/burst_paper/ds/all_burst_Pdynspec.npy'
ds_list = load_dict(savefile)

ds_dir = '/data/jrv/burst_paper/ds/' # where to save ds plots
if not os.path.exists(ds_dir):
    os.system('mkdir '+ds_dir)

obs = ds_list.keys()[0]
ds = ds_list[obs]

close('all')

### PLOT INDIVIDUAL OBSERVATIONS ###

srclist = ['ADLeo','UVCet','EQPeg']
#srclist = ['UVCet']

figure(figsize=(7,10))

n_rows = 5
n_cols = 3
gs = gridspec.GridSpec(n_rows, n_cols)
#ar0 = 0.55
ar0 = 0.9

proj = '15A-416'

for src in srclist:

    figname = ds_dir + src + '_' + proj + '_P.pdf'
    clf()
    sub_row = 1
    subplots_adjust(hspace=0.5,wspace=0)
    
    obslist = [s for s in ds_list.keys() if proj in s and src in s] # choose the obs names that have proj and src in them (eg all ADLeo 15A-416)
    obslist = sort(obslist)
    
    for obs in obslist:
        
        offset = (sub_row - 1) * n_cols
        
        print '\n', obs
        obsname,srcname = get_obsname(obs)
        ds = ds_list[obs]
        
        # flux limits for Stokes I
        smax = ds.get_rms('i')*6
        smin = -smax  # make colorbar symmetric about zero to be consistent with Stokes V
        
        # plot Stokes I real
        i = offset + 0
        subplot(gs[i])
        pp = {'pol':'i','smin':smin,'smax':smax,'trim_mask':False,'axis_labels':['cbar','ylabel'],'dy':0.05,'ar0':ar0}
        ds.plot_dynspec(plot_params=pp)
        if sub_row==1:
            title('Stokes I')
        
        # flux limits for Stokes U and V
        smax = ds.get_rms('v')*6
        smin = -smax  # make colorbar symmetric about zero to be consistent with Stokes V

        # plot Stokes U real
        i = offset + 1
        subplot(gs[i])
        pp = {'pol':'u','smin':smin,'smax':smax,'trim_mask':False,'axis_labels':['cbar','xlabel'],'ar0':ar0}
        ds.plot_dynspec(plot_params=pp)

        gca().yaxis.set_visible(False)
        gca().xaxis.set_label_coords(-0.15,-0.2)
        if sub_row==1:
            title('Stokes U')

        # plot Stokes V real
        i = offset + 2
        subplot(gs[i])
        pp = {'pol':'v','smin':smin,'smax':smax,'trim_mask':False,'axis_labels':['cbar','cbar_label'],'ar0':ar0}
        ds.plot_dynspec(plot_params=pp)
        gca().yaxis.set_visible(False)
        if sub_row==1:
            title('Stokes V')        
        
        sub_row += 1

    suptitle(src[0:2]+' '+src[2:5]+' 20'+proj[0:2],y=0.95)
    savefig(figname,bbox_inches='tight')
    clf()
    subrow=1

