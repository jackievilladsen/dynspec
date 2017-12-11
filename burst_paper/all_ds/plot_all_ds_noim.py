'''

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
          'axes.labelsize': 'x-small',
          'xtick.labelsize': 'xx-small',
          'ytick.labelsize': 'xx-small',
          'image.interpolation': 'none'}
rcParams.update(params)

savefile = '/data/jrv/burst_paper/all_burst_dynspec.npy'
ds_list = load_dict(savefile)

ds_dir = '/data/jrv/burst_paper/ds/' # where to save ds plots
if not os.path.exists(ds_dir):
    os.system('mkdir '+ds_dir)

obs = ds_list.keys()[0]
ds = ds_list[obs]

close('all')

### PLOT INDIVIDUAL OBSERVATIONS ###

srclist = [('ADLeo','15A-416'),('UVCet','13A-423'),('UVCet','15A-416'),('EQPeg','15A-416'),('YZCMi','15A-416')]
#srclist = [('UVCet','15A-416')]

for src,proj in srclist:

    figure(figsize=(7,10))

    n_rows = 5
    n_cols = 2
    if proj == '13A-423':
        n_rows = 3
    gs = gridspec.GridSpec(n_rows, n_cols)

    figname = ds_dir + src + '_' + proj + '.pdf'
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

        ar0 = 0.45
        if proj == '13A-423':
            ar0 *= 5./3.
        if obs == '/data/jrv/15A-416/UVCet/1' or obs == '/data/jrv/15A-416/UVCet/2':
            ar0 *= 87./38.
        
        if obs == '/data/jrv/15A-416/UVCet/1':
            ds = ds.expand_flims(fmax_new = 4.e9)

        # flux limits
        smax = max(percentile(real(ds.spec['i']),99),median(real(ds.spec['i']))*2)
        smin = -smax  # make colorbar symmetric about zero to be consistent with Stokes V

        # plot Stokes I real
        i = offset + 0
        subplot(gs[i])
        pp = {'pol':'i','smin':smin,'smax':smax,'trim_mask':False,'axis_labels':['ylabel','cbar'],'ar0':ar0}
        ds.plot_dynspec(plot_params=pp)
        cb = gca().images[-1].colorbar
        cb.remove()
        if obs == '/data/jrv/15A-416/UVCet/1' or obs == '/data/jrv/15A-416/UVCet/2':
            pos1 = gca().get_position()
            gca().set_position([pos1.x0-0.09,pos1.y0,pos1.width,pos1.height])
        if sub_row==1:
            title('Stokes I')
        
        # plot Stokes V real
        i = offset + 1
        subplot(gs[i])
        pp = {'pol':'v','smin':smin,'smax':smax,'trim_mask':False,'axis_labels':['cbar','xlabel','cbar_label'],'ar0':ar0}
        ds.plot_dynspec(plot_params=pp)
        gca().yaxis.set_visible(False)
        gca().xaxis.set_label_coords(-0.15,-0.2)
        if obs == '/data/jrv/15A-416/UVCet/1' or obs == '/data/jrv/15A-416/UVCet/2':
            gca().xaxis.set_label_coords(-1.0,-0.2)
            pos1 = gca().get_position()
            gca().set_position([pos1.x0-0.09,pos1.y0,pos1.width,pos1.height])
            cb = gca().images[-1].colorbar
            x = cb.ax.get_position()
            cb.ax.set_position([x.x0-0.09,x.y0,x.width,x.height])
        if sub_row==1:
            title('Stokes V')
        
        sub_row += 1

    suptitle(src[0:2]+' '+src[2:5]+' 20'+proj[0:2],y=0.95)
    savefig(figname,bbox_inches='tight')
    clf()
    subrow=1


