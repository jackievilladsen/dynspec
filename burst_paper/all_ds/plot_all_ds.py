'''
plot_all_ds.py - Load all Stokes I,V dynamic spectra from file produced by compile_ds.py, then produce
                 a series of figures showing all of the dynamic spectra - Re(vis) and Im(vis).
'''

from dynspec import load_dict
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

savefile = '/data/jrv/burst_paper/all_dynspec.npy'
ds_list = load_dict(savefile)

ds_dir = '/data/jrv/burst_paper/ds/' # where to save ds plots
if not os.path.exists(ds_dir):
    os.system('mkdir '+ds_dir)

obs = ds_list.keys()[0]
ds = ds_list[obs]

close('all')

### PLOT INDIVIDUAL OBSERVATIONS ###

srclist = [('ADLeo','13A-423'),('ADLeo','15A-416'),('UVCet','13A-423'),('UVCet','15A-416'),('EQPeg','15A-416'),('EVLac','15A-416'),('YZCMi','15A-416')]
srclist = [('YZCMi','15A-416')]

figure(figsize=(6,8.5))

n_rows = 5
n_cols = 4
gs = gridspec.GridSpec(n_rows, n_cols, width_ratios=[1,1,1.2,1])

for src,proj in srclist:

    figname = ds_dir + src + '_' + proj + '.pdf'
    clf()
    sub_row = 1
    subplots_adjust(hspace=1.0)

    cmd = 'ls -d /data/jrv/'+proj+'/'+src+'/*/*/*.tbavg.ms.dynspec | cut -d / -f 1-6 | uniq'
    obslist = subprocess.check_output(cmd, shell=True).rstrip().split('\n')
    
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
        pp = {'pol':'i','smin':smin,'smax':smax,'trim_mask':False,'axis_labels':['ylabel','cbar']}
        ds.plot_dynspec(plot_params=pp)
        cb = gca().images[-1].colorbar
        cb.remove()
        title('Stokes I')

        # plot Stokes I imaginary
        i = offset + 1
        subplot(gs[i])
        pp.update({'axis_labels':['cbar'],'func':imag})
        ds.plot_dynspec(plot_params=pp)
        gca().yaxis.set_visible(False)
        title('Imag(I)')

        # flux limits for Stokes V
        smax = ds.get_rms('v')*6
        smin = -smax  # make colorbar symmetric about zero to be consistent with Stokes V

        # plot Stokes V real
        i = offset + 2
        subplot(gs[i])
        pp = {'pol':'v','smin':smin,'smax':smax,'trim_mask':False,'axis_labels':['cbar','xlabel']}
        ds.plot_dynspec(plot_params=pp)
        cb = gca().images[-1].colorbar
        cb.remove()
        gca().yaxis.set_visible(False)
        gca().xaxis.set_label_coords(-0.5,-0.3)
        title('Stokes V')

        # plot Stokes V imaginary
        i = offset + 3
        subplot(gs[i])
        pp.update({'axis_labels':['cbar','cbar_label'],'func':imag})
        ds.plot_dynspec(plot_params=pp)
        gca().yaxis.set_visible(False)
        title('Imag(V)')
        
        sub_row += 1

    suptitle(src[0:2]+' '+src[2:5]+' 20'+proj[0:2])
    savefig(figname,bbox_inches='tight')
    clf()
    subrow=1

