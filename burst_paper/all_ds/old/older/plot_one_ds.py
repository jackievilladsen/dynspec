'''
plot_all_ds.py - Load all Stokes I,V dynamic spectra from file produced by compile_ds.py, then produce
                 a series of figures showing all of the dynamic spectra - Re(vis) and Im(vis).
'''

from dynspec import load_dict
from pylab import *
import os

rmsfac=3

def get_obsname(obsfile):
    # take a file directory such as '/data/jrv/15A-416/YZCMi/1' and
    # convert to obs name such as '15A-416_YZCMi_1' and srcname 'YZCMi'
    names = obsfile.split('/')
    srcname = names[4]
    obsname = names[3]+'_'+names[4]+'_'+names[5]
    return obsname,srcname

params = {'legend.fontsize': 'small',
          'axes.titlesize': 'medium',
          'axes.labelsize': 'small',
          'xtick.labelsize': 'x-small',
          'ytick.labelsize': 'x-small',
          'image.interpolation': 'hanning'}
rcParams.update(params)

savefile = '/data/jrv/burst_paper/all_dynspec.npy'
ds_list = load_dict(savefile)

ds_dir = '/data/jrv/burst_paper/ds/' # where to save individual ds plots
if not os.path.exists(ds_dir):
    os.system('mkdir '+ds_dir)

obs = ds_list.keys()[0]
ds = ds_list[obs]

close('all')

### PLOT INDIVIDUAL OBSERVATIONS ###

figure(figsize=(6,4))

for obs in ds_list:

    clf()
    print '\n', obs
    ds = ds_list[obs]

    # set up directory and filename to save plot
    obsname,srcname = get_obsname(obs)
    srcdir = ds_dir + srcname + '/'
    if not os.path.exists(srcdir):
        os.system('mkdir '+srcdir)
    plotfile = srcdir + obsname + '.pdf'

    # mask dynspec pixels with SNR < 3  (uses im(vis) to measure RMS --> minimal contamination from burst)
    ds_real_masked = ds.mask_SNR(rmsfac=rmsfac)
    ds_im_masked = ds.mask_SNR(rmsfac=rmsfac,func=imag)
    
    # plot Stokes I real and imaginary
    smin = 0
    smax = ds.get_rms('i')*6
    smin = -smax
    subplot(241)
    pp = {'pol':'i','smin':smin,'smax':smax,'trim_mask':False,'axis_labels':['ylabel','cbar']}
    ds.plot_dynspec(plot_params=pp)
    cb = gca().images[-1].colorbar
    cb.remove()
    title('Stokes I')
    subplot(242)
    pp['axis_labels']=['cbar']
    ds_real_masked.plot_dynspec(plot_params=pp)
    gca().yaxis.set_visible(False)
    cb = gca().images[-1].colorbar
    cb.remove()
    title('Stokes I SNR>3')
    subplot(243)
    pp['func'] = imag
    pp['axis_labels'] = ['cbar']
    ds.plot_dynspec(plot_params=pp)
    gca().yaxis.set_visible(False)
    cb = gca().images[-1].colorbar
    cb.remove()
    title('Imag(I)')
    subplot(244)
    pp['func'] = imag
    pp['axis_labels'] = ['cbar','cbar_label']
    ds_im_masked.plot_dynspec(plot_params=pp)
    gca().yaxis.set_visible(False)
    title('Imag(I) SNR>3')

    # plot Stokes V real and imaginary
    smax = ds.get_rms('v')*6
    smin = -smax
    pp = {'pol':'v','smin':smin,'smax':smax,'trim_mask':False,'axis_labels':['ylabel','cbar']}
    subplot(245)
    ds.plot_dynspec(plot_params=pp)
    cb = gca().images[-1].colorbar
    cb.remove()
    title('Stokes V')
    subplot(246)
    pp['axis_labels'] = ['xlabel','cbar']
    ds_real_masked.plot_dynspec(plot_params=pp)
    cb = gca().images[-1].colorbar
    cb.remove()
    gca().yaxis.set_visible(False)
    title('Stokes V SNR>3')
    subplot(247)
    pp['func'] = imag
    pp['axis_labels'] = ['cbar']
    ds.plot_dynspec(plot_params=pp)
    gca().yaxis.set_visible(False)
    cb = gca().images[-1].colorbar
    cb.remove()
    title('Imag(V)')
    subplot(248)
    pp['func'] = imag
    pp['axis_labels'] = ['cbar','cbar_label']
    ds_im_masked.plot_dynspec(plot_params=pp)
    gca().yaxis.set_visible(False)
    title('Imag(V) SNR>3')

    # save figure
    #subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.2, hspace=0.2)
    #tight_layout()
    savefig(plotfile,bbox_inches='tight')

