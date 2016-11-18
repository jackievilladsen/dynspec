'''
plot_all_ds.py - Load all Stokes I,V dynamic spectra from file produced by compile_ds.py, then produce
                 a series of figures showing all of the dynamic spectra - Re(vis) and Im(vis).
'''

from dynspec import load_dict
from dynspec.plot import Dynspec
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

savefile = '/data/jrv/burst_paper/ds/all_Pdynspec.npy'
ds_list = load_dict(savefile)

ds_dir = '/data/jrv/burst_paper/ds/' # where to save individual ds plots
if not os.path.exists(ds_dir):
    os.system('mkdir '+ds_dir)

obs = ds_list.keys()[0]
ds = ds_list[obs]

close('all')

### PLOT INDIVIDUAL OBSERVATIONS ###
'''
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
    plotfile = srcdir + obsname + '_P.pdf'

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
'''
### FOLLOW-UP PLOTS FOR INTERESTING BURSTS ###

close('all')

figure(figsize=(6,4))

pp_dict = {'ADLeo/3':{'nt':60,'nf':128},
           'ADLeo/5':{'tlims':array([170.,210.]),'nt':3,'nf':32,'smin':-.1,'smax':.1},
           'UVCet/3':{'nt':20,'nf':64,'smin':-.05,'smax':.05},
           'UVCet/4':{'tlims':array([100,160])},
           'UVCet/5':{'nt':30,'nf':64}}

pp_dict = {'EQPeg/1':{},
           'EQPeg/2':{}}

nt_def = 10
nf_def = 64

for s in pp_dict:
    pp = pp_dict[s]
    srcname, obsnum = s.split('/')
    ds_file = '/data/jrv/15A-416/'+s+'/P/'+srcname+'_'+obsnum+'P.tbavg.ms.dynspec'
    params={'filename':ds_file,'uniform':True}
    
    nt = pp.get('nt',nt_def)
    nf = pp.get('nf',nf_def)
    
    ds = Dynspec(params)
    ds.spec['i'] = (ds.spec['xx']+ds.spec['yy'])/2
    ds.spec['v'] = (ds.spec['xy']-ds.spec['yx'])/(2.j)
    ds.spec['u'] = (ds.spec['xy']+ds.spec['yx'])/2
    del ds.spec['xx']
    del ds.spec['yy']
    del ds.spec['xy']
    del ds.spec['yx']
    ds = ds.bin_dynspec(nt=nt,nf=nf)

    figfile = ds_dir + srcname + '_' + obsnum + 'P.pdf'
    
    smax = min(ds.rms_spec('v'))*5
    smin = -smax

    i_pp = {'pol':'i','smin':smin,'smax':smax,'trim_mask':False,'axis_labels':['ylabel','xlabel','cbar']}
    u_pp = {'pol':'u','smin':smin,'smax':smax,'trim_mask':False,'axis_labels':['cbar']}
    v_pp = {'pol':'v','smin':smin,'smax':smax,'trim_mask':False,'axis_labels':['cbar','cbar_label']}
    i_pp.update(pp)
    u_pp.update(pp)
    v_pp.update(pp)
    
    clf()
    subplot(131)
    ds.plot_dynspec(plot_params=i_pp)
    cb = gca().images[-1].colorbar
    cb.remove()
    gca().xaxis.set_label_coords(1.2,-0.15)
    title('Stokes I')
    subplot(132)
    ds.plot_dynspec(plot_params=u_pp)
    cb = gca().images[-1].colorbar
    cb.remove()
    title('Stokes U')
    subplot(133)
    ds.plot_dynspec(plot_params=v_pp)
    title('Stokes V')
    suptitle(srcname)
    savefig(figfile,bbox_inches='tight')
