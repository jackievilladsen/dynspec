'''

ms2dynspec.py

Goal: given an ms (w/ fully calibrated data in data column) and a model image, subtract the model and then produce a dynamic spectrum plot

for casa-prereleave version 5.0.0-141, ft does not work, but this will be fixed in next subversion of prerelease (142?)

'''

import os
from tbavg import tbavg, dyn_spec
from extract_dynspec import saveTxt
from pylab import *

params = {'legend.fontsize': 'small',
          'axes.titlesize': 'small',
          'axes.labelsize': 'x-small',
          'xtick.labelsize': 'xx-small',
          'ytick.labelsize': 'xx-small',
          'image.interpolation': 'hanning'}
mpl.rcParams.update(params)

def get_nterms(model):
    if type(model) is str:
        return 1
    else:
        return len(model)

def ms2dsfile(vis,model,dsdir):
    from tasks import *
    
    '''
    ms2dsfile subtracts model from vis, then runs tbavg and extracts the dynspec to a numpy data file.
    Returns location of dsfile.
    '''

    # subtract model from visibilities
    delmod(vis=vis)
    clearcal(vis=vis,addmodel=True)
    nterms = get_nterms(model)
    ft(vis=vis,model=model,nterms=nterms,usescratch=True)
    ft(vis=vis,model=model,nterms=nterms,usescratch=True)
    uvsub(vis=vis)

    # create dsdir, directory where stuff will be saved
    if dsdir[-1]!='/':
        dsdir+='/'
    tb = dsdir + 'tbavg.ms'
    if os.path.exists(dsdir):
        print 'Warning: ms2dsfile is now deleting and replacing the contents of', dsdir,'- remember not to try to create same ms filename w/ tbavg twice without restarting CASA'
        os.system('rm -rf '+dsdir)
    os.system('mkdir '+dsdir)
    
    # avg over all baselines and save dynspec to file
    try:
        tbavg(vis,tb,speed='fast',weight_mode='flat',datacolumn='corrected')
    except:
        print 'tbavg failed, probably b/c table', tb, 'is already open in the CASA cache - restart CASA to fix this problem'
    spec = dyn_spec(tb)
    dsfile = tb+'.dynspec'
    if os.path.exists(dsfile):
        os.system('rm -rf '+dsfile)
    saveTxt(spec,tb)
    
    return dsfile

    
def dsplot_Pband(dsfile,nt=50,nf=32,smax='auto',rmsfac=1.5):
    # go from vis and model to dynspec plot
    # the binning and plotting parameters are specific to the case of a P-band example data set
    
    from plot import Dynspec
        
    ds=Dynspec({'filename':dsfile,'uniform':True,'convert_stokes':True})
    ds.mask_RFI(rmsfac=rmsfac)
    ds = ds.bin_dynspec(nt=nt,nf=nf)
    
    if smax=='auto':
        smax = ds.get_rms('i') * 3
    smin = -smax
    
    figure()
    subplot(221)
    pp = {'pol':'i','smin':smin,'smax':smax,'trim_mask':False,'axis_labels':['ylabel','cbar'],'dy':0.05}
    ds.plot_dynspec(plot_params=pp)
    cb = gca().images[-1].colorbar
    #cb.remove()
    title('Stokes I')
    
    subplot(222)
    pp.update({'axis_labels':['cbar'],'func':imag})
    ds.plot_dynspec(plot_params=pp)
    gca().yaxis.set_visible(False)
    title('Imag(I)')
    
    subplot(223)
    pp = {'pol':'v','smin':smin,'smax':smax,'trim_mask':False,'axis_labels':['cbar','xlabel']}
    ds.plot_dynspec(plot_params=pp)
    cb = gca().images[-1].colorbar
    #cb.remove()
    gca().yaxis.set_visible(False)
    gca().xaxis.set_label_coords(-0.5,-0.3)
    title('Stokes V')
    
    subplot(224)
    pp.update({'axis_labels':['cbar','cbar_label'],'func':imag})
    ds.plot_dynspec(plot_params=pp)
    gca().yaxis.set_visible(False)
    title('Imag(V)')
    
    savefig(dsfile+'/dynspec.pdf',bbox_inches='tight')
    
    return ds
