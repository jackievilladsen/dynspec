'''

ms2dynspec.py

Goal: given an ms (w/ fully calibrated data in data column) and a model image, subtract the model and then produce a dynamic spectrum plot

for casa-prereleave version 5.0.0-141, ft does not work, but this will be fixed in next subversion of prerelease (142?)

'''

import os
from pylab import *
from taskinit import tbtool
try:
    from tasks import *
except:
    print 'You are trying to use ms2dynspec.py outside CASA, but it needs the CASA tasks module.'

from extract_dynspec import saveTxt
import tbavg
reload(tbavg)

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

def tbavg2dsfile(tbavg_ms_file,overwrite_mode='backup'):
    '''
    dsdir = tbavg2dsfile(tbavg_ms_file,overwrite_mode='backup')
    
    Read dynamic spectrum from tbavg_ms_file and save it as a set
    of text files in directory dsdir (=tbavg_ms_file+'.dynspec').
    
    overwrite_mode = 'backup' means that if there is an old dynspec
    file with the same name as dsfile, it is moved to [dsdir].old.
    Other options:
    - overwrite_mode = 'overwrite': just erase old dynspec file
    - overwrite_mode = 'none' (or anything else): return w/o creating new ds file
      if there is already an existing file
    '''
    # name of directory to put dynamic spectra text files in
    if tbavg_ms_file[-1] =='/':
        tbavg_ms_file = tbavg_ms_file[:-1]
    dsdir = tbavg_ms_file + '.dynspec'
    
    # check if directory already exists and if so, back it up, erase it,
    #   or exit w/o writing a new file (depending on overwrite_mode)
    if os.path.exists(dsdir):
        print 'Warning from tbavg2dsfile:',dsdir, 'already exists'
        if overwrite_mode=='backup':
            if os.path.exists(dsdir+'.old'):
                os.system('rm -rf '+dsdir+'.old')
            os.system('mv '+dsdir+' '+dsdir+'.old')
            print 'tbavg2dsfile is moving old',dsdir, 'to', dsdir+'.old before writing new dynspec to',dsdir, '(potentially overwriting old backup)'
        elif overwrite_mode=='overwrite':
            os.system('rm -rf '+dsdir)
            print 'tbavg2dsfile is erasing old',dsdir, 'before creating the new one'
        else:
            print 'tbavg2dsfile is not writing a new copy of',dsdir,'- returning old version'
            return dsdir
    
    # Read dynamic spectrum from tbavg.ms file
    from tbavg import dyn_spec
    spec = dyn_spec(tbavg_ms_file)
    
    # Save dynamic spectrum to text files in dsdir
    #  (name of dsdir is determined automatically by saveTxt)
    saveTxt(spec,tbavg_ms_file)
    return dsdir

def make_dsfile(vis,dsdir,datacolumn='corrected',weight_mode='flat'):
    # take post-uvsub ms and run tbavg, generate dsfile
        # create dsdir, directory where stuff will be saved
    
    from tbavg import tbavg, dyn_spec
    
    if dsdir[-1]!='/':
        dsdir+='/'
    tb = dsdir + 'tbavg.ms'
    if os.path.exists(dsdir):
        print 'Warning: ms2dsfile is now deleting and replacing the contents of', dsdir,'- remember not to try to create same ms filename w/ tbavg twice without restarting CASA'
        os.system('rm -rf '+dsdir)
    os.system('mkdir '+dsdir)
    
    # avg over all baselines and save dynspec to file
    try:
        tbavg(vis,tb,speed='fast',weight_mode=weight_mode,datacolumn=datacolumn)
    except:
        print 'tbavg failed, probably b/c table', tb, 'is already open in the CASA cache - restart CASA to fix this problem'
    spec = dyn_spec(tb)
    dsfile = tb+'.dynspec'
    if os.path.exists(dsfile):
        os.system('rm -rf '+dsfile)
    saveTxt(spec,tb)
    
    return dsfile

def copy_model_RRtoLL(vis):
    # copy RR column of model_data column to LL
    tb = tbtool()
    tb.open(vis,nomodify=False)
    model_vis = tb.getcol('MODEL_DATA')
    model_vis[3,:,:] = model_vis[0,:,:] # copy RR model column to LL model column
    tb.putcol('MODEL_DATA',model_vis)
    tb.unlock()
    tb.close()

def copy_data_RRtoLL(vis):
    # copy RR column of data column to LL
    # doesn't account for different flags on RR and LL
    tb = tbtool()
    tb.open(vis,nomodify=False)
    data_vis = tb.getcol('DATA')
    data_vis[3,:,:] = data_vis[0,:,:] # copy RR model column to LL model column
    tb.putcol('DATA',data_vis)
    tb.unlock()
    tb.close()

def mask_cen(model,pixel_rad=5):
    # create new model image(s) with circle of radius pixel_rad masked in center
    
    # if model image is a string (e.g., nterms=1), put it in a list so that we can iterate
    nterms = get_nterms(model)
    if nterms == 1:
        model = [model]
    
    # create mask with ones in center, which we can apply to model image
    cen_mask = 'center_mask'
    if os.path.exists(cen_mask):
        rmtables(cen_mask)
    ncen = str(int(imhead(imagename=model[0],mode='get',hdkey='crpix1')))   # center pixel of image (assuming square)
    rad = str(int(pixel_rad))
    circle_cen = 'circle[['+ncen+'pix,'+ncen+'pix],'+rad+'pix]' # describes a region of a circle in the image center with radius pixel_rad
    makemask(mode='copy',inpimage=model[0],output=cen_mask,inpmask=circle_cen,overwrite=True) # create mask with ones in center, of same size as model image
    
    # create new model images with center pixels masked out
    bgmodel = []
    for modelim in model:
        bgmodelim = modelim + '.cen_masked'
        immath(imagename=[modelim,cen_mask],outfile=bgmodelim,expr='IM0*(1-IM1)')
        bgmodel.append(bgmodelim)
    
    # if model image was a string, make bgmodel a string too
    if nterms == 1:
        bgmodel = bgmodel[0]
    
    return bgmodel

def ft_multi_im(vis,model_list,usescratch=True):
    '''
    ft_multi_im will ft multiple model images into the model_data column of vis.
    model_list should be of the form:
       nterms=1:  model_list = ['im1.model','im2.model','im3.model',...]
       nterms>1:  model_list = [['im1.model.tt0','im1.model.tt1'],['im2.model.tt0','im2.model.tt1'],...]
    model_list must be a list
    '''
    incremental=False
    for model in model_list:
        nterms = get_nterms(model)
        ft(vis=vis,model=model,nterms=nterms,usescratch=usescratch,incremental=incremental)
        incremental=True


def ms2dsfile(vis,model=[],dsdir='tbavg',reset_corrected=True,pop_model=True,weight_mode='flat'):    
    '''
    ms2dsfile subtracts model from vis, then runs tbavg and extracts the dynspec to a numpy data file.
    Returns location of dsfile.
    '''
    
    if model==[]:
        pop_model=False
    
    # subtract model from visibilities
    if reset_corrected:
        delmod(vis=vis)
        clearcal(vis=vis,addmodel=True)
    if pop_model:
        nterms = get_nterms(model)
        ft(vis=vis,model=model,nterms=nterms,usescratch=True)
        #ft(vis=vis,model=model,nterms=nterms,usescratch=True)
    uvsub(vis=vis)

    return make_dsfile(vis,dsdir,weight_mode=weight_mode)
    
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

def dsplot_Pband_fullstokes(dsfile,nt=50,nf=32,smax='auto',rmsfac=1.5):
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
    subplot(421)
    pp = {'pol':'i','smin':smin,'smax':smax,'trim_mask':False,'axis_labels':['ylabel','cbar'],'dy':0.05}
    ds.plot_dynspec(plot_params=pp)
    cb = gca().images[-1].colorbar
    #cb.remove()
    title('Stokes I')
    
    subplot(422)
    pp.update({'axis_labels':['cbar'],'func':imag})
    ds.plot_dynspec(plot_params=pp)
    gca().yaxis.set_visible(False)
    title('Imag(I)')
    
    subplot(423)
    pp = {'pol':'v','smin':smin,'smax':smax,'trim_mask':False,'axis_labels':['cbar','xlabel']}
    ds.plot_dynspec(plot_params=pp)
    cb = gca().images[-1].colorbar
    gca().yaxis.set_visible(False)
    gca().xaxis.set_label_coords(-0.5,-0.3)
    title('Stokes V')
    
    subplot(424)
    pp.update({'axis_labels':['cbar','cbar_label'],'func':imag})
    ds.plot_dynspec(plot_params=pp)
    gca().yaxis.set_visible(False)
    title('Imag(V)')

    subplot(425)
    pp = {'pol':'q','smin':smin,'smax':smax,'trim_mask':False,'axis_labels':['ylabel','cbar'],'dy':0.05}
    ds.plot_dynspec(plot_params=pp)
    cb = gca().images[-1].colorbar
    #cb.remove()
    title('Stokes Q')
    
    subplot(426)
    pp.update({'axis_labels':['cbar'],'func':imag})
    ds.plot_dynspec(plot_params=pp)
    gca().yaxis.set_visible(False)
    title('Imag(Q)')
    
    subplot(427)
    pp = {'pol':'u','smin':smin,'smax':smax,'trim_mask':False,'axis_labels':['cbar','xlabel']}
    ds.plot_dynspec(plot_params=pp)
    cb = gca().images[-1].colorbar
    #cb.remove()
    gca().yaxis.set_visible(False)
    gca().xaxis.set_label_coords(-0.5,-0.3)
    title('Stokes U')
    
    subplot(428)
    pp.update({'axis_labels':['cbar','cbar_label'],'func':imag})
    ds.plot_dynspec(plot_params=pp)
    gca().yaxis.set_visible(False)
    title('Imag(U)')
    
    savefig(dsfile+'/dynspec.pdf',bbox_inches='tight')
    
    return ds
