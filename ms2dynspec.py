'''

ms2dynspec.py

Goal: given an ms (w/ fully calibrated data in data column) and a model image, subtract the model and then produce a dynamic spectrum plot

for casa-prereleave version 5.0.0-141, ft does not work, but this will be fixed in next subversion of prerelease (142?)

'''

import os
from tasks import *
from tbavg import tbavg, dyn_spec
from extract_dynspec import saveTxt

def get_nterms(model):
    if type(model) is str:
        return 1
    else:
        return len(model)

def ms2dsfile(vis,model,dsdir):
    '''
    ms2dsfile subtracts model from vis, then runs tbavg and extracts the dynspec to a numpy data file
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
    if os.path.exists(tb+'.dynspec'):
        os.system('rm -rf '+tb+'.dynspec')
    saveTxt(spec,tb)

    
def ms2dsplot(vis,model,dsdir):
    # go from vis and model to dynspec plot
    # the binning and plotting parameters are specific to the case of a P-band example data set
    
    return 'hi'
