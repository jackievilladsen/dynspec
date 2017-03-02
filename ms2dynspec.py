'''

ms2dynspec.py

Goal: given an ms (w/ fully calibrated data in data column) and a model image, subtract the model and then produce a dynamic spectrum plot

'''

import os
from tasks import *

def get_nterms(model):
    if type(model) is str:
        return 1
    else:
        return len(model)


try:
    print 'Running ms2dynspec.py with vis='+vis+', model='+model+', dsdir='+dsdir
except:
    print ''
    pass

    # subtract model from visibilities
    delmod(vis=vis)
    clearcal(vis=vis,addmodel=True)
    nterms = get_nterms(model)
    ft(vis=vis,model=model,nterms=nterms,usescratch=True)
    uvsub(vis=vis)
    
    # create dsdir, directory where stuff will be saved
    if dsdir[-1]!='/':
        dsdir+='/'
    if os.path.exists(dsdir):
        rmtables(dsdir+'*')
        os.system('rm -rf '+dsdir)
    os.system('mkdir '+dsdir)
    
    # avg over all baselines and save dynspec to file
    tb = dsdir + 'tbavg.ms'
    try:
        tbavg(split,vis,tb,speed='fast',weight_mode='flat',datacolumn='corrected')
    except:
        print 'tbavg failed, probably b/c table', tb, 'is already open in the CASA cache - restart CASA to fix this problem'
