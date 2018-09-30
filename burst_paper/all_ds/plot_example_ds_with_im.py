

###########

# now plot just an example of ADLeo_3 with real and imaginary components

###########

import dynspec.plot
reload(dynspec.plot)

from dynspec import load_dict
from dynspec.plot import *
from pylab import *
import os, subprocess
import matplotlib.gridspec as gridspec


params = {'legend.fontsize': 'small',
          'axes.titlesize': 'small',
          'axes.labelsize': 'x-small',
          'xtick.labelsize': 'xx-small',
          'ytick.labelsize': 'xx-small',
          'image.interpolation': 'none'}
rcParams.update(params)

savefile = '/data/jrv/burst_paper/all_burst_dynspec.npy'
savedir='/data/jrv/burst_paper/adleo/'

nt = 60
nf = 32

ds_files = ['/data/jrv/15A-416/ADLeo/3/L/test_clean/ds_ap0_big_RR_n2_ms/tbavg.ms.dynspec','/data/jrv/15A-416/ADLeo/3/S/test_selfcal/ds_ap1_n3_bgsub_big/tbavg.ms.dynspec']
src = 'ADLeo'
'''
ds_obs = None
for f in ds_files:
    params={'filename':f,'uniform':True}
    ds = Dynspec(params)
    ds.spec['i'] = (ds.spec['rr']+ds.spec['ll'])/2
    ds.spec['v'] = (ds.spec['rr']-ds.spec['ll'])/2
    if ds_obs is None:
        ds_obs = deepcopy(ds)
    else:
        ds_obs.add_dynspec(ds)
    del ds
#ds_obs.mask_RFI(rmsfac=5.)

ds=ds_obs.bin_dynspec(nt=nt,nf=nf,mask_partial=0.9)
'''
figure(figsize=(6.5,6.5))

n_rows = 2
n_cols = 2
gs = gridspec.GridSpec(n_rows, n_cols)
ar0 = 1.0

clf()
sub_row = 1
subplots_adjust(hspace=0.2,wspace=0)
        
offset = 0

# flux limits
smax = percentile(real(ds.spec['i']),98)
smax = 0.015
smin = -smax  # make colorbar symmetric about zero to be consistent with Stokes V

# plot Stokes I real
i = offset + 0
subplot(gs[i])
pp = {'pol':'i','smin':smin,'smax':smax,'trim_mask':False,'axis_labels':['cbar','ylabel'],'ar0':ar0,'dy':0.5}
ds.plot_dynspec(plot_params=pp)
cb = gca().images[-1].colorbar
cb.remove()
gca().yaxis.set_label_coords(-0.05,-0.1)
title('Stokes I')

# plot Stokes I imag
i = offset + 1
subplot(gs[i])
pp = {'pol':'i','func':imag,'smin':smin,'smax':smax,'trim_mask':False,'axis_labels':['cbar'],'ar0':ar0}
ds.plot_dynspec(plot_params=pp)
gca().yaxis.set_visible(False)
cb = gca().images[-1].colorbar
cb.remove()
title('Imag(I)')

# plot Stokes V real
i = offset + 2
subplot(gs[i])
pp = {'pol':'v','smin':smin,'smax':smax,'trim_mask':False,'axis_labels':['cbar','xlabel'],'ar0':ar0,'dy':0.5}
ds.plot_dynspec(plot_params=pp)
#gca().yaxis.set_visible(False)
cb = gca().images[-1].colorbar
cb.remove()
gca().xaxis.set_label_coords(1.0,-0.1)
title('Stokes V')

# plot Stokes V imag
i = offset + 3
subplot(gs[i])
pp = {'pol':'v','func':imag,'smin':smin,'smax':smax,'trim_mask':False,'axis_labels':['cbar','cbar_label'],'ar0':ar0}
ds.plot_dynspec(plot_params=pp)
gca().yaxis.set_visible(False)
title('Imag(V)')
        
savefig(savedir+'ADLeo3_example_ds.pdf',bbox_inches='tight')

