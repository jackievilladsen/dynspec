"""

yzcmi1_plots.py - script to plot properties of YZ CMi 1 burst - rc vs. freq, tseries to check for harmonic structure...

"""

import dynspec.plot
reload(dynspec.plot)
from dynspec.plot import *
from pylab import *

# tested spline interpolation to smooth tseries
# From time series z (a masked array):
#    t = s.get_tlist()
#    ind = find(~s.mask)
#    ti = t[ind]
#    zi = z[ind]
#    us = scipy.interpolate.UnivariateSpline(ti,zi,s=0.5) #0.48 also works
#    z_smooth=us(ti)
# This gives us a pretty good smoothed version of the time series

def get_band(filename):
    # get band from dynspec filename
    return filename.split('/')[6]

params = {'legend.fontsize': 'small',
          'axes.titlesize': 'small',
          'axes.labelsize': 'small',
          'xtick.labelsize': 'x-small',
          'ytick.labelsize': 'x-small',
          'image.interpolation': 'none'}
mpl.rcParams.update(params)

savedir='/data/jrv/burst_paper/ds/'

nt = 5
nf = 32

ds_files = ['/data/jrv/15A-416/UVCet/3/L/sc/ds_bgsub_postcal_noflag/tbavg.ms.dynspec',
'/data/jrv/15A-416/UVCet/3/S/sc/ds_bgsub_postcal/tbavg.ms.dynspec']
src = 'UVCet'
'''
ds_obs = None
for f in ds_files:
    band = get_band(f)
    params={'filename':f,'uniform':True}
    ds = Dynspec(params)
    ds.spec['i'] = (ds.spec['rr']+ds.spec['ll'])/2
    ds.spec['v'] = (ds.spec['rr']-ds.spec['ll'])/2
    del ds.spec['ll']
    del ds.spec['rr']
    if ds_obs is None:
        ds_obs = deepcopy(ds)
    else:
        ds_obs.add_dynspec(ds)
    del ds

ds_obs2=deepcopy(ds_obs)
ds_obs2.mask_RFI_pixels(rmsfac=20.)
ds_obs2.mask_RFI(rmsfac=5.)
'''
# plot dynspec to file for paper
nt=15
nf=12
ds_bin=ds_obs2.bin_dynspec(nt=nt,nf=nf,mask_partial=0.9)
smax=0.07
pp = {'tlims':[20,120],'pol':'v','trim_mask':True,'dy':0.5,'dx':20.,'ar0':'auto','axis_labels':['cbar','cbar_label','xlabel','ylabel'],'scale':'linear','smin':-smax,'smax':smax}

close('all')
figure(figsize=(6.8,3.75))
axes([0.07,0.1,0.875,0.88])

ds_bin.plot_dynspec(plot_params=pp)
gca().set_axis_bgcolor('w')
#xlabel('Time (min) since '+ds_bin.time[0].iso[:-2])
cbar = gca().images[-1].colorbar
cbar.ax.tick_params(labelsize='x-small')
savefig(savedir+'UVCet3_dynspec.pdf') 
