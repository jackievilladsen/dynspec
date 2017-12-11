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

params = {'legend.fontsize': 'x-small',
          'axes.titlesize': 'small',
          'axes.labelsize': 'x-small',
          'xtick.labelsize': 'xx-small',
          'ytick.labelsize': 'xx-small',
          'image.interpolation': 'hanning'}
mpl.rcParams.update(params)

savedir='/data/jrv/burst_paper/uvcet/'

nt = 5
nf = 32

ds_files = ['/data/jrv/15A-416/UVCet/3/L/UVCet_3L.tbavg.ms.dynspec',
            '/data/jrv/15A-416/UVCet/3/S/UVCet_3S.tbavg.ms.dynspec']
src = 'UVCet'
'''
ds_obs = None
for f in ds_files:
    band = get_band(f)
    params={'filename':f,'uniform':True}
    ds = Dynspec(params)
    ds.spec['i'] = (ds.spec['rr']+ds.spec['ll'])/2
    ds.spec['v'] = (ds.spec['rr']-ds.spec['ll'])/2
    if ds_obs is None:
        ds_obs = deepcopy(ds)
    else:
        ds_obs.add_dynspec(ds)
    del ds
ds_obs.mask_RFI(rmsfac=5.)
'''
# plot dynspec to file for paper
nt=15
nf=8
ds_bin=ds_obs.bin_dynspec(nt=nt,nf=nf,mask_partial=0.8)
smax=0.1
pp = {'tlims':[20,120],'pol':'v','smin':-smax,'smax':smax,'trim_mask':False,'dy':0.5,'dx':60.}

params = {'legend.fontsize': 'x-small',
          'axes.titlesize': 'x-small',
          'axes.labelsize': 'xx-small',
          'xtick.labelsize': 'xx-small',
          'ytick.labelsize': 'xx-small',
          'image.interpolation': 'none'}
mpl.rcParams.update(params)

close('all')
figure(figsize=(3,2))
ds_bin.plot_dynspec(plot_params=pp)
gca().set_axis_bgcolor('w')
#xlabel('Time (min) since '+ds_bin.time[0].iso[:-2])
cbar = gca().images[-1].colorbar
cbar.ax.tick_params(labelsize='xx-small')
savefig(savedir+'UVCet3_dynspec.pdf',bbox_inches='tight')
