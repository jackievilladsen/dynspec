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

savedir='/data/jrv/burst_paper/yzcmi/'

nt = 5
nf = 32

ds_files = ['/data/jrv/15A-416/YZCMi/1/L/YZCMi_1L.tbavg.ms.dynspec', '/data/jrv/15A-416/YZCMi/1/S/YZCMi_1S.tbavg.ms.dynspec']
src = 'YZCMi'
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
ds_obs.mask_RFI(rmsfac=5.) # so far have only plotted YZ CMi dynspec w/ this - only makes a minor diff - maybe bigger when there is less binning?

# plot dynspec to narrow down on desired time range
ds_bin = ds_obs.bin_dynspec(nt=nt,nf=nf)
pp = {'tlims':[61.5,63.5],'pol':'ll','smin':0,'smax':0.03,'trim_mask':False,'dy':0.5,'dx':0.5}
clf()
ds_bin.plot_dynspec(plot_params=pp)
show()

ds_bin = ds_obs.bin_dynspec(nt=nt,nf=1)
tseries = ds_bin.tseries(mask_partial=0.8)
t = tseries.get_tlist()*tseries.dt()/60.

clf()
plot(t,real(tseries.spec['i']))
show()

#ds_bin = ds_obs.bin_dynspec(nt=1,nf=nf)
#spec = ds_bin.make_spectrum(tmin=61.5,tmax=63.5,trim_mask=False,mask_partial=0.8)
spec = ds_bin.make_spectrum(tmin=62.14,tmax=64.14,trim_mask=False,mask_partial=0.8)
rr = spec.spec['rr'] * 1.e3
ll = spec.spec['ll'] * 1.e3
#rc = real(rr-ll)/real(rr+ll)
rc = real(spec.spec['v'])/real(spec.spec['i'])
f = spec.f/1.e9

close('all')
figure(figsize=(3,4))
subplot(211)
plot(f,real(rr),'k',linewidth=1.8)
plot(f,real(ll),'k')
#plot(f,imag(rr),'b--')
#plot(f,imag(ll),'g--')
xlim([1,4])
legend(('RCP','LCP'))
gca().axhline(0,color='k')
ylabel('Time-averaged flux (mJy)')
#axis([1,4,-18,18])

subplot(212)
plot(f,rc,'k.')
gca().axhline(0,color='k')
gca().axhline(1,color='k',linestyle='--')
gca().axhline(-1,color='k',linestyle='--')
axis([1,4,-1.4,1.4])
xlabel('Frequency (GHz)')
ylabel('Degree of circular pol.')

print spec.time.iso

savefig(savedir+'YZCMi_burstspec.pdf',bbox_inches='tight')

# now plot tseries of RCP to check for harmonics

ds_bin2 = ds_obs.bin_dynspec(nt=nt,nf=1)
ds_bin=ds_bin2.clip(tmin=50.,tmax=80.,trim_mask=False)
tseries1 = ds_bin.tseries(fmin=1.e9,fmax=1.4e9,mask_partial=0.8)
tseries2 = ds_bin.tseries(fmin=2.e9,fmax=2.8e9,mask_partial=0.8)
tseries3 = ds_bin.tseries(fmin=1.5e9,fmax=2.e9,mask_partial=0.8)
tseries4 = ds_bin.tseries(fmin=3.e9,fmax=4.e9,mask_partial=0.8)
t = tseries1.get_tlist()*tseries1.dt()/60.

close('all')
figure(figsize=(3,4))
subplot(211)
plot(t,real(tseries1.spec['rr']),'k')
plot(t,real(tseries2.spec['rr']),'r')
legend(('1-1.4 GHz','2-2.8 GHz'))
ylabel('RCP Flux (mJy)')
subplot(212)
plot(t,real(tseries3.spec['ll']),'k')
plot(t,real(tseries4.spec['ll']),'r')
legend(('1.5-2 GHz','3-4 GHz'))
xlabel('Time (min) since '+tseries1.time[0].iso[:-2])
ylabel('LCP Flux (mJy)')
savefig(savedir+'YZCMi_tseries.pdf',bbox_inches='tight')
'''


# now plot dynspec to file for paper
nt=4
nf=16
ds_bin=ds_obs.bin_dynspec(nt=nt,nf=nf,mask_partial=0.8)
smax=0.08
pp = {'tlims':[58.,73.],'pol':'v','smin':-smax,'smax':smax,'trim_mask':False,'dy':0.5,'dx':5}

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
savefig(savedir+'YZCMi1_dynspec.pdf',bbox_inches='tight')
