# -*- coding: utf-8 -*-
"""
Last modified: August 7 2015

@author: jackie
"""

from pylab import *
from numpy import *
import matplotlib.pyplot as plt
import os
from astropy.time import Time
from copy import deepcopy
from dynspec.plot import *


##### USER-DEFINED PARAMETERS #####

# scaling factors for color plots
plot_scale_factor = 5

# enable or disable inset plot on binned dynamic spectrum figure
make_inset = False


mydir = os.getcwd()
tmp = mydir[1:].split('/')
proj = tmp[2]
src = tmp[3]
obs = tmp[4]
bands = ['L','S']

# parameter filelist is not currently used - will eventually be used to cycle through files to load
filelist={}
for band in bands:
    filename = band + '/' + src + '_' + obs + band + '.tbavg.ms.dynspec'
    filelist[band] = filename
print 'file list:', filelist

# font size for plots
rcParams.update({'font.size': 10})

# number of time and frequency channels to bin together
nt = 8
nf = 8

# properties for large dynspec plots (binned)
ar = 1.0 # aspect ratio setting for large dynspec (binned)
ar0 = ar

# properties for inset dynspec plots (unbinned for flare peak)
ar2 = ar * nt/nf # maintains same aspect ratio as large dynspec
l1min=log10(0.05) # different colorbar to highlight brightest features
l1max=log10(0.3)


##### END USER-DEFINED PARAMETERS #####
try:
    ds
except:
    ds = None
    for band in bands:
        params={'filename':filelist[band],'uniform':True}
        if band=='S' or band=='C':
            params['df']=2.e6
        ds_band = Dynspec(params)
        if band=='L':
            ds_band = ds_band.bin_dynspec(nt=1,nf=2) # bin L-band dynamic spectrum to 2-MHz channels
            ds_band.f -= 0.5e6   # this is a quick fix; figure out how to address this better later
        if ds is None:
            ds = deepcopy(ds_band)
        else:
            ds.add_dynspec(ds_band)
        ds_band = None

plot_params = {}
#plot_params = {'dx':30,'tlims':array([20.,120.]),'smin':0.01} # UVCet?
plot_params = {'flims':array([1.e9,2.1e9]),'dx':30} # ADLeo_4
#plot_params = {'tlims':array([60.,71.]),'dx':5,'smin':.030,'smax':.180} # YZCMi_1


nt = 4
nf = 4
ds_bin = ds.bin_dynspec(nt,nf)
ds_bin.mask_RFI(5.)

ds_bin.spec['v'] = ds_bin.spec['ll']-ds_bin.spec['rr']
pol_dict = {'ll':'LCP','rr':'RCP','v':'StokesV'}
for pol in ['ll','rr','v']:
    plot_params['pol'] = pol
    ds_bin.plot_dynspec(plot_params=plot_params)
    p = pol_dict[pol]
    title(p+' Dynamic Spectrum')
    plotfile = 'plots/'+p+'dynspec.png'
    if os.path.exists(plotfile):
        os.system('rm -f ' + plotfile)
    savefig(plotfile,bbox_inches='tight')

ll = ds_bin.spec['ll']
rms = std(ll,0) # rms for each frequency
tseries = average(ll,1,1/rms**2) # compute variance-weighted tseries (so RFI doesn't dominate)




'''
ds_bin.plot_dynspec(plot_params={'dx':30,'pol':'ll'})
title('UV Ceti LCP Dynamic Spectrum - 07/04/2015')
savefig('dynspecLCP.pdf',bbox_inches='tight')

ds_bin.spec['rc']=(ds_bin.spec['rr']-ds_bin.spec['ll'])/(ds_bin.spec['rr']+ds_bin.spec['ll'])
ds_bin.plot_dynspec(plot_params={'dx':30,'pol':'rc','scale':'linear','smin':-1.0,'smax':1.0})
title('UV Ceti Degree of Circular Poln - 07/04/2015')
savefig('rc.pdf',bbox_inches='tight')

nt=2
nf=2
ds_bin2 = ds_bin.bin_dynspec(nt,nf)
ar0 = 1.0
ds_bin2.plot_dynspec(plot_params={'dx':2,'tlims':array([30.5,40.]),'df':0.1})
title('UV Ceti RCP - Downwards Sweep')
savefig('shortburstup.pdf',bbox_inches='tight')

ds_bin2.plot_dynspec(plot_params={'dx':2,'tlims':array([50.,70.])})
title('UV Ceti RCP - Downwards Sweep')
savefig('shortburstdown.pdf',bbox_inches='tight')
close('all')







##### PLOT DYNAMIC SPECTRA #####



if make_inset:
    # clip out a portion of the dynspec to zoom in on
    Clims = [60,100,1,1.3]
    rrC,tC_min,fC_GHz=clip_dynspec(rr,Clims,t_min,f_GHz)
    llC,tC_min,fC_GHz=clip_dynspec(ll,Clims,t_min,f_GHz)
    
    # the commands below create an inset plot to highlight the flare peak or other features that are best
    #   displayed with different binning and/or a different color scale than the rest of the dynamic spectrum
    #
    # add yellow rectangle to large plot to highlight inset region
    ax.add_patch(Rectangle((imin,jmin),di,dj,edgecolor='#ffff00',fill=False,linewidth=0.5))
    # create inset plot
    ax=axes([0.5,0.5,0.25,0.3],axisbg='k')
    imshow(log10(rrC.T),aspect=ar2,vmin=l1min,vmax=l1max,origin='lower',cmap='seismic')
    title('Flare Peak',color='w',size=10)
    # add yellow rectangle around whole plot (matches rectangle marking inset region on large dynspec)
    ax.add_patch(Rectangle((0,0),shape(rrC)[0]-1,shape(rrC)[1]-1,edgecolor='#ffff00',fill=False,linewidth=2))
    # label axes
    tick_params(labelcolor='w',labelsize=8)
    tick_labels, tick_locs = make_tick_labels(arange(0.,80.,0.5),tC_min)
    xticks(tick_locs, tick_labels) # label x axis time in minutes
    tick_labels, tick_locs = make_tick_labels(arange(1.,4.,0.1),fC_GHz)
    yticks(tick_locs, tick_labels) # label y axis in GHz
    # add colorbar and change labels to show log scale
    cbar = colorbar()
    ticks = arange(l1min,l1max+0.25,0.25)
    cbar.set_ticks(ticks)
    ticklbls = (10**(ticks+3)).round().astype(int)
    cbar.set_ticklabels(ticklbls)
    cbar.ax.yaxis.set_tick_params(labelcolor='w',labelsize=8)

savefig('dynspecRCP.pdf',bbox_inches='tight')

### FIGURE: LCP DYNAMIC SPECTRUM ###
figure()

## Large plot (entire dynspec, binned) ##
ax=subplot(111,axisbg='k')
imshow(log10(ll_bin).T,aspect=ar,vmin=lmin,vmax=lmax,origin='lower',cmap='seismic')
title('UV Ceti Flare LCP Dynamic Spectrum')

# add yellow rectangle to highlight inset region
ax.add_patch(Rectangle((imin,jmin),di,dj,edgecolor='#ffff00',fill=False,linewidth=0.5))

# add colorbar and change labels to show log scale
cbar = colorbar()
cbar.set_label('LCP Flux (mJy)')
ticks = arange(lmin,lmax+0.5,0.5)
cbar.set_ticks(ticks)
ticklbls = (10**(ticks+3)).round(1)
cbar.set_ticklabels(ticklbls)

# label axes: x axis in minutes, y axis in GHz
tick_labels,tick_locs = make_tick_labels(xtick_lbls,t_bin/60.)
xticks(tick_locs, tick_labels)
xlabel('Time in minutes since '+tstring+' UTC')
tick_labels,tick_locs = make_tick_labels(ytick_lbls,f_bin/1.e9)
yticks(tick_locs, tick_labels)
ylabel('Frequency (GHz)')

## Inset plot (flare peak only, unbinned) ##
ax=axes([0.15,0.5,0.25,0.3],axisbg='k')
imshow(log10(llC.T),aspect=ar2,vmin=lmin,vmax=lmax,origin='lower',cmap='seismic')
title('Flare Peak',color='w',size=10)

# add yellow rectangle around whole plot (matches rectangle marking inset region on large dynspec)
ax.add_patch(Rectangle((0,0),shape(rrC)[0]-1,shape(rrC)[1]-1,edgecolor='#ffff00',fill=False,linewidth=2))

# label axes
tick_params(labelcolor='w',labelsize=8)
tick_labels, tick_locs = make_tick_labels(arange(0.,40.,0.5),tC_min)
xticks(tick_locs, tick_labels) # label x axis time in minutes
tick_labels, tick_locs = make_tick_labels(arange(4.,6.,0.1),fC_GHz)
yticks(tick_locs, tick_labels) # label y axis in GHz

# add colorbar and change labels to show log scale
cbar = colorbar()
ticks = arange(lmin,lmax+0.5,0.5)
cbar.set_ticks(ticks)
ticklbls = (10**(ticks+3)).round().astype(int)
cbar.set_ticklabels(ticklbls)
cbar.ax.yaxis.set_tick_params(labelcolor='w',labelsize=8)

savefig('dynspecLCP.pdf',bbox_inches='tight')



# compare linear and log color scales for inset dyn spec - which shows features better?
#  log plot won.
figure()
a = subplot(121)
a.imshow(log10(rrC.T),aspect='auto',vmin=log10(0.05),vmax=log10(0.8),origin='lower',cmap='seismic')
title('Log Scale')
a = subplot(122)
a.imshow(rrC.T*1000,aspect='auto',vmin=1,vmax=500,origin='lower',cmap='seismic')
title('Linear Scale')
savefig('compare.pdf',bbox_inches='tight')


##### FIGURE: PLOT TIME SERIES #####

print 'making plot 2'

# how many seconds to average over when making time series
n1 = 8
n2 = 8
n3 = 8

i_2GHz = ytick_locs[0]
i_4GHz = ytick_locs[2]
tseries_rr = rebin2d_ma(rr[:,i_2GHz:i_4GHz],(1,500))
tseries_ll = rebin2d_ma(ll[:,i_2GHz:i_4GHz],(1,500))
tseries_V = (tseries_rr - tseries_ll)*1000
t2 = rebin1d(t,n2)
V2 = rebin1d_ma(tseries_V[:,0],n2)
t3 = rebin1d(t,n3)
V3 = rebin1d_ma(tseries_V[:,1],n3)

tseries_rr = rebin2d_ma(rr[:,i_4GHz:],(1,978))
tseries_ll = rebin2d_ma(ll[:,i_4GHz:],(1,978))
tseries_V = (tseries_rr - tseries_ll)*1000
t4 = rebin1d(t,n4)
V4 = rebin1d_ma(tseries_V[:,0],n4)
t6 = rebin1d(t,n6)
V6 = rebin1d_ma(tseries_V[:,1],n6)

figure()
hold(True)
plot(t2/60.,V2,'c',label='2-3 GHz')
plot(t3/60.,V3,'g',label='3-4 GHz')
plot(t4[160:]/60.,V4[160:],'r',label='4-5 GHz')
plot(t6/60.,V6,'k',label='6-8 GHz')
xlabel('Time in minutes since '+tstring+' UTC')
ylabel('Flux (mJy)')
title('UV Ceti Flare Stokes V Lightcurve')
#axis([0,tlen/60.,-5,85])
legend(loc='upper left')
savefig('tseries.pdf',bbox_inches='tight')

'''

