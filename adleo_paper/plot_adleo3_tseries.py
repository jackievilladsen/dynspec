'''
plot_adleo_ds.py - Load ADLeo multi-band dynamic spectrum for a given epoch, bin to specified resolution, and plot to file
'''

#from dynspec.plot import Dynspec
from pylab import *
import pickle
#import os

n_sec_VLBA = 120 # must be multiple of 6
n_sec_VLA = 12 # must be multiple of 6
n_MHz = 16  # must be multiple of 2
mask_partial=0.75
weight_mode='rms'

src = 'ADLeo'
epoch='3'

params = {'legend.fontsize': 'small',
          'axes.titlesize': 'medium',
          'axes.labelsize': 'small',
          'xtick.labelsize': 'x-small',
          'ytick.labelsize': 'x-small',
          'image.interpolation': 'hanning'}
rcParams.update(params)

savedir = '/data/jrv/adleo_paper/'

ntVLBA = n_sec_VLBA/6  # number of integrations to bin together (current resolution is 6 sec)
ntVLA = n_sec_VLA/6  # number of integrations to bin together (current resolution is 6 sec)
nf = n_MHz/2  # number of channels to bin together (current resolution is 2 MHz)

close('all')
figure(figsize=(6,6))

# load multi-band dynspec
try:
    ds_dict
except:
    savefile = savedir + src + '_' + epoch + '.dynspec.pickle'
    ds_dict = pickle.load( open( savefile, "rb" ) )

ds = ds_dict['VLA']
dsVLBA = ds_dict['VLBA']

# Calculate VLBA time series
dsVLBA = dsVLBA.bin_dynspec(nt=nt,nf=1,mask_partial=mask_partial)
dsVLBA = dsVLBA.tseries(weight_mode='flat')
tVLBA = dsVLBA.get_tlist() * dsVLBA.dt()/60.
iflux = dsVLBA.spec['i']*1e3
vflux = dsVLBA.spec['v']*1e3
flux_err = std(imag(iflux))

# Calculate lower frequency tseries
ds_bin = ds.bin_dynspec(nt=nt,nf=1,mask_partial=mask_partial)

ds34 = ds_bin.tseries(fmin=3.e9,fmax=4.e9,weight_mode=weight_mode)
t34 = ds34.get_tlist() * ds34.dt()/60.
i34 = ds34.spec['i']*1e3
v34 = ds34.spec['v']*1e3
i_err34 = std(imag(i34))
v_err34 = std(imag(v34))

ds23 = ds_bin.tseries(fmin=2.4e9,fmax=3.e9,weight_mode=weight_mode)
t23 = ds23.get_tlist() * ds23.dt()/60.
i23 = ds23.spec['i']*1e3
v23 = ds23.spec['v']*1e3
i_err23 = std(imag(i23))
v_err23 = std(imag(v23))

ds12 = ds_bin.tseries(fmin=1.e9,fmax=1.5e9,weight_mode=weight_mode)
t12 = ds12.get_tlist() * ds12.dt()/60.
i12 = ds12.spec['i']*1e3
v12 = ds12.spec['v']*1e3
i_err12 = std(imag(i12))
v_err12 = std(imag(v12))

clf()
#subplot(221)
plot(tVLBA,iflux)
plot(t34,-v34)
plot(t23,-v23)
#plot(t12,i12)
plot(t12,-v12)
xlabel('Time in minutes')
ylabel('Flux (mJy)')
legend(('8.4 GHz I','3-4 GHz -V','2.4-3 GHz -V','1-1.5 GHz -V'))
gca().axhline(0,color='k')
savefig(savedir+'ADLeo3_tseries.pdf',bbox_inches='tight')
close('all')
'''
subplot(223)
plot(tVLBA,imag(iflux))
plot(t34,imag(i34))
plot(t23,imag(i23))
plot(t12,imag(i12))

subplot(222)
plot(tVLBA,vflux)
plot(t34,v34)
plot(t23,v23)
plot(t12,v12)
subplot(224)
plot(tVLBA,imag(vflux))
plot(t34,imag(v34))
plot(t23,imag(v23))
plot(t12,imag(v12))

'''
'''
    # plot tseries above dynspec (time axis aligned) and save figure
    # to do: set x-axis limits for time series
    clf()
    ax=axes([0,0,1,1])
    ds_bin.plot_dynspec(plot_params=pp)
    ax=axes([0,0.96,0.915,0.2]) #2: 0.96
    axhline(0,color='k')
    errorbar(t,iflux,flux_err)
    plot(t,vflux)
    xlim([min(t),max(t)])
    ylabel('VLBA Flux (mJy)')
    title(src + ' ' + epoch)
    gca().xaxis.set_visible(False)
    figfile = savedir + src + '_' + epoch + '_dynspec.pdf'
    savefig(figfile,bbox_inches='tight')


'''
