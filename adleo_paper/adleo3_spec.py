'''
adleo3_spec.py - Calculate time-averaged spectrum of AD Leo epoch 3
'''

#from dynspec import load_dict
#from dynspec.plot import Dynspec
from pylab import *
import pickle
from scipy.interpolate import interp1d
#import os

n_sec_P = 600 # must be multiple of 6
n_sec_VLA = 90
n_sec_VLBA = 150
n_MHz = 16  # must be multiple of 2
rmsfac = 10
smax = 0.08
smin = 0.002
scale = 'log'

src = 'ADLeo'
epoch='4'

params = {'legend.fontsize': 'small',
          'axes.titlesize': 'medium',
          'axes.labelsize': 'small',
          'xtick.labelsize': 'x-small',
          'ytick.labelsize': 'x-small',
          'image.interpolation': 'hanning'}
rcParams.update(params)

savedir = '/data/jrv/'+src.lower()+'_paper/'

close('all')
figure(figsize=(6,4))

# load multi-band dynspec
savefile = savedir + src + '_' + epoch + '.dynspec.pickle'
ds_dict = pickle.load( open( savefile, "rb" ) )

ds = ds_dict['VLA']
dsP = ds_dict['P']
dsVLBA = ds_dict['VLBA']

dsP.mask_RFI(rmsfac=1.5)

# Calculate number of dynspec pixels to bin together
nt_P = int(round(n_sec_P/dsP.dt()))  # number of integrations to bin together (current resolution is 6 sec)
nt_VLA = int(round(n_sec_VLA/ds.dt()))
nt_VLBA = int(round(n_sec_VLBA/dsVLBA.dt()))
nf = n_MHz/int(round(ds.df()/1e6))  # number of channels to bin together (current resolution is 2 MHz)
nf_VLBA = n_MHz/int(round(dsVLBA.df()/1e6))
    
# bin dynspec to improve signal-to-noise ratio
dsVLBA_bin = dsVLBA.bin_dynspec(nt=nt_VLBA,nf=nf_VLBA,mask_partial=0.5)
ds_bin = ds.bin_dynspec(nt_VLA,nf,mask_partial=0.5) # will flag if 3/4 of contributing pixels flagged --> 50% sensitivity
dsP_bin = dsP.bin_dynspec(nt_P,nf,mask_partial=0.75)
specVLA = ds_bin.bin_dynspec(nf=1,nt=len(ds_bin.get_tlist()),mask_partial=0.8)
specVLBA = dsVLBA_bin.bin_dynspec(nf=1,nt=len(dsVLBA_bin.get_tlist()),mask_partial=1)
specP = dsP_bin.bin_dynspec(nf=1,nt=len(dsP_bin.get_tlist()),mask_partial=1)
specI = ma.concatenate([specP.spec['i'],specVLA.spec['i'],specVLBA.spec['i']],1)[0,:]*1e3
specV = ma.concatenate([specP.spec['v'],specVLA.spec['v'],specVLBA.spec['v']],1)[0,:]*1e3
mask = specI.mask
rc = real(specV)/real(specI)

fVLA = ds_bin.f/1.e9
fVLBA = dsVLBA_bin.f/1.e9
fP = dsP_bin.f/1.e9
f = ma.concatenate([fP,fVLA,fVLBA])

subplot(211)
semilogx(f,specI,'b.')
semilogx(f,specV,'g.')
gca().axhline(0,color='k')
ylabel('Time-averaged flux (mJy)')
legend(('Stokes I','Stokes V'))
axis([0.2,9,-18,18])
subplot(212)
semilogx(f,rc,'k.')
gca().axhline(0,color='k')
gca().axhline(1,color='k',ls='--')
gca().axhline(-1,color='k',ls='--')
ylabel('Circ pol fraction')
axis([0.2,9,-1.4,1.4])
xlabel('Frequency (GHz)')
savefig(savedir+'ADLeo'+epoch+'_burstspec.pdf',bbox_inches='tight')

