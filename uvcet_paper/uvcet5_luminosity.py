'''
uvcet5_luminosity.py - Calculate time-averaged luminosity in UV Cet's radio aurora based on 2015 UVCet 5.
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
n_MHz = 64  # must be multiple of 2
rmsfac = 10
smax = 0.08
smin = 0.002
scale = 'log'

src = 'UVCet'
epoch='5'

params = {'legend.fontsize': 'small',
          'axes.titlesize': 'medium',
          'axes.labelsize': 'small',
          'xtick.labelsize': 'x-small',
          'ytick.labelsize': 'x-small',
          'image.interpolation': 'hanning'}
rcParams.update(params)

savedir = '/data/jrv/'+src.lower()+'_paper/'

close('all')
figure(figsize=(6,8))

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
dsVLBA_bin = dsVLBA.bin_dynspec(nt=nt_VLBA,nf=nf_VLBA,mask_partial=0.75)
ds_bin = ds.bin_dynspec(nt_VLA,nf,mask_partial=0.75) # will flag if 3/4 of contributing pixels flagged --> 50% sensitivity
dsP_bin = dsP.bin_dynspec(nt_P,nf,mask_partial=0.75)
'''
# calculate quiescent spectrum
t_quiet_0 = 0
t_quiet_f = 40 #min

ds_quiet = ds_bin.clip(tmin=t_quiet_0,tmax=t_quiet_f)
dsVLBA_quiet = dsVLBA_bin.clip(tmin=t_quiet_0,tmax=t_quiet_f)
dsP_quiet = dsP_bin.clip(tmin=t_quiet_0,tmax=t_quiet_f)

spec_quiet = ds_quiet.bin_dynspec(nf=8,nt=len(ds_quiet.get_tlist()))
specVLBA_quiet = dsVLBA_quiet.bin_dynspec(nf=1,nt=len(dsVLBA_quiet.get_tlist()))
specP_quiet = dsP_quiet.bin_dynspec(nf=1,nt=len(dsP_quiet.get_tlist()))
'''
# mask dynspec below certain flux
dsVLBA_bin=dsVLBA_bin.mask_SNR(rmsfac=6.)
ds_bin = ds_bin.mask_SNR(rmsfac=6.)
dsP_bin = dsP_bin.mask_SNR(rmsfac=1.4)
dsVLBA_bin.spec['v'].mask = dsVLBA_bin.spec['i'].mask
ds_bin.spec['v'].mask = ds_bin.spec['i'].mask
dsP_bin.spec['v'].mask = dsP_bin.spec['i'].mask

specVLA = ds_bin.bin_dynspec(nf=1,nt=len(ds_bin.get_tlist()),mask_partial=0.8)
specVLBA = dsVLBA_bin.bin_dynspec(nf=1,nt=len(dsVLBA_bin.get_tlist()),mask_partial=1)
specP = dsP_bin.bin_dynspec(nf=1,nt=len(dsP_bin.get_tlist()),mask_partial=1)
specI = ma.concatenate([specP.spec['i'],specVLA.spec['i'],specVLBA.spec['i']],1)[0,:]*1e3
specV = ma.concatenate([specP.spec['v'],specVLA.spec['v'],specVLBA.spec['v']],1)[0,:]*1e3
mask = specI.mask
rc = real(specV)/real(specI)

dtVLA = sum(1-ds_bin.spec['i'].mask,0) * ds_bin.dt() / 60.
dtVLBA = sum(1-dsVLBA_bin.spec['i'].mask,0) * dsVLBA_bin.dt() / 60.
dtP = sum(1-dsP_bin.spec['i'].mask,0) * dsP_bin.dt() / 60.
P = 5.45 * 60
dt = concatenate([dtP,dtVLA,dtVLBA])
dt = ma.masked_array(dt,mask=mask)

fVLA = ds_bin.f/1.e9
fVLBA = dsVLBA_bin.f/1.e9
fP = dsP_bin.f/1.e9
f = ma.concatenate([fP,fVLA,fVLBA])

E = specI * dt * 60
# currently E is in mJy * sec
# 1 mJy * sec = 1e-26 erg/s/Hz/cm^2 * s = 1e-26 erg/Hz/cm^2
# convert to erg/Hz/ster
# 1e-26 erg/Hz/cm^2 * 4 * pi (2.676 pc)^2 / (4*pi) = 6.82e11 erg/Hz/ster
E *= 6.82e11

P_rot = 5.45 * 60 # rotation period in minutes
theta = dt/P_rot * 2 * pi
theta_deg = theta * 180./pi
deg_factor = 180. * pi / P_rot
ster = pi * (theta/2)**2
E_ster = E * ster

E_iso = real(sum(E))*n_MHz*1.e6 *4 * pi
L_iso = E_iso / (P_rot * 60)
E_beamed = real(sum(E_ster))*n_MHz*1e6
L_beamed = E_beamed / (P_rot * 60)
print 'Energy radiated in burst (assuming isotropic emission):',E_iso,'erg'
print 'Luminosity (assuming isotropic):',L_iso,'erg/s =', L_iso/1.e7, 'W'
print 'Energy radiated in burst (assuming beamed):',E_beamed, 'erg'
print 'Luminosity (assuming beamed):', L_beamed, 'erg/s =', L_beamed/1.e7, 'W'

df = n_MHz/1.e3
f_even = arange(f[0],f[-1],df)
ind = ~E.mask
func_iso = interp1d(f[ind],E[ind])
func_beamed = interp1d(f[ind],E_ster[ind])
Esmooth = func_iso(f_even)
E_ster_smooth = func_beamed(f_even)

E_iso = real(sum(Esmooth))*n_MHz*1.e6 *4 * pi
L_iso = E_iso / (P_rot * 60)
E_beamed = real(sum(E_ster_smooth))*n_MHz*1e6
L_beamed = E_beamed / (P_rot * 60)
print 'Energy radiated in burst (interpolated, assuming isotropic emission):',E_iso,'erg'
print 'Luminosity (interpolated, assuming isotropic):',L_iso,'erg/s =', L_iso/1.e7, 'W'
print 'Energy radiated in burst (interpolated, assuming beamed):',E_beamed, 'erg'
print 'Luminosity (interpolated, assuming beamed):', L_beamed, 'erg/s =', L_beamed/1.e7, 'W'


subplot(411)
semilogx(f,specI,'b.')
semilogx(f,specV,'g.')
ylabel('Avg burst flux (mJy)')
legend(('Stokes I','Stokes V'))
axis([0.2,9,0,70])
subplot(412)
semilogx(f,rc,'k.')
ylabel('Circ pol fraction')
axis([0.2,9,0,1])
subplot(413)
semilogx(f,dt,'k.')
ylabel('Duration (min)')
#gca().tick_params(labeltop=False, labelright=True)
#ylabel('Degrees')
axis([0.2,9,0,100])
subplot(414)
semilogx(f,E*4*pi,'b.')
semilogx(f,E_ster,'g.')
semilogx(f_even,Esmooth*4*pi,'b')
semilogx(f_even,E_ster_smooth,'g')
ylabel('Burst energy (erg/Hz)')
xlabel('Frequency (GHz)')
legend(('Isotropic','Beamed'))
axis([0.2,9,0,ma.max(real(E*4*pi))*1.1])
savefig(savedir+'UVCet5_burstspec.pdf',bbox_inches='tight')


# plot masked dynspec so I can see which data are being used
for pol in ['i','v']:
    pp = {'pol':pol,'smin':smin,'smax':smax,'trim_mask':False,'ar0':0.8,'scale':scale,'axis_labels':['ylabel','cbar','cbar_label'],'dy':0.5}
    ppP = {'pol':pol,'smin':smin,'smax':smax,'ar0':0.2,'trim_mask':False,'scale':scale,'axis_labels':['xlabel','cbar'],'dy':0.1}
    ppX = {'pol':pol,'smin':smin,'smax':smax,'ar0':0.2,'trim_mask':False,'scale':scale,'axis_labels':['cbar'],'dy':0.1}

    clf()
    ax=axes([0,0.18,1,0.81])
    ds_bin.plot_dynspec(plot_params=pp)
    ax.xaxis.set_visible(False)
    ax = axes([0,0.97,1,0.2])
    dsVLBA_bin.plot_dynspec(plot_params=ppX)
    cb = gca().images[-1].colorbar
    cb.remove()
    ax.xaxis.set_visible(False)
    title(src + ' ' + epoch)
    ax.xaxis.set_visible(False)
    ax = axes([0,0,1,0.2])
    dsP_bin.plot_dynspec(plot_params=ppP)
    cb = gca().images[-1].colorbar
    cb.remove()
    figfile = savedir + 'mask'+pol.upper()+'.png'
    savefig(figfile,bbox_inches='tight')
