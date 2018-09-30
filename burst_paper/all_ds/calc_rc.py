###########

# calc_rc.py: load a dynamic spectrum, bin to an appropriate resolution, and clip in time (as needed) and frequency to just the burst,
#             then calculate rc = mean(V)/mean(I)

# for P band data, I took meanV, sig_meanV, meanI, sig_meanI from here, then ran rcirc.py to get a Bayesian confidence interval or lower limit on rc,
#        because the error bars on I are so large that you get problematic things like |V| > I

###########

import dynspec.plot
import dynspec.pipeline_utils
reload(dynspec.plot)
reload(dynspec.pipeline_utils)

from dynspec.pipeline_utils import load_burst_filelist
from dynspec import load_dict
from dynspec.plot import *
from pylab import *
import os, subprocess
import matplotlib.gridspec as gridspec
import scipy.odr.odrpack as odrpack

params = {'legend.fontsize': 'small',
          'axes.titlesize': 'small',
          'axes.labelsize': 'x-small',
          'xtick.labelsize': 'xx-small',
          'ytick.labelsize': 'xx-small',
          'image.interpolation': 'none'}
mpl.rcParams.update(params)

def get_srcname(filename):
    # get source name from dynspec filename
    return filename.split('/')[4]

def get_dist(src):
    # return distance in pc to src
    # (will use this to convert flux to flux at 1 pc)
    dstar = {'ADLeo':4.7,
             'UVCet':2.7,
             'EQPeg':6.2,
             'EVLac':5.1,
             'YZCMi':6.0}
    return dstar[src]

savedir = '/data/jrv/burst_paper/rc/'

ds_file_dict = load_burst_filelist()

obs_dict = { 'ADLeo_3L': {'fileroot':'15A-416/ADLeo/3','bands':['L'],'flims':[1.1e9,1.55e9]},
             'ADLeo_3P': {'fileroot':'15A-416/ADLeo/3','bands':['P'],'flims':[0.28e9,0.4e9]},
             'ADLeo_4L_long': {'fileroot':'15A-416/ADLeo/4','bands':['L'],'flims':[1e9,1.6e9]},
             'ADLeo_4L_short': {'fileroot':'15A-416/ADLeo/4','bands':['L','S'],'flims':[1.6e9,2.2e9],'tlims':[136.,137.],'nt':5}, 
             'ADLeo_4S_short': {'fileroot':'15A-416/ADLeo/4','bands':['S'],'flims':[2.8e9,4e9],'tlims':[133.,145.]},
             'ADLeo_5LS': {'fileroot':'15A-416/ADLeo/5','bands':['L','S'],'tlims':[120.,1000.],'flims':[1.e9,2.5e9]},
             'ADLeo_5P': {'fileroot':'15A-416/ADLeo/5','bands':['P'],'flims':[0.29e9,0.36e9],'tlims':[184.,190.],'nt':5.},
             'EQPeg_2': {'fileroot':'15A-416/EQPeg/2','bands':['L','S'],'flims':[1.e9,3.e9],'tlims':[0.,30.]},
             'EQPeg_2P':{'fileroot':'15A-416/EQPeg/2','bands':['P'],'flims':[0.35e9,0.385e9],'tlims':[0.,35.],'mask_partial':1.,'rmsfac':5.0,'nf':16,'nt':30}, # this is the one used for Table 3 in the paper, others are for exploring
             'EQPeg_2P_v2':{'fileroot':'15A-416/EQPeg/2','bands':['P'],'flims':[0.2e9,0.45e9],'tlims':[0.,60.],'mask_partial':1.,'rmsfac':5.0,'nf':16},
             'EQPeg_2P_shortB':{'fileroot':'15A-416/EQPeg/2','bands':['P'],'flims':[0.29e9,0.305e9],'tlims':[12.,25.],'mask_partial':1.,'rmsfac':5.0},
             'EQPeg_2P_shortC':{'fileroot':'15A-416/EQPeg/2','bands':['P'],'flims':[0.32e9,0.34e9],'tlims':[0.,25.],'mask_partial':1.,'rmsfac':5.0,'nf':8,'nt':20},
             'YZCMi_1_long': {'fileroot':'15A-416/YZCMi/1','bands':['S'],'flims':[3.e9,3.7e9]},
             'YZCMi_1_shortA': {'fileroot':'15A-416/YZCMi/1','bands':['L'],'flims':[1.e9,1.4e9],'tlims':[62.,65.]},
             'YZCMi_1_shortB': {'fileroot':'15A-416/YZCMi/1','bands':['L'],'flims':[1.6e9,1.8e9],'tlims':[61.,65.]},
             'YZCMi_1_shortC': {'fileroot':'15A-416/YZCMi/1','bands':['S'],'flims':[2e9,2.8e9],'tlims':[61.,63.5]},
             'YZCMi_1_short': {'fileroot':'15A-416/YZCMi/1','bands':['L','S'],'flims':[1e9,3.7e9],'tlims':[60.,70.]},
             'YZCMi_2': {'fileroot':'15A-416/YZCMi/2','bands':['S'],'flims':[2e9,2.7e9],'tlims':[0.,45.]},
             '2013_UVCet_1_long': {'fileroot':'13A-423/UVCet/1','bands':['L','S'],'tlims':[78.,1000.]},
             '2013_UVCet_1_short': {'fileroot':'13A-423/UVCet/1','bands':['L','S'],'tlims':[33.,38.],'flims':[1.e9,6.e9]},
             '2013_UVCet_2': {'fileroot':'13A-423/UVCet/2','bands':['L','S'],'tlims':[30.,90.],'flims':[1.e9,6.e9]},
             '2013_UVCet_2A': {'fileroot':'13A-423/UVCet/2','bands':['L','S'],'tlims':[40.,60.],'flims':[1.e9,2.2e9]},
             '2013_UVCet_2B': {'fileroot':'13A-423/UVCet/2','bands':['S'],'tlims':[78.,79.1],'flims':[3.2e9,4.e9]},
             '2015_UVCet_1': {'fileroot':'15A-416/UVCet/1','bands':['L'],'tlims':[48.,56.],'flims':[1.2e9,2.e9]},
             '2015_UVCet_2LS': {'fileroot':'15A-416/UVCet/2','bands':['L','S'],'tlims':[73.,1000.],'flims':[1.e9,3.7e9]},
             '2015_UVCet_2P': {'fileroot':'15A-416/UVCet/2','bands':['P'],'tlims':[73.,1000.],'nt':20,'flims':[0.3e9,0.5e9]},
             '2015_UVCet_3LS': {'fileroot':'15A-416/UVCet/3','bands':['L','S'],'tlims':[30.,110.]},

             '2015_UVCet_3P': {'fileroot':'15A-416/UVCet/3','bands':['P'],'tlims':[45.,110.],'flims':[0.29e9,0.5e9],'nt':20},
             '2015_UVCet_4LS_long': {'fileroot':'15A-416/UVCet/4','bands':['L','S'],'tlims':[184.,1000.]},
             '2015_UVCet_4LS_shortA': {'fileroot':'15A-416/UVCet/4','bands':['L','S'],'tlims':[68.5,71.]},
             '2015_UVCet_4LS_shortB': {'fileroot':'15A-416/UVCet/4','bands':['L','S'],'tlims':[89.,96.],'flims':[2.4e9,4.e9]},
             '2015_UVCet_4P_long': {'fileroot':'15A-416/UVCet/4','bands':['P'],'tlims':[194.,1000.],'flims':[0.35e9,0.5e9],'nt':20},
             '2015_UVCet_4P_short': {'fileroot':'15A-416/UVCet/4','bands':['P'],'tlims':[130.,140.],'flims':[0.28e9,0.37e9],'nt':10},
             '2015_UVCet_5LS': {'fileroot':'15A-416/UVCet/5','bands':['L','S'],'tlims':[68.,150.],'mask_partial':0.75},
             '2015_UVCet_5P': {'fileroot':'15A-416/UVCet/5','bands':['P'],'tlims':[68.,128.],'flims':[0.25e9,0.5e9],'nt':40,'nf':64}}

my_burst = 'YZCMi_2'

burst_dict = obs_dict[my_burst]

print '--------------------------------\n', my_burst

band_list = burst_dict['bands']
filepath = '/data/jrv/'+burst_dict['fileroot']
filelist_all_bands = ds_file_dict[filepath]
ds_files = [item for item in filelist_all_bands if item.split('/')[6] in band_list]

ds_obs = None
for f in ds_files:
    band = f.split('/')[6]
    params={'filename':f,'uniform':True}
    ds = Dynspec(params)
    if band=='P':
        ds.spec['i'] = (ds.spec['xx']+ds.spec['yy'])/2
        ds.spec['v'] = (ds.spec['xy']-ds.spec['yx'])/(2.0j)
        ds.spec['u'] = (ds.spec['xy']+ds.spec['yx'])/2
        #del ds.spec['xx']
        #del ds.spec['xy']
        #del ds.spec['yx']
        #del ds.spec['yy']
    else:
        ds.spec['i'] = (ds.spec['rr']+ds.spec['ll'])/2
        ds.spec['v'] = (ds.spec['rr']-ds.spec['ll'])/2
        del ds.spec['rr']
        del ds.spec['ll']
    if ds_obs is None:
        ds_obs = deepcopy(ds)
    else:
        ds_obs.add_dynspec(ds)
    del ds

ds2 = ds_obs
if 'tlims'in burst_dict.keys():
    tlims = burst_dict['tlims']
    ds2 = ds2.clip(tmin=tlims[0],tmax=tlims[1],trim_mask=False)
    print 'clipping ds by time:', tlims
if 'flims' in burst_dict.keys():
    flims = burst_dict['flims']
    ds2 = ds2.clip(fmin=flims[0],fmax=flims[1],trim_mask=False)
    print 'clipping ds by freq', flims

rmsfac = burst_dict.get('rmsfac',3.0)
ds2.mask_RFI(rmsfac=rmsfac)

nt0 = len(ds2.time)/30
nf0 = len(ds2.f)/30
nt = burst_dict.get('nt',nt0)
nf = burst_dict.get('nf',nf0)
mask_partial=burst_dict.get('mask_partial',0.5)
ds2=ds2.bin_dynspec(nt=nt,nf=nf,mask_partial=mask_partial) # need to bin in time and freq to get accurate estimate of sig_rc due to correlated noise

if band != 'P':
    dsI = ds2.spec['i']
    dsV = ds2.spec['v']
    
    new_mask = ma.mask_or(dsI.mask,dsV.mask)
    dsI.mask = new_mask
    dsV.mask = new_mask
    
    sI = dsI.compressed()
    sV = dsV.compressed()
    
    meanI = ma.mean(real(sI))
    meanV = ma.mean(real(sV))
    sigI = ma.std(imag(sI)) / sqrt(len(sI)-1)
    sigV = ma.std(imag(sV)) / sqrt(len(sV)-1)
    rc = meanV/meanI
    if meanI < 0:
        print 'Warning! meanI is less than 0:', meanI
    sig_rc = abs(rc) * sqrt( (sigI/meanI)**2 + (sigV/meanV)**2 )
    print my_burst, 'r_c:', round(rc,2), '+/-', round(sig_rc,3)
    print 'mean flux in V:', meanV*1e3, '+-', sigV*1e3, 'mJy --> SNR V:', abs(meanV)/sigV
    print 'mean flux in I:', meanI*1e3, '+-', sigI*1e3, 'mJy --> SNR I:', meanI/sigI

else:
    dsI = ds2.spec['i']
    dsV = ds2.spec['v']
    dsU = ds2.spec['u']
    
    new_mask = ma.mask_or(dsI.mask,dsV.mask) #,dsU.mask)
    dsI.mask = new_mask
    dsV.mask = new_mask
    dsU.mask = new_mask
    
    sI = dsI.compressed()
    sV = dsV.compressed()
    sU = dsU.compressed()
    
    meanI = ma.mean(real(sI))
    print 'meanI:', meanI
    meanV = ma.mean(real(sV))
    print 'meanV:', meanV
    meanU = ma.mean(real(sU))
    print 'meanU:', meanU

    sigI = ma.std(imag(sI)) / sqrt(len(sI)-1)
    sigV = ma.std(imag(sV)) / sqrt(len(sV)-1)
    sigU = ma.std(imag(sU)) / sqrt(len(sU)-1)
    meanUV = sqrt(meanV**2+meanU**2)
    sigUV = sqrt( sigU**2 * meanU**2 + sigV**2 * meanV**2 ) / meanUV
    rc = meanUV/meanI
    if meanI < 0:
        print 'Warning! meanI is less than 0:', meanI
    sig_rc = abs(rc) * sqrt( (sigI/meanI)**2 + (sigV/meanV)**2 + (sigU/meanU)**2 )
    print my_burst, 'r_c:', round(rc,2), '+/-', round(sig_rc,3)
    print 'mean flux in UV:', meanUV*1e3, '+-', sigUV*1e3, 'mJy --> SNR UV:', meanUV/sigUV
    print 'mean flux in I:', meanI*1e3, '+-', sigI*1e3, 'mJy --> SNR I:', meanI/sigI

# calculate total energy in the burst:
BW = len(ds2.f) * ds2.df() # in Hz
duration = len(ds2.time) * ds2.dt() # in seconds
srcname = get_srcname(filepath)
d_star = get_dist(srcname)
burst_energy = meanI / d_star**2 * (1.2e15) * BW * duration   # in erg
sig_energy = sigI / d_star**2 * (1.2e15) * BW * duration   # in erg
print 'Burst energy:', burst_energy/(1.e22), '+-', sig_energy/(1.e22), 'x 1e22 erg (for',srcname,'at',d_star,'pc)'




'''
#if band=='P':
dsP1squared = real ( ds2.spec['xy'] * ds2.spec['yx'].conj() )
dsP4squared = ds2.spec['xy'] * ds2.spec['yx'].conj()
dsP2 = (abs(ds2.spec['u'])**2 + abs(ds2.spec['v'])**2)**0.5
dsP3 = ((abs(ds2.spec['xy'])**2 + abs(ds2.spec['yx'])**2)/2)**0.5   # this is the same as P2

dsXY = ds2.spec['xy']
dsYX = ds2.spec['yx']
dsU = ds2.spec['u']
new_mask = ma.mask_or(dsI.mask,dsV.mask,dsP1.mask,dsP2.mask)
new_mask = ma.mask_or(new_mask,dsU.mask,dsXY.mask,dsYX.mask)
dsI.mask = new_mask
dsU.mask = new_mask
dsV.mask = new_mask
dsP1.mask = new_mask
dsP2.mask = new_mask
dsXY.mask = new_mask
dsYX.mask = new_mask
sI = dsI.compressed()
sV = dsV.compressed()
sP1squared = dsP1squared.compressed()
sP2 = dsP2.compressed()
sP3 = dsP3.compressed()
sXY = dsXY.compressed()
sYX = dsYX.compressed()
sU = dsU.compressed()
meanI = ma.mean(real(sI))
meanV = ma.mean(real(sV))
meanP1squared = ma.mean(real(sP1squared))
meanP2 = mean(sP2)
meanP3 = mean(sP3)
sigI = ma.std(imag(sI)) / sqrt(len(sI)-1)
sigV = ma.std(imag(sV)) / sqrt(len(sV)-1)
'''


# attempt to measure peak flux
signV = rc/abs(rc)
absV = real(sV) * signV
maxV = max(absV) * signV
peakV95 = percentile(absV,95) * signV
peakV98 = percentile(absV,98) * signV
peakV99 = percentile(absV,99) * signV
#print 'max V:', maxV
print '98% V:',int(round(peakV98*1e3,0)),'mJy'

maxI = max(real(sI))
peakI95 = percentile(real(sI),95)
peakI98 = percentile(real(sI),98)
peakI99 = percentile(real(sI),99)
#print 'max I:',maxI
print '98% I:',int(round(peakI98*1e3,0)),'mJy'

print 'V98/min(abs(rc),1):', int(round(peakV98/min(abs(rc),1) * signV * 1e3,0)), 'mJy\n'

if band != 'P':
    close('all')
    figure(figsize=(8,5))
    subplot(121)
    ds2.plot_dynspec({'pol':'i','smin':-peakI98,'smax':peakI98})
    title('Stokes I')
    subplot(122)
    ds2.plot_dynspec({'pol':'v','smin':-peakI98,'smax':peakI98})
    title('Stokes V')
    savefig(savedir+my_burst+'_rc_calc_ds.pdf')
else:
    close('all')
    figure(figsize=(10,10))
    subplot(221)
    ds2.plot_dynspec({'pol':'i','smin':-peakI98,'smax':peakI98})
    title('Stokes I')
    subplot(222)
    ds2.plot_dynspec({'pol':'u','smin':-peakI98,'smax':peakI98})
    title('Stokes U')
    subplot(223)
    ds2.plot_dynspec({'pol':'v','smin':-peakI98,'smax':peakI98})
    title('Stokes V')
    subplot(224)
    ds2.spec['uv'] = sqrt(abs(ds2.spec['u'])**2 + abs(ds2.spec['v'])**2)
    ds2.plot_dynspec({'pol':'uv','smin':-peakI98,'smax':peakI98})
    title('sqrt(U^2+V^2)')
    savefig(savedir+my_burst+'_rc_calc_ds.pdf')

show()
