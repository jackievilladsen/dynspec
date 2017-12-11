###########

# plot_rc.py: load a dynamic spectrum, bin to an appropriate resolution, and
#             then make a scatterplot of RCP vs. LCP to identify r_c from the slope

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
          'image.interpolation': 'hanning'}
rcParams.update(params)

savedir='/data/jrv/burst_paper/adleo/'

ds_file_dict = load_burst_filelist()

obs_dict = { 'ADLeo_3L': {'fileroot':'15A-416/ADLeo/3','bands':['L'],'flims':[1.1e9,1.6e9],'nt':600,'nf':64},
             'ADLeo_3P': {'fileroot':'15A-416/ADLeo/3','bands':['P'],'flims':[0.3e9,0.4e9],'nt':600,'nf':64},
             'ADLeo_5L': {'fileroot':'15A-416/ADLeo/5','bands':['L'],'nt':600,'nf':64}}
ds_files = ['/data/jrv/15A-416/ADLeo/3/L/test_clean/ds_ap0_big_RR_n2_ms/tbavg.ms.dynspec']
#            '/data/jrv/15A-416/ADLeo/3/S/ADLeo_3S.tbavg.ms.dynspec']
#ds_files = ['/data/jrv/15A-416/ADLeo/5/L/ADLeo_5L.tbavg.ms.dynspec']
#,'/data/jrv/15A-416/ADLeo/5/S/ADLeo_5S.tbavg.ms.dynspec']

src = 'ADLeo'

my_burst = 'ADLeo_5L'
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
        ds.spec['rr'] = (ds.spec['xx']+ds.spec['yy']-1.0j*(ds.spec['xy']-ds.spec['yx']))/2
        ds.spec['ll'] = (ds.spec['xx']+ds.spec['yy']+1.0j*(ds.spec['xy']-ds.spec['yx']))/2
    if ds_obs is None:
        ds_obs = deepcopy(ds)
    else:
        ds_obs.add_dynspec(ds)
    del ds
ds_obs.mask_RFI(rmsfac=1.5)

nt = burst_dict['nt']
nf = burst_dict['nf']
ds=ds_obs.bin_dynspec(nt=nt,nf=nf,mask_partial=0.5)
#ds2 = ds.mask_SNR(0.5)
ds2 = ds
if 'flims' in burst_dict.keys():
    flims = burst_dict['flims']
    ds2 = ds2.clip(fmin=flims[0],fmax=flims[1])

rr = ds2.spec['rr'].flatten()
ll = ds2.spec['ll'].flatten()
new_mask = rr.mask & ll.mask
rr.mask = new_mask
ll.mask = new_mask

x = real(rr).compressed()
y = real(ll).compressed()

def flin(C,x):
    return C[0]*x
linear = odrpack.Model(flin)
sx = std(imag(rr)) * ones(len(x))
sy = std(imag(ll)) * ones(len(y))
mydata = odrpack.RealData(x,y,sx=sx,sy=sy)
myodr = odrpack.ODR(mydata,linear,beta0=[1.]) #,0.])
myoutput = myodr.run()
C = myoutput.beta
Cerr = myoutput.sd_beta
print 'chi2, I think (\"residual variance\"):', myoutput.res_var

plot(x,y,'.')
xmax = max(ma.max(real(rr)),ma.max(real(ll)))
xmin = min(ma.min(real(rr)),ma.min(real(ll)))
dx = (xmax-xmin)/10
xmax += dx
xmin -= dx
ylin1 = C[0]*xmin #+ C[1]
ylin2 = C[0]*xmax #+ C[1]
plot([xmin,xmax],[ylin1,ylin2],'k-')
axis([xmin,xmax,xmin,xmax])
gca().set_aspect(1)

rc = (1-C[0])/(1+C[0])
print 'rc:', rc
xplt = sort(x)
mA = C[0] + Cerr[0]
mB = C[0] - Cerr[0]
yfit = C[0] * xplt #+ C[1]
yA = mA * xplt #+ C[1]
yB = mB * xplt #+ C[1]
plot(xplt,yA,'g--')
plot(xplt,yB,'r--')
grid(which='major')

rcA = (1-mA)/(1+mA)
rcB = (1-mB)/(1+mB)
print 'rc range:', rcA, rcB

r = mean(x)
l = mean(y)
rc_direct = (r-l)/(r+l)
print 'rc_direct:',rc_direct

show()
