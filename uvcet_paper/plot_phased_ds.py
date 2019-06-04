'''
plot_phased_ds.py - Plot 1-4 GHz Stokes I dynamic spectrum for each epoch as
                    a function of phase for a user-defined period.
    (I got the period by doing a Plavchan periodogram on 1-2 GHz RCP data.)

    Run this code in ipython.680
'''

# Issue: loading the dictionaries requires a different, older version of astropy (default in ipython)
#  but this older version of astropy does not have coordinates.EarthLocation.of_site()
#  (whereas the version that is default to use in CASA does have it)

#import sys
#sys.path.append("/data/jrv/casa_utils/astropy-1.3.2/")  # path to newer version of astropy

import dynspec.plot
reload(dynspec.plot)

from dynspec import load_dict
from dynspec.plot import *
from pylab import *
import os, subprocess
import matplotlib.gridspec as gridspec
from astropy import time, coordinates, units

def get_obsfile(obsname):
    # take an obs name such as '15A-416_YZCMi_1' and return srcname ('YZCMi')
    # and file directory ('/data/jrv/15A-416/YZCMi/1')
    names = obsname.split('_')
    srcname = names[1]
    obsfile = '/data/jrv/'+names[0]+'/'+names[1]+'/'+names[2]
    return obsfile, srcname

def mjd_to_phase(t_mjd,period_target='UVCet',t0_mjd = 0.0):
    '''
    mjd_to_phase_UVCet(t_mjd,period_target='UVCet',t0_mjd = 57207.5127083-0.12)
    - default t0_mjd is start time of 2015_UVCet_3 VLA dynspec minus roughly half a period
    
    Phase is scaled from 0 to 1 (where 1 is 360 degrees).
    
    This routine keeps the observation in time order - so if the observation starts at phase of 0.9,
    it may end at phase of 1.5.
    '''
    period_dict = {'UVCet':0.2268,'UVCetMod':5.44709/24,'BLCet':0.2430,'BLCetMod':0.2430-0.00052}
    try:
        Prot = period_dict[period_target]
    except KeyError:
        print 'period_target is not in period_dict - using period of 1 day'
        Prot = 1.0
    t_rel = t_mjd - t0_mjd
    n_periods_off = floor(t_rel[0] / Prot)
    phase = t_rel / Prot - n_periods_off
    return phase

def convert_ds_time_to_phase(ds,period_target='UVCet',t0_TDB = 57207.5127083-0.2,epoch='2'):
    t_mjd = ds.time.mjd
    
    # convert mjd time list from UTC to TDB (subbing off light travel time to get barycentric time)
    # not using the astropy functions b/c for some reason my data files
    #   will only load with an older version of astropy that doesn't have them,
    #   so I got these values (for light_travel_time diff at start of obs) from
    #   running make_uvcet_tseries.py in CASA with newer version of astropy
    light_time_dict = {'1':-0.00407680646506,'2':-0.00402307342939, '3':-0.000322840436988,'4':0.000897346506331,'5':0.00442142878382}
    t_TDB = ds.time.tdb.mjd + light_time_dict[epoch]
    
    # subtract off t0 and convert to phase
    phase = mjd_to_phase(t_TDB,period_target,t0_TDB)
    ds.phase = phase
    print 'Storing rotational phase in ds.phase using period of',period_target+'. WARNING:',  \
          'ds.phase does not update if you bin or truncate the time list, so make sure to do',\
          'this last before plotting.'

params = {'legend.fontsize': 'small',
          'axes.titlesize': 'small',
          'axes.labelsize': 'small',
          'xtick.labelsize': 'x-small',
          'ytick.labelsize': 'x-small',
          'image.interpolation': 'nearest'}
rcParams.update(params)

loadfile = '/data/jrv/burst_paper/all_burst_epoch_dynspec_LSband.npy'
ds_list = load_dict(loadfile)
loadfileP = '/data/jrv/burst_paper/all_burst_epoch_dynspec_Pband.npy'
dsP_list = load_dict(loadfileP)

close('all')

# params for ds plots
fig_params = {
    'tint_P':  300,
    'tint_LS': 60,
    'df_MHz_P':  16,
    'df_MHz_LS': 16,
    'smax_P':  None,
    'smax_LS': 0.1,
    'pixflag_sigfacP':  7.,
    'pixflag_sigfacLS': 10.,
    'chanflag_sigfacP':  3.,
    'chanflag_sigfacLS': 7.,
    'colorscale_P':'linear',
    'colorscale_LS':'linear',
    'maskpartial_P':0.5,
    'maskpartial_LS':0.5,
    'linthresh_P':None,
    'linthresh_LS':None}

### PLOT INDIVIDUAL OBSERVATIONS ###

#obs_list = ['15A-416_UVCet_1','15A-416_UVCet_2','15A-416_UVCet_3','15A-416_UVCet_4','15A-416_UVCet_5']
obs_list = ['15A-416_UVCet_2','15A-416_UVCet_3','15A-416_UVCet_4','15A-416_UVCet_5']

fig_max_width=6.5
fig_max_height=8.25

# set flux limits for LS band
smax = fig_params['smax_LS']
smin = -smax  # make colorbar symmetric about zero

# set axis ratio to 'auto' in order to fill specified subplot areas
# IMPORTANT: must not include 'cbar' and 'cbar_label' in axis_labels
ar0 = 'auto'

period_dict = {'UVCet':0.2268,'UVCetMod':5.44709/24,'BLCet':0.2430,'BLCetMod':0.2430-0.00052}

periodname='UVCetMod'
    
period = period_dict[periodname]

func = real

ds_dir = '/data/jrv/uvcet_paper/' # where to save ds plots
if not os.path.exists(ds_dir):
    os.system('mkdir '+ds_dir)

# create figure
close('all')
figname = ds_dir+'UVCet_Lband_phased_dynspec.pdf'

# figure will have one plot per row - each row is one epoch Stokes V dynspec
fig,axes = subplots(4,1,sharex=True, figsize=(6.5,8.25))
fig.subplots_adjust(hspace=0,right=0.9,left=0.08,top=0.97,bottom=0.05)

plot_n = 1

for obsname in obs_list:

    # load dynamic spectra for this observation
    print '\n-----', obsname, '-----'
    obsfile,srcname = get_obsfile(obsname)
    ds = ds_list[obsfile]

    # bin LS band dynamic spectrum to desired resolution
    # mask RFI pix and chans before binning, pix after binning
    ds.mask_RFI_pixels(rmsfac=fig_params['pixflag_sigfacLS'],func=imag)
    ds.mask_RFI(rmsfac=fig_params['chanflag_sigfacLS'])

    # calculate phase (do again after binning since it doesn't auto-update)
    epoch = obsname[-1]
    convert_ds_time_to_phase(ds,period_target=periodname,epoch=epoch)
    print 'ds.phase[0]:',ds.phase[0]

    # pad with zeros to cover phase of 0.35 to 1.6 rotation periods
    phase_min = 0.35
    phase_max = 1.6
    #delta_phase_before = 2*ds.phase[0] - ds.phase[1] - phase_min # when I used ds.phase[0], sometimes the first phase would end up as 0.99
    delta_phase_before = ds.phase[0] - phase_min # when I used ds.phase[0], sometimes the first phase would end up as 0.99
    delta_phase_after = phase_max - ds.phase[-1]
    delta_time_before = period * delta_phase_before * 24 * 3600
    delta_time_after = period * delta_phase_after * 24 * 3600
    print ds.spec['i'].shape
    print 'Expanding tlims...'
    ds2 = ds.expand_tlims(t_add_left=delta_time_before,t_add_right = delta_time_after)
    print ds2.spec['i'].shape
    convert_ds_time_to_phase(ds2,period_target=periodname,epoch=epoch)
    print 'ds2.phase[0]:',ds2.phase[0]


    nt = int(round(fig_params['tint_LS']/ds2.dt()))    # number of integrations to bin together
    nf = int(round(fig_params['df_MHz_LS']/(ds2.df()/1e6))) # number of channels to bin together
    ds2 = ds2.bin_dynspec(nt=nt,nf=nf,mask_partial=fig_params['maskpartial_LS'])
    ds2.mask_RFI_pixels(rmsfac=fig_params['pixflag_sigfacLS'],func=imag)

    # calculate phase (needs to be done after all binning)
    epoch = obsname[-1]
    convert_ds_time_to_phase(ds2,period_target=periodname,epoch=epoch)
    print 'ds2.phase[0]:',ds2.phase[0]

    # plot Stokes V real, LS band
    subplot(4,1,plot_n)
    pp = {'pol':'v','smin':smin,'smax':smax,'trim_mask':False,'axis_labels':[],'ar0':ar0,'dy':0.5,'scale':fig_params['colorscale_LS'],'func':func,'xaxis_type':'phase'}
    plt,cbar_ticks,cbar_ticklbls = ds2.plot_dynspec(plot_params=pp)
    date = ds.t0().split()[0]
    title(date,x=0.1,y=0.8) #,fontsize='small')

    if plot_n == 2:
        ylabel('Frequency (GHz)',y=-0.1)
    elif plot_n == 4:
        xlabel('Rotational Phase (0 to 1 is full period)')
        suptitle('Period = ' + str(period*24) + ' hours',y=0.99)

    plot_n += 1

cbar_ax = fig.add_axes([0.92, 0.05, 0.025, 0.92])
cbar = fig.colorbar(plt,cax=cbar_ax)
cbar.set_ticks(cbar_ticks)
cbar.set_ticklabels(cbar_ticklbls)
savefig(figname)
