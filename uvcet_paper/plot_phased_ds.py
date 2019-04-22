'''
plot_allband_ds.py - Load P,L,S band dynamic spectrum for a given epoch, bin to specified resolution, and plot to file
'''

import dynspec.plot
reload(dynspec.plot)

from dynspec import load_dict
from dynspec.plot import *
from pylab import *
import os, subprocess
import matplotlib.gridspec as gridspec

'''
def get_obsname(obsfile):
    # take a file directory such as '/data/jrv/15A-416/YZCMi/1' and
    # convert to obs name such as '15A-416_YZCMi_1' and srcname 'YZCMi'
    names = obsfile.split('/')
    srcname = names[4]
    obsname = names[3]+'_'+names[4]+'_'+names[5]
    return obsname,srcname
'''

def get_obsfile(obsname):
    # take an obs name such as '15A-416_YZCMi_1' and return srcname ('YZCMi')
    # and file directory ('/data/jrv/15A-416/YZCMi/1')
    names = obsname.split('_')
    srcname = names[1]
    obsfile = '/data/jrv/'+names[0]+'/'+names[1]+'/'+names[2]
    return obsfile, srcname

def mjd_to_phase(t_mjd,period_target='UVCet',t0_mjd = 57207.5127083-0.12):
    '''
    mjd_to_phase_UVCet(t_mjd,period_target='UVCet',t0_mjd = 57207.5127083-0.12)
    - default t0_mjd is start time of 2015_UVCet_3 VLA dynspec minus roughly half a period
    
    Phase is scaled from 0 to 1 (where 1 is 360 degrees).
    
    This routine keeps the observation in time order - so if the observation starts at phase of 0.9,
    it may end at phase of 1.5.
    '''
    period_dict = {'UVCet':0.2268,'UVCetMod':0.2268+0.00014,'BLCet':0.2430,'BLCetMod':0.2430-0.00052}
    try:
        Prot = period_dict[period_target]
    except KeyError:
        print 'period_target is not in period_dict - using period of 1 day'
        Prot = 1.0
    t_rel = t_mjd - t0_mjd
    n_periods_off = floor(t_rel[0] / Prot)
    phase = t_rel / Prot - n_periods_off
    return phase

def convert_ds_time_to_phase(ds,period_target='UVCet',t0_mjd = 57207.5127083-0.2):
    t_mjd = ds.time.mjd
    phase = mjd_to_phase(t_mjd,period_target,t0_mjd)
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

# note: throughout, "LS" can also include C band, I initially wrote this code for 2015 data (which only has LS band)
#  but it works for the 2013 data with LSC band

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

#obs_list = fig_params_dict.keys()
#obs_list = ['15A-416_UVCet_1','15A-416_UVCet_2','15A-416_UVCet_3','15A-416_UVCet_4','15A-416_UVCet_5']
obs_list = ['15A-416_UVCet_2','15A-416_UVCet_3','15A-416_UVCet_4','15A-416_UVCet_5']

fig_max_width=6.5
fig_max_height=8.25
#tot_w = 6.5
#tot_h = 8.25

# set flux limits for LS band
smax = fig_params['smax_LS']
smin = -smax  # make colorbar symmetric about zero

# set axis ratio to 'auto' in order to fill specified subplot areas
# IMPORTANT: must not include 'cbar' and 'cbar_label' in axis_labels
ar0 = 'auto'
#ar0 = 0.3

'''
# calculate horizontal positions of subplots in units from 0 to 1
# (0 is left edge)
dsplot_w = 4.5  #* frac_duration  # width of dynamic spectrum in inches
gap_l = 0.55    # width of x-axis blank space (left) in inches
gap_cbar = 0.45 # width of blank space between V plot & cbar in inches
gap_r = 0.57     # width of x-axis blank space (right) in inches
cbar_w = 0.13   # width of colorbar in inches
tot_w = dsplot_w + cbar_w + gap_l + gap_cbar + gap_r # total width in inches

print 'Total width of figure in inches:', tot_w, '(goal: <=6.8)'
x1 = gap_l/tot_w             # left edge of Stokes V dynspec
x2 = x1 + dsplot_w/tot_w     # right edge of Stokes V dynspec
x3 = x2 + gap_cbar/tot_w     # left edge of colorbar
x4 = x3 + cbar_w/tot_w         # right edge of colorbar

# calculate vertical positions of subplots in units from 0 to 1
# (0 is bottom edge)
ds_h = 8.5     # * frac_BW   # height of LS band dynspec in inches
gap_t = 0.43   # height of y-axis blank space at top (includes titles) in inches
gap_rows = 0.5 # heights of each gap between rows of dynspecs in inches
gap_b = 0.36   # height of y-axis blank space at bottom in inches
tot_h = 4*ds_h + gap_t + 3*gap_rows + gap_b # total height in inches
print 'Total height of figure in inches:', tot_h, '(goal: <=8.5)'
y1 = 1-(gap_t/tot_h)    # top edge of row 1 dynspec (also top of cbar)
y2 = y1 - ds_h/tot_h    # bottom edge of row 1 dynspec
y3 = y2 - gap_rows/tot_h  # top edge of row 2 dynspec
y4 = y3 - ds_h/tot_h     # bottom edge of row 2 dynspec
y5 = y4 - gap_rows/tot_h  # top edge of row 3 dynspec
y6 = y5 - ds_h/tot_h     # bottom edge of row 3 dynspec
y7 = y6 - gap_rows/tot_h  # top edge of row 4 dynspec
y8 = y7 - ds_h/tot_h     # bottom edge of row 4 dynspec (also bottom of cbar)
cbar_h = (4*ds_h + 3*gap_rows)/tot_h
'''

period_dict = {'UVCet':0.2268,'UVCetMod':0.2268+0.00014,'BLCet':0.2430,'BLCetMod':0.2430-0.00052}

for periodname in ['UVCet','UVCetMod','BLCet','BLCetMod']:
#for periodname in ['UVCetMod']:
    
    period = period_dict[periodname]
    
    func = real
    
    ds_dir = '/data/jrv/uvcet_paper/phased_ds/' # where to save ds plots
    if not os.path.exists(ds_dir):
        os.system('mkdir '+ds_dir)

    # create figure
    close('all')
    figname = ds_dir+periodname+'.pdf'
    fig=figure(figsize=(6.5,8.25))
    # figure will have one plot per row - each row is one epoch Stokes V dynspec
    
    plot_n = 1
    
    for obsname in obs_list:

        # load dynamic spectra for this observation
        print '\n-----', obsname, '-----'
        obsfile,srcname = get_obsfile(obsname)
        ds = ds_list[obsfile]
        
        '''
        # Duration of observation relative to 3h40m (max duration of any) - scale x-axis by this
        #   so they are all on the same time scale
        duration = ds.get_tlist()[-1]*ds.dt()
        print 'Duration:',duration,'sec'
        frac_duration = duration/(3*3600+40*60)
        print 'Fractional duration compared to 3h40m:', frac_duration
        
        # Bandwidth of >1 GHz data relative to 3 GHz (default for 2015) - scale y-axis of >1 GHz dynspec by this
        BW_LSC = max(ds.f)-min(ds.f)
        frac_BW = BW_LSC/3.e9
        print 'Fractional bandwidth of >1 GHz data compared to 3 GHz:',frac_BW
        '''
        
        # bin LS band dynamic spectrum to desired resolution
        # mask RFI pix and chans before binning, pix after binning
        ds.mask_RFI_pixels(rmsfac=fig_params['pixflag_sigfacLS'],func=imag)
        ds.mask_RFI(rmsfac=fig_params['chanflag_sigfacLS'])
        nt = int(round(fig_params['tint_LS']/ds.dt()))    # number of integrations to bin together
        nf = int(round(fig_params['df_MHz_LS']/(ds.df()/1e6))) # number of channels to bin together
        ds = ds.bin_dynspec(nt=nt,nf=nf,mask_partial=fig_params['maskpartial_LS'])
        ds.mask_RFI_pixels(rmsfac=fig_params['pixflag_sigfacLS'],func=imag)
        
        # calculate phase (needs to be done after all binning)
        convert_ds_time_to_phase(ds,period_target=periodname)
        print 'ds.phase[0]:',ds.phase[0]
        
        # pad with zeros to cover phase of 0 to 2 rotation periods
        phase_min = 0.0
        phase_max = 2.0
        delta_phase_before = 2*ds.phase[0] - ds.phase[1] - phase_min # when I used ds.phase[0], sometimes the first phase would end up as 0.99
        delta_phase_after = phase_max - ds.phase[-1]
        delta_time_before = period * delta_phase_before * 24 * 3600
        delta_time_after = period * delta_phase_after * 24 * 3600
        print ds.spec['i'].shape
        print 'Expanding tlims...'
        ds2 = ds.expand_tlims(t_add_left=delta_time_before,t_add_right = delta_time_after)
        print ds2.spec['i'].shape
        convert_ds_time_to_phase(ds2,period_target=periodname)
        print 'ds2.phase[0]:',ds2.phase[0]
        
        # Format for axes command is axes([x_left, y_bottom, width, height])
        # First row: y_bottom is y2, x_left is x1, x3, x5

        # plot Stokes V real, LS band

        #ax = axes([x1,y2,dsplot_w/tot_w,dsLS_h/tot_h])
        #ax.set_autoscale_on(False)
        subplot(4,1,plot_n)
        pp = {'pol':'v','smin':smin,'smax':smax,'trim_mask':False,'axis_labels':[],'ar0':ar0,'dy':0.5,'scale':fig_params['colorscale_LS'],'func':func,'xaxis_type':'phase'}
        plt,cbar_ticks,cbar_ticklbls = ds2.plot_dynspec(plot_params=pp)
        #gca().xaxis.set_visible(False)
        #gca().yaxis.set_label_coords(-0.2,0)
        date = ds.t0().split()[0]
        title(date,x=0.1,y=0.8) #,fontsize='small')
        #fig.text(0.01,0.5,'Frequency (GHz)',va='center',rotation='vertical',fontsize='small')

        if plot_n == 2:
            ylabel('Frequency (GHz)',y=-0.2)
        elif plot_n == 4:
            xlabel('Rotational Phase (0 to 1 is full period)')
            suptitle('Period = ' + str(period*24) + ' hours',y=0.93)

        plot_n += 1
        
        '''
        # plot Stokes V real, LS band
        ax=axes([x3,y2,dsplot_w/tot_w,dsLS_h/tot_h])
        pp = {'pol':'v','smin':smin,'smax':smax,'trim_mask':False,'axis_labels':['xlabel'],'ar0':ar0,'dy':0.5,'scale':fig_params['colorscale_LS'],'func':func,'xaxis_type':'phase'}
        if fig_params['linthresh_LS']:
            pp['linthresh']=fig_params['linthresh_LS']
        ds.plot_dynspec(plot_params=pp)
        gca().yaxis.tick_right()
        xlabel_text = ax.xaxis.get_label_text()
        ax.set_xlabel('')
        #gca().xaxis.set_visible(False)
        if dsP:
            title('Stokes V, 1-4 GHz')
        else:
            title('Stokes V')

        
        # plot LS band colorbar
        ax = axes([x5,y2,cbar_w/tot_w,dsLS_h/tot_h])
        cbar=colorbar(plt,cax=ax)
        cbar.set_ticks(cbar_ticks)
        cbar.set_ticklabels(cbar_ticklbls)
        ax = cbar.ax
        if dsP:
            cbar_label = '1-4 Flux Density (mJy)'
            ycbar = 0.75
        else:
            cbar_label = 'Flux Density (mJy)'
            ycbar=0.65
            if obsname=='15A-416_UVCet_1':
                ycbar=0.98
        ax.text(4.2,ycbar,cbar_label,rotation=90,fontsize='small')
        '''
        
    #fig.text(0.5,0.01,xlabel_text,ha='center',fontsize='small')
    savefig(figname,bbox_inches='tight')
