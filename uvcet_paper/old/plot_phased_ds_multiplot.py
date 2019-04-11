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

def mjd_to_phase(t_mjd,period_target='UVCet',t0_mjd = 57207.5127083):
    '''
    mjd_to_phase_UVCet(t_mjd,period_target='UVCet',t0_mjd = 57207.5127083)
    - default t0_mjd is start time of 2015_UVCet_3 VLA dynspec
    
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

def convert_ds_time_to_phase(ds,period_target='UVCet',t0_mjd = 57207.5127083):
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

# params that can be changed are listed in default_fig_params
default_fig_params = {
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

fig_params_dict = {
    '13A-423_UVCet_1':{},
    '13A-423_UVCet_2':{},
    '15A-416_UVCet_1':{},
    '15A-416_UVCet_2':{},
    '15A-416_UVCet_3':{},
    '15A-416_UVCet_4':{},
    '15A-416_UVCet_5':{}
}


### PLOT INDIVIDUAL OBSERVATIONS ###

obs_list = fig_params_dict.keys()
#obs_list = ['15A-416_UVCet_1','15A-416_UVCet_3'] # so I can work on just this event

fig_max_width=6.5
fig_max_height=8.25

for obsname in obs_list:

    for periodname in ['UVCet','UVCetMod','BLCet','BLCetMod']:
        func = real

        ds_dir = '/data/jrv/uvcet_paper/phased_ds/'+periodname+'/' # where to save ds plots
        if not os.path.exists(ds_dir):
            os.system('mkdir '+ds_dir)

        # load dynamic spectra for this observation
        print '\n-----', obsname, '-----'
        obsfile,srcname = get_obsfile(obsname)
        ds = ds_list[obsfile]
#        dsP = dsP_list.get(obsfile,None)

        # load custom parameters for plotting this epoch (binning, RFI flagging, color scale)
        fig_params = deepcopy(default_fig_params)
        fp_dict_temp = fig_params_dict.get(obsname,{})
        for k in fp_dict_temp:
             fig_params[k] = fp_dict_temp[k]

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

        '''
        if dsP:
            dsP.mask_RFI_pixels(rmsfac=fig_params['pixflag_sigfacP'])
            dsP.mask_RFI(rmsfac=fig_params['chanflag_sigfacP'])
            # bin P band dynamic spectrum to desired resolution
            nt = int(round(fig_params['tint_P']/dsP.dt()))    # number of integrations to bin together
            nf = int(round(fig_params['df_MHz_P']/(dsP.df()/1e6))) # number of channels to bin together
            dsP = dsP.bin_dynspec(nt=nt,nf=nf,mask_partial=fig_params['maskpartial_P'])
            dsP.mask_RFI(rmsfac=fig_params['chanflag_sigfacP'])
        '''
        
        # calculate horizontal positions of subplots in units from 0 to 1
        # (0 is left edge)
        dsplot_w = 3.2 * frac_duration  # width of dynamic spectrum in inches
        gap_l = 0.55    # width of x-axis blank space (left) in inches
        gap_c = 0.15    # width of x-axis blank space (center) in inches
        gap_cbar = 0.45 # width of blank space between V plot & cbar in inches
        gap_r = 0.57     # width of x-axis blank space (right) in inches
        cbar_w = 0.13   # width of colorbar in inches
        tot_w = 2*dsplot_w + cbar_w + gap_l + gap_c + gap_cbar + gap_r # total width in inches
        #if obs == '13A-423_UVCet_2':
        #    tot_w += gap_c + dsplot_w + gap_cbar + gap_r
        print 'Total width of figure in inches:', tot_w, '(goal: <=8.25)'
        x1 = gap_l/tot_w             # left edge of Stokes I dynspec
        x2 = x1 + dsplot_w/tot_w     # right edge of Stokes I dynspec
        x3 = x2 + gap_c/tot_w        # left edge of Stokes V dynspec
        x4 = x3 + dsplot_w/tot_w     # right edge of Stokes V dynspec
        x5 = x4 + gap_cbar/tot_w     # left edge of colorbar
        x6 = x5+cbar_w/tot_w         # right edge of colorbar
        #if obs == '13A-423_UVCet_2':
        #    x7 = x6 + (gap_r+gap_c)/tot_w   # left edge of second Stokes V dynspec
        #    x8 = x

        # calculate vertical positions of subplots in units from 0 to 1
        # (0 is bottom edge)
        dsLS_h = 3.2 * frac_BW   # height of LS band dynspec in inches
        dsP_h = 0.9    # height of P band dynspec in inches
        gap_t = 0.43   # height of y-axis blank space at top (includes titles) in inches
        gap_rows = 0.5 # heights of each gap between rows of dynspecs in inches
        gap_b = 0.36   # height of y-axis blank space at bottom in inches
        if dsP:
            tot_h = dsLS_h + 2*dsP_h + gap_t + 2*gap_rows + gap_b # total height in inches
        else:
            tot_h = gap_t + dsLS_h + gap_b    # total height in inches if no P band data
        print 'Total height of figure in inches:', tot_h, '(goal: <=6.8)'
        y1 = 1-(gap_t/tot_h)    # top edge of LS band dynspec
        y2 = y1 - dsLS_h/tot_h    # bottom edge of LS band dynspec
        y3 = y2 - gap_rows/tot_h  # top edge of P band I,V dynspecs
        y4 = y3 - dsP_h/tot_h     # bottom edge of P band I,V dynspecs
        y5 = y4 - gap_rows/tot_h  # top edge of P band U dynspec
        y6 = y5 - dsP_h/tot_h     # bottom edge of P band U dynspec
        cbarP_h = (2*dsP_h + gap_rows)/tot_h

        # create figure
        close('all')
        figname = ds_dir+obsname+'.pdf'
        if func == imag:
            figname = ds_dir+obsname+'_imag.pdf'
        fig=figure(figsize=(tot_w,tot_h))

        # First row of plots: Stokes I LS, Stokes V LS, colorbar LS
        # Format for axes command is axes([x_left, y_bottom, width, height])
        # First row: y_bottom is y2, x_left is x1, x3, x5

        # set flux limits for LS band
        smax = fig_params['smax_LS']
        if smax is None:
            smax = max(percentile(real(ds.spec['i']),99)*1.1,median(real(ds.spec['i']))*2)
        smin = -smax  # make colorbar symmetric about zero

        # set axis ratio to 'auto' in order to fill specified subplot areas
        # IMPORTANT: must not include 'cbar' and 'cbar_label' in axis_labels
        ar0 = 'auto'

        # plot Stokes I real, LS band
        ax = axes([x1,y2,dsplot_w/tot_w,dsLS_h/tot_h])
        #ax.set_autoscale_on(False)
        pp = {'pol':'i','smin':smin,'smax':smax,'trim_mask':False,'axis_labels':[],'ar0':ar0,'dy':0.5,'scale':fig_params['colorscale_LS'],'func':func,'xaxis_type':'phase'}
        if fig_params['linthresh_LS']:
            pp['linthresh']=fig_params['linthresh_LS']
        plt,cbar_ticks,cbar_ticklbls = ds.plot_dynspec(plot_params=pp)
        #gca().xaxis.set_visible(False)
        #gca().yaxis.set_label_coords(-0.2,0)
        if dsP:
            title('Stokes I, 1-4 GHz')
        else:
            title('Stokes I')
        fig.text(0.01,0.5,'Frequency (GHz)',va='center',rotation='vertical',fontsize='small')

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
        if dsP:
            # Second row of plots: Stokes I P, apparent Stokes V P
            # Format for axes command is axes([x_left, y_bottom, width, height])
            # Second row: y_bottom is y4, x_left is x1, x3    

            # set flux limits for P band
            smaxP = fig_params['smax_P']
            if smaxP is None:
                smaxP = dsP.get_rms('v')*6.
            sminP = -smaxP

            # plot Stokes I real, P band
            ax = axes([x1,y4,dsplot_w/tot_w,dsP_h/tot_h])
            pp = {'pol':'i','smin':sminP,'smax':smaxP,'trim_mask':False,'axis_labels':[],'dy':0.05,'ar0':ar0,'scale':fig_params['colorscale_P'],'func':func}
            if fig_params['linthresh_P']:
                pp['linthresh']=fig_params['linthresh_P']
            dsP.plot_dynspec(plot_params=pp)
            title('Stokes I, 0.2-0.5 GHz')

            # plot Stokes V real, P band
            ax = axes([x3,y4,dsplot_w/tot_w,dsP_h/tot_h])
            pp = {'pol':'v','smin':sminP,'smax':smaxP,'trim_mask':False,'axis_labels':[],'dy':0.05,'ar0':ar0,'scale':fig_params['colorscale_P'],'func':func}
            if fig_params['linthresh_P']:
                pp['linthresh']=fig_params['linthresh_P']
            plt,cbar_ticks,cbar_ticklbls=dsP.plot_dynspec(plot_params=pp)
            gca().yaxis.tick_right()
            title('Stokes V\', 0.2-0.5 GHz')

            # Third row of plots: [empty], apparent Stokes U P, P band colorbar (extra height)
            # Format for axes command is axes([x_left, y_bottom, width, height])
            # Third row: y_bottom is y6
            #            x_left is x3 (Stokes U), x5 (colorbar)
            #            height is dsP_h (Stokes U), 2*dsP_h+gap_rows (colorbar)

            # plot Stokes U real, P band
            ax = axes([x3,y6,dsplot_w/tot_w,dsP_h/tot_h])
            pp = {'pol':'u','smin':sminP,'smax':smaxP,'trim_mask':False,'axis_labels':[],'dy':0.05,'ar0':ar0,'scale':fig_params['colorscale_P'],'func':func}
            if fig_params['linthresh_P']:
                pp['linthresh']=fig_params['linthresh_P']
            dsP.plot_dynspec(plot_params=pp)
            gca().yaxis.tick_right()
            title('Stokes U\', 0.2-0.5 GHz')

            # plot P band colorbar
            ax = axes([x5,y6,cbar_w/tot_w,cbarP_h])
            cbar=colorbar(plt,cax=ax)
            cbar.set_ticks(cbar_ticks)
            cbar.set_ticklabels(cbar_ticklbls)
            ax = cbar.ax
            ax.text(4.2,0.9,'0.2-0.5 GHz Flux Density (mJy)',rotation=90,fontsize='small')
        '''
        
        fig.text(0.5,0.01,xlabel_text,ha='center',fontsize='small')
        date = ds.t0().split()[0]
        fig_title = srcname[0:2]+' '+srcname[2:5]+' - '+date
        if func == imag:
            fig_title += ' - Imag(vis)'
        suptitle(fig_title,y=0.99,fontsize='medium')
        savefig(figname)
