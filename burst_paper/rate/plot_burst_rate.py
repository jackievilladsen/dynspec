'''
plot_burst_rate.py

Purpose: Load burst rates from calc_burst_rate.py and generate a series of plots
         comparing the 3 transient surveys.
'''

from pylab import *

close('all')

savedir = '/data/jrv/burst_paper/rate/'

survey_list = ['VLASS','Lband','Pband']
leg_list = ['S band','L band','P band']
plot_file_root = 'band_'

survey_list = ['VAST','ThunderKAT','VLASS','Lband','Pband']
leg_list = ['L band/VAST','L band/ThunderKAT','S band/VLASS','L band','P band']
plot_file_root = 'band_'

survey_list = ['Phi','Llo','Lhi','Slo','Shi']
leg_list = ['340-480 MHz','1-1.4 GHz','1.4-2 GHz','2-2.8 GHz','2.8-4 GHz']
plot_file_root = 'halfband_'

params = {'legend.fontsize': 'small',
          'axes.titlesize': 'medium',
          'axes.labelsize': 'small',
          'xtick.labelsize': 'x-small',
          'ytick.labelsize': 'x-small'}
rcParams.update(params)

n = len(survey_list)
### PLOT 1: TIME SERIES ###

figure(figsize=(6.5,8))
j=0
for survey,survey_leg in zip(survey_list,leg_list):
    j+=1
    subplot(n,1,j)
    
    savefile = savedir + survey + '_burstrate.npy'
    survey_dict = load(savefile).reshape(1)[0]  # weird song and dance to get it to load a dictionary properly...
    times = survey_dict['times']
    flux = survey_dict['flux']
    smin = survey_dict['smin']
    t_int = survey_dict['t_int']
    print survey, 'smin (red line on tseries.png):', smin*1e3, 'mJy'
    src = survey_dict['src']
    srclist = unique(src)
    
    imag_offset = 1.0
    if survey[0]=='P':
        imag_offset = 1.0
    
    for srcname in srclist:
        ind = find([s == srcname for s in src])
        plot(ind*t_int,real(flux[ind]),'.')
    gca().set_prop_cycle(None)
    for srcname in srclist:
        ind = find([s == srcname for s in src])
        plot(ind*t_int,imag(flux[ind])-imag_offset,'.')
    if j==n:
        legend(srclist,loc='upper left',fontsize='x-small')
        xlabel('Time in seconds (zero point is arbitrary)')
    elif j==(n/2+1):
        ylabel('Flux (Jy) at 1 pc')
    title(survey_leg)
    axis([0,len(times)*t_int,-imag_offset-0.3,0.8])
    axhline(0.,color='k')
    axhline(-imag_offset,color='k')
    axhline(smin,color='r')
    axhline(smin-imag_offset,color='r')

tight_layout()
filename1 = savedir + plot_file_root+'tseries.pdf'
savefig(filename1,bbox_inches='tight')

### PLOT 2: Fraction of time 1-pc flux is > S ###

# load all data first so I can plot it in right order
dN_dS = {}
dN_dS_imag = {}
Srange={}
for survey in survey_list:
    savefile = savedir + survey + '_burstrate.npy'
    survey_dict = load(savefile).reshape(1)[0]  # weird song and dance to get it to load a dictionary properly...
    Srange[survey] = survey_dict['Srange']
    dN_dS[survey] = survey_dict['dN_dS']
    dN_dS_imag[survey] = survey_dict['dN_dS_imag']

figure(figsize=(6.5,3))
subplot(121)
for survey in survey_list:
    semilogx(Srange[survey],dN_dS[survey],'.-')
#legend(leg_list)
gca().set_prop_cycle(None)
#for survey in survey_list:
#    semilogx(Srange[survey],dN_dS_imag[survey],'-')
axis([smin,2,0,0.3])
xlabel('Flux S (Jy) at 1 pc')
ylabel('Fraction of time brighter than S')
#filename2 = savedir+plot_file_root+'fluxdist.pdf'
#savefig(filename2,bbox_inches='tight')

### PLOT 3: N(>S) vs. S for each survey band ###

# load data
N_S = {}
Neuclid = {}
for survey in survey_list:
    savefile = savedir + survey + '_burstrate.npy'
    survey_dict = load(savefile).reshape(1)[0]  # weird song and dance to get it to load a dictionary properly...
    Splot = survey_dict['Splot']*1e3  # same for all
    N_S[survey] = survey_dict['N_S']
    Neuclid[survey] = N_S[survey][0] * (Splot/Splot[0])**-1.5

#figure(figsize=(4,3))
subplot(122)
for survey in survey_list:
    loglog(Splot,Neuclid[survey])
gca().set_prop_cycle(None)
for survey in survey_list:
    loglog(Splot,N_S[survey],'.')
legend(leg_list)
xlabel('S (mJy)')
ylabel('N(>S) (per sq deg)')
axis([0.1,1e3,1e-5,1])
#filename3 = savedir + 'survey_NvS.pdf'
#savefig(filename3,bbox_inches='tight')
fname = savedir + plot_file_root+'burstrate.pdf'
tight_layout()
savefig(fname,bbox_inches='tight')


### PLOT 4: CONTRIBUTION VS. DISTANCE ###

Smin = {}
r_range = {}
f_burst = {}
n_burst = {}
for survey in survey_list:
    savefile = savedir + survey + '_burstrate.npy'
    survey_dict = load(savefile).reshape(1)[0]  # weird song and dance to get it to load a dictionary properly...
    Smin = survey_dict['Splot_contribution'] # same for all
    r_range[survey] = survey_dict['r_range']
    f_burst[survey] = survey_dict['fraction_bursting']
    n_burst[survey] = survey_dict['n_bursting_per_pc']

figure(figsize=(6.5,3))
for survey in survey_list:
    subplot(121)
    plot(r_range[survey],f_burst[survey])
    subplot(122)
    plot(r_range[survey],n_burst[survey])
subplot(121)
xlabel('Distance (pc)')
ylabel('Fraction of stars bursting >Smin')
legend(leg_list)
subplot(122)
xlabel('Distance (pc)')
ylabel('# bursting stars on sky (per pc)')
suptitle('Contribution vs. distance for Smin='+str(Smin*1e3)+' mJy',y=1.03)
tight_layout()
fname = savedir + plot_file_root+'contribution.pdf'
savefig(fname,bbox_inches='tight')

close('all')
