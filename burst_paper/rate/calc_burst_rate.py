'''
calc_burst_rate.py

Purpose: Calculate instantaneous # of bursts per solid angle vs. flux for various transient surveys.

Run compile_tseries.py before this.

'''

from pylab import *
import pickle

close('all')
savedir = '/data/jrv/burst_paper/rate/'

#survey_list = ['VAST','ThunderKAT','VLASS']
survey_list = ['Phi','Llo','Lhi','Slo','Shi'] #,'Clo']

for survey in survey_list:

    ### Load 1-pc flux time series ###

    savefile = savedir + survey + '_tseries.npy'
    survey_dict = load(savefile).reshape(1)[0]  # weird song and dance to get it to load a dictionary properly...
    flux = survey_dict['flux']
    flux_err = survey_dict['flux_err']


    ### Calculate fraction of time that 1-pc flux is > S ###

    # flux range to consider: from 3*RMS up to maximum 1-pc flux detected
    smin = median(flux_err) * 1.5 # minimum flux value for histogram
    smin = 0.1
    smax = max(real(flux)) * 1.1
    print 'Median 3-sigma for',survey,'time series:',smin*1e3,'mJy'

    # calculate # of integrations in time series with flux > S for a range of values of S
    n_int,Srange,patches = hist(real(flux),range=(smin,smax),bins=1000,cumulative=-1)
    n_int_imag,Srange,patches = hist(imag(flux),range=(smin,smax),bins=1000,cumulative=-1) # do it for imaginary tseries to check false alarm rate
    Srange = Srange[:-1]  # drop right edge of last bin so Srange contains the left edges of bins

    # convert to fraction of time bursting > S
    dN_dS = n_int / len(flux)
    dN_dS_imag = n_int_imag / len(flux)
    
    '''
    i = max(find(dN_dS>0.1))
    S10 = Srange[i]
    print 'In', survey, 'band, source flux at 1 pc spends 10% of time at >',S10*1e3,'mJy'
    '''

    ### Calculate N(>S) (per sq deg) vs. S ###

    # Volume density of stars like our sample
    n0 = 5.0/(4./3 * pi * 6.2**3) # conservative lower limit on volume density is our 5 stars divided by the volume containing them

    # Range of flux S for which we want to calculate N(>S)
    Splot = 10.0**arange(-4.,0.,0.1)
    Splot_contribution = array([0.0001,0.0003]) # make detailed plot showing which distances contribute most for this flux limit

    # Range of source distances we will consider
    dr = 0.1  # pc
    rmin = 0.1  # I tried 2.6 pc also - result depends on this only weakly
    rmax = sqrt(max(real(flux))/Splot[0]) # This is the furthest that the brightest burst could be while brighter than Splot[0]
    r_range = arange(rmin,rmax,dr)

    j=0
    for S in Splot_contribution:
        S_thresh = S * r_range**2  # S_thresh: 1-pc flux needed as a function of distance to be > S
        ind = [argmin(abs(Srange-St)) for St in S_thresh] # find index in Srange (the x-axis of the histogram of S at 1 pc) corresponding to S_thresh
        fraction_bursting = dN_dS[ind]  # fraction of stars like our sample that are bursting >S_thresh at a given instant
        n_bursting_in_shell = 4*pi*n0 * fraction_bursting * r_range**2 * dr  # number bursting >S_thresh per solid angle in (r,r+dr)
        N = sum(n_bursting_in_shell)  # total # of stars per solid angle bursting at > S
        subplot(311)
        plot(r_range,S_thresh)
        subplot(312)
        plot(r_range,fraction_bursting)
        subplot(313)
        plot(r_range,n_bursting_in_shell/N)
        if j==0:
            survey_dict['Splot_contribution'] = S
            survey_dict['fraction_bursting'] = fraction_bursting
            survey_dict['n_bursting_per_pc'] = n_bursting_in_shell/dr
            survey_dict['r_range'] = r_range
        j += 1
    subplot(311)
    legend(Splot_contribution*1e3,title="Smin (mJy)",loc='upper left')
    ylabel('S_1pc (Jy)')
    subplot(312)
    ylabel('Fraction >S_1pc')
    subplot(313)
    ylabel('Fractional contribution')
    xlabel('Distance (pc)')
    fname = savedir+survey+'_contribution.png'
    savefig(fname,bbox_inches='tight')

    N_S = []
    sqdeg_per_ster = (180./pi)**2
    for S in Splot:
        S_thresh = S * r_range**2  # flux needed as a function of distance to be > S
        ind = [argmin(abs(Srange-St)) for St in S_thresh] # find index in Srange (the x-axis of the histogram of S at 1 pc) corresponding to S_thresh
        fraction_bursting = dN_dS[ind]  # fraction of stars like our sample that are bursting >S at a given instant
        n_bursting_in_shell = n0 * fraction_bursting * r_range**2 * dr / sqdeg_per_ster  # number per solid angle in (r,r+dr)
        N = sum(n_bursting_in_shell)  # total # of stars per solid angle bursting at > S
        N_S.append(N)
    N_S = array(N_S)

    # save variables generated here that are needed for 
    survey_dict['smin'] = smin
    survey_dict['Srange'] = Srange
    survey_dict['dN_dS'] = dN_dS
    survey_dict['dN_dS_imag'] = dN_dS_imag
    survey_dict['Splot'] = Splot
    survey_dict['N_S'] = N_S

    savefile = savedir + survey + '_burstrate.npy'
    save(savefile,survey_dict)

    close('all')

