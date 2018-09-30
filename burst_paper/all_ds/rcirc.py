from pylab import *

obs_dict = {
    'ADLeo_3P': {'I':5.01, 'sigI':0.61, 'UV':4.84, 'sigUV': 0.21},
    'ADLeo_5P': {'I':15.8, 'sigI':4.9, 'UV':33.4, 'sigUV': 1.2},
    'UVCet_4P_short': {'I':31.7, 'sigI':2.6, 'UV':28.2, 'sigUV': 1.0},
    'EQPeg_2P': {'I':43.7, 'sigI':7.2, 'UV':39.0, 'sigUV': 5.4}}

# throughout this code:
# rc = degree of circular polarization
# S = Stokes I flux density ("true" value)
# V, I or {d} = measured flux in Stokes I & V, respectively
# sigV,sigI = rms noise on observations (assumed same for Stokes I & V)

# conf = 0.68 is 1-sigma confidence interval, 0.9 should give [rc_LL,rc_UL] if conf_UL is 0.95
conf = 0.68 # desired confidence level for confidence interval on rc
conf_UL = 0.68

for obs in obs_dict:
    
    print '\n-------',obs,'--------'
    
    I = obs_dict[obs]['I']
    V = obs_dict[obs]['UV']
    sigI = obs_dict[obs]['sigI']
    sigV = obs_dict[obs]['sigUV']

    # range of params to search
    rcmin = 0.0
    rcmax = 1.0
    nrc = 1001

    Smin = 0.0
    Smax = (max(abs(V),I)+sigI*3)*2# LL on a doesn't depend much on this - good
    nS = 1001 # LL on a doesn't depend much on this - good

    # 2-D posterior distribution
    # p(rc,S|{d})
    # A is an as-yet-unknown normalizing factor (use 1 until I figure it out)
    def posterior(rc,S,A):
        return A * exp( -1 * (V - S*rc)**2 / (2*sigV**2) - (I - S)**2 / (2*sigI**2) )

    # create vectors of param values for grid, then meshgrid them
    drc = (rcmax - rcmin) / (nrc - 1)
    rcrange = arange(rcmin,rcmax+drc,drc)

    dS = (Smax - Smin) / (nS - 1)
    Srange = arange(Smin,Smax+dS,dS)

    # rcmat varies horizontally (all rows are identical)
    # Smat varies vertically (all columns are identical)
    rcmat,Smat=meshgrid(rcrange,Srange)

    # calculate posterior (un-normalized)
    p2d_unnorm = posterior(rcmat,Smat,1)

    # calculate normalizing factor
    A = 1/sum(p2d_unnorm)

    # normalized posterior
    p2d = p2d_unnorm * A

    # 3-D mesh plot of the normalized 2-D posterior
    figure()
    subplot(1,2,1)
    imshow(p2d,extent=[rcmin,rcmax,Smin,Smax],aspect='auto',origin='lower') # figure out how to display axes later
    suptitle(obs)
    xlabel('Degree of circ pol')
    ylabel('Flux Density (mJy)')

    # marginalize over S2 to get p1d = p(rc|{d})
    p1d = sum(p2d,0)

    # line plot of p(rc|{d}) - posterior pdf for spectral index
    #figure()
    subplot(2,2,2)
    plot(rcrange,p1d)
    ymax = max(p1d)*1.1
    axis([rcmin,rcmax,0,ymax])
    title('Posterior PDF')

    # calculate CDF for spectral index: P(a|{d})
    cdf1d = cumsum(p1d)
    #subplot(2,1,2)
    #plot(rcrange,cdf1d)

    # figure out confidence interval
    # 1) find value of rc that gives CDF nearest to 2.5% (if conf=0.95)
    ind = argmin(abs(cdf1d-(1-conf)/2))
    rc_low = rcrange[ind]
    axvline(rc_low,color='k',linestyle='--')

    # 2) find value of rc that gives CDF nearest to 97.5% (for 95% conf)
    ind = argmin(abs(cdf1d-(1-(1-conf)/2)))
    rc_hi = rcrange[ind]
    axvline(rc_hi,color='k',linestyle='--')

    # line plot of CDF
    subplot(2,2,4)
    plot(rcrange,cdf1d)
    axvline(rc_low,color='k',linestyle='--')
    axvline(rc_hi,color='k',linestyle='--')
    axhline((1-conf)/2,color='k',linestyle=':')
    axhline(1-(1-conf)/2,color='k',linestyle=':')
    axis([rcmin,rcmax,0,1])
    xlabel('Degree of circ pol')
    title('Cumulative CDF')

    # find max of PDF
    ind = argmax(p1d)
    rc_peak = rcrange[ind]

    print str(conf*100)+'% confidence interval: ['+str(rc_low)+','+str(rc_hi)+']'
    diff_hi = rc_hi - rc_peak
    diff_low = rc_peak - rc_low
    print 'Other way:',rc_peak, '+', diff_hi, '-', diff_low
    if diff_hi < 0 or diff_low < 0:
        print 'Warning! Most probable rc is outside symmetric confidence interval, use lower or upper limit instead.'

    # find upper and lower limits for 95%
    ind_LL = argmin(abs(cdf1d-(1-conf_UL)))
    ind_UL = argmin(abs(cdf1d-conf_UL))
    rc_LL = rcrange[ind_LL]
    rc_UL = rcrange[ind_UL]
    print conf_UL*100,'% lower limit: >',rc_LL
    print conf_UL*100,'% upper limit: <',rc_UL

show()
