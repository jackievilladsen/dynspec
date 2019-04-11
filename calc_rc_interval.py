from pylab import *

# throughout this code:
# rc = degree of circular polarization
# S = Stokes I flux density ("true" value)
# V, I or {d} = measured flux in Stokes I & V, respectively
# sigV,sigI = rms noise on observations (assumed same for Stokes I & V)

# 2-D posterior distribution
# p(rc,S|{d})
# A is an as-yet-unknown normalizing factor (use 1 until I figure it out)
def posterior(rc,S,A,V,sigV,I,sigI):
    return A * exp( -1 * (V - S*rc)**2 / (2*sigV**2) - (I - S)**2 / (2*sigI**2) )

def rc_conf_int(IV_dict,conf = 0.68):
    I = IV_dict['I']
    V = IV_dict['V']
    sigI = IV_dict['Ierr']
    sigV = IV_dict['Verr']
    
    # range of params to search
    rcmin = -1.0
    rcmax = 1.0
    nrc = 1001
    
    Smin = 0.0
    Smax = (max(abs(V),I)+sigI*3)*2# LL on a doesn't depend much on this - good
    nS = 1001 # LL on a doesn't depend much on this - good
    
    # create vectors of param values for grid, then meshgrid them
    drc = (rcmax - rcmin) / (nrc - 1)
    rcrange = arange(rcmin,rcmax+drc,drc)

    dS = (Smax - Smin) / (nS - 1)
    Srange = arange(Smin,Smax+dS,dS)

    # rcmat varies horizontally (all rows are identical)
    # Smat varies vertically (all columns are identical)
    rcmat,Smat=meshgrid(rcrange,Srange)

    # calculate posterior (un-normalized)
    p2d_unnorm = posterior(rcmat,Smat,1,V,sigV,I,sigI)

    # calculate normalizing factor
    A = 1/sum(p2d_unnorm)

    # normalized posterior
    p2d = p2d_unnorm * A
    
    # marginalize over S2 to get p1d = p(rc|{d})
    p1d = sum(p2d,0)
    
    # calculate CDF for spectral index: P(a|{d})
    cdf1d = cumsum(p1d)
    
    # figure out confidence interval
    # 1) find value of rc that gives CDF nearest to 2.5% (if conf=0.95)
    ind = argmin(abs(cdf1d-(1-conf)/2))
    rc_low = rcrange[ind]

    # 2) find value of rc that gives CDF nearest to 97.5% (for 95% conf)
    ind = argmin(abs(cdf1d-(1-(1-conf)/2)))
    rc_hi = rcrange[ind]

    # find max of PDF
    ind = argmax(p1d)
    rc_peak = rcrange[ind]

    # if best-fit rc is outside CI range, change CI range to lower/upper limit
    if rc_peak > rc_hi:
        rc_hi = 1.0
        ind_LL = argmin(abs(cdf1d-(1-conf)))
        rc_low = rcrange[ind_LL]
    elif rc_peak < rc_low:
        rc_low = -1.0
        ind_UL = argmin(abs(cdf1d-conf))
        rc_hi = rcrange[ind_UL]

    return rc_low, rc_hi, rc_peak
