'''
burst_rate.py

Purpose: Calculate instantaneous # of bursts per solid angle vs. flux for various transient surveys.
'''

from pylab import *
'''
close('all')
savedir = '/data/jrv/burst_paper/'

survey_list = ['VAST','ThunderKAT','VLASS']
survey = 'VAST'

savefile = savedir + survey + '_tseries.npy'
survey_dict = load(savefile).reshape(1)[0]  # weird song and dance to get it to load a dictionary properly...
times = survey_dict['times']
flux = survey_dict['flux']
flux_err = survey_dict['flux_err']
src = survey_dict['src']

fig1=figure(figsize=(9,6.5))
srclist = unique(src)
for srcname in srclist:
    ind = find([s == srcname for s in src])
    subplot(211)
    plot(ind*30,real(flux[ind]),'.')
    subplot(212)
    plot(ind*30,imag(flux[ind]),'.')
legend(srclist)
xlabel('Time in seconds (zero point is arbitrary)')
ylabel('Imaginary component')
smin = min([min(flux),min(imag(flux))])*1.1
axis([0,len(times)*30,smin,max(real(flux))*1.1])
subplot(211)
axis([0,len(times)*30,smin,max(real(flux))*1.1])
ylabel('Flux (Jy) at 1 pc')
filename1 = savedir + survey + '_tseries.png'


smin = median(flux_err) * 3 # minimum flux value for histogram
smax = max(real(flux)) * 1.1
# calculate # of integrations in time series with flux < S for a range of values of S
n_int,Srange,patches = hist(real(flux),range=(smin,smax),bins=1000,cumulative=-1)
n_int_imag,Srange,patches = hist(imag(flux),range=(smin,smax),bins=1000,cumulative=-1)
dN_dS = n_int / len(times)
dN_dS_imag = n_int_imag / len(times)
Srange = Srange[:-1]
fig2 = figure()
semilogx(Srange,dN_dS,'.-')
semilogx(Srange,dN_dS_imag,'.-')
i = max(find(dN_dS>0.1))
S10 = Srange[i]
axvline(S10,color='r')
legend(('Real','Imag'))
xlabel('Flux S (Jy) at 1 pc')
ylabel('Fraction of time brighter than S')
filename2 = savedir+survey+'_fluxdist.png'
savefig(filename2,bbox_inches='tight')


figure(fig1.number)
subplot(211)
axhline(S10,color='r')
savefig(filename1,bbox_inches='tight')
'''
close('all')

### Calculate N(>S) (per solid angle) vs. S ###

# Volume density of stars like our sample
n0 = 5.0/(4./3 * pi * 6.2**3) # our lower limit on volume density is our 5 stars divided by the volume containing them

# Range of flux we want on the plot
Splot = 10.0**arange(-4.,0.,0.1)
#Splot_contribution = array([0.0001,0.001,0.01,0.1,1.0])
Splot_contribution = array([0.0001,0.0003])

# Range of source distances we will consider
dr = 0.1  # pc
rmin = 0.1  # play with this to see how it affects result
rmax = sqrt(max(real(flux))/Splot[0]) # This is the furthest that the brightest burst could be while brighter than Splot[0]
r_range = arange(rmin,rmax,dr)

for S in Splot_contribution:
    S_thresh = S * r_range**2  # flux needed as a function of distance to be > S
    ind = [argmin(abs(Srange-St)) for St in S_thresh] # find index in Srange (the x-axis of the histogram of S at 1 pc) corresponding to S_thresh
    fraction_bursting = dN_dS[ind]  # fraction of stars like our sample that are bursting >S at a given instant
    n_bursting_in_shell = n0 * fraction_bursting * r_range**2 * dr  # number per solid angle in (r,r+dr)
    N = sum(n_bursting_in_shell)  # total # of stars per solid angle bursting at > S
    subplot(311)
    plot(r_range,S_thresh)
    subplot(312)
    plot(r_range,fraction_bursting)
    subplot(313)
    plot(r_range,n_bursting_in_shell/N)
subplot(311)
legend(Splot_contribution*1e3,title="Smin (mJy)",loc='upper left')
ylabel('S_1pc (Jy)')
subplot(312)
ylabel('Fraction >S_1pc')
subplot(313)
ylabel('Fractional contribution')
xlabel('Distance (pc)')
filename3 = savedir+survey+'_contribution.png'
savefig(filename3,bbox_inches='tight')

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

# calculate N~S^-1.5 for comparison
Neuclid = N_S[0] * (Splot/Splot[0])**-1.5
dN = N_S - Neuclid

figure()
subplot(211)
loglog(Splot*1e3,N_S,'.')
loglog(Splot*1e3,Neuclid)
ylabel('N(>S) (per sq deg)')
subplot(212)
semilogx(Splot,dN/Neuclid * 100.,'.-')
ylabel('% deviation')
xlabel('S (mJy)')
filename4 = savedir + survey + '_NvS.png'
savefig(filename4,bbox_inches='tight')

close('all')
