'''
analyze_MonteCarlo_results.py:

Load the pgram and peak file outputs from the NASA Exoplanet Archive Periodogram Service,
for the Monte Carlo trials with different randomized observation times, then analyze and plot
the statistics of periodogram peaks to determine the significance of the peak in the periodogram
of the real data.
'''



##### USER DEFINED PARAMETERS #####

# set reload=False if I have already run this script once in ipython and don't want to reload the text files (saves time as I iterate in plotting/data analysis)
reload=True

# change this directory if I create different MonteCarlo trials in another directory
MC_dir = '/export/data_1/jvillads/UVCet_paper/pgram/new_tseries_RR/MonteCarlo/pgrams/'

# how many sigma around Barnes+2017 UV Cet period should we plot?
Nsig_barnes = 3

# number of rows and cols of subplots for plotting periodogram around Barnes UV Cet period
#nrows = 6
#ncols = 6

##### END USER DEFINED PARAMETERS #####

# parameters of periods from Barnes+2017 (units of hours)
period_uvcet_barnes = 5.4432
sigP_uvcet_barnes = 0.0072
period_blcet_barnes = 5.832
sigP_blcet_barnes = 0.012

from pylab import *
from numpy import *
from scipy import signal
import os, glob

# function to go between trial # (-1 for real data) and trial_string, which is appropriate title for trial
#   and key for dictionaries
def get_trial_string(i):
    if i < 0: # i = -1 is reserved for real data
        return 'Observed data'
    else:    # i >=0 are Monte-Carlo trials
        return 'Trial #'+str(i)


##### MAIN SCRIPT STARTS HERE #####

os.chdir(MC_dir)

# check how many trials have been run through the pgram and are in this directory
trial_pgram_filelist = glob.glob('trial*_pgram.tbl')
trial_n_list = sort([int(f.strip('trial').strip('_pgram.tb')) for f in trial_pgram_filelist])
n_trials = amax(trial_n_list)+1         # number of trials with saved data files in MC_dir

#n_trials = 35

if reload:
    
    # containers to compile pgram data for all trials once loaded
    pgram_power_dict = {}
    peak_dict = {}
    peak_power_list = []
    local_peak_dict = {}
    local_peak_power_list = []
    
    for i in range(-1,n_trials):

        trial_string = get_trial_string(i)
        if i < 0: # i = -1 is reserved for real data
            pgram_file = 'realdata_pgram.tbl'
            peak_file = 'realdata.dat.top'
        else:    # i >=0 are Monte-Carlo trials
            pgram_file = 'trial'+str(i)+'_pgram.tbl'
            peak_file = 'trial'+str(i)+'.dat.top'

        # load periodogram plot and list of top 50 peaks
        if not os.path.exists(pgram_file) or not os.path.exists(peak_file):
            print 'WARNING: One or more of these pgram files for', trial_string, 'does not exist, skipping this trial:',pgram_file+', '+peak_file
            continue
        print 'Loading data for '+trial_string+'...'
        period,power = loadtxt(pgram_file,skiprows=3,usecols=(1,2),unpack=True)
        period_peak,power_peak = loadtxt(peak_file,skiprows=22,usecols=(1,2),unpack=True)

        # plot periodogram in range w/in Nsig_barnes sigma of Barnes UV Cet period
        period_sig = (period - period_uvcet_barnes) / sigP_uvcet_barnes
        ind = find( abs(period_sig) < Nsig_barnes )

        # find peaks in this range
        peak_widths = arange(30,101,10)
        peak_ind = signal.find_peaks_cwt(power[ind],widths=peak_widths,max_distances=peak_widths,min_length=2)
        #peak_ind = signal.find_peaks(power[ind],distance=20)
        period_peak_local = period[ind][peak_ind]
        power_peak_local = power[ind][peak_ind]
        # this needs to be improved - I'd like it to take the local maximum, but find_peaks_cwt sometimes give points near
        #    the peak but not at the peak
        # find_peaks may be better but it is not available in the older scipy library on eris

        # save info from this trial to overall dictionaries
        pgram_power_dict[trial_string] = power   # only saving power b/c I'm assuming period list is the same for all!
        peak_dict[trial_string] = {'period': period_peak, 'power': power_peak}
        local_peak_dict[trial_string] = {'period': period_peak_local, 'power': power_peak_local,'peak_ind':peak_ind}
        if i < 0:
            real_peak_power = power_peak
            real_local_peak_power = power_peak_local
        else:
            for p in power_peak_local:
                local_peak_power_list.append(p)
            for p in power_peak:
                peak_power_list.append(p)


# figure to plot local periodogram around Barnes period
#  (iterate to plot all trials - black line is real data, trials are gray lines)
fig = figure(figsize=(18,12))
for i in range(-1,n_trials):
    trial_string = get_trial_string(i)
    power = pgram_power_dict[trial_string]
    power_peak_local = local_peak_dict[trial_string]['power']
    period_peak_local = local_peak_dict[trial_string]['period']
    peak_ind = local_peak_dict[trial_string]['peak_ind']
    
    # plot local periodogram in subplot
    #subplot(nrows,ncols,i+2)
    #title(trial_string, x = 0.05, y = 0.85, horizontalalignment='left')
    if i<0:
        plot_color = 'k'
    else:
        plot_color = '0.6'
    plot(period_sig[ind],power[ind],color=plot_color)
    if i>=0:
        plot(period_sig[ind][peak_ind],power_peak_local,'r*')
    #if i >= (n_trials - ncols):
    xlabel('Period (sigma away from Barnes+17)')
    #if (i+1) % ncols == 0:
    ylabel('Power')
    axis([-3,3,1,9])
plot(period_sig[ind],pgram_power_dict[get_trial_string(-1)][ind],color='k')
legend(('Observed data','Random Start Time Trials'))
savefig('pgram_near_barnes_uvcet.png',bbox_inches='tight')


# figure: histograms of 50 highest peaks from all trials, and peaks close to Barnes period
fig  = figure(figsize=(12,6))
subplot(1,2,1)
hist(peak_power_list,range=(0,10),bins=50,color='0.6',histtype='stepfilled',linewidth=0)
hist(real_peak_power,range=(0,10),bins=50,color='k',histtype='step',linewidth=2)
xlabel('Plavchan Periodogram Power')
ylabel('Number of Peaks')
title('50 Highest Peaks Between 4 and 7 h')

subplot(1,2,2)
hist(local_peak_power_list,range=(0,10),bins=20,color='0.6',histtype='stepfilled',linewidth=0)
hist(real_local_peak_power,range=(0,10),bins=20,color='k',histtype='step',linewidth=2)
xlabel('Plavchan Periodogram Power')
ylabel('Number of Peaks')
title(r'All Peaks within 3 sigma of Barnes')
legend(('Random Start Time Trials','Observed Data'))
savefig('power_histogram.png',bbox_inches='tight')

close('all')
