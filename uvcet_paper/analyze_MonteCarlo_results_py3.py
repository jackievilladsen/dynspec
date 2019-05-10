'''
analyze_MonteCarlo_results.py:

Load the pgram and peak file outputs from the NASA Exoplanet Archive Periodogram Service,
for the Monte Carlo trials with different randomized observation times, then analyze and plot
the statistics of periodogram peaks to determine the significance of the peak in the periodogram
of the real data.

# alternative to execfile in python 3:
filename = '/users/jvillads/casa_utils/dynspec/uvcet_paper/analyze_MonteCarlo_results_py3.py'
exec(open(filename).read())

'''

##### USER DEFINED PARAMETERS #####

# set reload=False if I have already run this script once in ipython and don't want to reload the text files (saves time as I iterate in plotting/data analysis)
reload=False

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

if reload:
    
    # containers to compile pgram data for all trials once loaded
    pgram_power_dict = {}
    peak_dict = {}
    
    for i in range(-1,n_trials): # i = -1 is for the real data

        trial_string = get_trial_string(i)
        if i < 0: # i = -1 is reserved for real data
            pgram_file = 'realdata_pgram.tbl'
            peak_file = 'realdata.dat.top'
        else:    # i >=0 are Monte-Carlo trials
            pgram_file = 'trial'+str(i)+'_pgram.tbl'
            peak_file = 'trial'+str(i)+'.dat.top'

        # load periodogram plot and list of top 50 peaks
        if not os.path.exists(pgram_file) or not os.path.exists(peak_file):
            print('WARNING: One or more of these pgram files for', trial_string, 'does not exist, skipping this trial:',pgram_file+', '+peak_file)
            continue
        print('Loading data for '+trial_string+'...')
        period,power = loadtxt(pgram_file,skiprows=3,usecols=(1,2),unpack=True)
        period_peak,power_peak = loadtxt(peak_file,skiprows=22,usecols=(1,2),unpack=True)
        
       # save info from this trial to overall dictionaries
        pgram_power_dict[trial_string] = power   # only saving power b/c I'm assuming period list is the same for all!
        peak_dict[trial_string] = {'period': period_peak, 'power': power_peak}
        if i < 0:
            real_peak_power = power_peak

# calculate indices and period-values w/in Nsig_barnes sigma of Barnes periods
period_sig = (period - period_uvcet_barnes) / sigP_uvcet_barnes
ind = find( abs(period_sig) < Nsig_barnes )
period_sig_blcet = (period - period_blcet_barnes) / sigP_blcet_barnes
ind_blcet = find( abs(period_sig_blcet) < Nsig_barnes )

# extract just highest power peak from all trials
peak1_list_trials_only = array([amax(peak_dict[t]['power']) for t in tstring_list if t[0]=='T'])
local_maxpower_list_trials_only = array([amax(pgram_power_dict[t][ind]) for t in tstring_list if t[0]=='T'])
local_blcet_maxpower_list_trials_only = array([amax(pgram_power_dict[t][ind_blcet]) for t in tstring_list if t[0]=='T'])

# figure to plot local periodogram around Barnes period
#  (iterate to plot all trials - black line is real data, trials are gray lines)
fig1 = figure(figsize=(8,6))
fig2 = figure(figsize=(8,6))
for i in range(-1,n_trials):
    trial_string = get_trial_string(i)
    power = pgram_power_dict[trial_string]
    
    # plot local periodogram in subplot
    if i<0:
        plot_color = 'k'
        label_line = 'Observed Data'
        label_star = None
    else:
        plot_color = '0.6'
        label_star = None
        if i==0:
            label_line='Random Trials'
            label_star = 'Peak Power'
        else:
            label_line=None
    figure(fig1.number)
    plot(period_sig[ind],power[ind],color=plot_color,label=label_line)
    i_max = argmax(power[ind])
    plot(period_sig[ind][i_max],power[ind][i_max],'r*',label=label_star)
    xlabel(r'Period ($\sigma$ away from B17 UV Cet)')
    ylabel('Power')
    axis([-3,3,1,9])

    figure(fig2.number)
    plot(period_sig_blcet[ind_blcet],power[ind_blcet],color=plot_color,label=label_line)
    i_max = argmax(power[ind_blcet])
    plot(period_sig_blcet[ind_blcet][i_max],power[ind_blcet][i_max],'r*',label=label_star)
    xlabel(r'Period ($\sigma$ away from B17 BL Cet)')
    ylabel('Power')
    axis([-3,3,1,9])

figure(fig1.number)
plot(period_sig[ind],pgram_power_dict[get_trial_string(-1)][ind],color='k',label=None)
legend()
savefig('pgram_near_barnes_uvcet.png',bbox_inches='tight')

figure(fig2.number)
plot(period_sig_blcet[ind_blcet],pgram_power_dict[get_trial_string(-1)][ind_blcet],color='k',label=None)
legend()
savefig('pgram_near_barnes_blcet.png',bbox_inches='tight')

# figure: histograms of highest peak from all trials (both for full 4 to 7 hour range, and for within 3 sigma of Barnes period)
fig  = figure(figsize=(5,4))
hist(peak1_list_trials_only,range=(0,10),bins=1000,color='0.6',histtype='stepfilled',linewidth=0, cumulative=-1, density=True)
#hist(local_blcet_maxpower_list_trials_only,range=(0,10),bins=1000,color='0.8',histtype='stepfilled',linewidth=2, cumulative=-1, density=True)
hist(local_maxpower_list_trials_only,range=(0,10),bins=1000,color='0.3',histtype='stepfilled',linewidth=2, cumulative=-1, density=True)
gca().vlines(amax(real_peak_power),0,2,color='k',linewidth=2)
xlabel('Plavchan Periodogram Highest Peak Power')
ylabel('Fraction of Trials Above This Power')
#legend(('Random Trials - 4 to 7 h',r'Random Trials - B17 BL Cet $\pm$ 3$\sigma$',r'Random Trials - B17 UV Cet $\pm$ 3$\sigma$','Observed Data'),framealpha=0.0,loc='upper center',bbox_to_anchor=(0.5,1.25))
legend(('Random Trials - 4 to 7 h',r'Random Trials - B17 UV Cet $\pm$ 3$\sigma$','Observed Data'),framealpha=0.0,loc='upper center',bbox_to_anchor=(0.5,1.25))
axis([0,10,0,1.0])

savefig('power_histogram.png',bbox_inches='tight')

close('all')
