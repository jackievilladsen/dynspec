# -*- coding: utf-8 -*-
"""
gen_random_uvcet_tseries.py: Like make_uvcet_tseries.py, this script creates a time series for the
  UV Ceti observations, as a text file usable by the NASA Periodogram Service: one column
  is time in TDB (Barycentric Dynamical Time), other column is flux density.

  However, this script randomizes the start time of each observation within the 4-month period that
  the observations occurred.  By generating many such time series and running the Plavchan periodogram
  on them, I can test whether the UV Cet period is significant or if it is just an effect of the sampling
  or small number of observations.
"""

from pylab import *
from numpy import *

# number of synthetic time series to generate
N_trials = 100

mydir = '/data/jrv/15A-416/UVCet/'
band = 'L'

# load original multi-obs tseries from file - this is much faster than re-generating from individual dynspecs
tseries_filename=mydir+'tseries/Lband_rr_doclipFalse_tseries.dat'
[tlist_original,fluxlist]= loadtxt(tseries_filename,unpack=True)

# get indices positions in tlist_original that correspond to starts of observations
dt_orig = tlist_original[1:]-tlist_original[:-1]
ind = zeros(5,dtype='int')  # need a 0 as the start index of the first observation
ind[1:5] = find(abs(dt_orig)>5.)+1
# observations are all 4 hours or less, so a gap of 5 hours or more is a break between obs
# the observations are not in order so need to take absolute value of time difference (but times are in order
#    within a single observation)

# NEXT STEP: USE START AND TIMES OF FIRST AND LAST OBS AS MAX AND MIN START TIME?  OR BEGINNING AND END OF THE
#   MONTHS THEY WERE IN?
tstart_orig = tlist_original[ind]

# figure out what range of start times to pick from
# Goal: median difference between earliest and latest start time is comparable to difference for real observations
# How much time range do we need to pick from? Plan: generate many trials of 5 numbers from zero to 1, and then see
#    what the median range is. Then multiple our real t5-t1 time 1/median to get the time range to generate times from.
'''
rand_5_nums = rand(10000000,5)
range_5 = amax(rand_5_nums,1)-amin(rand_5_nums,1)
time_range_factor = 1/median(range_5)
print 'Time range expanded by a factor of', time_range_factor
# typical value is 1.457
'''
time_range_factor = 1.457
tstart_min = amin(tstart_orig)
orig_tstart_diff = amax(tstart_orig)-tstart_min
tstart_range = orig_tstart_diff*time_range_factor
print 'Choosing 5 random start times within a period of roughly', tstart_range/24/30.44, 'months.'
print 'Conducting',N_trials, 'trials.'

# iterate once for each trial
tstart_list = zeros((N_trials,5))  # each row will contain the start times from one trial, so I can plot them
for i in range(N_trials):
    # generate a list of random start times for the trial,
    # with all separated by at least 6 hours
    min_tsep = 0.0
    while min_tsep < 6.0:
        tstart = rand(5) * tstart_range + tstart_min
        min_tsep = amax(tstart)-amin(tstart)
    tstart_list[i,:] = tstart # save start times to an array so I can plot them easily
    
    # for each observation, modify tlist to match that observation's start time
    tlist = tlist_original.copy()
    for j in range(5):  # cycle through the 5 observations
        # identify index range that corresponds to this observation in tlist
        ind_start = ind[j]
        if j==4:
            ind_stop = len(tlist)
        else:
            ind_stop = ind[j+1]
        
        # calculate difference between original start time and desired start time,
        # then add that to all tlist values for this observation so they are all
        # set relative to the desired start time
        tdiff = tstart[j]-tstart_orig[j]
        tlist[ind_start:ind_stop] += tdiff
    
    # save modified tlist and flux density list to file
    filename = mydir+'tseries/MonteCarlo/trial'+str(i)+'.dat'
    savetxt(filename,column_stack((tlist,fluxlist)))

# plot synthetic observation start times    
figure()
color_list = ['k','b','g','c','m']
for i in range(5):
    plot(tstart_list[:,i],'.',color=color_list[i])
ax = gca()
ax.hlines(tstart_orig,-1,N_trials,colors=color_list)
figfile = mydir+'tseries/MonteCarlo/start_times.png'
savefig(figfile,bbox_inches='tight')

# save list of synthetic start times to file for reference
filename = mydir+'tseries/MonteCarlo/start_times_list.dat'
savetxt(filename,tstart_list)

# calculate median range of synthetic start times
tstart_minmax_diff = amax(tstart_list,1)-amin(tstart_list,1)
median_tstart_diff = median(tstart_minmax_diff)
print 'Start time difference from first to last real observation: ~',orig_tstart_diff/24./30.44,'months'
print 'Median first/last start time difference of synthetic observations: ~',median_tstart_diff/24./30.44,'months'
