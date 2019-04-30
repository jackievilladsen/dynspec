# -*- coding: utf-8 -*-
"""
make_uvcet_tseries.py: Purpose of this script is to export a single time series concatenating the five 2015
  UV Ceti observations, as a text file usable by the NASA Periodogram Service: one column
  is time in TDB (Barycentric Dynamical Time), other column is flux density.

  Will write out a few different files, for combos of [I,V,RR,LL] with [L band, S band].
"""

from pylab import *
from numpy import *
import os
from dynspec.plot import *
from astropy import time, coordinates, units
from dynspec.pipeline_utils import load_band_filelist

# polarization product to use
pol = 'rr'

# if do_clip is True, will replace flux density values GREATER than clip_thresh_mJy with clip_thresh_mJy
# The purpose is to see if this helps with the box fitting periodogram algorithm since the bursts are variable
do_clip = True
clip_thresh_mJy = 15

# UV Cet coords in July 4, 2015 VLA observation: 1h39m5.103 -17d56m51.87
# (quite good agreement with SIMBAD coords for that date, probably because the SIMBAD
#   coords come from Gaia so the binary position is from a similar epoch to the observations)
# Should be fine to use same coords for all epochs - parallax should not have a huge effect on travel time
#  (compared to the time uncertainty introduced by variability in the bursts)
RA_hms_UVCet = '01:39:05.103'
dec_dms_UVCet = '-17:56:51.87'

def convert_to_bary_time(t_obs,RA_hms='00:00:00.0',dec_dms='00:00:00.0',site_name='vla'):
    '''
    convert_to_bary_time(t_obs,RA_hms='00:00:00.0',dec_dms='00:00:00.0',site_name='vla'):
    
    Convert from time at the observatory (t_obs - a Time object) to solar system
    barycentric time in TBD format by adding the light-travel time difference, which
    depends upon the coordinates of the target (RA_hms,dec_dms).

    Must run individually for each element in ds.time :/
    
    Run astropy.coordinates.EarthLocation.get_site_names() to see available site names.
    '''
    observatory_coords = coordinates.EarthLocation.of_site(site_name)
    target_coords = coordinates.SkyCoord(RA_hms,dec_dms,unit=(units.hourangle,units.deg),frame='icrs')
    dt_barycenter_observatory = t_obs.light_travel_time(target_coords,location=observatory_coords)
    t_barycentric = t_obs.tdb + dt_barycenter_observatory  # in TDB time
    return t_barycentric

close('all')

mydir = '/data/jrv/15A-416/UVCet/'
bands = ['S','L']

filelist=[]
times = {'L':[],'S':[]}
flux = {'L':[],'S':[]}
i=0

# maybe restructure to pull filelist from best_dynspec_file list
for band in bands:
    all_ds_filelist = load_band_filelist(band)
    filelist = [f for f in all_ds_filelist if 'UVCet' in f and '15A-416' in f]
    
    leg = []
    i+=1
    subplot(2,1,i)
    
    for dsfile in filelist:
        
        #filedir = filename[0:28] # where tseries and times list will be saved - is this really necessary?        
        
        # Load dsfile
        print 'Loading dynamic spectrum from', dsfile
        if pol in ['rr','ll']:
            convert_stokes=False
        else:
            convert_stokes=True
        params={'filename':dsfile,'uniform':True,'convert_stokes':convert_stokes}
        if band=='S':
            params['df']=2.e6
        ds_big = Dynspec(params)
        # decide whether to mask RFI here
        
        # bin ds in time before making tseries (can experiment with binning options)
        ds = ds_big.bin_dynspec(30,2) # rebin dynspec to 30 sec, 2 MHz resolution
        
        # make time series
        tseries = ds.tseries(weight_mode='rms',clipds=False)
        
        # clip values > clip_thresh_mJy and replace them with clip_thresh_mJy
        if do_clip:
            ind=find(real(tseries.spec[pol])>clip_thresh_mJy * 0.001)
            tseries.spec[pol][ind] = clip_thresh_mJy * 0.001 * (1+0j) * ones(len(ind))
        
        # convert time list from UTC to TDB
        t_TDB = [convert_to_bary_time(t,RA_hms=RA_hms_UVCet,dec_dms=dec_dms_UVCet,site_name='vla') for t in tseries.time]
        t_TDB_hours = array([t.value*24. for t in t_TDB])
            
        # select only the times that are not masked
        ind = find(~tseries.spec[pol].mask)
        
        # add those points to the times and flux lists for all obs together (for this band)
        times[band] += t_TDB_hours[ind].tolist()
        flux[band] += real(tseries.spec[pol])[ind].tolist()
        
        # plot this obs/band's tseries so I can see what they all look like
        # (time axis is relative to obs start time so they won't be phased)
        plot(t_TDB_hours[ind]-t_TDB_hours[0],real(tseries.spec[pol])[ind])
        obs = dsfile.split('/')[5]
        leg.append(obs+band)
    
    # label the plot of all time series from this band
    legend(leg)
    
    # save a time series of all observations from this band together
    # ! look at this file ! this may already be what I need.
    tseries_filename=mydir+'tseries/'+band+'band_'+pol+'_doclip'+str(do_clip)+'_tseries.dat'
    savetxt(tseries_filename,column_stack((times[band],flux[band])))

# save a figure of the time series in both bands (one subplot per band)
savefig(mydir+'tseries/tseries_band_'+pol+'_doclip'+str(do_clip)+'.pdf',bbox_inches='tight')


# plot the time series on an absolute time axis
#    --> shows the spacing between obs, can't see variation during individual epochs
figure()
for band in bands:
    plot(times[band],flux[band])
legend(bands)
savefig(mydir+'tseries/long_'+pol+'_doclip'+str(do_clip)+'.pdf',bbox_inches='tight')

