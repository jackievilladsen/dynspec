# -*- coding: utf-8 -*-
"""

plot.py: Tools for manipulating and plotting dynamic spectra, time series, etc.
This functionality is accessed by the user through the class Dynspec.

"""

from pylab import *
from numpy import *
import matplotlib.pyplot as plt
import os
from astropy.time import Time
from copy import deepcopy

class TimeSec(Time):
    # modify class Time to support using units of MJD in seconds (units of CASA's TIME column)
    def __init__(self,t,format='mjds'):
        if format=='mjds':
            Time.__init__(self,t/24./3600.,format='mjd',scale='utc')
        else:
            Time.__init__(self,t,format=format,scale='utc')
    def mjds(self):
        return self.mjd * 24. * 3600.

def rebin2d(a,wt,binsize):
    shape = tuple(array(a.shape)/array(binsize))
    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
    a1 = a.reshape(sh)
    wt1=wt.reshape(sh)
    tmp,wt2 = average(a1,len(sh)-1,wt1,True)
    return ma.average(tmp,1,wt2)

def rebin1d(a,binsize):
    l = floor(len(a)/binsize)*binsize
    a1 = a[0:l]
    sh = len(a1)/binsize,binsize
    a2 = a1.reshape(sh)
    return ma.average(a2,1)

def rebin1d_ma(a,binsize):
    l = floor(len(a)/binsize)*binsize
    a1 = a[0:l]
    sh = len(a1)/binsize,binsize
    a2 = a1.reshape(sh)
    return ma.average(a2,1)

def rebin2d_ma(b,binsize):
    nt,nf = binsize
    lt,lf = b.shape
    lt_new = floor(lt/nt)*nt
    lf_new = floor(lf/nf)*nf
    a = b[0:lt_new][:,0:lf_new]
    shape = tuple(array(a.shape)/array(binsize))
    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
    a1 = a.reshape(sh)
    tmp,wt2 = ma.average(a1,len(sh)-1,returned=True)
    return ma.average(tmp,1,wt2)

def make_ma(a):
    # make a into a masked array where all zero values are masked
    mask = logical_or.reduce((a==0,isinf(a),isnan(a)))
    return ma.masked_array(a,mask=mask)

def add_band(ma_big,t,f,ma_band,t_band,f_band):
    # add a band to our dynamic spectrum
    # use t_band and f_band to figure out what cells to put it into
    # if there is already data in the big dyn spec, then don't add it
    
    # t_band tells us the indices of the rows in ma_big (i.e., the times) for
    # which this band has data
    t_ind = t_band
    
    # add one col (one freq) at a time since there are some overlapping frequencies
    for i in range(len(f_band)):
        m = ma_band.mask[:,i]
        if not m.all():  # if not all values in this frequency channel are masked
            f_ind = find(f==f_band[i])[0]
            
            mask = ma_big.mask[t_ind,f_ind]
            
            # resulting cells are flagged only if both ma_big and ma_band are flagged
            #    in those cells
            ma_big.mask[t_ind,f_ind] *= ma_band.mask[:,i]
            
            # add data from ma_band if there is no data in the destined cell yet
            # mask = 0 when there is already data in a cell
            ma_big[t_ind,f_ind] += ma_band[:,i].data*mask
    return ma_big

def make_tick_labels(desired_ticks,x):
    # tick_labels,tick_locs = make_tick_labels(desired_ticks,x)
    # desired_ticks is a list of where to put tick marks but can include locations
    # that are not in the range of x (x is the range of values for an axis on an image).
    # This function identifies the values of desired_ticks and returns those as tick_labels,
    # and identifies in the indices in x that are closest to those values and returns those
    # as tick_locs.
    ind = find(logical_and(desired_ticks>=min(x),desired_ticks<=max(x)))
    tick_labels = desired_ticks[ind]
    tick_locs = [find(min(abs(x-t))==abs(x-t))[0] for t in tick_labels]    
    return tick_labels,tick_locs

def closest_ind(x,x0):
    # ind = closest_ind(x,x0)
    # For 1D array x and scalar x0, returns the index of the item in x that has the value closest to x0
    return argmin(abs(x-x0))

def trim_whitespace(dynspec,x=None,y=None):
    # dynspec1,x1,y1 = trim_whitespace(dynspec,x,y)
    # Trim off any rows and colummns on the outer edge of the dynspec that have no
    # unmasked values.
    sum0 = sum(dynspec,1)
    sum1 = sum(dynspec,0)
    if x is None:
        x = arange(len(sum0))
    if y is None:
        y = arange(len(sum1))
    try:
        i = find(sum0.mask==False)
        imin = i[0]
        imax = i[-1]+1
        j = find(sum1.mask==False)
        jmin = j[0]
        jmax = j[-1]+1
        dynspec1 = dynspec[imin:imax,jmin:jmax]
        x1 = x[imin:imax]
        y1 = y[jmin:jmax]
    except:
        print 'no unmasked values'
        dynspec1,x1,y1 = None,None,None
    return dynspec1,x1,y1

def clip_dynspec(dynspec,lims,x=None,y=None,trim_mask=True):
    # spec1,x1,y1 = clip_dynspec(dynspec,lims,x=None,y=None)
    # Given a 2D array dynspec, return a smaller 2D array cut at the indices or x,y values
    # in lims.  lims=[xmin,xmax,ymin,ymax] is assumed to be indices, unless arrays x and y
    # are provided giving the x and y values corresponding to each index, in which case the
    # cut is made at the indices corresponding to the values closest to lims.
    # If trim_mask=True (the default), trim off any rows and columns on the outer edge of the array
    # that have no unmasked values.
    (xlen,ylen) = shape(dynspec)
    [xmin,xmax,ymin,ymax] = lims
    if x is None:
        imin = xmin
        imax = xmax
    else:
        imin = closest_ind(x,xmin)
        imax = closest_ind(x,xmax)
    if y is None:
        jmin = ymin
        jmax = ymax
    else:
        jmin = closest_ind(y,ymin)
        jmax = closest_ind(y,ymax)
    spec1 = dynspec[imin:imax,jmin:jmax]
    x1 = x[imin:imax]
    y1 = y[jmin:jmax]
    if trim_mask:
        return trim_whitespace(spec1,x1,y1)
    return spec1,x1,y1

class Dynspec:
    ''' Dynspec: a class for manipulating and plotting dynamic spectra (using masked arrays) and
        keeping track of the frequencies and time lists corresponding to the array indices.
        
        Here is a list of all object attributes that may be defined by any class routines:
        - self.spec     : dictionary containing entries for each poln product (such as 'rr') - data are masked arrays
        - self.f        : list of frequencies corresponding to dynspec rows
        - self.time     : astropy.time.Time object containing list of times corresponding to dynspec columns        
        '''

    def __init__(self,params={}):
        # initiates a Dynspec object
        self.spec={}
        if 'filename' in params:
            self.load_dynspec(params)
    
    def read_params(self,params):
        # params is a dictionary with certain useful parameters:
        #   params['filename']: directory name to load dynspec from (must contain rr.dat, ll.dat, freq.dat, times.dat)
        #   params['uniform']: regrid to uniform time/frequency sampling after loading dynspec (default False)
        filename = params.get('filename','')
        uniform = params.get('uniform',False)
        return filename,uniform
    
    def load_dynspec(self,params):
        # if filename is a valid file, then loads dynspec (rr,ll,t,f) from that directory
        # self.spec['rr'] and self.spec['ll'] are loaded as masked arrays
        # optional parameter i is used to tell it to load the dynspec starting from time with index i
        # future modification: enable imax as well?
        filename,uniform=self.read_params(params)
        if not os.path.exists(filename):
            print 'Warning: bad dynspec filename:', filename
        else:
            for pol in ['rr','ll','xx','yy','xy','yx']:
                fname = filename + '/' + pol + '.npy'
                if os.path.exists(fname):
                    print 'loading', fname
                    self.spec[pol] = make_ma(load(fname))
                    print pol,'rms:', self.get_rms(pol)*1000, 'mJy'
                        
            self.f=array(loadtxt(filename+'/freq.dat'))       # units: Hz
            
            t=array(loadtxt(filename+'/times.dat'))  # units: MJD in seconds
            self.time = TimeSec(t,format='mjds')     # create Time object containing list of MJD times            
            
            if uniform:
                self.regrid_uniform()                # regrid to uniform time and frequency sampling
            
    
    def make_xlist(self,x,dx,x0=None):
        # xlist: for each element in x, count how many units of dx it is away from x0 (or x[0] if x0 is not defined)
        if x0 is None:
            x0 = x[0]
        diff = (x[1:]-x[:-1])/dx
        diff_int = diff.round().astype(int)
        diff0 = ((x[0]-x0)/dx).round().astype(int)
        xlist = diff0 * ones(shape(x)).astype(int)
        xlist[1:] += cumsum(diff_int)
        return xlist
    
    def make_full_indlist(self,xlist):
        # indlist: return a list counting from 0 to max(xlist)
        return arange(max(xlist)+1)
    
    def get_tlist(self):
        # return tlist: tlist is the amount of time (in units of integration time)
        # that each column in the dynamic spectrum is separated from the first integration (so tlist[0] is 0)
        t = self.time.mjds()
        tlist = self.make_xlist(t,self.dt())
        return tlist
        
    def get_flist(self,df=None):
        # return flist: flist is the number of frequency channels
        # that each row in the dynamic spectrum is separated from the first channel (so flist[0] is 0)
        # df, if given, MUST = self.df()/whole number (so we can add blank channels)
        if df is None:
            df = self.df()
        flist = self.make_xlist(self.f,df)
        return flist

    def gen_x(self,xlist,x0,dx):
        return x0 + xlist * dx
    
    def get_spacing(self,x):
        # return median spacing between elements of x
        # meant to help retrieve integration time or channel width
        return median(x[1:]-x[:-1])
    
    def dt(self):
        # return integration time (duration of dynspec pixels)
        return self.get_spacing(self.time.mjds())
    
    def df(self):
        # return channel width (bandwidth of dynspec pixels)
        return self.get_spacing(self.f)
    
    def set_time(self,tlist,t0):
        # given t0 and tlist (tlist is in units of integration times, t0 in MJD seconds), set self.time
        # as a Time object with a correct list of times in MJD seconds
        t = self.gen_x(tlist,t0,self.dt())
        self.time = TimeSec(t,format='mjds')
        
    def set_freq(self,flist,f0):
        # given f0 and flist (flist is in units of self.df(), f0 in Hz), set self.f as an array of frequencies
        # in units of Hz
        self.f = self.gen_x(flist,f0,self.df())
    
    def regrid_uniform(self,df=None):
        # regrid by adding blank rows and columns so that time and frequency sampling is uniform
        
        # make list of times that counts from 0 and includes times when there are no data
        tlist = self.get_tlist()
        t = self.make_full_indlist(tlist)
        tlen = len(t)
        
        # make list of frequencies that has no gaps
        flist = self.get_flist(df) 
        f = self.make_full_indlist(flist)
        flen = len(f)
        
        # cycle through all poln products in self.spec
        for pol in self.spec.keys():
            # create empty dynspec with regridded dimensions
            spec = ma.zeros((tlen,flen)) * 0j
            spec.mask = ones((tlen,flen))
        
            # use add_band to regrid onto new dimensions, overwriting self.spec[pol]
            print 'regridding', pol, 'onto uniform time-frequency grid'
            self.spec[pol] = add_band(spec,t,f,self.spec[pol],tlist,flist)
            spec = None    # so that it stops taking up memory
        
        # overwrite self.time and self.f with new values
        t0 = min(self.time.mjds())
        f0 = min(self.f)
        self.set_time(t,t0)
        self.set_freq(f,f0)
    
    def add_dynspec(self,dyn):
        # merge another Dynspec object dyn with this object
        # regrid etc as necessary to make them both fit in the dynamic spectrum
        # Current approach: dt must be the same, df must be integer multiples;
        #  new dynspec will have alternating blank channels in portion from more coarsely gridded dynspec
        
        dt = self.dt()
        t1 = self.time.mjds()
        t2 = dyn.time.mjds()
        t0 = min(concatenate([t1,t2]))
        tlist1 = self.make_xlist(t1,dt,x0=t0)
        tlist2 = self.make_xlist(t2,dt,x0=t0)
        tlen = max(concatenate([tlist1,tlist2]))+1
        tlist = arange(tlen)
        
        df = min(self.df(),dyn.df())
        # round self.f and dyn.f to nearest integer multiple of df (to ensure they are on same frequency grid)
        d_selff = df/2 - (self.f[0] + df/2) % df
        d_dynf = df/2 - (dyn.f[0] + df/2) % df
        selff = self.f + d_selff
        dynf = dyn.f + d_dynf
        if d_selff != 0 or d_dynf != 0:
            print 'add_dynspec is adding', d_selff/1e6, 'MHz to self.f and', d_dynf/1e6, 'MHz to dyn.f to make even frequency grid'
        ftemp = concatenate([selff,dynf])
        fmin = min(ftemp)
        fmax = max(ftemp)
        f = arange(fmin,fmax+df,df)
        flen = len(f)
        
        # cycle through all poln products in self.spec
        for pol in self.spec.keys():
            # create big empty masked array (with dimensions big enough to hold both dynspec)
            spec = ma.zeros((tlen,flen)) * 0j
            spec.mask = ones((tlen,flen))
            print 'merging dynspec:', pol
            # add our own dynspec to big dynspec
            spec = add_band(spec,tlist,f,self.spec[pol],tlist1,selff)
            # add new dynspec to big dynspec
            spec = add_band(spec,tlist,f,dyn.spec[pol],tlist2,dynf)
            # overwrite self.spec[pol] with new dynspec
            self.spec[pol] = spec
            spec = None
            
        # redefine self.time and self.f
        self.set_time(tlist,t0)
        self.f = f
        
    def t0(self):
        # return string format for min time in self.time
        return min(self.time).iso

    def fGHz(self):
        # return freq list in GHz
        return self.f/1.e9
    
    def get_rms(self,pol='',func=imag):
        # return rms of dynspec in Jy (calculates RMS per channel then takes median RMS)
        # default is RMS of Im(LL), but can use any pol that is a key in self.spec,
        # and func can be any function converting complex numbers to real numbers (such as imag, real, abs, angle)
        if pol == '':
            pol = self.spec.keys()[0]
        rms = ma.median(ma.std(func(self.spec[pol]),0)).data
        try:
            return rms[0] # in case rms is a single-valued array - this happens sometimes but not always, not sure why
        except:
            return rms
        
    def mask_RFI(self,rmsfac=5.):
        print 'Warning: mask_RFI handles the masking wrong, need to fix this!'
        print 'masking chans w/ rms >',rmsfac,'* median rms'
        for pol in self.spec.keys():
            rms = std(imag(self.spec[pol]),0)  # rms vs freq
            medrms = ma.median(rms)
            chanmask = (rms > rmsfac*medrms)
            n = sum(chanmask)
            mask0 = self.spec[pol].mask
            self.spec[pol].mask = ~ (~mask0 * ~chanmask)
            self.spec[pol] = self.spec[pol] * (1-chanmask)
            print n, 'channels masked for', pol, '- new RMS:', self.get_rms(pol)*1000, 'mJy'

    def mask_RFI_pixels(self,rmsfac=5.):
        print 'masking dynspec pixels >',rmsfac,'* median rms'
        for pol in self.spec.keys():
            rms = std(imag(self.spec[pol]),0)  # rms vs freq
            medrms = ma.median(rms)
            mask = abs(self.spec[pol]) > rmsfac*medrms
            n = sum(mask)
            self.spec[pol].mask = ma.mask_or(self.spec[pol].mask,mask)
            self.spec[pol] = self.spec[pol] * (1-mask)
            print n, 'pixels masked for', pol, '- new RMS:', self.get_rms(pol)*1000, 'mJy'

    def mask_SNR(self,rmsfac=3.,func=real):
        # returns dynspec with pixels masked that have SNR < rmsfac (where SNR is calculated relative to channel RMS in Imag component)
        ds = deepcopy(self)
        for pol in self.spec.keys():
            rms_spec = self.rms_spec(pol=pol)
            if pol in ['v','rc']:
                lowSNR = abs(func(self.spec[pol])) < (rmsfac * rms_spec)  # creates masked array of True/False values showing where to mask
            else:
                lowSNR = func(self.spec[pol]) < (rmsfac * rms_spec)  # creates masked array of True/False values showing where to mask
            ds.spec[pol] = ma.masked_where(lowSNR,self.spec[pol])
        return ds
        
    def bin_dynspec(self,nt,nf,mask_partial=1.):
        # returns a new dynspec object with binning in time or frequency
        # nt and nf are the number of time and frequency channels to bin together
        # mask_partial: if <1, mask pixels in new dynspec where a fraction of the contributing
        #               pixels in the original dynspec > mask_partial (e.g. 0.5) were masked
        # A lower value for mask_partial will flag more data
        
        # create empty Dynspec object for new binned dynamic spectrum
        ds = Dynspec()
        
        # calculate pixel duration and bandwidth for binned dynspec
        dt = nt * self.dt()
        df = nf * self.df()
        print 'binning dynamic spectrum to resolution of', dt, 'sec and', df/1.e6, 'MHz'
        
        # cycle through all poln products in self.spec
        for pol in self.spec.keys():
            # bin dynamic spectrum
            ds.spec[pol] = rebin2d_ma(self.spec[pol],(nt,nf))
            if mask_partial < 1.:
                frac_masked = rebin2d_ma(self.spec[pol].mask,(nt,nf))
                frac_too_high = frac_masked > mask_partial
                new_mask = ma.mask_or(ds.spec[pol].mask,frac_too_high)
                ds.spec[pol].mask = new_mask
        
        # create new time list with center times for each bin
        t = self.time.mjds()
        t_bin = rebin1d(t,nt)
        ds.time = TimeSec(t_bin,format='mjds')
        
        # create new frequency list with center times for each bin
        ds.f = rebin1d(self.f,nf)
        
        # print rms of new dynspec
        print 'binned dynspec rms:', ds.get_rms()*1000, 'mJy'
        
        return ds


    def plot_dynspec(self,plot_params={}):
        # create an imshow color plot of the dynamic spectrum
        
        func = plot_params.get('func',real) # part of complex visibility to plot (real, imag, abs, angle)
        pol = plot_params.get('pol','rr')
        rmspol = plot_params.get('rmspol',pol)
        
        # generate automatic plot limits (units: flux in Jy)
        smin = self.get_rms(rmspol)*3 # default minimum flux on color scale is 3*RMS (median channel-based RMS)
        smax = percentile(func(self.spec[pol]),99) # default max flux on color scale is (99th percentile flux)

        # plot params #
        scale = plot_params.get('scale','linear')    # options: log, linear       
        smin = plot_params.get('smin',smin)
        smax = plot_params.get('smax',smax)
        dx = plot_params.get('dx',0.)             # spacing between x axis tick marks - time in minutes (default: 0 --> auto)
        dy = plot_params.get('dy',0.)             # spacing between y axis tick marks - frequency in GHz (default: 0 --> auto)
        tlims = self.time[0].mjds()+plot_params.get('tlims',array([0,1e6]))*60.             # min and max time to plot (in min since beginning of obs)
        flims = plot_params.get('flims',array([min(self.f),max(self.f)+1]))                 # min and max frequencies to plot (in Hz)
        ar0 = plot_params.get('ar0',1.0)
        axis_labels = plot_params.get('axis_labels',['xlabel','ylabel','cbar','cbar_label'])
        trim_mask = plot_params.get('trim_mask',True) # whether to cut off fully masked edges when making plot
        
        # clip dynspec to match tlims, flims
        spec,t,f = clip_dynspec(func(self.spec[pol]),[tlims[0],tlims[1],flims[0],flims[1]],self.time.mjds(),self.f,trim_mask=trim_mask)
        t0 = TimeSec(t[0],format='mjds')
        
        ar = ar0*len(t)/len(f)
        
        ## Large plot (entire dynspec, binned) ##
        gca().set_axis_bgcolor('k')
        print func.func_name,pol,smin,smax
        if scale=='log':
            imshow(log10(spec).T,aspect=ar,vmin=log10(smin),vmax=log10(smax),origin='lower',cmap='seismic')
            ds = round(log10(smax)-log10(smin),1)/5         # spacing between colorbar ticks
            ticks = arange(log10(smin),log10(smax)+ds,ds)   # colorbar tick locations
            ticklbls = (10**(ticks+3)).round().astype(int)  # colorbar tick labels
        else: # scale = 'linear'
            imshow(spec.T,aspect=ar,vmin=smin,vmax=smax,origin='lower',cmap='seismic')
            ds = round((smax-smin),2)/5                            # spacing between colorbar ticks
            if ds == 0.0:
                ds = (smax-smin)/10.
            tickmin = sign(smin) * (abs(smin) - (abs(smin) % ds))
            tickmax = smax - (smax % ds)
            ticks = arange(tickmin,tickmax+ds,ds)                    # colorbar tick locations
            ticklbls = np.round(ticks*1000,1)                           # colorbar tick labels
        
        # add colorbar and change labels to show scale
        if 'cbar' in axis_labels:
            cbar = colorbar(fraction=0.046, pad=0.04) # fraction=0.046,pad=0.04 - these numbers magically make colorbar same size as plot
            if pol is 'rc':
                if 'cbar_label' in axis_labels:
                    cbar.set_label('Percent Circular Polarization')
                ticks = arange(-1.,1.1,0.2)
                ticklbls = arange(-100,101,20)
                cbar.set_ticks(ticks)
                cbar.set_ticklabels(ticklbls)
            else:
                if 'cbar_label' in axis_labels:
                    cbar.set_label('Flux (mJy)')
                cbar.set_ticks(ticks)
                cbar.set_ticklabels(ticklbls)

        # label x axis in minutes
        t_minutes = (t-t0.mjds())/60.  # get time since beginning of observation in minutes
        xmax = max(t_minutes)
        if dx == 0.:
            dx = ceil(xmax/40)*10
        xtick_lbls = arange(0,xmax+dx,dx)
        tick_labels,tick_locs = make_tick_labels(xtick_lbls,t_minutes)
        xticks(tick_locs, tick_labels)
        if 'xlabel' in axis_labels:
            xlabel('Time in minutes since '+ t0.iso +' UTC')

        # label y axis in GHz
        f_GHz = f/1.e9
        if dy == 0.:
            dy = ceil((max(f_GHz)-min(f_GHz))/4*2)/2
        ymin = round(min(f_GHz)/dy)*dy
        ymax = round(max(f_GHz)/dy)*dy
        ytick_lbls = arange(ymin,ymax+dy,dy)
        tick_labels,tick_locs = make_tick_labels(ytick_lbls,f_GHz)
        yticks(tick_locs, tick_labels)
        if 'ylabel' in axis_labels:
            ylabel('Frequency (GHz)')

    def clip(self,tmin=0,tmax=1e6,fmin=0,fmax=1.e12,trim_mask=True):
        # returns a new Dynspec object clipped to the time and frequency
        #  limits specified (tmin and tmax are time in minutes since beginning of obs,
        #  fmin and fmax are in Hz)
        
        ds = Dynspec()  # create an empty Dynspec object  - need to define ds.spec[pol], ds.time, ds.f
        
        # calculate time in minutes since beginning of obs (since this is what we're using to clip)
        mjds0 = self.time.mjds()[0]
        t_minutes = (self.time.mjds()-mjds0)/60.
        
        # clip each pol's spec
        for pol in self.spec.keys():
            ds.spec[pol],t,ds.f = clip_dynspec(self.spec[pol],[tmin,tmax,fmin,fmax],t_minutes,self.f,trim_mask=trim_mask)

        # calculate mjds times of new dynspec
        t_mjds = t*60. + mjds0
        ds.time = TimeSec(t_mjds,format='mjds')
        
        return ds

    def rms_spec(self,pol='i',func=imag):
        # return the RMS spectrum (RMS in each channel) in the complex func of the specified pol
        #  default is RMS of imag(I)
        return std(func(self.spec[pol]),0)

    def tseries(self,fmin=0,fmax=1.e12,weight_mode='rms',trim_mask=False):
        # return a Dynspec object that is a time series integrated from fmin to fmax
        # weight_mode: 'rms' --> weight by 1/rms^2; anything else --> no weights
        ds = self.clip(fmin=fmin,fmax=fmax,trim_mask=trim_mask)
        tseries = Dynspec()
        tseries.time = ds.time
        tseries.f = mean(ds.f)
        if weight_mode == 'rms':
            rms = ds.rms_spec()
            wt = 1/rms**2
        else:
            wt = ones(len(ds.f))
        for pol in ds.spec.keys():
            tseries.spec[pol] = ma.average(ds.spec[pol],1,wt)
        return tseries
