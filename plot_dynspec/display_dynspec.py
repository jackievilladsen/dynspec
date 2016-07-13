# -*- coding: utf-8 -*-
"""
Last modified: August 7 2015

@author: jackie
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
            Time.__init__(self,t/24./3600.,format='mjd')
        else:
            Time.__init__(self,t,format=format)    
    def mjds(self):
        return self.mjd * 24. * 3600.

def rebin2d(a,wt,binsize):
    shape = tuple(array(a.shape)/array(binsize))
    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
    a1 = a.reshape(sh)
    wt1=wt.reshape(sh)
    tmp,wt2 = average(a1,len(sh)-1,wt1,True)
    return average(tmp,1,wt2)

def rebin1d(a,binsize):
    l = floor(len(a)/binsize)*binsize
    a1 = a[0:l]
    sh = len(a1)/binsize,binsize
    a2 = a1.reshape(sh)
    return average(a2,1)

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
    return ma.masked_array(a,mask=(a==0))

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
        #   params['i']: load dynspec only starting from time with that index (default 0)
        #   params['uniform']: regrid to uniform time/frequency sampling after loading dynspec (default False)
        filename = params.get('filename','')
        i = params.get('i',0)
        uniform = params.get('uniform',False)
        return filename,i,uniform
    
    def load_dynspec(self,params):
        # if filename is a valid file, then loads dynspec (rr,ll,t,f) from that directory
        # self.spec['rr'] and self.spec['ll'] are loaded as masked arrays
        # optional parameter i is used to tell it to load the dynspec starting from time with index i
        # future modification: enable imax as well?
        filename,i,uniform=self.read_params(params)
        if os.path.exists(filename):
            print 'loading RR from', filename
            self.spec['rr'] = make_ma(loadtxt(filename+'/rr.dat'))[i:,:]
            print 'RCP rms:', self.get_rms('rr')*1000, 'mJy'
            
            print 'loading LL from', filename
            self.spec['ll'] = make_ma(loadtxt(filename+'/ll.dat'))[i:,:]
            print 'LCP rms:', self.get_rms('ll')*1000, 'mJy'
            
            self.f=array(loadtxt(filename+'/freq.dat'))       # units: Hz
            
            t=array(loadtxt(filename+'/times.dat')[i:])  # units: MJD in seconds
            self.time = TimeSec(t,format='mjds')         # create Time object containing list of MJD times            
            
            if uniform:
                self.regrid_uniform()                    # regrid to uniform time and frequency sampling
        else:
            print 'Warning: bad dynspec filename:', filename
    
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
        
    def get_flist(self):
        # return flist: flist is the number of frequency channels
        # that each row in the dynamic spectrum is separated from the first channel (so flist[0] is 0)
        flist = self.make_xlist(self.f,self.df())
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
    
    def regrid_uniform(self):
        # regrid by adding blank rows and columns so that time and frequency sampling is uniform
        
        # make list of times that counts from 0 and includes times when there are no data
        tlist = self.get_tlist()
        t = self.make_full_indlist(tlist)
        tlen = len(t)
        
        # make list of frequencies that has no gaps
        flist = self.get_flist()
        f = self.make_full_indlist(flist)
        flen = len(f)
        
        ll = ma.zeros((tlen,flen))
        ll.mask = ones((tlen,flen))
        
        # cycle through all poln products in self.spec
        for pol in self.spec.keys():
            # create empty dynspec with regridded dimensions
            spec = ma.zeros((tlen,flen))
            spec.mask = ones((tlen,flen))
        
            # use add_band to regrid onto new dimensions, overwriting self.spec[pol]
            print 'regridding', pol
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
        # assume for now that they both have the same df and dt
        # deal with fixing that later
        
        dt = self.dt()
        t1 = self.time.mjds()
        t2 = dyn.time.mjds()
        t0 = min(concatenate([t1,t2]))
        tlist1 = self.make_xlist(t1,dt,x0=t0)
        tlist2 = self.make_xlist(t2,dt,x0=t0)
        tlen = max(concatenate([tlist1,tlist2]))+1
        tlist = arange(tlen)
        
        df = self.df()
        ftemp = concatenate([self.f,dyn.f])
        fmin = min(ftemp)
        fmax = max(ftemp)
        f = arange(fmin,fmax+df,df)
        flen = len(f)
        
        # cycle through all poln products in self.spec
        for pol in self.spec.keys():
            # create big empty masked array (with dimensions big enough to hold both dynspec)
            spec = ma.zeros((tlen,flen))
            spec.mask = ones((tlen,flen))
            print 'merging dynspec:', pol
            # add our own dynspec to big dynspec
            spec = add_band(spec,tlist,f,self.spec[pol],tlist1,self.f)
            # add new dynspec to big dynspec
            spec = add_band(spec,tlist,f,dyn.spec[pol],tlist2,dyn.f)    
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
    
    def get_rms(self,pol='ll'):
        # return rms of dynspec in Jy
        # default pol product is 'll', but can use any pol that is a key in self.spec
        rms = ma.median(ma.std(self.spec[pol],0)).data
        try:
            return rms[0]
        except:
            return rms
        
    def bin_dynspec(self,nt,nf):
        # returns a new dynspec object with binning in time or frequency
        # nt and nf are the number of time and frequency channels to bin together
        
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
        
        # create new time list with center times for each bin
        t = self.time.mjds()
        t_bin = rebin1d(t,nt)
        ds.time = TimeSec(t_bin,format='mjds')
        
        # create new frequency list with center times for each bin
        ds.f = rebin1d(self.f,nf)
        
        # print rms of new dynspec
        print 'binned dynspec LCP rms:', ds.get_rms()*1000, 'mJy'
        
        return ds


    def plot_dynspec(self,plot_params={}):
        # create an imshow color plot of the dynamic spectrum
        # no binning
        
        pol = plot_params.get('pol','rr')
        rmspol = plot_params.get('rmspol','ll')
        # generate automatic plot limits (units: flux in Jy)
        smin = self.get_rms(rmspol)*1.5
        smax = percentile(self.spec[pol],98)*1.3 # get the 97th percentile flux

        # plot params #                                                        
        scale = plot_params.get('scale','log')    # options: log, linear       
        smin = plot_params.get('smin',smin)
        smax = plot_params.get('smax',smax)
        dx = plot_params.get('dx',5.)             # spacing between x axis tick marks - time in minutes (default: 5 min)
        dy = plot_params.get('dy',0.5)            # spacing between y axis tick marks - frequency in GHz (default: 0.5 GHz)
        tlims = self.time[0].mjds()+plot_params.get('tlims',array([0,1e6]))*60.             # min and max time to plot (in min since beginning of obs)
        flims = plot_params.get('flims',array([min(self.f),max(self.f)+1]))                 # min and max frequencies to plot (in Hz)
        
        # clip dynspec to match tlims, flims
        spec,t,f = clip_dynspec(self.spec[pol],[tlims[0],tlims[1],flims[0],flims[1]],self.time.mjds(),self.f)
        t0 = TimeSec(t[0],format='mjds')
        
        ar = ar0*len(t)/len(f)
        
        ## Large plot (entire dynspec, binned) ##
        figure()
        ax=subplot(111,axisbg='k')
        if scale=='log':
            imshow(log10(spec).T,aspect=ar,vmin=log10(smin),vmax=log10(smax),origin='lower',cmap='seismic')
            ds = round(log10(smax)-log10(smin),1)/5         # spacing between colorbar ticks                                                                                      
            ticks = arange(log10(smin),log10(smax)+ds,ds)   # colorbar tick locations                                                                                            
            ticklbls = (10**(ticks+3)).round().astype(int)  # colorbar tick labels                                                                                               
        else: # scale = 'linear'                                                                                                                                                 
            imshow(spec.T,aspect=ar,vmin=smin,vmax=smax,origin='lower',cmap='seismic')
            ds = round(smax,1)/5                            # spacing between colorbar ticks                                                                                    
            ticks = arange(0,smax+ds,ds)                    # colorbar tick locations                                                                                           
            ticklbls = ticks*1000                           # colorbar tick labels                                                                                              
            
        # add colorbar and change labels to show log scale
        cbar = colorbar()
        if pol is 'rc':
            cbar.set_label('Percent Circular Polarization')
            ticks = arange(-1.,1.1,0.2)
            ticklbls = arange(-100,101,20)
            cbar.set_ticks(ticks)
            cbar.set_ticklabels(ticklbls)
        else:
            cbar.set_label('Flux (mJy)')
            cbar.set_ticks(ticks)
            cbar.set_ticklabels(ticklbls)

        # label x axis in minutes                                                                                                                                             
        t_minutes = (t-t0.mjds())/60.  # get time since beginning of observation in minutes                                                                                   
        xmax = max(t_minutes)
        xtick_lbls = arange(0,xmax+dx,dx)
        tick_labels,tick_locs = make_tick_labels(xtick_lbls,t_minutes)
        xticks(tick_locs, tick_labels)
        xlabel('Time in minutes since '+ t0.iso +' UTC')

        # label y axis in GHz                                                                                                                                        
        f_GHz = f/1.e9
        ymin = round(min(f_GHz)/dy)*dy
        ymax = round(max(f_GHz)/dy)*dy
        ytick_lbls = arange(ymin,ymax+dy,dy)
        tick_labels,tick_locs = make_tick_labels(ytick_lbls,f_GHz)
        yticks(tick_locs, tick_labels)
        ylabel('Frequency (GHz)')



##### USER-DEFINED PARAMETERS #####

# scaling factors for color plots
plot_scale_factor = 5

# enable or disable inset plot on binned dynamic spectrum figure
make_inset = False

# parameter filelist is not currently used - will eventually be used to cycle through files to load
filelist = ['UVCet_4L.dynspec','UVCet_4S.dynspec']

# font size for plots
rcParams.update({'font.size': 10})

# number of time and frequency channels to bin together
nt = 8
nf = 8

# properties for large dynspec plots (binned)
ar0 = 1.0 # aspect ratio setting for large dynspec (binned)

# properties for inset dynspec plots (unbinned for flare peak)
ar2 = ar * nt/nf # maintains same aspect ratio as large dynspec
l1min=log10(0.05) # different colorbar to highlight brightest features
l1max=log10(0.3)


##### END USER-DEFINED PARAMETERS #####
'''
dynspecL = Dynspec(params={'filename':'UVCet_4L.dynspec','uniform':True})
dynspecS = Dynspec(params={'filename':'UVCet_4S.dynspec','df':2.e6})

dynspec_Lbin = dynspecL.bin_dynspec(nt=1,nf=2) # bin L-band dynamic spectrum to 2-MHz channels
dynspec_Lbin.f -= 0.5e6   # this is a quick fix; figure out how to address this better later
dynspec = deepcopy(dynspec_Lbin)
dynspec.add_dynspec(dynspecS)
dynspecL,dynspecS = None,None
'''
nf=8
ds_bin = dynspec.bin_dynspec(nt,nf)
ds_bin.plot_dynspec(plot_params={'dx':30})
title('UV Ceti RCP Dynamic Spectrum - 07/18/2015')
savefig('dynspecRCP.pdf',bbox_inches='tight')

ds_bin.plot_dynspec(plot_params={'dx':30,'pol':'ll'})
title('UV Ceti LCP Dynamic Spectrum - 07/18/2015')
savefig('dynspecLCP.pdf',bbox_inches='tight')

ds_bin.spec['rc']=(ds_bin.spec['rr']-ds_bin.spec['ll'])/(ds_bin.spec['rr']+ds_bin.spec['ll'])
ds_bin.plot_dynspec(plot_params={'dx':30,'pol':'rc','scale':'linear','smin':-1.0,'smax':1.0})
title('UV Ceti Degree of Circular Poln - 07/18/2015')
savefig('rc.pdf',bbox_inches='tight')

ar0 = 1.0
ds_bin.plot_dynspec(plot_params={'dx':2,'tlims':array([67.,100.])})
title('UV Ceti RCP - Short Burst - 07/18/2015')
savefig('shortburst.pdf',bbox_inches='tight')
close('all')


'''


##### PLOT DYNAMIC SPECTRA #####



if make_inset:
    # clip out a portion of the dynspec to zoom in on
    Clims = [60,100,1,1.3]
    rrC,tC_min,fC_GHz=clip_dynspec(rr,Clims,t_min,f_GHz)
    llC,tC_min,fC_GHz=clip_dynspec(ll,Clims,t_min,f_GHz)
    
    # the commands below create an inset plot to highlight the flare peak or other features that are best
    #   displayed with different binning and/or a different color scale than the rest of the dynamic spectrum
    #
    # add yellow rectangle to large plot to highlight inset region
    ax.add_patch(Rectangle((imin,jmin),di,dj,edgecolor='#ffff00',fill=False,linewidth=0.5))
    # create inset plot
    ax=axes([0.5,0.5,0.25,0.3],axisbg='k')
    imshow(log10(rrC.T),aspect=ar2,vmin=l1min,vmax=l1max,origin='lower',cmap='seismic')
    title('Flare Peak',color='w',size=10)
    # add yellow rectangle around whole plot (matches rectangle marking inset region on large dynspec)
    ax.add_patch(Rectangle((0,0),shape(rrC)[0]-1,shape(rrC)[1]-1,edgecolor='#ffff00',fill=False,linewidth=2))
    # label axes
    tick_params(labelcolor='w',labelsize=8)
    tick_labels, tick_locs = make_tick_labels(arange(0.,80.,0.5),tC_min)
    xticks(tick_locs, tick_labels) # label x axis time in minutes
    tick_labels, tick_locs = make_tick_labels(arange(1.,4.,0.1),fC_GHz)
    yticks(tick_locs, tick_labels) # label y axis in GHz
    # add colorbar and change labels to show log scale
    cbar = colorbar()
    ticks = arange(l1min,l1max+0.25,0.25)
    cbar.set_ticks(ticks)
    ticklbls = (10**(ticks+3)).round().astype(int)
    cbar.set_ticklabels(ticklbls)
    cbar.ax.yaxis.set_tick_params(labelcolor='w',labelsize=8)

savefig('dynspecRCP.pdf',bbox_inches='tight')

### FIGURE: LCP DYNAMIC SPECTRUM ###
figure()

## Large plot (entire dynspec, binned) ##
ax=subplot(111,axisbg='k')
imshow(log10(ll_bin).T,aspect=ar,vmin=lmin,vmax=lmax,origin='lower',cmap='seismic')
title('UV Ceti Flare LCP Dynamic Spectrum')

# add yellow rectangle to highlight inset region
ax.add_patch(Rectangle((imin,jmin),di,dj,edgecolor='#ffff00',fill=False,linewidth=0.5))

# add colorbar and change labels to show log scale
cbar = colorbar()
cbar.set_label('LCP Flux (mJy)')
ticks = arange(lmin,lmax+0.5,0.5)
cbar.set_ticks(ticks)
ticklbls = (10**(ticks+3)).round(1)
cbar.set_ticklabels(ticklbls)

# label axes: x axis in minutes, y axis in GHz
tick_labels,tick_locs = make_tick_labels(xtick_lbls,t_bin/60.)
xticks(tick_locs, tick_labels)
xlabel('Time in minutes since '+tstring+' UTC')
tick_labels,tick_locs = make_tick_labels(ytick_lbls,f_bin/1.e9)
yticks(tick_locs, tick_labels)
ylabel('Frequency (GHz)')

## Inset plot (flare peak only, unbinned) ##
ax=axes([0.15,0.5,0.25,0.3],axisbg='k')
imshow(log10(llC.T),aspect=ar2,vmin=lmin,vmax=lmax,origin='lower',cmap='seismic')
title('Flare Peak',color='w',size=10)

# add yellow rectangle around whole plot (matches rectangle marking inset region on large dynspec)
ax.add_patch(Rectangle((0,0),shape(rrC)[0]-1,shape(rrC)[1]-1,edgecolor='#ffff00',fill=False,linewidth=2))

# label axes
tick_params(labelcolor='w',labelsize=8)
tick_labels, tick_locs = make_tick_labels(arange(0.,40.,0.5),tC_min)
xticks(tick_locs, tick_labels) # label x axis time in minutes
tick_labels, tick_locs = make_tick_labels(arange(4.,6.,0.1),fC_GHz)
yticks(tick_locs, tick_labels) # label y axis in GHz

# add colorbar and change labels to show log scale
cbar = colorbar()
ticks = arange(lmin,lmax+0.5,0.5)
cbar.set_ticks(ticks)
ticklbls = (10**(ticks+3)).round().astype(int)
cbar.set_ticklabels(ticklbls)
cbar.ax.yaxis.set_tick_params(labelcolor='w',labelsize=8)

savefig('dynspecLCP.pdf',bbox_inches='tight')



# compare linear and log color scales for inset dyn spec - which shows features better?
#  log plot won.
figure()
a = subplot(121)
a.imshow(log10(rrC.T),aspect='auto',vmin=log10(0.05),vmax=log10(0.8),origin='lower',cmap='seismic')
title('Log Scale')
a = subplot(122)
a.imshow(rrC.T*1000,aspect='auto',vmin=1,vmax=500,origin='lower',cmap='seismic')
title('Linear Scale')
savefig('compare.pdf',bbox_inches='tight')


##### FIGURE: PLOT TIME SERIES #####

print 'making plot 2'

# how many seconds to average over when making time series
n1 = 8
n2 = 8
n3 = 8

i_2GHz = ytick_locs[0]
i_4GHz = ytick_locs[2]
tseries_rr = rebin2d_ma(rr[:,i_2GHz:i_4GHz],(1,500))
tseries_ll = rebin2d_ma(ll[:,i_2GHz:i_4GHz],(1,500))
tseries_V = (tseries_rr - tseries_ll)*1000
t2 = rebin1d(t,n2)
V2 = rebin1d_ma(tseries_V[:,0],n2)
t3 = rebin1d(t,n3)
V3 = rebin1d_ma(tseries_V[:,1],n3)

tseries_rr = rebin2d_ma(rr[:,i_4GHz:],(1,978))
tseries_ll = rebin2d_ma(ll[:,i_4GHz:],(1,978))
tseries_V = (tseries_rr - tseries_ll)*1000
t4 = rebin1d(t,n4)
V4 = rebin1d_ma(tseries_V[:,0],n4)
t6 = rebin1d(t,n6)
V6 = rebin1d_ma(tseries_V[:,1],n6)

figure()
hold(True)
plot(t2/60.,V2,'c',label='2-3 GHz')
plot(t3/60.,V3,'g',label='3-4 GHz')
plot(t4[160:]/60.,V4[160:],'r',label='4-5 GHz')
plot(t6/60.,V6,'k',label='6-8 GHz')
xlabel('Time in minutes since '+tstring+' UTC')
ylabel('Flux (mJy)')
title('UV Ceti Flare Stokes V Lightcurve')
#axis([0,tlen/60.,-5,85])
legend(loc='upper left')
savefig('tseries.pdf',bbox_inches='tight')

'''

