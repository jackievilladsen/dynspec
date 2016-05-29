# Module that assists with producing dynamic spectra in CASA.
# The primary routines are tbavg which averages all baselines together
# and dyn_spec which puts the visibilites into a numpy array.
#
# To be used from within CASA.
#
# Stephen Bourke
# Caltech, Oct 2013.
#
# Modified version:
# Jackie Villadsen
# 5/29/2016

'''
Notes:
    
1) WEIGHTS: 'flat' seems to work best. Ultimately I may want to write my own code to calculate weights per channel and
   baseline, where the weights are determined by std dev over time.  For this, it would
   be helpful to be able to specify distinct timeranges for calculating and applying
   the weights.  This might yield a better RMS when avg'ing over baselines than flat
   weights.
   
2) DISCREPANCIES VS PLOTMS AVG: There are discrepancies between the output of tbavg
   and plotms avg'ing over all baselines. The discrepancies are fractionally largest
   when the astronomical signal is small.  The result of plotms avg'ing over baselines
   does change if you change the weights, so plotms is using the weights somehow - but
   maybe not in the same way as tbavg?  The plotms results are a closer match to tbavg
   with weights than tbavg without weights.  WHY IS THERE A DIFFERENCE? Maybe plotms handles flags differently?

3) SPW AVG'ING: plotms avgspw=True no longer works (does not avg spw's together) on post-tbavg ms.  Why?

4) TBAVG vs TBAVG2: tbavg2, which uses split for avg'ing, runs ~100x faster than tbavg for large data sets,
   but it produces slightly less clean dynamic spectra - why? Use keepflags=False and see if this is fixed.
'''

import time as timemod
from taskinit import tbtool,ms
import numpy

def table(tablename, readonly=True):
    """Return a new table object with the specified table open"""
    t = tbtool()
    t.open(tablename, nomodify=readonly)
    return t

def subtable(tab, subname):
    try:
        label, path = tab.getkeyword(subname).split()
    except:
        raise KeyError, 'Subtable not found'
    if label != 'Table:':
        raise ValueError, 'Does not appear to be a sub table'
    return table(path)

def unique_col_values(tab, colname):
    """Return a numpy array of the unique values in the column"""
    return tab.query('', columns='DISTINCT '+colname).getcol(colname)

def plt(msname):
    intab = table(msname)
    x = []
    y = []
    for t in range(20):
        for spw in range(16):
            for v in intab.getcell('DATA', t * 16 + spw).ravel():
                x.append(t)
                y.append(v.real)
    return x,y

def tbavg(msname, savename, weight_mode='', datacolumn='DATA'):
    """
    Create a new MS with all baselines averaged together.
    
    Options for weight_mode are '' (use WEIGHT column for weights), 'spc' (use
    WEIGHT_SPECTRUM column for weights) and 'flat' (use no weights).
    
    Parameter datacolumn can be used to select which column to average over (default
    is DATA, the raw data column).  The column must have the format of visibility
    data: if they exist, columns CORRECTED_DATA and MODEL_DATA fit this criterion.
    No matter which column is selected, the results will be written to the DATA
    column in the output ms.
    """
    
    # track how long this takes (in seconds)
    t1=timemod.time()
    
    # create a table tool with the original ms open
    intab = table(msname)
 
    # create a table using first row from each unique combo of time,spw and save that
    # as the output ms
    # DATA_DESC_ID means spw
    subt = intab.taql('SELECT FROM %s ORDERBY DISTINCT TIME,DATA_DESC_ID' % msname)
    subt.copy(savename, True)
      
    # open the output ms with permission to write
    outtab = table(savename, readonly=False)
      
    # remove columns CORRECTED_DATA and MODEL_DATA from output ms if they exist
    # (since we do not want to spend extra time to average over those columns)
    for colname in ['CORRECTED_DATA','MODEL_DATA']:
        if colname in outtab.colnames():
            outtab.removecols(colname)  
    
    # number of rows in output ms
    outrows = outtab.nrows()
    
    # replace antenna and UVW columns in output ms so all baselines are 0-1
    # and all UVW is [0,0,0]
    outtab.putcol('ANTENNA1',numpy.zeros(outrows))
    outtab.putcol('ANTENNA2',numpy.ones(outrows))
    outtab.putcol('UVW',numpy.zeros([3,outrows]))
    
    # check to make sure that the column we are averaging over
    # (input paramater datacolumn) has the correct format (so we can average
    # over the CORRECTED_DATA coumn but not the ANTENNA1 column, for instance)
    data_shape = intab.getcell('DATA', 0).shape
    output_shape = intab.getcell(datacolumn, 0).shape
    if data_shape != output_shape:
        raise ValueError, 'DATA and %s shapes differ' % datacolumn
    
    # When we read in a column of data it will give us an N-d array - N depends
    # on whether the data have 0, 1, or 2 dimensions (0 dimensions would be for
    # only one corr and only one channel per spw, 2 dimensions is normal for VLA).
    # baseline_axis is the number of dimensions of data in a single cell, so that
    # we know that the next dimension is the one that iterates over baselines in
    # our subtables below.
    baseline_axis = len(data_shape)
    
    for i in range(outrows):
        if i % 100 == 0:
            print i, '/', outrows
        time = outtab.getcell('TIME', i)
        data_desc_id = outtab.getcell('DATA_DESC_ID', i)
        subt = intab.query('TIME=%.8f AND DATA_DESC_ID=%d AND FLAG_ROW=False' % (time, data_desc_id), columns='%s,WEIGHT,FLAG,WEIGHT_SPECTRUM' % datacolumn)
        if subt.nrows() == 0:
            outtab.putcell('FLAG_ROW', i, True)
            continue
        outtab.putcell('FLAG_ROW',i,False)
        data = subt.getcol(datacolumn)
        flag = subt.getcol('FLAG')
        good_mask = numpy.logical_not(flag)
        data *= good_mask
        outflag=flag.all(baseline_axis)
        outtab.putcell('FLAG',i,outflag)
        if weight_mode == 'flat':
            avg_data = numpy.average(data,baseline_axis)
            outtab.putcell('DATA', i, avg_data)
            # replace WEIGHT column with ones
            old_wt = outtab.getcell('WEIGHT',i)
            new_wt = numpy.ones(old_wt.shape)
            outtab.putcell('WEIGHT',i,new_wt)
        else:
            if weight_mode == 'spc':
                weights = subt.getcol('WEIGHT_SPECTRUM')
            else:
                weights = numpy.repeat(subt.getcol('WEIGHT'), data_shape[1], axis=0)
                weights.shape = data.shape
            weights *= good_mask
            avg_data, weights_sum = numpy.ma.average(data, baseline_axis, weights, True)
            outtab.putcell('DATA', i, avg_data.data)
            outtab.putcell('WEIGHT', i, weights_sum.data.sum(axis=min(2,len(weights_sum.data.shape))-1))
            if weight_mode == 'spc':
                outtab.putcell('WEIGHT_SPECTRUM', i, weights_sum.data)
    intab.unlock()
    outtab.unlock()
    subt.unlock()
    t2=timemod.time()
    print 'Duration:', t2-t1, 's'

def tbavg2(split, msname, savename, weight_mode='', datacolumn='DATA'):
    """
    tbavg2(split, msname, savename, weight_mode='', datacolumn='DATA')
    
    Create a new MS with all baselines averaged together.
    
    Options for weight_mode are '' (use WEIGHT column for weights) and 'flat' (use no weights).
    
    Parameter datacolumn can be used to select which column to average over (default
    is DATA, the raw data column).  The column (or columns) must have the format of visibility
    data: options are DATA, MODEL_DATA, CORRECTED_DATA, FLOAT_DATA, LAG_DATA, and/or all.

    Have to pass the function split to tbavg2 - I tried using ms.split but it doesn't give the same result.
    """
    
    # track how long this takes (in seconds)
    t1=timemod.time()
    
    # create a table tool with the original ms open
    intab = table(msname,readonly=False)

    # save copy of antenna columns of original ms
    ant1 = intab.getcol('ANTENNA1')
    ant2 = intab.getcol('ANTENNA2')
    
    # replace antenna columns in original ms so all baselines are 0-1
    nrows = intab.nrows()
    intab.putcol('ANTENNA1',numpy.zeros(nrows))
    intab.putcol('ANTENNA2',numpy.ones(nrows))
    
    # write over weights in original ms if weight_mode='flat'
    if weight_mode.lower()=='flat':
        wt = intab.getcol('WEIGHT')
        intab.putcol('WEIGHT',numpy.ones(wt.shape))

    # get minimum integration time and make sure our split averaging time is less than that
    try:
        interval = intab.getcol('INTERVAL')
        dt = min(interval) * 1e-3
    except:
        dt = 0.01
    timebin = str(dt)+'s'
#    timebin='0.1s'
    
    # split and avg over time < integration time
    split(vis=msname,outputvis=savename,datacolumn=datacolumn,timebin=timebin,keepflags=False)
    
    # put original antenna columns back in original ms
    intab.putcol('ANTENNA1',ant1)
    intab.putcol('ANTENNA2',ant2)
    
    # put original weights back in original ms if weight_mode='flat'
    if weight_mode.lower()=='flat':
        intab.putcol('WEIGHT',wt)
    
    intab.unlock()
    
    t2=timemod.time()
    print 'tbavg duration:', t2-t1, 's'

def interval(times):
    """Return a likely at the integration time"""
    t = numpy.sort(times)
    delta = numpy.around(t[1:] - t[:-1], 2)
    return numpy.median(delta)

def interpolate(times):
    """Return a regular array with missing elements filled in"""
    t = numpy.sort(times)
    delta = numpy.around(t[1:] - t[:-1], 2) # Intervals, 2 significant figures
    t_int = numpy.median(delta)
    add_positions = numpy.floor(numpy.fmax(delta / t_int, t_int))
    out_times = numpy.empty(shape=(add_positions.sum()+1,), dtype=times.dtype)
    i=0
    for j in range(len(add_positions)):
        for k in range(add_positions[j]):
            out_times[i] = t[j] + k*t_int
            i+=1
    out_times[i] = t[j+1]
    return out_times

def dyn_spec(msname, ddids=None, interpolate_times=False):
    tab = table(msname)
    if ddids is None:
        ddids = numpy.arange(subtable(tab, 'DATA_DESCRIPTION').nrows())
    else:
        ddids = numpy.array(ddids)
    times = unique_col_values(tab, 'TIME')
    if interpolate_times:
        out_times = interpolate(times)
    else:
        out_times = times
    tmp_data = tab.getcell('DATA', 0)
    pols, freqs = tmp_data.shape
    arr = numpy.zeros(shape=(len(out_times), pols, len(ddids)*freqs), dtype=tmp_data.dtype)
    for i in range(tab.nrows()):
        if tab.getcell('FLAG_ROW', i):
            continue
        ddid = tab.getcell('DATA_DESC_ID', i)
        time = tab.getcell('TIME', i)
        t_indx = numpy.where(out_times==time)[0][0]
        try:
            d_indx = numpy.where(ddids==ddid)[0][0]
        except IndexError:
            # This DATA_DESC is not requested in the output.
            continue
        if tab.iscelldefined('WEIGHT_SPECTRUM', i):
            arr[t_indx,:,d_indx*freqs:(d_indx+1)*freqs] = tab.getcell('DATA', i) * (tab.getcell('WEIGHT_SPECTRUM', i) > 0)
        else:
            arr[t_indx,:,d_indx*freqs:(d_indx+1)*freqs] = tab.getcell('DATA', i)

    # Generate array of frequency axis values
    ddesc = subtable(tab, 'DATA_DESCRIPTION')
    spw = subtable(tab, 'SPECTRAL_WINDOW')
    ftype = spw.getcell('CHAN_FREQ', 0).dtype
    freq_values = numpy.zeros(shape=(len(ddids)*freqs,), dtype=ftype)
    for i in range(len(ddids)):
        spw_id = ddesc.getcell('SPECTRAL_WINDOW_ID', ddids[i])
        freq_values[i*freqs:(i+1)*freqs] = spw.getcell('CHAN_FREQ', spw_id)

    tab.unlock()
    ddesc.unlock()
    spw.unlock()
    return {'data': arr, 'freqs': freq_values, 'times': out_times}

def dynspec_LR(msname):
    spec = dyn_spec(msname)
    nu = spec['freqs']
    t = spec['times']
    
    pol_list = ['RR','RL','LR','LL']
    pol_data = {}
    for i in range(4):
        pol_data[pol_list[i]] = spec['data'][:,i,:]
    return t,nu,pol_data

def threshold(arr, lsigma, usigma=None):
    """Clip the array to the specified multiple of the RMS."""
    arr = arr.copy()
    if usigma == None:
        usigma = lsigma
        lsigma *= -1
    rms = arr.std()
    arr[arr < lsigma * rms] = lsigma * rms
    arr[arr > usigma * rms] = usigma * rms
    return arr
