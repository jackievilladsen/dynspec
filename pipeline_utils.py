'''

pipeline_utils.py

Provides functionality used by various pipeline scripts in the dynspec package.

'''

import os
from glob import glob
from recipes import tec_maps
import analysisUtils as au

def get_fields(vis,vishead):
    # return dictionary f with keys 'src','gcal','bpcal' for a given vis
    # assumes that gcal starts with J, bpcal starts with 3C (not very sophisticated)
    fieldlist=vishead(vis=vis,mode='get',hdkey='field')[0]
    f = {}
    for field in fieldlist:
        if field[0:2]=='3C':
            f['bpcal']=field
        elif field[0]=='J':
            f['gcal']=field
        else:
            f['src']=field
    return f

def get_names():
    # read directory location and use it to return project (e.g. 13A-423), src (e.g. ADLeo - not necessarily
    # identical to the src name in the ms though), obs # (e.g. 1), band (e.g. L), and sb (e.g. ADLeo_1L)
    names = {}
    mydir = os.getcwd()
    tmp = mydir[1:].split('/')
    names['proj'] = tmp[2]
    names['src'] = tmp[3]
    names['obs'] = tmp[4]
    if len(tmp)>5:
        names['band'] = tmp[5]
        names['sb'] = names['src']+'_'+names['obs']+names['band']
    return names

def get_phasecen():
    # load phasecen dict from file phasecen.txt (in casa_utils) and find entry for this observation
    # format of one line of phasecen dict: 15A-416_YZCMi_1---J2000 7h44m39.810 3d33m1.91
    
    names = get_names()
    obs_str = names['proj']+'_'+names['src']+'_'+names['obs']
    
    fname = '/data/jrv/casa_utils/dynspec/phasecen.txt'
    f = open(fname)
    lines = f.readlines()
    f.close()
    pcen_dict = {}
    for l in lines:
        [obs_name,coords]=l.rstrip().split('---')
        pcen_dict[obs_name] = coords
    return pcen_dict[obs_str]

def get_plot_dir():
    # return path to current directory + '/plots/', creating this directory if it does not exist yet
    plotdir = os.getcwd() + '/plots/'
    if not os.path.exists(plotdir):
        os.system('mkdir '+plotdir)
    return plotdir

def field_flag_summary(vis,flagdata):
    # print % flagged for each field and return summary of flagged data
    summary = flagdata(vis=vis,mode='summary')
    axis='field'
    for value, stats in summary[axis].iteritems():
        print '%s %s: %5.1f percent flagged' % ( axis, value, 100.*stats['flagged']/stats['total'])
    return summary

def init_cal(vis,gencal,plotcal):
    # attempt to generate antpos, TEC, and requantizer caltables, and return
    #  a gaintable list (for applycal) of the caltable files that were successfully generated
    #  (if no antpos corrections exist, this table will not be generated)
    
    plotdir = get_plot_dir()
    
    # generate antpos corrections
    gencal(vis=vis,caltable=vis+'.antpos',caltype='antpos')
    ''' # currently not using TEC, waiting to hear back from NRAO helpdesk
    # generate TEC corrections
    print 'Generating TEC corrections. WARNING: TEC data are not available until 2 weeks after obs date.'
    tec_image, tec_rms_image = tec_maps.create(vis=vis)
    gencal(vis=vis,caltable=vis+'.tecim',caltype='tecim',infile=tec_image)
    plotcal(caltable=vis+'.tecim',xaxis='time',yaxis='tec',showgui=False,figfile=plotdir+'tec.png')
    print 'Check plots/tec.png for TEC vs. time.'
    '''
    # generate requantizer gain corrections
    gencal(vis=vis,caltype='rq',caltable=vis+'.rq')
    
    # put all existing caltables in gaintable list for use by applycal
    tablist = [vis+'.antpos',vis+'.tecim',vis+'.rq']
    gaintable = []
    for t in tablist:
        if os.path.exists(t):
            gaintable.append(t)
    return gaintable

def amp_plot_max(vis,visstat,datacolumn='data'):
    # returns 2*median visibility amplitude for single baseline/pol, which
    # is a rough estimate of a good maximum plot value for plotting raw visibilities
    # (may need to be more sophisticated for calibrated data)
    stat = visstat(vis,antenna='0&1',correlation='XX',datacolumn=datacolumn)
    return stat[datacolumn.upper()]['median']*2.

def im_params(vis,pblevel=0.05):
    # return a good cell size and imsize (in pixels) for the ms
    cell,imsize,fieldID=au.pickCellSize(vis,imsize=True,pblevel=0.05)
    pixelsize = str(cell) + 'arcsec'
    npixels = imsize[0]
    print 'Pixel size:', pixelsize, '/ Image size:', npixels, 'pixels'
    return pixelsize,npixels


'''

def get_fieldname():
    if proj == '13A-423':
        field_dict = {'UVCet':'V*UVCet','ADLeo':'adleo'}
        return field_dict[src]
    elif proj == '15A-416':
        return src

def get_clean_spws(band):
    # relatively RFI-free spw's that will be used for imaging and subtracting bg srcs
    clean_spw_dict = {'P':'','L':'6,7,11~13','S':'16,19~26','C':'4~13'}
    if proj=='15A-416' and src=='ADLeo' and obs=='4' and band=='L':
        return '10~15'
    return clean_spw_dict[band]

def get_clean_scans(band):
    # read clean (flare-free) scans from file clean_scans.txt in working directory
    # format of a line in this file should be:
    # L   0~13,15~17
    fname = mydir + '/clean_scans.txt'
    try:
        f = open(fname)
        lines = f.readlines()
        f.close()
    except:
        print 'no file exists at',fname,'- assuming all scans are flare-free'
        return ''
    clean_scan_dict = {}
    for l in lines:
        [band_name,scans]=l.rstrip().split()
        clean_scan_dict[band_name] = scans
    try:
        return clean_scan_dict[band]
    except KeyError:
        print 'no line for band', band, 'in file', fname, '- assuming all scans are flare-free'

'''
