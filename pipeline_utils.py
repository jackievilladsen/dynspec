'''

pipeline_utils.py

Provides functionality used by various pipeline scripts in the dynspec package.

'''

import os
from glob import glob
#from recipes import tec_maps
import subprocess
import numpy as np

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

def extract_ms(msfile):
    # give a tarfile named msfile+'.tar', this will run tar -xv and then
    # move the contents to msfile and msfile.flagversions

    tarfile = msfile+'.tar'
    flagfile = msfile+'.flagversions'
    
    if os.path.exists(msfile):
        print msfile, 'already exists, exiting extract_ms without re-extracting.', \
              'Delete',msfile,'before running extract_ms if you want a fresh version.'
        return
    
    print 'extracting ms from tarfile...'
    os.system('tar -xf '+tarfile)
    filelist = glob('lustre/*/*/*/*/*/*')
    for f in filelist:
        if f[-3:] == '.ms':
            os.system('mv '+f+' '+msfile)
        elif f[-13:]=='.flagversions':
            os.system('mv '+f+' '+flagfile)
    os.system('rm -rf lustre')    

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
        if l[0] != '#':
            [obs_name,coords]=[s.rstrip() for s in l.split('---')]
            pcen_dict[obs_name] = coords
    return pcen_dict[obs_str]

def get_phasecenP():
    # load phasecen dict from file phasecenP.txt (in casa_utils) and find entry for this observation
    # format of one line of phasecen dict: SB---Ipos---Vpos
    # e.g.:    15A-416_YZCMi_1---J2000 7h44m39.810 3d33m1.91---J2000 7h44m39.815 3d33m1.90
    
    names = get_names()
    obs_str = names['proj']+'_'+names['src']+'_'+names['obs']
    
    fname = '/data/jrv/casa_utils/dynspec/phasecenP.txt'
    f = open(fname)
    lines = f.readlines()
    f.close()
    pcenI_dict = {}
    pcenV_dict = {}
    for l in lines:
        if l[0] != '#':
            [obs_name,coordsI,coordsV]=[s.rstrip() for s in l.split('---')]
            pcenI_dict[obs_name] = coordsI
            pcenV_dict[obs_name] = coordsV
    return pcenI_dict.get(obs_str,''),pcenV_dict.get(obs_str,'')

def get_clean_scans():
    # read clean (flare-free) scans and spws from clean_scans.txt in dynspec code dir
    # format of a line in this file should be:
    # SB                 clean_scans  clean_spws
    # 15A-416_ADLeo_4L---0~13,15~17---10~15
    # clean_scans and clean_spws can both be blank, must separate by '---', but don't need 2nd entry, e.g.
    # 15A-416_ADLeo_4L--- ---10~15
    # 15A-416_ADLeo_4L---0~13,15~17
    
    names = get_names()
    obs_str = names['proj']+'_'+names['sb']
    
    fname = '/data/jrv/casa_utils/dynspec/clean_scans.txt'
    f = open(fname)
    lines = f.readlines()
    f.close()
    clean_scan_dict = {}
    clean_spw_dict = {}
    for l in lines:
        tmp=[s.rstrip() for s in l.split('---')]
        try:
            obs_name = tmp[0]
            clean_scan_dict[obs_name] = tmp[1]
            clean_spw_dict[obs_name] = tmp[2]
        except:
            pass
#    return clean_scan_dict, clean_spw_dict
    return clean_scan_dict.get(obs_str,''),clean_spw_dict.get(obs_str,'')

'''
def get_clean_spws(band):
    # relatively RFI-free spw's that will be used for imaging and subtracting bg srcs
    clean_spw_dict = {'P':'','L':'6,7,11~13','S':'16,19~26','C':'4~13'}
    if proj=='15A-416' and sb=='ADLeo_4L':
        return '10~15'
    return clean_spw_dict[band]
'''

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

def amp_plot_max(vis,visstat,datacolumn='data',field=''):
    # returns 2*median visibility amplitude for single baseline/pol, which
    # is a rough estimate of a good maximum plot value for plotting raw visibilities
    # (may need to be more sophisticated for calibrated data)
    stat = visstat(vis,antenna='0&1;2&3;4&5',correlation='XX',datacolumn=datacolumn,field=field)
    return stat[datacolumn.upper()]['median']*2.

def im_params(vis,pblevel=0.1):
    # return a good cell size and imsize (in pixels) for the ms
    import analysisUtils as au
    cell,imsize,fieldID=au.pickCellSize(vis,imsize=True,pblevel=pblevel,npix=3)
    if imsize[0]<5000:
        cell,imsize,fieldID=au.pickCellSize(vis,imsize=True,pblevel=pblevel,npix=5)
    pixelsize = str(cell) + 'arcsec'
    npixels = imsize[0]
    print 'Pixel size:', pixelsize, '/ Image size:', npixels, 'pixels'
    return pixelsize,npixels

def get_params(key,pipe_params):
    # helper function to load pipeline parameters for dynspec_pipeline.py
    pipe_params_default = {'overwrite':False,'makebig':False,'interactive':False}
    # overwrite: set this to True if you want existing .src.ms, .src.ms.tar, and .small.ms
    #            to be overwritten (otherwise pipeline skips steps for which files have already been produced)
    # makebig: set to True to run bg subtraction and tbavg on full-resolution src ms
    # interactive: set to True to make cleaning interactive
    try:
        return pipe_params[key]
    except:
        return pipe_params_default[key]

def get_ants(vis,msmd):
    # returns a list of antennas (ea-# style names) from measurement set vis
    # must pass msmd (the ms metadata tool)
    msmd.open(vis)
    antlist = msmd.antennanames()
    msmd.done()
    return antlist

def get_refant(vis,msmd):
    # return a refant near array center that is in subarray in vis]
    # must pass msmd (the ms metadata tool)
    preferred_refant = ['ea06','ea22','ea05','ea25'] # identified by hand
    antlist = get_ants(vis,msmd)
    for a in preferred_refant:
        if a in antlist:
            return a
    return antlist[0] # should not happen

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

def load_preferred_files():
    # returns a dictionary whose keys are observation names (e.g. '15A-416_UVCet_3S') and whose
    # entries are filenames pointing to the preferred dynspec file
    #
    # does not contain entries for observations that aren't in best_dynspec_file.txt (i.e. which
    # use the default dynspec from the pipeline)
    fname = '/data/jrv/casa_utils/dynspec/best_dynspec_file.txt'
    f = open(fname)
    lines = f.readlines()
    f.close()
    file_dict = {}
    for l in lines:
        if l[0] != '#':
            [obs_name,filename]=[s.rstrip() for s in l.split('---')]
            file_dict[obs_name] = filename
    return file_dict

def load_ds_filelist():
    # returns a dictionary whose keys are observation file names (e.g. '/data/jrv/15A-416/UVCet/3') and whose
    #  entries are lists of dynspec files (e.g. one for L band, one for S band)
    preferred_file_dict = load_preferred_files()
    filelist = {}
    cmd = 'ls -d /data/jrv/*/*/*/*/*.tbavg.ms.dynspec | cut -d / -f 1-6 | uniq'
    obslist = subprocess.check_output(cmd, shell=True).rstrip().split('\n') # get list of every obs that has dynspec files
    for obs in obslist:
        obsname = obs.split('/')[3] + '_' + obs.split('/')[4] + '_' + obs.split('/')[5]
        cmd = 'ls -d ' + obs + '/*/*.tbavg.ms.dynspec'
        flist = subprocess.check_output(cmd,shell=True).rstrip().split('\n')
        filelist[obs] = flist
        for i in np.arange(len(flist)):
            f = flist[i]
            band = f.split('/')[6]
            band_obsname = obsname+band
            flist[i] = preferred_file_dict.get(band_obsname,f)
    return filelist

def load_burst_filelist(band='',project=''):
    # returns a dictionary whose keys are observation file names (e.g. 'data/jrv/15A-416/UVCet/3') from
    # epochs where there was a burst, and whose entries are lists of dynspec files (e.g. one for L band, one for S band)
    # options: band='' --> uses burst_epochs.txt
    #          band='P' --> uses burst_epochs_Pband.txt (all these epochs are covered w/in
    #                       burst_epochs.txt as well since all have bursts at higher freq)
    #          band='withP' --> uses burst_epochs_with_Pband.txt
    #                          (all epochs with detected bursts that have
    #                           P band coverage, even if burst not detected
    #                           in P band)
    # project='': If project is set to '15A-416' or '13A-423', will only load 
    #             burst filenames from that project
    ds_filelist = load_ds_filelist()
    
    if band=='P':
        fname = '/data/jrv/casa_utils/dynspec/burst_epochs_Pband.txt'
    elif band=='withP':
        fname = '/data/jrv/casa_utils/dynspec/burst_epochs_with_Pband.txt'
    else:
        fname = '/data/jrv/casa_utils/dynspec/burst_epochs.txt'
    f = open(fname)
    lines = f.readlines()
    f.close()
    file_dict = {}
    for l in lines:
        l = l.rstrip()
        if l=='':
            continue
        if project not in l:
            continue
        if l[0] != '#':
            file_dict[l] = ds_filelist[l]
    return file_dict

def load_band_filelist(band='L'):
    # return list of filenames of all ds files for a certain band
    #
    ds_dict = load_ds_filelist()
    filelist = [item for sublist in ds_dict.values() for item in sublist if item.split('/')[6]==band.upper()]
    #flatten the list of lists and check if each item matches the desired band
    return filelist
        
