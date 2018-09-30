'''
refresh_dynspec_files.py

Purpose: Re-run the step of going from tbavg.ms to tbavg.ms.dynspec for all observations. This is
         useful because I found a bug in dyn_spec (which reads the dynspec out of tbavg.ms) and so
         want to redo just this step.
'''
import dynspec.ms2dynspec
reload(dynspec.ms2dynspec)
from dynspec.pipeline_utils import load_ds_filelist
from dynspec.ms2dynspec import tbavg2dsfile

ds_filelist = load_ds_filelist()

failed_list = []
for obs in ds_filelist:
    obs_filelist = ds_filelist[obs]
    for ds_file in obs_filelist:
        tbavg_ms_file = ds_file[:-8]
        if os.path.exists(tbavg_ms_file):
            ds_file = tbavg2dsfile(tbavg_ms_file)
        else:
            print tbavg_ms_file, 'does not exist'
            failed_list.append(tbavg_ms_file)

print 'Failed list:'
print failed_list
# failed on VLBA data (except for UVCet_3X...no it wasn't in ds_filelist) --> need to do those separately
# (although they shouldn't have much flagging so it shouldn't be as much of an issue)

# to do: find code used to create VLBA dynspecs in the first place
# find VLBA dynspecs
# write script to recreate VLBA dynspecs
ls -d /data/jrv/BV071/*/*/*tbavg.ms*
/data/jrv/BV071/ADLeo/3/ADLeo_3X.tbavg.ms          /data/jrv/BV071/UVCet/3/UVCet_3X.tbavg.ms
/data/jrv/BV071/ADLeo/3/ADLeo_3X.tbavg.ms.dynspec  /data/jrv/BV071/UVCet/3/UVCet_3X.tbavg.ms.dynspec
/data/jrv/BV071/ADLeo/4/ADLeo_4X.tbavg.ms          /data/jrv/BV071/UVCet/4/UVCet_4X.tbavg.ms
/data/jrv/BV071/ADLeo/4/ADLeo_4X.tbavg.ms.dynspec  /data/jrv/BV071/UVCet/4/UVCet_4X.tbavg.ms.dynspec
/data/jrv/BV071/ADLeo/5/ADLeo_5X.tbavg.ms          /data/jrv/BV071/UVCet/5/UVCet_5X.tbavg.ms
/data/jrv/BV071/ADLeo/5/ADLeo_5X.tbavg.ms.dynspec  /data/jrv/BV071/UVCet/5/UVCet_5X.tbavg.ms.dynspec

# Conclusion: I want to use selfcal'd versions of the data sets for dynspecs anyways, so don't worry about re-creating these dynspecs for the moment.

