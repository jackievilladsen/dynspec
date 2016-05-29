# starting point: raw data and pipeline products -need to apply pipeline caltables and flags - once that is done I can run the general pipeline
# run this file from the directory containing the ms file (e.g. /data/jrv/15A-416/YZCMi/1/L)

msfile = raw_input('Type in ms filename (e.g. YZCMi_1L.ms): ')

### APPLY PIPELINE ###

## APPLYCAL ##

# This is the final applycal command from the casa log from the pipeline.  Since selectdata = False, it will overwrite all previous applycals, so I only need to run this one.
print 'Applying final caltables from VLA CASA pipeline to ms file', msfile
gdir = 'final_caltables/'
gtab = ['gain_curves.g', 'opacities.g', 'requantizergains.g', 'finaldelay.k', 'finalBPcal.b', 'averagephasegain.g', 'finalampgaincal.g', 'finalphasegaincal.g']
for i in range(len(gtab)):
    gtab[i] = gdir + gtab[i]
applycal(vis=msfile,selectdata=False,calwt=[False, False, False, False, False, False, False, False],applymode="calflagstrict",gaintable=gtab)
# flag backup file should be before_applycal_1

## FLAGDATA ##

# based on pipeline casa log, last flag file (created with flagmanager after final applycal) is called "finalflags"
print 'Applying final flags from VLA CASA pipeline to ms file', msfile
flagmanager(vis=msfile,mode='restore',versionname='finalflags')

# inspect flagged data percentages
print 'Generating summary of flagged data'
tbl=flagdata(vis=msfile,mode='summary')
field=tbl['field']
for f in field:
    print f, field[f]['flagged']/field[f]['total']*100, '% flagged'
