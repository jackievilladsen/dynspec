# run in project directory (e.g., /data/jrv/15A-416)

import os

#srclist = ['ADLeo','UVCet','EQPeg','EVLac','YZCMi']
# for testing purposes
srclist = ['EVLac']
bandlist = ['S']

mydir = os.getcwd()
tmp = mydir[1:].split('/')
proj = tmp[2]

for src in srclist:
    for band in bandlist:
        os.chdir(src+'/'+band)
        obslist = 

        os.chdir('../..')
