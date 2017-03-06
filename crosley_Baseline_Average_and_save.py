import numpy
from datetime import datetime

timeStart = datetime.now()


File = 'target1'
Filename = File+'.ms'

FilePC = File + 'PhaseCorrected'
FilenamePC = FilePC + '.ms'

field = 'EQ_Peg'


########################################################################################################################
# Here is the general procedure that I suggest: https://casa.nrao.edu/docs/CasaRef/ms-Tool.html
#- use ms.msselect() to select a subset of the data (alternatively, check out ms.msseltoindex)
#- use data = ms.getdata() to get a dictionary for only that subset of the data
#- perform averaging on the relevant values in "data" (which will be a dictionary)
# (put all of the above into a loop, to loop through the data as necessary)

#If you wish to store the averaged data back into an ms, I think you'll want to do something like this admittedly kludge-y procedure:
#- use ms.putdata to put the baseline-averaged data back into the ms, pretending that it corresponds to a baseline
#         (so overwrite the data for one of the baselines). Use the shortest or longest uvrange baseline
#- use mstransform to split out that single baseline. mstransform doesn't allow a selection on baseline,
#          but it does allow a selection on uvrange which is why I suggested using one of the extremes
#########################################################################################################################


tb.open(FilenamePC)
scanList = tb.getcol("SCAN_NUMBER")  #getting a list of all the scan values
scanMax=scanList[-1]                #getting the maximum scan number  or use max(scanList) for slower, but garunteed accurate resulsts
tb.close()

tb.open(FilenamePC+'/DATA_DESCRIPTION')
windows = tb.getcol("SPECTRAL_WINDOW_ID")       #getting list of all the spw's
tb.close()

for j in xrange(scanMax+1):       #pick a certain scan number  +1 is needed to reach final scan
    if j in scanList:  #make sure its an observation scan,
        print 'Working on Scan: ', j
        for i in xrange(windows.size):  #pick a certain spw number
            print '    SPW: ', i
            ms.open(FilenamePC)
            ms.msselect({'field':field, 'spw':str(windows[i]) , 'scan' : str(j)})  #select the data cooresponding
            data = ms.getdata(["corrected_amplitude","axis_info","ut"],ifraxis=True)  #ifr = interferometer axis ifr_number, ifr_name, ifr_shortname & baseline
                ###axis_info = https://casa.nrao.edu/docs/CasaRef/ms.getdata.html
                ###ut adds the UT time information to the axis_info time_axis list
            timeCoord = data["axis_info"]["time_axis"]["MJDseconds"]  # UT for UTtime
            freqCoordMHz = data["axis_info"]["freq_axis"]["chan_freq"]/1.0E6  #units in MHz
            baselineCoord = data["axis_info"]["ifr_axis"]["baseline"]
            ampCoord = data["corrected_amplitude"][:,:,:,:]  #4, 128,325, 293   aka coorelations, channels, baselines, times
            npyAmp = numpy.array(ampCoord).mean(2)  #turn ampCoord into a npy array and then average that array along the baseline Axis
                ### all four amplitudes components: XX,XY,YX,YY for me

            npyData = numpy.zeros((freqCoordMHz.shape[0]*timeCoord.shape[0], 4))  # create a emptry matrix rows by 4 collums
            for l in xrange(freqCoordMHz.size): ## loop over freqeuncy
                for n in xrange(timeCoord.size): ## 185 for me
                    stokesI = 0.5 * (npyAmp[0,l,n] + npyAmp[3,l,n])  #1/2 XX + YY
                    stokesV = 0.5 * (npyAmp[1,l,n] - npyAmp[2,l,n])     # 1/2     XY - YX
                    npyData[l*len(xrange(timeCoord.size))+n]=[freqCoordMHz[l][0], timeCoord[n], stokesI, stokesV]  #freq, time, StokesI, stokesV

            numpy.save('./DYNAMIC_SPECTRA/'+'Scan'+str(j)+'_Spw'+str(i),npyData)


print 'Program Completed.  ',
timeEnd = datetime.now()
timeTot = timeEnd - timeStart
timeTot = timeTot.seconds / 60
print 'Total time: ' + str(timeTot) + ' minutes.'

#ms.msselect({'field':field, 'spw':str(1) , 'scan':str(6)})  #select the data cooresponding
