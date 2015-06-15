#!/usr/bin/python

#Usage: python read_bi_pulser.py <LEDcorrbool> <savebool>
# Defaults runlist used to ${UCNA_AUX}/UCNA Run Log 2012.txt.

import csv
import sys
from math import sqrt
#from os import path
import os
import string
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pd_led_linearity_combinedfit import getLEDdata, makeQualityCuts, makeChiCuts

sys.path.append('../../') # find RunPlotter.py
from RunPlotter import getTimeForRunlist
import matplotlib.dates as dates

class ChrisPulserFile:
    def __init__(self, run):
        infile = "/data1/saslutsky/OfficialReplayData/data/Monitors/Run_"
        infile += str(run)
        #    infile += sys.argv[1]
        infile += "/ChrisPulser.txt"

        if not os.path.exists(infile):
            #print "No file found for run " + str(run)
            self.array = 0
            #sys.exit(-1)
        else:
            with open(infile, 'r') as csvfile:
                pulsereader = csv.reader(csvfile, delimiter = '\t')
                pulsearray = list()
                for row in pulsereader:
                    pulsearray.append(row)
                    self.array = pulsearray

class ValwError:
    def __init__(self, val, err):
        self.value = val
        self.error = err
    
# makes an array with all center positions, regardless of tube
def getkey(keystring, array):
    keyarray = list()
    keyarray.append(keystring)

    for row in array:
        for element in row:
            tmp = element.split(' ')
            if tmp[0] == keystring:
                keyarray.append(tmp[2])
    print keyarray

# get all data for one tube
def get_tube_data(sideIn, tubeIn, array):
    tmp_tube_array = list()
    tube_array = list()
    #idstring = str(sideIn) + str(tubeIn)
    #tube_array.append(idstring)
    
    for row in array:
        for element in row:
            tmp = element.split(' ')
            if tmp[0] == 'side' and tmp[2] == sideIn:
                tmp_tube_array.append(row)
                break

   # print tmp_tube_array
    for row in tmp_tube_array:
        for element in row:
            tmpr = element.split(' ')
            if tmpr[0] == 'tube' and tmpr[2] == tubeIn:
                tube_array.append(row)
                break
            
    return tube_array

# get a value with error from a given row
def get_val_w_error(keystring, row):
    returnVal = ValwError(0.0, 0.0)
    tmpstr = keystring
    tmpstrerror = "d" + keystring
    
    for element in row: 
        tmp = element.split(' ') 
        if tmp[0] == tmpstr:
            returnVal.value = float(tmp[2])
        if tmp[0] == tmpstrerror:
            returnVal.error = float(tmp[2])
    
    return returnVal 

# get average value and error over all data from a given tube
def get_average_w_error(keystring, arrayIn):
    average = 0.0
    err = 0.0
    errsquared = 0.0
    cnt = 0;

    for row in arrayIn:
        cnt +=1 

        tmpvalwerror = get_val_w_error(keystring, row)            
        average += tmpvalwerror.value
        err  = tmpvalwerror.error
        errsquared += err*err

    average = average/cnt
    err = sqrt(errsquared)/cnt
    avgwerror = ValwError(average, err)
    
    return avgwerror

def GetBiRuns(runtype = False):
##    if len(sys.argv) < 2:
#        print "No run list.\nUsage: " +  "read_bi_pulser.py" + " <Run Log File>"
#        print "Using default runlist ../Aux/UCNA_LED_run_log.txt"

        # how does using this run list make any sense?
       # runstring = "../../Aux/UCNA_LED_run_log.txt" 
#        runstring = "../../Aux/UCNA_RUN_LOG_2012_renamed.txt" 
#        runstring = "../../Aux/UCNA Run Log 2012.txt" 
    runstring = os.environ['UCNA_AUX'] + "/UCNA Run Log 2012.txt" 
        
##    else:     
##        runstring = sys.argv[1]

    runlist = list()
    runfile = open(runstring, "r")
    for line in runfile:
        if line[0] == "*":
            splitted = string.lstrip(line,"*")
            if runtype == False:
                runlist.append(splitted[0:5])
            if runtype == True:
                runlist.append(splitted[6:9].rstrip())
            
    return runlist

def get_tube_average_w_error(sideIn, tubeIn, arrayIn, keystring):
    tube_data = get_tube_data(sideIn, tubeIn, arrayIn)
    tube_center = get_average_w_error(keystring, tube_data)
    
    return tube_center

def get_tube_numb(sideIn, tubeIn):
    boolp = 0
    if sideIn == "E":
        boopl = 0
    if sideIn == "W":
        boolp = 1
    if (sideIn != "W" and sideIn != "E"):
        print "Incorrect side specified"
        return 0
    
    tubeOut = int(tubeIn) + boolp*4
    return tubeOut 
    
def get_tube_separated_arrays(keyIn, runlistIn):
    # define tubes
    sidelist = ["E","W"]
    tubes = ['0', '1', '2', '3']
        
    #lists of parameters for each tube over all runs
    valarray = list()
    errarray = list()
    for side in sidelist:
        for tube in tubes:
            valarray.append([])
            errarray.append([])

    for run in runlistIn:
        cpf = ChrisPulserFile(run)
        for side in sidelist:
            for tube in tubes:
                #print side + tube
                usetube = get_tube_numb(side, tube)
                if cpf.array:
                    calc_avg_w_error = get_tube_average_w_error(side,
                                                                tube, cpf.array, keyIn)
                    valarray[usetube].append(calc_avg_w_error.value)
                    errarray[usetube].append(calc_avg_w_error.error)
                    #print str(tube_center.value) + " " + str(tube_center.error)
                else:
                    valarray[usetube].append(0.0)
                    errarray[usetube].append(0.0)
                    
    valerrarray = [valarray, errarray]
    return valerrarray


def print_key_arrays(keystringIn):
    # get list of good runs
    runlist = GetBiRuns()
    
    # form arrays 
    valerrarray = get_tube_separated_arrays(keystringIn, runlist)
    
    valarray = valerrarray[0]
    errarray = valerrarray[1]
        
    print valarray
    print errarray

    print "You're Clean!"

def return_key_vals(keystringIn, valerr):
    runlist = GetBiRuns()
    valerrarray = get_tube_separated_arrays(keystringIn, runlist)
    
    if valerr == "val":
        return valerrarray[0]
    if valerr == "err":
        return valerrarray[1]
    else: 
        print "Invalid option in read_bi_pulser::return_key_vals"
        return 0
    
def makeLEDcorrection(LEDdata, BiVal, BiErr, run, tube, corrbool):
    if BiVal == 0.0:
        return 0
    pars = ['p0', 'p1', 'p2']
    parvals = list()
    parerrs = list()

    condrun = LEDdata['Run'] == run
    condtube = LEDdata['Channel'] == tube
    for i in range(3):
        condpar = LEDdata['ParName'] == pars[i]
        condition = condpar & condrun & condtube
        if len(LEDdata[condition]['ParVal']) == 0:
            return 0
        parvals.append(LEDdata[condition]['ParVal'])
        parerrs.append(LEDdata[condition]['ParErr'])
    #    print parvals[i]

    #quadratic correction
    if corrbool == 1:
        correctedBiVal = parvals[0] + parvals[1]*BiVal + parvals[2]*BiVal*BiVal
    if corrbool == 2:
        correctedBiVal = parvals[0] + parvals[1]*BiVal

    #print "###"
#    print BiVal
#    print correctedBiVal
#    print "$$$"
    return correctedBiVal #, do errors later

def writeBiToFile(valarray, errarray, runlist):
    a = len(valarray[0])
    b = len(errarray[0])
    c = len(runlist)
    print a 
    print b
    print c
    if a != b or a != c:
        print "Incommensurate lists in writeBiToFile"
        return -1
    
    BiFile = open('BiPulser.txt','w')
    for i in range(a):
        string = str(runlist[i])
        
        for tube in range(8):
            string = string + " " + str(valarray[tube][i]) + " " + str(errarray[tube][i]) + " "
        
        string += "\n"
        BiFile.write(string)
        
    return 0

if __name__ == "__main__":
    LEDcorrbool = int(sys.argv[1])
    savebool = int(sys.argv[2])
    datebool = int(sys.argv[3])

    basename = "Bi_Pulser_21700_21950"
    if LEDcorrbool == 0:
        outputfile = PdfPages(basename + ".pdf")
    if LEDcorrbool == 1:
        outputfile = PdfPages(basename + "_LEDcorrected.pdf")
    if LEDcorrbool == 2:
        outputfile = PdfPages(basename + "_LEDcorrected_linonly.pdf")

    centvalarray = return_key_vals("center", "val")
    centerrarray = return_key_vals("center", "err")
    runlist = GetBiRuns()

    if savebool:
        writeBiToFile(centvalarray, centerrarray, runlist)
        
    if LEDcorrbool > 0:
        LEDdir = "/data1/saslutsky/LEDPulser/images_05_26_2015_16way_separate_wavelength_coeff_residuals_21650_21950/"
        LEDdata = getLEDdata(LEDdir)
        LEDdata = makeQualityCuts(LEDdata)   # cut runs found bad by eye
        LEDdata = makeChiCuts(LEDdata)       # cut runs with chisq == 0.0
        for i, run in enumerate(runlist):
            for tube in range(0, 8):
 #               print centvalarray[tube][i]
                #centvalarray[tube][i], centerrarray[tube][i] = makeLEDcorrection(LEDdata,  # fix errors later
                centvalarray[tube][i] = makeLEDcorrection(LEDdata, 
                                                          centvalarray[tube][i],
                                                          centerrarray[tube][i],
                                                          int(run),
                                                          tube,
                                                          LEDcorrbool)
#                print centvalarray[tube][i]
                
                
    figures = list()
    for i in range(0, 8):
        if i < 4:
            _chan = "E" 
            _chan = _chan + str(i)
        if i >3 and i < 8:
            _chan = "W" 
            _chan = _chan + str(i-4)
        
        fig, ax = plt.subplots()
        figures.append(fig)
        if not datebool:
            ax.errorbar(runlist, centvalarray[i], 
                        #yerr = centerrarray[i],
                        linestyle='None', 
                        marker= 'o', markersize = 4)       
            ax.set_xlabel("Run Number")

        if datebool: 
            timelist = getTimeForRunlist(runlist)
            ax.xaxis.set_major_formatter(dates.DateFormatter("%m/%d/%y %H:%M"))
            ax.plot_date(timelist, centvalarray[i],
                         #yerr = centerrarray[i],
                         linestyle='None', 
                         marker= 'o', markersize = 4)       
            ax.set_xlabel("Time")
             
        centvalarray_strip = [c for c in centvalarray[i] if c>0.0]
        ymin, ymax = min(centvalarray_strip), max(centvalarray[i]) 
#        ymin, ymax = max(centvalarray_strip) - 400, max(centvalarray[i]) 
        ax.set_ylim(ymin - 1, ymax + 1)
        ax.set_ylabel("PMT (ADC)")
        ax.set_title(" Bi Pulser " + _chan)

        if not datebool:
            ax.set_xlim(21703.5, 21950.5)
        
        if savebool == 1:
            outputfile.savefig(figures[i])

    outputfile.close()
    plt.show()

   # 
   # print runlist


    
