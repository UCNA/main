#!/usr/bin/python

#Usage: python read_bi_pulser.py <runlist>,
# where <runlist> is a formatted runlist as found in UCNA_AUX.
# If no runlist is supplied, default to ${UCNA_AUX}/UCNA Run Log 2012.txt.

import csv
import sys
from math import sqrt
#from os import path
import os
import string
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

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
    if len(sys.argv) < 2:
#        print "No run list.\nUsage: " +  "read_bi_pulser.py" + " <Run Log File>"
#        print "Using default runlist ../Aux/UCNA_LED_run_log.txt"

        # how does using this run list make any sense?
       # runstring = "../../Aux/UCNA_LED_run_log.txt" 
#        runstring = "../../Aux/UCNA_RUN_LOG_2012_renamed.txt" 
#        runstring = "../../Aux/UCNA Run Log 2012.txt" 
        runstring = os.environ['UCNA_AUX'] + "/UCNA Run Log 2012.txt" 
    else:     
        runstring = sys.argv[1]

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
    
if __name__ == "__main__":
    
    outputfile = PdfPages("BiPulser_21700_21950.pdf")

    centvalarray = return_key_vals("center", "val")
    centerrarray = return_key_vals("center", "err")
    runlist = GetBiRuns()
        
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
        ax.errorbar(runlist, centvalarray[i], 
                    yerr = centerrarray[i], linestyle='None', 
                    marker= 'o', markersize = 4)       

        centvalarray_strip = [c for c in centvalarray[i] if c>0.0]
        ymin, ymax = min(centvalarray_strip), max(centvalarray[i]) 
        ax.set_ylim(ymin - 1, ymax + 1)
        ax.set_xlabel("Run Number")
        ax.set_ylabel("PMT (ADC)")
        ax.set_title(" Bi Pulser " + _chan)

        ax.set_xlim(21703.5, 21950.5)
        outputfile.savefig(figures[i])

    outputfile.close()
    plt.show()

   # 
   # print runlist


    
