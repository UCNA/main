# 10/15/2014
# Script to plot p2/d(p2) for an estimate of how non-linear the LED fits are.

# dtype=[('myint', '<i8'), ('myfloat', '<f8'), ('mystring', '|S5')])
#dtype="i8,f8,S5",

import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from matplotlib.backends.backend_pdf import PdfPages
from math import sqrt
import sys
from scipy import signal as sig

sys.path.append('../../') # find RunPlotter.py
from RunPlotter import getTimeForRunlist
import matplotlib.dates as dates

def plotAllLEDLinParms_only405(minrun = 20970, sigthresh = 50): #21086
# NOT IN USE AT THE MOMENT!!!
#    linearitydata = np.genfromtxt("/data4/saslutsky/PulserComp/images_10_16_2014_allruns/FitResults.txt",
    linearitydata = np.genfromtxt("/data4/saslutsky/PulserComp/images_11_24_2014_10rangemin/FitResults.txt",
                                  delimiter = "\t", 
                                  names = ['Run','Channel', 'Wavelength','p0',
                                           'p0Err','p1','p1Err','p2','p2Err', 
                                           'Chi2', 'NDF'])
    fig = plt.figure()
    ax0 = fig.add_subplot(111)
    
    for tube in range(8):
        cutChannel = linearitydata['Channel'] == tube
        cutWave = linearitydata['Wavelength'] == 405
        cutCond = cutChannel & cutWave
        _data_cut = linearitydata[cutCond]
        
        runlist  = list()
        datalist = list()
        for row in _data_cut:
            run = row['Run']
            val = row['p2']
            err = row['p2Err']
            if val != 0.0:
                sigfactor = val/err
                if abs(sigfactor) < sigthresh:
                     if run > minrun: 
                         datalist.append(sigfactor)
                         runlist.append(run)
                else:
                    print "Run " + str(run) + ": abs(Sigfactor) > " + str(sigthresh) + " . Omitting."
            else: 
                print "Run " + str(run) + ": Parameter was 0"

        markie = "D"
        if tube < 7:
            markie = "o"
        mylabel = "PMT " + str(tube)
        ax0.errorbar(runlist, datalist,
                      linestyle = 'None', marker = markie, markersize = 4,
                      label = mylabel)

#    ax0.set_xlim([21000, 23500])
    ax0.legend(title = "405 nm LED")
    fig.text(0.5, 0.04, 'Run Number', ha='center', va='center')
    fig.text(0.06, 0.5, 'Non-linearity Significance', ha='center', va='center', rotation='vertical')

    plt.show()


def plotAllLEDLinParms_only405_separated(filename, plot_date = False, sigthresh = 50, minrun = 20970): #21086
#    linearitydata = np.genfromtxt("/data4/saslutsky/PulserComp/images_10_16_2014_allruns/FitResults.txt",

    rcParams['figure.figsize'] = 10, 12

    linearitydata = np.genfromtxt(filename,
                                  delimiter = "\t", 
                                  names = ['Run','Channel', 'Wavelength','p0',
                                           'p0Err','p1','p1Err','p2','p2Err', 
                                           'Chi2', 'NDF'])
 
    figE, (ax0, ax1, ax2, ax3) = plt.subplots(4, sharex = True, sharey=False)
    figW, (ax4, ax5, ax6, ax7) = plt.subplots(4, sharex = True, sharey=False)
    axes = [ax0, ax1, ax2, ax3, ax4, ax5, ax6, ax7]

    figEp2, (axp0, axp1, axp2, axp3) = plt.subplots(4, sharex = True, sharey=False)
    figWp2, (axp4, axp5, axp6, axp7) = plt.subplots(4, sharex = True, sharey=False)
    axesp2 = [axp0, axp1, axp2, axp3, axp4, axp5, axp6, axp7]

    for tube in range(8):
        cutChannel = linearitydata['Channel'] == tube
        cutWave = linearitydata['Wavelength'] == 405
        cutCond = cutChannel & cutWave
        _data_cut = linearitydata[cutCond]
        
        runlist  = list()
        datalist = list()

        p2list = list()
        p2errlist = list()

        for row in _data_cut:
            run = row['Run']
            val = row['p2']
            err = row['p2Err']
            if val != 0.0:
                sigfactor = val/err
                if abs(sigfactor) < sigthresh:
                     if run > minrun: 
                         if abs(val) < 0.5:
                             if err < 0.5:
                                 datalist.append(sigfactor)
                                 runlist.append(run)
                                 p2list.append(val)
                                 p2errlist.append(err)
                else:
                    print "Run " + str(run) + ": abs(Sigfactor) > " + str(sigthresh) + " . Omitting."
            else: 
                print "Run " + str(run) + ": Parameter was 0"

        markie = "D"
#        if tube < 7:
#        markie = "o"
        mylabel = "PMT " + str(tube)
#        datalist = sig.savgol_filter(datalist, 101, 2)
#        p2list = sig.savgol_filter(p2list, 101, 2)

        if not plot_date:
            axes[tube].errorbar(runlist, datalist,
                                linestyle = 'None', marker = markie, markersize = 4,
                                label = mylabel)
            axes[tube].legend()

            axesp2[tube].errorbar(runlist, p2list, yerr = p2errlist,
                                  linestyle = 'None', marker = markie, markersize = 4,
                                  label = mylabel)
            axesp2[tube].legend()
            xaxisLabel = "Run Number"
            
        if plot_date:
            timelist = getTimeForRunlist(runlist)
            axes[tube].xaxis.set_major_formatter(dates.DateFormatter("%m/%d/%y"))
            axes[tube].plot_date(timelist, datalist,
                                linestyle = 'None', marker = markie, markersize = 4,
                                label = mylabel)
            axes[tube].legend()

            axesp2[tube].xaxis.set_major_formatter(dates.DateFormatter("%m/%d/%y"))
            axesp2[tube].plot_date(timelist, p2list, 
                                   linestyle = 'None', marker = markie, markersize = 4,
                                   label = mylabel)
            axesp2[tube].legend()
            xaxisLabel = "Time"

#    ax0.set_xlim([21000, 23500])
    ax0.set_title("East Tubes, 405 nm LED")
    ax4.set_title("West Tubes, 405 nm LED")
    
    axp0.set_title("East Tubes, 405 nm LED")
    axp4.set_title("West Tubes, 405 nm LED")
    
  #  figE.text(0.5, 0.04, 'Run Number', ha='center', va='center')
    figE.text(0.5, 0.04, xaxisLabel, ha='center', va='center')
    figE.text(0.06, 0.5, 'Non-linearity Significance', ha='center', va='center', rotation='vertical')
#    figW.text(0.5, 0.04, 'Run Number', ha='center', va='center')
    figW.text(0.5, 0.04, xaxisLabel, ha='center', va='center')
    figW.text(0.06, 0.5, 'Non-linearity Significance', ha='center', va='center', rotation='vertical')

#    figEp2.text(0.5, 0.04, 'Run Number', ha='center', va='center')
    figEp2.text(0.5, 0.04, xaxisLabel, ha='center', va='center')
    figEp2.text(0.06, 0.5, 'p2 (Non-linear Term)', ha='center', va='center', rotation='vertical')
#    figWp2.text(0.5, 0.04, 'Run Number', ha='center', va='center')
    figWp2.text(0.5, 0.04, xaxisLabel, ha='center', va='center')
    figWp2.text(0.06, 0.5, 'p2 (Non-linear Term)', ha='center', va='center', rotation='vertical')

    plt.show()


def plotAllLEDLinParms():
    linearitydata = np.genfromtxt("/data4/saslutsky/PulserComp/images_10_16_2014_allruns/FitResults.txt",
                                  delimiter = "\t", 
                                  names = ['Run','Channel', 'Wavelength','p0',
                                           'p0Err','p1','p1Err','p2','p2Err', 
                                           'Chi2', 'NDF'])
    fig = plt.figure()
    ax0 = fig.add_subplot(211)
    ax1 = fig.add_subplot(212)
    
    for tube in range(8):
        cutChannel = linearitydata['Channel'] == tube
        for wave in (405, 465):
            cutWave = linearitydata['Wavelength'] == wave
            cutCond = cutChannel & cutWave
            _data_cut = linearitydata[cutCond]

            runlist  = list()
            datalist = list()
            for row in _data_cut:
                run = row['Run']
                val = row['p2']
                err = row['p2Err']
                if val != 0.0:
                    sigfactor = val/err
                    if sigfactor > -10000.0:
                        datalist.append(sigfactor)
                        runlist.append(run)
                    else:
                        print "Run " + str(run) + ": Sigfactor < -10000. Omitting."
                else: 
                    print "Run " + str(run) + ": Parameter was 0"

            axis = ax0 if wave == 405 else ax1
            markie = "o"
            if tube < 7:
                markie = "D"
            mylabel = "PMT " + str(tube)
            axis.errorbar(runlist, datalist,
                              linestyle = 'None', marker = markie, markersize = 5,
                              label = mylabel)

    ax0.set_xlim([21000, 23500])
    ax1.set_xlim([21000, 23500])
    ax0.legend(title = "405 nm LED")
    ax1.legend(title = "465 nm LED")
    fig.text(0.5, 0.04, 'Run Number', ha='center', va='center')
    fig.text(0.06, 0.5, 'Non-linearity Significance', ha='center', va='center', rotation='vertical')

    plt.show()

def plotAverageLEDLinParms():
    linearitydata = np.genfromtxt("../ELOGPics/AverageLEDLinearityParms_pol2.txt",
                                  delimiter = "\t", dtype = "i8,i8,i8,i8,S5,f8,f8,f8",
                                  names = ['Run_start', 'Run_end', 'tube', 'wave',
                                           'parm', 'val', 'meanerr', 'stddev'] )

    runsegments = [21086,21274, 21291,21679, 21704, 21914, 22215, 22294,22437]

    fig = plt.figure()
    ax0 = fig.add_subplot(211)
    ax1 = fig.add_subplot(212)

    p2data = linearitydata[linearitydata['parm'] == 'p2']
    for tube in range(8):
        cutChannel = p2data['tube'] == tube
        for wave in (405, 465):
            cutWave = p2data['wave'] == wave
            cutCond = cutChannel & cutWave
            _data_cut = p2data[cutCond]

            datalist = list()
            errlist = list()
            for run in runsegments:
                val      = _data_cut[_data_cut['Run_start'] == run]['val']
                meanerr  = _data_cut[_data_cut['Run_start'] == run]['meanerr']
                stddev   = _data_cut[_data_cut['Run_start'] == run]['stddev']
                
                # find significance factor
                quaderr = sqrt(meanerr**2 + stddev**2)
                sigfactor = val/quaderr
                datalist.append(sigfactor)

            axis = ax0 if wave == 405 else ax1
            markie = "o"
            if tube < 7:
                markie = "D"
            axis.errorbar(runsegments, datalist,
                              linestyle = 'None', marker = markie, markersize = 7, 
                              label = "PMT " + str(tube))

    ax0.set_xlim([21000, 23000])
    ax1.set_xlim([21000, 23000])
#    ax.set_xlim([21000, 23000])
    ax0.legend(title = "405 nm LED")
    ax1.legend(title = "465 nm LED")
    fig.text(0.5, 0.04, 'Run Number', ha='center', va='center')
    fig.text(0.06, 0.5, 'Non-linearity Significance', ha='center', va='center', rotation='vertical')


#    ax.set_xlabel("Run Number")
#    ax.set_ylabel("Non-linearity Significance")

    plt.show()


if __name__== "__main__":
#    filename = "/data4/saslutsky/PulserComp/images_10_22_2014_allruns/FitResults"
#    filename = "/data4/saslutsky/PulserComp/images_11_24_2014_10rangemin/FitResults"
#    filename = "/data4/saslutsky/PulserComp/images_12_15_2014_fixBeta/FitResults"
#    filename = "/data4/saslutsky/PulserComp/images_01_14_2015_scaleE/FitResults"
     filename = "/data1/saslutsky/PulserComp/images_09_18_2018_light_versus_ADCcounts_22534_22557/FitResults_Combined.txt
     
    try: 
        if sys.argv[1] == 'pe':
            filename += "_PE_PMT.txt"
    except IndexError: 
        filename += ".txt"
    
#    plotAllLEDLinParms_only405_separated(filename, True, 40)
    plotAllLEDLinParms_only405_separated(filename, False, 40)
