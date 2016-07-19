# pd_led_linearity.py
# Author: Simon Slutsky
# Created: 10/21/2013
# 
# Plot fitted gains output from pd_led_pmt_batch file 
# Useful syntax: carp = data[data['Channel'] == 4]

import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from matplotlib.backends.backend_pdf import PdfPages
from math import sqrt
#from ROOT import TCanvas
import sys

# date plotting stuff
sys.path.append('../../') # find RunPlotter.py
import matplotlib.dates as dates
from RunPlotter import getTimeForRunlist

def read_PD_LED_response():
#    imagedir = '/data4/saslutsky/PulserComp/imagesLEDDebug'
#    imagedir = '/data4/saslutsky/PulserComp/images_06_09_2014'
#    imagedir = '/data1/saslutsky/LEDPulser/images_06_10_2015_16way_separate_wavelength_coeff_20254_23173/'
    imagedir = '/data1/saslutsky/LEDPulser/images_08_26_2015_newPDGainFits_20970_23173/'
    filename = 'GainResults.txt'
    path = imagedir + "/" +filename
    
    outputfilename = 'GainsTogether.pdf'
    outputpath = imagedir + "/" +  outputfilename

    # import data. 
    # TODO: improve this to read the header string from the file
    data = np.genfromtxt(path, skip_header=1, 
                         delimiter = "\t", 
                         names = ['Run',
                                  'A', 'AErr','Mu','MuErr','Sigma','SigmaErr', #'Chi2', 
                                  'BiA', 'BiAErr','BiMu','BiMuErr','BiSigma','BiSigmaErr', 'Chi2'])#'BiChi2'])
    data.sort(order = 'Run') 

    return data

def readTemperatureData():
    data = np.genfromtxt("TavgFile.dat", delimiter = " ", 
                         names = ['Run', 'SCSW', 'SCSE', 'PPM'])

    data.sort(order = 'Run')
    return data
       
if __name__ == "__main__":

    saveBool = int(sys.argv[1])
    dateBool = int(sys.argv[2])

    plt.ion() #turn on interactive mode
    rcParams['figure.figsize'] = 10, 10     #Set default fig size
    
    data = read_PD_LED_response()
    tempdata = readTemperatureData()
    
    fig0, (ax0) = plt.subplots()
    fig1, (ax1) = plt.subplots()
    figures = [fig0, fig1]
    axes = [ax0, ax1]

    plt.rc('axes', color_cycle=['r', 'g', 'b', 'y'])

    ax0.set_title("A")
    ax1.set_title("Mean Response")

    marks = 4

    #cull data
    #   print data[np.where(data['BiAErr'] > data['BiA'])]
    #   data = np.delete(data, np.where(abs(data['BiAErr']) > abs(data['BiA'])), 0)

    dataMu = data['Mu']
    dataA = data['A']
    dataRun = data['Run']
    dataMuErr = data['MuErr']
    dataAErr = data['AErr']

    dataBiMu = data['BiMu']
    dataBiA = data['BiA']
    dataBiRun = data['Run']
    dataBiMuErr = data['BiMuErr']
    dataBiAErr = data['BiAErr']

    tempRun = tempdata['Run']
    tempSCSE = tempdata['SCSE']

    cutcond = np.where( (abs(data['MuErr']) > abs(data['Mu'])/50.) | (data['Mu'] < 0.00001) )
    dataMu =     np.delete(dataMu, cutcond, 0)
    dataA =      np.delete(dataA, cutcond, 0)
    dataMuErr =  np.delete(dataMuErr,  cutcond, 0)
    dataAErr =   np.delete(dataAErr,   cutcond, 0)
    dataRun =    np.delete(dataRun, cutcond, 0)

    Bicutcond = np.where( ( (abs(data['BiMuErr']) > abs(data['BiMu'])/50.) ) | (data['BiMu'] < 0.00001)| (data['BiA'] < 3) | (abs(data['BiAErr']) > abs(data['BiA']/10.) ) ) 
    dataBiMu =     np.delete(dataBiMu,  Bicutcond,  0)
    dataBiA =      np.delete(dataBiA,  Bicutcond,  0)
    dataBiMuErr =  np.delete(dataBiMuErr,  Bicutcond, 0)
    dataBiAErr =   np.delete(dataBiAErr, Bicutcond,  0)
    dataBiRun =    np.delete(dataBiRun, Bicutcond,  0)
    
    tempcutcond = np.where ( tempdata['SCSE'] < 0.00001 ) 
    tempRun = np.delete(tempRun, tempcutcond, 0)
    tempSCSE = np.delete(tempSCSE, tempcutcond, 0)
    print tempcutcond

    dataA = [da-10000 for da in dataA]
    BiMuScale = 1.25
    dataBiMu = [dmb/BiMuScale for dmb in dataBiMu]
    MuScale = 20
    dataMu = [dm+MuScale for dm in dataMu]
    tempSCSE = [tE + 273.15 for tE in tempSCSE]
    
    if not dateBool:
        ax0.errorbar(dataRun, dataA, yerr=dataAErr,
                     linestyle='None', marker='o', markersize=marks)
        ax1.errorbar(dataRun, dataMu, yerr=dataMuErr,
                     linestyle='None', marker='o', markersize=marks)
        ax0.errorbar(dataBiRun, dataBiA, yerr=dataBiAErr,
                     linestyle='None', marker='o', markersize=marks)
        ax1.errorbar(dataBiRun, dataBiMu, yerr=dataBiMuErr,
                     linestyle='None', marker='o', markersize=marks)
        ax1.errorbar(tempRun, tempSCSE, 
                     linestyle='None', marker='o', markersize=marks)

#    ax0.set_xlim([20800, 24000])
#    ax1.set_xlim([20800, 24000])
#    ax0.set_xlabel('Run Number')
#    ax1.set_xlabel('Run Number')

    if dateBool:
        timelist = getTimeForRunlist(dataRun)
        Bitimelist = getTimeForRunlist(dataBiRun)
        temptimelist = getTimeForRunlist(tempRun)

        ax0.errorbar(timelist, dataA, yerr=dataAErr,
                     linestyle='None', marker='o', markersize=marks)
        ax1.errorbar(timelist, dataMu, yerr=dataMuErr, 
                     label = "Normal PD LED Response + " + str(MuScale),
                     linestyle='None', marker='o', markersize=marks)
        ax0.errorbar(Bitimelist, dataBiA, yerr=dataBiAErr,
                     linestyle='None', marker='o', markersize=marks)
        ax1.errorbar(Bitimelist, dataBiMu, yerr=dataBiMuErr, 
                     label = "Bimodal PD LED Response/" + str(BiMuScale),
                     linestyle='None', marker='o', markersize=marks)
        ax1.errorbar(temptimelist, tempSCSE,
                     label = "SCS E Temperature (K)",
                     linestyle='None', marker='o', markersize=marks)

        for ax in axes:
            ax.set_xlabel('Time')


    leg1 = ax1.legend()

 #   ax0.set_ylim([0, 60])
#    ax1.set_ylim([0, 18])
  #  ax2.set_ylim([0.85, 1.1])
   # ax3.set_ylim([0.85, 1.05])
   # ax4.set_ylim([-0.0005, 0.0])
   # ax5.set_ylim([-0.00015, -0.00008])
   # ax8.set_ylim([0, 8])

    if saveBool:
        outputfile = PdfPages(outputpath)
        for f in range(0, len(figures)):
            outputfile.savefig(figures[f])
            outputfile.close()

    plt.show(block=True)     #block=True keeps the plot window open when in interactive mode 



# too hard to do it manually

#run = data['Run']
#_run0 = np.where(data['Channel'] == 0, run, 0]
#run0 = _run0[_run0.nonzero() #strip zeroes from the array

#p0 = data['p0']
#_p00 = np.where(data['Channel'] == 0, p0, 0]
#p00 = _p00[_p00.nonzero()

#for i in range (0, len(run)):
#    if y[i] > 1e8 or y[i] < -1e2:
#        x[i] = 0
#        y[i] = 0
#        print i
