# pd_led_pmt_gain.py
# Author: Simon Slutsky
# Created: 10/21/2013
# 
# Plot fitted pmt gains output from pd_led_pmt_batch file 

import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from matplotlib.backends.backend_pdf import PdfPages
from math import sqrt
#from ROOT import TCanvas
import sys

sys.path.append('../../') # find RunPlotter.py
from RunPlotter import getTimeForRunlist
import matplotlib.dates as dates

def ReadLEDFile():
    imagedir = '/data1/saslutsky/LEDPulser/images_06_10_2015_16way_separate_wavelength_coeff_20254_23173/'
    filename = 'PMTGainResults.txt'
    path = imagedir + "/" +filename
    
    outputfilename = 'GainsTogether.pdf'
    outputpath = imagedir + "/" +  outputfilename

    # import data. 
    # TODO: improve this to read the header string from the file
    readData = np.genfromtxt(path, skip_header=1, 
                         delimiter = "\t", 
                         names = ['Run','tube', 'A', 'AErr','Mu',
                                  'MuErr','Sigma','SigmaErr', 'Chi2'])
    readData.sort(order = 'Run')

    return readData, outputpath

def getLEDDataforTube(tubeIn):
    dat = ReadLEDFile()
    cutDat =  dat[dat['tube'] == tubeIn]

    return cutDat

# copy/pasted from pd_led_gain.py - needs updates
if __name__ == "__main__":

    savebool = sys.argv[1]
    datebool = sys.argv[2]

    plt.ion() #turn on interactive mode
    rcParams['figure.figsize'] = 10, 10     #Set default fig siz

    data, outputpath = ReadLEDFile()
    
    means = data['Mu']
    meanErrs = data['MuErr']
    average = sum(means[:]/len(means))
    print "AVERAGE = " + str(average)
    
    # calculate "gain" arbitrarily normalized to 260 
    gains = [m/average for m in means]
    gainsErr = [mErr/average for mErr in meanErrs]

    #fig = plt.figure()
    fig0, (ax0) = plt.subplots()
    fig1, (ax1) = plt.subplots()
    fig2, (ax2) = plt.subplots()
    fig3, (ax3) = plt.subplots()
    figChi, (ax4) = plt.subplots()
    figures = [fig0, fig1, fig2, fig3, figChi]
    
    plt.rc('axes', color_cycle=['r', 'g', 'b', 'y'])
    axes = [ax0, ax1, ax2, ax3, ax4]

    ax0.set_title("A")
    ax1.set_title("Mean Response")
    ax2.set_title("Relative Gain")    
    ax3.set_title("Sigma")
    ax4.set_title("Chi2")

    marks = 4
    if not datebool:
        ax0.errorbar(data['Run'], data['A'], yerr=data['AErr'],
                     linestyle='None', marker='o', markersize=marks)
        ax1.errorbar(data['Run'], data['Mu'], yerr=data['MuErr'],
                     linestyle='None', marker='o', markersize=marks)
        ax2.errorbar(data['Run'], gains, gainsErr,
                     linestyle='None', marker='o', markersize=marks)
        ax3.errorbar(data['Run'], data['Sigma'], yerr=data['SigmaErr'],
                     linestyle='None', marker='o', markersize=marks)
        ax4.errorbar(data['Run'], data['Chi2'], 
                 linestyle='None', marker='o', markersize=marks)
    
     #   ax0.set_xlim([20800, 24000])
     #   ax1.set_xlim([20800, 24000])
     #   ax2.set_xlim([20800, 24000])
     #   ax3.set_xlim([20800, 24000])
     #   ax4.set_xlim([20800, 24000])
        ax0.set_xlabel('Run Number')
        ax1.set_xlabel('Run Number')
        ax2.set_xlabel('Run Number')
        ax3.set_xlabel('Run Number')
        ax4.set_xlabel('Run Number')

    if datebool: 
        timelist = getTimeForRunlist(data['Run'])
        ax0.errorbar(timelist, data['A'], yerr=data['AErr'],
                     linestyle='None', marker='o', markersize=marks)
        ax1.errorbar(timelist, data['Mu'], yerr=data['MuErr'],
                     linestyle='None', marker='o', markersize=marks)
        ax2.errorbar(timelist, gains, gainsErr,
                     linestyle='None', marker='o', markersize=marks)
        ax3.errorbar(timelist, data['Sigma'], yerr=data['SigmaErr'],
                     linestyle='None', marker='o', markersize=marks)
        ax4.errorbar(timelist, data['Chi2'], 
                     linestyle='None', marker='o', markersize=marks)
        for ax in axes:
            ax.xaxis.set_major_formatter(dates.DateFormatter("%m/%d/%y"))
            ax.xaxis_date()
            ax.set_xlabel('Time')
            
    sepfigs = list()
    sepaxes = list()
    for i in range(0,8):
        tmpfig, (tmpax) = plt.subplots()
        if i < 4:
            chan = "E" + str(i)
        else: 
            chan = "W" + str(i%4)
        tmpax.set_title("PMT Response to LED Gain Pulse: " + chan)
        tmpax.set_ylabel("PMT (ADC)")
        sepfigs.append(tmpfig)
        sepaxes.append(tmpax)

        data_cut = data[[data['tube']  == i]]
#        print data_cut

        if not datebool:
            tmpax.errorbar(data_cut['Run'], data_cut['Mu'], yerr=data_cut['MuErr'],
                           linestyle='None', marker='o', markersize=marks)
            tmpax.set_xlabel('Run Number')
        if datebool:
            timelist = getTimeForRunlist(data_cut['Run'])
            tmpax.errorbar(timelist, data_cut['Mu'], yerr=data_cut['MuErr'],
                           linestyle='None', marker='o', markersize=marks)
            tmpax.set_xlabel('Time')

    if savebool:
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
