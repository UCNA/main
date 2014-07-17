# pd_led_linearity.py
# Author: Simon Slutsky
# Created: 10/01/2013
# 
# Plot fitted linearity constants output from pd_led_pmt_batch file 
# Useful syntax: carp = data[data['Channel'] == 4]


import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from matplotlib.backends.backend_pdf import PdfPages
from math import sqrt
#from ROOT import TCanvas

if __name__ == "__main__":

    plt.ion() #turn on interactive mode
    
    rcParams['figure.figsize'] = 10, 10     #Set default fig size

    imagedir = 'imagesPowerModFullRange'
    filename = 'FitResults.txt'
    path = imagedir + "/" +filename
    
    outputfilename = 'PlotsTogether.pdf'
    outputpath = imagedir + "/" +  outputfilename

    # import data. 
    # TODO: improve this to read the header string from the file
    data = np.genfromtxt(path, skip_header=1, 
                         delimiter = "\t", 
                         names = ['Run','Channel', 'Wavelength','p0',
                                  'p0Err','p1','p1Err','p2','p2Err', 'Chi2'])
    
    fitpars = ['p0', 'p1', 'p2']
    
    run = data['Run']
    data_cut_405 = list()
    data_cut_465 = list()
    data_cut_err_405 = list()
    data_cut_err_465 = list()
    ratioGain = list()
    ratioGainErr = list()
    for chan in range(0,8): # 8 channels
        for wave in (405, 465):
            cutChannel = data['Channel'] == chan
            cutWave = data['Wavelength'] == wave
            cutCond = cutChannel & cutWave
            _data_cut = data[cutCond]
            if wave == 405:
                data_cut_405.append(_data_cut)
            if wave == 465:
                data_cut_465.append(_data_cut)
        #Calculate ratio of Gain parameters and errors
        tmpRatio =  [dc405/dc465 for dc405, dc465 in zip(data_cut_405[chan]['p0'], 
                                                         data_cut_465[chan]['p0'])]
        ratioGain.append(tmpRatio)
        
        delta405 = [dc405Err/dc405 for dc405, dc405Err in zip(data_cut_405[chan]['p0'],
                                                            data_cut_405[chan]['p0Err'])]
        delta465 = [dc465Err/dc465 for dc465, dc465Err in zip(data_cut_465[chan]['p0'],
                                                            data_cut_465[chan]['p0Err'])]
        tmpRatioErr = [ratio*sqrt(d405*d405 + d465*d465) for ratio, d405, d465 in 
                       zip(tmpRatio,delta405, delta465)]
        ratioGainErr.append(tmpRatioErr)

    #fig = plt.figure()
    fig0, (ax0, ax1) = plt.subplots(nrows=2)
    fig1, (ax2, ax3) = plt.subplots(nrows=2)
    fig2, (ax4, ax5) = plt.subplots(nrows=2)
    figChi, (ax6, ax7) = plt.subplots(nrows=2)
    figRatio, ax8 = plt.subplots()
    figures = [fig0, fig1, fig2, figChi, figRatio]
    
    plt.rc('axes', color_cycle=['r', 'g', 'b', 'y'])
    #plt.rc('axes', color_cycle=['c', 'm', 'y', 'k']) # looks terrible
    #ax0.set_color_cycle(['c', 'm', 'y', 'k'])
    #ax1.set_color_cycle(['r', 'g', 'b', 'y'])
    
    ax0.set_title("A   405nm")
    ax1.set_title("A   465nm")    
    ax2.set_title("k   405nm")
    ax3.set_title("k   465nm")    
    ax4.set_title("b   405nm")
    ax5.set_title("b   465nm")   
    ax6.set_title("Chi^2   405nm")
    ax7.set_title("Chi^2   465nm")
    ax8.set_title("Gain Ratio (A_405nm/A_465nm)")

    ax1.set_xlabel('Run Number')
    ax3.set_xlabel('Run Number')
    ax5.set_xlabel('Run Number')
    ax7.set_xlabel('Run Number')
    ax8.set_xlabel('Run Number')
    #ax1.set_ylabel('your y label...')

    #chan = "Ch. "
    _chan = ""
    marks = 4
    for i in range(0,8):
        if i < 4:
            _chan = "E" 
            _chan = _chan + str(i)
        else:
            _chan = "W" 
            _chan = _chan + str(i-4)
        # 405nm
        ax0.errorbar(data_cut_405[i]['Run'], data_cut_405[i]['p0'], yerr=data_cut_405[i]['p0Err'],
                 linestyle='None', marker='o', markersize=marks, label=_chan)
        ax2.errorbar(data_cut_405[i]['Run'], data_cut_405[i]['p1'], yerr=data_cut_405[i]['p1Err'],
                 linestyle='None', marker='o', markersize=marks, label=_chan)
        ax4.errorbar(data_cut_405[i]['Run'], data_cut_405[i]['p2'], yerr=data_cut_405[i]['p2Err'],
                 linestyle='None', marker='o', markersize=marks, label=_chan) 
        ax6.plot(data_cut_405[i]['Run'], data_cut_405[i]['Chi2'], 
                 linestyle='None', marker='o', markersize=marks, label=_chan)
        # 465nm
        ax1.errorbar(data_cut_465[i]['Run'], data_cut_465[i]['p0'], yerr=data_cut_465[i]['p0Err'],
                 linestyle='None', marker='o', markersize=marks, label=_chan)
        ax3.errorbar(data_cut_465[i]['Run'], data_cut_465[i]['p1'], yerr=data_cut_465[i]['p1Err'], 
                 linestyle='None', marker='o', markersize=marks, label=_chan)
        ax5.errorbar(data_cut_465[i]['Run'], data_cut_465[i]['p2'], yerr=data_cut_465[i]['p2Err'],
                 linestyle='None', marker='o', markersize=marks, label=_chan)
        ax7.plot(data_cut_465[i]['Run'], data_cut_465[i]['Chi2'], 
                 linestyle='None', marker='o', markersize=marks, label=_chan)
        #ratio
        ax8.errorbar(data_cut_405[i]['Run'], ratioGain[i], yerr=ratioGainErr[i],
                 linestyle='None', marker='o', markersize=marks, label=_chan)
        
    
    leg0 = ax0.legend(title = "PMT")
    leg1 = ax1.legend(title = "PMT")
    leg2 = ax2.legend(title = "PMT")
    leg3 = ax3.legend(title = "PMT")
    leg4 = ax4.legend(title = "PMT")
    leg5 = ax5.legend(title = "PMT")
    leg6 = ax6.legend(title = "PMT")
    leg7 = ax7.legend(title = "PMT")
    leg8 = ax8.legend(title = "PMT")

    ax0.set_xlim([20800, 24000])
    ax1.set_xlim([20800, 24000])
    ax2.set_xlim([20800, 24000])
    ax3.set_xlim([20800, 24000])
    ax4.set_xlim([20800, 24000])
    ax5.set_xlim([20800, 24000])
    ax6.set_xlim([20800, 24000])
    ax7.set_xlim([20800, 24000])
    ax8.set_xlim([20800, 24000])

    ax0.set_ylim([0, 60])
    ax1.set_ylim([0, 18])
    ax2.set_ylim([0.85, 1.1])
    ax3.set_ylim([0.85, 1.05])
    ax4.set_ylim([-0.0005, 0.0])
    ax5.set_ylim([-0.00015, -0.00008])
    ax8.set_ylim([0, 8])

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
