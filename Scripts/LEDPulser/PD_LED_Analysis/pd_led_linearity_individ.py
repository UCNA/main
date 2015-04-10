# pd_led_linearity.py
# Author: Simon Slutsky
# Created: 10/01/2013
# 
# Plot fitted linearity constants output from pd_led_pmt_batch file 
# Useful syntax: carp = data[data['Channel'] == 4]


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pylab import *
from math import sqrt
#from ROOT import TCanvas

if __name__ == "__main__":

#    plt.ion() #turn on interactive mode
    
    rcParams['figure.figsize'] = 10, 10     #Set default fig size
    plt.rc('axes', color_cycle=['r', 'g', 'b', 'y'])

    # import data. 
    # TODO: improve this to read the header string from the file
#    data = np.genfromtxt("imagesPowerMod/FitResults.txt", skip_header=1, 
#                         delimiter = "\t", 
#                         names = ['Run','Channel', 'Wavelength','p0',
#                                  'p0Err','p1','p1Err','p2','p2Err', 'Chi2')]\

    basedir = "/data4/saslutsky/PulserComp/images_04_09_2015_21927_21939/"
    outputfile = PdfPages(basedir + "FitResultsIndividual.pdf")
#    outputfile = PdfPages('imagesPowerMod/ParmPlotsIndividual.pdf')
#    data = np.genfromtxt("/data4/saslutsky/PulserComp/images_04_09_2015_21927_21939/FitResults.txt", 
    data = np.genfromtxt(basedir + "FitResults.txt", 
                         skip_header=1, 
                         delimiter = "\t", 
                         names = ['Run','Channel', 'Wavelength','p0',
                                  'p0Err','p1','p1Err','p2','p2Err',
                                  'p3', 'p3Err', 'Chi2', 'NDF',
                                  'RangeStart', 'RangeEnd']  )
    
    #cut bad runs
#    cut1 = (data['Run']!=21599)
#    cut2 = (data['Run']!=21601)
#    cut3 = (data['Run']!=21602)
 #   cutsAll = cut1 & cut2 & cut3
    cutruns = [21599, 21601, 21602, 21933]
    cutsAll = True
    for i in cutruns:
        _cut = (data['Run']!=i)
        cutsAll = cutsAll & _cut
    data = data[cutsAll]
    #testing
#    print data[data['Run']==21933]
#    sys.exit()

    fitpars = ['p0', 'p1', 'p2']
    
    run = data['Run']
    data_cut_405 = list()
    data_cut_465 = list()
    data_cut_err_405 = list()
    data_cut_err_465 = list()
    ratioGain = list()
    ratioGainErr = list()
    axes = list()
    figures = list()
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
        
        tmpFig0, (tmpAx0, tmpAx1, tmpAx2, tmpAx3) = plt.subplots(nrows=4, sharex = True, sharey=False)
        tmpAxes = [tmpAx0, tmpAx1, tmpAx2, tmpAx3]
        axes.append(tmpAxes)
        figures.append(tmpFig0)
        
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
#        for j in range(0, 3):
        for j in range(0, 1):
#            fitstring = fitpars[j]
            fitstring = 'p2'
            errstring = fitstring + 'Err'
            axes[i][j].errorbar(data_cut_405[i]['Run'], data_cut_405[i][fitstring], yerr=data_cut_405[i][errstring],
                         linestyle='None', marker='o', markersize=marks, label=_chan)
#            axes[i][j+3].errorbar(data_cut_465[i]['Run'], data_cut_465[i][fitstring], yerr=data_cut_465[i][errstring],
            axes[i][j+2].errorbar(data_cut_465[i]['Run'], data_cut_465[i][fitstring], yerr=data_cut_465[i][errstring],
                         linestyle='None', marker='o', markersize=marks, label=_chan)
#        axes[i][6].plot(data_cut_405[i]['Run'], data_cut_405[i]['Chi2'], 
        axes[i][1].plot(data_cut_405[i]['Run'], data_cut_405[i]['Chi2'], 
                        linestyle='None', marker='o', markersize=marks, label=_chan)
 #       axes[i][7].plot(data_cut_465[i]['Run'], data_cut_465[i]['Chi2'], 
        axes[i][3].plot(data_cut_465[i]['Run'], data_cut_465[i]['Chi2'], 
                        linestyle='None', marker='o', markersize=marks, label=_chan)
        #ratio
   #     axes[i][8].errorbar(data_cut_405[i]['Run'], ratioGain[i], yerr=ratioGainErr[i],
#                 linestyle='None', marker='o', markersize=marks, label=_chan)

 #       axes[i][0].set_ylabel("p0   405nm")
 #       axes[i][3].set_ylabel("p0   465nm")    
 #       axes[i][1].set_ylabel("p1   405nm")
 #       axes[i][4].set_ylabel("p1   465nm")    
 #       axes[i][2].set_ylabel("p2   405nm")
 #       axes[i][5].set_ylabel("p2   465nm")   
 #       axes[i][6].set_ylabel("Chi^2/nu   405nm")
 #       axes[i][7].set_ylabel("Chi^2/nu  465nm")
 #       axes[i][8].set_ylabel("Gain Ratio (A_405nm/A_465nm)")
  
        axes[i][0].set_ylabel("p2   405nm")
        axes[i][2].set_ylabel("p2   465nm")
        axes[i][1].set_ylabel("Chi^2/nu  405nm")
        axes[i][3].set_ylabel("Chi^2/nu  465nm")

 #       axes[i][0].set_yscale("log",nonposy='clip')
 #       axes[i][3].set_yscale("log",nonposy='clip')
        
        axes[i][0].set_title(_chan)
#        axes[i][6].set_title(_chan)

 #       for m in range(0, 9):
 #           axes[i][m].set_xlim([21300, 24000])
        

        xmin, xmax = axes[i][0].get_xlim()
        axes[i][0].set_xlim([xmin - 0.5, xmax + 0.5])
        
    #    axes[i][0].set_ylim([0.0, 60.]) 
    #    axes[i][3].set_ylim([0.0, 18.])   
    #    axes[i][1].set_ylim([0.85, 1.1])
    #    axes[i][4].set_ylim([0.85, 1.05])
    #    axes[i][2].set_ylim([-0.0005, 0.0])
    #    axes[i][5].set_ylim([-0.00015, -0.00008])
    #    axes[i][8].set_ylim([0, 8])

#        axes[i][5].set_xlabel('Run Number')
        axes[i][3].set_xlabel('Run Number')
#        axes[i][8].set_xlabel('Run Number')
 
#    for f in range (0, 16):
    for f in range (0, 8):
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
