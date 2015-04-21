# pd_led_linearity.py
# Author: Simon Slutsky
# Created: 04/20/2015
# 
# Plot fitted linearity constants output from pd_led_pmt_combinedfit 
# Useful syntax: carp = data[data['Channel'] == 4]

# usage: python pd_led_linearity_combinedfit.py <show> <save> <LED>

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pylab import *
from math import sqrt
#from ROOT import TCanvas
import sys

if __name__ == "__main__":
    
    showbool = int(sys.argv[1]) # boolean to control plotting
    savebool = int(sys.argv[2]) # boolean to control saving
    LEDbool = int(sys.argv[3]) # control which LED is plotted 0: 405, 1: 465, 2: both

    rcParams['figure.figsize'] = 10, 10     #Set default fig size
    plt.rc('axes', color_cycle=['r', 'g', 'b', 'y'])

    # import data. 

    basedir = "/data4/saslutsky/PulserComp/images_04_17_2015_8waysimulfit_21927_21939/"
    outputfile = PdfPages(basedir + "FitResultsCombined.pdf")
    data = np.genfromtxt(basedir + "FitResults_Combined.txt", 
                         skip_header=0, 
                         delimiter = "\t", 
                         names = ['Run','Channel', 'Wavelength','ParName',
                                  'ParVal', 'ParErr'],
                         dtype = "int, int, int, S4, float, float" )
    
    #cut bad runs
    cutruns = [21599, 21601, 21602, 21933, 21295, 21298,21312, 21318, 21288, 21933]
    cutsAll = True
    for i in cutruns:
        _cut = (data['Run']!=i)
        cutsAll = cutsAll & _cut
    data = data[cutsAll]

    run = data['Run']
    axes = list()
    figures = list()
    for i in range(0,9): # 8 channels + PD 
        data_cut_405 = list()
        data_cut_465 = list()
        for wave in (0, 1):
            for parm in range (0, 3):
                cutChannel = data['Channel'] == int(i)
                cutWave = data['Wavelength'] == int(wave)
                if i < 8:
                    cutParm = data['ParName'] == "t" + str(i) + "p" + str(parm)
                else: 
                    cutParm = data['ParName'] == "PD" + "p" + str(parm)
                cutCond = cutChannel & cutWave & cutParm
                _data_cut = data[cutCond]
                if wave == 0:
                    data_cut_405.append(_data_cut)
                if wave == 1:
                    data_cut_465.append(_data_cut)
        tmpFig0, (tmpAx0, tmpAx1, tmpAx2, tmpAx3) = plt.subplots(nrows=4, sharex = True, sharey=False)
#        tmpAxes = [tmpAx0, tmpAx1, tmpAx2, tmpAx3]
        axes = [tmpAx0, tmpAx1, tmpAx2, tmpAx3]
        figures.append(tmpFig0)
   
        _chan = ""
        marks = 8
        if i < 4:
            _chan = "E" 
            _chan = _chan + str(i)
        else:
            _chan = "W" 
            _chan = _chan + str(i-4)
        if i == 8:
            _chan = "PD" 
        for j in range(0, 3): #        for j in range(0, 4): # add chi2 later
            if (LEDbool == 0 or LEDbool == 2):
                axes[j].errorbar(data_cut_405[j]['Run'], 
                                 data_cut_405[j]['ParVal'], 
                                 yerr=data_cut_405[j]['ParErr'], 
                                 linestyle='None', marker='^',
                                 markersize=marks, label=_chan, color = 'Black')
            if (LEDbool == 1 or LEDbool == 2):
                axes[j].errorbar(data_cut_465[j]['Run'],
                                 data_cut_465[j]['ParVal'], 
                                 yerr=data_cut_465[j]['ParErr'],
                                 linestyle='None', marker='v', 
                                 markersize=marks, label=_chan, color = 'Blue')
            axes[j].set_ylabel("p" + str(j))

        axes[3].set_ylabel("Chi^2/nu")        
        axes[0].set_title(_chan)
        xmin, xmax = data_cut_405[j]['Run'][0], data_cut_405[j]['Run'][-1] 
        axes[0].set_xlim([xmin - 0.5, xmax + 0.5])
        axes[3].set_xlabel('Run Number')
   
    if savebool: 
        for f in range (0, 9):
            outputfile.savefig(figures[f])
    
    outputfile.close()
    if showbool:
        plt.show(block=True)     #block=True keeps the plot window open when in interactive mode 
