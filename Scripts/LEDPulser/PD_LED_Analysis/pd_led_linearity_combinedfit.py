# pd_led_linearity.py
# Author: Simon Slutsky
# Created: 04/20/2015
# 
# Plot fitted linearity constants output from pd_led_pmt_combinedfit 
# Useful syntax: carp = data[data['Channel'] == 4]

# usage: python pd_led_linearity_combinedfit.py <show> <save>

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
    for i in range(0,8): # 8 channels
        data_cut_405 = list()
        data_cut_465 = list()
        for wave in (0, 1):
            for parm in range (0, 3):
                cutChannel = data['Channel'] == int(i)
                cutWave = data['Wavelength'] == int(wave)
                cutParm = data['ParName'] == "t" + str(i) + "p" + str(parm)
                cutCond = cutChannel & cutWave & cutParm
                _data_cut = data[cutCond]
                if wave == 0:
                    data_cut_405.append(_data_cut)
                if wave == 1:
                    data_cut_465.append(_data_cut)

       # print data_cut_405[i]

        tmpFig0, (tmpAx0, tmpAx1, tmpAx2, tmpAx3) = plt.subplots(nrows=4, sharex = True, sharey=False)
#        tmpAxes = [tmpAx0, tmpAx1, tmpAx2, tmpAx3]
        axes = [tmpAx0, tmpAx1, tmpAx2, tmpAx3]
#        axes.append(tmpAxes)
        figures.append(tmpFig0)
        
#        print data_cut_405[2]#['ParVal']
#        sys.exit()
   
        _chan = ""
        marks = 8
        if i < 4:
            _chan = "E" 
            _chan = _chan + str(i)
        else:
            _chan = "W" 
            _chan = _chan + str(i-4)
        for j in range(0, 3): #        for j in range(0, 4): # add chi2 later
        #    print "_________________" + str(j) + "    " + str(i)
#            print data_cut_405[i]
            axes[j].errorbar(data_cut_405[j]['Run'], 
                             data_cut_405[j]['ParVal'], 
                             yerr=data_cut_405[j]['ParErr'], 
                             linestyle='None', marker='^',
                             markersize=marks, label=_chan, color = 'Black')
            for tl in axes[j].get_yticklabels():
                tl.set_color('Black')
     
            twinax = axes[j].twinx()
            twinax.errorbar(data_cut_465[j]['Run'],
                            data_cut_465[j]['ParVal'], 
                            yerr=data_cut_465[j]['ParErr'],
                            linestyle='None', marker='v', 
                            markersize=marks, label=_chan, color = 'Blue')
            for tl in twinax.get_yticklabels():
                tl.set_color('Blue')
                
            axes[j].set_ylabel("p" + str(j))

      #  axes[i][3].plot(data_cut_405[i]['Run'], data_cut_405[i]['Chi2'], 
#                        linestyle='None', marker='o', markersize=marks, label=_chan)
     #   axes[i][3].plot(data_cut_465[i]['Run'], data_cut_465[i]['Chi2'], 
 #                       linestyle='None', marker='o', markersize=marks, label=_chan)

 #       axes[i][0].set_ylabel("p0   405nm")
 #       axes[i][3].set_ylabel("p0   465nm")    
 #       axes[i][1].set_ylabel("p1   405nm")
 #       axes[i][4].set_ylabel("p1   465nm")    
 #       axes[i][2].set_ylabel("p2   405nm")
 #       axes[i][5].set_ylabel("p2   465nm")   
 #       axes[i][6].set_ylabel("Chi^2/nu   405nm")
 #       axes[i][7].set_ylabel("Chi^2/nu  465nm")
 #       axes[i][8].set_ylabel("Gain Ratio (A_405nm/A_465nm)")
  
 #       axes[i][0].set_ylabel("p2   405nm")
 #       axes[i][2].set_ylabel("p2   465nm")

 #       axes[i][3].set_ylabel("Chi^2/nu  465nm")

 #       axes[i][0].set_yscale("log",nonposy='clip')
 #       axes[i][3].set_yscale("log",nonposy='clip')
        

#        axes[i][6].set_title(_chan)

 #       for m in range(0, 9):
 #           axes[i][m].set_xlim([21300, 24000])
        
        axes[3].set_ylabel("Chi^2/nu")        
        axes[0].set_title(_chan)
#        xmin, xmax = axes[0].get_xlim()
        xmin, xmax = data_cut_405[j]['Run'][0], data_cut_405[j]['Run'][-1] 
        axes[0].set_xlim([xmin - 0.5, xmax + 0.5])
        axes[3].set_xlabel('Run Number')

        
    #    axes[i][0].set_ylim([0.0, 60.]) 
    #    axes[i][3].set_ylim([0.0, 18.])   
    #    axes[i][1].set_ylim([0.85, 1.1])
    #    axes[i][4].set_ylim([0.85, 1.05])
    #    axes[i][2].set_ylim([-0.0005, 0.0])
    #    axes[i][5].set_ylim([-0.00015, -0.00008])
    #    axes[i][8].set_ylim([0, 8])

#        axes[i][5].set_xlabel('Run Number')

#        axes[i][8].set_xlabel('Run Number')
 
#    for f in range (0, 16):
    if savebool: 
        for f in range (0, 8):
            outputfile.savefig(figures[f])
    
    outputfile.close()
    if showbool:
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
