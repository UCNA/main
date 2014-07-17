# Run QuickBiPeakStudy for all runs on Sources List
# typcial command: 
#./QuickBiPeakStudy /data/ucnadata/2011/rootfiles/full20818.root
# use standard output to dump tube Bi peak positions to a text file
# Calls compiled root script, eg:
#### ./QuickBiPeakStudy /data/ucnadata/2011/rootfiles/full20818.root#"

#### USAGE: python QuickBiPeakBatch.py >>  txt.txt

## Simon Slutsky 05/28/2014

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def createSourceRunList(filename):
    sourcerunlist = list()
    f = open(filename, "r")
    bi_in_use = 0 #boolean flag to parse Run Log
    for line in f:      
        if line[0:8] == "@sources": # toggle to use run or not
            if line.find("Bi") > 0:
                bi_in_use = 1
            else:
                bi_in_use = 0
        if line.find("SourcesCal") > 0:
            if line[0] == "*" and bi_in_use == 1:    # only take "good" runs, in which Bi source is in use
                run = line[1:6] # parse run list
                sourcerunlist.append(run)

    return sourcerunlist

def makeAllBiPeakPlot(filename, showb=0):
    data = np.genfromtxt(filename, delimiter = "\t", 
                  names = ['Run', 'm1', 'e1', 'm2', 'e2', 'm3', 'e3', 'm4', 'e4',
                                  'm5', 'e5', 'm6', 'e6', 'm7', 'e7', 'm8', 'e8'])
    
    fig = plt.figure()
    plt.subplot(111)
    colors = ('b', 'g', 'r', 'c', 'm', 'y', 'k')
    runlist = data['Run']
    for i in range(1,9):
        means = data['m' + str(i)]
        errs  = data['e' + str(i)]
        plt.errorbar(runlist, means, yerr = errs, color = colors[i%7], 
                     linestyle = "None", marker = 'o', markersize =4 )
        ax = plt.gca()
    
    ax.set_ylim(0, 2000)
    if showb:
        plt.show()
    return fig
 
def makeBiTubePlot(filename, tube, showb=0):
    data = np.genfromtxt(filename, delimiter = "\t", 
                  names = ['Run', 'm1', 'e1', 'm2', 'e2', 'm3', 'e3', 'm4', 'e4',
                                  'm5', 'e5', 'm6', 'e6', 'm7', 'e7', 'm8', 'e8'])
    
    fig = plt.figure()
    plt.subplot(111)
    runlist = data['Run']
    means = data['m' + str(tube)]
    errs  = data['e' + str(tube)]
    plt.errorbar(runlist, means, yerr = errs, linestyle = "None", marker = 'o', markersize =4 )
    ax = plt.gca()
    
    ax.set_ylim(0, 3000)
    if showb:
        plt.show()
    return fig
   
if __name__ == "__main__":
    basepath = "/home/saslutsky/UCNA/UCNAReplay_052214/UCNA/"
    filepath = "/Aux/UCNA Run Log 2012.txt"
    
    sourcerunlist = createSourceRunList(basepath+filepath)

    basecmdstring = "./QuickBiPeakStudy /data/ucnadata/2011/rootfiles/full"
    
    for run in sourcerunlist:
        cmdstring = basecmdstring + str(run) + ".root"
#        print cmdstring
        os.system(cmdstring)

