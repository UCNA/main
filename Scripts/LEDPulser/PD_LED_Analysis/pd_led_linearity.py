# pd_led_linearity.py
# Author: Simon Slutsky
# Created: 10/01/2013
# 
# Plot fitted linearity constants output from pd_led_pmt_batch file 
# Useful syntax: carp = data[data['Channel'] == 4]

# 09/23/2014: Updated to work with constrained fits file, too
# usage: python pd_led_linearity.py <bool> 
# for bool = 0 --> normal fit
#            1 --> constrained fit
# 10/01/2014: Kill constrained fit boolean
# add a save bool to control plot output to file
# ---------  usage: python pd_led_linearity.py <savebool> 

import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from matplotlib.backends.backend_pdf import PdfPages
from math import sqrt
#from ROOT import TCanvas
import sys

sys.path.append('/home/saslutsky/UCNA/UCNAReplay_052214/UCNA/Scripts/LEDPulser/BetaSpectrumFitting')
#from BetaEnergyScaleBatch import truncate_runs_before
from RunLogUtilities import identifySourceSegments

def find_segment_average(runlist, datalist, errorlist, segment_start, segment_end):
# find all runs bookended by segment_start, segment_end and find the weighted average data value
# wt. mean = sum(x_i*dx_i^-2)/sum(dx_i^-2) 
# err = 1/sqrt(sum(dx_i^-2)) --> This is only for uncorrelated measurements
### track avg statistical error and stddev  (max stat error is not realistic)

#    for i in range(len(L)):    ===   for i, item in enumerate(L):
 #       item = L[i]

    meansum = 0.0
    wtsum = 0.0
    unwtmeansum = 0.0  # use for stddev
    cnt = 0
#    maxerr = 0.0 
    errsum = 0.0 
    for i, run in enumerate(runlist):
        if run > (segment_start - 1): # include segment_start
            if run < (segment_end): # segment_end is start of next segment
                val = datalist[i]
                err = errorlist[i] 
                
                unwtmeansum += val   # get unwtmean in first pass for stddev
                errsum += err
   
                if err == 0:
                    continue
                weight = 1/(err*err)
                meansum += val*weight
                wtsum += weight
                
 #               if err > maxerr: 
 #                   maxerr = err

                cnt += 1
                print str(run) + " " + str(val) + "+-" + str(err)
                print unwtmeansum 
                print meansum
                print errsum
                

    wtmean = meansum/wtsum
    unwtmean = unwtmeansum/cnt
    avgerr = errsum/cnt

    stddevsum = 0.0
    for j, run in enumerate(runlist): # 2nd pass for stddev
        if run > (segment_start - 1): 
            if run < (segment_end): 
                val = datalist[j]
                stddevsum += (val - unwtmean)**2
    
    stddev = sqrt( stddevsum / (cnt - 1) )
   
#    wtmeanerr = 1/sqrt(wtsum) # --> This is only for uncorrelated measurements
    #return wtmean, wtmeanerr

#    return wtmean, maxerr, stddev
    return wtmean, avgerr, stddev

def find_segment_average_by_parm(dataset, parameter, tube, segment_start, segment_end):
    print '--------'
    print parameter
    print '--------'
    runlist = dataset[tube]['Run']
    datalist = dataset[tube][parameter]
    errorname = parameter + 'Err' 
    errorlist = dataset[tube][errorname]
    wtmean, maxerr, stddev = find_segment_average(runlist, datalist,
                                             errorlist, segment_start, segment_end)
    print '--------'
    print wtmean
    print '--------'
    return wtmean, maxerr, stddev

def print_segment_averages(runlist, datalist, errorlist):
    segmentlist = identifySourceSegments(1)
    for i, run in enumerate(segmentlist):
        if i < len(segmentlist) - 2: # stop before last run
            startrun = segmentlist[i]
            endrun = segmentlist[i+1]
            print "Average for runs " + str(startrun) + " to " + str(endrun) + ":"
            wtmean, maxerr, stddev = find_segment_average(runlist, datalist, errorlist, 
                                                     startrun, endrun)
            print str(wtmean) + " +/- " + str(maxerr) + " (stat)  +/- " + str(stddev) + " (stddev)" 

    return 0

def print_segment_averages_by_parm(dataset, parameter, tube):
    segmentlist = identifySourceSegments(1)
    for i, run in enumerate(segmentlist):
        if i < len(segmentlist) - 2: # stop before last run
            startrun = segmentlist[i]
            endrun = segmentlist[i+1]
            print "Average " + parameter + " for runs " + str(startrun) + " to " + str(endrun) + ":"
            wtmean, maxerr, stddev = find_segment_average_by_parm(dataset, parameter, tube,
                                                             startrun, endrun)
            print str(wtmean) + " +/- " + str(maxerr) + " (stat)  +/- " + str(stddev) + " (stddev)" 

    return 0
    
def print_segment_averages_by_run(dataset, tube, parmlist):
    segmentlist = identifySourceSegments(1)
    for i, run in enumerate(segmentlist):
        if i < len(segmentlist) - 2: # stop before last run
            startrun = segmentlist[i]
            endrun = segmentlist[i+1]
            print "Startrun: " + str(startrun) + ", Endrun: " + str(endrun)
            for parm in parmlist:
                wtmean, maxerr, stddev = find_segment_average_by_parm(dataset, parm, tube,
                                                                      startrun, endrun)
                print parm + ": " + str(wtmean) + " +/- " + str(maxerr) + " (stat)  +/- " + str(stddev) + " (stddev)" 
            
    return 0 
    

if __name__ == "__main__":
    constrainbool = 1
#    constrainbool = int(sys.argv[1])
    savebool = int(sys.argv[1])
    plotbool = int(sys.argv[2])

    plt.ion() #turn on interactive mode
    
    rcParams['figure.figsize'] = 10, 10     #Set default fig size

#    imagedir = '/data4/saslutsky/PulserComp/images_10_02_2014'
    imagedir = '/data4/saslutsky/PulserComp/images_10_22_2014_allruns/'
#    imagedir = '/data4/saslutsky/PulserComp/images_11_24_2014_10rangemin/'
#    imagedir = '/data4/saslutsky/PulserComp/images_12_15_2014_fixBeta'
#    imagedir = '/data4/saslutsky/PulserComp/images_01_14_2015_scaleE'
#    imagedir = '/data4/saslutsky/PulserComp/images_01_23_2015_widerrange'
#    imagedir = '/data4/saslutsky/PulserComp/images_04_08_2015_21596_21605'
        
#    if ~constrainbool:
    filename = 'FitResults.txt'
#    if constrainbool:
#        filename = 'FitResults_PE_PMT.txt'
#        filename = 'FitResults_constrained.txt'
        
    path = imagedir + "/" +filename
    
    outputfilename = 'PlotsTogether.pdf'
    outputpath = imagedir + "/" +  outputfilename

    # import data. 
    # TODO: improve this to read the header string from the file
    if ~constrainbool:
        data = np.genfromtxt(path, skip_header=1, 
                             delimiter = "\t", 
                             names = ['Run','Channel', 'Wavelength','p0',
                                      'p0Err','p1','p1Err','p2','p2Err', 
                                      'Chi2'])
    if  constrainbool:
        data = np.genfromtxt(path, skip_header=1, 
                             delimiter = "\t", 
                             names = ['Run','Channel', 'Wavelength','p0',
                                      'p0Err','p1','p1Err','p2','p2Err',
                                      'p3', 'p3Err', 'Chi2', ])
    fitpars = ['p0', 'p1', 'p2']
    if constrainbool:
        fitpars = ['p0', 'p1', 'p2', 'p3']

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

    if plotbool:
    #fig = plt.figure()
        fig0, (ax0, ax1) = plt.subplots(nrows=2)
        fig1, (ax2, ax3) = plt.subplots(nrows=2)
        fig2, (ax4, ax5) = plt.subplots(nrows=2)
        figChi, (ax6, ax7) = plt.subplots(nrows=2)
        figRatio, ax8 = plt.subplots()
        figures = [fig0, fig1, fig2, figChi, figRatio]
        
        plt.rc('axes', color_cycle=['r', 'g', 'b', 'y'])
    #plt.rc('axes', color_cycle=['c', 'm', 'y', 'k']) # looks terrible
        
        ax0.set_title("p0   405nm")
        ax1.set_title("p0   465nm")    
        ax2.set_title("p1   405nm")
        ax3.set_title("p1   465nm")    
        ax4.set_title("p2   405nm")
        ax5.set_title("p2   465nm")   
        ax6.set_title("Chi^2/nu   405nm")
        ax7.set_title("Chi^2/nu   465nm")
        ax8.set_title("Gain Ratio (A_405nm/A_465nm)")
        
        ax1.set_xlabel('Run Number')
        ax3.set_xlabel('Run Number')
        ax5.set_xlabel('Run Number')
        ax7.set_xlabel('Run Number')
        ax8.set_xlabel('Run Number')
        
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
                
#            if i != 2:
#                continue
 
       # 405nm
#        print_segment_averages(data_cut_405[i]['Run'], data_cut_405[i]['p0'], 
 #                              data_cut_405[i]['p0Err'])
 
# Testing
#        print "Tube: " + str(i)
#        print "405 nm"
#        print_segment_averages_by_run(data_cut_405, i, fitpars) 
#        print "465 nm" 
#        print_segment_averages_by_run(data_cut_465, i, fitpars) 

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
        #ax8.errorbar(data_cut_405[i]['Run'], ratioGain[i], yerr=ratioGainErr[i],
        #         linestyle='None', marker='o', markersize=marks, label=_chan)
   
        leg0 = ax0.legend(title = "PMT")
        leg1 = ax1.legend(title = "PMT")
        leg2 = ax2.legend(title = "PMT")
        leg3 = ax3.legend(title = "PMT")
        leg4 = ax4.legend(title = "PMT")
        leg5 = ax5.legend(title = "PMT")
        leg6 = ax6.legend(title = "PMT")
        leg7 = ax7.legend(title = "PMT")
        leg8 = ax8.legend(title = "PMT")

        ax0.set_xlim([21000, 23500])
        ax1.set_xlim([21000, 23500])
        ax2.set_xlim([21000, 23500])
        ax3.set_xlim([21000, 23500])
        ax4.set_xlim([21000, 23500])
        ax5.set_xlim([21000, 23500])
        ax6.set_xlim([21000, 23500])
        ax7.set_xlim([21000, 23500])
        ax8.set_xlim([21000, 23500])

        ax0.xaxis.set_major_formatter(FormatStrFormatter('%d'))
        ax1.xaxis.set_major_formatter(FormatStrFormatter('%d'))
        ax2.xaxis.set_major_formatter(FormatStrFormatter('%d'))
        ax3.xaxis.set_major_formatter(FormatStrFormatter('%d'))
        ax4.xaxis.set_major_formatter(FormatStrFormatter('%d'))
        ax5.xaxis.set_major_formatter(FormatStrFormatter('%d'))
        ax6.xaxis.set_major_formatter(FormatStrFormatter('%d'))
        ax7.xaxis.set_major_formatter(FormatStrFormatter('%d'))
        ax8.xaxis.set_major_formatter(FormatStrFormatter('%d'))

#    makeVLinesforSegments(ax0)
#    makeVLinesforSegments(ax2)
#    makeVLinesforSegments(ax4)
#    makeVLinesforSegments(ax6)


        ax0.set_ylim([-1000, 1000])
        ax1.set_ylim([-1000., 1000.])
        ax2.set_ylim([-10, 10])
        ax3.set_ylim([-10, 10.0])
        ax4.set_ylim([-1, 1])
        ax5.set_ylim([-1, 1])
        ax6.set_ylim([-0.06, 0.02])
        ax7.set_ylim([0.0, 25.0])
        ax8.set_ylim([0, 8])

        if savebool:
            outputfile = PdfPages(outputpath)
            for f in range(0, len(figures)):
                outputfile.savefig(figures[f])
            outputfile.close()

        plt.show(block=True)     #block=True keeps the plot window open when in interactive mode 


    # write to file the averages
    if savebool:
        outfilename = imagedir + "/AverageLEDLinearityParms.txt"
        outfile = open(outfilename, "w")
        segmentlist = identifySourceSegments(1)
        print segmentlist
        datasets = [data_cut_405, data_cut_465]
        for s in range(2):
            wave = ""
            if s == 0:
                wave = "405"
            if s == 1: 
                wave = "465"
            for parm in fitpars: # number of fitted parameters
                for tube in range(8): #tubes
                    for i, run in enumerate(segmentlist): #run segments
                        if i < len(segmentlist) - 2: # stop before last run
                            startrun = segmentlist[i]
                            endrun = segmentlist[i+1]                        
                            wtmean, maxerr, stddev = find_segment_average_by_parm(datasets[s], parm, tube,
                                                                                  startrun, endrun)
#                            print str(startrun) + " to  " + str(endrun) + ": tube " + str(tube) + ", " + parm
#                            print wtmean
                            writestring = str(startrun) + "\t" + str(endrun) + "\t"
                            writestring += str(tube) + "\t" + wave + "\t" + parm + "\t"
                            writestring += str(wtmean) + "\t" + str(maxerr) + "\t" + str(stddev)
                            writestring += "\n"
                            outfile.write(writestring)
                    
        outfile.close()


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
