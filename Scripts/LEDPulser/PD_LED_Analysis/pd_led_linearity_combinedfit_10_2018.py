# pd_led_linearity_combinedfit_10_2018.py
# Author: Simon Slutsky
# Created: 10/08/2018
### Adapted (mostly simplified) from pd_led_linearity_combinedfit_09_2018.py for legibility.

# usage: python pd_led_linearity_combinedfit.py <start run> <end run>

import numpy as np
from pylab import *
from math import sqrt
import sys
import weightedaverages as wavg
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

def makeQualityCuts(data):
    cutruns = [18424, # explicit bimodal change in the middle of the run (only E?) 
               18429, 18430, 18520, 18559, 18644, 18646, 18648, 18787,  # LEDs swapped, bad cycle finding
               18456, 21517, 21825, 22278, 23159, 23163, 22547, # bad fits
               18505, 18506, 18524, 18532, 18533, # no data 
               18577, 18969, 19322, # bad cylce finding
               18746, 18751, 18760, 18761, 
               21186, 22647, 22037, # just a bad run
               21288, 21295, 21298, 21312, 21318, 
               21378,  # bad run
               21599, 21601, 21602, 
               21674, 21686, 
               21792, # refuses to fit
               21903, 21932, 21933, 
               21777, 21891, 21948, 21911, #runs including a bad cycle
               21429, 21430, 21431, 21432, 22261, # no PMT data
               20782, 20889, 20685, 20287]  # haven't/should investigate by eye
 
    bimodalranges = [range(18462, 18477), 
                     range(18539, 18544),
                     [18731], 
                   range(20988, 20994), 
                     range(20994, 21041), # no PMT data
                     range(21045, 21048),
                     range(21053, 21056), range(21070, 21114), range(21127, 21166),
                     range(21170, 21192), range(21194, 21232), range(21238, 21253),
                     range(21290, 21343), range(21353, 21377), range(21412, 21419),
                     range(21629, 21633), range(21640, 21704), range(21827, 21842),
                     range(21865, 21877), range(21993, 22003), range(22084, 22093),
                     [22111],
                     range(22123, 22163), range(22311, 22339), range(22344, 22350),
                     range(22358, 22389), range(22412, 22437), range(22441, 22480),
                     range(22483, 22511), range(22733, 22737), range(22746, 22753),
                     range(22793, 22803), range(22829, 22834), range(22873, 22896),
                     range(22947, 22953), range(22955, 22980), range(23010, 23017),
                     range(23081, 23084)]
    for bmr in bimodalranges:
        for j in bmr:
            cutruns.append(j)

    ledcullrange0 = range(0, 17544) # cull runs before 17544, the first run with linear ramp
    # not implying 17544-17646 are useful, but are not obviously cullable
    ledcullrange1 = range(17647, 18422) # cull runs 17647-18421, LEDs swapped, bad cycle finding
    ledcullrange2 = range(18651, 18661) # cull 18651-18660, no sweep data?
    ledcullrange3 = range(18684, 18711) # cull 18684-18710, no or little data
    ledcullrange4 = range(18769, 18781) # cull 18769-18780, 
    ledcullrange4a = range(18794, 18801) # cull 18794-18800, no sweep data or bad data
    ledcullrange4b = range(19002, 19024 ) # cull 19002-19023, no sweep data
    ledcullrange4c = range(19047, 19317 ) # cull 19047-19316, no sweep data

#########    # 20970 seems to be first good run, runs before it have confusingly large offset (p0)
    ledcullrange5 = range(20254, 20970) 
    
###############################
    ledcullrange6 = [22222, 22301, 22302]#, 22448, 22450] 10/02/18 22448, 22450 have fine linearity fits but some bimodal stuff may be affecting... shadowing shows up on linearity sweeps	 
    ledcullrange7 = range(22453, 22463) # cull 22453-22462 10/02/18 ditto
    ledcullrange8 = [22778]
    ledcullranges = [ledcullrange0, ledcullrange1, 
                     ledcullrange2, ledcullrange3, 
                     ledcullrange4, ledcullrange4a,
                     ledcullrange4b,ledcullrange4c,
                     ledcullrange5, ledcullrange6,
                     ledcullrange7, ledcullrange8 ]
    for lcr in ledcullranges:
        for k in lcr:
            cutruns.append(k)
  #  print cutruns
    cutsAll = True
    for i in cutruns:
        _cut = (data['Run']!=i)
        cutsAll = cutsAll & _cut
    _data = data[cutsAll]
    return _data

def checkChiSq(data, run):
    runCut = data['Run'] == run
    chiCut = data['ParName'] == "Chisq/nvars"
    cutCond = runCut & chiCut
    chi_data_cut = data[cutCond]
    retarray = chi_data_cut['ParVal'] != 0.0  # cut runs with chisq = 0
#    print retarray
    if retarray[0]:
        return run
    else:
        return 0

def makeChiCuts(data):
    cutlist = []
    runlist = data['Run']
    runlist_short = list()
    for run in runlist:
        runinshort = 0
        for run_short in runlist_short:
            if run == run_short:
                runinshort = 1
        if runinshort == 0:
            runlist_short.append(run)

    for run in runlist_short:
#        print "Checking ChiSq for run " + str(run)
        cutrun = checkChiSq(data, run)
        cutlist.append(cutrun)
        
#    _data = [row for row in data if row['Run'] not in cutlist]

    # Black magic, 
    # see http://stackoverflow.com/questions/1962980/selecting-rows-from-a-numpy-ndarray
    _data = data[np.logical_or.reduce([data['Run'] == x for x in cutlist])]

    return _data

def makeErrCuts(data, threshold=10):
    cutcond = np.where( (data['ParErr'] > 1e5) | (abs(data['ParErr']) > abs(threshold*data['ParVal'])) | abs((data['ParVal']) > 1e5) )
   # for row in data[cutcond]:
   #     print row
    newdata = np.delete(data, cutcond, 0)

    return newdata

def getLEDdata(basedir):
    # import data. 

#    basedir = "/data1/saslutsky/LEDPulser/images_05_26_2015_16way_separate_wavelength_coeff_residuals_21650_21950/

    data = np.genfromtxt(basedir + "FitResults_Combined.txt", 
                         skip_header=0, 
                         delimiter = "\t", 
                         usecols = range(0,6),
                         names = ['Run','Channel', 'ParName',
                                  'ParVal', 'ParErr', 
                                  'Rmin', 'Rmax', 'Rmin2', 'Rmax2'],
                         dtype = "int, int, S12, float, float, float" )
    return data


def takeAverages(dataarray, stripbool):
    segments = catSourceSegments() # find run segment boundaries; use for averaging data
    
    dataruns = dataarray['Run']
    datadata = dataarray['ParVal']

    avglist = list()
    stdlist = list()
    seglist = list()
    for i, segrun in enumerate(segments):
        segavg = 0
        segstd = 0
        segdatalist = list()
        for j, drun in enumerate(dataruns):
            if i == len(segments) - 1: #special case for end of array   #            if drun > segments[-1]: #special case for runs after end of array
                if drun > segrun:
                    if drun == dataruns[j-1]:
                        continue
                    segdatalist.append(datadata[j])
            elif i < len(segments) - 1:
                if drun > segrun and drun < segments[i+1]: # start with initial segment 
                    if drun == dataruns[j-1]:
                        continue
                    segdatalist.append(datadata[j])
            else:
                print "Waaah? That shouldn't happen."    
         
        print "SEGRUN: " + str(segrun)
        print "LIST: " + str(segdatalist)
        #sometimes it cuts out entire arrays?
        if stripbool:
            stripsegdatalist = stripOutlier(segdatalist) # see below
            if len(stripsegdatalist) != 0:
                segavg = np.mean(stripsegdatalist)
                segstd = np.std(stripsegdatalist)
            else:
                segavg = 0
                segstd = 0
        else:
            if len(segdatalist) != 0:
                segavg = np.mean(segdatalist)
                segstd = np.std(segdatalist)
            else:
                segavg = 0
                segstd = 0
        print "AVERAGE: " + str(segavg)
        avglist.append(segavg)
        stdlist.append(segstd)

    return avglist, stdlist, segments

def stripOutlier(array, thresh1 = 5, thresh2 = 2):  # strip values more than thresh2*stddev from mean
    print "Stripping array " 
    print array

    astddev = np.std(array)
    amean = np.mean(array)

    # check for outliers
    candidateoutliers = [a for a in array if abs(a - amean) > thresh1*astddev]
    print candidateoutliers
    # make subarray of values that are not potentially outliers
    trimarray1 = [a for a in array if abs(a - amean) < thresh1*astddev]

    astddev_trim = np.std(trimarray1)
    amean_trim = np.mean(trimarray1)
    print "STD: " + str(astddev_trim)
    print "MEAN: " + str(amean_trim)
    

    # strip outliers far from trimmed mean
    trimarray2 = [a for a in array if abs(a - amean_trim) < thresh2*astddev_trim]
    outliers2 = [a for a in array if abs(a - amean_trim) > thresh2*astddev_trim]
    print "Cutting Outliers"
    print outliers2

    print "Returning array " 
    print trimarray2
    return trimarray2
    
if __name__ == "__main__":
    
    if len(sys.argv) < 5:
        print "\n Usage: python pd_led_linearity_combinedfit.py <showbool> <savebool> <datebool> <averagebool> <stripbool>"
        print " Please set flags\n"
        sys.exit()

    showbool = int(sys.argv[1]) # boolean to control plotting
    savebool = int(sys.argv[2]) # boolean to control saving
    datebool = int(sys.argv[3]) # boolean to toggle run#/date for x-axis
    avgebool = int(sys.argv[4]) # boolean to average data points over a run segment
    stripbool = int(sys.argv[5]) # boolean to remove outliers from averages

    rcParams['figure.figsize'] = 10, 10     #Set default fig size
    plt.rc('axes', color_cycle=['r', 'g', 'b', 'y'])

    basedir = "/data1/saslutsky/LEDPulser/images_10_01_2018_SourcePeriod23_22767_22791/"

    data = getLEDdata(basedir)

    runlist = data['Run']
#cut bad runs
    print "Making Quality Cuts"
    data = makeQualityCuts(data)   # cut runs found bad by eye
    print "Making Chi^2 cuts"
    data = makeChiCuts(data)       # cut runs with chisq == 0.0
    data = makeErrCuts(data)

    # open output file
    outputfile = PdfPages(basedir + "FitResultsCombined_" + 
                          str(min(runlist)) + "_" +
                          str(max(runlist)) + ".pdf")
    p1file = PdfPages(basedir + "FitResultsCombined_" + 
                          str(min(runlist)) + "_" +
                          str(max(runlist)) + "_p1.pdf")
    

    npars = 4   #p0, p1, p2, nlambda for each tube
    marks = 8
    parnames = ["p0", "p1", "p2", "nlambda"]
    pdparnames = ["PDq1", "PDq2", "PDq3"]

    axes = list()
    figures = list()
    p1Figures = list()
   
    pd_data_cut = list()
    print "making pdpar cuts"
    for k in range(0, len(pdparnames)):
        _pd_data_cut = data[data['ParName'] == pdparnames[k]]
        pd_data_cut.append(_pd_data_cut)

    tmpFig1, (ax0, ax1, ax2) = plt.subplots(nrows = 3, sharex = True, sharey = False)
    pdaxes = [ax0, ax1, ax2]
#    figures.append(tmpFig1)

    for k in range(0, len(pdparnames)): 
        if datebool:
            pdtimelist = getTimeForRunlist(pd_data_cut[k]['Run'])
            pdaxes[k].xaxis.set_major_formatter(dates.DateFormatter("%m/%d/%y"))
        
        if not datebool:
            pdaxes[k].errorbar(pd_data_cut[k]['Run'], 
                               pd_data_cut[k]['ParVal'], 
                               yerr=pd_data_cut[k]['ParErr'], 
                               linestyle='None', marker='o',
                               markersize=marks, color = 'Black')
            xmin, xmax = min(pd_data_cut[k]['Run']), max(pd_data_cut[k]['Run']) 
            pdaxes[0].set_xlim([xmin - 0.5, xmax + 0.5])

        if datebool:
            pdaxes[k].errorbar(pdtimelist, 
                               pd_data_cut[k]['ParVal'], 
                               yerr=pd_data_cut[k]['ParErr'], 
                               linestyle='None', marker='o',
                               markersize=marks, color = 'Black')
        
        pdaxes[k].set_ylabel("q" + str(k+1))
        pdaxes[0].set_title("PD nonlinearity parameters")
        if not datebool:
            pdaxes[2].set_xlabel('Run Number')
        if datebool:
            pdaxes[2].set_xlabel('Time')

    txtoutfile = open('NonLinAverages.txt', 'a')
    wtdaveragefile = 
        
    for i in range(0,8): # 8 channels 
        data_cut = list()
        avgdata = list() # one entry for each parval
        stddata = list()
        for j in range (0, npars): 
            # Select the channel and parameter, will plot vs run #
            cutChannel = data['Channel'] == int(i)
            cutParm = data['ParName'] == parnames[j]
            #cutErr = data['ParErr'] > abs(data['ParVal']) # this wrongly cuts values near 0
            cutCond = cutChannel & cutParm
            _data_cut = data[cutCond]
            if avgebool: # average the data now if requested
                avglist, stdlist, seglist = takeAverages(_data_cut, stripbool)
                avgdata.append(avglist)
                stddata.append(stdlist)
                segdata.append(seglist)
   
   	    data_cut.append(_data_cut) # npars arrays of ParName (for all runs)    return _data_cut
            # write to file
            if savebool and avgebool:
                if parnames[j] != "nlambda":
                    for k, segment in enumerate(seglist):
                        txtoutfile.write(str(segment) + "\t")  # segment start run number
                        txtoutfile.write(str(i) + "\t")        # tube
                        txtoutfile.write(parnames[j] + "\t")   # parameter name
                        txtoutfile.write(str(avglist[k])+ "\t")    # average parameter
                        txtoutfile.write(str(stdlist[k]) + "\n")  # average error
                        
        for avg in avgdata:
            print "---"
            print avg


        # Print Weighted Averages to screen
	print "---------"
        print "Weighted average for p2:"
        print wavg.takeWeightedAverage(data_cut[2]['ParVal'], data_cut[2]['ParErr'])
	print "---------"

        ## Prepare Plots

        tmpFig0, (tmpAx0, tmpAx1, tmpAx2, tmpAx3) = plt.subplots(nrows=4, sharex = True, sharey=False)
        p1Fig, p1Ax = plt.subplots(nrows = 1)
        axes = [tmpAx0, tmpAx1, tmpAx2, tmpAx3]
        figures.append(tmpFig0)
        p1Figures.append(p1Fig)
        _chan = ""
        if i < 4:
            _chan = "E" 
            _chan = _chan + str(i)
        if i >3 and i < 8:
            _chan = "W" 
            _chan = _chan + str(i-4)
        if i == 10:
            _chan = "PD" 
        for j in range(0, npars): 
            if not datebool:
                axes[j].errorbar(data_cut[j]['Run'], 
                                 data_cut[j]['ParVal'], 
                                 yerr=data_cut[j]['ParErr'], 
                                 linestyle='None', marker='^',
                                 markersize=marks, label=_chan, color = 'Black')
                xmin, xmax = min(data_cut[j]['Run']), max(data_cut[j]['Run']) 
                axes[0].set_xlim([xmin - 0.5, xmax + 0.5])
                axes[3].set_xlabel('Run Number')
                p1Ax.set_xlim([xmin - 0.5, xmax + 0.5])
                x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
                p1Ax.xaxis.set_major_formatter(x_formatter)
                ##
              #  print avgdata[0]
                if avgebool:
                    numsegments = len(seglist)
                    print "NUM SEGMENTS " +  str(numsegments)
                    print "PARNAME " + parnames[j]
                    for segment in range(numsegments):
                        yline = avgdata[j][segment]
                        xxminx = seglist[segment]
                        if segment < numsegments - 1:
                            xxmaxx = seglist[segment+1]
                        else: 
                            xxmaxx = seglist[segment] + 50
                        axes[j].hlines(yline, xxminx, xxmaxx, linestyles='solid',color = 'blue', linewidth = 3.0 )

            if datebool:
                timelist = getTimeForRunlist(data_cut[j]['Run'])
                axes[j].xaxis.set_major_formatter(dates.DateFormatter("%m/%d/%y"))
                axes[j].errorbar(timelist,
                                 data_cut[j]['ParVal'], 
                                 yerr=data_cut[j]['ParErr'], 
                                 linestyle='None', marker='^',
                                 markersize=marks, label=_chan, color = 'Black')
                axes[j].xaxis_date()
                axes[3].set_xlabel('Time')
                mintime, maxtime = findMinMaxTime(timelist) # from RunPlotter.py
                axes[3].set_xlim(mintime, maxtime)
                if avgebool:  
                    numsegments = len(seglist)
                    for segment in range(numsegments):
                        yline = avgdata[j][segment]
                        xxminx = seglist[segment]
                        timemin = getTimeForRunlist([xxminx])
                        if segment < numsegments - 1:
                            xxmaxx = seglist[segment+1]
                        else: 
                            xxmaxx = seglist[segment] + 50
                        timemax = getTimeForRunlist([xxmaxx])
                        axes[j].hlines(yline, timemin, timemax, linestyles='solid',color = 'blue', linewidth = 3.0 )
               #         axes[0].vlines(timemin, -10000, 10000, linestyles = 'dashed', color = 'red', linewidth = 1.0)
                    axstd = np.std(avgdata[j])
                    axmean = np.mean(avgdata[j])
                    axthresh = 4
                    axes[j].set_ylim(axmean - axthresh*axstd, axmean + axthresh*axstd)

            axes[j].set_ylabel("p" + str(j))
            
        ip1Ax = 2
        if not datebool:
            p1Ax.errorbar(data_cut[ip1Ax]['Run'], 
                          data_cut[ip1Ax]['ParVal'], 
                          yerr=data_cut[ip1Ax]['ParErr'], 
                          linestyle='None', marker='^',
                          markersize=marks, label=_chan, color = 'Black')
            p1Ax.set_xlabel('Run Number')
            if avgebool:
                for segment in range(numsegments):
                    yline = avgdata[ip1Ax][segment]
                    xxminx = seglist[segment]
                    if segment < numsegments - 1:
                        xxmaxx = seglist[segment+1]
                    else: 
                        xxmaxx = seglist[segment] + 50
                    p1Ax.hlines(yline, xxminx, xxmaxx, linestyles='solid',color = 'blue', linewidth = 3.0 )
            
        if datebool:
            timelist = getTimeForRunlist(data_cut[ip1Ax]['Run'])
            p1Ax.errorbar(timelist, 
                          data_cut[ip1Ax]['ParVal'], 
                          yerr=data_cut[ip1Ax]['ParErr'], 
                          linestyle='None', marker='^',
                          markersize=marks, label=_chan, color = 'Black')
            p1Ax.xaxis_date()
            p1Ax.set_xlabel('Time')
            mintime, maxtime = findMinMaxTime(timelist) # from RunPlotter.py
            p1Ax.set_xlim(mintime, maxtime)
            if avgebool:  
                numsegments = len(seglist)
                for segment in range(numsegments):
                    yline = avgdata[ip1Ax][segment]
                    xxminx = seglist[segment]
                    timemin = getTimeForRunlist([xxminx])
                    if segment < numsegments - 1:
                        xxmaxx = seglist[segment+1]
                    else: 
                        xxmaxx = seglist[segment] + 50
                    timemax = getTimeForRunlist([xxmaxx])
                    p1Ax.hlines(yline, timemin, timemax, linestyles='solid',color = 'blue', linewidth = 3.0 )

            
        axes[3].set_ylabel("$\eta_{\lambda}$")        
        axes[0].set_title(_chan)
        

        p1Ax.set_ylabel("p2")
        p1Ax.set_title(_chan)

    if savebool: 
        for f in range (8):
            outputfile.savefig(figures[f])
            p1file.savefig(p1Figures[f])
    
    txtoutfile.close()
    outputfile.close()
    p1file.close()
    if showbool:
        plt.show(block=True)     #block=True keeps the plot window open when in interactive mode 


