# Ties together Beta Endpoint Analysis with LED/PD and Bi Pulser analysis
# Makes plots comparing the 3 calibrations

import sys
sys.path.append('/home/saslutsky/UCNA/UCNAReplay_052214/UCNA/Scripts/LEDPulser/PD_LED_Analysis')
sys.path.append('/home/saslutsky/UCNA/UCNAReplay_052214/UCNA/Scripts/LEDPulser/BetaSpectrumFitting')

#from led_bi_compare import fetch_bi_data
from led_bi_compare import make_bi_arrays, normalize_array_to_avg, process_LED_PD_corr, find_array_overlap
#from led_bi_compare import process_LED_PD_corr
from BetaEnergyScaleBatch import makeBetaEndpointPlot, truncate_runs_before
from pd_led_pmt_gain import getLEDDataforTube
from RunLogUtilities import identifySourceSegments, makeVLinesforSegments

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

def makeTripleComparePlot(tubeIn):
#    ax_bi = plt.axes()
#    ax_beta = plt.axes()
#    ax_LED = plt.axes()
#    fig = plt.figure()

    fig, (ax_bi, ax_led_pmt) = plt.subplots(nrows=2)

    ax_bi = make_bi_ax(ax_bi, tubeIn, 1)
    ax_beta = makeBetaEndpointPlot(tubeIn)

    return fig

def findAVG(arrayIn):
    avg = sum(arrayIn)/float(len(arrayIn))
    return avg

def getBetaData(tubeIn, filename):
    print "Making Beta Axis"
    beta_runlist, beta_data =  makeBetaEndpointPlot(tubeIn, 1, 1, 1, filename, 0)
    print "Beta Avg: " + str(findAVG(beta_data))
    
    return beta_runlist, beta_data

def getBiData(tubeIn):
    # Get Bi Runs and Data
    print "Making Bi arrays"
    #bi_arrays = make_bi_arrays(tubeIn, 1)  # 1 = normalize
    bi_arrays = make_bi_arrays(tubeIn, 0)  # 0 = don't normalize
    bi_data    = bi_arrays[0]
    bi_runlist = bi_arrays[2]
    bi_runlist, bi_data = truncate_runs_before(bi_runlist, bi_data, 20976.0)
    bi_data = normalize_array_to_avg(bi_data)
    print "Bi Avg: " + str(findAVG(bi_data))
    
    return bi_runlist, bi_data

def getLEDData(tubeIn, make_PD_LED_corr = 1):
    # Get LED Runs and Data
    print "Getting LED Data"
    LED_arrays  = getLEDDataforTube(tubeIn)
    if make_PD_LED_corr:
        LED_arrays = process_LED_PD_corr(LED_arrays)
    LED_runlist = LED_arrays['Run']
    LED_data    = LED_arrays['Mu']
    LED_runlist, LED_data = truncate_runs_before(LED_runlist, LED_data, 20976.0)
    LED_data = normalize_array_to_avg(LED_data)
    print "LED Avg: " + str(findAVG(LED_data))
    
    return LED_runlist, LED_data

def makeTripleComparePlot2(tubeIn, zoombool = 0, make_PD_LED_corr = 1):
    fig = plt.figure(figsize = (14, 12))
    markSize = 5
    
    # Get Beta Runs and Data
#    filename = "/home/saslutsky/UCNA/UCNAReplay_052214/UCNA/Scripts/LEDPulser/BetaSpectrumFitting/BetaLinEndpoints_FixedFitRange.txt"
    filename = "/home/saslutsky/UCNA/UCNAReplay_052214/UCNA/Scripts/LEDPulser/BetaSpectrumFitting/BetaKurieEndpoints_clean.txt"

    beta_runlist, beta_data = getBetaData(tubeIn, filename)
    bi_runlist, bi_data = getBiData(tubeIn)
    LED_runlist, LED_data = getLEDData(tubeIn, make_PD_LED_corr)

    print "Plotting Bi"
    plt.errorbar(bi_runlist, bi_data, color = 'k',
                 linestyle='None', marker='o', markersize=markSize, label = "Bi")

    print "Plotting Beta"
    plt.errorbar(beta_runlist, beta_data, color = 'b', 
                 linestyle = 'None', marker = 'o', markersize=markSize, label = "Beta")

    print "Plotting LED"
    plt.errorbar(LED_runlist, LED_data, color = 'g', 
                 linestyle = 'None', marker = 'o', markersize=markSize, label = "LED")

    plt.legend( loc = 2, title = ("W" if tubeIn/4 else "E") + str(tubeIn%4) )
    
    ax = plt.gca()
    if zoombool:
        ax.set_ylim(0.5, 1.5)
    
#    makeVLinesforSegments(ax, 20976)  # add vertical lines for each run cycle
    
    return fig

def makeTripleCompareDiffs(tubeIn, zoombool = 0, make_PD_LED_corr = 1):
    fig = plt.figure(figsize = (14, 12))
    markSize = 5
    filename = "/home/saslutsky/UCNA/UCNAReplay_052214/UCNA/Scripts/LEDPulser/BetaSpectrumFitting/BetaKurieEndpoints_clean.txt"
    
    beta_runlist, beta_data = getBetaData(tubeIn, filename)
    bi_runlist, bi_data = getBiData(tubeIn)
    LED_runlist, LED_data = getLEDData(tubeIn, make_PD_LED_corr)
   
    overlap_runlist = find_array_overlap(beta_runlist, LED_runlist) # LED list is smallest
    overlap_runlist = find_array_overlap(overlap_runlist, bi_runlist)

    beta_data = sync_arrays(overlap_runlist, beta_runlist, beta_data)
    bi_data   = sync_arrays(overlap_runlist, bi_runlist, bi_data)
    LED_data  = sync_arrays(overlap_runlist, LED_runlist, LED_data)

    offset = 0.25
    beta_led_diff = [beta - LED for beta, LED in zip(beta_data, LED_data)]
    beta_bi_diff  = [beta - bi + offset for beta, bi in zip(beta_data, bi_data)]
    led_bi_diff   = [led  - bi - offset for led,  bi in zip(LED_data, bi_data)]
    
    # plots of diff vs runlist
    print "Plotting Beta-LED diff"
    plt.errorbar(overlap_runlist, beta_led_diff, color = 'k', 
                 linestyle = 'None', marker = 'o', markersize = markSize, label = "Beta minus LED")
    print "Plotting Beta-Bi diff"
    plt.errorbar(overlap_runlist, beta_bi_diff, color = 'b', 
                 linestyle = 'None', marker = 'o', markersize = markSize, label = "Beta minus Bi + " + str(offset) )
    print "Plotting LED-Bi diff"
    plt.errorbar(overlap_runlist, led_bi_diff, color = 'r', 
                 linestyle = 'None', marker = 'o', markersize = markSize, label = "LED  minus Bi - " + str(offset) )

    plt.legend( loc = 2, title = ("W" if tubeIn/4 else "E") + str(tubeIn%4) )

    ax = plt.gca() 
    if zoombool:
        ax.set_ylim(0.5, 1.5)
        
#    return beta_led_diff, beta_bi_diff
    return fig

def sync_arrays(overlap_runlist, runlist, datalist):  #remove data from list to match a given runlist
    if len(runlist) != len(datalist):
        print "Bad input"
        return 0
    outdata_list = list()
    for i in range(len(overlap_runlist)):
        data_index = runlist.index(overlap_runlist[i])
        outdata_list.append(datalist[data_index])
    print len(outdata_list)
    print len (overlap_runlist)

    return outdata_list

def LED_PD_corr_test(tubeIn, norm = 1):
    fig = plt.figure()
    LED_arrays  = getLEDDataforTube(tubeIn)
    LED_runlist = LED_arrays['Run']
    LED_data    = LED_arrays['Mu']
    LED_data_1 = LED_data
    if norm:
        LED_data_1 = normalize_array_to_avg(LED_data)
    LED_runlist_1 = LED_runlist
    LED_runlist_1, LED_data_1 = truncate_runs_before(LED_runlist_1, LED_data_1, 20976.0)
    plt.errorbar(LED_runlist_1, LED_data_1, color = 'g', 
                 linestyle = 'None', marker = 'o', markersize=6)
    
    LED_arrays = process_LED_PD_corr(LED_arrays)
    LED_runlist_2 = LED_arrays['Run']
    LED_data_2 = LED_arrays['Mu']
    LED_runlist_2, LED_data_2 = truncate_runs_before(LED_runlist_2, LED_data_2, 20976.0)
    if norm:
        LED_data_2 = normalize_array_to_avg(LED_data_2)
    diff = [LED1 - LED2 for LED1, LED2 in zip(LED_data_1, LED_data_2)]
#    return LED_data_1, LED_data_2, diff

    plt.errorbar(LED_runlist_2, LED_data_2, color = 'r', 
                 linestyle = 'None', marker = 'o', markersize=6)

    plt.show()

    return 0

def makeAllTriplePlots(savebool = 0):
    outputfile = PdfPages("./Figures/TripleComparePlots_Kurie.pdf")
    outputfilezoom = PdfPages("./Figures/TripleComparePlots_zoom_Kurie.pdf")
    for i in range(0, 8):
        print "--------------------------"
        print "Making Plot " + str(i)
        print "--------------------------"
        fig = makeTripleComparePlot2(i, 0)
        figzoom = makeTripleComparePlot2(i, 1)
        if savebool:
            outputfile.savefig(fig)
            outputfilezoom.savefig(figzoom)
    outputfile.close()
    outputfilezoom.close()

    return 0

def makeAllTriplePlotsDiffs(savebool = 0):
    outputfile = PdfPages("./Figures/TripleCompareDiffs_Kurie.pdf")
    outputfilezoom = PdfPages("./Figures/TripleCompareDiffs_zoom_Kurie.pdf")
    for i in range(0, 8):
        print "--------------------------"
        print "Making Plot " + str(i)
        print "--------------------------"
        fig = makeTripleCompareDiffs(i, 0)
        figzoom = makeTripleCompareDiffs(i, 1)
        if savebool:
            outputfile.savefig(fig)
            outputfilezoom.savefig(figzoom)
    outputfile.close()
    outputfilezoom.close()
    

if __name__ == "__main__":
    #makeTripleComparePlot2(0)

#    makeAllTriplePlots(1)
    makeAllTriplePlotsDiffs(1)
   #plt.show()
