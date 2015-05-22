#led_bi_compare.py
# Simon Slutsky
# 11/06/2013
#
# compare led pulser and Bi pulser time history

from pd_led_pmt_gain import * 
from read_bi_pulser import * 
from pd_led_gain import * 
from LEDutils import *

def fetch_bi_data(keystringIn, valerr):
    # from read_bi_pulser.py
    bi_array = return_key_vals(keystringIn, valerr)
    
    return bi_array

def fetch_LED_PD_response():
    LED_PD_data = read_PD_LED_response()  # Get PD response to LED, see pd_led_gain.py
    LED_PD_means = LED_PD_data['Mu']

    return LED_PD_means

def fetch_LED_PD_errs():
    LED_PD_data = read_PD_LED_response()  # Get PD response to LED, see pd_led_gain.py
    LED_PD_errs = LED_PD_data['MuErr']

    return LED_PD_errs

def fetch_LED_PD_keystring(keystring):
    LED_PD_data = read_PD_LED_response()  # Get PD response to LED, see pd_led_gain.py
    LED_PD_keystring = LED_PD_data[keystring]

    return LED_PD_keystring

def correct_LED_with_PD(LED_PMT_array, LED_PD_array):
    # divide PMT response by PD response
    # works for errors too, since dominated by PMT errors
    # i.e., d(PMT/PD) = (dPMT)/PD
  
    # Runs should be sync'd already by the time this is called, but just in case
    if len(LED_PMT_array['Run']) != len(LED_PD_array['Run']):
        print "---------------------------------------------------"
        print "PMT data and PD data have different numbers of runs"
        print len(LED_PMT_array['Run'])
        print len(LED_PD_mu['Run'])
        print "---------------------------------------------------"
        return LED_PMT_array

    numbruns = len(LED_PMT_array)

    # check rows are identical 
    for i in arange(numbruns):
        #print str(LED_PMT_array[i]['Run']) + " " + str(LED_PD_array[i]['Run']) 
        if float(LED_PMT_array[i]['Run']) !=  float(LED_PD_array[i]['Run']):
            print "Runs aren't synced"
            return LED_PMT_array
    
    print "Correcting using PD values with runlist of " + str(len(LED_PMT_array['Run'])) + " runs"

    for i in arange(numbruns):
        if LED_PD_array[i]['Mu'] != 0.0:       # one last check
            LED_PMT_array[i]['Mu'] = LED_PMT_array[i]['Mu']/LED_PD_array[i]['Mu']
            LED_PMT_array[i]['MuErr'] = LED_PMT_array[i]['MuErr']/LED_PD_array[i]['Mu']
        else: 
            print "PD data = 0 for run " + str(LED_PMT_array[i]['Run'])
               
    return LED_PMT_array
 
def check_runs(totalBool = False, nonBool = False): # for testing only
    bi_runlist = GetBiRuns()
    led_pmt_data = getLEDDataforTube(0)
    led_pmt_runlist = led_pmt_data['Run']
    led_pd_data = read_PD_LED_response()
    led_pd_runlist = led_pd_data['Run']

#    overlap_runlist = list()
#    for birun in bi_runlist:
#        for ledrun in led_runlist:
#            if float(birun) == ledrun:
#                overlap_runlist.append(ledrun)
#    print len(overlap_runlist)
    tmp = find_array_overlap(led_pmt_runlist, led_pd_runlist, totalBool, nonBool)
    
    print "PD - LED:"
    print len(led_pd_data)
    print "PMT - LED:"
    print len(led_pmt_data)
    print "PMT - Bi:"
    print len(bi_runlist)

    return tmp

def find_array_overlap(list_1, list_2, total = False, non = False):
    overlap_list = list()
    total_list = list()
    non_list_1 = list()
    non_list_2 = list()
    for i in list_1:
        found = False
        for j in list_2:
            if float(i) == float(j):
                overlap_list.append(float(i))
                total_list.append("Y")
                found = True
        if found == False:
            total_list.append("N")
            non_list_2.append(i)
    for j in list_2:
        found = False
        for i in list_1:
            if float(i) == float(j):
                found = True
        if found == False: 
            non_list_1.append(j)
    non_list = zip(non_list_1, non_list_2)

    print str(len(overlap_list)) + " overlapping runs"
    if total == False:
        if non == False:
            return overlap_list
        if non == True:
            return non_list
    else:
        return total_list

def find_good_runs(dataset):
    # scan the list for runs where Peak Position ('Mu') > 0
    good_runlist = list()
    for row in dataset:
        if float(row['Mu']) > 0.0:
            good_runlist.append(row['Run'])
    
    return good_runlist
    
def sync_list(in_data, in_runlist):
    cond = [] # create boolean list indicating good runs
    for row in in_data:
        cond.append(row['Run'] in in_runlist)
    indices = np.where(cond)
    # create truncated list of only good runs
    truncated_list = in_data[indices]
    
    return truncated_list

def make_bi_ax(ax, tubeIn, normBool=1):
    markSize = 4
    bi_runlist = GetBiRuns()
    bi_data = fetch_bi_data("center", "val")
    bi_errs = fetch_bi_data("center", "err")

    bi_plot_data = bi_data[tubeIn]
    bi_plot_errs = bi_errs[tubeIn]

    if normBool:   # plot normalized PMT responses
        bi_avg_err = get_average_w_error(bi_plot_data, bi_plot_errs)
        bi_avg = bi_avg_err.value
        bi_plot_data = normalize_array(bi_plot_data, bi_avg)
        bi_plot_errs = normalize_array(bi_plot_errs, bi_avg)

    Biavg = get_avg(bi_plot_data)
    if normBool: 
        ax.set_ylabel("PMT (Normalized)")
 #       ax.set_ylim(Biavg - 0.2, Biavg + 0.2)
        ax.set_ylim(0.5, 1.3)
    else:
        ax.set_ylim(Biavg-500, Biavg+500)
        ax.set_ylabel("PMT ADC")

    ax.errorbar(bi_runlist, bi_plot_data, yerr = bi_plot_errs, 
                   linestyle='None', marker='o', markersize=markSize)

    return ax

def make_led_pmt_diff_array(tubeIn, div):
    led_pmt_data = getLEDDataforTube(tubeIn)
    led_pmt_data = process_LED_PD_corr(led_pmt_data)
    led_pmt_plot_data = led_pmt_data['Mu'] 
    led_pmt_plot_errs = led_pmt_data['MuErr']
    
    bi_runlist = GetBiRuns()
    bi_data = fetch_bi_data("center", "val")
    bi_errs = fetch_bi_data("center", "err")
    bi_plot_data = bi_data[tubeIn]
    bi_plot_errs = bi_errs[tubeIn]

    #Normalize if subtracting
    if div == 0:
        bi_avg_err = get_average_w_error(bi_plot_data, bi_plot_errs)
        bi_avg = bi_avg_err.value
        bi_plot_data = normalize_array(bi_plot_data, bi_avg)
        led_avg_err = get_average_w_error(led_pmt_plot_data, led_pmt_plot_errs)
        led_avg = led_avg_err.value
        led_pmt_plot_data = normalize_array(led_pmt_plot_data, led_avg)

    print len(bi_runlist) 
    print len(bi_plot_data)

    if len(bi_runlist) != len(bi_plot_data):
        print "Bismuth Run list and data don't match, exiting"
        return -1

    diff_list = list()
    for j in arange(len(bi_runlist)):
        hit = 0 
        for k in arange(len(led_pmt_data)):
            if float(led_pmt_data[k]['Run']) == float(bi_runlist[j]): # compare run No.
                #only main array led_pmt_data has run info
                #Use main array to identify which row of normalized array to use.
                if div == 0:
                    led_pmt_diff = led_pmt_plot_data[k] - bi_plot_data[j]
                if div == 1: 
                    if bi_plot_data[j] > 0.0:
                        led_pmt_diff = led_pmt_plot_data[k]/bi_plot_data[j]
                    else:
                        led_pmt_diff = 0.0
                hit = 1
                diff_list.append(led_pmt_diff)
                break
                
        if hit == 0:
            diff_list.append(0.0)
                           
    return diff_list

def make_led_pmt_diff_plot(tubeIn, div = 0):
    markSize = 4

    diff_array = make_led_pmt_diff_array(tubeIn, div)
    bi_runlist = GetBiRuns()    
    
    fig, ax_led_pmt_diff = plt.subplots(nrows=1)
    ax_led_pmt_diff.errorbar(bi_runlist, diff_array,
                    linestyle='None', marker='o', markersize=markSize)

    ax_led_pmt_diff.set_title("Bi Minus LED, Tube " + str(tubeIn))

    if div > 0:
        ax_led_pmt_diff.set_ylim(-0.002, 0.002)
    if div == 0:
        ax_led_pmt_diff.set_ylim(-0.2, 0.2)
    return fig
    
def make_tube_to_tube_diff_plot(tube1, tube2, axis):
    markSize = 4
    
    diffArray1 = make_led_pmt_diff_array(tube1, 0)
    diffArray2 = make_led_pmt_diff_array(tube2, 0)

    bi_runlist = GetBiRuns()

    diffArrayCompare = [d1 - d2 
                        for (d1, d2) in 
                        zip(diffArray1, diffArray2)]                        
    
    axis.errorbar(bi_runlist, diffArrayCompare, 
                  linestyle='None', marker='o', markersize=markSize)

    return axis

def make_all_tube_diffs_for_tube(tubeIn):
    #figE, (ax1, ax2, ax3, ax4) = plt.subplots(nrows = 4)
    #figW, (ax5, ax6, ax7, ax8) = plt.subplots(nrows = 4)
#    axlist = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]

#    for i in arange(len(axlist)): 
    for i in range(8):
        fig, ax = plt.subplots()
        ax = make_tube_to_tube_diff_plot(tubeIn, i, ax)
    #        axlist[i].set_ylim(-0.2, 0.2)
        ax.set_ylim(-0.2, 0.2)

    rcParams['figure.figsize'] = 10, 10

    #output = [figE, figW]

    plt.show()

    #return output
    return 0


def make_led_pmt_ax(ax_led_pmt, tubeIn, normBool, corrBool):
#def make_led_pmt_ax(tubeIn, normBool, corrBool): ### Testing only
    markSize = 4
    
    ### Testing only
    #fig, ax_led_pmt = plt.subplots(nrows=1)
    ###

    led_pmt_data = getLEDDataforTube(tubeIn)

    # Subject to revision if corrBool == 1
#    led_pmt_plot_data = led_pmt_data['Mu']
#    led_pmt_plot_errs = led_pmt_data['MuErr']
#    led_pmt_plot_runlist = led_pmt_data['Run']

    if corrBool:    # compensate for long term LED variations using PD response
        led_pmt_data = process_LED_PD_corr(led_pmt_data) # correction function
        
    # Define plot array now that potential corrections are made.    
    led_pmt_plot_data = led_pmt_data['Mu'] 
    led_pmt_plot_errs = led_pmt_data['MuErr']
    led_pmt_plot_runlist = led_pmt_data['Run']
    
    print "and plot array have"
    print str(len(led_pmt_plot_runlist)) + " Runs"
    print str(len(led_pmt_plot_data)) + " Data Points"
    print str(len(led_pmt_plot_errs)) + " Error Points"

    if normBool:   # plot normalized PMT responses
        led_avg_err = get_average_w_error(led_pmt_plot_data, led_pmt_plot_errs)
        led_avg = led_avg_err.value
        led_pmt_plot_data = normalize_array(led_pmt_plot_data, led_avg)
        led_pmt_plot_errs = normalize_array(led_pmt_plot_errs, led_avg)
   
    ax_led_pmt.errorbar(led_pmt_plot_runlist, led_pmt_plot_data, yerr=led_pmt_plot_errs,
                    linestyle='None', marker='o', markersize=markSize)
  
    LEDavg = get_avg(led_pmt_plot_data)

    if normBool:
        ax_led_pmt.set_ylabel("PMT (Normalized)")
        #ax_led_pmt.set_ylim(LEDavg - 0.2, LEDavg + 0.2)
        ax_led_pmt.set_ylim(0.5, 1.3)
    else:
        ax_led_pmt.set_ylim(LEDavg-300, LEDavg+300)
        ax_led_pmt.set_ylabel("PMT ADC")

    return ax_led_pmt

def normalize_array(array, normval):
    retarray = [arraycell/normval for arraycell in array]
    return retarray

def process_LED_PD_corr(led_pmt_data):
        # Find runs where PD response exists and truncate led_pmt_data to match
        led_pd_data = read_PD_LED_response()  # Get array with PD response to LED, see pd_led_gain.py
        led_pd_runlist = led_pd_data['Run']
        print "led_pd_runlist = " + str(len(led_pd_runlist)) + " runs"
#        print "led_pmt_plot_runlist = " + str(len(led_pmt_plot_runlist)) + " runs"
        
        # Remove bad runs from each run list
        led_pmt_good_runlist = find_good_runs(led_pmt_data)
        led_pd_good_runlist = find_good_runs(led_pd_data)
        print "led_pd_good_runlist = " + str(len(led_pd_good_runlist)) + " runs"
        print "led_pmt_good_runlist = " + str(len(led_pmt_good_runlist)) + " runs"
        
        # find overlap in good runs
        overlap_runlist = find_array_overlap(led_pmt_good_runlist, led_pd_good_runlist)

        # Remove bad runs from data sets
        led_pmt_data = sync_list(led_pmt_data, overlap_runlist) 
        led_pd_data = sync_list(led_pd_data, overlap_runlist)
        print "led_pmt_runlist truncated to"
        print str(len(led_pmt_data['Run'])) + " entries" 
        print "led_pd_runlist truncated to"
        print str(len(led_pd_data['Run'])) + " entries" 

        # Carry out correction using truncated lists
        led_pmt_data = correct_LED_with_PD(led_pmt_data, led_pd_data) 
        print "led_pmt_data has " + str(len(led_pmt_data)) + " runs after correction"
                                  
        return led_pmt_data


def make_plot(tubeIn, normBool=1, corrBool=1, fitBool=0):
    fig, (ax_bi, ax_led_pmt) = plt.subplots(nrows=2)

    ax_bi = make_bi_ax(ax_bi, tubeIn, normBool)
    ax_led_pmt = make_led_pmt_ax(ax_led_pmt, tubeIn, normBool, corrBool)

    ax_led_pmt.set_title("LED Pulser, Tube " + str(tubeIn))
    ax_bi.set_title("Bi Pulser, Tube " + str(tubeIn))
    ax_led_pmt.set_xlabel("Run Number")
    
#    ax_bi.set_xlim([21250, 24000])
#    ax_led_pmt.set_xlim([21250, 24000])
    ax_bi.set_xlim([20000, 24000])
    ax_led_pmt.set_xlim([20000, 24000])

#    if fitBool:
#        led_lin_fit, led_cov = np.polyfit(
#            led_data['Run'], led_plot_data, 1)
        # , full=True, cov=True) #, w = bi_plot_errs) need numpy update
#        print led_lin_fit, led_cov
#        led_lin_curve = np.poly1d(led_lin_fit)
#        ax_led.plot(led_lin_curve)

    return fig

# Broken since bi data and led data don't exist for same runs
def make_plot_corr(tubeIn):
    fig, ax = plt.subplots()

    bi_runlist = GetBiRuns()
    bi_data = fetch_bi_data("center", "val")
  
    led_data = getLEDDataforTube(tubeIn)

    markSize = 4
    ax.errorbar(bi_data[tubeIn], led_data['Mu'], 
                   linestyle='None', marker='o', markersize=markSize)

    ax.set_title("Bi Pulser:LED Pulser, Tube " + str(tubeIn))
    
    ax.set_xlabel("PMT response to Bi")
    ax.set_ylabel("PMT response to LED")
    
    return fig

#def make_plot2(tubeIn):
#    fig, ax = plt.subplots()

#    bi_runlist = GetBiRuns()
#    bi_data = fetch_bi_data("center", "val")
#    bi_errs = fetch_bi_data("center", "err")

#    led_data = getLEDDataforTube(tubeIn)
#    led_data_scaled_mu = scale_array(led_data['Mu'], 2000)

#    markSize = 4
#    ax.errorbar(bi_runlist, bi_data[tubeIn], yerr = bi_errs[tubeIn], 
#                   linestyle='None', marker='o', markersize=markSize)

#    ax.errorbar(led_data['Run'], led_data_scaled_mu, yerr=led_data['MuErr'],
#                 linestyle='None', marker='o', markersize=markSize)
  
    #ax.set_title("Bi Pulser, Tube " + str(tubeIn))
    #ax.set_title("LED Pulser, Tube " + str(tubeIn))
#    ax.set_title("Bi and LED Pulser, Tube " + str(tubeIn))

#    ax.set_xlabel("Run Number")
    
#    ax.set_xlim([21250, 24000])

#    return fig

def scale_array(arrayIn, scaler):
    arrayOut = list()
    for element in arrayIn:
        tmp = float(element) + scaler
        arrayOut.append(tmp)
    return arrayOut
    
# This is nice code but currently no runs
# have both Bi and LED pulser data so it does nothing
def make_ratio_plot(tubeIn):
    fig_ratio, ax_ratio = plt.subplots()

    bi_runlist = GetBiRuns()
    bi_data_all = fetch_bi_data("center", "val")
    bi_errs_all = fetch_bi_data("center", "err")

    bi_data = bi_data_all[tubeIn]
    bi_errs = bi_errs_all[tubeIn]

    led_all_data = getLEDDataforTube(tubeIn)
    led_runlist_ndarray = led_all_data['Run']
    led_data = led_all_data['Mu']
    led_errs = led_all_data['MuErr']
    led_runlist = list()
    for run in led_runlist_ndarray:
        led_runlist.append(run)
    
    bi_led_ratio_list = list()
    bi_led_ratio_err_list = list()

    for run in bi_runlist:
        bi_index = bi_runlist.index(run)
        bi_val = bi_data[bi_index]
        bi_err = bi_errs[bi_index]
        try:
            led_index = led_runlist.index(run)
            led_val = led_data[led_index]
            led_err = led_errs[led_index]
            if led_val:
                bi_led_ratio_list.append(bi_val/led_val)
                errorprop = find_ratio_err(bi_val, 
                                           led_val, bi_err, led_err)
                bi_led_ratio_err_list.append(errorprop)
            else: 
                bi_led_ratio_list.append(0.0)
        except ValueError:
            bi_led_ratio_list.append(0.0)
            bi_led_ratio_err_list.append(0.0)

#    ax_ratio.errorbar(bi_runlist, bi_led_ratio_list, yerr=bi_led_ratio_err_list,
#                      linestyle='None', marker='o', markersize=4)
    
    plt.show()

def find_ratio_err(num, den, num_err, den_err):
    # d(x/y) = (x/y)*sqrt[(dx/x)^2 + (dy/y)^2]
    ratio = num/den
    radical = sqrt( (num_err/num)*(num_err/num) + 
                    (den_err/den)*(den_err/den) )
    ratio_err = ratio*radical
    return ratio_err

def get_avg(arrayIn):
    sum = 0
    cnt = 0
    for element in arrayIn:
        sum += element
        if element:
            cnt += 1
    avg = 0
    if cnt: 
        avg = sum/cnt   

    return avg
    


def make_all_plots(save, corr, filename):
    if save:
        #outputfile = PdfPages("/data4/saslutsky/PulserComp/imagesLEDDebug/BiLEDcompare_PDcorrected_longrange.pdf")
        path = "/data4/saslutsky/PulserComp/imagesPMTGain/"
        filepath = path + filename + ".pdf"
        outputfile = PdfPages(filepath)

    for tube in range(8):
        print "Making Plot " + str(tube)
        figr = make_plot(tube, 1, corr, 0) #tube/norm/corr/fit
        if save:
            outputfile.savefig(figr)
    if save:
        outputfile.close()

def make_all_diff_plots(div = 0, save = 0):
    if save:
        outputfile = PdfPages("/data4/saslutsky/PulserComp/imagesLEDDebug/BiLEDcompate_diff.pdf")
    for tube in range(8):
        print "Making Diff Plot " + str(tube)
        figr = make_led_pmt_diff_plot(tube, div)
        if save:
            outputfile.savefig(figr)
    if save: 
        outputfile.close()


if __name__ == "__main__":

    #for tube in range(8):
     #   figr = make_plot(tube)
    rcParams['figure.figsize'] = 10, 10
    plt
#    make_all_plots(1, "corr")
    make_all_plots(1, 1)

    plt.show()

#    make_tube_plots()
       
