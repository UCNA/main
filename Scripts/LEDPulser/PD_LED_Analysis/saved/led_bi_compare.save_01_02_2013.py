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

def correct_LED_with_PD(LED_PMT_array, LED_PMT_runlist, errBool = False):
    LED_PD_runs = fetch_LED_PD_keystring("Run")
    LED_PD_mu = fetch_LED_PD_keystring("Mu")

   # print LED_PD_mu

    # divide PMT response by PD response
    # works for errors too, since dominated by PMT errors
    # i.e., d(PMT/PD) = (dPMT)/PD
  
    # If runs are out of order, this doesn't work!
    if len(LED_PMT_array) != len(LED_PD_mu):
        print "---------------------------------------------------"
        print "PMT data and PD data have different numbers of runs"
        print len(LED_PMT_array)
        print len(LED_PD_mu)
        print "---------------------------------------------------"
        return zip(LED_PMT_array, LED_PMT_runlist)
    
    #doesn't catch zeroes
    #LED_PMT_corr = [mPMT/mPD for mPMT, mPD in 
    #                zip(LED_PMT_array, LED_PD_array)]
    
    print "Correcting with runlist of " + str(len(LED_PMT_runlist)) + " runs"

    if errBool == False: 
        LED_PMT_corr = list()
        good_run_list = list()
        print "Processing " + str(len(LED_PD_mu)) + " entries"
        for i in range(len(LED_PD_mu)):
            iPD = LED_PD_mu[i]
            iPMT = LED_PMT_array[i]
            rPD = LED_PD_runs[i]
            rPMT = LED_PMT_runlist[i]
            if float(rPD) == float(rPMT): 
                if LED_PD_mu[i] != 0:
                    LED_PMT_corr.append(iPMT/iPD)
                    good_run_list.append(rPD)
            else:
                print "---------------------------------------------------"
                print " PD and PMT runs do not correspond for "
                print " Entry " + i + "(PD run " + str(rPD) + ", PMT run " + str(rPMT) + ")"
                print "---------------------------------------------------"
                
        print "Found " + str(len(LED_PMT_corr)) + " good data points and"
        print str(len(good_run_list)) + " good runs"

        q = zip(LED_PMT_corr, good_run_list)
        ####
        return q
    
    if errBool == True: # errors may not be 0 at the same time as values
        LED_PMT_err_corr = list()
        print "Processing " + str(len(LED_PMT_runlist)) + " entries"
        for i in range(len(LED_PMT_runlist)):
            iPD = LED_PD_mu[i]
            iPMT = LED_PMT_array[i]
            cnt = 0
            if iPD == 0:
           #     print LED_PMT_runlist[i]
            #    print LED_PD_mu[i]
                LED_PMT_err_corr.append(iPMT/iPD)
        ####
        print "Found " + str(len(LED_PMT_err_corr)) + " good error points" 
        return LED_PMT_err_corr 
                
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

def sync_pmt_list(in_data, in_runlist):
    cond = [] # create boolean list indicating good runs
    for row in in_data:
        cond.append(row['Run'] in in_runlist)
    indices = np.where(cond)
    # create truncated list of only good runs
    truncated_list = in_data[indices]
    
    return truncated_list

def make_plot(tubeIn, normBool=1, corrBool=1, fitBool=0):
    markSize = 4
    fig, (ax_bi, ax_led_pmt) = plt.subplots(nrows=2)

    bi_runlist = GetBiRuns()
    bi_data = fetch_bi_data("center", "val")
    bi_errs = fetch_bi_data("center", "err")

    led_pmt_data = getLEDDataforTube(tubeIn)

    bi_plot_data = bi_data[tubeIn]
    bi_plot_errs = bi_errs[tubeIn]

    # Subject to revision if corrBool == 1
    led_pmt_plot_data = led_pmt_data['Mu']
    led_pmt_plot_errs = led_pmt_data['MuErr']
    led_pmt_plot_runlist = led_pmt_data['Run']

    if corrBool:    # compensate for long term LED variations using PD response
        # Find runs where PD response exists and truncate led_pmt_data to match
        led_pd_runlist = fetch_LED_PD_keystring('Run')
        print "led_pd_runlist = " + str(len(led_pd_runlist)) + " runs"
        print "led_pmt_plot_runlist = " + str(len(led_pmt_plot_runlist)) + " runs"
        led_pmt_plot_runlist = find_array_overlap(
            led_pmt_plot_runlist, led_pd_runlist)
        print "Correcting using PD response."
        print "led_pmt_runlist truncated to"
        print str(len(led_pmt_plot_runlist)) + " entries" 
        led_pmt_data = sync_pmt_list(led_pmt_data, led_pmt_plot_runlist) 

        # Redefine led_pmt_plot_data with truncated runs
        led_pmt_plot_data = led_pmt_data['Mu']
        led_pmt_plot_errs = led_pmt_data['MuErr']
        print "led_pmt_data truncated to"
        print str(len(led_pmt_plot_data)) + " entries"
        
        # Carry out correction using truncated lists
        # correction returns a zip of corrected data and
        # a valid runlist where PD data is not 0
        q1 = correct_LED_with_PD(
            led_pmt_plot_data, led_pmt_plot_runlist, False) # correct data
        unzip_q1 = zip(*q1)
        led_pmt_plot_data = unzip_q1[0]
        print "correct_LED_with_PD returns " + str(len(led_pmt_plot_data)) + " entries"
        led_pmt_plot_runlist = unzip_q1[1] #updated to exclude bad PD runs
        
        led_pmt_plot_errs = correct_LED_with_PD(
            led_pmt_plot_errs, led_pmt_plot_runlist, True) # correct errs
        print "correct_LED_with_PD returns " + str(len(led_pmt_plot_errs)) + " entries"

        print "Using PD to correct PMT responses"
        print str(len(led_pmt_plot_runlist)) + " Runs"
        print str(len(led_pmt_plot_data)) + " Data Points"
        print str(len(led_pmt_plot_errs)) + " Error Points"

    if normBool:   # plot normalized PMT responses
        bi_avg_err = get_average_w_error(bi_plot_data, bi_plot_errs)
        bi_avg = bi_avg_err.value
        bi_plot_data = [bi_point/bi_avg for bi_point in bi_plot_data]
        bi_plot_errs = [bi_err/bi_avg for bi_err in bi_plot_errs]

        led_avg_err = get_average_w_error(led_pmt_plot_data, led_pmt_plot_errs)
        led_avg = led_avg_err.value
        led_pmt_plot_data = [led_point/led_avg for led_point in led_pmt_plot_data]
        led_pmt_plot_errs = [led_err/led_avg for led_err in led_pmt_plot_errs]
   
#    ax_bi.errorbar(bi_runlist, bi_data[tubeIn], yerr = bi_errs[tubeIn] ,
    ax_bi.errorbar(bi_runlist, bi_plot_data, yerr = bi_plot_errs, 
                   linestyle='None', marker='o', markersize=markSize)

#    ax_led.errorbar(led_data['Run'], led_data['Mu'], yerr=led_data['MuErr'],
#    ax_led.errorbar(led_data['Run'], led_plot_data, yerr=led_plot_errs,
    print str(len(led_pmt_plot_runlist)) + " Runs"
    print str(len(led_pmt_plot_data)) + " Data Points"
    print str(len(led_pmt_plot_errs)) + " Error Points"
    ax_led_pmt.errorbar(led_pmt_plot_runlist, led_pmt_plot_data, yerr=led_pmt_plot_errs,
                    linestyle='None', marker='o', markersize=markSize)
  
    ax_bi.set_title("Bi Pulser, Tube " + str(tubeIn))
    ax_led_pmt.set_title("LED Pulser, Tube " + str(tubeIn))
    
  #  ax_bi.set_xlabel("Run Number")
    ax_led_pmt.set_xlabel("Run Number")
    
    ax_bi.set_xlim([21250, 24000])
    ax_led_pmt.set_xlim([21250, 24000])

    LEDavg = get_avg(led_pmt_plot_data)
    Biavg = get_avg(bi_plot_data)

    if normBool:
        ax_bi.set_ylabel("PMT (Normalized)")
        ax_led_pmt.set_ylabel("PMT (Normalized)")

     #   ax_bi.set_ylim(0.8,1.2)
     #   ax_led_pmt.set_ylim(0.8,1.2)
        ax_bi.set_ylim(Biavg - 0.2, Biavg + 0.2)
        ax_led_pmt.set_ylim(LEDavg - 0.2, LEDavg + 0.2)
                    
    else:
        ax_bi.set_ylim(Biavg-500, Biavg+500)
        ax_led_pmt.set_ylim(LEDavg-300, LEDavg+300)
        
        ax_bi.set_ylabel("PMT ADC")
        ax_led_pmt.set_ylabel("PMT ADC")

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

def make_plot2(tubeIn):
    fig, ax = plt.subplots()

    bi_runlist = GetBiRuns()
    bi_data = fetch_bi_data("center", "val")
    bi_errs = fetch_bi_data("center", "err")

    led_data = getLEDDataforTube(tubeIn)
    led_data_scaled_mu = scale_array(led_data['Mu'], 2000)

    markSize = 4
    ax.errorbar(bi_runlist, bi_data[tubeIn], yerr = bi_errs[tubeIn], 
                   linestyle='None', marker='o', markersize=markSize)

    ax.errorbar(led_data['Run'], led_data_scaled_mu, yerr=led_data['MuErr'],
                 linestyle='None', marker='o', markersize=markSize)
  
    #ax.set_title("Bi Pulser, Tube " + str(tubeIn))
    #ax.set_title("LED Pulser, Tube " + str(tubeIn))
    ax.set_title("Bi and LED Pulser, Tube " + str(tubeIn))

    ax.set_xlabel("Run Number")
    
    ax.set_xlim([21250, 24000])

    return fig

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
    
def make_all_plots(save, corr):
    if save:
        outputfile = PdfPages("/data4/saslutsky/PulserComp/imagesPMTGain/BiLEDcompare_PDcorrected.pdf")
    for tube in range(8):
        print "Making Plot " + str(tube)
        figr = make_plot(tube, 1, corr, 0) #tube/norm/corr/fit
        if save:
            outputfile.savefig(figr)
    if save:
        outputfile.close()

if __name__ == "__main__":

    #for tube in range(8):
     #   figr = make_plot(tube)
    rcParams['figure.figsize'] = 10, 10
    
    make_all_plots(1, "corr")

    plt.show()

#    make_tube_plots()
       
