# small script to run pd_led_pmt_combinedfit.cc over a series of runs

import os
import sys

if __name__== "__main__":
    startrun = int(sys.argv[1])
    endrun = int(sys.argv[2]) + 1
    
    print "Commands:"
    for run in range(startrun, endrun):
        print "./pd_led_pmt_combinedfit_analysis "+ str(run)

    for run in range(startrun, endrun):
        print "Doing:"
        print "./pd_led_pmt_combinedfit_analysis "+ str(run)
        os.system("./pd_led_pmt_combinedfit_analysis "+ str(run))


