# weightedaverages.py
# Author: Simon Slutsky
# Created: 10/09/2018
### Simple Script to take a weightedaverage (for pd_led_linearity_combinedfit_10_2018.py)

import numpy as np
from pylab import *
from math import sqrt
import sys

def takeWeightedAverage(valarray, errorarray):
    if len(valarray) != len(errorarray):
    	print "Incompatible array lengths"
    	return -999
    	
    _weight = 0	
    _val = 0
    _weightedval = 0
    _sumofweightedvals = 0
    _sumofweights = 0
    print "length of array = " + str(len(valarray))
    print " --- " 
    for i in range(len(valarray)):		
	if errorarray[i] == 0:
	    return 0
	_weight = 1/errorarray[i]
        _val = valarray[i]
        _weightedval = _weight*_val
        _sumofweightedvals += _weightedval
        _sumofweights += _weight
#        print "i = " + str(i)
#        print _sumofweightedvals
#        print _sumofweights

    if _sumofweights == 0:
    	return -999 
	
    return _sumofweightedvals/_sumofweights
           
if __name__ == "__main__":
    arr = (1,2,3,4,5)
    err = (0.3, 0.4, 0.23, 0.1, .99)
    print takeWeightedAverage(arr, err)
    
