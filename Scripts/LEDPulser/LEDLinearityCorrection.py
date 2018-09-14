# Apply a correction to a peak energy based on the PMT - PD linearity curves
# Idea is to remove quadratic term: 
# PMT = p0 + p1*PD + p2*PD*PD (1)
# so linearized PMT value PMT' = p0 + p1*PD = PMT - p2*PD*PD
# solve (1) for PD(PMT), PD signal as function of PMT, and 
# rewrite PMT' in terms of new function:
# PMT' = PMT - p2*PD(PMT)*PD(PMT) 
# Now we have a function that can be applied to any PMT signal, 
# independent of PD data.

# 11/19/2014, SS
# 1/15/2015, SS: changed correction function
###############   now fitting about PD_Beta, best beta value for pd endpoint

from numpy import genfromtxt
from math import sqrt 

def readFitParms(filename):
    linearitydata = genfromtxt(filename, delimiter = "\t", 
                                  dtype = "i8,i8,i8,i8,S5,f8,f8,f8",
                                  names = ['Run_start', 'Run_end', 'tube', 'wave',
                                           'parm', 'val', 'meanerr', 'stddev'] )

    return linearitydata

def readFitResultsTxt(filename):
    linearitydata = genfromtxt(filename, delimiter = "\t", 
                                  dtype = "i8,i8,i8,f8,f8,f8,f8,f8,f8,f8,i8",
                                  names = ['run', 'tube', 'wave',
                                           'p0', 'p0Err', 'p1', 'p1Err', 
                                           'p2', 'p2Err', 'Chi2', 'dof'] )
                               
    return linearitydata

def getConversionFactor(filename, run, tube, wave=405):
    conversiondata = genfromtxt(filename, delimiter = "\t",
                                dtype = "i8, i8, i8, f8",
                                names = ["run","tube", "wave", "val"])

    cutRun = conversiondata['run'] == run
    cutTube = conversiondata['tube'] == tube
    cutWave = conversiondata['wave'] == wave
    cutCondition = cutRun & cutTube & cutWave

    _cutdata = conversiondata[cutCondition]
    convfactor = _cutdata[val]
    return convfactor

def findGainFactor(run, tube):
    BetaEndpoint = 782.
    # Beta endpoint in ADC counts by tube and run segment
    # run < 20500
    BetaADC_below_20500 = [500., 500., 850., 500., 700., 550., 700., 1050.]
    # 20500 < run < 21250
    BetaADC_20500_21250 = [550., 400., 850., 500., 450., 600., 500., 1050.]
    # run > 21250
    BetaADC_21250_above = [650., 650., 1200., 550., 750., 700., 900., 1000.]

    gainfactor = 0
    if run < 20501:
        gainfactor = BetaADC_below_20500[tube]/BetaEndpoint
    if run > 20500 and run < 21251:
        gainfactor = BetaADC_20500_21250[tube]/BetaEndpoint   
    if run > 21250:
        gainfactor = BetaADC_21250_above[tube]/BetaEndpoint
    print "Gain factor = " + str(gainfactor)
    
    return gainfactor

def makeLEDLinCorrByRun(value, run, tube, wave, parm_filename):
    print "making correction for run, tube:"
    print str(run) + ", " + str(tube)
 
    data = readFitResultsTxt(parm_filename)    
    
    if run < min(data['run']) or run > max(data['run']):
gain        print "Run out of range"
        return -1

    cutTube = data['tube'] == tube
    cutRun = data['run'] == run
    cutWave = data['wave'] == wave
    cutCondition = cutTube & cutRun  & cutWave

    _data_cut = data[cutCondition] 
    
    if len(_data_cut) == 0:
        print "No data found, exiting" 
        return -1

    linearityparms = list()
    linearityerrs = list()
    parms = ['p0', 'p1', 'p2']
    for p in parms:
        _data_parm = _data_cut[p]
        _data_err  = _data_cut[p+'Err'] 
        linearityparms.append( float(_data_parm) )
        linearityerrs.append( float(_data_err) )
    print linearityparms

    gainfactor = findGainFactor(run, tube)
#    new_value = value - correctionFunctionQuadratic(value, linearityparms)
    new_value = value - correctionFunctionQuadratic2(value, linearityparms, gainfactor)
 
    return new_value

def find_segment_average(data, parameter, tube, wave=405, start_run=0, end_run=99999):
    print '--------'
    print parameter
    print '--------'

    cutTube = data['tube'] == tube
    cutWave = data['wave'] == wave
    cutCondition = cutTube & cutWave
    _data_cut = data[cutCondition] 
    _dcp = _data_cut[parameter]
    _dcp = [dcpi for dcpi in _dcp if abs(dcpi) < 1e4] # 1e4 should be larger than any parm val
    
    if len(_dcp) == 0:
        print "Averaging zero terms in find_segment_average"
        return 0

    avg = sum(_dcp)/float(len(_dcp))
#    print "Average = " + str(avg)
    return avg

# add run_segment functionality soon 01/16/2015 SS
def makeLEDLinCorrByAvg(value, run, tube, wave, parm_filename, start_run=0, end_run = 99999):
    print "making correction by segment average for run, tube:"
    print str(run) + ", " + str(tube)
 
    data = readFitResultsTxt(parm_filename)    
    
    if run < min(data['run']) or run > max(data['run']):
        print "Run out of range"
        return -1

    cutTube = data['tube'] == tube
    cutRun = data['run'] == run
    cutWave = data['wave'] == wave
    cutCondition = cutTube & cutRun  & cutWave

    _data_cut = data[cutCondition] 
    
    if len(_data_cut) == 0:
        print "No data found, exiting" 
        return -1

    linearityparms = list()
    linearityerrs = list()
    parms = ['p0', 'p1', 'p2']
    for p in parms:
        _data_parm = find_segment_average(data, p, tube,wave)
        _data_err  = _data_cut[p+'Err'] 
        linearityparms.append( float(_data_parm) )
        linearityerrs.append( float(_data_err) )
    print linearityparms
 
    gainfactor = findGainFactor(run, tube)
#    new_value = value - correctionFunctionQuadratic(value, linearityparms)
    new_value = value - correctionFunctionQuadratic2(value, linearityparms, gainfactor)
 
    return new_value

 ### Averaged value washes out data 
def makeLEDLinCorr(value, run, tube, wave, parm_filename="../ELOGPics/AverageLEDLinearityParms_pol2.txt"):
    data = readFitParms(parm_filename)

    if run < min(data['Run_start']) or run > max(data['Run_end']):
        print "Run out of range"
        return -1

    # make cut conditions
    cutTube = data['tube'] == tube
    cutBegin = data['Run_start'] < run
    cutEnd = data['Run_end'] >= run # I think it's correct to lump with the following set but needs to be checked
    cutWave = data['wave'] == wave
    cutCondition = cutTube & cutBegin & cutEnd & cutWave
    
    _data_cut = data[cutCondition] 

    if len(_data_cut) == 0:
        print "Exiting" 
        return -1

    parms = ['p0', 'p1', 'p2']

    linearityparms = list()
    linearityerrs = list()
    for p in parms:
        cutParm = _data_cut['parm'] == p
        _data_parm = _data_cut[cutParm]
#        print "\nData_parm: " 
#        print _data_parm
#        print "\n"
        linearityparms.append( float(_data_parm['val']) )
        linearityerrs.append( float(_data_parm['meanerr']) )

    #    convert_factor = getConversionFactor(factor_filename, run, tube, wave)
    # this conversion factor is not relevant for PMT only correction -- 11/21/2014 SS

#    new_value = value - correctionFunctionQuadratic(value, linearityparms)
    new_value = value - correctionFunctionQuadratic2(value, linearityparms)

    print value
    print new_value
    
    return new_value

def correctionFunctionQuadratic(PMTval, parms):
    # 'parms' takes array of 3 parameters, [p0, p1, p2], where PMT = p0 + p1*PD * p2*PD*PD
    p0 = parms[0]
    p1 = parms[1]
    p2 = parms[2]
    arg = p1*p1 + 4*(PMTval - p0)*p2
    if arg < 0.0:
        print "negative sign under radical in correction function."
        print "returning 0"
        return 0

    correction = (1 / (2 * p2)) * ( p1*p1 + 2*(PMTval - p0)*p2 - p1*sqrt(p1*p1 + 4*(PMTval - p0)*p2) )
#    print str(p0) + " + " + str(p1) + "*PMT + " + str(p2) + "*PMT^2"

    return correction

def correctionFunctionQuadratic2(PMTval, parms, gainfactor =1, fitorigin = 782):  # 782 keV = beta endpt
    # 'parms' takes array of 3 parameters, [p0, p1, p2], where PMT = p0 + p1*(PD-fitorigin) + p2*(PD-fitorigin)^2
    # will take PMT' = PMT -  p2*(PMT-fitorigin)^2, ideally eliminating the quadratic behavior, assuming PD ~ PMT
    p2 = parms[2]
    correction = p2*(PMTval/gainfactor - fitorigin)**2

    return correction

if __name__ == "__main__":
#### testing
#    value = 782
#    print makeCorrections(value, 22000, 0, 405, "../ELOGPics/AverageLEDLinearityParms_pol2.txt")
#    for i in range(50):
#        print str(i*20) + "  " + str(i*20 - correctionFunctionQuadratic(i*20, [19.42, 11.42, -0.0115]))
#        if i > 0:
#            print correctionFunctionQuadratic(i*20, [19.42, 11.42, -0.0115])/float(i*i)
    value = 200
#    print makeLEDLinCorr(value, 22767, 0, 405, "../ELOGPics/AverageLEDLinearityParms_pol2.txt")


#    print makeLEDLinCorrByRun(value, 22767, 0, 405, "/data4/saslutsky/PulserComp/images_11_24_2014_10rangemin/FitResults.txt")


#    print makeLEDLinCorrByRun(value, 21274, 2, 405, "/data4/saslutsky/PulserComp/images_01_14_2015_scaleE/FitResults.txt")

    print makeLEDLinCorrByAvg(value, 21274, 2, 405, "/data4/saslutsky/PulserComp/images_01_14_2015_scaleE/FitResults.txt")
