import os
import sys
sys.path.append("~/MBpython")
import MButils

#taking care of bad runs when run through simulated Xe to fit for isotope composition
if 0:
    runs = [19890]
    pcmd = "./UCNAnalyzer pmap sim %i %i 12 x x"
    for run in runs:
        print "Running " + pcmd%(run, run)
        os.system(pcmd%(run, run))


################### Set the year #################
year = 2011
fit='quad'
runlow=[]
runhigh=[]
missingRuns=[] #list of runs which seem to be missing from data
rerunAllSimMatching=True #Whether or not to automatically rerun matching of 
                          #Simulation to data to determine isotopic comp.
                          #This should be true unless trying to track down
                          #a problem!

if 1:
#NOTE: REPLAYING XENON RUNS ISN'T NECESSARY EVERY ITERATION
#Replaying Xenon runs and extracting section-by-section spectra for individual runs
# runs 21596 and beyond are from 2012/2013 and should be replayed again once new sims have been run
# make sure you change the environment path of the srcsims and source the .bashrc file
    if 1:
        if year==2011:
            runlow = [17570, 18081, 18390, 18712, 19589, 19873]
            runhigh = [17600, 18090, 18413, 18744, 19606, 19898]
        else:
            runlow=[21596, 21966, 22961]           
            runhigh = [21605, 22003, 22979]
        #pcmd1 = "cd Scripts/launchers; ./OfficialReplayManager.py -x --rmin=%i --rmax=%i"
        pcmd2 = "cd Scripts/launchers; ./ReplayManager.py -x --rmin=%i --rmax=%i"
    
        for i in range(0, len(runlow), 1):
            #os.system(pcmd1%(runlow[i], runhigh[i]))
            os.system(pcmd2%(runlow[i], runhigh[i]))


# Merging the individual runs into one file. Range may not match that above
# based on looking at Erecon 
    if 1:
        if year==2011:
            runlow = [17570, 18082, 18391, 18712, 19873]
            runhigh = [17590, 18090, 18413, 18744, 19898]
        elif year==2012:
            runlow=[21596, 21966, 22962]           
            runhigh = [21605, 22003, 22979]
                                        
        pcmd1 = "./UCNAnalyzer pmap gen %i %i 12 x x"
    
        for i in range(0, len(runlow), 1):
                os.system(pcmd1%(runlow[i], runhigh[i]))

# determine isotopic composition via simulation
    if 1:
        redoRuns=[]
        failedRuns=[]
        path=' '
        if fit=='lin':
            path = '/extern/mabrow05/ucna/Analysis_Output/PositionMaps/SingleRunsSim/'
        elif fit=='quad':
            path = '/extern/mabrow05/ucna/Analysis_Output_quad/PositionMaps/SingleRunsSim/'
        if year==2011:
            runlow = [17570, 18081, 18390, 18712, 19873]
            runhigh = [17590, 18090, 18413, 18744, 19898]       
        elif year==2012:
            runlow=[21596, 21966, 22961]           
            runhigh = [21605, 22003, 22979]
        
        pcmd1 = "cd Scripts/launchers; ./ReplayManager.py -X --rmin=%i --rmax=%i"

        for i in range(0, len(runlow), 1):
            if rerunAllSimMatching:
                os.system(pcmd1%(runlow[i], runhigh[i]))
        
            for run in range(runlow[i], runhigh[i]+1, 1):
                if not os.path.exists(path+'/Xenon_%i_12_52'%run) or MButils.fileIsOlderThan(path+'/Xenon_%i_12_52'%run,2):
                    redoRuns.append(run)

        print redoRuns
        pcmd = "./UCNAnalyzer pmap sim %i %i 12 x x"
        for redo in redoRuns:
            print "Running " + pcmd%(redo, redo)
            os.system(pcmd%(redo, redo))
            if not os.path.exists(path+'/Xenon_%i_12_52'%redo):
                missingRuns.append(redo)
                continue
            if MButils.fileIsOlderThan(path+'/Xenon_%i_12_52'%redo,2):
                failedRuns.append(redo)

        if len(failedRuns)>0:
            print "These runs failed" 
            print failedRuns
            for r in failedRuns:
                tries = 0
                while MButils.fileIsOlderThan(path+'/Xenon_%i_12_52'%r,2):
                    os.system(pcmd%(r,r))
                    tries+=1
                    if tries==5:
                        exit(0) #quit after so many tries with same run...
            


#Merge simulated runs. These ranges should match that of the data which was merged
    if 1:
        if year==2011:
            runlow = [17570, 18082, 18391, 18712, 19873]
            runhigh = [17590, 18090, 18413, 18744, 19898]
        elif year==2012:
            runlow=[21596, 21966, 22962]           
            runhigh = [21605, 22003, 22979]
        
        pcmd1 = "./UCNAnalyzer pmap sim %i %i 12 x x"
    
        for i in range(0, len(runlow), 1):
            os.system(pcmd1%(runlow[i], runhigh[i]))

#make maps and plots
    if 1:
        if year==2011:
            runlow = [17570, 18082, 18391, 18712, 19873]
            runhigh = [17590, 18090, 18413, 18744, 19898]
        elif year==2012:
            runlow=[21596, 21966, 22962]           
            runhigh = [21605, 22003, 22979]
        
        #pmap = [275,277,279,281,283]
        pcmd1 = "./UCNAnalyzer pmap comp %i %i 12 x x"
        pcmd2 = "./UCNAnalyzer pmap plot %i x x"
    
        for i in range(0, len(runlow), 1):
            os.system(pcmd1%(runlow[i], runhigh[i]))
            #os.system(pcmd2%(pmap[i]))


 
    print 'These Runs seem to be missing from data...'
    print missingRuns

       
         
