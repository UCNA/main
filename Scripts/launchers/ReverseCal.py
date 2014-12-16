import os
from math import *

runs = [18433, 18436, 18438, 18441, 18445, 18448, 18450, 18453]

octetNum = 25
fiducialRad = 50.

outputDir = "/extern/mabrow05/ucna/SimAnalysisOutput/"

octetFile = open(outputDir+"BetaRunsInOctet_2011-2012.txt","r")
octets = octetFile.readlines()
#print octets[octetNum]
octetRuns = [rn for rn in octets[octetNum].split(" ") if rn!='\n']
#print octetRuns

os.system("mkdir %s/hists/Octet_%i_%s-%s"%(outputDir,octetNum,octetRuns[0],octetRuns[len(octetRuns)-1]))

for run in octetRuns:
    cmd = "./MC_Plugin_Analyzer beta %s %s %f" %(outputDir,run,fiducialRad)
    os.system("cd ../../; %s"%(cmd))
    #print cmd
    os.system("cd %s/hists/; mv %s_sim.root Octet_%i_%s-%s/%s_sim_fid%.0f.root"%(outputDir,run,octetNum,octetRuns[0],octetRuns[len(octetRuns)-1],run,fiducialRad))





