#!/usr/bin/python

import os
import sys

runs2011 = [(17359,17387),(17517,17527),(17871,17941),(18020,18055),(18357,18386),
        (18617,18640),(18745,18768),(19203,19239),(19347,19377),(19823,19863)]
runs2012 = [(20515,20527),(20818,20830),(21087,21099),(21299,21328),(21679,21718),
        (21914,21939),(22215,22238),(22437,22462),(22767,22793)]

#Runs replay of source runs
#for i in range(0,len(runs2011),1):
    #os.system("./OfficialReplayManager.py -s --rmin=%i --rmax=%i"%(runs2011[i][0],runs2011[i][1]))

#Creates Data to monte carlo comparisons and fills DB with peak information for making error envelope
for i in range(0,len(runs2011),1):
    os.system("./ReplayManager.py -s --rmin=%i --rmax=%i"%(runs2011[i][0],runs2011[i][1]))

#for i in range(0,len(runs2012),1):
 #   os.system("./ReplayManager.py -s --rmin=%i --rmax=%i"%(runs2012[i][0],runs2012[i][1]))

