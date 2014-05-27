#!/usr/bin/python
# nohup ./ReplayManager.py --mwpccal_sim < /dev/null > scriptlog.txt 2>&1 &
# nohup ./ReplayManager.py -s --rmin=20818 --rmax=20828 < /dev/null > scriptlog.txt 2>&1 &
# nohup ./ReplayManager.py -o < /dev/null > scriptlog.txt 2>&1 &

# nohup ./ReplayManager.py -x --rmin=19873 --rmax=19898 --nrings=12 < /dev/null > scriptlog.txt 2>&1 &
# nohup ./ReplayManager.py -X --rmin=17561 --rmax=17650 --nrings=12 < /dev/null > scriptlog.txt 2>&1 &

import os
import time
from optparse import OptionParser
import sys
sys.path.append("..")
from ucnacore.EncalDB import *


anaBinDir = "../../"

def processOctets(sim,omin,omax):
	pcmd = "cd "+anaBinDir+"; ./UCNAnalyzer pr oct %i x x\n"
	freplaylist = open("oct_replaylist.txt","w")
	for r in range(100)[omin:omax+1]:
		if sim:
			freplaylist.write(pcmd%(-r-1));
		else:
			freplaylist.write(pcmd%r);
	freplaylist.close()
	os.system("cat oct_replaylist.txt")
	if sim:
		os.system("nice -n 5 parallel -P 3 < oct_replaylist.txt")
	else:
		os.system("nice -n 5 parallel -P 4 < oct_replaylist.txt")
	os.system("rm oct_replaylist.txt")
	if sim:
		os.system("cd "+anaBinDir+"; ./UCNAnalyzer pr oct -1000 x x\n");
	else:
		os.system("cd "+anaBinDir+"; ./UCNAnalyzer pr oct 1000 x x\n");
	

def processSources(rmin,rmax):
		pcmd = "cd "+anaBinDir+"; ./UCNAnalyzer sr %i %i x\n"
		freplaylist = open("source_replaylist.txt","w")
		for r in getRunType(open_connection(),"SourceCalib",rmin,rmax):
			freplaylist.write(pcmd%(r,r))
		freplaylist.close()
		os.system("cat source_replaylist.txt")
		os.system("nice -n 10 parallel -P 5 < source_replaylist.txt")
		os.system("rm source_replaylist.txt")

def processXeMap(rmin,rmax,nr):
		pcmd = "cd "+anaBinDir+"; ./UCNAnalyzer pmap gen %i %i %i x x\n"
		freplaylist = open("xenon_replaylist.txt","w")
		for r in range(rmin,rmax+1):
			freplaylist.write(pcmd%(r,r,nr))
		freplaylist.close()
		os.system("cat xenon_replaylist.txt")
		os.system("nice -n 15 parallel -P 4 < xenon_replaylist.txt")
		os.system("rm xenon_replaylist.txt")
		os.system(pcmd%(rmin,rmax,nr))

def processXeSim(rmin,rmax,nr):
		pcmd = "cd "+anaBinDir+"; ./UCNAnalyzer pmap sim %i %i %i x x\n"
		freplaylist = open("xenon_simlist.txt","w")
		for r in range(rmin,rmax+1):
			freplaylist.write(pcmd%(r,r,nr))
		freplaylist.close()
		os.system("cat xenon_simlist.txt")
		nproc = 5
		if nr > 15:
			nproc = 3
		os.system("nice -n 15 parallel -P %i < xenon_simlist.txt"%nproc)
		os.system("rm xenon_simlist.txt")
		os.system(pcmd%(rmin,rmax,0,nr))

def betaOctetPositions(dosim):
	pcmd = "cd "+anaBinDir+"; ./BetaOctetPositions %i\n"
	freplaylist = open("parallel_cmd.txt","w")
	for r in range(60):
		if dosim:
			freplaylist.write(pcmd%(-r-1))
		else:
			freplaylist.write(pcmd%r)
	freplaylist.close()
	os.system("cat parallel_cmd.txt")
	os.system("nice -n 15 parallel -P 3 < parallel_cmd.txt")
	os.system("rm parallel_cmd.txt")
	if dosim:
		os.system(pcmd%(-1000))
	else:
		os.system(pcmd%1000)

def process_MWPCCal_octets(sim,omin,omax):
	pcmd = "cd "+anaBinDir+"; ./MWPC_Energy_Cal oct %%i %s 8 x\n"%{True:"sim", False:"dat"}[sim]
	freplaylist = open("mwpccal_replaylist.txt","w")
	for r in range(60)[omin:omax+1]:
		freplaylist.write(pcmd%r);
	freplaylist.close()
	os.system("cat mwpccal_replaylist.txt")
	os.system("nice -n 15 parallel -P 4 < mwpccal_replaylist.txt")
	os.system("rm mwpccal_replaylist.txt")


if __name__ == "__main__":
	
	parser = OptionParser()
	parser.add_option("-k", "--kill", dest="kill", action="store_true", default=False, help="kill running replays")
	parser.add_option("-o", "--octs", dest="octs", action="store_true", default=False, help="process octets")
	parser.add_option("-z", "--simocts", dest="simocts", action="store_true", default=False, help="re-simulate octets")
	parser.add_option("-x", "--xenon", dest="xenon", action="store_true", default=False, help="process Xenon runs for position map")
	parser.add_option("-X", "--xesim", dest="xesim", action="store_true", default=False, help="simulate Xenon runs for position map")
	parser.add_option("-s", "--sources", dest="sources", action="store_true", default=False, help="process calibration source data")
	parser.add_option("--mwpccal", dest="mwpccal", action="store_true", default=False, help="MWPC calibration octet analysis")
	parser.add_option("--mwpccal_sim", dest="mwpccal_sim", action="store_true", default=False, help="MWPC calibration octet analysis, simulation")
	parser.add_option("--rmin", type="int", dest="rmin", default=0)
	parser.add_option("--rmax", type="int", dest="rmax", default=100000)
	parser.add_option("--nrings", type="int", dest="nrings", default=12, help="number of rings for position map")
	parser.add_option("--bops", dest="bops", action="store_true", default=False, help="beta octet positions")
	
	options, args = parser.parse_args()
	if options.kill:
		os.system("killall -9 parallel")
		os.system("killall -9 UCNAnalyzer")
		os.system("killall -9 BetaOctetPositions")
		os.system("killall -9 ReplayManager.py")
		exit(0)
	
	if len(os.popen("ps -a | grep ReplayManager").readlines()) > 1:
		print "Already running! I die!"
		exit(1)
	
	if options.xenon:
		processXeMap(options.rmin,options.rmax,options.nrings)
		exit(0)
	
	if options.xesim:
		processXeSim(options.rmin,options.rmax,options.nrings)
		exit(0)
		
	if options.octs:
		processOctets(False,options.rmin,options.rmax)
		
	if options.simocts:
		processOctets(True,options.rmin,options.rmax)
		
	if options.sources:
		processSources(options.rmin,options.rmax)
		exit(0)

	if options.bops:
		betaOctetPositions(False)
		exit(0)

	if options.mwpccal:
		process_MWPCCal_octets(False,0,65)
		exit(0)

	if options.mwpccal_sim:
		process_MWPCCal_octets(True,0,65)
		exit(0)
