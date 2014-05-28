#!/usr/bin/python

# nohup ./OfficialReplayManager.py -s --rmin=22767 --rmax=22793 < /dev/null > officreplaylog.txt 2>&1 &
# nohup ./OfficialReplayManager.py -r --rmin=20000 --rmax=30000 < /dev/null > officreplaylog.txt 2>&1 &

# nohup ./OfficialReplayManager.py -x --rmin=19589 --rmax=19606 < /dev/null > officreplaylog.txt 2>&1 &

import os
import time
from optparse import OptionParser
import sys
sys.path.append("..")
from ucnacore.EncalDB import *

def replay_betas(rmin,rmax):	
	
	pcmd = "cd ../../OfficialReplay; ./ucnaDataAnalyzer11b %i cutbeam ledtree\n"
	rlist = getRunType(open_connection(),"Asymmetry",rmin,rmax)
	rlist.sort()
	freplaylist = open("officialreplaylist.txt","w")
	
	for r in rlist:
		freplaylist.write(pcmd%r)
	freplaylist.close()
	print "------- Replay List -------"
	os.system("cat officialreplaylist.txt")
	os.system("nice -n 10 parallel < officialreplaylist.txt")
	os.system("rm officialreplaylist.txt")


def replay_new(rmin,rmax):
	
	rlist = range(rmin,rmax+1)
	pcmd = "cd ../../OfficialReplay; ./ucnaDataAnalyzer11b %i\n"
		
	freplaylist = open("officialreplaylist.txt","w")
	for r in rlist:
		freplaylist.write(pcmd%r)
	freplaylist.close()
	print "------- Replay List -------"
	os.system("cat officialreplaylist.txt")
	os.system("nice -n 10 parallel < officialreplaylist.txt")
	os.system("rm officialreplaylist.txt")


def replay_backscatters(rmin,rmax):
	rlist = getRunType(open_connection(),"Asymmetry",rmin,rmax)
	rlist.sort()
	pcmd = "cd ../../OfficialReplay; ./TriggerTree %i\n"
	freplaylist = open("backscatterreplaylist.txt","w")
	for r in rlist:
		freplaylist.write(pcmd%r)
	freplaylist.close()
	print "------- Replay List -------"
	os.system("cat backscatterreplaylist.txt")
	os.system("nice -n 10 parallel < backscatterreplaylist.txt")
	os.system("rm backscatterreplaylist.txt")


def replay_sources(rmin,rmax,doXe=False):
	
	rlist = getRunType(open_connection(),"SourceCalib",rmin,rmax)
	if doXe:
		rlist = getRunType(open_connection(),"Xenon",rmin,rmax)
	pcmd = "cd ../../OfficialReplay; ./ucnaDataAnalyzer11b %i\n"
			
	freplaylist = open("officialreplaylist_sources.txt","w")
	
	rlist.sort()
	for r in rlist:
		freplaylist.write(pcmd%r)
	freplaylist.close()
	print "------- Replay List -------"
	os.system("cat officialreplaylist_sources.txt")
	os.system("nice -n 10 parallel -P 4 < officialreplaylist_sources.txt")
	os.system("rm officialreplaylist_sources.txt")
	
	
if __name__ == "__main__":
		
	parser = OptionParser()
	parser.add_option("-k", "--kill", dest="kill", action="store_true", default=False, help="kill running replays")
	parser.add_option("-r", "--replay", dest="replay", action="store_true", default=False, help="run official replay for beta decay data")
	parser.add_option("-s", "--sources", dest="sources", action="store_true", default=False, help="run official replay for source data")
	parser.add_option("-x", "--xenon", dest="xenon", action="store_true", default=False, help="run official replay for xenon data")
	parser.add_option("-n", "--new", dest="newruns", action="store_true", default=False, help="run official replay for new data")
	parser.add_option("-b", "--bksc", dest="backscatter", action="store_true", default=False, help="run backscatter event extractor")
	parser.add_option("--rmin", type="int", dest="rmin", default=0)
	parser.add_option("--rmax", type="int", dest="rmax", default=100000)
	
	options, args = parser.parse_args()
	
	if options.kill:
		os.system("killall -9 parallel")
		os.system("killall -9 ucnaDataAnalyzer11b")
		os.system("killall -9 OfficialReplayManager.py")
		exit(0)

	if len(os.popen("ps -a | grep OfficialReplay").readlines()) > 1:
		print "Already running! I die!"
		exit(1)
		
	if options.xenon:
		replay_sources(options.rmin,options.rmax,doXe=True)
	if options.sources:
		replay_sources(options.rmin,options.rmax)
	if options.replay:
		replay_betas(options.rmin,options.rmax)
	if options.newruns:
		replay_new(options.rmin,options.rmax)
	if options.backscatter:
		replay_backscatters(options.rmin,options.rmax)
		