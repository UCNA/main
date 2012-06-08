#!/usr/bin/python
import os
import time
from EncalDB import *
from optparse import OptionParser

def replay_betas(rmin,rmax):	
	
	pcmd = "cd ../OfficialReplay; ./ucnaDataAnalyzer11b %i cutbeam ledtree\n"
	rlist = getRunType(open_connection(),"Asymmetry",rmin,rmax)	
	freplaylist = open("officialreplaylist.txt","w")
	rlist.sort()
	for r in rlist:
		freplaylist.write(pcmd%r)
	freplaylist.close()
	print "------- Replay List -------"
	os.system("cat officialreplaylist.txt")
	os.system("nice -n 10 parallel < officialreplaylist.txt")
	os.system("rm officialreplaylist.txt")
	
def replay_new(rmin,rmax):
	
	rlist = range(rmin,rmax+1)
	pcmd = "cd ../OfficialReplay; ./ucnaDataAnalyzer11b %i\n"	
		
	freplaylist = open("officialreplaylist.txt","w")
	
	rlist.sort()
	for r in rlist:
		freplaylist.write(pcmd%r)
	freplaylist.close()
	print "------- Replay List -------"
	os.system("cat officialreplaylist.txt")
	os.system("nice -n 10 parallel < officialreplaylist.txt")
	os.system("rm officialreplaylist.txt")

	
def replay_sources(rmin,rmax,doXe=False):
	
	rlist = getRunType(open_connection(),"SourceCalib",rmin,rmax)
	if doXe:
		rlist = getRunType(open_connection(),"Xenon",rmin,rmax)
	pcmd = "cd ../OfficialReplay; ./ucnaDataAnalyzer11b %i\n"
			
	freplaylist = open("officialreplaylist_sources.txt","w")
	
	rlist.sort()
	for r in rlist:
		freplaylist.write(pcmd%r)
	freplaylist.close()
	print "------- Replay List -------"
	os.system("cat officialreplaylist_sources.txt")
	os.system("nice -n 10 parallel < officialreplaylist_sources.txt")
	os.system("rm officialreplaylist_sources.txt")
	
	
if __name__ == "__main__":
		
	parser = OptionParser()
	parser.add_option("-k", "--kill", dest="kill", action="store_true", default=False, help="kill running replays")
	parser.add_option("-r", "--replay", dest="replay", action="store_true", default=False, help="run official replay for beta decay data")
	parser.add_option("-s", "--sources", dest="sources", action="store_true", default=False, help="run official replay for source data")
	parser.add_option("-x", "--xenon", dest="xenon", action="store_true", default=False, help="run official replay for xenon data")
	parser.add_option("-n", "--new", dest="newruns", action="store_true", default=False, help="run official replay for new data")
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
		