#!/usr/bin/python

#'official analyzer' tube numbering	
jTubes = {'E':[3,4,1,2,5],'W':[1,2,3,4,5]}

# listing of all PMTs
allPMTs = [ [(s,t) for t in range(4)] for s in ['E','W'] ]
allPMTs = allPMTs[0]+allPMTs[1]

# various sets of runs
betagroups = [	(13905,13964), (14077,14100), (14127,14261), (14356,14380), (14397,14421), (14432,14467),
				(14535,14667), (14688,14782), (14888,14994), (15084,15150), (15172,15356), (15448,15639), (15667,15915),
				(15943,15966), (16097,16216),
				(18801,18994),
				(19023,19046),
				(19050,19201)]
sourcegroups = [(13883,13894), (14104,14118), (14383,14394), (14516,14530), (14736,14746), (15645,15662), (15916,15939),
				#(15357,15371),	# badly tilted sources
				(16240,16257),	# 2010 end of year sources
				(17233,17249),
				(17359,17387),
				(17517,17527),
				(17871,17922),
				(17925,17956),	# long centered runs
				(18020,18055),	# old and new Cd sources
				(18357,18386),	# new In source
				(18617,18640),	#
				(18745,18768),	#
				(19203,19239),	# Jan 23
				(19347,19377),	# Feb 10
				(19505,19544)	# Feb 14, Cd/In only
				]
thecalruns = [13890,14111,14390,14524,14743,15653,15931]
xegroups = [	(14264,14347),	# 2010 #1
				(15991,16077),	# 2010 #2
				(16983,17078),	# 2011 start; PMT W2 dead
				(17224,17230),	# old Xenon
				(17561,17734),	# long, W0 pulser dead
				(18081,18090),	# everything working; not much Xe produced; short set
				(18390,18413),	# everything working; W0 pulser very low (& badly gain calibrated)
				(18712,18744),	# 2012 start; W0 pulser still low
				(19589,19606)	# Feb. 16, short-ish high rates run
				]
				
def expandRanges(rgs,rmin=0,rmax=100000):
	rns = []
	for rg in rgs:
		rns += [r for r in range(rg[0],rg[1]+1) if rmin <= r <= rmax]
	return rns
	