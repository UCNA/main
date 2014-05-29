#!/sw/bin/python2.7

import sys
sys.path.append("..")
import os
from ucnacore.EncalDB import *
from ucnacore.QFile import *

import os
import string

ppSegments = [["A1","A2","A4","A5"],["A7","A9","A10","A12"],["B1","B2","B4","B5"],["B7","B9","B10","B12"]]
hoSegments = [ppSegments[0]+ppSegments[1],ppSegments[2]+ppSegments[3]]
octSegments = [hoSegments[0]+hoSegments[1],]

##########################
#                        #
#   A1   Background Off  #
#   A2   Beta Off        #
#   A3   Depol Off->On   #
#   A4   Background On   #
#   A5   Beta On         #
#   A6   Depol On->Off   #
#   A7   Beta On         #
#   A8   Depol On->Off   #
#   A9   Background On   #
#   A10  Beta Off        #
#   A11  Depol Off->On   #
#   A12  Background Off  #
#                        #
#   B1   Background On   #
#   B2   Beta On         #
#   B3   Depol On->Off   #
#   B4   Background Off  #
#   B5   Beta Off        #
#   B6   Depol Off->On   #
#   B7   Beta Off        #
#   B8   Depol Off->On   #
#   B9   Background Off  #
#   B10  Beta On         #
#   B11  Depol On->Off   #
#   B12  Background On   #
#                        #
##########################

fgRuns = ["A2","A5","A7","A10","B2","B5","B7","B10"]
bgRuns = ["A1","A4","A9","A12","B1","B4","B9","B12"]

def load_octetruns(fname, tplist = octSegments[0]):
	rnlist = []
	for m in QFile(fname).dat.get("Octet",[]):
		for k in m.dat:
				if k in tplist:
					for rns in m.dat[k]:
						for rn in rns.split(","):
							rnlist.append(int(rn))
	return rnlist

def collectPlots(pattern, outbase, rlist, inDir=os.environ['UCNA_ANA_PLOTS']+"/figures/"):
		
	fname = pattern.split("/")[-1]	# copy target file name
	outdir = outbase+"/"+string.join(pattern.split("/")[:-1],'/')+'/'+fname+'/' # output directory
	os.system("mkdir -p "+outdir)

	for r in rlist:
		targfile = "%s/run_%i/%s"%(inDir,r,pattern)
		if os.path.exists(targfile):
			cmd = "cp %s %s/%i_%s"%(targfile,outdir,r,fname)
			print cmd
			os.system(cmd)
		else:
			print "Missing plot for run",r
		
if __name__=="__main__":
	conn = open_connection()
	#rlist = getRunType(conn,"Asymmetry",20542,23173)
	rlist = load_octetruns(os.environ["UCNA_AUX"]+"/OctetList_20122013.txt",fgRuns)
	#for t in range(4):
	#	collectPlots("PMTs/BiPulser_W%i.pdf"%t,"/Users/michael/Desktop",rlist)
	collectPlots("Rates/EventRates.pdf","/Users/michael/Desktop",rlist)

