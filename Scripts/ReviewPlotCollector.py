#!/sw/bin/python2.7
import os
from EncalDB import *

import os
import string
					
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
	rlist = getRunType(conn,"Asymmetry",16500,23200)
	collectPlots("PMTs/Erecon_Type_1.pdf","/Users/michael/Desktop",rlist)
