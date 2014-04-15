#!/usr/bin/python

import os

if __name__ == "__main__":
	
	simfldir = os.environ["G4OUTDIR"]+"/AFP_Fringe_n1_f_n"
	baseout = os.environ["UCNA_ANA_PLOTS"]+"/test/SimPosOff"
	cmds = open("cmdlist.txt","w")
	
	pcmd = "../../MC_Plugin_Analyzer posoff "
	for f in os.listdir(simfldir):
		if f.split("_")[0]=="analyzed":
			outdir = baseout+"/"+f[:-5]
			os.system("mkdir -p "+outdir)
			flist = open(outdir+"/simlist.txt","w")
			flist.write("%s/%s\n"%(simfldir,f))
			flist.close()
			cmds.write(pcmd+outdir+"\n")
			
	cmds.close();
	os.system("cat cmdlist.txt")
	os.system("nice -n 20 parallel < cmdlist.txt")
