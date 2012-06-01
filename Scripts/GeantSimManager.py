#!/usr/bin/python
import os
import time
from math import *

# killall -9 GeantSimManager.py; killall -9 parallel; killall -9 ucnG4_prod

class GeantSimManager:
	
	def __init__(self, simName, vacuum="1.e-5 torr",
				 fmap=None, geometry="C", sourceHolderPos = None):
		
		self.podsh = False				
		self.settings = {}
		
		self.settings["simName"] = simName
		self.settings["geometry"] = geometry
		self.settings["vacuum"] = vacuum
		self.settings["sourceholderpos"] = sourceHolderPos
		self.settings["sourceRadius"] = "1.5 mm"
		self.settings["sourceScan"] = 0.
		self.settings["particle"] = "e-"
		self.settings["fieldmapcmd"] = "#/detector/fieldmapfile UNUSED"
		if fmap:
			self.settings["fieldmapcmd"] = "/detector/fieldmapfile "+fmap
		self.settings["physlist"] = "livermore"
		self.settings["scintstep"] = "1.0 mm"
		self.settings["extra_cmds"] = ""
		
		self.settings["vis_cmd"] = ""
				
	def enable_vis(self):
		
		self.settings["vis_cmd"] = "/vis/open HepRepFile\n"
		#self.settings["vis_cmd"] = "/vis/open OGLSX\n"
		self.settings["vis_cmd"] += "/vis/viewer/set/viewpointThetaPhi 90 0\n"
		self.settings["vis_cmd"] += "/vis/viewer/panTo 2.2 0\n"
		self.settings["vis_cmd"] += "/vis/viewer/zoom 10\n"
		self.settings["vis_cmd"] += "/vis/viewer/set/viewpointThetaPhi 10 20\n"
		self.settings["vis_cmd"] += "/vis/drawVolume\n"
		self.settings["vis_cmd"] += "/vis/viewer/flush\n"
		self.settings["vis_cmd"] += "/vis/modeling/trajectories/create/drawByCharge myTrackVis\n"
		self.settings["vis_cmd"] += "/vis/modeling/trajectories/myTrackVis/default/setDrawStepPts true\n"
		self.settings["vis_cmd"] += "/vis/modeling/trajectories/myTrackVis/default/setDrawAuxPts true\n"
		self.settings["vis_cmd"] += "/vis/modeling/trajectories/select myTrackVis\n"
		self.settings["vis_cmd"] += "/vis/scene/add/trajectories\n"
		#self.settings["vis_cmd"] += "/vis/scene/add/trajectories rich\n"
		self.settings["vis_cmd"] += "/vis/scene/add/hits\n"
		
		#self.settings["vis_cmd"] += "/tracking/verbose 2\n"


	def set_generator(self,generator,forcePositioner=None):
		
		self.settings["generator"] = generator
		self.settings["gunenergy"] = 0
		
		# whether to construct the source holder and use source positioning
		self.settings["positioner"] = "DecayTrapFiducial" # DecayTrapUniform
		self.settings["gunpos_mm"] = [0.,0.,0.]
		self.settings["magf"] = "on"
		if self.settings["generator"] in ["In114E","In114W"]:
			self.settings["makeinfoil"] = "y"
		else:
			self.settings["makeinfoil"] = "n"
		needsHolder = False
		
		if self.settings["generator"][:2] in ["Bi","Ce","Sn","Cd","In","Sr","Cs"]:
			needsHolder = True
			self.settings["positioner"] = "SourceDrop"
		if self.settings["generator"][:2] in ["Ce","Cd"]:
			self.settings["sourceRadius"] = "1.25 mm"
		if self.settings["generator"] == "In114E":
			self.settings["gunpos_mm"] = [0.,0.,-.004975]
		if self.settings["generator"] == "In114W":
			self.settings["gunpos_mm"] = [0.,0.,.004975]
		if self.settings["generator"] == "eGun":
			needsHolder = True
			self.settings["positioner"] = "Fixed"
			self.settings["gunpos_mm"] = [0.,0.,-1000]
			self.settings["magf"] = "off"
		if self.settings["generator"][:2] == "Xe":
			self.settings["positioner"] = "UniformRadialGasFill"
		if forcePositioner:
			self.settings["positioner"] = forcePositioner
		
		if not self.settings["sourceholderpos"]:	
			if needsHolder:
				self.settings["sourceholderpos"] = "0 0	0 m"
			else:
				self.settings["sourceholderpos"] = "0 0.5 0 m"
		
		self.settings["run_num"] = 0
	
	
	#def throw_lines(self,nLines,eStart,eStop,logSpace):
	#
	#	energyLines = [0.5*(eStart+eStop),]
	#	if nLines > 1:
	#		if logSpace:
	#			energyLines = [ exp(log(eStart)+i*log(float(eStop)/eStart)/(nLines-1)) for i in range(nLines) ]
	#		else:
	#			energyLines = [eStart+i*(eStop-eStart)/float(nLines-1) for i in range(nLines) ]
	#	
	#	for self.settings["gunenergy"] in energyLines:
	
	def set_dirs(self):
		self.g4_workdir = os.environ["G4WORKDIR"]
		self.g4_bindir = os.environ["G4BINDIR"]
		if self.podsh:
			g4_workdir = "~/geant4/"
		self.type_dir = self.settings["simName"]+"_%s"%(self.settings["generator"].replace("/","_"))
		if self.settings["gunenergy"]:
			self.type_dir += "_%.1fkeV"%self.settings["gunenergy"]
		self.g4_out_dir = self.g4_workdir+"/output/%s/"%self.type_dir
		self.g4_log_dir = self.g4_workdir+"/logs/%s/"%self.type_dir
		self.g4_macro_dir = self.g4_workdir+"/macros/%s/"%self.type_dir
		self.g4_out_name = "%s/g4_run_%%s.root"%self.g4_out_dir
	
	def launch_sims(self,nEvents,nClusters=6,hours_old=0):
		
		nruns = 0
		nperclust = 0
		if self.podsh:
			nperclust = 8
		else:
			import multiprocessing
			nperclust = multiprocessing.cpu_count()	# number of simulation files to produce (generated in parallel)
		nruns = nperclust*nClusters
		if not nruns:
			nruns = 1
		oldtime = time.time() - hours_old*3600
		
		self.set_dirs()
		parallel_jobfile = "%s/jobs.txt"%self.g4_macro_dir
		
		os.system("mkdir -p %s"%self.g4_macro_dir)
		os.system("mkdir -p %s"%self.g4_out_dir)
		os.system("mkdir -p %s"%self.g4_log_dir)
		
		# account for low Bi conversion efficiency, ~14.3% as many electrons thrown as events run, mostly Auger
		ineffic_mul = 1
		if self.settings["generator"] == "Bi207":
			ineffic_mul = 4
		if self.settings["generator"] == "Ce139":
			ineffic_mul = 2
		
		# main simulations
		os.system("rm -r %s/*"%self.g4_macro_dir)
		jobsout = None
		if not self.podsh:
			jobsout = open(parallel_jobfile,"w")
		ucnG4_prod = self.g4_bindir+"/ucnG4_prod"
		onejob = ""
		subcmds = []
		
		# set up macros for each job
		for rn in range(nruns):
			self.settings["run_num"] += 1
			self.settings["jobname"] = self.settings["simName"]+"_%i"%self.settings["run_num"]
			self.settings["outfile"]=self.g4_out_name%str(self.settings["run_num"])
			self.settings["nevt"]=(ineffic_mul*nEvents)/nruns
			self.settings["joblog"] = "%s/gen_macro_%i.txt"%(self.g4_log_dir,self.settings["run_num"])
			g4_sub_file = "%s/geantjob_%i.sub"%(self.g4_macro_dir,self.settings["run_num"])
			
			# source position scan
			if self.settings["sourceScan"]:
				xpos = (((rn%nperclust)*nClusters+(rn/nperclust))/float(nruns-1)-0.5)*self.settings["sourceScan"]
				self.settings["gunpos_mm"][0] = xpos;
				if self.settings["sourceholderpos"] != "0 0.5 0 m":
					self.settings["sourceholderpos"] = "%g 0 0 mm"%xpos
			self.settings["gunpos"] = "%g %g %g mm"%tuple(self.settings["gunpos_mm"])
			
			# estimate wall time
			twall = int((60.0 + 0.2*self.settings["nevt"])/60)
			self.settings["walltime"]="%i:%02i:00"%(twall/60,twall%60)
			
			if os.path.exists(self.g4_out_name%str(self.settings["run_num"])) and os.stat(self.g4_out_name%str(self.settings["run_num"])).st_mtime > oldtime:
				continue;
			onejob = ucnG4_prod + " %s/geantgen_%i.mac %s"%(self.g4_macro_dir,self.settings["run_num"],self.settings["physlist"])
			# geant macro
			open(os.path.expanduser("%s/geantgen_%i.mac"%(self.g4_macro_dir,self.settings["run_num"])),"w").write(open("GeantGenMacroTemplate.mac","r").read()%self.settings)
			# submission script
			self.settings["jobcmd"] =  "mkdir -p %s\n"%self.g4_out_dir
			self.settings["jobcmd"] += "mkdir -p %s\n"%self.g4_log_dir
			self.settings["jobcmd"] += onejob+"\n"
			#self.settings["jobcmd"] += "rsync -av -e ssh %s 'mmendenhall@ucn.krl.caltech.edu:%s/'\n"%(self.settings["outfile"],g4_out_dir)
			#self.settings["jobcmd"] += "rm %s\n"%self.settings["outfile"]
			#self.settings["jobcmd"] += "rsync -av -e ssh %s 'mmendenhall@ucn.krl.caltech.edu:%s/'\n"%(self.settings["joblog"],g4_log_dir)
			#self.settings["jobcmd"] += "rm %s\n"%self.settings["joblog"]
			open(os.path.expanduser(g4_sub_file),"w").write(open("GeantJobTemplate.sub","r").read()%self.settings)
			subcmds.append("podsh submit --stageout='%s':'%s' %s"%(self.settings["joblog"],self.settings["joblog"],g4_sub_file))
			
			if jobsout:
				jobsout.write(onejob+" > %s 2>&1\n"%self.settings["joblog"])
		
		if jobsout:
			jobsout.close()
		
		if self.podsh:
			# send executable, macro to remote system
			print "Staging in macro files..."
			os.system("podsh stagein --files=~/geant4/bin:'~/geant4',~/geant4/macros:'~/geant4',~/geant4/logs:'~/geant4'")
			# submit each job
			for cmd in subcmds:
				print cmd
				os.system(cmd)
		else:	
			print "Running simulation jobs..."
			os.system("cat "+parallel_jobfile)
			if nruns > 1:
				os.system("nice -n 20 parallel < %s"%parallel_jobfile)
			else:
				os.system(onejob)
			os.system("rm "+parallel_jobfile)
	
	
	def launch_postanalyzer(self,nMin=0,nMax=100000):
		print "Running post analyzer..."
		self.set_dirs()
		resim_jobfile = "%s/resim_jobs.txt"%self.g4_macro_dir
		jobsout = open(resim_jobfile,"w")
		anafiles = [ (int(f[:-5].split("_")[-1]),self.g4_out_dir+"/"+f) for f in os.listdir(self.g4_out_dir) if f[:7]=="g4_run_"]
		anafiles.sort()
		nanalyzed = 0
		self.settings["analyzer"]="UCNA_MC_Analyzer"
		if self.settings["geometry"]=="siDet":
			self.settings["analyzer"]="SiDet_Analyzer"
		while anafiles:
			outlist_name = self.g4_out_dir+"outlist_%i.txt"%nanalyzed
			fout = open(outlist_name,"w")
			for f in anafiles[:6]:
				fout.write(f[1]+"\n")
			fout.close()
			anafiles = anafiles[6:]
			print "\n----- %s ------"%outlist_name
			os.system("cat "+outlist_name)
			allpts = ""
			if self.settings["generator"] in ["neutronBetaUnpol","eGunRandMomentum","eGun"]:
				allpts = " y"
			analyzer_bin = self.g4_bindir+"/"+self.settings["analyzer"]
			if nMin <= nanalyzed <= nMax:
				jobsout.write("%s %s %s/analyzed_%i.root%s\n"%(analyzer_bin,outlist_name,self.g4_out_dir,nanalyzed,allpts))
			nanalyzed += 1
		jobsout.close()
		print "\n----- %s ------"%resim_jobfile
		os.system("cat "+resim_jobfile)
		print
		os.system("nice -n 10 parallel < %s"%resim_jobfile)
		os.system("rm %s/outlist_*.txt"%self.g4_out_dir)
		os.system("rm "+resim_jobfile)


if __name__ == "__main__":
	
	####################				
	# neutrons
	####################
	
	# unpolarized beta baseline: 5e7 in 52 clusters
	if 0:
		betaSim = GeantSimManager("Test")
		betaSim.settings["physlist"]="livermore"
		betaSim.set_generator("neutronBetaUnpol")
		betaSim.launch_sims(nEvents=2e4,nClusters=1,hours_old=0)
		betaSim.launch_postanalyzer()
	
	# beta decay in magnetic field wiggles, 1e-3 vacuum: 1e7 in 52 clusters
	if 0:
		betaSim = GeantSimManager("WideKev",vacuum="1.e-3 torr",fmap="/home/mmendenhall/UCNA/Aux/Fieldmap_20101028_b.txt")
		betaSim.settings["physlist"]="livermore"
		betaSim.set_generator("neutronBetaUnpol")
		betaSim.launch_sims(nEvents=1e7,nClusters=52,hours_old=1000)
		betaSim.launch_postanalyzer(16)
	
	####################				
	# calibration sources
	####################
	
	# sources ["Bi207","Sn113","Ce139","Cd109","In114E","In114W","Cd113m"] 1e6 each
	if 0:
		for g in ["Bi207","Sn113","Ce139","Cd109"]:
			sourceSim = GeantSimManager("PosLine",fmap="/home/mmendenhall/UCNA/Aux/Fieldmap_20101028_b.txt")
			sourceSim.settings["physlist"]="livermore"
			sourceSim.settings["sourceScan"]=80.
			sourceSim.set_generator(g)
			sourceSim.launch_sims(nEvents=1e6,nClusters=12,hours_old=0)
			sourceSim.launch_postanalyzer()
	
	if 0:
		for g in ["Bi207","Sn113","Ce139"]:
			sourceSim = GeantSimManager("DeepCrinkle",fmap="/home/mmendenhall/UCNA/Aux/Fieldmap_20101028_b.txt")
			sourceSim.settings["physlist"]="livermore"
			sourceSim.settings["sourceScan"]=80.
			sourceSim.settings["extra_cmds"] += "/detector/foilcrinkle 1.5708\n"
			sourceSim.set_generator(g)
			#sourceSim.enable_vis()
			sourceSim.launch_sims(nEvents=1e6,nClusters=12,hours_old=0)
			sourceSim.launch_postanalyzer()


	####################				
	# Xenon ["Xe135_3-2+","Xe125_1-2+","Xe133_3-2+"], 3e7 for main and 3e6 for secondary
	####################
	if 1:
		for g in ["Xe133_3-2+"]:
			sourceSim = GeantSimManager("WideKev",vacuum="1.e-3 torr")
			sourceSim.settings["physlist"]="livermore"
			sourceSim.set_generator(g)
			sourceSim.launch_sims(nEvents=3e6,nClusters=18,hours_old=0)
			sourceSim.launch_postanalyzer()

	
	####################				
	# silicon detector
	####################
	
	# Silicon detector test
	if 0:
		siDet = GeantSimManager("SiDet",geometry="siDet")
		siDet.set_generator("Cs137")
		siDet.settings["extra_cmds"] += "/sourceholder/windowthickness 1.5 mm\n"
		siDet.launch_sims(nEvents=1e6,nClusters=6,hours_old=0)
		siDet.launch_postanalyzer()

	# Alpha particles through foil
	if 0:
		for th in range(11):
			siDet = GeantSimManager("AlphaFoil_%i"%th,geometry="siDet")
			siDet.set_generator("eGun")
			siDet.settings["extra_cmds"] += "/sourceholder/windowthickness %i um\n"%th
			siDet.settings["particle"] = "alpha"
			siDet.settings["gunenergy"] = 5485.56
			siDet.launch_sims(nEvents=1e4,nClusters=1,hours_old=0)
			siDet.launch_postanalyzer()

# visualization test
#launch_simulations(generators = ["eGunRandMomentum"], forcePositioner="SourceDrop",
#					nEvents = 99, nClusters=0, folderPrefix = "VisTest",
#					nMonoLines=1, eStart=60.0, eStop=60.0, hours_old = 0, resimOnly=False)
