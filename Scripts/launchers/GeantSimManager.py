#!/usr/bin/python
import os
import time
from math import *
from optparse import OptionParser

# nohup ./GeantSimManager.py --xesrcs --evtgen < /dev/null > g4log.txt 2>&1 &
# nohup ./GeantSimManager.py --xesrcs --sim --ana < /dev/null > g4log.txt 2>&1 &
# nohup ./GeantSimManager.py --calsrcs --sim --ana < /dev/null > g4log.txt 2>&1 &

class GeantSimManager:
	
	def __init__(self, simName, vacuum="1.e-5 torr",
				 fmap=None, geometry="C", sourceHolderPos = None):
		
		self.settings = {}
		
		self.settings["run_num"] = 0
		self.settings["simName"] = simName
		self.settings["geometry"] = geometry
		self.settings["vacuum"] = vacuum
		self.settings["sourceholderpos"] = sourceHolderPos
		self.settings["sourceRadius"] = "0 mm"
		self.settings["sourceScan"] = 0.
		self.settings["fieldmapcmd"] = "#/field/mapfile UNUSED"
		if fmap:
			self.settings["fieldmapcmd"] = "/field/mapfile "+fmap
		#self.settings["ana_args"] = "undead cathodes"
		self.settings["ana_args"] = "cathodes"
		
		self.settings["extra_cmds"] = ""
		self.settings["extra_cmds"] += "/detector/MWPCBowing 5 mm\n"
		self.settings["extra_cmds"] += "/detector/sourcefoilthick 9.5 um\n"
		
		self.settings["extra_post_cmds"] = ""
		
		self.settings["vis_cmd"] = ""
		
		self.g4_out_dir_base = "/extern/mabrow05/ucna/geant4work/SrcSims"
		
		self.anagroup = 10 # number of files to group together for final analyzer result
				
	def enable_vis(self):
		
		self.settings["vis_cmd"] = "/vis/open HepRepFile\n"
		#self.settings["vis_cmd"] = "/vis/open OGLIX\n"
		
		#self.settings["vis_cmd"] = "/vis/open OGLIXm\n"
		#self.settings["vis_cmd"] = "/vis/open OGLIQt\n"
		
		self.settings["vis_cmd"] += "/vis/drawVolume\n"
		
		self.settings["vis_cmd"] += "/vis/viewer/set/viewpointThetaPhi 90 0\n"
		self.settings["vis_cmd"] += "/vis/viewer/panTo 2.2 0\n"
		self.settings["vis_cmd"] += "/vis/viewer/zoom 10\n"
		#self.settings["vis_cmd"] += "/vis/viewer/set/viewpointThetaPhi 10 20\n"
		
		self.settings["vis_cmd"] += "/vis/viewer/set/auxiliaryEdge true\n"
		self.settings["vis_cmd"] += "/vis/viewer/set/style surface"
		
		#self.settings["vis_cmd"] += "/vis/viewer/set/hiddenEdge 1"
		
		if 1:
			self.settings["vis_cmd"] += "/vis/modeling/trajectories/create/drawByCharge myTrackVis\n"
			#self.settings["vis_cmd"] += "/vis/modeling/trajectories/myTrackVis/default/setDrawStepPts true\n"
			#self.settings["vis_cmd"] += "/vis/modeling/trajectories/myTrackVis/default/setDrawAuxPts true\n"
			self.settings["vis_cmd"] += "/vis/modeling/trajectories/select myTrackVis\n"
			self.settings["vis_cmd"] += "/vis/scene/add/trajectories\n"
			self.settings["vis_cmd"] += "/vis/scene/add/trajectories rich\n"
			self.settings["vis_cmd"] += "/vis/scene/add/hits\n"
		
		#self.settings["vis_cmd"] += "/tracking/verbose 2\n"
		
		self.settings["vis_cmd"] += "/vis/viewer/flush\n"
	
	def set_evtsrc(self,evtsrc):
	
		self.settings["evtsrc"] = evtsrc
		self.settings["evtsrc_name"] = evtsrc
		self.needsHolder = False
		
		# special settings for asymmetric In114 source
		if self.settings["evtsrc"] in ["In114E","In114W"]:
			self.settings["extra_cmds"] += "/detector/infoil true\n"
			if self.settings["evtsrc"][-1] == "E":
				self.settings["extra_post_cmds"] += "/generator/offset 0 0 -4.975 um\n"
			else:
				self.settings["extra_post_cmds"] += "/generator/offset 0 0 4.975 um\n"
			self.settings["evtsrc"] = "In114m"
			
		# special setting for calibration source generators
		if self.settings["evtsrc"][:2] in ["Bi","Ce","Sn","Cd","In","Sr","Cs"]:
			srd = "1.5 mm"
			if self.settings["evtsrc"][:2] in ["Ce","Cd"]:
				srd = "1.25 mm"
			self.settings["extra_post_cmds"] += "/generator/sourceRadius %s\n"%srd
			self.needsHolder = True
			self.settings["evtsrc"] += "_o_n"
		
		# special settings for xenon source generators
		if self.settings["evtsrc"][:2] == "Xe":
			self.settings["evtsrc"] += "_g_n"
			
		# set location of source holder
		if not self.settings["sourceholderpos"]:
			if self.needsHolder:
				self.settings["sourceholderpos"] = "0 0	0 m"
				self.settings["extra_post_cmds"] += "/generator/relholder true\n"
			else:
				self.settings["sourceholderpos"] = "0 0.5 0 m"


	def set_dirs(self):
		self.g4_workdir = os.environ["G4WORKDIR"]
		self.g4_bindir = os.environ["G4BINDIR"]
		self.g4_evtsdir = os.environ["G4EVTDIR"]+"/"+self.settings["evtsrc"]
		
		#self.type_dir = self.settings["simName"]+"_"+self.settings["evtsrc_name"]
		self.type_dir = self.settings["evtsrc_name"]
		if "gunenergy" in self.settings:
			self.type_dir += "_%.1fkeV"%self.settings["gunenergy"]
		
		if not self.g4_out_dir_base:
			self.g4_out_dir_base = self.g4_workdir+"/output/"
		self.g4_out_dir = self.g4_out_dir_base+"/%s/"%self.type_dir
		self.g4_log_dir = self.g4_workdir+"/logs/%s/"%self.type_dir
		self.g4_macro_dir = self.g4_workdir+"/macros/%s/"%self.type_dir
		self.g4_out_name = "%s/g4_run_%%s.root"%self.g4_out_dir

	def set_detector_offsets(self):
		betaSim.settings["extra_cmds"] += "/detector/rotation 0.037\n"
		betaSim.settings["extra_cmds"] += "/detector/offset -3.98 0.44 0 mm\n"

	def launch_sims(self, maxIn=100000, hours_old=0):
		
		self.set_dirs()
		parallel_jobfile = "%s/jobs.txt"%self.g4_macro_dir
		
		os.system("mkdir -p %s"%self.g4_macro_dir)
		os.system("mkdir -p %s"%self.g4_out_dir)
		os.system("mkdir -p %s"%self.g4_log_dir)
		
		inflist = [f for f in os.listdir(self.g4_evtsdir) if f[:5]=="Evts_"]
		inflist.sort()
		inflist = inflist[:maxIn]
		nruns = len(inflist)
		
		oldtime = time.time() - hours_old*3600
		
		# main simulations
		os.system("rm -r %s/*"%self.g4_macro_dir)
		jobsout = open(parallel_jobfile,"w")
		ucnG4_prod = self.g4_bindir+"/ucnG4_prod"
		onejob = ""
		
		# set up macros for each job
		for rn in range(nruns):
			self.settings["run_num"] += 1
			self.settings["jobname"] = self.settings["simName"]+"_%i"%self.settings["run_num"]
			self.settings["outfile"]=self.g4_out_name%str(self.settings["run_num"])
			self.settings["evtfile"]=self.g4_evtsdir+"/"+inflist[rn]

			#self.settings["nevt"] = 1000 # assume this many events per input file... TODO something more elegant

			self.settings["nevt"] = 1000000000	# actual limit is number of events in file, nominally 10000

			self.settings["joblog"] = "%s/gen_macro_%i.txt"%(self.g4_log_dir,self.settings["run_num"])
			g4_sub_file = "%s/geantjob_%i.sub"%(self.g4_macro_dir,self.settings["run_num"])
			
			# source position scan
			if self.settings["sourceScan"]:
				xpos = (((rn*self.anagroup + rn/(nruns/self.anagroup))%nruns)/float(nruns-1)-0.5)*self.settings["sourceScan"]
				if self.settings["sourceholderpos"] != "0 0.5 0 m":
					self.settings["sourceholderpos"] = "%g 0 0 mm"%xpos
						
			# skip recently-run jobs
			if os.path.exists(self.g4_out_name%str(self.settings["run_num"])) and os.stat(self.g4_out_name%str(self.settings["run_num"])).st_mtime > oldtime:
				continue;
			
			# generate macro file
			open(os.path.expanduser("%s/geantgen_%i.mac"%(self.g4_macro_dir,self.settings["run_num"])),"w").write(open("GeantGenMacroTemplate.mac","r").read()%self.settings)
			# single job execution command, appended to batch job file
			onejob = ucnG4_prod + " %s/geantgen_%i.mac"%(self.g4_macro_dir,self.settings["run_num"])
			jobsout.write(onejob+" > %s 2>&1\n"%self.settings["joblog"])
		
		jobsout.close()
		
		print "Running simulation jobs..."
		os.system("cat "+parallel_jobfile)
		if nruns > 1:
			os.system("nice -n 20 parallel -P 6 < %s"%parallel_jobfile)
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
			for f in anafiles[:self.anagroup]:
				fout.write(f[1]+"\n")
			fout.close()
			anafiles = anafiles[self.anagroup:]
			print "\n----- %s ------"%outlist_name
			os.system("cat "+outlist_name)
			analyzer_bin = self.g4_bindir+"/"+self.settings["analyzer"]
			if nMin <= nanalyzed <= nMax:
				jobsout.write("%s %s %s/analyzed_%i.root %s\n"%(analyzer_bin,outlist_name,self.g4_out_dir,nanalyzed,self.settings["ana_args"]))
			nanalyzed += 1
		jobsout.close()
		print "\n----- %s ------"%resim_jobfile
		os.system("cat "+resim_jobfile)
		print
		os.system("nice -n 10 parallel -P 4 < %s"%resim_jobfile)
		os.system("rm %s/outlist_*.txt"%self.g4_out_dir)
		os.system("rm "+resim_jobfile)


if __name__ == "__main__":
	
	parser = OptionParser()
	parser.add_option("-k", "--kill", dest="kill", action="store_true", default=False, help="kill running replays")
	parser.add_option("--evtgen", dest="evtgen", action="store_true", default=False, help="run event generators")
	parser.add_option("--calsrcs", dest="calsrcs", action="store_true", default=False, help="simulate sealed calibration sources")
	parser.add_option("--xesrcs", dest="xesrcs", action="store_true", default=False, help="simulate xenon sources")
	parser.add_option("--sim", dest="sim", action="store_true", default=False, help="run initial Geant4 simulation")
	parser.add_option("--ana", dest="ana", action="store_true", default=False, help="run post-analyzer")
	options, args = parser.parse_args()
	
	if options.kill:
		os.system("killall -9 parallel")
		os.system("killall -9 ucnG4_prod")
		os.system("killall -9 UCNA_MC_Analyzer")
		os.system("killall -9 GeantSimManager.py")
		exit(0)
	
	######################
	# calibration sources.
	# sources ["Bi207","Sn113","Ce139","Cd109","Cs137","In114E","In114W","Cd113m"] 1e6 each
	######################
	if options.calsrcs:
		
		if options.evtgen:
			for g in ["Bi207","Sn113","Ce139","Cd109","Cs137","In114m","Cd113m"]:
				os.system("rm -rf %s/%s_o_n"%(os.environ["G4EVTDIR"],g))
				os.system("../../MC_EventGen run %s %s o n 10000 100 x"%(g,os.environ["G4EVTDIR"]))

		for g in ["Bi207","Sn113","Ce139","Cd109","Cs137","In114E","In114W"]: #,"Cd113m"]:
			sourceSim = GeantSimManager("thinfoil",
										fmap="/home/mmendenhall/UCNA/Aux/Fieldmap_20101028_b.txt", # sources simulatable with full field map
										geometry="thinFoil")
			sourceSim.settings["sourceScan"] = 80.	# move source holder across detector by 80mm
			sourceSim.set_evtsrc(g)
			if options.sim:
				sourceSim.launch_sims(maxIn=100)
			if options.ana:
				sourceSim.launch_postanalyzer()




	####################				
	# Xenon (Xe135_3-2+ is 915keV beta endpoint; simulate extra events)
	####################
	XeIsots =  [	"Xe135_3-2+","Xe133_3-2+",
					"Xe129_11-2-","Xe131_11-2-","Xe133_11-2-",
					"Xe135_11-2-","Xe137_7-2-","Xe127_1-2+","Xe125_1-2+" ]
	if options.xesrcs:
		
		if options.evtgen:
			for g in XeIsots:
				os.system("rm -rf %s/%s_g_n"%(os.environ["G4EVTDIR"],g))
				os.system("../../MC_EventGen run %s %s g n 10000 1000 x"%(g,os.environ["G4EVTDIR"]))

		for g in XeIsots:
			sourceSim = GeantSimManager("thinfoil", geometry="thinFoil")
			sourceSim.set_evtsrc(g)
			if options.sim:
				maxIn = {"Xe135_3-2+": 300}.get(g,100)
				sourceSim.launch_sims(maxIn=maxIn)
			if options.ana:
				sourceSim.launch_postanalyzer()





if False:

	# self.settings["ana_args"] += " saveall"
	
	####################				
	# neutrons: 2010
	####################
	
	# unpolarized beta baseline: 5e7 in 520 clusters
	if 0:
		betaSim = GeantSimManager("20120823")
		#betaSim.set_generator("neutronBetaUnpol")
		betaSim.set_evtsrc("neutronBetaUnpol")
		betaSim.set_detector_offsets()
		betaSim.g4_out_dir_base = "/data2/mmendenhall/G4Out/2010/"
		### DON'T RE-SIMULATE ### betaSim.launch_sims(nEvents=5e7,nClusters=520,hours_old=10*24)
		betaSim.launch_postanalyzer()
		exit(0)
		
	# beta decay in magnetic field wiggles, 1e-3 vacuum: 1e7 in 104 clusters
	if 0:
		betaSim = GeantSimManager("20120824_MagF",vacuum="1.e-3 torr",fmap="/home/mmendenhall/UCNA/Aux/Fieldmap_20101028_b.txt")
		betaSim.set_generator("neutronBetaUnpol")
		betaSim.settings["extra_cmds"] += "/detector/rotation 0.037\n"
		betaSim.settings["extra_cmds"] += "/detector/offset -3.98 0.44 0 mm\n"
		#betaSim.launch_sims(nEvents=1e7,nClusters=104,hours_old=0)
		betaSim.launch_postanalyzer()
	
	# thin foils
	if 0:
		betaSim = GeantSimManager("endcap_180_150",geometry="thinFoil")
		betaSim.set_generator("neutronBetaUnpol")
		betaSim.settings["extra_cmds"] += "/detector/rotation 0.037\n"
		betaSim.settings["extra_cmds"] += "/detector/offset -3.98 0.44 0 mm\n"
		#betaSim.launch_sims(nEvents=5e7,nClusters=520,hours_old=100*24)
		betaSim.launch_postanalyzer()


	# 2010 sources, **postanalyzer only**
	if 0:
		for g in ["Bi207","Sn113","Ce139"]:
			sourceSim = GeantSimManager("20120823",fmap="/home/mmendenhall/UCNA/Aux/Fieldmap_20101028_b.txt",geometry="C")
			sourceSim.g4_out_dir_base = "/data2/mmendenhall/G4Out/2010/"
			sourceSim.set_evtsrc(g)
			sourceSim.launch_postanalyzer()
		exit(0)
		
	# sources ["Bi207","Sn113","Ce139","Cd109","Cs137","In114E","In114W","Cd113m"] 1e6 each
	if 0:
		for g in ["In114E","In114W"]:
			sourceSim = GeantSimManager("thinfoil",fmap="/home/mmendenhall/UCNA/Aux/Fieldmap_20101028_b.txt",geometry="thinFoil")
			sourceSim.settings["physlist"]="livermore"
			sourceSim.settings["sourceScan"]=80.
			sourceSim.settings["extra_cmds"] += "/detector/sourcefoilthick 9.5 um\n"
			sourceSim.settings["extra_cmds"] += "/detector/MWPCBowing 5 mm\n"
			sourceSim.settings["ana_args"] += " undead"
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
	# Xenon [	"Xe135_3-2+","Xe133_3-2+",
	#			"Xe129_11-2-","Xe131_11-2-","Xe133_11-2-",
	#			"Xe135_11-2-","Xe137_7-2-","Xe127_1-2+","Xe125_1-2+"	]
	# 3M for most; do lots more for important Xe135_3-2+
	####################
	if 1:
		start = time.time()
		for g in [ "Xe135_3-2+","Xe133_3-2+",
				"Xe129_11-2-","Xe131_11-2-","Xe133_11-2-",
				"Xe135_11-2-","Xe137_7-2-","Xe127_1-2+","Xe125_1-2+" ]:
			sourceSim = GeantSimManager("reanalyze",geometry="2011/2012", vacuum="1.e-3 torr")
			sourceSim.settings["extra_cmds"] += "/detector/MWPCBowing 5 mm\n"
			sourceSim.settings["ana_args"] += " saveall"
			sourceSim.settings["ana_args"] += " cathodes"
			sourceSim.settings["physlist"]="livermore"
			sourceSim.set_evtsrc(g)
			#sourceSim.launch_sims(maxIn=50000000, hours_old=0)
			sourceSim.launch_postanalyzer()
		Time = time.time()-start
		print time.time()-start
	
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


	####################				
	# isotropic line
	####################
	if 0:
		for l in [56.7,]:
			iline = GeantSimManager("IsotLine")
			iline.set_generator("eGunRandMomentum")
			iline.settings["positioner"] = "Fixed"
			iline.settings["gunenergy"] = l
			iline.settings["ana_args"] += " saveall"
			iline.settings["sourceholderpos"] = "0 0 0 m"
			iline.launch_sims(nEvents=1e6,nClusters=36,hours_old=0)
			iline.launch_postanalyzer()


