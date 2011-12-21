#!/usr/bin/python

import os
import numpy
from numpy import *
import numpy.linalg as linalg
#from Histogram import *
#from LinFitter import *
from QFile import *

sourceNames = { 1:"Sn113", 2:"Bi207", 3:"Sr85", 4:"Cd109", 7:"Sr90", 8:"In114", 9:"Ce139", 10:"Monosmuth", 11:"500keV" }
fancySourceNames = { 1:"$^{113}$Sn", 2:"$^{207}$Bi", 3:"$^{85}$Sr", 4:"$^{109}$Cd", 7:"$^{90}$Sr", 8:"$^{114}$In", 9:"$^{139}$Ce", 10:"Monosmuth", 11:"500keV" }		
peakNames = { 8:"Bi207 1", 9:"Bi207 2", 11:"Sn113", 12:"Sr85", 13:"Cd109", 14:"In114", 15:"Ce139", 16:"Monosmuth", 17:"500keV" }  
betaSensorNames = { 'E':["ADCE3Beta","ADCE4Beta","ADCE1Beta","ADCE2Beta","ADCRefELED"], 'W':["ADCW1Beta","ADCW2Beta","ADCW3Beta","ADCW4Beta","ADCRefWLED"] }
class PeakNum:
	Sn = 11
	Ce = 15
	Bi_1 = 8
	Bi_2 = 9
	
#'official analyzer' tube numbering	
jTubes = {'E':[3,4,1,2,5],'W':[1,2,3,4,5]}

# listing of all PMTs
allPMTs = [ [(s,t) for t in range(4)] for s in ['E','W'] ]
allPMTs = allPMTs[0]+allPMTs[1]

# various sets of runs
betagroups = [	(13905,13964), (14077,14100), (14127,14261), (14356,14380), (14397,14421), (14432,14467),
				(14535,14667), (14688,14782), (14888,14994), (15084,15150), (15172,15356), (15448,15639), (15667,15915),
				(15943,15966), (16097,16216) ]
sourcegroups = [(13883,13894), (14104,14118), (14383,14394), (14516,14530), (14736,14746), (15645,15662), (15916,15939),
				#(15357,15371),	# badly tilted sources
				(16240,16257),	# 2010 end of year sources
				(17233,17249),
				(17359,17387),
				(17517,17527),
				(17871,17922),
				(17925,17956),	# long centered runs
				(18020,18055),	# old and new Cd sources
				(18357,18386)	# new In source
				]
thecalruns = [13890,14111,14390,14524,14743,15653,15931]
xegroups = [	(14264,14347), (15991,16077), (16983,17078),
				(17224,17230),	# old Xenon
				(17561,17734) ]
				
def expandRanges(rgs,rmin=0,rmax=100000):
	rns = []
	for rg in rgs:
		rns += [r for r in range(rg[0],rg[1]+1) if rmin <= r <= rmax]
	return rns
	
	
# electron mass, keV/c^2
m_e = 511.0
# beta, for electron by default
electron_beta = (lambda ke, m = m_e: sqrt(ke**2+2*m*abs(ke))/(m+abs(ke)) )
# neutron decay beta spectrum
def betaSpectrumUnpol(ke):
	betaEP = 782.347 
	if ke < 0 or ke > betaEP:
		return 0
	return  (ke+m_e)*sqrt(ke*ke+2*m_e*ke) * (betaEP-ke)**2
	
# get geometry for run number
def geomName(runNum):
	gn = "X"
	if runNum >= 7000:
		gn = "A"
	if runNum >= 9220:
		gn = "B"
	if runNum >= 10400:
		gn = "C"
	if runNum >= 11390:
		gn = "D'"
	if gn == "X":
		exit(1)
	return gn
		
# beta run range for each geometry
geomRanges = {'A':(7662,9189),'B':(9220,10375),'C':(10404,11110),"C'":(12425,12721),
				0:(13905,16216),1:(13905,13964),2:(14077,14380),3:(14397,14667),4:(14535,14667),5:(14688,15356),6:(15448,15915),7:(15943,16216)}
	
# source from KVM
class Source:
	def __init__(self,kvm):
		self.runNum = 0
		self.side = kvm.getFirst("side")
		if not self.side:
			self.side = 'N'
		self.name = kvm.getFirst("name")
		self.x = kvm.getFirstF("x")
		self.y = kvm.getFirstF("y")
		self.wx = kvm.getFirstF("wx")
		self.wy = kvm.getFirstF("wy")
		self.type = kvm.getFirst("type")
		self.sID = int(kvm.getFirst("sID"))
		self.etrue = kvm.getFirstF("Etrue")
		self.eta = [kvm.getFirstF("eta_%i"%t) for t in range(4)]
				
	def __repr__(self):
		return "<Source '%s' [%s:%i] at (%.1f,%.1f)>"%(self.name, self.side, self.sID, self.x,self.y)	
		
# source spectrum line from KVM
class SourceLine:
	def __init__(self,kvm,nm):
		self.runNum = 0
		self.name = nm
		self.tube = 4
		if nm[-2] == '_':
			self.tube = int(nm[-1]);
		self.center = 	kvm.getFirstF("center")
		self.dcenter = 	kvm.getFirstF("dcenter")
		self.width = 	kvm.getFirstF("width")
		self.encenter = kvm.getFirstF("energy")
		self.dencenter = kvm.getFirstF("denergy")
		self.enwidth = 	kvm.getFirstF("energywidth")
		self.height = 	kvm.getFirstF("height")
		self.area = 	self.height*self.enwidth*sqrt(2*pi)
		self.integral = kvm.getFirstF("integral","0")
		self.nPE = 		kvm.getFirstF("nPE","0")
		self.sID = 		int(kvm.getFirst("sID"))
		self.type = 	int(kvm.getFirst("type"))
		self.x = 		kvm.getFirstF("x")
		self.y = 		kvm.getFirstF("y")
		self.wx = self.wy = 0
		self.simulated = not (kvm.getFirst("simulated") == "no")
		self.side = ''
	
	def uid(self):
		return (self.runNum,self.sID,self.tube,self.type)
	
	def radius(self):
		return sqrt(self.x**2+self.y**2)
		
	def __repr__(self):
			return "<SourceLine for %i at (%.1f,%.1f) = %.1f>"%(self.sID,self.x,self.y,self.center)

class ManualLine:
	def __init__(self,s,t,type,adc,evisExp):
		self.side = s
		self.tube = t
		self.type = type
		self.center = adc
		self.evisExp = evisExp
		self.eta = 1.0
		self.gms0 = 1.0
		self.runNum = 0
		self.sID = 0

# result of a kurie plot fit
class kurieFit(KVMap):
	
	def __init__(self,kvm):
		self.dat = kvm.dat
		self.side = kvm.getFirst("side")
		self.tube = kvm.getFirstF("tube")
		self.ep = kvm.getFirstF("endpoint")
		self.dep = kvm.getFirstF("dendpoint")
		self.simep = kvm.getFirstF("simendp")
		self.fitrange = (kvm.getFirstF("fitStart"),kvm.getFirstF("fitEnd"))		
		self.afp = kvm.getFirst("afp")
		self.type = kvm.getFirst("type")
		
# muon peak fit
class muonFit:
	
	def __init__(self,kvm):
		
		self.side = kvm.getFirst("side")
		self.tube = int(kvm.getFirstF("tube"))
		self.adc = kvm.getFirstF("peakADC")
		self.dadc = kvm.getFirstF("dADC")
		self.fiterr = int(kvm.getFirstF("fiterr"))
		self.n = int(kvm.getFirstF("segment"))
		self.sigma = kvm.getFirstF("peakSigma")
		self.dsigma = kvm.getFirstF("dSigma")
		self.time = kvm.getFirstF("time")
		self.dtime = kvm.getFirstF("dtime")
		self.counts = int(kvm.getFirstF("counts"))
		
		
# gms info data
class gmsInfo:
	
	def __init__(self,kvm):
		
		self.run = kvm.getFirstF("run")
		self.tstart = kvm.getFirstF("tstart")
		
		self.gms0 = {}
		self.gms = {}
		self.ledadc = {}
		self.gmsfact_Muon = {}
		self.gmsfact_LED = {}
		self.gmsfact_Kurie = {}
		self.ledbright = {}
		self.deltaL = {}
		
		for s in ['E','W']:
			self.gms0[s] = []
			self.gms[s] = []
			self.ledadc[s] = []
			self.gmsfact_LED[s] = []
			self.gmsfact_Muon[s] = []
			self.gmsfact_Kurie[s] = []
			self.deltaL[s] = []
			self.ledbright[s] = kvm.getFirstF(s+"_ledbright")
			
			for t in range(4):
				self.gms0[s].append(kvm.getFirstF("%s%i_gms0"%(s,t)))
				self.gms[s].append(kvm.getFirstF("%s%i_gms"%(s,t)))
				self.ledadc[s].append(kvm.getFirstF("%s%i_ledadc"%(s,t)))
				self.gmsfact_Muon[s].append(kvm.getFirstF("%s%i_gmsfactor_Muon"%(s,t)))
				self.gmsfact_LED[s].append(kvm.getFirstF("%s%i_gmsfactor_LED"%(s,t)))
				self.gmsfact_Kurie[s].append(kvm.getFirstF("%s%i_gmsfactor_Kurie"%(s,t)))
				self.deltaL[s].append(kvm.getFirstF("%s%i_deltaL"%(s,t)))
				
			self.ledadc[s].append(kvm.getFirstF(s+"_ref_ledadc"))
			
			
			
# read in linearity correction data
class LinearityDB:
	def __init__(self):
		pass


class betaTube:
	def __init__(self,s,t,kvm=None):
		self.side = s
		self.tube = t
		self.load_kvm(kvm)
	
	def load_kvm(self,kvm):
		if not kvm:
			return
		self.gms0 = kvm.getFirstF("gms0_%i"%self.tube)
		self.gmsRel = kvm.getFirstF("gmsRel_%i"%self.tube)
		self.pmt_thresh = [float(x) for x in kvm.getFirst("trig_effic_%i"%self.tube,"0,0").split(',')]

class counter:
	def __init__(self,kvm):
		self.counts = {}
		k = kvm.dat.keys()
		k.sort()
		assert len(k)%5 == 0
		for i in range(len(k)/5):
			s0 = k[5*i]
			self.counts[(s0,'T')]=kvm.getFirstF(s0)
			for s in ['E','W','B','N']:
				self.counts[(s0,s)]=kvm.getFirstF(s0+"_"+s)
				
			
# RunData QFile reader
class RunData(QFile):

	def __init__(self,rn,basePath="../RunData/"):
		QFile.__init__(self,basePath+"/Run_%i.txt"%rn)
		
		self.runNum = int(self.getItem("RunInfo","Run_Number","0"))
		self.afpOn = bool(self.getItem("RunInfo","AFP_On","0"))
		self.vLed = [float(self.getItem("RunInfo","V_LED_E","0")),float(self.getItem("RunInfo","V_LED_W","0"))]
		self.scs = float(self.getItem("RunInfo","SCS_Field","0"))
		self.role = self.getItem("RunInfo","Role_Name")
		self.pmt_thresh = {}
		self.tubes = {}
		if basePath=="../RunData/":
			self.runTime = self.getItemF("Trigger","runTime",1.0)
			self.nEvents = self.getItemF("Trigger","nEvents",0.0)
			z = self.getFirst("EventCounts")
			if z:
				self.cntr = counter(z)
		else:
			R0 = RunData(rn)
			self.runTime = R0.runTime
			self.nEvents = R0.nEvents
		self.rate = self.nEvents/(self.runTime+0.00001)
		for s in ['E','W']:
			self.tubes[s] = []
			for t in range(4):
				self.tubes[s].append( betaTube(s,t,self.getFirst("BetaSc"+s)) )				
		self.gmsMatchedRun = {'E':int(self.getItem("BetaScE","gmsMatchedRun",str(rn))),'W':int(self.getItem("BetaScW","gmsMatchedRun",str(rn)))}
		self.runTime = float(self.getItem("Trigger","runTime","0"))
	
	# total rate (after beam cuts)
	def totalRate(self):
		return (self.cntr.counts[("Total","T")]-self.cntr.counts[("Crud","T")])/self.runTime
	# 'correct' events rate
	def betaRate(self):
		return (self.cntr.counts[("Correct","T")]+self.cntr.counts[("TypeI","T")]+self.cntr.counts[("TypeII-III","T")])/self.runTime
	# monitor events rate
	def monRate(self):
		return self.getItemF("EventCounts","UCNMon",0)/self.runTime
	# gms events rate
	def gmsRate(self):
		return self.getItemF("EventCounts","GMS",0)/self.runTime
	# vetoed events rate (muons + gammas)
	def vetoedRate(self):
		return (self.getItemF("EventCounts","Muon")+self.getItemF("EventCounts","Unknown"))/self.runTime
	# get run time, possibly accounting for deadtime
	def getRunTime(self,fixdt=False):
		if fixdt:
			return self.runTime - (self.cntr.counts[("Total","T")]-self.cntr.counts[("Crud","T")])*(12.0e-6)
		return self.runTime
		
	# check whether this is a self-calibrated reference run
	def isReference(self):
		if self.role != "SourcesCal":
			return False
		for t in range(4):
			if self.gmsMatching['E'][t] != 1.0 or self.gmsMatching['W'][t] != 1.0:
				return False
		return True
	
	# return list of sources
	def getSources(self):
		sources = self.getMatching("object","Source")
		l = []
		for k in sources.dat.keys():
			for j in sources.dat[k]:
				l.append(Source(j));
				l[-1].runNum = self.runNum
		return l
		
	# return list of source lines
	def getSourceLines(self):
		sourcelines = self.getMatching("object","SpectrumPeak")
		sources = {}
		for s in self.getSources():
			sources[s.sID] = s
		l = []
		for k0 in sourcelines.dat:
			for k in sourcelines.dat[k0]:
				if not int(k.getFirst("sID","0")):
					continue
				
				l.append(SourceLine(k,k0))
				l[-1].runNum = self.runNum

				mysid = 1*(l[-1].sID)
				mysrc = sources[mysid]
				
				l[-1].x = mysrc.x
				l[-1].y = mysrc.y
				l[-1].wx = mysrc.wx
				l[-1].wy = mysrc.wy
				l[-1].side = mysrc.side
				l[-1].runRate = self.rate
				if l[-1].tube < 4:
					if mysrc.eta[l[-1].tube]:
						l[-1].eta = mysrc.eta[l[-1].tube]
					else:
						l[-1].eta = 0
				else:
					if mysrc.eta[0]:
						l[-1].eta = 0.25*sum(mysrc.eta)
					else:
						l[-1].eta = 0
				if l[-1].tube < 4 and l[-1].side in ['E','W']:
					l[-1].gmsRel = self.tubes[l[-1].side][l[-1].tube].gmsRel
					l[-1].gms0 = self.tubes[l[-1].side][l[-1].tube].gms0
				if l[-1].center < 0:
					l[-1].center = 0
					
				if not l[-1].eta:
					print "*** Warning: missing eta for",l[-1].uid()
					
		return l
		
# reader for RunsDB info				
class RunsDB:
	
	def __init__(self,fin="../SummaryData/UCNA Run Log.txt",rmin=0,rmax=1000000,basepath="../RunData"):
		self.runs = {}
		for l in [l.strip() for l in open(fin).readlines() if l[0] in ['*','@']]:
			if l[0] == '*':
				#try:
					rn = int(l.split()[0][1:])
					if rn < rmin or rn > rmax:
						continue
					if os.path.exists(basepath+"/Run_%i.txt"%rn):
						self.runs[rn] = RunData(rn,basepath)
				#except:
				#	print "Bad line or unable to process RunData file:",l
					
		print "Loaded RunsDB with data for %i runs."%len(self.runs)

# load a histogram from a stringmap
def hFromStringmap(m):
	nbins = m.getFirstF("nbins")
	
	h = histogram()
				
	
