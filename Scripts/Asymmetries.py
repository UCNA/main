#!/usr/bin/python

import os
from LinFitter import *
from PyxUtils import *
from QFile import *
from EncalDB import *
from Histogram import *
try:
	from scipy import stats
except:
	stats = None
	
# runs in pulse pair, half octet, and octet groups			  		  
#ppSegments = [["A1","A2","A3","A4","A5","A6"],["A7","A8","A9","A10","A11","A12"],["B1","B2","B3","B4","B5","B6"],["B7","B8","B9","B10","B11","B12"]]
ppSegments = [["A1","A2","A4","A5"],["A7","A9","A10","A12"],["B1","B2","B4","B5"],["B7","B9","B10","B12"]]
hoSegments = [ppSegments[0]+ppSegments[1],ppSegments[2]+ppSegments[3]]
octSegments = [hoSegments[0]+hoSegments[1],]
divSegments = [octSegments,hoSegments,ppSegments]
unitNames = {0:"Octet",1:"Half Octet",2:"Pulse Pair",3:"Run"}

scols = {"East":rgb.red,"West":rgb.blue}
afpSymbs = {"On":symbol.circle,"Off":symbol.triangle}
afpLines = {"Off":style.linestyle.dashed,"On":style.linestyle.dotted}


def addQuad(a,b):
	return sqrt(a**2+b**2)
	
class asymmetryFit(KVMap):
	def __init__(self,m=KVMap()):
		self.dat = m.dat
		self.loadFloats(["fitMin","fitMax","A0_fit","dA0"])
		self.A0 = self.A0_fit
	def __repr__(self):
		return "<Asym: %.2f +- %.2f, (%g,%g)>"%(self.A0,self.dA0,self.fitMin,self.fitMax)

class instAsym(KVMap):
	def __init__(self,m=KVMap()):
		self.dat = m.dat
		self.loadFloats(["fitMin","fitMax","IA","dIA"])
	def __repr__(self):
		return "<Asym: %.2f +- %.2f, (%g,%g)>"%(self.A0,self.dA0,self.fitMin,self.fitMax)


class eventRate(KVMap):
	def __init__(self,m):
		self.dat = m.dat
		self.loadFloats(["counts","d_counts","rate","d_rate"])
		self.loadStrings(["side","afp","fg","name"])
	def dRate(self):
		if not self.counts:
			return 0
		return self.rate/sqrt(self.counts)
		
class kurieFit(KVMap):
	def __init__(self,m):
		self.dat = m.dat
		self.loadFloats(["endpoint","dendpoint","tube"])
		self.loadStrings(["side","afp","type"])
		self.ep = self.endpoint
		self.dep = self.dendpoint

class runCal(KVMap):
	def __init__(self,m):
		self.dat = m.dat
		self.run = int(self.getFirstF("run"))
		self.eOrig={}
		self.eFinal={}
		self.gms={}
		self.gms0={}
		for s in ['E','W']:
			for t in range(5):
				self.eOrig[(s,t)]=self.getFirstF("%s%i_eOrig"%(s,t))
				self.eFinal[(s,t)]=self.getFirstF("%s%i_eFinal"%(s,t))
				if t<4:
					self.gms[(s,t)]=self.getFirstF("%s%i_gms"%(s,t))
					self.gms0[(s,t)]=self.getFirstF("%s%i_gms0"%(s,t))
					
	def getGMSTweak(self,side,tube):
		return self.eFinal[(side,tube)]/self.eOrig[(side,tube)]

class bgSubtr(KVMap):
	def __init__(self,m):
		self.dat = m.dat
		self.loadFloats(["nBG","d_nBG","xs","d_xs","eMin","eMax"])
		self.loadStrings(["side","afp","name"])

class asymData:
	"""Asymmetry data extracted from a group of runs"""
	def __init__(self):
		self.octruns = {}
		self.asyms = []
		
	def getRuns(self):
		rns = []
		for k in self.octruns:
			rns += self.octruns[k]
		rns.sort()
		return rns

	def whichSegment(self,divlev):
		sn = 0
		smax = 0
		ntot = len(self.octruns)
		for (n,sg) in enumerate(divSegments[divlev]):
			sc = 0
			for k in self.octruns:
				if k in sg:
					sc += len(self.octruns[k])
			if sc > smax:
				smax = sc
				sn = n
		print self.octruns.keys(),smax,"of",ntot,"in div",sn
		return sn
		
	def addRun(self,rtype,rn):
		self.octruns.setdefault(rtype,[]).append(int(rn))
	
	def getAsym(self,fitMin=225,fitMax=675):
		matchAsyms = [a for a in self.asyms if a.fitMin==fitMin and a.fitMax==fitMax]
		if len(matchAsyms)==1:
			return matchAsyms[0]
		assert not len(matchAsyms)
		return None
	
		
		
class AsymmetryFile(QFile,asymData):

	def __init__(self,fname):
		QFile.__init__(self,fname)
		asymData.__init__(self)
		self.asyms = [asymmetryFit(m) for m in self.dat.get("asymmetry",[])]
		self.kuries = [kurieFit(m) for m in self.dat.get("kurieFit",[])]
		self.rates = [eventRate(m) for m in self.dat.get("rate",[])]
		self.bgs = dict([((b.side,b.afp,b.type),b) for b in [bgSubtr(m) for m in self.dat.get("bg_subtr_fit",[])]])
		for m in self.dat.get("Octet",[]):
			for k in m.dat:
				if k in octSegments[0]:
					for rns in m.dat[k]:
						for rn in rns.split(","):
							self.addRun(k,rn)
		if not self.getRuns():
			rns = self.fname.split("/")[-1].split("_")[0].split("-")
			if len(rns)==2:
				self.addRun("start",int(rns[0]))
				self.addRun("end",int(rns[1]))
		self.runcals = dict([(r.run,r) for r in [runCal(m) for m in self.dat.get("runcal",[])]])
		
	def kepDelta(self,side,type="0"):
		"""Calculate difference between On/Off endpoint"""
		kOn = self.getKurie(side,"On",type)
		kOff = self.getKurie(side,"Off",type)
		return [kOn.ep-kOff.ep,addQuad(kOn.dep,kOff.dep)]
	
	def getKurie(self,side,afp,type,tube=4):
		for k in self.kuries:
			if k.side==side and k.afp==afp and k.type==type and k.tube==tube:
				return k
		print "*** Can't find kurie in",side,afp,type,self.fname
		return None
	
	def getRate(self,side,afp,fg,name):
		hname = name
		if hname[-1] == "_":
			hname = "%s%s_%s"%(name,side,afp)
		for r in self.rates:
			if r.side==side and r.afp==afp and r.fg==fg and r.name==hname:
				return r
		print "Missing rate for",hname,fg
		return None
	
	def getGMSTweak(self,side,tube):
		tweaks = [self.runcals[r].getGMSTweak(side,tube) for r in self.runcals]
		return sum(tweaks)/len(tweaks)
		
	def getInstAsym(self):
		return instAsym(self.dat.get("instasym",[KVMap()])[0])
				
def collectAsymmetries(basedir,depth):
	print "Collecting asymmetry data from",basedir,"at depth",depth
	asyms = []
	fls = os.listdir(basedir)
	fls.sort()
	for f in [ (basedir+'/'+f+'/',basedir+'/'+f+'/'+f+'.txt') for f in fls if os.path.isdir(basedir+'/'+f) and '-' in f and f[:3]!="139"]:
		if not depth and os.path.exists(f[1]):
			asyms.append(AsymmetryFile(f[1]))
		else:
			asyms += collectAsymmetries(f[0],depth-1)
	return asyms
	
def plot_octet_asymmetries(basedir,depth=0):
	
	##############
	# collect data
	##############
	gdat = []
	iadat = []
	ndat = []
	bgRateDat = {'E':{"Off":[],"On":[]},'W':{"Off":[],"On":[]}}
	muRateDat = {'E':{"Off":[],"On":[]},'W':{"Off":[],"On":[]}}
	kdat = {'E':{"Off":[],"On":[]},'W':{"Off":[],"On":[]}}
	kepDelta = {'E':[],'W':[],'C':[]}
	n=0
	print "--------------------- Division",depth,"-------------------------"
	for af in collectAsymmetries(basedir,depth):
		# asymmetry data
		a = af.getAsym()
		if a:
			if abs(a.A0+0.115) > 0.05 or abs(a.A0+0.115)/a.dA0 > 3:
				print "*** Funny asym",a.A0,"+-",a.dA0,"in",af.fname
			gdat.append([n,a.A0,a.dA0,af.getRuns()[0]])
		else:
			print "**** Missing asym in",af.fname
		ia = af.getInstAsym()
		iadat.append([n,ia.IA,ia.dIA])
		ndat.append([n,0])
		if len(ndat)>1:
			ndat[-1][-1] += ndat[-2][-1]
		# kurie data
		for s in kdat:
			kepDelta[s].append([n,]+af.kepDelta(s))
			if not (-15 < kepDelta[s][-1][1] < 15 and 0 < kepDelta[s][-1][2]):
				print "*** Funny delta",kepDelta[s][-1][1],kepDelta[s][-1][2],"in",af.fname
			for afp in kdat[s]:
				k = af.getKurie(s,afp,"0")
				if k:
					kdat[s][afp].append([n,k.ep,k.dep]+[af.getKurie(s,afp,"0",t) for t in range(4)])
					if not 782 < k.ep < 810:
						print "*** Funny Endpoint",k.ep,k.side,k.afp,"in",af.fname
				
				rt = af.getRate(s,afp,'1',"hEnergy_Type_0_%s_%s"%(s,afp))
				if rt:
					ndat[-1][-1] += rt.counts
				
				rt = af.getRate(s,afp,'0',"hEnergy_Type_0_%s_%s"%(s,afp))
				if rt:
					if not 0.10 < rt.rate < 0.30:
						print "*** Funny BG rate",rt.rate,rt.side,rt.afp,"in",af.fname
					if not 0.001 < rt.dRate():
						print "**** Funny dRate",rt.dRate(),rt.side,rt.afp,"in",af.fname
						continue
					bgRateDat[s][afp].append([n,rt.rate,rt.dRate()])
				else:
					print "*** Can't find rate for",n,s,afp
				
				murt = af.getRate(s,afp,'1',"hMuonNoSub_%s_%s"%(s,afp))
				if murt:
					if not 0.85 < murt.rate < 1.15:
						print "*** Funny muons",murt.rate,murt.side,murt.afp,"in",af.fname
					muRateDat[s][afp].append([n,murt.rate,murt.dRate()])
				else:
					print "*** Can't find muon rate for",n,s,afp

		kepDelta['C'].append([n,0.5*(kepDelta['E'][-1][1]-kepDelta['W'][-1][1]),0.5*addQuad(kepDelta['E'][-1][2],kepDelta['W'][-1][2])])
		n+=1
	print

	##############
	# set up graphs
	##############
	unitName=unitNames[depth]
	
	gAsyms=graph.graphxy(width=25,height=8,
				x=graph.axis.lin(title=unitName,min=0,max=gdat[-1][0]),
				y=graph.axis.lin(title="Asymmetry",min=-0.2,max=0),
				key = graph.key.key(pos="bl"))
	setTexrunner(gAsyms)

	gCount=graph.graphxy(width=15,height=15,
						 x=graph.axis.lin(title=unitName,min=0,max=gdat[-1][0]),
						 y=graph.axis.lin(title="Type 0 Counts",min=0),
						 key = None)
	setTexrunner(gCount)
			
	gBgRate=graph.graphxy(width=25,height=8,
				x=graph.axis.lin(title=unitName,min=0,max=gdat[-1][0]),
				y=graph.axis.lin(title="Background Rate [Hz]",min=0,max=0.35),
				key = graph.key.key(pos="bl"))
	setTexrunner(gBgRate)

	gMuRate=graph.graphxy(width=25,height=8,
						  x=graph.axis.lin(title=unitName,min=0,max=gdat[-1][0]),
						  y=graph.axis.lin(title="Muon Rate [Hz]"),
						  key = graph.key.key(pos="bl"))
	setTexrunner(gMuRate)
		
	gdEp=graph.graphxy(width=25,height=8,
				x=graph.axis.lin(title=unitName,min=0,max=gdat[-1][0]),
				y=graph.axis.lin(title="Endpoint On/Off Difference [keV]",min=-30,max=20),
				key = graph.key.key(pos="bl"))
	setTexrunner(gdEp)
	
	##############
	# event count
	##############
	gCount.plot(graph.data.points(ndat,x=1,y=2),[graph.style.line(lineattrs=[rgb.red,])])
	gCount.writetofile(basedir+"/OctetCounts_%i.pdf"%depth)
		
	##############
	# asymmetry plots
	##############
	gAsyms.plot(graph.data.points(gdat,x=1,y=2,dy=3,title=None),
				[graph.style.symbol(symbol.circle,size=0.2,symbolattrs=[rgb.red,]),
				graph.style.errorbar(errorbarattrs=[rgb.red,])])
			
	LF = LinearFitter(terms=[polyterm(0)])
	
	rbreak = 14800 # 15400 = magnet ramp, calibration changes
	gdat_A = [g for g in gdat if 14000 < g[3] < rbreak]
	gdat_B = [g for g in gdat if rbreak < g[3] ]
	
	if gdat_A:
		LF.fit(gdat_A,cols=(0,1,2),errorbarWeights=True)
		chi2 = LF.ssResids()
		ndf = len(gdat_A)-len(LF.coeffs)
		gtitle = "$A=%.5f\\pm%.5f$, $\\chi^2/\\nu = %.1f/%i$"%(LF.coeffs[0],1.0/sqrt(LF.sumWeights()),chi2,ndf)
		if stats:
			gtitle += " $(p=%.2f)$"%stats.chisqprob(chi2,ndf)
		gAsyms.plot(graph.data.points(LF.fitcurve(gdat_A[0][0],gdat_A[-1][0]),x=1,y=2,title=gtitle),[graph.style.line()])
	
	if gdat_B:
		LF.fit(gdat_B,cols=(0,1,2),errorbarWeights=True)
		chi2 = LF.ssResids()
		ndf = len(gdat_B)-len(LF.coeffs)
		gtitle = "$A=%.5f\\pm%.5f$, $\\chi^2/\\nu = %.1f/%i$"%(LF.coeffs[0],1.0/sqrt(LF.sumWeights()),chi2,ndf)
		if stats:
			gtitle += " $(p=%.2f)$"%stats.chisqprob(chi2,ndf)
		gAsyms.plot(graph.data.points(LF.fitcurve(gdat_B[0][0],gdat_B[-1][0]),x=1,y=2,title=gtitle),
					[graph.style.line([style.linestyle.dashed,])])
				
	gAsyms.writetofile(basedir+"/OctetAsym_%i.pdf"%depth)
	
	##############
	# instrumental asymmetry
	##############
	gIA=graph.graphxy(width=25,height=8,
				x=graph.axis.lin(title=unitName,min=0,max=iadat[-1][0]),
				y=graph.axis.lin(title="Instrumental Asymmetry",min=-0.04,max=0.01),
				key = graph.key.key(pos="bl"))
	setTexrunner(gIA)
	gIA.plot(graph.data.points(iadat,x=1,y=2,dy=3,title=None),
				[graph.style.symbol(symbol.circle,size=0.2,symbolattrs=[rgb.red,]),
				graph.style.errorbar(errorbarattrs=[rgb.red,])])
	gIA.writetofile(basedir+"/OctetInstAsym_%i.pdf"%depth)

	##############
	# endpoint, background, muons plots
	##############
	sideCols = {'E':rgb.red,'W':rgb.blue,'C':rgb(0.,0.7,0.)}
	for s in ['E','C','W']:
		
		LF.fit([kd for kd in kepDelta[s] if kd[2]>0],cols=(0,1,2),errorbarWeights=True)
		chi2 = LF.ssResids()
		ndf = len(gdat)-len(LF.coeffs)
		sname = "$%s_{\\rm on}-%s_{\\rm off}$"%(s,s)
		if s == 'C':
			sname = "$(E_{\\rm on}-E_{\\rm off}-W_{\\rm on}+W_{\\rm off})/2$"
		ch2str = "$\\chi^2/\\nu = %.1f/%i$"%(chi2,ndf)
		err = LF.rmsDeviation()/sqrt(len(gdat)) #1.0/sqrt(LF.sumWeights())
		print "dEp:",s,LF.coeffs[0],LF.rmsDeviation(),err
					
		gdEp.plot(graph.data.points(kepDelta[s],x=1,y=2,dy=3,
					title="%s: $\\mu=%.2f \pm %.3f$"%(sname,LF.coeffs[0],err)),
						[graph.style.symbol(symbol.circle,size=0.2,symbolattrs=[sideCols[s],]),
						graph.style.errorbar(errorbarattrs=[sideCols[s],])])
		gdEp.plot(graph.data.points(LF.fitcurve(0,gdat[-1][0],),x=1,y=2,title=None),[graph.style.line([sideCols[s],])])
		
		if s in bgRateDat:
			for afp in bgRateDat[s]:
				if len(bgRateDat[s][afp])<2:
					continue;
				LF.fit(bgRateDat[s][afp],cols=(0,1,2),errorbarWeights=True)
				chi2 = LF.ssResids()
				print chi2
				ndf = len(bgRateDat[s][afp])-len(LF.coeffs)
				gtitle = s+" Side, AFP=%s: $%.4f$Hz $\\chi^2/\\nu = %.1f/%i$"%(afp,LF.coeffs[0],chi2,ndf)
				if stats:
					gtitle += " $(p=%.2f)$"%stats.chisqprob(chi2,ndf)
				gBgRate.plot(graph.data.points(bgRateDat[s][afp],x=1,y=2,dy=3,title=gtitle),
							[graph.style.symbol(afpSymbs[afp],size=0.2,symbolattrs=[sideCols[s],]),
							graph.style.errorbar(errorbarattrs=[sideCols[s],])])
				gBgRate.plot(graph.data.points(LF.fitcurve(0,gdat[-1][0],),x=1,y=2,title=None),[graph.style.line([sideCols[s],])])

		if s in muRateDat:
			for afp in muRateDat[s]:
				if len(muRateDat[s][afp])<2:
					continue;
				gMuRate.plot(graph.data.points(muRateDat[s][afp],x=1,y=2,dy=3,title=None),
							 [graph.style.symbol(afpSymbs[afp],size=0.2,symbolattrs=[sideCols[s],]),
							  graph.style.errorbar(errorbarattrs=[sideCols[s],])])
								 #gMuRate.plot(graph.data.points(LF.fitcurve(0,gdat[-1][0],),x=1,y=2,title=None),[graph.style.line([sideCols[s],])])
					 
	gdEp.writetofile(basedir+"/Octet_dEP_%i.pdf"%depth)
	gBgRate.writetofile(basedir+"/Octet_BgRate_%i.pdf"%depth)
	gMuRate.writetofile(basedir+"/Octet_MuRate_%i.pdf"%depth)
	
			







					
def get_gain_tweak(conn,rn,s,t):
	conn.execute("SELECT e_orig,e_final FROM gain_tweak WHERE start_run <= %i AND %i <= end_run AND side = '%s'\
	 AND quadrant = %i ORDER BY end_run-start_run LIMIT 1"%(rn,rn,s,t))
	tweak = conn.fetchall()
	if not tweak:
		return (500.,500.)
	return tweak[0]
	
def delete_gain_tweak(conn,rn,s,t):
	conn.execute("DELETE FROM gain_tweak WHERE start_run <= %i AND %i <= end_run AND side = '%s' AND quadrant = %i"%(rn,rn,s,t))

def upload_gain_tweak(conn,rns,s,t,e0,e1):
	cmd = "INSERT INTO gain_tweak(start_run,end_run,side,quadrant,e_orig,e_final) VALUES (%i,%i,'%s',%i,%f,%f)"%(rns[0],rns[-1],s,t,e0,e1)
	print cmd
	conn.execute(cmd)


class MC_Comparator:
	"""Class for comparing data and MC results"""

	def __init__(self,basedir,simdir,depth=0):
		def sortByRuns(asyms):
			return dict([((af.getRuns()[0],af.getRuns()[-1]),af) for af in asyms])
		self.basedir = basedir
		self.simdir = simdir
		self.depth = abs(depth)
		self.datAsyms = sortByRuns(collectAsymmetries(basedir,depth))
		self.simAsyms = sortByRuns(collectAsymmetries(simdir,depth))
		self.rungrps = self.datAsyms.keys()
		self.rungrps.sort()
		for rns in self.rungrps:
			if rns not in self.simAsyms:
				print "**** Missing simulations for ",rns
				self.datAsyms.pop(rns)
		self.rungrps = self.datAsyms.keys()
		self.rungrps.sort()
		
	def endpoint_gain_tweak(self,conn=None):
		"""Set gain tweak factors to match spectrum endpoints between data and simulation"""
		# collect data
		gdat = {}
		for (n,rns) in enumerate(self.rungrps):
			print n,(rns[0],rns[-1]),"-"*80
			for s in ["East","West"]:
				for t in range(4):
					kdat = 0.5*(self.datAsyms[rns].getKurie(s[0],"Off","0",t).ep+self.datAsyms[rns].getKurie(s[0],"On","0",t).ep)
					ksim = 0.5*(self.simAsyms[rns].getKurie(s[0],"Off","0",t).ep+self.simAsyms[rns].getKurie(s[0],"On","0",t).ep)
					oldtweak = self.datAsyms[rns].getGMSTweak(s[0],t)
					gdat.setdefault((s,t),[]).append([n,kdat,ksim,oldtweak])
					print "\t%s %i:\t%f\t%f\t%f\tOld: %f"%(s,t,kdat,ksim,ksim/kdat,oldtweak)
					if conn:
						delete_gain_tweak(conn,rns[0],s,t)
						upload_gain_tweak(conn,rns,s,t,kdat/oldtweak,ksim)
		# plot
		tcols = [rgb.red,rgb.green,rgb.blue,rgb(1,1,0)]
		for s in ["East","West"]:
			gEp=graph.graphxy(width=30,height=10,
				x=graph.axis.lin(title=unitNames[self.depth],min=0,max=len(self.rungrps)-1),
				y=graph.axis.lin(title="Endpoint [keV]",min=700,max=780),
				key = graph.key.key(pos="tc",columns=2))
			setTexrunner(gEp)
			
			gTwk=graph.graphxy(width=15,height=15,
				x=graph.axis.lin(title="\\% Gain Tweak",min=-4.,max=4.),
				y=graph.axis.lin(title="Counts",min=0,max=15),
				key = graph.key.key(pos="tl"))
			setTexrunner(gTwk)
			
			gTwkOct=graph.graphxy(width=30,height=10,
								  x=graph.axis.lin(title=unitNames[self.depth],min=0,max=len(self.rungrps)-1),
								  y=graph.axis.lin(title="\\% Gain Tweak"),
								  key = graph.key.key(pos="bc",columns=4))
			setTexrunner(gTwkOct)

			
			for t in range(4):
				LF = LinearFitter(terms=[polyterm(0)])
				LFsim = LinearFitter(terms=[polyterm(0)])
				LF.fit(gdat[(s,t)],cols=(0,1))
				LFsim.fit(gdat[(s,t)],cols=(0,2))
				err = LF.rmsDeviation()
				gtitle = "%s %i: $\\mu=%.1f$, $\\sigma=%.1f$; $\\mu_{\\rm sim}=%.1f$"%(s,t,LF.coeffs[0],err,LFsim.coeffs[0])
				gEp.plot(graph.data.points(gdat[(s,t)],x=1,y=2,title=gtitle),[graph.style.symbol(symbol.circle,symbolattrs=[tcols[t],])])
				gEp.plot(graph.data.points(gdat[(s,t)],x=1,y=3,title=None),[graph.style.line([tcols[t]])])
				
				hTweak = histogram(40,-4.,4.)
				for g in gdat[(s,t)]:
					g.append(100*(g[2]/g[1]*g[3]-1.))
					hTweak.fill(g[-1],1.0,False)
				
				gtitle = "%s %i: $\\mu=%+.2f$\\%%, $\\sigma=%.2f$"%(s,t,hTweak.avg(),hTweak.rms())
				gTwk.plot(graph.data.points(hTweak.lineData(),x=1,y=2,title=gtitle),[graph.style.line([tcols[t]])])
				gtitle = "%s %i"%(s,t)
				gTwkOct.plot(graph.data.points(gdat[(s,t)],x=1,y=5,title=gtitle),[graph.style.symbol(symbol.circle,symbolattrs=[tcols[t],])])
							
			gEp.writetofile(self.basedir+"/TubeEp_%s_%i.pdf"%(s,self.depth))
			#gTwk.writetofile(self.basedir+"/TubeEp_%s_%i_Tweak.pdf"%(s,self.depth))
			gTwkOct.writetofile(self.basedir+"/TubeEp_%s_%i_Tweak.pdf"%(s,self.depth))
	
	def plot_endpoints(self):
		"""Compare data and MC spectrum endpoints"""
		gEp=graph.graphxy(width=30,height=10,
				x=graph.axis.lin(title=unitNames[self.depth],min=0,max=len(self.rungrps)-1),
				y=graph.axis.lin(title="Endpoint [keV]",min=760,max=820),
				key = graph.key.key(pos="bc",columns=2))
		setTexrunner(gEp)

			
		sideCols = {'E':rgb.red,'W':rgb.blue}
		LF = LinearFitter(terms=[polyterm(0)])
		LFsim = LinearFitter(terms=[polyterm(0)])
		
		for s in sideCols:
			for afp in afpSymbs:
				gdat = [(n,
							self.datAsyms[rns].getKurie(s,afp,"0",4).ep,
							self.simAsyms[rns].getKurie(s,afp,"0",4).ep)  for (n,rns) in enumerate(self.rungrps)]
				
				LF.fit(gdat,cols=(0,1))
				LFsim.fit(gdat,cols=(0,2))
				err = LF.rmsDeviation()
				gtitle = "%s$_{\\rm %s}$: $\\mu=%.1f$, $\\sigma=%.1f$; $\\mu_{\\rm sim}=%.1f$"%(s,afp,LF.coeffs[0],err,LFsim.coeffs[0])
				gEp.plot(graph.data.points(gdat,x=1,y=2,dy=3,title=gtitle),
							[graph.style.symbol(afpSymbs[afp],size=0.2,symbolattrs=[sideCols[s],])])
				gEp.plot(graph.data.points(gdat,x=1,y=3,title=None),[graph.style.line([sideCols[s],afpLines[afp]])])
		
		gEp.writetofile(self.basedir+"/OctetEP_%i.pdf"%self.depth)
					
							
	def backscatter_fractions(self):
		"""Plot data and simulation backscatter fractions"""
		nmax = len(self.rungrps)-1
		gScatter=graph.graphxy(width=30,height=15,
			x=graph.axis.lin(title=unitNames[self.depth],min=-nmax/4,max=nmax),
			y=graph.axis.lin(title="Backscatter Fraction (\\% of Type 0)",min=0,max=4.0),
			key = graph.key.key(pos="ml"))
		setTexrunner(gScatter)
		
		scols = {'E':rgb.red,'W':rgb.blue}
		afpFills = {"Off":[],"On":[deco.filled]}
		afpThick = {"Off":style.linewidth.thick,"On":style.linewidth.THick}
		typeSymbs = {1:symbol.circle,2:symbol.triangle,3:symbol.square}
		typeLines = {1:style.linestyle.solid,2:style.linestyle.dashed,3:style.linestyle.dotted}
		
		for s in ["E","W"]:
			for afp in ["Off","On"]:
				gdat = []
				for (n,rns) in enumerate(self.rungrps):
					dat = [ self.datAsyms[rns].getRate(s,afp,'1',"hEnergy_Type_%i_"%t).counts for t in range(4) ]
					sim = [ self.simAsyms[rns].getRate(s,afp,'1',"hEnergy_Type_%i_"%t).counts for t in range(4) ]
					gdat.append([n,rns[0]] + [100.*dat[t]/dat[0] for t in range(4)[1:]] + [100.*sim[t]/sim[0] for t in range(4)[1:]])
					print n,rns
					print gdat[-1]
				
				for t in range(4)[1:]:
					gScatter.plot(graph.data.points(gdat,x=1,y=2+t,title="%s$_{\\rm %s}$ Type "%(s,afp)+"I"*t),
								[graph.style.symbol(typeSymbs[t],symbolattrs=[scols[s],]+afpFills[afp]),])
					gScatter.plot(graph.data.points(gdat,x=1,y=5+t,title=None),
								[graph.style.line([typeLines[t],scols[s],afpThick[afp]]),])
					
		gScatter.writetofile(self.basedir+"/BackscatterFrac_%i.pdf"%self.depth)
		
		
		
def backscatterFracTableLine(fname):
	AF = AsymmetryFile(fname)
	counts = [0,0,0,0]
	for s in ["E","W"]:
		for afp in ["Off","On"]:
			for t in range(4):
				counts[t] += AF.getRate(s,afp,'1',"hEnergy_Type_%i_"%t).counts
	return tuple([100.*counts[1]/counts[0],100.*counts[2]/counts[0],100.*counts[3]/counts[0]])

def backscatterFracTable(simV = "OctetAsym_Offic_SimMagF"):
	print "\\begin{table} \\centering \\begin{tabular}{|c||c|c|c|} \\hline"
	print "& Type I & Type II & Type III \\\\ \\hline \\hline"
	
	
	sOct = None
	#sOct = "15764-15789_0"
	
	datf = None
	simf = None
	if sOct:
		datf = backscatterFracTableLine(os.environ["UCNA_ANA_PLOTS"]+"/OctetAsym_Offic/%s/%s.txt"%(sOct,sOct))
		simf = backscatterFracTableLine(os.environ["UCNA_ANA_PLOTS"]+"/%s/%s/%s.txt"%(simV,sOct,sOct))
	else:
		datf = backscatterFracTableLine(os.environ["UCNA_ANA_PLOTS"]+"/OctetAsym_Offic/OctetAsym_Offic.txt")
		simf = backscatterFracTableLine(os.environ["UCNA_ANA_PLOTS"]+"/%s/%s.txt"%(simV,simV))
		
	print "Data & %.2f\\%% & %.2f\\%% & %.2f\\%% \\\\ \\hline"%datf
	print "MC & %.2f\\%% & %.2f\\%% & %.2f\\%% \\\\ \\hline"%simf
	print "Error & %.0f\\%% & %.0f\\%% & %.0f\\%% \\\\ \\hline"%tuple([ 100*(simf[t]-datf[t])/datf[t] for  t in range(3)])
	print "\\end{tabular}\\caption{Backscatter events as fraction of Type 0}\\end{table}"



if __name__=="__main__":
	
	#backscatterFracTable("OctetAsym_Offic_SimMagF")
	#exit(0)
	
	if 0:
		MCC = MC_Comparator(os.environ["UCNA_ANA_PLOTS"]+"/OctetAsym_Offic/",os.environ["UCNA_ANA_PLOTS"]+"/OctetAsym_Offic_SimMagF/")
		#MCC.backscatter_fractions()
		MCC.plot_endpoints()
		#exit(0)
		
	
		conn=open_connection()
		#conn=None
		MCC.endpoint_gain_tweak(conn)
		exit(0)
	
	for i in range(3):				
		plot_octet_asymmetries(os.environ["UCNA_ANA_PLOTS"]+"/OctetAsym_Offic/",i)
		#plot_octet_asymmetries(os.environ["UCNA_ANA_PLOTS"]+"/OctetAsym_Offic_Sim_MagF_2",2-i)
			