#!/sw/bin/python2.6

from LinFitter import *
from PyxUtils import *
from QFile import *
from scipy import stats
	
# runs in pulse pair, half octet, and octet groups			  		  
#ppSegments = [["A1","A2","A3","A4","A5","A6"],["A7","A8","A9","A10","A11","A12"],["B1","B2","B3","B4","B5","B6"],["B7","B8","B9","B10","B11","B12"]]
ppSegments = [["A1","A2","A4","A5"],["A7","A9","A10","A12"],["B1","B2","B4","B5"],["B7","B9","B10","B12"]]
hoSegments = [ppSegments[0]+ppSegments[1],ppSegments[2]+ppSegments[3]]
octSegments = [hoSegments[0]+hoSegments[1],]
divSegments = [octSegments,hoSegments,ppSegments]

def addQuad(a,b):
	return sqrt(a**2+b**2)
	
class asymmetryFit:
	def __init__(self,m=KVMap()):
		self.fitMin = m.getFirstF("fitMin",0)
		self.fitMax = m.getFirstF("fitMax",0)
		self.A0 = m.getFirstF("A0_fit",0)
		self.dA0 = m.getFirstF("dA0",0)
	def __repr__(self):
		return "<Asym: %.2f +- %.2f, (%g,%g)>"%(self.A0,self.dA0,self.fitMin,self.fitMax)

class eventRate(KVMap):
	def __init__(self,m):
		self.dat = m.dat
		self.loadFloats(["counts","rate"])
		self.loadStrings(["side","afp","fg","name"])
	def dRate(self):
		if not self.counts:
			return 0
		return self.rate/sqrt(self.counts)
		
class kurieFit(KVMap):
	def __init__(self,m):
		self.dat = m.dat
		self.loadFloats(["endpoint","dendpoint"])
		self.loadStrings(["side","afp","type"])
		self.ep = self.endpoint
		self.dep = self.dendpoint

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
	
	def addRun(self,rtype,rn):
		self.octruns.setdefault(rtype,[]).append(int(rn))
	
	def getAsym(self,fitMin=225,fitMax=675):
		for a in self.asyms:
			if a.fitMin==fitMin and a.fitMax==fitMax:
				return a
		return None
	
		
		
class AsymmetryFile(QFile,asymData):

	def __init__(self,fname):
		QFile.__init__(self,fname)
		asymData.__init__(self)
		self.asyms = [asymmetryFit(m) for m in self.dat.get("asymmetry",[])]
		self.kuries = [kurieFit(m) for m in self.dat.get("kurieFit",[])]
		self.rates = [eventRate(m) for m in self.dat.get("rate",[])]
		for m in self.dat.get("Octet"):
			for k in m.dat:
				if k in octSegments[0]:
					for rns in m.dat[k]:
						for rn in rns.split(","):
							self.addRun(k,rn)
		
	def kepDelta(self,side,type="0"):
		"""Calculate difference between On/Off endpoint"""
		kOn = self.getKurie(side,'1',type)
		kOff = self.getKurie(side,'0',type)
		return [kOn.ep-kOff.ep,addQuad(kOn.dep,kOff.dep)]
	
	def getKurie(self,side,afp,type):
		for k in self.kuries:
			if k.side==side and k.afp==afp and k.type==type:
				return k
		print "*** Can't find kurie in",side,afp,type,self.fname
		return None
	
	def getRate(self,side,afp,fg,name):
		for r in self.rates:
			if r.side==side and r.afp==afp and r.fg==fg and r.name==name:
				return r
		return None
		
def collectAsymmetries(basedir,depth):
	asyms = []
	for f in [ (basedir+'/'+f+'/',basedir+'/'+f+'/'+f+'.txt') for f in os.listdir(basedir) if os.path.isdir(basedir+'/'+f) and '-' in f and f[:3]!="139"]:
		if not depth and os.path.exists(f[1]):
			asyms.append(AsymmetryFile(f[1]))
		else:
			asyms += collectAsymmetries(f[0],depth-1)
	return asyms
	
def plot_octet_asymmetries(basedir="../PostPlots/OctetAsym_div0",depth=0):
	
	gdat = []
	bgRateDat = {'E':{'0':[],'1':[]},'W':{'0':[],'1':[]}}
	kdat = {'E':{'0':[],'1':[]},'W':{'0':[],'1':[]}}
	kepDelta = {'E':[],'W':[],'C':[]}
	n=0
	afpNames = {'0':"Off",'1':"On"}
	print "--------------------- Division",depth,"-------------------------"
	for af in collectAsymmetries(basedir,depth):
		a = af.getAsym()
		if a:
			if abs(a.A0+0.115) > 0.05 or abs(a.A0+0.115)/a.dA0 > 3:
				print "*** Funny asym",a.A0,"+-",a.dA0,"in",af.fname
			gdat.append([n,a.A0,a.dA0,1.0/a.dA0**2])
		else:
			print "**** Missing asym in",af.fname
		for s in kdat:
			kepDelta[s].append([n,]+af.kepDelta(s))
			if not (-15 < kepDelta[s][-1][1] < 15 and 0 < kepDelta[s][-1][2]):
				print "*** Funny delta",kepDelta[s][-1][1],kepDelta[s][-1][2],"in",af.fname
			for afp in kdat[s]:
				k = af.getKurie(s,afp,"0")
				if k:
					kdat[s][afp].append([n,k.ep,k.dep])
					if not 780 < k.ep < 840:
						print "*** Funny Endpoint",k.ep,k.side,k.afp,"in",af.fname
				rt = af.getRate(s,afp,'0',"hEnergy_Type_0_%s_%s"%(s,afpNames[afp]))
				if rt:
					if not 0.10 < rt.rate < 0.25:
						print "*** Funny rate",rt.rate,rt.side,rt.afp,"in",af.fname
					if not 0.001 < rt.dRate():
						print "**** Funny dRate",rt.dRate(),rt.side,rt.afp,"in",af.fname
						continue
					bgRateDat[s][afp].append([n,rt.rate,rt.dRate()])
				else:
					print "*** Can't find rate for",n,s,afp
		kepDelta['C'].append([n,0.5*(kepDelta['E'][-1][1]-kepDelta['W'][-1][1]),0.5*addQuad(kepDelta['E'][-1][2],kepDelta['W'][-1][2])])
		n+=1
	print

	unitName = "Octet"
	if depth == 1:
		unitName = "Half Octet"
	if depth == 2:
		unitName = "Pulse Pair"
		
	gAsyms=graph.graphxy(width=25,height=8,
				x=graph.axis.lin(title=unitName,min=0,max=gdat[-1][0]),
				y=graph.axis.lin(title="Asymmetry",min=-0.2,max=0),
				key = graph.key.key(pos="bl"))
	gAsyms.texrunner.set(lfs='foils17pt')
	
	gBgRate=graph.graphxy(width=25,height=8,
				x=graph.axis.lin(title=unitName,min=0,max=gdat[-1][0]),
				y=graph.axis.lin(title="Background Rate [Hz]",min=0,max=0.25),
				key = graph.key.key(pos="bl"))
	gBgRate.texrunner.set(lfs='foils17pt')
	
	gEp=graph.graphxy(width=25,height=8,
				x=graph.axis.lin(title=unitName,min=0,max=gdat[-1][0]),
				y=graph.axis.lin(title="Endpoint [keV]",min=750,max=850),
				key = graph.key.key(pos="bl"))
	gEp.texrunner.set(lfs='foils17pt')
	
	gdEp=graph.graphxy(width=25,height=8,
				x=graph.axis.lin(title=unitName,min=0,max=gdat[-1][0]),
				y=graph.axis.lin(title="Endpoint On/Off Difference [keV]",min=-30,max=20),
				key = graph.key.key(pos="bl"))
	gdEp.texrunner.set(lfs='foils17pt')

	
	gAsyms.plot(graph.data.points(gdat,x=1,y=2,dy=3,title=None),
				[graph.style.symbol(symbol.circle,size=0.2,symbolattrs=[rgb.red,]),
				graph.style.errorbar(errorbarattrs=[rgb.red,])])
	
	LF = LinearFitter(terms=[polyterm(0)])
	LF.fit(gdat,cols=(0,1,2),errorbarWeights=True)
	chi2 = LF.ssResids()
	ndf = len(gdat)-len(LF.coeffs)
	gAsyms.plot(graph.data.points(LF.fitcurve(0,gdat[-1][0],),x=1,y=2,
				title="$A=%.5f\\pm%.5f$, $\\chi^2/\\nu = %.1f/%i$ $(p=%.2f)$"%(LF.coeffs[0],1.0/sqrt(LF.sumWeights()),chi2,ndf,stats.chisqprob(chi2,ndf))),
				[graph.style.line()])
	
	gAsyms.writetofile(basedir+"/OctetAsym_%i.pdf"%depth)
	
	sideCols = {'E':rgb.red,'W':rgb.blue,'C':rgb(0.,0.7,0.)}
	afpSymbs = {'0':symbol.circle,'1':symbol.triangle}
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
		
		if s in kdat:
			for afp in kdat[s]:
				LF.fit([kd for kd in kdat[s][afp] if 0<kd[2]],cols=(0,1,2),errorbarWeights=True)

				gEp.plot(graph.data.points(kdat[s][afp],x=1,y=2,dy=3,title=s+" Side, AFP=%s: $\\mu=%.1f$, $\\sigma=%.1f$"%(afp,LF.coeffs[0],LF.rmsDeviation())),
							[graph.style.symbol(afpSymbs[afp],size=0.2,symbolattrs=[sideCols[s],]),
							graph.style.errorbar(errorbarattrs=[sideCols[s],])])
				gEp.plot(graph.data.points(LF.fitcurve(0,gdat[-1][0],),x=1,y=2,title=None),[graph.style.line([sideCols[s],])])
		
		if s in bgRateDat:
			for afp in bgRateDat[s]:
				if len(bgRateDat[s][afp])<2:
					continue;
				LF.fit(bgRateDat[s][afp],cols=(0,1,2),errorbarWeights=True)
				chi2 = LF.ssResids()
				print chi2
				ndf = len(bgRateDat[s][afp])-len(LF.coeffs)
				gBgRate.plot(graph.data.points(bgRateDat[s][afp],x=1,y=2,dy=3,
								title=s+" Side, AFP=%s: $%.4f$Hz $\\chi^2/\\nu = %.1f/%i$ $(p=%.2f)$"%(afp,LF.coeffs[0],chi2,ndf,stats.chisqprob(chi2,ndf))),
							[graph.style.symbol(afpSymbs[afp],size=0.2,symbolattrs=[sideCols[s],]),
							graph.style.errorbar(errorbarattrs=[sideCols[s],])])
				gBgRate.plot(graph.data.points(LF.fitcurve(0,gdat[-1][0],),x=1,y=2,title=None),[graph.style.line([sideCols[s],])])

				
	gEp.writetofile(basedir+"/OctetEP_%i.pdf"%depth)
	gdEp.writetofile(basedir+"/Octet_dEP_%i.pdf"%depth)
	gBgRate.writetofile(basedir+"/Octet_BgRate_%i.pdf"%depth)
	
if __name__=="__main__":
	for i in range(3):
		plot_octet_asymmetries("../../PostPlots/OctetAsym_Offic/",2-i)
		#plot_octet_asymmetries("../PostPlots/OctetAsym_Offic_Simulated",2-i)
	