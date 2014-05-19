from QFile import *
from PyxUtils import *

# runs in pulse pair, half octet, and octet groups			  		  
#ppSegments = [["A1","A2","A3","A4","A5","A6"],["A7","A8","A9","A10","A11","A12"],["B1","B2","B3","B4","B5","B6"],["B7","B8","B9","B10","B11","B12"]]
ppSegments = [["A1","A2","A4","A5"],["A7","A9","A10","A12"],["B1","B2","B4","B5"],["B7","B9","B10","B12"]]
hoSegments = [ppSegments[0]+ppSegments[1],ppSegments[2]+ppSegments[3]]
octSegments = [hoSegments[0]+hoSegments[1],]
divSegments = [octSegments,hoSegments,ppSegments]
unitNames = {0:"Octet",1:"Half Octet",2:"Pulse Pair",3:"Run"}

class runCal(KVMap):
	def __init__(self,m):
		self.dat = m.dat
		self.loadInts(["run","tstart","tend"])
		self.eOrig = self.get_tubesparam("eOrig")
		self.eFinal = self.get_tubesparam("eFinal")
		self.gms = self.get_tubesparam("gms")

	def get_tubesparam(self,pname):
		return dict( [ ((s,t),self.getFirstF("%s%i_%s"%(s,t,pname))) for s in ['E','W'] for t in range(5) ] )

	def getGMSTweak(self,side,tube):
		return self.eFinal[(side,tube)]/self.eOrig[(side,tube)]

	def midTime(self):
		return 0.5*(self.tend+self.tstart)

class eventRate(KVMap):
	def __init__(self,m):
		self.dat = m.dat
		self.loadFloats(["counts","d_counts","rate","d_rate"])
		self.loadStrings(["side","afp","fg","name"])
	def dRate(self):
		if not self.counts:
			return 0
		return self.rate/sqrt(self.counts)

# base class for files produced by RunAccumulators
class RunAccumulatorFile(QFile):
	def __init__(self,fname):
		QFile.__init__(self,fname)

		# included runs in file
		self.runcals = dict([(r.run,r) for r in [runCal(m) for m in self.dat.get("runcal",[])]])
		self.runCounts = dict([(int(i[0]),float(i[1][0])) for i in self.getFirst("runCounts").dat.items()])
		self.runTimes = dict([(int(i[0]),float(i[1][0])) for i in self.getFirst("runTimes").dat.items()])
		
		# histogram rate summaries
		self.rates = [eventRate(m) for m in self.dat.get("rate",[])]
		
		# list of included runs by octet
		self.octruns = {}
		for m in self.dat.get("Octet",[]):
			for k in m.dat:
				if k in octSegments[0]:
					for rns in m.dat[k]:
						for rn in rns.split(","):
							self.addOctRun(k,rn)

	def addOctRun(self,rtype,rn):
		"""add run to runs-by-octet listing"""
		self.octruns.setdefault(rtype,[]).append(int(rn))

	def getRuns(self):
		"""get list of runs included in file"""
		rns = self.runCounts.keys()
		rns.sort()
		return rns

	def getRate(self,side,afp,fg,name):
		"""get rate info for named histogram"""
		hname = name
		if hname[-1] == "_":
			hname = "%s%s_%s"%(name,side,afp)
		for r in self.rates:
			if r.side==side and r.afp==afp and r.fg==fg and r.name==hname:
				return r
		print "Missing rate for",hname,fg
		return None

def collectOctetFiles(basedir,depth,baseClass=RunAccumulatorFile):
	"""Collect octet-based files"""
	print "Collecting octet data from",basedir,"at depth",depth
	datfls = []
	fls = os.listdir(basedir)
	fls.sort()
	for f in [ (basedir+'/'+f+'/',basedir+'/'+f+'/'+f+'.txt') for f in fls if os.path.isdir(basedir+'/'+f) and '-' in f]:
		if not depth and os.path.exists(f[1]):
			datfls.append(baseClass(f[1]))
		else:
			datfls += collectOctetFiles(f[0],depth-1)
	return datfls

def make_runaxis(rmin,rmax):
	
	tckdist = [5,1]
	if rmax-rmin > 100:
		tckdist = [10,1]
	if rmax-rmin > 500:
		tckdist = [100,10]	
	if rmax-rmin > 1000:
		tckdist = [100,20]
				
	return graph.axis.lin(title="Run Number",min=rmin,max=rmax,
							parter=graph.axis.parter.linear(tickdists=tckdist),
							texter = graph.axis.texter.rational(),
							painter=graph.axis.painter.regular(labeldist=0.1,labeldirection=graph.axis.painter.rotatetext(135)))
