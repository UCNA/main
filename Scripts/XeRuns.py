#!/usr/bin/python

from Asymmetries import *

class XeFit(KVMap):
	def __init__(self,m=KVMap()):
		self.dat = m.dat
		self.loadFloats(["xe_hi","d_xe_hi","xe_lo","d_xe_lo","tube"])
		self.tube = int(self.tube)
		self.loadStrings(["side"])

class XeDecomp(KVMap):
	def __init__(self,m=KVMap()):
		self.dat = m.dat
		self.nIsots = int(self.getFirstF("nTerms"))
		self.names = self.getFirst("isots").split(",")
		self.counts = [float(x) for x in self.getFirst("counts").split(",")]
		self.terms = [float(x) for x in self.getFirst("terms").split(",")]
		self.err = [float(x) for x in self.getFirst("errs").split(",")]
							  
class XeFile(QFile):
	def __init__(self,fname):
		QFile.__init__(self,fname)
		self.fits = {}
		for m in [XeFit(m) for m in self.dat.get("tuben",())]:
			self.fits[(m.side,m.tube)] = m
		self.runcals = dict([(r.run,r) for r in [runCal(m) for m in self.dat.get("runcal",[])]])
		scomps = self.dat.get("spectrumComp",[])
		if scomps:
			self.xecomp = XeDecomp(scomps[0])
		self.runtimes = {}
		for rt in self.dat.get("runTimes",[]):
			for r in rt.dat:
				self.runtimes[int(r)] = rt.getFirstF(r)
							
def XeGainTweak(rn,conn):
	datpath = os.environ["UCNA_ANA_PLOTS"]+"/PositionMaps/SingleRuns/"
	simpath = os.environ["UCNA_ANA_PLOTS"]+"/PositionMaps/SingleRunsSim/"
	rname = "Xenon_%i_12_52"%rn
	
	xdat = XeFile(datpath+"/"+rname+"/"+rname+".txt")
	xsim = XeFile(simpath+"/"+rname+"/"+rname+".txt")
	
	print "----------",rn,"----------"
	for s in ["East","West"]:
		for t in range(4):
			ldat = xdat.fits[(s[0],t)].xe_hi
			lsim = xsim.fits[(s[0],t)].xe_hi
			oldtweak = xdat.runcals[rn].getGMSTweak(s[0],t)
			print "\t",s,t,ldat,"->",lsim,"\t Err =",100.0*(lsim-ldat)/lsim,"\tOld =",100.*(oldtweak-1)
			if(conn):
				delete_gain_tweak(conn,rn,s,t)
				upload_gain_tweak(conn,[rn],s,t,ldat/oldtweak,lsim)
							
							
def XeTimeEvolution(rmin,rmax):
	simpath = os.environ["UCNA_ANA_PLOTS"]+"/PositionMaps/SingleRunsSim/"
	conn = open_connection()
	isotdat = {}						
	
	# collect data
	tmin = 1e100
	for rn in range(rmin,rmax+1):
		try:
			rname = "Xenon_%i_12_52"%rn
			xsim = XeFile(simpath+"/"+rname+"/"+rname+".txt")
			trange = (getRunStartTime(conn,rn),getRunEndTime(conn,rn))
			rtime = xsim.runtimes[rn]
			print rn,trange
			if trange[0]<tmin:
				tmin = trange[0]
			for i in range(xsim.xecomp.nIsots):
				dpt = [trange[0]+0.5*rtime,xsim.xecomp.counts[i]/rtime]
				dpt.append(dpt[-1]*xsim.xecomp.err[i])
				isotdat.setdefault(xsim.xecomp.names[i],[]).append(dpt)
		except:
			print "*** Fail on run",rn,"***"
	
	# plot
	gIA=graph.graphxy(width=15,height=15,
					  x=graph.axis.lin(title="Time [h]",min=0),
					  y=graph.axis.log(title="Decay rate [Hz]",min=20,max=300),
					  key = graph.key.key(pos="tr"))
	setTexrunner(gIA)
	icols = rainbowDict(isotdat)
	ks = isotdat.keys()
	ks.sort()
	for k in ks:
		for d in isotdat[k]:
			d[0] = (d[0]-tmin)/3600.
		LF = LogYer(terms=[polyterm(0),polyterm(1)])
		LF.fit(isotdat[k],cols=(0,1))
		thalf = 0
		if LF.coeffs[1]:
			thalf = -log(2)/LF.coeffs[1]
		gtitle = k.replace("_"," ")
		gtitle = "$^{"+gtitle[2:5]+"}$Xe"+gtitle[5:-1].replace("-","/")+"$^{"+gtitle[-1]+"}$: $T_{1/2}$ = "
		if abs(thalf) < 24:
			gtitle += "%.1f h"%thalf
		else:
			gtitle += "%.1f d"%(thalf/24)
		gIA.plot(graph.data.points(isotdat[k],x=1,y=2,dy=3,title=gtitle),
					[graph.style.symbol(symbol.circle,size=0.2,symbolattrs=[icols[k],]),
					graph.style.errorbar(errorbarattrs=[icols[k],])])
		gIA.plot(graph.data.points(LF.fitcurve(0,30),x=1,y=2,title=None),[graph.style.line([icols[k],])])
	gIA.writetofile(os.environ["UCNA_ANA_PLOTS"]+"/test/XeDecomp/DecompHistory_%i-%i.pdf"%(rmin,rmax))
							
if __name__ == "__main__":
	
	XeTimeEvolution(14283,14333)
	exit(0)
	
	conn = open_connection()
	#conn = None
	for rn in range(14282,14347+1):
		XeGainTweak(rn,conn)

