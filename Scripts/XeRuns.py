#!/usr/bin/python

from Asymmetries import *
import sys
import traceback

class SectorCutter(KVMap):
	def __init__(self,m=KVMap()):
		self.dat = m.dat
		self.loadFloats(["nRings","nSectors","radius"])
		self.nRings = int(self.nRings)
		self.nSectors = int (self.nSectors)

class XeFit(KVMap):
	def __init__(self,m=KVMap()):
		self.dat = m.dat
		self.loadFloats(["xe_hi","d_xe_hi","xe_lo","d_xe_lo","xe_lo_w","d_xe_lo_w","tube","m","eta"])
		self.tube = int(self.tube)
		self.m = int(self.m)
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
		
		self.sects = SectorCutter(self.getFirst("SectorCutter_Xe"))
		
		self.tuben = {}
		for m in [XeFit(m) for m in self.dat.get("sectDat",[]) if m.dat.get("m")==["%i"%self.sects.nSectors]]:
			self.tuben[(m.side,m.tube)] = m
		
		self.sectdat = {}
		for m in [XeFit(m) for m in self.dat.get("sectDat",[])]:
			self.sectdat[(m.side,m.tube,m.m)] = m
				
		self.runcals = dict([(r.run,r) for r in [runCal(m) for m in self.dat.get("runcal",[])]])
		scomps = self.dat.get("spectrumComp",[])
		if scomps:
			self.xecomp = XeDecomp(scomps[0])
		self.runtimes = {}
		for rt in self.dat.get("runTimes",[]):
			for r in rt.dat:
				self.runtimes[int(r)] = rt.getFirstF(r)
							
def XeGainTweak(rn,conn,nrings):
	datpath = os.environ["UCNA_ANA_PLOTS"]+"/PositionMaps/SingleRuns/"
	simpath = os.environ["UCNA_ANA_PLOTS"]+"/PositionMaps/SingleRunsSim/"
	rname = "Xenon_%i_%i_50"%(rn,nrings)
	
	xdat = XeFile(datpath+"/"+rname+"/"+rname+".txt")
	xsim = XeFile(simpath+"/"+rname+"/"+rname+".txt")
	
	print "----------",rn,"----------"
	for s in ["East","West"]:
		for t in range(4):
			ldat = xdat.tuben[(s[0],t)].xe_hi
			lsim = xsim.tuben[(s[0],t)].xe_hi
			oldtweak = xdat.runcals[rn].getGMSTweak(s[0],t)
			if not lsim or not ldat or not abs((lsim-ldat)/lsim) < 0.10:
				print "\t***** BAD FIT",s,t,ldat,"->",lsim
				continue
			print "\t",s,t,"%.2f -> %.2f\t\tErr = %+.2f%%\tOld = %+.2f%%"%(ldat,lsim,100.0*(lsim-ldat)/lsim,100.*(oldtweak-1))
			if(conn):
				delete_gain_tweak(conn,rn,s,t)
				upload_gain_tweak(conn,[rn],s,t,ldat/oldtweak,lsim)
							
							
def XeTimeEvolution(rmin,rmax,nrings):
	simpath = os.environ["UCNA_ANA_PLOTS"]+"/PositionMaps/SingleRunsSim/"
	conn = open_connection()
	isotdat = {}						
	
	# collect data
	tmin = 1e100
	for rn in range(rmin,rmax+1):
		try:
			rname = "Xenon_%i_%i_50"%(rn,nrings)
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
			#traceback.print_exc(file=sys.stdout)
			#exit(1)
	
	# plot
	gIA=graph.graphxy(width=15,height=15,
					  x=graph.axis.lin(title="Time [h]",min=0),
					  y=graph.axis.log(title="Decay rate [Hz]",min=5),
					  key = graph.key.key(pos="tr"))
	setTexrunner(gIA)
	icols = rainbowDict(isotdat)
	ks = isotdat.keys()
	ks.sort()
	for k in ks:
		for d in isotdat[k]:
			d[0] = (d[0]-tmin)/3600.
		LF = LogYer(terms=[polyterm(0),polyterm(1)])
		LF.fit([d for d in isotdat[k] if d[1]>20],cols=(0,1))
		thalf = 0
		if LF.coeffs[1]:
			thalf = -log(2)/LF.coeffs[1]
		gtitle = k.replace("_"," ")
		gtitle = "$^{"+gtitle[2:5]+"}$Xe"+gtitle[5:-1].replace("-","/")+"$^{"+gtitle[-1]+"}$: $T_{1/2}$ = "
		if abs(thalf) < 1.0:
			gtitle += "%.1f m"%(60*thalf)
		elif abs(thalf) < 24:
			gtitle += "%.1f h"%thalf
		else:
			gtitle += "%.1f d"%(thalf/24)
		gIA.plot(graph.data.points(isotdat[k],x=1,y=2,dy=3,title=gtitle),
					[graph.style.symbol(symbol.circle,size=0.2,symbolattrs=[icols[k],]),
					graph.style.errorbar(errorbarattrs=[icols[k],])])
		gIA.plot(graph.data.points(LF.fitcurve(0,30),x=1,y=2,title=None),[graph.style.line([icols[k],])])
	gIA.writetofile(os.environ["UCNA_ANA_PLOTS"]+"/test/XeDecomp/DecompHistory_%i-%i.pdf"%(rmin,rmax))


def ep_v_eta(fname):
	
	basepath = os.environ["UCNA_ANA_PLOTS"]+"/PositionMaps/"+fname+"/"
	xdat = XeFile(basepath+"/"+fname+".txt")
				  
	gEp=graph.graphxy(width=20,height=10,
					  x=graph.axis.lin(title="eta",min=0,max=2.5),
					  y=graph.axis.lin(title="Endpoint",min=800,max=1000),
					  key = graph.key.key(pos="tr"))
	setTexrunner(gEp)
	
	gHvL=graph.graphxy(width=20,height=20,
					  x=graph.axis.lin(title="Low Peak",min=0),
					  y=graph.axis.lin(title="Endpoint",min=0),
					  key = graph.key.key(pos="tl"))
	setTexrunner(gHvL)

	
	tcols = rainbow(4)
	ssymbs = {"E":symbol.circle,"W":symbol.triangle}
	for s in ["E","W"]:
		for t in range(4):
			print s,t
			gdat = [xdat.sectdat[(s,t,m)] for m in range(xdat.sects.nSectors)]
			#        0      1        2          3              4         
			gdat = [(g.eta, g.xe_hi, g.d_xe_hi, g.xe_lo*g.eta, g.xe_hi*g.eta) for g in gdat]
			gEp.plot(graph.data.points(gdat,x=1,y=2,dy=3,title=None),
					 [graph.style.symbol(ssymbs[s],size=0.2,symbolattrs=[tcols[t],]),
					  graph.style.errorbar(errorbarattrs=[tcols[t],])])
			gHvL.plot(graph.data.points(gdat,x=4,y=5,title=None),
					[graph.style.symbol(ssymbs[s],size=0.2,symbolattrs=[tcols[t],])])
			LF = LinearFitter(terms = [polyterm(0),polyterm(1)])
			LF.fit(gdat,cols=(3,4))
			gHvL.plot(graph.data.points(LF.fitcurve(0,250),x=1,y=2,title=s+" %i: $%s$"%(t,LF.toLatex())),
					 [graph.style.line(lineattrs=[tcols[t],])])

	gEp.writetofile(basepath+"/Ep_v_Eta.pdf")
	gHvL.writetofile(basepath+"/Hi_v_Lo_Raw.pdf")


def data_v_sim(rmin,rmax,nrings):
	
	basepath = os.environ["UCNA_ANA_PLOTS"]+"/PositionMaps/"
	#rname = "%i-%i_%i"%(rmin,rmax,nrings)
	rname = "%i-%i"%(rmin,rmax)
	xdat = XeFile(basepath+"/Xenon_"+rname+"/Xenon_"+rname+".txt")
	xsim = XeFile(basepath+"/SimXe_"+rname+"/SimXe_"+rname+".txt")

	gHvL=graph.graphxy(width=20,height=20,
					   x=graph.axis.lin(title="$\eta$ from low peak",min=0,max=2.5),
					   y=graph.axis.lin(title="$\eta$ from endpoint",min=0,max=2.5),
					   key = graph.key.key(pos="tl"))
	setTexrunner(gHvL)

	tcols = rainbow(4)
	ssymbs = {"E":symbol.circle,"W":symbol.triangle}
	for s in ["E","W"]:
		for t in range(4):
			print s,t
			gdat = [ (xdat.sectdat[(s,t,m)],xsim.sectdat[(s,t,m)]) for m in range(xdat.sects.nSectors)]
			gdat = [ (g[0].eta*g[0].xe_lo/g[1].xe_lo,g[0].eta*g[0].xe_hi/g[1].xe_hi) for g in gdat ]
			gHvL.plot(graph.data.points(gdat,x=1,y=2,title=None),
					  [graph.style.symbol(ssymbs[s],size=0.2,symbolattrs=[tcols[t],])])
			LF = LinearFitter(terms = [polyterm(0),polyterm(1)])
			LF.fit(gdat,cols=(0,1))
			gHvL.plot(graph.data.points(LF.fitcurve(0,2.5),x=1,y=2,title=s+" %i: $%s$"%(t,LF.toLatex())),
					  [graph.style.line(lineattrs=[tcols[t],])])
	
	gHvL.writetofile(basepath+"/Xenon_"+rname+"/Hi_v_Lo.pdf")

if __name__ == "__main__":
	
	#XeTimeEvolution(14283,14333,15)
	#XeTimeEvolution(15992,16077,15)
	#exit(0)
	
	#data_v_sim(14282,14347,12)
	#exit(0)
	
	#ep_v_eta("Xenon_14282-14347")
	#ep_v_eta("SimXe_14282-14347")
	#exit(0)
	
	conn = open_connection()
	#conn = None
	for rn in range(14282,14347+1):
	#for rn in range(15991,16077+1):
		XeGainTweak(rn,conn,15)

