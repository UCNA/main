#!/usr/bin/python

import sys
sys.path.append("..")

from ucnacore.AnaDB import *
from Asymmetries import *
from pyx import *
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
		self.runtimes = {}
		for rt in self.dat.get("runTimes",[]):
			for r in rt.dat:
				self.runtimes[int(r)] = rt.getFirstF(r)
							
def XeGainTweak(rn,conn,aconn=open_anadb_connection()):
	"""Add gain tweak to xenon runs to match data to MC endpoint"""
	
	ADBL = AnaDBLocator()
	ADBL.req["name"] = "XeEndpt"
	ADBL.req["grouping"] = "run"
	ADBL.req["start_run"] = ADBL.req["end_run"] = "%i"%rn
	
	ddat = dict([ ((p.side,p.n),p) for p in ADBL.find(aconn) ])
	
	ADBL.req["name"] = "prevGainTweak"
	dprev = dict([ ((p.side,p.n),p) for p in ADBL.find(aconn) ])
	
	ADBL.req["source"] = os.environ["UCNA_ANA_AUTHOR"]+"_Sim"
	ADBL.req["name"] = "XeEndpt"
	dsim = dict([ ((p.side,p.n),p) for p in ADBL.find(aconn) ])
	
	print "----------",rn,"----------"
	for s in ["East","West"]:
		for t in range(4):
		
			if (s,t) not in ddat:
				print "Missing endpoint data",s,t
				continue
			if (s,t) not in dsim:
				print "Missing endpoint sim",s,t
				continue
		
			ldat = ddat[(s,t)].value
			lsim = dsim[(s,t)].value
			
			oldtweak = 1.0
			if (s,t) in dprev:
				oldtweak = dprev[(s,t)].value
			
			if ldat < 300 or lsim < 300:
				print "*** Crazy value! not uploading.",ldat,lsim
				continue
				
			print "\t",s,t,"%.2f -> %.2f\t\tErr = %+.2f%%\t(old: %+.2f%%)"%(ldat, lsim, 100.0*(lsim-ldat)/lsim, 100.*(oldtweak-1))

			if(conn):
				delete_gain_tweak(conn,rn,s,t)
				upload_gain_tweak(conn, [rn], s, t, ldat / oldtweak, lsim)
							
							
def XeTimeEvolution(rmin, rmax, yrange=(20,1000), showFits=True):
	"""Spectrum-fit isotope composition as a function of time for a series of Xe runs"""

	aconn=open_anadb_connection()
	conn=open_connection()
	
	ADBL = AnaDBLocator()
	ADBL.req["grouping"] = "run"
	ADBL.req["source"] = os.environ["UCNA_ANA_AUTHOR"]+"_Sim"
	ADBL.xcond = "start_run = end_run AND %i <= start_run AND end_run <= %i AND name LIKE 'XeComp_%%'"%(rmin,rmax)
	
	# xenon compositions
	xecomps = ADBL.find(aconn)
	
	# run timing info
	rtimes = dict([r,(getRunStartTime(conn,r),getRunEndTime(conn,r),getRunLiveTime(conn,r))] for r in range(rmin,rmax+1))
	t0 = min([t[0] for t in rtimes.values() if t[0]])
	for x in xecomps:
		x.t_mid = ((rtimes[x.start_run][0]+rtimes[x.start_run][1])*0.5-t0)/3600.
		x.dt = rtimes[x.start_run][2]
		x.rate = x.value / x.dt
		x.drate = x.err / x.dt

	# plot setup
	gIA=graph.graphxy(width=25,height=16,
					  x=graph.axis.lin(title="time [h]"),
					  y=graph.axis.log(title="event rate [Hz]",min=yrange[0], max=yrange[1]),
					  key = graph.key.key(pos="br",columns=2))
	setTexrunner(gIA)

	ks    = ['Xe135_3-2+', 'Xe125_1-2+', 'Xe133_3-2+', 'Xe131_11-2-', 'Xe129_11-2-', 'Xe133_11-2-', 'Xe137_7-2-', 'Xe135_11-2-']
	kshort = ['Xe137_7-2-', 'Xe135_11-2-']
	ksymb = [symbol.circle, symbol.triangle, symbol.square, symbol.plus, symbol.cross, symbol.diamond, symbol.circle, symbol.triangle]
	
	icols = rainbowDict(ks)
	# black-and-white version
	#for k in isotdat:
	#	icols[k] = rgb.black

	for (n,k) in enumerate(ks):
	
		isotdat = [ (x.t_mid,x.rate,x.drate) for x in xecomps if x.name[7:]==k and x.dt and x.rate > 0.1 and x.t_mid < 100 ]
		isotdat.sort()
		if not isotdat:
			continue
				
		gtitle = k.replace("_"," ")
		gtitle = "$^{"+gtitle[2:5]+"}$Xe"+gtitle[5:-1].replace("-","/")+"$^{"+gtitle[-1]+"}$"

		LF = None
		fdat = [d for d in isotdat[1:] if d[1]>1 and not (k not in kshort and d[0]<0.15) and not (k=='Xe125_1-2+' and d[0]>10)]
		if showFits and len(fdat) > 2:
			LF = LogYer(terms=[polyterm(0),polyterm(1)])
			LF.fit(fdat,cols=(0,1))
			thalf = 0
			if LF.coeffs[1]:
				thalf = -log(2)/LF.coeffs[1]

			gtitle += ": $T_{1/2}$ = "
			if abs(thalf) < 1.0:
				gtitle += "%.1f m"%(60*thalf)
			elif abs(thalf) < 24:
				gtitle += "%.1f h"%thalf
			else:
				gtitle += "%.1f d"%(thalf/24)

		sfill = []
		if k in kshort:
			sfill = [deco.filled]

		gIA.plot(graph.data.points(isotdat,x=1,y=2,dy=3,title=gtitle),
					[graph.style.symbol(ksymb[n],size=0.2,symbolattrs=[icols[k],]+sfill),
					graph.style.errorbar(errorbarattrs=[icols[k],])])
		if LF:
			gIA.plot(graph.data.points(LF.fitcurve(isotdat[0][0],isotdat[-1][0]),x=1,y=2,title=None),[graph.style.line([icols[k]])])

	gIA.writetofile(os.environ["UCNA_ANA_PLOTS"]+"/test/XeDecomp/DecompHistory_%i-%i.pdf"%(rmin,rmax))


def ep_v_eta(fname):
	
	basepath = os.environ["UCNA_ANA_PLOTS"]+"/PositionMaps/"+fname+"/"
	writepath = os.environ["UCNA_ANA_PLOTS"]+"/PositionMaps/Endpts/"
        xdat = XeFile(basepath+"/"+fname+".txt")
	
	gEp=graph.graphxy(width=20,height=10,
					  x=graph.axis.lin(title="eta",min=0,max=2.5),
					  y=graph.axis.lin(title="Endpoint (keV)"),
					  key = graph.key.key(pos="tr"))
        #setTexrunner(gEp)
	
	gHvL=graph.graphxy(width=20,height=20,
					  x=graph.axis.lin(title="Low Peak (keV)",min=0),
					  y=graph.axis.lin(title="Endpoint (keV)",min=0),
					  key = graph.key.key(pos="tl"))
	#setTexrunner(gHvL)

	
	tcols = rainbow(4)
	ssymbs = {"E":symbol.circle,"W":symbol.triangle}
	for s in ["E","W"]:
		for t in range(4):
			print s,t
                        print "Processing Endpoint data..."
			gdat = [xdat.sectdat[(s,t,m)] for m in range(xdat.sects.nSectors)]
			#        0      1        2          3              4         
			gdat = [(g.eta, g.xe_hi, g.d_xe_hi, g.xe_lo*g.eta, g.xe_hi*g.eta) for g in gdat]
			gEp.plot(graph.data.points(gdat,x=1,y=2,dy=3,title=s+"%i"%t),
					 [graph.style.symbol(ssymbs[s],size=0.2,symbolattrs=[tcols[t],]),
					  graph.style.errorbar(errorbarattrs=[tcols[t],])])
			gEp.text(gEp.width/2, gEp.height + 0.2, "Runs 19891-19898", 
       					[text.halign.center, text.valign.bottom, text.size.Large])
			gHvL.plot(graph.data.points(gdat,x=4,y=5,title=None),
					[graph.style.symbol(ssymbs[s],size=0.2,symbolattrs=[tcols[t],])])
			LF = LinearFitter(terms = [polyterm(0),polyterm(1)])
			LF.fit(gdat,cols=(3,4))
			gHvL.plot(graph.data.points(LF.fitcurve(0,250),x=1,y=2,title=s+" %i: $%s$"%(t,LF.toLatex())),
					 [graph.style.line(lineattrs=[tcols[t],])])
                        gHvL.text(gHvL.width/2, gHvL.height + 0.2, "Runs 19891-19898",
                                        [text.halign.center, text.valign.bottom, text.size.Large])

	gEp.writetofile(writepath+"/Ep_v_Eta_"+fname+".pdf")
	gHvL.writetofile(writepath+"/Hi_v_Lo_Raw_"+fname+".pdf")


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
	
	XeSets = [(17562,17650),	# subset of (17561,17734)
			  (18081,18090),
			  (18390,18413),
			  (18712,18744),
			  (19873,19898)]
	
	XeRange = XeSets[-1]
	
	XeTimeEvolution(XeRange[0], XeRange[1], yrange = (0.5,1000)); exit(0)
	
	#data_v_sim(14282,14347,12); exit(0)
	#delete_gain_tweak_range(conn,17651,17734); exit(0);
	
	#ep_v_eta("Xenon_19891-19898_12")
	#ep_v_eta("SimXe_14282-14347")
	#exit(0)
	
	conn = open_connection()
	conn = None	# to display changes without uploading
	for rn in range(XeRange[0],XeRange[1]+1):
		XeGainTweak(rn,conn)

