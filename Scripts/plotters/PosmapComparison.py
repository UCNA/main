#!/usr/bin/python

from ucnacore.PyxUtils import *
from ucnacore.EncalDB import *
from ucnacore.QFile import *
from math import *
from ucnacore.LinFitter import *

class etaPoint(KVMap):
	def __init__(self,m):
		KVMap.__init__(self)
		self.dat = m.dat
		self.loadFloats(["eta","nPE","sector","tube","x","y"])
		self.loadStrings(["side",])
		self.sector = int(self.sector)
		self.tube = int(self.tube)


def loadEtaFile(rn):
	Q = QFile(os.environ["UCNA_ANA_PLOTS"]+"/PosmapDump/Posmap_%i.txt"%rn)
	etapts = {"E":[[] for t in range(5)],"W":[[] for t in range(5)]}
	for d in Q.dat["PosmapPoint"]:
		ep=etaPoint(d)
		etapts[ep.side][ep.tube].append(ep)
	return etapts

def PosmapArrowPlot(pnum):

	etapts = loadEtaFile(pnum)
	tcols = [rgb.red,rgb.blue,rgb(0.8,0,0.8),rgb(0,0.7,0)]
	
	for s in ["E","W"]:
	
		gwid = 15
		rwid = 55
		lscale = gwid/(2.*rwid)*unit.v_cm
		gPosmap=graph.graphxy(width=gwid,height=gwid,
					x=graph.axis.lin(title="x position [mm]",min=-rwid,max=rwid),
					y=graph.axis.lin(title="y position [mm]",min=-rwid,max=rwid),
					key = None)
		setTexrunner(gPosmap)
		
		for (t,tdir) in enumerate([[1.,1.],[-1.,1.],[-1.,-1.],[1.,-1.]]):
				print s,t,tdir
				peScale = 1/60.
				tdir = [x/(2*sqrt(2)) for x in tdir]
				
				gcenter = [ [ep.x,ep.y,ep.nPE*peScale,180*atan2(tdir[1],tdir[0])/3.141592] for ep in etapts[s][t] ]
				goff = [ [g[0]+tdir[0]*g[2],g[1]+tdir[1]*g[2]] + g[2:] for g in gcenter ]
				gPosmap.plot(graph.data.points(goff,x=1,y=2,size=3,angle=4),
						[graph.style.arrow(linelength=lscale,arrowsize=0.5*lscale,lineattrs=[style.linewidth.Thick,tcols[t]])])
				if t==3:
					gPosmap.plot(graph.data.points(gcenter,x=1,y=2), [graph.style.symbol(symbol.circle,size=0.07,symbolattrs=[deco.filled])])
				
				gCircs=graph.graphxy(width=gwid,height=gwid,
					x=graph.axis.lin(title="x position [mm]",min=-rwid,max=rwid),
					y=graph.axis.lin(title="y position [mm]",min=-rwid,max=rwid),
					key = None)
				setTexrunner(gCircs)
				gcenter = [ g[:2]+[g[2]*0.7] for g in gcenter]
				gCircs.plot(graph.data.points(gcenter,x=1,y=2,size=3),[varCircle()])
				gCircs.writetofile(os.environ["UCNA_ANA_PLOTS"]+"/PosmapDump/nPE_%i_%s%i.pdf"%(pnum,s,t))
				
		gPosmap.writetofile(os.environ["UCNA_ANA_PLOTS"]+"/PosmapDump/nPE_%i_%s.pdf"%(pnum,s))

def SmearCorrectingPlot(pmid):
	conn = open_connection()
	pts = getPosmapPoints(conn,pmid)

	for s in ["East","West"]:
		for t in range(4):
			if (s,t) not in pts:
				continue
			gSmear=graph.graphxy(width=15,height=10,
					x=graph.axis.lin(title="relative light transport $\\eta$"),
					y=graph.axis.lin(title="resolution correction [\\%]"),
					key = None)
			setTexrunner(gSmear)

			gdat = [ (p.sig/p.norm, 100*(p.norm-860)/860) for p in pts[(s,t)]]
			gSmear.plot(graph.data.points(gdat,x=1,y=2), [graph.style.symbol(symbol.circle,size=0.15)])

			print s,t,musigma([g[1] for g in gdat])
			gSmear.writetofile(os.environ["UCNA_ANA_PLOTS"]+"/PositionMaps/Posmap_%i/Smear_%s%i.pdf"%(pmid,s,t))


def compareEndpoints(rn1,rn2):
	
	eta1=loadEtaFile(rn1)
	eta2=loadEtaFile(rn2)
			
	for s in ['E','W']:
		for t in range(5):
			
			d1 = [ (e.sector,e) for e in eta1[s][t]]
			d1.sort()
			d2 = [ (e.sector,e) for e in eta2[s][t]]
			d2.sort()
			assert len(d1)==len(d2)
			gdat = [ [d1[n][1].nPE,d2[n][1].nPE] for n in range(len(d1))]
			
						
			# plot
			xrange = (0,1.1*max([d[1].nPE for d in d1]))
			yrange = (0,1.1*max([d[1].nPE for d in d2]))
			
			gResid=graph.graphxy(width=10,height=2,
					x=graph.axis.lin(title="Run %i nPE"%rn1,min=xrange[0],max=xrange[1]),
					y=graph.axis.lin(title="\\% Resid",min=-10,max=10))
			#gResid.texrunner.set(lfs='foils17pt')

			gComp=graph.graphxy(width=10,height=10,ypos=gResid.height+0.5,
					x=graph.axis.linkedaxis(gResid.axes["x"]),
					y=graph.axis.lin(title="Run %i nPE"%rn2,min=yrange[0],max=yrange[1]),
					key = graph.key.key(pos="tl"))
			#gComp.texrunner.set(lfs='foils17pt')
			
			cnvs = canvas.canvas()
			cnvs.insert(gResid)
			cnvs.insert(gComp)

			gComp.plot(graph.data.points(gdat,x=1,y=2,title=None), [graph.style.symbol(symbol.circle)])
			#fit comparison
			LF = LinearFitter(terms=[polyterm(1)])
			LF.fit(gdat,cols=(0,1))
			gComp.plot(graph.data.points(LF.fitcurve(xrange[0],xrange[1]),x=1,y=2,title=None),[graph.style.line([rgb.blue,])])
			
			gdat = [ g+[100.0*(g[1]-LF(g[0]))/LF(g[0]),] for g in gdat]
			gResid.plot(graph.data.points(gdat,x=1,y=3,title=None), [graph.style.symbol(symbol.circle)])
			gResid.plot(graph.data.function("y(x)=0",title=None),[graph.style.line()])
			
			# calculate discrepancy
			rms = sqrt(sum([g[2]**2 for g in gdat])/len(gdat))
			print s,t,"RMS =",rms

			gComp.text(7.5,1+gResid.height+0.5,"%s%i"%(s,t+1))
			gComp.text(7.5,0.3+gResid.height+0.5,"RMS = %.1f\\%%"%rms)
		
			
			outdir = os.environ["UCNA_ANA_PLOTS"]+"/nPE/Comparisons/%i_%i/"%(rn1,rn2)
			os.system("mkdir -p "+outdir)
			cnvs.writetofile(outdir+"%s_%i.pdf"%(s,t))


class CircInterp:
	def __init__(self,dpts):
		self.dat = dpts

	def __call__(self,th):
		while th>=2*pi:
			th -= 2*pi
		n = len(self.dat)
		i = int(len(self.dat)*th/(2*pi))%n
		y = n*th/(2*pi)-i
		p0 = self.dat[i-1]
		p1 = self.dat[i]
		p2 = self.dat[(i+1)%n]
		p3 = self.dat[(i+2)%n]
		A = -0.5
		return ( A*p0*(1-y)*(1-y)*y + p1*(1-y)*(1-y*((2+A)*y-1)) - p2*y*(A*(1-y)*(1-y)+y*(2*y-3)) + A*p3*(1-y)*y*y );

def sc_rings():
	cumdivs = [0,1]
	ndivs = [1]
	for i in range(20)[1:]:
		ndivs.append(int(ceil(2*pi*i)))
		cumdivs.append(cumdivs[-1]+ndivs[-1])
	return cumdivs

def interpPlot():

	rn = 16000
	eta = loadEtaFile(rn)
	gInterp=graph.graphxy(width=15,height=10,
					x=piaxis(title="$\\theta$",min=0,max=2*pi),
					y=graph.axis.lin(title="eta"))
	setTexrunner(gInterp)
	
	rngs = sc_rings()
	print rngs
	for i in range(len(rngs)-1):
		pts = [ (p.sector, p.nPE) for p in eta["W"][2] if rngs[i] <= p.sector <= rngs[i+1]-1 ]
		pts.sort()
		if not pts:
			continue
			
		CI = CircInterp([p[1] for p in pts])
		

		gdat = [(x,CI(x)) for x in unifrange(0,2*pi,400)]
		gInterp.plot(graph.data.points(gdat,x=1,y=2,title=None), [graph.style.line()])
		pts.append(pts[0])
		gdat = [(float(n)*2*pi/(len(pts)-1),p[1]) for (n,p) in enumerate(pts)]
		gInterp.plot(graph.data.points(gdat,x=1,y=2,title=None), [graph.style.symbol(symbol.circle,size=0.1,symbolattrs=[deco.filled([rgb.white])])])
		
	gInterp.writetofile(os.environ["UCNA_ANA_PLOTS"]+"/test/Interpolation/Posmap_%i.pdf"%rn)

if __name__=="__main__":
	
	#PosmapArrowPlot(16000)
	#SmearCorrectingPlot(164)
	interpPlot()
	#compareEndpoints(14300,16000)
	#compareEndpoints(16000,18362)
	#compareEndpoints(18362,18750)
	#compareEndpoints(18081,18362)