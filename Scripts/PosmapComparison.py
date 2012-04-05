#!/usr/bin/python

from PyxUtils import *
from EncalDB import *
from QFile import *
from math import *
from LinFitter import *

class etaPoint(KVMap):
	def __init__(self,m):
		KVMap.__init__(self)
		self.dat = m.dat
		self.loadFloats(["eta","nPE","sector","t","x","y"])
		self.loadStrings(["side",])
		self.sector = int(self.sector)
		self.t = int(self.t)

def loadEtaFile(rn):
	Q = QFile(os.environ["UCNA_ANA_PLOTS"]+"/nPE/Run_%i/NPE.txt"%rn)
	etapts = {"E":[[] for t in range(5)],"W":[[] for t in range(5)]}
	for d in Q.dat["pmap"]:
		ep=etaPoint(d)
		etapts[ep.side][ep.t].append(ep)
	return etapts

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

if __name__=="__main__":
	compareEndpoints(14300,16000)
	#compareEndpoints(16000,18362)
	#compareEndpoints(18362,18750)
	#compareEndpoints(18081,18362)