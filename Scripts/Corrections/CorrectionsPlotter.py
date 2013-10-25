#!/usr/bin/python

from math import *
import os
import sys
sys.path.append("..")
from ucnacore.PyxUtils import *
from EnergyErrors import *

class bg_graphxy(graph.graphxy):
    def path_between(self, f, a=None, b=None):
        """Return the path f between x=a and x=b"""
        if a == None: a = self.axes['x'].axis.min
        if b == None: b = self.axes['x'].axis.max
        xa = self.xgridpath(a)
        xb = self.xgridpath(b)
        f_at_xa = f.intersect(xa)[0][0]
        f_at_xb = f.intersect(xb)[0][0]
        return f.split([f_at_xa, f_at_xb])[1]

    def area_between(self, f0, f1, a=None, b=None):
        """Return the area between the two plot paths f0 and f1, from x=a to
        x=b."""
        f0 = self.path_between(f0, a, b)
        f1 = self.path_between(f1, a, b).reversed()
        return f0 << f1

    def move_to_back(self, n=1):
        """Cycle the previous n items to the bottom (background)"""
        self.items = self.items[-n:] + self.items[:-n]


def coerce(y,ymn,ymx):
	if y > ymx:
		return ymx
	if y < ymn:
		return ymn
	return y

def errorband_area(g,gdat,cols=(0,1,2)):

	ymx = g.axes['y'].axis.max
	ymn = g.axes['y'].axis.min
	xmn = g.axes['x'].axis.min
	xmx = g.axes['x'].axis.max
	gup = [(p[cols[0]],coerce(p[cols[1]]+p[cols[2]],ymn,ymx)) for p in gdat]
	gdn = [(p[cols[0]],coerce(p[cols[1]]-p[cols[2]],ymn,ymx)) for p in gdat]
	
	xlo = g.xgridpath(xmn)
	xhi = g.xgridpath(xmx)
	ylo = g.ygridpath(ymn)
	yhi = g.ygridpath(ymx)
	
	xlo = xlo.split([xlo.intersect(ylo)[0][0],xlo.intersect(yhi)[0][0]])[1].reversed()
	xhi = xhi.split([xhi.intersect(ylo)[0][0],xhi.intersect(yhi)[0][0]])[1].reversed()
	
	u = g.plot(graph.data.points(gup,x=1,y=2,title=None),[graph.style.line(None)])
	d =	g.plot(graph.data.points(gdn,x=1,y=2,title=None),[graph.style.line(None)])
	g.finish()
	#return g.area_between(u.path,d.path)
	area = u.path << xhi << d.path.reversed() << xlo.reversed()
	#area = d.path.reversed() << u.path
	area.append(path.closepath())
	#area[-1].close()
	return area


class CorrFile:
	def __init__(self,fname):
		self.dat = [[float(x) for x in l.split()] for l in open(fname,"r").readlines() if len(l)>10 and l[0] != '#']
baseCorrPath = "../../Aux/Corrections/"

def PlotCorrections():
	
	clist = {
			None:"NGBG.txt",
		#"Muon Veto":"MuonEffic.txt",
		#	"Linearity":"EnergyLinearityUncertainty_2010.txt",
		#	"Gain Flucts":"GainFlucts.txt",
		#	"Ped Shifts":"PedShifts.txt",
		#	"Recoil Order":"RecoilOrder.txt",
		#	"Radiative":"Radiative_h-g.txt"
			}
	cxns = dict([(k,CorrFile(baseCorrPath+clist[k])) for k in clist])
	
	gCx=bg_graphxy(width=20,height=12,
						  x=graph.axis.lin(title="Energy [keV]",min=0,max=800),
						  y=graph.axis.lin(title="Uncertainty on $A$ [\\%]",min=-0.05,max=0.05),
						  key = graph.key.key(pos="tl"))
	setTexrunner(gCx)

	cxcols = rainbowDict(cxns)
	for cx in cxns:
		gdat = [ [0.5*(d[0]+d[1]),100*d[2],abs(100*d[3])] for d in cxns[cx].dat][2:]
		gCx.plot(graph.data.points(gdat,x=1,y=2,dy=3,title=cx),[graph.style.line([style.linewidth.THick])])
		gCx.plot(graph.data.function("y(x)=0",title=None),[graph.style.line([style.linewidth.thick,style.linestyle.dashed])])
		eband = errorband_area(gCx,gdat)
		gCx.fill(eband,[pattern.hatched(0.1,45)])


	gCx.writetofile(os.environ["UCNA_ANA_PLOTS"]+"/test/NGBGCorrections.pdf")
				 

def PlotUncerts():
	
	clist = { None:"MuonEffic.txt"}
	cxns = dict([(k,CorrFile(baseCorrPath+clist[k])) for k in clist])
	
	gCx=graph.graphxy(width=20,height=12,
					  x=graph.axis.lin(title="Energy [keV]",min=0,max=800),
					  y=graph.axis.lin(title="Uncertainty in A [\\%]",min=0,max=0.5),
					  key = graph.key.key(pos="tl"))
	setTexrunner(gCx)
	
	cxcols = {None:rgb.black} #rainbowDict(cxns)
	for cx in cxns:
		gdat = [ [0.5*(d[0]+d[1]),100*d[3]] for d in cxns[cx].dat if 20 < d[0] < 700]
		print gdat
		gCx.plot(graph.data.points(gdat,x=1,y=2,title=cx),
				 [
					#graph.style.symbol(symbol.circle,size=0.1,symbolattrs=[cxcols[cx]]),
				  	graph.style.line([style.linewidth.THick]),
					#graph.style.errorbar(errorbarattrs=[cxcols[cx]])])
				  ])
	
	gCx.writetofile(os.environ["UCNA_ANA_PLOTS"]+"/test/AsymUncerts.pdf")




if __name__=="__main__":
	PlotCorrections()
	#PlotUncerts()

