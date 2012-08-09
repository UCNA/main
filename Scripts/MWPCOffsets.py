#!/usr/bin/python

from QFile import *
from PyxUtils import *
from LinFitter import *
import os
from math import *

class posOffset(KVMap):
	def __init__(self,m):
		KVMap.__init__(self)
		self.dat = m.dat
		self.loadFloats(["x","y","dx","dy","d_dx","d_dy","wx","wy","d_wx","d_wy"])
		self.dr = sqrt(self.dx**2+self.dy**2);
		self.East = (self.x+0.5*self.dx,self.y+0.5*self.dy)
		self.West = (self.x-0.5*self.dx,self.y-0.5*self.dy)

if __name__ == "__main__":
	
	poffs = [posOffset(p) for p in QFile(os.environ["UCNA_ANA_PLOTS"]+"/OctetAsym_Offic/OctetAsym_Offic.txt").dat["posOffset"]]
	
	gwid = 15
	rwid = 50
	gOffs=graph.graphxy(width=gwid,height=gwid,
					   x=graph.axis.lin(title = "x [mm]",min=-rwid,max=rwid),
					   y=graph.axis.lin(title="y [mm]",min=-rwid,max=rwid),
					   key = None)
	setTexrunner(gOffs)

	gdat = [(p.x,p.y,p.dr,atan2(p.dy,p.dx)*180/pi) for p in poffs]
	lscale = gwid/(2.*rwid)*unit.v_cm
	gOffs.plot(graph.data.points(gdat,x=1,y=2,size=3,angle=4),
			   [graph.style.arrow(linelength=lscale,arrowsize=0.7*lscale,lineattrs=[style.linewidth.THIck])])

	#fdat = [ (p.East,p.x,p.y) for p in poffs]
	fdat = [ (p.West,p.x,p.y) for p in poffs]
	t0 = (lambda x: 1)
	tx = (lambda x: x[0])
	ty = (lambda x: x[1])
	LF = LinearFitter(terms = [t0,tx,ty])
	LF.fit(fdat,cols=(0,1))
	print LF.coeffs
	LF.fit(fdat,cols=(0,2))
	print LF.coeffs
	
	if False:
		gdat =[(p.x-0.5*p.dx,p.y-0.5*p.dy) for p in poffs]
		gOffs.plot(graph.data.points(gdat,x=1,y=2),[graph.style.symbol(symbol.circle)])
		gdat =[(p.x+0.5*p.dx,p.y+0.5*p.dy) for p in poffs]
		gOffs.plot(graph.data.points(gdat,x=1,y=2),[graph.style.symbol(symbol.circle)])
	
	gOffs.writetofile(os.environ["UCNA_ANA_PLOTS"]+"/test/WirechamberOffsets.pdf")
