#!/usr/bin/python

from QFile import *
from PyxUtils import *
from LinFitter import *
from Asymmetries import *
import os
from math import *

class posOffset(KVMap):
	def __init__(self,m):
		KVMap.__init__(self)
		self.dat = m.dat
		self.loadFloats(["x","y","dx","dy","d_dx","d_dy","wx","wy","d_wx","d_wy"])
		self.calcRTheta()
		self.East = (self.x+0.5*self.dx,self.y+0.5*self.dy)
		self.West = (self.x-0.5*self.dx,self.y-0.5*self.dy)

	def calcRTheta(self):
		self.dr = self.dr = sqrt(self.dx**2+self.dy**2)
		self.theta = atan2(self.dy,self.dx)
			
	def applyFit(self,LFx,LFy):
		self.dx = 2*(LFx((self.x,self.y))-self.x)
		self.dy = 2*(LFy((self.x,self.y))-self.y)
		self.calcRTheta()

def fitOffsets(poffs):
	
	#fdat = [ (p.East,p.x,p.y) for p in poffs]
	fdat = [ (p.West,p.x,p.y) for p in poffs]
	
	t0 = (lambda x: 1)
	tx = (lambda x: x[0])
	ty = (lambda x: x[1])
	
	LFx = LinearFitter(terms = [t0,tx,ty])
	LFx.fit(fdat,cols=(0,1))
	xShift = LFx.coeffs[0]
	sinTh = LFx.coeffs[2]
	
	LFy = LinearFitter(terms = [t0,tx,ty])
	LFy.fit(fdat,cols=(0,2))
	yShift = LFy.coeffs[0]
	sinTh = 0.5*(sinTh-LFy.coeffs[1])
	
	LFx.coeffs[1] = sqrt(1-sinTh**2)
	LFx.coeffs[2] = sinTh
	
	LFy.coeffs[1] = -sinTh
	LFy.coeffs[2] = sqrt(1-sinTh**2)
	
	return (xShift,yShift,sinTh,LFx,LFy)
	
def offsetPlot(basedir,fname):
	
	poffs = [posOffset(p) for p in QFile(basedir+"/"+fname).dat["posOffset"]]
	(xShift,yShift,sinTh,LFx,LFy) = fitOffsets(poffs)
	
	print "Average offset:",(xShift,yShift),sinTh
	
	gwid = 15
	rwid = 50
	lscale = gwid/(2.*rwid)*unit.v_cm
	gOffs=graph.graphxy(width=gwid,height=gwid,
						x=graph.axis.lin(title = "x [mm]",min=-rwid,max=rwid),
						y=graph.axis.lin(title="y [mm]",min=-rwid,max=rwid),
						key = None)
	setTexrunner(gOffs)
	
	gdat = [(p.x,p.y,p.dr,p.theta*180/pi) for p in poffs]
	
	for p in poffs:
		p.applyFit(LFx,LFy);
	fdat = [(p.x,p.y,p.dr,p.theta*180/pi) for p in poffs]
	
	gOffs.plot(graph.data.points(fdat,x=1,y=2,size=3,angle=4),
			   [graph.style.arrow(linelength=lscale,arrowsize=0.7*lscale,lineattrs=[style.linewidth.THIck,rgb.red])])
	gOffs.plot(graph.data.points(gdat,x=1,y=2,size=3,angle=4),
			   [graph.style.arrow(linelength=lscale,arrowsize=0.7*lscale,lineattrs=[style.linewidth.THIck])])
	
	
	gOffs.writetofile(basedir+"/MWPCOffsets.pdf")

def offsetHistory(basedir,depth=0):

	shiftHist = []
	n = 0
	for AF in collectAsymmetries(basedir,depth):
		
		print n,AF.getRuns()
		
		poffs = [posOffset(p) for p in AF.dat["posOffset"]]
		(xShift,yShift,sinTh,LFx,LFy) = fitOffsets(poffs)
		shiftHist.append((n,xShift,yShift,sinTh))
		n += 1
		
		if False:
			gdat =[(p.x-0.5*p.dx,p.y-0.5*p.dy) for p in poffs]
			gOffs.plot(graph.data.points(gdat,x=1,y=2),[graph.style.symbol(symbol.circle)])
			gdat =[(p.x+0.5*p.dx,p.y+0.5*p.dy) for p in poffs]
			gOffs.plot(graph.data.points(gdat,x=1,y=2),[graph.style.symbol(symbol.circle)])
	
	unitName=unitNames[depth]
	gShift=graph.graphxy(width=15,height=10,
						 x=graph.axis.lin(title=unitName,min=0,max=shiftHist[-1][0]),
						 y=graph.axis.lin(title="Detector offset [mm]"),
						 key = graph.key.key(pos="tr",columns=2))
	setTexrunner(gShift)
	gRot=graph.graphxy(width=15,height=10,
					   x=graph.axis.lin(title=unitName,min=0,max=shiftHist[-1][0]),
					   y=graph.axis.lin(title="Detector rotation [radians]"),
					   key = None)
	setTexrunner(gRot)
	
	
	gShift.plot(graph.data.points(shiftHist,x=1,y=2,title="x"),[graph.style.symbol(symbol.cross,symbolattrs=[rgb.red])])
	gShift.plot(graph.data.points(shiftHist,x=1,y=3,title="y"),[graph.style.symbol(symbol.circle,symbolattrs=[rgb.blue])])
	gShift.writetofile(basedir+"/MWPCShift.pdf")
	
	gRot.plot(graph.data.points(shiftHist,x=1,y=4),[graph.style.symbol(symbol.circle)])
	gRot.writetofile(basedir+"/MWPCRot.pdf")


if __name__ == "__main__":
	
	offsetPlot(os.environ["UCNA_ANA_PLOTS"]+"/test/SimReplay/Test_neutronBetaUnpol/","Test_neutronBetaUnpol.txt")
	
	if 0:
		basedir = os.environ["UCNA_ANA_PLOTS"]+"/OctetAsym_Offic/"
		
		offsetPlot(basedir,"OctetAsym_Offic.txt")
		offsetHistory(basedir)
	