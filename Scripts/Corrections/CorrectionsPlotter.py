#!/usr/bin/python

from math import *
import os
import sys
sys.path.append("..")
from PyxUtils import *

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
	
	gCx=graph.graphxy(width=20,height=12,
						  x=graph.axis.lin(title="Energy [keV]",min=0,max=800),
						  y=graph.axis.lin(title="Correction [\\%]",min=-0.05,max=0.2),
						  key = graph.key.key(pos="tl"))
	setTexrunner(gCx)

	cxcols = rainbowDict(cxns)
	for cx in cxns:
		gdat = [ [0.5*(d[0]+d[1]),100*d[2],abs(100*d[3])] for d in cxns[cx].dat]
		gCx.plot(graph.data.points(gdat,x=1,y=2,dy=3,title=cx),
				[graph.style.symbol(symbol.circle,size=0.1,symbolattrs=[cxcols[cx]]),
				graph.style.errorbar(errorbarattrs=[cxcols[cx]])])

	gCx.writetofile(os.environ["UCNA_ANA_PLOTS"]+"/test/AsymCorrections.pdf")
				 

def PlotUncerts():
	
	clist = { None:"NGBG.txt" }
	cxns = dict([(k,CorrFile(baseCorrPath+clist[k])) for k in clist])
	
	gCx=graph.graphxy(width=20,height=12,
					  x=graph.axis.lin(title="Energy [keV]",min=0,max=800),
					  y=graph.axis.lin(title="Error in A [\\%]",min=-0.1,max=0.1),
					  key = graph.key.key(pos="tl"))
	setTexrunner(gCx)
	
	cxcols = rainbowDict(cxns)
	for cx in cxns:
		gdat = [ [0.5*(d[0]+d[1]),100*d[3]] for d in cxns[cx].dat]
		gCx.plot(graph.data.points(gdat,x=1,y=2,title=cx),
				 [graph.style.symbol(symbol.circle,size=0.1,symbolattrs=[cxcols[cx]]),
				  graph.style.errorbar(errorbarattrs=[cxcols[cx]])])
	
	gCx.writetofile(os.environ["UCNA_ANA_PLOTS"]+"/test/AsymUncerts.pdf")

if __name__=="__main__":
	PlotCorrections()
	#PlotUncerts()
