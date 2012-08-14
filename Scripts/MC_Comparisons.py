#!/usr/bin/python

from QFile import *
from PyxUtils import *
from math import *
import os

class enhist(KVMap):
	def __init__(self,m=KVMap()):
		self.dat = m.dat
		self.loadFloats(["counts","mean","normcounts","rms","sim","type"])
		self.type = int(self.type)
		self.sim = int(self.sim)

class MCCompFile(QFile):
	def __init__(self,fname):
		QFile.__init__(self,fname)
		
		self.evis = {}
		for m in [enhist(m) for m in self.dat.get("evis",[])]:
			self.evis[(m.sim,m.type)] = m

		self.equench = {}
		for m in [enhist(m) for m in self.dat.get("equench",[])]:
			self.equench[(m.sim,m.type)] = m




def G4_Pen_Compare():
	basepath = outpath = os.environ["UCNA_ANA_PLOTS"]+"/test/MC_Compare/"
	mccols = {"Geant4":rgb.red,"Penelope":rgb.blue}
	mcsymbs = {"Geant4":symbol.circle,"Penelope":symbol.triangle}
	tpname = {0:"0",1:"I",2:"II-III"}
	for tp in range(3):
		print "Type",tp,"events"
		gdat = {}
		for mc in range(2):
			for l in [50,100,150,200,300,400,600,800]:
				mcf = MCCompFile(basepath+"%i_keV/MC_Compare.txt"%l)
				if not mcf.evis:
					continue
				counts = mcf.evis[(mc,tp)].counts
				if not counts:
					counts = 1
				ncounts = mcf.evis[(mc,tp)].normcounts
				euncert = mcf.evis[(mc,tp)].rms/sqrt(counts)
				equncert = mcf.equench[(mc,tp)].rms/sqrt(counts)
				evis = mcf.evis[(mc,tp)].mean
				d = (l,l-evis,euncert,100*ncounts,100*ncounts/sqrt(counts),evis-mcf.equench[(mc,tp)].mean,equncert)
				gdat.setdefault(mc,[]).append(d)

		gEvis=graph.graphxy(width=10,height=10,
						  x=graph.axis.log(title="primary energy [keV]",min=40,max=1000),
						  y=graph.axis.lin(title="mean energy loss [keV]",min=0),
						  key = graph.key.key(pos="bl"))
		setTexrunner(gEvis)
		
		gEquench=graph.graphxy(width=10,height=10,
						x=graph.axis.log(title="primary energy [keV]",min=40,max=1000),
						y=graph.axis.lin(title="quenching loss [keV]",min=0),
						key = graph.key.key(pos="bl"))
		setTexrunner(gEquench)
		
		gScatter=graph.graphxy(width=10,height=10,
							x=graph.axis.log(title="primary energy [keV]",min=40,max=1000),
							y=graph.axis.lin(title="Type %s backscatter [\\%% of Type 0]"%tpname[tp],min=0,max=5),
							key = graph.key.key(pos="tl"))
		setTexrunner(gScatter)
		
		for (n,mc) in enumerate(["Geant4","Penelope"]):
			gEvis.plot(graph.data.points([g for g in gdat[n] if g[3]],x=1,y=2,dy=3,title=mc),
					 [graph.style.symbol(mcsymbs[mc],size=0.2,symbolattrs=[mccols[mc],]),
					  graph.style.errorbar(errorbarattrs=[mccols[mc],])])
			gScatter.plot(graph.data.points(gdat[n],x=1,y=4,dy=5,title=mc),
					   [graph.style.symbol(mcsymbs[mc],size=0.2,symbolattrs=[mccols[mc],]),
						graph.style.errorbar(errorbarattrs=[mccols[mc],])])
			gEquench.plot(graph.data.points([g for g in gdat[n] if g[3]],x=1,y=6,dy=7,title=mc),
						 [graph.style.symbol(mcsymbs[mc],size=0.2,symbolattrs=[mccols[mc],]),
						  graph.style.errorbar(errorbarattrs=[mccols[mc],])])
					
		gEvis.writetofile(basepath+"/Evis_Type_%i.pdf"%tp)
		gEquench.writetofile(basepath+"/Equench_Type_%i.pdf"%tp)
		gScatter.writetofile(basepath+"/Frac_Type_%i.pdf"%tp)


if __name__=="__main__":
	G4_Pen_Compare()
