#!/usr/bin/python

import sys
sys.path.append("..")

from review.Asymmetries import *
from ucnacore.PyxUtils import *
from ucnacore.DecayPhysics import *
import copy

class AsymCorrFile(QFile):
	def __init__(self,fname):
		QFile.__init__(self,fname)
		
		self.cxns = {}
		for p in self.dat.get("corrPoint",[]):
			if p.getFirst("anChoice") != 'C':
				continue
			p.loadFloats(["KE","cTot","cErr"])
			self.cxns.setdefault("tot",[]).append((p.KE,p.cTot,p.cErr))
			for k in p.dat:
				if p.dat.has_key("d_"+k):
					self.cxns.setdefault(k,[]).append((p.KE,p.getFirstF(k),p.getFirstF("d_"+k)))

def compare_anchoices(datf, simf):
	gdat = []
	ancs = ["A","B","C","D","E","F","G","H","I","J","K"]
	
	pol_corr = 1.0067 # polarization correction
	
	for (n,a) in enumerate(ancs):
		aDat = AsymmetryFile(datf+"/Anchoice_"+a+"/Anchoice_"+a+".txt").getAsym(225,675)
		aSim = AsymmetryFile(simf+"/Anchoice_"+a+"/Anchoice_"+a+".txt").getAsym(225,675)
		gdat.append([n+1,aDat.A0*pol_corr,aDat.dA0*pol_corr,aSim.A0,aSim.dA0,aDat.A0*pol_corr*A0_PDG/aSim.A0,aDat.dA0*pol_corr*abs(A0_PDG/aSim.A0)])

			
	myticks = [ graph.axis.tick.tick(n+1,label=a) for (n,a) in enumerate(ancs) ]

	gA=graph.graphxy(width=10,height=10,
				  x=graph.axis.lin(title="Analysis Choice",manualticks=myticks,min=0.5,max=len(ancs)+0.5,parter=None),
				  y=graph.axis.lin(title="Asymmetry"),
				  key = graph.key.key(pos="tl"))
	setTexrunner(gA)

	gA.plot(graph.data.points(gdat,x=1,y=4,dy=5,title="MC"),
			[graph.style.symbol(symbol.triangle,size=0.2,symbolattrs=[rgb.blue,]),
			graph.style.errorbar(errorbarattrs=[rgb.blue,])])
	gA.plot(graph.data.points(gdat,x=1,y=2,dy=3,title="Data"),
		 [graph.style.symbol(symbol.circle,size=0.2,symbolattrs=[rgb.red,]),
		  graph.style.errorbar(errorbarattrs=[rgb.red,])])

	gA.writetofile(simf+"/Anchoices.pdf")

	ancs = ancs[:5]
	myticks = myticks[:5]
	gdat = gdat[:5]
	#ancs = ancs[5:]
	#myticks = myticks[5:]
	#gdat = gdat[5:]
	gdA=graph.graphxy(width=10,height=10,
					 x=graph.axis.lin(title="Analysis Choice",manualticks=myticks,min=0.5,max=len(ancs)+0.5,parter=None),
					 y=graph.axis.lin(title="Asymmetry",min=-0.125,max=-0.118),
					 key = graph.key.key(pos="bl"))
	setTexrunner(gdA)
	gdA.plot(graph.data.points(gdat,x=1,y=2,dy=3,title="Uncorrected"),
			[graph.style.symbol(symbol.triangle,size=0.2,symbolattrs=[deco.filled(),]),
			 graph.style.errorbar(errorbarattrs=[])])
	gdA.plot(graph.data.points(gdat,x=1,y=6,dy=7,title="MC Corrected"),
			[graph.style.symbol(symbol.circle,size=0.2,symbolattrs=[deco.filled(),]),
			 graph.style.errorbar(errorbarattrs=[])])
	gdA.plot(graph.data.function("y(x)=%f"%A0_PDG,title=None),
			[graph.style.line(lineattrs=[style.linestyle.dashed,])])
			 
	gdA.writetofile(simf+"/Anchoices_Delta.pdf")


class anchoiceFit(KVMap):
	def __init__(self,m=KVMap()):
		self.dat = m.dat
		self.loadFloats(["A0","dA0","chi2","ndf"])
		self.loadStrings(["anChoice"])
		self.n = {'A':1,'B':2,'C':3,'D':4,'E':5,'F':6,'G':7,'H':8,'I':9,'J':10,'K':11}[self.anChoice]

def compare_corrections(basedir,fin):
	q = QFile(basedir+"/"+fin)
	rawdat = [anchoiceFit(f) for f in q.dat["rawFit"]]
	corrdat = [anchoiceFit(f) for f in q.dat["correctedFit"]]
	
	ancs = ["A","D","B","C","E","F","G","H","I","J","K"]
	ancs = ancs[1:4]
	anclabs = {"B":"0+I","C":"0+I+II/III","D":"0"}
	myticks = [ graph.axis.tick.tick(n+1,label=anclabs[a]) for (n,a) in enumerate(ancs) ]
	gdA=graph.graphxy(width=10,height=10,
					  x=graph.axis.lin(title="Analysis Choice",manualticks=myticks,min=0.5,max=len(ancs)+0.5,parter=None),
					  y=graph.axis.lin(title="Extracted $A_0$",min=-.125,max=-.115),
					  key = graph.key.key(pos="bl"))
	setTexrunner(gdA)

	gdat = [(f.n,f.A0,f.dA0) for f in rawdat]
	gdA.plot(graph.data.points(gdat,x=1,y=2,dy=3,title="Raw"),
			 [graph.style.symbol(symbol.triangle,size=0.2,symbolattrs=[rgb.blue,]),
			  graph.style.errorbar(errorbarattrs=[rgb.blue,])])
	gdat = [(f.n,f.A0,f.dA0) for f in corrdat]
	gdA.plot(graph.data.points(gdat,x=1,y=2,dy=3,title="MC Corrected"),
			 [graph.style.symbol(symbol.circle,size=0.2,symbolattrs=[rgb.red,]),
			  graph.style.errorbar(errorbarattrs=[rgb.red,])])

	gdA.writetofile(basedir+"/Asymmetries.pdf")

	A0s = [g[1] for g in gdat]
	mu,sigma = mu_sigma(A0s)
	dmax = max([abs(a-mu) for a in A0s])
	print "mu,sigma,dmax=",(mu,sigma,dmax)

def plot_A_corrs(basedir,fin):
	
	AC = AsymCorrFile(basedir+"/"+fin)
	
	ygridpainter = graph.axis.painter.regular(gridattrs=[attr.changelist([style.linestyle.dashed, style.linestyle.dotted])])
	
	for ctp in ["Delta_2","Delta_3"]:
		
		ccol = ctp
		ctit = "Backscatter"
		if ctp=="Delta_3":
			ccol += "_C"
			ctit = "Angle"
		
		gdA=graph.graphxy(width=12,height=8,
					   x=graph.axis.lin(title="Energy [keV]",min=0,max=800),
					   y=graph.axis.lin(title=ctit+" Correction [\\%]",min=-5,max=5,painter=ygridpainter),
					   key = graph.key.key(pos="bl"))
		setTexrunner(gdA)
		
		
		gdat = [(p[0],100.*p[1],abs(100.*p[2])) for p in AC.cxns[ccol] if 75 < p[0] < 780]
		gdat.sort()
		area = errorBand(gdA,gdat,0,1,2)
		gdA.fill(area, [deco.filled([color.rgb(0.0,1.0,0.0),color.transparency(0.5)])])
		
		gdA.plot(graph.data.points(gdat,x=1,y=2,dy=3,title=None),
				 [graph.style.line(lineattrs=[style.linewidth.THick]), ])

		gdA.writetofile(basedir+"/"+ctp+".pdf")


if __name__=="__main__":
	
	#plot_A_corrs(os.environ["UCNA_ANA_PLOTS"]+"/OctetAsym_Offic/Range_0-1000/CorrectAsym/","CorrectedAsym.txt")
	#plot_A_corrs(os.environ["UCNA_ANA_PLOTS"]+"/test/MCCors/thinfoil/","test.txt")
	#exit(0)
	
	compare_anchoices(os.environ["UCNA_ANA_PLOTS"]+"/test/Anchoices/",os.environ["UCNA_ANA_PLOTS"]+"/test/Anchoices_SimMagF/")

	#compare_corrections(os.environ["UCNA_ANA_PLOTS"]+"/test/CorrectAsym_Sim0823_4x/","CorrectedAsym.txt")
	#compare_corrections(os.environ["UCNA_ANA_PLOTS"]+"/OctetAsym_Offic/Range_0-16/CorrectAsym/","CorrectedAsym.txt")
	#compare_corrections(os.environ["UCNA_ANA_PLOTS"]+"/OctetAsym_Offic/Range_17-1000/CorrectAsym/","CorrectedAsym.txt")
	#compare_corrections(os.environ["UCNA_ANA_PLOTS"]+"/OctetAsym_Offic/Range_0-1000/CorrectAsym/","CorrectedAsym.txt")

	#compare_corrections(os.environ["UCNA_ANA_PLOTS"]+"/test/CorrectAsym_SimPen/","CorrectedAsym.txt")

	
