#!/usr/bin/python

from Asymmetries import *



def compare_anchoices(datf, simf):
	gdat = []
	ancs = ["A","B","C","D","E","F","G","H","I","J","K"]
	
	for (n,a) in enumerate(ancs):
		aDat = AsymmetryFile(datf+"/Anchoice_"+a+"/Anchoice_"+a+".txt").getAsym(200,675)
		aSim = AsymmetryFile(simf+"/Anchoice_"+a+"/Anchoice_"+a+".txt").getAsym(200,675)
		gdat.append([n+1,aDat.A0,aDat.dA0,aSim.A0,aSim.dA0,aDat.A0/aSim.A0-1.,aDat.dA0/abs(aSim.A0)])

			
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
	gdA=graph.graphxy(width=10,height=10,
					 x=graph.axis.lin(title="Analysis Choice",manualticks=myticks,min=0.5,max=len(ancs)+0.5,parter=None),
					 y=graph.axis.lin(title="${\\rm Data}/{\\rm MC}-1$ Asymmetry Difference",min=0,max=0.025),
					 key = graph.key.key(pos="tl"))
	setTexrunner(gdA)
	gdA.plot(graph.data.points(gdat,x=1,y=6,dy=7,title=None),
			[graph.style.symbol(symbol.triangle,size=0.2,symbolattrs=[rgb.blue,]),
			 graph.style.errorbar(errorbarattrs=[rgb.blue,])])
	
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
	
	ancs = ["A","B","C","D","E","F","G","H","I","J","K"]
	ancs = ancs[:4]
	myticks = [ graph.axis.tick.tick(n+1,label=a) for (n,a) in enumerate(ancs) ]	
	gdA=graph.graphxy(width=10,height=10,
					  x=graph.axis.lin(title="Analysis Choice",manualticks=myticks,min=0.5,max=len(ancs)+0.5,parter=None),
					  y=graph.axis.lin(title="Extracted $A_0$",min=-.124,max=-.118),
					  key = graph.key.key(pos="bl"))
	setTexrunner(gdA)

	gdat = [(f.n,f.A0,f.dA0) for f in rawdat]
	gdA.plot(graph.data.points(gdat,x=1,y=2,dy=3,title="Raw"),
			 [graph.style.symbol(symbol.triangle,size=0.2,symbolattrs=[rgb.blue,]),
			  graph.style.errorbar(errorbarattrs=[rgb.blue,])])
	gdat = [(f.n,f.A0,f.dA0) for f in corrdat]
	gdA.plot(graph.data.points(gdat,x=1,y=2,dy=3,title="Corrected"),
			 [graph.style.symbol(symbol.circle,size=0.2,symbolattrs=[rgb.red,]),
			  graph.style.errorbar(errorbarattrs=[rgb.red,])])

	gdA.writetofile(basedir+"/Asymmetries.pdf")

if __name__=="__main__":
	#compare_anchoices(os.environ["UCNA_ANA_PLOTS"]+"/test/Anchoices/",os.environ["UCNA_ANA_PLOTS"]+"/test/Anchoices_SimMagF/")
	compare_corrections(os.environ["UCNA_ANA_PLOTS"]+"/test/CorrectAsym_SimPen/","CorrectedAsym.txt")
	compare_corrections(os.environ["UCNA_ANA_PLOTS"]+"/test/CorrectAsym/","CorrectedAsym.txt")