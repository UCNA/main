#!/usr/bin/python

from QFile import *
from PyxUtils import *
from pyx import pattern
from math import *
from LinFitter import *
import os

class enhist(KVMap):
	def __init__(self,m=KVMap()):
		self.dat = m.dat
		self.loadFloats(["counts","mean","normcounts","rms","sim","type"])
		self.type = int(self.type)
		self.sim = int(self.sim)

class asymplot(KVMap):
	def __init__(self,m=KVMap()):
		self.dat = m.dat
		self.loadFloats(["A0","dA0","corr","dcorr","ssbg","ssfg","ssMC","KE"])

class MCCompFile(QFile):
	def __init__(self,fname):
		QFile.__init__(self,fname)
		
		self.evis = {}
		for m in [enhist(m) for m in self.dat.get("evis",[])]:
			self.evis[(m.sim,m.type)] = m

		self.equench = {}
		for m in [enhist(m) for m in self.dat.get("equench",[])]:
			self.equench[(m.sim,m.type)] = m


def errorBand(gdat,nx,ny,ndy):
	gwid = 12
	g0 = graph.graphxy(width=gwid,height=0.2*gwid,
				    x=graph.axis.lin(title="Reconstructed Energy [keV]",min=0,max=800),
				    y=graph.axis.lin(title="Decay Rate [Hz/keV]",min=-5,max=5),)
	d1 = g0.plot(graph.data.points([(g[nx],g[ny]-g[ndy]) for g in gdat],x=1,y=2), [graph.style.line()])
	d2 = g0.plot(graph.data.points([(g[nx],g[ny]+g[ndy]) for g in gdat],x=1,y=2), [graph.style.line()])
	g0.finish()
	p1 = d1.path
	p2 = d2.path
	
	area = (p1 << p2.reversed())
	#area[-1].close()

	return area



def paperFig1(inColor=False):
	pdat = [asymplot(a) for a in QFile(os.environ["UCNA_ANA_PLOTS"]+"/Paper/PlotData.txt").dat["asymplot"]]

	xrange = (0,800)
	gwid = 12
	
	gCorr=graph.graphxy(width=gwid,height=0.2*gwid,
							  x=graph.axis.lin(title="Reconstructed Energy [keV]",min=xrange[0],max=xrange[1]),
							  y=graph.axis.lin(title="\\% Corr.",min=-5,max=5),
						key = graph.key.key(pos="tl"))
	setTexrunner(gCorr)
	
	gAsym=graph.graphxy(width=gwid,height=0.4*gwid,ypos=gCorr.height+0.5,
						x=graph.axis.linkedaxis(gCorr.axes["x"]),
						y=graph.axis.lin(title="$A_0$",min=-0.14,max=-0.10),
						key = graph.key.key(pos="tr"))
	setTexrunner(gAsym)
	
	gSpec=graph.graphxy(width=gwid,height=0.5*gwid,ypos=gCorr.height+gAsym.height+1.0,
							 x=graph.axis.linkedaxis(gCorr.axes["x"]),
							 y=graph.axis.lin(title="Decay Rate [Hz/keV]",min=0,max=0.065),
							 key = graph.key.key(pos="tr"))
	setTexrunner(gSpec)
	
	cnvs = canvas.canvas()
	cnvs.insert(gCorr)
	cnvs.insert(gAsym)
	cnvs.insert(gSpec)
	
	
	specCircAttrs = None
	bgCircAttrs = [deco.filled]
	if inColor:
		specCircAttrs = [color.rgb.blue]
		bgCircAttrs += [color.rgb(0,0.7,0)]
		
	gSpec.plot(graph.data.points([(p.KE,p.ssfg) for p in pdat],x=1,y=2,title="Data"),
			   [graph.style.symbol(symbol.circle,size=0.1,symbolattrs=specCircAttrs),])
	gSpec.plot(graph.data.points([(p.KE,p.ssMC) for p in pdat],x=1,y=2,title="Geant4 MC"),
			 [graph.style.line(),])
	gSpec.plot(graph.data.points([(p.KE,p.ssbg) for p in pdat],x=1,y=2,title="Background"),
			   [graph.style.symbol(symbol.circle,size=0.1,symbolattrs=bgCircAttrs),])
	
	pdat = [p for p in pdat if p.KE > 60]
	
	area = errorBand([(p.KE,100.*p.corr,100.*p.dcorr) for p in pdat],0,1,2)
	asLineAttrs = [style.linewidth.THick]
	if inColor:
		gCorr.fill(area, [deco.filled([color.rgb(0.5,0.5,1)])])
		asLineAttrs += [color.rgb.green]
	else:
		gCorr.fill(area, [pattern.hatched(0.10,-45)])
		gCorr.fill(area, [pattern.hatched(0.10,45)])
	
	gCorr.plot(graph.data.points([(p.KE,100.*p.corr,100.*p.dcorr) for p in pdat],x=1,y=2,dy=3,title=None),
			 [graph.style.line(asLineAttrs), ]) #[graph.style.symbol(symbol.circle,size=0.1), ]) #graph.style.errorbar()])
	for y in [-2.5,0,2.5]:
		gCorr.plot(graph.data.points([(0,y),(800,y)],x=1,y=2,title=None), [graph.style.line([style.linestyle.dashed])])

	if False:
		d = gCorr.plot(graph.data.points([(p.KE,100.*p.dcorr) for p in pdat],x=1,y=2,title=None),
				 [graph.style.line()])
		
		gCorr.finish()
		p = d.path

		pa = gCorr.xgridpath(100)
		pb = gCorr.xgridpath(750)
		p0 = gCorr.ygridpath(0)
		(splita,), (splitpa,) = p.intersect(pa)
		(splitb,), (splitpb,) = p.intersect(pb)
		splitp0a = p0.intersect(pa)[1][0]
		splitp0b = p0.intersect(pb)[1][0]
		
		area = (pa.split([splitp0a,splitpa])[1] <<
			   p.split([splita, splitb])[1] <<
			   pb.split([splitp0b,splitpb])[1].reversed())
		area[-1].close()
		gCorr.stroke(area, [deco.filled([color.gray(0.8)])])
	
	
	asymBarAttrs = None
	asymCircAttrs = [deco.filled([color.gray(1.0)])]
	asymFitAttrs = [style.linewidth.THIck]
	if inColor:
		asymBarAttrs = [color.rgb.blue]
		asymCircAttrs += [color.rgb.blue]
		asymFitAttrs += [color.rgb(0,0.7,0)]
		
	gAsym.plot(graph.data.points([(p.KE,p.A0,p.dA0) for p in pdat],x=1,y=2,dy=3,title=None),
			   [graph.style.errorbar(errorbarattrs=asymBarAttrs), graph.style.symbol(symbol.circle,size=0.1,symbolattrs=asymCircAttrs)])
	A0 = -0.11954
	gAsym.plot(graph.data.points([(220,A0),(670,A0)],x=1,y=2,title=None), [graph.style.line(asymFitAttrs)])
	
	# different fit ranges
	for (e0,e1) in [(220,670),]+[(200-10*x,700+10*x) for x in range(11)]:
		LF = LinearFitter(terms=[polyterm(0)])
		fdat = [(p.KE,p.A0,p.dA0) for p in pdat if e0 < p.KE < e1]
		LF.fit(fdat,cols=(0,1,2),errorbarWeights=True)
		print "<A0> [%i,%i] ="%(e0,e1),LF.coeffs[0],"chi^2/ndf =",LF.ssResids(),len(fdat)-1

	cnvs.writetofile(os.environ["UCNA_ANA_PLOTS"]+"/Paper/Fig1_new.pdf")

	if False:
		fgtot = bgtot = 0
		for p in pdat:
				if 100 < p.KE and p.KE < 800:
				    print p.KE,p.ssbg,p.ssfg,p.ssfg/p.ssbg
				if True or (220 < p.KE and p.KE < 670):
					fgtot += p.ssfg
					bgtot += p.ssbg
		print "\n\n",fgtot,bgtot,fgtot/bgtot


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
	#G4_Pen_Compare()
	paperFig1(True)
