#!/usr/bin/python

from PyxUtils import *
from QFile import *
from math import *
import os

class spectrumCorrs(KVMap):
	def __init__(self,m):
		KVMap.__init__(self,m)
		self.loadFloats(self.dat.keys())
		
def product1(l):
	p = 1.
	for x in l:
		p *= 1.+x
	return p

def tensorCoeff(E):
	m_e = 511.
	l = 1.2694
	return 2/(1+3*l)**2/938000*(m_e**2./E*l**2*(1-l)+782*(l**2+2./3.*l-1./3.)+E*2./3.*(1+l+3*l**2+3*l**3))*(1+3.*l**2)/(2*l*(1-l))

def gluckx2E(x):
	return x*(782.6+511.)-511.
	
	
def plotWilkinsonCorrs(fin,outpath):

	os.system("mkdir -p %s"%outpath)
	os.system("cd ..; ./UCNAnalyzer rc x");
	
	corrs = [spectrumCorrs(m) for m in QFile(fin).dat["spectrumPoint"]]
		
	##############
	# individual spectrum shape corrections
	##############
	
	gWC=graph.graphxy(width=24,height=10,
			x=graph.axis.lin(title="Energy [keV]",min=0,max=1000),
			y=graph.axis.log(title="Spectrum Correction",min=5e-6,max=0.5),
			key = graph.key.key(pos="tr"))
	setTexrunner(gWC)

	gdat = [ [	c.energy,		# 0
				c.F0m1,			# 1
				c.L0m1,			# 2
				-c.Qm1,			# 3
				c.RVm1,			# 4
				c.RAm1,			# 5
				abs(c.BiRWM),	# 6
				-c.VCm1,		# 7
				-c.ACm1,		# 8
				c.g,			# 9
				0.000401,		# 10
				0.0004*(1+0.2*sin(c.energy/30)) ] for c in corrs if c.energy<783]
				
	gWC.plot(graph.data.points(gdat,x=1,y=2,title="$F_0-1$"),[graph.style.line([style.linewidth.Thick])])
	gWC.plot(graph.data.points(gdat,x=1,y=10,title="g"),[graph.style.line([rgb.green,style.linestyle.dashdotted,style.linewidth.Thick])])
	gWC.plot(graph.data.points(gdat,x=1,y=6,title="$R^A-1$"),[graph.style.line([rgb.red,style.linestyle.dotted,style.linewidth.Thick])])
	gWC.plot(graph.data.points(gdat,x=1,y=5,title="$R^V-1$"),[graph.style.line([rgb.red,style.linestyle.dashed,style.linewidth.Thick])])
	gWC.plot(graph.data.points(gdat,x=1,y=7,title="$|R+{\\rm WM}|$"),[graph.style.line([rgb.red,style.linewidth.Thick])])
	
	
	gWC.plot(graph.data.points(gdat,x=1,y=11,title="$J(Z)-1$"),[graph.style.line([rgb(0.7,0.,0.7),style.linewidth.Thick])])
	gWC.plot(graph.data.points(gdat,x=1,y=12,title="$O[ K(\\alpha) ]$"),[graph.style.line([rgb(0.7,0.,0.7),style.linestyle.dotted,style.linewidth.Thick])])
	
	gWC.plot(graph.data.points(gdat,x=1,y=4,title="$L_0-1$"),[graph.style.line([style.linestyle.dashed,style.linewidth.Thick])])
	gWC.plot(graph.data.points(gdat,x=1,y=8,title="$1-{^VC}$"),[graph.style.line([rgb.blue,style.linestyle.dashed,style.linewidth.Thick])])
	gWC.plot(graph.data.points(gdat,x=1,y=9,title="$1-{^AC}$"),[graph.style.line([rgb.blue,style.linestyle.dotted,style.linewidth.Thick])])
	gWC.plot(graph.data.points(gdat,x=1,y=4,title="$1-Q$"),[graph.style.line([style.linestyle.dotted,style.linewidth.Thick])])

	gWC.writetofile(outpath+"/WilkinsonCorrs.pdf")
	
	##############
	# total spectrum shape corrections
	##############
	
	gWCTot=graph.graphxy(width=20,height=10,
			x=graph.axis.lin(title="Energy [keV]",min=0,max=783),
			y=graph.axis.lin(title="Total spectrum correction factor",min=1.0,max=1.2),
			key = graph.key.key(pos="tr"))
	setTexrunner(gWCTot)
	gdat = [[c.energy,c.S/c.S0] for c in corrs if c.S0]
	gdat.sort()
	gWCTot.plot(graph.data.points(gdat,x=1,y=2,title=None),[graph.style.line([style.linewidth.Thick])])
	gWCTot.writetofile(outpath+"/WilkinsonCorrsTot.pdf")
	
	###############
	# spectrum shapes
	###############

	gSpec=graph.graphxy(width=10,height=8,
			x=graph.axis.lin(title="Energy [keV]",min=0,max=800),
			y=graph.axis.lin(title="Decay Rate",min=0),
			key = graph.key.key(pos="tr"))
	setTexrunner(gSpec)
	gdat = [[c.energy,c.S0,c.S] for c in corrs if c.energy<783]
	gdat.sort()
	n = len(gdat)/2
	snorm = gdat[n][1]/gdat[n][2]
	gdat = [g+[g[2]*snorm,] for g in gdat]
	gSpec.plot(graph.data.points(gdat,x=1,y=3,title="Corrected"),[graph.style.line([style.linestyle.dashed,rgb.blue,style.linewidth.Thick])])
	gSpec.plot(graph.data.points(gdat,x=1,y=4,title="Normalized"),[graph.style.line([rgb.red,style.linewidth.Thick])])
	gSpec.plot(graph.data.points(gdat,x=1,y=2,title="Plain"),[graph.style.line([style.linestyle.dotted,style.linewidth.Thick])])
	gSpec.writetofile(outpath+"/BetaSpectrum.pdf")

	################
	# Asymmetry corrections
	################
	
	gWCA=graph.graphxy(width=20,height=10,
			x=graph.axis.lin(title="Energy [keV]",min=0,max=800),
			y=graph.axis.lin(title="Fractional Correction to $A$",min=0.0000),
			key = graph.key.key(pos="tl"))
	setTexrunner(gWCA)
	gWCA.plot(graph.data.points([(c.energy,c.RWM) for c in corrs],x=1,y=2,title="Wilkinson recoil-order"),[graph.style.line([style.linewidth.Thick])])
	gWCA.plot(graph.data.points([(c.energy,c.hmg) for c in corrs],x=1,y=2,title="Radiative $h - g$"),[graph.style.line([style.linewidth.Thick,rgb.red])])
	if 0:
		gluck_da = [ [gluckx2E(0.4),0.03/1.66],[gluckx2E(0.5),0.01/6.55],[gluckx2E(0.6),0.01/8.08],
						[gluckx2E(0.7),0.01/8.90],[gluckx2E(0.8),0.01/9.41],[gluckx2E(0.9),0.01/9.76]]
		#gluck_da = [ [gluckx2E(0.4),0.003],[gluckx2E(0.5),0.001],[gluckx2E(0.6),0.001],
		#				[gluckx2E(0.7),0.001],[gluckx2E(0.8),0.001],[gluckx2E(0.9),0.001]]
		gWCA.plot(graph.data.points(gluck_da,x=1,y=2,title='Gl\\"uck 1992 $\\delta \\alpha_e$ table'),[graph.style.symbol(symbolattrs=[rgb.blue])])
	gWCA.writetofile(outpath+"/A_Rad_Corrections.pdf")

	################
	# Asymmetry curves
	################
	
	gAsym=graph.graphxy(width=10,height=8,
			x=graph.axis.lin(title="Energy [keV]",min=0,max=800),
			y=graph.axis.lin(title="Asymmetry $A(E)\\beta/2$"),
			key = graph.key.key(pos="tr"))
	setTexrunner(gAsym)
	gdat = [[c.energy,c.A0,c.A] for c in corrs if c.energy<783]
	gAsym.plot(graph.data.points(gdat,x=1,y=2,title="$A(E)=A_0$"),[graph.style.line([style.linewidth.Thick])])
	gAsym.plot(graph.data.points(gdat,x=1,y=3,title="Corrected"),[graph.style.line([rgb.red,style.linestyle.dashdotted,style.linewidth.Thick])])
	gAsym.writetofile(outpath+"/Asymmetry.pdf")

	
	
	if 0:
		# size of tensor corrections
		gT=graph.graphxy(width=20,height=10,
				x=graph.axis.lin(title="Energy [keV]",min=0,max=800),
				y=graph.axis.lin(title="$\\Delta A/A$ for $g_{II}/g_V=1$",min=-0.01,max=0.01),
				key = graph.key.key(pos="tr"))
		setTexrunner(gT)
		gT.plot(graph.data.points(gdat,x=1,y=12,title=None),[graph.style.line([style.linewidth.Thick])])
		gT.writetofile(outpath+"/GardnerTensor.pdf")
		
		
		
	
if __name__ == "__main__":
	plotWilkinsonCorrs(os.environ["UCNA_ANA_PLOTS"]+"/SpectrumCorrection/SpectrumCorrection.txt",os.environ["UCNA_ANA_PLOTS"]+"/SpectrumCorrection/")
	