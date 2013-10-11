#!/usr/bin/python

from ucnacore.PyxUtils import *
from ucnacore.QFile import *
from math import *
import os
import sys
sys.path.append("Corrections")
from EnergyErrors import *

class spectrumCorrs(KVMap):
	def __init__(self,m):
		KVMap.__init__(self,m)
		self.loadFloats(self.dat.keys())
		
		
def plotWilkinsonCorrs(fin,outpath):

	#os.system("mkdir -p %s"%outpath)
	#os.system("cd ..; ./UCNAnalyzer rc x");
	
	Q = QFile(fin)
	A = Q.getItemF("decayInfo","A")
	Z = Q.getItemF("decayInfo","Z")
	ep = Q.getItemF("decayInfo","endpt")
	corrs = [spectrumCorrs(m) for m in Q.dat["spectrumPoint"]]
		
	##############
	# individual spectrum shape corrections
	##############
	
	gWC=graph.graphxy(width=18,height=12,
			x=graph.axis.lin(title="Energy [keV]",min=0,max=ep), #max=ep*1.25),
			#y=graph.axis.log(title="Spectrum Correction",min=5e-6,max=0.5),
			y=graph.axis.log(title="Spectrum Correction",min=5e-6),
			key = graph.key.key(pos="tc",columns=4))
	setTexrunner(gWC)

	gdat = [ [	c.energy,			# 0
				c.F0m1,			# 1
				abs(c.L0m1),		# 2
				-c.Qm1,			# 3
				c.RVm1,			# 4
				c.RAm1,			# 5
				abs(c.BiRWM),		# 6
				-c.VCm1,			# 7
				-c.ACm1,			# 8
				c.g,				# 9
				0.000401,			# 10
				0.0004*(1+0.2*sin(c.energy/30)), # 11
				-c.Cm1,			# 12
				c.Rm1			# 13
			] for c in corrs]
			
	gWC.plot(graph.data.points([(c.energy,c.F0m1) for c in corrs],x=1,y=2,title="$F_0-1$"),[graph.style.line([style.linewidth.THICk])])
	gWC.plot(graph.data.points([(c.energy,abs(c.g)) for c in corrs],x=1,y=2,title="$|g \\alpha/2\\pi|$"),[graph.style.line([style.linestyle.dashed,style.linewidth.THick])])
	gWC.plot(graph.data.points([(c.energy,abs(c.BiRWM)) for c in corrs],x=1,y=2,title="$|R+{\\rm WM}|$"),[graph.style.line([style.linestyle.dotted,style.linewidth.THIck])])
	#gWC.plot(graph.data.points([(c.energy,-c.ACm1) for c in corrs],x=1,y=2,title="$1-{^AC}$"),[graph.style.line([style.linewidth.thick])])
	#gWC.plot(graph.data.points([(c.energy,-c.VCm1) for c in corrs],x=1,y=2,title="$1-{^VC}$"),[graph.style.line([style.linestyle.dashdotted,style.linewidth.thick])])
	gWC.plot(graph.data.points([(c.energy,-c.Cm1) for c in corrs],x=1,y=2,title="$1-C$"),[graph.style.line([style.linestyle.dashdotted,style.linewidth.thick])])
	gWC.plot(graph.data.points([(c.energy,-c.L0m1) for c in corrs],x=1,y=2,title="$1-L_0$"),[graph.style.line([style.linestyle.dashed,style.linewidth.thick])])
	gWC.plot(graph.data.points([(c.energy,-c.Qm1) for c in corrs],x=1,y=2,title="$1-Q$"),[graph.style.line([style.linestyle.dotted,style.linewidth.thick])])
	

	
	#gWC.plot(graph.data.points(gdat,x=1,y=10,title="$g$"),[graph.style.line([rgb.green,style.linestyle.dashdotted,style.linewidth.THick])])
	#gWC.plot(graph.data.points(gdat,x=1,y=6,title="$R^A-1$"),[graph.style.line([rgb.red,style.linestyle.dotted,style.linewidth.THick])])
	#gWC.plot(graph.data.points(gdat,x=1,y=5,title="$R^V-1$"),[graph.style.line([rgb.red,style.linestyle.dashed,style.linewidth.THick])])
	#gWC.plot(graph.data.points(gdat,x=1,y=14,title="$R-1$"),[graph.style.line([rgb.red,style.linestyle.dashdotted,style.linewidth.thick])])
	#gWC.plot(graph.data.points(gdat,x=1,y=7,title="$|R+{\\rm WM}|$"),[graph.style.line([rgb.red,style.linewidth.THick])])
	
	#gWC.plot(graph.data.points(gdat,x=1,y=11,title="$J(Z)-1$"),[graph.style.line([rgb(0.7,0.,0.7),style.linestyle.dashdotted,style.linewidth.thick])])
	#gWC.plot(graph.data.points(gdat,x=1,y=12,title="$O[ K(\\alpha) ]$"),[graph.style.line([rgb(0.7,0.,0.7),style.linestyle.dotted,style.linewidth.thick])])
	
	#gWC.plot(graph.data.points(gdat,x=1,y=4,title="$1-Q$"),[graph.style.line([style.linestyle.dotted,style.linewidth.Thick])])
	#gWC.plot(graph.data.points(gdat,x=1,y=8,title="$1-{^VC}$"),[graph.style.line([rgb.blue,style.linestyle.dashed,style.linewidth.Thick])])
	#gWC.plot(graph.data.points(gdat,x=1,y=9,title="$1-{^AC}$"),[graph.style.line([rgb.blue,style.linestyle.dotted,style.linewidth.Thick])])
	#gWC.plot(graph.data.points(gdat,x=1,y=13,title="$1-C$"),[graph.style.line([rgb.blue,style.linestyle.dashdotted,style.linewidth.thick])])
	#gWC.plot(graph.data.points(gdat,x=1,y=3,title="$L_0-1$"),[graph.style.line([style.linestyle.dashed,style.linewidth.Thick])])
	
	gWC.writetofile(outpath+"/WilkinsonCorrs_%i_%i_%f.pdf"%(A,Z,ep))
	return
	
	###############
	# spectrum shapes
	###############
	
	gSpec=graph.graphxy(width=10,height=8,
						x=graph.axis.lin(title="Energy [keV]",min=0,max=ep),
						y=graph.axis.lin(title="Decay Rate",min=0),
						key = graph.key.key(pos="tr"))
	setTexrunner(gSpec)
	gdat = [[c.energy,c.S0,c.S] for c in corrs]
	gdat.sort()
	print "Uncorrected beta energy:",sum([g[0]*g[1] for g in gdat])/sum([g[1] for g in gdat])
	print "Average beta energy:",sum([g[0]*g[2] for g in gdat])/sum([g[2] for g in gdat])
	snorm = sum([g[1] for g in gdat])/sum([g[2] for g in gdat])
	gdat = [g+[g[2]*snorm,] for g in gdat]
	if 1/1.2 < snorm < 1.2:
		gSpec.plot(graph.data.points(gdat,x=1,y=3,title="Corrected"),[graph.style.line([style.linestyle.dashed,rgb.blue,style.linewidth.Thick])])
		gSpec.plot(graph.data.points(gdat,x=1,y=4,title="Normalized"),[graph.style.line([rgb.red,style.linewidth.Thick])])
	else:
		gSpec.plot(graph.data.points(gdat,x=1,y=4,title="Corrected"),[graph.style.line([style.linestyle.dashed,rgb.blue,style.linewidth.Thick])])
	gSpec.plot(graph.data.points(gdat,x=1,y=2,title="Plain"),[graph.style.line([style.linestyle.dotted,style.linewidth.Thick])])
	gSpec.writetofile(outpath+"/BetaSpectrum_%i_%i_%f.pdf"%(A,Z,ep))

	if A != 1:
		return
	
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
	
	################
	# Asymmetry corrections
	################
	
	gWCA=graph.graphxy(width=20,height=10,
			x=graph.axis.lin(title="Energy [keV]",min=0,max=800),
			y=graph.axis.lin(title="Fractional Correction to $A$",min=0.0000),
			key = graph.key.key(pos="tl"))
	setTexrunner(gWCA)
	gWCA.plot(graph.data.points([(c.energy,c.RWM) for c in corrs],x=1,y=2,title="recoil-order (Wilkinson)"),[graph.style.line([style.linewidth.Thick])])
	gWCA.plot(graph.data.points([(c.energy,c.hmg) for c in corrs],x=1,y=2,title="Radiative $h - g$ (Shann/Sirlin)"),[graph.style.line([style.linewidth.Thick,rgb.red])])
	if 1:
		def gluckx2E(x):
			return x*(782.6+511.)-511.
		gluck_da = [ [gluckx2E(0.4),0.03/1.66],[gluckx2E(0.5),0.01/6.55],[gluckx2E(0.6),0.01/8.08],
						[gluckx2E(0.7),0.01/8.90],[gluckx2E(0.8),0.01/9.41],[gluckx2E(0.9),0.01/9.76]]

		gluck_hmg = [ [gluckx2E(0.4),0.005/1.731],[gluckx2E(0.5),0.014/6.961],[gluckx2E(0.6),0.013/8.553],
				   [gluckx2E(0.7),0.011/9.383],[gluckx2E(0.8),0.009/9.884],[gluckx2E(0.9),0.008/10.214]]

		#gWCA.plot(graph.data.points(gluck_da,x=1,y=2,title='Gl\\"uck 1992 $\\delta \\alpha_e$ table'),[graph.style.symbol(symbolattrs=[rgb.green])])
		gWCA.plot(graph.data.points(gluck_hmg,x=1,y=2,title='Gl\\"uck private communication'),[graph.style.symbol(symbolattrs=[rgb.green])])
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

	
	
	################
	# shape of tensor coupling corrections to A
	################
	
	def tensorCoeff(E):
		m_e = 511.
		l = 1.2694
		return 2/(1+3*l)**2/938000*(m_e**2./E*l**2*(1-l)+782*(l**2+2./3.*l-1./3.)+E*2./3.*(1+l+3*l**2+3*l**3))*(1+3.*l**2)/(2*l*(1-l))
	gT=graph.graphxy(width=20,height=10,
			x=graph.axis.lin(title="Energy [keV]",min=0,max=800),
			y=graph.axis.lin(title="$\\Delta A/A$ for $g_{II}/g_V=1$",min=-0.01,max=0.01),
			key = graph.key.key(pos="tr"))
	setTexrunner(gT)
	gT.plot(graph.data.points([(c.energy,tensorCoeff(c.energy)) for c in corrs],x=1,y=2,title=None),[graph.style.line([style.linewidth.Thick])])
	gT.writetofile(outpath+"/GardnerTensor.pdf")


def plot_hg_radiative(fin,outpath):
	"""h-g radiative corrections for comparison with Albert's numbers"""
	Q = QFile(fin)
	ep = Q.getItemF("decayInfo","endpt")
	corrs = [spectrumCorrs(m) for m in Q.dat["spectrumPoint"]]
	ghg=graph.graphxy(width=24,height=16,
			   x=graph.axis.lin(title="Energy [keV]",min=500,max=1300),
			   y=graph.axis.lin(title="Relative Asymmetry Correction",min=0,max=0.003),
			   key = graph.key.key(pos="mr"))
	setTexrunner(ghg)

	ghg.plot(graph.data.points([(c.energy+511,c.hmg) for c in corrs],x=1,y=2,title="Radiative $h - g$ (Shann/Sirlin)"),
		    [graph.style.line([style.linewidth.Thick,rgb.blue])])

	ghg.writetofile(outpath+"/Radiative_h-g.pdf")



def thesisCorrsPlot(fin,outpath):
	
	Q = QFile(fin)
	A = Q.getItemF("decayInfo","A")
	Z = Q.getItemF("decayInfo","Z")
	ep = Q.getItemF("decayInfo","endpt")
	corrs = [spectrumCorrs(m) for m in Q.dat["spectrumPoint"]]
		
	##############
	# Fine spectrum shape corrections
	##############
	
	gWC=graph.graphxy(width=16,height=8,
			x=graph.axis.lin(title="Energy [keV]",min=0,max=ep),
			y=graph.axis.lin(title="Spectrum Correction"),
			key = graph.key.key(pos="bc",columns=3))
	setTexrunner(gWC)
				
	gWC.plot(graph.data.points([(c.energy,c.g) for c in corrs],x=1,y=2,title="$g \\alpha/2\\pi$"),[graph.style.line([style.linestyle.dashdotted,style.linewidth.THick])])

	gWC.plot(graph.data.points([(c.energy,c.ACm1+c.L0m1) for c in corrs],x=1,y=2,title="${^AC}L_0-1$"),[graph.style.line([style.linewidth.Thick])])
	gWC.plot(graph.data.points([(c.energy,c.VCm1+c.L0m1) for c in corrs],x=1,y=2,title="${^VC} L_0-1$"),[graph.style.line([style.linestyle.dotted,style.linewidth.Thick])])
	
	gWC.writetofile(outpath+"/ThesisSpecCorr_%i_%i_%f.pdf"%(A,Z,ep))

	###############
	# spectrum shapes
	###############
	
	gSpec=graph.graphxy(width=10,height=8,
						x=graph.axis.lin(title="Energy [keV]",min=0,max=ep),
						y=graph.axis.lin(title="Decay Rate",min=0,max=1.1),
						key = graph.key.key(pos="bc"))
	setTexrunner(gSpec)

	gdat = [[c.energy,c.S0,c.S0*(1+c.F0m1)*(1+c.L0m1)*(1+c.VCm1)*(1+c.g)] for c in corrs]
	gdat.sort()
	print "Uncorrected beta energy:",sum([g[0]*g[1] for g in gdat])/sum([g[1] for g in gdat])
	print "Average beta energy:",sum([g[0]*g[2] for g in gdat])/sum([g[2] for g in gdat])
	vnorm = 1./max([g[1] for g in gdat])
	snorm = sum([g[1] for g in gdat])/sum([g[2] for g in gdat])*vnorm
	gdat = [ (g[0],g[1]*vnorm,g[2]*snorm) for g in gdat]
	
	gSpec.plot(graph.data.points(gdat,x=1,y=2,title="plain phase space"),[graph.style.line([style.linestyle.dotted,style.linewidth.Thick])])
	gSpec.plot(graph.data.points(gdat,x=1,y=3,title="corrected"),[graph.style.line([style.linewidth.Thick])])
	gSpec.writetofile(outpath+"/Thesis_Spectrum_%i_%i_%f.pdf"%(A,Z,ep))


if __name__ == "__main__":						  
	plotWilkinsonCorrs(os.environ["UCNA_ANA_PLOTS"]+"/SpectrumCorrection/SpectrumCorrection_1_1_782.347.txt",
					   os.environ["UCNA_ANA_PLOTS"]+"/SpectrumCorrection/")
	#plotWilkinsonCorrs(os.environ["UCNA_ANA_PLOTS"]+"/SpectrumCorrection/SpectrumCorrection_137_55_513.97.txt",
	#				   os.environ["UCNA_ANA_PLOTS"]+"/SpectrumCorrection/")

	#plot_hg_radiative(os.environ["UCNA_ANA_PLOTS"]+"/SpectrumCorrection/SpectrumCorrection_1_1_782.347.txt",
	#		   				   os.environ["UCNA_ANA_PLOTS"]+"/SpectrumCorrection/")

	#thesisCorrsPlot(os.environ["UCNA_ANA_PLOTS"]+"/SpectrumCorrection/SpectrumCorrection_135_55_915.txt",
	#				   os.environ["UCNA_ANA_PLOTS"]+"/SpectrumCorrection/")
