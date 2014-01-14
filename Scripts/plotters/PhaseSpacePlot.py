#!/usr/bin/python

import sys
sys.path.append("..")

from ucnacore.PyxUtils import *
from ucnacore.LinFitter import *
import os

PDG2010 = {"year":2010, "lambda":1.2694, "d_lambda":0.0028, "tau":885.7, "d_tau":0.8, "Vud":0.97425, "d_Vud":0.00022}
PDG2012 = {"year":2012, "lambda":1.2701, "d_lambda":0.0025, "tau":880.1, "d_tau":1.1, "Vud":0.97425, "d_Vud":0.00022}

# PDG 2012 values
lambda_PDG = 1.2701
d_lambda_PDG = 0.0025
lambda_Mund = 1.2755	 # Mund 2012 eq. 8, combined result
d_lambda_Mund = 0.0013
tau_PDG = 880.1
d_tau_PDG = 1.1
tau_SBV = 878.5
d_tau_SBV = 0.7
Vud_PDG = 0.97425
d_Vud_PDG = 0.00022

# Serebrov 2005 lifetime
tau_sbv = 878.5
d_tau_sbv = 0.8

# Nico 2005 NIST beam lifetime
tau_nico = 886.3
d_tau_nico = 3.4

# Yue 2013 NIST beam lifetime
tau_yue = 887.7
d_tau_yue = 2.25

# Perkeo II 2002 (Abele)
lambda_PII02 = 1.2739
d_lambda_PII02 = .0019

def Vud(gA,tau):
	return sqrt(4908.7/(tau*(1+3*gA**2)))
def dVdTau(gA,tau):
	return 0.5*sqrt(4908.7/(tau**3*(1+3*gA**2)))


class PhaseSpacePlotter:
	def __init__(self):
		self.l0,self.l1 = (1.265,1.285)
		self.V0,self.V1 = (.969,.979)
		self.npts = 100

		self.gPS=graph.graphxy(width=12,height=12,
			   x=graph.axis.lin(title="$|\\lambda| \\equiv |g_A/g_V|$",min=self.l0,max=self.l1,
							parter=graph.axis.parter.linear(tickdists=[0.005,0.001])),
			   y=graph.axis.lin(title="$V_{ud}$",min=self.V0,max=self.V1),
			   key = graph.key.key(pos="bl"))
		setTexrunner(self.gPS)

	def Vud_area(self,Vud,d_Vud):
		gdat = [(l,Vud,d_Vud) for l in unifrange(self.l0,self.l1,self.npts)]
		return errorBand(self.gPS,gdat,0,1,2)

	def lambda_area(self,lm,d_lm,dV0=0):
		gdat = [(lm,v,d_lm) for v in unifrange(self.V0+dV0,self.V1,self.npts)]
		return errorBand(self.gPS,gdat,0,1,2,True)

	def tau_area(self,tau,d_tau):
		gdat = [(l,Vud(l,tau),dVdTau(l,tau)*d_tau) for l in unifrange(self.l0,self.l1,self.npts)]
		return errorBand(self.gPS,gdat,0,1,2)

	def PDG_bands(self,pdg):
		return self.Vud_area(pdg["Vud"],pdg["d_Vud"]),self.lambda_area(pdg["lambda"],pdg["d_lambda"]),self.tau_area(pdg["tau"],pdg["d_tau"])


	def ucnaPRL2012(self, drawMode = "hatched"):
	
		area_Vud,area_PDGgA,area_PDGtau = self.PDG_bands(PDG2012)
		area_MundgA = self.lambda_area(lambda_Mund,d_lambda_Mund,0.002)
		area_UCNA = self.lambda_area(1.2756,0.0030)
		area_YueTau = self.tau_area(tau_yue,d_tau_yue)
		area_NicoTau = self.tau_area(tau_nico,d_tau_nico)
		
		if False:
			gdat = [(l,Vud(l,tau_SBV),dVdTau(l,tau_SBV)*d_tau_SBV) for l in unifrange(l0,l1,npts)]
			area = errorBand(gPS,gdat,0,1,2)
			self.gPS.fill(area, [deco.filled([color.rgb(0.0,1.0,0.0),color.transparency(0.5)])])
		
		if drawMode == "color":
			self.gPS.fill(area_Vud, [deco.filled([color.rgb(0.0,0.0,1.0),color.transparency(0.5)])])
			#gPS.stroke(area_PDGgA, [style.linestyle.solid, color.rgb.blue,
			#				    deco.filled([color.rgb(0.0,1.0,1.0),color.transparency(0.5)])])
			self.gPS.stroke(area_MundgA, [style.linestyle.solid, color.rgb.blue,
								deco.filled([color.rgb(0.0,1.0,1.0),color.transparency(0.5)])])
			self.gPS.stroke(area_UCNA, [style.linestyle.dashed, color.rgb.red,
							   deco.filled([color.rgb(1.0,0.0,0.0),color.transparency(0.5)])])
			self.gPS.fill(area_PDGtau, [deco.filled([color.rgb(0.0,1.0,0.0),color.transparency(0.5)])])

		elif drawMode == "grey":
			#self.gPS.stroke(area_PDGgA, [style.linestyle.solid, color.rgb.blue,
			#				    deco.filled([color.rgb(0.0,1.0,1.0),color.transparency(0.5)])])
			self.gPS.stroke(area_UCNA, [style.linestyle.dashed,
							   deco.filled([color.rgb(0.7,0.7,0.7),color.transparency(0.5)])])
			self.gPS.stroke(area_MundgA, [style.linestyle.solid,
								deco.filled([color.rgb(0.4,0.4,0.4),color.transparency(0.5)])])
			self.gPS.fill(area_PDGtau, [deco.filled([color.rgb(0.3,0.3,0.3),color.transparency(0.2)])])
			self.gPS.fill(area_Vud, [deco.filled([color.rgb(0.5,0.5,0.5),color.transparency(0.2)])])

		elif drawMode == "hatched":
			self.gPS.fill(area_Vud, [pattern.hatched(0.07,0)])
			self.gPS.fill(area_MundgA, [pattern.hatched(0.15,-45)])
			self.gPS.fill(area_UCNA, [pattern.hatched(0.10,45)])
			self.gPS.fill(area_PDGtau, [pattern.hatched(0.07,90)])
			#self.gPS.fill(area_NicoTau, [pattern.hatched(0.07,75)])
			self.gPS.fill(area_YueTau, [pattern.hatched(0.07,75)])

		self.gPS.text(9.0,6.8,"$0^+ \\rightarrow 0^+$")

		tbox = text.text(0,0,"$\\tau_n$ PDG")
		self.gPS.insert(tbox,[trafo.rotate(-52),trafo.translate(9.0,3.55)])
		#gPS.text(9.5,2.5,"$\\tau_n$ PDG")

		tbox = text.text(0,0,"$\\tau_n$ NIST")
		self.gPS.insert(tbox,[trafo.rotate(-52),trafo.translate(2.2,3.55)])
		
		tbox = text.text(5.55, 1, "UCNA")
		tpath = tbox.bbox().enlarged(0.2).path()
		self.gPS.stroke(tpath,[deco.filled([color.rgb.white])])
		self.gPS.insert(tbox)

		tbox = text.text(0, 0, "PERKEO II '12")
		tpath = tbox.bbox().enlarged(0.2).path()
		c = canvas.canvas()
		c.stroke(tpath,[deco.filled([color.rgb.white])])
		c.insert(tbox)
		self.gPS.insert(c,[trafo.rotate(90),trafo.scale(0.9),trafo.translate(6.50,8.0)])
		
		if False:
			tbox = text.text(1.95, 1, r"PDG '12")
			tpath = tbox.bbox().enlarged(0.2).path()
			self.gPS.stroke(tpath,[deco.filled([color.rgb.white])])
			self.gPS.insert(tbox)

		self.gPS.writetofile(os.environ["UCNA_ANA_PLOTS"]+"/Paper/PhaseSpace2012_%s.pdf"%drawMode)


	def thesis_conflicted(self):

		a_Vud,a_PDGgA,a_PDGtau = self.PDG_bands(PDG2010)
		a_SBVtau = self.tau_area(tau_sbv,d_tau_sbv)
		a_PII02 = self.lambda_area(lambda_PII02,d_lambda_PII02,0.00)
		
		self.gPS.fill(a_Vud, [pattern.hatched(0.07,0)])
		self.gPS.fill(a_PDGgA, [pattern.hatched(0.15,-45)])
		self.gPS.fill(a_PDGtau, [pattern.hatched(0.07,90)])
		self.gPS.fill(a_SBVtau, [pattern.hatched(0.07,90)])
		self.gPS.fill(a_PII02, [pattern.hatched(0.10,45)])
		
		self.gPS.text(9.0,6.8,"$0^+ \\rightarrow 0^+$")

		tbox = text.text(0,0,"$\\tau_n$ PDG")
		self.gPS.insert(tbox,[trafo.rotate(-52.5),trafo.translate(7.0,2.25)])

		tbox = text.text(0,0,"$\\tau_n$ Serebrov")
		self.gPS.insert(tbox,[trafo.rotate(-52.5),trafo.translate(9.2,4.25)])
			
		tbox = text.text(0, 0, "PERKEO II '02")
		tpath = tbox.bbox().enlarged(0.2).path()
		c = canvas.canvas()
		c.stroke(tpath,[deco.filled([color.rgb.white])])
		c.insert(tbox)
		self.gPS.insert(c,[trafo.rotate(90),trafo.translate(5.55,7.5)])
		
		tbox = text.text(0, 0, "PDG 2010")
		tpath = tbox.bbox().enlarged(0.2).path()
		c = canvas.canvas()
		c.stroke(tpath,[deco.filled([color.rgb.white])])
		c.insert(tbox)
		self.gPS.insert(c,[trafo.rotate(90),trafo.translate(2.75,1.5)])

		self.gPS.writetofile(os.environ["UCNA_ANA_PLOTS"]+"/Paper/PhaseSpace2010.pdf")


def A_meas_history():

	gdat = [	(1975,-0.113,0.006),					# Krohn 75 Argonne
				(1986,-0.1146,0.0019),					# Bopp 86
				(1995,-0.1160,sqrt(.0009**2+.0012**2)),	# TPC Liaud '95
				(1997,-0.1189,0.0012),					# Abele '97
				(2002,-0.1189,0.0007),					# Abele '02
				#(2012.3,-0.11972,0.00065)				# Mund '13
			]

	gYeroz = [
				(1990.5,-0.1116,0.0014),				# Yeroz. 1990, 1991
				(1996,-0.1135,0.0014),					# Yeroz. '97 - revised
			]
				
	gUCNA = [	(2009,-0.1138,sqrt(.0046**2+.0021**2)),
				(2010,-.11966,sqrt(.00089**2+.0014**2)),
				#(2013.1,-0.11954,sqrt(.00055**2+.00098**2))
			]
	
	gA = graph.graphxy(width=15,height=10,
			   x=graph.axis.lin(title="Year",max=2015),
			   y=graph.axis.lin(title="$A_0$"),
			   key = graph.key.key(pos="bl"))
	setTexrunner(gA)
	
	gA.plot(graph.data.points(gdat,x=1,y=2, dy=3, title="CN experiments"),
			   [graph.style.errorbar(),graph.style.symbol(symbol.circle,size=0.22,symbolattrs=[deco.filled()]),])
	gA.plot(graph.data.points(gYeroz,x=1,y=2, dy=3, title=None),
			   [graph.style.errorbar(),graph.style.symbol(symbol.circle,size=0.22,symbolattrs=[deco.filled()]),graph.style.line()])
	gA.plot(graph.data.points(gUCNA,x=1,y=2, dy=3, title="UCNA"),
			   [graph.style.errorbar(),graph.style.symbol(symbol.diamond,size=0.22,symbolattrs=[deco.filled([rgb.white])]),])
			   
	gA.writetofile(os.environ["UCNA_ANA_PLOTS"]+"/Paper/AsymHistory_2010.pdf")

if __name__=="__main__":
	
	#A_meas_history()
	
	PSP = PhaseSpacePlotter()
	PSP.ucnaPRL2012()
	
	#PSP = PhaseSpacePlotter()
	#PSP.thesis_conflicted()