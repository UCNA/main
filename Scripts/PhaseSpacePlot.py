#!/usr/bin/python

from PyxUtils import *
from LinFitter import *
import os

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

def Vud(gA,tau):
	return sqrt(4908.7/(tau*(1+3*gA**2)))
def dVdTau(gA,tau):
	return 0.5*sqrt(4908.7/(tau**3*(1+3*gA**2)))

def lambdaVudPlot(inColor = True):
	
	V0 = Vud_PDG-0.005
	V1 = Vud_PDG+0.005
	l0 = lambda_PDG-0.005
	l1 = lambda_PDG+0.015
	npts = 100
	gPS=graph.graphxy(width=12,height=12,
			   x=graph.axis.lin(title="$|\\lambda| \\equiv |g_A/g_V|$",min=l0,max=l1),
			   y=graph.axis.lin(title="$V_{ud}$",min=V0,max=V1),
			   key = graph.key.key(pos="bl"))
	setTexrunner(gPS)

	# PDG Vud
	gdat = [(l,Vud_PDG,d_Vud_PDG) for l in unifrange(l0,l1,npts)]
	area_Vud = errorBand(gPS,gdat,0,1,2)
	
	# PDG gA
	gdat = [(lambda_PDG,v,d_lambda_PDG) for v in unifrange(V0,V1,npts)]
	area_PDGgA = errorBand(gPS,gdat,0,1,2,True)
		
	# PerkeoII 2012 gA
	gdat = [(lambda_Mund,v,d_lambda_Mund) for v in unifrange(V0+0.002,V1,npts)]
	area_MundgA = errorBand(gPS,gdat,0,1,2,True)
	
	# UCNA
	gdat = [(1.2756,v,0.0030) for v in unifrange(V0,V1,25)]
	area_UCNA = errorBand(gPS,gdat,0,1,2,True)
		
	# PDG lifetime
	gdat = [(l,Vud(l,tau_PDG),dVdTau(l,tau_PDG)*d_tau_PDG) for l in unifrange(l0,l1,npts)]
	area_PDGtau = errorBand(gPS,gdat,0,1,2)
	
	if False:
		gdat = [(l,Vud(l,tau_SBV),dVdTau(l,tau_SBV)*d_tau_SBV) for l in unifrange(l0,l1,npts)]
		area = errorBand(gPS,gdat,0,1,2)
		gPS.fill(area, [deco.filled([color.rgb(0.0,1.0,0.0),color.transparency(0.5)])])
	
	if inColor:
		gPS.fill(area_Vud, [deco.filled([color.rgb(0.0,0.0,1.0),color.transparency(0.5)])])
		#gPS.stroke(area_PDGgA, [style.linestyle.solid, color.rgb.blue,
		#				    deco.filled([color.rgb(0.0,1.0,1.0),color.transparency(0.5)])])
		gPS.stroke(area_MundgA, [style.linestyle.solid, color.rgb.blue,
						    deco.filled([color.rgb(0.0,1.0,1.0),color.transparency(0.5)])])
		gPS.stroke(area_UCNA, [style.linestyle.dashed, color.rgb.red,
						   deco.filled([color.rgb(1.0,0.0,0.0),color.transparency(0.5)])])
		gPS.fill(area_PDGtau, [deco.filled([color.rgb(0.0,1.0,0.0),color.transparency(0.5)])])

	else:
		#gPS.stroke(area_PDGgA, [style.linestyle.solid, color.rgb.blue,
		#				    deco.filled([color.rgb(0.0,1.0,1.0),color.transparency(0.5)])])
		gPS.stroke(area_UCNA, [style.linestyle.dashed,
						   deco.filled([color.rgb(0.7,0.7,0.7),color.transparency(0.5)])])
		gPS.stroke(area_MundgA, [style.linestyle.solid,
							deco.filled([color.rgb(0.4,0.4,0.4),color.transparency(0.5)])])
		gPS.fill(area_PDGtau, [deco.filled([color.rgb(0.3,0.3,0.3),color.transparency(0.2)])])
		gPS.fill(area_Vud, [deco.filled([color.rgb(0.5,0.5,0.5),color.transparency(0.2)])])

		if 0:
			gPS.fill(area_Vud, [pattern.hatched(0.07,0)])
			gPS.fill(area_PDGgA, [pattern.hatched(0.15,-45)])
			gPS.fill(area_UCNA, [pattern.hatched(0.10,45)])
			gPS.fill(area_PDGtau, [pattern.hatched(0.07,90)])

	gPS.text(9.0,6.5,"$0^+ \\rightarrow 0^+$")

	tbox = text.text(0,0,"$\\tau_n$ PDG")
	gPS.insert(tbox,[trafo.rotate(-52),trafo.translate(9.0,3.5)])
	#gPS.text(9.5,2.5,"$\\tau_n$ PDG")

	tbox = text.text(5.5, 1, "UCNA")
	tpath = tbox.bbox().enlarged(0.2).path()
	gPS.stroke(tpath,[deco.filled([color.rgb.white])])
	gPS.insert(tbox)

	tbox = text.text(0, 0, "PERKEO II")
	tpath = tbox.bbox().enlarged(0.2).path()
	c = canvas.canvas()
	c.stroke(tpath,[deco.filled([color.rgb.white])])
	c.insert(tbox)
	gPS.insert(c,[trafo.rotate(90),trafo.translate(6.45,8.5)])
	
	if False:
		tbox = text.text(1.9, 1, r"PDG '12")
		tpath = tbox.bbox().enlarged(0.2).path()
		gPS.stroke(tpath,[deco.filled([color.rgb.white])])
		gPS.insert(tbox)

	if inColor:
		gPS.writetofile(os.environ["UCNA_ANA_PLOTS"]+"/Paper/PhaseSpace.pdf")
	else:
		gPS.writetofile(os.environ["UCNA_ANA_PLOTS"]+"/Paper/Fig2.pdf")

if __name__=="__main__":
	lambdaVudPlot(False)