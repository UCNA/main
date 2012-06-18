#!/usr/bin/python

import sys
sys.path.append("..")
from math import *
from PyxUtils import *
from ErrorEnvelope import *

def bin_edges(w=10,n=100):
	return [w*i for i in range(n+1)]

baseOutPath = "../../Aux/Corrections/"

def writeUncertaintyTable(fout,dat):
	"""Write uncertainties table to file in correct format"""
	fout.write("#E_lo\tE_hi\tcorrection\tuncertainty\n")
	for d in dat:
		fout.write("%i\t%i\t%f\t%f\n"%tuple(d))

def neutronBeta(KE, m = 511.):
	return sqrt(KE**2+2*m*KE)/(m+KE)

def ObsAsymApprox(KE):
	"""Phenomenological fit to MC observed asymmetry"""
	A0 = .1172
	p2 = .966648
	p3 = .0001174
	p4 = 56.5529
	p5 = 12.8861
	p6 = 2.18673
	return A0*neutronBeta(KE)*0.5*p2*(1+p3*KE)*(1+p6/(1+exp((KE-p4)/p5)))

def linearityUncertaintyTable():
	"""Uncertainty due to energy calibration errors; see EnergyErrorsRevis.pdf"""
	edges = bin_edges()
	dat = []
	errmax = 0
	for i in range(len(edges)-1)[-1::-1]:
		c = 0.5*(edges[i]+edges[i+1])
		Eprim = c+calEnvelope(c)
		err = ObsAsymApprox(Eprim)/ObsAsymApprox(c)-1
		if err > errmax:
			errmax = err
		dat.append((edges[i],edges[i+1],0.0,errmax))
	dat = dat[-1::-1]
	fout = open(baseOutPath+"/EnergyLinearityUncertainty_2010.txt","w")
	fout.write("# Uncertainty from energy calibration linearity envelope for 2010 data\n")
	writeUncertaintyTable(fout,dat)
	
def gainFluctsUncertaintyTable():
	"""Uncertainty due to gain fluctuations; see"""
	edges = bin_edges()
	dat = []
	for i in range(len(edges)-1):
		c = 0.5*(edges[i]+edges[i+1])
		dat.append((edges[i],edges[i+1],0.0,0.0))
	writeUncertaintyTable("gainfluctsUncertainty.txt",dat)
							
if __name__=="__main__":
	linearityUncertaintyTable()
	#gainFluctsUncertaintyTable()
