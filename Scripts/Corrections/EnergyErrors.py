#!/usr/bin/python

import sys
sys.path.append("..")
from PyxUtils import *
from ErrorEnvelope import *

def bin_edges(w=10,n=100):
	return [w*i for i in range(n+1)]

baseOutPath = "../../Aux/Corrections/"

def writeUncertaintyTable(fout,dat):
	"""Write uncertainties table to file in correct format"""
	fout.write("#E_lo\tE_hi\tcorrection\tuncertainty\n")
	for d in dat:
		fout.write("%i\t%i\t%f\t%f\n"%d)

def linearityUncertaintyTable():
	"""Uncertainty due to energy calibration errors; see EnergyErrorsRevis.pdf"""
	edges = bin_edges()
	dat = []
	def neutronBeta(KE, m = 511.):
		return sqrt(KE**2+2*m*KE)/(m+KE)
	for i in range(len(edges)-1):
		c = 0.5*(edges[i]+edges[i+1])
		Eprim = c+calEnvelope(c)
		dat.append((edges[i],edges[i+1],0.0,neutronBeta(Eprim)/neutronBeta(c)-1))
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