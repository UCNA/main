#!/usr/bin/python

from EnergyErrors import *
from BetaSpectrumCorrections import *

def GenRecoilCorrTable():
	"""Generate Recoil/Weak Magnetism corrections table"""
	fin = os.environ["UCNA_ANA_PLOTS"]+"/SpectrumCorrection/SpectrumCorrection_1_1_782.347.txt"
	corrs = [spectrumCorrs(m) for m in QFile(fin).dat["spectrumPoint"]]
	RWMs = [(c.energy,1/(1+c.RWM)-1) for c in corrs]
	RWMs.sort()
	RWMs = RWMs[5::10]
	de = RWMs[1][0]-RWMs[0][0]
	be = bin_edges(w=de,n=len(RWMs))
	
	fout = open(baseOutPath+"/RecoilOrder.txt","w")
	fout.write("# Recoil and weak magnetism theory corrections\n#\n")
	writeUncertaintyTable(fout,[[be[n],be[n+1],RWMs[n][1],0.0] for n in range(len(RWMs))])

def GenRadiativeCorrections():
	fin = os.environ["UCNA_ANA_PLOTS"]+"/SpectrumCorrection/SpectrumCorrection_1_1_782.347.txt"
	corrs = [spectrumCorrs(m) for m in QFile(fin).dat["spectrumPoint"]]
	RCs = [(c.energy,1/(1+c.hmg)-1) for c in corrs]
	RCs.sort()
	RCs = RCs[5::10]
	de = RCs[1][0]-RCs[0][0]
	be = bin_edges(w=de,n=len(RCs))
	
	fout = open(baseOutPath+"/Radiative_h-g.txt","w")
	fout.write('# Shann Radiative "h-g" corrections\n#\n')
	writeUncertaintyTable(fout,[[be[n],be[n+1],RCs[n][1],0.0] for n in range(len(RCs))])

if __name__ == "__main__":
	GenRecoilCorrTable()
	GenRadiativeCorrections()