#!/sw/bin/python2.7

import re
from string import *
from math import *
import os

from QFile import *

def fmt_err(xdx):
	x,dx = xdx
	x = float(x)
	dx = float(dx)
	
	if not x and not dx:
		return "---"
	
	ndig = floor(log10(dx)-0.5)
	if ndig < -5:
		ndig = -5
	idx = int(dx*10**(-ndig))
	fstr = "%%.%if"%(-ndig)
	if idx:
		return (fstr+"(%i)")%(x,idx)
	return fstr%x

def divErr(n,dn,d,dd):
	return n/d,sqrt((dn/d)**2+(n*dd/d**2)**2)


def sumWeightedErr(l):
	e = sum([x[0]*x[1] for x in l])
	w = sum([x[1] for x in l])
	e /= w
	serr = sum([((x[0]-e)*x[2])**2 for x in l])
	return e,sqrt(serr)/w
	

class AtomBindingTable:
	def __init__(self):
		self.dat = dict([(int(x.getFirst("Z")),x) for x in QFile(os.environ["UCNA_AUX"]+"/NuclearDecays/ElectronBindingEnergy.txt").dat["binding"]])
	def binding(self,Z,shell):
		if shell=="K1":
			shell = "K"
		if not Z in self.dat:
			print "Binding Z",Z,"not found!"
			return 0
		e = self.dat[Z].getFirstF(shell)
		if e is None:
			print "Binding",Z,shell,"missing!"
		return e

class BrIccDat:
	def __init__(self,Z):
		self.Z = Z
		self.E = self.dE = ''
		self.mx = self.dmx = ''
		
	def display(self):
		print "BrIcc Data:",self.Z,self.Egamma,self.M,"%s(%s)"%(self.mx,self.dmx)
		
	def Egamma(self):
		return "%s(%s)"%(self.E,self.dE)

	def Multipolarity(self):
		mm = self.M
		if self.mx:
			mm += " $\\delta=%s(%s)$"%(self.mx,self.dmx)
		return mm
		
	def ShellFrac(self,S):
		for r in self.rows:
			if r[0]==S:
				return float(r[self.IccCol]),float(r[self.dIccCol])
		return 0,0

	def ShellRat(self,Sn,Sd):
		n,dn = self.ShellFrac(Sn)
		d,dd = self.ShellFrac(Sd)
		if not d:
			return 0,0
		return divErr(n,dn,d,dd)

	def calcEnergy(self,ABT):
		lns = []
		shells = [("K",1),("L",3),("M",5),("N",7),("O",5),("P",2)]
		for s in shells:
			for i in range(s[1]):
				shnm = s[0]
				if s[1]>1:
					shnm += "%i"%(i+1)
				x,dx = self.ShellFrac(shnm)
				if not x:
					continue
				b = ABT.binding(self.Z,shnm)
				if not b:
					print "Missing",self.Z,shnm
					continue
				lns.append((float(self.E)-b/1000.,x,dx))
		return sumWeightedErr(lns)







def nextHeader(f):
	"""Find and load next block header"""
	hdrStart = "BrIcc v2.3a"
	l = ''
	x = None
	while not x:
		l = f.readline()
		print "###",l[:-1]
		if not l:
			break
		x = re.match(' BrIcc v2.3a.*Z=\s*(\d*)\s*Egamma=\s*(\d*\.*\d*)\s+(\d*).*Multipolarity=\s*(\S*)',l)
	if not x:
		return None
	Z,E,dE,M = x.groups()
	B = BrIccDat(int(Z))
	B.E,B.dE,B.M = E,dE,M

	x = re.search('Mixing ratio=\s*(\S*)\s*(\S*)',f.readline())
	if x:
		B.mx, B.dmx = x.groups()

	B.display()
	return B

def readShelltab(f):
	l = ' '
	while len(l)<2 or l[1]==' ':
		l = f.readline()
		print "###",l[:-1]
	headers = l.split()
	print "Headers:",headers
	if len(headers) < 3:
		return [],[]
	rows = []
	while "Compare" not in l:
		l = f.readline()
		s = l.split()
		if len(s) == len(headers):
			rows.append(s)
			print s
	return headers,rows


def Read_BrIcc_lst(fname):
	"""Read info from BrIcc.lst file"""
	f = open(fname,"r")
	ldat = []
	while True:
		B = nextHeader(f)
		if not B:
			break
		B.headers,B.rows = readShelltab(f)
		if len(B.headers)<3:
			continue
		B.IccCol = 1
		if "Icc" in B.headers:
			B.IccCol = B.headers.index("Icc")
		B.dIccCol = B.headers.index("dIcc")
		ldat.append(B)
	return ldat

def GigantoTable(tcols):

	znms = {47:"Ag",49:"In",56:"Ba",57:"La",82:"Pb"}
	ABT = AtomBindingTable()
	
	lx = "{|c|"+"|l"*len(tcols)+'|} \hline'+"\n"
	lx += "\tZ\t" + join(["& %i %s"%(c.Z,znms[c.Z]) for c in tcols],"\t")+' \\\\ '+"\n"
	lx += "\t$E_\\gamma$\t" + join(["& "+c.Egamma() for c in tcols],"\t")+" \\\\\n"
	lx += "\ttype\t" + join(["& "+c.Multipolarity() for c in tcols],"\t")+" \\\\\n"
	lx += "\\hline \\hline\n"
	lx += "\tK\t" + join(["& "+fmt_err(c.ShellFrac("K")) for c in tcols],"\t")+' \\\\ \hline'+"\n"
	
	shells = [("L",3),("M",5),("N",7),("O",5),("P",2)]
	for s in shells:
		lx += "\t%s$_\\text{tot}$\t"%s[0] + join(["& "+fmt_err(c.ShellFrac("%s-tot"%s[0])) for c in tcols],"\t")+' \\\\'+"\n"
		for i in range(s[1])[1:]:
			lx += "\t%s$_%i$/%s$_1$\t"%(s[0],i+1,s[0]) + join(["& "+fmt_err(c.ShellRat("%s%i"%(s[0],i+1),"%s1"%s[0])) for c in tcols],"\t")+' \\\\'+"\n"
		lx += "\\hline\n"
	lx += "\hline\t$\\alpha_\\text{tot}$\t" + join(["& "+fmt_err(c.ShellFrac("Tot")) for c in tcols],"\t")+' \\\\'+"\n"
	lx += "\t$\\langle E_e \\rangle$, keV\t" + join(["& "+fmt_err(c.calcEnergy(ABT)) for c in tcols],"\t")+' \\\\'+"\n"
	lx += "\\hline\n"
	print
	print
	print lx
	return lx

#\begin{tabular}{|c||c|ccc|ccc|}

def BETable():
	ABT = AtomBindingTable()
	shells = [("K",1),("L",3),("M",5)]
	
	atoms = [(47,"Ag"),(49,"In"),(56,"Ba"),(57,"La"),(82,"Pb")]
	
	lx = "{|c||"
	lx2 = "\tZ "
	lxa = ["\t%i %s"%a for a in atoms]
	
	for s in shells:
		lx += 'c'*s[1]+'|'
		for i in range(s[1]):
			lx2 += "\t& "+s[0]
			if s[1]>1:
				lx2 += "$_%i$"%(i+1)
			for (n,a) in enumerate(atoms):
				lxa[n] += "\t& %.2f"%(ABT.binding(a[0],s[0]+"%i"%(i+1))/1000.)
	lx += "} \\hline \n"+lx2+" \\\\ \\hline \\hline \n"
	for l in lxa:
		lx += l + " \\\\ \n"
	print
	print
	print lx


if __name__=="__main__":
	LBi = Read_BrIcc_lst("/Users/michael/Documents/UCNA/Papers/RefArticles/RadSources/ENSDF/BrIcc_207Bi-207PB.txt")
	LSn = Read_BrIcc_lst("/Users/michael/Documents/UCNA/Papers/RefArticles/RadSources/ENSDF/BrIcc_113Sn-113IN.txt")
	LCd = Read_BrIcc_lst("/Users/michael/Documents/UCNA/Papers/RefArticles/RadSources/ENSDF/BrIcc_109Cd-109AG.txt")
	LCe = Read_BrIcc_lst("/Users/michael/Documents/UCNA/Papers/RefArticles/RadSources/ENSDF/BrIcc_139Ce-139LA.txt")
	LIn = Read_BrIcc_lst("/Users/michael/Documents/UCNA/Papers/RefArticles/RadSources/ENSDF/BrIcc_114mIn-114IN.txt")
	LCs = Read_BrIcc_lst("/Users/michael/Documents/UCNA/Papers/RefArticles/RadSources/ENSDF/BrIcc_137Cs-137BA.txt")
	#GigantoTable(LCs)
	GigantoTable(LCd+LCe+LIn+[LSn[1],LBi[1],LCs[0],LBi[3]])

	BETable()

