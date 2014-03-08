#!/usr/bin/python

import sys
sys.path.append("..")
from ucnacore.EncalDB import *
from ucnacore.RunAccumulatorFile import *
from review.Asymmetries import *
from Sources import SourceDatDirectory
from math import *
from ucnacore.QFile import *
from ucnacore.PyxUtils import *
import numpy
from numpy import zeros,matrix,linalg
import os

class wireSpec:
	def __init__(self,pos,nm,n):
		self.pos = pos
		self.nm = nm
		self.n = n
		self.norm = 1.0

class cathnorm_fits(KVMap):
	"""cathnorm_fits from CathodeGainPlugin analysis, indicating shape of cathode charge distribution"""
	def __init__(self,m=KVMap()):
		self.dat = m.dat
		self.loadFloats(["center","d_center","height","d_height","width","d_width","cathode","prev_gain"])
		self.cathode = int(self.cathode)
		for i in range(3):
			self.loadFloats(["cx_lo_a%i"%i,"d_cx_lo_a%i"%i,"cx_hi_a%i"%i,"d_cx_hi_a%i"%i])
		self.loadStrings(["side","plane"])
	def cathkey(self):
		return (self.side,self.plane,self.cathode)
	def getY(self,x):
		if -0.6 <= x <= -0.4:
			return self.cx_lo_a0 + self.cx_lo_a1*x + self.cx_lo_a2*x**2
		if 0.4 <= x <= 0.6:
			return self.cx_hi_a0 + self.cx_hi_a1*x + self.cx_hi_a2*x**2
		return None
	def getDeriv(self,x):
		if -0.6 <= x <= -0.4:
			return self.cx_lo_a1 + 2*self.cx_lo_a2*x
		if 0.4 <= x <= 0.6:
			return self.cx_hi_a1 + 2*self.cx_hi_a2*x
		return None

class CathnormFile(QFile):
	"""QFile containing 'cathnorm_fits' keys"""
	def __init__(self,fname):
		QFile.__init__(self,fname)
		self.cathnorms = dict([(c.cathkey(),c) for c in [cathnorm_fits(c) for c in self.dat.get("cathnorm_fits",[])]])

class cathseg(KVMap):
	"""CathodeSeg struct created by CathodeTweakPlugin analysis, indicating counts vs data and derivatives at boundary"""
	def __init__(self,m=KVMap()):
		self.dat = m.dat
		self.loadFloats(["position","n_exp","n_obs","dndx_lo","dndx_hi","i"])
		self.i = int(self.i)
		self.loadStrings(["side","plane"])
	def __eval__(self,x):
		return self.height*exp(-(x-self.center)**2/(2*self.width**2))
	def deriv(self,x):
		return -self(x)*x/self.width**2
	def expansionFactor(self,x0,x1):
		return exp(-(x0**2-x1**2)/(2*self.width**2))
	def renorm(self):
		return self.expansionFactor(0.5,0.5/self.fill_frac)
	
	def cathkey(self):
		return (self.side,self.plane,self.i)

class hitdist(KVMap):
	"""Fourier decomposition of hit distribution around cathode"""
	def __init__(self,m=KVMap()):
		self.dat = m.dat
		self.loadFloats(["cathode","ehi","elo","eavg"])
		self.cathode = int(self.cathode)
		self.loadStrings(["side","plane"])
		self.terms = [float(x) for x in self.getFirst("terms","").split(",")]
		self.dterms = [float(x) for x in self.getFirst("dterms","").split(",")]
	def cathkey(self):
		return (self.side,self.plane,self.cathode)

class mwpcGainCal(KVMap):
	"""MWPC gain calibration fit data"""
	def __init__(self,m=KVMap()):
		self.dat = m.dat
		self.loadFloats(["fit_gain","d_fit_gain","g0_avg","g0_d_avg","type"])
		self.type = int(self.type)
		self.loadStrings(["side"])

class CathFillFile(QFile):
	"""QFile containing cathseg entries with fill fraction vs MC"""
	def __init__(self,fname):
		QFile.__init__(self,fname)
		self.cathsegs = dict([(c.cathkey(),c) for c in [cathseg(c) for c in self.dat.get("cathseg",[])]])
		# cathode normalization from runcal
		self.cnorms = {}
		for s in ["E","W"]:
			for d in ["x","y"]:
				self.cnorms[(s,d)] = [float(x) for x in self.getItem("runcal","cnorm_"+s+d).split(",")]
		for k in self.cathsegs:
			self.cathsegs[k].cnorm = self.cnorms[(k[0],k[1])][k[2]]
	def getCathkeys(self,side,plane):
		ckeys = [ k for k in  self.cathsegs if k[0] == side and k[1] == plane]
		ckeys.sort()
		return ckeys

class CathFile(QFile):
	def __init__(self,fname):
		QFile.__init__(self,fname)
		self.cathsegs = {}
		for m in self.dat.get("cathseg",[]):
			c = cathseg(m)
			self.cathsegs[(c.side,c.plane,c.i)] = c
		self.hitdists = {}
		for m in self.dat.get("hitdist",[]):
			h = hitdist(m)
			self.hitdists.setdefault((h.side,h.plane,h.cathode),[]).append(h)
		self.cnorms = {}
		for s in ["E","W"]:
			for d in ["x","y"]:
				self.cnorms[(s,d)] = [float(x) for x in self.getItem("runcal","cnorm_"+s+d).split(",")]
		
class count23(KVMap):
	def __init__(self,m=KVMap()):
		self.dat = m.dat
		self.loadFloats(["e0","e1","II","IIx","III","IIIx"])
		self.loadStrings(["side"])
		self.total = self.II+self.IIx+self.III+self.IIIx


###############
# DB commands
###############

def delete_shape_graphs(conn,csid):
	print "Deleting shape graphs for",csid
	conn.execute("SELECT graph_id FROM cathshape_graphs WHERE cathseg_id = %i"%csid)
	for gid in [r[0] for r in conn.fetchall()]:
		delete_graph(conn,gid)
	conn.execute("DELETE FROM cathshape_graphs WHERE cathseg_id = %i"%csid)

def delete_cathseg_cal(conn,csid):
	print "Deleting cathseg_cal",csid
	delete_shape_graphs(conn,csid)
	conn.execute("DELETE FROM cathseg_cal WHERE cathseg_id = %i"%csid)
	
def delete_cathcal_set(conn,ccsid):
	print "Deleting cathode calibration set",ccsid
	conn.execute("SELECT cathseg_id FROM cathseg_cal WHERE cathcal_set_id = %i"%ccsid)
	for csid in [r[0] for r in conn.fetchall()]:
		delete_cathseg_cal(conn,csid)
	conn.execute("DELETE FROM cathcal_set WHERE cathcal_set_id = %i"%ccsid)

def new_cathcal_set(conn,side,plane,r0,r1):
	conn.execute("INSERT INTO cathcal_set(side,plane,start_run,end_run) VALUES ('%s','%s',%i,%i)"%(side,plane,r0,r1))
	conn.execute("SELECT LAST_INSERT_ID()")
	ccsid = int(conn.fetchone()[0])
	print "Generating cathcal set for",side,plane,r0,r1,":",ccsid
	return ccsid

def find_cathcal_sets(conn,r0,r1):
	conn.execute("SELECT cathcal_set_id FROM cathcal_set WHERE start_run = %i AND end_run = %i"%(r0,r1))
	return [r[0] for r in conn.fetchall()]

def new_cathseg_cal(conn,ccsid,wspec):
	conn.execute("INSERT INTO cathseg_cal(cathcal_set_id,position,sensor_name,norm) VALUES (%i,%f,'%s',%f)"
				 %(ccsid,wspec.pos,wspec.nm,wspec.norm))
	conn.execute("SELECT LAST_INSERT_ID()")
	return int(conn.fetchone()[0])

def set_cathshape_graph(conn,csid,gid):
	conn.execute("INSERT INTO cathshape_graphs(graph_id,cathseg_id) VALUES (%i,%i)"%(gid,csid))

###############
# wires info
###############

def getWires(rn,s,d):
	"""Get list of wires active for run/side/plane"""	
	assert s in ["East","West"] and d in ["X","Y"] and rn >= 13000
	nWires = 16
	wireSpacing = 4*2.54*sqrt(0.6)
	#								0    1    2    3    4    5    6    7    8    9    10   11   12   13   14   15
	padc_nums = {
					("East","X"): ( 231, 230, 229, 228, 227, 226, 225, 224, 223, 222, 221, 220, 219, 218, 217, 216 ),
					("East","Y"): ( 215, 214, 213, 212, 211, 210, 29,  28,  27,  26,  25,  24,  23,  22,  21,  20  ),
					("West","X"): ( 16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31  ),
					("West","Y"): ( 15,  14,  13,  12,  11,  10,  9,   8,   7,   6,   5,   4,   3,   2,   1,   0 ) }
	wires = []
	for i in range(nWires):
		pos = (-nWires*0.5+0.5+i)*wireSpacing
		nm = "Pdc%i"%padc_nums[(s,d)][i]
		if s=="West":
			nm = "Padc%i"%padc_nums[(s,d)][i]
		wires.append(wireSpec(pos,nm,i))
	return wires


###############
# cathode calibrations
###############

def gen_cathcal_set(conn,r0,r1,cr0,cr1):
	"""Plain cathcal set with no shape corrections"""
	
	basepath = os.environ["UCNA_ANA_PLOTS"]+"/WirechamberCal/%i-%i"%(cr0,cr1)
	cathnorms = PlotCathNorm(basepath)
	corrcurves = PlotCathShapes(basepath)
	
	if not conn:
		return
	
	for ccsid in find_cathcal_sets(conn,r0,r1):
		delete_cathcal_set(conn,ccsid)
	
	for s in ["East","West"]:
		for d in ["X","Y"]:
			ccsid = new_cathcal_set(conn,s,d,r0,r1)
			for w in getWires(r0,s,d):
				w.norm = cathnorms[(s[0],d.lower(),w.n)]
				csid = new_cathseg_cal(conn,ccsid,w)
				for c in corrcurves[(s[0],d.lower(),w.n)]:
					gid = upload_graph(conn,"Cathode %s%s %i"%(s,d,w.n),c)
					set_cathshape_graph(conn,csid,gid)

def upload_cathcal_set(conn,r0,r1,cathnorms,corrcuves=None):
	
	if not corrcuves:
		corrcurves = dict([(k,[]) for k in cathnorms])
	
	for ccsid in find_cathcal_sets(conn,r0,r1):
		delete_cathcal_set(conn,ccsid)
	
	for s in ["East","West"]:
		for d in ["X","Y"]:
			ccsid = new_cathcal_set(conn,s,d,r0,r1)
			for w in getWires(r0,s,d):
				w.norm = cathnorms[(s[0],d.lower(),w.n)]
				csid = new_cathseg_cal(conn,ccsid,w)
				for c in corrcurves[(s[0],d.lower(),w.n)]:
					gid = upload_graph(conn,"Cathode %s%s %i"%(s,d,w.n),c)
					set_cathshape_graph(conn,csid,gid)


# Cathode hit distribution correction factors
def PlotCathShapes(basepath):
	CF = CathFile(basepath+"/MWPCCal.txt")
	caths = CF.hitdists.keys()
	caths.sort()
	corrcurves = {}
	
	lwidths = {0:[style.linewidth.thin],1:[style.linewidth.Thick]}
	#lstyles = {0:[],1:[style.linestyle.dashed],2:[style.linestyle.dotted]}
	lstyles = {0:[style.linestyle.dotted],1:[],3:[style.linestyle.dashed]}
	#lcols = {0:[rgb.blue],1:[rgb.red]}
	lcols = {0:[],1:[]}
	
	for k in caths:
		print k
		gdat = [ [h.eavg,]+h.terms[1:]+h.dterms[1:]+[h,] for h in CF.hitdists[k]]
		gdat.sort()
		nterms = (len(gdat[0])-1)/2
		tcols = rainbow(nterms)
						  
		gHitDist=graph.graphxy(width=15,height=10,
						x=graph.axis.lin(title="Scintillator Visible Energy [keV]",min=0,max=1000),
						y=graph.axis.lin(title="%s%s-%i Distortion Coeffiecient"%k,min=-0.2,max=0.2),
						key = graph.key.key(pos="bc",columns=3))
		setTexrunner(gHitDist)
		
		corrcurves[k] = []
		for n in range(nterms):
			tdat = [(g[0],0,g[n+1],g[nterms+n+1]) for g in gdat]
			# intermediate points
			zdat = [(gdat[0][-1].elo,0,tdat[0][2],0),tdat[0]]
			for i in range(len(tdat))[:-1]:
				zdat += [(gdat[i][-1].ehi,0,0.5*(tdat[i][2]+tdat[i+1][2]),0),tdat[i+1]]
			zdat += [(gdat[-1][-1].ehi,0,tdat[-1][2],0)]
			corrcurves[k].append(zdat)
			gtitle = "$\\sin(%i\\pi x)$"%(2*(n/2+1))
			if n%2:
				gtitle = "$\\cos(%i\\pi x)$"%(2*(n/2+1))
			gHitDist.plot(graph.data.points(tdat,x=1,y=3,dy=4,title=None),
						  [graph.style.errorbar(errorbarattrs=lcols[n%2])])
			gHitDist.plot(graph.data.points(zdat,x=1,y=3,dy=4,title=gtitle),
						  [graph.style.line(lineattrs=lstyles[n/2]+lwidths[n%2]+lcols[n%2])])
		if False:
			try:
				gHitDist.writetofile(basepath+"/HitDist_%s%s_%i.pdf"%k)
			except:
				print "******** PLOTTING ERROR! *********"

	# replacement for end wires; average distribution
	avgcurves = {}
	cavg = {}
	for s in ['E','W']:
		for d in ['x','y']:
						
			gAvgDist=graph.graphxy(width=15,height=10,
								   x=graph.axis.lin(title="Scintillator Visible Energy [keV]",min=0,max=1000),
								   y=graph.axis.lin(title="%s%s Average Coefficient"%(s,d),min=-0.1,max=0.2),
								   key = graph.key.key(pos="bc",columns=3))
			setTexrunner(gAvgDist)

			avgcurves[(s,d)] = []
			cavg[(s,d)] = []
			for nt in range(len(corrcurves[(s,d,8)])):
				avgcurves[(s,d)].append([])
				for ne in range(len(corrcurves[(s,d,8)][nt])):
					e = corrcurves[(s,d,8)][0][ne][0]
					#(mu,sg) = musigma([corrcurves[(s,d,i)][nt][ne][2] for i in range(16)[2:14]])
					(mu,sg) = musigma([corrcurves[(s,d,i)][nt][ne][2] for i in range(16)[3:13]])
					avgcurves[(s,d)][-1].append([e,0,mu,sg*(ne%2)])
				cavg[(s,d)].append(musigma([g[2] for g in avgcurves[(s,d)][-1]])[0])
				gtitle = "$\\sin(%i\\pi x)$"%(2*(nt/2+1))
				if nt%2:
					gtitle = "$\\cos(%i\\pi x)$"%(2*(nt/2+1))
				if nt not in [0,1,3]:
					continue
				gAvgDist.plot(graph.data.points(avgcurves[(s,d)][-1],x=1,y=3,dy=4,title=gtitle),
							  #[graph.style.line(lineattrs=lstyles[nt/2]+lwidths[nt%2]+lcols[nt%2]),
							  [graph.style.line(lineattrs=lstyles[nt]+lwidths[nt%2]+lcols[nt%2]),
							   graph.style.errorbar(errorbarattrs=lcols[nt%2])])

			gAvgDist.writetofile(basepath+"/HitDist_%s%s.pdf"%(s,d))
		
			corrcurves[(s,d,0)]=[]
			corrcurves[(s,d,15)]=[]
	
			corrcurves[(s,d,1)]=avgcurves[(s,d)]
			corrcurves[(s,d,14)]=avgcurves[(s,d)]
	
			if 0:
				for i in range(16)[1:15]:
					c1avgDel = musigma([g[2] for g in corrcurves[(s,d,i)][1]])[0]-cavg[(s,d)][1]
					print s,d,i,"Corrected by",c0avgDel
					corrcurves[(s,d,i)] = avgcurves[(s,d)]
					corrcurves[(s,d,i)][1] = [(g[0],g[1],g[2]+c1avgDel,g[3]) for g in avgcurves[(s,d)][1]]
				
	return corrcurves

def CathNormCorrection(basePath, normFile, fillFile):

	CN = CathnormFile(normFile)
	CF = CathFillFile(fillFile)
	
	# fill fractions
	gFill=graph.graphxy(width=15,height=10,
						 x=graph.axis.lin(title="cathode number",min=0,max=15,
										  parter=graph.axis.parter.linear(tickdists=[5,1])),
						 y=graph.axis.lin(title="observed/expected counts",min=0.9,max=1.12),
						 key = graph.key.key(pos="bc",columns=2))
	setTexrunner(gFill)
	
	# cathode norms after correction
	gCNorm=graph.graphxy(width=15,height=10,
						   x=graph.axis.lin(title="cathode number",min=0,max=15,
											parter=graph.axis.parter.linear(tickdists=[5,1])),
						   y=graph.axis.lin(title="gain change $dg$"),
						   key = graph.key.key(pos="tc",columns=2))
	setTexrunner(gCNorm)
	
	lstyles = {"x":[],"y":[style.linestyle.dashed]}
	lcols = {"E":[rgb.red],"W":[rgb.blue]}
	ssymbs = {"E":symbol.circle,"W":symbol.triangle}
	
	# all cathode gain normalizations recommended
	cathnorms = {}
	numpy.set_printoptions(precision=3, linewidth=250)
	
	for s in ["E","W"]:
		for d in ["x","y"]:
			print "-----------",s,d,"------------"
			
			# collect data
			ckeys = CF.getCathkeys(s,d)[1:-1]
			nCaths = len(ckeys)
			cathSegs = [CF.cathsegs[k] for k in ckeys]
			cathNorms = [CN.cathnorms[k] for k in ckeys]
			n_obs = sum([c.n_obs for c in cathSegs])
			n_exp = sum([c.n_exp for c in cathSegs])
			print "Total events Expected:",n_exp,"Observed:",n_obs
			gdat = [ (c.i,c.n_obs/c.n_exp*n_exp/n_obs) for c in cathSegs]
			print "Fill factors",gdat
			
			gtitle = "%s%s"%(s,d)
			gFill.plot(graph.data.points(gdat,x=1,y=2,title=gtitle),
						[graph.style.line(lineattrs=lstyles[d]+lcols[s]),
						 graph.style.symbol(ssymbs[s],symbolattrs=lcols[s])])

			# check crossover consistency
			print "Cx consistency",  [cathNorms[i].getY(0.5)-cathNorms[i+1].getY(-0.5) for i in range(nCaths-1)]
			# dx/dgain at each crossing
			dxdg = [0.5*(cathNorms[i].getY(0.5)+cathNorms[i+1].getY(-0.5))/(-cathNorms[i].getDeriv(0.5)+cathNorms[i+1].getDeriv(-0.5)) for i in range(nCaths-1)]
			print "Crossing dx/dg",dxdg
			# dn/dx at each crossing, averaged between cathodes
			dndx = [0.5*(cathSegs[i].dndx_hi+cathSegs[i+1].dndx_lo) for i in range(nCaths-1)]
			print "Crossing dn/dx",dndx
			# desired delta n expected-observed
			deltaN = matrix([ c.n_exp*n_obs/n_exp - c.n_obs for c in cathSegs]).transpose()
			print "Desired delta N"
			print deltaN
			
			# matrix: deltaN sensitivity to gain changes
			gmatrix = matrix(zeros((nCaths,nCaths)))
			for i in range(nCaths):
				dndg = 0
				if i>0:
					gmatrix[i-1,i] -= dxdg[i-1]*dndx[i-1]
					dndg -= gmatrix[i-1,i]
				if i<nCaths-1:
					gmatrix[i+1,i] -= dxdg[i]*dndx[i]
					dndg -= gmatrix[i+1,i]
				gmatrix[i,i] += dndg
			print "gain sensitivity [counts]"
			print gmatrix
		
			# matrix: transform gmatrix into local imbalance system
			ncomb = matrix(zeros((nCaths-1+nCaths,nCaths)))
			for i in range(nCaths-1):
				#ncomb[i,i] = -1.0/cathSegs[i].n_exp
				#ncomb[i,i+1] = 1.0/cathSegs[i+1].n_exp
				ncomb[i,i] = -1.0/n_exp
				ncomb[i,i+1] = 1.0/n_exp
			print "combining matrix"
			print ncomb
			
			gmatrix = ncomb*gmatrix
			deltaN = ncomb*deltaN
			#gmatrix[nCaths-1,3] = gmatrix[nCaths,12] = gmatrix[nCaths+1,7] = 10.0
			for i in range(nCaths):
				gmatrix[nCaths-1+i,i] = 0.005
			print "gain equation matrix"
			print gmatrix
			
			(gain,ssr,matrixRank,singularValues) = linalg.lstsq(gmatrix,deltaN)
			print "recommended gain change"
			print gain
			print "rms dev =",sqrt(ssr/len(deltaN))

			gdat = [(n+1,c) for (n,c) in enumerate(gain)]
			gCNorm.plot(graph.data.points(gdat,x=1,y=2,title=gtitle),
						[graph.style.line(lineattrs=lstyles[d]+lcols[s]),
						 graph.style.symbol(ssymbs[s],symbolattrs=lcols[s])])

			for (n,k) in enumerate(ckeys):
				cathnorms[k] = 1.0+gain[n]
			cathnorms[(s,d,0)] = cathnorms[(s,d,15)] = 1.0

	print cathnorms
	
	gFill.writetofile(basePath+"/CathFill.pdf")
	gCNorm.writetofile(basePath+"/CathNorm.pdf")

	return cathnorms

###############
# anode gain calibration
###############

class anodeCalAvg(KVMap):
	def __init__(self,m=KVMap()):
		self.dat = m.dat
		self.loadFloats(["avg","d_avg","type"])
		self.type = int(self.type)
		self.loadStrings(["side"])

class anodeGaincorr(KVMap):
	def __init__(self,m=KVMap()):
		self.dat = m.dat
		self.loadFloats(["avg","d_avg"])
		self.loadStrings(["side"])

def getAnodeCalAvgs(AF):
	return dict([((m.side,m.type),m) for m in [anodeCalAvg(m) for m in AF.dat["anodeCalAvg"]]])

class Anode_MC_Comparator(MC_Comparator):
	
	def __init__(self,basedir,simdir,depth=0):
		MC_Comparator.__init__(self,basedir,simdir,depth)

	def betaAnodeGainCal(self,conn=None,pmap=None):
		"""Set anode calibration to match data and simulation"""
		# collect data
		gdat = {}
		for (n,rns) in enumerate(self.rungrps):
			print n,(rns[0],rns[-1]),"-"*80
			# previous anode gain corrections
			prevagc = dict([(m.side,m) for m in [anodeGaincorr(m) for m in self.datAsyms[rns].dat["anodeGaincorr"]]])
			# data and simulation average values
			datAvgs = getAnodeCalAvgs(self.datAsyms[rns])
			simAvgs = getAnodeCalAvgs(self.simAsyms[rns])
			# plottable data by side, event type
			for s in ["E","W"]:
				for tp in range(3):
					adat = datAvgs[(s,tp)].avg
					asim = simAvgs[(s,tp)].avg
					gOld = prevagc[s].avg
					gdat.setdefault((s,tp),[]).append([n,adat,asim,asim/adat*gOld])
					print "\t%s%i:\t%f\t%f\t%f\tgOld: %f"%(s,tp,adat,asim,gdat[(s,tp)][-1][-1],gOld)
			# upload correction factors
			if conn:
				uploadAnodeCal(conn,rns[0],rns[-1],pmap,gdat[('E',0)][-1][-1],gdat[('W',0)][-1][-1])
	
		# plot
		gAGC=graph.graphxy(width=25,height=10,
						  x=graph.axis.lin(title=unitNames[self.depth],min=0,max=len(self.rungrps)-1),
						  y=graph.axis.lin(title="Anode calibration [keV/channel]"),
						  key = graph.key.key(pos="bl"))
		setTexrunner(gAGC)
		tpsymbs = {0:symbol.circle,1:symbol.triangle,2:symbol.plus}
		sideStyles = {"East":[scols["East"]], "West":[scols["West"],deco.filled]}
		for s in ["East","West"]:
			for tp in range(2):
				gtitle = s + " Type " + {0:"0",1:"I",2:"II/III"}[tp]
				gAGC.plot(graph.data.points(gdat[(s[0],tp)],x=1,y=-1,title=gtitle),[graph.style.symbol(tpsymbs[tp],symbolattrs=sideStyles[s])])
		print "Saving gain calibration plot..."
		gAGC.writetofile(self.basedir+"/AnodeGainCal.pdf")

		# side average values
		LF = LinearFitter(terms=[polyterm(0)])
		savg = {}
		for s in ['E','W']:
			LF.fit(gdat[(s,0)],cols=(0,-1))
			savg[s] = LF.coeffs[0]
		return savg
	
def anodeGainCal(conn=None):
	if(conn):
		deleteAllAnodeCals(conn)
	
	pmap = 145
			
	AMC = Anode_MC_Comparator(os.environ["UCNA_ANA_PLOTS"]+"/OctetAsym_Offic/",os.environ["UCNA_ANA_PLOTS"]+"/OctetAsym_Offic_SimMagF/",0)
	savg = AMC.betaAnodeGainCal(conn,pmap)
	if conn:
		uploadAnodeCal(conn,7000,100000,pmap,savg['E'],savg['W'])

#####################
# II/III Separation #
#####################

def sep23effic(basedir=os.environ["UCNA_ANA_PLOTS"]+"/test/23Separation/"):
	cdat = [ count23(c) for c in QFile(basedir+"/23Separation.txt").dat["23counts"]]
	
	for s in ["East","West"]:
		gdat = [ (c.e0,c.e1,c.total,c.II/c.total,c.IIIx/c.total,c.III/c.total,c.IIx/c.total,0.5*(c.e0+c.e1)) for c in cdat if c.side==s ]
		gdat.sort()
		
		g23=graph.graphxy(width=10,height=7,
						   x=graph.axis.lin(title="Energy [keV]",min=0,max=700),
						   y=graph.axis.lin(title="Fraction of II+III",min=0,max=1),
						   key = graph.key.key(pos="tc",columns=2))
		setTexrunner(g23)

		print gdat
		
		iiCol,iiiCol = rgb.red,rgb.blue
		iiCol,iiiCol = rgb.black,rgb.black
		
		g23.plot(graph.data.points(gdat,x=-1,y=4,title="II"),[graph.style.line(lineattrs=[iiCol,style.linestyle.dashed])])
		g23.plot(graph.data.points(gdat,x=-1,y=5,title="III*"),[graph.style.line(lineattrs=[iiiCol,style.linestyle.dashdotted])])
		g23.plot(graph.data.points(gdat,x=-1,y=6,title="III"),[graph.style.line(lineattrs=[iiiCol,style.linestyle.solid])])
		g23.plot(graph.data.points(gdat,x=-1,y=7,title="II*"),[graph.style.line(lineattrs=[iiCol,style.linestyle.dotted])])
		
		g23.writetofile(basedir+"frac23_%s.pdf"%s)


################################
# Wirechamber calibration sets #
################################

# make new MWPC position maps from energy and charge signal maps
def make_mwpc_posmap(conn,charge_pmid,energy_pmid):
	c_pset = getPosmapSet(conn,charge_pmid)
	e_pset = getPosmapSet(conn,energy_pmid)
	
	pinfo = c_pset["East",0].info
	pinfo.descrip = "MWPC Position Response Q/E=[%i]/[%i]"%(charge_pmid,energy_pmid)
	newPosmap(conn,pinfo)
	
	for s in ["East","West"]:
		c_pts = c_pset[(s,0)].get_pts_sorted()
		#e_pts = e_pset[(s,0)].get_pts_sorted()
		#assert len(c_pts)==len(e_pts)
		
		# average and standard deviation of energy fits
		emu,esig = musigma(e_pset[(s,0)].get_pt_vals())
		print "Energy",s,"mu =",emu," rms =",esig
		
		m = posmap(pinfo)
		m.side = s
		m.quadrant = 0
		for (i,c) in enumerate(c_pts):
			c.sig /= c.norm;
			c.norm = emu #e_pts[i].sig/e_pts[i].norm # use energy position map average
			m.add_pt(c)
		uploadPosmap(conn,m)
	
	return pinfo

# Wirechamber gain/position calibrations class
class MWPC_Ecal:
	def __init__(self):
		self.charge_meas="ccloud"
		self.gain=None
		self.side="BOTH"
		self.start_run = self.end_run = 0
		self.priority = 0
		self.pmid = None

	def upload(self,conn):
		if self.side == "BOTH":
			for self.side in ["East","West"]:
				self.upload(conn)
			self.side = "BOTH"
			return

		# delete old entries for same run range
		wh = "WHERE start_run = %i AND end_run = %i AND side = '%s' AND charge_meas='%s'"%(self.start_run,self.end_run,self.side,self.charge_meas)
		conn.execute("SELECT count(*) FROM mwpc_ecal "+wh)
		nold = conn.fetchone()[0]
		if nold:
			print "Deleting",nold,"previous entries."
			cmd = "DELETE FROM mwpc_ecal "+wh
			print cmd
			conn.execute(cmd)
		
		if self.pmid is None:
			return
		
		print "Setting up MWPC calibrations for",self.start_run,"-",self.end_run,self.side,"using",self.charge_meas
		if self.gain is None: # synthesize gain from position map
			self.gain = 1.0/getPosmapSet(conn,self.pmid)[(self.side,0)].avg_val()
			print "Applying average gain factor",self.gain,"from position map",self.pmid
		else:
			print "Gain factor",self.gain,"and position map",self.pmid

		cmd = "INSERT INTO mwpc_ecal (start_run,end_run,priority,side,charge_meas,gain_posmap_id,gain_factor) "
		cmd += "VALUES (%i,%i,%i,'%s','%s',%i,%g)"%(self.start_run,self.end_run,self.priority,self.side,self.charge_meas,self.pmid,self.gain)
		print cmd
		conn.execute(cmd)

	def __repr__(self):
		return "<MWPC Ecal %i-%i %s %s %i %g p%i>"%(self.start_run,self.end_run,self.side,self.charge_meas,self.pmid,self.gain,self.priority)


def delete_cath_ccloud_scale(conn,ccsid):
	print "Deleting cathode cloud scaling set",ccsid
	conn.execute("SELECT gain_graph_id FROM cath_ccloud_scale WHERE cath_ccloud_scale_id = %i"%ccsid);
	for r in conn.fetchall():
		delete_graph(conn,r[0])
	conn.execute("DELETE FROM cath_ccloud_scale WHERE cath_ccloud_scale_id = %i"%ccsid);

# Cathode CCloud scaling factor g_i*g_CC data class
class CthCclScl:
	def __init__(self,srun,erun,side,plane):
		self.srun = srun
		self.erun = erun
		self.side = side
		self.plane = plane
		self.cfacts = None
	def delete_from_db(self,conn):
		conn.execute("SELECT cath_ccloud_scale_id FROM cath_ccloud_scale WHERE start_run = %i AND end_run = %i AND side = '%s' AND plane = '%s'"%(self.srun,self.erun,self.side,self.plane))
		for r in conn.fetchall():
			delete_cath_ccloud_scale(conn,r[0])
	def upload_to_db(self,conn):
		assert self.cfacts is not None
		self.gid = upload_graph(conn,"Cathode CCloud Scaling %s %s"%(self.side,self.plane),self.cfacts)
		conn.execute("INSERT INTO cath_ccloud_scale (start_run,end_run,side,plane,gain_graph_id) VALUES (%i,%i,'%s','%s',%i)"%(self.srun,self.erun,self.side,self.plane,self.gid))


##################
# Functions for updating from Cal. DB prior to full cathode gain calibrations system

# transfer deprecated anode calibration table to new system
def assign_MWPC_calib_from_anode_table(conn):
	conn.execute("SELECT start_run,end_run,anode_posmap_id,calfactor_E,calfactor_W FROM anode_cal")
	for r in conn.fetchall():
		M = MWPC_Ecal()
		M.start_run,M.end_run = r[0],r[1]
		M.pmid = r[2]
		M.charge_meas = "anode"
		for M.side in ["East","West"]:
			M.gain = {"East":r[3],"West":r[4]}[M.side]
			M.upload(conn)

# set approximate default cathode scale factors in DB
def set_default_cathode_scalefactors(conn):
	for s in ["East","West"]:
		for p in ["x","y"]:
			C = CthCclScl(0,100000,s,p)
			C.delete_from_db(conn)
			C.cfacts = [ [i,0,0.32,0] for i in range(16)]
			C.upload_to_db(conn)

##################
# New wirechamber calibrations

def MWPC_calib_from_RunAccumulatorFile(conn,RAcF,M):
	"""Generate adjusted MWPC energy calibrations from file"""
	calruns = RAcF.getRuns()
	gcals = dict([ ((c.side,c.type),c) for c in [mwpcGainCal(c) for c in RAcF.dat.get("mwpcGainCal",[])]])
	if ("East",0) not in gcals or ("West",0) not in gcals or not calruns:
		print calruns
		print "No data found in",RAcF.fname
		M.gain = None
		return None
		
	M.start_run = calruns[0]
	M.end_run = calruns[-1]
	
	sgains = {}
	for M.side in ["East","West"]:
		M.gain = gcals[(M.side,0)].g0_avg / gcals[(M.side,0)].fit_gain
		sgains[M.side] = M.gain
		print M
		if conn:
			M.upload(conn)
	return sgains

def MWPC_calib_for_beta_octets(conn, basedir, runAxis=None, gGC=None):
	"""Set MWPC ccloud gain factors for runs in each beta decay octet"""
	depth = 0
	
	M = MWPC_Ecal()
	M.priority = 10
	M.pmid = 181
	M.charge_meas = "ccloud"
	
	gdat = []
	for (n,f) in enumerate(collectOctetFiles(basedir,0)):
		gdat.append( (n, MWPC_calib_from_RunAccumulatorFile(conn,f,M), 0.5*(f.getRuns()[0]+f.getRuns()[-1])) )

	# plot
	if gGC is None:
		xaxis = graph.axis.lin(title=unitNames[depth],min=0,max=len(gdat)-1)
		if runAxis:
			xaxis = runAxis
		gGC=graph.graphxy(width=25,height=10,
						  x=xaxis,
						  y=graph.axis.lin(title="MWPC gain calibration $g_Q$",min=2e-4,max=3.5e-4),
						  key = graph.key.key(pos="bl"))
		setTexrunner(gGC)
	tpsymbs = {0:symbol.circle,1:symbol.triangle,2:symbol.plus}
	sideStyles = {"East":[scols["East"]], "West":[scols["West"],deco.filled]}
	for s in ["East","West"]:
		gtitle = s + " $\\beta$ octets"		# + " Type " + {0:"0",1:"I",2:"II/III"}[tp]
		gpts = [(g[0],g[1][s]) for g in gdat]
		if runAxis:
			gpts = [(g[-1],g[1][s]) for g in gdat]
		gGC.plot(graph.data.points(gpts,x=1,y=2,title=gtitle),[graph.style.symbol(tpsymbs[0],symbolattrs=sideStyles[s])])

	#gGC.writetofile(basedir+"/MWPC_GainCal_%i.pdf"%depth)
	return gGC

def MWPC_calib_for_source_runs(conn,basedir,gGC=None):
	"""Set MWPC ccloud gain factors for each calibration source run"""
	
	M = MWPC_Ecal()
	M.priority = 10
	M.pmid = 181
	M.charge_meas = "ccloud"

	SDD = SourceDatDirectory(basedir)
	gdat = [(rn,MWPC_calib_from_RunAccumulatorFile(conn,SDD.getQFile(rn,RunAccumulatorFile),M)) for rn in SDD.getRunlist()]
		
	# plot
	if gGC is None:
		gGC=graph.graphxy(width=25,height=10,
						  x=make_runaxis(gdat[0][0],gdat[-1][0]),
						  y=graph.axis.lin(title="MWPC gain calibration $g_Q$"),
						  key = graph.key.key(pos="bl"))
		setTexrunner(gGC)
	tpsymbs = {0:symbol.circle,1:symbol.triangle,2:symbol.plus}
	#sideStyles = {"East":[scols["East"]], "West":[scols["West"],deco.filled]}
	sideStyles = {"East":[rgb.red], "West":[rgb.red,deco.filled]}
	for s in ["East","West"]:
		gtitle = s + " source calibrations" # + " Type " + {0:"0",1:"I",2:"II/III"}[tp]
		gGC.plot(graph.data.points([(g[0],g[1][s]) for g in gdat if g[1]],x=1,y=2,title=gtitle),[graph.style.symbol(tpsymbs[0],symbolattrs=sideStyles[s])])
	gGC.writetofile(basedir+"/MWPC_GainCal.pdf")

	return gGC



def MWPC_Cathode_CCloud_Scaling(conn,datfName,simfName,start_run=0,end_run=100000):
	"""Set MWPC individual cathode to charge cloud size scaling factors"""
	
	dat = CathnormFile(datfName)
	sim = CathnormFile(simfName)
	
	for s in ["East","West"]:
		for p in ["x","y"]:
			C = CthCclScl(start_run,end_run,s,p.upper())
			C.cfacts = [ [i,0,0.32,0] for i in range(16)]
						
			print "Cathode gain factor changes",s,p
			for c in range(16)[1:-1]:
				delta = dat.cathnorms[(s[0],p,c)].height/sim.cathnorms[(s[0],p,c)].height
				C.cfacts[c][2] = delta * dat.cathnorms[(s[0],p,c)].prev_gain
				print c,delta

			# end cathodes lacking independent data... fill in with neighbor's values
			C.cfacts[0][2] = C.cfacts[1][2]
			C.cfacts[-1][2] = C.cfacts[-2][2]

			print C.cfacts
			print musigma([c[2] for c in C.cfacts[1:-1]])

			C.delete_from_db(conn)
			C.upload_to_db(conn)



###############
#             #
###############

if __name__ == "__main__":

	conn = open_connection()
	conn = None
	ppath = os.environ["UCNA_ANA_PLOTS"]
	
	#########
	# transition from DB before mwpc calibrations
	#########
	#assign_MWPC_calib_from_anode_table(conn)
	#set_default_cathode_scalefactors(conn)
	#exit(0)
	
	#########
	# 2010 cathode-based calibration
	#########
	
	#make_mwpc_posmap(conn,167,169) # only run this once
	#assign_MWPC_calib(conn,0,100000,10,181,"ccloud") # default calibration for all runs
	
	gGC = graph.graphxy(width=25,height=10,
					x=make_runaxis(14000,16400),
					y=graph.axis.lin(title="MWPC gain calibration $g_Q$",min=2e-4,max=3.5e-4),
					key = graph.key.key(pos="bl"))
	setTexrunner(gGC)

	MWPC_calib_for_beta_octets(conn, ppath+"/MWPC_ECal_8_Sim0823/", True, gGC)
	MWPC_calib_for_source_runs(conn, ppath+"/SourceFitsSim2010", gGC)
	
	#MWPC_Cathode_CCloud_Scaling(conn, ppath+"/MWPC_ECal_8/MWPC_ECal_8.txt", ppath+"/MWPC_ECal_8_Sim0823/MWPC_ECal_8_Sim0823.txt") # cathode to charge-cloud scaling
	
	exit(0)
	
	
	
	
	#anodeGainCal(None)
	
	#sep23effic(os.environ["UCNA_ANA_PLOTS"]+"/test/23Separation_Sim0823_4x/")
	
	#conn = open_connection()
	#for ccsid in find_cathcal_sets(conn,14077,16216):
	#	delete_cathcal_set(conn,ccsid)
	#anodeGainCal(conn)
	#conn = None
	#gen_cathcal_set(conn,13000,100000,14264,16077)

	wpath = ppath+"/WirechamberCal/14077-16216"
	cnorms = CathNormCorrection(wpath,ppath+"BetaOctetPositions_NewNorm/BetaOctetPositions_NewNorm.txt",wpath+"/MWPCCal.txt")
	#upload_cathcal_set(conn,14077,16216,cnorms)
