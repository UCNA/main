#!/usr/bin/python

from EncalDB import *
from math import *
from QFile import *
from PyxUtils import *
from Asymmetries import *
import os

class wireSpec:
	def __init__(self,pos,nm,n):
		self.pos = pos
		self.nm = nm
		self.n = n
		self.norm = 1.0

class cathseg(KVMap):
	def __init__(self,m=KVMap()):
		self.dat = m.dat
		self.loadFloats(["center","d_center","height","d_height","width","d_width","i","position","max","d_max","fill_frac"])
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

class hitdist(KVMap):
	def __init__(self,m=KVMap()):
		self.dat = m.dat
		self.loadFloats(["cathode","ehi","elo","eavg"])
		self.cathode = int(self.cathode)
		self.loadStrings(["side","plane"])
		self.terms = [float(x) for x in self.getFirst("terms","").split(",")]
		self.dterms = [float(x) for x in self.getFirst("dterms","").split(",")]

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

def deleteAllAnodeCals(conn):
	conn.execute("DELETE FROM anode_cal WHERE 1")

def uploadAnodeCal(conn,rmin,rmax,pmap,ecor,wcor):
	cmd = "INSERT INTO anode_cal(start_run,end_run,anode_posmap_id,calfactor_E,calfactor_W) VALUES (%i,%i,%i,%g,%g)"%(rmin,rmax,pmap,ecor,wcor)
	print cmd
	conn.execute(cmd)

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
	
	#return
	
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

def musigma(l):
	mu = sum(l)/len(l)
	s = sqrt(sum([x*x for x in l])/len(l)-mu*mu)/sqrt(len(l))
	return (mu,s)

def PlotCathShapes(basepath):
	CF = CathFile(basepath+"/MWPCCal.txt")
	caths = CF.hitdists.keys()
	caths.sort()
	corrcurves = {}
	
	lwidths = {0:[style.linewidth.thin],1:[style.linewidth.Thick]}
	lstyles = {0:[],1:[style.linestyle.dashed],2:[style.linestyle.dotted]}
	lcols = {0:[rgb.blue],1:[rgb.red]}
	
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
								   y=graph.axis.lin(title="%s%s Average Coefficient"%(s,d),min=-0.2,max=0.2),
								   key = graph.key.key(pos="bc",columns=nterms/2))
			setTexrunner(gAvgDist)

			avgcurves[(s,d)] = []
			cavg[(s,d)] = []
			for nt in range(len(corrcurves[(s,d,8)])):
				avgcurves[(s,d)].append([])
				for ne in range(len(corrcurves[(s,d,8)][nt])):
					e = corrcurves[(s,d,8)][0][ne][0]
					(mu,sg) = musigma([corrcurves[(s,d,i)][nt][ne][2] for i in range(16)[2:14]])
					avgcurves[(s,d)][-1].append([e,0,mu,sg*(ne%2)])
				cavg[(s,d)].append(musigma([g[2] for g in avgcurves[(s,d)][-1]])[0])
				gtitle = "$\\sin(%i\\pi x)$"%(2*(nt/2+1))
				if nt%2:
					gtitle = "$\\cos(%i\\pi x)$"%(2*(nt/2+1))
				gAvgDist.plot(graph.data.points(avgcurves[(s,d)][-1],x=1,y=3,dy=4,title=gtitle),
							  [graph.style.line(lineattrs=lstyles[nt/2]+lwidths[nt%2]+lcols[nt%2]),
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

def PlotCathNorm(basepath):
	CF = CathFile(basepath+"/MWPCCal.txt")
	caths = CF.hitdists.keys()
	caths.sort()
	cathnorms = {}
	
	gCNorm=graph.graphxy(width=15,height=10,
						   x=graph.axis.lin(title="cathode number",min=0,max=15,
											parter=graph.axis.parter.linear(tickdists=[5,1])),
						   y=graph.axis.lin(title="normalization factor",min=0.8,max=1.2),
						   key = graph.key.key(pos="bc",columns=2))
	setTexrunner(gCNorm)
	
	gFill=graph.graphxy(width=15,height=10,
						 x=graph.axis.lin(title="cathode number",min=0,max=15,
										  parter=graph.axis.parter.linear(tickdists=[5,1])),
						 y=graph.axis.lin(title="fill fraction",min=0.75,max=1.25),
						 key = graph.key.key(pos="bc",columns=2))
	setTexrunner(gFill)
	
	lstyles = {"x":[],"y":[style.linestyle.dashed]}
	lcols = {"E":[rgb.red],"W":[rgb.blue]}
	ssymbs = {"E":symbol.circle,"W":symbol.triangle}
	
	for s in ["E","W"]:
		for d in ["x","y"]:
			print "-----------",s,d,"------------"
			gdat = [ CF.cathsegs[k] for k in caths if k[0]==s and k[1]==d ]
			gdat = [ (g.i,g.renorm()*CF.cnorms[(s,d)][g.i],g.fill_frac,g) for g in gdat ]
			gdat.sort()
			gdat = gdat[1:-1]
			gnorm = sum([g[1] for g in gdat])/len(gdat)
			for g in gdat:
				print s,d,g[0],g[1],g[-1].fill_frac,g[-1].renorm(),g[-1].width
			gdat = [ (g[0],g[1]/gnorm,g[2]) for g in gdat ]
			gtitle = "%s%s"%(s,d)
			gCNorm.plot(graph.data.points(gdat,x=1,y=2,title=gtitle),
			  [graph.style.line(lineattrs=lstyles[d]+lcols[s]),
			   graph.style.symbol(ssymbs[s],symbolattrs=lcols[s])])
			gFill.plot(graph.data.points(gdat,x=1,y=3,title=gtitle),
						[graph.style.line(lineattrs=lstyles[d]+lcols[s]),
						 graph.style.symbol(ssymbs[s],symbolattrs=lcols[s])])
			
			for g in gdat:
				cathnorms[(s,d,g[0])] = g[1]
			#cathnorms[(s,d,0)]=cathnorms[(s,d,1)]
			#cathnorms[(s,d,15)]=cathnorms[(s,d,14)]
			cathnorms[(s,d,0)] = 0.85
			cathnorms[(s,d,15)] = 0.85
	
	cathnorms[("W","x",1)] = 0.9
	cathnorms[("E","x",14)] = 0.9
				
	gCNorm.writetofile(basepath+"/CathNorm.pdf")
	gFill.writetofile(basepath+"/CathFill.pdf")
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
		for s in ["East","West"]:
			for tp in range(2):
				gtitle = s + " Type " + {0:"0",1:"I",2:"II/III"}[tp]
				gAGC.plot(graph.data.points(gdat[(s[0],tp)],x=1,y=-1,title=gtitle),[graph.style.symbol(tpsymbs[tp],symbolattrs=[scols[s],])])
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
		
		g23=graph.graphxy(width=15,height=10,
						   x=graph.axis.lin(title="Energy [keV]"),
						   y=graph.axis.lin(title="Fraction of II+III",min=0,max=1),
						   key = graph.key.key(pos="tc",columns=2))
		setTexrunner(g23)

		print gdat
		
		g23.plot(graph.data.points(gdat,x=-1,y=4,title="II"),[graph.style.line(lineattrs=[rgb.red,])])
		g23.plot(graph.data.points(gdat,x=-1,y=5,title="III*"),[graph.style.line(lineattrs=[rgb.blue,style.linestyle.dashed])])
		g23.plot(graph.data.points(gdat,x=-1,y=6,title="III"),[graph.style.line(lineattrs=[rgb.blue,])])
		g23.plot(graph.data.points(gdat,x=-1,y=7,title="II*"),[graph.style.line(lineattrs=[rgb.red,style.linestyle.dashed])])
		
		g23.writetofile(basedir+"frac23_%s.pdf"%s)

###############
#             #
###############

if __name__ == "__main__":
	
	sep23effic(os.environ["UCNA_ANA_PLOTS"]+"/test/23Separation_SimPen/")
	
	#conn = open_connection()
	#anodeGainCal(conn)
	#gen_cathcal_set(conn,13000,100000,14264,16077)
