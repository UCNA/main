#!/usr/bin/python

from EncalDB import *
from math import *
from QFile import *
from PyxUtils import *
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
		self.loadFloats(["center","d_center","height","d_height","width","d_width","i","position"])
		self.i = int(self.i)
		self.loadStrings(["side","plane"])

class hitdist(KVMap):
	def __init__(self,m=KVMap()):
		self.dat = m.dat
		self.loadFloats(["cathode","energy"])
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
					("East","X"): ( 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231 ),
					("East","Y"): ( 20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  210, 211, 212, 213, 214, 215 ),
					("West","X"): ( 31,  30,  29,  28,  27,  26,  25,  24,  23,  22,  21,  20,  19,  18,  17,  16 ),
					("West","Y"): ( 0,   1,   2,   3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  13,  14,  15 ) }
	wires = []
	for i in range(nWires):
		pos = (nWires*0.5-0.5-i)*wireSpacing
		nm = "Pdc%i"%padc_nums[(s,d)][i]
		if s=="West":
			nm = "Padc%i"%padc_nums[(s,d)][i]
		wires.append(wireSpec(pos,nm,i))
	return wires

def gen_cathcal_set(conn,r0,r1,cr0,cr1):
	"""Plain cathcal set with no shape corrections"""
	
	basepath = os.environ["UCNA_ANA_PLOTS"]+"/WirechamberCal/%i-%i"%(cr0,cr1)
	corrcurves = PlotCathShapes(basepath)
	cathnorms = PlotCathNorm(basepath)
	
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
		gdat = [ [h.energy,]+h.terms[1:]+h.dterms[1:] for h in CF.hitdists[k]]
		gdat.sort()
		nterms = (len(gdat[0])-1)/2
		tcols = rainbow(nterms)
						  
		gHitDist=graph.graphxy(width=15,height=10,
						x=graph.axis.lin(title="Scintillator Visible Energy [keV]",min=0,max=800),
						y=graph.axis.lin(title="%s%s-%i Distortion Coeffiecient"%k,min=-0.2,max=0.2),
						key = graph.key.key(pos="bc",columns=3))
		setTexrunner(gHitDist)
		
		corrcurves[k] = []
		for n in range(nterms):
			tdat = [(g[0],0,g[n+1],g[nterms+n+1]) for g in gdat]
			corrcurves[k].append(tdat)
			gtitle = "$\\sin(%i\\pi x)$"%(2*(n/2+1))
			if n%2:
				gtitle = "$\\cos(%i\\pi x)$"%(2*(n/2+1))
			gHitDist.plot(graph.data.points(tdat,x=1,y=3,dy=4,title=gtitle),
					 [graph.style.line(lineattrs=lstyles[n/2]+lwidths[n%2]+lcols[n%2]),
					  graph.style.errorbar(errorbarattrs=lcols[n%2])])
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
								   x=graph.axis.lin(title="Scintillator Visible Energy [keV]",min=0,max=800),
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
					avgcurves[(s,d)][-1].append([e,0,mu,sg])
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
						   y=graph.axis.lin(title="max normalized signal",min=0,max=4),
						   key = graph.key.key(pos="bc",columns=2))
	setTexrunner(gCNorm)
	
	lstyles = {"x":[],"y":[style.linestyle.dashed]}
	lcols = {"E":[rgb.red],"W":[rgb.blue]}
	ssymbs = {"E":symbol.circle,"W":symbol.triangle}
	
	for s in ["E","W"]:
		for d in ["x","y"]:
			gdat = [ CF.cathsegs[k] for k in caths if k[0]==s and k[1]==d ]
			gdat = [ (g.i,g.height,g.d_height) for g in gdat ]
			gdat.sort()
			gdat = gdat[1:-1]
			gtitle = "%s%s"%(s,d)
			gCNorm.plot(graph.data.points(gdat,x=1,y=2,dy=3,title=gtitle),
			  [graph.style.line(lineattrs=lstyles[d]+lcols[s]),
			   graph.style.symbol(ssymbs[s],symbolattrs=lcols[s]),
			   graph.style.errorbar(errorbarattrs=lcols[s])])
			
			gnorm = sum([g[1] for g in gdat])/len(gdat)
			for g in gdat:
				cathnorms[(s,d,g[0])] = gnorm/g[1]
			cathnorms[(s,d,0)]=cathnorms[(s,d,1)]
			cathnorms[(s,d,15)]=cathnorms[(s,d,14)]
						
	gCNorm.writetofile(basepath+"/CathNorm.pdf")
	return cathnorms
		



###############
#             #
###############

if __name__ == "__main__":
	conn = open_connection()
	gen_cathcal_set(conn,13000,100000,14282,14347)
