#!/usr/bin/python

from Sources import *
	
class pmtCorr(KVMap):
	def __init__(self,m):
		self.dat = m.dat
		self.loadFloats(["E0","E1","corr","dcorr","t1","t2","sID"])
		self.t1 = int(self.t1)
		self.t2 = int(self.t2)
		self.sID = int(self.sID)
		self.loadStrings(["name","side","simulated"])
		self.simulated = (self.simulated == "yes")
							
def plotCorrelations(conn,rlist,outdir):

	# gather data
	rlist.sort()
	SDD = SourceDatDirectory()
	corrdat = {}
	for rn in rlist:
		for c in [pmtCorr(c) for c in SDD.getKey(rn,"correlation")]:
			c.run = rn
			corrdat.setdefault((c.E0,c.E1),{}).setdefault((c.t1,c.t2),[]).append(c)

	# plot
	tcols = rainbow(4)
	for kr in corrdat.keys():
		for s in ['E','W']:
		
			# set up graph
			tckdist = [5,1]
			if rlist[-1]-rlist[0] > 100:
				tckdist = [10,1]
			runaxis = graph.axis.lin(title="Run Number",min=rlist[0]-5,max=rlist[-1]+1,
								parter=graph.axis.parter.linear(tickdists=tckdist),
								texter = graph.axis.texter.rational(),
								painter=graph.axis.painter.regular(labeldist=0.1,labeldirection=graph.axis.painter.rotatetext(135)))
			gRuns=graph.graphxy(width=30,height=15,
				x2=runaxis,
				y=graph.axis.lin(title="PMT Correlation Slope"),
				key = graph.key.key(pos="tl"))
			try:
				gRuns.texrunner.set(lfs='foils17pt')
			except:
				pass
			
			for t1 in range(4):
				for t2 in range(4):
					if t1 == t2:
						continue
					gsim = [ (c.run,c.corr,c.dcorr) for c in corrdat[kr][(t1,t2)] if c.side == s and c.simulated ]	
					gsim.sort()
					gRuns.plot(graph.data.points(gsim,x=1,y=2,dy=3,title=None),
								[graph.style.line([tcols[t1],style.linewidth.THick]),
								 graph.style.line([tcols[t2],style.linewidth.thick,style.linestyle.dashed]),
								 graph.style.errorbar(errorbarattrs=[tcols[t1]])])

			for t1 in range(4):
				for t2 in range(4):
					if t1 == t2:
						continue
					gdat = [ (c.run,c.corr,c.dcorr) for c in corrdat[kr][(t1,t2)] if c.side == s and not c.simulated ]
					gdat.sort()
					gRuns.plot(graph.data.points(gdat,x=1,y=2,dy=3,title="%i v %i"%(t1,t2)),
								[graph.style.symbol(symbol.circle,symbolattrs=[deco.filled,tcols[t1],deco.stroked([style.linewidth.THick,tcols[t2],])]),
								 graph.style.line([tcols[t1],style.linewidth.thick]),
								 graph.style.errorbar(errorbarattrs=[tcols[t1]])])
					

			os.system("mkdir -p "+outdir)
			gname = outdir+"/%i%s_%i-%i.pdf"%(rlist[0],s,kr[0],kr[1])
			print gname
			gRuns.writetofile(gname)


if __name__=="__main__":
	
	outpath = os.environ["UCNA_ANA_PLOTS"]+"/Sources/"
	conn = open_connection()
	
	for c in cal_2010:
	
		rlist = range(c[0],c[1]+1)
		plotCorrelations(conn,rlist,outpath+"/PMTCorr")
		