#!/usr/bin/python

from Sources import *
from Histogram import *

limdat = {2008:[(0,5.0),(250,5.0),(500,500*0.013),(900,900*0.025),(1000,1000*0.025),(1200,1200*0.025)],
			2010:[(0,2.5),(200,200*0.0125),(500,500*0.0125),(1000,500*0.0125)],
			2011:[(0,2.5),(200,200*0.0125),(500,500*0.0125),(1000,500*0.0125)] }

def calEnvelope(E,year=2010):	
	i = 0
	while E > limdat[year][i+1][0]:
		i+=1
	l = (E-limdat[year][i][0])/(limdat[year][i+1][0]-limdat[year][i][0])
	return (1-l)*limdat[year][i][1]+l*limdat[year][i+1][1]

def plot_Cal_Uncertainty(g,title=None,st=[graph.style.line([style.linestyle.dashed])],year=2010):
	"""Plot energy uncertainty envelope."""
	g.plot(graph.data.points(limdat[year],x=1,y=2,title=title),st)
	g.plot(graph.data.points([ (x[0],-x[1]) for x in limdat[year]],x=1,y=2,title=None),st)
	
# plot linearity calibration source errors
def plotAllErrors(outpath,year,s="Both",t=4):
	
	yrange = 15

	g=graph.graphxy(width=45,height=15,
			x=graph.axis.lin(title="Expected E$_{vis}$ [keV]",min=0,max=1300),
			y=graph.axis.lin(title="Reconstructed E$_{vis}$ Error [keV]",min=-yrange,max=yrange),
			key = graph.key.key(pos="tr"))
	if s=="Both":
		g.texrunner.set(lfs='foils17pt')
	
	# plot error envelope
	plot_Cal_Uncertainty(g,"Provisional Envelope 2010",[graph.style.line([style.linestyle.dashed,style.linewidth.THick])],year=year)

	# gather source data from calibration runs
	conn = open_connection()
	rlist = []
	if year==2010:
		rlist = range(14000,14746+1)+range(15645,15939+1)
		#rlist = range(13883,14746+1)+range(15645,15939+1)
	if year==2011:
		rlist = range(17233,100000)
	slines = gather_peakdat(conn,rlist)
	srs = sort_by_type(slines)
	scols = rainbowDict(srs)
	
	# plot
	for k in srs:
		gdat = [(l.sim.erecon,l.erecon-l.sim.erecon,l) for l in srs[k] if l.tube==t and l.src.radius()<50. and (s=="Both" or l.side==s)]
		if not gdat:
			continue
		print k,peakNames[k]
		for l in gdat:
			if abs(l[-1].erecon-l[-1].sim.erecon) > 15:
				print "****************"
			if abs(l[-1].erecon-l[-1].sim.erecon) > 6:
				print "Large error",l[-1].src.run,l[-1].uid,l[-1].erecon-l[-1].sim.erecon,"from expected"
		g.plot(graph.data.points(gdat,x=1,y=2,title=None), [graph.style.symbol(symbol.circle,size=0.2,symbolattrs=[scols[k]]),])
		x0 = gdat[0][0]
		hErr = histogram(int(100/(5*(x0+200)/1000)),-yrange-1,yrange+1)
		for l in gdat:
			hErr.fill(l[1],10)
		g.plot(graph.data.points(hErr.lineData(yoff=x0),x=2,y=1,title="%s: $%.1f\\pm%.1f$keV"%(peakNames[k],hErr.avg(),hErr.rms())), [graph.style.line([scols[k],style.linewidth.THIck]),])
		g.plot(graph.data.points([[x0,-30],[x0,30]],x=1,y=2,title=None), [graph.style.line([scols[k],style.linestyle.dashed]),])
	
	if s != "Both":
		g.writetofile(outpath+"/CalErrors_%i_%s_%i.pdf"%(year,s,t))
	else:
		g.writetofile(outpath+"/CalErrors_%i.pdf"%year)
	
	
	
	
if __name__=="__main__":

	# set up output paths
	outpath = os.environ["UCNA_ANA_PLOTS"]+"/Sources/ErrorEnvelope/"
	os.system("mkdir -p %s"%outpath)
		
	plotAllErrors(outpath,2010)
	for s in ["East","West"]:
		for t in range(5):
			plotAllErrors(outpath,2010,s,t)
			
	#plotAllErrors(outpath,2011)
	