#!/usr/bin/python


import sys
sys.path.append("..")
import os
from calib.Sources import *
from ucnacore.Histogram import *

limdat = {2008:[(0,5.0),(250,5.0),(500,500*0.013),(900,900*0.025),(1000,1000*0.025),(1200,1200*0.025)],
			2010:[(0,2.5),(200,200*0.0125),(500,500*0.0125),(1000,500*0.0125)],
			2011:[(0,0.02*180.3),(180.3,180.3*0.02),(418.4938,0.018*418.4938),(1200,1200.*0.012)], #Parallel Analyzer
	                20112:[(0,6.),(200,6.),(500,10.),(1000,10.)], #Conservative envelope
	                #2011:[(0,200.*0.02),(200,200*0.02),(500,500*0.0150),(1000,500*0.0150)],#liberal envelope
			2012:[(0,0.02*180.3),(180.3,180.3*0.02),(418.4938,0.018*418.4938),(1200,1200.*0.011)]} #Parallel Analyzer 

badruns = {2011:[(17359,17360),(17361,17361),(17376,17376),(17520,17520),(17874,17874),(17891,17891),(18037,18037),(18370,18386),(18620,18620),(19357,19357),(17233,17249)],#,(19823,19863)],#,(18617,19863)],
	   #2012:[(20818,20818),(20829,20829),(21094,21094),(21322,21322),(22451,22451),(22782,22782)] }

#The badruns below are to mimic precisely the runs which are not being used in the parallel analyzer. For 
# consistency with what Dave has said to ignore, use above list.
#badruns = {2011:[(17383,17383),(17385,17385),(17382,17382),(17886,17886),(17912,17912),(19347,19364),(17233,17249),(17871,18055),(19899,20000)],#,(19823,19863)],#,(18617,19863)],
	   2012:[(20818,20818),(20829,20829),(21094,21094),(21322,21322),(22451,22451),(22782,22782)] }
#This is for MB personal test of running energy cal
#badruns = {2011:[(17233,17249)],2012:[]}
badIn114Runs = {2011:[],#[(19823,19863)],
		2012:[]}
RmSrcTypeList = ["Cd109","Cd113m","In114","In114E","In114W","Cs137"]

# Uncertainties for linear fit
#gausUnc = {2011:[(4.089,5.516),(-1.851,4.756),(-1.026,3.253),(-0.284,1.182),(-1.639,1.921),(-2.298,1.861)],
#	   2012:[]}
 # Uncertainties for quadratic fits
gausUnc = {2011:[(-2.575,7.265),(-2.187,6.367),(-0.8229,5.208),(-0.9457,1.419),(-2.894,2.386),(-1.276,3.772)],
	   2012:[]}

def calEnvelope(E,year=2011):	
	i = 0
	while E > limdat[year][i+1][0]:
		i+=1
	l = (E-limdat[year][i][0])/(limdat[year][i+1][0]-limdat[year][i][0])
	return (1-l)*limdat[year][i][1]+l*limdat[year][i+1][1]

def in_cal_list(rn,clist,year):
	for c in clist:
		if c[0] <= rn <= c[1]:
			for pair in badruns[year]: 
				for r in range(pair[0],pair[1]+1,1):
					if rn==r:
						return False
			 
			return True
	return False;

def calrun_ranges(year):
	if year==2010:
		return range(14000,14746+1)+range(15645,15939+1)
		#rlist = range(13883,14746+1)+range(15645,15939+1)
	if year==2011:
		return [r for r in range(17233,20000) if in_cal_list(r,cal_2011,year) ]
	if year==2012:
		return [r for r in range(20000,22300) if in_cal_list(r,cal_2012,year) ]

def plot_Cal_Uncertainty(g,title=None,st=[graph.style.line([style.linestyle.dotted])],year=2010):
	"""Plot energy uncertainty envelope."""
	g.plot(graph.data.points(limdat[year],x=1,y=2,title=title),st)
	g.plot(graph.data.points([ (x[0],-x[1]) for x in limdat[year]],x=1,y=2,title=None),st)

srcYearDat = {}
def getSourceYearDat(year):
	if year not in srcYearDat:
		conn = open_connection()
		rlist = calrun_ranges(year)
		SDC= SourceDataCollector(conn)
		SDC.gather_peakdat(rlist, "AND peak_num<100")
		srcYearDat[year] = SDC;
	return srcYearDat[year]

# plot linearity calibration source errors
def plotAllErrors(outpath,year,s="Both",t=4,gaussFiles=False,gaussUnc=False):
	
	# gather source data from calibration runs
	srs = sort_by_type( getSourceYearDat(year).slines )
	
	scols = rainbowDict(srs)
	if "PUBLICATION_PLOTS" in os.environ:
			for k in scols:
				scols[k] = rgb.black


	yrange = 30
	if t != 4 or year > 2010:
		yrange = 30
		
	g=graph.graphxy(width=30,height=10,
			#x=graph.axis.lin(title="Expected E$_{\\rm vis}$ [keV]",min=0,max=1500),
			x=graph.axis.lin(title="E$_{\\rm recon}$ [keV]",min=0,max=1500),
			y=graph.axis.lin(title="reconstructed energy error [keV]",min=-yrange,max=yrange),
			key = graph.key.key(pos="br"))
	setTexrunner(g)
	
	# plot error envelope
	if t == 4:
		plot_Cal_Uncertainty(g,"uncertainty envelope",[graph.style.line([style.linestyle.dotted,style.linewidth.THick])],year=year)
	else:
		g.plot(graph.data.function("y(x)=0",title=None), [graph.style.line([style.linestyle.dotted,style.linewidth.THick])])

	gausInc=0
	#making array of bad indium runs for whatever year is chosen
	In114List = []
	for pair in badIn114Runs[year]:
		for r in range(pair[0],pair[1]+1,1):
			In114List.append(r)

	# plot
	for k in srs:
		
		gdat = [(l.sim.erecon+50,l.erecon-l.sim.erecon,l) for l in srs[k] if l.tube==t and l.src.radius()>0. and l.src.radius()<45. and (s=="Both" or l.side==s) and l.src.run not in In114List and l.src.type not in RmSrcTypeList]

		gdat = [p for p in gdat if abs(p[1])<yrange]
		if not gdat:
			continue
		if gaussFiles:
			ofile = open("%s/%s_%i.txt"%(os.environ["UCNA_ANA_PLOTS"],peakNames[k][8:],year),'w')
		print k,peakNames[k]
		for l in gdat:
			if abs(l[-1].erecon-l[-1].sim.erecon) > 15:
				print "****************"
			#if abs(l[-1].erecon-l[-1].sim.erecon) > 6:
				#print "Large error",l[-1].src.run,l[-1].uid,l[-1].erecon-l[-1].sim.erecon,"from expected"
		x0 = gdat[0][0]
		hErr = histogram(int(100/(5*(x0+200)/1000)),-yrange-1,yrange+1)
		for l in gdat:
			hErr.fill(l[1],6)
			if gaussFiles:	
				ofile.write('%f\n'%l[1])
		if gaussFiles:	
			ofile.close
		if not gaussUnc:
			gtitle = "%s: $%.1f\\pm%.1f$keV"%(peakNames.get(k,k),hErr.avg(),hErr.rms())
		else:
			gtitle = "%s: $%.1f\\pm%.1f$keV"%(peakNames.get(k,k),gausUnc[year][gausInc][0],gausUnc[year][gausInc][1])
		#g.plot(graph.data.points(gdat,x=1,y=2,title=gtitle), [graph.style.symbol(peakSymbs.get(k,symbol.circle),size=0.2,symbolattrs=[scols[k]]),])
		g.plot(graph.data.points(hErr.lineData(yoff=x0),x=2,y=1,title=None), [graph.style.line([scols[k],style.linewidth.THick]),])
		g.plot(graph.data.points([[x0,-30],[x0,30]],x=1,y=2,title=None), [graph.style.line([scols[k],style.linestyle.dotted]),])
		if not gaussUnc:
			g.plot(graph.data.points([[x0,hErr.avg(),hErr.rms()]],x=1,y=2,dy=3,title=gtitle),
			[
				graph.style.errorbar(errorbarattrs=[style.linewidth.THick,scols[k]]),
				graph.style.symbol(peakSymbs.get(k,symbol.circle),size=0.3,symbolattrs=[scols[k],deco.filled([rgb.white]),style.linewidth.THick]),
			])
		else:
			g.plot(graph.data.points([[x0,gausUnc[year][gausInc][0],gausUnc[year][gausInc][1]]],x=1,y=2,dy=3,title=gtitle),
			[
				graph.style.errorbar(errorbarattrs=[style.linewidth.THick,scols[k]]),
				graph.style.symbol(peakSymbs.get(k,symbol.circle),size=0.3,symbolattrs=[scols[k],deco.filled([rgb.white]),style.linewidth.THick]),
			])
		gausInc+=1
	
	if s != "Both":
		g.writetofile(outpath+"/CalErrors_%i_%s_%i.pdf"%(year,s,t))
	else:
		g.writetofile(outpath+"/CalErrors_%i.pdf"%year)
	
	
def plotAllWidths(outpath,year,s="Both",t=4):

	# gather source data from calibration runs
	srs = getSourceYearDat(year)
	scols = rainbowDict(srs)
	if "PUBLICATION_PLOTS" in os.environ:
			for k in scols:
				scols[k] = rgb.black

	# set up graph
	yrange = 20
	xmin,xmax = 25,175
	gwid = 15
	if t==4:
		yrange = 10
		xmin,xmax = 20,80
		gwid = 30
	g=graph.graphxy(width=gwid,height=8,
			x=graph.axis.lin(title="expected peak width [keV]",min=xmin,max=xmax),
			y=graph.axis.lin(title="observed $-$ expected width [keV]",min=-yrange,max=yrange),
			key = graph.key.key(pos="bc",columns=2))
	setTexrunner(g)

	g.plot(graph.data.function("y(x)=0",title=None), [graph.style.line([style.linestyle.dotted,style.linewidth.THick])])

	# plot data for each source
	for k in srs.keys()[::-1]:
		gdat = [(l.sim.enwidth,l.enwidth-l.sim.enwidth,sqrt(l.denwidth**2+l.sim.denwidth**2)) for l in srs[k] if l.tube==t and l.src.radius()<50. and (s=="Both" or l.side==s)]
		if not gdat:
			continue
		mu,sigma = musigma([x[1] for x in gdat if abs(x[1])<yrange])
		gtitle = "%s: $%.1f\\pm%.1f$keV"%(peakNames.get(k,k),mu,sigma)
		print gtitle
		g.plot(graph.data.points(gdat,x=1,y=2,dy=3,title=gtitle), [ graph.style.symbol(peakSymbs.get(k,symbol.circle),size=0.2,symbolattrs=[scols[k]])] )

	if s != "Both":
		g.writetofile(outpath+"/WidthErrors_%i_%s_%i.pdf"%(year,s,t))
		#g.writetofile(outpath+"/%s_%i.pdf"%(s,t))
	else:
		g.writetofile(outpath+"/WidthErrors_%i.pdf"%year)


	
if __name__=="__main__":

	# set up output paths
	#outpath = os.environ["UCNA_ANA_PLOTS"]+"/Sources/ErrorEnvelope/FitTypeStudy/"
	outpath = os.environ["UCNA_ANA_PLOTS"]+"/Sources/ErrorEnvelope/MBTest/"
	os.system("mkdir -p %s"%outpath)
	
	if 0:
		plotAllWidths(outpath,2010)
		plotAllErrors(outpath,2010)
		exit(0)
		
		for s in ["East","West"]:
			for t in range(5):
				plotAllErrors(outpath,2010,s,t)
				#plotAllWidths(outpath,2010,s,t)

	plotAllErrors(outpath,2011,gaussFiles=False,gaussUnc=False)
	#for s in ["East","West"]:
	#	plotAllErrors(outpath,2011,s,gaussFiles=False,gaussUnc=False)
	#for s in ["East","West"]:
	#		for t in range(4):
	#			plotAllErrors(outpath,2011,s,t,gaussFiles=False,gaussUnc=False)
				#plotAllWidths(outpath,2010,s,t)
		
