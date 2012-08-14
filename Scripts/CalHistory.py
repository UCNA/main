#!/usr/bin/python

from LinFitter import *
from PyxUtils import *
from EncalDB import *
import os

def make_runaxis(rmin,rmax):
	
	tckdist = [5,1]
	if rmax-rmin > 100:
		tckdist = [10,1]
	if rmax-rmin > 500:
		tckdist = [100,10]	
	if rmax-rmin > 1000:
		tckdist = [100,20]
				
	return graph.axis.lin(title="Run Number",min=rmin,max=rmax,
							parter=graph.axis.parter.linear(tickdists=tckdist),
							texter = graph.axis.texter.rational(),
							painter=graph.axis.painter.regular(labeldist=0.1,labeldirection=graph.axis.painter.rotatetext(135)))


def plot_run_monitor(rlist,sname,tp="pedestal",outpath=None):

	rlist.sort()
	conn = open_connection()
	
	gdat = []
	for rn in rlist:
		if not rn%10:
			print rn
		cgwgid = getRunMonitorGIDs(conn,rn,sname,tp)
		if not cgwgid:
			print "*** Can't find run monitor for",rn,"***"
			continue	
		cgid,wgid = cgwgid
		cg = getGraph(conn,cgid)
		wg = getGraph(conn,wgid)
		if cg and wg:
			gdat.append([rn,cg[0][2],cg[0][3],wg[0][2],wg[0][3]])
	
	grmon=graph.graphxy(width=25,height=8,
				x=make_runaxis(rlist[0]-1,rlist[-1]+1),
				y=graph.axis.lin(title=sname),
				key = None)
	#grmon.texrunner.set(lfs='foils17pt')

	grmon.plot(graph.data.points(gdat,x=1,y=2,dy=4,title=None),
				[graph.style.errorbar(errorbarattrs=[rgb.blue,]),graph.style.symbol(symbol.circle,size=0.10,symbolattrs=[rgb.red,])])

	if outpath:
		grmon.writetofile(outpath+"/%s.pdf"%sname)
		
	return grmon
	
def plot_all_pedestals(rmin,rmax):

	rlist = getRunType(open_connection(),"Asymmetry",rmin,rmax)
	
	outpath = os.environ["UCNA_ANA_PLOTS"]+"/Pedestals/%i-%i/"%(rlist[0],rlist[-1])
	os.system("mkdir -p %s"%outpath)
	
	for s in ['E','W']:
		for t in range(4):
			sname = "ADC%s%iBeta"%(s,t+1)
			grmon = plot_run_monitor(rlist,sname,"pedestal",outpath)
		
		sname = "MWPC%sAnode"%s
		grmon = plot_run_monitor(rlist,sname,"pedestal",outpath)
				
		for p in ['x','y']:
			for c in range(16):
				sname = "MWPC%s%s%i"%(s,p,c+1)
				grmon = plot_run_monitor(rlist,sname,"pedestal",outpath)
	
			
def plot_trigeff_history(rmin,rmax):
	
	rlist = getRunType(open_connection(),"Asymmetry",rmin,rmax)
	outpath = os.environ["UCNA_ANA_PLOTS"]+"/Trigeff/%i-%i/"%(rlist[0],rlist[-1])
	os.system("mkdir -p %s"%outpath)
	conn = open_connection()
	
	for s in ['E','W']:
		for t in range(4):
			gdat = []
			for rn in rlist:
				tparms = list(getTrigeffParams(conn,rn,s,t))
				tparms.sort()
				try:
					gdat.append([rn,tparms[0][2],tparms[1][2],tparms[3][2],tparms[3][3]])
				except:
					print "***** Missing data for",rn,s,t
					continue
				if not rn%50:
					print rn
				
			
			gthresh=graph.graphxy(width=25,height=8,
						x=make_runaxis(rlist[0]-1,rlist[-1]+1),
						y=graph.axis.lin(title="Trigger ADC Threshold",min=0,max=100),
						key = None)
			gthresh.plot(graph.data.points(gdat,x=1,y=2,dy=3,title=None),
						[graph.style.errorbar(errorbarattrs=[rgb.blue,]),graph.style.symbol(symbol.circle,size=0.10,symbolattrs=[rgb.red,])])
			gthresh.writetofile(outpath+"/%s%i.pdf"%(s,t))
	
			geff=graph.graphxy(width=25,height=8,
						x=make_runaxis(rlist[0]-1,rlist[-1]+1),
						y=graph.axis.lin(title="Trigger Efficiency"),
						key = None)
			geff.plot(graph.data.points(gdat,x=1,y=4,dy=5,title=None),
						[graph.style.errorbar(errorbarattrs=[rgb.blue,]),graph.style.symbol(symbol.circle,size=0.10,symbolattrs=[rgb.red,])])
			geff.writetofile(outpath+"/Eff_%s%i.pdf"%(s,t))
	
if __name__ == "__main__":
	
	#plot_all_pedestals(14000,16300)
	#plot_all_pedestals(16500,19900)
	
	plot_trigeff_history(14000,16300)
	#plot_trigeff_history(16500,19900)