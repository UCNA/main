#!/sw/bin/python2.7

import sys
sys.path.append("..")

import os
from review.Asymmetries import *
from ucnacore.PyxUtils import *
from ucnacore.LinFitter import *

def collectRecursive(basedir,depth):
	print "Collecting asymmetry data from",basedir,"at depth",depth
	asyms = []
	fls = os.listdir(basedir)
	fls.sort()
	for f in [ (basedir+'/'+f+'/',basedir+'/'+f+'/'+f+'.txt') for f in fls if os.path.isdir(basedir+'/'+f) and '_' in f]:
		if not depth and os.path.exists(f[1]):
			asyms.append(AsymmetryFile(f[1]))
		else:
			asyms += collectRecursive(f[0],depth-1)
	return asyms

def fmt_llog(coeffs,xpt=4):
	a = exp(coeffs[0])*10**(xpt*coeffs[1])
	return r"$%.2f \cdot (N/10^%i)^{%.2f}$"%(a,xpt,coeffs[1])

if __name__=="__main__":

	outdir = os.environ["UCNA_ANA_PLOTS"]+"/test/SimStats/"
	os.system("mkdir -p %s"%outdir)

	gK=graph.graphxy(width=15,height=15,
				  x=graph.axis.log(title="Simulated events",min=1e3,max=2e5),
				  y=graph.axis.log(title="RMS endpoint spread [keV]",min=0.5,max=50),
				  key = graph.key.key(pos="tr"))
	setTexrunner(gK)
	
	gA=graph.graphxy(width=15,height=15,
				  x=graph.axis.log(title="Simulated events",min=1e3,max=2e5),
				  y=graph.axis.log(title="RMS asymmetry spread [\\% of $A$]",min=.05,max=5),
				  key = graph.key.key(pos="tr"))
	setTexrunner(gA)
	
	scols = [rgb.blue,rgb.red]
	snames = ["pseudo-random","quasi-random"]
	
	for (nstatf,statf) in enumerate(["SimStats_PR","SimStats_QR"]):
	
		gdat = []
		for d in [1,2,3,4,5,6,7,8]:
			fls = collectRecursive(os.environ["UCNA_ANA_PLOTS"]+"/test/%s/"%statf,d)
			asyms = []
			kuries = []
			for af in fls:
				kuries.append(af.getKurie("E","On","0").ep)
				asyms.append(af.getAsym().A0)
			n = len(fls)
			sscl = 1./sqrt(2.*(n-1))
			k_m,k_s = musigma(kuries)
			a_m,a_s = musigma(asyms)
			a_rel = 100.*a_s/abs(a_m)
			gdat.append([d,1000*(2**(8-d)),len(kuries),k_s,a_rel,k_s*sscl,a_rel*sscl])

		print gdat

		LL = LogLogger(terms=[polyterm(0),polyterm(1)])
		LL.fit(gdat,(1,3,5),errorbarWeights=True)
		
		gK.plot(graph.data.points(gdat,x=2,y=4,dy=6,title=snames[nstatf]+": "+fmt_llog(LL.coeffs)),
				[graph.style.symbol(symbol.triangle,size=0.2,symbolattrs=[scols[nstatf],]), graph.style.errorbar(errorbarattrs=[scols[nstatf]])])

		gK.plot(graph.data.points(LL.fitcurve(1e3,2e5,100),x=1,y=2,title=None),
				[graph.style.line(lineattrs=[scols[nstatf],])])
				

		LL = LogLogger(terms=[polyterm(0),polyterm(1)])
		LL.fit(gdat,(1,4,6),errorbarWeights=True)
		
		gA.plot(graph.data.points(gdat,x=2,y=5,dy=7,title=snames[nstatf]+": "+fmt_llog(LL.coeffs)),
				[graph.style.symbol(symbol.triangle,size=0.2,symbolattrs=[scols[nstatf],]), graph.style.errorbar(errorbarattrs=[scols[nstatf]])])

		gA.plot(graph.data.points(LL.fitcurve(1e3,2e5,100),x=1,y=2,title=None),
				[graph.style.line(lineattrs=[scols[nstatf],])])


	gK.writetofile(outdir+"/KurieSpread.pdf")
	gA.writetofile(outdir+"/AsymSpread.pdf")