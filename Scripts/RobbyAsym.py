#!/sw/bin/python2.6

from Asymmetries import *

def Robby_asymmetryFit(l):
	"""Parse an asymmetry line from Robby's output file"""
	l = l.split("=")[1].split("+/-")
	a = asymmetryFit()
	a.A0 = float(l[0])
	a.dA0 = float(l[1])
	a.fitMin=200
	a.fitMax=675
	return a

def load_rfile(fname):
	"""Load asymmetries for each level of octet division from Robby's file"""
	fin = open(fname,"r")
	asyms = [[asymData()],[asymData(),asymData()],[asymData(),asymData(),asymData(),asymData()]]
	
	# load runs list
	fin.readline()
	fin.readline()
	l = fin.readline()
	while l[0] != '=':
		l = l.split()
		rn = int(l[5])
		if rn:
			for i in range(3):
				for j in range(len(divSegments[i])):
					if l[0] in divSegments[i][j]:
						asyms[i][j].addRun(l[0],rn)
		l = fin.readline()
		l = fin.readline()
	
	# load asymmetries
	fin.readline()
	for i in range(4):
		asyms[2][i].asyms.append(Robby_asymmetryFit(fin.readline()))
	
	fin.readline()
	fin.readline()
	for i in range(2):
		fin.readline()
		asyms[1][i].asyms.append(Robby_asymmetryFit(fin.readline()))
		
	fin.readline()
	fin.readline()
	fin.readline()
	asyms[0][0].asyms.append(Robby_asymmetryFit(fin.readline()))
																
	return asyms
	
def collectRobbyAsymmetries(basedir,depth):
	print "Loading Robby's asymmetries for depth",depth
	asyms = []
	for f in [ f for f in os.listdir(basedir) if f[-4:]==".txt"]:
		rn = int(f.split("_")[-1][:-4])
		if rn > 13000:
			print rn
			asyms += load_rfile(basedir+"/"+f)[depth]
	return asyms

def pairAsyms(as1,as2):
	"""Pair up asymmetries for the same run ranges"""
	
	print "Pairing",len(as1),"and",len(as2),"asymmetries."
	
	rd1 = {}
	for a in as1:
		rd1[tuple(a.getRuns())] = a
	rd2 = {}
	for a in as2:
		rd2[tuple(a.getRuns())] = a
	
	pairs = {}
	rkeys = rd1.keys()
	rkeys.sort()
	for k in rkeys:
		if k not in rd2:
			print k,"in 1 not 2"
		else:
			pairs[k] = (rd1[k],rd2[k])
	
	print "------------------"
	rkeys = rd2.keys()
	rkeys.sort()
	for k in rkeys:
		if k not in rd1:
			print k,"in 2 not 1"

	print "------------------\nFound",len(pairs),"pairs."
	
	return pairs
	

def compareMichaelRobbyAsyms(depth):
	rasyms = collectRobbyAsymmetries("/Users/michael/Desktop/octet_debugging/",depth)
	masyms = collectAsymmetries("../../PostPlots/OctetAsym_Offic/",depth)
	pairs = pairAsyms(rasyms,masyms)

	gdat = [ (p,pairs[p][0].getAsym(200,675),pairs[p][1].getAsym(200,675)) for p in pairs]
	gdat = [ (g[0][0],-abs(g[1].A0),g[1].dA0,g[2].A0,g[2].dA0) for g in gdat if g[1] and g[2]]

	# Robby's A
	LFr = LinearFitter(terms=[polyterm(0)])
	LFr.fit(gdat,cols=(0,1,2),errorbarWeights=True)
	chi2_r = LFr.ssResids()
	nu_r = len(gdat)-len(LFr.coeffs)
	dA_r = 1.0/sqrt(LFr.sumWeights())
	# Michael's A
	LFm = LinearFitter(terms=[polyterm(0)])
	LFm.fit(gdat,cols=(0,3,4),errorbarWeights=True)
	chi2_m = LFm.ssResids()
	nu_m = len(gdat)-len(LFm.coeffs)
	dA_m = 1.0/sqrt(LFm.sumWeights())

	#stats.chisqprob(chi2,ndf);

	gAsyms=graph.graphxy(width=15,height=15,
				x=graph.axis.lin(title="Robby's $A=%.5f(%i)$, $\\chi^2/\\nu=%.1f/%i$"%(LFr.coeffs[0],int(100000*dA_r),chi2_r,nu_r),min=-0.06,max=-0.035),
				y=graph.axis.lin(title="Michael's $A=%.5f(%i)$, $\\chi^2/\\nu=%.1f/%i$"%(LFm.coeffs[0],int(100000*dA_m),chi2_m,nu_m),min=-0.06,max=-0.035),
				key = graph.key.key(pos="bl"))
	gAsyms.texrunner.set(lfs='foils17pt')
	
		
	gAsyms.plot(graph.data.points(gdat,x=2,dx=3,y=4,dy=5,title=None),
				[graph.style.errorbar(errorbarattrs=[rgb.blue,]),
				graph.style.symbol(symbol.circle,size=0.2,symbolattrs=[rgb.red,]),
				])
	gAsyms.plot(graph.data.function("y(x)=x",title=None), [graph.style.line(lineattrs=[style.linestyle.dashed,]),])
	
	gAsyms.writetofile("/Users/michael/Desktop/Robby_v_Michael_%i.pdf"%depth)
				
if __name__=="__main__":
	for i in range(3):
		compareMichaelRobbyAsyms(i)
	
	