#!/usr/bin/python

from Asymmetries import *

def plot_bgSubtrHist(basedir,depth=2):
	
	print "--------------------- Division",depth,"-------------------------"
	bgs = {}
	n = 0
	conn = open_connection()
	for af in collectAsymmetries(basedir,depth):
		for k in af.bgs:
			bgs.setdefault(k,[]).append((getRunStartTime(conn,af.getRuns()[0]),af.bgs[k],af))
		n += 1
	
	gBgs=graph.graphxy(width=25,height=8,
					   #x=graph.axis.lin(title=unitNames[depth],min=0,max=n-1),
					   x=graph.axis.lin(title="Time [days]"),
					   y=graph.axis.lin(title="Residual Background [$\\mu$Hz/keV]"),
					   key = graph.key.key(pos="bl"))
	setTexrunner(gBgs)
	
	sideCols = {'E':rgb.red,'W':rgb.blue}
	segCols = {0:rgb.red,1:rgb.blue}
	segNames = {0:"BG Before",1:"BG After"}
	s = 'E'
	LF = LinearFitter(terms=[polyterm(0)])
	
	
	for seg in [0,1]:
		a0 = [b for b in bgs[(s,"Off","0")] if b[-1].whichSegment(depth)%2 == seg ] 
		a1 = [b for b in bgs[(s,"On","0")] if b[-1].whichSegment(depth)%2 == seg ]
		t0 = a0[0][0]
		
		gdat = [ (n,(b[0]-t0)/(24*3600.),(b[1].nrate+a1[n][1].nrate)*5e5,sqrt(b[1].d_nrate**2+a1[n][1].d_nrate**2)*5e5) for (n,b) in enumerate(a0)]
		LF.fit(gdat,cols=(0,2,3),errorbarWeights=True)
		err = 1.0/sqrt(LF.sumWeights())
		gBgs.plot(graph.data.points(gdat,x=2,y=3,dy=4,title="%s: $%.1f \\pm %.1f$"%(segNames[seg],LF.coeffs[0],err)),
				  [graph.style.symbol(symbol.circle,size=0.2,symbolattrs=[segCols[seg],]),
				   graph.style.errorbar(errorbarattrs=[segCols[seg],])])
	
	gBgs.writetofile(basedir+"/BGResid_%i.pdf"%depth)



def plot_TypeIV_resid(basedir,depth=2):
	print "--------------------- Division",depth,"-------------------------"
	
	rts = {}
	n = -1
	hname = "ExcessTheta"
	for af in collectAsymmetries(basedir,depth):
		n += 1
		if [badrun for badrun in [14166,14888,15518] if badrun in af.getRuns()]:
			continue
		for s in ["E","W"]:
			for afp in ["On","Off"]:
				rts.setdefault((s,afp),[]).append((n,af.getRate(s,afp,"1",hname+"_")))
	
	
	gRts=graph.graphxy(width=30,height=10,
					   x=graph.axis.lin(title=unitNames[depth],min=0,max=n),
					   y=graph.axis.lin(title="Excess Rate [Hz]"),
					   key = graph.key.key(pos="bc",columns=2))
	setTexrunner(gRts)
	
	scols = {"E":rgb.red,"W":rgb.blue}
	LF = LinearFitter(terms=[polyterm(0)])
	
	for s in scols:
		for afp in afpSymbs:
			gdat = [ [n,r.rate,r.d_rate] for (n,r) in rts[(s,afp)] ]
			LF.fit(gdat,cols=(0,1,2),errorbarWeights=True)
			err = 1.0/sqrt(LF.sumWeights())
			chi2 = LF.ssResids()
			ndf = len(gdat)-len(LF.coeffs)
			gtitle = "%s %s: $%.3f \\pm %.3f$, $\\chi^2/\\nu = %.1f/%i$"%(s,afp,LF.coeffs[0],err,chi2,ndf)
			print gtitle
			gRts.plot(graph.data.points(gdat,x=1,y=2,dy=3,title=gtitle),
					  [graph.style.symbol(afpSymbs[afp],size=0.2,symbolattrs=[scols[s],]),
					   graph.style.errorbar(errorbarattrs=[scols[s]])])
	
	gRts.writetofile(basedir+"/Rt_%s_%i.pdf"%(hname,depth))


def excess_table(dset,ltxt="",dscale=1.0):
	
	print "\n------------",dset,"--------------\n"
	
	AF = AsymmetryFile(os.environ["UCNA_ANA_PLOTS"]+"/"+dset)
	xs = [bgSubtr(b) for b in AF.dat["bg_subtr_xs"]]

	tsets = [[("Type 0","hEnergy_Type_0",1000),("Type I","hEnergy_Type_1",1000),("Type II/III","hEnergy_Type_2",1000),
			  ("$\\mu$ Backing","hBackMu",1000),("gamma","hEnergy_Type_4",1000),],
			 [("$\\beta$ 1--2.2MeV","ExcessE",1000),("$\\beta > 2.2$MeV","ExcessE",2200),
			  ("$\\gamma$ 0.2--1MeV","ExcessGamma",200),("$\\gamma$ 1--2.2MeV","ExcessGamma",1000),("$\\gamma >2.2$MeV","ExcessGamma",2200)]]
	
	for cols in tsets:
		
		tbl = "\\begin{table} \\centering \\begin{tabular}{| c ||"
		for c in cols:
			tbl += " c |"
		tbl += "}\\hline\n "
		for c in cols:
			tbl += "& %s\t"%c[0]
		tbl += "\\\\ \\hline \hline\n"
				
		for s in ["E","W"]:
			for afp in ["Off","On"]:
				if ltxt:
					if afp == "On":
						continue
					tbl += ltxt+" "
				else:
					tbl += "%s %s "%(s,afp)
				cdat = dict([(c,[x for x in xs if x.side==s and x.afp==afp and x.name==c[1] and x.eMin==c[2]]) for c in cols])
				
				for c in cols:
					bxs = cdat[c]
					assert len(bxs)==1
					bxs = bxs[0]
					tbl += "& $%i \\pm %i$\t"%(bxs.xs*dscale,bxs.d_xs*dscale)
				
				if ltxt:
					tbl += "\\\\ \\hline\n"
					continue
						
				tbl += "\\\\\nBG "
				for c in cols:
					bxs = cdat[c][0]
					tbl += "& $%i $\t"%(bxs.nBG)
				tbl += "\\\\ \\hline\n"

		tbl += "\\end{tabular}\\caption{Excess event counts after BG subtraction and size of subtracted BG}\\end{table}"

		print
		print tbl


if __name__=="__main__":
	
	excess_table("OctetAsym_Offic/OctetAsym_Offic.txt")
	
	excess_table("NGBG/ScintFace_nCaptH/ScintFace_nCaptH.txt","Scintillator",0.1)
	excess_table("NGBG/DetAl_nCaptAl/DetAl_nCaptAl.txt","Al $\\beta+\gamma$",0.1)
	excess_table("NGBG/DetAl_nCaptAlGamma/DetAl_nCaptAlGamma.txt","Al $\\gamma$",0.1)
	excess_table("NGBG/Combined/Combined.txt","Combined",1.0)
	exit(0)
	
	#plot_bgSubtrHist(os.environ["UCNA_ANA_PLOTS"]+"/OctetAsym_Annulus/")
	for i in range(3):
		plot_TypeIV_resid(os.environ["UCNA_ANA_PLOTS"]+"/OctetAsym_Offic/",i)
