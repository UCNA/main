#!/usr/bin/python

from review.Asymmetries import *
from ucnacore.PyxUtils import *

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
			chi2 = LF.chisquared()
			ndf = LF.nu()
			gtitle = "%s %s: $%.3f \\pm %.3f$, $\\chi^2/\\nu = %.1f/%i$"%(s,afp,LF.coeffs[0],err,chi2,ndf)
			print gtitle
			gRts.plot(graph.data.points(gdat,x=1,y=2,dy=3,title=gtitle),
					  [graph.style.symbol(afpSymbs[afp],size=0.2,symbolattrs=[scols[s],]),
					   graph.style.errorbar(errorbarattrs=[scols[s]])])
	
	gRts.writetofile(basedir+"/Rt_%s_%i.pdf"%(hname,depth))


class fltErr:
	def __init__(self,x,dx):
		self.x = x
		self.dx = dx
	def __repr__(self):
		return "%g~%g"%(self.x,self.dx)
	def __add__(self,other):
		if type(other)==type(self):
			return fltErr(self.x+other.x,sqrt(self.dx**2+other.dx**2))
		return fltErr(self.x+other,self.dx)

class bgSubtr(KVMap):
	def __init__(self,m):
		self.dat = m.dat
		self.loadFloats(["nBG","d_nBG","xs","d_xs","eMin","eMax"])
		self.nBG = fltErr(self.nBG,self.d_nBG)
		self.xs = fltErr(self.xs,self.d_xs)
		self.loadStrings(["side","afp","name"])

class XSFile(AsymmetryFile):
	def __init__(self,fname):
		AsymmetryFile.__init__(self,fname)
		self.bgs = dict([((b.side,b.afp,b.type),b) for b in [bgSubtr(m) for m in self.dat.get("bg_subtr_fit",[])]])
		self.XS = dict([((b.side,b.afp,b.name,b.eMin),b) for b in [bgSubtr(b) for b in self.dat["bg_subtr_xs"]]])
	
	def afpSumXS(self,side,name,emin):
		return self.XS[(side,"Off",name,emin)].xs+self.XS[(side,"On",name,emin)].xs

def excess_table(dset,ltxt="",dscale=1.0):
	
	print "\n------------",dset,"--------------\n"
	
	XS = XSFile(os.environ["UCNA_ANA_PLOTS"]+"/"+dset)
	
	cols = [("Type 0","hEnergy_Type_0",1000), ("Type I","hEnergy_Type_1",1000), ("Type II/III","hEnergy_Type_2",1000),
			("$\\beta$ 1--2.2MeV","ExcessE",1000), ("$\\beta > 2.2$MeV","ExcessE",2200),
			("$\\gamma$ 0.2--1MeV","ExcessGamma",200), ("$\\gamma$ 1--2.2MeV","ExcessGamma",1000) ,("$\\gamma >2.2$MeV","ExcessGamma",2200)]
			
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
			
			for c in cols:
				bxs = XS.XS[(s,afp,c[1],c[2])]
				tbl += "& $%i \\pm %i$\t"%(bxs.xs.x*dscale,bxs.xs.dx*dscale)
			
			if ltxt:
				tbl += "\\\\ \\hline\n"
				continue
					
			tbl += "\\\\\nBG "
			for c in cols:
				bxs = XS.XS[(s,afp,c[1],c[2])]
				tbl += "& %i\t\t"%(bxs.nBG.x)
			tbl += "\\\\ \\hline\n"

	tbl += "\\end{tabular}\\caption{Excess event counts after BG subtraction and size of subtracted BG}\\end{table}"

	print
	print tbl


def NGBG_combo_plot():

	# comparison categories
	cats = [("Type 0 $(r<50)$","hEnergy_Type_0",1000), ("Type I $(r<50)$","hEnergy_Type_1",1000),
			("$\\beta$ 1--2.2MeV","ExcessE",1000), ("$\\beta > 2.2$MeV","ExcessE",2200),
			("$\\gamma$ 1--2.2MeV","ExcessGamma",1000)]
			
	myticks = [ graph.axis.tick.tick(n,label=c[0]) for (n,c) in enumerate(cats) ]
	catAxis = graph.axis.lin(title=None,min=-0.5,max=len(cats)-0.5,parter=None,manualticks=myticks,painter=graph.axis.painter.regular(labeldist=0.1,labeldirection=graph.axis.painter.rotatetext(135)))
	
	#
	# individual MC attributes
	#
	
	XSalb =	XSFile(os.environ["UCNA_ANA_PLOTS"]+"/NGBG/DetAl_nCaptAl/DetAl_nCaptAl.txt")
	XSalg = XSFile(os.environ["UCNA_ANA_PLOTS"]+"/NGBG/DetAl_nCaptAlGamma/DetAl_nCaptAlGamma.txt")
	XSsch = XSFile(os.environ["UCNA_ANA_PLOTS"]+"/NGBG/ScintFace_nCaptH/ScintFace_nCaptH.txt")
	
	sims = [(XSalb,symbol.triangle,"Al $\\beta+\\gamma$"),(XSalg,symbol.square,"Al $\\gamma$"),(XSsch,symbol.circle,"H $\\gamma$")]

	
	gXS=graph.graphxy(width=10,height=7,
				x=catAxis,
				y=graph.axis.lin(title="events per $10^6$ captures",min=0),
				key = graph.key.key(pos="tl"))
	setTexrunner(gXS)

	for (sm,ssymb,sname) in sims:
		gdat = [ (n,sm.afpSumXS("E",c[1],c[2]).x*0.1) for (n,c) in enumerate(cats)]
		print gdat
		gXS.plot(graph.data.points(gdat,x=1,y=2,title=sname),[graph.style.symbol(ssymb),graph.style.line([style.linestyle.dotted])])

	gXS.writetofile(os.environ["UCNA_ANA_PLOTS"]+"/NGBG/NGBGSimAttrs.pdf")


	#
	# compare combined MC and data
	#
	
	# East data counts for each category
	#XSdat = XSFile(os.environ["UCNA_ANA_PLOTS"]+"/OctetAsym_Offic/OctetAsym_Offic.txt")
	Eoff = [fltErr(2157,503),fltErr(1213,454),fltErr(7435,781),fltErr(1318,478),fltErr(3522,1767)]
	Eon =  [fltErr(1237,567),fltErr(876,435),fltErr(3636,879),fltErr(456,476),fltErr(5791,1930)]
	Ecomb = [Eoff[n]+Eon[n] for n in range(len(Eoff))]
	
	# combined MC
	XSComb = XSFile(os.environ["UCNA_ANA_PLOTS"]+"/NGBG/Combined/Combined.txt")
	
	gComp=graph.graphxy(width=10,height=7,
				x=catAxis,
				y=graph.axis.lin(title="excess counts ($\\times 10^3$)",min=0),
				key = graph.key.key(pos="tl"))
	setTexrunner(gComp)

	gdat = [(n-0.1,Ecomb[n].x/1000.,Ecomb[n].dx/1000.) for (n,c) in enumerate(cats)]
	print gdat
	gComp.plot(graph.data.points(gdat,x=1,y=2,dy=3,title="East data"),[graph.style.symbol(symbol.circle),graph.style.errorbar()])
	gdat = [ (n+0.1,XSComb.afpSumXS("E",c[1],c[2]).x*0.001) for (n,c) in enumerate(cats)]
	print gdat
	gComp.plot(graph.data.points(gdat,x=1,y=2,title="Combo MC"),[graph.style.symbol(symbol.triangle,symbolattrs=[deco.filled])])
	gComp.writetofile(os.environ["UCNA_ANA_PLOTS"]+"/NGBG/NGBGDatSimCompare.pdf")

if __name__=="__main__":
	
	#NGBG_combo_plot()
	#exit(0)
	
	excess_table("OctetAsym_Offic/OctetAsym_Offic.txt")
	
	#excess_table("NGBG/ScintFace_nCaptH/ScintFace_nCaptH.txt","Scintillator",0.1)
	#excess_table("NGBG/DetAl_nCaptAl/DetAl_nCaptAl.txt","Al $\\beta+\gamma$",0.1)
	#excess_table("NGBG/DetAl_nCaptAlGamma/DetAl_nCaptAlGamma.txt","Al $\\gamma$",0.1)
	#excess_table("NGBG/Combined/Combined.txt","Combined",1.0)
	exit(0)
	
	
	#plot_bgSubtrHist(os.environ["UCNA_ANA_PLOTS"]+"/OctetAsym_Annulus/")
	for i in range(3):
		plot_TypeIV_resid(os.environ["UCNA_ANA_PLOTS"]+"/OctetAsym_Offic/",i)
