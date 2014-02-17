#!/usr/bin/python

import sys
sys.path.append("..")
from ucnacore.LinFitter import *
from ucnacore.PyxUtils import *
from ucnacore.EncalDB import *
from ucnacore.QFile import *

def alpha_eloss_plot(basePath = os.environ["UCNA_ANA_PLOTS"]+"/test/"):
	
	dat = QFile(basePath+"/AlphaEnergyLossSim.txt").dat["alpha_eloss"]
	for d in dat:
		d.loadFloats()

	e0 = 5485.56
	gdat = [(2*d.thickness,(5485.6-d.mean)/1000) for d in dat]

	g=graph.graphxy(width=15,height=15,
			   x=graph.axis.lin(title="Mylar foil thickness [$\\mu$m]",parter=graph.axis.parter.linear(tickdists=[5,1])),
			   y=graph.axis.lin(title="$^{241}$Am $\\alpha$ energy loss [MeV]"),
			   key = graph.key.key(pos="tl"))
	setTexrunner(g)

	g.plot(graph.data.points(gdat,x=1,y=2,title="Geant4 prediction"),[graph.style.symbol(symbol.circle,size=0.15)])

	LF = LinearFitter(terms=[polyterm(i) for i in range(3)])
	LF.fit(gdat)
	g.plot(graph.data.points(LF.fitcurve(0,20,200),x=1,y=2,title="$\\Delta E = %.2f + %.3f \\cdot l + %.4f \\cdot l^2$"%tuple(LF.coeffs)),[graph.style.line([style.linestyle.dashed])])

	g.plot(graph.data.points([[6,0.690]],x=1,y=2,title="Nominal 6$\\mu$m foil, 0.69 MeV measured"),[graph.style.symbol(symbol.plus,size=0.3,symbolattrs=[style.linewidth.Thick])])

	x0 = inverseFunction(LF,(1,15))(1.120)
	g.plot(graph.data.points([[x0,LF(x0)]],x=1,y=2,title="Source foil: 1.12 MeV loss $\\Rightarrow$ $l = %.1f$ $\\mu$m"%x0),[graph.style.symbol(symbol.diamond,size=0.3,symbolattrs=[style.linewidth.Thick])])


	print "6 um ->",LF(6)
	
	g.writetofile(basePath+"/AlphaEnergyLossSim.pdf")




if __name__ == "__main__":
	alpha_eloss_plot()
