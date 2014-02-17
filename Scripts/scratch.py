#!/sw/bin/python2.7

from PyxUtils import *
from LinFitter import *

def foilPlots():
	gdat = [[0, 5454.613657, 5.315878], [1, 5236.656481, 13.680301], [2, 5011.568849, 19.121518], [3, 4779.527688, 23.611639], [4, 4538.714897, 27.468423], [5, 4288.680534, 31.241567], [6, 4027.857506, 35.217324], [7, 3754.420640, 39.238778], [8, 3466.925528, 43.227529], [9, 3161.792778, 48.056306], [10, 2835.220217, 52.903874]]
	
	
	gFoo=graph.graphxy(width=15,height=15,
					   x=graph.axis.lin(title="Foil thickness [$\\mu$m]",min=0,parter=graph.axis.parter.linear(tickdists=[2,1])),
					   y=graph.axis.lin(title="Energy loss [MeV]",min=0),
					   key = graph.key.key(pos="tl"))
	setTexrunner(gFoo)
	
	e0 = 5485.56
	gdat = [(2*g[0],(e0-g[1])/1000.,g[2]/1000.) for g in gdat]
	gFoo.plot(graph.data.points(gdat,x=1,y=2,dy=3,title="MC Results"),
			  [graph.style.symbol(symbol.circle,size=0.1),
			   graph.style.errorbar()])
	
	LF = LinearFitter(terms=[polyterm(0),polyterm(1),polyterm(2)])
	LF.fit(gdat,cols=(0,1))
	gFoo.plot(graph.data.points(LF.fitcurve(0,gdat[-1][0]),x=1,y=2,title="$y=%s$"%LF.toLatex('x',".3g")),
			  [graph.style.line(lineattrs=[style.linestyle.dashed,rgb.red]),])
		
	tnom = 7.2
	gFoo.plot(graph.data.points([[tnom,LF(tnom)]],x=1,y=2,title="Nominal: %.1f$\\mu$m $\\Rightarrow$ %.2fMeV"%(tnom,LF(tnom))),
			  [graph.style.symbol(symbol.plus,size=0.2,symbolattrs=[rgb.blue])])
	
	tobs = isolve(LF,1.125,6,12)
	gFoo.plot(graph.data.points([[tobs,LF(tobs)]],x=1,y=2,title="Observed: %.3fMeV $\\Rightarrow$ %.2f$\\mu$m"%(LF(tobs),tobs)),
			  [graph.style.symbol(symbol.triangle,size=0.2,symbolattrs=[rgb.blue])])
	
	gFoo.writetofile("/Users/michael/Desktop/foils.pdf")


if __name__ == "__main__":
	
	
	foilPlots()
	exit(0)
	
	gFoo=graph.graphxy(width=10,height=10,
						  x=graph.axis.lin(title="Field Peak [Gauss]"),
						  y=graph.axis.lin(title="Instrumental Asymmetry",min=0))
	setTexrunner(gFoo)

	b0 = 310
	gdat = [(0,0),(335-b0,.49),(435-b0,1.10),(600-b0,1.71),(600-b0,1.48)]
	gFoo.plot(graph.data.points(gdat,x=1,y=2,title=None),
			  [graph.style.symbol(symbol.circle,size=0.1,symbolattrs=[rgb.red,])])
	gFoo.plot(graph.data.function("y(x)=0.0171/sqrt(%f)*100*sqrt(x*1e-5)"%((600-b0)*1e-5),title=None),
			  [graph.style.line()])
	
	gFoo.writetofile("/Users/michael/Desktop/asym.pdf")
					   
