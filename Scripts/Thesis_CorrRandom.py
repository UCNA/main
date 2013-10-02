#!/sw/bin/python2.7

from math import *
from PyxUtils import *
from LinFitter import *
import os

def calc_ab(c,n):
	a = sqrt(1+(n-1)**2-(n-1)*(n-2)*c + 2*(n-1)*sqrt((1-c)*(1+(n-1)*c)))/n
	b = ( 2*c*(n-1) + (n-2)*(sqrt((1-c)*(1+(n-1)*c))-1) ) / (n**2 * a)
	#print c,n,(a,b),"a^2+(n-1)b^2-1 =",a**2+(n-1)*b**2-1,"\t2ab+(n-2)b^2-c =",2*a*b+(n-2)*b**2-c
	return a,b

def calc_de(a,b,n):
	x = (a-b)*(a+(n-1)*b)
	d = (a+(n-2)*b)/x
	e = -b/x
	#print (a,b),(d,e)
	return d,e

def calc_abde(c,n):
	a,b = calc_ab(c,n)
	d,e = calc_de(a,b,n)
	#print (c,n),(a,b),(d,e)
	return a,b,d,e

if __name__=="__main__":

	
	g=graph.graphxy(width=15,height=8,
			   x=graph.axis.lin(title="correlation $c$",min=-1./3.,max=1),
			   y=graph.axis.lin(title="$n=4$ correlating matrix terms",min=-1,max=3),
			   key = graph.key.key(pos="tc",columns=2))
	setTexrunner(g)

	gdat = [ [c,]+list(calc_abde(c,4)) for c in unifrange(-0.333,0.999,500)]
	g.plot(graph.data.points(gdat,x=1,y=2,title="a"),[graph.style.line([style.linewidth.Thick])])
	g.plot(graph.data.points(gdat,x=1,y=3,title="b"),[graph.style.line([style.linestyle.dashdotted,style.linewidth.Thick])])
	g.plot(graph.data.points(gdat,x=1,y=4,title="d"),[graph.style.line([style.linestyle.dotted])])
	g.plot(graph.data.points(gdat,x=1,y=5,title="e"),[graph.style.line([style.linestyle.dashed])])
	g.writetofile(os.environ["UCNA_ANA_PLOTS"]+"/Thesis/DecorrelationMatrix.pdf")

