#!/sw/bin/python2.6

from ucnacore.PyxUtils import *
from math import *
from LinFitter import *

class fieldMap:

	def __init__(self):
		self.flatranges = []
		self.interpranges = []
		
	def addFlat(self,x0,x1,m):
		self.flatranges.append((x0,x1,m))
		if len(self.flatranges) > 1:
			self.interpranges.append((self.flatranges[-2][1],self.flatranges[-1][0],self.flatranges[-2][2],self.flatranges[-1][2]))
		
	def __call__(self,z):
		for r in self.flatranges:
			if r[0] <= z <= r[1]:
				return r[2]
		for r in self.interpranges:
			if r[0] <= z <= r[1]:
				base = 0.5*(r[2]+r[3])
				amp = 0.5*(r[2]-r[3])
				dz = (r[1]-r[0])
				l = (z-r[0])/dz
				return base + amp*cos(l*3.1415926535)
				
		return 0

	def points(self):
		p = []
		for r in self.flatranges:
			p += [(r[0],r[2]),(r[1],r[2])]
		return p
	
	def max(self):
		return max([r[2] for r in self.flatranges])
	
	def average_dips(self):
		bmax = self.max()
		n = 0.0
		xs = 0.0
		rxs = 0.0
		for x in unifrange(-1.5,1.5,200):
			n += 1
			xs += (bmax - self(x))/bmax
			rxs += sqrt((bmax - self(x))/bmax)
		print "Field max:",bmax,"Average dip:",xs/n,"Trapped:",rxs/n
		
def gen_fmap(mapdat,mapname=""):
	
	farfield = sum(mapdat)/len(mapdat)*0.6
	
	M = fieldMap()
	M.addFlat(-3.0,-2.2,farfield)
	
	fout = None
	if mapname:
		fout = open(mapname,"w")
		
		fout.write("z		B_z(tesla)\n")
		fout.write("-3.0		%f\n"%farfield)
		fout.write("-2.2		%f\n"%farfield)
	
	
	xrange = 1.5
	for (n,m) in enumerate(mapdat):
		x0 = 2.0*n*xrange/(len(mapdat)-1.0)-xrange
		M.addFlat(x0,x0,m)
		if fout:
			fout.write("%.2f		%.5f\n"%(x0,m))
	
	if fout:
		fout.write("2.2		%f\n"%farfield)
		fout.write("3.0		%f\n"%farfield)
	M.addFlat(2.2,3.0,farfield)
	
	return M
	
	
def load_fmaps_tsv(fname):
	return [ gen_fmap([float(x)/10. for x in l.split("\t")[1:16]]) for l in open(fname).readlines()[1:] ]
	
def plot_ramp(fmaps):

	print len(fmaps),"Maps..."
	
	gF=graph.graphxy(width=20,height=30,
				x=graph.axis.lin(title="Position [m]",min=-2,max=2),
				y=graph.axis.lin(title="Field [T]",min=0),
				key = graph.key.key(pos="tl"))
	gF.texrunner.set(lfs='foils17pt')

	mcols = rainbow(len(fmaps))
	
	for n,M in enumerate(fmaps):
		print n
		gF.plot(graph.data.points(M.points(),x=1,y=2,title=None), [graph.style.symbol(symbol.circle,size=0.1,symbolattrs=[mcols[n]]),])
		gF.plot(graph.data.points([(x,M(x)) for x in unifrange(-3,3,100)],x=1,y=2,title=None), [ graph.style.line([mcols[n]]),])

	gF.writetofile("MagField/Ramp.pdf")
	
def plot_fmaps(fmaps,mapnames):
	
	gF=graph.graphxy(width=20,height=10,
				x=graph.axis.lin(title="Position [m]",min=-2,max=2),
				y=graph.axis.lin(title="Field [T]",min=0.955,max=0.97),
				key = graph.key.key(pos="tl"))
	gF.texrunner.set(lfs='foils17pt')
	
	mcols = [rgb.red,rgb.blue]
	
	for (n,M) in enumerate(fmaps):
		M.average_dips()
		gF.plot(graph.data.points(M.points(),x=1,y=2,title=mapnames[n]), [graph.style.symbol(symbol.circle,symbolattrs=[mcols[n]]),])
		gF.plot(graph.data.points([(x,M(x)) for x in unifrange(-3,3,500)],x=1,y=2,title=None), [ graph.style.line([mcols[n]]),])
		gF.plot(graph.data.function("y(x)=%f"%M.max(),min=-3,max=3,title=None), [ graph.style.line([style.linestyle.dashed,mcols[n]]),])
		
		gRoot=graph.graphxy(width=20,height=4,
				x=graph.axis.lin(title="Position [m]",min=-1.5,max=1.5),
				y=graph.axis.lin(title="$\\sqrt{\\Delta B/B}$",min=0,max=0.1),
				key = graph.key.key(pos="tl"))
		gRoot.texrunner.set(lfs='foils17pt')
		bmax = M.max()
		gRoot.plot(graph.data.points([(x,sqrt(bmax-M(x))) for x in unifrange(-3,3,500)],x=1,y=2,title=None), [ graph.style.line([style.linewidth.THIck,mcols[n]]),])
		gRoot.writetofile("MagField/FRoot_%s.pdf"%mapnames[n])
		
	gF.writetofile("MagField/FMap.pdf")

fdat_20101028b = [9.633498413, 9.622409042, 9.622984363, 9.613266744, 9.61011692, 9.628310038, 9.615143825,
					9.582097653, 9.601583413, 9.624877538, 9.626377334, 9.625702384, 9.630312998, 9.626244183, 9.624789519]
fdat_20101028b = [m/10.0 for m in fdat_20101028b]
fdat_20101121 = [9.652568663, 9.641401899, 9.643352467, 9.635215422, 9.63018593, 9.634213194, 9.639542063,
					9.615458507, 9.626148018, 9.645075848, 9.648127796, 9.648129616, 9.652001974, 9.64841926, 9.646634086]
fdat_20101121 = [m/10.0 for m in fdat_20101121]

if __name__ == "__main__":
		
	m1028 = gen_fmap(fdat_20101028b[-1::-1],"../Aux/Fieldmap_20101028_b.txt")
		
	#m1121 = gen_fmap("../SummaryData/Fieldmap_20101121.txt",fdat_20101121)	
	#plot_fmaps([m1028,m1121],["October 28","November 21"])
	
	#plot_ramp(load_fmaps_tsv("MagField/10_15_2011.tsv")[:-30])