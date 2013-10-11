from math import *
from bisect import bisect

# class for filling/plotting a histogram
class histogram:
	def __init__(self,nbins,xlo,xhi):
		if nbins:
			self.setBinEdges([ xlo+(xhi-xlo)*float(n)/float(nbins) for n in range(nbins+1)])
	
	def setBinEdges(self,bedge):
		self.binEdges = bedge
		self.nbins = len(bedge)-1
		if self.nbins <= 0:
			return
		self.binCenters = [ (self.binEdges[n]+self.binEdges[n+1])*0.5 for n in range(self.nbins) ]
		self.data = [ 0.0 for n in range(self.nbins) ]
		self.xlo = self.binEdges[0]
		self.xhi = self.binEdges[-1]
		self.xsum = 0
		self.xxsum = 0
		self.wsum = 0

	def fill(self,x,w=1.0,avgOutliers=True):
		b = bisect(self.binEdges,x)
		if avgOutliers or not (b<1 or b>self.nbins):
			self.xsum += w*x
			self.xxsum += w*x*x
			self.wsum += w
		if b<1 or b>self.nbins:
			return
		self.data[b-1] += w
		
	def xyData(self):
		return [ (self.binCenters[n],self.data[n]) for n in range(self.nbins) ]
		
	def lineData(self,xoff=0,yoff=0):
		l = [(self.xlo,0.0)]
		for n in range(self.nbins):
			l += [ (self.binEdges[n]+xoff,self.data[n]+yoff),(self.binEdges[n+1]+xoff,self.data[n]+yoff) ]
		l.append((self.xhi,0.0))
		return l
		
	def avg(self):
		if not self.wsum:
			return 0
		return self.xsum/self.wsum
		
	def rms(self):
		if not self.wsum:
			return 0
		d = self.xxsum/self.wsum-(self.xsum/self.wsum)**2
		if not d>0:
			return 0
		return sqrt(d)
