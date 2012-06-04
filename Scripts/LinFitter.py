#!/usr/bin/python

from numpy import array
import numpy.linalg as linalg
from numpy import dot
from math import *
from bisect import bisect

# uniformly spaced points
def unifrange(xmin,xmax,npts):
	return [ xmin + float(i)/float(npts-1)*(xmax-xmin) for i in range(npts)]

# LaTeX formatting for polynomial term
def latexPoly(varname,n):
	if n==0:
		return "1"
	if n==1:
		return varname
	if varname[0]=='{':
		return "\\left[ %s \\right]^{%i}"%(varname,n)
	return "%s^{%i}"%(varname,n)

# zero function
def returnZero(x):
	return 0

# one function
def returnOne(x):
	return 1
	
# polynomial term for fit (use in place of lambda for better picklability)
class polyterm:
	def __init__(self,n):
		self.n = n
	def __call__(self,x):
		return x**self.n
	def toLatex(self,varname='x'):
		return latexPoly(varname,self.n)
		
# polynomial term turning on after a certain point
class halfpolyterm:
	def __init__(self,n,x0):
		self.n = n
		self.x0= x0
	def __call__(self,x):
		return (x-self.x0)**self.n * (x>self.x0)
	def toLatex(self,varname='x'):
		clhs = varname
		crhs = "%.4g"%self.x0
		vname2 ="{"+varname+"%+.4g"%(-self.x0)+"}"
		if varname[:5] == "{\\ln ":
			clhs = "{"+varname[5:]
			crhs = "%.4g"%exp(self.x0)
			vname2 = "{\\ln \\left("+varname[5:-1]+"/%.4g"%(exp(self.x0))+"\\right)}"
		return "???"
		#return "\\left( \\begin{array}{l|l} %s & %s \\geq %s \\\\ 0 & %s < %s \\end{array} \\right)"%(latexPoly(vname2,self.n),clhs,crhs,clhs,crhs)
		
class tagWrapper:
	def __init__(self,f):
		self.f = f
	def __call__(self,x):
		return (self.f(x[0]),x[1])
	
# switching wrapper for split-key fits
class idSwitcher:
	def __init__(self,f,t0=None):
		self.t0 = t0
		self.f = f
	def __call__(self,x):
		if self.t0 is None or x[1]==self.t0:
			return self.f(x[0])
		return 0
			
# easy linear fitter
class LinearFitter:

	#  constructor with terms and values for any fixed terms
	def __init__(self, terms=[polyterm(0),polyterm(1)], fixparams = {}):
		self.terms = terms
		self.fixparams = fixparams
		self.coeffs = [0.0 for t in terms]
		for n in fixparams:
			self.coeffs[n] = fixparams[n]
	
	# perform linear fit with on (x,y) pair data (using specified columns);			
	def fit(self,xydat,cols=(0,1),errorbarWeights=False):
	
		self.xdat = array([ x[cols[0]] for x in xydat])
		self.ydat = array([ x[cols[1]] for x in xydat])
		
		self.weights = []
		if len(cols) < 3:
			self.weights = array([1.0 for x in self.xdat])
		else:
			if errorbarWeights:
				self.weights = array([ 1./x[cols[2]]**2 for x in xydat])
			else:
				self.weights = array([ x[cols[2]] for x in xydat])
		
		# non-fixed terms to use in fit
		varterms = [ n for n in range(len(self.terms)) if n not in self.fixparams ]
		
		# contribution from fixed terms
		fixedshift = array([ sum([self.fixparams[n]*self.terms[n](x) for n in self.fixparams])  for x in self.xdat])
		
		# fill arrays and perform fit
		a = array([ [self.terms[n](x)*self.weights[i] for n in varterms] for (i,x) in enumerate(self.xdat) ])
		b = (array(self.ydat)-fixedshift)*self.weights
		(varcoeffs,ssr,matrixRank,singularValues) = linalg.lstsq(a,b)
		
		# merge variable and fixed terms
		self.coeffs = range(len(self.terms))
		for (n,t) in enumerate(varterms):
			self.coeffs[t] = varcoeffs[n]
		for n in self.fixparams:
			self.coeffs[n] = self.fixparams[n]

	
	# evaluate fit function at position x
	def __call__(self,x):
		return sum([f(x)*self.coeffs[n] for (n,f) in enumerate(self.terms)])
	
	# get a list of the fitted points along with the fit value
	def fittedpoints(self,dosort = True):
		fp = [ [self.xdat[i],self.ydat[i],self(self.xdat[i])] for i in range(len(self.xdat))]
		if dosort:
			fp.sort()
		return fp
		
	# sum of square of residuals (weighted = sum(chi^2))
	def ssResids(self):
		return sum([ (self.ydat[i]-self(self.xdat[i]))**2*self.weights[i] for i in range(len(self.xdat))])
	
	# sum of weights
	def sumWeights(self):
		return sum(self.weights)
	
	# rms deviation		
	def rmsDeviation(self):
		return (self.ssResids()/sum(self.weights))**0.5
	
	# get a list of npts uniformly spaced points between xmin and xmax
	def unifPoints(self,xmin,xmax,npts):
		return unifrange(xmin,xmax,npts)
	
	# plottable curve produced by fit
	def fitcurve(self,xmin,xmax,npts = 100):
		return [ [x,self(x)] for x in self.unifPoints(xmin,xmax,npts) ]

	# print coefficients
	def setCoeffs(self,c):
		self.coeffs = c
	
	# latex printable form
	def toLatex(self,varname='x',cfmt=".4g"):
		s = ""
		for (n,f) in enumerate(self.terms):
		
			coeffstr = f.toLatex(varname)
				
			if self.coeffs[n]==0:
				continue
			if self.coeffs[n]==1:
				if n==0:
					s += coeffstr
				else:
					s += "+"+coeffstr
				continue
			if self.coeffs[n]==-1:
				s += "-"+coeffstr
				continue
				
			fmt = "%+"+cfmt
			if n==0:
				fmt = "%"+cfmt
			if coeffstr == "1":
				coeffstr = ""
			else:
				coeffstr = "\\cdot "+coeffstr
			s += fmt%self.coeffs[n] + coeffstr
		return s
		
# wrap a linear fitter with axis transforms
class TransformedFitter:

	# constructor with transform and inverse-transform functions
	def __init__(self,xtrans,ytrans,xinv,yinv,**kwargs):
		self.LF = LinearFitter(**kwargs)
		self.xtrans = xtrans
		self.ytrans = ytrans
		self.xinv = xinv
		self.yinv = yinv
		self.coeffs = self.LF.coeffs
	
	# perform fit
	def fit(self,xydat,cols=(0,1)):
		if len(cols) < 3:
			self.LF.fit([(self.xtrans(x[cols[0]]),self.ytrans(x[cols[1]])) for x in xydat])
		else:
			self.LF.fit([(self.xtrans(x[cols[0]]),self.ytrans(x[cols[1]]),x[cols[2]]) for x in xydat],cols=(0,1,2))
		self.coeffs = self.LF.coeffs
	
	# evaluate at x
	def __call__(self,x):
		return self.yinv(self.LF(self.xtrans(x)))

	# return fitted points plus fit value
	def fittedpoints(self,dosort = True):
		fp = [ (self.xinv(p[0]),self.yinv(p[1]),self.yinv(p[2])) for p in self.LF.fittedpoints(False) ]
		if dosort:
			fp.sort()
		return fp
	
	# list of points uniformly spaced in transform space
	def unifPoints(self,xmin,xmax,npts):
		return [ self.xinv(x) for x in self.LF.unifPoints(self.xtrans(xmin),self.xtrans(xmax),npts) ]

	# plottable fit curve
	def fitcurve(self,xmin,xmax,npts = 100):
		return [ [x,self(x)] for x in self.unifPoints(xmin,xmax,npts) ]

	# set coefficients
	def setCoeffs(self,c):
		self.LF.setCoeffs(c)
		
# log x vs. linear y transformed fitter		
class LogXer(TransformedFitter):
	def __init__(self,**kwargs):
		TransformedFitter.__init__(self,log,polyterm(1),exp,polyterm(1),**kwargs)

# log y vs. linear x transformed fitter
class LogYer(TransformedFitter):
	def __init__(self,**kwargs):
		TransformedFitter.__init__(self,polyterm(1),log,polyterm(1),exp,**kwargs)
		
# log-log plot transformed fitter
class LogLogger(TransformedFitter):
	def __init__(self,**kwargs):
		TransformedFitter.__init__(self,log,log,exp,exp,**kwargs)
	
	# latex printable form
	def toLatex(self,varname='x'):
		return "\\exp \\left[ %s \\right]"%self.LF.toLatex("{\\ln %s}"%varname)

def nderiv(f,x,h=0.001):
	return (f(x+h)-f(x-h))/(2*h)
	
def logderiv(f,x,h=0.001):
	return nderiv(f,x,h)*x/f(x)
	
# composed chain of functions
class functionChain:

	def __init__(self):
		self.fcns = []
		self.myDomain = [0,0]
		
	def __call__(self,x):
		for f in self.fcns:
			x = f(x)
		return x

	# add new outer layer of composition
	def compose(self,f):
		if not self.fcns:
			if hasattr(f,"myDomain"):
				self.myDomain = list(f.myDomain)
		self.fcns.append(f)
		return self
	
	# add a new inner layer of composition
	def prepend(self,f):
		self.fcns.insert(0,f)
		return self		

class piecewiseExtender:
	
	def __init__(self,f0,xmin,xmax,loglog=True,loderiv=None,hideriv=None,scale=1.0):
		self.f0 = f0
		self.xmin = xmin
		self.xmax = xmax
		self.range = [0,0]
		self.scale = scale
		self.offset = 0.0
		
		if loglog:
			if xmin>0:
				self.flow = LogLogger()
				if loderiv is None:
					loderiv = logderiv(f0,xmin)
				self.flow.setCoeffs( [log(f0(xmin))-loderiv*log(xmin),loderiv] )
			else:
				self.flow = None
				
			self.fhi = LogLogger()
			if hideriv is None:
				hideriv = logderiv(f0,xmax)
			self.fhi.setCoeffs( [log(f0(xmax))-hideriv*log(xmax),hideriv] )
			
	def __call__(self,x):
		y = 0
		if x < self.xmin:
			y = self.flow(x)
		elif x <= self.xmax:
			y = self.f0(x)
		else:
			y = self.fhi(x)
		if hasattr(self,"scale"):
			y *= self.scale
		if hasattr(self,"offset"):
			y += self.offset
		return y
	
	def fitcurve(self,xmin,xmax,npts = 100):
		return [ [x,self(x)] for x in unifrange(xmin,xmax,npts) ]
	
	def unifPoints(self,xmin,xmax,npts):
		return unifrange(xmin,xmax,npts)
		
	def toLatex(self,varname='x'):
		return "\\left. %s \\right|_{%s < %i}"%(self.f0.toLatex(varname),varname,self.xmax)

class yMulShifter:
	def __init__(self,f0,x0,y0):
		self.f0 = f0
		self.mul = y0/f0(x0)
	def __call__(self,x):
		return self.f0(x)*self.mul

# function inverse class
class inverseFunction:

	def __init__(self,f,xrange,abstol = 0.001,reltol = 0.001):
		self.f = f				# function to invert
		self.abstol = abstol	# absolute error tolerance
		self.reltol = reltol	# relative error tolerance
		self.xrange = xrange	# x range for inversion
		
		self.make_invtable(10)
		
	def make_invtable(self,npts):
		invtable = [ (self.f(x),x) for x in unifrange(self.xrange[0],self.xrange[1],npts) ]
		invtable.sort()
		self.xtable = [ x[1] for x in invtable ]
		self.ytable = [ x[0] for x in invtable ]
		
	def __call__(self,y):
	
		# initial guess
		b = bisect(self.ytable,y)
		if b < 1:
			b=1
		if b >= len(self.xtable):
			b = len(self.xtable)-1
		x0 = self.xtable[b-1]
		x1 = self.xtable[b]
		y0 = self.f(x0)
		y1 = self.f(x1)
		
		# refine until tolerances met
		while abs(y1-y)>self.abstol and abs(y1-y)/y>self.reltol:
			xnew = ( (x1-x0)*y + x0*y1-x1*y0 )/(y1-y0)
			if abs(x0-xnew)<abs(x1-xnew):
				x1 = xnew
				y1 = self.f(xnew)
			else:
				x0 = xnew
				y0 = self.f(xnew)
		
		return x1

# fit several groups of data with unknown normalization	
def multiGroupFit(datagroups,fitter,normlGuess=None):
	
	# default normalization guess
	if not normlGuess:
		normlGuess = [1.0 for d in datagroups]
	assert(len(normlGuess)==len(datagroups))
	
	# fit combined, rescaled data
	combodat = []
	for (n,d) in enumerate(datagroups):
		combodat += [ (p[0]*normlGuess[n],p[1]) for p in d]
	fitter.fit(combodat)
	
	# find normalization factor for each group
	normlfit = LogLogger(terms=[inverseFunction(fitter,(10,4000))])
	for (n,d) in enumerate(datagroups):
		normlfit.fit(d)
		normlGuess[n] = exp(normlfit.LF.coeffs[0])
	
	# normalize normalizations
	navg = sum(normlGuess)/len(normlGuess)
	normlGuess = [ x/navg for x in normlGuess ]
	return normlGuess
	
					
# solve for inverse point on a function
def isolve(f,y,x0,x1,abstol = 0.001,reltol = 0.001):

	y0 = f(x0)
	y1 = f(x1)
	
	while abs(y1-y)>abstol and abs(y1-y)/y>reltol:
		xnew = ( (x1-x0)*y + x0*y1-x1*y0 )/(y1-y0)
		if abs(x0-xnew)<abs(x1-xnew):
			x1 = xnew
			y1 = f(xnew)
		else:
			x0 = xnew
			y0 = f(xnew)
	
	return x1
	
	
# bilinear term for class below
class bilinearTerm:
	def __init__(self,xki,xkj):
		self.xki=xki
		self.xkj=xkj
	def __call__(self,k):
		k = int(k)
		assert 0 <= k < len(self.xki)
		return self.xki[k]*self.xkj[k]
		
# set up bilinear matrix best fit y_k = x_ki M_ij x_kj^T for (lower triangular) M_ij
class bilinearFactory:
	def __init__(self,xki):
		self.xki = [1 for x in xki[0]] + xki
		self.n = len(self.xki)
		self.terms = []
		for i in range(self.n):
			for j in range(i+1):
				self.terms.append(bilinearTerm([x[i] for x in self.xki],[x[j] for x in self.xki]))
				
	def fit(self,yk):
		LF = LinearFitter(terms=self.terms)
		LF.fit([(k,y) for (k,y) in enumerate(yk)])
		self.coeffs = LF.coeffs
		
	def __call__(self,xi):
		s = 0
		n = 0
		row = [1,]+xi
		for i in range(self.n):
			for j in range(i+1):
				s += self.coeffs[n]*row[i]*row[j]
				n += 1
		return s
	
	