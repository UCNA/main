from pyx import *
from pyx import style
from pyx import graph
from pyx.color import rgb, hsb
from pyx.graph.style import symbol
import math

def setTexrunner(g,opt='foils17pt'):
	try:
		g.texrunner.set(lfs=opt)
	except:
		pass

# error band path
def errorBand(g,gdat,nx,ny,ndy,xerr=False):
	
	g0=graph.graphxy(width=g.width,height=g.height,
				  x=graph.axis.linkedaxis(g.axes["x"]),
				  y=graph.axis.linkedaxis(g.axes["y"]))
	d1,d2 = (None,None)
	if xerr:
		d1 = g0.plot(graph.data.points([(g[nx]-g[ndy],g[ny]) for g in gdat],x=1,y=2), [graph.style.line()])
		d2 = g0.plot(graph.data.points([(g[nx]+g[ndy],g[ny]) for g in gdat],x=1,y=2), [graph.style.line()])
	else:
		d1 = g0.plot(graph.data.points([(g[nx],g[ny]-g[ndy]) for g in gdat],x=1,y=2), [graph.style.line()])
		d2 = g0.plot(graph.data.points([(g[nx],g[ny]+g[ndy]) for g in gdat],x=1,y=2), [graph.style.line()])
	g0.finish()
	p1 = d1.path
	p2 = d2.path
	
	area = (p1 << p2.reversed())
	
	return area
		
# variable size circles
class varCircle(graph.style.symbol):

	def __init__(self, sizecolumnname="size",
					   symbol=graph.style.symbol.circle,
					   symbolattrs=[deco.stroked([color.gray.black])],
					   **kwargs):
		self.sizecolumnname = sizecolumnname
		graph.style.symbol.__init__(self, symbol=symbol, symbolattrs=symbolattrs, **kwargs)

	def columnnames(self, privatedata, sharedata, agraph, columnnames):
		# register the new column names
		if self.sizecolumnname not in columnnames:
			raise ValueError("column '%s' missing" % self.sizecolumnname)
		return ([self.sizecolumnname,] +
				graph.style.symbol.columnnames(self, privatedata,
											   sharedata, agraph, columnnames))

	def drawpoint(self, privatedata, sharedata, graph, point):
		# replace the original drawpoint method by a slightly revised one
		if sharedata.vposvalid and privatedata.symbolattrs is not None:
			x_pt, y_pt = graph.vpos_pt(*sharedata.vpos)
			privatedata.symbol(privatedata.symbolcanvas, x_pt, y_pt,
							   privatedata.size_pt*point[self.sizecolumnname],
							   privatedata.symbolattrs)


# variable color and size circles
class changesymbol(graph.style.symbol):

    def __init__(self, sizecolumnname="size", colorcolumnname="color",
                       gradient=color.gradient.Rainbow,
                       symbol=graph.style.symbol.circle,
                       symbolattrs=[deco.filled, deco.stroked([color.gray.black])],
                       **kwargs):
        # add some configuration parameters and modify some other
        self.sizecolumnname = sizecolumnname
        self.colorcolumnname = colorcolumnname
        self.gradient = gradient
        graph.style.symbol.__init__(self, symbol=symbol, symbolattrs=symbolattrs, **kwargs)

    def columnnames(self, privatedata, sharedata, agraph, columnnames):
        # register the new column names
        if self.sizecolumnname not in columnnames:
            raise ValueError("column '%s' missing" % self.sizecolumnname)
        if self.colorcolumnname not in columnnames:
            raise ValueError("column '%s' missing" % self.colorcolumnname)
        return ([self.sizecolumnname, self.colorcolumnname] +
                graph.style.symbol.columnnames(self, privatedata,
                                               sharedata, agraph, columnnames))

    def drawpoint(self, privatedata, sharedata, graph, point):
        # replace the original drawpoint method by a slightly revised one
        if sharedata.vposvalid and privatedata.symbolattrs is not None:
            x_pt, y_pt = graph.vpos_pt(*sharedata.vpos)
            color = self.gradient.getcolor(point[self.colorcolumnname])
            privatedata.symbol(privatedata.symbolcanvas, x_pt, y_pt,
                               privatedata.size_pt*point[self.sizecolumnname],
                               privatedata.symbolattrs + [color])
							   

if 0:
	# variable color text
	class colortext(graph.style.text):

		def __init__(self, colorcolumnname="color",
						   gradient=color.gradient.Rainbow,
						   textattrs=[color.gray.black,],
						   **kwargs):
			# add some configuration parameters and modify some other
			self.colorcolumnname = colorcolumnname
			self.gradient = gradient
			graph.style.text.__init__(self, textattrs=textattrs, **kwargs)

		def columnnames(self, privatedata, sharedata, agraph, columnnames):
			# register the new column names
			if self.colorcolumnname not in columnnames:
				raise ValueError("column '%s' missing" % self.colorcolumnname)
			return ([self.colorcolumnname,] +
					graph.style.symbol.columnnames(self, privatedata,
												   sharedata, agraph, columnnames))

		def drawpoint(self, privatedata, sharedata, graph, point):
			# replace the original drawpoint method by a slightly revised one
			if sharedata.vposvalid and privatedata.symbolattrs is not None:
				x_pt, y_pt = graph.vpos_pt(*sharedata.vpos)
				color = self.gradient.getcolor(point[self.colorcolumnname])
				privatedata.text(privatedata.symbolcanvas, x_pt, y_pt, privatedata.size_pt, privatedata.symbolattrs + [color])
							   
							   
def rainbow(n, b=1.0):
	return [ hsb((1.0*x)/n,1,b) for x in range(n) ]
	
def rainbowDict(keys, b=1.0):
	n = len(keys)
	knew = [k for k in keys]
	knew.sort()
	return dict([ (k,hsb((1.0*x)/n,1,b)) for (x,k) in enumerate(knew) ])

# axis in fractions of pi
class piaxis(graph.axis.linear):
    def __init__(self, divisor=math.pi,
                 texter=graph.axis.texter.rational(suffix="\pi"), **kwargs):
        graph.axis.linear.__init__(self, divisor=divisor, texter=texter, **kwargs)



#-------------------------------------------------------------------------------------------------------------------

import time
import pyx

class timetexter:
	__implements__ = pyx.graph.axis.texter._Itexter

	def __init__(self):
		#label format strings, for strftime(), also easily overridden in derived classes
		#format of start-of-day string
		self.dateformat=r"$\rm{%a} \; %m%d$"
		#format of time string without seconds
		self.minuteformat=r"$%H \! \! : \! \!%M$"
		#format of time string with seconds
		self.secondsformat=r"$%H \! \! : \! \!%M \! \! : \! \!%S$"
	
	def createtext(self, tick):
		timeval=time.gmtime(float(tick.num)/tick.denom) #convert back to seconds and to time object
		seconds=int(((float(tick.num)/float(tick.denom))%86400 + 0.5))
		if seconds== 0: #tick 0 of a day gets date
			tick.label=time.strftime(self.dateformat, timeval) #at day turnover, print date
		elif seconds%60==0:
			tick.label=time.strftime(self.minuteformat, timeval)
		else:
			tick.label=time.strftime(self.secondsformat, timeval)

	def labels(self, ticks):
		labeledticks = []
		maxdecprecision = 0
		for tick in ticks:
			if tick.label is None and tick.labellevel is not None:
				labeledticks.append(tick)
				self.createtext(tick)
