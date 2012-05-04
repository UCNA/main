from pyx import *
from pyx import style
from pyx import graph
from pyx.color import rgb, hsb
from pyx.graph.style import symbol

def setTexrunner(g,opt='foils17pt'):
	try:
		g.texrunner.set(lfs=opt)
	except:
		pass
		
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
							   
							   
def rainbow(n):
	return [ hsb((1.0*x)/n,1,1) for x in range(n) ]
	
def rainbowDict(keys):
	n = len(keys)
	knew = [k for k in keys]
	knew.sort()
	return dict([ (k,hsb((1.0*x)/n,1,1)) for (x,k) in enumerate(knew) ])
