#!/usr/bin/python

import sys
sys.path.append("..")
from ucnacore.LinFitter import *
from ucnacore.PyxUtils import *
from ucnacore.EncalDB import *
from ucnacore.QFile import *
import os

peakNames = { 8:"$^{207}$Bi 1", 9:"$^{207}$Bi 2", 19:"$^{207}$Bi $1+2$", 11:"$^{113}$Sn", 12:"Sr85", 13:"$^{109}$Cd", 14:"$^{114}$In", 15:"$^{139}$Ce", 20:"$^{137}$Cs" }
peakSymbs = { 8:symbol.circle, 9:symbol.square, 11:symbol.triangle, 15:symbol.diamond }

####################
# CalDB operations
####################

def delete_calibration(conn,ecid):
	"""Delete calibration info with given calibration ID."""
	print "Deleting calibration",ecid
	conn.execute("SELECT linearity_graph FROM tube_calibration WHERE ecal_id = %i"%ecid);
	for f in conn.fetchall():
		delete_graph(conn,f[0])
	conn.execute("DELETE FROM tube_calibration WHERE ecal_id = %i"%ecid)
	conn.execute("DELETE FROM energy_calibration WHERE ecal_id = %i"%ecid)

def delete_ALL_calibrations(conn):
	"""Clear all calibration data from DB"""
	conn.execute("SELECT ecal_id FROM energy_calibration WHERE 1")
	for f in conn.fetchall():
		delete_calibration(conn,f[0])

def makeCalset(conn,r0,r1,rgms,posmap,replace=False):
	"""Set calibration definition for range of runs in calibration DB."""
	
	# check (and clear) existing calibration
	conn.execute("SELECT ecal_id FROM energy_calibration WHERE start_run = %i AND end_run = %i"%(r0,r1))
	res = conn.fetchall()
	if res and not replace:
		print "Calibration for",(r0,r1),"already loaded."
		return 0
	for r in res:
		delete_calibration(conn,r[0])
	
	# generate new ecal_id
	cmd = "INSERT INTO energy_calibration(start_run,end_run,gms_run,posmap_set_id) VALUES (%i,%i,%i,%i)"%(r0,r1,rgms,posmap)
	print cmd
	conn.execute(cmd)
	conn.execute("SELECT LAST_INSERT_ID()")
	ecid = int(conn.fetchone()[0])
	print "Added new calibration set",ecid
	return ecid

####################
# Source data collection
####################

class SourcePos:
	"""Position of a source for a run"""
	def __init__(self,sID):
		self.sID = sID
	def radius(self):
		return sqrt(self.x**2+self.y**2)

class SourceLine:
	"""Measured source spectrum line"""
	def __init__(self):
		pass

class SourceRate(KVMap):
	def __init__(self,m):
		KVMap.__init__(self)
		self.dat = m.dat
		self.loadFloats(["counts","rate","type","type0frac","sID"])
		self.loadStrings(["side","name","simulated"])

def get_run_sources(conn,rn):
		"""Get all sources present in listed run number."""
		conn.execute("SELECT source_id,side,x_pos,y_pos,x_width,y_width,counts,sourcetype FROM sources WHERE run_number = %i"%rn)
		srcs = []
		for r in conn.fetchall():
			src = SourcePos(r[0])
			src.run = rn
			src.side = r[1]
			src.x = r[2]
			src.y = r[3]
			src.wx = r[4]
			src.wy = r[5]
			src.counts = r[6]
			src.type = r[7]
			srcs.append(src)
		return srcs

def gather_sourcedat(conn,rlist):
	"""Gather sources for a list of run numbers"""
	srcs = []
	for rn in rlist:
		srcs += get_run_sources(conn,rn)
	print "Located",len(srcs),"sources."
	return srcs

def get_source_lines(conn,src,xquery=""):
	"""Get source peaks measured for source."""
	conn.execute("SELECT side,tube,peak_num,peak_data,adc,dadc,adcwidth,erecon,derecon,ereconwidth,dereconwidth,eta,gms,nPE FROM sourcepeaks WHERE source_id = %i %s"%(src.sID,xquery))
	slines = []
	for r in conn.fetchall():
		sline = SourceLine()
		sline.src = src
		sline.side = r[0]
		sline.tube = r[1]
		sline.type = r[2]
		sline.simulation = (r[3]=='simulation')
		sline.adc = r[4]
		sline.dadc = r[5]
		sline.adcwidth = r[6]
		sline.erecon = r[7]
		sline.derecon = r[8]
		sline.enwidth = r[9]
		sline.denwidth = r[10]
		sline.eta = r[11]
		sline.gms = r[12]
		sline.nPE = r[13]
		
		# unscramble occaisionally swapped points
		if sline.type == 109 and sline.erecon > 1000:
			sline.type = 119
		elif sline.type == 119 and sline.erecon < 1000:
			sline.type = 109

		sline.uid = (sline.src.sID,sline.side,sline.tube,sline.type)
		slines.append(sline)
	return slines

def sort_by_type(slines):
	"""Sort out source peaks by peak type"""
	pks = {}
	for l in slines:
		pks.setdefault(l.type,[]).append(l)
	return pks

class SourceDataCollector:

	def __init__(self,conn):
		self.conn = conn
	
	def gather_peakdat(self,rlist,xquery=""):
		"""Collect source peak data for runs in list."""
		
		# load all lines
		self.srcs = gather_sourcedat(self.conn,rlist)
		
		self.alllines = []
		for src in self.srcs:
			self.alllines += get_source_lines(self.conn, src, xquery)

		# filter out blatantly crazy fits
		self.slines = [l for l in self.alllines if 10 < l.erecon < 2000 and 5 < l.enwidth < 1000]
		
		# connect data and simulations	
		sims = dict([(l.uid,l) for l in self.slines if l.simulation])
		for l in self.slines:
			if l.uid in sims:
				l.sim = sims[l.uid]
			else:
				print "** Missing simulation for",l.src.run,l.uid
		self.slines = [l for l in self.slines if l.uid in sims and not l.simulation]
		
		print "\twith",len(self.slines),"lines."
		return self.slines

	def getTubeLines(self,side,tube):
		return [l for l in self.slines if l.side==side and l.tube==tube]

	def getRunGMS(self,rn,s,t):
		"""Infer GMS for run from any available lines"""
		gmsi = [l.gms for l in self.alllines if l.src.run == rn and l.side==s and l.tube==t]
		if not gmsi:
			print "*** GMS data not found for",rn,s,t
			return 1
		gmsi.sort()
		return gmsi[len(gmsi)/2]

	def rebase_gms(self,rn):
		"""Shift GMS relative to specified run"""
		for s in ["East","West"]:
			for t in range(4):
				gms0 = self.getRunGMS(rn,s,t)
				print "Resetting GMS base",s,t,gms0
				for l in self.slines:
					l.gms /= gms0
				for l in self.alllines:
					l.gms /= gms0


####################
# Linearity/calibration fitting
####################
		
class LinearityCurve:
	"""Linearity curve produced from a set of source peaks."""
	def __init__(self,side,tube,SDC):
		self.side = side
		self.tube = tube
		self.SDC = SDC
		
		self.prefitter = LinearFitter(terms=[polyterm(i) for i in range(2)])
		self.cnvs = None
		self.uselist = None
		
		fterms = [polyterm(i) for i in range(2)]
		
		if 0:
			if not (tube in (2,3) and side=="West"):
				fterms += [expterm(-0.03),expterm(-0.015)]
			if (tube,side) == (2,"West"):
				fterms.append(expterm(-0.03))
			if (tube,side) == (1,"East"):
				fterms = fterms[:-2] + [expterm(-0.02), expterm(-0.007)]

		self.fitter = LinearFitter(terms=fterms)
		
		#self.fitter = LinearFitter(terms=[polyterm(i) for i in range(2)])
		self.axisType = graph.axis.log
	
				
	##
	# Linearity plot & fit
	##
	def fitLinearity(self):
	
		self.slines = self.SDC.getTubeLines(self.side,self.tube)
	
		if self.uselist:
			print "\nFallback straight-line fits to specified sources",self.uselist
		
			self.slines = [l for l in self.slines if l.src.sID in self.uselist]
			self.fitter = LinearFitter(terms=[polyterm(1)])
			combodat = [ (l.adc*l.gms, l.sim.erecon*l.eta, l.enwidth, l.sim.enwidth) for l in self.slines ]
			
			print combodat
			
			self.fitter = LinearFitter(terms=[polyterm(1)])
			self.LFwid = LinearFitter(terms=[polyterm(1)])
			self.datrange = (100,2000)
			
			if len(combodat):
					self.fitter.fit(combodat, cols=(0,1))
					self.LFwid.fit(combodat,cols=(2,3))
			else:
				print "*** NO DATA FOUND *** Fallback to terrible calibration..."
				self.fitter.fit([(10,10),(100,100)])
				self.LFwid.fit([(10,10),(100,100)])
			
			return
	
	
		# estimated position map fractional error
		etaErr = 0.014
	
		pks = sort_by_type(self.slines)
		if not pks:
			print "\n\n*********",self.side,self.tube,"NO DATA FOUND!! ************\n\n"
			return
		cP = rainbowDict(pks.keys())
		if "PUBLICATION_PLOTS" in os.environ:
			for k in cP:
				cP[k] = rgb.black
		adcmax = max([l.adc for l in self.slines])

		##
		# Set up plots
		##
		xrange = [0,1500]
		yrange = [0,1500]
		if adcmax > 1500:
			xrange = (0,2000)
		if adcmax > 2000:
			xrange = (0,3000)
		if self.axisType == graph.axis.log:
			xrange = [10,4000]
			yrange = [10,2000]
			
		self.gResid=graph.graphxy(width=10,height=2,
				x=self.axisType(title="PMT ADC",min=xrange[0],max=xrange[1]),
				y=graph.axis.lin(title="\\% Resid",min=-10,max=10))
		setTexrunner(self.gResid)

		self.gEvis=graph.graphxy(width=10,height=10,ypos=self.gResid.height+0.5,
				x=graph.axis.linkedaxis(self.gResid.axes["x"]),
				y=self.axisType(title="Expected Light $\\eta\\cdot E_{\\rm vis}$",min=yrange[0],max=yrange[1]),
				key = graph.key.key(pos="tl"))
		setTexrunner(self.gEvis)
		
		self.cnvs = canvas.canvas()
		self.cnvs.insert(self.gResid)
		self.cnvs.insert(self.gEvis)
		
		##
		# Plot data
		##
		combodat = []
		for k in pks:
			gdat = [ (l.adc*l.gms, l.sim.erecon*l.eta, l.dadc*l.gms, l.sim.erecon*l.eta*etaErr, l) for l in pks[k] if 0 < l.adc < 3500 and 5 < l.sim.erecon*l.eta < 2500 ]
			gdat = [ g for g in gdat if xrange[0] < g[0] < xrange[1] and yrange[0] < g[1] < yrange[1] and 0 < g[3] < 100 ]
			combodat += [g for g in gdat if g[-1].src.radius() <= 45. and k != 11 and  g[-1].sim.erecon* g[-1].eta > 20]
			if not gdat:
				continue
			self.gEvis.plot(graph.data.points(gdat,x=1,y=2,dy=4,title=peakNames.get(k,k)),
				[graph.style.symbol(peakSymbs.get(k,symbol.circle),size=0.2,symbolattrs=[cP[k],]),graph.style.errorbar(errorbarattrs=[cP[k],])])
		
		self.gEvis.text(7.5,3.5,"%s %i"%(self.side,self.tube+1))
		
		##
		# Fit data
		##
		if not combodat:
			print "****** No data found!"
			self.cnvs=None
			return
		self.datrange = ( min([p[0] for p in combodat]), max([p[0] for p in combodat]) )
		self.prefitter.fit([p for p in combodat if p[-1].type in  [8,9,11,15]],cols=(0,1))
		trimcdat = []
		for p in combodat:
			if not (1/1.3 < p[1]/self.prefitter(p[0]) < 1.3 or abs(p[1]-self.prefitter(p[0])) < 25):
				print "--> Check fit",p[-1].src.run,p[-1].uid
			else:
				trimcdat.append(p)
		if not trimcdat:
			print "********* DATA IFFY! *********",self.side,self.tube
			trimcdat = combodat
		self.fitter.fit(trimcdat,cols=(0,1,3),errorbarWeights=True)
		print "Fit",s,t,":",self.fitter.toLatex()
		print "chi^2/nu =",self.fitter.chisquared(),"/",self.fitter.nu()
		self.fitter.displayCoeffErrCorr()
		print
		self.gEvis.plot(graph.data.points(self.make_lcurve(),x=1,y=2,title=None), [graph.style.line(),])
			
		##
		# residuals plotting
		##
		for k in pks:
			gdat = [ (l.adc*l.gms, self.fitter(l.adc*l.gms), l.sim.erecon*l.eta, l.sim.erecon*l.eta*etaErr) for l in pks[k] if l.adc > 0]
			gdat = [ (x,100.0*(y-yexp)/yexp,100*dy/yexp) for (x,yexp,y,dy) in gdat ]
			gdat = [ g for g in gdat if xrange[0] < g[0] < xrange[1] and -100 < g[1] < 100 ]
			if not gdat:
				continue
			self.gResid.plot(graph.data.points(gdat,x=1,y=2,dy=3,title=None),
				[graph.style.symbol(peakSymbs.get(k,symbol.circle),size=0.2,symbolattrs=[cP[k],]), graph.style.errorbar(errorbarattrs=[cP[k]])])
		self.gResid.plot(graph.data.function("y(x)=0.0",title=None), [graph.style.line(lineattrs=[style.linestyle.dashed])])
		
		
	##
	# Reconstructed energy plot with residuals
	##
	def plot_erecon(self):
	
		if self.tube == 4:
			self.slines = self.SDC.getTubeLines(self.side,self.tube)
		
		pks = sort_by_type(self.slines)
		if not pks:
			print "\n\n*********",self.side,self.tube,"NO DATA FOUND!! ************\n\n"
			self.cnvs=None
			return
		cP = rainbowDict(pks.keys())
		if "PUBLICATION_PLOTS" in os.environ:
			for k in cP:
				cP[k] = rgb.black
		
		# set up graphs
		title = "Tube %i"%(t+1)
		csize = 0.20
		if t == 4:
			csize = 0.30
			title = "Combined"
			
		rwidth = 30
		
		self.gRes=graph.graphxy(width=15,height=3,xpos=rwidth+0.5,
			x=graph.axis.lin(title="Expected Visible Energy",min=0,max=1200),
			y=graph.axis.lin(title="\\% Error",min=-10,max=10))
		setTexrunner(self.gRes)
		
		rmin = min([l.src.run for l in self.slines])
		rmax = max([l.src.run for l in self.slines])
		tckdist = [5,1]
		if rmax-rmin > 100:
			tckdist = [10,1]
		runaxis = graph.axis.lin(title="Run Number",min=rmin-1,max=rmax+1,
							parter=graph.axis.parter.linear(tickdists=tckdist),
							texter = graph.axis.texter.rational(),
							painter=graph.axis.painter.regular(labeldist=0.1,labeldirection=graph.axis.painter.rotatetext(135)))

		self.gRuns=graph.graphxy(width=rwidth,height=15,ypos=self.gRes.height+0.5,
				x2=runaxis,
				y=graph.axis.lin(title="Observed Visible Energy",min=0,max=1200),
				key = graph.key.key(pos="tl"))
		setTexrunner(self.gRuns)
					
		self.gEn=graph.graphxy(width=15,height=15,ypos=self.gRes.height+0.5,xpos=rwidth+0.5,
				x=graph.axis.linkedaxis(self.gRes.axes["x"]),
				y=graph.axis.linkedaxis(self.gRuns.axes["y"]),
				key = graph.key.key(pos="tl"))
		setTexrunner(self.gEn)
		
		maxWidth = 175
		if self.tube==4:
			maxWidth = 80
		self.gWidth=graph.graphxy(width=15,height=15,
				x=graph.axis.lin(title="Expected Width [keV]",min=0,max=maxWidth),
				y=graph.axis.lin(title="Observed Width [keV]",min=0,max=maxWidth),
				key = graph.key.key(pos="tl"))
		setTexrunner(self.gWidth)
		
					
		self.cnvs = canvas.canvas()
		self.cnvs.insert(self.gRes)
		self.cnvs.insert(self.gEn)
		self.cnvs.insert(self.gRuns)
		
		# plot
		combodat = []
		for k in pks:
			#		  0          1             2         3                                           4      5          6
			gdat = [ (q.src.run, q.sim.erecon, q.erecon, 100.0*(q.erecon-q.sim.erecon)/q.sim.erecon, q.eta, q.enwidth, q.denwidth,
			#           7              8               9                                     -1
						q.sim.enwidth, q.sim.denwidth, sqrt(q.sim.denwidth**2+q.denwidth**2), q) for q in pks[k]]
			gdat = [ g for g in gdat if 0 < g[5] < 1000 and 0 < g[7] < 1000]
			combodat += gdat
			if not gdat:
				continue
			self.gEn.plot(graph.data.points(gdat,x=2,y=3,title=peakNames.get(k,k)), [graph.style.symbol(peakSymbs.get(k,symbol.circle),size=csize,symbolattrs=[cP[k]]),])
			self.gRes.plot(graph.data.points(gdat,x=2,y=4,title=None), [graph.style.symbol(peakSymbs.get(k,symbol.circle),size=csize,symbolattrs=[cP[k]]),])
			self.gRuns.plot(graph.data.points(gdat,x=1,y=3,size=5,title=None), [ varCircle(symbolattrs=[cP[k]]),])
			self.gRuns.plot(graph.data.points(gdat,x=1,y=2,title=None), [ graph.style.line([style.linestyle.dashed,cP[k]]),])
			self.gWidth.plot(graph.data.points(gdat,x=8,y=6,dy=7,dx=9,title=peakNames.get(k,k)),
							 [graph.style.symbol(peakSymbs.get(k,symbol.circle),size=csize,symbolattrs=[cP[k]]),
							  graph.style.errorbar(errorbarattrs=[cP[k]])])
			
		self.gEn.plot(graph.data.function("y(x)=x",title=None), [graph.style.line(lineattrs=[]),])
		self.gEn.text(11,1.5,"%s Reconstructed Energy"%self.side)
		self.gRes.plot(graph.data.function("y(x)=0",title=None), [graph.style.line(lineattrs=[style.linestyle.dashed,]),])
		#self.gWidth.plot(graph.data.function("y(x)=x",title="$y=x$"), [graph.style.line(lineattrs=[style.linestyle.dashed,]),])
		
		######
		# Fit widths
		######
		cselect = [g for g in combodat if 1/1.25 < g[5]/g[7] < 1.25 and g[-1].src.radius() <= 45. ]
		if not cselect:
			cselect = [g for g in combodat if 1/10 < g[5]/g[7] < 10 and g[-1].src.radius() <= 45. ]
		self.LFwid = LinearFitter(terms=[polyterm(1)])
		self.LFwid.fit(cselect,cols=(7,5))
		wxmax = max([g[7] for g in cselect])
		for g in cselect:
			if not abs(g[5]-self.LFwid(g[7]))/g[9] < 4:
				print "--> Check width",g[-1].src.run,g[-1].uid

		wdat = [g for g in cselect if abs(g[5]-self.LFwid(g[7]))/g[9] < 4]
		if len(wdat)>2:
			self.LFwid.fit(wdat,cols=(7,5,9),errorbarWeights=True)
		else:
			self.LFwid.fit([(10,10),(50,50),(100,100)])
		print "Width Fit",s,t,":",self.LFwid.toLatex()
		print "chi^2/nu =",self.LFwid.chisquared(),"/",self.LFwid.nu()
		self.LFwid.displayCoeffErrCorr()
		print
		# re-fit without errorbars
		#self.LFwid.fit([g for g in cselect if 1/1.2 < g[5]/self.LFwid(g[7]) < 1.2],cols=(7,5))

		self.gWidth.plot(graph.data.points(self.LFwid.fitcurve(0,maxWidth),x=1,y=2,title="$y=%.3f \\cdot x$"%self.LFwid.coeffs[0]),
			[graph.style.line(lineattrs=[style.linestyle.dashed]),])
				
	def make_lcurve(self):
		"""Construct linearity curve points for upload or plotting"""
								
		# generate points with logarithmic spacing, upload as linearity curve
		lg = LogLogger(terms=[(lambda x: x)])
		lindat = [ (x,self.fitter(x)) for x in lg.unifPoints(10,4100,100) if self.datrange[0] < x < self.datrange[1] ]
		
		lindat = [ p for (n,p) in enumerate(lindat[:-1]) if 0 < p[1] < lindat[n+1][1] ]	# increasing only!
		if len(lindat) < 2:
			return [(0,0),(4000,4000)]
		d0 = nderiv(self.fitter, lindat[0][0])
		d1 = nderiv(self.fitter, lindat[-1][0])
		lindat = [(10, lindat[0][1] - (lindat[0][0]-10)*d0)] + lindat + [(3999, lindat[-1][1] + (3999-lindat[-1][0])*d1)]
		
		return lindat
	
	##
	# Upload linearity to calibrations DB
	##
	def dbUpload(self,conn,ecid,refline_id):
		"""Upload PMT calibration curves to DB for given energy calibration ID."""
	
		lindat = self.make_lcurve()
		lgid = upload_graph(conn,"Tube Linearity %s %i ID=%i"%(self.side,self.tube,ecid),lindat)
		
		# reference line for anchoring energy, resolution
		reflines = [l for l in self.slines if l.src.sID==refline_id]
		refline = SourceLine()
		widthscale = 1.0
		if reflines and self.uselist:
			refline = reflines[0]
		else:
			if reflines:
				refline = reflines[0]
				# rescaled width to true nPE width based on simulation
				simApparentnPE = (refline.sim.erecon/refline.sim.enwidth)**2
				widthscale = sqrt(simApparentnPE/refline.sim.nPE)
				print "MC width scale",widthscale
				widthscale *= self.LFwid.coeffs[0]
				print "Corrected for sources average",widthscale
				widthscale *= self.LFwid(refline.sim.enwidth)/refline.enwidth
				print "Corrected for reference source offset",widthscale
			else:
				print "\n\n******** Reference source",refline_id,"not found!!! Using defaults!\n"
				raw_input("Press enter to acknowledge and continue...")
				refline.adc=500
				refline.adcwidth=200
				
		print t,"width",refline.adcwidth,"Corrected width by",widthscale
		
		# upload to DB
		try:
			caldat = (ecid, self.side, self.tube, lgid, refline.adc, refline.adcwidth*widthscale)
			print caldat
			conn.execute("""INSERT INTO tube_calibration
							(ecal_id, side, quadrant, linearity_graph, noisecal_adc,noisecal_width)
							VALUES (%i,'%s',%i,%i, %.2f, %.2f)""" % caldat)
		except:
			print "Tube calibration upload failed!"
			delete_graph(conn,lgid)
			raw_input("Press enter to acknowledge and continue...")



####################
# Other plots/analysis
####################
	

def plotSourcePositions(conn,rlist):
	"""Source positions in source scan"""

	sdat = gather_sourcedat(conn,rlist)
	stypes = sort_by_type(sdat)
	
	sstyle = {	"Ce139":[style.linewidth.thick],
				"Sn113":[style.linewidth.thick,style.linestyle.dotted],
				"Bi207":[style.linewidth.thick,style.linestyle.dashed] }
				
	sname = {	"Ce139":"$^{139}$Ce",
				"Sn113":"$^{113}$Sn",
				"Bi207":"$^{207}$Bi",
				"Cd109":"$^{109}$Cd",
				"Cs137":"$^{137}$Cs",
				"In114E":"$^{114}$In (E)",
				"In114W":"$^{114}$In (W)" }
	
	for s in ['East','West']:
	
		gwid = 13
		gSourcepos=graph.graphxy(width=gwid,height=gwid,
			x=graph.axis.lin(title="x position [mm]",min=-60,max=60),
			y=graph.axis.lin(title="y position [mm]",min=-60,max=60),
			key = graph.key.key(pos="tc", columns=3))
		setTexrunner(gSourcepos)
		
		gSourcepos.dolayout()
		
		gx0,gy0 = gSourcepos.pos(0,0)
		gx1,gy1 = gSourcepos.pos(1,1)
		gscale = sqrt(0.5*( (gx0-gx1)**2 + (gy0-gy1)**2 ))
		
		# draw 50mm radius fiducial volume
		gSourcepos.stroke(path.circle(gx0, gy0, 50*gscale), [style.linestyle.dashdotted, style.linewidth.thin])
		
		for styp in stypes:
			print styp
			gSourcepos.plot(graph.data.points([(-1000,-1000)],x=1,y=2,title=sname.get(styp,styp)),[graph.style.line(lineattrs=sstyle.get(styp,[]))])
			for sc in stypes[styp]:
				if not sc.side == s:
					continue
				x,y = gSourcepos.pos(sc.x,sc.y)
				gSourcepos.stroke(path.circle(x, y, 0.5*gscale*(sc.wx+sc.wy)) ,sstyle.get(styp,[]))
				
		gSourcepos.writetofile(os.environ["UCNA_ANA_PLOTS"]+"/Sources/Positions/Pos_%i_%s.pdf"%(rlist[0],s))

def backscatterEnergy(conn, rlist):

	rlist.sort()
	slines = gather_peakdat(conn, rlist, "AND peak_num > 100 AND tube=4")
	slines0 = gather_peakdat(conn, rlist, "AND peak_num < 100 AND tube=4")
	
	# connect Type I to Type 0 lines
	tp0s = dict([(l.uid,l) for l in slines0])
	for l in slines:
		u = (l.src.sID, l.side, l.tube, l.type-100)
		l.tp0 = tp0s.get(u, None)
	
	pks = sort_by_type(slines)
	
	gB = graph.graphxy(width=15,height=15,
			x=graph.axis.lin(title="expected Type 0 - Type I difference [keV]", min=0, max=60),
			y=graph.axis.lin(title="observed Type 0 - Type I difference [keV]", min=0, max=60),
			key = graph.key.key(pos="br"))
	setTexrunner(gB)

	cP = rainbowDict(pks.keys(), 0.6)
	
	for k in pks:
		#		  0             1         2              3
		#gdat = [ (q.sim.erecon, q.erecon, q.sim.derecon, q.derecon) for q in pks[k] ]
		
		gdat = [ (q.tp0.sim.erecon - q.sim.erecon, q.tp0.erecon - q.erecon, q.sim.derecon, q.derecon) for q in pks[k] if q.tp0 ]
		
		if not gdat:
			continue
			
		gB.plot(graph.data.points(gdat,x=1,y=2,dx=3,dy=4, title=peakNames.get(k-100,k)),
			[graph.style.symbol(peakSymbs.get(k,peakSymbs.get(k-100,symbol.circle)), size=0.2, symbolattrs=[cP[k]]),
			 graph.style.errorbar(errorbarattrs=[cP[k]])])

	gB.plot(graph.data.function("y(x)=x",title=None), [graph.style.line(lineattrs=[style.linestyle.dashed,]),])
	gB.writetofile(os.environ["UCNA_ANA_PLOTS"]+"/Sources/Backscatter/Energy_%i.pdf"%(rlist[0]))




# calibration definitions:
#				source runs;	gms;	calibrated range; 	E,W ref sources;	posmap
cal_2010 = [
			(	13883,	13894,	13890,	13879,	13964,		94,		97,			161 ),	# 63 0 first usable? data + little Xe
			(	14104,	14116,	14111,	14077,	14380,		144,	147,		161 ),	# 63 1	Columbus Day weekend + big Xe
			(	14383,	14394,	14390,	14383,	14507,		212,	215,		161 ),	# 63 2 Oct. 15-21 week
			(	14516,	14530,	14524,	14513,	14667,		268,	271,		161 ),	# 63 3 Oct. 22-24 weekend
			(	14736,	14746,	14743,	14688,	14994,		330,	333,		161 ),	# 63 4 Oct. 27-29 weekend; Nov. 12-14, including isobutane running and tilted sources
			(	15645,	15662,	15653,	15084,	15915,		437,	440,		164 ),	# 65/151 5 Nov. 22-29 Thanksgiving Week
			(	15916,	15939,	15931,	15916,	100000,		553,	555,		164 )	# 65 6 Post-Thanksgiving
			]
			
cal_2011 = [
			(	17233,	17249,	17238,	16983,	17297,		678,	681,		55	),	# 0 New Sn, Ce sources; Xenon, Betas, Dead PMT W2
			(	17359,	17387,	17371,	17359,	17439,		1348,	1351,		55	),	# 1 Beta Decay; PMT W0 missing pulser
			(	17517,	17527,	17522,	17440,	17734,		1125,	1128,		55	),	# 2 Calibrations for Xe; W0 pulser still dead
			(	17871,	17922,	17892,	17735,	17955,		807,	810,		55	),	# 3 Big Scan; W0 pulser still dead
			(	18020,	18055,	18039,	18020,	18055,		1018,	1021,		55	),	# 4 Old and new Cd Source; self-calibration; W0 pulser still dead
			(	18357,	18386,	18362,	18081,	18413,		1469,	1472,		55	),	# 5 Beta decay, new In source, Xe; everything working now
			(	18617,	18640,	18622,	18432,	18683,		1894,	1897,		55	),	# 6 Beta decay; PMT W4 Bi pulser very low
			(	18745,	18768,	18750,	18712,	18994,		2113,	2116,		59	),	# 7 Start of 2012; PMT W4 pulser still low
			(	19203,	19239,	19233,	19023,	19239,		2338,	2341,		59	),	# 8 W4 Pulser now higher... drifty
			(	19347,	19377,	19359,	19347,	19544,		2387,	2390,		61	),	# 9 W4 Pulser now low...
			#(	19505,	19544	),														# Feb. 14, Cd/In only; not used for calibration
			(	19823,	19863,	19858,	19583,	20000,		2710,	2713,		61)		# 10 Feb. 16-24 Xe, Betas, long sources
			]

cal_2012 = [
			(	21087,	21098,	21094,	21087,	100000,		3146,	3149,		61),	# Bi/Ce/Sn only
			(	21299,	21328,	21314,	21274,	100000,		3302,	3305,		61),	# Thanksgiving; first Cs137
			(	21679,	21718,	21687,	21679,	100000,		4149,	4151,		61),	# Dec. 6 weekend betas
			(	21914,	21939,	21921,	21914,	100000,		4193,	4196,		61),	# Dec. 12
			(	22215,	22238,	22222,	22004,	100000,		4292,	4295,		61),
			(	22294,	22306),															# Jan. 11 Bi/Ce/Sn
			]

			
	
if __name__=="__main__":

	# set up output paths
	outpath = os.environ["UCNA_ANA_PLOTS"]+"/Sources/"
	os.system("mkdir -p %s/Linearity"%outpath)
	os.system("mkdir -p %s/Erecon"%outpath)
	os.system("mkdir -p %s/Widths"%outpath)	
	os.system("mkdir -p %s/Positions"%outpath)
	os.system("mkdir -p %s/Backscatter"%outpath)
	
	conn = open_connection() # connection to calibrations DB
	replace = False 	# whether to replace previous calibration data
	makePlots = True
	#delete_calibration(conn,8466); exit(0)


	fCalSummary = open(os.environ["UCNA_ANA_PLOTS"]+"/Sources/CalSummary.txt","w")
	
	for c in cal_2011[10:11]:
	
		#print "./ReplayManager.py -s --rmin=%i --rmax=%i < /dev/null > scriptlog.txt 2>&1 &\n"%(c[0],c[1])
		#continue
		
		rlist = range(c[0],c[1]+1)
		fCalSummary.write("\n--------- %i-%i ----------\n"%(rlist[0],rlist[-1]))
		
		# gather source data from calibration runs
		SDC= SourceDataCollector(conn)
		SDC.gather_peakdat(rlist, "AND peak_num<100")
		if len(c) >= 3:
			SDC.rebase_gms(c[2])
		
		#plotSourcePositions(conn,rlist)
		#backscatterEnergy(conn, rlist)
		#plotBackscatters(conn,rlist).writetofile(outpath+"/Backscatter/Backscatter_%i.pdf"%(rlist[0]))
		#continue
		
		# make new calibrations set
		ecid = None
		if len(c)>=8:
			ecid = makeCalset(conn,c[3],c[4],c[2],c[7],replace)
		
		# fit linearity curves for each PMT
		for (sn,s) in enumerate(["East","West"]):
			for t in range(5):
				LC = LinearityCurve(s,t,SDC)
				#LC.uselist = [1348,1351]
				if t<4:
					LC.fitLinearity()
					if LC.cnvs:
						fCalSummary.write("%s %i\t%s\n"%(s,t,LC.fitter.toLatex()))
					else:
						fCalSummary.write("%s %i\t*** CAL MISSING ***\n"%(s,t))
					if LC.cnvs and makePlots:
						LC.cnvs.writetofile(outpath+"/Linearity/ADC_v_Light_%i_%s%i.pdf"%(rlist[0],s[0],t))
				if makePlots and not LC.uselist:
					LC.plot_erecon()
			
				if ecid and t<4:
					LC.dbUpload(conn,ecid,c[5+sn])
				if not (makePlots and LC.cnvs):
						continue
				try:
					LC.cnvs.writetofile(outpath+"/Erecon/Erecon_v_Etrue_%i_%s%i.pdf"%(rlist[0],s[0],t))
					LC.gWidth.writetofile(outpath+"/Widths/Widths_%i_%s%i.pdf"%(rlist[0],s[0],t))
				except:
					print "Plotting failure for",sn,s,t,"!"

		print "source replay command:"
		print "./ReplayManager.py -s --rmin=%i --rmax=%i < /dev/null > scriptlog.txt 2>&1 &\n"%(c[0],c[1])
