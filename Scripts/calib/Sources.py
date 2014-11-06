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

def get_source_lines(conn, srclist, xquery=""):
	"""Get source peaks measured for sources in list."""
	
	sdict = dict([(s.sID,s) for s in srclist])
	
	#             0    1    2        3         4   5    6        7      8       9           10           11  12  13  14
	cmd = "SELECT side,tube,peak_num,peak_data,adc,dadc,adcwidth,erecon,derecon,ereconwidth,dereconwidth,eta,gms,nPE,source_id FROM sourcepeaks WHERE source_id IN ("
	cmd += join(["%i"%s for s in sdict.keys()],',') + ") %s"%xquery
	conn.execute(cmd)
	
	slines = []
	for r in conn.fetchall():
		sline = SourceLine()
		sline.src = sdict[r[14]]
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
		sline.nPE_per_keVeta = sline.nPE/(sline.erecon*sline.eta)
		
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
		self.alllines = get_source_lines(self.conn, self.srcs, xquery)

		# filter out blatantly crazy fits
		self.slines = [l for l in self.alllines if 10 < l.erecon < 2000 and 5 < l.enwidth < 1000 and 0 < l.denwidth]
		
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

	def getTubeLines(self, side, tube, linelist = None):
		if linelist is not None:
			return [l for l in self.slines if l.side==side and l.tube==tube and l.type in linelist]
		return [l for l in self.slines if l.side==side and l.tube==tube]

	def getRunGMS(self,rn,s,t):
		"""Infer GMS for run from any available lines"""
		gmsi = [l.gms for l in self.alllines if l.src.run == rn and l.side==s and l.tube==t]
		if not gmsi:
			print "*** GMS data not found for",rn,s,t
			return 1
		gmsi.sort()
		return gmsi[len(gmsi)/2]
	
	def get_median_nPE_per_keVeta(self,s,t):
		"""Get typical nPE/MeV used to simulate lines"""
		l = [l.nPE_per_keVeta for l in self.getTubeLines(s,t)]
		l.sort()
		return l[len(l)/2]

	def rebase_gms(self,rn):
		"""Shift GMS relative to specified run"""
		for s in ["East","West"]:
			for t in range(4):
				gms0 = self.getRunGMS(rn,s,t)
				print "Resetting GMS base",s,t,gms0
				if not gms0:
					continue
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
		self.slines = None
		self.intlin = None
		self.datrange = (100,2000)
		self.width_adc = 300
		self.width_dadc = 45
		
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
	
	def load_lines(self):
		if self.slines is None:
			self.slines = self.SDC.getTubeLines(self.side,self.tube)
			self.pks = sort_by_type(self.slines)
			self.cP = rainbowDict(self.pks.keys())
			if "PUBLICATION_PLOTS" in os.environ:
				for k in self.cP:
					self.cP[k] = rgb.black
				
	##
	# Linearity/Width plot & fit
	##
	def fitLinearity(self):
	
		self.load_lines()
			
		if self.uselist is not None:
			print "\nFallback straight-line fits to specified sources",self.uselist
		
			self.slines = [l for l in self.slines if l.src.sID in self.uselist]
			self.fitter = LinearFitter(terms=[polyterm(1)])
			combodat = [ (l.adc*l.gms, l.sim.erecon*l.eta, l.enwidth, l.sim.enwidth) for l in self.slines ]
			
			print combodat
			
			self.fitter = LinearFitter(terms=[polyterm(1)])
			self.LFwid = LinearFitter(terms=[polyterm(1)])
			
			if len(combodat):
					self.fitter.fit(combodat, cols=(0,1))
					self.LFwid.fit(combodat,cols=(2,3))
			else:
				print "*** NO DATA FOUND *** Fallback to terrible calibration..."
				self.fitter.fit([(10,10),(100,100)])
				self.LFwid.fit([(10,10),(100,100)])
			
			return True
	
	
		# estimated position map fractional error
		etaErr = 0.014
	
		if not self.pks:
			print "\n\n*********",self.side,self.tube,"NO DATA FOUND!! ************\n\n"
			return False
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
		#setTexrunner(self.gResid)

		self.gEvis=graph.graphxy(width=10,height=10,ypos=self.gResid.height+0.5,
				x=graph.axis.linkedaxis(self.gResid.axes["x"]),
				y=self.axisType(title="Expected Light $\\eta\\cdot E_{\\rm vis}$",min=yrange[0],max=yrange[1]),
				key = graph.key.key(pos="tl"))
		#setTexrunner(self.gEvis)
		
		self.cnvs = canvas.canvas()
		self.cnvs.insert(self.gResid)
		self.cnvs.insert(self.gEvis)
		
		##
		# Plot data
		##
		combodat = []
		for k in self.pks:
			gdat = [ (l.adc*l.gms, l.sim.erecon*l.eta, l.dadc*l.gms, l.sim.erecon*l.eta*etaErr, l) for l in self.pks[k] if 0 < l.adc < 3500 and 5 < l.sim.erecon*l.eta < 2500 ]
			gdat = [ g for g in gdat if xrange[0] < g[0] < xrange[1] and yrange[0] < g[1] < yrange[1] and 0 < g[3] < 100 ]
			combodat += [g for g in gdat if g[-1].src.radius() <= 45. and k != 11 and  g[-1].sim.erecon* g[-1].eta > 20]
			if not gdat:
				continue
			self.gEvis.plot(graph.data.points(gdat,x=1,y=2,dy=4,title=peakNames.get(k,k)),
				[graph.style.symbol(peakSymbs.get(k,symbol.circle),size=0.2,symbolattrs=[self.cP[k],]),graph.style.errorbar(errorbarattrs=[self.cP[k],])])
		
		#self.gEvis.text(7.5,3.5,"%s %i"%(self.side,self.tube+1))
		
		##
		# Fit data
		##
		if not combodat:
			print "****** No data found!"
			self.cnvs=None
			return False
		try:
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
			#print "chi^2/nu =",self.fitter.chisquared(),"/",self.fitter.nu()
			self.fitter.displayCoeffErrCorr()
			print
			self.gEvis.plot(graph.data.points(self.make_lcurve(),x=1,y=2,title=None), [graph.style.line(),])
		except:
			print "**** Fit failure!"
			self.fitter = None
			return False
		refvolt = []
                sourcedev = []
                for k in self.pks:
                        gdat = [ (l.adc*l.gms, self.fitter(l.adc*l.gms), l.sim.erecon*l.eta, l.sim.erecon*l.eta*etaErr) for l in self.pks[k] if l.adc > 0]
                        gdat = [ (x,100.0*(y-yexp)/yexp,100*dy/yexp,y-yexp) for (x,yexp,y,dy) in gdat ]
                        gdat = [ g for g in gdat if xrange[0] < g[0] < xrange[1] and -100 < g[1] < 100 ]
                        refvolt.append(max(abs(g[0]) for g in gdat))
                        sourcedev.append(max(abs(g[3]) for g in gdat))
                        if not gdat:
                                continue
                delv = max(sourcedev)
                vmax = max(refvolt)
                self.intlin = 100*delv/3000 # using estimated vmax for whole experiment
		#self.intlin = 100*delv/vmax # using vmax calculated for run group 
		self.gEvis.text(5.0,3.5,"%s %i INL = %.3f"%(self.side,self.tube+1,self.intlin))
			
		##
		# residuals plotting
		##
		for k in self.pks:
			gdat = [ (l.adc*l.gms, self.fitter(l.adc*l.gms), l.sim.erecon*l.eta, l.sim.erecon*l.eta*etaErr) for l in self.pks[k] if l.adc > 0]
			gdat = [ (x,100.0*(y-yexp)/yexp,100*dy/yexp) for (x,yexp,y,dy) in gdat ]
			gdat = [ g for g in gdat if xrange[0] < g[0] < xrange[1] and -100 < g[1] < 100 ]
			if not gdat:
				continue
			self.gResid.plot(graph.data.points(gdat,x=1,y=2,dy=3,title=None),
				[graph.style.symbol(peakSymbs.get(k,symbol.circle),size=0.2,symbolattrs=[self.cP[k],]), graph.style.errorbar(errorbarattrs=[self.cP[k]])])
		self.gResid.plot(graph.data.function("y(x)=0.0",title=None), [graph.style.line(lineattrs=[style.linestyle.dashed])])
		
		return True

	def fitWidths(self):
		
		######
		# Fit/plot widths
		######
		
		if not self.fitter:
			print "Bad linearity fitter; exiting width calculation."
			return False
			
		self.load_lines()
		
		maxWidth = 200
		if self.tube==4:
			maxWidth = 100
		self.gWidth=graph.graphxy(width=15,height=15,
				x=graph.axis.lin(title="Expected Width [keV]",min=0,max=maxWidth),
				y=graph.axis.lin(title="Observed Width [keV]",min=0,max=maxWidth),
				key = graph.key.key(pos="tl"))
		#setTexrunner(self.gWidth)
		
		csize = 0.20
		if self.tube == 4:
			csize = 0.30
			
		combodat = []
		for k in self.pks:
			#		  0          1           2              3               4                                       -1
			gdat = [ (q.enwidth, q.denwidth, q.sim.enwidth, q.sim.denwidth, sqrt(q.denwidth**2 + q.sim.denwidth**2), q) for q in self.pks[k] if 0 < q.enwidth < maxWidth and 0 < q.sim.enwidth < maxWidth ]
			combodat += gdat
			if not gdat:
				continue
			self.gWidth.plot(graph.data.points(gdat, x=3, y=1, dx=4, dy=2,title=peakNames.get(k,k)),
							 [graph.style.symbol(peakSymbs.get(k,symbol.circle),size=csize,symbolattrs=[self.cP[k]]),
							  graph.style.errorbar(errorbarattrs=[self.cP[k]])])
		
		# filter points for width fit
		cselect = [g for g in combodat if 1/1.25 < g[0]/g[2] < 1.25 and g[-1].src.radius() <= 45. ]
		if not cselect:
			cselect = [g for g in combodat if 1/10 < g[0]/g[2] < 10 and g[-1].src.radius() <= 45. ]
		
		# fit with 'a*x'; remove outliers
		self.LFwid = LinearFitter(terms=[polyterm(1)])
		self.LFwid.fit(cselect,cols=(2,0))
		wxmax = max([g[2] for g in cselect])
		for g in cselect:
			if not abs(g[0]-self.LFwid(g[2]))/g[4] < 4:
				print "--> Check width",g[-1].src.run,g[-1].uid

		wdat = [g for g in cselect if abs(g[0]-self.LFwid(g[2]))/g[4] < 4]
		if len(wdat):
			self.LFwid.fit(wdat,cols=(2,0,4),errorbarWeights=True)
		else:
			print "*** Insufficient width fit data! Keeping defaults."
			return False
		print "Width Fit",s,t,":",self.LFwid.toLatex()
		#print "chi^2/nu =",self.LFwid.chisquared(),"/",self.LFwid.nu()
		self.LFwid.displayCoeffErrCorr()
		print

		self.gWidth.plot(graph.data.points(self.LFwid.fitcurve(0,maxWidth),x=1,y=2,title="$y=%.3f \\cdot x$"%self.LFwid.coeffs[0]),
			[graph.style.line(lineattrs=[style.linestyle.dashed]),])
		
		# calculate updated energy resolution at width_adc channels
		if self.tube < 4:
			old_pE_per_keV = self.SDC.get_median_nPE_per_keVeta(self.side,self.tube)
			pE_per_keV = old_pE_per_keV/self.LFwid.coeffs[0]**2
			print "Updating from %.2f to %.2f photoelectrons per MeV"%(1000*old_pE_per_keV,1000*pE_per_keV)
			width_light = self.fitter(self.width_adc)
			try:
				width_dlight = sqrt(width_light/pE_per_keV)
				self.width_dadc = width_dlight/nderiv(self.fitter, self.width_adc)
			except:
				self.width_dadc = None
			if not 30 < self.width_dadc < 90:
				print "Crazy width result",self.width_dadc,"--- returning to default"
				self.width_dadc = 45
			print "Energy resolution +/-%.2f channels at %.2f channels."%(self.width_dadc,self.width_adc)
		
		return True


		
	##
	# Reconstructed energy plot with residuals
	##
	def plot_erecon(self):
	
		self.load_lines()
		
		if not self.pks:
			print "\n\n*********",self.side,self.tube,"NO DATA FOUND!! ************\n\n"
			self.cnvs=None
			return False
		
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
		#setTexrunner(self.gRes)
		
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
		#setTexrunner(self.gRuns)
					
		self.gEn=graph.graphxy(width=15,height=15,ypos=self.gRes.height+0.5,xpos=rwidth+0.5,
				x=graph.axis.linkedaxis(self.gRes.axes["x"]),
				y=graph.axis.linkedaxis(self.gRuns.axes["y"]),
				key = graph.key.key(pos="tl"))
		#setTexrunner(self.gEn)
		
		self.cnvs = canvas.canvas()
		self.cnvs.insert(self.gRes)
		self.cnvs.insert(self.gEn)
		self.cnvs.insert(self.gRuns)
		
		# plot
		combodat = []
		for k in self.pks:
			#		  0          1             2         3                                           4     	-1
			gdat = [ (q.src.run, q.sim.erecon, q.erecon, 100.0*(q.erecon-q.sim.erecon)/q.sim.erecon, q.eta,	q) for q in self.pks[k]]
			combodat += gdat
			if not gdat:
				continue
			self.gEn.plot(graph.data.points(gdat,x=2,y=3,title=peakNames.get(k,k)), [graph.style.symbol(peakSymbs.get(k,symbol.circle),size=csize,symbolattrs=[self.cP[k]]),])
			self.gRes.plot(graph.data.points(gdat,x=2,y=4,title=None), [graph.style.symbol(peakSymbs.get(k,symbol.circle),size=csize,symbolattrs=[self.cP[k]]),])
			#self.gRuns.plot(graph.data.points(gdat,x=1,y=3,size=5,title=None), [ varCircle(symbolattrs=[self.cP[k]]),])
			self.gRuns.plot(graph.data.points(gdat,x=1,y=3,title=None), [graph.style.symbol(peakSymbs.get(k,symbol.circle),size=csize,symbolattrs=[self.cP[k]]),])
			self.gRuns.plot(graph.data.points(gdat,x=1,y=2,title=None), [ graph.style.line([style.linestyle.dashed,self.cP[k]]),])

		self.gEn.plot(graph.data.function("y(x)=x",title=None), [graph.style.line(lineattrs=[]),])
		self.gEn.text(11,1.5,"%s Reconstructed Energy"%self.side)
		self.gRes.plot(graph.data.function("y(x)=0",title=None), [graph.style.line(lineattrs=[style.linestyle.dashed,]),])
		
								
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
	def dbUpload(self,conn,ecid):
		"""Upload PMT calibration curves to DB for given energy calibration ID."""
	
		if not ecid:
			return
			
		lindat = self.make_lcurve()
		lgid = upload_graph(conn,"Tube Linearity %s %i ID=%i"%(self.side,self.tube,ecid),lindat)
		caldat = (ecid, self.side, self.tube, lgid, self.width_adc, self.width_dadc)
		cmd = """INSERT INTO tube_calibration (ecal_id, side, quadrant, linearity_graph, noisecal_adc, noisecal_width) VALUES (%i,'%s',%i, %i, %.2f, %.2f)""" % caldat
		print cmd
		conn.execute(cmd)



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
		#setTexrunner(gSourcepos)
		
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
				if not (0.6 < sc.wx < 6 and 0.6 < sc.wy < 6):
					print "Strange size",sc.wx,sc.wy,sc.sID
				gSourcepos.stroke(path.circle(x, y, 0.5*gscale*(sc.wx+sc.wy)) ,sstyle.get(styp,[]))
				
		gSourcepos.writetofile(os.environ["UCNA_ANA_PLOTS"]+"/Sources/Positions/Pos_%i_%s.pdf"%(rlist[0],s))

def backscatterEnergy(conn, rlist):

	rlist.sort()
	SDC = SourceDataCollector(conn)
	SDC.gather_peakdat(rlist, "AND peak_num > 100 AND tube=4")
	slines = SDC.slines
	SDC.gather_peakdat(rlist, "AND peak_num < 100 AND tube=4")
	slines0 = SDC.slines
	
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
	#setTexrunner(gB)

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
#				source runs;	gms;	calibrated range;	posmap
cal_2010 = [
			(	13883,	13894,	13890,	13879,	13964,		161 ),	# 63 0 first usable? data + little Xe
			(	14104,	14116,	14111,	14077,	14380,		161 ),	# 63 1	Columbus Day weekend + big Xe
			(	14383,	14394,	14390,	14383,	14507,		161 ),	# 63 2 Oct. 15-21 week
			(	14516,	14530,	14524,	14513,	14667,		161 ),	# 63 3 Oct. 22-24 weekend
			(	14736,	14746,	14743,	14688,	14994,		161 ),	# 63 4 Oct. 27-29 weekend; Nov. 12-14, including isobutane running and tilted sources
			(	15645,	15662,	15653,	15084,	15915,		164 ),	# 65/151 5 Nov. 22-29 Thanksgiving Week
			(	15916,	15939,	15931,	15916,	100000,		164 )	# 65 6 Post-Thanksgiving
			]
			
cal_2011 = [
#			(	17233,	17249,	17238,	16983,	17297,		233	),	# 0 New Sn, Ce sources; Xenon, Betas, Dead PMT W2
			(	17359,	17387,	17371,	17359,	17439,		233	),	# 1 Beta Decay; PMT W0 missing pulser
			(	17517,	17527,	17522,	17440,	17734,		233	),	# 2 Calibrations for Xe; W0 pulser still dead
			(	17871,	17941,	17892,	17735,	17955,		233	),	# 3 Big Scan; W0 pulser still dead; originally to 17922
			(	18020,	18055,	18039,	18020,	18055,		233	),	# 4 Old and new Cd Source; self-calibration; W0 pulser still dead
			(	18357,	18386,	18362,	18081,	18386,		223	),	# 5 Beta decay, new In source, Xe; PMT W4 pulser low & drifty
			(	18617,	18640,	18622,	18390,	18683,		225	),	# 6 Beta decay; PMT W4 Bi pulser very low
			(	18745,	18768,	18750,	18712,	18994,		227	),	# 7 Start of 2012; PMT W4 pulser still low
			(	19203,	19239,	19233,	19023,	19239,		227	),	# 8 W4 Pulser now higher... drifty
			(	19347,	19377,	19359,	19347,	19544,		229	),	# 9 W4 Pulser now low...
			#(	19505,	19544, 19539	),					# Feb. 14, Cd2In only; not used for calibration
			(	19823,	19863,	19858,	19583,	20000,		229)	# 10 Feb. 16-24 Xe, Betas, long sources
			]

cal_2012 = [
			(	20515,	20527,	20522,	20121,	20741,		61),	# 0 Oct. 12, Bi/Ce/Sn
			(	20818,	20830,	20823,	20782,	20837,		61),	# 1 Oct. 20, Bi/Ce/Sn, calibrates 1 octet
			(	21087,	21099,	21092,	21087,	21237,		61),	# 2 Nov. 16; Bi/Ce/Sn
			(	21299,	21328,	21314,	21274,	21623,		61),	# 3 Nov. 20, Thanksgiving; Cd, In, +first Cs137
			(	21679,	21718,	21687,	21625,	21863,		61),	# 4 Dec. 6 weekend betas
			(	21914,	21939,	21921,	21890,	22118,		61),	# 5 Dec. 12
			(	22215,	22238,	22222,	22124,	22238,		61),	# 6 Dec. 18
			#(	22294,	22306),										#   Jan. 11; Bi/Ce/Sn only
			(	22437,	22462,	22442,	22294,	22630,		61),	# 7 Jan. 14
			(	22767,	22793,	22772,	22631,	100000,		61)		# 8 Jan. 25
			]


# some useful DB commands:
# UPDATE energy_calibration SET posmap_set_id=213 WHERE posmap_set_id = 211;

if __name__=="__main__":

	# set up output paths
	outpath = os.environ["UCNA_ANA_PLOTS"]+"/Sources/"
	os.system("mkdir -p %s/Linearity"%outpath)
	os.system("mkdir -p %s/Erecon"%outpath)
	os.system("mkdir -p %s/Widths"%outpath)	
	os.system("mkdir -p %s/Positions"%outpath)
	os.system("mkdir -p %s/Backscatter"%outpath)
	
	conn = open_connection() # connection to calibrations DB
	replace = False		# whether to replace previous calibration data
	makePlots = True	# whether to output linearity/width/Erecon plots
	#delete_calibration(conn,9752); exit(0)


	fCalSummary = open(os.environ["UCNA_ANA_PLOTS"]+"/Sources/CalSummary.txt","w")
	allintl = [[[] for i in range(4)]for j in range(2)]
	allruns = []
	for c in cal_2011:
		
		#print "./ReplayManager.py -s --rmin=%i --rmax=%i"%(c[0],c[1]); continue
		
		rlist = range(c[0],c[1]+1)
		calname = "\n--------- %i-%i ----------\n"%(rlist[0],rlist[-1])
		fCalSummary.write(calname)
		print calname
		
		#plotSourcePositions(conn,rlist); continue
		#backscatterEnergy(conn, rlist); continue
		#plotBackscatters(conn,rlist).writetofile(outpath+"/Backscatter/Backscatter_%i.pdf"%(rlist[0])); continue
		
		# gather source data from calibration runs
		SDC= SourceDataCollector(conn)
		SDC.gather_peakdat(rlist, "AND peak_num<100")
		if len(c) >= 3:
			SDC.rebase_gms(c[2])
		
				
		# fit linearity curves for each PMT
		calib_OK = True
		LC = {}
		k = 0
		for s in ["East","West"]:
			for t in range(5):
				print "\n-----",s,t,rlist[0],"-----"
				LC[(s,t)] = LinearityCurve(s,t,SDC)
				# LC[(s,t)].uselist = []	# set this for "emergency" defaults calibration, hopefully good enough to bootstrap to real calibration
				if t<4:
					if LC[(s,t)].fitLinearity():
						allintl[k][t].append(LC[(s,t)].intlin)
						fCalSummary.write("%s %i\t%s\n"%(s,t,LC[(s,t)].fitter.toLatex()))
						if makePlots and LC[(s,t)].uselist is None:
							LC[(s,t)].cnvs.writetofile(outpath+"/Linearity/ADC_v_Light_%i_%s%i.pdf"%(rlist[0],s[0],t))
					else:
						fCalSummary.write("%s %i\t*** CAL MISSING ***\n"%(s,t))
						calib_OK = False
				
				if LC[(s,t)].fitWidths():
					if makePlots and LC[(s,t)].uselist is None:
						LC[(s,t)].gWidth.writetofile(outpath+"/Widths/Widths_%i_%s%i.pdf"%(rlist[0],s[0],t))
				else:
					pass
					#calib_OK = False
					
				# reconstructed energy and widths plots
				if makePlots and LC[(s,t)].uselist is None:
					LC[(s,t)].plot_erecon()
					try:
						LC[(s,t)].cnvs.writetofile(outpath+"/Erecon/Erecon_v_Etrue_%i_%s%i.pdf"%(rlist[0],s[0],t))
					except:
						print "Plotting failure for",s,t,"!"
			k += 1
		allruns.append(rlist[0])
		# make new calibrations set; upload to DB
		if not calib_OK:
			print "\n*** Calibration curve generation failure; DB not updated!"
		if calib_OK and len(c) >= 6:
			print "\n --- Uploading calibrations... ---"
			ecid = makeCalset(conn,c[3],c[4],c[2],c[5],replace)
			klist = LC.keys()
			klist.sort()
			for k in klist:
				if k[1] < 4:
					LC[k].dbUpload(conn,ecid)

		print "\nsource replay command:"
		print "nohup ./ReplayManager.py -s --rmin=%i --rmax=%i < /dev/null > scriptlog.txt 2>&1 &\n"%(c[0],c[1])
	gheight = 3
	goff = 0
	cIntLin = canvas.canvas()	
	
	for i in range(2): 
		if (i==0): compass = 'East'
		else: compass = 'West'
		for t in range(4):
			gIntLin=graph.graphxy(width=15,height=gheight,ypos=goff,
                                x=graph.axis.lin(title="Run Number",min=17000,max=20000),
                                y=graph.axis.lin(title="INL(\\%)",min=0,max=8))
			cIntLin.insert(gIntLin)
			gIntLin.plot(graph.data.values(x=allruns,y=allintl[i][t]),
                                [graph.style.symbol(symbol.circle,size=0.15,symbolattrs=[rgb.red,deco.filled])])
			gIntLin.plot(graph.data.function("y(x)=5.0",title=None), [graph.style.line(lineattrs=[style.linestyle.dashed])])
			gIntLin.text(0.7,goff+2.2,"%s %i"%(compass,t))
			goff += gheight+1.4
		goff += 0.5
	cIntLin.writetofile(outpath+"/Linearity/integrallin.pdf")
