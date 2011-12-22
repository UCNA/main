#!/usr/bin/python

from LinFitter import *
from PyxUtils import *
from EncalDB import *
import os

sourceNames = { 1:"Sn113", 2:"Bi207", 3:"Sr85", 4:"Cd109", 7:"Sr90", 8:"In114", 9:"Ce139", 10:"Monosmuth", 11:"500keV" }
fancySourceNames = { 1:"$^{113}$Sn", 2:"$^{207}$Bi", 3:"$^{85}$Sr", 4:"$^{109}$Cd", 7:"$^{90}$Sr", 8:"$^{114}$In", 9:"$^{139}$Ce", 10:"Monosmuth", 11:"500keV" }		
peakNames = { 8:"$^{207}$Bi 1", 9:"$^{207}$Bi 2", 11:"$^{113}$Sn", 12:"Sr85", 13:"$^{109}$Cd", 14:"$^{114}$In", 15:"$^{139}$Ce", 16:"Monosmuth", 17:"500keV" }  
				
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
	
def get_source_lines(conn,src):
	"""Get source peaks measured for source."""
	conn.execute("SELECT side,tube,peak_num,peak_data,adc,dadc,adcwidth,erecon,derecon,ereconwidth,eta,gms,nPE FROM sourcepeaks WHERE source_id = %i"%src.sID)
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
		sline.eta = r[10]
		sline.gms = r[11]
		sline.nPE = r[12]
		sline.uid = (sline.src.sID,sline.side,sline.tube,sline.type)
		slines.append(sline)
	return slines

def gather_peakdat(conn,rlist):
	"""Collect source peak data for runs in list."""
	
	# load all lines
	srcs = []
	for rn in rlist:
		srcs += get_run_sources(conn,rn)
	print "Located",len(srcs),"sources."
	slines = []
	for src in srcs:
		slines += get_source_lines(conn,src)
		
	# connect data and simulations	
	sims = dict([(l.uid,l) for l in slines if l.simulation])
	for l in slines:
		if l.uid in sims:
			l.sim = sims[l.uid]
		else:
			print "** Missing simulation for",l.src.run,l.uid
	slines = [l for l in slines if l.uid in sims and not l.simulation]
	
	print "\twith",len(slines),"lines."
	return slines
	
def sort_peaks_by_type(slines):
	"""Sort out source peaks by peak type"""
	pks = {}
	for l in slines:
		pks.setdefault(l.type,[]).append(l)
	return pks

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
			
			
			
			
			
			
class LinearityCurve:
	"""Linearity curve produced from a set of source peaks."""
	def __init__(self,side,tube):
		self.side = side
		self.tube = tube
		self.fitter = LinearFitter(terms=[polyterm(0),polyterm(1)])
		self.axisType = graph.axis.lin
		
	def fitLinearity(self,slines):
	
		self.slines = [l for l in slines if l.side==self.side and l.tube==self.tube]
		pks = sort_peaks_by_type(self.slines)
		cP = rainbowDict(pks.keys())
		adcmax = max([l.adc for l in self.slines])

		##
		# Set up plots
		##
		xrange = [0,1500]
		yrange = [0,1500]
		if adcmax > 1500:
			xrange = (0,3000)
		if self.axisType == graph.axis.log:
			xrange = yrange[0] = [20,4000]
			
		self.gResid=graph.graphxy(width=10,height=2,
				x=self.axisType(title="PMT ADC",min=xrange[0],max=xrange[1]),
				y=graph.axis.lin(title="\\% Resid",min=-10,max=10))
		#self.gResid.texrunner.set(lfs='foils17pt')

		self.gEvis=graph.graphxy(width=10,height=10,ypos=self.gResid.height+0.5,
				x=graph.axis.linkedaxis(self.gResid.axes["x"]),
				y=self.axisType(title="Expected Light $\\eta\\cdot E_Q$",min=yrange[0],max=yrange[1]),
				key = graph.key.key(pos="tl"))
		#self.gEvis.texrunner.set(lfs='foils17pt')
		
		self.cnvs = canvas.canvas()
		self.cnvs.insert(self.gResid)
		self.cnvs.insert(self.gEvis)
		
		##
		# Plot data
		##
		combodat = []
		for k in pks:
			gdat = [ (l.adc*l.gms,l.sim.erecon*l.eta,l.dadc*l.gms,l) for l in pks[k] if l.adc > 0]
			combodat += [g for g in gdat if g[-1].src.type != "Sn113"]
			if not gdat:
				continue
			self.gEvis.plot(graph.data.points(gdat,x=1,y=2,dx=3,title=peakNames[k]),
				[graph.style.symbol(symbol.circle,size=0.1,symbolattrs=[cP[k],]),graph.style.errorbar(errorbarattrs=[cP[k],])])
		
		self.gEvis.text(7.5,3.5,"%s %i"%(self.side,self.tube+1))
		
		##
		# Fit data
		##
		if not combodat:
			print "****** No data found!"
			return
		dmin,dmax = min([p[0] for p in combodat]),max([p[0] for p in combodat])
		dmin=min(dmin,100)
		self.fitter.fit(combodat,cols=(0,1,2),errorbarWeights=True)
		for p in combodat:
			if not 1/1.2 < p[1]/self.fitter(p[0]) < 1.2:
				print "--> Check fit",p[-1].src.run,p[-1].uid
		self.fitter.fit([p for p in combodat if 1/1.2 < p[1]/self.fitter(p[0]) < 1.2])
		print "Fit",s,t,":",self.fitter.toLatex()
		self.gEvis.plot(graph.data.points([ [x,self.fitter(x)] for x in self.fitter.unifPoints(xrange[0]+0.1,xrange[1],100)],x=1,y=2,title=None),
			[graph.style.line(),])
			
		##
		# residuals plotting
		##
		for k in pks:
			gdat = [ (l.adc*l.gms,self.fitter(l.adc*l.gms),l.sim.erecon*l.eta) for l in pks[k] if l.adc > 0]
			gdat = [ (x,100.0*(y-yexp)/yexp) for (x,yexp,y) in gdat ]
			self.gResid.plot(graph.data.points(gdat,x=1,y=2,title=None),
				[graph.style.symbol(symbol.circle,size=0.1,symbolattrs=[cP[k],])])
		self.gResid.plot(graph.data.function("y(x)=0.0",title=None), [graph.style.line(lineattrs=[style.linestyle.dashed])])
		
	def plot_erecon(self,slines):
	
		self.slines = [l for l in slines if l.side==self.side and l.tube==self.tube]
		pks = sort_peaks_by_type(self.slines)
		cP = rainbowDict(pks.keys())
		
		# set up graphs
		title = "Tube %i"%(t+1)
		csize = 0.10
		if t == 4:
			csize = 0.30
			title = "Combined"
			
		rwidth = 30
		
		self.gRes=graph.graphxy(width=15,height=3,xpos=rwidth+0.5,
			x=graph.axis.lin(title="Expected Visible Energy",min=0,max=1200),
			y=graph.axis.lin(title="\\% Error",min=-10,max=10))
		#self.gRes.texrunner.set(lfs='foils17pt')
		
		rmin = min([l.src.run for l in slines])
		rmax = max([l.src.run for l in slines])
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
		#self.gRuns.texrunner.set(lfs='foils17pt')
					
		self.gEn=graph.graphxy(width=15,height=15,ypos=self.gRes.height+0.5,xpos=rwidth+0.5,
				x=graph.axis.linkedaxis(self.gRes.axes["x"]),
				y=graph.axis.linkedaxis(self.gRuns.axes["y"]),
				key = graph.key.key(pos="tl"))
		#self.gEn.texrunner.set(lfs='foils17pt')
		
		self.gWidth=graph.graphxy(width=15,height=15,
				x=graph.axis.lin(title="Expected Width [keV]",min=0),
				y=graph.axis.lin(title="Observed Width [keV]",min=0),
				key = graph.key.key(pos="tl"))
		#self.gWidth.texrunner.set(lfs='foils17pt')
		
		self.cnvs = canvas.canvas()
		self.cnvs.insert(self.gRes)
		self.cnvs.insert(self.gEn)
		self.cnvs.insert(self.gRuns)
		
		# plot
		for k in pks.keys():
			gdat = [ (q.src.run,q.sim.erecon,q.erecon,100.0*(q.erecon-q.sim.erecon)/q.sim.erecon,q.eta,q.enwidth,q.sim.enwidth) for q in pks[k]]
			self.gEn.plot(graph.data.points(gdat,x=2,y=3,title=peakNames[k]), [graph.style.symbol(symbol.circle,size=csize,symbolattrs=[cP[k]]),])
			self.gRes.plot(graph.data.points(gdat,x=2,y=4,title=None), [graph.style.symbol(symbol.circle,size=csize,symbolattrs=[cP[k]]),])
			self.gRuns.plot(graph.data.points(gdat,x=1,y=3,size=5,title=None), [ varCircle(symbolattrs=[cP[k]]),])
			self.gRuns.plot(graph.data.points(gdat,x=1,y=2,title=None), [ graph.style.line([style.linestyle.dashed,cP[k]]),])
			self.gWidth.plot(graph.data.points(gdat,x=7,y=6,title=peakNames[k]), [graph.style.symbol(symbol.circle,size=csize,symbolattrs=[cP[k]]),])
			
		self.gEn.plot(graph.data.function("y(x)=x",title=None), [graph.style.line(lineattrs=[style.linestyle.dashed,]),])
		self.gEn.text(11,1.5,"%s Reconstructed Energy"%self.side)
		self.gRes.plot(graph.data.function("y(x)=0",title=None), [graph.style.line(lineattrs=[style.linestyle.dashed,]),])
		self.gWidth.plot(graph.data.function("y(x)=x",title="$y=x$"), [graph.style.line(lineattrs=[style.linestyle.dashed,]),])
		
			
	def dbUpload(self,conn,ecid,refline_id):
		"""Upload PMT calibration curves to DB for given energy calibration ID."""
				
		# generate points with logarithmic spacing, upload as linearity curve
		lg = LogLogger(terms=[(lambda x: x)])
		lindat = [(0,0),]+[ (x,self.fitter(x)) for x in lg.unifPoints(10,4000,50) ]
		lgid = upload_graph(conn,"Tube Linearity %s %i ID=%i"%(self.side,self.tube,ecid),lindat)
		
		# reference line for anchoring energy, resolution
		refline = [l for l in self.slines if l.src.sID==refline_id][0]
		# rescaled width to true nPE width based on simulation
		simApparentnPE = (refline.sim.erecon/refline.sim.enwidth)**2
		widthscale = sqrt(simApparentnPE/refline.sim.nPE)
		print t,"width",refline.adcwidth,"Corrected width by",widthscale
		
		# upload to DB
		try:
			refEnergy = self.fitter(refline.adc*refline.gms)/refline.eta # synthesize energy for reference source based on ADC->Evis*eta curve
			print (ecid,self.side,self.tube,lgid,refline.src.x,refline.src.y,refline.adc*refline.gms,refEnergy,refline.adc,refline.adcwidth*widthscale)
			conn.execute("""INSERT INTO tube_calibration
							(ecal_id,side,quadrant,linearity_graph,encal_xpos,encal_ypos,encal_adc,encal_evis,noisecal_adc,noisecal_width)
							VALUES (%i,'%s',%i,%i,%.2f,%.2f,%.1f,%.1f,%.1f,%.1f)"""
							% (ecid,self.side,self.tube,lgid,refline.src.x,refline.src.y,refline.adc*refline.gms,refEnergy,refline.adc,refline.adcwidth*widthscale))	
		except:
			print "Tube calibration upload failed!"
			delete_graph(conn,lgid)
		
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
			
						
if __name__=="__main__":

	# set up output paths
	outpath = os.environ["UCNA_ANA_PLOTS"]+"/Sources/"
	os.system("mkdir -p %s/Linearity"%outpath)
	os.system("mkdir -p %s/Erecon"%outpath)
	os.system("mkdir -p %s/Widths"%outpath)	
	os.system("mkdir -p %s/Positions"%outpath)
	
	# calibration definitions:
	#				source runs;	gms;	calibrated range; 	E,W ref sources;	posmap
	cal_2010 = [
				(	13883,	13894,	13890,	13879,	13964,		94,		97,			2	),	# 0 first usable data + little Xe
				(	14104,	14116,	14111,	14077,	14380,		144,	147,		2	),	# 1	Columbus Day weekend + big Xe	
				(	14383,	14394,	14390,	14383,	14507,		212,	215,		2	),	# 2 Oct. 15-21 week
				(	14516,	14530,	14524,	14513,	14667,		268,	271,		2	),	# 3 Oct. 22-24 weekend
				(	14736,	14746,	14743,	14688,	14994,		330,	333,		2	),	# 4 Oct. 27-29 weekend; Nov. 12-14, including isobutane running and tilted sources
				(	15645,	15662,	15653,	15084,	15915,		437,	440,		1	),	# 5 Nov. 22-29 Thanksgiving Week
				(	15916,	15939,	15931,	15916,	100000,		553,	555,		1	)	# 6 Post-Thanksgiving
				]
				
	cal_2011 = [
				(	17233,	17249,	17238,	16983,	17279,		678,	681,		1	),	# New Sn, Ce sources; Xenon, Betas, Dead PMT W2
				#(17368,17359,17509,x,x),		# Beta Decay; PMT W0 missing pulser
				(	17517,	17527,	17522,	17517,	17734,		1125,	1128,		1	),	# Calibrations for Xe; W0 pulser still dead
				#(	17871,	17922,	17876	),	# Big Scan; W0 pulser still dead
				#(17903,17735,17956,x,x),		# Beta decay, long source runs
				(	18020,	18055,	18039,	18020,	18055,		1018,	1021,		1	),	# Old and new Cd Source; W0 pulser still dead
				(	18357,	18386,	18362,	18081,	18413,		1469,	1472,		1	),	# Beta decay, new In source, Xe; everything working now
				(	18617,	18640,	18622,	18432,	100000,		1894,	1897,		1	)	# Beta decay; PMT W4 Bi pulser very low
				]

		
	conn = open_connection() # connection to calibrations DB
	replace = True	# whether to replace previous calibration data
	
	for c in cal_2011[4:]:
	
		# make new calibrations set
		ecid = makeCalset(conn,c[3],c[4],c[2],c[7],replace)
			
		# gather source data from calibration runs
		rlist = range(c[0],c[1]+1)
		slines = gather_peakdat(conn,rlist)
		
		# fit linearity curves for each PMT
		for (sn,s) in enumerate(["East","West"]):
			for t in range(5):
				if t<4:
					LC = LinearityCurve(s,t)
					LC.fitLinearity(slines)
					LC.cnvs.writetofile(outpath+"/Linearity/ADC_v_Light_%i_%s%i.pdf"%(rlist[0],s[0],t))
					if ecid:
						LC.dbUpload(conn,ecid,c[5+sn])
				LC.plot_erecon(slines)
				LC.cnvs.writetofile(outpath+"/Erecon/Erecon_v_Etrue_%i_%s%i.pdf"%(rlist[0],s[0],t))
				LC.gWidth.writetofile(outpath+"/Widths/Widths_%i_%s%i.pdf"%(rlist[0],s[0],t))
				