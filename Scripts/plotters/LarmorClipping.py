#!/sw/bin/python2.7

import sys
sys.path.append("..")

from ucnacore.PyxUtils import *
from math import *
from ucnacore.LinFitter import *
#from UCNAUtils import *
from bisect import bisect
from calib.FieldMapGen import *

def clip_function(y,rho,h,R):
	sqd = sqrt(rho**2-y**2)
	if sqd==0:
		sqd = 1e-10
	return h*rho**2/R*atan(y/sqd)+2*sqd/(3*R)*(3*h*y/2+rho**2-y**2)

def survival_fraction(h,rho,R):
	d = R-h
	if d < -rho:
		return 1
	if h <= -rho:
		return 0
	c1 = 0
	if  d < rho:
		sqd = sqrt(rho**2-d**2)
		c1 = pi/2*rho**2-d*sqd-rho**2*atan(d/sqd)
	return ( c1 + clip_function(min(h,rho),rho,h,R)
				- clip_function(max(h-R,-rho),rho,h,R))/(pi*rho**2)

def radial_clip_function(r,rho,h,R):
	return r**2*(3*h-2*r)/(6*R**2)
	
def radial_survival_fraction(h,rho,R):
	d = h-R
	if d > rho:
		return 1
	if h <= 0:
		return 0
	c1 = 0
	if d > 0:
		c1 = (h-R)**2
		
	return ( c1 + radial_clip_function(min(h,rho),rho,h,R) - radial_clip_function(max(d,0),rho,h,R) )/(rho**2)


class rot3:
	def __init__(self,t1,t2,t3,s=1.0):
		self.c1,self.s1 = cos(t1),sin(t1)
		self.c2,self.s2 = cos(t2),sin(t2)
		self.c3,self.s3 = cos(t3),sin(t3)
		self.s = s
	def __call__(self,(x,y,z)):
		x,y = self.c1*x+self.s1*y,self.c1*y-self.s1*x
		y,z = self.c2*y+self.s2*z,self.c2*z-self.s2*y
		z,x = self.c3*z+self.s3*x,self.c3*x-self.s3*z
		return self.s*x,self.s*y,self.s*z


class path3d:

	def __init__(self):
		self.pts = []
		self.sty = []
		self.endsty = []
		self.breakunder = False
		self.nopatch = False
		
	def addpt(self,(x,y,z),s=1):
		self.pts.append((x*s,y*s,z*s))
		
	def apply(self,transf):
		self.pts = [transf(xyz) for xyz in self.pts]
	
	def finish(self):
		self.p = path.path()
		self.p.append(path.moveto(self.pts[0][0],self.pts[0][1]))
		for g in self.pts[1:]:	
			self.p.append(path.lineto(g[0],g[1]))
		self.patchpts = []
		self.underpts = []

	def nearestpt(self,(x,y)):
		d0 = 1e20
		n = None
		for i in range(len(self.pts)):
			d1 = (self.pts[i][0]-x)**2+(self.pts[i][1]-y)**2
			if d1 < d0:
				d0 = d1
				n = i
		return n
		
	def znear(self,(x,y)):
		return self.pts[self.nearestpt((x,y))][2]
	
	def znearc(self,c):
		x,y = self.p.at(c)
		x,y = 100*x.t,100*y.t
		return self.znear((x,y))
	
	def addPatch(self,c,z):
		self.patchpts.append((c,z))
	
	def drawto(self,cnvs):
		cnvs.stroke(self.p,self.sty)
		
		
		
	
def interleave(p3d1,p3d2):
	print "Finding intersection points..."
	is1,is2 = p3d1.p.intersect(p3d2.p)
	print "determining patch z..."
	assert len(is1)==len(is2)
	for i in range(len(is1)):
		z1 = p3d1.znearc(is1[i])
		z2 = p3d2.znearc(is2[i])
		if z1>z2:
			p3d1.addPatch(is1[i],z1)
			p3d2.underpts.append(is2[i])
		else:
			p3d2.addPatch(is2[i],z2)
			p3d1.underpts.append(is1[i])
	print "done."
			
def drawInterleaved(c,ps):
	print "Drawing base curves..."
	for p in ps:
		p.p = p.p.normpath()
		if p.breakunder:
			splits = []
			for s in p.underpts:
				splits += [s-p.breakunder*0.5,s+p.breakunder*0.5]
			psplit = p.p.split(splits)
			for seg in psplit[0::2]:
				c.stroke(seg,p.sty)
		else:
			c.stroke(p.p,p.sty+p.endsty)
	print "Preparing patches..."
	patches = []
	for (pn,p) in enumerate(ps):
		if p.nopatch:
			continue
		p.patchpts.sort()
		splits = []
		for s in p.patchpts:
			splits += [s[0]-0.05,s[0]+0.05]
		psplit = p.p.split(splits)
		patches += [ (patch[1],pn,psplit[2*n+1]) for n,patch in enumerate(p.patchpts) ]
	patches.sort()
	print "Patching intersections..."
	for p in patches:
		c.stroke(p[2],ps[p[1]].sty)
	print "Done."

def fieldPath(fmap,z0,z1,c,cmax,npts=50):
	pfield = path3d()
	for z in unifrange(z0,z1,npts):
		Bdens = c/sqrt(fmap(z)+0.0001)
		if abs(Bdens) < cmax:
			pfield.addpt((0,Bdens,z))
	return pfield

def larmor_unif(fT,theta,KE,t):
	b = electron_beta(KE)
	z = t*b*cos(theta)*3e8				# m
	r = 3.3e-6*b*(KE+511)*sin(theta)/fT	# m
	f = 2.8e10*fT						# Hz
	return r*cos(2*pi*f*t),r*sin(2*pi*f*t),z
	
def larmor_step(p,pt2_per_B,fT):
	nu = 2.8e10*fT*2*pi				# angular frequency, Hz
	pt = sqrt(fT*pt2_per_B)			# transverse momentum component, keV
	if p<=pt:
		return 0,nu
	pl = sqrt(p**2-pt**2)			# longitudinal momentum, keV
	vz = pl/sqrt(p*p+511*511)*3e8;	# z velocity, m/s
	
	return vz,nu
	
def larmorPath(fmap,p,pt2_per_B,z0,z1,dt,theta=0):
	lpath = path3d()
	z = z0
	vz = 1
	while z0 <= z <= z1 and vz>0:
		fT = fmap(z)						# magnetic field, T
		r = 3.3e-6*sqrt(pt2_per_B/fT)		# larmor radius, m
		lpath.addpt((r*cos(theta),r*sin(theta),z))
		# step to next point
		vz,nu = larmor_step(p,pt2_per_B,fmap(z))
		theta += nu*dt
		z += vz*dt
	return lpath
		
					
def plot_larmor_trajectory():
	
	fmap = fieldMap()
	fmap.addFlat(-1.0,0.01,1.0)
	fmap.addFlat(0.015,1.0,0.6)
	#fmap.addFlat(-1.0,0.01,0.6)
	#fmap.addFlat(0.08,1.0,1.0)
	
	fT = fmap(0)
	theta = 1.4
	KE = 511.
	#rot = rot3(0,0.0,-pi/2-0.2,500)
	rot = rot3(0,0.0,-pi/2+0.2,500)
	tm = 1e-9
	doFinal = True
	
	plarmor = larmorPath(fmap,500,495**2/fmap(0),0,0.02,5e-13,3*pi/4)
	plarmor.apply(rot)
	#plarmor.sty = [style.linewidth.thick,rgb.red]
	plarmor.sty = [style.linewidth.thick]
	plarmor.endsty = [deco.earrow()]
	plarmor.finish()
	x0,y0 = plarmor.p.at(plarmor.p.begin())
	
	fieldlines = []
	w = 0.0025
	cmagf = canvas.canvas()
	for o in unifrange(-w,w,20):
		pf = fieldPath(fmap,-0.002,0.022,o,1.02*w)
		if len(pf.pts) < 10:
			continue
		pf.apply(rot)
		pf.finish()
		pf.breakunder = 0.07
		pf.nopatch = True
		#pf.sty=[style.linewidth.thin,rgb.blue]
		pf.sty=[style.linewidth.thin]	# field line color/style
		fieldlines.append(pf)
		pf.drawto(cmagf)
		if doFinal:
			interleave(plarmor,pf)
	#cmagf.stroke(path.circle(x0,y0,0.07),[deco.filled([rgb.green])])
	cmagf.stroke(path.circle(x0,y0,0.07),[deco.filled([rgb.white]),style.linewidth.Thick])
	cmagf.writetofile("/Users/michael/Desktop/Bfield.pdf")	
	
			
	c = canvas.canvas()
	if doFinal:
		drawInterleaved(c,[plarmor,]+fieldlines)
	else:
		plarmor.drawto(c)
		for pf in fieldlines:
			pf.drawto(c)
			
	
	#c.stroke(path.circle(x0,y0,0.07),[deco.filled([rgb.green])])
	c.stroke(path.circle(x0,y0,0.07),[deco.filled([rgb.white]),style.linewidth.Thick])
	
	c.writetofile("/Users/michael/Desktop/larmor_spiral.pdf")
	

def plot_spectrometer_field():

	fmap = fieldMap()
	fmap.addFlat(-3,-2.8,0.01)
	fmap.addFlat(-2.3,-2.1,0.6)
	fmap.addFlat(-1.6,1.6,1.0)
	fmap.addFlat(2.1,2.3,0.6)
	fmap.addFlat(2.8,3,0.01)
	
	rot = rot3(0.0,0.0,-pi/2.,10.)
	
	w = 0.25
	cmagf = canvas.canvas()
	for o in unifrange(-w,w,20):
		pf = fieldPath(fmap,-2.6,2.6,o,w,400)
		pf.apply(rot)
		#if len(pf.pts) < 10:
		#	continue
		pf.finish()
		#pf.sty=[style.linewidth.thin,rgb.blue]
		pf.sty=[style.linewidth.thin]	# field line color/style
		pf.drawto(cmagf)
	cmagf.writetofile("/Users/michael/Desktop/Bfield.pdf")



	
def larmor_clipping_plot():
	
	gSurv=graph.graphxy(width=20,height=10,
				x=graph.axis.lin(title="Source offset [mm]"),
				y=graph.axis.lin(title="",min=0,max=1),
				key = graph.key.key(pos="bl"))
	gSurv.texrunner.set(lfs='foils17pt')
	
	rho = 1.5
	h0 = 9.5
	
	gdat = [ [h0-h,survival_fraction(h,rho,2*3.3),survival_fraction(h,rho,2*3.3/2)] for h in unifrange(h0-10,h0,100) ]
	gdat = [ g+[0.5*(g[2]<=1e-3)+(g[2]>1e-3)*(g[1]/(g[2]+1e-6)),] for g in gdat]
	gSurv.plot(graph.data.points(gdat,x=1,y=3,title="500keV line survival"),[graph.style.line([style.linewidth.Thick,rgb.blue])])
	gSurv.plot(graph.data.points(gdat,x=1,y=2,title="1MeV line survival"),[graph.style.line([style.linewidth.Thick,rgb.red])])
	gSurv.plot(graph.data.points(gdat,x=1,y=4,title="1MeV:500keV survival ratio"),[graph.style.line([style.linewidth.Thick])])
	
	gSurv.writetofile("/Users/michael/Desktop/survival_%g.pdf"%rho)
	
	
	

def radial_clipping_plot():
	
	gSurv=graph.graphxy(width=20,height=10,
				x=graph.axis.lin(title="Source spot radius [mm]",min=0,max=9.5),
				y=graph.axis.lin(title="",min=0,max=1),
				key = graph.key.key(pos="bl"))
	gSurv.texrunner.set(lfs='foils17pt')
	
	h = 9.5
	
	gdat = [ [rho,radial_survival_fraction(h,rho,3.3),radial_survival_fraction(h,rho,3.3/2.0)] for rho in unifrange(0.,9.5,200) ]
	gdat = [ g+[0.5*(g[2]<=1e-3)+(g[2]>1e-3)*(g[1]/(g[2]+1e-6)),] for g in gdat]
	gSurv.plot(graph.data.points(gdat,x=1,y=3,title="500keV line survival"),[graph.style.line([style.linewidth.Thick,rgb.blue])])
	gSurv.plot(graph.data.points(gdat,x=1,y=2,title="1MeV line survival"),[graph.style.line([style.linewidth.Thick,rgb.red])])
	gSurv.plot(graph.data.points(gdat,x=1,y=4,title="1MeV:500keV survival ratio"),[graph.style.line([style.linewidth.Thick])])
		
	gSurv.writetofile("/Users/michael/Desktop/survival_radial.pdf")
	
	
	
if __name__ == "__main__":
	#larmor_clipping_plot()
	#radial_clipping_plot()
	#plot_larmor_trajectory()
	plot_spectrometer_field()