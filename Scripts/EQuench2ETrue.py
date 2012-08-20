#!/usr/bin/python

from QFile import *
from LinFitter import *
from PyxUtils import *
from EncalDB import *
import os

def deleteEQ2ET(conn,eqid):
	print "Deleting Equench to Etrue",eqid
	conn.execute("SELECT conversion_curve_id FROM evis_conversion WHERE evis_conversion_id = %i"%eqid)
	for i in conn.fetchall():
		if i[0]:
			delete_graph(conn,i[0])
	conn.execute("DELETE FROM evis_conversion WHERE evis_conversion_id = %i"%eqid)

def delete_all_EQ2ET(conn):
	print "Deleting ALL Equench->Etrue conversions"
	conn.execute("SELECT evis_conversion_id FROM evis_conversion")
	for i in conn.fetchall():
		deleteEQ2ET(conn,i[0])


def uploadEQ2ET(conn,tp,LF,rmin=13000,rmax=100000):
	tp = int(tp)
	print "Loading Equench->Etrue curve for type",tp,"Runs",rmin,"-",rmax
	curvedat = [ (x,LF(x)) for x in LF.unifPoints(1,1000,400) ]	
	for s in ["East","West"]:
		lgid = upload_graph(conn,"Equench to Etrue %s Type %i"%(s,tp),curvedat)
		conn.execute("""INSERT INTO evis_conversion
						(start_run,end_run,side,type,conversion_curve_id)
						VALUES (%i,%i,'%s',%i,%i)"""
						% (rmin,rmax,s,tp,lgid))

class EnergyPoint(KVMap):
	def __init__(self,m):
		self.dat = m.dat
		self.loadFloats([k for k in self.dat.keys() if k!="side"])
		self.loadStrings(["side"])
		
def EQ2ET(fbase, conn=None):
	
	f = QFile(fbase+"/Evis2ETrue.txt")
	epts = [EnergyPoint(m) for m in f.dat.get("spectrumInfo",[])]
	
	scols = {"East":rgb.red,"West":rgb.blue}
	typeSymbs = {0:symbol.circle,1:symbol.triangle,2:symbol.square,3:symbol.cross}
	typeLines = {0:style.linestyle.solid,1:style.linestyle.dashed,2:style.linestyle.dotted,3:style.linestyle.dashdotted}
	typeNames = {0:"0",1:"I",2:"II",3:"III"}
	
	for forward in [True,]:
	
		xtitle = "Mean Observed Quenched Energy [keV]"
		ytitle = "True Initial Energy [keV]"
		if not forward:
			xtitle = "Mean True Energy [keV]"
			ytitle = "Observed Quenched Energy [keV]"

		gResid=graph.graphxy(width=7.5,height=1.0,
				x=graph.axis.lin(title=xtitle,min=0,max=800),
				y=graph.axis.lin(title="Residuals",min=-10,max=10))
		gResid.texrunner.set(lfs='10ptex')
		
		gEn=graph.graphxy(width=7.5,height=7.5,ypos=gResid.height+0.5,
				x=graph.axis.linkedaxis(gResid.axes["x"]),
				y=graph.axis.lin(title=ytitle,min=0,max=800),
				key = graph.key.key(pos="tl"))
		gEn.texrunner.set(lfs='10ptex')
		
		cnvs = canvas.canvas()
		cnvs.insert(gResid)
		cnvs.insert(gEn)
		
		types = typeSymbs.keys()
		types.sort()
		
		for tp in types:
			
			combodat = []
			
			for s in scols:
			
				gdat = []
				if forward:
					gdat = [ (p.evis_avg,p.etrue_input,p.evis_rms/sqrt(p.evis_counts)) for p in epts if p.side==s and p.type==tp and p.evis_counts > 50 ]
				else:
					gdat = [ (p.etrue_avg,p.evis_input,p.etrue_rms/sqrt(p.etrue_counts)) for p in epts if p.side==s and p.type==tp and p.etrue_counts > 50 ]
				combodat += gdat
		
				gEn.plot(graph.data.points(gdat,x=1,y=2,dx=3,title=None),
					[graph.style.symbol(typeSymbs[tp],size=0.07,symbolattrs=[scols[s],]),graph.style.errorbar(errorbarattrs=[scols[s]])])
			
			combodat.sort()		
			LF = LinearFitter(terms=[polyterm(0),polyterm(1),polyterm(-1),polyterm(-2)])
			LF.fit([c for c in combodat],cols=(0,1))
			print s,LF.toLatex()
			LF = piecewiseExtender(LF,combodat[1][0],900)
			
			if forward and conn:
				uploadEQ2ET(conn,tp,LF)
				
			gEn.plot(graph.data.points(LF.fitcurve(1,1000,400),x=1,y=2,title=None),
				[graph.style.line([typeLines[tp],])])
			
			ctitle = "Type %s"%typeNames[tp]
			#ctitle += ": $E_T = %s$"%LF.toLatex("E_Q")
			gEn.plot(graph.data.points([(-1,-1),],x=1,y=2,title=ctitle),
				[graph.style.line([typeLines[tp],]),graph.style.symbol(typeSymbs[tp],size=0.10)])
				
			gResid.plot(graph.data.points([ (g[0],g[1]-LF(g[0])) for g in combodat],x=1,y=2,title=None),
				[graph.style.symbol(typeSymbs[tp],size=0.07)])
		
		if forward:
			cnvs.writetofile(fbase+"/EvisToETrue.pdf")
		else:
			cnvs.writetofile(fbase+"/EtrueToEvis.pdf")
			
if __name__ == "__main__":
	conn = open_connection()
	#conn = None
	if conn:
		delete_all_EQ2ET(conn)
	EQ2ET(os.environ["UCNA_ANA_PLOTS"]+"/Evis2ETrue/20120810/",conn)
	