#!/usr/bin/python

import os
import time
from EncalDB import *

class RunInfo:
	def __init__(self,rn,rtype="Other",ccyc="Unknown"):
		self.runNum = rn			# run number
		self.runType = rtype		# A1...B12
		self.run_class = "Other"	# generic class, Asym/Depol/etc
		self.calcycle = ccyc		# cal cycle name
		self.slowdaq = 0			# slowdaq number
		self.ledv = [0,0]
		self.scs = 0.0
		self.startTime = 0			# start timestamp
		self.endTime = 0			# end timestamp
		self.comments = ""			# comments from my runlog
		self.title = ""				# title given when run
		self.gate_valve="Other"
		self.flipper="Other"
	
	def display(self):
		print self.runNum,self.calcycle,"'%s'"%self.runType,self.ledv,self.scs,self.comments

	# set run parameters from run type
	def setType(self):
		if self.runType[0] in ('A','B') and self.runType[1:].isdigit():

			if self.runType in ('A1','A12','B4','B9'):
				self.run_class = 'Asymmetry'
				self.gate_valve = 'Closed'
				self.flipper = 'Off'
			
			elif self.runType in ('A2','A10','B5','B7'):
				self.run_class = 'Asymmetry'
				self.gate_valve = 'Open'
				self.flipper = 'Off'
			
			elif self.runType in ('A4','A9','B1','B12'):
				self.run_class = 'Asymmetry'
				self.gate_valve = 'Closed'
				self.flipper = 'On'
			
			elif self.runType in ('A5','A7','B2','B10'):
				self.run_class = 'Asymmetry'
				self.gate_valve = 'Open'
				self.flipper = 'On'

			elif self.runType in ('A3','A11','B6','B8'):
				self.run_class = 'Depol'
				self.gate_valve = 'Other'
				self.flipper = 'Off2On'
				
			elif self.runType in ('A6','A8','B3','B11'):
				self.run_class = 'Depol'
				self.gate_valve = 'Other'
				self.flipper = 'On2Off'
			
		else:
		
			self.run_class = "Other"
			if self.runType == 'LEDCal':
				self.run_class = 'LEDCalib'
			elif self.runType == 'SourcesCal':
				self.run_class = 'SourceCalib'
			elif self.runType == 'BgOn':
				self.gate_valve = 'Closed'
				self.flipper = 'On'
			elif self.runType == 'BgOff':
				self.gate_valve = 'Closed'
				self.flipper = 'Off'
			elif self.runType == 'BetaOn':
				self.gate_valve = 'Open'
				self.flipper = 'On'
			elif self.runType == 'BetaOff':
				self.gate_valve = 'Open'
				self.flipper = 'Off'
			elif self.runType == 'Xe':
				self.run_class = 'Xenon'
			else:
				print "Unknown run type"
				self.display()

		
		
		
#read in auto generated runlog.txt
def load_runlog(fname="../SummaryData/runlog_2011.txt"):
	runs={}
	#0          1       2           3             4     5           6   7   8   9        10      11      12
	#Tue Aug 17	13552	20:08:10	20:08:34	  2.2k	daq test	0	0	0	open	 0.00T	 0.00T	 0.00T
	for l in [l.split('\t') for l in open(fname,"r").readlines() ]:
		rn = int(l[1])
		R = RunInfo(rn)
		R.startTime = time.mktime(time.strptime("%s %s 2010"%(l[0],l[2]),"%a %b %d %H:%M:%S %Y"))
		R.endTime = time.mktime(time.strptime("%s %s 2010"%(l[0],l[3]),"%a %b %d %H:%M:%S %Y"))
		if R.endTime < R.startTime:
			R.endTime += 24*3600
		R.title = l[5]
		R.gate_valve = "Closed"
		if int(l[6]):
			R.gate_valve = "Open"
		R.flipper = "Off"
		if int(l[7]):
			R.flipper = "On"
		R.scs_field = float(l[10][:-1])
		R.afp_field = float(l[11][:-1])
		R.ppm_field = float(l[12][:-2])
		runs[rn] = R
	return runs
		
#read in my UCNA run log		
def load_mylog(fname="../Aux/UCNA Run Log.txt"):
	f = [l.split()+[l,] for l in open(fname,"r").readlines() if l[0] in ['*','@']]
	runs = {}
	ledv = [0,0]
	calcycle = 0
	scs = 0
	cname = "Unknown"
	sources={'E':[],'W':[]}
	
	for l in f:
		if l[0] == '@cal':
			calcycle += 1
		if l[0] == '@LED':
			#ledv = [float(l[1]),float(l[2])]
			pass
		if l[0] == '@scs':
			scs = float(l[1])
		if l[0] == '@cycle':
			cname = l[-1][6:].strip()
		if l[0] == '@sources':
			pass
		if l[0][0] == '*':
			rn = int(l[0][1:])
			rt = l[1]
			R = RunInfo(rn,rt,calcycle)
			R.ledv = ledv
			R.scs = scs
			R.cname = cname
			if rt[0] in ['A','B'] and l[2].isdigit():
				R.slowdaq = int(l[2])
			
			#if rt == 'LEDCal':
			#	if len(l) > 3:
			#		R.ledv = [float(l[2]),float(l[3])]
				
			if l[-1].find('#')>0:
				R.comments = l[-1][l[-1].find('#')+1:].strip()
			if rn in runs:
				print "******** WARNING: Duplicate Run",rn,"in Run Log!"
			runs[rn] = R
				
	return runs




# merge auto and custom runlogs
def merge_runlogs(mylog,autolog):
	runs = {}
	for r in autolog:
		R = autolog[r]
		if r in mylog:
			R.comments = mylog[r].comments
			R.calcycle = mylog[r].calcycle
			R.runType = mylog[r].runType
			R.setType()
			
			runs[r] = R
			
	for r in mylog:
		if r not in autolog:
			print "*** Missing run",r,"from Midas log!"
					
	return runs


# update runs DB with run list
def fillRunsDB(runs,rmin=0,rmax=100000):

		conn = open_connection()
		
		for rn in runs:
			if not rmin <= rn <= rmax:
				continue
			R = runs[rn]
			R.comments = R.comments.replace("'","\\'")
			R.title = R.title.replace("'","\\'")
			
			rdata = {}
			rdata["run_number"]=R.runNum
			rdata["slow_run_number"]=R.slowdaq
			rdata["asym_oct"]="'%s'"%R.runType
			rdata["run_type"]="'%s'"%R.run_class
			rdata["gate_valve"]="'%s'"%R.gate_valve
			rdata["flipper"]="'%s'"%R.flipper
			rdata["scs_field"]=R.scs
			rdata["comments"]="'%s'"%R.comments
			rdata["title"]="'%s'"%R.title
			rdata["geometry"]="'%s'"%"D"
			rdata["start_time"]=time.strftime("'%Y-%m-%d %H:%M:%S'",time.localtime(R.startTime))
			rdata["end_time"]=time.strftime("'%Y-%m-%d %H:%M:%S'",time.localtime(R.endTime))
			
			cmd = "SELECT COUNT(*) FROM run WHERE run_number=%i"%rn
			conn.execute(cmd)
	
			if conn.fetchone()[0]:
				cmd = "UPDATE run SET"
				for k in rdata.keys():
					cmd += " %s=%s,"%(k,str(rdata[k]).strip())
				cmd = cmd[:-1]+" WHERE run_number = %i"%R.runNum				
			else:
				cmd = "INSERT INTO run ("
				fields = rdata.keys()
				for f in fields:
					cmd += f+","
				cmd = cmd[:-1]+")\n\tVALUES ("
				for f in fields:
					cmd += "'%s',"%str(rdata[f]).strip()
				cmd = cmd[:-1]+")"
				
			print cmd
			conn.execute (cmd)

				
def fillRunGroups():
	conn = open_connection()
	conn.execute("DELETE FROM run_group WHERE 1")
	cycs = {}
	runs = load_mylog()
	for rn in runs:
		cycs.setdefault(runs[rn].cname,[]).append(runs[rn].runNum)
	for c in cycs:
		cycs[c].sort()
		print c,cycs[c]
		conn.execute("INSERT INTO run_group(start_run,end_run,name) VALUES (%i,%i,'%s')"%(cycs[c][0],cycs[c][-1],c))
	
if __name__=="__main__":   
	al = load_runlog("/data/ucnadata/midfiles_2010/runlog.txt")
	ml = load_mylog()
	runs = merge_runlogs(ml,al)
	fillRunsDB(runs,rmin=13500,rmax=16300)
	fillRunGroups()
	