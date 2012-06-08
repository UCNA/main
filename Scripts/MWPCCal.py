#!/usr/bin/python

from EncalDB import *
from math import *

class wireSpec:
	def __init__(self,pos,nm):
		self.pos = pos
		self.nm = nm
		self.norm = 1.0


###############
# DB commands
###############

def delete_shape_graphs(conn,csid):
	conn.execute("SELECT graph_id FROM cathshape_graphs WHERE cathseg_id = %i"%csid)
	for gid in [r[0] for r in conn.fetchall()]:
		delete_graph(conn,gid)
	conn.execute("DELETE FROM cathshape_graphs WHERE cathseg_id = %i"%csid)

def delete_cathcal_set(conn,ccsid):
	print "Deleting cathode calibration set",ccsid
	conn.execute("SELECT cathseg_id FROM cathseg_cal WHERE cathcal_set_id = %i"%ccsid)
	for csid in [r[0] for r in conn.fetchall()]:
		delete_shape_graphs(conn,csid)
	conn.execute("DELETE FROM cathcal_set WHERE cathcal_set_id = %i"%ccsid)

def new_cathcal_set(conn,side,plane,r0,r1):
	conn.execute("INSERT INTO cathcal_set(side,plane,start_run,end_run) VALUES ('%s','%s',%i,%i)"%(side,plane,r0,r1))
	conn.execute("SELECT LAST_INSERT_ID()")
	ccsid = int(conn.fetchone()[0])
	print "Generating cathcal set for",side,plane,r0,r1,":",ccsid
	return ccsid

def find_cathcal_sets(conn,r0,r1):
	conn.execute("SELECT cathcal_set_id FROM cathcal_set WHERE start_run = %i AND end_run = %i"%(r0,r1))
	return [r[0] for r in conn.fetchall()]

def new_cathseg_cal(conn,ccsid,wspec):
	conn.execute("INSERT INTO cathseg_cal(cathcal_set_id,position,sensor_name,norm) VALUES (%i,%f,'%s',%f)"
				 %(ccsid,wspec.pos,wspec.nm,wspec.norm))
	conn.execute("SELECT LAST_INSERT_ID()")
	return int(conn.fetchone()[0])

def set_cathshape_graph(conn,csid,gid):
	conn.execute("INSERT INTO cathshape_graphs(graph_id,cathseg_id) VALUES (%i,%i)"%(gid,csid))

###############
# wires info
###############

def getWires(rn,s,d):
	"""Get list of wires active for run/side/plane"""	
	assert s in ["East","West"] and d in ["X","Y"] and rn >= 13000
	nWires = 16
	wireSpacing = 4*2.54*sqrt(0.6)
	#								0    1    2    3    4    5    6    7    8    9    10   11   12   13   14   15
	padc_nums = {
					("East","X"): ( 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231 ),
					("East","Y"): ( 20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  210, 211, 212, 213, 214, 215 ),
					("West","X"): ( 31,  30,  29,  28,  27,  26,  25,  24,  23,  22,  21,  20,  19,  18,  17,  16 ),
					("West","Y"): ( 0,   1,   2,   3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  13,  14,  15 ) }
	wires = []
	for i in range(nWires):
		pos = (nWires*0.5-0.5-i)*wireSpacing
		nm = "Pdc%i"%padc_nums[(s,d)][i]
		if s=="West":
			nm = "Padc%i"%padc_nums[(s,d)][i]
		wires.append(wireSpec(pos,nm))
	return wires

def gen_cathcal_set(conn,r0,r1):
	"""Plain cathcal set with no shape corrections"""
	for ccsid in find_cathcal_sets(conn,r0,r1):
		delete_cathcal_set(conn,ccsid)
	for s in ["East","West"]:
		for d in ["X","Y"]:
			ccsid = new_cathcal_set(conn,s,d,r0,r1)
			for w in getWires(r0,s,d):
				new_cathseg_cal(conn,ccsid,w)


###############
#             #
###############

if __name__ == "__main__":
	conn = open_connection()
	gen_cathcal_set(conn,13000,100000)
