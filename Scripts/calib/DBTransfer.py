#!/usr/bin/python

import sys
sys.path.append("..")
from ucnacore.EncalDB import *
from string import join


def dbInputStr(i):
	if i is None or i=='':
		return "Null"
	return "'%s'"%str(i)

class DBTransferManager:
	"""Class for transferring calibration information from one version of DB to another"""
	
	def __init__(self, connIn, connOut):
		self.connIn = connIn
		self.connOut = connOut
		self.tblcols = {}

	def transfer_lines(self, tbl, qry, xargs = {}, subgraphs = [], withPrimary = False):
		
		# identify table field names
		if tbl not in self.tblcols:
			cmd = "SHOW COLUMNS FROM `%s`"%tbl
			if not withPrimary:
				cmd += " WHERE `Key`!='PRI'"
			self.connIn.execute(cmd)
			s1 = frozenset([f[0] for f in self.connIn.fetchall()])
			if self.connOut:
				self.connOut.execute(cmd)
				s2 = frozenset([f[0] for f in self.connOut.fetchall()])
			else:
				s2 = s1
			self.tblcols[tbl] = s1.intersection(s2)
			print "Table",tbl,"has columns",self.tblcols[tbl]
		fields = [f for f in self.tblcols[tbl] if f not in xargs]

		# select items to transfer
		datlist = []
		if fields:
			cmd = "SELECT " + join(["`%s`"%str(i) for i in fields],",") + " FROM `" + tbl + "` WHERE " + qry
			print cmd
			self.connIn.execute(cmd);
			datlist = self.connIn.fetchall()
		else:
			cmd = "SELECT COUNT(*) FROM `" + tbl + "` WHERE " + qry
			print cmd
			self.connIn.execute(cmd);
			datlist = [[] for i in range(self.connIn.fetchone()[0])]
		
		if not datlist:
			return None
		
		# copy over to other DB
		cmd = None
		for z in datlist:
		
			fdict = dict([(f,z[n]) for (n,f) in enumerate(fields)])
			
			# move and re-number any graph-type subobjects
			for g in subgraphs:
				fdict[g] = self.transfer_graph(fdict[g])
			
			fdict.update(xargs)
			
			if not cmd:
				cmd = "INSERT INTO `" + tbl + "` (" + join(["`%s`"%i for i in fdict.keys()],',') + ") VALUES"
			cmd += "\n\t(" + join([dbInputStr(i) for i in fdict.values()],',') + "),"
		
		cmd = cmd[:-1]
		print cmd
		if self.connOut is not None:
			self.connOut.execute(cmd)
			self.connOut.execute("SELECT LAST_INSERT_ID()")
			return self.connOut.fetchone()[0]

		return None

	def transfer_graph(self,gid):
		print
		print "Relocating old graph",gid
		iid = self.transfer_lines("graphs","graph_id = %i"%gid)
		self.transfer_lines("graph_points","graph_id = %i"%gid, {"graph_id":iid})
		return iid

	def transfer_posmap(self,pmid):
		iid = self.transfer_lines("posmap_set", "posmap_set_id = %i"%pmid)
		if not iid:
			print "No such position map!"
			return None
		self.transfer_lines("posmap_points", "posmap_set_id = %i"%pmid, {"posmap_set_id":iid})
		return iid

	def transfer_energy_calibration(self,ecid,pmid):
		iid = self.transfer_lines("energy_calibration", "ecal_id = %i"%ecid, {"posmap_set_id":pmid})
		self.transfer_lines("tube_calibration", "ecal_id = %i"%ecid, {"ecal_id":iid}, subgraphs=["linearity_graph"])

	def transfer_mwpc_ecal(self,ecid,pmid):
		return self.transfer_lines("mwpc_ecal", "mwpc_ecal_id = %i"%ecid, {"gain_posmap_id":pmid})

	def transfer_evis_conversion(self,evcid):
		return self.transfer_lines("evis_conversion", "evis_conversion_id = %i"%evcid, subgraphs=["conversion_curve_id"])

	def transfer_cathscale(self,csid):
		return self.transfer_lines("cath_ccloud_scale", "cath_ccloud_scale_id = %i"%csid, subgraphs=["gain_graph_id"])

	def transfer_cathcal_set(self, csid, transferShape=False):
		ncsid = self.transfer_lines("cathcal_set", "cathcal_set_id = %i"%csid)
		self.connIn.execute("SELECT cathseg_id FROM cathseg_cal WHERE cathcal_set_id=%i"%csid)
		ocsgids = [l[0] for l in self.connIn.fetchall()]
		print "Transferring cathseg calibrations",ocsgids
		for ocsgid in ocsgids:
			ncsgid = self.transfer_lines("cathseg_cal", "cathseg_id = %i"%ocsgid, {"cathcal_set_id":ncsid})
			if transferShape:
				self.transfer_lines("cathshape_graphs", "cathseg_id = %i"%ocsgid, {"cathseg_id":ncsgid}, subgraphs=["graph_id"], withPrimary=True)
		return ncsid

	def transfer_runlist(self,rmin,rmax):
		self.transfer_lines("run", "%i <= run_number AND run_number <= %i"%(rmin,rmax), withPrimary=True)

	def transfer_rungroups(self):
		self.transfer_lines("run_group","1")




def delete_all_entries(conn):
	"""Print commands to empty all tables. Be careful!!"""
	conn.execute("SHOW TABLES")
	for tbl in conn.fetchall():
		print "TRUNCATE TABLE `%s`;"%tbl[0]



if __name__ == "__main__":

	connIn = open_connection(db="mpm_debug")
	connOut = open_connection(db="cal_starter_DB")
	dbtm = DBTransferManager(connIn, connOut)
	
	#delete_all_entries(connOut); exit(0)
	
	if 1:
		pmid = dbtm.transfer_posmap(213)
		dbtm.transfer_energy_calibration(9726,pmid)
		if connOut is not None:
			connOut.execute("UPDATE energy_calibration SET start_run=13000, end_run=100000")
		
		pmid = dbtm.transfer_posmap(181)
		dbtm.transfer_mwpc_ecal(939,pmid)
		dbtm.transfer_mwpc_ecal(940,pmid)
		
		for i in range(79,87):
			dbtm.transfer_evis_conversion(i)
		for i in range(21,25):
			dbtm.transfer_cathscale(i)
		for i in range(198,202):
			dbtm.transfer_cathcal_set(i)

		dbtm.transfer_runlist(16500,100000)
		dbtm.transfer_rungroups()

####
# mysqldump --skip-lock-tables -u root -p cal_starter_DB > cal_starter_DB_v1.sql

