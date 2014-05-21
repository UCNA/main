from EncalDB import *
from os import *
import string

# SELECT DISTINCT name FROM analysis_numbers;

# mysql> describe analysis_runset;
#	+--------------------+------------------------------------------------------+------+-----+---------+----------------+
#	| Field              | Type                                                 | Null | Key | Default | Extra          |
#	+--------------------+------------------------------------------------------+------+-----+---------+----------------+
#	| analysis_runset_id | bigint(20) unsigned                                  | NO   | PRI | NULL    | auto_increment |
#	| start_run          | int(10) unsigned                                     | YES  | MUL | NULL    |                |
#	| end_run            | int(10) unsigned                                     | YES  |     | NULL    |                |
#	| grouping           | enum('run','fgbg','ppair','quartet','octet','range') | YES  |     | NULL    |                |
#	| gate_valve         | enum('Open','Closed','Other')                        | YES  | MUL | NULL    |                |
#	| afp                | enum('On','Off','Other','On2Off','Off2On')           | YES  |     | NULL    |                |
#	+--------------------+------------------------------------------------------+------+-----+---------+----------------+

# mysql> describe analysis_numbers;
#	+--------------------+-----------------------------------+------+-----+---------+----------------+
#	| Field              | Type                              | Null | Key | Default | Extra          |
#	+--------------------+-----------------------------------+------+-----+---------+----------------+
#	| analysis_number_id | bigint(20) unsigned               | NO   | PRI | NULL    | auto_increment |
#	| analysis_runset_id | int(10) unsigned                  | YES  | MUL | NULL    |                |
#	| source             | varchar(64)                       | YES  | MUL | NULL    |                |
#	| name               | varchar(128)                      | YES  |     | NULL    |                |
#	| date               | datetime                          | YES  |     | NULL    |                |
#	| side               | enum('East','West','Both','None') | YES  |     | NULL    |                |
#	| event_type         | set('0','I','II','III')           | YES  |     | NULL    |                |
#	| n                  | int(11)                           | YES  |     | NULL    |                |
#	| value              | double                            | YES  |     | NULL    |                |
#	| err                | double                            | YES  |     | NULL    |                |
#	+--------------------+-----------------------------------+------+-----+---------+----------------+



analysis_runset_table_fields = [
								"analysis_runset_id",
								"start_run",
								"end_run",
								"grouping",
								"gate_valve",
								"afp" ]
analysis_numbers_table_fields = [
								"analysis_number_id",
								"analysis_runset_id",
								"source",
								"name",
								"date",
								"side",
								"event_type",
								"n",
								"value",
								"err" ]


def open_anadb_connection():
	return open_connection(usr="ucn",db="analysis_results",pw=os.environ["UCNADBPASS"][:-2])

class anaDbResult:
	"""AnalysisDB entry"""
	def __init__(self,r,fields):
		self.dat = dict([(f,r[n]) for (n,f) in enumerate(fields)])
		for f in self.dat:
			self.__dict__[f] = self.dat[f]
		self.rrange = (self.start_run,self.end_run)

	def __repr__(self):
		return "(%i-%i %g~%g)"%(self.start_run,self.end_run,self.value,self.err)
		
	def display(self):
		print self.rrange,self.date,self.author,self.source,self.grouping,self.type,self.event_type,self.side,self.afp,self.gate_valve,"%.3f~%.3f"%(self.value,self.err)

class AnaDBLocator:
	"""Class for constructing queries to find AnalysisDB results"""
	
	def __init__(self):
		self.fields = analysis_numbers_table_fields[2:] + analysis_runset_table_fields[1:]
		self.req = {"source":os.environ["UCNA_ANA_AUTHOR"]+"_Data"}
		self.xcond = None
		
	def find(self,conn):
		cmd = "SELECT "+string.join(self.fields,",")+" FROM analysis_numbers,analysis_runset WHERE analysis_runset.analysis_runset_id = analysis_numbers.analysis_runset_id"
		for k in self.req:
			cmd += " AND %s='%s'"%(k,str(self.req[k]))
		if self.xcond:
			cmd += " AND " + self.xcond
		print cmd
		conn.execute(cmd)
		res = [anaDbResult(r,self.fields) for r in conn.fetchall()]
		print "Located",len(res),"entries"
		return res
