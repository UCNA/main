from EncalDB import *
import string

# mysql> describe analysis_results;
#	+----------------------+--------------------------------------------+------+-----+---------+----------------+
#	| Field                | Type                                       | Null | Key | Default | Extra          |
#	+----------------------+--------------------------------------------+------+-----+---------+----------------+
#	| analysis_results_id | bigint(20) unsigned                        | NO   | PRI | NULL    | auto_increment |
#	| author               | enum('MPM','RWP','BZ')                     | YES  | MUL | NULL    |                |
#	| date                 | datetime                                   | YES  | MUL | NULL    |                |
#	| type                 | enum('Asymmetry','Counts')                 | YES  | MUL | NULL    |                |
#	| source               | enum('Data','G4','Pen')                    | YES  |     | NULL    |                |
#	| start_run            | int(10) unsigned                           | YES  | MUL | NULL    |                |
#	| end_run              | int(10) unsigned                           | YES  |     | NULL    |                |
#	| event_type           | set('0','I','II','III')                    | YES  |     | NULL    |                |
#	| ana_choice           | enum('A','B','C','D','E','F','G','H')      | YES  |     | NULL    |                |
#	| side                 | enum('East','West','Both')                 | YES  |     | NULL    |                |
#	| afp                  | enum('On','Off','Other','On2Off','Off2On') | YES  |     | NULL    |                |
#	| value                | double                                     | YES  |     | NULL    |                |
#	| err                  | double                                     | YES  |     | NULL    |                |
#	| cut_spec_id          | int(10) unsigned                           | YES  | MUL | NULL    |                |
#	| grouping             | enum('run','pair','quartet','octet')       | YES  | MUL | NULL    |                |
#	| gate_valve           | enum('Open','Closed')                      | YES  | MUL | NULL    |                |
#	+----------------------+--------------------------------------------+------+-----+---------+----------------+
#
# mysql> describe cut_spec;
#	+-------------+-------------------------+------+-----+---------+----------------+
#	| Field       | Type                    | Null | Key | Default | Extra          |
#	+-------------+-------------------------+------+-----+---------+----------------+
#	| cut_spec_id | bigint(20) unsigned     | NO   | PRI | NULL    | auto_increment |
#	| energy_min  | double                  | YES  | MUL | NULL    |                |
#	| energy_max  | double                  | YES  |     | NULL    |                |
#	| radius      | double                  | YES  |     | NULL    |                |
#	| positioning | enum('plain','rotated') | YES  |     | NULL    |                |
#	+-------------+-------------------------+------+-----+---------+----------------+


analysis_results_table_fields = [
									"analysis_results_id",
									"author",
									"date",
									"type",
									"source",
									"start_run",
									"end_run",
									"event_type",
									"ana_choice",
									"side",
									"afp",
									"value",
									"err",
									"grouping",
									"gate_valve",
									"cut_spec_id" ]
cut_spec_table_fields = [
							"cut_spec_id",
							"energy_min",
							"energy_max",
							"radius",
							"positioning" ]


def open_anadb_connection():
	#return open_connection(usr="ucn",db="analysis_results_offic_2010",pw=os.environ["UCNADBPASS"][:-2])
	return open_connection(usr="ucn",db="analysis_results",pw=os.environ["UCNADBPASS"][:-2])

class anaDbResult:
	"""AnalysisDB entry"""
	def __init__(self,r,fields):
		self.dat = dict([(f,r[n]) for (n,f) in enumerate(fields)])
		for f in self.dat:
			self.__dict__[f] = self.dat[f]
		self.rrange = (self.start_run,self.end_run)
		self.erange = (self.energy_min,self.energy_max)
		
	def __repr__(self):
		return "(%i-%i %g~%g)"%(self.start_run,self.end_run,self.value,self.err)
		
	def display(self):
		print self.rrange,self.date,self.author,self.source,self.grouping,self.type,"%i-%ikeV"%self.erange,self.event_type,self.side,self.afp,self.gate_valve,"%.3f~%.3f"%(self.value,self.err)

class AnaDBLocator:
	"""Class for constructing queries to find AnalysisDB results"""
	
	def __init__(self):
		self.fields = analysis_results_table_fields[:-1]+cut_spec_table_fields[1:]
		self.req = {"source":"Data", "author":"MPM", "type":"Asymmetry", "gate_valve":"Open", "ana_choice":"C", "radius":50, "energy_min":230, "energy_max":660}
		
	def find(self,conn):
		cmd = "SELECT "+string.join(self.fields,",")+" FROM analysis_results,cut_spec WHERE cut_spec.cut_spec_id = analysis_results.cut_spec_id"
		for k in self.req:
			cmd += " AND %s='%s'"%(k,str(self.req[k]))
		print cmd
		conn.execute(cmd)
		res = [anaDbResult(r,self.fields) for r in conn.fetchall()]
		print "Located",len(res),"entries"
		return res


