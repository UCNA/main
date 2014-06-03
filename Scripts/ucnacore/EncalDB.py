
import os
import MySQLdb
import getpass
from string import join

# open connection to calibration DB
def open_connection(db=None,usr=None,pw=None,host=None):
	if not host:
		host = os.environ["UCNADBADDRESS"]
	if not db:
		db = os.environ["UCNADB"]
	if not usr:
		usr = os.environ["UCNADBUSER"]
	if not pw:
		pw = os.environ["UCNADBPASS"]
	print "Connecting to database",db,"with user",usr,"on",host
	conn = None
	nfails = 0
	while nfails < 3 and not conn:
		try:
			passwd = pw
			if not pw:
				passwd = getpass.getpass("password for %s? "%usr)
			conn = MySQLdb.connect(host = host, user = usr, passwd = passwd, db = db)
		except:
			print "Connection failed! Try again?"
			nfails += 1
	if not conn:
			passwd = pw
			if not pw:
				passwd = getpass.getpass("password for %s? "%usr)
			conn = MySQLdb.connect(host = host, user = usr, passwd = passwd, db = db)
	

	curs = conn.cursor()
	curs.execute("SELECT VERSION()")
	r=curs.fetchone()
	if r:
		print "Connected to server version:",r[0]
		conn.autocommit(True)
		return curs
	else:
		print "CONNECTION FAILED"
		exit()

# generate a new (empty) graph ID
def newgraph(conn,gdescrip):
	conn.execute("INSERT INTO graphs (text_description) VALUES ('%s')"%gdescrip)
	conn.execute("SELECT LAST_INSERT_ID()")
	return int(conn.fetchone()[0])

# upload a new graph, formatted [[x,dx,y,dy],...]; return graph ID
def upload_graph(conn,gdescrip,gdata):
	gid=newgraph(conn,gdescrip)
	print "Loading new graph id",gid
	if len(gdata[0]) >= 4:
		conn.execute("INSERT INTO graph_points (graph_id, x_value, x_error, y_value, y_error) VALUES " + join(["(%i,%f,%f,%f,%f)"%(gid,d[0],d[1],d[2],d[3]) for d in gdata],','))
	elif len(gdata[0]) >= 2:
		conn.execute("INSERT INTO graph_points (graph_id, x_value, y_value) VALUES " + join(["(%i,%f,%f)"%(gid,d[0],d[1]) for d in gdata],','))
	else:
		print "Bad dimensions for data array!"
		exit(-1)
	return gid

# delete graph by ID
def delete_graph(conn,gid):
	print "Deleting graph",gid
	conn.execute("DELETE FROM graph_points WHERE graph_id = %i"%gid)
	conn.execute("DELETE FROM graphs WHERE graph_id = %i"%gid)

# get graph data by ID
def getGraph(conn,gid):
	conn.execute("SELECT x_value,x_error,y_value,y_error FROM graph_points WHERE graph_id = %i"%gid)
	return conn.fetchall()

# get beta runs in range
def getRunType(conn,rtype="Asymmetry",rmin=0,rmax=1e6):
	rlist = []
	if type(rtype) != type([]):
		rtype = [rtype]
	for t in rtype:
		conn.execute("SELECT run_number FROM run WHERE run_type = '%s' AND %i <= run_number AND run_number <= %i"%(t,rmin,rmax))
		rlist += [r[0] for r in conn.fetchall()]
	rlist.sort()
	return rlist

# get start time for run
def getRunStartTime(conn,rn):
	conn.execute("SELECT UNIX_TIMESTAMP(start_time) FROM run WHERE run_number = %i"%rn)
	st = conn.fetchone()
	if not st:
		print "** unknown start time for",rn
		return 0
	return st[0]

# get end time for run
def getRunEndTime(conn,rn):
	conn.execute("SELECT UNIX_TIMESTAMP(end_time) FROM run WHERE run_number = %i"%rn)
	st = conn.fetchone()
	if not st:
		print "** unknown end time for",rn
		return 0
	return st[0]

# get table of all run start/end times
def getRunTimeTable(conn):
	conn.execute("SELECT run_number,UNIX_TIMESTAMP(start_time),UNIX_TIMESTAMP(end_time) FROM run")
	return dict([(r[0],(r[1],r[2])) for r in conn.fetchall()])

# get live time from analyzer
def getRunLiveTime(conn,rn):
	conn.execute("SELECT live_time FROM analysis WHERE run_number = %i"%rn)
	st = conn.fetchone()
	if not st:
		print "** unknown live time for",rn
		return 0
	return st[0]

# get list of all GMS runs
def getGMSruns(conn):
	conn.execute("SELECT gms_run FROM energy_calibration WHERE 1 ORDER BY gms_run ASC")
	return [ r[0] for r in conn.fetchall() ]

# get run start times for GMS runs
def gmsRunTimes(conn):
	conn.execute("SELECT run_number,UNIX_TIMESTAMP(start_time) FROM run WHERE run_number IN (SELECT gms_run FROM energy_calibration WHERE 1 ORDER BY gms_run ASC)")
	return conn.fetchall()
	
# get sensor ID number for name
sensorIDtable = {}
def getSensorID(conn,sname):
	if sname in sensorIDtable:
		return sensorIDtable[sname]
	conn.execute("SELECT sensor_id FROM sensors WHERE sensor_name = '%s'"%sname)
	sid = conn.fetchone()[0]
	sensorIDtable[sname] = sid
	return sid
	
# get graph IDs for run monitor
def getRunMonitorGIDs(conn,rn,sname,tp):
	if type(sname) != type(123):
		sname = getSensorID(conn,sname)
	conn.execute("SELECT center_graph_id,width_graph_id FROM run_monitors WHERE run_number = %i AND sensors_sensor_id = %i AND monitor_type = '%s'"%(rn,sname,tp))
	return conn.fetchone()

# get trigger efficiency params
def getTrigeffParams(conn, rn, s, t):
	if s in ['E','W']:
		s = {'E':"East",'W':"West"}[s]
	conn.execute("SELECT params_graph FROM mpm_trigeff WHERE run_number = %i AND side = '%s' AND quadrant = %i"%(rn,s,t))
	return getGraph(conn,conn.fetchone()[0])

# delete run monitor by ID
def deleteRunMonitor(conn,rmid):
	print "Deleting run monitor",rmid
	conn.execute("SELECT center_graph_id,width_graph_id FROM run_monitors WHERE monitor_id = %i"%rmid)
	for r in conn.fetchall():
		delete_graph(conn,r[0])
		delete_graph(conn,r[1])
	conn.execute("DELETE FROM run_monitors WHERE monitor_id = %i"%rmid)

# delete all run monitors in range
def deleteRunMonitors(conn,rmin,rmax,mtype='pedestal'):
	conn.execute("SELECT monitor_id FROM run_monitors WHERE %i <= run_number AND run_number <= %i and monitor_type = '%s'"%(rmin,rmax,mtype))
	for r in conn.fetchall():
		deleteRunMonitor(conn,r[0])


################
# Position maps
################

class posmap_point:
	pass
class posmap_info:
	def __repr__(self):
		return "[Posmap %i: '%s' %i, %g]"%(self.posmap_set_id,self.descrip,self.n_rings,self.radius)

# one position map
class posmap:
	def __init__(self,info):
		self.pts = {}
		self.info = info
	def add_pt(self,pt):
		self.pts[pt.n] = pt
	def get_pts_sorted(self):
		l = self.pts.keys()
		l.sort()
		return [self.pts[k] for k in l]
	def get_pt_vals(self):
		return [p.sig/p.norm for p in self.get_pts_sorted()]
	def avg_val(self):
		return sum(self.get_pt_vals())/len(self.pts)

# get information about a posmap
def getPosmapInfo(conn,pmid):
	conn.execute("SELECT descrip,n_rings,radius FROM posmap_set WHERE posmap_set_id=%i"%pmid)
	p = conn.fetchone()
	info = posmap_info()
	info.posmap_set_id = pmid
	info.descrip = p[0]
	info.n_rings = p[1]
	info.radius = p[2]
	return info

# get set of all position maps corresponding to given ID number
def getPosmapSet(conn,pmid):
	pinfo = getPosmapInfo(conn,pmid)
	conn.execute("SELECT side,quadrant,pixel_id,`signal`,norm,center_x,center_y FROM posmap_points WHERE posmap_set_id = %i"%pmid)
	pts = {}
	for p in conn.fetchall():
		x = posmap_point()
		x.n = p[2]
		x.sig = p[3]
		x.norm = p[4]
		x.x = p[5]
		x.y = p[6]
		pts.setdefault((p[0],p[1]),posmap(pinfo)).add_pt(x)
	for k in pts:
		pts[k].side = k[0]
		pts[k].quadrant = k[1]
	return pts

# generate new position map ID from information (assigning new posmap set ID to info)
def newPosmap(conn,pinfo):
	conn.execute("INSERT INTO posmap_set (descrip,n_rings,radius) VALUES ('%s',%i,%g)"%(pinfo.descrip,pinfo.n_rings,pinfo.radius))
	conn.execute("SELECT LAST_INSERT_ID()")
	pinfo.posmap_set_id = int(conn.fetchone()[0])
	print "Uploaded",pinfo

# upload posmap points
def uploadPosmap(conn,pmap):
	print "Uploading points for",pmap.info
	cmd = "INSERT INTO posmap_points (posmap_set_id, side, quadrant, pixel_id, center_x, center_y, `signal`, norm) VALUES "
	cmd += join(["\n\t(%i,'%s',%i,%i,%g,%g,%g,%g)"%(pmap.info.posmap_set_id,pmap.side,pmap.quadrant, p.n,p.x,p.y,p.sig,p.norm) for p in pmap.get_pts_sorted()],',')
	print cmd
	conn.execute(cmd)
		