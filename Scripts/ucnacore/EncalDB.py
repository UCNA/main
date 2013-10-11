#!/sw/bin/python2.6
import os
import MySQLdb
import getpass

# open write connection to calibration DB
def open_connection():
	host = os.environ["UCNADBADDRESS"]
	db = os.environ["UCNADB"]
	usr = os.environ["UCNADBUSER"]
	pw = os.environ["UCNADBPASS"]
	print "Connecting to database",db,"on",host
	conn = None
	nfails = 0
	while nfails < 3 and not conn:
		try:
			passwd = pw
			if not pw:
				passwd = getpass.getpass("password for %s? "%usr)
			conn = MySQLdb.connect(host = host, user = usr, passwd = passwd, db = db).cursor()
		except:
			print "Connection failed! Try again?"
			nfails += 1
	if not conn:
			passwd = pw
			if not pw:
				passwd = getpass.getpass("password for %s? "%usr)
			conn = MySQLdb.connect(host = host, user = usr, passwd = passwd, db = db).cursor()
			
	conn.execute("SELECT VERSION()")
	r=conn.fetchone()
	if r:
		print "Connected to server version:",r[0]
		return conn
	else:
		print "CONNECTION FAILED"
		exit()

# generate a new (empty) graph ID
def newgraph(conn,gdescrip):
	conn.execute("INSERT INTO graphs (text_description) VALUES ('%s')"%gdescrip)
	conn.execute("SELECT LAST_INSERT_ID()")
	return int(conn.fetchone()[0])

# upload a new graph
def upload_graph(conn,gdescrip,gdata):
	gid=newgraph(conn,gdescrip)
	print "Loading new graph id",gid
	if len(gdata[0]) >= 4:
		conn.executemany("INSERT INTO graph_points (graph_id, x_value, x_error, y_value, y_error) VALUES (%s,%s,%s,%s,%s)",[(gid,d[0],d[1],d[2],d[3]) for d in gdata])
	elif len(gdata[0]) >= 2:
		conn.executemany("INSERT INTO graph_points (graph_id, x_value, y_value) VALUES (%s,%s,%s)",[(gid,d[0],d[1]) for d in gdata])
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
	return conn.fetchone()[0]

# get end time for run
def getRunEndTime(conn,rn):
	conn.execute("SELECT UNIX_TIMESTAMP(end_time) FROM run WHERE run_number = %i"%rn)
	return conn.fetchone()[0]

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
