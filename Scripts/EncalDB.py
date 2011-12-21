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

# delete a posmap
def delete_posmap(conn,pmid):
	print "Deleting posmap",pmid
	conn.execute("SELECT linearity_graph FROM posmap_info WHERE posmap_set_id = %i"%pmid)
	for i in conn.fetchall():
		if i[0]:
			delete_graph(conn,i[0])
	conn.execute("DELETE FROM posmap_info WHERE posmap_set_id = %i"%pmid)
	conn.execute("DELETE FROM posmap_points WHERE posmap_set_id = %i"%pmid)

# delete all posmaps
def delete_all_posmaps(conn):
	conn.execute("SELECT DISTINCT posmap_set_id FROM posmap_info")
	for i in conn.fetchall():
		delete_posmap(conn,i[0])

# delete a run monitor		
def delete_run_monitor(conn,rmid):
	conn.execute("SELECT center_graph_id,width_graph_id FROM run_monitors WHERE monitor_id = %i"%rmid)	
	for i in conn.fetchall():
		delete_graph(conn,i[0])
		delete_graph(conn,i[1])
	conn.execute("DELETE FROM run_monitors WHERE monitor_id = %i"%rmid)
	
def get_sensor_id(conn,sname):
	if sname == "UNKNOWN":
		return 0
	conn.execute("SELECT sensor_id FROM sensors WHERE sensor_name = '%s'"%sname)
	srow = conn.fetchone()
	if not srow:
		return 0
	return srow[0]
	
def get_miscinfo(conn,s):
	conn.execute("SELECT misc_info_id,int_1,int_2,double_1,double_2 FROM misc_info WHERE text_description = '%s'"%s)
	return conn.fetchall();

def set_miscinfo(conn,s,i1=None,i2=None,d1=None,d2=None):
	conn.execute("INSERT INTO misc_info(text_description) VALUES ('%s')"%s);
	conn.execute("SELECT LAST_INSERT_ID()")
	miid = conn.fetchone()[0]
	if i1 is not None:
		conn.execute("UPDATE misc_info SET int_1 = %i WHERE misc_info_id = %i"%(i1,miid))
	if i2 is not None:
		conn.execute("UPDATE misc_info SET int_2 = %i WHERE misc_info_id = %i"%(i2,miid))
	if d1 is not None:
		conn.execute("UPDATE misc_info SET double_1 = %f WHERE misc_info_id = %i"%(d1,miid))
	if d2 is not None:
		conn.execute("UPDATE misc_info SET double_2 = %f WHERE misc_info_id = %i"%(d2,miid))
	return miid		

# get list of rundata files for given run range
def get_rundata_files(rmin=7000,rmax=12000,basepath = "../RunData/"):
	return [ r for r in [ (r,int(r.split('_')[1].split('.')[0])) for r in os.listdir(basepath) if r[:4]=="Run_"] if rmin <= r[1] <= rmax]

# get plots folders for runs in range
def get_plot_folders(rmin,rmax,basepath = "../Plots/"):
	rdirs = []
	for cd in os.listdir(basepath):
		if not os.path.isdir(basepath+cd):
			continue
		rdirs += [r for r in [ (int(d.split()[0].split('_')[0]),basepath+cd+'/'+d) for d in os.listdir(basepath+cd) if d[0].isdigit() ] if rmin <= r[0] <= rmax]
	rdirs.sort()
	return rdirs

					
def delete_LED(rn0,rn1):
	
	# open connection
	conn = open_connection()
	
	# find deletables
	conn.execute("SELECT center_graph_id,width_graph_id,monitor_id,run_number FROM run_monitors WHERE run_number >= %i AND run_number <= %i AND monitor_type = 'GMS_peak'"%(rn0,rn1))
	gids = conn.fetchall()
	for g in gids:
		print g[3]
		conn.executemany("DELETE FROM graph_points WHERE graph_id = %s",[(g[0]),(g[1])])
		conn.executemany("DELETE FROM graphs WHERE graph_id = %s",[(g[0]),(g[1])])
		conn.execute("DELETE FROM run_monitors WHERE monitor_id = %i"%g[2])

# get GMS run for given run
def getGMSrun(conn,rn):
	conn.execute("SELECT gms_run FROM energy_calibration WHERE start_run <= %i AND %i <= end_run ORDER BY end_run-start_run ASC LIMIT 1"%(rn,rn))
	g = conn.fetchone()
	if g:
		return g[0]
	return 0

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

# get list of all GMS runs
def getGMSruns(conn):
	conn.execute("SELECT gms_run FROM energy_calibration WHERE 1 ORDER BY gms_run ASC")
	return [ r[0] for r in conn.fetchall() ]

# get run start times for GMS runs
def gmsRunTimes(conn):
	conn.execute("SELECT run_number,UNIX_TIMESTAMP(start_time) FROM run WHERE run_number IN (SELECT gms_run FROM energy_calibration WHERE 1 ORDER BY gms_run ASC)")
	return conn.fetchall()

					
# delete all Co60 data
def clearCo60():
	conn = open_connection()
	for side in ['E','W']:
		sid = get_sensor_id(conn,"ADCRef%sCo60"%side)
		for pk in range(2):
			conn.execute("SELECT center_graph_id,width_graph_id FROM continuous_monitors WHERE monitor_type = 'Co60_Peak_%i' AND sensors_sensor_id = %i"%(pk+1,sid))
			for gid in conn.fetchone():
				delete_graph(conn,gid)			

# get run monitor graph for sensor name	
def getPeakMonitor(rNum,mName,mtype="GMS_peak",conn=None):
	if not conn:
		conn = open_connection()
	print "Fetching",mtype,mName,rNum
	sid = mName
	if type(mName) == type(1):
		sid = get_sensor_id(conn,mName)
		if not sid:
			print "*** Unknown sensor:",mName
			return []
	conn.execute("SELECT center_graph_id FROM run_monitors WHERE run_number = %i AND sensors_sensor_id = %i AND monitor_type = '%s'"%(rNum,sid,mtype))
	return getGraph(conn,conn.fetchone()[0])

