#include "CalDBSQL.hh"
#include "PathUtils.hh"
#include <utility>

CalDBSQL* CalDBSQL::getCDB(bool readonly) {
	if(readonly) {
		static CalDBSQL* CDBr = NULL;
		if(!CDBr) CDBr = new CalDBSQL(getEnvSafe("UCNADB"),getEnvSafe("UCNADBADDRESS"),
									  getEnvSafe("UCNADBUSER_READONLY"),getEnvSafe("UCNADBPASS_READONLY"));
		return CDBr;
	} else {
		static CalDBSQL* CDBw = NULL;
		if(!CDBw) CDBw = new CalDBSQL(getEnvSafe("UCNADB"),getEnvSafe("UCNADBADDRESS"),
									  getEnvSafe("UCNADBUSER"),getEnvSafe("UCNADBPASS"));
		return CDBw;		
	}
}

CalDBSQL::~CalDBSQL() {
	for(unsigned int i=0; i<pcors.size(); i++)
		delete pcors[i];
}

TGraphErrors* CalDBSQL::getRunMonitor(RunNum rn, const std::string& sensorName, const std::string& monType, bool centers) {
	if(centers)
		sprintf(query,"SELECT center_graph_id FROM run_monitors,sensors WHERE run_number = %i \
				AND sensors_sensor_id = sensor_id AND sensor_name = '%s' AND monitor_type = '%s'",rn,sensorName.c_str(),monType.c_str());
	else
		sprintf(query,"SELECT width_graph_id FROM run_monitors,sensors WHERE run_number = %i \
				AND sensors_sensor_id = sensor_id AND sensor_name = '%s' AND monitor_type = '%s'",rn,sensorName.c_str(),monType.c_str());
	TSQLRow* r = getFirst();
	if(!r) {
		printf("Failed to locate requested monitor:\n\t%s\n",query);
		return NULL;
	}
	TGraphErrors* tg = getGraph(fieldAsInt(r,0));
	delete(r);
	if(!tg)
		printf("*** Warning: no graph for run monitor %i %s (%s)\n",rn,sensorName.c_str(),monType.c_str());
	return tg;
}

float CalDBSQL::getRunMonitorStart(RunNum rn, const std::string& sensorName, const std::string& monType) {
	sprintf(query,"SELECT center_graph_id FROM run_monitors,sensors WHERE run_number = %i \
			AND sensors_sensor_id = sensor_id AND sensor_name = '%s' AND monitor_type = '%s'",rn,sensorName.c_str(),monType.c_str());
	TSQLRow* r = getFirst();
	if(!r) {
		printf("Failed to locate requested monitor start:\n\t%s\n",query);
		return 0;
	}
	sprintf(query,"SELECT y_value FROM graph_points WHERE graph_id = %i ORDER BY x_value ASC LIMIT 1",fieldAsInt(r,0));
	delete(r);
	r = getFirst();
	if(!r) {
		printf("Failed to locate requested monitor start:\n\t%s\n",query);
		return 0;
	}
	float z = fieldAsFloat(r,0);
	delete(r);
	return z;
}


TGraphErrors* CalDBSQL::getContinuousMonitor(const std::string& sensorName, const std::string& monType, RunNum rn, bool centers) {
	sprintf(query,"SELECT center_graph_id,width_graph_id FROM continuous_monitors,sensors \
			WHERE sensors_sensor_id = sensor_id AND sensor_name = '%s' AND monitor_type = '%s'",sensorName.c_str(),monType.c_str());
	TSQLRow* r = getFirst();
	if(!r) {
		printf("Failed to locate requested monitor:\n\t%s\n",query);
		return NULL;
	}
	TGraphErrors* tg; 
	if(centers)
		tg = getGraph(fieldAsInt(r,0),rn);
	else
		tg = getGraph(fieldAsInt(r,1),rn);
	delete(r);
	if(tg && !tg->GetN()) {
		delete(tg);
		return NULL;
	}
	return tg;
}


float CalDBSQL::getKurieEnergy(RunNum rn, Side s, unsigned int t) {
	sprintf(query,"SELECT energy FROM kurie_cal WHERE run_number = %i AND side = %s AND quadrant = %i",rn,dbSideName(s),t);
	TSQLRow* r = getFirst();		
	if(!r)
		return 0;
	float ken = fieldAsFloat(r,0);
	delete(r);
	return ken;	
}
float CalDBSQL::getKurieADC(RunNum rn, Side s, unsigned int t) {
	sprintf(query,"SELECT adc FROM kurie_cal WHERE run_number = %i AND side = %s AND quadrant = %i",rn,dbSideName(s),t);
	TSQLRow* r = getFirst();		
	if(!r)
		return 0;
	float kadc = fieldAsFloat(r,0);
	delete(r);
	return kadc;		
}



PositioningCorrector* CalDBSQL::getPositioningCorrector(RunNum rn) {
	return getPositioningCorrectorByID(getCalSetInfo(rn,"posmap_set_id"));
}


float CalDBSQL::getAnodeCalInfo(RunNum R, const char* field) {
	sprintf(query,"SELECT %s FROM anode_cal WHERE start_run <= %i AND %i <= end_run ORDER BY end_run-start_run LIMIT 1",field,R,R);
	TSQLRow* r = getFirst();
	if(!r) {
		sprintf(query,"SELECT %s,start_run,end_run FROM anode_cal ORDER BY pow(1.0*end_run-%i,2)+pow(1.0*start_run-%i,2) LIMIT 1",field,R,R);
		r = getFirst();
		assert(r);
		unsigned int r0 = fieldAsInt(r,1);
		unsigned int r1 = fieldAsInt(r,2);
		printf("**** WARNING: No matching anode calibration found for %i; using nearest range (%i,%i)...\n",R,r0,r1);
	}
	float f = fieldAsFloat(r,0);
	delete(r);
	return f;
}

PositioningCorrector* CalDBSQL::getAnodePositioningCorrector(RunNum rn) {	
	return getPositioningCorrectorByID(int(getAnodeCalInfo(rn,"anode_posmap_id")));
}

float CalDBSQL::getAnodeGain(RunNum rn, Side s) {
	return getAnodeCalInfo(rn,sideSubst("calfactor_%c",s).c_str());
}


PositioningCorrector* CalDBSQL::getPositioningCorrectorByID(unsigned int psid) {
	
	std::map<unsigned int,PositioningCorrector*>::iterator it = pcors.find(psid);
	if(it != pcors.end())
		return it->second;
	
	printf("Loading positioning corrector %i...\n",psid);
	std::vector<PosmapInfo> pinf;
	TSQLRow* r;
	for(Side s = EAST; s<=WEST; s=nextSide(s)) {
		for(unsigned int t=0; t<=nBetaTubes; t++) {
			sprintf(query,"SELECT n_rings,radius FROM posmap_info WHERE posmap_set_id = %i AND side = %s AND quadrant = %i",psid,dbSideName(s),t);
			r = getFirst();
			if(!r)
				continue;
			pinf.push_back(PosmapInfo());
			pinf.back().s = s;
			pinf.back().t = t;
			pinf.back().nRings = fieldAsInt(r,0);
			pinf.back().radius = fieldAsFloat(r,1);
			delete(r);
			sprintf(query,
					"SELECT endpoint_adc,smear_correction FROM posmap_points WHERE posmap_set_id = %i AND side = %s AND quadrant = %i ORDER BY pixel_id ASC",
					psid,dbSideName(s),t);
			Query();
			if(!res)
				return NULL;
			while((r = res->Next())) {
				pinf.back().adc.push_back(fieldAsFloat(r,0));
				pinf.back().energy.push_back(fieldAsFloat(r,1));
				delete(r);
			}
			// delete entries with no data points
			if(!pinf.back().adc.size())
				pinf.pop_back();
		}
	}
	assert(pinf.size() || IGNORE_DEAD_DB);
	if(!pinf.size()) return NULL;
	
	pcors.insert(std::make_pair(psid,new PositioningCorrector(pinf)));
	return getPositioningCorrectorByID(psid);
}

unsigned int CalDBSQL::getCalSetInfo(RunNum R, const char* field) {
	sprintf(query,"SELECT %s,start_run FROM energy_calibration WHERE start_run <= %i AND %i <= end_run ORDER BY end_run - start_run LIMIT 1",field,R,R);
	TSQLRow* r = getFirst();
	if(!r)
		return 0;
	int i = fieldAsInt(r,0);
	if(fieldAsInt(r,1) == 1)
		printf("**** WARNING: Catchall Calibration selected for Run %i\n",R);
	delete(r);
	return i;
}

unsigned int CalDBSQL::getTubecalID(RunNum R, Side s, unsigned int t) {
	unsigned int csid = getCalSetInfo(R,"ecal_id");
	sprintf(query,"SELECT tubecal_id FROM tube_calibration WHERE ecal_id = %i AND side = %s AND quadrant = %i",csid,dbSideName(s),t);
	TSQLRow* r = getFirst();
	if(!r) {
		printf("Failed to locate Tube Data:\n\t%s\n;",query);
		return 0;
	}
	int tsid = fieldAsInt(r); 
	delete(r);
	return tsid;		
}

float CalDBSQL::getTubecalData(RunNum rn, Side s, unsigned int t, const char* field) {
	unsigned int tsid = getTubecalID(rn, s, t);
	assert(tsid || IGNORE_DEAD_DB);
	sprintf(query,"SELECT %s FROM tube_calibration WHERE tubecal_id = %i",field,tsid);
	TSQLRow* r = getFirst();
	if(!r)
		return 0;
	float x = fieldAsFloat(r); 
	delete(r);
	return x;
}

int CalDBSQL::getTubecalInt(RunNum rn, Side s, unsigned int t, const char* field) {
	unsigned int tsid = getTubecalID(rn, s, t);
	assert(tsid || IGNORE_DEAD_DB);
	sprintf(query,"SELECT %s FROM tube_calibration WHERE tubecal_id = %i",field,tsid);
	TSQLRow* r = getFirst();
	if(!r)
		return 0;
	int x = fieldAsInt(r); 
	delete(r);
	return x;
}


const char* CalDBSQL::dbSideName(Side s) const {
	if(s==EAST)
		return "'East'";
	if(s==WEST)
		return "'West'";
	return "NULL";
}

int CalDBSQL::startTime(RunNum rn, int t0) {
	sprintf(query,"SELECT UNIX_TIMESTAMP(start_time)-%i FROM run WHERE run_number = %u",t0,rn);
	TSQLRow* row = getFirst();
	if(!row)
		return -t0;
	t0 = fieldAsInt(row,0,-t0);
	delete(row);
	return t0;
}

int CalDBSQL::endTime(RunNum rn, int t0) {
	sprintf(query,"SELECT UNIX_TIMESTAMP(end_time)-%i FROM run WHERE run_number = %u",t0,rn);
	TSQLRow* row = getFirst();
	if(!row)
		return -t0;
	t0 = fieldAsInt(row,0,-t0);
	delete(row);
	return t0;	
}

BlindTime CalDBSQL::fiducialTime(RunNum rn) {
	BlindTime b = 0;
	sprintf(query,"SELECT live_time_e,live_time_w,live_time FROM analysis WHERE run_number = %u",rn);
	TSQLRow* row = getFirst();
	if(!row)
		return 0;
	b.t[EAST] = fieldAsFloat(row,0);
	b.t[WEST] = fieldAsFloat(row,1);
	b.t[BOTH] = fieldAsFloat(row,2);
	delete(row);
	return b;		
}

float CalDBSQL::totalTime(RunNum rn) {
	sprintf(query,"SELECT total_time FROM analysis WHERE run_number = %u",rn);
	TSQLRow* row = getFirst();
	if(!row)
		return 0;
	float t = fieldAsFloat(row,0);
	delete(row);
	return t;	
}

std::string CalDBSQL::getGroupName(RunNum rn) {
	sprintf(query,"SELECT name,start_run FROM run_group WHERE start_run <= %i AND %i <= end_run ORDER BY end_run-start_run LIMIT 1",rn,rn);
	TSQLRow* r = getFirst();
	if(!r)
		return "0 Unknown";
	std::string nm = fieldAsString(r,0);
	int r0 = fieldAsInt(r,1);
	delete(r);
	return itos(r0)+" "+nm;
}

RunInfo CalDBSQL::getRunInfo(RunNum r) {
	
	RunInfo R = RunInfo(r);
	R.startTime = startTime(r);
	R.groupName = getGroupName(r);
	
	//                    0               1        2        3          4       5
	sprintf(query,"SELECT slow_run_number,run_type,asym_oct,gate_valve,flipper,scs_field FROM run WHERE run_number = %u",r);
	TSQLRow* row = getFirst();
	assert(row);
	
	R.roleName = fieldAsString(row,2,"Other");
	R.octet = 0;
	if(R.roleName[0] == 'A' || R.roleName[0] == 'B') {
		R.octet += atoi(R.roleName.c_str()+1);
		if(R.roleName[0] == 'B')
			R.octet += 12;
	}
	if(R.octet)
		R.triad = TriadType((R.octet-1)/3+1);
	
	R.slowDaq = fieldAsInt(row,0,0);
	std::string s = fieldAsString(row,1,"Other");
	if(s=="Asymmetry")
		R.type = ASYMMETRY_RUN;
	else if(s=="LEDCalib")
		R.type = LED_RUN;
	else if(s=="SourceCalib") {
		R.type = SOURCE_RUN;
		R.roleName = "SourcesCal";
	} else if(s=="Xe") {
		R.type = XE_RUN;
		R.roleName = "Xe";
	} else
		R.type = UNKNOWN_RUN;
	
	s = fieldAsString(row,3,"Other");
	if(s=="Open")
		R.gvState = GV_OPEN;
	else if(s=="Closed")
		R.gvState = GV_CLOSED;
	else
		R.gvState = GV_OTHER;
	
	s=fieldAsString(row,4,"Other");
	if(s=="On")
		R.afpState = AFP_ON;
	else if (s=="Off")
		R.afpState = AFP_OFF;
	else if (s=="On2Off")
		R.afpState = AFP_ON2OFF;
	else if (s=="Off2On")
		R.afpState = AFP_OFF2ON;
	else
		R.afpState = AFP_OTHER;
	
	R.scsField = fieldAsFloat(row,5,0);
	
	delete(row);
	
	return R;
}


std::vector<RunNum> CalDBSQL::findRuns(const char* whereConditions) {
	sprintf(query,"SELECT run_number FROM run WHERE %s",whereConditions);
	Query();
	std::vector<RunNum> v;
	if(!res)
		return v;
	TSQLRow* row;
	while((row = res->Next())) {
		v.push_back((RunNum)fieldAsInt(row,0,0));
		delete(row);
	}
	return v;
}

EfficCurve* CalDBSQL::getTrigeff(RunNum rn, Side s, unsigned int t) {
	sprintf(query,"SELECT params_graph FROM mpm_trigeff WHERE run_number = %i AND side = %s AND quadrant = %i",rn,dbSideName(s),t);
	Query();
	if(!res)
		return NULL;
	TSQLRow* r = res->Next();
	if(!r)
		return NULL;
	int gid = fieldAsInt(r,0);
	delete(r);
	sprintf(query,"SELECT y_value FROM graph_points WHERE graph_id = %i ORDER BY x_value ASC",gid);
	Query();
	if(!res)
		return NULL;
	std::vector<double> v;
	while((r = res->Next())) {
		v.push_back(fieldAsFloat(r,0));
		delete(r);
	}
	if(v.size() != 4)
		return NULL;
	EfficCurve* C = new EfficCurve();
	for(unsigned int i=0; i<4; i++)
		C->params[i] = v[i];
	return C;
}

TGraph* CalDBSQL::getEvisConversion(RunNum rn, Side s, EventType tp) {
	sprintf(query,"SELECT conversion_curve_id FROM evis_conversion WHERE side = %s AND type = %i \
			AND start_run <= %i AND %i <= end_run ORDER BY end_run - start_run LIMIT 1",dbSideName(s),tp,rn,rn);
	TSQLRow* r = getFirst();
	if(!r)
		return NULL;
	int gid = fieldAsInt(r,0);
	delete(r);
	return getGraph(gid);
}

TGraphErrors* CalDBSQL::getGraph(unsigned int gid) {
	std::vector<float> gdata[4];
	sprintf(query,"SELECT x_value,x_error,y_value,y_error FROM graph_points WHERE graph_id = %i ORDER BY x_value ASC",gid);
	Query();
	if(!res) {
		printf("*** Warning: no graph found for <%i>!\n",gid);
		return NULL;
	}
	TSQLRow* r;
	while((r = res->Next())) {
		for(unsigned int i=0; i<4; i++)
			gdata[i].push_back(fieldAsFloat(r,i));
		delete(r);
	}
	unsigned int npts = gdata[0].size();
	if(!npts) {
		printf("*** Warning: no points found for graph <%i>!\n",gid);
		return NULL;
	}
	if(npts == 1) {
		printf("Notice: only 1 graph point found for <%i>; extending to 2.\n",gid);
		for(unsigned int i=0; i<4; i++)
			gdata[i].push_back(gdata[i].back());
		gdata[0].back() += 10.0;
		npts++;
	}		
	TGraphErrors* tg = new TGraphErrors(npts);
	for(unsigned int i=0; i<npts; i++) {
		tg->SetPoint(i,gdata[0][i],gdata[2][i]);
		tg->SetPointError(i,gdata[1][i],gdata[3][i]);
	}
	return tg;
}

TGraphErrors* CalDBSQL::getGraph(unsigned int gid, RunNum rn) {
	
	// determine start time
	sprintf(query,"SELECT x_value FROM graph_points WHERE graph_id = %i \
			AND x_value < %i ORDER BY %i-x_value ASC LIMIT 1",gid,startTime(rn),startTime(rn));
	Query();
	if(!res) {
		printf("*** Warning: no graph found for <%i>!\n",gid);
		return NULL;
	}
	TSQLRow* r = res->Next();
	if(!r) {
		sprintf(query,"SELECT x_value FROM graph_points WHERE graph_id = %i ORDER BY x_value ASC LIMIT 1",gid);
		Query();
		r = res->Next();
	}
	if(!r) {
		printf("*** Warning: no points found for graph <%i>!\n",gid);
		return NULL;
	}
	float tstart = fieldAsFloat(r);
	delete(r);
	
	// determine end time
	sprintf(query,"SELECT x_value FROM graph_points WHERE graph_id = %i \
			AND x_value > %i ORDER BY x_value-%i ASC LIMIT 1",gid,endTime(rn),endTime(rn));
	Query();
	r = res->Next();
	if(!r) {
		sprintf(query,"SELECT x_value FROM graph_points WHERE graph_id = %i ORDER BY x_value DESC LIMIT 1",gid);
		Query();
		r = res->Next();
	}
	if(!r)
		return NULL;
	float tend = fieldAsFloat(r);
	delete(r);
	
	// gather data
	sprintf(query,"SELECT x_value-%i,x_error,y_value,y_error FROM graph_points WHERE graph_id = %i \
			AND x_value >= %f AND x_value <= %f ORDER BY x_value ASC",startTime(rn),gid,tstart-10,tend+10);
	Query();
	std::vector<float> gdata[4];
	while((r = res->Next())) {
		for(unsigned int i=0; i<4; i++)
			gdata[i].push_back(fieldAsFloat(r,i));
		delete(r);
	}
	
	// compile graph
	unsigned int npts = gdata[0].size();
	if(!npts)
		return NULL;
	if(npts == 1) {
		printf("Notice: only 1 graph point found for <%i;%i>; extending to 2.\n",gid,rn);
		for(unsigned int i=0; i<4; i++)
			gdata[i].push_back(gdata[i].back());
		gdata[0].back() += 10.0;
		npts++;
	}
	TGraphErrors* tg = new TGraphErrors(npts);
	for(unsigned int i=0; i<npts; i++) {
		tg->SetPoint(i,gdata[0][i],gdata[2][i]);
		tg->SetPointError(i,gdata[1][i],gdata[3][i]);
	}
	return tg;
}


unsigned int CalDBSQL::newGraph(const std::string& description) {
	sprintf(query,"INSERT INTO graphs (text_description) VALUES ('%s')",description.c_str());
	execute();
	sprintf(query,"SELECT LAST_INSERT_ID()");
	Query();
	TSQLRow* r = getFirst();
	assert(r);
	int gid = fieldAsInt(r,0);
	delete(r);
	return gid;
}

void CalDBSQL::deleteGraph(unsigned int gid) {
	printf("Deleting graph %i...\n",gid);
	sprintf(query,"DELETE FROM graph_points WHERE graph_id = %i",gid);
	execute();
	sprintf(query,"DELETE FROM graphs WHERE graph_id = %i",gid);
	execute();
}

unsigned int CalDBSQL::uploadGraph(const std::string& description, std::vector<double> x, std::vector<double> y,
								   std::vector<double> dx, std::vector<double> dy) {
	assert(x.size()==y.size()||!x.size());
	unsigned int gid = newGraph(description);
	for(unsigned int i=0; i<y.size(); i++) {
		double pdx = dx.size()>i?dx[i]:0;
		double pdy = dy.size()>i?dy[i]:0;
		double xi = x.size()==y.size()?x[i]:i;
		sprintf(query,"INSERT INTO graph_points (graph_id, x_value, y_value, x_error, y_error) VALUES (%i,%f,%f,%f,%f)",gid,xi,y[i],pdx,pdy);
		execute();
	}
	printf("Uploaded graph '%s' to %i (%i points)\n",description.c_str(),gid,(int)y.size());
	return gid;
}

void CalDBSQL::uploadTrigeff(RunNum rn, Side s, unsigned int t, std::vector<double> params, std::vector<double> dparams) {
	unsigned int pgid = uploadGraph("Trigger Efficiency Params",std::vector<double>(),params,std::vector<double>(),dparams);
	sprintf(query,"INSERT INTO mpm_trigeff(run_number,side,quadrant,params_graph) VALUES (%i,'%s',%i,%i)",
			rn,sideWords(s),t,pgid);
	execute();
}

void CalDBSQL::deleteTrigeff(RunNum rn, Side s, unsigned int t) {
	sprintf(query,"SELECT params_graph FROM mpm_trigeff WHERE run_number = %i AND side = '%s' AND quadrant = %i",rn,sideWords(s),t);
	Query();
	TSQLRow* r;
	std::vector<int> gids;
	while((r = res->Next())) {
		gids.push_back(fieldAsInt(r,0));
		delete(r);
	}
	while(gids.size()) {
		deleteGraph(gids.back());
		gids.pop_back();
	}
	sprintf(query,"DELETE FROM mpm_trigeff WHERE run_number = %i AND side = '%s' AND quadrant = %i",rn,sideWords(s),t);
	execute();
}


void CalDBSQL::deleteRunMonitor(RunNum rn, const std::string& sensorName, const std::string& monType) {
	sprintf(query,"SELECT center_graph_id, width_graph_id,monitor_id FROM run_monitors,sensors WHERE run_number = %i \
			AND sensors_sensor_id = sensor_id AND sensor_name = '%s' AND monitor_type = '%s'",rn,sensorName.c_str(),monType.c_str());
	Query();
	TSQLRow* r;
	std::vector<int> rmids;
	std::vector<int> gids;
	while((r = res->Next())) {
		gids.push_back(fieldAsInt(r,0));
		gids.push_back(fieldAsInt(r,1));
		rmids.push_back(fieldAsInt(r,2));
		delete(r);
	}
	while(gids.size()) {
		deleteGraph(gids.back());
		gids.pop_back();
	}
	while(rmids.size()) {
		sprintf(query,"DELETE FROM run_monitors WHERE monitor_id = %i",rmids.back());
		execute();
		rmids.pop_back();
	}
}

unsigned int CalDBSQL::getSensorID(const std::string& sname) {
	sprintf(query,"SELECT sensor_id FROM sensors WHERE sensor_name = '%s'",sname.c_str());
	Query();
	TSQLRow* r = getFirst();
	if(!r) {
		printf("Failed to locate requested sensor:\n\t%s\n",query);
		return 0;
	}
	unsigned int sid = fieldAsInt(r,0);
	delete(r);
	return sid;
}

void CalDBSQL::addRunMonitor(RunNum rn, const std::string& sensorName, const std::string& monType, unsigned int cgid, unsigned int wgid) {
	unsigned int sid = getSensorID(sensorName);
	sprintf(query,"INSERT INTO run_monitors (run_number,sensors_sensor_id,monitor_type,center_graph_id,width_graph_id) \
			VALUES (%i,%i,'%s',%i,%i)",rn,sid,monType.c_str(),cgid,wgid);
	execute();
}

