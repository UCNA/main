#include "AnalysisDB.hh"

AnaResult::AnaResult(const std::string& auth): arid(0), author(auth), timestamp(time(NULL)),
startRun(0), endRun(0), anach(ANCHOICE_C), s(BOTH), afp(AFP_OTHER), gv(GV_OPEN), value(0), err(0), csid(0) { }

std::string AnaResult::atypeWord(AnaType t) { return t==ANA_ASYM?"Asymmetry":"Counts"; }
std::string AnaResult::dsourceWord(DataSource d) { return d==REAL_DATA?"Data":d==G4_DATA?"G4":"Pen"; }
std::string AnaResult::groupWord(RunGrouping g) { return g==GROUP_RUN?"run":g==GROUP_PAIR?"pair":g==GROUP_QUARTET?"quartet":"octet"; }
AnaResult::AnaType AnaResult::strToAtype(const std::string& str) { return str=="Asymmetry"?ANA_ASYM:ANA_COUNTS; }
AnaResult::DataSource AnaResult::strToDsource(const std::string& str) { return str=="Data"?REAL_DATA:str=="G4"?G4_DATA:PEN_DATA; }

Stringmap AnaResult::toStringmap() const {
	Stringmap m;
	m.insert("author",author);
	m.insert("timestamp",itos(timestamp));
	m.insert("startRun",startRun);
	m.insert("endRun",endRun);
	m.insert("grouping",groupWord(grouping));
	m.insert("evtps",typeSetString());
	m.insert("side",sideWords(s));
	m.insert("afp",afpWords(afp));
	m.insert("value",value);
	m.insert("err",err);
	m.insert("csid",csid);
	m.insert("type",anatp==AnaResult::ANA_ASYM?"Asymmetry":"Counts");
	m.insert("source",datp==AnaResult::REAL_DATA?"Data":datp==AnaResult::G4_DATA?"G4":"Pen");
	return m;
}

std::string AnaResult::typeSetString() const {
	std::string ts;
	for(std::set<EventType>::const_iterator it = etypes.begin(); it != etypes.end(); it++) {
		if(it != etypes.begin()) ts += ",";
		if((*it)<=TYPE_III_EVENT)
			ts += (*it)==TYPE_0_EVENT?"0":(*it)==TYPE_I_EVENT?"I":(*it)==TYPE_II_EVENT?"II":"III";
	}
	return ts;
}


AnalysisDB* AnalysisDB::getADB() {
	static AnalysisDB* ADB = NULL;
	if(!ADB) ADB = new AnalysisDB();
	return ADB;
}

unsigned int AnalysisDB::uploadCutSpec(AnaCutSpec& c) {
	sprintf(query,"INSERT INTO cut_spec(energy_min,energy_max,radius,positioning) VALUES (%f,%f,%f,'%s')",
			c.emin,c.emax,c.radius,c.postp == AnaCutSpec::POS_PLAIN?"plain":"rotated");
	execute();
	c.csid = getInsertID();
	return c.csid;
}

unsigned int AnalysisDB::uploadAnaResult(AnaResult& r) {
	printf("Uploading analysis result:\n");
	r.toStringmap().display();
	r.value = r.value==r.value && fabs(r.value) < FLT_MAX?r.value:0;
	r.err = r.err==r.err && fabs(r.err) < FLT_MAX?r.err:0;
	sprintf(query,"INSERT INTO analysis_results(author,date,type,source,start_run,end_run,grouping,event_type,ana_choice,side,afp,gate_valve,value,err,cut_spec_id) \
			VALUES ('%s',FROM_UNIXTIME(%u),'%s','%s',%i,%i,'%s','%s','%c',%s,'%s','%s',%f,%f,%u)",
			r.author.c_str(),
			(unsigned int)r.timestamp,
			AnaResult::atypeWord(r.anatp).c_str(),
			AnaResult::dsourceWord(r.datp).c_str(),
			r.startRun,
			r.endRun,
			AnaResult::groupWord(r.grouping).c_str(),
			r.typeSetString().c_str(),
			choiceLetter(r.anach),
			dbSideName(r.s),
			afpWords(r.afp).c_str(),
			gvWords(r.gv).c_str(),
			r.value,
			r.err,
			r.csid);
	execute();
	r.arid = getInsertID();
	return r.arid;
}

void AnalysisDB::deleteAnaResult(unsigned int arid) {
	printf("Deleting analysis result %i\n",arid);
	sprintf(query,"SELECT cut_spec_id FROM analysis_results WHERE analysis_results_id = %i",arid);
	Query();
	TSQLRow* r = getFirst();
	if(!r) return;
	int csid = fieldAsInt(r); 
	delete(r);
	
	sprintf(query,"DELETE FROM analysis_results WHERE analysis_results_id = %i",arid);
	execute();
	
	sprintf(query,"SELECT COUNT(*) FROM analysis_results WHERE cut_spec_id = %i",csid);
	Query();
	r = getFirst();
	assert(r);
	if(!fieldAsInt(r))
		deleteCutSpec(csid);
	delete(r);
}

void AnalysisDB::deleteCutSpec(unsigned int csid) {
	printf("Deleting CutSpec %i\n",csid);
	sprintf(query,"DELETE FROM cut_spec WHERE cut_spec_id = %i",csid);
	execute();
}

AnaResult AnalysisDB::getAnaResult(unsigned int arid) {
	//                    0      1    2    3      4         5       6          7          8    9   10    11  12
	sprintf(query,"SELECT author,date,type,source,start_run,end_run,event_type,ana_choice,side,afp,gate_valve,value,err,cut_spec_id \
			FROM analysis_results WHERE analysis_results_id = %i",arid);
	TSQLRow* r = getFirst();
	if(!r) {
		SMExcept e("MissingAnaResult");
		e.insert("analysis_results_id",arid);
		throw(e);
	}
	AnaResult a(fieldAsString(r,0));
	a.arid = arid;
	a.timestamp = fieldAsInt(r,1);
	a.anatp = AnaResult::strToAtype(fieldAsString(r,2));
	a.datp = AnaResult::strToDsource(fieldAsString(r,3));
	a.startRun = fieldAsInt(r,4);
	a.endRun = fieldAsInt(r,5);
	std::vector<std::string> tps = split(fieldAsString(r,6));
	for(std::vector<std::string>::iterator it = tps.begin(); it != tps.end(); it++)
		a.etypes.insert((*it)=="0"?TYPE_0_EVENT:EventType(TYPE_0_EVENT+it->size()));
	a.anach = AnalysisChoice(fieldAsString(r,7)[0]-'A'+1);
	a.s = strToSide(fieldAsString(r,8));
	a.afp = strToAfp(fieldAsString(r,9));
	a.gv = strToGV(fieldAsString(r,10));
	a.value = fieldAsFloat(r,11);
	a.err = fieldAsFloat(r,12);
	a.csid = fieldAsInt(r,13);
	delete(r);
	return a;
}

AnaCutSpec AnalysisDB::getCutSpec(unsigned int csid) {
	AnaCutSpec c;
	sprintf(query,"SELECT energy_min,energy_max,radius,positioning FROM cut_spec WHERE cut_spec_id = %i",csid);
	TSQLRow* r = getFirst();
	if(!r) {
		SMExcept e("MissingCutSpec");
		e.insert("cut_spec_id",csid);
		throw(e);
	}
	c.emin = fieldAsFloat(r,0);
	c.emax = fieldAsFloat(r,1);
	c.radius = fieldAsFloat(r,2);
	c.postp = fieldAsString(r,3)=="plain"?AnaCutSpec::POS_PLAIN:AnaCutSpec::POS_ROTATED;
	delete(r);
	return c;
}

std::vector<AnaResult> AnalysisDB::findMatching(const AnaResult& A) {
	std::string qry = "SELECT analysis_results_id FROM analysis_results WHERE author = '"+A.author+"'";
	qry += " AND type = '"+AnaResult::atypeWord(A.anatp)+"'";
	qry += " AND source = '"+AnaResult::dsourceWord(A.datp)+"'";
	qry += " AND ana_choice = '"+ctos(choiceLetter(A.anach))+"'";
	qry += " AND side = "; qry += dbSideName(A.s);
	qry += " AND afp = '"; qry += afpWords(A.afp)+"'";
	qry += " AND gate_valve = '"; qry += gvWords(A.gv)+"'";
	qry += " AND event_type = '"+A.typeSetString()+"'";
	if(A.startRun)
		qry += " AND start_run = "+itos(A.startRun);
	if(A.endRun)
		qry += " AND end_run = "+itos(A.endRun);
	sprintf(query,"%s",qry.c_str());
	Query();
	std::vector<unsigned int> arids;
	while(TSQLRow* r = res->Next()) {
		arids.push_back(fieldAsInt(r,0));
		delete(r);
	}
	std::vector<AnaResult> v;
	for(unsigned int i=0; i<arids.size(); i++)
		v.push_back(getAnaResult(arids[i]));
	return v;
}
