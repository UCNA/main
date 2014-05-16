#include "AnalysisDB.hh"
#include <time.h>

bool AnalysisDB::disableADB = false;

std::string typeSetString(const std::set<EventType>& etypes)  {
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
	if(disableADB) return NULL;
	if(!ADB) ADB = new AnalysisDB();
	return ADB;
}

void AnalysisDB::locateRunset(AnaRunset& AR) {
	std::string qry =  ( "SELECT analysis_runset_id FROM analysis_runset WHERE start_run = " + itos(AR.start_run)
						+ " AND end_run = " + itos(AR.end_run) + " AND grouping = '" + groupWords(AR.grouping)
						+ "' AND gate_valve = '" + gvWords(AR.gate_valve) + "' AND afp = '" + afpWords(AR.afp) + "'");
	sprintf(query,"%s",qry.c_str());
	printf("%s\n",query);
	TSQLRow* r = getFirst();
	if(r) {
		AR.rsid = fieldAsInt(r,0);
		delete(r);
		return;
	}
	qry = ( "INSERT INTO analysis_runset(start_run,end_run,grouping,gate_valve,afp) VALUES ("
		   + itos(AR.start_run) + "," + itos(AR.end_run) + ",'" + groupWords(AR.grouping) + "','"
		   + gvWords(AR.gate_valve) + "','" + afpWords(AR.afp) +"')");
	sprintf(query,"%s",qry.c_str());
	Query();
	sprintf(query,"SELECT LAST_INSERT_ID()");
	r = getFirst();
	if(!r) {
		SMExcept e("MissingRunset");
		throw(e);
	}
	AR.rsid = fieldAsInt(r,0);
	delete(r);
}

std::vector<AnaNumber> AnalysisDB::findMatching(const AnaNumber& AN) {
	std::string qry = ( "SELECT analysis_number_id,date,value,err FROM analysis_numbers WHERE analysis_runset_id = "
					   + itos(AN.rsid) + " AND source = '" + AN.source + "' AND name = '" + AN.name + "' AND side = " + dbSideName(AN.s)
					   + " AND event_type = '" + typeSetString(AN.etypes) + "' AND n = " + itos(AN.n) + " ORDER BY date DESC");
	sprintf(query,"%s",qry.c_str());
	Query();
	std::vector<AnaNumber> v;
	while(TSQLRow* r = res->Next()) {
		v.push_back(AN);
		v.back().anid = fieldAsInt(r,0);
		v.back().value = fieldAsFloat(r,2);
		v.back().err = fieldAsFloat(r,3);
		delete(r);
	}
	return v;
}

void AnalysisDB::uploadAnaNumber(AnaNumber& AN, bool replace) {
	
	if(!AN.date) AN.date = time(NULL);
	
	if(replace) {
		std::vector<AnaNumber> v = findMatching(AN);
		if(v.size()) {
			AN.anid = v[0].anid;
			std::string qry = ("UPDATE analysis_numbers SET date=FROM_UNIXTIME(" + itos((unsigned int)AN.date)
							   +"), value=" + dtos(AN.value,"NULL") +", err= " + dtos(AN.err,"NULL")
							   + " WHERE analysis_number_id = " + itos(AN.anid));
			sprintf(query,"%s",qry.c_str());
			Query();
			return;
		}
	}
	
	std::string qry = ("INSERT INTO analysis_numbers(analysis_runset_id,source,name,date,side,event_type,n,value,err) VALUES ("
					   + itos(AN.rsid) + ",'" + AN.source + "','" + AN.name + "',FROM_UNIXTIME(" + itos((unsigned int)AN.date)
					   + ")," + dbSideName(AN.s) + ",'" + typeSetString(AN.etypes) + "'," + itos(AN.n)
					   + "," + dtos(AN.value,"NULL") + "," + dtos(AN.err,"NULL") +")");
	sprintf(query,"%s",qry.c_str());
	Query();
	sprintf(query,"SELECT LAST_INSERT_ID()");
	TSQLRow* r = getFirst();
	if(!r) {
		SMExcept e("MissingAnaNumber");
		throw(e);
	}
	AN.anid = fieldAsInt(r,0);
	delete(r);
}
