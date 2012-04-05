#include "SourceDBSQL.hh"
#include "CalDBSQL.hh"

SourceDBSQL* SourceDBSQL::getSourceDBSQL() {
	static SourceDBSQL* SDB = NULL;
	if(!SDB) SDB = new SourceDBSQL(getEnvSafe("UCNADB"),getEnvSafe("UCNADBADDRESS"),
								   getEnvSafe("UCNADBUSER"),getEnvSafe("UCNADBPASS"),atoi(getEnvSafe("UCNADBPORT","3306").c_str()));
	return SDB;
}

std::vector<Source> SourceDBSQL::runSources(RunNum rn, Side s) {	
	//                    0         1          2    3     4     5       6       7      8
	if(s==EAST || s==WEST)
		sprintf(query,"SELECT source_id,run_number,side,x_pos,y_pos,x_width,y_width,counts,sourcetype FROM sources WHERE run_number = %i and side = '%s' ORDER BY x_pos",
				rn,s==EAST?"East":"West");
	else 
		sprintf(query,"SELECT source_id,run_number,side,x_pos,y_pos,x_width,y_width,counts,sourcetype FROM sources WHERE run_number = %i ORDER BY x_pos",rn);
	Query();
	TSQLRow* r;
	std::vector<Source> srcs;
	while((r = res->Next())) {
		srcs.push_back(Source(fieldAsString(r,8), (fieldAsString(r,2)=="East")?EAST:WEST, fieldAsInt(r,0)));
		srcs.back().myRun = fieldAsInt(r,1);
		srcs.back().x = fieldAsFloat(r,3);
		srcs.back().y = fieldAsFloat(r,4);
		srcs.back().wx = fieldAsFloat(r,5);
		srcs.back().wy = fieldAsFloat(r,6);
		srcs.back().nCounts = fieldAsFloat(r,7);
		delete(r);
	}
	return srcs;
}

void SourceDBSQL::clearSources(RunNum rn) {
	sprintf(query,"DELETE sourcepeaks FROM sourcepeaks,sources WHERE sources.source_id = sourcepeaks.source_id and sources.run_number = %i",rn);
	sprintf(query,"DELETE FROM sources WHERE run_number = %i",rn);
	execute();
}

void SourceDBSQL::addSource(const Source& src) {
	if(src.sID) {
		sprintf(query,"UPDATE sources SET x_pos=%f, y_pos=%f, x_width=%f, y_width=%f, counts=%f, sourcetype='%s' WHERE source_id = %i",
				src.x,src.y,src.wx,src.wy,src.nCounts,src.t.c_str(),src.sID);
	} else {
		sprintf(query,"INSERT INTO sources(run_number,side,x_pos,y_pos,x_width,y_width,counts,sourcetype) VALUES (%i,'%s',%f,%f,%f,%f,%f,'%s')",
				src.myRun,src.mySide==EAST?"East":"West",src.x,src.y,src.wx,src.wy,src.nCounts,src.t.c_str());
	}
	execute();
}

void SourceDBSQL::clearPeaks(unsigned int sID) {
	sprintf(query,"DELETE FROM sourcepeaks WHERE source_id = %i",sID);
	execute();
}

void SourceDBSQL::addPeak(const SpectrumPeak& pk) {
	sprintf(query,"INSERT INTO sourcepeaks%s",sm2insert(pk.toStringmap()).c_str());
	execute();
}
