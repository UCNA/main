#include "SQL_Utils.hh"
#include <stdlib.h>
#include "strutils.hh"
#include "SMExcept.hh"

SQLHelper::SQLHelper(const std::string& dbnm,
					 const std::string& dbAddress,
					 const std::string& dbUser,
					 const std::string& dbPass,
					 unsigned int port,
					 unsigned int ntries): db(NULL), res(NULL), dbName(dbnm) {
	
	std::string dbAddressFull = "mysql://"+dbAddress+":"+itos(port)+"/"+dbnm;
	
	while(!db) {
		ntries--;
		db = TSQLServer::Connect(dbAddressFull.c_str(),dbUser.c_str(),dbPass.c_str());
		if(!db) {
			if(!ntries || IGNORE_DEAD_DB)
				break;
			printf("** DB Connection %s@%s failed... retrying...\n",dbUser.c_str(),dbAddressFull.c_str());
			sleep(2);
		} else {
			printf("Connected to DB server: %s\n", db->ServerInfo());
			return;
		}
	}
	SMExcept e("DBConnectFail");
	e.insert("dbAddress",dbAddressFull);
	e.insert("dbUser",dbUser);
	throw(e);
}

void SQLHelper::execute() {
	if(res) delete(res);
	res = NULL;
	db->Exec(query);
}

void SQLHelper::Query() { 
	if(!db) {
		res = NULL;
	} else {
		if(res) delete(res);
		res = db->Query(query);
	}
}

std::string SQLHelper::fieldAsString(TSQLRow* row, unsigned int fieldnum, const std::string& dflt) {
	const char* s = row->GetField(fieldnum);
	isNullResult = !s;
	if(isNullResult)
		return dflt;
	return std::string(s);
}

int SQLHelper::fieldAsInt(TSQLRow* row, unsigned int fieldnum, int dflt) {
	std::string s = fieldAsString(row,fieldnum);
	if(isNullResult)
		return dflt;
	return atoi(s.c_str());
}

float SQLHelper::fieldAsFloat(TSQLRow* row, unsigned int fieldnum, float dflt) {
	std::string s = fieldAsString(row,fieldnum);
	if(isNullResult)
		return dflt;
	return atof(s.c_str());
}

int SQLHelper::getInsertID() {
	sprintf(query,"SELECT LAST_INSERT_ID()");
	Query();
	TSQLRow* r = getFirst();
	if(!r)
		throw(SMExcept("failedInsert"));
	int rid = fieldAsInt(r,0);
	delete(r);
	if(!rid)
		throw(SMExcept("failedInsert"));
	return rid;	
}

void SQLHelper::printResult() {
	TSQLRow* row;
	while( (row = res->Next()) ) {
		printf("----------------\n");
		for(int i=0; i<res->GetFieldCount(); i++)
			printf("%s:\t%s\n",res->GetFieldName(i),fieldAsString(row,i,"NULL").c_str());
		delete(row);
	}		
}

TSQLRow* SQLHelper::getFirst() {
	Query();
	if(!res)
		return NULL;
	return res->Next();
}

std::string sm2insert(const Stringmap& m) {
	std::string svars = "(";
	std::string svals = "VALUES (";
	for(std::map<std::string,std::string>::const_iterator it = m.dat.begin(); it != m.dat.end(); it++) {
		if(it != m.dat.begin()) {
			svars += ",";
			svals += ",";
		}
		svars += it->first;
		svals += it->second;
	}
	return svars+") "+svals+")";
}
