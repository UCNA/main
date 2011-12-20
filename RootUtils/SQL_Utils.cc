#include "SQL_Utils.hh"
#include <stdlib.h>
#include "strutils.hh"

SQLHelper::SQLHelper(const std::string& dbnm,
					 const std::string& dbAddress,
					 const std::string& dbUser,
					 const std::string& dbPass,
					 unsigned int port,
					 unsigned int ntries): db(NULL), res(NULL), dbName(dbnm) {
	
	std::string dbAddressFull = std::string("mysql://")+dbAddress+":"+itos(port)+"/"+dbnm;
	
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
	printf("\n\n************* FAILED TO CONNECT TO ANALYSIS DB! *******************\n\n\n");
	assert(false || IGNORE_DEAD_DB);
}

void SQLHelper::execute() {
	assert(db || IGNORE_DEAD_DB);
	if(res) delete(res);
	res = NULL;
	db->Exec(query);
}

void SQLHelper::Query() { 
	assert(db || IGNORE_DEAD_DB);
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
