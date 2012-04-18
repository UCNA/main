#ifndef SQL_UTILS_HH
#define SQL_UTILS_HH 1


#include <TSQLServer.h>
#include <TSQLRow.h>
#include <TSQLResult.h>
#include <cassert>
#include <vector>
#include <string>
#include <map>
#include "QFile.hh"

#define IGNORE_DEAD_DB false

/// convenient wrapper class for DB access
class SQLHelper {
public:
	
	/// constructor
	SQLHelper(const std::string& dbName,
			  const std::string& dbAddress,
			  const std::string& dbUser,
			  const std::string& dbPass,
			  unsigned int port = 3306,
			  unsigned int ntries = 3);
	
	/// destructor
	virtual ~SQLHelper() { if(res) delete(res); if(db) db->Close(); }
	
	/// get name of DB in use
	std::string getDBName() const { return dbName; }
	
	bool isNullResult;	//< whether field query returned NULL
	
	char query[9182];			//< buffer space for SQL query strings
	/// execute a non-info-returning query
	void execute();
	
protected:
	/// use current query string, return first row
	TSQLRow* getFirst();
	/// get field as string (with default for NULL)
	std::string fieldAsString(TSQLRow* row, unsigned int fieldnum=0, const std::string& dflt = "");	
	/// get field as integer (with default for NULL)
	int fieldAsInt(TSQLRow* row, unsigned int fieldnum=0, int dflt = 0);
	/// get field as float (with default for NULL)
	float fieldAsFloat(TSQLRow* row, unsigned int fieldnum=0, float dflt = 0);
	/// display a query result to stdout
	void printResult();
	/// get most recently inserted row ID
	int getInsertID();
	
	/// check if table exists in DB
	bool checkTable(const std::string& tname) { return db && db->HasTable(tname.c_str()); }
	
	/// execute an info-returning query
	void Query();

protected:	
	TSQLServer* db;				//< DB server connection
	TSQLResult* res;			//< result of most recent query
	std::string dbName;			//< name of DB in use
};

/// convert a stringmap to "(vars,...) VALUES (vals,...)" for DB insert query
std::string sm2insert(const Stringmap& m);

#endif
