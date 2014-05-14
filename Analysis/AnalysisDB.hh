#ifndef ANALYSISDB_HH
#define ANALYSISDB_HH

#include "SQL_Utils.hh"
#include "Enums.hh"
#include "Types.hh"
#include "PathUtils.hh"
#include "SMExcept.hh"
#include <string>
#include <set>
#include <time.h>

/// description of analysis run set
class AnaRunset {
public:
	/// constructor
	AnaRunset(): rsid(0), start_run(0), end_run(0), grouping(GROUP_RANGE), gate_valve(GV_OTHER), afp(AFP_OTHER) {}
	
	unsigned int rsid;		///< run set ID number
	RunNum start_run;		///< star run
	RunNum end_run;			///< end run
	RunGrouping grouping;	///< type of run grouping
	GVState gate_valve;		///< gate valve state
	AFPState afp;			///< AFP state
};

/// generic analysis number entry
class AnaNumber {
public:
	/// constructor
	AnaNumber(const std::string& nm = "", const std::string& sr = "Unknown"):
	anid(0), rsid(0), source(sr), name(nm), date(0), s(BOTH), n(0), value(0), err(0) {}
	
	unsigned int anid;			///< analysis number ID
	unsigned int rsid;			///< run set ID
	std::string source;			///< data source/author
	std::string name;			///< entry name
	time_t date;				///< result generation timestamp
	Side s;						///< applicable result side
	std::set<EventType> etypes;	///< event types considered
	int n;						///< result integer enumerator
	double value;				///< result value
	double err;					///< result error
};

/// string for type set
std::string typeSetString(const std::set<EventType>& tps);

/// DB interface for sharing analysis results
class AnalysisDB: public SQLHelper, NoCopy {
public:
	/// globally available AnalysisDB
	static AnalysisDB* getADB();
	/// disable analysis DB access
	static bool disableADB;
	
	/// locate existing or generate new AnaRunset specifier
	void locateRunset(AnaRunset& AR);
	/// locate matching analysis numbers (same run set, author, source, name, side, event types, side, n)
	std::vector<AnaNumber> findMatching(const AnaNumber& AN);
	/// upload analysis result; sets timestamp and analysis number ID; optionally, over-write earlier results
	void uploadAnaNumber(AnaNumber& AN, bool replace = true);
		
protected:
	/// constructor (use AnalysisDB::getADB() if you need access to DB)
	AnalysisDB(const std::string& dbAddress = getEnvSafe("UCNAANARESDBADDRESS"),
			   const std::string& dbUser =  getEnvSafe("UCNADBUSER"),
			   const std::string& dbPass = getEnvSafe("UCNADBPASS"),
			   unsigned int port = atoi(getEnvSafe("UCNADBPORT","3306").c_str())
			   ): SQLHelper("analysis_results",dbAddress,dbUser,dbPass,port) {}
};

#endif
