#ifndef ANALYSISDB_HH
#define ANALYSISDB_HH 1

#include "SQL_Utils.hh"
#include "Enums.hh"
#include "Types.hh"
#include "PathUtils.hh"
#include "SMExcept.hh"
#include <string>
#include <set>
#include <time.h>

/// description of an analysis cut
class AnaCutSpec {
public:
	/// constructor
	AnaCutSpec(): csid(0), postp(AnaCutSpec::POS_PLAIN) {}
	
	unsigned int csid;	//< ID number
	double emin;		//< minimum energy
	double emax;		//< maximum energy
	double radius;		//< maximum radius
	enum PosType {
		POS_PLAIN,		//< plain positions
		POS_ROTATED		//< positions corrected for rotation
	} postp;			//< positioning method used
};

/// description of an analysis result
class AnaResult {
public:
	/// constructor
	AnaResult(const std::string& auth = getEnvSafe("UCNA_ANA_AUTHOR"));
	/// return as Stringmap
	Stringmap toStringmap() const;
	/// string representation of types set for DB
	std::string typeSetString() const;
	
	unsigned int arid;		//< analysis result ID number in database
	std::string	author;		//< author of analysis result
	time_t timestamp;		//< timestamp of analysis result
	enum AnaType {
		ANA_ASYM,			//< asymmetry analysis
		ANA_COUNTS			//< counts analysis
	} anatp;				//< type of analysis result
	enum DataSource {
		REAL_DATA,			//< actual replay data
		G4_DATA,			//< Geant4 MC
		PEN_DATA			//< Penelope MC
	} datp;					//< data/MC source for analysis
	enum RunGrouping {
		GROUP_RUN=0,		//< single run (or grouped runs with same purpose)
		GROUP_PAIR=1,		//< AFP On/Off pulse pair
		GROUP_QUARTET=2,	//< A or B pair of pulse pairs
		GROUP_OCTET=3		//< A->B or B->A octet
	} grouping;
	RunNum startRun;		//< start of run range analyzed
	RunNum endRun;			//< end of run range analyzed
	std::set<EventType> etypes;	//< event types considered
	AnalysisChoice anach;	//< ``analysis choice'' used
	Side s;					//< side for result
	AFPState afp;			//< AFP state for result
	double value;			//< value of result
	double err;				//< uncertainty on result
	unsigned int csid;		//< cut specification ID
	
	/// words for analysis type
	static std::string atypeWord(AnaType t);
	/// analysis type from word
	static AnaType strToAtype(const std::string& str);
	/// words for data source
	static std::string dsourceWord(DataSource d);
	/// words for octet grouping
	static std::string groupWord(RunGrouping g);
	/// data type from word
	static DataSource strToDsource(const std::string& str);
};

/// DB interface for sharing analysis results
class AnalysisDB: public SQLHelper, NoCopy {
public:
	/// globally available AnalysisDB
	static AnalysisDB* getADB();
	
	/// upload a cut spec and return ID
	unsigned int uploadCutSpec(AnaCutSpec& c);
	/// upload an analysis result and return ID
	unsigned int uploadAnaResult(AnaResult& r);
	/// delete analysis result by ID
	void deleteAnaResult(unsigned int arid);
	/// delete cut spec by ID
	void deleteCutSpec(unsigned int csid);
	/// get analysis result by ID
	AnaResult getAnaResult(unsigned int arid);
	/// get cut spec by ID
	AnaCutSpec getCutSpec(unsigned int csid);
	
	/// find analysis results matching prototype result
	std::vector<AnaResult> findMatching(const AnaResult& A);
	
protected:
	/// constructor (use AnalysisDB::getADB() if you need access to DB)
	AnalysisDB(const std::string& dbAddress = getEnvSafe("UCNAANARESDBADDRESS"),
			   const std::string& dbUser =  getEnvSafe("UCNADBUSER"),
			   const std::string& dbPass = getEnvSafe("UCNADBPASS"),
			   unsigned int port = atoi(getEnvSafe("UCNADBPORT","3306").c_str())
			   ): SQLHelper("analysis_results",dbAddress,dbUser,dbPass,port) {}
	/// string for type set
	static std::string typeSetString(const std::set<EventType>& tps);
};

#endif
