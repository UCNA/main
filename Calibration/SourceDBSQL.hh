#ifndef SOURCEDBSQL_HH
#define SOURCEDBSQL_HH 1

#include "Enums.hh"
#include "SQL_Utils.hh"
#include "PathUtils.hh"
#include "Source.hh"
#include "SpectrumPeak.hh"
#include <vector>

class SourceDBSQL: protected SQLHelper {
public:
	/// get singleton instance
	static SourceDBSQL* getSourceDBSQL();
	/// get list of sources for run (optionally for one side), sorted by x
	std::vector<Source> runSources(RunNum rn, Side s = BOTH);
	/// clear all run sources
	void clearSources(RunNum rn);
	/// add source to DB
	void addSource(const Source& src);
	/// delete peaks for source
	void clearPeaks(unsigned int sID);
	/// add peak to DB
	void addPeak(const SpectrumPeak& pk);
protected:
	/// constructor (use static GetSourceDBSQL for active instance)
	SourceDBSQL(const std::string& dbName,
				const std::string& dbAddress,
				const std::string& dbUser,
				const std::string& dbPass,
				unsigned int port): SQLHelper(dbName,dbAddress,dbUser,dbPass,port) {}
};

#endif
