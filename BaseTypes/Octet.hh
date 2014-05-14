#ifndef OCTET_HH
#define OCTET_HH

#include "QFile.hh"
#include "strutils.hh"
#include <vector>
#include "Enums.hh"
#include <algorithm>

/// get AFP state for given octet type
AFPState afpForOctet(OctetRole t);
/// get GV state for given octet type
GVState gvForOctet(OctetRole t);
/// get run type for given octet type
RunType runTypeForOctet(OctetRole t);
/// get name for octet segment
std::string nameForOctet(OctetRole t);

/// class for octets and their subsets
class Octet {
public:
	/// constructor
	Octet(const Stringmap& m = Stringmap());
	/// split into (completed) sub-groupings
	std::vector<Octet> getSubdivs(RunGrouping divlvl, bool onlyComplete = false) const;
	/// add run to octet listing
	void addRun(RunNum rn, OctetRole t);
	/// get runs of given type
	std::vector<RunNum> getRuns(OctetRole t) const;
	/// get asymmetry runs list
	std::vector<RunNum> getAsymRuns(bool foreground) const;
	/// get sorted list of all runs
	std::vector<RunNum> getAllRuns() const;
	/// get number of runs
	unsigned int getNRuns() const { return runs.size(); }
	
	/// load all octets from QFile
	static std::vector<Octet> loadOctets(const QFile& q);
	/// load numbered octet from QFile
	static Octet loadOctet(const QFile& q, unsigned int n);
	/// convert to stringMap
	Stringmap toStringmap() const;
	
	/// get AFP state of octet, AFP_OTHER if mixed
	AFPState octAFPState() const;
	/// get name for this octet
	std::string octName() const;
	/// get first run, for sorting Octets
	RunNum getFirstRun() const;
	
	RunGrouping grouping;	///< dividion level of octet
	
protected:
	
	std::map< OctetRole,std::vector<RunNum> > runs;	///< list of octet runs
};

#endif
