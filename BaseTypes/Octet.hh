#ifndef OCTET_HH
#define OCTET_HH 1

#include "QFile.hh"
#include "strutils.hh"
#include <vector>
#include "Enums.hh"
#include <algorithm>

/// get AFP state for given octet type
AFPState afpForOctet(OctetType t);
/// get GV state for given octet type
GVState gvForOctet(OctetType t);
/// get run type for given octet type
RunType runTypeForOctet(OctetType t);
/// get triad grouping for given octet type
TriadType triadForOctet(OctetType t);
/// get name for octet segment
std::string nameForOctet(OctetType t);

/// class for octets and their subsets
class Octet {
public:
	/// constructor
	Octet(const Stringmap& m = Stringmap());
	/// split into (completed) sub-groupings
	std::vector<Octet> getSubdivs(unsigned int divlvl = 3, bool onlyComplete = false) const;
	/// add run to octet listing
	void addRun(RunNum rn, OctetType t);
	/// get runs of given type
	std::vector<RunNum> getRuns(OctetType t) const;
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
	
	unsigned int divlevel;	//< dividion level of octet
	
protected:
	
	std::map< OctetType,std::vector<RunNum> > runs;	//< list of octet runs
};

#endif
