#ifndef OCTET_HH
#define OCTET_HH 1

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

/// subdivision level for octet list
enum Octet_Division_t {
	DIV_OCTET	= 0,
	DIV_HALFOCT	= 1,
	DIV_PP		= 2,
	DIV_TRIAD	= 3,
	DIV_RUN		= 4
};
/// increment subdivision depth
inline Octet_Division_t& operator++(Octet_Division_t& d) { return d = Octet_Division_t(d+1); }
/// return next division depth
inline Octet_Division_t nextDiv(Octet_Division_t d) { Octet_Division_t d2 = d; return ++d2; }

/// class for octets and their subsets
class Octet {
public:
	/// constructor
	Octet(const Stringmap& m = Stringmap());
	/// split into (completed) sub-groupings
	std::vector<Octet> getSubdivs(Octet_Division_t divlvl = DIV_TRIAD, bool onlyComplete = false) const;
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
	
	Octet_Division_t divlevel;	//< dividion level of octet
	
protected:
	
	std::map< OctetRole,std::vector<RunNum> > runs;	//< list of octet runs
};

#endif
