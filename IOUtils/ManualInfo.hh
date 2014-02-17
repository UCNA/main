#ifndef MANUALINFO_HH
#define MANUALINFO_HH 1

#include "QFile.hh"
#include "Types.hh"

/// class for looking up manually-stored miscellaneous analysis info
class ManualInfo: public QFile {
public:
	/// constructor
	ManualInfo(std::string fname): QFile(fname) {}
	
	/// get (double start,double end) pairs for key (or any other named pairs)
	std::vector< std::pair<double,double> > getRanges(const std::string& key, const std::string& k1="start", const std::string& k2="end") const;
	
	/// get matching keys for item in range
	std::vector<Stringmap> getInRange(const std::string& key,
									  const double x,
									  const std::string& k1="runStart",
									  const std::string& k2="runEnd") const;
	
	static ManualInfo MI;	//< static global instance to use
};

/// simple class for cuts from Stringmap
class RangeCut {
public:
	/// constructor
	RangeCut(const Stringmap& m = Stringmap());
	/// constructor with start, end times
	RangeCut(double s, double e): start(s), end(e) {}
	
	/// check if value is in range
	inline bool inRange(double x) const { return start <= x && x <= end; }
	
	double start;	//< cut minimum
	double end;		//< cut maximum
};

/// simple class for value + cuts range
class CutVariable {
public:
	/// constructor
	CutVariable(std::string sn=""): sname(sn) {}
	/// check if in range
	inline bool inRange() const { return R.inRange(val); }
	std::string sname;	//< sensor name (for pedestal subtraction)
	Float_t val;		//< stored value
	RangeCut R;			//< cuts range
};

Stringmap loadCut(RunNum rn, const std::string& cutName);

void loadRangeCut(RunNum rn, CutVariable& c, const std::string& cutName);

#endif
