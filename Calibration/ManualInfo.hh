#ifndef MANUALINFO_HH
#define MANUALINFO_HH 1

#include "QFile.hh"

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

#endif
