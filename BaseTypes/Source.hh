#ifndef SOURCE_HH
#define SOURCE_HH 1

#include "SpectrumPeak.hh"
#include "QFile.hh"
#include <vector>

/// Class for representing a radioactive calibration source
class Source: public StringmapProvider {
public:
		
	/// constructor
	Source(std::string T = "", Side s = NOSIDE, unsigned int sid=0): t(T), mySide(s), x(0), y(0), wx(0), wy(0), nCounts(0), sID(sid) {}
	/// constructor from a StringMap
	Source(Stringmap S);
	
	/// get a list of peaks associated with this source type
	std::vector<SpectrumPeak> getPeaks() const;
		
	/// printable name for this source
	std::string name() const { return sID?t+"_"+itos(sID):t; }
	
	/// check whether a point is within this source's range
	inline bool inSourceRegion(float xx, float yy, float nSigma=3.0) const { return (xx-x)*(xx-x)/(wx*wx)+(yy-y)*(yy-y)/(wy*wy) < nSigma*nSigma; }
	
	std::string t;		//< this source's type
	Side mySide;		//< side on which this source appears
	RunNum myRun;		//< run number for this source
	float x;			//< x coordinate
	float y;			//< y coordinate
	float wx;			//< width in x
	float wy;			//< width in y
	float nCounts; 		//< number of counts in cut
	unsigned int sID;	//< unique ID number when source is constructed

protected:
	
	virtual Stringmap getProperties() const;
};

#endif
