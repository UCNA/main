#ifndef POSITIONRESPONSE_HH
#define POSITIONRESPONSE_HH 1

#include "Enums.hh"
#include "Types.hh"
#include "Interpolator.hh"
#include "SectorCutter.hh"
#include "QFile.hh"
#include "strutils.hh"
#include <vector>
#include <cassert>
#include <TGraph.h>
#include <fstream>

/// info for building a positioning interpolator
struct PosmapInfo {
	unsigned int nRings;		//< number of rings
	float radius;				//< radius of sampled area
	Side s;						//< side this calibration applies to
	unsigned int t;				//< tube number for this calibration
	std::vector<float> signal;	//< light signal at each position
	std::vector<float> norm;	//< normalization at each position (accounting for energy resolution smearing, etc.)
};

/// r-theta interpolator for beta tube position response map
class PositioningInterpolator {
public:
	
	/// constructor, from input data
	PositioningInterpolator(const PosmapInfo& PMI);
	/// destructor
	~PositioningInterpolator();	
	/// forbid copying
	PositioningInterpolator& operator=(PositioningInterpolator&) { assert(false); return *this; }
	/// evaluate from interpolation table
	double eval(double x, double y);	
	
	SectorCutter S;
	
protected:	
	InterpoSequence sRadial;
	CubiTerpolator L;
	//Interpolator L;
	std::vector<DoubleSequence*> phiSeqs;
	std::vector<CubiTerpolator*> phiInterps;
	//std::vector<Interpolator*> phiInterps;
};

/// positioning interpolators for each tube
class PositioningCorrector: private NoCopy {
public:
	/// constructor, from input data
	PositioningCorrector(std::vector<PosmapInfo>& indat);
	/// constructor, from QFile
	PositioningCorrector(QFile& qin);
	/// destructor
	~PositioningCorrector();
	
	/// get sector cutter
	SectorCutter& getSectors(Side s, unsigned int t) { assert((s==EAST||s==WEST) && t<tubes[s].size() && tubes[s][t]); return tubes[s][t]->S; }
	/// get positioning correction for given tube
	double eval(Side s, unsigned int t, double x, double y, bool normalize = false) const;
	/// set normalization to c at center
	void setNormCenter(double c = 1.0);
	/// set normalization to average of c
	void setNormAvg(double c = 1.0);
	
	/// get number of maps available on side
	unsigned int getNMaps(Side s) const { assert(s==EAST||s==WEST); return tubes[s].size(); }
	
private:
	/// init positioning interpolators
	void initPIs(std::vector<PosmapInfo>& indat);
	
	std::vector<PositioningInterpolator*> tubes[2];	//< interpolated position response maps for each tube
	std::vector<float> neta[2];						//< position map center normalization
};

#endif
