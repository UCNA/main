#ifndef POSITIONRESPONSE_HH
#define POSITIONRESPONSE_HH 1

#include "Enums.hh"
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
	std::vector<float> adc;		//< adc of Kurie plot endpoint
	std::vector<float> energy;	//< estimated energy of Kurie plot endpoint
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
class PositioningCorrector {
public:
	/// constructor, from input data
	PositioningCorrector(std::vector<PosmapInfo>& indat);
	/// constructor, from QFile
	PositioningCorrector(QFile& qin);
	/// get sector cutter
	SectorCutter& getSectors(Side s, unsigned int t) { assert((s==EAST||s==WEST) && t<tubes[s].size() && tubes[s][t]); return tubes[s][t]->S; }
	/// destructor
	~PositioningCorrector();	
	/// forbid copying
	PositioningCorrector& operator=(PositioningCorrector&) { assert(false); return *this; }
	/// get positioning correction for given tube
	double eval(Side s, unsigned int t, double x, double y, bool normalize = false) const;
	/// generate positioning information for each line in a file
	void processFile(const std::string& fInName, const std::string& fOutName) const;
	
private:
	/// init positioning interpolators
	void initPIs(std::vector<PosmapInfo>& indat);
	
	std::vector<PositioningInterpolator*> tubes[2];	//< interpolated position response maps for each tube
	std::vector<float> neta[2];						//< position map center normalization
};

#endif
