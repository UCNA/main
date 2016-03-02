#ifndef POSITIONRESPONSE_HH
#define POSITIONRESPONSE_HH

#include "Enums.hh"
#include "Types.hh"
#include "Interpolator.hh"
#include "SectorCutter.hh"
#include "QFile.hh"
#include "strutils.hh"
#include <vector>
#include "SMExcept.hh"
#include <TGraph.h>
#include <fstream>

/// info for building a positioning interpolator
struct PosmapInfo {
	unsigned int nRings;		///< number of rings
	float radius;				///< radius of sampled area
	Side s;						///< side this calibration applies to
	unsigned int t;				///< tube number for this calibration
	std::vector<float> signal;	///< light signal at each position
	std::vector<float> norm;	///< normalization at each position (accounting for energy resolution smearing, etc.)
};

/// StringMap info for SectorCutter
Stringmap SCtoSM(const SectorCutter& SC);

/// r-theta interpolator for beta tube position response map
class PositioningInterpolator {
public:
	
	/// constructor, from input data
	PositioningInterpolator(const PosmapInfo& PMI,
		Interpolator* (*phiInterp)(DataSequence*, double, double) = &CubiTerpolator::newCubiTerpolator,
		Interpolator* (*rInterp)(DataSequence*, double, double) = &CubiTerpolator::newCubiTerpolator);
	/// destructor
	~PositioningInterpolator();	
	/// forbid copying
	PositioningInterpolator& operator=(PositioningInterpolator&) { smassert(false); return *this; }
	/// evaluate from interpolation table
	double eval(double x, double y);
        ///evaluate straight from position map
        double evalNoInterp(double x, double y);
	
        const PosmapInfo pmi;
	SectorCutter S;
        
	
protected:	
	InterpoSequence sRadial;
	Interpolator* L;
	std::vector<DoubleSequence*> phiSeqs;
	std::vector<Interpolator*> phiInterps;
	//std::vector<Interpolator*> phiInterps;
};

/// positioning interpolators for each tube
class PositioningCorrector: private NoCopy {
public:
	/// constructor
	PositioningCorrector(): interpType(PositioningCorrector::defaultInterpType) { }
	/// destructor
	~PositioningCorrector();
	
	/// load input data to define interpolators
	void loadData(const std::vector<PosmapInfo>& indat);
	/// load input data from QFile
	void loadData(const QFile& qin);
	/// get raw posmap data
	const std::vector<PosmapInfo>& getData() const { return myData; }
	
	/// get sector cutter
	const SectorCutter& getSectors(Side s, unsigned int t) const { smassert((s==EAST||s==WEST) && t<tubes[s].size() && tubes[s][t]); return tubes[s][t]->S; }
	/// get positioning correction for given tube
	double eval(Side s, unsigned int t, double x, double y, bool normalize = false) const;
        /// get positioning correction for given tube not interpolated
	double evalNoInterp(Side s, unsigned int t, double x, double y, bool normalize = false) const;
	/// set normalization to c at center
	void setNormCenter(double c = 1.0);
	/// set normalization to average of c
	void setNormAvg(double c = 1.0);
	
	/// get number of maps available on side
	unsigned int getNMaps(Side s) const { smassert(s==EAST||s==WEST); return tubes[s].size(); }
	
	/// interpolation type
	Interpolator* (*interpType)(DataSequence*, double, double);
	/// default interpolation type
	static Interpolator* (*defaultInterpType)(DataSequence*, double, double);
	
private:
	/// delete existing interpolators
	void deleteInterpolators();
	std::vector<PosmapInfo> myData;					///< position map building data
	std::vector<PositioningInterpolator*> tubes[2];	///< interpolated position response maps for each tube
	std::vector<float> neta[2];						///< position map center normalization
};

#endif
