#ifndef SIMNONLINEARITY_HH
#define SIMNONLINEARITY_HH 1

#include "strutils.hh"
#include "Enums.hh"
#include "QFile.hh"
#include <TGraph.h>

/// class for simulating nonlinear PMT signals
class SimNonlinearizer {
public:
	/// constructor
	SimNonlinearizer();
	
	/// apply nonlinearity to observed light
	float delinearize(Side s, unsigned int t, float l0) const;
	
	/// stringmap description
	Stringmap toStringmap() const;
	
	/// fix calibration point by stretching
	void fixCalPoint(Side s, unsigned int t, float l0, float l1 = 10000);
	
	/// make random error within limits
	void makeRanderr(bool correlate = false);
	
	/// set energy offset nonlinearity
	void setOffset(float dl, Side s0 = BOTH);
	/// generate random PMT offsets
	void randomOffsets(float sigma, float mu = 0, bool gaussStats=true);
	/// generate random gain fluctuations
	void randomGainflucts(float sigma, float mu = 0, bool gaussStats=true);
	/// random n-term sum of relative fluctuations
	void randomRelcurve(float sigmaRel, unsigned int nTerms, bool gaussStats=true);
	/// random n-term sum of absolute fluctuations
	void randomAbscurve(float sigma, unsigned int nTerms, bool gaussStats=true);
	/// add noise to gain factors
	void addGainNoise(float sigmaRel);
	/// set relative nonlinearity
	void setRelerr(float relerr, Side s0 = BOTH);
	/// set resolution error
	void setReserr(float reserr, Side s0 = BOTH);
	/// notify of new run start
	void startNewRun(AFPState afp);
	/// set all PMTs on each side to same error
	void unifyErrors();
	/// correlate errors between sides
	void correlateErrors(bool correlate = true);
	
	/// check whether passes error limits
	bool checkLimits() const;
	
	/// get maximum energy error over energy range
	float maxError(Side s, unsigned int t, float emin, float emax) const;
	/// get maximum energy error over energy range, all sides/tubes
	float maxError(float emin, float emax) const;
	
	float gainfactor[2][nBetaTubes];				//< final gain correction coefficients
	float newrunGainerr;							//< gain error to apply at start of new run
	float afpCorrGain;								//< AFP-correlated gain factor
	float resErr[2][nBetaTubes];					//< energy resolution error factor

protected:
	float range;									//< nonlinearity function range
	TGraph* errlim;									//< limits on acceptable error
	std::vector<float> abscoeffs[2][nBetaTubes];	//< coefficients for absolute offsets
	std::vector<float> relcoeffs[2][nBetaTubes];	//< coefficients for relative offsets
};

#endif
