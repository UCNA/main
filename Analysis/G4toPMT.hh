#ifndef G4TOPMT_HH
#define G4TOPMT_HH 1

#include "Sim2PMT.hh"
#define N_SD 24

/// converts Geant 4 simulation results to PMT spectra
class G4toPMT: public Sim2PMT {
public:
	/// constructor
	G4toPMT(bool ext = false): Sim2PMT("anaTree"), extended(ext) { }
	/// unit conversions
	virtual void doUnits();
		
protected:
	/// set read points for input tree
	virtual void setReadpoints();
	
	bool extended;			//< whether to read in additional variables
	double eDepSD[N_SD];	//< energy deposition array
	double thetaInSD[N_SD];	//< entrance angle array
	double thetaOutSD[N_SD];//< exit angle array
	double keInSD[N_SD];	//< entrance energy array
	double keOutSD[N_SD];	//< exit energy array
};

/// For consistency checks, swaps E/W sides on Geant4 sim data
class G4toPMT_SideSwap: public G4toPMT {
public:
	/// constructor
	G4toPMT_SideSwap(): G4toPMT() { }
	/// unit conversions
	virtual void doUnits();
};

/// multiply Geant4 data to fill all segments of SectionCutter
class G4SegmentMultiplier: public G4toPMT, SimPositioner {
public:
	/// constructor
	G4SegmentMultiplier(const SectorCutter& S);
	/// start scan
	virtual void startScan(bool startRandom = false);
	/// overrides G4toPMT::nextPoint to re-use points
	virtual bool nextPoint();
	/// unit conversions, only done once per rotation sequence
	virtual void doUnits();
	/// calculate spectrum re-weighting factor
	virtual void calcReweight();
	/// apply rotational offset
	virtual void applyOffset(Sim2PMT& S);
protected:
	/// rotate a point
	void rotpt(double& x0, double& y0);
	
	SectorCutter SC;			//< SectorCutter to determine event multiplication
	unsigned int nrots;			//< number of remaining point rotations
	unsigned int rcurrent;		//< current radial ring
	std::vector<double> vc,vs;	//< pre-calculated rotation matrix cosines, sines for each ring
	bool morePts;				//< whether the data still has more points to come
};

#endif
