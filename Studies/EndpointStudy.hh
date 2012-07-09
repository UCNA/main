#ifndef ENDPOINTSTUDY_HH
#define ENDPOINTSTUDY_HH 1

#include "SectorCutter.hh"
#include "RunAccumulator.hh"
#include "WirechamberCalibrator.hh"


/// Struct for cathode calibration data
struct CathodeSeg {
	Side s;					//< side
	AxisDirection d;		//< plane measuring direction
	unsigned int i;			//< cathode number
	float_err height;		//< normalized height
	float_err width;		//< width
	float_err center;		//< measured center
	float_err max;			//< maximum normalized value
	float fill_frac;		//< proportion of expected events
	float pos;				//< cathode position
};

/// wirechamber position/calibration analysis
class WirechamberAnalyzer: public RunAccumulator {
public:
	/// constructor
	WirechamberAnalyzer(OutputManager* pnt, const std::string& nm, const std::string& infl="");
	
	/// cloning generator: just return another of the same subclass (with any settings you want to preserve)
	virtual SegmentSaver* makeAnalyzer(const std::string& nm,
									   const std::string& inflname) { return new WirechamberAnalyzer(this,nm,inflname); }
	
	/// process a data point into position histograms
	virtual void fillCoreHists(ProcessedDataScanner& PDS, double weight);
	/// make output plots
	virtual void makePlots();
	
	fgbgPair* hitPos[2];						//< hit positions on each side
	fgbgPair* hitPosRaw[2];						//< uncorrected hit positions on each side
	fgbgPair* cathHitpos[2][2][kMaxCathodes];	//< raw position distribution around each cathode by energy
	fgbgPair* cathNorm[2][2][kMaxCathodes];		//< cathode normalization histograms
};

/// data collected for each sector
struct SectorDat {
	Side s;
	unsigned int t;
	unsigned int m;
	float eta;
	float_err low_peak;
	float_err low_peak_width;
	float_err xe_ep;
};

/// class for position-segmented analysis
class PositionBinner: public WirechamberAnalyzer {
public:
	/// constructor
	PositionBinner(OutputManager* pnt, const std::string& nm, unsigned int nr, const std::string& infl="");
	
	/// cloning generator: just return another of the same subclass (with any settings you want to preserve)
	virtual SegmentSaver* makeAnalyzer(const std::string& nm,
									   const std::string& inflname) { return new PositionBinner(this,nm,sects.n,inflname); }
	
	/// process a data point into position histograms
	virtual void fillCoreHists(ProcessedDataScanner& PDS, double weight);
	
	/// overall spectrum info
	virtual void calculateResults();
	/// fit enpoints in each sector
	void fitSectors();
	/// fit peak and endpoint for energy spectrum
	void fitSpectrum(TH1* hSpec,SectorDat& sd);
	/// make output plots
	virtual void makePlots();
	/// virtual routine for MC/Data comparison plots/calculations
	virtual void compareMCtoData(RunAccumulator& OAdata);
	
	fgbgPair* energySpectrum;							//< total energy spectrum
	TH1F* hTuben[2][nBetaTubes];						//< type 0 PMT energy for gain matching
	std::vector<SectorDat> sectDat[2][nBetaTubes];		//< processed data for each sector
	SectorCutter sects;									//< sector cutter
	bool sectorPlots;									//< whether to make plots for each sector
	
	static double fidRadius;							//< analysis fiducial radius
	
protected:
	std::vector<fgbgPair*> sectEnergy[2][nBetaTubes];	//< visible light histograms by position, PMT
	std::vector<fgbgPair*> sectAnode[2];				//< anode histograms by position
};

/// process xenon runs
void process_xenon(RunNum r0, RunNum r1, unsigned int nrings);

/// generate xenon simulation runs
void simulate_xenon(RunNum r0, RunNum r1, RunNum rsingle=0, unsigned int nRings =12);

/// comparison between data, simulated Xenon to produce position map
void xenon_posmap(RunNum r0, RunNum r1, unsigned int nRings);


/// process wirechamber calibraiton data
void processWirechamberCal(WirechamberAnalyzer& WCdat, WirechamberAnalyzer& WCsim);
/// process wirechamber calibraiton data specified by run range, nrings
void processWirechamberCal(RunNum r0, RunNum r1, unsigned int nrings);

#endif
