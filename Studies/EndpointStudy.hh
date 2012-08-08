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

/// analyzer plugin for wirechamber position/calibration analysis
class WirechamberAnalyzer: public AnalyzerPlugin {
public:
	/// constructor
	WirechamberAnalyzer(RunAccumulator* RA);
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

/// analyzer plugin for position-segmented analysis
class PositionBinner: public AnalyzerPlugin {
public:
	/// constructor
	PositionBinner(RunAccumulator* RA, unsigned int nr = 0);
	
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
	virtual void compareMCtoData(AnalyzerPlugin* AP);
	
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

/// analyzer for xenon data
class XenonAnalyzer: public RunAccumulator {
public:
	/// constructor
	XenonAnalyzer(OutputManager* pnt, const std::string& nm = "XenonAnalyzer", unsigned int nr = 0, const std::string& inflName = ""):
	RunAccumulator(pnt,nm,inflName) {
		addPlugin(myEndpt = new PositionBinner(this,nr));
		addPlugin(myWC = new WirechamberAnalyzer(this));
	}
	/// create a new instance of this object (cloning self settings) for given directory
	virtual SegmentSaver* makeAnalyzer(const std::string& nm, const std::string& inflname) {
		return new XenonAnalyzer(this,nm,myEndpt->sects.n,inflname); }
	
	PositionBinner* myEndpt;	//< xenon endpoint map plugin
	WirechamberAnalyzer* myWC;	//< wirechamber positions plugin
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
