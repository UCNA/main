#ifndef ENDPOINTSTUDY_HH
#define ENDPOINTSTUDY_HH 1

#include "SectorCutter.hh"
#include "RunAccumulator.hh"

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
class PositionBinner: public RunAccumulator {
public:
	/// constructor
	PositionBinner(OutputManager* pnt, const std::string& nm, float r, unsigned int nr, const std::string& infl="");
	
	/// cloning generator: just return another of the same subclass (with any settings you want to preserve)
	virtual SegmentSaver* makeAnalyzer(const std::string& nm,
									   const std::string& inflname) { return new PositionBinner(this,nm,sects.r,sects.n,inflname); }
	
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
	
	fgbgPair energySpectrum;							//< total energy spectrum
	fgbgPair hitPos[2];									//< hit positions on each side
	TH1F* hTuben[2][nBetaTubes];						//< type 0 PMT energy for gain matching
	std::vector<SectorDat> sectDat[2][nBetaTubes];		//< processed data for each sector
	SectorCutter sects;									//< sector cutter
	bool sectorPlots;									//< whether to make plots for each sector
	
protected:
	std::vector<fgbgPair> sectEnergy[2][nBetaTubes];	//< visible light histograms by position, PMT
	std::vector<fgbgPair> sectAnode[2];					//< anode histograms by position
};

/// process xenon runs
void process_xenon(RunNum r0, RunNum r1, unsigned int nrings);

/// compare xenon to simulation
void simulate_xenon(RunNum r0, RunNum r1, RunNum rsingle=0, unsigned int nRings =12);
 
#endif
