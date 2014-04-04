#ifndef XENONANALYZER_HH
#define XENONANALYZER_HH

#include "PositionBinnedPlugin.hh"
#include "WirechamberGainMapPlugins.hh"
#include "CathodeTuningAnalyzer.hh"
#include "WirechamberEnergyPlugins.hh"

/// xenon data collected for each sector
struct SectorDat {
	Side s;
	unsigned int t;
	unsigned int m;
	float eta;
	float_err low_peak;
	float_err low_peak_width;
	float_err xe_ep;
};
/// convert Stringmap to SectorDat
SectorDat sm2sd(const Stringmap& m);
/// convert SectorDat to Stringmap
Stringmap sd2sm(const SectorDat& sd);

/// analyzer plugin for Xenon spectrum shape fitting and position map
class XenonSpectrumPlugin: public PositionBinnedPlugin {
public:
	/// constructor
	XenonSpectrumPlugin(RunAccumulator* RA, unsigned int nr = 0);
	
	/// process a data point into position histograms
	virtual void fillCoreHists(ProcessedDataScanner& PDS, double weight);
	/// overall spectrum info
	virtual void calculateResults();
	/// generate position map from data endpoints
	void compareMCtoData(AnalyzerPlugin* AP);
	/// fit endpoint in each sector
	void fitSectors();
	
	fgbgPair* energySpectrum;							//< overall energy spectrum for isotope decomposition
	std::vector<fgbgPair*> sectEnergy[2][nBetaTubes+1];	//< visible light histograms by position, PMT
	std::vector<SectorDat> sectDat[2][nBetaTubes+1];	//< processed data for each sector
	
protected:
	
	/// fit a xenon spectrum
	void fitSpectrum(TH1* hSpec,SectorDat& sd);
};

/// analyzer for xenon data
class XenonAnalyzer: public CathodeTuningAnalyzer {
public:
	/// constructor
	XenonAnalyzer(OutputManager* pnt, const std::string& nm, const std::string& inflName = "", unsigned int nrE = 0, unsigned int nrA = 12);
	/// create a new instance of this analyzer
	virtual SegmentSaver* makeAnalyzer(const std::string& nm, const std::string& inflname) {
		return new XenonAnalyzer(this,nm,inflname,myXeSpec->sects.n,myAnode->sects.n); }
	
	XenonSpectrumPlugin* myXeSpec;	//< position-binned Xenon spectrum analysis
	AnodeGainMapPlugin* myAnode;	//< anode position gain
	MWPCGainPlugin* myWG;			//< MWPC energy calibration
};

/// analyzer for simulated xenon data
class SimXenonAnalyzer: public CathodeTuningAnalyzer {
public:
	/// constructor
	SimXenonAnalyzer(OutputManager* pnt, const std::string& nm, const std::string& inflName = "", unsigned int nrE = 0);
	
	XenonSpectrumPlugin* myXeSpec;	//< position-binned Xenon spectrum analysis
	MWPCGainPlugin* myWG;			//< MWPC energy calibration
};

/// process xenon runs
void process_xenon(RunNum r0, RunNum r1, unsigned int nrings);
/// generate xenon simulation runs
void simulate_xenon(RunNum r0, RunNum r1, RunNum rsingle=0, unsigned int nRings =12);
/// comparison between data, simulated Xenon to produce position map
void xenon_posmap(RunNum r0, RunNum r1, unsigned int nRings);

#endif
