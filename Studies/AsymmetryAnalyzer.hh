#ifndef ASYMMETRYANALYZER_HH
#define ASYMMETRYANALYZER_HH 1

#include "OctetAnalyzer.hh"

/// primary octet data analysis class
class AsymmetryAnalyzer: public OctetAnalyzer {
public:
	/// constructor
	AsymmetryAnalyzer(OutputManager* pnt, const std::string& nm, const std::string& inflname = "");
	
	/// cloning generator
	virtual SegmentSaver* makeAnalyzer(const std::string& nm,
									   const std::string& inflname) { return new AsymmetryAnalyzer(this,nm,inflname); }
	
	/// MC/data comparison
	void compareMCtoData(RunAccumulator& OAdata);
	
	/// calculate super-ratio asymmetry from anode spectra
	virtual void calculateResults();
	/// output plot generation
	virtual void makePlots();
	
	/// location of already-processed data (after first run) for errorbar estimation
	virtual std::string estimatorHistoLocation() const { return AsymmetryAnalyzer::processedLocation; }
	static std::string processedLocation;	//< set location here for already-processed files
	
	TH1F* hAsym;					//< asymmetry
	TH1F* hInstAsym;				//< instrumental asymmetry
	TH1F* hSuperSum;				//< super-sum spectrum
	TH1F* hEvtSS[TYPE_III_EVENT];	//< super-sum for each event type
	
protected:
	
	/// fill from scan data point
	virtual void fillCoreHists(ProcessedDataScanner& PDS, double weight);
	
	/// fit asymmetry over given range
	void fitAsym(float fmin, float fmax, unsigned int color, bool avg = false);
	/// fit instrumental asymmetry over given range
	void fitInstAsym(float fmin=200, float fmax=600, unsigned int color = 2);
	/// various beta spectrum endpoint fits
	void endpointFits();
	/// anode calibration fits
	void anodeCalFits();
	
	static TF1 asymmetryFit;		//< fit function for asymmetry
	static TF1 averagerFit;			//< pol0 line fit for averaging
	static AnalysisChoice anChoice;	//< asymmetry analysis choice
	
	quadHists qEnergySpectra[2][nBetaTubes+1][TYPE_IV_EVENT+1];	//< energy spectra quad hists for [side][tube][event type]
	TH1F* hEnergySpectra[2][nBetaTubes+1][TYPE_IV_EVENT+1];		//< energy spectra write point for [side][tube][event type]
	quadHists qPositions[2][TYPE_III_EVENT+1];		//< event positions quad hists for [side][type]
	TH2F* hPositions[2][TYPE_III_EVENT+1];			//< event positions write point for [side][type]
	quadHists qAnodeCal[2];							//< anode calibration spectrum (Type 0, Erecon>225)
	TH1F* hAnodeCal[2];								//< anode cal write point
};

#endif
