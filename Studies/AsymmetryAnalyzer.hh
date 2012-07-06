#ifndef ASYMMETRYANALYZER_HH
#define ASYMMETRYANALYZER_HH 1

#include "OctetAnalyzer.hh"
#include "AnalysisDB.hh"

/// primary octet data analysis class
class AsymmetryAnalyzer: public OctetAnalyzer {
public:
	/// constructor
	AsymmetryAnalyzer(OutputManager* pnt, const std::string& nm, const std::string& inflname = "");
	
	/// cloning generator
	virtual SegmentSaver* makeAnalyzer(const std::string& nm,
									   const std::string& inflname) { return new AsymmetryAnalyzer(this,nm,inflname); }
	
	/// calculate super-ratio asymmetry from anode spectra
	virtual void calculateResults();
	/// output plot generation
	virtual void makePlots();
	/// upload results to analysis results DB
	virtual void uploadAnaResults();
	
	/// location of already-processed data (after first run) for errorbar estimation
	virtual std::string estimatorHistoLocation() const { return AsymmetryAnalyzer::processedLocation; }
	static std::string processedLocation;	//< set location here for already-processed files
	
	/// get base AnaResult to fill in
	AnaResult getResultBase() const;
	
	TH1F* hAsym;					//< asymmetry
	TH1F* hTpAsym[TYPE_III_EVENT];	//< asymmetry by event type
	TH1F* hInstAsym;				//< instrumental asymmetry
	TH1F* hSuperSum;				//< super-sum spectrum
	TH1F* hEvtSS[TYPE_III_EVENT];	//< super-sum for each event type
	
	quadHists qEnergySpectra[2][nBetaTubes+1][TYPE_IV_EVENT+1];	//< energy spectra quad hists for [side][tube][event type]
	quadHists qPositions[2][TYPE_III_EVENT+1];					//< event positions quad hists for [side][type]
	quadHists qAnodeCal[2];										//< anode calibration spectrum (Type 0, Erecon>225)
	
	quadHists qTotalSpectrum[2];								//< total spectrum based on analysis choice
	
protected:
	
	/// fill from scan data point
	virtual void fillCoreHists(ProcessedDataScanner& PDS, double weight);
	
	/// fit asymmetry over given range
	void fitAsym(float fmin, float fmax, unsigned int color, AnaResult AR, bool avg = false);
	/// fit instrumental asymmetry over given range
	void fitInstAsym(float fmin=200, float fmax=600, unsigned int color = 2);
	/// various beta spectrum endpoint fits
	void endpointFits();
	/// anode calibration fits
	void anodeCalFits();
	
	std::vector<AnaResult> asymFits;	//< list of asymmetry fits available for upload
	std::vector<AnaCutSpec> asymCuts;	//< list of cut specifications for asymmetry fits
	
	static TF1 asymmetryFit;		//< fit function for asymmetry
	static TF1 averagerFit;			//< pol0 line fit for averaging
	static AnalysisChoice anChoice;	//< asymmetry analysis choice
	
	TH1F* hEnergySpectra[2][nBetaTubes+1][TYPE_IV_EVENT+1];		//< energy spectra write point for [side][tube][event type]
	TH2F* hPositions[2][TYPE_III_EVENT+1];			//< event positions write point for [side][type]
	TH1F* hAnodeCal[2];								//< anode cal write point
};

#endif
