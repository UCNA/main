#ifndef ASYMMETRYANALYZER_HH
#define ASYMMETRYANALYZER_HH 1

#include "OctetAnalyzer.hh"
#include "MuonAnalyzer.hh"
#include "AnalysisDB.hh"

/// energy spectra and asymmetry analysis class
class AsymmetryAnalyzer: public OctetAnalyzerPlugin {
public:
	/// constructor
	AsymmetryAnalyzer(OctetAnalyzer* OA);
	
	/// fill from scan data point
	virtual void fillCoreHists(ProcessedDataScanner& PDS, double weight);
	/// calculate super-ratio asymmetry from anode spectra
	virtual void calculateResults();
	/// combine events used to form asymmetry, super-sums
	void calcSuperCombos();
	/// output plot generation
	virtual void makePlots();
	/// upload results to analysis results DB
	virtual void uploadAnaResults();
	/// MC/data comparison
	virtual void compareMCtoData(AnalyzerPlugin* AP);
	
	/// get base AnaResult to fill in
	AnaResult getResultBase() const;
	
	unsigned int nEnergyBins;		//< number of bins for energy histograms
	double energyMax;				//< energy range for energy histograms
	
	TH1F* hAsym;						//< asymmetry
	TH1F* hTpAsym[TYPE_III_EVENT+1];	//< asymmetry by event type
	TH1F* hInstAsym;					//< instrumental asymmetry
	TH1F* hSuperSum;					//< super-sum spectrum
	TH1F* hEvtSS[TYPE_III_EVENT+1];		//< super-sum for each event type
	
	quadHists* qEnergySpectra[2][nBetaTubes+1][TYPE_IV_EVENT+1];	//< energy spectra quad hists for [side][tube][event type]
	quadHists* q23ProbCut[2][TYPE_III_EVENT+1];						//< Type II/III spectra by probability cut for [side][type]
	quadHists* qTotalSpectrum[2];									//< total spectrum based on analysis choice
	
	AnalysisChoice anChoice;	//< asymmetry analysis choice

protected:
	
	/// fit asymmetry over given range
	void fitAsym(float fmin, float fmax, unsigned int color, AnaResult AR, bool avg = false);
	/// fit instrumental asymmetry over given range
	void fitInstAsym(float fmin=200, float fmax=600, unsigned int color = 2);
	/// various beta spectrum endpoint fits
	void endpointFits();
	
	AnaResult ARtot;					//< temporary analysis result holder
	std::vector<AnaResult> asymFits;	//< list of asymmetry fits available for upload
	std::vector<AnaCutSpec> asymCuts;	//< list of cut specifications for asymmetry fits
	
	static TF1 asymmetryFit;		//< fit function for asymmetry
	static TF1 averagerFit;			//< pol0 line fit for averaging
};

#endif
