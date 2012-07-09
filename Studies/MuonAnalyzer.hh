#ifndef MUONANALYZER_HH
#define MUONANALYZER_HH 1

#include "OctetAnalyzer.hh"

/// muon data analyzer class
class MuonAnalyzer: public OctetAnalyzer {
public:
	/// constructor
	MuonAnalyzer(OutputManager* pnt, const std::string& nm, const std::string& inflname = "");
	
	/// cloning generator
	virtual SegmentSaver* makeAnalyzer(const std::string& nm,
									   const std::string& inflname) { return new MuonAnalyzer(this,nm,inflname); }
	
	/// output plot generation
	virtual void makePlots();
	
	/// location of already-processed data (after first run) for errorbar estimation
	virtual std::string estimatorHistoLocation() const { return MuonAnalyzer::processedLocation; }
	static std::string processedLocation;	//< set location here for already-processed files
	
	quadHists qMuonSpectra[2][2];			//< muon-veto event energy for [side][subtracted]
	quadHists qBackMuons[2][2];				//< back-veto tagged muon spectrum for [side][subtracted]
	fgbgPair* pMuonPos[2];					//< muon event positions
	fgbgPair* pBackMuPos[2];				//< backing veto muon event positions
	
protected:
	
	/// fill from scan data point
	virtual void fillCoreHists(ProcessedDataScanner& PDS, double weight);
	
	TH1F* hMuonSpectra[2][2];				//< muon spectra write point for [side][subtracted]
	TH1F* hBackMuons[2][2];					//< backing veto muon spectra write point for [side][subtracted]
};

#endif
