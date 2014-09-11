#ifndef SIMASYMMETRYANALYZER_HH
#define SIMASYMMETRYANALYZER_HH

#include "OctetAnalyzer.hh"
#include "AsymmetryPlugin.hh"
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>

/// plugin with extra plots for simulated runs
class SimAsymmetryPlugin: public OctetAnalyzerPlugin {
public:
	/// constructor
	SimAsymmetryPlugin(OctetAnalyzer* OA);
        /// destructor
  //~SimAsymmetryPlugin() {SimOut->Print(); SimOut->Delete();}
	/// fill from scan data point
	virtual void fillCoreHists(ProcessedDataScanner& PDS, double weight);
	/// output plot generation
	virtual void makePlots();
	/// calculate corrections given actual data; return asymmetry correction stages
	std::vector<TH1*> calculateCorrections(AsymmetryPlugin& Adat, AsymmetryPlugin& Asim);
	/// "old style" corrections calculation
	std::vector<TH1*> calculateCorrectionsOld(AsymmetryPlugin& Asim);
	
	quadHists* qBCT[2][TYPE_IV_EVENT+1];		///< average beta cos theta TProfiles by [side][event type]
	quadHists* qCosth[2][TYPE_IV_EVENT+1];		///< average cos theta TProfiles by [side][type]
	quadHists* qBeta[2][TYPE_IV_EVENT+1];		///< average energy TProfiles by [side][type]
	quadHists* qAsymAcc[2][TYPE_IV_EVENT+1];	///< average asymmetry acceptance
	quadHists* qWrongSide[2][TYPE_III_EVENT+1];	///< Energy spectra of events ID'd on wrong side by [side][type]
	quadHists* qMissedSpectrum;					///< energy spectrum of missed events


  void makeSimTree();
  void fillSimTree(ProcessedDataScanner& PDS, double weight);
  //void printSimTree() { SimOut->Print(); }
  //Tree to hold simulated data which looks like real analyzed data
  TFile *SimFile;
  TTree *SimOut;
  TBranch *MWPCE, *MWPCW, *Erecon, *Type, *MWPCPosE, *MWPCPosW, *Event;
  Float_t mwpce, mwpcw, erecon;
  Float_t mwpcposE[2], mwpcposW[2];
  Int_t type;
};
#endif
