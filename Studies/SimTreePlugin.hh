#ifndef SIMTREEANALYZER_HH
#define SIMTREEANALYZER_HH


#include "OctetAnalyzer.hh"
#include "AsymmetryPlugin.hh"
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>

// Data structures for the two different branches

// Variables which come directly from the simulation (NOTE: 0 = E, 1 = W in the array index of the leaves)
struct simData {
  Double_t eQ[2];						///< Scintillator quenched energy [keV]
  Double_t eDep[2];						///< Scintillator deposited energy [keV]
  Double_t eW[2];						///< Wirechamber active volume deposited energy [keV]
  Double_t scintPos[2][3];	///< hit position in scintillator [mm in decay trap]
  Double_t mwpcPos[2][3];	///< hit position in MWPC
  Double_t primPos[4];			///< primary event vertex position (4=radius)
  Double_t time[2];						///< hit time [s] in each scintillator
  Double_t costheta;						///< primary event cos pitch angle
  Double_t ePrim;							///< primary event energy
  Double_t edepFoils[2];			///< energy deposition [keV] in decay trap foils
  Double_t edepWinOut[2];		///< energy deposition [keV] in outer wirechamber window
  Double_t edepWinIn[2];			///< energy deposition [keV] in inner wirechamber window
  Double_t edepDeadMWPC[2];		///< energy deposition [keV] in MWPC dead volume
  Double_t edepKevlar[2];		///< energy deposition [keV] in kevlar strings
  Double_t edepWires[2];			///< energy deposition [keV] in wire planes
  Double_t edepDeadScint[2];		///< energy deposition [keV] in dead scintillator
};
// Variables which come after reverse calibration
struct simCalData {
  Double_t Erecon;                       ///< reconstructed energy 
  Double_t physicsWeight;               ///< Physics weight for each event
  Float_t mwpcEnergyE, mwpcEnergyW;	///< calibrated wirechamber energy deposition on each side
  Int_t type;                   ///< Event type after reverse calibration
  Int_t side;                   ///< Event side after reverse calibration
  Int_t PID;                    ///< Event PID after reverse calibration
  
};

/// plugin with extra plots for simulated runs
class SimTreePlugin: public OctetAnalyzerPlugin {
public:
	/// constructor
	SimTreePlugin(OctetAnalyzer* OA);
        /// destructor
        ~SimTreePlugin() {SimOut->Print(); SimFile->Write();}
        /// make simulation tree
        virtual void makeSimTree(int rn);
	/// fill from scan data point
	//void fillSimTree(ProcessedDataScanner& PDS, double weight);
        virtual void fillCoreHists(ProcessedDataScanner& PDS, double weight);
	/// output plot generation
	//virtual void makePlots();
	

  //void printSimTree() { SimOut->Print(); }
  //Tree to hold simulated data which looks like real analyzed data
  TFile *SimFile;
  TTree *SimOut;
  TBranch *SimDataBranch, *RevCalSimDataBranch; 

  simData SD;
  simCalData SCD;

  //std::ofstream ofile;
};
#endif
