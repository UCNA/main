///This is a simple class that will handle the MC-ROOT interface.
///The main point is to create a global interface object which 
///renders the flexibility for user "rooting" at various level 
///in the MC, e.g. detector construction, tracking, etc.
///Created based on IcaG4 from http://root.cern.ch/root/HowtoMC.html
// Jianglai 05/26/2006 --- Modified MPM 9/2011
#ifndef bmAnalysisManager_h
#define bmAnalysisManager_h 1

class G4VPhysicalVolume;
class G4Event;
class G4Run;
class G4Track;
class G4Step;
class G4PrimaryVertex;
class G4PrimaryParticle;

#include "globals.hh"
#include <TH1.h>
#include <TTree.h>
#include <TFile.h>
#include <TStopwatch.h>
#include <TString.h>

#include "bmEventAction.hh"
#include "bmSteppingAction.hh"
#include "bmTrackerHit.hh"
#include "bmTrackerSD.hh"
#include "bmMCEvent.hh"
#include <vector>

class bmAnalysisManager;

extern bmAnalysisManager *gbmAnalysisManager; // global bmAnalysisManager

/// handles storing ROOT data for each event
class bmAnalysisManager {
	
public:
	/// constructor
	bmAnalysisManager();
	/// destructor
	~bmAnalysisManager() {
		if (gbmAnalysisManager == this)
			gbmAnalysisManager = (bmAnalysisManager *)0;
	}
	/// get global analysis manager
	static bmAnalysisManager* GetAnalysisManager() { return gbmAnalysisManager; }
	
	/// open ROOT output file
	void OpenFile(const G4String filename);
	/// write and close ROOT output file
	void CloseFile();
	/// clear event info
	void Clear(){ mcEvent.ClearEvent(); }
	/// set up output tree
	void CreateTrees();
	
	void StoreHitCollectionIDs();
	
	void SetRunNumber(const Int_t runno) { fRunNumber = runno; }
	Int_t GetRunNumber(void) const {return fRunNumber;}
	
	/// convert primary events data to ROOT form
	void FillPrimaryData(const G4Event* evt_in, const long);
	/// convert tracking data to ROOT form
	void FillTrackerData(const G4Event *evt);
	/// fill output tree
	void FillEventTree();
	
	bmMCEvent mcEvent;		//< saveable event
	bmMCEvent* pMcEvent;	//< pointer to said event for TTree branch setup
	
	/// store a sensitive detector name to record tracking info for
	void SaveSDName(const G4String name){ fSDNames.push_back(name); }
		
private:
	
	TFile *fROOTOutputFile;		//< ROOT output file
	TTree *fEventTree;			//< ROOT output TTree
	Int_t fRunNumber;  			//< MC run number
	UInt_t fSeed;				//< MC random seed
	vector<G4int> detectorIDs;	//< list of SD ID numbers
	vector<G4String> fSDNames;	//< list of SD names corresponding to ID numbers
};

#endif
