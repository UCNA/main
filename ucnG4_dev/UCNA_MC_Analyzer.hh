#include "AnalyzerBase.hh"
#include "Enums.hh"

/// number of sensitive detector regions
#define N_SD 24

class UCNA_MC_Analyzer: public ucnG4_analyzer {
public:
	/// constructor
	UCNA_MC_Analyzer(const std::string& outfname): ucnG4_analyzer(outfname), saveAllEvents(false) {}
	
	bool saveAllEvents;				//< whether to save non-energy-depositing events to file
	
	Double_t Edep[2];				//< scintillator deposited energy
	Double_t EdepQ[2];				//< quenched energy in scintillator
	Double_t MWPCPos[2][3];			//< MWPC deposited energy weighted position
	Double_t ScintPos[2][3];		//< scintillator quenched energy weighted position	
	Double_t MWPCPosSigma[2][3];	//< MWPC deposited energy weighted position variance
	Double_t ScintPosSigma[2][3];	//< scintillator quenched energy weighted position variance
	
	Double_t fMWPCEnergy[2];		//< MWPC deposited energy
	Double_t EdepSD[N_SD];			//< array for energy deposition in all SDs
	Double_t thetaInSD[N_SD];		//< entrance angle in each sensitive detector
	Double_t thetaOutSD[N_SD];		//< exit angle for each sensitive detector
	Double_t keInSD[N_SD];			//< kinetic energy entering each sensitive detector
	Double_t keOutSD[N_SD];			//< kinetic energy exiting each sensitive detector
	Int_t hitCountSD[N_SD];			//< count of primary tracks in each volume
	Double_t EdepAll;				//< total edep in all SDs
	Double_t hitTime[2];			//< timing info for hits on each side
	Double_t hitTimeSD[N_SD];		//< earliest hit time in each SD
	Double_t trapMonTime[2];		//< timing for trap monitor hits
	
protected:
	
	/// add additional branches to output tree
	virtual void setupOutputTree();
	
	/// reset analysis values for new event
	virtual void resetAnaEvt();
	/// process current track segment
	virtual void processTrack();
	/// final whole-event processing
	virtual void processEvent();
	/// determine whether an event should be saved to output file
	virtual bool saveEvent() { return saveAllEvents || Edep[EAST]+Edep[WEST]>0; }
};
