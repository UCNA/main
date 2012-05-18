#include "AnalyzerBase.hh"
#include "Enums.hh"

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
	Double_t EdepSD[19];			//< array for energy deposition in all SDs
	Double_t EdepAll;				//< total edep in all SDs
	Double_t hitTime[2];			//< timing info for hits on each side
	Double_t trapMonTime[2];		//< timing for trap monitor hits
	Double_t thetaIn[2];			//< particle angle entering wirechamber
	Double_t thetaOut[2];			//< particle angle exiting wirechamber
	Double_t thetaInDF[2];			//< particle angle entering decay trap foil
	Double_t thetaOutDF[2];			//< particle angle exiting decay trap foil
	Double_t kEIn[2];				//< particle energy entering scintillator
	Double_t kEInTrapMon[2];		//< particle energy entering decay trap monitor volume
	Double_t kEOut[2];				//< particle energy leaving scintillator
	
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
