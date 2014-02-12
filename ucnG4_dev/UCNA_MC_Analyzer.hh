#include "AnalyzerBase.hh"
#include "Enums.hh"
#include "WirechamberCalibrator.hh"

/// number of sensitive detector regions
#define N_SD 24
/// number of points to integrate around circle for cathode response
#define N_CHG_CIRC_PTS 7

class UCNA_MC_Analyzer: public ucnG4_analyzer {
public:
	/// constructor
	UCNA_MC_Analyzer(const std::string& outfname);
	
	bool saveAllEvents;			//< whether to save non-energy-depositing events to file
	bool undeadLayer;			//< whether to ignore scintillator dead layer
	bool calcCathCharge;		//< whether to calculate cathode charge distribution
	
	Double_t Edep[BOTH];						//< scintillator deposited energy
	Double_t EdepQ[BOTH];						//< quenched energy in scintillator
	Double_t MWPCPos[BOTH][Z_DIRECTION+1];		//< MWPC deposited energy weighted position
	Double_t ScintPos[BOTH][Z_DIRECTION+1];		//< scintillator deposited energy weighted position
	Double_t MWPCPosSigma[BOTH][Z_DIRECTION+1];	//< MWPC deposited energy weighted position variance
	Double_t ScintPosSigma[BOTH][Z_DIRECTION+1];//< scintillator quenched energy weighted position variance
	
	Double_t fMWPCEnergy[BOTH];		//< MWPC deposited energy
	Double_t EdepSD[N_SD];			//< array for energy deposition in all SDs
	Double_t thetaInSD[N_SD];		//< entrance angle in each sensitive detector
	Double_t thetaOutSD[N_SD];		//< exit angle for each sensitive detector
	Double_t keInSD[N_SD];			//< kinetic energy entering each sensitive detector
	Double_t keOutSD[N_SD];			//< kinetic energy exiting each sensitive detector
	Int_t hitCountSD[N_SD];			//< count of primary tracks in each volume
	Double_t EdepAll;				//< total edep in all SDs
	Double_t hitTime[BOTH];			//< timing info for hits on each side
	Double_t hitTimeSD[N_SD];		//< earliest hit time in each SD
	Double_t trapMonTime[BOTH];		//< timing for trap monitor hits
	
	Double_t cathWirePos[kWiresPerCathode*kMaxCathodes];		//< individual cathode wire positions, [cm] in physical (non-projected) position
	Float_t cathCharge[BOTH][Y_DIRECTION+1][kMaxCathodes];		//< estimated cathode segment charge signal
	
protected:

	Double_t circ_pts[Y_DIRECTION+1][N_CHG_CIRC_PTS];
	
	/// add additional branches to output tree
	virtual void setupOutputTree();
	
	/// reset analysis values for new event
	virtual void resetAnaEvt();
	/// process current track segment
	virtual void processTrack();
	/// final whole-event processing
	virtual void processEvent();
	/// determine whether an event should be saved to output file
	virtual bool saveEvent() { return saveAllEvents || Edep[EAST]+Edep[WEST] > 0; }
};
