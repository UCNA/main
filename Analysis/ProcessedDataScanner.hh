#ifndef PROCESSEDDATASCANNER_HH
#define PROCESSEDDATASCANNER_HH 1

#include "Types.hh"
#include "Enums.hh"

#include "RunSetScanner.hh"
#include "EnergyCalibrator.hh"

#include "QFile.hh"
#include "PMTGenerator.hh"
#include "CalDBSQL.hh"
#include "TagCounter.hh"

#include <map>

/// Generic class for processed data TChains
class ProcessedDataScanner: public RunSetScanner {
public:
	/// constructor
	ProcessedDataScanner(const std::string& treeName, bool withCalibrators = false);
	
	/// radius squared of event
	virtual float radius2(Side s) const;
	/// radius of event
	virtual float radius(Side s) const { return sqrt(radius2(s)); }
	/// event energy
	virtual float getEnergy() const { return scints[EAST].energy.x + scints[WEST].energy.x; }
	/// get event true (reconstructed) energy
	virtual float getEtrue();
	/// re-calibrate tube energy of currently loaded event
	virtual void recalibrateEnergy();
	/// whether event passes fiducia/position cut on side
	virtual bool passesPositionCut(Side s);
	/// get info about current event
	virtual Stringmap evtInfo();
	
	bool redoPositions;			//< whether to re-calibrate positions

	ScintEvent scints[2];		//< readout point for scintillator data
	Float_t led_pd[2];			//< readout point for reference photodiode
	wireHit wires[2][2];		//< readout point for wirechamber data [side][plane]
	float cathodes[2][2][kMaxCathodes];	//< readout point for pedestal-subtracted cathode values, [side][plane][cathode]
	MWPCevent mwpcs[2];			//< readout point for mwpc data (anode & cathode sum)
	Float_t mwpcEnergy[2];		//< calibrated wirechamber energy deposition on each side
	BlindTime runClock;			//< time of current event since run start
	
	PID fPID;					//< analysis particle ID
	EventType fType;			//< analysis event type
	Float_t fProbIII;			//< probability of Type III backscatter
	Side fSide;					//< analysis event side
	UInt_t fTaggedBack[2];		//< whether event was tagged by the muon backing veto on each side
	double physicsWeight;		//< event spectrum re-weighting factor
	
	AnalysisChoice anChoice;	//< which analysis choice to use in identifying event types
	float fiducialRadius;		//< radius for position cut
		
protected:
	
	/// generate event classification flags
	virtual void calcEventFlags() {}
};

#endif
