#ifndef PROCESSEDDATASCANNER_HH
#define PROCESSEDDATASCANNER_HH 1

#include "Types.hh"
#include "Enums.hh"

#include "RunSetScanner.hh"
#include "EnergyCalibrator.hh"
#include "EventClassifier.hh"

#include "QFile.hh"
#include "PMTGenerator.hh"
#include "CalDBSQL.hh"
#include "TagCounter.hh"

#include <map>
#include "SMExcept.hh"

/// Generic class for processed data TChains
class ProcessedDataScanner: public RunSetScanner, public EventClassifier {
public:
	/// constructor
	ProcessedDataScanner(const std::string& treeName, bool withCalibrators = false);
	
	/// radius squared of event
	virtual float radius2(Side s) const;
	/// radius of event
	virtual float radius(Side s) const { return sqrt(radius2(s)); }
	/// MWPC charge cloud volume
	virtual float mwpcCharge(Side s) const { return wires[s][X_DIRECTION].ccloud_size() + wires[s][Y_DIRECTION].ccloud_size(); }
	/// MWPC cathode maximum sum
	virtual float mwpcCathMaxSum(Side s) const { return wires[s][X_DIRECTION].maxValue + wires[s][Y_DIRECTION].maxValue; }
	/// event energy
	virtual float getEnergy() const { return scints[EAST].energy.x + scints[WEST].energy.x; }
	/// get event true (reconstructed) energy
	virtual float getErecon() const;
	/// re-calibrate tube energy of currently loaded event
	virtual void recalibrateEnergy();
	/// whether event passes fiducia/position cut on side
	virtual bool passesPositionCut(Side s);
	/// whether event was simulated as triggering the given side
	virtual bool Sis00_2fold(Side s) { return SIS00 & (1<<s); }
	/// get trigger probability for individual PMT
	virtual float probTrig(Side s, unsigned int t);
	/// get info about current event
	virtual Stringmap evtInfo();
	
	/// Type I initial hit side determination --- not available here
	virtual Side getFirstScint() const { smassert(false); return BOTH; }
	/// Type II/III separation probability
	virtual float getProbIII() const { return WirechamberCalibrator::sep23Prob(fSide,getEnergy(),mwpcEnergy[fSide]); }
	
	static bool redoPositions;		//< whether to re-calibrate positions

	ScintEvent scints[BOTH];	//< readout point for scintillator data
	Float_t led_pd[BOTH];		//< readout point for reference photodiode
	wireHit wires[BOTH][2];		//< readout point for wirechamber data [side][plane]
	float cathodes[BOTH][2][kMaxCathodes];	//< readout point for pedestal-subtracted cathode values, [side][plane][cathode]
	MWPCevent mwpcs[BOTH];		//< readout point for mwpc data (anode & cathode sum)
	Float_t mwpcEnergy[BOTH];	//< calibrated wirechamber energy deposition on each side
	BlindTime runClock;			//< time of current event since run start
		
	double physicsWeight;		//< event spectrum re-weighting factor
	
	AnalysisChoice anChoice;	//< which analysis choice to use in identifying event types
	float fiducialRadius;		//< radius for position cut
};

#endif
