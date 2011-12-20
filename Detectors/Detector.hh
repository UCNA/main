/// \file Detector.hh \brief The entire UCNA experiment apparatus
#ifndef DETECTOR_HH
/// make sure this file is only included once
#define DETECTOR_HH 1

#include "Subsystem.hh"
#include "Wirechamber.hh"
#include "EventParser.hh"
#include "MWPC.hh"
#include "BetaScint.hh"
#include "MuonVeto.hh"
#include "GraphicsUtils.hh"
#include <TProfile2D.h>
#include <TGraph2D.h>
#include <TF1.h>
#include <TF2.h>
#include <TGraphErrors.h>


/// Complete UCNA Experiment detector, combining all subsystems
class Detector: public EventParser, public Subsystem {
public:
	/// constructor
	Detector(RunManager* T, Trigger* tg, MWPC* mwpce, MWPC* mwpcw, BetaScint* bse, BetaScint* bsw, MuonVeto* mv);
	/// destructor
	virtual ~Detector() { if(events) delete(events); }
	
	MWPC* MWPCs[2];			//< E and W MWPCs
	BetaScint* Scints[2];	//< E and W Beta Scintillators
	MuonVeto* muVeto;		//< muon veto
	Trigger* TG;			//< event triggers
	
	/// generate plots
	void genHistograms();
	/// draw cuts for fiducial volume, sources
	void draw2Dcuts(Side s);
	
	/// whether an event is a valid trigger for this subsystem
	virtual bool triggered(UInt_t e) const { 
		if(thisRun->RI.type == LED_RUN)
			return isGMS(e);
		return whichScintSide(e) != NONE;
	}
	
	void addOutBranches(TTree* T);
	void fillEvent(UInt_t e);

protected:
	/// identify each event
	void idEvents();
	/// analyze sources and energy calibration spectra
	void calibrationSpectra();
	
	BetaHit anEvent;
	BetaHit* events;
	std::vector<int> idflags[2];
};

#endif

