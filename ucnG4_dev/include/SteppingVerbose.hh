#ifndef SteppingVerbose_h
#define SteppingVerbose_h

#include <G4SteppingVerbose.hh>

/// class for printing out information at every step
class SteppingVerbose : public G4SteppingVerbose {
public:
	/// constructor
	SteppingVerbose() {}
	/// display info for step
	void StepInfo();
	/// display info at start of track
	void TrackingStarted();
};

#endif
