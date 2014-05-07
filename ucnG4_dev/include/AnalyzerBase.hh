#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TVector3.h>
#include <Rtypes.h>

#include <vector>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <cassert>

#include "TrackInfo.hh"
#include "MCEvent.hh"
#include "PrimaryInfo.hh"

using namespace std;

/// event primary vertex information
struct vtxInfo {
	Int_t fTrapped;				///< whether this particle is trapped in mag field (lives too long)
	Double_t fCompTime;			///< computer time needed to track this event
	Long_t seed;				///< random seed used for this event
	
	Double_t primKE;			///< primary event kinetic energy
	Double_t primTheta;			///< primary event emission angle
	Double_t primWeight;		///< primary event generator weight
	Double_t primPos[4];		///< primary event position (4th coordinate = radius)
};

/// base class for analyzing UCNA Geant4 MC simulations
class ucnG4_analyzer {
public:
	/// constructor
	ucnG4_analyzer(const std::string& outfname);
	/// destructor
	virtual ~ucnG4_analyzer();
	
	/// analyze a file
	void analyzeFile(const string& fname);
	/// analyze all files in list file
	void analyzeFileList(const string& flist);
	
	vtxInfo vtx;				///< event vertex data write point
	unsigned int pass_number;	///< event analysis pass number
	int pID;					///< track pID
	int detectorID;				///< track detector ID
	int trackID;				///< track segment number
	
protected:
	
	TTree* anaTree;			///< analysis results output tree
	TFile* outf;			///< output file
	MCEvent* myevt;			///< current event being analyzed
	TrackInfo* trackinfo;	///< current track info
	PrimaryInfo* priminfo;	///< current primary info
	
	/// add additional branches to output tree
	virtual void setupOutputTree() {}
	
	/// reset analysis values for new event
	virtual void resetAnaEvt() {}
	/// individual track segment processing
	virtual void processTrack() {}
	/// whole-event processing after scanning individual track segments
	virtual void processEvent() {}
	/// write results to output TTree
	virtual void writeEvent() { anaTree->Fill(); }
};
