#ifndef __TrackInfo_hh__
#define __TrackInfo_hh__

#include <Rtypes.h>
#include <TObject.h>

using namespace std;

/// information for a section of a particle track within one detector
/// converted to ROOT-friendly form from TrackerHit
/// we will save a list of these for each event
class TrackInfo: public TObject {
public:
	/// constructor
	TrackInfo() {}
		
	Int_t trackID;			///< ID number for this track
	Int_t hcID;				///< sensitive detector ID number for this track
	Double_t hitTime;		///< time of track start [ns]
	Double_t KE;			///< KE at start of track [keV]
	Double_t Edep;			///< accumulated deposited energy along track (not counting secondaries), [keV]
	Double_t EdepQuenched;	///< accumulated "quenched" energy along track [keV]
	Double_t pIn[3];		///< momentum at entry to volume
	Double_t pOut[3];		///< momentum at exit to volume

	Double_t edepPos[3];	///< hit position weighted by deposited energy, [keV*cm]
	Double_t edepPos2[3];	///< hit position^2 weighted by deposited energy, [kev*cm^2]
	
	Double_t inPos[3];		///< entry position to volume [cm]
	Double_t vertexPos[3];	///< track vertex position [cm]
	
	Int_t pID;				///< particle creating this track (PDG code)
	bool isEntering;		///< whether this is initial track entering a volume
	
	ClassDef(TrackInfo, 2);
};

#endif
