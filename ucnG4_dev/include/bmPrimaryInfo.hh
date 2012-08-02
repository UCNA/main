#ifndef __bmPrimaryInfo_hh__
#define __bmPrimaryInfo_hh__

#include "Rtypes.h"
#include "TObject.h"
#include "TObjString.h"

/// ROOT-friendly info on a primary particle
/// we will save a list of these for each event
class bmPrimaryInfo: public TObject {
public:
	/// constructor
	bmPrimaryInfo() {}
	
	Float_t vertex[3];	//< primary vertex position [m]
	Float_t p[3];		//< primary momentum [keV]
	
	Float_t KE;			//< primary KE
	Float_t weight;		//< event generator weight
	Long_t  seed;		//< primary random seed
	
	ClassDef(bmPrimaryInfo, 1);
};

#endif
