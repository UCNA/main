#include "AnalyzerBase.hh"
#include "Enums.hh"

class SiDet_Analyzer: public ucnG4_analyzer {
public:
	/// constructor
	SiDet_Analyzer(const std::string& outfname): ucnG4_analyzer(outfname) {}
	
	Double_t Edep;				//< deposited energy in Si detector
	Double_t Pos[3];			//< average position in Si detector
		
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
	virtual bool saveEvent() { return Edep>0; }
};
