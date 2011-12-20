#ifndef WIRECHAMBER_HH
#define WIRECHAMBER_HH 1

#include <stdlib.h>
#include "Subsystem.hh"
#include "Trigger.hh"
#include "WirechamberReconstruction.hh"
#include <vector>
#include <TLine.h>

/// Class for a wirechamber plane detector
class Wirechamber: public Subsystem {
public:
	/// constructor
	Wirechamber(RunManager* T, Trigger* tg, Side S, AxisDirection d);
	/// destructor
	virtual ~Wirechamber() { if(hits) delete[](hits); }
	
	/// Brad Plaster's reconstruction algorithm
	Float_t bradPosition(UInt_t e);
	/// get wirechamber acceptance lower bound
	Float_t getLowerEdge() const { return lowerEdge; }
	/// get wirechamber acceptance lower bound
	Float_t getUpperEdge() const { return upperEdge; }

	/// whether the wirechamber reconstructed position is reasonable
	bool passesPositionCut(UInt_t e) const { return getLowerEdge() < hits[e].center && hits[e].center < getUpperEdge(); }
	
	/// draw wires onto 2D plot
	void drawWires(Float_t low, Float_t high, Int_t c = 7) const {
		for(unsigned int i=0; i<nWires; i++) {
			TLine* l;
			if(myDirection==X_DIRECTION)
				l = new TLine(wirePositions[i],low,wirePositions[i],high);
			else
				l = new TLine(low,wirePositions[i],high,wirePositions[i]);
			l->SetLineColor(c);
			l->Draw();
		}
	}
		
	wireHit anEvent;	//< output tree readout point
	
	/// add relevant branches to output tree
	virtual void addOutBranches(TTree* T) { 
		T->Branch(name.c_str(),&anEvent,
				  "center/F:width/F:maxValue/F:cathodeSum/F:maxWire/I:nClipped/I:multiplicity/I:err/I:gausCenter/F:avgCenter/F:parabCenter/F"); 
	}
	
	/// fill event in output tree
	virtual void fillEvent(UInt_t e) { anEvent = hits[e]; }
	
	
	
	Float_t csumPedestal;	//< cathode sum pedestal
	Float_t csumPedwidth;	//< cathode sum width
	Trigger* TG;			//< event triggers
	wireHit* hits;			//< reconstructed info about each hit
	std::vector<float*> cathdat;	//< shortcut to cathode data
	AxisDirection myDirection;		//< whether the wires track in the x (y) direction
	Float_t cathMaxCut;		//< cathode maximum trigger cut
	Float_t lowerEdge;		//< minimum position tracked by detector
	Float_t upperEdge;		//< maximum position tracked by detector
	
	/// generate relevant histograms
	void genHistograms();

	unsigned int nWires;	//< number of wires in wirechamber
	Float_t wireSpacing;	//< spacing between wires
	std::vector<float> wirePositions;	//< positions of each wire
	int* adcNums;			//< ADC numbers of each wire
	float cathNorm[16];		//< cathode gain normalization
	
protected:
	/// run-specific actions
	void specialize();
};

#endif

