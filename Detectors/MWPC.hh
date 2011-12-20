#ifndef MWPC_HH
#define MWPC_HH 1

#include "Subsystem.hh"
#include "Wirechamber.hh"
#include "BetaScint.hh"
#include <vector>

/// Class for an MWPC wire positioning chamber (2 wire planes + anode)
class MWPC: public Subsystem {
public:
	/// constructor from component detectors
	MWPC(RunManager* T, Side s, Wirechamber* x, Wirechamber* y, BetaScint* bs);
	
	Wirechamber* xPlane;	//< wire plane for x positioning
	Wirechamber* yPlane;	//< wire plane for y positioning
	std::vector<MWPCevent> events;		//< event-by-event information
	
	/// generate relevant histograms
	void genHistograms();
	
	/// get list of points falling in specified ellipse
	std::vector<unsigned int> cutEllipse(Source p, Float_t nsigma = 1.0, bool* selected = NULL) const;
	
	/// whether an event passes the maximum wire value cut
	bool passesMaxwireCut(UInt_t e) const { return xPlane->triggers[e] && yPlane->triggers[e]; }
	/// whether an event passes the cathode sum cut
	bool passesCathodesumCut(UInt_t e) const { return events[e].cathodeSum > cathodeCut; }
	/// whether an event passes the anode cut
	bool passesAnodeCut(UInt_t e) const { return events[e].anode > anodeCut; }
	/// whether an event passes the positioning cut
	bool passesPositionCut(UInt_t e) const { return xPlane->passesPositionCut(e) && yPlane->passesPositionCut(e); }
	/// event radius
	float radius(UInt_t e) const { return sqrt(xPlane->hits[e].center*xPlane->hits[e].center + yPlane->hits[e].center*yPlane->hits[e].center); }
	
	MWPCevent anEvent; //< event data readout point for output tree
	
	/// add relevant branches to output tree
	virtual void addOutBranches(TTree* T) { T->Branch(name.c_str(),&anEvent,"cathodeSum/F:anode/F:sID/I:pos[2]/F:errflags[2]/I"); }
	
	/// fill event data for output tree
	virtual void fillEvent(UInt_t e) { anEvent = events[e]; }
	
	/// draw wires onto 2D plot
	void drawWires(Int_t c = 46) const {
		xPlane->drawWires(yPlane->getLowerEdge(), yPlane->getUpperEdge(), c);
		yPlane->drawWires(xPlane->getLowerEdge(), xPlane->getUpperEdge(), c);
	}
	
	/// blank histogram of correct dimensions for wirechamber position plots, registered with run
	TH2F* blankHisto(std::string name, std::string title, unsigned int ngrid = 100) {
		return registeredTH2F(name,title,ngrid,xPlane->getLowerEdge(),xPlane->getUpperEdge(),ngrid,yPlane->getLowerEdge(),yPlane->getUpperEdge());
	}
		
	BetaScint* BS;		//< Beta scintillator behind this wirechamber
	
private:
	
	/// run-specific configuration
	void specialize();
	float anodeCut;		//< anode cut level
	float cathodeCut;	//< cathode cut level
};

#endif

