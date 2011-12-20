#ifndef TRIGGER_HH
#define TRIGGER_HH 1


#include "Subsystem.hh"
#include "Enums.hh"
#include "TH1F.h"
#include <stdlib.h>

/// blip in data (useless segment to be discarded)
struct Blip {
	unsigned int start;
	unsigned int end;
	float time;
	int counts() const { return int(end)-int(start); }
	void display() const { printf("Blip [%i-%i] %i counts, t=%fs\n",start,end,counts(),time); }
};

/// UCNA data triggers
class Trigger: public Subsystem {
public:
	/// constructor
	Trigger(RunManager* T);
	
	/// get sis00 trigger flags for an event
	unsigned int sis00(unsigned int e) const { return (int)getData("Sis00")[e]; }
	/// if this event is quality cut (during beam, no beam, etc.)
	bool isCrud(UInt_t e) const { return !triggers[e]; }	
	
	/// if this is a UCN monitor event
	bool isUCNMon(UInt_t e) const { return sis00(e) & ((1<<2) | (1<<8) | (1<<9) | (1<<10) | (1<<11)); }
	/// if this is a UCN monitor event
	bool isMon1(UInt_t e) const { return sis00(e) & (1<<8); }
	/// if this is a UCN monitor event
	bool isMon2(UInt_t e) const { return sis00(e) & (1<<9); }
	/// if this is a UCN monitor event
	bool isMon3(UInt_t e) const { return sis00(e) & (1<<10); }
	/// if this is a UCN monitor event
	bool isMon4(UInt_t e) const { return sis00(e) & (1<<11); }
	
	/// if this is a GMS event
	bool gmsTrigger(UInt_t e) const { return sis00(e) & (1<<7); }
	/// which side this GMS event is on
	Side gmsSide(UInt_t e) const {
		if(sis00(e) & (1<<7)) return BOTH;
		return NONE; 
	}
	
	/// if this is a PMT Pulser event
	bool pulserTrigger(UInt_t e) const { return nFiring(e,EAST)+nFiring(e,WEST) == 1; }
	
	/// if this is a normal decay event with no extra triggers
	bool isDataEvent(UInt_t e) const { return sis00(e) <= 3; }
	/// if this is during beam noise
	bool isBeamnoise(UInt_t e) const { return getData("beamclock")[e] < 0.08; }
	/// if this is during a beam outage
	bool isBeamout(UInt_t e) const { return getData("beamclock")[e] > 20.0; }
	
	/// if this is a reference Co60 source decay on the given side
	bool gmsCo(UInt_t e, Side s) const { return (s == EAST && sis00(e) == 32) || (s == WEST && sis00(e) == 64); }
	/// if this is a GMS LED firing
	bool gmsLED(UInt_t e) const { return sis00(e) & (1<<7) || ((sis00(e) & (1<<5)) && (sis00(e) & (1<<6))); }
	/// if 2 of 4 Beta Scintillator PMTs fired on the given side
	bool beta2of4(UInt_t e, Side s) const { return getData("2of4_tdc",s)[e] > 5.0; }	
	/// which side the scintillator trigger was on
	Side beta2of4Side(UInt_t e) const { return sideCombo(beta2of4(e,EAST), beta2of4(e,WEST)); }
	/// whether a particular PMT triggered
	bool betaTrigger(UInt_t e, Side s, unsigned int t) const { return events[e].trigflags & 1<<(nBetaTubes*s+t); }
	/// number of PMTs triggering for event
	int nFiring(UInt_t e, Side s) const {
		int nf = 0;
		for(unsigned int t=0; t<nBetaTubes; t++)
			if(betaTrigger(e,s,t))
				nf++;
		return nf;
	}
	
	/// whether the given side was a primary scintillator event
	bool primaryScint(UInt_t e, Side s) const {
		assert(s<=WEST);
		return scintTDCcuts[s] < getData("2of4_tdc",s)[e] && getData("2of4_tdc",s)[e] < 3990; 
	}
	/// whether this was the second scintillator to fire during an event
	bool secondaryScint(UInt_t e, Side s) const {
		assert(s<=WEST);
		return 5 < getData("2of4_tdc",s)[e] && getData("2of4_tdc",s)[e] < scintTDCcuts[s]; 
	}
	/// whether this was a non-firing scintillator during an event
	bool nonfiringScint(UInt_t e, Side s) const {
		assert(s<=WEST);
		return getData("2of4_tdc",s)[e] < 5; 
	}
		
	/// time of this event in s since beginning of the run (with blinding clocks)
	Float_t eventTime(UInt_t e) const { return getData("runclock")[e]-getData("runclock")[0]; }
	
	/// time of useful data in this run in s, with blinding clock
	Float_t runTime(Side s = BOTH) const { return runtimes[s]; }
	
	/// total time of this run
	Float_t totalTime() const { return getData("runclock")[nEvents-1] - getData("runclock")[0]; }
	
	
	/// set up output tree
	void addOutBranches(TTree* T);
	/// fill output tree
	void fillEvent(UInt_t e);
	
	/// Count missed triggers due to dead time
	void checkDeadtime();
	/// Make rate histograms
	void makeHistos();
	/// output items for RunData file
	virtual Stringmap finalWords() { 
		forParent.insert("runTime",runTime());
		forParent.insert("runTime_E",runTime(EAST));
		forParent.insert("runTime_W",runTime(WEST));
		forParent.insert("nEvents",itos(nEvents));
		forParent.insert("totalTime",totalTime());
		return forParent;
	}
	
	std::vector<TrigInfo> events;
	
private:
	
	/// run-specific configuration
	void specialize();
	
	/// fix time wrap-around problem, convert times to seconds
	void fixTimes();
	/// calculate actual running time after cuts
	void calcRuntime();
	
	/// characterize beam pulse events and make relevant cuts
	void examineBeamStructure();
	/// characterize rapid event clusters and make relevant cuts
	void examineClusterTiming();
	
	/// generate trigger flags
	TrigFlags getTrigFlags(unsigned int e) const {
		unsigned int f = 0;
		for(Side s = EAST; s <= WEST; s = nextSide(s))
			for(unsigned int t=0; t<nBetaTubes; t++)
				if(getData("tube_tdc",nBetaTubes*s+t)[e] > 5)
					f |= 1<<(nBetaTubes*s+t);
		return TrigFlags(f);
	}
	
	/// make blips from a list of trigger events
	std::vector<Blip> makeBlips(const bool* trg) const;
	
	std::vector<Blip> blips;	//< data quality blips
	TrigInfo anEvent;		//< output tree filling point
	float scintTDCcuts[2];	//< timing cut positions
	float runtimes[4];		//< (blinded) useful running times after cuts
};

#endif


