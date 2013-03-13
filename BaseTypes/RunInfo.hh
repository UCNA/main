#ifndef RUNINFO_HH
#define RUNINFO_HH 1

#include <map>
#include <vector>
#include <string>
#include <Rtypes.h>

#include "QFile.hh"
#include "Enums.hh"

/// information about a particular UCNA data run
class RunInfo: public StringmapProvider {
public:	
	
	/// constructor
	RunInfo(RunNum r = 0): StringmapProvider(), runNum(r), slowDaq(0), groupName("0 Unknown"),
	type(UNKNOWN_RUN), octet(0), triad(TRIAD_NONE), gvState(GV_OTHER), afpState(AFP_OTHER), runGeometry(GEOMETRY_OTHER), runQuality(DQ_UNKNOWN),
	roleName("Unknown"), startTime(0), liveTime(0), wallTime(0), scsField(0), comments("") { }
	
	/// print summary info
	void display() const {
		printf("Run=%i AFP=%i GV=%i Type=%i Oct=%i Triad=%i\n",runNum,afpState,gvState,type,octet,triad);
	}
	
	
	RunNum runNum;				//< run number
	Int_t slowDaq;				//< SlowDAQ run number
	std::string groupName;		//< name for this run's group
	
	RunType type;			//< what type of run this is
	OctetType octet;			//< what octet stage this run is
	TriadType triad;			//< what type of triad this belongs to
	GVState gvState;			//< state of gate valve during run
	AFPState afpState;			//< whether spin-flipper is on
	RunGeometry runGeometry;	//< geometry configuration for this run
	DataQuality runQuality;		//< quality of data for this run
	std::string roleName;		//< name of this run's role
	Float_t startTime;			//< absolute start time of run
	BlindTime liveTime;	//< fiducial time for run
	BlindTime wallTime;	//< wall clock time for run
	
	Float_t scsField;			//< field in SCS spectrometer
	std::string comments;		//< run log comments about this run
	
protected:
	
	/// properties for output Stringmap
	virtual Stringmap getProperties() const {
		Stringmap ri;
		ri.insert("Run_Number",runNum);
		ri.insert("Role_Name",roleName);
		ri.insert("SlowDaq",slowDaq);
		ri.insert("SCS_Field",scsField);
		ri.insert("AFP_State",int(afpState));
		ri.insert("GV_State",int(gvState));
		ri.insert("Octet",int(octet));
		ri.insert("Triad",int(triad));
		ri.insert("liveTime",liveTime[BOTH]);
		ri.insert("wallTime",wallTime[BOTH]);
		return ri;
	}
};

#endif
