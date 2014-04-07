#ifndef TrackerHit_h
#define TrackerHit_h

#include "TrackInfo.hh"

#include <G4VHit.hh>
#include <G4THitsCollection.hh>
#include <G4Allocator.hh>
#include <G4ThreeVector.hh>
#include <G4String.hh>


/// accumuates segment-by-segment information for a track in an SD
class TrackerHit : public G4VHit {
public:
	/// constructor
	TrackerHit();
	/// special new()
	inline void* operator new(size_t);
	/// special delete()
	inline void  operator delete(void*);
	
	/// print out info about this hit
	void Print();
	
	/// copy info into ROOT-friendly TrackInfo
	void fillTrackInfo(TrackInfo& h) const;
	
	void SetTrackID(G4int track) { trackID = track; }
	void SetIncidentEnergy (G4double energy) { incidentEnergy = energy; }
	void SetPos (G4ThreeVector xyz) { hitPosition = xyz; }
	void SetHitTime(G4double time) { hitTime = time; }
	void AddEdep(G4double edep, G4ThreeVector xyz) {
		eDepSoFar += edep;
		edepWeightedPosition += xyz*edep;
		for(unsigned int i=0; i<3; i++)
			edepWeightedPosition2[i] += xyz[i]*xyz[i]*edep;
	}
	void AddEdepQuenched(G4double edep) { eDepQuenchedSoFar += edep; }
	void SetIncidentMomentum(G4ThreeVector pin) { incidentMomentum=pin; }
	void SetExitMomentum(G4ThreeVector pout) { exitMomentum=pout; }
	void SetPID(G4int p) { pID=p; }
	void SetProcessName(G4String sproc) { processName = sproc; }
	void SetVolumeName(G4String svol) { volumeName = svol; }
	void SetVertex(G4ThreeVector xyz) { vertex = xyz; };
	void SetCreatorVolumeName(G4String sname) { creatorVolumeName = sname; }
	
	G4int GetTrackID() const { return trackID; };
	G4double GetIncidentEnergy() const { return incidentEnergy; };      
	G4ThreeVector GetPos() const { return hitPosition; };
	G4double GetEdep() const { return eDepSoFar; };
	G4double GetEdepQuenched() const { return eDepQuenchedSoFar; };
	G4ThreeVector GetEdepPos() const { return edepWeightedPosition; }
	G4ThreeVector GetEdepPos2() const { return edepWeightedPosition2; }
	G4double GetHitTime() const { return hitTime; };
	G4ThreeVector GetIncidentMomentum() const {return incidentMomentum;}
	G4ThreeVector GetExitMomentum() const {return exitMomentum;}
	G4int GetPID() const { return pID; }
	G4String GetProcessName() const { return processName; }
	G4String GetVolumeName() const { return volumeName; }
	G4ThreeVector GetVertex() const { return vertex; };
	G4String GetCreatorVolumeName() const { return creatorVolumeName; }
	
	G4double		originEnergy;			//< energy at split from "originating" track for EdepQ tracking
	unsigned int	nSecondaries;			//< number of secondaries produced along track
private:
	
	G4int			trackID;				//< ID number for this track
	G4double		incidentEnergy;			//< incident energy at start of track
	G4double		eDepSoFar;				//< accumulator for energy deposition
	G4double		eDepQuenchedSoFar;		//< accumulator for quenched energy
	G4double		hitTime;				//< entry time into volume
	G4ThreeVector	hitPosition;			//< position where this track entered volume
	G4ThreeVector	edepWeightedPosition;	//< track position weighted by deposited energy
	G4ThreeVector	edepWeightedPosition2;	//< track position^2 weighted by deposited energy
	G4ThreeVector	incidentMomentum;		//< momentum when entering volume
	G4ThreeVector	exitMomentum;			//< momentum when exiting volume
	G4int 			pID;					//< particle creating track (PDG code)
	G4String processName;					//< name of process creating track
	G4String volumeName;					//< name of volume where track is
	G4ThreeVector	vertex;					//< track vertex position
	G4String creatorVolumeName;				//< volume where track was created
};

typedef G4THitsCollection<TrackerHit> TrackerHitsCollection;

extern G4Allocator<TrackerHit> TrackerHitAllocator;

inline void* TrackerHit::operator new(size_t) {
	void *aHit;
	aHit = (void *) TrackerHitAllocator.MallocSingle();
	return aHit;
}

inline void TrackerHit::operator delete(void *aHit) {
	TrackerHitAllocator.FreeSingle((TrackerHit*) aHit);
}

#endif
