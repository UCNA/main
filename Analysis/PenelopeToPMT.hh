#ifndef PENELOPETOPMT_HH
#define PENELOPETOPMT_HH 1

#include "Sim2PMT.hh"

/// converts Robby's Penelope data to PMT spectra
class PenelopeToPMT: public Sim2PMT {
public:
	/// constructor
	PenelopeToPMT(): Sim2PMT("h34") { }	
	/// unit conversions
	virtual void doUnits();
	
	float fEprim;			//< float version for primary energy
	float fEdep[2];			//< float version for scintillator energy
	float fEW[2];			//< float version of wirechamber energy
	float fMWPCpos[2][2];	//< float version of MWPC position
	float fPrimPos[3];		//< float version of primary position
	float fTime[2];			//< float version of time
	float fCostheta;		//< float version of cos theta
	
protected:
	virtual void setReadpoints();
};

#endif
