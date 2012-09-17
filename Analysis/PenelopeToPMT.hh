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
	/// calculate physics weighting factor (undo pre-baked asymmetry)
	void calcReweight();
	
	float fEprim;				//< float version for primary energy
	float fEdep[2];				//< float version for scintillator energy
	float fEquench[2];			//< float version for quenched energy
	float fEW[2];				//< float version of wirechamber energy
	float fMWPCpos[2][2];		//< float version of MWPC position
	float fEWd[2];				//< float version of MWPC dead gas energy
	float fEMWPCDead[2][2];		//< float version of MWPC deead gas energy for each [side][front/back]
	float fEMWPCWires[2][2];	//< float version of wire plane energy 
	
	float fedepFoils[2];
	float fedepWinOut[2];
	float fedepWinIn[2];
	float fedepDeadMWPC[2];                 
	float fedepKevlar[2];           
	float fedepWires[2];
	float fedepDeadScint[2];
	
	float fcosThetaInFoils[2];
	float fcosThetaInWinOut[2];
	float fcosThetaInWinIn[2];
	float fcosThetaInScint[2];
	
	float fcosThetaOutFoils[2];             
	float fcosThetaOutWinOut[2];     
	float fcosThetaOutWinIn[2];
	float fcosThetaOutScint[2];
	
	float fPrimPos[3];			//< float version of primary position
	float fTime[2];				//< float version of time
	float fCostheta;			//< float version of cos theta
	
protected:
	virtual void setReadpoints();
};

#endif
