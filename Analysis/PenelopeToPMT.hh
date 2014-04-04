#ifndef PENELOPETOPMT_HH
#define PENELOPETOPMT_HH

#include "Sim2PMT.hh"

/// converts Robby's Penelope data to PMT spectra
class PenelopeToPMT: public Sim2PMT {
public:
	/// constructor
	PenelopeToPMT(): Sim2PMT("tree") { }	
	/// unit conversions
	virtual void doUnits();
	/// calculate physics weighting factor (undo pre-baked asymmetry)
	void calcReweight();
	
	double fEprim;				//< float version for primary energy
	double fEdep[2];				//< float version for scintillator energy
	double fEquench[2];			//< float version for quenched energy
	double fEW[2];				//< float version of wirechamber energy
	double fMWPCpos[2][2];		//< float version of MWPC position
	double fSCINTpos[2][2];           //< float version of Scintillator absorption position
	double fEWd[2];				//< float version of MWPC dead gas energy
	double fEMWPCDead[2][2];		//< float version of MWPC deead gas energy for each [side][front/back]
	double fEMWPCWires[2][3];	//< float version of wire plane energy
	
	double fedepFoils[2];
	double fedepWinOut[2];
	double fedepWinIn[2];
	double fedepDeadMWPC[2];                 
	double fedepKevlar[2];           
	double fedepWires[2];
	double fedepDeadScint[2];
	
	double fcosThetaInFoils[2];
	double fcosThetaInWinOut[2];
	double fcosThetaInWinIn[2];
	double fcosThetaInScint[2];
	
	double fcosThetaOutFoils[2];             
	double fcosThetaOutWinOut[2];     
	double fcosThetaOutWinIn[2];
	double fcosThetaOutScint[2];
	
	double fPrimPos[3];			//< float version of primary position
	double fTime[2];				//< float version of time
	double fCostheta;			//< float version of cos theta
	
protected:
	virtual void setReadpoints();
};

#endif
