#include "PenelopeToPMT.hh"
#include "Types.hh"
#include "strutils.hh"
#include "BetaSpectrum.hh"

void PenelopeToPMT::setReadpoints() {
	
	for(Side s = EAST; s <= WEST; ++s) {
		SetBranchAddress(sideSubst("Ep%c",s,true),&fEdep[s]);
		SetBranchAddress(sideSubst("nph%c",s,true),&fEquench[s]);
		SetBranchAddress(sideSubst("Ep%cd",s,true),&fedepDeadScint[s]);
		
		SetBranchAddress(sideSubst("Eg%ca",s,true),&fEW[s]);
		SetBranchAddress(sideSubst("Eg%cdf",s,true),&fEMWPCDead[s][0]);
		SetBranchAddress(sideSubst("Eg%cdb",s,true),&fEMWPCDead[s][1]);
		
		SetBranchAddress(sideSubst("E%cc1",s,true),&fEMWPCWires[s][0]);
		SetBranchAddress(sideSubst("E%cc2",s,true),&fEMWPCWires[s][1]);
		SetBranchAddress(sideSubst("E%can",s,true),&fEMWPCWires[s][2]);
		
		SetBranchAddress(sideSubst("Em%c1",s,true),&fedepWinIn[s]);
		SetBranchAddress(sideSubst("Em%c2",s,true),&fedepWinOut[s]);
		
		for(AxisDirection d=X_DIRECTION; d<=Y_DIRECTION; ++d)
			SetBranchAddress(sideSubst("%cpos",s,true)+(d==X_DIRECTION?"x":"y"),&fMWPCpos[s][d]);
		
		SetBranchAddress(sideSubst("Efl%c",s,true),&fedepFoils[s]);
		SetBranchAddress(sideSubst("t%c",s,true),&fTime[s]);
	}
	
	for(AxisDirection d=X_DIRECTION; d<=Z_DIRECTION; ++d)
		SetBranchAddress(d==X_DIRECTION?"X":d==Y_DIRECTION?"Y":"Z",&fPrimPos[d]);
	
	
	SetBranchAddress("W1",&fcosThetaInFoils[EAST]);
	SetBranchAddress("W3",&fcosThetaInWinIn[EAST]);
	SetBranchAddress("W4",&fcosThetaInWinOut[EAST]);
	SetBranchAddress("W6",&fcosThetaInScint[EAST]);
	
	SetBranchAddress("W7",&fcosThetaOutFoils[EAST]);
	SetBranchAddress("W9",&fcosThetaOutWinIn[EAST]);
	SetBranchAddress("W10",&fcosThetaOutWinOut[EAST]);
	SetBranchAddress("W12",&fcosThetaOutScint[EAST]);
	
	SetBranchAddress("E",&fEprim);
	SetBranchAddress("W",&fCostheta);
}

void PenelopeToPMT::calcReweight() {
	physicsWeight = basePhysWeight;
	//Sim2PMT::calcReweight();
	// remove baked-in asymmetry and spectrum corrections
	//if(afp != AFP_OTHER)
	//	physicsWeight /= (1.0+correctedAsymmetry(ePrim,-costheta))*spectrumCorrectionFactor(ePrim);
}

void PenelopeToPMT::doUnits() {
	const double posConversion = 10.0;
	for(AxisDirection d=X_DIRECTION; d<=Z_DIRECTION; ++d) {
		primPos[d] = fPrimPos[d]*posConversion;
		for(Side s = EAST; s <= WEST; ++s)
			mwpcPos[s][d] = fMWPCpos[s][d]*posConversion*sqrt(0.6);
	}
	for(Side s = EAST; s <= WEST; ++s) {
		eDep[s] = fEdep[s]*0.001;
		eQ[s] = fEquench[s];
		eW[s] = (fEW[s]+0.5*(fEMWPCDead[s][0]+fEMWPCDead[s][1]))*0.001;
		time[s] = fTime[s]*1e-9; // convert from nanoseconds to seconds
		for(AxisDirection d=X_DIRECTION; d<=Y_DIRECTION; ++d)
			scintPos[s][d] = mwpcPos[s][d];	// fake scintillator pos from MWPC
	}
	ePrim = fEprim*0.001;
	costheta = -fCostheta;
	for(Side s = EAST; s <= WEST; ++s) {
		edepWires[s] = (fEMWPCWires[s][0]+fEMWPCWires[s][1]+fEMWPCWires[s][2])*0.001;
		edepDeadMWPC[s] = fEMWPCDead[s][0]*0.001 + fEMWPCDead[s][1]*0.001;
		edepFoils[s]  = fedepFoils[s]*0.001;
		edepWinOut[s] = fedepWinOut[s]*0.001;
		edepWinIn[s]  = fedepWinIn[s]*0.001;
		edepDeadScint[s] = fedepDeadScint[s]*0.001;
		
		cosThetaInFoils[s] = fcosThetaInFoils[s];
		cosThetaInWinOut[s] = fcosThetaInWinOut[s];
		cosThetaInWinIn[s] = fcosThetaInWinIn[s];
		cosThetaInScint[s] = fcosThetaInScint[s];
		
		cosThetaOutFoils[s] = fcosThetaOutFoils[s];
		cosThetaOutWinOut[s] = fcosThetaOutWinOut[s];
		cosThetaOutWinIn[s] = fcosThetaOutWinIn[s];
		cosThetaOutScint[s] = fcosThetaOutScint[s];
	}
}
