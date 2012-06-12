#include "PenelopeToPMT.hh"

void PenelopeToPMT::setReadpoints() {
	Tch->SetBranchAddress("Epe",&fEdep[EAST]);
	Tch->SetBranchAddress("Epw",&fEdep[WEST]);
	Tch->SetBranchAddress("Egea",&fEW[EAST]);
	Tch->SetBranchAddress("Egwa",&fEW[WEST]);
	
	Tch->SetBranchAddress("eposx",&fMWPCpos[EAST][X_DIRECTION]);
	Tch->SetBranchAddress("eposy",&fMWPCpos[EAST][Y_DIRECTION]);
	Tch->SetBranchAddress("wposx",&fMWPCpos[WEST][X_DIRECTION]);
	Tch->SetBranchAddress("wposy",&fMWPCpos[WEST][Y_DIRECTION]);
	
	Tch->SetBranchAddress("te",&fTime[EAST]);
	Tch->SetBranchAddress("tw",&fTime[WEST]);
	
	Tch->SetBranchAddress("X",&fPrimPos[X_DIRECTION]);
	Tch->SetBranchAddress("Y",&fPrimPos[Y_DIRECTION]);
	Tch->SetBranchAddress("Z",&fPrimPos[Z_DIRECTION]);
	
	Tch->SetBranchAddress("E",&fEprim);
	Tch->SetBranchAddress("W",&fCostheta);
}

void PenelopeToPMT::doUnits() {
	const double wcPosConversion = 10.0;
	for(unsigned int i=0; i<3; i++)
		primPos[i] = fPrimPos[i]*wcPosConversion;
	for(Side s = EAST; s <= WEST; ++s) {
		eQ[s] = fEdep[s]*0.001;	// really Edep and not Equenched
		eW[s] = fEW[s]*0.001;
		time[s] = fTime[s];
		for(AxisDirection d=X_DIRECTION; d<=Y_DIRECTION; ++d) {			
			scintPos[s][d] = primPos[d];	// fake scintillator pos from primary
			mwpcPos[s][d]  = primPos[d]; 	// same position in scintillator
			wires[s][d].center = mwpcPos[s][d];
		}
	}
	ePrim = fEprim*0.001;
	costheta = fCostheta;
}
