#include "G4toPMT.hh"
#include "PathUtils.hh"
#include <cassert>
#include <cmath>

void G4toPMT::setReadpoints() {
	Tch->SetBranchAddress("EdepQ",eQ);
	Tch->SetBranchAddress("Edep",eDep);
	Tch->SetBranchAddress("MWPCEnergy",eW);
	Tch->SetBranchAddress("ScintPos",scintPos);
	Tch->SetBranchAddress("MWPCPos",mwpcPos);
	Tch->SetBranchAddress("time",time);
	Tch->SetBranchAddress("primTheta",&costheta);
	Tch->SetBranchAddress("primKE",&ePrim);
	if(Tch->GetBranch("primPos"))
		Tch->SetBranchAddress("primPos",primPos);
	else
		primPos[0] = primPos[1] = primPos[2] = primPos[3] = 0;
}

void G4toPMT::doUnits() {
	for(unsigned int i=0; i<4; i++)
		primPos[i] *= 1000.0;	// convert m to mm
	// wirechamber position projection plus empirical window-diameter-matching fudge factor
	const double wcPosConversion = 10.0*sqrt(0.6)*(51.96/52.8);
	for(Side s = EAST; s <= WEST; ++s) {
		for(AxisDirection d=X_DIRECTION; d<=Y_DIRECTION; ++d) {
			int flip = (d==X_DIRECTION && s==EAST)?-1:1;
			mwpcPos[s][d] *= wcPosConversion*flip;
			scintPos[s][d] *= wcPosConversion*flip;
		}
	}
	costheta=cos(costheta);
}

//-------------------------------------------

void G4toPMT_SideSwap::doUnits() {
	std::swap(eQ[EAST],eQ[WEST]);
	std::swap(eDep[EAST],eDep[WEST]);
	std::swap(eW[EAST],eW[WEST]);
	for(AxisDirection d=X_DIRECTION; d<=Y_DIRECTION; ++d) {
		std::swap(mwpcPos[EAST][d],mwpcPos[WEST][d]);
		std::swap(scintPos[EAST][d],scintPos[WEST][d]);
	}
	G4toPMT::doUnits();
}

//-------------------------------------------

G4SegmentMultiplier::G4SegmentMultiplier(const SectorCutter& S):
G4toPMT(), SimPositioner(), SC(S), nrots(0), rcurrent(0) {
	SP = this;
	for(unsigned int i=0; i<SC.n; i++) {
		double th = 2*M_PI/SC.ndivs[i];
		vc.push_back(cos(th));
		vs.push_back(sin(th));
	}
}

void G4SegmentMultiplier::calcReweight() {
	G4toPMT::calcReweight();
	physicsWeight *= primRadius()*SC.ndivs[SC.n-1]/(SC.r*SC.ndivs[rcurrent]);
}

void G4SegmentMultiplier::startScan(bool startRandom) {
	nrots = 0;
	G4toPMT::startScan(startRandom);
}

bool G4SegmentMultiplier::nextPoint() {
	if(!nrots) {
		morePts = ProcessedDataScanner::nextPoint();
		reverseCalibrate();
		rcurrent = SC.getRing(SC.sector(primPos[X_DIRECTION],primPos[Y_DIRECTION]));
		rcurrent = rcurrent<SC.n?rcurrent:SC.n-1;
		nrots = SC.ndivs[rcurrent]-1;
		calcReweight();
		nSimmed++;
		nCounted+=simEvtCounts();
	} else {
		reverseCalibrate();
		nrots--;
		nSimmed++;
		nCounted+=simEvtCounts();
	}
	return morePts || nrots;
}

void G4SegmentMultiplier::rotpt(double& x0, double& y0) {
	double x = x0*vc[rcurrent]-y0*vs[rcurrent];
	y0 = x0*vs[rcurrent]+y0*vc[rcurrent];
	x0 = x;
}

void G4SegmentMultiplier::applyOffset(Sim2PMT& S) {
	if(!nrots) return;
	for(Side s = EAST; s <= WEST; ++s) {
		rotpt(S.scintPos[s][X_DIRECTION],S.scintPos[s][Y_DIRECTION]);
		rotpt(S.mwpcPos[s][X_DIRECTION],S.mwpcPos[s][Y_DIRECTION]);
	}
	rotpt(S.primPos[X_DIRECTION],S.primPos[Y_DIRECTION]);
}

void G4SegmentMultiplier::doUnits() {
	if(!nrots) G4toPMT::doUnits();
}

