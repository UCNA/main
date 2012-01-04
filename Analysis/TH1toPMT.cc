#include "TH1toPMT.hh"
#include <TRandom3.h>
#include <climits>

TRandom3 TH1toPMTrand;

TH1toPMT::TH1toPMT(TH1* h, PosGen* P): Sim2PMT(""), mySpectrum(h), PG(P), nToSim(0) {
	assert(PG);
	ePrim = costheta = eW[EAST] = eW[WEST] = 0;
}

void TH1toPMT::doUnits() {
	Side offside = otherSide(genside);
	// set event position
	PG->next();
	for(AxisDirection d=X_DIRECTION; d<=Y_DIRECTION; ++d) {
		primPos[d] = scintPos[genside][d] = mwpcPos[genside][d] = wires[genside][d].center = PG->pos[d];
		scintPos[offside][d] = mwpcPos[offside][d] = wires[offside][d].center = 0;
	}
	// set event energy randomly selected from histogram
	if(!nToSim) { 
		eDep[genside] = eQ[genside] = mySpectrum->GetRandom();
	} else {
		assert(false); // TODO deterministic energy selections
	}
	eW[genside] = 1.0;
	eW[offside] = eDep[offside] = eQ[offside] = 0;
}
