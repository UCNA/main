#include "TH1toPMT.hh"
#include <TRandom3.h>
#include <climits>

TRandom3 TH1toPMTrand;

TH1toPMT::TH1toPMT(TH1* h): Sim2PMT(""), mySpectrum(h), randomPositionRadius(-1), nToSim(0) {
	genpos[X_DIRECTION]=genpos[Y_DIRECTION]=0;
	ePrim = costheta = eW[EAST] = eW[WEST] = 0;
}

void TH1toPMT::doUnits() {
	Side offside = otherSide(genside);
	if(randomPositionRadius>=0) {
		// select random event position in circle
		genpos[X_DIRECTION]=genpos[Y_DIRECTION]=randomPositionRadius;
		while(genpos[X_DIRECTION]*genpos[X_DIRECTION]+genpos[Y_DIRECTION]*genpos[Y_DIRECTION] > randomPositionRadius*randomPositionRadius) {
			genpos[X_DIRECTION] = TH1toPMTrand.Uniform(-randomPositionRadius,randomPositionRadius);
			genpos[Y_DIRECTION] = TH1toPMTrand.Uniform(-randomPositionRadius,randomPositionRadius);
		}
	}
	for(AxisDirection d=X_DIRECTION; d<=Y_DIRECTION; ++d) {
		// set event position
		primPos[d] = scintPos[genside][d] = mwpcPos[genside][d] = wires[genside][d].center = genpos[d];
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
