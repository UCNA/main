/// \file DataScannerExample.cc example code for using MC simulation data
#include "G4toPMT.hh"
#include "CalDBSQL.hh"

int main(int argc, char *argv[]) {
	
	// Geant4 MC data scanner object
	G4toPMT G2P;
	// use data from these MC files (most recent unpolarized beta decay, includes Fermi function spectrum correction)
	// note wildcard * in filename; MC output is split up over many files, but G2P will TChain them together
	G2P.addFile("/home/mmendenhall/geant4/output/LivPhys_495_neutronBetaUnpol_geomC/analyzed_*.root");
	
	// PMT Calibrator loads run-specific energy calibrations info for selected run (14111)
	RunNum rn = 14111;
	PMTCalibrator PCal(rn);
	// set the data scanner to use this calibrator
	G2P.setCalibrator(PCal);
	
	unsigned int nToSim = 50;	// how many triggering events to simulate
	G2P.resetSimCounters();		// reset counters for how many (triggering) events have been simulated
	// start a scan over the data. Argument "true" means start at random offset in file instead of at beginning
	// if you really want this to be random, you will need to seed rand() with something other than default
	// note that it can take many seconds to load the first point of a scan (loading file segment into memory), but will go much faster afterwards.
	G2P.startScan(true);
	
	while(true) {
		// load next point. If end of data is reached, this will loop back and start at the beginning again.
		G2P.nextPoint();
		
		// do whatever analysis you want on this event.
		
		// check the event characteristics on each side
		for(Side s = EAST; s <= WEST; ++s) {
			// get event classification type. TYPE_IV_EVENT means the event didn't trigger this side.
			EventType tp = G2P.fType;
			// skip non-triggering events, or those outside 50mm position cut (you could add your own custom cuts here, if you cared)
			if(tp>=TYPE_IV_EVENT || !G2P.passesPositionCut(s) || G2P.fSide != s)
				continue;
			
			// print out event info, (simulated) reconstructed true energy and position, comparable to values in data
			printf("Event on side %c: type=%i, Etrue=%g @ position (%g,%g)\n",
				   sideNames(s),tp,G2P.getEtrue(),G2P.wires[s][X_DIRECTION].center,G2P.wires[s][Y_DIRECTION].center);
			// print out event primary info, only available in simulation
			printf("\tprimary KE=%g, cos(theta)=%g\n",G2P.ePrim,G2P.costheta);
		}
		
		// break when enough data has been generated.
		if(G2P.nSimmed>=nToSim)
			break;
	}
	
	
	
	return 0;
}
