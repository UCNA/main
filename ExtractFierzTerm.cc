/// \file DataScannerExample.cc example code for using MC simulation data
#include "G4toPMT.hh"
#include "CalDBSQL.hh"
#include <TH1F.h>

int main(int argc, char *argv[]) {
	
	// Geant4 MC data scanner object
	G4toPMT G2P;
	// use data from these MC files (most recent unpolarized beta decay, includes Fermi function spectrum correction)
	// note wildcard * in filename; MC output is split up over many files, but G2P will TChain them together
	G2P.addFile("/home/mmendenhall/geant4/output/Baseline_20110826_neutronBetaUnpol_geomC/analyzed_*.root");
	
	// PMT Calibrator loads run-specific energy calibrations info for selected run (14111)
	// and uses default Calibrations DB connection to most up-to-date though possibly unstable "mpm_debug"
	RunNum rn = 14111;
	PMTCalibrator PCal(rn,CalDBSQL::CDB);
	
	// Energy simulators for both sides using same PMT Calibrator
	PMTGenerator PGenE;
	PGenE.addCalibrator(&PCal);
	PMTGenerator PGenW;
	PGenW.addCalibrator(&PCal);
	// set the data scanner to use these PMT Calibrators
	G2P.setGenerators(&PGenE,&PGenW);
	
    double electron_mass = 510.9989; // needed for the physics pf Fierz interference
	unsigned int nToSim = 10000000;	// how many triggering events to simulate
	unsigned int nSimmed = 0;	// counter for how many (triggering) events have been simulated
    TH1F *fierz_histogram = new TH1F("fierz_histogram", "", 128, 0, 900);
    TH1F *sm_histogram = new TH1F("standard_model_histogram", "", 128, 0, 900);

	// start a scan over the data. Argument "true" means start at random offset in file instead of at beginning
	// if you really want this to be random, you will need to seed rand() with something other than default
	// note that it can take many seconds to load the first point of a scan (loading file segment into memory), but will go much faster afterwards.
	G2P.startScan(true);
	


	while(true) {
		// load next point. If end of data is reached, this will loop back and start at the beginning again.
		G2P.nextPoint();
		// perform energy calibrations/simulations to fill class variables with correct values for this simulated event
		G2P.recalibrateEnergy();
		
		// do whatever analysis you want on this event.
        
		// check the event characteristics on each side
		for(Side s = EAST; s <= WEST; s = nextSide(s)) {
			// get event classification type. TYPE_IV_EVENT means the event didn't trigger this side.
			EventType tp = G2P.getEventType(s);
			// skip non-triggering events, or those outside 50mm position cut (you could add your own custom cuts here, if you cared)
			if(tp>=TYPE_IV_EVENT || !G2P.passesPositionCut(s))
				continue;
			
			// print out event info, (simulated) reconstructed true energy and position, comparable to values in data
			printf("Event on side %c: type=%i, Etrue=%g @ position (%g,%g)\n",sideNames(s),tp,G2P.getEtrue(s),G2P.mwpcs[s].pos[0],G2P.mwpcs[s].pos[1]);
			// print out event primary info, only available in simulation
			printf("\tprimary KE=%g, cos(theta)=%g\n",G2P.ePrim,G2P.costheta);

            double energy = G2P.ePrim + electron_mass;
            double fierz_weight = electron_mass / energy;
            if (nSimmed % 2)
                fierz_histogram->Fill(G2P.getEtrue(s), fierz_weight);
            else
                sm_histogram->Fill(G2P.getEtrue(s), 1);

			nSimmed++;
		}
		
		// break when enough data has been generated.
		if(nSimmed>=nToSim)
			break;
	}
    fierz_histogram->Scale(1/fierz_histogram->Integral());
    sm_histogram->Scale(1/sm_histogram->Integral());
    fierz_histogram->SetLineColor(2);
    sm_histogram->SetLineColor(3);
	
    TCanvas *canvas = new TCanvas("fierz_canvas", "Fierz component of energy spectrum");
    fierz_histogram->Draw("");
    sm_histogram->Draw("Same");
    TString gif_filename = "/data/kevinh/mc/fierz_test.gif";
    canvas->SaveAs(gif_filename);
	
	return 0;
}
