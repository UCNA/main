/// \file DataScannerExample.cc example code for using MC simulation data
#include "G4toPMT.hh"
#include "CalDBSQL.hh"
#include <TH1F.h>

static double electron_mass = 510.9989; // needed for the physics pf Fierz interference

class FierzHistogram {
  public: 
    double minBin;
    double maxBin;
    unsigned int nBins;
    TH1F *fierz_histogram;
    TH1F *sm_histogram;

    FierzHistogram( double _minBin, double _maxBin, unsigned int _nBins) {
        minBin = _minBin;
        maxBin = _maxBin;
        nBins = _nBins;
        fierz_histogram = new TH1F("fierz_histogram", "", nBins, minBin, maxBin);
        sm_histogram = new TH1F("standard_model_histogram", "", nBins, minBin, maxBin);
    }

/*
    double evaluate(double *x, double*p) {
        double rv = 0;
        rv += p[0] * sm_histogram->GetBinContent(sm_histogram->FindBin(x[0]));        
        rv += p[1] * fierz_histogram->GetBinContent(fierz_histogram->FindBin(x[0]));
        return rv;
    }
*/

    void normalizeHistogram(TH1F* hist) {
        hist->Scale(1/(hist->GetBinWidth(2)*hist->Integral()));
    }
};

FierzHistogram betaSpectrum(0,1000,40);

double evaluate(double *x, double*p) {
    double rv = 0;
    rv += (1 - 0.65 * p[0]) * betaSpectrum.sm_histogram->GetBinContent(betaSpectrum.sm_histogram->FindBin(x[0]-p[1]));        
    rv += p[0] * 0.65 * betaSpectrum.fierz_histogram->GetBinContent(betaSpectrum.fierz_histogram->FindBin(x[0]-p[1]));
    return rv;
}

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

	unsigned int nToSim = 1E6;	// how many triggering events to simulate
	unsigned int nSimmed = 0;	// counter for how many (triggering) events have been simulated
	
    /*
    double minBin = 0;
    double maxBin = 1000;
    unsigned int nBins = 40;
    TH1F *fierz_histogram = new TH1F("fierz_histogram", "", nBins, minBin, maxBin);
    TH1F *sm_histogram = new TH1F("standard_model_histogram", "", nBins, minBin, maxBin);
    */
    //FierzHistogram betaSpectrum(0,1000,40);


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
                betaSpectrum.fierz_histogram->Fill(G2P.getEtrue(s), fierz_weight);
            else
                betaSpectrum.sm_histogram->Fill(G2P.getEtrue(s), 1);

			nSimmed++;
		}
		
		// break when enough data has been generated.
		if(nSimmed>=nToSim)
			break;
	}
    betaSpectrum.normalizeHistogram(betaSpectrum.fierz_histogram);
    betaSpectrum.normalizeHistogram(betaSpectrum.sm_histogram);
    betaSpectrum.fierz_histogram->SetLineColor(3);
    betaSpectrum.sm_histogram->SetLineColor(4);
	
    TCanvas *canvas = new TCanvas("fierz_canvas", "Fierz component of energy spectrum");
    betaSpectrum.fierz_histogram->Draw("");
    betaSpectrum.sm_histogram->Draw("Same");

    /*
        If you want a fast way to get the combined data spectrums for comparison, I already have it extracted as a ROOT histogram for my own data/MC comparisons.
        You can read the ROOT TFile at:
        /home/mmendenhall/mpmAnalyzer/PostPlots/OctetAsym_div0/Combined/Combined.root
        and get the TH1F histograms named:
        Combined_Events_<S><afp>1<type>
        where <S>='E' or 'W' is the side, <afp>='0' or '1' for AFP Off/On, and <type>='0','1','2' for type 0, I, or II/III backscatter events.
        (the '1' in the name after <afp> indicates foreground runs; set to '0' if you want to see the background data).
    */

    TFile *ucna_data_tfile = new TFile(
        "/home/mmendenhall/mpmAnalyzer/PostPlots/OctetAsym_div0/Combined/Combined.root");

    TH1F *ucna_data_histogram = (TH1F*)ucna_data_tfile->Get("Combined_Events_E010");
    ucna_data_histogram->Scale(1/(ucna_data_histogram->GetBinWidth(2)*ucna_data_histogram->Integral()));
    ucna_data_histogram->Draw("Same");
    //printf("Number of bins in data %d\n", ucna_data_histogram->GetNbinsX());

    TF1 *fit = new TF1("fierz_fit", evaluate, 0, 1000, 2);
    fit->SetParameter(0,1);
    fit->SetParameter(1,1);
    ucna_data_histogram->Fit("fierz_fit");
    fit->SetLineColor(6);
    fit->Draw("Same");

    TString pdf_filename = "/data/kevinh/mc/fierz_test.pdf";
    canvas->SaveAs(pdf_filename);

	return 0;
}
