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
    //TH1F *ucna_data_histogram;
    //TH1F *backgrund_histogram;

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

FierzHistogram betaSpectrum(0,1000,400);

/**
 * x[0] : kenetic energy
 * p[0] : b, fierz term
 */
double theoretical_fierz_spectrum(double *x, double*p) {
    double rv = 0;
    //unsigned n = betaSpectrum.sm_histogram->FindBin(p[3]*x[0]*x[0] + p[2]*x[0] - p[1]);        
    unsigned n = betaSpectrum.sm_histogram->FindBin(p[2]*x[0] - p[1]);        
    //unsigned n = betaSpectrum.sm_histogram->FindBin(x[0]);        
    double b = p[0];
    double norm = 1 + 0.65 * b;
    rv += betaSpectrum.sm_histogram->GetBinContent(n) / norm;
    rv += b * 0.65 * betaSpectrum.fierz_histogram->GetBinContent(n) / norm;
    return rv;
}

unsigned deg = 4;
double mc_model(double *x, double*p) {
    double _exp = 0;
    double _x = x[0] / electron_mass;
    for (int i = deg; i >= 0; i--)
        _exp = p[i] + _x * _exp;
    return TMath::Exp(_exp);
}

void normalize(TH1F* hist) {
    hist->Scale(1/(hist->GetBinWidth(2)*hist->Integral()));
}

// S = (r[0][0] * r[1][1]) / (r[0][1] * r[1][0]);
TH1F* compute_super_ratio(TH1F* rate_histogram[2][2] ) {
    TH1F *super_ratio_histogram = new TH1F(*(rate_histogram[0][0]));
    int bins = super_ratio_histogram->GetNbinsX();
    for (int bin = 0; bin < bins; bin++) {
        double r[2][2];
        for (int side = 0; side < 2; side++)
            for (int spin = 0; spin < 2; spin++)
                r[side][spin] = rate_histogram[side][spin]->GetBinContent(bin);
        double super_ratio = (r[0][0] * r[1][1]) / (r[0][1] * r[1][0]);
        super_ratio_histogram->SetBinContent(bin, super_ratio);
    }
    return super_ratio_histogram;
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

	unsigned int nToSim = 1E5;	// how many triggering events to simulate
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
		
		// check the event characteristics on each side
		for(Side s = EAST; s <= WEST; s = nextSide(s)) {
			// get event classification type. TYPE_IV_EVENT means the event didn't trigger this side.
			EventType tp = G2P.fType;

			// skip non-triggering events, or those outside 50mm position cut (you could add your own custom cuts here, if you cared)
			if(tp>=TYPE_IV_EVENT || !G2P.passesPositionCut(s) || G2P.fSide != s)
				continue;
			
			// print out event info, (simulated) reconstructed true energy and position, comparable to values in data
			printf("Event on side %c: type=%i, Etrue=%g @ position (%g,%g)\n",
				   sideNames(s), tp, G2P.getEtrue(), G2P.wires[s][X_DIRECTION].center, G2P.wires[s][Y_DIRECTION].center);

			// print out event primary info, only available in simulation
			printf("\tprimary KE=%g, cos(theta)=%g\n",G2P.ePrim,G2P.costheta);

            double energy = G2P.ePrim + electron_mass;
            double fierz_weight = electron_mass / energy;
            if (nSimmed % 2)
                betaSpectrum.fierz_histogram->Fill(G2P.getEtrue(), fierz_weight);
            else
                betaSpectrum.sm_histogram->Fill(G2P.getEtrue(), 1);

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

    // fit a smooth model to the mc
    TF1 *mc_fit = new TF1("fierz_mc_fit", mc_model, 0, 1000, deg+1);
    mc_fit->SetParameter(0,-0.5);
    mc_fit->SetParameter(1,1.0);
    mc_fit->SetParameter(2,1.0);
    mc_fit->SetParameter(3,1.0);
    mc_fit->SetParameter(4,1.0);
    //mc_fit->SetParameter(5,1.0);
    //mc_fit->SetParameter(6,1.0);
    betaSpectrum.sm_histogram->Fit("fierz_mc_fit");

    TString pdf_filename = "/data/kevinh/mc/fierz_models_to_fit.pdf";
    canvas->SaveAs(pdf_filename);

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

    TH1F *ucna_data_histogram[2][2] = {
        {
            (TH1F*)ucna_data_tfile->Get("Combined_Events_E010"),
            (TH1F*)ucna_data_tfile->Get("Combined_Events_E110")
        }, {
            (TH1F*)ucna_data_tfile->Get("Combined_Events_W010"),
            (TH1F*)ucna_data_tfile->Get("Combined_Events_W110")
        }
    };
    //printf("Number of bins in data %d\n", ucna_data_histogram->GetNbinsX());

    /* Already background subtracted...
        TH1F *background_histogram = (TH1F*)ucna_data_tfile->Get("Combined_Events_E000");
        ucna_data_histogram->Add(background_histogram,-1);
        // normalize after background subtraction
        background_histogram->Draw("");
    */
    normalize(ucna_data_histogram[0][0]);
    ucna_data_histogram[0][0]->Draw("Same");

    TF1 *fit = new TF1("fierz_fit", theoretical_fierz_spectrum, 0, 1000, 3);
    fit->SetParameter(0,0.0);
    fit->SetParameter(1,0.0);
    fit->SetParameter(2,1.0);
    ucna_data_histogram[0][0]->Fit("fierz_fit");
    double chisq = fit->GetChisquare();
    double N = fit->GetNDF();
    printf("Chi^2 / ( N - 1) = %f / %f = %f\n",chisq, N-1, chisq/(N-1));

    TString fit_pdf_filename = "/data/kevinh/mc/fierz_fit_data.pdf";
    canvas->SaveAs(fit_pdf_filename);

    TH1F *super_ratio_histogram = compute_super_ratio(ucna_data_histogram);

	return 0;
}
