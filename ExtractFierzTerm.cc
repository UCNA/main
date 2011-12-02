/// \file DataScannerExample.cc example code for using MC simulation data
#include "G4toPMT.hh"
#include "CalDBSQL.hh"
#include <TH1F.h>
#include <TLegend.h>
//#include <TFitResult.h> // v5.27
#include <TF1.h>

static double electron_mass = 510.9989; // needed for the physics of Fierz interference
static unsigned nToSim = 1E6;	// how many triggering events to simulate

class FierzHistogram {
  public: 
    double minBin;
    double maxBin;
    unsigned int nBins;
    TH1F *fierz_super_sum_histogram;
    TH1F *sm_super_sum_histogram;
    TH1F* fierz_histogram[2][2];
    TH1F* sm_histogram[2][2];

    FierzHistogram( double _minBin, double _maxBin, unsigned int _nBins) {
        minBin = _minBin;
        maxBin = _maxBin;
        nBins = _nBins;
        fierz_super_sum_histogram = new TH1F("fierz_histogram", "", nBins, minBin, maxBin);
        sm_super_sum_histogram = new TH1F("standard_model_histogram", "", nBins, minBin, maxBin);
        for (int side = 0; side < 2; side++)
            for (int spin = 0; spin < 2; spin++) {
                fierz_histogram[side][spin] = new TH1F("fierz_super_sum", "", nBins, minBin, maxBin);
                sm_histogram[side][spin] = new TH1F("standard_model_super_sum", "", nBins, minBin, maxBin);
            }
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

// ug. needs to be static
FierzHistogram mc(0,1000,100);

/**
 * x[0] : kenetic energy
 * p[0] : b, fierz term
 */
double theoretical_fierz_spectrum(double *x, double*p) {
    double rv = 0;
    //unsigned n = mc.sm_histogram->FindBin(p[3]*x[0]*x[0] + p[2]*x[0] - p[1]);        
    unsigned n = mc.sm_super_sum_histogram->FindBin(p[2]*x[0] - p[1]);        
    //unsigned n = mc.sm_histogram->FindBin(x[0]);        
    double b = p[0];
    double norm = 1 + 0.65 * b;
    rv += mc.sm_super_sum_histogram->GetBinContent(n) / norm;
    rv += b * 0.65 * mc.fierz_super_sum_histogram->GetBinContent(n) / norm;
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
    for (int bin = 1; bin < bins+2; bin++) {
        double r[2][2];
        for (int side = 0; side < 2; side++)
            for (int spin = 0; spin < 2; spin++)
                r[side][spin] = rate_histogram[side][spin]->GetBinContent(bin);
        double super_ratio = (r[0][0] * r[1][1]) / (r[0][1] * r[1][0]);
        super_ratio_histogram->SetBinContent(bin, super_ratio);
        super_ratio_histogram->SetBinError(bin, 0.01);
    }
    return super_ratio_histogram;
}

// S = (r[0][0] * r[1][1]) / (r[0][1] * r[1][0]);
TH1F* compute_super_sum(TH1F* rate_histogram[2][2] ) {
    TH1F *super_sum_histogram = new TH1F(*(rate_histogram[0][0]));
    int bins = super_sum_histogram->GetNbinsX();
    for (int bin = 1; bin < bins+2; bin++) {
        double r[2][2];
        for (int side = 0; side < 2; side++)
            for (int spin = 0; spin < 2; spin++)
                r[side][spin] = rate_histogram[side][spin]->GetBinContent(bin);
        double super_sum = TMath::Sqrt(r[0][0] * r[1][1]) + TMath::Sqrt(r[0][1] * r[1][0]);
        if ( TMath::IsNaN(super_sum) ) 
            super_sum = 0;

        printf("Setting bin content for super sum bin %d, to %f\n", bin, super_sum);
        super_sum_histogram->SetBinContent(bin, super_sum);
        super_sum_histogram->SetBinError(bin, TMath::Sqrt(super_sum));
    }
    return super_sum_histogram;
}

TH1F* compute_asymmetry(TH1F* rate_histogram[2][2] ) {
    TH1F *asymmetry_histogram = new TH1F(*(rate_histogram[0][0]));
    int bins = asymmetry_histogram->GetNbinsX();
    for (int bin = 1; bin < bins+2; bin++) {
        double r[2][2];
        for (int side = 0; side < 2; side++)
            for (int spin = 0; spin < 2; spin++)
                r[side][spin] = rate_histogram[side][spin]->GetBinContent(bin);
        double sqrt_super_ratio = TMath::Sqrt((r[0][0] * r[1][1]) / (r[0][1] * r[1][0]));
        if ( TMath::IsNaN(sqrt_super_ratio) ) 
            sqrt_super_ratio = 0;
        double asymmetry = (1 - sqrt_super_ratio) / (1 + sqrt_super_ratio);
        asymmetry_histogram->SetBinContent(bin, asymmetry);
        printf("Setting bin content for super sum bin %d, to %f\n", bin, asymmetry);
        asymmetry_histogram->SetBinError(bin, 0.0);
    }
    return asymmetry_histogram;
}

TH1F* compute_rate_function(TH1F* rate_histogram[2][2], 
                            double (*rate_function)(double r[2][2]))
{
    TH1F *out_histogram = new TH1F(*(rate_histogram[0][0]));
    int bins = out_histogram->GetNbinsX();
    for (int bin = 1; bin < bins+2; bin++) {
        double r[2][2];
        for (int side = 0; side < 2; side++)
            for (int spin = 0; spin < 2; spin++)
                r[side][spin] = rate_histogram[side][spin]->GetBinContent(bin);

        double value = rate_function(r);
        out_histogram->SetBinContent(bin, value);
    }
    return out_histogram;
}

TH1F* compute_rate_function(TH1F* rate_histogram[2][2], 
                            double (*rate_function)(double r[2][2]),
                            double (*error_function)(double r[2][2])) 
{
    TH1F *out_histogram = new TH1F(*(rate_histogram[0][0]));
    int bins = out_histogram->GetNbinsX();
    for (int bin = 1; bin < bins+2; bin++) {
        double r[2][2];
        double e[2][2];
        for (int side = 0; side < 2; side++)
            for (int spin = 0; spin < 2; spin++) {
                r[side][spin] = rate_histogram[side][spin]->GetBinContent(bin);
                e[side][spin] = rate_histogram[side][spin]->GetBinError(bin);
            }
        double value = 0;
        double error = 0; 

        if (rate_function)
            value = rate_function(r);
        if (error_function)
            error = error_function(e);

        out_histogram->SetBinContent(bin, value);
        out_histogram->SetBinError(bin, error);
    }
    return out_histogram;
}

/*
TH1F* compute_rate_error_function(TH1F* rate_histogram[2][2], 
                                  double (*rate_error_function)(double r[2][2])) 
{
    TH1F *out_histogram = new TH1F(*(rate_histogram[0][0]));
    int bins = out_histogram->GetNbinsX();
    for (int bin = 1; bin < bins+2; bin++) {
        double sr[2][2];
        for (int side = 0; side < 2; side++)
            for (int spin = 0; spin < 2; spin++)
                sr[side][spin] = rate_histogram[side][spin]->GetBinError(bin);
    }
    return out_histogram;
}
*/

double bonehead_sum(double r[2][2]) {
    return r[0][0] + r[0][1] + r[1][0] + r[1][1];
}

double bonehead_asymmetry(double r[2][2]) {
    return (r[0][0] - r[0][1])/(r[1][0] + r[1][1]);
}

double super_ratio_asymmetry(double r[2][2]) {
    double super_ratio = (r[0][0] * r[1][1]) / (r[0][1] * r[1][0]);
    double sqrt_super_ratio = TMath::Sqrt(super_ratio);
    if ( TMath::IsNaN(sqrt_super_ratio) ) 
        sqrt_super_ratio = 0;
    return (1 - sqrt_super_ratio) / (1 + sqrt_super_ratio);
}

double super_sum_error(double r[2][2]) {
    double super_ratio = (r[0][0] * r[1][1]) / (r[0][1] * r[1][0]);
    double sqrt_super_ratio = TMath::Sqrt(super_ratio);
    if ( TMath::IsNaN(sqrt_super_ratio) ) 
        sqrt_super_ratio = 0;
    //return (1 - sqrt_super_ratio) / (1 + sqrt_super_ratio);
}

int main(int argc, char *argv[]) {
	
	// Geant4 MC data scanner object
	G4toPMT G2P;
	// use data from these MC files (most recent unpolarized beta decay, includes Fermi function spectrum correction)
	// note wildcard * in filename; MC output is split up over many files, but G2P will TChain them together
	G2P.addFile("/home/mmendenhall/geant4/output/Livermore_neutronBetaUnpol_geomC/analyzed_1.root");
	//G2P.addFile("/home/mmendenhall/geant4/output/Baseline_20110826_neutronBetaUnpol_geomC/analyzed_*.root");
	//G2P.addFile("/home/mmendenhall/mpmAnalyzer/PostPlots/OctetAsym_Offic_10keV_bins/Combined");
    //G2P.addFile("/home/mmendenhall/mpmAnalyzer/PostPlots/OctetAsym_10keV_Bins/Combined");
	
	// PMT Calibrator loads run-specific energy calibrations info for selected run (14111)
	// and uses default Calibrations DB connection to most up-to-date though possibly unstable "mpm_debug"
	RunNum run_number = 14111;
	PMTCalibrator PCal(run_number, CalDBSQL::getCDB());
	
	// Energy simulators for both sides using same PMT Calibrator
	/*
	PMTGenerator PGenE;
	PGenE.setCalibrator(&PCal);
	PMTGenerator PGenW;
	PGenW.setCalibrator(&PCal);
	*/
	//PMTGenerator PGen;
	//PGen.setCalibrator(&PCal);
	// set the data scanner to use these PMT Calibrators
	//G2P.setGenerators(&PGenE,&PGenW);
	G2P.setCalibrator(PCal);

	unsigned int nSimmed = 0;	// counter for how many (triggering) events have been simulated
	
    /*
    double minBin = 0;
    double maxBin = 1000;
    unsigned int nBins = 40;
    TH1F *fierz_histogram = new TH1F("fierz_histogram", "", nBins, minBin, maxBin);
    TH1F *sm_histogram = new TH1F("standard_model_histogram", "", nBins, minBin, maxBin);
    */
    //FierzHistogram mc(0,1000,40);

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
			#ifdef DEBUG
			printf("Event on side %c: type=%i, Etrue=%g @ position (%g,%g), %d\n",
				   sideNames(s), tp, G2P.getEtrue(), G2P.wires[s][X_DIRECTION].center, 
				   G2P.wires[s][Y_DIRECTION].center, (unsigned)G2P.getAFP());

			// print out event primary info, only available in simulation
			printf("\tprimary KE=%g, cos(theta)=%g\n",G2P.ePrim,G2P.costheta);
			#endif 

            /*
            double energy = G2P.ePrim + electron_mass;
            double fierz_weight = electron_mass / energy;
            if (nSimmed % 2)
                mc.fierz_histogram->Fill(G2P.getEtrue(), fierz_weight);
            else
                */
            
            //if (G2P.afp == AFP_ON)
            if (nSimmed % 100 > 40) // redo with real loading eff.
                mc.sm_histogram[s][0]->Fill(G2P.getEtrue(), 1);
            else
                mc.sm_histogram[s][1]->Fill(G2P.getEtrue(), 1);

			nSimmed++;
		}
		
		// break when enough data has been generated.
		if(nSimmed>=nToSim)
			break;
	}
    
    mc.sm_super_sum_histogram = compute_super_sum(mc.sm_histogram);
    normalize(mc.sm_super_sum_histogram);

    mc.fierz_super_sum_histogram = compute_super_sum(mc.fierz_histogram);
    normalize(mc.fierz_super_sum_histogram);

    for (int side = 0; side < 2; side++)
        for (int spin = 0; spin < 2; spin++) {
            normalize(mc.fierz_histogram[side][spin]);
            normalize(mc.sm_histogram[side][spin]);
        }


    TCanvas *canvas = new TCanvas("fierz_canvas", "Fierz component of energy spectrum");

	mc.fierz_super_sum_histogram->SetStats(0);
    mc.fierz_super_sum_histogram->SetLineColor(3);
    mc.fierz_super_sum_histogram->Draw("");
    mc.sm_super_sum_histogram->SetLineColor(1);
    mc.sm_super_sum_histogram->Draw("Same");

    // fit a smooth model to the mc
    /*
    TF1 *mc_fit = new TF1("fierz_mc_fit", mc_model, 0, 1000, deg+1);
    mc_fit->SetParameter(0,-0.5);
    mc_fit->SetParameter(1,1.0);
    mc_fit->SetParameter(2,1.0);
    mc_fit->SetParameter(3,1.0);
    mc_fit->SetParameter(4,1.0);
    //mc_fit->SetParameter(5,1.0);
    //mc_fit->SetParameter(6,1.0);
    mc.sm_histogram->Fit("fierz_mc_fit");

    TString pdf_filename = "/data/kevinh/mc/fierz_models_to_fit.pdf";
    canvas->SaveAs(pdf_filename);
    */

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
        //"/home/mmendenhall/mpmAnalyzer/PostPlots/OctetAsym_div0/Combined/Combined.root");
	    //"/home/mmendenhall/mpmAnalyzer/PostPlots/OctetAsym_Offic_10keV_bins/Combined.root");
        //"/home/mmendenhall/mpmAnalyzer/PostPlots/OctetAsym_10keV_Bins/Combined");
		"/home/mmendenhall/mpmAnalyzer/PostPlots/OctetAsym_10keV_Bins/OctetAsym_10keV_Bins.root");
	if (ucna_data_tfile->IsZombie())
	{
		//printf("File "+beta_filename+"not found.\n");
		std::cout << "File not found." << std::endl;
		exit(1);
	}

    TH1F *ucna_data_histogram[2][2] = {
        {
            //(TH1F*)ucna_data_tfile->Get("Combined_Events_E010"),
            //(TH1F*)ucna_data_tfile->Get("Combined_Events_E110")
            (TH1F*)ucna_data_tfile->Get("hEnergy_Type_0_E_Off"),
            (TH1F*)ucna_data_tfile->Get("hEnergy_Type_0_E_On")
        }, {
            //(TH1F*)ucna_data_tfile->Get("Combined_Events_W010"),
            //(TH1F*)ucna_data_tfile->Get("Combined_Events_W110")
            (TH1F*)ucna_data_tfile->Get("hEnergy_Type_0_W_Off"),
            (TH1F*)ucna_data_tfile->Get("hEnergy_Type_0_W_On")
        }
    };
    //printf("Number of bins in data %d\n", ucna_data_histogram->GetNbinsX());

    /* Already background subtracted...
        TH1F *background_histogram = (TH1F*)ucna_data_tfile->Get("Combined_Events_E000");
        ucna_data_histogram->Add(background_histogram,-1);
        // normalize after background subtraction
        background_histogram->Draw("");
    */
    //normalize(ucna_data_histogram[0][0]);
    if (ucna_data_histogram[0][0] == NULL)
	{
		puts("histogram is null. Aborting...");
		exit(1);
	}

	
    /*
    TF1 *fit = new TF1("fierz_fit", theoretical_fierz_spectrum, 0, 1000, 3);
    fit->SetParameter(0,0.0);
    fit->SetParameter(1,0.0);
    fit->SetParameter(2,1.0);
    ucna_data_histogram[0][0]->Fit("fierz_fit");
    double chisq = fit->GetChisquare();
    double N = fit->GetNDF();
    printf("Chi^2 / ( N - 1) = %f / %f = %f\n",chisq, N-1, chisq/(N-1));
    */

    TString fit_pdf_filename = "/data/kevinh/mc/fierz_fit_data.pdf";
    canvas->SaveAs(fit_pdf_filename);

    // compute and plot the super ratio
    TH1F *super_ratio_histogram = compute_super_ratio(ucna_data_histogram);
    super_ratio_histogram->Draw();
    TString super_ratio_pdf_filename = "/data/kevinh/mc/super_ratio_data.pdf";
    canvas->SaveAs(super_ratio_pdf_filename);

    // compute and plot the super ratio asymmetry 
    TH1F *asymmetry_histogram = compute_asymmetry(ucna_data_histogram);
    asymmetry_histogram->Draw();
    TString asymmetry_pdf_filename = "/data/kevinh/mc/asymmetry_data.pdf";
    canvas->SaveAs(asymmetry_pdf_filename);

    // Compute the super sums
    TH1F *super_sum_histogram = compute_super_sum(ucna_data_histogram);
    normalize(super_sum_histogram);
    super_sum_histogram->SetLineColor(2);
    super_sum_histogram->Draw("");

    // Compute the bonehead sum 
	/* TH1F *bonehead_sum_histogram = compute_rate_function(ucna_data_histogram, &bonehead_sum);
    normalize(bonehead_sum_histogram);
    bonehead_sum_histogram->SetLineColor(45);
    bonehead_sum_histogram->Draw("same"); */

    // Draw Monte Carlo
    mc.sm_super_sum_histogram->SetLineColor(1);
    mc.sm_super_sum_histogram->SetMarkerStyle(4);
    mc.sm_super_sum_histogram->Draw("same p0");

    // lets make a pretty legend
    TLegend * legend = new TLegend(0.6,0.8,0.7,0.6);
    legend->AddEntry(super_sum_histogram, "Super sum", "l");
    legend->AddEntry(mc.sm_super_sum_histogram, "Monte Carlo", "l");
    //legend->AddEntry(bonehead_sum_histogram, "Bonehead sum", "l");
    legend->SetTextSize(0.04);
    legend->SetBorderSize(0);
    legend->Draw();

    // save the data and Mote Carlo plots
    TString super_sum_pdf_filename = "/data/kevinh/mc/super_sum_data.pdf";
    canvas->SaveAs(super_sum_pdf_filename);

    // compute little b factor
    TH1F *fierz_ratio_histogram = new TH1F(*super_sum_histogram);
    fierz_ratio_histogram->Divide(super_sum_histogram, mc.sm_super_sum_histogram);
    fierz_ratio_histogram->GetYaxis()->SetRangeUser(0.6,1.6); // Set the range
    fierz_ratio_histogram->SetTitle("Ratio to Monte Carlo");
	//fierz_ratio_histogram->Fit("pol2", "r", "", 100, 800);
	//TFitResult* fit = ((TFitResultPtr)fierz_ratio_histogram->Fit("1++511/(511+x)", "Sr", "", 100, 800)).Get(); // works in v5.27 ?
	fierz_ratio_histogram->Fit("1++(511/(511+x)-0.65)", "Sr", "", 100, 700);
    //TF1 *mc_fit = new TF1("fierz_mc_fit", mc_model, 0, 1000, deg+1);
    //mc_fit->SetParameter(0,-0.5);

/*
    double chisq = fit->GetChisquare();
    double N = fit->GetNDF();
    printf("Chi^2 / ( N - 1) = %f / %f = %f\n",chisq, N-1, chisq/(N-1));
*/

	fierz_ratio_histogram->SetStats(0);
    fierz_ratio_histogram->Draw();
    TString fierz_ratio_pdf_filename = "/data/kevinh/mc/fierz_ratio.pdf";
    canvas->SaveAs(fierz_ratio_pdf_filename);

	return 0;
}
