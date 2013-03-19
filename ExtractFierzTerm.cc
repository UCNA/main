// UCNA includes
#include "G4toPMT.hh"
#include "PenelopeToPMT.hh"
#include "CalDBSQL.hh"

// ROOT includes
#include <TH1F.h>
#include <TLegend.h>
//#include <TFitResult.h> // v5.27
#include <TF1.h>

// c++ includes
#include <iostream>
#include <fstream>
#include <string>

// c includes
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


///
/// physical constants
///
const double    pi          = 3.1415926535;         /// pi
const double    alpha       = 1/137.035999;         /// fine-structure constant
const double    a2pi        = alpha/(2*pi);         /// alpha over 2pi
const double    m_e         = 510.9988484;          /// electron mass       (14)
const double    m_p         = 938272.046;           /// proton mass         (21)
const double    m_n         = 939565.378;           /// neutron mass        (21)
const double    mu_p        = 2.792846;             /// proton moment       (7)
const double    mu_n        = -1.9130427;           /// neutron moment      (5)
const double    muV         = mu_p - mu_n;          /// nucleon moment      (9)
const double    Q           = 782.344;              /// end point KE        (30)
const double    E0          = m_e + Q;              /// end point E         (30)
const double    L           = -1.27590;             /// gA/gV               (450)
const double    M           = 1 + 3*L*L;            /// matrix element
const double    I_0         = 1.63632;              /// 0th moment          (25)
const double    I_1         = 1.07017;              /// 1st moment          (15)
const double    x_1         = I_1/I_0;              /// first m/E moment    (9)


///
/// settings
///
double min_E = 220;
double max_E = 670;
static double expected_fierz = 0.6540;		/// full range (will get overwritten) 
//static double expected_fierz = 0.6111;	/// for range 150 - 600
//static double expected_gluck = 11.8498;   /// for range 150 - 600
static unsigned nToSim = 1E4;				/// how many triggering events to simulate
static double loading_prob = 40; 			/// ucn loading probability (percent)
static int bins = 150;						/// replace with value from data or smoothing fit
static int integral_size = 1000;
//double scale_x = 1.015;
double scale_x = 1.0;


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

    void normalizeHistogram(TH1F* hist, double min, double max) {
		int _min = hist->FindBin(min);
		int _max = hist->FindBin(max);
        hist->Scale(1/(hist->GetBinWidth(2)*hist->Integral(_min, _max)));
    }
};


// ug. needs to be static
FierzHistogram mc(0,1500,bins);

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
    double norm = 1 + expected_fierz * b;
    rv += mc.sm_super_sum_histogram->GetBinContent(n) / norm;
    rv += b * expected_fierz * mc.fierz_super_sum_histogram->GetBinContent(n) / norm;
    return rv;
}


/// beta spectrum with little b term
double beta_spectrum(const double *val, const double *par) {
	const double K = val[0];                    /// kinetic energy
	const double b = par[0];                    /// Fierz parameter
	if (K <= 0 or K >= Q)
		return 0;                               /// zero outside range

	const double E = K + m_e;                   /// electron energy
	const double e = Q - K;                     /// neutrino energy
	const double p = sqrt(E*E - m_e*m_e);       /// electron momentum
	const double x = m_e/E;                     /// Fierz term
	const double f = (1 + b*x)/(1 + b*x_1);     /// Fierz factor
	const double k = 1.3723803E-11/Q;           /// normalization factor
	const double P = k*p*e*e*E*f;               /// the output PDF value

	return P;
}


/// beta spectrum with little b term
double fierz_beta_spectrum(const double *val, const double *par) {
	const double K = val[0];                    /// kinetic energy
	if (K <= 0 or K >= Q)
		return 0;                               /// zero outside range

	const double E = K + m_e;                   /// electron energy
	const double e = Q - K;                     /// neutrino energy
	const double p = sqrt(E*E - m_e*m_e);       /// electron momentum
	const double x = m_e/E;                     /// Fierz term
	const double k = 1.3723803E-11/Q;           /// normalization factor
	const double P = k*p*e*e*E*x;               /// the output PDF value

	return P;
}


unsigned deg = 4;
double mc_model(double *x, double*p) {
    double _exp = 0;
    double _x = x[0] / m_e;
    for (int i = deg; i >= 0; i--)
        _exp = p[i] + _x * _exp;
    return TMath::Exp(_exp);
}


void normalize(TH1F* hist) {
    hist->Scale(1/(hist->GetBinWidth(2)*hist->Integral()));
}


void normalize(TH1F* hist, double min, double max) {
	int _min = hist->FindBin(min);
	int _max = hist->FindBin(max);
	hist->Scale(1/(hist->GetBinWidth(2)*hist->Integral(_min, _max)));
}


double evaluate_expected_fierz(double min, double max) {
    TH1D *h1 = new TH1D("beta_spectrum", "Beta Spectrum", integral_size, min, max);
    TH1D *h2 = new TH1D("beta_spectrum", "Beta Spectrum", integral_size, min, max);
	for (int i = 0; i < integral_size; i++)
	{
		double K = min + double(i)*(max-min)/integral_size;
		double b = 0;
		double y1 = fierz_beta_spectrum(&K, 0);
		double y2 = beta_spectrum(&K, &b);
		h1->SetBinContent(K, y1);
		h2->SetBinContent(K, y2);
	}
	return h1->Integral(0, integral_size) / h2->Integral(0, integral_size);
}


/// S = (r[0][0] * r[1][1]) / (r[0][1] * r[1][0]);
TH1F* compute_super_ratio(TH1F* rate_histogram[2][2] ) {
    TH1F *super_ratio_histogram = new TH1F(*(rate_histogram[0][0]));
    int bins = super_ratio_histogram->GetNbinsX();
	std::cout << "Number of bins " << bins << std::endl;
    for (int bin = 1; bin < bins+2; bin++) {
        double r[2][2];
        for (int side = 0; side < 2; side++)
            for (int spin = 0; spin < 2; spin++)
                r[side][spin] = rate_histogram[side][spin]->GetBinContent(bin);
        double super_ratio = (r[0][0] * r[1][1]) / (r[0][1] * r[1][0]);
        super_ratio_histogram->SetBinContent(bin, super_ratio);
        super_ratio_histogram->SetBinError(bin, 0.01);   // TODO compute correctly!!
    }
    return super_ratio_histogram;
}


TH1F* compute_super_sum(TH1F* rate_histogram[2][2]) {
    TH1F *super_sum_histogram = new TH1F(*(rate_histogram[0][0]));
    int bins = super_sum_histogram->GetNbinsX();
    for (int bin = 1; bin < bins+2; bin++) {
        double r[2][2];
        for (int side = 0; side < 2; side++)
            for (int spin = 0; spin < 2; spin++)
                r[side][spin] = rate_histogram[side][spin]->GetBinContent(bin);
        double super_sum = TMath::Sqrt(r[0][0] * r[1][1]) + TMath::Sqrt(r[0][1] * r[1][0]);
        double rel_error = TMath::Sqrt( 1/(r[0][0] + r[1][0]) + 1/(r[1][1] * r[0][1]));
        if ( TMath::IsNaN(super_sum)) 
            super_sum = 0;

        if (TMath::IsNaN(rel_error)) 
			rel_error = 0;

        printf("Setting bin content for super sum bin %d, to %f\n", bin, super_sum);
        super_sum_histogram->SetBinContent(bin, super_sum);
        //super_sum_histogram->SetBinError(bin, TMath::Sqrt(super_sum));
        super_sum_histogram->SetBinError(bin, super_sum*rel_error);
    }
    return super_sum_histogram;
}


TH1F* compute_asymmetry(TH1F* rate_histogram[2][2] ) {
    TH1F *asymmetry_histogram = new TH1F(*(rate_histogram[0][0]));
    int bins = asymmetry_histogram->GetNbinsX();
    for (int bin = 1; bin < bins+2; bin++) 
	{
        double r[2][2];
        for (int side = 0; side < 2; side++)
            for (int spin = 0; spin < 2; spin++)
                r[side][spin] = rate_histogram[side][spin]->GetBinContent(bin);
        double sqrt_super_ratio = TMath::Sqrt((r[0][0] * r[1][1]) / (r[0][1] * r[1][0]));
        if ( TMath::IsNaN(sqrt_super_ratio) ) 
            sqrt_super_ratio = 0;
		double denom = 1 + sqrt_super_ratio;
        double asymmetry = (1 - sqrt_super_ratio) / denom;
		double sqrt_inverse_sum = TMath::Sqrt(1/r[0][0] + 1/r[1][1] + 1/r[0][1] + 1/r[1][0]);
		double asymmetry_error = sqrt_inverse_sum * sqrt_super_ratio / (denom * denom);  
        asymmetry_histogram->SetBinContent(bin, asymmetry);
        asymmetry_histogram->SetBinError(bin, asymmetry_error);
        //printf("Setting bin content for asymmetry bin %d, to %f\n", bin, asymmetry);
    }
    return asymmetry_histogram;
}


TH1F* compute_corrected_asymmetry(TH1F* rate_histogram[2][2] ) {
    TH1F *asymmetry_histogram = new TH1F(*(rate_histogram[0][0]));
    int bins = asymmetry_histogram->GetNbinsX();
    for (int bin = 1; bin < bins+2; bin++) 
	{
        double r[2][2];
        for (int side = 0; side < 2; side++)
            for (int spin = 0; spin < 2; spin++)
                r[side][spin] = rate_histogram[side][spin]->GetBinContent(bin);
        double sqrt_super_ratio = TMath::Sqrt((r[0][0] * r[1][1]) / (r[0][1] * r[1][0]));
        if ( TMath::IsNaN(sqrt_super_ratio) ) 
            sqrt_super_ratio = 0;
		double denom = 1 + sqrt_super_ratio;
        double asymmetry = (1 - sqrt_super_ratio) / denom;
		double sqrt_inverse_sum = TMath::Sqrt(1/r[0][0] + 1/r[1][1] + 1/r[0][1] + 1/r[1][0]);
		double asymmetry_error = sqrt_inverse_sum * sqrt_super_ratio / (denom * denom);  
		double K = asymmetry_histogram->GetBinCenter(bin);
		double E = K + m_e;                   /// electron energy
		double p = sqrt(E*E - m_e*m_e);       /// electron momentum
		double beta = p / E;				  /// v/c
        asymmetry_histogram->SetBinContent(bin, -2*asymmetry/beta);
        asymmetry_histogram->SetBinError(bin, 2*asymmetry_error/beta);
        //printf("Setting bin content for corrected asymmetry bin %d, to %f\n", bin, asymmetry);
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


/*
double super_sum_error(double r[2][2]) {
    double super_ratio = (r[0][0] * r[1][1]) / (r[0][1] * r[1][0]);
    double sqrt_super_ratio = TMath::Sqrt(super_ratio);
    if ( TMath::IsNaN(sqrt_super_ratio) ) 
        sqrt_super_ratio = 0;
    //return (1 - sqrt_super_ratio) / (1 + sqrt_super_ratio);
}
*/


using namespace std;
void output_histogram(string filename, TH1F* h, double ax, double ay)
{
	using namespace std;
	ofstream ofs;
	ofs.open(filename.c_str());
	for (int i = 1; i < h->GetNbinsX(); i++)
	{
		double x = ax * h->GetBinCenter(i);
		//double sx = h->GetBinWidth(i);
		double r = ay * h->GetBinContent(i);
		double sr = ay * h->GetBinError(i);
		ofs << x << '\t' << r << '\t' << sr << endl;
	}
	ofs.close();
}


int main(int argc, char *argv[]) 
{
	expected_fierz = evaluate_expected_fierz(min_E, max_E);
	std::cout << "Expected Fierz " << expected_fierz << "\n";
	
	// Geant4 MC data scanner object
	G4toPMT G2P;
	//PenelopeToPMT G2P;

	// use data from these MC files (most recent unpolarized beta decay, includes Fermi function spectrum correction)
	// note wildcard * in filename; MC output is split up over many files, but G2P will TChain them together
	//G2P.addFile("/home/ucna/penelope_output/ndecay_10/event_*.root"); // standard final Penelope
	//G2P.addFile("/home/mmendenhall/geant4/output/20120824_MagF_neutronBetaUnpol/analyzed_*.root"); // magnetic wiggles Monte Carlo
	//G2P.addFile("/home/mmendenhall/geant4/output/20120823_neutronBetaUnpol/analyzed_*.root"); // standard final Monte Carlo 
	//G2P.addFile("/home/mmendenhall/geant4/output/20120810_neutronBetaUnpol/analyzed_*.root");
	//G2P.addFile("/home/mmendenhall/geant4/output/Livermore_neutronBetaUnpol_geomC/analyzed_*.root");
	//G2P.addFile("/home/mmendenhall/geant4/output/Baseline_20110826_neutronBetaUnpol_geomC/analyzed_*.root");
	//G2P.addFile("/home/mmendenhall/mpmAnalyzer/PostPlots/OctetAsym_Offic_10keV_bins/Combined");
    //G2P.addFile("/home/mmendenhall/mpmAnalyzer/PostPlots/OctetAsym_10keV_Bins/Combined");
	G2P.addFile("/data2/mmendenhall/G4Out/2010/20120823_neutronBetaUnpol/analyzed_*.root");
	
	// PMT Calibrator loads run-specific energy calibrations info for selected run (14111)
	// and uses default Calibrations DB connection to most up-to-date though possibly unstable "mpm_debug"
	RunNum run_number = 14111;
	PMTCalibrator PCal(run_number);
	
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
	G2P.startScan(false);
	std::cout << "Finished scan." << std::endl;

	srand ( time(NULL) );

	TChain* tchain = G2P.getChain();  // can be used to for GetEntries()
	int n = tchain->GetEntries();
	std::cout << "Total number of Monte Carlo entries without cuts: " << n << std::endl;

	tchain->SetBranchStatus("*",0);
	tchain->SetBranchStatus("EdepQ",1);
	tchain->SetBranchStatus("Edep",1);
	tchain->SetBranchStatus("MWPCEnergy",1);
	tchain->SetBranchStatus("ScintPos",1);
	tchain->SetBranchStatus("MWPCPos",1);
	tchain->SetBranchStatus("time",1);
	tchain->SetBranchStatus("primTheta",1);
	tchain->SetBranchStatus("primKE",1);
	tchain->SetBranchStatus("primPos",1);
	tchain->SetBranchStatus("EdepSD",1);
	tchain->SetBranchStatus("thetaInSD",1);
	tchain->SetBranchStatus("thetaOutSD",1);
	tchain->SetBranchStatus("keInSD",1);
	tchain->SetBranchStatus("keOutSD",1);
	//G2P.setReadpoints();

	unsigned int nSimmed = 0;	// counter for how many (triggering) events have been simulated
	while(G2P.nextPoint()) { // will stop 
		// load next point. If end of data is reached, this will loop back and start at the beginning again.
		//G2P.nextPoint();

		// perform energy calibrations/simulations to fill class variables with correct values for this simulated event
		//G2P.recalibrateEnergy();
		
		// check the event characteristics on each side
		for(Side s = EAST; s <= WEST; ++s) {
			// get event classification type. TYPE_IV_EVENT means the event didn't trigger this side.
			EventType tp = G2P.fType;

			// skip non-triggering events, or those outside 50mm position cut (you could add your own custom cuts here, if you cared)
			//if(tp>=TYPE_I_EVENT || !G2P.passesPositionCut(s) || G2P.fSide != s)
			if (tp >= TYPE_I_EVENT or !G2P.passesPositionCut(s) or G2P.fSide != s)
				continue;
			
			// print out event info, (simulated) reconstructed true energy and position, comparable to values in data
			#ifdef DEBUG
			printf("Event on side %c: type=%i, Etrue=%g @ position (%g,%g), %d\n",
				   sideNames(s), tp, G2P.getEtrue(), G2P.wires[s][X_DIRECTION].center, 
				   G2P.wires[s][Y_DIRECTION].center, (unsigned)G2P.getAFP());

			// print out event primary info, only available in simulation
			printf("\tprimary KE=%g, cos(theta)=%g\n", G2P.ePrim, G2P.costheta);
			#endif 

			// fill with loading efficiency 
			bool load = (nSimmed % 100 < loading_prob);

			// calculate the energy with a distortion factor
			double energy = scale_x * G2P.getEtrue();
            mc.sm_histogram[s][load]->Fill(energy, 1);

			nSimmed++;
		}
		
		// break when enough data has been generated.
		if(nSimmed >= nToSim)
			break;
	}
    
	std::cout << "Total number of Monte Carlo entries with cuts: " << nSimmed << std::endl;

    mc.sm_super_sum_histogram = compute_super_sum(mc.sm_histogram);
    normalize(mc.sm_super_sum_histogram, min_E, max_E);

    mc.fierz_super_sum_histogram = compute_super_sum(mc.fierz_histogram);
    normalize(mc.fierz_super_sum_histogram, min_E, max_E);

    for (int side = 0; side < 2; side++)
        for (int spin = 0; spin < 2; spin++)
		{
            normalize(mc.fierz_histogram[side][spin], min_E, max_E);
            normalize(mc.sm_histogram[side][spin], min_E, max_E);
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
		//"/home/mmendenhall/mpmAnalyzer/PostPlots/OctetAsym_10keV_Bins/OctetAsym_10keV_Bins.root");
		//"/home/mmendenhall/UCNA/PostPlots/OctetAsym_Offic/OctetAsym_Offic.root");
		"/home/mmendenhall/Plots/OctetAsym_Offic/OctetAsym_Offic.root");
	if (ucna_data_tfile->IsZombie())
	{
		//printf("File "+beta_filename+"not found.\n");
		std::cout << "File not found." << std::endl;
		exit(1);
	}

	#define EVENT_TYPE 0 
    TH1F *ucna_data_histogram[2][2] = {
	#if EVENT_TYPE == 0
        {
            (TH1F*)ucna_data_tfile->Get("hEnergy_Type_0_E_Off"),
            (TH1F*)ucna_data_tfile->Get("hEnergy_Type_0_E_On")
        }, {
            (TH1F*)ucna_data_tfile->Get("hEnergy_Type_0_W_Off"),
            (TH1F*)ucna_data_tfile->Get("hEnergy_Type_0_W_On")
        }
	#endif
	#if EVENT_TYPE == -1
        {
            (TH1F*)ucna_data_tfile->Get("Combined_Events_E010"),
            (TH1F*)ucna_data_tfile->Get("Combined_Events_E110")
        }, {
            (TH1F*)ucna_data_tfile->Get("Combined_Events_W010"),
            (TH1F*)ucna_data_tfile->Get("Combined_Events_W110")
        }
	#endif
    };
    //printf("Number of bins in data %d\n", ucna_data_histogram->GetNbinsX());

    /* Already background subtracted...
        TH1F *background_histogram = (TH1F*)ucna_data_tfile->Get("Combined_Events_E000");
        ucna_data_histogram->Add(background_histogram,-1);
        // normalize after background subtraction
        background_histogram->Draw("");
    */

	for (int side = 0; side < 2; side++)
		for (int spin = 0; spin < 2; spin++)
		{
			std::cout << "Number of entries in (" 
					  << side << ", " << spin << ") is "
					  << (int)ucna_data_histogram[side][spin]->GetEntries() << std::endl;
			if (ucna_data_histogram[side][spin] == NULL)
			{
				puts("histogram is null. Aborting...");
				exit(1);
			}
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
    TH1F *asymmetry_histogram = compute_corrected_asymmetry(ucna_data_histogram);

	// fit the Fierz term from the asymmetry
	char A_fit_str[1024];
    sprintf(A_fit_str, "[0]/(1+[1]*(%f/(%f+x)))", m_e, m_e);
    TF1 *A_fierz_fit = new TF1("A_fierz_fit", A_fit_str, min_E, max_E);
    A_fierz_fit->SetParameter(-.12,0);
    A_fierz_fit->SetParameter(1,0);
	asymmetry_histogram->Fit(A_fierz_fit, "Sr");
	asymmetry_histogram->SetMaximum(0);
	asymmetry_histogram->SetMinimum(-0.2);

	// compute chi squared
    double chisq = A_fierz_fit->GetChisquare();
    double N = A_fierz_fit->GetNDF();
	char A_str[1024];
	sprintf(A_str, "A = %1.3f #pm %1.3f", A_fierz_fit->GetParameter(0), A_fierz_fit->GetParError(0));
	char A_b_str[1024];
	sprintf(A_b_str, "b = %1.3f #pm %1.3f", A_fierz_fit->GetParameter(1), A_fierz_fit->GetParError(1));
	char A_b_chisq_str[1024];
    printf("Chi^2 / ( N - 1) = %f / %f = %f\n", chisq, N-1, chisq/(N-1));
	sprintf(A_b_chisq_str, "#frac{#chi^{2}}{n-1} = %f", chisq/(N-1));

	// draw the ratio plot
	asymmetry_histogram->SetStats(0);
    asymmetry_histogram->Draw();

	// draw a legend on the plot
    TLegend* asym_legend = new TLegend(0.3,0.85,0.6,0.65);
    asym_legend->AddEntry(asymmetry_histogram, "Asymmetry data", "l");
    asym_legend->AddEntry(A_fierz_fit, "Fierz term fit", "l");
    asym_legend->AddEntry((TObject*)0, A_str, "");
    asym_legend->AddEntry((TObject*)0, A_b_str, "");
    asym_legend->AddEntry((TObject*)0, A_b_chisq_str, "");
    asym_legend->SetTextSize(0.03);
    asym_legend->SetBorderSize(0);
    asym_legend->Draw();

    TString asymmetry_pdf_filename = "/data/kevinh/mc/asymmetry_data.pdf";
    canvas->SaveAs(asymmetry_pdf_filename);

    // Compute the super sums
    TH1F *super_sum_histogram = compute_super_sum(ucna_data_histogram);
	std::cout << "Number of super sum entries " << (int)super_sum_histogram->GetEntries() << std::endl;
    normalize(super_sum_histogram, min_E, max_E);
    super_sum_histogram->SetLineColor(2);
	super_sum_histogram->SetStats(0);
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

    // make a pretty legend
    TLegend * legend = new TLegend(0.6,0.8,0.7,0.6);
    legend->AddEntry(super_sum_histogram, "Type 0 supersum", "l");
    legend->AddEntry(mc.sm_super_sum_histogram, "Monte Carlo supersum", "p");
    //legend->AddEntry(bonehead_sum_histogram, "Bonehead sum", "l");
    legend->SetTextSize(0.03);
    legend->SetBorderSize(0);
    legend->Draw();

    // save the data and Mote Carlo plots
    TString super_sum_pdf_filename = "/data/kevinh/mc/super_sum_data.pdf";
    canvas->SaveAs(super_sum_pdf_filename);

    // compute little b factor
    TH1F *fierz_ratio_histogram = new TH1F(*super_sum_histogram);
    //fierz_ratio_histogram->Divide(super_sum_histogram, mc.sm_super_sum_histogram);
    int bins = super_sum_histogram->GetNbinsX();
    for (int bin = 1; bin < bins+1; bin++) {
		double X = super_sum_histogram->GetBinContent(bin);
		double Y = mc.sm_super_sum_histogram->GetBinContent(bin);
		double Z = X/Y;
		if (Y > 0)
			fierz_ratio_histogram->SetBinContent(bin, X/Y);
		else
			fierz_ratio_histogram->SetBinContent(bin, 0);

		double x = super_sum_histogram->GetBinError(bin);
		double y = mc.sm_super_sum_histogram->GetBinError(bin);
		fierz_ratio_histogram->SetBinError(bin, Z*TMath::Sqrt(x*x/X/X + y*y/Y/Y));
	}
    fierz_ratio_histogram->GetYaxis()->SetRangeUser(0.9,1.1); // Set the range
    fierz_ratio_histogram->SetTitle("Ratio of UCNA data to Monte Carlo");
	std::cout << super_sum_histogram->GetNbinsX() << std::endl;
	std::cout << mc.sm_super_sum_histogram->GetNbinsX() << std::endl;

	// fit the Fierz ratio 
	char fit_str[1024];
    sprintf(fit_str, "1+[0]*(%f/(%f+x)-%f)", m_e, m_e, expected_fierz);
    TF1 *fierz_fit = new TF1("fierz_fit", fit_str, min_E, max_E);
    fierz_fit->SetParameter(0,0);
	fierz_ratio_histogram->Fit(fierz_fit, "Sr");

	// A fit histogram for output to gnuplot
    TH1F *fierz_fit_histogram = new TH1F(*super_sum_histogram);
	for (int i = 0; i < fierz_fit_histogram->GetNbinsX(); i++)
		fierz_fit_histogram->SetBinContent(i, fierz_fit->Eval(fierz_fit_histogram->GetBinCenter(i)));

	// compute chi squared
    chisq = fierz_fit->GetChisquare();
    N = fierz_fit->GetNDF();
	char b_str[1024];
	sprintf(b_str, "b = %1.3f #pm %1.3f", fierz_fit->GetParameter(0), fierz_fit->GetParError(0));
	char chisq_str[1024];
    printf("Chi^2 / ( N - 1) = %f / %f = %f\n", chisq, N-1, chisq/(N-1));
	sprintf(chisq_str, "#frac{#chi^{2}}{n-1} = %f", chisq/(N-1));

	// draw the ratio plot
	fierz_ratio_histogram->SetStats(0);
    fierz_ratio_histogram->Draw();

	// draw a legend on the plot
    TLegend* ratio_legend = new TLegend(0.3,0.85,0.6,0.65);
    ratio_legend->AddEntry(fierz_ratio_histogram, "Data ratio to Monte Carlo (Type 0)", "l");
    ratio_legend->AddEntry(fierz_fit, "Fierz term fit", "l");
    ratio_legend->AddEntry((TObject*)0, "1+b(#frac{m_{e}}{E} - #LT #frac{m_{e}}{E} #GT)", "");
    ratio_legend->AddEntry((TObject*)0, b_str, "");
    ratio_legend->AddEntry((TObject*)0, chisq_str, "");
    ratio_legend->SetTextSize(0.03);
    ratio_legend->SetBorderSize(0);
    ratio_legend->Draw();

	// output for root
    TString fierz_ratio_pdf_filename = "/data/kevinh/mc/fierz_ratio.pdf";
    canvas->SaveAs(fierz_ratio_pdf_filename);

	// output for gnuplot
	output_histogram("/data/kevinh/mc/super-sum-data.dat", super_sum_histogram, 1, 1000);
	output_histogram("/data/kevinh/mc/super-sum-mc.dat", mc.sm_super_sum_histogram, 1, 1000);
	output_histogram("/data/kevinh/mc/fierz-ratio.dat", fierz_ratio_histogram, 1, 1);
	output_histogram("/data/kevinh/mc/fierz-fit.dat", fierz_fit_histogram, 1, 1);

	return 0;
}
