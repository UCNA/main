/// UCNA includes
#include "G4toPMT.hh"
#include "PenelopeToPMT.hh"
#include "CalDBSQL.hh"
#include "FierzFitter.hh"

/// ROOT includes
#include <TH1.h>
#include <TLegend.h>
#include <TF1.h>
#include <TVirtualFitter.h>
#include <TList.h>
#include <TStyle.h>
#include <TApplication.h>
#include <TMatrixD.h>
#include <TNtuple.h>
#include <TLeaf.h>
#include <TString.h>


/// C++ includes
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <limits>
#include <ctime>

/// C includes
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

<<<<<<< Updated upstream
/// name spaces
using std::setw;
using std::cout;
using namespace TMath;
double infity = std::numeric_limits<double>::infinity();

/// cuts and settings
double KEmin = 0;                 /// min kinetic energy for plots
double KEmax = 800;               /// max kinetic range for plots
int    KEbins=(KEmax-KEmin)/10;   /// number of bins to use fit spectral plots
double fit_min = 150;             /// min kinetic energy for plots
//ucna.fidcut2 = fedutial_cut*fedutial_cut;
double fit_max = 630;             /// max kinetic range for plots
int fit_bins=(fit_max-fit_min)/10;/// number of bins to use fit spectral plots
double fedutial_cut = 50;         /// radial cut in millimeters TODO!! HARD CODED IN MODEL
int file_number_start = 1;
int file_number_stop = 3;//99;
int file_number_batch = 1;//10;

/// set up free fit parameters with best guess
#if 0
static TString FIT_TYPE = "AbN";
static const int nPar = 3;
double afp_ratio = 0.40;
TString paramNames[3] = {"A", "b", "N"};
double paramInits[3] = {-0.12, 0.0, 1};
#endif
#if 0
static TString FIT_TYPE = "bN";
static const int nPar = 2;
double afp_ratio = 0.40;
TString paramNames[2] = {"b", "N"};
double paramInits[2] = {0.0, 1000000};
int A_index=-1, b_index=0, N_index=1;
#endif
#if 1
static TString FIT_TYPE = "b";
static const int nPar = 1;
double afp_ratio = 0.40;
TString paramNames[1] = {"b"};
double paramInits[1] = {0};
int A_index=-1, b_index=0, N_index=-1;
#endif
#if 0
static TString FIT_TYPE = "Ab";
static const int nPar = 2;
double afp_ratio = 0.40;
TString paramNames[2] = {"A", "b"};
double paramInits[2] = {1, 0};
int A_index=0, b_index=1, N_index=-1;
#endif
=======
/// settings
double expected_fierz = 0.6540;				/// full range (will get overwritten) 
//static double expected_fierz = 0.6111;	/// for range 150-600
//static double expected_gluck = 11.8498;   /// for range 150-600
static unsigned nToSim = 5e7;				/// how many triggering events to simulate
static double afp_on_prob = 0.68/1.68; 		/// afp on probability per neutron
static double afp_off_prob = 1/1.68; 	    /// afp off probability per neutron
static unsigned rand_bins = 100; 	        /// probability granularity
int large_prime = 179430331;                /// a large prime number
static int bins = 150;						/// replace with value from data or smoothing fit
double scale_x = 1.00000;

double min_E = 220;                         /// min energy from the 2013 paper
double max_E = 670;                         /// max range from the 2013 paper
double fedutial_cut = 50;                   /// radial cut in millimeters 
double fidcut2 = 50*50;                     /// mm^2 radial cut
>>>>>>> Stashed changes

double expected[3][3];                      /// expected values (based on the energy range)
                                            /// need to be visible to the chi^2 code

<<<<<<< Updated upstream
/// BASIC DEFAULT FILENAME LANDSCAPE

/// path to experiment data files
TString data_dir = "/media/hickerson/boson/Data"; 
TString data_filename = "OctetAsym_Offic.root";

/// path to Monte Carlo files
TString mc_dir = "/home/xuansun/Documents/SimProcessedFiles/100mill_beta";
TString mc_sys_dir = "/home/xuansun/Documents/SimProcessedFiles/100mill_beta";
TString mc_name = "SimAnalyzed";
TString mc_filename_stem = "SimAnalyzed-beta-";
TString mc_sys_filename_stem = "SimAnalyzed-beta-twittled-";
=======
// ug. Needs to be static
FierzHistogram mc(0,1500,bins);
>>>>>>> Stashed changes

/// path to save output plots
TString plots_dir = "/home/hickerson/Root";

/// path to save output root structures 
TString root_output_dir = "/home/hickerson/Root";



/// GLOBAL MODELS

/// This needs to be static and global for MINUIT to work
UCNAFierzFitter* global_ff = 0;

void combined_chi2(Int_t & /*n*/, Double_t * /*grad*/ , Double_t &chi2, Double_t *p, Int_t /*iflag */  )
{
    assert(global_ff);
    if (FIT_TYPE=="AbN") {
        double A=p[0], b=p[1], N=p[2];
        chi2 = global_ff->combined_chi2(A,b,N);
    }
    else if (FIT_TYPE=="bN") {
        double b=p[0], N=p[1];
        chi2 = global_ff->supersum_chi2(b,N);
    }
    else if (FIT_TYPE=="b") {
        double b=p[0];
        chi2 = global_ff->supersum_chi2(b,1/(1+b*x_1));
        //chi2 = global_ff->supersum_chi2(b,1);
    }
}

int cl = 14;
void output_matrix(TString title, TMatrixD matrix, TString colNames[], TString rowNames[])
{
    cout<<"\n "<<title;
    int cols = matrix.GetNcols();
    int rows = matrix.GetNrows();
    if (cols>1) {
        cout<<"\n     ";
        for (int i=0; i<cols; i++)
            cout<<setw(cl)<<colNames[i];
    }
    cout<<"\n";
    for (int i=0; i<cols; i++) {
        cout<<"    "<<rowNames[i];
        for (int j=0; j<rows; j++)
            cout<<setw(cl)<<matrix[i][j];
        cout<<"\n";
    }
}

void output_matrix(TString title, TMatrixD matrix)
{
    output_matrix(title, matrix, paramNames, paramNames);
}

/// FITTING 
TF1* fit(UCNAFierzFitter *ff) {
    assert(ff);
    /// Set up the fit functions parameters.
    TF1* asymmetry_func = new TF1("asymmetry_fit_func", &asymmetry_fit_func, fit_min, fit_max, nPar);
    asymmetry_func->SetParameters(paramInits);
    for (int i=0; i<nPar; i++)
        asymmetry_func->SetParName(i, paramNames[i]);

    /*TF1 supersum_func("supersum_fit_func", &supersum_fit_func, fit_min, fit_max, nPar);
    supersum_func.SetParameters(paramInits);
    for (int i=0; i<nPar; i++)
        supersum_func->SetParName(i, paramNames[i]);*/

    /// Actually do the combined fitting.
    global_ff = ff; /// TODO cannot be called in a member
    TMatrixD cov(nPar,nPar);
    ff->combined_fit(cov, asymmetry_func, /*supersum_func,*/ &combined_chi2); /// TODO no member
    //ff->compute_fit(&asymmetry_func);
    global_ff = 0;

    /// Set up reasonable guesses 
    double A = -0.12;
    double b = 0;
    double N = 1;

    double *p = asymmetry_func->GetParameters();
    if (FIT_TYPE=="AbN") {
        assert(nPar == 3);
        A = p[0]; b = p[1]; N = p[2];
    }
    else if (FIT_TYPE=="bN") {
        assert(nPar == 2);
        b = p[0]; N = p[1];
    }
    else if (FIT_TYPE=="b") {
        assert(nPar == 1);
        b = p[0];
        //N = 1/(1+b*x_1));
    }

    /// Compute a fit from the parameters
    ff->compute_fit(A,b,N);

    /// PRINT OUT REPORT OF FITS, CORRELATIONS AND ERRORS

<<<<<<< Updated upstream
    /// Output the data and cut info.
    cout<<" NUMBER OF EVENTS WITH ENERGY RANGE CUTS:\n";
    N = ff->comupte_sizes();
    ff->print_sizes();
    
    /// Set all expectation values for this range.
    double ex_cos = 0.5; // TODO evaluate_expected_cos_theta(fit_min,fit_max);
    double ec2 = ex_cos*ex_cos; // TODO evaluate_expected_cos_theta(fit_min,fit_max);
    const int nSpec = 4;
    TMatrixD ex(nSpec,nSpec);
    TString colNames[nSpec], rowNames[nSpec];
    for (int m=0; m<nSpec; m++) {       /// beta exponent
        colNames[m].Form("(v/c)^%d",m);
        for (int n=0; n<nSpec; n++) {  /// m_e/E exponent
            rowNames[n].Form("(m/E)^%d",n); 
            ex[m][n] = evaluate_expected_fierz(m,n,fit_min,fit_max);
            //cout<<"Expected value for beta exponent: "<<m
            //    <<" and m/E exponent "<<n<<" is "<<ex[m][n]<<".\n";
        }
    }
    output_matrix("APPROXIMATE ANALYTICAL EXPECTATION MATRIX", ex, colNames, rowNames);
    
    /// Calculate the predicted inverse covariance matrix for this range.
    TMatrixD est_cov_inv(nPar,nPar);
    for (int i=0; i<nPar; i++)
        for (int j=0; j<nPar; j++)
            est_cov_inv[i][j] = 0;

    if (A_index >= 0)
        est_cov_inv[A_index][A_index] += N*ec2*ex[2][0];

    if (b_index >= 0) {
        est_cov_inv[b_index][b_index] += N*(ex[0][2] - ex[0][1]*ex[0][1]);
        if (A_index >= 0) {
            est_cov_inv[A_index][b_index] = 
            est_cov_inv[b_index][A_index] = -N*A*ec2*ex[2][1];
            est_cov_inv[b_index][b_index] += N*A*A*ec2*ex[2][2];
        }
=======
/**
 * Si := (r[0][0] r[1][1]) / (r[0][1] * r[1][0])
 */
TH1F* compute_super_ratio(TH1F* rate_histogram[2][2] ) 
{
    TH1F *super_ratio_histogram = new TH1F(*(rate_histogram[0][0]));
    int bins = super_ratio_histogram->GetNbinsX();
	std::cout << "Number of bins " << bins << std::endl;
    for (int bin = 1; bin < bins+2; bin++) {
        double r[2][2];
        for (int side = 0; side < 2; side++)
            for (int spin = 0; spin < 2; spin++)
                r[side][spin] = rate_histogram[side][spin]->GetBinContent(bin);
        double super_ratio = r[0][0]*r[1][1]/r[0][1]/r[1][0];
        super_ratio_histogram->SetBinContent(bin, super_ratio);
        super_ratio_histogram->SetBinError(bin, 0.01);   // TODO compute correctly!!
    }
    return super_ratio_histogram;
}


/**
 * S := r[0][0] + r[1][1] + r[0][1] + r[1][0]
 */
TH1F* compute_super_sum(TH1F* rate_histogram[2][2]) 
{
    TH1F *super_sum_histogram = new TH1F(*(rate_histogram[0][0]));
    int bins = super_sum_histogram->GetNbinsX();
    for (int bin = 1; bin < bins+2; bin++) {
        double r[2][2];
        for (int side = 0; side < 2; side++)
            for (int spin = 0; spin < 2; spin++)
                r[side][spin] = rate_histogram[side][spin]->GetBinContent(bin);
        double super_sum = TMath::Sqrt(r[0][0] * r[1][1]) + TMath::Sqrt(r[0][1] * r[1][0]);
        double rel_error = TMath::Sqrt( 1/r[0][0] + 1/r[1][0] + 1/r[1][1] + 1/r[0][1]);
        if ( TMath::IsNaN(super_sum)) 
            super_sum = 0;

        if (TMath::IsNaN(rel_error)) 
			rel_error = 0;

        printf("Setting bin content for super sum bin %d, to %f\n", bin, super_sum);
        super_sum_histogram->SetBinContent(bin, super_sum);
        super_sum_histogram->SetBinError(bin, super_sum*rel_error);
>>>>>>> Stashed changes
    }

    if (N_index >= 0) {
        est_cov_inv[N_index][N_index] = N;
        if (A_index >= 0) {
            est_cov_inv[A_index][N_index] =
            est_cov_inv[N_index][A_index] = N*ex[0][1];
        }
    }

    /* if (b_index > -1)
        est_cov_inv[0][0] = N*ex_cos*ex[2][0];
    if (nPar > 1) {
        est_cov_inv[1][0] = 
        est_cov_inv[0][1] = -N*A*ex_cos*ex[2][1];
        //est_cov_inv[1][1] =  N*(A*A*ex[2][2]/4 + ex[0][2] - ex[0][1]*ex[0][1]);
        est_cov_inv[1][1] = N*(A*A*ex_cos*ex[2][2] + ex[0][2] - ex[0][1]*ex[0][1]);
    }
    if (nPar > 2) {
        est_cov_inv[1][2] =
        est_cov_inv[2][1] = N*ex[0][1];
        est_cov_inv[2][2] = N;
    } */

    /// Output the fit covariance matrix from the fit.
    output_matrix("FIT COVARIANCE MATRIX", cov);

    /// Compute the estimated covariance matrix by inverse.
    double det = 0;
    TMatrixD est_cov = est_cov_inv;
    output_matrix("ESTIMATED COVARIANCE MATRIX", est_cov.Invert(&det));

    /// Output the estimated covariance details.
    TMatrixD ratio_cov(nPar,nPar);
    for (int i=0; i<nPar; i++)
        for (int j=0; j<nPar; j++)
            ratio_cov[i][j] = est_cov[i][j]/cov[i][j];
    output_matrix("ESTIMATED/ACTUAL COVARIANCE RATIOS", ratio_cov);

    /// Output the fit covariance inverse matrix time 2.
    TMatrixD cov_inv = cov;
    output_matrix("FIT HALF CHI^2 DERIVATIVE MATRIX", cov_inv.Invert(&det));

    /// Compute the estimated covariance matrix by inverse.
    output_matrix("ESTIMATED HALF CHI^2 DERIVATIVE MATRIX", est_cov_inv);

    /// Output the estimated covariance details.
    TMatrixD ratio_cov_inv(nPar,nPar);
    for (int i=0; i<nPar; i++)
        for (int j=0; j<nPar; j++)
            ratio_cov_inv[i][j] = est_cov_inv[i][j]/cov_inv[i][j];
    output_matrix("ESTIMATED/ACTUAL HALF CHI^2 DERIVATIVE RATIOS", ratio_cov_inv);

#if 0
    /// Compute independent standard errors.
    cout<<"\n FOR UNCOMBINED FITS:\n";
    cout<<"    Expected independent statistical sigma and error:\n";
    cout<<"    "<<setw(1)<<" " 
                <<setw(cl)<<"value"
                <<setw(cl)<<"sigma" 
                <<setw(cl+1)<<"error\n";
    for (int i=0; i<nPar; i++) {
        TString name = paramNames[i];
        mc_dir+"/SimAnalyzed_Beta_7.root",
        double value = asymmetry_func->GetParameter(i);
        double sigma = 1/Sqrt(est_cov_inv[i][i]);
        double error = 100*sigma/value;
        cout<<"    "<<setw(1) <<name
                    <<setw(cl)<<value
                    <<setw(cl)<<sigma
                    <<setw(cl-1)<<error<<"%\n";
    }
#endif

    /// Compare predicted and actual standard errors.
    cout<<"\n FOR COMBINED FITS:\n";
    cout<<"    Actual and estimated combined statistical sigmas and errors:\n";
    cout<<"    "<<setw(1)<<" "
                <<setw(cl)<<"value" 
                <<setw(3*cl/2)<<"sigma (act./est)"
                <<setw(3*cl/2)<<"factor (act./est.)"
                <<setw(3*cl/2)<<"error (act./est.)" 
                <<setw(cl+1)<<"est./actual\n";
    for (int i=0; i<nPar; i++) {
        TString param = paramNames[i];
        double scale = Sqrt(N);
        double value = asymmetry_func->GetParameter(i);
        double sigma = Sqrt(cov[i][i]);
        double error = sigma*100/value;
        double factor = scale*sigma;
        double est_sigma = Sqrt(est_cov[i][i]);
        double est_error = est_sigma*100/value;
        double est_factor = scale*est_sigma;
        double ratio = est_sigma/sigma;
        cout<<"    "<<setw(1)   <<param 
                    <<setw(cl)  <<value 
                    <<setw(cl-1)<<sigma<<"/"
                    <<setw(cl/2)<<est_sigma 
                    <<setw(cl-1)<<factor<<"/"
                    <<setw(cl/2)<<est_factor 
                    <<setw(cl-2)<<error <<"%/"
                    <<setw(cl/2-1)<<est_error<<"%"
                    <<setw(cl)  <<ratio<<"\n";
    }

    /// Compare predicted and actual correlations.
    if (nPar > 1)
        cout<<"\n CORRELATIONS FACTORS FOR COMBINED FITS:\n"
            <<"    "<<setw(8)<<" "
                    <<setw(cl)<<"actual"
                    <<setw(cl)<<"estimate" 
                    <<setw(cl)<<"act. sq(N)"
                    <<setw(cl)<<"est. sq(N)" 
                    <<setw(cl+1)<<"est./actual\n";
    for (int i=0; i<nPar; i++) {
        TString name_i = paramNames[i];
        for (int j = i+1; j<nPar; j++) {
            TString name_j = paramNames[j];
            TString cor_name = "cor("+name_i+","+name_j+")";
            double cor_ij = cov[j][i]/Sqrt(cov[i][i]*cov[j][j]);
            double est_cor_ij = est_cov[j][i]/Sqrt(est_cov[i][i]*est_cov[j][j]);
            double ratio = est_cor_ij/cor_ij;
            cout<<"    "<<setw(8)<<cor_name
                        <<setw(cl)<<cor_ij
                        <<setw(cl)<<est_cor_ij
                        <<setw(cl)<<cor_ij*Sqrt(N)
                        <<setw(cl)<<est_cor_ij*Sqrt(N)
                        <<setw(cl)<<ratio<<"\n";
        }
    }
    return asymmetry_func;
}

void setenvvar(TString & var, const TString env) {
    if (getenv(env))
        var = getenv(env);
    else 
        cout<<"Warning: Environmental variable "<<env<<"is not set.\n"
            <<"Using default "<<var<<" instead.\n";
}


/*
TString timstr() {
    time_t t = time(0);   // get time now
    struct tm * now = localtime(&t);
    tm = *localtime(&t);
    char buffer[256];
    strftime(buffer, sizeof(buffer), "-%a-%b-%d-%H:%M:%S-%Y-", tm);
    string time = put_time(&tm, "-%a-%b-%d-%H:%M:%S-%Y-");
    //filename+="report-"+(now->tm_year + 1900)+"-"
    //        +(now->tm_mon+1)+"-"+now->tm_mday+"-"+now->tm_hour+":"+now->tm_min;

    return TString(time);
}
*/


#define TYPES 2
#define MAXRUNS 512
TString type_name[] = {"TYPE 0", "TYPE 1", "TYPE 2", "TYPE 3", "All"};
TString type_upr[] = {"Type_0", "Type_1", "Type_2", "Type_3", "All"};
TString type_lwr[] = {"type_0", "type_1", "type_2", "type_3", "all"};



/// MAIN APPLICATION
int main(int argc, char *argv[])
{
    //TApplication app("Extract Combined A+b Fitter", &argc, argv);
    srand( time(NULL) );    /// set this to make random or repeatable

    UCNAFierzFitter* ucna[TYPES][MAXRUNS];
    //UCNAFierzFitter ucna("ucna", "UCNA", KEbins, KEmin, KEmax, fit_bins, fit_min, fit_max);

    /// Find file paths from environment
    setenvvar(data_dir, "UCNA_DATA_DIR");
    setenvvar(mc_dir, "UCNA_MONTE_CARLO_DIR");
    setenvvar(mc_sys_dir, "UCNA_SYSTEMATICS_DIR");
    setenvvar(plots_dir, "UCNA_PLOTS_DIR");
    /*
    if (getenv("UCNA_DATA_DIR"))
        data_dir = getenv("UCNA_DATA_DIR");
    else 
        cout<<"Warning: Environmental variable UCNA_DATA_DIR is not set.\n"
            <<"Using default "<<data_dir<<" instead.\n";

    if (getenv("UCNA_MONTE_CARLO_DIR"))
        mc_dir = getenv("UCNA_MONTE_CARLO_DIR");
    else 
        cout<<"Warning: Environmental variable UCNA_MONTE_CARLO_DIR is not set.\n"
            <<"Using default "<<mc_dir<<" instead.\n";

    if (getenv("UCNA_SYSTEMATICS_DIR"))
        mc_sys_dir = getenv("UCNA_SYSTEMATICS_DIR");
    else 
        cout<<"Warning: Environmental variable UCNA_SYSTEMATICS_DIR is not set.\n"
            <<"Using default "<<mc_sys_dir<<" instead.\n";
            */

    if (FIT_TYPE == "b" or FIT_TYPE == "bN") {
    }
    else {
        cout<<"Can run b or bN mode right now.\n";
        exit(1);
    }

    if (argc < 2) {
        cout<<"Requires at least one argument.\n";
        cout<<"Usage: "<<argv[0]<<" <filename | file pattern | dataset> [<filename | twiddled-n | file pattern>]\n";
        exit(1);
    }

    TString dataset = argv[1];
    int arg = 2, run = 0;
    while (arg < argc) {
        cout<<"arg "<<arg<<" is "<<argv[arg]<<"\n";
        /// loops through types
        for (int type=0; type<=TYPES; type++) {
            /// Construct a fierz fitter for this run.
            ucna[type][run] = new UCNAFierzFitter(
                "ucna_"+type_lwr[type]+"_run_"+to_string(run), 
                "UCNA "+type_name[type]+" Run "+to_string(run), 
                KEbins, KEmin, KEmax, fit_bins, fit_min, fit_max);

            /// Load Monte Carlo files .
            TString file_stem = mc_dir+"/"+mc_filename_stem;
            cout<< "    Loading Monte Carlo files - "<<type_name[type]<<" for run "<<run<<".\n";
            assert(ucna[type][run]);
            cout<<file_stem<<" "<<mc_name<<" "<<type<<" "<<afp_ratio<<" "<<argv[arg]<<".\n";
            ucna[type][run]->fill( // TODO file_stem+"%s-%04d.root",
                file_stem+"vector-"+argv[arg]+".root",
                // file_stem+"axial-"+argv[arg]+".root",
                file_stem+"fierz-"+argv[arg]+".root",
                file_number_start, file_number_stop, 
                mc_name, type, 0.4);
            ucna[type][run]->save(plots_dir+"/monte_carlo_"+type_lwr[type]+"_run_"+to_string(run)+"_");

            /// Load data set files that already contain super sum histograms.
            if (dataset == "2010") {
                cout<<" LOADING REAL EVENTS FROM "<<dataset<<" UCNA DATASET - "
                    <<type_name[type]<<" - RUN "<<run<<":\n";
                ucna[type][run]->data.super_sum.fill(
                    data_dir+"/OctetAsym_Offic_2010_FINAL/"+data_filename,
                    "SuperSum_"+type_upr[type],
                    dataset+" final official supersum - "+type_name[type]);
                ucna[type][run]->data.super_sum.save(
                    plots_dir+"/"+dataset+"_data_supersum_"+type_lwr[type]+"_run_"+to_string(run)+".txt");
                // TODO ucna[type][run]->back.super_sum.fill(
                //    ...
                //    dataset+" final official supersum background");
            }
            else {
                /// TODO construct fake from Monte Carlo data
                cout<<"Error: check this code before trusting the results.\n";
                exit(1);

                /// Load Monte Carlo simulated Fierz as fake events.
                /// Pick fake values of asymmetry and Fierz terms.
                //double A = -0.12;
                //double b = 0.00;
                //double N = 1;

                /// Generate fake signal curve from different simulated spectra.
                //fake.compute_fit(A,b,N);
                //ucna[type][run]->data = fake.fit;
            }
            ucna[type][run]->comupte_sizes();
            ucna[type][run]->print_sizes();
        }
        run++;
        arg++;
    }
<<<<<<< Updated upstream

    if(argc <3) { /// TODO change to checking run
        cout << "Nothing to compare to... Done.\n";
        exit(0);    /// TODO temp hack
    }
=======
};

struct UCNAModel {
    int          bins;
    TH1F*        raw[2][2];
    STLhistogram counts[2][2];
    STLhistogram super_ratio;
    STLhistogram super_sum;
    STLhistogram asymmetry;

    void init(int bins) {
        this->bins = bins;
        // TODO init TH1Fs
        for (int i=0; i<2; i++)
            for (int j=0; j<2; j++)
                counts[i][j].init(bins);
        super_ratio.init(bins);
        super_sum.init(bins);
        asymmetry.init(bins);
   }
};

struct UCNAEvent {
    double EdepQ;
    double Edep;
    double MWPCEnergy;
    double ScintPos;
    double MWPCPos;
    double time;
    double primMomentum;
};

// data need to be globals to be visible by func 
/*
vector<double> asymmetry_energy;        
vector<double> asymmetry_values;        
vector<double> asymmetry_errors;        

vector<double> super_ratio_energy;        
vector<double> super_ratio_values;        
vector<double> super_ratio_errors;        

vector<double> super_sum_energy;        
vector<double> super_sum_values;        
vector<double> super_sum_errors;        
*/

UCNAModel ucna_data; // Need construction.
UCNAModel ucna_sm_mc; // Need construction.
UCNAModel ucna_fierz_mc; // Need construction.

void combined_chi2(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *p, Int_t /*iflag */  )
{
	double chi2 = 0; 
	double chi,	E; 

	int n = ucna_data.asymmetry.energy.size();
	//int n = asymmetry_energy.size();
	for (int i = 0; i < n; i++)
	{
		double par[2] = {p[0],p[1]}; // A, b
		E = ucna_data.asymmetry.energy[i];
		//chi = (asymmetry_values[i] - asymmetry_fit_func(&E,par)) / asymmetry_errors[i];
		double Y = ucna_data.asymmetry.values[i];
        double f = asymmetry_fit_func(&E,par);
        double eY = ucna_data.asymmetry.errors[i];
		chi = (Y - f)/eY;
		chi2 += chi*chi; 
	}

	//n = fierzratio_energy.size();
	for (int i = 0; i < n; i++) { 
		//double par[2] = {p[1], expected[0][1]};
		E = ucna_data.super_sum.energy[i];
		//chi = (fierzratio_values[i] - fierzratio_fit_func(&E,par)) / fierzratio_errors[i];
		double Y =      ucna_data    .super_sum.values[i];
        double f = p[1]*ucna_sm_mc   .super_sum.values[i] 
                 + p[2]*ucna_fierz_mc.super_sum.values[i];
        double eY =     ucna_data    .super_sum.errors[i];
		chi = (Y - f)/eY;
		chi2 += chi*chi; 
	}
	fval = chi2; 
}



#if 1
static const int nPar = 3;

TF1* combined_fit(TH1F* asymmetry, TH1F* super_sum, double cov[nPar][nPar]) 
{ 
	/// set up free fit parameters with best guess
	double iniParams[nPar] = {-0.12, 0, 1e6};
	const char * iniParamNames[nPar] = {"A", "b", "N"};

	/// create fit function
	TF1 * func = new TF1("func", asymmetry_fit_func, min_E, max_E, nPar);
	func->SetParameters(iniParams);
	for (int i=0; i<nPar; i++)
		func->SetParName(i, iniParamNames[i]);

	/// fill data structure for fit (coordinates + values + errors) 
	std::cout << "Do global fit" << std::endl;
    	ucna_data.asymmetry.fill(asymmetry);

	/// set up the minuit fitter
	TVirtualFitter::SetDefaultFitter("Minuit");
	TVirtualFitter *minuit = TVirtualFitter::Fitter(0,nPar);
	for (int i=0; i<nPar; ++i)
		minuit->SetParameter(i, func->GetParName(i), func->GetParameter(i), 1, 0, 0);
	minuit->SetFCN(combined_chi2);
	minuit->SetErrorDef(1);	        /// 1 for chi^2

	/// set print level
	double arglist[100];
	arglist[0] = 0;
	minuit->ExecuteCommand("SET PRINT", arglist, 1);

	/// minimize
	arglist[0] = 50;    /// number of function calls
	arglist[1] = 0.1;   /// tolerance
	minuit->ExecuteCommand("MIGRAD",arglist,nPar);

	/// extract results from minuit
	double minParams[nPar];
	double parErrors[nPar];
	for (int i=0; i<nPar; ++i) {  
		minParams[i] = minuit->GetParameter(i);
		parErrors[i] = minuit->GetParError(i);
	}
	double chi2, edm, errdef; 
	int nvpar, nparx;
	minuit->GetStats(chi2,edm,errdef,nvpar,nparx);

	func->SetParameters(minParams);
	func->SetParErrors(parErrors);
	func->SetChisquare(chi2);
	int ndf = ucna_data.asymmetry.energy.size() 
            	+ ucna_data.super_sum.energy.size() - nvpar;
	func->SetNDF(ndf);
    
	TMatrixD matrix( nPar, nPar, minuit->GetCovarianceMatrix() );
	for (int i=0; i<nPar; i++)
		for (int j=0; j<nPar; j++)
			cov[i][j] = minuit->GetCovarianceMatrixElement(i,j);
	
	cout << endl;
	cout << "    chi^2 = " << chi2 << ", ndf = " << ndf << ", chi^2/ndf = " << chi2/ndf << endl;
	cout << endl;

	return func; 
}
#endif
>>>>>>> Stashed changes

    /// TODO loop
    int type = 0;
    run = 0;

<<<<<<< Updated upstream
    /// Take a look at the state of the fits
    ucna[type][run]->data.super_sum.snapshot();
    ucna[type][run]->vector.super_sum.snapshot();
    //ucna[type][run]->axial.super_sum.snapshot();
    ucna[type][run]->fierz.super_sum.snapshot();

    /// fitting the data using the Monte Carlo spectra
    TF1* fit_func = fit(ucna[type][run]);
    if (not fit_func) {
        cout<<" Error: No fitting function returned.\n";
        exit(1);
    }
    cout<<" Successfully fit the Monte Carlo to the data.\n";
    ucna[type][run]->fit.super_sum.snapshot();
    ucna[type][run]->display(plots_dir+"/fit-plot-");

    /// extract most critical results from fit functions and report
    ofstream rofs;
    // time_t t = time(0);   // get time now
    // struct tm * now = localtime(&t);
    //tm = *localtime(&t);
    // char buffer[256];
    // strftime(buffer, sizeof(buffer), "-%a-%b-%d-%H:%M:%S-%Y-", tm);
    //string time = put_time(&tm, "-%a-%b-%d-%H:%M:%S-%Y-");
    //filename+="report-"+(now->tm_year + 1900)+"-"
    //        +(now->tm_mon+1)+"-"+now->tm_mday+"-"+now->tm_hour+":"+now->tm_min;
    string outfilename(plots_dir+"/report.dat");
    rofs.open(outfilename, std::fstream::app | std::fstream::out);
    cout << " FINAL REPORT\n";
    cout<<setw(cl)<<"chi^2"<<setw(cl)<<"ndf"<<setw(cl)<<"chi^2/ndf"<<setw(cl)<<"p(chi2,ndf)";
    for (int i=0; i<nPar; ++i) {  
        TString par_name = fit_func->GetParName(i);
        cout<<setw(cl)<<par_name
            <<setw(cl-5)<<par_name<<" err."
            <<setw(cl-3)<<par_name<<" ex"
            <<setw(cl-8)<<par_name<<" ex err.";
    }
    cout<<"\n";

    double chi2 = fit_func->GetChisquare(); /// Get the chi squared for this fit
    double ndf = fit_func->GetNDF();        /// Get the degrees of freedom in fit
    double prob = Prob(chi2,ndf);           /// Get the probability of this chi2 with this ndf
    
    cout<<setw(cl)<<chi2<<setw(cl)<<ndf
        <<setw(cl)<<chi2/ndf<<setw(cl)<<prob;
    rofs<<setw(cl)<<chi2<<setw(cl)<<ndf
        <<setw(cl)<<chi2/ndf<<setw(cl)<<prob;

    for (int i=0; i<nPar; ++i) {  
        double param = fit_func->GetParameter(i);
        double error = fit_func->GetParError(i);
        cout<<setw(cl)<<param<<setw(cl)<<error
            <<setw(cl)<<param*prob<<setw(cl)<<error*prob;
        rofs<<setw(cl)<<param<<setw(cl)<<error
            <<setw(cl)<<param*prob<<setw(cl)<<error*prob;
=======
int fill_simulation(string filename, string title, THD1* histogram, THD1* super_sum)
{
    TChain *chain = (TChain*)fierz_mc_tfile->Get("SimAnalyzed");
    if (not chain) {
        printf("Error getting Fierz Monte Carlo.\n");
        exit(1);
>>>>>>> Stashed changes
    }
    cout<<"\n";
    rofs<<"\n";
    rofs.close();

<<<<<<< Updated upstream
    //app.Run();
    return 0;
}
=======
	int nevents = chain->GetEntries();
	std::cout << "Total number of Monte Carlo entries without cuts: " << nevents << std::endl;
	chain->SetBranchStatus("*",false);
	chain->SetBranchStatus("PID",true);
	chain->SetBranchStatus("side",true);
	chain->SetBranchStatus("type",true);
	chain->SetBranchStatus("Erecon",true);
	chain->SetBranchStatus("primMomentum",true);
    chain->SetBranchStatus("ScintPosAdjusted",true);
    chain->SetBranchStatus("ScintPosAdjusted",true);

	/*
	TFile* mc_tfile = new TFile("Fierz/mc.root", "recreate");
	if (mc_tfile->IsZombie())
	{
		std::cout << "Can't make MC file.\n";
		exit(1);
	}
	TNtuple* tntuple = new TNtuple("mc_ntuple", "MC NTuple", "s:load:energy");
    */

	unsigned int nSimmed = 0;	/// events have been simulated after cuts
    int PID, side, type;
    double energy;
    double mwpcPosW[3], mwpcPosE[3], primMomentum[3];

    chain->SetBranchAddress("PID",&PID);
    chain->SetBranchAddress("side",&side);
    chain->SetBranchAddress("type",&type);
    chain->SetBranchAddress("Erecon",&energy);
	chain->SetBranchAddress("primMomentum",primMomentum);
    chain->GetBranch("ScintPosAdjusted")->GetLeaf("ScintPosAdjE")->SetAddress(mwpcPosE);
    chain->GetBranch("ScintPosAdjusted")->GetLeaf("ScintPosAdjW")->SetAddress(mwpcPosW);

    for (int evt=0; evt<nevents; evt++) {
        chain->GetEvent(evt);

        /// cut out bad events
        if (PID!=1) 
            continue;

        /// radial fiducial cut
        double radiusW = mwpcPosW[0]*mwpcPosW[0] + mwpcPosW[1]*mwpcPosW[1]; 
        double radiusE = mwpcPosE[0]*mwpcPosE[0] + mwpcPosE[1]*mwpcPosE[1]; 
        if (radiusW > fidcut2 or radiusE > fidcut2) 
            continue;

        /// Type 0, Type I, Type II/III events 
        if (type<4) { 
            /// fill with loading efficiency 
			double p = (large_prime*nSimmed % rand_bins)/rand_bins;
			double afp = (p < afp_off_prob)? -1 : +1;
            //bool spin = (afp * primMomentum[2]) < 0? EAST : WEST;
            bool spin = (afp < 0)? EAST : WEST;
            histogram[side][spin]->Fill(energy, 1);
			tntuple->Fill(side, spin, energy);
			nSimmed++;
        }
        /*  hEreconALL->Fill(Erecon);
        if (type==0) 
            hErecon0->Fill(Erecon);
        else if (type==1) 
            hErecon1->Fill(Erecon);
        else if (type==2 or type==3) 
            hErecon23->Fill(Erecon); */

		/// break when enough data has been generated.
		if(nSimmed >= nToSim)
			break;
    }    
     
    #if 0
	while (false) {
		/// perform energy calibrations/simulations to fill class 
        /// variables with correct values for this simulated event
		//G2P.recalibrateEnergy();
		
		/// check the event characteristics on each side
		for(Side s=EAST; s<=WEST; ++s) {
			/// get event classification type. 
            /// TYPE_IV_EVENT means the event didn't trigger this side.
			/// TODO EventType tp = G2P.fType;

			/// skip non-triggering events, or those outside 50mm position 
            /// cut (you could add your own custom cuts here, if you cared)
            /// TODO place all cuts here!
			//if (tp >= TYPE_I_EVENT or !G2P.passesPositionCut(s) or G2P.fSide != s)
			//	continue;
			
			/// print out event info, (simulated) reconstructed true energy
            /// and position, comparable to values in data
			#ifdef DEBUG
			printf("Event on side %c: type=%i, Erecon=%g @ position (%g,%g), %d\n",
				   sideNames(s), tp, G2P.getErecon(), G2P.wires[s][X_DIRECTION].center, 
				   G2P.wires[s][Y_DIRECTION].center, (unsigned)G2P.getAFP());

			/// print out event primary info, only available in simulation
			printf("\tprimary KE=%g, cos(theta)=%g\n", G2P.ePrim, G2P.costheta);
			#endif 

			/// fill with loading efficiency 
			/// NOT the way to do this anymore 
            //bool load = (nSimmed % 100 < loading_prob);

			/// calculate the energy with a distortion factor
			/// double energy = scale_x * G2P.getErecon();
            double energy = 0;  /// TODO XXXXXXXX
            histogram[s][load]->Fill(energy, 1);
			tntuple->Fill(s, load, energy);
			nSimmed++;
		}
		
		/// break when enough data has been generated.
		if(nSimmed >= nToSim)
			break;
	}
    #endif
    
	std::cout << "Total number of Monte Carlo entries with cuts: " << nSimmed << std::endl;

	//tntuple->SetDirectory(mc_tfile);
	//tntuple->Write();

	/// compute and normalize super sum
    super_sum = compute_super_sum(histogram);
    normalize(histogram, min_E, max_E);

    for (int side=0; side<2; side++)
        for (int spin=0; spin<2; spin++)
            normalize(histogram[side][spin], min_E, max_E);
}


int main(int argc, char *argv[])
{
	TApplication app("Extract Combined A + b Fitter", &argc, argv);
	//TH1::AddDirectory(kFALSE);
	/// Geant4 MC data scanner object
	/// G4toPMT G2P;
	/// use data from these MC files (most recent unpolarized beta decay, 
    /// includes Fermi function spectrum correction)
	/// note wildcard * in filename; 
    /// MC output is split up over many files, but G2P will TChain them together
	//G2P.addFile("/data2/mmendenhall/G4Out/2010/20120823_neutronBetaUnpol/analyzed_*.root");
	
	/// PMT Calibrator loads run-specific energy calibrations info for selected run
	/// and uses default Calibrations DB connection to most up-to-date though possibly unstable "mpm_debug"
	
	/// If you really want this to be random, 
    /// you will need to seed rand() with something 
    /// other than default.
	srand( time(NULL) );

	/// load the files that contain data histograms
    TFile *asymmetry_data_tfile = new TFile(
        "/media/hickerson/boson/Data/OctetAsym_Offic_2010_FINAL/"
        "Range_0-1000/CorrectAsym/CorrectedAsym.root");
	if (asymmetry_data_tfile->IsZombie())
	{
		std::cout << "File not found." << std::endl;
		exit(1);
	}

    TFile *ucna_data_tfile = new TFile(
        "/media/hickerson/boson/Data/OctetAsym_Offic_2010_FINAL/"
		"OctetAsym_Offic.root");
	if (ucna_data_tfile->IsZombie())
	{
		std::cout << "File not found." << std::endl;
		exit(1);
	}

	/// extract the histograms from the files
    TH1F *asymmetry_histogram = 
            (TH1F*)asymmetry_data_tfile->Get("hAsym_Corrected_C");
    if (not asymmetry_histogram) {
        printf("Error getting asymmetry histogram.\n");
        exit(1);
    }

    TH1F *super_sum_histogram = 
            (TH1F*)ucna_data_tfile->Get("Total_Events_SuperSum");
    if (not super_sum_histogram) {
        printf("Error getting super sum histogram.\n");
        exit(1);
    }

    fill_simulation("/home/xuansun/Documents/SimData_Beta/"
                    "SimAnalyzed_Beta.root",
                    "Monte Carlo Standard Model beta spectrum",
                    mc.sm_histogram,
					mc.sm_super_sum_hisotgram);

    fill_simulation("/home/xuansun/Documents/SimData_Beta/"
                    "SimAnalyzed_Beta_fierz.root",
                    "Monte Carlo Fierz beta spectrum",
                    mc.fierz_histogram,
					mc.fierz_super_sum_hisotgram);

	//histogram->SetDirectory(mc_tfile);
	//histogram->Write();

	//mc.sm_super_sum_histogram->SetDirectory(mc_tfile);
	//mc.sm_super_sum_histogram->Write();

	//mc.fierz_super_sum_histogram->SetDirectory(mc_tfile);
	//mc.fierz_super_sum_histogram->Write();

	//mc_tfile->Close();

/*
    for (int side=0; side<2; side++)
        for (int spin=0; spin<2; spin++)
		{
            normalize(mc.fierz_histogram[side][spin], min_E, max_E);
            normalize(mc.sm_histogram[side][spin], min_E, max_E);
        }*/

    TCanvas *canvas = new TCanvas("fierz_canvas", "Fierz component of energy spectrum");
>>>>>>> Stashed changes


    /*
    /// LOAD PRECOMPUTED HISTOGRAMS AND OVERWRITE 
    /// Load the files that already contain data asymmetry histogram.
    ucna.data.asymmetry.fill(
        data_dir_2010+"Range_0-1000/CorrectAsym/CorrectedAsym.root",
        "hAsym_Corrected_C",
        "2010 final official asymmetry");
    /// Load the files that already contain data super histogram.
    ucna.data.super_sum.fill(
        data_dir_2010+"OctetAsym_Offic.root",
        "Total_Events_SuperSum",
        "2010 final official supersum");
    */
    /*
    /// Load the files that already contain data super histogram.
    bool use_real_data = false;
    if (use_real_data) {
        for (int side=EAST; side<=WEST; side++) {
            for (int afp=EAST; afp<=WEST; afp++) {
                TString s = side? "W":"E", a = afp? "On":"Off";
                TString title = "2010 final official "+s+" afp "+a;
                TString name = "hTotalEvents_"+s+"_"+a+";1";
                int entries = real.data.counts[side][afp]->fill(data_dir_2010+"OctetAsym_Offic.root", name, title);
                if (entries) {
                    cout<<"Status: Number of entries in "<<(side? "west":"east")
                        <<" side with afp "<<a<<" is "<<entries<<".\n";
                }
                else
                    cout<<"Error: found no events in "<<title<<".\n";
            }
        }
    }
    */
<<<<<<< Updated upstream
=======

    TString fit_pdf_filename = "mc/fierz_fit_data.pdf";
    canvas->SaveAs(fit_pdf_filename);

    // compute and plot the super ratio
    TH1F *super_ratio_histogram = compute_super_ratio(ucna_data.raw);
    super_ratio_histogram->Draw();
    TString super_ratio_pdf_filename = "mc/super_ratio_data.pdf";
    canvas->SaveAs(super_ratio_pdf_filename);

    // compute and plot the super ratio asymmetry 
    //TH1F *asymmetry_histogram = compute_corrected_asymmetry(ucna_data.raw, ucna_correction_histogram);

	// fit the Fierz term from the asymmetry
	/*
	char A_fit_str[1024];
    sprintf(A_fit_str, "[0]/(1+[1]*(%f/(%f+x)))", m_e, m_e);
    TF1 *A_fierz_fit = new TF1("A_fierz_fit", A_fit_str, min_E, max_E);
    A_fierz_fit->SetParameter(-.12,0);
    A_fierz_fit->SetParameter(1,0);
	asymmetry_histogram->Fit(A_fierz_fit, "Sr");
	asymmetry_histogram->SetMaximum(0);
	asymmetry_histogram->SetMinimum(-0.2);
	*/

	// compute chi squared
	/*
    double chisq = A_fierz_fit->GetChisquare();
    double NDF = A_fierz_fit->GetNDF();
	char A_str[1024];
	sprintf(A_str, "A = %1.3f #pm %1.3f", A_fierz_fit->GetParameter(0), A_fierz_fit->GetParError(0));
	char A_b_str[1024];
	sprintf(A_b_str, "b = %1.3f #pm %1.3f", A_fierz_fit->GetParameter(1), A_fierz_fit->GetParError(1));
	char A_b_chisq_str[1024];
    printf("Chi^2 / ( NDF - 1) = %f / %f = %f\n", chisq, NDF-1, chisq/(NDF-1));
	sprintf(A_b_chisq_str, "#frac{#chi^{2}}{n-1} = %f", chisq/(NDF-1));
	*/

	// draw the ratio plot
	//asymmetry_histogram->SetStats(0);
    //asymmetry_histogram->Draw();

	// draw a legend on the plot
	/*
    TLegend* asym_legend = new TLegend(0.3,0.85,0.6,0.65);
    asym_legend->AddEntry(asymmetry_histogram, "Asymmetry data", "l");
    asym_legend->AddEntry(A_fierz_fit, "Fierz term fit", "l");
    asym_legend->AddEntry((TObject*)0, A_str, "");
    asym_legend->AddEntry((TObject*)0, A_b_str, "");
    asym_legend->AddEntry((TObject*)0, A_b_chisq_str, "");
    asym_legend->SetTextSize(0.03);
    asym_legend->SetBorderSize(0);
    asym_legend->Draw();
	*/

    TString asymmetry_pdf_filename = "mc/asymmetry_data.pdf";
    canvas->SaveAs(asymmetry_pdf_filename);

    /// Compute the super sums
    TH1F *super_sum_histogram = compute_super_sum(ucna_data.raw);
	std::cout << "Number of super sum entries " << (int)super_sum_histogram->GetEntries() << std::endl;
    normalize(super_sum_histogram, min_E, max_E);
    super_sum_histogram->SetLineColor(2);
	super_sum_histogram->SetStats(0);
    super_sum_histogram->Draw("");

    /// Draw Monte Carlo
    mc.sm_super_sum_histogram->SetLineColor(1);
    mc.sm_super_sum_histogram->SetMarkerStyle(4);
    mc.sm_super_sum_histogram->Draw("same p0");

    /// make a pretty legend
    TLegend * legend = new TLegend(0.6,0.8,0.7,0.6);
    legend->AddEntry(super_sum_histogram, "Type 0 super_sum", "l");
    legend->AddEntry(mc.sm_super_sum_histogram, "Monte Carlo super_sum", "p");
    //legend->AddEntry(bonehead_sum_histogram, "Bonehead sum", "l");
    legend->SetTextSize(0.03);
    legend->SetBorderSize(0);
    legend->Draw();

    /// save the data and Mote Carlo plots
    TString super_sum_pdf_filename = "mc/super_sum_data.pdf";
    canvas->SaveAs(super_sum_pdf_filename);
	TFile* ratio_tfile = new TFile("Fierz/ratio.root", "recreate");
	if (ratio_tfile->IsZombie())
	{
		std::cout << "Can't recreate MC file" << std::endl;
		exit(1);
	}

    /// compute little b factor
    TH1F *fierz_ratio_histogram = new TH1F(*super_sum_histogram);
	fierz_ratio_histogram->SetName("fierz_ratio_histogram");
    //fierz_ratio_histogram->Divide(super_sum_histogram, mc.sm_super_sum_histogram);
    int bins = super_sum_histogram->GetNbinsX();
    for (int bin = 1; bin < bins+1; bin++) {
		double X = super_sum_histogram->GetBinContent(bin);
		double Y = mc.sm_super_sum_histogram->GetBinContent(bin);
		double Z = Y > 0 ? X/Y : 0;

		fierz_ratio_histogram->SetBinContent(bin, Z);

		double x = super_sum_histogram->GetBinError(bin);
		double y = mc.sm_super_sum_histogram->GetBinError(bin);
		fierz_ratio_histogram->SetBinError(bin, Z*TMath::Sqrt(x*x/X/X + y*y/Y/Y));
	}
    fierz_ratio_histogram->GetYaxis()->SetRangeUser(0.9,1.1); // Set the range
    fierz_ratio_histogram->SetTitle("Ratio of UCNA data to Monte Carlo");
	std::cout << super_sum_histogram->GetNbinsX() << std::endl;
	std::cout << mc.sm_super_sum_histogram->GetNbinsX() << std::endl;
	

	/// fit the Fierz ratio 
	char fit_str[1024];
    sprintf(fit_str, "1+[0]*(%f/(%f+x)-%f)", m_e, m_e, expected_fierz);
    TF1 *fierz_fit = new TF1("fierz_fit", fit_str, min_E, max_E);
    fierz_fit->SetParameter(0,0);
	fierz_ratio_histogram->Fit(fierz_fit, "Sr");

	/// A fit histogram for output to gnuplot
    TH1F *fierz_fit_histogram = new TH1F(*super_sum_histogram);
	for (int i = 0; i < fierz_fit_histogram->GetNbinsX(); i++)
		fierz_fit_histogram->SetBinContent(i, fierz_fit->Eval(fierz_fit_histogram->GetBinCenter(i)));

	/// compute chi squared
    double chisq = fierz_fit->GetChisquare();
    double NDF = fierz_fit->GetNDF();
	char b_str[1024];
	sprintf(b_str, "b = %1.3f #pm %1.3f", fierz_fit->GetParameter(0), fierz_fit->GetParError(0));
	char chisq_str[1024];
    printf("Chi^2 / (NDF-1) = %f / %f = %f\n", chisq, NDF-1, chisq/(NDF-1));
	sprintf(chisq_str, "#frac{#chi^{2}}{n-1} = %f", chisq/(NDF-1));

	/// draw the ratio plot
	fierz_ratio_histogram->SetStats(0);
    fierz_ratio_histogram->Draw();

	/// draw a legend on the plot
    TLegend* ratio_legend = new TLegend(0.3,0.85,0.6,0.65);
    ratio_legend->AddEntry(fierz_ratio_histogram, "Data ratio to Monte Carlo (Type 0)", "l");
    ratio_legend->AddEntry(fierz_fit, "Fierz term fit", "l");
    ratio_legend->AddEntry((TObject*)0, "1+b(#frac{m_{e}}{E} - #LT #frac{m_{e}}{E} #GT)", "");
    ratio_legend->AddEntry((TObject*)0, b_str, "");
    ratio_legend->AddEntry((TObject*)0, chisq_str, "");
    ratio_legend->SetTextSize(0.03);
    ratio_legend->SetBorderSize(0);
    ratio_legend->Draw();

	/// output for root
    TString fierz_ratio_pdf_filename = "mc/fierz_ratio.pdf";
    canvas->SaveAs(fierz_ratio_pdf_filename);

	/// output for gnuplot
	output_histogram("mc/super-sum-data.dat", super_sum_histogram, 1, 1000);
	output_histogram("mc/super-sum-mc.dat", mc.sm_super_sum_histogram, 1, 1000);
	output_histogram("mc/fierz-ratio.dat", fierz_ratio_histogram, 1, 1);
	output_histogram("mc/fierz-fit.dat", fierz_fit_histogram, 1, 1);


	fierz_ratio_histogram->SetDirectory(ratio_tfile);
	fierz_ratio_histogram->Write();
	//ratio_tfile->Close();


    /// CODE BREAK


	double cov[nPar][nPar]; 
	double entries = super_sum_histogram->GetEffectiveEntries();
	double N = GetEntries(super_sum_histogram, min_E, max_E);

	/// set all expectation values for this range
	for (int m=0; m<=2; m++)
		for (int n=0; n<=2; n++)
			expected[m][n] = evaluate_expected_fierz(m, n, min_E, max_E);
	
	/// find the predicted inverse covariance matrix for this range
	double A = -0.12;
	TMatrixD p_cov_inv(2,2);
	p_cov_inv[0][0] =  N/4*expected[2][0];
	p_cov_inv[1][0] = 
    p_cov_inv[0][1] = -N*A/4*expected[2][1];
	p_cov_inv[1][1] =  N*(expected[0][2] - expected[0][1]*expected[0][1]);

	/// find the covariance matrix
	double det = 0;
	double p_var_A = 1 / p_cov_inv[0][0];
	double p_var_b = 1 / p_cov_inv[1][1];
	TMatrixD p_cov = p_cov_inv.Invert(&det);

	/// do the fitting
	TF1* func = combined_fit(asymmetry_histogram, super_sum_histogram, cov);

	/// output the data info
	cout << " ENERGY RANGE:" << endl;
	cout << "    Energy range is " << min_E << " - " << max_E << " keV" << endl;
	cout << "    Number of counts in full data is " << (int)entries << endl;
	cout << "    Number of counts in energy range is " <<  (int)N << endl;
	cout << endl;

	/// output the details	
	cout << " FIT COVARIANCE MATRIX cov(A,b) =\n";
	for (int i = 0; i < nPar; i++) {
		for (int j = 0; j < nPar; j++) {
			cout << "\t" << cov[i][j];
		}
		cout << "\n";
	}

	double sig_A = sqrt(cov[0][0]);
	double sig_b = sqrt(cov[1][1]);

	cout << endl;
	cout << " PREDICTED COVARIANCE MATRIX cov(A,b) =\n";
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			cout << "\t" << p_cov[i][j];
		}
		cout << "\n";
	}

	double p_sig_A = sqrt(p_cov[0][0]);
	double p_sig_b = sqrt(p_cov[1][1]);

	cout << endl;
	cout << " FOR UNCOMBINED FITS:\n";
	cout << "    Expected statistical error for A in this range is without b is " 
					<< sqrt(p_var_A) << endl;
	cout << "    Expected statistical error for b in this range is without A is "
					<< sqrt(p_var_b) << endl;
	cout << endl;
	cout << " FOR COMBINED FITS:\n";
	cout << "    Expected statistical error for A in this range is " << p_sig_A << endl;
	cout << "    Actual statistical error for A in this range is " << sig_A << endl;
	cout << "    Ratio for A error is " << sig_A / p_sig_A << endl;
	cout << "    Expected statistical error for b in this range is " << p_sig_b << endl;
	cout << "    Actual statistical error for b in this range is " << sig_b << endl;
	cout << "    Ratio for b error is " << sig_b / p_sig_b << endl;
	cout << "    Expected cor(A,b) = " << p_cov[1][0] / p_sig_A / p_sig_b << endl;
	cout << "    Actual cor(A,b) = " << cov[1][0] / sqrt(cov[0][0] * cov[1][1]) << endl;

	/*
	// A fit histogram for output to gnuplot
    TH1F *fierz_fit_histogram = new TH1F(*asymmetry_histogram);
	for (int i = 0; i < fierz_fit_histogram->GetNbinsX(); i++)
		fierz_fit_histogram->SetBinContent(i, fierz_fit->Eval(fierz_fit_histogram->GetBinCenter(i)));

	

	/// output for root
    TString pdf_filename = "/data/kevinh/mc/asymmetry_fierz_term_fit.pdf";
    canvas->SaveAs(pdf_filename);

	/// output for gnuplot
	//output_histogram("/data/kevinh/mc/super-sum-data.dat", super_sum_histogram, 1, 1000);
	//output_histogram("/data/kevinh/mc/super-sum-mc.dat", mc.sm_super_sum_histogram, 1, 1000);
	//output_histogram("/data/kevinh/mc/fierz-ratio.dat", fierz_ratio_histogram, 1, 1);
	//output_histogram("/data/kevinh/mc/fierz-fit.dat", fierz_fit_histogram, 1, 1);

	*/

	// Create a new canvas.
	TCanvas * c1 = new TCanvas("c1","Two HIstogram Fit example",100,10,900,800);
	c1->Divide(2,2);
	gStyle->SetOptFit();
	gStyle->SetStatY(0.6);

	c1->cd(1);
	asymmetry_histogram->Draw();
	func->SetRange(min_E, max_E);
	func->DrawCopy("cont1 same");
	/*
	c1->cd(2);
	asymmetry->Draw("lego");
	func->DrawCopy("surf1 same");
	*/
	c1->cd(3);
	func->SetRange(min_E, max_E);
	///TODO fierzratio_histogram->Draw();
	func->DrawCopy("cont1 same");
	/*
	c1->cd(4);
	fierzratio->Draw("lego");
	gPad->SetLogz();
	func->Draw("surf1 same");
	*/

	app.Run();

	return 0;
}
>>>>>>> Stashed changes
