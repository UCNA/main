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

    if(argc <3) { /// TODO change to checking run
        cout << "Nothing to compare to... Done.\n";
        exit(0);    /// TODO temp hack
    }

    /// TODO loop
    int type = 0;
    run = 0;

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
    }
    cout<<"\n";
    rofs<<"\n";
    rofs.close();

    //app.Run();
    return 0;
}


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
