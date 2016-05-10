#include "FierzFitter.hh"


double UCNAhistogram::normalize() {
    return normalize(GetMinimum(), GetMaximum());
}

double UCNAhistogram::normalize(double min, double max)
{
	int bin_min = GetMinimumBin();
	int bin_max = GetMaximumBin();
    double integrand = 0;
    integrand +=(GetBinCenter(bin_min) - min)
               / GetBinWidth(bin_min)
               / GetBinWidth(bin_min)
               * GetBinContent(bin_min)
               +(max - GetBinCenter(bin_max))
               / GetBinWidth(bin_max)
               / GetBinWidth(bin_max)
               * GetBinContent(bin_max);

    for (int bin=bin_min+1; bin<bin_max; bin++)
        integrand += GetBinContent(bin)
                   / GetBinWidth(bin);

	Scale(1/integrand);
    return integrand;
}


//TF1* UCNAFierzFitter::combined_fit(TH1D* asymmetry, TH1D* super_sum, TMatrixD &cov, TF1* func)
TF1* UCNAFierzFitter::combined_fit(
        TMatrixD &cov, TF1 *func, 
        void (*global_fcn_ptr)(Int_t&, Double_t*, Double_t&, Double_t*, Int_t) ) { 
    int nPar = func->GetNpar();
    if (not func) {
        cout<<"Fit function not set.\n";
        exit(1);
    }

	/// set up the Minuit fitter
	TVirtualFitter::SetDefaultFitter("Minuit");
	TVirtualFitter *minuit = TVirtualFitter::Fitter(0,nPar);
	for (int i=0; i<nPar; ++i)
		minuit->SetParameter(i, func->GetParName(i), func->GetParameter(i), 1, 0, 0);
	minuit->SetFCN(global_fcn_ptr);
	minuit->SetErrorDef(1);	    /// 1 for chi^2

	/// set print level
	double arglist[100];
	arglist[0] = 0;
	minuit->ExecuteCommand("SET PRINT", arglist, 1);

	/// minimize
	arglist[0] = 50;            /// number of function calls
	arglist[1] = 0.1;           /// tolerance
	minuit->ExecuteCommand("MIGRAD", arglist, nPar);

	/// extract results from Minuit
	double chi2, edm, errdef; 
	int nvpar, nparx;
	minuit->GetStats(chi2,edm,errdef,nvpar,nparx);
	func->SetChisquare(chi2);

	for (int i=0; i<nPar; ++i) {  
		double param = minuit->GetParameter(i);
		double error = minuit->GetParError(i);
        func->SetParameter(i,param);
        func->SetParError(i,error);
	}

	int ndf = data.asymmetry.GetNbinsX() + data.super_sum.GetNbinsX() - nvpar;
	func->SetNDF(ndf);
    
	TMatrixD matrix( nPar, nPar, minuit->GetCovarianceMatrix() );
	for (int i=0; i<nPar; i++)
		for (int j=0; j<nPar; j++)
			cov[i][j] = minuit->GetCovarianceMatrixElement(i,j);
	
	cout<<"\n    chi^2 = "<<chi2<<", ndf = "<<ndf<<", chi^2/ndf = "<<chi2/ndf<<".\n\n";

	return func; 
}


double UCNAFierzFitter::asymmetry_chi2(double A, double b) {
        return data.asymmetry_chi2(A,b);
}

double UCNAmodel::asymmetry_chi2(double A, double b) {
    /*if (not asymmetry) {
        cout<<"Error: Asymmetry histogram is not constructed.\n";
        exit(1);
    }*/

	double chi2 = 0;
    int n = asymmetry.GetNbinsX();
	for (int i = 0; i < n; i++)
	{
		double E = asymmetry.GetBinCenter(i);
        if (E < min or E > max)
            continue;
		double Y = asymmetry.GetBinContent(i);
		double eY = asymmetry.GetBinError(i);
        double f = A/(1 + b*m_e/(E+m_e)); 
        //double f = asymmetry_fit_func(&E,p);
        if (eY > 0) {
            double chi = (Y-f)/eY;
            chi2 += chi*chi; 
        }
	}
    return chi2;
}


double UCNAFierzFitter::supersum_chi2(double b, double N)
{
	double chi2 = 0;
    int n = data.asymmetry.GetNbinsX();
	for (int i = 0; i < n; i++)
	{
		double E = data.asymmetry.GetBinCenter(i);
        if (E < min or E > max)
            continue;
		double Y = data.super_sum.GetBinContent(i);
		double eY = data.super_sum.GetBinError(i);
        double f  = N*sm.super_sum.GetBinContent(i) 
                  + N*b*fierz.super_sum.GetBinContent(i);
        if (eY > 0) {
            double chi = (Y-f)/eY;
            chi2 += chi*chi; 
        }
	}
    return chi2;
}

double UCNAFierzFitter::combined_chi2(double A, double b, double N)
{
	double chi2 = asymmetry_chi2(A,b) + supersum_chi2(b,N);
    return chi2;
}

void UCNAFierzFitter::compute_fit(/*double A,*/ double b, double N)
{
    int n = data.asymmetry.GetNbinsX();
	for (int i = 0; i < n; i++)
	{
		double E = data.asymmetry.GetBinCenter(i);
        if (E < min or E > max)
            continue;
        double pSM = sm.super_sum.GetBinContent(i);
        double pF = fierz.super_sum.GetBinContent(i);
        double f  = N*(pSM + b*pF);
        fit.super_sum.SetBinContent(i,f);

        double eSM = sm.super_sum.GetBinContent(i);
        double eF = fierz.super_sum.GetBinContent(i);
        double ef = N*Sqrt(eSM*eSM + b*b+eF*eF);
        fit.super_sum.SetBinError(i,ef);
	}
}

double asymmetry_fit_func(double *, double *)
{
    //return f = A/(1 + b*m_e/(E+m_e)); 
    return 0;
}

/** 
 * UCNAhistogram::fill
 * load the files that contain data histograms
 */
int UCNAhistogram::fill(TString filename) {
    return fill(filename, GetName(), GetTitle());
}

int UCNAhistogram::fill(TString filename, 
                        TString name, 
                        TString title) {
    //, TH1D* histogram)
	TFile* tfile = new TFile(filename);
	if (tfile->IsZombie()) {
		cout<<"Error loading "<<title<<":\n";
		cout<<"File not found: "<<filename<<".\n";
		return 0;
	}

    /*if (histogram) {
		cout<<"Warning: Histogram "<<title<<" already exists and is being ignored.\n";
        cout<<histogram<<endl;
        delete histogram;
    }*/

    // TODO warn if overwriting data...
    TH1D* histogram = (TH1D*)tfile->Get(name);
    if (not histogram) {
		cout<<"Error: In file "<<filename<<":\n";
		cout<<"       Error getting "<<title<<".\n";
		cout<<"       Cannot find histogram named "<<name<<".\n";
        return 0;
    }
    *(TH1D*)this = *histogram;

	int entries = GetEntries();
	cout<<"Number of entries in "<<title<<" is "<<entries<<".\n";
    return entries;
}


/**
 * UCNAhistogram::save() and UCNAmodel::save()
 * Save root data
 */
void UCNAmodel::save(TString filename)
     //TH1D* counts[2][2], TH1D* super_sum, TH1D* asymmetry)
{
	TFile* tfile = new TFile(filename, "recreate");
	if (tfile->IsZombie()) {
		cout<<"Error: Problem saving "<<title<<":\n";
		cout<<"       Cannot create file "<<filename<<".\n";
        cout<<"Aborting...\n";
		exit(1);
	}

    //if (test_construction(counts, super_sum)) {
        for (int side=0; side<2; side++)
            for (int spin=0; spin<2; spin++) {
                counts[side][spin]->SetDirectory(tfile);
                counts[side][spin]->Write();
            }
    /*} else {
        cout<<"Error: Rates or super sum for "<<name<<":\n";
        cout<<"       Histogram not constructed.\n";
        cout<<"Aborting...\n";
        exit(1);
    }

    //if (super_sum) {*/
        super_sum.SetDirectory(tfile);
        super_sum.Write();
    /*} else {
        cout<<"Error: Super sum for "<<name<<":\n";
        cout<<"       Histogram not constructed.\n";
        cout<<"Aborting...\n";
        exit(1);
    }

    if (asymmetry) {*/
        asymmetry.SetDirectory(tfile);
        asymmetry.Write();
    /*} else {
        cout<<"Error: Asymmetry for "<<name<<":\n";
        cout<<"       Histogram not constructed.\n";
        cout<<"Aborting...\n";
        exit(1);
    }*/

    if (ntuple) {
        ntuple->SetDirectory(tfile);
        ntuple->Write();
    } else {
        cout<<"Warning: Ntuple not set. Can't save data.\n";
    }

	tfile->Close();
}

/**
 * UCNAhistogram::fill() and UCNAmodel::fill()
 * Fill in asymmetry and super_ratio, and super sums from simulation data.
 * Use wild card * in filename where data is split up over many files
 * and they get Tchained together.
 */
int UCNAmodel::fill(TString filename, TString name, TString title)
        //TH1D* counts[2][2], TH1D* super_sum, TH1D* asymmetry)
{
    /*
    if (not test_construction(counts, super_sum)) {
        cout<<"Error: Rates or super sum for "<<name<<":\n";
        cout<<"       Histogram not constructed.\n";
        cout<<"Aborting...\n";
        exit(1);
    }

    if (not asymmetry) {
        cout<<"Error: Asymmetry for "<<name<<":\n";
        cout<<"       Histogram not constructed.\n";
        cout<<"Aborting...\n";
        exit(1);
    }
    */

	TFile* tfile = new TFile(filename);
	if (tfile->IsZombie()) {
		cout<<"Error: Problem filling "<<title<<":\n";
		cout<<"       File "<<filename<<" not found\n";
        cout<<"Aborting...\n";
		exit(1);
	}

    TChain *chain = (TChain*)tfile->Get(name);
    if (not chain) {
		cout<<"Error: In file "<<filename<<":\n";
		cout<<"       Cannot get "<<title<<":\n";
		cout<<"       Cannot find chain or tree named "<<name<<".\n";
        cout<<"Aborting...\n";
        exit(1);
    }

	int nEvents = chain->GetEntries();
	chain->SetBranchStatus("*",false);
	chain->SetBranchStatus("PID",true);
	chain->SetBranchStatus("side",true);
	chain->SetBranchStatus("type",true);
	chain->SetBranchStatus("Erecon",true);
	chain->SetBranchStatus("primMomentum",true);
    chain->SetBranchStatus("ScintPosAdjusted",true);
    chain->SetBranchStatus("ScintPosAdjusted",true);

	//TNtuple* tntuple = new TNtuple("mc_ntuple", "MC NTuple", "s:load:energy"); */
	TNtuple* tntuple = 0;

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

    for (int evt=0; evt<nEvents; evt++) {
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
            double p = rand.Uniform(1);
			double afp = (p < afp_off_prob)? -1 : +1;
            bool spin = (afp < 0)? EAST : WEST;
            //cout<<"energy: "<<energy<<" side: "<<"side: "<<side
            //         <<" spin: "<<spin<<" afp: "<<afp<<" p: "<<p<<".\n";
            counts[side][spin]->Fill(energy, 1);
            if (tntuple)
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
     
	cout<<"Total number of Monte Carlo entries:\n";
	cout<<"      Entries without cuts:  "<<nEvents<<endl;
	cout<<"      Entries with cuts:     "<<nSimmed<<endl;
	cout<<"      Efficiencies of cuts:  "<<(100.0*nSimmed)/double(nEvents)<<"%\n";

	/// compute and normalize super sum
    compute_super_sum();
    compute_asymmetry();

    super_sum.normalize();
    //normalize(asymmetry, min_E, max_E);
    //for (int side=0; side<2; side++)
    //    for (int spin=0; spin<2; spin++)
    //        normalize(counts[side][spin], min_E, max_E);

    return nSimmed;
}


/**
 * test_counts
 * Tests that all histograms are constructed.
    bool test_counts(TH1D* counts[2][2])
 */
bool UCNAmodel::test_counts() {
    for (int side=0; side<2; side++)
        for (int spin=0; spin<2; spin++)
            if (not counts[side][spin]) {
                cout<<"Error: rate histogram on the "
                    <<(side? "west":"east")<<" side with afp "
                    <<(spin? "on":"off")<<" is not constructed.\n";
                return false;
            }
    return true;
}


/**
 * test_minimum
 * Tests the range properties of a histogram. 
 * Is the top and bottom energies aligned with bins, for example.
 * 
 * bin = 0          underflow bin
 * bin = 1          first bin with low-edge INCLUDED
 * bin = bins       last bin with upper-edge EXCLUDED
 * bin = bins + 1   overflow bin
 *
    int test_min(TH1D* histogram, double min = 0) {
        if (not histogram) {
            cout<<"Error: No histogram to test range on.\n";
            return false;
        }
 */
bool UCNAhistogram::test_min() {
    double min = GetMinimum();
    return test_min(min);
}

bool UCNAhistogram::test_min(double min) {
    TAxis *axis = GetXaxis();
    if (not axis) {
        cout<<"Error: No axis in histogram to test range on.\n";
        return false;
    }

    double bin_min = axis->FindBin(min);
    double lower = axis->GetBinLowEdge(bin_min);
    if (min and min != lower) {
        cout<<"Error: Minimum does not align with a bin minimum.\n";
        cout<<"       Minimum is "<<min<<" and bin minimum is "<<lower<<".\n";
        return false;
    }

    return true;
}


/**
 * test_max
 * Tests the max range properties of a histogram. 
 * 
 * bin = 0          underflow bin
 * bin = 1          first bin with low-edge INCLUDED
 * bin = bins       last bin with upper-edge EXCLUDED
 * bin = bins + 1   overflow bin
 *
    int UCNAhistogram::test_max(TH1D* histogram, double max = 0) 
    {
        TAxis *axis = histogram->GetXaxis();
        if (not axis) {
            cout<<"Error: No axis in histogram to test range on.\n";
            return false;
        }
 */
bool UCNAhistogram::test_max() {
    double max = GetMaximum();
    return test_max(max);
}

bool UCNAhistogram::test_max(double max) {
    TAxis *axis = GetXaxis();
    if (not axis) {
        cout<<"Error: No axis in histogram to test range on.\n";
        return false;
    }

    double bin_max = axis->FindBin(max);
    double upper = axis->GetBinUpEdge(bin_max-1); // TODO check this works in all cases
    if (max and max != upper) {
        cout<<"Error: Maximum does not align with a bin maximum.\n";
        cout<<"       Maximum is "<<max<<" and bin maximum is "<<upper<<".\n";
        return false;
    }

    return bin_max;
}


/**
 * test_range
 * Tests the rang properties of a histogram. 
 * Is the top and bottom energies aligned with bins, for example.
 * 
 * bin = 0          underflow bin
 * bin = 1          first bin with low-edge INCLUDED
 * bin = bins       last bin with upper-edge EXCLUDED
 * bin = bins+1     overflow bin
 * 
    bool test_range(TH1D* histogram, double min = 0, double max = 0) 
    {
        if (not histogram) {
            cout<<"Error: No histogram to test ranges on.\n";
            return false;
        }
 */
bool UCNAhistogram::test_range() {
    double min = GetMinimum();
    double max = GetMaximum();
    return test_range(min, max);
}

bool UCNAhistogram::test_range(double min, double max) {
    TAxis *axis = GetXaxis();
    if (not axis) {
        cout<<"Error: No axis in histogram to test ranges on.\n";
        return false;
    }

    double bin_min = axis->FindBin(min);
    double lower = axis->GetBinLowEdge(bin_min);
    if (min and min != lower) {
        cout<<"Error: Minimum does not align with a bin minimum.\n";
        cout<<"       Minimum is "<<min<<" and bin minimum is "<<lower<<".\n";
        return false;
    }

    double bin_max = axis->FindBin(max);
    double upper = axis->GetBinUpEdge(bin_max-1);
    if (max and max != upper) {
        cout<<"Error: Maximum does not align with a bin maximum.\n";
        cout<<"       Maximum is "<<max<<" and bin maximum is "<<upper<<".\n";
        return false;
    }

    if (lower < 0) {
        cout<<"Error: Minimum is negative.\n";
        return false;
    }

    if (upper <= lower) {
        cout<<"Error: Maximum is not greater than minimum.\n";
        return false;
    }

    return true;
}


/**
 * test_construction
 * Tests the rang properties of a histogram. 
 * Is the top and bottom energies aligned with bins, for example.
 *
    bool test_construction(TH1D* counts[2][2]) 
 */
bool UCNAmodel::test_construction() {
    if (not test_counts())
        return false;

    for (int side=0; side < 2; side++)
        for (int spin = 0; spin < 2; spin++) {
            int rate_bins = counts[side][spin]->GetNbinsX();
            if (bins != rate_bins) {
                cout<<"Error: Extracted and side spin histogram sizes don't match.\n";
                cout<<"       Rate histogram on the "
                    <<(side? "west":"east")<<" side with afp "
                    <<(spin? "on":"off")<<" has "<<rate_bins<<".\n";
                cout<<"       Extracted histogram has "<<bins<<" bins.\n";
                return false;
            }
            if (rate_bins <= 0) {
                cout<<"Error: Bad bin number.\n";
                cout<<"       Rate histogram on the "
                    <<(side? "west":"east")<<" side with afp "
                    <<(spin? "on":"off")<<" has "<<rate_bins<<".\n";
                return false;
            }
        }

    return true;
}


/// Gets the counts in bin
void UCNAmodel::get_counts(int bin, double n[2][2]) {
    for (int side = 0; side < 2; side++)
        for (int spin = 0; spin < 2; spin++)
            n[side][spin] = counts[side][spin]->GetBinContent(bin);
}


/// Gets the counts and errors in bin
void UCNAmodel::get_counts(int bin, double n[2][2], double e[2][2]) {
    for (int side = 0; side < 2; side++)
        for (int spin = 0; spin < 2; spin++) {
            n[side][spin] = counts[side][spin]->GetBinContent(bin);
            e[side][spin] = counts[side][spin]->GetBinError(bin);
        }
}



/**
 * S := (n[0,0] n[1,1]) / (n[0,1] n[1,0])
 *
    TH1D& compute_super_ratio(TH1D* counts[2][2], TH1D* super_ratio_histogram = 0) 
        if (not test_construction(counts, super_ratio_histogram)) {
            cout<<"Error: computing super ratio.\nAborting.\n";
            exit(1);
        }

        if (not super_ratio_histogram) {
            cout<<"Warning: No super ratio histogram. Copying rate[0][0].\n";
            super_ratio_histogram = new TH1D(*(counts[0][0]));
        }
double UCNAmodel::compute_super_ratio(double n[2][2]) {
    double Y0 = Sqrt(n[0][0]*n[1][1]);
    double Y1 = Sqrt(n[0][1]*n[1][0]);
    // TODO check check positive number
    double S = Y0/Y1;
    return S;
}

double UCNAmodel::compute_super_ratio_error(double e[2][2]) {
    // TODO fix 
    double SE = n[0][0]*n[1][1] / n[0][1]/n[1][0];
}

TH1D& UCNAmodel::compute_super_ratio() {
    int bins = super_ratio->GetNbinsX();
	cout<<"Number of bins "<<bins<<endl;
    // TODO copy other method of bin checking 
    for (int bin = 1; bin < bins; bin++) {
        double n[2][2];
        for (int side = 0; side < 2; side++)
            for (int spin = 0; spin < 2; spin++)
                n[side][spin] = counts[side][spin]->GetBinContent(bin);
        double S = n[0][0]*n[1][1]/n[0][1]/n[1][0];
        if (IsNaN(S)) {
            cout<<"Warning: Super ratio in bin "<<bin<<" is not a number:\n"
                <<"         Was "<<S<<". Setting to zero and continuing.\n";
            S = 0;
        }
        S.SetBinContent(bin, super_ratio);
        S.SetBinError(bin, 0.01);   // TODO compute correctly!!
        cout<<"Warning: Super ratio is not computed correctly.\n";
    }
    return super_ratio
}
 */


/**
 * Y+ := Sqrt(r[0,0] n[1,1])
 * Y- := Sqrt(n[0,1] n[1,0])
 * D := (Y+ - Y-)/2
 * S := (Y+ + Y-)/2
 *
    TH1D* compute_super_sum(TH1D* counts[2][2], TH1D* super_sum = 0, double min = 0, double max = 0) 
        if (not test_construction(counts, super_sum)) {
            cout<<"Error: Problem constructing super sum histogram.\n Aborting.\n";
            exit(1);
        }

        if (not super_sum) {
            cout<<"Warning: Super sum histogram is not constructed.\n";
            super_sum = new TH1D(*(counts[0][0]));
        }
 */
double UCNAmodel::compute_super_sum(double n[2][2]) {
    double Y0 = Sqrt(n[0][0]*n[1][1]);
    double Y1 = Sqrt(n[0][1]*n[1][0]);
    // TODO check if these are numbers.
    double S = (Y0 + Y1)/2;
    if (IsNaN(S)) {
        cout<<"Warning: Super sum is not a number.\n";
        S = -1;
    }
    return S;
}

double UCNAmodel::compute_super_sum(double n[2][2], double e[2][2], double& S, double& sigmaS) {
    double varS = 0;
    S = compute_super_sum(n);
    for (int side=0; side<2; side++) {
        for (int spin=0; spin<2; spin++) {
            double counts = n[side][spin];
            double sigma = e[side][spin];
            if (sigma > 0) {
                double alter = n[side?0:1][spin?0:1];   
                if (counts > 0) {
                    double var = sigma*sigma;
                    varS += var*alter/counts;  /// errors already set
                }
                else
                    varS += alter;             /// using Poisson stats
            } else 
                cout<<"Warning: Negative number in super sum error.\n";
        }
    }

    sigmaS = Sqrt(varS)/2;
    if (IsNaN(sigmaS)) {
        cout<<"Warning: Super sum error multiplier is not a number.\n";
        sigmaS = -1;
    }
    return S;
}

double UCNAmodel::compute_super_sum(int bin, double& S, double& sigmaS) {
    double n[2][2], e[2][2];
    get_counts(bin,n,e);
    compute_super_sum(n,e,S,sigmaS);
    super_sum.SetBinContent(bin,S);
    super_sum.SetBinError(bin,sigmaS);

    if (bin % 10 == 0) {
        double KE = super_sum.GetBinCenter(bin);
        cout<<"Status: Super sum bin "<<bin<<" and KE "<<KE<<" to "<<S<<"("<<sigmaS<<").\n";
    }

    return S;
}

double UCNAmodel::compute_super_sum(int bin) {
    double count, error;
    return compute_super_sum(bin,count,error);
}

TH1D& UCNAmodel::compute_super_sum() {
    int min_bin, max_bin;
    return compute_super_sum(min, max, min_bin, max_bin);
}

TH1D& UCNAmodel::compute_super_sum(double min, double max, 
                                   int& min_bin, int& max_bin) {
    min_bin = super_sum.FindBin(min);
    max_bin = super_sum.FindBin(max);
    if (not super_sum.test_range(min, max)) {
        cout<<"Error: Problem with ranges in super sum histogram.\n";
        exit(1);
    }
    return compute_super_sum(min_bin, max_bin);
}

TH1D& UCNAmodel::compute_super_sum(int min_bin, int max_bin) 
{
    for (int bin = min_bin; bin <= max_bin; bin++)
        compute_super_sum(bin);
    return super_sum;
}


/**
 * Y+ := Sqrt(r[0,0] n[1,1])
 * Y- := Sqrt(n[0,1] n[1,0])
 * D := (Y+ - Y-)/2
 * S := (Y+ + Y-)/2
 * A := D / S
 *
    TH1D* compute_asymmetry(TH1D* counts[2][2], TH1D* asymmetry, double min = 0, double max = 0) 
    {
        if (not test_construction(counts, asymmetry)) {
            cout<<"Error: Bad input histograms for computing asymmetry.\nAborting.\n";
            exit(1);
        }
 */
double UCNAmodel::compute_asymmetry(double n[2][2]) {
    double Y0 = Sqrt(n[0][0]*n[1][1]);
    double Y1 = Sqrt(n[0][1]*n[1][0]);
    double D = (Y0 - Y1)/2;
    double S = (Y0 + Y1)/2;
    double A = D / S;
    if (IsNaN(A)) {
        cout<<"Warning: Asymmetry is not a number.\n";
        S = -1;
    }
    return S;
}

double UCNAmodel::compute_asymmetry(double n[2][2], double e[2][2], double& A, double& sigmaA) {
    double S, sigmaS;
    compute_super_sum(n,e,S,sigmaS);
    A = compute_asymmetry(n);
    sigmaA = Sqrt(1 + A*A) * sigmaS/S;
    if (IsNaN(sigmaA)) {
        cout<<"Warning: Asymmetry error multiplier is not a number.\n";
        sigmaA = -1;
    }
    return A;
}

double UCNAmodel::compute_asymmetry(int bin, double& A, double& sigmaA) {
    double n[2][2], e[2][2];
    get_counts(bin,n,e);
    compute_asymmetry(n,e,A,sigmaA);
    asymmetry.SetBinContent(bin,A);
    asymmetry.SetBinError(bin,sigmaA);

    if (bin % 10 == 0) {
        double KE = super_sum.GetBinCenter(bin);
        cout<<"Status: Asymmetry bin "<<bin<<" and KE "
            <<KE<<" set to "<<A<<"("<<sigmaA<<").\n";
    }

    return A;
}


double UCNAmodel::compute_asymmetry(int bin) {
    double count, error;
    return compute_asymmetry(bin,count,error);
}


TH1D& UCNAmodel::compute_asymmetry() {
    int min_bin, max_bin;
    return compute_asymmetry(min, max, min_bin, max_bin);
}


TH1D& UCNAmodel::compute_asymmetry(double min, double max, 
                                   int& min_bin, int& max_bin) {
    min_bin = asymmetry.FindBin(min);
    max_bin = asymmetry.FindBin(max);
    if (not asymmetry.test_range(min, max)) {
        cout<<"Error: Problem with ranges in asymmetry histogram.\n";
        exit(1);
    }
    return compute_asymmetry(min_bin, max_bin);
}

TH1D& UCNAmodel::compute_asymmetry(int min_bin, int max_bin) 
{
    for (int bin = min_bin; bin <= max_bin; bin++)
        compute_asymmetry(bin);
    return asymmetry;
}

#if 0
TH1D& UCNAmodel::compute_asymmetry(double min, double max) 
{
    int bin_min = test_min(asymmetry, min);
    int bin_max = test_max(asymmetry, max);
    if (not test_range(asymmetry, min, max)) {
        cout<<"Error: Problem with ranges in super sum histogram.\nAborting.\n";
        exit(1);
    }

    for (int bin = bin_min; bin <= bin_max; bin++)
	{
        double n[2][2];
        for (int side = 0; side < 2; side++)
            for (int spin = 0; spin < 2; spin++)
                n[side][spin] = counts[side][spin]->GetBinContent(bin);

        double Lambda = Sqrt(n[0][0]*n[1][1]/(n[0][1]*n[1][0]));
        if (IsNaN(Lambda)) {
            cout<<"Warning: While computing asymmetry, super ratio in bin "
                <<bin<<" is not a number:\n";
                <<"         Was "<<super_ratio
                <<". Setting to zero and continuing.\n";
            Lambda = 0;
        }

		double norm = 1 + super_ratio;  assert(norm > 0);
        double asymmetry = (1 - super_ratio) / norm;
        asymmetry->SetBinContent(bin, asymmetry);

		double inv_sum = Sqrt(1/n[0][0] + 1/n[1][1] + 1/n[0][1] + 1/n[1][0]) / norm;
		double asymmetry_error = super_ratio * inv_sum / norm;  
        if (IsNaN(asymmetry_error) or asymmetry_error <= 0) {
            cout<<"Warning: Asymmetry error in bin "<<bin<<" is not a number:\n";
            cout<<"         Was "<<asymmetry_error<<". Setting to 0.01 and continuing.\n";
            asymmetry_error = 0.01;
        }
        asymmetry->SetBinError(bin, asymmetry_error);
        //printf("Setting bin content for asymmetry bin %d, to %f\n", bin, asymmetry);
    }
    return asymmetry;
}
#endif


/// END UCNAobject codes...

#if 0
//void UCNAFierzFitter::combined_chi2(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *p, Int_t /*iflag */  )
void UCNAFierzFitter::combined_chi2(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *p, Int_t /*iflag */  )
{
	double chi2=0, chi;
    int n = data.asymmetry.bins;
    double A=p[0], b=p[1], N=p[2]; // TODO make nPar correct here
	//double par[3] = {p[0],p[1],p[2]}; // A, b, N
	for (int i = 0; i < n; i++)
	{
		double E = data.asymmetry.GetBinCenter(i);
		double Y = data.asymmetry.GetBinContent(i);
        //double f = asymmetry_fit_func(&E,p);
        double f = A/(1 + b*m_e/(E+m_e));
		//double eY = data.asymmetry.GetBinError(i);
		double eY = 0.1;
        if (eY > 0) {
            chi = (Y-f)/eY;
            chi2 += chi*chi; 
        }
	}

    n = data.super_sum.bins;
	//double par[2] = {p[1], expected[0][1]};
	for (int i=0; i<n; i++) { 
		//E = data.super_sum.energy[i];
		/*double Y =      data .super_sum.values[i];
        double f = p[1]*sm   .super_sum.values[i] 
                 + p[2]*fierz.super_sum.values[i];
        double eY =     data .super_sum.errors[i];*/
		//double E      = data .super_sum.GetBinCenter(i);
		//chi = (fierzratio_values[i] - fierzratio_fit_func(&E,par)) / fierzratio_errors[i];
		double Y  = data .super_sum.GetBinContent(i);
        double eY = data .super_sum.GetBinError(i);
        double f  = N*sm .super_sum.GetBinContent(i) 
                  + N*b*fierz.super_sum.GetBinContent(i);
        if (eY > 0) {
		    chi = (Y-f)/eY;
		    chi2 += chi*chi; 
        }
	}
	fval = chi2; 
}
#endif


double GetEntries(TH1D* histogram, double min, double max)
{
	double entries = histogram->GetEffectiveEntries();
	double part_int = histogram->Integral(
					  histogram->FindBin(min),
					  histogram->FindBin(max));
	double full_int = histogram->Integral();
	double N = entries * part_int / full_int;

	return N;
}


double evaluate_expected_fierz(double min, double max, double integral_size = 1234) 
{
    TH1D *h1 = new TH1D("beta_spectrum_fierz", 
                        "Beta spectrum with Fierz term", 
                        integral_size, min, max);
    TH1D *h2 = new TH1D("beta_spectrum",
                        "Beta Spectrum", 
                        integral_size, min, max);
	for (int i = 0; i < integral_size; i++)
	{
		double K = min + double(i)*(max-min)/integral_size;
		double par1[2] = {0, 1};
		double par2[2] = {0, 0};
		double y1 = fierz_beta_spectrum(&K, par1);
		double y2 = fierz_beta_spectrum(&K, par2);
		h1->SetBinContent(K, y1);
		h2->SetBinContent(K, y2);
	}
	return h1->Integral(0, integral_size) / h2->Integral(0, integral_size);
}


/// beta spectrum with little b term
double fierz_beta_spectrum(const double *val, const double *par)
{
	const double K = val[0];                    /// kinetic energy
	if (K <= 0 or K >= Q)
		return 0;                               /// zero outside range

	const double b = par[0];                    /// Fierz parameter
	const int n = par[1];                    	/// Fierz exponent
	const double E = K + m_e;                   /// electron energy
	const double e = Q - K;                     /// neutrino energy
	const double p = sqrt(E*E - m_e*m_e);       /// electron momentum
	const double x = pow(m_e/E,n);              /// Fierz term
	const double f = (1 + b*x)/(1 + b*x_1);     /// Fierz factor
	const double k = 1.3723803E-11/Q;           /// normalization factor
	const double P = k*p*e*e*E*f*x;             /// the output PDF value

	return P;
}


/// beta spectrum with expected x^-n and beta^m
double beta_spectrum(const double *val, const double *par)
{
	const double K = val[0];                    	///< kinetic energy
	if (K <= 0 or K >= Q)
		return 0;                               	///< zero beyond endpoint

	const double m = par[0];                    	///< beta exponent
	const double n = par[1];                    	///< Fierz exponent
	const double E = K + m_e;                   	///< electron energy
	const double B = pow(1-m_e*m_e/E/E,(1+m)/2);  	///< beta power factor
	const double x = E / m_e;                   	///< reduced electron energy
	const double y = (Q - K) / m_e;             	///< reduced neutrino energy
	const double z = pow(x,2-n);          			///< Fierz power term
	const double k = 1.3723803E-11/Q;           	///< normalization factor
	const double P = k*B*z*y*y;             		///< the output PDF value

	return P;
}


double evaluate_expected_fierz(double m, double n, double min, double max, double integral_size = 1234)
{
    // TODO These don't need to be pointers.
    TH1D *h1 = new TH1D("beta_spectrum_fierz", 
                        "Beta spectrum with Fierz term", 
                        integral_size, min, max);
    TH1D *h2 = new TH1D("beta_spectrum", 
                        "Beta Spectrum", 
                        integral_size, min, max);
	for (int i = 0; i < integral_size; i++)
	{
		double K = min + double(i)*(max-min)/integral_size;
		double par1[2] = {m, n};
		double par2[2] = {0, 0};
		double y1 = beta_spectrum(&K, par1);
		double y2 = beta_spectrum(&K, par2);
		h1->SetBinContent(K, y1);
		h2->SetBinContent(K, y2);
	}
	double rv = h1->Integral(0, integral_size) / h2->Integral(0, integral_size);
    delete h1;
    delete h2;
    return rv;
}


double evaluate_expected_fierz(double min, double max)
{
	return evaluate_expected_fierz(0, 1, min, max, 1234);
}
