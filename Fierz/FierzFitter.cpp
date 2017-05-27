#include "FierzFitter.hh"
#include <regex>


UCNAmodel & UCNAmodel::operator=(const UCNAmodel & other)
{
    /** TODO check that if name and title are already assigned, to NOT copy them
    name(other.name), title(other.title),
    bins(other.bins), min(other.min), max(other.max),
    rand(0)
    super_ratio(other.super_ratio),
    super_sum(other.super_sum),
    asymmetry(other.asymmetry) */
    if (name == "")
        name = other.name + "_copy";
    if (title == "")
        title = other.title + " Copy";
    if (min > other.min)
        min = other.min;
    if (max < other.max)
        max = other.max;
    // TODO check if bins are bigger
        bins = other.bins;
    super_ratio = other.super_ratio;
    super_ratio.SetName(name+"_super_ratio");
    super_ratio.SetTitle(title+" Super Ratio");

    super_sum = other.super_sum;
    super_sum.SetName(name+"_super_sum");
    super_sum.SetTitle(title+" Super Sum");

    asymmetry = other.asymmetry;
    asymmetry.SetName(name+"_asymmetry");
    asymmetry.SetTitle(title+" Asymmetry");

    for (int side = 0; side < 2; side++) {
        TString sub_name = name;    /// keep name the same
        TString sub_title = title;  /// keep title the same
        if (not side) {
            sub_name += "_E";
            sub_title += " East";
        } else {
            sub_name += "_W";
            sub_title += " West";
        }
        for (int spin = 0; spin < 2; spin++) {
            if (not spin) {
                sub_name += "_off";
                sub_title += " AFP Off";
            } else {
                sub_name += "_on";
                sub_title += " AFP On";
            }
            *counts[side][spin] = *other.counts[side][spin];
            counts[side][spin]->SetName(sub_name);
            counts[side][spin]->SetTitle(sub_title);
        }
    }
    ntuple = other.ntuple;
    return *this;
}

bool UCNAhistogram::test_compatable(UCNAhistogram & other)
{
    int error = 0;
    double min = GetXaxis()->GetXmin();
    double max = GetXaxis()->GetXmax();
    double other_min = other.GetXaxis()->GetXmin();
    double other_max = other.GetXaxis()->GetXmax();
    double inner_min = min > other_min? min : other_min;
    double inner_max = max < other_max? max : other_max;
    int bin_min = FindBin(inner_min);
    int bin_max = FindBin(inner_max);
    /*
    if (min != other_min or max != other_max) {
        cout<<"Warning: Overall bins are not identical. They still might be compatable.\n";
        cout<<"         This min is "<<min<<" and the other min is "<<other_min<<".\n";
        cout<<"         This max is "<<max<<" and the other max is "<<other_max<<".\n";
    }
    */
    int other_bin_min = other.FindBin(inner_min);
    //int other_bin_max = other.FindBin(inner_max);
    for (int bin=bin_min; bin<bin_max; bin++) {
        int other_bin = bin - bin_min + other_bin_min;
        double center = GetBinCenter(bin);
        double other_center = other.GetBinCenter(other_bin);
        if (center != other_center) {
            cout<<"Error: Centers in bin "<<bin<<" do not match.\n";
            cout<<"       This center is "<<center<<" and the other's is "<<other_center<<".\n";
            error++;
        }

        double width = GetBinWidth(bin);
        double other_width = other.GetBinWidth(bin);
        if (width != other_width) {
            cout<<"Error: Bin widths in bin "<<bin<<" do not match.\n";
            cout<<"       This width is "<<width<<" and the other's is "<<other_width<<".\n";
            error++;
        }
    }
    return (error > 0);
}

/// Normalize to units of probability per MeV over default range
double UCNAhistogram::normalize() {
    double min = GetXaxis()->GetXmin();
    double max = GetXaxis()->GetXmax();
    test_range();
    return normalize(min,max);
}

/// Normalize to units of probability per MeV
double UCNAhistogram::normalize(double min, double max)
{
    //Sumw2();    /// set up sum of squares structure.
    assert(min < max);
	int bin_min = FindBin(min);
	int bin_max = FindBin(max);
    assert(bin_min < bin_max);
    double integrand = 0;
    /* TODO ???
    integrand +=(GetBinCenter(bin_min) - min)
               / GetBinWidth(bin_min)
               / GetBinWidth(bin_min)
               * GetBinContent(bin_min)
               +(max - GetBinCenter(bin_max))
               / GetBinWidth(bin_max)
               / GetBinWidth(bin_max)
               * GetBinContent(bin_max);
               */

    //for (int bin=bin_min+1; bin<bin_max; bin++)
    for (int bin=bin_min; bin<bin_max; bin++)
    {
        double area = GetBinContent(bin)
                    * GetBinWidth(bin);
        if (IsNaN(area))
            cout<<"Warning: Differential area is not a number.\n";
        else
            integrand += area;
    }
    if (integrand > 0)
	    Scale(GetEntries()/integrand);
    else 
        cout<<"Warning: Integrand is zero. Leaving unscaled.\n";

    return integrand;
}


TF1* UCNAFierzFitter::combined_fit(
        TMatrixD &cov, TF1 *func, 
        void (*global_fcn_ptr)(Int_t&, Double_t*, Double_t&, Double_t*, Int_t) )
{
    int nPar = func->GetNpar();
    if (not func) {
        cout<<"Fit function not set.\n";
        exit(1);
    }

	/// initialize the Minuit fitter
	TVirtualFitter::SetDefaultFitter("Minuit");
	TVirtualFitter *minuit = TVirtualFitter::Fitter(0,nPar);
	for (int i=0; i<nPar; i++) {
		minuit->SetParameter(i, func->GetParName(i), func->GetParameter(i), 1e-3, 0, 0);
    }

    /// set up chi^2 function to minimize
	minuit->SetErrorDef(1);	        /// 1 for chi^2
	minuit->SetFCN(global_fcn_ptr); /// needs to be a global function

	/// set print level
	double arglist[100];
	arglist[0] = 0;
	minuit->ExecuteCommand("SET PRINT", arglist, 1);

	/// minimization parameters
	arglist[0] = 5000;     /// number of function calls
	arglist[1] = 0.01;     /// tolerance

    /// execute minimization.
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

    /// find degrees of freedom in fit
	//int ndf = fit.asymmetry.GetNbinsX() + fit.super_sum.GetNbinsX() - nvpar;
	int ndf = fit.super_sum.GetNbinsX() - nvpar;
	func->SetNDF(ndf);
    
	TMatrixD matrix( nPar, nPar, minuit->GetCovarianceMatrix() );
	for (int i=0; i<nPar; i++)
		for (int j=0; j<nPar; j++)
			cov[i][j] = minuit->GetCovarianceMatrixElement(i,j);
	
	cout<<"\n    chi^2 = "<<chi2<<", ndf = "<<ndf<<", chi^2/ndf = "<<chi2/ndf<<".\n\n";
	return func; 
}


/* TODO some of the details of this are still wrong...
double UCNAhistogram::chi2(const UCNAhistogram & fit) {
	double chi2 = 0;
    int n = fit.FindBin(min);
    int delta = FindBin(fit_min) - 1;
	for (int bin=1; bin<n; bin++)
	{
		double Y = GetBinContent(bin + delta);
		double eY = GetBinError(bin + delta);
        double F = fit.GetBinContent(bin);
        double eF = fit.GetBinError(bin);
        double eYF2 = eY*eY + eF*eF;
        if (eYF2 > 0)
            chi2 += (Y-F)*(Y-F)/eYF2; 
	}
    return chi2;
}
*/


double UCNAhistogram::chi2(const UCNAhistogram &fit) {
	double chi2 = 0;
    int n = fit.GetNbinsX();
    double min = fit.GetXaxis()->GetXmin();
    int delta = FindBin(min) - 1;
	for (int bin=1; bin<n; bin++)
	{
		double Y = GetBinContent(bin+delta);
		double eY = GetBinError(bin+delta);
        double F = fit.GetBinContent(bin);
        double eF = fit.GetBinError(bin);
        double eYF2 = eY*eY + eF*eF;
        if (eYF2 > 0)
            chi2 += (Y-F)*(Y-F)/eYF2; 
        else
        {
            if (eYF2 == 0) {
                cout<<"Warning: division-by-zero while summing chi2.\n";
            }
            if (eYF2 < 0) {
                cout<<"Warning: negative square error while summing chi2: "<<eYF2<<".\n";
            }
            if (IsNaN(eYF2)) {
                cout<<"Warning: error is not a number while summing chi2.\n";
            }
            //exit(0);
        }
	}
    if (chi2 <= 0)
        cout<<"Error: Non-positive total chi^2 fitting "<<fit.GetTitle()<<" to "<<GetTitle()<<".\n";
    return chi2;
}


double UCNAFierzFitter::asymmetry_chi2(double A, double b) {
    assert(false);
    compute_asymmetry_fit(A,b);
    return data.asymmetry.chi2(fit.asymmetry);
    /*
	double chi2 = 0;
    int n = fit.asymmetry.GetNbinsX();
    int delta = data.asymmetry.FindBin(fit_min) - 1;
	for (int bin=1; bin<n; bin++)
	{
		double Y = data.asymmetry.GetBinContent(bin + delta);
		double eY = data.asymmetry.GetBinError(bin + delta);
        double F = fit.asymmetry.GetBinContent(bin);
        double eF = fit.asymmetry.GetBinError(bin);
        double eYF2 = eY*eY + eF*eF;
        if (eYF2 > 0)
            chi2 += (Y-F)*(Y-F)/eYF2; 
        else
            cout<<"Warning: Encountered non-positive error in asymmetry fit.\n";
	}
    if (chi2 <= 0)
        cout<<"Error: Non-positive total chi^2 in asymmetry fit.\n";
    return chi2;
    */
}


double UCNAmodel::asymmetry_chi2(double A, double b) {
    assert(false); /*
	double chi2 = 0;
    int n = asymmetry.GetNbinsX();
	for (int i=1; i<n; i++)
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
    */
}


double UCNAFierzFitter::supersum_chi2(double b, double N)
{
    assert(N > 0);
    compute_supersum_fit(b,N);
	double chi2 = 0;
    int n = fit.super_sum.GetNbinsX();
    int delta = data.super_sum.FindBin(fit_min) - 1;
	for (int bin=1; bin<n; bin++)
	{
		double Y = data.super_sum.GetBinContent(bin + delta);
		double eY = data.super_sum.GetBinError(bin + delta);
        double F = fit.super_sum.GetBinContent(bin);
        double eF = fit.super_sum.GetBinError(bin);
        double eYF2 = eF*eF + eY*eY;
        if (eYF2 > 0)
            chi2 += (Y-F)*(Y-F)/eYF2; 
	}
    return chi2;
}

double UCNAFierzFitter::combined_chi2(double A, double b, double N)
{
	double chi2 = //asymmetry_chi2(A,b) + 
    supersum_chi2(b,N);
    return chi2;
}

void UCNAFierzFitter::compute_fit(/*MatrixD &cov,*/ TF1 *func) {
    assert(false);
    double A = func->GetParameter(0);
    double b = func->GetParameter(1);
    double N = func->GetParameter(2);
    compute_fit(A,b,N);
}

void UCNAFierzFitter::compute_fit(double A, double b, double N) {
    assert(N > 0);
    //compute_asymmetry_fit(A,b);
    compute_supersum_fit(b,N);
}

void UCNAFierzFitter::compute_data(double A, double b, double N) {
    assert(N > 0);
    //compute_asymmetry_fit(A,b);
    compute_supersum_fit(b,N);
    data.super_sum = fit.super_sum;
}


void UCNAFierzFitter::compute_supersum_fit(double b, double N)
{
    assert(N > 0);
    int n = fit.super_sum.GetNbinsX();
    // TODO add axial supersum?
    int deltaV = vector.super_sum.FindBin(fit_min) - 1;
    int deltaF = fierz.super_sum.FindBin(fit_min) - 1;
    double exV = vector.super_sum.GetEffectiveEntries(fit_min, fit_max);
    double exF = fierz.super_sum.GetEffectiveEntries(fit_min, fit_max);
    //double e_x = evaluate_expected_fierz(fit_min,fit_max); // not in the right range
    //double m_E = evaluate_expected_fierz(0,Q); // not in the right range
    double norm = exV + b*x_1*exF;
    cout <<"norm = exV + b*x_1*exF : "<<norm<<" = "<<exV<<"+"<<b<<"*"<<x_1<<"*"<<exF<<".\n";
	for (int bin=1; bin<n; bin++)
	{
		double KE = fit.super_sum.GetBinCenter(bin);
        if (KE < fit_min or KE > fit_max) {
            cout<<"Error: Found energy out of bounds in fit\n";
            exit(1);
        }
        double pV = vector.super_sum.GetBinContent(bin + deltaV);
        assert(pV >= 0);
        double pF = fierz.super_sum.GetBinContent(bin + deltaF);
        assert(pF >= 0);
        double f  = N*(pV + b*x_1*pF)/norm;
        assert(f >= 0);
        fit.super_sum.SetBinContent(bin,f);

        double eV = vector.super_sum.GetBinError(bin + deltaV)/norm;
        double eF = fierz.super_sum.GetBinError(bin + deltaF)/norm;
        double ef = N*Sqrt(eV*eV + b*b*x_1*x_1*eF*eF);
        fit.super_sum.SetBinError(bin,ef);

        //cout<<"Super sum fit bin: "<<bin<<"\tKE: "<<KE<<"\tSS: "<<f<<"("<<ef<<")\n";
	}
}

void UCNAFierzFitter::fill(TString vector_pattern, 
                           TString axial_pattern, 
                           // TODO TString axial_up_pattern, 
                           // TODO TString axial_down_pattern, 
                           TString fierz_pattern,
                           int min, int max, /// TODO read filename pattern
                           TString name,
                           int type,
                           double flip)
{
    if (flip == -1)
        flip = spin_ratio;

    /// Fill Standard Model vector events.
    if (vector_pattern != "")
        vector.fill(
            vector_pattern, 
            min, max,     /// TODO read from filename pattern
            name,
            "Standard Model vector current", 
            type, flip);

    /// Fill Standard Model axial-vector events
    /// 
    /// TODO fill A+ and A- axial vector events 
    /// axial[0].fill( ...
    /// axial[1].fill( ...
    /// -- or --
    /// for (i;1,2) axial[i].fill( ...
    if (axial_pattern != "")
        axial.fill(
            axial_pattern,
            min, max,     /// TODO read from filename pattern
            name,
            "Standard Model axial-vector current", 
            type, flip);

    /// Fill BSM Fierz events
    if (fierz_pattern != "")
        fierz.fill(
            fierz_pattern,
            min, max,     /// TODO read from filename pattern
            name,
            "BSM Fierz current", 
            type, flip);
}


/// DISPLAYING AND OUTPUTTING
void UCNAFierzFitter::display(TString &plots_base) {
    TCanvas canvas("fierz_fitter_canvas",
                   "Combined Fierz component of energy spectrum");
    TLegend legend(0.55,0.65,0.85,0.85);

    /* TODO Uncomment when ready for fitting asymmetry
    data.asymmetry.draw(
                   "data_asymmetry", 
                   title+" fit #Lambda(E)", 
                   &canvas,&legend,"",1,0);
    canvas.SaveAs(plots_base+"data_asymmetry.pdf");
    fit.asymmetry.draw(
                   "fit_asymmetry", 
                   title+" fit #Lambda(E)", 
                   &canvas,&legend,"Same",6,0);
    canvas.SaveAs(plots_base+"data_fit_asymmetry.pdf"); */

    vector.super_sum.draw(
                   "vector_supersum", 
                   title+" vector Monte Carlo #Sigma(E)", 
                   &canvas,&legend,"",4,0);
    /* TODO Uncomment when ready for fitting asymmetry
    axial.super_sum.draw(
                   "axial_supersum", 
                   title+" axial-vector Monte Carlo #Sigma(E)", 
                   &canvas,&legend,"Same",8,0); */
    fierz.super_sum.draw(
                   "fierz_supersum", 
                   title+" Fierz Monte Carlo #Sigma(E)", 
                   &canvas,&legend,"Same",6,0);
    canvas.SaveAs(plots_base+"monte_carlo.pdf");

    data.super_sum.draw(
                   "data_supersum", 
                   title+" #Sigma(E)", 
                   &canvas,&legend,"",1,0);
    canvas.SaveAs(plots_base+"data_supersum.pdf");
    fit.super_sum.draw(
                   "fit_supersum", 
                   title+" fit #Sigma(E)", 
                   &canvas,&legend,"Same",6,0);
    canvas.SaveAs(plots_base+"data_fit_supersum.pdf");

    /* /// Output for gnuplot
	//save(plots_base+"super-sum-data.dat", data.super_sum, 1, 1000);
	//save(plots_base+"super-sum-mc.dat", sm.super_sum, 1, 1000);
	//save(plots_base+"fierz-ratio.dat", fierz_ratio_histogram, 1, 1);
	//save(plots_base+"fierz-fit.dat", fierz_fit_histogram, 1, 1); */
}

/// compute effective size
double UCNAFierzFitter::comupte_sizes() {
    Nsim_data = data.super_sum.GetEntries();
    Nall_data = data.super_sum.GetEffectiveEntries(min, max);
    Nfit_data = data.super_sum.GetEffectiveEntries(fit_min, fit_max);
    Nfit_vector = vector.super_sum.GetEffectiveEntries(fit_min, fit_max);
    // TODO Nfit_axial = axial.super_sum.GetEffectiveEntries(fit_min, fit_max);
    // TODO Nfit_fierz = fierz.super_sum.GetEffectiveEntries(fit_min, fit_max);
    //double N = Nfit_data*Nfit_vector/Sqrt(Nfit_data*Nfit_data + Nfit_vector*Nfit_vector);
    Neff = Nfit_data*Nfit_vector/(Nfit_data + Nfit_vector);
    return Neff;
}


/// Output the number of events for data, mcs and fits and with cut info.
double UCNAFierzFitter::print_sizes() {
    cout<<setprecision(5);
    //cout<<" NUMBER OF EVENTS WITH ENERGY RANGE CUTS:\n";
    cout<<"    Full energy range is "<<min<<" - "<<max<<" keV.\n";
    cout<<"    Fit energy range cut is "<<fit_min<<" - "<<fit_max<<" keV.\n";
    cout<<"    Number of actual counts in data is "<<int(Nsim_data)<<".\n";
    cout<<"    Entries full energy range is "<<int(Nsim_data)<<".\n";
    cout<<"    Effective number of counts in full energy range is "<<int(Nall_data)<<".\n";
    cout<<"    Effective number of counts in fit energy range cut is "<<int(Nfit_data)<<".\n";
    cout<<"    Efficiency of energy cut is "<<int(Nfit_data/Nall_data*10000)/100<<"%.\n";
    //cout<<"    Effective N "<<int(Neff*100)/100<<".\n";
    return Neff;
}


/* XXX old formula method
void UCNAFierzFitter::compute_asymmetry_fit(double A, double b)
{
    int n = fit.asymmetry.GetNbinsX();
	for (int bin=1; bin<n; bin++)
	{
		double KE = fit.asymmetry.GetBinCenter(bin);
        if (KE < fit_min or KE > fit_max) {
            cout<<"Error: Found energy out of bounds in fit\n";
            exit(1);
        }
        double E = KE + m_e;
        double f = A/(1 + b*(m_e/E));
        double ef = 0;  // TODO get from MC asymmetry (for 2011-2013 data)
        fit.asymmetry.SetBinContent(bin,f);
        fit.asymmetry.SetBinError(bin,ef); 
        //cout<<"Asymmetry fit bin: "<<bin<<"\tKE: "<<KE<<"\tA: "<<f<<"("<<ef<<")\n";
	}
}
*/

 void UCNAFierzFitter::compute_asymmetry_fit(double A, double b/*, double N*/)
{
    cout<<"Function compute_asymmetry_fit.\n";
    exit(0);
    int n = fit.asymmetry.GetNbinsX();
    int deltaA = axial.asymmetry.FindBin(fit_min) - 1;
    int deltaV = vector.asymmetry.FindBin(fit_min) - 1;
    int deltaF = fierz.asymmetry.FindBin(fit_min) - 1;
    //double ex_cos = evaluate_expected_fierz(fit_min,fit_max);
    //double ex_x = evaluate_expected_fierz(fit_min,fit_max);
    double ex_cos = 1/2;
	for (int bin=1; bin<n; bin++)
	{
		double KE = fit.asymmetry.GetBinCenter(bin);
        if (KE < fit_min or KE > fit_max) {
            cout<<"Error: Found energy out of bounds in fit.\n";
            exit(1);
        }
        //double E = KE + m_e;
        //double x = m_e/E;
        //double beta_cos = Sqrt(1 - x*x)/2;
        double pA = axial.asymmetry.GetBinContent(bin + deltaA);
        double pV = vector.super_sum.GetBinContent(bin + deltaV);
        double pF = fierz.super_sum.GetBinContent(bin + deltaF);
        //double xFR = pF/pV;
        //double f = A*pA*beta_cos/(1 + b*x);
        //double f = A*pA/(1 + xFRMC);
        double f = A*pA/(pV + b*pF)/ex_cos;
        fit.asymmetry.SetBinContent(bin,f);

        double eA = axial.asymmetry.GetBinError(bin + deltaA);
        double eV = vector.super_sum.GetBinError(bin + deltaV);
        double eF = fierz.super_sum.GetBinError(bin + deltaF);
        //double ef = N*Sqrt(eV*eV + e_x*e_x*beF*beF);
        //double eFRMC = Sqrt(ebF*ebF + eV*eV);
        double ef = A*pF/pA*Sqrt(eA*eA + eV*eV + b*b*eF*eF)/ex_cos;
        fit.asymmetry.SetBinError(bin,ef);

        //cout<<"Asymmetry fit: bin: "<<bin<<"\tKE: "<<KE<<"\tA: "<<f<<"("<<ef<<")\n";
	}
}

double asymmetry_fit_func(double *x, double *p)
{
    exit(0);
    double E = x[0] + m_e;
    double A = p[0];
    double b = p[1];
    return A/(1 + b*m_e/E); 
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
	TFile* tfile = new TFile(filename);
	if (tfile->IsZombie()) {
		cout<<"Error: While loading "<<title<<":\n";
		cout<<"       File "<<filename<<" not found.\n";
		return 0;
	}

    // TODO warn if overwriting data...
    TH1D* histogram = (TH1D*)tfile->Get(name);
    if (not histogram) {
		cout<<"Error: In file "<<filename<<":\n";
		cout<<"       Problem getting "<<title<<".\n";
		cout<<"       Cannot find histogram named "<<name<<".\n";
        return 0;
    }
    *(TH1D*)this = *histogram;

	int entries = GetEntries();
	cout<<"    Number of entries in "<<title<<" is "<<entries<<".\n";
    return entries;
}


/** 
 * UCNAhistogram::draw
 * overload TH1's draw to set up certain features
 */
void UCNAhistogram::draw(TString name, TString title,
          TCanvas* canvas, TLegend* legend, 
          TString draw = "", int color = 0, int marker = 0)
{
    if (not canvas)
        canvas = new TCanvas(name + "_canvas", title + " Canvas");
    if (not canvas) {
        cout<<"Error: Can't construct a new canvas for "<<name<<".\n";
        exit(1);
    }

    /// Draw a histogram.
	SetStats(0);
    SetLineColor(color);
    SetMarkerStyle(marker);
    Draw(draw);

    /// Make a pretty legend.
    if (legend) {
        if (draw == "")
            legend->Clear();
        legend->SetTextSize(0.033);
        legend->SetBorderSize(0);
        legend->SetFillColor(0);
        legend->AddEntry(this, title, "lep");
        legend->Draw();
    }

    /// Save the data and Mote Carlo plots.
    //TString filename = plots_base + name + ".pdf";
    //canvas->SaveAs(filename);
}




/// DISPLAYING AND OUTPUTTING
/**
 *
 */
void UCNAhistogram::save(TString filename, double ax, double ay)
{
    ofstream ofs;
	ofs.open(filename);
	for (int i=1; i<=GetNbinsX(); i++)
	{
		double x = ax * GetBinCenter(i);
		double sx = GetBinWidth(i);
		double r = ay * GetBinContent(i);
		double sr = ay * GetBinError(i);
		ofs<<x<<'\t'<<r<<'\t'<<sx<<'\t'<<sr<<endl;
	}
	ofs.close();
}


void UCNAFierzFitter::save(TString filename)
{
    vector.super_sum.save(filename+"-vector.dat");
    fierz.super_sum.save(filename+"-fierz.dat");
    data.super_sum.save(filename+"-data.dat");
    fit.super_sum.save(filename+"-fit.dat");
}


/**
 * UCNAhistogram::save() and UCNAmodel::save()
 * Save root data
 */
void UCNAmodel::save(TString filename)
{
	TFile* tfile = new TFile(filename, "recreate");
	if (tfile->IsZombie()) {
		cout<<"Error: Problem saving "<<title<<":\n";
		cout<<"       Cannot create file "<<filename<<".\n";
        cout<<"Aborting...\n";
		exit(1);
	}

    for (int side=0; side<2; side++)
        for (int spin=0; spin<2; spin++) {
            counts[side][spin]->SetDirectory(tfile);
            counts[side][spin]->Write();
        }

    super_sum.SetDirectory(tfile);
    super_sum.Write();

    asymmetry.SetDirectory(tfile);
    asymmetry.Write();

    /* if (ntuple) {
        ntuple->SetDirectory(tfile);
        ntuple->Write();
    } else {
        cout<<"Warning: Ntuple not set. Can't save data.\n";
    } */
	tfile->Close();
}


/** DEPRECATED *******
 * UCNAhistogram::fill() and UCNAmodel::fill()
 * Fill in asymmetry and super_ratio, and super sums from simulation data.
 * Use wild card * in filename where data is split up over many files
 * and they get Tchained together.
int UCNAmodel::fill(TString filename, TString name, TString title) {
    return fill(filename, name, title, 0.68/1.68, -0.12, 0.0);
} */


int UCNAmodel::fill(TString pattern, int first, int last, 
                    TString name, TString title, 
                    int type, double flip)
{
    assert(flip > 0);
    /* TODO add pattern matching for file numbers
    cout <<pattern.SubString("[",3)<<"\n";
    cout <<pattern.SubString("[]",0)<<"\n";
    //cout <<pattern.Contains("[7771234]")<<"\n";
    cout <<pattern(TRegexp("[*]"))<<"\n";
    static regex ex(pattern);
    sregex_iterator next(ex
    */
    TChain *chains = 0;
    int added = 0;
    for (int i = first; i <= last; i++) {
        TString filename(pattern);
        TString number; number.Form("%d",i); 
        filename.ReplaceAll("[*]",number);
        TFile* tfile = new TFile(filename);
        if (not tfile->IsZombie()) {
            cout<<"Loading filename "<<filename<<".\n";
            if (not chains) {
                assert(added==0);
                chains = new TChain(name, title);
            }
            chains->Add(filename);
            added++;
        }
        else {
            cout<<"Error: Problem filling "<<title<<":\n";
            cout<<"       File "<<filename<<" not found\n";
            cout<<"       from "<<pattern<<".\n";
            cout<<"Skipping...\n";
        }
    }
    if (added > 0 and chains) {
        return fill(chains, type, flip); 
    }
    else {
		cout<<"Error: In file pattern "<<pattern<<":\n";
		cout<<"       Cannot get "<<title<<":\n";
		cout<<"       Cannot find any chain or tree named "<<name<<" in any file with pattern "<<pattern<<".\n";
        cout<<"Aborting...\n";
        exit(1);
    }
}


int UCNAmodel::fill(TString filename, TString name, TString title, 
                    int type, double flip)
{
    /// TODO check if filename is a pattern, and if so call chain version above
    if (filename.MaybeWildcard()) {
        /// find * wild card 
        cout << filename.SubString("*",1);
        // TODO call chain version
        //fill(filename,1,1,name,title,type,flip);
    }
	TFile* tfile = new TFile(filename);
	if (tfile->IsZombie()) {
		cout<<"Error: Problem filling "<<title<<":\n";
		cout<<"       File "<<filename<<" not found\n";
        return 0;
        cout<<"Aborting...\n";
		exit(1);
	}
    TChain *chain = new TChain(name, title);
    chain->Add(filename);
    //TChain *chain = (TChain*)tfile->Get(name);
    if (not chain) {
		cout<<"Error: In file "<<filename<<":\n";
		cout<<"       Cannot get "<<title<<":\n";
		cout<<"       Cannot find chain or tree named "<<name<<".\n";
        return 0;
        cout<<"Aborting...\n";
        exit(1);
    }
    else {
        SetAllBranches(chain);
        int Nsim = fill(chain, type, flip);
        if (Nsim <= 0) {
            cout<<"Error: No events were filled.\n";
            return Nsim;
            cout<<"Aborting...\n";
            exit(1);
        }
        else
            return Nsim;
    }
}

void UCNAmodel::SetAllBranches(TChain *chain)
{
	chain->SetBranchStatus("*",false);
	chain->SetBranchStatus("PID",true);
	chain->SetBranchStatus("side",true);
	chain->SetBranchStatus("type",true);
	chain->SetBranchStatus("Erecon",true);
	chain->SetBranchStatus("primMomentum",true);
    chain->SetBranchStatus("ScintPosAdjusted",true);
    chain->SetBranchStatus("ScintPosAdjusted",true);
}

int UCNAmodel::fill(TChain *chain, int type, double flip)
{
    assert(flip > 0);

    if (not chain) {
		cout<<"Error: Null TChain pointer.\n";
    }
	double Nevents = chain->GetEntries();
    SetAllBranches(chain);

	//TNtuple* tntuple = new TNtuple("mc_ntuple", "MC NTuple", "s:load:energy"); */
	TNtuple* tntuple = 0;

	//int Noff=0, Non=0;	/// events have been simulated after cuts
	int Nspin[]={0,0};	    /// events have been simulated after cuts
    int PID, side, event_type;
    double energy;
    double mwpcPosW[3], mwpcPosE[3], p[3];
    //double mwpcPos[2][3], p[3]; 
    //double ex_cos = 4;  // TODO compute ex+cos correctly

    chain->SetBranchAddress("PID",&PID);
    chain->SetBranchAddress("side",&side);
    chain->SetBranchAddress("type",&event_type);
    chain->SetBranchAddress("Erecon",&energy);
	chain->SetBranchAddress("primMomentum",p); /// not momentum; normalized direction!
    chain->GetBranch("ScintPosAdjusted")->GetLeaf("ScintPosAdjE")->SetAddress(mwpcPosE);
    chain->GetBranch("ScintPosAdjusted")->GetLeaf("ScintPosAdjW")->SetAddress(mwpcPosW);

    for (int evt=0; evt<Nevents; evt++) {
        chain->GetEvent(evt);

        /// cut out bad events
        if (PID!=1) 
            continue;

        /// radial fiducial cut
        double fidcut2 = 50*50; // TODO !!!! VERY HARD CODED 
        assert(fidcut2 > 0);
        double radius2W = mwpcPosW[0]*mwpcPosW[0] + mwpcPosW[1]*mwpcPosW[1]; 
        double radius2E = mwpcPosE[0]*mwpcPosE[0] + mwpcPosE[1]*mwpcPosE[1]; 
        if (radius2W > fidcut2 or radius2E > fidcut2) 
            continue;

        /// fill with loading efficiency 
        double load = rand.Uniform(1);
        double afp = (load < flip)? +1 : -1; // TODO don't need load as var..
        //double axial = rand.Uniform(2)-1;
        double cos = p[2]/Sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
        assert(cos = p[2]);
        //bool spin = ((afp < 0) xor (cos/ex_cos < A)) ? EAST : WEST;
        //double E = energy + m_e;        /// relativistic energy
        //double x = m_e/E;
        bool spin = (afp < 0)? EAST : WEST;
        Nspin[spin]++;
        /*
        if (b > 0.5) {
            double d = 1 + b*x;
            if (afp < 0 ) {
                //Noff++;
                if (axial > (1 + A*beta*cos/d))
                    spin = EAST;
                else 
                    spin = WEST;
            }
            else {
                //Non++;
                if (axial > (1 - A*beta*cos/d))
                    spin = WEST;
                else
                    spin = EAST;
            }
        } */
        /*cout<<"energy: "<<energy<<" side: "<<side<<" spin: "<<spin 
            <<" afp: "<<afp<<" p: ("<<p[0]<<","<<p[1]<<","<<p[2]<<")" 
            <<" cos: "<<cos<<" load: "<<load<<".\n";*/
        if (event_type==type) 
            counts[side][spin]->Fill(energy, 1);
        /*
        if (type<4) 
        if (type==2 or type==3)
            counts[side][spin]->Fill(energy, 1);
        counts[side][spin]->Fill(energy, 1); */
        if (tntuple)
            tntuple->Fill(side, spin, energy);

		/// break when enough data has been generated.
        /*
		if(nSimmed >= nToSim)
		break;
        */
    }    
     
    int Nsim = Nspin[0] + Nspin[1];
	cout<<"Total number of entries in "<<name<<":\n";
	cout<<"      Entries without cuts:      "<<int(Nevents)<<endl;
	cout<<"      Entries with cuts:         "<<Nsim<<endl;
	cout<<"      Efficiencies of cuts:      "<<(100.0*Nsim)/double(Nevents)<<"%\n";
	cout<<"      Entries with AFP off cut:  "<<Nspin[0]<<" ("<<(100.0*Nspin[0])/double(Nsim)<<"%)\n";
	cout<<"      Entries with AFP on cut:   "<<Nspin[1]<<" ("<<(100.0*Nspin[1])/double(Nsim)<<"%)\n";

	/// compute and normalize super sum
    compute_super_sum();
    //compute_asymmetry();
    //compute_super_diff();
    //super_sum.normalize();
    if (Nsim > 0)
        super_sum.Scale(1/double(Nsim));
    //asymmetry.Scale(1/A);
    //super_sum.Scale(1/double(Nevents));
    //super_sum.Scale(double(nSimmed)/double(Nevents));
    //normalize(asymmetry, min_E, max_E);
    //for (int side=0; side<2; side++)
    //    for (int spin=0; spin<2; spin++)
    //        normalize(counts[side][spin], min_E, max_E);

    return Nsim;
}

/*
int UCNAmodel::fill(TString filename, TString name, TString title, 
                    double flip, double A, double b)
                    //double flip = 0.68/1.68, double A = -0.12, double b = 0)
{
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

	double Nevents = chain->GetEntries();
	chain->SetBranchStatus("*",false);
	chain->SetBranchStatus("PID",true);
	chain->SetBranchStatus("side",true);
	chain->SetBranchStatus("type",true);
	chain->SetBranchStatus("Erecon",true);
	chain->SetBranchStatus("primMomentum",true);
    chain->SetBranchStatus("ScintPosAdjusted",true);
    chain->SetBranchStatus("ScintPosAdjusted",true);

	//TNtuple* tntuple = new TNtuple("mc_ntuple", "MC NTuple", "s:load:energy");
	TNtuple* tntuple = 0;

	//int Noff=0, Non=0;	/// events have been simulated after cuts
	int Nspin[]={0,0};	    /// events have been simulated after cuts
    int PID, side, type;
    double energy;
    double mwpcPosW[3], mwpcPosE[3], p[3];
    //double mwpcPos[2][3], p[3]; 
    //double ex_cos = 4;  // TODO compute ex+cos correctly

    chain->SetBranchAddress("PID",&PID);
    chain->SetBranchAddress("side",&side);
    chain->SetBranchAddress("type",&type);
    chain->SetBranchAddress("Erecon",&energy);
	chain->SetBranchAddress("primMomentum",p); /// not momentum; normalized direction!
    chain->GetBranch("ScintPosAdjusted")->GetLeaf("ScintPosAdjE")->SetAddress(mwpcPosE);
    chain->GetBranch("ScintPosAdjusted")->GetLeaf("ScintPosAdjW")->SetAddress(mwpcPosW);

    for (int evt=0; evt<Nevents; evt++) {
        chain->GetEvent(evt);

        /// cut out bad events
        if (PID!=1) 
            continue;

        /// radial fiducial cut
        double fidcut2 = 50*50; // TODO !!!! VERY HARD CODED 
        assert(fidcut2 > 0);
        double radius2W = mwpcPosW[0]*mwpcPosW[0] + mwpcPosW[1]*mwpcPosW[1]; 
        double radius2E = mwpcPosE[0]*mwpcPosE[0] + mwpcPosE[1]*mwpcPosE[1]; 
        if (radius2W > fidcut2 or radius2E > fidcut2) 
            continue;

        /// Type 0, Type I, Type II/III events 
        if (type<4) { 
        //if (type==0) { 
            /// fill with loading efficiency 
            double load = rand.Uniform(1);
			double afp = (load < flip)? +1 : -1; // TODO don't need load as var..
            double axial = rand.Uniform(2)-1;
            double cos = p[2]/Sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
            assert(cos = p[2]);
            //bool spin = ((afp < 0) xor (cos/ex_cos < A)) ? EAST : WEST;
            double E = energy + m_e;        /// relativistic energy
            double x = m_e/E;
            bool spin = (afp < 0)? EAST : WEST;
            Nspin[spin]++;
            if (b > 0.5) {
                double d = 1 + b*x;
                double beta = Sqrt(1 - x*x);
                if (afp < 0 ) {
                    //Noff++;
                    if (axial > (1 + A*beta*cos/d))
                        spin = EAST;
                    else 
                        spin = WEST;
                }
                else {
                    //Non++;
                    if (axial > (1 - A*beta*cos/d))
                        spin = WEST;
                    else
                        spin = EAST;
                }
            }
            //cout<<"energy: "<<energy<<" side: "<<side<<" spin: "<<spin 
                <<" afp: "<<afp<<" p: ("<<p[0]<<","<<p[1]<<","<<p[2]<<")" 
                <<" cos: "<<cos<<" load: "<<load<<".\n";
            counts[side][spin]->Fill(energy, 1);
            if (tntuple)
			    tntuple->Fill(side, spin, energy);
        }
            hEreconALL->Fill(Erecon);
        if (type==0) 
            hErecon0->Fill(Erecon);
        else if (type==1) 
            hErecon1->Fill(Erecon);
        else if (type==2 or type==3) 
            hErecon23->Fill(Erecon); 

		/// break when enough data has been generated.
    }    
     
    int Nsim = Nspin[0] + Nspin[0];
	cout<<"Total number of Monte Carlo entries:\n";
	cout<<"      Entries without cuts:      "<<int(Nevents)<<endl;
	cout<<"      Entries with cuts:         "<<Nsim<<endl;
	cout<<"      Efficiencies of cuts:      "<<(100.0*Nsim)/double(Nevents)<<"%\n";
	cout<<"      Entries with AFP off cut:  "<<Nspin[0]<<" ("<<(100.0*Nspin[0])/double(Nsim)<<"%)\n";
	cout<<"      Entries with AFP on cut:   "<<Nspin[1]<<" ("<<(100.0*Nspin[1])/double(Nsim)<<"%)\n";

	/// compute and normalize super sum
    compute_super_sum();
    compute_asymmetry();
    //compute_super_diff();
    super_sum.normalize();
    //super_sum.Scale(1/Nsim);
    asymmetry.Scale(1/A);
    //super_sum.Scale(1/double(Nevents));
    //super_sum.Scale(double(nSimmed)/double(Nevents));
    //normalize(asymmetry, min_E, max_E);
    //for (int side=0; side<2; side++)
    //    for (int spin=0; spin<2; spin++)
    //        normalize(counts[side][spin], min_E, max_E);

    return Nsim;
}
*/


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
    double min = GetXaxis()->GetXmin();
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
    double max = GetXaxis()->GetXmax();
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
    double min = GetXaxis()->GetXmin();
    double max = GetXaxis()->GetXmax();
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

    if (upper < lower) {
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
            double count = n[side][spin];
            double sigma = e[side][spin];
            if (count > 0) {
                double alter = n[side?0:1][spin?0:1];   
                if (sigma > 0) {
                    double var = sigma*sigma;
                    varS += var*alter/count;  /// errors already set
                }
                else
                    varS += alter;             /// using Poisson stats
            }        /*
            } else {
                cout<<"Warning: Negative number in super sum error.\n";
                cout<<"         Counts for side "<<side<<" and spin "<<spin
                    <<" are "<<count<<" and sigma is "<<sigma<<".\n";
                sigma = 0.3;
                cout<<"         Setting to "<<sigma<<".\n";
            }
            */
        }
    }

    /// The 1/4 factor comes from the definition of S=(Y0+Y1)/2
    sigmaS = Sqrt(varS)/4;
    if (IsNaN(sigmaS)) {
        cout<<"Warning: Super sum error is not a number.\n";
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

    /*
    if (bin % 10 == 0) {
        double KE = super_sum.GetBinCenter(bin);
        cout<<"Status: Super sum bin "<<bin<<" and KE "<<KE<<" to "<<S<<"("<<sigmaS<<").\n";
    }
    */

    return S;
}

void UCNAhistogram::snapshot(int every)
{
    double min = GetXaxis()->GetXmin();
    double max = GetXaxis()->GetXmax();
    int bin_min = FindBin(min);
    int bin_max = FindBin(max);
    TString name = GetName();
    TString title = GetTitle();
    cout<<" HISTOGRAM SNAPSHOT:\n"
        <<"    Name: "<<name<<"   Title: "<<title<<"\n";
    for (int bin = bin_min; bin < bin_max; bin += every) {
        double KE = GetBinCenter(bin);
        double H = GetBinContent(bin);
        double sigmaH = GetBinError(bin);
        cout<<"     bin: "<<bin<<" KE: "<<KE<<" is "<<H<<"("<<sigmaH<<").\n";
    }
    cout<<"\n FULL HISTOGRAM:\n";
    for (int bin = bin_min; bin < bin_max; bin++) {
        double KE = GetBinCenter(bin);
        cout<<"\t"<<KE;
    }
    cout<<"\n";
    for (int bin = bin_min; bin < bin_max; bin++) {
        double H = GetBinContent(bin);
        cout<<"\t"<<H;
    }
    cout<<"\n";
    for (int bin = bin_min; bin < bin_max; bin++) {
        double sigmaH = GetBinError(bin);
        cout<<"\t"<<sigmaH;
    }
    cout<<"\n";
}

double UCNAmodel::compute_super_sum(int bin) {
    double count, error;
    return compute_super_sum(bin,count,error);
}

TH1D& UCNAmodel::compute_super_sum() {
    int bin_min, bin_max;
    return compute_super_sum(min, max, bin_min, bin_max);
}

TH1D& UCNAmodel::compute_super_sum(double min, double max, 
                                   int& bin_min, int& bin_max) {
    bin_min = super_sum.FindBin(min);
    bin_max = super_sum.FindBin(max);
    if (not super_sum.test_range(min, max)) {
        cout<<"Error: Problem with ranges in super sum histogram.\n";
        exit(1);
    }

/*  double lower = super_sum.GetXaxis()->GetBinLowEdge(bin_max);
    double bin_max_epsi = super_sum.FindBin(max+0.0001);
    cout<<"max="<<max<<" bin_max="<<bin_max<<" lower edge="<<lower<<" epsi="<<bin_max_epsi<<".\n";
    exit(1); */

    return compute_super_sum(bin_min, bin_max);
}

TH1D& UCNAmodel::compute_super_sum(int bin_min, int bin_max) 
{
    for (int bin = bin_min; bin < bin_max; bin++)
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
        if (S == 0)
            cout<<"         Division by zero super sum.\n"; 
        cout<<"         Setting asymmetry to zero.\n";

    }
    return A;
}

double UCNAmodel::compute_asymmetry(double n[2][2], double e[2][2], double& A, double& sigmaA) {
    //cout<<"Error: Called compute_asymmetry()\n";
    //exit(1);
    double S, sigmaS;
    compute_super_sum(n,e,S,sigmaS);
    A = compute_asymmetry(n);
    sigmaA = Sqrt(1 + A*A) * sigmaS/S;
    if (IsNaN(sigmaA)) {
        cout<<"Warning: Asymmetry error is not a number.\n"
            <<"         Setting asymmetry error to 1.\n";
        sigmaA = 1;
    }
    return A;
}

double UCNAmodel::compute_asymmetry(int bin, double& A, double& sigmaA) {
    double n[2][2], e[2][2];
    get_counts(bin,n,e);
    compute_asymmetry(n,e,A,sigmaA);
    if (A != 0 and (IsNaN(A) or IsNaN(sigmaA) or sigmaA < 0)) {
        double KE = asymmetry.GetBinCenter(bin);
        cout<<"Warning: Asymmetry error in bin "<<bin<<" and KE "<<KE<<" where\n         ";
        for (int side=0; side<2; side++) {
            for (int spin=0; spin<2; spin++) {
                double count = n[side][spin];
                cout<<"n["<<side<<","<<spin<<"] = "<<count<<";    ";
            }
        }
        cout<<"\n";
        //cout<<"\n         Skipping this data point.\n";
    }
    else {
        asymmetry.SetBinContent(bin,A);
        asymmetry.SetBinError(bin,sigmaA);
    }

    return A;
}


double UCNAmodel::compute_asymmetry(int bin) {
    double count, error;
    return compute_asymmetry(bin,count,error);
}


TH1D& UCNAmodel::compute_asymmetry() {
    int bin_min, bin_max;
    return compute_asymmetry(min, max, bin_min, bin_max);
}


TH1D& UCNAmodel::compute_asymmetry(double min, double max, 
                                   int& bin_min, int& bin_max) {
    bin_min = asymmetry.FindBin(min);
    bin_max = asymmetry.FindBin(max);
    if (not asymmetry.test_range(min,max)) {
        cout<<"Error: Problem with ranges in asymmetry histogram.\n";
        cout<<"       Computing asymmetry in bin range "<<bin_min<<"-"<<bin_max<<"\n";
        cout<<"       and energy range "<<min<<" keV - "<<max<<" keV.\n";
        exit(1);
    }
    return compute_asymmetry(bin_min,bin_max);
}

TH1D& UCNAmodel::compute_asymmetry(int bin_min, int bin_max) 
{
    for (int bin=bin_min; bin<bin_max; bin++)
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

    for (int bin=bin_min; bin<bin_max; bin++)
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

double UCNAhistogram::GetEffectiveEntries(double min, double max)
{
    /*
	double eff_ent = GetEffectiveEntries();
    double bin_min = FindBin(min);
	double bin_max = FindBin(max);
	double part_int = Integral(bin_min, bin_max);
	double full_int = Integral();
	double N = eff_ent * part_int / full_int;

	return N;*/
	double Neff = 0;
    double bin_min = FindBin(min);
	double bin_max = FindBin(max);
	for (int bin=bin_min; bin<bin_max; bin++)
	{
		double Y = GetBinContent(bin);
		double eY = GetBinError(bin);
        if (eY > 0) {
            double sqrtn = Y/eY;
            Neff += sqrtn*sqrtn;
        }
	}
    return Neff;
}



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




/*
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
*/


/// relativistic beta factor
double beta_factor(const double KE)
{
	if (KE <= 0 or KE >= Q)             ///< kinetic energy
		return 0;                       ///< zero beyond endpoint

	double E = KE + m_e;                ///< electron energy
	double x = m_e / E;                 ///< reduced electron energy
	double B = sqrt(1 - x*x);  	        ///< beta power factor

	return B;
}


/// relativistic beta factor
/*
double cos_theta(const double p[3])
{
	if (KE <= 0 or KE >= Q)                     ///< kinetic energy
		return 0;                               ///< zero beyond endpoint

	const double E = KE + m_e;                  ///< electron energy
	const double x = m_e / E;                   ///< reduced electron energy
	const double B = sqrt(1 - x*x);  	        ///< beta power factor

	return B;
}
*/


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


double evaluate_expected_fierz(double _min, double _max) {
    return evaluate_expected_fierz(0,1,_min,_max);
}

double evaluate_expected_fierz(double ym, double xn, double _min, double _max)
{
    int integral_size = 5512;
    // TODO These don't need to be pointers.
    /*
    TH1D *h1 = new TH1D("beta_spectrum_fierz", 
                        "Beta spectrum with Fierz term", 
                        integral_size, _min, _max);
    TH1D *h2 = new TH1D("beta_spectrum", 
                        "Beta Spectrum", 
                        integral_size, _min, _max);
	for (int i = 1; i < integral_size; i++)
	{
		double K = _min + double(i)*(_max-_min)/integral_size;
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
    */
    TH1D hmn("beta_spectrum_fierz", 
             "Beta spectrum with Fierz term", 
             integral_size, _min, _max);
    TH1D h10("beta_spectrum", 
             "Beta Spectrum", 
             integral_size, _min, _max);
    double parmn[2] = {ym, xn};
    double par10[2] = {0, 0};
	for (int bin = 1; bin < integral_size; bin++)
	{
		//double K = _min + double(i)*(_max-_min)/integral_size;
		double KE[1] = {hmn.GetBinCenter(bin)};
		double ymn = beta_spectrum(KE, parmn);
		double y10 = beta_spectrum(KE, par10);
		hmn.SetBinContent(bin, ymn);
		h10.SetBinContent(bin, y10);
	}
	 return hmn.Integral(0, integral_size) / h10.Integral(0, integral_size);
}


/*
double evaluate_expected_fierz(double _min, double _max)
{
	return evaluate_expected_fierz(0, 1, _min, _max);
}
*/


TH1D* compute_rate_function(TH1D* rate_histogram[2][2], 
                            double (*rate_function)(double r[2][2]));
TH1D* compute_rate_function(TH1D* rate_histogram[2][2], 
                            double (*rate_function)(double r[2][2]),
                            double (*error_function)(double r[2][2]));

/*
TH1D* compute_rate_error_function(TH1D* rate_histogram[2][2], 
                                  double (*rate_error_function)(double r[2][2])) 
{
    TH1D *out_histogram = new TH1D(*(rate_histogram[0][0]));
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
double bonehead_sum(double r[2][2]);
double bonehead_asymmetry(double r[2][2]);
double super_ratio_asymmetry(double r[2][2]);


/*
double super_sum_error(double r[2][2]) {
    double super_ratio = (r[0][0] * r[1][1]) / (r[0][1] * r[1][0]);
    double sqrt_super_ratio = Sqrt(super_ratio);
    if ( IsNaN(sqrt_super_ratio) ) 
        sqrt_super_ratio = 0;
    //return (1 - sqrt_super_ratio) / (1 + sqrt_super_ratio);
}
*/
