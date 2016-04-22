#include "FierzFitter.hh"


double UCNAhistogram::normalize() 
{
	int bin_min = histogram->FindBin(min);
	int bin_max = histogram->FindBin(max);
    double integrand = 1;
    integrand +=(histogram->GetBinCenter(bin_min) - min)
               / histogram->GetBinWidth(bin_min)
               / histogram->GetBinWidth(bin_min)
               * histogram->GetBinContent(bin_min)
               +(max - histogram->GetBinCenter(bin_max))
               / histogram->GetBinWidth(bin_max)
               / histogram->GetBinWidth(bin_max)
               * histogram->GetBinContent(bin_max);

    for (int bin=bin_min+1; bin<bin_max; bin++)
        integrand += histogram->GetBinContent(bin)
                   / histogram->GetBinWidth(bin);

	histogram->Scale(1/integrand);
    return integrand;
}


//TF1* UCNAFierzFitter::combined_fit(TH1D* asymmetry, TH1D* super_sum, TMatrixD &cov, TF1* func)
TF1* UCNAFierzFitter::combined_fit(
TMatrixD &cov, TF1 *func, 
        void (*global_fcn_ptr)(Int_t&, Double_t*, Double_t&, Double_t*, Int_t) )
{ 
    int nPar = func->GetNpar();
    if (not func) {
        cout<<"Fit function not set.\n";
        exit(1);
    }

	/// set up the minuit fitter
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

	/// extract results from minuit
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

    TH1D* asymmetry = data.asymmetry.histogram; 
    TH1D* super_sum = data.super_sum.histogram; 
	int ndf = asymmetry->GetNbinsX() + super_sum->GetNbinsX() - nvpar;
	func->SetNDF(ndf);
    
	TMatrixD matrix( nPar, nPar, minuit->GetCovarianceMatrix() );
	for (int i=0; i<nPar; i++)
		for (int j=0; j<nPar; j++)
			cov[i][j] = minuit->GetCovarianceMatrixElement(i,j);
	
	cout<<"\n    chi^2 = "<<chi2<<", ndf = "<<ndf<<", chi^2/ndf = "<<chi2/ndf<<".\n\n";

	return func; 
}


double UCNAFierzFitter::asymmetry_chi2(double A, double b)
{
	double chi2 = 0;
    int n = data.asymmetry.bins;
	for (int i = 0; i < n; i++)
	{
		double E = data.asymmetry.histogram->GetBinCenter(i);
        if (E < min or E > max)
            continue;
		double Y = data.asymmetry.histogram->GetBinContent(i);
		double eY = data.asymmetry.histogram->GetBinError(i);
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
    int n = data.asymmetry.bins;
	for (int i = 0; i < n; i++)
	{
		double E = data.asymmetry.histogram->GetBinCenter(i);
        if (E < min or E > max)
            continue;
		double Y = data.asymmetry.histogram->GetBinContent(i);
		double eY = data.asymmetry.histogram->GetBinError(i);
        double f  = N*sm.super_sum.histogram->GetBinContent(i) 
                  + N*b*fierz.super_sum.histogram->GetBinContent(i);
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
		double E = data.asymmetry.histogram->GetBinCenter(i);
		double Y = data.asymmetry.histogram->GetBinContent(i);
        //double f = asymmetry_fit_func(&E,p);
        double f = A/(1 + b*m_e/(E+m_e));
		//double eY = data.asymmetry.histogram->GetBinError(i);
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
		//double E      = data .super_sum.histogram->GetBinCenter(i);
		//chi = (fierzratio_values[i] - fierzratio_fit_func(&E,par)) / fierzratio_errors[i];
		double Y  = data .super_sum.histogram->GetBinContent(i);
        double eY = data .super_sum.histogram->GetBinError(i);
        double f  = N*sm .super_sum.histogram->GetBinContent(i) 
                  + N*b*fierz.super_sum.histogram->GetBinContent(i);
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

/*
void normalize(double min, double max) 
{
	int bin_min = hist->FindBin(min);
	int bin_max = hist->FindBin(max);
	hist->Scale(1/(hist->GetBinWidth(2)*hist->Integral(_min, _max)));
}
*/


double evaluate_expected_fierz(double min, double max, double integral_size) 
{
    TH1D *h1 = new TH1D("beta_spectrum_fierz", "Beta spectrum with Fierz term", integral_size, min, max);
    TH1D *h2 = new TH1D("beta_spectrum", "Beta Spectrum", integral_size, min, max);
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

#if 0

/**
 * Si := (r[0][0] r[1][1]) / (r[0][1] * r[1][0])
 */
TH1D* compute_super_ratio(TH1D* rate_histogram[2][2] ) 
{
    TH1D *super_ratio_histogram = new TH1D(*(rate_histogram[0][0]));
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
TH1D* compute_super_sum(TH1D* histogram[2][2], TH1D* super_sum_histogram) 
{
    //TH1D *super_sum_histogram = new TH1D(*(rate_histogram[0][0]));
    for (int side=0; side<2; side++)
        for (int spin=0; spin<2; spin++)
            if (not histogram[side][spin]) {
                std::cout << "Error: histogram is not constructed.\n";
                std::cout << "Side: "<< side << " Spin: " << spin <<".\n";
                exit(1);
            }

    if (not super_sum_histogram) {
        std::cout << "Error: super sum histogram is not constructed.\n";
        exit(1);
    }

    int bins = super_sum_histogram->GetNbinsX();
    for (int side = 0; side < 2; side++)
        for (int spin = 0; spin < 2; spin++) {
            int ss_bins = histogram[side][spin]->GetNbinsX();
            if (bins != ss_bins and ss_bins == 0) {
                std::cout << "Error: super sum and side spin histogram sizes don't match.\n";
                std::cout << "super sum bins: " << bins << "\n";
                exit(1);
            }
        }

    for (int bin = 1; bin < bins+2; bin++) {
        double r[2][2];
        for (int side = 0; side < 2; side++)
            for (int spin = 0; spin < 2; spin++)
                r[side][spin] = histogram[side][spin]->GetBinContent(bin);

        double super_sum = TMath::Sqrt(r[0][0] * r[1][1]) + TMath::Sqrt(r[0][1] * r[1][0]);
        double rel_error = TMath::Sqrt( 1/r[0][0] + 1/r[1][0] + 1/r[1][1] + 1/r[0][1]);
        if ( TMath::IsNaN(super_sum)) {
            super_sum = 0;
            std::cout << "Warning: super sum: division by zero.\n";
        }

        if (TMath::IsNaN(rel_error)) {
			rel_error = 0;
            std::cout << "Warning: super sum error: division by zero.\n";
        }
        
        if (bin % 10 == 0)
            printf("Setting bin content for super sum bin %d, to %f\n", bin, super_sum);

        super_sum_histogram->SetBinContent(bin, super_sum);
        super_sum_histogram->SetBinError(bin, super_sum*rel_error);
    }
    return super_sum_histogram;
}


TH1D* compute_asymmetry(TH1D* rate_histogram[2][2]) 
{
    TH1D *asymmetry_histogram = new TH1D(*(rate_histogram[0][0]));
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


TH1D* compute_corrected_asymmetry(TH1D* rate_histogram[2][2], TH1D* correction) 
{
    TH1D *asymmetry_histogram = new TH1D(*(rate_histogram[0][0]));
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
        asymmetry_histogram->Multiply(correction);
        //printf("Setting bin content for corrected asymmetry bin %d, to %f\n", bin, asymmetry);
    }
    return asymmetry_histogram;
}

#endif
