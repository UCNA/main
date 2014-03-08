/// 

#include "BetaSpectrum.hh"
#include <stdio.h>
#include <TMath.h>
#include <vector>
#include <map>
#include "SMExcept.hh"

/// hyperbolic sine
double my_sinh(double x) { return (exp(x)-exp(-x))*0.5; }

/// inverse hyperbolic tangent
//double atanh(double x) throw() { return 0.5*log((1.+x)/(1.-x)); }

/// struct for 1-index coefficients
struct coeff1 {
	coeff1(int ii, double cc): i(ii), c(cc) {}
	int i;
	double c;
};

/// sum coeff3 power series
double sumCoeffs(const std::vector<coeff1>& coeffs, double x=1.0) {
	double s = 0;
	for(std::vector<coeff1>::const_iterator it = coeffs.begin(); it != coeffs.end(); it++)
		s += it->c*pow(x,it->i);
	return s;
}

/// struct for 3-index coefficients
struct coeff2 {
	coeff2(int ii, int jj, double cc): i(ii), j(jj), c(cc) {}
	int i;
	int j;
	double c; 
};

/// sum coeff3 power series
double sumCoeffs(const std::vector<coeff2>& coeffs, double x=1.0, double y=1.0) {
	double s = 0;
	for(std::vector<coeff2>::const_iterator it = coeffs.begin(); it != coeffs.end(); it++)
		s += it->c*pow(x,it->i)*pow(y,it->j);
	return s;
}

/// struct for 3-index coefficients
struct coeff3 {
	coeff3(int ii, int jj, int kk, double cc): i(ii), j(jj), k(kk), c(cc) {}
	int i;
	int j;
	int k;
	double c; 
};

/// sum coeff3 power series
double sumCoeffs(const std::vector<coeff3>& coeffs, double x=1.0, double y=1.0, double z=1.0) {
	double s = 0;
	for(std::vector<coeff3>::const_iterator it = coeffs.begin(); it != coeffs.end(); it++)
		s += it->c*pow(x,it->i)*pow(y,it->j)*pow(z,it->k);
	return s;
}

/// definition of \gamma in [1]
double WilkinsonGamma(double Z = 1.) { return sqrt(1-(alpha*Z)*(alpha*Z)); }

double WilkinsonF_PowerSeries(double Z, double W, double R) {
	// set up coeffs
	static std::vector<coeff3> coeffs;
	if(!coeffs.size()) {
		coeffs.push_back(coeff3(0,0,0,1.));
		
		coeffs.push_back(coeff3(1,1,0,M_PI));
		
		coeffs.push_back(coeff3(2,0,0,0.577216));
		coeffs.push_back(coeff3(2,0,1,-1.));
		coeffs.push_back(coeff3(2,2,0,3.289868));
		
		coeffs.push_back(coeff3(3,1,0,1.813376));
		coeffs.push_back(coeff3(3,1,1,-M_PI));
	}
	
	double p = sqrt(W*W-1);
	double gm = WilkinsonGamma(Z);
	double z = TMath::Gamma(2.*gm+1.);
	return 2.*(gm+1.)/(z*z)*pow(2*R,2*(gm-1))*sumCoeffs(coeffs,alpha*Z,W/p,log(p));
}

/// approximation to |Gamma(gm+iaZW/p)|^2, as per [3] eq. 1
double WilkinsonGammaMagSquaredApprox(double Z, double W, unsigned int N) {
	double s = 0;
	double gm = WilkinsonGamma(Z);
	double y = alpha*Z*W/sqrt(W*W-1);
	double Ngm = N+gm;
	double a = (N+1.)/Ngm;
	double y1 = a*y;
	for(unsigned int n=0; n<N; n++)
		s += log((n*n+y1*y1)/((n+gm)*(n+gm)+y*y));
	return exp(s + log(M_PI*(N*N+y1*y1)/(y1*my_sinh(M_PI*y1)))
			   +(1.-gm)*(2.-log(Ngm*Ngm+y*y) + 2.*y/Ngm*atan(y/Ngm)
						 +1./(Ngm*Ngm+y*y)/(6.*a))
			   -(2.*N+1.)*log(a));	
}

double F_approx(double W, double Z, double R, bool fullterms) {
	if(W<=1) return 1;
	double F = 1;
	double p = sqrt(W*W-1.);
	double azn = alpha*Z;
	double x = 11./4.-gamma_euler-log(2*p*R);
	F += azn*(M_PI*W/p);
	if(!fullterms) return F;
	azn *= alpha*Z;
	F += azn*(x+M_PI*M_PI*W*W/(3.*p*p));
	azn *= alpha*Z;
	F += azn*(M_PI*W/p*x);
	return F<20?F:20;
}

// F0 = 2/(1+WilkinsonGamma) * F
double WilkinsonF0(double Z, double W, double R, unsigned int N) {
	if(W<=1) return 0;
	double gm = WilkinsonGamma(Z);
	double GMi = 1./TMath::Gamma(2*gm+1);
	double p = sqrt(W*W-1.);
	double F0 = 4.*pow(2.*p*R,2.*gm-2.)*GMi*GMi*exp(M_PI*Z*alpha*W/p)*WilkinsonGammaMagSquaredApprox(Z,W,N);
	return F0<1e3?F0:0;
}

double WilkinsonRV(double W, double W0, double M) {
	return (1.+W0*W0/(2.*M*M)-11./(6.*M*M)
			+W0/(3*M*M)/W
			+(2./M-4.*W0/(3.*M*M))*W
			+16./(3.*M*M)*W*W);
}

double WilkinsonRA(double W, double W0, double M) {
	return (1.+2.*W0/(3.*M)-W0*W0/(6.*M*M)-77./(18.*M*M)
			+(-2./(3.*M)+7*W0/(9.*M*M))/W
			+(10./(3.*M)-28.*W0/(9.*M*M))*W
			+88./(9.*M*M)*W*W);
}

double CombinedR(double W, double M2_F, double M2_GT, double W0, double M) {
	return (M2_F*WilkinsonRV(W,W0,M)+lambda*lambda*M2_GT*WilkinsonRA(W,W0,M))/(M2_F+lambda*lambda*M2_GT);
}

double Bilenkii59_RWM(double W) {
	const double mu = 2.792847356-(-1.91304273);
	return (-2.*lambda*(lambda+mu)*beta_W0
			+2.*(5.*lambda*lambda+2.*lambda*mu+1.)*W
			-2.*lambda*(mu+lambda)/W)/(1.+3*lambda*lambda)/proton_M0;
}

double WilkinsonL0(double Z, double W, double R) {
	// set up coeffs
	static std::vector<coeff1> ai[6];
	static std::vector<coeff1> aminus1;
	static std::map<double,std::vector<coeff1> > aiZ;
	static std::map<double,double> aminus1Z;
	
	if(!ai[0].size()) {
		aminus1.push_back(coeff1(1,0.115));
		aminus1.push_back(coeff1(2,-1.8123));
		aminus1.push_back(coeff1(3,8.2498));
		aminus1.push_back(coeff1(4,-11.223));
		aminus1.push_back(coeff1(5,-14.854));
		aminus1.push_back(coeff1(6,32.086));
		
		ai[0].push_back(coeff1(1,-0.00062));
		ai[1].push_back(coeff1(1,0.02482));
		ai[2].push_back(coeff1(1,-0.14038));
		ai[3].push_back(coeff1(1,0.008152));
		ai[4].push_back(coeff1(1,1.2145));
		ai[5].push_back(coeff1(1,-1.5632));
		
		ai[0].push_back(coeff1(2,0.007165));
		ai[1].push_back(coeff1(2,-0.5975));
		ai[2].push_back(coeff1(2,3.64953));
		ai[3].push_back(coeff1(2,-1.15664));
		ai[4].push_back(coeff1(2,-23.9931));
		ai[5].push_back(coeff1(2,33.4192));
		
		ai[0].push_back(coeff1(3,0.01841));
		ai[1].push_back(coeff1(3,4.84199));
		ai[2].push_back(coeff1(3,-38.8143));
		ai[3].push_back(coeff1(3,49.9663));
		ai[4].push_back(coeff1(3,149.9718));
		ai[5].push_back(coeff1(3,-255.1333));
		
		ai[0].push_back(coeff1(4,-0.53736));
		ai[1].push_back(coeff1(4,-15.3374));
		ai[2].push_back(coeff1(4,172.1368));
		ai[3].push_back(coeff1(4,-273.711));
		ai[4].push_back(coeff1(4,-471.2985));
		ai[5].push_back(coeff1(4,938.5297));
		
		ai[0].push_back(coeff1(5,1.2691));
		ai[1].push_back(coeff1(5,23.9774));
		ai[2].push_back(coeff1(5,-346.708));
		ai[3].push_back(coeff1(5,657.6292));
		ai[4].push_back(coeff1(5,662.1909));
		ai[5].push_back(coeff1(5,-1641.2845));
		
		ai[0].push_back(coeff1(6,-1.5467));
		ai[1].push_back(coeff1(6,-12.6534));
		ai[2].push_back(coeff1(6,288.7873));
		ai[3].push_back(coeff1(6,-603.7033));
		ai[4].push_back(coeff1(6,-305.6804));
		ai[5].push_back(coeff1(6,1095.358));
	}
	
	if(!aiZ.count(Z)) {
		std::vector<coeff1> aiZi;
		for(unsigned int i=0; i<6; i++)
			aiZi.push_back(coeff1(i,sumCoeffs(ai[i],alpha*Z)));
		aiZ.insert(std::make_pair(Z,aiZi));
		aminus1Z.insert(std::make_pair(Z,sumCoeffs(aminus1,alpha*Z)));
	}
	
	if(W<=1)
		return 0;
	
	double gm = WilkinsonGamma(Z);
	double L0 = (1.+13.*(alpha*Z)*(alpha*Z)/60.-W*R*alpha*Z*(41.-26.*gm)/(15.*(2.*gm-1.))
				 -alpha*Z*R*gm*(17.-2.*gm)/(30.*W*(2.*gm-1.))
				 +aminus1Z[Z]*R/W+sumCoeffs(aiZ[Z],W*R)
				 +0.41*(R-0.0164)*pow(alpha*Z,4.5) );
	
	return L0==L0?L0*2./(1.+gm):0;
}

double WilkinsonVC(double Z, double W, double W0, double R) {
	double gm = WilkinsonGamma(Z);
	return (1.-233.*alpha*Z*alpha*Z/630.-W0*R*W0*R/5.-6*W0*R*alpha*Z/35.
			+(-13.*R*alpha*Z/35.+4.*W0*R*R/15.)*W
			+(2.*gm*W0*R*R/15.+gm*R*alpha*Z/70.)/W
			-4*R*R/15.*W*W);
}

double WilkinsonAC(double Z, double W, double W0, double R) {
	return (1.-233.*alpha*Z*alpha*Z/630.-W0*R*W0*R/5.+2*W0*R*alpha*Z/35.
			+(-21.*R*alpha*Z/35.+4.*W0*R*R/9.)*W
			-4.*R*R/9.*W*W);
}

double CombinedC(double Z, double W, double M2_F, double M2_GT, double W0, double R) {
	return (M2_F*WilkinsonVC(Z,W,W0,R)+lambda*lambda*M2_GT*WilkinsonAC(Z,W,W0,R))/(M2_F+lambda*lambda*M2_GT);
}

double WilkinsonQ(double, double W, double W0, double M) {
	double B = (1.-lambda)/(1.+3.*lambda*lambda);
	return 1.-M_PI*alpha/(M*sqrt(W*W-1))*(1+B*(W0-W)/(3.*W));
}



double SpenceL(double x, unsigned int N=20) {
	smassert(-1.<x && x<=1.);
	double s = 0;
	double xk = x;
	for(unsigned int k=1; k<=N; k++) {
		s += xk/(k*k);
		xk *= x;
	}
	return -s;
}

double Sirlin_g_a2pi(double KE,double KE0,double m) {
	if(KE<=0 || KE>=KE0)
		return 0;
	double b = beta(KE,m);
	double E = KE+m;
	double E0 = KE0+m;
	double athb = atanh(b);
	return (3.*log(m_p/m)-3./4.
			+4.*(athb/b-1.)*((E0-E)/(3.*E)-3./2.+log(2.*(E0-E)/m))
			+4./b*SpenceL(2.*b/(1.+b))
			+athb/b*(2.*(1.+b*b)+(E0-E)*(E0-E)/(6.*E*E)-4.*athb)
			)*alpha/(2.*M_PI);
}

double shann_h_a2pi(double KE, double KE0, double m) {
	if(KE<=0 || KE>=KE0)
		return 0;
	double b = beta(KE,m);
	double E = KE+m;
	double E0 = KE0+m;
	double athb = atanh(b);
	return (3.*log(m_p/m)-3./4.
			+ 4.*(athb/b-1.)*((E0-E)/(3.*E*b*b)+(E0-E)*(E0-E)/(24*E*E*b*b)-3./2.+log(2.*(E0-E)/m))
			+ 4./b*SpenceL(2.*b/(1.+b))+4*athb/b*(1-athb)
			)*alpha/(2.*M_PI);
}

double Wilkinson_g_a2pi(double W, double W0, double M) {
	if(W>=W0 || W<=1)
		return 0;
	double b = sqrt(W*W-1)/W;
	double athb = atanh(b);
	double g = (3.*log(M)-3./4.
				+4.*(athb/b-1.)*((W0-W)/(3.*W)-3./2.+log(2))
				+4./b*SpenceL(2.*b/(1.+b))
				+athb/b*(2.*(1.+b*b)+(W0-W)*(W0-W)/(6.*W*W)-4.*athb)
				)*alpha/(2.*M_PI)+pow((W0-W),2*alpha/M_PI*(athb/b-1.))-1.;
	return g==g?g:0;
}


double neutronSpectrumCorrectionFactor(double KE) {
	double W = (KE+m_e)/m_e;
	double c = WilkinsonF0(1,W,neutron_R0);				// Fermi function Coulomb
	c *= WilkinsonL0(1,W,neutron_R0);					// Nonzero charge radius
	c *= CombinedC(1,W,1.,3.,beta_W0,neutron_R0);		// electron/nucleon nonzero size wavefunction convolution
	c *= WilkinsonQ(1,W,beta_W0,proton_M0);				// Coulomb effect on recoil
	c *= (1.+Wilkinson_g_a2pi(W,beta_W0));				// outer radiative corrections
	c *= (1.+Bilenkii59_RWM(W));						// recoil + WM for free neutron
	return c;
}

double neutronCorrectedBetaSpectrum(double KE) {
	double W = (KE+m_e)/m_e;
	return plainPhaseSpace(W,beta_W0)*neutronSpectrumCorrectionFactor(KE);
}

double Davidson_C1T(double W, double W0, double Z, double R) {
	double p = sqrt(W*W-1);
	double y = alpha*Z*W/p;
	double a2Z2 = alpha*alpha*Z*Z;
	double S0 = sqrt(1-a2Z2);
	double S1 = sqrt(4-a2Z2);
	const double C = pow(TMath::Gamma(0.25),2)/sqrt(8*M_PI*M_PI*M_PI); // "Gauss number"
	double sm = 0;
	for(unsigned int n=1; n<10; n++) sm += 1/(n*(n*n+y*y));
	double A = ( (S1+2)/(2*S0+2) * pow(12*TMath::Gamma(2.*S0+1.)/TMath::Gamma(2.*S1+1.),2) *
			  pow(2*p*R,a2Z2/2) * (pow(1-a2Z2/4,2)+y*y) * (1-a2Z2*C/2+a2Z2*y*y*sm/2) );
	
	return (1+S0)*((W0-W)*(W0-W)+A*(W*W-1))/24;
}

double Langer_Cs137_C2T(double W, double W0) {
	const double k = 0.030;
	return (W0-W)*(W0-W)+k*(W*W-1);
}

double Behrens_l2(double W, double W0, double Z, double R) {
	double p = sqrt(W*W-1);
	double y = alpha*Z*W/p;
	double a2Z2 = alpha*alpha*Z*Z;
	double S0 = sqrt(1-a2Z2);
	double S1 = sqrt(4-a2Z2);
	const double C = pow(TMath::Gamma(0.25),2)/sqrt(8*M_PI*M_PI*M_PI); // "Gauss number"
	double sm = 0;
	for(unsigned int n=1; n<10; n++) sm += 1/(n*(n*n+y*y));
	return ( (S1+2)/(2*S0+2) * pow(12*TMath::Gamma(2.*S0+1.)/TMath::Gamma(2.*S1+1.),2) *
		   pow(2*p*R,a2Z2/2) * (pow(1-a2Z2/4,2)+y*y) * (1-a2Z2*C/2+a2Z2*y*y*sm/2) );
}

double Behrens_Cs137_C(double W, double W0) {
	double q2 = (W0-W)*(W0-W);
	double p2 = W*W-1;
	double l2 = Behrens_l2(W,W0,56,pow(137,1./3.)*neutron_R0);
	double a1 = 0.000346865*q2 + 0.00331725*l2*p2 - 0.000050327*q2*W + 0.000155636*l2*p2*W + 0.000114834*q2/W;
	double a2 = -0.00427141*q2 - 0.00645269*l2*p2 + 0.000063321*q2*W - 0.000913829*l2*p2*W - 0.000566409*q2/W + 0.0000576232*l2*p2/W;
	double a3 =   0.0131499*q2 + 0.00313793*l2*p2 + 0.00151806*q2*W  + 0.000741562*l2*p2*W - 0.000865957*q2/W - 0.000272219*l2*p2/W;
	double x = 1.07;
	return a1+a2*x+a3*x*x;
}

//-----------------------------------------------------//

BetaSpectrumGenerator::BetaSpectrumGenerator(double a, double z, double ep): A(a), Z(z), EP(ep),
W0((EP+m_e)/m_e), R(pow(A,1./3.)*neutron_R0), M0(fabs(Z)*proton_M0+(A-fabs(Z))*neutron_M0),
forbidden(0), M2_F(0), M2_GT(1) { }

double BetaSpectrumGenerator::spectrumCorrectionFactor(double W) const {
	double c = WilkinsonF0(Z,W,R);			// Fermi function Coulomb
	c *= WilkinsonL0(Z,W,R);				// Nonzero charge radius effect on Coulomb correction
	c *= CombinedC(Z,W,M2_F,M2_GT,W0,R);	// electron/nucleon nonzero size wavefunction convolution
	c *= WilkinsonQ(Z,W,W0,M0);				// Recoil effect on Coulomb correction (tiny tiny!)
	c *= (1.+Wilkinson_g_a2pi(W,W0,M0));	// outer radiative corrections
	if(A==1 && Z==1) {
		c *= (1.+Bilenkii59_RWM(W));		// recoil + WM for free neutron
	} else {
		c *= CombinedR(W,M2_F,M2_GT,W0,M0);	// recoil effect on phase space
	}
	// first forbidden Axial-Vector decays
	if(forbidden==1 && M2_GT>0 && M2_F==0)
		c *= Davidson_C1T(W, W0, Z, R);
	// Cs137 second-forbidden decay
	if(forbidden==2 && A==137)
		c *= Behrens_Cs137_C(W, W0);
	
	return c;
}

double BetaSpectrumGenerator::decayProb(double KE) const {
	double W = (KE+m_e)/m_e;
	return plainPhaseSpace(W,W0)*spectrumCorrectionFactor(W);
}

//-----------------------------------------------------//


double shann_h_minus_g_a2pi(double W, double W0) {
	if(W>=W0 || W<=1)
		return 0;
	double b = sqrt(W*W-1)/W;
	double athb = atanh(b);
	return ( 4.*(athb/b-1.)*(1/(b*b)-1)*(W0-W)/(3*W)*(1+(W0-W)/(8*W))
			+athb/b*(2.-2.*b*b) - (W0-W)*(W0-W)/(6.*W*W) )*alpha/(2.*M_PI);
}

double WilkinsonACorrection(double W) {
	const double W0 = neutronBetaEp/m_e+1.;
	const double mu = 2.792847356-(-1.91304273);	// mu_p - mu_n = 2.792847356(23) - -1.91304273(45) PDG 2010
	const double A_uM = (lambda+mu)/(lambda*(1.-lambda)*(1.+3.*lambda*lambda)*m_p/m_e);
	const double A_1 = lambda*lambda+2.*lambda/3.-1./3.;
	const double A_2 = -lambda*lambda*lambda-3.*lambda*lambda-5.*lambda/3.+1./3.;
	const double A_3 = 2.*lambda*lambda*(1.-lambda);
	return A_uM*(A_1*W0+A_2*W+A_3/W);
}
