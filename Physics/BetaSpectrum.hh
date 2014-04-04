/// \file BetaSpectrum.hh calculations for the beta spectrum and various related corrections
#ifndef BETASPECTRUM_HH
#define BETASPECTRUM_HH

// references:
// [0] Wilkinson, Analysis of Neutron Beta-Decay, Nucl. Phys. A 377 (1982) 474-504
// [1] Wilkinson, Evaluation of Beta-Decay I,   NIM A 275 (1989) 378-386
// [2] Wilkinson, Evaluation of Beta-Decay II,  NIM A 290 (1990) 509-515
// [3] Wilkinson, Evaluation of Beta-Decay III, NIM A 335 (1995) 305-309
// [4] Wilkinson, Evaluation of Beta-Decay IV,  NIM A 365 (1995) 203-207
// [5] Wilkinson, Evaluation of Beta-Decay V,   NIM A 365 (1995) 497-507
	
#include <math.h>
// useful physics constants; several in hbar=m_e=c=1 "natural units"
const double neutronBetaEp = 782.347;			//< neutron beta decay endpoint, keV
const double m_e = 511.00;						//< electron mass, keV/c^2
const double m_p = 938272.046;					//< proton mass, keV/c^2
const double m_n = m_p+m_e+neutronBetaEp;		//< neutron mass, keV/c^2
const double alpha = 1./137.036;				//< fine structure constant
const double lambda = fabs(-1.2694);			//< +/-0.0028, PDG 2010 value, Wilkinson sign convention
const double A0_PDG = -0.1173;					//< +/-0.0013, PDG 2010 value
const double beta_W0 = (neutronBetaEp+m_e)/m_e;	//< beta spectrum endpoint, ``natural'' units
const double neutron_R0 = 0.0025896*1.2;		//< neutron and proton radius approximation, in "natural" units (1.2fm)/(hbar/m_e*c)
const double proton_M0 = m_p/m_e;				//< proton mass, ``natural'' units
const double neutron_M0 = m_n/m_e;				//< neutron mass, ``natural'' units
const double gamma_euler = 0.577215;			//< Euler's constant

// NOTE: functions of W are using Wilkinson's ``natural'' units for energy, W=(KE+m_e)/m_e

//-------------- Spectrum corrections ------------------

/// beta decay phase space without corrections
inline double plainPhaseSpace(double W, double W0=beta_W0) { return (1.<W && W<W0)?sqrt(W*W-1)*W*(W0-W)*(W0-W):0; }
/// beta for particle with given KE
inline double beta(double KE, double m = m_e) { return sqrt(KE*KE+2*m*KE)/(m+KE); }

/// lowest order approximation of F
inline double crudeF(double Z, double W) { return 1+M_PI*alpha*Z*W/sqrt(W*W-1.); }
/// power series approximation of F(Z,W;R) in [1]
double WilkinsonF_PowerSeries(double Z, double W, double R=neutron_R0);
/// Wilkinson's F0(Z,W;R) as in [0],[1],[2],[3]; using complex gamma approximation to N terms
double WilkinsonF0(double Z, double W, double R = neutron_R0, unsigned int N = 3);

/// L_0(Z,W) as parametrized in [2], correction to point-like point charge used for F0(Z,W)
double WilkinsonL0(double Z, double W, double R = neutron_R0);

/// R(W,W0,M) as parametrized in [2], phase space correction for nuclear recoil, Vector part
double WilkinsonRV(double W, double W0=beta_W0, double M=proton_M0);
/// R(W,W0,M) as parametrized in [2], phase space correction for nuclear recoil, Axial Vector part
double WilkinsonRA(double W, double W0=beta_W0, double M=proton_M0);
/// Combined Vector/Axial-Vector nuclear recoil correction to spectrum shape
double CombinedR(double W, double M2_F, double M2_GT, double W0=beta_W0, double M=proton_M0);

/// correction to spectrum shape from recoil + weak magnetism according to Bilenkii 1959, eq. 11, after factoring out (1+3*lambda^2)
double Bilenkii59_RWM(double W);

/// Wilkinson ^VC(Z,W) as in [2] nucleon/lepton wavefunction convolution correction, Vector part
double WilkinsonVC(double Z, double W, double W0=beta_W0, double R=neutron_R0);
/// Wilkinson ^AC(Z,W) as in [2], nucleon/lepton wavefunction convolution correction, axial part
double WilkinsonAC(double Z, double W, double W0=beta_W0, double R=neutron_R0);
/// Combined Vector/Axial-Vector C
double CombinedC(double Z, double W, double M2_F, double M2_GT, double W0=beta_W0, double R=neutron_R0);

/// Wilkinson Q(Z,W,M) as in [0], nucleon recoil effect on Coulomb corrections
double WilkinsonQ(double Z,double W,double W0=beta_W0,double M=proton_M0);

/// Sirlin 1967 g * alpha/2pi radiative corrections to order alpha, also in [5]
double Sirlin_g_a2pi(double KE,double KE0,double m=m_e);
/// Wilkinson g * alpha/2pi: Sirlin g + fix for logarithm divergence [5]
double Wilkinson_g_a2pi(double W,double W0=beta_W0, double M=proton_M0);

/// shape factor for first forbidden Tensor/Axial decays, per J. Davidson, Phys. Rev. 82(1) p. 48, 1951
double Davidson_C1T(double W, double W0, double Z, double R);
/// shape factor for Cs137 second forbidden beta decay, per L.M. Langer and R.J.D. Moffat,  Phys. Rev. 82(5), p. 635, 1951
double Langer_Cs137_C2T(double W, double W0);
/// shape factor for Cs137 second forbidden beta decay, per H. Behrens and P. Christmas, Nucl. Phys. A399, pp. 131-140, 1983
double Behrens_Cs137_C(double W, double W0);

/// combined spectrum correction factor for unpolarized neutron beta decay
double neutronSpectrumCorrectionFactor(double KE);
/// corrected beta spectrum for unpolarized neutron beta decay
double neutronCorrectedBetaSpectrum(double KE);

/// beta decay spectrum calculating class
class BetaSpectrumGenerator {
public:
	/// constructor
	BetaSpectrumGenerator(double a, double z, double ep);
	
	/// shape correction to basic phase space
	double spectrumCorrectionFactor(double W) const;
	
	/// decay probability at given KE
	double decayProb(double KE) const;
	
	double A;				//< number of nucleons
	double Z;				//< number of protons
	double EP;				//< endpoint kinetic energy, keV
	double W0;				//< endpoint total energy, m_e*c^2
	double R;				//< effective nuclear radius
	double M0;				//< nuclear mass, m_e*c^2
	unsigned int forbidden;	//< "forbidden" level of decay
	double M2_F;			//< |M_F|^2 Fermi decay matrix element
	double M2_GT;			//< |M_GT|^2 Gamov-Teller decay matrix element
};


//-------------- A corrections ------------------

/// uncorrected asymmetry as a function of kinetic energy [keV]
inline double plainAsymmetry(double KE, double costheta=0.5) { return A0_PDG*beta(KE)*costheta; }

/// Shann's h * alpha/2pi radiative correction
double shann_h_a2pi(double KE, double KE0=neutronBetaEp, double m = m_e);
/// (h-g) * alpha/2pi radiative correction to A
double shann_h_minus_g_a2pi(double W, double W0=beta_W0);
/// Wilkinson weak magnetism + g_V*g_A interference + recoil correction to A [0]
double WilkinsonACorrection(double W);

/// combined order-alpha asymmetry corrections
inline double asymmetryCorrectionFactor(double KE) { double W = (KE+m_e)/m_e; return 1.0+WilkinsonACorrection(W)+shann_h_minus_g_a2pi(W); }
/// corrected asymmetry
inline double correctedAsymmetry(double KE, double costheta=0.5) { return plainAsymmetry(KE,costheta)*asymmetryCorrectionFactor(KE); }

#endif
