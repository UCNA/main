/// \file BetaSpectrum.hh calculations for the beta spectrum and various related corrections
#ifndef BETASPECTRUM_HH
#define BETASPECTRUM_HH 1

// references:
// [0] Wilkinson, Analysis of Neutron Beta-Decay, Nucl. Phys. A 377 (1982) 474-504
// [1] Wilkinson, Evaluation of Beta-Decay I, NIM A 275 (1989) 378-386
// [2] Wilkinson, Evaluation of Beta-Decay II, NIM A 290 (1990) 509-515
// [3] Wilkinson, Evaluation of Beta-Decay III

#include <math.h>
// useful physics constants
const double neutronBetaEp = 782.347;	//< neutron beta decay endpoint, keV
const double m_e = 511.00;				//< electron mass, keV/c^2
const double m_p = 938272.046;			//< proton mass, keV/c^2
const double m_n = m_p-m_e-neutronBetaEp;	//< neutron mass, keV/c^2
const double alpha = 1./137.036;		//< fine structure constant
const double lambda = fabs(-1.2694);	//< +/-0.0028, PDG 2010 value, Wilkinson sign convention
const double A0_PDG = -0.1173;			//< +/-0.0013, PDG 2010 value
const double beta_W0 = (neutronBetaEp+m_e)/m_e;	//< beta spectrum endpoint, ``natural'' units
const double neutron_R0 = 0.0025896*1.2;		//< neutron and proton radius approximation in ``natural'' units
const double proton_M0 = m_p/m_e;				//< proton mass, ``natural'' units
const double neutron_M0 = m_n/m_e;				//< neutron mass, ``natural'' units

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
/// correction to spectrum shape from recoil, weak magnetism according to Bilenkii 1958, eq. 11
double Bilenkii_1958_11(double W);
/// Wilkinson ^VC(Z,W) as in [2] nucleon/lepton wavefunction convolution correction, Vector part
double WilkinsonVC(double Z, double W, double W0=beta_W0, double R=neutron_R0);
/// Wilkinson ^AC(Z,W) as in [2], nucleon/lepton wavefunction convolution correction, axial part
double WilkinsonAC(double Z, double W, double W0=beta_W0, double R=neutron_R0);

/// Wilkinson Q(Z,W,M) as in [0], nucleon recoil effect on Coulomb corrections
double WilkinsonQ(double Z,double W,double W0=beta_W0,double M=proton_M0);

/// Sirlin 1967 g * alpha/2pi radiative corrections to order alpha, also in [5]
double Sirlin_g(double KE,double KE0,double m=m_e);
/// Wilkinson g: Sirlin g + fix for logarithm divergence [5]
double Wilkinson_g(double W,double W0=beta_W0);

/// combined spectrum correction factor
double spectrumCorrectionFactor(double KE,int A = 1, int Z = 1, double ep = neutronBetaEp);
/// corrected beta spectrum
double correctedBetaSpectrum(double KE, int A = 1, int Z = 1, double ep = neutronBetaEp);

//-------------- A corrections ------------------

/// uncorrected asymmetry as a function of kinetic energy [keV]
inline double plainAsymmetry(double KE, double costheta=0.5) { return A0_PDG*beta(KE)*costheta; }

/// h*alpha/2pi radiative correction
double shann_h(double KE, double KE0=neutronBetaEp, double m = m_e);
/// (h-g)*alpha/2pi radiative correction to A
double shann_h_minus_g(double W, double W0=beta_W0);
/// Wilkinson weak magnetism + g_V*g_A interference + recoil correction to A [0]
double WilkinsonACorrection(double W);

/// combined order-alpha asymmetry corrections
inline double asymmetryCorrectionFactor(double KE) { double W = (KE+m_e)/m_e; return 1.0+WilkinsonACorrection(W)+shann_h_minus_g(W); }
/// corrected asymmetry
inline double correctedAsymmetry(double KE, double costheta=0.5) { return plainAsymmetry(KE,costheta)*asymmetryCorrectionFactor(KE); }

#endif
