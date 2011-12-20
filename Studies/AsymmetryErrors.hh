#ifndef ASYMMETRYERRORS_HH
#define ASYMMETRYERRORS_HH 1

#include "OutputManager.hh"
#include "PMTGenerator.hh"
#include "SRAsym.hh"
#include "SimNonlinearity.hh"
#include "Types.hh"
#include <TPad.h>
#include <TH1F.h>
#include "KurieFitter.hh"

/// data related to asymmetry-producing spectra
class AsymData {
public:
	/// constructor
	AsymData(): AsymSR(NULL), AsymBoner(NULL) {}
	/// write entries to a QFile
	void write(QFile& qf, std::string pfx) const;
	
	SRAsym* AsymSR;		//< super-ratio asymmetry
	SRAsym* AsymBoner;	//< bonehead asymmetry
	TH1F* hSimulated[2][2][nBetaTubes+1];	//< simulated smeared spectra
	TH1F* hSpinavg[2][nBetaTubes+1];		//< spin-averaged simulated spectra
	float_err kurieEp[2][2][nBetaTubes+1];	//< spectra endpoints
	float_err kurieEpAvg[2][nBetaTubes+1];	//< spin-averaged spectra endpoints
};

/// class for exploring effects of linearity errors on measured A
class AsymErrorExplorer {
public:
	
	/// constructor, generates reference spectra
	AsymErrorExplorer(RunGeometry g = GEOMETRY_OTHER, RunNum rn = 15916);	
	/// simulate distorted histograms for each side / flipper state
	void simulateSpinstates(OutputManager* OMz, std::string pfx, SimNonlinearizer* SNL, unsigned int nToSim, AsymData& ad);
	
	OutputManager* OM;		//< base directory output manager
	CalDBSQL CDB;			//< calibration DB connection
	PMTGenerator PGen;		//< event simulator
	TH1F* hInput[2][2];		//< input spectrum histograms
	float n0;				//< normalization for spectrum histograms
	
	AsymData refSpectra;	//< asymmetry data for reference spectra
	
	unsigned int nBetas;	//< number of betas to simulate in each state
	float refMultiple;		//< count multiple factor for 'reference' spectrum
	float flipon_effic;		//< reduced efficiency in flipper on state
	unsigned int nBinsE;	//< number of energy bins
	float eMax;				//< max energy to plot
	float fitMin;			//< fit range min
	float fitMax;			//< fit range max
	
	/// "repair" nonlinearity to match Kurie endpoints
	void gainmatchSNL(SimNonlinearizer& SNL, const AsymData& obs);	
	/// process a trial nonlinearity, compare to 'correct' results
	AsymData processTrial(OutputManager& OMz, SimNonlinearizer& SNL, std::string pfx);
};

/// main routine for asymmetry energy errors simulations
void AsymmetryErrors();

#endif
