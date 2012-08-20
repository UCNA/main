#ifndef RESOURCE_HH
#define RESOURCE_HH 1

#include "Source.hh"
#include "Types.hh"
#include "EnergyCalibrator.hh"
#include "OutputManager.hh"
#include "ProcessedDataScanner.hh"
#include <vector>
#include <TH1F.h>
#include <TProfile.h>

/// class for fitting source data
class ReSourcer {
public:
	/// constructor
	ReSourcer(OutputManager* O, const Source& s, PMTCalibrator* P = NULL);
	
	OutputManager* OM;		//< output manager to write results to
	Source mySource;		//< source being fit for
	PMTCalibrator* PCal;	//< PMTCalibrator for estimating energy scale
	bool dbgplots;			//< whether to print extra "debugging" plots
	bool simMode;			//< whether the data is from simulation
	
	TH1F* hTubes[nBetaTubes+1][TYPE_III_EVENT+1];		//< calibrated PMT energy spectra by event type
	TH1F* hTubesRaw[nBetaTubes+1];						//< raw ADC energy spectra
	std::vector<SpectrumPeak> tubePeaks[nBetaTubes+1];	//< found peaks for each PMT
	TH1F* hitPos[2];									//< x,y hit profiles
				 
	unsigned int nBins;		//< histogram binning
	float eMin;				//< histogram lower energy
	float eMax;				//< histogram upper energy
	float pkMin;			//< minimum peak value (avoid Bi auger peak false alarms)
	float nSigma;			//< number of sigma to fit peaks
	
	TProfile* pPMTCorr[nBetaTubes][nBetaTubes];	//< correlation histograms between PMTs
	/// add correlation fit window
	void addCorrFit(float E0, float E1) { corrFitE0.push_back(E0); corrFitE1.push_back(E1); }
	/// calculate inter-PMT correlations
	void calcPMTcorr();
	
	/// fill histograms from source data; return whether point filled
	unsigned int fill(const ProcessedDataScanner& P);
	
	/// fit for source peaks, attach output to given subsystem, return tube spectrum histograms
	void findSourcePeaks(float runtime = 1.0);
	
	/// get number of counts
	unsigned int counts() const { assert(hTubes[nBetaTubes][TYPE_0_EVENT]); return (unsigned int)(hTubes[nBetaTubes][TYPE_0_EVENT]->GetEntries()); }

protected:
	std::vector<float> corrFitE0;	//< visible energy correlation measurement window start
	std::vector<float> corrFitE1;	//< visible energy correlation measurement window end

};

/// re-generate source fits / plots
void reSource(RunNum rn);

/// upload sources from run log
void uploadRunSources();


#endif
