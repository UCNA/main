#ifndef ENDPOINTSTUDY_HH
#define ENDPOINTSTUDY_HH 1

#include "SpectrumHistos.hh"
#include "GraphicsUtils.hh"
#include "SectorCutter.hh"
#include "EnergyCalibrator.hh"
#include "PMTGenerator.hh"
#include "DataSource.hh"
#include "KurieFitter.hh"
#include "RunAccumulator.hh"

/// class for position-segmented analysis
class PositionBinner: public RunAccumulator {
public:
	/// constructor
	PositionBinner(OutputManager* pnt, const std::string& nm, float r, unsigned int nr, const std::string& infl="");
	
	/// process a data point into position histograms
	virtual void fillCoreHists(ProcessedDataScanner& PDS, double weight);
	
	/// determine what section a point belongs to
	virtual unsigned int getSector(float x, float y) { return sects.sector(x,y); }
	/// get number of sectors
	virtual unsigned int getNSectors() const { return sects.nSectors(); }
	
	/// fit enpoints in each sector
	void fitEndpoints();
	
protected:
	std::vector<fgbgPair> sectEnergy[2][nBetaTubes];	//< visible light histograms by position, PMT
	fgbgPair energySpectrum[2];						//< combined energy spectrum, each side
	std::vector<fgbgPair> sectAnode[2];				//< anode histograms by position
	SectorCutter sects;								//< sector cutter
};

#endif
