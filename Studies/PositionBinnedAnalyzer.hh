#ifndef POSITIONBINNEDANALYZER_HH
#define POSITIONBINNEDANALYZER_HH 1

#include "RunAccumulator.hh"
#include "SectorCutter.hh"

/// analyzer plugin for position-segmented analysis
class PositionBinnedAnalyzer: public AnalyzerPlugin {
public:
	/// constructor
	PositionBinnedAnalyzer(RunAccumulator* RA, const std::string& nm, unsigned int nr);
	/// allocate histograms for each position segment
	std::vector<fgbgPair*> allocateSegmentHistograms(TH1& hTemplate, AFPState a = AFP_OTHER, Side s = BOTH);
	
	SectorCutter sects;			//< sector cutter for position binning
	static double fidRadius;	//< analysis fiducial radius
};

#endif
