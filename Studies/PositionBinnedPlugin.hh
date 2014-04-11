#ifndef POSITIONBINNEDANALYZER_HH
#define POSITIONBINNEDANALYZER_HH

#include "RunAccumulator.hh"
#include "SectorCutter.hh"

/// analyzer plugin for position-segmented analysis
class PositionBinnedPlugin: public AnalyzerPlugin {
public:
	/// constructor
	PositionBinnedPlugin(RunAccumulator* RA, const std::string& nm, unsigned int nr);
	/// allocate histograms for each position segment
	std::vector<fgbgPair*> allocateSegmentHistograms(TH1& hTemplate, AFPState a = AFP_OTHER, Side s = BOTH);
	
	SectorCutter sects;			///< sector cutter for position binning
	static double fidRadius;	///< analysis fiducial radius
};

#endif
