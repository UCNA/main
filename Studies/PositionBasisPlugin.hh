#ifndef POSITIONBASISANALYZER_HH
#define POSITIONBASISANALYZER_HH

#include "RunAccumulator.hh"

/// analyzer plugin for position-segmented analysis
class PositionBasisPlugin: public AnalyzerPlugin {
public:
	/// constructor
	PositionBasisPlugin(RunAccumulator* RA, const std::string& nm, unsigned int nr, double r);
	/// allocate histograms for each position segment
	//std::vector<fgbgPair*> allocateSegmentHistograms(TH1& hTemplate, AFPState a = AFP_OTHER, Side s = BOTH);
	
	// SectorCutter sects;		///< sector cutter for position binning
};

#endif
