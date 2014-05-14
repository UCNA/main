#ifndef GRAVITYSPECTROMETER_HH
#define GRAVITYSPECTROMETER_HH

#include "OctetAnalyzer.hh"

/// event y positions, by event time
class GravitySpectrometerPlugin: public OctetAnalyzerPlugin {
public:
	/// constructor
	GravitySpectrometerPlugin(OctetAnalyzer* OA);
	
	/// fill from scan data point
	virtual void fillCoreHists(ProcessedDataScanner& PDS, double weight);
	/// calculate offset info
	virtual void calculateResults();
	/// output plot generation
	virtual void makePlots();
	
	quadHists* qTime;			///< event positions, by height and time
	quadHists* qHeight;			///< event positions height
	
	TH1* hAFPRat;				///< AFP off/on ratio vs position
	TH1* timeProf[AFP_ON+1];	///< time profiles by [afp]
	TH1* hTimeRat;				///< afp off/on counts vs time
	std::vector<TH1F*> timeSlices[AFP_ON+1];
	TGraphErrors* timeEvol[AFP_ON+1];
	
};

#endif
