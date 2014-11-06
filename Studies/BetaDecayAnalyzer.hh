#ifndef BETADECAYANALYZER_HH
#define BETADECAYANALYZER_HH

#include "MuonPlugin.hh"
#include "PositionsPlugin.hh"
#include "CathodeTuningAnalyzer.hh"
#include "AsymmetryPlugin.hh"
#include "SimAsymmetryPlugin.hh"
#include "SimTreePlugin.hh"
#include "HighEnergyExcessPlugin.hh"

/// analyzer for beta decay data
class BetaDecayAnalyzer: public OctetAnalyzer {
public:
	/// constructor
	BetaDecayAnalyzer(OutputManager* pnt, const std::string& nm = "BetaDecayAnalyzer", const std::string& inflName = "");
	/// create a new instance of this object (cloning self settings) for given directory
	virtual SegmentSaver* makeAnalyzer(const std::string& nm, const std::string& inflname) { return new BetaDecayAnalyzer(this,nm,inflname); }
	
	MuonPlugin* myMuons;			///< muons plugin
	PositionsPlugin* myPos;			///< positions plugin
	HighEnergyExcessPlugin* myHEE;	///< high energy excess events, indicating neutron generated backgrounds
	
	AsymmetryPlugin* myAsym;		///< asymmetry plugin
};

/// analyzer for beta decay simulation
class SimBetaDecayAnalyzer: public OctetAnalyzer {
public:
	/// constructor
	SimBetaDecayAnalyzer(OutputManager* pnt, const std::string& nm = "BetaDecayAnalyzer", const std::string& inflName = "");
	/// create a new instance of this object (cloning self settings) for given directory
	virtual SegmentSaver* makeAnalyzer(const std::string& nm, const std::string& inflname);
	
	MuonPlugin* myMuons;			///< muons plugin
	PositionsPlugin* myPos;			///< positions plugin
	AsymmetryPlugin* myAsym;		///< asymmetry plugin

	SimAsymmetryPlugin* mySimAsym;	///< simulated asymmetry plugin
};

class SimTreeBetaDecayAnalyzer: public OctetAnalyzer {
public:
	/// constructor
        SimTreeBetaDecayAnalyzer(OutputManager* pnt, const std::string& nm = "BetaDecayAnalyzer", const std::string& inflName = "");
        /// create a new instance of this object (cloning self settings) for given directory
	virtual SegmentSaver* makeAnalyzer(const std::string& nm, const std::string& inflname);

	SimTreePlugin* mySimTree;	///< simulated tree plugin
};

#endif
