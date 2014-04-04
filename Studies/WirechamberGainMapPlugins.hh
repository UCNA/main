#ifndef WIRECHAMBERGAINMAPPLUGINS_HH
#define WIRECHAMBERGAINMAPPLUGINS_HH

#include "PositionBinnedPlugin.hh"

/// base class for wirechamber position mapping functions
class WirechamberGainMapPluginBase: public PositionBinnedPlugin {
public:
	/// constructor
	WirechamberGainMapPluginBase(RunAccumulator* RA, unsigned int nr, const std::string& nm);
	
	/// process a data point into position histograms
	virtual void fillCoreHists(ProcessedDataScanner& PDS, double weight);
	
	/// identify histogram to book event in; return false if event discarded
	bool evtLoc(const ProcessedDataScanner& PDS, Side& s, unsigned int& m, float& x, float& y) const;
	
	/// generate position map, upload to CalDB
	void genPosmap(const std::string& pmapNameBase) const;
	
	std::vector<fgbgPair*> sectHists[BOTH];	//< reconstructed energy histograms for each sector
	fgbgPair* sectGains[BOTH];				//< prior average gains in each sector
	
protected:
	ChargeProxyType myChgPrx;				//< what charge measure is being used
};

/// simulated wirechamber energy deposition position map
class WirechamberEdepMapPlugin: public WirechamberGainMapPluginBase {
	public:
	WirechamberEdepMapPlugin(RunAccumulator* RA, unsigned int nr): WirechamberGainMapPluginBase(RA,nr,"MWPC_SimEdepPos") {}
	/// process a data point into position histograms
	virtual void fillCoreHists(ProcessedDataScanner& PDS, double weight);
};

/// HACK: stub edep by position plugin for data, to fake out SectorCutter info for simulation
class WirechamberNullEdepMapPlugin: public PositionBinnedPlugin {
public:
	/// constructor
	WirechamberNullEdepMapPlugin(RunAccumulator* RA, unsigned int nr): PositionBinnedPlugin(RA,"MWPC_SimEdepPos",nr) {}
	/// fill histograms
	virtual void fillCoreHists(ProcessedDataScanner&, double) {}
};


/// anode signal position map
class AnodeGainMapPlugin: public WirechamberGainMapPluginBase {
public:
	/// constructor
	AnodeGainMapPlugin(RunAccumulator* RA, unsigned int nr): WirechamberGainMapPluginBase(RA,nr,"MWPC_AnodePos") { myChgPrx = CHARGE_PROXY_ANODE; }
};

/// charge cloud size position map
class CCloudGainMapPlugin: public WirechamberGainMapPluginBase {
public:
	/// constructor
	CCloudGainMapPlugin(RunAccumulator* RA, unsigned int nr): WirechamberGainMapPluginBase(RA,nr,"MWPC_CCloudPos") { myChgPrx = CHARGE_PROXY_CCLOUD; }
};

#endif
