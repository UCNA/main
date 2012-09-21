#ifndef OCTETANALYZER_HH
#define OCTETANALYZER_HH 1

#include "Enums.hh"
#include "SegmentSaver.hh"
#include "RunAccumulator.hh"
#include "ProcessedDataScanner.hh"
#include "G4toPMT.hh"
#include "Octet.hh"
#include <map>

/// histograms for [flipper][fg/bg]
class quadHists: private NoCopy {
public:
	/// constructor
	quadHists(const std::string& nm="", const std::string& ttl="", Side s = BOTH): name(nm), title(ttl), mySide(s), fillPoint(NULL) {}
	
	std::string name;	//< histogram base name
	std::string title;	//< histogram base title
	Side mySide;		//< side for these histograms (for blind time rates)
	TH1* fillPoint;		//< pointer to currently active fill histogram
	fgbgPair* fgbg[2];	//< background-subtracted pair for each flipper state
	
	/// add another quadHists
	void operator+=(const quadHists& p);
	/// scale by a constant
	void operator*=(double c);
	/// set fill point to point to appropriate histograms
	void setFillPoint(AFPState afp, GVState gv);
	/// check if another quadHists is equivalent (same size histograms)
	bool isEquivalent(const quadHists& qh) const;
	/// set background subtraction option
	void setSubtraction(bool b);
	/// set time scaling option
	void setTimeScaling(bool b);
	/// set draw range y minimum/maximum for all histograms
	void setDrawRange(double y, bool maximum);
	/// set draw range for specified axis
	void setRangeUser(double mn, double mx, AxisDirection d = X_DIRECTION);
	/// set axis title
	void setAxisTitle(AxisDirection d, const std::string& ttl);
	/// naming convention for each histogram
	std::string getHistoName(AFPState afp, bool fg) const { return fgbg[afp]->getHistoName(fg); }
	/// get (side-mangled) name
	std::string getName() const { return name+((mySide==EAST||mySide==WEST)?(mySide==EAST?"_E":"_W"):""); }
};

/// class for managing histogram inputs for asymmetry
class OctetAnalyzer: public RunAccumulator {
public:
	/// constructor, optionally with input filename
	OctetAnalyzer(OutputManager* pnt, const std::string& nm = "OctetAnalyzer", const std::string& inflName = "");
	/// destructor
	virtual ~OctetAnalyzer();
	
	/// create or load a FG/BG,OFF/ON histogram set
	quadHists* registerCoreHist(const std::string& hname, const std::string& title,
							   unsigned int nbins, float xmin, float xmax, Side s = BOTH);
	/// create or load a FG/BG,OFF/ON histogram set based on a template TH1
	quadHists* registerCoreHist(const TH1& hTemplate, Side s = BOTH);
	/// clone (and register) a quadHists
	quadHists* cloneQuadHist(const quadHists* qh, const std::string& newName, const std::string& newTitle);
	/// get core histogram by name
	quadHists* getCoreHist(const std::string& qname);
	/// get core histogram by name, const version
	const quadHists* getCoreHist(const std::string& qname) const;
	/// set all quadHists fill points
	void setFillPoints(AFPState afp, GVState gv);
	
	/// load data from ProcessedDataScanner
	void loadProcessedData(AFPState afp, ProcessedDataScanner& FG, ProcessedDataScanner& BG);
	/// load simulation data
	virtual void loadSimData(Sim2PMT& simData, unsigned int nToSim = 0, bool countAll = false);
	/// make a simulation clone (using simulation data from simData) of analyzed data in directory basedata; return number of cloned pulse-pairs
	unsigned int simuClone(const std::string& basedata, Sim2PMT& simData, double simfactor = 1., double replaceIfOlder = 0., bool doPlots = true);
	
	int depth;				//< octet division depth
	bool simPerfectAsym;	//< whether to simulate "perfect" asymmetry by re-using simulation events
	
private:
	
	std::map<std::string,quadHists*> coreHists;	//< core histograms for merging, BG subtraction
};

/// plug-in for OctetAnalyzer
class OctetAnalyzerPlugin: public AnalyzerPlugin {
public:
	/// constructor
	OctetAnalyzerPlugin(OctetAnalyzer* OA, const std::string& nm): AnalyzerPlugin(OA,nm), myA(OA) { }
	
	/// create or load a FG/BG,OFF/ON histogram set
	quadHists* registerCoreHist(const std::string& hname, const std::string& title,
								unsigned int nbins, float xmin,
								float xmax, Side s = BOTH) { return myA->registerCoreHist(hname,title,nbins,xmin,xmax,s); }
	/// create or load a FG/BG,OFF/ON histogram set based on a template TH1
	quadHists* registerCoreHist(const TH1& hTemplate, Side s = BOTH) { return myA->registerCoreHist(hTemplate,s); }
	
	// ---- some utility routines for common analysis/output operations ---- //
	
	/// calculate (blinded) super-ratio from quadHists for each side (optionally asymmetry of background and/or instrumental)
	TH1* calculateSR(const std::string& hname, const quadHists* qEast, const quadHists* qWest, bool fg=true, bool instr=false);
	/// calculate super-sum from quadHists for each side (optionally super-sum of background)
	TH1* calculateSuperSum(const std::string& hname, const quadHists* qEast, const quadHists* qWest, bool fg=true);
	/// draw whole quadHists of histograms
	void drawQuad(quadHists* qh, const std::string& subfolder = ".", const char* opt = "");
	/// draw East/West pair of quadHists together, optionally also drawing AFP states together
	void drawQuadSides(quadHists* qhE, quadHists* qhW, bool combineAFP = false,  const std::string& subfolder = ".", const std::string& opt = "");
	
	OctetAnalyzer* myA;	//< OctetAnalyzer with which this plugin is associated
};

/// process one pulse-pair worth of data
unsigned int processPulsePair(OctetAnalyzer& OA, const Octet& PP);

/// process a set of octets; return number of processed pulse-pairs
unsigned int processOctets(OctetAnalyzer& OA, const std::vector<Octet>& O, double replaceIfOlder = 0, bool doPlots = true);
	
#endif
