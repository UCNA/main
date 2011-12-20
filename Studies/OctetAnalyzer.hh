#ifndef OCTETANALYZER_HH
#define OCTETANALYZER_HH 1

#include "Enums.hh"
#include "SegmentSaver.hh"
#include "RunAccumulator.hh"
#include "ProcessedDataScanner.hh"
#include "DataSource.hh"
#include "G4toPMT.hh"
#include "Octet.hh"
#include <map>

/// histograms for [flipper][fg/bg]
class quadHists {
public:
	/// constructor
	quadHists(std::string nm="", Side s = BOTH): name(nm), mySide(s), fillPoint(NULL) {}
	
	std::string name;	//< histogram base name
	Side mySide;		//< side for these histograms (for blind time rates)
	TH1** fillPoint;	//< location of currently active histogram
	fgbgPair fgbg[2];	//< background-subtracted pair for each flipper state
	
	/// add another quadHists
	void operator+=(const quadHists& p);
	/// scale by a constant
	void operator*=(double c);
	/// set fill point to point to appropriate histograms
	void setFillPoint(AFPState afp, GVState gv);
	/// check if another quadHists is equivalent (same size histograms)
	bool isEquivalent(const quadHists& qh) const;
	/// set draw range minimum for all histograms
	void setDrawMinimum(double y);
	/// naming convention for each histogram
	std::string getHistoName(AFPState afp, bool fg) const { return fgbg[afp].getHistoName(fg); }
	/// get (side-mangled) name
	std::string getName() const { return name+((mySide==EAST||mySide==WEST)?(mySide==EAST?"_E":"_W"):""); }
};

/// class for managing histogram inputs for asymmetry
class OctetAnalyzer: public RunAccumulator {
public:
	/// constructor, optionally with input filename
	OctetAnalyzer(OutputManager* pnt, const std::string& nm = "OctetAnalyzer", const std::string& inflName = "");
	/// destructor
	virtual ~OctetAnalyzer() {}
	
	/// create or load a FG/BG,OFF/ON histogram set
	quadHists registerCoreHist(const std::string& hname, const std::string& title,
							   unsigned int nbins, float xmin, float xmax,
							   Side s = BOTH, TH1F** fillPoint = NULL);
	/// create or load a FG/BG,OFF/ON histogram set based on a template TH1
	quadHists registerCoreHist(const TH1& hTemplate, Side s = BOTH, TH1** fillPoint = NULL);
	/// clone (and register) a quadHists
	quadHists cloneQuadHist(const quadHists& qh, const std::string& newName);
	/// get core histogram by name
	quadHists& getCoreHist(const std::string& qname);
	/// get core histogram by name, const version
	const quadHists& getCoreHist(const std::string& qname) const;
	
	/// load data from ProcessedDataScanner
	void loadProcessedData(AFPState afp, ProcessedDataScanner& FG, ProcessedDataScanner& BG);
	/// load simulation data
	virtual void loadSimData(Sim2PMT& simData, unsigned int nToSim);
		
	// ---- some utility routines for common analysis/output operations ---- //
	
	/// calculate (blinded) super-ratio from quadHists for each side (optionally asymmetry of background)
	TH1* calculateSR(const std::string& hname, const quadHists& qEast, const quadHists& qWest, bool fg=true);
	/// calculate super-sum from quadHists for each side (optionally super-sum of background)
	TH1* calculateSuperSum(const std::string& hname, const quadHists& qEast, const quadHists& qWest, bool fg=true);
	/// draw whole quadHists of histograms
	void drawQuad(quadHists& qh, const std::string& subfolder = ".", const char* opt = "");
	/// draw East/West pair of quadHists together, optionally also drawing AFP states together
	void drawQuadSides(quadHists& qhE, quadHists& qhW, bool combineAFP = false,  const std::string& subfolder = ".", const std::string& opt = "");
	
	int depth;		//< octet division depth
		
protected:
	
	/// set all quadHists fill points
	void setFillPoints(AFPState afp, GVState gv);
	
private:
	
	std::map<std::string,quadHists> coreHists;	//< core histograms for merging, BG subtraction
};

/// process one pulse-pair worth of data
unsigned int processPulsePair(OctetAnalyzer& OA, const Octet& PP);

/// process a set of octets; return number of processed pulse-pairs
unsigned int processOctets(OctetAnalyzer& OA, const std::vector<Octet>& O, double replaceIfOlder = 0);

/// make a simulation clone (using simulation data from simData) of analyzed data in directory basedata; return number of cloned pulse-pairs
unsigned int simuClone(const std::string& basedata, OctetAnalyzer& OA, Sim2PMT& simData, double simfactor = 1.0, double replaceIfOlder = 0);
	
#endif
