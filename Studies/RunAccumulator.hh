#ifndef RUNACCUMULATOR_HH
#define RUNACCUMULATOR_HH 1

#include "SegmentSaver.hh"
#include "Enums.hh"
#include "Types.hh"
#include "TagCounter.hh"
#include "ProcessedDataScanner.hh"
#include "G4toPMT.hh"
#include <TRandom3.h>

class AnalyzerPlugin;

/// background-subtracting histograms pair
class fgbgPair: private NoCopy {
public:
	/// constructor
	fgbgPair(const std::string& nm="", const std::string& ttl="", AFPState a = AFP_OTHER, Side s = BOTH);

	/// get (side,afp-mangled) name
	std::string getName() const { return baseName+(mySide<=WEST?(mySide==EAST?"_E":"_W"):"")+(afp<=AFP_ON?(afp?"_On":"_Off"):""); }
	/// get (side,afp-mangled) title
	std::string getTitle() const { return (mySide<=WEST?sideSubst(baseTitle,mySide):baseTitle)+(afp<=AFP_ON?(afp?" AFP On":" AFP Off"):""); }
	/// get histogram name
	std::string getHistoName(bool fg) const { return getName()+(fg==GV_CLOSED?"_BG":""); }
	/// get histogram title
	std::string getHistoTitle(bool fg) const { return getTitle()+(fg==GV_CLOSED?" Background":""); }
	/// set axis title
	void setAxisTitle(AxisDirection d, const std::string& ttl);
	/// perform background subtraction
	void bgSubtract(BlindTime tFG, BlindTime tBG);
	/// add another fgbgPair
	void operator+=(const fgbgPair& p);
	/// scale by a constant
	void operator*=(double c);
	
	std::string baseName;	//< base naming convention
	std::string baseTitle;	//< title naming convention
	TH1* h[2];				//< background, foreground pair
	AFPState afp;			//< AFP state for data (determines which time to use for BG subtraction)
	Side mySide;			//< side for data
	bool doSubtraction;		//< whether to do background subtraction
	bool doTimeScale;		//< whether to scale the BG by relative time for subtraction
	bool isSubtracted;		//< whether this pair is already background-subtracted
};

class RunAccumulator: public SegmentSaver, private NoCopy {
public:
	/// constructor
	RunAccumulator(OutputManager* pnt, const std::string& nm = "RunAccumulator", const std::string& inflName = "");
	/// destructor
	virtual ~RunAccumulator();
	
	/// add histograms from another RunAccumulator of the same type
	virtual void addSegment(const SegmentSaver& S);
	/// zero out run times
	void zeroCounters();
	/// scale all saved histograms by a factor
	virtual void scaleData(double s);
	
	/// create or load a FG/BG TH1F* set
	fgbgPair* registerFGBGPair(const std::string& hname, const std::string& title,
							  unsigned int nbins, float xmin, float xmax, AFPState a = AFP_OTHER, Side s = BOTH);
	/// create or load a FG/BG,OFF/ON histogram set based on a template TH1
	fgbgPair* registerFGBGPair(const TH1& hTemplate, AFPState a = AFP_OTHER, Side s = BOTH);
	/// check whether this has the named FGBGPair
	bool hasFGBGPair(const std::string& qname) const;
	/// get fgbg histogram by name
	fgbgPair& getFGBGPair(const std::string& qname);
	/// get fgbg histogram by name, const version
	const fgbgPair& getFGBGPair(const std::string& qname) const;
	/// make a new (unregistered) fgbgPair copy
	fgbgPair* cloneFGBGPair(const fgbgPair& p, const std::string& newName, const std::string& newTitle);

	/// get total run time for given state
	BlindTime getTotalTime(AFPState afp, bool fg) const;
	/// get total counts for given state
	float getTotalCounts(AFPState afp, GVState gv) const { return totalCounts[afp][gv]; }
	/// get time contributed from given rumber
	float getRunTime(RunNum rn) const { return runTimes[rn]; }
	/// save rates summary for core histograms to qOut
	void makeRatesSummary();
	/// write to QFile
	virtual void write(std::string outName = "");
	
	/// fill data from a ProcessedDataScanner
	virtual void loadProcessedData(AFPState afp, GVState gv, ProcessedDataScanner& PDS);
	/// fill data from simulations
	virtual void loadSimData(Sim2PMT& simData, unsigned int nToSim = 0, bool countAll = false);
	/// load single current event from simulator
	virtual void loadSimPoint(Sim2PMT& simData);
	/// load sim data to match a given run
	void simForRun(Sim2PMT& simData, RunNum rn, unsigned int nToSim, bool countAll = false);
	/// simulate for many runs, possibly apportioning a fixed number of counts among them; return number of events thrown
	unsigned int simMultiRuns(Sim2PMT& simData, const TagCounter<RunNum>& runReqs, unsigned int nCounts = 0);
	/// simulate background fluctuations based on "reference" data
	void simBgFlucts(const RunAccumulator& RefOA, double simfactor, bool addFluctCounts = true);
	/// perform background subtraction
	virtual void bgSubtractAll();
	
	AFPState currentAFP;			//< current state of AFP during data scanning
	GVState currentGV;				//< current foreground/background status during data scanning
	bool needsSubtraction;			//< whether background subtraction is pending
	bool isSimulated;				//< flag for whether this is based on simulated data
	
	TagCounter<RunNum> runCounts;	//< type-0 event counts by run, for re-simulation
	
	/// set current AFP, GV state
	void setCurrentState(AFPState afp, GVState gv);
	/// fill core histograms in plugins from data point
	virtual void fillCoreHists(ProcessedDataScanner& PDS, double weight);
	/// calculate results from filled histograms
	virtual void calculateResults();
	/// upload results to analysis DB
	virtual void uploadAnaResults();
	/// make plots from each plugin
	virtual void makePlots();
	/// MC/Data comparison plots/calculations from each plugin
	virtual void compareMCtoData(RunAccumulator& OAdata);
	/// add an analyzer plugin
	void addPlugin(AnalyzerPlugin* AP);
	/// get plugin by name
	AnalyzerPlugin* getPlugin(const std::string& nm);
	
	/// location of errorbar estimates for low-rate histograms
	virtual std::string estimatorHistoLocation() const { return ""; }
	
	std::map<std::string,fgbgPair*> fgbgHists;	//< background-subtractable quantities
	float totalCounts[AFP_OTHER+1][2];			//< total type-0 event counts by [flipper][fg/bg], for re-simulation
	BlindTime totalTime[AFP_OTHER+1][2]; 		//< total time for [flipper][fg/bg]
	TagCounter<RunNum> runTimes;				//< time spent on each run
	
protected:
	
	std::map<std::string,AnalyzerPlugin*> myPlugins;	//< analysis plugins
	static TRandom3 rnd_source;							//< random number source
	
	/// get matching RunAccumulator with "master" histograms for estimating error bars on low-counts bins
	RunAccumulator* getErrorEstimator();
};

/// generic analyzer plug-in class
class AnalyzerPlugin: private NoCopy {
public:
	/// constructor
	AnalyzerPlugin(RunAccumulator* RA, const std::string& nm): name(nm), myA(RA) { }
	/// destructor
	virtual ~AnalyzerPlugin() {}
	/// create or load a FG/BG TH1F* set
	fgbgPair* registerFGBGPair(const std::string& hname, const std::string& title,
							   unsigned int nbins, float xmin, float xmax,
							   AFPState a = AFP_OTHER, Side s = BOTH) { return myA->registerFGBGPair(hname,title,nbins,xmin,xmax,a,s); }
	/// create or load a FG/BG,OFF/ON histogram set based on a template TH1
	fgbgPair* registerFGBGPair(const TH1& hTemplate, AFPState a = AFP_OTHER, Side s = BOTH) { return myA->registerFGBGPair(hTemplate,a,s); }
	/// save canvas image
	void printCanvas(std::string fname, std::string suffix=".pdf") const { myA->printCanvas(fname,suffix); }
	
	std::string name;				//< plugin name
	RunAccumulator* myA;			//< RunAccumulator with which this plugin is associated
	AFPState currentAFP;			//< current state of AFP during data scanning
	GVState currentGV;				//< current foreground/background status during data scanning
	
	/// virtual routine for filling core histograms from data point
	virtual void fillCoreHists(ProcessedDataScanner& PDS, double weight) {}
	
	/// generate output plots
	virtual void makePlots() {}
	/// generate calculated hists
	virtual void calculateResults() {}
	/// upload results to analysis DB
	virtual void uploadAnaResults() {}
	/// virtual routine for MC/Data comparison plots/calculations
	/// NOTE: this MUST NOT change the contents of saved histograms (calculated ones are OK)
	virtual void compareMCtoData(AnalyzerPlugin* AP) {}
};

#endif
