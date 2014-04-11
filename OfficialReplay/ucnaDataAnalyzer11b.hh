#ifndef UCNADATAANALYZER11B
#define UCNADATAANALYZER11B 1

#include "OutputManager.hh"
#include "TChainScanner.hh"
#include "Enums.hh"
#include "Types.hh"
#include "EnergyCalibrator.hh"
#include "CalDBSQL.hh"
#include "ManualInfo.hh"
#include "RollingWindow.hh"
#include "EventClassifier.hh"

const size_t kNumModules = 5;	///< number of DAQ modules for internal event header checks
const size_t kNumUCNMons = 4;	///< number of UCN monitors

enum UCN_MON_ID {
	UCN_MON_GV = 0,
	UCN_MON_SW = 1,
	UCN_MON_FE = 2,
	UCN_MON_SCS = 3
};

/// cut blip in data
struct Blip {
	Blip(): start(0), end(0) {}
	BlindTime start;
	BlindTime end;
	BlindTime length() const { return end-start; }
};

/// new clean-ish re-write of data analyzer; try to be backward compatible with old output tree
class ucnaDataAnalyzer11b: public TChainScanner, public OutputManager, public EventClassifier {
public:
	/// constructor
	ucnaDataAnalyzer11b(RunNum R, std::string bp, CalDB* CDB);
	/// destructor
	virtual ~ucnaDataAnalyzer11b() {}
	
	/// run analysis
	void analyze();
	/// set output DB connection
	inline void setOutputDB(CalDBSQL* CDB = NULL) { CDBout = CDB; }
	/// set to ignore beam cuts
	inline void setIgnoreBeamOut(bool ibo) { ignore_beam_out = ibo; }
	
	/// beam + data cuts
	bool passesBeamCuts();
	/// figure out whether this is a Bi pulser trigger
	bool isPulserTrigger();
	/// whether one PMT fired
	inline bool pmtFired(Side s, unsigned int t) const { return fScint_tdc[s][t].inRange(); }
	/// whether 2-of-4 trigger fired
	inline bool trig2of4(Side s) const { return fScint_tdc[s][nBetaTubes].val > 5; }
	/// calculate number of PMTs firing on given side
	unsigned int nFiring(Side s) const;
	
	/// set whether to separately analyze LED events
	bool analyzeLED;
	/// set whether pedestals need to be measured
	bool needsPeds;
	/// set whether to produce color or black-and-white plots
	bool colorPlots;
	
	//---- event classifier subclass functions
	/// event type classification
	virtual void classifyEvent();
	/// Type I initial hit side determination --- needed in subclass
	virtual Side getFirstScint() const { return fScint_tdc[WEST][nBetaTubes].val < ScintSelftrig[WEST].start ? EAST : WEST; }
	/// Type II/III separation probability --- needed in subclass
	virtual float getProbIII() const { return WirechamberCalibrator::sep23Prob(fSide, sevt[EAST].energy.x + sevt[WEST].energy.x, fEMWPC[fSide]); }
	
protected:
	// read variables
	Float_t r_Sis00;							///< Sis00 trigger flags
	Float_t r_TriggerNumber;					///< event trigger number
	BlindTime r_Clk;							///< event time blinded clock
	Float_t r_BClk;								///< proton beam clock
	Float_t r_Delt0;							///< time since last event
	Float_t r_AbsTime;							///< absolute time clock
	Float_t r_MonADC[kNumUCNMons];				///< UCN monitor ADCs
	Float_t r_PMTADC[BOTH][nBetaTubes];			///< PMT ADCs
	Float_t r_PMTTDC[BOTH][nBetaTubes+1];		///< PMT TDCs
	Float_t r_MWPC_caths[BOTH][2][kMaxCathodes];///< cathodes on [side][xplane][wire]
	Float_t r_MWPC_anode[BOTH];					///< MWPC anode on each side
	Float_t r_Backing_TDC[BOTH];				///< Backing Veto TDC
	Float_t r_Drift_TAC[BOTH];					///< Drift tubes TAC
	Float_t r_Backing_ADC[BOTH];				///< Backing Veto ADC
	Float_t r_Top_TDC[BOTH];					///< Top Veto TDC (for East only)
	Float_t r_Top_ADC[BOTH];					///< Top Veto ADC (for East only)
	Float_t r_Evnb[kNumModules];				///< header and footer counters per module
	Float_t r_Bkhf[kNumModules];				///< header and footer counters per module
	
	// whole run variables
	RunNum rn;									///< run number for file being processed
	PMTCalibrator PCal;							///< PMT Calibrator for this run
	CalDBSQL* CDBout;							///< output database connection
	std::vector<std::string> cathNames[BOTH][2];///< cathode sensor names on each [side][xplane]
	Float_t fAbsTimeStart;						///< absolute start time of run
	Float_t fAbsTimeEnd;						///< absolute end time of run
	BlindTime deltaT;							///< time scaler wraparound fix
	BlindTime totalTime;						///< total running time (accumulated from last event); after scanning, total ``live'' time (less global cuts)
	Float_t wallTime;							///< initial estimate of run time before actually scanning events; after scanning, total run time
	unsigned int nLiveTrigs;					///< number of triggers not removed by cuts
	bool ignore_beam_out;						///< whether to ignore long beam outages (e.g. for a source run)
	Float_t nFailedEvnb;						///< total Evnb failures
	Float_t nFailedBkhf;						///< total Bkhf failures
	RangeCut ScintSelftrig[BOTH];				///< self-trigger range cut for each scintillator (used for Type I origin side determination)
	std::vector< std::pair<double,double> > manualCuts;	///< manually cut time segments
	std::vector<Blip> cutBlips;							///< keep track of cut run time
	
	
	// event-by-event calibrated variables
	int iTriggerNumber;							///< trigger number
	BlindTime fTimeScaler;						///< absolute event time scaler, blinded E, W, and unblinded
	CutVariable fBeamclock;						///< time since last beam pulse scaler
	Float_t fDelt0;								///< time since previous event
	CutVariable fWindow;						///< time window between previous and next event
	CutVariable fScint_tdc[BOTH][nBetaTubes+1];	///< TDC readout for each PMT and side
	ScintEvent sevt[BOTH];						///< scintillator event, for reconstructing energy
	CutVariable fMWPC_anode[BOTH];				///< anode ADC
	float fMWPC_caths[BOTH][2][kMaxCathodes];	///< cathodes on [side][xplane][wire]
	CutVariable fBacking_tdc[BOTH];				///< muon backing veto TDC
	Float_t fBacking_adc[BOTH];					///< muon backing veto ADC
	CutVariable fDrift_tac[BOTH];				///< muon veto drift tubes TAC
	CutVariable fTop_tdc[BOTH];					///< top veto TDCs (only East)
	Float_t fTop_adc[BOTH];						///< top veto ADCs (only East)
	CutVariable fMonADC[kNumUCNMons];			///< UCN monitor ADCs = GV, Sw, Fe, SCS
	RollingWindow gvMonChecker;					///< rolling window check on gv monitor rate
	bool prevPassedCuts;						///< whether passed cuts on previous event
	bool prevPassedGVRate;						///< whether passed GV rate on previous event
	MWPCevent mwpcs[BOTH];						///< MWPC information
	
	/// load all cuts for run
	void loadCuts();
	
	/// set read points for input TChain
	virtual void setReadpoints();
	/// setup output tree, read points
	void setupOutputTree();
	
	// additional event variables for output tree
	TTree* TPhys;								///< physics events output tree
	TTree* TLED;								///< LED pulser events output tree
	Int_t fPassedGlobal;						///< whether this event passed global beam/misc cuts on each side (counts for total run time)
	wireHit wirePos[BOTH][2];					///< wire positioning data for [side][direction]
	CutVariable fCathSum[BOTH];					///< combined x+y cathode sum on each side
	CutVariable fCathMax[BOTH];					///< min(max cathode each plane) for each side
	CutVariable fCathMaxSum[BOTH];				///< sum of max cathode from each plane for each side
	Float_t fEMWPC[BOTH];						///< reconstructed energy deposition in wirechamber
	Float_t fEtrue;								///< event reconstructed true energy
	
	/// pre-scan data to extract pedestals
	void pedestalPrePass();
	/// fit pedestals
	void monitorPedestal(const std::vector<float>& vdata, const std::vector<float>& vtime,
						 const std::string& mon_name, double graphWidth, float tmin=60., unsigned int cmin=3000, bool isPed=true);
	
	/*--- event processing loop ---*/
	/// process current event raw->phys
	bool processEvent();
	/// move readin variables to appropriate locations for processing
	void convertReadin();
	/// check event headers for errors
	void checkHeaderQuality();
	/// fix scaler overflows, convert times to seconds
	void calibrateTimes();
	/// reconstruct wirechamber positions
	void reconstructPosition();
	/// apply PMT calibrations to get visible energy
	void reconstructVisibleEnergy();
	/// classify muon veto response
	void checkMuonVetos();
	/// reconstruct true energy based on event type
	void reconstructTrueEnergy();
	
	/*--- end of processing ---*/
	/// trigger efficiency curves
	void calcTrigEffic();
	/// Bi pulser gain stabilizer
	void processBiPulser();
	/// estimate muon veto accidentals rates
	void muonVetoAccidentals();
	/// tally total run time
	void tallyRunTime();
	/// output replay summary info, optionally into AnalysisDB
	void replaySummary();
	/// print info to replace old "quick analyzer"
	void quickAnalyzerSummary() const;
	
	
	/*--- histograms ---*/
	/// set up summary histograms
	void setupHistograms();
	/// fill summary histograms with few calibration dependencies
	void fillEarlyHistograms();
	/// fill summary histograms from event data
	void fillHistograms();
	/// locate sources on positions histogram
	void locateSourcePositions();
	/// output histograms of interest
	void plotHistos();
	/// draw cut range lines
	void drawCutRange(const RangeCut& r, Int_t c=4);
	/// draw regions excluded by blip cuts
	void drawExclusionBlips(Int_t c=4);
	/// set drawing style for second histogram in group
	void setSecondaryStyle(TH1* h);
	/// set drawing style based on detector side
	void setSideStyle(TH1* h, Side s);
	/// scale histogram to rate
	void scaleToRate(TH1* h, double rscale = 1.0);
	// histograms
	TH1F* hCathMax[BOTH][2];						///< cathode max, [side][cut]
	TH1F* hCathMaxSum[BOTH][2];						///< cathode max sum, [side][cut]
	TH1F* hAnode[BOTH][2];							///< anode, [side][cut]
	TH1F* hCathSum[BOTH][2];						///< cathode sum, [side][cut]
	TH1F* hBackTDC[BOTH];							///< backing TDC
	TH1F* hBackADC[BOTH][2];						///< backing ADC, [side][cut]
	TH1F* hDriftTAC[BOTH];							///< drift TAC [side]
	TH1F* hTopTDC[BOTH];							///< top TDC [side] (East only)
	TH1F* hTopADC[BOTH][2];							///< top ADC [side][cut] (East only)
	TH1F* hScintTDC[BOTH][nBetaTubes+1];			///< scintillator 2-of-4 TDC by [side][tube]
	TH1F* hEtrue[BOTH][TYPE_IV_EVENT+1];			///< true energy, by [side][type]
	TH1F* hTuben[BOTH][nBetaTubes];					///< individual PMT visible energy
	TH1F* hMonADC[kNumUCNMons];						///< UCN Monitor ADCs
	TH1F* hMonRate[kNumUCNMons];					///< UCM Monitor rates
	TH1F* hTypeRate[TYPE_IV_EVENT+1];				///< rate of event type for betas
	TH1F* hSideRate[BOTH][2];						///< rate for [side][muon/beta]
	TH1F* hBkhfFailRate;							///< rate of bad Bkhf events
	TH1F* hEvnbFailRate;							///< rate of bad Evnb events
	TH1F* hHitsProfile[BOTH][2];					///< 1D hit position histograms [side][plane]
	TH2F* hHitPos[BOTH];							///< hit position on each side, 2D
	TH1F* hTrigEffic[BOTH][nBetaTubes][2];			///< trigger efficiency for [side][tube][all/trig]
	std::vector<TH1*> hBiPulser[BOTH][nBetaTubes];	///< Bi puser for [side][tube]
	TH1F* hClusterTiming[2];						///< event cluster timing for [all/beta]
};


#endif
