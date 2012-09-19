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

const size_t kNumModules = 5;	//< number of DAQ modules for internal event header checks
const size_t kNumUCNMons = 4;	//< number of UCN monitors

enum UCN_MON_ID {
	UCN_MON_GV = 0,
	UCN_MON_SW = 1,
	UCN_MON_FE = 2,
	UCN_MON_SCS = 3
};

/// simple class for cuts from Stringmap
class RangeCut {
public:
	/// constructor
	RangeCut(const Stringmap& m = Stringmap());
	/// constructor with start, end times
	RangeCut(double s, double e): start(s), end(e) {}
	
	/// check if value is in range
	inline bool inRange(double x) const { return start <= x && x <= end; }
	
	double start;	//< cut minimum
	double end;		//< cut maximum
};

/// simple class for value + cuts range
class CutVariable {
public:
	/// constructor
	CutVariable(std::string sn=""): sname(sn) {}
	/// check if in range
	inline bool inRange() const { return R.inRange(val); }
	std::string sname;	//< sensor name (for pedestal subtraction)
	Float_t val;		//< stored value
	RangeCut R;			//< cuts range
};

/// cut blip in data
struct Blip {
	Blip(): start(0), end(0) {}
	BlindTime start;
	BlindTime end;
	BlindTime length() const { return end-start; }
};

/// new clean-ish re-write of data analyzer; try to be backward compatible with old output tree
class ucnaDataAnalyzer11b: public TChainScanner, public OutputManager {
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
	/// overall wirechamber cut
	inline bool passedMWPC(Side s) const { return fPassedCathMaxSum[s]; }
	/// scintillator TDC cut
	inline bool passedCutTDC(Side s) const { return fScint_tdc[s][nBetaTubes].inRange(); }
	/// overall wirechamber+scintillator cut
	inline bool Is2fold(Side s) const { return passedMWPC(s) && passedCutTDC(s); }
	/// overall muon tag check
	inline bool taggedMuon() const { return fTaggedBack[EAST] || fTaggedBack[WEST] || fTaggedDrift[EAST] || fTaggedDrift[WEST] || fTaggedTop[EAST]; }
	/// sis-tagged LED events
	inline bool isLED() const {return iSis00 & (1<<7); }
	/// sis-tagged UCN Mon events
	inline bool isUCNMon() const { return iSis00 & ((1<<2) | (1<<8) | (1<<9) | (1<<10) | (1<<11)); }
	/// specific UCN monitors
	inline bool isUCNMon(unsigned int n) const { return (iSis00 & (1<<2)) && (iSis00 & (1<<(8+n))); }
	/// sis-tagged scintillator trigger events
	inline bool isScintTrigger() const { return iSis00 & 3; }
	/// figure out whether this is a Bi pulser trigger
	bool isPulserTrigger();
	/// whether one PMT fired
	inline bool pmtFired(Side s, unsigned int t) const { return fScint_tdc[s][t].inRange(); }
	/// whether 2-of-4 trigger fired
	inline bool trig2of4(Side s) const { return fScint_tdc[s][nBetaTubes].val > 5; }
	/// calculate number of PMTs firing on given side
	unsigned int nFiring(Side s) const;
	/// qadc sum
	inline double qadcSum(Side s) const { double q = 0; for(unsigned int t=0; t<nBetaTubes; t++) q+=sevt[s].adc[t]; return q; }
	
	/// set whether to separately analyze LED events
	bool analyzeLED;
	/// set whether pedestals need to be measured
	bool needsPeds;
	
protected:
	// read variables
	Float_t r_Sis00;							//< Sis00 trigger flags
	Float_t r_TriggerNumber;					//< event trigger number
	BlindTime r_Clk;							//< event time blinded clock
	Float_t r_BClk;								//< proton beam clock
	Float_t r_Delt0;							//< time since last event
	Float_t r_AbsTime;							//< absolute time clock
	Float_t r_MonADC[kNumUCNMons];				//< UCN monitor ADCs
	Float_t r_PMTADC[2][nBetaTubes];			//< PMT ADCs
	Float_t r_PMTTDC[2][nBetaTubes+1];			//< PMT TDCs
	Float_t r_MWPC_caths[2][2][kMaxCathodes];		//< cathodes on [side][xplane][wire]
	Float_t r_MWPC_anode[2];					//< MWPC anode on each side
	Float_t r_Backing_TDC[2];					//< Backing Veto TDC
	Float_t r_Drift_TAC[2];						//< Drift tubes TAC
	Float_t r_Backing_ADC[2];					//< Backing Veto ADC
	Float_t r_Top_TDC[2];						//< Top Veto TDC (for East only)
	Float_t r_Top_ADC[2];						//< Top Veto ADC (for East only)
	Float_t r_Evnb[kNumModules];				//< header and footer counters per module
	Float_t r_Bkhf[kNumModules];				//< header and footer counters per module
	
	// whole run variables
	RunNum rn;									//< run number for file being processed
	PMTCalibrator PCal;							//< PMT Calibrator for this run
	CalDBSQL* CDBout;							//< output database connection
	std::vector<std::string> cathNames[2][2];	//< cathode sensor names on each [side][xplane]
	Float_t fAbsTimeStart;						//< absolute start time of run
	Float_t fAbsTimeEnd;						//< absolute end time of run
	BlindTime deltaT;							//< time scaler wraparound fix
	BlindTime totalTime;						//< total running time (accumulated from last event); after scanning, total ``live'' time (less global cuts)
	Float_t wallTime;							//< initial estimate of run time before actually scanning events; after scanning, total run time
	unsigned int nLiveTrigs;					//< number of triggers not removed by cuts
	bool ignore_beam_out;						//< whether to ignore long beam outages (e.g. for a source run)
	Float_t nFailedEvnb;						//< total Evnb failures
	Float_t nFailedBkhf;						//< total Bkhf failures
	RangeCut ScintSelftrig[2];					//< self-trigger range cut for each scintillator (used for Type I origin side determination)
	std::vector< std::pair<double,double> > manualCuts;	//< manually cut time segments
	std::vector<Blip> cutBlips;							//< keep track of cut run time
	
	
	// event-by-event calibrated variables
	int iTriggerNumber;						//< trigger number
	int iSis00;								//< Sis00 flags
	BlindTime fTimeScaler;					//< absolute event time scaler, blinded E, W, and unblinded
	CutVariable fBeamclock;					//< time since last beam pulse scaler
	Float_t fDelt0;							//< time since previous event
	CutVariable fWindow;					//< time window between previous and next event
	CutVariable fScint_tdc[2][nBetaTubes+1];//< TDC readout for each PMT and side
	ScintEvent sevt[2];						//< scintillator event, for reconstructing energy
	CutVariable fMWPC_anode[2];				//< anode ADC
	std::vector<float> fMWPC_caths[2][2];	//< cathodes on [side][xplane][wire]
	CutVariable fBacking_tdc[2];			//< muon backing veto TDC
	Float_t fBacking_adc[2];				//< muon backing veto ADC
	CutVariable fDrift_tac[2];				//< muon veto drift tubes TAC
	CutVariable fTop_tdc[2];				//< top veto TDCs (only East)
	Float_t fTop_adc[2];					//< top veto ADCs (only East)
	CutVariable fMonADC[kNumUCNMons];		//< UCN monitor ADCs = GV, Sw, Fe, SCS
	RollingWindow gvMonChecker;				//< rolling window check on gv monitor rate
	bool prevPassedCuts;					//< whether passed cuts on previous event
	bool prevPassedGVRate;					//< whether passed GV rate on previous event
	
	/// load a cut range for a CutVariable
	void loadCut(CutVariable& c, const std::string& cutName);
	/// load all cuts for run
	void loadCuts();
	
	/// set read points for input TChain
	virtual void setReadpoints();
	/// setup output tree, read points
	void setupOutputTree();
	
	// additional event variables for output tree
	TTree* TPhys;			//< physics events output tree
	TTree* TLED;			//< LED pulser events output tree
	Int_t fEvnbGood;		//< DAQ data quality checks
	Int_t fBkhfGood;		//< DAQ data quality checks
	Int_t fPassedAnode[2];	//< whether passed anode cut on each side 
	Int_t fPassedCath[2];	//< whether passed cathode sum cut on each side
	Int_t fPassedCathMax[2];//< whether passed cathode max cut on each side
	Int_t fPassedCathMaxSum[2];	//< whether passed cathode max sum on each side
	Int_t fPassedGlobal;	//< whether this event passed global beam/misc cuts on each side (counts for total run time)
	wireHit wirePos[2][2];	//< wire positioning data for [side][direction]
	CutVariable fCathSum[2];//< combined x+y cathode sum on each side
	CutVariable fCathMax[2];//< min(max cathode each plane) for each side
	CutVariable fCathMaxSum[2];	//< sum of max cathode from each plane for each side
	Float_t fEMWPC[2];		//< reconstructed energy deposition in wirechamber
	UInt_t fTaggedBack[2];	//< whether event was tagged by the muon backing veto on each side
	UInt_t fTaggedDrift[2];	//< whether event was tagged by muon veto drift tubes on each side
	UInt_t fTaggedTop[2];	//< whether event was tagged by top veto on each side (only meaningful on East)
	Side fSide;				//< event primary scintillator side
	EventType fType;		//< event backscatter type
	PID fPID;				//< event particle ID
	Float_t fProbIII;		//< event estimated probability of being a Type III backscatter
	Float_t fEtrue;			//< event reconstructed true energy
	
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
	/// classify event type
	void classifyEventType();
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
	// histograms
	TH1F* hCathMax[2][2];				//< cathode max, [side][cut]
	TH1F* hCathMaxSum[2][2];			//< cathode max sum, [side][cut]
	TH1F* hAnode[2][2];					//< anode, [side][cut]
	TH1F* hCathSum[2][2];				//< cathode sum, [side][cut]
	TH1F* hBackTDC[2];					//< backing TDC
	TH1F* hBackADC[2][2];				//< backing ADC, [side][cut]
	TH1F* hDriftTAC[2];					//< drift TAC [side]
	TH1F* hTopTDC[2];					//< top TDC [side] (East only)
	TH1F* hTopADC[2][2];				//< top ADC [side][cut] (East only)
	TH1F* hScintTDC[2][nBetaTubes+1];	//< scintillator 2-of-4 TDC by [side][tube]
	TH1F* hEtrue[2][5];					//< true energy, by [side][type]
	TH1F* hTuben[2][nBetaTubes];		//< individual PMT visible energy
	TH1F* hMonADC[kNumUCNMons];			//< UCN Monitor ADCs
	TH1F* hMonRate[kNumUCNMons];		//< UCM Monitor rates
	TH1F* hTypeRate[TYPE_IV_EVENT+1];	//< rate of event type for betas
	TH1F* hSideRate[2][2];				//< rate for [side][muon/beta]
	TH1F* hBkhfFailRate;				//< rate of bad Bkhf events
	TH1F* hEvnbFailRate;				//< rate of bad Evnb events
	TH1F* hHitsProfile[2][2];			//< 1D hit position histograms
	TH2F* hHitPos[2];					//< hit position on each side, 2D
	TH1F* hTrigEffic[2][nBetaTubes][2];			//< trigger efficiency for [side][tube][all/trig]
	std::vector<TH1*> hBiPulser[2][nBetaTubes];	//< Bi puser for [side][tube]
	TH1F* hClusterTiming[2];			//< event cluster timing for [all/beta]
};


#endif
