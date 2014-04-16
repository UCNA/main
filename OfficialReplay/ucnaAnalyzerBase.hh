#ifndef UCNAANALYZERBASE_HH
#define UCNAANALYZERBASE_HH

#include "OutputManager.hh"
#include "TChainScanner.hh"
#include "Enums.hh"
#include "Types.hh"
#include "EnergyCalibrator.hh"
#include "CalDBSQL.hh"
#include "ManualInfo.hh"

const size_t kNumModules = 5;	///< number of DAQ modules for internal event header checks
const size_t kNumUCNMons = 4;	///< number of UCN monitors

enum UCN_MON_ID {
	UCN_MON_GV = 0,
	UCN_MON_SW = 1,
	UCN_MON_FE = 2,
	UCN_MON_SCS = 3
};

/// new clean-ish re-write of data analyzer; try to be backward compatible with old output tree
class ucnaAnalyzerBase: public TChainScanner, public OutputManager {
public:
	/// constructor
	ucnaAnalyzerBase(RunNum R, const std::string& bp, const std::string& nm, CalDB* CDB);
	/// destructor
	virtual ~ucnaAnalyzerBase() {}
	
	/// set to ignore beam cuts
	void setIgnoreBeamOut(bool ibo) { ignore_beam_out = ibo; }
	
	/// whether one PMT fired
	inline bool pmtFired(Side s, unsigned int t) const { return fScint_tdc[s][t].inRange(); }
	/// whether 2-of-4 trigger fired
	inline bool trig2of4(Side s) const { return fScint_tdc[s][nBetaTubes].val > 5; }
	/// calculate number of PMTs firing on given side
	unsigned int nFiring(Side s) const;
	
	/// whether event passes beam time cuts (+ manual time cuts)
	bool passesBeamCuts() const;
		
protected:

	// read variables
	Float_t r_Sis00;							///< Sis00 trigger flags
	Float_t r_TriggerNumber;					///< event trigger number
	BlindTime r_Clk;							///< event time blinded clock
	Float_t r_BClk;								///< proton beam clock
	Float_t r_Delt0;							///< time since last event
	Float_t r_AbsTime;							///< absolute time clock
	Float_t r_PMTADC[BOTH][nBetaTubes];			///< PMT ADCs
	Float_t r_PMTTDC[BOTH][nBetaTubes+1];		///< PMT TDCs
	Float_t r_MWPC_caths[BOTH][2][kMaxCathodes];///< cathodes on [side][xplane][wire]
	Float_t r_MWPC_anode[BOTH];					///< MWPC anode on each side
	Float_t r_MonADC[kNumUCNMons];				///< UCN monitor ADCs
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
	Float_t fAbsTimeStart;						///< absolute start time of run
	Float_t fAbsTimeEnd;						///< absolute end time of run
	BlindTime totalTime;						///< total running time (accumulated from last event); after scanning, total ``live'' time (less global cuts)
	BlindTime deltaT;							///< time scaler wraparound fix
	RangeCut ScintSelftrig[BOTH];				///< self-trigger range cut for each scintillator (used for Type I origin side determination)
	std::vector< std::pair<double,double> > manualCuts;	///< manually cut time segments
	std::vector<std::string> cathNames[BOTH][2];///< cathode sensor names on each [side][xplane]

	bool ignore_beam_out;						///< whether to ignore long beam outages (e.g. for a source run)

	// event-by-event calibrated variables and cuts
	int iTriggerNumber;							///< trigger number
	BlindTime fTimeScaler;						///< absolute event time scaler, blinded E, W, and unblinded
	CutVariable fBeamclock;						///< time since last beam pulse scaler
	Float_t fDelt0;								///< time since previous event
	CutVariable fWindow;						///< time window between previous and next event
	CutVariable fScint_tdc[BOTH][nBetaTubes+1];	///< TDC readout for each PMT and side
	CutVariable fMWPC_anode[BOTH];				///< anode ADC
	CutVariable fBacking_tdc[BOTH];				///< muon backing veto TDC
	CutVariable fDrift_tac[BOTH];				///< muon veto drift tubes TAC
	CutVariable fTop_tdc[BOTH];					///< top veto TDCs (only East)
	CutVariable fMonADC[kNumUCNMons];			///< UCN monitor ADCs = GV, Sw, Fe, SCS
	CutVariable fCathSum[BOTH];					///< combined x+y cathode sum on each side
	CutVariable fCathMax[BOTH];					///< min(max cathode each plane) for each side
	CutVariable fCathMaxSum[BOTH];				///< sum of max cathode from each plane for each side

	/// set to read in timing variables
	void readInTiming();
	/// set to read in wirechamber variables
	void readInWirechambers();
	/// set to read in PMT ADC variables
	void readInPMTADC();
	/// set to read in UCN monitors
	void readInUCNMon();
	/// set to read in muon veto data
	void readInMuonVetos();
	/// set to read in header check variables
	void readInHeaderChecks();

	/// process timing variables
	virtual void calibrateTimes();
};

#endif
