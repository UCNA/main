#ifndef LEDSCANS_HH
#define LEDSCANS_HH 1

#include "ProcessedDataScanner.hh"
#include "RollingWindow.hh"
#include "Enums.hh"
#include <deque>
#include <vector>
#include <TGraphErrors.h>

/// Class for processing LED scan data files
class LEDScanScanner: public ProcessedDataScanner {
public:
	/// constructor
	LEDScanScanner(): ProcessedDataScanner("LED",true) {}
	/// find path to processed run .root file
	virtual std::string locateRun(RunNum r);
	
protected:
	/// set up tree read points
	virtual void setReadpoints();
};




class SweepAverager;

/// Class for caching points to compare to averaged values
class AveragerUnit {
public:
	/// constructor
	AveragerUnit(unsigned int n): npts(n), RW(2*n+1) {}
	/// destructor
	virtual ~AveragerUnit() {}
	
	/// get current cached point
	double getP() const { return pts.front(); }
	/// get average
	double getAvg() const { return RW.getAvg(); }
	/// get point-excluded average
	double getAvgExcl() const { return RW.getAvgExcl(getP()); }
	
	double x;	//< point to add to average
	

protected:

	friend class SweepAverager;
	
	/// clear averaging data
	void clear();
	
	const unsigned int npts;	//< half window size
	RollingWindow RW;			//< rolling window with average
	std::deque<double> pts; 	//< cached points for later evaluation
};

/// Class for handling and evaluating multiple averagers
class SweepAverager {
public:
	/// constructor
	SweepAverager(unsigned int n): npts(n), ncached(0) {}
	/// destructor
	virtual ~SweepAverager() { while(AUs.size()) { delete(AUs.back()); AUs.pop_back(); } }
	
	/// start averaging
	void start();
	/// evaluate all units to process point
	void addPt();
	/// stop averaging
	void stop();
	
	/// add a new AveragerUnit
	AveragerUnit* addUnit() { AUs.push_back(new AveragerUnit(npts)); return AUs.back(); }
	
	/// subclass this for actions on averaged data
	virtual void doWithAverage() {}
	
protected:
	const unsigned int npts;		//< half window width
	unsigned int ncached;			//< counter for cached points
	std::vector<AveragerUnit*> AUs;	//< individual averaging units
};






class LEDAnalyzer: public OutputManager, SweepAverager {
public:
	/// constructor
	LEDAnalyzer(std::string nm, std::string bp);
	/// destructor
	virtual ~LEDAnalyzer() {}
	
	/// collect data from specified scanner
	void ScanData(ProcessedDataScanner& LSS);
	/// calculate correlations from collected data
	void CalcCorrelations();
	/// make trigger efficiency curves
	void CalcTrigEffic();
	/// calculate light balance
	void CalcLightBal();
	/// set up MLP learning set for trigger efficiency
	void buildPerceptronTree();
	
protected:
	/// locate jumps indicating scan starts
	std::vector<unsigned int> findJumps(ProcessedDataScanner& LSS, float emin=100., float emax=1000., Side s = EAST);
	/// calculate correlations between PMTs based on LED
	void PMT_LED_Correlations();
	
	AveragerUnit* Eavg;						//< total energy averager
	AveragerUnit* clockAvg;					//< access to event time data
	AveragerUnit* Tavg[2][nBetaTubes];		//< PMT energy averager
	/// fill histograms from average
	virtual void doWithAverage();
	
	float emax;											//< maximum energy for plot ranged
	float wmax;											//< maximum deviation from mean to consider
	
	unsigned int pass;					//< analysis pass
	std::vector<unsigned int> jumps;	//< start/end of scan segments
	/// perform data pass over segments
	void dataPass(ProcessedDataScanner& LSS);
	/// scan one range of points
	void ScanSeg(ProcessedDataScanner& LSS, unsigned int eStart, unsigned int eEnd);
	/// make plot of one segment
	void TimePlot(ProcessedDataScanner& LSS, unsigned int jn);
	
	TProfile* pE0;						//< master energy correction profile
	TProfile* pEi[2][nBetaTubes];		//< individual PMT rolling average deviations
	TH1F* hPed;							//< low-energy histogram for locating pedestal energy
	float EZero;						//< actual pedestal location on energy scale
	TH1F* hAvgEnergy;					//< high-resolution averaged energy variable
	TGraphErrors* gE0;					//< Interpolatable TGraph version of profile
	TGraphErrors* gEi[2][nBetaTubes];	//< Interpolatable TGraph version of profile
	
	TH2F* hE8;											//< combined energy spread
	TProfile* pE8;										//< profile for combined energy spread
	TH2F* corrs[2][nBetaTubes+1][nBetaTubes+1];			//< correlations histogram between each PMT pair
	TProfile* corrsProf[2][nBetaTubes+1][nBetaTubes+1];	//< TProfile of correlations on same side
	TProfile* corrsProfEW[nBetaTubes+1][nBetaTubes+1];	//< TProfile of correlations between sides
	
	// time domain plot
	TGraph* gRawVT;
	TGraph* gAvgVT;
	TGraph* gCorVT;
	unsigned int nVT;
				
	TH1F* hLowE[2][nBetaTubes+1][2];					//< low-energy events for each [side][PMT][triggered?]
	TH1F* hLightBal[2][nBetaTubes];						//< PMT ``light balance'' from LED
	TH1F* h2foldProb[2][2];								//< calculated 2-fold probability for [side][triggered]
	TTree* MLPdat[2];									//< tree with data for MLP for triggers on each side
	TriggerProb TP[2];									//< plain trigger probability calculator
	Float_t ewt;										//< event weight for trigger tree
	Float_t evtE[2];									//< event energy on each side
	Float_t is2fold[2];									//< whether event was 2-fold trigger, for training MLP output
};

/// tests on spectrum generation
void spectrumGenTest();

/// Trigger efficiency neural network fit
void makeMLPfit(OutputManager& OM);


#endif
