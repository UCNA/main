#ifndef CALDBFake_HH
#define CALDBFake_HH 1

#include "CalDB.hh"
#include "QFile.hh"
#include "strutils.hh"
#include <TGraphErrors.h>
#include <TF1.h>
#include <map>

/// class for retrieving calibration data from DB
class CalDBFake: public CalDB {
public:
	/// constructor
	CalDBFake(): calfile("../SummaryData/Cal_Default.txt") {}
	/// destructor
	~CalDBFake() {}
	
	/// check if valid data available for run
	bool isValid(RunNum rn) { return true; }
	
	/// get linearity correction data
	virtual TGraph* getLinearity(RunNum rn, Side s, unsigned int t) { return constGraph(0.0,1.0); }
	
	/// get noise estimate calibration point ADC width
	virtual float getNoiseWidth(RunNum rn, Side s, unsigned int t)  { return calfile.getFirst("PMT_"+sideNames(s)+itos(t)).getDefault("noiseWidth",50); }
	/// get noise estimate calibration point raw ADC value
	virtual float getNoiseADC(RunNum rn, Side s, unsigned int t)  { return calfile.getFirst("PMT_"+sideNames(s)+itos(t)).getDefault("noiseADC",400); }
	
	/// get energy calibration point ADC value
	virtual float getEcalADC(RunNum rn, Side s, unsigned int t)  { return calfile.getFirst("PMT_"+sideNames(s)+itos(t)).getDefault("ecalADC",500); }
	/// get energy calibration point visible energy
	virtual float getEcalEvis(RunNum rn, Side s, unsigned int t)  { return calfile.getFirst("PMT_"+sideNames(s)+itos(t)).getDefault("ecalEvis",500); }
	/// get energy calibration point x position
	virtual float getEcalX(RunNum rn, Side s)  { return calfile.getFirst(sideSubst("PMT_%c",s)+"1").getDefault("ecalX",0); }
	/// get energy calibration point y position
	virtual float getEcalY(RunNum rn, Side s) { return calfile.getFirst(sideSubst("PMT_%c",s)+"1").getDefault("ecalY",0); }
	/// get GMS Cal run number
	virtual RunNum getGMSRun(RunNum rn) { return 0; } //(int)calfile.getFirst("PMT_E0").getDefault("gmsRun",0); }
	/// get a run monitor graph
	TGraphErrors* getRunMonitor(RunNum rn, const std::string& sensorName, const std::string& monType, bool centers = true) { return NULL; }
	/// get start of run monitor graph
	float getRunMonitorStart(RunNum rn, const std::string& sensorName, const std::string& monType) { return 1.0; }
	/// get a continuous monitor graph for a run
	TGraphErrors* getContinuousMonitor(const std::string& sensorName, const std::string& monType, RunNum rn, bool centers) { return NULL; }
	
	/// get pedestals for named sensor
	virtual TGraph* getPedestals(RunNum rn, const std::string& sensorName) { return NULL; }
	/// get pedestal widths for named sensor
	virtual TGraph* getPedwidths(RunNum rn, const std::string& sensorName) { return NULL; }
	
	/// get kurie endpoint calibration energy for given run, side, PMT
	virtual float getKurieEnergy(RunNum rn, Side s, unsigned int t) { return 782.0; }
	/// get kurie endpoint calibration energy for given run, side, PMT
	virtual float getKurieADC(RunNum rn, Side s, unsigned int t) { return 782.0; }
	
	/// get positioning corrector for given run
	virtual PositioningCorrector* getPositioningCorrector(RunNum rn) { return getPositioningCorrectorByID(0); }
	/// get anode positioning corrector for given run
	virtual PositioningCorrector* getAnodePositioningCorrector(RunNum rn) { return getPositioningCorrectorByID(0); }
	/// get anode gain correction factor
	virtual float getAnodeGain(RunNum rn, Side) { return 1.0; }
	/// get GMS gain tweaking factors
	virtual void getGainTweak(RunNum rn, Side s, unsigned int t, float& orig, float& final) { orig=final=500.0; }
	
	/// get positioning corrector by ID number
	PositioningCorrector* getPositioningCorrectorByID(unsigned int psid) {
		printf("*** Warning: loading default fake position map for %i!\n",psid);
		QFile qin("../SummaryData/Posmap_Default.txt");
		return new PositioningCorrector(qin);	
	}
	
	/// get trigger efficiency function
	EfficCurve* getTrigeff(RunNum rn, Side s, unsigned int t) { return NULL; /*TODO*/ }
	/// Evis to Etrue conversion
	TGraph* getEvisConversion(RunNum rn, Side s, EventType tp) { return NULL; /*TODO*/ }
	/// get list of cathode segment calibrators (caller is responsible for deletion)
	virtual std::vector<CathSegCalibrator*> getCathSegCalibrators(RunNum rn, Side s, AxisDirection d) { return std::vector<CathSegCalibrator*>(); /*TODO*/ }	
	
	/// run start time
	int startTime(RunNum rn, int t0 = 0) { return t0; }
	/// run end time
	int endTime(RunNum rn, int t0 = 0) { return t0; }
	/// run time after cuts, blinded
	BlindTime fiducialTime(RunNum rn) { return BlindTime(600.0); }
	/// get RunInfo for given run
	RunInfo getRunInfo(RunNum r) { return RunInfo(r); }
	
	/// get name
	virtual std::string getName() const { return "Fake DB!"; }
	
protected:
	
	/// make a constant-value TGraphErrors
	TGraphErrors* constGraph(float y, float dy = 0) const {
		TGraphErrors* g = new TGraphErrors(2);
		g->SetPoint(0,0,y);
		g->SetPoint(1,1,y+dy);
		return g;
	}
	
	QFile calfile;	//< calibration defaults file
	
};

#endif

