#ifndef CALDBSQL_HH
#define CALDBSQL_HH 1

#include "SQL_Utils.hh"
#include "CalDB.hh"
#include <TGraphErrors.h>
#include <TF1.h>
#include <map>
#include "PathUtils.hh"
#include "Types.hh"

/// class for retrieving calibration data from DB
class CalDBSQL: public SQLHelper, public CalDB, private NoCopy {
public:
	/// destructor
	~CalDBSQL();
	
	/// check if valid data available for run
	bool isValid(RunNum rn) { return getCalSetInfo(rn,"ecal_id"); }
	
	/// get linearity correction data
	virtual TGraph* getLinearity(RunNum rn, Side s, unsigned int t) { return getGraph( getTubecalInt(rn,s,t,"linearity_graph") ); }
	
	/// get noise estimate calibration point ADC width
	virtual float getNoiseWidth(RunNum rn, Side s, unsigned int t)  { return getTubecalData(rn,s,t,"noisecal_width"); }
	/// get noise estimate calibration point raw ADC value
	virtual float getNoiseADC(RunNum rn, Side s, unsigned int t)  { return getTubecalData(rn,s,t,"noisecal_adc"); }
	
	/// get energy calibration point ADC value
	virtual float getEcalADC(RunNum rn, Side s, unsigned int t)  { return getTubecalData(rn,s,t,"encal_adc"); }
	/// get energy calibration point visible energy
	virtual float getEcalEvis(RunNum rn, Side s, unsigned int t)  { return getTubecalData(rn,s,t,"encal_evis"); }
	/// get energy calibration point x position
	virtual float getEcalX(RunNum rn, Side s)  { return getTubecalData(rn,s,0,"encal_xpos"); }
	/// get energy calibration point y position
	virtual float getEcalY(RunNum rn, Side s) { return getTubecalData(rn,s,0,"encal_ypos"); }
	/// get GMS Cal run number
	virtual RunNum getGMSRun(RunNum rn) { return getCalSetInfo(rn,"gms_run"); }
	/// get a run monitor graph
	TGraphErrors* getRunMonitor(RunNum rn, const std::string& sensorName, const std::string& monType, bool centers = true);
	/// get initial value for run monitor
	float getRunMonitorStart(RunNum rn, const std::string& sensorName, const std::string& monType);
	/// get a continuous monitor graph for a run
	TGraphErrors* getContinuousMonitor(const std::string& sensorName, const std::string& monType, RunNum rn, bool centers);
	
	/// get pedestals for named sensor
	virtual TGraph* getPedestals(RunNum rn, const std::string& sensorName) { return getRunMonitor(rn,sensorName,"pedestal",true); }
	/// get pedestal widths for named sensor
	virtual TGraph* getPedwidths(RunNum rn, const std::string& sensorName) { return getRunMonitor(rn,sensorName,"pedestal",false); }
	
	/// get kurie endpoint calibration energy for given run, side, PMT
	virtual float getKurieEnergy(RunNum rn, Side s, unsigned int t);
	/// get kurie endpoint calibration energy for given run, side, PMT
	virtual float getKurieADC(RunNum rn, Side s, unsigned int t);
	
	/// get positioning corrector for given run
	virtual PositioningCorrector* getPositioningCorrector(RunNum rn);
	/// get positioning corrector by ID number
	PositioningCorrector* getPositioningCorrectorByID(unsigned int psid);
	
	/// get anode positioning corrector for given run
	virtual PositioningCorrector* getAnodePositioningCorrector(RunNum rn);
	/// get anode gain correction factor
	virtual float getAnodeGain(RunNum rn, Side s);
	/// get GMS gain tweaking factors
	virtual void getGainTweak(RunNum rn, Side s, unsigned int t, float& orig, float& final);
	
	/// get trigger efficiency function
	EfficCurve* getTrigeff(RunNum rn, Side s, unsigned int t);
	/// get E_vis -> E_true parametrization
	TGraph* getEvisConversion(RunNum rn, Side s, EventType tp);
	/// get list of cathode segment calibrators (caller is responsible for deletion)
	virtual std::vector<CathSegCalibrator*> getCathSegCalibrators(RunNum rn, Side s, AxisDirection d);
	
	/// run start time
	int startTime(RunNum rn, int t0 = 0);
	/// run end time
	int endTime(RunNum rn, int t0 = 0);
	/// run time after cuts, blinded
	BlindTime fiducialTime(RunNum rn);
	/// total run time
	virtual float totalTime(RunNum rn);
	/// get RunInfo for given run
	RunInfo getRunInfo(RunNum r);
	/// get a list of run numbers matching the given conditions
	std::vector<RunNum> findRuns(const std::string& whereConditions, RunNum r0, RunNum r1);
	
	/// get name
	virtual std::string getName() const { return getDBName(); }
	
	/// globally available CalDB
	static CalDBSQL* getCDB(bool readonly = true);
	
	/// create an ID for a new graph
	unsigned int newGraph(const std::string& description);
	/// delete graph with given ID
	void deleteGraph(unsigned int gid);
	/// upload new graph
	unsigned int uploadGraph(const std::string& description, std::vector<double> x, std::vector<double> y,
							 std::vector<double> dx = std::vector<double>(), std::vector<double> dy = std::vector<double>());
	/// delete a run monitor
	void deleteRunMonitor(RunNum rn, const std::string& sensorName, const std::string& monType);
	/// upload a run monitor
	void addRunMonitor(RunNum rn, const std::string& sensorName, const std::string& monType, unsigned int cgid, unsigned int wgid);
	/// upload a trigger efficiency curve
	void uploadTrigeff(RunNum rn, Side s, unsigned int t, std::vector<double> params, std::vector<double> dparams);	
	/// delete a trigger efficiency curve
	void deleteTrigeff(RunNum rn, Side s, unsigned int t);
	/// list position maps
	void listPosmaps();
	/// add new position map set
	unsigned int newPosmap(const std::string& descrip, unsigned int nrings, double radius);
	/// upload position map point
	void addPosmapPoint(unsigned int pmid, Side s, unsigned int t, unsigned int n, double sig, double norm, double x, double y);
	/// delete a position map
	void deletePosmap(unsigned int pmid);

protected:
	/// constructor (use CalDBSQL::getCDB() if you need access to DB)
	CalDBSQL(const std::string& dbName = getEnvSafe("UCNADB"),
			 const std::string& dbAddress = getEnvSafe("UCNADBADDRESS"),
			 const std::string& dbUser =  getEnvSafe("UCNADBUSER_READONLY"),
			 const std::string& dbPass = getEnvSafe("UCNADBPASS_READONLY"),
			 unsigned int port = atoi(getEnvSafe("UCNADBPORT","3306").c_str())
			 ): SQLHelper(dbName,dbAddress,dbUser,dbPass,port) {}
	/// get sensor ID for given name
	unsigned int getSensorID(const std::string& sname);
	/// get calibration set table info for given run
	unsigned int getCalSetInfo(RunNum R, const char* field);
	/// get tube calibrations table tubecal_id
	unsigned int getTubecalID(RunNum R, Side s, unsigned int t);
	/// get anode calibrations ID
	float getAnodeCalInfo(RunNum R, const char* field);
 	/// get tube calibrations table entry
	virtual float getTubecalData(RunNum rn, Side s, unsigned int t, const char* field);
	/// get tube calibrations table (int) entry
	virtual int getTubecalInt(RunNum rn, Side s, unsigned int t, const char* field);
	/// get graph by ID number
	TGraphErrors* getGraph(unsigned int gid);	
	/// get graph for run time range
	TGraphErrors* getGraph(unsigned int gid, RunNum rn);
	/// get run group name for run
	std::string getGroupName(RunNum rn);
	
	std::map<unsigned int,PositioningCorrector*> pcors;	//< cached positioning correctors
};

#endif

