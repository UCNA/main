#ifndef LEDSCANS_HH
#define LEDSCANS_HH 1

#include "RunSetScanner.hh"
#include "Enums.hh"

/// Class for processing LED scan data files
class LEDScanScanner: public RunSetScanner {
public:
	/// constructor
	LEDScanScanner(): RunSetScanner("LED",true) {}
	
	/// find path to processed run .root file
	virtual std::string locateRun(RunNum r);
	
	/// locate jumps indicating scan starts
	std::vector<unsigned int> findJumps(float emin=100., float emax=1000., Side s = EAST);
	
	ScintEvent scints[2];		//< readout point for scintillator data
	Float_t led_pd[2];			//< readout point for reference photodiode
	Float_t runClock;			//< time of current event since run start
	Float_t anode[2];			//< wirechamber anode on each side
	Float_t cathMax[2];			//< wirechamber cathode on each side
	
	std::vector<float> ledAvg[2][nBetaTubes+1];	//< LED smoothed average
	std::vector<float> ledRMS[2][nBetaTubes+1];	//< LED rolling RMS
	/// generate rolling averages for events over event range
	void analyzeSegment(unsigned int estart, unsigned int eend, unsigned int w = 50);
	
protected:
	/// set up tree read points
	virtual void setReadpoints();
};

/// study correlations between PMTs using LED data
void PMT_LED_Correlations(OutputManager& OM, LEDScanScanner& LSS);

/// tests on spectrum generation
void spectrumGenTest();

#endif
