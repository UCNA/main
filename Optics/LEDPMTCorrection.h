#include "WireChamberResponse.h"
#include <string>
//#include <sys/stat.h>
#include <fstream>
#include <climits>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <sstream>
#include <unistd.h> 
#include <math.h>
#include <cmath> 
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
//#include "TTreeReader.h"
//#include "TTreeReaderValue.h"
#include "TMath.h"
#include "TRandom3.h"
//#include "/home/cmswank/Documents/ucna/main/Scripts/lgprobmap/pmtprobstuff.h"
//#include "/home/cmswank/Documents/ucna/main/Scripts/lgprobmap/lgpmtTools.h"
#include "pmtprobstuff.h"
#include "lgpmtTools.h"
#include "TLatex.h"
#include "TFormula.h"
#include "TAxis.h"
#include "TPave.h"
#include "TPaveStats.h" //Loading MPM's TPhys TTree requires a lot of includes... maybe there is a TPhys.h?
#include "TBox.h"      //
#include "TText.h"    //   
#include "TH2.h"
#include "TF1.h"
#include "TH2F.h"
#include "TPaveText.h"
#include "FlattenPMTMap.hh"
#include <gsl/gsl_mode.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_sf_ellint.h> 

	//this class is designed to take a maxtrix with LED quadratic fits (in simon's format but turned into a matrix) and apply them to the data;
			

class LEDPMTApp{
private:
        std::string privatestring="HOW DARE YOU! THIS IS PRIVATE!!";
public:
	char * ORD=getenv ("UCNAOUTPUTDIR"); //Official Replay Data 
	Int_t runnum; //current run number.	
	WireChamberResponse *WCR =new WireChamberResponse();  //we need this to fill in the wire chamber response. 
	Float_t led0[8];
	Float_t led1[8];
	Float_t led2[8];
	Float_t lederror0[8];
	Float_t lederror1[8];
	Float_t lederror2[8];
	Float_t PMT[8];
	Float_t PMTLED[8];
	Float_t PMTLEDError[8];
	TFile * physfile;	
	TTree * phys;
		///////////////////This is the TTree with the applied PMT data.////////////////////////////////	
	TFile * LEDphysfile;
	TTree * LEDphys; //new tree with corrected PMT data! branches: Runnum,PMTLED,PMTLEDError,PMT,WCResponse,WCMaxind,WCMax. 
		//////////////////////////////////////////////////
	Float_t scintE[16];    //scintillation assumed LED corrected. 
   	Float_t scintW[16];		
	Int_t entnum;       
	bool import=false;    //did we import from text? or use a root file?
	TFile * LEDfitfile;
		//////////////////This is the fit tree. 
	TTree * LEDfit;  //TTree to store the LED fit matrix in for easy(ish) access. 
	std::string matrixfile;	

	void ApplyCorrection(int runnumber);  //results are found in PMTLED[8] after function is run.

	void SetMatrixFile(std::string filename);
	bool FilExists(const std::string& name);
	void SetPhysTree(int rnum);
	int ImportLED(std::string filename);   ///1 is a successful import, 0 is a failure. also it yells at you if it fails. 
	void LoadLEData(std::string filename);

	
};
