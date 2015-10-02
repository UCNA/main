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

	//this class to create a spectrum for all the different WC response types grid. 
			

class PMTSpectrumMap{
private:
        std::string privatestring="Stop reading my private strings. Computers have rights too, you know.";
public:
	char * ORD=getenv ("UCNAOUTPUTDIR"); //Official Replay Data 
	Int_t rnum; //desired run number
	Int_t runnum;//TTree Value of LEDphys runnumber //current run number.	
	Int_t entnum;	//number entries in LEDphys. 
	
	
	int xindsize=16;
	int yindsize=16;
	int typeindsize=7;
	int mapentries=16*16*7*7*8;
	TH1F * specmap[16*16*7*7*8];	
	TH1F * spec0[7*7*8];
	
	TFile * SpecMapfile;
	TTree * SpecMap;
	TFile * LEDphysfile;
	TTree * LEDphys;
	Float_t specmax=1000;
	Float_t fitrangehi=1000;
	Float_t fitrangelow=0;
	Float_t pmt[8];
	Float_t pmtled[8];
	Float_t pmtlederror[8];
	Float_t wcmax[16]; 
	Int_t wcind[16];
	Float_t wcresp[4]; 
	Int_t wcpos[4]; 


		 ///write to  position map (flattened array)
	void FillSpectrumMap();
	bool FilExists(const std::string& name);
	void SetPhysTree(int rnum);
	void FitSpectrum(TH1F * spec0,TH1F * spectemp);
	void FindSpecZero();	
	int ArrayInd(int xind,int yind, int typex, int typey, int pmt);
};



	

