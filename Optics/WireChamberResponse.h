
#include <fstream>
#include <string>
#include <iostream>
#include <stdlib.h>
//#include <stdio.h>
//#include <unistd.h> 
#include <math.h>
#include <climits>
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
#include <gsl/gsl_sf_ellint.h> //double gsl_sf_ellint_E(double phi, double k, gsl_mode_t mode);


class WireChamberResponse{
private:
        std::string privatestring="DO NOT, I repeat DO NOT read this string, its private";
public:
	//look at all of these public (government funded) types!
	char * ORD=getenv ("UCNAOUTPUTDIR"); //Official Replay Data 
	Float_t threshold=120;    //what is the event threshold on the wires. 
	Float_t threshold2=82;   //what is the event threshold on the neighbors. 	
	Float_t cathwx[16];   //cathode data
	Float_t platfrac=0.9;
	Float_t trifrac=1.5;	
	Float_t cathwy[16];
	Float_t cathex[16];
	Float_t cathey[16];
	Float_t tempcath[16];
	Float_t scintE[16];    //scintillation assumed LED corrected. 
   	Float_t scintW[16];    //West Scint (these are not currently used)
	int entnum;	 //number of entries in the tree. 
	TFile *physfile;    //TFile
	TTree *phys;     //MPMs data tree. 
	int resptype[4];   //Response classification 	
	Float_t max;	   //single value max
	Float_t maxind;      //single value max
	int triind[3];      //multi value max
	Float_t trimax[3];    //I pooped a little just now.
	int quadind[4];     //what if its a 4 fold platue??
	Float_t quadmax[4];   //then we should store it as well. but 5 is getting thrown out. 
	Float_t quind[4];
	Int_t wcpos;	    //sorted quad ind. 	


	int ResponseType(Float_t cath[]);   //main function of the class

	//sub functions. used to make the main function. 
	bool FilExists(const std::string& name);	
	void SetPhysTree(int rnum);
	int MaxInd(Float_t cath[]);   //finds max and ind, used by tri ind and quad ind. 
	void SetCaths(int cathdex);
	void SetTempCaths(Float_t cath[]); //sets the cathode
	void TriMax(Float_t cath[]);  //sets the trimax and triind values, it is not used anywhere,
	void QuadMax(Float_t cath[]); //sets the quadmax and quadind values, 
};














