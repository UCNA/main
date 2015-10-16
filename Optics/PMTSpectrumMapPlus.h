#include <TSpline.h>
#include <TCanvas.h>
#include <TGraph2D.h>
#include "Math/Interpolator.h"
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
	char * ORD=getenv ("UCNAOUTPUTDIR"); //Official Replay Data Folder! 
	Int_t rnum; //desired run number
	Int_t runnum;//TTree Value of LEDphys runnumber //current run number.	
	Int_t entnum;	//number entries in LEDphys. 
	
	const double * ROOTsucks;
	const double * ROOTblows;
	unsigned int numpoint=300;
	int xindsize=16;
	int yindsize=16;
	int typeindsize=9;
	int mapentries=16*16*9*9*8;
	TH1F * specmap[16*16*9*9*8];	
	TH1F * specmaperror[16*16*9*9*8];	
	TGraph* tempgraph;
	Float_t PMTERRORMAP[16*16*9*9*8];
	Float_t pmtfiterror;	
	Float_t PMTMAP[16*16*9*9*8];

	Float_t pmttemp[300];
	Float_t pmttemperror[300];
	Float_t pmttempsum=0;
	
	Float_t spec11sum=0;
	Float_t spec12sum=0;
	Float_t spec21sum=0;
	Float_t spec22sum=0;
	Float_t spec13sum=0;
	Float_t spec31sum=0;	
	Float_t spec23sum=0;
	Float_t spec32sum=0;
	Float_t spec33sum=0;

	int	xit[9],yit[9],txit[9],tyit[9],ti[9];//array position of spec0; 	
	//std::vector<double> spec0;  // double vector of spec0 values
	//std::vector<double> spec0x; //x-axis for spec0;
	Double_t spec11[300];
	Double_t spec12[300];
	Double_t spec21[300];
	Double_t spec22[300];
	Double_t spec13[300];
	Double_t spec31[300];
	Double_t spec23[300];
	Double_t spec32[300];
	Double_t spec33[300];
	Double_t spec0x[300];

	TSpline3 *spec11inter;
	TSpline3 *spec12inter;
	TSpline3 *spec21inter;
	TSpline3 *spec22inter;
	TSpline3 *spec13inter;
	TSpline3 *spec31inter;
	TSpline3 *spec23inter;
	TSpline3 *spec32inter;
	TSpline3 *spec33inter;

	//Broke as a joke! ROOT sucks !interpolation class for spec0. 
	//ROOT::Math::Interpolator* spec0inter= new ROOT::Math::Interpolator(300);
	
	
	TFile * SpecMapfile;
	TTree * SpecMap;
	TFile * LEDphysfile;
	TTree * LEDphys;
	Float_t min;
	Float_t minind;
	Float_t specmax=1000;
	Float_t fitrangehi=1000;
	Float_t fitrangelow=0;
	Float_t pmt[8];
	Float_t pmtled[8];
	Float_t pmtlederror[80];
	Float_t wcmax[16]; 
	Int_t wcind[16];
	Float_t wcresp[4]; 
	Int_t wcpos[4]; 
	Float_t chisquared;

		 
	void CreateSpectrumError();	
	void FillSpectrumMap();
	void FillPMTMap();
	void FindSpec0();   ///Find Spec0. 
	void ConstructSpectrumMap();
	void LoadSpectrumMap(std::string filename);
	bool FilExists(const std::string& name);
	void SetPhysTree(int rnum);

	Float_t MinAlpha(int flatind,int resolution,TSpline3 *spec0inter);
	void FindSpecZero();	
	int ArrayInd(int xind,int yind, int typex, int typey, int pmt);
	Float_t ChiSquaredBrah(Float_t alpha,TSpline3 *spec0inter);
};



	

