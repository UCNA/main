//g++ -o lgprob LGProbError.cpp pmtprobstuff.cpp lgpmtTools.cpp `root-config --cflags --glibs`
//#include <fstream>
#include <iostream>
#include <stdlib.h>
//#include <stdio.h>
//#include <unistd.h> 
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



using namespace std;

void LGProbError(string infile, string outfile, double LGFitParamE[],double LGFitParamW[])
	{
   const double pi = 3.1415926535897;  //strawberry
   Double_t lgEprob[12];
   Double_t lgWprob[12];
   Double_t lgEerr[12];
   Double_t lgWerr[12];
   Double_t fp1E[12];
   Double_t fp1W[12];
   Double_t fp2E[12];
   Double_t fp2W[12];
   Double_t fp3E[12];
   Double_t fp3W[12];
   Double_t fp4E[12];
   Double_t fp4W[12];
   Double_t offsetE[2];
   Double_t offsetW[2];
   Float_t xEpos[50];//needs to have a lot of space for some reason or another.
   Float_t yEpos[50];
   Float_t xWpos[50];
   Float_t yWpos[50];
   Float_t scintE[4];   //screw adding tree friends, this might make the code less sea sick. 
   Float_t scintW[4];   //phew I'm exhausted 


//open mpm's analyzed root file, then load a tree
   TFile *myFile2 = TFile::Open(infile.c_str());
   TTree *phys = (TTree*)myFile2->Get("phys");
    ////get number of entries. 
   Int_t linum=(Int_t)phys->GetEntries();
		
	//set position branch addres to position array   
	phys->SetBranchAddress("xEmpm",&xEpos);
	phys->SetBranchAddress("yEmpm",&yEpos);
	phys->SetBranchAddress("xWmpm",&xWpos);
	phys->SetBranchAddress("yWmpm",&yWpos);
	phys->SetBranchAddress("ScintE",&scintE);
	phys->SetBranchAddress("ScintW",&scintW);
	
  
//create new file for the map to go into 
   TFile myFile(outfile.c_str(),"recreate");
	//Light guide and PMT Tree
   TTree *pmtTree=new TTree("pmtTree","Light guide probability and qADC for each event");
 
	//pmtTree Branches
    TBranch *LGEprob = pmtTree->Branch("LGEprob",&lgEprob,"LGEprob[12]/D");
    TBranch *LGWprob = pmtTree->Branch("LGWprob",&lgWprob,"LGWprob[12]/D");
    TBranch *LGEerror = pmtTree->Branch("LGEerror",&lgEerr,"LGEerror[12]/D");
    TBranch *LGWerror = pmtTree->Branch("LGWerror",&lgWerr,"LGWerror[12]/D");
    TBranch *ScintE = pmtTree->Branch("ScintE",&scintE,"ScintE[4]/F");
    TBranch *ScintW = pmtTree->Branch("ScintW",&scintW,"ScintW[4]/F");
    
	///Fit parameters Tree
    TTree *fitTree=new TTree("fitTree","LightGuide Map Fit Parameters");
	//Fit Parameter Branches
    TBranch *FitParam1E = fitTree->Branch("FitParamE1",&fp1E,"FitParamE1[12]/D");
    TBranch *FitParam1W = fitTree->Branch("FitParamW1",&fp1W,"FitParamW1[12]/D");
    TBranch *FitParam2E = fitTree->Branch("FitParamE2",&fp2E,"FitParamE2[12]/D");
    TBranch *FitParam2W = fitTree->Branch("FitParamW2",&fp2W,"FitParamW2[12]/D");
    TBranch *FitParam3E = fitTree->Branch("FitParamE3",&fp3E,"FitParamE3[12]/D");
    TBranch *FitParam3W = fitTree->Branch("FitParamW3",&fp3W,"FitParamW3[12]/D");
    TBranch *FitParam4E = fitTree->Branch("FitParamE4",&fp4E,"FitParamE4[12]/D");
    TBranch *FitParam4W = fitTree->Branch("FitParamW4",&fp4W,"FitParamW4[12]/D");
    TBranch *OffsetE = fitTree->Branch("OffsetE",&offsetE,"OffsetE[2]/D");
    TBranch *OffsetW = fitTree->Branch("OffsetW",&offsetW,"OffsetW[2]/D");

	//define initial fit parameters	
    offsetE[0]=LGFitParamE[0];
    offsetE[1]=LGFitParamE[1];
    offsetW[0]=LGFitParamW[0];
    offsetW[1]=LGFitParamW[1];
	for(int i =0; i<12;++i){	
	    fp1E[i]=LGFitParamE[2+i];
	    fp2E[i]=LGFitParamE[2+12+i];
	    fp3E[i]=LGFitParamE[2+24+i];	
       	    fp4E[i]=LGFitParamE[2+36+i];	
  	    fp1W[i]=LGFitParamW[2+i];
	    fp2W[i]=LGFitParamW[2+12+i];
	    fp3W[i]=LGFitParamW[2+24+i];	
	    fp4W[i]=LGFitParamW[2+36+i];
	   }
	
	//fill initial fit parameters to TTree 
	fitTree->Fill();
  	
	//double tlgp[linum]; // for probability checking

    	int ii = 0;
     while (ii<linum)
    	{   
	phys->GetEntry(ii);	

	    
	///get light guide probability value. 
		for (int i=0; i<12; ++i){ //east side 		
		lgEprob[i]=lightguideprob(xEpos[0]-LGFitParamE[0],yEpos[0]-LGFitParamE[1],i);
		lgEerr[i]=LGerror(xEpos[1],yEpos[1],xEpos[0]-LGFitParamE[0],yEpos[0]-LGFitParamE[1],i);
		}
	    		
	        for (int i=0; i<12; ++i) { //west side
		  lgWprob[i]=lightguideprob(xWpos[0]-LGFitParamW[0],yWpos[0]-LGFitParamW[1],i);
      		  lgWerr[i]=LGerror(xWpos[1],yWpos[1],xWpos[0]-LGFitParamW[0],yWpos[0]-LGFitParamW[1],i);
		}
		
		/*
		///PMT Total  Probablility check!  Uncomment the non-comment lines to check.
		tlgp[ii]=PMTprob(xpos[ii],ypos[ii],LGFitParam,1)+PMTprob(xpos[ii],ypos[ii],LGFitParam,,2)+PMTprob(xpos[ii],ypos[ii],0,0,1,1,1,3)+PMTprob(xpos[ii],ypos[ii],0,0,1,1,1,4); 
		if (tlgp[ii]>1.005)
		{      //all errors tend to be less than 1% e.g. prob=1.00X where X is < 5.
			// still its pretty big, thats what we get for using atan2?...
		  cout<<"Warning, Probabililty>1, something isn't right here.\n";				
		}
		*/
		pmtTree->Fill(); //fill the branch  (How does it know what entry I'm on? Magic!?!? it probably just appends entries.)
		//if (ii==19) {cout<<xEpos[0]<<"  "<<yEpos[0]<<"\n";
		//	     cout<<xWpos[0]<<"  "<<yWpos[0]<<"\n";}  //TESTING
		
		++ii;	 
	
  	   }   
 	pmtTree->Write(); //write the file. 
   	fitTree->Write();               
   return;
}
	

int main()
{
  double LGFitParamE[50]; //2 offsets + 12*4 coupling coefficients 
  double LGFitParamW[50]; //2 offsets + 12*4 coupling coefficients
///this is bad. eventually write a Class that writes this info to a file maybe, maybe a TTree?

	 
	for(int i = 0; i<50;++i) {
	LGFitParamE[i]=0;
	LGFitParamW[i]=0;
	}   //are c++ arrays zero by default?? or is it random bits 
        
	for (int i = 0; i <4; ++i){
	LGFitParamE[2+i*15]=1;    /// 15 because the pmt location rotates but the lg designation does not.  
	LGFitParamE[3+i*15]=1;
	LGFitParamE[4+i*15]=1;
	LGFitParamW[2+i*15]=1;
	LGFitParamW[3+i*15]=1;
	LGFitParamW[4+i*15]=1;
	}

//pass the name and location of the file you want to work and the name and location of the probablilty you want to save, with the intention the new TTree file is added to phys Tree as a friend when calibration of the sources is done. 	
  LGProbError("$UCNAOUTPUTDIR/hists/spec_22770.root","$UCNAOUTPUTDIR/hists/pmtprob_22770.root",LGFitParamE,LGFitParamW);

  return 0;
}


