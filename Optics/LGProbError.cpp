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
#include "TPaveStats.h"
#include "TBox.h"
#include "TText.h"
#include "TH2.h"
#include "TF1.h"
#include "TH2F.h"
#include "TPaveText.h"



using namespace std;

void LGProbError(string infile, string outfile, double LGFitParamE[],double LGFitParamW[])
	{
   const double pi = 3.1415926535897;  //strawberry
   Double_t pmtprob[8];
   Double_t pmterr[8];

//open mpm's analyzed root file, then load a tree
   TFile *myFile2 = TFile::Open(infile.c_str());
   TTree *phys = (TTree*)myFile2->Get("phys");
  
//create new file for the map to go into 
   TFile myFile(outfile.c_str(),"recreate");
   TTree *pmtTree=new TTree("pmtTree","pmt probabiliity for each event");
 
	//add 	 Branch
   TBranch *PMTmap = pmtTree->Branch("PMTmap",&pmtprob,"PMTmap[8]/D");
   TBranch *PMTerr = pmtTree->Branch("PMTerr",&pmterr,"PMTerr[8]/D");
  
	//define position array (center is [0] width is [1])
   Float_t xEpos[10];
   Float_t yEpos[10];
   Float_t xWpos[10];
   Float_t yWpos[10];
 
    ////get number of entries. 
   Int_t linum=(Int_t)phys->GetEntries();
		
	//set position branch addres to position array   
	phys->SetBranchAddress("xEmpm",&xEpos);
	phys->SetBranchAddress("yEmpm",&yEpos);
	phys->SetBranchAddress("xWmpm",&xWpos);
	phys->SetBranchAddress("yWmpm",&yWpos);

	
  	
	//double tlgp[linum];  

    	int ii = 0;
    
    	while (ii<linum)
    	{   
	phys->GetEntry(ii);	
	///what if its both?? Not ignoring for now. 
	    
	       //clear pmtprob (stupid way)
	for (int i=0; i<8;++i) {pmtprob[i]=0;pmterr[i]=0;}    

	///get PMTmap value. uses function PMTprob... 
		for (int i=0; i<4; ++i){ //east side 		
		pmtprob[i]=PMTprob(xEpos[0],yEpos[0],LGFitParamE,i+1);
		pmterr[i]=PMTerror(xEpos[1],yEpos[1],xEpos[0],yEpos[0],LGFitParamE,i+1);
		}
	    		
	        for (int i=0; i<4; ++i) { //west side
		  pmtprob[4+i]=PMTprob(xWpos[0],yWpos[0],LGFitParamW,i+1);
      		pmterr[4+i]=PMTerror(xWpos[1],yWpos[1],xWpos[0],yWpos[0],LGFitParamW,i+1);
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
		pmtTree->Fill(); //fill the branch  (How does it know what entry I'm on? Magic!?!?)
		++ii;	 
	
  	   }   
 	pmtTree->Write(); //write the file. 
   	               
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

//pass the name and location of the file you want to work and the name and location of the probablilty you want to save, with the intention the new TTree file is added to the anaTree as a friend when calibration of the sources is done. 	
  LGProbError("$UCNAOUTPUTDIR/hists/spec_22770.root","pmtprob_1.root",LGFitParamE,LGFitParamW);

  return 0;
}


