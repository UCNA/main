//#include <fstream>
#include <iostream>
#include <stdlib.h>
//#include <stdio.h>
//#include <unistd.h> 
#include <math.h>
#include <cmath> 
#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"
#include "TTree.h"
//#include "TTreeReader.h"          /// Later versions of root
//#include "TTreeReaderValue.h"     /// Later versions of root
#include "TRandom3.h"
//#include "/home/cmswank/Documents/ucna/main/Scripts/lgprobmap/pmtprobstuff.h"
//#include "/home/cmswank/Documents/ucna/main/Scripts/lgprobmap/lgpmtTools.h"
#include "pmtprobstuff.h"
#include "lgpmtTools.h"


using namespace std;

	
////PMTMapGen operations
//import wirechamber position from root file.
//create lightguide probability from positions.
//export probability to new root file. 
 	
void PMTProbGen(string infile, string outfile, double LGFitParamE[],double LGFitParamW[]){	
   const double pi = 3.1415926535897;  //mine is pumpkin
   Double_t pmtprob[8];

//open mpm's analyzed root file, then load a tree
   TFile *myFile2 = TFile::Open(infile.c_str());
   TTree *anaTree = (TTree*)myFile2->Get("anaTree");
//create new file for the map to go into 
   TFile myFile(outfile.c_str(),"recreate");
   TTree *pmtTree=new TTree("pmtTree","pmt map");

	//add 	 Branch
   TBranch *PMTmap = pmtTree->Branch("PMTmap",&pmtprob,"PMTmap[8]/D");
	//define position array
   Double_t hitpos[6];
    ////get number of entries. 
   Int_t linum=(Int_t)anaTree->GetEntries();
	//set position branch addres to position array   
	anaTree->SetBranchAddress("MWPCPos",&hitpos);
	
  	double xpos;
	double ypos;
	//double tlgp[linum];  
	double both;	
    	int ii = 0;
         int west = 0;
    	while (ii<linum)
    	{   
	anaTree->GetEntry(ii);	
	///what if its both?? Not ignoring for now. 
	    
	    //west side or east side?
	    if (hitpos[0]!=0 && hitpos[3]!=0){
		cout<<ii<<"\n";//this happens a lot.  cout<<"Both PMT triggered something weird happened";
		both=1;
		}
	    if (hitpos[3]!=0){  
	    west = 1; //west side
     	    xpos = 10*hitpos[3];  //change units from cm to mm (10*hitpos). 
  	    ypos = 10*hitpos[4];
	    }	
   	    if (hitpos[0]!=0)
            {
	    west=0;//east side
	    xpos = 10*hitpos[0];
  	    ypos = 10*hitpos[1];	
             }
	for (int i=0; i<8;++i) pmtprob[i]=0.000;         //clear pmtprob (stupid way)

	///get PMTmap value. uses function PMTprob... 
	if (west==1 || both==1){	
	        for (int i=0; i<4; ++i) pmtprob[4+i]=PMTprob(xpos,ypos,LGFitParamW,i+1); 
	      }
         if(west==0 || both==1){ 
		for (int i=0; i<4; ++i) pmtprob[i]=PMTprob(xpos,ypos,LGFitParamE,i+1);
		}
	    both=0;
		
		///////
		///PMT Total  Probablility check!  Uncomment the non-comment lines to check.
		//tlgp[ii]=PMTprob(xpos[ii],ypos[ii],LGFitParam,1)+PMTprob(xpos[ii],ypos[ii],LGFitParam,,2)+PMTprob(xpos[ii],ypos[ii],0,0,1,1,1,3)+PMTprob(xpos[ii],ypos[ii],0,0,1,1,1,4); 
		//if (tlgp[ii]>1.005)
		//{      //all errors tend to be less than 1% e.g. prob=1.00X where X is < 5.
			// still its pretty big, thats what we get for using atan2?...
		  //cout<<"Warning, Probabililty>1, something isn't right here.\n";				
		//}
		/////
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
	LGFitParamE[2+i*15]=1;    /// 15 not 12, the pmt location rotates w.r.t the lg designation.  
	LGFitParamE[3+i*15]=1;
	LGFitParamE[4+i*15]=1;
	LGFitParamW[2+i*15]=1;
	LGFitParamW[3+i*15]=1;
	LGFitParamW[4+i*15]=1;
	}

//pass the name and location of the file you want to work and the name and location of the probablilty you want to save, with the intention the new TTree file is added to the anaTree as a friend when calibration of the sources is done. 	
  PMTProbGen("/home/cmswank/G4Work/output/thinfoil_Xe135_3-2+/analyzed_1.root","pmtprob_1.root",LGFitParamE,LGFitParamW);

  return 0;
}

