#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h> 
#include <stdlib.h>
#include <math.h>
#include <cmath> 
#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "/home/cmswank/Documents/ucna/lgprobmap/pmtprobstuff.h"
#include "/home/cmswank/Documents/ucna/lgprobmap/lgpmtTools.h"

using namespace std;
 
//This is a collection of functions that are useful for light guide map creation. 
///For an example of use look in main. 



	
////the main function:
//import wirechamber position from root file.
//create lightguide probability from positions.
//export probability to new root file. 
 	
int main()
{	
	const double pi = 3.1415926535897;  //pi
	Double_t pmtprob[8];
	
	

//open mpm's analyzed root file, then load a tree
   TFile *myFile2 = TFile::Open("/home/cmswank/G4Work/output/thinfoil_Xe135_3-2+/analyzed_1.root");
   TTree *anaTree = (TTree*)myFile2->Get("anaTree");
//create new file for the map to go into 
   TFile myFile("pmtprob_1.root","recreate");
   TTree *pmtTree=new TTree("pmtTree","pmt map");

	//add 	 Branch
   TBranch *PMTmap = pmtTree->Branch("PMTmap",&pmtprob,"PMTmap[8]/D");
	//define position array
   Double_t hitpos[6];
    ////get number of entries. 
   Int_t linum=(Int_t)anaTree->GetEntries();
	//set position branch addres to position array   
	anaTree->SetBranchAddress("MWPCPos",&hitpos);



	
	string lgps;
  	double xpos[linum];
	double ypos[linum];
	double tlgp[linum];  	
    	int ii = 0;
         int west = 0;
    	while (ii<linum)
    	{   
	    anaTree->GetEntry(ii);
	    //west side or east side?
	    if (hitpos[3]==0){  
	    west = 1; //west side
     	    xpos[ii] = 10*hitpos[4];  //change units from cm to mm (10*hitpos). 
  	    ypos[ii] = 10*hitpos[5];}	
   	    else
            {
	    west=0;//east side
	    xpos[ii] = 10*hitpos[1];
  	    ypos[ii] = 10*hitpos[2];	
             }
	for (int i=0; i<8;++i) pmtprob[i]=0.000;         //clear pmtprob (stupid way)

///eventually I will need to use coupling coefficients/offsets.
	for (int i=0; i<4; ++i) pmtprob[4*west+i]=PMTprob(xpos[ii],ypos[ii],0,0,1,1,1,i+1); //
            //lgp[ii] = lightguideprob(xpos[ii],ypos[ii],0); ///LG Probability function
	    tlgp[ii]=PMTprob(xpos[ii],ypos[ii],0,0,1,1,1,1)+PMTprob(xpos[ii],ypos[ii],0,0,1,1,1,2)+PMTprob(xpos[ii],ypos[ii],0,0,1,1,1,3)+PMTprob(xpos[ii],ypos[ii],0,0,1,1,1,4); ///PMT Total Probablility check!
		
		if (tlgp[ii]>1.05)
		{      //all errors tend to be less than 1% e.g. prob=1.00X where X is < 5.
			// still its pretty big, thats what we get for using atan2?...
		  cout<<"Warning, Probabililty>1, something isn't right here.\n";
	          				
		}
		pmtTree->Fill(); //fill the branch  (How does it know what entry I'm on? Magic!?!?)
		++ii;	 	
  	   }   
 	pmtTree->Write(); //write the file. 
 	                 



   return 0;

}
