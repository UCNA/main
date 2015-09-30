//g++ -o WCR WireChamberResponse.cpp `root-config --cflags --glibs`
//g++ -o lgprob LGProbPMTQadc.cpp pmtprobstuff.cpp lgpmtTools.cpp FlattenPMTMap.cc -lgsl -lgslcblas -lm `root-config --cflags --glibs`


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


using namespace std;

//Set File Name       
void WireChamberResponse::SetPhysTree(string filename){

  TFile* myFile = TFile::Open(filename.c_str());
  phys = (TTree*)myFile->Get("phys");
  entnum=phys->GetEntries();
  phys->SetBranchAddress("Cathodes_Wx",&cathwx);
  phys->SetBranchAddress("ScintE",&scintE);
  phys->SetBranchAddress("ScintW",&scintW);
  phys->SetBranchAddress("Cathodes_Ey",&cathey);
  phys->SetBranchAddress("Cathodes_Ex",&cathex);
  phys->SetBranchAddress("Cathodes_Wy",&cathwy);
 
  return;
}

//Find Max index
int WireChamberResponse::MaxInd(Float_t cath[])
{   
	Float_t tempcath[16];
	for(int i =0; i< (16-sizeof(cath)*CHAR_BIT/32); i++)
	{
	tempcath[i]=cath[i];
	}
	
	this->max=-3.40282e+38; //barely floats
	for (int i=0; i<16; ++i)
	{
                    
		if  (tempcath[i]>this->max)
		{
                
		this->max=tempcath[i];
		maxind=i;
		}
	}
	
	return maxind;
}


void WireChamberResponse::TriMax(Float_t cath[]){
	this->SetTempCaths(cath);	
	for (int i = 0; i<3;i++){	
	triind[i]=MaxInd(tempcath); //This function exports the index and stores the maximum internally
	trimax[i]=this->max;	//this is where the maximum is located in the class. 
	tempcath[triind[i]]=0;
	}
	return;
}

void WireChamberResponse::QuadMax(Float_t cath[]){
	this->SetTempCaths(cath);
	for (int i = 0; i<4;i++){	
	quadind[i]=MaxInd(tempcath); //This function exports the index and stores the maximum internally
	quadmax[i]=max;	//this is where the maximum is located in the class. 
	tempcath[quadind[i]]=0;
	}
	return;
	
}





 void WireChamberResponse::SetCaths(int cathdex){
  

	phys->GetEntry(cathdex);

	return;
}

void WireChamberResponse::SetTempCaths(Float_t cath[]){
   
	for(int i = 0; i<16; i++){
	tempcath[i]=cath[i];	
	}
	return;
}




//This is were the rubber meets the road. 
int WireChamberResponse::ResponseType(Float_t cath[]){
	this->QuadMax(cath);
	if (quadmax[0]>threshold){	
	//Platue Response. 	
		int platnum=1;
		for (int i =0; i<3; i++){	
		if ((0.9*quadmax[0])<quadmax[i+1]) {platnum++;}
		}
	
		if (platnum>1) {return (2+platnum);}
	
	//Triangle Response Type
		else{
		//CLASS 0  single point. 
			if(quadmax[1]<threshold2){return 0; }	
			else{
				if( (quadind[0]-quadind[1])^2==1){
			   	     if((quadmax[1]>1.5*quadmax[2] ||  (quadind[0]-quadind[2])^2>1 ) && (quadind[0]-quadind[1])==1 ){return 2; }
			   	    else if((quadmax[1]>1.5*quadmax[2] ||  (quadind[0]-quadind[2])^2>1 ) && (quadind[0]-quadind[1])==-1 ){return 3; }
			            else if ((quadind[0]-quadind[2])^2==1){return 1;}	   
			         }
		                 else{return 0;}			
			    }	
		    }
	
	}	
	else {return 7;}	 	
}

/*

//testing grounds
int main()
{
	
  	WireChamberResponse *WCR =new WireChamberResponse();
	
	//optional definitions; default values are shown here. 	
	WCR->threshold=120;  //adc threshold (is it really an event?)	
	WCR->threshold2=82;  //triangle daughter threshold. 		
	WCR->platfrac=0.90;
	WCR->trifrac=1.5;	
	WCR->SetPhysTree("$UCNAOUTPUTDIR/hists/spec_22000.root");

  	
	
	for(int i =10; i<20; i++){	
	
	WCR->phys->GetEntry(i);
	
	cout<<WCR->ResponseType(WCR->cathwx)<<"\n";
	cout<<WCR->quadmax[0]<<" "<<WCR->quadmax[1]<<" "<<WCR->quadmax[2]<<" "<<WCR->quadmax[3]<<"\n";	
	cout<<WCR->quadind[0]<<" "<<WCR->quadind[1]<<" "<<WCR->quadind[2]<<" "<<WCR->quadind[3]<<"\n";	

	}
	

	
  return 0;
}


	

/*  //Make Histograms! if 2D use TH2F.. etc.
TH1F * gainhist = new TH1F("gain","gain hist",30,-1,4);//,30,0,10);
 
for (int i = 0; i< linum; ++i){
 //pmtTree->GetEntry(i);
    if(xpos[i]*xpos[i]+ypos[i]*ypos[i] < 25 && gain0[i]<2 && gain0[i]>0){ //&& LG02error[i]>0.0001){
      
  gainhist->Fill(gain0[i]);
  
      }
     }

gainhist->Draw("colz");
*/



