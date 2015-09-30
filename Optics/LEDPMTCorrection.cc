//g++ -o LEDPMT LEDPMTCorrection.cc WireChamberResponse.cpp `root-config --cflags --glibs`
#include "WireChamberResponse.h"
#include "LEDPMTCorrection.h"
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

	//this class is designed to take a maxtrix with LED quadratic fits (in simon's format but turned into a matrix) and apply them to the data;
	//it will have a feature to export a new Ttree with the new data. 		


     //phys tree is michaels output, I don't use analyzed data (i think) just using basic data. 
void LEDPMTApp::SetPhysTree(string filename){

  	TFile* myFile2 = TFile::Open(filename.c_str());
  	phys = (TTree*)myFile2->Get("phys");
  	entnum=phys->GetEntries();
  	phys->SetBranchAddress("ScintE",&scintE);
  	phys->SetBranchAddress("ScintW",&scintW);					
  return;
}


//this function is used if the LEData root file exists already
void LEDPMTApp::LoadLEData(string filename){
  	this->myFile = TFile::Open(filename.c_str());
  	LEDfit = (TTree*)myFile->Get("LEDfit");
  	LEDfit->SetBranchAddress("LED0",&led0);
  	LEDfit->SetBranchAddress("LEDerror0",&lederror0);
	LEDfit->SetBranchAddress("LED1",&led1);
  	LEDfit->SetBranchAddress("LEDerror1",&lederror1);
	LEDfit->SetBranchAddress("LED2",&led2);
  	LEDfit->SetBranchAddress("lederror2",&lederror2);						
  return;
}

int LEDPMTApp::ImportLED(string filename){
	SetMatrixFile(filename);
	string fname;
	string dot=".";
	std::size_t dotpos;	
	dotpos=filename.find(dot);
	fname=filename.substr(0,dotpos);
	fname.append(".root");	
	this->myFile= new TFile(fname.c_str(),"recreate");
	this->LEDfit= new TTree("LEDfit", "Matrix of LED fit data");  
	TBranch *LED0 = LEDfit->Branch("LED0",&led0,"led0[8]/F");
	TBranch *LEDerror0 = LEDfit->Branch("LEDerror0",&lederror0,"lederror0[8]/F");
	TBranch *LED1 = LEDfit->Branch("LED1",&led1,"led1[8]/F");
	TBranch *LEDerror1=LEDfit->Branch("LEDerror1",&lederror1,"lederror1[8]/F");
	TBranch *LED2= LEDfit->Branch("LED2",&led2,"led2[8]/F"); 	
	TBranch *LEDerror2=LEDfit->Branch("LEDerror2",&lederror2,"lederror2[8]/F");
	string ledin;
  	std::string delim (" ");
	
	
	if (this->FilExists(filename.c_str())){  	
	ifstream myfile(filename.c_str());
  	string stupid="stupid";
  	
    		while ( getline (myfile,ledin) )
    			{	

			
			   
			//cycle throughthe first few rounds (of useless data (for our purposes (at the moment)))
			ledin.append(stupid);
			std::size_t found = ledin.find(delim);
			std::size_t dumb  = ledin.find(stupid);
      			std::string local = ledin.substr (0,found);
			std::string ledintemp = ledin.substr(found+1,dumb);
			runnum=atoi(local.c_str());
			for(int i = 0;i<4;i++){
						ledintemp.append(stupid);
						found=ledintemp.find(delim);
						dumb=ledintemp.find(stupid);
						local = ledintemp.substr (0,found);
						ledintemp=ledintemp.substr(found+1,dumb);					
						}
      			
			///starting this loop we are on the 6th entry. This is the start of the LED fit data. 
			for(int i = 0; i<64; i++){
      						ledintemp.append(stupid);
						found=ledintemp.find(delim);
						dumb=ledintemp.find(stupid);
						local = ledintemp.substr (0,found);
						ledintemp=ledintemp.substr(found+1,dumb);
      						if(i%8==0) this->led0[i/8]=atof(local.c_str());
						if(i%8==1) this->lederror0[i/8]=atof(local.c_str());
						if(i%8==2) this->led1[i/8]=atof(local.c_str());
						if(i%8==3) this->lederror1[i/8]=atof(local.c_str());
						if(i%8==4) this->led2[i/8]=atof(local.c_str());
						if(i%8==5) this->lederror2[i/8]=atof(local.c_str());
						cout<<i/8<<"\n";
						}
			Float_t rootsucks=led1[0];
			cout<<rootsucks<<"\n";
			this->LEDfit->Fill();
			//Float_t rootsucks=led1[0];
			//cout<<rootsucks<<"\n";
    			}
    			myfile.close();
			return 1;
  	}

  	else {cout <<"\n"<<"Unable to open file"<<"\n"<<"\n"; return 0;}

}


  bool LEDPMTApp::FilExists(const std::string& name) {
    ifstream f(name.c_str());
    if (f.good()) {
        f.close();
        return true;
    } else {
        f.close();
        return false;
    }   
}

void LEDPMTApp::SetMatrixFile(string filename){
	char* ORD;
  	ORD = getenv ("UCNAOUTPUTDIR");
	string matrixfile;
	this->matrixfile.append(ORD);
	this->matrixfile.append("/");
	this->matrixfile.append(filename.c_str());
}





int main(){
	
	//Initialize a wire chamber response!!!...
	WireChamberResponse *WCR =new WireChamberResponse();
	//optional values (default shown)	
	WCR->threshold=120;  //adc threshold (is it really an event?)	
	WCR->threshold2=82;  //triangle daughter threshold. 		
	WCR->platfrac=0.90; //platue criteria 90% of highest value is partof the platue
	WCR->trifrac=1.5;  // a triangle hasto be worth 1.5 times the lower triangle to be a left or right leaning triangle.
	//load file. 	
	WCR->SetPhysTree("$UCNAOUTPUTDIR/hists/spec_22000.root");	
	
	LEDPMTApp *LED =new LEDPMTApp();
	LED->SetPhysTree("$UCNAOUTPUTDIR/hists/spec_22000.root");
	

	//ONLY THE FILE NAME. must be in official replay data). 
	LED->ImportLED("MatrixLEData.txt");

	//LED->phys->Show(13); 
	//LED->phys->Show(12);//WORKS!!!
	LED->LEDfit->Write();  //WORKS, numbers are in the right spots and are the right values also. 
	
}

