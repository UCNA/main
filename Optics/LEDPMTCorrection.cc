//g++ -o LEDPMT LEDPMTCorrection.cc WireChamberResponse.cpp `root-config --cflags --glibs`
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
	

     //phys tree is Michael M's output, I don't use analyzed data (i think, atleast its not too analyzed) just using basic data. 
void LEDPMTApp::SetPhysTree(int rnum){
	stringstream fname;
 	fname<<this->ORD<<"hists/spec_"<<rnum<<".root";
	if(this->FilExists(fname.str())){
	string filename=fname.str();
  	this->physfile = TFile::Open(filename.c_str());
  	phys = (TTree*)physfile->Get("phys");
  	entnum=phys->GetEntries();
  	phys->SetBranchAddress("ScintE",&scintE);
  	phys->SetBranchAddress("ScintW",&scintW);					
  return;}
	else{cout<<"\nspec_"<<rnum<<".root does not exist in:\n"<<ORD<<"/hists/ \n";return;}
}


//this function is used if the LEData root file exists already
void LEDPMTApp::LoadLEData(string filename){
	if(!import){  	
	this->LEDfitfile = TFile::Open(filename.c_str());
  	LEDfit = (TTree*)LEDfitfile->Get("LEDfit");
	}
	LEDfit->SetBranchAddress("Runnum",&runnum);  	
	LEDfit->SetBranchAddress("LED0",&led0);
  	LEDfit->SetBranchAddress("LEDerror0",&lederror0);
	LEDfit->SetBranchAddress("LED1",&led1);
  	LEDfit->SetBranchAddress("LEDerror1",&lederror1);
	LEDfit->SetBranchAddress("LED2",&led2);
  	LEDfit->SetBranchAddress("LEDerror2",&lederror2);						
  return;
}

int LEDPMTApp::ImportLED(string filename){
	SetMatrixFile(filename);
	string fname;
	string dot=".";
	std::size_t dotpos;	
	dotpos=filename.find(dot);
	fname.append(ORD);
	fname.append("/");
	fname.append(filename.substr(0,dotpos));
	fname.append(".root");	
	
	this->LEDfitfile= new TFile(fname.c_str(),"recreate");
	this->LEDfit= new TTree("LEDfit", "Matrix of LED fit data");  
	TBranch *Runnum = LEDfit->Branch("Runnum",&runnum,"Runnum/I");	
	TBranch *LED0 = LEDfit->Branch("LED0",&led0,"led0[8]/F");
	TBranch *LEDerror0 = LEDfit->Branch("LEDerror0",&lederror0,"lederror0[8]/F");
	TBranch *LED1 = LEDfit->Branch("LED1",&led1,"led1[8]/F");
	TBranch *LEDerror1=LEDfit->Branch("LEDerror1",&lederror1,"lederror1[8]/F");
	TBranch *LED2= LEDfit->Branch("LED2",&led2,"led2[8]/F"); 	
	TBranch *LEDerror2=LEDfit->Branch("LEDerror2",&lederror2,"lederror2[8]/F");
	string ledin;
  	std::string delim (" ");
	
	
	if (this->FilExists(filename.c_str())){  	
	import=true;	
	ifstream myfile(filename.c_str());
  	string stupid="stupid";
  	
    		while ( getline (myfile,ledin) )
    			{	

			
			   
			//cycle throughthe first few rounds (of useless data (for our purposes (at the moment)))^20
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
						}
			
			this->LEDfit->Fill();
		
    			}
			LoadLEData(" "); //file name not needed after import.  but a string is. 
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
	string matrixfile;
	this->matrixfile.append(ORD);
	this->matrixfile.append(filename.c_str());
}


void LEDPMTApp::ApplyCorrection(int rnum){
	
	stringstream fname,fname2;
	fname2<<this->ORD<<"hists/spec_"<<rnum<<".root";
	if(this->FilExists(fname2.str())){
	Float_t wctemp[16]; Int_t wcindtemp[16];Float_t wcresptemp[4]; Int_t wcpostemp[4];		
 	fname<<this->ORD<<"hists/LEDspec_"<<rnum<<".root";
	string filename=fname.str();	
		
	

	isthere=false;	
	for(int iii =0; iii<this->LEDfit->GetEntries(); iii++){	 	
			this->LEDfit->GetEntry(iii);
		
	 		if(this->runnum==rnum){
					this->LEDphysfile = new TFile(filename.c_str(),"recreate");
					this->LEDphys= new TTree("LEDphys", "LED fit Corrected PMT to photon proportional value");  
					TBranch *Runnum = LEDphys->Branch("Runnum",&rnum,"Runnum/I");	
					TBranch *PMTLED = LEDphys->Branch("PMTLED",&this->PMTLED,"PMTLED[8]/F");
					TBranch *PMTLEDError = LEDphys->Branch("PMTLEDError",&this->PMTLEDError,"PMTLEDError[8]/F");
					TBranch *PMT = LEDphys->Branch("PMT",&this->PMT,"PMT[8]/F");
					TBranch *WCResponse = LEDphys->Branch("WCResponse",&wcresptemp,"WCResponse[4]/F");
					TBranch *WCPosition = LEDphys->Branch("WCPosition",&wcpostemp,"WCPosition[4]/I");	
					TBranch *WCMaxind = LEDphys->Branch("WCMaxind",&wcindtemp,"WCMaxind[16]/I");	
					TBranch *WCMax = LEDphys->Branch("WCMax",&wctemp,"WCMax[16]/F");
					isthere=true;			
			cout<<"\nDanger: Filling run "<<rnum<<"'s TTree with LED corrected PMT values.\nKeep all arms and legs inside the cart.\n";
				for(int ii = 0; ii< this->phys->GetEntries();ii++){
				   if(ii%10000==0)cout<<"*";
				   this->WCR->phys->GetEntry(ii);
				   this->phys->GetEntry(ii);
				   for(int i =0; i<4;i++){
					//collecting and setting all the values.
					this->PMT[i]=scintE[i];
					this->PMT[i+4]=scintW[i];
					this->PMTLED[i]=this->led0[i]+this->led1[i]*scintE[i]+this->led2[i]*pow(scintE[i],2);
					this->PMTLED[i+4]=this->led0[i+4]+this->led1[i+4]*scintW[i]+this->led2[i+4]*pow(scintW[i],2);
					this->PMTLEDError[i]=this->lederror0[i]+this->lederror1[i]*scintE[i]+this->lederror2[i]*pow(scintE[i],2);
					this->PMTLEDError[i+4]=this->lederror0[i+4]+this->lederror1[i+4]*scintW[i]+this->lederror2[i+4]*pow(scintW[i],2);
					wcresptemp[0]=this->WCR->ResponseType(this->WCR->cathex);
					wcpostemp[0]=this->WCR->wcpos;					
					wctemp[i]=this->WCR->quadmax[i];wcindtemp[i]=this->WCR->quadind[i];					
					wcresptemp[1]=this->WCR->ResponseType(this->WCR->cathey);
					wcpostemp[1]=this->WCR->wcpos;					
					wctemp[i+4]=this->WCR->quadmax[i];wcindtemp[i+4]=this->WCR->quadind[i];	
				    wcresptemp[2]=this->WCR->ResponseType(this->WCR->cathwx);
					wcpostemp[2]=this->WCR->wcpos;					
					wctemp[i+8]=this->WCR->quadmax[i];wcindtemp[i+8]=this->WCR->quadind[i];	
					wcresptemp[3]=this->WCR->ResponseType(this->WCR->cathwy);
					wcpostemp[3]=this->WCR->wcpos;					
					wctemp[i+12]=this->WCR->quadmax[i];wcindtemp[i+12]=this->WCR->quadind[i];
											
				}
				
				this->LEDphys->Fill();
				
			}
		}
	}
	if(!isthere) cout<<"\nThere is no LED fit data for run "<<rnum<<". Sorry (well, I'm a computer, so I'm not sorry.)\n";
	cout<<"\n"; return;
	}
	else {cout<<"\n"<<"File "<<fname2.str()<<" does not exist.\nSorry (not really.)\n";return;}
}


///EXAMPLES!!!
int main(){
	stringstream fname;

	LEDPMTApp *LED =new LEDPMTApp();
	
	//file must must be in $UCNAOUPPUTDIR). 
	LED->ImportLED("MatrixLEData.txt");	
	LED->LEDfit->Write();	
	for(int runnum=20900;runnum<23713;runnum++){ 
	
		
		fname<<LED->ORD<<"hists/spec_"<<runnum<<".root";
	
		//int runnum=21966;
		LED->SetPhysTree(runnum);
		//Default values shown. 	
		LED->WCR->threshold=100;  //adc threshold (is it really an event?)	
		LED->WCR->threshold2=62;  //triangle daughter threshold. 		
		LED->WCR->platfrac=0.90; //platue criteria 90% of highest value is partof the platue
		LED->WCR->trifrac=1.5;  // a triangle hasto be worth 1.5 times the lower triangle to be a left or right leaning 		triangle.
		LED->WCR->SetPhysTree(runnum);
		/*for(int i= 1; i<100; i++){	
		LED->WCR->phys->GetEntry(i);
		cout<<"Cathode Type: "<<LED->WCR->ResponseType(LED->WCR->cathwx)<<"\n";
		cout<<"At wire position: "<<LED->WCR->wcpos<<"\n";
		}*/
	
		LED->ApplyCorrection(runnum); 
	

	//there is no controll to not seg fault on the TTree subclasses if the file does not exist. only the LED and WCR classes. 
		if(LED->FilExists(fname.str())){	
			LED->physfile->Close();
			LED->WCR->physfile->Close();
			if(LED->isthere){			
			LED->LEDphys->Write();	
			LED->LEDphysfile->Close();	
			}
	} 	
	//always close the files at the end of a use
	

	}
	//same LEDfits use same TFile for all run numbers, close after loop. 
	LED->LEDfitfile->Close();
	
	///Create Spectrum Map
	
	PMTSpectrumMap * PSM = new PMTSpectrumMap();
	PSM->ConstructSpectrumMap();
 
	for(int runnum=21966;runnum<22003;runnum++){
		fname<<LED->ORD<<"hists/LEDspec_"<<runnum<<".root";
		PSM->SetPhysTree(runnum);		
		PSM->FillSpectrumMap();
		if(PSM->FilExists(fname.str())) PSM->LEDphysfile->Close();
	}

  


}

