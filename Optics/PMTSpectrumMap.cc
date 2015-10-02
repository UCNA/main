//g++ -o MapSpec PMTSpectrumMap.cc `root-config --cflags --glibs`
#include "PMTSpectrumMap.h"
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



void PMTSpectrumMap::SetPhysTree(int rnum){
	stringstream fname;
	this->rnum=rnum;
 	fname<<this->ORD<<"hists/LEDspec_"<<rnum<<".root";
	if(this->FilExists(fname.str())){
	string filename=fname.str();
  	this->LEDphysfile = TFile::Open(filename.c_str());
  	LEDphys = (TTree*)LEDphysfile->Get("LEDphys");
  	entnum=LEDphys->GetEntries();
	LEDphys->SetBranchAddress("Runnum",&this->runnum);	
	LEDphys->SetBranchAddress("PMTLED",&this->pmtled);
	LEDphys->SetBranchAddress("PMTLEDError",&this->pmtlederror);
	LEDphys->SetBranchAddress("PMT",&this->pmt);
	LEDphys->SetBranchAddress("WCResponse",&this->wcresp);
	LEDphys->SetBranchAddress("WCPosition",&this->wcpos);
	LEDphys->SetBranchAddress("WCMaxind",&this->wcind);
	LEDphys->SetBranchAddress("WCMax",&this->wcmax);


					
  return;}
	else{cout<<"\nspec_"<<rnum<<".root does not exist in:\n"<<ORD<<"/hists/ \n";return;}
}


 bool PMTSpectrumMap::FilExists(const std::string& name) {
    ifstream f(name.c_str());
    if (f.good()) {
        f.close();
        return true;
    } else {
        f.close();
        return false;
    }   
}

void PMTSpectrumMap::FillSpectrumMap(){
		

	for(int i = 0; i<this->mapentries; i++){
		char *histname = new char[30];
    	sprintf(histname, "%d",i);
		this->specmap[i]= new TH1F(histname,"",300,0,1000);
		if(i%10000==0) cout<<"Creating histograms..."<<i<<"\n";
	}
				cout<<"About to fill "<<this->entnum<<" entries into histograms. Good luck!\n";
				for(int i = 0; i<this->entnum; i++){
				this->LEDphys->GetEntry(i);	
				if(i%this->entnum/10==0) cout<<"Filling Histograms.."<<i<<"\n";			
				  for(int ii = 0; ii<8; ii++){
				if(ArrayInd(wcpos[2*(ii>3)],wcpos[1+2*(ii>3)],wcresp[2*(ii>3)],wcresp[1+2*(ii>3)],ii)<mapentries){
					//flattened spec map. a 16x16x7x7x8 array! maybe it should be a 16x7x16x7x8 array. 
			 specmap[ArrayInd(wcpos[2*(ii>3)],wcpos[1+2*(ii>3)],wcresp[2*(ii>3)],wcresp[1+2*(ii>3)],ii)]->Fill(pmtled[ii]);
						}
					}								
				}	
	cout<<"seems like everything is ok!";	
	return;
}
//flatten to un flatten array
int PMTSpectrumMap::ArrayInd(int xind,int yind, int typex, int typey, int pmt){	
	int aa=mapentries/8*pmt;
	int bb=xindsize*yindsize*typeindsize*typey;
	int cc=xindsize*yindsize*typex;
	int dd=xindsize*yind;		
	return (aa+bb+cc+dd+xind);
}

int main(){
	int runnum=21968;	
	PMTSpectrumMap * PSM = new PMTSpectrumMap();
	PSM->SetPhysTree(runnum);	
	PSM->FillSpectrumMap();
	cout<<"hello!"<<"\n";
	PSM->specmap[PSM->ArrayInd(8,8,0,0,4)]->Draw("colz");

}










