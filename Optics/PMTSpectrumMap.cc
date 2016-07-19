//g++ -o MapSpec PMTSpectrumMap2.cc `root-config --cflags --glibs`
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
#include <gsl/gsl_mode.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_sf_ellint.h> 
#include <TApplication.h>

using namespace std;



void PMTSpectrumMap::SetPhysTree(int runum){	
	stringstream fname;
	this->rnum=runum;
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
	else{cout<<"\nLEDspec_"<<runum<<".root does not exist in:\n"<<ORD<<"/hists/ \n";return;}
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
void PMTSpectrumMap::ConstructSpectrumMap(){
	for(int i = 0; i<this->mapentries; i++){
		char *histname = new char[30];
    	sprintf(histname, "pmtmap%d",i);
		this->specmap[i]= new TH1F(histname,"",300,0,1000);
		if(i%(mapentries/3)==0) cout<<"Creating histograms..."<<i<<"\n";
	}
}

void PMTSpectrumMap::LoadSpectrumMap(string filename){  //slighlty faster than simply remaking?? is it????!?!?!!!
		stringstream fname;
		fname<<this->ORD<<filename;
		if(this->FilExists(fname.str())){
			filename=fname.str();			
			SpecMapfile =  TFile::Open(filename.c_str());
			for(int i = 0; i<this->mapentries; i++){
				char *histname = new char[30];
    			sprintf(histname, "pmtmap%d",i);
				this->specmap[i]= (TH1F*)SpecMapfile->Get(histname);
				if(i%(mapentries/10)==0 && i>0) cout<<"Loading histograms..."<<i<<"/"<<mapentries<<"\n";

			}
		return;		
		}
		else{ cout<<"Spectrum Map File Not Found: "<<fname<<"\n"; return;}


	}

void PMTSpectrumMap::FillSpectrumMap(){
	stringstream fname;
 	fname<<this->ORD<<"hists/LEDspec_"<<rnum<<".root";
	if(this->FilExists(fname.str())){


				cout<<"About to fill "<<this->entnum<<" entries into histograms. Good luck!\n";
				for(int i = 0; i<this->entnum; i++){
				this->LEDphys->GetEntry(i);	
				if(i%(this->entnum/3)==0) cout<<"Filling... completed "<<i<<" Histograms for run "<<rnum<<"\n";			
					
				  for(int ii = 0; ii<8; ii++){
	
				
									
			//flattened spec map. a 16x16x8x8x8 array! maybe it should be a 16x7x16x7x8 array. 
			 specmap[ArrayInd(wcpos[2*(ii>3)],wcpos[1+2*(ii>3)],wcresp[2*(ii>3)],wcresp[1+2*(ii>3)],ii)]->Fill(pmtled[ii]);
						
					}								
				}	
	cout<<"seems like it worked and everything is ok!\n";	
	return;}
	else cout<<"There doesn't appear to be a TTree for run "<<rnum<<"\n";
	return;
}
// un flatten to flatten array
int PMTSpectrumMap::ArrayInd(int xind,int yind, int typex, int typey, int pmt){	

	return (mapentries/xindsize*xind + yind*typeindsize*typeindsize*8 + typex*typeindsize*8 + typey*8 + pmt);
}





//find the best local spectrum for fitting. 
void PMTSpectrumMap::FindSpec0(){
	spec0sum=0;
	
	/// loop over total array (human readable) version. 
	for(int xi = 2; xi < xindsize-2; xi++){// exclude rows too close to the edge. 
	for(int yi =2; yi<yindsize-2; yi++){
	for(int txi=0;txi<typeindsize;txi++){
	for(int tyi=0; tyi<typeindsize;tyi++){
	for(int i=0; i<8;i++){
    	
			float sum=0;	
			for(int ii=0; ii<300; ii++){
				sum+=this->specmap[ArrayInd(xi,yi, txi, tyi,i)]->GetBinContent(ii);
				}				
			if(sum>spec0sum){
				spec0sum=sum;
				xit=xi;yit=yi;txit=txi;tyit=tyi;ti=i;			
			}
		  
						
			
	}}}}}
		 /*	float sum=0;	
			for(int ii=0; ii<300; ii++){
				spec0sum+=this->specmap[ArrayInd(8,8,2,2,4)]->GetBinContent(ii);
				}	*/

	//  ArrayInd(xit,yit,txit,tyit,ti)]   ArrayInd(8,8,2,2,4)
	cout<<"Spec0 is found in index "<<xit<<" "<<yit<<" "<<txit<<" "<<tyit<<" "<<ti<<"\n";
	for(int i = 0; i<300; i++){
		spec0x[i]=(double)((double)i)*1000.000/300.000;		
		if(i>3){//ignore pedistal. 
			spec0[i]=(double)((specmap[ArrayInd(xit,yit,txit,tyit,ti)]->GetBinContent(i))/spec0sum); 
			
			}
		else{
			spec0[i]=0;
			}
		}
		//TSpline3 *sp = new TSpline3("Cubic Spline", Point_x, Point_y, n, "b2e2", 0, 0);
		//spec0inter->SetData(spec0x,spec0);
		spec0inter = new TSpline3("Cubic Spline", spec0x, spec0, 300, "b2e2", 0, 0);
	return;
}

void PMTSpectrumMap::FillPMTMap(){
	stringstream fname;
	string filename;
	Float_t alphatemp;
	Double_t tempspec[300];	
	fname<<this->ORD<<"/PMT_SPEC_MAP.root";
	filename=fname.str();	
	this->SpecMapfile = new TFile(filename.c_str(),"recreate");
	this->SpecMap= new TTree("SpecMap", "PMT Spectrum Maps");  
	//TBranch *SpectrumMap = LEDphys->Branch("SpectrumMap",this->specmap,"Runnum/TH1F[165888]");
	TBranch *PmtMap = SpecMap->Branch("PmtMap", &this->PMTMAP,"PmtMap/F[165888]");
	TBranch *PmtError = SpecMap->Branch("PmtError", &this->PMTERRORMAP,"PmtError/F[165888]");
	
	
	//loop through human readable array index 
	for(int xi = 0; xi < xindsize; xi++){
	for(int yi =0; yi<yindsize; yi++){ 
	for(int txi=0;txi<typeindsize;txi++){  
	for(int tyi=0; tyi<typeindsize;tyi++){
	for(int i=0; i<8;i++){
	 alphatemp=MinAlpha(ArrayInd(xi,yi,txi,tyi,i),3);
	PMTMAP[ArrayInd(xi,yi,txi,tyi,i)]=alphatemp;
	PMTERRORMAP[ArrayInd(xi,yi,txi,tyi,i)]=pmtfiterror;
	//if(pmtfiterror>0) cout<<"Finished minimization correction is "<<PMTMAP[ArrayInd(xi,yi,txi,tyi,i)]<<" with error "<<pmtfiterror<<" for index "<<xi<<" "<<yi<<" "<<txi<<" "<<tyi<<" "<<i<<"\n";
	/*if (txi==2 && tyi==2){	
		if(pmttempsum>3000){
		for(int ii = 0; ii<300;ii++){pmttemp[ii]=pmttemp[ii]/pmttempsum*alphatemp;tempspec[ii]=(Double_t)spec0x[ii]*alphatemp;}
		tempgraph=new TGraph(300,tempspec,(Double_t*)pmttemp);
		tempgraph->Draw("same");}
		}*/
	}}}}}
	SpecMap->Fill();
	SpecMap->Write();
	//specmap->Write();
	SpecMapfile->Write();
	SpecMapfile->Close();	
	return;
}




Float_t PMTSpectrumMap::MinAlpha(int flatind,int resolution){
	pmttempsum=0;
	if(resolution<=0)resolution=1;
	Float_t alpha=0;
	Float_t alphaprime=100;
	Float_t alphatemp;
	Float_t error;
	Float_t temperror;
	
	for (int i =0; i<300; ++i){
		pmttemp[i]=specmap[flatind]->GetBinContent(i);
		pmttempsum+=pmttemp[i]; 		
		}
	
	if (pmttempsum>3000){
		error=1;
		for(int iii=0; iii<(resolution+1);iii++){
			for(int ii = 0;ii<20;ii++){
			temperror=ChiSquaredBrah(alpha+pow(0.1,iii)*((float)(ii+1-(iii>0)*10)));
				if(temperror<error) {error=temperror; alpha=alpha+pow(0.1,iii)*((float)(ii+1-(iii>0)*10));
				}
			}
		}	
	this->pmtfiterror=error;	
	return alpha;	
	}
	else{this->pmtfiterror=0; return 0;}
}


Float_t PMTSpectrumMap::ChiSquaredBrah(Float_t alpha){//TODO: Propagate LED Error to make real Chi^2!
	
	Float_t chitemp=0;
	if (alpha>0){
	for(int i = 0; i<300;i++){
		if(spec0x[i]/alpha<1000) chitemp+=pow(pmttemp[i]/pmttempsum-1/alpha*(float)spec0inter->Eval((Double_t)spec0x[i]/alpha),2);
		else chitemp+=pow(pmttemp[i]/pmttempsum,2);
		
	} 

	
     return sqrt(chitemp);  //THIS IS FRACTIONAL ERROR,  not a chi^2 because LED error is not implemented.  
	}
	else return 1e30; //big number... hopefully. 
}



//examples!!!!!

int main(int argc, char **argv){
	stringstream fname;	
	string filename;			
	PMTSpectrumMap * PSM = new PMTSpectrumMap();
	fname<<PSM->ORD<<"PMT_SPEC_MAP.root";
	filename=fname.str();
	PSM->SpecMapfile = new TFile(filename.c_str(),"recreate");
	PSM->ConstructSpectrumMap();
	
	//If you would rather load the spectrum map Do this!!!, assuming the file hasn't broke! Which it just did for me! DAMN!
	//PSM->SpecMapfile =  TFile::Open(filename.c_str());SpecMapFile
	//PSM->LoadSpectrumMap(filename);

 for(int runnum=21966;runnum<22003;runnum++){
	fname<<PSM->ORD<<"hists/LEDspec_"<<runnum<<".root";
	PSM->SetPhysTree(runnum);		
	PSM->FillSpectrumMap();
	if(PSM->FilExists(fname.str())) PSM->LEDphysfile->Close();
	}
	
	PSM->SpecMapfile->Write();
    

	PSM->FindSpec0();
	/* Double_t x[300],y[300];

	for(int i =0; i<300; i++){
	//PSM->spec0inter->Eval((Double_t)spec0x[i],2);
	x[i]=(Double_t)PSM->spec0x[i]/1.3;
	y[i]=(Double_t)PSM->spec0inter->Eval((Double_t)PSM->spec0x[i]);
	}
	TGraph *spec0spline= new TGraph(300,x,y);
	TGraph *spec0graph = new TGraph(300,PSM->spec0x,PSM->spec0);
	*/
	//TCanvas *c1d= new TCanvas();   	
	
	//PSM->specmap[PSM->ArrayInd(8,8,2,2,4)]->Draw("colz");
   // TCanvas *c2d= new TCanvas("c2d","PMT Map",600,600);
	//c2d->cd();
	//PSM->specmap[PSM->ArrayInd(8,8,2,2,4)]->Draw("colz");
	//spec0spline->Draw();
	//spec0graph->SetLineColor(2);
	//spec0graph->Draw("same");	
	
	//PSM->FindSpec0();
	PSM->FillPMTMap();
	//app.Run();
	//PSM->SpecMapfile->Close();
	

 	TApplication app("proportional to photon spectrum", &argc, argv);
	

	
	Float_t x,y,z;
	TGraph2D * pmap = new TGraph2D();
	int ii=0;

	for(int i = 0; i<15; i++){
		for(int j=0; j<15; j++){
		x=(float)(i+1);
		y=(float)(j+1);		
		z=PSM->PMTERRORMAP[PSM->ArrayInd(i,j,3,3,4)];
		//cout<<x<<" "<<y<<" "<<z<<"\n";
		pmap->SetPoint(ii,(double)x,(double)y,(double)z);		
		ii++;
	}}
	pmap->Draw("colz");
	
	TCanvas * c2=new TCanvas();	
	c2->cd();
	Float_t x2,y2,z2;
	TGraph2D * pmap2 = new TGraph2D();
	ii=0;
	for(int i = 0; i<15; i++){
		for(int j=0; j<15; j++){
		x2=(float)(i+1);
		y2=(float)(j+1);		
		z2=PSM->PMTMAP[PSM->ArrayInd(i,j,3,3,4)];
		//cout<<x<<" "<<y<<" "<<z<<"\n";
		pmap2->SetPoint(ii,(double)x2,(double)y2,(double)z2);		
		ii++;
	}}
	pmap2->Draw("colz");

app.Run(); 
	PSM->SpecMapfile->Close();
}

/*
XenonRuns 16581-16583, 16983-17078, 17224-17230, 17561-17734,
18081-18090, 18390-18413, 18712-18744, 19589-19606, 19873-19898,
21596-21605, 21966-22003, 22961-22979
*/











