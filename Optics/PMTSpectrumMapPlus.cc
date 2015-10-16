//g++ -o MapSpec PMTSpectrumMapPlus.cc `root-config --cflags --glibs`
#include "PMTSpectrumMapPlus.h"
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

///This is PMTSPectrumMapPlus, it finds the spectrum of each type and 

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
		char *histnameerr = new char[30];
    	sprintf(histname, "pmtmap%d",i);
		sprintf(histnameerr, "pmterrormap%d",i);			
		this->specmap[i]= new TH1F(histname,"",300,0,1000);
		this->specmaperror[i]= new TH1F(histnameerr,"",300,0,1000);
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
				char *histnameerr = new char[30];
    			sprintf(histname, "pmtmap%d",i);
				sprintf(histnameerr, "pmterrormap%d",i);				
				this->specmap[i]= (TH1F*)SpecMapfile->Get(histname);
				this->specmaperror[i]=(TH1F*)SpecMapfile->Get(histnameerr);
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
		specmaperror[ArrayInd(wcpos[2*(ii>3)],wcpos[1+2*(ii>3)],wcresp[2*(ii>3)],wcresp[1+2*(ii>3)],ii)]->Fill(pmtled[ii],pmtlederror[ii]);
					//if(pmtlederror[ii]>1000)cout<<"found broke ass root bullshit, "<<ArrayInd(wcpos[2*(ii>3)],wcpos[1+2*(ii>3)],wcresp[2*(ii>3)],wcresp[1+2*(ii>3)],ii)<<"\n"; seems some are really broken!!				
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


void PMTSpectrumMap::CreateSpectrumError(){
	Float_t number1[300];
	Float_t number2[300];	
	for(int i =0; i<mapentries;i++){
		for(int ii=0; ii<300; ii++){
			number1[ii]=specmaperror[i]->GetBinContent(ii);
			number2[ii]=specmap[i]->GetBinContent(ii);
			//if(number2[ii]>0) cout<<number1[ii]<<"  "<<number2[ii]<<"\n";
		}
		specmaperror[i]->Reset();
		for(int ii=0; ii<300;ii++){
			if(number2[ii]>0) specmaperror[i]->AddBinContent(ii,number1[ii]/number2[ii]);	
			else specmaperror[i]->AddBinContent(ii,0);
		}
	}
}



//find the best local spectrum for fitting. 
void PMTSpectrumMap::FindSpec0(){
	spec11sum=spec12sum=spec13sum=spec21sum=spec22sum=spec23sum=spec31sum=spec32sum=spec33sum=0;
	
	
	/// loop over total array (human readable) version. 
	for(int xi = 2; xi < (xindsize-2); xi++){ // exclude rows too close to the edge. 
	for(int yi =2; yi<(yindsize-2); yi++){
	for(int txi=1;txi<4;txi++){
	for(int tyi=1; tyi<4;tyi++){
	for(int i=0; i<8;i++){
    	
			float sum=0;	
			for(int ii=0; ii<300; ii++){
				sum+=this->specmap[ArrayInd(xi,yi, txi, tyi,i)]->GetBinContent(ii);
				}		
					
			if(txi==1 && tyi==1 && sum>spec11sum){
				spec11sum=sum; 
				xit[0]=xi;yit[0]=yi;txit[0]=txi;tyit[0]=tyi;ti[0]=i;			
			}
		  	else if(txi==2 && tyi==1 && sum>spec21sum){
				spec21sum=sum;
				xit[1]=xi;yit[1]=yi;txit[1]=txi;tyit[1]=tyi;ti[1]=i;			
			}
			else if(txi==1 && tyi==2 && sum>spec12sum){
				spec12sum=sum;
				xit[2]=xi;yit[2]=yi;txit[2]=txi;tyit[2]=tyi;ti[2]=i;			
			}
			else if(txi==2 && tyi==2 && sum>spec22sum){
				spec22sum=sum;
				xit[3]=xi;yit[3]=yi;txit[3]=txi;tyit[3]=tyi;ti[3]=i;			
			}	
			else if(txi==3 && tyi==1 && sum>spec31sum){
				spec31sum=sum; 
				xit[4]=xi;yit[4]=yi;txit[4]=txi;tyit[4]=tyi;ti[4]=i;			
			}		
			else if(txi==1 && tyi==3 && sum>spec13sum){
				spec13sum=sum;
				xit[5]=xi;yit[5]=yi;txit[5]=txi;tyit[5]=tyi;ti[5]=i;			
			}						
			else if(txi==3 && tyi==2 && sum>spec32sum){
				spec32sum=sum;
				xit[6]=xi;yit[6]=yi;txit[6]=txi;tyit[6]=tyi;ti[6]=i;			
			}			
			else if(txi==2 && tyi==3 && sum>spec23sum){
				spec23sum=sum;
				xit[7]=xi;yit[7]=yi;txit[7]=txi;tyit[7]=tyi;ti[7]=i;			
			}
			else if(txi==3 && tyi==3 && sum>spec33sum){
				spec33sum=sum;
				xit[8]=xi;yit[8]=yi;txit[8]=txi;tyit[8]=tyi;ti[8]=i;			
			}

	}}}}}

		 /*	float sum=0;	
			for(int ii=0; ii<300; ii++){
				spec0sum+=this->specmap[ArrayInd(8,8,2,2,4)]->GetBinContent(ii);
				}	*/

	//  ArrayInd(xit,yit,txit,tyit,ti)]   ArrayInd(8,8,2,2,4)
	//cout<<"Spec0 is found in index "<<xit<<" "<<yit<<" "<<txit<<" "<<tyit<<" "<<ti<<"\n";
	for(int i = 0; i<300; i++){
		spec0x[i]=(double)((double)i)*1000.000/300.000;		
		if(i>3){//ignore pedistal. 
			spec11[i]=(double)((specmap[ArrayInd(xit[0],yit[0],txit[0],tyit[0],ti[0])]->GetBinContent(i))/spec11sum); 			
			spec21[i]=(double)((specmap[ArrayInd(xit[1],yit[1],txit[1],tyit[1],ti[1])]->GetBinContent(i))/spec21sum); 
			spec12[i]=(double)((specmap[ArrayInd(xit[2],yit[2],txit[2],tyit[2],ti[2])]->GetBinContent(i))/spec12sum); 
			spec22[i]=(double)((specmap[ArrayInd(xit[3],yit[3],txit[3],tyit[3],ti[3])]->GetBinContent(i))/spec22sum); 
			spec31[i]=(double)((specmap[ArrayInd(xit[4],yit[4],txit[4],tyit[4],ti[4])]->GetBinContent(i))/spec31sum); 
			spec13[i]=(double)((specmap[ArrayInd(xit[5],yit[5],txit[5],tyit[5],ti[5])]->GetBinContent(i))/spec13sum); 
			spec32[i]=(double)((specmap[ArrayInd(xit[6],yit[6],txit[6],tyit[6],ti[6])]->GetBinContent(i))/spec32sum); 
			spec23[i]=(double)((specmap[ArrayInd(xit[7],yit[7],txit[7],tyit[7],ti[7])]->GetBinContent(i))/spec23sum); 
			spec33[i]=(double)((specmap[ArrayInd(xit[8],yit[8],txit[8],tyit[8],ti[8])]->GetBinContent(i))/spec33sum); 

			}
		else{
			spec11[i]=spec12[i]=spec13[i]=spec21[i]=spec22[i]=spec23[i]=spec31[i]=spec32[i]=spec33[i]=0;
			}
		}
		//TSpline3 *sp = new TSpline3("Cubic Spline", Point_x, Point_y, n, "b2e2", 0, 0);
		//spec0inter->SetData(spec0x,spec0);
		spec11inter = new TSpline3("Cubic Spline", spec0x, spec11, 300, "b2e2", 0, 0);
		spec12inter = new TSpline3("Cubic Spline", spec0x, spec12, 300, "b2e2", 0, 0);
		spec13inter = new TSpline3("Cubic Spline", spec0x, spec13, 300, "b2e2", 0, 0);
		spec21inter = new TSpline3("Cubic Spline", spec0x, spec21, 300, "b2e2", 0, 0);
		spec22inter = new TSpline3("Cubic Spline", spec0x, spec22, 300, "b2e2", 0, 0);
		spec23inter = new TSpline3("Cubic Spline", spec0x, spec23, 300, "b2e2", 0, 0);
		spec31inter = new TSpline3("Cubic Spline", spec0x, spec31, 300, "b2e2", 0, 0);
		spec32inter = new TSpline3("Cubic Spline", spec0x, spec32, 300, "b2e2", 0, 0);
		spec33inter = new TSpline3("Cubic Spline", spec0x, spec33, 300, "b2e2", 0, 0);
	return;
}

void PMTSpectrumMap::FillPMTMap(){
	TSpline3 *spec0inter;	
	stringstream fname;
	string filename;
	Float_t alphatemp;
	Double_t tempspec[300];	
	fname<<this->ORD<<"/PMT_SPEC_MAP.root";
	filename=fname.str();	
	this->SpecMapfile = new TFile(filename.c_str(),"recreate");
	this->SpecMap= new TTree("SpecMap", "PMT Spectrum Maps");  
	//TBranch *SpectrumMap = LEDphys->Branch("SpectrumMap",this->specmap,"Runnum/TH1F[165888]");
	TBranch *PmtMap = SpecMap->Branch("PmtMap", &this->PMTMAP,"PmtMap[165888]/F");
	TBranch *PmtError = SpecMap->Branch("PmtError", &this->PMTERRORMAP,"PmtError[165888]/F");
	
	
	//loop through human readable array index 
	for(int xi = 0; xi < xindsize; xi++){
	for(int yi =0; yi<yindsize; yi++){ 
	for(int txi=1;txi<4;txi++){  
	for(int tyi=1; tyi<4;tyi++){
		if(txi==1 && tyi==1) spec0inter=spec11inter;
		if(txi==1 && tyi==2) spec0inter=spec12inter;
		if(txi==1 && tyi==3) spec0inter=spec13inter;
		if(txi==2 && tyi==1) spec0inter=spec21inter;
		if(txi==2 && tyi==2) spec0inter=spec22inter;
		if(txi==2 && tyi==3) spec0inter=spec23inter;
		if(txi==3 && tyi==1) spec0inter=spec31inter;
		if(txi==3 && tyi==2) spec0inter=spec32inter;
		if(txi==3 && tyi==3) spec0inter=spec33inter;
	for(int i=0; i<8;i++){
    
	 alphatemp=MinAlpha(ArrayInd(xi,yi,txi,tyi,i),3,spec0inter);
	PMTMAP[ArrayInd(xi,yi,txi,tyi,i)]=alphatemp;
	PMTERRORMAP[ArrayInd(xi,yi,txi,tyi,i)]=pmtfiterror;
	//if(pmtfiterror<100) cout<<"Finished minimization correction is "<<PMTMAP[ArrayInd(xi,yi,txi,tyi,i)]<<" with error "<<pmtfiterror<<" for index "<<xi<<" "<<yi<<" "<<txi<<" "<<tyi<<" "<<i<<"\n";
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




Float_t PMTSpectrumMap::MinAlpha(int flatind,int resolution,TSpline3 *spec0inter){
	pmttempsum=0;
	if(resolution<=0)resolution=1;
	Float_t alpha=0;
	Float_t alphaprime=100;
	Float_t alphatemp;
	Float_t error;
	Float_t temperror;
	
	for (int i =0; i<300; ++i){
		pmttemp[i]=specmap[flatind]->GetBinContent(i);
		pmttemperror[i]=specmaperror[flatind]->GetBinContent(i);
		pmttempsum+=pmttemp[i]; 		
		}
	
	if (pmttempsum>3000){
		error=1E29;
		for(int iii=0; iii<(resolution+1);iii++){
			for(int ii = 0;ii<20;ii++){
			temperror=ChiSquaredBrah(alpha+pow(0.1,iii)*((float)(ii+1-(iii>0)*10)),spec0inter);
				if(temperror<error) {error=temperror; alpha=alpha+pow(0.1,iii)*((float)(ii+1-(iii>0)*10));
				}
			}
		}	
	this->pmtfiterror=error/300;	

	return alpha;
	}
	else{this->pmtfiterror=0; return 1;}
}



Float_t PMTSpectrumMap::ChiSquaredBrah(Float_t alpha,TSpline3 *spec0inter){
	Float_t dpmt;	//derivative 
	Float_t chitemp=0;

	Float_t derror; //funny sounding
	if (alpha>0){
	for(int i = 0; i<300;i++){
		dpmt=(pmttemp[i+1*(i<299)]-pmttemp[i-1*(i>0)])/pmttempsum*(1000.000/300.000*(float)(2-1*(i==299)-1*(i==0)));
		derror=sqrt(pow(dpmt,2)*pow(pmttemperror[i],2));	
	
		if(spec0x[i]/alpha < 1000 && derror>0){
			chitemp+=pow((pmttemp[i]/pmttempsum - 1/alpha*(float)spec0inter->Eval((Double_t)spec0x[i]/alpha))/derror,2);
					
			}
		else if(derror==0) chitemp+=0; //?? do I need to put this?
		else chitemp+=pow(pmttemp[i]/pmttempsum/derror,2);    
	} 


     return (chitemp);  // error again, chisquared is super fucked up right now, no idea why cause erros and such look perfect.  
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

//Loop to fill all entries into the specturm map. 
 for(int runnum=21966;runnum<22003;runnum++){
	fname<<PSM->ORD<<"hists/LEDspec_"<<runnum<<".root";
	PSM->SetPhysTree(runnum);		
	PSM->FillSpectrumMap();
	if(PSM->FilExists(fname.str())) PSM->LEDphysfile->Close();
	}
	PSM->CreateSpectrumError(); //we were perfect before all this error we created!	
	PSM->SpecMapfile->Write();
    
	cout<<"We made it past specmap write!\n";
	PSM->FindSpec0();
	cout<<"We Made it past find Spec0!\n";
	
	//Graphing spec 0. 	
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
	
	
	
	//Fit and Fill Map!
	PSM->FillPMTMap();
	cout<<"We Made it past FIllPMTMap();\n";
	//app.Run();
	//PSM->SpecMapfile->Close();
	

	//plotting maps. 
 	TApplication app("proportional to photon spectrum", &argc, argv);
	
	Float_t x,y,z;
	TGraph2D * pmap = new TGraph2D();
	int ii=0;

	for(int i = 0; i<15; i++){
		for(int j=0; j<15; j++){
		x=(float)(i+1);
		y=(float)(j+1);		
		z=PSM->PMTERRORMAP[PSM->ArrayInd(i,j,2,2,4)];
		//cout<<x<<" "<<y<<" "<<z<<"\n";
		pmap->SetPoint(ii,(double)x,(double)y,(double)z);		
		ii++;
	}}
	pmap->Draw("surf1");
	
	TCanvas * c2=new TCanvas();	
	c2->cd();
	Float_t x2,y2,z2;
	TGraph2D * pmap2 = new TGraph2D();
	ii=0;
	for(int i = 0; i<15; i++){
		for(int j=0; j<15; j++){
		x2=(float)(i+1);
		y2=(float)(j+1);		
		z2=PSM->PMTMAP[PSM->ArrayInd(i,j,2,2,4)];
		//cout<<x<<" "<<y<<" "<<z<<"\n";
		pmap2->SetPoint(ii,(double)x2,(double)y2,(double)z2);		
		ii++;
	}}
	pmap2->Draw("surf1");

app.Run(); 
	PSM->SpecMapfile->Close();
}

/*
XenonRuns 16581-16583, 16983-17078, 17224-17230, 17561-17734,
18081-18090, 18390-18413, 18712-18744, 19589-19606, 19873-19898,
21596-21605, 21966-22003, 22961-22979
*/










