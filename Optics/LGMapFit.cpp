//g++ -o PMTfit LGMapFit.cpp `root-config --cflags --glibs`
#include <string>
#include <ctype.h>
#include <algorithm>
#include <functional>

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
//#include "pmtprobstuff.h"
//#include "lgpmtTools.h"
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
#include "TH3F.h"
#include "TPaveText.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TError.h"

#include <TString.h>
#include <TCut.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TApplication.h>


using namespace std;



vector<double> PMTadc;vector<double> gain0;vector<double> xpos;vector<double> ypos;
//a dozen LG probability vectors.        
vector<double> LG01prob;vector<double> LG02prob;vector<double> LG03prob;vector<double> LG04prob;vector<double> LG05prob;vector<double> LG06prob;vector<double> LG07prob;vector<double> LG08prob;vector<double> LG09prob;vector<double> LG10prob;vector<double> LG11prob;vector<double> LG12prob;  
//a dozen LG error vectors 
vector<double> LG01error;vector<double> LG02error;vector<double> LG03error;vector<double> LG04error;vector<double> LG05error;
vector<double> LG06error;vector<double> LG07error;vector<double> LG08error;vector<double> LG09error;vector<double> LG10error;vector<double> LG11error;vector<double> LG12error;




void SetLGMap (TTree * pmtTree, TTree * flatpmtTree,int east,int pmt) {

  Double_t lgEprob[12];
  Double_t lgWprob[12];
  Double_t lgEerr[12];
  Double_t lgWerr[12];
  Float_t scintE[4];  ///
  Float_t scintW[4];
  Float_t Wpos[2];
 pmtTree->SetBranchAddress("LGEprob",&lgEprob);
 pmtTree->SetBranchAddress("LGWprob",&lgWprob);
 pmtTree->SetBranchAddress("LGEerror",&lgEerr);
 pmtTree->SetBranchAddress("LGWerror",&lgWerr);
 //flatpmtTree->SetBranchAddress("ScintEflat",&scintE);
 //flatpmtTree->SetBranchAddress("ScintWflat",&scintW);
 pmtTree->SetBranchAddress("ScintW",&scintW);
 pmtTree->SetBranchAddress("Wpos",&Wpos);
 
 int linum=pmtTree->GetEntries();
   
//define vector lengths
     PMTadc.reserve(linum);gain0.reserve(linum);LG01prob.reserve(linum); LG02prob.reserve(linum); LG03prob.reserve(linum); LG04prob.reserve(linum); LG05prob.reserve(linum); LG06prob.reserve(linum); LG07prob.reserve(linum); LG08prob.reserve(linum); LG09prob.reserve(linum); LG10prob.reserve(linum); LG11prob.reserve(linum); LG12prob.reserve(linum);LG01error.reserve(linum); LG02error.reserve(linum); LG03error.reserve(linum); LG04error.reserve(linum); LG05error.reserve(linum); LG06error.reserve(linum); LG07error.reserve(linum); LG08error.reserve(linum); LG09error.reserve(linum); LG10error.reserve(linum); LG11error.reserve(linum); LG12error.reserve(linum);
xpos.reserve(linum); ypos.reserve(linum); 

if (east == 0) {   cout<<linum<<"\n";
   for( int i =0; i<linum; i++){
        pmtTree->GetEntry(i);
        flatpmtTree->GetEntry(i);
       PMTadc.push_back(scintW[pmt]);
       xpos.push_back(Wpos[0]);
       ypos.push_back(Wpos[1]);
       gain0.push_back(4*scintW[0]/(scintW[0]+scintW[1]+scintW[2]+scintW[3]));
         
       LG01prob.push_back(lgWprob[0]);     //why no std::maxtrix<double>? why c++ WHY?!?!!
       LG02prob.push_back(lgWprob[1]);	     //preprocessor loop?
       LG03prob.push_back(lgWprob[2]);	
       LG04prob.push_back(lgWprob[3]);	
       LG05prob.push_back(lgWprob[4]);	
       LG06prob.push_back(lgWprob[5]);	
       LG07prob.push_back(lgWprob[7]);	
       LG08prob.push_back(lgWprob[8]);	
       LG09prob.push_back(lgWprob[9]);	
       LG10prob.push_back(lgWprob[10]);	
       LG11prob.push_back(lgWprob[11]);	
       LG01error.push_back(lgWerr[0]);	
       LG02error.push_back(lgWerr[1]);	
       LG03error.push_back(lgWerr[2]);	
       LG04error.push_back(lgWerr[3]);	
       LG05error.push_back(lgWerr[4]);	
       LG06error.push_back(lgWerr[5]);	
       LG07error.push_back(lgWerr[7]);	
       LG08error.push_back(lgWerr[8]);	
       LG09error.push_back(lgWerr[9]);	
       LG10error.push_back(lgWerr[10]);	
       LG11error.push_back(lgWerr[11]);
      }

}

    if (east == 1) {
      for( int i =0; i<linum; i++){
       pmtTree->GetEntry(i); 
       PMTadc[i]=scintE[pmt];       
       LG01prob[i]=lgEprob[0];	     //why no std::maxtrix<double>? why c++ EHY?!?!!
       LG02prob[i]=lgEprob[1];	
       LG03prob[i]=lgEprob[2];	
       LG04prob[i]=lgEprob[3];	
       LG05prob[i]=lgEprob[4];	
       LG06prob[i]=lgEprob[5];	
       LG07prob[i]=lgEprob[7];	
       LG08prob[i]=lgEprob[8];	
       LG09prob[i]=lgEprob[9];	
       LG10prob[i]=lgEprob[10];	
       LG11prob[i]=lgEprob[11];	
       LG01error[i]=lgEerr[0];	
       LG02error[i]=lgEerr[1];	
       LG03error[i]=lgEerr[2];	
       LG04error[i]=lgEerr[3];	
       LG05error[i]=lgEerr[4];	
       LG06error[i]=lgEerr[5];	
       LG07error[i]=lgEerr[7];	
       LG08error[i]=lgEerr[8];	
       LG09error[i]=lgEerr[9];	
       LG10error[i]=lgEerr[10];	
       LG11error[i]=lgEerr[11];
       }
    }
}


double ChiSquaredBrah(const double *LGFitParam){
	const Double_t lg0=LGFitParam[0];
	const Double_t lg1=LGFitParam[1];
	const Double_t lg2=LGFitParam[2];
	const Double_t lg3=LGFitParam[3];
	const Double_t lg4=LGFitParam[4];
	const Double_t lg5=LGFitParam[5];
	const Double_t lg6=LGFitParam[6];
	const Double_t lg7=LGFitParam[7];
	const Double_t lg8=LGFitParam[8];
	const Double_t lg9=LGFitParam[9];
	const Double_t lg10=LGFitParam[10];
	const Double_t lg11=LGFitParam[11];
	//const Double_t lg12=LGFitParam[12];  pow(lg12,2)*

        Double_t chiBrah;

        int entnum = LG01prob.size();
	
     for (int i =0; i<entnum; ++i){
	if(xpos[i]*xpos[i]+ypos[i]*ypos[i] < 25 && LG02error[i]>0.0001){//gain calibration use fiducial circle of 5mm radius. 
       chiBrah+=pow(gain0[i]-(lg0*LG01prob[i]+lg1*LG02prob[i]+lg2*LG03prob[i]+lg3*LG04prob[i]+lg4*LG05prob[i]+lg5*LG06prob[i]+lg6*LG07prob[i]+lg7*LG08prob[i]+lg8*LG09prob[i]+lg9*LG10prob[i]+lg10*LG11prob[i]+lg11*LG12prob[i]),2)/entnum/(pow(lg0*LG01error[i],2.0)+pow(lg1*LG02error[i],2.0)+pow(lg2*LG03error[i],2.0)+pow(lg3*LG04error[i],2.0)+pow(lg4*LG05error[i],2.0)+pow(lg5*LG06error[i],2.0)+pow(lg6*LG07error[i],2.0)+pow(lg7*LG08error[i],2.0)+pow(lg8*LG09error[i],2.0)+pow(lg9*LG10error[i],2.0)+pow(lg10*LG11error[i],2.0)+pow(lg11*LG12error[i],2.0));
	     
		}
	} 

	
     return chiBrah;  //chi^2 

}



void PMTMapFit(TTree * pmtTree, TTree * fitTree,TTree * flatpmtTree){
 int linum=pmtTree->GetEntries();
 	
 TH2F * qmap = new TH2F("qmap","adc",30,-60,60,30,-60,60);//,30,0,10);
  Float_t Wpos[2];
  Float_t scintW[4];  
 //; 
 pmtTree->SetBranchAddress("ScintW",&scintW);
//flatpmtTree->SetBranchAddress("ScintWflat",&scintW);
 pmtTree->SetBranchAddress("Wpos",&Wpos);
  vector<double> xWpos;
  vector<double> yWpos;
  vector<double> ScintWflat;
  xWpos.reserve(linum);
  yWpos.reserve(linum);
  ScintWflat.reserve(linum);
 
 // Fill vector values
  for (int i = 0; i< linum; ++i){
 pmtTree->GetEntry(i);
 flatpmtTree->GetEntry(i);
 xWpos.push_back(Wpos[0]);
 yWpos.push_back(Wpos[1]);
 ScintWflat.push_back(scintW[0]);
}

std::vector<double> fidIN(gain0.size());
 double gainsum=0; 
double gainsize=0;
  int pmtnum=0; //0-3 
  int east = 0;
   Double_t gain=0;
  SetLGMap(pmtTree,flatpmtTree,east,pmtnum);
	//int gainsize=0;
     for(int i =0; i < linum; i++){
       if(xpos[i]*xpos[i]+ypos[i]*ypos[i] < 25 && gain0[i]<2 && gain0[i]>0){ //&& LG02error[i]>0.0001){
       fidIN.push_back(1);gainsum+=gain0[i];++gainsize;}
        else fidIN.push_back(0);
      }
    
     


double gainmean = gainsum / gainsize;
double gainvar=0;
for(int i =0; i < linum; i++){
     if(xpos[i]*xpos[i]+ypos[i]*ypos[i] < 25 && gain0[i]<2 && gain0[i]>0){ //&& LG02error[i]>0.0001){
       gainvar+=pow(gain0[i]-gainmean,2);}
        
      }

double gainstdev=sqrt(gainvar/gainsize);
 cout<<gainmean<<"   "<<gainstdev<<"  "<<gainsize<<"\n";   ///test how good this is. 

TH1F * gainhist = new TH1F("gain","gain hist",30,-1,4);//,30,0,10);
 
for (int i = 0; i< linum; ++i){
 //pmtTree->GetEntry(i);
    if(xpos[i]*xpos[i]+ypos[i]*ypos[i] < 25 && gain0[i]<2 && gain0[i]>0){ //&& LG02error[i]>0.0001){
      
  gainhist->Fill(gain0[i]);
  
      }
     }

gainhist->Draw("colz");










   const double pi = 3.1415926535897;  //strawberry
  
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



 
 fitTree->SetBranchAddress("FitParamE1",&fp1E);
 fitTree->SetBranchAddress("FitParamW1",&fp1W);
 fitTree->SetBranchAddress("FitParamE2",&fp2E);
 fitTree->SetBranchAddress("FitParamW2",&fp2W);
 fitTree->SetBranchAddress("FitParamE3",&fp3E);
 fitTree->SetBranchAddress("FitParamW3",&fp3W);
 fitTree->SetBranchAddress("FitParamE4",&fp4E);
 fitTree->SetBranchAddress("FitParamW4",&fp4W);
 fitTree->SetBranchAddress("OffsetE",&offsetE);
 fitTree->SetBranchAddress("OffsetW",&offsetW);

	//get fit parameters. or else...
 fitTree->GetEntry(0);

    //create minimizer   Minuit 2 with migrad
 ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad");

   min->SetMaxFunctionCalls(10000000); // for Minuit/Minuit2 
   min->SetMaxIterations(1000);  // for GSL 
   min->SetTolerance(1);
   min->SetPrintLevel(1);
	

  
  
  
 	unsigned int dimnum=12;//12;

  ROOT::Math::Functor f(&ChiSquaredBrah,dimnum);

  //step size 

  double step[12]={.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01,.01};  
  const double * FP1W=fp1W;

  
///////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////

/////////////////////////////////////////////
/* 
  min->SetFunction(f);
  min->SetVariable(0,"lg0",.1,0.01);
  min->SetVariable(1,"lg1",1,1);
  min->SetVariable(2,"lg2",8,1);
  min->SetVariable(3,"lg3",9,1);//fp1W[3],step[3]);
  min->SetVariable(4,"lg4",2,1);//fp1W[4],step[4]);
  min->SetVariable(5,"lg5",2,1);//fp1W[5],step[5]);
  min->SetVariable(6,"lg6",2,1);//fp1W[6],step[6]);
  min->SetVariable(7,"lg7",0.1,1);//fp1W[7],step[7]);
  min->SetVariable(8,"lg8",0.1,1);//fp1W[8],step[8]);
  min->SetVariable(9,"lg9",0.1,1);//fp1W[9],step[9]);
min->SetVariable(10,"lg10",0.1,1);//fp1W[10],step[10]);
min->SetVariable(11,"lg11",0.1,1);//fp1W[11],step[11]);
  //min->SetVariable(12,"lg12",1,1);  
 
/* min->SetVariable(0,"lg0",fp1W[0],step[0]);
  min->SetVariable(1,"lg1",fp1W[1],step[1]);
  min->SetVariable(2,"lg2",fp1W[2],step[2]);
  min->SetVariable(3,"lg3",fp1W[3],step[3]);*/

   //////////////////////////////////////////////////////////////////////////////////////////                                               
 /*
   
 min->Minimize();
  
  const double *lgcoeff = min->X();
  std::cout << "Minimum: f(" << lgcoeff[0] << "," << lgcoeff[0] << "...): " << min->MinValue()  << std::endl;
  
  
 for (int i = 0; i< linum; ++i){
 //pmtTree->GetEntry(i);
   if(xpos[i]*xpos[i]+ypos[i]*ypos[i] < 25 && LG02error[i]>0.0001){
      
  qmap->Fill(xWpos[i],yWpos[i],ScintWflat[i]);
  
      }
     }

qmap->Draw("colz");






TH2F * qfit = new TH2F("qfit","fit adc",30,-60,60,30,-60,60);//,30,0,10);
 
 
 //; 
 
  TCanvas * c2= new TCanvas();                       //c2->cd();
 // Fill Histogram.
  for (int i = 0; i< linum; ++i){
 //pmtTree->GetEntry(i);
    if(xpos[i]*xpos[i]+ypos[i]*ypos[i] < 25 && LG02error[i]>0.0001){
      qfit->Fill(xWpos[i],yWpos[i],(lgcoeff[0]*LG01prob[i]+lgcoeff[1]*LG02prob[i]+lgcoeff[2]*LG03prob[i]+lgcoeff[3]*LG04prob[i]+lgcoeff[4]*LG05prob[i]+lgcoeff[5]*LG06prob[i]+lgcoeff[6]*LG07prob[i]+lgcoeff[7]*LG08prob[i]+lgcoeff[8]*LG09prob[i]+lgcoeff[9]*LG10prob[i]+lgcoeff[10]*LG11prob[i]+lgcoeff[11]*LG12prob[i]));
  
  
      }
     }


c2->cd();
qfit->Draw("colz");

 


 TH2F * qdiff = new TH2F("qdiff","Qadc - Fit adc",30,-60,60,30,-60,60);//,30,0,10);
 
  TCanvas * c3= new TCanvas();                       //c2->cd();
 // Fill Histogram.
  for (int i = 0; i< linum; ++i){
 //pmtTree->GetEntry(i);
   if(xpos[i]*xpos[i]+ypos[i]*ypos[i] < 25 && LG02error[i]>0.0001){
      qdiff->Fill(xWpos[i],yWpos[i],ScintWflat[i]-(lgcoeff[0]*LG01prob[i]+lgcoeff[1]*LG02prob[i]+lgcoeff[2]*LG03prob[i]+lgcoeff[3]*LG04prob[i]+lgcoeff[4]*LG05prob[i]+lgcoeff[5]*LG06prob[i]+lgcoeff[6]*LG07prob[i]+lgcoeff[7]*LG08prob[i]+lgcoeff[8]*LG09prob[i]+lgcoeff[9]*LG10prob[i]+lgcoeff[10]*LG11prob[i]+lgcoeff[11]*LG12prob[i]));
     }
   }
c3->cd();
qdiff->Draw("colz");

*/
 

//std::cout<<ChiSquaredBrah(FP1W)<<"\n";
}

 

int main(int argc, char **argv)
{
   
   TFile *myFile = TFile::Open("$UCNAOUTPUTDIR/hists/pmtprob_22000.root");
   TTree *pmtTree = (TTree*)myFile->Get("pmtTree");
   TTree *fitTree = (TTree*)myFile->Get("fitTree");
   TTree *flatpmtTree = (TTree*)myFile->Get("flatpmtTree");
     TApplication app("Qadc Stuff", &argc, argv);

   PMTMapFit(pmtTree,fitTree,flatpmtTree);
  
   app.Run();

}












