//g++ -o PMTfit LGMapFit.cpp `root-config --cflags --glibs`
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
#include "TPaveText.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TError.h"


using namespace std;



vector<double> PMTadc;
//a dozen LG probability vectors.        
vector<double> LG01prob;vector<double> LG02prob;vector<double> LG03prob;vector<double> LG04prob;vector<double> LG05prob;vector<double> LG06prob;vector<double> LG07prob;vector<double> LG08prob;vector<double> LG09prob;vector<double> LG10prob;vector<double> LG11prob;vector<double> LG12prob;  
//a dozen LG error vectors 
vector<double> LG01error;vector<double> LG02error;vector<double> LG03error;vector<double> LG04error;vector<double> LG05error;
vector<double> LG06error;vector<double> LG07error;vector<double> LG08error;vector<double> LG09error;vector<double> LG10error;vector<double> LG11error;vector<double> LG12error;




void SetLGMap (TTree * pmtTree,int east,int pmt) {

  double lgEprob[12];
  double lgWprob[12];
  double lgEerr[12];
  double lgWerr[12];
  double scintE[4];
  double scintW[4];
 pmtTree->SetBranchAddress("LGEprob",&lgEprob);
 pmtTree->SetBranchAddress("LGWprob",&lgWprob);
 pmtTree->SetBranchAddress("LGEerror",&lgEerr);
 pmtTree->SetBranchAddress("LGWerror",&lgWerr);
 pmtTree->SetBranchAddress("ScintE",&scintE);
 pmtTree->SetBranchAddress("ScintW",&scintW);
 int linum=pmtTree->GetEntries();
 //PMTadc=vector<double>();
//define vector lengths
     PMTadc.reserve(linum);LG01prob.reserve(linum); LG02prob.reserve(linum); LG03prob.reserve(linum); LG04prob.reserve(linum); LG05prob.reserve(linum); LG06prob.reserve(linum); LG07prob.reserve(linum); LG08prob.reserve(linum); LG09prob.reserve(linum); LG10prob.reserve(linum); LG11prob.reserve(linum); LG12prob.reserve(linum);LG01error.reserve(linum); LG02error.reserve(linum); LG03error.reserve(linum); LG04error.reserve(linum); LG05error.reserve(linum); LG06error.reserve(linum); LG07error.reserve(linum); LG08error.reserve(linum); LG09error.reserve(linum); LG10error.reserve(linum); LG11error.reserve(linum); LG12error.reserve(linum);

  if (east = 0) {
   for( int i =0; i<linum; i++){
   /*    pmtTree->GetEntry(i);
       PMTadc[i]=scintW[pmt];
      LG01prob[i]=lgWprob[0];	     //why no std::maxtrix<double>? why c++ WHY?!?!!
       LG02prob[i]=lgWprob[1];	     //preprocessor loop?
       LG03prob[i]=lgWprob[2];	
       LG04prob[i]=lgWprob[3];	
       LG05prob[i]=lgWprob[4];	
       LG06prob[i]=lgWprob[5];	
       LG07prob[i]=lgWprob[7];	
       LG08prob[i]=lgWprob[8];	
       LG09prob[i]=lgWprob[9];	
       LG10prob[i]=lgWprob[10];	
       LG11prob[i]=lgWprob[11];	
       LG01error[i]=lgWerr[0];	
       LG02error[i]=lgWerr[1];	
       LG03error[i]=lgWerr[2];	
       LG04error[i]=lgWerr[3];	
       LG05error[i]=lgWerr[4];	
       LG06error[i]=lgWerr[5];	
       LG07error[i]=lgWerr[7];	
       LG08error[i]=lgWerr[8];	
       LG09error[i]=lgWerr[9];	
       LG10error[i]=lgWerr[10];	
       LG11error[i]=lgWerr[11];	*/
      }

}

    if (east = 1) {
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
	
        Double_t chisquaredBrah;
        int entnum = std::distance(LG01prob.begin(),LG01prob.end());
     for (int i =0; i<entnum; ++i){
	
       chisquaredBrah+=pow((PMTadc[i]-(lg0*LG01prob[i]+lg1*LG02prob[i]+lg2*LG03prob[i]+lg3*LG04prob[i]+lg4*LG05prob[i]+lg5*LG06prob[i]+lg6*LG07prob[i]+lg7*LG08prob[i]+lg8*LG09prob[i]+lg9*LG10prob[i]+lg10*LG11prob[i]+lg11*LG12prob[i]))/(pow(lg0*LG01error[i],2.0)+pow(lg1*LG02error[i],2.0)+pow(lg2*LG03error[i],2.0)+pow(lg3*LG04error[i],2.0)+pow(lg4*LG05error[i],2.0)+pow(lg5*LG06error[i],2.0)+pow(lg6*LG07error[i],2.0)+pow(lg7*LG08error[i],2.0)+pow(lg8*LG09error[i],2.0)+pow(lg9*LG10error[i],2.0)+pow(lg10*LG11error[i],2.0)+pow(lg11*LG12error[i],2.0)),2.0);   
  }
     return chisquaredBrah;

}



void PMTMapFit(TTree * pmtTree, TTree * fitTree){
   
  int linum=pmtTree->GetEntries();
	cout<<linum<<"\n";
  int pmtnum=0; //0-3 
  int east = 0;
   
  SetLGMap(pmtTree,east,pmtnum);


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

   min->SetMaxFunctionCalls(100); // for Minuit/Minuit2 
  // min->SetMaxIterations(1000);  // for GSL 
   min->SetTolerance(0.001);
   min->SetPrintLevel(1);
	

  
  
  ////import pmtTree, 
  SetLGMap (pmtTree,east,pmtnum);
 	unsigned int dimnum=12;

  ROOT::Math::Functor f(&ChiSquaredBrah,dimnum);

  //step size (should be automatic in a good minimizer)

  double step[12]={0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01};  
 
  min->SetFunction(f);
  min->SetVariable(0,"lg0",fp1W[0],step[0]);
  min->SetVariable(1,"lg1",fp1W[1],step[1]);
  min->SetVariable(2,"lg2",fp1W[2],step[2]);
  min->SetVariable(3,"lg3",fp1W[3],step[3]);
  min->SetVariable(4,"lg4",fp1W[4],step[4]);
  min->SetVariable(5,"lg5",fp1W[5],step[5]);
  min->SetVariable(6,"lg6",fp1W[6],step[6]);
  min->SetVariable(7,"lg7",fp1W[7],step[7]);
  min->SetVariable(8,"lg8",fp1W[8],step[8]);
  min->SetVariable(9,"lg9",fp1W[9],step[9]);
  min->SetVariable(10,"lg10",fp1W[10],step[10]);
  min->SetVariable(11,"lg11",fp1W[11],step[11]);
 
  min->Minimize();

  const double *xs = min->X();
  std::cout << "Minimum: f(" << xs[0] << "," << xs[1] << "...): " << min->MinValue()  << std::endl;

}

 

int main()
{
   TFile *myFile = TFile::Open("$UCNAOUTPUTDIR/hists/pmtprob_22770.root");
   TTree *pmtTree = (TTree*)myFile->Get("pmtTree");
   TTree *fitTree = (TTree*)myFile->Get("fitTree");
 

   PMTMapFit(pmtTree,fitTree);
  

}












