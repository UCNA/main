//g++ -o lgprob LGProbError.cpp pmtprobstuff.cpp lgpmtTools.cpp `root-config --cflags --glibs`
//#include <fstream>
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
#include "pmtprobstuff.h"
#include "lgpmtTools.h"
#include "TLatex.h"
#include "TFormula.h"
#include "TAxis.h"
#include "TPave.h"
#include "TPaveStats.h"
#include "TBox.h"
#include "TText.h"
#include "TH2.h"
#include "TF1.h"
#include "TH2F.h"
#include "TPaveText.h"


using namespace std;


//a bit much for this. root probably has a function for this. 
int sign(double value) {
    return ( value > 0) - (value < 0);
}


//This probably gives a probability to a pmt given the offset and coupling coefficients, and pmt number.
double PMTprob(double xpos,double ypos, double LGFitParam[],int pmtnum){
	const double pi = 3.1415926535897;  //pi
	double xp, yp; 
	double pmtprob=0;
		
	if (pmtnum>0 && pmtnum<5){
	xp=xpos-LGFitParam[0];///center LGs w.r.t. wirechamber. 
	yp=ypos-LGFitParam[1]; //
	}	
	
	else {cout<<"enter a valid pmt# (they are 1,2,3 or 4)\n"; return 0;}
	for (int i = 0; i<12; ++i) pmtprob+=LGFitParam[i+2+12*(pmtnum-1)]*lightguideprob(xp,yp,i);  //sum over all lightguides 
												    //with LG FIt Parameter in the 													   //right spots
	return pmtprob;
}

//generate pmt uncertianty 
double PMTerror(double xerr, double yerr,double xpos,double ypos, double LGFitParam[],int pmtnum){
	
	//random number class
	
	TRandom3 * alThor = new TRandom3(0);
       
	// generate a gaussian distributed number with mu=xpos,ypos,  sigma=xerr,yerr (default values)

	double PMTtot[10];  //ten points enough??? NO!!!
	
	for(int i = 0; i<10; ++i){
	
	
	PMTtot[i]=PMTprob(alThor->Gaus(xpos,xerr),alThor->Gaus(ypos,yerr),LGFitParam,pmtnum);	
	}

	return (TMath::RMS(&PMTtot[0],&PMTtot[9]));
}

double LGerror(double xerr, double yerr,double xpos,double ypos,int lgnum){
	
	//random number class
	
	TRandom3 * alThor = new TRandom3(0);
       
	// generate a gaussian distributed number with mu=xpos,ypos,  sigma=xerr,yerr (default values)

	double PMTtot[10];  //ten points enough??? NO!!!
	
	for(int i = 0; i<100; ++i){
	
	
	PMTtot[i]=lightguideprob(alThor->Gaus(xpos,xerr),alThor->Gaus(ypos,yerr),lgnum);	
	}

	return (TMath::RMS(&PMTtot[0],&PMTtot[9]));
}


