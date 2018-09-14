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
#include <gsl/gsl_mode.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_sf_ellint.h> //double gsl_sf_ellint_E(double phi, double k, gsl_mode_t mode);

using namespace std;

///This function returns the probability of entering a lightguide given a position, 
	///lgnum is the light guide number, 0 to 11, think analog clock hour positions!! 
double lightguideprob(double xpos, double ypos,int lgnum){
	double lightprob, thet1, thet2, xp, yp, r, dist, lgsorb;
	const double Radius=75*sqrt(.6); //what radius??? scintilator or wire chamber.	                       
	const double pi = 3.1415926535897;  //pi
	const double p1x=0;
	const double p1y=Radius;
	const double p2x=Radius*sin(-pi/6);
	const double p2y=p1y*cos(-pi/6);
	const double thetaC = 40*pi/180; //critical angle in radians 
	///My model started with theta designating the angle from the vertical, 
    /// therefore x=r*sin(theta), etc... whoops.	
    if  (sqrt(pow(xpos,2) + pow(ypos,2)) > Radius){
        //cout<<"Impact position missed Scintillator, try again"<<"\n";
        return 0.00; 
    }
    if (lgnum>=0 && lgnum<12){ //analog clock hour positions are lightguide positions!! 
        
        double theta;
        r=sqrt(pow(xpos,2)+pow(ypos,2));  
        theta=-atan2(xpos,ypos);         
        xp=r*sin(theta+lgnum*pi/6+pi/12);   //rotate to the lg number!
        yp=r*cos(theta+lgnum*pi/6+pi/12); 
        thet1=-(atan2((p1x-xp),(p1y-yp)));
        thet2=-(atan2((p2x-xp),(p2y-yp)));
    }
    else{ cout<<"please choose a valid light guide number 0-11.\n";
        return 0;
    }
    // computerified heaviside in Kevin's integral: if theta gets bigger than thetaC then make it equal, mind the sign. 
    if (fabs(thet1)>thetaC){
        thet1=sign(thet1)*thetaC;
    }
    if (fabs(thet2)>thetaC){
        thet2=sign(thet2)*thetaC;
    }
    dist=sqrt(pow(Radius-yp,2)+pow(xp,2));
    lightprob=fabs((thet1-thet2)/2/pi);///I can't remember which way it integrates, so I take abs value. 
    
    //turns out we need to include some absorption. SO...
  ///r0a
     ///mine has a minus because my dTheta is opposite of kevins. 
   // lgsorb=-r*(sin(thet2)-sin(thet1))-Radius*(sign(cos(thet2))*gsl_sf_ellint_E(thet2, r/Radius, 2)-sign(cos(thet1))*gsl_sf_ellint_E(thet1,r/Radius, 2)); //0 double precision, 1 signle precision, 2 approx. 
 
  return lightprob*exp(-dist/28);		//1/75 is placeholder scattering rate (1/distance)
}	



