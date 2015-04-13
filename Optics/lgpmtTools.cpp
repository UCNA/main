//#include <fstream>
#include <iostream>
#include <stdlib.h>
//#include <stdio.h>
//#include <unistd.h> 
#include <math.h>
#include <cmath> 
#include <TFile.h>
#include <TH1F.h>
#include <TTreeReader.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TTreeReaderValue.h>
//#include "/home/cmswank/Documents/ucna/main/Scripts/lgprobmap/pmtprobstuff.h"
//#include "/home/cmswank/Documents/ucna/main/Scripts/lgprobmap/lgpmtTools.h"
#include "pmtprobstuff.h"
#include "lgpmtTools.h"

using namespace std;

///This function returns the probability of entering a lightguide given a position, 
	///lgnum is the light guide number, 0 to 11, think analog clock hour positions!! 
double lightguideprob(double xpos, double ypos,int lgnum){
	double lightprob, thet1, thet2;
	const double Radius=75;//*sqrt(.6); what radius??? scintilator or wire chamber.	                       
	const double pi = 3.1415926535897;  //pi
	const double p1x=0;
	const double p1y=Radius;
	const double p2x=Radius*sin(pi/6);
	const double p2y=p1y*cos(pi/6);
	const double thetaC = 40*pi/180; //critical angle in radians 
	///My model started with theta designating the angle from the vertical, 
    /// therefore x=r*sin(theta), etc... whoops.	
    if  (sqrt(pow(xpos,2) + pow(ypos,2)) > Radius){
        //cout<<"Impact position missed Scintillator, try again"<<"\n";
        return 0.00; 
    }
    if (lgnum>=0 && lgnum<12){ //analog clock hour positions are lightguide positions!! 
        
        double xp, yp, r, theta;
        r=sqrt(pow(xpos,2)+pow(ypos,2));  
        theta=-atan2(xpos,ypos);         
        xp=r*sin(theta+lgnum*pi/6);   //rotate to the lg to the right
        yp=r*cos(theta+lgnum*pi/6); 
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

    lightprob=fabs((thet1-thet2)/2/pi);///I can't remember which way it integrates, so I take abs value. 
    return lightprob;		
}	



