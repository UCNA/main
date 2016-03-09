#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h> 
#include <stdlib.h>
#include <math.h>
#include <cmath> 
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
//#include "TTreeReader.h"          /// only in newer Root
//#include "TTreeReaderValue.h"     /// only in newer Root
using namespace std;
 
//This is a collection of functions that are useful for light guide map creation. 
///For an example of use look in main. 





void twoDvec2txt(const char* filename,double number1[],double number2[], double value[],int length){     
 
  ofstream myfile;
  myfile.open (filename);
 	for(int i=0; i<length; i++) {
  myfile << number1[i] <<" "<< number2[i] << " " << value[i]<<"\n"; 
  }
  myfile.close();
}

class PostextToArray   //dumb class, I only could get linenumber to work right. 
{
 	
  string line;
  
  public:
  int linenumber(string str);
  
   //// double *lgpos(const int ln);  /doesn't work in c++, cause i hate c++, 
  //void setpos(string str);        /doesn't work this way in c++ . 
  };

int PostextToArray::linenumber(string str) { 

  ifstream file (str.c_str());
  if (file.is_open())
  {
    int i = 0;
    while ( getline (file,line) )
    {
      ++i; 	
     }
    file.close();
    return i-1;	
     }

  else {
	cout << "Unable to open file, Sorry :o( \n";
	return 0;
	}  
}

/// sign function for double input  sign(>0)=1, sign(<0)=-1, sign(0)=0;
int sign(double value) {
    return ( value > 0) - (value < 0);
}




///This function returns the probability of entering a lightguide given a position, 
	///lgnum is the light guide number, 0 is for center, 1 is for the left, 2 is for the right. 
double lightguideprob(double xpos, double ypos,int lgnum){
	
	
	double lightprob;
	double thet1;
	double thet2;
	const double Radius=75;//*sqrt(.6); what radius??? scintilator or wire chamber.	                       
	const double pi = 3.1415926535897;  //pi
	const double p1x=Radius*sin(-pi/12);
	const double p1y=Radius*cos(pi/12);
	const double p2x=-p1x;
	const double p2y=p1y;
	const double thetaC = 40*pi/180;	//critical angle in radians 
	
	
	
		
		if  (sqrt(pow(xpos,2) + pow(ypos,2)) > Radius){
		cout<<"Impact position missed Scintillator, try again"<<"\n";
		return 0.00; 
		}

if (lgnum==2){//%on da right
	double xright;
	double yright;
	double r;
	double theta;
	r=sqrt(pow(xpos,2)+pow(ypos,2));  
	theta=-atan2(xpos,ypos);
	xright=r*sin(theta+pi/6);   //rotate lg to the right
	yright=r*cos(theta+pi/6); 
	thet1=-(atan2((p1x-xright),(p1y-xright)));
	thet2=-(atan2((p2x-xright),(p2y-yright)));
}
else if (lgnum==1){//%relative to lightguide on the left;
	
	double xleft;
	double yleft;
	double r;
	double theta;
	r=sqrt(pow(xpos,2)+pow(ypos,2));  
	theta=-atan2(xpos,ypos);
	xleft=r*sin(theta-pi/6); ///rotate lg to the left
	yleft=r*cos(theta-pi/6); 
	thet1=-(atan2((p1x-xleft),(p1y-yleft)));
	thet2=-(atan2((p2x-xleft),(p2y-yleft)));
}
	
	//Kevin's equation for phi1 and phi2 (find all accepted angles) 
else if (lgnum==0){	
	thet1=-(atan2((p1x-xpos),(p1y-ypos)));
	thet2=-(atan2((p2x-xpos),(p2y-ypos)));
}
else{ 
	cout<<"please choose a valid light guide number.\n";
	return 0;
 	 }
	//computerified heaviside in Kevin's integral: if theta gets bigger than thetaC then make it equal, mind the sign. 
		if (abs(thet1)>thetaC){
		thet1=sign(thet1)*thetaC;
		}
		if (abs(thet2)>thetaC){
		thet2=sign(thet2)*thetaC;
		}

	lightprob=fabs((thet1-thet2)/2/pi);///I can't remember which way it integrates, so I take abs value. 
	return lightprob;		
}	


//This gives a probability to a pmt given the offset and coupling coefficients, and pmt number.
double PMTprob(double xpos,double ypos, double offsetx,double offsety, double lga, double lgb, double lgc,int pmtnum){
	const double pi = 3.1415926535897;  //pi
	double xp, yp, r, theta, pmtprob;
	
	r=sqrt(pow(xpos-offsetx,2)+pow(ypos-offsety,2));  //r and theta in terms of scintillator center. 
	theta=-atan2(xpos-offsetx,ypos-offsety);
	
	if (pmtnum>0 && pmtnum<5){
	xp=r*sin(theta-3*pi/12-(pmtnum-1)*pi/2); ///rotate to correct orientation of LG's w.r.t. wirechamber. 
	yp=r*cos(theta-3*pi/12-(pmtnum-1)*pi/2); 
	}	
	else {cout<<"enter a valid pmt# (they are 1,2,3 or 4)\n";}

	pmtprob=lga*lightguideprob(xp,yp,1)+lgb*lightguideprob(xp,yp,0)+lgc*lightguideprob(xp,yp,2); 
	
	return pmtprob;

}






	
////example of use..
//import wirechamber position from root.
//create lightguide probability from positions.
//export positions and probability to txt file. 
 	
int main()
{	
	const double pi = 3.1415926535897;  //pi
	///The following code gets the data from root file and exports it to a text file... or maybe I'll export it to a new root file with a branch name thats really hard to spell. 

   TFile *myFile2 = TFile::Open("/home/cmswank/G4Work/output/thinfoil_Xe135_3-2+/analyzed_1.root");
   TTree *anaTree = (TTree*)myFile2->Get("anaTree");

   // Create a TTreeReader for the tree, for instance by passing the
   // TTree's name and the TDirectory / TFile it is in.
     //TTreeReader myReader("ntuple", myFile);
   //TTreeReader myReader2("anaTree", myFile2); //TTreeReader doesn't seem to work with nicely with g++ 4.8 
   
 
   Double_t hitpos[6];
    
   Int_t linum=anaTree->GetEntries();
   anaTree->SetBranchAddress("MWPCPos",&hitpos);



	
	string lgps;
  	double xpos[linum];
	double ypos[linum];
	double lgp[linum];  	
    	int ii = 0;
    	while (ii<linum)
    	{   
	    anaTree->GetEntry(ii);
	    if (hitpos[2]==0){  //west side or east side? this is important. and WRONG VERY WRONG 
	    
     	    xpos[ii] = 10*hitpos[4];  //change units from cm to mm (10*hitpos). 
  	    ypos[ii] = 10*hitpos[5];}	
   	    else
            {
	    xpos[ii] = 10*hitpos[1];
  	    ypos[ii] = 10*hitpos[2];	
             }
	    

            //lgp[ii] = lightguideprob(xpos[ii],ypos[ii],0); ///LG Probability function
	    lgp[ii]=PMTprob(xpos[ii],ypos[ii],0,0,1,1,1,1)+PMTprob(xpos[ii],ypos[ii],0,0,1,1,1,2)+PMTprob(xpos[ii],ypos[ii],0,0,1,1,1,3)+PMTprob(xpos[ii],ypos[ii],0,0,1,1,1,4); ///PMT Probablility function
		
		if (lgp[ii]>1.0000)
		{      //all errors tend to be less than 1% e.g. prob=1.00X where X is < 5.
			if (lgp[ii]>1.05){ cout<<"Warning, Probabililty>1, something isn't right here.\n";}
		 lgp[ii]=1.000;			
		}
		++ii;	 	
  	   }   
 	//
 


//Export function.
    twoDvec2txt("lightprob.txt",xpos,ypos,lgp,linum);

   return 0;

}
