#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h> 
#include <stdlib.h>
#include <math.h>
#include <cmath> 
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
	const double Radius=75*sqrt(.6);	                       
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
	thet1=(atan((p1x-xright)/(p1y-xright)));
	thet2=(atan((p2x-xright)/(p2y-yright)));
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
	thet1=(atan((p1x-xleft)/(p1y-yleft)));
	thet2=(atan((p2x-xleft)/(p2y-yleft)));
}
	
	//Kevin's equation for phi1 and phi2 (find all entrence angles) atan works becuase we are limited in theta. atan2 is used later.  
else if (lgnum==0){	
	thet1=(atan((p1x-xpos)/(p1y-ypos)));
	thet2=(atan((p2x-xpos)/(p2y-ypos)));
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

	lightprob=fabs((thet1-thet2)/2/pi);
	return lightprob;	

	
}	
	
////example of use..
//import impact position from txt.
//create lightguide probability from positions.
//export positions and probability to txt file. 
 	
int main()
{	
	PostextToArray pos; ///class died a birth. only gets line numbers, kinda dumb. 
/// Get line numbers;			   
 	int linum=pos.linenumber("hitpos.txt");  //see this is all that works. 

//import hit positions from text file. 	I tried to include this in a class, but it doesn't save any space because you need the array size specified before the class is called. 
  	
	string lgps;
  	double xpos[linum];
	double ypos[linum];
	double lgp[linum];
  	ifstream file2 ("hitpos.txt");	  
  	if (file2.is_open())
  	{
    	int ii = 0;
    	while ( getline (file2,lgps) )
    	{   
     	    xpos[ii] = atof(lgps.c_str());
   	    std::size_t pos = lgps.find(" ");	
  	    ypos[ii] = atof(lgps.substr(pos).c_str());	
   	    lgp[ii] = lightguideprob(xpos[ii],ypos[ii],0); ///Probability function
  	    //cout<<xpos[ii]<<" "<<ypos[ii]<<"\n";
	    ++ii; 	
  	   }   
 	  file2.close();
 	 }
	else {cout<<"please read from the same file\n";}


//Export function.
    twoDvec2txt("lightprob.txt",xpos,ypos,lgp,linum);//,xpos,ypos,lightprob,possize);

   return 0;

}
