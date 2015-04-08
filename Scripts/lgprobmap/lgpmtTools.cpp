


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
		if (fabs(thet1)>thetaC){
		thet1=sign(thet1)*thetaC;
		}
		if (fabs(thet2)>thetaC){
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

