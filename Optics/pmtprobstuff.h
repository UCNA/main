
#ifndef PMTPROBSTUFF_H
#define PMTPROBSTUFF_H

//sign function
int sign(double value);

double PMTprob(double xpos,double ypos, double LGFitParam[],int pmtnum);

double PMTerror(double xerr, double yerr,double xpos,double ypos, double LGFitParam[],int pmtnum);

double LGerror(double xerr, double yerr,double xpos,double ypos,int lgnum);

#endif
