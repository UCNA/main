//#include "G4toPMT.hh"
//#include "CalDBSQL.hh"
//#include <TH1F.h>
//#include <TLegend.h>
//#include <TFitResult.h> // v5.27
//#include <TF1.h>

#include <iostream>
#include <fstream>
#include <string>
#include "Studies/BetaSpectrum.hh"


const double    Q           = 782.344;              /// end point KE        			(30)
const int 		n 			= 1E4;					/// number of points to simulate


/// beta phase space integral
//double PhaseSpaceIntegral(const double *val, const double *par) {
double PhaseSpaceIntegral() {
	double integral = 0;
	for (double KE = 0; KE < Q; KE += Q/(double)n)
	{
		double P = correctedBetaSpectrum(KE) / n * 1.2078377907569861626768670;
		integral += P;
		printf("%f\t%.12f\n", KE, P);
	}

	//printf("%16.25f\n", 1/integral);
	return integral;
}



int main (int argc, char ** argv) {
	PhaseSpaceIntegral();
}
