//#include "G4toPMT.hh"
//#include "CalDBSQL.hh"
//#include <TH1F.h>
//#include <TLegend.h>
//#include <TFitResult.h> // v5.27
//#include <TF1.h>

#include <iostream>
#include <fstream>
#include <string>
<<<<<<< HEAD
#include "Studies/BetaSpectrum.hh" 
//#include "ucnG4_dev/include/bmPrimaryGeneratorAction.hh"
=======
#include "BetaSpectrum.hh"
>>>>>>> c40651173c2f4f0e350b4812aee2e1c2f305f846


const double    Q           = 782.344;              /// end point KE        			(30)
const int 		n 			= 1E4;					/// number of points to simulate


/// beta phase space integral
//double PhaseSpaceIntegral(const double *val, const double *par) {
double PhaseSpaceIntegral() {
	double integral = 0;
	for (double KE = 0; KE < Q; KE += Q/(double)n)
		integral += correctedBetaSpectrum(KE) / n;

	printf("%16.25f\n", integral);
	return integral;
}

void PhaseSpaceCurve() {
	for (double KE = 0; KE < Q; KE += Q/(double)n)
	{
		double P = correctedBetaSpectrum(KE) / n / 1.1356741907175162964449555;
		printf("%f\t%.12f\n", KE, P);
	}
}


int main (int argc, char ** argv) {
	//PhaseSpaceIntegral();
	PhaseSpaceCurve();
}
