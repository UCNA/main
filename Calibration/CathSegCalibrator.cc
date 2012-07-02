#include "CathSegCalibrator.hh"
#include <cmath>

double CathSegCalibrator::adjustPos(double x0, double E) const {
	double x = x0;
	E = E>0?(E<1000?E:1000):0;
	for(unsigned int i=0; i<pcoeffs.size(); i++) {
		unsigned int n = 1+(i/2);
		if(i%2)
			x += pcoeffs[i]->Eval(E)*sin(2*M_PI*n*x0)/(2*M_PI*n);
		else
			x += -pcoeffs[i]->Eval(E)*(cos(2*M_PI*n*x0)+1)/(2*M_PI*n);
	}
	return x;
}

CathSegCalibrator::~CathSegCalibrator() {
	for(unsigned int i=0; i<pcoeffs.size(); i++)
		delete(pcoeffs[i]);
}
