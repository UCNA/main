#include "CathSegCalibrator.hh"
#include <cmath>

double CathSegCalibrator::adjustPos(double x0, double E) const {
	double x = x0;
	E = E>50?(E<1000?E:1000):50;
	for(unsigned int i=0; i<pcoeffs.size(); i++)
		x += (i+1)*pcoeffs[i]->Eval(E)*cos(M_PI*x0);
	return x;
}

CathSegCalibrator::~CathSegCalibrator() {
	for(unsigned int i=0; i<pcoeffs.size(); i++)
		delete(pcoeffs[i]);
}
