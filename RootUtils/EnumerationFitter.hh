#ifndef ENUMERATIONFITTER_HH
#define ENUMERATIONFITTER_HH 1

#include <vector>
#include <string>
#include <TF1.h>
#include <TGraphErrors.h>

class EnumerationFitter {
public:
	/// constructor
	EnumerationFitter(): fitter(NULL) {}
	/// destructor
	~EnumerationFitter() { if(!fitter) delete fitter; }
	/// add a fit terms set
	void addTerm(const std::vector<double>& t);
	/// fit evaluation from sum of terms
	double Evaluate(double *x, double *p);
	/// get number of fit parameters
	unsigned int getNParams() const { return fterms.size(); }
	/// get fitter
	TF1* getFitter();
	/// load fittable data and terms from a file
	TGraphErrors* loadFitFile(const std::string& fname);
	
protected:
	
	std::vector< std::vector<double> > fterms;	//< fit term sets
	TF1* fitter;								//< fitter based on these terms
};


#endif
