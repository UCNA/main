#include "EnumerationFitter.hh"
#include "SMExcept.hh"
#include "strutils.hh"
#include "PathUtils.hh"
#include <cassert>



double EnumerationFitter::Evaluate(double *x, double *p) {
	assert(x);
	int i = (int)(*x);
	double s = 0;
	for(unsigned int n=0; n<fterms.size(); n++) {
		if(i>=0 && i<(int)fterms[n].size());
			s += fterms[n][i]*p[n];
	}
	return s;
}

void EnumerationFitter::addTerm(const std::vector<double>& t) {
	if(fitter) {
		delete fitter;
		fitter = NULL;
	}
	fterms.push_back(t);
}

TF1* EnumerationFitter::getFitter() {
	if(!fitter)
		fitter = new TF1("fEnumFit",this,&EnumerationFitter::Evaluate,0,1,fterms.size());
	return fitter;
}

TGraphErrors* EnumerationFitter::loadFitFile(const std::string& fname) {
	if(!fileExists(fname)) {
		SMExcept e("fileUnreadable");
		e.insert("filename",fname);
		throw(e);
	}
	
	std::ifstream fin(fname.c_str());
	std::string s;
	std::vector<double> datenum;
	std::vector<double> dat;
	std::vector<double> daterr;
	fterms.clear();
	if(fitter) {
		delete fitter;
		fitter = NULL;
	}
	
	printf("Loading data from '%s'...\n",fname.c_str());
	while (fin.good()) {
		std::getline(fin,s);
		s = strip(s);
		if(!s.size() || s[0]=='#')
			continue;
		std::vector<double> v = sToDoubles(s," ,\t");
		if(v.size() < 2) continue;
		datenum.push_back(0.5+dat.size());
		dat.push_back(v[0]);
		daterr.push_back(v[1]);
		for(unsigned int i=2; i<v.size(); i++) {
			while(fterms.size()<i-1)
				fterms.push_back(std::vector<double>());
			fterms[i-2].push_back(v[i]);
		}
	}
	fin.close();
	
	printf("Loaded %i fit points and %i parameters\n",(int)dat.size(),(int)fterms.size());
	
	return new TGraphErrors(dat.size(),&datenum[0],&dat[0],NULL,&daterr[0]);
}
