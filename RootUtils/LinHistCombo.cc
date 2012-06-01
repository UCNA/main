#include "LinHistCombo.hh"
#include "strutils.hh"

unsigned int LinHistCombo::nFitters = 0;

int LinHistCombo::Fit(TH1* h, double xmin, double xmax, const std::string& fitopt) {
	if(myFit && myFit->GetNpar() != (int)terms.size()) { delete(myFit); myFit = NULL; }
	if(!myFit) myFit = new TF1((std::string("fCombo")+itos(nFitters++)).c_str(),
							   this,&LinHistCombo::Evaluate,xmin,xmax,terms.size());
	myFit->SetRange(xmin,xmax);
	int err = h->Fit(myFit,fitopt.c_str());
	coeffs.clear();
	dcoeffs.clear();
	for(unsigned int i=0; i<terms.size(); i++) {
		coeffs.push_back(myFit->GetParameter(i));
		dcoeffs.push_back(myFit->GetParError(i));
	}
	return err;
}

double LinHistCombo::Evaluate(double* x, double* p) {
	double s = 0;
	for(unsigned int i=0; i<terms.size(); i++) {
		int bn = terms[i]->FindBin(*x);
		if(bn < 1 || bn >= terms[i]->GetNbinsX()-1) continue;
		s += terms[i]->GetBinContent(bn)*p[i];
	}
	return s;
}
