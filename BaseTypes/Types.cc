#include "Types.hh"
#include <cassert>

Double_t poiscdf(const Double_t *x, const Double_t *par) {
	float w0 = par[2]*par[2];
	float x0 = (x[0]-par[0])/par[1]*par[2] + w0 + 1.0/3.0 - 0.015/w0;
	if(x0<0)
		return 0;
	return (1.0-TMath::Gamma(x0,w0))*par[3];
}

Double_t fancyfish(const Double_t *x, const Double_t *par) {
	float adc50 = par[0];						// ADC at 50% efficiency
	float w = par[1];							// width of efficiency transition
	float n50 = par[2]*adc50/w;					// "number of photoelectrons" at 50% efficiency
	float n0 = n50+sqrt(n50)*(x[0]-adc50)/w;	// "number of photoelectrons" at current point
	if(n0 < 0)
		return 0;
	return TMath::Gamma(n50+1/3-0.015/n50,n0)*par[3];
}

double binterpolate(const TAxis* ax, double binloc) {
	return ax->GetBinCenter(1) + (ax->GetBinCenter(2)-ax->GetBinCenter(1))*binloc;
}

void zero(TH1* h) {
	assert(h);
	h->Scale(0.0);
	h->SetEntries(0);
}

unsigned int totalBins(const TH1* h) {
	assert(h);
	Int_t nbinsx = h->GetNbinsX()+2;
	Int_t nbinsy = h->GetDimension()<2?1:h->GetNbinsY()+2;
	Int_t nbinsz = h->GetDimension()<3?1:h->GetNbinsZ()+2;
	return nbinsx*nbinsy*nbinsz;
}


BlindTime::BlindTime(const Stringmap& m) {
	for(Side s = EAST; s != BADSIDE; ++s)
		t[s]=m.getDefault(sideWords(s),0.);
}

Stringmap BlindTime::toStringmap() const {
	Stringmap m;
	for(Side s = EAST; s != BADSIDE; ++s)
		m.insert(sideWords(s),t[s]);
	return m;
}

BlindTime operator+(const BlindTime& a, const BlindTime& b) {
	BlindTime c = a;
	c += b;
	return c;
}

BlindTime operator-(const BlindTime& a, const BlindTime& b) {
	BlindTime c = a;
	c -= b;
	return c;
}

BlindTime operator*(double x, const BlindTime& a) {
	BlindTime b = a;
	for(Side s = EAST; s != BADSIDE; ++s)
		b.t[s] *= x;
	return b;
}

