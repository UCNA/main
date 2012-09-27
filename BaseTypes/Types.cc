#include "Types.hh"
#include <cassert>

double binterpolate(const TAxis* ax, double binloc) {
	return ax->GetBinCenter(1) + (ax->GetBinCenter(2)-ax->GetBinCenter(1))*binloc;
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
		b[s] *= x;
	return b;
}

