#include "SpectrumPeak.hh"
#include <TH1F.h>

Float_t SpectrumPeak::energy() const {
	if(type==CD109_PEAK)
		return 75.2;
	if(type==SN_PEAK)
		return 364;
	if(type==SR85_PEAK)
		return 500.1;
	if(type==BI_PEAK_1)
		return 502.6;
	if(type==BI_PEAK_2)
		return 1000.0;
	if(type==M000_PEAK)
		return 502.6;
	if(type==M500_PEAK)
		return 500;
	if(type==CE139_PEAK)
		return 131;
	if(type==IN114_PEAK)
		return 150;
	if(type==LINE_PEAK)
		return energyCenter.x;
	return 0;
}

const char* SpectrumPeak::name() const {
	
	if(type==UNKNOWN_PEAK)
		return "Unknown";
	if(type==SN_PEAK)
		return "Sn";
	if(type==BI_PEAK_1)
		return "Bi_1";
	if(type==BI_PEAK_2)
		return "Bi_2";
	if(type==M000_PEAK)
		return "M000";
	if(type==M500_PEAK)
		return "M500";
	if(type==SR85_PEAK)
		return "Sr85";
	if(type==CD109_PEAK)
		return "Cd109";
	if(type==IN114_PEAK)
		return "In114";
	if(type==CE139_PEAK)
		return "Ce139";
	if(type==REF_CO60_1)
		return "Ref_Co60_1";
	if(type==REF_CO60_2)
		return "Ref_Co60_2";
	if(type==REF_LED)
		return "Ref_LED";
	if(type==TUBE1_LED)
		return "Tube1_LED";
	if(type==TUBE2_LED)
		return "Tube2_LED";
	if(type==TUBE3_LED)
		return "Tube3_LED";
	if(type==TUBE4_LED)
		return "Tube4_LED";
	if(type==LINE_PEAK)
		return "Line_Peak";
	
	return "UnnamedPeak";
}

void SpectrumPeak::print() const { printf("[%s] ADC %f [%f] Energy %f [%f]",name(),center.x,width.x,energyCenter.x,energyWidth.x); }

Stringmap SpectrumPeak::toStringmap() const {
	Stringmap M;
	M.insert("object","SpectrumPeak");
	M.insert("type",type);
	M.insert("center",center.x);
	M.insert("dcenter",center.err);
	M.insert("energy",energyCenter.x);
	M.insert("denergy",energyCenter.err);
	M.insert("width",width.x);
	M.insert("dwidth",width.err);
	M.insert("energywidth",energyWidth.x);
	M.insert("denergywidth",energyWidth.err);
	M.insert("height",h.x);
	M.insert("dheight",h.err);
	M.insert("sID",sID);
	M.insert("x",xpos);
	M.insert("y",ypos);
	M.insert("integral",integral);
	M.insert("nPE",nPE);
	M.insert("simulated",simulated?"yes":"no");
	return M;
}

SpectrumPeak SpectrumPeak::fromStringmap(const Stringmap& M) {
	SpectrumPeak p;
	p.type = PeakType(int(M.getDefault("type",0)));
	p.center.x = M.getDefault("center",0);
	p.center.err = M.getDefault("dcenter",0);
	p.energyCenter.x = M.getDefault("energy",0);
	p.energyCenter.err = M.getDefault("denergy",0);
	p.width.x = M.getDefault("width",0);
	p.width.err = M.getDefault("dwidth",0);
	p.energyWidth.x = M.getDefault("energywidth",0);
	p.energyWidth.err = M.getDefault("denergywidth",0);
	p.h.x = M.getDefault("height",0);
	p.h.err = M.getDefault("dheight",0);
	p.sID = (unsigned int)M.getDefault("sID",0);
	p.xpos = M.getDefault("x",0);
	p.ypos = M.getDefault("y",0);
	return p;
}

void SpectrumPeak::fromGaussian(const TF1* f, unsigned int n) {
	h.x = f->GetParameter(0+3*n);
	h.err = f->GetParError(0+3*n);
	center.x = f->GetParameter(1+3*n);
	center.err = f->GetParError(1+3*n);
	width.x = f->GetParameter(2+3*n);
	width.err = f->GetParError(2+3*n); 
}

SpectrumPeak gausFitter(std::vector<float> dat, Float_t mu0, Float_t sigma0, int nSigma) {
	
	TH1F* foo = new TH1F("foo","Meanfinder",200,mu0-nSigma*sigma0,mu0+nSigma*sigma0);
	for(UInt_t i=0; i<dat.size(); i++)
		foo->Fill(dat[i]);
	
	TF1* g1 = new TF1("g1","gaus",mu0-nSigma*sigma0, mu0 + nSigma*sigma0); 
	foo->Fit(g1,"QR");
	
	SpectrumPeak p;
	p.fromGaussian(g1);

	delete(foo);
	delete(g1);
	
	return p;
}
