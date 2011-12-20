#ifndef SRASYM_HH
#define ARASYM_HH 1

#include "QFile.hh"
#include "Types.hh"
#include <TF1.h>
#include <TH1F.h>

/// super ratio asymmetry calculator/plotter
class SRAsym {
public:
	/// constructor from combinations of spin state and side
	SRAsym(TH1F* eOff, TH1F* wOff, TH1F* eOn, TH1F* wOn, bool bonehead = false);
	
	/// fit and plot asymmetry
	void fitAsym(float emin, float emax);	
	/// output writable stringmap
	Stringmap toStringmap() const;
	/// draw asymmetry spectra
	void drawSpectra(int lcolor = 0,bool same = false,std::string outname = "");

	TH1F* hAsym;	//< asymmetry histogram
	TH1F* hAnorm;	//< asymmetry / beta
	TH1F* hSum;		//< sum rate histogram
	bool boneit;	//< whether to use bonehead asymmetry
		
protected:
	
	/// calculate super-ratio asymmetry from given histograms
	void calcAsym();
	/// normalize asymetry by beta
	void normAsym();
	/// calculate sum histogram
	void calcSum();
	/// calculate bonehead asymmetry
	void calcBonehead();
	
	TF1* afit;		//< fit to asymmetry
	TF1* nafit;		//< fit to normalized asymmetry
	TH1F* hOff[2];	//< flipper off spectra
	TH1F* hOn[2];	//< flipper on spectra
};

#endif
