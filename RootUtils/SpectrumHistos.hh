#ifndef SPECTRUMHISTOS_HH
#define SPECTRUMHISTOS_HH 1

#include <string>
#include <vector>
#include <TPad.h>
#include <TH1.h>
#include "RunManager.hh"
#include "OutputManager.hh"

/// Class for generating background-subtracted histograms
class SpectrumHistos: public OutputManager {
public:
	/// constructor
	SpectrumHistos(std::string nm, std::string bp): OutputManager(nm,bp), fg(NULL), bg(NULL), chName("UnNamed") { }
	
	/// fill background-subtracted spectrum histograms for given points
	void fill(float fgNorm, float bgNorm,
			  const std::vector<float>& fgPoints, const std::vector<float>& bgPoints,
			  const std::string& hName, unsigned int nbins = 400, float xmin=-200, float xmax=4000);
	
	/// fill background-subtracted 2D histograms for given points
	void fill(float fgNorm, float bgNorm,
			  const std::vector<float>& fgPointsX, const std::vector<float>& fgPointsY,
			  const std::vector<float>& bgPointsX, const std::vector<float>& bgPointsY,
			  const std::string& hName,
			  unsigned int nbinsX=200, float xmin=-60, float xmax=60,
			  unsigned int nbinsY=200, float ymin=-60, float ymax=60);
	
	/// setup 1D histos for filling
	void setupHistos(const std::string& hName, unsigned int nbins, float xmin, float xmax);
	/// setup 2D histos for filling
	void setup2D( const std::string& hName,
				 unsigned int nbinsX=200, float xmin=-60, float xmax=60,
				 unsigned int nbinsY=200, float ymin=-60, float ymax=60);

	/// finalize 1D + BG subtraction
	void finishHistos(float fgNorm, float bgNorm);
	/// finalize 2D + BG subtraction
	void finish2D(float fgNorm, float bgNorm);
	
	/// fill 1D point to fg or bg
	void fill(float x, bool isFg);
	/// fill 2D histogram
	void fill(float x, float y, bool isFg);

			  
	/// draw spectrum histograms
	void draw(bool logScale = true, const char* xlabel = NULL, const char* ylabel = NULL, const char* title = NULL);
	
	/// bg-subtracted rate
	float diffRate() const { return diff->Integral()*diff->GetBinWidth(1); }
	
	TH1F* fg;		//< foreground histogram
	TH1F* bg;		//< background histogram
	TH1F* diff;		//< background-subtracted histogram
	
	TH2F* fg2;		//< foreground histogram
	TH2F* bg2;		//< background histogram
	TH2F* diff2;	//< background-subtracted histogram
		
	float fgCount;		//< foreground run counts
	float bgCount;		//< background run counts
	
	std::string chName;	//< current histograms name
	
protected:
	
	static unsigned int nHistos;	//< unique incrementing id number for generating histogram names
	static unsigned int nHistos2;	//< unique incrementing id number for generating histogram names
};

void setLabels(TH1* h, const char* xlabel = NULL, const char* ylabel = NULL, const char* title = NULL);

#endif
