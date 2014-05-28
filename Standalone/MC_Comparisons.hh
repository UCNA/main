#include "G4toPMT.hh"
#include "PenelopeToPMT.hh"
#include "OutputManager.hh"
#include "PathUtils.hh"
#include "GraphicsUtils.hh"
#include <TROOT.h>
#include <TStyle.h>
#include <TLegend.h>
#include <iostream>

	// Declare simulations.
        std::vector<Sim2PMT*> sps;

        // set up lists for histograms/data from each simulation
        std::vector<TH1F*> hEvis[TYPE_II_EVENT+1];
        std::vector<TH1F*> hEvisd[TYPE_II_EVENT+1];
        std::vector<TH1F*> hEquench[TYPE_II_EVENT+1];
        std::vector<TH1F*> hEMWPC[TYPE_II_EVENT+1];
        std::vector<TH1F*> hEMWPCd[TYPE_II_EVENT+1];
        std::vector<TH1F*> hEWires[TYPE_II_EVENT+1];
        std::vector<TH1F*> hEFoils[TYPE_II_EVENT+1];
        std::vector<TH1F*> hEMylarf[TYPE_II_EVENT+1];
        std::vector<TH1F*> hEMylarb[TYPE_II_EVENT+1];
        std::vector<TH1F*> hCTScint;
        std::vector<TH1F*> hCTMWPC;
        std::vector<TH1F*> hCTMWPCo;
        std::vector<TH1F*> hCTFoil;
        std::vector<TH1F*> hPScintX;
        std::vector<TH1F*> hPMWPCX;
        std::vector<TH1F*> hPScintY;
        std::vector<TH1F*> hPMWPCY;
        std::vector<double> t0norm;
        std::vector<double> mwpcnorm;
        std::vector<double> foilnorm;
        std::vector<double> mwpcoutnorm;
        std::vector<double> scintnorm;

	Sim2PMT* SP;
