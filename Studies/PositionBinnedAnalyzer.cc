#include "PositionBinnedAnalyzer.hh"
#include "SMExcept.hh"

double PositionBinnedAnalyzer::fidRadius = 50.;

PositionBinnedAnalyzer::PositionBinnedAnalyzer(RunAccumulator* RA, const std::string& nm, unsigned int nr):
AnalyzerPlugin(RA,nm), sects(nr,PositionBinnedAnalyzer::fidRadius) {
	
	// load sector cutter
	if(myA->fIn) {
		QFile qOld(myA->inflname+".txt");
		Stringmap sct = qOld.getFirst("SectorCutter_"+name);
		sects = SectorCutter(int(sct.getDefault("nRings",0)),sct.getDefault("radius",0));
		if(!(sects.n && sects.r)) {
			SMExcept e("MissingSectorInfo");
			e.insert("name",name);
			e.insert("file",myA->inflname+".txt");
			throw e;
		}
	}
	
	// save sector cutter
	Stringmap ms;
	ms.insert("nRings",sects.n);
	ms.insert("radius",sects.r);
	ms.insert("nSectors",sects.nSectors());
	myA->qOut.insert("SectorCutter_"+name,ms);
}

std::vector<fgbgPair*> PositionBinnedAnalyzer::allocateSegmentHistograms(TH1& hTemplate, AFPState a, Side s) {
	std::vector<fgbgPair*> segHists;
	std::string hname0 = hTemplate.GetName();
	for(unsigned int m=0; m<=sects.nSectors(); m++) {
		if(m<sects.nSectors())
			hTemplate.SetName((hname0+"_"+itos(m)).c_str());
		else
			hTemplate.SetName(hname0.c_str());
		segHists.push_back(registerFGBGPair(hTemplate,a,s));
	}
	return segHists;
}
