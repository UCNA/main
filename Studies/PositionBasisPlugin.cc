#include "PositionBasisPlugin.hh"
#include "SMExcept.hh"

PositionBasisPlugin::PositionBasisPlugin(RunAccumulator* RA, const std::string& nm, unsigned int nr, double r):
AnalyzerPlugin(RA,nm) { /// , sects(nr,r) {
	
	// load sector cutter
    /*
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
	Stringmap ms = SCtoSM(sects);
	myA->qOut.insert("SectorCutter_"+name,ms);
    */
}

/*
std::vector<fgbgPair*> PositionBasisPlugin::allocateSegmentHistograms(TH1& hTemplate, AFPState a, Side s) {
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
*/
