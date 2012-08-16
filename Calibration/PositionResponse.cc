#include "PositionResponse.hh"
#include <stdlib.h>

PositioningInterpolator::PositioningInterpolator(const PosmapInfo& PMI):
S(PMI.nRings,PMI.radius), sRadial(BC_DERIVCLAMP_ZERO), L(&sRadial, PMI.radius*(1.0+1.0/(2*PMI.nRings-1.0)), 0) {
	
	assert(PMI.signal.size() == S.nSectors());
	
	// set up sequences/interpolators
	for(unsigned int n=0; n<PMI.nRings; n++) {
		phiSeqs.push_back(new DoubleSequence(BC_CYCLIC));
		phiInterps.push_back(new CubiTerpolator(phiSeqs.back(),2.0*PI,PI/float(S.getNDivs(n))));
		//phiInterps.push_back(new Interpolator(phiSeqs.back(),2.0*PI,PI/float(S.getNDivs(n))));
		sRadial.addPoint(phiInterps.back());
	}
	
	// load data points
	for(unsigned int i=0; i<S.nSectors(); i++)
		phiSeqs[S.getRing(i)]->addPoint(PMI.signal[i]/PMI.norm[i]);
	
}

PositioningInterpolator::~PositioningInterpolator() {
	for(unsigned int i=0; i<phiSeqs.size(); i++)
		delete(phiSeqs[i]);
	for(unsigned int i=0; i<phiInterps.size(); i++)
		delete(phiInterps[i]);
}

double PositioningInterpolator::eval(double x, double y) {
	double xp[2] = {sqrt(x*x+y*y),atan2(y,x)};
	return L.eval(xp);
}

void PositioningCorrector::initPIs(std::vector<PosmapInfo>& indat) {
	for(std::vector<PosmapInfo>::iterator it = indat.begin(); it != indat.end(); it++) {
		assert(it->s==EAST || it->s==WEST);
		while(tubes[it->s].size()<=it->t)
			tubes[it->s].push_back(NULL);
		assert(!tubes[it->s][it->t]);
		tubes[it->s][it->t] = new PositioningInterpolator(*it);
	}
}

void PositioningCorrector::setNormCenter(double c) {
	for(Side s = EAST; s <= WEST; ++s) {
		neta[s].clear();
		for(unsigned int t=0; t<tubes[s].size(); t++)
			neta[s].push_back(eval(s,t,0.0,0.0,false)/c);
	}
}

void PositioningCorrector::setNormAvg(double c) {
	for(Side s = EAST; s <= WEST; ++s) {
		neta[s].clear();
		for(unsigned int t=0; t<tubes[s].size(); t++) {
			SectorCutter& tsects = getSectors(s,t);
			float x,y;
			double avg = 0;
			for(unsigned int m=0; m<tsects.nSectors(); m++) {
				tsects.sectorCenter(m,x,y);
				avg += eval(s,t,x,y,false)*tsects.sectorArea(m);
			}
			avg /= tsects.totalArea();
			neta[s].push_back(avg/c);
		}
	}
}

PositioningCorrector::PositioningCorrector(std::vector<PosmapInfo>& indat) {
	initPIs(indat);
	setNormCenter();
}

PositioningCorrector::PositioningCorrector(QFile& qin) {
	
	// init PosmapInfo vector
	int nRings = atoi(qin.getDefault("SectorCutter","nRings","0").c_str());
	float radius = atof(qin.getDefault("SectorCutter","radius","0").c_str());
	assert(nRings && radius);
	SectorCutter S(nRings,radius);
	std::vector<PosmapInfo> pinf;
	pinf.resize(2*nBetaTubes);
	for(Side s = EAST; s <= WEST; ++s) {
		for(unsigned int t=0; t<nBetaTubes; t++) {
			PosmapInfo pmi;
			pmi.s = s;
			pmi.t = t;
			pmi.nRings = nRings;
			pmi.radius = radius;
			pmi.signal.resize(S.nSectors());
			pmi.norm.resize(S.nSectors());
			pinf.push_back(pmi);
		}
	}
	// read in data points
	std::vector<Stringmap> dpts =  qin.retrieve("PosmapPoint");
	for(std::vector<Stringmap>::iterator it = dpts.begin(); it != dpts.end(); it++) {
		std::string ss = it->getDefault("side","N");
		assert(ss=="E" || ss=="W");
		Side s = ss=="E"?EAST:WEST;
		unsigned int t = (unsigned int)it->getDefault("tube",nBetaTubes);
		unsigned int n = (unsigned int)it->getDefault("sector",S.nSectors());
		float z = it->getDefault("light",1.0);
		float z0 = it->getDefault("energy",1.0);
		assert(t<nBetaTubes && n<S.nSectors());
		pinf[nBetaTubes*s+t].signal[n] = z;
		pinf[nBetaTubes*s+t].norm[n] = z0;
	}
	
	initPIs(pinf);
	setNormCenter();
}

PositioningCorrector::~PositioningCorrector() {
	for(Side s = EAST; s <= WEST; ++s)
		for(unsigned int i=0; i<tubes[s].size(); i++)
			if(tubes[s][i]) delete tubes[s][i];
}

double PositioningCorrector::eval(Side s, unsigned int t, double x, double y, bool normalize) const {
	if(s>WEST || t>=tubes[s].size() || !tubes[s][t])
		return 0;
	if(normalize)
		return tubes[s][t]->eval(x,y)/neta[s][t];
	return tubes[s][t]->eval(x,y);
}
