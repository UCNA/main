#include "NuclEvtGen.hh"
#include "SMExcept.hh"
#include "strutils.hh"
#include <cassert>
#include <stdlib.h>
#include <float.h>
#include <algorithm>

double rand_dbl() { return double(rand())/double(RAND_MAX); }

unsigned int PSelector::select() const {
	std::vector<double>::const_iterator itsel = std::upper_bound(cumprob.begin(),cumprob.end(),cumprob.back()*rand_dbl());
	unsigned int selected = (unsigned int)(itsel-cumprob.begin()-1);
	assert(selected<cumprob.size()-1);
	return selected;
}

void PSelector::scale(double s) {
	for(std::vector<double>::iterator it = cumprob.begin(); it != cumprob.end(); it++)
		(*it) *= s;
}

double PSelector::getProb(unsigned int n) const {
	assert(n<cumprob.size()-1);
	return (cumprob[n+1]-cumprob[n])/cumprob.back();
}

//-----------------------------------------

NucLevel::NucLevel(const Stringmap& m): fluxIn(0), fluxOut(0) {
	name = m.getDefault("nm","0.0.0");
	std::vector<std::string> v = split(name,".");
	assert(v.size()==3);
	A = atoi(v[0].c_str());
	Z = atoi(v[1].c_str());
	n = atoi(v[2].c_str());
	E = m.getDefault("E",0);
	hl = m.getDefault("hl",0);
	if(hl < 0) hl = DBL_MAX;
	jpi = m.getDefault("jpi","");
}

//-----------------------------------------

ConversionGamma::ConversionGamma(NucLevel& f, NucLevel& t, const Stringmap& m): TransitionBase(f,t), BET(NULL) {
	Egamma = from.E - to.E;
	Igamma = m.getDefault("Igamma",0.0)/100.0;
	// load conversion electron and subshell probabilities
	unsigned int nshells = BindingEnergyTable::shellnames.size();
	for(unsigned int i=0; i<nshells; i++) {
		std::vector<double> v = sToDoubles(m.getDefault(std::string("CE_")+ctos(BindingEnergyTable::shellnames[i]),""));
		if(!v.size()) break;
		if(v.size()==1) v.push_back(1.);
		shells.addProb(v[0]);
		subshells.push_back(PSelector());
		for(unsigned int s=1; s<v.size(); s++)
			subshells.back().addProb(v[s]);
	}
	// assign remaining probability for gamma
	Itotal = Igamma*(1.0+shells.getCumProb());
	shells.addProb(1.0-shells.getCumProb());
}

void ConversionGamma::run(std::vector<NucDecayEvent>& v) {
	assert(BET);
	shell = (int)shells.select();
	if(shell < (int)subshells.size())
		subshell = (int)subshells[shell].select();
	else
		shell = subshell = -1;
	NucDecayEvent evt;
	evt.E = Egamma;
	if(shell<0) {
		evt.d = D_GAMMA;
	} else {
		evt.d = D_ELECTRON;
		evt.E -= BET->getSubshellBinding(shell,subshell);
	}	
	v.push_back(evt);
}

double ConversionGamma::shellAverageE(unsigned int n) const {
	assert(n<subshells.size());
	double e = 0;
	double w = 0;
	for(unsigned int i=0; i<subshells[n].getN(); i++) {
		double p = subshells[n].getProb(i);
		w += p;
		e += (Egamma-BET->getSubshellBinding(n,i))*p;
	}
	return e/w;
}

double ConversionGamma::averageE() const {
	double e = 0;
	double w = 0;
	for(unsigned int n=0; n<subshells.size(); n++) {
		double p = shells.getProb(n);
		e += shellAverageE(n)*p;
		w += p;
	}
	return e/w;
}

//-----------------------------------------

BetaDecayTrans::BetaDecayTrans(NucLevel& f, NucLevel& t):
TransitionBase(f,t), betaTF1((f.name+"-"+t.name+"_Beta").c_str(),this,&BetaDecayTrans::evalBeta,0,1,0) {
	betaTF1.SetNpx(1000);
	betaTF1.SetRange(0,from.E-to.E-m_e);
	W0 = (from.E-to.E)/m_e;
	R = pow(to.A,1./3.)*neutron_R0;
	positron = false;
}

void BetaDecayTrans::run(std::vector<NucDecayEvent>& v) {
	NucDecayEvent evt;
	evt.d = D_ELECTRON;
	evt.E = betaTF1.GetRandom();
	v.push_back(evt);
}

double BetaDecayTrans::evalBeta(double* x, double*) {
	double W = (x[0]+m_e)/m_e;
	return 1<W && W<W0 ? plainPhaseSpace(W,W0)*WilkinsonF0(to.Z*(positron?-1:1),W,R)*WilkinsonL0(to.Z,W,R)*(1.+Wilkinson_g(W,W0)) : 0;
}


//-----------------------------------------


AugerManager::AugerManager(const Stringmap& m, unsigned int s): shell(s), Iauger(0), Ikxr(0), IshellOpen(0) {
	for(std::multimap< std::string, std::string >::const_iterator it = m.dat.begin(); it != m.dat.end(); it++) {
		assert(it->first.size());
		if(it->first[0]=='a') Iauger += atof(it->second.c_str())/100.0;
		else if(it->first[0]=='k') Ikxr += atof(it->second.c_str())/100.0;
	}
}

void AugerManager::calculate() {
	pInitCapt = Iauger+Ikxr-IshellOpen;
	pAuger = Iauger/(Iauger+Ikxr);
	if(!Iauger) pInitCapt = pAuger = 0;
}

//-----------------------------------------

bool sortLevels(const NucLevel& a, const NucLevel& b) { return (a.E < b.E); }

NucDecaySystem::NucDecaySystem(const QFile& Q, const BindingEnergyLibrary& B): BEL(B), AM(Q.getFirst("AugerK"),0) {
	// load levels data
	std::vector<Stringmap> levs = Q.retrieve("level");
	for(std::vector<Stringmap>::iterator it = levs.begin(); it != levs.end(); it++) {
		levels.push_back(NucLevel(*it));
		transIn.push_back(std::vector<unsigned int>());
		transOut.push_back(std::vector<unsigned int>());
		levelDecays.push_back(PSelector());
	}
	std::sort(levels.begin(),levels.end(),sortLevels);
	for(std::vector<NucLevel>::iterator it = levels.begin(); it != levels.end(); it++) {
		assert(levelIndex.find(it->name) == levelIndex.end());
		it->n = it-levels.begin();
		levelIndex.insert(std::make_pair(it->name,it->n));
	}
	
	// set up internal conversions
	printf("Loading internal converions...\n");
	std::vector<Stringmap> gammatrans = Q.retrieve("gamma");
	for(std::vector<Stringmap>::iterator it = gammatrans.begin(); it != gammatrans.end(); it++) {
		ConversionGamma* CG = new ConversionGamma(levels[levIndex(it->getDefault("from",""))],levels[levIndex(it->getDefault("to",""))],*it);
		CG->BET = &BEL.getBindingTable(CG->to.Z);
		CG->display();
		addTransition(CG);
	}
	
	printf("Loading beta decays...\n");
	std::vector<Stringmap> betatrans = Q.retrieve("beta");
	for(std::vector<Stringmap>::iterator it = betatrans.begin(); it != betatrans.end(); it++) {
		BetaDecayTrans* BD = new BetaDecayTrans(levels[levIndex(it->getDefault("from",""))],levels[levIndex(it->getDefault("to",""))]);
		BD->Itotal = it->getDefault("I",0);
		BD->positron = it->getDefault("positron",0);
		BD->display();
		addTransition(BD);
	}
	
	// set up electron captures
	printf("Loading electron captures...\n");
	std::vector<Stringmap> ecapts = Q.retrieve("ecapt");
	for(std::vector<Stringmap>::iterator it = ecapts.begin(); it != ecapts.end(); it++) {
		NucLevel& Lorig = levels[levIndex(it->getDefault("from",""))];
		std::string to = it->getDefault("to","AUTO");
		if(to == "AUTO") {
			for(std::vector<NucLevel>::iterator Ldest = levels.begin(); Ldest != levels.end(); Ldest++) {
				if(Ldest->A == Lorig.A && Ldest->Z+1 == Lorig.Z && Ldest->E < Lorig.E) {
					double missingFlux = Ldest->fluxOut - Ldest->fluxIn;
					if(missingFlux <= 0) continue;
					ECapture* EC = new ECapture(Lorig,*Ldest);
					EC->Itotal = missingFlux;
					EC->display();
					addTransition(EC);
				}
			}
		} else {
			NucLevel& Ldest = levels[levIndex(to)];
			assert(Ldest.A == Lorig.A && Ldest.Z+1 == Lorig.Z && Ldest.E < Lorig.E);
			ECapture* EC = new ECapture(Lorig,Ldest);
			EC->Itotal = it->getDefault("I",0.);
			EC->display();
			addTransition(EC);
		}
	}
	
	printf("Calculating for Augers...\n");	
	AM.calculate();
	setCutoff(DBL_MAX);
}

void NucDecaySystem::addTransition(TransitionBase* T) {
	AM.IshellOpen += T->getPVacant(0)*T->Itotal;
	transIn[T->to.n].push_back(transitions.size());
	transOut[T->from.n].push_back(transitions.size());
	levelDecays[T->from.n].addProb(T->Itotal);
	T->from.fluxOut += T->Itotal;
	T->to.fluxIn += T->Itotal;
	transitions.push_back(T);
}

void NucDecaySystem::setCutoff(double t) {
	tcut = t;
	lStart = PSelector();
	for(unsigned int n=0; n<levels.size(); n++) {
		levelDecays[n] = PSelector();
		for(unsigned int t = 0; t < transOut[n].size(); t++)
			levelDecays[n].addProb(transitions[transOut[n][t]]->Itotal);
		
		double pStart = (n+1==levels.size());
		if(levels[n].hl > tcut && transOut[n].size())
			for(unsigned int t = 0; t < transIn[n].size(); t++)
				pStart += transitions[transIn[n][t]]->Itotal;
		lStart.addProb(pStart);
	}
}

void NucDecaySystem::display() const {
	printf("---- Nuclear Level System ----\n");
	unsigned int lnum = 0;
	for(std::vector<NucLevel>::const_iterator it = levels.begin(); it != levels.end(); it++) {
		unsigned int n = it-levels.begin();
		printf("%s %.1f\t[%.3g; %.3g]\t%.3g <<",it->name.c_str(),it->E,levelDecays[n].getCumProb(),lStart.getProb(n),it->hl);
		for(std::vector<unsigned int>::const_iterator it2 = transIn[n].begin(); it2 != transIn[n].end(); it2++)
			printf(" %s(%.2g)",transitions[*it2]->from.name.c_str(),transitions[*it2]->Itotal);
		printf(" >>");
		for(std::vector<unsigned int>::const_iterator it2 = transOut[n].begin(); it2 != transOut[n].end(); it2++)
			printf(" %s(%.2g)",transitions[*it2]->to.name.c_str(),transitions[*it2]->Itotal);
		printf("\n");
		lnum++;
	}
	printf("-- %.2g initial K-shell capture, %.2g Auger K probability --\n",AM.pInitCapt,AM.pAuger);
}

unsigned int NucDecaySystem::levIndex(const std::string& s) const {
	std::map<std::string,unsigned int>::const_iterator n = levelIndex.find(s);
	if(n==levelIndex.end()) {
		SMExcept e("UnknownLevel");
		e.insert("name",s);
		throw(e);
	}
	return n->second;
}

void NucDecaySystem::genDecayChain(std::vector<NucDecayEvent>& v, unsigned int n) {
	bool init = n>=levels.size();
	if(init) n = lStart.select();
	if(!levels[n].fluxOut || (!init && levels[n].hl > tcut)) return;
	
	TransitionBase* T = transitions[transOut[n][levelDecays[n].select()]];
	T->run(v);
	
	unsigned int nAugerK = init*(rand_dbl()<AM.pInitCapt)+T->nVacant(0);
	while(nAugerK--) {
		if(rand_dbl() > AM.pAuger) continue;
		const BindingEnergyTable& BT = BEL.getBindingTable(T->to.Z);
		NucDecayEvent evt;
		evt.d = D_ELECTRON;
		evt.E = BT.getSubshellBinding(0,0) - BT.getSubshellBinding(1,0) - BT.getSubshellBinding(1,1);
		v.push_back(evt);
	}
	
	genDecayChain(v,T->to.n);
}


