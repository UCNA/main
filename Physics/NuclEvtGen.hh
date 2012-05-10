#ifndef NUCLEVTGEN_HH
#define NUCLEVTGEN_HH 1

#include "ElectronBindingEnergy.hh"
#include "BetaSpectrum.hh"
#include <TF1.h>
#include <vector>
#include <map>
#include <climits>
#include <float.h>
#include <stdio.h>

/// random event selector
class PSelector {
public:
	/// constructor
	PSelector() { cumprob.push_back(0); }
	/// add a probability
	void addProb(double p) { cumprob.push_back(p+cumprob.back()); }
	/// randomly generate according to probability partitions
	unsigned int select() const;
	/// get cumulative probability
	double getCumProb() const { return cumprob.back(); }
	/// get number of items
	unsigned int getN() const { return cumprob.size()-1; }
	/// get probability of numbered item
	double getProb(unsigned int n) const;
	/// scale all probabilities
	void scale(double s);
		
protected:
	std::vector<double> cumprob; //< cumulative probabilites
};

/// Nuclear energy level
class NucLevel {
public:
	/// constructor
	NucLevel(const Stringmap& m);
	
	std::string name;	//< name for this level
	unsigned int A;		//< nucleus A
	unsigned int Z;		//< nucleus Z
	unsigned int n;		//< level number
	double E;			//< energy
	double hl;			//< half-life
	std::string jpi;	//< spin/parity
	double fluxIn;		//< net flux into level
	double fluxOut;		//< net flux out of level
};

enum DecayType {
	D_GAMMA,
	D_ELECTRON,
	D_POSITRON,
	D_NEUTRINO,
	D_NONEVENT
};

/// specification for a decay particle
struct NucDecayEvent {
	double E;		//< particle energy
	DecayType d;	//< particle type
	double t;		//< time of event
};

/// Transition base class
class TransitionBase {
public:
	/// constructor
	TransitionBase(NucLevel& f, NucLevel& t): from(f), to(t), Itotal(0) {}
	/// destructor
	virtual ~TransitionBase() {}
	/// display transition line info
	virtual void display() const { printf("%s->%s (p=%.2g)\n",from.name.c_str(),to.name.c_str(),Itotal); }
	
	/// select transition outcome
	virtual void run(std::vector<NucDecayEvent>&) { }
	
	/// get probability of removing an electron from a given shell
	virtual double getPVacant(unsigned int) const { return 0; }
	/// how many of given electron type were knocked out
	virtual unsigned int nVacant(unsigned int) const { return 0; }
	
	NucLevel& from;	//< level this transition is from
	NucLevel& to;	//< level this transition is to
	double Itotal;	//< total transition intensity
};


/// Gamma transition with possible conversion electrons
class ConversionGamma: public TransitionBase {
public:
	/// constructor
	ConversionGamma(NucLevel& f, NucLevel& t, const Stringmap& m);
	/// select transition outcome
	virtual void run(std::vector<NucDecayEvent>& v);
	/// display transition line info
	virtual void display() const { printf("Gamma %.1f / CE %.1f ",Egamma,averageE()); TransitionBase::display(); }
	
	/// get probability of knocking conversion electron from a given shell
	virtual double getPVacant(unsigned int n) const { return n<shells.getN()-1?shells.getProb(n):0; }
	/// get whether said electron was knocked out
	virtual unsigned int nVacant(unsigned int n) const { return int(n)==shell; }
	/// shell weighted average energy
	double shellAverageE(unsigned int n) const;
	/// line weighted average
	double averageE() const;
	
	double Egamma;		//< gamma energy
	int shell;			//< selected conversion electron shell
	int subshell;		//< selected conversion electron subshell
	double Igamma;		//< total gamma intensity
	BindingEnergyTable const* BET;		//< binding energy table
	
protected:
	PSelector shells;					//< conversion electron shells
	std::vector<PSelector> subshells;	//< subshell choices for each shell
};

class ECapture: public TransitionBase {
public:
	/// constructor
	ECapture(NucLevel& f, NucLevel& t): TransitionBase(f,t), pKCapt(0) {}
	/// select transition outcome
	virtual void run(std::vector<NucDecayEvent>&);
	/// display transition line info
	virtual void display() const { printf("Ecapture "); TransitionBase::display(); }
	/// get probability of removing an electron from a given shell
	virtual double getPVacant(unsigned int n) const { return n==0?pKCapt:0; }
	/// get whether said electron was knocked out
	virtual unsigned int nVacant(unsigned int n) const { return n==0?isKCapt:0; }
	
	double pKCapt;	//< probability of k-shell capture
	bool isKCapt;	//< whether transition was a K capture
};

/// beta decay transitions
class BetaDecayTrans: public TransitionBase {
public:
	/// constructor
	BetaDecayTrans(NucLevel& f, NucLevel& t);
	/// select transition outcome
	virtual void run(std::vector<NucDecayEvent>& v);
	/// display transition line info
	virtual void display() const { printf("Beta(%.1f) ",(W0-1)*m_e); TransitionBase::display(); }
	
	bool positron;		//< whether this is positron decay
	
protected:
	/// evaluate beta spectrum probability
	double evalBeta(double* x, double*);
	/// TF1 for beta spectrum
	TF1 betaTF1;
	
	double R;			//< estimated nucleus radius
	double W0;			//< endpoint, "natural" units
};

/// class for determining probability of Auger electron emission
class AugerManager {
public:
	/// constructor
	AugerManager(const Stringmap& m, unsigned int s);
	/// calculate derived values
	void calculate();
	
	unsigned int shell;	//< which shell to consider augers from
	double Iauger;		//< intensity of Auger electron emissions
	double Ikxr;		//< intensity of k X-Ray emissions
	double IshellOpen;	//< intensity of known shell openings
	double pInitCapt;	//< probability of initial e-capture from shell
	double pAuger;		//< probability of auger given any opening
};

/// Decay system
class NucDecaySystem {
public:
	/// constructor from specification file
	NucDecaySystem(const QFile& Q, const BindingEnergyLibrary& B, double t = DBL_MAX);
	/// set cutoff lifetime for intermediate states
	void setCutoff(double t);
	/// display transitions summary
	void display() const;
	/// get index for named level
	unsigned int levIndex(const std::string& s) const;
	/// generate a chain of decay events starting from level n
	void genDecayChain(std::vector<NucDecayEvent>& v, unsigned int n = UINT_MAX);
protected:
	/// add a transition
	void addTransition(TransitionBase* T);
	
	BindingEnergyLibrary const&  BEL;					//< electron binding energy info
	AugerManager AM;									//< Auger emission probability
	double tcut;
	std::vector<NucLevel> levels;						//< levels, enumerated
	std::map<std::string,unsigned int> levelIndex;		//< energy levels by name
	PSelector lStart;									//< selector for starting level (for breaking up long decays)
	std::vector<PSelector> levelDecays;					//< probabilities for transitions from each level
	std::vector<TransitionBase*> transitions;			//< transitions, enumerated
	std::vector< std::vector<unsigned int> > transIn;	//< transitions into each level
	std::vector< std::vector<unsigned int> > transOut;	//< transitions out of each level
};


#endif
