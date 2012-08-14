#ifndef NUCLEVTGEN_HH
#define NUCLEVTGEN_HH 1

#include "ElectronBindingEnergy.hh"
#include "BetaSpectrum.hh"
#include "FloatErr.hh"
#include <TF1.h>
#include <vector>
#include <map>
#include <set>
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
	/// print info
	void display(bool verbose = false) const;
	/// scale probabilities
	void scale(double s) { fluxIn *= s; fluxOut *= s; }
	
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

/// Atom/electron information
class DecayAtom {
public:
	/// constructor
	DecayAtom(BindingEnergyTable const* B);
	/// load Auger data from Stringmap
	void load(const Stringmap& m);
	/// generate Auger K probabilistically
	void genAuger(std::vector<NucDecayEvent>& v);
	/// display info
	void display(bool verbose = false) const;
	
	BindingEnergyTable const* BET;	//< binding energy table
	double Eauger;		//< Auger K energy
	double Iauger;		//< intensity of Auger electron emissions
	double Ikxr;		//< intensity of k X-Ray emissions
	double ICEK;		//< intensity of CE K events
	double IMissing;	//< intensity of un-accounted for Augers (from initial capture)
	double pAuger;		//< probability of auger given any opening
};

/// Transition base class
class TransitionBase {
public:
	/// constructor
	TransitionBase(NucLevel& f, NucLevel& t): from(f), to(t), Itotal(0) {}
	/// destructor
	virtual ~TransitionBase() {}
	/// display transition line info
	virtual void display(bool verbose = false) const;
	
	/// select transition outcome
	virtual void run(std::vector<NucDecayEvent>&) { }
	
	/// scale probability
	virtual void scale(double s) { Itotal *= s; }
	
	/// get probability of removing an electron from a given shell
	virtual double getPVacant(unsigned int) const { return 0; }
	/// how many of given electron type were knocked out
	virtual unsigned int nVacant(unsigned int) const { return 0; }
	
	DecayAtom* toAtom;	//< final state atom info
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
	virtual void display(bool verbose = false) const;
	/// get total conversion efficiency
	double getConversionEffic() const;
	/// get probability of knocking conversion electron from a given shell
	virtual double getPVacant(unsigned int n) const { return n<shells.getN()-1?shells.getProb(n):0; }
	/// get whether said electron was knocked out
	virtual unsigned int nVacant(unsigned int n) const { return int(n)==shell; }
	/// shell weighted average energy
	double shellAverageE(unsigned int n) const;
	/// line weighted average
	float_err averageE() const;
	/// scale probability
	virtual void scale(double s);
	
	double Egamma;		//< gamma energy
	int shell;			//< selected conversion electron shell
	int subshell;		//< selected conversion electron subshell
	double Igamma;		//< total gamma intensity
	
protected:
	PSelector shells;					//< conversion electron shells
	std::vector<float> shellUncert;		//< uncertainty on shell selection probability
	std::vector<PSelector> subshells;	//< subshell choices for each shell
};

/// electron capture transitions
class ECapture: public TransitionBase {
public:
	/// constructor
	ECapture(NucLevel& f, NucLevel& t): TransitionBase(f,t) {}
	/// select transition outcome
	virtual void run(std::vector<NucDecayEvent>&);
	/// display transition line info
	virtual void display(bool verbose = false) const { printf("Ecapture "); TransitionBase::display(verbose); }
	/// get probability of removing an electron from a given shell
	virtual double getPVacant(unsigned int n) const { return n==0?toAtom->IMissing:0; }
	/// get whether said electron was knocked out
	virtual unsigned int nVacant(unsigned int n) const { return n==0?isKCapt:0; }
	
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
	virtual void display(bool verbose = false) const { printf("Beta(%.1f) ",from.E-to.E); TransitionBase::display(verbose); }
	
	bool positron;		//< whether this is positron decay
	
protected:
	/// evaluate beta spectrum probability
	double evalBeta(double* x, double*);
	/// TF1 for beta spectrum
	TF1 betaTF1;
};

/// Decay system
class NucDecaySystem {
public:
	/// constructor from specification file
	NucDecaySystem(const QFile& Q, const BindingEnergyLibrary& B, double t = DBL_MAX);
	/// destructor
	~NucDecaySystem();
	/// set cutoff lifetime for intermediate states
	void setCutoff(double t);
	/// display transitions summary
	void display(bool verbose = false) const;
	/// display list of levels
	void displayLevels(bool verbose = false) const;
	/// display list of transitions
	void displayTransitions(bool verbose = false) const;
	/// display list of atoms
	void displayAtoms(bool verbose = false) const;
	/// generate a chain of decay events starting from level n
	void genDecayChain(std::vector<NucDecayEvent>& v, unsigned int n = UINT_MAX);
	/// rescale all probabilities
	void scale(double s);
	
	/// LaTeX name for generator
	std::string fancyname;
	
protected:
	/// get index for named level
	unsigned int levIndex(const std::string& s) const;
	/// get atom info for given Z
	DecayAtom* getAtom(unsigned int Z);
	/// add a transition
	void addTransition(TransitionBase* T);
	
	BindingEnergyLibrary const&  BEL;					//< electron binding energy info
	double tcut;
	std::vector<NucLevel> levels;						//< levels, enumerated
	std::map<std::string,unsigned int> levelIndex;		//< energy levels by name
	PSelector lStart;									//< selector for starting level (for breaking up long decays)
	std::vector<PSelector> levelDecays;					//< probabilities for transitions from each level
	std::map<unsigned int, DecayAtom*> atoms;			//< atom information
	std::vector<TransitionBase*> transitions;			//< transitions, enumerated
	std::vector< std::vector<unsigned int> > transIn;	//< transitions into each level
	std::vector< std::vector<unsigned int> > transOut;	//< transitions out of each level
};

/// manager for loading decay event generators
class NucDecayLibrary {
public:
	/// constructor
	NucDecayLibrary(const std::string& datp, double t = DBL_MAX);
	/// destructor
	~NucDecayLibrary();
	/// check if generator is available
	bool hasGenerator(const std::string& nm);
	/// get decay generator by name
	NucDecaySystem& getGenerator(const std::string& nm);
	
	std::string datpath;	//< path to data folder
	double tcut;			//< event generator default cutoff time
	BindingEnergyLibrary  BEL;	//< electron binding energy info

protected:
	std::map<std::string,NucDecaySystem*> NDs;	//< loaded decay systems
	std::set<std::string> cantdothis;			//< list of decay systems that can't be loaded
};

/// class for throwing from large list of gammas
class GammaForest {
public:
	/// constructor
	GammaForest(const std::string& fname, double E2keV = 1000);
	/// get total cross section
	double getCrossSection() const { return gammaProb.getCumProb(); }
	/// generate cluster of gamma decays
	void genDecays(std::vector<NucDecayEvent>& v, double n = 1.0);
protected:
	std::vector<double> gammaE;	//< gamma energies
	PSelector gammaProb;		//< gamma probabilities selector
};

#endif
