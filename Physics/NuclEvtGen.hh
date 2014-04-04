#ifndef NUCLEVTGEN_HH
#define NUCLEVTGEN_HH

#include "ElectronBindingEnergy.hh"
#include "BetaSpectrum.hh"
#include "FloatErr.hh"
#include "GraphUtils.hh"
#include "TChainScanner.hh"
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
	/// select partition for given input (random if not specified); re-scale input to partition range to pass along to sub-selections
	unsigned int select(double* x = NULL) const;
	/// get cumulative probability
	double getCumProb() const { return cumprob.back(); }
	/// get number of items
	unsigned int getN() const { return cumprob.size()-1; }
	/// get probability of numbered item
	double getProb(unsigned int n) const;
	/// scale all probabilities
	void scale(double s);
		
protected:
	std::vector<double> cumprob;	//< cumulative probabilites
};

/// generate an isotropic random direction, from optional random in [0,1]^2
void randomDirection(double& x, double& y, double& z, double* rnd = NULL);

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

/// primary particle types
/// using PDG numbering scheme, http://pdg.lbl.gov/2012/reviews/rpp2012-rev-monte-carlo-numbering.pdf
enum DecayType {
	D_GAMMA = 22,
	D_ELECTRON = 11,
	D_POSITRON = -11,
	D_NEUTRINO = -12, // anti-nu_e
	D_NONEVENT = 0
};

/// string name of particle types
std::string particleName(DecayType t);
/// decay type from particle name
DecayType particleType(const std::string& s);


/// specification for a decay particle
class NucDecayEvent {
public:
	/// constructor
	NucDecayEvent(): eid(0), E(0), d(D_NONEVENT), t(0), w(1.) {}
	/// randomize momentum direction
	void randp(double* rnd = NULL) { randomDirection(p[0],p[1],p[2],rnd); }
	
	unsigned int eid;	//< event ID number
	double E;			//< particle energy
	double p[3];		//< particle momentum direction
	double x[3];		//< vertex position
	DecayType d;		//< particle type
	double t;			//< time of event
	double w;			//< weighting for event
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
	virtual void run(std::vector<NucDecayEvent>&, double* = NULL) { }
	
	/// return number of continuous degrees of freedom needed to specify transition
	virtual unsigned int getNDF() const { return 2; }
	
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
	virtual void run(std::vector<NucDecayEvent>& v, double* rnd = NULL);
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
	virtual void run(std::vector<NucDecayEvent>&, double* rnd = NULL);
	/// display transition line info
	virtual void display(bool verbose = false) const { printf("Ecapture "); TransitionBase::display(verbose); }
	/// get probability of removing an electron from a given shell
	virtual double getPVacant(unsigned int n) const { return n==0?toAtom->IMissing:0; }
	/// get whether said electron was knocked out
	virtual unsigned int nVacant(unsigned int n) const { return n==0?isKCapt:0; }
	
	/// return number of continuous degrees of freedom needed to specify transition
	virtual unsigned int getNDF() const { return 0; }
	
	bool isKCapt;	//< whether transition was a K capture
};

/// beta decay transitions
class BetaDecayTrans: public TransitionBase {
public:
	/// constructor
	BetaDecayTrans(NucLevel& f, NucLevel& t, bool pstrn = false, unsigned int forbidden = 0);
	/// destructor
	~BetaDecayTrans();
	/// select transition outcome
	virtual void run(std::vector<NucDecayEvent>& v, double* rnd = NULL);
	/// display transition line info
	virtual void display(bool verbose = false) const { printf("Beta(%.1f) ",from.E-to.E); TransitionBase::display(verbose); }
	
	/// return number of continuous degrees of freedom needed to specify transition
	virtual unsigned int getNDF() const { return 3; }
	
	bool positron;		//< whether this is positron decay
	BetaSpectrumGenerator BSG;	//< spectrum shape generator

protected:
	/// evaluate beta spectrum probability
	double evalBeta(double* x, double*);
	TF1 betaTF1;					//< TF1 for beta spectrum shape
	TF1_Quantiles* betaQuantiles;	//< inverse CDF of beta spectrum shape for random point selection
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
	void genDecayChain(std::vector<NucDecayEvent>& v, double* rnd = NULL, unsigned int n = UINT_MAX);
	/// rescale all probabilities
	void scale(double s);
	
	/// return number of degrees of freedom needed to specify decay from given level
	unsigned int getNDF(unsigned int n = UINT_MAX) const;
	
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
	std::vector< std::vector<TransitionBase*> > transIn;	//< transitions into each level
	std::vector< std::vector<TransitionBase*> > transOut;	//< transitions out of each level
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

//------------------------------------------------------------------------------

/// base class for generating event positions
class PositionGenerator {
public:
	/// constructor
	PositionGenerator() {}
	/// destructor
	virtual ~PositionGenerator() {}
	
	/// get number of random DF consumed
	virtual unsigned int getNDF() const { return 3; }
	/// generate vertex position
	virtual void genPos(double* v, double* rnd = NULL) const = 0;
};

/// map unit square onto circle of specified radius
void square2circle(double& x, double& y, double r = 1.0);

/// class for generating positions in a cylinder
class CylPosGen: public PositionGenerator {
public:
	/// constructor
	CylPosGen(double zlength, double radius): dz(zlength), r(radius) {}
	/// generate vertex position
	virtual void genPos(double* v, double* rnd = NULL) const;
	
	double dz;	//< length of cylinder
	double r;	//< radius of cylinder
};

/// uniform cube [0,1]^3 positions, for later transform
class CubePosGen: public PositionGenerator {
public:
	/// constructor
	CubePosGen() {}
	/// generate vertex position
	virtual void genPos(double* v, double* rnd = NULL) const;
};

/// generate fixed event position
class FixedPosGen: public PositionGenerator {
public:
	/// constructor
	FixedPosGen(double x0=0, double y0=0, double z0=0): x(x0), y(y0), z(z0) {}
	/// generate vertex position
	virtual void genPos(double* v, double* = NULL) const { v[0]=x; v[1]=y; v[2]=z; }
	/// get number of random DF consumed
	virtual unsigned int getNDF() const { return 0; }
	
	double x;
	double y;
	double z;
};

//------------------------------------------------------------------------------


/// class for reading events trees
class EventTreeScanner: protected TChainScanner {
public:
	/// constructor
	EventTreeScanner(): TChainScanner("Evts") {}
	/// add a file to the TChain
	virtual int addFile(const std::string& filename);
	/// load next event into vector; return number of primaries
	unsigned int loadEvt(std::vector<NucDecayEvent>& v);

	bool firstpass;	//< whether read is on first pass through data
	
protected:
	/// set tree readpoints
	virtual void setReadpoints();

	NucDecayEvent evt;	//< event readpoint
	unsigned int prevN;	//< previous event number
};

#endif
