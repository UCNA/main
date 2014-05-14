#ifndef ELECTRONBINDINGENERGY_HH
#define ELECTRONBINDINGENERGY_HH

#include "QFile.hh"
#include "SMExcept.hh"
#include <vector>
#include <map>

/// table of electron binding energies
class BindingEnergyTable {
public:
	/// constructor from Stringmap
	BindingEnergyTable(const Stringmap& m);
	/// get subshell binding energies for given shell
	const std::vector<double>& getShellBinding(unsigned int n) const;
	/// get subshell binding energy
	double getSubshellBinding(unsigned int n, unsigned int m) const;
	/// display summary of binding energies
	void display() const;
	/// get Z
	unsigned int getZ() const { return Z; }
	/// get element name
	std::string getName() const { return nm; }
	
	static const std::string shellnames;
	
protected:
	unsigned int Z;									///< element number
	std::string nm;									///< element name abbrev.
	std::vector< std::vector<double> > eBinding;	///< binding energy by shell and subshell
};

/// catalog of many BindingEnergyTables
class BindingEnergyLibrary {
public:
	/// constructor from QFile containing element tables
	BindingEnergyLibrary(const QFile& Q);
	/// destructor
	~BindingEnergyLibrary();
	/// get BindingEnergyTable for specified element
	const BindingEnergyTable* getBindingTable(unsigned int Z) const;
	/// display contents
	void display() const;
protected:
	std::map<unsigned int,BindingEnergyTable*> tables;	///< binding energy tables by element
};

#endif
