#ifndef KDATA_HH
#define KDATA_HH 1

#include <string>
#include <vector>
#include <map>
#include <cassert>
#include "strutils.hh"

/// base class for recursive string-tagged data
class RData {
public:
	/// constructor
	RData() {}
	/// destructor
	virtual ~RData() {}
	
	// implement these in subclasses
	
	/// get list of all keys
	virtual std::vector<std::string> getKeys() const { return std::vector<std::string>(); }
	/// get first key
	virtual std::string getFirstKey(std::string dflt = "") const { return dflt; }
	/// get subdata for key
	virtual std::vector<RData*> getSubdata(const std::string& key) { return std::vector<RData*>(); } 
	/// get first subdata for a key
	virtual RData* getFirst(const std::string& key) { return RData::NullRData; }
	/// insert a key --- base class is read-only
	virtual RData* insert(const std::string& key) { assert(false); return NULL; }
	/// return number of keys
	virtual unsigned int size() { return 0; }

	// these should already work for subclasses
	
	/// get first key as float
	virtual float getFirstD(double dflt) const;
	/// get first subdata, creating missing key if necessary
	virtual RData* getForced(const std::string& key);
	/// insert a numerical key
	virtual RData* insert(double x) { return insert(dtos(x)); }
	/// generate a string representation of data
	virtual std::string toString(bool displaymode = false, const std::string& pfx = "");
	/// check whether this is the null data object
	bool isNull() const { return this==NullRData; }
	/// recursively search for matching subdata, returning all
	std::vector<RData*> getSubdata(const std::vector<std::string>& keys, unsigned int kdepth = 0);
	/// recursively search for matching subdata, returning first found
	RData* getFirst(const std::vector<std::string>& keys, unsigned int kdepth = 0);
	/// get first by "path" style string
	RData* getFirstByPath(const std::string& path, const std::string& pathsep = "/");
	/// get all matching "path" style string
	std::vector<RData*> getSubdataByPath(const std::string& path, const std::string& pathsep = "/");
protected:
	static RData* NullRData;	//< "null" return result for keys not found
};

/// memory-resident RData
class RDataMem: public RData {
public:
	/// constructor
	RDataMem(): RData() {}
	/// copy constructor from RData
	RDataMem(RData* KD);
	/// destructor
	~RDataMem();
	
	/// get list of all keys
	virtual std::vector<std::string> getKeys() const;
	/// get first key
	virtual std::string getFirstKey(std::string dflt = "") const;
	/// get subdata for key
	virtual std::vector<RData*> getSubdata(const std::string& key);
	/// get first subdata for a key
	virtual RData* getFirst(const std::string& key);
	/// insert a key
	virtual RData* insert(const std::string& key);
	/// return number of keys
	virtual unsigned int size() { return dat.size(); }
	
	/// insert a copy from a RData for key
	RData* insert(const std::string& key, RData* KD);
	
protected:
	std::multimap< std::string, RDataMem* > dat;	//< data
};

#endif
