#ifndef QFILE_HH
#define QFILE_HH 1

#include <map>
#include <vector>
#include <string>

/// wrapper for multimap<std::string,std::string> with useful functions
class Stringmap {
public:
	
	/// constructor
	Stringmap(const std::string& s = "");
	/// copy constructor from another Stringmap
	Stringmap(const Stringmap& m);
	/// destructor
	~Stringmap() {}
	
	
	/// insert key/(string)value pair
	void insert(const std::string& s, const std::string& v);
	/// insert key/(double)value
	void insert(const std::string& s, double d);
	/// retrieve key values
	std::vector<std::string> retrieve(const std::string& s) const;	
	/// get first key value (string) or default
	std::string getDefault(const std::string& s, const std::string& d) const;
	/// return number of elements
	unsigned int size() const { return dat.size(); }
	/// serialize to a string
	std::string toString() const;
	
	/// get first key value (double) or default
	double getDefault(const std::string& s, double d) const;
	/// retrieve key values as doubles
	std::vector<double> retrieveDouble(const std::string& s) const;
	/// remove a key
	void erase(const std::string& s);
	
	/// display to screen
	void display(std::string linepfx = "") const;
	
	/// merge data from another stringmap
	void operator+=(const Stringmap& S) { S.mergeInto(*this); }
	
	/// convert to RData format
	//RData* toRData() const;
	
	std::multimap< std::string, std::string > dat;	//< key-value multimap
	
protected:
	
	/// merge data into another stringmap
	void mergeInto(Stringmap& S) const;
};

/// base class for objects that provide stringmaps
class StringmapProvider {
public:
	/// constructor
	StringmapProvider(): Sxtra() {}
	/// destructor
	virtual ~StringmapProvider() {}
	
	/// insert key/(string)value pair
	void insert(const std::string& s, const std::string& v) { Sxtra.insert(s,v); }
	/// insert key/(double)value
	void insert(const std::string& s, double d) { Sxtra.insert(s,d); }
	
	/// provide stringmap from self properties
	Stringmap toStringmap() const {
		Stringmap m = getProperties();
		m += Sxtra;
		return m;
	}
	
	/// display
	void display(std::string linepfx = "") const { toStringmap().display(linepfx); }
	
protected:
	Stringmap Sxtra;
	virtual Stringmap getProperties() const { return Stringmap(); }
};

/// wrapper for multimap<std::string,Stringmap> with useful functions
class QFile {
public:
	
	/// constructor given a string
	QFile(const std::string& s = "", bool readit = true);
	
	/// insert key/(string)value pair
	void insert(const std::string& s, const Stringmap& v);
	/// remove a key
	void erase(const std::string& s);
	/// retrieve values for key
	std::vector<Stringmap> retrieve(const std::string& s) const;
	/// retrieve first value for key
	Stringmap getFirst(const std::string& s, const Stringmap& dflt = Stringmap()) const;
	/// retrieve all sub-key values
	std::vector<std::string> retrieve(const std::string& k1, const std::string& k2) const;
	/// retreive sub-key with default
	std::string getDefault(const std::string& k1, const std::string& k2, const std::string& d) const;
	/// retrieve sub-key as double with default
	double getDefault(const std::string& k1, const std::string& k2, double d) const;
	/// retrieve all sub-key values as doubles
	std::vector<double> retrieveDouble(const std::string& k1, const std::string& k2) const;	
	/// return number of elements
	unsigned int size() const { return dat.size(); }
	/// transfer all data for given key from other QFile
	void transfer(const QFile& Q, const std::string& k);
	
	/// set output file location
	void setOutfile(std::string nm) { name = nm; }
	/// commit data to file
	void commit(std::string outname = "") const;
	/// display to stdout
	void display() const;
	
	/// convert to RData format
	//RData* toRData() const;

protected:
	
	std::string name;								//< name for this object
	std::multimap< std::string, Stringmap > dat;	//< key-value multimap

};

#endif
