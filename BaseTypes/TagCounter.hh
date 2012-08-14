#ifndef TAGCOUNTER_HH
#define TAGCOUNTER_HH 1

#include <map>
#include <iostream>
#include <istream>
#include <sstream>
#include "QFile.hh"

template<typename T>
class TagCounter {
public:
	/// constructor
	TagCounter(Stringmap m = Stringmap());
	/// destructor
	~TagCounter() {}
	/// add counts
	void add(const T& itm, double c);
	/// add another counter
	void operator+=(const TagCounter<T>& c);
	/// multiply all counts
	void scale(double s);
	/// make into Stringmap
	Stringmap toStringmap();
	/// get number of counted items
	unsigned int nTags() const { return counts.size(); }
	/// get total counts on all objects
	double total() const;
	/// get count for given item
	double operator[](const T& itm) const;

	std::map<T,double> counts;	//< counts per object
};

template<typename T>
void TagCounter<T>::add(const T& itm, double c) {
	counts[itm] += c;
}

template<typename T>
void TagCounter<T>::operator+=(const TagCounter<T>& c) {
	for(typename std::map<T,double>::const_iterator it = c.counts.begin(); it != c.counts.end(); it++)
		add(it->first,it->second);
}

template<typename T>
void TagCounter<T>::scale(double s) {
	if(s==1) return;
	for(typename std::map<T,double>::iterator it = counts.begin(); it != counts.end(); it++)
		it->second *= s;
}

template<typename T>
Stringmap TagCounter<T>::toStringmap() {
	Stringmap m;
	for(typename std::map<T,double>::const_iterator it = counts.begin(); it != counts.end(); it++) {
		std::ostringstream s;
		s << (*it).first;
		m.insert(s.str(),dtos(it->second));
	}
	return m;
}

template<typename T>
double TagCounter<T>::total() const {
	double d = 0;
	for(typename std::map<T,double>::const_iterator it = counts.begin(); it != counts.end(); it++)
		d += it->second;
	return d;
}

template<typename T>
double TagCounter<T>::operator[](const T& itm) const {
	typename std::map<T,double>::const_iterator it = counts.find(itm);
	if(it==counts.end()) return 0;
	return it->second;
}

#endif
