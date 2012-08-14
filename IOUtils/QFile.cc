#include "QFile.hh"

#include <iostream>
#include <sstream>
#include <fstream>
#include <utility>
#include "strutils.hh"
#include "PathUtils.hh"
#include "SMExcept.hh"

Stringmap::Stringmap(const std::string& s) {
	std::vector<std::string> pairs = split(s,"\t");
	for(std::vector<std::string>::const_iterator it = pairs.begin(); it!=pairs.end(); it++) {
		std::vector<std::string> keyval = split(*it,"=");
		if(keyval.size() != 2)
			continue;
		dat.insert(std::make_pair(strip(keyval[0]),strip(keyval[1])));
	}
}

Stringmap::Stringmap(const Stringmap& m) {
	for(std::multimap< std::string, std::string >::const_iterator it = m.dat.begin(); it!=m.dat.end(); it++)
		dat.insert(std::make_pair(it->first,it->second));
}

void Stringmap::insert(const std::string& s, const std::string& v) {
	dat.insert(std::make_pair(s,v));
}

void Stringmap::insert(const std::string& s, double d) {
	insert(s,dtos(d));
}

void Stringmap::erase(const std::string& s) { dat.erase(s); }

std::vector<std::string> Stringmap::retrieve(const std::string& s) const {
	std::vector<std::string> v;
	for(std::multimap<std::string,std::string>::const_iterator it = dat.lower_bound(s); it != dat.upper_bound(s); it++)
		v.push_back(it->second);
	return v;
}

std::string Stringmap::getDefault(const std::string& s, const std::string& d) const {
	std::multimap<std::string,std::string>::const_iterator it = dat.find(s);
	if(it == dat.end())
		return d;
	return it->second;
}

std::string Stringmap::toString() const {
	std::string s;
	for(std::multimap<std::string,std::string>::const_iterator it = dat.begin(); it != dat.end(); it++)
		s += "\t" + it->first + " = " + it->second;
	return s;
}

void Stringmap::display(std::string linepfx) const {
	for(std::multimap<std::string,std::string>::const_iterator it = dat.begin(); it != dat.end(); it++)
		std::cout << linepfx <<	it->first << ": " << it->second << "\n";
}


double Stringmap::getDefault(const std::string& k, double d) const {
	std::string s = getDefault(k,"");
	if(!s.size())
		return d;
	std::istringstream ss(s);
	ss >> d;
	return d;
}

std::vector<double> Stringmap::retrieveDouble(const std::string& k) const {
	std::vector<std::string> vs = retrieve(k);
	std::vector<double> v;
	double d;
	for(std::vector<std::string>::const_iterator it = vs.begin(); it != vs.end(); it++) {
		std::istringstream s(*it);
		s >> d;
		v.push_back(d);
	}
	return v;
}

void Stringmap::mergeInto(Stringmap& S) const {
	for(std::multimap<std::string,std::string>::const_iterator it = dat.begin(); it != dat.end(); it++)
		S.insert(it->first,it->second);	
}

/*
RData* Stringmap::toRData() const {
	RDataMem* RM = new RDataMem;
	for(std::multimap<std::string,std::string>::const_iterator it = dat.begin(); it != dat.end(); it++)
		RM->getForced(it->first)->insert(it->second);
	return RM;
}
*/

//----------------------------------------------------------------------------------------------




QFile::QFile(const std::string& fname, bool readit) {
	name = fname;
	if(!readit || name=="")
		return;
	if(!fileExists(fname)) {
		SMExcept e("fileUnreadable");
		e.insert("filename",fname);
		throw(e);
	}
	std::ifstream fin(fname.c_str());
	std::string s;
	while (fin.good()) {
		std::getline(fin,s);
		s = strip(s);
		size_t n = s.find(':');
		if(n==std::string::npos || s[0]=='#') continue;
		std::string key = s.substr(0,n);
		std::string vals = s.substr(n+1);
		vals=strip(vals);
		while(vals.size() && vals[vals.size()-1]=='\\') {
			vals.erase(vals.size()-1);
			std::getline(fin,s);
			s = strip(s);
			vals += '\t';
			vals += s;
		}
		insert(key,Stringmap(vals));
	}
	fin.close();
}

void QFile::insert(const std::string& s, const Stringmap& v) {
	dat.insert(std::make_pair(s,v));
}

void QFile::erase(const std::string& s) { dat.erase(s); }

std::vector<Stringmap> QFile::retrieve(const std::string& s) const {
	std::vector<Stringmap> v;
	for(std::multimap<std::string,Stringmap>::const_iterator it = dat.lower_bound(s); it != dat.upper_bound(s); it++)
		v.push_back(it->second);
	return v;
}

void QFile::transfer(const QFile& Q, const std::string& k) {
	std::vector<Stringmap> v = Q.retrieve(k);
	for(std::vector<Stringmap>::iterator it = v.begin(); it != v.end(); it++)
		insert(k,*it);
}

void QFile::display() const {
	for(std::multimap<std::string, Stringmap>::const_iterator it = dat.begin(); it != dat.end(); it++) {
		std::cout << "--- " << it->first << " ---:\n";
		it->second.display();
	}
}

void QFile::commit(std::string outname) const {
	if(outname=="")
		outname = name;
	makePath(outname,true);
	std::ofstream fout(outname.c_str());
	if(!fout.good()) {
		SMExcept e("fileUnwriteable");
		e.insert("filename",outname);
		throw(e);
	}
	printf("Writing File '%s'.\n",outname.c_str());
	for(std::multimap<std::string, Stringmap>::const_iterator it = dat.begin(); it != dat.end(); it++)
		fout << it->first << ":\t" << it->second.toString() << "\n";
	fout.close();
}

std::vector<std::string> QFile::retrieve(const std::string& k1, const std::string& k2) const {
	std::vector<std::string> v1;
	for(std::multimap<std::string,Stringmap>::const_iterator it = dat.lower_bound(k1); it != dat.upper_bound(k1); it++) {
		std::vector<std::string> v2 = it->second.retrieve(k2);
		for(std::vector<std::string>::const_iterator it2 = v2.begin(); it2 != v2.end(); it2++)
			v1.push_back(*it2);
	}
	return v1;
}

std::vector<double> QFile::retrieveDouble(const std::string& k1, const std::string& k2) const {
	std::vector<double> v1;
	for(std::multimap<std::string,Stringmap>::const_iterator it = dat.lower_bound(k1); it != dat.upper_bound(k1); it++) {
		std::vector<double> v2 = it->second.retrieveDouble(k2);
		for(std::vector<double>::const_iterator it2 = v2.begin(); it2 != v2.end(); it2++)
			v1.push_back(*it2);
	}
	return v1;
}

std::string QFile::getDefault(const std::string& k1, const std::string& k2, const std::string& d) const {
	for(std::multimap<std::string,Stringmap>::const_iterator it = dat.lower_bound(k1); it != dat.upper_bound(k1); it++) {
		std::vector<std::string> v2 = it->second.retrieve(k2);
		if(v2.size())
			return v2[0];
	}
	return d;
}

double QFile::getDefault(const std::string& k1, const std::string& k2, double d) const {
	for(std::multimap<std::string,Stringmap>::const_iterator it = dat.lower_bound(k1); it != dat.upper_bound(k1); it++) {
		std::vector<double> v2 = it->second.retrieveDouble(k2);
		if(v2.size())
			return v2[0];
	}
	return d;
}

Stringmap QFile::getFirst(const std::string& s, const Stringmap& dflt) const {
	std::multimap<std::string,Stringmap>::const_iterator it = dat.find(s);
	if(it == dat.end())
		return dflt;
	return it->second;
}

/*
RData* QFile::toRData() const {
	RDataMem* RM = new RDataMem;
	for(std::multimap<std::string,Stringmap>::const_iterator it = dat.begin(); it != dat.end(); it++) {
		RData* rsub = it->second.toRData();
		RM->insert(it->first,rsub);
		delete(rsub);
	}
	return RM;
}
*/

