#include "RData.hh"
#include <sstream>
#include <utility>

RData* RData::NullRData = new RData();

std::string RData::toString(bool displaymode, const std::string& pfx) {
	if(!size()) {
		if(displaymode)
			return "";
		else
			return "{}";
	}
	
	std::string closelist = "]";
 	std::string dictfirst = "{";
	std::string dictnext = ",\n"+pfx;
	std::string closedict = "}";
	std::string qt = "\"";
	std::string listfirst = ": [";
	std::string listnext = ",\n\t"+pfx;
	std::string ppfx = "\t";
	if(displaymode) {
		qt = "";
		dictfirst = pfx+"+";
		dictnext = pfx+"|";
		closedict = closelist = "";
		listfirst = "\n";
		listnext = "";
		ppfx = "|\t";
	}
	
	std::string s = dictfirst;
	std::vector<std::string> keys = getKeys();
	for(std::vector<std::string>::iterator it = keys.begin(); it != keys.end(); it++) {
		if(it != keys.begin())
			s += dictnext;
		s += qt + *it + qt;
		std::vector<RData*> subdat = getSubdata(*it);
		for(std::vector<RData*>::iterator it2 = subdat.begin(); it2 != subdat.end(); it2++) {
			if(it2 == subdat.begin())
				s += listfirst;
			else
				s += listnext;
			s += (*it2)->toString(displaymode, ppfx + pfx);
		}
		s += closelist;
	}
	s += closedict;
	return s;
}

float RData::getFirstD(double dflt) const {
	std::string s = getFirstKey("");
	if(!s.size())
		return dflt;
	std::istringstream ss(s);
	ss >> dflt;
	return dflt;
}

RData* RData::getForced(const std::string& key) {
	RData* RD = getFirst(key);
	if(!RD->isNull())
		return RD;
	return insert(key);
}

std::vector<RData*> RData::getSubdata(const std::vector<std::string>& keys, unsigned int kdepth) {
	std::vector<RData*> v;
	if(kdepth >= keys.size())
		return v;
	std::vector<RData*> submatches = getSubdata(keys[kdepth]);
	if(kdepth == keys.size()-1)
		return submatches;
	for(std::vector<RData*>::iterator it = submatches.begin(); it != submatches.end(); it++) {
		std::vector<RData*> v2 = (*it)->getSubdata(keys,kdepth+1);
		while(v2.size()) {
			v.push_back(v2.back());
			v2.pop_back();
		}
	}
	return v;
}

RData* RData::getFirst(const std::vector<std::string>& keys, unsigned int kdepth) {
	if(kdepth >= keys.size())
		return RData::NullRData;
	if(kdepth == keys.size()-1)
		return getFirst(keys[kdepth]);
	std::vector<RData*> submatches = getSubdata(keys[kdepth]);
	for(std::vector<RData*>::iterator it = submatches.begin(); it != submatches.end(); it++) {
		RData* RD = (*it)->getFirst(keys,kdepth+1);
		if(!RD->isNull())
			return RD;
	}
	return NullRData;	
}
	
RData* RData::getFirstByPath(const std::string& path, const std::string& pathsep) {
	return getFirst(split(strip(path,pathsep),pathsep));
}

std::vector<RData*> RData::getSubdataByPath(const std::string& path, const std::string& pathsep) {
	return getSubdata(split(strip(path,pathsep),pathsep));
}

RDataMem::RDataMem(RData* KD): RData() {
	std::vector<std::string> keys = KD->getKeys();
	for(std::vector<std::string>::iterator it = keys.begin(); it != keys.end(); it++) {
		std::vector<RData*> subdat = KD->getSubdata(*it);
		for(std::vector<RData*>::iterator it2 = subdat.begin(); it2 != subdat.end(); it2++)
			insert(*it,*it2);
	}
}

RDataMem::~RDataMem() {
	for(std::multimap<std::string,RDataMem*>::iterator it = dat.begin(); it != dat.end(); it++)
		delete(it->second);
	dat.clear();
}

std::vector<std::string> RDataMem::getKeys() const {
	std::vector<std::string> v;
	for(std::multimap<std::string,RDataMem*>::const_iterator it = dat.begin(); it != dat.end(); it++)
		if(!v.size() || it->first != v.back())
			v.push_back(it->first);
	return v;
}

std::string RDataMem::getFirstKey(std::string dflt) const {
	if(dat.size())
		return dat.begin()->first;
	return dflt;
}

RData* RDataMem::insert(const std::string& key) {
	RDataMem* KD = new RDataMem();
	dat.insert(std::make_pair(key,KD));
	return KD;
}

RData* RDataMem::insert(const std::string& key, RData* KD) {
	RDataMem* nKD = new RDataMem(KD);
	dat.insert(std::make_pair(key,nKD));
	return nKD;
}

std::vector<RData*> RDataMem::getSubdata(const std::string& key) {
	std::vector<RData*> v;
	for(std::multimap<std::string,RDataMem*>::const_iterator it = dat.lower_bound(key); it != dat.upper_bound(key); it++)
		v.push_back(it->second);
	return v;
}

RData* RDataMem::getFirst(const std::string& key) {
	std::multimap<std::string,RDataMem*>::const_iterator it = dat.find(key);
	if(it != dat.end())
		return it->second;
	return RData::NullRData;
}
