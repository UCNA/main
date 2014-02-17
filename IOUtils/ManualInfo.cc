#include "ManualInfo.hh"
#include "SMExcept.hh"
#include "PathUtils.hh"
#include <cfloat>
#include <cmath>
#include <algorithm>
#include <stdio.h>

ManualInfo ManualInfo::MI(getEnvSafe("UCNA_AUX")+"/ManualInfo.txt");

struct rangeSorter {
	rangeSorter(const std::string& key1, const std::string& key2): k1(key1), k2(key2) {}
	const std::string k1,k2;
	bool operator() (const Stringmap& m1, const Stringmap& m2) const {
		return fabs(m1.getDefault(k2,DBL_MAX) - m1.getDefault(k1,0)) < fabs(m2.getDefault(k2,DBL_MAX) - m2.getDefault(k1,0));
	}
};

std::vector< std::pair<double,double> > ManualInfo::getRanges(const std::string& key, const std::string& k1, const std::string& k2) const {
	std::vector< std::pair<double,double> > v;
	for(std::multimap<std::string,Stringmap>::const_iterator it = dat.lower_bound(key); it != dat.upper_bound(key); it++)
		v.push_back(std::make_pair(it->second.getDefault(k1,0.),it->second.getDefault(k2,0.)));
	return v;
}

std::vector<Stringmap> ManualInfo::getInRange(const std::string& key,
											  const double x,
											  const std::string& k1,
											  const std::string& k2) const {
	std::vector<Stringmap> v;
	for(std::multimap<std::string,Stringmap>::const_iterator it = dat.lower_bound(key); it != dat.upper_bound(key); it++)
		if(it->second.getDefault(k1,DBL_MAX) <= x && x <= it->second.getDefault(k2,0.))
			v.push_back(it->second);
	rangeSorter RS(k1,k2);
	std::sort(v.begin(),v.end(),RS);
	return v;
}

RangeCut::RangeCut(const Stringmap& m): start(m.getDefault("start",0.0)), end(m.getDefault("end",0.0)) {}

Stringmap loadCut(RunNum rn, const std::string& cutName) {
	std::vector<Stringmap> v = ManualInfo::MI.getInRange(cutName,rn);
	if(!v.size()) {
		SMExcept e("missingCut");
		e.insert("cutName",cutName);
		e.insert("runNum",rn);
		throw(e);
	}
	return v[0];
}

void loadRangeCut(RunNum rn, CutVariable& c, const std::string& cutName) {
	c.R = RangeCut(loadCut(rn,cutName));
	printf("Loaded cut %s/%i = (%g,%g)\n",cutName.c_str(),rn,c.R.start,c.R.end);
}
