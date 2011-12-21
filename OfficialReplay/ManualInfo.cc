#include "ManualInfo.hh"
#include <cfloat>

ManualInfo ManualInfo::MI("../Aux/ManualInfo.txt");

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
	return v;
}
