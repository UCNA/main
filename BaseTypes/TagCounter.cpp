#include "TagCounter.hh"
#include <stdlib.h>

template<>
TagCounter<int>::TagCounter(Stringmap m) {
	for(std::multimap< std::string, std::string >::iterator it = m.dat.begin(); it != m.dat.end(); it++)
		add(atoi(it->first.c_str()),atof(it->second.c_str()));
}

template<>
TagCounter<unsigned int>::TagCounter(Stringmap m) {
	for(std::multimap< std::string, std::string >::iterator it = m.dat.begin(); it != m.dat.end(); it++)
		add(atoi(it->first.c_str()),atof(it->second.c_str()));
}

template<>
TagCounter<std::string>::TagCounter(Stringmap m) {
	for(std::multimap< std::string, std::string >::iterator it = m.dat.begin(); it != m.dat.end(); it++)
		add(it->first,atof(it->second.c_str()));
}

