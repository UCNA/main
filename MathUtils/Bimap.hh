#ifndef BIMAP_HH
#define BIMAP_HH 1

#include <multimap>

/// class for mapping both directions between two sets of values in lists 'A' and 'B'
template<typename T1, typename T2>
class Bimap {
public:
	/// constructor
	Bimap() {}
	/// destructor
	~Bimap() {}
	/// insert a pair of mapped values
	void insert(T1 a, T2 b) {
		m1.insert(std::pair<T1,T2>(a,b));
		m2.insert(std::pair<T2,T1>(b,a));
	}
	/// locate entries in 'B' mapped to from given value in A
	std::pair< std::multimap<T1,T2>::iterator, std::multimap<T1,T2>::iterator > findA(T1 a) { return m1.equal_range(a); }
	/// locate entries in 'A' mapped from given value in B
	std::pair< std::multimap<T2,T1>::iterator, std::multimap<T2,T1>::iterator > findB(T2 b) { return m2.equal_range(b); }
	/// get list of Bs mapped from given A
	std::vector<T2> get_Bs_for_A(T1 a) {
		std::vector<T2> v;
		std::pair< std::multimap<T1,T2>::iterator, std::multimap<T1,T2>::iterator > its = findA(a);
		for(std::multimap<T1,T2>::iterator it = its.first; it != its.second; it++)
			v.push_back(it.second);
		return v;
	}
	/// get list of As mapped from given B
	std::vector<T1> get_As_for_B(T2 b) {
		std::vector<T1> v;
		std::pair< std::multimap<T2,T1>::iterator, std::multimap<T2,T1>::iterator > its = findB(b);
		for(std::multimap<T2,T1>::iterator it = its.first; it != its.second; it++)
			v.push_back(it.second);
		return v;
	}
	
	/// delete entry in A and its mappings
	void delete_A(T1 a);
	/// delete entry in B and its mappings
	void delete_B(T2 b);
	
protected:
	std::multimap<T1,T2> m1;	//< internal map from A to B
	std::multimap<T2,T1> m2;	//< internal map from B to A
}

#endif
