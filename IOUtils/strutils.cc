#include "strutils.hh"
#include <stdlib.h>
#include <math.h>

std::string dtos(double d, const std::string& badnum) {
	char c[16];
	if(badnum.size() && !(d==d && !isinf(d))) sprintf(c,"%s",badnum.c_str());
	else sprintf(c,"%g",d);
	return std::string(c);
}

std::string itos(int i) {
	char c[32];
	sprintf(c,"%i",i);
	return std::string(c);	
}

std::string itosRN(int i) {
	if(!i) return "0";
	std::string s;
	if(i<0) { s += "-"; i = -i; }
	while(i>=1000) { s += "M"; i-=1000; }
	while(i>=900) { s += "CM"; i-=900; }
	while(i>=500) { s += "D"; i-=500; }
	while(i>=400) { s += "CD"; i-=400; }
	while(i>=100) { s += "C"; i-=100; }
	while(i>=90) { s += "XC"; i-=90; }
	while(i>=50) { s += "L"; i-=50; }
	while(i>=40) { s += "XL"; i-=40; }
	while(i>=10) { s += "X"; i-=10; }
	while(i>=9) { s += "IX"; i-=9; }
	while(i>=5) { s += "V"; i-=5; }
	while(i>=4) { s += "IV"; i-=4; }
	while(i>0) { s += "I"; i--; }
	return s;
}

std::string vtos(const double* st, const double* en, std::string sep) {
	std::string s = "";
	if(st==en)
		return s;
	s = dtos(*st);
	for(const double* it = st+1; it != en; it++)
		s += sep+dtos(*it);
	return s;
}

std::string vtos(const std::vector<double>& ds,std::string sep) { return vtos(&*ds.begin(),&*ds.end(),sep); }

std::string vtos(const float* st, const float* en, std::string sep) {
	std::string s = "";
	if(st==en)
		return s;
	s = dtos(*st);
	for(const float* it = st+1; it != en; it++)
		s += sep+dtos(*it);
	return s;
}

std::string vtos(const std::vector<float>& ds,std::string sep) { return vtos(&*ds.begin(),&*ds.end(),sep); }

std::string vtos(const int* st, const int* en, std::string sep) {
	std::string s = "";
	if(st==en)
		return s;
	s = itos(*st);
	for(const int* it = st+1; it != en; it++)
		s += sep+itos(*it);
	return s;
}

std::string vtos(const std::vector<int>& ds,std::string sep) { return vtos(&*ds.begin(),&*ds.end(),sep); }

std::string ctos(char c) {
	char ch[3];
	sprintf(ch,"%c",c);
	return std::string(ch);		
}

std::string lower(std::string s) {
	std::transform(s.begin(), s.end(), s.begin(), (int(*)(int))tolower);
	return s;
}

std::string upper(std::string s) {
	std::transform(s.begin(), s.end(), s.begin(), (int(*)(int))toupper);
	return s;
}

std::string replace(std::string s, char o, char n) {
	std::string::size_type found = s.find_first_of(o);
	while( found != std::string::npos ) {
		s[found] = n;
		found = s.find_first_of(o,found+1);
	}
	return s;
}

bool startsWith(const std::string& a, const std::string& b) { return a.substr(0,b.size()) == b; }

std::vector<std::string> split(const std::string& s, const std::string splitchars) {
	std::vector<std::string> v;
	size_t p = 0;
	while(p<s.size()) {
		size_t wstart = s.find_first_not_of(splitchars,p);
		if(wstart == std::string::npos)
			break;
		p = s.find_first_of(splitchars,wstart);
		if(p == std::string::npos)
			p = s.size();
		v.push_back(s.substr(wstart,p-wstart));
	}
	return v;
}

std::string join(const std::vector<std::string>& ss, const std::string& sep) {
	std::string s = "";
	if(!ss.size())
		return s;
	s = ss[0];
	for(std::vector<std::string>::const_iterator it = ss.begin()+1; it < ss.end(); it++)
		s += sep + *it;
	return s;
}

std::string strip(const std::string& s, const std::string stripchars) {
	size_t wstart = s.find_first_not_of(stripchars);
	if(wstart == std::string::npos)
		return "";
	size_t wend = s.find_last_not_of(stripchars);
	return s.substr(wstart,wend-wstart+1);
}

std::vector<double> sToDoubles(const std::string& s, const std::string splitchars) {
	std::vector<double> v;
	std::vector<std::string> words = split(s,splitchars);
	for(unsigned int i=0; i<words.size(); i++)
		v.push_back(atof(words[i].c_str()));
	return v;
}

std::vector<float> sToFloats(const std::string& s, const std::string splitchars) {
	std::vector<float> v;
	std::vector<std::string> words = split(s,splitchars);
	for(unsigned int i=0; i<words.size(); i++)
		v.push_back(atof(words[i].c_str()));
	return v;
}

std::vector<int> sToInts(const std::string& s, const std::string splitchars) {
	std::vector<int> v;
	std::vector<std::string> words = split(s,splitchars);
	for(unsigned int i=0; i<words.size(); i++)
		v.push_back(atoi(words[i].c_str()));
	return v;
}

std::vector< std::vector<float> > readArray(std::ifstream& fin, unsigned int minitems, const std::string splitchars) {
	std::vector< std::vector<float> > a;
	std::string s;
	while (fin.good()) {
		std::getline(fin,s);
		std::vector<float> v = sToFloats(s,splitchars);
		if(v.size() >= minitems)
			a.push_back(v);
	}
	return a;
}
