#include "FloatErr.hh"
#include "strutils.hh"
#include <math.h>

float_err::float_err(const std::string& s): x(0), err(0) {
	std::vector<std::string> v = split(s,"~");
	if(v.size()) x = atof(v[0].c_str());
	if(v.size()>1) err = atof(v[1].c_str());
}

std::string float_err::toString() const {
	return dtos(x)+"~"+dtos(err);
}

float_err operator+(float_err a, float_err b) {
	return float_err(a.x+b.x,sqrt(a.err*a.err + b.err*b.err));
}

float_err operator*(float a, float_err b) {
	return float_err(a*b.x,a*b.err);
}

float_err weightedSum(unsigned int n, const float_err* d) {
	float_err sum;
	sum.x = sum.err = 0;
	for(unsigned int i=0; i<n; i++) {
		sum.x += d[i].x/(d[i].err*d[i].err);
		sum.err += 1/(d[i].err*d[i].err);
	}
	sum.x /= sum.err;
	sum.err = 1/sqrt(sum.err);
	return sum;
}

float proximity(unsigned int n, const float_err* d, float_err c) {
	float psum = 0;
	for(unsigned int i=0; i<n; i++)
		psum += (d[i].x-c.x)*(d[i].x-c.x)/(d[i].err*d[i].err + c.err*c.err);
	return sqrt(psum/n);
}
