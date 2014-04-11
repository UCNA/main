#ifndef FLOATERR_HH
#define FLOATERR_HH

#include <string>

/// float value with an error estimate
class float_err {
public:
	/// constructor
	float_err(float c=0, float dc=0): x(c), err(dc) {}
	/// constructor from std::string
	float_err(const std::string& s);
	/// converto to a string
	std::string toString() const;
	
	float x;	///< central value
	float err;	///< error
};

/// add float_errs assuming independent statistics
float_err operator+(float_err a, float_err b);
/// multiply float_err by float
float_err operator*(float a, float_err b);
/// statistically weighted sum of n values with errors, useful for combining PMT results
float_err weightedSum(unsigned int n, const float_err* d);
/// measure of combined statistical proximity of points d to central value c
float proximity(unsigned int n, const float_err* d, float_err c);

#endif
