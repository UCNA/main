#ifndef ROLLINGWINDOW_HH
#define ROLLINGWINDOW_HH

#include <utility>
#include <deque>
#include <cfloat>
#include <cmath>

/// rolling window averager with length and time limit
class RollingWindow {
public:
	/// constructor
	RollingWindow(unsigned int n, double l=FLT_MAX): nMax(n), lMax(l), sw(0), sww(0) {}
	
	/// introduce next element
	void addCount(double t, double w=1.);
	/// remove item off back
	void popExcess();
	/// update leading time limit without adding items
	void moveTimeLimit(double t);	
	/// get items sum
	inline double getSum() const { return sw; }
	/// get items count
	inline unsigned int getCount() const { return itms.size(); }
	/// get average
	inline double getAvg() const { return sw/itms.size(); }
	/// get RMS deviation
	inline double getRMS() const { return sqrt((sww-sw*sw)/itms.size()); }
	/// get average excluding value
	inline double getAvgExcl(double x) const { return (sw-x)/(itms.size()-1); }
	/// get RMS excluding value
	inline double getRmsExcl(double x) const { return sqrt((sww-x*x-(sw-x)*(sw-x))/(itms.size()-1)); }
	/// clear all data
	void clear() { sw = sww = 0; itms.clear(); }
	
	unsigned int nMax;	//< maximum number of items to track
	double lMax;		//< maximum time span to track from leading object
	
protected:
	std::deque< std::pair<double,double> > itms;	//< items in window
	double sw;										//< sum of weights
	double sww;										//< sum of weights squared
};


#endif
