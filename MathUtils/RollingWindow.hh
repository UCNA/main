#ifndef ROLLINGWINDOW_HH
#define ROLLINGWINDOW_HH 1

#include <utility>
#include <deque>
#include <cfloat>

/// rolling window averager with length and time limit
class RollingWindow {
public:
	/// constructor
	RollingWindow(unsigned int n, double l=FLT_MAX): nMax(n), lMax(l), sw(0) {}
	
	/// introduce next element
	void addCount(double t, double w=1.);
	/// update leading time limit without adding items
	void moveTimeLimit(double t);	
	/// get items sum
	inline double getSum() const { return sw; }
	/// get items count
	inline unsigned int getCount() const { return itms.size(); }
	
	unsigned int nMax;	//< maximum number of items to track
	double lMax;		//< maximum time span to track from leading object
	
protected:
	/// remove item off back
	void popExcess();
	std::deque< std::pair<double,double> > itms;	//< items in window
	double sw;										//< sum of weights
};


#endif
