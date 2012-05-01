#include "RollingWindow.hh"
#include <utility>
#include <cassert>

void RollingWindow::addCount(double t, double w) {
	itms.push_front(std::make_pair(t,w));
	sw += w;
	sww += w*w;
	while(itms.size()>nMax)
		popExcess();
	moveTimeLimit(t);
}

void RollingWindow::moveTimeLimit(double t) {
	while(itms.size() && t - itms.back().first > lMax)
		popExcess();
}

void RollingWindow::popExcess() {
	assert(itms.size());
	double w = itms.back().second;
	sw -= w;
	sww -= w*w;
	itms.pop_back();
	if(!itms.size())
		sw=sww=0;
}
