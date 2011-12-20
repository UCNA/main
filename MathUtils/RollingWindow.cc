#include "RollingWindow.hh"
#include <utility>

void RollingWindow::addCount(double t, double w) {
	itms.push_front(std::make_pair(t,w));
	sw += w;
	while(itms.size()>nMax)
		popExcess();
	moveTimeLimit(t);
}

void RollingWindow::moveTimeLimit(double t) {
	while(itms.size() && t - itms.back().first > lMax)
		popExcess();
}

void RollingWindow::popExcess() {
	sw -= itms.back().second;
	itms.pop_back();
	if(!itms.size())
		sw=0;
}
