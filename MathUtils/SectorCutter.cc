#include "SectorCutter.hh"
#include <algorithm>
#include <math.h>
#include <stdlib.h>

SectorCutter::SectorCutter(unsigned int N, float R): n(N), r(R) {
	cumdivs.push_back(0);
	for(unsigned int i=0; i<n; i++) {
		if(i==0)
			ndivs.push_back(1);
		else
			ndivs.push_back((unsigned int)ceil(2*3.1415926535*i));
		cumdivs.push_back(cumdivs[i]+ndivs[i]);
	}
}

unsigned int SectorCutter::sector(float x, float y) const {
	
	// which ring this belongs in
	unsigned int rng = 0;
	float rrel = sqrt(x*x+y*y)/r*(n-0.5);
	if(n>1 && rrel>=0.5)
		rng=(unsigned int)(rrel+0.5);
	if(rng>=n || (n==1 && rrel>1))
		return nSectors();
	
	// which phi this belongs in
	int ph = int(((atan2(-y,-x)+3.141592653589)/6.28318531)*ndivs[rng]);
	if(ph<0) ph = 0;
	if(ph==(int)ndivs[rng]) ph--;
	
	return cumdivs[rng]+ph;
}

unsigned int SectorCutter::getRing(unsigned int s) const {
	std::vector<unsigned int>::const_iterator it = std::upper_bound(cumdivs.begin(),cumdivs.end(),s);
	return int(it-cumdivs.begin())-1;
}

void SectorCutter::sectorBounds(unsigned int s, float& r0, float& r1, float& ph0, float& ph1) const {
	if(s>=nSectors()) {
		r0=r1=r;
		ph0=0;
		ph1=6.28318531;
		return;
	}
	unsigned int rng = getRing(s);
	r0 = rng>0?ringRadius(rng-1):0;
	r1 = ringRadius(rng);
	s -= cumdivs[rng];
	ph0 = 6.28318531*float(s)/float(ndivs[rng]);
	ph1 = 6.28318531*float(s+1)/float(ndivs[rng]);
}

void SectorCutter::sectorCenter(unsigned int s, float& x, float& y) const {
	if(s>=nSectors()||!s) {
		x=y=0;
		return;
	}
	float r0,r1,ph0,ph1;
	sectorBounds(s,r0,r1,ph0,ph1);
	float rc = 0.5*(r0+r1);
	float phc = 0.5*(ph0+ph1);
	x = rc*cos(phc);
	y = rc*sin(phc);
}

float SectorCutter::sectorCenterRadius(unsigned int s) const {
	if(s>=nSectors()||!s) return 0;
	float r0,r1,ph0,ph1;
	sectorBounds(s,r0,r1,ph0,ph1);
	return 0.5*(r0+r1);
}

float SectorCutter::sectorArea(unsigned int s) const {
	if(s>=nSectors()) return 0;
	unsigned int rng = getRing(s);
	float r1 = ringRadius(rng);
	if(rng==0)
		return 3.141592653589*r1*r1;
	float r2 = ringRadius(rng-1);
	return 3.141592653589*(r1-r2)*(r1+r2)/getNDivs(rng);
}

double randunif(double x0, double x1) {
	return x0+(x1-x0)*double(rand())/double(RAND_MAX);
}

void SectorCutter::randPos(unsigned int s, float& x, float& y) const {
	float r0,r1,ph0,ph1;
	sectorBounds(s,r0,r1,ph0,ph1);
	double pr = sqrt(randunif(r0*r0,r1*r1));
	double pth = randunif(ph0,ph1);
	x = pr*cos(pth);
	y = pr*sin(pth);
}

