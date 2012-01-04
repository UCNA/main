#include "SectorCutter.hh"
#include <algorithm>

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
	float rrel = sqrt(x*x+y*y)/r*n;
	unsigned int m = 0;
	if(n>1 && rrel>=0.5)
		m=(unsigned int)(rrel+0.5);
	if(m>=n || (n==1 && rrel>1))
		return nSectors();
	
	// which phi this belongs in
	int ph = int(((atan2(-y,-x)+3.141592653589)/6.28318531)*ndivs[m]);
	if(ph<0) ph = 0;
	if(ph==(int)ndivs[m]) ph--;
	
	return cumdivs[m]+ph;
}

unsigned int SectorCutter::getRing(unsigned int s) const {
	std::vector<unsigned int>::const_iterator it = std::upper_bound(cumdivs.begin(),cumdivs.end(),s);
	return int(it-cumdivs.begin())-1;
}

void SectorCutter::sectorCenter(unsigned int s, float& x, float& y) const {
	if(s>=nSectors()) {
		x=y=0;
		return;
	}
	unsigned int n0 = getRing(s);
	if(n0==0) {
		x=y=0;
		return;
	}
	float r0 = r*float(n0)/(float(n)-0.5);
	s -= cumdivs[n0];
	float phi = 6.28318531*float(s+0.5)/float(ndivs[n0]);
	x = r0*cos(phi);
	y = r0*sin(phi);
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
