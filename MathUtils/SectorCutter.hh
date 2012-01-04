#ifndef SECTORCUTTER_HH
#define SECTORCUTTER_HH 1

#include <vector>
#include <math.h>

/// class for dividing a circular region into smaller radial/angular patches
class SectorCutter {
public:
	/// constructor
	SectorCutter(unsigned int N, float R);	
	/// return the total number of sectors
	unsigned int nSectors()  const { return cumdivs[n]; }
	/// identify sector for a given point
	unsigned int sector(float x, float y) const;
	/// identify the ring of this sector
	unsigned int getRing(unsigned int s) const;	
	/// get number of sectors in ring
	unsigned int getNDivs(unsigned int s) const { return s<n?ndivs[s]:0; }
	/// get coordinates for the center of the given sector
	void sectorCenter(unsigned int s, float& x, float& y) const;	
	/// get area of specified sector
	float sectorArea(unsigned int s) const;	
	/// get total area
	float totalArea() const { return 3.141592653589*r*r; }
	/// get outer radius of ring
	float ringRadius(unsigned int rng) const { return rng<n?r*(rng+0.5)/(n-0.5):r; }
	
	unsigned int n;						//< number of concentric rings
	float r;							//< radius of outermost ring
	std::vector<unsigned int> ndivs;	//< number of phi divisions in each ring
	std::vector<unsigned int> cumdivs;	//< cumulative number of divisions in lower rings
};

#endif
