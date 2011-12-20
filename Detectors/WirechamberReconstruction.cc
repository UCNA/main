#include "WirechamberReconstruction.hh"
#include "strutils.hh"
#include <cfloat>

std::vector<bool> getLiveWires(RunNum rn, Side s, AxisDirection d) {
	
	std::vector<bool> isLive;
	unsigned int nWires = 16;
	
	for(unsigned int i=0; i<nWires; i++) {
		
		isLive.push_back(false);
		
		// 2008 data, dead East wire 8
		if(s==EAST && i==8) {
			if( 7200 < rn && rn < 7318 && d == Y_DIRECTION)
				continue;
			if(7318 <= rn && rn <= 7551 && d == X_DIRECTION)
				continue;
			if(7553 <= rn && rn <= 7588 && d == Y_DIRECTION)
				continue;
			if(7716 <= rn && rn <= 7732 && d == X_DIRECTION)
				continue;
			if(7837 <= rn && rn <= 7840 && d == Y_DIRECTION)
				continue;
		}
		// 2008 data, dead East edge wires
		if(rn < 8460 && s==EAST && i==15)
			continue;
		// 2009 data
		if(s==EAST && d == X_DIRECTION && i==8) {
			if(11984 <= rn && rn <= 11991) //missing
				continue;
			if(12388 <= rn && rn <= 12426) //iffy
				continue;
			if(12430 <= rn && rn <= 12432) //missing
				continue;
			if(12433 <= rn && rn <= 12534) //iffy
				continue;
			if(12535 <= rn && rn <= 12588) //missing
				continue;
		}
		
		isLive.back() = true;
	}
	
	return isLive;
}

std::vector<std::string> getCathodeNames(RunNum rn, Side s, AxisDirection d) {
	std::vector<bool> liveWires = getLiveWires(rn,s,d);
	std::vector<std::string> v;
	std::string sensname = std::string("MWPC")+(s==EAST?"E":"W")+(d==X_DIRECTION?"x":"y");
	for(unsigned int i=0; i<liveWires.size(); i++)
		if(liveWires[i])
			v.push_back(sensname+itos(i+1));
	return v;
}

std::vector<unsigned int> getPadcNumbers(RunNum rn, Side s, AxisDirection d) {
	
	//                              0    1    2    3    4    5    6    7    8    9    10   11   12   13   14   15
	unsigned int EX_Padc_nums[] = { 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231 };
	unsigned int EY_Padc_nums[] = { 20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  210, 211, 212, 213, 214, 215 };
	unsigned int WX_Padc_nums[] = { 31,  30,  29,  28,  27,  26,  25,  24,  23,  22,  21,  20,  19,  18,  17,  16 };
	unsigned int WY_Padc_nums[] = { 0,   1,   2,   3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  13,  14,  15 };	
	
	std::vector<bool> liveWires = getLiveWires(rn,s,d);
	std::vector<unsigned int> padcNums;
	
	for(unsigned int i=0; i<liveWires.size(); i++) {
		if(!liveWires[i])
			continue;
		if(s == EAST && d == X_DIRECTION)
			padcNums.push_back(EX_Padc_nums[i]);
		if(s == EAST && d == Y_DIRECTION)
			padcNums.push_back(EY_Padc_nums[i]);
		if(s == WEST && d == X_DIRECTION)
			padcNums.push_back(WX_Padc_nums[i]);
		if(s == WEST && d == Y_DIRECTION)
			padcNums.push_back(WY_Padc_nums[i]);
	}
	
	return padcNums;
}

std::vector<float> calcWirePositions(RunNum rn, Side s, AxisDirection d, float wireSpacing) {
	
	std::vector<bool> liveWires = getLiveWires(rn,s,d);
	unsigned int nWires = liveWires.size();
	
	std::vector<float> wirepos;
	for(unsigned int i=0; i<nWires; i++)
		if(liveWires[i])
			wirepos.push_back( (nWires-1)*wireSpacing*0.5 - i*wireSpacing );
	
	return wirepos;
}

int sign(float v) { return v > 0 ? 1 : (v < 0 ? -1 : 0); }

wireHit mpmGaussianPositioner(const std::vector<float>& wirepos, float* wireValues, const float* wirePeds) {
	
	float x1,x2,x3,y1,y2,y3,d;
	const double sigma0 = 0.75;	// TODO is this really a sigma? re-understand this
	unsigned int nWires = wirepos.size();
	
	std::vector<float> xs;
	std::vector<float> ys;
	int maxn = -1;
	
	// initialize values in h
	wireHit h;
	h.nClipped = 0;
	h.maxWire = 0;
	h.cathodeSum = 0;
	h.multiplicity = 0;
	h.errflags = 0;
	h.maxValue = -1;
	
	for(unsigned int c=0; c<nWires; c++) {
		// values above 3950 count as "clipped"
		bool isClipped = wireValues[c] > 3950;
		// pedestal subtract wire readout
		wireValues[c] -= wirePeds[c];
		if(isClipped) {
			h.nClipped++;
		} else {
			// record in usable wires
			xs.push_back(wirepos[c]);
			ys.push_back(wireValues[c]);
		}
		// values above 70 count towards multiplicity
		if(wireValues[c]>70)
			h.multiplicity++;
		h.cathodeSum += wireValues[c];
		// check if this is the maximum value found so far
		if(wireValues[c] >= h.maxValue && !isClipped) {
			// if wires are tied for max value, randomly choose which is labeled as maxWire
			if(wireValues[c] == h.maxValue && rand()%2)
				continue;
			h.maxValue = wireValues[c];
			h.maxWire = c;
			maxn = xs.size()-1;
		}
	}
	
	if(h.nClipped > 0)
		h.errflags |= WIRES_CLIPPED;
	if(h.nClipped > 1)
		h.errflags |= WIRES_MULTICLIP;
	if(h.maxValue <= 0)
		h.errflags |= WIRES_NONE;
	if(maxn==0 || maxn==int(xs.size())-1)
		h.errflags |= WIRES_EDGE;
	
	// no usable wires? There's no hit to reconstruct.
	if(h.maxValue <= 0 || !xs.size()) {
		h.center = h.width = 0;
		return h;
	}
	// only one usable wire?
	if(xs.size()==1) {
		h.errflags |= WIRES_SINGLET;
		h.center = xs[0];
		return h;
	}
	
	// special case for hits near edge (assumes fixed width)
	if(h.errflags & WIRES_EDGE) {
		
		y1 = ys[maxn];
		x1 = xs[maxn];
		if(maxn==0)
			maxn++;
		else
			maxn--;
		y2 = ys[maxn];
		x2 = xs[maxn];
		
		h.avgCenter = (x1*y1+x2*y2)/(y1+y2);
		
		// next wire over is negative:
		// isolated edge wire firing
		if(y2 <= 0) {
			h.errflags |= WIRES_SINGLET;
			y2 = 1;
		}
		
		h.width = sigma0*fabs(x2-x1);
		float l = 0.5 + (log(y2)-log(y1))*sigma0*sigma0;
		h.center = (1-l)*x1 + l*x2;
		h.parabCenter = -1000;
		return h;
	}
	
	// collect x and y values
	y1 = ys[maxn-1];
	y2 = ys[maxn];
	y3 = ys[maxn+1];
	x1 = xs[maxn-1];
	x2 = xs[maxn];
	x3 = xs[maxn+1];
	float dx1 = x2-x1;
	float dx2 = x3-x2;
	
	// "charge center" average from main wires
	h.avgCenter = (x1*y1+x2*y2+x3*y3)/(y1+y2+y3);
	
	if(fabs(dx1-dx2)<0.001) {
		// normal case: uniformly spaced wires; convert from parabola center to position based on assumed sigma
		float x = (y3-y1)/(4.0*y2-2.0*(y1+y3));
		h.parabCenter = x2+dx1*x;
		float sigma2 = 0.65*0.65;
		float gx = sign(x)*sigma2*log(( 2.0*exp(1.0/(2.0*sigma2))*fabs(x) + sqrt( 1.0 + 4.0*( exp(1.0/sigma2) - 1.0 )*x*x ) )/(1.0+2.0*fabs(x)));
		h.center = x2+dx1*gx;
	} else {
		// nonuniform wire spacing (or wires missing due to clipping)
		h.errflags |= WIRES_NONUNIF;
		h.center = -1000;
		h.parabCenter = -1000;
	}
	
	// isolated singlet: wires on both sides negative
	if(y1 <= 0 && y3 <= 0) {
		h.errflags |= WIRES_SINGLET;
		h.width = 0;
		h.center = x2;
		return h;
	}
	
	// isolated doublet (wires on one side negative): fixed sigma reconstruction
	if(y1 <= 0 || y3 <= 0) {
		if(y3 <= 0) {
			y3 = y1;
			x3 = x1;
		}
		h.errflags |= WIRES_DOUBLET;
		h.width = sigma0*fabs(x3-x2);
		float l = 0.5 + (log(y3)-log(y2))*sigma0*sigma0;
		h.center = (1-l)*x2 + l*x3;
		return h;
	}
	
	// calculate gaussian center for non-uniform wires, width for all wires
	y1 = log(y1);
	y2 = log(y2);
	y3 = log(y3);
	d = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2);
	h.width = sqrt((x1-x2)*(x2-x3)*(x3-x1)/(-2*d));
	if(h.center == -1000)
		h.center = 0.5*(x1*x1*(y2-y3)  + x2*x2*(y3-y1) + x3*x3*(y1-y2)) / d;
	
	return h;
}
