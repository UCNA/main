#include "WirechamberCalibrator.hh"
#include "SMExcept.hh"
#include <cmath>
#include <cassert>
#include <algorithm>

#define kUndefinedPosition 666

WirechamberCalibrator::WirechamberCalibrator(RunNum rn, CalDB* cdb): anodeP(cdb->getAnodePositioningCorrector(rn)) {
	assert(anodeP);
	anodeP->setNormAvg();
	for(Side s = EAST; s <= WEST; ++s) {
		anodeGainCorr[s]=cdb->getAnodeGain(rn,s);
		for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d) {
			cathsegs[s][d] = cdb->getCathSegCalibrators(rn,s,d);
			if(cathsegs[s][d].size() <= 1) {
				SMExcept e("missingCathodeCalibrators");
				e.insert("run",rn);
				e.insert("side",sideWords(s));
				e.insert("plane",d==X_DIRECTION?"X":"Y");
				throw(e);
			}
			domains[s][d].push_back(0);
			for(unsigned int i=0; i<cathsegs[s][d].size(); i++) {
				wirePos[s][d].push_back(cathsegs[s][d][i]->pos);
				if(i)
					domains[s][d].push_back(0.5*(wirePos[s][d][i]+wirePos[s][d][i-1]));
			}
			domains[s][d][0] = 2*wirePos[s][d][0]-domains[s][d][1];
			domains[s][d].push_back(2*wirePos[s][d].back()-domains[s][d].back());
		}
	}
}

WirechamberCalibrator::~WirechamberCalibrator() {
	for(Side s = EAST; s <= WEST; ++s)
		for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d)
			for(unsigned int i=0; i<cathsegs[s][d].size(); i++)
				delete(cathsegs[s][d][i]);
}

float WirechamberCalibrator::wirechamberGainCorr(Side s, float) const {  
	assert(s==EAST||s==WEST);
	return anodeGainCorr[s];
}

float WirechamberCalibrator::calibrateAnode(float adc, Side s, float x, float y, float t) const {
	assert(s==EAST||s==WEST);
	return adc*wirechamberGainCorr(s,t)/anodeP->eval(s,0,x,y,true);
}

void WirechamberCalibrator::tweakPosition(Side s, AxisDirection d, wireHit& h, double E) const {
	if(h.rawCenter == kUndefinedPosition) return;
	unsigned int n;
	float c;
	toLocal(s,d,h.rawCenter,n,c);
	h.center = fromLocal(s,d,n,cathsegs[s][d][n]->adjustPos(c,E));
}

std::vector<std::string> WirechamberCalibrator::getCathChans(Side s, AxisDirection d) const {
	assert(s<=WEST && d<=Y_DIRECTION);
	std::vector<std::string> chans;
	for(std::vector<CathSegCalibrator*>::const_iterator it = cathsegs[s][d].begin(); it != cathsegs[s][d].end(); it++)
		chans.push_back((*it)->channel);
	return chans;
}

float WirechamberCalibrator::sep23Cut(Side, float Escint) {
	// magic numbers formula from simulation
	//return 2.68 + 4.17*exp(-Escint/146);		// 20120810 baseline simulation
	//return 3.69 + 4.80*exp(-Escint/135.3);	// 20120822 small cuts + 50% dead contrib
	//return 4.16 + 5.94*exp(-Escint/125.6);	// 20120823 +5mm window bowing
	return 3.99 + 5.39*exp(-Escint/148.8); 		// 20120824 MagF
}

float WirechamberCalibrator::normMWPC(Side s, float Escint, float Emwpc) {
	double c = sep23Cut(s,Escint); // hard cut location
	return (Emwpc-c)/c; // normalized MWPC
}

float WirechamberCalibrator::sep23Prob(Side s, float Escint, float Emwpc) {
	double m = normMWPC(s,Escint,Emwpc);
	
	// asymptotic value
	double asympt = ((m<0)?
					 0.053+0.000301*Escint+1.27e-07*Escint*Escint :
					 0.962-0.285*exp(-Escint/149.5) ); // 20120824 MagF
					 //0.0456+3.138e-4*Escint+8.240e-8*Escint*Escint :
					 //0.9505-0.2827*exp(-Escint/132.3) );
	// falloff scale
	double m0 = (m<0)? -0.15 : 0.20;
	
	return asympt+(0.5-asympt)*exp(-m/m0);
}

void WirechamberCalibrator::printSummary() {
	printf("Wirechamber Calibrator for %i,%i, %i,%i cathodes\n",
		   (int)cathsegs[EAST][X_DIRECTION].size(),(int)cathsegs[EAST][Y_DIRECTION].size(),
		   (int)cathsegs[WEST][X_DIRECTION].size(),(int)cathsegs[WEST][Y_DIRECTION].size());
	for(Side s = EAST; s <= WEST; ++s) {
		for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d) {
			printf("%s-%s ",sideWords(s),d==X_DIRECTION?"x":"y");
			for(unsigned int i=0; i<cathsegs[s][d].size(); i++)
				printf("[%+.1f,%i]",(cathsegs[s][d][i]->norm-1)*100,(int)cathsegs[s][d][i]->pcoeffs.size());
			printf("\n");
		}
	}
	printf("Anode 1000:\tE=%.2f\tW=%.2f\n",calibrateAnode(1000,EAST,0,0,0),calibrateAnode(1000,WEST,0,0,0));
}

void WirechamberCalibrator::toLocal(Side s, AxisDirection d, float x, unsigned int& n, float& c) const {
	if(x<=domains[s][d][0]) n=0;
	else if(x>=domains[s][d].back()) n = wirePos[s][d].size()-1;
	else n = std::upper_bound(domains[s][d].begin(),domains[s][d].end(),x)-domains[s][d].begin()-1;
	c = (x-domains[s][d][n])/(domains[s][d][n+1]-domains[s][d][n])-0.5;
}

float WirechamberCalibrator::fromLocal(Side s, AxisDirection d, unsigned int n, float c) const {
	assert(n>=0 && n<=wirePos[s][d].size()-1);
	return domains[s][d][n]*(0.5-c)+domains[s][d][n+1]*(0.5+c);
}

float WirechamberCalibrator::getCathNorm(Side s, AxisDirection d, unsigned int c) const {
	assert(s<=WEST && d<=Y_DIRECTION && c<cathsegs[s][d].size());
	return cathsegs[s][d][c]->norm;
}

Stringmap WirechamberCalibrator::wirecalSummary() const {
	Stringmap m;
	for(Side s = EAST; s <= WEST; ++s) {
		m.insert(sideSubst("anodegain_%c",s),anodeGainCorr[s]);
		for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d) {
			std::string pname = sideSubst("%c",s)+(d==X_DIRECTION?"x":"y");
			std::vector<double> cnorm;
			for(unsigned int c = 0; c < cathsegs[s][d].size(); c++) cnorm.push_back(cathsegs[s][d][c]->norm);
			m.insert("cnorm_"+pname,vtos(cnorm));
		}
	}
	return m;
}

wireHit WirechamberCalibrator::calcHitPos(Side s, AxisDirection d,
										  std::vector<float>& wireValues, std::vector<float>& wirePeds) const {
	
	assert(s<=WEST && d<=Y_DIRECTION);
	const unsigned int nWires = wirePos[s][d].size();
	assert(wireValues.size()>=nWires);
	float x1,x2,x3,y1,y2,y3;
	const double sigma0 = 0.75;	//< typical event width in sire spacings
	
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
	
	h.center = h.rawCenter = kUndefinedPosition;
	
	for(unsigned int c=0; c<nWires; c++) {
		// values above 3950 count as "clipped"
		bool isClipped = wireValues[c] > 3950;
		// pedestal subtract wire readout
		if(wirePeds.size()>=nWires)
			wireValues[c] -= wirePeds[c];
		// cathode normalization
		float cnm = cathsegs[s][d][c]->norm;
		if(isClipped) {
			wireValues[c] = 100000;
			h.nClipped++;
		} else {
			// record in usable wires
			xs.push_back(wirePos[s][d][c]);
			ys.push_back(wireValues[c]*cnm);
		}
		// values above 70 count towards multiplicity
		if(wireValues[c]>70)
			h.multiplicity++;
		h.cathodeSum += wireValues[c]*cnm;
		// check if this is the maximum value found so far
		if(wireValues[c]*cnm >= h.maxValue && !isClipped) {
			// if wires are tied for max value, randomly choose which is labeled as maxWire
			if(wireValues[c]*cnm == h.maxValue && rand()%2)
				continue;
			h.maxValue = wireValues[c]*cnm;
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
	
	
	//---------------
	// special cases: edge wires
	//---------------
	
	// no usable wires? There's no hit to reconstruct.
	if(h.maxValue <= 0 || !xs.size()) {
		h.center = h.width = 0;
		return h;
	}
	// only one usable wire?
	if(xs.size()==1) {
		h.errflags |= WIRES_SINGLET;
		h.rawCenter = h.center = xs[0];
		return h;
	}
	// edge wire; reconstruct assuming fixed sigma0
	if(h.errflags & WIRES_EDGE) {
		y1 = ys[maxn];
		x1 = xs[maxn];
		if(maxn==0)
			maxn++;
		else
			maxn--;
		y2 = ys[maxn];
		x2 = xs[maxn];
				
		// next wire over is negative:
		// isolated edge wire firing
		if(y2 <= 0) {
			h.errflags |= WIRES_SINGLET;
			y2 = 1;
		}
		
		h.width = sigma0*fabs(x2-x1);
		float l = 0.5 + (log(y2)-log(y1))*sigma0*sigma0;
		h.rawCenter = h.center = (1-l)*x1 + l*x2;
		return h;
	}
	
	//---------------
	// interior wires
	//---------------
	
	// collect x and y values
	y1 = ys[maxn-1];
	y2 = ys[maxn];
	y3 = ys[maxn+1];
	x1 = xs[maxn-1];
	x2 = xs[maxn];
	x3 = xs[maxn+1];
	float dx1 = x2-x1;
	float dx2 = x3-x2;
	
	// isolated singlet: wires on both sides negative
	if(y1 <= 0 && y3 <= 0) {
		h.errflags |= WIRES_SINGLET;
		h.width = 0;
		h.rawCenter = h.center = x2;
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
		h.rawCenter = h.center = (1-l)*x2 + l*x3;
		return h;
	}
	
	// normal case: uniformly spaced wires; convert from parabola center to position based on fixed sigma0
	if(fabs(dx1-dx2)<0.001) {
		float x = (y3-y1)/(4.0*y2-2.0*(y1+y3));
		float sigma2 = sigma0*sigma0;
		float gx = (x<0?-1:x==0?0:1)*sigma2*log(( 2.0*exp(1.0/(2.0*sigma2))*fabs(x) + sqrt( 1.0 + 4.0*( exp(1.0/sigma2) - 1.0 )*x*x ) )/(1.0+2.0*fabs(x)));
		h.rawCenter = h.center = x2+dx1*gx;
	} else {
		h.errflags |= WIRES_NONUNIF;
	}
	
	// calculate gaussian center for non-uniform wires, width for all wires
	y1 = log(y1);
	y2 = log(y2);
	y3 = log(y3);
	float denom = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2);
	h.width = sqrt((x1-x2)*(x2-x3)*(x3-x1)/(-2*denom));
	float gcenter = 0.5*(x1*x1*(y2-y3)  + x2*x2*(y3-y1) + x3*x3*(y1-y2)) / denom;
	if(h.center == kUndefinedPosition)
		h.rawCenter = h.center = gcenter;
	
	return h;
}

