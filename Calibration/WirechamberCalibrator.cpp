#include "WirechamberCalibrator.hh"
#include "SMExcept.hh"
#include "GraphicsUtils.hh"
#include <cmath>
#include "SMExcept.hh"
#include <algorithm>

#define kUndefinedPosition 666

bool WirechamberCalibrator::calibrateCathodes = true;

WirechamberCalibrator::WirechamberCalibrator(RunNum rn, CalDB* cdb): sigma(5.90) {
	for(Side s = EAST; s <= WEST; ++s) {
		
		// load energy calibration methods
		std::vector<MWPC_Ecal_Spec> ecalList = cdb->get_MWPC_Ecals(rn,s);
		smassert(ecalList.size());
		for(std::vector<MWPC_Ecal_Spec>::iterator it = ecalList.begin(); it != ecalList.end(); it++) {
			smassert(it->pcorr);
			it->pcorr->setNormAvg();
			ecalMethods[s].insert(std::pair<ChargeProxyType,MWPC_Ecal_Spec>(it->charge_meas,*it));
		}
		
		// set primary energy calibration method
		MWPC_Ecal_Spec mySpec = ecalList[0];
		mwpcGainCorr[s] = mySpec.gain_factor;
		myChargeProxy[s] = mySpec.charge_meas;
		chargeP[s] = mySpec.pcorr;
		
		// organize other calibration methods, synthesizing approximate stand-ins where needed
		for(ChargeProxyType c = CHARGE_PROXY_ANODE; c <= CHARGE_PROXY_CCLOUD; ++c) {
			if(ecalMethods[s].find(c) != ecalMethods[s].end()) continue;
			MWPC_Ecal_Spec sp = mySpec;
			sp.charge_meas = c;
			sp.gain_factor = (c==CHARGE_PROXY_ANODE?0.008:.00022);
			ecalMethods[s].insert(std::pair<ChargeProxyType,MWPC_Ecal_Spec>(c,sp));
		}
		
		MWPC_Ecal_Spec ccalSpec = getAltEcal(s,CHARGE_PROXY_CCLOUD);
		ccloudGainCorr[s] = ccalSpec.gain_factor;
		ccloud_eta[s] = ccalSpec.pcorr;
		
		for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d) {
			cath_ccloud_gains[s][d] = cdb->getCathCCloudGains(rn,s,d);
			cathsegs[s][d] = cdb->getCathSegCalibrators(rn,s,d);
			nCaths[s][d] = (unsigned int)cathsegs[s][d].size();
			if(nCaths[s][d] <= 1 || cath_ccloud_gains[s][d].size() < nCaths[s][d]) {
				SMExcept e("missingCathodeCalibrators");
				e.insert("run",rn);
				e.insert("side",sideWords(s));
				e.insert("plane",d==X_DIRECTION?"X":"Y");
				throw(e);
			}
			domains[s][d].push_back(0);
	
 			for(unsigned int i=0; i<nCaths[s][d]; i++) {
				cathNames[s][d].push_back(sideSubst("MWPC%c",s)+(d==X_DIRECTION?"x":"y")+itos(i+1));
				cathseg_energy_norm[s][d].push_back(cath_ccloud_gains[s][d][i]/ccloudGainCorr[s]/cathsegs[s][d][i]->norm);
				
				if(!calibrateCathodes) {
					cathsegs[s][d][i]->norm = 1.0;
					cathsegs[s][d][i]->pcoeffs.clear();
				}
				wirePos[s][d].push_back(cathsegs[s][d][i]->pos);
				if(i)
					domains[s][d].push_back(0.5*(wirePos[s][d][i]+wirePos[s][d][i-1]));
			}
			domains[s][d][0] = 2*wirePos[s][d][0]-domains[s][d][1];
			domains[s][d].push_back(2*wirePos[s][d].back()-domains[s][d].back());
		}
		
		// load cuts
		loadRangeCut(rn,fCathMax[s],sideSubst("Cut_MWPC_%c_CathMax",s));
		loadRangeCut(rn,fCathSum[s],sideSubst("Cut_MWPC_%c_CathSum",s));
		loadRangeCut(rn,fCathMaxSum[s],sideSubst("Cut_MWPC_%c_CathMaxSum",s));
	}
}

WirechamberCalibrator::~WirechamberCalibrator() {
	for(Side s = EAST; s <= WEST; ++s)
		for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d)
			for(unsigned int i=0; i<cathsegs[s][d].size(); i++)
				delete(cathsegs[s][d][i]);
}

float WirechamberCalibrator::chargeProxy(Side s, ChargeProxyType c, const wireHit& x_wires, const wireHit& y_wires, const MWPCevent& mwpc) const {
	if(c==CHARGE_PROXY_ANODE) return mwpc.anode;
	if(c==CHARGE_PROXY_CCLOUD) return x_wires.ccloud_size()+y_wires.ccloud_size();
	return 0;
}

float WirechamberCalibrator::wirechamberEnergy(Side s, const wireHit& x_wires, const wireHit& y_wires, const MWPCevent& mwpc) const {
	smassert(s==EAST||s==WEST);
	float eta_w = chargeP[s]->eval(s,0,x_wires.center, y_wires.center,true);
	float Q = chargeProxy(s,myChargeProxy[s],x_wires,y_wires,mwpc);
	return mwpcGainCorr[s]*Q/eta_w;
}

const MWPC_Ecal_Spec& WirechamberCalibrator::getAltEcal(Side s, ChargeProxyType tp) const {
	std::map<ChargeProxyType,MWPC_Ecal_Spec>::const_iterator it = ecalMethods[s].find(tp);
	smassert(it != ecalMethods[s].end());
	return it->second;
}


void WirechamberCalibrator::tweakPosition(Side s, AxisDirection d, wireHit& h, double E) const {
	if(h.rawCenter == kUndefinedPosition) return;
	unsigned int n;
	float c;
	toLocal(s,d,h.rawCenter,n,c);
	h.center = fromLocal(s,d,n,cathsegs[s][d][n]->adjustPos(c,E));
}

std::vector<std::string> WirechamberCalibrator::getCathChans(Side s, AxisDirection d) const {
	smassert(s<=WEST && d<=Y_DIRECTION);
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
		   nCaths[EAST][X_DIRECTION],nCaths[EAST][Y_DIRECTION],
		   nCaths[WEST][X_DIRECTION],nCaths[WEST][Y_DIRECTION]);
	for(Side s = EAST; s <= WEST; ++s) {
		for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d) {
			printf("%s-%s ",sideWords(s),d==X_DIRECTION?"x":"y");
			for(unsigned int i=0; i<nCaths[s][d]; i++)
				printf("[%+.1f,%i]",(cathsegs[s][d][i]->norm-1)*100,(int)cathsegs[s][d][i]->pcoeffs.size());
			printf("\n\t");
			for(unsigned int i=0; i<nCaths[s][d]; i++)
				printf("[%03i]",(int)cathPeds0[s][d][i]);
			printf(" ped\n\t");
			for(unsigned int i=0; i<nCaths[s][d]; i++)
				printf("[%03i]",(int)cathPedW0[s][d][i]);
			printf(" pedw\n\t");
			for(unsigned int i=0; i<nCaths[s][d]; i++)
				printf("[%03i]",(int)cathseg_energy_norm[s][d][i]);
			printf(" enorm\n");
		}
		MWPC_Ecal_Spec mySpec = getAltEcal(s,myChargeProxy[s]);
		printf("Energy via %s, gain factor %.4g (runs %i-%i)\n",chargeProxyName(myChargeProxy[s]).c_str(),mySpec.gain_factor,mySpec.start_run,mySpec.end_run);
	}
}

void WirechamberCalibrator::toLocal(Side s, AxisDirection d, float x, unsigned int& n, float& c) const {
	if(x<=domains[s][d][0]) n=0;
	else if(x>=domains[s][d].back()) n = wirePos[s][d].size()-1;
	else n = std::upper_bound(domains[s][d].begin(),domains[s][d].end(),x)-domains[s][d].begin()-1;
	c = (x-domains[s][d][n])/(domains[s][d][n+1]-domains[s][d][n])-0.5;
}

float WirechamberCalibrator::fromLocal(Side s, AxisDirection d, unsigned int n, float c) const {
	smassert(n<=wirePos[s][d].size()-1);
	return domains[s][d][n]*(0.5-c)+domains[s][d][n+1]*(0.5+c);
}

float WirechamberCalibrator::getCathNorm(Side s, AxisDirection d, unsigned int c) const {
	smassert(s<=WEST && d<=Y_DIRECTION && c<nCaths[s][d]);
	return cathsegs[s][d][c]->norm;
}

float WirechamberCalibrator::getCathCCloudGain(Side s, AxisDirection d, unsigned int c) const {
	smassert(s<=WEST && d<=Y_DIRECTION);
	smassert(c<nCaths[s][d]);
	return cath_ccloud_gains[s][d][c];
}

Stringmap WirechamberCalibrator::wirecalSummary() const {
	Stringmap m;
	for(Side s = EAST; s <= WEST; ++s) {
		m.insert(sideSubst("mwpc_gain_%c",s),mwpcGainCorr[s]);
		m.insert(sideSubst("mwpc_calmethod_%c",s),chargeProxyName(myChargeProxy[s]));
		for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d) {
			std::string pname = sideSubst("%c",s)+(d==X_DIRECTION?"x":"y");
			std::vector<double> cnorm;
			for(unsigned int c = 0; c < nCaths[s][d]; c++) cnorm.push_back(cathsegs[s][d][c]->norm);
			m.insert("cnorm_"+pname,vtos(cnorm));
		}
	}
	return m;
}

void WirechamberCalibrator::calcDoubletHitPos(wireHit& h, float x0, float x1, float y0, float y1) const {
		float sclamp = exp(-2*(x1-x0)*(x1-x0)/(2*sigma*sigma));
		if(y1<y0*sclamp) y1 = y0*sclamp;
		
		h.width = sigma;
		h.rawCenter = h.center = (x0+x1)/2 + sigma*sigma/(x1-x0)*log(y1/y0);
		const double sigma0 = sigma/fabs(x1-x0);
		if(!(h.errflags & WIRES_SINGLET))
			h.height = exp( pow(log(y1/y0)*sigma0,2)/2 + log(y0*y1)/2 + 1/(8*sigma0*sigma0) );
}

wireHit WirechamberCalibrator::calcHitPos(Side s, AxisDirection d, const float* cathADC, const float* cathPed) const {
	
	smassert(s<=WEST && d<=Y_DIRECTION);
	const unsigned int nWires = nCaths[s][d];
	float x1,x2,x3,y1,y2,y3;
	
	std::vector<float> xs;
	std::vector<float> ys;
	int maxn = -1;
	
	// initialize values in h
	wireHit h;
	h.nClipped = 0;
	h.maxWire = 0;
	h.cathodeSum = 0;
	h.multiplicity = 0;
	h.height = 0;
	h.errflags = 0;
	h.maxValue = -1;
	
	h.center = h.rawCenter = kUndefinedPosition;
	
	float wireValues[kMaxCathodes];
	
	for(unsigned int c=0; c<nWires; c++) {
		// values above 3950 count as "clipped"
		bool isClipped = (cathADC[c] + (cathPed?cathPed[c]:0)) > 3950;
		// normalized wire value
		wireValues[c] = cathADC[c] * cathsegs[s][d][c]->norm;
		
		if(isClipped) {
			h.nClipped++;
		} else {
			// record in usable wires
			xs.push_back(wirePos[s][d][c]);
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
			maxn = (int)xs.size()-1;
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
		h.center = h.width = h.height = 0;
		return h;
	}
	// only one usable wire?
	if(xs.size()==1) {
		h.errflags |= WIRES_SINGLET;
		h.rawCenter = h.center = xs[0];
		h.height = ys[0];
		return h;
	}
	// edge wire; reconstruct assuming fixed sigma
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
			h.height = y1;
			y2 = 1;
		}
		
		calcDoubletHitPos(h,x1,x2,y1,y2);
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
		h.height = y2;
		return h;
	}
	
	// isolated doublet (wires on one side negative): fixed sigma reconstruction
	if(y1 <= 0 || y3 <= 0) {
		if(y3 <= 0) {
			y3 = y1;
			x3 = x1;
		}
		h.errflags |= WIRES_DOUBLET;
		calcDoubletHitPos(h,x2,x3,y2,y3);
		return h;
	}
	
	// normal case: uniformly spaced wires; convert from parabola center to position based on fixed sigma0
	if(fabs(dx1-dx2)<0.001) {
		float x = (y3-y1)/(4.0*y2-2.0*(y1+y3));
		float sigma0 = sigma/fabs(dx1);
		float sigma2 = sigma0*sigma0;
		float gx = (x<0?-1:x==0?0:1)*sigma2*log(( 2.0*exp(1.0/(2.0*sigma2))*fabs(x) + sqrt( 1.0 + 4.0*( exp(1.0/sigma2) - 1.0 )*x*x ) )/(1.0+2.0*fabs(x)));
		h.rawCenter = h.center = x2+dx1*gx;
	} else {
		h.errflags |= WIRES_NONUNIF;
	}
	
	// calculate gaussian center for non-uniform wires, width and height for all wires
	y1 = log(y1);
	y2 = log(y2);
	y3 = log(y3);
	float denom = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2);
	h.width = sqrt((x1-x2)*(x2-x3)*(x3-x1)/(-2*denom));
	float gcenter = 0.5*(x1*x1*(y2-y3)  + x2*x2*(y3-y1) + x3*x3*(y1-y2)) / denom;
	h.height = exp( ((x3-x2)*x3*x2*y1 - (x3-x1)*x3*x1*y2 + (x2-x1)*x2*x1*y3)/((x2-x1)*(x3-x1)*(x3-x2)) + gcenter*gcenter/(2*h.width*h.width) );
	if(h.center == kUndefinedPosition)
		h.rawCenter = h.center = gcenter;
	
	return h;
}

void WirechamberCalibrator::drawWires(Side s, AxisDirection p, TVirtualPad* C, Int_t color, AxisDirection onAxis) const {
	smassert(C);
	smassert(s<=WEST && p<=Y_DIRECTION);
	for(std::vector<double>::const_iterator it = wirePos[s][p].begin(); it != wirePos[s][p].end(); it++)
		if(onAxis==X_DIRECTION) drawVLine(*it,C,color);
		else drawHLine(*it,C,color);
}

