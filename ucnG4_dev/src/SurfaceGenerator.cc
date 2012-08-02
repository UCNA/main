#include "SurfaceGenerator.hh"
#include <Randomize.hh>
#include <algorithm>
#include <cmath>

ConeFrustum::ConeFrustum(const G4ThreeVector& O, double R0, double R1, double DZ, bool outn):
SurfaceSeg(O), outwardNormal(outn), r0(R0), r1(R1), dz(fabs(DZ)) {
	double dr = r1-r0;
	double l = sqrt(dr*dr+dz*dz);
	sn0 = dz/l*(outn?1:-1);
	snorm[2] = -dr/l*(outn?1:-1);
}

G4ThreeVector ConeFrustum::_getSurfRandom() {
	double phi = G4UniformRand()*2*M_PI;
	// solve for l: a*l^2 + b*l = c
	// inverting cumulative distribution
	double c = G4UniformRand()*(r1+r0)*0.5;
	double a = (r1-r0)*0.5;
	double l = a?(-r0+sqrt(r0*r0+4*a*c))/(2*a):c/r0;
	double r = r0*(1-l)+r1*l;
	double cp = cos(phi);
	double sp = sin(phi);
	snorm[0] = sn0*cp;
	snorm[1] = sn0*sp;
	return G4ThreeVector(r*cp,r*sp,dz*(l-0.5));
}

bool ConeFrustum::_intersectsSurface(G4ThreeVector p0, G4ThreeVector d) const {
	// special case: constant z direction
	if(!d[2]) {
		if(!(-dz/2. < p0[2] && p0[2] < dz/2)) return false;
		double rp2 = p0[0]*p0[0]+p0[1]*p0[1];
		double l = (p0[2]+dz/2)/dz;
		double rf = (1-l)*r0+l*r1;
		if(rp2 < rf*rf) return true;
		//TODO
	}
		
	// special case: flat disc
	if(!dz) {
		double k = -p0[2]/d[2];
		if(k<0) return false;
		double rp2 = pow(p0[0]+k*d[0],2)+pow(p0[1]+k*d[1],2);
		return rp2 <= r0*r0 ? rp2 > r1*r1 : rp2 <= r1*r1;
	}
		
	// quadratic equation coefficients for intersection z coordinates
	//G4cout "From " << d << "in direction " << p0 << G4endl;
	double a = d.mag2()-pow(d[2]/(dz/2)*(r1-r0),2);
	double b = 2*d[2]*p0.dot(d)-d[2]*d[2]*(r1-r0)*2/dz;
	double c = pow(d[2]*p0[0]-p0[2]*d[0],2)+pow(d[2]*p0[1]-p0[2]*d[1],2)-pow((r0+r1)*d[2]/2,2);
	//G4cout "From " << d << "in direction " << p0 << ": a,b,c = " << a << ", " << b << ", " << c << G4endl;
	
	// quadratic equation solutions
	double delta = b*b-4*a*c;
	if(delta<0) return false;
	delta = sqrt(delta);
	double z1 = (-b-delta)/(2*a);
	double z2 = (-b+delta)/(2*a);
	//G4cout << "  z1,z2 = " << z1 << ", " << z2 << G4cout;
	
	// check if intersections fall within range and in correct direction from p0
	int sgn = d[2]>0?1:-1;
	return (-dz/2. <= z1 && z1 <= dz/2. && z1*sgn > p0[2]*sgn) || (-dz/2. <= z2 && z2 <= dz/2. && z2*sgn > p0[2]*sgn);
}

//-----------------------------------------------------------------

void SurfaceAssembly::addSegment(SurfaceSeg* s, double w) {
	assert(s && s != this);
	G4cout << "Adding subsurface " << cumweights.size() << ": Area " << s->getArea()/(cm*cm) << " cm^2 * " << w << G4endl;
	subsegs.push_back(s);
	weights.push_back(w);
	cumweights.push_back(cumweights.back()+s->getArea()*w);
}
	
G4ThreeVector SurfaceAssembly::SurfaceAssembly::_getSurfRandom() {
	assert(cumweights.back());
	std::vector<double>::const_iterator itsel = std::upper_bound(cumweights.begin(),cumweights.end(),cumweights.back()*G4UniformRand());
	unsigned int selected = int(itsel-cumweights.begin())-1;
	G4ThreeVector p = subsegs[selected]->getSurfRandom();
	snorm = subsegs[selected]->snorm;
	return p;
}

bool SurfaceAssembly::_intersectsSurface(G4ThreeVector p0, G4ThreeVector d) const {
	for(unsigned int i=0; i<subsegs.size(); i++)
		if(subsegs[i]->intersectsSurface(p0,d) && G4UniformRand()<weights[i])
			return true;
	return false;
}

