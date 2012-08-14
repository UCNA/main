#ifndef SURFACEGNERATOR_HH
#define SURFACEGNERATOR_HH 1

#include <G4ThreeVector.hh>
#include <vector>
#include <cmath>

/// base class for portion of a surface
class SurfaceSeg {
public:
	/// constructor
	SurfaceSeg(const G4ThreeVector& O): o(O) {}
	/// destructor
	virtual ~SurfaceSeg() {}
	
	/// get surface area
	virtual double getArea() const = 0;
	
	/// get random point on surface
	G4ThreeVector getSurfRandom() { return o+_getSurfRandom(); }
	
	/// check whether a direction vector from a given point intersects this surface
	bool intersectsSurface(G4ThreeVector p0, G4ThreeVector d) const { return _intersectsSurface(p0-o,d); }
	
	G4ThreeVector snorm;	//< surface normal for most recently generated point
	
protected:
	/// get random point in local coordinates
	virtual G4ThreeVector _getSurfRandom() = 0;
	/// check whether a vector from a given point, in local coordinates, intersects this surface
	virtual bool _intersectsSurface(G4ThreeVector p0, G4ThreeVector d) const = 0;
		
	G4ThreeVector o;	//< center positioning of surface
};

/// conical frustum surface
class ConeFrustum: public SurfaceSeg {
public:
	/// constructor
	ConeFrustum(const G4ThreeVector& O, double R0, double R1, double DZ, bool outn = true);
	
	/// get surface area
	virtual double getArea() const { return M_PI*(r0+r1)*sqrt((r0-r1)*(r0-r1)+dz*dz); }
	
protected:
	/// get random point in local coordinates
	virtual G4ThreeVector _getSurfRandom();
	/// check whether a vector from a given point, in local coordinates, intersects this surface
	virtual bool _intersectsSurface(G4ThreeVector p0, G4ThreeVector d) const;
	
	bool outwardNormal;
	double r0,r1,dz,sn0;
};

/// (weighted) assembly of surfaces
class SurfaceAssembly: public SurfaceSeg {
public:
	/// constructor
	SurfaceAssembly(const G4ThreeVector& O): SurfaceSeg(O) { cumweights.push_back(0); }
	
	/// get surface area
	virtual double getArea() const { return cumweights.back(); }
	/// add sub-segment
	void addSegment(SurfaceSeg* s, double w = 1.0);
	
protected:
	/// get random point in local coordinates
	virtual G4ThreeVector _getSurfRandom();
	/// check whether a vector from a given point, in local coordinates, intersects this surface
	virtual bool _intersectsSurface(G4ThreeVector p0, G4ThreeVector d) const;
	
	std::vector<SurfaceSeg*> subsegs;	//< sub-segments
	std::vector<double> weights;		//< weights assigned to each subsegment
	std::vector<double> cumweights;		//< cumulative weights with areas for each segment
};


#endif
