#ifndef INTERPOLATOR_HH
#define INTERPOLATOR_HH 1

#include <math.h>
#include <vector>
#include <cassert>
#include <stdio.h>

/// boundary conditions for interpolation
enum BoundaryCondition {
	BC_INFINITE,		//< extend endpoint values out to infinity
	BC_CYCLIC,			//< cyclic boundary conditions
	BC_DERIVCLAMP_ZERO	//< BC for clamping derivative to 0 at 0 for bicubic interpolation
};

/// infinite, bi-directional sequence of points
class DataSequence {
public:
	/// constructor
	DataSequence(BoundaryCondition b = BC_CYCLIC): npts(0), bc(b) {}
	/// destructor
	virtual ~DataSequence() {}
	/// data retreival
	virtual double valueAt(int i, void* xopts) = 0;
	/// get number of internal points
	virtual int getNpts() { return npts; }
	/// display state
	virtual void display() const { printf("<Generic DataSequence>\n"); }
protected:
	/// coerce data point into actual data range
	virtual unsigned int coerce(int i) {
		if(bc == BC_CYCLIC) {
			if(i>=0)
				return i%npts;
			return ((i%npts)+npts)%npts;
		} else if(bc == BC_INFINITE) {
			if(i<=0) return 0;
			if(i>=npts) return npts-1;
			return (unsigned int)i;
		} else if (bc == BC_DERIVCLAMP_ZERO) {
			assert(npts>=2);
			if(i<0) return 1;
			if(i>=npts) return npts-1;
			return (unsigned int)i;
		}
		return 0;
	}
	int npts;				//< number of internal gridpoints
	BoundaryCondition bc;	//< boundary conditions for generating sequence
};

/// data sequence based on internal array of doubles
class DoubleSequence: public DataSequence {
public:
	DoubleSequence(BoundaryCondition b = BC_CYCLIC): DataSequence(b) {}
	/// value at given grid location
	virtual double valueAt(int i, void*) { return pts[coerce(i)]; }
	/// add data point
	void addPoint(double p) { pts.push_back(p); npts++; }
	/// display state
	virtual void display() const {
		printf("<DoubleSequence> "); 
		for(unsigned int i = 0; i<pts.size(); i++)
			printf("%f, ",pts[i]);
		printf("(%i)\n",npts);
	}
protected:	
	std::vector<double> pts;	//< internal list of points
};

/// generic interpolator for intermediate points in a sequence
class Interpolator {
public:
	/// constructor, with input scale factor s, offset o (where interpolator '0' point is placed on input axis)
	Interpolator(DataSequence* L, double s = 1.0, double o = 0.0): myData(L), scale(s), offset(o) {}
	/// destructor
	virtual ~Interpolator() {}
	/// subclass interpolation method here; nearest-neighbor interpolation used by default
	virtual double eval(double* x) { 
		double y;
		int i = locate(*x,&y);
		if(y<=0.5)
			return myData->valueAt(i,(void*)(x+1));
		return myData->valueAt(i+1,(void*)(x+1));
	}
protected:
	/// locate position in data coordinates, remainder in [0,1)
	virtual int locate(double x, double* remainder = NULL) {
		double l = (x-offset)*myData->getNpts()/scale;
		if(remainder)
			*remainder = l-floor(l);
		return int(floor(l));
	}
	DataSequence* myData;	//< sequence to be interpolated
	double scale;			//< internal length scale
	double offset;			//< zero point coordinate offset
};

/// linear interpolator
class LinTerpolator: public Interpolator {
public:
	/// constructor, with scale factor s, offset o
	LinTerpolator(DataSequence* L, double s = 1.0, double o = 0.0): Interpolator(L,s,o) {}
	/// evaluation by linear interpolation
	virtual double eval(double* x) {
		double y;
		int i = locate(*x,&y);
		return myData->valueAt(i,(void*)(x+1))*(1-y)+myData->valueAt(i+1,(void*)(x+1))*y;
	}
};

/// cubic interpolator
class CubiTerpolator: public Interpolator {
public:
	/// constructor, with scale factor s, offset o, 'sharpening' a
	CubiTerpolator(DataSequence* L, double s = 1.0, double o = 0.0, double a = -0.5): Interpolator(L,s,o), A(a) { }
	/// evaluation by linear interpolation
	virtual double eval(double* x) {
		double y;
		int i = locate(*x,&y);
		double p0 = myData->valueAt(i-1,(void*)(x+1));
		double p1 = myData->valueAt(i,(void*)(x+1));
		double p2 = myData->valueAt(i+1,(void*)(x+1));
		double p3 = myData->valueAt(i+2,(void*)(x+1));
		return ( A*p0*(1-y)*(1-y)*y
				+p1*(1-y)*(1-y*((2+A)*y-1))
				-p2*y*(A*(1-y)*(1-y)+y*(2*y-3))
				+A*p3*(1-y)*y*y );

	}
protected:
	double A;	//< "sharpening" coefficient, default = -0.5
};

/// a sequence of interpolators for multi-dimensional interpolation between interpolators
class InterpoSequence: public DataSequence {
public:
	/// constructor
	InterpoSequence(BoundaryCondition b = BC_CYCLIC): DataSequence(b) {}
	/// add data point
	void addPoint(Interpolator* i) { myInterpolators.push_back(i); npts++; }
	/// data retreival
	virtual double valueAt(int i, void* xopts) { return myInterpolators[coerce(i)]->eval((double*)xopts); }
protected:
	std::vector<Interpolator*> myInterpolators;	//< interpolators in sequence
};


#endif
