#ifndef WIRECHAMBERRECONSTRUCTION_HH
#define WIRECHAMBERRECONSTRUCTION_HH 1

#include "Enums.hh"
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <string>

/// enumeration for various types of wirechamber reconstruction errors
enum reconstructionErrors {
	WIRES_GOOD		= 0,	//< perfectly good wirechamber events
	WIRES_NONE		= 1<<0,	//< no wires fired
	WIRES_EDGE		= 1<<1,	//< reconstructed event on edge of wire chamber
	WIRES_SINGLET	= 1<<2,	//< single wire fired
	WIRES_DOUBLET	= 1<<3,	//< 2 wires fired
	WIRES_CLIPPED	= 1<<4,	//< clipping in wirechamber signal
	WIRES_MULTICLIP	= 1<<5,	//< multiple wire clipping
	WIRES_NONUNIF	= 1<<6	//< non-uniform wire spacing required in calculations
};

/// struct for reconstruction information about a wirechamber hit
struct wireHit {
	float center;				//< reconstructed center of charge cloud
	float width;				//< reconstructed 1-sigma width of charge cloud
	float maxValue;				//< maximum wire signal for event
	float cathodeSum;			//< sum of cathode signals for event
	unsigned int maxWire;		//< wire with largest signal
	unsigned int nClipped;		//< number of wires clipped at maximum value
	unsigned int multiplicity;	//< number of wires firing above threshold for this event
	int errflags;				//< reconstruction error warnings
	float gausCenter;			//< center based on 3-pt gaussian fit
	float avgCenter;			//< center based on 3-wire weighted average
	float parabCenter;			//< center of 3-point parabola fit
};

/// get std::vector of which wires are 'live' for a given run/side/plane
std::vector<bool> getLiveWires(RunNum rn, Side s, AxisDirection d);

/// get std::vector of Padc channel numbers for live wires
std::vector<unsigned int> getPadcNumbers(RunNum rn, Side s, AxisDirection d);

/// get standardized "sensor name" for each live wire
std::vector<std::string> getCathodeNames(RunNum rn, Side s, AxisDirection d);

/// get wire positions for given run number, side, and plane
std::vector<float> calcWirePositions(RunNum rn, Side s, AxisDirection d, float wireSpacing = 10.16*sqrt(0.6));

/// reconstruct hit positions based on analytic gaussian model
wireHit mpmGaussianPositioner(const std::vector<float>& wirepos, float* wireValues, const float* wirePeds);

#endif
