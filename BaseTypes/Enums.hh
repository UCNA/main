/// \file Enums.hh \brief various enumerations
#ifndef ENUMS_HH
/// make sure this file is only included once
#define ENUMS_HH 1

#include <math.h>
#include <string>

#define nBetaTubes 4
#define PI 3.1415926535

/// detector side specification
enum Side {
	EAST = 0,	//< East side of detector
	WEST = 1,	//< West side of detector
	BOTH = 2,	//< Both sides of detector
	NONE = 3,	//< Neither side of detector
	BADSIDE = 4
};

/// iteration to next side
inline Side& operator++(Side& s) { return s = Side(s+1); }

/// run number type
typedef unsigned int RunNum;
/// get letter for side names
char sideNames(Side s);
/// get full word for side names
const char* sideWords(Side s);
/// substitute side names into string expression
std::string sideSubst(std::string instr, Side s);

/// opposite side to given side
Side otherSide(Side s);

/// directions, to clean up wireplane direction confusion
enum AxisDirection {
	X_DIRECTION = 0,
	Y_DIRECTION = 1,
	Z_DIRECTION = 2,
	T_DIRECTION = 3
};

/// float value with an error estimate
struct float_err {
	float_err(float c=0, float dc=0): x(c), err(dc) {}
	float x;
	float err;
	//const float_err& operator=(const float_err& rhs) { x = rhs.x; err = rhs.err; return *this; }
};

/// add float_errs assuming independent statistics
float_err operator+(float_err a, float_err b);
/// multiply float_err by float
float_err operator*(float a, float_err b);
/// statistically weighted sum of n values with errors, useful for combining PMT results
float_err weightedSum(unsigned int n, const float_err* d);
/// measure of combined statistical proximity of points d to central value c
float proximity(unsigned int n, const float_err* d, float_err c);

/// types of UCNA data runs
enum RunType {
	UNKNOWN_RUN,	//< unknown/invalid run
	ASYMMETRY_RUN,	//< background or beta decay asymmetry run
	DEPOL_RUN,		//< depolarization run
	LED_RUN,		//< LED linearity run
	SOURCE_RUN,		//< source calibration run
	XE_RUN			//< activated Xenon calibration run
};

/// octet type specification
/// 0 -> not in octet
/// 1...12 -> A1...A12
/// 13...24 -> B1...B12
typedef unsigned int OctetType;

/// types of beta/background/depol triads
enum TriadType {
	TRIAD_NONE		= 0,		//< not part of a triad
	TRIAD_A010203	= 1,	//< A1-A3 triad
	TRIAD_A040506	= 2,	//< A4-A6 triad
	TRIAD_A070809	= 3,	//< A7-A9 triad
	TRIAD_A101112	= 4,	//< A10-A12 triad
	TRIAD_B010203	= 5,	//< B1-B3 triad
	TRIAD_B040506	= 6,	//< B4-B6 triad
	TRIAD_B070809	= 7,	//< B7-B9 triad
	TRIAD_B101112	= 8		//< B10-B12 triad
};

/// state of AFP during run
enum AFPState {
	AFP_OFF		=	0,	//< AFP is off
	AFP_ON		=	1,	//< AFP is on
	AFP_OTHER	=	2,	//< AFP state unspecified
	AFP_OFF2ON	=	3,	//< (Depolarization run) AFP off->on
	AFP_ON2OFF	=	4	//< (Depolarization run) AFP on->off
};

/// state of gate valve during a run
enum GVState {
	GV_CLOSED	= 0,	//< GV closed (background)
	GV_OPEN		= 1,	//< GV open (beta decay)
	GV_OTHER	= 2		//< GV state unspecified
};

/// iteration to next GV state
inline GVState& operator++(GVState& a) { return a = GVState(a+1); }

/// SCS geometry for a run
enum RunGeometry {
	GEOMETRY_OTHER,	//< unspecified geometry
	GEOMETRY_2007,	//< 2007 Geometry: 2.5u mylar, 0.3u Be, 25u wirechamber windows
	GEOMETRY_A,		//< 'A' 0.7u mylar, 0.3u Be, 25u wirechamber windows
	GEOMETRY_B,		//< 'B' 0.7+12.5u mylar, 0.3u Be, 25u wirechamber windows
	GEOMETRY_C,		//< 'C' 0.7u mylar, 0.3u Be, 6u wirechamber windows
	GEOMETRY_D		//< 'D' No endcaps, 6u wirechamber windows
};

enum DataQuality {
	DQ_UNKNOWN		=	0,		//< unknown data quality
	DQ_CHECKED		=	1<<0,	//< data quality has been checked
	DQ_TRASH		=	1<<1,	//< run is garbage, DO NOT USE
	DQ_CORRUPT_END	=	1<<2,	//< data at end of run is corrupted
	DQ_BEAM_ISSUES	=	1<<3,	//< beam quality dubious during run
};

/// get the geometry for a given run
RunGeometry whichGeometry(RunNum rn);
/// get the name (letter) for a run geometry
std::string geomName(RunGeometry g);
/// get hardware PMT number for given side/quadrant
unsigned int pmtHardwareNum(Side s, unsigned int quadrant);

/// enumeration for Jianglai asymmetry analysis choices
enum AnalysisChoice {
	ANCHOICE_NONE=0,
	ANCHOICE_A,		//< Assign 0, I, II+III (unseparated) to primary trigger side
	ANCHOICE_B,		//< Only use 0, I
	ANCHOICE_C,		//< Separate II/III with 4keV MWPC cut
	ANCHOICE_D,		//< Only use 0
	ANCHOICE_E,		//< Separate II/III with likelihood
	ANCHOICE_F,
	ANCHOICE_G
};

/// enumeration for Jianglai particle ID
enum PID {
	PID_SINGLE = 0,	//< "single" gamma event
	PID_BETA = 1,	//< beta event
	PID_MUON = 2,	//< muon event
	PID_LED = 3,		//< tagged LED event
	PID_PULSER = 4	//< tagged Bi pulser event
};

//--------------------------------------------------------------------------
// 1	"Type I" -- Both scintillators and wirechambers
// 2	"Type II" -- Scatters off wirechamber, goes through other side
// 3	"Type III" -- Scatters off scintillator, dies in other wirechamber
// 5	"Type IV" -- Undetectable scatters off wirechamber surface
//--------------------------------------------------------------------------

/// event type enum
enum EventType {
	TYPE_0_EVENT	= 0,
	TYPE_I_EVENT	= 1,
	TYPE_II_EVENT	= 2,
	TYPE_III_EVENT	= 3,
	TYPE_IV_EVENT	= 4
};

std::string typeWords(EventType tp);

#endif
