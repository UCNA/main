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

/// run number type
typedef unsigned int RunNum;
/// get letter for side names
char sideNames(Side s);
/// get full word for side names
const char* sideWords(Side s);
/// substitute side names into string expression
std::string sideSubst(std::string instr, Side s);

/// combined Side
Side sideCombo(bool E, bool W);
/// iterate to next Side
Side nextSide(Side s);
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
	GV_CLOSED,	//< GV closed (background)
	GV_OPEN,	//< GV open (beta decay)
	GV_OTHER	//< GV state unspecified
};

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


/// basic trigger data
enum TrigFlags {
	TRIG_E1_TDC = 1<<0,		//< E1 TDC not timed out
	TRIG_E2_TDC = 1<<1,		//< E2 TDC not timed out
	TRIG_E3_TDC = 1<<2,		//< E3 TDC not timed out
	TRIG_E4_TDC = 1<<3,		//< E4 TDC not timed out
	TRIG_W1_TDC = 1<<4,		//< W1 TDC not timed out
	TRIG_W2_TDC = 1<<5,		//< W2 TDC not timed out
	TRIG_W3_TDC = 1<<6,		//< W3 TDC not timed out
	TRIG_W4_TDC = 1<<7,		//< W4 TDC not timed out	
};

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

/// event id flag bits
enum EventFlags {
	
	IS_CORRECT			= 1<<0, //<	correct or "Type IV" -- Undetectable scatters off wirechamber surface
	IS_TYPE_I			= 1<<1, //<	Type I hits scintillators and WCs on both sides
	IS_TYPE_II			= 1<<2, //<	Type II scatters off wirechambers, ends up in other scintillator
	IS_TYPE_III			= 1<<3,	//< Type III scatters off scintillator, ends up in other wirechamber
	
	IS_BEAMCRUD			= 1<<4, //<	Beam Crud (too close to beam pulse)
	IS_GMS				= 1<<5, //<	GMS LED event
	IS_MUON				= 1<<6, //<	tagged muon background event
	IS_UCN_MON			= 1<<7,	//<	UCN monitor trigger event
	IS_PMT_PULSER		= 1<<8,	//<	Individual PMT pulser event
	
	HIT_WIRES			= 1<<9, //<	tracked in MWPC on this side
	HIT_SCINT			= 1<<10,//<	seen in scintillator
	HIT_MUON_BACK		= 1<<11,//< hit muon backing veto on this side
	HIT_SIDE_FIRST		= 1<<12,//<	whether this side was hit first (in TDC zero-time peak or not)
	
	PRIMARY_DIRECTION	= 1<<13,//< whether the (simulation) particle was originally heading in this direction
	
	IS_NOT_BETADATA		= IS_BEAMCRUD | IS_GMS | IS_MUON | IS_UCN_MON | IS_PMT_PULSER,
	IS_CLASSIFIED_HIT	= IS_CORRECT | IS_TYPE_I | IS_TYPE_II | IS_TYPE_III
};

#endif
