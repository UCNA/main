/// \file Enums.hh \brief various enumerations
#ifndef ENUMS_HH
/// make sure this file is only included once
#define ENUMS_HH

#include <math.h>
#include <string>

#define nBetaTubes 4

/// detector side specification
enum Side {
	EAST = 0,	///< East side of detector
	WEST = 1,	///< West side of detector
	BOTH = 2,	///< Both sides of detector
	NOSIDE = 3,	///< Neither side of detector
	BADSIDE = 4 ///< Improperly defined side
};

/// iteration to next side
inline Side& operator++(Side& s) { return s = Side(s+1); }

/// run number type
typedef unsigned int RunNum;
/// get letter for side names
char sideNames(Side s, bool clower = false);
/// get full word for side names
const char* sideWords(Side s);
/// get database enum side names (in single quotes)
const char* dbSideName(Side s);
/// convert string back to side
Side strToSide(const std::string& s);
/// substitute side names into string expression
std::string sideSubst(const std::string& instr, Side s, bool clower = false);

/// opposite side to given side
Side otherSide(Side s);

/// "sign" for side, East = -1, West = +1
inline int ssign(Side s) { return s==EAST?-1:s==WEST?1:0; }

/// directions, to clean up wireplane direction confusion
enum AxisDirection {
	X_DIRECTION = 0,
	Y_DIRECTION = 1,
	Z_DIRECTION = 2,
	T_DIRECTION = 3
};
/// iteration to next axis
inline AxisDirection& operator++(AxisDirection& d) { return d = AxisDirection(d+1); }

/// types of UCNA data runs
enum RunType {
	UNKNOWN_RUN,	///< unknown/invalid run
	ASYMMETRY_RUN,	///< background or beta decay asymmetry run
	DEPOL_RUN,		///< depolarization run
	LED_RUN,		///< LED linearity run
	SOURCE_RUN,		///< source calibration run
	XE_RUN			///< activated Xenon calibration run
};

/// octet type specification
/// 0 -> not in octet
/// 1...12 -> A1...A12
/// 13...24 -> B1...B12
enum OctetRole {
	OCTR_UNKNOWN	= 0,	///< unknown run type
	OCTR_A1			= 1,	///< start of A octet
	OCTR_A12		= 12,	///< end of A octet
	OCTR_B1			= 13,	///< start of B octet
	OCTR_B12		= 24,	///< end of B octet
	OCTR_BG			= 25,	///< generic background run, undefined AFP
	OCTR_FG			= 26	///< generic foreground run, undefined AFP
};
inline OctetRole& operator++(OctetRole& d) { return d = OctetRole(d+1); }

/// state of AFP during run
enum AFPState {
	AFP_OFF		=	0,	///< AFP is off
	AFP_ON		=	1,	///< AFP is on
	AFP_OTHER	=	2,	///< AFP state unspecified
	AFP_OFF2ON	=	3,	///< (Depolarization run) AFP off->on
	AFP_ON2OFF	=	4	///< (Depolarization run) AFP on->off
};

inline AFPState& operator++(AFPState& a) { return a = AFPState(a+1); }
/// string for AFP state
std::string afpWords(AFPState afp);
/// convert string back to AFPState
AFPState strToAfp(const std::string& s);

/// state of gate valve during a run
enum GVState {
	GV_CLOSED	= 0,	///< GV closed (background)
	GV_OPEN		= 1,	///< GV open (beta decay)
	GV_OTHER	= 2		///< GV state unspecified
};

/// run grouping level
enum RunGrouping {
	GROUP_RUN		= 0,	///< single run
	GROUP_FGBG		= 1,	///< foreground/background pair
	GROUP_PPAIR		= 2,	///< AFP On/Off pulse pair
	GROUP_QUARTET	= 3,	///< A or B pair of pulse pairs
	GROUP_OCTET		= 4,	///< A->B or B->A octet
	GROUP_RANGE		= 5		///< arbitrary run range
};
/// words for octet grouping
std::string groupWords(RunGrouping g);
/// return next lower division depth
inline RunGrouping subdivide(RunGrouping g) { return RunGrouping(g-1); }


/// iteration to next GV state
inline GVState& operator++(GVState& a) { return a = GVState(a+1); }
/// string for GV state
std::string gvWords(GVState gv);
/// convert string to GV state
GVState strToGV(const std::string& s);

/// SCS geometry for a run
enum RunGeometry {
	GEOMETRY_OTHER,	///< unspecified geometry
	GEOMETRY_2007,	///< 2007 Geometry: 2.5u mylar, 0.3u Be, 25u wirechamber windows
	GEOMETRY_A,		///< 'A' 0.7u mylar, 0.3u Be, 25u wirechamber windows
	GEOMETRY_B,		///< 'B' 0.7+12.5u mylar, 0.3u Be, 25u wirechamber windows
	GEOMETRY_C,		///< 'C' 0.7u mylar, 0.3u Be, 6u wirechamber windows
	GEOMETRY_D		///< 'D' No endcaps, 6u wirechamber windows
};

enum DataQuality {
	DQ_UNKNOWN		=	0,		///< unknown data quality
	DQ_CHECKED		=	1<<0,	///< data quality has been checked
	DQ_TRASH		=	1<<1,	///< run is garbage, DO NOT USE
	DQ_CORRUPT_END	=	1<<2,	///< data at end of run is corrupted
	DQ_BEAM_ISSUES	=	1<<3	///< beam quality dubious during run
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
	ANCHOICE_A=1,	///< 0, I, II+III (unseparated) to primary trigger side
	ANCHOICE_B=2,	///< 0, I
	ANCHOICE_C=3,	///< 0, I, Separate II/III with 4keV MWPC cut
	ANCHOICE_D=4,	///< 0
	ANCHOICE_E=5,	///< 0, I, Separate II/III with likelihood
	ANCHOICE_F=6,	///< I
	ANCHOICE_G=7,	///< II+III unseparated on trigger side
	ANCHOICE_H=8,	///< II/III hard cut
	ANCHOICE_I=9,	///< II/III probability assigned
	ANCHOICE_J=10,	///< II hard cut
	ANCHOICE_K=11,	///< III hard cut
	ANCHOICE_Z=26	///< placeholder
};
/// iteration to next analysis choice
inline AnalysisChoice& operator++(AnalysisChoice& ac) { return ac = AnalysisChoice(ac+1); }
/// choice name letter
char choiceLetter(AnalysisChoice a);

/// enumeration for Jianglai particle ID
enum PID {
	PID_SINGLE = 0,		///< "single" gamma event
	PID_BETA = 1,		///< beta event
	PID_MUON = 2,		///< muon event
	PID_LED = 3,		///< tagged LED event
	PID_PULSER = 4,		///< tagged Bi pulser event
	PID_UNKNOWN = 666	///< unknown/non-event
};

/// names for particle ID
std::string pidWords(PID p);

/// event backscattering types
enum EventType {
	TYPE_0_EVENT	= 0,	///< "Correct" events arrive on only one side of detector
	TYPE_I_EVENT	= 1,	///< Scatters to hit both scintillators and wirechambers
	TYPE_II_EVENT	= 2,	///< Scatters off wirechamber, hits scintillator and wirechamber on opposite side
	TYPE_III_EVENT	= 3,	///< Triggers and scatters off scintillator; triggers only wirechamber on opposite side
	TYPE_IV_EVENT	= 4,	///< Undetectable backscatters (also used for "unclassifiable" events)
	TYPE_NONEVENT	= 666	///< Unclassified discard event
};

/// iteration to next event type
inline EventType& operator++(EventType& tp) { return tp = EventType(tp+1); }

/// names for event types
std::string typeWords(EventType tp);

/// enumeration for what signal to use for wirechamber energy reconstruction
enum ChargeProxyType {
	CHARGE_PROXY_NONE	= 0,	///< charge signal undefined
	CHARGE_PROXY_ANODE	= 1,	///< use anode ADC for charge signal
	CHARGE_PROXY_CCLOUD = 2		///< use cathode "charge cloud size" for charge signal
};
/// iteration to next charge type
inline ChargeProxyType& operator++(ChargeProxyType& tp) { return tp = ChargeProxyType(tp+1); }
/// string name for charge proxy
std::string chargeProxyName(ChargeProxyType tp);
/// charge proxy from string
ChargeProxyType strToChgPrx(const std::string& s);

#endif
