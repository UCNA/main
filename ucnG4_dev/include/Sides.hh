/// \file Sides.hh \brief detector side enumeration
#ifndef SIDES_HH
/// make sure this file is only included once
#define SIDES_HH 1

#include <string>

/// detector side specification
enum Side {
	EAST = 0,		//< East side of detector
	WEST = 1,		//< West side of detector
	BOTH = 2,		//< Both sides of detector
	NO_SIDE = 3,	//< Neither side of detector
	BADSIDE = 4
};

/// iteration to next side
inline Side& operator++(Side& s) { return s = Side(s+1); }

/// get letter for side names
char sideNames(Side s);
/// get full word for side names
const char* sideWords(Side s);
/// substitute side names into string expression
std::string sideSubst(std::string instr, Side s);

/// opposite side to given side
Side otherSide(Side s);

/// "sign" for side, East = +1, West = -1
int sign(Side s);

#endif
