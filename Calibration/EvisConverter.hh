#ifndef EVISCONVERTER_HH
#define EVISCONVERTER_HH 1

#include "Enums.hh"
#include "CalDB.hh"

/// Evis to Etrue conversion class
class EvisConverter {
public:
	/// constructor
	EvisConverter(RunNum rn, CalDB* CDB);	
	/// destructor
	virtual ~EvisConverter() { /*TODO*/ }
	/// get true energy for side given Evis on each side
	float Etrue(Side s, EventType tp, float EvisE, float EvisW) const;			
protected:
	TGraph* conversions[2][3];	//< energy conversion curves
};

#endif
