#ifndef EVISCONVERTER_HH
#define EVISCONVERTER_HH

#include "Enums.hh"
#include "CalDB.hh"

/// Evis to Erecon conversion class
class EvisConverter {
public:
	/// constructor
	EvisConverter(RunNum rn, CalDB* CDB);	
	/// destructor
	virtual ~EvisConverter();
	/// get true energy for side given Evis on each side
	float Erecon(Side s, EventType tp, float EvisE, float EvisW) const;			
protected:
	TGraph* conversions[2][TYPE_III_EVENT+1];	///< energy conversion curves by [side][type]
};

#endif
