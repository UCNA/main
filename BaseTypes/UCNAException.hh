/// \file UCNAException.hh \brief exception handling class
#ifndef UCNAEXCEPTION_HH
/// make sure this file is only included once
#define UCNAEXCEPTION_HH 1

#include <exception>
#include "QFile.hh"

/// exception class for error handling
class UCNAException: public std::exception, public Stringmap {
public:
	/// constructor
	UCNAException(std::string tp);
	/// destructor
	~UCNAException() throw() {}
	/// display error
	virtual const char* what() const throw();
	/// string for holding error message
	mutable std::string msg;
};

#endif
