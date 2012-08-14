/// \file SMExcept.hh \brief exception handling class
#ifndef UCNAEXCEPTION_HH
/// make sure this file is only included once
#define UCNAEXCEPTION_HH 1

#include <exception>
#include "QFile.hh"

/// exception class for error handling
class SMExcept: public std::exception, public Stringmap {
public:
	/// constructor
	SMExcept(std::string tp);
	/// destructor
	~SMExcept() throw() {}
	/// display error
	virtual const char* what() const throw();
	/// string for holding error message
	mutable std::string msg;
};

#endif
