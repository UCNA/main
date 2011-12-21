#include "UCNAException.hh"

UCNAException::UCNAException(std::string tp): std::exception(), Stringmap() {
	insert("type",tp);
}

const char* UCNAException::what() const throw() { 
	return toString().c_str(); 
}
