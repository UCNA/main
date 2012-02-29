#include "UCNAException.hh"

UCNAException::UCNAException(std::string tp): std::exception(), Stringmap() {
	insert("type",tp);
}

const char* UCNAException::what() const throw() { 
	msg = toString();
	return msg.c_str(); 
}
