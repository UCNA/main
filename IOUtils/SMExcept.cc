#include "SMExcept.hh"

SMExcept::SMExcept(std::string tp): std::exception(), Stringmap() {
	insert("type",tp);
}

const char* SMExcept::what() const throw() { 
	msg = toString();
	return msg.c_str(); 
}
