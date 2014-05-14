#include "SMExcept.hh"

SMExcept::SMExcept(const std::string& tp): std::exception(), Stringmap() {
	insert("type",tp);
}

const char* SMExcept::what() const throw() { 
	msg = toString();
	return msg.c_str(); 
}

void smassert(bool b, const std::string& tp, const Stringmap& m) {
	if(!b) {
		SMExcept e(tp);
		e += m;
		throw e;
	}
}

