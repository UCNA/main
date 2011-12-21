#include "Sides.hh"
#include <stdio.h>

char sideNames(Side s) {
	const char snm[] = {'E','W','B','N'};
	return snm[s];
}

const char* sideWords(Side s) {
	const char* swd[] = {"East","West","Both","None"};
	return swd[s];
}

// TODO: make this robust
std::string sideSubst(std::string instr, Side s) {
	char tmp[1024];
	if(instr.find("%c") != instr.npos)
		sprintf(tmp,instr.c_str(),sideNames(s));
	if(instr.find("%s") != instr.npos)
		sprintf(tmp,instr.c_str(),sideWords(s));
	return std::string(tmp);
}

Side otherSide(Side s) {
	switch(s) {
		case EAST:
			return WEST;
		case WEST:
			return EAST;
		case BOTH:
			return NONE;
		default:
			return BOTH;
	}
}

int sign(Side s) {
	if(s==EAST)
		return 1;
	if(s==WEST)
		return -1;
	return 0;
}
