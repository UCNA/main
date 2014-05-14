#include "Enums.hh"
#include "strutils.hh"
#include <math.h>
#include <stdio.h>
#include "SMExcept.hh"

//---------------------------------------

char sideNames(Side s, bool clower) {
	smassert(s<BADSIDE);
	const char snm[] = {'E','W','B','N'};
	const char snml[] = {'e','w','b','n'};
	return clower?snml[s]:snm[s];
}

const char* sideWords(Side s) {
	smassert(s<=BADSIDE);
	const char* swd[] = {"East","West","Both","None","BadSide"};
	return swd[s];
}

std::string pidWords(PID p) {
	const char* twd[] = {"gamma","beta","muon","LED","BiPulser"};
	if(p<=PID_PULSER) return twd[p];
	return "unknown";
}

std::string typeWords(EventType tp) {
	smassert(tp<=TYPE_IV_EVENT);
	const char* twd[] = {"Type0","TypeI","TypeII","TypeIII","TypeIV"};
	return twd[tp];
}

const char* dbSideName(Side s) {
	smassert(s<=BADSIDE);
	const char* swd[] = {"'East'","'West'","'Both'","'None'","'BadSide'"};
	return swd[s];
}

std::string sideSubst(const std::string& instr, Side s, bool clower) {
	std::vector<std::string> segs = split(instr,"%");
	if(!segs.size()) return "";
	std::string sout = (instr[0]=='%')?"":segs[0];
	for(unsigned int i=(instr[0]!='%'); i<segs.size(); i++) {
		if(!segs[i].size()) continue;
		if(segs[i][0]=='c') sout += ctos(sideNames(s,clower));
		else if(segs[i][0]=='s') sout += sideWords(s);
		else sout += ctos(segs[i][0]);
		sout += segs[i].substr(1,segs[i].size());
	}
	return sout;
}

Side otherSide(Side s) {
	return s==EAST?WEST : s==WEST?EAST : s==BOTH?NOSIDE : BOTH;
}

unsigned int pmtHardwareNum(Side s, unsigned int quadrant) {
	return s==EAST?(quadrant+2)%4+1 : s==WEST?quadrant+1 : 0;
}

Side strToSide(const std::string& s) {
	return s=="East"?EAST:s=="West"?WEST:s=="Both"?BOTH:s=="None"?NOSIDE:BADSIDE;
}

//---------------------------------------

std::string afpWords(AFPState afp) {
	if(afp>AFP_ON2OFF) afp = AFP_OTHER;
	const char* awd[] = {"Off","On","Other","Off2On","On2Off"};
	return awd[afp];
}

AFPState strToAfp(const std::string& s) {
	return s=="On"?AFP_ON:s=="Off"?AFP_OFF:s=="On2Off"?AFP_ON2OFF:s=="Off2On"?AFP_OFF2ON:AFP_OTHER;
}

std::string gvWords(GVState gv) {
	return gv==GV_OPEN?"Open":gv==GV_CLOSED?"Closed":"Other";
}

GVState strToGV(const std::string& s) {
	return s=="Open"?GV_OPEN:s=="Closed"?GV_CLOSED:GV_OTHER;
}

//---------------------------------------

std::string groupWords(RunGrouping g) {
	return (g==GROUP_RUN?"run" :
			g==GROUP_FGBG?"fgbg" :
			g==GROUP_PPAIR?"ppair" :
			g==GROUP_QUARTET?"quartet" :
			g==GROUP_OCTET?"octet" :
			"range");
}

//---------------------------------------

RunGeometry whichGeometry(RunNum rn) {
	if( rn < 1000 )
		return GEOMETRY_OTHER;
	if( rn < 7000 )
		return GEOMETRY_2007;
	if( rn < 9220)
		return GEOMETRY_A;
	if( rn < 10400)
		return GEOMETRY_B;
	if( rn < 11390)
		return GEOMETRY_C;
	if( rn < 11501)
		return GEOMETRY_D;
	return GEOMETRY_C;
}

std::string geomName(RunGeometry g) {
	if( g == GEOMETRY_2007 )
		return "2007";
	if( g == GEOMETRY_A )
		return "A";
	if( g == GEOMETRY_B )
		return "B";
	if( g == GEOMETRY_C )
		return "C";
	if( g == GEOMETRY_D )
		return "D";
	return "OTHER";
}

char choiceLetter(AnalysisChoice a) {
	return (ANCHOICE_A <= a && a <= ANCHOICE_Z)?'A'+(a-ANCHOICE_A):'0';
}

//----------------------------------------

std::string chargeProxyName(ChargeProxyType tp) {
	if(tp==CHARGE_PROXY_ANODE) return "anode";
	else if(tp==CHARGE_PROXY_CCLOUD) return "ccloud";
	return "None";
}

ChargeProxyType strToChgPrx(const std::string& s) {
	if(s=="anode") return CHARGE_PROXY_ANODE;
	else if(s=="ccloud") return CHARGE_PROXY_CCLOUD;
	return CHARGE_PROXY_NONE;
}


