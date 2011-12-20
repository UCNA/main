#include "EvisConverter.hh"

EvisConverter::EvisConverter(RunNum rn, CalDB* CDB) {
	bool hasConverters = true;
	for(Side s = EAST; s <= WEST; ++s) {
		for(unsigned int tp = TYPE_0_EVENT; tp <= TYPE_II_EVENT; tp++) {
			conversions[s][tp] = CDB->getEvisConversion(rn,s,EventType(tp));
			hasConverters = hasConverters && conversions[s][tp];
		}
	}
	if(!hasConverters)
		printf("******* WARNING: missing Evis->Etrue conversion for run %i ********\n",rn);
}

float EvisConverter::Etrue(Side s, EventType tp, float EvisE, float EvisW) const {
	assert((s==EAST || s==WEST) && tp >= TYPE_0_EVENT);
	if(tp>=TYPE_III_EVENT)
		tp=TYPE_II_EVENT;
	float Evis = (tp==TYPE_I_EVENT)? EvisE+EvisW:(s==EAST?EvisE:EvisW);
	if(!conversions[s][tp]) return Evis;
	return conversions[s][tp]->Eval(Evis);
}
