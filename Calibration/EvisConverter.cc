#include "EvisConverter.hh"

EvisConverter::EvisConverter(RunNum rn, CalDB* CDB) {
	bool hasConverters = true;
	for(Side s = EAST; s <= WEST; ++s) {
		for(unsigned int tp = TYPE_0_EVENT; tp <= TYPE_III_EVENT; tp++) {
			conversions[s][tp] = CDB->getEvisConversion(rn,s,EventType(tp));
			hasConverters = hasConverters && conversions[s][tp];
		}
	}
	if(!hasConverters)
		printf("******* WARNING: missing Evis->Erecon conversion for run %i ********\n",rn);
}

EvisConverter::~EvisConverter() {
	for(Side s = EAST; s <= WEST; ++s)
		for(unsigned int tp = TYPE_0_EVENT; tp <= TYPE_III_EVENT; tp++)
			if(conversions[s][tp])
				delete conversions[s][tp];
}

float EvisConverter::Erecon(Side s, EventType tp, float EvisE, float EvisW) const {
	smassert((s==EAST || s==WEST));
	float Evis = (tp==TYPE_I_EVENT)? EvisE+EvisW:(s==EAST?EvisE:EvisW);
	if(tp>TYPE_III_EVENT || !conversions[s][tp]) return Evis;
	return conversions[s][tp]->Eval(Evis);
}
