#include "Source.hh"

Source::Source(Stringmap S) {
	sID = (int)S.getDefault("sID", 0);
	t = S.getDefault("type","");
	x = S.getDefault("x",0);
	y = S.getDefault("y",0);
	wx = S.getDefault("wx",0);
	wy = S.getDefault("wy",0);
	nCounts = S.getDefault("nCounts",0);
	std::string sd = S.getDefault("side", "N");
	if(sd[0] == sideNames(EAST))
		mySide = EAST;
	else if(sd[0] == sideNames(WEST))
		mySide = WEST;
	else
		mySide = NOSIDE;
}

Stringmap Source::getProperties() const {
	Stringmap M;
	M.insert("object","Source");
	M.insert("x",x);
	M.insert("y",y);
	M.insert("wx",wx);
	M.insert("wy",wy);
	M.insert("nCounts",nCounts);
	M.insert("type",t);
	M.insert("sID",sID);
	M.insert("side",ctos(sideNames(mySide)));
	M.insert("name",name());
	
	return M;
}

std::vector<SpectrumPeak> Source::getPeaks() const {
	std::vector<SpectrumPeak> v;
	if(t=="Sn113")
		v.push_back(SpectrumPeak(SN_PEAK,sID,mySide));
	if(t=="Cd109")
		v.push_back(SpectrumPeak(CD109_PEAK,sID,mySide));
	if(t=="Cs137")
		v.push_back(SpectrumPeak(CS137_PEAK,sID,mySide));
	if(t=="In114" || t=="In114E" || t=="In114W")
		v.push_back(SpectrumPeak(IN114_PEAK,sID,mySide));
	if(t=="Ce139")
		v.push_back(SpectrumPeak(CE139_PEAK,sID,mySide));
	if(t=="Bi207") {
		v.push_back(SpectrumPeak(BI_PEAK_1,sID,mySide));
		v.push_back(SpectrumPeak(BI_PEAK_2,sID,mySide));
	}
	return v;
}
