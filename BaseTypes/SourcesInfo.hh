#ifndef SOURCESINFO_HH
#define SOURCESINFO_HH 1

#include "Source.hh"
#include <vector>
#include <string>

/// class for storing source information for a run
class SourcesInfo {
public:	
	/// constructor
	SourcesInfo(RunNum r = 0): sourceCounter(0) {
		
		if(!r) return;
		
		// attempt to read source holder contents
		
		RunNum rn;
		std::vector<std::string> sources[2];
		std::string l;
		
		printf("Loading run log...\n");
		std::ifstream fin("Aux/UCNA Run Log.txt");
		
		while (fin.good()) {
			
			std::getline(fin,l);
			if(!l.size()) continue;
			std::vector<std::string> words = split(l);
			if(!words.size()) continue;
			
			if(l[0]=='*') {
				if(words.size() < 2)
					continue;
				if(!sscanf(words[0].c_str(),"*%i",&rn))
					continue;
				if(rn == r) {
					sourceHolder[EAST] = sources[EAST];
					sourceHolder[WEST] = sources[WEST];
					break;
				}
			} else if(words[0]=="@sources") {
				
				bool notSide[2];
				notSide[EAST] = notSide[WEST] = false;
				if(words.size() >= 2) {
					if(words[1]=="E")
						notSide[WEST] = true;
					else if(words[1]=="W")
						notSide[EAST] = true;
				}
				for(Side s = EAST; s<=WEST; ++s) {
					if(notSide[s])
						continue;
					sources[s].clear();
					for(std::vector<std::string>::const_iterator it = words.begin(); it != words.end(); it++) {
						if(*it == "Sn")
							sources[s].push_back("Sn113");
						else if(*it == "Bi")
							sources[s].push_back("Bi207");
						else if(*it == "Cd")
							sources[s].push_back("Cd109");
						else if(*it == "Sr85")
							sources[s].push_back("Sr85");
						else if(*it == "Sr90")
							sources[s].push_back("Sr90");
						else if(*it == "In")
							sources[s].push_back("In114");
						else if(*it == "Ce")
							sources[s].push_back("Ce139");
					}
				}
			}
		}
		fin.close();
	}
	/// destructor
	~SourcesInfo() {}
	
	/// check whether this is a source run
	bool isSourceRun() const { return sources[EAST].size() + sources[WEST].size() > 0; }
	
	std::vector<std::string> sourceHolder[2];	//< sources in the source holder, left to right visible on each side
	
	/// add identified source to this runs' source list
	Source addSource(Side s, Source src) { src.sID = ++sourceCounter; src.mySide = s; sources[s].push_back(src); return src; }
	/// get start of sources list
	std::vector<Source>::const_iterator sourcesBegin(Side s) const { return sources[s].begin(); }
	/// get end of sources list
	std::vector<Source>::const_iterator sourcesEnd(Side s) const { return sources[s].end(); }
	
protected:
	std::vector<Source> sources[2];			//< position of radioactive sources on each side
	unsigned int sourceCounter;				//< counter for sequential source numbering
};

#endif
