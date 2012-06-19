#include "ControlMenu.hh"
#include <utility>
#include <cassert>

inputRequester::inputRequester(std::string d, void (*f)(std::deque<std::string>&,std::stack<std::string>&), unsigned int n):
namedInteractor(d), myFunc(f) {
	for(unsigned int i=0; i<n; i++)
		addArg("Option "+itos(i+1),"",NULL);
}

void inputRequester::doIt(std::deque<std::string>& args, std::stack<std::string>& stack) {
		
	for(unsigned int i=0; i<inputFilters.size(); i++) {
		if(inputFilters[i]) {
			inputFilters[i]->doIt(args,stack);
		} else if(!args.size()) {
			char indata[64];
			printf("%s",argNames[i].c_str());
			if(defaultArgs[i].size())
				printf(" [%s]",defaultArgs[i].c_str());
			printf(" > ");
			assert(fgets(indata,64,stdin));
			std::string argin = strip(indata);
			if(!argin.size())
				stack.push(defaultArgs[i]);
			else
				stack.push(argin);
		} else {
			stack.push(popString(args));
		}
	}	
		
	if(myFunc) (*myFunc)(args,stack);
}

void inputRequester::addArg(std::string s, std::string dflt, namedInteractor* filter) {
	argNames.push_back(s);
	defaultArgs.push_back(dflt);
	inputFilters.push_back(filter);
}

void inputRequester::setArgOpts(unsigned int i, std::string s, std::string dflt, namedInteractor* filter) { 
	assert(i<argNames.size()); 
	argNames[i] = s;
	defaultArgs[i] = dflt;
	inputFilters[i] = filter;
}

std::string inputRequester::getArgname(unsigned int i) const { 
	assert(i<argNames.size()); 
	return argNames[i];
}

std::string inputRequester::getDescription() {
	std::string s = name;
	if(argNames.size() > 0) {
		s += " (";
		for(unsigned int i=0; i<argNames.size(); i++) {
			if(inputFilters[i] && !getArgname(i).size()) {
				s += inputFilters[i]->getDescription();
			} else {
				s += getArgname(i);
				if(defaultArgs[i].size())
					s += " = "+defaultArgs[i];
			}
			if(i<argNames.size()-1)
				s += +", ";
		}
		s += + ")";
	}
	return s;
}

bool nameselector_default_softmatch(const std::string& a, const std::string& b) { return lower(a) == lower(b.substr(0,a.size())); }
std::string NameSelector::barf_control = "\033_BARF";
std::string NameSelector::exit_control = "\033_EXIT";

NameSelector::NameSelector(std::string t, std::string promptval, bool persist): inputRequester(t), catchAll(false), isPersistent(persist) {
	inputRequester::addArg(promptval);
	setSoftmatch();
}

void NameSelector::addChoice(std::string d, std::string nm, Selector_Option_Flags o, std::string mname, streamInteractor* action) {
	if(!nm.size())
		nm = itos(choiceNames.size()+1);
	assert(nameMap.find(nm) == nameMap.end());
	if(!mname.size())
		mname = nm;
	nameMap.insert(std::make_pair(nm,choiceNames.size()));
	choiceNames.push_back(nm);
	choiceDescrips.push_back(d);
	choiceOut.push_back(mname);
	oflags.push_back(o);
	actions.push_back(action);
}

void NameSelector::addSynonym(std::string arg0, std::string syn) {
	std::map<std::string,unsigned int>::iterator it = nameMap.find(arg0);
	assert(it != nameMap.end());
	NameSelector::addChoice(choiceDescrips[it->second],syn,
			  Selector_Option_Flags(oflags[it->second] | SELECTOR_SYNONYM | SELECTOR_HIDDEN),
			  choiceOut[it->second],actions[it->second]);
}

void NameSelector::displayOptions() {
	printf("%s:\n---------------------------\n",name.c_str());
	for(unsigned int i=0; i<choiceNames.size(); i++)
		if(!(oflags[i] & SELECTOR_HIDDEN))
			printf("%s\t%s\n",choiceNames[i].c_str(),choiceDescrips[i].c_str());
}		

std::string NameSelector::getDescription() {
	std::string s = name;
	if(defaultArgs[0].size())
		s += " = "+defaultArgs[0];
	return s; 
}

void NameSelector::doIt(std::deque<std::string>& args, std::stack<std::string>& stack) {
	bool forceBreak = false;
	do {
		
		// display options if choice not already made
		if(!args.size()) {
			displayOptions();
			printf("---------------------------\n");
		}
		
		// input request loop
		while(1) {
			
			inputRequester::doIt(args,stack);
			
			std::string myArg = popString(stack);
			
			// break for control sequences
			if(startsWith(myArg, exit_control) || startsWith(myArg, barf_control) ) {
				if(startsWith(myArg, barf_control))
					args.push_front(myArg);
				forceBreak = true;
				break;
			}
			
			// try again on null input
			if(!myArg.size())
				continue;
			
			// find matching selection name
			std::map<std::string,unsigned int>::iterator it = nameMap.find(myArg);
			
			// if no direct match, apply soft matching when available
			if(it == nameMap.end() && softmatch) {
				
				// find soft matches
				std::vector<std::map<std::string,unsigned int>::iterator> matches;
				for(std::map<std::string,unsigned int>::iterator i = nameMap.begin(); i != nameMap.end(); i++)
					if(!(oflags[i->second] & SELECTOR_DISABLED) && softmatch(myArg,i->first))
						matches.push_back(i);
				
				// exactly one soft match is just right
				if(matches.size() == 1) {
					myArg = matches.back()->first;
					it = matches.back();
				}
				// too many ambiguous soft matches
				else if(matches.size() > 1) {
					printf("Error: ambiguous selection from:\n");
					for(unsigned int i=0; i<matches.size(); i++)
						printf("\t%s\n",matches[i]->first.c_str());
					continue;
				}
			}
			
			if(it == nameMap.end() || (oflags[it->second] & SELECTOR_DISABLED)) {
				if(catchAll) {
					stack.push(myArg);
					catchAll->doIt(args,stack);
					break;
				} else {
					printf("Error: unknown selection '%s'\n",myArg.c_str());
				}
			} else {
				if(actions[it->second])
					actions[it->second]->doIt(args,stack);
				else
					stack.push(choiceOut[it->second]);
				break;
			}
		}
		
	} while(isPersistent && !forceBreak);
}


void OptionsMenu::addChoice(namedInteractor* M, std::string nm, Selector_Option_Flags o) {
	NameSelector::addChoice(M->getDescription(),nm,o,"",M);
}

void menutils_PrintQue(std::deque<std::string>& args, std::stack<std::string>&) {
	printf("[ ");
	for(std::deque<std::string>::iterator it = args.begin(); it != args.end(); it++)
		printf("'%s' ",it->c_str());
	printf("]\n");
}
void menutils_PrintStack(std::deque<std::string>&, std::stack<std::string>& stack) {
	std::stack<std::string> tmp;
	while(stack.size())
		tmp.push(streamInteractor::popString(stack));
	printf("[ ");
	while(tmp.size()) {
		stack.push(streamInteractor::popString(tmp));
		printf("'%s' ",stack.top().c_str());
	}
	printf("]\n");
}
void menutils_StackSize(std::deque<std::string>&, std::stack<std::string>& stack) { stack.push(itos(stack.size())); }
bool menutils_CheckStackSize(unsigned int n, std::deque<std::string>& args, std::stack<std::string>& stack) {
	menutils_StackSize(args,stack);
	unsigned int nstack = streamInteractor::popInt(stack);
	if(nstack < n) {
		args.push_front(NameSelector::barf_control+" Insufficient Arguments ["+itos(n-nstack)+"]");
		return false;
	}
	return true;
}
void menutils_Drop(std::deque<std::string>& args, std::stack<std::string>& stack) {
	if(!menutils_CheckStackSize(1,args,stack))
		return;
	stack.pop();
}
void menutils_Dup(std::deque<std::string>& args, std::stack<std::string>& stack) {
	if(!menutils_CheckStackSize(1,args,stack))
		return;
	stack.push(stack.top());
}
void menutils_DropN(std::deque<std::string>& args, std::stack<std::string>& stack) {
	if(!menutils_CheckStackSize(1,args,stack))
		return;
	unsigned int n = streamInteractor::popInt(stack);
	if(!menutils_CheckStackSize(n,args,stack))
		return;
	while(n--)
		stack.pop();
}
void menutils_ClearStack(std::deque<std::string>&, std::stack<std::string>& stack) { while(!stack.empty()) stack.pop(); }
void menutils_Swap(std::deque<std::string>& args, std::stack<std::string>& stack) {
	if(!menutils_CheckStackSize(2,args,stack))
		return;
	std::string a = streamInteractor::popString(stack);
	std::string b = streamInteractor::popString(stack);
	stack.push(a);
	stack.push(b);
}
void menutils_Rot(std::deque<std::string>& args, std::stack<std::string>& stack) {
	if(!menutils_CheckStackSize(1,args,stack))
		return;
	unsigned int n = streamInteractor::popInt(stack);
	if(!n || !menutils_CheckStackSize(n,args,stack))
		return;
	for(unsigned int i=0; i<n-1; i++)
		args.push_front(streamInteractor::popString(stack));
	std::string s = streamInteractor::popString(stack);
	for(unsigned int i=0; i<n-1; i++)
		stack.push(streamInteractor::popString(args));
	stack.push(s);
}
void menutils_Select(std::deque<std::string>& args, std::stack<std::string>& stack) {
	if(!menutils_CheckStackSize(3,args,stack))
		return;
	std::string c = streamInteractor::popString(stack);
	std::string b = streamInteractor::popString(stack);
	std::string a = streamInteractor::popString(stack);
	if(c == "true" || atof(c.c_str()))
		stack.push(a);
	else
		stack.push(b);
}
void menutils_Exec(std::deque<std::string>& args, std::stack<std::string>& stack) {
	if(!menutils_CheckStackSize(1,args,stack))
		return;
	std::vector<std::string> v = split(streamInteractor::popString(stack));
	while(v.size()) {
		args.push_front(v.back());
		v.pop_back();
	}
}
void menutils_Barf(std::deque<std::string>& args, std::stack<std::string>&) { args.push_front(NameSelector::barf_control); }
void menutils_Exit(std::deque<std::string>& args, std::stack<std::string>&) { args.push_front(NameSelector::exit_control); }


