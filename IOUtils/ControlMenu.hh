#ifndef CONTROLMENU_HH
#define CONTROLMENU_HH 1

#include "strutils.hh"
#include <cassert>
#include <algorithm>
#include <deque>
#include <stack>
#include <stdio.h>
#include <stdlib.h>
#include <map>

/// base class for interacting with a deque of strings
class streamInteractor {
public:
	/// constructor
	streamInteractor() {}
	/// destructor
	virtual ~streamInteractor() {}
	/// do something! (put something useful here in subclasses)
	virtual void doIt(std::deque<std::string>&, std::stack<std::string>&) {}
	
	/// pop string off front of deque
	static std::string popString(std::deque<std::string>& args) { std::string s = args.front(); args.pop_front(); return s; }
	/// pop int off front of deque
	static int popInt(std::deque<std::string>& args) { return atoi(popString(args).c_str()); }
	/// pop float off front of deque
	static float popFloat(std::deque<std::string>& args) { return atof(popString(args).c_str()); }
	
	/// pop int off stack
	static std::string popString(std::stack<std::string>& stack) { std::string s = stack.top(); stack.pop(); return s; }
	/// pop int off stack
	static int popInt(std::stack<std::string>& stack) { return atoi(popString(stack).c_str()); }
	/// pop float off stack
	static float popFloat(std::stack<std::string>& stack) { return atof(popString(stack).c_str()); }
};


/// stream interactor with a name for screen display
class namedInteractor: public streamInteractor {
public:
	/// constructor
	namedInteractor(std::string nm): name(nm) {}
	/// get my name/description
	virtual std::string getDescription() { return name; }
protected:
	std::string name;
};


/// stream interactor with named/numbered arguments that can prompt for necessary input
class inputRequester: public namedInteractor {
public:
	/// constructor
	inputRequester(std::string d, void (*f)(std::deque<std::string>&,std::stack<std::string>&) = NULL, unsigned int n = 0);
	/// action when selected
	virtual void doIt(std::deque<std::string>& args, std::stack<std::string>& stack);
	/// add new argument
	virtual void addArg(std::string s, std::string dflt = "", namedInteractor* filter = NULL);
	/// set argument parameters
	virtual void setArgOpts(unsigned int i, std::string s, std::string dflt = "", namedInteractor* filter = NULL);
	/// get option name
	virtual std::string getArgname(unsigned int i) const;
	///get my name/description
	virtual std::string getDescription();
	
protected:
	std::vector<std::string> argNames;				//< names/descriptions of arguments
	std::vector<std::string> defaultArgs;			//< default values for arguments
	std::vector<namedInteractor*> inputFilters;	//< input data filters
	void (*myFunc)(std::deque<std::string>&,std::stack<std::string>&);		//< function to do when selected
};


/// option display/activity flags
enum Selector_Option_Flags {
	SELECTOR_NORMAL = 0,		//< normal mode
	SELECTOR_HIDDEN = 1<<0,		//< option is hidden in menu
	SELECTOR_DISABLED = 1<<1,	//< option is inactive in menu
	SELECTOR_SYNONYM = 1<<2		//< option is a synonym for another option
};

/// default soft-matching routine
bool nameselector_default_softmatch(const std::string& a, const std::string& b);

/// class for selecting an item from a list and applying an associated action
class NameSelector: public inputRequester {
public:
	/// constructor
	NameSelector(std::string t, std::string promptval = "Selection", bool persist = false);
	/// add selection choice
	virtual void addChoice(std::string d, std::string nm = "", Selector_Option_Flags o = SELECTOR_NORMAL, std::string mname = "", streamInteractor* action = NULL);
	/// display available options
	virtual void displayOptions();
	/// get choice from input queue, return selected item on stack
	virtual void doIt(std::deque<std::string>& args, std::stack<std::string>& stack);
	/// get my name/description
	virtual std::string getDescription();
	/// set default choice
	virtual void setDefault(std::string s) { defaultArgs[0] = s; }
	/// set catchall action
	virtual void setCatchall(streamInteractor* I) { catchAll = I; }
	/// prevent adding arguments (doesn't make sense in this context)
	virtual void addArg(std::string, std::string, namedInteractor*) { assert(false); }
	/// add a synonym for an existing argument
	virtual void addSynonym(std::string arg0, std::string syn);
	/// set soft-matching function (set to NULL to disable soft matching)
	void setSoftmatch(bool (*f)(const std::string& a, const std::string& b) = &nameselector_default_softmatch) { softmatch = f; }
	static std::string barf_control;
	static std::string exit_control;
	
protected:
	std::map<std::string,unsigned int> nameMap;		//< map from choice names to selected content
	std::vector<std::string> choiceNames;			//< choice names
	std::vector<std::string> choiceDescrips;		//< choice descriptions
	std::vector<std::string> choiceOut;				//< output for each choice
	std::vector<streamInteractor*> actions;			//< output action for each choice
	std::vector<Selector_Option_Flags> oflags;		//< option display flags
	streamInteractor* catchAll;						//< catch-all action for unidentified choices
	bool isPersistent;								//< whether menu is persistent (repeats after selection is made)
	bool (*softmatch)(const std::string& a, const std::string& b);	//< function for "soft matching" comparison of selections
};

/// Text menu of selectable items (simplified NameSelector specifically for actions)
class OptionsMenu: public NameSelector {
public:
	/// constructor
	OptionsMenu(std::string t, bool persist = true): NameSelector(t,"Selection",persist) { }
	/// destructor
	virtual ~OptionsMenu() {}
	/// add choice to selections list
	virtual void addChoice(namedInteractor* M, std::string nm = "", Selector_Option_Flags o = SELECTOR_NORMAL);
	/// prevent adding choices through base mechanism
	virtual void addChoice(std::string, std::string, Selector_Option_Flags, std::string, streamInteractor*) { assert(false); }
};



/// print contents of control queue
void menutils_PrintQue(std::deque<std::string>& args, std::stack<std::string>&);
/// print contents of control queue
void menutils_PrintStack(std::deque<std::string>&, std::stack<std::string>& stack);
/// push stack size onto stack
void menutils_StackSize(std::deque<std::string>&, std::stack<std::string>& stack);
/// check that there are at least the specified number of items on the stack
bool menutils_CheckStackSize(unsigned int n, std::deque<std::string>& args, std::stack<std::string>& stack);
/// drop top stack item
void menutils_Drop(std::deque<std::string>& args, std::stack<std::string>& stack);
/// duplicate stack item
void menutils_Dup(std::deque<std::string>& args, std::stack<std::string>& stack);
/// drop top n stack items
void menutils_DropN(std::deque<std::string>& args, std::stack<std::string>& stack);
/// clear stack
void menutils_ClearStack(std::deque<std::string>&, std::stack<std::string>& stack);
/// swap top two stack items
void menutils_Swap(std::deque<std::string>& args, std::stack<std::string>& stack);
/// rotate n stack items, bringing nth item to top
void menutils_Rot(std::deque<std::string>& args, std::stack<std::string>& stack);
/// select c?a:b
void menutils_Select(std::deque<std::string>& args, std::stack<std::string>& stack);
/// move string on stack to command stream for execution
void menutils_Exec(std::deque<std::string>& args, std::stack<std::string>& stack);
/// add error code ("BARF") to control queue
void menutils_Barf(std::deque<std::string>& args, std::stack<std::string>&);
/// add exit code ("EXIT") to control queue
void menutils_Exit(std::deque<std::string>& args, std::stack<std::string>&);

#endif
