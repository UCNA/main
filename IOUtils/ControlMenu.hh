#ifndef CONTROLMENU_HH
#define CONTROLMENU_HH

#include "strutils.hh"
#include <cassert>
#include <algorithm>
#include <deque>
#include <stack>
#include <stdio.h>
#include <stdlib.h>
#include <map>

/// base class for interacting with a deque of strings
class StreamInteractor {
public:
	/// constructor
	StreamInteractor(): mydeque(NULL), mystack(NULL)  {}
	/// destructor
	virtual ~StreamInteractor() {}
	/// do something! (put something useful here in subclasses)
	virtual void doIt() {}
	
	/// check that there are at least the specified number of items on the stack
	bool menutils_CheckStackSize(unsigned int n);

	/// pop string off front of deque
	std::string popStringD() { assert(mydeque); std::string s = mydeque->front(); mydeque->pop_front(); return s; }
	/// pop int off front of deque
	int popIntD() { return atoi(popStringD().c_str()); }
	/// pop float off front of deque
	float popFloatD() { return atof(popStringD().c_str()); }
	
	/// pop int off stack
	std::string popString() { assert(mystack); std::string s = mystack->top(); mystack->pop(); return s; }
	/// pop int off stack
	int popInt() { return atoi(popString().c_str()); }
	/// pop float off stack
	float popFloat() { return atof(popString().c_str()); }
	
	std::deque<std::string>* mydeque;	///< command arguments deque
	std::stack<std::string>* mystack;	///< working space stack
};

/// stream interactor with a name for screen display
class NamedInteractor: public StreamInteractor {
public:
	/// constructor
	NamedInteractor(std::string nm): name(nm) {}
	/// get my name/description
	virtual std::string getDescription() { return name; }

	std::string name;	///< name for this interactor
};


/// stream interactor with named/numbered arguments that can prompt for necessary input
class InputRequester: public NamedInteractor {
public:
	/// constructor
	InputRequester(std::string d, void (*f)(StreamInteractor*) = NULL, StreamInteractor* fObj = NULL);
	/// action when selected
	virtual void doIt();
	/// add new argument
	virtual void addArg(const std::string& s, const std::string& dflt = "", const std::string& descrip = "", NamedInteractor* filter = NULL);
	/// add new argument, assuming descriptions come from filter
	virtual void addArg(NamedInteractor* filter, const std::string& s = "") { addArg(s,"","",filter); }
	/// set argument parameters
	virtual void setArgOpts(unsigned int i, std::string s, std::string dflt = "", NamedInteractor* filter = NULL);
	/// get option name
	virtual std::string getArgname(unsigned int i) const;
	///get my name/description
	virtual std::string getDescription();
	
	static InputRequester exitMenu;
	
protected:
	std::vector<std::string> argNames;				///< names of arguments
	std::vector<std::string> argDescrips;			///< extended descriptions of arguments
	std::vector<std::string> defaultArgs;			///< default values for arguments
	std::vector<NamedInteractor*> inputFilters;		///< input data filters
	StreamInteractor* myFuncObject;					///< object called with function
	void (*myFunc)(StreamInteractor*);				///< function to do when selected
};


/// option display/activity flags
enum Selector_Option_Flags {
	SELECTOR_NORMAL = 0,		///< normal mode
	SELECTOR_HIDDEN = 1<<0,		///< option is hidden in menu
	SELECTOR_DISABLED = 1<<1,	///< option is inactive in menu
	SELECTOR_SYNONYM = 1<<2		///< option is a synonym for another option
};

/// default soft-matching routine
bool nameselector_default_softmatch(const std::string& a, const std::string& b);

/// class for selecting an item from a list and applying an associated action
class NameSelector: public InputRequester {
public:
	/// constructor
	NameSelector(std::string t, std::string promptval = "Selection", bool persist = false);
	/// add selection choice
	virtual void addChoice(std::string d, std::string nm = "", Selector_Option_Flags o = SELECTOR_NORMAL, std::string mname = "", StreamInteractor* action = NULL);
	/// display available options
	virtual void displayOptions();
	/// get choice from input queue, return selected item on stack
	virtual void doIt();
	/// get my name/description
	virtual std::string getDescription();
	/// set default choice
	virtual void setDefault(std::string s) { defaultArgs[0] = s; }
	/// set catchall action
	virtual void setCatchall(StreamInteractor* SI) { catchAll = SI; }
	/// prevent adding arguments (doesn't make sense in this context)
	virtual void addArg(const std::string&, const std::string& = "", const std::string& = "", NamedInteractor* = NULL) { assert(false); }
	/// prevent adding arguments (doesn't make sense in this context)
	virtual void addArg(NamedInteractor*, const std::string& = "") { assert(false); }
	/// add a synonym for an existing argument
	virtual void addSynonym(std::string arg0, std::string syn);
	/// set soft-matching function (set to NULL to disable soft matching)
	void setSoftmatch(bool (*f)(const std::string& a, const std::string& b) = &nameselector_default_softmatch) { softmatch = f; }
	static std::string barf_control;
	static std::string exit_control;
	
protected:
	std::map<std::string,unsigned int> nameMap;		///< map from choice names to selected content
	std::vector<std::string> choiceNames;			///< choice names
	std::vector<std::string> choiceDescrips;		///< choice descriptions
	std::vector<std::string> choiceOut;				///< output for each choice
	std::vector<StreamInteractor*> actions;			///< output action for each choice
	std::vector<Selector_Option_Flags> oflags;		///< option display flags
	StreamInteractor* catchAll;						///< catch-all action for unidentified choices
	bool isPersistent;								///< whether menu is persistent (repeats after selection is made)
	bool (*softmatch)(const std::string& a, const std::string& b);	///< function for "soft matching" comparison of selections
};

/// Text menu of selectable items (simplified NameSelector specifically for actions)
class OptionsMenu: public NameSelector {
public:
	/// constructor
	OptionsMenu(std::string t, bool persist = true): NameSelector(t,"Selection",persist) { }
	/// destructor
	virtual ~OptionsMenu() {}
	/// add choice to selections list
	virtual void addChoice(NamedInteractor* M, std::string nm = "", Selector_Option_Flags o = SELECTOR_NORMAL);
	/// prevent adding choices through base mechanism
	virtual void addChoice(std::string, std::string, Selector_Option_Flags, std::string, StreamInteractor*) { assert(false); }
};



// standard utility functions
	
/// print contents of control queue
void menutils_PrintQue(StreamInteractor*);
/// print contents of control queue
void menutils_PrintStack(StreamInteractor*);
/// push stack size onto stack
void menutils_StackSize(StreamInteractor*);
/// drop top stack item
void menutils_Drop(StreamInteractor*);
/// duplicate stack item
void menutils_Dup(StreamInteractor*);
/// drop top n stack items
void menutils_DropN(StreamInteractor*);
/// clear stack
void menutils_ClearStack(StreamInteractor*);
/// swap top two stack items
void menutils_Swap(StreamInteractor*);
/// rotate n stack items, bringing nth item to top
void menutils_Rot(StreamInteractor*);
/// select c?a:b
void menutils_Select(StreamInteractor*);
/// move string on stack to command stream for execution
void menutils_Exec(StreamInteractor*);
/// add error code ("BARF") to control queue
void menutils_Barf(StreamInteractor*);
/// add exit code ("EXIT") to control queue
void menutils_Exit(StreamInteractor*);

#endif
