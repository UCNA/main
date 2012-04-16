#ifndef OUTPUTMANAGER_HH
#define OUTPUTMANAGER_HH 1

#include "QFile.hh"
#include <cassert>
#include <TObject.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>

/// collection of saved TObjects
class TObjCollector {
public:
	/// destructor
	virtual ~TObjCollector() { clearItems(); }
	virtual void writeItems() {
		printf("Saving registered objects...");
		fflush(stdout);
		for(std::vector<TObject*>::iterator i = rootItems.begin(); i != rootItems.end(); i++)
			(*i)->Write();
		printf(" Done.\n");
	}
	void clearItems() {
		for(std::vector<TObject*>::iterator i = rootItems.begin(); i != rootItems.end(); i++)
			delete(*i);
		rootItems.clear();
		for(std::vector<TObject*>::iterator i = deleteItems.begin(); i != deleteItems.end(); i++)
			delete(*i);
		deleteItems.clear();
	}
	/// register a root object for output (and eventual deletion)
	virtual TObject* addObject(TObject* o, bool noWrite=false) { if(noWrite) deleteItems.push_back(o); else rootItems.push_back(o); return o; }
	std::vector<TObject*> rootItems;	//< objects held until deleted
	std::vector<TObject*> deleteItems;	//< other objects never written to file
};

/// indicator for level of problem with analysis
enum WarningLevel {
	BENIGN_WARNING,		//< just FYI for troubleshooting
	MODERATE_WARNING,	//< might be a problem to look at
	SEVERE_WARNING,		//< almost certainly something is very wrong
	FATAL_WARNING		//< data is corrupted and cannot be analyzed
};

/// manages output directory for grouping related information; manages a canvas, output QFile, output ROOT file, recursive subdirectories
class OutputManager: public TObjCollector {
public:
	/// constructor for top-level
	OutputManager(std::string nm, std::string bp);
	/// constructor for nested
	OutputManager(std::string nm, OutputManager* pnt);
	
	/// destructor
	virtual ~OutputManager() {
		if(writeRootOnDestruct) writeROOT();
		clearItems();
		if(rootOut) rootOut->Close();
		if(defaultCanvas && !parent) delete(defaultCanvas); 
	}
	
	/// attach a sub-output to this one
	virtual void attach(OutputManager* S) {
		printf("OutputManager '%s' attached to '%s'\n",S->name.c_str(),name.c_str());
		subouts.push_back(S);
	}
		
	/// data entries to include in parent's output
	virtual Stringmap finalWords() { return forParent; }
	
	/// generate a TH1F registered with this runs output objects list
	TH1F* registeredTH1F(std::string hname, std::string htitle, unsigned int nbins, float x0, float x1);
	/// generate a TH2F registered with this runs output objects list
	TH2F* registeredTH2F(std::string hname, std::string htitle, unsigned int nbinsx, float x0, float x1, unsigned int nbinsy, float y0, float y1);
	/// print current canvas
	virtual void printCanvas(std::string fname, std::string suffix=".pdf") const;

	/// put a data quality warning in the parent output file
	void warn(WarningLevel l, std::string descrip, Stringmap M = Stringmap());
	
	/// write output QFile
	virtual void write(std::string outName = "");
	/// open output ROOT file for writing (useful if output tree is being created and needs a home)
	void openOutfile();
	/// set whether to write ROOT output when destructed
	void setWriteRoot(bool w) { writeRootOnDestruct = w; }

	QFile qOut;					//< QFile output
	TFile* rootOut;				//< ROOT file output
	Stringmap forParent;		//< items to include in parent's output
	TCanvas* defaultCanvas;		//< canvas for drawing plots
	OutputManager* parent;		//< parent output manager
	std::string basePath;		//< general output path
	std::string plotPath;		//< specific output path for plots
	std::string dataPath;		//< specific output path for output data
	std::string rootPath;		//< specific output path for ROOT files
	std::string name;			//< name for this subsystem			
	bool writeRootOnDestruct;	//< whether to write ROOT file when destructed
	
protected:
	
	/// set name of this detector
	virtual void setName(std::string nm);
	/// write output ROOT file; WARNING: THIS DELETES ALL REGISTERED ITEMS; do last if you reference these.
	void writeROOT();
	
	std::vector<OutputManager*> subouts;	//< output subdirectories
};

#endif
