#include "OutputManager.hh"
#include "PathUtils.hh"
#include <TH1.h>

OutputManager::OutputManager(std::string nm, std::string bp): rootOut(NULL), defaultCanvas(new TCanvas()),
parent(NULL), writeRootOnDestruct(false) {
	TH1::AddDirectory(kFALSE);
	// set up output canvas	
	defaultCanvas->SetFillColor(0);
	defaultCanvas->SetCanvasSize(300,300);
	
	basePath = plotPath = dataPath = rootPath = bp;
	setName(nm);
}

OutputManager::OutputManager(std::string nm, OutputManager* pnt):
rootOut(NULL), defaultCanvas(NULL), parent(pnt), writeRootOnDestruct(false) {
	TH1::AddDirectory(kFALSE);
	if(parent)
		defaultCanvas = parent->defaultCanvas;
	setName(nm);
}

void OutputManager::setName(std::string nm) {
	name = nm;
	if(parent) {
		plotPath = dataPath = basePath = rootPath = parent->basePath+"/"+name+"/";
	}
}

void OutputManager::warn(WarningLevel l, std::string descrip, Stringmap M) {
	
	M.insert("description",descrip);
	M.insert("subsystem",name);
	
	if(l==BENIGN_WARNING) {
		printf("* Warning: %s\n",descrip.c_str());
		M.insert("level","benign");
	} else if(l==MODERATE_WARNING) {
		printf("\n*** WARNING: %s\n",descrip.c_str());
		M.insert("level","moderate");
	} else if(l==SEVERE_WARNING) {
		printf("\n\n*******************\n* SEVERE WARNING: %s\n*******************\n",descrip.c_str());
		M.insert("level","severe");
	} if(l==FATAL_WARNING) {
		printf("\n\n*******************\n*******************\n** FATAL WARNING: %s\n*******************\n*******************\n",descrip.c_str());
		M.insert("level","fatal");
	}
	
	M.display("\t");
	printf("\n");
	
	if(parent)
		parent->qOut.insert("Warning",M);
	else
		qOut.insert("Warning",M);
}

void OutputManager::write(std::string outName) {
	
	// record subout info
	for(std::vector<OutputManager*>::iterator it = subouts.begin(); it != subouts.end(); it++) {
		Stringmap m = (*it)->finalWords();
		if(m.size())
			qOut.insert((*it)->name,m);
	}
	
	// write run data file
	if(qOut.size()) {
		makePath(dataPath+"/"+outName,true);
		if(outName.size())
			qOut.setOutfile(dataPath+"/"+outName);
		else
			qOut.setOutfile(dataPath+"/"+name+".txt");
		qOut.commit();
	}
}

void OutputManager::openOutfile() {
	if(rootOut) { rootOut->Close(); }
	makePath(rootPath);
	rootOut = new TFile((rootPath+"/"+name+".root").c_str(),"RECREATE");
	rootOut->cd();	
}

void OutputManager::writeROOT() {
	printf("\n--------- Building output .root file... ----------\n");
	if(!rootOut) openOutfile();
	rootOut->cd();
	writeItems();
	clearItems();
	rootOut->Close();
	rootOut = NULL;
	printf("---------          Done.          ----------\n");
}


TH1F* OutputManager::registeredTH1F(std::string hname, std::string htitle, unsigned int nbins, float x0, float x1) {
	if(rootOut) rootOut->cd();
	return (TH1F*)addObject(new TH1F(hname.c_str(),htitle.c_str(),nbins,x0,x1));
}

TH2F* OutputManager::registeredTH2F(std::string hname, std::string htitle, unsigned int nbinsx, float x0, float x1, unsigned int nbinsy, float y0, float y1) {
	if(rootOut) rootOut->cd();
	return (TH2F*)addObject(new TH2F(hname.c_str(),htitle.c_str(),nbinsx,x0,x1,nbinsy,y0,y1));
}

void OutputManager::printCanvas(std::string fname, std::string suffix) const {
	makePath(plotPath+"/"+fname+suffix,true);
	printf("Printing canvas '%s' in '%s'\n",(fname+suffix).c_str(), plotPath.c_str());
	defaultCanvas->Print((plotPath+"/"+fname+suffix).c_str());
}

