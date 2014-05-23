#include "G4toPMT.hh"
#include "PenelopeToPMT.hh"
#include "OutputManager.hh"
#include "PathUtils.hh"
#include "GraphicsUtils.hh"
#include "MC_Comparisons.hh"
#include <TROOT.h>
#include <TStyle.h>
#include <TLegend.h>
#include <iostream>

void mc_compare_plots(OutputManager& OM, Sim2PMT& SP1, Sim2PMT& SP2, double emax) {
	
	// put two simulations in a list for easy iteration
	sps.push_back(&SP1);
	sps.push_back(&SP2);
	// process each simulation to fill histograms
	for(unsigned int i=0; i<sps.size(); i++) {
		SP = sps[i];
		// define histograms
		for(EventType t=TYPE_0_EVENT; t <= TYPE_II_EVENT; ++t) {
			hEvis[t].push_back(OM.registeredTH1F("hEvis_"+itos(t)+"_"+itos(i),
												 "Type "+itos(t)+" Energy Spectrum",
												 100,0,emax));
			hEvisd[t].push_back(OM.registeredTH1F("hEvisd_"+itos(t)+"_"+itos(i),
                                                                                                 "Type "+itos(t)+" Raw energy loss in scint. dead region",
                                                                                                 100,0,emax));
			hEquench[t].push_back(OM.registeredTH1F("hEquench_"+itos(t)+"_"+itos(i),
												 "Type "+itos(t)+" Quenched Energy Spectrum",
												 100,0,emax));
			hEMWPC[t].push_back(OM.registeredTH1F("hEMWPC_"+itos(t)+"_"+itos(i),
                                                                                                 "Type "+itos(t)+" Energy loss in MWPC gas",
                                                                                                 100,0,150.0));
			hEMWPCd[t].push_back(OM.registeredTH1F("hEMWPCd_"+itos(t)+"_"+itos(i),
                                                                                                 "Type "+itos(t)+" Energy loss in MWPC gas dead regions",
                                                                                                 100,0,100.0));
			hEMylarf[t].push_back(OM.registeredTH1F("hEMylarf_"+itos(t)+"_"+itos(i),
                                                                                                 "Type "+itos(t)+" Energy loss in front MWPC mylar",
                                                                                                 100,0,150.0));
			hEMylarb[t].push_back(OM.registeredTH1F("hEMylarb_"+itos(t)+"_"+itos(i),
                                                                                                 "Type "+itos(t)+" Energy loss in back MWPC mylar",
                                                                                                 100,0,150.0));
			hEFoils[t].push_back(OM.registeredTH1F("hEFoils_"+itos(t)+"_"+itos(i),
                                                                                                 "Type "+itos(t)+" Energy loss in decay tube endcap foils",
                                                                                                 100,0,50));
			hEWires[t].push_back(OM.registeredTH1F("hEWires_"+itos(t)+"_"+itos(i),
                                                                                                 "Type "+itos(t)+" Energy deposited in wire planes",
                                                                                                 25,0,50.0));
		}
		hCTScint.push_back(OM.registeredTH1F("hCTScint_"+itos(i),"cos(theta) entering scintillator",100,0.0,1.0));
		hCTMWPC.push_back(OM.registeredTH1F("hCTMWPC_"+itos(i),"cos(theta) entering MWPC",100,0.0,1.0));
		hPScintX.push_back(OM.registeredTH1F("hPScintX_"+itos(i),"scintillator trigger position (x)",100,-100.0,100.0));
                hPMWPCX.push_back(OM.registeredTH1F("hPMWPCX_"+itos(i),"MWPC x-position",100,-80.0,80.0));
		hPScintY.push_back(OM.registeredTH1F("hPScintY_"+itos(i),"scintillator trigger position (y)",100,-100.0,100.0));
                hPMWPCY.push_back(OM.registeredTH1F("hPMWPCY_"+itos(i),"MWPC y-position",100,-80.0,80.0));
		hCTMWPCo.push_back(OM.registeredTH1F("hCTMWPCo_"+itos(i),"cos(theta) exiting MWPC",100,0.0,1.0));		
		hCTFoil.push_back(OM.registeredTH1F("hCTFoil_"+itos(i),"cos(theta) entering decay trap foil",100,0.0,1.0));		

		// fill histograms
		SP->startScan();
		while(SP->nextPoint()) {
			if(SP->fType == TYPE_III_EVENT) SP->fType = TYPE_II_EVENT; // merge Type II/III
			if(SP->fPID != PID_BETA || SP->fSide != EAST) continue;
			if(SP->fType <= TYPE_II_EVENT) {
				hEvis[SP->fType].back()->Fill(SP->eDep[EAST]+SP->eDep[WEST],SP->physicsWeight);
				hEvisd[SP->fType].back()->Fill(SP->edepDeadScint[EAST]+SP->edepDeadScint[WEST],SP->physicsWeight);
				hEquench[SP->fType].back()->Fill(SP->eQ[EAST]+SP->eQ[WEST],SP->physicsWeight);
				hEMWPC[SP->fType].back()->Fill(SP->eW[EAST]+SP->eW[WEST],SP->physicsWeight);
				hEMWPCd[SP->fType].back()->Fill(SP->edepDeadMWPC[EAST]+SP->edepDeadMWPC[WEST],SP->physicsWeight);
				hEMylarf[SP->fType].back()->Fill(SP->edepWinIn[EAST]+SP->edepWinIn[WEST],SP->physicsWeight);
				hEMylarb[SP->fType].back()->Fill(SP->edepWinOut[EAST]+SP->edepWinOut[WEST],SP->physicsWeight);
				hEFoils[SP->fType].back()->Fill(SP->edepFoils[EAST]+SP->edepFoils[WEST],SP->physicsWeight);
				hEWires[SP->fType].back()->Fill(SP->edepWires[EAST]+SP->edepWires[WEST],SP->physicsWeight);
				hPScintX.back()->Fill(SP->scintPos[EAST][X_DIRECTION],SP->physicsWeight);
                       	 	hPMWPCX.back()->Fill(SP->mwpcPos[EAST][X_DIRECTION],SP->physicsWeight);
                        	hPScintY.back()->Fill(SP->scintPos[EAST][Y_DIRECTION],SP->physicsWeight);
                        	hPMWPCY.back()->Fill(SP->mwpcPos[EAST][Y_DIRECTION],SP->physicsWeight);
			}
			if(SP->fSide <= WEST) {
				if(SP->cosThetaInScint[SP->fSide] != 0){
					hCTScint.back()->Fill(SP->cosThetaInScint[SP->fSide],SP->physicsWeight);
				}
				if(SP->cosThetaInWinOut[SP->fSide] != 0){
					hCTMWPCo.back()->Fill(SP->cosThetaInWinOut[SP->fSide],SP->physicsWeight);
				}
				if(SP->cosThetaInFoils[SP->fSide] != 0){
					hCTFoil.back()->Fill(SP->cosThetaInFoils[SP->fSide],SP->physicsWeight);
				}
				if(SP->cosThetaInWinIn[SP->fSide] != 0){
					hCTMWPC.back()->Fill(SP->cosThetaInWinIn[SP->fSide],SP->physicsWeight);
				}
			}
		}
		
		// record normalization counts
		t0norm.push_back(hEvis[TYPE_0_EVENT].back()->Integral());
		scintnorm.push_back(hCTScint.back()->Integral());
		mwpcoutnorm.push_back(hCTMWPCo.back()->Integral());
		mwpcnorm.push_back(hCTMWPC.back()->Integral());
                foilnorm.push_back(hCTFoil.back()->Integral());
	}
	
	/////////////
	// make plots, process data
	/////////////
	
	// normaliztion factors to Type 0 counts
	for(unsigned int i=0; i<t0norm.size(); i++){
		t0norm[i] = 1.0/t0norm[i];
		mwpcnorm[i] = 1.0/mwpcnorm[i];
		foilnorm[i] = 1.0/foilnorm[i];
		mwpcoutnorm[i] = 1.0/mwpcoutnorm[i];
		scintnorm[i] = 1.0/scintnorm[i];
	}
	std::vector<TH1*> hToPlot;
	OM.defaultCanvas->cd();
        //gPad->SetLogy();	
	for(EventType t=TYPE_0_EVENT; t <= TYPE_II_EVENT; ++t) {
		hToPlot.clear();
		for(unsigned int i=0; i<sps.size(); i++) {
			Stringmap m;
			m.insert("type",t);
			m.insert("sim",i);
			m.insert("mean",hEvis[t][i]->GetMean());
			m.insert("rms",hEvis[t][i]->GetRMS());
			m.insert("counts",hEvis[t][i]->Integral());
			m.insert("normcounts",hEvis[t][i]->Integral()*t0norm[i]);
			gPad->SetLogy(1);
			hEvis[t][i]->GetXaxis()->SetTitle("E (keV)");
			hEvis[t][i]->SetLineColor(2+2*i);
			hEvis[t][i]->Scale(t0norm[i]/hEvis[t][i]->GetBinWidth(1));
			hToPlot.push_back(hEvis[t][i]);
			OM.qOut.insert("evis",m);
		}
		drawSimulHistos(hToPlot);
		TLegend *evleg = new TLegend(0.13,0.82,0.4,0.9);
                evleg->AddEntry(hEvis[t][0],"Geant4","l");
                evleg->AddEntry(hEvis[t][1],"Penelope","l");
		evleg->Draw();
		OM.printCanvas("Evis_Type_"+itos(t));
		
		hToPlot.clear();
                for(unsigned int i=0; i<sps.size(); i++) {
                        Stringmap m;
                        m.insert("type",t);
                        m.insert("sim",i);
                        m.insert("mean",hEvisd[t][i]->GetMean());
                        m.insert("rms",hEvisd[t][i]->GetRMS());
                        m.insert("counts",hEvisd[t][i]->Integral());
                        m.insert("normcounts",hEvisd[t][i]->Integral()*t0norm[i]);
                        hEvisd[t][i]->GetXaxis()->SetTitle("E (keV)");
                        hEvisd[t][i]->SetLineColor(2+2*i);
                        hEvisd[t][i]->Scale(t0norm[i]/hEvisd[t][i]->GetBinWidth(1));
                        hToPlot.push_back(hEvisd[t][i]);
                        OM.qOut.insert("evisdead",m);
                }
                drawSimulHistos(hToPlot);
                TLegend *evdleg = new TLegend(0.13,0.82,0.4,0.9);
                evdleg->AddEntry(hEvisd[t][0],"Geant4","l");
                evdleg->AddEntry(hEvisd[t][1],"Penelope","l");
                evdleg->Draw();
                OM.printCanvas("Evisdead_Type_"+itos(t));

		hToPlot.clear();
		for(unsigned int i=0; i<sps.size(); i++) {
			Stringmap m;
			m.insert("type",t);
			m.insert("sim",i);
			m.insert("mean",hEquench[t][i]->GetMean());
			m.insert("rms",hEquench[t][i]->GetRMS());
			m.insert("counts",hEquench[t][i]->Integral());
			m.insert("normcounts",hEquench[t][i]->Integral()*t0norm[i]);
			hEquench[t][i]->GetXaxis()->SetTitle("E (keV)");
			hEquench[t][i]->SetLineColor(2+2*i);
			hEquench[t][i]->Scale(t0norm[i]/hEquench[t][i]->GetBinWidth(1));
			hToPlot.push_back(hEquench[t][i]);
			OM.qOut.insert("equench",m);
		}
		drawSimulHistos(hToPlot);
		TLegend *eqleg = new TLegend(0.13,0.82,0.4,0.9);
                eqleg->AddEntry(hEquench[t][0],"Geant4","l");
                eqleg->AddEntry(hEquench[t][1],"Penelope","l");
                eqleg->Draw();
		OM.printCanvas("Equench_Type_"+itos(t));

		hToPlot.clear();
                for(unsigned int i=0; i<sps.size(); i++) {
                        Stringmap m;
                        m.insert("type",t);
                        m.insert("sim",i);
                        m.insert("mean",hEMWPC[t][i]->GetMean());
                        m.insert("rms",hEMWPC[t][i]->GetRMS());
                        m.insert("counts",hEMWPC[t][i]->Integral());
                        m.insert("normcounts",hEMWPC[t][i]->Integral()*t0norm[i]);
                        hEMWPC[t][i]->GetXaxis()->SetTitle("E (keV)");
                        hEMWPC[t][i]->SetLineColor(2+2*i);
                        hEMWPC[t][i]->Scale(t0norm[i]/hEMWPC[t][i]->GetBinWidth(1));
                        hToPlot.push_back(hEMWPC[t][i]);
                        OM.qOut.insert("emwpc",m);
                }
                drawSimulHistos(hToPlot);
                TLegend *emleg = new TLegend(0.13,0.82,0.4,0.9);
                emleg->AddEntry(hEMWPC[t][0],"Geant4","l");
                emleg->AddEntry(hEMWPC[t][1],"Penelope","l");
                emleg->Draw();
                OM.printCanvas("EMWPCgas_Type_"+itos(t));

		hToPlot.clear();
                for(unsigned int i=0; i<sps.size(); i++) {
                        Stringmap m;
                        m.insert("type",t);
                        m.insert("sim",i);
                        m.insert("mean",hEMWPCd[t][i]->GetMean());
                        m.insert("rms",hEMWPCd[t][i]->GetRMS());
                        m.insert("counts",hEMWPCd[t][i]->Integral());
                        m.insert("normcounts",hEMWPCd[t][i]->Integral()*t0norm[i]);
                        hEMWPCd[t][i]->GetXaxis()->SetTitle("E (keV)");
                        hEMWPCd[t][i]->SetLineColor(2+2*i);
                        hEMWPCd[t][i]->Scale(t0norm[i]/hEMWPCd[t][i]->GetBinWidth(1));
                        hToPlot.push_back(hEMWPCd[t][i]);
                        OM.qOut.insert("emwpcdead",m);
                }
                drawSimulHistos(hToPlot);
                TLegend *emdleg = new TLegend(0.13,0.82,0.4,0.9);
                emdleg->AddEntry(hEMWPCd[t][0],"Geant4","l");
                emdleg->AddEntry(hEMWPCd[t][1],"Penelope","l");
                emdleg->Draw();
                OM.printCanvas("EMWPCgasd_Type_"+itos(t));

		hToPlot.clear();
                for(unsigned int i=0; i<sps.size(); i++) {
                        Stringmap m;
                        m.insert("type",t);
                        m.insert("sim",i);
                        m.insert("mean",hEMylarf[t][i]->GetMean());
                        m.insert("rms",hEMylarf[t][i]->GetRMS());
                        m.insert("counts",hEMylarf[t][i]->Integral());
                        m.insert("normcounts",hEMylarf[t][i]->Integral()*t0norm[i]);
                        hEMylarf[t][i]->GetXaxis()->SetTitle("E (keV)");
                        hEMylarf[t][i]->SetLineColor(2+2*i);
                        hEMylarf[t][i]->Scale(t0norm[i]/hEMylarf[t][i]->GetBinWidth(1));
                        hToPlot.push_back(hEMylarf[t][i]);
                        OM.qOut.insert("emylarf",m);
                }
                drawSimulHistos(hToPlot);
                TLegend *emyfleg = new TLegend(0.13,0.82,0.4,0.9);
                emyfleg->AddEntry(hEMylarf[t][0],"Geant4","l");
                emyfleg->AddEntry(hEMylarf[t][1],"Penelope","l");
                emyfleg->Draw();
                OM.printCanvas("EMylarfront_Type_"+itos(t));

		hToPlot.clear();
                for(unsigned int i=0; i<sps.size(); i++) {
                        Stringmap m;
                        m.insert("type",t);
                        m.insert("sim",i);
                        m.insert("mean",hEMylarb[t][i]->GetMean());
                        m.insert("rms",hEMylarb[t][i]->GetRMS());
                        m.insert("counts",hEMylarb[t][i]->Integral());
                        m.insert("normcounts",hEMylarb[t][i]->Integral()*t0norm[i]);
                        hEMylarb[t][i]->GetXaxis()->SetTitle("E (keV)");
                        hEMylarb[t][i]->SetLineColor(2+2*i);
                        hEMylarb[t][i]->Scale(t0norm[i]/hEMylarb[t][i]->GetBinWidth(1));
                        hToPlot.push_back(hEMylarb[t][i]);
                        OM.qOut.insert("emylarb",m);
                }
                drawSimulHistos(hToPlot);
                TLegend *emybleg = new TLegend(0.13,0.82,0.4,0.9);
                emybleg->AddEntry(hEMylarb[t][0],"Geant4","l");
                emybleg->AddEntry(hEMylarb[t][1],"Penelope","l");
                emybleg->Draw();
                OM.printCanvas("EMylarback_Type_"+itos(t));

		hToPlot.clear();
                for(unsigned int i=0; i<sps.size(); i++) {
                        Stringmap m;
                        m.insert("type",t);
                        m.insert("sim",i);
                        m.insert("mean",hEFoils[t][i]->GetMean());
                        m.insert("rms",hEFoils[t][i]->GetRMS());
                        m.insert("counts",hEFoils[t][i]->Integral());
                        m.insert("normcounts",hEFoils[t][i]->Integral()*t0norm[i]);
                        hEFoils[t][i]->GetXaxis()->SetTitle("E (keV)");
                        hEFoils[t][i]->SetLineColor(2+2*i);
                        hEFoils[t][i]->Scale(t0norm[i]/hEFoils[t][i]->GetBinWidth(1));
                        hToPlot.push_back(hEFoils[t][i]);
                        OM.qOut.insert("efoils",m);
                }
                drawSimulHistos(hToPlot);
                TLegend *eefleg = new TLegend(0.13,0.82,0.4,0.9);
                eefleg->AddEntry(hEFoils[t][0],"Geant4","l");
                eefleg->AddEntry(hEFoils[t][1],"Penelope","l");
                eefleg->Draw();
                OM.printCanvas("EFoils_Type_"+itos(t));

		hToPlot.clear();
                for(unsigned int i=0; i<sps.size(); i++) {
                        Stringmap m;
                        m.insert("type",t);
                        m.insert("sim",i);
                        m.insert("mean",hEWires[t][i]->GetMean());
                        m.insert("rms",hEWires[t][i]->GetRMS());
                        m.insert("counts",hEWires[t][i]->Integral());
                        m.insert("normcounts",hEWires[t][i]->Integral()*t0norm[i]);
			gPad->SetLogy(1);
			hEWires[t][i]->GetXaxis()->SetTitle("E (keV)");
                        hEWires[t][i]->SetLineColor(2+2*i);
                        hEWires[t][i]->Scale(t0norm[i]/hEWires[t][i]->GetBinWidth(1));
                        //gPad->SetLogy(1);
			hToPlot.push_back(hEWires[t][i]);
                        OM.qOut.insert("ewires",m);
                }
                drawSimulHistos(hToPlot);
                TLegend *ewleg = new TLegend(0.13,0.82,0.4,0.9);
                ewleg->AddEntry(hEWires[t][0],"Geant4","l");
                ewleg->AddEntry(hEWires[t][1],"Penelope","l");
                ewleg->Draw();
                OM.printCanvas("EWires_Type_"+itos(t));
	}
	
	hToPlot.clear();
	for(unsigned int i=0; i<sps.size(); i++) {
		hCTScint[i]->SetLineColor(2+2*i);
		hCTScint[i]->Scale(scintnorm[i]/hCTScint[i]->GetBinWidth(1));
		hToPlot.push_back(hCTScint[i]);
	}
	drawSimulHistos(hToPlot);
	TLegend *ctsleg = new TLegend(0.13,0.82,0.4,0.9);
        ctsleg->AddEntry(hCTScint[0],"Geant4","l");
        ctsleg->AddEntry(hCTScint[1],"Penelope","l");
        ctsleg->Draw();
	OM.printCanvas("Scint_Costheta");
	
	hToPlot.clear();
	for(unsigned int i=0; i<sps.size(); i++) {
		hCTMWPC[i]->SetLineColor(2+2*i);
		hCTMWPC[i]->Scale(mwpcnorm[i]/hCTMWPC[i]->GetBinWidth(1));
		hToPlot.push_back(hCTMWPC[i]);
	}
	drawSimulHistos(hToPlot);
 	TLegend *ctmleg = new TLegend(0.13,0.82,0.4,0.9);
        ctmleg->AddEntry(hCTMWPC[0],"Geant4","l");
        ctmleg->AddEntry(hCTMWPC[1],"Penelope","l");
        ctmleg->Draw();
	OM.printCanvas("MWPC_In_Costheta");

	hToPlot.clear();
        for(unsigned int i=0; i<sps.size(); i++) {
                hCTMWPCo[i]->SetLineColor(2+2*i);
                hCTMWPCo[i]->Scale(mwpcoutnorm[i]/hCTMWPCo[i]->GetBinWidth(1));
                hToPlot.push_back(hCTMWPCo[i]);
        }
        drawSimulHistos(hToPlot);
        TLegend *ctmlego = new TLegend(0.13,0.82,0.4,0.9);
        ctmlego->AddEntry(hCTMWPCo[0],"Geant4","l");
        ctmlego->AddEntry(hCTMWPCo[1],"Penelope","l");
        ctmlego->Draw();
        OM.printCanvas("MWPC_Out_Costheta");

	hToPlot.clear();
        for(unsigned int i=0; i<sps.size(); i++) {
                hPScintX[i]->GetXaxis()->SetTitle("Pos (mm)");
		hPScintX[i]->SetLineColor(2+2*i);
                hPScintX[i]->Scale(t0norm[i]/hPScintX[i]->GetBinWidth(1));
                hToPlot.push_back(hPScintX[i]);
        }
        drawSimulHistos(hToPlot);
        TLegend *psxleg = new TLegend(0.13,0.82,0.4,0.9);
        psxleg->AddEntry(hPScintX[0],"Geant4","l");
        psxleg->AddEntry(hPScintX[1],"Penelope","l");
        psxleg->Draw();
        OM.printCanvas("TrigposX_Scint");
	
	hToPlot.clear();
        for(unsigned int i=0; i<sps.size(); i++) {
		hPMWPCX[i]->GetXaxis()->SetTitle("Pos (mm)");
                hPMWPCX[i]->SetLineColor(2+2*i);
                hPMWPCX[i]->Scale(t0norm[i]/hPMWPCX[i]->GetBinWidth(1));
                hToPlot.push_back(hPMWPCX[i]);
        }
        drawSimulHistos(hToPlot);
        TLegend *pmxleg = new TLegend(0.13,0.82,0.4,0.9);
        pmxleg->AddEntry(hPMWPCX[0],"Geant4","l");
        pmxleg->AddEntry(hPMWPCX[1],"Penelope","l");
        pmxleg->Draw();
        OM.printCanvas("TrigposX_MWPC");

	hToPlot.clear();
        for(unsigned int i=0; i<sps.size(); i++) {
		hPScintY[i]->GetXaxis()->SetTitle("Pos (mm)");
                hPScintY[i]->SetLineColor(2+2*i);
                hPScintY[i]->Scale(t0norm[i]/hPScintY[i]->GetBinWidth(1));
                hToPlot.push_back(hPScintY[i]);
        }
        drawSimulHistos(hToPlot);
        TLegend *psyleg = new TLegend(0.13,0.82,0.4,0.9);
        psyleg->AddEntry(hPScintY[0],"Geant4","l");
        psyleg->AddEntry(hPScintY[1],"Penelope","l");
        psyleg->Draw();
        OM.printCanvas("TrigposY_Scint");

        hToPlot.clear();
        for(unsigned int i=0; i<sps.size(); i++) {
		hPMWPCY[i]->GetXaxis()->SetTitle("Pos (mm)");
                hPMWPCY[i]->SetLineColor(2+2*i);
                hPMWPCY[i]->Scale(t0norm[i]/hPMWPCY[i]->GetBinWidth(1));
                hToPlot.push_back(hPMWPCY[i]);
        }
        drawSimulHistos(hToPlot);
        TLegend *pmyleg = new TLegend(0.13,0.82,0.4,0.9);
        pmyleg->AddEntry(hPMWPCY[0],"Geant4","l");
        pmyleg->AddEntry(hPMWPCY[1],"Penelope","l");
        pmyleg->Draw();
        OM.printCanvas("TrigposY_MWPC");

	hToPlot.clear();
        for(unsigned int i=0; i<sps.size(); i++) {
                hCTFoil[i]->SetLineColor(2+2*i);
                hCTFoil[i]->Scale(foilnorm[i]/hCTFoil[i]->GetBinWidth(1));
                hToPlot.push_back(hCTFoil[i]);
        }
        drawSimulHistos(hToPlot);
        TLegend *ctfleg = new TLegend(0.13,0.82,0.4,0.9);
        ctfleg->AddEntry(hCTFoil[0],"Geant4","l");
        ctfleg->AddEntry(hCTFoil[1],"Penelope","l");
        ctfleg->Draw();
        OM.printCanvas("Foils_Costheta");
}

int main(int, char **) {

	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1);
	gStyle->SetNumberContours(255);
	
	// list of simulated energies
	//int enlist[] = {50,100,150,200,300,400,600,800};
	
	//for(int i = 1; i < 8; i++) {
	//	int l = enlist[i];
		
		// set up output directories
		//OutputManager OM("MC_Compare",getEnvSafe("UCNA_ANA_PLOTS")+"/test/MC_Compare/"+itos(l)+"_keV");
		OutputManager OM("MC_Compare",getEnvSafe("UCNA_ANA_PLOTS")+"/test/MC_Compare/Xenon131");
                Stringmap mcdat;
		//mcdat.insert("energy",l);
		mcdat.insert("n_MC",2);
		
		// make PMT response look like run 15925
		PMTCalibrator PCal(19891);
	
		// load Geant4 simulation
		G4toPMT g2p(true);
		//std::string fname = getEnvSafe("G4OUTDIR")+"/IsotLine_eGunRandMomentum_"+itos(l)+".0keV/analyzed_*";
		std::string fname = getEnvSafe("G4OUTDIR")+"/geant4_data_MB/NEW_Xe131_11-2-/analyzed*";
                mcdat.insert("MC_1",fname);
		g2p.addFile(fname);
		//if(!g2p.getnFiles()) continue;
		g2p.setCalibrator(PCal);
	
		// load Penelope simulation
                PenelopeToPMT p2p;
                //fname = "/home/ucna/penelope_output/iso_line_sources/event_"+itos(l/50-1)+" _*.root";
                std::string fpname = getEnvSafe("PENOUTDIR")+"/Xenon131_001*";
                mcdat.insert("MC_2",fpname);
                p2p.addFile(fpname);
                //if(!p2p.getnFiles()) continue;
                p2p.setCalibrator(PCal);
	
		// make comparison plots between simulations
		//mc_compare_plots(OM,g2p,p2p,2*l);
		mc_compare_plots(OM,g2p,p2p,200);
                OM.qOut.insert("mcdat",mcdat);
		OM.write();
		OM.setWriteRoot(true);
	//}
	return 0;
}
