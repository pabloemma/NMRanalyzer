/*
 * NMR_DST.h
 *
 *  Created on: Feb 2, 2017
 *      Author: klein
 */
#ifndef _NMR_DST_H
#define _NMR_DST_H

#include <iostream>
#include <vector>
#include <TH1D.h>
#include <TString.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>



#include <Math/Factory.h>
#include <Math/Minimizer.h>
#include <Math/Functor.h>

class NMR_DST
{
public:

NMR_DST();
~NMR_DST();

	TTree *OpenFile();
	TTree *CreateTree();
	void AddBranches();
	void CloseFile();
	TTree * WriteTree();
	void FillTree(Double_t ,Long64_t, TH1D *);

	Double_t area;    TFile *f1;
    TTree *Dtree;


private:
	Int_t BranchCounter; // this keeps track of the number of branches we will have
    std::string DST_pr = "DST> ";
    Long64_t time;
    TH1D *histo;


};

NMR_DST::NMR_DST(){

}

NMR_DST::~NMR_DST(){

}


TTree	*NMR_DST::OpenFile(){
		Dtree = CreateTree();
		AddBranches();
        return Dtree;

		}


TTree  *  NMR_DST::CreateTree(){
    	//
    	 //this creates the DST tree
	    cout<<DST_pr<<"at producing tree"<<endl;
    	Dtree = new TTree("Dtree","analyzed variables from NMRanalyzer");
    	//Dtree->Print();
    	return Dtree;

     }

void	NMR_DST::AddBranches()	{
	Dtree->Branch("area",&area,"area/D");
	Dtree->Branch("time",&time,"time/L");
	   histo = new TH1D("histo","Real TimeSignal Fit background subtracted",401,212.6,213.4);
	Dtree->Branch("histo","TH1D",&histo,128000,0);

		//Dtree->Print();

	}
void	NMR_DST::FillTree(Double_t ar ,Long64_t timel,TH1D *histol){
		area = ar;
		time = timel;
		histo = histol;
		Dtree->Fill();
		//Dtree->Show(0);

		//cout<<DST_pr<<area<<endl;

	}
void	NMR_DST::CloseFile(){
		f1->Close();

	}
TTree	* NMR_DST::WriteTree(){
		//f1->cd();
		//f1->ls();
	    cout<<DST_pr<<"writing tree"<<endl;
	    //Dtree->Print();
	    //Dtree->Show(50);
	    //Dtree->Write();
	   // Dtree->Write(f1->GetName());
		// Int_t test1 = f1->Write();
	    //f1->Close();
		//cout<<DST_pr<<test1<<"  "<<endl;
        return Dtree;
	}

#endif
