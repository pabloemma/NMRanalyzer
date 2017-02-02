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


#include <Math/Factory.h>
#include <Math/Minimizer.h>
#include <Math/Functor.h>

class NMR_DST
{
public:

NMR_DST();
~NMR_DST();

	void OpenFile(TString);
	void CreateTree();
	void AddBranches();
	void CloseFile();
	void WriteTree();
	void FillTree();

	Int_t area;


private:
	Int_t BranchCounter; // this keeps track of the number of branches we will have
    std::string DST_pr = "DST> ";
    TFile *f;
    TTree *Dtree;


};

void	NMR_DST::OpenFile(TString DSTfile){

		     cout<<DST_pr<<"opening root outputfile "<<DSTfile<<"\n";

		     f = new TFile(DSTfile);
		     if(!f->IsOpen()){
		    	 cout <<DST_pr<<"Cannot open DST File "<<DSTfile<<" \n  will exit \n";
		    	 exit(EXIT_FAILURE);
		     }


		}


void    NMR_DST::CreateTree(){
    	//
    	 //this creates the DST tree
    	Dtree = new TTree("DTree","analyzed variables from NMRanalyzer");

     }

void	NMR_DST::AddBranches()	{
		Dtree->Branch("area",&area,"area/I");

	}
void	NMR_DST::FillTree(){
		Dtree->Fill();
		cout<<DST_pr<<area<<endl;

	}
void	NMR_DST::CloseFile(){
		Dtree->Write();
		f->Close();

	}

#endif
