/*
 * NMRana.h
 *
 *  Created on: May 10, 2016
 *      Author: klein
 */

//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue May 10 11:02:55 2016 by ROOT version 5.34/34
// from TTree NMRtree/NMR
// found on file: /Volumes/FastDisk/NMR/april22/ana_folder/NMR3544192508.root
//////////////////////////////////////////////////////////

#ifndef NMRana_h
#define NMRana_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>


// Header file for the classes stored in the TTree if any.
#include <vector>
#include <iostream>

using namespace std;

// Fixed size dimensions of array or collections stored in the TTree if any.

class NMRana {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           IntScanPoints;
   Double_t        FreqStep;
   Double_t        FreqCenter;
   Double_t        Temperature;
   Double_t        ScanPoints;
   Double_t        TuneV;
   Double_t        Offset;
   Double_t        ControllerV;
   Char_t          time[11];
   std::vector<double>  *array;

   Double_t MinFreq ; // limits of frequency sweep
   Double_t MaxFreq ;
   std::vector<TString> RootFileArray ; // if there is a list of input files it will put them into vector



   // List of branches
   TBranch        *b_IntScanPoints;   //!
   TBranch        *b_FreqStep;   //!
   TBranch        *b_FreqCenter;   //!

   TBranch        *b_Temperature;   //!
   TBranch        *b_ScanPoints;   //!
   TBranch        *b_TuneV;   //!
   TBranch        *b_Offset;   //!
   TBranch        *b_ControllerV;   //!
   TBranch        *b_time;   //!
   TBranch        *b_array;   //!
   TFile		  *f;

   TTree 		   *tree;
   TH1D 	 *NMR1; // Signal histogram
   TCanvas	 *c1;
   TChain    *NMRchain; // if we have more than one root file





   NMRana();
   virtual ~NMRana();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual Int_t 	OpenFile(TString);
   virtual void	    SetupHistos();
   virtual void		CloseFile();
   virtual void		DrawHistos();
   virtual int      OpenChain(std::vector<TString> );
};

#endif

#ifdef NMRana_cxx
NMRana::NMRana(){

}

int NMRana::OpenFile(TString rootfile){

	// oepn file and initialize tree
     cout<<"opening file "<<rootfile<<"\n";

     f = new TFile(rootfile);

     f->GetObject("NMRtree",tree);
     Init(tree);

   return 0;

}
int NMRana::OpenChain(std::vector<TString> RootFileArray){

	// This creates a chain fo trees instead of just one

	 NMRchain = new TChain("NMRtree");
	 // Now loop over all the rootfiles we have
		for(Int_t pos = 0 ; pos < RootFileArray.size() ; pos++)
		{
			cout<<RootFileArray[pos]<<"   filename \n";
			NMRchain->Add(RootFileArray[pos]);
		}



     //NMRchain->GetObject("NMRtree",tree);
     Init(NMRchain);

   return 0;

}

void NMRana::CloseFile(){
	f->Close();
}
NMRana::~NMRana()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t NMRana::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t NMRana::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void NMRana::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   array = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("IntScanPoints", &IntScanPoints, &b_IntScanPoints);
   fChain->SetBranchAddress("FreqStep", &FreqStep, &b_FreqStep);
   fChain->SetBranchAddress("Temperature", &Temperature, &b_Temperature);
   fChain->SetBranchAddress("ScanPoints", &ScanPoints, &b_ScanPoints);
   fChain->SetBranchAddress("ControllerV", &ControllerV, &b_ControllerV);
   fChain->SetBranchAddress("TuneV", &TuneV, &b_TuneV);
   fChain->SetBranchAddress("Offset", &Offset, &b_Offset);
//    fChain->SetBranchAddress("ControllerV", &ControllerV, &b_ControllerV);
   fChain->SetBranchAddress("time", time, &b_time);
   fChain->SetBranchAddress("array", &array, &b_array);
   Notify();
}

Bool_t NMRana::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void NMRana::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t NMRana::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
void NMRana::SetupHistos(){
// Here we setup histos if needed
// first we get the first entry to calculate the limits
/*	   if (fChain == 0) return;

	   Long64_t nentries = fChain->GetEntriesFast();

	   Long64_t nbytes = 0, nb = 0;

	   Long64_t ientry = LoadTree(0);

	   fChain->GetEntry(jentry);

	   fchain->Show(0);
*/
	   GetEntry(0);
	   Show(0);

	   // histo limits
	   MinFreq = FreqCenter - (ScanPoints-1)/2 * FreqStep;
	   MaxFreq = FreqCenter + (ScanPoints-1)/2 * FreqStep;

	   NMR1 = new TH1D("NMR1","Signal histogram",IntScanPoints,MinFreq,MaxFreq);


}

void NMRana::DrawHistos(){
// draw histos, mainly for debug purpose
	c1 = new TCanvas("c1","NMR signal",200,200,1000,800);
	NMR1->Draw();
	c1->Update();
}
void NMRana::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L NMRana.C
//      Root > NMRana t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      //	now fill histogram
      // reset freq_temp to lower bound
      Double_t freq_temp = MinFreq;

      for (UInt_t j = 0; j < array->size(); ++j) {
          NMR1->Fill(freq_temp,array->at(j));
          freq_temp = freq_temp+FreqStep;
      	  }

   }
}

#endif // #ifdef NMRana_cxx













