/*
 * TEana.h
 *
 *  Created on: Jun 7, 2016
 *      Author: klein
 *      * came from  Makeclass
 *      * will calculate the calibration constant
 */

#ifndef TEANA_H_
#define TEANA_H_

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.

class TEana {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        FreqCenter;
   Double_t        FreqStep;
   Double_t        Temperature;
   Double_t        ScanPoints;
   Double_t        TuneV;
   Double_t        Offset;
   Double_t        ControllerV;
   vector<double>  *array;
   Long64_t        timel;
   Int_t           IntScanPoints;

   // List of branches
   TBranch        *b_FreqCenter;   //!
   TBranch        *b_FreqStep;   //!
   TBranch        *b_Temperature;   //!
   TBranch        *b_ScanPoints;   //!
   TBranch        *b_TuneV;   //!
   TBranch        *b_Offset;   //!
   TBranch        *b_ControllerV;   //!
   TBranch        *b_array;   //!
   TBranch        *b_timel;   //!
   TBranch        *b_IntScanPoints;   //!

   TEana();
   virtual ~TEana();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual Int_t 	OpenFile(TString);
};









#endif /* TEANA_H_ */
#ifdef TEana_cxx



TEana::TEana()
{

		cout<<"**************       Initialize TE analyzer **************** \n\n\n";
}

int TEana::OpenFile(TString rootfile){

	// oepn file and initialize tree
     cout<<"opening file "<<rootfile<<"\n";

     f = new TFile(rootfile);

     f->GetObject("NMRtree",tree);
     Init(tree);

   return 0;

}

TEana::~TEana()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t TEana::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TEana::LoadTree(Long64_t entry)
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

void TEana::Init(TTree *tree)
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

   fChain->SetBranchAddress("FreqCenter", &FreqCenter, &b_FreqCenter);
   fChain->SetBranchAddress("FreqStep", &FreqStep, &b_FreqStep);
   fChain->SetBranchAddress("Temperature", &Temperature, &b_Temperature);
   fChain->SetBranchAddress("ScanPoints", &ScanPoints, &b_ScanPoints);
   fChain->SetBranchAddress("TuneV", &TuneV, &b_TuneV);
   fChain->SetBranchAddress("Offset", &Offset, &b_Offset);
   fChain->SetBranchAddress("ControllerV", &ControllerV, &b_ControllerV);
   fChain->SetBranchAddress("array", &array, &b_array);
   fChain->SetBranchAddress("timel", &timel, &b_timel);
   fChain->SetBranchAddress("IntScanPoints", &IntScanPoints, &b_IntScanPoints);
   Notify();
}

Bool_t TEana::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}



void TEana::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L TEana.C
//      Root > TEana t
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
   }
}

#endif // #ifdef TEana_cxx


