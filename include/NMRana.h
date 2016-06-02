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
#include <TDatime.h>
#include <TStyle.h>  // so we can use gStyle
#include <TTimeStamp.h>


// Header file for the classes stored in the TTree if any.
#include <vector>
#include <iostream>
#include <string>


#include "Ana.h"

using namespace std;

// Fixed size dimensions of array or collections stored in the TTree if any.

class NMRana : public Ana {  //inherit from analysis

private:
	Double_t LowArea_X; // lower bound for area calculation
	Double_t HighArea_X;// upper bound for area claculation
	Int_t low_id;		// determines lower and upper index of array for integration
	Int_t high_id;


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
   Long64_t         timel; // note this time is down to 100musec, in order to deal only on the second level, strip the last 4 digits
   Long64_t	        time_offset;
   std::vector<double>  *array;

   Double_t MinFreq ; // limits of frequency sweep
   Double_t MaxFreq ;
   Double_t low_x;  //area limits
   Double_t high_x;
   Double_t SignalArea ; // the summed area of every signal // not normalized yet
   Double_t SignalAreaNormalized ; // aka polarization
   std::vector<TString> RootFileArray ; // if there is a list of input files it will put them into vector
   time_t TimeStamp;  // timestamp in seconds on UNIX time 0
   Long64_t TimeStamp_usec;  // timestamp in 100 musec
   timespec root_time; // timespec for root time stamp

   	   char *timel_ptr; // because  Root stored the string as a charcater array

   // List of branches
   TBranch        *b_IntScanPoints;   //!
   TBranch        *b_FreqStep;   //!
   TBranch        *b_FreqCenter;   //!

   TBranch        *b_Temperature;   //!
   TBranch        *b_ScanPoints;   //!
   TBranch        *b_TuneV;   //!
   TBranch        *b_Offset;   //!
   TBranch        *b_ControllerV;   //!
   TBranch        *b_timel;   //!
   TBranch        *b_array;   //!
   TFile		  *f;

   TTree 		   *tree;
   // the block for the histograms
   TH1D 	 *NMR1; // Signal histogram
   TH1D		 *PolTime; // polarization vs time
   //
   TCanvas	 *GeneralCanvas;
   TCanvas	 *StripCanvas;
   TChain    *NMRchain; // if we have more than one root file
   TString timestring; // from labview time
   TDatime *td ; // Datetime for time histogram
   TTimeStamp *RootTimeStamp;






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
   virtual void		AreaSetLimits(Double_t , Double_t);
   virtual Double_t CalculateArea(std::vector<Double_t> *);
   virtual void		PrintTime();  // prints time from Labview time stamp
   virtual void		SetupCanvas();
   virtual TH1D * 	SetupStripChart(TString);
   virtual void 	GetTimeStamp();


};
#endif

#ifdef NMRana_cxx


NMRana::NMRana(){
	time_offset = 2082852019 ; // unix offset in seconds


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

   fChain->SetBranchAddress("FreqStep", &FreqStep, &b_FreqStep);
   fChain->SetBranchAddress("FreqCenter", &FreqCenter, &b_FreqCenter);
   fChain->SetBranchAddress("Temperature", &Temperature, &b_Temperature);
   fChain->SetBranchAddress("ScanPoints", &ScanPoints, &b_ScanPoints);
   fChain->SetBranchAddress("ControllerV", &ControllerV, &b_ControllerV);
   fChain->SetBranchAddress("TuneV", &TuneV, &b_TuneV);
   fChain->SetBranchAddress("Offset", &Offset, &b_Offset);
//    fChain->SetBranchAddress("ControllerV", &ControllerV, &b_ControllerV);
   fChain->SetBranchAddress("timel", &timel, &b_timel);
   fChain->SetBranchAddress("array", &array, &b_array);
   fChain->SetBranchAddress("IntScanPoints", &IntScanPoints, &b_IntScanPoints);
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
	   //gROOT->cd(); //this prevents from histos being deleted once files are closed
	   if (fChain == 0) return;

	   Long64_t nentries = fChain->GetEntriesFast();

	   Long64_t nbytes = 0, nb = 0;

	   Long64_t ientry = LoadTree(0);

	   fChain->GetEntry(0);

	   fChain->Show(0);


	   //GetEntry(0);
	   //Show(0);


	   // histo limits
	   MinFreq = (FreqCenter) - (ScanPoints-1.)/2. * FreqStep;
	   MaxFreq = FreqCenter + (ScanPoints-1.)/2. * FreqStep;

	   NMR1 = new TH1D("NMR1","Signal histogram",IntScanPoints,MinFreq,MaxFreq);

	   // Determine the Integration or summation limits for peak in terms of channels.
	   low_id = NMR1->GetXaxis()->FindBin(LowArea_X);
	   high_id = NMR1->GetXaxis()->FindBin(HighArea_X);

       //




	   // Now setup the time versus polarization histo
	   // the time is in 100 microseconds, which is much more precise than what we need
	   // so I didvide by 10000  with an integer division, also the time and date is now in UNIX time
	   TimeStamp_usec = timel -time_offset; //
	   //cout<<"timestamp   "<<TimeStamp<<"\n";
	   time_t timm = Int_t(TimeStamp_usec/10000);
	   td = new TDatime(timm);
	   // we set the time stamp to 0
	   gStyle->SetTimeOffset(td->Convert());
	   PolTime = SetupStripChart("Polarization vs time");
	   PolTime->SetMaximum(100.);
	   PolTime->SetMinimum(-100.);
	   PolTime->SetLineColor(2);




	   // Print original time
	   PrintTime();


	   Int_t timewindow = 644;
	   PolTime = new TH1D("PolTime","Polarization vs time",timewindow,(TimeStamp),(TimeStamp)+timewindow);
	   //PolTime = new TH1D("PolTime","Polarization vs time",timewindow,0,timewindow);
	   PolTime->GetXaxis()->SetTimeDisplay(1);

}
void NMRana::SetupCanvas(){
	// creates the different Canvas
	// master canvas for all histograms
	GeneralCanvas = new TCanvas("GeneralCanvas","NMR signal",200,200,1000,800);


	//canvas for all strip charts
	StripCanvas =  new TCanvas("StripCanvas","NMR strip charts",800,50,1000,600);
	StripCanvas->SetFillColor(42);
	StripCanvas->SetFrameFillColor(33);
	StripCanvas->SetGrid();

}

void NMRana::DrawHistos(){

	// analyze spectra
	FindPeak(NMR1);
// draw histos, mainly for debug purpose
	GeneralCanvas->Divide(1,2);
	GeneralCanvas->cd(1);
	NMR1->Draw();
    FitSpectrum(NMR1,1);
	GeneralCanvas->cd(2);
	PolTime->Draw();
	GeneralCanvas->Update();




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


   	   // go to strip chart
   StripCanvas->cd();

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t time_prev = 0;
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

//sum the peak area
      SignalArea = CalculateArea(array);
      // file polarization vs time
/*      Long64_t timediff =timel - time_prev;
      time_prev = timel;
      cout<<"time difference"<<timediff<< "\n";
//      TimeStamp = time_t((timel)/10000 -2082844800);
      TimeStamp_usec = timel - time_offset;
      TimeStamp = TimeStamp_usec/10000;
//      TimeStamp = timel-2082844800;
      cout<<"Timestamp"<<TimeStamp<<"\n";
      */
	  GetTimeStamp();
	  // RootTimeStamp->Print();
      //PolTime->Fill((TimeStamp),SignalArea*100.);
      //PolTime->Fill(jentry,SignalArea*100.);
	  PolTime->SetBinContent(jentry,SignalArea*100.);
	  PolTime->GetXaxis()->SetRange(jentry-500,jentry+500);
	  StripCanvas->Clear();
	  PolTime->Draw();
	  StripCanvas->Modified();
	  StripCanvas->Update();



//      cout<<timel<<"another one \n";

   }// end of entry loop
}

Double_t NMRana::CalculateArea(std::vector<Double_t> *array){
	// this function calculates the area of the NMR perak by simpy summing it in
	//  limits given by AreaSetLimits in the main program
	// the area is calculated as sum i_low to i_high (a(i)*binwidth)

	Double_t sum = 0.;
    for (Int_t j = low_id; j < high_id; ++j) {
         sum = sum + array->at(j);
     	  }
    return sum * FreqStep;




}

void NMRana::AreaSetLimits(Double_t low_x, Double_t high_x){
		LowArea_X = low_x;
		HighArea_X = high_x;
}

void NMRana::GetTimeStamp(){
	// creates time stamp
	   TimeStamp = timel/10000 -time_offset; //
	   TimeStamp_usec = timel-time_offset*10000;
	   // now split the time in sec and musec
	   long nsec = Int_t(TimeStamp_usec % 10000)*100000;// cpnverts 100 musec into nsek
	   // fell the timespect tructure
	   root_time.tv_sec = TimeStamp;
	   root_time.tv_nsec = nsec;


	   RootTimeStamp = new TTimeStamp(root_time);
}

void NMRana::PrintTime(){
	// prints time from time stamp
	   TimeStamp = timel/10000 -time_offset; //
	   TimeStamp_usec = timel-time_offset*10000;
	   // now split the time in sec and musec
	   long nsec = Int_t(TimeStamp_usec % 10000)*100000;// cpnverts 100 musec into nsek
	   // fell the timespect tructure
	   root_time.tv_sec = TimeStamp;
	   root_time.tv_nsec = nsec;
	   RootTimeStamp=new TTimeStamp(root_time);


//	   cout<<timel<<"   "<<nsec<<"   "<<TimeStamp_usec<<"   \n";
	   time_t timm = Int_t(TimeStamp);

	   tm *ltm = localtime(&timm);
	    cout<<" \n \n ******************************************\n\n";
	    cout << "Year: "<< 1900 + ltm->tm_year << endl;
	       cout << "Month: "<< 1 + ltm->tm_mon<< endl;
	       cout << "Day: "<<  ltm->tm_mday << endl;
	       cout << "Time: "<< 1 + ltm->tm_hour << ":";
	       cout << 1 + ltm->tm_min << ":";
	       cout << 1 + ltm->tm_sec << endl;

	      cout<<asctime(ltm)<< " \n";
	      cout<<" \n \n ******************************************\n\n";

/*       cout << "Time: "<< 1 + ltm->tm_hour << ":";
       cout << 1 + ltm->tm_min << ":";
       cout << 1 + ltm->tm_sec << endl;
*/
}

TH1D *NMRana::SetupStripChart(TString Title){
	// this routine sets up a strip chart for time vs value
	Float_t bintime =1.;
	TH1D *ht = new TH1D("ht",Title,10,0,10*bintime);
	ht->SetStats(0);
	ht->SetLineColor(2);
	ht->GetXaxis()->SetTimeDisplay(1);
	ht->GetYaxis()->SetNdivisions(520);

	return ht;



}



#endif // #ifdef NMRana_cxx













