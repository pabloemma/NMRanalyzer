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
#include <TString.h>


// Header file for the classes stored in the TTree if any.
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <map>



#include "Ana.h"
#include "TEana.h"

using namespace std;

// Fixed size dimensions of array or collections stored in the TTree if any.

class NMRana : public Ana {  //inherit from analysis

private:
	Int_t Control;  			// determines the time reading (when has the file been created and what readout mechanism are we using

protected:    // for TEana inheritance
	Double_t LowArea_X; // lower bound for area calculation
	Double_t HighArea_X;// upper bound for area claculation
	Int_t low_id;		// determines lower and upper index of array for integration
	Int_t high_id;


public :
	std::string NMR_pr ="NMR_ana> ";
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
   Double_t        Phase_Voltage;
   Double_t        Peak_Area;
   Double_t        Pol_Calib_Const;
   Double_t        Gain;
   Double_t        Pol_Sign;
   Double_t        Log_Channel;
   Double_t 	   *gr_freq;
   Double_t		   *gr_amp; // needed for creating and filling the background graph

   Double_t			Qamp; // amplifier setting for Qcurve
   Long64_t         timel; // note this time is down to 100musec, in order to deal only on the second level, strip the last 4 digits
   Long64_t	        time_offset;
   std::vector<double>  *array;
   std::vector<double>  Qcurve_array; // this will contain the QCurve normalized to the sweep number

   Double_t MinFreq ; // limits of frequency sweep
   Double_t MaxFreq ;
   Double_t low_x;  //area limits
   Double_t high_x;
   Double_t SignalArea ; // the summed area of every signal // not normalized yet
   Double_t SignalAreaNormalized ; // aka polarization

   Double_t fit_x1, fit_x2, fit_x3, fit_x4;// limits for fitting widows
   Double_t CalibConstant; // is calculated from area of TE peak and pressure of measuring CalibConstant = pol_calculated/area



   Int_t	DEBUG;// level of debugging information  // currently 0: nodebug
   	   	   	   	   	   	   	   	   	   	   	   	     //1: regular debug
 	   	   	   	   	     //2:
 	   	   	   	   	     //3: You really do not want this




   Bool_t TEmeasurement; // if set to true we have TEmeasurement, gets determined from filename


   std::vector<TString> RootFileArray ; // if there is a list of input files it will put them into vector
   time_t TimeStamp;  // timestamp in seconds on UNIX time 0
   Long64_t TimeStamp_usec;  // timestamp in 100 musec
   timespec root_time; // timespec for root time stamp

   Int_t TimeControl;// this is set to a value depending on the time of the file openin
   	   	   	     // it is used to change the data format according to the evolution of the NMR files
   	   	        // anything after June15 is Control=1, before it is control=0

   	   char *timel_ptr; // because  Root stored the string as a charcater array
   Bool_t QC; // set true for qcurve subtraction
   Bool_t QC_DISP; // If we have a Qcurve we can display it with setting this switch to 1
   Bool_t FitLimit;
   TTree *QCUtree;
   std::map<std::string,std::string> Parameters; // input parameters






   // List of branches
   TBranch        *b_IntScanPoints;   //!
   TBranch        *b_FreqStep;   //!
   TBranch        *b_FreqCenter;   //!

   TBranch        *b_Temperature;   //!
   TBranch        *b_ScanPoints;   //!
   TBranch        *b_TuneV;   //!
   TBranch        *b_Offset;   //!
   TBranch        *b_ControllerV;   //!
   TBranch        *b_Phase_Voltage;   //!
   TBranch        *b_Peak_Area;   //!
   TBranch        *b_Pol_Calib_Const;   //!
   TBranch        *b_Gain;   //!
   TBranch        *b_Pol_Sign;   //!
   TBranch        *b_Log_Channel;   //!
   TBranch        *b_timel;   //!
   TBranch        *b_array;   //!
   TFile		  *f;

   TTree 		   *tree;
   // the block for the histograms
   TH1D 	 *NMR1; // Signal histogram
   TH1D		 *NMR1_NoQ ; // Signal without Qcurve subtraction
   TH1D      *NMR_RT;// real time display of the NMR signal
   TH1D      *NMR_RT_Corr;// real time display of the NMR signal, with background fit subtracted
   TH1D		 *PolTime; // polarization vs time
   TH1D		 *CalibTime; // TE calibration constant vs time
   TH1D		 *Qcurve_histo; // Displays Qcurve if it will be subtracted
   TGraph    *Background; // this is the background I determine from the signal thorugh a spline
   //
   TCanvas	 *GeneralCanvas;  // has the signal and polarization vs time on it
   TCanvas	 *StripCanvas;   // shows polarization vs time
   TCanvas	 *StripCanvas_1; // for TE measurement shows the pressure over time.
   TCanvas	 *AuxCanvas;   // all the auxiliary plots, like Qcurve
   TCanvas	 *RTCanvas ;    // Life time display canvas, get used in Loop

   TChain    *NMRchain; // if we have more than one root file
   TString timestring; // from labview time
   TDatime *td ; // Datetime for time histogram
   TTimeStamp *RootTimeStamp;
   TFile * QcurveFile ; // Qcurve file name


	TEana TE;




   NMRana();
   virtual ~NMRana();
   virtual Int_t    Cut(Long64_t entry);
   virtual void     Show(Long64_t entry = -1);
   virtual int      OpenChain(std::vector<TString> );
   virtual void		AreaSetLimits(Double_t , Double_t);
   // virtual Double_t CalculateArea(std::vector<Double_t> *);
    virtual Double_t CalculateArea(TH1D *);
   virtual void		PrintTime();  // prints time from Labview time stamp
   virtual void 	GetTimeStamp();
   virtual void     GetQcurve(std::string );
   virtual void		FillQcurveArray();
   virtual void	    SetTimeControl(int);

   // memebre which will be inhertied from TEana
   void     Loop(); // TEana inherits
   void		SetupCanvas();
   TH1D * 	SetupStripChart(TString);
   void	    SetupHistos();
   void		ReadParameterFile(TString );
   void     Init(TTree *tree);
   Bool_t   Notify();
   Int_t    GetEntry(Long64_t entry);
   Long64_t LoadTree(Long64_t entry);
   Int_t 	OpenFile(TString);
   void		CloseFile();
   void		DrawHistos();
   TString  GetDate(TString input);


};
#endif

#ifdef NMRana_cxx


NMRana::NMRana(){
	time_offset = 2082852019 ; // unix offset in seconds
	TimeControl = 1; // always assume the file is from the newest generation // if not use TIMEC = value in parameter file
	QC_DISP = false;



}
void NMRana::ReadParameterFile(TString ParameterFile){
	// reads in parameters for running the NMRanalyzer
	// needs the Qcurve file
	// calibration constants from TE measurements
	// line length is a maximum of 80

	// the format of the file is  name , value
	// example: Qcurve   Qcurve.root
	char temp_string[132];
	std::string temp;
	std::string string1, string2;
	std::string temp_file;

	Double_t lownmr(212.8);
	Double_t hinmr(213.3);

	ifstream ParFile; // create instream file
	ParFile.open(ParameterFile);
	// check if found and opened
	if(!ParFile.is_open()){
		exit(EXIT_FAILURE);
	}
	// read the lines


	do{
		ParFile >>string1 >> string2;
		if (ParFile.eof()) break;  // get out of the loop
		if(string1.find("#") == std::string::npos) Parameters[string1] = string2;  // check for comments

	}while(ParFile.good());


	// Now print out parameter map
	for( std::map<string,string>::iterator pos=Parameters.begin(); pos !=  Parameters.end(); ++pos){
		cout<<NMR_pr<<"parameters from file  :"<<pos->first<<"\t"<<pos->second <<"\n";


		if(pos->first.find("QCurve")!= std::string::npos){
			QC=true;
			temp_file = pos->second;
		}


		if(pos->first.find("QAMP")!= std::string::npos){
			// amplifier setting for QCurve
			Qamp = std::stod(pos->second);
			cout<<NMR_pr<<" test qamp"<<string2<<"\n";

		}
		if(pos->first.find("TIMEC")!= std::string::npos){
			// amplifier setting for QCurve
			TimeControl = std::stoi(pos->second);

		}
		if(pos->first.find("NMR_LOW")!= std::string::npos){
			// amplifier setting for QCurve
			lownmr = std::stod(pos->second);

		}
		if(pos->first.find("NMR_HI")!= std::string::npos){
			// amplifier setting for QCurve
			hinmr = std::stod(pos->second);

		}
		if(pos->first.find("QC_DISPLAY")!= std::string::npos){
			// amplifier setting for QCurve
			if(std::stoi(pos->second)==1) QC_DISP=true ;

		}
		if(pos->first.find("FITX1")!= std::string::npos){
			// amplifier setting for QCurve
		fit_x1 = std::stod(pos->second);
			FitLimit=true;
		}
		if(pos->first.find("FITX2")!= std::string::npos){
			// amplifier setting for QCurve
		 fit_x2 = std::stod(pos->second);

		}
		if(pos->first.find("FITX3")!= std::string::npos){
			// amplifier setting for QCurve
		fit_x3 = std::stod(pos->second);

		}
		if(pos->first.find("FITX4")!= std::string::npos){
			// amplifier setting for QCurve
		fit_x4 = std::stod(pos->second);

		}
		if(pos->first.find("DEBUG")!= std::string::npos){
			// amplifier setting for QCurve
		DEBUG = std::stoi(pos->second);

		}


	}
	if(QC) GetQcurve(temp_file);


		// set limits for are either default or what comes in from parameter file
	if(lownmr > hinmr){
		cout<<NMR_pr<<"the integration limits are wrong Low limit is larger than high limit  "<<lownmr<<"  "<<hinmr<<"\n";
	}
	else  AreaSetLimits(lownmr,hinmr);
	if(FitLimit)SetFitLimits(fit_x1,fit_x2,fit_x3,fit_x4);


}







int NMRana::OpenFile(TString rootfile){

	// oepn file and initialize tree
     cout<<NMR_pr<<"opening file "<<rootfile<<"\n";

     if(rootfile.Contains("TER")){
    	 TEmeasurement = true;
    	 cout <<NMR_pr<< "\n\n this is a TE measurement \n\n\n";
     }

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
			cout<<NMR_pr<<RootFileArray[pos]<<"   filename \n";
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
   // add the new branches for the later data format
   if(TimeControl == 1){
	   fChain->SetBranchAddress("Phase_Voltage", &Phase_Voltage, &b_Phase_Voltage);
	   fChain->SetBranchAddress("Peak_Area", &Peak_Area, &b_Peak_Area);
	   fChain->SetBranchAddress("Pol_Calib_Const", &Pol_Calib_Const, &b_Pol_Calib_Const);
	   fChain->SetBranchAddress("Gain", &Gain, &b_Gain);
	   fChain->SetBranchAddress("Pol_Sign", &Pol_Sign, &b_Pol_Sign);
	   fChain->SetBranchAddress("Log_Channel", &Log_Channel, &b_Log_Channel);
	    }


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
	   NMR_RT = new TH1D("NMR_RT","Real TimeSignal histogram",IntScanPoints,MinFreq,MaxFreq);
	   NMR_RT->SetLineColor(kSpring-2);
	   NMR_RT->SetLineWidth(4);
	   NMR_RT_Corr = new TH1D("NMR_RT_corr","Real TimeSignal background subtracted",IntScanPoints,MinFreq,MaxFreq);
	   NMR_RT_Corr->SetLineColor(kBlue-2);
	   NMR_RT_Corr->SetLineWidth(4);

	   // Determine the Integration or summation limits for peak in terms of channels.
	   low_id = NMR1->GetXaxis()->FindBin(LowArea_X);
	   high_id = NMR1->GetXaxis()->FindBin(HighArea_X);

       //




	   // Now setup the time versus polarization histo
	   // the time is in 100 microseconds, which is much more precise than what we need
	   // so I didvide by 10000  with an integer division, also the time and date is now in UNIX time
	   //####################Note the time is only correct if the updates take about 1 second.
	   // if it is longer or shorter it is not correct, since tickmaks assue 1 second.######
	   GetTimeStamp();
	   td = new TDatime(TimeStamp);
	   td->Print();
	   // we set the time stamp to 0
	   //gStyle->SetTimeOffset(td->Convert());
	   PolTime = SetupStripChart("Polarization vs time");
	   PolTime->SetMaximum(100.);
	   PolTime->SetMinimum(-100.);
	   PolTime->SetLineColor(2);

	   //setup calibration constant histogram
	   if(TEmeasurement){
		   CalibTime = SetupStripChart("Calibration Constant vs time");
		   CalibTime->SetMaximum(100.);
		   CalibTime->SetMaximum(0.);

	   }




	   // Print original time
	   PrintTime();


	   Int_t timewindow = 1;
	   PolTime = new TH1D("PolTime","Polarization vs time",timewindow,0,10*timewindow);
	   PolTime->GetXaxis()->SetTimeDisplay(1);
	   PolTime->GetXaxis()->SetTimeOffset(TimeStamp);
	   PolTime->GetXaxis()->SetTimeFormat("%d\ %m\ %H\:%M \:%S");
	   PolTime->GetXaxis()->SetNdivisions(405) ;
	   // Now setup a Canvas for Qcurve

	   if(Qcurve_array.size()!=0){
		   Qcurve_histo = new TH1D("Qcurve_hist","Normalized Qcurve histogram",IntScanPoints,MinFreq,MaxFreq);
		   NMR1_NoQ = new TH1D("NMR1_NoQ","Signal without QCurve subtraction",IntScanPoints,MinFreq,MaxFreq);
		   }

	   // Do the graph for the histo
	   gr_freq = new Double_t[IntScanPoints];
	   gr_amp = new Double_t[IntScanPoints];
	  // Background = new TGraph(IntScanPoints,gr_freq,gr_amp);


}
void NMRana::SetupCanvas(){
	// creates the different Canvas
	// master canvas for all histograms
	GeneralCanvas = new TCanvas("GeneralCanvas","NMR signal",200,50,800,800);


	//canvas for all strip charts
	StripCanvas =  new TCanvas("StripCanvas","NMR strip charts",1020,50,1000,600);
	StripCanvas->SetGrid();
	StripCanvas->SetFillColor(42);
	StripCanvas->SetFrameFillColor(33);
    // only do a strip chart for pressure if we have a TE measurement.
	if(TEmeasurement){
		StripCanvas_1 =  new TCanvas("StripCanvas_1","Calibration Constant strip charts",1020,550,1000,600);
		StripCanvas_1->SetGrid();
		StripCanvas_1->SetFillColor(40);
		StripCanvas_1->SetFrameFillColor(30);

	}

	RTCanvas =  new TCanvas("RTCanvas","Real Time charts",100,1000,600,400);
	RTCanvas->SetGrid();
	RTCanvas->SetFillColor(23);
	RTCanvas->SetFrameFillColor(16);
	RTCanvas->Divide(1,2);





	if(Qcurve_array.size()!=0){

	AuxCanvas = new TCanvas("AuxCanvas","Auxiliary plots",1020,700,1000,600);

	if(QC_DISP)AuxCanvas->Divide(1,2); // if we want to see the Qcurve as well
	else AuxCanvas->Divide(1,1);

	}
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

	if(Qcurve_array.size()!=0){
		if(QC_DISP){
			AuxCanvas->cd(2);
			Qcurve_histo->Draw();
			AuxCanvas->cd(1);
			NMR1_NoQ->Draw();

		AuxCanvas->Update();
		}
		else {
			AuxCanvas->cd(1);
			NMR1_NoQ->Draw();
			AuxCanvas->Update();
		}

	}


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
   if(TEmeasurement)StripCanvas_1->cd();

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t time_prev = 0;
   Long64_t nbytes = 0, nb = 0;

   RTCanvas->cd(1);
   NMR_RT->Draw();
   RTCanvas->cd(2);
   NMR_RT_Corr->Draw();
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
	   NMR_RT->Reset();
	   NMR_RT_Corr->Reset();
	   Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

//	now fill histogram
// reset freq_temp to lower bound
      Double_t freq_temp = MinFreq;


      for (UInt_t j = 0; j < array->size(); ++j) {
    	  // subtract QCurve if existing

    	  if(Qcurve_array.size()!=0) {
              NMR1_NoQ->Fill(freq_temp,array->at(j));
             array->at(j) = array->at(j) - Qcurve_array.at(j);
    	  }

          NMR1->Fill(freq_temp,array->at(j));
          NMR_RT->Fill(freq_temp,array->at(j));
          NMR_RT_Corr->Fill(freq_temp,array->at(j));
          // for backgroiund graph
          gr_freq[j] = freq_temp;
          gr_amp[j] = array->at(j);
          freq_temp = freq_temp+FreqStep;
      	  }
//		fill the background graph and go to determine the spline
//      Background = new TGraph(IntScanPoints,gr_freq,gr_amp);
//      BackSpline(Background);
  	  //FindPeak(NMR1);
      FitBackground(NMR_RT_Corr);
//sum the peak area
      StripCanvas->cd();
      SignalArea = CalculateArea(NMR_RT_Corr);
// now for every point in a TE we will calculate the polarization from the pressure
// the ratio of calculated polarization/ area gives the calibration constant calib
// so that calib*area = polarization of the reL SIGNAL
// at the end we will calculate an average caibration constant with a deviation





 	  GetTimeStamp();
	  PolTime->SetBinContent(jentry,SignalArea);
	  PolTime->GetXaxis()->SetRange(jentry-10000,jentry+5);
	  StripCanvas->Clear();
	  PolTime->Draw();
	  StripCanvas->Modified();
	  StripCanvas->Update();

	  if(TEmeasurement){
		  if(jentry ==0)cout<<NMR_pr<<"!!!!!!!!!!!!!!!need to change the call to TE polarization calculation!!!!!!!!!!!\n";
		  StripCanvas_1->cd();
		  CalibConstant = TE.CalculateTEP("proton",.5,5.,.117) ; // needs to change
		  CalibConstant = CalibConstant/SignalArea;
		  CalibTime->SetBinContent(jentry,CalibConstant);
		  CalibTime ->GetXaxis()->SetRange(jentry-10000,jentry+5);
		  StripCanvas_1->Clear();
		  CalibTime->Draw();
		  StripCanvas_1->Modified();
		  StripCanvas_1->Update();

	  }

// draw the signal histogram
	  RTCanvas->cd(1);
	  NMR_RT->Draw();
	  RTCanvas->cd(2);
	  NMR_RT_Corr->Draw();
	  RTCanvas->Modified();
	  RTCanvas->Update();

	  if(DEBUG ==2)cout<<NMR_pr<<timel<<"another one \n";

   }// end of entry loop
   // Now fill the QCurve histo if it is used
	  if(Qcurve_array.size()!=0){
	      Double_t freq_temp = MinFreq;


	      for (UInt_t j = 0; j < Qcurve_array.size(); ++j) {
	    	  // subtract QCurve if existing

	          Qcurve_histo->Fill(freq_temp,Qcurve_array.at(j));
	          freq_temp = freq_temp+FreqStep;
	      	  }

	  }
	  //Fill qcurve histogram if desired
/*	  if(QC_DISP){
	      Double_t freq_temp = MinFreq;


	      for (UInt_t j = 0; j < Qcurve_array.size(); ++j) {
	    	  // subtract QCurve if existing

	          QCC->Fill(freq_temp,Qcurve_array.at(j));
	          freq_temp = freq_temp+FreqStep;
	      	  }

	  }
*/
	  	  	 delete [] gr_freq;
	  	  	 delete [] gr_amp;
}

Double_t NMRana::CalculateArea(TH1D *histo){
	// this function calculates the area of the NMR perak by simpy summing it in
	//  limits given by AreaSetLimits in the main program
	// the area is calculated as sum i_low to i_high (a(i)*binwidth)

	Double_t sum1 = 0.;
    for (Int_t j = low_id; j < high_id; ++j) {
         sum1 += array->at(j);
     	  }

		//Double_t sum = histo->Integral(histo->FindBin(fit_x2),histo->FindBin(fit_x3));
    Double_t sum11 =0;
    Double_t sum = histo->Integral(histo->FindBin(low_id),histo->FindBin(high_id));
		for(Int_t k=low_id;k<high_id;k++){
		if(DEBUG ==3)cout<<NMR_pr<<histo->GetBinContent(k)<<"  \n";
		sum11 += histo->GetBinContent(k);
		}
		if(DEBUG==3){
			cout<<NMR_pr<<sum1<<"  "<<sum<<"  "<<sum11<<"\n";

		   cout<<NMR_pr<<low_id<<"   "<<high_id<<" \n";
		}
    	    return sum ;
   // return sum * FreqStep;




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


//	   cout<<NMR_pr<<timel<<"   "<<nsec<<"   "<<TimeStamp_usec<<"   \n";
	   time_t timm = Int_t(TimeStamp);

	   tm *ltm = localtime(&timm);
	    cout<<NMR_pr<<" \n \n ******************************************\n\n";
	    cout<<NMR_pr << "Year: "<< 1900 + ltm->tm_year << endl;
	       cout<<NMR_pr << "Month: "<< 1 + ltm->tm_mon<< endl;
	       cout<<NMR_pr << "Day: "<<  ltm->tm_mday << endl;
	       cout<<NMR_pr << "Time: "<< 1 + ltm->tm_hour << ":";
	       cout<<NMR_pr << 1 + ltm->tm_min << ":";
	       cout<<NMR_pr << 1 + ltm->tm_sec << endl;

	      cout<<NMR_pr<<asctime(ltm)<< " \n";
	      cout<<NMR_pr<<" \n \n ******************************************\n\n";

/*       cout<<NMR_pr << "Time: "<< 1 + ltm->tm_hour << ":";
       cout<<NMR_pr << 1 + ltm->tm_min << ":";
       cout<<NMR_pr << 1 + ltm->tm_sec << endl;
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

void NMRana::GetQcurve(std::string filename) {

	cout<<NMR_pr<<" QCurve file "<<filename<<"\n";

	QcurveFile = new TFile(filename.c_str());
	// now create a new tree, so that we have a different name
	QCUtree = (TTree*)QcurveFile->Get("NMRtree");
	Init(QCUtree);
	FillQcurveArray();
}

void NMRana::FillQcurveArray(){
	// this function fill the Qcurve array and normalizes it
	   Long64_t nentries = QCUtree->GetEntriesFast();
	   Long64_t nbytes = 0, nb = 0;
       if(DEBUG ==1) cout<<NMR_pr<<"FillQcurve>  "<<" nentries"<< nentries<<" \n";
       Qcurve_array.clear();
	   for (Long64_t jentry=0; jentry<nentries;jentry++) {
	      Long64_t ientry = LoadTree(jentry);
	      if (ientry < 0) break;
	      nb = QCUtree->GetEntry(jentry);   nbytes += nb;
	      Double_t freq_temp = MinFreq;


	      for (UInt_t j = 0; j < array->size(); ++j) {
	          //NMR1->Fill(freq_temp,array->at(j));
	    	  // I have to decide if I want to subtract ROOT histos or arrays.
	    	  // I think arrays is better.
	    	  // so create new array vetcor and continuosly add and ifinally mormalize it to number of
	    	  //sweps.
	          freq_temp = freq_temp+FreqStep;
	          if(jentry ==0){  // the first time we create the array
        	  Qcurve_array.push_back(array->at(j));

	          }
	          else {
	        	  Qcurve_array.at(j) = Qcurve_array.at(j)+array->at(j);
	          }
	      	  }



	   }// end of entry loop

	   //now we need to normalize the QCurve to the number of sweeps, which is nentries
	   if(DEBUG==1)cout<<NMR_pr<<"FillQcurve>  "<<" qamp"<<Qamp<<"\n";
	   for (UInt_t j = 0; j < array->size(); ++j) {

		   Qcurve_array.at(j)= Qcurve_array.at(j)/nentries/Qamp;

	   }



}
void NMRana::SetTimeControl(Int_t main_input){
	TimeControl = main_input;
}

TString NMRana::GetDate(TString input) {


	TString timestring = input;
    time_t time_test = timestring.Atoi()-2082844800; // calculated with offset since stupid labview uses jan-1-1904
    tm *ltm = localtime(&time_test);          //and unix uses jan-1-1970

    cout<<NMR_pr<<" \n \n ******************************************\n\n";
    cout<<NMR_pr << "Year: "<< 1900 + ltm->tm_year << endl;
       cout<<NMR_pr << "Month: "<< 1 + ltm->tm_mon<< endl;
       cout<<NMR_pr << "Day: "<<  ltm->tm_mday << endl;
       cout<<NMR_pr << "Time: "<< 1 + ltm->tm_hour << ":";
       cout<<NMR_pr << 1 + ltm->tm_min << ":";
       cout<<NMR_pr << 1 + ltm->tm_sec << endl;

      cout<<NMR_pr<<asctime(ltm)<<"  "<<time_test<<"   "<<"    \n";
      cout<<NMR_pr<<" \n \n ******************************************\n\n";
      if(Int_t(time_test) > 1465948800) TimeControl =1; // this gives a control value for which time the polarization file is from
      return  asctime(ltm);
}

#endif // #ifdef NMRana_cxx

