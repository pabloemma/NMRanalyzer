/*
 * TEana.h
 *
 *  Created on: Jun 7, 2016
 *      Author: klein
 *      * came from  Makeclass
 *      * will calculate the calibration constant
 */

#ifndef TEana_h
#define TEana_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TDatime.h>
#include <TStyle.h>  // so we can use gStyle
#include <TTimeStamp.h>
#include <TMath.h>
#include <TGraph.h>


// Header file for the classes stored in the TTree if any.
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>





// Fixed size dimensions of array or collections stored in the TTree if any.

class TEana {
private:
	Double_t deuteron_g;
	Double_t nucleon_mag_moment;
	Double_t proton_g;
	Double_t neutron_g;
	Double_t boltzmann ; // J/K
	Double_t proton_fact;
	Double_t deuteron_fact;
	Int_t TimeControl;


	// coefficiewnts for helium pressute to T calculation
	Double_t aLowT[10];
	Double_t bLow;
	Double_t cLow;
	Double_t aHighT[10];
	Double_t bHigh;
	Double_t cHigh;
	TGraph *lowT;  // used for temp calulation at low pressure


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
   Double_t        Phase_Voltage;
   Double_t        Peak_Area;
   Double_t        Pol_Calib_Const;
   Double_t        Gain;
   Double_t        Pol_Sign;
   Double_t        Log_Channel;
   vector<double>  *array;
   Long64_t        timel;
   Int_t           IntScanPoints;
   TFile		   *f;
   TTree           *tree;
   TChain    	   *NMRchain; // if we have more than one root file




   // List of branches
   TBranch        *b_FreqCenter;   //!
   TBranch        *b_FreqStep;   //!
   TBranch        *b_Temperature;   //!
   TBranch        *b_ScanPoints;   //!
   TBranch        *b_TuneV;   //!
   TBranch        *b_Offset;   //!
   TBranch        *b_ControllerV;   //!
   TBranch        *b_Peak_Area;   //!
   TBranch		  *b_Phase_Voltage;
   TBranch        *b_Pol_Calib_Const;   //!
   TBranch        *b_Gain;   //!
   TBranch        *b_Pol_Sign;   //!
   TBranch        *b_Log_Channel;   //!
   TBranch        *b_array;   //!
   TBranch        *b_timel;   //!
   TBranch        *b_IntScanPoints;   //!

   TEana();
   virtual ~TEana();
 //  virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
//   virtual void     Show(Long64_t entry = -1);
   virtual Int_t 	OpenFile(TString);
   virtual int      OpenChain(std::vector<TString> );
   virtual void CloseFile();
   virtual Double_t  CalculateTEP(std::string, Double_t,Double_t , Double_t );
   virtual Double_t	CalcT(Double_t); // calculates temperature from pressure, input in TORR
   virtual Double_t CalculateT(Double_t *, Double_t , Double_t, Double_t);

};









#endif /* TEana_h */


#ifdef TEana_cxx



TEana::TEana(){
		cout<<"**************       Initialize TE analyzer **************** \n\n\n";
		proton_g    =  5.585694702;
		neutron_g   = -3.82608545;
		deuteron_g = 0.85741; // needs to be multiplied with nucleon mag moment (kind of a g factor)
		nucleon_mag_moment  = 5.05078353e-27;  //SI J/T
		boltzmann = 1.38064852e-23 ; // SI units J/K
		proton_fact = nucleon_mag_moment * proton_g /boltzmann;
		deuteron_fact =  nucleon_mag_moment * deuteron_g /boltzmann;

		// setup the coefficients for P to T calculation

		aLowT[0] =1.392408;
		aLowT[1] =0.527153;
		aLowT[2] =0.166756;
		aLowT[3] =0.050988;
		aLowT[4] =0.026514;
		aLowT[5] =0.001975;
		aLowT[6] =-.017976;
		aLowT[7] =0.005409;
		aLowT[8] =0.013259;
		aLowT[9] =0.0;
		bLow = 5.6;
		cLow = 2.9;



		aHighT[0] =3.146631;
		aHighT[1] =1.357655;
		aHighT[2] =0.413923;
		aHighT[3] =0.091159;
		aHighT[4] =0.016349;
		aHighT[5] =0.001826;
		aHighT[6] =-.004325;
		aHighT[7] =-.004973;
		aHighT[8] =0.0;
		aHighT[9] =0.0;
		bHigh = 10.3;
		cHigh = 1.9;

		// setup TGparh for our later interpolation
		Double_t tt[13] = {.650,.7,.75,.8,.85,.9,.95,1.0,1.05,1.1,1.15,1.2,1.25};
		Double_t pp[13] = {1.101e-01,2.923e-01,6.893e-01,1.475,2.914,5.380,9.381,15.58,24.79,38.02,56.47,81.52,114.7};
		lowT = new TGraph(13,pp,tt);

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

int TEana::OpenChain(std::vector<TString> RootFileArray){

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

void TEana::CloseFile(){
	f->Close();
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
   if(TimeControl == 1){
	   fChain->SetBranchAddress("Phase_Voltage", &Phase_Voltage, &b_Phase_Voltage);
	   fChain->SetBranchAddress("Peak_Area", &Peak_Area, &b_Peak_Area);
	   fChain->SetBranchAddress("Pol_Calib_Const", &Pol_Calib_Const, &b_Pol_Calib_Const);
	   fChain->SetBranchAddress("Gain", &Gain, &b_Gain);
	   fChain->SetBranchAddress("Pol_Sign", &Pol_Sign, &b_Pol_Sign);
	   fChain->SetBranchAddress("Log_Channel", &Log_Channel, &b_Log_Channel);
	    }
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
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();


   	  Double_t press = .05;
   	  Double_t p_step=.01;
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;


      cout<<"pressure  "<<press<<"   temp:  "<<CalcT(press)<<"   Polarization:  "<<CalculateTEP("proton",.5,5.,press)<<" \n";
      press = press+p_step;
      // if (Cut(ientry) < 0) continue;
   }
}

Double_t TEana::CalculateTEP(std::string particle ,Double_t spin, Double_t field, Double_t pressure){
	// this routine calculates the TE polarization for a pressure P and a give field
	// It calculates the Brouillouiin  function for p and d currently
	Double_t TEPol;

	//
	Double_t fact; // calculates the constants
	Double_t brill; //Brillouin function
	Double_t Temp = CalcT(pressure); // get temperature from calculation
	Double_t arg1 = (2.*spin+1.)/(2.*spin);
	Double_t arg2 = (1.)/(2.*spin);
	if(particle.find("proton") != std::string::npos){
		fact = proton_fact *spin *field /Temp;
	}
	else if(particle.find("deuteron") != std::string::npos){
		fact = deuteron_fact *spin *field /Temp;
	}
	else {
		cout<<"TEana> No particle found for TE calculation************\n\n\n ";
	}

	TEPol = arg1*cosh(arg1*fact)/sinh(arg1*fact)-arg2*cosh(arg2*fact)/sinh(arg2*fact);



	return TEPol;
}

Double_t TEana::CalcT(Double_t pressure){
	// Calculates Temp as a function of pressure in Torr
	// follows Journal of Physical and Chemical Ref. Data 27, 1217 (1998)
	//
	Double_t temp;
	Double_t pa = pressure*133.322; // convert pressure into pascal
	if(pressure>.0009 && pressure<.826)
		temp = lowT->Eval(pa,0,"S");  //cubic spline

	else if(pressure>=.826 && pressure<37.82)
		temp = CalculateT(aLowT,bLow,cLow,pa);

	else if(pressure>=72.82 && pressure<1471.)
		temp = CalculateT(aHighT,bHigh,cHigh,pa);


	return temp;
}


Double_t TEana::CalculateT(Double_t *a, Double_t b, Double_t c, Double_t press){
	// calculate temperature
	Double_t con = (TMath::Log(press)-b)/c;
	Double_t t_i = 0.0;
	for(Int_t k=0 ;k<10;k++){
		t_i = t_i+a[k]*pow(con,k);
	}
	return t_i;

}


#endif // #ifdef TEana_cxx


