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
/*
 * histos:
 * NMR1: Signal Histogram (all signals in a final histogram) withb Qcurve subtrtacted
 * NMR1_Qfit: Signal Histogram with Qfit subtracted
 * NMR1_B ; // created by Ana::Fitbackground and is the background fitted and subtracted NMR1 spectrum
 * NMR1_Qfit_B ; // created by Ana::Fitbackground and is the background fitted and subtracted NMR1_Qfit spectrum
 *
 * NMR_RT: Real_time signal histo ; the 20 sweep histograms // reset everye sweep
 * NMR_RT_Corr; NMR_RT - Qcurve histogram reset every sweep;  from this we calculate the polarization
 * * however careful; I am doing a linear background subtraction det
 * *
 * Qcurve_histo: Normazlied Qcurve histogram
 * Qcurve_histo_shifted: shifted QCurve hist  by number of channels, so that the min of the data and the qcurve agree
 *                        only used at the moment when we do the aux plots
 * NMR1_NoQ: NMR signal without Qcurve subtraction
 *
 * raw_histo: raw signal; only used in Debug mode
 * raw_histo_QC: signal histo with QC subtracted.
 *
 * Stripcharts:
 * PolTime: polarization vs time (for regular measurement)
 * PolTime: TE signal area vs time for TE measurement
 * CalibTime: Calibration constant vs time; gets calculated from polarziation and temp
 * PressTime: Pressure vs time (or convesrly temp)
 * SysTempTime: System temeperature vs time (from NMR system boards)
 * PatTemp: This is the difference of the backgrounds left - right and an indiaction of temeperature drift.
 *
 */

#ifndef NMRana_h
#define NMRana_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TDatime.h>
#include <TStyle.h>  // so we can use gStyle
#include <TTimeStamp.h>
#include <TString.h>


// Header file for the classes stored in the TTree if any.
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdio.h>
#include <map>



#include "AnaGraph1.h"
#include "TEana.h"
#include "NMRFastAna.h"  //Kun's Qcurve hifting class
#include "NMR_DST.h"
using namespace std;

// Fixed size dimensions of array or collections stored in the TTree if any.

class NMRana : public Ana {  //inherit from analysis

private:
	Int_t Control;  			// determines the time reading (when has the file been created and what readout mechanism are we using
    Int_t NumberOfStripCharts;
    Double_t PatLowArea;
    Double_t PatHighArea;   // the lower and higher areas in thye NMR signal
    Double_t PatDiffArea;  		// the difference of PatLow and PatHighArea
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
   Double_t			ScanNumber;
   Double_t        Pol_Sign;
   Double_t        Log_Channel;
   Double_t 		Peak_Amp;//
   Double_t 		NMRchan; // whci coil
   Double_t 		PeakCenter ;//
   Double_t 		BeamOn;//
   Double_t 		RFlevel;//
   Double_t 		IFatt; //

   Double_t 		HeT;
   Double_t			HeP;



   Double_t			Qamp; // amplifier setting for Qcurve
   Long64_t         timel; // note this time is down to 100musec, in order to deal only on the second level, strip the last 4 digits
   Long64_t	        time_offset;
   std::vector<double>  *array = 0;
   std::vector<int>  xoffset ; // vector of offsert for Qcurve, filled by NMRfastana
   std::vector<double>  yoffset ; // vector of y offsert for Qcurve, filled by NMRfastana

   Double_t MinFreq ; // limits of frequency sweep for histos
   Double_t MaxFreq ;
   Double_t MinFreqR ; // real limits of frequency sweep
   Double_t MaxFreqR ;
   Double_t low_x;  //area limits
   Double_t high_x;

   Double_t low_fit_x;  // center of lower fit window; fit_x2+fit_x1/2
   Double_t high_fit_x; //width of higher fit window; used by subtract linear to remove any slope


   Double_t SignalArea ; // the summed area of every signal // not normalized yet
   Double_t SignalAreaNormalized ; // aka polarization

   Double_t fit_x1, fit_x2, fit_x3, fit_x4;// limits for fitting widows
   Int_t Ifit_x1, Ifit_x2, Ifit_x3, Ifit_x4;// corresponding channel limits for fitting widows
   Double_t CalibConstant; // is calculated from area of TE peak and pressure of measuring CalibConstant = pol_calculated/area
   Double_t CalConst;// Calibration constant read from parameter file
   Double_t gain_array[3];
   Double_t NumberOfSweeps;  // total number of sweeps
   Double_t TotalEntries; //running tally of all the entries from the different runs, used to renormaize the histograms NMR1, NMR1-NoQ and NMR1_Qfit

   Double_t CurveOffset; // this is the offset which is caluclated as an averag between the left and the right window
   	   	   	   	   	   	   // it is subratced from the signal and then the calibration constant is calculated and the area

   // The variables asociated with the Qcurve fit from QCana
   Double_t QcurveCoil; // the corresponding coil
   Double_t QcurveMin; // the min found in the fit
   Double_t QcurveTune; // The tune volgae of the QCurve. This will be compared with the tune voltage of the data to make sure
   	   	   	   	   	     // they agree. if they don;'t give a warning, and l;ater correct with firts
   Double_t QcurveScale;// This is the number of sweeps in the Qcuvre root file which was used to do the fit. This gives the scale
   	   	   	   	   	   	   // how much we have to renormalize: It is really fitpar/Qcurve_scale
   Double_t QcurveGain; // gain for QCurve measurement
   Double_t QfitPar[5]; // the parameters determined from the QCurve fit and written to the Qcurve.txt; these are normalized by QcurveScale !! so they are really
   	   	   	   	   	   	 // QfitPar/QcurveScale;
   // end of Qcurve parameters

   // Variables associated with the real Qcurve file
   Long64_t QcurveEntries; // the number of entries in the measured Qcurve
   std::vector<double>  Qcurve_array; // this  QCurve array normalized to the sweep number and the gain.
   Double_t QCtune; // used to make sure that Qcurve tune voltages are the same as from data taking
   Double_t QCcoil; // used to make sure that Qcurve tune voltages are the same as from data taking

   const char * analysis_text_file; // if cosen in parameter file, redirects output to analysis text file
   	   Int_t StripLength ; // how many points in the stripchart


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
   Bool_t QCshift; // if true calcuate the QCurve shift and subtract the shifted Curve
   Bool_t TXT ;//output to file
   Bool_t Ave ; // we average over TE sweeps
   Bool_t Reached_Average;// flag we set when we have formed an average
   	   	   	   	   	      // this is needed for caluclating the calib constant in stripper

   TTree *QCUtree;
   TTree *mytr;
   std::map<std::string,std::string> Parameters; // input parameters
   std::vector<Double_t> CalibConstantVector; // this hold the calculated calibration constants and will calculate mean and error
   Double_t CalibConstantMean;
   Double_t CalibConstantMeanError;

   Int_t AverageNumber ; // number of TE sweeps to average over



   std::string  NMR_ROOT ; // top directory of NMR system, needs to be defined thorugh enviro variable $NMR_ROOT
   TString DstFile;
   // define structure for qcurve parameters
   // needed for reading in the csv file.







   // List of branches
   TBranch        *b_IntScanPoints;   //!
   TBranch        *b_FreqStep;   //!
   TBranch        *b_FreqCenter;   //!

   TBranch        *b_Temperature;   //!
   TBranch        *b_ScanPoints;   //!
   TBranch        *b_ScanNumber;   //!
   TBranch        *b_TuneV;   //!
   TBranch        *b_Offset;   //!
   TBranch        *b_ControllerV;   //!
   TBranch        *b_Phase_Voltage;   //!
   TBranch        *b_Peak_Area;   //!
   TBranch        *b_Pol_Calib_Const;   //!
   TBranch        *b_QcurveAmp;   //!

   TBranch        *b_Gain;   //!
   TBranch        *b_Pol_Sign;   //!
   TBranch        *b_Log_Channel;   //!
   TBranch        *b_Peak_Amp;   //!
   TBranch        *b_HeT;   //!
   TBranch        *b_HeP;   //!
	TBranch 	  *b_NMRchan; // whci coil
	TBranch 		*b_PeakCenter ;//
	TBranch 		*b_BeamOn;//
	TBranch 		*b_RFlevel;//
	TBranch 		*b_IFatt; //

   TBranch        *b_timel;   //!
   TBranch        *b_array;   //!
   TFile		  *f;

   TTree 		   *tree;
   // the block for the histograms
   TH1D 	 *NMR1; // Signal histogram
   TH1D		 *NMR1_NoQ ; // Signal without Qcurve subtraction
   TH1D		 *NMR1_Qfit ; // Signal Qfit subtraction
   TH1D 	 *NMR1_B ; // created by Ana::Fitbackground and is the background fitted and subtracted NMR1 spectrum
   TH1D 	 *NMR1_Qfit_B ; // created by Ana::Fitbackground and is the background fitted and subtracted NMR1 spectrum


   TH1D      *NMR_RT;// real time display of the NMR signal
   TH1D      *NMR_RT_Corr;// real time display of the NMR signal, with qcurve subtracted
   TH1D      *NMR_RT_Corr_Fit;// real time display of the NMR signal, with qcurrve subtracted and background fit subtracted
   TH1D		 *PolTime; // polarization vs time
   TH1D		 *CalibTime; // TE calibration constant vs time
   TH1D		 *PressTime ; // Pressure vs time for TE measurement
   TH1D		 *SysTempTime ; // Pressure vs time for TE measurement
   TH1D		 *PatTemp; // pats way to check for temeparture shifts, uses low and high integral of signal and stakes the difference
   TH1D		*CrateTemp; // temperature of crate vs time
   TH1D		 *Qcurve_histo; // Displays Qcurve if it will be subtracted
   TH1D		 *Qcurve_histo_shifted; // Displays Qcurve shifted if it will be subtracted
   TH1D	     *ht[20];// Number of strip charts
   TH1D		 *raw_histo;// this is the histogram filled by the raw numbers//
   TH1D		 *raw_histo_QC;// this is the histogram filled by the raw numbers and at the end QC subtracted

   TGraph    *Background; // this is the background I determine from the signal thorugh a spline
   TF1		 *Qfit; // The Qcurve Fit function determin ed from QCana and read in from Qcurve.txt filk
   //
   TGraph    *FitGraph;  // new background fit attempt
   TCanvas	 *GeneralCanvas;  // has the signal and polarization vs time on it
   TCanvas	 *StripCanvas;   // shows polarization vs time
   TCanvas	 *StripCanvas_1; // for TE measurement shows the calibration constant over time.
   TCanvas	 *StripCanvas_2; // for TE measurement shows the pressure over time.
   TCanvas	 *StripCanvas_3; // for TE measurement shows the system tempearture over time.
   TCanvas	 *StripCanvas_4; //  background over time.
   TCanvas	 *StripCanvas_5; //  Crate temperature

   TCanvas	 *AuxCanvas;   // all the auxiliary plots, like Qcurve
   TCanvas	 *RTCanvas ;    // Life time display canvas, get used in Loop
   TCanvas	 *DebugCanvas;

   TChain    *NMRchain; // if we have more than one root file
   TString timestring; // from labview time
   TDatime *td ; // Datetime for time histogram
   TTimeStamp *RootTimeStamp;
   TFile *QcurveFile; //
   NMRFastAna* fastAna; // the fast system
   NMR_DST *DST;// handles output
   TFile *mydst; //DST file
   std::string QcurveFileName ; // Qcurve file name
   std::ofstream teout ;   // global teoutput file
   std::ofstream teout1 ;   // global teoutput file for temporary spectrum


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
   virtual void	    PrintWarnings();
   virtual void 	Finish(); // used to clean up memory
   virtual void		SetEnvironment(std::string);
   virtual void		ReadQcurveParFile(std::string);
   virtual void     DrawFitHisto(); // draws signal hist with subtracted fitted background
   virtual void		InitFastAna(); // initializes the fast Qcurve treatment

   // memebre which will be inhertied from TEana
   void     Loop(); // TEana inherits
   void		SetupCanvas();
   TH1D * 	SetupStripChart(TString);
   TH1D *    CopyHisto(const char *, const char *,TH1D *); // Instead of cloning a histo really does a copy
   void	    SetupHistos();
   void		ReadParameterFile(TString );
   void     Init(TTree *tree);
   Bool_t   Notify();
   Int_t    GetEntry(Long64_t entry);
   Long64_t LoadTree(Long64_t entry);
   Int_t 	OpenFile(TString);
   void		CloseFile();
   void		DrawHistos();
   void		Stripper(Long64_t); // darws the stri charts
   Double_t CalculatePatArea(TH1D *); // calculates the difference of left and right background areas
   static Double_t FitFcn2(Double_t * , Double_t *);

   TString  GetDate(TString input);


};
#endif

#ifdef NMRana_cxx


NMRana::NMRana(){
	time_offset = 2082852019 ; // unix offset in seconds
	TimeControl = 2; // always assume the file is from the newest generation // if not use TIMEC = value in parameter file
	QC_DISP = true;
    NumberOfStripCharts=-1; // then start with 0 in loop
    CalConst = 1; // set calibration to default 1
    //StripLength =  68400000; for very long runs.
    //slows down program
    StripLength =  2000;
    PrintWarnings();
    //currently set all the gain arrays to 1
    //pat thinks we only will run on one setting
    gain_array[0]=1.;
    gain_array[1]=1.;
    gain_array[2]=1.;
    QfitPar[0]=0.;
    QfitPar[1]=0.;
    QfitPar[2]=0.;
    QfitPar[3]=0.;
    QfitPar[4]=0.;
    TotalEntries = 0;
    QCshift = true ; // default calculate Qcurve shift.
    Ave = true; // Default we do average over  measurement
    AverageNumber =1;




}


//!!!!!!!!!!!!!!!!!!!!!!!!!start of main loop!!!!!!!!!!!!!!!!!!

void NMRana::Loop()
{


	if (fChain == 0) return;


   	   // go to strip chart
   StripCanvas->cd();
   if(TEmeasurement){
	   StripCanvas_1->cd();
	   StripCanvas_2->cd();
	   StripCanvas_4->cd();
	   StripCanvas_5->cd();

   }

   Long64_t nentries = fChain->GetEntries(); //
   Long64_t time_prev = 0;
   Long64_t nbytes = 0, nb = 0;
   // insert Kun's fst analyzer to get the QCurve offset
   if(QC){
	for(int i = 0; i < nentries; ++i)   //This is equivalent to NMRana::Loop()
	{
		fChain->GetEntry(i);

		TH1D* teHist = new TH1D("teHist", "teHist", ScanPoints, MinFreq, MaxFreq);
		for(int j = 0; j < ScanPoints; ++j)
		{
			teHist->Fill(MinFreqR+j*FreqStep, array->at(j)/gain_array[int(Gain+0.01)]);
		}
		//cout << "Loop " << i << ": xOffset = " << fastAna->getXOffset(teHist) << ", yOffset = " << fastAna->getYOffset() << endl;
		xoffset.push_back(fastAna->getXOffset(teHist)); // fill vector of xoffsets
		yoffset.push_back(fastAna->getYOffset()); // fill vector of xoffsets
		delete teHist;
	}
} // end of fastana

   RTCanvas->cd(1);
   NMR_RT->Draw("HIST P");
   RTCanvas->cd(2);
   Int_t counter_per=0;
   NMR_RT_Corr->Draw("HIST P");
   RTCanvas->cd(3);
   Int_t help_counter = 0; // this counter is reset whenever we reach the averagenumber for averaging over TE signal
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
   //for (Long64_t jentry=0; jentry<20;jentry++) {
	   //if(jentry ==3)break; // only do 3 loops
	   //tell how much has been analyzed
	   Reached_Average = false;
	   Int_t modnum = 100;
	   if(jentry%modnum == 0){

		   Int_t percent= counter_per*modnum*100/nentries;
		   cout<<NMR_pr<<percent <<"  percent of total data analyzed  \n";
		   counter_per++;
	   }

	   // reset the histos every
	   if(Ave && help_counter ==0){
		   NMR_RT->Reset();
		   NMR_RT_Corr->Reset();
		   NMR_RT_Corr_Fit->Reset();
		   NMR1_Qfit->Reset();
	   }
	   //FitGraph->Clear();

	   Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

//	Consistency checks:
// first test that tune voltage between measurement and Qcurve are the same
      if(QC){
      if(NMRchan != QCcoil || NMRchan != QcurveCoil) {
    	  cout<<NMR_pr<<" !!error QCurve coil not the same as TE: Data coil"<< NMRchan<<"    Qcurve voil : "<<QCcoil<< "Fit Qcurve coil :"<<QcurveCoil<<endl;
    	  cout<<NMR_pr<<"!!!!!!!!!!!sever failure!!!!!!!"<<endl;
    	  //exit(EXIT_FAILURE);
         }
      }
      if(TuneV != QcurveTune || TuneV != QCtune) {
    	  cout<<NMR_pr<<" !!warning Tune voltages are not the same Data tuneV :"<< TuneV<<"    Qcurve data tunev : "<<QCtune<< "Fit Qcurve tune :"<<QcurveTune<<endl;
      }




//	now fill histogram
// reset freq_temp to lower bound
      Double_t freq_temp = MinFreqR;

      Double_t DataTemp = 0.; // holds the point of data at position j of data array
      Double_t QcurTemp; // holds the point of data at position j of Qcurev array
      for (UInt_t j = 0; j < array->size(); ++j) {
    	  // subtract QCurve if existing
    	  //renormailze signal by amplifier setting
       	  DataTemp = array->at(j);
          if(QCshift){
           	  // now shift the qcurve. the j has to be shifted by xoffset(jentry); however make sure we do not go beyond the bumdary
              Int_t shift = j-xoffset.at(jentry);
           	  // now shift the qcurve. the j has to be shifted by xoffset(jentry); however make sure we do not go beyond the bumdary
               if(shift<Qcurve_array.size() and shift >=0 ){ QcurTemp = Qcurve_array.at(shift);
               DataTemp = DataTemp-yoffset.at(jentry); //subtract the fitted offset
              }
               else QcurTemp = Qcurve_array.at(j); ;
               //cout<<"problem with shift"<<shift<<"  "<< "entry  "<<j<< "\n";
          }
          else if(Qcurve_array.size()!=0)QcurTemp = Qcurve_array.at(j);


        	  if(DEBUG==1)raw_histo->Fill(freq_temp,DataTemp);

    	  if(Qcurve_array.size()!=0) {
    		 NMR1_NoQ->Fill(freq_temp,DataTemp);

             // cout<<" data  "<<DataTemp<<"  Qcurve  "<<QcurTemp<<"\n";
    		 NMR1->Fill(freq_temp,DataTemp- QcurTemp);

    		 // correct for the Qcurve shift
    		 if(QCshift) NMR1_Qfit->Fill(freq_temp, DataTemp-Qfit->Eval(freq_temp-xoffset.at(jentry)*FreqStep));
       		     else NMR1_Qfit->Fill(freq_temp, DataTemp-Qfit->Eval(freq_temp));
           	 NMR_RT_Corr->Fill(freq_temp,(DataTemp- QcurTemp));
       		 NMR_RT_Corr_Fit->Fill(freq_temp,(DataTemp- QcurTemp));
     		 // now take background out
     	  	  }
    	  else{
    		  NMR1->Fill(freq_temp,DataTemp);
    		  NMR_RT_Corr->Fill(freq_temp,DataTemp);
        	  NMR_RT_Corr_Fit->Fill(freq_temp,DataTemp);


    	  	  }
          NMR_RT->Fill(freq_temp,DataTemp);

          freq_temp = freq_temp+FreqStep;


    	  // for backgroiund graph
      	  }
      StripCanvas->cd();

      	 // herwe calculate the average left and right of the peak for Patarea
      PatDiffArea = CalculatePatArea(NMR_RT_Corr);




// draw the signal histogram
      // only do this for the average number
      if(Ave && help_counter == AverageNumber-1) {

    	  // rescale histos by averagenumber



	   RTCanvas->cd(1);
	   NMR_RT->Draw("HIST P");
	   RTCanvas->cd(2);
	  // copy histo so that we can fit the background without affecting the original hist
           for(Int_t ihelp=1;ihelp<= NMR_RT_Corr_Fit->GetNbinsX(); ihelp++){
        	   NMR_RT_Corr_Fit->SetBinError(ihelp, 0.1);
           }
			  TF1 * FitBckFunc = NewFitBackground(NMR_RT_Corr_Fit);




		 if(TEmeasurement){ SignalArea = CalculateArea(NMR_RT_Corr_Fit);
		 //cout<<NMR_pr<<"signal area    "<<SignalArea<<endl;
		 //cout<<CalculateArea(NMR_RT_Corr_Fit)<<" temp    "<<CalculateArea(NMR_RT_Corr)<<"\n";

		 }
	      else  SignalArea = CalculateArea(NMR_RT_Corr_Fit);



	   NMR_RT_Corr_Fit->Draw("HIST P");
	   RTCanvas->cd(3);
	   NMR_RT_Corr->Draw("HIST P");

	   FitBckFunc->Draw("SAME");



	   RTCanvas->Modified();
	   RTCanvas->Update();

	   // scale histos
	   NMR_RT_Corr_Fit->Scale(1./Double_t(AverageNumber));
	   NMR_RT_Corr->Scale(1./Double_t(AverageNumber));
	   NMR1_Qfit->Scale(1./Double_t(AverageNumber));

	   mytr->Fill();


      // reset the help_counter
	   help_counter = -1;
	   Reached_Average = true;


	   // fill the csv file
	   if(TEmeasurement) {
			  Double_t temp_pol =TE.CalculateTEP("proton",.5,5.0027,HeP);
			  teout<<timel/10000<<","<< ScanNumber<<","<<HeT<<",1.,1.,"<<CalculateArea(NMR_RT_Corr_Fit)*FreqStep<<","<<CalculateArea(NMR_RT_Corr)*FreqStep<<","<<temp_pol<<"\n";
		  }

      }
      help_counter++;
	     //warninghook
	//end warninghook

	  // write out csv file

	      // Convert to polarization

	            SignalArea *=CalConst;

	// now for every point in a TE we will calculate the polarization from the pressure
	// the ratio of calculated polarization/ area gives the calibration constant calib
	// so that calib*area = polarization of the reL SIGNAL
	// at the end we will calculate an average caibration constant with a deviation



	  if(DEBUG ==2)cout<<NMR_pr<<timel<<"another one \n";
	   Stripper(jentry); // draw all the strip charts


  }// end of entry loop


	  	  	 if(TEmeasurement) TE.ShowDistribution(CalibConstantVector);
	  	  	 // now we need to put the sweep number away, so that we know how to rescale the Qcurve
	  	  	 // distribution for display.
	  	  	 // the Qcurve should be subratced on a sweep to seep basis.
		  	 NumberOfSweeps = nentries;

   //scale NMR1 and NMR1_NoQ to one sweep. This might be a problem for having many files
		  	TotalEntries += nentries;
}



void NMRana::DrawHistos(){

	// analyze spectra
	// FindPeak(NMR1);  // temporarily removed, does only seem to be a waste of time
// draw histos, mainly for debug purpose
	if(Qcurve_array.size()!=0){
		if(QC_DISP){

			// now find how much to shift the qcurve
			// determine the average of all offsets
			Int_t temp_shift = 0;

			for (Int_t k=0;k < xoffset.size() ; k++){
				temp_shift +=xoffset.at(k);
			}
			temp_shift = temp_shift/xoffset.size();




			temp_shift = -21;
			cout<< " temporary shift "<<temp_shift<<endl;
	 		 	for (Int_t k=0 ; k < Qcurve_histo->GetNbinsX(); k++){

	 		 		Qcurve_histo_shifted->Fill(Qcurve_histo->GetBinCenter(k)-temp_shift*FreqStep,Qcurve_histo->GetBinContent(k));
	 		 	}


			Qcurve_histo->SetLineColor(3);
			//Qcurve_histo->SetMarkerStyle(21);
			Qcurve_histo_shifted->SetLineColor(4);
			//Qcurve_histo_shifted->SetMarkerStyle(22);
			NMR1_NoQ->SetLineColor(1);
			//NMR1_NoQ->SetMarkerStyle(23);

			AuxCanvas->cd(2);
			Qcurve_histo->GetXaxis()->SetRangeUser(fit_x1,fit_x4);
			Qcurve_histo->Draw("HIST ");
			NMR1_NoQ->Draw("SAME HIST ");
			Qcurve_histo_shifted->Draw("SAME HIST ");

			AuxCanvas->cd(1);
			NMR1_NoQ->GetXaxis()->SetRangeUser(fit_x1,fit_x4);
			NMR1_NoQ->Draw("HIST ");
			AuxCanvas->cd(3);
			if(DEBUG == 1){
				cout<<MinFreq<<"   "<<MaxFreq<<endl;
				raw_histo_QC = CopyHisto("raw_histo_QC","Qcurve subtracted signal",NMR1_NoQ);
				//raw_histo_QC=(TH1D*)NMR1_NoQ->Clone("raw_histo_QC");
				// now subtract the Qcurve
				Double_t NE = 1/TotalEntries;

				raw_histo_QC->Scale(NE);
				raw_histo_QC->Add(Qcurve_histo_shifted,-1.);
				raw_histo_QC->GetXaxis()->SetRangeUser(fit_x1,fit_x4);
				raw_histo_QC->Draw("HIST");
			}

		AuxCanvas->Update();
		}
		else {
			AuxCanvas->cd(1);
			NMR1_NoQ->Draw("HIST P");
			AuxCanvas->Update();
		}

	}

	GeneralCanvas->Divide(1,3);
	// renormalize the histos by the total number of entries we read in
	Double_t NE = 1/TotalEntries;
	NMR1->Scale(NE);
	if(QC){
		NMR1_Qfit->Scale(NE);

	NMR1_NoQ->Scale(NE);
    }

	GeneralCanvas->cd(1);

	NMR1->GetXaxis()->SetRangeUser(FreqCenter*.9986,FreqCenter*1.0014); //set it to same axis as the next histogram
	NMR1->Draw("HIST ");
	//
		if(QC) {
			Qcurve_histo_shifted->Draw("HIST  SAME"); //
		 }



		//	BckFct1->Draw();
//    FitSpectrum(NMR1,1);
	GeneralCanvas->cd(2);
	NMR1_Qfit->Draw("HIST ");


	GeneralCanvas->cd(3);
	PolTime->Draw();
	GeneralCanvas->Update();






   // DrawFitHisto();
	// print out the different integrals;
	// they are unnormalized
	cout<<NMR_pr<< "*************   IntegrGetdateal of histo NMR1 "<<NMR1->Integral(low_id,high_id)<<endl;
	if(QC)cout<<NMR_pr<< "*************   Integral of histo NMR1_NoQ  "<<NMR1_NoQ->Integral(low_id,high_id)<<endl;

}


//!!!!!!!!!!!!!!!!!!!!!!!!!end of main loop!!!!!!!!!!!!!!!!!!







//!!!!!!!!!!!!!!!!! Initialization block !!!!!!!!!!!!!!!!!!!!!!!!!


void NMRana::SetEnvironment(std::string environment){
	NMR_ROOT = environment ;
}

void NMRana::InitFastAna(){
	// initializes the fast ana system, which is used for the QCurve shifting
    cout<<NMR_pr<< "Initializing NMRFastAna" <<endl;

	fastAna = new NMRFastAna();
	fastAna->setFreqRange(212.62, 213.38);
	fastAna->setExclusionWin(213.0, 0.1);
	fastAna->setSampleRate(1);
	fastAna->setXOffsetRange(-50, 50);
	fastAna->setYOffsetRange(-1., 1.);
	fastAna->init();



	//Set the QCurve to fastAna
	fastAna->setQCurve(Qcurve_histo);




}

void NMRana::ReadParameterFile(TString ParameterFile){
	// reads in parameters for running the NMRanalyzer
	// needs the Qcurve file
	// calibration constants from TE measurements
	// line length is a maximum of 80

	// the format of the file is  name , value
	// example: Qcurve   Qcurve.root
	std::string temp;
	std::string string1, string2;
	std::string temp_file;

	Double_t lownmr(212.8);
	Double_t hinmr(213.3);

	ifstream ParFile; // create instream file
	ParFile.open(ParameterFile);
	// check if found and opened
	if(!ParFile.is_open()){
		cout<<"error with parameter file, most likely wrong path   "<<ParameterFile<<" \n";
		cout<< " will terminate \n";
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

		if(pos->first.find("ANATEXT")!= std::string::npos){
			// amplifier setting for QCurve
		std::string filenamestring = pos->second;
		analysis_text_file = filenamestring.c_str();
		cout<<endl;
		cout<<endl;
		cout<<endl;
		cout<<NMR_pr<<"**************************************************"<<endl;
		cout<<NMR_pr<<"!!!!!!!!!!input redirected to "<< analysis_text_file<<" !!!!!!!"<<endl;
		cout<<NMR_pr<<"**************************************************"<<endl;
		cout<<endl;
		cout<<endl;
		cout<<endl;

		freopen(analysis_text_file,"w",stdout); // redirect cout

		TXT=true;

		}

		cout<<NMR_pr<<"parameters from file  :"<<pos->first<<"\t"<<pos->second <<"\n";



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
		if(pos->first.find("CALCONST")!= std::string::npos){
			// amplifier setting for QCurve
		CalConst = std::stod(pos->second);

		}

		if(pos->first.find("QCURVE")!= std::string::npos){
			// amplifier setting for QCurve
		QcurveFileName = pos->second;
		QC=true;

		cout<<"Qcurve File "<<QcurveFileName<<endl;

		}
		if(pos->first.find("AVERAGE")!= std::string::npos){
			// amplifier setting for QCurve
		AverageNumber = std::stoi(pos->second);
		Ave = true;

		cout<<"We are averaging over  "<<AverageNumber<<endl;

		}
		if(pos->first.find("QCSHIFT")!= std::string::npos){
			// shift qcurve or not
		Int_t help = std::stod(pos->second);
		if(help == 1) QCshift = true;
		else QCshift = false;
		}


	}
	if(QC){
		GetQcurve(QcurveFileName);
		//ReadQcurveParFile("QCurveParDec16.csv");
		ReadQcurveParFile("QCurve.txt");
	}


		// set limits for are either default or what comes in from parameter file
	if(lownmr > hinmr){
		cout<<NMR_pr<<"the integration limits are wrong Low limit is larger than high limit  "<<lownmr<<"  "<<hinmr<<"\n";
	}
	else  AreaSetLimits(lownmr,hinmr);
	if(FitLimit){
		SetFitLimits(fit_x1,fit_x2,fit_x3,fit_x4);
		low_fit_x = (fit_x2 + fit_x1)/2.;
		high_fit_x = (fit_x4 + fit_x3)/2.;
	}


}



void NMRana::ReadQcurveParFile(std::string name){
	// this routine reads the QCurve parameter file which has been generated by QCana
	// the format is: QCurve filename
	// par [0]-par[5]
	// min found of fit fucntion
	// tune voltage
	// Note if you ever open the file with microsoft excel, make sure you DO NOT save it
	// MS puts line breaks in, which getline won't understand.

	// open tune file
	std::string file;
	Double_t par[5];
	std::string dummy;
	std::string filename = NMR_ROOT+"/QC_files/"+name;




	std::ifstream Qpar(filename);
	if(!Qpar.is_open()){
		cout<<NMR_pr<<"Cannot find the QCurve parameter file   "<<filename<<"   from QCana \n, will abort \n ";
		exit(EXIT_FAILURE);
	}

    int cols = 11; // number of columns
	std::string line; // problem with line ending is that ffffing excel uses carriage return instead of line feed for files
	while(getline(Qpar,line)){
		size_t len = line.length();
		size_t tmp =0;
		if(line.empty()) {continue;}

		const char* last_start =&line[0];
		int num_parts = 0 ;

			// now check for comma
			while(tmp < len)
			{
				if(line[tmp]== ',' || line[tmp] == '\n' || line[tmp] == '\r' )
				{

					line[tmp] = '\0';
					if(num_parts == cols) {break;}
					else if(num_parts == 0) {file = last_start;}
					else if(num_parts == 1) {QcurveCoil = atoi(last_start);}
					else if(num_parts == 7) {QcurveMin = atof(last_start);}
					else if(num_parts == 8) {QcurveTune = atof(last_start);}
					else if(num_parts == 9) {QcurveScale = atof(last_start);}
					else if(num_parts == 10) {QcurveGain = atof(last_start);}
					else {par[num_parts - 2] = atof(last_start);}
					tmp++;
					num_parts++;
					last_start = &line[tmp];
				}
				tmp++;


			}


		// if we find the correct qcurve  we extract the parameters and exit the loop.
		if(file == QcurveFileName)
		{
			cout<<NMR_pr<< " We are using the following  parameters for the fit function \n\n";
			cout<<file<<" ";
			for (int l =0 ; l<5;l++){
				cout<<"  "<<par[l]<<" ";
				QfitPar[l] = par[l]/QcurveScale/gain_array[int(QcurveGain+.01)]; // normalize to gain
			} // normalzied to number of sweeps


		 //		 	 cout<<Min<<" "<<Qcurve_tune<<" "<<Qcurve_scale<<endl;
		 	 break;


		}
		if(Qpar.eof()) break;
	}


}

int NMRana::OpenFile(TString rootfile){

	// oepn file and initialize tree
     cout<<NMR_pr<<"opening file "<<rootfile<<"\n";

     if(rootfile.Contains("TER") || rootfile.Contains("TEQ") ){
    	 TEmeasurement = true;
    	 cout <<NMR_pr<< "\n\n this is a TE measurement \n\n\n";
    	 std::string tefile = NMR_ROOT+"/TE/TEglobal.csv";
    	 cout<< NMR_pr<<"Opening TE global csv  file"<<tefile<< "\n";
    	 std::string tefile1 = NMR_ROOT+"/TE/TEglobal_spectrum.csv";
    	 cout<< NMR_pr<<"Opening TE global csv  file"<<tefile<< "\n";

    	 teout.open(tefile,ios::out | ios::app) ;
    	 teout1.open(tefile1,ios::out | ios::app) ;

    	 // perform the TEfile read as well, so that we have the map
    	 // this is only used if we use the Yurov file and mesaurement
    	//  TE.ReadTE();
     }
     // make the DST file
     Int_t posc =rootfile.First('.');
     DstFile = rootfile;
     DstFile.Insert(posc,'D',3);
     DstFile.Insert(posc+1,'S',3);
     DstFile.Insert(posc+2,'T',3);

     cout<<posc<<"  "<<DstFile<<endl;
	 mydst = new TFile(DstFile,"recreate");
	 mytr = new TTree("Dtree","analyzed variables from NMRanalyzer");



     f = new TFile(rootfile);
     if(!f->IsOpen()){
    	 cout <<NMR_pr<<"Cannot open ROOT File "<<rootfile<<" \n  will exit \n";
    	 exit(EXIT_FAILURE);
     }

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



		Int_t posc =RootFileArray[0].First('.');




		DstFile =RootFileArray[0] ;
	     DstFile.Insert(posc,'D',3);
	     DstFile.Insert(posc+1,'S',3);
	     DstFile.Insert(posc+2,'T',3);
	     DstFile.Insert(posc+3,'G',3);

	     cout<<posc<<"  "<<DstFile<<endl;
		 mydst = new TFile(DstFile,"recreate");
		 mytr = new TTree("Dtree","analyzed variables from NMRanalyzer");




	     if(RootFileArray[0].Contains("TER") || RootFileArray[0].Contains("TEQ") ){

	    	 TEmeasurement = true;
	    	 cout <<NMR_pr<< "\n\n this is a TE measurement \n\n\n";
	    	 // perform the TEfile read as well, so that we have the map
	    	 //TE.ReadTE();
	    	 std::string tefile = NMR_ROOT+"/TE/TEglobal.csv";
	    	 cout<< NMR_pr<<"Opening TE global csv  file"<<tefile<< "\n";
	    	 std::string tefile1 = NMR_ROOT+"/TE/TEglobal_spectrum.csv";
	    	 cout<< NMR_pr<<"Opening TE global csv  file"<<tefile<< "\n";

	    	 teout.open(tefile,ios::out | ios::app) ;
	    	 teout1.open(tefile,ios::out | ios::app) ;

	     }


     //NMRchain->GetObject("NMRtree",tree);
     Init(NMRchain);

   return 0;

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
   fChain->SetBranchAddress("ScanNumber", &ScanNumber, &b_ScanNumber);
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
   if(TimeControl == 2){
	   fChain->SetBranchAddress("Phase_Voltage", &Phase_Voltage, &b_Phase_Voltage);
	   fChain->SetBranchAddress("Peak_Area", &Peak_Area, &b_Peak_Area);
	   fChain->SetBranchAddress("Pol_Calib_Const", &Pol_Calib_Const, &b_Pol_Calib_Const);
	   fChain->SetBranchAddress("Gain", &Gain, &b_Gain);
	   fChain->SetBranchAddress("Pol_Sign", &Pol_Sign, &b_Pol_Sign);
	   fChain->SetBranchAddress("Log_Channel", &Log_Channel, &b_Log_Channel);
	   fChain->SetBranchAddress("Peak_Amp", &Peak_Amp, &b_Peak_Amp);
	 	 fChain->SetBranchAddress("NMRchan",&NMRchan,&b_NMRchan);
	 	 fChain->SetBranchAddress("PeakCenter",&PeakCenter,&b_PeakCenter);
	 	 fChain->SetBranchAddress("BeamOn",&BeamOn,&b_BeamOn);
	 	 fChain->SetBranchAddress("RFlevel",&RFlevel,&b_RFlevel);
	 	 fChain->SetBranchAddress("IFatt",&IFatt,&b_IFatt);
		   fChain->SetBranchAddress("HeT", &HeT, &b_HeT);
		   fChain->SetBranchAddress("HeP", &HeP, &b_HeP);
	    }


//    fChain->SetBranchAddress("ControllerV", &ControllerV, &b_ControllerV);
   fChain->SetBranchAddress("timel", &timel, &b_timel);
   fChain->SetBranchAddress("array", &array, &b_array);
   fChain->SetBranchAddress("IntScanPoints", &IntScanPoints, &b_IntScanPoints);
   Notify();
}


//!!!!!!!!!!!!!!!end initialization block !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Double_t NMRana::FitFcn2(Double_t *x , Double_t *par){
	//  polynomial fit function from QCana

	Double_t y = x[0]-213.0;
     Double_t result = (par[0]+par[1]*x[0]+par[2]*pow(x[0],2.))-(par[3]*pow(y,3)-y*par[4]);
	return result;
}


//!!!!!!!!!!!!!!! finish everything up !!!!!!!!!!!!!!!!!!

void NMRana::CloseFile(){
	f->Close();
	DST->CloseFile();
}

NMRana::~NMRana()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();

}
void NMRana::Finish(){
	// clean up memory
		cout<<NMR_pr<< "***************************************************************"<<endl;
		cout<<NMR_pr<<" *"<<endl;
		cout<<NMR_pr<<"                  Analysis finished"<<endl;
		cout<<NMR_pr<<" *"<<endl;
		cout<<NMR_pr<< "***************************************************************"<<endl;
		//DST->WriteTree();
		cout<<NMR_pr<< "Everything initialized"<<endl;
		//mydst = new TFile(DstFile,"recreate");


		mydst->cd();
	    mytr->Write();
		mytr->Print();
		NMR1_NoQ->Write();
		NMR1->Write();
		NMR1_Qfit->Write();
		Qcurve_histo->Write();
		NMR_RT_Corr_Fit->Write();
		mydst->ls();
		mydst->Write();
		mydst->Write();
		 Long64_t help = mydst->GetFileBytesWritten();
		cout<<" tree file writing "<<help<<endl;

		mydst->Close();

		teout.close(); // close te file
		teout1.close();

}


//!!!!!!!!!!!!!!! end of finish everything up !!!!!!!!!!!!!!!!!!



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


	   Long64_t ientry = LoadTree(0);

	   fChain->GetEntry(0);

	   fChain->Show(0);


	   //GetEntry(0);
	   //Show(0);


	   // histo limits
	   MinFreq = (FreqCenter) - (ScanPoints-1.)/2. * FreqStep-FreqStep/2.;
	   MinFreqR = (FreqCenter) - (ScanPoints-1.)/2. * FreqStep;
	   MaxFreq = FreqCenter + (ScanPoints-1.)/2. * FreqStep+FreqStep/2.;
	   MaxFreqR = (FreqCenter) + (ScanPoints-1.)/2. * FreqStep;

	   NMR1 = new TH1D("NMR1","Signal histogram",IntScanPoints,MinFreq,MaxFreq);
	   NMR1->Sumw2();
	   NMR1_Qfit = new TH1D("NMR1_Qfit","Signal histogram with Qcurve fit subtracted",IntScanPoints,FreqCenter*.9986,FreqCenter*1.0014);
	   NMR1_Qfit->Sumw2();
	   NMR_RT = new TH1D("NMR_RT","Real TimeSignal histogram",IntScanPoints,MinFreq,MaxFreq);
	   NMR_RT->Sumw2();
	   NMR_RT->SetLineColor(kSpring-2);
	   NMR_RT->SetLineWidth(4);
	   NMR_RT->GetXaxis()->SetRangeUser(fit_x1,fit_x4);
	   NMR_RT_Corr = new TH1D("NMR_RT_Corr","Real TimeSignal Qcurve subtracted",IntScanPoints,MinFreq,MaxFreq);
	   NMR_RT_Corr->Sumw2();
	   NMR_RT_Corr->SetLineColor(kBlue-2);
	   NMR_RT_Corr->SetLineWidth(4);
	   NMR_RT_Corr->GetXaxis()->SetRangeUser(fit_x1,fit_x4);

	   NMR_RT_Corr_Fit = new TH1D("NMR_RT_Corr_Fit","Real TimeSignal Fit background subtracted",IntScanPoints,MinFreq,MaxFreq);
	   NMR_RT_Corr_Fit->Sumw2();
	   NMR_RT_Corr_Fit->SetLineColor(kBlue-2);
	   NMR_RT_Corr_Fit->SetLineWidth(4);
	   NMR_RT_Corr_Fit->GetXaxis()->SetRangeUser(fit_x1,fit_x4);
//	   Double_t xgr[1000], ygr[1000];
//	   FitGraph = new TGraph(401,xgr,ygr);

	   // Determine the Integration or summation limits for peak in terms of channels.
	   	  // low_id = NMR1->GetXaxis()->FindBin(LowArea_X);
	   	  // high_id = NMR1->GetXaxis()->FindBin(HighArea_X);
	   	   low_id = NMR_RT_Corr_Fit->GetXaxis()->FindBin(fit_x2);
	   	   high_id = NMR_RT_Corr_Fit->GetXaxis()->FindBin(fit_x3);
	   Ifit_x1 = NMR_RT_Corr_Fit->GetXaxis()->FindBin(fit_x1);
	   Ifit_x2 = NMR_RT_Corr_Fit->GetXaxis()->FindBin(fit_x2);
	   Ifit_x3 = NMR_RT_Corr_Fit->GetXaxis()->FindBin(fit_x3);
	   Ifit_x4 = NMR_RT_Corr_Fit->GetXaxis()->FindBin(fit_x4);

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
	   if(!TEmeasurement){
		   PolTime = SetupStripChart("Polarization vs time");

	   PolTime->SetMaximum(100.);
	   PolTime->SetMinimum(-100.);
	   PolTime->SetLineColor(2);
	   }
	   //setup calibration constant histogram
	   if(TEmeasurement){
		   PolTime = SetupStripChart("TE area vs time");

		   PolTime->SetMaximum(100.);
		   PolTime->SetMinimum(-100.);
		   PolTime->SetLineColor(2);

		   CalibTime = SetupStripChart("Calibration Constant (pol over TEarea) vs time");
		   CalibTime->SetMaximum(100.);
		   CalibTime->SetMinimum(0.);
		   PressTime = SetupStripChart("TE pressure [Torr] vs time");
		   PressTime->SetMaximum(20.);
		   PressTime->SetMinimum(0.);
		   SysTempTime = SetupStripChart("TEmperature vs time");
		   SysTempTime->SetMaximum(20.);
		   SysTempTime->SetMinimum(0.);
		   PatTemp = SetupStripChart("Background difference vs time");
		   PatTemp->SetMaximum(20.);
		   PatTemp->SetMinimum(-20.);
		   CrateTemp = SetupStripChart("Crate Temp vs time");
		   CrateTemp->SetMaximum(30.);
		   CrateTemp->SetMinimum(20.);

	   }




	   // Print original time
	   PrintTime();


	   Int_t timewindow = 1;
	   if(!TEmeasurement){
		   PolTime = new TH1D("PolTime","Polarization vs time",timewindow,0,10*timewindow);

	   PolTime->GetXaxis()->SetTimeDisplay(1);
	   PolTime->GetXaxis()->SetTimeOffset(TimeStamp);
	   PolTime->GetXaxis()->SetTimeFormat("%d %m %H:%M :%S");
	   PolTime->GetXaxis()->SetNdivisions(405) ;
	   }
	   if(TEmeasurement){
		   PolTime = new TH1D("PolTime","TE area vs time",timewindow,0,10*timewindow);

	   PolTime->GetXaxis()->SetTimeDisplay(1);
	   PolTime->GetXaxis()->SetTimeOffset(TimeStamp);
	   PolTime->GetXaxis()->SetTimeFormat("%d %m %H:%M :%S");
	   PolTime->GetXaxis()->SetNdivisions(405) ;
		   CalibTime = new TH1D("CalibTime","Polarization over TEarea  vs time",timewindow,0,10*timewindow);
		   CalibTime->GetXaxis()->SetTimeDisplay(1);
		   CalibTime->GetXaxis()->SetTimeOffset(TimeStamp);
		   CalibTime->GetXaxis()->SetTimeFormat("%d %m %H :%M :%S");
		   CalibTime->GetXaxis()->SetNdivisions(405) ;
		   CalibTime->SetLineColor(kRed+2);

		   PressTime = new TH1D("PressTime","Pressure  vs time",timewindow,0,10*timewindow);
		   PressTime->GetXaxis()->SetTimeDisplay(1);
		   PressTime->GetXaxis()->SetTimeOffset(TimeStamp);
		   PressTime->GetXaxis()->SetTimeFormat("%d %m %H :%M :%S");
		   PressTime->GetXaxis()->SetNdivisions(405) ;
		   SysTempTime = new TH1D("SysTempTime","System temeperature  vs time",timewindow,0,10*timewindow);
		   SysTempTime->GetXaxis()->SetTimeDisplay(1);
		   SysTempTime->GetXaxis()->SetTimeOffset(TimeStamp);
		   SysTempTime->GetXaxis()->SetTimeFormat("%d %m %H :%M :%S");
		   SysTempTime->GetXaxis()->SetNdivisions(405) ;
		   SysTempTime->SetMaximum(35.);
		   SysTempTime->SetMinimum(25.);

		   PatTemp = new TH1D("PatTemp","System temeperature  vs time",timewindow,0,10*timewindow);
		   PatTemp->GetXaxis()->SetTimeDisplay(1);
		   PatTemp->GetXaxis()->SetTimeOffset(TimeStamp);
		   PatTemp->GetXaxis()->SetTimeFormat("%d %m %H :%M :%S");
		   PatTemp->GetXaxis()->SetNdivisions(405) ;

	   }

	   // Now setup a histo for Qcurve




	   if(Qcurve_array.size()!=0){
		   Qcurve_histo = new TH1D("Qcurve_histo","Normalized Qcurve histogram",IntScanPoints,MinFreq,MaxFreq);
		   Qcurve_histo->Sumw2();
		   Qcurve_histo_shifted = new TH1D("Qcurve_histo_shifted","Normalized Qcurve histogram shifted ",IntScanPoints,MinFreq,MaxFreq);
		   Qcurve_histo_shifted->Sumw2();
		   NMR1_NoQ = new TH1D("NMR1_NoQ","Signal without QCurve subtraction",IntScanPoints,MinFreq,MaxFreq);
		   NMR1_NoQ->Sumw2();
			// Now fill the QCurve histo if it is used
				  if(Qcurve_array.size()!=0){
				      Double_t freq_temp = MinFreqR; //shift everything by 22 steps


				      for (UInt_t j = 0; j < Qcurve_array.size(); ++j) {
				    	  // subtract QCurve if existing

				          Qcurve_histo->Fill(freq_temp,Qcurve_array.at(j));
				          freq_temp = freq_temp+FreqStep;
				      	  }

				  }
			if(QfitPar[0] !=0)	{   // we have a qfit from the file
				// Now Create Qfit

				Qfit= new TF1("Qfit",FitFcn2,FreqCenter*.9986,FreqCenter*1.0014,5);
					//	FitH2->SetParameters(7e7,6e5,1500,1116,107);
					 	Qfit->SetParameters(QfitPar);

					 	//cout<<"************************** value for 212.71 :"<<Qfit->Eval(212.71)<<"value for 213.26 : "<<  Qfit->Eval(213.26)<<"\n\n" <<endl;

				}



	   }

	   // Do the graph for the histo
	  // Background = new TGraph(IntScanPoints,gr_freq,gr_amp);

	   	   // histo for debugging
	   if(DEBUG==1)	 {
		   raw_histo = new TH1D("raw_histo","Raw Signal histogram",IntScanPoints,MinFreq,MaxFreq);
		   raw_histo_QC = new TH1D("raw_histo_QC"," Signal histogram with QCurve subtracted",IntScanPoints,MinFreq,MaxFreq);//
	   }

	   // Initialze the QCurve Fast Ana system, but only if there is aQcurve

	   if(QC) InitFastAna();


		mytr->Branch("SignalArea",&SignalArea,"area/D");
		mytr->Branch("timel",&timel,"timel/L");
		mytr->Branch("NMR_RT_Corr_Fit","TH1D",&NMR_RT_Corr_Fit,128000,0);
		mytr->Branch("NMR_RT_Corr","TH1D",&NMR_RT_Corr,128000,0);
		mytr->Branch("NMR_RT","TH1D",&NMR_RT,128000,0);
		mytr->Branch("NMR1_Qfit","TH1D",&NMR1_Qfit,128000,0);



}
void NMRana::SetupCanvas(){
	// creates the different Canvas
	// master canvas for all histograms
	GeneralCanvas = new TCanvas("GeneralCanvas","NMR signal",500,50,400,400);

    Int_t strip_w = 400;
    Int_t strip_x =1000;
    Int_t strip_y =50;
    Int_t strip_h = 100;
	//canvas for all strip charts
	if(!TEmeasurement) StripCanvas =  new TCanvas("StripCanvas","NMR strip charts",strip_x,strip_y,strip_w,strip_h);
	else 		StripCanvas =  new TCanvas("StripCanvas","NMR strip charts",strip_x,strip_y,strip_w,strip_h);
	StripCanvas->SetGrid();
	StripCanvas->SetFillColor(42);
	StripCanvas->SetFrameFillColor(33);
    // only do a strip chart for pressure if we have a TE measurement.
	if(TEmeasurement){
		strip_y += strip_h;
		//StripCanvas_1 =  new TCanvas("StripCanvas_1","Calibration Constant strip charts",1200,350,1200,200);
		StripCanvas_1 =  new TCanvas("StripCanvas_1","Calibration Constant strip charts",strip_x,strip_y,strip_w,strip_h);
		StripCanvas_1->SetGrid();
		StripCanvas_1->SetFillColor(40);
		StripCanvas_1->SetFrameFillColor(30);
		strip_y += strip_h;
		StripCanvas_2 =  new TCanvas("StripCanvas_2","Pressure strip charts",strip_x,strip_y,strip_w,strip_h);
		StripCanvas_2->SetGrid();
		StripCanvas_2->SetFillColor(38);
		StripCanvas_2->SetFrameFillColor(28);
		strip_y += strip_h;
		StripCanvas_3 =  new TCanvas("StripCanvas_3","System Temperature strip charts",strip_x,strip_y,strip_w,strip_h);
		StripCanvas_3->SetGrid();
		StripCanvas_3->SetFillColor(36);
		StripCanvas_3->SetFrameFillColor(26);
		strip_y += strip_h;

		StripCanvas_4 =  new TCanvas("StripCanvas_4","Left Righ Backrond difference strip charts",strip_x,strip_y,strip_w,strip_h);
		StripCanvas_4->SetGrid();
		StripCanvas_4->SetFillColor(34);
		StripCanvas_4->SetFrameFillColor(24);
		strip_y += strip_h;
		StripCanvas_5 =  new TCanvas("StripCanvas_5","Crate temperature",strip_x,strip_y,strip_w,strip_h);
		StripCanvas_5->SetGrid();
		StripCanvas_5->SetFillColor(34);
		StripCanvas_5->SetFrameFillColor(24);

	}

	//RTCanvas =  new TCanvas("RTCanvas","Real Time charts",100,1000,400,300);
	// For debugging purposes
	RTCanvas =  new TCanvas("RTCanvas","Real Time charts",100,200,800,800);
	RTCanvas->SetGrid();
	RTCanvas->SetFillColor(23);
	RTCanvas->SetFrameFillColor(16);
	RTCanvas->Divide(1,3);





	if(Qcurve_array.size()!=0){

	AuxCanvas = new TCanvas("AuxCanvas","Auxiliary plots",500,500,200,200);

	if(QC_DISP)AuxCanvas->Divide(1,3); // if we want to see the Qcurve as well
	else AuxCanvas->Divide(1,1);

	}
	if(DEBUG ==1)	DebugCanvas = new TCanvas("DebugCanvas","Debugging plots",50,700,200,200);

}

void NMRana::DrawFitHisto(){
	// this creates at the end a new function from the original QCurve fit but scaled by the number of entries
	// then it subtracts this new function from the uncorrected histo gram and draws it
	// this way we make sure we have the same limits as the function
	// originally this is just a test to compare the fitted functuion with the actual subtracted qcurve
	       TF1 *Qfit1= new TF1("Qfit1",FitFcn2,FreqCenter*.9986,FreqCenter*1.0014,5);
		//	FitH2->SetParameters(7e7,6e5,1500,1116,107);
	   	   Long64_t nentries = fChain->GetEntriesFast();

		 	Qfit1->SetParameters(QfitPar[0]*nentries,QfitPar[1]*nentries,QfitPar[2]*nentries,QfitPar[3]*nentries,QfitPar[4]*nentries);
		 	TCanvas *c10 = new TCanvas();
		 	TH1D *newhisto = new TH1D("newhisto","test hist of signal with fitted background subtracted ",IntScanPoints,FreqCenter*.9986,FreqCenter*1.0014);
		 	newhisto->Sumw2();
		 	// Now fille the new histo with the signal noq histo gram
		 	// first get number of channels in NMR1_NoQ

		 	for (Int_t k=0 ; k < NMR1_NoQ->GetNbinsX(); k++){
		 		//cout<<k<<"  "<<NMR1_NoQ->GetBinCenter(k)<<"  "<<NMR1_NoQ->GetBinContent(k)<<endl;
		 		newhisto->Fill(NMR1_NoQ->GetBinCenter(k),NMR1_NoQ->GetBinContent(k));
		 	}

		 	c10->Divide(1,2);
		 	//newhisto->Add(Qfit,-1.);
		 	c10->cd(1);
		 	newhisto->Draw("HIST P");
		 	Qfit->SetLineColor(3);
		 	Qfit->Draw("SAME");
		 	c10->cd(2);
		 	Qcurve_histo->Draw();
		 	Qfit->Draw("SAME");
		 	c10->Update();


}


TH1D * NMRana::CopyHisto(const char *hiname , const char *hititle,TH1D *temphist){
	// instead of cloning copies a histo

		Int_t chan = temphist->GetNbinsX();
		Double_t xl = temphist->GetXaxis()->GetXmin();
		Double_t xh = temphist->GetXaxis()->GetXmax();

		cout<< "copy histo"<< xl<< "   "<<xh<<"   "<<chan<<endl;
		TH1D * hi =  new TH1D(hiname,hititle,chan,xl,xh)	;

	 	for (Int_t k=0 ; k < temphist->GetNbinsX(); k++){

	 		hi->Fill(temphist->GetBinCenter(k),temphist->GetBinContent(k));
	 	}
     return hi;


}

Double_t NMRana::CalculateArea(TH1D *histo){
	// this function calculates the area of the NMR perak by simpy summing it in
	//  limits given by AreaSetLimits in the main program
	// the area is calculated as sum i_low to i_high (a(i)*binwidth)

	   // determine sum of channels from low channel to high channel as determined from read in.
    Double_t sum = histo->Integral(low_id,high_id);
    //cout<< histo->Integral(low_id,high_id) <<"  "<<(high_id-low_id+1)*CurveOffset<<"  "<<sum<<endl;

    return sum ;
   // return sum * FreqStep;




}

Double_t NMRana::CalculatePatArea(TH1D *histo){
	// inetgrates lwer and higher background and takes the difference
		Int_t Left_ch = Ifit_x2-Ifit_x1;
		Int_t Right_ch = Ifit_x4-Ifit_x3;
		Double_t low = histo->Integral(Ifit_x1,Ifit_x2)/(Left_ch);
		Double_t high = histo->Integral(Ifit_x3,Ifit_x4)/(Right_ch);
		// calculate the offset shift which is the average of the left window and the right window
		CurveOffset = (high+low)/2.;
		return (high - low);
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
	NumberOfStripCharts++;


	Float_t bintime =1.;
	ht[NumberOfStripCharts] = new TH1D(Form("ht_%d",NumberOfStripCharts),Title,10,0,10*bintime);
	ht[NumberOfStripCharts]->SetStats(0);
	ht[NumberOfStripCharts]->SetLineColor(2);
	ht[NumberOfStripCharts]->GetXaxis()->SetTimeDisplay(1);
	ht[NumberOfStripCharts]->GetYaxis()->SetNdivisions(520);
	return ht[NumberOfStripCharts];



}

void NMRana::GetQcurve(std::string filename) {

	cout<<NMR_pr<<" QCurve file "<<filename<<"\n";
	// Now add the envieronmnet to it
	filename = NMR_ROOT+"/QC_files/"+filename;

	QcurveFile = new TFile(filename.c_str());
	// now create a new tree, so that we have a different name
	QCUtree = (TTree*)QcurveFile->Get("NMRtree");
	Init(QCUtree);
	FillQcurveArray();



}

void NMRana::FillQcurveArray(){
	// this function fill the Qcurve array and normalizes it
	   Long64_t QcurveEntries = QCUtree->GetEntriesFast();
	   Long64_t nbytes = 0, nb = 0;
       if(DEBUG ==1) cout<<NMR_pr<<"FillQcurve>  "<<" QcurveEntries"<< QcurveEntries<<" \n";
       Qcurve_array.clear();
	   for (Long64_t jentry=0; jentry<QcurveEntries;jentry++) {
	      Long64_t ientry = LoadTree(jentry);
	      if (ientry < 0) break;
	      nb = QCUtree->GetEntry(jentry);   nbytes += nb;
	      Double_t freq_temp = MinFreqR;
	      QCtune = TuneV;
	      QCcoil = NMRchan;

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
	   // historicif(DEBUG==1)cout<<NMR_pr<<"FillQcurve>  "<<" qamp"<<Qamp<<"\n";

	   if(DEBUG==1)cout<<NMR_pr<<"FillQcurve>  "<<" qamp"<<gain_array[int(Gain+.01)]<<"\n";

	   for (UInt_t j = 0; j < array->size(); ++j) {

		   Qcurve_array.at(j)= Qcurve_array.at(j)/QcurveEntries/gain_array[int(Gain+.01)];

	   }
	   cout<<" nmr ana  "<<gain_array[int(Gain+.01)]<<"  "<<QcurveEntries<<" "<<ScanNumber<<endl;



}
void NMRana::SetTimeControl(Int_t main_input){
	TimeControl = main_input;
}

TString NMRana::GetDate(TString input) {


	// TString timestring = input;
	timestring = input;
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
      if(Int_t(time_test) > 1465948800 && Int_t(time_test) <= 1468972800 ) TimeControl =1; // this gives a control value for which time the polarization file is from
      if(Int_t(time_test) > 1468972800 ) TimeControl =2; // this gives a control value for which time the polarization file is from

      return  asctime(ltm);
}

void NMRana::PrintWarnings(){
	// this is a routine which prints out current warnings and problems in the code

	// First it prints the version number
	// and then an array of warnings
	// these warnings have to be corrected for later versions:
	   std::string version = "1.01";
	   std::string bold = "\e[1m";
	   std::string nonbold ="\e[0m";
	   std::vector<std::string> warning;

	   std::string NMRwarning = "NMR warning> ";
	   // Block of warning messages
	   warning.push_back(" In NMRana loop, the polarization currently gets calculated from the uncorrected signal ");
	   warning.push_back(" In NMRana loop, The qcurve subtraction is off; needs to be turned on for new NMR signals ");
	   warning.push_back(" In NMRana loop, currently we are calculating no background subtraction for the real signal, this needs to change line 779");
	   warning.push_back(" In NMRana loop, the Control does not work for stripper and read in");


	   // print out
	   cout<<"**********************************************************************************************\n\n";
	   cout<<"welcome to NMRanalyzer version "<<bold<<version<<nonbold<<"\n\n";
	   cout<<"**********************************************************************************************\n\n\n";

	   // Now print out the warning messages
	   for(auto a =warning.begin() ; a != warning.end(); a++){
		   cout<<NMRwarning<<bold<<*a<<nonbold<<"\n";
	   }
	   cout<<"\n\n\n";
}

void NMRana::Stripper(Long64_t jentry){
	// called by Loop and does the strip charts

	  GetTimeStamp();
	  PolTime->SetBinContent(jentry,SignalArea);
	  PolTime->GetXaxis()->SetRange(jentry-StripLength,jentry+20);
	  StripCanvas->Clear();
	  PolTime->Draw();
	  StripCanvas->Modified();
	  StripCanvas->Update();

	  if(TEmeasurement){
		  //if(jentry ==0)cout<<NMR_pr<<"!!!!!!!!!!!!!!!need to change the call to TE polarization calculation!!!!!!!!!!!\n";
		  StripCanvas_1->cd();
		  // temporary fix for separate time file
		  //replace press_help with pressure variable from ROOT file//
		  // give a pressure if there is none
		  //if( HeP==0.0 ) HeP = 5.5 ; //pressure in Torr
		  // take care of the averaging
		  // we can only calculate the calib constant when we have an average
		  // otherwise signal area will be 0
		  // this is because of the averaging we have to do, so we only fill when we havae reached a full average.

		  if(Reached_Average){

			  CalibConstant = TE.CalculateTEP("proton",.5,5.0027,HeP) ; // needs to change to ROOTfile pressure
//			  cout<<" calib constant"<<CalibConstant<<"   area"<<SignalArea<<"\n";
			  CalibConstant = CalibConstant/SignalArea*AverageNumber;
			  CalibConstantVector.push_back(CalibConstant);
			  CalibTime->SetBinContent(jentry,CalibConstant);
			  CalibTime ->GetXaxis()->SetRange(jentry-StripLength,jentry+20);
			  StripCanvas_1->Clear();
			  CalibTime->Draw();
			  StripCanvas_1->Modified();
			  StripCanvas_1->Update();
		  }

		  StripCanvas_2->cd();
		  PressTime->SetBinContent(jentry,HeP);
		  PressTime ->GetXaxis()->SetRange(jentry-StripLength,jentry+20);
		  StripCanvas_2->Clear();
		  PressTime->Draw();
		  StripCanvas_2->Modified();
		  StripCanvas_2->Update();

		  StripCanvas_3->cd();
		  SysTempTime->SetBinContent(jentry,Temperature);
		  SysTempTime ->GetXaxis()->SetRange(jentry-StripLength,jentry+20);
		  StripCanvas_3->Clear();
		  SysTempTime->Draw();
		  StripCanvas_3->Modified();
		  StripCanvas_3->Update();

		  StripCanvas_4->cd();
		  PatTemp->SetBinContent(jentry,PatDiffArea);
		  PatTemp ->GetXaxis()->SetRange(jentry-StripLength,jentry+20);
		  StripCanvas_4->Clear();
		  PatTemp->Draw();
		  StripCanvas_4->Modified();
		  StripCanvas_4->Update();

	  }




}


#endif // #ifdef NMRana_cxx

