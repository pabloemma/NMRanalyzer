/*
 * Ana.h
 *
 *  Created on: May 20, 2016 on an airplane
 *      Author: klein
 */

#ifndef ANA_H_
#define ANA_H_
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <iomanip>

#include <TROOT.h>
#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TString.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TSpectrum.h>
#include <TLegend.h>
#include <TMath.h>
#include <TGraph.h>
#include <TGraph.h>
#include "Math/MinimizerOptions.h"

using namespace std;

class Ana
{
private:
	Bool_t peakfind;
	 Double_t  fit_low_overall; // fit limits
	 Double_t fit_high_overall;
     Double_t fit_limit[4]; ; // these are the limits for the background fit
     Double_t reject1, reject2;// the lower and upper limits of what backrgound fit rejects.
public:
	std::string Ana_pr ="Ana_ana> ";


	TH1D *Spectrum;  // the histogram we want to find the spectrum of
	TSpectrum *spec;
	TF1 *FitFcn;  // convoluted fit function
	TF1 *FitGaus; // first Gauss fit
	TF1 *FitBck; // background fit
	TF1 *FitBckCopy; // copy of background fit
	TF1 *BckFct; // function from background fit
	TF1 *BckFct1; // copy of function from background fit
	TF1 *ff1;


	Int_t npeaks;	//number of peaks to find

	Double_t *xpeaks; // xpositions of peaks
	Float_t sigma; // sigma for peak search
	Double_t amplitude;// amplitude from peak serach
	Float_t Offset; // Offset which is added so that we do not have any negative counts
					// this offset gets determned from the loweest histogram entry and then added if that one is negative
	Int_t NumberOfSpectra; // how many spectra are in the NMR sweep.


	Double_t par[6]; // the fit parameters
	Double_t gaus_par[3];// for initial gaus fit
	Double_t bck_par[4]; // background polynomial (3rd deg polynomial






	Ana();
	virtual ~Ana();


//	 static Double_t BackSpline(Double_t *, Double_t *,Double_t *, Double_t *, TGraph *); // This will determine a spline back ground
	 static Double_t BackSpline( TGraph *);
	 TH1D* FitBackground(TH1D*); // this one gets called at every sweep event
	 void SetFitLimits(Double_t, Double_t, Double_t, Double_t);
	 Double_t Background(Double_t *, Double_t *);
	 Double_t Background_1(Double_t *, Double_t *);
	 static Double_t CopyBackground(Double_t *, Double_t *);
	 TF1 *NewFitBackground(TH1D *);
	 void makeTF1();
	 void makeTF11();



}; //end of class definition

Ana::Ana(){
	//initialize class
	cout<< Ana_pr<<"************* initializing analyzer  ************** \n";
}

Ana::~Ana(){
	cout<<Ana_pr<<"Done with analysis \n\n\n";
}



// Quadratic background function
Double_t Ana::Background(Double_t *x, Double_t *par) {
	// new version with point rejection, the idea being that
	// I will not fit background in the peak area but on left and right side of it
	// see fit descrition in Root reference manual
	// 3rd degree polynomial

		if(x[0]>reject1  && x[0] < reject2){
	     TF1::RejectPoint();
	     return 0;
	 }

   return (par[0] + par[1]*x[0] + par[2]*pow(x[0],2));
		//return (par[0] + par[1]*x[0] + par[2]*pow(x[0],2)+par[3]*pow(x[0],2));
		//return (par[0] + par[1]*x[0]);


}
Double_t Ana::Background_1(Double_t *x, Double_t *par) {
	// new version with point rejection, the idea being that
	// I will not fit background in the peak area but on left and right side of it
	// see fit descrition in Root reference manual
	// 3rd degree polynomial

		if(x[0]>reject1  && x[0] < reject2){
	     TF1::RejectPoint();
	     return 0;
	 }

   return (par[0] + par[1]*x[0] + par[2]*pow(x[0],2));
		//return (par[0] + par[1]*x[0] + par[2]*pow(x[0],2)+par[3]*pow(x[0],2));
		//return (par[0] + par[1]*x[0]);



}
Double_t Ana::CopyBackground(Double_t *x, Double_t *par) {
	// new version with point rejection, the idea being that
	// I will not fit background in the peak area but on left and right side of it
	// see fit descrition in Root reference manual
	// 3rd degree polynomial, this is just a copy of the real backgound but with no rejection
	// this way I can subtract the function from the histogram.




	//return (par[0] + par[1]*x[0] + par[2]*pow(x[0],2)+par[3]*pow(x[0],2));
	return (par[0] + par[1]*x[0] + par[2]*pow(x[0],2));
	//return (par[0] + par[1]*x[0]);


}


void Ana::makeTF1(){
	// this is so that I can pass parameters to the background fit
	//FitBck = new TF1("FitBck",this,&Ana::Background,fit_limit[0],fit_limit[3],4);
	FitBck = new TF1("FitBck",this,&Ana::Background,fit_limit[0],fit_limit[3],3);
}
void Ana::makeTF11(){
	// this is so that I can pass parameters to the background fit
	//ff1 = new TF1("ff1",this,&Ana::Background,fit_limit[0],fit_limit[3],4);
	ff1 = new TF1("ff1",this,&Ana::Background_1,fit_limit[0],fit_limit[3],3);

}


TF1 *Ana::NewFitBackground(TH1D *spectrum){

// this is using a TGraph to fit the spectrum
// first we create a TGraph from a histo by filling x, y vector
// then we fit the background
// and finally we return this graph

	// get info on histo
    Double_t dataAmp[1000] , freq[1000];
	Int_t nchan = spectrum->GetNbinsX();
	for(Int_t k=1; k<=nchan;k++){
		dataAmp[k]=(spectrum->GetBinContent(k));
		freq[k]=(spectrum->GetBinCenter(k));
	}
	//Now create Graph
	reject1 = fit_limit[1];
	reject2 = fit_limit[2];
    makeTF11();
	//ff1 = new TF1("ff1",Background,fit_limit[0],fit_limit[3],3);
	TGraph *gr1 = new TGraph(nchan,freq,dataAmp);
//

	//Now comes the fit part

	ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2","Minimize");
	ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);
	//ROOT::Math::MinimizerOptions::SetDefaultTolerance(.0001);
	ROOT::Math::MinimizerOptions::SetDefaultTolerance(.001);
	ROOT::Math::MinimizerOptions::SetDefaultPrecision(1.e-10);
	ROOT::Math::MinimizerOptions::SetDefaultMaxIterations(100000);
	ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1000000);
	ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(0);

	gr1->Fit("ff1","REWF0");
	return ff1;
}
	/*


	makeTF1();
	FitBck->SetNpx(1000);
	FitBck->SetLineWidth(4);
	FitBck->SetLineColor(kYellow);

	// first we try to get starting pararemets by fitting the background only
		// to the lower background area
		reject1 = fit_limit[0];
		reject2 = fit_limit[1];

		// the upper limit of the lower background window we choose
		gr1->Fit(FitBck,"RE0Q");
		FitBck->GetParameters(&bck_par[0]);
		// now we fit the whole spectrum with the new variables
		reject1 = fit_limit[1];
		reject2 = fit_limit[2];
		gr1->Fit(FitBck,"RE0Q"); // parameter Q for quiet
		FitBck->GetParameters(&bck_par[0]);
		FitBckCopy = new TF1("FitBckCopy",CopyBackground,fit_limit[0],fit_limit[3],3);
		FitBckCopy->SetParameters(bck_par);


    return FitBckCopy;


} // end of newfit
*/

//TH1D* Ana::FitBackground(TH1D *spectrum1){
	TH1D* Ana::FitBackground(TH1D *spectrum){
	// this just determines the backgtound parameters of the spectrum
	// Currently it is a simple 2degree polynomial

   //check for offset
/*	Double_t offset = spectrum->GetBinContent(spectrum->GetMinimumBin());
    if(offset<0.){
    for(Int_t k=0;k<spectrum->GetNbinsX();k++){
    	spectrum->AddBinContent(k,-1.*offset);
   	   }
    }
*/
		// check if spectrum exists, if yes delete so that we do not get mem leaks

	// First copy spectrum into new histo, so we do not overwrite stuff

	/*
	Int_t nchan = spectrum1->GetNbinsX();
 	TH1D *spectrum = new TH1D("spectrum","Spectrum for background subtraction ",nchan,spectrum1->GetBinCenter(0),spectrum1->GetBinCenter(nchan-1));
 	// delete spectrum and recreate it so no memery leak; stupid way

 	spectrum->Sumw2();
 	// Now fille the new histo with the signal noq histo gram
 	// first get number of channels in NMR1_NoQ

 	for (Int_t k=0 ; k < spectrum1->GetNbinsX(); k++){

 		spectrum->Fill(spectrum1->GetBinCenter(k),spectrum1->GetBinContent(k));
 	}
*/


	 //ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Fumili2"); // faster minimizer

		ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2","Minimize");
		ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);
		//ROOT::Math::MinimizerOptions::SetDefaultTolerance(.0001);
		ROOT::Math::MinimizerOptions::SetDefaultTolerance(.001);
		ROOT::Math::MinimizerOptions::SetDefaultPrecision(1.e-10);
		ROOT::Math::MinimizerOptions::SetDefaultMaxIterations(100000);
		ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1000000);
		ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(0);

//		extern int gErrorIgnoreLevel;
//		gErrorIgnoreLevel = 1001;



	makeTF1();
	FitBck->SetNpx(1000);
	FitBck->SetLineWidth(4);
	FitBck->SetLineColor(kYellow);

// first we try to get starting pararemets by fitting the background only
	// to the lower background area
	reject1 = fit_limit[1];
	reject2 = fit_limit[2];

	// the upper limit of the lower background window we choose
	spectrum->Fit(FitBck,"RE0Q");
	FitBck->GetParameters(&bck_par[0]);
	// now we fit the whole spectrum with the new variables
/*	reject1 = fit_limit[1];
	reject2 = fit_limit[2];
	spectrum->Fit(FitBck,"RE0Q"); // parameter Q for quiet
	FitBck->GetParameters(&bck_par[0]);
*/

	// Create new function to subtract from spectrum
	// need to do this since otherwise it only subtracts in the range defined by the fit range
	FitBckCopy = new TF1("FitBckCopy",CopyBackground,fit_limit[0],fit_limit[3],3);
	FitBckCopy->SetParameters(bck_par);
// 	now subtract the background
	spectrum->Add(FitBckCopy,-1.);


	//TF1 *fhelp = new TF1("fhelp","[0]+[1]*x+[2]*x*x",fit_limit[0],fit_limit[3]);
//	TF1 *BckFct = new TF1("BckFct","[0]+[1]*x+[2]*x*x + [3]*pow(x,3)",	spectrum->GetXaxis()->GetXmin(),	spectrum->GetXaxis()->GetXmax());
//	BckFct->SetParameters(bck_par);

//    BckFct1 = (TF1*)BckFct->Clone("BckFct1");


    return spectrum;
}

Double_t Ana::BackSpline(TGraph *gr1){
	TCanvas * chelp = new TCanvas();
	gr1->Draw();
	chelp->Update();
	return 1.;

}
void  Ana::SetFitLimits(Double_t x1, Double_t x2 , Double_t x3, Double_t x4){
	// this sets the lower and upper window of the quadratic background fit for the background subtraction
	fit_limit[0] = x1;
	fit_limit[1] = x2;
	fit_limit[2] = x3;
	fit_limit[3] = x4;
	}

#endif /* ANA_H_ */
