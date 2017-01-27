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
 virtual void FindPeak(TH1D* );
 virtual void FindOffset(TH1D* );

	 static Double_t CombinedFit(Double_t *, Double_t *); // Gaus folded with Lorentzian distribution
		Int_t FitSpectrum(TH1D *, Int_t  );

	 static Double_t GausPeak(Double_t *, Double_t *);
	 static Double_t Background2(Double_t *, Double_t *,Double_t *, Double_t *);
//	 static Double_t BackSpline(Double_t *, Double_t *,Double_t *, Double_t *, TGraph *); // This will determine a spline back ground
	 static Double_t BackSpline( TGraph *);
	 TH1D* FitBackground(TH1D*); // this one gets called at every sweep event
	 TH1D* FitBackground1(TH1D*); // this one gets called at the end by the signal histo
	 void SubtractLinear(TH1D*,Int_t,Int_t,Int_t,Int_t); //
	 void SetFitLimits(Double_t, Double_t, Double_t, Double_t);
	 Double_t Background(Double_t *, Double_t *);
	 static Double_t CopyBackground(Double_t *, Double_t *);
	 void makeTF1();



}; //end of class definition

Ana::Ana(){
	//initialize class
	cout<< Ana_pr<<"************* initializing analyzer  ************** \n";
}

Ana::~Ana(){
	cout<<Ana_pr<<"Done with analysis \n\n\n";
}
void Ana::FindPeak(TH1D * Spectrum){
	// this is from peaks.C at Cern
	// Create new spectrum,we only assume one peak for the moment
	//

	npeaks = 3;
	peakfind = true;
	spec= new TSpectrum(npeaks,1.);
	Int_t nfound = spec->Search(Spectrum,sigma=.02,"nobackground new",.5);
	//cout<<"\n\n\n*******************************************\n";
	//cout<<Ana_pr<<"Number of peaks found "<<nfound<<"\n";
	// fill array with peak posistions

	xpeaks = (Double_t *)(spec->GetPositionX());

	//fit_low_overall = *xpeaks -.2;;
	//fit_high_overall = *xpeaks +.2;
	//fit_low_overall = *xpeaks -.2;;
	//fit_high_overall = *xpeaks +.2;


	if(*xpeaks < fit_limit[1] - .1 || *xpeaks > fit_limit[2] + .1 ){
		cout<<Ana_pr<< "TSpectrum failed, assign peak to 212.98 \n";
		*xpeaks = 212.982;
		sigma = .018;
		peakfind = false;
	}


}

void Ana::FindOffset(TH1D *Spectrum ){
	//this routine determines if there is a negative channel and moves this up by that amount
	// so that a fittin would still be okay
	// Determine if any of the cnannels is negative and calculate the offset
	Offset = Spectrum->GetBinContent(Spectrum->GetMinimumBin());

	if (Offset < 0.){
		cout<<Ana_pr<<"Add offset to spectrum  "<< Offset<<"   "<<Spectrum->GetNbinsX()<<"\n";
		for (Int_t k = 0;  k< Spectrum->GetNbinsX();k++){
			//cout<<"adding bin"<<k<<"\n";
			Spectrum->AddBinContent(k,-Offset);  // add offset
		}
	}
}

Double_t Ana::CombinedFit(Double_t *x, Double_t *par){
	// this function is a convolution of Lorentzian and Gaus
	  Double_t Lo = (0.5*par[0]*par[1]/TMath::Pi()) /
	    TMath::Max( 1.e-10,(x[0]-par[2])*(x[0]-par[2])
	   + .25*par[1]*par[1]);
		Double_t arg = 0;
		if(par[5] !=0) arg = (x[0]-par[4])/par[5];
	  Double_t Ga = par[3]*TMath::Exp(-.5*arg*arg);

      return Lo*Ga;
}

Double_t Ana::GausPeak(Double_t *x, Double_t *par){
	Double_t arg = 0;
	if(par[2] !=0) arg = (x[0]-par[1])/par[2];
	return par[0]*TMath::Exp(-.5*arg*arg);
}

// Quadratic background function
Double_t Ana::Background(Double_t *x, Double_t *par) {
	// new version with point rejection, the idea being that
	// I will not fit background in the peak area but on left and right side of it
	// see fit descrition in Root reference manual
	// 3rd degree polynomial

		if(x[0]>reject1  && x[0] < reject2){
	     TF1::RejectPoint();
	     //return 0;
	 }

   return (par[0] + par[1]*x[0] + par[2]*pow(x[0],2)+par[3]*pow(x[0],3));
;

}
Double_t Ana::CopyBackground(Double_t *x, Double_t *par) {
	// new version with point rejection, the idea being that
	// I will not fit background in the peak area but on left and right side of it
	// see fit descrition in Root reference manual
	// 3rd degree polynomial, this is just a copy of the real backgound but with no rejection
	// this way I can subtract the function from the histogram.



   return (par[0] + par[1]*x[0] + par[2]*pow(x[0],2)+par[3]*pow(x[0],3));
;

}

int Ana::FitSpectrum(TH1D * Spectrum,Int_t NumberOfSpectra){
	// function to fit spectrum with this needs to be called after Findpeak
	// since it relies on finding the parameters.



	   // correct for negative
//		FindOffset(Spectrum);
		FitBackground(Spectrum);

    fit_low_overall = fit_limit[1];
    fit_high_overall = fit_limit[2];

	FitFcn =  new TF1("FitFcn",CombinedFit,fit_low_overall,fit_high_overall,6);

	FitFcn->SetNpx(1000);
	FitFcn->SetLineWidth(4);
	FitFcn->SetLineColor(kMagenta);

	FitGaus = new TF1("FitGaus",GausPeak,fit_low_overall,fit_high_overall,3);
    FitGaus->SetParameters(1.,xpeaks[0],1.);



	Spectrum->Fit("FitGaus","RE","C");

	// lets fit the background
	//FitBackground(Spectrum);


	FitFcn->SetParameters(1,1,1,FitGaus->GetParameter(0),FitGaus->GetParameter(1),FitGaus->GetParameter(2));




	Spectrum->Fit("FitFcn","MRE+","C");
	//Spectrum->Draw();
	FitFcn->Draw("SAME");
	FitBck->Draw("SAME");
	BckFct1->Draw("SAME");






	   // draw the legend
	   TLegend *legend=new TLegend(0.6,0.65,0.88,0.85);
	   legend->SetTextFont(72);
	   legend->SetTextSize(0.04);
	   legend->AddEntry(Spectrum,"NMR corrected","lpe");
	   legend->AddEntry(FitBck,"Background fit","l");
	   legend->AddEntry(FitGaus,"Gaus fit","l");
	   legend->AddEntry(FitFcn,"Global Fit","l");
	   legend->Draw();







	return 1;
}

void Ana::makeTF1(){
	// this is so that I can pass parameters to the background fit
	FitBck = new TF1("FitBck",this,&Ana::Background,fit_limit[0],fit_limit[3],4);
}






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
		ROOT::Math::MinimizerOptions::SetDefaultTolerance(.0001);
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
	reject1 = fit_limit[0];
	reject2 = fit_limit[1];

	// the upper limit of the lower background window we choose
	spectrum->Fit(FitBck,"RE0Q");
	FitBck->GetParameters(&bck_par[0]);
	// now we fit the whole spectrum with the new variables
	reject1 = fit_limit[1];
	reject2 = fit_limit[2];
	spectrum->Fit(FitBck,"RE0Q"); // parameter Q for quiet
	FitBck->GetParameters(&bck_par[0]);


	// Create new function to subtract from spectrum
	// need to do this since otherwise it only subtracts in the range defined by the fit range
	FitBckCopy = new TF1("FitBckCopy",CopyBackground,fit_limit[0],fit_limit[3],4);
	FitBckCopy->SetParameters(bck_par);
// 	now subtract the background
	spectrum->Add(FitBckCopy,-1.);


	//TF1 *fhelp = new TF1("fhelp","[0]+[1]*x+[2]*x*x",fit_limit[0],fit_limit[3]);
//	TF1 *BckFct = new TF1("BckFct","[0]+[1]*x+[2]*x*x + [3]*pow(x,3)",	spectrum->GetXaxis()->GetXmin(),	spectrum->GetXaxis()->GetXmax());
//	BckFct->SetParameters(bck_par);

//    BckFct1 = (TF1*)BckFct->Clone("BckFct1");


    return spectrum;
}
TH1D* Ana::FitBackground1(TH1D *spectrum1){
		//TH1D* Ana::FitBackground(TH1D *spectrum){
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


		Int_t nchan = spectrum1->GetNbinsX();
	 	TH1D *spectrum = new TH1D("spectrum","Spectrum for background subtraction ",nchan,spectrum1->GetBinCenter(0),spectrum1->GetBinCenter(nchan-1));
	 	// delete spectrum and recreate it so no memery leak; stupid way

	 	spectrum->Sumw2();
	 	// Now fille the new histo with the signal noq histo gram
	 	// first get number of channels in NMR1_NoQ

	 	for (Int_t k=0 ; k < spectrum1->GetNbinsX(); k++){

	 		spectrum->Fill(spectrum1->GetBinCenter(k),spectrum1->GetBinContent(k));
	 	}



		 ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Fumili2"); // faster minimizer
		makeTF1();
		FitBck->SetNpx(1000);
		FitBck->SetLineWidth(4);
		FitBck->SetLineColor(kYellow);

	// first we try to get starting pararemets by fitting the background only
		// to the lower background area
		reject1 = fit_limit[0];
		reject2 = fit_limit[1];

		// the upper limit of the lower background window we choose
		spectrum->Fit(FitBck,"RE0Q");
		FitBck->GetParameters(&bck_par[0]);
		// now we fit the whole spectrum with the new variables
		reject1 = fit_limit[1];
		reject2 = fit_limit[2];
		spectrum->Fit(FitBck,"RE0Q"); // parameter Q for quiet
		FitBck->GetParameters(&bck_par[0]);


		// Create new function to subtract from spectrum
		// need to do this since otherwise it only subtracts in the range defined by the fit range
		FitBckCopy = new TF1("FitBckCopy",CopyBackground,fit_limit[0],fit_limit[3],4);
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
void Ana::SubtractLinear(TH1D *spec,Int_t Ifit_x1,Int_t Ifit_x2,Int_t Ifit_x3,Int_t Ifit_x4){
	// determines a linear fit on the QCurve subtracted histogram
	// first we take the integral on the low x and the high x
	// take lower range and integrate and higher range
	// mostly pseudocode right now



    Double_t ylow = spec->Integral(Ifit_x1,Ifit_x2)/(Ifit_x2-Ifit_x1+1);
	Double_t yhigh = spec->Integral(Ifit_x3,Ifit_x4)/(Ifit_x4-Ifit_x3+1);
	Double_t xlow = (spec->GetBinCenter(Ifit_x2)+spec->GetBinCenter(Ifit_x1))/(Ifit_x2-Ifit_x1+1);
	Double_t xhigh = (spec->GetBinCenter(Ifit_x4)+spec->GetBinCenter(Ifit_x3))/(Ifit_x4-Ifit_x3+1);



	Double_t slope = (yhigh-ylow)/(xhigh-xlow);
	Double_t offset = ylow - slope*xlow;

    for (Int_t k =0; k < spec->GetNbinsX(); k++){
    	Double_t newval = offset +slope*spec->GetBinCenter(k);
      	Double_t temp = spec->GetBinContent(k)-newval;
    	spec->SetBinContent(k,temp);
    }


}


#endif /* ANA_H_ */
