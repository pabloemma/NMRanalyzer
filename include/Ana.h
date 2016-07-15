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

using namespace std;

class Ana
{
private:
	Bool_t peakfind;
	 Double_t  fit_low_overall; // fit limits
	 Double_t fit_high_overall;
     Double_t fit_limit[4]; ; // these are the limits for the background fit
public:
	std::string Ana_pr ="Ana_ana> ";


	TH1D *Spectrum;  // the histogram we want to find the spectrum of
	TSpectrum *spec;
	TF1 *FitFcn;  // convoluted fit function
	TF1 *FitGaus; // first Gauss fit
	TF1 *FitBck; // background fit
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
	Double_t bck_par[3]; // background polynomial






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
	 void FitBackground(TH1D*);
	 void SetFitLimits(Double_t, Double_t, Double_t, Double_t);
	 Double_t Background(Double_t *, Double_t *);
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
	Int_t nfound = spec->Search(Spectrum,sigma,"nobackground new",.01);
//	cout<<"\n\n\n*******************************************\n";
//	cout<<Ana_pr<<"Number of peaks found "<<nfound<<"\n";
	// fill array with peak posistions

	xpeaks = (spec->GetPositionX());

	fit_low_overall = *xpeaks-.2;
	fit_high_overall = *xpeaks +.2;


	if(*xpeaks < 212.982 - .3 || *xpeaks > 212.982 + .3 ){
		cout<<Ana_pr<< "TSpectrum failed, assign peak to 212.98 \n";
		*xpeaks = 212.982;
		sigma = .018;
		peakfind = false;
	}

	cout<<"peak_find \n";

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
//	if(x[0]> 212.85  && x[0] < 213.26){
		if(x[0]>fit_limit[1]  && x[0] < fit_limit[2]){
	     TF1::RejectPoint();
	     //return 0;
	 }

   return (par[0] + par[1]*x[0] + par[2]*x[0]*x[0]);
   //return par[0] + par[1]*x[0] ;

}
// Quadratic background function
/*Double_t Ana::Background2(Double_t low_x, Double_t high_x,Double_t *x, Double_t *par) {
	// new version with point rejection, the idea being that
	// I will not fit background in the peak area but on left and right side of it
	// see fit descrition in Root reference manual
	if(x[0]> low_x  && x[0] < high_x){
	     TF1::RejectPoint();
	 }

   return (par[0] + par[1]*x[0] + par[2]*x[0]*x[0]);
   //return par[0] + par[1]*x[0] ;

}
*/

int Ana::FitSpectrum(TH1D * Spectrum,Int_t NumberOfSpectra){
	// function to fit spectrum with this needs to be called after Findpeak
	// since it relies on finding the parameters.



	   // correct for negative
		FindOffset(Spectrum);
		FitBackground(Spectrum);



	FitFcn =  new TF1("FitFcn",CombinedFit,fit_low_overall,fit_high_overall,6);

	FitFcn->SetNpx(1000);
	FitFcn->SetLineWidth(4);
	FitFcn->SetLineColor(kMagenta);

	FitGaus = new TF1("FitGaus",GausPeak,fit_low_overall,fit_high_overall,3);
    FitGaus->SetParameters(1.,213.,1.);



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
	FitBck = new TF1("FitBck",this,&Ana::Background,fit_limit[0],fit_limit[3],3);
}






void Ana::FitBackground(TH1D *spectrum){
	// this just determines the backgtound parameters of the spectrum
	// Currently it is a simple 2degree polynomial

	makeTF1();
	FitBck->SetNpx(1000);
	FitBck->SetLineWidth(4);
	FitBck->SetLineColor(kYellow);

	FitBck->SetParameters(1.,1.,1.); // initialze parameters

	spectrum->Fit(FitBck,"RQM");
	FitBck->GetParameters(bck_par);
	// Create new function to subtract from spectrum
	// need to do this since otherwise it only subtracts in the range defined by the fit range

	//TF1 *fhelp = new TF1("fhelp","[0]+[1]*x+[2]*x*x",fit_limit[0],fit_limit[3]);
	TF1 *BckFct = new TF1("BckFct","[0]+[1]*x+[2]*x*x",	spectrum->GetXaxis()->GetXmin(),	spectrum->GetXaxis()->GetXmax());
	BckFct->SetParameters(bck_par);
    BckFct1 = (TF1*)BckFct->Clone("BckFct1");
/*	TCanvas *cc1 = new TCanvas();
	cc1->cd();
	spectrum->Draw();
	fhelp->Draw("SAME");
    spectrum->Add(fhelp,-1.);
    spectrum->Draw("SAME");
    cc1->Update(); */


    // now check for offset below 0 and then shift this offset up
    Double_t offset = spectrum->GetBinContent(spectrum->GetMinimumBin());
    if(offset<0.){
    for(Int_t k=0;k<spectrum->GetNbinsX();k++){
    	spectrum->AddBinContent(k,-1.*offset);
   	   }
    }

    //return 1;
}

Double_t Ana::BackSpline(TGraph *gr1){
	TCanvas * chelp = new TCanvas();
	gr1->Draw();
	chelp->Update();
	return 1.;

}
void  Ana::SetFitLimits(Double_t x1, Double_t x2 , Double_t x3, Double_t x4){
	// this sets the lower and upper window of the quadratic background fit for the backgriound subtraction
	fit_limit[0] = x1;
	fit_limit[1] = x2;
	fit_limit[2] = x3;
	fit_limit[3] = x4;
	}

#endif /* ANA_H_ */
