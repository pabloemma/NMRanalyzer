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

using namespace std;

class Ana
{
private:
	Bool_t peakfind;
	Double_t  fit_low_overall; // fit limits
	Double_t fit_high_overall;

public:
	std::string Ana_pr ="NMR_ana> ";


	TH1D *Spectrum;  // the histogram we want to find the spectrum of
	TSpectrum *spec;
	TF1 *FitFcn;  // convoluted fit function
	TF1 *FitGaus; // first Gauss fit
	Int_t npeaks;	//number of peaks to find

	Float_t *xpeaks; // xpositions of peaks
	Float_t sigma; // sigma for peak search
	Double_t amplitude;// amplitude from peak serach
	Float_t Offset; // Offset which is added so that we do not have any negative counts
					// this offset gets determned from the loweest histogram entry and then added if that one is negative
	Int_t NumberOfSpectra; // how many spectra are in the NMR sweep.


	Double_t par[6]; // the fit parameters
	Double_t gaus_par[3];// for initial gaus fit






	Ana();
	virtual ~Ana();
	virtual void FindPeak(TH1D* );
	virtual void FindOffset(TH1D* );

	static Double_t CombinedFit(Double_t *, Double_t *); // Gaus folded with Lorentzian distribution
	Int_t FitSpectrum(TH1D *, Int_t  );
	static Double_t GausPeak(Double_t *, Double_t *);


}; //end of class definition

Ana::Ana(){
	//initialize class
	cout<< Ana_pr<<"************* initializing analyzer  ************** \n";
}

Ana::~Ana(){
	cout<<Ana_pr<<"Done with analysis \n";
}
void Ana::FindPeak(TH1D * Spectrum){
	// this is from peaks.C at Cern
	// Create new spectrum,we only assume one peak for the moment
	//


	peakfind = true;
	spec= new TSpectrum(npeaks,1.);

	Int_t nfound = spec->Search(Spectrum,sigma,"nobackground new",.01);
	cout<<"\n\n\n*******************************************\n";
	cout<<Ana_pr<<"Number of peaks found "<<nfound<<"\n";
	// fill array with peak posistions

	xpeaks = (spec->GetPositionX());

	if(*xpeaks < 212.982 - .2 || *xpeaks > 212.982 + .2 ){
		cout<< "TSpectrum failed, assign peak to 212.98 \n";
		*xpeaks = 212.982;
		sigma = .018;
		peakfind = false;
	}
	else
	{
		cout<<Ana_pr<<" Tspectrum found peak at :  "<<*xpeaks<<" \n\n";
	}


}

void Ana::FindOffset(TH1D *Spectrum ){
	//this routine determines if there is a negative channel and moves this up by that amount
	// so that a fittin would still be okay
	// Determine if any of the cnannels is negative and calculate the offset
	Offset = Spectrum->GetBinContent(Spectrum->GetMinimumBin());

	if (Offset < 0.){
		cout<<"Add offset to spectrum  "<< Offset<<"   "<<Spectrum->GetNbinsX()<<"\n";
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



int Ana::FitSpectrum(TH1D * Spectrum,Int_t NumberOfSpectra){
	// function to fit spectrum with this needs to be called after Findpeak
	// since it relies on finding the parameters.




	fit_low_overall = *xpeaks-.2;
	fit_high_overall = *xpeaks +.2;
	FitFcn =  new TF1("FitFcn",CombinedFit,fit_low_overall,fit_high_overall,6);

	FitFcn->SetNpx(1000);
	FitFcn->SetLineWidth(4);
	FitFcn->SetLineColor(kMagenta);

	FitGaus = new TF1("FitGaus",GausPeak,fit_low_overall,fit_high_overall,3);
    FitGaus->SetParameters(1.,213.,1.);



	Spectrum->Fit("FitGaus","RE","C");


	FitFcn->SetParameters(1,1,1,FitGaus->GetParameter(0),FitGaus->GetParameter(1),FitGaus->GetParameter(2));




	Spectrum->Fit("FitFcn","MRE+","C");
	//Spectrum->Draw();
	FitFcn->Draw("SAME");






	   // draw the legend
	   TLegend *legend=new TLegend(0.6,0.65,0.88,0.85);
	   legend->SetTextFont(72);
	   legend->SetTextSize(0.04);
	   legend->AddEntry(Spectrum,"NMR corrected","lpe");
//	   legend->AddEntry(BackGround,"Background fit","l");
	   legend->AddEntry(FitGaus,"Gaus fit","l");
	   legend->AddEntry(FitFcn,"Global Fit","l");
	   legend->Draw();







	return 1;
}




#endif /* ANA_H_ */
