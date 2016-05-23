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

public:
	std::string Ana_pr ="NMR_ana> ";


	TH1D *Spectrum;  // the histogram we want to find the spectrum of
	TSpectrum *spec;
	Int_t npeaks;	//number of peaks to find

	Float_t *xpeaks; // xpositions of peaks
	Float_t sigma; // sigma for peak search
	Double_t amplitude;// amplitude from peak serach
	Float_t Offset; // Offset which is added so that we do not have any negative counts
					// this offset gets determned from the loweest histogram entry and then added if that one is negative
	Int_t NumberOfSpectra; // how many spectra are in the NMR sweep.


	Double_t par[6]; // the fit parameters






	Ana();
	virtual ~Ana();
	virtual void FindPeak(TH1D* );
	virtual void FindOffset(TH1D* );

	static Double_t CombinedFit(Double_t *, Double_t *); // Gaus folded with Lorentzian distribution


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


Double_t Ana::CombinedFit(Double_t *x, Double_t *par){
	// this function is a convolution of Lorentzian and Gaus
	  return (0.5*par[0]*par[1]/TMath::Pi()) /
	    TMath::Max( 1.e-10,(x[0]-par[2])*(x[0]-par[2])
	   + .25*par[1]*par[1])*TMath::Exp(-.5*pow((x[3]-par[4])/par[5],2));


}


}


#endif /* ANA_H_ */
