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
	Int_t FitSpectrum(TH1D *, Int_t  );

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
	  return (0.5*par[0]*par[1]/TMath::Pi()) /
	    TMath::Max( 1.e-10,(x[0]-par[2])*(x[0]-par[2])
	   + .25*par[1]*par[1])*TMath::Exp(-.5*pow((x[3]-par[4])/par[5],2));


}

int Ana::FitSpectrum(TH1D * Spectrum,Int_t NumberOfSpectra){
	// function to fit spectrum with this needs to be called after Findpeak
	// since it relies on finding the parameters.




	fit_low_overall = *xpeaks-.1;
	fit_high_overall = *xpeaks +.1;
	FitFcn =  new TF1("FitFcn",CombinedFit,fit_low_overall,fit_high_overall,6);

	FitFcn->SetNpx(1000);
	FitFcn->SetLineWidth(4);
	FitFcn->SetLineColor(kMagenta);
/*	if(FitF.Contains("gaus")){

	FitFcn->SetParameters(1,1,1,PeakFitPar[0],PeakFitPar[1],PeakFitPar[2]);
	FitFcn->FixParameter(4,PeakFitPar[1]);
	FitFcn->FixParameter(5,PeakFitPar[2]);
	}
	else {
		FitFcn->SetParameters(1,1,1,PeakFitPar[0],PeakFitPar[1],PeakFitPar[2]);
		FitFcn->FixParameter(4,PeakFitPar[1]);
		FitFcn->FixParameter(5,PeakFitPar[2]);
	}
*/


	Spectrum->Fit("FitFcn","+REM","C");
	//Spectrum->Draw();
	//FitFcn->Draw("SAME");

	// Now get the single fits for background and signal

/*	   // improve the picture:
	   BackGround = new TF1("BackGround",Background,fit_low_overall,fit_high_overall,3);
	   BackGround->SetLineColor(kRed);
		if(FitF.Contains("gaus"))
		{
			Signal = new TF1("Signal",GausPeak,fit_low_overall,fit_high_overall,3);
		}
		else
		{
			Signal = new TF1("Signal",LorentzianPeak,fit_low_overall,fit_high_overall,3);

		}
	   Signal->SetLineColor(kBlue);
	   Signal->SetNpx(500);
	   Double_t par[6];
	   // writes the fit results into the par array
	   FitFcn->GetParameters(par);
	   BackGround->SetParameters(par);
	   BackGround->Draw("same");
	   Float_t low =0.;
	   Float_t high = 0.;
	   if(FitF.Contains("gaus")){

			Signal->SetParameters(&par[3]);
			Signal->Draw("same");
			 low = Signal->GetParameter(1)-3*(Signal->GetParameter(2));
			 high = Signal->GetParameter(1)+3*(Signal->GetParameter(2));
		}
		else
		{
			   Signal->SetParameters(&par[3]);
			   Signal->Draw("same");
			    low = Signal->GetParameter(2)-1.2*(Signal->GetParameter(1));
			    high = Signal->GetParameter(2)+1.2*(Signal->GetParameter(1));
		}

	   cout<< "\n\n\n $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ \n\n\n";
	   cout<<"low integration limit  "<<low<<"   high integration limit  "<<high<<"   \n";
	   NMRArea = Signal->Integral(low,high);
       cout<<"integral of fit  corrected histogram  "<<NMRArea/NumberOfSpectra<<" \n";
	   cout<< "\n\n\n $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ \n\n\n";
	   //fprintf(fp_out,"%lf \n", NMRArea/NumberOfSpectra);

	   // Now integrate the whole peak to get the NMR value
	   //
*/
	   // draw the legend
	   TLegend *legend=new TLegend(0.6,0.65,0.88,0.85);
	   legend->SetTextFont(72);
	   legend->SetTextSize(0.04);
	   legend->AddEntry(Spectrum,"NMR corrected","lpe");
//	   legend->AddEntry(BackGround,"Background fit","l");
//	   legend->AddEntry(Signal,"Signal fit","l");
	   legend->AddEntry(FitFcn,"Global Fit","l");
	   legend->Draw();







	return 1;
}




#endif /* ANA_H_ */
