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
     Double_t fit_limit[4]; ; // these are the limits for the background fit
     Double_t reject0,reject1, reject2,reject3;// the lower and upper limits of what backrgound fit rejects.
public:
	std::string Ana_pr ="Ana_ana> ";


	TH1D *Spectrum;  // the histogram we want to find the spectrum of
	TF1 *FitBck; // background fit
	TF1 *ff1, *ff2,*ff3,*ff4;


	Int_t npeaks;	//number of peaks to find


	Double_t par[3] , bck_par[3]; // the fit parameters






	Ana();
	virtual ~Ana();


//	 static Double_t BackSpline(Double_t *, Double_t *,Double_t *, Double_t *, TGraph *); // This will determine a spline back ground
	 void SetFitLimits(Double_t, Double_t, Double_t, Double_t);
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




void Ana::makeTF1(){
	// this is so that I can pass parameters to the background fit
	//FitBck = new TF1("FitBck",this,&Ana::Background,fit_limit[0],fit_limit[3],4);
	//FitBck = new TF1("FitBck",this,&Ana::Background,fit_limit[0],fit_limit[3],3);
	FitBck = new TF1("FitBck", "pol2", fit_limit[0], fit_limit[3]);
}
void Ana::makeTF11(){
	// this is so that I can pass parameters to the background fit
	//ff1 = new TF1("ff1",this,&Ana::Background,fit_limit[0],fit_limit[3],4);
	//ff1 = new TF1("ff1",this,&Ana::Background_1,fit_limit[0],fit_limit[3],3);
	ff1 = new TF1("pol2", "pol2", fit_limit[0], fit_limit[3]);

}


TF1 *Ana::NewFitBackground(TH1D *spectrum){

// this is using a TGraph to fit the spectrum
// first we create a TGraph from a histo by filling x, y vector
// then we fit the background
// and finally we return this graph

	// get info on histo
    Double_t dataAmp[1000] , freq[1000];
	Int_t nchan = spectrum->GetNbinsX();
	reject0 = fit_limit[0];
	reject1 = fit_limit[1];
	reject2 = fit_limit[2];
	reject3 = fit_limit[3];
    Int_t counter =0;
	for(Int_t k=1; k<=nchan;k++){
		Double_t freq_temp = (spectrum->GetBinCenter(k));
		if( (freq_temp > reject0 && freq_temp < reject1) || (freq_temp >reject2 && freq_temp < reject3) ){
		     freq[counter]=freq_temp;
		     dataAmp[counter]=(spectrum->GetBinContent(k));
		     counter++;
		}
	}
	//Now create Graph
    makeTF11();
	TGraph *gr1 = new TGraph(counter,freq,dataAmp);
//

	//Now comes the fit part

	ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2","Minimize");
	ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);
	//ROOT::Math::MinimizerOptions::SetDefaultTolerance(.0001);
	ROOT::Math::MinimizerOptions::SetDefaultTolerance(.0001);
	ROOT::Math::MinimizerOptions::SetDefaultPrecision(1.e-10);
	ROOT::Math::MinimizerOptions::SetDefaultMaxIterations(100000);
	ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1000000);
	//ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(0);



	//gr1->Fit(ff1,"");
	gr1->Fit(ff1,"Q");


	ff1->GetParameters(&bck_par[0]);
	// Create new function to subtract from spectrum
	// need to do this since otherwise it only subtracts in the range defined by the fit range
	ff2 = new TF1("ff2", "pol2", fit_limit[0], fit_limit[3]);
	ff2->SetParameters(bck_par);

	// Now do background subtraction
	spectrum->Add(ff2,-1.);

	return ff2;
}


void  Ana::SetFitLimits(Double_t x1, Double_t x2 , Double_t x3, Double_t x4){
	// this sets the lower and upper window of the quadratic background fit for the background subtraction
	fit_limit[0] = x1;
	fit_limit[1] = x2;
	fit_limit[2] = x3;
	fit_limit[3] = x4;
	}

#endif /* ANA_H_ */
