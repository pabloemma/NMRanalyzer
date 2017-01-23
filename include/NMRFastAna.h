#ifndef _NMRFastAna_H
#define _NMRFastAna_H

#include <iostream>
#include <vector>
#include <TH1D.h>
#include <TString.h>
#include <TGraph.h>
#include <TCanvas.h>

#include <Math/Factory.h>
#include <Math/Minimizer.h>
#include <Math/Functor.h>

class NMRFastAna
{
public:
	NMRFastAna();
	~NMRFastAna();

	//!Initialization
	void init();

	//!Set fit frequency range
	void setFreqRange(double min, double max) { freqMin = min; freqMax = max; } 
	void setExclusionWin(double center, double win) { freqCenter = center; freqWin = win; }
	void setSampleRate(int rate) { sampleRate = rate; }

	//!Set fit parameter range
	void setXOffsetRange(int min, int max) { xoffsetMin = min; xoffsetMax = max; }
	void setYOffsetRange(double min, double max) { yoffsetMin = min; yoffsetMax = max; }

	//!Set the Q curve
	void setQCurve(TH1D* ref);

	//!main external call
	int getXOffset(TH1D* data);

	//!additional getter
	int getXOffset() { return xoffset; }
	double getYOffset() { return yoffset; }

	//!make a plot of fitting results
	void plot(TString name);

	//!chisq calcuation -- internal use
	double chisq(const double* par);

private:
	//!Minimizer and functor
	//@{
	ROOT::Math::Minimizer* minimizer;
	ROOT::Math::Functor fcn;
	//@}

	//!Fit configurations
	//@{
	double freqMin;
	double freqMax;
	double freqCenter;
	double freqWin;
	int sampleRate;
	//@}

	//!parameter ranges
	//@{
	int xoffsetMin;
	int xoffsetMax;
	double yoffsetMin;
	double yoffsetMax;
	//@}

	//!internal array to store data
	//@{
	unsigned int nScanPoints;
	std::vector<double> freq;
	std::vector<double> dataAmp;
	std::vector<double> refAmp;
	//@}

	//fit parameters
	//@{
	int xoffset;    //number of bins to shift
	double yoffset; //offset in y
	//@}
};

NMRFastAna::NMRFastAna()
{
	//default configurations
	freqCenter = 213.;
	freqWin = 0.1;
	freqMin = 213. - 0.5;
	freqMax = 213. + 0.5;
	sampleRate = 4;

	xoffsetMin = -50;
	xoffsetMax = 50;
	yoffsetMin = -0.5;
	yoffsetMax = 0.5;

	//set the pointers to NULL
	minimizer = NULL;
}

NMRFastAna::~NMRFastAna()
{
	if(minimizer != NULL) delete minimizer;
}

void NMRFastAna::init()
{
	fcn = ROOT::Math::Functor(this, &NMRFastAna::chisq, 2);

	//Initialize Fitter
	minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Combined");
	minimizer->SetMaxFunctionCalls(1000000);
	minimizer->SetMaxIterations(10000);
	minimizer->SetTolerance(0.0001);
	minimizer->SetFunction(fcn);
	minimizer->SetPrintLevel(0);
}

double NMRFastAna::chisq(const double* par)
{
	xoffset = int(par[0]);
	yoffset = par[1];

	double chisq = 0.;
	unsigned int refIndex = xoffset > 0 ? 0 : abs(xoffset);
	unsigned int datIndex = xoffset > 0 ? abs(xoffset) : 0;
	while(refIndex < nScanPoints && datIndex < nScanPoints)
	{
		if(freq[datIndex] > freqMin && freq[datIndex] < freqMax && fabs(freq[datIndex] - freqCenter) > freqWin) 
			chisq = chisq + (dataAmp[datIndex] - yoffset - refAmp[refIndex])*(dataAmp[datIndex] - yoffset - refAmp[refIndex]);
		
		refIndex += sampleRate;
		datIndex += sampleRate;
	}

	//std::cout << xoffset << "  " << yoffset << " : " << chisq << std::endl;
	return chisq;
}

void NMRFastAna::setQCurve(TH1D* ref)
{
	//Load the data from hitstogram to a vector
	nScanPoints = ref->GetNbinsX();

	freq.clear();
	refAmp.clear();
	for(unsigned i = 1; i <= nScanPoints; ++i)
	{
		freq.push_back(ref->GetBinCenter(i));
		refAmp.push_back(ref->GetBinContent(i));
	}
}

int NMRFastAna::getXOffset(TH1D* data)
{
	//initialize parameters
	xoffset = 0;
	yoffset = 0.;

	dataAmp.clear();
	for(unsigned i = 1; i <= nScanPoints; ++i)
	{
		dataAmp.push_back(data->GetBinContent(i));
	}

	//Set initial value and range
	minimizer->SetLimitedVariable(0, "xoffset", 0., 1.,    xoffsetMin, xoffsetMax);
	minimizer->SetLimitedVariable(1, "yoffset", 0., 0.001, yoffsetMin, yoffsetMax);

	minimizer->Minimize();

	//Get the result
	xoffset = int(minimizer->X()[0]);
	yoffset = minimizer->X()[1];

	return xoffset;
}

void NMRFastAna::plot(TString name)
{
	const int nPointsMax = 1000;
	double x[nPointsMax], yBkg[nPointsMax], yAll[nPointsMax];

	int nPointsLeft = 0;
	int refIndex = xoffset > 0 ? 0 : abs(xoffset);
	int datIndex = xoffset > 0 ? abs(xoffset) : 0;
	while(refIndex < nScanPoints && datIndex < nScanPoints)
	{
		//std::cout << refIndex << "  " << ampRef[refIndex] << std::endl;
		x[nPointsLeft] = freq[datIndex];
		yBkg[nPointsLeft] = yoffset + refAmp[refIndex];
		yAll[nPointsLeft] = dataAmp[datIndex];
		
		++refIndex;
		++datIndex;
		++nPointsLeft;
	}

	TCanvas c1;
	TGraph gBkg(nPointsLeft, x, yBkg);
	TGraph gAll(nPointsLeft, x, yAll);

	gBkg.SetLineColor(kBlue);
	gBkg.SetLineWidth(2);
	//gBkg.GetYaxis()->SetRangeUser(0, 0.0015);
	gAll.SetLineColor(kRed);
	gAll.SetLineWidth(2);

	c1.cd();
	gBkg.Draw("AC");
	gAll.Draw("Csame");

	c1.SaveAs(name.Data());
}

#endif
