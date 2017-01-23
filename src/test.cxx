#include <iostream>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TStopwatch.h>

#include "NMRFastAna.h"

using namespace std;

int main(int argc, char* argv[])
{
	//consts
	double gainarray[3] = {1., 20., 200.};

	//Load the Q and TE data
	double FreqCenter;
	double FreqStep;
	double ScanPoints;
	double Gain;
	vector<double>* array; 

	//Load Q-curve
	TFile* qcFile = new TFile(argv[1], "READ");
	TTree* qcTree = (TTree*)qcFile->Get("NMRtree");

	qcTree->SetBranchAddress("array", &array);
	qcTree->SetBranchAddress("FreqCenter", &FreqCenter);
	qcTree->SetBranchAddress("FreqStep", &FreqStep);
	qcTree->SetBranchAddress("ScanPoints", &ScanPoints);
	qcTree->SetBranchAddress("Gain", &Gain);

	cout << "Loading Q-curve from file " << argv[1] << endl;

	qcTree->GetEntry(0);
	int nScanPoints = int(ScanPoints+0.1);
	double minFreq = FreqCenter - (ScanPoints-1.)/2.*FreqStep;
	double maxFreq = FreqCenter + (ScanPoints-1.)/2.*FreqStep;

	vector<double> qcurve_arr(nScanPoints, 0.);
	vector<double> freq_arr(nScanPoints, 0.);
	for(int i = 0; i < qcTree->GetEntries(); ++i)
	{
		qcTree->GetEntry(i);

		for(int j = 0; j < nScanPoints; ++j)
		{
			freq_arr[j] = minFreq + j*FreqStep;
			qcurve_arr[j] = qcurve_arr[j] + array->at(j)/gainarray[int(Gain+0.01)];
		}
	}

	TH1D* qcHist = new TH1D("qcHist", "qcHist", nScanPoints, minFreq, maxFreq);
	qcHist->Sumw2(); 
	for(int i = 0; i < nScanPoints; ++i) 
	{
		qcurve_arr[i] = qcurve_arr[i]/qcTree->GetEntries();
		qcHist->Fill(freq_arr[i], qcurve_arr[i]);

	}
	//for(int i = 0; i < nScanPoints; ++i) cout << i << "  " << qcurve_arr[i] << "  " << freq_arr[i] << endl;

	//Load TE
	TFile* teFile = new TFile(argv[2], "READ");
	TTree* teTree = (TTree*)teFile->Get("NMRtree");

	teTree->SetBranchAddress("array", &array);
	teTree->SetBranchAddress("FreqCenter", &FreqCenter);
	teTree->SetBranchAddress("FreqStep", &FreqStep);
	teTree->SetBranchAddress("ScanPoints", &ScanPoints);
	teTree->SetBranchAddress("Gain", &Gain);

	cout << "Loading TE data from " << argv[2] << endl;

	//Initialize NMRFastAna, the parameters are not very carefully tuned
	NMRFastAna* fastAna = new NMRFastAna();
	fastAna->setFreqRange(212.7, 213.3);
	fastAna->setExclusionWin(213.0, 0.1);
	fastAna->setSampleRate(1);
	fastAna->setXOffsetRange(-50, 50);
	fastAna->setYOffsetRange(-1., 1.);
	fastAna->init();

	//Set the QCurve to fastAna
	fastAna->setQCurve(qcHist);

	//Loop over all TE entries
	TStopwatch timer;
	timer.Start();
	for(int i = 0; i < teTree->GetEntries(); ++i)   //This is equivalent to NMRana::Loop()
	{
		teTree->GetEntry(i);

		TH1D* teHist = new TH1D("teHist", "teHist", nScanPoints, minFreq, maxFreq);
		for(int j = 0; j < nScanPoints; ++j)
		{
			teHist->Fill(freq_arr[j], array->at(j)/gainarray[int(Gain+0.01)]);
		}
		cout << "Loop " << i << ": xOffset = " << fastAna->getXOffset(teHist) << ", yOffset = " << fastAna->getYOffset() << endl;
		//fastAna->plot(Form("res_%d.pdf", i));

		delete teHist;
	}
	timer.Stop();
	timer.Print();

	delete fastAna;
	return 0;
}