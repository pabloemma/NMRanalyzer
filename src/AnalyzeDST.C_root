/*
 * AnalyzeDST.C
 *
 *  Created on: Feb 27, 2018
 *      Author: klein
 */


#include <unistd.h>
void AnalyzeDST(TString DSTfile){

// this analyzes DST files produced by NMRana
// Currently the DST file has a few hist saved at the end and one TRee
// the teee consists of the time stamp,the signal area and
//the background corrected signal

TFile *f1 = new TFile(DSTfile,"r"); // only open for read

// bocl kf definitions
TH1D * NMR_RT_Corr_Fit =0;

TH1D * NMR_RT_Corr =0;
TH1D * NMR_RT =0;
TH1D * NMR_Copy =0;
TH1D * NMR1_Qfit = 0;
TH1D *integ1= new TH1D("integ1","area distribution",200,0.0002,.001);
TH1D *integ2 = new TH1D("integ2","area distribution",200,0.0002,.001);
TH1D *integ3 = new TH1D("integ3","area distribution Qfit",200,0.0002,.001);


Double_t *SignalArea =0; // calculated signal area
Long64_t *timel =0; // time stamp


// Now get the bloody tree from the file
TTree *Dtree = (TTree*)f1->Get("Dtree");

// set branch addresses first
Dtree->SetBranchAddress("NMR_RT_Corr_Fit",&NMR_RT_Corr_Fit);
Dtree->SetBranchAddress("NMR_RT_Corr",&NMR_RT_Corr);
Dtree->SetBranchAddress("NMR_RT",&NMR_RT);
Dtree->SetBranchAddress("NMR1_Qfit",&NMR1_Qfit);
//Dtree->SetBranchAddress("timel",&timel);
//Dtree->SetBranchAddress("SignalArea",&SignalArea);

// Get number of entries
Long64_t nentries = Dtree->GetEntries();
//cout<<nentries<< "   we have so many entries \n";

TCanvas * c1 = new TCanvas("c1","Analysis of DST",500,50,800,800);
	c1->Divide(2,2);


//start Loop
	for(int i = 0; i < nentries; ++i){
	//for(int i = 0; i < 50; ++i){
	//for(int i = 0; i < 10; ++i){
	 NMR_RT_Corr_Fit->Reset();
	 NMR_RT_Corr->Reset();
	 NMR_RT->Reset();
	 NMR1_Qfit->Reset();
	 Dtree->GetEntry(i);
//	 cout<<i<<endl;
	 c1->cd(1);
	 NMR_RT_Corr_Fit->Draw("HIST P");
	 c1->cd(2);
	 NMR_RT_Corr->Draw("HIST P");


//	 NMR_RT->Draw("HIST P");
         NMR_Copy = (TH1D*)NMR_RT_Corr->Clone("NMR_Copy");
	 NMR_Copy->Add(NMR_RT_Corr_Fit,-1.);
	 NMR_Copy->Draw("HIST P SAME");
	 	 c1->cd(3);
		 NMR1_Qfit->SetMarkerStyle(24);
		 NMR1_Qfit->SetMarkerSize(.6);
	 NMR1_Qfit->Draw("HIST P ");
	 c1->cd(4);
	 NMR_RT_Corr->SetMarkerStyle(22);
	 NMR_RT_Corr->SetMarkerSize(.6);
	 NMR_RT_Corr->SetMarkerColor(2);
	 NMR_RT_Corr_Fit->SetMarkerColor(4);
	 NMR_RT_Corr_Fit->SetMarkerStyle(23);
	 NMR_RT_Corr_Fit->SetMarkerSize(.6);

	 NMR_RT_Corr_Fit->Draw("HIST P");
	 NMR_RT_Corr->Draw("HIST P SAME");
	 //NMR_RT->Draw("HIST P SAME");
	 //
	 //now intergrate the area and multiply with channel width
         integ1->Fill((NMR_RT_Corr_Fit->Integral(119,301))*.002);
         integ2->Fill((NMR_RT_Corr->Integral(119,301))*.002);
         integ3->Fill((NMR1_Qfit->Integral(119,301))*.002);

//	 cout<<"integral"  <<temp<<endl;

	  c1->Modified();
	  c1->Update();
	  //sleep(1);

	 }
TCanvas * c2 = new TCanvas("c2","Integral of DST",10,50,800,800);
	c2->Divide(1,2);
	c2->cd(1);
		 integ1->SetMarkerStyle(22);
	 integ1->SetMarkerSize(.6);
	 integ2->SetMarkerSize(.6);
	 integ2->SetLineColor(3);

	integ1->Draw("HIST ");
	c2->cd(2);
	integ2->Draw("HIST ");
	integ3->SetLineColor(1);
	integ3->Draw("HIST SAME");
	  c2->Modified();
	  c2->Update();


} // end of analyzer



