/*
 * NMRana.cpp
 *
 *  Created on: May 10, 2016
 *      Author: klein
 *      analyzes the NMR rootfile
 *      modified to work on linux
 */

#define NMRana_cxx
#define TEana_cxx


#include <unistd.h>  // to get gnu getoption
#include <ctime> // to manipulate time
#include <stdlib.h>
#include  <vector>
#include <string>
#include <cstdlib> // contaings getenv

#include <iostream>
#include "NMRana.h"
#include "NMRana_main.h"

#include "TEana.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TObjArray.h>
#include <TApplication.h>
#include <TStopwatch.h>

//#include "gperftools/profiler.h"


#ifdef __APPLE__
 std::string OS = "OSX";
#elif __linux__
 std::string OS = "OSX";
 #elif _WIN64
 std::string OS = "WIN";
 // add windows includes
#elif _WIN32
 std::string OS = "WIN";
 // add windows includes
#endif

 NMRana SIG;			// create object for signal
// TEana TE;			// create object for TE
 NMRana QCU;			// create object QCU



using namespace std;

#ifdef __APPLE__
int main(Int_t argc,char **argv) {
#endif

#ifdef __linux__
int main(Int_t argc,char **argv) {
#endif

#ifdef __WIN64
int main(Int_t argc,char *argv[],char *envp[] ) {
#endif

#ifdef __WIN32
int main(Int_t argc,char *argv[],char *envp[] ) {
#endif




/*     Get and handle command line options. Arguments handled are removed
 from the argument array. The following arguments are handled:
    -b : run in batch mode without graphics
    -x : exit on exception
    -n : do not execute logon and logoff macros as specified in .rootrc
    -q : exit after processing command line macro files
    -l : do not show splash screen
 The last three options are only relevant in conjunction with TRint.
 The following help and info arguments are supported:
    -?       : print usage
    -h       : print usage
    --help   : print usage
    -config  : print ./configure options
    -memstat : run with memory usage monitoring
 In addition to the above options the arguments that are not options,
 i.e. they don't start with - or + are treated as follows (and also removed
 from the argument array):
   <dir>       is considered the desired working directory and available
               via WorkingDirectory(), if more than one dir is specified the
               first one will prevail
   <file>      if the file exists its added to the InputFiles() list
   <file>.root are considered ROOT files and added to the InputFiles() list,  **this only works in TRint
               the file may be a remote file url
   <macro>.C   are considered ROOT macros and also added to the InputFiles() list
 In TRint we set the working directory to the <dir>, the ROOT files are
 connected, and the macros are executed. If your main TApplication is not
 TRint you have to decide yourself what to do whith these options.
 All specified arguments (also the ones removed) can always be retrieved
 via the TApplication::Argv() method.
*/


	// Get the environment, make sure it is defined. If not exit
			if(const char* env_p = std::getenv("NMR_ROOT")){
				NMR_ROOT = env_p;
				cout<< "NMR_main" <<" your environment is "<<NMR_ROOT<<endl;
				SIG.SetEnvironment(NMR_ROOT);
			}
			else{
				cout<< "NMR_main"<<" your environment is not defined; most likely you need to set NMR_ROOT environment variable \n ";
				cout<<"NMR_main"<<" This means that the program does not find the directory where the par files and qc curve etc are located \n";
				cout<<"NMR_main"<<" You should have a directory, which has the following subdirectories: NMR_par,QC_files, TE_files \n";
				cout<<"NMR_main"<<" I define one for you, but most likely this will bomb and you have to go into the main code and change the line around line 118 \n";
				//std::string env_p1="/Volumes/macsmall/nmrwork/";
				//NMR_ROOT= env_p1;
				cout<<"NMR_main"<<NMR_ROOT<<"   this is what I define for you \n";
						//cout<< " You must define it now; Good bye  \n" <<endl;
						//exit(EXIT_FAILURE);
				SIG.SetEnvironment(NMR_ROOT);


			}




	TApplication *theApp = new TApplication("theApp",&argc,argv);  // problem with the two commandline args getting in conflict
	theApp->GetOptions(&argc,argv);


	extern Bool_t gPrintViaErrorHandler;

	gPrintViaErrorHandler = true;


	//	extern int gErrorIgnoreLevel;
//	gErrorIgnoreLevel = 9999;


	InputRootDirectory = theApp->WorkingDirectory();

	if(InputRootDirectory =="") InputRootDirectory = "/home/klein/NMRanalysis/LanlData/root/";
	cout<<"root data dir"<<InputRootDirectory<<"\n";

    Timer = new TStopwatch; // create a stopwatch
    Timer->Start();
	// First we look at the filenames and put thenm into a TObjArray
	// so that we can then work thorugh them in order.
	// If the rootname starts with QCR we deal with a Qcurve
	//TER is TE measurement
	// and NMR is signal.


	for (Int_t k=0;k < theApp->Argc();k++){  // find the filenames
	TString temp = std::string(theApp->Argv(k));  //Argv with index is char*
		if(temp.Contains("-f")){
			parameter_file = std::string(theApp->Argv(k+1));
			cout<<"parameter file"<<parameter_file<<"\n";

		}

		if(temp.Contains(".root")) {
			// Now check what kind of run it is
			if(temp.BeginsWith("POL")|| temp.BeginsWith("TER") || temp.BeginsWith("NMR")|| temp.BeginsWith("TEQ")|| temp.BeginsWith("QCV")) {
				InputSignalFile.push_back(InputRootDirectory+temp);
				// create datestring
				datestring =temp;
				datestring = datestring.Remove(13,5);
				datestring = datestring.Remove(0,3);
			}
			else if(temp.BeginsWith("QCR") ) {
				InputQcurveFile.push_back(InputRootDirectory+temp);
			}
			else {
				cout<<" the file identifier in "<<temp<< " is not recognized , fatal error , will exit \n";
				return 1;
				}
		}
	}

// Here starts the main activity
// check if the particular type of run is in the file list
// First any possible QCurve run

	if(!InputQcurveFile.empty() ){

// Open ROOT file
			if(InputQcurveFile.size()>1){
				QCU.OpenChain(InputQcurveFile);  // we will create a TChain
			}
			else {
				QCU.OpenFile(InputQcurveFile[0]); // just one spectrum
			}
// Do stuff with it
			QCU.SetupHistos();
			QCU.Loop();
			QCU.DrawHistos();
			QCU.Finish();
// Close the file if only one
	}








	// start profiling with google profuiler
//	ProfilerStart("/home/plm/profile.dat");
	if(!InputSignalFile.empty() || !InputTEFile.empty() ){
			SIG.ReadParameterFile(parameter_file);
			SIG.GetDate(datestring); // this will set time control;
// Open ROOT file
			if(InputSignalFile.size()>1){
				SIG.OpenChain(InputSignalFile);  // we will create a TChain
			}
			else {
				SIG.OpenFile(InputSignalFile[0]); // just one spectrum
			}
// Do stuff with it
			//!!!!!!!!!!!warning, has to be set differently, right now for testing hard coded limits
//			SIG.AreaSetLimits(212.83,213.21);
			//!!!!!!!!!!!!!!!!!!


			SIG.SetupCanvas();

			SIG.SetupHistos();
			SIG.Loop();
			SIG.DrawHistos();
			SIG.Finish();
// Close the file if only one


	}
//	ProfilerStop();


	cout<< "Finished working \n";
	Timer->Print();


    theApp->Run();

    theApp->Terminate();

    // close file
     // Close the file if only one

	if(InputQcurveFile.size() == 1)QCU.CloseFile();
	if(InputSignalFile.size() == 1 || InputTEFile.size() == 1)SIG.CloseFile();

    return 0;




} //end of main loop
