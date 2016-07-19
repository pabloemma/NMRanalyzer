/*
 * NMRana.cpp
 *
 *  Created on: May 10, 2016
 *      Author: klein
 *      analyzes the NMR rootfile
 */

#define NMRana_cxx
#define TEana_cxx


#include <unistd.h>  // to get gnu getoption
#include <ctime> // to manipulate time
#include <stdlib.h>
#include  <vector>
#include <string>

#include <iostream>

#include "NMRana.h"
#include "NMRana_main.h"
#include "TEana.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TObjArray.h>
#include <TApplication.h>


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



	TApplication *theApp = new TApplication("theApp",&argc,argv);  // problem with the two commandline args getting in conflict
	theApp->GetOptions(&argc,argv);



	InputRootDirectory = theApp->WorkingDirectory();


	// First we look at the filenames and put thenm into a TObjArray
	// so that we can then work thorugh them in order.
	// If the rootname starts with QCR we deal with a Qcurve
	//TER is TE measurement
	// and NMR is signal.


	for (Int_t k=0;k < theApp->Argc();k++){  // find the filenames
	TString temp = std::string(theApp->Argv(k));  //Argv with index is char*
		if(temp.Contains("-f")){
			parameter_file = std::string(theApp->Argv(k+1));
		}

		if(temp.Contains(".root")) {
			// Now check what kind of run it is
			if(temp.BeginsWith("POL")|| temp.BeginsWith("TER") || temp.BeginsWith("NMR")) {
				InputSignalFile.push_back(InputRootDirectory+temp);
			}
			else if(temp.BeginsWith("QCR")) {
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
// Close the file if only one
	}









	if(!InputSignalFile.empty() || !InputTEFile.empty() ){
			SIG.ReadParameterFile(parameter_file);
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
// Close the file if only one


	}





    theApp->Run();

    theApp->Terminate();

    // close file
     // Close the file if only one

	if(InputQcurveFile.size() == 1)QCU.CloseFile();
	if(InputSignalFile.size() == 1 || InputTEFile.size() == 1)SIG.CloseFile();

    return 0;




} //end of main loop
