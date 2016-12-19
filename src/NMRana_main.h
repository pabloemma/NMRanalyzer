/*
 * NMRana_main.h
 *
 *  Created on: May 10, 2016
 *      Author: klein
 */

#ifndef NMRANA_MAIN_H_
#define NMRANA_MAIN_H_

TString InputRootDirectory ; // the directory where the root input file is located
TString root_infile; // the current root file


std::vector<TString> InputSignalFile ; // if there is a list of input files it will put them into vector
std::vector<TString> InputQcurveFile ; // if there is a list of input files it will put them into anarry
std::vector<TString> InputTEFile ; // if there is a list of input files it will put them into anarry

TString parameter_file; //
TStopwatch *Timer;
TString datestring;

std::string  NMR_ROOT ; // top directory of NMR system, needs to be defined thorugh enviro variable $NMR_ROOT



#endif /* NMRANA_MAIN_H_ */
