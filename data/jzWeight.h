// this info is for 

#ifndef JZWEIGHT_H
#define JZWEIGHT_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <TROOT.h>


// units of nb
const double JZ0_cs_PythDijet = 6.7890E+07;
const double JZ0_fltEff_PythDijet = 9.9713E-01 ;

// Note: filt ef is avg between 2.8306 and 2.8748
const double JZ1_cs_PythDijet     = 6.7890E+07;
const double JZ1_fltEff_PythDijet = 2.8527E-03 ;

// Note: filt ef is avg between 4.2774 and 4.2952
const double JZ2_cs_PythDijet     = 6.3997E+05;
const double JZ2_fltEff_PythDijet = 4.2863E-03 ;

// Note: filt ef is avg between 5.2851 and 5.2994
const double JZ3_cs_PythDijet     = 4.7194E+03;
const double JZ3_fltEff_PythDijet = 5.29225E-03;

// Note: filt ef is avg between 4.5834 and 4.5901
const double JZ4_cs_PythDijet     = 2.6602E+01;
const double JZ4_fltEff_PythDijet = 4.58675E-03;

// Note: filt ef is avg between 2.1826 and 2.1846
const double JZ5_cs_PythDijet     = 2.2476E-01;
const double JZ5_fltEff_PythDijet = 2.1836E-03 ;

const double nEvents_jz0_PythDijet  = 1000000;
const double nEvents_jz1_PythDijet  = 7999000;
const double nEvents_jz2_PythDijet  = 7998000;
const double nEvents_jz3_PythDijet  = 7999000;
const double nEvents_jz4_PythDijet  = 8000000;
const double nEvents_jz5_PythDijet  = 7999000;


// Powheg event weight defined as: event weight = (JZ_xsection * JZ_filterEff * Powheg_jetWeight) / (num_events*sum_jet_weights)
// note: jet weight is applied per jet. Corrections are stored in the ntuple and included in the reweighting section of the analysis scripts
// Pythia/herwig event weight defined as: event weight = (JZ_xsection * JZ_filterEff) / (num_events)

double perEvtWgtJZ0_PythDijet = (JZ0_cs_PythDijet * JZ0_fltEff_PythDijet ) / nEvents_jz0_PythDijet ; 
double perEvtWgtJZ1_PythDijet = (JZ1_cs_PythDijet * JZ1_fltEff_PythDijet ) / nEvents_jz1_PythDijet ;
double perEvtWgtJZ2_PythDijet = (JZ2_cs_PythDijet * JZ2_fltEff_PythDijet ) / nEvents_jz2_PythDijet ; 
double perEvtWgtJZ3_PythDijet = (JZ3_cs_PythDijet * JZ3_fltEff_PythDijet ) / nEvents_jz3_PythDijet ; 
double perEvtWgtJZ4_PythDijet = (JZ4_cs_PythDijet * JZ4_fltEff_PythDijet ) / nEvents_jz4_PythDijet ;
double perEvtWgtJZ5_PythDijet = (JZ5_cs_PythDijet * JZ5_fltEff_PythDijet ) / nEvents_jz5_PythDijet ; 

#endif 
