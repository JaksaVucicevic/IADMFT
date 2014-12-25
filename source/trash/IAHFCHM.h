//*********************************************************//
// Imaginary Axis Half Filled Clean Hubbard Model Solver   //
//                d=inf, Bethe Lattice                     //
//                                                         //
//                       by Jaksa Vucicevic January 2010   //
//*********************************************************//

#include <cstdio>
#include <cmath> 
#include <iostream>
#include <complex>
#include "FFT.h"

class IAGRID;

using namespace std;

class IAHFCHM
{
  public:
    bool Initialized;
    
    //---options---
    bool UseBethe; 
    bool print_final;
    bool print_initial;
    int print_intermediate;   
    bool print_lambda;
    bool show_status;  

   
    //--parameters--
    int N;                    //number of points in arrays
    double Temp;              //temperature 
    double U;                 //interaction
    int metal;               //start from
    double t;                 //hopping
    double min_accr;          //accr demand

    //Green's functions and Sigmas
    complex<double>* G_f;
    complex<double>* G_0_f;
    complex<double>* G_t;
    complex<double>* Sigma_t;
    complex<double>* Sigma_f;
   
    //Fast Fourier Transform object
    FFT fft;
    IAGRID* iagrid;
    double * omega;
    double * tau;

    //--calculation routines--
    void Init_G_f();
    void get_G_0_f();
    void get_Sigma_t();
    void get_G_f();
    
    //---output---

    void ShowStatus(int msg);
    void ShowStatus(int it, double accr);

  public:
    //--user interface--
    IAHFCHM();
    ~IAHFCHM();
    void Initialize(int N, 
                    double Temp, 
                    double U, 
                    int metal,
                    double t,
                    double min_accr);
    void SetMetal(int metal);
    void SetUT(double U, double Temp);
    void SetOptions(bool print_final, 
                    int print_intermediate, 
                    bool print_initial, 
                    bool print_lambda, 
                    bool show_status);

    double Run(int ITS, bool INIT_GF);
    void GetResults(double omega_n[], complex<double> G_f[]);

    double DoubleOccupancy();
    complex<double> InternalEnergy();
    double ChargeSusceptibility();
    //complex<double> OpticalConductivity(int n);
    complex<double> OpticalConductivity(int n, double (*v)(double)=NULL);
    complex<double> Lambda(int n, double (*v)(double)=NULL); //current-current correlator
    complex<double> OpticalConductivity(double nu);
    complex<double> OpticalConductivity2(double nu);

};

