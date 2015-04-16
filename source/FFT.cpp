#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <complex>
#include "IAGRID.h"
#include "FFT.h"
#include "routines.h"


using namespace std;


//------------IMPLEMENTATION------------------------------------//
//--------------------------------------------------------------//

//--------------------General stuff---------------------

double FFT::get_a_n(int n)
{
  return (n+0.5)*pi/N;
}

double FFT::get_A_n(int n)
{ 
  return sqr( sin(a_n[n])/a_n[n] );
}

double FFT::get_C(complex<double> G_f_max)
{ 
  return real(complex<double>(1.0, 0.0) / G_f_max);
}

double FFT::get_Omega()
{
  return (N+0.5)*domega;
}

//--------------Little Algebra-----------------------------------
 
complex<double> FFT::Stage_1(int n)
{ 
  return exp(complex<double> (0.0,-(n+1.0)*pi/N));
}

complex<double> FFT::Stage_3(int m)
{
  return exp(complex<double> (0.0,-pi*(m-0.5)/N)); 
}

//---------------------F to T-----------------------------------------//
//--------------------------------------------------------------------//


//---------F to T: Correction Terms--------------------------------

complex<double> FFT::CorrTerm_FtoT_1(int m)
{
  return -1.0 + ( SI((2.0*(m+1)-1.0)*pi) + SI((2.0*N-2.0*(m+1)+1.0)*pi) )/pi;
}

complex<double> FFT::CorrTerm_FtoT_2(int m)
{
  return tau_m[m]/2.0 - cos(tau_m[m]*Omega)/(Omega*pi) 
                      - tau_m[m]*SI(tau_m[m]*Omega)/pi  ;
}

complex<double> FFT::CorrTerm_FtoT_3(int m)
{
  return  -sqr(tau_m[m])/4.0 +  tau_m[m]*cos(tau_m[m]*Omega)/(Omega*2.0*pi)  
                                + sin( tau_m[m] * Omega ) / (sqr(Omega)*2.0*pi)
                                + sqr( tau_m[m] ) * SI( tau_m[m] * Omega ) / (2.0*pi)  ;
}

complex<double> FFT::Correction_FtoT(int m, complex<double> G_f_max)
{
  return             CT_FtoT_1[m]
          +     C  * CT_FtoT_2[m]
          + sqr(C) * CT_FtoT_3[m];
}

//---------------F to T: Calculate G_0 and G_beta------------------

void FFT::Prepare_G_0_and_beta(complex<double> G_t[])
{
  double sum = 0.0;

  if (not HalfFilling)
  { G_0 = CubicFrom4points(G_t, tau_m);
    G_beta = - 1.0 - G_0;
  }
  else
  { G_0 = -0.5;
    G_beta = - 0.5;
  }
  
  //ovo treba videti da li je dobro
/*  for (int n=0; n<N; n++)
  {
    sum+=real(G_f[n]);
  }
  G_0 = 2.0*Temp*sum - 0.5;
  G_beta = -2.0*Temp*sum - 0.5;
*/
}

//--------------------------------------------------------------------



//-----------------------T to F----------------------------------------//
//---------------------------------------------------------------------//


//-----------T to F: Correction Terms-------------------------

complex<double> FFT::get_dG(int n, complex<double> G_t[], int one_or_N)
{ 
  return complex<double>( dtau * A_n[n], 0.0 ) 
         * G_t[one_or_N-1]
         * exp( complex<double>( 0.0, omega_n[n]*tau_m[one_or_N-1] ) );
}

complex<double> FFT::CorrTerm_TtoF_1(int n, complex<double> G_t[], int one_or_N)
{ 
  complex<double> xi = (one_or_N==1) ? complex<double>(1.0 ,0.0) : complex<double>(-1.0 ,0.0);
  complex<double> G_0_or_beta = (one_or_N==1) ? G_0 : G_beta;
  
  return   complex<double>( 0.0, -B_1_n[n] )
              * G_t[one_or_N-1] * exp( complex<double> ( 0.0, real(xi)*a_n[n] ) )
         + complex<double>( 0.0, B_1_n[n] )
              * G_0_or_beta 
         + xi * complex<double>( B_2_n[n], 0.0 ) 
              * ( G_t[one_or_N-1] - G_0_or_beta ) 
              * ( exp( complex<double> (0.0, real(xi)*a_n[n]) ) - 1.0 );
}

complex<double> FFT::CorrTerm_TtoF_2(int n, complex<double> G_t[], int one_or_N)
{  
  complex<double> xi = (one_or_N==1) ? complex<double>(1.0 ,0.0) : complex<double>(-1.0 ,0.0);
  complex<double> G_0_or_beta = (one_or_N==1) ? G_0 : G_beta;
  
  return   get_dG(n,G_t,one_or_N)  
         + G_t[one_or_N-1] * exp( complex<double>(0.0, real(xi) * a_n[n]) ) 
                         * complex<double> ( -real(xi) * B_3_n[n], B_1_n[n] )
         + xi* G_t[one_or_N-1] * exp( complex<double>(0.0, -real(xi)*a_n[n]) )
                             * complex<double> (B_3_n[n], 0.0);
}


complex<double> FFT::Correction_TtoF(int n, complex<double> G_t[])
{
  return   CorrTerm_TtoF_1(n, G_t, 1)
         + CorrTerm_TtoF_1(n, G_t, N)
         + CorrTerm_TtoF_2(n, G_t, 1)
         + CorrTerm_TtoF_2(n, G_t, N);
}

//-------FAST FOURIER TRANSFORM ROUTINE----------------------//
//                from Numerical recipes                     //
//-----------------------------------------------------------//

//--------------F to T: Prepare input data-----------

void FFT::PrepareData_FtoT(complex<double> G_f[], double data[])
{ 
  for (int n=0; n<N; n++)
  {
    data[2*n+1] = real(G_f[n]*ST_1[n] );
    data[2*n+2] = imag(G_f[n]*ST_1[n] );
  } 
}

//--------------T to F: Prepare input data-----------

void FFT::PrepareData_TtoF(complex<double> G_t[], double data[])
{
  for (int m=0; m<N; m++) 
  {
    data[2*m+1] = real( G_t[m] * conj(ST_1[m]) );
    data[2*m+2] = imag( G_t[m] * conj(ST_1[m]) );
  } 
}

//----------Fast Fourier-----------------------------
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void FFT::FastFourier(double data[], unsigned long N, int isign)
{
	unsigned long n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta;
	double tempr,tempi;

	n = N << 1;
	j=1;
	for (i=1;i<n;i+=2) {
		if (j > i) {
			SWAP(data[j],data[i]);
			SWAP(data[j+1],data[i+1]);
		}
		m=n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax=2;
	while (n > mmax) {
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=n;i+=istep) {
				j=i+mmax;
				tempr=wr*data[j]-wi*data[j+1];
				tempi=wr*data[j+1]+wi*data[j];
				data[j]=data[i]-tempr;
				data[j+1]=data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}
#undef SWAP
//----------------------------------------------------------------.//


//---------------------------------------------------------------------//
//                  USER INTERFACE FUNCTIONS                           //
//---------------------------------------------------------------------//


//-----------Precalculate expressions------------------


void FFT::ReleaseMemory()
{
  delete [] a_n;
  delete [] A_n;
  delete [] ST_1;
  delete [] ST_3;
  delete [] CT_FtoT_1;
  delete [] CT_FtoT_2;
  delete [] CT_FtoT_3;
  delete [] B_1_n;
  delete [] B_2_n;
  delete [] B_3_n;
}

void FFT::Initialize(int N, double T, double* omega, double* tau)
{
  if (Initialized) ReleaseMemory();

  HalfFilling = false;
  this->N = N;
  Temp=T;
  omega_n=omega;
  tau_m=tau;
  domega = omega[1]-omega[0];
  dtau = tau[1]-tau[0];

  //----general stuff----
  a_n = new double [N];
  A_n = new double [N];
  
  //---little algebra------
  ST_1 = new complex<double> [N];
  ST_3 = new complex<double> [N];
  
  //------F to T---------
  Omega = get_Omega();
  CT_FtoT_1 = new complex<double> [N];
  CT_FtoT_2 = new complex<double> [N];
  CT_FtoT_3 = new complex<double> [N]; 

  //------T to F---------
  B_1_n = new double [N];
  B_2_n = new double [N];
  B_3_n = new double [N];

  //---precalculated arrays-----
  for (int i=0; i<N; i++)
  {   
    a_n[i] = get_a_n(i);
    A_n[i] = get_A_n(i);
    
    CT_FtoT_1[i]=CorrTerm_FtoT_1(i);  
    CT_FtoT_2[i]=CorrTerm_FtoT_2(i); 
    CT_FtoT_3[i]=CorrTerm_FtoT_3(i); 

    B_1_n[i] = dtau / (2*a_n[i]);
    B_2_n[i] = B_1_n[i] / a_n[i];
    B_3_n[i] = B_2_n[i] / 2.0;

    ST_1[i] = Stage_1(i);
    ST_3[i] = Stage_3(i);
  }
  Initialized = true;
}

void FFT::SetHalfFilling(bool HalfFilling)
{
  this->HalfFilling = HalfFilling;
}


//----------------Freq to Time ---------------------------



void FFT::FtoT(complex<double> G_f[], complex<double> G_t[])
{ 
  if (!Initialized)
  {
    cout << "FFT not initialized!" << endl;
    exit(1);
  }  

  C = get_C(G_f[N-1]);
  
  //------fast four---------
  double* data = new double [2*N + 1];
  PrepareData_FtoT(G_f,data);
  
  FastFourier(data, N, -1);
  
  //--extracting data--
  
  for (int m=0; m<N; m++)   
    G_t[m] = 2.0 * Temp 
                 * real( complex<double>( data[2*m+1], data[2*m+2] ) 
                         * ST_3[m] )
             + Correction_FtoT(m, G_f[N]);
      
  delete [] data;
  
}

//---------------Time to Freq ----------------------------



void FFT::TtoF(complex<double> G_t[], complex<double> G_f[])
{ 
  if (!Initialized)
  {
    cout << "FFT not initialized!" << endl;
    exit(1);
  }

  //---fast fourier------------
  double* data = new double [2*N + 1];
  PrepareData_TtoF(G_t,data);
 
  FastFourier(data, N, +1);
 
  //--extracting data--
  for (int n=0; n<N; n++)  
  {
    G_f[n] = dtau * A_n[n] * ( complex<double>( data[2*n+1], data[2*n+2] ) 
                               * conj( ST_3[n] )
                               - exp( complex<double>( 0.0, omega_n[n]*tau_m[0] ) )
                                 * G_t[0]
                               - exp( complex<double>( 0.0, omega_n[n]*tau_m[N-1] ) )
                                 * G_t[N-1] );
  }

  Prepare_G_0_and_beta(G_t);
  for (int n=0; n<N; n++)  
    G_f[n] += Correction_TtoF(n, G_t);

  delete [] data;
}

//--------default constructor-----------------------

FFT::FFT()
{
  Initialized = false;
}

//--------------------------------------------------------------//
