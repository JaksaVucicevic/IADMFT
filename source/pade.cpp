#include <cstdio>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <gmpxx.h>
#include "pade.h"

using namespace std;

const double pi = 3.14159265359;
const complex<double> ii = complex<double>(0.0, 1.0);
typedef std::complex<double> dcomplex;


//--------------- GMP complex data type, taken from TRIS ---------------//

struct gmp_complex {
 mpf_class re, im;
 gmp_complex operator* (const gmp_complex &rhs){ return { rhs.re*re-rhs.im*im, rhs.re*im+rhs.im*re }; }
 friend gmp_complex inverse (const gmp_complex &rhs){ mpf_class d=rhs.re*rhs.re + rhs.im*rhs.im; return { rhs.re/d, -rhs.im/d}; }
 gmp_complex operator/ (const gmp_complex &rhs){ return (*this)*inverse(rhs); }
 gmp_complex operator+ (const gmp_complex &rhs){ return { rhs.re + re, rhs.im + im }; }
 gmp_complex operator- (const gmp_complex &rhs){ return { re - rhs.re, im - rhs.im }; }
 friend mpf_class real(const gmp_complex &rhs) { return rhs.re; }
 friend mpf_class imag(const gmp_complex &rhs) { return rhs.im; }
 gmp_complex& operator= (const std::complex<double> &rhs) { re = real(rhs); im = imag(rhs); return *this; }
 friend std::ostream & operator << (std::ostream & out,gmp_complex const & r) { return out << " gmp_complex("<<r.re<<","<<r.im<<")"<<std::endl ;}
};

mpf_class abs(mpf_class a)
{
  if ( a > 0.00 )
    return  a;
  else return -a;

}

namespace pade_local
{

void PrintFunc(const char* FileName, int N, complex<double>* Y, double* X)
{
  FILE *f;
  f = fopen(FileName, "w");
  for (int i=0; i<N; i++)
    fprintf(f,"%.15le %.15le %.15le\n", X[i], real(Y[i]), imag(Y[i]));
  fclose(f);
}

void ReadFunc(const char* FileName, int &N, complex<double>* &Y, double* &X)
{ // reads file formatted like 
  //  X  ReY(X) ImY(X)
  //----open file---//
  FILE *f;
  f = fopen(FileName, "r");
  if (f==NULL) perror ("Error opening file");
  
  //----count rows----//
  int i=0;
  char str[1000];  
  while (!feof(f))
  {  
    fgets ( str, 1000, f ); 
    if (str[0]=='\n') continue;
    i++;
  }
  N=i-1;

  printf("N: %d \n", N);
  fclose(f);
 
  X = new double[N];
  Y = new complex<double>[N];

  f = fopen(FileName, "r");
   
  for (int i=0; i<N; i++)
  { double Dummy1,Dummy2,Dummy3;
    fscanf(f, "%le", &Dummy1);
    fscanf(f, "%le", &Dummy2);
    fscanf(f, "%le", &Dummy3);
    X[i]=Dummy1;
    Y[i]=complex<double>(Dummy2,Dummy3);
  }
  fclose(f);
}

}
//------------------------ PADE -------------------------------------//

void pade( int M, double* iw, complex<double>* Giw, 
           int N, double* w, complex<double>* Gw )
{
  double eta = 0.0; 	//broadening (G gets calculated in (w+i*eta) points)
  double accr = 5e-7;	//minimum accuracy required
  int Mmin = 500;	//minimal number of matsubara points to use

  int GMP_default_prec = 256; 
  unsigned long old_prec = mpf_get_default_prec();
  mpf_set_default_prec(GMP_default_prec);

  int P = M;
  gmp_complex** g = new gmp_complex*[P];
  for(int p=0; p<P; p++)
    g[p] = new gmp_complex[M];   
    
  for(int n=0; n<M; n++)
    g[0][n] = Giw[n];

  for(int p=1; p<P; p++)   
    for(int n=p; n<M; n++)   
    {  gmp_complex gmp_one = {1.0, 0.0}; 
       gmp_complex y = {0.0, iw[n] - iw[p-1]};  
       gmp_complex x = g[p-1][p-1]/g[p-1][n] - gmp_one ;
       g[p][n] =  x / y ;
    }

  gmp_complex* A = new gmp_complex[M+1];
  gmp_complex* B = new gmp_complex[M+1];
    
  A[0].re = 0.0;
  A[0].im = 0.0;
  A[1] = g[0][0];
  B[0].re = 1.0;
  B[0].im = 0.0;
  B[1].re = 1.0;
  B[1].im = 0.0;
  
  for(int i=0; i<N; i++)
  {  
     gmp_complex G_old = {0.0,0.0}; 
     for(int n=2; n<=M; n++)
     {  gmp_complex y = {w[i], -iw[n-2]+eta}; 
        A[n] = A[n-1] + y * g[n-1][n-1] * A[n-2];
        B[n] = B[n-1] + y * g[n-1][n-1] * B[n-2];

        gmp_complex G = A[n]/B[n];

        mpf_class diff = abs( G_old.im - G.im );
        if ( (n==M)
             or
             ( ( diff < accr ) and ( n > Mmin ) ) 
           ) 
        { //printf("--- Pade done!!! i: %d w[i]: %.3f ---- n used: %d DIFF:%.2le\n", i, w[i], n, diff.get_d() );
          Gw[i] = complex<double>( G.re.get_d(), G.im.get_d());
          break;
        }
        else G_old = G;
     }
  }

  for(int p=0; p<P; p++)
    delete [] g[p];
  delete [] g;
  delete [] A;
  delete [] B;

  mpf_set_default_prec(old_prec);
}

void PadeToFile( int Mmax, complex<double>* Giw,  double* iw, const char* outputFN, int N, double wmax )
{ 
  complex<double>* Gw = new complex<double>[N];
  double* w = new double[N];
  for(int i=0; i<N; i++)
    w[i] = - wmax + i * 2.0*wmax/(N-1.0); 

  pade( Mmax, iw, Giw, 
        N,     w, Gw );

  pade_local::PrintFunc(outputFN, N, Gw, w);

  delete [] Gw;
  delete [] w;
}

int PadeFromFile(const char* FN, int Mmax, double wmax, int N)
{ 
  int M;  
  complex<double>* Giw;
  double* iw;
  pade_local::ReadFunc(FN, M, Giw, iw);

  complex<double>* Gw = new complex<double>[N];
  double* w = new double[N];
  for(int i=0; i<N; i++)
    w[i] = - wmax + i * 2.0*wmax/(N-1.0); 

  pade( Mmax, iw, Giw, 
        N,     w, Gw );

  char outputFN[200];
  sprintf(outputFN,"PadeOut.%s",FN);
  pade_local::PrintFunc(outputFN, N, Gw, w);

  delete [] Giw;
  delete [] iw;
  delete [] Gw;
  delete [] w;
}

complex<double> pade( int M, double* iw, complex<double>* Giw, double w )
{
  double* ws = new double[1];
  ws[0] = w;
  complex<double>* Gws = new complex<double>[1];

  pade( M, iw, Giw, 
        1, ws, Gws );

  complex<double> res = Gws[0];

  delete [] ws;
  delete [] Gws;
 
  return res;
}
