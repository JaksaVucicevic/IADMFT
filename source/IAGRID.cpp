#include "IAGRID.h"
#include "routines.h"

IAGRID::IAGRID(int N, int M, double T)
{
  this->N = N;
  this->M = M;
  this->T = T;
  
  beta = 1.0/T;
}

IAGRID::~IAGRID()
{
  //add something maybe?
}

double IAGRID::get_tau(int m)
{ 
  return ( (double) m + 0.5 )
         /
         (  ((double) N) * T  );
}

double IAGRID::get_omega(int n)
{
  return 2.0*(n+0.5)*pi*T;
}

void IAGRID::assign_omega(double* omega)
{
  for(int n=0; n<N; n++)
    omega[n] = get_omega(n);
} 

void IAGRID::assign_tau(double* tau)
{
  for(int m=0; m<M; m++)
    tau[m] = get_tau(m);
} 


