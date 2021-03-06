#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "IAResult.h"
#include "IAGRID.h"
#include "routines.h"
#include "IBZ.h"
#include <limits>
#include <cmath>
#include <iostream>
using namespace std;

IAResult::IAResult()
{
}

IAResult::IAResult(IAGRID* g)
{
  Initialize(g);
}

IAResult::IAResult(const IAResult &r)
{
  Initialize(r.iagrid);
  CopyFrom(r);
}

IAResult::~IAResult()
{
  ReleaseMemory();
}

void IAResult::Reset()
{
  ReleaseMemory();
  Initialize(iagrid);
}

void IAResult::Reset(IAGRID* g)
{
  ReleaseMemory();
  Initialize(g);
}

void IAResult::Initialize(IAGRID* g)
{
  this->iagrid = g;
  
  int N = iagrid->get_N();

  tau = new double[N];
  iagrid->assign_tau(tau);

  G0_tau = new complex<double>[N];
  SOCSigma_tau = new complex<double>[N];

  omega = new double[N];
  iagrid->assign_omega(omega);

  Delta = new complex<double>[N];
  G0 = new complex<double>[N];
  SOCSigma = new complex<double>[N];
  Sigma = new complex<double>[N];
  G = new complex<double>[N];

  for(int i=0; i<N; i++)
  {    
    G0_tau[i] = 0.0; 				
    SOCSigma_tau[i] = 0.0;
    Delta[i] = 0.0;			
    G0[i] = 0.0; 				
    SOCSigma[i] = 0.0; 	
    Sigma[i] = 0.0;			
    G[i] = 0.0;	
  }
  n=0.0;
  n0=0.0;
  mu=0.0;
  mu0=0.0;
}

void IAResult::ReleaseMemory()
{
  delete [] tau;
  
  delete [] G0_tau;
  delete [] SOCSigma_tau;

  delete [] omega;

  delete [] Delta;
  delete [] G0;  
  delete [] SOCSigma;
  delete [] Sigma;
  delete [] G;
}

void IAResult::PrintResult(const char* ResultFN)
{ 
  FILE *f;
  f = fopen(ResultFN, "w");
  
  fprintf(f,"# n = %le n0 = %le mu = %le mu0=%le\n",n,n0,mu,mu0);   

  int N = iagrid->get_N();
  int i;
  for (i=0; i<N; i++)
  { 
     // loop through and store the numbers into the file
    fprintf(f, "%.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le\n", 
                   tau[i],						//1
                   real(G0_tau[i]), imag(G0_tau[i]),			//2 3
                   real(SOCSigma_tau[i]), imag(SOCSigma_tau[i]), 	//4 5
                   omega[i],						//6 
                   real(Delta[i]), imag(Delta[i]),			//7 8
                   real(G0[i]), imag(G0[i]), 				//9 10
                   real(SOCSigma[i]), imag(SOCSigma[i]), 		//11 12
                   real(Sigma[i]), imag(Sigma[i]),			//13 14
                   real(G[i]), imag(G[i]));				//15 16                                
  }
  fclose(f);
}

void IAResult::PrintShort(const char* ResultFN)
{ 
  FILE *f;
  f = fopen(ResultFN, "w");
  
  fprintf(f,"# n = %le n0 = %le mu = %le mu0=%le\n",n,n0,mu,mu0);   

  int N = iagrid->get_N();
  int i;
  for (i=0; i<N; i++)
  { 
     // loop through and store the numbers into the file
    fprintf(f, "%.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le\n", 
                   tau[i],						//1
                   real(G0_tau[i]),					//2
                   real(SOCSigma_tau[i]),			 	//3
                   omega[i],						//4 
                   real(Sigma[i]), imag(Sigma[i]),			//5 6
                   real(G[i]), imag(G[i]));				//7 8                                
  }
  fclose(f);
}

void IAResult::PrintMinimal(const char* ResultFN, double cutoff)
{ 
  FILE *f;
  f = fopen(ResultFN, "w");
  
  fprintf(f,"# n = %le n0 = %le mu = %le mu0=%le\n",n,n0,mu,mu0);   

  int N = iagrid->get_N();
  int i;
  for (i=0; i<N; i++)
  { if (omega[i]>cutoff) break;
     // loop through and store the numbers into the file
    fprintf(f, "%.15le %.15le %.15le %.15le %.15le %.15le %.15le\n", 
                   omega[i],						//1
                   real(Delta[i]), imag(Delta[i]),			//2 3 
                   real(Sigma[i]), imag(Sigma[i]),			//4 5
                   real(G[i]), imag(G[i]));				//6 7                                
  }
  fclose(f);
}

void IAResult::PrintUberMinimal(const char* ResultFN, double cutoff)
{ 
  FILE *f;
  f = fopen(ResultFN, "w");
  
  fprintf(f,"# n = %le n0 = %le mu = %le mu0=%le\n",n,n0,mu,mu0);   

  int N = iagrid->get_N();
  int i;
  for (i=0; i<N; i++)
  { if (omega[i]>cutoff) break;
     // loop through and store the numbers into the file
    fprintf(f, "%.15le %.15le %.15le %.15le %.15le\n", 
                   omega[i],						//1
                   real(Sigma[i]), imag(Sigma[i]),			//2 3
                   real(G[i]), imag(G[i]));				//4 5                                
  }
  fclose(f);
}


bool IAResult::ReadFromFile(const char* ResultFN)
{ 
  FILE *f;
  f = fopen(ResultFN, "r");
  if (f==NULL) { return false; }

  char rstLine[1000];
  fgets ( rstLine, 1000, f );

  char * pch;
  printf ("rstline: %s\n",rstLine);
  pch = strtok (rstLine,"=");
  int counter=1;
  while (pch != NULL)
  { 
    switch (counter)
    {
      case 2: n=atof(pch); printf ("%d: %s, atof: %le \n",counter,pch,n); break;
      case 3: n0=atof(pch); printf ("%d: %s, atof: %le \n",counter,pch,n0); break;
      case 4: mu=atof(pch); printf ("%d: %s, atof: %le \n",counter,pch,mu); break;
      case 5: mu0=atof(pch); printf ("%d: %s, atof: %le \n",counter,pch,mu0); break;
    }
          
    pch = strtok (NULL, "=");
    counter++; 
  }

  int N = iagrid->get_N();
  int i;
  for (i=0; i<N; i++)
  { double t,rG0t, iG0t, rSSt, iSSt, o, rD, iD, rG0, iG0, rSS, iSS, rS, iS, rG, iG; 
     // loop through and store the numbers into the file
    fscanf(f, "%le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le\n", 
                   &t,				//1
                   &rG0t, &iG0t,		//2 3
                   &rSSt, &iSSt,		//4 5
                   &o,				//6 
                   &rD, &iD,			//7 8
                   &rG0, &iG0, 			//9 10
                   &rSS, &iSS, 			//11 12
                   &rS, &iS,			//13 14
                   &rG, &iG);			//15 16    

    tau[i] = t;
    G0_tau[i] = complex<double>(rG0t,iG0t); 				
    SOCSigma_tau[i] = complex<double>(rSSt,iSSt);
    omega[i] = o;
    Delta[i] = complex<double>(rD,iD);			
    G0[i] = complex<double>(rG0,iG0); 				
    SOCSigma[i] = complex<double>(rSS,iSS); 	
    Sigma[i] = complex<double>(rS,iS);			
    G[i] = complex<double>(rG,iG);			
  }
  fclose(f);
  
  return true;
}

bool IAResult::ReadFromMinimal(const char* ResultFN)
{ 
  FILE *f;
  f = fopen(ResultFN, "r");
  if (f==NULL) { return false; }

  char rstLine[1000];
  fgets ( rstLine, 1000, f );

  char * pch;
  printf ("rstline: %s\n",rstLine);
  pch = strtok (rstLine,"=");
  int counter=1;
  while (pch != NULL)
  { 
    switch (counter)
    {
      case 2: n=atof(pch); printf ("%d: %s, atof: %le \n",counter,pch,n); break;
      case 3: n0=atof(pch); printf ("%d: %s, atof: %le \n",counter,pch,n0); break;
      case 4: mu=atof(pch); printf ("%d: %s, atof: %le \n",counter,pch,mu); break;
      case 5: mu0=atof(pch); printf ("%d: %s, atof: %le \n",counter,pch,mu0); break;
    }
          
    pch = strtok (NULL, "=");
    counter++; 
  }

  int N = iagrid->get_N();
  int i;
  for (i=0; i<N; i++)
  { double o, rS, iS, rG, iG, rD, iD; 
     // loop through and store the numbers into the file
    if (feof(f)) 
      o = rS = iS = rG = iG = rD = iD = std::numeric_limits<double>::quiet_NaN();
    else       
      fscanf(f, "%le %le %le %le %le %le %le\n", 
                 &o,				//1
                 &rD, &iD,			//2 3
                 &rS, &iS,			//4 5
                 &rG, &iG);			//6 7        

    tau[i] =  std::numeric_limits<double>::quiet_NaN();
    G0_tau[i] =  std::numeric_limits<double>::quiet_NaN(); 				
    SOCSigma_tau[i] = std::numeric_limits<double>::quiet_NaN();
    omega[i] = o;
    Delta[i] = complex<double>(rD,iD);			
    G0[i] =  std::numeric_limits<double>::quiet_NaN(); 				
    SOCSigma[i] =  std::numeric_limits<double>::quiet_NaN(); 	
    Sigma[i] = complex<double>(rS,iS);			
    G[i] = complex<double>(rG,iG);			
  }
  fclose(f);
  
  return true;
}

bool IAResult::ReadFromUberMinimal(const char* ResultFN)
{ 
  FILE *f;
  f = fopen(ResultFN, "r");
  if (f==NULL) { return false; }

  char rstLine[1000];
  fgets ( rstLine, 1000, f );

  char * pch;
  printf ("rstline: %s\n",rstLine);
  pch = strtok (rstLine,"=");
  int counter=1;
  while (pch != NULL)
  { 
    switch (counter)
    {
      case 2: n=atof(pch); printf ("%d: %s, atof: %le \n",counter,pch,n); break;
      case 3: n0=atof(pch); printf ("%d: %s, atof: %le \n",counter,pch,n0); break;
      case 4: mu=atof(pch); printf ("%d: %s, atof: %le \n",counter,pch,mu); break;
      case 5: mu0=atof(pch); printf ("%d: %s, atof: %le \n",counter,pch,mu0); break;
    }
          
    pch = strtok (NULL, "=");
    counter++; 
  }

  int N = iagrid->get_N();
  int i;
  for (i=0; i<N; i++)
  { double o, rS, iS, rG, iG; 
     // loop through and store the numbers into the file
    if (feof(f)) 
      o = rS = iS = rG = iG = std::numeric_limits<double>::quiet_NaN();
    else       
      fscanf(f, "%le %le %le %le %le\n", 
                 &o,				//1
                 &rS, &iS,			//2 3
                 &rG, &iG);			//4 5        

    tau[i] =  std::numeric_limits<double>::quiet_NaN();
    G0_tau[i] =  std::numeric_limits<double>::quiet_NaN(); 				
    SOCSigma_tau[i] = std::numeric_limits<double>::quiet_NaN();
    omega[i] = o;
    Delta[i] =  std::numeric_limits<double>::quiet_NaN();			
    G0[i] =  std::numeric_limits<double>::quiet_NaN(); 				
    SOCSigma[i] =  std::numeric_limits<double>::quiet_NaN(); 	
    Sigma[i] = complex<double>(rS,iS);			
    G[i] = complex<double>(rG,iG);			
  }
  fclose(f);
  
  return true;
}

/*
bool IAResult::ReadFromFile(const char* ResultFN, int Ncolumns, 
                            int Nrqs, const double* rqs[], const int ris[], 
                            int Ncqs, const complex<double>* cqs [], const int cis[])
{ 
  FILE *f;
  f = fopen(ResultFN, "r");
  if (f==NULL) { return false; }

  char rstLine[1000];
  fgets ( rstLine, 1000, f );

  char * pch;
  printf ("rstline: %s\n",rstLine);
  pch = strtok (rstLine,"=");
  int counter=1;
  while (pch != NULL)
  { 
    switch (counter)
    {
      case 2: n=atof(pch); printf ("%d: %s, atof: %le \n",counter,pch,n); break;
      case 3: n0=atof(pch); printf ("%d: %s, atof: %le \n",counter,pch,n0); break;
      case 4: mu=atof(pch); printf ("%d: %s, atof: %le \n",counter,pch,mu); break;
      case 5: mu0=atof(pch); printf ("%d: %s, atof: %le \n",counter,pch,mu0); break;
    }
          
    pch = strtok (NULL, "=");
    counter++; 
  }

  int N = iagrid->get_N();
  int i;
  for (i=0; i<N; i++)
  { 
    for(int j=0; j>Ncolumns; j++)
    { double dum; 
      fscanf(f, "%le", &dum);
      for(int k=0; k<Nrqs; k++)
        if (ris[k]==j) rqs[k][i]=dum;
      for(int k=0; k<Ncqs; k++)
        if (cis[k]==j) cqs[k][i]=dum;       
      for(int k=0; k<Ncqs; k++)
        if (cis[k]==j-1) cqs[k][i]+=ii*dum;
    }
  }
  fclose(f);
  
  return true;
}
*/

void IAResult::CopyFrom(const IAResult &r)
{
  Reset(r.iagrid);

  int N = r.iagrid->get_N();

  for (int i=0; i<N; i++)
  {
    //printf("m"); 
    tau[i] = r.tau[i];
    G0_tau[i] = r.G0_tau[i];
    SOCSigma_tau[i] = r.SOCSigma_tau[i];
    omega[i] = r.omega[i];
    Delta[i] = r.Delta[i];
    G0[i] = r.G0[i];
    SOCSigma[i] = r.SOCSigma[i];
    Sigma[i] = r.Sigma[i];
    G[i] = r.G[i];

  }
  n = r.n;
  n0 = r.n0;
  mu = r.mu;
  mu0 = r.mu0;
}


complex<double> IAResult::OpticalConductivity(int n, IBZ* ibz)
{
  return (1.0/(2.0*pi*n))* (Lambda(n, ibz) - Lambda(0, ibz));
}

complex<double> IAResult::Lambda(int n, IBZ* ibz)
{
  //-- bosonic frequency
 
  int N = iagrid->get_N();

  double* iw_large = new double[2*N];
  complex<double>* Sig_large = new complex<double>[2*N];
  for(int i = 0; i < N; i++)
  {
    iw_large[N+i]    =  omega[i];
    iw_large[N-1-i]  = -omega[i];        
    Sig_large[N+i]   =  Sigma[i];
    Sig_large[N-1-i] = conj(Sigma[i]); // CONJ: G(iw)=G*(-iw) 
  }

  complex<double> sum = 0.0;
  int Nsum = 200;
  //for (int m=0; m<2*N-n; m++)
  for (int m=N-1-Nsum; m<N+Nsum; m++)
  { 
    for(int i=0; i<ibz->get_Nx(); i++)
    for(int j=0; j<ibz->get_Ny(); j++)
    { 
      double v = ibz->velocity[i][j];
      double eps = ibz->epsilon[i][j];
      ibz->summand[i][j] =   
                sqr(v) 
              * 2.0 
              * ( 1.0
                  / ( ( ii*iw_large[m] + mu - eps - Sig_large[m] ) 
                       * 
                      ( ii*iw_large[m+n] + mu - eps - Sig_large[m+n] ) 
                    )
                );
    }
    sum += ibz->sum();
    /*  if ((m==N)and(n==1)) 
      { 
        char FN[300];
        sprintf(FN,"summand.n%d",n);
        ibz->PrintToFile(FN);
      }
    */
  }

  delete [] iw_large;
  delete [] Sig_large;

  return sum;
}


complex<double> IAResult::OpticalConductivity(int n, double (*v)(double))
{
  return (1.0/(2.0*pi*n))* (Lambda(n, v) - Lambda(0, v));
}

complex<double> IAResult::Lambda(int n, double (*v)(double))
{
  double t = 0.5; //!!!!!!!!!!!!!!!
  int N = iagrid->get_N();

  double* iw_large = new double[2*N];
  complex<double>* Sig_large = new complex<double>[2*N];
  for(int i = 0; i < N; i++)
  {
    iw_large[N+i]    =  omega[i];
    iw_large[N-1-i]  = -omega[i];        
    Sig_large[N+i]   =  Sigma[i];
    Sig_large[N-1-i] = conj(Sigma[i]); //!!!!!!!!!!!!!!!!!!! CONJ
  }

  int Neps = 2000;
  double* eps = new double[Neps];
  for(int i = 0; i < Neps; i++) 
    eps[i] = -2.0*t + i * 4.0 * t / ((double)(Neps-1));
  complex<double>* g = new complex<double>[Neps];
  complex<double>* summand = new complex<double>[2*N];
  complex<double> sum = 0.0;
  for (int m=0; m<2*N-n; m++)
  { 
    //#pragma omp parallel for
    for(int i =0; i<Neps; i++)
    {  g[i] = DOS(DOStypes::SemiCircle, t, eps[i]) 
              * ( (v!=NULL) ? v(eps[i]) : 1.0 ) 
              * 2.0 
              * ( 1.0
                  / ( ( ii*iw_large[m] + mu - eps[i] - Sig_large[m] ) 
                       * 
                      ( ii*iw_large[m+n] + mu - eps[i] - Sig_large[m+n] ) 
                    )
                );
    }
//    char integFN[300];
//    sprintf(integFN,"integ.n%d.m%d",n,m);
//    if ((m<10)and(n<10)) PrintFunc(integFN,Neps,g,eps);  

    summand[m] = TrapezIntegral(Neps, g, eps);
    sum += summand[m]; //(1.0/nu) * Temp 
  }

  delete [] iw_large;
  delete [] Sig_large;
  delete [] summand;
  delete [] eps;
  delete [] g;

  return sum;

/*  char summandFN[200];
  sprintf(summandFN,"summand.n%.4d",n);
  PrintFunc(summandFN,N,summand,omega); 
*/

}



//---------------------------------------------------------------------------//

complex<double> IAResult::InternalEnergy(double T, IBZ* ibz)
{
  int N = iagrid->get_N();

  double* iw_large = new double[2*N];
  complex<double>* Sig_large = new complex<double>[2*N];
  for(int i = 0; i < N; i++)
  {
    iw_large[N+i]    =  omega[i];
    iw_large[N-1-i]  = -omega[i];        
    Sig_large[N+i]   =  Sigma[i];
    Sig_large[N-1-i] = conj(Sigma[i]);
  }

  complex<double> sum = 0.0;
  int Nsum = 5000;
  //for (int m=0; m<2*N-n; m++)
  for (int m=N-1-Nsum; m<N+Nsum; m++)
  { 
    for(int i=0; i<ibz->get_Nx(); i++)
    for(int j=0; j<ibz->get_Ny(); j++)
    { 
      double eps = ibz->epsilon[i][j];
      ibz->summand[i][j] =            (eps + ii*iw_large[m])  
                           / ( ii*iw_large[m] + mu - eps - Sig_large[m] )  
                           - 1.0; 
                  
    }
    sum += T*ibz->sum();
    
  }

  delete [] iw_large;
  delete [] Sig_large;

  return sum;
}



complex<double> IAResult::InternalEnergy(double T)
{
  double t = 0.5; //!!!!!!!!!!!!!!!
  int N = iagrid->get_N();

  double* iw_large = new double[2*N];
  complex<double>* Sig_large = new complex<double>[2*N];
  for(int i = 0; i < N; i++)
  {
    iw_large[N+i]    =  omega[i];
    iw_large[N-1-i]  = -omega[i];        
    Sig_large[N+i]   =  Sigma[i];
    Sig_large[N-1-i] = conj(Sigma[i]);
  }

  int Neps = 2000;
  double* eps = new double[Neps];
  for(int i = 0; i < Neps; i++) 
    eps[i] = -2.0*t + i * 4.0 * t / ((double)(Neps-1));
 
  complex<double> sum = 0.0;
  complex<double>* g = new complex<double>[2*N];
  for (int m=0; m<2*N; m++)
  { for(int i = 0; i<Neps; i++)
    {  g[i] =  DOS(DOStypes::SemiCircle, t, eps[i]) 
               *
               (                  (eps[i]+ii*iw_large[m])
                  / ( ii*iw_large[m] + mu - eps[i] - Sig_large[m] ) 
                  -1.0 
               );
    }   
    sum += TrapezIntegral(Neps, g, eps);; 
  }
  return T*sum + mu*n;


  delete [] iw_large;
  delete [] Sig_large;
  delete [] eps;
  delete [] g;
}

complex<double> IAResult::SmartInternalEnergy(double T)
{
  double t = 0.5; //!!!!!!!!!!!!!!!
  int N = iagrid->get_N();

  double* iw_large = new double[2*N];
  complex<double>* G_large = new complex<double>[2*N];
  complex<double>* Delta_large = new complex<double>[2*N];
  for(int i = 0; i < N; i++)
  {
    iw_large[N+i]    =  omega[i];
    iw_large[N-1-i]  = -omega[i];        
    G_large[N+i]   =  G[i];
    G_large[N-1-i] = conj(G[i]);
    Delta_large[N+i]   =  Delta[i];
    Delta_large[N-1-i] = conj(Delta[i]);
  }

  complex<double> sum = 0.0;
  for (int m=0; m<2*N; m++)
    sum += (Delta_large[m]+ii*iw_large[m])*G_large[m] - 1.0;

  delete [] iw_large;
  delete [] G_large;
  delete [] Delta_large;

  
  return T*sum + mu*n;
}

void IAResult::PatchAtomicLimitSigma(double AtomicCutoff, double U)
{
  for (int i=0; i<iagrid->get_N(); i++) 
    if (omega[i]>AtomicCutoff) Sigma[i] = AtomicLimitSigma(omega[i], mu, n, U);
}

void IAResult::PatchAtomicLimitG(double AtomicCutoff, double U)
{
  for (int i=0; i<iagrid->get_N(); i++) 
    if (omega[i]>AtomicCutoff) G[i] = AtomicLimitG(omega[i], mu, n, U);
}

void IAResult::PatchDelta(double AtomicCutoff)
{
  int i;
  for (i=0; i<iagrid->get_N(); i++) 
    if (omega[i]>AtomicCutoff) break;
  
  double A = real(Delta[i])*sqr(omega[i]);
  double B = imag(Delta[i])*omega[i];

  for (int j=i+1; j<iagrid->get_N(); j++) 
    Delta[j] = complex<double>(A/sqr(omega[j]),B/omega[j]);
}

/*
void IAResult::InitG(complex<double>* G, double a, double t, int type)
{
  switch (type)
  {  case 0:
     {
       for(int n=0; n<N; n++)
       {  complex<double> z = complex<double>(a,get_omega(n));
          complex<double> sq  = sqrt( complex<double>( z*z - 4.0*sqr(t) ) );
          double sign = abs( imag(sq) ) / imag(sq) ;
          complex<double> G0 = ( z - sign*sq ) / (2.0*sqr(t)) ;
          G[n] =  1.0/( ii*get_omega(n) - sqr(t)*G0 );
       } 
       break;  
     } 
     case 1:
     {
       for(int n=0; n<N; n++) G[n] = 0.5/(ii*get_omega(n)+5.0) + 0.5/(ii*get_omega(n)-5.0);
       break;  
     } 
     default:
     {  
       printf("Initial G type not implemented!\n");
       exit(1);
     }
  }
} 


void IAGRID::InitDelta(complex<double>* Delta, double a, double t, int type)
{
  InitG(Delta,a,t,type);
  for(int n=0; n<N; n++) Delta[n] = ii*get_omega(n) - 1.0/Delta[n];  
}
*/
/*
double IAResult::TriangularConductivity(double T, int Nkx, int Nky, int Nnu, const char * integrandFN) //bool excludeSmallOmega)
//--------------------- DC CONDUCTIVITY --------------------------// Neps,Nnu ~ 400 or 800
{
  FILE* integrandFile;
  if (integrandFN!=NULL)
    integrandFile = fopen(integrandFN,"w");


  IBZ ibz(IBZtypes::TriangularLattice, Nkx, Nky );
  printf("-----ibz ready\n");
  complex<double> sum = 0.0;
  
  double k = 20.0;		//determines the nu range

  double nu_start = -k*T;
  double nu_end = k*T;

  double dnu = (nu_end-nu_start)/ (double) Nnu;
 
  for (double nu = nu_start; nu < nu_end; nu += dnu )
  { printf("-- nu: %.3f\n",nu);
    //---------- ONLY FOR INSULATOR -----------//
    //if ((abs(nu)<0.1)and(excludeSmallOmega))
    //  continue;
    //-----------------------------------------//

    complex<double> Sigma_nu=grid->interpl(Sigma,nu);
    double fprim = 1.0 / ( 4.0 * T * sqr( cosh( nu/(2.0*T) 
                                              ) 
                                        )   
                           );

    for(int i=0; i<Nkx; i++)
    for(int j=0; j<Nky; j++)
    { 
      double v = ibz.velocity[i][j];
 
      complex<double> G = 1.0/( nu + mu - ibz.epsilon[i][j] - Sigma_nu);
      double rho = - imag(G) / pi;

      ibz.summand[i][j] = sqr(rho) * sqr(v) * fprim;

      //if (integrandFN!=NULL) fprintf(integrandFile,"%.15le %.15le %.15le %.15le %.15le %.15le\n", eps, nu, integrand, rho, fprim, rho0);  
    }
    sum+=ibz.sum();

    if (integrandFN!=NULL) fprintf(integrandFile,"\n");  
  }
  if (integrandFN!=NULL) fclose(integrandFile);

  return 2.0 * pi * real(sum) * dnu ;
}
*/




/*
double IAResult::Conductivity(double w, double T, double mu, double Neps, double Nnu, double Nscan)
{
//  double n = 100.0;
//  double k = 10.0;
  double sum = 0.0;
  //char integrandFN[100];
  //sprintf(integrandFN, "integrand.U%.3f.T%.3f",2.0*mu,T);
  //FILE* integrandFile = fopen(integrandFN,"w");
  double k=5.0;
  double nu = k*T/2.0;
  double eps_max = 0.0;
  double rho_max = 0.0;
  for (double eps = -1.0; eps < 1.0; eps += 2.0 / Nscan )
  {
    complex<double> G = 1.0/( nu + mu - eps - grid->interpl(Sigma,nu) );
    double rho = - imag(G) / pi;
    
    if (rho_max < rho) { rho_max = rho; eps_max = eps; }
  }
  double halfWidth=0;
  for (double eps = -1.0; eps < 1.0; eps += 2.0 / Nscan )
  {
    complex<double> G = 1.0/( nu + mu - eps - grid->interpl(Sigma,nu) );
    double rho = - imag(G) / pi;
    
    if (rho_max < 100*rho) { halfWidth = abs(eps_max-eps); break; } 
  }

  double angle = nu/eps_max;

  if (angle<0) k = 12.0; 
  printf("rho_max: %f eps_max: %f angle: %f halfWidth: %f\n",rho_max, eps_max, angle, halfWidth);

  k=10.0; 
  double dnu = 2.0 * k * T / Nnu;
  for (double nu = -k*T; nu < k*T; nu += dnu )
  { 
    double eps_start = -1.0;
    double eps_end = 1.0;
    //double eps_start = nu/angle - halfWidth;
    //double eps_end = nu/angle + halfWidth;
    //if ((eps_start < -1.0)or(angle<0)) eps_start = -1.0;
    //if ((eps_end > 1.0)or(angle<0)) eps_end = 1.0;
    double deps = (eps_end-eps_start) / Neps ;
    for (double eps = eps_start; eps < eps_end; eps += deps )
    { 
      double v = sqrt( (1.0-sqr(eps))/3.0 );
      double rho0 = grid->interpl(NIDOS,eps);
 
      complex<double> G = 1.0/( nu + mu - eps - grid->interpl(Sigma,nu) );
      complex<double> Gw = 1.0/( nu + w + mu - eps - grid->interpl(Sigma,nu) );
      double rho = - imag(G) / pi;
      double rhow = - imag(Gw) / pi;
      double fnu = 1.0 / (1.0 + exp(nu/T));
      double fnuw = 1.0 / (1.0 + exp((nu+w)/T));
      double integrand = rho0 * rho * rhow * sqr(v) * (fnu-fnuw)/w;
      sum += integrand * deps;
      //fprintf(integrandFile,"%.15le %.15le %.15le %.15le %.15le\n", eps, nu, integrand, rho, (fnu-fnuw)/w);  
    }
    //fprintf(integrandFile,"\n");  
  }
  //fclose(integrandFile);

  return 2.0 * pi * sum * dnu ;
}*/



