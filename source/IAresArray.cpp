#include "IAresArray.h"
#include "IAResult.h"
#include "routines.h"
#include "IAGRID.h"
#include "pade.h"
#include <cstdio>
#include <limits>

#ifdef _OMP
#include <omp.h>
#endif

using namespace std;

IAresArray::IAresArray(int N, IAGRID* g)
{
  this->N = N;
  
  r = new IAResult[N];
   
  for(int id=0; id<N; id++)
    r[id].Initialize(g);  

  totalDelta = new complex<double>[N*g->get_N()];
}

IAresArray::IAresArray(IAresArray &a)
{
  this->N = a.get_N();

  r = new IAResult[N];
   
  for(int id=0; id<N; id++)
  { r[id].Initialize(a.r[0].iagrid);  
    r[id].CopyFrom(a.r[id]); 
  }

  int Nw = a.r[0].iagrid->get_N();
  totalDelta = new complex<double>[Nw*N];
  for(int i=0; i<N*Nw; i++)
    totalDelta[i] = a.totalDelta[i];
}
 
IAresArray::~IAresArray()
{
  //this is NOT NEEDED, delete calls the destructor
  //for(int id=0; id<N; id++)
  //  r[id].~IAResult(); /
  delete [] r;
  delete [] totalDelta;
}

void IAresArray::Set_n(double n)
{
  for(int id=0; id<N; id++)
    r[id].n = n;
}
void IAresArray::Set_n0(double n0)
{
  for(int id=0; id<N; id++)
    r[id].n0 = n0;
}

void IAresArray::Set_mu(double mu)
{
  for(int id=0; id<N; id++)
    r[id].mu = mu;
}

void IAresArray::Set_mu0(double mu0)
{
  for(int id=0; id<N; id++)
    r[id].mu0 = mu0;
}
//---------------- get arrays -------------------//

double* IAresArray::Get_ns()
{
  double* ns = new double[N];
  for(int id=0; id<N; id++)
    ns[id] = r[id].n;

  return ns;
}

double* IAresArray::Get_n0s()
{
  double* n0s = new double[N];
  for(int id=0; id<N; id++)
    n0s[id] = r[id].n0;

  return n0s;
}

double* IAresArray::Get_mus()
{
  double* mus = new double[N];
  for(int id=0; id<N; id++)
    mus[id] = r[id].mu;

  return mus;
}

double* IAresArray::Get_mu0s()
{
  double* mu0s = new double[N];
  for(int id=0; id<N; id++)
    mu0s[id] = r[id].mu0;

  return mu0s;
}

double* IAresArray::Get_Aw0s()
{
  double* Aw0s = new double[N];
  for(int id=0; id<N; id++)
    Aw0s[id] = -(1.0/pi)*imag( CubicFrom4points(r[id].G,  r[id].omega) );

  return Aw0s;
}

double* IAresArray::Get_taus() //scattering rate
{
  double* taus = new double[N];
  for(int id=0; id<N; id++)
  {
    double inv_Sigma[4];
    for(int i=0; i<4; i++)
      inv_Sigma[i] = -1.0/imag(r[id].Sigma[i]); 
    taus[id] = CubicFrom4points(inv_Sigma,  r[id].omega);
  }
  return taus;
}

double* IAresArray::Get_Zs() //quasiparticle weight
{
  int M = GetActualN();
  if (M>50) M=50;
  double* Zs = new double[N];
  for(int id=0; id<N; id++)
  { 
    double ReSigw0 = real(pade( M, r[id].omega, r[id].Sigma, 0.0 ));
    double eta = 1e-1;
    double ReSigweta = real(pade( M, r[id].omega, r[id].Sigma, eta ));
   
    Zs[id] = 1.0/( 1.0 - (ReSigweta - ReSigw0)/eta );
    /*Zs[id] =  1.0 
              / ( 1.0 
                  - imag(r[id].Sigma[0])
                    / r[id].omega[0]
                 );*/
  }  


  return Zs;
}

double* IAresArray::Get_Zs_RASC(complex<double>* Delta) //quasiparticle weight
{
  int M = GetActualN();
  if (M>50) M=50;
  double* Zs = new double[N];
  for(int id=0; id<N; id++)
  { double eta = 1e-1;
    complex<double> Gw0 = pade( M, r[id].omega, r[id].G, 0.0 );
    complex<double> Deltaw0 = pade( M, r[id].omega, ((Delta!=NULL)?Delta:r[id].Delta), 0.0 );
    complex<double> Gweta = pade( M, r[id].omega, r[id].G, eta );
    complex<double> Deltaweta = pade( M, r[id].omega, ((Delta!=NULL)?Delta:r[id].Delta), eta );
    complex<double> Sigw0 = r[id].mu - Deltaw0 - 1.0/Gw0;
    complex<double> Sigweta = eta + r[id].mu - Deltaweta - 1.0/Gweta;
   
    Zs[id] = 1.0/( 1.0 - (real(Sigweta) - real(Sigw0))/eta );
  }  
  return Zs;
}

double* IAresArray::Get_taus_RASC(complex<double>* Delta) //quasiparticle weight
{
  int M = GetActualN();
  if (M>50) M=50;
  double* taus = new double[N];
  for(int id=0; id<N; id++)
  { complex<double> Gw0 = pade( M, r[id].omega, r[id].G, 0.0 );
    complex<double> Deltaw0 = pade( M, r[id].omega, ((Delta!=NULL)?Delta:r[id].Delta), 0.0 );
    complex<double> Sigw0 = r[id].mu - Deltaw0 - 1.0/Gw0;
   
    taus[id] = 1.0/imag(Sigw0);
  }  
  return taus;
}

double* IAresArray::Get_renormalized_mus() //renormalized mu
{
  double* rmus = new double[N];
  for(int id=0; id<N; id++)
  {
    double ReSigma[4];
    for(int i=0; i<4; i++)
      ReSigma[i] = real(r[id].Sigma[i]); 
    rmus[id] = -CubicFrom4points(ReSigma,  r[id].omega)+r[id].mu;
  }
  return rmus;
}

double* IAresArray::Get_ReSigw0s() //energy shift from self-energy
{
  double* s = new double[N];
  for(int id=0; id<N; id++)
  {
    double ReSigma[4];
    for(int i=0; i<4; i++)
      ReSigma[i] = real(r[id].Sigma[i]); 
    s[id] = CubicFrom4points(ReSigma,  r[id].omega);
  }
  return s;
}

double* IAresArray::Get_pade_Aw0s()
{
  double* Aw0s = new double[N];
  for(int id=0; id<N; id++)
  { int M = GetActualN();
    if (M>200) M=200;
    Aw0s[id] = -(1.0/pi)*imag( pade( M, r[id].omega, r[id].G, 0.0 ) );
  }

  return Aw0s;
}

double* IAresArray::Get_pade_taus() //scattering rate
{
  double* taus = new double[N];
  for(int id=0; id<N; id++)
  { int M = GetActualN();
    if (M>200) M=200;
    taus[id] = -1.0/imag( pade( M, r[id].omega, r[id].Sigma, 0.0 ) );
  }
  return taus;
}

//------------get sorted arrays------------------------//

/*
double* IAresArray::GetSorted_ns()
{
  IAResult** sorted_rs = GetSorted();
  double* ns = new double[N];
  for(int id=0; id<N; id++)
    ns[id] = sorted_rs[id]->n;
  delete [] sorted_rs;
  return ns; 
}

double* IAresArray::GetSorted_n0s()
{
  IAResult** sorted_rs = GetSorted();
  double* n0s = new double[N];
  for(int id=0; id<N; id++)
    n0s[id] = sorted_rs[id]->n0;
  delete [] sorted_rs;
  return n0s;
}

double* IAresArray::GetSorted_mus()
{
  IAResult** sorted_rs = GetSorted();
  double* mus = new double[N];
  for(int id=0; id<N; id++)
    mus[id] = sorted_rs[id]->mu;
  delete [] sorted_rs;
  return mus;
}

double* IAresArray::GetSorted_mu0s()
{
  IAResult** sorted_rs = GetSorted();
  double* mu0s = new double[N];
  for(int id=0; id<N; id++)
    mu0s[id] = sorted_rs[id]->mu0;
  delete [] sorted_rs;
  return mu0s;
}

double* IAresArray::GetSorted_Aw0s()
{
  IAResult** sorted_rs = GetSorted();
  double* Aw0s = new double[N];
  for(int id=0; id<N; id++)
    Aw0s[id] = -(1.0/pi)*imag( CubicFrom4points(sorted_rs[id]->G,  sorted_rs[id]->omega) );
  delete [] sorted_rs;
  return Aw0s;
}

double* IAresArray::GetSorted_taus() //scattering rate
{ 
  IAResult** sorted_rs = GetSorted(); 
  double* taus = new double[N];
  for(int id=0; id<N; id++)
  {
    double inv_Sigma[4];
    for(int i=0; i<4; i++)
      inv_Sigma[i] = -1.0/imag(sorted_rs[id]->Sigma[i]);
    taus[id] = CubicFrom4points(inv_Sigma, sorted_rs[id]->omega);
  }
  return taus;
}
*/

//-------------------------------------------------------//

int IAresArray::GetActualN()
{
  int M = r[0].iagrid->get_N();
  for(int i=0; i<r[0].iagrid->get_N(); i++)
    if (r[0].G[i] != r[0].G[i]) { M = i; break; }
  return M;
}

double IAresArray::Global_n()
{
  double n = 0;
  for(int id=0; id<N; id++)
    n += r[id].n;
  return n/((double)N);
}

void IAresArray::GetAverageGandSigma(complex<double>* G, complex<double>* Sigma)
{
  int M = GetActualN();

  #pragma omp parallel for
  for(int i=0; i<r[0].iagrid->get_N(); i++) 
  { 
    if (i>M)
    { G[i] = std::numeric_limits<double>::quiet_NaN();; 
      Sigma[i] = std::numeric_limits<double>::quiet_NaN();;
    }
    else
    {
      G[i]=0.0; 
      Sigma[i] = 0.0;
      int rejects = 0;
      for(int j=0; j<N; j++)
      {  if ( abs(r[j].mu0) > 5.0 ) { rejects++; continue; } //precautionary measure - if IASIAM failed to find a physical mu0, reject the result
         G[i]+=r[j].G[i]; 
         Sigma[i]+=r[j].Sigma[i]; 
      }
      G[i] /= (double) (N-rejects);
      Sigma[i] /= (double) (N-rejects);
    }
  }  
}

complex<double>* IAresArray::GetGtypical()
{
  int M = GetActualN();
 
  complex<double>* G = new complex<double>[r[0].iagrid->get_N()];
  #pragma omp parallel for
  for(int i=0; i<r[0].iagrid->get_N(); i++) 
  { 
    if (i>M)
      G[i] = std::numeric_limits<double>::quiet_NaN();
    else
    {
      double expo = 0.0;
      int rejects = 0;
      for(int j=0; j<N; j++)
      {  if ( abs(r[j].mu0) > 5.0 ) { rejects++; continue; } //precautionary measure - if IASIAM failed to find a physical mu0, reject the result
         expo+=log(-imag(r[j].G[i])); 
      }
      expo /= (double) (N-rejects);
      G[i] = exp(expo);
    }
  }  

  return G;
}

 
void IAresArray::WriteTotalDelta()
{
  int Nfreq = r[0].iagrid->get_N();
  for(int id=0; id<N; id++)
  for(int i=0; i<Nfreq; i++)
    totalDelta[id*Nfreq + i] = r[id].Delta[i];
}

void IAresArray::ReadTotalDelta()
{
  int Nfreq = r[0].iagrid->get_N();
  for(int id=0; id<N; id++)
  for(int i=0; i<Nfreq; i++)
    r[id].Delta[i] = totalDelta[id*Nfreq + i];
}
 
void IAresArray::PrintAll(const char* bareFN)
{
  for(int id=0; id<N; id++)
  { char FN[300];
    sprintf(FN, "%s/%d", bareFN, id);
    r[id].PrintResult(FN);  
  }
}

void IAresArray::PrintAllShort(const char* bareFN)
{
  for(int id=0; id<N; id++)
  { char FN[300];
    sprintf(FN, "%s/%d", bareFN, id);
    r[id].PrintShort(FN);  
  }
}

void IAresArray::PrintAllMinimal(const char* bareFN, double cutoff)
{
  for(int id=0; id<N; id++)
  { char FN[300];
    sprintf(FN, "%s/%d", bareFN, id);
    r[id].PrintMinimal(FN,cutoff);  
  }
}

void IAresArray::PrintAllUberMinimal(const char* bareFN, double cutoff)
{
  for(int id=0; id<N; id++)
  { char FN[300];
    sprintf(FN, "%s/%d", bareFN, id);
    r[id].PrintUberMinimal(FN,cutoff);  
  }
}

IAResult** IAresArray::GetSortedResults()
{
  // copy result array
  IAResult** sorted_rs = new IAResult*[N];
  for(int id=0; id<N; id++) 
    sorted_rs[id] = &(r[id]);
  
  // sort it
  for(int i=0; i<N-1; i++) 
  for(int j=i+1; j<N; j++) 
    if (sorted_rs[i]->mu>sorted_rs[j]->mu)
    { IAResult* temp = sorted_rs[i];
      sorted_rs[i] = sorted_rs[j];
      sorted_rs[j] = temp;
    }

  return sorted_rs;
}

double* IAresArray::GetSorted(double* X)
{
  // make a copy of X
  double* X_sorted = new double[N];
  for(int id=0; id<N; id++) 
    X_sorted[id] = X[id];
  
  //get mus
  double* mus = Get_mus();

  // sort mus and X simultaneously
  for(int i=0; i<N-1; i++) 
  for(int j=i+1; j<N; j++) 
    if (mus[i] > mus[j])
    { double temp = mus[i];
      mus[i] = mus[j];
      mus[j] = temp;

      temp = X_sorted[i];
      X_sorted[i] = X_sorted[j];
      X_sorted[j] = temp;
    }

  delete [] mus;

  return X_sorted;
}



void IAresArray::PrintSortedMinimal(const char* bareFN, double cutoff)
{
  // copy result array
  IAResult** sorted_rs = GetSortedResults();
  for(int i=0; i<N; i++)
  { char FN[300];
    sprintf(FN, "%s/mu%.5f", bareFN, sorted_rs[i]->mu);
    sorted_rs[i]->PrintMinimal(FN,cutoff);  
  }

  //print out in files named with mu and put the sorted mus in an array
  double* mus = Get_mus();
  double* mus_sorted = GetSorted(mus);
  char FN[300];
  sprintf(FN, "%s/mu_list", bareFN);
  PrintList(FN, N, mus_sorted); //prints it in one line, 5 decimal digits

  //print mus
  //later in gnuplot plot sorted with
  // list = "`echo $(ls mus)`"
  // splot for [i in list] 'mu'.i using 1:(i):7 w lp notitle

  //release memory   
  delete [] mus;
  delete [] mus_sorted;
  delete [] sorted_rs;
}

void IAresArray::Print_ns_and_mus(const char* bareFN)
{
    double* ns = Get_ns();
    double* mus = Get_mus();
    double* sorted_ns =  GetSorted(ns);
    double* sorted_mus = GetSorted(mus);

    char FN[300];
    sprintf(FN, "%s/ns", bareFN);
    PrintFunc(FN, N, ns, mus);
    
    sprintf(FN, "%s/ns_sorted", bareFN);
    PrintFunc(FN, N, sorted_ns, sorted_mus);

    delete [] sorted_ns;
    delete [] sorted_mus;
    delete [] ns;
    delete [] mus;
}

void IAresArray::Print_Histogram(const char* FN, double* X, bool Logarithmic)
{
  int Nbins = 12;
  double* x = new double[Nbins];
  double* P = new double[Nbins];
 
  Histogram(N, X, Nbins, x, P, Logarithmic);
  
  PrintFunc(FN, Nbins, P, x);

  delete [] x;
  delete [] P;
}

void IAresArray::Get_Histogram(double* X, int Nbins, double* &x, double* &P, bool Logarithmic)
{
  x = new double[Nbins];
  P = new double[Nbins];
 
  Histogram(N, X, Nbins, x, P, Logarithmic);
}


bool IAresArray::ReadFromFiles(const char* bareFN)
{
  for(int id=0; id<N; id++)
  { char FN[300];
    sprintf(FN, "%s/%d", bareFN, id);
    if (not FileExists(FN)) return false;
    r[id].ReadFromFile(FN);  
  }
  return true;
}

bool IAresArray::ReadFromMinimalFiles(const char* bareFN)
{
  for(int id=0; id<N; id++)
  { char FN[300];
    sprintf(FN, "%s/%d", bareFN, id);
    if (not FileExists(FN)) return false;
    r[id].ReadFromMinimal(FN);  
  }
  return true;
}

bool IAresArray::ReadFromUberMinimalFiles(const char* bareFN)
{
  for(int id=0; id<N; id++)
  { char FN[300];
    sprintf(FN, "%s/%d", bareFN, id);
    if (not FileExists(FN)) return false;
    r[id].ReadFromUberMinimal(FN);  
  }
  return true;
}

void IAresArray::CopyFrom(IAresArray &a)
{
  for(int id=0; id<N; id++)
    r[id].CopyFrom(a.r[id]);  
}

