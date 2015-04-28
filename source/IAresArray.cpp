#include "IAresArray.h"
#include "IAResult.h"
#include "routines.h"
#include "IAGRID.h"
#include <cstdio>

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

double IAresArray::Global_n()
{
  double n = 0;
  for(int id=0; id<N; id++)
    n += r[id].n;
  return n/((double)N);
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
    sprintf(FN, "%s%d", bareFN, id);
    r[id].PrintResult(FN);  
  }
}

void IAresArray::PrintAllShort(const char* bareFN)
{
  for(int id=0; id<N; id++)
  { char FN[300];
    sprintf(FN, "%s%d", bareFN, id);
    r[id].PrintShort(FN);  
  }
}

bool IAresArray::ReadFromFiles(const char* bareFN)
{
  for(int id=0; id<N; id++)
  { char FN[300];
    sprintf(FN, "%s%d", bareFN, id);
    if (not FileExists(FN)) return false;
    r[id].ReadFromFile(FN);  
  }
  return true;
}

void IAresArray::CopyFrom(IAresArray &a)
{
  for(int id=0; id<N; id++)
    r[id].CopyFrom(a.r[id]);  
}

