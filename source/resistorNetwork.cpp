#include<sys/stat.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include<complex>
// #include <algorithm>
#include <mkl.h>
#include "resistorNetwork.h"
#include "routines.h"

using namespace std;

void PrintMatrix(double* C, long long int N)
{
  FILE *f=fopen("Matrica.dat", "w");
  for (long long int i=0; i<N; ++i)
  {
    for (long long int j=0; j<N; ++j)
      fprintf(f,"%lf\n", C[i*N+j]);
    printf("%s%5.5lf   ", "", C[i*(N+1)]);
    if (i%20 == 0) printf("\n");
  }
  fclose(f);
  printf("\n");
}

void LUDijag(double* A, long long int N, double *Aii)
{
  long long int INFO=0;
  char UL = 'L'; //  !!!!!!!! mkl routine accepts transposed matrix (U is L)
  dpotrf_(&UL, &N, A, &N, &INFO);
  for (long long int i=0; i<N; ++i)
    Aii[i] = A[i*(N+1)];
}

double Det(double* A, long long int N)
{
  double *M;
  M = new double [N*N];
  for (long long int i=0; i<N; ++i)
  for (long long int j=i; j<N; ++j)
    M[i*N+j] = A[i*N+j];
  long long int INFO=0;
  char UL = 'L'; //  !!!!!!!! mkl routine accepts transposed matrix (U is L)
  double det=1.0;
  dpotrf_(&UL, &N, M, &N, &INFO);
  for (long long int i=0; i<N; ++i)
    det *= M[i*(N+1)];
  delete [] M;
  return sqr(det);
}

void Minor_in(long long int in, double* Cout, double* Cin, long long int N)
{
  long long int N1, k;
  N1 = N - 1;
  k=0;
  for (long long int i=0; i<N; ++i)
  {
    if (i != in)
    {
/*-----Half or full matrix-----*/
      if (i > in)
        k--;
      for (long long int j=i*(N+1); j<(i+1)*N; ++j)
//       for (long long int j=i*N; j<(i+1)*N; j++)
/*-----------------------------*/
      {
        if (j-i*N != in)
          Cout[k+j] = Cin[j];
        else
          k--;
      }
    }
    else
      k -= N;
  }
}

double RhoEquiv(long long int in, long long int out, double* C, long long int N)
{
  if (in == out)
  {
    printf("Cvorovi moraju da budu razliciti!\n");
    exit(1);
  }
  long long int N1, N2;
  N1 = N-1;
  N2 = N-2;
  double *Ci, *Ci_dij;
  double *Cij, *Cij_dij;
  double Req=1.0;
  Ci = new double [N1*N1];
  Cij = new double [N2*N2];
  Ci_dij = new double [N1];
  Cij_dij = new double [N2];
  if (in < out)
    out--;
  Minor_in(in, Ci, C, N);
  Minor_in(out, Cij, Ci, N1);
  LUDijag(Ci, N1, Ci_dij);
  LUDijag(Cij, N2, Cij_dij);
  for (long long int i=0; i<N2; ++i)
    Req *= Cij_dij[i]/Ci_dij[i];
  Req /= Ci_dij[N1-1];
  Req = sqr(Req);
  delete [] Ci;
  delete [] Cij;
  delete [] Ci_dij;
  delete [] Cij_dij;
  return Req;
}

void neighbors(long long int i, vector<long long int>& j, long long int Nx, long long int Ny, long long int Nz, long long int BoundaryX, long long int BoundaryY, long long int BoundaryZ, bool SelfNeighbor)
// j is vector of neighbours of i
{
  bool   nexist;
  long long int nx, nxy, nxyz, Nxyz, Nxyz2;
  nx = (Nx-1)*BoundaryX; nxy = Nx*(Ny-1)*BoundaryY; nxyz = Nx*Ny*(Nz-1)*BoundaryZ;
  Nxyz=Nx*Ny*Nz;
  Nxyz2=Nxyz*Nxyz;
  long long int    narr[6];       // array of indexes movement to obtain nearest neighbours [xm, xp, ym, yp, zm, zp]
  narr[0] = -1;  narr[1] = 1;  narr[2] = -Nx;  narr[3] = Nx;  narr[4] = -Nx*Ny;  narr[5] = Nx*Ny;
  if (i < Nx*Ny) narr[4] = nxyz;                //    else
  if (i >= Nx*Ny*(Nz-1)) narr[5] = -nxyz;
  if ((i/Nx) % Ny == 0) narr[2] = nxy;          //    else
  if ((i/Nx) % Ny == Ny-1) narr[3] = -nxy;
  if (i % Nx == 0) narr[0] = nx;                //    else
  if (i % Nx == Nx-1) narr[1] = -nx;
  for (long long int br=5; br>=0; br--)
  {
    nexist = true;
    for (long long int ibr = 0; ibr < j.size(); ibr++)
      if ((narr[br]+i == j[ibr]) ) nexist = false;
    if ((not SelfNeighbor) && (narr[br] == 0)) nexist = false;
    if (nexist)
      j.push_back(narr[br]+i);
  }
}

bool comp (long long int i, long long int j) { return (i>j); }

bool elem (long long int a, long long int *La, long long int NLa)
{
  for (long long int i=0; i<NLa; i++) if (a==La[i]) return true;
  return false;
}

//=================================================================================================================//

double Resistivity(int Direction, //0: x, 1: y, 2: z, 3: average
                   long long int Nx, 
                   long long int Ny,
                   long long int Nz,
                      long long int BoundaryX,
                      long long int BoundaryY,
                      long long int BoundaryZ,
                   double* rhos,
                   bool SelfNeighbor)
{
  if (Direction==3)
    return (   Resistivity(0, Nx, Ny, Nz, BoundaryX, BoundaryY, BoundaryZ, rhos, SelfNeighbor)
             + Resistivity(1, Nx, Ny, Nz, BoundaryX, BoundaryY, BoundaryZ, rhos, SelfNeighbor)
             + Resistivity(2, Nx, Ny, Nz, BoundaryX, BoundaryY, BoundaryZ, rhos, SelfNeighbor)
           )/3.0;
   
  long long int N = Nx*Ny*Nz;
  long long int NN = N*N;

  long long int NLead;

  switch (Direction)
  { //x
    case 0: NLead = Ny*Nz;
            break;
    //y
    case 1: NLead = Nx*Nz;
            break;
    //z
    case 2: NLead = Nx*Ny;
            break;
    default: printf("Direction not implemented!\n"); 
             return 0;
  }

  long long int* inList = new long long int[NLead];
  long long int* outList = new long long int[NLead];

  int counter = 0;
  switch (Direction)
  { //x
    case 0: for (long long int i=0; i<NLead; i++)
            { inList[i]=i*Nx;
              outList[i]=inList[i]+Nx/2;
              //printf("%.2d---- in: %lld out: %lld\n", counter++, inList[i], outList[i]);
            }
            break;
    //y
    case 1: for (long long int i=0; i<N; i+=Nx*Ny)
              for (long long int j=0; j<Nx; j++) 
              { inList[counter]=i+j;
                outList[counter]=inList[counter]+Nx*Ny/2;
                //printf("%.2d---- in: %lld out: %lld\n", counter, inList[counter], outList[counter]);
                counter++;
              }
            break;
    //z
    case 2: for (long long int i=0; i<NLead; i++)
            { inList[i]=i;
              outList[i]=inList[i]+N/2;
              //printf("%.2d---- in: %lld out: %lld\n", counter++, inList[i], outList[i]);
            }
            break;
    default: printf("Direction not implemented!\n"); 
             return 0;
  }

  double* A = new double[N];
  for (long long int i=0; i<N; ++i)
    A[i] = rhos[i];
  
  vector<long long int> *Neighbors;
  Neighbors = new vector<long long int>[N];
  for (long long int i=0; i<N; i++)
    neighbors(i, Neighbors[i], Nx, Ny, Nz, BoundaryX, BoundaryY, BoundaryZ, SelfNeighbor);

  double *CND;
  double *CNDc;
  CND = new double[NN];
  CNDc = new double[NN];
  for (long long int i=0; i<N; ++i)
  {
    CND[i*(N+1)] = 0.0;
    CNDc[i*(N+1)] = Neighbors[i].size();
    for (long long int j=0; j<Neighbors[i].size(); ++j)
    {
      CND[i*(N+1)] += 2.0/(A[i]+A[Neighbors[i][j]]);
      if (Neighbors[i][j]!=i)
      {
        CND[i*N+Neighbors[i][j]]= -2.0/(A[i]+A[Neighbors[i][j]]);
        CNDc[i*N+Neighbors[i][j]]= -1.0;
      }
    }
  }

// /*  Conductance between two leads
// ----------------------------------------- -------- -----------------------------------//
  double *CNDr, *CNDcr;
  long long int Nr=N-2*NLead+2, NNr=Nr*Nr;
  if (NLead > 1)
  {
    CNDr = new double[NNr];
    CNDcr = new double[NNr];
    for (long long int i=0; i<NNr; i++)
    {
      CNDr[i] = 0.0;
      CNDcr[i] = 0.0;
    }
    long long int *Index=new long long int[N];
    for (long long int i=0; i<N; ++i) Index[i]=-1;
    for (long long int i=0; i<NLead; ++i) Index[inList[i]] = 0;
    for (long long int i=0; i<NLead; ++i) Index[outList[i]] = Nr-1;
    long long int li=0;
    for (long long int i=0; i<N; ++i)
      if (Index[i]==-1)
      {
        ++li;
        Index[i] = li;
      }
    if (li != Nr-2) { printf("\n\nLOSE PREBROJANO!!!\nli = %lld\nNr-2 = %lld\n\n", li, Nr-2); exit(1); }
//     else printf("\n\nDOBRO PREBROJANO!!!\nli = %d\nNr-3 = %d\n\n", li, Nr-2);
    for (long long int i=0; i<N; ++i)
    {
      for (long long int j=0; j<Neighbors[i].size(); ++j)
      {
          CNDr[Index[i]*Nr+Index[Neighbors[i][j]]] += CND[i*N+Neighbors[i][j]];
          CNDcr[Index[i]*Nr+Index[Neighbors[i][j]]] += CNDc[i*N+Neighbors[i][j]];
      }
      CNDr[Index[i]*(Nr+1)] += CND[i*(N+1)];
      CNDcr[Index[i]*(Nr+1)] += CNDc[i*(N+1)];
    }
    delete [] Index;
  }
  else
  {
    printf("Provodnost izmedju dva cvora:\n\n");
    CNDr = CND;
    CNDcr = CNDc;
  }

  double Rho, NormRho, RhoN;
  NormRho = RhoEquiv(0, Nr-1, CNDcr, Nr);
//   PrintMatrix(CNDr, Nr);
  Rho = RhoEquiv(0, Nr-1, CNDr, Nr);
  RhoN = Rho/NormRho;
  printf("\n\nNormRho = %le\n\n", NormRho);
  printf("\n\nRhoNonNorm = %le\n\n", Rho);
  printf("\n\nRho_eff = %le\n\n", Rho/NormRho);

  if (NLead > 1)
  {
    delete [] CNDr;
    delete [] CNDcr;
  }

  delete [] Neighbors;
  delete [] A;
  delete [] CND;
  delete [] CNDc;


  return Rho/NormRho;
}
