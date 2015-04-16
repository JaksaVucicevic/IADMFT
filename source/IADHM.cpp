#include <cstdio>
#include <cstdlib>
#include "IADHM.h"
#include "IAGRID.h"
#include "IAResult.h"
#include "IAresArray.h"
#include "Input.h"
#include "IASIAM.h"
#include "routines.h"
#include "arrayInitializers.h"
#include "mkl_invertMatrix.h"

#ifdef _OMP
#include <omp.h>
#endif

class IASIAM;

void IADHM::Defaults() 
{
    U = 2.0;
    W = 1.0; 
    T = 0.05;

    //UseBethe = true;
    //PHSymmetricCase = true;
}

IADHM::IADHM() : Loop()
{
  Defaults();
}

IADHM::IADHM(const char* ParamsFN) : Loop(ParamsFN)
{
  Defaults();
  this->ParamsFN.assign(ParamsFN);
  cout << "-- INFO -- CHM: Params File Name set to:" << this->ParamsFN << endl;

  Input input(ParamsFN);
  
  input.ReadParam(U,"IADHM::U");
  input.ReadParam(T,"IADHM::T");
  input.ReadParam(W,"IADHM::W");
//  input.ReadParam(UseBethe,"IACHM::UseBethe");
//  input.ReadParam(UseFixedMuSIAMRun,"IACHM::UseFixedMuSIAMRun");

}

void IADHM::ReleaseMemory()
{
  printf("CHM release\n");
}

IADHM::~IADHM()
{
  ReleaseMemory();
}

void IADHM::SetParams(int Nsites, double U, double W, double T)
{
  this->N = Nsites;
  this->U = U;
  this->T = T;
  this->W = W;
}
/*
void IADHM::SetNIDOS(int Nw, double* NIDOS, double* w)
{
  this->Nw = Nw;
  this->NIDOS = NIDOS;
  this->w = w;
}

void IADHM::SetUseBethe(bool UseBethe)
{
  this->UseBethe = UseBethe;
}
*/
//-------------------------------------------------------//
  
void IADHM::AssignMus( IAresArray &a, double mu, unsigned int seed )
{
  srand (seed);

  for(int id=0; id<N; id++)
    a.r[id].mu = mu - W/2.0 + W * ((double) rand()) / ((double) RAND_MAX);
}

void IADHM::GetMus( IAresArray &a, double* mus )
{
  for(int id=0; id<N; id++)
    mus[id] = a.r[id].mu;
}


bool IADHM::SolveSIAM()
{
  iasiams = new IASIAM[N];
  #pragma omp parallel for
  for(int id=0; id<N; id++)
  { iasiams[id].SetUTepsilon(U,T,0.0);
    //iasiam->PHSymmetricCase = PHSymmetricCase;
    iasiams[id].fft = &(fft[id]);
    iasiams[id].Run(&(a->r[id]));
  }
  delete [] iasiams;
  return false;
}

void IADHM::CalcDelta()
{
  int Nfreq = a->r[0].iagrid->get_N();
  #pragma omp parallel for 
  for(int n=0; n<Nfreq; n++)
  { 
    printf("IADMFT::CalcDelta: working iw_%d\n",n);
    complex<double>** invG = Array2D< complex<double> >(N, N);

    for(int i=0; i<N; i++)
    for(int j=0; j<N; j++)
      invG[i][j] = H0[i][j];

    for(int id=0; id<N; id++)
      invG[id][id] += ii*a->r[id].omega[n] + a->r[id].mu - a->r[id].Sigma[n];   

    //if(n==0) PrintMatrix("invG", N, N, invG);

    complex<double>** G = Array2D< complex<double> >(N, N);   
    InvertSymmetricMatrix(N, invG, G);

    //if(n==0) PrintMatrix("G", N, N, G);

    for(int id=0; id<N; id++)
      a->r[id].Delta[n] = ii*a->r[id].omega[n] + a->r[id].mu - a->r[id].Sigma[n] - 1.0/G[id][id];

    FreeArray2D< complex<double> >(invG, N);
    FreeArray2D< complex<double> >(G, N);
  }

/*
  if (PrintIntermediate)
  { char FN[300];
    sprintf( FN, "IADHM.U%.3f.T%.3f.it%d", U, T, Iteration);
    r->PrintResult(FN);
  }  
*/  
}
