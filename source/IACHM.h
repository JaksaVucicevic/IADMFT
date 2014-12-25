#include <iostream>
#include "Loop.h"
#include "FFT.h"
using namespace std;

class IASIAM;

class IACHM: public Loop
{
  protected:
    void Defaults();
    string ParamsFN;

    double U;
    double T;
    double t;   //needed only when Bethe self-consistency is used (see below)

    IASIAM* iasiam;


    bool UseBethe; //Delta=t^2 G
                   //... or  ..
    int Nw;        //fill in the non-interacting DOS. and then G_latt= int deps rho_0(eps) (iw+mu-eps-Sigma)^{-1}
    double* w;     //eps
    double* NIDOS; //rho_0

    int SiamNt;
 
    virtual bool SolveSIAM();
    virtual void CalcDelta();  
   
    virtual void ReleaseMemory();

  public:
    IACHM();
    IACHM(const char* ParamsFN);   
    ~IACHM();

    FFT fft;
  
    void SetParams(double U, double T, double t);
    void SetNIDOS(int Nw, double* NIDOS, double* w);
    bool UseFixedMuSIAMRun;
    bool PHSymmetricCase;

    void SetUseBethe(bool UseBethe);

    double get_U() { return U; };
    double get_T() { return T; };
};


