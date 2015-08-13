#include "ctqmcDHM.h"


class ctqmcCPA: public ctqmcDHM
{
  public:
    ctqmcCPA();
    ~ctqmcCPA(); 

    double t;
    bool UseBethe;
    bool UseAverageSigma;
    bool PHSymmetry;

    int Nw;        //fill in the non-interacting DOS. and then G_latt= int deps rho_0(eps) (iw+mu-eps-Sigma)^{-1}
    double* w;     //eps
    double* NIDOS; //rho_0

    void DiscretizeMus(IAresArray &a, double mu);
    double mu; //this is necessary for CalcDelta!

  protected:
    virtual void CalcDelta();

};
