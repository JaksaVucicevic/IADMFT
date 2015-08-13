#include <iostream>
#include "IACHM.h"
#include "Broyden.h"

class ctqmcSIAM;


using namespace std;

class ctqmcCHM: public IACHM
{
  private:

    int LastSign;

  public:
    ctqmcCHM();
    ~ctqmcCHM();

    ctqmcSIAM* ctqmcsiam;

    double c;

    bool PrintMuHistory;
    bool UseBroydenForMu;
    bool UseSmart_c;
    bool UseSmart_M;
    int Smart_it;
    unsigned long int Smart_smallM;
    unsigned long int Smart_largeM; 

    double BroydenStartDiff;
    int BroydenStatus;
    Broyden B;

    bool UseIPT;




  protected:
    virtual bool SolveSIAM();

};
