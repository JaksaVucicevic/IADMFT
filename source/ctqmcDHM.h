#include "IADHM.h"

class ctqmcSIAM;

using namespace std;

class ctqmcDHM: public IADHM
{
  public:
    ctqmcDHM();
    ~ctqmcDHM(); 
 
    ctqmcSIAM* ctqmcsiam;

    bool UseSeparateFolders; 

    bool UseIPT;

  protected:
    virtual bool SolveSIAM();

};
