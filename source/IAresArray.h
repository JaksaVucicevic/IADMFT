#include <complex>
using namespace std;

class IAResult;
class IAGRID;

class IAresArray
{
  public:
    IAresArray(int N, IAGRID* g);
    IAresArray(IAresArray &a); 
    ~IAresArray();

    void Set_n(double n);
    void Set_n0(double n0);
    void Set_mu(double mu);
    void Set_mu0(double mu0);
    double Global_n();
    
    IAResult* r;
    complex<double>* totalDelta;
    void WriteTotalDelta();
    void ReadTotalDelta();

    void PrintAll(const char* bareFN);
    void PrintAllShort(const char* bareFN);
    bool ReadFromFiles(const char* bareFN);
    void CopyFrom(IAresArray &a);
 
    int get_N() {return N;};
  private:
    int N;
};


