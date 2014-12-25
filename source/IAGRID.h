#include <cstdio>
#include <complex>
#include <iostream>
#include <cmath>

using namespace std;

namespace IntialGtypes
{
  const int Metallic = 0;
  const int Insulating = 1;
}

class IAGRID
{ 
  private:  
    int N;
    int M;
    double T;
    double beta;

    double get_tau(int m);
    double get_omega(int n);

  public:
    IAGRID(int N, int M, double T);
    ~IAGRID();
    int get_N() { return N; };
    int get_M() { return M; };
     
    void assign_omega(double* omega);
    void assign_tau(double* tau);
};
    
   
