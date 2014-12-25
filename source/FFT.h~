/*****************************************************************/
//     Discrete Fourier Transform for Green's Functions          //
//                     with corrections                          //
//                            by Jaksa Vucicevic, January 2010   //
/*****************************************************************/

#include <complex>

using namespace std;

class IAGRID;

//-----------------Fast Fourier Transform Class---------------------//
//------------------------------------------------------------------//

class FFT
{ 
  private:
    bool Initialized;
    
    //-----general stuff-----
    IAGRID* iagrid; 
   
    double Temp;
    int N; 
    double* tau_m;
    double dtau;
    double* omega_n;
    double domega;
    

    double* a_n;
    double* A_n;
    double get_a_n(int n);
    double get_A_n(int n);

    //---little algebra----------
    complex<double>* ST_1;
    complex<double>* ST_3; 
    complex<double> Stage_1(int n);
    complex<double> Stage_3(int m);

    //-----F to T------
    double C;
    double Omega;
    double get_C(complex<double> G_f_max);
    double get_Omega(); 

    complex<double>* CT_FtoT_1;
    complex<double>* CT_FtoT_2;
    complex<double>* CT_FtoT_3;
    complex<double> CorrTerm_FtoT_1(int m);
    complex<double> CorrTerm_FtoT_2(int m);
    complex<double> CorrTerm_FtoT_3(int m);
    complex<double> Correction_FtoT(int m, complex<double> G_f_max); 
       
    complex<double> G_0;
    complex<double> G_beta;
  
    void Prepare_G_0_and_beta(complex<double> G_f[]);

    //-----T to F-----
    complex<double> get_dG(int n, complex<double> G_t[], int one_or_N);
    
    double* B_1_n;
    double* B_2_n;
    double* B_3_n;
  
    complex<double> CorrTerm_TtoF_1(int n, complex<double> G_t[], int one_or_N);
    complex<double> CorrTerm_TtoF_2(int n, complex<double> G_t[], int one_or_N);
    complex<double> Correction_TtoF(int n, complex<double> G_t[]);
    
    //-------FastFourier------------
    void PrepareData_FtoT(complex<double> G_f[], double data[]);
    void PrepareData_TtoF(complex<double> G_t[], double data[]);

    void FastFourier(double data[], unsigned long N, int isign);
    
  public:
    //----------user interface functions------------------------ 
    FFT();
    void Initialize(IAGRID* iagrid);
    void ReleaseMemory();
    void FtoT(complex<double> G_f[], complex<double> G_t[]);
    void TtoF(complex<double> G_t[], complex<double> G_f[]);
};








