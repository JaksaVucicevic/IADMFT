//**********************************************************//
//              IASIAM at Arbitrary Filling                   //
//                    by Jaksha Vuchichevicc                //
//**********************************************************//

#include <iostream>
#include <complex>

class IAResult;
class IAGRID;

using namespace std;

//======================= IASIAM Class ==========================================//

class IASIAM
{
  private:

    void Defaults();
    string ParamsFN;

    IAResult* r;

    //--impurity parameters--//
    double U;			//on-site repulsion
    double T;			//temperature
    double epsilon;		//impurity energy level
   
    //--bath parameters--// 
    double mu;			//global chemical potential
    double mu0;			//fictious chemical potential
    
    //--don't touch this---//
    bool HalfFilling;

    //-- Broyden solver options--//
    double Accr;
    int MAX_ITS; 

    //-- mu0 search --//
    bool UseBroydenFormu0;	//set this to false if only Amoeba method is tu be used
    int max_tries;		//number of tries (with different initial guesses) of broyden search before Amoeba is used
                                //set to a large number (say 100) if MPT corrections are used. Amoeba does only the mu0 search
    				
   
    //--storage arrays--//
    GRID* grid;
    int N;

     //--get functions--//
    double get_n(complex<double> X[]);

    //--get procedures--//
    void get_G0();
    void get_G0(complex<double>* V);
    void get_SOCSigma();
    double get_b();
    void get_Sigma();
    void get_G();
    void get_G(complex<double>* V); //used by broyden in solving systems of equations
    void get_G_CHM();
    void get_G_CHM(complex<double>* V); //used by broyden in solving systems of equations

    //--- IASIAM solver ---//
    void SolveSiam(complex<double>* V);
    void Amoeba(double accr, complex<double>* V);	//amoeba method for mu0 search.
    void Amoeba_CHM(double accr, complex<double>* V);

    double AmoebaScanStart;	//before amoeba starts, the equation is solved roughly 
                                //  (with accuracy AmobeScanStep) by scanning from AmoebaScanStart to AmoebaScanEnd.
    double AmoebaScanEnd; 	//make sure AmoebaScanStart and AmoebaScanEnd are sufficiently far apart (when U or W is large).
    double AmoebaScanStep;
    int AmoebaMaxIts;		//maximum number of Amoeba iterations
    bool AmoebaForceScanAndPrintOut;	//output n, n0, n-n0 and result when scanning for mu0 candidate
  
  public:
    //------ OPTIONS -------//
    void SetBroydenParameters(int MAX_ITS, double Accr);
    void SetT(double T);
    void SetU(double U);    
    void SetEpsilon(double epsilon);
    void SetAmoebaParams(double AmoebaScanStart, double AmoebaScanEnd, double AmoebaScanStep);

    //--Constructors/destructors--//
    IASIAM();  
    IASIAM(const char* ParamsFN);
    ~IASIAM();
       
    //--------RUN IASIAM--------//
    
    bool Run(IAResult* r); 
    bool Run_CHM(IAResult* r); 

  //---- FRIENDS -----//
  //function that will be calling private member functions in case of solving (systems of) equations
  friend bool UseBroyden(int, int, double, void (IASIAM::*)(complex<double>*), IASIAM*, complex<double>*);
   
};
