#include <cstdlib>
#include "IAHFCHM.h"
#include "IAGRID.h"
#include "routines.h"
#include "omp.h"


//---------------------CALCULATION ROUTINES---------------------------//

void IAHFCHM::Init_G_f()
{ 
  iagrid->InitG(G_f,0.0,t,metal);
}

void IAHFCHM::get_G_0_f()
{
  if (UseBethe)
    for (int n=0; n<N; n++)
      G_0_f[n] = 1.0/(ii * omega[n] - sqr(t) * G_f[n]) ;
  else
    for (int n=0; n<N; n++)
      G_0_f[n] =  1.0 / 
                  ( Sigma_f[n] + 1.0/G_f[n] );
    
}

void IAHFCHM::get_Sigma_t()
{
  for (int m=0; m<N; m++)
    Sigma_t[m] = sqr(U) * G_t[m] * G_t[m] * G_t[m];
  
}

void IAHFCHM::get_G_f()
{
  if (UseBethe)
    for (int n=0; n<N; n++)
      G_f[n] = 1.0 / ( 1.0/G_0_f[n] - Sigma_f[n] );
  else
  {
    int Neps = 2000;
    double* eps = new double[Neps];
    for(int i =0; i<Neps; i++) 
      eps[i] = -2.0*t + i * 4.0 * t / ((double)(Neps-1));
    complex<double>* g = new complex<double>[Neps];

    for (int n=0; n<N; n++)
    { for(int i =0; i<Neps; i++)
        g[i] =   DOS(DOStypes::SemiCircle, t, eps[i]) 
               / ( ii*omega[n] - eps[i] - Sigma_f[n] );
      G_f[n] = TrapezIntegral(Neps, g, eps);
    }
    delete [] eps;
    delete [] g;
  }  
}


//----------------------------OUTPUT------------------------------//

void IAHFCHM::ShowStatus(int msg)
{ 
  if (msg == 0)
    cout << "Error: IAHFCHM not INITIALIZED!!!"<< endl;
  if (show_status) 
    switch (msg)
    {
      case 1: cout << "----------" <<
                      "IAHFCHM Run for T=" << Temp
                              << " U=" << U 
                   << "----------" << endl; break;
      case 2: cout << "Gf initialized!" << endl; break;
      case 3: cout << "Converged!!!" << endl; break;
      case 4: cout << "----------   Done!   ----------" << endl; break;
      case 5: cout << "----------   FFT done   ----------" << endl; break;

    }
}

void IAHFCHM::ShowStatus(int it, double accr)
{
  if (show_status)
    cout << "Iteration: " << it << ", accr: " << accr <<endl;
}

//----------------------------------------------------------------//
//                 USER INTERFACE FUNCTIONS                       //
//----------------------------------------------------------------//

void IAHFCHM::Initialize(int N, 
                     double Temp, 
                     double U, 
                     int metal, 
                     double t, 
                     double min_accr)
{
  this->N = N;
  this->Temp = Temp;  
  this->U = U;
  this->metal = metal;
  this->t = t;
  this->min_accr = min_accr;
  
  G_f = new complex<double> [N];
  G_0_f = new complex<double> [N];
  G_t = new complex<double> [N];
  Sigma_t = new complex<double> [N];
  Sigma_f = new complex<double> [N];
  for (int n=0; n<N; n++)
   Sigma_f[n] = 0.0; 
  
   
  iagrid = new IAGRID(N,N,Temp);
  iagrid->GetGrid(omega,tau);
  fft.Initialize(iagrid);
    
  Initialized = true;
}

void IAHFCHM::SetMetal(int metal)
{
  this->metal = metal;
}

void IAHFCHM::SetUT(double U, double Temp)
{
  this->Temp = Temp;  
  this->U = U;
 
  if (Initialized) iagrid->~IAGRID();
  iagrid = new IAGRID(N,N,Temp);
  iagrid->GetGrid(omega,tau);
  fft.Initialize(iagrid);
}

void IAHFCHM::SetOptions(bool print_final, 
                     int print_intermediate, 
                     bool print_initial, 
                     bool print_lambda, 
                     bool show_status)
{
  this->print_final = print_final;
  this->print_intermediate = print_intermediate;
  this->print_initial = print_initial;
  this->print_lambda = print_lambda;
  this->show_status = show_status;
}

double IAHFCHM::Run(int MAX_ITS, bool INIT_GF)
{ 
 
  if (!Initialized) 
  { 
    ShowStatus(0);
    exit(1); //if IAHFCHM not initialized, exit
  }
  ShowStatus(1);  
  
  //------initialize G_f------
  if (INIT_GF) Init_G_f();
  if (print_initial) 
    PrintFunc("gfin", N, G_f, omega);
  
  ShowStatus(2);
  
  //======== start iterations =========//
  int it;
  double* Sigma_f_old = new double[N];
  double* G_f_old = new double[N];
  double lndGw0_it5;
  double lndGw0_last;

  double lambda=0.0;

  FILE* lndGw0File;
  if (print_lambda)
  {
    char lndGw0FN[50];
    sprintf(lndGw0FN,"results/lndGw0.U%.3f.T%.3f",U,Temp);
    lndGw0File = fopen(lndGw0FN,"w");
  }

  for (it=1; it<=MAX_ITS; it++)
  { 

/*    for (int i=0; i<N; i++)
    {  Sigma_f_old[i] = imag(Sigma_f[i]);
       G_f_old[i] = imag(G_f[i]);
    }    
*/
    //-----------------------------------------// 
    get_G_0_f();  

    fft.FtoT(G_0_f, G_t);

   // ShowStatus(5);

    get_Sigma_t();
   
     //save old Sigma_f[1] to check convergence later 
   
    fft.TtoF(Sigma_t, Sigma_f); 

    //ShowStatus(5);

    for (int i=0; i<N; i++) Sigma_f[i]= complex<double>(0.0, imag(Sigma_f[i]));
    get_G_f();
    for (int i=0; i<N; i++) G_f[i]= complex<double>(0.0, imag(G_f[i]));
    /*
    if (it==1) 
      G_f_old = imag(G_f[0]);
    else
    {
      if (print_lambda)
      { double lndGfw0 = log(abs(imag(G_f[0])-G_f_old));
        if (it==5) lndGw0_it5 = lndGfw0;
        lndGw0_last = lndGfw0;
        fprintf(lndGw0File,"%d %.15le\n",it,lndGfw0);
      }
      G_f_old = imag(G_f[0]);
    }
    */
    //------------------------------------------//
  
    //--check convergence--
    /*double diffG = 0.0;
    double diffSigma = 0.0;
    for (int i=0; i<N; i++)
    {
      if ( diffSigma < abs( imag(Sigma_f[i]) - Sigma_f_old[i] ) )
        diffSigma = abs( imag(Sigma_f[i]) - Sigma_f_old[i] );
      if ( diffG < abs( imag(G_f[i]) - G_f_old[i] ) )
        diffG = abs( imag(G_f[i]) - G_f_old[i] );
    }
    double diff = (diffG>diffSigma) ? diffG : diffSigma;*/
 

    if ((print_lambda)and(it>1))
    { double lndGfw0 = log(abs(imag(G_f[0])-G_f_old[0]));
      ShowStatus(it, lndGfw0);
      fprintf(lndGw0File,"%d %.15le\n",it,lndGfw0);
    }
    for (int i=0; i<N; i++) G_f_old[i] = imag(G_f[i]);      
/*
    if (( diff < min_accr ) and (it>5))
    { 
      ShowStatus(3); 
     
      if (print_lambda)
      { 
        lambda = (lndGw0_last - lndGw0_it5)/(it-5);

        char lambdaFN[50];
        sprintf(lambdaFN,"lambda.T%.3f",Temp);
        FILE* lambdaFile = fopen(lambdaFN,"a");
        fprintf(lambdaFile,"%.15le %.15le\n",U,lambda);
        fclose(lambdaFile);
      }
     
      break; //if min_accr reached, exit the loop
    }
*/        
       
    //--intermediate print out--
    if (print_intermediate>0)
      if (it % print_intermediate == 0)
      {  
        char FN[200];

        sprintf(FN,"results/IAHFCHM.Gf.U%.3f.T%.5f.%s.it%.3d",U,Temp,(metal==0) ? "FromMet" : "FromIns", it);
        PrintFunc(FN,N,G_f,omega);

        sprintf(FN,"results/IAHFCHM.G0f.U%.3f.T%.5f.%s.it%.3d",U,Temp,(metal==0) ? "FromMet" : "FromIns", it);
        PrintFunc(FN,N,G_0_f,omega);

        sprintf(FN,"results/IAHFCHM.Gt.U%.3f.T%.5f.%s.it%.3d",U,Temp,(metal==0) ? "FromMet" : "FromIns", it);
        PrintFunc(FN,N,G_t,tau);

        sprintf(FN,"results/IAHFCHM.Sigmaf.U%.3f.T%.5f.%s.it%.3d",U,Temp,(metal==0) ? "FromMet" : "FromIns", it);
        PrintFunc(FN,N,Sigma_f,omega);

        sprintf(FN,"results/IAHFCHM.Sigmat.U%.3f.T%.5f.%s.it%.3d",U,Temp,(metal==0) ? "FromMet" : "FromIns", it);
        PrintFunc(FN,N,Sigma_t,tau);
      }
     
  } 

  if (print_lambda) fclose(lndGw0File);
 
  //========== end iterations ===========//

  if (print_final) 
  {  char FN[50];
     sprintf(FN,"results/IAHFCHM.Gf.U%.3f.T%.5f.%s",U,Temp,(metal==0) ? "FromMet" : "FromIns");
     PrintFunc(FN,N,G_f,omega);
     sprintf(FN,"results/IAHFCHM.Sigma.U%.3f.T%.5f.%s",U,Temp,(metal==0) ? "FromMet" : "FromIns");
     PrintFunc(FN,N,Sigma_f,omega);
  } 
  
  ShowStatus(4);

  delete [] Sigma_f_old;
  delete [] G_f_old;
  
  return lambda;
}


double IAHFCHM::DoubleOccupancy()
{
  double eta = 1e-5;
  complex<double> ii = complex<double>(0.0,1.0);
  complex<double> d = 0;
  for (int n = N-1; n>=0; n--)
    d += (ii * omega[n] * G_f[n] - 1.0) * exp( ii * omega[n] * eta) 
         - sqr(t) * G_f[n] * G_f[n] ;
  return real(d) * Temp / U;
}

complex<double> IAHFCHM::InternalEnergy()
{
  double eta = 1e-5;
  complex<double> ii = complex<double>(0.0,1.0);
  complex<double> E = 0;
  
  int Neps = 200;
  double* eps = new double[Neps];
  double* epsdos = new double[Neps];

  double deps = 2.0/((double) Neps - 1.0);
  for (int i=0; i<Neps; i++)
  {  eps[i] = -1.0 + i*deps;
     epsdos[i] = eps[i] * DOS(DOStypes::SemiCircle,0.5,eps[i]);
  }

  complex<double> sum = 0.0;
  for (int n = N-1; n>=0; n--)
  {
    complex<double>* intg = new complex<double>[Neps];
    complex<double>* negintg = new complex<double>[Neps];
    for (int i=0; i<Neps; i++)
    {  intg[i] = epsdos[i]
                 / (ii * omega[n] - Sigma_f[n] - eps[i]);

       negintg[i] = epsdos[i]
                    / (- ii * omega[n] + Sigma_f[n] - eps[i]);

    }
    complex<double> term1 = TrapezIntegral(Neps, intg, eps);
    complex<double> term2 = TrapezIntegral(Neps, negintg, eps);
    complex<double> term3 = Sigma_f[n] * G_f[n];
    sum += term1 + term2 + term3;
    delete [] intg;
    delete [] negintg;
  }
  delete [] eps;
  delete [] epsdos;
  return Temp*sum ;
}

double IAHFCHM::ChargeSusceptibility()
{
  double eta = 1e-4;
  complex<double> ii = complex<double>(0.0,1.0);
  double sum = 0.0;
  complex<double>* Gprobno = new complex<double>[N];


  int Neps = 200;
  double* eps = new double[Neps];
  double* dos = new double[Neps];
  for (int i=0; i<Neps; i++)
  {  eps[i] = -1.0 + i * 2.0/(Neps-1.0);
     dos[i] = DOS(DOStypes::SemiCircle,0.5,eps[i]);
  }
  //PrintFunc("dos",Neps,dos,eps);
  for(int n=N-1; n>=0; n--)
  {
    complex<double>* g = new complex<double>[Neps];
    for (int i=0; i<Neps; i++)
      g[i] =  dos[i] 
              / (ii * omega[n] + eta - Sigma_f[n] - eps[i]) ;
    
    Gprobno[n] = TrapezIntegral(Neps, g, eps);       
    
    delete [] g;
    //Gprobno[n] =                            1.0
    //               / ( ii * omega[n] + eta - Sigma_f[n] - sqr(t) * G_f[n] ) ;
    sum += real (    Gprobno[n] 
                );
  }
  
  //char fn[100];
  //sprintf(fn,"Gprobno.U%.3f.T%.3f", U, Temp );
  //PrintFunc(fn,N,Gprobno,omega);
  delete [] Gprobno;      
  delete [] eps;
  delete [] dos;
  printf("n: %f dn/dmu: %f\n",  2.0 * Temp * sum + 0.5, (2.0 * Temp * sum)/eta );
  return (2.0 * Temp * sum)/eta  ;

}


int intabs(int x) { if (x<0) return -x; else return x; };
//complex<double> sqr(complex<double> x) {return x*x; };
/*
complex<double> IAHFCHM::OpticalConductivity(int n)
{
  //-- bosonic frequency
  double nu = 2.0*pi*n*Temp; 
 

  double* iw_large = new double[2*N];
  complex<double>* Sig_large = new complex<double>[2*N];
  for(int i = 0; i < N; i++)
  {
    iw_large[N+i]    =  omega[i];
    iw_large[N-1-i]  = -omega[i];        
    Sig_large[N+i]   =  Sigma_f[i];
    Sig_large[N-1-i] = -Sigma_f[i];
  }

  int Neps = 2000;
  double* eps = new double[Neps];
  for(int i = 0; i < Neps; i++) 
    eps[i] = -2.0*t + i * 4.0 * t / ((double)(Neps-1));
  complex<double>* g = new complex<double>[Neps];
  complex<double>* summand = new complex<double>[2*N];
  complex<double> sum = 0.0;
  for (int m=0; m<2*N-n; m++)
  { for(int i =0; i<Neps; i++)
    {  g[i] = DOS(DOStypes::SemiCircle, t, eps[i]) 
             * (2.0/(2.0*pi*(double)n)) *
             (
               ( 1.0
                 / ( ( ii*iw_large[m] - eps[i] - Sig_large[m] ) 
                      * 
                     ( ii*iw_large[m+n] - eps[i] - Sig_large[m+n] ) 
                   )
               ) 
               -
               ( 1.0 / sqr( ii*iw_large[m] - eps[i] - Sig_large[m] ) 
               )
             );
      
    }
//    char integFN[300];
//    sprintf(integFN,"integ.n%d.m%d",n,m);
//    if ((m<10)and(n<10)) PrintFunc(integFN,Neps,g,eps);  

    summand[m] = TrapezIntegral(Neps, g, eps);
    sum +=  summand[m];
  }
  return sum;

  char summandFN[200];
  sprintf(summandFN,"summand.n%.4d",n);
  PrintFunc(summandFN,N,summand,omega); 

  delete [] iw_large;
  delete [] Sig_large;
  delete [] summand;
  delete [] eps;
  delete [] g;
}
*/
complex<double> IAHFCHM::OpticalConductivity(int n, double (*v)(double))
{
  return (1.0/2.0*pi*n)* (Lambda(n, v) - Lambda(0, v));
}

complex<double> IAHFCHM::Lambda(int n, double (*v)(double))
{
  //-- bosonic frequency
  double nu = 2.0*pi*n*Temp; 
 

  double* iw_large = new double[2*N];
  complex<double>* Sig_large = new complex<double>[2*N];
  for(int i = 0; i < N; i++)
  {
    iw_large[N+i]    =  omega[i];
    iw_large[N-1-i]  = -omega[i];        
    Sig_large[N+i]   =  Sigma_f[i];
    Sig_large[N-1-i] = -Sigma_f[i];
  }

  int Neps = 2000;
  double* eps = new double[Neps];
  for(int i = 0; i < Neps; i++) 
    eps[i] = -2.0*t + i * 4.0 * t / ((double)(Neps-1));
  complex<double>* g = new complex<double>[Neps];
  complex<double>* summand = new complex<double>[2*N];
  complex<double> sum = 0.0;
  for (int m=0; m<2*N-n; m++)
  { for(int i =0; i<Neps; i++)
    {  g[i] = DOS(DOStypes::SemiCircle, t, eps[i]) 
              * ( (v!=NULL) ? v(eps[i]) : 1.0 ) 
              * 2.0 
              * ( 1.0
                  / ( ( ii*iw_large[m] - eps[i] - Sig_large[m] ) 
                       * 
                      ( ii*iw_large[m+n] - eps[i] - Sig_large[m+n] ) 
                    )
                );
    }
//    char integFN[300];
//    sprintf(integFN,"integ.n%d.m%d",n,m);
//    if ((m<10)and(n<10)) PrintFunc(integFN,Neps,g,eps);  

    summand[m] = TrapezIntegral(Neps, g, eps);
    sum += /*(1.0/nu) * Temp */ summand[m];
  }
  return sum;

  char summandFN[200];
  sprintf(summandFN,"summand.n%.4d",n);
  PrintFunc(summandFN,N,summand,omega); 

  delete [] iw_large;
  delete [] Sig_large;
  delete [] summand;
  delete [] eps;
  delete [] g;
}

int find_index(int N, double* X, double x)
{
  for(int i = 0; i < N; i++)
  {
    if (X[i]>x) return i;   
  }  
  return -1;
}




complex<double> IAHFCHM::OpticalConductivity(double nu)
{

  double* iw_large = new double[2*N];
  complex<double>* Sig_large = new complex<double>[2*N];

  for(int i = 0; i < N; i++)
  {
    iw_large[N+i]    =  omega[i];
    iw_large[N-1-i]  = -omega[i];        
    Sig_large[N+i]   =  Sigma_f[i];
    Sig_large[N-1-i] = -Sigma_f[i];
  }
//  PrintFunc("Sig_large",2*N, Sig_large, iw_large);

/*
  FILE* f = fopen("Sig_large.interpl","w");
  for(double w = -100.0; w < 100.0; w+=0.3)
  {
    fprintf(f,"%.15le %.15le %.15le\n",w,real(interpl(2*N,Sig_large,iw_large,w)), imag(interpl(2*N,Sig_large,iw_large,w)));
  } 
  fclose(f);
*/

  int Neps = 2000;
  double* eps = new double[Neps];
  for(int i = 0; i < Neps; i++) 
    eps[i] = -2.0*t + i * 4.0 * t / ((double)(Neps-1));
  complex<double>* g = new complex<double>[Neps];
  complex<double>* summand = new complex<double>[2*N];
  complex<double> sum = 0.0;
  for (int m=0; m<2*N; m++)
  { 
    summand[m] = 0.0;
    double iw_plus_nu=iw_large[m]+nu;
    int nun = find_index(2*N, iw_large, iw_large[m]+nu);
    //if ((nun==m+1)and(m%1000==0)) printf("nun=m+1!!! m=%d, nu=%f\n",m,nu);
    if (nun==0) printf("nun=0!!! m=%d, nu=%f\n",m,nu);
    if (nun<0) { printf("nun<0!!! m=%d, nu=%f\n",m,nu); continue;}
    if ((nun>=2*N-2)or(nun==1)) { printf("nun+1 out of range. m=%d, nu=%f\n",m,nu); continue;}
    #pragma omp parallel for
    for(int i =0; i<Neps; i++)
    {
       complex<double> G_m = 1.0 / ( ii*iw_large[m] - eps[i] - Sig_large[m] );                 


       complex<double> G_nun_minus2 = 1.0 / ( ii*iw_large[nun-2] - eps[i] - Sig_large[nun-2] );
       complex<double> G_nun_minus1 = 1.0 / ( ii*iw_large[nun-1] - eps[i] - Sig_large[nun-1] );
       complex<double> G_nun = 1.0 / ( ii*iw_large[nun] - eps[i] - Sig_large[nun] );
       complex<double> G_nun_plus1 = 1.0 / ( ii*iw_large[nun+1] - eps[i] - Sig_large[nun+1] );
       complex<double> Gs[4] = {G_nun_minus2-G_m,G_nun_minus1-G_m,G_nun-G_m,G_nun_plus1-G_m};
       double iws[4] = {iw_large[nun-2]-iw_plus_nu, iw_large[nun-1]-iw_plus_nu, iw_large[nun]-iw_plus_nu, iw_large[nun+1]-iw_plus_nu};
       
       //complex<double> G_iw_plus_nu_minus_Gm = G_nun_minus1 + ((iw_large[m]+nu-iw_large[nun-1])/(iw_large[nun]-iw_large[nun-1]))*(G_nun-G_nun_minus1) - G_m;
       complex<double> G_iw_plus_nu_minus_Gm = CubicFrom4points(Gs, iws);
       //complex<double> G_iw_plus_nu_minus_Gm = ParabolaFrom3points(Gs+1, iws+1);

       //if ((m%4000==0)and(i%500==0))
       //  printf("iw: %f espilon: %f nu: %f nun: %d G_nun: (%f,%.15le) G_nun_minus1: (%f,%.15le) G_iw_plus_nu: (%f,%.15le)\n",
       //     iw_large[m],eps[i], nu, nun, real(G_nun), imag(G_nun), real(G_nun_minus1), imag(G_nun_minus1), real(G_iw_plus_nu), imag(G_iw_plus_nu) );
       //  printf("nu: %.2le iw: %.2le espilon: %.2le nun: %d Gs: %.5le,%.5le,%.5le,%.5le iws: %.5le,%.5le,%.5le,%.5le iw+nu: %.15le G_iw_plus_nu: %.15le\n",
       //       nu, iw_large[m], eps[i], nun,  imag(Gs[0]),imag(Gs[1]),imag(Gs[2]),imag(Gs[3]), iws[0],iws[1],iws[2],iws[3], iw_large[m]+nu, imag(G_iw_plus_nu) );
       g[i] = DOS(DOStypes::SemiCircle, t, eps[i]) 
             * (2.0*Temp/nu) * 
             (  G_m * G_iw_plus_nu_minus_Gm );      
    }
//    char integFN[300];
//    sprintf(integFN,"integ.n%d.m%d",n,m);
//    if ((m<10)and(n<10)) PrintFunc(integFN,Neps,g,eps);  

    summand[m] = TrapezIntegralMP(Neps, g, eps);
  }
  for (int m=0; m<2*N; m++)  sum += summand[m];

  char summandFN[200];
  sprintf(summandFN,"summand.nu%.3f",nu);
  PrintFunc(summandFN,2*N,summand,iw_large); 

  delete [] iw_large;
  delete [] Sig_large;
  delete [] summand;
  delete [] eps;
  delete [] g;

  return sum;

}

complex<double> IAHFCHM::OpticalConductivity2(double nu)
{

  double* iw_large = new double[2*N];
  complex<double>* Sig_large = new complex<double>[2*N];

  for(int i = 0; i < N; i++)
  {
    iw_large[N+i]    =  omega[i];
    iw_large[N-1-i]  = -omega[i];        
    Sig_large[N+i]   =  Sigma_f[i];
    Sig_large[N-1-i] = -Sigma_f[i];
  }
//  PrintFunc("Sig_large",2*N, Sig_large, iw_large);

/*
  FILE* f = fopen("Sig_large.interpl","w");
  for(double w = -100.0; w < 100.0; w+=0.3)
  {
    fprintf(f,"%.15le %.15le %.15le\n",w,real(interpl(2*N,Sig_large,iw_large,w)), imag(interpl(2*N,Sig_large,iw_large,w)));
  } 
  fclose(f);
*/

  int Neps = 2000;
  double* eps = new double[Neps];
  for(int i = 0; i < Neps; i++) 
    eps[i] = -2.0*t + i * 4.0 * t / ((double)(Neps-1));
  complex<double>* g = new complex<double>[Neps];
  complex<double>* summand = new complex<double>[2*N];
  complex<double> sum = 0.0;
  for (int m=0; m<2*N; m++)
  { 
    summand[m] = 0.0;
    double iw_plus_nu=iw_large[m]+nu;
    int nun = find_index(2*N, iw_large, iw_large[m]+nu);
    //if ((nun==m+1)and(m%1000==0)) printf("nun=m+1!!! m=%d, nu=%f\n",m,nu);
    if (nun==0){ printf("nun=0!!! m=%d, nu=%f\n",m,nu); continue;}
    if (nun<0) { printf("nun<0!!! m=%d, nu=%f\n",m,nu); continue;}
    if ((nun>=2*N-2)or(nun==1)) { printf("nun+1 out of range. m=%d, nu=%f\n",m,nu); continue;}

    complex<double> Sigs[4] = {Sig_large[nun-2],Sig_large[nun-1],Sig_large[nun],Sig_large[nun+1]};
    double iws[4] = {iw_large[nun-2]-iw_plus_nu, iw_large[nun-1]-iw_plus_nu, iw_large[nun]-iw_plus_nu, iw_large[nun+1]-iw_plus_nu};
       
    complex<double> Sig_iw_plus_nu = CubicFrom4points(Sigs, iws);

    #pragma omp parallel for
    for(int i =0; i<Neps; i++)
    {
       complex<double> G_m = 1.0 / ( ii*iw_large[m] - eps[i] - Sig_large[m] );                 
       complex<double> G_iw_plus_nu = 1.0 / ( ii*iw_plus_nu - eps[i] - Sig_iw_plus_nu );                 

       //if ((m%4000==0)and(i%500==0))
       //  printf("iw: %f espilon: %f nu: %f nun: %d G_nun: (%f,%.15le) G_nun_minus1: (%f,%.15le) G_iw_plus_nu: (%f,%.15le)\n",
       //     iw_large[m],eps[i], nu, nun, real(G_nun), imag(G_nun), real(G_nun_minus1), imag(G_nun_minus1), real(G_iw_plus_nu), imag(G_iw_plus_nu) );
       //  printf("nu: %.2le iw: %.2le espilon: %.2le nun: %d Gs: %.5le,%.5le,%.5le,%.5le iws: %.5le,%.5le,%.5le,%.5le iw+nu: %.15le G_iw_plus_nu: %.15le\n",
       //       nu, iw_large[m], eps[i], nun,  imag(Gs[0]),imag(Gs[1]),imag(Gs[2]),imag(Gs[3]), iws[0],iws[1],iws[2],iws[3], iw_large[m]+nu, imag(G_iw_plus_nu) );
       g[i] = DOS(DOStypes::SemiCircle, t, eps[i]) 
             * (2.0*Temp/nu) * 
             (  G_m * ( G_iw_plus_nu - G_m ) );      
    }
//    char integFN[300];
//    sprintf(integFN,"integ.n%d.m%d",n,m);
//    if ((m<10)and(n<10)) PrintFunc(integFN,Neps,g,eps);  

    summand[m] = TrapezIntegralMP(Neps, g, eps);
  }
  for (int m=0; m<2*N; m++)  sum += summand[m];

  char summandFN[200];
  sprintf(summandFN,"summand.nu%.3f",nu);
  PrintFunc(summandFN,N,summand,omega); 

  delete [] iw_large;
  delete [] Sig_large;
  delete [] summand;
  delete [] eps;
  delete [] g;

  return sum;

}

/*


complex<double> IAHFCHM::OpticalConductivity(double nu)
{

  double* iw_large = new double[2*N];
  complex<double>* Sig_large = new complex<double>[2*N];

  for(int i = 0; i < N; i++)
  {
    iw_large[N+i]    =  omega[i];
    iw_large[N-1-i]  = -omega[i];        
    Sig_large[N+i]   =  Sigma_f[i];
    Sig_large[N-1-i] = -Sigma_f[i];
  }
//  PrintFunc("Sig_large",2*N, Sig_large, iw_large);


//  FILE* f = fopen("Sig_large.interpl","w");
//  for(double w = -100.0; w < 100.0; w+=0.3)
//  {
//    fprintf(f,"%.15le %.15le %.15le\n",w,real(interpl(2*N,Sig_large,iw_large,w)), imag(interpl(2*N,Sig_large,iw_large,w)));
//  } 
//  fclose(f);


  int Neps = 2000;
  double* eps = new double[Neps];
  for(int i = 0; i < Neps; i++) 
    eps[i] = -2.0*t + i * 4.0 * t / ((double)(Neps-1));
  complex<double>* g = new complex<double>[Neps];
  complex<double>* summand = new complex<double>[2*N*10];
  complex<double> sum = 0.0;
  for (int m=0; m<2*N; m++)
  { if (iw_large[m]+nu>iw_large[2*N-2]) break;
    for(double w=iw_large[m]; w<iw_large[m+1]-((iw_large[m+1]-iw_large[m])/20.0); w+=(iw_large[m+1]-iw_large[m])/10.0)
    {
      complex<double> Sig_iw_plus_nu = interpl(2*N,Sig_large,iw_large,w+nu);
      complex<double> Sig_iw = interpl(2*N,Sig_large,iw_large,w);

      #pragma omp parallel for
      for(int i =0; i<Neps; i++)
      {  g[i] = DOS(DOStypes::SemiCircle, t, eps[i]) 
               * (2.0*Temp/nu) *
               (
                 ( 1.0
                   / ( ( ii*w - eps[i] - Sig_iw ) 
                        * 
                       ( ii*(w+nu) - eps[i] - Sig_iw_plus_nu ) 
                     )
                 ) 
                 -
                 ( 1.0 / sqr( ii*w - eps[i] - Sig_iw ) 
                 )
               );
       
      }
//    char integFN[300];
//    sprintf(integFN,"integ.n%d.m%d",n,m);
//    if ((m<10)and(n<10)) PrintFunc(integFN,Neps,g,eps);  

      summand[m] = TrapezIntegralMP(Neps, g, eps);
      sum +=  summand[m]; //   (1.0/nu) * Temp 
    }
  }
  return sum;

  char summandFN[200];
  sprintf(summandFN,"summand.nu%.3f",nu);
  PrintFunc(summandFN,N,summand,omega); 

  delete [] iw_large;
  delete [] Sig_large;
  delete [] summand;
  delete [] eps;
  delete [] g;
}

*/



//-----------------------default Constructor----------------------//

IAHFCHM::IAHFCHM()
{
  Initialized = false;
 
  //--default options--
  print_initial = true;
  print_intermediate = 0;
  print_lambda = false; 
  print_final = true;
  show_status = true;
}

IAHFCHM::~IAHFCHM()
{
  iagrid->~IAGRID();
  //delete [] iagrid;
  delete [] G_f;
  delete [] G_0_f;
  delete [] G_t;
  delete [] Sigma_t;
  delete [] Sigma_f;
}


