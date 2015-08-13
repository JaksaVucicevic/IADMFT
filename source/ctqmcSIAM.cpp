#include <complex>
#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "ctqmcSIAM.h"
#include "IAResult.h"
#include "IAGRID.h"
#include "routines.h"

using namespace std;

ctqmcSIAM::ctqmcSIAM()
{ printf("ctqmcSIAM::ctqmcSIAM()\n");
  Defaults();
}

//------------------- copy constructor -----------------------------//
/*ctqmcSIAM::ctqmcSIAM(const ctqmcSIAM &s)
{
  printf("COPY constructor ctqmcSIAM::ctqmcSIAM(s)\n");

  this->DoPrintParamsFile = s.DoPrintParamsFile;

  this->execName = s.execName; //full path!
  this->runFolderName = s.runFolderName; //full path!
  this->runParallel = s.runParallel;

  this->inputf   = s.inputf;      // filename of the parameter file
  this->fDelta   = s.fDelta; // input bath function Delta(iom)
  this->fcix     = s.fcix;   // input file with atom eigenvalues and eigenvectors (exact diagonalization results)
  this->fGf      = s.fGf;    // output file for green's function G(iom)
  this->fSig     = s.fSig;   // output file for self-energy Sig(iom)
  this->mu       = s.mu;           // chemical potential
  this->beta     = s.beta;         // inverse temperature
  this->U        = s.U;           // Coulomb repulsion (should be used only in single site DMFT)
  this->M	 = s.M;     // number of all qmc steps
  this->Ntau     = s.Ntau;         // number of imaginary time points for Delta(tau)
  this->Nmax     = s.Nmax;        // maximum number of kinks
  this->nom      = s.nom;         // number of frequency points sample in qmc
  this->nomb     = s.nomb;         //number of bosonic frequencies for susceptibilities
  this->PChangeOrder 	= s.PChangeOrder;     // probability to change the order (as compared to move a kink)
  this->tsample  	= s.tsample;          // how often to record measurements (each tsample steps)
  this->warmup   	= s.warmup;      // how many qmc steps should be ignored
  this->CleanUpdate 	= s.CleanUpdate;      // clean update is sometimes useful
  this->minM     = s.minM;       // trace shuld always be larger than this minimal value when accepting the step
  this->minD     = s.minD;       // determinant of hybrodization should always be larger than this number when accepting the step
  this->Ncout    = s.Ncout;       // screen output for basic information every Ncout steps
  this->Naver    = s.Naver;     // much more information about simulation every Naver steps
  this->TwoKinks = s.TwoKinks;          // probability to add two kinks
  this->GlobalFlip	= s.GlobalFlip;         // global flip is tried after GlobalFlip qmc steps
  this->treshold 	= s.treshold;       // energy to throw away atomic states when creating new.cix for next iteration
  this->SampleGtau 	= s.SampleGtau;          // How often to sample for Gtau (if at all)
  this->SampleVertex 	= s.SampleVertex;          // How often to sample two particle vertex
  this->ReduceComm 	= s.ReduceComm;           // only part of Gtau is transfered between nodes (if slow communication)
  this->Ncorrect 	= s.Ncorrect;          // baths with index higher than Ncorrect should not be added a high-frequency tail (if -1, all tails are added)
  this->aom 		= s.aom;           // number of frequency points to find Sigma(om_last_sampled)
  this->som 		= s.som;           // number of frequency points to find Susc(om_last_sampled)
  this->PreciseP 	= s.PreciseP;           // computes probabilities more precisely
  this->sderiv 		= s.sderiv;         // maximum mismatch when matching high and low frequency of imaginary part of Sigmas
  this->minDeltat	= s.minDeltat;        // Delta(tau) is sometimes not causal due to numerical error. In this case we set it to small value minDeltat.
  this->SampleSusc 	= s.SampleSusc;       // If spin and charge dynamic susceptibility should be sampled during simulation
  this->nomv 		= s.nomv;       // number of fermionic frequencies for computing vertex
  this->nOm 		= s.nOm;     
}
*/
ctqmcSIAM::~ctqmcSIAM()
{
 // whateva
}


void ctqmcSIAM::PrintFunc(const char* FileName, int N, complex<double>* Y, double* X)
{
  FILE *f;
  f = fopen(FileName, "w");
  fprintf(f,"#whateva\n");
  for (int i=0; i<N; i++)
    fprintf(f,"%.15le %.15le %.15le\n", X[i], real(Y[i]), imag(Y[i]));
  fclose(f);
}

void ctqmcSIAM::ReadXY(const char* FileName, int N, complex<double>* Y, double* X)
{ 
  FILE *f;
  f = fopen(FileName, "r");
  char str[1000];
  fgets(str,1000,f);
  for (int i=0; i<N; i++)
  { 
    double x, Rey, Imy;
    //fscanf(f, "%e %e %e", &x, &Rey, &Imy);
    fscanf(f, "%le", &x);
    fscanf(f, "%le", &Rey);
    fscanf(f, "%le", &Imy);

    Y[i]=complex<double>(Rey, Imy);
    X[i] = x;
  }
  fclose(f);
}

void ctqmcSIAM::Defaults()
{
  printf("ctqmcSIAM::Defaults()\n");
  
  DoPrintParamsFile = true;
  UseSmartNom = true;
  FreqCutoff = 10.0;

  PatchTailWithAtomicLimit=false;
  AtomicCutoff=12.0;

  execName = "ctqmc"; //full path!
  runFolderName = "run_ctqmc"; //full path!
  runParallel = false;

  inputf = "PARAMS";      // filename of the parameter file
  fDelta   = "Delta"; // input bath function Delta(iom)
  fcix     = "icix.cix";   // input file with atom eigenvalues and eigenvectors (exact diagonalization results)
  fGf      = "Gf";    // output file for green's function G(iom)
  fSig     = "Sigma";   // output file for self-energy Sig(iom)
  mu       = 0;           // chemical potential
  beta     = 100;         // inverse temperature
  U        = 0;           // Coulomb repulsion (should be used only in single site DMFT)
  M = 50000000;     // number of all qmc steps
  Ntau     = 500;         // number of imaginary time points for Delta(tau)
  Nmax     = 1024;        // maximum number of kinks
  nom      = 50;         // number of frequency points sample in qmc
  nomb     = 50;          // number of bosonic frequencies for susceptibilities
  PChangeOrder = 0.9;     // probability to change the order (as compared to move a kink)
  tsample  = 10;          // how often to record measurements (each tsample steps)
  warmup   = 500000;      // how many qmc steps should be ignored
  CleanUpdate = 100000;      // clean update is sometimes useful
  minM     = 1e-10;       // trace shuld always be larger than this minimal value when accepting the step
  minD     = 1e-10;       // determinant of hybrodization should always be larger than this number when accepting the step
  Ncout    = 1000000000;       // screen output for basic information every Ncout steps
  Naver    = 1000000000;     // much more information about simulation every Naver steps
  TwoKinks = 0.;          // probability to add two kinks
  GlobalFlip= -1;         // global flip is tried after GlobalFlip qmc steps
  treshold = 1e-10;       // energy to throw away atomic states when creating new.cix for next iteration
  SampleGtau = 100000;          // How often to sample for Gtau (if at all)
  SampleVertex =-1;          // How often to sample two particle vertex
  ReduceComm = 0;           // only part of Gtau is transfered between nodes (if slow communication)
  Ncorrect = -1;          // baths with index higher than Ncorrect should not be added a high-frequency tail (if -1, all tails are added)
  aom = 3;           // number of frequency points to find Sigma(om_last_sampled)
  som = 3;           // number of frequency points to find Susc(om_last_sampled)
  PreciseP = 1;           // computes probabilities more precisely
  sderiv = 0.1;         // maximum mismatch when matching high and low frequency of imaginary part of Sigmas
  minDeltat= 1e-7;        // Delta(tau) is sometimes not causal due to numerical error. In this case we set it to small value minDeltat.
  SampleSusc = true;       // If spin and charge dynamic susceptibility should be sampled during simulation
  nomv = 10000;       // number of fermionic frequencies for computing vertex
  nOm = 1;           // number of bosonic frequencies for vertex calculation
}

void ctqmcSIAM::PrintParamsFile()
{
  char FN[300];
  sprintf(FN, "%s/%s", runFolderName.c_str(), inputf.c_str() );
  printf("PARAMS FILE: %s, inputf: %s\n",FN, inputf.c_str());
  FILE* f = fopen(FN,"w");
  
  fprintf(f, "Delta %s\n",fDelta.c_str());			// input bath function Delta(iom)
  fprintf(f, "cix %s\n",fcix.c_str());				// input file with atom eigenvalues and eigenvectors (exact diagonalization results)
  fprintf(f, "Gf %s\n",fGf.c_str());				// output file for green's function G(iom)
  fprintf(f, "Sig %s\n",fSig.c_str());				// output file for self-energy Sig(iom)
  fprintf(f, "mu %f\n",mu);				// chemical potential
  fprintf(f, "beta %f\n",beta);				// inverse temperature
  fprintf(f, "U %f\n",U);				// Coulomb repulsion (should be used only in single site DMFT)
  fprintf(f, "M %u\n",M);				// number of all qmc steps
  fprintf(f, "Ntau %d\n",Ntau);				// number of imaginary time points for Delta(tau)
  fprintf(f, "Nmax %d\n",Nmax);				// maximum number of kinks
  fprintf(f, "nom %d\n",nom);				// number of frequency points sample in qmc
  fprintf(f, "nomb %d\n",nomb);				// number of bosonic frequencies for susceptibilities
  fprintf(f, "PChangeOrder %f\n",PChangeOrder);		// probability to change the order (as compared to move a kink)
  fprintf(f, "tsample %d\n",tsample);			// how often to record measurements (each tsample steps)
  fprintf(f, "warmup %d\n",warmup);			// how many qmc steps should be ignored
  fprintf(f, "CleanUpdate %d\n",CleanUpdate);		// clean update is sometimes useful
  fprintf(f, "minM %le\n",minM);			// trace shuld always be larger than this minimal value when accepting the step
  fprintf(f, "minD %le\n",minD);			// determinant of hybrodization should always be larger than this number when accepting the step
  fprintf(f, "Ncout %d\n",Ncout);			// screen output for basic information every Ncout steps
  fprintf(f, "Naver %d\n",Naver);			// much more information about simulation every Naver steps
  fprintf(f, "TwoKinks %f\n",TwoKinks);			// probability to add two kinks
  fprintf(f, "GlobalFlip %d\n",GlobalFlip);		// global flip is tried after GlobalFlip qmc steps
  fprintf(f, "treshold %le\n",treshold);		// energy to throw away atomic states when creating new.cix for next iteration
  fprintf(f, "SampleGtau %d\n",SampleGtau);		// How often to sample for Gtau (if at all)
  fprintf(f, "SampleVertex %d\n",SampleVertex);		// How often to sample two particle vertex
  fprintf(f, "ReduceComm %d\n",ReduceComm);		// only part of Gtau is transfered between nodes (if slow communication)
  fprintf(f, "Ncorrect %d\n",Ncorrect);			// baths with index higher than Ncorrect should not be added a high-frequency tail (if -1, all tails are added)
  fprintf(f, "aom %d\n",aom);				// number of frequency points to find Sigma(om_last_sampled)
  fprintf(f, "som %d\n",som);				// number of frequency points to find Susc(om_last_sampled)
  fprintf(f, "PreciseP %d\n",PreciseP);			// computes probabilities more precisely
  fprintf(f, "sderiv %f\n",sderiv);			// maximum mismatch when matching high and low frequency of imaginary part of Sigmas
  fprintf(f, "minDeltat %le\n",minDeltat);		// Delta(tau) is sometimes not causal due to numerical error. In this case we set it to small value minDeltat.
  fprintf(f, "SampleSusc %d\n",(SampleSusc) ? 1 : 0);	// If spin and charge dynamic susceptibility should be sampled during simulation
  fprintf(f, "nomv %d\n",nomv);				// number of fermionic frequencies for computing vertex
  fprintf(f, "nOm %d\n",nOm);				// number of bosonic frequencies for vertex calculation

  fclose(f); 
}

void ctqmcSIAM::PrintCixFile()
{
  char FN[300];
  sprintf(FN, "%s/%s", runFolderName.c_str(), fcix.c_str() );
  FILE * CixFile = fopen(FN,"w");

  fprintf(CixFile, "# Cix file for cluster DMFT with CTQMC\n");
  fprintf(CixFile, "# cluster_size, number of states, number of baths\n");
  fprintf(CixFile, "1 4 2 1\n");
  fprintf(CixFile, "# baths, dimension, symmetry\n");
  fprintf(CixFile, "0       1 0 0\n");
  fprintf(CixFile, "1       1 0 0\n");
  fprintf(CixFile, "# cluster energies for unique baths, eps[k]\n");
  fprintf(CixFile, "0 0\n");
  fprintf(CixFile, "#   N   K   Sz size\n");
  fprintf(CixFile, "1   0   0    0   1     2  3     0   0\n");
  fprintf(CixFile, "2   1   0 -0.5   1     0  4     0   0.5\n");
  fprintf(CixFile, "3   1   0  0.5   1     4  0     0   0.5\n");
  fprintf(CixFile, "4   2   0    0   1     0  0     0   0\n");
  fprintf(CixFile, "# matrix elements\n");
  fprintf(CixFile, "1  2  1  1    1\n");
  fprintf(CixFile, "1  3  1  1    1\n");
  fprintf(CixFile, "2  0  0  0\n");
  fprintf(CixFile, "2  4  1  1   -1\n");
  fprintf(CixFile, "3  4  1  1    1\n");
  fprintf(CixFile, "3  0  0  0\n");
  fprintf(CixFile, "4  0  0  0\n");
  fprintf(CixFile, "4  0  0  0\n");
  fprintf(CixFile, "HB1\n");
  fprintf(CixFile, "# number of operators needed\n");
  fprintf(CixFile, "0\n");
  fprintf(CixFile, "# Data for HB1\n");
  fprintf(CixFile, "1 4 2 1\n");
  fprintf(CixFile, "#  ind  N   K   Sz size\n");
  fprintf(CixFile, "1  1   0   0    0   1     2  3     0   0\n");
  fprintf(CixFile, "2  2   1   0 -0.5   1     0  4     0   0.5\n");
  fprintf(CixFile, "3  3   1   0  0.5   1     4  0     0   0.5\n");
  fprintf(CixFile, "4  4   2   0    0   1     0  0     0   0\n");
  fprintf(CixFile, "# matrix elements\n");
  fprintf(CixFile, "1  2  1  1    1\n");
  fprintf(CixFile, "1  3  1  1    1\n");
  fprintf(CixFile, "2  0  0  0\n");
  fprintf(CixFile, "2  4  1  1   -1\n");
  fprintf(CixFile, "3  4  1  1    1\n");
  fprintf(CixFile, "3  0  0  0\n");
  fprintf(CixFile, "4  0  0  0\n");
  fprintf(CixFile, "4  0  0  0\n");
  
  fclose(CixFile);
}



double ctqmcSIAM::GetOccupancy()
{
  char occFN[300];
  sprintf(occFN, "%s/observables",runFolderName.c_str());          
  FILE* occFile = fopen(occFN,"r");
  double occ;
  fscanf(occFile,"%le",&occ);
  fclose(occFile);

  return occ/2.0;
}


bool ctqmcSIAM::Run(IAResult* r)
{
  printf("ctqmcSIAM::Run\n");

  if (UseSmartNom)
  { nom = (int) ( 0.5 * ( FreqCutoff/(pi/beta) - 1.0 ) );
    if (nom == 0) nom = 1;
    if (nom > 100) nom = 100;
  }
 
  char cmd[300];
  sprintf(cmd, "mkdir %s",runFolderName.c_str());
  system(cmd);

  //sprintf(cmd, "cp %s %s/", execName.c_str(), runFolderName.c_str());
  //system(cmd);

  if (DoPrintParamsFile) PrintParamsFile();
  PrintCixFile();

  char FN[300];
  
  sprintf(FN,"%s/%s",runFolderName.c_str(), fDelta.c_str());
  PrintFunc(FN, r->iagrid->get_N(), r->Delta, r->omega);

 // sprintf(cmd, "cd %s",runFolderName.c_str());
 // system(cmd);

  if (runParallel)
    sprintf(cmd, "cd %s ; rm observables ; $MPI_OPENMPI_MPIEXEC %s",runFolderName.c_str(), execName.c_str()); 
  else 
    sprintf(cmd, "cd %s ; rm observables ; %s",runFolderName.c_str(), execName.c_str());
  system(cmd);

  sprintf(FN,"%s/%s",runFolderName.c_str(), fGf.c_str());
  ReadXY(FN, r->iagrid->get_N(), r->G, r->omega);

  sprintf(FN,"%s/%s",runFolderName.c_str(), fSig.c_str());
  ReadXY(FN, r->iagrid->get_N(), r->Sigma, r->omega);
 
  if (PatchTailWithAtomicLimit)
  {  double n = r->n;
     r->n = GetOccupancy(); //this is now in ctqmcCHM::SolveSiam
     r->PatchAtomicLimitG(AtomicCutoff, U);
     r->PatchAtomicLimitSigma(AtomicCutoff, U);
     r->n = n;
  } 
}
