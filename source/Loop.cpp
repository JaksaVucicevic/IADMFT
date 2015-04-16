#include "Broyden.h"
#include "Mixer.h"
#include "IAResult.h"
#include "IAresArray.h"
#include "Input.h"
#include "IAGRID.h"
#include "Loop.h"
#include "LambdaCalculator.h"

void Loop::Defaults()
{
    printf("\n\n\n\nLOOP DEFAULTS <<<<<<<<<<<<<<<<<<<<\n\n");
    UseBroyden = true;
    ForceBroyden = false;
    BroydenStartDiff =5e-3;

    //---- Mixer Options ------//
    NtoMix = 2;
    Coefs = new int[NtoMix];
    Coefs[0] = 1;
    Coefs[1] = 0;
    
    //---- Loop Options -------//
    MAX_ITS = 300;
    MIN_ITS = 5;
    Accr = 5e-5;

    //---- PrintOut/Debugging optins----//
    PrintIntermediate = false;
    HaltOnIterations = false;
    ForceSymmetry = false;

    LC = new LambdaCalculator;
}

Loop::Loop()
{
  Defaults();
}

Loop::Loop(const char* ParamsFN)
{
  Defaults();
  
  this->ParamsFN.assign(ParamsFN);
  cout << "-- INFO -- Loop: Params File Name set to:" << this->ParamsFN << endl;

  Input input(ParamsFN);
  
  input.ReadParam(UseBroyden,"Loop::UseBroyden");
  input.ReadParam(ForceBroyden,"Loop::ForceBroyden");
  input.ReadParam(BroydenStartDiff,"Loop::BroydenStartDiff");
  input.ReadParam(NtoMix,"Loop::NtoMix");
  delete [] Coefs;
  Coefs = new int[NtoMix];
  Coefs[0] = 1;
  printf("-------- LOOP:: C0 = %d, C1 = %d\n",Coefs[0],Coefs[1]);   
  input.ReadArray(NtoMix, Coefs, "Loop::Coefs");
  input.ReadParam(MAX_ITS,"Loop::MAX_ITS");
  input.ReadParam(MIN_ITS,"Loop::MIN_ITS");
  input.ReadParam(Accr,"Loop::Accr");
  input.ReadParam(PrintIntermediate,"Loop::PrintIntermediate");
  input.ReadParam(HaltOnIterations,"Loop::HaltOnIterations");
  input.ReadParam(ForceSymmetry,"Loop::ForceSymmetry");
  
  printf("HOI: %s PI: %s\n", (HaltOnIterations) ? "yes" : "no", (PrintIntermediate) ? "yes" : "no");

  LC = new LambdaCalculator(ParamsFN);
}

void Loop::ReleaseMemory()
{
  printf("Loop release\n");
  delete [] Coefs;
  LC->~LambdaCalculator(); 
}

Loop::~Loop()
{
  ReleaseMemory();
}

//------- Option Stetters --------//

void Loop::SetGrid(IAGRID* g)
{
  iagrid = g;
  N = iagrid->get_N();
}

void Loop::SetMixerOptions(int NtoMix, const int * Coefs)
{
  delete [] this->Coefs;
  this->NtoMix = NtoMix;
  this->Coefs = new int[NtoMix];
  for (int i=0; i<NtoMix; i++) this->Coefs[i] = Coefs[i];
}

void Loop::SetBroydenOptions(bool UseBroyden, bool ForceBroyden, double BroydenStartDiff)
{
  this->UseBroyden = UseBroyden;
  this->ForceBroyden = ForceBroyden;
  this->BroydenStartDiff = BroydenStartDiff;
}

void Loop::SetLoopOptions(int MAX_ITS, double Accr)
{
  this->MAX_ITS = MAX_ITS;
  this->Accr = Accr;
}

void Loop::SetPrintOutOptions(bool PrintIntermediate, bool HaltOnIterations)
{
  this->PrintIntermediate = PrintIntermediate;
  this->HaltOnIterations = HaltOnIterations;
  if (HaltOnIterations) printf("-- INFO -- Loop: Halt on iterations set ON\n");
}

void Loop::SetForceSymmetry(bool FS)
{
  ForceSymmetry = FS;

}

//=======================================================================================//
// Run with single result
//=======================================================================================//

bool Loop::Run(IAResult* r)
{
  this->r = r;
  this->iagrid = r->iagrid;
  N = r->iagrid->get_N();
  
  //Initialize mixer
  printf("|||||||||||||||||||||||||| LOOP:: C0 = %d, C1 = %d\n",Coefs[0],Coefs[1]);
  Mixer< complex<double> > mixer(N, NtoMix, Coefs, (UseBroyden) ? BroydenStartDiff : Accr);
  mixer.Mix(r->Delta);

  //initialize broyden
  Broyden B;
  B.SetParameters(N, MAX_ITS, 1.0, 0.01, Accr);
  int BroydenStatus = 0;
  // Broyden status: 0 - Waiting for mixer to reach BroydenStartDiff
  //                 1 - Running
  //                 2 - Suspended

  //-------LambdaCalculatorPrepare------//  
  LC->ResetCounter();
  LC->SetOmega(r->omega); 
  LC->SetN(N);
  LC->SetOffset(0);


  //Halt on first iteration if HaltOnIterations
  int Halt = (HaltOnIterations) ? 1 : 0; 

  bool converged = false;
  //------------ DMFT loop-------------//
  for (int it = 1; it<=MAX_ITS; it++)
  {  printf("--- DMFT Loop Iteration %d ---\n", it);
     Iteration = it;
    
     //----- solve SIAM ------//
     if ( SolveSIAM() ) return true;
     //-----------------------//

     //halt
     if (it==Halt)
     {
       r->PrintResult("intermediate");
       printf("Next stop: ");
       cin >> Halt; 
     }

     //print out intermediate results
     if (PrintIntermediate)
     {  char FN[50];
        sprintf(FN,"intermediate.%d", it);
        //sprintf(FN,"intermediate");
        r->PrintResult(FN);
     }

     //--- self-consistency ---// 
     CalcDelta(); 
     //------------------------//

     // check for nans
     for (int i = 0; i < N; i++) if ( r->Delta[i] != r->Delta[i] ) { printf("nan in Delta!!!!\n"); return true; }

     // clip off
     bool ClippingDelta= false;

     // force symmetry
     if (ForceSymmetry)
       for (int i = 0; i < N; i++)
       { r->Delta[i] =complex<double>( 0, imag(r->Delta[i]) );
         r->Sigma[i] =complex<double>( 0, imag(r->Sigma[i]) );
         r->SOCSigma[i] =complex<double>( 0, imag(r->SOCSigma[i]) );
         r->G[i] =complex<double>( 0, imag(r->G[i]) );
       }

     printf("   Loop: mixing and checking convergence...\n");
     // now mix and check if converged
     int conv = 0;
     if (BroydenStatus == 1) 
       conv = B.CalculateNew(r->Delta,it);
     else
     { if (mixer.Mix(r->Delta))
         if ((UseBroyden)and(BroydenStatus == 0)) 
         { B.TurnOn(it); //switch to broyden if mixer converged
           BroydenStatus++;
           B.CurrentDiff = mixer.CurrentDiff;
         } 
         else conv = 1;
     }

     LC->CalculateLambda(r->Delta);
     
     if ((conv==1)and(it>MIN_ITS)) { converged = true; break; }
  }
  //-----------------------------------//
  
  if (BroydenStatus == 1) B.TurnOff();
  return !converged;
}

//=======================================================================================//
// Run with IAresArray - applicable for StatDMFT
//=======================================================================================//


bool Loop::Run(IAresArray* a)
{
  this->a = a;
  int Nsites = a->get_N();
  this->iagrid = a->r[0].iagrid;
  N = iagrid->get_N();
  
  //Initialize mixer
  printf("|||||||||||||||||||||||||| LOOP:: C0 = %d, C1 = %d\n",Coefs[0],Coefs[1]);
  Mixer< complex<double> > mixer(N*Nsites, NtoMix, Coefs, (UseBroyden) ? BroydenStartDiff : Accr);
  mixer.Mix(a->totalDelta);

  //initialize broyden
  Broyden B;
  B.SetParameters(N*Nsites, MAX_ITS, 1.0, 0.01, Accr);
  int BroydenStatus = 0;
  // Broyden status: 0 - Waiting for mixer to reach BroydenStartDiff
  //                 1 - Running
  //                 2 - Suspended

  //-------LambdaCalculatorPrepare------//  
  LC->ResetCounter();
  LC->SetOmega(a->r[0].omega); 
  LC->SetN(N*Nsites);
  LC->SetOffset(0);


  //Halt on first iteration if HaltOnIterations
  int Halt = (HaltOnIterations) ? 1 : 0; 

  bool converged = false;
  //------------ DMFT loop-------------//
  for (int it = 1; it<=MAX_ITS; it++)
  {  printf("--- DMFT Loop Iteration %d ---\n", it);
     Iteration = it;
    
     //----- solve SIAM ------//
     if ( SolveSIAM() ) return true;
     //-----------------------//

     //halt
     if (it==Halt)
     { char FN[300];
       sprintf(FN,"intermediate/");
       char command[300];
       sprintf(command,"mkdir %s",FN);  
       a->PrintAll(FN);
       printf("Next stop: ");
       cin >> Halt; 
     }

     //print out intermediate results
     if (PrintIntermediate)
     { char FN[300];
       sprintf(FN,"intermediate.%d/", it);
       char command[300];
       sprintf(command,"mkdir %s",FN);  
       a->PrintAll(FN);
     }

     //--- self-consistency ---// 
     CalcDelta(); 
     //------------------------//

     // check for nans
     for (int id = 0; id<Nsites; id++) 
     for (int i = 0; i < N; i++) 
       if ( a->r[id].Delta[i] != a->r[id].Delta[i] ) 
       { printf("nan in Delta!!!!\n"); return true; }
   
     if (ForceSymmetry)
       for (int id = 0; id<Nsites; id++) 
       for (int i = 0; i < N; i++) 
         a->r[id].Delta[i] =complex<double>( 0, imag(a->r[id].Delta[i]) );

     a->WriteTotalDelta();

     printf("   Loop: mixing and checking convergence...\n");
     // now mix and check if converged
     int conv = 0;
     if (BroydenStatus == 1) 
       conv = B.CalculateNew(a->totalDelta,it);
     else
     { if (mixer.Mix(a->totalDelta))
         if ((UseBroyden)and(BroydenStatus == 0)) 
         { B.TurnOn(it); //switch to broyden if mixer converged
           BroydenStatus++;
           B.CurrentDiff = mixer.CurrentDiff;
         } 
         else conv = 1;
     }

     LC->CalculateLambda(a->totalDelta);

     a->ReadTotalDelta();
     
     if ((conv==1)and(it>MIN_ITS)) { converged = true; break; }
  }
  //-----------------------------------//
  
  if (BroydenStatus == 1) B.TurnOff();
  return !converged;
}



bool Loop::SolveSIAM()
{
  printf("SS Loop"); 
  return false;
}

void Loop::CalcDelta()
{
  printf("CD Loop");
}
