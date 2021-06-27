// =========================================================================================
// Evaluate the parameters for the CHAMP iteration for partitioned schemes.
// =========================================================================================
#include "Cgad.h"
#include "Interface.h"  


// #include "CompositeGridOperators.h"
// #include "ParallelUtility.h"
// //#include "AdamsPCData.h"
// //#include "AdvanceOptions.h"
// //#include "OGPulseFunction.h"
// //#include "OGTrigFunction.h"
// #include "Oges.h"

// #include "SparseRep.h"

// #include "gridFunctionNorms.h"


// =========================================================================================
/// \brief Evaluate the parameters for the CHAMP iteration for partitioned schemes.
///
/// Output:
///      champParameters(0,side,axis,grid)=pl;     // optimized Scwartz Parameter for side 1 : Sl = pl/dx
///      champParameters(1,side,axis,grid)=pr;     // optimized Scwartz Parameter for side 2 : Sr = pr/dx 
///      champParameters(2,side,axis,grid)=theta;  // K1/K2
///      champParameters(3,side,axis,grid)=beta;   // D1/D2        
// =========================================================================================
int getChampParameters( int grid, int side, int axis, real dt, real dx, Parameters & parameters, RealArray & champParameters )
{

   // Retrieve the parameters from the opposite side of the interface:
  Parameters & paramRight = getInterfaceParameters( grid, side, axis, parameters);

  // This only works if the opposite side is a Cgad object *** FIX ME ***
  std::vector<real> & kappaLeft  = parameters.dbase.get<std::vector<real> >("kappa");
  std::vector<real> & kappaRight = paramRight.dbase.get<std::vector<real> >("kappa");
  const Real & kThermalLeft  = parameters.dbase.get<Real>("thermalConductivity");
  const Real & kThermalRight = paramRight.dbase.get<Real>("thermalConductivity");

  Real theta = kThermalLeft/kThermalRight;
  Real beta = kappaLeft[0]/kappaRight[0];

  Real lambdaD = kappaLeft[0]*dt/(dx*dx); 

  // do this for now 
  Real pl =1.;
  Real pr =1.; 
  // printF("***** fillChampBC: (side,axis,grid)=(%d,%d,%d) : theta=Kleft/Kright=%g, beta=Dleft/Dright=%g, Sl=%g, Sr=%g\n ****",side,axis,grid,theta,beta,Sl,Sr);

  champParameters(0,side,axis,grid)=pl;     // optimized Scwartz Parameter for side 1 : Sl = pl/dx
  champParameters(1,side,axis,grid)=pr;     // optimized Scwartz Parameter for side 2 : Sr = pr/dx 
  champParameters(2,side,axis,grid)=theta;  // K1/K2
  champParameters(3,side,axis,grid)=beta;   // D1/D2         


  return 0;
}