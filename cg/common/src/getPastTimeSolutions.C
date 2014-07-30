#include "DomainSolver.h"
// #include "CompositeGridOperators.h"
// #include "GridCollectionOperators.h"
// #include "interpPoints.h"
// #include "ExposedPoints.h"
// #include "PlotStuff.h"
// #include "InterpolateRefinements.h"
// #include "Regrid.h"
// #include "Ogen.h"
// #include "MatrixTransform.h"
// #include "updateOpt.h"
// #include "App.h"
// #include "ParallelUtility.h"
// #include "Oges.h"
// #include "AdamsPCData.h"



// ===========================================================================================
/// \brief Determine a past time solution and grids (needed for predictor corrector schemes)
///
/// \param current (input) : index into the gf array of the current (t=0) solution
/// \param numberOfPast  (input) : number of past time solutions to compute
/// \param 
// ===========================================================================================
int DomainSolver::
getPastTimeSolutions( int current, int numberOfPast, int *previous  )
{

  printF("--DS-- DomainSolver::getPastTimeSolution: current=%i (t=%8.2e) numberOfPast=%i\n",
	 current,gf[current].t, numberOfPast );

  for( int past=0; past<numberOfPast; past ++ )
  {
    const int prev = previous[past];
    
    real tPast = gf[prev].t;

    printF("--DS-- DomainSolver::getPastTimeSolution: construct past time solution: prev=%i, tPast=%8.2e\n",prev,tPast);
    
    if( movingGridProblem() )
    {
      // For moving grid problems we generate an overlapping grid at the previous time

      CompositeGrid & cg = gf[prev].cg;
      
      MovingGrids & movingGrids = parameters.dbase.get<MovingGrids >("movingGrids");
      
      movingGrids.getPastTimeGrid( gf[prev] );

      printF(" --DS-- getPastTimeSolutions: REGENERATE THE PAST-TIME OVERLAPPING GRID  ---\n\n");
      parameters.regenerateOverlappingGrid( cg , cg, true );

    }

    // Assign the "initial" conditions 

    assignInitialConditions(prev);


    if( true && parameters.dbase.get<GUIState* >("runTimeDialog")!=NULL  )
    {
      // -- optionally plot the solution and grid --
      // optionIn: 0=wait, 1=plot-and-wait, 2=plot-but-don't-wait
      int optionIn = 1;
      real tFinal=tPast+1;
      plot(tPast, optionIn, tFinal, prev ); 
    }
    

  }
  

  return 0;
}
