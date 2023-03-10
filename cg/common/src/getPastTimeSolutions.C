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

  printF("\n --DS-- GET-PAST-TIME SOLUTIONS getPastTimeSolution: current=%i (t=%8.2e) numberOfPast=%i\n",
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

      printF("\n --DS-- getPastTimeSolutions: REGENERATE THE PAST-TIME OVERLAPPING GRID t=%8.2e ---\n",tPast);
      if( debug() & 4 )
      {
        FILE *& debugFile =parameters.dbase.get<FILE* >("debugFile");
        fPrintF(debugFile," --DS-- getPastTimeSolutions: REGENERATE THE PAST-TIME OVERLAPPING GRID  t=%9.3e---\n\n",
                gf[prev].t );
      }
      
      parameters.regenerateOverlappingGrid( cg , cg, true );

      if( debug() & 64 )
      {
        FILE *& debugFile =parameters.dbase.get<FILE* >("debugFile");
        ::displayMask(cg[0].mask(),"Past time grid - mask on grid 0",debugFile);
      }

      if( true )
      {
        // -- *new* Jan 15, 2019 **  who computes the grid velocity ???
        printF("Compute grid-velocity at past time t=%9.3e\n",gf[prev].t);
        getGridVelocity( gf[prev], gf[prev].t );

        // ::display(gf[prev].getGridVelocity(1),"gridVelocity","%8.1e ");
      }
    
    } // end if movingGrid
    

    // -- Assign the "initial" conditions --
    assignInitialConditions(prev);

    if( debug() & 8 ) 
    {
      gf[prev].u.display(sPrintF("Past time solution prev=%i t=%9.3e BEFORE PROJECT \n",prev,gf[prev].t),"%6.2f ");
    }
    


    // -- compute the pressure on moving grids when the pressure and body accelerations are coupled --
    if( true ) // parameters.dbase.get<bool>("projectInitialConditions") )
    {
      projectInitialConditionsForMovingGrids(prev);

      if( debug() & 8) 
      {
        gf[prev].u.display(sPrintF("Past time solution prev=%i t=%9.3e AFTER PROJECT \n",prev,gf[prev].t),"%6.2f ");
      }
    }
    
    if( (false || debug() & 16)  && parameters.dbase.get<GUIState* >("runTimeDialog")!=NULL  )
    {
      // -- optionally plot the solution and grid --
      // optionIn: 0=wait, 1=plot-and-wait, 2=plot-but-don't-wait
      int optionIn = 1;
      real tFinal=tPast+1;
      printF(" --DS-- getPastTimeSolutions: plot past time solution at t=%9.3e\n",tPast);
      
      plot(tPast, optionIn, tFinal, prev ); 
    }
    

  }
  

  return 0;
}

