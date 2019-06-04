#include "ReferenceSolution.h"

#include "ShowFileReader.h"
#include "InterpolatePointsOnAGrid.h"
#include "ParallelUtility.h"
#include "PlotStuff.h"
 
// =================================================================================
/// \brief Constructor.
// =================================================================================
ReferenceSolution::
ReferenceSolution()
{
  showFileReader = NULL;
  uTarget        = NULL;
  interpolatePointsOnAGrid = NULL;

  currentSolution=-1;
  currentTime    =-1e10;
  tPlot          =-1.;
}


// =================================================================================
/// \brief Destructor.
// =================================================================================
ReferenceSolution::
~ReferenceSolution()
{
  delete showFileReader;
  delete uTarget;
  delete [] interpolatePointsOnAGrid;
}

// =================================================================================
/// \brief Set the name of the reference show file.
// =================================================================================
int ReferenceSolution::
setShowFileName( const aString & name )
{
  nameOfShowFile=name;
  return 0;
}


// =================================================================================
/// \brief Return the reference solution at time t, interpolated onto grid cgTarget.
// =================================================================================
realCompositeGridFunction & 
ReferenceSolution::getSolution( real t, CompositeGrid & cgTarget, Range & C )
{
  int debug=1;
  
  const int interpolationWidth=max(cgTarget.interpolationWidth());
  
  if( debug &1 )
    printF("XXXXXXX ReferenceSolution::getSolution: numberOfDomains=%i, interpWidth=%i\n",cgTarget.numberOfDomains(),
           interpolationWidth);
  
  const real tEps=REAL_EPSILON*100.;
  if( fabs(t-currentTime)/max(1.,t) < tEps )
  {
    // We have already found this solution 
    assert( uTarget != NULL );
    return *uTarget;
  }

  CompositeGrid cgRef;              // hold reference grid from the show file 
  realCompositeGridFunction uRef;   // hold reference solution from the show file 

  if( showFileReader==NULL )
  {
    showFileReader = new ShowFileReader(nameOfShowFile);

    const int numberOfFrames = showFileReader->getNumberOfFrames();
    const int numberOfSolutions = showFileReader->getNumberOfSolutions();
    printF("ReferenceSolution::INFO: open show file=[%s], numberOfFrames=%i, numberOfSolutions=%i\n",
           (const char*)nameOfShowFile,numberOfFrames,numberOfSolutions);

    // Compute "tPlot" :  the time between solutions in the show file.
    real time[2];
    for( int solution=0; solution<2; solution++ )
    {
      // Fine the times for the first two frames and subtract to get tPlot
      int solutionNumber = solution+1;
      // showFileReader->getASolution(solutionNumber,cgRef,uRef);  // show file solutions are base 1 
      HDF_DataBase & db = *(showFileReader->getFrame(solution+1)); 
      db.get(time[solution],"time");  
    }
    tPlot = time[1]-time[0];
    printF("ReferenceSolution:: time[0]=%9.3e, time[1]=%9.3e : tPlot=%9.3e\n",time[0],time[1],tPlot);
   
  }

	    
  int solutionNumber = 1 + int( t/tPlot + .5);   // This is our guess for the reference solution we need
    
  printF(" showFileReader->getASolution ...\n");
  showFileReader->getASolution(solutionNumber,cgRef,uRef);        // read in a grid and solution
  printF(" showFileReader->getASolution ...done\n");

  HDF_DataBase & db = *(showFileReader->getFrame(solutionNumber));
  real t0=-1.;
  db.get(t0,"time");
  if( debug & 1 )
    printF(" **** ReferenceSolution: t=%9.3e, get solutionNumber=%i, timeFromFile=%9.3e\n",t,solutionNumber,t0);
  if( fabs(t-t0)/max(1.,t) > tEps )
  {
    printF("\n XXXXXXXXX ReferenceSolution::ERROR: times do NOT match! t=%9.3e, timeFromFile=%9.3e, diff=%9.3e  XXXXXXXXXX \n",
           t,t0,t-t0);
    printF(" XXX comparing: fabs(t-t0)/max(1.,t)=%9.3e >? tEps=%9.3e\n",fabs(t-t0)/max(1.,t),tEps);
    
  }
  currentTime=t0;
  



  // The reference solution may live on a larger domain so we interpolate to the targetDomain
  Range all;
  uTarget = new realCompositeGridFunction(cgTarget,all,all,all,C);
  for( int c=C.getBase(); c<=C.getBound(); c++ )
  {
    uTarget->setName(sPrintF("E%i",c),c);
  }


  cgTarget.update(MappedGrid::THEmask | GridCollection::THEdomain  );

  cgRef.update(MappedGrid::THEmask | GridCollection::THEdomain  );
    
  assert( cgTarget.numberOfDomains() == cgRef.numberOfDomains() );

  // We need to keep track of domains so we only interpolate from the same domain
  const int numberOfDomains = cgTarget.numberOfDomains();

  if( interpolatePointsOnAGrid==NULL )
  {
    interpolatePointsOnAGrid = new  InterpolatePointsOnAGrid [numberOfDomains];
    for( int domain=0; domain<numberOfDomains; domain++ )
    {
      
      InterpolatePointsOnAGrid & interpolator = interpolatePointsOnAGrid[domain];

      interpolator.setInfoLevel( 1 );
      interpolator.setInterpolationWidth(interpolationWidth);
      // Set the number of valid ghost points that can be used when interpolating from a grid function: 
      int numGhostToUse=1;
      interpolator.setNumberOfValidGhostPoints( numGhostToUse );
      
      // Assign all points, extrapolate pts if necessary:
      interpolator.setAssignAllPoints(true);
    }
    
  }
  real time0=getCPU();
  for( int domain=0; domain<numberOfDomains; domain++ )
  {
    // --- Interpolate domain by domain ---- (code from comp.C)

    bool gridsMatch = false; // domain==0;   // **FIX ME ***

    if( gridsMatch )
    {
      printF("Just copy solution for domain=%i\n",domain);
      for( int grid=0; grid<cgRef.numberOfComponentGrids(); grid++ )
      {
        if( cgRef.domainNumber(grid)==domain )
        {
          (*uTarget)[grid] = uRef[grid];
        }
        
      }
      
    }
    else
    {

      CompositeGrid & cgRefd    = cgRef.domain[domain];
      CompositeGrid & cgTargetd = cgTarget.domain[domain];

      InterpolatePointsOnAGrid & interpolator = interpolatePointsOnAGrid[domain];

      // compute the master grid number for ref grid: 
      //   cgRefd[grid] = cgRef[masterGrid(k)]
      IntegerArray masterGrid(cgRefd.numberOfComponentGrids());
      int k=0;
      for( int grid=0; grid<cgRef.numberOfComponentGrids(); grid++ )
      {
        if( cgRef.domainNumber(grid)==domain )
        {
          masterGrid(k)=grid; k++;   // this grid in master collection is in the current domain
        }
      }
      assert( k==cgRefd.numberOfComponentGrids() );


      realCompositeGridFunction uRefd(cgRefd,all,all,all,C);       uRefd=0.;
      realCompositeGridFunction uTargetd(cgTargetd,all,all,all,C); uTargetd=0.; 
	
      for( int grid=0; grid<cgRefd.numberOfComponentGrids(); grid++ )
      {
        // printF("domain=%i: Set uRefd[%i]=uRef[%i]\n",domain,grid,masterGrid(grid));
      
        uRefd[grid]=uRef[masterGrid(grid)];  // copy master to domain grid-function
      }


      int numGhost=0;  // number of ghost to interpolate ** FIX ME **

      int num=interpolator.interpolateAllPoints( uRefd,uTargetd,C,C,numGhost);    // interpolate uTargetd from uRefd


      masterGrid.redim(cgTargetd.numberOfComponentGrids());
      k=0;
      for( int grid=0; grid<cgTarget.numberOfComponentGrids(); grid++ )
      {
        if( cgTarget.domainNumber(grid)==domain )
        {
          masterGrid(k)=grid; k++;   // this grid in master collection is in the current domain
        }
      }
      assert( k==cgTargetd.numberOfComponentGrids() );

      for( int grid=0; grid<cgTargetd.numberOfComponentGrids(); grid++ )
      {
        // printF("domain=%i: Set uTarget[%i]=uTargetd[%i]\n",domain,masterGrid(grid),grid);

        (*uTarget)[masterGrid(grid)]=uTargetd[grid];  // copy interpolated values into master 
      }
    }
    

  } // end for domain

  
  real time=getCPU()-time0;
  time=ParallelUtility::getMaxValue(time);

  if( debug & 1 )
    printF("ReferenceSolution:: time to interpolate to the target grid = %8.2e(s)\n",time);

  if( false )
  {
    uTarget->interpolate();
  }
  

  if( false )
  {
    // plot the solutions 
    GenericGraphicsInterface & gi = *Overture::getGraphicsInterface();
    PlotStuffParameters psp;
    psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);
    aString buff;

    gi.erase();
    psp.set(GI_TOP_LABEL,sPrintF(buff,"uRef (from show file) t=%9.3e",t));
    PlotIt::contour(gi,uRef,psp);  

    gi.erase();
    psp.set(GI_TOP_LABEL,sPrintF(buff,"uTarget (interpolated) t=%9.3e",t));
    PlotIt::contour(gi,*uTarget,psp);  
  }
  

  return *uTarget;

}

