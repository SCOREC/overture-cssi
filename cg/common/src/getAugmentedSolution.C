#include "DomainSolver.h"
#include "CompositeGridOperators.h"
#include "ParallelUtility.h"
#include "App.h"
#include "GenericGraphicsInterface.h"

// int 
// getAmrErrorFunction(realCompositeGridFunction & u, 
//                     real t,
//                     Parameters & parameters,
//                     realCompositeGridFunction & error,
//                     bool computeOnFinestLevel =false );


//\begin{>>CompositeGridSolverInclude.tex}{\subsection{getAugmentedSolution}} 
realCompositeGridFunction & DomainSolver::
getAugmentedSolution( GridFunction & gf0, realCompositeGridFunction & v )
// ========================================================================================
// /Description:
//    Return a grid function for plotting and for output that may contain extra variables, such
// as errors or such as the pressure for the compressible NS.
// /gf0 (input) : input grid function.
// /v (input) : grid function to hold the result, IF extra variables are added.
// /Return values: The possibly augmented solution.
//\end{CompositeGridSolverInclude.tex}  
// ========================================================================================
{
  checkArrayIDs(sPrintF("getAugmentedSolution: start") );

  // No need to compute the augmented solution if graphics plotting is off
  if( (parameters.dbase.get<GenericGraphicsInterface* >("ps")!=NULL && 
      !parameters.dbase.get<GenericGraphicsInterface* >("ps")->graphicsIsOn()) ||
      parameters.dbase.get<int>("simulateGridMotion")>0  )
  {
    return gf0.u;
  }
  // if( parameters.dbase.get<GenericGraphicsInterface* >("ps")!=NULL && !parameters.dbase.get<GenericGraphicsInterface* >("ps")->graphicsIsOn() ) return gf0.u;

  const int numberOfComponents = parameters.dbase.get<int >("numberOfComponents");
  const bool & twilightZoneFlow = parameters.dbase.get<bool >("twilightZoneFlow"); 

  CompositeGrid & cg = gf0.cg;
  realCompositeGridFunction & u = gf0.u;

  bool plotMoreComponents = false;
  int extra=0;
  int dc=0;   // save divergence here
  int kc=0;   // save known solution here 
  if( parameters.isAdaptiveGridProblem() && parameters.dbase.get<int >("showAmrErrorFunction") )
    extra+=1;


  if( parameters.dbase.get<int >("showResiduals") || parameters.dbase.get<Parameters::TimeSteppingMethod >("timeSteppingMethod")==Parameters::steadyStateRungeKutta )
  {
    // plot residuals too
    extra+=numberOfComponents;
    dc=2*numberOfComponents;
  }
  
  bool plotErrorsInKnownSolution = parameters.dbase.get<Parameters::KnownSolutionsEnum >("knownSolution")!=Parameters::noKnownSolution;
  // plotErrorsInKnownSolution=false; // *** TEMP ****
  // printF("$$$$$$$$$$$$$ getAugmentedSolution: plotErrorsInKnownSolution = %d $$$$$$$$$$$$$$$\n",(int)plotErrorsInKnownSolution);

  if( false )
    printP(" ###### DomainSolver::getAugmentedSolution at t=%9.2e, plotErrorsInKnownSolution=%d #####\n",gf0.t,plotErrorsInKnownSolution);

  if( plotErrorsInKnownSolution )
  {
    assert( !twilightZoneFlow );
    kc = numberOfComponents;
    dc = 2*numberOfComponents;
  }


  if( twilightZoneFlow || plotErrorsInKnownSolution )
  {
    Range all;
    v.updateToMatchGrid(cg,all,all,all,numberOfComponents*2+extra); 
    for( int n=0; n<numberOfComponents; n++ )
    {
      v.setName(u.getName(n),n);
      v.setName(u.getName(n)+" (error)",n+numberOfComponents);
    }
  }

  if( twilightZoneFlow )
  {
    // For twilight zone flow we also plot the errors
    plotMoreComponents=true;
    dc+=numberOfComponents;
    
    Range all;
    Index I1,I2,I3;
    // v.updateToMatchGrid(cg,all,all,all,numberOfComponents*2+extra);  // extra space for errors
    Range N(0,numberOfComponents-1);
    Range N2(numberOfComponents,2*numberOfComponents-1);

    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
      MappedGrid & mg = cg[grid];
      getIndex(mg.dimension(),I1,I2,I3);

      realArray & x= mg.center();
      OV_GET_SERIAL_ARRAY(real,x,xLocal);
      OV_GET_SERIAL_ARRAY(real,u[grid],uLocal);
      OV_GET_SERIAL_ARRAY(real,v[grid],vLocal);
      const int includeGhost=1;
      bool ok = ParallelUtility::getLocalArrayBounds(u[grid],uLocal,I1,I2,I3,includeGhost); 

      // bool ok=true;
      // #ifdef USE_PPP
      //   realSerialArray xLocal; getLocalArrayWithGhostBoundaries(x,xLocal);
      //   realSerialArray uLocal; getLocalArrayWithGhostBoundaries(u[grid],uLocal);
      //   realSerialArray vLocal; getLocalArrayWithGhostBoundaries(v[grid],vLocal);
      //   const int includeGhost=1;
      //   ok = ParallelUtility::getLocalArrayBounds(u[grid],uLocal,I1,I2,I3,includeGhost); 
      // #else
      //   const realSerialArray & xLocal = x;
      //   const realSerialArray & uLocal = u[grid];
      //   const realSerialArray & vLocal = v[grid];
      // #endif

      if( ok )
      {
        vLocal(I1,I2,I3,N)=uLocal(I1,I2,I3,N);

        bool isRectangular=false;
        realSerialArray u0(I1,I2,I3);
        for( int n=N.getBase(), n2=N2.getBase(); n<=N.getBound(); n++, n2++ )
        {
          parameters.dbase.get<OGFunction* >("exactSolution")->gd( u0,xLocal,mg.numberOfDimensions(),isRectangular,0,0,0,0,I1,I2,I3,n,gf0.t);
          vLocal(I1,I2,I3,n2)=vLocal(I1,I2,I3,n)-u0;
        }
      }

//       v[grid](I1,I2,I3,N)=u[grid](I1,I2,I3,N);
//       v[grid](I1,I2,I3,N2)=v[grid](I1,I2,I3,N)-(*parameters.dbase.get<OGFunction* >("exactSolution"))(cg[grid],I1,I2,I3,N,gf0.t);

    }

  }

  if( plotErrorsInKnownSolution )
  {     
    // printF("*************** getAugmented: SAVE ERRORS IN KNOWN *********************\n");
    if( !parameters.dbase.has_key("errorGridFunction") )
    {
      printF("getAugmentedSolution: known solution is there but no errors have been computed for plotting !?\n");

      // added *wdh* Nov 27, 2021 for cgmp so that it includes error components in option menu
      printF("getAugmentedSolution: I am calling determineErrors ....\n");
      determineErrors( gf0 );  
    } 

    if( true )
    {
      plotMoreComponents=true;

      realCompositeGridFunction & error = parameters.dbase.get<realCompositeGridFunction>("errorGridFunction");
      Range N=numberOfComponents, all; 
      for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
      {
        OV_GET_SERIAL_ARRAY(real,error[grid],errLocal);
        OV_GET_SERIAL_ARRAY(real,v[grid],vLocal);
        OV_GET_SERIAL_ARRAY(real,u[grid],uLocal);
        vLocal(all,all,all,N)    = uLocal(all,all,all,N);
        vLocal(all,all,all,N+kc) = errLocal(all,all,all,N);
      }
    }
  }
  
  if( parameters.isAdaptiveGridProblem() && parameters.dbase.get<int >("showAmrErrorFunction") )
  {
    // add on the AMR error estimate
    Range all;
    if( !plotMoreComponents )
    {
      v.updateToMatchGrid(cg,all,all,all,numberOfComponents+extra); 
      plotMoreComponents=true;
    }
    
    int ec = v.getComponentBound(0);
    #ifndef USE_PPP
      realCompositeGridFunction error;
      error.link(v,Range(ec,ec));
      bool computeOnFinestLevel=true;
      getAmrErrorFunction(u,gf0.t,error,computeOnFinestLevel);
    #else
      realCompositeGridFunction error(cg);
      bool computeOnFinestLevel=true;
      getAmrErrorFunction(u,gf0.t,error,computeOnFinestLevel);
      for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
      {
        assign(v[grid],all,all,all,ec, u[grid],all,all,all, 0 );
      }
      
    #endif

    v.setName("error estimate",ec);
  }
  
  if( parameters.dbase.get<int >("showResiduals") || parameters.dbase.get<Parameters::TimeSteppingMethod >("timeSteppingMethod")==Parameters::steadyStateRungeKutta )
  {
    if( !plotMoreComponents )
    {
      Range all;
      v.updateToMatchGrid(cg,all,all,all,numberOfComponents+extra); 
      plotMoreComponents=true;

      Range N=numberOfComponents;
      for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        v[grid](all,all,all,N)=u[grid](all,all,all,N);    // copy existing components
      for( int n=0; n<numberOfComponents; n++ )
        v.setName(u.getName(n),n);
    }

    Range N(0,numberOfComponents-1);
    int offset=numberOfComponents;
    if( parameters.dbase.get<bool >("twilightZoneFlow") )
      offset+=numberOfComponents;

    realCompositeGridFunction & residual = fn[0];   // pointer to du/dt

    // printf(" getAug: *** residual.numberOfComponentGrids=%i **\n",residual.numberOfComponentGrids());
    // printf(" getAug: *** v.numberOfComponentGrids=%i **\n",v.numberOfComponentGrids());
    
    Index I1,I2,I3;
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
      getIndex(cg[grid].dimension(),I1,I2,I3);
      v[grid](I1,I2,I3,N+offset)=residual[grid](I1,I2,I3,N);
      where( cg[grid].mask()(I1,I2,I3)<=0 )
      {
        for( int n=N.getBase(); n<=N.getBound(); n++ )
          v[grid](I1,I2,I3,n+offset)=0.;   // set residual to zero at unused and interp points
      }
      
    }
    for( int n=0; n<numberOfComponents; n++ )
    {
      v.setName(u.getName(n)+" (residual)",n+offset);
    }
  }
  
  checkArrayIDs(sPrintF("getAugmentedSolution: near end") ); 

  realCompositeGridFunction & ur = plotMoreComponents ? v : gf[current].u;
//   if( true ) // *wdh* 060228
//   {
//     gf0.u.display("\n\n getAugmentedSolution: gf0.u","%5.2f ");
//     ur.display("\n\n getAugmentedSolution: ur","%5.2f ");
//   }
  return ur;
}
