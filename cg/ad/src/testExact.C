// ---------------------------------------------------------------------------
//
// Test exact solutions for testExact
//
// ---------------------------------------------------------------------------

#include "CompositeGridOperators.h"
#include "PlotStuff.h"
// #include "SquareMapping.h"
// #include "AnnulusMapping.h"
// #include "MatrixTransform.h"
// #include "DataPointMapping.h"

// #include "OGTrigFunction.h"
// #include "OGPolyFunction.h"
#include "display.h"

#include "AdExactSolutions.h"

#include "ParallelUtility.h"
#include "Oges.h"
#include "CgSolverUtil.h"

#include "gridFunctionNorms.h"

#define ForBoundary(side,axis)   for( int axis=0; axis<mg.numberOfDimensions(); axis++ ) \
                                 for( int side=0; side<=1; side++ )
int 
getLineFromFile( FILE *file, char s[], int lim);

// void display(realArray & u )
// {
//   printF("u.getlength(0)=%i\n",u.getLength(0));
  
//   ::display(u,"u");
// }



int 
main(int argc, char *argv[])
{
  Overture::start(argc,argv);
  // Use this to avoid un-necessary communication: 
  Optimization_Manager::setForceVSG_Update(Off);
  const int myid=Communication_Manager::My_Process_Number;

  // This macro will initialize the PETSc solver if OVERTURE_USE_PETSC is defined.
  INIT_PETSC_SOLVER();

  int plotOption=true;
  bool smartRelease=false;
  bool reportMemory=false;
  bool loadBalance=false;
  int numberOfParallelGhost=2;

  aString caseName = "doubleAnnulusTT"; // 
  // aString caseName = "hollowCylinderTT";

  
  aString commandFileName="";
  if( argc > 1 )
  { // look at arguments for "noplot" or some other name
    int len=0;
    aString line;
    for( int i=1; i<argc; i++ )
    {
      line=argv[i];
      if( line=="-noplot" || line=="noplot" )
        plotOption=false;
      else if( line=="-nopause" || line=="-abortOnEnd" || line=="-nodirect" ||
               line=="-readCollective" || line=="-writeCollective" ||
               line=="nopause" || line=="abortOnEnd" || line=="nodirect" )
        continue; // these commands are processed by getGraphicsInterface below 
      else if( line=="memory" )
      {
        reportMemory=true;
        Diagnostic_Manager::setTrackArrayData(TRUE);
      }
      else if( line=="loadBalance" || line=="-loadBalance" )
      {
        loadBalance=true;
      }
      else if( len=line.matches("-numberOfParallelGhost=") )
      {
        sScanF(line(len,line.length()-1),"%i",&numberOfParallelGhost);
        if( numberOfParallelGhost<0 || numberOfParallelGhost>10 )
        {
          printF("ERROR: numberOfParallelGhost=%i is no valid!\n",numberOfParallelGhost);
          OV_ABORT("error");
        }
        printF("Setting numberOfParallelGhost=%i\n",numberOfParallelGhost);
      }
      else if( len=line.matches("-caseName=") )
      {
        caseName = line(len,line.length()-1);
        printF("Setting caseName=[%s]\n",(const char*)caseName);
      }      
      else if( line=="release" )
      {
        smartRelease=true;
        printF("*** turn on smart release of memory ***\n");
        Diagnostic_Manager::setSmartReleaseOfInternalMemory( ON );
      }
      else if( commandFileName=="" )
      {
        commandFileName=line;    
        printF("testExact: reading commands from file [%s]\n",(const char*)commandFileName);
      }
      
    }
  }
  else
    printF("Usage: `testExact [options][file.cmd]' \n"
            "     options:                            \n" 
            "          -noplot:   run without graphics \n" 
            "          -caseName: defines the exact solution \n" 
            "          -abortOnEnd: abort if command file ends \n" 
            "          -numberOfParallelGhost=<num> : number of parallel ghost lines \n" 
            "          memory:   run with A++ memory tracking\n" 
            "          release:  run with A++ smart release of memory\n"
            "     file.cmd: read this command file \n");

  GenericGraphicsInterface & ps = *Overture::getGraphicsInterface("testExact",false,argc,argv);

  // By default start saving the command file called "testExact.cmd"
  aString logFile="testExact.cmd";
  ps.saveCommandFile(logFile);
  printF("User commands are being saved in the file `%s'\n",(const char *)logFile);

  ps.appendToTheDefaultPrompt("testExact>");

  // read from a command file if given
  if( commandFileName!="" )
  {
    printF("read command file =%s\n",(const char*)commandFileName);
    ps.readCommandFile(commandFileName);
  }


  CompositeGrid cg;
  const int maxWidthExtrapInterpNeighbours=4;  // This means we support 3rd-order extrap, (1,-3,3,-1)
  aString nameOfGridFile="";
  nameOfGridFile = readOrBuildTheGrid(ps, cg, loadBalance, numberOfParallelGhost,maxWidthExtrapInterpNeighbours );


  // Interpolant & interpolant = *new Interpolant(cg); interpolant.incrementReferenceCount();
  // interpolant.setImplicitInterpolationMethod(Interpolant::iterateToInterpolate);

  const bool evalDoubleAnnulus = caseName=="doubleAnnulusCase4" ||
                                 caseName=="doubleAnnulusCase5";

  const bool evalDoubleRectangle = caseName=="doubleRectangleCase1" ||
                                   caseName=="doubleRectangleCase";                                 

  // const bool evalSphere = caseName=="hollowSphereDD" || 
  //                         caseName=="hollowSphereDD" ||
  //                         caseName=="solidSphereD"   ||
  //                         caseName=="solidSphereT";

  // printF("&&&& evalCylinder=%d, evalSphere=%d\n",(int)evalCylinder,(int)evalSphere);

  aString exactSolution;
  if( evalDoubleAnnulus )
  {
     exactSolution="doubleAnnulus";
  } 
  else if( evalDoubleRectangle )
  {
    exactSolution = "doubleRectangle";
  }
  else
  {
    printF("testExact: unknown caseName=[%s]\n",(const char*)caseName);
    OV_ABORT("error");
  }

  cg.update(GridCollection::THEdomain);

  AdExactSolutions daExact[2];
  for( int d=0; d<cg.numberOfDomains(); d++ )
  {
    // int domain=1-d; // NOTE: 
    // daExact[d].setDomain( domain );

    daExact[d].initialize( cg.domain[d], exactSolution, caseName );
  }

  int numberOfComponents=1;
  Range all;
  realCompositeGridFunction u(cg,all,all,all,numberOfComponents);
  u.setName("u",0);
  u=0.;

  Index I1,I2,I3;
  Index Ib1,Ib2,Ib3;
  real t=0.;
  int numberOfTimeDerivatives = 0;
  for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
  {
     MappedGrid & mg = cg[grid];
     getIndex( mg.dimension(),I1,I2,I3);

     const int d = cg.domainNumber(grid); // this grid is in this domain
     int gridInDomain =0; // fix me 
     if( evalDoubleAnnulus )
       daExact[d].evalSolution( t, cg.domain[d], gridInDomain, u[grid], I1,I2,I3, numberOfTimeDerivatives );

  }


  // for( int d=0; d<cg.numberOfDomains(); d++ )
  // {
  //   CompositeGrid & cgd = cg.domain[d];
  //   realCompositeGridFunction & ud = u.domain[d];

  //   for( int grid=0; grid<cgd.numberOfComponentGrids(); grid++ )
  //   {
  //      MappedGrid & mg = cgd[grid];
  //      getIndex( mg.dimension(),I1,I2,I3);

  //      if( evalDoubleAnnulus )
  //        daExact[d].evalSolution( t, cgd, grid, ud[grid], I1,I2,I3, numberOfTimeDerivatives );

  //   }
  // }

  // Real sr,si;
  // if( evalDoubleAnnulus )
  // {
  //   daExact.getParameter( "sr",sr );
  // }

  // printF(">>>>> [sr,si]=[%12.4e,%12.4e] for the exact solution\n",sr,si);

  PlotStuffParameters psp;
  psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);
  psp.set(GI_TOP_LABEL,"Solution");
  PlotIt::contour( ps,u,psp );

/* ---
  CompositeGridOperators cgop(cg);

  const int numberOfDimensions = cg.numberOfDimensions();
  const int uc=0;

  // Put the residual here:
  realCompositeGridFunction res(cg,all,all,all,numberOfComponents);
  res.setName("u",0);

  for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
  {
    MappedGrid & mg = cg[grid];
    MappedGridOperators & op = cgop[grid];

    mg.update(MappedGrid::THEvertex | MappedGrid::THEcenter | MappedGrid::THEinverseVertexDerivative );

    OV_GET_SERIAL_ARRAY(real,u[grid],uLocal);
    OV_GET_SERIAL_ARRAY(real,res[grid],resLocal);

    resLocal = 0.;

    int extra=-1; // skip boundaries -- *fix* me for traction
    getIndex( mg.gridIndexRange(),I1,I2,I3,extra ); 
     
    RealArray uLap(I1,I2,I3,C);
    op.derivative( MappedGridOperators::laplacianOperator,uLocal, uLap,I1,I2,I3,uc );

    RealArray ux(I1,I2,I3), uy(I1,I2,I3), uz(I1,I2,I3);
    op.derivative(MappedGridOperators::xDerivative,uLocal,ux  ,I1,I2,I3,uc);
    op.derivative(MappedGridOperators::yDerivative,uLocal,uy  ,I1,I2,I3,uc);

    resLocal(I1,I2,I3,u1c) = SQR(omega)*uLocal(I1,I2,I3,u1c) + mu*uLap(I1,I2,I3,u1c) - px(I1,I2,I3); 
    resLocal(I1,I2,I3,u2c) = SQR(omega)*uLocal(I1,I2,I3,u2c) + mu*uLap(I1,I2,I3,u2c) - py(I1,I2,I3); 
    if( numberOfDimensions==3 )
    {
     op.derivative(MappedGridOperators::zDerivative,uLocal,pz  ,I1,I2,I3,pc);
     resLocal(I1,I2,I3,u3c) = SQR(omega)*uLocal(I1,I2,I3,u3c) + mu*uLap(I1,I2,I3,u3c) - pz(I1,I2,I3); 
    }

    resLocal(I1,I2,I3,pc) = uLap(I1,I2,I3,pc); 

    RealArray u1x(I1,I2,I3), u2y(I1,I2,I3), u3z(I1,I2,I3);
    op.derivative( MappedGridOperators::xDerivative,uLocal, u1x,I1,I2,I3,u1c );
    op.derivative( MappedGridOperators::yDerivative,uLocal, u2y,I1,I2,I3,u2c );
    op.derivative( MappedGridOperators::zDerivative,uLocal, u3z,I1,I2,I3,u3c );

    resLocal(I1,I2,I3,nDiv) = u1x + u2y + u3z;

    ForBoundary(side,axis) 
    {
      if( mg.boundaryCondition(side,axis)>0 )
      {
        getBoundaryIndex(mg.gridIndexRange(),side,axis,Ib1,Ib2,Ib3);
        OV_GET_SERIAL_ARRAY(real,u[grid],uLocal);

        if( caseName == "hollowCylinderDD" || 
            caseName == "solidCylinderD"   ||
            caseName == "hollowSphereDD"   ||
            caseName == "solidSphereD"   )
        {
          // check displacement BCs

          Real u1Err = max(fabs(uLocal(Ib1,Ib2,Ib3,u1c)));
          Real u2Err = max(fabs(uLocal(Ib1,Ib2,Ib3,u2c)));
          Real u3Err = max(fabs(uLocal(Ib1,Ib2,Ib3,u3c)));
          Real dispErr = max(fabs(uLocal(Ib1,Ib2,Ib3,Uc)));
          printF(" grid=%d, (side,axis)=(%d,%d), error in u=0 BC: u1=%8.2e, u2=%8.2e, u3=%8.2e, max=%8.2e\n",
               grid,side,axis,u1Err,u2Err,u3Err,dispErr);
        }
      }
    }
  }

  // compute the residual norms
  printF("============= caseName = %s, grid=%s =============\n",(const char*)caseName,(const char*)nameOfGridFile);
  printF("============= residual norms: omega=%12.6e =============\n",omega);
  for( int m=0; m<=numberOfComponents; m++ )
  {
    Real maxRes = maxNorm(res,m);
    Real l2Res = l2Norm(res,m);

    printF("testExact: equation %d:  max-res=%9.3e, l2-res=%9.3e  (%s)\n",m,maxRes,l2Res,(const char*)res.getName(m));

  }


  // --- plot the residuals ----
  ps.erase();
  psp.set(GI_TOP_LABEL,"Residuals");
  psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);
  PlotIt::contour( ps,res,psp );

--- */

  Overture::finish();          
  return 0;
}
