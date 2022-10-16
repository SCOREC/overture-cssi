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

#include "incompressible/SmCylinderExactSolution.h"
#include "incompressible/SmSphereExactSolution.h"
#include "incompressible/SmRectangleExactSolution.h"

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

  aString caseName = "hollowCylinderDD"; // 
  // aString caseName = "hollowCylinderTT";

  
  int mx=2, my=1, mz=1; // mode numbers for rectangle

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
      else if( len=line.matches("-mx=") )
      {
        sScanF(line(len,line.length()-1),"%i",&mx);
        printF("Setting mx=%i\n",mx);   
      } 
      else if( len=line.matches("-my=") )
      {
        sScanF(line(len,line.length()-1),"%i",&my);
        printF("Setting my=%i\n",my);   
      } 
      else if( len=line.matches("-mz=") )
      {
        sScanF(line(len,line.length()-1),"%i",&mz);
        printF("Setting mz=%i\n",mz);   
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

  const bool evalCylinder = caseName=="hollowCylinderDD" || 
                            caseName=="hollowCylinderDD" ||
                            caseName=="solidCylinderD"   ||
                            caseName=="solidCylinderT";

  const bool evalSphere = caseName=="hollowSphereDD" || 
                          caseName=="hollowSphereDD" ||
                          caseName=="solidSphereD"   ||
                          caseName=="solidSphereT";

  const bool evalRectangle = caseName=="rectangleDD" || 
                             caseName=="rectangleTD" ||
                             caseName=="rectangleTT";

  printF("&&&& evalCylinder=%d, evalSphere=%d, evalRectangle=%d,\n",(int)evalCylinder,(int)evalSphere,(int)evalRectangle);

  assert( evalCylinder || evalSphere || evalRectangle );


  SmCylinderExactSolution cylExact;
  if( evalCylinder )
  {
    cylExact.initialize( cg, caseName );
  }

  SmSphereExactSolution sphereExact;
  if( evalSphere )
  {
    sphereExact.initialize( cg, caseName );  
  }

  SmRectangleExactSolution rectangleExact;
  if( evalRectangle )
  {
    rectangleExact.initialize( cg, caseName );  
    rectangleExact.setParameter( "mx",mx );
    rectangleExact.setParameter( "my",my );
    rectangleExact.setParameter( "mz",mz );    
  }  

  const int numberOfDimensions = cg.numberOfDimensions();
  int numberOfComponents= numberOfDimensions+1;
  Range all;
  realCompositeGridFunction u(cg,all,all,all,numberOfComponents);
  int mc=0; 
  u.setName("u1",mc); mc++; 
  u.setName("u2",mc); mc++; 
  if( numberOfDimensions>2 ){ u.setName("u3",mc); mc++; }
  u.setName("p", mc); mc++;

  u=0.;

  Index I1,I2,I3;
  Index Ib1,Ib2,Ib3;
  real t=0.;
  int numberOfTimeDerivatives = 0;
  for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
  {
     MappedGrid & mg = cg[grid];
     getIndex( mg.dimension(),I1,I2,I3);

     if( evalCylinder )
       cylExact.evalSolution( t, cg, grid, u[grid], I1,I2,I3, numberOfTimeDerivatives );
     else if( evalSphere )
       sphereExact.evalSolution( t, cg, grid, u[grid], I1,I2,I3, numberOfTimeDerivatives );
     else
       rectangleExact.evalSolution( t, cg, grid, u[grid], I1,I2,I3, numberOfTimeDerivatives );

  }

  Real omega;
  if( evalCylinder )
  {
    cylExact.getParameter( "omega",omega );
  }
  if( evalSphere )
  {
    sphereExact.getParameter( "omega",omega );
  } 
  if( evalRectangle )
  {
    rectangleExact.getParameter( "omega",omega );
  }   

  printF(">>>>> omega=%12.4e for the exact solution\n",omega);

  PlotStuffParameters psp;
  psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);
  psp.set(GI_TOP_LABEL,"Solution");
  PlotIt::contour( ps,u,psp );

  CompositeGridOperators cgop(cg);

  const int u1c=0; 
  const int u2c=1; 
  const int u3c=2; 
  const int pc = u1c + numberOfDimensions;

  Real mu = 1.; // ** FIX ME 

  Range Uc(u1c,u1c+numberOfDimensions-1);  // displacement components
  Range C(u1c,pc);                         // all components 

  // Put the residual here:
  realCompositeGridFunction res(cg,all,all,all,numberOfComponents+1);
  mc=0; 
  res.setName("u1",mc); mc++; 
  res.setName("u2",mc); mc++; 
  if( numberOfDimensions>2 ){ res.setName("u3",mc); mc++; }
  res.setName("p", mc); mc++;
  res.setName("div",mc);
  const int nDiv=mc; // position of div in res

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
    op.derivative( MappedGridOperators::laplacianOperator,uLocal, uLap,I1,I2,I3,C );

    RealArray px(I1,I2,I3), py(I1,I2,I3), pz(I1,I2,I3);
    op.derivative(MappedGridOperators::xDerivative,uLocal,px  ,I1,I2,I3,pc);
    op.derivative(MappedGridOperators::yDerivative,uLocal,py  ,I1,I2,I3,pc);

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
    if( numberOfDimensions==2 )
    {
      resLocal(I1,I2,I3,nDiv) = u1x + u2y;
    }
    else
    {
      op.derivative( MappedGridOperators::zDerivative,uLocal, u3z,I1,I2,I3,u3c );
      resLocal(I1,I2,I3,nDiv) = u1x + u2y + u3z;
    }

    printF("============= Check boundary conditions =============\n");
    ForBoundary(side,axis) 
    {
      if( mg.boundaryCondition(side,axis)>0 )
      {
        getBoundaryIndex(mg.gridIndexRange(),side,axis,Ib1,Ib2,Ib3);
        OV_GET_SERIAL_ARRAY(real,u[grid],uLocal);

        if( caseName == "hollowCylinderDD" || 
            caseName == "solidCylinderD"   ||
            caseName == "hollowSphereDD"   ||
            caseName == "solidSphereD"     ||
            caseName == "rectangleDD"      ||
           (caseName == "rectangleTD" && axis==0 && side==1 )
          )
        {
          // check displacement BCs

          Real u1Err = max(fabs(uLocal(Ib1,Ib2,Ib3,u1c)));
          Real u2Err = max(fabs(uLocal(Ib1,Ib2,Ib3,u2c)));
          Real u3Err=0;
          if( numberOfDimensions==3 ) 
            u3Err = max(fabs(uLocal(Ib1,Ib2,Ib3,u3c)));
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


  int movieMode=1;
  if( movieMode )
  {
    Real dt=.1; 
    int numSteps=100;
    for( int n=0; n<numSteps; n++ )
    {
      real t= n*dt;

      int numberOfTimeDerivatives = 0;
      for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
      {
         MappedGrid & mg = cg[grid];
         getIndex( mg.dimension(),I1,I2,I3);

         if( evalCylinder )
           cylExact.evalSolution( t, cg, grid, u[grid], I1,I2,I3, numberOfTimeDerivatives );
         else if( evalSphere )
           sphereExact.evalSolution( t, cg, grid, u[grid], I1,I2,I3, numberOfTimeDerivatives );
         else
           rectangleExact.evalSolution( t, cg, grid, u[grid], I1,I2,I3, numberOfTimeDerivatives );
      } 

      psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);
      psp.set(GI_TOP_LABEL,sPrintF("Solution: t=%5.2f",t));
      ps.erase();
      PlotIt::contour( ps,u,psp );
    }
  }

  Overture::finish();          
  return 0;
}
