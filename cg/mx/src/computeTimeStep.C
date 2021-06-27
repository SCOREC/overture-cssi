#include "Maxwell.h"
#include "CompositeGridOperators.h"
#include "display.h"
#include "UnstructuredMapping.h"
#include "ParallelUtility.h"
#include "GridStatistics.h"
#include "DispersiveMaterialParameters.h"

#include "ULink.h"

extern bool verifyUnstructuredConnectivity( UnstructuredMapping &umap, bool verbose );

#define FOR_3D(i1,i2,i3,I1,I2,I3) \
int I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase();  \
int I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(); \
for(i3=I3Base; i3<=I3Bound; i3++) \
for(i2=I2Base; i2<=I2Bound; i2++) \
for(i1=I1Base; i1<=I1Bound; i1++)

#define FOR_3(i1,i2,i3,I1,I2,I3) \
I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase();  \
I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(); \
for(i3=I3Base; i3<=I3Bound; i3++) \
for(i2=I2Base; i2<=I2Bound; i2++) \
for(i1=I1Base; i1<=I1Bound; i1++)

// ======================================================================================
//  Return the max value of a scalar over all processors
//  /processor: return the result to this processor (-1 equals all processors)
// ======================================================================================
// real Maxwell::
// getMaxValue(real value, int processor /* = -1 */)
// {
//   real maxValue=value;
//   #ifdef USE_PPP 
//   if( processor==-1 )
//     MPI_Allreduce(&value, &maxValue, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
//   else
//     MPI_Reduce        (&value, &maxValue, 1, MPI_DOUBLE, MPI_MAX, processor, MPI_COMM_WORLD);
//   #endif
//   return maxValue;
// }

// int Maxwell::
// getMaxValue(int value, int processor /* = -1 */)
// {
//   int maxValue=value;
//   #ifdef USE_PPP 
//   if( processor==-1 )
//     MPI_Allreduce(&value, &maxValue, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
//   else
//     MPI_Reduce        (&value, &maxValue, 1, MPI_INT, MPI_MAX, processor, MPI_COMM_WORLD);
//   #endif
//   return maxValue;
// }

// real Maxwell::
// getMinValue(real value, int processor /* = -1 */ )
// {
//   real minValue=value;
//   #ifdef USE_PPP 
//   if( processor==-1 )
//     MPI_Allreduce(&value, &minValue, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
//   else
//     MPI_Reduce        (&value, &minValue, 1, MPI_DOUBLE, MPI_MIN, processor, MPI_COMM_WORLD);
//   #endif
//   return minValue;
// }

// int Maxwell::
// getMinValue(int value, int processor /* = -1 */)
// {
//   int minValue=value;
//   #ifdef USE_PPP 
//   if( processor==-1 )
//     MPI_Allreduce(&value, &minValue, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
//   else
//     MPI_Reduce        (&value, &minValue, 1, MPI_INT, MPI_MIN, processor, MPI_COMM_WORLD);
//   #endif
//   return minValue;
// }

// =============================================================================================
/// \brief Determine the time step
// =============================================================================================
int Maxwell::
computeTimeStep()
{
  if( method==bamx )
    return computeTimeStepBA();
  

  real time0=getCPU();

  assert( cgp!=NULL );
  CompositeGrid & cg= *cgp;
  const int numberOfComponentGrids = cg.numberOfComponentGrids();
  const int numberOfDimensions = cg.numberOfDimensions();

  const real & dtMax = dbase.get<real>("dtMax");

  deltaT=REAL_MAX*.01;
  
//    RealArray dtGrid(numberOfComponentGrids);  // time step for each grid by itself
//    dtGrid=REAL_MAX;
  
  // ====== SOSUP STABILITY REGIONS =======
  // SOSUP: dt depends on the order of accuracy
  // Approximate stability regions:
  //      (c*dt/dx)^sp + (c*dt/dy)^sp = lambda^sp
  // Then
  //     dt =  (lambda/c) / [  (1/dx)^sigma + (1/dy)^sigma )^(1/sigma) ]
  // 
  // where sp and lambda depend on the orderOfAccuracyInSpace:
  //   sp     = sosupPower[orderOfAccuracyInSpace], 
  //   lambda = sosupLambda[orderOfAccuracyInSpace], 
  const int maxOrderOfAccuracy=10;
  assert( orderOfAccuracyInSpace<=maxOrderOfAccuracy );
  real sosupPower2d[maxOrderOfAccuracy]  = { 2., 2., 2., 2., 2., 2., 2., 2., 2., 2. }; // note: some entries not used
  real sosupLambda2d[maxOrderOfAccuracy] = { 1., 1., 1., 1., 1., 1., 1., 1., 1., 1. }; // 
  real sosupPower3d[maxOrderOfAccuracy]  = { 2., 2., 2., 2., 2., 2., 2., 2., 2., 2. }; // 
  real sosupLambda3d[maxOrderOfAccuracy] = { 1., 1., 1., 1., 1., 1., 1., 1., 1., 1. }; // 
    
  // From Jeff Banks:
  //  2nd order: sigma=1.35,   b=.605
  //  4th order: sigma=2.175,  b=1.075
  //  6th order: sigma=1.6,    b=1.275
  //  6th order: sigma=1.5,    b=1.55 (Sept 2015)
  sosupPower2d[2]=1.35;  sosupLambda2d[2]=.605;   // 2nd order 2D
  sosupPower2d[4]=2.175; sosupLambda2d[4]=1.075;  // 4th order 2D, beta = 1
  //sosupPower2d[4]=1.6;   sosupLambda2d[4]=1.4;    // 4th order 2D, beta = 0.8
  //sosupPower2d[6]=1.6;   sosupLambda2d[6]=1.275;  // 6th order 2D
  sosupPower2d[6]=1.5;   sosupLambda2d[6]=1.55;  // 6th order 2D (Sept 2015)
    
  // *finish me for 3D:*
  sosupPower3d[2]=1.35;  sosupLambda3d[2]=.605;   // 2nd order 3D

  // sosupPower3d[4]=2.175; sosupLambda3d[4]=1.075;  // 4th order 3D, beta = 1
  // *wdh* July 1, 2016 -- reduced cfl seems to be needed (dieletric sphere G2)
  // should be able to run at cfl=.95 with this: 
  sosupPower3d[4]=2.175; sosupLambda3d[4]=1.075*.75;  // 4th order 3D, beta = 1 *wdh* reduce by .75 
             
  //sosupPower3d[4]=1.6;   sosupLambda3d[4]=1.4;    // 4th order 3D, beta = 0.8
  //sosupPower3d[6]=1.6;   sosupLambda3d[6]=1.275;  // 6th order 3D
  sosupPower3d[6]=1.5;   sosupLambda3d[6]=1.55;  // 6th order 3D (Sept 2015)


  // ============== FD + SOSUP Dissipation ============
  const int & useSosupDissipation = parameters.dbase.get<int>("useSosupDissipation");
  const real & sosupParameter = parameters.dbase.get<real>("sosupParameter");    // scaling of sosup dissipation

  // ** FINISH ME **
  

  real cMax=max(cGrid);
  if( numberOfMaterialRegions>1 && method==yee )
  { // Compute maximum c for variable eps and mu
    assert( numberOfComponentGrids==1 );
    cMax = sqrt(1./min(epsv*muv));
    // cMax = ParallelUtility::getMaxValue(cMax);
    printF("computeTimeStep: numberOfMaterialRegions=%i cMax=%9.3e\n",numberOfMaterialRegions,cMax);
  }
  
  real dtBA=0;
  if( method==bamx )
  {
    // If there is only one material, the materialRegionParameters vector may not have been created.
    // Do this here for now. Is this the right place?
    if( !dbase.has_key("materialRegionParameters") )
    {
      std::vector<DispersiveMaterialParameters> & dmpVector =
        dbase.put<std::vector<DispersiveMaterialParameters> >("materialRegionParameters");

      // Material "0" is the background material: 
      const int grid=0;
      const int domain = cg.domainNumber(grid);
      const DispersiveMaterialParameters & dmp0 = getDomainDispersiveMaterialParameters(domain);
      dmpVector.push_back(dmp0);
      // aString label="dmpVector[0]";
      // dmpVector[0].display(stdout,label);
    }


    // --- Fill in the material matrix for the embedded material regions ----
    std::vector<DispersiveMaterialParameters> & dmpVector = 
      dbase.get<std::vector<DispersiveMaterialParameters> >("materialRegionParameters");

    cMax=0.;
    for( int mr=0; mr<numberOfMaterialRegions; mr++ )
    {
      DispersiveMaterialParameters & dmp = dmpVector[mr];  
      real reLambda, imLambda;
      dmp.evalBianisotropicEigenValues( numberOfDimensions, reLambda, imLambda ); 
      // imLambda better be zero or the problem is ill-posed !
      c = reLambda;  // maximum wave speed
      cMax=max(cMax,c);

      printF("CgMx:computeTimeStep: BA material region %d : c=%9.3e (imLambda=%9.3e) cMax=%9.3e\n",mr,c,
             imLambda,cMax );

    }
    
  }
  
   
  // *wdh* Oct. 1, 2020 -- more accurate time step (e.g. nonSquare now matches square)
  // Nov 23, 2020 -- something is wrog with the new way for twoSquaresInterface left-grid, eps=2
  // Old: dt=2.8e-3
  // New: dt=3.8e-3 !!
  const bool useNewTimeStep=true; // false; // true;

  for( int grid=0; grid<numberOfComponentGrids; grid++ )
  {
    MappedGrid & mg = cg[grid];

    const int numberOfDimensions=mg.numberOfDimensions();
    real c = cGrid(grid);   
    // eps = epsGrid(grid);
    // mu = muGrid(grid);

    if( method==bamx )
    {
      c = cMax;

      // const int domain = cg.domainNumber(grid);
      // DispersiveMaterialParameters & dmp = getDomainDispersiveMaterialParameters(domain);
      // real reLambda, imLambda;
      // dmp.evalBianisotropicEigenValues( numberOfDimensions, reLambda, imLambda ); 
      // // imLambda better be zero or the problem is ill-posed !
      // c = reLambda;  // maximum wave speed
      // printF("CgMx:computeTimeStep: grid=%d - BA max wave speed c=%9.3e (reLambda=%9.3e, imLambda=%9.3e)\n",grid,c,
      //             reLambda,imLambda);

      // // OV_ABORT("stop here for now");
      
    }
    

    if( method==yee && numberOfMaterialRegions>1 )
      c=cMax;
        

    real dtg=REAL_MAX*.01;
    if( mg.getGridType()==MappedGrid::structuredGrid )
    {
      real dx[3];
      if( mg.isRectangular() )
      {
        mg.getDeltaX(dx);

        if( method==nfdtd || method==yee  )
        {
          if( numberOfDimensions==2 )
            dtg=cfl*1./( c*sqrt( 1./(dx[0]*dx[0])+1./(dx[1]*dx[1]) ) );  
          else
            dtg=cfl*1./( c*sqrt( 1./(dx[0]*dx[0])+1./(dx[1]*dx[1])+1./(dx[2]*dx[2]) ) ); 
        }
        else if( method==bamx )
        {
          // **Check me** New Oct 2019
          if( numberOfDimensions==2 )
            dtg=cfl*1./( c*( 1./dx[0] + 1./dx[1] ) );  
          else
            dtg=cfl*1./( c*( 1./dx[0] +1./dx[1] +1./dx[2] ) ); 
        }
        else if( method==sosup )
        {
          // SOSUP: dt depends on the order of accuracy: 
          if( numberOfDimensions==2 )
          {
            const real lambda = sosupLambda2d[orderOfAccuracyInSpace], sp = sosupPower2d[orderOfAccuracyInSpace];
            dtg = cfl*(lambda/c)/( pow( pow(1./dx[0],sp) + pow(1./dx[1],sp) , 1./sp ) );
          }
          else
          {
            const real lambda = sosupLambda3d[orderOfAccuracyInSpace], sp = sosupPower3d[orderOfAccuracyInSpace];
            dtg = cfl*(lambda/c)/( pow( pow(1./dx[0],sp) + pow(1./dx[1],sp) + pow(1./dx[2],sp) , 1./sp ) );
          }
                               
        }
        else
        {
          OV_ABORT("computeTimeSTep::ERROR: unknown method");
        }
        
        dxMinMax(grid,0)=numberOfDimensions==2 ? min(dx[0],dx[1]) : min(dx[0],dx[1],dx[2]);
        dxMinMax(grid,1)=numberOfDimensions==2 ? max(dx[0],dx[1]) : max(dx[0],dx[1],dx[2]);
        
        printf(" computeTimeStep: grid=%d dx=%8.2e dy=%8.2e c=%8.2e, dtg=%8.2e\n",grid,dx[0],dx[1],c,dtg);
      }
      else // curvilinear grids 
      {  // curvilinear grids 

        mg.update(MappedGrid::THEinverseVertexDerivative);
        const realArray & rx = mg.inverseVertexDerivative();
        const intArray & mask = mg.mask();
        
      
        Index I1,I2,I3;
        getIndex( mg.indexRange(),I1,I2,I3);

        // Grid spacings on unit square:
        real dr1 = mg.gridSpacing(axis1);
        real dr2 = mg.gridSpacing(axis2);
        real dr3 = mg.gridSpacing(axis3);

        // parallel version here --- also broadcast max error in forcing.bC *************************
        #ifdef USE_PPP
          realSerialArray rxLocal; getLocalArrayWithGhostBoundaries(rx,rxLocal);
          intSerialArray maskLocal; getLocalArrayWithGhostBoundaries(mask,maskLocal);
        #else
          const realSerialArray & rxLocal = rx;
          const intSerialArray & maskLocal = mask;
        #endif

        real *rxp = rxLocal.Array_Descriptor.Array_View_Pointer3;
        const int rxDim0=rxLocal.getRawDataSize(0);
        const int rxDim1=rxLocal.getRawDataSize(1);
        const int rxDim2=rxLocal.getRawDataSize(2);
        const int rxDim3=mg.numberOfDimensions();   // note
#undef RX
#define RX(i0,i1,i2,i3,i4) rxp[i0+rxDim0*(i1+rxDim1*(i2+rxDim2*(i3+rxDim3*(i4))))]

        const int *maskp = maskLocal.Array_Descriptor.Array_View_Pointer2;
        const int maskDim0=maskLocal.getRawDataSize(0);
        const int maskDim1=maskLocal.getRawDataSize(1);
        const int md1=maskDim0, md2=md1*maskDim1; 
#define MASK(i0,i1,i2) maskp[(i0)+(i1)*md1+(i2)*md2]


        int includeGhost=0;
        bool ok = ParallelUtility::getLocalArrayBounds(rx,rxLocal,I1,I2,I3,includeGhost);

        int i1,i2,i3;
        real a11Min=REAL_MAX*.001;
        real a11Max=-a11Min;
        //  **** this is a guess **** check this.
        const real alpha0=1.;

        dxMinMax(grid,0)=REAL_MAX*.01; 
        dxMinMax(grid,1)=0.;
        dtg = REAL_MAX*.01;

        real a11,a12,a22;
        if( ok )
        {
          if( numberOfDimensions==2 )
          {

            if( method!=sosup )
            {
              FOR_3D(i1,i2,i3,I1,I2,I3)
              {
              
                if( MASK(i1,i2,i3)>0 )
                {
                  a11 = ( RX(i1,i2,i3,0,0)*RX(i1,i2,i3,0,0) + RX(i1,i2,i3,0,1)*RX(i1,i2,i3,0,1) );
                  a12 = ( RX(i1,i2,i3,0,0)*RX(i1,i2,i3,1,0) + RX(i1,i2,i3,0,1)*RX(i1,i2,i3,1,1) )*2.;
                  a22 = ( RX(i1,i2,i3,1,0)*RX(i1,i2,i3,1,0) + RX(i1,i2,i3,1,1)*RX(i1,i2,i3,1,1) );

                  // we could save work by delaying the sqrt to after the loop
                  a11=1./sqrt( a11 *(1./(alpha0*dr1*dr1)) 
                               +abs(a12)*(.25/(alpha0*dr1*dr2))
                               +a22 *(1./(alpha0*dr2*dr2)) 
                    );

                  a11Min=min(a11Min,a11);
                  a11Max=max(a11Max,a11);
          
                }

              }
            }
            else
            {
              // sosup:

              // FIX dxMin dxMax !
              const real lambda = sosupLambda2d[orderOfAccuracyInSpace], sp = sosupPower2d[orderOfAccuracyInSpace];
              const real spBy2=sp*.5;
              FOR_3D(i1,i2,i3,I1,I2,I3)
              {
              
                if( MASK(i1,i2,i3)>0 )
                {
                  a11 = ( RX(i1,i2,i3,0,0)*RX(i1,i2,i3,0,0) + RX(i1,i2,i3,0,1)*RX(i1,i2,i3,0,1) );
                  a12 = ( RX(i1,i2,i3,0,0)*RX(i1,i2,i3,1,0) + RX(i1,i2,i3,0,1)*RX(i1,i2,i3,1,1) )*2.;
                  a22 = ( RX(i1,i2,i3,1,0)*RX(i1,i2,i3,1,0) + RX(i1,i2,i3,1,1)*RX(i1,i2,i3,1,1) );

                  // we could save work by delaying the outer pow to after the loop
                  a11=lambda/pow( pow(a11 *(1./(alpha0*dr1*dr1)),spBy2) +
                                  pow(abs(a12)*(.25/(alpha0*dr1*dr2)),spBy2) +
                                  pow(a22 *(1./(alpha0*dr2*dr2)),spBy2), 1./sp );

                  a11Min=min(a11Min,a11);
                  a11Max=max(a11Max,a11);
          
                }
              }
            }
            
            
          }
          else  //   ***** 3D ********
          { //   ***** 3D ********

#define rxDotRx(axis,dir) (RX(i1,i2,i3,axis,0)*RX(i1,i2,i3,dir,0) \
                         + RX(i1,i2,i3,axis,1)*RX(i1,i2,i3,dir,1) \
                         + RX(i1,i2,i3,axis,2)*RX(i1,i2,i3,dir,2))
      
            if( method!=sosup )
            {
              // There would be a factor of 4 for the worst case plus/minus wave but we also
              // divide by a factor of 4 for the 2nd-order time stepping.
              FOR_3D(i1,i2,i3,I1,I2,I3)
              {
                if( MASK(i1,i2,i3)>0 )
                {
                 // we could save work by delaying the sqrt to after the loop
                  a11=1./sqrt(   rxDotRx(0,0) *(1./(dr1*dr1)) 
                                 +rxDotRx(1,1) *(1./(dr2*dr2))
                                 +rxDotRx(2,2) *(1./(dr3*dr3))
                                 +abs(rxDotRx(1,0))*(.5/(dr2*dr1))  
                                 +abs(rxDotRx(2,0))*(.5/(dr3*dr1)) 
                                 +abs(rxDotRx(2,1))*(.5/(dr3*dr2)) );

                  // ** a11 =  pow(a11,-.5);
                
                  a11Min=min(a11Min,a11);
                  a11Max=max(a11Max,a11);

                }
              }
            }
            else
            { // sosup: 
              const real lambda = sosupLambda3d[orderOfAccuracyInSpace], sp = sosupPower3d[orderOfAccuracyInSpace];
              const real spBy2=sp*.5;
              FOR_3D(i1,i2,i3,I1,I2,I3)
              {
                if( MASK(i1,i2,i3)>0 )
                {
                  // we could save work by delaying the outer pow to after the loop
                  a11=lambda/pow( pow(rxDotRx(0,0) *(1./(dr1*dr1)),spBy2) +
                                  pow(rxDotRx(1,1) *(1./(dr2*dr2)),spBy2) +
                                  pow(rxDotRx(2,2) *(1./(dr3*dr3)),spBy2) +
                                  pow(abs(rxDotRx(1,0))*(.5/(dr2*dr1)),spBy2) + 
                                  pow(abs(rxDotRx(2,0))*(.5/(dr3*dr1)),spBy2) + 
                                  pow(abs(rxDotRx(2,1))*(.5/(dr3*dr2)),spBy2), 1./sp );

                  // ** a11 =  pow(a11,-.5);
                
                  a11Min=min(a11Min,a11);
                  a11Max=max(a11Max,a11);

                }
              }
            }
            
#undef rxDotRx
          } // end if 3D
          
          

        }
        // end if OK 
        
        if( useNewTimeStep )
        {
          // *wdh* Oct. 1, 2020 -- more accurate time step (e.g. nonSquare now matches square)
          real dsMin[3],dsAve[3],dsMax[3];
          // Note: there is communication here: 
          GridStatistics::getGridSpacing( mg, dsMin, dsAve, dsMax );

          assert( method!=sosup );

          // Nov 23, 2020 *fixed* formula
          dx[0]=dsMin[0]; dx[1]=dsMin[1]; dx[2]=dsMin[2];
          
          if( numberOfDimensions==2 )
            dtg=cfl*1./( c*sqrt( 1./(dx[0]*dx[0])+1./(dx[1]*dx[1]) ) );  
          else
            dtg=cfl*1./( c*sqrt( 1./(dx[0]*dx[0])+1./(dx[1]*dx[1])+1./(dx[2]*dx[2]) ) ); 

          dxMinMax(grid,0) = numberOfDimensions == 2 ? min(dsMin[0],dsMin[1]) : min(dsMin[0],dsMin[1],dsMin[2]);
          dxMinMax(grid,1) = numberOfDimensions == 2 ? max(dsMax[0],dsMax[1]) : max(dsMax[0],dsMax[1],dsMax[2]);
        }
        else
        {
          dxMinMax(grid,0)=a11Min;  // *wdh* Sept 30, 2020 -- use actual grid spacing below:
          dxMinMax(grid,1)=a11Max;
  
          dtg = (cfl/c) * dxMinMax(grid,0); 
        }
          



      }    // end curvilinear grid 

      // if( debug & 1 )
      // {
      //   printf("grid=%d: dtg=%e debug=%d(before)\n",grid,dtg,debug);
      //   fprintf(pDebugFile,"grid=%d: dtg=%e (before)\n",grid,dtg);
      //   ::display(dxMinMax,"dxMinMax",pDebugFile);
      //   fflush(pDebugFile);
      // }
      
      dtg             =ParallelUtility::getMinValue(dtg);  // compute min over all processors
      dxMinMax(grid,0)=ParallelUtility::getMinValue(dxMinMax(grid,0));
      dxMinMax(grid,1)=ParallelUtility::getMaxValue(dxMinMax(grid,1));
      
      // if( debug & 1 )
      // {
      //   fprintf(pDebugFile,"grid=%d: dtg=%e (after)\n",grid,dtg);
      //   fflush(pDebugFile);
      // }
      
      // Cartesian grids use: artificialDissipation
      // Curvilinear grids use: artificialDissipationCurvilinear
      const real artDiss = mg.isRectangular() ? artificialDissipation : artificialDissipationCurvilinear;
       
      if( useSosupDissipation!=0 && artDiss !=0. )
      {
        printF("--MX-- getTimeStep: ERROR: useSosupDissipaton but normal artificial dissipation is also on!\n");
        OV_ABORT("error");
      }
      

      if( artDiss>0. )
      {
        // Here is the correction for artificial dissipation
        //
        // The equation for dt looks like
        //   dt*dt *c*c*(  1/dx^2 + 1/dy^2 ) = 1 - beta*dt
        //
        real gamma = dtg*dtg;
        real beta;
        // const real adc = artDiss*SQR(cMax); // scale dissipation by c^2 *wdh* 041103
        const real adc = c*artDiss; // do this now *wdh* 090602

        beta = .5*adc*( numberOfDimensions*pow(2.,real(orderOfArtificialDissipation)) );
        real factor=2.;  // safety factor
        beta *=factor;

        dtg = sqrt( gamma + pow(beta*gamma*.5,2.) ) - beta*gamma*.5;

        if( debug & 4 )
          fprintf(pDebugFile," getTimeStep: Correct for art. dissipation: new dt=%9.3e (old = %9.3e, new/old=%4.2f) myid=%i\n",
                 dtg,sqrt(gamma),dtg/sqrt(gamma),myid);

        if( true )
          printF("***** getTimeStep: Correct for art. dissipation: new dt=%9.3e (old = %9.3e, new/old=%4.2f)\n",
                 dtg,sqrt(gamma),dtg/sqrt(gamma));
        
      }
      
      if( timeSteppingMethod==modifiedEquationTimeStepping )
      {
        if( true || orderOfAccuracyInTime==2 || orderOfAccuracyInTime==4 )
        {
          dtg*=1.; // Check this for 3D
        }
        else
        {
          printF("getTimeStep:ERROR: modifiedEquationTimeStepping -- orderOfAccuracyInTime=%i ??\n",
          orderOfAccuracyInTime);
          
          Overture::abort("getTimeStep:ERROR: modifiedEquationTimeStepping -- orderOfAccuracyInTime?? ");
        }
        
      }
      else if( timeSteppingMethod==rungeKutta )
      {
        // MOL RK stability bounds -- Oct 20, 2019 

        const int & orderOfRungeKutta = dbase.get<int>("orderOfRungeKutta");
        if( orderOfRungeKutta == 1 || orderOfRungeKutta==2 )
        {
          // formally unstable with no-dissipation 
        }
        else if( orderOfRungeKutta == 3 )
        {
          if( !useSosupDissipation && artDiss==0. )
          {
            dtg *= 1.7;   //  approximate RK3-SPP, maybe 1.7 on axis 
          }
          else
          {
             dtg *= 2.3; //  approximate RK3-SPP with dissipation -- **FIX ME** 
          }
          
          
          if( orderOfAccuracyInSpace==4 )
            dtg = dtg/1.4;   // approximate for four-order first derivative *check me*

        }
        else if( orderOfRungeKutta == 4 )
        {
          const real rk4ImBound = 2.8;  // RK stability bound on the imaginary axis 


          if( orderOfAccuracyInSpace==2 )
          {     
            dtg *= rk4ImBound; 
          }
          else if( orderOfAccuracyInSpace==4 )
          {
            const real firstDerivSymbolBound = 1.4; // approximate bound on the symbol of D0( I - h^2/6 D+D-) operator 
            dtg = rk4ImBound*dtg/firstDerivSymbolBound;  
          }
          else
          {
            OV_ABORT("finish me");
          }
          
          
        }
        else
        {
          OV_ABORT("ERROR: getTimeStep: unexpected orderOfRungeKutta");
        }
        

      }
      else

      {

        if( orderOfAccuracyInSpace==2 )
        {
        }
        else if( orderOfAccuracyInSpace==4 )
          dtg*=sqrt(3./4.);
        else if( orderOfAccuracyInSpace==6 )
          dtg*=sqrt(.6618);
        else if( orderOfAccuracyInSpace==8 )
          dtg*=sqrt(.6152);
        else
        {
          Overture::abort("getTimeStep:ERROR: modifiedEquationTimeStepping -- orderOfAccuracyInSpace?? ");
        }

        if( orderOfAccuracyInTime==4 )
        {
          dtg*=1.41/2.;
        }
        else if( orderOfAccuracyInTime==6 )
        {
          dtg*=.84/2.;
        }
        else if( orderOfAccuracyInTime==8 )
        {
          dtg*=.46/2.;
        }
        else if( orderOfAccuracyInTime==3 && method==dsi )
        {
          dtg*=(12./7.)/2.;   // ABS3
        }
        else if( orderOfAccuracyInTime!=2 )
        {
          Overture::abort();
        }
      }
      
      
      printF(" computeTimeStep: grid=%i c=%8.2e, dtg=%8.2e min-dx=%8.2e max-dx=%8.2e\n",grid,c,dtg,dxMinMax(grid,0),dxMinMax(grid,1));
      // printf(" computeTimeStep: grid=%i c=%8.2e, dtg=%8.2e min-dx=%8.2e max-dx=%8.2e myid=%i\n",
      //          grid,c,dtg,dxMinMax(grid,0),dxMinMax(grid,1),myid);

#undef RX
    
    }
    else
    {
      // unstructured grid.

//        UnstructuredMapping & map = (UnstructuredMapping &) mg.mapping().getMapping();
//        const int numberOfElements = map.getNumberOfElements();
//        const intArray & faces = map.getFaces();
//        const intArray & elementFace = map.getElementFaces();
//        const realArray & nodes = map.getNodes();

//        real dsMax=0.;
//        int e;
//        for( e=0; e<numberOfElements; e++ )
//        {
//      const int numFacesThisElement=map.getNumberOfFacesThisElement(e);      

//      // select two faces from the element:
//      int f1=elementFace(e,0);
//      int f2=elementFace(e,numFacesThisElement-1);
      
//      int n0=faces(f1,0), n1=faces(f1,1);
//      real ds1=SQR( nodes(n1,0)-nodes(n0,0) )+SQR( nodes(n1,1)-nodes(n0,1) );
//      n0=faces(f2,0), n1=faces(f2,1);
//      real ds2=SQR( nodes(n1,0)-nodes(n0,0) )+SQR( nodes(n1,1)-nodes(n0,1) );

//      dsMax = max( dsMax, 1./ds1+1./ds2 );
//        }
    
//        if( map.getNumberOfFacesThisElement(0)==3 )
//        {
//      dsMax*=4.;  // **** fudge for triangular grids -- fix this ---
//        }
    
#if 0
      mg.update(MappedGrid::THEminMaxEdgeLength);
      real dsMax =8./(mg.minimumEdgeLength(0)*mg.minimumEdgeLength(0));

      printF(" computeTimeStep: dsMax=%e c=%e\n",dsMax,c);

      dtg=cfl*1./( c*sqrt(dsMax) );
#else

      UnstructuredMapping::EntityTypeEnum cellType = UnstructuredMapping::EntityTypeEnum(cg.numberOfDimensions());
      UnstructuredMapping::EntityTypeEnum cellBdyType = UnstructuredMapping::EntityTypeEnum(cg.numberOfDimensions()-1);
      UnstructuredMapping::EntityTypeEnum faceType = UnstructuredMapping::Face;
      UnstructuredMapping::EntityTypeEnum edgeType = UnstructuredMapping::Edge;

      UnstructuredMapping &umap = (UnstructuredMapping &)mg.mapping().getMapping();

      int rDim = umap.getRangeDimension();
      int dDim = umap.getDomainDimension();

      bool vCent = mg.isAllVertexCentered();
      realArray &cFArea = vCent ? mg.centerArea() : mg.faceArea();
      const realArray &cFNorm = vCent ? mg.centerNormal() : mg.faceNormal();
      const realArray &cEArea = vCent ? mg.faceArea() : mg.centerArea();
      const realArray &cENorm = vCent ? mg.faceNormal() : mg.centerNormal();

      realArray cellBdyCenters, cellCenters;
      getCenters( mg, cellBdyType, cellBdyCenters);
      getCenters( mg, cellType, cellCenters);

      const realArray &nodes  = umap.getNodes();
      const intArray &edges = umap.getEntities(edgeType);

      UnstructuredMappingIterator iter,iter_end;
      UnstructuredMappingAdjacencyIterator aiter, aiter_end;

      iter_end = umap.end(UnstructuredMapping::Face);
      real minLen = REAL_MAX;
      real minFEL = REAL_MAX;
      for ( iter = umap.begin(UnstructuredMapping::Face); iter!=iter_end; iter++ )
      {
        int f=*iter;
        real area = cFArea(f,0,0);
        real lsum = 0.;
        aiter_end = umap.adjacency_end(iter, UnstructuredMapping::Edge);
        for ( aiter = umap.adjacency_begin(iter, UnstructuredMapping::Edge); aiter!=aiter_end; aiter++ )
        {
          int e = *aiter;
          real edgeL = 0;
          int v1 = edges(e,0);
          int v2 = edges(e,1);
          for ( int a=0; a<rDim; a++ )
            edgeL += ( nodes(v2,a)-nodes(v1,a) )*( nodes(v2,a)-nodes(v1,a) );
          //          lsum += sqrt(edgeL);
          lsum = max(lsum,edgeL);
        }

        lsum = sqrt(lsum);

        minLen = min(minLen, area/lsum);
        minFEL = min(lsum,minFEL);
      }

      iter_end = umap.end(cellBdyType);
      real minCCL = REAL_MAX;
      for ( iter=umap.begin(cellBdyType); iter!=iter_end; iter++ )
      {
        int e = *iter;
        aiter_end = umap.adjacency_end(iter,cellType);
        real edgeL = 0;
        for ( aiter=umap.adjacency_begin(iter,cellType); aiter!=aiter_end; aiter++ )
        {
          real L = 0;
          int c = *aiter;
          for ( int a=0; a<rDim; a++ )
            L += (cellCenters(c,a)-cellBdyCenters(e,a))*(cellCenters(c,a)-cellBdyCenters(e,a));
          edgeL += sqrt(L);
        }

        minCCL = min(minCCL, edgeL);
        if ( rDim==2 )
          minLen = min(edgeL,minLen);
        else
        {
          aiter_end = umap.adjacency_end(iter,edgeType);
          for ( aiter=umap.adjacency_begin(iter,edgeType); aiter!=aiter_end; aiter++ )
          {
            real area=0;
            for ( int a=0; a<rDim; a++ )
              area+=edgeAreaNormals(*aiter,0,0,a)*edgeAreaNormals(*aiter,0,0,a);
            minLen = min(minLen, sqrt(area)/(2*edgeL));
          }
        }
          
      }

      minLen = min(minLen,minCCL);
      // *wdh* 100827 -- commented these next lines out since MappedGrid no longer has minMaxEdgeLength. Is this ok?
      // mg.update(MappedGrid::THEminMaxEdgeLength);
      // cout<<"minFEL = "<<minFEL<<", minCCL = "<<minCCL<<", minLen = "<<minLen<<", minEdge = "<<mg.minimumEdgeLength(0)<<endl;
      // minLen = min(minLen,mg.minimumEdgeLength(0));
      dtg=cfl*minLen/c;

#endif
      
    }
    
    
    deltaT=min(deltaT,dtg);
    
  } // end for grid

  // if( debug & 1 )
  // {
  //   fprintf(pDebugFile,"computeDT: deltaT=%e (before)\n",deltaT);
  //   fflush(pDebugFile);
  // }
  

  deltaT = ParallelUtility::getMinValue(deltaT);  // min value over all processors
  deltaT = min(dtMax,deltaT);
  
  // if( debug & 1 )
  // {
  //   fprintf(pDebugFile,"computeDT: deltaT=%e (after)\n",deltaT);
  //   fflush(pDebugFile);
  // }

  printF("==== computeTimeStep: deltaT=%8.2e\n",deltaT);
  if( debug & 4 )
    fprintf(pDebugFile,"==== computeTimeStep: deltaT=%8.2e, myid=%i\n",deltaT,myid);

  timing(timeForComputingDeltaT)+=getCPU()-time0;
  return 0;
}



// =============================================================================================
/// \brief Determine the time step for the BA Maxwell equations
// =============================================================================================
int Maxwell::
computeTimeStepBA()
{
  real time0=getCPU();

  assert( method==bamx );
  
  assert( cgp!=NULL );
  CompositeGrid & cg= *cgp;
  const int numberOfComponentGrids = cg.numberOfComponentGrids();
  const int numberOfDimensions = cg.numberOfDimensions();

  const real & dtMax = dbase.get<real>("dtMax");

  deltaT=REAL_MAX*.01;
  
  
  real dtBA=0;
  // If there is only one material, the materialRegionParameters vector may not have been created.
  // Do this here for now. Is this the right place?
  if( !dbase.has_key("materialRegionParameters") )
  {
    std::vector<DispersiveMaterialParameters> & dmpVector =
      dbase.put<std::vector<DispersiveMaterialParameters> >("materialRegionParameters");

    // Material "0" is the background material: 
    const int grid=0;
    const int domain = cg.domainNumber(grid);
    const DispersiveMaterialParameters & dmp0 = getDomainDispersiveMaterialParameters(domain);
    dmpVector.push_back(dmp0);
    // aString label="dmpVector[0]";
    // dmpVector[0].display(stdout,label);
  }

  real imStabilityBound=-1;
  if( timeSteppingMethod==modifiedEquationTimeStepping )
  {
    imStabilityBound=1;
  }
  else if( timeSteppingMethod==rungeKutta )
  {
    const int & orderOfRungeKutta = dbase.get<int>("orderOfRungeKutta");
    if( orderOfRungeKutta == 1 || orderOfRungeKutta==2 )
    {
      // formally unstable with no-dissipation 
     imStabilityBound=1;
    }
    else if( orderOfRungeKutta == 3 )
    {
      if( artificialDissipation==0. )
      {
        imStabilityBound = 1.7;   //  approximate RK3-SPP, maybe 1.7 on axis 
      }
      else
      {
        imStabilityBound = 2.; //  approximate RK3-SPP with dissipation -- **check me** 
      }
    }
    else if( orderOfRungeKutta == 4 )
    {
      imStabilityBound = 2.8;  // RK stability bound on the imaginary axis 
    }
    

  }
  else
  {
    OV_ABORT("finish me");
  }
  assert( imStabilityBound>0. );

  std::vector<DispersiveMaterialParameters> & dmpVector = 
    dbase.get<std::vector<DispersiveMaterialParameters> >("materialRegionParameters");

  for( int grid=0; grid<numberOfComponentGrids; grid++ )
  {
    MappedGrid & mg = cg[grid];
    for( int mr=0; mr<numberOfMaterialRegions; mr++ )
    {
      DispersiveMaterialParameters & dmp = dmpVector[mr];  

      real reLambda, imLambda;
      dmp.evalBianisotropicTimeSteppingEigenvalue( mg, orderOfAccuracyInSpace, reLambda, imLambda );
      
      dtBA = cfl*imStabilityBound/imLambda;
      printF("CgMx: computeTimeStepBA: grid=%d, region=%d: estimated dt= cfl*imStabilityBound/imLambda = %9.3e\n",
             grid,mr,dtBA);

      deltaT=min(deltaT,dtBA);

    }
  }
  

  deltaT = ParallelUtility::getMinValue(deltaT);  // min value over all processors
  deltaT = min(dtMax,deltaT);
  
  printF("==== computeTimeStepBA: deltaT=%8.2e =========\n\n",deltaT);
  if( debug & 4 )
    fprintf(pDebugFile,"==== computeTimeStep: deltaT=%8.2e, myid=%i\n",deltaT,myid);

  timing(timeForComputingDeltaT)+=getCPU()-time0;
  return 0;
}


