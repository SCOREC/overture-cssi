// This file automatically generated from advanceIncompressible.bC with bpp.
//
//     INCOMPRESSIBLE ELASTICITY
// 


#include "Cgsm.h"
#include "display.h"
#include "CompositeGridOperators.h"
#include "ParallelUtility.h"
// #include "ParallelOverlappingGridInterpolator.h"
#include "SmParameters.h"
// #include "GridFunctionFilter.h"
#include "GridMaterialProperties.h"

// ===================================================================================
//   This macro extracts the boundary data arrays
// ===================================================================================

// ===============================================================================================
// This macro determines the pointers to the variable material properties that are
// used when calling fortran routines.
// ===============================================================================================

#define advIsm EXTERN_C_NAME(advism)
extern "C"
{
    void advIsm(const int&nd,
                            const int&n1a,const int&n1b,const int&n2a,const int&n2b,const int&n3a,const int&n3b,
                            const int&nd1a,const int&nd1b,const int&nd2a,const int&nd2b,const int&nd3a,const int&nd3b,
                            const int&nd4a,const int&nd4b,
                            const int&mask,const real&rx, const real&xy, 
                            const real&um, const real&u, real&un, const real&f,
                            const real&vt, const real&vn,const real&v,
                            const real&v1, const real&v2, const real&v3, const real&v4,
                            const real&vt1,const real&vt2,const real&vt3,const real&vt4,
                            const int&gridIndexRange, const int&boundaryCondition, const int&ipar, const real&rpar, int&ierr );

}


// fourth order dissipation 2D: ***** NOTE: this is minus of the 4th difference:  -(D+D-)^2 *********
#define FD4_2D(u,i1,i2,i3,c) (    -( u(i1-2,i2,i3,c)+u(i1+2,i2,i3,c)+u(i1,i2-2,i3,c)+u(i1,i2+2,i3,c) )   +4.*( u(i1-1,i2,i3,c)+u(i1+1,i2,i3,c)+u(i1,i2-1,i3,c)+u(i1,i2+1,i3,c) ) -12.*u(i1,i2,i3,c) )

// fourth order dissipation 3D:
#define FD4_3D(u,i1,i2,i3,c) (    -( u(i1-2,i2,i3,c)+u(i1+2,i2,i3,c)+u(i1,i2-2,i3,c)+u(i1,i2+2,i3,c)+u(i1,i2,i3-2,c)+u(i1,i2,i3+2,c) )   +4.*( u(i1-1,i2,i3,c)+u(i1+1,i2,i3,c)+u(i1,i2-1,i3,c)+u(i1,i2+1,i3,c)+u(i1,i2,i3-1,c)+u(i1,i2,i3+1,c) ) -18.*u(i1,i2,i3,c) )

#define FOR_3(i1,i2,i3,I1,I2,I3) for( i3=I3.getBase(); i3<=I3.getBound(); i3++ )  for( i2=I2.getBase(); i2<=I2.getBound(); i2++ )  for( i1=I1.getBase(); i1<=I1.getBound(); i1++ )  

#define FN(m) fn[m+numberOfFunctions*(grid)]




// =======================================================================================
// Macro: Advance the solution to the next step
//          Call the optimized advance routine
//  option: 0 = update solution
//          1 = add upwind dissipation
// =======================================================================================

// =======================================================================================
// Macro: Interpolate and apply boundary conditions
// =======================================================================================

// ===============================================================================================================
//  Get any forcing 
// ===============================================================================================================





// =============================================================================
/// \brief Advance the solution for
///       Incompressible Linear Elastic Equations
///
///     >>>>>>>>   MODIFIED EQUATION METHOD : Taylor Time-Stepping  <<<<<<<<
///
// \param ut : if not NULL, compute u.tt and return in this grid-function. Otherwise
//             return the solution gf[next] at t=t+dt.
// =============================================================================
void Cgsm::
advanceIncompressible( int current, real t, real dt,
                                              RealCompositeGridFunction *ut /* = NULL */, 
                                              real tForce /* = 0. */ )
{   

    int & globalStepNumber = parameters.dbase.get<int >("globalStepNumber");
    if( t<2*dt || debug() & 4 )
        printF("\n @@@@@@@ Cgsm::advanceIncompressible *NEW* Modified Equation Method called at globalStepNumber=%d, t=%9.3e\n",globalStepNumber,t);

    checkArrays("advanceIncompressible:start");


    FILE *& debugFile  =parameters.dbase.get<FILE* >("debugFile");
    FILE *& logFile    =parameters.dbase.get<FILE* >("logFile");
    FILE *& pDebugFile =parameters.dbase.get<FILE* >("pDebugFile");
    
    const int numberOfDimensions           = cg.numberOfDimensions();
    const int numberOfComponentGrids       = cg.numberOfComponentGrids();
    const int & numberOfComponents         = parameters.dbase.get<int >("numberOfComponents");
    const int & uc                         = parameters.dbase.get<int >("uc");
    const int & vc                         = parameters.dbase.get<int >("vc");
    const int & wc                         = parameters.dbase.get<int >("wc");
    const int & rc                         = parameters.dbase.get<int >("rc");
    const int & tc                         = parameters.dbase.get<int >("tc");
    const int & pc                         = parameters.dbase.get<int >("pc");
    const int & orderOfAccuracyInSpace     = parameters.dbase.get<int>("orderOfAccuracy");
    const int & orderOfTimeAccuracy        = parameters.dbase.get<int>("orderOfTimeAccuracy");
    const int & numberOfCorrections        = parameters.dbase.get<int>("numberOfCorrections"); 
    const int & skipLastPressureSolve      = parameters.dbase.get<int>("skipLastPressureSolve");  // For ILE predictor-corrector schemes  
    const int & projectLinearModeForILE    = parameters.dbase.get<int>("projectLinearMode"); 
    const int & projectLinearModeFrequency = parameters.dbase.get<int>("projectLinearModeFrequency"); 

    SmParameters::TimeSteppingMethodSm & timeSteppingMethodSm = 
                                                                      parameters.dbase.get<SmParameters::TimeSteppingMethodSm>("timeSteppingMethodSm");
    RealArray & timing = parameters.dbase.get<RealArray >("timing");

    Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
    Range C=numberOfComponents;
    const int prev = (current-1+numberOfTimeLevels) % numberOfTimeLevels;
    const int next = (current+1)                    % numberOfTimeLevels;

    real & rho=parameters.dbase.get<real>("rho");
    real & mu = parameters.dbase.get<real>("mu");
    real & lambda = parameters.dbase.get<real>("lambda");
    RealArray & muGrid = parameters.dbase.get<RealArray>("muGrid");
    RealArray & lambdaGrid = parameters.dbase.get<RealArray>("lambdaGrid");
    bool & gridHasMaterialInterfaces = parameters.dbase.get<bool>("gridHasMaterialInterfaces");
    int & debug = parameters.dbase.get<int >("debug");

    const int & upwindSOS = parameters.dbase.get<int>("upwindSOS"); 

    const real cMax=max(lambdaGrid+muGrid)/rho;

    const int computeUt = ut != NULL;  // compute u.tt


  // const bool extrapInterpolationNeighbours = upwindSOS;
  // if( extrapInterpolationNeighbours )
  // {
  //   // -- Extrapolate interpolation neighbours for the artificial dissipation ---
  //   // This may also be done in advanceSOS
  //   if( debug & 4 )
  //     printF("advanceIncompressible: Extrapolate interpolation neighbours for upwind dissipation, t=%9.3e\n",t);

  //   extrapolateInterpolationNeighbours( gf[current], C );
  // }
    
    const int useWhereMask = numberOfComponentGrids>1;  

    sizeOfLocalArraysForAdvance=0.;


    if( timeSteppingMethodSm != SmParameters::modifiedEquationTimeStepping )
    {
        printF("advanceIncompressible: not implemented for timeSteppingMethodSm=%d\n",(int)timeSteppingMethodSm);
        printF(" advanceAdamsPredictorCorrector2=%d\n",(int)SmParameters::adamsPredictorCorrector2);
        printF(" advanceAdamsPredictorCorrector4=%d\n",(int)SmParameters::adamsPredictorCorrector4);
        OV_ABORT("error"); 
    }


    for( int correction=0; correction<numberOfCorrections; correction++ )
    {
        if( numberOfCorrections>1 && t<=2.*dt )
        {
            printF(" >>>>>>> advanceIncompressible: t=%9.3e : correction=%d <<<<<<<\n",t,correction);
        }

    // ----------------- ADVANCE EQUATIONS -----------------
    // --------------- START LOOP OVER GRIDS -------------------
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            real time0=getCPU();

            MappedGrid & mg = cg[grid];
            MappedGridOperators & mgop = (*cgop)[grid];

            getIndex(mg.gridIndexRange(),I1,I2,I3);
        
      // getBoundsForPML( mg,Iv );

            realMappedGridFunction & fieldPrev    =gf[prev].u[grid];
            realMappedGridFunction & fieldCurrent =gf[current].u[grid];
            realMappedGridFunction & fieldNext    =gf[next].u[grid];


            realArray & um = fieldPrev;
            realArray & u  = fieldCurrent;
            realArray & un = fieldNext;

            if( debug & 4 )
            {
                fPrintF(debugFile," **** start of advance, t=%8.2e\n",t);
                fprintf(pDebugFile," **** start of advance, t=%8.2e\n",t);
                
                if( debug & 8 )
                {
                    display(um,sPrintF("um start of advance, t=%8.2e",t),debugFile,"%8.2e ");
                    display(u,sPrintF("u start of advance, t=%8.2e",t),debugFile,"%8.2e ");
                }
            }

            lambda = lambdaGrid(grid);
            mu = muGrid(grid);
            c1=(mu+lambda)/rho, c2= mu/rho;
            
            if( numberOfStepsTaken<1 ) 
                printF(" advanceIncompressible:INFO mu=%8.2e for grid=%i (%s) \n",mu,grid,(const char*)cg[grid].getName());
            
      // const real dtsq=dt*dt; 
      // const real adc=artificialDissipation*SQR(cMax); // scale dissipation by c^2 *wdh* 041103


            bool useOpt=true; // true;
            const bool isRectangular=mg.isRectangular();
            if( !isRectangular )
            {
         // *** PUT THIS IN SETUP SOMEWHERE ***** 
                real timea=getCPU();
                mg.update( MappedGrid::THEinverseVertexDerivative | MappedGrid::THEinverseCenterDerivative );
                timea=getCPU()-timea;
                timing(parameters.dbase.get<int>("timeForInitialize"))+=timea;

                time0+=timea;  // do not include with time for curvilinear
            }

      //  ::display(lap,"lap","%8.1e ");

      // --- get any forcing ----
                realArray f;
                const bool addForcing = forcingOption!=noForcing && 
                                                                forcingOption!=planeWaveBoundaryForcing && 
                                                                forcingOption!=twilightZoneForcing;
                if( addForcing )
                {
                    Index D1,D2,D3;
                    getIndex(mg.dimension(),D1,D2,D3);
                    f.partition(mg.getPartition());
                    f.redim(D1,D2,D3,C);  // could use some other array for work space ??
                    int option=1;  
                    if( addForcing )
                    {
            //kkc getForcing is called from advance but is also called from assignIC and getErrors.
            //    we have to add the timing in external to getForcing to avoid double counting the time
            //    in assignIC and getErrors
                        real timef = getCPU();
                        getForcing( next, grid,f,t,dt,option );
                        timing(parameters.dbase.get<int>("timeForForcing")) += getCPU()-timef;
                    }
                }
                OV_GET_SERIAL_ARRAY(real,f ,fLocal);

            OV_GET_SERIAL_ARRAY(real,um,umLocal);
            OV_GET_SERIAL_ARRAY(real,u ,uLocal);
            OV_GET_SERIAL_ARRAY(real,un,unLocal);
      

      // ------ advance the interior points -------
            {
                real timeAdv=getCPU();
                int gridType = isRectangular? 0 : 1;
                const real adc=artificialDissipation; // do not scale *wdh* 090216
                real *umptr=umLocal.getDataPointer();
                real *uptr =uLocal.getDataPointer();
                real *unptr=unLocal.getDataPointer();
                real *vtptr = uptr; 
                real *vnptr = uptr; 
                real *vptr  = uptr; 
                real *v1ptr = uptr; 
                real *v2ptr = uptr; 
                real *v3ptr = uptr; 
                real *v4ptr = uptr; 
                real *vt1ptr = uptr; 
                real *vt2ptr = uptr; 
                real *vt3ptr = uptr; 
                real *vt4ptr = uptr; 
                    real *fptr   = addForcing ? fLocal.getDataPointer() : uptr;
        // assert( !useVariableDissipation || variableDissipation!=NULL );
        // real *pVarDis = useVariableDissipation ? varDis.getDataPointer() : uptr;
                const intArray & mask = mg.mask();
                OV_GET_SERIAL_ARRAY_CONST(int,mask,maskLocal)
        // #ifdef USE_PPP
        //     intSerialArray maskLocal;  getLocalArrayWithGhostBoundaries(mask,maskLocal);
        // #else
        //     const intSerialArray & maskLocal = mask; 
        // #endif
                real *rxptr;
                if( isRectangular )
                {
                    rxptr=uptr;
                }
                else
                {
                #ifdef USE_PPP
                            realSerialArray rxLocal; getLocalArrayWithGhostBoundaries(mg.inverseVertexDerivative(),rxLocal);
                #else
                            const realSerialArray & rxLocal=mg.inverseVertexDerivative();
                #endif
                    rxptr = rxLocal.getDataPointer();
                }
                const bool centerNeeded=forcingOption==twilightZoneForcing || (forcingOption==planeWaveBoundaryForcing); // **************** fix this 
                #ifdef USE_PPP
                    realSerialArray xy;
                    if( centerNeeded ) getLocalArrayWithGhostBoundaries(mg.center(),xy);
                #else
                    const realSerialArray & xy = centerNeeded ? mg.center() : umLocal;
                #endif
                real *xyptr=xy.getDataPointer();
                int maskNull[1];
                int *maskptr = maskLocal.getDataPointer();
        // realSerialArray *dis = NULL;
        // real *pdis=uptr;
        // Macro to extract the pointers to the variable material property arrays
         // --- Variable material properies ---
                  GridMaterialProperties::MaterialFormatEnum materialFormat = GridMaterialProperties::constantMaterialProperties;
                  int ndMatProp=1;  // for piecewise constant materials, this is the leading dimension of the matVal array
                  int *matIndexPtr=maskptr;  // if not used, point to mask
                  real*matValPtr=uptr;       // if not used, point to u
                  if( parameters.dbase.get<int>("variableMaterialPropertiesOption")!=0 )
                  {
           // Material properties do vary 
                      std::vector<GridMaterialProperties> & materialProperties = 
                                parameters.dbase.get<std::vector<GridMaterialProperties> >("materialProperties");
                      GridMaterialProperties & matProp = materialProperties[grid];
                      materialFormat = matProp.getMaterialFormat();
                      if( materialFormat==GridMaterialProperties::piecewiseConstantMaterialProperties )
                      {
                                IntegerArray & matIndex = matProp.getMaterialIndexArray();
                          matIndexPtr = matIndex.getDataPointer();
                      }
                      RealArray & matVal = matProp.getMaterialValuesArray();
                      matValPtr = matVal.getDataPointer();
                      ndMatProp = matVal.getLength(0);  
           // ::display(matVal,"matVal");
                  }
                bool twilightZone = forcingOption==twilightZoneForcing; 
                int combineDissipationWithAdvance=0;  // not used anymore
                real tForce = t;   // evaluate forcing at this time 
                    int option=0; //  0=update solution, 1=add upwind dissipation
                int ipar[]={option,
                                        gridType,
                                        orderOfAccuracyInSpace,
                                        orderOfTimeAccuracy,
                                        (int)addForcing,
                                        (int)twilightZone,
                                        uc,
                                        vc,
                                        wc,
                                        useWhereMask,
                                        (int)timeSteppingMethodSm,
                                        (int)useVariableDissipation,
                                        (int)useConservative,           // 12 
                                        combineDissipationWithAdvance,
                                        debug,
                                        computeUt,
                                        materialFormat,                 // 16 
                                        myid,
                                        pc,                             // 18 
                                        upwindSOS,
                                        correction
                                        };    
                real dx[3]={1.,1.,1.};
                if( isRectangular )
                    mg.getDeltaX(dx);
                real rpar[30];
                rpar[ 0]=dt;
                rpar[ 1]=dx[0];
                rpar[ 2]=dx[1];
                rpar[ 3]=dx[2];
                rpar[ 4]=adc;
                rpar[ 5]=mg.gridSpacing(0);
                rpar[ 6]=mg.gridSpacing(1);
                rpar[ 7]=mg.gridSpacing(2);
                rpar[ 8]=rho;
                rpar[ 9]=mu;
                rpar[10]=kx; // for plane wave scattering
                rpar[11]=ky;
                rpar[12]=kz;
                rpar[13]=(real &)parameters.dbase.get<OGFunction* >("exactSolution");  // twilight zone pointer, ep
                rpar[14]=tForce;
                rpar[15]=parameters.dbase.get<real>("dtOld");  // dt used on the previous step
                rpar[16]=rho;
                rpar[17]=mu;
                rpar[18]=lambda;
                rpar[20]=0.;  // return cpu for dissipation
        // printF("** AdvanceIncompressible: gridType=%i, isRectangular=%i\n",gridType,(int)isRectangular);
                bool ok = ParallelUtility::getLocalArrayBounds(u,uLocal,I1,I2,I3);
                int ierr=0;
                if( ok )
                {
          // if( adc>0. && !combineDissipationWithAdvance )
          // {
          //   // create a temp array to hold the artificial dissipation
          //   dis = new realSerialArray(uLocal.dimension(0),uLocal.dimension(1),uLocal.dimension(2),uLocal.dimension(3));
          //   pdis = dis->getDataPointer();
          //   assert( pdis!=NULL );
          //   sizeOfLocalArraysForAdvance=max(sizeOfLocalArraysForAdvance,(double)(dis->elementCount()*sizeof(real)));
          // }
          // real timeAdv=getCPU();
          // if( combineDissipationWithAdvance )
          // {
          //   assert( umptr!=unptr );
          // }
                    advIsm(mg.numberOfDimensions(),
                                  I1.getBase(),I1.getBound(),I2.getBase(),I2.getBound(),I3.getBase(),I3.getBound(),
                                  uLocal.getBase(0),uLocal.getBound(0),
                                  uLocal.getBase(1),uLocal.getBound(1),
                                  uLocal.getBase(2),uLocal.getBound(2),
                                  uLocal.getBase(3),uLocal.getBound(3),
                                  *maskptr,*rxptr,*xyptr,  
                                  *umptr,*uptr,*unptr, *fptr, 
                                  *vtptr, *vnptr, *vptr,
                                  *v1ptr, *v2ptr, *v3ptr, *v4ptr,
                                  *vt1ptr,*vt2ptr,*vt3ptr,*vt4ptr,
                                  mg.gridIndexRange(0,0), mg.boundaryCondition(0,0),  ipar[0], rpar[0], ierr );
                }
                timeAdv=getCPU()-timeAdv;
                    timing(parameters.dbase.get<int>("timeForAdvOpt"))+=timeAdv;
                if( debug & 8 )
                {
                    display(unLocal,sPrintF("unLocal after advIncompressibleSolidMechanics, processor=%i before BC's t=%8.2e",
                                                                    Communication_Manager::My_Process_Number,t),pDebugFile,"%8.2e ");
                    display(un,sPrintF("un after advIncompressibleSolidMechanics, before BC's t=%8.2e",t),debugFile,"%8.2e ");
                }
                if( isRectangular )   
                    timing(parameters.dbase.get<int>("timeForAdvanceRectangularGrids"))+=getCPU()-time0;
                else
                    timing(parameters.dbase.get<int>("timeForAdvanceCurvilinearGrids"))+=getCPU()-time0;
                if( debug & 8 )
                {
                    display(u,sPrintF("u after advOpt and updateGhost, t=%8.2e",t),debugFile,"%8.2e ");
                }
                if( debug & 16 )
                {
                    if( addForcing ) 
                        ::display(f,sPrintF("  *** advIncompressibleSolidMechanics: Here is the forcing f grid=%i t=%9.3e ********\n",grid,t),
                                            pDebugFile,"%7.4f ");
                } 
            }

        
        } // end grid
        
        if( debug & 4 )
        {
            getErrors( next,t+dt,dt,sPrintF("\n *** advanceIncompressible: Errors after advance, before BC, t=%9.3e ******\n",t+dt) );
        }

        #undef FN

        
        gf[next].t=t+dt;

        if( correction==0 && orderOfAccuracyInSpace > 2  )
        {
      // ---- Extrapolate pressure in time ----
            if( t<= 2*dt )
                printF("Extrapolate p in time for BCs...\n"); 
            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
            {
                MappedGrid & mg = cg[grid];
                OV_GET_SERIAL_ARRAY(real,gf[prev].u[grid]    ,umLocal);
                OV_GET_SERIAL_ARRAY(real,gf[current].u[grid] , uLocal);
                OV_GET_SERIAL_ARRAY(real,gf[next].u[grid]    ,unLocal);
                getIndex(mg.dimension(),I1,I2,I3);

        // *** NOTE *** THIS ONLY NEEDS TO BE DONE NEAR CERTAIN BOUNDARIES ****************************
                if( orderOfTimeAccuracy==2 )
                {
          // unLocal(I1,I2,I3,pc) = 2.*uLocal(I1,I2,I3,pc) - umLocal(I1,I2,I3,pc);
          // assume unLocal holds p(t-2*dt)
                    unLocal(I1,I2,I3,pc) = 3.*uLocal(I1,I2,I3,pc) - 3.*umLocal(I1,I2,I3,pc) + unLocal(I1,I2,I3,pc);
                }
                else
                {
                    OV_ABORT("finish me -- extrap p in time"); 
                }

        // OV_ABORT("finish me");     
            }
        }


        bool applyBC=true; // ** FIX ME **
            if( applyBC )
            {
                if( true ||   // *wdh* 091205 -- interpolate will call periodicUpdate and updateGhost 
                        cg.numberOfComponentGrids()>1 )
                {
                    real timei=getCPU();
                    if( debug & 4 )
                        gf[next].u.display(sPrintF("Cgsm::advanceIncompressible: gf[next].u before interpolate, t=%8.2e",t+dt),debugFile,"%8.2e ");
          // --- Note: interpolate performs a periodicUpdate and updateGhostBoundaries even if there is only one grid
                    gf[next].u.interpolate();
                    if( debug & 4 )
                        gf[next].u.display(sPrintF("Cgsm::advanceIncompressible: gf[next].u after interpolate, t=%8.2e",t+dt),debugFile,"%8.2e ");
                    if( debug & 4 )
                    {
                        getErrors( next,t+dt,dt,sPrintF("\n ************** advance Errors after interpolate t=%9.3e ******\n",t+dt));
                    }
                    timing(parameters.dbase.get<int>("timeForInterpolate"))+=getCPU()-timei;
                }
        // ============= Boundary Conditions =============
                int option=0; // not used.
                applyBoundaryConditions( option, dt, next,current ); // apply BC to "next" (current=previous time step)
                if( debug & 8 )  // & 64
                {
                    gf[next].u.display(sPrintF("Cgsm::advanceIncompressible: gf[next].u after applyBC, t=%8.2e",gf[next].t),debugFile,"%8.2e ");
                }
                if( debug & 4 )
                  {
                      getErrors( next,t+dt,dt,sPrintF("\n ************** advanceIncompressible Errors after applyBC t=%9.3e ******\n",t+dt));
                  }    
        // ---- assign values at material interfaces ------
        // *** this does nothing currently ***
                assignInterfaceBoundaryConditions( next, t, dt );  // is this the right place to do this?
            }





    // ----------------- ADD UPWIND DISSIPATION -----------------
        if( upwindSOS && globalStepNumber>=1 )
        {  
      // --------------- START LOOP OVER GRIDS -------------------
            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
            {
                real time0=getCPU();

                MappedGrid & mg = cg[grid];
                MappedGridOperators & mgop = (*cgop)[grid];

                getIndex(mg.gridIndexRange(),I1,I2,I3);
            
        // getBoundsForPML( mg,Iv );

                realMappedGridFunction & fieldPrev    =gf[prev].u[grid];
                realMappedGridFunction & fieldCurrent =gf[current].u[grid];
                realMappedGridFunction & fieldNext    =gf[next].u[grid];


                realArray & um = fieldPrev;
                realArray & u  = fieldCurrent;
                realArray & un = fieldNext;

        // if( debug & 4 )
        // {
        //   fPrintF(debugFile," **** start of advance, t=%8.2e\n",t);
        //   fprintf(pDebugFile," **** start of advance, t=%8.2e\n",t);
                    
        //   if( debug & 8 )
        //   {
        //     display(um,sPrintF("um start of advance, t=%8.2e",t),debugFile,"%8.2e ");
        //     display(u,sPrintF("u start of advance, t=%8.2e",t),debugFile,"%8.2e ");
        //   }
        // }

                lambda = lambdaGrid(grid);
                mu = muGrid(grid);
                c1=(mu+lambda)/rho, c2= mu/rho;
                
                if( numberOfStepsTaken<1 ) 
                    printF(" advanceIncompressible:INFO mu=%8.2e for grid=%i (%s) \n",mu,grid,(const char*)cg[grid].getName());
                

                bool useOpt=true; // true;
                const bool isRectangular=mg.isRectangular();

        // ----- Add upwind dissipation ----

                OV_GET_SERIAL_ARRAY(real,um,umLocal);
                OV_GET_SERIAL_ARRAY(real,u ,uLocal);
                OV_GET_SERIAL_ARRAY(real,un,unLocal);
      
                {
                    real timeAdv=getCPU();
                    int gridType = isRectangular? 0 : 1;
                    const real adc=artificialDissipation; // do not scale *wdh* 090216
                    real *umptr=umLocal.getDataPointer();
                    real *uptr =uLocal.getDataPointer();
                    real *unptr=unLocal.getDataPointer();
                    real *vtptr = uptr; 
                    real *vnptr = uptr; 
                    real *vptr  = uptr; 
                    real *v1ptr = uptr; 
                    real *v2ptr = uptr; 
                    real *v3ptr = uptr; 
                    real *v4ptr = uptr; 
                    real *vt1ptr = uptr; 
                    real *vt2ptr = uptr; 
                    real *vt3ptr = uptr; 
                    real *vt4ptr = uptr; 
                        const bool addForcing =false;
                        real *fptr  =  uptr; // not used for upwind
          // assert( !useVariableDissipation || variableDissipation!=NULL );
          // real *pVarDis = useVariableDissipation ? varDis.getDataPointer() : uptr;
                    const intArray & mask = mg.mask();
                    OV_GET_SERIAL_ARRAY_CONST(int,mask,maskLocal)
          // #ifdef USE_PPP
          //     intSerialArray maskLocal;  getLocalArrayWithGhostBoundaries(mask,maskLocal);
          // #else
          //     const intSerialArray & maskLocal = mask; 
          // #endif
                    real *rxptr;
                    if( isRectangular )
                    {
                        rxptr=uptr;
                    }
                    else
                    {
                    #ifdef USE_PPP
                                realSerialArray rxLocal; getLocalArrayWithGhostBoundaries(mg.inverseVertexDerivative(),rxLocal);
                    #else
                                const realSerialArray & rxLocal=mg.inverseVertexDerivative();
                    #endif
                        rxptr = rxLocal.getDataPointer();
                    }
                    const bool centerNeeded=forcingOption==twilightZoneForcing || (forcingOption==planeWaveBoundaryForcing); // **************** fix this 
                    #ifdef USE_PPP
                        realSerialArray xy;
                        if( centerNeeded ) getLocalArrayWithGhostBoundaries(mg.center(),xy);
                    #else
                        const realSerialArray & xy = centerNeeded ? mg.center() : umLocal;
                    #endif
                    real *xyptr=xy.getDataPointer();
                    int maskNull[1];
                    int *maskptr = maskLocal.getDataPointer();
          // realSerialArray *dis = NULL;
          // real *pdis=uptr;
          // Macro to extract the pointers to the variable material property arrays
           // --- Variable material properies ---
                      GridMaterialProperties::MaterialFormatEnum materialFormat = GridMaterialProperties::constantMaterialProperties;
                      int ndMatProp=1;  // for piecewise constant materials, this is the leading dimension of the matVal array
                      int *matIndexPtr=maskptr;  // if not used, point to mask
                      real*matValPtr=uptr;       // if not used, point to u
                      if( parameters.dbase.get<int>("variableMaterialPropertiesOption")!=0 )
                      {
             // Material properties do vary 
                          std::vector<GridMaterialProperties> & materialProperties = 
                                    parameters.dbase.get<std::vector<GridMaterialProperties> >("materialProperties");
                          GridMaterialProperties & matProp = materialProperties[grid];
                          materialFormat = matProp.getMaterialFormat();
                          if( materialFormat==GridMaterialProperties::piecewiseConstantMaterialProperties )
                          {
                                    IntegerArray & matIndex = matProp.getMaterialIndexArray();
                              matIndexPtr = matIndex.getDataPointer();
                          }
                          RealArray & matVal = matProp.getMaterialValuesArray();
                          matValPtr = matVal.getDataPointer();
                          ndMatProp = matVal.getLength(0);  
             // ::display(matVal,"matVal");
                      }
                    bool twilightZone = forcingOption==twilightZoneForcing; 
                    int combineDissipationWithAdvance=0;  // not used anymore
                    real tForce = t;   // evaluate forcing at this time 
                        int option=1; //  0=update solution, 1=add upwind dissipation
                        if( false )
                            printF(">>>> advanceIncompressible: add upwind dissipation at globalStepNumber=%d, t=%9.3e <<<<<\n",globalStepNumber,t);
                    int ipar[]={option,
                                            gridType,
                                            orderOfAccuracyInSpace,
                                            orderOfTimeAccuracy,
                                            (int)addForcing,
                                            (int)twilightZone,
                                            uc,
                                            vc,
                                            wc,
                                            useWhereMask,
                                            (int)timeSteppingMethodSm,
                                            (int)useVariableDissipation,
                                            (int)useConservative,           // 12 
                                            combineDissipationWithAdvance,
                                            debug,
                                            computeUt,
                                            materialFormat,                 // 16 
                                            myid,
                                            pc,                             // 18 
                                            upwindSOS,
                                            correction
                                            };    
                    real dx[3]={1.,1.,1.};
                    if( isRectangular )
                        mg.getDeltaX(dx);
                    real rpar[30];
                    rpar[ 0]=dt;
                    rpar[ 1]=dx[0];
                    rpar[ 2]=dx[1];
                    rpar[ 3]=dx[2];
                    rpar[ 4]=adc;
                    rpar[ 5]=mg.gridSpacing(0);
                    rpar[ 6]=mg.gridSpacing(1);
                    rpar[ 7]=mg.gridSpacing(2);
                    rpar[ 8]=rho;
                    rpar[ 9]=mu;
                    rpar[10]=kx; // for plane wave scattering
                    rpar[11]=ky;
                    rpar[12]=kz;
                    rpar[13]=(real &)parameters.dbase.get<OGFunction* >("exactSolution");  // twilight zone pointer, ep
                    rpar[14]=tForce;
                    rpar[15]=parameters.dbase.get<real>("dtOld");  // dt used on the previous step
                    rpar[16]=rho;
                    rpar[17]=mu;
                    rpar[18]=lambda;
                    rpar[20]=0.;  // return cpu for dissipation
          // printF("** AdvanceIncompressible: gridType=%i, isRectangular=%i\n",gridType,(int)isRectangular);
                    bool ok = ParallelUtility::getLocalArrayBounds(u,uLocal,I1,I2,I3);
                    int ierr=0;
                    if( ok )
                    {
            // if( adc>0. && !combineDissipationWithAdvance )
            // {
            //   // create a temp array to hold the artificial dissipation
            //   dis = new realSerialArray(uLocal.dimension(0),uLocal.dimension(1),uLocal.dimension(2),uLocal.dimension(3));
            //   pdis = dis->getDataPointer();
            //   assert( pdis!=NULL );
            //   sizeOfLocalArraysForAdvance=max(sizeOfLocalArraysForAdvance,(double)(dis->elementCount()*sizeof(real)));
            // }
            // real timeAdv=getCPU();
            // if( combineDissipationWithAdvance )
            // {
            //   assert( umptr!=unptr );
            // }
                        advIsm(mg.numberOfDimensions(),
                                      I1.getBase(),I1.getBound(),I2.getBase(),I2.getBound(),I3.getBase(),I3.getBound(),
                                      uLocal.getBase(0),uLocal.getBound(0),
                                      uLocal.getBase(1),uLocal.getBound(1),
                                      uLocal.getBase(2),uLocal.getBound(2),
                                      uLocal.getBase(3),uLocal.getBound(3),
                                      *maskptr,*rxptr,*xyptr,  
                                      *umptr,*uptr,*unptr, *fptr, 
                                      *vtptr, *vnptr, *vptr,
                                      *v1ptr, *v2ptr, *v3ptr, *v4ptr,
                                      *vt1ptr,*vt2ptr,*vt3ptr,*vt4ptr,
                                      mg.gridIndexRange(0,0), mg.boundaryCondition(0,0),  ipar[0], rpar[0], ierr );
                    }
                    timeAdv=getCPU()-timeAdv;
                        timing(parameters.dbase.get<int>("timeForUpwindDissipation"))+=timeAdv;
                    if( debug & 8 )
                    {
                        display(unLocal,sPrintF("unLocal after advIncompressibleSolidMechanics, processor=%i before BC's t=%8.2e",
                                                                        Communication_Manager::My_Process_Number,t),pDebugFile,"%8.2e ");
                        display(un,sPrintF("un after advIncompressibleSolidMechanics, before BC's t=%8.2e",t),debugFile,"%8.2e ");
                    }
                    if( isRectangular )   
                        timing(parameters.dbase.get<int>("timeForAdvanceRectangularGrids"))+=getCPU()-time0;
                    else
                        timing(parameters.dbase.get<int>("timeForAdvanceCurvilinearGrids"))+=getCPU()-time0;
                    if( debug & 8 )
                    {
                        display(u,sPrintF("u after advOpt and updateGhost, t=%8.2e",t),debugFile,"%8.2e ");
                    }
                }


            } // end grid

                if( applyBC )
                {
                    if( true ||   // *wdh* 091205 -- interpolate will call periodicUpdate and updateGhost 
                            cg.numberOfComponentGrids()>1 )
                    {
                        real timei=getCPU();
                        if( debug & 4 )
                            gf[next].u.display(sPrintF("Cgsm::advanceIncompressible: gf[next].u before interpolate, t=%8.2e",t+dt),debugFile,"%8.2e ");
            // --- Note: interpolate performs a periodicUpdate and updateGhostBoundaries even if there is only one grid
                        gf[next].u.interpolate();
                        if( debug & 4 )
                            gf[next].u.display(sPrintF("Cgsm::advanceIncompressible: gf[next].u after interpolate, t=%8.2e",t+dt),debugFile,"%8.2e ");
                        if( debug & 4 )
                        {
                            getErrors( next,t+dt,dt,sPrintF("\n ************** advance Errors after interpolate t=%9.3e ******\n",t+dt));
                        }
                        timing(parameters.dbase.get<int>("timeForInterpolate"))+=getCPU()-timei;
                    }
          // ============= Boundary Conditions =============
                    int option=0; // not used.
                    applyBoundaryConditions( option, dt, next,current ); // apply BC to "next" (current=previous time step)
                    if( debug & 8 )  // & 64
                    {
                        gf[next].u.display(sPrintF("Cgsm::advanceIncompressible: gf[next].u after applyBC, t=%8.2e",gf[next].t),debugFile,"%8.2e ");
                    }
                    if( debug & 4 )
                      {
                          getErrors( next,t+dt,dt,sPrintF("\n ************** advanceIncompressible Errors after applyBC t=%9.3e ******\n",t+dt));
                      }    
          // ---- assign values at material interfaces ------
          // *** this does nothing currently ***
                    assignInterfaceBoundaryConditions( next, t, dt );  // is this the right place to do this?
                }
        } // end if upwindSOS



        if( correction<=(numberOfCorrections-1) || !skipLastPressureSolve)
        {
      // -----------------------------------------------------
      // ---------------- PRESSURE SOLVE ---------------------
      // -----------------------------------------------------
            solveForPressure( current, t+dt, dt );

            if( orderOfAccuracyInSpace==4 )
            { // re-apply BC's that depend on p 
                    if( applyBC )
                    {
                        if( true ||   // *wdh* 091205 -- interpolate will call periodicUpdate and updateGhost 
                                cg.numberOfComponentGrids()>1 )
                        {
                            real timei=getCPU();
                            if( debug & 4 )
                                gf[next].u.display(sPrintF("Cgsm::advanceIncompressible: gf[next].u before interpolate, t=%8.2e",t+dt),debugFile,"%8.2e ");
              // --- Note: interpolate performs a periodicUpdate and updateGhostBoundaries even if there is only one grid
                            gf[next].u.interpolate();
                            if( debug & 4 )
                                gf[next].u.display(sPrintF("Cgsm::advanceIncompressible: gf[next].u after interpolate, t=%8.2e",t+dt),debugFile,"%8.2e ");
                            if( debug & 4 )
                            {
                                getErrors( next,t+dt,dt,sPrintF("\n ************** advance Errors after interpolate t=%9.3e ******\n",t+dt));
                            }
                            timing(parameters.dbase.get<int>("timeForInterpolate"))+=getCPU()-timei;
                        }
            // ============= Boundary Conditions =============
                        int option=0; // not used.
                        applyBoundaryConditions( option, dt, next,current ); // apply BC to "next" (current=previous time step)
                        if( debug & 8 )  // & 64
                        {
                            gf[next].u.display(sPrintF("Cgsm::advanceIncompressible: gf[next].u after applyBC, t=%8.2e",gf[next].t),debugFile,"%8.2e ");
                        }
                        if( debug & 4 )
                          {
                              getErrors( next,t+dt,dt,sPrintF("\n ************** advanceIncompressible Errors after applyBC t=%9.3e ******\n",t+dt));
                          }    
            // ---- assign values at material interfaces ------
            // *** this does nothing currently ***
                        assignInterfaceBoundaryConditions( next, t, dt );  // is this the right place to do this?
                    }
            }
        }
        else
        {
            printF("advanceIncompressible: SKIP last pressure solve\n");
        }


    } // end correction


    if( projectLinearModeForILE &&
            globalStepNumber>0 && (globalStepNumber % projectLinearModeFrequency) == 0 )
    {
        projectLinearMode( current, t, dt );
    }


    parameters.dbase.get<real>("dtOld")=dt;  // set dtOld   

    checkArrays("advanceIncompressible:end");
    
}


// =============================================================================
/// \brief Advance the solution for
///       Incompressible Linear Elastic Equations
/// using
///       METHOD OF LINES 
///
// =============================================================================
void Cgsm::
advanceIncompressibleMethodOfLines( int current, real t, real dt,
                                                                      RealCompositeGridFunction *ut /* = NULL */, 
                                                                      real tForce /* = 0. */ )
{   

    int & globalStepNumber = parameters.dbase.get<int >("globalStepNumber");
    if( t<2*dt || debug() & 4 )
        printF("\n @@@@@@@ Cgsm::advanceIncompressibleMethodOfLines called at globalStepNumber=%d, t=%9.3e\n",globalStepNumber,t);

    checkArrays("advanceIncompressibleMethodOfLines:start");


    FILE *& debugFile  =parameters.dbase.get<FILE* >("debugFile");
    FILE *& logFile    =parameters.dbase.get<FILE* >("logFile");
    FILE *& pDebugFile =parameters.dbase.get<FILE* >("pDebugFile");
    
    const int numberOfDimensions       = cg.numberOfDimensions();
    const int numberOfComponentGrids   = cg.numberOfComponentGrids();
    const int & numberOfComponents     = parameters.dbase.get<int >("numberOfComponents");
    const int & uc                     = parameters.dbase.get<int >("uc");
    const int & vc                     = parameters.dbase.get<int >("vc");
    const int & wc                     = parameters.dbase.get<int >("wc");
    const int & rc                     = parameters.dbase.get<int >("rc");
    const int & tc                     = parameters.dbase.get<int >("tc");
    const int & pc                     = parameters.dbase.get<int >("pc");
    const int & orderOfAccuracyInSpace = parameters.dbase.get<int>("orderOfAccuracy");
    const int & orderOfTimeAccuracy    = parameters.dbase.get<int>("orderOfTimeAccuracy");
    const int & numberOfCorrections    = parameters.dbase.get<int>("numberOfCorrections"); 
    const int & skipLastPressureSolve  = parameters.dbase.get<int>("skipLastPressureSolve");  // For ILE predictor-corrector schemes
  

    SmParameters::TimeSteppingMethodSm & timeSteppingMethodSm = 
                                                                      parameters.dbase.get<SmParameters::TimeSteppingMethodSm>("timeSteppingMethodSm");
    RealArray & timing = parameters.dbase.get<RealArray >("timing");


  // --- Here is where we save the velocity v, and acceleration vt ------
    realCompositeGridFunction *& vgf  = parameters.dbase.get<realCompositeGridFunction*>("vgf");
    realCompositeGridFunction *& vtgf = parameters.dbase.get<realCompositeGridFunction*>("vtgf");

    const int & numberOfVelocityFunctions     = parameters.dbase.get<int>("numberOfVelocityFunctions");
    const int & numberOfAccelerationFunctions = parameters.dbase.get<int>("numberOfAccelerationFunctions");


    Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
    Range C=numberOfComponents;

    const int prev = (current-1+numberOfTimeLevels) % numberOfTimeLevels;
    const int next = (current+1)                    % numberOfTimeLevels;

  // Here are the index pointers into the velocity and acceleration arrays 
    int & currentVelocity   = parameters.dbase.get<int>("currentVelocity");
    const int nextVelocity  = (currentVelocity +1 + numberOfVelocityFunctions) % numberOfVelocityFunctions;
    const int prevVelocity  = (currentVelocity -1 + numberOfVelocityFunctions) % numberOfVelocityFunctions;
    const int prev2Velocity = (currentVelocity -2 + 2*numberOfVelocityFunctions) % numberOfVelocityFunctions;
    const int prev3Velocity = (currentVelocity -3 + 2*numberOfVelocityFunctions) % numberOfVelocityFunctions;

  // Here are the index pointers into the pressure arrays (for time extrapolation)
    realCompositeGridFunction *& pgf      = parameters.dbase.get<realCompositeGridFunction*>("pgf");
    const int & numberOfPressureFunctions = parameters.dbase.get<int>("numberOfPressureFunctions");
    int & currentPressure                 = parameters.dbase.get<int>("currentPressure");
    const int prevPressure                = (currentPressure - 1 + numberOfPressureFunctions) % max(1,numberOfPressureFunctions);
    const int nextPressure                = (currentPressure + 1 + numberOfPressureFunctions) % max(1,numberOfPressureFunctions);

    real & rho=parameters.dbase.get<real>("rho");
    real & mu = parameters.dbase.get<real>("mu");
    real & lambda = parameters.dbase.get<real>("lambda");
    RealArray & muGrid = parameters.dbase.get<RealArray>("muGrid");
    RealArray & lambdaGrid = parameters.dbase.get<RealArray>("lambdaGrid");
    bool & gridHasMaterialInterfaces = parameters.dbase.get<bool>("gridHasMaterialInterfaces");
    int & debug = parameters.dbase.get<int >("debug");

    const int & upwindSOS = parameters.dbase.get<int>("upwindSOS"); 

    const real cMax=max(lambdaGrid+muGrid)/rho;

    const int computeUt = ut != NULL;  // compute u.tt


    
    const int useWhereMask = numberOfComponentGrids>1;  

    sizeOfLocalArraysForAdvance=0.;


    if( timeSteppingMethodSm != SmParameters::adamsPredictorCorrector2 &&
            timeSteppingMethodSm != SmParameters::adamsPredictorCorrector4 )
    {
        printF("advanceIncompressibleMethodOfLines: not implemented for timeSteppingMethodSm=%d\n",(int)timeSteppingMethodSm);
        printF("Currently the following MOL schemes are supported:\n");
        printF("   advanceAdamsPredictorCorrector2=%d\n",(int)SmParameters::adamsPredictorCorrector2);
        printF("   advanceAdamsPredictorCorrector4=%d\n",(int)SmParameters::adamsPredictorCorrector4);
        OV_ABORT("error"); 
    }


  // -------------- LOOP FOR PREDICTOR AND CORRECTOR ---------------------
    for( int correction=0; correction<=numberOfCorrections; correction++ )
    {

    // ----------------- ADVANCE SOLUTION -----------------

    // --------------- START LOOP OVER GRIDS -------------------
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            real time0=getCPU();

            MappedGrid & mg = cg[grid];
            MappedGridOperators & mgop = (*cgop)[grid];

            getIndex(mg.gridIndexRange(),I1,I2,I3);
        
      // getBoundsForPML( mg,Iv );

            realMappedGridFunction & fieldPrev    =gf[prev].u[grid];
            realMappedGridFunction & fieldCurrent =gf[current].u[grid];
            realMappedGridFunction & fieldNext    =gf[next].u[grid];


            realArray & um = fieldPrev;
            realArray & u  = fieldCurrent;
            realArray & un = fieldNext;

            if( debug & 4 )
            {
                fPrintF(debugFile," **** start of advance, t=%8.2e\n",t);
                fprintf(pDebugFile," **** start of advance, t=%8.2e\n",t);
                
                if( debug & 8 )
                {
                    display(um,sPrintF("um start of advance, t=%8.2e",t),debugFile,"%8.2e ");
                    display(u,sPrintF("u start of advance, t=%8.2e",t),debugFile,"%8.2e ");
                }
            }

            lambda = lambdaGrid(grid);
            mu = muGrid(grid);
            c1=(mu+lambda)/rho, c2= mu/rho;
            
            if( numberOfStepsTaken<1 ) 
                printF(" advanceIncompressibleMethodOfLines:INFO mu=%8.2e for grid=%i (%s) \n",mu,grid,(const char*)cg[grid].getName());
            
      // const real dtsq=dt*dt; 
      // const real adc=artificialDissipation*SQR(cMax); // scale dissipation by c^2 *wdh* 041103


            bool useOpt=true; // true;
            const bool isRectangular=mg.isRectangular();
            if( !isRectangular )
            {
         // *** PUT THIS IN SETUP SOMEWHERE ***** 
                real timea=getCPU();
                mg.update( MappedGrid::THEinverseVertexDerivative | MappedGrid::THEinverseCenterDerivative );
                timea=getCPU()-timea;
                timing(parameters.dbase.get<int>("timeForInitialize"))+=timea;

                time0+=timea;  // do not include with time for curvilinear
            }

      //  ::display(lap,"lap","%8.1e ");

      // --- get any forcing ----
                realArray f;
                const bool addForcing = forcingOption!=noForcing && 
                                                                forcingOption!=planeWaveBoundaryForcing && 
                                                                forcingOption!=twilightZoneForcing;
                if( addForcing )
                {
                    Index D1,D2,D3;
                    getIndex(mg.dimension(),D1,D2,D3);
                    f.partition(mg.getPartition());
                    f.redim(D1,D2,D3,C);  // could use some other array for work space ??
                    int option=1;  
                    if( addForcing )
                    {
            //kkc getForcing is called from advance but is also called from assignIC and getErrors.
            //    we have to add the timing in external to getForcing to avoid double counting the time
            //    in assignIC and getErrors
                        real timef = getCPU();
                        getForcing( next, grid,f,t,dt,option );
                        timing(parameters.dbase.get<int>("timeForForcing")) += getCPU()-timef;
                    }
                }
                OV_GET_SERIAL_ARRAY(real,f ,fLocal);

            OV_GET_SERIAL_ARRAY(real,um,umLocal);
            OV_GET_SERIAL_ARRAY(real,u ,uLocal);
            OV_GET_SERIAL_ARRAY(real,un,unLocal);
      

      // ------ getVt  -------
            {
                real timeAdv=getCPU();
                int gridType = isRectangular? 0 : 1;
                const real adc=artificialDissipation; // do not scale *wdh* 090216
                real *umptr=umLocal.getDataPointer();
                real *uptr =uLocal.getDataPointer();
                real *unptr=unLocal.getDataPointer();
                real *vtptr = uptr; 
                real *vnptr = uptr; 
                real *vptr  = uptr; 
                real *v1ptr = uptr; 
                real *v2ptr = uptr; 
                real *v3ptr = uptr; 
                real *v4ptr = uptr; 
                real *vt1ptr = uptr; 
                real *vt2ptr = uptr; 
                real *vt3ptr = uptr; 
                real *vt4ptr = uptr; 
                    real *fptr   = addForcing ? fLocal.getDataPointer() : uptr;
        // assert( !useVariableDissipation || variableDissipation!=NULL );
        // real *pVarDis = useVariableDissipation ? varDis.getDataPointer() : uptr;
                const intArray & mask = mg.mask();
                OV_GET_SERIAL_ARRAY_CONST(int,mask,maskLocal)
        // #ifdef USE_PPP
        //     intSerialArray maskLocal;  getLocalArrayWithGhostBoundaries(mask,maskLocal);
        // #else
        //     const intSerialArray & maskLocal = mask; 
        // #endif
                real *rxptr;
                if( isRectangular )
                {
                    rxptr=uptr;
                }
                else
                {
                #ifdef USE_PPP
                            realSerialArray rxLocal; getLocalArrayWithGhostBoundaries(mg.inverseVertexDerivative(),rxLocal);
                #else
                            const realSerialArray & rxLocal=mg.inverseVertexDerivative();
                #endif
                    rxptr = rxLocal.getDataPointer();
                }
                const bool centerNeeded=forcingOption==twilightZoneForcing || (forcingOption==planeWaveBoundaryForcing); // **************** fix this 
                #ifdef USE_PPP
                    realSerialArray xy;
                    if( centerNeeded ) getLocalArrayWithGhostBoundaries(mg.center(),xy);
                #else
                    const realSerialArray & xy = centerNeeded ? mg.center() : umLocal;
                #endif
                real *xyptr=xy.getDataPointer();
                int maskNull[1];
                int *maskptr = maskLocal.getDataPointer();
        // realSerialArray *dis = NULL;
        // real *pdis=uptr;
        // Macro to extract the pointers to the variable material property arrays
         // --- Variable material properies ---
                  GridMaterialProperties::MaterialFormatEnum materialFormat = GridMaterialProperties::constantMaterialProperties;
                  int ndMatProp=1;  // for piecewise constant materials, this is the leading dimension of the matVal array
                  int *matIndexPtr=maskptr;  // if not used, point to mask
                  real*matValPtr=uptr;       // if not used, point to u
                  if( parameters.dbase.get<int>("variableMaterialPropertiesOption")!=0 )
                  {
           // Material properties do vary 
                      std::vector<GridMaterialProperties> & materialProperties = 
                                parameters.dbase.get<std::vector<GridMaterialProperties> >("materialProperties");
                      GridMaterialProperties & matProp = materialProperties[grid];
                      materialFormat = matProp.getMaterialFormat();
                      if( materialFormat==GridMaterialProperties::piecewiseConstantMaterialProperties )
                      {
                                IntegerArray & matIndex = matProp.getMaterialIndexArray();
                          matIndexPtr = matIndex.getDataPointer();
                      }
                      RealArray & matVal = matProp.getMaterialValuesArray();
                      matValPtr = matVal.getDataPointer();
                      ndMatProp = matVal.getLength(0);  
           // ::display(matVal,"matVal");
                  }
                bool twilightZone = forcingOption==twilightZoneForcing; 
                int combineDissipationWithAdvance=0;  // not used anymore
                real tForce = t;   // evaluate forcing at this time 
          // MOL
                    int option=2; // compute vt for MOL schemes
                    if( correction==0 )
                    {
                        if( false )
                            printF("++++++ GetVt for MOL: PREDICTOR: currentVelocity=%d for t=%9.3e, tForce=%9.3e\n",currentVelocity,t,tForce);
                        vtptr = vtgf[currentVelocity][grid].getDataPointer();  // evaluate v.t(t)  
                    }
                    else
                    {
                        tForce = t+dt;
                        if( false )
                            printF("++++++ GetVt for MOL: CORRECTOR: nextVelocity=%d for t+dt=%9.3e, tForce=%9.3e\n",nextVelocity,t+dt,tForce);
                        uptr = unptr; // eval v.t using u(t+dt)
                        vtptr = vtgf[nextVelocity][grid].getDataPointer();     // evaluate v.t(t+dt)
                    }
                int ipar[]={option,
                                        gridType,
                                        orderOfAccuracyInSpace,
                                        orderOfTimeAccuracy,
                                        (int)addForcing,
                                        (int)twilightZone,
                                        uc,
                                        vc,
                                        wc,
                                        useWhereMask,
                                        (int)timeSteppingMethodSm,
                                        (int)useVariableDissipation,
                                        (int)useConservative,           // 12 
                                        combineDissipationWithAdvance,
                                        debug,
                                        computeUt,
                                        materialFormat,                 // 16 
                                        myid,
                                        pc,                             // 18 
                                        upwindSOS,
                                        correction
                                        };    
                real dx[3]={1.,1.,1.};
                if( isRectangular )
                    mg.getDeltaX(dx);
                real rpar[30];
                rpar[ 0]=dt;
                rpar[ 1]=dx[0];
                rpar[ 2]=dx[1];
                rpar[ 3]=dx[2];
                rpar[ 4]=adc;
                rpar[ 5]=mg.gridSpacing(0);
                rpar[ 6]=mg.gridSpacing(1);
                rpar[ 7]=mg.gridSpacing(2);
                rpar[ 8]=rho;
                rpar[ 9]=mu;
                rpar[10]=kx; // for plane wave scattering
                rpar[11]=ky;
                rpar[12]=kz;
                rpar[13]=(real &)parameters.dbase.get<OGFunction* >("exactSolution");  // twilight zone pointer, ep
                rpar[14]=tForce;
                rpar[15]=parameters.dbase.get<real>("dtOld");  // dt used on the previous step
                rpar[16]=rho;
                rpar[17]=mu;
                rpar[18]=lambda;
                rpar[20]=0.;  // return cpu for dissipation
        // printF("** AdvanceIncompressible: gridType=%i, isRectangular=%i\n",gridType,(int)isRectangular);
                bool ok = ParallelUtility::getLocalArrayBounds(u,uLocal,I1,I2,I3);
                int ierr=0;
                if( ok )
                {
          // if( adc>0. && !combineDissipationWithAdvance )
          // {
          //   // create a temp array to hold the artificial dissipation
          //   dis = new realSerialArray(uLocal.dimension(0),uLocal.dimension(1),uLocal.dimension(2),uLocal.dimension(3));
          //   pdis = dis->getDataPointer();
          //   assert( pdis!=NULL );
          //   sizeOfLocalArraysForAdvance=max(sizeOfLocalArraysForAdvance,(double)(dis->elementCount()*sizeof(real)));
          // }
          // real timeAdv=getCPU();
          // if( combineDissipationWithAdvance )
          // {
          //   assert( umptr!=unptr );
          // }
                    advIsm(mg.numberOfDimensions(),
                                  I1.getBase(),I1.getBound(),I2.getBase(),I2.getBound(),I3.getBase(),I3.getBound(),
                                  uLocal.getBase(0),uLocal.getBound(0),
                                  uLocal.getBase(1),uLocal.getBound(1),
                                  uLocal.getBase(2),uLocal.getBound(2),
                                  uLocal.getBase(3),uLocal.getBound(3),
                                  *maskptr,*rxptr,*xyptr,  
                                  *umptr,*uptr,*unptr, *fptr, 
                                  *vtptr, *vnptr, *vptr,
                                  *v1ptr, *v2ptr, *v3ptr, *v4ptr,
                                  *vt1ptr,*vt2ptr,*vt3ptr,*vt4ptr,
                                  mg.gridIndexRange(0,0), mg.boundaryCondition(0,0),  ipar[0], rpar[0], ierr );
                }
                timeAdv=getCPU()-timeAdv;
                    timing(parameters.dbase.get<int>("timeForAdvOpt"))+=timeAdv;
                if( debug & 8 )
                {
                    display(unLocal,sPrintF("unLocal after advIncompressibleSolidMechanics, processor=%i before BC's t=%8.2e",
                                                                    Communication_Manager::My_Process_Number,t),pDebugFile,"%8.2e ");
                    display(un,sPrintF("un after advIncompressibleSolidMechanics, before BC's t=%8.2e",t),debugFile,"%8.2e ");
                }
                if( isRectangular )   
                    timing(parameters.dbase.get<int>("timeForAdvanceRectangularGrids"))+=getCPU()-time0;
                else
                    timing(parameters.dbase.get<int>("timeForAdvanceCurvilinearGrids"))+=getCPU()-time0;
                if( debug & 8 )
                {
                    display(u,sPrintF("u after advOpt and updateGhost, t=%8.2e",t),debugFile,"%8.2e ");
                }
                if( debug & 16 )
                {
                    if( addForcing ) 
                        ::display(f,sPrintF("  *** advIncompressibleSolidMechanics: Here is the forcing f grid=%i t=%9.3e ********\n",grid,t),
                                            pDebugFile,"%7.4f ");
                } 
            }

      // ------ update MOL solution  -------
            {
                real timeAdv=getCPU();
                int gridType = isRectangular? 0 : 1;
                const real adc=artificialDissipation; // do not scale *wdh* 090216
                real *umptr=umLocal.getDataPointer();
                real *uptr =uLocal.getDataPointer();
                real *unptr=unLocal.getDataPointer();
                real *vtptr = uptr; 
                real *vnptr = uptr; 
                real *vptr  = uptr; 
                real *v1ptr = uptr; 
                real *v2ptr = uptr; 
                real *v3ptr = uptr; 
                real *v4ptr = uptr; 
                real *vt1ptr = uptr; 
                real *vt2ptr = uptr; 
                real *vt3ptr = uptr; 
                real *vt4ptr = uptr; 
                    real *fptr   = addForcing ? fLocal.getDataPointer() : uptr;
        // assert( !useVariableDissipation || variableDissipation!=NULL );
        // real *pVarDis = useVariableDissipation ? varDis.getDataPointer() : uptr;
                const intArray & mask = mg.mask();
                OV_GET_SERIAL_ARRAY_CONST(int,mask,maskLocal)
        // #ifdef USE_PPP
        //     intSerialArray maskLocal;  getLocalArrayWithGhostBoundaries(mask,maskLocal);
        // #else
        //     const intSerialArray & maskLocal = mask; 
        // #endif
                real *rxptr;
                if( isRectangular )
                {
                    rxptr=uptr;
                }
                else
                {
                #ifdef USE_PPP
                            realSerialArray rxLocal; getLocalArrayWithGhostBoundaries(mg.inverseVertexDerivative(),rxLocal);
                #else
                            const realSerialArray & rxLocal=mg.inverseVertexDerivative();
                #endif
                    rxptr = rxLocal.getDataPointer();
                }
                const bool centerNeeded=forcingOption==twilightZoneForcing || (forcingOption==planeWaveBoundaryForcing); // **************** fix this 
                #ifdef USE_PPP
                    realSerialArray xy;
                    if( centerNeeded ) getLocalArrayWithGhostBoundaries(mg.center(),xy);
                #else
                    const realSerialArray & xy = centerNeeded ? mg.center() : umLocal;
                #endif
                real *xyptr=xy.getDataPointer();
                int maskNull[1];
                int *maskptr = maskLocal.getDataPointer();
        // realSerialArray *dis = NULL;
        // real *pdis=uptr;
        // Macro to extract the pointers to the variable material property arrays
         // --- Variable material properies ---
                  GridMaterialProperties::MaterialFormatEnum materialFormat = GridMaterialProperties::constantMaterialProperties;
                  int ndMatProp=1;  // for piecewise constant materials, this is the leading dimension of the matVal array
                  int *matIndexPtr=maskptr;  // if not used, point to mask
                  real*matValPtr=uptr;       // if not used, point to u
                  if( parameters.dbase.get<int>("variableMaterialPropertiesOption")!=0 )
                  {
           // Material properties do vary 
                      std::vector<GridMaterialProperties> & materialProperties = 
                                parameters.dbase.get<std::vector<GridMaterialProperties> >("materialProperties");
                      GridMaterialProperties & matProp = materialProperties[grid];
                      materialFormat = matProp.getMaterialFormat();
                      if( materialFormat==GridMaterialProperties::piecewiseConstantMaterialProperties )
                      {
                                IntegerArray & matIndex = matProp.getMaterialIndexArray();
                          matIndexPtr = matIndex.getDataPointer();
                      }
                      RealArray & matVal = matProp.getMaterialValuesArray();
                      matValPtr = matVal.getDataPointer();
                      ndMatProp = matVal.getLength(0);  
           // ::display(matVal,"matVal");
                  }
                bool twilightZone = forcingOption==twilightZoneForcing; 
                int combineDissipationWithAdvance=0;  // not used anymore
                real tForce = t;   // evaluate forcing at this time 
          // ----  MOL update the solution ----
                    int option=3; // update MOL solution using vt
                    if( false )
                          printF("++++++ Update MOL: Fill un[next=%d], u[cur=%d], AND vn[nextV=%d] from v[curV=%d] t=%9.3e\n",next,current,nextVelocity,currentVelocity,t);
                    vnptr  = vgf[nextVelocity   ][grid].getDataPointer();
                    vptr   = vgf[currentVelocity][grid].getDataPointer();
                    if( correction==0 )
                    { // predictor:
            //   un = u +  dt*( c1*v1  + c2*v2  + ... )
            //   vn = v +  dt*( c1*vt1 + c2*vt2 + ... )
                        if( false )
                              printF("++++++ Update MOL: PREDICTOR currentVelocity=%d, prevVelocity=%d, prev2V=%d, prev3V=%d, t=%9.3e\n",
                                              currentVelocity,prevVelocity,prev2Velocity,prev3Velocity,t);
                        v1ptr  = vgf[currentVelocity][grid].getDataPointer();
                        v2ptr  = vgf[prevVelocity   ][grid].getDataPointer();
                        v3ptr  = vgf[prev2Velocity  ][grid].getDataPointer();
                        v4ptr  = vgf[prev3Velocity  ][grid].getDataPointer();
                        vt1ptr = vtgf[currentVelocity][grid].getDataPointer();
                        vt2ptr = vtgf[prevVelocity   ][grid].getDataPointer();
                        vt3ptr = vtgf[prev2Velocity  ][grid].getDataPointer();
                        vt4ptr = vtgf[prev3Velocity  ][grid].getDataPointer();
                    }
                    else
                    { // corrector 
            //   un = u +  dt*( c1*v[next] + c2*v[cur] + ... )
                        if( false )
                            printF("++++++ Update MOL: CORRECTOR correction=%d, nextVelocity=%d,, currentVelocity=%d, prevV=%d, prev2=%d, t=%9.3e\n",
                                          correction,nextVelocity,currentVelocity,prevVelocity,prev2Velocity,t);
                        v1ptr  = vgf[nextVelocity   ][grid].getDataPointer();
                        v2ptr  = vgf[currentVelocity][grid].getDataPointer();
                        v3ptr  = vgf[prevVelocity   ][grid].getDataPointer();
                        v4ptr  = vgf[prev2Velocity  ][grid].getDataPointer();
                        vt1ptr = vtgf[nextVelocity   ][grid].getDataPointer();
                        vt2ptr = vtgf[currentVelocity][grid].getDataPointer();
                        vt3ptr = vtgf[prevVelocity   ][grid].getDataPointer();
                        vt4ptr = vtgf[prev2Velocity  ][grid].getDataPointer();
                    }
                int ipar[]={option,
                                        gridType,
                                        orderOfAccuracyInSpace,
                                        orderOfTimeAccuracy,
                                        (int)addForcing,
                                        (int)twilightZone,
                                        uc,
                                        vc,
                                        wc,
                                        useWhereMask,
                                        (int)timeSteppingMethodSm,
                                        (int)useVariableDissipation,
                                        (int)useConservative,           // 12 
                                        combineDissipationWithAdvance,
                                        debug,
                                        computeUt,
                                        materialFormat,                 // 16 
                                        myid,
                                        pc,                             // 18 
                                        upwindSOS,
                                        correction
                                        };    
                real dx[3]={1.,1.,1.};
                if( isRectangular )
                    mg.getDeltaX(dx);
                real rpar[30];
                rpar[ 0]=dt;
                rpar[ 1]=dx[0];
                rpar[ 2]=dx[1];
                rpar[ 3]=dx[2];
                rpar[ 4]=adc;
                rpar[ 5]=mg.gridSpacing(0);
                rpar[ 6]=mg.gridSpacing(1);
                rpar[ 7]=mg.gridSpacing(2);
                rpar[ 8]=rho;
                rpar[ 9]=mu;
                rpar[10]=kx; // for plane wave scattering
                rpar[11]=ky;
                rpar[12]=kz;
                rpar[13]=(real &)parameters.dbase.get<OGFunction* >("exactSolution");  // twilight zone pointer, ep
                rpar[14]=tForce;
                rpar[15]=parameters.dbase.get<real>("dtOld");  // dt used on the previous step
                rpar[16]=rho;
                rpar[17]=mu;
                rpar[18]=lambda;
                rpar[20]=0.;  // return cpu for dissipation
        // printF("** AdvanceIncompressible: gridType=%i, isRectangular=%i\n",gridType,(int)isRectangular);
                bool ok = ParallelUtility::getLocalArrayBounds(u,uLocal,I1,I2,I3);
                int ierr=0;
                if( ok )
                {
          // if( adc>0. && !combineDissipationWithAdvance )
          // {
          //   // create a temp array to hold the artificial dissipation
          //   dis = new realSerialArray(uLocal.dimension(0),uLocal.dimension(1),uLocal.dimension(2),uLocal.dimension(3));
          //   pdis = dis->getDataPointer();
          //   assert( pdis!=NULL );
          //   sizeOfLocalArraysForAdvance=max(sizeOfLocalArraysForAdvance,(double)(dis->elementCount()*sizeof(real)));
          // }
          // real timeAdv=getCPU();
          // if( combineDissipationWithAdvance )
          // {
          //   assert( umptr!=unptr );
          // }
                    advIsm(mg.numberOfDimensions(),
                                  I1.getBase(),I1.getBound(),I2.getBase(),I2.getBound(),I3.getBase(),I3.getBound(),
                                  uLocal.getBase(0),uLocal.getBound(0),
                                  uLocal.getBase(1),uLocal.getBound(1),
                                  uLocal.getBase(2),uLocal.getBound(2),
                                  uLocal.getBase(3),uLocal.getBound(3),
                                  *maskptr,*rxptr,*xyptr,  
                                  *umptr,*uptr,*unptr, *fptr, 
                                  *vtptr, *vnptr, *vptr,
                                  *v1ptr, *v2ptr, *v3ptr, *v4ptr,
                                  *vt1ptr,*vt2ptr,*vt3ptr,*vt4ptr,
                                  mg.gridIndexRange(0,0), mg.boundaryCondition(0,0),  ipar[0], rpar[0], ierr );
                }
                timeAdv=getCPU()-timeAdv;
                    timing(parameters.dbase.get<int>("timeForAdvOpt"))+=timeAdv;
                if( debug & 8 )
                {
                    display(unLocal,sPrintF("unLocal after advIncompressibleSolidMechanics, processor=%i before BC's t=%8.2e",
                                                                    Communication_Manager::My_Process_Number,t),pDebugFile,"%8.2e ");
                    display(un,sPrintF("un after advIncompressibleSolidMechanics, before BC's t=%8.2e",t),debugFile,"%8.2e ");
                }
                if( isRectangular )   
                    timing(parameters.dbase.get<int>("timeForAdvanceRectangularGrids"))+=getCPU()-time0;
                else
                    timing(parameters.dbase.get<int>("timeForAdvanceCurvilinearGrids"))+=getCPU()-time0;
                if( debug & 8 )
                {
                    display(u,sPrintF("u after advOpt and updateGhost, t=%8.2e",t),debugFile,"%8.2e ");
                }
                if( debug & 16 )
                {
                    if( addForcing ) 
                        ::display(f,sPrintF("  *** advIncompressibleSolidMechanics: Here is the forcing f grid=%i t=%9.3e ********\n",grid,t),
                                            pDebugFile,"%7.4f ");
                } 
            }

        
        } // end grid
        
        if( debug & 4 )
        {
            getErrors( next,t+dt,dt,sPrintF("\n *** advanceIncompressible: Errors after advance, before BC, t=%9.3e ******\n",t+dt) );
        }

        #undef FN

        
        gf[next].t=t+dt;

        if( correction==0 && orderOfAccuracyInSpace > 2  )
        {
      // --------------------------------------
      // ---- Extrapolate pressure in time ----
      // --------------------------------------

            if( t<= 2*dt )
                printF("Extrapolate p in time for BCs...\n"); 
            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
            {
                MappedGrid & mg = cg[grid];
                OV_GET_SERIAL_ARRAY(real,gf[prev].u[grid]    ,umLocal);
                OV_GET_SERIAL_ARRAY(real,gf[current].u[grid] , uLocal);
                OV_GET_SERIAL_ARRAY(real,gf[next].u[grid]    ,unLocal);
                getIndex(mg.dimension(),I1,I2,I3);

        // *** NOTE *** THIS ONLY NEEDS TO BE DONE NEAR CERTAIN BOUNDARIES ****************************
                if( orderOfTimeAccuracy==2 )
                {
          // unLocal(I1,I2,I3,pc) = 2.*uLocal(I1,I2,I3,pc) - umLocal(I1,I2,I3,pc);
          // assume unLocal holds p(t-2*dt)
                    unLocal(I1,I2,I3,pc) = 3.*uLocal(I1,I2,I3,pc) - 3.*umLocal(I1,I2,I3,pc) + unLocal(I1,I2,I3,pc);
                }
                else
                {
          // int numberOfPressureFunctions    = parameters.dbase.get<int>("numberOfPressureFunctions");
          // realCompositeGridFunction *& pgf = parameters.dbase.get<realCompositeGridFunction*>("pgf");  // holds p for time extrap

                    assert( numberOfGridFunctionsToUse==numberOfTimeLevels ); // this is assumed here I think

                    OV_GET_SERIAL_ARRAY(real,pgf[currentPressure][grid],pLocal); 
                    RealArray pSave(I1,I2,I3);
                    pSave(I1,I2,I3) = unLocal(I1,I2,I3,pc); // hold p(t-2*dt)

                    unLocal(I1,I2,I3,pc) = 4.*uLocal(I1,I2,I3,pc) - 6.*umLocal(I1,I2,I3,pc) + 4.*pSave(I1,I2,I3) - pLocal(I1,I2,I3);

          // Update saved value of p to the value at t-2*dt
                    assert( numberOfPressureFunctions==1 );
                    pLocal(I1,I2,I3) = pSave(I1,I2,I3); 

                }

        // OV_ABORT("finish me");     
            }
        }


        bool applyBC=true; // ** FIX ME **
            if( applyBC )
            {
                if( true ||   // *wdh* 091205 -- interpolate will call periodicUpdate and updateGhost 
                        cg.numberOfComponentGrids()>1 )
                {
                    real timei=getCPU();
                    if( debug & 4 )
                        gf[next].u.display(sPrintF("Cgsm::advanceIncompressible: gf[next].u before interpolate, t=%8.2e",t+dt),debugFile,"%8.2e ");
          // --- Note: interpolate performs a periodicUpdate and updateGhostBoundaries even if there is only one grid
                    gf[next].u.interpolate();
                    if( debug & 4 )
                        gf[next].u.display(sPrintF("Cgsm::advanceIncompressible: gf[next].u after interpolate, t=%8.2e",t+dt),debugFile,"%8.2e ");
                    if( debug & 4 )
                    {
                        getErrors( next,t+dt,dt,sPrintF("\n ************** advance Errors after interpolate t=%9.3e ******\n",t+dt));
                    }
                    timing(parameters.dbase.get<int>("timeForInterpolate"))+=getCPU()-timei;
                }
        // ============= Boundary Conditions =============
                int option=0; // not used.
                applyBoundaryConditions( option, dt, next,current ); // apply BC to "next" (current=previous time step)
                if( debug & 8 )  // & 64
                {
                    gf[next].u.display(sPrintF("Cgsm::advanceIncompressible: gf[next].u after applyBC, t=%8.2e",gf[next].t),debugFile,"%8.2e ");
                }
                if( debug & 4 )
                  {
                      getErrors( next,t+dt,dt,sPrintF("\n ************** advanceIncompressible Errors after applyBC t=%9.3e ******\n",t+dt));
                  }    
        // ---- assign values at material interfaces ------
        // *** this does nothing currently ***
                assignInterfaceBoundaryConditions( next, t, dt );  // is this the right place to do this?
            }


    // ----------------- ADD UPWIND DISSIPATION -----------------
        if( upwindSOS && globalStepNumber>=1 )
        {  
      // --------------- START LOOP OVER GRIDS -------------------
            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
            {
                real time0=getCPU();

                MappedGrid & mg = cg[grid];
                MappedGridOperators & mgop = (*cgop)[grid];

                getIndex(mg.gridIndexRange(),I1,I2,I3);
            
        // getBoundsForPML( mg,Iv );

                realMappedGridFunction & fieldPrev    =gf[prev].u[grid];
                realMappedGridFunction & fieldCurrent =gf[current].u[grid];
                realMappedGridFunction & fieldNext    =gf[next].u[grid];


                realArray & um = fieldPrev;
                realArray & u  = fieldCurrent;
                realArray & un = fieldNext;

                lambda = lambdaGrid(grid);
                mu = muGrid(grid);
                c1=(mu+lambda)/rho, c2= mu/rho;
                
                if( numberOfStepsTaken<1 ) 
                    printF(" advanceIncompressible:INFO mu=%8.2e for grid=%i (%s) \n",mu,grid,(const char*)cg[grid].getName());
                

                bool useOpt=true; // true;
                const bool isRectangular=mg.isRectangular();

        // ----- Add upwind dissipation ----

                OV_GET_SERIAL_ARRAY(real,um,umLocal);
                OV_GET_SERIAL_ARRAY(real,u ,uLocal);
                OV_GET_SERIAL_ARRAY(real,un,unLocal);
      
                {
                    real timeAdv=getCPU();
                    int gridType = isRectangular? 0 : 1;
                    const real adc=artificialDissipation; // do not scale *wdh* 090216
                    real *umptr=umLocal.getDataPointer();
                    real *uptr =uLocal.getDataPointer();
                    real *unptr=unLocal.getDataPointer();
                    real *vtptr = uptr; 
                    real *vnptr = uptr; 
                    real *vptr  = uptr; 
                    real *v1ptr = uptr; 
                    real *v2ptr = uptr; 
                    real *v3ptr = uptr; 
                    real *v4ptr = uptr; 
                    real *vt1ptr = uptr; 
                    real *vt2ptr = uptr; 
                    real *vt3ptr = uptr; 
                    real *vt4ptr = uptr; 
                        const bool addForcing =false;
                        real *fptr  =  uptr; // not used for upwind
          // assert( !useVariableDissipation || variableDissipation!=NULL );
          // real *pVarDis = useVariableDissipation ? varDis.getDataPointer() : uptr;
                    const intArray & mask = mg.mask();
                    OV_GET_SERIAL_ARRAY_CONST(int,mask,maskLocal)
          // #ifdef USE_PPP
          //     intSerialArray maskLocal;  getLocalArrayWithGhostBoundaries(mask,maskLocal);
          // #else
          //     const intSerialArray & maskLocal = mask; 
          // #endif
                    real *rxptr;
                    if( isRectangular )
                    {
                        rxptr=uptr;
                    }
                    else
                    {
                    #ifdef USE_PPP
                                realSerialArray rxLocal; getLocalArrayWithGhostBoundaries(mg.inverseVertexDerivative(),rxLocal);
                    #else
                                const realSerialArray & rxLocal=mg.inverseVertexDerivative();
                    #endif
                        rxptr = rxLocal.getDataPointer();
                    }
                    const bool centerNeeded=forcingOption==twilightZoneForcing || (forcingOption==planeWaveBoundaryForcing); // **************** fix this 
                    #ifdef USE_PPP
                        realSerialArray xy;
                        if( centerNeeded ) getLocalArrayWithGhostBoundaries(mg.center(),xy);
                    #else
                        const realSerialArray & xy = centerNeeded ? mg.center() : umLocal;
                    #endif
                    real *xyptr=xy.getDataPointer();
                    int maskNull[1];
                    int *maskptr = maskLocal.getDataPointer();
          // realSerialArray *dis = NULL;
          // real *pdis=uptr;
          // Macro to extract the pointers to the variable material property arrays
           // --- Variable material properies ---
                      GridMaterialProperties::MaterialFormatEnum materialFormat = GridMaterialProperties::constantMaterialProperties;
                      int ndMatProp=1;  // for piecewise constant materials, this is the leading dimension of the matVal array
                      int *matIndexPtr=maskptr;  // if not used, point to mask
                      real*matValPtr=uptr;       // if not used, point to u
                      if( parameters.dbase.get<int>("variableMaterialPropertiesOption")!=0 )
                      {
             // Material properties do vary 
                          std::vector<GridMaterialProperties> & materialProperties = 
                                    parameters.dbase.get<std::vector<GridMaterialProperties> >("materialProperties");
                          GridMaterialProperties & matProp = materialProperties[grid];
                          materialFormat = matProp.getMaterialFormat();
                          if( materialFormat==GridMaterialProperties::piecewiseConstantMaterialProperties )
                          {
                                    IntegerArray & matIndex = matProp.getMaterialIndexArray();
                              matIndexPtr = matIndex.getDataPointer();
                          }
                          RealArray & matVal = matProp.getMaterialValuesArray();
                          matValPtr = matVal.getDataPointer();
                          ndMatProp = matVal.getLength(0);  
             // ::display(matVal,"matVal");
                      }
                    bool twilightZone = forcingOption==twilightZoneForcing; 
                    int combineDissipationWithAdvance=0;  // not used anymore
                    real tForce = t;   // evaluate forcing at this time 
                        int option=1; //  0=update solution, 1=add upwind dissipation
                        if( false )
                            printF(">>>> advanceIncompressible: add upwind dissipation at globalStepNumber=%d, t=%9.3e <<<<<\n",globalStepNumber,t);
                    int ipar[]={option,
                                            gridType,
                                            orderOfAccuracyInSpace,
                                            orderOfTimeAccuracy,
                                            (int)addForcing,
                                            (int)twilightZone,
                                            uc,
                                            vc,
                                            wc,
                                            useWhereMask,
                                            (int)timeSteppingMethodSm,
                                            (int)useVariableDissipation,
                                            (int)useConservative,           // 12 
                                            combineDissipationWithAdvance,
                                            debug,
                                            computeUt,
                                            materialFormat,                 // 16 
                                            myid,
                                            pc,                             // 18 
                                            upwindSOS,
                                            correction
                                            };    
                    real dx[3]={1.,1.,1.};
                    if( isRectangular )
                        mg.getDeltaX(dx);
                    real rpar[30];
                    rpar[ 0]=dt;
                    rpar[ 1]=dx[0];
                    rpar[ 2]=dx[1];
                    rpar[ 3]=dx[2];
                    rpar[ 4]=adc;
                    rpar[ 5]=mg.gridSpacing(0);
                    rpar[ 6]=mg.gridSpacing(1);
                    rpar[ 7]=mg.gridSpacing(2);
                    rpar[ 8]=rho;
                    rpar[ 9]=mu;
                    rpar[10]=kx; // for plane wave scattering
                    rpar[11]=ky;
                    rpar[12]=kz;
                    rpar[13]=(real &)parameters.dbase.get<OGFunction* >("exactSolution");  // twilight zone pointer, ep
                    rpar[14]=tForce;
                    rpar[15]=parameters.dbase.get<real>("dtOld");  // dt used on the previous step
                    rpar[16]=rho;
                    rpar[17]=mu;
                    rpar[18]=lambda;
                    rpar[20]=0.;  // return cpu for dissipation
          // printF("** AdvanceIncompressible: gridType=%i, isRectangular=%i\n",gridType,(int)isRectangular);
                    bool ok = ParallelUtility::getLocalArrayBounds(u,uLocal,I1,I2,I3);
                    int ierr=0;
                    if( ok )
                    {
            // if( adc>0. && !combineDissipationWithAdvance )
            // {
            //   // create a temp array to hold the artificial dissipation
            //   dis = new realSerialArray(uLocal.dimension(0),uLocal.dimension(1),uLocal.dimension(2),uLocal.dimension(3));
            //   pdis = dis->getDataPointer();
            //   assert( pdis!=NULL );
            //   sizeOfLocalArraysForAdvance=max(sizeOfLocalArraysForAdvance,(double)(dis->elementCount()*sizeof(real)));
            // }
            // real timeAdv=getCPU();
            // if( combineDissipationWithAdvance )
            // {
            //   assert( umptr!=unptr );
            // }
                        advIsm(mg.numberOfDimensions(),
                                      I1.getBase(),I1.getBound(),I2.getBase(),I2.getBound(),I3.getBase(),I3.getBound(),
                                      uLocal.getBase(0),uLocal.getBound(0),
                                      uLocal.getBase(1),uLocal.getBound(1),
                                      uLocal.getBase(2),uLocal.getBound(2),
                                      uLocal.getBase(3),uLocal.getBound(3),
                                      *maskptr,*rxptr,*xyptr,  
                                      *umptr,*uptr,*unptr, *fptr, 
                                      *vtptr, *vnptr, *vptr,
                                      *v1ptr, *v2ptr, *v3ptr, *v4ptr,
                                      *vt1ptr,*vt2ptr,*vt3ptr,*vt4ptr,
                                      mg.gridIndexRange(0,0), mg.boundaryCondition(0,0),  ipar[0], rpar[0], ierr );
                    }
                    timeAdv=getCPU()-timeAdv;
                        timing(parameters.dbase.get<int>("timeForUpwindDissipation"))+=timeAdv;
                    if( debug & 8 )
                    {
                        display(unLocal,sPrintF("unLocal after advIncompressibleSolidMechanics, processor=%i before BC's t=%8.2e",
                                                                        Communication_Manager::My_Process_Number,t),pDebugFile,"%8.2e ");
                        display(un,sPrintF("un after advIncompressibleSolidMechanics, before BC's t=%8.2e",t),debugFile,"%8.2e ");
                    }
                    if( isRectangular )   
                        timing(parameters.dbase.get<int>("timeForAdvanceRectangularGrids"))+=getCPU()-time0;
                    else
                        timing(parameters.dbase.get<int>("timeForAdvanceCurvilinearGrids"))+=getCPU()-time0;
                    if( debug & 8 )
                    {
                        display(u,sPrintF("u after advOpt and updateGhost, t=%8.2e",t),debugFile,"%8.2e ");
                    }
                }


            } // end grid

                if( applyBC )
                {
                    if( true ||   // *wdh* 091205 -- interpolate will call periodicUpdate and updateGhost 
                            cg.numberOfComponentGrids()>1 )
                    {
                        real timei=getCPU();
                        if( debug & 4 )
                            gf[next].u.display(sPrintF("Cgsm::advanceIncompressible: gf[next].u before interpolate, t=%8.2e",t+dt),debugFile,"%8.2e ");
            // --- Note: interpolate performs a periodicUpdate and updateGhostBoundaries even if there is only one grid
                        gf[next].u.interpolate();
                        if( debug & 4 )
                            gf[next].u.display(sPrintF("Cgsm::advanceIncompressible: gf[next].u after interpolate, t=%8.2e",t+dt),debugFile,"%8.2e ");
                        if( debug & 4 )
                        {
                            getErrors( next,t+dt,dt,sPrintF("\n ************** advance Errors after interpolate t=%9.3e ******\n",t+dt));
                        }
                        timing(parameters.dbase.get<int>("timeForInterpolate"))+=getCPU()-timei;
                    }
          // ============= Boundary Conditions =============
                    int option=0; // not used.
                    applyBoundaryConditions( option, dt, next,current ); // apply BC to "next" (current=previous time step)
                    if( debug & 8 )  // & 64
                    {
                        gf[next].u.display(sPrintF("Cgsm::advanceIncompressible: gf[next].u after applyBC, t=%8.2e",gf[next].t),debugFile,"%8.2e ");
                    }
                    if( debug & 4 )
                      {
                          getErrors( next,t+dt,dt,sPrintF("\n ************** advanceIncompressible Errors after applyBC t=%9.3e ******\n",t+dt));
                      }    
          // ---- assign values at material interfaces ------
          // *** this does nothing currently ***
                    assignInterfaceBoundaryConditions( next, t, dt );  // is this the right place to do this?
                }
        } // end if upwindSOS


    // -----------------------------------------------------
    // ---------------- PRESSURE SOLVE ---------------------
    // -----------------------------------------------------
        if( correction<=(numberOfCorrections-1) || !skipLastPressureSolve)
        {
            solveForPressure( current, t+dt, dt );

            if( orderOfAccuracyInSpace==4 )
            { // re-apply BC's that depend on p 
                    if( applyBC )
                    {
                        if( true ||   // *wdh* 091205 -- interpolate will call periodicUpdate and updateGhost 
                                cg.numberOfComponentGrids()>1 )
                        {
                            real timei=getCPU();
                            if( debug & 4 )
                                gf[next].u.display(sPrintF("Cgsm::advanceIncompressible: gf[next].u before interpolate, t=%8.2e",t+dt),debugFile,"%8.2e ");
              // --- Note: interpolate performs a periodicUpdate and updateGhostBoundaries even if there is only one grid
                            gf[next].u.interpolate();
                            if( debug & 4 )
                                gf[next].u.display(sPrintF("Cgsm::advanceIncompressible: gf[next].u after interpolate, t=%8.2e",t+dt),debugFile,"%8.2e ");
                            if( debug & 4 )
                            {
                                getErrors( next,t+dt,dt,sPrintF("\n ************** advance Errors after interpolate t=%9.3e ******\n",t+dt));
                            }
                            timing(parameters.dbase.get<int>("timeForInterpolate"))+=getCPU()-timei;
                        }
            // ============= Boundary Conditions =============
                        int option=0; // not used.
                        applyBoundaryConditions( option, dt, next,current ); // apply BC to "next" (current=previous time step)
                        if( debug & 8 )  // & 64
                        {
                            gf[next].u.display(sPrintF("Cgsm::advanceIncompressible: gf[next].u after applyBC, t=%8.2e",gf[next].t),debugFile,"%8.2e ");
                        }
                        if( debug & 4 )
                          {
                              getErrors( next,t+dt,dt,sPrintF("\n ************** advanceIncompressible Errors after applyBC t=%9.3e ******\n",t+dt));
                          }    
            // ---- assign values at material interfaces ------
            // *** this does nothing currently ***
                        assignInterfaceBoundaryConditions( next, t, dt );  // is this the right place to do this?
                    }
            }
        }
        else
        {
            if( t<=2*dt )
            {
                  printF("\n\n &&&&&&&&&& SKIP LAST PRESSURE SOLVE correction=%d numberOfCorrections=%d &&&&&&&&&&\n\n",correction,numberOfCorrections);
            }
        }

    } // end correction

    currentVelocity = nextVelocity; // update velocity index 
    currentPressure = nextPressure; // update pressure index 

    parameters.dbase.get<real>("dtOld")=dt;  // set dtOld   

    checkArrays("advanceIncompressibleMethodOfLines:end");
    
}







// // =============================================================================
// //   *** OLD VERSION AUg 26, 2021 --> upwind with no BC's afterward
// ///
// /// \brief Advance the solution for
// ///       Incompressible Linear Elastic Equations
// ///
// // \param ut : if not NULL, compute u.tt and return in this grid-function. Otherwise
// //             return the solution gf[next] at t=t+dt.
// // =============================================================================
// void Cgsm::
// advanceIncompressibleOld( int current, real t, real dt,
//                           RealCompositeGridFunction *ut /* = NULL */, 
//                          real tForce /* = 0. */ )
// {   

//   int & globalStepNumber = parameters.dbase.get<int >("globalStepNumber");
//   if( debug() & 4 )
//     printF("\n @@@@@@@ Cgsm::advanceIncompressible called at globalStepNumber=%d, t=%9.3e\n",globalStepNumber,t);

//   checkArrays("advanceIncompressible:start");


//   FILE *& debugFile  =parameters.dbase.get<FILE* >("debugFile");
//   FILE *& logFile    =parameters.dbase.get<FILE* >("logFile");
//   FILE *& pDebugFile =parameters.dbase.get<FILE* >("pDebugFile");
    
//   const int numberOfDimensions       = cg.numberOfDimensions();
//   const int numberOfComponentGrids   = cg.numberOfComponentGrids();
//   const int & numberOfComponents     = parameters.dbase.get<int >("numberOfComponents");
//   const int & uc                     = parameters.dbase.get<int >("uc");
//   const int & vc                     = parameters.dbase.get<int >("vc");
//   const int & wc                     = parameters.dbase.get<int >("wc");
//   const int & rc                     = parameters.dbase.get<int >("rc");
//   const int & tc                     = parameters.dbase.get<int >("tc");
//   const int & pc                     = parameters.dbase.get<int >("pc");
//   const int & orderOfAccuracyInSpace = parameters.dbase.get<int>("orderOfAccuracy");
//   const int & orderOfTimeAccuracy    = parameters.dbase.get<int>("orderOfTimeAccuracy");

//   SmParameters::TimeSteppingMethodSm & timeSteppingMethodSm = 
//                                    parameters.dbase.get<SmParameters::TimeSteppingMethodSm>("timeSteppingMethodSm");
//   RealArray & timing = parameters.dbase.get<RealArray >("timing");

//   Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
//   Range C=numberOfComponents;
//   const int prev = (current-1+numberOfTimeLevels) % numberOfTimeLevels;
//   const int next = (current+1)                    % numberOfTimeLevels;

//   real & rho=parameters.dbase.get<real>("rho");
//   real & mu = parameters.dbase.get<real>("mu");
//   real & lambda = parameters.dbase.get<real>("lambda");
//   RealArray & muGrid = parameters.dbase.get<RealArray>("muGrid");
//   RealArray & lambdaGrid = parameters.dbase.get<RealArray>("lambdaGrid");
//   bool & gridHasMaterialInterfaces = parameters.dbase.get<bool>("gridHasMaterialInterfaces");
//   int & debug = parameters.dbase.get<int >("debug");

//   const int & upwindSOS = parameters.dbase.get<int>("upwindSOS"); 

//   const real cMax=max(lambdaGrid+muGrid)/rho;

//   const int computeUt = ut != NULL;  // compute u.tt


//   const bool extrapInterpolationNeighbours = upwindSOS;
//   if( extrapInterpolationNeighbours )
//   {
//     // -- Extrapolate interpolation neighbours for the artificial dissipation ---
//     // This may also be done in advanceSOS
//     if( debug & 4 )
//       printF("advanceIncompressible: Extrapolate interpolation neighbours for upwind dissipation, t=%9.3e\n",t);

//     extrapolateInterpolationNeighbours( gf[current], C );
//   }
    
//    // // In some cases we combine the artificial dissipation loop with the main loop
//    //  int combineDissipationWithAdvance = isRectangular && 
//    //    !useVariableDissipation &&
//    //    timeSteppingMethodSm==SmParameters::modifiedEquationTimeStepping &&
//    //    orderOfAccuracyInSpace==4 && orderOfTimeAccuracyInTime==4;

//   const int useWhereMask = numberOfComponentGrids>1;  

//   sizeOfLocalArraysForAdvance=0.;
//   // --------------- START LOOP OVER GRIDS -------------------
//   for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
//   {
//     real time0=getCPU();

//     MappedGrid & mg = cg[grid];
//     MappedGridOperators & mgop = (*cgop)[grid];

//     getIndex(mg.gridIndexRange(),I1,I2,I3);
    
//     getBoundsForPML( mg,Iv );

//     realMappedGridFunction & fieldPrev    =gf[prev].u[grid];
//     realMappedGridFunction & fieldCurrent =gf[current].u[grid];
//     realMappedGridFunction & fieldNext    =gf[next].u[grid];


//     realArray & um = fieldPrev;
//     realArray & u  = fieldCurrent;
//     realArray & un = fieldNext;

//     if( debug & 4 )
//     {
//       fPrintF(debugFile," **** start of advance, t=%8.2e\n",t);
//       fprintf(pDebugFile," **** start of advance, t=%8.2e\n",t);
            
//       if( debug & 8 )
//       {
//         display(um,sPrintF("um start of advance, t=%8.2e",t),debugFile,"%8.2e ");
//         display(u,sPrintF("u start of advance, t=%8.2e",t),debugFile,"%8.2e ");
//       }
//     }

//     lambda = lambdaGrid(grid);
//     mu = muGrid(grid);
//     c1=(mu+lambda)/rho, c2= mu/rho;
        
//     if( numberOfStepsTaken<1 ) 
//       printF(" advanceIncompressible:INFO mu=%8.2e for grid=%i (%s) \n",mu,grid,(const char*)cg[grid].getName());
        
//     // const real dtsq=dt*dt; 
//     // const real adc=artificialDissipation*SQR(cMax); // scale dissipation by c^2 *wdh* 041103


//     bool useOpt=true; // true;
//     const bool isRectangular=mg.isRectangular();
//     if( !isRectangular )
//     {
//        // *** PUT THIS IN SETUP SOMEWHERE ***** 
//       real timea=getCPU();
//       mg.update( MappedGrid::THEinverseVertexDerivative | MappedGrid::THEinverseCenterDerivative );
//       timea=getCPU()-timea;
//       timing(parameters.dbase.get<int>("timeForInitialize"))+=timea;

//       time0+=timea;  // do not include with time for curvilinear
//     }

//     //  ::display(lap,"lap","%8.1e ");

//     // --- get any forcing ----
//     getForcingMacro(); 

//     // We new 3 levels of solution to apply the upwinding
//     // t=0 -> dt    : globalStepNumber=0
//     // t=dt -> 2*dt : globalStepNumber=1
//     if( upwindSOS && globalStepNumber>=1 )
//     {
//       // ----- Add upwind dissipation ----
//       int numberOfTimeLevelsStored = numberOfTimeLevels; // **CHECK ME**
//       const int prev2= (current-2+numberOfTimeLevelsStored) % numberOfTimeLevelsStored;

//       // -- rename these ---
//       realMappedGridFunction & fieldPrev    =gf[prev2 ].u[grid];
//       realMappedGridFunction & fieldCurrent =gf[prev  ].u[grid];
//       realMappedGridFunction & fieldNext    =gf[current].u[grid];


//       realArray & um = fieldPrev;
//       realArray & u  = fieldCurrent;
//       realArray & un = fieldNext;      
    
//       OV_GET_SERIAL_ARRAY(real,gf[prev2  ].u[grid],umLocal);  // u(t-2*dt)
//       OV_GET_SERIAL_ARRAY(real,gf[prev   ].u[grid],uLocal);   // u(t-dt)
//       OV_GET_SERIAL_ARRAY(real,gf[current].u[grid],unLocal);  // u(t)

//       if( false )
//         printf("\n >>>>Upwind at step=%d, t=%9.3e: prev2=%d (t=%9.3e) prev=%d (t=%9.3e) current=%d (t=%9.3e) (next=%d)\n",
//                globalStepNumber,t, prev2,gf[prev2].t, prev,gf[prev].t, current,gf[current].t,next);
//       advanceSolutionMacro(UPWIND);

//       fieldNext.periodicUpdate();          

//     }

//     OV_GET_SERIAL_ARRAY(real,um,umLocal);
//     OV_GET_SERIAL_ARRAY(real,u ,uLocal);
//     OV_GET_SERIAL_ARRAY(real,un,unLocal);
  

//     // ------ advance the interior points -------
//     advanceSolutionMacro(UPDATE);

    
//   } // end grid
    
//   if( debug & 4 )
//   {
//     getErrors( next,t+dt,dt,sPrintF("\n *** advanceIncompressible: Errors after advance, before BC, t=%9.3e ******\n",t+dt) );
//   }

// #undef FN

    
//   gf[next].t=t+dt;


//   bool applyBC=true; // ** FIX ME **
//   interpolateAndApplyBoundaryConditionsMacro();



//   // -----------------------------------------------------
//   // ---------------- PRESSURE SOLVE ---------------------
//   // -----------------------------------------------------
//   solveForPressure( current, t+dt, dt );

//   parameters.dbase.get<real>("dtOld")=dt;  // set dtOld 

//   checkArrays("advanceIncompressible:end");
    
// }
