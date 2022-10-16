// This file automatically generated from applyBoundaryConditions.bC with bpp.
#include "Cgad.h"
#include "AdParameters.h"
#include "App.h"
#include "ParallelUtility.h"
#include "ParallelGridUtility.h"
#include "Interface.h"  

// ===================================================================================
//   This macro extracts the boundary data arrays
//
//  *wdh* 110312 THIS WAS COPIED FROM cg/sm/src -- FIX ME ---
// ===================================================================================


#define ForBoundary(side,axis)   for( int axis=0; axis<numberOfDimensions; axis++ ) for( int side=0; side<=1; side++ )

#define bcOptAd EXTERN_C_NAME(bcoptad)
extern "C"
{

void bcOptAd( const int&nd, 
                            const int&nd1a,const int&nd1b,const int&nd2a,const int&nd2b,const int&nd3a,const int&nd3b,
                            const int&gridIndexRange, const int& dimRange, const int &isPeriodic, real&u, const int&mask,
                            const real&rsxy, const real&xy, const int&boundaryCondition, 
                            const int&addBoundaryForcing, const int&interfaceType, const int&ndb, const real&bcData,  
                            const int&dim, const real&bcf0, const int64_t & bcfOffset, 
                            const int&ipar, const real&rpar, int&ierr );

}

//    Mixed-derivative BC for component i: 
//          mixedCoeff(i)*u(i) + mixedNormalCoeff(i)*u_n(i) = mixedRHS(i)
#define mixedRHS(component,side,axis,grid)         bcData(component+numberOfComponents*(0),side,axis,grid)
#define mixedCoeff(component,side,axis,grid)       bcData(component+numberOfComponents*(1),side,axis,grid)
#define mixedNormalCoeff(component,side,axis,grid) bcData(component+numberOfComponents*(2),side,axis,grid)

// =======================================================================================
// Macro: get the parameters for calling the optimized fortran BC routine
// =======================================================================================

int Cgad::
applyBoundaryConditions(const real & t, realMappedGridFunction & u, 
                                                realMappedGridFunction & gridVelocity,
                                                const int & grid,
                                                const int & option /* =-1 */,
                                                realMappedGridFunction *puOld /* =NULL */,  
                                                realMappedGridFunction *pGridVelocityOld /* =NULL */,
                                                const real & dt /* =-1. */ )
//=========================================================================================
// /Description:
//   Apply boundary conditions for the advection diffusion equations
// 
// /t (input):
// /u (input/output) : apply to this grid function.
// /gridIsMoving (input) : true if this grid is moving.
// /gridVelocity (input) : the grid velocity if gridIsMoving==true.
// /variableBoundaryData (input) : true if there is boundary data that depends on the position along the boundary.
// /boundaryData (input) : boundary data used if variableBoundaryData==true.
// grid (input) : the grid number if this MappedGridFunction is part of a CompositeGridFunction.
// option (input): not used here.
//
// /Note:
// ***Remember to also change the BC routine for implicit time stepping if changes are made here
// applyBoundaryConditionsForImplicitTimeStepping
//
// 
//\end{MappedGridSolverInclude.tex}  
//=========================================================================================
{
    real time0=getCPU();

  // WE SHOULD BE ABLE TO JUST RETURN FOR IMPLICIT GRIDS  -- *wdh* April 24, 2020
    if( parameters.getGridIsImplicit(grid) )
        return 0;


    const int & myid = parameters.dbase.get<int>("myid");
    const int np= max(1,Communication_Manager::numberOfProcessors());

    const int & multiDomainProblem = parameters.dbase.get<int>("multiDomainProblem"); 
    const bool twilightZoneFlow = parameters.dbase.get<bool>("twilightZoneFlow");
    const int & gridIsImplicit = parameters.getGridIsImplicit(grid);
    const int applyChampInterfaceConditions = parameters.dbase.get<int>("applyChampInterfaceConditions");


    const int orderOfAccuracy=parameters.dbase.get<int >("orderOfAccuracy");
    assert( orderOfAccuracy==2 || orderOfAccuracy==4 );

    BoundaryConditionParameters bcParams;
    const int numberOfComponents=parameters.dbase.get<int >("numberOfComponents");
    Range C(0,numberOfComponents-1);

    const RealArray & bcData = parameters.dbase.get<RealArray>("bcData");
    BoundaryData::BoundaryDataArray & pBoundaryData = parameters.getBoundaryData(grid);

    const bool & gridIsMoving = parameters.gridIsMoving(grid);

    MappedGrid & mg = *u.getMappedGrid();
    const int numberOfDimensions = mg.numberOfDimensions();
    const bool isRectangular=mg.isRectangular();

  // -- evaluate the known solution where needed ----
    const Parameters::KnownSolutionsEnum & knownSolution = parameters.dbase.get<Parameters::KnownSolutionsEnum >("knownSolution");
    const bool & assignKnownSolutionAtBoundaries = parameters.dbase.get<bool>("assignKnownSolutionAtBoundaries");

    if( knownSolution == Parameters::userDefinedKnownSolution && assignKnownSolutionAtBoundaries )
    {
    // Fill in Dirichlet BCs with the known solution
        OV_GET_SERIAL_ARRAY(real,u,uLocal);
        RealArray ua; 

        Index Ib1,Ib2,Ib3;
        ForBoundary(side,axis)
        {
            if( mg.boundaryCondition(side,axis) == AdParameters::dirichletBoundaryCondition )
            {
                printF("))))))))) CGAD:applyExplicitBC:  Set Known solution at dirichlet BC (side,axis,grid)=(%dm%dm%d) t=%9.3e  (((((((((\n",side,axis,grid,t);

                getBoundaryIndex(mg.gridIndexRange(),side,axis,Ib1,Ib2,Ib3);
                ua = parameters.getKnownSolution(t, grid, Ib1,Ib2,Ib3 );  
                uLocal(Ib1,Ib2,Ib3,C) = ua(Ib1,Ib2,Ib3,C); 


        //  getUserDefinedKnownSolution(real t, CompositeGrid & cg, int grid, RealArray & ua, 
        // const Index & I1, const Index &I2, const Index &I3, int numberOfTimeDerivatives = 0 );        

            }

        }

    }


  // New optimized version -- currently for explicit TS, and 1 component only  *wdh* April 17, 2021
    bool useOpt= multiDomainProblem==0 && numberOfComponents==1 && !(parameters.getGridIsImplicit(grid));  // **FIX ME**
    
    if( t<=2*dt )
        printF("CGAD: applyBC: useOpt=%d, numberOfComponents=%d\n",(int)useOpt,numberOfComponents);
    if( useOpt )
    {
    // ----- Optimized Boundary Conditions ------
    // *wdh* April 17, 2021
        const IntegerArray & interfaceType = parameters.dbase.get<IntegerArray >("interfaceType");

        OV_GET_SERIAL_ARRAY(real,u,uLocal);
        OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);

    // get parameters for calling fortran
            IntegerArray indexRangeLocal(2,3), dimLocal(2,3), bcLocal(2,3);
            ParallelGridUtility::getLocalIndexBoundsAndBoundaryConditions( u,indexRangeLocal,dimLocal,bcLocal );
            real dx[3]={1.,1.,1.};
            if( isRectangular )
                mg.getDeltaX(dx);
            int gridType = isRectangular ? 0 : 1;
            int tc=0; // first component
            int ipar[] = {
                tc                  ,            // ipar( 0)
                numberOfComponents  ,            // ipar( 1)
                grid                ,            // ipar( 2)
                gridType            ,            // ipar( 3)
                orderOfAccuracy     ,            // ipar( 4)
                gridIsMoving        ,            // ipar( 5)
                gridIsImplicit      ,            // ipar( 6)
                twilightZoneFlow    ,            // ipar( 7)
                np                  ,            // ipar( 8)
                debug()             ,            // ipar( 9)
                myid                ,            // ipar(10)
                applyChampInterfaceConditions,   // ipar(11)
                assignKnownSolutionAtBoundaries  // ipar(12)
                                      };
            real rpar[] = {
                t                , //  rpar( 0)
                dt               , //  rpar( 1)
                dx[0]            , //  rpar( 2)
                dx[1]            , //  rpar( 3)
                dx[2]            , //  rpar( 4)
                mg.gridSpacing(0), //  rpar( 5)
                mg.gridSpacing(1), //  rpar( 6)
                mg.gridSpacing(2), //  rpar( 7)
                (real &)(parameters.dbase.get<OGFunction* >("exactSolution")) ,        //  rpar( 8) ! pointer for exact solution -- new : 110311 
                REAL_MIN           //  rpar( 9)
                                        };
      // The next macro extracts the boundary data arrays for use in Fortran
                int pdbc[2*3*2*3];
                #define dbc(s,a,side,axis) (pdbc[(s)+2*((a)+3*((side)+2*(axis)))])
                int pAddBoundaryForcing[6];
                #define addBoundaryForcing(side,axis) (pAddBoundaryForcing[(side)+2*(axis)])
                real *pbcf[2][3];
        // long int pbcfOffset[6];
        // We need an 8 byte integer so we can pass to fortran: int64_t is in stdint.h 
                int64_t pbcfOffset[6];
                #define bcfOffset(side,axis) pbcfOffset[(side)+2*(axis)]
                for( int axis=0; axis<=2; axis++ )
                {
                    for( int side=0; side<=1; side++ )
                    {
            // *** for now make sure the boundary data array is allocated on all sides
                        if( false &&   // We do NOT need to always allocate the boundaryDataArray for INS *wdh* 110313
                                ( pBoundaryData[side][axis]==NULL || parameters.isAdaptiveGridProblem() ) && 
                                mg.boundaryCondition(side,axis)>0 )
                        {
                            parameters.getBoundaryData(side,axis,grid,mg);
              // RealArray & bd = *pBoundaryData[side][axis]; // this is now done in the above line *wdh* 090819
              // bd=0.;
                        }
                        if( pBoundaryData[side][axis]!=NULL )
                        {
                            addBoundaryForcing(side,axis)=true;
                            RealArray & bd = *pBoundaryData[side][axis];
                            pbcf[side][axis] = bd.getDataPointer();
              // if( debug & 8 )
              // ::display(bd,sPrintF(" ++++ Cgsm: Here is bd (side,axis)=(%i,%i) ++++",side,axis),"%4.2f ");
                            for( int a=0; a<=2; a++ )
                            {
                                dbc(0,a,side,axis)=bd.getBase(a);
                                dbc(1,a,side,axis)=bd.getBound(a);
                            }
                        }
                        else
                        {
                            addBoundaryForcing(side,axis)=false;
                            pbcf[side][axis] =bcData.getDataPointer();  // should not be used in this case 
                            for( int a=0; a<=2; a++ )
                            {
                                dbc(0,a,side,axis)=0;
                                dbc(1,a,side,axis)=0;
                            }
                        }
            // for now we save the offset in a 4 byte int (double check that this is ok)
                        int64_t offset = pbcf[side][axis]- pbcf[0][0];
      //       if( offset > INT_MAX )
      //       {
      //      printF("ERROR: offset=%li INT_MAX=%li \n",offset,(long int)INT_MAX);
      //       }
      //       assert( offset < INT_MAX );
                        bcfOffset(side,axis) = offset;
            // bcfOffset(side,axis) = pbcf[side][axis]- pbcf[0][0];
            // cout << " **** bcfOffset= " << bcfOffset(side,axis) << endl;
                    }
                }
            real *pu = uLocal.getDataPointer();
            int *pmask = maskLocal.getDataPointer();
            real temp, *pxy=&temp, *prsxy=&temp;
            if( !isRectangular )
            {
                #ifdef USE_PPP
                  prsxy=mg.inverseVertexDerivative().getLocalArray().getDataPointer();
                #else
                  prsxy=mg.inverseVertexDerivative().getDataPointer();
                #endif    
            }
            bool vertexNeeded = twilightZoneFlow;
            if( vertexNeeded )
            {
                #ifdef USE_PPP
                  pxy=mg.vertex().getLocalArray().getDataPointer();
                #else
                  pxy=mg.vertex().getDataPointer();
                #endif    
            }

        int ierr=0;
        bcOptAd(mg.numberOfDimensions(),
                      uLocal.getBase(0),uLocal.getBound(0),uLocal.getBase(1),uLocal.getBound(1),
                      uLocal.getBase(2),uLocal.getBound(2),
                      indexRangeLocal(0,0), dimLocal(0,0), mg.isPeriodic(0),
                      *pu, *pmask, *prsxy, *pxy,  bcLocal(0,0),  
                      *pAddBoundaryForcing, *interfaceType.getDataPointer(), 
                      bcData.getLength(0), *bcData.getDataPointer(), 
                      *pdbc, *pbcf[0][0], pbcfOffset[0], ipar[0],rpar[0], ierr );

    // ...swap periodic edges 
        u.periodicUpdate();
        u.updateGhostBoundaries();
        return 0;

    }

    if( parameters.getGridIsImplicit(grid) )
    {
    // *wdh* 080829 -- only apply some BC's for implicit time stepping grids

        if( !assignKnownSolutionAtBoundaries )
            u.applyBoundaryCondition(C,BCTypes::dirichlet,Parameters::dirichletBoundaryCondition,bcData,pBoundaryData,t,bcParams,grid); 

        u.applyBoundaryCondition(C,BCTypes::neumann,Parameters::neumannBoundaryCondition,bcData,pBoundaryData,t,bcParams,grid);

    // assigned extended boundaries: *wdh* 080829 -- treat the case when a dirichlet-BC is next to an interface
        if( false )
        {
            bcParams.extraInTangentialDirections= orderOfAccuracy==2 ? 1 : 2;
            u.applyBoundaryCondition(C,BCTypes::dirichlet,Parameters::dirichletBoundaryCondition,bcData,pBoundaryData,t,
                                                              bcParams,grid); 
            bcParams.extraInTangentialDirections=0;
        }
        
        return 0;
    }
    

    if( debug() & 8 )
    {
        printF(">>>>> Cgad::applyBoundaryConditions <<<<<<<\n");
    }
    checkArrayIDs("  Cgad::applyBoundaryConditions: start"); 

    const IntegerArray & interfaceType = parameters.dbase.get<IntegerArray >("interfaceType");


    for( int axis=0; axis<mg.numberOfDimensions(); axis++ ) 
    {
        for( int side=0; side<=1; side++ )
        {
            int bc = mg.boundaryCondition(side,axis);
            if( !(bc==AdParameters::dirichletBoundaryCondition ||
                        bc==AdParameters::neumannBoundaryCondition ||
                        bc==AdParameters::mixedBoundaryCondition ||
                        bc == AdParameters::interfaceBoundaryCondition ||
                        bc == AdParameters::axisymmetric ||
                        bc <=0 ) )
            {
                printF("Cgad:applyBoundaryConditions:ERROR: unexpected boundary condition bc=%i for (grid,side,axis)=(%i,%i,%i)\n",
                              bc,grid,side,axis);
                Overture::abort("error");
            }
      // --- check that interface boundaries use the correct bc ---
            if( parameters.dbase.get<int>("applyInterfaceBoundaryConditions")==1 &&
                    interfaceType(side,axis,grid)!=Parameters::noInterface &&
                    bc!=AdParameters::mixedBoundaryCondition )
            {
                printP("Cgad:applyBoundaryConditions:ERROR:the interface on (side,axis,grid)=(%i,%i,%i)\n"
                              " should have a mixed boundary condition associated with it, but bc=%i\n",
                              side,axis,grid,bc);
                Overture::abort("error");
            }
        }
    }
    


    u.applyBoundaryCondition(C,BCTypes::dirichlet,Parameters::dirichletBoundaryCondition,bcData,pBoundaryData,t,
                                                      bcParams,grid); 

    BoundaryConditionParameters extrapParams;
    extrapParams.orderOfExtrapolation=parameters.dbase.get<int >("orderOfAccuracy")+1;
  // *wdh* 071125 u.applyBoundaryCondition(C,BCTypes::extrapolate,BCTypes::allBoundaries,0.,t,extrapParams);
    u.applyBoundaryCondition(C,BCTypes::extrapolate,Parameters::dirichletBoundaryCondition,0.,t,extrapParams);


    u.applyBoundaryCondition(C,BCTypes::neumann,Parameters::neumannBoundaryCondition,bcData,pBoundaryData,t,bcParams,grid);

  // An interface could have a mixed-BC which is really a dirichlet BC -- we need to check this 
  // u.applyBoundaryCondition(C,BCTypes::neumann,Parameters::mixedBoundaryCondition,bcData,pBoundaryData,t,bcParams,grid);
    for( int axis=0; axis<mg.numberOfDimensions(); axis++ ) 
    {
        for( int side=0; side<=1; side++ )
        {
      // printP("applyBC: interfaceType(side=%i,axis=%i,grid=%i)=%i (noInterface=%i)\n",side,axis,grid,interfaceType(side,axis,grid),
      //             (int)Parameters::noInterface );
            
            if( mg.boundaryCondition(side,axis)==AdParameters::mixedBoundaryCondition )
            {
                if( interfaceType(side,axis,grid)!=Parameters::noInterface )
                { // This is an interface between domains

                    if( parameters.dbase.get<int>("applyInterfaceBoundaryConditions")==0 )
                    {
                        printP("Cgad:applyBC:skip interface bc: interfaceType(side=%i,axis=%i,grid=%i)=%i\n",
                                      side,axis,grid,interfaceType(side,axis,grid));
                        continue;
                    }

          // what about BC's applied at t=0 before the boundary data is set ??
          // if( parameters.dbase.get<int >("globalStepNumber") < 2 ) continue; // ********************* TEMP *****

          // if this is an interface we should turn off the TZ forcing for the boundary condition since we want
          // to use the boundary data instead.
                    if( debug() & 4 )
                        printP("Cgad:applyBC: apply a mixed BC on the interface (side,axis,grid)=(%i,%i,%i) %f*T + %f*T.n \n",
                                      side,axis,grid,mixedCoeff(0,side,axis,grid),mixedNormalCoeff(0,side,axis,grid));

          // RealArray *&pBoundaryData = boundaryData[grid].boundaryData[side][axis];
          // if( pBoundaryData!=NULL )
                    if( pBoundaryData[side][axis]!=NULL )
                        u.getOperators()->setTwilightZoneFlow( false );
                    else
                    {
                        printP("$$$$ Cgad:applyBC:INFO:interface: pBoundaryData is NULL for [side,axis,grid]=[%i,%i,%i], "
                                      "t=%9.3e.\n",
                                      side,axis,grid,t );
                    }
                    
                }
                for( int c=0; c<numberOfComponents; c++ )
                {
                    if( mixedNormalCoeff(c,side,axis,grid)!=0. ) // coeff of T.n is non-zero
                    {
                        real a0=mixedCoeff(c,side,axis,grid);
                        real a1=mixedNormalCoeff(c,side,axis,grid);
                        bcParams.a.redim(3);
                        bcParams.a(0)=a0;
                        bcParams.a(1)=a1;
                        bcParams.a(2)=0.;
                        u.applyBoundaryCondition(c,BCTypes::mixed,BCTypes::boundary(side,axis),bcData,pBoundaryData,
                                                                          t,bcParams,grid);
                    }
                    else
                    {
                        printP("Cgad:applyBC: apply a dirichlet BC on the interface (side,axis,grid)=(%i,%i,%i)\n",side,axis,grid);
                        real a0=mixedCoeff(c,side,axis,grid);
                        assert( a0==1. );
                        u.applyBoundaryCondition(c,BCTypes::dirichlet,BCTypes::boundary(side,axis),bcData,pBoundaryData,
                                                                          t,bcParams,grid);

            // do this for now:  -- we could do better --
                        u.applyBoundaryCondition(c,BCTypes::extrapolate,BCTypes::boundary(side,axis),0.,t,extrapParams);
                    }
                }
                if( interfaceType(side,axis,grid)!=Parameters::noInterface )
                { // reset TZ
                    u.getOperators()->setTwilightZoneFlow( parameters.dbase.get<bool >("twilightZoneFlow") );
                }
            }
        }
    }
    


    u.applyBoundaryCondition(C,BCTypes::neumann,Parameters::axisymmetric,0.,t);

  // --- Extrapolate the second ghost line for fourth-order -----
    if( orderOfAccuracy > 2 )
    {
    // **check me *** April 16, 2021
        extrapParams.ghostLineToAssign=2;
        u.applyBoundaryCondition(C,BCTypes::extrapolate,BCTypes::allBoundaries,0.,t,extrapParams);
    }

    u.finishBoundaryConditions();

    parameters.dbase.get<RealArray>("timing")(parameters.dbase.get<int>("timeForBoundaryConditions"))+=getCPU()-time0;

    return 0;
}

// ====================================================================================
// Macro: Add twilght-zone corrections
// ====================================================================================

// ====================================================================================
// Macro: Add twilght-zone corrections for CHAMP
// ====================================================================================

// // ====================================================================================
// // Macro: Add twilght-zone corrections for CHAMP
// // ====================================================================================
// #beginMacro addTwilightZoneCorrectionForChamp()              

//   OGFunction & e = *(parameters.dbase.get<OGFunction* >("exactSolution"));

//   const RealArray & champParameters = parameters.dbase.get<RealArray>("champParameters");
//   // const real Sl    = champParameters(0,side,axis,grid);
//   // const real theta = champParameters(2,side,axis,grid);
//   // const real beta  = champParameters(3,side,axis,grid);  
//   const Real pl    = champParameters(0,side,axis,grid);    // optimized Scwartz Parameter for side 1
//   const Real pr    = champParameters(1,side,axis,grid);    // optimized Scwartz Parameter for side 2
//   const Real theta = champParameters(2,side,axis,grid);    // K1/K2
//   const Real beta  = champParameters(3,side,axis,grid);    // D1/D2    
//   const Real Sl    = champParameters(4,side,axis,grid);    // pl/dxs; 

//   printP("**** ADD TZ Correction to results from getData for CHAMP:  Sl=%g, theta=%g, beta=%g *********************************\n",Sl,theta,beta);                  

//   OV_GET_VERTEX_BOUNDARY_NORMAL(mg,side,axis,normal);

//   RealArray ue(Ib1,Ib2,Ib3,N), uex(Ib1,Ib2,Ib3,N), uey(Ib1,Ib2,Ib3,N), uexx(Ib1,Ib2,Ib3,N), ueyy(Ib1,Ib2,Ib3,N);
//   int rectangular=0;
//   e.gd( ue  ,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,0,Ib1,Ib2,Ib3,N,t);
//   e.gd( uex ,xLocal,mg.numberOfDimensions(),rectangular,0,1,0,0,Ib1,Ib2,Ib3,N,t);
//   e.gd( uey ,xLocal,mg.numberOfDimensions(),rectangular,0,0,1,0,Ib1,Ib2,Ib3,N,t);
//   e.gd( uexx,xLocal,mg.numberOfDimensions(),rectangular,0,2,0,0,Ib1,Ib2,Ib3,N,t);
//   e.gd( ueyy,xLocal,mg.numberOfDimensions(),rectangular,0,0,2,0,Ib1,Ib2,Ib3,N,t);


//   if( isRectangular )
//   {

//     // ONLY VALID FOR CARTESIAN AND  axis==0 
//     assert( axis==0 );

//     const Real dxs = dx[axis]; 

//     const real a0 = Sl;
//     const real a1 = theta + Sl*dxs*theta;


//     const real a2 = dxs*( beta      ) + Sl*( .5*SQR(dxs)*beta      );
//     const real a3 = dxs*( (beta-1.) ) + Sl*( .5*SQR(dxs)*(beta-1.) );
//     const real a4 = a3; // for 3D 

//     int n=0;
//     assert( N.getLength()==1 );
//     interfaceData.u(Ib1,Ib2,Ib3,n) += a0*ue + a1*( normal(Ib1,Ib2,Ib3,0)*uex + normal(Ib1,Ib2,Ib3,1)*uey ) + a2*uexx + a3*ueyy;
//   }
//   else
//   {
//     // -------- CURVILINEAR GRID ---------

//     const Real dxs = mg.gridSpacing(axis);

//     const real a0 = Sl;
//     const real a1 = theta + Sl*dxs*theta;


//     const real a2 = dxs*( beta      ) + Sl*( .5*SQR(dxs)*beta      );
//     const real a3 = dxs*( (beta-1.) ) + Sl*( .5*SQR(dxs)*(beta-1.) );
//     const real a4 = a3; // for 3D 

//     printP("**** ADD TZ Correction to results from getData for CHAMP:  Sl=%g, theta=%g, beta=%g **************** FINISH ME *****************\n",Sl,theta,beta);  
//     int n=0;
//     assert( N.getLength()==1 );
//     interfaceData.u(Ib1,Ib2,Ib3,n) += a0*ue + a1*( normal(Ib1,Ib2,Ib3,0)*uex + normal(Ib1,Ib2,Ib3,1)*uey ) + a2*uexx + a3*ueyy;


//   }
// #endMacro


// =======================================================================================================
// / \brief  Fill in the RHS values for implicit time-stepping.
///
/// \param u (output) : the RHS
/// \param  (input) : for linearized equations
/// 
/// \note See ins/src/implicit.C for another example of this function is written.
// =======================================================================================================
int Cgad::
applyBoundaryConditionsForImplicitTimeStepping(realMappedGridFunction & u, 
                                                                                              realMappedGridFunction & uL,
                                                                                              realMappedGridFunction & uOld, // *wdh* Dec 20, 2017 -- added uOld
                                                                                              realMappedGridFunction & gridVelocity,
                                                                                              real t,
                                                                                              int scalarSystem,
                                                                                              int grid )

{
//  Overture::abort("Cgad::applyBoundaryConditionsForImplicitTimeStepping:ERROR: not implemented");
//  return 0;

    MappedGrid & mg = *u.getMappedGrid();
    const int & numberOfDimensions = mg.numberOfDimensions();

    const int orderOfAccuracy=parameters.dbase.get<int >("orderOfAccuracy");
    assert( orderOfAccuracy==2 || orderOfAccuracy==4 );  

    const bool twilightZoneFlow = parameters.dbase.get<bool >("twilightZoneFlow");

    Range N = parameters.dbase.get<Range >("Rt");   // time dependent variables

    const int & multiDomainProblem = parameters.dbase.get<int>("multiDomainProblem"); 
    const int applyChampInterfaceConditions = parameters.dbase.get<int>("applyChampInterfaceConditions");

    typedef int BoundaryCondition;
    
    const BoundaryCondition & dirichletBoundaryCondition = Parameters::dirichletBoundaryCondition;
    const BoundaryCondition & neumannBoundaryCondition = Parameters::neumannBoundaryCondition;
    const Parameters::BoundaryCondition & interfaceBoundaryCondition= Parameters::interfaceBoundaryCondition;
    
    BCTypes::BCNames dirichlet             = BCTypes::dirichlet,
                                      neumann               = BCTypes::neumann,
  //   extrapolate           = BCTypes::extrapolate,
  //   allBoundaries         = BCTypes::allBoundaries,
        normalComponent       = BCTypes::normalComponent;

    const RealArray & bcData = parameters.dbase.get<RealArray >("bcData");
    const int & numberOfComponents = parameters.dbase.get<int >("numberOfComponents");

#define mixedRHS(component,side,axis,grid)         bcData(component+numberOfComponents*(0),side,axis,grid)
#define mixedCoeff(component,side,axis,grid)       bcData(component+numberOfComponents*(1),side,axis,grid)
#define mixedNormalCoeff(component,side,axis,grid) bcData(component+numberOfComponents*(2),side,axis,grid)

  // *** assign boundary conditions for the implicit method ***** 

    if( parameters.getGridIsImplicit(grid) )
    {
    // ** Note that we are assigning the RHS for the implicit solve ***
        if( debug() & 2 )
            printP("AD: applyBCForImplicitTimeStepping: fill RHS for implicit solve, orderOfAccuracy=%d, t=%9.3e\n",orderOfAccuracy,t);

    // The RHS for the interface equations is placed in the BoundaryDataArray
        BoundaryData::BoundaryDataArray & pBoundaryData = parameters.getBoundaryData(grid);


        const Parameters::KnownSolutionsEnum & knownSolution = parameters.dbase.get<Parameters::KnownSolutionsEnum >("knownSolution");
        const bool & assignKnownSolutionAtBoundaries = parameters.dbase.get<bool>("assignKnownSolutionAtBoundaries");
        if( knownSolution == Parameters::userDefinedKnownSolution && assignKnownSolutionAtBoundaries  )
        {
      // ----- Fill in Dirichlet BCs with the known solution ------
            OV_GET_SERIAL_ARRAY(real,u,uLocal);
            RealArray ua; 

            Index Ib1,Ib2,Ib3;
            ForBoundary(side,axis)
            {
                if( mg.boundaryCondition(side,axis) == AdParameters::dirichletBoundaryCondition )
                {
                    getBoundaryIndex(mg.gridIndexRange(),side,axis,Ib1,Ib2,Ib3);
                    ua = parameters.getKnownSolution(t, grid, Ib1,Ib2,Ib3 );  
                    uLocal(Ib1,Ib2,Ib3,N) = ua(Ib1,Ib2,Ib3,N); 


          //  getUserDefinedKnownSolution(real t, CompositeGrid & cg, int grid, RealArray & ua, 
          // const Index & I1, const Index &I2, const Index &I3, int numberOfTimeDerivatives = 0 );        

                }

            }

        }
        else
        {
            u.applyBoundaryCondition(N,dirichlet,dirichletBoundaryCondition,bcData,pBoundaryData,t,
                                                              Overture::defaultBoundaryConditionParameters(),grid);
        }

    // ** Note: neumann BC is applied below since we fill in the RHS 
    // u.applyBoundaryCondition(N,neumann,neumannBoundaryCondition,bcData,t,
    //                          Overture::defaultBoundaryConditionParameters(),grid);
    // u.applyBoundaryCondition(N,dirichlet,interfaceBoundaryCondition,bcData,t,
    //                       Overture::defaultBoundaryConditionParameters(),grid);
        

        const bool isRectangular=mg.isRectangular();
        if( !isRectangular || twilightZoneFlow ) // we need this now
            mg.update(MappedGrid::THEvertexBoundaryNormal); 

        OV_GET_SERIAL_ARRAY(real,u,uLocal);
    // #ifdef USE_PPP
    //   const realSerialArray & uLocal = u.getLocalArray();
    // #else
    //   realSerialArray & uLocal = u; 
    // #endif


        real dx[3]={1.,1.,1.};
        if( isRectangular )
            mg.getDeltaX(dx);

        const bool rectangular= mg.isRectangular() && !twilightZoneFlow;
            
        realArray & x= mg.center();
        #ifdef USE_PPP
            realSerialArray xLocal; 
            if( !rectangular || twilightZoneFlow ) 
                getLocalArrayWithGhostBoundaries(x,xLocal);
        #else
            const realSerialArray & xLocal = x;
        #endif
        
        const IntegerArray & interfaceType = parameters.dbase.get<IntegerArray >("interfaceType");

        const bool interfaceBoundaryConditionsAreSpecified=parameters.dbase.has_key("interfaceCondition");
        IntegerArray & interfaceCondition = (interfaceBoundaryConditionsAreSpecified ? 
                                                                                  parameters.dbase.get<IntegerArray>("interfaceCondition") :
                                                                                  Overture::nullIntArray() );

        int side,axis;
        Index Ibv[3], &Ib1=Ibv[0], &Ib2=Ibv[1], &Ib3=Ibv[2];
        Index Igv[3], &Ig1=Igv[0], &Ig2=Igv[1], &Ig3=Igv[2];

        Index Jbv[3], &Jb1=Jbv[0], &Jb2=Jbv[1], &Jb3=Jbv[2];
        Index Jgv[3], &Jg1=Jgv[0], &Jg2=Jgv[1], &Jg3=Jgv[2];

        Index Ip1,Ip2,Ip3;
        for( axis=0; axis<mg.numberOfDimensions(); axis++ ) 
        {
            for( side=0; side<=1; side++ )
            {

                int bc = mg.boundaryCondition(side,axis);
        // if( bc == Parameters::interfaceBoundaryCondition && interfaceBoundaryConditionsAreSpecified )
                if( interfaceType(side,axis,grid) != Parameters::noInterface )
                {
          // this face is on a domain interface
                    bc = interfaceCondition(side,axis,grid);
                    assert( bc==Parameters::dirichletInterface || bc==Parameters::neumannInterface );
                }


        // **NEW** April 8, 2021 -- ADD OPTION FOR REQUEST DATA ON DEMAND --   ... FINISH ME 
                const Parameters::InterfaceCommunicationModeEnum & interfaceCommunicationMode = 
                                                    parameters.dbase.get<Parameters::InterfaceCommunicationModeEnum>("interfaceCommunicationMode");
                if( interfaceCommunicationMode==Parameters::requestInterfaceDataWhenNeeded && 
                        interfaceType(side,axis,grid) != Parameters::noInterface )
                {

                    if( debug() & 4 && interfaceCommunicationMode==Parameters::requestInterfaceDataWhenNeeded  )
                    {
                        printP("applyBCIMP:INFO: interfaceCommunicationMode==requestInterfaceDataWhenNeeded\n");
                    }

                    bool testNew=true;
                    if( true && (testNew || interfaceCommunicationMode==Parameters::requestInterfaceDataWhenNeeded) )
                    {
            // *** TEST NEW WAY **** April 8, 2021 
                        if( debug() & 4 )
                        {
                            printP("applyBCImp:REQUEST DATA for (grid,side,axis)=(%d,%d,%d) interfaceType(side,axis,grid)=%d\n",grid,side,axis,interfaceType(side,axis,grid));
                            if( bc==Parameters::dirichletInterface )
                                printP("applyBCImp: dirichletInterface: TEST:request heatFluxInterface interface data, t=%9.3e\n",t);
                            else
                                printP("applyBCImp: neumannInterface: TEST:request heatFluxInterface data, t=%9.3e\n",t);
                        }

                        InterfaceData interfaceData;
            // Range Rx=numberOfDimensions;
                        getBoundaryIndex(mg.gridIndexRange(),side,axis,Ib1,Ib2,Ib3);
                        getGhostIndex(   mg.gridIndexRange(),side,axis,Ig1,Ig2,Ig3);
            // *wdh* April 16, 2022
            //getBoundaryIndex(mg.indexRange(),side,axis,Ib1,Ib2,Ib3);
            //getGhostIndex(   mg.indexRange(),side,axis,Ig1,Ig2,Ig3);

                        interfaceData.u.redim(Ib1,Ib2,Ib3,N); // heat flux is returned here 
                        interfaceData.t=t;
                        interfaceData.u=0;
            // ::display(interfaceData.u(Ib1,Ib2,Ib3,N),"BEFORE REQUEST DATA FOR INTERFACE","%7.4f ");

            // Boundary data should go here
                        RealArray & bd = parameters.getBoundaryData(side,axis,grid,mg);
            // bd=0;

                        if( debug() & 4 )
                        {
                            printF("\n");
                            printP(" >>>> REQUEST INTERFACE DATA for (grid,side,axis)=(%d,%d,%d) t=%9.3e\n",grid,side,axis,t);
                        }
                        
            // We pass the coefficients in the dirichlet, neumann or mixed BC

            // From ins/src/addedMassImplicitBoundaryConditions.bC: 
            //Parameters & bulkSolidParams = getInterfaceParameters( grid,side,axis,parameters );

                        GridFaceDescriptor & myGFD = getInterfaceGridFaceDescriptor( grid, side, axis, parameters );
                        real *sourceMixedBC = myGFD.dbase.get<real[2]>("sourceMixedBC");



                        GridFaceDescriptor targetInfo(-1,grid,side,axis);
                        const int n=0; // component -- FIX ME
                        const real a0=mixedCoeff(n,side,axis,grid), a1=mixedNormalCoeff(n,side,axis,grid);
            // targetInfo.a[0]=a0; targetInfo.a[1]=-a1; targetInfo.a[2]=0.;  // **NOTE*** flip sign of a1 for normal for other side, *check me*

                        targetInfo.a[0]=sourceMixedBC[0]; targetInfo.a[1]=sourceMixedBC[1]; targetInfo.a[2]=0.; 
            // const real a0 = targetInfo.a[0], a1= targetInfo.a[1]; 

                        if( debug() & 2 )
                            printF("--- applyBC: After getInterfaceGridFaceDescriptor: [a0,a1]=[%g,%g] myGFD.a=[%g,%g], sourceMixedBC=[%g,%g]\n",
                                              a0,a1,myGFD.a[0],myGFD.a[1],sourceMixedBC[0],sourceMixedBC[1]);


                        if( applyChampInterfaceConditions )
                        { 
               // CHAMP CONDITIONS
              // Evaluate   S*uR(dx) + n.grad( uR(dx) ) 
              // for the solution on the other side ("right") of the interface, at ONE grid line inside the domain
                            const RealArray & champParameters = parameters.dbase.get<RealArray>("champParameters");
                            const Real Sl    = champParameters(4,side,axis,grid);    // pl/dxs; 
              // const real Sl    = champParameters(0,side,axis,grid);
              // const real theta = champParameters(2,side,axis,grid);
              // const real beta  = champParameters(3,side,axis,grid);  

                            targetInfo.a[0] =  Sl;  // optimized Schwartz parameter for this domain ("left")
                            targetInfo.a[1] = -1.;  // flip sign of normal
              // Do this for now: 
                            targetInfo.dbase.put<int>("getChampData")=1;
                        }
                        else
                        {
              // What should we do here ?
                        }

                        int interfaceDataOptions = Parameters::heatFluxInterfaceData;
                        bool saveTimeHistory=true;

                        if( debug() & 8 )
                        {
                            printF("\n");
                            printP(" >>>> REQUEST INTERFACE DATA for (grid,side,axis)=(%d,%d,%d) t=%9.3e. source: (a0,a1)=(%g,%g),  target: (a0,a1)=(%g,%g) (champ=%d)\n",
                                          grid,side,axis,t,a0,a1,targetInfo.a[0],targetInfo.a[1],(int)applyChampInterfaceConditions);
                        }

                        getInterfaceData( t, grid, side, axis, 
                                                            interfaceDataOptions,
                                                            interfaceData.u,
                                                            parameters,saveTimeHistory,
                                                            &targetInfo );

                        if( twilightZoneFlow  )
                        {
                            if( applyChampInterfaceConditions )
                            {
                                    Index Ib1,Ib2,Ib3; 
                                    getBoundaryIndex(mg.indexRange(),side,axis,Ib1,Ib2,Ib3);
                                    OGFunction & e = *(parameters.dbase.get<OGFunction* >("exactSolution"));
                                    const RealArray & champParameters = parameters.dbase.get<RealArray>("champParameters");
                  // const real Sl    = champParameters(0,side,axis,grid);
                  // const real theta = champParameters(2,side,axis,grid);
                  // const real beta  = champParameters(3,side,axis,grid);  
                                    const Real pl    = champParameters(0,side,axis,grid);    // optimized Scwartz Parameter for side 1
                                    const Real pr    = champParameters(1,side,axis,grid);    // optimized Scwartz Parameter for side 2
                                    const Real theta = champParameters(2,side,axis,grid);    // K1/K2
                                    const Real beta  = champParameters(3,side,axis,grid);    // D1/D2    
                                    const Real Sl    = champParameters(4,side,axis,grid);    // pl/dxs; 
                                    const Real KLR = theta;
                                    const Real DLR = beta;
                  // printF("**** ADD TZ Correction to results from getData for CHAMP:  Sl=%g, theta=%g, beta=%g *********************************\n",Sl,theta,beta);                  
                                    OV_GET_VERTEX_BOUNDARY_NORMAL(mg,side,axis,normal);
                                    RealArray ue(Ib1,Ib2,Ib3,N), uex(Ib1,Ib2,Ib3,N), uey(Ib1,Ib2,Ib3,N), uexx(Ib1,Ib2,Ib3,N), ueyy(Ib1,Ib2,Ib3,N);
                                    int rectangular=0;
                                    e.gd( ue  ,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,0,Ib1,Ib2,Ib3,N,t);
                                    e.gd( uex ,xLocal,mg.numberOfDimensions(),rectangular,0,1,0,0,Ib1,Ib2,Ib3,N,t);
                                    e.gd( uey ,xLocal,mg.numberOfDimensions(),rectangular,0,0,1,0,Ib1,Ib2,Ib3,N,t);
                                    e.gd( uexx,xLocal,mg.numberOfDimensions(),rectangular,0,2,0,0,Ib1,Ib2,Ib3,N,t);
                                    e.gd( ueyy,xLocal,mg.numberOfDimensions(),rectangular,0,0,2,0,Ib1,Ib2,Ib3,N,t);
                                    if( isRectangular )
                                    {
                    // const Real dxs = dx[axis];   // *** FIX ME need dx[opposite-side]
                    // We use dx fro the opposite side
                                        const Real dxs = champParameters(5,side,axis,grid); // value saved in champBoundaryConditions.bC 
                                        if( orderOfAccuracy==2 )
                                        {
                                            const real a0 = Sl;
                                            const real a1 = theta + Sl*dxs*theta;
                      // const real a2 = dxs*( beta      ) + Sl*( .5*SQR(dxs)*beta      );
                      // const real a3 = dxs*( (beta-1.) ) + Sl*( .5*SQR(dxs)*(beta-1.) );
                      // const real a4 = a3; // for 3D 
                                            real axx, ayy, azz;
                                            if( axis==0 )
                                            {
                                                axx = (beta   )*( dxs + Sl*( .5*SQR(dxs) ) );
                                                ayy = (beta-1.)*( dxs + Sl*( .5*SQR(dxs) ) );
                                                azz = ayy; // for 3D 
                                            }
                                            else if( axis==1 )
                                            {
                                                ayy = (beta   )*( dxs + Sl*( .5*SQR(dxs) ) );
                                                axx = (beta-1.)*( dxs + Sl*( .5*SQR(dxs) ) );
                                                azz = axx; // for 3D       
                                            }
                                            assert( numberOfDimensions==2 );
                                            int n=0;
                                            assert( N.getLength()==1 );
                                            interfaceData.u(Ib1,Ib2,Ib3,n) += a0*ue + a1*( normal(Ib1,Ib2,Ib3,0)*uex + normal(Ib1,Ib2,Ib3,1)*uey ) + axx*uexx + ayy*ueyy;
                                        }
                                        else if( orderOfAccuracy==4 )
                                        {
                                            if( true )
                                                printF("Add TZ correction for CHAMP4 at t=%9.3e, KLR=%g\n",t,KLR);
                                            int n=0;
                                            assert( N.getLength()==1 );
                                            if( axis==0 )
                                            {
                                                RealArray uexxx(Ib1,Ib2,Ib3,N), uexyy(Ib1,Ib2,Ib3,N), uexxxx(Ib1,Ib2,Ib3,N), uexxyy(Ib1,Ib2,Ib3,N), ueyyyy(Ib1,Ib2,Ib3,N);
                                                e.gd( uexxx ,xLocal,mg.numberOfDimensions(),rectangular,0,3,0,0,Ib1,Ib2,Ib3,N,t);
                                                e.gd( uexyy ,xLocal,mg.numberOfDimensions(),rectangular,0,1,2,0,Ib1,Ib2,Ib3,N,t);
                                                e.gd( uexxxx,xLocal,mg.numberOfDimensions(),rectangular,0,4,0,0,Ib1,Ib2,Ib3,N,t);
                                                e.gd( uexxyy,xLocal,mg.numberOfDimensions(),rectangular,0,2,2,0,Ib1,Ib2,Ib3,N,t);
                                                e.gd( ueyyyy,xLocal,mg.numberOfDimensions(),rectangular,0,0,4,0,Ib1,Ib2,Ib3,N,t);
                                                RealArray L1(Ib1,Ib2,Ib3,N),L2(Ib1,Ib2,Ib3,N),L3(Ib1,Ib2,Ib3,N),L4(Ib1,Ib2,Ib3,N);
                                                const Real nSign = 2*side-1;
                                                L1 = (KLR*nSign)*uex;
                                                L2 = DLR*uexx + (DLR-1.)*ueyy;
                                                L3 = (KLR*DLR*nSign)*( uexxx ) + ((KLR*DLR-KLR)*nSign)*( uexyy );
                                                L4 = (DLR*DLR)*( uexxxx + 2.*uexxyy + ueyyyy ) - ueyyyy - (2.*DLR)*( uexxyy + ueyyyy) + 2.*ueyyyy;
                                                interfaceData.u(Ib1,Ib2,Ib3,n) += Sl*( ue + dxs*L1 + (dxs*dxs*.5)*L2 + (dxs*dxs*dxs/6.)*L3 + (dxs*dxs*dxs*dxs/24.)*L4 )
                                                                                                                          + L1 + dxs*L2 + (dxs*dxs*.5)*L3 + (dxs*dxs*dxs/6.)*L4;
                                            }
                                            else if( axis==1 )
                                            {
                                                OV_ABORT("finish me");
                                            }
                                        }
                                        else
                                        {
                                            OV_ABORT("error -- orderOfAccuracy");
                                        }
                                    }
                                    else
                                    {
                    // -------- CURVILINEAR GRID ---------
                                        const Real dxs = mg.gridSpacing(axis);
                                        const real a0 = Sl;
                                        const real a1 = theta + Sl*dxs*theta;
                                        const real a2 = dxs*( beta      ) + Sl*( .5*SQR(dxs)*beta      );
                                        const real a3 = dxs*( (beta-1.) ) + Sl*( .5*SQR(dxs)*(beta-1.) );
                                        const real a4 = a3; // for 3D 
                    //  printF("**** ADD TZ Correction to results from getData for CHAMP:  Sl=%g, theta=%g, beta=%g **************** FINISH ME *****************\n",Sl,theta,beta); 
                    // NCoeff(m) = theta*( an1L*xCoeff(m,i1,i2,i3) + an2L*yCoeff(m,i1,i2,i3) ) - b2R*r2Coeff(m,j1,j2,j3); 
                    // // NOTE: Here we assume the tangential r derivatives are the same on both sides! (see champ4/notes) 
                    // LCoeff(m) = (b1R*b1R/c11R)*( beta*lapCoeff(m,i1,i2,i3) 
                    //                             - (c12R+c21R)*( theta*( (b1L/b1R)*r1r2Coeff(m,i1,i2,i3) + (b2L/b1R)*r2r2Coeff(m,j1,j2,j3) + b2Lb1Rr2*r2Coeff(m,j1,j2,j3) )
                    //                                             + b2Rb1Rr2*r2Coeff(m,j1,j2,j3) + (b2R/b1R)*r2r2Coeff(m,j1,j2,j3) )
                    //                             - c22R*r2r2Coeff(m,j1,j2,j3) 
                    //                             + (c1R/b1R)*NCoeff(m) 
                    //                             + c2R*r2Coeff(m,j1,j2,j3) );
                                        const Real dr1 = mg.gridSpacing(axis); 
                    // -- compute tangential derivatives w.r.t r2  ---
                                        const int axisp1 = (axis+1) % mg.numberOfDimensions();
                                        const Real dr2 = mg.gridSpacing(axisp1); 
                                        int jsv[3]={0,0,0}, &js1=jsv[0], &js2=jsv[1], &js3=jsv[2];
                                        jsv[axisp1]=1; 
                                        int n=0;
                                        assert( N.getLength()==1 );
                    // interfaceData.u(Ib1,Ib2,Ib3,n) += a0*ue + a1*( normal(Ib1,Ib2,Ib3,0)*uex + normal(Ib1,Ib2,Ib3,1)*uey ) + a2*uexx + a3*ueyy;
                    // Some coefficients in the CHAMP conditions were saved when the matrix was created (champBoundaryConditions.bC)
                                        RealArray & cc = myGFD.dbase.get<RealArray>("ccChamp");
                                        if( orderOfAccuracy==2 )
                                        {
                                            const RealArray & dra = mg.gridSpacing(); 
                                            real dr[3]={dra(0),dra(1),dra(2)};
                      // evaluate the exact solution at points on and near the boundary
                                            Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
                                            I1=Ib1; I2=Ib2; I3=Ib3;
                                            for( int dir=0; dir<numberOfDimensions; dir++ )
                                            {
                                                Iv[dir] = Range( Iv[dir].getBase()-2, Iv[dir].getBound()+2 );  // we have a 5-pt stencil
                                            }
                                            RealArray ue(I1,I2,I3,N);
                                            e.gd( ue,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,0,I1,I2,I3,N,t);
                                            RealArray    uer(Ib1,Ib2,Ib3),    ues(Ib1,Ib2,Ib3);
                                            RealArray   uerr(Ib1,Ib2,Ib3),   uers(Ib1,Ib2,Ib3),   uess(Ib1,Ib2,Ib3);
                                            uer  = ((ue(Ib1+1,Ib2,Ib3)-ue(Ib1-1,Ib2,Ib3))/(2.*dr[0]));
                                            ues  = ((ue(Ib1,Ib2+1,Ib3)-ue(Ib1,Ib2-1,Ib3))/(2.*dr[1]));
                                            uerr = ((ue(Ib1+1,Ib2,Ib3)-2.*ue(Ib1,Ib2,Ib3)+ue(Ib1-1,Ib2,Ib3))/(SQR(dr[0])));
                                            uers = ((((ue(Ib1+1,Ib2+1,Ib3)-ue(Ib1+1,Ib2-1,Ib3))/(2.*dr[1]))-((ue(Ib1-1,Ib2+1,Ib3)-ue(Ib1-1,Ib2-1,Ib3))/(2.*dr[1])))/(2.*dr[0]));
                                            uess = ((ue(Ib1,Ib2+1,Ib3)-2.*ue(Ib1,Ib2,Ib3)+ue(Ib1,Ib2-1,Ib3))/(SQR(dr[1])));
                      // interfaceData.u(Ib1,Ib2,Ib3,n) += Sl*ue(Ib1,Ib2,Ib3); 
                      // interfaceData.u(Ib1,Ib2,Ib3,n) +=  cc(Ib1,Ib2,Ib3, 0)*ue(Ib1,Ib2,Ib3);
                                            interfaceData.u(Ib1,Ib2,Ib3,n) +=  cc(Ib1,Ib2,Ib3, 0)*ue(Ib1,Ib2,Ib3)
                                                                                                                +cc(Ib1,Ib2,Ib3, 1)*uer
                                                                                                                +cc(Ib1,Ib2,Ib3, 2)*ues
                                                                                                                +cc(Ib1,Ib2,Ib3, 3)*uerr
                                                                                                                +cc(Ib1,Ib2,Ib3, 4)*uers
                                                                                                                +cc(Ib1,Ib2,Ib3, 5)*uess;
                      // RealArray uep(Ib1+js1,Ib2+js2,Ib3+js3,N), uem(Ib1-js1,Ib2-js2,Ib3-js3,N);
                      // e.gd( uep,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,0,Ib1+js1,Ib2+js2,Ib3+js3,N,t);
                      // e.gd( uem,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,0,Ib1-js1,Ib2-js2,Ib3-js3,N,t);
                      // RealArray uer2(Ib1,Ib2,Ib3,N), uer2r2(Ib1,Ib2,Ib3,N);
                      // uer2   = ( uep - uem )*(1./(2.*dr2));
                      // uer2r2 = ( uep -2.*ue + uem )*(1./(dr2*dr2));
                      // // -- compute the mixed derivative w.r.t r1 and r2 ----
                      // RealArray uepp(Ib1+1,Ib2+1,Ib3,N), uemm(Ib1-1,Ib2-1,Ib3,N), uepm(Ib1+1,Ib2-1,Ib3,N), uemp(Ib1-1,Ib2+1,Ib3,N);
                      // e.gd( uemm,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,0,Ib1-1,Ib2-1,Ib3,N,t);
                      // e.gd( uepm,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,0,Ib1+1,Ib2-1,Ib3,N,t);
                      // e.gd( uemp,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,0,Ib1-1,Ib2+1,Ib3,N,t);
                      // e.gd( uepp,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,0,Ib1+1,Ib2+1,Ib3,N,t);
                      // RealArray uer1r2(Ib1,Ib2,Ib3,N);
                      // uer1r2 = ( uepp - uemp -  uepm + uemm )*(1./(4.*dr1*dr2));   // D0r1 * D0r2       
                      // RealArray Nc(Ib1,Ib2,Ib3,N);
                      // Nc = theta*( normal(Ib1,Ib2,Ib3,0)*uex + normal(Ib1,Ib2,Ib3,1)*uey ) + cc(Ib1,Ib2,Ib3,0)*uer2; 
                      // interfaceData.u(Ib1,Ib2,Ib3,n) += Sl*ue   
                      //                                   + cc(Ib1,Ib2,Ib3,1)*( uexx + ueyy ) 
                      //                                   + cc(Ib1,Ib2,Ib3,2)*uer1r2 
                      //                                   + cc(Ib1,Ib2,Ib3,3)*uer2r2
                      //                                   + cc(Ib1,Ib2,Ib3,4)*uer2
                      //                                   + cc(Ib1,Ib2,Ib3,5)*Nc;
                                        }
                                        else
                                        {
                                            const RealArray & dra = mg.gridSpacing(); 
                                            real dr[3]={dra(0),dra(1),dra(2)};
                      // evaluate the exact solution at points on and near the boundary
                                            Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
                                            I1=Ib1; I2=Ib2; I3=Ib3;
                                            for( int dir=0; dir<numberOfDimensions; dir++ )
                                            {
                                                Iv[dir] = Range( Iv[dir].getBase()-2, Iv[dir].getBound()+2 );  // we have a 5-pt stencil
                                            }
                                            RealArray ue(I1,I2,I3,N);
                                            e.gd( ue,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,0,I1,I2,I3,N,t);
                                            RealArray    uer(Ib1,Ib2,Ib3),    ues(Ib1,Ib2,Ib3);
                                            RealArray   uerr(Ib1,Ib2,Ib3),   uers(Ib1,Ib2,Ib3),   uess(Ib1,Ib2,Ib3);
                                            RealArray  uerrr(Ib1,Ib2,Ib3),  uerrs(Ib1,Ib2,Ib3),  uerss(Ib1,Ib2,Ib3),  uesss(Ib1,Ib2,Ib3);
                                            RealArray uerrrr(Ib1,Ib2,Ib3), uerrrs(Ib1,Ib2,Ib3), uerrss(Ib1,Ib2,Ib3), uersss(Ib1,Ib2,Ib3), uessss(Ib1,Ib2,Ib3);
                                            uer  = (((1./12.)*ue(Ib1-2,Ib2,Ib3)-(2./3.)*ue(Ib1-1,Ib2,Ib3)+(2./3.)*ue(Ib1+1,Ib2,Ib3)-(1./12.)*ue(Ib1+2,Ib2,Ib3))/(dr[0]));
                                            ues  = (((1./12.)*ue(Ib1,Ib2-2,Ib3)-(2./3.)*ue(Ib1,Ib2-1,Ib3)+(2./3.)*ue(Ib1,Ib2+1,Ib3)-(1./12.)*ue(Ib1,Ib2+2,Ib3))/(dr[1]));
                                            uerr = ((-1.*ue(Ib1-2,Ib2,Ib3)+16.*ue(Ib1-1,Ib2,Ib3)-30.*ue(Ib1,Ib2,Ib3)+16.*ue(Ib1+1,Ib2,Ib3)-1.*ue(Ib1+2,Ib2,Ib3))/(12.*SQR(dr[0])));
                                            uers = (((1./12.)*(((1./12.)*ue(Ib1-2,Ib2-2,Ib3)-(2./3.)*ue(Ib1-2,Ib2-1,Ib3)+(2./3.)*ue(Ib1-2,Ib2+1,Ib3)-(1./12.)*ue(Ib1-2,Ib2+2,Ib3))/(dr[1]))-(2./3.)*(((1./12.)*ue(Ib1-1,Ib2-2,Ib3)-(2./3.)*ue(Ib1-1,Ib2-1,Ib3)+(2./3.)*ue(Ib1-1,Ib2+1,Ib3)-(1./12.)*ue(Ib1-1,Ib2+2,Ib3))/(dr[1]))+(2./3.)*(((1./12.)*ue(Ib1+1,Ib2-2,Ib3)-(2./3.)*ue(Ib1+1,Ib2-1,Ib3)+(2./3.)*ue(Ib1+1,Ib2+1,Ib3)-(1./12.)*ue(Ib1+1,Ib2+2,Ib3))/(dr[1]))-(1./12.)*(((1./12.)*ue(Ib1+2,Ib2-2,Ib3)-(2./3.)*ue(Ib1+2,Ib2-1,Ib3)+(2./3.)*ue(Ib1+2,Ib2+1,Ib3)-(1./12.)*ue(Ib1+2,Ib2+2,Ib3))/(dr[1])))/(dr[0]));
                                            uess = ((-1.*ue(Ib1,Ib2-2,Ib3)+16.*ue(Ib1,Ib2-1,Ib3)-30.*ue(Ib1,Ib2,Ib3)+16.*ue(Ib1,Ib2+1,Ib3)-1.*ue(Ib1,Ib2+2,Ib3))/(12.*SQR(dr[1])));
                                            uerrr = (((-1./2.)*ue(Ib1-2,Ib2,Ib3)+1.*ue(Ib1-1,Ib2,Ib3)-1.*ue(Ib1+1,Ib2,Ib3)+(1./2.)*ue(Ib1+2,Ib2,Ib3))/(pow(dr[0],3)));
                                            uerrs = ((-1.*(((1./12.)*ue(Ib1-2,Ib2-2,Ib3)-(2./3.)*ue(Ib1-2,Ib2-1,Ib3)+(2./3.)*ue(Ib1-2,Ib2+1,Ib3)-(1./12.)*ue(Ib1-2,Ib2+2,Ib3))/(dr[1]))+16.*(((1./12.)*ue(Ib1-1,Ib2-2,Ib3)-(2./3.)*ue(Ib1-1,Ib2-1,Ib3)+(2./3.)*ue(Ib1-1,Ib2+1,Ib3)-(1./12.)*ue(Ib1-1,Ib2+2,Ib3))/(dr[1]))-30.*(((1./12.)*ue(Ib1,Ib2-2,Ib3)-(2./3.)*ue(Ib1,Ib2-1,Ib3)+(2./3.)*ue(Ib1,Ib2+1,Ib3)-(1./12.)*ue(Ib1,Ib2+2,Ib3))/(dr[1]))+16.*(((1./12.)*ue(Ib1+1,Ib2-2,Ib3)-(2./3.)*ue(Ib1+1,Ib2-1,Ib3)+(2./3.)*ue(Ib1+1,Ib2+1,Ib3)-(1./12.)*ue(Ib1+1,Ib2+2,Ib3))/(dr[1]))-1.*(((1./12.)*ue(Ib1+2,Ib2-2,Ib3)-(2./3.)*ue(Ib1+2,Ib2-1,Ib3)+(2./3.)*ue(Ib1+2,Ib2+1,Ib3)-(1./12.)*ue(Ib1+2,Ib2+2,Ib3))/(dr[1])))/(12.*SQR(dr[0])));
                                            uerss = (((1./12.)*((-1.*ue(Ib1-2,Ib2-2,Ib3)+16.*ue(Ib1-2,Ib2-1,Ib3)-30.*ue(Ib1-2,Ib2,Ib3)+16.*ue(Ib1-2,Ib2+1,Ib3)-1.*ue(Ib1-2,Ib2+2,Ib3))/(12.*SQR(dr[1])))-(2./3.)*((-1.*ue(Ib1-1,Ib2-2,Ib3)+16.*ue(Ib1-1,Ib2-1,Ib3)-30.*ue(Ib1-1,Ib2,Ib3)+16.*ue(Ib1-1,Ib2+1,Ib3)-1.*ue(Ib1-1,Ib2+2,Ib3))/(12.*SQR(dr[1])))+(2./3.)*((-1.*ue(Ib1+1,Ib2-2,Ib3)+16.*ue(Ib1+1,Ib2-1,Ib3)-30.*ue(Ib1+1,Ib2,Ib3)+16.*ue(Ib1+1,Ib2+1,Ib3)-1.*ue(Ib1+1,Ib2+2,Ib3))/(12.*SQR(dr[1])))-(1./12.)*((-1.*ue(Ib1+2,Ib2-2,Ib3)+16.*ue(Ib1+2,Ib2-1,Ib3)-30.*ue(Ib1+2,Ib2,Ib3)+16.*ue(Ib1+2,Ib2+1,Ib3)-1.*ue(Ib1+2,Ib2+2,Ib3))/(12.*SQR(dr[1]))))/(dr[0]));
                                            uesss = (((-1./2.)*ue(Ib1,Ib2-2,Ib3)+1.*ue(Ib1,Ib2-1,Ib3)-1.*ue(Ib1,Ib2+1,Ib3)+(1./2.)*ue(Ib1,Ib2+2,Ib3))/(pow(dr[1],3)));
                                            uerrrr = ((ue(Ib1-2,Ib2,Ib3)-4.*ue(Ib1-1,Ib2,Ib3)+6.*ue(Ib1,Ib2,Ib3)-4.*ue(Ib1+1,Ib2,Ib3)+ue(Ib1+2,Ib2,Ib3))/(pow(dr[0],4)));
                                            uerrrs = (((-1./2.)*(((1./12.)*ue(Ib1-2,Ib2-2,Ib3)-(2./3.)*ue(Ib1-2,Ib2-1,Ib3)+(2./3.)*ue(Ib1-2,Ib2+1,Ib3)-(1./12.)*ue(Ib1-2,Ib2+2,Ib3))/(dr[1]))+1.*(((1./12.)*ue(Ib1-1,Ib2-2,Ib3)-(2./3.)*ue(Ib1-1,Ib2-1,Ib3)+(2./3.)*ue(Ib1-1,Ib2+1,Ib3)-(1./12.)*ue(Ib1-1,Ib2+2,Ib3))/(dr[1]))-1.*(((1./12.)*ue(Ib1+1,Ib2-2,Ib3)-(2./3.)*ue(Ib1+1,Ib2-1,Ib3)+(2./3.)*ue(Ib1+1,Ib2+1,Ib3)-(1./12.)*ue(Ib1+1,Ib2+2,Ib3))/(dr[1]))+(1./2.)*(((1./12.)*ue(Ib1+2,Ib2-2,Ib3)-(2./3.)*ue(Ib1+2,Ib2-1,Ib3)+(2./3.)*ue(Ib1+2,Ib2+1,Ib3)-(1./12.)*ue(Ib1+2,Ib2+2,Ib3))/(dr[1]))/(pow(dr[0],3))));
                                            uerrss = ((-1.*((-1.*ue(Ib1-2,Ib2-2,Ib3)+16.*ue(Ib1-2,Ib2-1,Ib3)-30.*ue(Ib1-2,Ib2,Ib3)+16.*ue(Ib1-2,Ib2+1,Ib3)-1.*ue(Ib1-2,Ib2+2,Ib3))/(12.*SQR(dr[1])))+16.*((-1.*ue(Ib1-1,Ib2-2,Ib3)+16.*ue(Ib1-1,Ib2-1,Ib3)-30.*ue(Ib1-1,Ib2,Ib3)+16.*ue(Ib1-1,Ib2+1,Ib3)-1.*ue(Ib1-1,Ib2+2,Ib3))/(12.*SQR(dr[1])))-30.*((-1.*ue(Ib1,Ib2-2,Ib3)+16.*ue(Ib1,Ib2-1,Ib3)-30.*ue(Ib1,Ib2,Ib3)+16.*ue(Ib1,Ib2+1,Ib3)-1.*ue(Ib1,Ib2+2,Ib3))/(12.*SQR(dr[1])))+16.*((-1.*ue(Ib1+1,Ib2-2,Ib3)+16.*ue(Ib1+1,Ib2-1,Ib3)-30.*ue(Ib1+1,Ib2,Ib3)+16.*ue(Ib1+1,Ib2+1,Ib3)-1.*ue(Ib1+1,Ib2+2,Ib3))/(12.*SQR(dr[1])))-1.*((-1.*ue(Ib1+2,Ib2-2,Ib3)+16.*ue(Ib1+2,Ib2-1,Ib3)-30.*ue(Ib1+2,Ib2,Ib3)+16.*ue(Ib1+2,Ib2+1,Ib3)-1.*ue(Ib1+2,Ib2+2,Ib3))/(12.*SQR(dr[1]))))/(12.*SQR(dr[0])));
                                            uersss = (((1./12.)*(((-1./2.)*ue(Ib1-2,Ib2-2,Ib3)+1.*ue(Ib1-2,Ib2-1,Ib3)-1.*ue(Ib1-2,Ib2+1,Ib3)+(1./2.)*ue(Ib1-2,Ib2+2,Ib3))/(pow(dr[1],3)))-(2./3.)*(((-1./2.)*ue(Ib1-1,Ib2-2,Ib3)+1.*ue(Ib1-1,Ib2-1,Ib3)-1.*ue(Ib1-1,Ib2+1,Ib3)+(1./2.)*ue(Ib1-1,Ib2+2,Ib3))/(pow(dr[1],3)))+(2./3.)*(((-1./2.)*ue(Ib1+1,Ib2-2,Ib3)+1.*ue(Ib1+1,Ib2-1,Ib3)-1.*ue(Ib1+1,Ib2+1,Ib3)+(1./2.)*ue(Ib1+1,Ib2+2,Ib3))/(pow(dr[1],3)))-(1./12.)*(((-1./2.)*ue(Ib1+2,Ib2-2,Ib3)+1.*ue(Ib1+2,Ib2-1,Ib3)-1.*ue(Ib1+2,Ib2+1,Ib3)+(1./2.)*ue(Ib1+2,Ib2+2,Ib3))/(pow(dr[1],3))))/(dr[0]));
                                            uessss = ((ue(Ib1,Ib2-2,Ib3)-4.*ue(Ib1,Ib2-1,Ib3)+6.*ue(Ib1,Ib2,Ib3)-4.*ue(Ib1,Ib2+1,Ib3)+ue(Ib1,Ib2+2,Ib3))/(pow(dr[1],4)));
                      // printF("TZ correction: champ4: Sl=%9.3e \n",Sl);
                      // ::display(cc(Ib1,Ib2,Ib3, 0),"cc(Ib1,Ib2,Ib3, 0)","%9.3e");
                      // printF("TZ correction: ||uer||=%8.2e\n",max(fabs(uer)));
                      // printF("TZ correction: ||ues||=%8.2e\n",max(fabs(ues)));
                      // printF("TZ correction: ||uerr||=%8.2e\n",max(fabs(uerr)));
                      // printF("TZ correction: ||uers||=%8.2e\n",max(fabs(uers)));
                      // printF("TZ correction: ||uess||=%8.2e\n",max(fabs(uess)));
                      // printF("TZ correction: ||uerrr||=%8.2e\n",max(fabs(uerrr)));
                      // printF("TZ correction: ||uerrs||=%8.2e\n",max(fabs(uerrs)));
                      // printF("TZ correction: ||uerss||=%8.2e\n",max(fabs(uerss)));
                      // printF("TZ correction: ||uesss||=%8.2e\n",max(fabs(uesss)));
                      // printF("TZ correction: ||uerrrr||=%8.2e\n",max(fabs(uerrrr)));
                      // printF("TZ correction: ||uerrrs||=%8.2e\n",max(fabs(uerrrs)));
                      // printF("TZ correction: ||uerrss||=%8.2e\n",max(fabs(uerrss)));
                      // printF("TZ correction: ||uersss||=%8.2e\n",max(fabs(uersss)));
                      // printF("TZ correction: ||uessss||=%8.2e\n",max(fabs(uessss)));
                      // interfaceData.u(Ib1,Ib2,Ib3,n) += Sl*ue(Ib1,Ib2,Ib3); 
                      // interfaceData.u(Ib1,Ib2,Ib3,n) +=  cc(Ib1,Ib2,Ib3, 0)*ue(Ib1,Ib2,Ib3);
                                            interfaceData.u(Ib1,Ib2,Ib3,n) +=  cc(Ib1,Ib2,Ib3, 0)*ue(Ib1,Ib2,Ib3)
                                                                                                                +cc(Ib1,Ib2,Ib3, 1)*uer
                                                                                                                +cc(Ib1,Ib2,Ib3, 2)*ues
                                                                                                                +cc(Ib1,Ib2,Ib3, 3)*uerr
                                                                                                                +cc(Ib1,Ib2,Ib3, 4)*uers
                                                                                                                +cc(Ib1,Ib2,Ib3, 5)*uess
                                                                                                                +cc(Ib1,Ib2,Ib3, 6)*uerrr
                                                                                                                +cc(Ib1,Ib2,Ib3, 7)*uerrs
                                                                                                                +cc(Ib1,Ib2,Ib3, 8)*uerss
                                                                                                                +cc(Ib1,Ib2,Ib3, 9)*uesss
                                                                                                                +cc(Ib1,Ib2,Ib3,10)*uerrrr
                                                                                                                +cc(Ib1,Ib2,Ib3,11)*uerrrs
                                                                                                                +cc(Ib1,Ib2,Ib3,12)*uerrss
                                                                                                                +cc(Ib1,Ib2,Ib3,13)*uersss
                                                                                                                +cc(Ib1,Ib2,Ib3,14)*uessss;
                                        }
                                    }
                // OV_ABORT("Cgad::applyBoundaryConditionsForImplicitTimeStepping get RHS for CHAMP -- STOP HERE FOR NOW");  
                            }
                            else
                            {
                                    printF("**** applyBC: ADD TZ Correction to results from getData a0=%g, a1=%g **CHECK ME**\n",a0,a1);
                  // fix me for parallel
                                    OGFunction & e = *(parameters.dbase.get<OGFunction* >("exactSolution"));
                                    int rectangular=0;
                                    if( a0!=0. )
                                    {
                                        RealArray ue(Ib1,Ib2,Ib3,N);
                                        e.gd( ue,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,0,Ib1,Ib2,Ib3,N,t);
                                        interfaceData.u(Ib1,Ib2,Ib3,N) += a0*ue;
                                    }
                                    if( a1!=0. )
                                    {
                                        OV_GET_VERTEX_BOUNDARY_NORMAL(mg,side,axis,normal);
                    // #ifdef USE_PPP
                    //   const RealArray & normal = mg.vertexBoundaryNormalArray(side,axis);
                    // #else
                    //   const RealArray & normal = mg.vertexBoundaryNormal(side,axis);
                    // #endif
                                        RealArray uex(Ib1,Ib2,Ib3,N), uey(Ib1,Ib2,Ib3,N);
                                        e.gd( uex,xLocal,mg.numberOfDimensions(),rectangular,0,1,0,0,Ib1,Ib2,Ib3,N,t);
                                        e.gd( uey,xLocal,mg.numberOfDimensions(),rectangular,0,0,1,0,Ib1,Ib2,Ib3,N,t);
                                        if( mg.numberOfDimensions()==2 )
                                        {
                                            for( int n=N.getBase(); n<=N.getBound(); n++ )
                                                interfaceData.u(Ib1,Ib2,Ib3,n) += a1*( normal(Ib1,Ib2,Ib3,0)*uex(Ib1,Ib2,Ib3,n) + normal(Ib1,Ib2,Ib3,1)*uey(Ib1,Ib2,Ib3,n) );
                                        }
                                        else
                                        {
                                            RealArray uez(Ib1,Ib2,Ib3,N);
                                            e.gd( uez,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,1,Ib1,Ib2,Ib3,N,t);
                                            for( int n=N.getBase(); n<=N.getBound(); n++ )
                                                interfaceData.u(Ib1,Ib2,Ib3,n) += a1*( normal(Ib1,Ib2,Ib3,0)*uex(Ib1,Ib2,Ib3,n) + normal(Ib1,Ib2,Ib3,1)*uey(Ib1,Ib2,Ib3,n) + normal(Ib1,Ib2,Ib3,2)*uez(Ib1,Ib2,Ib3,n) );
                                        }
                                    }
                            }
        
                        }

                        if( debug() & 8 )
                        {
                            printP("AFTER REQUEST DATA FOR INTERFACE + TZ FIX (grid,side,axis)=(%d,%d,%d) t=%9.3e\n",grid,side,axis,t);
                            ::display(interfaceData.u(Ib1,Ib2,Ib3,N),"AFTER REQUEST DATA FOR INTERFACE + TZ FIX","%7.4f ");
                            ::display(bd,"bd: existing boundary data","%7.4f ");
                            printP("MAX DIFF getInterface - bd = %8.2e\n",max(fabs(bd(Ib1,Ib2,Ib3,N)-interfaceData.u(Ib1,Ib2,Ib3,N))));
                        }

                        if( true && interfaceCommunicationMode==Parameters::requestInterfaceDataWhenNeeded  )
                        {
                            if( multiDomainProblem )
                            {  // FILL IN RHS FOR implicit CHAMP condition
                                uLocal(Ig1,Ig2,Ig3,N) = interfaceData.u(Ib1,Ib2,Ib3,N);
                                bd(Ib1,Ib2,Ib3,N)     = interfaceData.u(Ib1,Ib2,Ib3,N);      // *wdh* March 25, 2022
                                if( debug() & 8 )
                                {
                                    ::display(uLocal(Ib1,Ib2,Ib3,N),"RHS ON BOUNDARY FOR CHAMP face uLocal(Ib1,Ib2,Ib3,N)","%7.4f ");
                                }
                            }
                            else
                            { // FOR TESTING IN A SINGLE DOMAIN 
                                bd(Ib1,Ib2,Ib3,N) = interfaceData.u(Ib1,Ib2,Ib3,N);
                            }
                        }
            // bd.redim(Ib1,Ib2,Ib3,N);
            // bd=interfaceData.u;

                    } // end if interfaceCommunicationMode==Parameters::requestInterfaceDataWhenNeeded

                } // end if interfaceCommunicationMode==Parameters::requestInterfaceDataWhenNeeded


                if( (bc==Parameters::neumannInterface || bc==Parameters::dirichletInterface ) && applyChampInterfaceConditions && ( multiDomainProblem == 0 || false ) )
                {
          // **** ASSIGN RHS FOR FAKE CHAMP CONDITIONS *****

                    getBoundaryIndex(mg.gridIndexRange(),side,axis,Ib1,Ib2,Ib3);
                    getGhostIndex(mg.gridIndexRange(),side,axis,Ig1,Ig2,Ig3);

          // OV_ABORT("CHAMP: RHS - FINISH ME ");

                    OV_GET_VERTEX_BOUNDARY_NORMAL(mg,side,axis,normal);
                    if( twilightZoneFlow )
                    {

                        OGFunction & e = *(parameters.dbase.get<OGFunction* >("exactSolution"));
                        RealArray ue(Ib1,Ib2,Ib3,N), uex(Ib1,Ib2,Ib3,N), uey(Ib1,Ib2,Ib3,N), uexx(Ib1,Ib2,Ib3,N), ueyy(Ib1,Ib2,Ib3,N);
                        int rectangular=0;
                        e.gd( ue  ,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,0,Ib1,Ib2,Ib3,N,t);
                        e.gd( uex ,xLocal,mg.numberOfDimensions(),rectangular,0,1,0,0,Ib1,Ib2,Ib3,N,t);
                        e.gd( uey ,xLocal,mg.numberOfDimensions(),rectangular,0,0,1,0,Ib1,Ib2,Ib3,N,t);
                        e.gd( uexx,xLocal,mg.numberOfDimensions(),rectangular,0,2,0,0,Ib1,Ib2,Ib3,N,t);
                        e.gd( ueyy,xLocal,mg.numberOfDimensions(),rectangular,0,0,2,0,Ib1,Ib2,Ib3,N,t);

                        const RealArray & champParameters = parameters.dbase.get<RealArray>("champParameters");
                        const Real pl    = champParameters(0,side,axis,grid);    // optimized Scwartz Parameter for side 1
                        const Real pr    = champParameters(1,side,axis,grid);    // optimized Scwartz Parameter for side 2
                        const Real theta = champParameters(2,side,axis,grid);    // K1/K2
                        const Real beta  = champParameters(3,side,axis,grid);    // D1/D2    
                        const Real Sl    = champParameters(4,side,axis,grid);    // pl/dxs; 

                        const Real dxs = dx[axis]; 

            // const real Sl    = champParameters(0,side,axis,grid);
            // const real theta = champParameters(2,side,axis,grid);
            // const real beta  = champParameters(3,side,axis,grid);  

                        printP("APPLY FAKE CHAMP BC's for TESTING Sl=%g, theta=%g, beta=%g *************************************************************\n",Sl,theta,beta);                  

                        const real a0 = Sl;
                        const real a1 = theta + Sl*dxs*theta;
            // ONLY VALID FOR CARTESIAN : axis==0 
                        const real a2 = dxs*( beta      ) + Sl*( .5*SQR(dxs)*beta      );
                        const real a3 = dxs*( (beta-1.) ) + Sl*( .5*SQR(dxs)*(beta-1.) );
                        const real a4 = a3; // for 3D 

                        int n=0;
                        assert( N.getLength()==1 );
                        uLocal(Ig1,Ig2,Ig3,n)= a0*ue + a1*( normal(Ib1,Ib2,Ib3,0)*uex + normal(Ib1,Ib2,Ib3,1)*uey ) + a2*uexx + a3*ueyy;

                    }
                }
                else if( bc==Parameters::dirichletInterface  && !applyChampInterfaceConditions )
                {
          // DIRICHLET INTERFACE (NO CHAMP)
                    getBoundaryIndex(mg.gridIndexRange(),side,axis,Ib1,Ib2,Ib3);
                    const int includeGhost=1;
                    bool ok = ParallelUtility::getLocalArrayBounds(u,uLocal,Ib1,Ib2,Ib3,includeGhost);
                    if( !ok ) continue;
                        
                    if( debug() & 4 )
                        printP("** Cgad:applyImpBC: setting RHS for a dirichlet interface bc(%i,%i,%i)=%i\n",side,axis,grid,bc);
          // u.applyBoundaryCondition(N,dirichlet,BCTypes::boundary(side,axis),bcData,pBoundaryData,t,
          //                       Overture::defaultBoundaryConditionParameters(),grid);

                    RealArray & bd = parameters.getBoundaryData(side,axis,grid,mg);
                    
                    if( debug() & 4 )
                    {
                        ::display(uLocal(Ib1,Ib2,Ib3,N),"RHS for dirichlet BC: uLocal(Ib1,Ib2,Ib3,N)","%7.4f ");

                        if( twilightZoneFlow )
                        {
                            
                            OGFunction & e = *(parameters.dbase.get<OGFunction* >("exactSolution"));
                            RealArray ue(Ib1,Ib2,Ib3,N);
                            int rectangular=0;
                            e.gd( ue,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,0,Ib1,Ib2,Ib3,N,t);

                            ::display(ue(Ib1,Ib2,Ib3,N),"Cgad:applyImpBC: EXACT RHS for dirichlet BC: ue(Ib1,Ib2,Ib3,N)","%7.4f ");

                            ::display(fabs(bd(Ib1,Ib2,Ib3,N)-ue(Ib1,Ib2,Ib3,N)),
                                                      "Cgad:applyImpBC: ERROR in RHS for Dirichlet BC: bd-ue","%7.4f ");


              // printP("Cgad:applyImpBC: set RHS to exact\n");
              // bd(Ib1,Ib2,Ib3,N)=ue(Ib1,Ib2,Ib3,N);  // *****************************************************************

                        }
                    }

                    uLocal(Ib1,Ib2,Ib3,N)=bd(Ib1,Ib2,Ib3,N);
                    
                }
                else if( bc==neumannBoundaryCondition || 
                                  bc==AdParameters::mixedBoundaryCondition ||
                                  ( bc==Parameters::neumannInterface && !applyChampInterfaceConditions ) ||
                                  bc==Parameters::axisymmetric )
                {
          // ****** neumann boundary condition or neumann interface (but not CHAMP) *****
                    getBoundaryIndex(mg.gridIndexRange(),side,axis,Ib1,Ib2,Ib3);
                    getGhostIndex(mg.gridIndexRange(),side,axis,Ig1,Ig2,Ig3);

                    #ifdef USE_PPP
                      RealArray uLocal; getLocalArrayWithGhostBoundaries(u,uLocal);
                      const RealArray & normal = mg.vertexBoundaryNormalArray(side,axis);
                    #else
                      const RealArray & normal = mg.vertexBoundaryNormal(side,axis);
                    #endif

                    const int includeGhost=1;
                    bool ok = ParallelUtility::getLocalArrayBounds(u,uLocal,Ib1,Ib2,Ib3,includeGhost);
                    ok = ParallelUtility::getLocalArrayBounds(u,uLocal,Ig1,Ig2,Ig3,includeGhost);
                    if( !ok ) continue;
                        
                    if( bc==Parameters::neumannInterface )
                    {
            // The RHS for the interface equations is placed in the BoundaryDataArray
            // (This includes TZ forcing)
                        if( debug() & 4 )
                            printP("** Cgad:applyImpBC: setting RHS for a neumann interface bc(%i,%i,%i)=%i\n",side,axis,grid,bc);
                        RealArray & bd = parameters.getBoundaryData(side,axis,grid,mg);




                        if( debug() & 4 && twilightZoneFlow )
                        {
                            ::display(bd(Ib1,Ib2,Ib3,N),"Cgad:applyImpBC: RHS for Neumann BC: bd(Ib1,Ib2,Ib3,N)","%7.4f ");

                            OGFunction & e = *(parameters.dbase.get<OGFunction* >("exactSolution"));
                            RealArray ue(Ib1,Ib2,Ib3,N), uex(Ib1,Ib2,Ib3,N), uey(Ib1,Ib2,Ib3,N);
                            int rectangular=0;
                            e.gd( ue ,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,0,Ib1,Ib2,Ib3,N,t);
                            e.gd( uex,xLocal,mg.numberOfDimensions(),rectangular,0,1,0,0,Ib1,Ib2,Ib3,N,t);
                            e.gd( uey,xLocal,mg.numberOfDimensions(),rectangular,0,0,1,0,Ib1,Ib2,Ib3,N,t);

                            if( mg.numberOfDimensions()==2 )
                            {
                                for( int n=N.getBase(); n<=N.getBound(); n++ )
                                {
                                    real a0=mixedCoeff(n,side,axis,grid), a1=mixedNormalCoeff(n,side,axis,grid);
                                    ue(Ib1,Ib2,Ib3,n)=a1*(uex(Ib1,Ib2,Ib3,n)*normal(Ib1,Ib2,Ib3,0)+
                                                                                uey(Ib1,Ib2,Ib3,n)*normal(Ib1,Ib2,Ib3,1))+ a0*ue(Ib1,Ib2,Ib3,n);
                                    ::display(ue(Ib1,Ib2,Ib3,n),"Cgad:applyImpBC: EXACT RHS for Neumann BC: bd(Ib1,Ib2,Ib3,N)","%7.4f ");
                                    ::display(fabs(bd(Ib1,Ib2,Ib3,n)-ue(Ib1,Ib2,Ib3,n)),
                                                      "Cgad:applyImpBC: ERROR in RHS for Neumann BC: bd(Ib1,Ib2,Ib3,N)","%7.4f ");

                  // ** bd(Ib1,Ib2,Ib3,N) *= 1./a1;
                                    
                                }
                            }
              // printP("Cgad:applyImpBC: set RHS to exact\n");
              // bd(Ib1,Ib2,Ib3,N)=ue(Ib1,Ib2,Ib3,N);  // ***********************************************************
                        }
                        
                        uLocal(Ig1,Ig2,Ig3,N)=bd(Ib1,Ib2,Ib3,N);

            // If the adjacent face is a dirichlet BC then we change the right-hand-side on the extended
            // boundary point: 
            // 
            //   bc=dirichlet
            //       ------c--o   <-  adjust RHS to neuman interface here 
            //             |
            //       ------+--+   <- interface side 
            //             |
            // Since the interior equation is not applied at the corner pt "c" we can not just impose
            // [ k T.n ]=0 since we need another equation to determine the 2 ghost point values (there is
            // another ghost pt value on the other side of the interface ). Therefore we set 
            //  k*T.n = given at point "c" 

            // loop over adjacent sides
                        for( int dir=1; dir<mg.numberOfDimensions(); dir++ ) for( int side2=0; side2<=1; side2++ )
                        {
                            int dir2 = (axis+dir) % mg.numberOfDimensions();
                            if( mg.boundaryCondition(side2,dir2)==dirichletBoundaryCondition )
                            {
                                Jb1=Ib1, Jb2=Ib2, Jb3=Ib3;
                                Jg1=Ig1, Jg2=Ig2, Jg3=Ig3;
                // check for parallel: 
                                if( Jbv[dir2].getBase()  > mg.gridIndexRange(side2,dir2) ||
                                        Jbv[dir2].getBound() < mg.gridIndexRange(side2,dir2) )
                                {  
                                    ok=false;
                                    continue;
                                }
                                Jbv[dir2]=mg.gridIndexRange(side2,dir2);
                                Jgv[dir2]=mg.gridIndexRange(side2,dir2);

                                if( debug() & 8 )
                                {
                                    printP("Cgad:applyImpBC: set RHS to exact where interface meets adj-dirichlet "
                                                  " (side,axis)=(%i,%i) (side2,dir2)=(%i,%i) Jv=[%i,%i][%i,%i][%i,%i]\n",
                                                  side,axis,side2,dir2,Jb1.getBase(),Jb1.getBound(),Jb2.getBase(),Jb2.getBound(),
                                                  Jb3.getBase(),Jb3.getBound() );
                                }
                                
                                if( twilightZoneFlow )
                                {
                                    OGFunction & e = *(parameters.dbase.get<OGFunction* >("exactSolution"));
                                    RealArray ue(Jb1,Jb2,Jb3,N), uex(Jb1,Jb2,Jb3,N), uey(Jb1,Jb2,Jb3,N);
                                    int rectangular=0;
                                    e.gd( ue ,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,0,Jb1,Jb2,Jb3,N,t);
                                    e.gd( uex,xLocal,mg.numberOfDimensions(),rectangular,0,1,0,0,Jb1,Jb2,Jb3,N,t);
                                    e.gd( uey,xLocal,mg.numberOfDimensions(),rectangular,0,0,1,0,Jb1,Jb2,Jb3,N,t);

                                    if( mg.numberOfDimensions()==2 )
                                    {
                                        for( int n=N.getBase(); n<=N.getBound(); n++ )
                                        {
                                            real a0=mixedCoeff(n,side,axis,grid), a1=mixedNormalCoeff(n,side,axis,grid);
                                            ue(Jb1,Jb2,Jb3,n)=a1*(uex(Jb1,Jb2,Jb3,n)*normal(Jb1,Jb2,Jb3,0)+
                                                                                        uey(Jb1,Jb2,Jb3,n)*normal(Jb1,Jb2,Jb3,1))+ a0*ue(Jb1,Jb2,Jb3,n);
                                        }
                                        uLocal(Jg1,Jg2,Jg3,N)=ue(Jb1,Jb2,Jb3,N);
                                    }
                                    else 
                                    {
                                        RealArray uez(Jb1,Jb2,Jb3,N);
                                        e.gd( uez,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,1,Jb1,Jb2,Jb3,N,t);
                                        for( int n=N.getBase(); n<=N.getBound(); n++ )
                                        {
                                            real a0=mixedCoeff(n,side,axis,grid), a1=mixedNormalCoeff(n,side,axis,grid);
                                            ue(Jb1,Jb2,Jb3,n)=a1*(uex(Jb1,Jb2,Jb3,n)*normal(Jb1,Jb2,Jb3,0)+
                                                                                        uey(Jb1,Jb2,Jb3,n)*normal(Jb1,Jb2,Jb3,1)+
                                                                                        uez(Jb1,Jb2,Jb3,n)*normal(Jb1,Jb2,Jb3,2))+ a0*ue(Jb1,Jb2,Jb3,n);
                                        }
                                        uLocal(Jg1,Jg2,Jg3,N)=ue(Jb1,Jb2,Jb3,N);
                                    }
                                    
                                }
                                else
                                {
                  //MappedGridOperators & op = *(u.getOperators());
                  //op.derivative(MappedGridOperators::xDerivative,uLocal,ux  ,Jb1,Jb2,Jb3,N);
                  //op.derivative(MappedGridOperators::yDerivative,uLocal,uy  ,Jb1,Jb2,Jb3,N);
                  //op.derivative(MappedGridOperators::zDerivative,uLocal,uz  ,Jb1,Jb2,Jb3,N);

                  // Here we assume the Dirichlet BC is a constant value so that the normal derivative  *fix me*
                  // of the solution along the Dirichlet BC is zero: 
                  // *NOTE* if we fix this for the case of a variable Dirichlet BC then we must also adjust the 
                  //   interface getRHS to eval   k*u.n - k*ue.n so that the residuals in the interface equations 
                  //   will go to zero. 
                                    for( int n=N.getBase(); n<=N.getBound(); n++ )
                                    {
                                        real a0=mixedCoeff(n,side,axis,grid), a1=mixedNormalCoeff(n,side,axis,grid);
                    // note: uLocal(Jb1,Jb2,Jb3,n) has been set above by the Dirichlet BC 
                                        uLocal(Jg1,Jg2,Jg3,n)=a0*uLocal(Jb1,Jb2,Jb3,n);
                                    }
                                        
                                }
                            }
                        }
                            


                    }
                    else if( !twilightZoneFlow )
                    {
                        if( bc==neumannBoundaryCondition || 
                                bc==AdParameters::mixedBoundaryCondition )
                        {
                            for( int n=N.getBase(); n<=N.getBound(); n++ )
                            {
                                uLocal(Ig1,Ig2,Ig3,n)=mixedRHS(n,side,axis,grid);  // set ghost value to RHS for neumann oe mixed BC
                            }
                        }
                        else
                        {
                            uLocal(Ig1,Ig2,Ig3,N)=0.;                          // axisymmetric
                        }
                    }
                    else
                    {
            // Twilight-zone forcing:

                        OGFunction & e = *(parameters.dbase.get<OGFunction* >("exactSolution"));
                        RealArray ue(Ib1,Ib2,Ib3,N), uex(Ib1,Ib2,Ib3,N), uey(Ib1,Ib2,Ib3,N);
                        int rectangular=0;
                        e.gd( ue ,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,0,Ib1,Ib2,Ib3,N,t);
                        e.gd( uex,xLocal,mg.numberOfDimensions(),rectangular,0,1,0,0,Ib1,Ib2,Ib3,N,t);
                        e.gd( uey,xLocal,mg.numberOfDimensions(),rectangular,0,0,1,0,Ib1,Ib2,Ib3,N,t);

                        if( mg.numberOfDimensions()==2 )
                        {
                            for( int n=N.getBase(); n<=N.getBound(); n++ )
                            {
                                real a0=mixedCoeff(n,side,axis,grid), a1=mixedNormalCoeff(n,side,axis,grid);
                // printF(" ***** BC: a0=%8.2e a1=%8.2e\n",a0,a1);
                                
                                uLocal(Ig1,Ig2,Ig3,n)=a1*(uex(Ib1,Ib2,Ib3,n)*normal(Ib1,Ib2,Ib3,0)+
                                                                                    uey(Ib1,Ib2,Ib3,n)*normal(Ib1,Ib2,Ib3,1)) + a0*ue(Ib1,Ib2,Ib3,n);
                            }
                        }
                        else
                        {
                            RealArray uez(Ib1,Ib2,Ib3,N);
                            e.gd( uez,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,1,Ib1,Ib2,Ib3,N,t);
                            for( int n=N.getBase(); n<=N.getBound(); n++ )
                            {
                                real a0=mixedCoeff(n,side,axis,grid), a1=mixedNormalCoeff(n,side,axis,grid);
                                uLocal(Ig1,Ig2,Ig3,n)=a1*(uex(Ib1,Ib2,Ib3,n)*normal(Ib1,Ib2,Ib3,0)+
                                                                                    uey(Ib1,Ib2,Ib3,n)*normal(Ib1,Ib2,Ib3,1)+
                                                                                    uez(Ib1,Ib2,Ib3,n)*normal(Ib1,Ib2,Ib3,2)) + a0*ue(Ib1,Ib2,Ib3,n);
                            }
                        }

                    }
                }
            }
        } // end for axis
                
    // ************ try this for interfaces ********* 080909
        u.updateGhostBoundaries();


      if( debug() & 8 )
      {
          display(uLocal,"RHS AFTER FILLING IN BCs uLocal(Ib1,Ib2,Ib3,N)","%7.4f ");
      }

    } // end if grid is implicit
    

    return 0;


}

