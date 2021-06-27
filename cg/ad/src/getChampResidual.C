// This file automatically generated from getChampResidual.bC with bpp.
// ================================================================================================
//  Evaluate the resdiual in the CHAMP interface conditions
// ================================================================================================

#include "Cgad.h"
#include "CompositeGridOperators.h"
#include "ParallelUtility.h"
#include "ParallelGridUtility.h"
#include "Interface.h"  

// #include "Oges.h"

// #include "SparseRep.h"

// #include "gridFunctionNorms.h"

#define ForBoundary(side,axis)   for( int axis=0; axis<mg.numberOfDimensions(); axis++ ) for( int side=0; side<=1; side++ )

#define champResidualOpt EXTERN_C_NAME(champresidualopt)
extern "C"
{
// subroutine champResidualOpt( nd,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,gridIndexRange,dimRange,isPeriodic,boundaryCondition,mask,xy,rsxy,u,f,ipar,rpar,ierr )
void champResidualOpt( const int & nd,
                                              const int & nd1a,const int & nd1b,const int & nd2a,const int & nd2b,const int & nd3a,const int & nd3b,const int & nd4a,const int & nd4b,
                                              const int & md1a,const int & md1b,const int & md2a,const int & md2b,const int & md3a,const int & md3b,const int & md4a,const int & md4b,
                                              const int & ld1a,const int & ld1b,const int & ld2a,const int & ld2b,const int & ld3a,const int & ld3b,const int & ld4a,const int & ld4b,
                                              const int & gridIndexRange,const int & dimRange,const int & isPeriodic,const int & boundaryCondition,
                                              const int & mask,const real & xy,const real & rsxy,const real & u,const real & f, const real & coeff, const int & ipar,const real & rpar,const int & ierr );

}

// ====================================================================================
// Macro: Add twilght-zone corrections for CHAMP
// ====================================================================================

// // =========================================================================================
// // Macro: Add twilght-zone corrections for CHAMP  **THIS IS COPIED FROM applyBoundaryConditions.bC ****** FIX ME ****
// // =========================================================================================
// #beginMacro addTwilightZoneCorrectionForChamp()              

//   OGFunction & e = *(parameters.dbase.get<OGFunction* >("exactSolution"));
//   RealArray ue(Ib1,Ib2,Ib3,N), uex(Ib1,Ib2,Ib3,N), uey(Ib1,Ib2,Ib3,N), uexx(Ib1,Ib2,Ib3,N), ueyy(Ib1,Ib2,Ib3,N);
//   int rectangular=0;
//   e.gd( ue  ,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,0,Ib1,Ib2,Ib3,N,t);
//   e.gd( uex ,xLocal,mg.numberOfDimensions(),rectangular,0,1,0,0,Ib1,Ib2,Ib3,N,t);
//   e.gd( uey ,xLocal,mg.numberOfDimensions(),rectangular,0,0,1,0,Ib1,Ib2,Ib3,N,t);
//   e.gd( uexx,xLocal,mg.numberOfDimensions(),rectangular,0,2,0,0,Ib1,Ib2,Ib3,N,t);
//   e.gd( ueyy,xLocal,mg.numberOfDimensions(),rectangular,0,0,2,0,Ib1,Ib2,Ib3,N,t);

//   const RealArray & champParameters = parameters.dbase.get<RealArray>("champParameters");
//   // const real Sl    = champParameters(0,side,axis,grid);
//   // const real theta = champParameters(2,side,axis,grid);
//   // const real beta  = champParameters(3,side,axis,grid);  
//   const Real pl    = champParameters(0,side,axis,grid);    // optimized Scwartz Parameter for side 1
//   const Real pr    = champParameters(1,side,axis,grid);    // optimized Scwartz Parameter for side 2
//   const Real theta = champParameters(2,side,axis,grid);    // K1/K2
//   const Real beta  = champParameters(3,side,axis,grid);    // D1/D2    
//   const Real Sl    = champParameters(4,side,axis,grid);    // pl/dxs; 

//   const Real dxs = dx[axis]; 

//   printF("**** ADD TZ Correction to results from getData for CHAMP:  Sl=%g, theta=%g, beta=%g *********************************\n",Sl,theta,beta);                  

//   const real a0 = Sl;
//   const real a1 = theta + Sl*dxs*theta;
//   // ONLY VALID FOR CARTESIAN : axis==0 
//   assert( isRectangular );

//   const real a2 = dxs*( beta      ) + Sl*( .5*SQR(dxs)*beta      );
//   const real a3 = dxs*( (beta-1.) ) + Sl*( .5*SQR(dxs)*(beta-1.) );
//   const real a4 = a3; // for 3D 

//   OV_GET_VERTEX_BOUNDARY_NORMAL(mg,side,axis,normal);

//   int n=0;
//   assert( N.getLength()==1 );
//   interfaceData.u(Ib1,Ib2,Ib3,n) += a0*ue + a1*( normal(Ib1,Ib2,Ib3,0)*uex + normal(Ib1,Ib2,Ib3,1)*uey ) + a2*uexx + a3*ueyy;
// #endMacro

// =========================================================================================
// Macro: Call optimized residual routine
// =========================================================================================

// ================================================================================================
/// \brief Evaluate the resdiual in the CHAMP interface conditions
/// \return the maximum residual
// ================================================================================================
real
getChampResidual( realCompositeGridFunction & u, Parameters & parameters, Real t, Real dt )
{

    CompositeGrid & cg = *u.getCompositeGrid();
    const int numberOfDimensions = cg.numberOfDimensions();
    int numberOfComponents = 1;
    Range N = parameters.dbase.get<Range >("Rt");   // time dependent variables
    
  // const int eq1=0, eq2=1;   // equation numbers
  // const int uc=0, vc=1;     // component numbers
    
  // const std::vector<real> & kappa  = parameters.dbase.get<std::vector<real> >("kappa");
  // const real & thermalConductivity = parameters.dbase.get<real>("thermalConductivity");

    const int & multiDomainProblem = parameters.dbase.get<int>("multiDomainProblem"); 
    const bool twilightZoneFlow    = parameters.dbase.get<bool>("twilightZoneFlow");
    const int & debug              = parameters.dbase.get<int>("debug");
    const int & myid               = parameters.dbase.get<int>("myid");
    const int & orderOfAccuracy    = parameters.dbase.get<int>("orderOfAccuracy");
    assert( orderOfAccuracy==2 || orderOfAccuracy==4 );



    RealArray & champParameters = parameters.dbase.get<RealArray>("champParameters");

    const IntegerArray & interfaceType = parameters.dbase.get<IntegerArray >("interfaceType");

  // printF("\n -- CGAD-- champBoundaryConditions: multiDomainProblem=%d,  kappa=%g, thermalConductivity=%g\n",multiDomainProblem,kappa[0],thermalConductivity);

    
    Index Ibv[3], &Ib1=Ibv[0], &Ib2=Ibv[1], &Ib3=Ibv[2];
    Index Igv[3], &Ig1=Igv[0], &Ig2=Igv[1], &Ig3=Igv[2];

    int i1,i2,i3, j1,j2,j3, i1m,i2m,i3m, m1,m2,m3;
    int isv[3], &is1=isv[0], &is2=isv[1], &is3=isv[2];

    real maxInterfaceResidual=0.; // max residual over all interfaces 
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {

        MappedGrid & mg = cg[grid];
        OV_GET_SERIAL_ARRAY(real,u[grid],uLocal);

        const int isRectangular=mg.isRectangular();

        real dx[3]={1.,1.,1.};
        if( isRectangular )
            mg.getDeltaX(dx);
        else
            mg.update(MappedGrid::THEinverseVertexDerivative);

    // if( !isRectangular )
    // {
    //   OV_ABORT("FINISH ME");
    // }
        

    // --- loop over faces ----      
        ForBoundary(side,axis)
        {
            if( interfaceType(side,axis,grid) == Parameters::heatFluxInterface && mg.boundaryCondition(side,axis)>0 )
            {
        // printF("+++++ getChampResidual for (grid,side,axis)=(%d,%d,%d) (heatFluxInterface)\n",grid,side,axis);

        // Set the index-shift for this side
                is1=is2=is3=0;
                isv[axis]=1-2*side;


        // const Real dxs = dx[axis]; 

        // const Real pl    = champParameters(0,side,axis,grid);    // optimized Scwartz Parameter for side 1
        // const Real pr    = champParameters(1,side,axis,grid);    // optimized Scwartz Parameter for side 2
                const Real theta = champParameters(2,side,axis,grid);    // K1/K2
                const Real beta  = champParameters(3,side,axis,grid);    // D1/D2    
                const Real Sl    = champParameters(4,side,axis,grid); 



                OV_GET_VERTEX_BOUNDARY_NORMAL(mg,side,axis,normal); 

                InterfaceData interfaceData;
        // Range Rx=numberOfDimensions;
                getBoundaryIndex(mg.gridIndexRange(),side,axis,Ib1,Ib2,Ib3);
                getGhostIndex(   mg.gridIndexRange(),side,axis,Ig1,Ig2,Ig3);

                interfaceData.u.redim(Ib1,Ib2,Ib3,N); // heat flux is returned here 
                interfaceData.t=t;
                interfaceData.u=0;
        // ::display(interfaceData.u(Ib1,Ib2,Ib3,N),"BEFORE REQUEST DATA FOR INTERFACE","%7.4f ");

        // // Boundary data should go here
        // RealArray & bd = parameters.getBoundaryData(side,axis,grid,mg);
        // // bd=0;

                if( debug & 8 )
                {
                    printF("\n");
                    printF(" >>>> getChampResidual:: REQUEST INTERFACE DATA for (grid,side,axis)=(%d,%d,%d) t=%9.3e\n",grid,side,axis,t);
                }
                        
        // We pass the coefficients in the dirichlet, neumann or mixed BC
        // *** THIS IS NOT CORRECT -- FIX ME ***

        // From ins/src/addedMassImplicitBoundaryConditions.bC: 
        //Parameters & bulkSolidParams = getInterfaceParameters( grid,side,axis,parameters );

                GridFaceDescriptor & myGFD = getInterfaceGridFaceDescriptor( grid, side, axis, parameters );
        // real *sourceMixedBC = myGFD.dbase.get<real[2]>("sourceMixedBC");

        // if( debug & 8 )
        //   printF("getChampResidual:: After getInterfaceGridFaceDescriptor: myGFD.a=[%g,%g], sourceMixedBC=[%g,%g]\n",myGFD.a[0],myGFD.a[1],sourceMixedBC[0],sourceMixedBC[1]);



                GridFaceDescriptor targetInfo(-1,grid,side,axis);

        // const int n=0; // component -- FIX ME
        // const real a0=mixedCoeff(n,side,axis,grid), a1=mixedNormalCoeff(n,side,axis,grid);
        // // targetInfo.a[0]=a0; targetInfo.a[1]=-a1; targetInfo.a[2]=0.;  // **NOTE*** flip sign of a1 for normal for other side, *check me*

        // targetInfo.a[0]=sourceMixedBC[0]; targetInfo.a[1]=sourceMixedBC[1]; targetInfo.a[2]=0.; 
        // // const real a0 = targetInfo.a[0], a1= targetInfo.a[1]; 

  
        // CHAMP CONDITIONS
        // Evaluate   S*uR(dx) + n.grad( uR(dx) ) 
  
                targetInfo.a[0] =  Sl;  // optimized Schwartz parameter for this domain ("left")
                targetInfo.a[1] = -1.;  // flip sign of normal
        // Do this for now: 
                targetInfo.dbase.put<int>("getChampData")=1;

                int interfaceDataOptions = Parameters::heatFluxInterfaceData;
                bool saveTimeHistory=true;

                if( debug & 8 )
                {
                    printF("\n");
                    printF(" >>>> getChampResidual:: REQUEST INTERFACE DATA for (grid,side,axis)=(%d,%d,%d) t=%9.3e. target: (a0,a1)=(%g,%g)\n",
                                  grid,side,axis,t,targetInfo.a[0],targetInfo.a[1]);
                }

                getInterfaceData( t, grid, side, axis, 
                                                    interfaceDataOptions,
                                                    interfaceData.u,
                                                    parameters,saveTimeHistory,
                                                    &targetInfo );

                if( twilightZoneFlow  )
                {
                    OV_GET_SERIAL_ARRAY(real,mg.vertex(),xLocal);    
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
                        printF("**** ADD TZ Correction to results from getData for CHAMP:  Sl=%g, theta=%g, beta=%g *********************************\n",Sl,theta,beta);                  
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
              // ONLY VALID FOR CARTESIAN AND  axis==0 
                            assert( axis==0 );
                            const Real dxs = dx[axis]; 
                            const real a0 = Sl;
                            const real a1 = theta + Sl*dxs*theta;
                            const real a2 = dxs*( beta      ) + Sl*( .5*SQR(dxs)*beta      );
                            const real a3 = dxs*( (beta-1.) ) + Sl*( .5*SQR(dxs)*(beta-1.) );
                            const real a4 = a3; // for 3D 
                            int n=0;
                            assert( N.getLength()==1 );
                            interfaceData.u(Ib1,Ib2,Ib3,n) += a0*ue + a1*( normal(Ib1,Ib2,Ib3,0)*uex + normal(Ib1,Ib2,Ib3,1)*uey ) + a2*uexx + a3*ueyy;
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
                            printF("**** ADD TZ Correction to results from getData for CHAMP:  Sl=%g, theta=%g, beta=%g **************** FINISH ME *****************\n",Sl,theta,beta); 
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
                            RealArray uep(Ib1+js1,Ib2+js2,Ib3+js3,N), uem(Ib1-js1,Ib2-js2,Ib3-js3,N);
                            e.gd( uep,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,0,Ib1+js1,Ib2+js2,Ib3+js3,N,t);
                            e.gd( uem,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,0,Ib1-js1,Ib2-js2,Ib3-js3,N,t);
                            RealArray uer2(Ib1,Ib2,Ib3,N), uer2r2(Ib1,Ib2,Ib3,N);
                            uer2   = ( uep - uem )*(1./(2.*dr2));
                            uer2r2 = ( uep -2.*ue + uem )*(1./(dr2*dr2));
              // -- compute the mixed derivative w.r.t r1 and r2 ----
                            RealArray uepp(Ib1+1,Ib2+1,Ib3,N), uemm(Ib1-1,Ib2-1,Ib3,N), uepm(Ib1+1,Ib2-1,Ib3,N), uemp(Ib1-1,Ib2+1,Ib3,N);
                            e.gd( uemm,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,0,Ib1-1,Ib2-1,Ib3,N,t);
                            e.gd( uepm,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,0,Ib1+1,Ib2-1,Ib3,N,t);
                            e.gd( uemp,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,0,Ib1-1,Ib2+1,Ib3,N,t);
                            e.gd( uepp,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,0,Ib1+1,Ib2+1,Ib3,N,t);
                            RealArray uer1r2(Ib1,Ib2,Ib3,N);
                            uer1r2 = ( uepp - uemp -  uepm + uemm )*(1./(4.*dr1*dr2));   // D0r1 * D0r2 
                            int n=0;
                            assert( N.getLength()==1 );
              // interfaceData.u(Ib1,Ib2,Ib3,n) += a0*ue + a1*( normal(Ib1,Ib2,Ib3,0)*uex + normal(Ib1,Ib2,Ib3,1)*uey ) + a2*uexx + a3*ueyy;
              // Some coefficients in the CHAMP conditions were saved when the matrix was created (champBoundaryConditions.bC)
                            RealArray & cc = myGFD.dbase.get<RealArray>("ccChamp");
              // RealArray cc(Ib1,Ib2,Ib3,6); // need to save these when matrix is formed.  **FIX ME***
              // cc=0.;
              // cc(Ib1,Ib2,Ib3,1) = beta; 
              // cc(Ib1,Ib2,Ib3,3) = -1.; 
              // OV_ABORT("finish me");
                            RealArray Nc(Ib1,Ib2,Ib3,N);
                            Nc = theta*( normal(Ib1,Ib2,Ib3,0)*uex + normal(Ib1,Ib2,Ib3,1)*uey ) + cc(Ib1,Ib2,Ib3,0)*uer2; 
                            interfaceData.u(Ib1,Ib2,Ib3,n) += Sl*ue   
                                                                                                + cc(Ib1,Ib2,Ib3,1)*( uexx + ueyy ) 
                                                                                                + cc(Ib1,Ib2,Ib3,2)*uer1r2 
                                                                                                + cc(Ib1,Ib2,Ib3,3)*uer2r2
                                                                                                + cc(Ib1,Ib2,Ib3,4)*uer2
                                                                                                + cc(Ib1,Ib2,Ib3,5)*Nc;
                        }
          // OV_ABORT("Cgad::applyBoundaryConditionsForImplicitTimeStepping get RHS for CHAMP -- STOP HERE FOR NOW");  
  
                }

        // --- evaluate the residual -----
                    int option = 0;    // eval residual
                    int gridType = isRectangular ? 0 : 1;
                    int tc=0; // first component
                    int ipar[] = {
                        option              ,            // ipar( 0)
                        side                ,            // ipar( 1)
                        axis                ,            // ipar( 2)
                        grid                ,            // ipar( 3)
                        gridType            ,            // ipar( 4)
                        orderOfAccuracy     ,            // ipar( 5)
                        twilightZoneFlow    ,            // ipar( 6)
                        tc                  ,            // ipar( 7)
                        debug               ,            // ipar( 8)
                        myid                             // ipar( 9)
                                              };
         // option             = ipar( 0)
         // side               = ipar( 1)
         // axis               = ipar( 2)
         // grid               = ipar( 3)
         // gridType           = ipar( 4)
         // orderOfAccuracy    = ipar( 5)
         // twilightZone       = ipar( 6)
         // tc                 = ipar( 7)
         // debug              = ipar( 8)
         // myid               = ipar( 9)
                    real maxRes=0., l2Res=0.; 
                    real rpar[] = {
                        t                , //  rpar( 0)
                        dt               , //  rpar( 1)
                        dx[0]            , //  rpar( 2)
                        dx[1]            , //  rpar( 3)
                        dx[2]            , //  rpar( 4)
                        mg.gridSpacing(0), //  rpar( 5)
                        mg.gridSpacing(1), //  rpar( 6)
                        mg.gridSpacing(2), //  rpar( 7)
                        (real &)(parameters.dbase.get<OGFunction*>("exactSolution")) ,        //  rpar( 8) ! pointer for exact solution -- new : 110311 
                        theta,             //  rpar( 9)
                        beta,              //  rpar(10)
                        Sl,                //  rpar(11)
                        REAL_MIN,          //  rpar(12)
                        maxRes,            //  rpar(13)  ... returned 
                        l2Res              //  rpar(14)  ... returned
                                                };
         // t        = rpar( 0)
         // dt       = rpar( 1)
         // dx(0)    = rpar( 2)
         // dx(1)    = rpar( 3)
         // dx(2)    = rpar( 4)
         // dr(0)    = rpar( 5)
         // dr(1)    = rpar( 6)
         // dr(2)    = rpar( 7)
         // ep       = rpar( 8)  ! for TZ
         // theta    = rpar( 9)  ! ratio of thermal conductivities 
         // beta     = rpar(10)  ! ratio od thermal diffusivities 
         // Sl       = rpar(11)  ! Optimized Schwartz parameter for the current side
         // REAL_MIN = rpar(12)
                    IntegerArray indexRangeLocal(2,3), dimLocal(2,3), bcLocal(2,3);
                    ParallelGridUtility::getLocalIndexBoundsAndBoundaryConditions( u[grid],indexRangeLocal,dimLocal,bcLocal );                
                    OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
                    int *pmask = maskLocal.getDataPointer();
                    real temp, *pxy=&temp, *prsxy=&temp;
                    if( !isRectangular )
                    {
                        OV_GET_SERIAL_ARRAY(real,mg.inverseVertexDerivative(),rxLocal);
                        prsxy = rxLocal.getDataPointer();
                    }
                    bool vertexNeeded = twilightZoneFlow;
                    if( vertexNeeded )
                    {
                        OV_GET_SERIAL_ARRAY(real,mg.vertex(),xLocal);
                        pxy=xLocal.getDataPointer();
                    }
                    real *pu = uLocal.getDataPointer();
                    RealArray & f = interfaceData.u;  // Hold RHS of champ condition for residual computation
                    real *pf = f.getDataPointer();
                    if( false )
                        ::display(f,"f: RHS TO CHAMP from getData","%8.2e ");
          // -- Find the copy of the CHAMP coefficient matrix for computing the residual in champResidualOpt.bf90 ---
          // bool sameSide=true; 
          // GridFaceDescriptor & myGFD = getInterfaceGridFaceDescriptor( grid, side, axis, parameters, sameSide );
                    RealArray & coeffChamp = myGFD.dbase.get<RealArray>("coeffChamp");
                    int ierr=0;
                    champResidualOpt( numberOfDimensions,
                                                        uLocal.getBase(0),uLocal.getBound(0),uLocal.getBase(1),uLocal.getBound(1),uLocal.getBase(2),uLocal.getBound(2),uLocal.getBase(3),uLocal.getBound(3),
                                                        f.getBase(0),f.getBound(0),f.getBase(1),f.getBound(1),f.getBase(2),f.getBound(2),f.getBase(3),f.getBound(3),
                                                        coeffChamp.getBase(0),coeffChamp.getBound(0),coeffChamp.getBase(1),coeffChamp.getBound(1),coeffChamp.getBase(2),coeffChamp.getBound(2),coeffChamp.getBase(3),coeffChamp.getBound(3),
                                                        indexRangeLocal(0,0), dimLocal(0,0), mg.isPeriodic(0),bcLocal(0,0), 
                                                        *pmask, *pxy, *prsxy, *pu, *pf, *coeffChamp.getDataPointer(), ipar[0], rpar[0], ierr );
                    maxRes = rpar[13]; //   ... returned 
                    l2Res  = rpar[14]; //   ... returned
        //   champResidualOpt( const int & nd,
        //                        const int & nd1a,const int & nd1b,const int & nd2a,const int & nd2b,const int & nd3a,const int & nd3b,const int & nd4a,const int & nd4b,
        //                        const int & gridIndexRange,const int & dimRange,const int & isPeriodic,const int & boundaryCondition,
        //                        const int & mask,const real & xy,const real & rsxy,const real & u,const real & f,const int & ipar,const real & rpar,const int & ierr );
            // int ierr=0;
            // bcOptAd(mg.numberOfDimensions(),
            //        uLocal.getBase(0),uLocal.getBound(0),uLocal.getBase(1),uLocal.getBound(1),
            //        uLocal.getBase(2),uLocal.getBound(2),
            //        indexRangeLocal(0,0), dimLocal(0,0), mg.isPeriodic(0),
            //        *pu, *pmask, *prsxy, *pxy,  bcLocal(0,0),  
            //        *pAddBoundaryForcing, *interfaceType.getDataPointer(), 
            //        bcData.getLength(0), *bcData.getDataPointer(), 
            //        *pdbc, *pbcf[0][0], pbcfOffset[0], ipar[0],rpar[0], ierr );

                if( debug & 1 )
                    printF("+++++ getChampResidual for (grid,side,axis)=(%d,%d,%d) maxRes=%8.2e, l2Res=%8.2e\n",grid,side,axis,maxRes,l2Res);

                maxInterfaceResidual = max(maxInterfaceResidual,maxRes);

        // OV_ABORT("finish me");

            }
        }
    }
                



    return maxInterfaceResidual;
}


// =======================================================================================================================
/// \brief Optionally evaluate the max residual in any interface conditions.
///
/// \details Normally Cgmp will compute the residual in the heat flux jump conditions. The CHAMP conditions, however,
///    do not exactly satisfy the normal jump conditions in Temperaure and heat-flux. In this case each Domain solver
///    is called to compute the residual in the CHAMP conditions.
///
/// \return Return a value of 0 if NO interface residual is available, Return 1 if the residual was computed
// =======================================================================================================================
int Cgad::
getInterfaceResidual( real t, real dt, GridFunction & cgf, Real & residual )
{
    const int & multiDomainProblem = parameters.dbase.get<int>("multiDomainProblem"); 
    const int applyChampInterfaceConditions = parameters.dbase.get<int>("applyChampInterfaceConditions");
    if( multiDomainProblem && applyChampInterfaceConditions )
    {

        const Real maxRes = getChampResidual( cgf.u, parameters, t, dt );

        printP("CGAD::getInterfaceResidual: t=%9.3e: Max residual in champ interface conditions = %8.2e\n",t,maxRes);
        residual = maxRes;

        return 1; // 1 means we have computed a residual
    }




    OV_ABORT("finish me");
    return 0;
}
