// This file automatically generated from addedMassImplicitBoundaryConditions.bC with bpp.
#include "Cgins.h"
#include "Parameters.h"
#include "MappedGridOperators.h"
#include "ParallelUtility.h"
#include "ParallelGridUtility.h"
// #include "Oges.h"
// #include "SparseRep.h"
// #include "App.h"
#include "GridMaterialProperties.h"
#include "Interface.h"
#include "DeformingBodyMotion.h"

#define ForBoundary(side,axis)   for( int axis=0; axis<mg.numberOfDimensions(); axis++ ) for( int side=0; side<=1; side++ )

#define FOR_3D(i1,i2,i3,I1,I2,I3) int I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase();  int I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(); for(i3=I3Base; i3<=I3Bound; i3++) for(i2=I2Base; i2<=I2Bound; i2++) for(i1=I1Base; i1<=I1Bound; i1++)

#define FOR_3(i1,i2,i3,I1,I2,I3) I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase();  I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(); for(i3=I3Base; i3<=I3Bound; i3++) for(i2=I2Base; i2<=I2Bound; i2++) for(i1=I1Base; i1<=I1Bound; i1++)

#define FOR_3IJD(i1,i2,i3,I1,I2,I3,j1,j2,j3,J1,J2,J3) int I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase();  int I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(); int J1Base =J1.getBase(),   J2Base =J2.getBase(),  J3Base =J3.getBase();  for(int i3=I3Base,j3=J3Base; i3<=I3Bound; i3++,j3++) for(int i2=I2Base,j2=J2Base; i2<=I2Bound; i2++,j2++) for(int i1=I1Base,j1=J1Base; i1<=I1Bound; i1++,j1++)


#define FOR_M(m,M)int mBase=M.getBase(), mBound=M.getBound(); for(m=mBase; m<=mBound; m++)

#define U(c)     u(I1,I2,I3,c)   
#define UU(c)   uu(I1,I2,I3,c)
#define UX(c)   ux(I1,I2,I3,c)
#define UY(c)   uy(I1,I2,I3,c)
#define UZ(c)   uz(I1,I2,I3,c)




//    Mixed-derivative BC for component i: 
//          mixedCoeff(i)*u(i) + mixedNormalCoeff(i)*u_n(i) = mixedRHS(i)
#define mixedRHS(component,side,axis,grid)         bcData(component+nc*(0),side,axis,grid)
#define mixedCoeff(component,side,axis,grid)       bcData(component+nc*(1),side,axis,grid)
#define mixedNormalCoeff(component,side,axis,grid) bcData(component+nc*(2),side,axis,grid)


#define AMP(n,m) (delta(n,m) - an[n]*an[m])
// Here is the "P" operator, (1-alpha) n n^T 
#define AMGP(n,m) ((1.-alpha)*an[n]*an[m])
// Here is I-P 
#define AMG(n,m) (delta(n,m) - (1.-alpha)*an[n]*an[m])

// ===============================================================================================
// Macro: evalRightHandSide: evaluate the right-hand-side to the AMP velocity BCs
//
// The velocity BC's can be combined to give two vector equations:
//   (I)  (div(v))*n + (I-n n^T)(tauv.n - zs*v  )/mu = RHS    (Eqn for boundary pt)
//   (II)  v - theta*dt*nu*(I+P) Delta v = RHS                (Eqn for ghost pt)
//  where     P = (alpha-1) n n^T 
// ===============================================================================================


// ======================================================================================
/// \brief Return some parameters for the AMP bulk solid interface conditions
/// \param grid,side,axis,dt (input):
/// \param zs,zp,zf,alpha (output) :
// ======================================================================================
int Cgins::
getBulkSolidAmpParameters( MappedGrid & mg, const int grid, const int side, const int axis, const real dt,
                                                      real & zs, real & zp, real & zf, real & alpha ) 
{
    Parameters & bulkSolidParams = getInterfaceParameters( grid,side,axis,parameters );
    real rhoSolid=bulkSolidParams.dbase.get<real>("rho");
    real lambdaSolid=bulkSolidParams.dbase.get<real>("lambda");
    real muSolid=bulkSolidParams.dbase.get<real>("mu");
    real cp=sqrt((lambdaSolid+2.*muSolid)/rhoSolid);
    real cs=sqrt(muSolid/rhoSolid);
                    
    zp=rhoSolid*cp;
    zs=rhoSolid*cs;
  // fluid impedance = rho*H/dt 
    const real & fluidDensity = parameters.dbase.get<real >("fluidDensity");
    assert( dt>0. && fluidDensity>0. );

    const real nu = parameters.dbase.get<real>("nu");
    const real mu = parameters.dbase.get<real >("mu");
    assert( fabs(fluidDensity*nu/mu-1.) < REAL_EPSILON*100. );

    const real fluidAddedMassLengthScale =  parameters.dbase.get<real>("fluidAddedMassLengthScale");
    const real zfOld=fluidDensity*fluidAddedMassLengthScale/dt; 

    const real zfNew = 2.*(fluidDensity*mu)/(zp*dt);

  // Need dn -- grid spacing in the normal direction. 

    bool useFluidImpedanceMuByH=false;  // ++++++++++++++++++++++++++   TEST ++++++++++++++

    const real & zfMuByH = parameters.dbase.get<real>("zfMuByH");  
    const real & zfMono = parameters.dbase.get<real>("zfMono");  
    const real & zfMuRhoByZpDt = parameters.dbase.get<real>("zfMuRhoByZpDt");  

    real zfmu=1.,dn=1.;

    bool useMonolithicImpedance=zfMono>0.; // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    real zfm=1.;
    if( true || useFluidImpedanceMuByH || useMonolithicImpedance )
    {
    //  --------------- Compute the normal distance of the first grid line  ---------------------
    //   *** SAVE THIS VALUE TO AVOID RECOMPUTING ****
        parameters.getNormalGridSpacing( mg,side,axis, dn );
        zfmu = zfMuByH*mu/dn; 


    // --- Impedance from the monolithic algorithm ---
        real rho = fluidDensity;

        const int axisp1 = (axis + 1) % mg.numberOfDimensions();
        const real ds = mg.gridSpacing(axisp1);  

        real k = 1./dn;
    // real k = 1./(2.*ds);

        real beta=sqrt( rho/(mu*dt) + k*k );
    // real zfm1;
        zfm = (rho*rho + ( rho*beta*zp + 4.*rho*mu*k*k )*dt - ( 4.*(beta-k)*mu*mu*k*k*k )*dt*dt )/( (rho+ zp*dt*(beta-k) )*dt*k );

<<<<<<< HEAD
        printF("GGGGGG getBulkSolidAmpParameters: ");
        printF("old zfm=%.2e, ",zfm);

    // this is asymptotically equivalent to old zfm
        zfm = 2.*mu*k+1./(k*dt);
        printF("new zfm=%.2e\n",zfm);

        printF("GGGGGG getBulkSolidAmpParameters: rhoSolid=%.2e, zfNew=(%.2g)*(rho*mu)/(zp*dt) = %.2e,  zfmu=(%.2e)*mu/h=%.2e, zfMonolithic=%.2e (zfMono=%.2e)\n",
                      rhoSolid,zfMuRhoByZpDt,zfNew, zfMuByH, zfmu,  zfm,zfMono);
=======
        if( false )
            printF("GGGGGG getBulkSolidAmpParameters: rhoSolid=%.2e, zfNew=(%.2g)*(rho*mu)/(zp*dt) = %.2e,  zfmu=(%.2e)*mu/h=%.2e, zfMonolithic=%.2e (zfMono=%.2e)\n",
                          rhoSolid,zfMuRhoByZpDt,zfNew, zfMuByH, zfmu,  zfm,zfMono);
>>>>>>> wdh: changes to cgmx for 3D GDM interfaces
    
    }
    

    bool useNew=true;  // ************
    if( useMonolithicImpedance )
    {
        zf = zfMono*zfm;

    // 8/28/18: some cases may be fixed by uncommenting the following line
    // zf = zfNew;
    }
    else if( useFluidImpedanceMuByH )
    {
        zf=zfmu;
    }
    else
    {
      
        if( useNew )
        {
      // -- new way ---
            zf=zfNew;
        }
        else
        {
            zf=zfOld;
        }
    }
    
    alpha = zf/(zf+zp);

    if( debug() & 2 )
    {
        FILE *&debugFile = parameters.dbase.get<FILE* >("debugFile");
        fPrintF(debugFile,"GGGGGG getBulkSolidAmpParameters: useFluidImpedanceMuByH=%i, useNew=%i, zfOld=rho*L/dt = %.2e  zfNew=2*(rho*mu)/(zp*dt) = %.2e, zfmu=(%.2e)*mu/h=%.2e, alpha=%.2e "
                        " rho=%.2e mu=%.2e zp=%.2e dt=%.2e, zf=%.2e\n",(int)useFluidImpedanceMuByH,(int)useNew,zfOld,zfNew,
                        zfMuByH,zfmu,alpha,fluidDensity,mu,zp,dt,zf);

    }
    return 0;
}



int Cgins::
addedMassImplicitBoundaryConditions(int option, 
                                                                        realMappedGridFunction & u, 
                                                                        realMappedGridFunction &uL,
                                                                        realMappedGridFunction &uOld,  
                                                                        realMappedGridFunction & gridVelocity,
                                                                        real t,
                                                                        int scalarSystem,
                                                                        int grid, 
                                                                        real dt0 )
// ======================================================================================
/// \brief Assign the right-hand side (or eval the residuals) to the Added-Mass Implicit Boundary conditions
///       for INS + BULK ELASTIC SOLID
///
/// \param option (input/output) : 
///          option = 0 : fill in the right-hand-side (before the implicit solve)
///          option = 1 : evaluate the residual in the added-mass equations (after the implicit solve)
/// \param u (input/output) : apply boundary conditions to this grid function.
/// \param gridVelocity (input) : for BC's on moving grids.
/// \param t (input) : time (new time)
/// \param scalarSystem (input) : 
/// \param grid (input) : component grid number.
///
/// \Return : 1= BC's were applied, 0=BC's were not applied
// ==========================================================================================
{
    int applied=0; // set to 1 if AMP BC's are applied

    const bool & useAddedMassAlgorithm = parameters.dbase.get<bool>("useAddedMassAlgorithm");
    const bool & projectAddedMassVelocity = parameters.dbase.get<bool>("projectAddedMassVelocity");
    const int initialConditionsAreBeingProjected = parameters.dbase.get<int>("initialConditionsAreBeingProjected");
    const bool & useImplicitAmpBCs = parameters.dbase.get<bool>("useImplicitAmpBCs");
    const bool & predictedBoundaryPressureNeeded = parameters.dbase.get<bool>("predictedBoundaryPressureNeeded");

    bool applyAddedMass = ( useAddedMassAlgorithm && 
                                                    projectAddedMassVelocity && 
                                                    parameters.gridIsMoving(grid)
                                                    && !initialConditionsAreBeingProjected 
                                                    && t!=0. );

  // For testing with Cgins we may fill in the implicit AMP BCs even when running Cgins alone
    applyAddedMass = applyAddedMass || (useAddedMassAlgorithm && useImplicitAmpBCs);
    

    if( !applyAddedMass )
    {
        return applied;
    } 
    

    MovingGrids & movingGrids = parameters.dbase.get<MovingGrids >("movingGrids");
    if( movingGrids.getNumberOfDeformingBodies()==0 && !useImplicitAmpBCs )
        return applied;

    real dt=dt0;
    if( dt<= 0. )
        dt = parameters.dbase.get<real>("dt");  // *wdh* 2017/05/31
  // const real & dt = parameters.dbase.get<real>("dt");
    const real epsT = REAL_EPSILON*10.;
    assert( dt> epsT );
    
    FILE *&debugFile = parameters.dbase.get<FILE* >("debugFile");
    FILE *&pDebugFile = parameters.dbase.get<FILE* >("pDebugFile");

    if( t <= 2.*dt && debug() &4 )
    {
        fPrintF(debugFile,"\n"
                      "--------------------------------------------------------------------------------------------------\n"
                      " --INS-- addedMassImplicitBoundaryConditions: ADDED MASS ALGORITHM t=%8.2e\n"
                      "--------------------------------------------------------------------------------------------------\n"
                        ,t);
    }
    
    const real nu = parameters.dbase.get<real>("nu");
    const real rho = parameters.dbase.get<real>("rho");
    const real & fluidDensity = parameters.dbase.get<real >("fluidDensity");
    assert( rho==fluidDensity );
    
    const real mu = rho*nu;
    const real implicitFactor = parameters.dbase.get<real >("implicitFactor");
    const real theta = implicitFactor;

    const bool & projectNormalComponentOfAddedMassVelocity =
                              parameters.dbase.get<bool>("projectNormalComponentOfAddedMassVelocity");
    const bool & projectVelocityOnBeamEnds = parameters.dbase.get<bool>("projectVelocityOnBeamEnds"); 

    const int addedMassVelocityBC = parameters.dbase.get<int>("addedMassVelocityBC");
    const int correctionStage = parameters.dbase.get<int>("correctionStage");
    real pcSwitch = correctionStage>0 ? 1. : 0.;

  // pcSwitch=0.;
    

    const bool twilightZoneFlow = parameters.dbase.get<bool >("twilightZoneFlow");

    assert(  useAddedMassAlgorithm && projectAddedMassVelocity && (useImplicitAmpBCs || parameters.gridIsMoving(grid)) );
    
    MappedGrid & mg = *u.getMappedGrid();
    const IntegerArray & gid = mg.gridIndexRange();
    const int numberOfDimensions = mg.numberOfDimensions();
    Range Rx=numberOfDimensions;

    const bool gridIsMoving = parameters.gridIsMoving(grid);

  // get components of solution
    const int uc = parameters.dbase.get<int >("uc");
    const int vc = parameters.dbase.get<int >("vc");
    const int wc = parameters.dbase.get<int >("wc");
    const int tc = parameters.dbase.get<int >("tc");
    const int pc = parameters.dbase.get<int >("pc");
    const int & nc = parameters.dbase.get<int >("nc");

    const int orderOfAccuracy=min(4,parameters.dbase.get<int >("orderOfAccuracy"));
    Range V = Range(uc,uc+numberOfDimensions-1);

    BoundaryData::BoundaryDataArray & pBoundaryData = parameters.getBoundaryData(grid); // this will create the BDA if it is not there
    std::vector<BoundaryData> & boundaryDataArray =parameters.dbase.get<std::vector<BoundaryData> >("boundaryData");
    BoundaryData & bd = boundaryDataArray[grid];
            
    const Parameters::InterfaceCommunicationModeEnum & interfaceCommunicationMode= 
        parameters.dbase.get<Parameters::InterfaceCommunicationModeEnum>("interfaceCommunicationMode");

    const Parameters::KnownSolutionsEnum & knownSolution = 
        parameters.dbase.get<Parameters::KnownSolutionsEnum >("knownSolution"); 

  // if( FALSE && !parameters.gridIsMoving(grid) )
  // {
  //   // ---- this is a test run using non-moving grids ----
  //   ForBoundary(side,axis)
  //   {
  //     if( mg.boundaryCondition(side,axis)==Parameters::noSlipWall )
  //     {
  //       printF(" XXXX addedMassImplicitBoundaryConditions: TEST RUN using non-moving grids XXXX \n\n");

  //       applied=1;  // we have applied an AMP BC

  //       Index Ib1,Ib2,Ib3, Ig1,Ig2,Ig3, Ip1,Ip2,Ip3;
  //       getBoundaryIndex(mg.gridIndexRange(),side,axis,Ib1,Ib2,Ib3);
  //       getGhostIndex(mg.gridIndexRange(),side,axis,Ig1,Ig2,Ig3,+1);  // first ghost line 

  //       RealArray solidTraction(Ib1,Ib2,Ib3,numberOfDimensions), vSolidLocal(Ib1,Ib2,Ib3,numberOfDimensions);

  //       // Do this for now
  //       real zf=1., zp=1., zs=1., alpha=.5; // defaults when tetsing cgins alone
  //       solidTraction=0.;
  //       vSolidLocal=0.;

  //       OV_GET_SERIAL_ARRAY(real,u   ,uLocal);
  //       OV_GET_SERIAL_ARRAY(real,uOld,uOldLocal);
                
  //       if( projectNormalComponentOfAddedMassVelocity )
  //         mg.update(MappedGrid::THEvertexBoundaryNormal);

  //       OV_GET_VERTEX_BOUNDARY_NORMAL(mg,side,axis,normal);         

  //       // evaluate the right-hand-side to the AMP velocity BCs
  //       evalAmpRightHandSide();

                
  //     }
            
  //   }
        
  //   // OV_ABORT("addedMassImplicitBoundaryConditions: STOP HERE FOR NOW");
        
  //   return applied;
  // }
    

  // ========== REAL DEFORMING GRID CASE ==========

   // --- DEFORMING BULK SOLID ----
   //   For a deforming bulk solid we only apply moving no-slip wall Bc's to the face that is actually moving.
   // 
   //           +------------------------+
   //           |                        |
   //           |  deforming grid        |
   //  noslip   |                        |  noSlip
   //   u=0     |                        |   u=0
   //           |                        |
   //           |    moving interface    |
   //           +------------------------+
   //                u=gridVelocity

  // -- extract parameters from any deforming solids ---
    
    if( bd.dbase.has_key("deformingBodyNumber") )
    {
        const real & fluidDensity = parameters.dbase.get<real >("fluidDensity");
        assert( fluidDensity>0. );

        const real fluidAddedMassLengthScale =  parameters.dbase.get<real>("fluidAddedMassLengthScale");

        int (&deformingBodyNumber)[2][3] = bd.dbase.get<int[2][3]>("deformingBodyNumber");
        Index Ib1,Ib2,Ib3, Ig1,Ig2,Ig3, Ip1,Ip2,Ip3;
        for( int side=0; side<=1; side++ )
        {
            for( int axis=0; axis<numberOfDimensions; axis++ )
            {

      	if( deformingBodyNumber[side][axis]>=0 )
      	{
                    if( false )
                        printF("WWWWW implicit: apply MOVING noSlipWall BC to a DEFORMING BULK SOLID, grid=%i (side,axis)=(%i,%i), t=%.2e\n",
                                      grid,side,axis,t);

        	  int body=deformingBodyNumber[side][axis];
        	  if( t<=0. )
          	    printF("--INS-- grid=%i, (side,axis)=(%i,%i) belongs to deforming body %i\n",grid,side,axis,body);

        	  DeformingBodyMotion & deform = movingGrids.getDeformingBody(body);

        	  getBoundaryIndex(mg.gridIndexRange(),side,axis,Ib1,Ib2,Ib3);
                    getGhostIndex(mg.gridIndexRange(),side,axis,Ig1,Ig2,Ig3,+1);  // first ghost line 
                    getGhostIndex(mg.gridIndexRange(),side,axis,Ip1,Ip2,Ip3,-1);  // first line in

        	  realArray vSolid(Ib1,Ib2,Ib3,Rx); // holds velocity of solid on the boundary
                    #ifndef USE_PPP
                    deform.getVelocityBC( t, side,axis,grid, mg, Ib1,Ib2,Ib3, vSolid );
                    #else
                        OV_ABORT("finish me");
                    #endif

          // OV_GET_SERIAL_ARRAY(real,gridVelocity,gridVelocityLocal);
        	  OV_GET_SERIAL_ARRAY(real,u   ,uLocal);
                    OV_GET_SERIAL_ARRAY(real,uOld,uOldLocal);
        	  OV_GET_SERIAL_ARRAY(real,vSolid,vSolidLocal);

                    if( projectNormalComponentOfAddedMassVelocity )
                        mg.update(MappedGrid::THEvertexBoundaryNormal);

                    OV_GET_VERTEX_BOUNDARY_NORMAL(mg,side,axis,normal);
        	  
                    if( !deform.isBulkSolidModel() && !deform.isBeamModel() )
        	  {

          	    OV_ABORT("addedMassImplicitBoundaryConditions::ERROR: un-expected deformation type");

        	  }
        	  else if( deform.isBulkSolidModel() )
        	  {

            // **********************************************************************
            // *************** ADDED MASS BC FOR A BULK SOLID ***********************
            // **********************************************************************

                        applied=1;  // we have applied an AMP BC

                        RealArray solidTraction;
                        if( interfaceCommunicationMode==Parameters::requestInterfaceDataWhenNeeded )
                        {
                            if( debug() & 4 )
                                fPrintF(debugFile,"--INS-- PIV: REQUEST interface traction at t=%9.3e\n",t);

                            InterfaceData interfaceData;
                            Range Rx=numberOfDimensions;
                            interfaceData.u.redim(Ib1,Ib2,Ib3,Rx); // traction is returned here 
                            interfaceData.t=t;
                            interfaceData.u=0;

                            int interfaceDataOptions = Parameters::tractionInterfaceData;
                            bool saveTimeHistory=true;
                            getInterfaceData( t, grid, side, axis, 
                                                                interfaceDataOptions,
                                                                interfaceData.u,
                                                                parameters,saveTimeHistory );

                            solidTraction=interfaceData.u;
                            if( t <= 2.*dt && debug() & 4 ) {
                                ::display(solidTraction(Ib1,Ib2,Ib3,1),"--INS-- AM-IMP-BC: Here is the SOLID TRACTION (I1,I2,I3,1)",debugFile,"%6.3f ");
                            }
                        }
                        else
                        {
                            OV_ABORT("finish me");
                        }

            // *new* April 1, 2018
                        real zs,zp,zf,alpha;
                        getBulkSolidAmpParameters( mg,grid,side,axis,dt, zs,zp,zf,alpha );

	    // real zpOld;
            // // old way:
            // deform.getBulkSolidParameters( zpOld );

            // // new way:
            // // Retrieve the parameters from the bulk solid
            // // FIX ME -- lookup first time and then save locally 
            // Parameters & bulkSolidParams = getInterfaceParameters( grid,side,axis,parameters );
            // real rhoSolid=bulkSolidParams.dbase.get<real>("rho");
            // real lambdaSolid=bulkSolidParams.dbase.get<real>("lambda");
            // real muSolid=bulkSolidParams.dbase.get<real>("mu");
            // real cp=sqrt((lambdaSolid+2.*muSolid)/rhoSolid);
            // real cs=sqrt(muSolid/rhoSolid);
                    

            // real zp=rhoSolid*cp;
            // real zs=rhoSolid*cs;
                        
            // if( t<=3.*dt && debug() & 4  )
            // {
            //   fPrintF(debugFile,"--INS-- AMP-IMP-BC: rhoSolid=%9.3e cp=%9.3e cs=%9.3e zp=%9.3e (old: zp=%9.3e) zs=%9.3e\n",
            //          rhoSolid,cp,cs,zp,zpOld,zs);
            //   // printF("  fluidAddedMassLengthScale=%9.3e\n",fluidAddedMassLengthScale);
            // }
                        
	    // // const real & fluidDensity = parameters.dbase.get<real >("fluidDensity");
          	    
            // // fluid impedance = rho*H/dt 
	    // assert( dt>0. );
                        
            // const real zfOld=fluidDensity*fluidAddedMassLengthScale/dt; 
            // // *new* June 30, 2018: 
            // real zf= 2.*(fluidDensity*mu)/(zp*dt);
            // const real alpha = zf/(zf+zp);
            // printF("BBBBBB addedMassImpBC: zf=rho*L/dt = %.2e  zf=2*(rho*mu)/(zp*dt) = %.2e, alpha=%.2e "
            //        " rho=%.2e mu=%.2e zp=%.2e dt=%.2e\n",zfOld,zf,alpha,rho,mu,zp,dt);
            // // zf=zfOld;



                        bool useExactSolidVelocity=false;
                        bool useExactSolidTraction=false;
                        
                        if( useExactSolidVelocity && knownSolution==Parameters::userDefinedKnownSolution )
                        {
                            printF("\n--INS-- AMP-IMP-BC: **TEST** set exact KNOWN-SOLUTION values for solid velocity ***TEMP***, t=%9.3e\n\n",t);
                            int body=0;
                            parameters.getUserDefinedDeformingBodyKnownSolution( body,Parameters::boundaryVelocity,
                                                                                                                                      t, grid, mg, Ib1,Ib2,Ib3,Rx,vSolid );
                        }
                        if( useExactSolidTraction && knownSolution==Parameters::userDefinedKnownSolution )
                        {
                            printF("\n--INS-- AMP-IMP-BC: **TEST** set exact KNOWN-SOLUTION values for solid traction ***TEMP***, t=%9.3e\n\n",t);
                            int body=0;
                            parameters.getUserDefinedDeformingBodyKnownSolution( body,Parameters::boundaryTraction,
                                                                                                                                      t, grid, mg, Ib1,Ib2,Ib3,Rx,solidTraction );
                        }


            // evaluate the right-hand-side to the AMP velocity BCs
                        {
                            if( debug() & 4 || t<=3.*dt ) 
                                fPrintF(debugFile,"--INS-- AMP-BC-IMP FOR BULK SOLID MODEL, t=%9.3e dt=%9.3e zf=%9.3e zp=%9.3e zs=%9.3e "
                                              "alpha=%9.2e **new**\n",t,dt,zf,zp,zs,alpha);
                            if( predictedBoundaryPressureNeeded==0 )
                            {
                                printF("--INS-- AMP-BC-IMP FOR BULK SOLID MODEL: ERROR: predicted pressure near boundary needed but "
                                              "predictedBoundaryPressureNeeded=0 !\n");
                                OV_ABORT("error");
                            }
                            OV_GET_SERIAL_ARRAY(real,gridVelocity,gridVelocityLocal);
                            if( debug() & 4 )
                            {
                                real maxDiff = max(fabs(gridVelocityLocal(Ib1,Ib2,Ib3,Rx)-vSolidLocal(Ib1,Ib2,Ib3,Rx)));
                                printF("--INS--AMP-IMP-BC: t=%.3e, max diff [gridVelocity - vSolidLocal] =%.2e\n",t,maxDiff);
                            }
                            const real beta = (mu*theta*dt/rho);
              // -- evaluate derivatives that appear in the INS equations ---       
                            MappedGridOperators & opOld = *(uOld.getOperators()); 
                            realSerialArray uOldxx(Ib1,Ib2,Ib3,V), uOldyy(Ib1,Ib2,Ib3,V);
                            realSerialArray pOldx(Ib1,Ib2,Ib3), pOldy(Ib1,Ib2,Ib3), pOldz;
                            opOld.derivative(MappedGridOperators::xxDerivative ,uOldLocal,uOldxx,Ib1,Ib2,Ib3,V);
                            opOld.derivative(MappedGridOperators::yyDerivative ,uOldLocal,uOldyy,Ib1,Ib2,Ib3,V);
                            opOld.derivative(MappedGridOperators::xDerivative ,uOldLocal,pOldx,Ib1,Ib2,Ib3,pc);
                            opOld.derivative(MappedGridOperators::yDerivative ,uOldLocal,pOldy,Ib1,Ib2,Ib3,pc);
                            if( numberOfDimensions==3 )
                            {
                                pOldz.redim(Ib1,Ib2,Ib3);
                                opOld.derivative(MappedGridOperators::zDerivative ,uOldLocal,pOldz,Ib1,Ib2,Ib3,pc);
                            }
              // -- We assume a predicted value for "p" at the new time
                            MappedGridOperators & op = *(u.getOperators()); 
                            realSerialArray px(Ib1,Ib2,Ib3), py(Ib1,Ib2,Ib3),pz;
                            op.derivative(MappedGridOperators::xDerivative ,uLocal,px,Ib1,Ib2,Ib3,pc);
                            op.derivative(MappedGridOperators::yDerivative ,uLocal,py,Ib1,Ib2,Ib3,pc);
                            if( numberOfDimensions==3 )
                            {
                                pz.redim(Ib1,Ib2,Ib3);
                                op.derivative(MappedGridOperators::zDerivative ,uLocal,pz,Ib1,Ib2,Ib3,pc);
                            }
                            RealArray nSigmaFluidN, nSigmaSolidN;
                            if( pcSwitch>0. )
                            {
                // --- evaluate the normal component of the fluid traction during correction stages ----
                                nSigmaSolidN.redim(Ib1,Ib2,Ib3);
                                nSigmaFluidN.redim(Ib1,Ib2,Ib3);
                                realSerialArray ux(Ib1,Ib2,Ib3,V), uy(Ib1,Ib2,Ib3,V);
                                op.derivative(MappedGridOperators::xDerivative ,uLocal,ux,Ib1,Ib2,Ib3,V);
                                op.derivative(MappedGridOperators::yDerivative ,uLocal,uy,Ib1,Ib2,Ib3,V);
                // Remove divergence from Tau
                // n.Sigma.n = - p + n.Tau.n 
                // nTaun = mu*( an1*( 2.*ux*an1 + (uy+vx)*an2) + an2*( (uy+vx)*an1 +2.*vy*an2 ) )
                                nSigmaFluidN = -uLocal(Ib1,Ib2,Ib3,pc) 
                                            + mu*( normal(Ib1,Ib2,Ib3,0)*( 
                                                                    +2.*ux(Ib1,Ib2,Ib3,uc)*normal(Ib1,Ib2,Ib3,0) + 
                                                                          (uy(Ib1,Ib2,Ib3,uc)+ux(Ib1,Ib2,Ib3,vc))*normal(Ib1,Ib2,Ib3,1) ) + 
                                                          normal(Ib1,Ib2,Ib3,1)*( 
                                                                          (uy(Ib1,Ib2,Ib3,uc)+ux(Ib1,Ib2,Ib3,vc))*normal(Ib1,Ib2,Ib3,0) 
                                                                    +2.*uy(Ib1,Ib2,Ib3,vc)*normal(Ib1,Ib2,Ib3,1) ) );
                // **FLIP SIGN -- should FLIP sign on solidTraction to start with *fix me*
                                nSigmaSolidN = -(normal(Ib1,Ib2,Ib3,0)*solidTraction(Ib1,Ib2,Ib3,0) + normal(Ib1,Ib2,Ib3,1)*solidTraction(Ib1,Ib2,Ib3,1));
                                assert( numberOfDimensions==2 );
                            }
                            RealArray delta(3,3); // Dirac delta
                            delta=0.; delta(0,0)=1.; delta(1,1)=1.; delta(2,2)=1.;
                            RealArray solidChar(Ib1,Ib2,Ib3,Rx);   // solid characteristic variable 
                            RealArray fluidRhs(Ib1,Ib2,Ib3,Rx);    // fluid interior update without implicit viscous term
              // Solid characteristic variable used in the RHS 
              // Use MINUS of solid traction since the solid normal is in the opposite direction to the fluid
                            solidChar(Ib1,Ib2,Ib3,Rx) = -solidTraction(Ib1,Ib2,Ib3,Rx) + zs*vSolidLocal(Ib1,Ib2,Ib3,Rx);
                            if( debug() & 4 )
                            {
                                ::display(solidChar,"--INS-- IMP-AMP-V-BC: solid-characteristic: sigmas.nv + zs*vs",debugFile,"%.3e ");
                            }
              // fluidRhs(Ib1,Ib2,Ib3,dir) = fluid interior update without implicit viscous term 
                            for( int dir=0; dir<numberOfDimensions; dir++ )
                            {
                // pd = px, py or pz 
                                RealArray & pd    = dir==0 ? px    : dir==1? py    : pz;
                                RealArray & pOldd = dir==0 ? pOldx : dir==1? pOldy : pOldz;
                // fluidRhs = v^{n-1} - (theta*dt/rho)*grad(p)^n + ((1-theta)*dt/rho)*( -grad(p)^{n-1} + mu*Delta(v)^{n-1} )
                                fluidRhs(Ib1,Ib2,Ib3,dir)= uOldLocal(Ib1,Ib2,Ib3,uc+dir) 
                                    + (-implicitFactor*dt/rho)*( pd(Ib1,Ib2,Ib3) )
                                    + ((1.-implicitFactor)*dt/rho)*( 
                                                mu*( uOldxx(Ib1,Ib2,Ib3,uc+dir)+uOldyy(Ib1,Ib2,Ib3,uc+dir) ) - pOldd(Ib1,Ib2,Ib3) 
                                                                                                  );
                            }
                            if( twilightZoneFlow )
                            {
                // ------- ADD TWILIGHT ZONE --------
                                OGFunction & e = *(parameters.dbase.get<OGFunction* >("exactSolution"));
                                const bool rectangular=false;
                                OV_GET_SERIAL_ARRAY(real,mg.vertex(),xLocal);
                                if( debug() & 4 )
                                {
                                    RealArray pex(Ib1,Ib2,Ib3), pey(Ib1,Ib2,Ib3);
                                    e.gd( pex,xLocal,mg.numberOfDimensions(),rectangular,0,1,0,0,Ib1,Ib2,Ib3,pc,t); // t=new time 
                                    e.gd( pey,xLocal,mg.numberOfDimensions(),rectangular,0,0,1,0,Ib1,Ib2,Ib3,pc,t); // t=new time 
                                    real maxErrPx = max(fabs(px-pex));
                                    real maxErrPy = max(fabs(py-pey));
                                    fprintf(debugFile,">>>INS IMP-AMP-V-BC: error in predicted (px,py)=(%9.3e,%9.3e) at t=%9.3e\n",maxErrPx,maxErrPy,t);
                                    ::display(fabs(px-pex),"error in px",debugFile,"%8.2e ");
                                    ::display(fabs(py-pey),"error in py",debugFile,"%8.2e ");
                                }
                                realSerialArray ue(Ib1,Ib2,Ib3,V), uex(Ib1,Ib2,Ib3,V), uey(Ib1,Ib2,Ib3,V);
                                RealArray tractione(Ib1,Ib2,Ib3,Rx);
                                e.gd( ue  ,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,0,Ib1,Ib2,Ib3,V,t); // t=new time 
                                e.gd( uex ,xLocal,mg.numberOfDimensions(),rectangular,0,1,0,0,Ib1,Ib2,Ib3,V,t); 
                                e.gd( uey ,xLocal,mg.numberOfDimensions(),rectangular,0,0,1,0,Ib1,Ib2,Ib3,V,t); 
                // tau_ij = mu*( partial_j u_i + partial_i u_j )
                                if( numberOfDimensions==2 )
                                {
                  // ====== Get traction in two dimensions =====
                                    RealArray taue(Ib1,Ib2,Ib3,3); // hold viscous stress tensor tau
                                    taue(Ib1,Ib2,Ib3,0) = (2.*mu)*( uex(Ib1,Ib2,Ib3,uc) );                        // tau_11
                                    taue(Ib1,Ib2,Ib3,1) = (   mu)*( uey(Ib1,Ib2,Ib3,uc) + uex(Ib1,Ib2,Ib3,vc) );  // tau_12 = tau_21
                                    taue(Ib1,Ib2,Ib3,2) = (2.*mu)*( uey(Ib1,Ib2,Ib3,vc) );                        // tau_22
                  //  t1 = tau11*n1 + tau12*n2 
                                    tractione(Ib1,Ib2,Ib3,0) = taue(Ib1,Ib2,Ib3,0)*normal(Ib1,Ib2,Ib3,0) + taue(Ib1,Ib2,Ib3,1)*normal(Ib1,Ib2,Ib3,1);
                  //  t2 = tau21*n1 + tau22*n2 
                                    tractione(Ib1,Ib2,Ib3,1) = taue(Ib1,Ib2,Ib3,1)*normal(Ib1,Ib2,Ib3,0) + taue(Ib1,Ib2,Ib3,2)*normal(Ib1,Ib2,Ib3,1);
                                }
                                else 
                                {
                  // ====== Get traction in three dimensions =====
                                    RealArray uez(Ib1,Ib2,Ib3,V);
                                    e.gd( uez ,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,1,Ib1,Ib2,Ib3,V,t); 
                                    RealArray taue(Ib1,Ib2,Ib3,6);
                                    taue(Ib1,Ib2,Ib3,0) = (2.*mu)*( uex(Ib1,Ib2,Ib3,uc) );                        // tau_11
                                    taue(Ib1,Ib2,Ib3,1) = (   mu)*( uey(Ib1,Ib2,Ib3,uc) + uex(Ib1,Ib2,Ib3,vc) );  // tau_12 = tau_21
                                    taue(Ib1,Ib2,Ib3,2) = (   mu)*( uez(Ib1,Ib2,Ib3,uc) + uex(Ib1,Ib2,Ib3,wc) );  // tau_13 = tau_31
                                    taue(Ib1,Ib2,Ib3,3) = (2.*mu)*( uey(Ib1,Ib2,Ib3,vc) );                        // tau_22
                                    taue(Ib1,Ib2,Ib3,4) = (   mu)*( uez(Ib1,Ib2,Ib3,vc) + uey(Ib1,Ib2,Ib3,wc) );  // tau_23 = tau_32
                                    taue(Ib1,Ib2,Ib3,5) = (2.*mu)*( uez(Ib1,Ib2,Ib3,wc) );                        // tau_33
                  //  t1 = tau11*n1 + tau12*n2 + tau13*n3 
                                    tractione(Ib1,Ib2,Ib3,0) = ( taue(Ib1,Ib2,Ib3,0)*normal(Ib1,Ib2,Ib3,0) + 
                                                                                              taue(Ib1,Ib2,Ib3,1)*normal(Ib1,Ib2,Ib3,1) +
                                                                                              taue(Ib1,Ib2,Ib3,2)*normal(Ib1,Ib2,Ib3,2) );
                  //  t2 = tau21*n1 + tau22*n2 + tau23*n3 
                                    tractione(Ib1,Ib2,Ib3,1) = ( taue(Ib1,Ib2,Ib3,1)*normal(Ib1,Ib2,Ib3,0) + 
                                                                                              taue(Ib1,Ib2,Ib3,3)*normal(Ib1,Ib2,Ib3,1) +
                                                                                              taue(Ib1,Ib2,Ib3,4)*normal(Ib1,Ib2,Ib3,2) );
                  //  t3 = tau31*n1 + tau32*n2 + tau33*n3 
                                    tractione(Ib1,Ib2,Ib3,2) = ( taue(Ib1,Ib2,Ib3,2)*normal(Ib1,Ib2,Ib3,0) + 
                                                                                              taue(Ib1,Ib2,Ib3,4)*normal(Ib1,Ib2,Ib3,1) +
                                                                                              taue(Ib1,Ib2,Ib3,5)*normal(Ib1,Ib2,Ib3,2) );
                                }
                                solidChar(Ib1,Ib2,Ib3,Rx) = tractione(Ib1,Ib2,Ib3,Rx) + zs*ue(Ib1,Ib2,Ib3,V);
                // fluidRhs = v^{n-1} + (theta*dt/rho)*grad(p)^n + ((1-theta)*dt/rho)*( -grad(p)^{n-1} + mu*Delta(v)^{n-1} )
                //          = [ I - (mu*theta*dt/rho)*Delta ]* v^n 
                                RealArray & uexx = uex; RealArray & ueyy =uey;  // Reuse arrays
                                e.gd( uexx ,xLocal,mg.numberOfDimensions(),rectangular,0,2,0,0,Ib1,Ib2,Ib3,V,t); 
                                e.gd( ueyy ,xLocal,mg.numberOfDimensions(),rectangular,0,0,2,0,Ib1,Ib2,Ib3,V,t); 
                                vSolidLocal(Ib1,Ib2,Ib3,Rx) = ue(Ib1,Ib2,Ib3,V);
                                if( numberOfDimensions==2 )
                                {
                  //  [ I - (mu*theta*dt/rho)*Delta ]* v^n
                                    printF("addedMassImpBC: beta=%9.3e, max(uexx=%8.2e)\n",beta,max(fabs(uexx)));
                                    fluidRhs(Ib1,Ib2,Ib3,Rx)=ue(Ib1,Ib2,Ib3,V) - beta*( uexx+ueyy ); 
                                }
                                else
                                {
                                    RealArray & uezz = tractione; // Reuse arrays
                                    e.gd( uezz ,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,2,Ib1,Ib2,Ib3,V,t); 
                  //  [ I - (mu*theta*dt/rho)*Delta ]* v^n
                  // *** THIS IS WRONG ***
                                    fluidRhs(Ib1,Ib2,Ib3,Rx)=ue(Ib1,Ib2,Ib3,V) - (mu*theta*dt/rho)*( uexx+ueyy+uezz ); 
                                }
                                if( pcSwitch>0. )
                                {
                                    RealArray pe(Ib1,Ib2,Ib3);
                                    e.gd( pe,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,0,Ib1,Ib2,Ib3,pc,t); // t=new time 
                                    nSigmaFluidN -= -pe + normal(Ib1,Ib2,Ib3,0)*tractione(Ib1,Ib2,Ib3,0) + normal(Ib1,Ib2,Ib3,1)*tractione(Ib1,Ib2,Ib3,1);
                                    assert( numberOfDimensions==2 );
                                }
                            }
            //  -- "vector" versions of macros that use full normal vector : 
                        #define AMPV(n,m) (delta(n,m) - normal(Ib1,Ib2,Ib3,n)*normal(Ib1,Ib2,Ib3,m) )
            // Here is the "P" operator, (1-alpha) n n^T 
                        #define AMGPV(n,m) ((1.-alpha)*normal(Ib1,Ib2,Ib3,n)*normal(Ib1,Ib2,Ib3,m))
            // Here is I-P 
                        #define AMGV(n,m) (delta(n,m) - (1.-alpha)*normal(Ib1,Ib2,Ib3,n)*normal(Ib1,Ib2,Ib3,m))
                            for( int n=0; n<numberOfDimensions; n++ )
                            {
                // --- RHS FOR GHOST EQUATION: ---
                //  RHS = (1/mu)*( I - n n^T)( solidTraction + zs*vs )
                                uLocal(Ig1,Ig2,Ig3,uc+n) = (1./mu)*( AMPV(n,0)*solidChar(Ib1,Ib2,Ib3,0) + AMPV(n,1)*solidChar(Ib1,Ib2,Ib3,1) );
                // --- RHS FOR BOUNDARY EQUATION: ---
                // printF("AM-IMP_BC: n=%i,  AMGV=%9.3e %9.3e %9.3e,  AMGPV=%9.3e %9.3e\n",n,
                // RHS = (I-P)*( v - beta*Delta(v) ) + P*vs 
                                uLocal(Ib1,Ib2,Ib3,uc+n) = ( AMGV(n,0) *fluidRhs(Ib1,Ib2,Ib3,0)    + AMGV(n,1) *fluidRhs(Ib1,Ib2,Ib3,1) + 
                                                                                          AMGPV(n,0)*vSolidLocal(Ib1,Ib2,Ib3,0) + AMGPV(n,1)*vSolidLocal(Ib1,Ib2,Ib3,1) );
                                if( true && pcSwitch>0. )
                                {
                  // *** add jump in normal traction ***
                                    uLocal(Ib1,Ib2,Ib3,uc+n) +=  (1./(zf+zp))*normal(Ib1,Ib2,Ib3,n)*( nSigmaSolidN - nSigmaFluidN );
                                }
                            }
                            if( debug() & 4 )
                            {
                // ::display(solidChar,"--INS-- IMP-AMP-V-BC: solid-characteristic: sigmas.nv + zs*vs (AGAIN)",debugFile,"%.3e ");
                                ::display(uLocal(Ig1,Ig2,Ig3,V),"--INS-- IMP-AMP-V-BC: RHS for GHOST",debugFile,"%.3e ");
                            }
              // ------ corners with dirichlet sides are treated specially -------
                            int iv[3], &i1=iv[0], &i2=iv[1], &i3=iv[2];
                            i3=gid(0,2); // set default
                            int axisp1= (axis+1) % numberOfDimensions;
                            for( int sidep1=0; sidep1<=1; sidep1++ )
                            {
                                if( mg.boundaryCondition(sidep1,axisp1)==Parameters::dirichletBoundaryCondition )
                                {
                                    if( numberOfDimensions==3 )
                                    {
                                        OV_ABORT("AMP IMP BC CORNER -- finish me for 3D");
                                    }
                                    for( int ghost=0; ghost<=1; ghost++ )  // set boundary and ghost values
                                    {
                                        iv[axis  ]=gid(side  ,axis  ) + ghost*(2*side-1);
                                        iv[axisp1]=gid(sidep1,axisp1);
                                        printf("Set AMP -IMP RHS for pt=(%i,%i,%i)\n",i1,i2,i3);
                                        uLocal(i1,i2,i3,V)=0.;
                                        if( twilightZoneFlow )
                                        {
                      // ------- FIX FOR TWILIGHT ZONE --------
                                            OGFunction & e = *(parameters.dbase.get<OGFunction* >("exactSolution"));
                                            const bool rectangular=false;
                                            Index J1=i1, J2=i2, J3=i3;
                                            OV_GET_SERIAL_ARRAY(real,mg.vertex(),xLocal);
                                            realSerialArray ue(J1,J2,J3,V);
                                            e.gd( ue  ,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,0,J1,J2,J3,V,t); // t=new time 
                                            uLocal(i1,i2,i3,V)=ue;
                                        }
                                    }
                                }
                            }
            // if( FALSE )
            // {
            //   // *** OLD WAY 
            //   real an[3], val[3], eqnVal[3];
            //   FOR_3IJD(i1,i2,i3,Ib1,Ib2,Ib3,ig1,ig2,ig3,Ig1,Ig2,Ig3)
            //   {
            //     if( numberOfDimensions==2 )
            //     {
            //       // ------ TWO DIMENSIONS -------
            //       for( int dir=0; dir<numberOfDimensions; dir++ )
            //       {
            //         an[dir]=normal(i1,i2,i3,dir);
            //         val[dir] = solidTraction(i1,i2,i3,dir) - zs*vSolidLocal(i1,i2,i3,dir);
            //         eqnVal[dir] = uOldLocal(i1,i2,i3,uc+dir) + 
            //           (-implicitFactor*dt/rho)*( px(i1,i2,i3) )
            //           + ((1.-implicitFactor)*dt/rho)*( 
            //             mu*(uOldxx(i1,i2,i3,uc+dir)+uOldyy(i1,i2,i3,uc+dir)) - pOldx(i1,i2,i3) );
            //       }
            //       for( int n=0; n<numberOfDimensions; n++ )
            //       {
            //         // --- RHS FOR GHOST EQUATION: ---
            //         //  RHS = (1/mu)*( I - n n^T)( solidTraction - zs*vs )
            //         u(ig1,ig2,ig3,uc+n) = (1./mu)*( AMP(n,0)*val[0] + AMP(n,1)*val[1] );
            //         // --- RHS FOR BOUNDARY EQUATION: ---
            //         u(i1,i2,i3,uc+n) = ( AMG(n,0)*eqnVal[0] + AMG(n,1)*eqnVal[1] + 
            //                              AMGP(n,0)*vSolidLocal(i1,i2,i3,0) + AMGP(n,1)*vSolidLocal(i1,i2,i3,1) );
            //       }
            //       OV_ABORT("finish me");
            //     }
            //     else
            //     {
            //       // ------ THREE DIMENSIONS -------
            //       for( int n=0; n<numberOfDimensions; n++ )
            //       {
            //         u(ig1,ig2,ig3,uc+n) = (1./mu)*( AMP(n,0)*val[0] + AMP(n,1)*val[1]+ AMP(n,2)*val[2] );
            //       }
            //       OV_ABORT("addedMassImplicitBoundaryConditions: 3D: finish me");
            //     }
            //   } // end for ijd
            // } // end if FALSE
                        }

                    }
                    
                }
            }
        }
    }
    

  // --- Now assign the RHS for adjacent no-slip walls ----
  // Note: do this last so corners are zero 
  // Note: Other boundary condition right-hand-sides are handled in applyBoundaryConditionsForImplicitTimeStepping
    if( applied )
    {
        const IntegerArray & interfaceType = parameters.dbase.get<IntegerArray >("interfaceType");
        ForBoundary(side,axis)
        {
            if( mg.boundaryCondition(side,axis)==Parameters::noSlipWall && interfaceType(side,axis,grid)==Parameters::noInterface )
            {
        // This must be an adajent no-slip wall
                if( false )
                    printF("WWWWW addedMassImpBC: apply HOMOGENEOUS noSlipWall BC to a DEFORMING BULK SOLID, grid=%i (side,axis)=(%i,%i), t=%.2e\n",
                                  grid,side,axis,t);
                u.applyBoundaryCondition(V,BCTypes::dirichlet,BCTypes::boundary(side,axis),0.,t);
            }
      // else if( mg.boundaryCondition(side,axis)>0  && interfaceType(side,axis,grid)==Parameters::noInterface )
      // {
      //   printF("INS: addedMassImpBC: ERROR: adjacent face to a deforming interface is not an noSlipWall -- fix me\n"
      //          " grid=%i (side,axis)=(%i,%i) bc=%i\n",grid,side,axis,mg.boundaryCondition(side,axis));
      //   OV_ABORT("error");
      // }
        
        }
    }
    

    return applied;
}







// ===============================================================================================
// Macro: checkAmpBoundaryConditions: check the solution from the implicit solve to see that
//   the AMP boundary conditions were all assigned properly
//
// The velocity BC's can be combined to give two vector equations:
//   (I)  (div(v))*n + (I-n n^T)(tauv.n - zs*v  )/mu = RHS    (Eqn for boundary pt)
//   (II)  v - theta*dt*nu*(I+P) Delta v = RHS                (Eqn for ghost pt)
//  where     P = (alpha-1) n n^T 
// ===============================================================================================




int Cgins::
checkAddedMassImplicitBoundaryConditions(realMappedGridFunction & u, 
                                                                                  realMappedGridFunction &uOld,  
                                                                                  real t,
                                                                                  int grid, 
                                                                                  real dt0 )
// ======================================================================================
// \brief Check the result from the implicit solve that the AMP BC's were properly assigned
///       for INS + BULK ELASTIC SOLID
///
/// \param u (input/output) : solution from solve
/// \param uOld (input/output) : solution from previous time
/// \param t (input) : time (new time)
/// \param grid (input) : component grid number.
///
// ==========================================================================================
{

    const bool & useAddedMassAlgorithm = parameters.dbase.get<bool>("useAddedMassAlgorithm");
    const bool & projectAddedMassVelocity = parameters.dbase.get<bool>("projectAddedMassVelocity");
    const int initialConditionsAreBeingProjected = parameters.dbase.get<int>("initialConditionsAreBeingProjected");
    const bool & useImplicitAmpBCs = parameters.dbase.get<bool>("useImplicitAmpBCs");
    const bool & predictedBoundaryPressureNeeded = parameters.dbase.get<bool>("predictedBoundaryPressureNeeded");

    bool applyAddedMass = ( useAddedMassAlgorithm && 
                                                    projectAddedMassVelocity && 
                                                    parameters.gridIsMoving(grid)
                                                    && !initialConditionsAreBeingProjected 
                                                    && t!=0. );

  // For testing with Cgins we may fill in the implicit AMP BCs even when running Cgins alone
    applyAddedMass = applyAddedMass || (useAddedMassAlgorithm && useImplicitAmpBCs);
    

    if( !applyAddedMass )
    {
        return 0;
    } 
    
    if( FALSE )
    {
        u.periodicUpdate(); // **TEST**
    }
    

    const int correctionStage = parameters.dbase.get<int>("correctionStage");

    MovingGrids & movingGrids = parameters.dbase.get<MovingGrids >("movingGrids");
    if( movingGrids.getNumberOfDeformingBodies()==0 && !useImplicitAmpBCs )
        return 0;

    real dt=dt0;
    if( dt<= 0. )
        dt = parameters.dbase.get<real>("dt");  // *wdh* 2017/05/31
  // const real & dt = parameters.dbase.get<real>("dt");
    const real epsT = REAL_EPSILON*10.;
    assert( dt> epsT );
    
    FILE *&debugFile = parameters.dbase.get<FILE* >("debugFile");
    FILE *&pDebugFile = parameters.dbase.get<FILE* >("pDebugFile");

    if( t <= 2.*dt && debug() &4 )
    {
        fPrintF(debugFile,"\n"
                      "--------------------------------------------------------------------------------------------------\n"
                      " --INS-- checkAddedMassImplicitBoundaryConditions: ADDED MASS ALGORITHM t=%8.2e\n"
                      "--------------------------------------------------------------------------------------------------\n"
                        ,t);
    }
    
    const real nu = parameters.dbase.get<real>("nu");
    const real rho = parameters.dbase.get<real>("rho");
    const real mu = rho*nu;
    const real implicitFactor = parameters.dbase.get<real >("implicitFactor");
    const real theta = implicitFactor;

    const bool & projectNormalComponentOfAddedMassVelocity =
                              parameters.dbase.get<bool>("projectNormalComponentOfAddedMassVelocity");
    const bool & projectVelocityOnBeamEnds = parameters.dbase.get<bool>("projectVelocityOnBeamEnds"); 

    const bool twilightZoneFlow = parameters.dbase.get<bool >("twilightZoneFlow");

    assert(  useAddedMassAlgorithm && projectAddedMassVelocity && (useImplicitAmpBCs || parameters.gridIsMoving(grid)) );
    
    MappedGrid & mg = *u.getMappedGrid();
    const IntegerArray & gid = mg.gridIndexRange();
    const int numberOfDimensions = mg.numberOfDimensions();
    Range Rx=numberOfDimensions;

    const bool gridIsMoving = parameters.gridIsMoving(grid);

  // get components of solution
    const int uc = parameters.dbase.get<int >("uc");
    const int vc = parameters.dbase.get<int >("vc");
    const int wc = parameters.dbase.get<int >("wc");
    const int tc = parameters.dbase.get<int >("tc");
    const int pc = parameters.dbase.get<int >("pc");
    const int & nc = parameters.dbase.get<int >("nc");

    const int orderOfAccuracy=min(4,parameters.dbase.get<int >("orderOfAccuracy"));
    Range V = Range(uc,uc+numberOfDimensions-1);

    BoundaryData::BoundaryDataArray & pBoundaryData = parameters.getBoundaryData(grid); // this will create the BDA if it is not there
    std::vector<BoundaryData> & boundaryDataArray =parameters.dbase.get<std::vector<BoundaryData> >("boundaryData");
    BoundaryData & bd = boundaryDataArray[grid];
            
    const Parameters::InterfaceCommunicationModeEnum & interfaceCommunicationMode= 
        parameters.dbase.get<Parameters::InterfaceCommunicationModeEnum>("interfaceCommunicationMode");

    const Parameters::KnownSolutionsEnum & knownSolution = 
        parameters.dbase.get<Parameters::KnownSolutionsEnum >("knownSolution"); 

    if( !parameters.gridIsMoving(grid) )
    {
    // ---- this is a test run using non-moving grids ----
        ForBoundary(side,axis)
        {
            if( mg.boundaryCondition(side,axis)==Parameters::noSlipWall )
            {
                printF(" XXXX addedMassImplicitBoundaryConditions: TEST RUN using non-moving grids XXXX \n\n");


                Index Ib1,Ib2,Ib3, Ig1,Ig2,Ig3, Ip1,Ip2,Ip3;
                getBoundaryIndex(mg.gridIndexRange(),side,axis,Ib1,Ib2,Ib3);
                getGhostIndex(mg.gridIndexRange(),side,axis,Ig1,Ig2,Ig3,+1);  // first ghost line 

                RealArray solidTraction(Ib1,Ib2,Ib3,numberOfDimensions), vSolidLocal(Ib1,Ib2,Ib3,numberOfDimensions);

        // Do this for now
                real zf=1., zp=1., zs=1., alpha=.5; // defaults when tetsing cgins alone
                solidTraction=0.;
                vSolidLocal=0.;

                OV_GET_SERIAL_ARRAY(real,u   ,uLocal);
                OV_GET_SERIAL_ARRAY(real,uOld,uOldLocal);
                
                if( projectNormalComponentOfAddedMassVelocity )
                    mg.update(MappedGrid::THEvertexBoundaryNormal);

                OV_GET_VERTEX_BOUNDARY_NORMAL(mg,side,axis,normal);         

        // evaluate the right-hand-side to the AMP velocity BCs
                if( t <= 3.*dt && debug() & 1 )
                {
          // Macro: 
                    {
                        if( debug() & 4 || t<=3.*dt ) 
                        {
                            fPrintF(debugFile,"--INS-- CHECK-AMP-BC-IMP FOR BULK SOLID MODEL, t=%9.3e dt=%9.3e zf=%9.3e zp=%9.3e zs=%9.3e "
                                          "alpha=%9.2e **new**\n",t,dt,zf,zp,zs,alpha);
                            printF("\n--INS-- CHECK-AMP-BC-IMP FOR BULK SOLID MODEL, t=%9.3e dt=%9.3e zf=%9.3e zp=%9.3e zs=%9.3e "
                                          "alpha=%9.2e **new**\n",t,dt,zf,zp,zs,alpha);
                        }
                        if( FALSE )
                        {
                            uOld.periodicUpdate(); // **TEST***
                        }
                        const real beta = (mu*theta*dt/rho);
            // -- evaluate derivatives that appear in the INS equations ---       
                        MappedGridOperators & opOld = *(uOld.getOperators()); 
                        realSerialArray uOldxx(Ib1,Ib2,Ib3,V), uOldyy(Ib1,Ib2,Ib3,V);
                        realSerialArray pOldx(Ib1,Ib2,Ib3), pOldy(Ib1,Ib2,Ib3), pOldz;
                        opOld.derivative(MappedGridOperators::xxDerivative ,uOldLocal,uOldxx,Ib1,Ib2,Ib3,V);
                        opOld.derivative(MappedGridOperators::yyDerivative ,uOldLocal,uOldyy,Ib1,Ib2,Ib3,V);
                        opOld.derivative(MappedGridOperators::xDerivative ,uOldLocal,pOldx,Ib1,Ib2,Ib3,pc);
                        opOld.derivative(MappedGridOperators::yDerivative ,uOldLocal,pOldy,Ib1,Ib2,Ib3,pc);
                        if( numberOfDimensions==3 )
                        {
                            pOldz.redim(Ib1,Ib2,Ib3);
                            opOld.derivative(MappedGridOperators::zDerivative ,uOldLocal,pOldz,Ib1,Ib2,Ib3,pc);
                        }
            // -- We assume a predicted value for "p" at the new time
                        MappedGridOperators & op = *(u.getOperators()); 
                        realSerialArray px(Ib1,Ib2,Ib3), py(Ib1,Ib2,Ib3),pz;
                        op.derivative(MappedGridOperators::xDerivative ,uLocal,px,Ib1,Ib2,Ib3,pc);
                        op.derivative(MappedGridOperators::yDerivative ,uLocal,py,Ib1,Ib2,Ib3,pc);
                        if( numberOfDimensions==3 )
                        {
                            pz.redim(Ib1,Ib2,Ib3);
                            op.derivative(MappedGridOperators::zDerivative ,uLocal,pz,Ib1,Ib2,Ib3,pc);
                        }
                        realSerialArray ux(Ib1,Ib2,Ib3,V), uy(Ib1,Ib2,Ib3,V),uz;
                        op.derivative(MappedGridOperators::xDerivative ,uLocal,ux,Ib1,Ib2,Ib3,V);
                        op.derivative(MappedGridOperators::yDerivative ,uLocal,uy,Ib1,Ib2,Ib3,V);
                        if( numberOfDimensions==3 )
                        {
                            uz.redim(Ib1,Ib2,Ib3,V);
                            op.derivative(MappedGridOperators::zDerivative ,uLocal,uz,Ib1,Ib2,Ib3,V);
                        }
                        RealArray div(Ib1,Ib2,Ib3);
                        div = ux(Ib1,Ib2,Ib3,uc)+uy(Ib1,Ib2,Ib3,vc);
                        if( numberOfDimensions==3 )
                            div += uz(Ib1,Ib2,Ib3,wc);
                        real maxDiv = max(fabs(div));
            // printF("--INS--check-AMP-IMP-BC: max divergence on boundary =%.2e (dt=%.2e)\n",maxDiv,dt);
                        fPrintF(debugFile,"--INS--check-AMP-IMP-BC: max divergence on boundary =%.2e (dt=%.2e)\n",maxDiv,dt);
                        if( debug() & 4 )
                            ::display(div,"Divergence on the boundary","%.2e ");
                        RealArray delta(3,3); // Dirac delta
                        delta=0.; delta(0,0)=1.; delta(1,1)=1.; delta(2,2)=1.;
                        RealArray solidChar(Ib1,Ib2,Ib3,Rx);   // solid characteristic variable 
                        RealArray fluidRhs(Ib1,Ib2,Ib3,Rx);    // fluid interior update without implicit viscous term
            // Solid characteristic variable used in the RHS 
            // Use MINUS of solid traction since the solid normal is in the opposite direction to the fluid
                        solidChar(Ib1,Ib2,Ib3,Rx) = -solidTraction(Ib1,Ib2,Ib3,Rx) + zs*vSolidLocal(Ib1,Ib2,Ib3,Rx);
            // fluidRhs(Ib1,Ib2,Ib3,dir) = fluid interior update without implicit viscous term 
                        for( int dir=0; dir<numberOfDimensions; dir++ )
                        {
              // pd = px, py or pz 
                            RealArray & pd    = dir==0 ? px    : dir==1? py    : pz;
                            RealArray & pOldd = dir==0 ? pOldx : dir==1? pOldy : pOldz;
              // fluidRhs = v^{n-1} - (theta*dt/rho)*grad(p)^n + ((1-theta)*dt/rho)*( -grad(p)^{n-1} + mu*Delta(v)^{n-1} )
                            fluidRhs(Ib1,Ib2,Ib3,dir)= uOldLocal(Ib1,Ib2,Ib3,uc+dir) 
                                + (-implicitFactor*dt/rho)*( pd(Ib1,Ib2,Ib3) )
                                + ((1.-implicitFactor)*dt/rho)*( 
                                            mu*( uOldxx(Ib1,Ib2,Ib3,uc+dir)+uOldyy(Ib1,Ib2,Ib3,uc+dir) ) - pOldd(Ib1,Ib2,Ib3) 
                                                                                              );
                        }
                        if( twilightZoneFlow )
                        {
              // ------- ADD TWILIGHT ZONE --------
                            OGFunction & e = *(parameters.dbase.get<OGFunction* >("exactSolution"));
                            const bool rectangular=false;
                            OV_GET_SERIAL_ARRAY(real,mg.vertex(),xLocal);
                            if( debug() & 4 )
                            {
                                RealArray pex(Ib1,Ib2,Ib3), pey(Ib1,Ib2,Ib3);
                                e.gd( pex,xLocal,mg.numberOfDimensions(),rectangular,0,1,0,0,Ib1,Ib2,Ib3,pc,t); // t=new time 
                                e.gd( pey,xLocal,mg.numberOfDimensions(),rectangular,0,0,1,0,Ib1,Ib2,Ib3,pc,t); // t=new time 
                                real maxErrPx = max(fabs(px-pex));
                                real maxErrPy = max(fabs(py-pey));
                                fprintf(debugFile,">>>INS IMP-AMP-V-BC: error in predicted (px,py)=(%9.3e,%9.3e) at t=%9.3e\n",maxErrPx,maxErrPy,t);
                                ::display(fabs(px-pex),"error in px",debugFile,"%8.2e ");
                                ::display(fabs(py-pey),"error in py",debugFile,"%8.2e ");
                            }
                            realSerialArray ue(Ib1,Ib2,Ib3,V), uex(Ib1,Ib2,Ib3,V), uey(Ib1,Ib2,Ib3,V);
                            RealArray tractione(Ib1,Ib2,Ib3,Rx);
                            e.gd( ue  ,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,0,Ib1,Ib2,Ib3,V,t); // t=new time 
                            e.gd( uex ,xLocal,mg.numberOfDimensions(),rectangular,0,1,0,0,Ib1,Ib2,Ib3,V,t); 
                            e.gd( uey ,xLocal,mg.numberOfDimensions(),rectangular,0,0,1,0,Ib1,Ib2,Ib3,V,t); 
              // tau_ij = mu*( partial_j u_i + partial_i u_j )
                            if( numberOfDimensions==2 )
                            {
                // ====== Get traction in two dimensions =====
                                RealArray taue(Ib1,Ib2,Ib3,3); // hold viscous stress tensor tau
                                taue(Ib1,Ib2,Ib3,0) = (2.*mu)*( uex(Ib1,Ib2,Ib3,uc) );                        // tau_11
                                taue(Ib1,Ib2,Ib3,1) = (   mu)*( uey(Ib1,Ib2,Ib3,uc) + uex(Ib1,Ib2,Ib3,vc) );  // tau_12 = tau_21
                                taue(Ib1,Ib2,Ib3,2) = (2.*mu)*( uey(Ib1,Ib2,Ib3,vc) );                        // tau_22
                //  t1 = tau11*n1 + tau12*n2 
                                tractione(Ib1,Ib2,Ib3,0) = taue(Ib1,Ib2,Ib3,0)*normal(Ib1,Ib2,Ib3,0) + taue(Ib1,Ib2,Ib3,1)*normal(Ib1,Ib2,Ib3,1);
                //  t2 = tau21*n1 + tau22*n2 
                                tractione(Ib1,Ib2,Ib3,1) = taue(Ib1,Ib2,Ib3,1)*normal(Ib1,Ib2,Ib3,0) + taue(Ib1,Ib2,Ib3,2)*normal(Ib1,Ib2,Ib3,1);
                            }
                            else 
                            {
                // ====== Get traction in three dimensions =====
                                RealArray uez(Ib1,Ib2,Ib3,V);
                                e.gd( uez ,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,1,Ib1,Ib2,Ib3,V,t); 
                                RealArray taue(Ib1,Ib2,Ib3,6);
                                taue(Ib1,Ib2,Ib3,0) = (2.*mu)*( uex(Ib1,Ib2,Ib3,uc) );                        // tau_11
                                taue(Ib1,Ib2,Ib3,1) = (   mu)*( uey(Ib1,Ib2,Ib3,uc) + uex(Ib1,Ib2,Ib3,vc) );  // tau_12 = tau_21
                                taue(Ib1,Ib2,Ib3,2) = (   mu)*( uez(Ib1,Ib2,Ib3,uc) + uex(Ib1,Ib2,Ib3,wc) );  // tau_13 = tau_31
                                taue(Ib1,Ib2,Ib3,3) = (2.*mu)*( uey(Ib1,Ib2,Ib3,vc) );                        // tau_22
                                taue(Ib1,Ib2,Ib3,4) = (   mu)*( uez(Ib1,Ib2,Ib3,vc) + uey(Ib1,Ib2,Ib3,wc) );  // tau_23 = tau_32
                                taue(Ib1,Ib2,Ib3,5) = (2.*mu)*( uez(Ib1,Ib2,Ib3,wc) );                        // tau_33
                //  t1 = tau11*n1 + tau12*n2 + tau13*n3 
                                tractione(Ib1,Ib2,Ib3,0) = ( taue(Ib1,Ib2,Ib3,0)*normal(Ib1,Ib2,Ib3,0) + 
                                                                                          taue(Ib1,Ib2,Ib3,1)*normal(Ib1,Ib2,Ib3,1) +
                                                                                          taue(Ib1,Ib2,Ib3,2)*normal(Ib1,Ib2,Ib3,2) );
                //  t2 = tau21*n1 + tau22*n2 + tau23*n3 
                                tractione(Ib1,Ib2,Ib3,1) = ( taue(Ib1,Ib2,Ib3,1)*normal(Ib1,Ib2,Ib3,0) + 
                                                                                          taue(Ib1,Ib2,Ib3,3)*normal(Ib1,Ib2,Ib3,1) +
                                                                                          taue(Ib1,Ib2,Ib3,4)*normal(Ib1,Ib2,Ib3,2) );
                //  t3 = tau31*n1 + tau32*n2 + tau33*n3 
                                tractione(Ib1,Ib2,Ib3,2) = ( taue(Ib1,Ib2,Ib3,2)*normal(Ib1,Ib2,Ib3,0) + 
                                                                                          taue(Ib1,Ib2,Ib3,4)*normal(Ib1,Ib2,Ib3,1) +
                                                                                          taue(Ib1,Ib2,Ib3,5)*normal(Ib1,Ib2,Ib3,2) );
                            }
                            solidChar(Ib1,Ib2,Ib3,Rx) = tractione(Ib1,Ib2,Ib3,Rx) + zs*ue(Ib1,Ib2,Ib3,V);
              // fluidRhs = v^{n-1} + (theta*dt/rho)*grad(p)^n + ((1-theta)*dt/rho)*( -grad(p)^{n-1} + mu*Delta(v)^{n-1} )
              //          = [ I - (mu*theta*dt/rho)*Delta ]* v^n 
                            RealArray & uexx = uex; RealArray & ueyy =uey;  // Reuse arrays
                            e.gd( uexx ,xLocal,mg.numberOfDimensions(),rectangular,0,2,0,0,Ib1,Ib2,Ib3,V,t); 
                            e.gd( ueyy ,xLocal,mg.numberOfDimensions(),rectangular,0,0,2,0,Ib1,Ib2,Ib3,V,t); 
                            vSolidLocal(Ib1,Ib2,Ib3,Rx) = ue(Ib1,Ib2,Ib3,V);
                            if( numberOfDimensions==2 )
                            {
                //  [ I - (mu*theta*dt/rho)*Delta ]* v^n
                                fPrintF(debugFile,"addedMassImpBC: beta=%9.3e, max(uexx=%8.2e)\n",beta,max(fabs(uexx)));
                                fluidRhs(Ib1,Ib2,Ib3,Rx)=ue(Ib1,Ib2,Ib3,V) - beta*( uexx+ueyy ); 
                            }
                            else
                            {
                                RealArray & uezz = tractione; // Reuse arrays
                                e.gd( uezz ,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,2,Ib1,Ib2,Ib3,V,t); 
                //  [ I - (mu*theta*dt/rho)*Delta ]* v^n
                // *** THIS IS WRONG ***
                                fluidRhs(Ib1,Ib2,Ib3,Rx)=ue(Ib1,Ib2,Ib3,V) - (mu*theta*dt/rho)*( uexx+ueyy+uezz ); 
                            }
                        }
                        RealArray tau(Ib1,Ib2,Ib3,3); // holds viscous stress tensor tau
                        RealArray traction(Ib1,Ib2,Ib3,3);
                        RealArray res(Ib1,Ib2,Ib3);
                        if( numberOfDimensions==2 )  
                        {
                            tau(Ib1,Ib2,Ib3,0) = (2.*mu)*( ux(Ib1,Ib2,Ib3,uc) );                        // tau_11
                            tau(Ib1,Ib2,Ib3,1) = (   mu)*( uy(Ib1,Ib2,Ib3,uc) + ux(Ib1,Ib2,Ib3,vc) );   // tau_12 = tau_21
                            tau(Ib1,Ib2,Ib3,2) = (2.*mu)*( uy(Ib1,Ib2,Ib3,vc) );                        // tau_22
              //  t1 = tau11*n1 + tau12*n2 
                            traction(Ib1,Ib2,Ib3,0) = tau(Ib1,Ib2,Ib3,0)*normal(Ib1,Ib2,Ib3,0) + tau(Ib1,Ib2,Ib3,1)*normal(Ib1,Ib2,Ib3,1);
              //  t2 = tau21*n1 + tau22*n2 
                            traction(Ib1,Ib2,Ib3,1) = tau(Ib1,Ib2,Ib3,1)*normal(Ib1,Ib2,Ib3,0) + tau(Ib1,Ib2,Ib3,2)*normal(Ib1,Ib2,Ib3,1);
              // Check BC:
              //    tv.[  tau.nv + zs*v - (sigmaSolid.nv + zs*vSolid ) ]=0 
              // In 2D the tangent-vector is tv = [ -n_y, n_x ]
                            res = ( -normal(Ib1,Ib2,Ib3,1)*( traction(Ib1,Ib2,Ib3,0) + zs*uLocal(Ib1,Ib2,Ib3,uc) - solidChar(Ib1,Ib2,Ib3,0) )
                                            +normal(Ib1,Ib2,Ib3,0)*( traction(Ib1,Ib2,Ib3,1) + zs*uLocal(Ib1,Ib2,Ib3,vc) - solidChar(Ib1,Ib2,Ib3,1) ) );
                        }
                        else
                        {
                            OV_ABORT("finish me for 3D");
                        }
                        real maxChar = max(fabs(res));
                        fPrintF(debugFile,"--INS--check-AMP-IMP-BC: max error in tv.( tau.nv + zs*v )= tv.solidChar =%.2e\n",maxChar);
                        if( debug() & 8 )
                            ::display(res,"Residual in  tv.( tau.nv + zs*v )= ... on the boundary",debugFile,"%.2e ");
            // ============= Check tangential component of interior eqn on the boundary =============
                        realSerialArray uxx(Ib1,Ib2,Ib3,V), uyy(Ib1,Ib2,Ib3,V), uzz;
                        op.derivative(MappedGridOperators::xxDerivative ,uLocal,uxx,Ib1,Ib2,Ib3,V);
                        op.derivative(MappedGridOperators::yyDerivative ,uLocal,uyy,Ib1,Ib2,Ib3,V);
                        if( numberOfDimensions==3 )
                        {
                            uzz.redim(Ib1,Ib2,Ib3,V);
                            op.derivative(MappedGridOperators::zzDerivative ,uLocal,uzz,Ib1,Ib2,Ib3,V);
                        }
                        RealArray unp1(Ib1,Ib2,Ib3,V);
                        unp1(Ib1,Ib2,Ib3,V) = fluidRhs(Ib1,Ib2,Ib3,Rx) + (implicitFactor*mu*dt/rho)*( uxx + uyy );
                        if( numberOfDimensions==2 )  
                        {
              // In 2D the tangent-vector is tv = [ -n_y, n_x ]
                            res = ( -normal(Ib1,Ib2,Ib3,1)*( unp1(Ib1,Ib2,Ib3,uc)-uLocal(Ib1,Ib2,Ib3,uc) )
                                            +normal(Ib1,Ib2,Ib3,0)*( unp1(Ib1,Ib2,Ib3,vc)-uLocal(Ib1,Ib2,Ib3,vc) ) );
                        }
                        real maxRes = max(fabs(res));
                        fPrintF(debugFile,"--INS--check-AMP-IMP-BC: max error in tv.( Eqn on boundary ) =%.2e\n",maxRes);
                        if( debug() & 8 )
                            ::display(res,"Residual in  tv.( Eqn on boundary )= ... on the boundary",debugFile,"%.2e ");
            // ============= Check normal component AMP condition =============
                        if( numberOfDimensions==2 )  
                        {
                            res = ( normal(Ib1,Ib2,Ib3,0)*( uLocal(Ib1,Ib2,Ib3,uc) - (alpha*unp1(Ib1,Ib2,Ib3,uc)+(1.-alpha)*vSolidLocal(Ib1,Ib2,Ib3,0)))
                                          +normal(Ib1,Ib2,Ib3,1)*( uLocal(Ib1,Ib2,Ib3,vc) - (alpha*unp1(Ib1,Ib2,Ib3,vc)+(1.-alpha)*vSolidLocal(Ib1,Ib2,Ib3,1))) );
                            if( correctionStage>0 )
                            {
                // jump in normal traction 
                                RealArray nSigmaSolidN, nSigmaFluidN;
                // **use minus for solidTraction **
                                nSigmaSolidN= -(normal(Ib1,Ib2,Ib3,0)*solidTraction(Ib1,Ib2,Ib3,0) +normal(Ib1,Ib2,Ib3,1)*solidTraction(Ib1,Ib2,Ib3,1) );
                // n.sigmaFluid.n = - p + n.Tau.n 
                                nSigmaFluidN = -uLocal(Ib1,Ib2,Ib3,pc) +
                                    normal(Ib1,Ib2,Ib3,0)*traction(Ib1,Ib2,Ib3,0) + normal(Ib1,Ib2,Ib3,1)*traction(Ib1,Ib2,Ib3,1);
                // ** fix me for TZ **
                                real tractionJump = max(fabs( nSigmaSolidN - nSigmaFluidN));
                                fPrintF(debugFile,"Correction stage: Max jump in traction = %.2e\n",tractionJump);
                                res -=  (1./(zf+zp))*( nSigmaSolidN - nSigmaFluidN );
                            }
                        }
                        else
                        {
                            OV_ABORT("finish me...");
                        }
                        maxRes = max(fabs(res));
                        fPrintF(debugFile,"--INS--check-AMP-IMP-BC: max error in nv.( vI - [alpha*vf + (1-alpha)*vs ] ) =%.2e\n",maxRes);
                        if( debug() & 8 )
                            ::display(res,"Residual in  nv.( vI - [alpha*vf + (1-alpha)*vs ] )",debugFile,"%.2e ");
                        if( false )
                        {
                            ::display(px,"px","%.9e ");
                            ::display(py,"py","%.9e ");
                            ::display(solidChar(Ib1,Ib2,Ib3,Rx),"solidChar","%.9e ");
                            ::display(fluidRhs(Ib1,Ib2,Ib3,Rx),"fluidRhs","%.9e ");
                            ::display(uxx(Ib1,Ib2,Ib3,V),"uxx","%.9e ");
                            ::display(uyy(Ib1,Ib2,Ib3,V),"uyy","%.9e ");
                            ::display(unp1(Ib1,Ib2,Ib3,V),"unp1","%.9e ");
                            ::display(uLocal(Ib1,Ib2,Ib3,V),"uLocal","%.9e ");
                            ::display(vSolidLocal(Ib1,Ib2,Ib3,Rx),"vSolid","%.9e ");
                            ::display(normal(Ib1,Ib2,Ib3,Rx),"normalVector","%.6e ");
                            const realArray & rx = mg.inverseVertexDerivative();
                            Range all;
                            ::display(rx(Ib1,Ib2,Ib3,all),"rx(Ib1,Ib2,Ib3,all)","%.6e ");
                            OV_ABORT("stop here for now");
                        }
          //  -- "vector" versions of macros that use full normal vector : 
                    #define AMPV(n,m) (delta(n,m) - normal(Ib1,Ib2,Ib3,n)*normal(Ib1,Ib2,Ib3,m) )
          // Here is the "P" operator, (1-alpha) n n^T 
                    #define AMGPV(n,m) ((1.-alpha)*normal(Ib1,Ib2,Ib3,n)*normal(Ib1,Ib2,Ib3,m))
          // Here is I-P 
                    #define AMGV(n,m) (delta(n,m) - (1.-alpha)*normal(Ib1,Ib2,Ib3,n)*normal(Ib1,Ib2,Ib3,m))
            // for( int n=0; n<numberOfDimensions; n++ )
            // {
            //   // --- RHS FOR GHOST EQUATION: ---
            //   //  RHS = (1/mu)*( I - n n^T)( solidTraction + zs*vs )
            //   uLocal(Ig1,Ig2,Ig3,uc+n) = (1./mu)*( AMPV(n,0)*solidChar(Ib1,Ib2,Ib3,0) + AMPV(n,1)*solidChar(Ib1,Ib2,Ib3,1) );
            //   // --- RHS FOR BOUNDARY EQUATION: ---
            //   if( false )
            //   {
            //     // TEST
            //   }
            //   else
            //   {
            //     // printF("AM-IMP_BC: n=%i,  AMGV=%9.3e %9.3e %9.3e,  AMGPV=%9.3e %9.3e\n",n,
            //     // RHS = (I-P)*( v - beta*Delta(v) ) + P*vs 
            //     uLocal(Ib1,Ib2,Ib3,uc+n) = ( AMGV(n,0) *fluidRhs(Ib1,Ib2,Ib3,0)    + AMGV(n,1) *fluidRhs(Ib1,Ib2,Ib3,1) + 
            //                                  AMGPV(n,0)*vSolidLocal(Ib1,Ib2,Ib3,0) + AMGPV(n,1)*vSolidLocal(Ib1,Ib2,Ib3,1) );
            //   }
            // }
            // // ------ corners with dirichlet sides are treated specially -------
            // int iv[3], &i1=iv[0], &i2=iv[1], &i3=iv[2];
            // i3=gid(0,2); // set default
            // int axisp1= (axis+1) % numberOfDimensions;
            // for( int sidep1=0; sidep1<=1; sidep1++ )
            // {
            //   if( mg.boundaryCondition(sidep1,axisp1)==Parameters::dirichletBoundaryCondition )
            //   {
            //     if( numberOfDimensions==3 )
            //     {
            //       OV_ABORT("AMP IMP BC CORNER -- finish me for 3D");
            //     }
            //     for( int ghost=0; ghost<=1; ghost++ )  // set boundary and ghost values
            //     {
            //       iv[axis  ]=gid(side  ,axis  ) + ghost*(2*side-1);
            //       iv[axisp1]=gid(sidep1,axisp1);
            //       printf("Set AMP -IMP RHS for pt=(%i,%i,%i)\n",i1,i2,i3);
            //       uLocal(i1,i2,i3,V)=0.;
            //       if( twilightZoneFlow )
            //       {
            //         // ------- FIX FOR TWILIGHT ZONE --------
            //         OGFunction & e = *(parameters.dbase.get<OGFunction* >("exactSolution"));
            //         const bool rectangular=false;
            //         Index J1=i1, J2=i2, J3=i3;
            //         OV_GET_SERIAL_ARRAY(real,mg.vertex(),xLocal);
            //         realSerialArray ue(J1,J2,J3,V);
            //         e.gd( ue  ,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,0,J1,J2,J3,V,t); // t=new time 
            //         uLocal(i1,i2,i3,V)=ue;
            //       }
            //     }
            //   }
            // }
                    }
                }
                
                
            }
            
        }
        
    // OV_ABORT("addedMassImplicitBoundaryConditions: STOP HERE FOR NOW");
        
        return 0;
    }
    

  // ========== REAL DEFORMING GRID CASE ==========

  // -- extract parameters from any deforming solids ---
    
    if( bd.dbase.has_key("deformingBodyNumber") )
    {
        const real & fluidDensity = parameters.dbase.get<real >("fluidDensity");
        assert( fluidDensity>0. );

        const real fluidAddedMassLengthScale =  parameters.dbase.get<real>("fluidAddedMassLengthScale");

        int (&deformingBodyNumber)[2][3] = bd.dbase.get<int[2][3]>("deformingBodyNumber");
        Index Ib1,Ib2,Ib3, Ig1,Ig2,Ig3, Ip1,Ip2,Ip3;
        for( int side=0; side<=1; side++ )
        {
            for( int axis=0; axis<numberOfDimensions; axis++ )
            {
      	if( deformingBodyNumber[side][axis]>=0 )
      	{
        	  int body=deformingBodyNumber[side][axis];
        	  if( t<=0. )
          	    printF("--INS-- grid=%i, (side,axis)=(%i,%i) belongs to deforming body %i\n",grid,side,axis,body);

        	  DeformingBodyMotion & deform = movingGrids.getDeformingBody(body);

        	  getBoundaryIndex(mg.gridIndexRange(),side,axis,Ib1,Ib2,Ib3);
                    getGhostIndex(mg.gridIndexRange(),side,axis,Ig1,Ig2,Ig3,+1);  // first ghost line 
                    getGhostIndex(mg.gridIndexRange(),side,axis,Ip1,Ip2,Ip3,-1);  // first line in

        	  realArray vSolid(Ib1,Ib2,Ib3,Rx); // holds velocity of solid on the boundary
                    #ifndef USE_PPP
                    deform.getVelocityBC( t, side,axis,grid, mg, Ib1,Ib2,Ib3, vSolid );
                    #else
                        OV_ABORT("finish me");
                    #endif

          // OV_GET_SERIAL_ARRAY(real,gridVelocity,gridVelocityLocal);
        	  OV_GET_SERIAL_ARRAY(real,u   ,uLocal);
                    OV_GET_SERIAL_ARRAY(real,uOld,uOldLocal);
        	  OV_GET_SERIAL_ARRAY(real,vSolid,vSolidLocal);

                    if( projectNormalComponentOfAddedMassVelocity )
                        mg.update(MappedGrid::THEvertexBoundaryNormal);

                    OV_GET_VERTEX_BOUNDARY_NORMAL(mg,side,axis,normal);
        	  
                    if( !deform.isBulkSolidModel() && !deform.isBeamModel() )
        	  {

          	    OV_ABORT("addedMassImplicitBoundaryConditions::ERROR: un-expected deformation type");

        	  }
        	  else if( deform.isBulkSolidModel() )
        	  {

            // **********************************************************************
            // *************** ADDED MASS BC FOR A BULK SOLID ***********************
            // **********************************************************************

                        RealArray solidTraction;
                        if( interfaceCommunicationMode==Parameters::requestInterfaceDataWhenNeeded )
                        {
                            if( debug() & 4 )
                                fPrintF(debugFile,"--INS-- PIV: REQUEST interface traction at t=%9.3e\n",t);

                            InterfaceData interfaceData;
                            Range Rx=numberOfDimensions;
                            interfaceData.u.redim(Ib1,Ib2,Ib3,Rx); // traction is returned here 
                            interfaceData.t=t;
                            interfaceData.u=0;

                            int interfaceDataOptions = Parameters::tractionInterfaceData;
                            bool saveTimeHistory=true;
                            getInterfaceData( t, grid, side, axis, 
                                                                interfaceDataOptions,
                                                                interfaceData.u,
                                                                parameters,saveTimeHistory );

                            solidTraction=interfaceData.u;
                            if( t <= 2.*dt && debug() & 4 ) {
                                ::display(solidTraction(Ib1,Ib2,Ib3,1),"--INS-- AM-IMP-BC: Here is the SOLID TRACTION (I1,I2,I3,1)",debugFile,"%6.3f ");
                            }
                        }
                        else
                        {
                            OV_ABORT("finish me");
                        }

            // *new* April 1, 2018
                        real zs,zp,zf,alpha;
                        getBulkSolidAmpParameters( mg,grid,side,axis,dt, zs,zp,zf,alpha );

	    // real zpOld;
            // // old way:
            // deform.getBulkSolidParameters( zpOld );

            // // new way:
            // // Retrieve the parameters from the bulk solid
            // // FIX ME -- lookup first time and then save locally 
            // Parameters & bulkSolidParams = getInterfaceParameters( grid,side,axis,parameters );
            // real rhoSolid=bulkSolidParams.dbase.get<real>("rho");
            // real lambdaSolid=bulkSolidParams.dbase.get<real>("lambda");
            // real muSolid=bulkSolidParams.dbase.get<real>("mu");
            // real cp=sqrt((lambdaSolid+2.*muSolid)/rhoSolid);
            // real cs=sqrt(muSolid/rhoSolid);
                    

            // real zp=rhoSolid*cp;
            // real zs=rhoSolid*cs;
            // if( t<=3.*dt && debug() & 4  )
            // {
            //   fPrintF(debugFile,"--INS-- AMP-IMP-BC: rhoSolid=%9.3e cp=%9.3e cs=%9.3e zp=%9.3e (old: zp=%9.3e) zs=%9.3e\n",
            //          rhoSolid,cp,cs,zp,zpOld,zs);
            //   // printF("  fluidAddedMassLengthScale=%9.3e\n",fluidAddedMassLengthScale);
            // }
                        
	    // // const real & fluidDensity = parameters.dbase.get<real >("fluidDensity");
          	    
            // // fluid impedance = rho*H/dt 
	    // assert( dt>0. );
            // const real zf=fluidDensity*fluidAddedMassLengthScale/dt; 
            // const real alpha = zf/(zf+zp);


                        bool useExactSolidVelocity=false;
                        bool useExactSolidTraction=false;
                        
                        if( useExactSolidVelocity && knownSolution==Parameters::userDefinedKnownSolution )
                        {
                            printF("\n--INS-- AMP-IMP-BC: **TEST** set exact KNOWN-SOLUTION values for solid velocity ***TEMP***, t=%9.3e\n\n",t);
                            int body=0;
                            parameters.getUserDefinedDeformingBodyKnownSolution( body,Parameters::boundaryVelocity,
                                                                                                                                      t, grid, mg, Ib1,Ib2,Ib3,Rx,vSolid );
                        }
                        if( useExactSolidTraction && knownSolution==Parameters::userDefinedKnownSolution )
                        {
                            printF("\n--INS-- AMP-IMP-BC: **TEST** set exact KNOWN-SOLUTION values for solid traction ***TEMP***, t=%9.3e\n\n",t);
                            int body=0;
                            parameters.getUserDefinedDeformingBodyKnownSolution( body,Parameters::boundaryTraction,
                                                                                                                                      t, grid, mg, Ib1,Ib2,Ib3,Rx,solidTraction );
                        }


            // evaluate the right-hand-side to the AMP velocity BCs
                        {
                            if( debug() & 4 || t<=3.*dt ) 
                            {
                                fPrintF(debugFile,"--INS-- CHECK-AMP-BC-IMP FOR BULK SOLID MODEL, t=%9.3e dt=%9.3e zf=%9.3e zp=%9.3e zs=%9.3e "
                                              "alpha=%9.2e **new**\n",t,dt,zf,zp,zs,alpha);
                                printF("\n--INS-- CHECK-AMP-BC-IMP FOR BULK SOLID MODEL, t=%9.3e dt=%9.3e zf=%9.3e zp=%9.3e zs=%9.3e "
                                              "alpha=%9.2e **new**\n",t,dt,zf,zp,zs,alpha);
                            }
                            if( FALSE )
                            {
                                uOld.periodicUpdate(); // **TEST***
                            }
                            const real beta = (mu*theta*dt/rho);
              // -- evaluate derivatives that appear in the INS equations ---       
                            MappedGridOperators & opOld = *(uOld.getOperators()); 
                            realSerialArray uOldxx(Ib1,Ib2,Ib3,V), uOldyy(Ib1,Ib2,Ib3,V);
                            realSerialArray pOldx(Ib1,Ib2,Ib3), pOldy(Ib1,Ib2,Ib3), pOldz;
                            opOld.derivative(MappedGridOperators::xxDerivative ,uOldLocal,uOldxx,Ib1,Ib2,Ib3,V);
                            opOld.derivative(MappedGridOperators::yyDerivative ,uOldLocal,uOldyy,Ib1,Ib2,Ib3,V);
                            opOld.derivative(MappedGridOperators::xDerivative ,uOldLocal,pOldx,Ib1,Ib2,Ib3,pc);
                            opOld.derivative(MappedGridOperators::yDerivative ,uOldLocal,pOldy,Ib1,Ib2,Ib3,pc);
                            if( numberOfDimensions==3 )
                            {
                                pOldz.redim(Ib1,Ib2,Ib3);
                                opOld.derivative(MappedGridOperators::zDerivative ,uOldLocal,pOldz,Ib1,Ib2,Ib3,pc);
                            }
              // -- We assume a predicted value for "p" at the new time
                            MappedGridOperators & op = *(u.getOperators()); 
                            realSerialArray px(Ib1,Ib2,Ib3), py(Ib1,Ib2,Ib3),pz;
                            op.derivative(MappedGridOperators::xDerivative ,uLocal,px,Ib1,Ib2,Ib3,pc);
                            op.derivative(MappedGridOperators::yDerivative ,uLocal,py,Ib1,Ib2,Ib3,pc);
                            if( numberOfDimensions==3 )
                            {
                                pz.redim(Ib1,Ib2,Ib3);
                                op.derivative(MappedGridOperators::zDerivative ,uLocal,pz,Ib1,Ib2,Ib3,pc);
                            }
                            realSerialArray ux(Ib1,Ib2,Ib3,V), uy(Ib1,Ib2,Ib3,V),uz;
                            op.derivative(MappedGridOperators::xDerivative ,uLocal,ux,Ib1,Ib2,Ib3,V);
                            op.derivative(MappedGridOperators::yDerivative ,uLocal,uy,Ib1,Ib2,Ib3,V);
                            if( numberOfDimensions==3 )
                            {
                                uz.redim(Ib1,Ib2,Ib3,V);
                                op.derivative(MappedGridOperators::zDerivative ,uLocal,uz,Ib1,Ib2,Ib3,V);
                            }
                            RealArray div(Ib1,Ib2,Ib3);
                            div = ux(Ib1,Ib2,Ib3,uc)+uy(Ib1,Ib2,Ib3,vc);
                            if( numberOfDimensions==3 )
                                div += uz(Ib1,Ib2,Ib3,wc);
                            real maxDiv = max(fabs(div));
              // printF("--INS--check-AMP-IMP-BC: max divergence on boundary =%.2e (dt=%.2e)\n",maxDiv,dt);
                            fPrintF(debugFile,"--INS--check-AMP-IMP-BC: max divergence on boundary =%.2e (dt=%.2e)\n",maxDiv,dt);
                            if( debug() & 4 )
                                ::display(div,"Divergence on the boundary","%.2e ");
                            RealArray delta(3,3); // Dirac delta
                            delta=0.; delta(0,0)=1.; delta(1,1)=1.; delta(2,2)=1.;
                            RealArray solidChar(Ib1,Ib2,Ib3,Rx);   // solid characteristic variable 
                            RealArray fluidRhs(Ib1,Ib2,Ib3,Rx);    // fluid interior update without implicit viscous term
              // Solid characteristic variable used in the RHS 
              // Use MINUS of solid traction since the solid normal is in the opposite direction to the fluid
                            solidChar(Ib1,Ib2,Ib3,Rx) = -solidTraction(Ib1,Ib2,Ib3,Rx) + zs*vSolidLocal(Ib1,Ib2,Ib3,Rx);
              // fluidRhs(Ib1,Ib2,Ib3,dir) = fluid interior update without implicit viscous term 
                            for( int dir=0; dir<numberOfDimensions; dir++ )
                            {
                // pd = px, py or pz 
                                RealArray & pd    = dir==0 ? px    : dir==1? py    : pz;
                                RealArray & pOldd = dir==0 ? pOldx : dir==1? pOldy : pOldz;
                // fluidRhs = v^{n-1} - (theta*dt/rho)*grad(p)^n + ((1-theta)*dt/rho)*( -grad(p)^{n-1} + mu*Delta(v)^{n-1} )
                                fluidRhs(Ib1,Ib2,Ib3,dir)= uOldLocal(Ib1,Ib2,Ib3,uc+dir) 
                                    + (-implicitFactor*dt/rho)*( pd(Ib1,Ib2,Ib3) )
                                    + ((1.-implicitFactor)*dt/rho)*( 
                                                mu*( uOldxx(Ib1,Ib2,Ib3,uc+dir)+uOldyy(Ib1,Ib2,Ib3,uc+dir) ) - pOldd(Ib1,Ib2,Ib3) 
                                                                                                  );
                            }
                            if( twilightZoneFlow )
                            {
                // ------- ADD TWILIGHT ZONE --------
                                OGFunction & e = *(parameters.dbase.get<OGFunction* >("exactSolution"));
                                const bool rectangular=false;
                                OV_GET_SERIAL_ARRAY(real,mg.vertex(),xLocal);
                                if( debug() & 4 )
                                {
                                    RealArray pex(Ib1,Ib2,Ib3), pey(Ib1,Ib2,Ib3);
                                    e.gd( pex,xLocal,mg.numberOfDimensions(),rectangular,0,1,0,0,Ib1,Ib2,Ib3,pc,t); // t=new time 
                                    e.gd( pey,xLocal,mg.numberOfDimensions(),rectangular,0,0,1,0,Ib1,Ib2,Ib3,pc,t); // t=new time 
                                    real maxErrPx = max(fabs(px-pex));
                                    real maxErrPy = max(fabs(py-pey));
                                    fprintf(debugFile,">>>INS IMP-AMP-V-BC: error in predicted (px,py)=(%9.3e,%9.3e) at t=%9.3e\n",maxErrPx,maxErrPy,t);
                                    ::display(fabs(px-pex),"error in px",debugFile,"%8.2e ");
                                    ::display(fabs(py-pey),"error in py",debugFile,"%8.2e ");
                                }
                                realSerialArray ue(Ib1,Ib2,Ib3,V), uex(Ib1,Ib2,Ib3,V), uey(Ib1,Ib2,Ib3,V);
                                RealArray tractione(Ib1,Ib2,Ib3,Rx);
                                e.gd( ue  ,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,0,Ib1,Ib2,Ib3,V,t); // t=new time 
                                e.gd( uex ,xLocal,mg.numberOfDimensions(),rectangular,0,1,0,0,Ib1,Ib2,Ib3,V,t); 
                                e.gd( uey ,xLocal,mg.numberOfDimensions(),rectangular,0,0,1,0,Ib1,Ib2,Ib3,V,t); 
                // tau_ij = mu*( partial_j u_i + partial_i u_j )
                                if( numberOfDimensions==2 )
                                {
                  // ====== Get traction in two dimensions =====
                                    RealArray taue(Ib1,Ib2,Ib3,3); // hold viscous stress tensor tau
                                    taue(Ib1,Ib2,Ib3,0) = (2.*mu)*( uex(Ib1,Ib2,Ib3,uc) );                        // tau_11
                                    taue(Ib1,Ib2,Ib3,1) = (   mu)*( uey(Ib1,Ib2,Ib3,uc) + uex(Ib1,Ib2,Ib3,vc) );  // tau_12 = tau_21
                                    taue(Ib1,Ib2,Ib3,2) = (2.*mu)*( uey(Ib1,Ib2,Ib3,vc) );                        // tau_22
                  //  t1 = tau11*n1 + tau12*n2 
                                    tractione(Ib1,Ib2,Ib3,0) = taue(Ib1,Ib2,Ib3,0)*normal(Ib1,Ib2,Ib3,0) + taue(Ib1,Ib2,Ib3,1)*normal(Ib1,Ib2,Ib3,1);
                  //  t2 = tau21*n1 + tau22*n2 
                                    tractione(Ib1,Ib2,Ib3,1) = taue(Ib1,Ib2,Ib3,1)*normal(Ib1,Ib2,Ib3,0) + taue(Ib1,Ib2,Ib3,2)*normal(Ib1,Ib2,Ib3,1);
                                }
                                else 
                                {
                  // ====== Get traction in three dimensions =====
                                    RealArray uez(Ib1,Ib2,Ib3,V);
                                    e.gd( uez ,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,1,Ib1,Ib2,Ib3,V,t); 
                                    RealArray taue(Ib1,Ib2,Ib3,6);
                                    taue(Ib1,Ib2,Ib3,0) = (2.*mu)*( uex(Ib1,Ib2,Ib3,uc) );                        // tau_11
                                    taue(Ib1,Ib2,Ib3,1) = (   mu)*( uey(Ib1,Ib2,Ib3,uc) + uex(Ib1,Ib2,Ib3,vc) );  // tau_12 = tau_21
                                    taue(Ib1,Ib2,Ib3,2) = (   mu)*( uez(Ib1,Ib2,Ib3,uc) + uex(Ib1,Ib2,Ib3,wc) );  // tau_13 = tau_31
                                    taue(Ib1,Ib2,Ib3,3) = (2.*mu)*( uey(Ib1,Ib2,Ib3,vc) );                        // tau_22
                                    taue(Ib1,Ib2,Ib3,4) = (   mu)*( uez(Ib1,Ib2,Ib3,vc) + uey(Ib1,Ib2,Ib3,wc) );  // tau_23 = tau_32
                                    taue(Ib1,Ib2,Ib3,5) = (2.*mu)*( uez(Ib1,Ib2,Ib3,wc) );                        // tau_33
                  //  t1 = tau11*n1 + tau12*n2 + tau13*n3 
                                    tractione(Ib1,Ib2,Ib3,0) = ( taue(Ib1,Ib2,Ib3,0)*normal(Ib1,Ib2,Ib3,0) + 
                                                                                              taue(Ib1,Ib2,Ib3,1)*normal(Ib1,Ib2,Ib3,1) +
                                                                                              taue(Ib1,Ib2,Ib3,2)*normal(Ib1,Ib2,Ib3,2) );
                  //  t2 = tau21*n1 + tau22*n2 + tau23*n3 
                                    tractione(Ib1,Ib2,Ib3,1) = ( taue(Ib1,Ib2,Ib3,1)*normal(Ib1,Ib2,Ib3,0) + 
                                                                                              taue(Ib1,Ib2,Ib3,3)*normal(Ib1,Ib2,Ib3,1) +
                                                                                              taue(Ib1,Ib2,Ib3,4)*normal(Ib1,Ib2,Ib3,2) );
                  //  t3 = tau31*n1 + tau32*n2 + tau33*n3 
                                    tractione(Ib1,Ib2,Ib3,2) = ( taue(Ib1,Ib2,Ib3,2)*normal(Ib1,Ib2,Ib3,0) + 
                                                                                              taue(Ib1,Ib2,Ib3,4)*normal(Ib1,Ib2,Ib3,1) +
                                                                                              taue(Ib1,Ib2,Ib3,5)*normal(Ib1,Ib2,Ib3,2) );
                                }
                                solidChar(Ib1,Ib2,Ib3,Rx) = tractione(Ib1,Ib2,Ib3,Rx) + zs*ue(Ib1,Ib2,Ib3,V);
                // fluidRhs = v^{n-1} + (theta*dt/rho)*grad(p)^n + ((1-theta)*dt/rho)*( -grad(p)^{n-1} + mu*Delta(v)^{n-1} )
                //          = [ I - (mu*theta*dt/rho)*Delta ]* v^n 
                                RealArray & uexx = uex; RealArray & ueyy =uey;  // Reuse arrays
                                e.gd( uexx ,xLocal,mg.numberOfDimensions(),rectangular,0,2,0,0,Ib1,Ib2,Ib3,V,t); 
                                e.gd( ueyy ,xLocal,mg.numberOfDimensions(),rectangular,0,0,2,0,Ib1,Ib2,Ib3,V,t); 
                                vSolidLocal(Ib1,Ib2,Ib3,Rx) = ue(Ib1,Ib2,Ib3,V);
                                if( numberOfDimensions==2 )
                                {
                  //  [ I - (mu*theta*dt/rho)*Delta ]* v^n
                                    fPrintF(debugFile,"addedMassImpBC: beta=%9.3e, max(uexx=%8.2e)\n",beta,max(fabs(uexx)));
                                    fluidRhs(Ib1,Ib2,Ib3,Rx)=ue(Ib1,Ib2,Ib3,V) - beta*( uexx+ueyy ); 
                                }
                                else
                                {
                                    RealArray & uezz = tractione; // Reuse arrays
                                    e.gd( uezz ,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,2,Ib1,Ib2,Ib3,V,t); 
                  //  [ I - (mu*theta*dt/rho)*Delta ]* v^n
                  // *** THIS IS WRONG ***
                                    fluidRhs(Ib1,Ib2,Ib3,Rx)=ue(Ib1,Ib2,Ib3,V) - (mu*theta*dt/rho)*( uexx+ueyy+uezz ); 
                                }
                            }
                            RealArray tau(Ib1,Ib2,Ib3,3); // holds viscous stress tensor tau
                            RealArray traction(Ib1,Ib2,Ib3,3);
                            RealArray res(Ib1,Ib2,Ib3);
                            if( numberOfDimensions==2 )  
                            {
                                tau(Ib1,Ib2,Ib3,0) = (2.*mu)*( ux(Ib1,Ib2,Ib3,uc) );                        // tau_11
                                tau(Ib1,Ib2,Ib3,1) = (   mu)*( uy(Ib1,Ib2,Ib3,uc) + ux(Ib1,Ib2,Ib3,vc) );   // tau_12 = tau_21
                                tau(Ib1,Ib2,Ib3,2) = (2.*mu)*( uy(Ib1,Ib2,Ib3,vc) );                        // tau_22
                //  t1 = tau11*n1 + tau12*n2 
                                traction(Ib1,Ib2,Ib3,0) = tau(Ib1,Ib2,Ib3,0)*normal(Ib1,Ib2,Ib3,0) + tau(Ib1,Ib2,Ib3,1)*normal(Ib1,Ib2,Ib3,1);
                //  t2 = tau21*n1 + tau22*n2 
                                traction(Ib1,Ib2,Ib3,1) = tau(Ib1,Ib2,Ib3,1)*normal(Ib1,Ib2,Ib3,0) + tau(Ib1,Ib2,Ib3,2)*normal(Ib1,Ib2,Ib3,1);
                // Check BC:
                //    tv.[  tau.nv + zs*v - (sigmaSolid.nv + zs*vSolid ) ]=0 
                // In 2D the tangent-vector is tv = [ -n_y, n_x ]
                                res = ( -normal(Ib1,Ib2,Ib3,1)*( traction(Ib1,Ib2,Ib3,0) + zs*uLocal(Ib1,Ib2,Ib3,uc) - solidChar(Ib1,Ib2,Ib3,0) )
                                                +normal(Ib1,Ib2,Ib3,0)*( traction(Ib1,Ib2,Ib3,1) + zs*uLocal(Ib1,Ib2,Ib3,vc) - solidChar(Ib1,Ib2,Ib3,1) ) );
                            }
                            else
                            {
                                OV_ABORT("finish me for 3D");
                            }
                            real maxChar = max(fabs(res));
                            fPrintF(debugFile,"--INS--check-AMP-IMP-BC: max error in tv.( tau.nv + zs*v )= tv.solidChar =%.2e\n",maxChar);
                            if( debug() & 8 )
                                ::display(res,"Residual in  tv.( tau.nv + zs*v )= ... on the boundary",debugFile,"%.2e ");
              // ============= Check tangential component of interior eqn on the boundary =============
                            realSerialArray uxx(Ib1,Ib2,Ib3,V), uyy(Ib1,Ib2,Ib3,V), uzz;
                            op.derivative(MappedGridOperators::xxDerivative ,uLocal,uxx,Ib1,Ib2,Ib3,V);
                            op.derivative(MappedGridOperators::yyDerivative ,uLocal,uyy,Ib1,Ib2,Ib3,V);
                            if( numberOfDimensions==3 )
                            {
                                uzz.redim(Ib1,Ib2,Ib3,V);
                                op.derivative(MappedGridOperators::zzDerivative ,uLocal,uzz,Ib1,Ib2,Ib3,V);
                            }
                            RealArray unp1(Ib1,Ib2,Ib3,V);
                            unp1(Ib1,Ib2,Ib3,V) = fluidRhs(Ib1,Ib2,Ib3,Rx) + (implicitFactor*mu*dt/rho)*( uxx + uyy );
                            if( numberOfDimensions==2 )  
                            {
                // In 2D the tangent-vector is tv = [ -n_y, n_x ]
                                res = ( -normal(Ib1,Ib2,Ib3,1)*( unp1(Ib1,Ib2,Ib3,uc)-uLocal(Ib1,Ib2,Ib3,uc) )
                                                +normal(Ib1,Ib2,Ib3,0)*( unp1(Ib1,Ib2,Ib3,vc)-uLocal(Ib1,Ib2,Ib3,vc) ) );
                            }
                            real maxRes = max(fabs(res));
                            fPrintF(debugFile,"--INS--check-AMP-IMP-BC: max error in tv.( Eqn on boundary ) =%.2e\n",maxRes);
                            if( debug() & 8 )
                                ::display(res,"Residual in  tv.( Eqn on boundary )= ... on the boundary",debugFile,"%.2e ");
              // ============= Check normal component AMP condition =============
                            if( numberOfDimensions==2 )  
                            {
                                res = ( normal(Ib1,Ib2,Ib3,0)*( uLocal(Ib1,Ib2,Ib3,uc) - (alpha*unp1(Ib1,Ib2,Ib3,uc)+(1.-alpha)*vSolidLocal(Ib1,Ib2,Ib3,0)))
                                              +normal(Ib1,Ib2,Ib3,1)*( uLocal(Ib1,Ib2,Ib3,vc) - (alpha*unp1(Ib1,Ib2,Ib3,vc)+(1.-alpha)*vSolidLocal(Ib1,Ib2,Ib3,1))) );
                                if( correctionStage>0 )
                                {
                  // jump in normal traction 
                                    RealArray nSigmaSolidN, nSigmaFluidN;
                  // **use minus for solidTraction **
                                    nSigmaSolidN= -(normal(Ib1,Ib2,Ib3,0)*solidTraction(Ib1,Ib2,Ib3,0) +normal(Ib1,Ib2,Ib3,1)*solidTraction(Ib1,Ib2,Ib3,1) );
                  // n.sigmaFluid.n = - p + n.Tau.n 
                                    nSigmaFluidN = -uLocal(Ib1,Ib2,Ib3,pc) +
                                        normal(Ib1,Ib2,Ib3,0)*traction(Ib1,Ib2,Ib3,0) + normal(Ib1,Ib2,Ib3,1)*traction(Ib1,Ib2,Ib3,1);
                  // ** fix me for TZ **
                                    real tractionJump = max(fabs( nSigmaSolidN - nSigmaFluidN));
                                    fPrintF(debugFile,"Correction stage: Max jump in traction = %.2e\n",tractionJump);
                                    res -=  (1./(zf+zp))*( nSigmaSolidN - nSigmaFluidN );
                                }
                            }
                            else
                            {
                                OV_ABORT("finish me...");
                            }
                            maxRes = max(fabs(res));
                            fPrintF(debugFile,"--INS--check-AMP-IMP-BC: max error in nv.( vI - [alpha*vf + (1-alpha)*vs ] ) =%.2e\n",maxRes);
                            if( debug() & 8 )
                                ::display(res,"Residual in  nv.( vI - [alpha*vf + (1-alpha)*vs ] )",debugFile,"%.2e ");
                            if( false )
                            {
                                ::display(px,"px","%.9e ");
                                ::display(py,"py","%.9e ");
                                ::display(solidChar(Ib1,Ib2,Ib3,Rx),"solidChar","%.9e ");
                                ::display(fluidRhs(Ib1,Ib2,Ib3,Rx),"fluidRhs","%.9e ");
                                ::display(uxx(Ib1,Ib2,Ib3,V),"uxx","%.9e ");
                                ::display(uyy(Ib1,Ib2,Ib3,V),"uyy","%.9e ");
                                ::display(unp1(Ib1,Ib2,Ib3,V),"unp1","%.9e ");
                                ::display(uLocal(Ib1,Ib2,Ib3,V),"uLocal","%.9e ");
                                ::display(vSolidLocal(Ib1,Ib2,Ib3,Rx),"vSolid","%.9e ");
                                ::display(normal(Ib1,Ib2,Ib3,Rx),"normalVector","%.6e ");
                                const realArray & rx = mg.inverseVertexDerivative();
                                Range all;
                                ::display(rx(Ib1,Ib2,Ib3,all),"rx(Ib1,Ib2,Ib3,all)","%.6e ");
                                OV_ABORT("stop here for now");
                            }
            //  -- "vector" versions of macros that use full normal vector : 
                        #define AMPV(n,m) (delta(n,m) - normal(Ib1,Ib2,Ib3,n)*normal(Ib1,Ib2,Ib3,m) )
            // Here is the "P" operator, (1-alpha) n n^T 
                        #define AMGPV(n,m) ((1.-alpha)*normal(Ib1,Ib2,Ib3,n)*normal(Ib1,Ib2,Ib3,m))
            // Here is I-P 
                        #define AMGV(n,m) (delta(n,m) - (1.-alpha)*normal(Ib1,Ib2,Ib3,n)*normal(Ib1,Ib2,Ib3,m))
              // for( int n=0; n<numberOfDimensions; n++ )
              // {
              //   // --- RHS FOR GHOST EQUATION: ---
              //   //  RHS = (1/mu)*( I - n n^T)( solidTraction + zs*vs )
              //   uLocal(Ig1,Ig2,Ig3,uc+n) = (1./mu)*( AMPV(n,0)*solidChar(Ib1,Ib2,Ib3,0) + AMPV(n,1)*solidChar(Ib1,Ib2,Ib3,1) );
              //   // --- RHS FOR BOUNDARY EQUATION: ---
              //   if( false )
              //   {
              //     // TEST
              //   }
              //   else
              //   {
              //     // printF("AM-IMP_BC: n=%i,  AMGV=%9.3e %9.3e %9.3e,  AMGPV=%9.3e %9.3e\n",n,
              //     // RHS = (I-P)*( v - beta*Delta(v) ) + P*vs 
              //     uLocal(Ib1,Ib2,Ib3,uc+n) = ( AMGV(n,0) *fluidRhs(Ib1,Ib2,Ib3,0)    + AMGV(n,1) *fluidRhs(Ib1,Ib2,Ib3,1) + 
              //                                  AMGPV(n,0)*vSolidLocal(Ib1,Ib2,Ib3,0) + AMGPV(n,1)*vSolidLocal(Ib1,Ib2,Ib3,1) );
              //   }
              // }
              // // ------ corners with dirichlet sides are treated specially -------
              // int iv[3], &i1=iv[0], &i2=iv[1], &i3=iv[2];
              // i3=gid(0,2); // set default
              // int axisp1= (axis+1) % numberOfDimensions;
              // for( int sidep1=0; sidep1<=1; sidep1++ )
              // {
              //   if( mg.boundaryCondition(sidep1,axisp1)==Parameters::dirichletBoundaryCondition )
              //   {
              //     if( numberOfDimensions==3 )
              //     {
              //       OV_ABORT("AMP IMP BC CORNER -- finish me for 3D");
              //     }
              //     for( int ghost=0; ghost<=1; ghost++ )  // set boundary and ghost values
              //     {
              //       iv[axis  ]=gid(side  ,axis  ) + ghost*(2*side-1);
              //       iv[axisp1]=gid(sidep1,axisp1);
              //       printf("Set AMP -IMP RHS for pt=(%i,%i,%i)\n",i1,i2,i3);
              //       uLocal(i1,i2,i3,V)=0.;
              //       if( twilightZoneFlow )
              //       {
              //         // ------- FIX FOR TWILIGHT ZONE --------
              //         OGFunction & e = *(parameters.dbase.get<OGFunction* >("exactSolution"));
              //         const bool rectangular=false;
              //         Index J1=i1, J2=i2, J3=i3;
              //         OV_GET_SERIAL_ARRAY(real,mg.vertex(),xLocal);
              //         realSerialArray ue(J1,J2,J3,V);
              //         e.gd( ue  ,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,0,J1,J2,J3,V,t); // t=new time 
              //         uLocal(i1,i2,i3,V)=ue;
              //       }
              //     }
              //   }
              // }
                        }

                    }
                    
                }
            }
        }
    }
    

    return 0;
}


