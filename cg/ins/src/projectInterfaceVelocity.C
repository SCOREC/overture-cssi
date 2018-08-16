// This file automatically generated from projectInterfaceVelocity.bC with bpp.
#include "Cgins.h"
#include "Parameters.h"
#include "turbulenceModels.h"
#include "Insbc4WorkSpace.h"
#include "App.h"
#include "ParallelUtility.h"
#include "ParallelGridUtility.h"
#include "DeformingBodyMotion.h"
#include "BeamModel.h"
#include "BeamFluidInterfaceData.h"
#include "Interface.h"

#define ForBoundary(side,axis)   for( axis=0; axis<mg.numberOfDimensions(); axis++ ) for( side=0; side<=1; side++ )

#define FOR_3D(i1,i2,i3,I1,I2,I3) int I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase();  int I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(); for(i3=I3Base; i3<=I3Bound; i3++) for(i2=I2Base; i2<=I2Bound; i2++) for(i1=I1Base; i1<=I1Bound; i1++)

#define FOR_3(i1,i2,i3,I1,I2,I3) I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase();  I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(); for(i3=I3Base; i3<=I3Bound; i3++) for(i2=I2Base; i2<=I2Bound; i2++) for(i1=I1Base; i1<=I1Bound; i1++)

// lapack routines
#ifdef OV_USE_DOUBLE
    #define GECO EXTERN_C_NAME(dgeco)
    #define GESL EXTERN_C_NAME(dgesl)
#else
    #define GECO EXTERN_C_NAME(sgeco)
    #define GESL EXTERN_C_NAME(sgesl)
#endif

#define insExplicitAMPVelocityBC EXTERN_C_NAME(insexplicitampvelocitybc)

extern "C"
{
void insExplicitAMPVelocityBC(const int&bcOption,const int&nd,
                                                            const int&nd1a,const int&nd1b,const int&nd2a,const int&nd2b,
                                                            const int&nd3a,const int&nd3b,const int&nd4a,const int&nd4b,
                                                            const int&ipar,const real&rpar, real&u, real&un, const int&mask, const real&x,const real&rx, 
                                                            const int&nds1a,const int&nds1b,const int&nds2a,const int&nds2b, const int&nds3a,const int&nds3b,
                                                            const real &solidVelocity, const real&solidTraction, const real &fluidVelocity, 
                                                            const int&bc, const int&indexRange, const int&ierr );
}


extern "C"
{
    void GECO( real & b, const int & nbd, const int & nb, int & ipvt, real & rcond, real & work );
    void GESL( real & a, const int & lda,const int & n,const int & ipvt, real & b, const int & job);
}


// These next two pointers are used in evalUserDefinedDeformingBodyKnownSolution 
static Parameters *parametersPointer=NULL;
static MappedGrid *mappedGridPointer=NULL;

#define evalUserDefinedDeformingBodyKnownSolution EXTERN_C_NAME(evaluserdefineddeformingbodyknownsolution)
extern "C"
{
// Function to allow the deforming body known solution to be called from Fortran
// The two global pointers parametersPointer and mappedGridPointer must be set
void evalUserDefinedDeformingBodyKnownSolution( int & body, int & stateChoice, real & t, int & grid, 
                                                                                              int & i1, int & i2, int & i3, int & n, real *val )
{
    assert( parametersPointer!=NULL && mappedGridPointer!=NULL );

    Parameters & parameters = *parametersPointer;
    MappedGrid & mg = *mappedGridPointer;

    Parameters::DeformingBodyStateOptionEnum stateOption = (Parameters::DeformingBodyStateOptionEnum)stateChoice;
    Index Ib1=i1, Ib2=i2, Ib3=i3;
    RealArray value(Ib1,Ib2,Ib3,5);
    parameters.getUserDefinedDeformingBodyKnownSolution( body,stateOption,t, grid, mg, Ib1,Ib2,Ib3,n,value );
    for( int dir=0; dir<mg.numberOfDimensions(); dir++ )
    {
        val[dir]=value(i1,i2,i3,dir);
    }
    
}

}

// ================================================================================================
/// Macro: 
///  Assign the AMP velocity BC's for explicit time-stepping
// ================================================================================================


// ==================================================================================
// Macro: Project the interface velocity for the beam model
// ==================================================================================

// =================================================================================================
/// Macro: Project the interface velocity for the bulk solid case -- linearized version for testing
// =================================================================================================

// =========================================================================================
// Macro: Project the interface velocity for the bulk solid case 
//
// =========================================================================================


// =========================================================================================
// Macro: Project the interface velocity for the bulk solid case 
// Interface conditions that define values on the boundary and one ghost :
//  (1)  tv.tau.n + zs* tv.v = tv.sigmaBar.n + zs tv.vs
//  (2)  div( v ) = 0
//  (3)  nv.v = alpha nv.ve + (1-alpha) nv.vs
//  (4)  tv.v - tv.ve = 0 
//      (v -vn)/dt - (nu/2)*( Delta(v) ) = (3/2)Fn - (1/2)Fnm1 + (nu/2)*Delta(vn) 
//
//      
//      v - (nu*dt/2)*( Delta(v) ) = vn + (3/2)*dt*Fn - (1/2)*dt*Fnm1 + (nu*dt/2)*Delta(vn) 
// =========================================================================================

// =======================================================================================
/// \brief Project the velocity on the interface for FSI problems. 
/// \return value: 1=an interface was assigned, 0=no interface was assigned
///
/// NOTE:
//     To project the fluid interface velocity we set the gridVelocity on the boundary
//   as this will later be used to set the values on the boundary.
// =======================================================================================
int Cgins::
projectInterfaceVelocity(const real & t, realMappedGridFunction & u, 
                                                  realMappedGridFunction & uOld,
                                                  realMappedGridFunction & gridVelocity,
                                                  const int & grid,
                                                  const real & dt0 /* =-1. */)
    
{
    int interfaceWasAssigned=0;
    
    MovingGrids & movingGrids = parameters.dbase.get<MovingGrids >("movingGrids");
    if( movingGrids.getNumberOfDeformingBodies()==0 )
        return interfaceWasAssigned;

    real dt=dt0;
    if( dt<= 0. )
        dt = parameters.dbase.get<real>("dt");  // *wdh* 2017/05/31
  // const real & dt = parameters.dbase.get<real>("dt");
    const real epsT = REAL_EPSILON*10.;
    assert( dt> epsT );
    
    FILE *&debugFile = parameters.dbase.get<FILE* >("debugFile");
    FILE *&pDebugFile = parameters.dbase.get<FILE* >("pDebugFile");
    
    const real nu = parameters.dbase.get<real>("nu");
    const real rho = parameters.dbase.get<real>("rho");
    const real mu = rho*nu;

    const bool & useAddedMassAlgorithm = parameters.dbase.get<bool>("useAddedMassAlgorithm");
    const bool & projectAddedMassVelocity = parameters.dbase.get<bool>("projectAddedMassVelocity");
    const bool & projectNormalComponentOfAddedMassVelocity =
                              parameters.dbase.get<bool>("projectNormalComponentOfAddedMassVelocity");
    const bool & projectVelocityOnBeamEnds = parameters.dbase.get<bool>("projectVelocityOnBeamEnds"); 
    const int & correctionStage = parameters.dbase.get<int>("correctionStage");
    assert( correctionStage>=0 );

    const Parameters::KnownSolutionsEnum & knownSolution = 
        parameters.dbase.get<Parameters::KnownSolutionsEnum >("knownSolution"); 

    assert(  useAddedMassAlgorithm && projectAddedMassVelocity && parameters.gridIsMoving(grid) );
    
    MappedGrid & mg = *u.getMappedGrid();
    const int numberOfDimensions = mg.numberOfDimensions();

    const bool gridIsMoving = parameters.gridIsMoving(grid);

  // These next two pointers are used in evalUserDefinedDeformingBodyKnownSolution 
    parametersPointer=&parameters;
    mappedGridPointer=&mg;


    if( t <= 2.*dt )
    {
        printF("\n"
                      "------------------------------------------------------------------------------------------------------\n"
                      " --INS-- projectInterfaceVelocity: ADDED MASS ALGORITHM - PROJECT VELOCITY at t=%8.2e correction=%i\n"
                      "------------------------------------------------------------------------------------------------------\n"
                        ,t,correctionStage);

        if( debug() & 4 )
            fPrintF(debugFile,"\n"
                            "------------------------------------------------------------------------------------------------------\n"
                            " --INS-- projectInterfaceVelocity: ADDED MASS ALGORITHM - PROJECT VELOCITY at t=%8.2e correction=%i\n"
                            "------------------------------------------------------------------------------------------------------\n"
                            ,t,correctionStage);
    }
    


  // get components of solution
    const int uc = parameters.dbase.get<int >("uc");
    const int vc = parameters.dbase.get<int >("vc");
    const int wc = parameters.dbase.get<int >("wc");
    const int tc = parameters.dbase.get<int >("tc");
    const int pc = parameters.dbase.get<int >("pc");
    const int & nc = parameters.dbase.get<int >("nc");

    const bool twilightZoneFlow = parameters.dbase.get<bool >("twilightZoneFlow");

    const int orderOfAccuracy=min(4,parameters.dbase.get<int >("orderOfAccuracy"));
    Range V = Range(uc,uc+numberOfDimensions-1);

    BoundaryData::BoundaryDataArray & pBoundaryData = parameters.getBoundaryData(grid); // this will create the BDA if it is not there
    std::vector<BoundaryData> & boundaryDataArray =parameters.dbase.get<std::vector<BoundaryData> >("boundaryData");
    BoundaryData & bd = boundaryDataArray[grid];
            
    const Parameters::InterfaceCommunicationModeEnum & interfaceCommunicationMode= 
        parameters.dbase.get<Parameters::InterfaceCommunicationModeEnum>("interfaceCommunicationMode");

  // get grid spacing
    const bool isRectangular = mg.isRectangular();

    mg.update(MappedGrid::THEinverseVertexDerivative);
    const real dr[3]={mg.gridSpacing(0),mg.gridSpacing(1),mg.gridSpacing(2)};
    const realArray & rsxy = mg.inverseVertexDerivative();
#define RSXY(I1,I2,I3,m1,m2) rsxy(I1,I2,I3,(m1)+numberOfDimensions*(m2))
#define RSXYR(I1,I2,I3,m1,m2) (rsxy(I1+1,I2  ,I3,(m1)+numberOfDimensions*(m2))-rsxy(I1-1,I2  ,I3,(m1)+numberOfDimensions*(m2)))*(1./(2.*dr[0]));
#define RSXYS(I1,I2,I3,m1,m2) (rsxy(I1  ,I2+1,I3,(m1)+numberOfDimensions*(m2))-rsxy(I1  ,I2-1,I3,(m1)+numberOfDimensions*(m2)))*(1./(2.*dr[1]));

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
            
                    interfaceWasAssigned=1;
                    
        	  int body=deformingBodyNumber[side][axis];
        	  if( t<=0. )
          	    printF("--INS-- grid=%i, (side,axis)=(%i,%i) belongs to deforming body %i\n",grid,side,axis,body);

        	  DeformingBodyMotion & deform = movingGrids.getDeformingBody(body);

        	  getBoundaryIndex(mg.gridIndexRange(),side,axis,Ib1,Ib2,Ib3);
                    getGhostIndex(mg.gridIndexRange(),side,axis,Ig1,Ig2,Ig3,+1);  // first ghost line 
                    getGhostIndex(mg.gridIndexRange(),side,axis,Ip1,Ip2,Ip3,-1);  // first line in

        	  Range Rx=numberOfDimensions;
        	  realArray vSolid(Ib1,Ib2,Ib3,Rx); // holds velocity of solid on the boundary
                    #ifndef USE_PPP
                    deform.getVelocityBC( t, side,axis,grid, mg, Ib1,Ib2,Ib3, vSolid );
                    #else
                        OV_ABORT("finish me");
                    #endif

        	  OV_GET_SERIAL_ARRAY(real,gridVelocity,gridVelocityLocal);
        	  OV_GET_SERIAL_ARRAY(real,u   ,uLocal);
                    OV_GET_SERIAL_ARRAY(real,uOld,uP    );
        	  OV_GET_SERIAL_ARRAY(real,vSolid,vSolidLocal);

                    if( projectNormalComponentOfAddedMassVelocity )
                        mg.update(MappedGrid::THEvertexBoundaryNormal);

                    OV_GET_VERTEX_BOUNDARY_NORMAL(mg,side,axis,normal);
        	  
                    if( !deform.isBulkSolidModel() && !deform.isBeamModel() )
        	  {

          	    OV_ABORT("projectInterfaceVelocity::ERROR: un-expected deformation type");

        	  }
        	  else if( deform.isBulkSolidModel() )
        	  {

            // *******************************************************************
            // *************** PROJECT VELOCITY BULK SOLID ***********************
            // *******************************************************************

                        RealArray solidTraction;
                        if( interfaceCommunicationMode==Parameters::requestInterfaceDataWhenNeeded )
                        {
                            if( true )
                            {
                                printF("--INS-- PIV: REQUEST interface traction at t=%9.3e\n",t);
                                if( debug() & 4 ) fPrintF(debugFile,"--INS-- PIV: REQUEST interface traction at t=%9.3e\n",t);
                            }
                            
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
                                ::display(solidTraction(Ib1,Ib2,Ib3,1),"--INS-- PIV: Here is the SOLID TRACTION (I1,I2,I3,1)",debugFile,"%6.3f ");
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
            //   fPrintF(debugFile,"--INS-- PIV: rhoSolid=%9.3e cp=%9.3e cs=%9.3e zp=%9.3e (old: zp=%9.3e) zs=%9.3e\n",
            //          rhoSolid,cp,cs,zp,zpOld,zs);
            //   // printF("  fluidAddedMassLengthScale=%9.3e\n",fluidAddedMassLengthScale);
            // }
                        
	    // // const real & fluidDensity = parameters.dbase.get<real >("fluidDensity");
          	    
            // // fluid impedance = rho*H/dt 
	    // assert( dt>0. );
            // const real zf=fluidDensity*fluidAddedMassLengthScale/dt; 

            // *NEW* June 29, 2018
            // const real zf= (fluidDensity*mu)/(zp*dt);

            // const real alpha = zf/(zf+zp);

                        const real theta = 0.5;

                        if( t<=3.*dt )
                            printF("--INS-- PIV PROJECT INTERFACE VELOCITY FOR BULK SOLID MODEL, t=%9.3e dt=%9.3e zf=%9.3e zp=%9.3e fluidAddedMassLengthScale=%.2e alpha=%.2e\n",
                                          t,dt,zf,zp,fluidAddedMassLengthScale,alpha);

          	    if( t<=3.*dt && debug() &4 ) 
                        {
            	      fPrintF(debugFile,"--INS-- PIV PROJECT INTERFACE VELOCITY FOR BULK SOLID MODEL, t=%9.3e dt=%9.3e zf=%9.3e zp=%9.3e "
                                          "alpha=%9.2e **new**\n",t,dt,zf,zp,alpha);
                        }
                        
            // check if uP is a null pointer
                        bool setAsSolidVelocity = false;
                        if (uP.getLength(0) == 0) {
            // if (&uP == NULL) {
                            printF("projectInterfaceVelocity: u at previous time is a null  pointer ***FIX ME IMPLICIT***\n");
                            setAsSolidVelocity = true;
                        }

                        if( false || setAsSolidVelocity ) // set this to true for debugging
                        {
              // ** do this for now****
              // --- set the gridVelocity to the solid velocity
                            printF("--INS--PIV: SETTING FLUID VELOCITY AS THE SOLID VELOCITY ***TEMP***\n");
                            gridVelocityLocal(Ib1,Ib2,Ib3,Rx)= vSolidLocal(Ib1,Ib2,Ib3,Rx);
                            uLocal(Ib1,Ib2,Ib3,Rx) = vSolidLocal(Ib1,Ib2,Ib3,Rx);
                            return 0;
                        }

                        


            // ::display(solidTraction(Ib1,Ib2,Ib3,0),"solid traction - x component");
            // ::display(solidTraction(Ib1,Ib2,Ib3,1),"solid traction - y component");
            // OV_ABORT("test traction now");

            // ************************************************************
            // ************* PROJECT THE INTERFACE VELOCITY ***************
            // ************************************************************

            // Macro: 
                        if(  FALSE ) // For testing -- 
                        {
              // Use linearized conditions and assume that interface normal is [0,-1]
                            printF("***TEMP*** use linearized AMP velocity BC (explicit) bulk solid\n");
                            {
                // *************** START LINEARIZED ***************
                // Use linearized conditions and assume that interface normal is [0,-1]
                                assert( axis == 1 && side == 0);
                // solve interior equation + mixed BC for lines 0 and ghost
                // 
                // Interior: 
                //   v1(t) = v1(t-dt) + theta*(dt/rho)*( -p.x +mu*(v1.xx + v1.yy ) ) + (1-theta)*( -p.x(-dt) +mu*(v1.xy(-dt) + v1.yy(-dt) ) ) +dt*f
                //
                //   v1 - theta*dt*(mu/rho)*( v1.xx + v1.yy ) = G1 = v1(t-dt) + theta*(dt/rho)*( -p.x  ) 
                //                                                    + (1-theta)*(dt/rho)*( -p.x(-dt) +mu*(v1.xx(-dt) + v1.yy(-dt) ) ) +dt*f
                //   zs*v1 - mu*( v1_y )                      = G2 = mu*( v2_x ) + zs*vs1 - sigmas12 
                //
                //   Diagonal term from approximation to v1.xx is on LHS, off-diagonal terms are on the RHS
                //   Mixed equation assumes that n = [0,-1]
                //
                //   a11*v1 + a12*v1g  = G1p
                //   a21*v1 + a22*v1g  = G2p
                // grid is approximately cartesian, for now use dx
                                real dx[3];
                                dx[0] = dr[0]/RSXY(0,0,0,0,0);
                                dx[1] = dr[1]/RSXY(0,0,0,1,1);
                                printF("dr = (%f,%f)\n",dr[0],dr[1]);
                                printF("dx = (%f,%f)\n",dx[0],dx[1]);
                // apply impedance weighted average
                // gridVelocityLocal(Ib1,Ib2,Ib3,Rx)= alpha*uLocal(Ib1,Ib2,Ib3,V) + (1.-alpha)*vSolidLocal(Ib1,Ib2,Ib3,Rx);
                // uLocal(Ib1,Ib2,Ib3,vc)= alpha*uLocal(Ib1,Ib2,Ib3,vc) + (1.-alpha)*vSolidLocal(Ib1,Ib2,Ib3,1);
                // ::display(alpha*uLocal(Ib1,Ib2,Ib3,vc) + (1.-alpha)*vSolidLocal(Ib1,Ib2,Ib3,1),"v2 on boundary (linearized)");
                // solve the interior equation + mixed BC
                                real a11,a12,a21,a22,det;
                                RealArray g1(Ib1,Ib2,Ib3), g2(Ib1,Ib2,Ib3);  // forcings for BC 
                                real a = ((theta*dt*mu)/rho)/SQR(dx[1]);
                                real b = ((theta*dt*mu)/rho)/SQR(dx[0]);
                                a11=1.+2.*(a+b);   a12= -a;
                                a21=zs;            a22= +mu/(2.*dx[1]);
                                det = a11*a22-a21*a12;
                // TODO: add any RHS forcing needed (eg twilight zone)
                                g1(Ib1,Ib2,Ib3)=
                                    ( uP(Ib1,Ib2,Ib3,uc)  
                                        +(theta*dt/rho)*(
                                            -(uLocal(Ib1+1,Ib2,Ib3,pc)-uLocal(Ib1-1,Ib2,Ib3,pc) )*(1./(2.*dx[0]))  // -p.x(t)
                                            +(uLocal(Ib1+1,Ib2,Ib3,uc)+uLocal(Ib1-1,Ib2,Ib3,uc) )*(mu/SQR(dx[0])) // part of mu*v1.xx
                                            +(                         uLocal(Ib1,Ib2+1,Ib3,uc) )*(mu/SQR(dx[1])) // part of mu*v1.yy
                                            )  
                                        +((1.-theta)*dt/rho)*(
                                            -(uP(Ib1+1,Ib2,Ib3,pc)-uP(Ib1-1,Ib2,Ib3,pc))*(1./(2.*dx[0]))  // -p.x(t-dt)
                                            +(uP(Ib1+1,Ib2,Ib3,uc)-2.*uP(Ib1,Ib2,Ib3,uc)+uP(Ib1-1,Ib2,Ib3,uc))*(mu/SQR(dx[0])) // mu*v1.xx(t-dt)
                                            +(uP(Ib1,Ib2+1,Ib3,uc)-2.*uP(Ib1,Ib2,Ib3,uc)+uP(Ib1,Ib2-1,Ib3,uc))*(mu/SQR(dx[1])) // mu*v1.yy(t-dt)
                                            )
                                        );
                                g2(Ib1,Ib2,Ib3)=( +(mu/(2.*dx[0]))*( uLocal(Ib1+1,Ib2,Ib3,vc)- uLocal(Ib1-1,Ib2,Ib3,vc))
                                                                    +(mu/(2.*dx[1]))*( uLocal(Ib1,Ib2+1,Ib3,uc) )  // part of mu*v1.y 
                                                                    -solidTraction(Ib1,Ib2,Ib3,0) + zs*vSolidLocal(Ib1,Ib2,Ib3,0) 
                                    );
                                RealArray UB(Ib1,Ib2,Ib3), UG(Ib1,Ib2,Ib3);  // for checking
                                UB = ( a22/det)*g1 + (-a12/det)*g2;
                                UG = (-a21/det)*g1 + ( a11/det)*g2;
                // ::display(UB,"v1 on boundary (linearized)");
                // ::display(UG,"v1 on ghost    (linearized)");
                // test the residual of the equations
                                RealArray res1(Ib1,Ib2,Ib3), res2(Ib1,Ib2,Ib3);  // forcings for BC 
                                res1 = (UB-uP(Ib1,Ib2,Ib3,uc))
                                    -(theta*dt/rho)*(
                                        -(uLocal(Ib1+1,Ib2,Ib3,pc)-uLocal(Ib1-1,Ib2,Ib3,pc) )*(1./(2.*dx[0])) // -p.x(t)
                                        +(uLocal(Ib1+1,Ib2  ,Ib3,uc)-2.*UB+uLocal(Ib1-1,Ib2  ,Ib3,uc) )*(mu/SQR(dx[0])) // mu*v1.xx
                                        +(uLocal(Ib1  ,Ib2+1,Ib3,uc)-2.*UB+UG )*(mu/SQR(dx[1])) // mu*v1.yy
                                        )  
                                    -((1.-theta)*dt/rho)*(
                                        -(uP(Ib1+1,Ib2,Ib3,pc)-uP(Ib1-1,Ib2,Ib3,pc))*(1./(2.*dx[0]))  // -p.x(t-dt)
                                        +(uP(Ib1+1,Ib2,Ib3,uc)-2.*uP(Ib1,Ib2,Ib3,uc)+uP(Ib1-1,Ib2,Ib3,uc))*(mu/SQR(dx[0])) // mu*v1.xx(t-dt)
                                        +(uP(Ib1,Ib2+1,Ib3,uc)-2.*uP(Ib1,Ib2,Ib3,uc)+uP(Ib1,Ib2-1,Ib3,uc))*(mu/SQR(dx[1])) // mu*v1.yy(t-dt)
                                        );
                                res2 = zs*UB
                                    - mu*(uLocal(Ib1  ,Ib2+1,Ib3,uc)-UG)*(1./(2.*dx[1]))
                                    - mu*(uLocal(Ib1+1,Ib2  ,Ib3,uc)-uLocal(Ib1-1,Ib2  ,Ib3,uc))*(1./(2.*dx[0]))
                                    -(zs*vSolidLocal(Ib1,Ib2,Ib3,0)-solidTraction(Ib1,Ib2,Ib3,0));
                // ::display(res1,sPrintF("FluidBC: residual for equation 1, t=%9.3e",t),"%8.2e ");
                // ::display(res2,sPrintF("FluidBC: residual for equation 2, t=%9.3e",t),"%8.2e ");
                // ::display(uLocal,sPrintF("FluidBC: uLocal before mixed robin update, t=%9.3e",t),"%8.2e ");
                // uLocal(Ib1,Ib2,Ib3,uc) = ( a22/det)*g1 + (-a12/det)*g2;
                // uLocal(Ig1,Ig2,Ig3,uc) = (-a21/det)*g1 + ( a11/det)*g2;
                // ::display(uLocal,sPrintF("FluidBC: uLocal after mixed robin update, t=%9.3e",t),"%8.2e ");
                // apply div.u = 0
                // uLocal(Ib1,Ib2-1,Ib3,vc) = uLocal(Ib1,Ib2+1,Ib3,vc) + (dx[1]/dx[0])*(uLocal(Ib1+1,Ib2,Ib3,uc)-uLocal(Ib1-1,Ib2,Ib3,uc));
                // gridVelocityLocal(Ib1,Ib2,Ib3,0)= uLocal(Ib1,Ib2,Ib3,uc);
                            }
                            
                        } 
                        else
                        {
                            if( true )
                            {
                // *wdh* newer version
                                {
                                    assert( numberOfDimensions==2 );
                                    bool useHeavySolidLimit=false;
                                    bool useLightSolidLimit=true;
                                    const bool & predictedBoundaryPressureNeeded = parameters.dbase.get<bool>("predictedBoundaryPressureNeeded");
                                    if( !predictedBoundaryPressureNeeded)
                                    {
                                        printF("--INS-- ERROR: Explicit AMP velocity BCs need a predicted pressure\n"
                                                      " You should turn on predictedBoundaryPressureNeeded\n");
                                        OV_ABORT("ERROR");
                                    }
                                    if( useHeavySolidLimit )
                                    {
                    // heavy solid limit 
                                        printF("\n--INS-- projectVelocityBulkSolidNew HEAVY SOLID LIMIT **TEMP** t=%.2e  zp*dt=%.2e\n\n",t,zp*dt);
                                        uLocal(Ib1,Ib2,Ib3,V)=vSolidLocal(Ib1,Ib2,Ib3,Rx);
                                        if( true && knownSolution==Parameters::userDefinedKnownSolution )
                                        {
                                            int body=0;
                                            RealArray vExact(Ib1,Ib2,Ib3,V);
                                            parameters.getUserDefinedDeformingBodyKnownSolution( body,Parameters::boundaryVelocity,
                                                                                                                                                      t, grid, mg, Ib1,Ib2,Ib3,V,vExact );
                                            real diff = max(fabs(uLocal(Ib1,Ib2,Ib3,V)-vExact));
                                            printF("--INS-- projectVelocityBulkSolidNew: max diff |vSolid-vExact|=%.2e, max(vExact)=%.2e t=%.2e\n",
                                                          diff,max(fabs(vExact)),t);
                                            if( false )
                                            {
                                                printF("--INS-- PIV: **TEST** set exact KNOWN-SOLUTION values for velocity on the interface+ghost ***TEMP***, t=%9.3e\n",t);
                                                fPrintF(debugFile,"--INS-- PIV: **TEST** set exact KNOWN-SOLUTION values for velocity on the interface+ghost ***TEMP***, t=%9.3e\n",t);
                                                uLocal(Ib1,Ib2,Ib3,V)=vExact;
                                            }
                      // parameters.getUserDefinedDeformingBodyKnownSolution( body,Parameters::boundaryVelocity,
                      //                                                      t, grid, mg, Ig1,Ig2,Ig3,V,uLocal );
                                        }
                                        BoundaryConditionParameters extrapParams;
                                        extrapParams.orderOfExtrapolation = 3;
                                        u.applyBoundaryCondition(V,BCTypes::extrapolate,Parameters::noSlipWall,0.,t,extrapParams);
                                        u.applyBoundaryCondition(V,BCTypes::generalizedDivergence,Parameters::noSlipWall,0.,t);
                    // -- set the grid velocity on the boundary too
                                        gridVelocityLocal(Ib1,Ib2,Ib3,Rx)= uLocal(Ib1,Ib2,Ib3,V);
                                    }
                                    if( useLightSolidLimit )
                                    {
                    // **THIS IS NOW THE GENERAL CASE WITH "PROPER" AMP BC's ***
                                        printF("\n--INS-- projectVelocityBulkSolidNew LIGHT SOLID LIMIT **TEMP** t=%.2e  zp*dt=%.2e\n\n",t,zp*dt);
                    // we should NOT extrapolate ghost on the corrector step 
                    // Now we extrapolate in time in the predictor
                                        if( false && correctionStage==0 )
                                        {
                      // Extrpolate ghost to give predicted values
                                            BoundaryConditionParameters extrapParams;
                                            extrapParams.orderOfExtrapolation = 3;
                      // u.applyBoundaryCondition(V,BCTypes::extrapolate,Parameters::noSlipWall,0.,t,extrapParams);
                      // Extrapolate all boundaries -- we need to evaluate Delta(v) at corners too
                                            u.applyBoundaryCondition(V,BCTypes::extrapolate,BCTypes::allBoundaries,0.,t,extrapParams);
                                        }
                    // Macro: 
                                        RealArray fluidVelocity(vSolidLocal.dimension(0),vSolidLocal.dimension(1),vSolidLocal.dimension(2),Rx); // make a copy of the predicted fluid values on the boundary 
                                      fluidVelocity(Ib1,Ib2,Ib3,Rx)=uLocal(Ib1,Ib2,Ib3,V);
                                      if( FALSE )
                                      {
                     // **TEMP** FIX ME -- initial projection
                                          uLocal(Ib1,Ib2,Ib3,V)=alpha*uLocal(Ib1,Ib2,Ib3,V)+(1.-alpha)*vSolidLocal(Ib1,Ib2,Ib3,Rx);
                                          u.periodicUpdate();
                                      }
                                      int numIts=1; // 3 
                                      for( int it=0; it<numIts; it++ )
                                      {
                                          {
                                              if( debug() & 2 )
                                                  fPrintF(debugFile,"++ assign explicit AMP Velocity BCs at t=%.3e\n",t);
                      // twilight always needs the vertex: 
                                              bool vertexNeeded = !isRectangular || parameters.dbase.get<bool >("twilightZoneFlow");
                                              OV_GET_SERIAL_ARRAY(real,u,uLocal);
                                              real *pu = uLocal.getDataPointer();
                                              OV_GET_SERIAL_ARRAY(real,uOld,uOldLocal);
                                              real *pun = uOldLocal.getDataPointer();
                                              OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
                                              const int *pmask = maskLocal.getDataPointer();
                                              real temp;
                                              real *pxy=&temp, *prsxy=&temp;
                                              if( !isRectangular )
                                              {
                                                  OV_GET_SERIAL_ARRAY(real,mg.inverseVertexDerivative(),rxLocal);
                                                  prsxy=rxLocal.getDataPointer();
                                              }
                                              if( vertexNeeded )
                                              {
                                                  OV_GET_SERIAL_ARRAY(real,mg.vertex(),xLocal);
                                                  pxy=xLocal.getDataPointer();
                                              }
                       // Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
                       // getIndex(mg.gridIndexRange(),I1,I2,I3);
                       // int includeGhost=1;
                       // bool ok = ParallelUtility::getLocalArrayBounds(u0,uLocal,I1,I2,I3,includeGhost);
                                              OV_GET_SERIAL_ARRAY_CONDITIONAL(real,gridVelocity,gridVelocityLocal,parameters.gridIsMoving(grid));
                                              realSerialArray gtt;  // holds the boundary acceleration
                       // if( parameters.gridIsMoving(grid) )
                       // {
                       //   // Moving grid problem: compute the acceleration on the boundary
                       //   MovingGrids & movingGrids = parameters.dbase.get<MovingGrids >("movingGrids");
                       //   int boundaryAccelerationOption=4; // return gtt
                       //   if( ok )
                       //   {
                       //     // note: insbc4 assumes gtt has the same first 3 dimensions as uLocal
                       //     gtt.redim(uLocal.dimension(0),uLocal.dimension(1),uLocal.dimension(2),numberOfDimensions);
                       //     gtt=0.;
                       //     movingGrids.getBoundaryAcceleration( mg, gtt, grid, t, boundaryAccelerationOption );
                       //   }
                       // }
                                              real *pSolidVelocity = vSolidLocal.getDataPointer();   // pointer to the solid velocity
                                              real *pSolidTraction = solidTraction.getDataPointer();   // pointer to the solid traction
                                              real *pFluidVelocity = fluidVelocity.getDataPointer();   // pointer to the solid traction
                       // check -- is it ok to use gid instead of ir?
                                              IntegerArray indexRangeLocal(2,3), dimLocal(2,3), bcLocal(2,3);
                                              ParallelGridUtility::getLocalIndexBoundsAndBoundaryConditions( u,indexRangeLocal,dimLocal,bcLocal ); 
                       // *wdh* 110311 - Add Boussinesq terms and boundary conditions for T
                                              const InsParameters::PDEModel & pdeModel = parameters.dbase.get<InsParameters::PDEModel >("pdeModel");
                                              const bool assignTemperature = pdeModel==InsParameters::BoussinesqModel ||
                                                                                                            pdeModel==InsParameters::viscoPlasticModel;
                                              const int np = max(Communication_Manager::Number_Of_Processors,1);
                                              int useWhereMask=false; // **NOTE** for  moving grids we may need to evaluate at more points than just mask >0 
                                              real ajs=1.;  // is this used?
                       //   real ajs=getSignForJacobian(mg);
                                              real thermalExpansivity=1.;
                       // parameters.dbase.get<ListOfShowFileParameters >("pdeParameters").getParameter("thermalExpansivity",thermalExpansivity)
                                              real dx[3]={1.,1.,1.};
                                              bool isRectangular=mg.isRectangular();
                                              if( isRectangular )
                                              {
                                                  mg.getDeltaX(dx);
                                              }
                                              int gridType = isRectangular ? 0 : 1;
                                              int iparam[] ={parameters.dbase.get<int >("pc"),
                                                       		 parameters.dbase.get<int >("uc"),
                                                       		 parameters.dbase.get<int >("vc"),
                                                       		 parameters.dbase.get<int >("wc"),
                                                       		 parameters.dbase.get<int >("sc"),
                                                       		 grid,
                                                       		 gridType,
                                                       		 min(4,parameters.dbase.get<int >("orderOfAccuracy")),
                                                       		 (int)parameters.gridIsMoving(grid),
                                                       		 useWhereMask,
                                                       		 parameters.getGridIsImplicit(grid),
                                                       		 (int)parameters.dbase.get<Parameters::ImplicitMethod >("implicitMethod"),
                                                       		 (int)parameters.dbase.get<Parameters::ImplicitOption >("implicitOption"),
                                                       		 (int)parameters.isAxisymmetric(),
                                                       		 (int)parameters.dbase.get<bool >("twilightZoneFlow"),
                                                       		 np,
                                                                            parameters.dbase.get<int >("debug"),
                                                       		 parameters.dbase.get<int >("myid"),
                                                                            (int)assignTemperature,
                                                                            parameters.dbase.get<int >("tc"),
                                                                            side,
                                                                            axis,
                                                                            knownSolution,
                                                                            parameters.dbase.get<int>("addedMassVelocityBC"), // *wdh* June 24, 2018
                                                                            correctionStage
                                              };
                                              real gravity[3];
                                              real rparam[]={dx[0],dx[1],dx[2],
                                                                            mg.gridSpacing(0),
                                                                            mg.gridSpacing(1),
                                                                            mg.gridSpacing(2),
                                                                            parameters.dbase.get<real >("nu"),
                                                                            t,
                                                                            zp, 
                                                                            zs, 
                                                                            alpha, // added-mass ratio in velocity projection
                                                                            mu,
                                                                            ajs,
                                                                            gravity[0],
                                                                            gravity[1],
                                                                            gravity[2],
                                                                            thermalExpansivity,
                                                                            (real &)(parameters.dbase.get<OGFunction* >("exactSolution")), // pointer to TZ
                                                                            dt,
                                                                            REAL_MIN,
                                                                            fluidDensity
                                              };
                                              int ierr=0;
                                              int bcOption=0;
                                              insExplicitAMPVelocityBC(bcOption,
                                                                                                mg.numberOfDimensions(),
                                                                                                uLocal.getBase(0),uLocal.getBound(0),uLocal.getBase(1),uLocal.getBound(1),
                                                                                                uLocal.getBase(2),uLocal.getBound(2),uLocal.getBase(3),uLocal.getBound(3),
                                                                                                iparam[0],rparam[0], *pu, *pun, *pmask, *pxy, *prsxy, 
                                                                                                solidTraction.getBase(0),solidTraction.getBound(0),
                                                                                                solidTraction.getBase(1),solidTraction.getBound(1),
                                                                                                solidTraction.getBase(2),solidTraction.getBound(2),
                                                                                                *pSolidVelocity, *pSolidTraction, *pFluidVelocity,
                                                                                                bcLocal(0,0), indexRangeLocal(0,0), ierr ) ;
                                          }
                                          u.periodicUpdate();
                     // **TEST**
                     // u.applyBoundaryCondition(V,BCTypes::generalizedDivergence,Parameters::noSlipWall,0.,t);
                                      }
                    // -- set the grid velocity on the boundary too
                                        gridVelocityLocal(Ib1,Ib2,Ib3,Rx)= uLocal(Ib1,Ib2,Ib3,V);
                                    }
                  // OV_GET_SERIAL_ARRAY(real,cg[grid].vertex(),xLocal);
                  // assert( u.getOperators()!=NULL );
                  // MappedGridOperators & op = *(u.getOperators());
                  // op.derivative(MappedGridOperators::xDerivative ,nuT,nuTx,I1,I2,I3,0);
                  // // Set the normal component from (3)
                  // RealArray nDotv(Ib1,Ib2,Ib3);
                  // // Set n.v by subtracting current n.v and adding desired value:
                  // // nDotV = nv^T( alpha*( ve ) + (1-alpha)*vSolid ) - nv^T ve
                  // nDotV = ( (alpha-1.)*( normal(Ib1,Ib2,Ib3,0)*u(Ib1,Ib2,Ib3,uc)     + normal(Ib1,Ib2,Ib3,1)*u(Ib1,Ib2,Ib3,vc) ) +
                  //           (1.-alpha)*( normal(Ib1,Ib2,Ib3,0)*vSolid(Ib1,Ib2,Ib3,0) + normal(Ib1,Ib2,Ib3,1)*vSolid(Ib1,Ib2,Ib3,1) ));
                  // u(Ib1,Ib2,Ib3,uc) += normal(Ib1,Ib2,Ib3,0)*nDotV;
                  // u(Ib1,Ib2,Ib3,vc) += normal(Ib1,Ib2,Ib3,1)*nDotV;
                  // // Set divergence 
                  // u.applyBoundaryCondition(V,generalizedDivergence,dirichletBoundaryCondition,0.,t);
                  // Set t.v(0) and t.v(-1) from (1) and (4)
                  //  [ a11  a12 ][ t.v(0)  ] = [ b1 ]
                  //  [ a21  a22 ][ t.v(-1) ]   [ b2 ]
                                }
                            }
                            else
                            {
                                {
                                    assert( axis == 1 && side == 0);
                                    assert( numberOfDimensions == 2 );
                  // solve the system:
                  // n (div.u) + 1/mu (I-n.n') (tau n + zs v) = 1/mu (I-n.n') ( Ts + zs vs )
                  // v - (I + (alpha - 1) n.n') theta dt nu Delta v
                  //      = (I + (alpha - 1) n.n') (vp + theta dt / rho (- grad p)
                  //                                +(1-theta) dt / rho (-grad pp + mu Delta vp))
                  //        (1 - alpha) n.n' vs
                  //
                  // *** See Dropbox/CG/fibr/codes/curvilinearBC/* for maple file ***
                                    const int numberOfEquations=4;
                                    Range R(0,numberOfEquations-1);
                  // solve: A.x = b
                                    RealArray A(R,R), b(R), x(R);
                                    intArray ipvt(R);
                                    RealArray work(R);
                                    real rcond;
                                    int job=0;
                                    real rx,ry,sx,sy, rxr,ryr,sxr,syr, rxs,rys,sxs,sys;
                                    real n1,n2;
                  // predicted values for derivatives
                                    real v1r, v1rr, v1s, v1ss, v1rs; 
                                    real v2r, v2rr, v2s, v2ss, v2rs; 
                  // pressure gradient at current and old time level
                                    real pr, ps, ppr, pps; 
                                    real px, py, ppx, ppy; 
                  // second derivatives at old time level
                                    real vp1r, vp1rr, vp1s, vp1ss, vp1rs; 
                                    real vp2r, vp2rr, vp2s, vp2ss, vp2rs; 
                                    real vp1xx, vp1yy, vp2xx, vp2yy;
                  // solid values for velocity and traction
                                    real vs1, vs2, Ts1, Ts2;
                  // boundary ghost and interior values
                                    real v1B, v1G, v1I, v2B, v2G, v2I, vp1, vp2; 
                  // RHS of equations
                                    real g1,g2,g3,g4;
                                    real s = (side==0) ? 1.0 : -1.0;
                  // arrays for the boundary and ghost lines
                                    Range rAssign(0,1);
                                    RealArray UB(Ib1,Ib2,Ib3,rAssign), UG(Ib1,Ib2,Ib3,rAssign);  
                                    int i1,i2,i3;
                                    FOR_3D(i1,i2,i3,Ib1,Ib2,Ib3) {
                    // print
                                        printF("***projectVelocityBulkSolid (i1,i2,i3)=(%d,%d,%d) zs=%.6e alpha=%.6e, mu=%.3e \n",i1,i2,i3,zs,alpha,mu);
                    // get metric terms
                                        rx  = RSXY (i1,i2,i3,0,0); ry  = RSXY (i1,i2,i3,0,1); sx  = RSXY (i1,i2,i3,1,0); sy  = RSXY (i1,i2,i3,1,1);
                                        rxr = RSXYR(i1,i2,i3,0,0); ryr = RSXYR(i1,i2,i3,0,1); sxr = RSXYR(i1,i2,i3,1,0); syr = RSXYR(i1,i2,i3,1,1);
                                        rxs = RSXYS(i1,i2,i3,0,0); rys = RSXYS(i1,i2,i3,0,1); sxs = RSXYS(i1,i2,i3,1,0); sys = RSXYS(i1,i2,i3,1,1);
                    // get normal
                                        n1 = normal(i1,i2,i3,0); n2 = normal(i1,i2,i3,1);
                                        printF("  normal: (%f,%f)\n",n1,n2);
                    // printF("  (rx ,ry ,sx ,sy ) = (%f,%f,%f,%f)\n",rx ,ry ,sx ,sy );
                    // printF("  (rxr,ryr,sxr,syr) = (%f,%f,%f,%f)\n",rxr,ryr,sxr,syr);
                    // printF("  (rxs,rys,sxs,sys) = (%f,%f,%f,%f)\n",rxs,rys,sxs,sys);
                    // printF("  dr=%e, ds=%e, mu=%e, zs=%e\n\n",dr[0],dr[1],mu,zs);
                    // calculate predicted derivatives
                                        v1r  = (uLocal(i1+1,i2,i3,uc)-uLocal(i1-1,i2,i3,uc))*(1./(2.*dr[0]));
                                        v2r  = (uLocal(i1+1,i2,i3,vc)-uLocal(i1-1,i2,i3,vc))*(1./(2.*dr[0]));
                                        v1rr = (uLocal(i1+1,i2,i3,uc)-2.*uLocal(i1,i2,i3,uc)+uLocal(i1-1,i2,i3,uc))*(1./(SQR(dr[0])));
                                        v2rr = (uLocal(i1+1,i2,i3,vc)-2.*uLocal(i1,i2,i3,vc)+uLocal(i1-1,i2,i3,vc))*(1./(SQR(dr[0])));
                                        v1rs  = (uLocal(i1+1,i2+1,i3,uc)-uLocal(i1-1,i2+1,i3,uc)-uLocal(i1+1,i2-1,i3,uc)+uLocal(i1-1,i2-1,i3,uc))*(1./(4.*dr[0]*dr[1]));
                                        v2rs  = (uLocal(i1+1,i2+1,i3,vc)-uLocal(i1-1,i2+1,i3,vc)-uLocal(i1+1,i2-1,i3,vc)+uLocal(i1-1,i2-1,i3,vc))*(1./(4.*dr[0]*dr[1]));
                    // previous time level
                                        vp1 = uP(i1,i2,i3,uc); vp2 = uP(i1,i2,i3,vc);
                                        vp1r  = (uP(i1+1,i2,i3,uc)-uP(i1-1,i2,i3,uc))*(1./(2.*dr[0]));
                                        vp2r  = (uP(i1+1,i2,i3,vc)-uP(i1-1,i2,i3,vc))*(1./(2.*dr[0]));
                                        vp1rr = (uP(i1+1,i2,i3,uc)-2.*uP(i1,i2,i3,uc)+uP(i1-1,i2,i3,uc))*(1./(SQR(dr[0])));
                                        vp2rr = (uP(i1+1,i2,i3,vc)-2.*uP(i1,i2,i3,vc)+uP(i1-1,i2,i3,vc))*(1./(SQR(dr[0])));
                                        vp1s  = (uP(i1,i2+1,i3,uc)-uP(i1,i2-1,i3,uc))*(1./(2.*dr[1]));
                                        vp2s  = (uP(i1,i2+1,i3,vc)-uP(i1,i2-1,i3,vc))*(1./(2.*dr[1]));
                                        vp1ss = (uP(i1,i2+1,i3,uc)-2.*uP(i1,i2,i3,uc)+uP(i1,i2-1,i3,uc))*(1./(SQR(dr[1])));
                                        vp2ss = (uP(i1,i2+1,i3,vc)-2.*uP(i1,i2,i3,vc)+uP(i1,i2-1,i3,vc))*(1./(SQR(dr[1])));
                                        vp1rs  = (uP(i1+1,i2+1,i3,uc)-uP(i1-1,i2+1,i3,uc)-uP(i1+1,i2-1,i3,uc)+uP(i1-1,i2-1,i3,uc))*(1./(4.*dr[0]*dr[1]));
                                        vp2rs  = (uP(i1+1,i2+1,i3,vc)-uP(i1-1,i2+1,i3,vc)-uP(i1+1,i2-1,i3,vc)+uP(i1-1,i2-1,i3,vc))*(1./(4.*dr[0]*dr[1]));
                                        vp1xx = SQR(rx)*vp1rr+rx*rxr*vp1r+2*rx*sx*vp1rs+rx*sxr*vp1s+rxs*sx*vp1r+SQR(sx)*vp1ss+sx*sxs*vp1s;
                                        vp2xx = SQR(rx)*vp2rr+rx*rxr*vp2r+2*rx*sx*vp2rs+rx*sxr*vp2s+rxs*sx*vp2r+SQR(sx)*vp2ss+sx*sxs*vp2s;
                                        vp1yy = SQR(ry)*vp1rr+ry*ryr*vp1r+2*ry*sy*vp1rs+ry*syr*vp1s+rys*sy*vp1r+SQR(sy)*vp1ss+sy*sys*vp1s;
                                        vp2yy = SQR(ry)*vp2rr+ry*ryr*vp2r+2*ry*sy*vp2rs+ry*syr*vp2s+rys*sy*vp2r+SQR(sy)*vp2ss+sy*sys*vp2s;
                    // get solid values
                                        vs1 = vSolidLocal(i1,i2,i3,0);   vs2 = vSolidLocal(i1,i2,i3,1);
                                        Ts1 = -solidTraction(i1,i2,i3,0); Ts2 = -solidTraction(i1,i2,i3,1);
                    // pressure gradient terms
                                        pr  = (uLocal(i1+1,i2  ,i3,pc)-uLocal(i1-1,i2  ,i3,pc))*(1./(2.*dr[0]));
                                        ps  = (uLocal(i1  ,i2+1,i3,pc)-uLocal(i1  ,i2-1,i3,pc))*(1./(2.*dr[1]));
                                        ppr  = (uP(i1+1,i2  ,i3,pc)-uP(i1-1,i2  ,i3,pc))*(1./(2.*dr[0]));
                                        pps  = (uP(i1  ,i2+1,i3,pc)-uP(i1  ,i2-1,i3,pc))*(1./(2.*dr[1]));
                                        px = rx*pr+sx*ps;
                                        py = ry*pr+sy*ps;
                                        ppx = rx*ppr+sx*pps;
                                        ppy = ry*ppr+sy*pps;
                    // interior points
                                        v1I = uLocal(i1,i2+1,i3,uc); v2I = uLocal(i1,i2+1,i3,vc);
                                        printF(" px,py=%.5e,%5e, old: ppx,ppy=%.5e,%.5e\n",px,py,ppx,ppy);
                                        printF(" vSolid=%.5e,%5e, solidTraction=%.5e,%.5e\n",vSolidLocal(i1,i2,i3,0),vSolidLocal(i1,i2,i3,1),
                                                                solidTraction(i1,i2,i3,0),solidTraction(i1,i2,i3,1));
                                        printF("Old: uxx,uyy=%.5e,%5e, vxx,vyy=%.5e,%.5e\n",vp1xx,vp1yy,vp2xx,vp2yy);
                    // printF("v1I=%e, v2I=%e, v1r=%e, v2r=%e\n\n",v1I,v2I,v1r,v2r);
                    // form A and b
                    // x = [v1B, v2B, v1G, v2G]
                    // *** Equation 1 ***
                                        A(0,0) = (1.0/mu)*(1.0-SQR(n1))*zs; 
                                        A(0,1) = -n1*n2*zs*(1.0/mu); 
                                        A(0,2) = -.5*sx*s*n1/dr[1]
                                            -s*(1.0-SQR(n1))*(sx*n1/dr[1]+0.5*sy*n2/dr[1])
                                            +0.5*SQR(n1)*n2*sy*s/dr[1];
                                        A(0,3) = -.5*sy*s*n1/dr[1]
                                            -0.5*(1-SQR(n1))*sx*s*n2/dr[1]
                                            +s*n1*n2*(0.5*sx*n1/dr[1]+sy*n2/dr[1]);
                                        g1 = (1.-SQR(n1))*(vs1*zs+Ts1)/mu-n1*n2*(vs2*zs+Ts2)/mu; // RHS of equation
                                        b(0) = -((rx*v1r+.5*sx*s*v1I/dr[1]+ry*v2r+.5*sy*s*v2I/dr[1])*n1 // terms from div.u
                                                          +(1.0-SQR(n1))*(2.*(rx*v1r+.5*sx*s*v1I/dr[1])*n1
                                                                                          +(ry*v1r+.5*sy*s*v1I/dr[1]+rx*v2r+.5*sx*s*v2I/dr[1])*n2)
                                                          -n1*n2*        (2.*(ry*v2r+.5*sy*s*v2I/dr[1])*n2
                                                                                          +(ry*v1r+.5*sy*s*v1I/dr[1]+rx*v2r+.5*sx*s*v2I/dr[1])*n1)) // terms from fluid traction
                                            +g1;
                    // printF("b0=%e, b1=%e",-10.*v1I-1.0*v2r,-v1r-10.*v2I);
                    // *** Equation 2 ***
                                        A(1,0) = -n1*n2*zs*(1.0/mu); 
                                        A(1,1) = (1.0/mu)*(1.0-SQR(n2))*zs; 
                                        A(1,2) = -0.5*sx*s*n2/dr[1]
                                            -0.5*s*(1.0-SQR(n2))*sy*n1/dr[1]
                                            +s*n1*n2*(0.5*sy*n2/dr[1]+sx*n1/dr[1]);
                                        A(1,3) = -0.5*sy*s*n2/dr[1]
                                            -s*(1.0-SQR(n2))*(sy*n2/dr[1]+0.5*sx*n1/dr[1])
                                            +0.5*s*SQR(n2)*n1*sx/dr[1]; 
                                        g2 = -n1*n2*(vs1*zs+Ts1)/mu+(1-SQR(n2))*(vs2*zs+Ts2)/mu; // RHS of equation
                                        b(1) = -((rx*v1r+.5*sx*s*v1I/dr[1]+ry*v2r+.5*sy*s*v2I/dr[1])*n2 // terms from div.u
                                                          +(1.0-SQR(n2))*(2.*(ry*v2r+.5*sy*s*v2I/dr[1])*n2
                                                                                          +(ry*v1r+.5*sy*s*v1I/dr[1]+rx*v2r+.5*sx*s*v2I/dr[1])*n1)
                                                          -n1*n2*        (2.*(rx*v1r+.5*sx*s*v1I/dr[1])*n1
                                                                                          +(ry*v1r+.5*sy*s*v1I/dr[1]+rx*v2r+.5*sx*s*v2I/dr[1])*n2)) // terms from fluid traction
                                            +g2;
                    // *** Equation 3 ***
                                        A(2,0) = 1.
                                            -(1.+(alpha-1.)*SQR(n1))
                                            *theta*dt*nu*(-2.*SQR(sx/dr[1])
                                                                        -2.*SQR(sy/dr[1])); 
                                        A(2,1) = -(alpha-1.)*n1*n2
                                            *theta*dt*nu*(-2.*SQR(sx/dr[1])
                                                                        -2.*SQR(sy/dr[1]));
                                        A(2,2) = -(1.+(alpha-1.)*SQR(n1))
                                            *theta*dt*nu*(-.5*rx*sxr*s/dr[1]
                                                                        +SQR(sx/dr[1])
                                                                        -.5*sx*sxs*s/dr[1]
                                                                        -.5*ry*syr*s/dr[1]
                                                                        +SQR(sy/dr[1])
                                                                        -.5*sy*sys*s/dr[1]); 
                                        A(2,3) = -(alpha-1.)*n1*n2
                                            *theta*dt*nu*(-.5*rx*sxr*s/dr[1]
                                                                        +SQR(sx/dr[1])
                                                                        -.5*sx*sxs*s/dr[1]
                                                                        -.5*ry*syr*s/dr[1]
                                                                        +SQR(sy/dr[1])
                                                                        -.5*sy*sys*s/dr[1]);
                                        g3 = (1.+(alpha-1.)*SQR(n1))*(vp1-theta*dt*px/rho+(1.-theta)*dt*(mu*(vp1xx+vp1yy)-ppx)/rho)
                                            +(alpha-1.)*n1*n2*(vp2-theta*dt*py/rho+(1.-theta)*dt*(mu*(vp2xx+vp2yy)-ppy)/rho)
                                            -(alpha-1.)*SQR(n1)*vs1
                                            -(alpha-1.)*n1*n2*vs2; // RHS of equation
                                        b(2) = -(-(1.+(alpha-1.)*SQR(n1))
                                                          *theta*dt*nu*(SQR(rx)*v1rr+rx*rxr*v1r+2.*rx*sx*v1rs+.5*rx*sxr*s*v1I/dr[1]+rxs*sx*v1r
                                                                                      +SQR(sx)*v1I/SQR(dr[1])+.5*sx*sxs*s*v1I/dr[1]+SQR(ry)*v1rr+ry*ryr*v1r+2.*ry*sy*v1rs
                                                                                      +.5*ry*syr*s*v1I/dr[1]+rys*sy*v1r+SQR(sy)*v1I/SQR(dr[1])+.5*sy*sys*s*v1I/dr[1])
                                                          -(alpha-1.)*n1*n2
                                                          *theta*dt*nu*(SQR(rx)*v2rr+rx*rxr*v2r+2.*rx*sx*v2rs+.5*rx*sxr*s*v2I/dr[1]+rxs*sx*v2r
                                                                                      +SQR(sx)*v2I/SQR(dr[1])+.5*sx*sxs*s*v2I/dr[1]+SQR(ry)*v2rr+ry*ryr*v2r+2.*ry*sy*v2rs
                                                                                      +.5*ry*syr*s*v2I/dr[1]+rys*sy*v2r+SQR(sy)*v2I/SQR(dr[1])+.5*sy*sys*s*v2I/dr[1]))
                                            +g3;
                    // *testing* try setting to predicted values
                    // A(2,0) = 1.;
                    // A(2,1) = 0.;
                    // A(2,2) = 0.;
                    // A(2,3) = 0.;
                    // b(2) = uLocal(i1,i2,i3,uc)+(alpha-1)*(SQR(n1)*(uLocal(i1,i2,i3,uc)-vs1)+n1*n2*(uLocal(i1,i2,i3,vc)-vs2));
                    // *** Equation 4 ***
                                        real ds = dr[1];
                                        A(3,0) = -(alpha-1.)*n1*n2
                                            *theta*dt*nu*(-2.*SQR(sx/dr[1])
                                                                        -2.*SQR(sy/dr[1]));
                                        A(3,1) = 1.
                                            -(1.+(alpha-1.)*SQR(n2))*theta*dt*nu*(-2.*SQR(sx/dr[1])
                                                                                                                        -2.*SQR(sy/dr[1]));
                                        A(3,2) = -(alpha-1.)*n1*n2
                                            *theta*dt*nu*(-.5*rx*sxr*s/dr[1]
                                                                        +SQR(sx/dr[1])
                                                                        -.5*sx*sxs*s/dr[1]
                                                                        -.5*ry*syr*s/dr[1]
                                                                        +SQR(sy/dr[1])-.5*sy*sys*s/dr[1]);
                                        A(3,3) = -(1.+(alpha-1.)*SQR(n2))
                                            *theta*dt*nu*(-.5*rx*sxr*s/dr[1]
                                                                        +SQR(sx/dr[1])
                                                                        -.5*sx*sxs*s/dr[1]
                                                                        -.5*ry*syr*s/dr[1]
                                                                        +SQR(sy/dr[1])
                                                                        -.5*sy*sys*s/dr[1]);
                                        g4 = (alpha-1.)*n1*n2*(vp1-theta*dt*px/rho+(1.-theta)*dt*(mu*(vp1xx+vp1yy)-ppx)/rho)
                                            +(1.+(alpha-1.)*SQR(n2))*(vp2-theta*dt*py/rho+(1.-theta)*dt*(mu*(vp2xx+vp2yy)-ppy)/rho)
                                            -(alpha-1.)*n1*n2*vs1
                                            -(alpha-1.)*SQR(n2)*vs2;
                    // g4 = (alpha-1)*n1*n2*(vp1-theta*dt*px/rho+(1-theta)*dt*(mu*(vp1xx+vp1yy)-ppx)/rho)
                    //   +(1+(alpha-1)*SQR(n2))*(vp2-theta*dt*py/rho+(1-theta)*dt*(mu*(vp2xx+vp2yy)-ppy)/rho)
                    //   -(alpha-1)*n1*n2*vs1-(alpha-1)*SQR(n2)*vs2;
                    // printF("%e\n",alpha);
                    // printF("%e,%e,%e,%e\n",(alpha-1)*n1*n2,
                    //        +(1.+(alpha-1.)*SQR(n2)),
                    //        -(alpha-1.)*n1*n2,
                    //        -(alpha-1.)*SQR(n2));
                    // printF("%e,%e\n",(vp2-theta*dt*py/rho+(1.-theta)*dt*(mu*(vp2xx+vp2yy)-ppy)/rho),vs2);
                                        b(3) = -(-(alpha-1.)*n1*n2
                                                          *theta*dt*nu*(SQR(rx)*v1rr+rx*rxr*v1r+2.*rx*sx*v1rs+.5*rx*sxr*s*v1I/dr[1]+rxs*sx*v1r
                                                                                      +SQR(sx)*v1I/SQR(dr[1])+.5*sx*sxs*s*v1I/dr[1]+SQR(ry)*v1rr+ry*ryr*v1r+2.*ry*sy*v1rs
                                                                                      +.5*ry*syr*s*v1I/dr[1]+rys*sy*v1r+SQR(sy)*v1I/SQR(dr[1])+.5*sy*sys*s*v1I/dr[1])
                                                          -(1.+(alpha-1.)*SQR(n2))
                                                          *theta*dt*nu*(SQR(rx)*v2rr+rx*rxr*v2r+2.*rx*sx*v2rs+.5*rx*sxr*s*v2I/dr[1]+rxs*sx*v2r
                                                                                      +SQR(sx)*v2I/SQR(dr[1])+.5*sx*sxs*s*v2I/dr[1]+SQR(ry)*v2rr+ry*ryr*v2r+2.*ry*sy*v2rs
                                                                                      +.5*ry*syr*s*v2I/dr[1]+rys*sy*v2r+SQR(sy)*v2I/SQR(dr[1])+.5*sy*sys*s*v2I/dr[1]))
                                            +g4;
                    // *testing* try setting to predicted values
                    // A(3,0) = 0.;
                    // A(3,1) = 1.;
                    // A(3,2) = 0.;
                    // A(3,3) = 0.;
                    // b(3) = uLocal(i1,i2,i3,vc)+(alpha-1)*(n1*n2*(uLocal(i1,i2,i3,uc)-vs1)+SQR(n2)*(uLocal(i1,i2,i3,vc)-vs2));
                    // printF("  A = \n");
                    // A.display("");
                    // printF("\n");
                    // printF("  b = \n");
                    // b.display("");
                    // printF("\n");
                    // printF("g1=%e, g2=%e, g3=%e, g4=%e\n",g1,g2,g3,g4);
                    // factor
                                        GECO(A(0,0), numberOfEquations, numberOfEquations, ipvt(0), rcond, work(0));
                    // printF("  rcond = %f\n\n",rcond);
                    // solve
                                        x=b;
                                        GESL(A(0,0), numberOfEquations, numberOfEquations, ipvt(0), x(0), job);
                    // printF("  x = \n");
                    // x.display("");
                    // printF("\n");
                    // *** check ***
                                        v1B = x(0);
                                        v2B = x(1);
                                        v1G = x(2);
                                        v2G = x(3);
                                        v1s = (v1I-v1G)*(1./(2*dr[1]));
                                        v2s = (v2I-v2G)*(1./(2*dr[1]));
                                        v1ss = (v1I-2*v1B+v1G)*(1./(SQR(dr[1])));
                                        v2ss = (v2I-2*v2B+v2G)*(1./(SQR(dr[1])));
                    // calculate residual
                                        real resEq1 = (rx*v1r+ry*v2r+sx*v1s+sy*v2s)*n1
                                            +(1-SQR(n1))*(2*(rx*v1r+sx*v1s)*n1
                                                                        +(rx*v2r+ry*v1r+sx*v2s+sy*v1s)*n2
                                                                        +zs*v1B/mu)
                                            -n1*n2      *(2*(ry*v2r+sy*v2s)*n2
                                                                        +(rx*v2r+ry*v1r+sx*v2s+sy*v1s)*n1
                                                                        +zs*v2B/mu) 
                                            -g1;
                                        real resEq2 = (rx*v1r+ry*v2r+sx*v1s+sy*v2s)*n2
                                            +(1-SQR(n2))*(2*(ry*v2r+sy*v2s)*n2
                                                                        +(rx*v2r+ry*v1r+sx*v2s+sy*v1s)*n1
                                                                        +zs*v2B/mu)
                                            -n1*n2      *(2*(rx*v1r+sx*v1s)*n1
                                                                        +(rx*v2r+ry*v1r+sx*v2s+sy*v1s)*n2
                                                                        +zs*v1B/mu)
                                            -g2;
                                        real resEq3 = v1B-
                                            (1.+(alpha-1.)*SQR(n1))
                                            *theta*dt*nu*(SQR(rx)*v1rr+rx*rxr*v1r+2*rx*sx*v1rs+rx*sxr*v1s+rxs*sx*v1r
                                                                        +SQR(ry)*v1rr+ry*ryr*v1r+2*ry*sy*v1rs+ry*syr*v1s+rys*sy*v1r
                                                                        +SQR(sx)*v1ss+sx*sxs*v1s+SQR(sy)*v1ss+sy*sys*v1s)
                                            -(alpha-1.)*n1*n2
                                            *theta*dt*nu*(SQR(rx)*v2rr+rx*rxr*v2r+2*rx*sx*v2rs+rx*sxr*v2s+rxs*sx*v2r
                                                                        +SQR(ry)*v2rr+ry*ryr*v2r+2*ry*sy*v2rs+ry*syr*v2s+rys*sy*v2r
                                                                        +SQR(sx)*v2ss+sx*sxs*v2s+SQR(sy)*v2ss+sy*sys*v2s)
                                            -g3;
                                        real resEq4 = v2B-
                                            (alpha-1.)*n1*n2
                                            *theta*dt*nu*(SQR(rx)*v1rr+rx*rxr*v1r+2*rx*sx*v1rs+rx*sxr*v1s+rxs*sx*v1r
                                                                        +SQR(ry)*v1rr+ry*ryr*v1r+2*ry*sy*v1rs+ry*syr*v1s+rys*sy*v1r
                                                                        +SQR(sx)*v1ss+sx*sxs*v1s+SQR(sy)*v1ss+sy*sys*v1s)
                                            -(1.+(alpha-1.)*SQR(n2))
                                            *theta*dt*nu*(SQR(rx)*v2rr+rx*rxr*v2r+2*rx*sx*v2rs+rx*sxr*v2s+rxs*sx*v2r
                                                                        +SQR(ry)*v2rr+ry*ryr*v2r+2*ry*sy*v2rs+ry*syr*v2s+rys*sy*v2r
                                                                        +SQR(sx)*v2ss+sx*sxs*v2s+SQR(sy)*v2ss+sy*sys*v2s)
                                            -g4;
                                        printF("New: u,v=%e,%e, ghost=%e,%e\n",v1B,v2B,v1G,v2G);
                                        printF("Old: u,v=%e,%e, ghost=%e,%e\n",uP(i1,i2,i3,uc),uP(i1,i2,i3,vc),uP(i1,i2-1,i3,uc),uP(i1,i2-1,i3,vc));
                                        printF("res1=%e, res2=%e, res3=%e, res4=%e \n",resEq1,resEq2,resEq3,resEq4);
                    // printF("uFp=%e, uS=%e, uF=%e\n",uLocal(i1,i2,i3,uc),v1s,v1B);
                    // printF("vFp=%e, vS=%e, vF=%e\n",uLocal(i1,i2,i3,vc),v2s,v2B);
                    // assign
                                        UB(i1,i2,i3,0) = v1B;
                                        UB(i1,i2,i3,1) = v2B;
                                        UG(i1,i2,i3,0) = v1G;
                                        UG(i1,i2,i3,1) = v2G;
                                    }
                  // ::display(UB(Ib1,Ib2,Ib3,1),"v2 on boundary (not linearized)");
                  // ::display(UB(Ib1,Ib2,Ib3,0),"v1 on boundary (not linearized)");
                  // ::display(UG(Ib1,Ib2,Ib3,0),"v1 on ghost (not linearized)");
                                    uLocal(Ib1,Ib2,Ib3,uc) = UB(Ib1,Ib2,Ib3,0);
                                    uLocal(Ib1,Ib2,Ib3,vc) = UB(Ib1,Ib2,Ib3,1);
                                    uLocal(Ib1,Ib2-1,Ib3,uc) = UG(Ib1,Ib2,Ib3,0);
                                    uLocal(Ib1,Ib2-1,Ib3,vc) = UG(Ib1,Ib2,Ib3,1);
                                    gridVelocityLocal(Ib1,Ib2,Ib3,0)= uLocal(Ib1,Ib2,Ib3,uc);
                                    gridVelocityLocal(Ib1,Ib2,Ib3,1)= uLocal(Ib1,Ib2,Ib3,vc);
                  // check divergence
                  // check traction
                                    OV_ABORT("stop here for now");
                                }
                            }
                            
                        }
                        
                        bool assignExactVelocity=FALSE; // *****************************************
                        if( assignExactVelocity &&  twilightZoneFlow )
                        {
              // -- for testing -- see the exact values 
                            fPrintF(debugFile,"--INS-- PIV: **TEMP** set exact TZ values for velocity on the interface+ghost, ***TEMP*** t=%9.3e\n",t);

                            OV_GET_SERIAL_ARRAY(real,cg[grid].vertex(),xLocal);
                            OGFunction & e = *(parameters.dbase.get<OGFunction* >("exactSolution"));
                            bool isRectangular=false;

                            RealArray ue(Ib1,Ib2,Ib3,V);
                            e.gd( ue,xLocal,numberOfDimensions,isRectangular,0,0,0,0,Ib1,Ib2,Ib3,V,t); 
                            uLocal(Ib1,Ib2,Ib3,V)=ue;  // boundary values 
                            RealArray ueg(Ig1,Ig2,Ig3,V);
                            e.gd( ueg,xLocal,numberOfDimensions,isRectangular,0,0,0,0,Ig1,Ig2,Ig3,V,t);
                            uLocal(Ig1,Ig2,Ig3,V)=ueg;  // ghost values 

                        }
                        if( assignExactVelocity && knownSolution==Parameters::userDefinedKnownSolution )
                        {
                            printF("--INS-- PIV: **TEST** set exact KNOWN-SOLUTION values for velocity on the interface+ghost ***TEMP***, t=%9.3e\n",t);
                            fPrintF(debugFile,"--INS-- PIV: **TEST** set exact KNOWN-SOLUTION values for velocity on the interface+ghost ***TEMP***, t=%9.3e\n",t);
                            int body=0;
                            parameters.getUserDefinedDeformingBodyKnownSolution( body,Parameters::boundaryVelocity,
                                                                                                                                      t, grid, mg, Ib1,Ib2,Ib3,V,uLocal );
                            parameters.getUserDefinedDeformingBodyKnownSolution( body,Parameters::boundaryVelocity,
                                                                                                                                      t, grid, mg, Ig1,Ig2,Ig3,V,uLocal );
                        }

                        

                        if( FALSE )
                        {
                            printF("\n >>>>>>>> PIV **TEMP** TESTING set v_1=0 \n");
                            gridVelocityLocal(Ib1,Ib2,Ib3,0)= 0.; 
                        }
                        
                        if( debug() & 8 )
                        {
                            ::display(gridVelocityLocal(Ib1,Ib2,Ib3,Rx),sPrintF("PROJECTED interface velocity t=%9.3e",t),"%6.2f ");
                            ::display(uLocal(Ib1,Ib2,Ib3,V),sPrintF(" -> fluid portion t=%9.3e",t),"%7.3f ");
                            ::display(vSolidLocal(Ib1,Ib2,Ib3,Rx),sPrintF("--> solid velocity t=%9.3e",t),"%7.3f ");
                            ::display(fabs(uLocal(Ib1,Ib2,Ib3,V)-vSolidLocal(Ib1,Ib2,Ib3,Rx)),sPrintF("--> DIFF fluid-solid velocity t=%9.3e",t),"%8.2e ");
                        }
                        

          	    continue;
        	  }
        	  else if( deform.isBeamModel() )
        	  {

            // **********************************************************
	    // ************ PROJECT VELOCITY BEAM MODEL ******************                
            // **********************************************************
                        {
                        #ifndef USE_PPP
                            BeamModel & beamModel = deform.getBeamModel();
                            real beamMassPerUnitLength=-1.;
                            beamModel.getMassPerUnitLength( beamMassPerUnitLength );
                            real alpha = 1./( 1. + beamMassPerUnitLength/(fluidDensity*fluidAddedMassLengthScale) );
              // alpha=0.; // ***************
                            if( t<=0. )
                                printF("--PIV-- alpha=%8.2e, beamMassPerUnitLength = %8.2e, fluidDensity=%8.2e hf=%8.2e\n",
                                              alpha,beamMassPerUnitLength,fluidDensity,fluidAddedMassLengthScale);
              // --- Extract the "weight" array for weighting the velocity projection ---
              //  This is used when we do not project the velocity on the ends of the beam
                            RealArray *pWeight= &Overture::nullRealArray(); // set pWeight to a default value if it is not used.
                            if( !projectVelocityOnBeamEnds )
                            {
                                DeformingBodyMotion & deformingBody = movingGrids.getDeformingBody(body);
                                DataBase & deformingBodyDataBase = deformingBody.deformingBodyDataBase;
                                const int & numberOfFaces = deformingBodyDataBase.get<int>("numberOfFaces");
                                const IntegerArray & boundaryFaces = deformingBodyDataBase.get<IntegerArray>("boundaryFaces");
                                BeamFluidInterfaceData &  beamFluidInterfaceData = 
                                    deformingBodyDataBase.get<BeamFluidInterfaceData>("beamFluidInterfaceData");
                                int face=-1;
                                for( int face0=0; face0<numberOfFaces; face0++ )
                                {
                                    const int side0=boundaryFaces(0,face0);
                                    const int axis0=boundaryFaces(1,face0);
                                    const int grid0=boundaryFaces(2,face0);
                                    if( grid==grid0 && side==side0 && axis==axis0 )
                                    { 
                                        face=face0;
                                        break;
                                    }
                                }
                                assert( face>=0 && face<numberOfFaces);
                                RealArray *& weightArray = beamFluidInterfaceData.dbase.get<RealArray*>("weightArray");
                                pWeight = &(weightArray[face]);
                            }
                            RealArray & weight = *pWeight;
                            if( (false && t<=max(0.,dt))  || debug() & 4 )
                            {
                                ::display(u(Ib1,Ib2,Ib3,V),sPrintF("---PIV-- fluid velocity at t=%9.3e",t),"%7.2e ");
                                ::display(vSolid(Ib1,Ib2,Ib3,Rx),sPrintF("---PIV-- beam  velocity at t=%9.3e",t),"%7.2e ");
                                ::display(gridVelocity(Ib1,Ib2,Ib3,Rx),sPrintF("---PIV-- grid velocity at t=%9.3e",t),"%7.2e ");
                                if( !projectVelocityOnBeamEnds )
                                    ::display(weight(Ib1,Ib2,Ib3),sPrintF("---PIV-- weight at t=%9.3e",t),"%5.2f ");
                            }
                            if( true ) // project v and adjust the grid velocity
                            {
                                if( projectNormalComponentOfAddedMassVelocity )
                                {
                  // --- only project the normal component of the fluid velocity ---
                  // Project the normal component by subtracting the current normal component and then adding the new
                  //  vp = AMP projected velocity
                  //   v = v - (n.v)n + (n.vp)n
                  //     = v - (n.(vp-v))n 
                  //
                                    if( t <= 10.*dt )
                                        printF("--PIV--: project NORMAL component of velocity only, t=%9.3e\n",t);
                                    if( true )
                                    {
                                        RealArray vp(Ib1,Ib2,Ib3,Rx), nDotV(Ib1,Ib2,Ib3);
                    // vp=( alpha*uLocal(Ib1,Ib2,Ib3,V) + (1.-alpha)*vSolidLocal(Ib1,Ib2,Ib3,Rx) -gridVelocityLocal(Ib1,Ib2,Ib3,Rx) );
                    // vp= (alpha-1.)*uLocal(Ib1,Ib2,Ib3,V) + (1.-alpha)*vSolidLocal(Ib1,Ib2,Ib3,Rx);
                                        if( projectVelocityOnBeamEnds )
                                        {
                                            vp= alpha*uLocal(Ib1,Ib2,Ib3,V) + (1.-alpha)*vSolidLocal(Ib1,Ib2,Ib3,Rx);
                                        }
                                        else
                                        {
                      // use solid velocity when the weight is zero: 
                                            for( int dir=0; dir<numberOfDimensions; dir++ )
                                                vp(Ib1,Ib2,Ib3,dir) = ( (   alpha*weight(Ib1,Ib2,Ib3))*uLocal(Ib1,Ib2,Ib3,uc+dir) + 
                                                                                                (1.-alpha*weight(Ib1,Ib2,Ib3))*vSolidLocal(Ib1,Ib2,Ib3,dir) );
                                        }
                                        nDotV = (normal(Ib1,Ib2,Ib3,0)*vp(Ib1,Ib2,Ib3,0)+
                                                          normal(Ib1,Ib2,Ib3,1)*vp(Ib1,Ib2,Ib3,1) );
                                        if( numberOfDimensions==3 )
                                            nDotV += normal(Ib1,Ib2,Ib3,2)*vp(Ib1,Ib2,Ib3,2);
                                        if( TRUE ) // *WDH* try this 2015/03/06
                                        {
                      // t.v = t.vs 
                                            gridVelocityLocal(Ib1,Ib2,Ib3,Rx)=vSolidLocal(Ib1,Ib2,Ib3,Rx);  // set all components equal to vs 
                                            nDotV -= (normal(Ib1,Ib2,Ib3,0)*vSolidLocal(Ib1,Ib2,Ib3,0)+
                                                                normal(Ib1,Ib2,Ib3,1)*vSolidLocal(Ib1,Ib2,Ib3,1) );
                                            if( numberOfDimensions==3 )
                                                nDotV -= normal(Ib1,Ib2,Ib3,2)*vSolidLocal(Ib1,Ib2,Ib3,2);
                                            for( int dir=0; dir<numberOfDimensions; dir++ )
                                                gridVelocityLocal(Ib1,Ib2,Ib3,dir) += nDotV*normal(Ib1,Ib2,Ib3,dir);  // n.v = n.vp 
                      // gridVelocityLocal(Ib1,Ib2,Ib3,Rx)=vSolidLocal(Ib1,Ib2,Ib3,Rx);  // **********
                                        }
                                        else
                                        {
                      // t.v= ZERO
                      // gridVelocityLocal(Ib1,Ib2,Ib3,Rx)=vSolidLocal(Ib1,Ib2,Ib3,Rx);
                      // gridVelocityLocal(Ib1,Ib2,Ib3,Rx)=0.;
                                            for( int dir=0; dir<numberOfDimensions; dir++ )
                                                gridVelocityLocal(Ib1,Ib2,Ib3,dir) = nDotV*normal(Ib1,Ib2,Ib3,dir);
                                        }
                                    }
                                    else // ** FALSE ***
                                    {
                    // new way
                                        RealArray vp(Ib1,Ib2,Ib3,Rx), nDotV(Ib1,Ib2,Ib3);
                    // vp - v  ( note (alpha-1.) in first term )
                                        if( projectVelocityOnBeamEnds )
                                        {
                      // vp= (alpha-1.)*uLocal(Ib1,Ib2,Ib3,V) + (1.-alpha)*vSolidLocal(Ib1,Ib2,Ib3,Rx);
                                            vp= (alpha)*uLocal(Ib1,Ib2,Ib3,V) + (1.-alpha)*vSolidLocal(Ib1,Ib2,Ib3,Rx);
                                        }
                                        else
                                        {
                      // use solid velocity when the weight is zero: 
                                            for( int dir=0; dir<numberOfDimensions; dir++ )
                                                vp(Ib1,Ib2,Ib3,dir) = ( (   alpha*weight(Ib1,Ib2,Ib3))*uLocal(Ib1,Ib2,Ib3,uc+dir) + 
                                                                                                (1.-alpha*weight(Ib1,Ib2,Ib3))*vSolidLocal(Ib1,Ib2,Ib3,dir) );
                                        }
                    // vp=( alpha*uLocal(Ib1,Ib2,Ib3,V) + (1.-alpha)*vSolidLocal(Ib1,Ib2,Ib3,Rx) 
                    // 		 -gridVelocityLocal(Ib1,Ib2,Ib3,Rx) );
                                        nDotV = (normal(Ib1,Ib2,Ib3,0)*vp(Ib1,Ib2,Ib3,0)+
                                                          normal(Ib1,Ib2,Ib3,1)*vp(Ib1,Ib2,Ib3,1) );
                                        if( numberOfDimensions==3 )
                                            nDotV += normal(Ib1,Ib2,Ib3,2)*vp(Ib1,Ib2,Ib3,2);
                    // if( !projectVelocityOnBeamEnds ) 
                    // {
                    //   nDotV *= weight(Ib1,Ib2,Ib3);  // turn off projection near beam ends 
                    // }
                    // gridVelocityLocal(Ib1,Ib2,Ib3,Rx)=vSolidLocal(Ib1,Ib2,Ib3,Rx);
                    // -- set the normal component of the velocity ---
                    // -- set tangential component of the velocity to zero --
                                        for( int dir=0; dir<numberOfDimensions; dir++ )
                                            gridVelocityLocal(Ib1,Ib2,Ib3,dir) = nDotV*normal(Ib1,Ib2,Ib3,dir);
                                    }
                                }
                                else
                                {
                                    if( projectVelocityOnBeamEnds ) 
                                    {
                                        gridVelocityLocal(Ib1,Ib2,Ib3,Rx) = alpha*uLocal(Ib1,Ib2,Ib3,V) + (1.-alpha)*vSolidLocal(Ib1,Ib2,Ib3,Rx);
                                    }
                                    else
                                    {
                                        for( int dir=0; dir<numberOfDimensions; dir++ )
                                        {
                      // use solid velocity when the weight is zero
                                            gridVelocityLocal(Ib1,Ib2,Ib3,dir) = ( alpha*weight(Ib1,Ib2,Ib3)*uLocal(Ib1,Ib2,Ib3,uc+dir) + 
                                                                                                                          (1.-alpha*weight(Ib1,Ib2,Ib3))*vSolidLocal(Ib1,Ib2,Ib3,dir) );
                                        }
                                    }
                                }
                            } // end if false
                            if( true )
                            {
                                if( t<=0. )
                                    printF("--PIV-- ****TEST*** set gridVelocity=0 on ends\n");
                                Index I1,I2,I3;
                                getIndex(mg.gridIndexRange(),I1,I2,I3);
                                const int axisp1 = (axis+1) % numberOfDimensions;
                                assert( axisp1==0 );
                                for( int sidea=0; sidea<=1; sidea++ )
                                {
                  // *** FINISH ME ***
                                    if( mg.boundaryCondition(sidea,axisp1)==Parameters::noSlipWall )
                                    {
                                        int i1 = sidea==0 ? Ib1.getBase() : Ib1.getBound();
                                        gridVelocityLocal(i1,I2,I3,Rx)=0.;  // set values on WHOLE FACE
                                    }
                                    else if( mg.boundaryCondition(sidea,axisp1)==Parameters::slipWall )
                                    {
                                        int i1 = sidea==0 ? Ib1.getBase() : Ib1.getBound();
                                        gridVelocityLocal(i1,I2,I3,0)=0.; // set values on WHOLE FACE
                                    }
                                    else if( mg.boundaryCondition(sidea,axisp1)==InsParameters::inflowWithPressureAndTangentialVelocityGiven ||
                                                      mg.boundaryCondition(sidea,axisp1)==InsParameters::inflowWithVelocityGiven ) // **FINISH ME**
                                    {
                                        int i1 = sidea==0 ? Ib1.getBase() : Ib1.getBound();
                                        gridVelocityLocal(i1,Ib2,Ib3,Rx)=0.; // set values on end point
                                    }
                                }
                            }
                        #else
                            OV_ABORT("FINISH ME FOR PARALLEL");
                        #endif
                        }
                        

        	  
        	  } // end deform.isBeamModel() ************* END PROJECT BEAM MODEL ******************

          // *********** FIX ME -- THIS IS DUPLICATED ************
          // -- Add a fourth-order filter to interface velocity --
        	  const bool & smoothInterfaceVelocity = parameters.dbase.get<bool>("smoothInterfaceVelocity");
        	  const int numberOfInterfaceVelocitySmooths=parameters.dbase.get<int>("numberOfInterfaceVelocitySmooths");
        	  if( smoothInterfaceVelocity )
        	  {
          	    const real omega=1.; // .5;
	    // real omega=.125; 
                        if( t <= 10.*dt )
                        {
            	      printF("--PIV--: smooth interface velocity, numberOSmooths=%i (4th order filter, omega=%g) grid=%i t=%9.3e...\n",
                 		     numberOfInterfaceVelocitySmooths,omega,grid,t);
                            if( debug()& 4 )
                            {
                                fPrintF(debugFile,"--PIV--: smooth interface velocity, numberOSmooths=%i (4th order filter, omega=%g) grid=%i t=%9.3e...\n",
                 		     numberOfInterfaceVelocitySmooths,omega,grid,t);
                            }
                            
                        }
                        
                        int extra=-1;
                        getBoundaryIndex(mg.gridIndexRange(),side,axis,Ib1,Ib2,Ib3,extra); // leave off end points
          	    Range Rx=numberOfDimensions;

          	    assert( numberOfDimensions==2 );  // *FIX ME for 3D*
          	    
          	    int isv[3], &is1=isv[0], &is2=isv[1], &is3=isv[2];
          	    is1=is2=is3=0;
                        const int axisp1 = (axis+1) % numberOfDimensions;
          	    isv[axisp1]=1;
                        realArray & gv = gridVelocity;
          	    for( int smooth=0; smooth<numberOfInterfaceVelocitySmooths; smooth++ )
          	    {
              // ADJACENT boundary conditions **FINISH ME**
            	      if( true )
            	      {
                                assert( axisp1==0 );
            		int i1a=mg.gridIndexRange(0,0), i1b=mg.gridIndexRange(1,0);
                // -- extrapolate ghost points ---
            		gv(i1a-1,Ib2,Ib3,Rx)=3.*gv(i1a,Ib2,Ib3,Rx)-3.*gv(i1a+1,Ib2,Ib3,Rx)+gv(i1a+2,Ib2,Ib3,Rx);
            		gv(i1b+1,Ib2,Ib3,Rx)=3.*gv(i1b,Ib2,Ib3,Rx)-3.*gv(i1b-1,Ib2,Ib3,Rx)+gv(i1b-2,Ib2,Ib3,Rx);
            	      }
            	      

	      // smooth interface values
              // NOTE: for now we smooth all components of the velocity
                            gv(Ib1,Ib2,Ib3,Rx)= gv(Ib1,Ib2,Ib3,Rx) + 
            		(omega/16.)*( -   gv(Ib1-2*is1,Ib2-2*is2,Ib3,Rx) 
                        			      +4.*gv(Ib1-  is1,Ib2-  is2,Ib3,Rx) 
                        			      -6.*gv(Ib1,      Ib2      ,Ib3,Rx) 
                        			      +4.*gv(Ib1+  is1,Ib2+  is2,Ib3,Rx) 
                        			      -   gv(Ib1+2*is1,Ib2+2*is2,Ib3,Rx) );
          	    } // end smooths
        	  } // end smoothSurface
        	  


      	}
            }
        }
      	

    } // end if bd.dbase.has_key("deformingBodyNumber") )


    return interfaceWasAssigned;
}



//============================================================================================
/// \brief This function is called from applyBoundaryConditions to assign some 
/// interface conditions (e.g. velocity projection for beams) that are not handled by cgmp. 
//============================================================================================
int Cgins::
assignInterfaceBoundaryConditions(GridFunction & cgf,
                          				  const int & option /* =-1 */,
                          				  int grid_ /* = -1 */,
                          				  GridFunction *puOld /* =NULL */, 
                          				  const real & dt /* =-1. */ )
{
    const bool & useAddedMassAlgorithm = parameters.dbase.get<bool>("useAddedMassAlgorithm");
    const bool & projectAddedMassVelocity = parameters.dbase.get<bool>("projectAddedMassVelocity");
    const int initialConditionsAreBeingProjected = parameters.dbase.get<int>("initialConditionsAreBeingProjected");
    const bool & projectBeamVelocity = parameters.dbase.get<bool>("projectBeamVelocity");
    if( useAddedMassAlgorithm && projectAddedMassVelocity && !initialConditionsAreBeingProjected 
            && projectBeamVelocity 
            && cgf.t>0. )
    {

    // --- project the velocity of the beam to match that from the fluid ---
        MovingGrids & movingGrids = parameters.dbase.get<MovingGrids >("movingGrids");
        const int numberOfDeformingBodies= movingGrids.getNumberOfDeformingBodies();

    // ----- We only project defomring body interfaces for now ----
        if( numberOfDeformingBodies==0 )
            return 0;


        if( cgf.t <= 5*dt ) 
            printF("--INS-- assignInterfaceBoundaryConditions: PROJECT-INTERFACE-VELOCITY at t=%8.2e\n",cgf.t);

        movingGrids.projectInterfaceVelocity( cgf );

        if( true  ) // this could be an option : reprojectFluidVelocity
        {
      // ---- After projecting the solid  (beam) velocity we need to make sure the
      //   the fluid velocity on both sides of the beam matches the new beam velocity.
      //   This ensures the fluid velocity on opposite sides of the beam are now consistent with the single beam velocity,
      //   otherwise the fluid velocity on opposite sides of the beam could get out of sync.
            CompositeGrid & cg = cgf.cg;
            const int numberOfDimensions = cg.numberOfDimensions();
            const int uc = parameters.dbase.get<int >("uc");
            Range V = Range(uc,uc+numberOfDimensions-1);

            for( int body=0; body<numberOfDeformingBodies; body++ )
            {
      	DeformingBodyMotion & deformingBody = movingGrids.getDeformingBody(body);
                const bool beamModelHasFluidOnTwoSides = deformingBody.beamModelHasFluidOnTwoSides();
      	if( beamModelHasFluidOnTwoSides )
      	{
          // --- this body is a beam with fluid on two sides ---
        	  DataBase & deformingBodyDataBase = deformingBody.deformingBodyDataBase;
        	  const int & numberOfFaces = deformingBodyDataBase.get<int>("numberOfFaces");
        	  const IntegerArray & boundaryFaces = deformingBodyDataBase.get<IntegerArray>("boundaryFaces");
          	    

        	  BeamFluidInterfaceData &  beamFluidInterfaceData = 
          	    deformingBodyDataBase.get<BeamFluidInterfaceData>("beamFluidInterfaceData");
        	  Index Ib1,Ib2,Ib3;
        	  for( int face=0; face<numberOfFaces; face++ )
        	  {
          	    const int side=boundaryFaces(0,face);
          	    const int axis=boundaryFaces(1,face);
          	    const int grid=boundaryFaces(2,face);
          	    MappedGrid & mg = cg[grid];
          	    OV_GET_SERIAL_ARRAY(real,cgf.u[grid],uLocal);
          	    getBoundaryIndex(mg.gridIndexRange(),side,axis,Ib1,Ib2,Ib3); // boundary index's for mg

          	    Range Rx=numberOfDimensions;
          	    RealArray vSolid(Ib1,Ib2,Ib3,Rx); // holds velocity of solid on the boundary
          	    deformingBody.getVelocityBC( cgf.t, side,axis,grid, mg, Ib1,Ib2,Ib3, vSolid );

          	    uLocal(Ib1,Ib2,Ib3,V)=vSolid(Ib1,Ib2,Ib3,Rx);


            // *********** FIX ME -- THIS IS DUPLICATED ************
	    // -- Add a fourth-order filter to interface velocity --
          	    const bool & smoothInterfaceVelocity = parameters.dbase.get<bool>("smoothInterfaceVelocity");
          	    const int numberOfInterfaceVelocitySmooths=parameters.dbase.get<int>("numberOfInterfaceVelocitySmooths");
          	    if( smoothInterfaceVelocity )
          	    {
            	      const real omega=1.; // .5 
	      // real omega=.125; 
            	      if( cgf.t <= 10.*dt )
            		printF("--IBC--: smooth interface velocity, numberOSmooths=%i (4th order filter, omega=%g) grid=%i t=%9.3e...\n",
                   		       numberOfInterfaceVelocitySmooths,omega,grid,cgf.t);
          	    
            	      int extra=-1;
            	      getBoundaryIndex(mg.gridIndexRange(),side,axis,Ib1,Ib2,Ib3,extra); // leave off end points

            	      assert( numberOfDimensions==2 );  // *FIX ME for 3D*
          	    
            	      int isv[3], &is1=isv[0], &is2=isv[1], &is3=isv[2];
            	      is1=is2=is3=0;
            	      const int axisp1 = (axis+1) % numberOfDimensions;
            	      isv[axisp1]=1;
            	      RealArray & gv = uLocal;
            	      for( int smooth=0; smooth<numberOfInterfaceVelocitySmooths; smooth++ )
            	      {
		// ADJACENT boundary conditions **FINISH ME**
            		if( true )
            		{
              		  assert( axisp1==0 );
              		  int i1a=mg.gridIndexRange(0,0), i1b=mg.gridIndexRange(1,0);
		  // -- extrapolate ghost points ---
              		  gv(i1a-1,Ib2,Ib3,V)=3.*gv(i1a,Ib2,Ib3,V)-3.*gv(i1a+1,Ib2,Ib3,V)+gv(i1a+2,Ib2,Ib3,V);
              		  gv(i1b+1,Ib2,Ib3,V)=3.*gv(i1b,Ib2,Ib3,V)-3.*gv(i1b-1,Ib2,Ib3,V)+gv(i1b-2,Ib2,Ib3,V);
            		}
            	      

		// smooth interface values
		// NOTE: for now we smooth all components of the velocity
            		gv(Ib1,Ib2,Ib3,V)= gv(Ib1,Ib2,Ib3,V) + 
              		  (omega/16.)*( -   gv(Ib1-2*is1,Ib2-2*is2,Ib3,V) 
                        				+4.*gv(Ib1-  is1,Ib2-  is2,Ib3,V) 
                        				-6.*gv(Ib1,      Ib2      ,Ib3,V) 
                        				+4.*gv(Ib1+  is1,Ib2+  is2,Ib3,V) 
                        				-   gv(Ib1+2*is1,Ib2+2*is2,Ib3,V) );
            	      } // end smooths
          	    } // end smoothSurface


        	  }
        	  
      	} // end if beamModelHasFluidOnTwoSides
      	
            } // end for body
        }
        
    }
    
    return 0;
}
