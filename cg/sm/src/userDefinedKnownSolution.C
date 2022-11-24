// This file automatically generated from userDefinedKnownSolution.bC with bpp.
#include "SmParameters.h"
#include "FlowSolutions.h"
#include "GenericGraphicsInterface.h"
#include "FluidPiston.h"
#include "PistonMotion.h"
#include "ParallelUtility.h"
#include "TimeFunction.h"
#include "BeamModel.h"

#include "incompressible/SmCylinderExactSolution.h"
#include "incompressible/SmSphereExactSolution.h"
#include "incompressible/SmRectangleExactSolution.h"


#define rotatingDiskSVK   EXTERN_C_NAME(rotatingdisksvk)
#define evalFibShearSolid EXTERN_C_NAME(evalfibshearsolid)
#define evalRadialFibShearSolid EXTERN_C_NAME(evalradialfibshearsolid)
#define evalFibShearSolidFull EXTERN_C_NAME(evalfibshearsolidfull)
#define evalFibCartWaveSolid EXTERN_C_NAME(evalfibcartwavesolid)
#define evalFibRadialWaveSolid EXTERN_C_NAME(evalfibradialwavesolid)

extern "C"
{
  // rotating disk (SVK) exact solution:
    void rotatingDiskSVK( const real & t, const int & numberOfGridPoints, real & uDisk, real & param,
                                                const int & nrwk, real & rwk );

  // exact fsi solution for shearing solid
    void evalFibShearSolid( const real & ksr, const real & ksi,
                                                    const real & ar, const real & ai,
                                                    const real & br, const real & bi,
                                                    const real & y, 
                                                    real & ur, real & ui, real & uyr, real & uyi);

    void evalFibShearSolidFull( const real & ksr, const real & ksi,
                                                            const real & ar, const real & ai,
                                                            const real & br, const real & bi,
                                                            const real & y, const real & t,
                                                            real & ur, real & ui, 
                                                            real & vr, real & vi, 
                                                            real & uyr, real & uyi,
                                                            const real & omegar, const real & omegai);

    void evalRadialFibShearSolid( const real & ksr, const real & ksi,
                                                                const real & omegar, const real & omegai,
                                                                const real & cr, const real & ci,
                                                                const real & r, const real & t,
                                                                real & ut, real & vt, real & at, 
                                                                real & utr, real & vtr);

    void evalFibCartWaveSolid( const real & omegar, const real & omegai,
                                                          const real & k, const real & mubar,
                                                          const real & rhobar, const real & lambdabar,
                                                          const real & k1r, const real & k1i,
                                                          const real & k2r, const real & k2i,
                                                          const real & amp, const real & x,
                                                          const real & y, const real & t,
                                                          const real & Hbar, 
                                                          real & u1barr, real & u2barr, real & v1barr, 
                                                          real & v2barr, real & s11barr, 
                                                          real & s12barr, real & s22barr);

    void evalFibRadialWaveSolid( const real & omegar, const real & omegai,
                                                              const real & d1r, const real & d1i,
                                                              const real & d2r, const real & d2i,
                                                              const real & d3r, const real & d3i,
                                                              const real & d4r, const real & d4i,
                                                              const real & n, const real & rhobar,
                                                              const real & muBar, const real & lambdaBar,
                                                              const real & r,  const real & theta,
                                                              const real & t,
                                                              real & ur, real & ut, real & vr, real & vt, 
                                                              real & ar, real & at,
                                                              real & srr, real & srt, real & stt, 
                                                              real & sdrr, real & sdrt, real & sdtt);
    
}

#define FOR_3D(i1,i2,i3,I1,I2,I3) int I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase();  int I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(); for(int i3=I3Base; i3<=I3Bound; i3++) for(int i2=I2Base; i2<=I2Bound; i2++) for(int i1=I1Base; i1<=I1Bound; i1++)

#define FOR_3(i1,i2,i3,I1,I2,I3) I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase();  I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(); for(int i3=I3Base; i3<=I3Bound; i3++) for(int i2=I2Base; i2<=I2Bound; i2++) for(int i1=I1Base; i1<=I1Bound; i1++)  

// Macro to get the vertex array
#define GET_VERTEX_ARRAY(x)                                     mg.update(MappedGrid::THEvertex | MappedGrid::THEcenter);       OV_GET_SERIAL_ARRAY_CONST(real,mg.vertex(),x);


// --------------------------------------------------------------------------------
// Assign the incompressible annulus solution
// Output: uNorm : Return the norm of the solution
// --------------------------------------------------------------------------------


// ----------------------------------------------------------------
// Evaluate the polar cooridinates from the cartesian coordinates 
// ----------------------------------------------------------------


// --------------------------------------------------------------------------------
// Assign the incompressible disk solution
// --------------------------------------------------------------------------------

// =======================================================================
// Macro: INITIALIZE variables for the method of lines approach
// This should only be done at the INITIAL TIME (identified by t<0)
//      velocity
//      acceleration
//      pressure 
// =======================================================================


// =======================================================================
// Macro: INITIALIZE variables for the method of lines approach
// This should only be done at the INITIAL TIME (identified by t<0)
//      velocity
//      acceleration
//      pressure 
// =======================================================================


// =======================================================================
// Macro: INITIALIZE variables for the method of lines approach
// This should only be done at the INITIAL TIME (identified by t<0)
//      velocity
//      acceleration
//      pressure 
// =======================================================================



// =======================================================================
// Macro: INITIALIZE variables for the method of lines approach
// This should only be done at the INITIAL TIME (identified by t<0)
//      velocity
//      acceleration
//      pressure 
// =======================================================================

int SmParameters::
getUserDefinedKnownSolution(real t, CompositeGrid & cg, int grid, RealArray & ua, 
                                                        const Index & I1, const Index &I2, const Index &I3, 
                                                        int numberOfTimeDerivatives /* = 0 */  )
// ==========================================================================================
///  \brief Evaluate a user defined known solution.
// ==========================================================================================
{
    MappedGrid & mg = cg[grid];
    const int numberOfDimensions = cg.numberOfDimensions();

    if( ! dbase.get<DataBase >("modelData").has_key("userDefinedKnownSolutionData") )
    {
        printf("getUserDefinedKnownSolution:ERROR: sub-directory `userDefinedKnownSolutionData' not found!\n");
        Overture::abort("error");
    }
    DataBase & db =  dbase.get<DataBase >("modelData").get<DataBase>("userDefinedKnownSolutionData");

    const aString & userKnownSolution = db.get<aString>("userKnownSolution");

    real *rpar = db.get<real[20]>("rpar");
    int *ipar = db.get<int[20]>("ipar");

    const real dt = dbase.get<real>("dt");

    const real & rho    = dbase.get<real>("rho");
    const real & mu     = dbase.get<real>("mu");
    const real & lambda = dbase.get<real>("lambda");

    const real cp=sqrt((2.*mu+lambda)/rho);
    const real cs =sqrt(mu/rho);

  // Here are the comopnents for displacement velocity and stress
    const int v1c = dbase.get<int >("v1c");
    const int v2c = dbase.get<int >("v2c");
    const int v3c = dbase.get<int >("v3c");

    const int u1c = dbase.get<int >("u1c");
    const int u2c = dbase.get<int >("u2c");
    const int u3c = dbase.get<int >("u3c");

    const int s11c = dbase.get<int >("s11c");
    const int s12c = dbase.get<int >("s12c");
    const int s21c = dbase.get<int >("s21c");
    const int s22c = dbase.get<int >("s22c");

    const bool assignVelocity = v1c>=0 ;
    const bool assignStress   = s11c>=0 ;

    SmParameters::TimeSteppingMethodSm & timeSteppingMethodSm   = dbase.get<SmParameters::TimeSteppingMethodSm>("timeSteppingMethodSm");
    SmParameters::CompressibilityTypeEnum & compressibilityType = dbase.get<SmParameters::CompressibilityTypeEnum>("compressibilityType");   

    bool computeVelocity=false;
    realCompositeGridFunction *vgf  = NULL;
    realCompositeGridFunction *vtgf = NULL;
    realCompositeGridFunction *pgf  = NULL;
    int numberOfVelocityFunctions =0;
    int numberOfPressureFunctions =0;

  // ----- NOTE: only do this when t<0 -----
    if( t<0. && compressibilityType==SmParameters::incompressibleSolid && timeSteppingMethodSm!=SmParameters::modifiedEquationTimeStepping )
    {
    // compute the velocity for MOL schemes
        computeVelocity=true; 
        vgf                       = dbase.get<realCompositeGridFunction*>("vgf");  // holds v 
        vtgf                      = dbase.get<realCompositeGridFunction*>("vtgf"); // holds vt
        numberOfVelocityFunctions = dbase.get<int>("numberOfVelocityFunctions");

        pgf                       = dbase.get<realCompositeGridFunction*>("pgf");  // holds p for time extrap
        numberOfPressureFunctions = dbase.get<int>("numberOfPressureFunctions");

        const int & numberOfAccelerationFunctions = dbase.get<int>("numberOfAccelerationFunctions");
        assert( numberOfVelocityFunctions==numberOfAccelerationFunctions );

        printF("\n ##########  getUserDefinedKnownSolution: EVAL VELOCITY AND ACCELERATION: t=%9.3e, dt=%9.3e, numberOfVelocityFunctions=%d,"
                      " numberOfPressureFunctions=%d ##########\n",
                      t,dt,numberOfVelocityFunctions,numberOfPressureFunctions);
    }


  //  assert( numberOfTimeDerivatives==0 );  // for now we don't use this in Cgsm

    if( userKnownSolution=="rotatingDisk" )
    {
    // ---- return the exact solution for the rotating disk ---
        printF(" userDefinedKnownSolution: rotatingDisk: t=%9.3e\n",t);

        assert( v1c>=0 && u1c>=0 && s22c >=0 );

        MappedGrid & mg = cg[grid];
        mg.update(MappedGrid::THEvertex | MappedGrid::THEcenter);

        const realArray & center = mg.center();
        RealArray & u = ua;

        
    // tDisk : solution has been computed at this time
        real & tDisk = db.get<real>("tDisk");

        const int & numberOfGridPoints = db.get<int>("numberOfGridPointsDisk");

        const real & omega = db.get<real>("omegaDisk");

        const real & innerRadius = db.get<real>("innerRadiusDisk");
        const real & outerRadius = db.get<real>("outerRadiusDisk");


        const real dr = (outerRadius-innerRadius)/(numberOfGridPoints-1.);

        RealArray & uDisk = db.get<RealArray>("uDisk"); // exact solution is stored here

        if( t!=tDisk )
        {
      // Compute the exact solution at time t
            tDisk=t;


            int numberOfComponents=8;  // number of components needed in 1D solution: we compute u, v, and sigma 
            if( uDisk.getLength(0)!=numberOfGridPoints )
            {
                uDisk.redim(numberOfGridPoints,numberOfComponents);
            }
            
      // call the fortran routine here to evaluate the solution at time t and return the result in uDisk     

            int nrwk=10*numberOfGridPoints+6;
            RealArray rwk(nrwk);
            int npar=10;
            RealArray param(npar);

            param(0)=innerRadius;
            param(1)=outerRadius;
            param(2)=omega;
            param(3)=1.;   // this is lambda
            param(4)=1.;   // this mu

            rotatingDiskSVK( t, numberOfGridPoints, *uDisk.getDataPointer(), *param.getDataPointer(), nrwk, *rwk.getDataPointer() );

        }
        

        const real twopi=6.283185307179586;
        const real x0=0., y0=0.;  // center of the disk 

        int i1,i2,i3;
        FOR_3D(i1,i2,i3,I1,I2,I3)
        {
      // Reference coordinates:
            real x= center(i1,i2,i3,0);
            real y= center(i1,i2,i3,1);

            real r = sqrt( SQR(x-x0) + SQR(y-y0) );

      // closest point in uDisk less than or equal to r:
            int i = int( (r-innerRadius)/dr );  
            i = max(0,min(i,numberOfGridPoints-2));
            
      // linear interpolation
            real alpha = (r-innerRadius)/dr - i;
            real c1=1.-alpha, c2=alpha;

      // radial and angular displacements
            real w=c1*uDisk(i,0)+c2*uDisk(i+1,0);
            real p=c1*uDisk(i,1)+c2*uDisk(i+1,1);

      // angular position
            real theta=0.;
            if (r>0.)
            {
                if (y<0.)
                {
                    theta=twopi-acos((x-x0)/r);
                }
                else
                {
                    theta=acos((x-x0)/r);
                }
            }

//      printF("rotatingDisk: i=%i, uDisk(i,0)=%9.3e, uDisk(i+1,0)=%9.3e, uDisk(i,1)=%9.3e, uDisk(i+1,1)=%9.3e\n",i,uDisk(i,0),uDisk(i+1,0),uDisk(i,1),uDisk(i+1,1));
//      printF("rotatingDisk: i=%i, alpha=%9.3e, c1=%9.3e, c2=%9.3e, w=%9.3e, p=%9.3e, theta=%9.3e\n",i,alpha,c1,c2,w,p,theta);

      // displacements
            real thbar=theta+p;
            u(i1,i2,i3,u1c)=(r+w)*cos(thbar)-x;
            u(i1,i2,i3,u2c)=(r+w)*sin(thbar)-y;

      // velocities
            real wt=c1*uDisk(i,2)+c2*uDisk(i+1,2);
            real pt=c1*uDisk(i,3)+c2*uDisk(i+1,3);
            u(i1,i2,i3,v1c)=wt*cos(thbar)-pt*(r+w)*sin(thbar);
            u(i1,i2,i3,v2c)=wt*sin(thbar)+pt*(r+w)*cos(thbar);

      // stresses
            real pbar11=c1*uDisk(i,4)+c2*uDisk(i+1,4);
            real pbar12=c1*uDisk(i,5)+c2*uDisk(i+1,5);
            real pbar21=c1*uDisk(i,6)+c2*uDisk(i+1,6);
            real pbar22=c1*uDisk(i,7)+c2*uDisk(i+1,7);
            real pt11=pbar11*cos(thbar)-pbar12*sin(thbar);
            real pt12=pbar11*sin(thbar)+pbar12*cos(thbar);
            real pt21=pbar21*cos(thbar)-pbar22*sin(thbar);
            real pt22=pbar21*sin(thbar)+pbar22*cos(thbar);
            u(i1,i2,i3,s11c)=cos(theta)*pt11-sin(theta)*pt21;
            u(i1,i2,i3,s12c)=cos(theta)*pt12-sin(theta)*pt22;
            u(i1,i2,i3,s21c)=sin(theta)*pt11+cos(theta)*pt21;
            u(i1,i2,i3,s22c)=sin(theta)*pt12+cos(theta)*pt22;

        }


    }

    else if ( userKnownSolution=="incompressibleHollowCylinder" ||
                        userKnownSolution=="incompressibleSolidCylinder" ) 
    {
        SmCylinderExactSolution & cylExact = dbase.get<SmCylinderExactSolution>("cylExact");

        printF("UDKS: evaluate hollow or solid cylinder exact solution at t=%9.3e\n",t);

        cylExact.evalSolution( t, cg, grid, ua, I1,I2,I3, numberOfTimeDerivatives );
        if( computeVelocity )
        {
            printF("\n $$$$$$$$$$ userDefinedKnownSolution: EVALUATE VELOCITY, ACCEL and PRESSURE AT PAST TIMES t=%9.3e dt=%16.6e $$$$$$$\n\n",t,dt);
            assert( dt>0. );
            
            OV_ABORT("UDKS: FINISH ME");

        }    

    }

    else if ( userKnownSolution=="incompressibleHollowSphere" ||
                        userKnownSolution=="incompressibleSolidSphere" ) 
    {
        SmSphereExactSolution & sphereExact = dbase.get<SmSphereExactSolution>("sphereExact");

        if( t<2*dt )
            printF("UDKS: evaluate hollow or solid sphere exact solution at t=%9.3e\n",t);

        sphereExact.evalSolution( t, cg, grid, ua, I1,I2,I3, numberOfTimeDerivatives );
        if( computeVelocity )
        {
            printF("\n $$$$$$$$$$ userDefinedKnownSolution: EVALUATE VELOCITY, ACCEL and PRESSURE AT PAST TIMES t=%9.3e dt=%16.6e $$$$$$$\n\n",t,dt);
            assert( dt>0. );
            
            OV_ABORT("UDKS: FINISH ME");

        }    

    } 

    else if ( userKnownSolution=="incompressibleRectangleEigenmode" )
    {
        SmRectangleExactSolution & rectangleExact = dbase.get<SmRectangleExactSolution>("rectangleExact");

        if( t<2*dt )
            printF("UDKS: evaluate rectangle exact solution at t=%9.3e\n",t);

        rectangleExact.evalSolution( t, cg, grid, ua, I1,I2,I3, numberOfTimeDerivatives );

            if( computeVelocity )
            {
                printF("\n $$$$$$$$$$ userDefinedKnownSolution: EVALUATE VELOCITY, ACCEL and PRESSURE AT PAST TIMES for rectangleExact "
                              "t=%9.3e dt=%16.6e $$$$$$$\n\n",t,dt);
                assert( dt>0. );
                assert( vgf!=NULL );
                const int currentVelocity = dbase.get<int>("currentVelocity");
                assert( currentVelocity>=0 && currentVelocity<numberOfVelocityFunctions );
        // fill in t0=0 and some past time levels
                for( int level=0; level<numberOfVelocityFunctions; level++ )
                {
                    real tc = 0. - level*dt;
                    const int prevVelocity = (currentVelocity - level + numberOfVelocityFunctions) % numberOfVelocityFunctions;
          // eval velocity 
                    int numberOfTimeDerivatives=1; 
                    OV_GET_SERIAL_ARRAY(real,vgf[prevVelocity][grid],vgfLocal);
                    rectangleExact.evalSolution( tc, cg, grid, vgfLocal, I1,I2,I3, numberOfTimeDerivatives );
          // rectangleExact.evalSolution( tc, cg, grid, vgf[prevVelocity][grid], I1,I2,I3, numberOfTimeDerivatives );
          // eval acceleration      
                    numberOfTimeDerivatives=2; 
                    OV_GET_SERIAL_ARRAY(real,vtgf[prevVelocity][grid],vtgfLocal);
                    rectangleExact.evalSolution( tc, cg, grid, vtgfLocal, I1,I2,I3, numberOfTimeDerivatives );
          // rectangleExact.evalSolution( tc, cg, grid, vtgf[prevVelocity][grid], I1,I2,I3, numberOfTimeDerivatives );
                }
        // --- save pressure for time extrapolation ---
                const int currentPressure = dbase.get<int>("currentPressure");
                assert( currentPressure>=0 && currentPressure<max(1,numberOfPressureFunctions) );
                for( int level=0; level<numberOfPressureFunctions; level++ )
                {
                    real tc = 0. - (level+3)*dt;  // -- save pressure starting at t-3*dt 
                    const int prevPressure = (currentPressure - level + numberOfPressureFunctions) % numberOfPressureFunctions;
                    assert( pgf !=NULL );
                    int numberOfTimeDerivatives=-1; // this means only compute the pressure 
                    OV_GET_SERIAL_ARRAY(real,pgf[prevPressure][grid],pgfLocal);
                    rectangleExact.evalSolution( tc, cg, grid, pgfLocal, I1,I2,I3, numberOfTimeDerivatives ); 
          // rectangleExact.evalSolution( tc, cg, grid, pgf[prevPressure][grid], I1,I2,I3, numberOfTimeDerivatives ); 
                }
            }
        

    }     

    else if( userKnownSolution=="incompressibleSurfaceWave" )
    {
    // --- incompressible surface wave -----
        const int pc = dbase.get<int >("pc");
        assert( pc>=0 );

        MappedGrid & mg = cg[grid];
        mg.update(MappedGrid::THEvertex | MappedGrid::THEcenter);
        OV_GET_SERIAL_ARRAY_CONST(real,mg.vertex(),vertex);
        RealArray & u = ua;

        int & iswCase = db.get<int>("iswCase");
        int & mPeriod = db.get<int>("mPeriod");  // used in case 1

        if( t<=2*dt )
            printF("UDKS: evaluate the incompressible surface wave: iswCase=%i, t=%9.3e\n",iswCase,t);

    // This next file will define some variables and macros ueSurfWaveCase1, veSurfWaveCase1, peSurfWaveCase1
        Real D1 =1.; // defines the scale
        if( iswCase==1 )
        {

            if( numberOfDimensions==2 )
            {
        // #Include "incompressible/ismStripTractionDirichlet.h"
// ----------------------------------------------------------------------------------------------------
//     SurfWaveCase1 : Exact solution for the incompressible elasticity equations.
//    Periodic strip [0,1]x[0:Ly] :  x=0: traction BC, x=1 : displacement BC.
//    File writtten by exactSolutionsIncompressibleElasticity.mw (maple)
//    Set D1 to scale the solution. Set mPeriod to set periods in y direction.
// ----------------------------------------------------------------------------------------------------
Real rho, c, mu, omega, Pi, k, Ly, LyPeriod1, A, B, C; 
rho = 1;
c = 1;
mu = 1;
omega = 4.9434667681945157224;
Pi=3.1415926535897932385;
k = 4.9434667681945157224;
LyPeriod1 = 1.2710078982637464465;
Ly = 2*Pi*mPeriod/k;
A = -.70285260952300159138e-2*D1;
B = -1.9929714739047699842*D1;
C = -.80346163570342302261*D1;
                int i1,i2,i3;
                FOR_3D(i1,i2,i3,I1,I2,I3)
                {
          // Reference coordinates:
                    Real x = vertex(i1,i2,i3,0);
                    Real y = vertex(i1,i2,i3,1);

                    u(i1,i2,i3,u1c) = -D1*(.80346163570342302253-x+.14217807916603065948e-2*exp(4.9434667681945157224*x)-.40315259864337181788*exp(-4.9434667681945157224*x))*(cos(4.9434667681945157224*t)*cos(omega*y)+sin(4.9434667681945157224*t)*sin(omega*y));;
                    u(i1,i2,i3,u2c) = D1*(-1+.70285260952300159139e-2*exp(4.9434667681945157224*x)+1.9929714739047699841*exp(-4.9434667681945157224*x))*(cos(4.9434667681945157224*t)*sin(omega*y)-sin(4.9434667681945157224*t)*cos(omega*y))/omega;;
                    u(i1,i2,i3,pc ) = -D1*(.70285260952300159138e-2*exp(4.9434667681945157224*x)+1.9929714739047699842*exp(-4.9434667681945157224*x))*(cos(4.9434667681945157224*t)*cos(omega*y)+sin(4.9434667681945157224*t)*sin(omega*y));;
                } 

        // At the initial time we initialize some MOL variables
                    if( computeVelocity )
                    {
                        printF("\n $$$$$$$$$$ userDefinedKnownSolution: EVALUATE VELOCITY, ACCEL and PRESSURE AT PAST TIMES t=%9.3e dt=%16.6e $$$$$$$\n\n",t,dt);
                        assert( dt>0. );
                        assert( vgf!=NULL );
                        const int currentVelocity = dbase.get<int>("currentVelocity");
                        assert( currentVelocity>=0 && currentVelocity<numberOfVelocityFunctions );
            // fill in t0=0 and some past time levels
                        for( int level=0; level<numberOfVelocityFunctions; level++ )
                        {
                            real tc = 0. - level*dt;
                            const int prevVelocity = (currentVelocity - level + numberOfVelocityFunctions) % numberOfVelocityFunctions;
                            OV_GET_SERIAL_ARRAY(real,vgf[prevVelocity][grid],vLocal);
                            OV_GET_SERIAL_ARRAY(real,vtgf[prevVelocity][grid],vtLocal);
                            FOR_3D(i1,i2,i3,I1,I2,I3)
                            {
                                Real x = vertex(i1,i2,i3,0);
                                Real y = vertex(i1,i2,i3,1);
                                vLocal(i1,i2,i3,u1c)  = -4.9434667681945157224*D1*(.80346163570342302253-x+.14217807916603065948e-2*exp(4.9434667681945157224*x)-.40315259864337181788*exp(-4.9434667681945157224*x))*(cos(4.9434667681945157224*tc)*sin(omega*y)-sin(4.9434667681945157224*tc)*cos(omega*y));;
                                vLocal(i1,i2,i3,u2c)  = -4.9434667681945157224*D1*(-1+.70285260952300159139e-2*exp(4.9434667681945157224*x)+1.9929714739047699841*exp(-4.9434667681945157224*x))*(cos(4.9434667681945157224*tc)*cos(omega*y)+sin(4.9434667681945157224*tc)*sin(omega*y))/omega;;
                                vtLocal(i1,i2,i3,u1c) = 24.437863688243529843*D1*(.80346163570342302253-x+.14217807916603065948e-2*exp(4.9434667681945157224*x)-.40315259864337181788*exp(-4.9434667681945157224*x))*(cos(4.9434667681945157224*tc)*cos(omega*y)+sin(4.9434667681945157224*tc)*sin(omega*y));;
                                vtLocal(i1,i2,i3,u2c) = -24.437863688243529843*D1*(-1+.70285260952300159139e-2*exp(4.9434667681945157224*x)+1.9929714739047699841*exp(-4.9434667681945157224*x))*(cos(4.9434667681945157224*tc)*sin(omega*y)-sin(4.9434667681945157224*tc)*cos(omega*y))/omega;;        
                            }         
                        }
            // --- save pressure for time extrapolation ---
                        const int currentPressure = dbase.get<int>("currentPressure");
                        assert( currentPressure>=0 && currentPressure<max(1,numberOfPressureFunctions) );
                        for( int level=0; level<numberOfPressureFunctions; level++ )
                        {
                            real tc = 0. - (level+3)*dt;  // -- save pressure starting at t-3*dt 
                            const int prevPressure = (currentPressure - level + numberOfPressureFunctions) % numberOfPressureFunctions;
                            assert( pgf !=NULL );
                            OV_GET_SERIAL_ARRAY(real,pgf[prevPressure][grid],pLocal);
                            FOR_3D(i1,i2,i3,I1,I2,I3)
                            {
                                Real x = vertex(i1,i2,i3,0);
                                Real y = vertex(i1,i2,i3,1);
                                pLocal(i1,i2,i3) = -D1*(.70285260952300159138e-2*exp(4.9434667681945157224*x)+1.9929714739047699842*exp(-4.9434667681945157224*x))*(cos(4.9434667681945157224*tc)*cos(omega*y)+sin(4.9434667681945157224*tc)*sin(omega*y));; 
                            }
                        }
                    }
            }
            else
            {
        // 3D 
                assert( u1c>=0 && u2c>=0 && u3c>=0 );

                Real u1ex,u2ex,u3ex,pex, v1ex,v2ex,v3ex, a1ex,a2ex,a3ex;

// ----------------------------------------------------------------------------------------------------
//     caseName : Exact solution for the incompressible elasticity equations.
//  Period strip in y and z.
//    File written by ismStri3d.maple
// ----------------------------------------------------------------------------------------------------
Real rho, c, mu, omega; 
rho = 1;
c = 1;
mu = 1;
omega = 8.51405494367735310088807080924;
Real t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49
        ,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,t86,t87,t88,t89,t90,t91,t92,t93,t94,t95,t96,t97,t98,t99
        ,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t130,t131,t132,t133,t134,t135,t136,t137,t138,t139,t140,t141,t142,t143,t144,t145,t146,t147,t148,t149
        ,t150,t151,t152,t153,t154,t155,t156,t157,t158,t159,t160,t161,t162,t163,t164,t165,t166,t167,t168,t169,t170,t171,t172,t173,t174,t175,t176,t177,t178,t179,t180,t181,t182,t183,t184,t185,t186,t187,t188,t189,t190,t191,t192,t193,t194,t195,t196,t197,t198,t199
        ,t200,t201,t202,t203,t204,t205,t206,t207,t208,t209,t210,t211,t212,t213,t214,t215,t216,t217,t218,t219,t220,t221,t222,t223,t224,t225,t226,t227,t228,t229,t230,t231,t232,t233,t234,t235,t236,t237,t238,t239,t240,t241,t242,t243,t244,t245,t246,t247,t248,t249
        ,t250,t251,t252,t253,t254,t255,t256,t257,t258,t259,t260,t261,t262,t263,t264,t265,t266,t267,t268,t269,t270,t271,t272,t273,t274,t275,t276,t277,t278,t279,t280,t281,t282,t283,t284,t285,t286,t287,t288,t289,t290,t291,t292,t293,t294,t295,t296,t297,t298,t299
        ,t300,t301,t302,t303,t304,t305,t306,t307,t308,t309,t310,t311,t312,t313,t314,t315,t316,t317,t318,t319,t320,t321,t322,t323,t324,t325,t326,t327,t328,t329,t330,t331,t332,t333,t334,t335,t336,t337,t338,t339,t340,t341,t342,t343,t344,t345,t346,t347,t348,t349
        ,t350,t351,t352,t353,t354,t355,t356,t357,t358,t359,t360,t361,t362,t363,t364,t365,t366,t367,t368,t369,t370,t371,t372,t373,t374,t375,t376,t377,t378,t379,t380,t381,t382,t383,t384,t385,t386,t387,t388,t389,t390,t391,t392,t393,t394,t395,t396,t397,t398,t399
        ,t400,t401,t402,t403,t404,t405,t406,t407,t408,t409,t410,t411,t412,t413,t414,t415,t416,t417,t418,t419,t420,t421,t422,t423,t424,t425,t426,t427,t428,t429,t430,t431,t432,t433,t434,t435,t436,t437,t438,t439,t440,t441,t442,t443,t444,t445,t446,t447,t448,t449
        ,t450,t451,t452,t453,t454,t455,t456,t457,t458,t459,t460,t461,t462,t463,t464,t465,t466,t467,t468,t469,t470,t471,t472,t473,t474,t475,t476,t477,t478,t479,t480,t481,t482,t483,t484,t485,t486,t487,t488,t489,t490,t491,t492,t493,t494,t495,t496,t497,t498,t499
        ,t500;


                int i1,i2,i3;
                FOR_3D(i1,i2,i3,I1,I2,I3)
                {
          // Reference coordinates:
                    const Real x = vertex(i1,i2,i3,0);
                    const Real y = vertex(i1,i2,i3,1);
                    const Real z = vertex(i1,i2,i3,2);
                    t1 = 0.888576587631673249403176198012e1 * x;
                    t2 = exp(-t1);
                    t4 = exp(t1);
                    t6 = 0.254316802920257171663099173637e1 * x;
                    t7 = exp(t6);
                    t9 = exp(-t6);
                    t11 = -0.119361657176847896968407609330e1 * t2 + 0.194198977888135501364797419613e-4 * t4 - 0.248134004570264724810053450765e-1 * t7 + 0.223126730498585545802558754804e1 * t9;
                    t12 = 0.628318530717958647692528676656e1 * z;
                    t13 = cos(t12);
                    t18 = -0.954571178312983497665745806050e-29 * t4 + 0.2320e-28 * t7 - 0.1270e-26 * t9;
                    t19 = sin(t12);
                    t20 = t19 * t18;
                    t22 = 0.628318530717958647692528676656e1 * y;
                    t23 = cos(t22);
                    t25 = sin(t22);
                    t32 = 0.851405494367735310088807080924e1 * t;
                    t33 = cos(t32);
                    t36 = sin(t32);
                    t44 = t36 * t13;
                    t46 = t36 * t25;
                    u1ex = t33 * (t23 * (t11 * t13 + t20) + t13 * t25 * t18 - t25 * t19 * t11) + t23 * (t11 * t19 * t36 - t13 * t18 * t36) + t44 * t25 * t11 + t46 * t20;
                    t51 = -0.674983753310343661719230574477e-29 * t4 - 0.1110290e-26 * t7 + 0.32880e-26 * t9;
                    t57 = -0.137319414164197011028089002269e-4 * t4 + 0.502170823021904425967534773136e-2 * t7 + 0.451561381466894154031760665905e0 * t9 - 0.844014372034130858066736159642e0 * t2;
                    t58 = t19 * t57;
                    u2ex = t33 * (t23 * (t13 * t51 + t58) + t13 * t25 * t57 - t25 * t19 * t51) + t23 * (-t13 * t36 * t57 + t19 * t36 * t51) + t44 * t25 * t51 + t46 * t58;
                    t81 = -0.674983753310343661719230574480e-29 * t4 + 0.111968038010069201967379861790e-26 * t7 - 0.277395764103970409544292048561e-26 * t9;
                    t87 = -0.137319414164197011028089002268e-4 * t4 + 0.502170823021904425967534865720e-2 * t7 + 0.451561381466894154031760695134e0 * t9 - 0.844014372034130858066736159640e0 * t2;
                    t88 = t19 * t87;
                    u3ex = t33 * (t23 * (t13 * t81 + t88) + t13 * t25 * t87 - t25 * t19 * t81) + t23 * (-t13 * t36 * t87 + t19 * t36 * t81) + t44 * t25 * t81 + t46 * t88;
                    t112 = t19 * t25;
                    t114 = t25 * t13;
                    t118 = t36 * t19;
                    t129 = t2 * t13;
                    pex = t4 * (t33 * (t23 * (0.158425457721357387311960260192e-3 * t13 - 0.7787290e-28 * t19) - 0.158425457721357387311960260192e-3 * t112 - 0.7787290e-28 * t114) + t23 * (0.158425457721357387311960260192e-3 * t118 + 0.7787290e-28 * t44) + 0.158425457721357387311960260192e-3 * t36 * t114 - 0.7787290e-28 * t36 * t112) + t33 * (0.973739685875924704361563447422e1 * t23 * t129 - 0.973739685875924704361563447422e1 * t19 * t25 * t2) + 0.973739685875924704361563447422e1 * t118 * t23 * t2 + 0.973739685875924704361563447422e1 * t46 * t129;
                    t145 = 0.812727145980757329702715098298e-28 * t4 - 0.197526074693314591940603242774e-27 * t7 + 0.108128497784702384381278499277e-25 * t9;
                    t151 = 0.165342076774556904631055919516e-3 * t4 - 0.211262654830592130851901386399e0 * t7 + 0.189971324286804670354286206957e2 * t9 - 0.101625170737206325095890229578e2 * t2;
                    t152 = t19 * t151;
                    v1ex = t33 * (t23 * (t13 * t145 + t152) + t13 * t25 * t151 - t25 * t19 * t145) + t23 * (-t13 * t151 * t36 + t145 * t19 * t36) + t44 * t25 * t145 + t46 * t152;
                    t176 = 0.718598473675192854892334279278e1 * t2 + 0.116914503702755950636577879930e-3 * t4 - 0.427550997832017053933186857275e-1 * t7 - 0.384461841225198526496981706061e1 * t9;
                    t181 = -0.574684876177382640444316486014e-28 * t4 - 0.945307006341552837438501613879e-26 * t7 + 0.279942126548111369957199768208e-25 * t9;
                    t182 = t19 * t181;
                    v2ex = t33 * (t23 * (t13 * t176 + t182) + t13 * t25 * t181 - t25 * t19 * t176) + t23 * (-t13 * t36 * t181 + t36 * t19 * t176) + t44 * t25 * t176 + t46 * t182;
                    t206 = 0.718598473675192854892334279272e1 * t2 + 0.116914503702755950636577879929e-3 * t4 - 0.427550997832017053933186936106e-1 * t7 - 0.384461841225198526496981730946e1 * t9;
                    t211 = -0.574684876177382640444316486018e-28 * t4 + 0.953302027553483470528909588326e-26 * t7 - 0.236176277672456611229535425342e-25 * t9;
                    t212 = t19 * t211;
                    v3ex = t33 * (t23 * (t13 * t206 + t212) + t13 * t25 * t211 - t25 * t19 * t206) + t23 * (-t13 * t211 * t36 + t19 * t206 * t36) + t44 * t25 * t206 + t46 * t212;
                    t236 = 0.865242287317166590720655495383e2 * t2 - 0.140773152616029667891934083527e-2 * t4 + 0.179870185077480517364788155132e1 * t7 - 0.161742629270100291889087646183e3 * t9;
                    t241 = 0.691960357509825277834485997503e-27 * t4 - 0.168174785274779720976568633067e-26 * t7 + 0.920611971116251058794147258600e-25 * t9;
                    t242 = t19 * t241;
                    a1ex = t33 * (t23 * (t13 * t236 + t242) + t13 * t25 * t241 - t25 * t19 * t236) + t23 * (-t13 * t241 * t36 + t19 * t236 * t36) + t44 * t25 * t236 + t46 * t242;
                    t265 = 0.489289861107465219758955163844e-27 * t4 + 0.804839579063513691392562015552e-25 * t7 - 0.238344264648049880418516235140e-24 * t9;
                    t271 = 0.995416508238033505861871474849e-3 * t4 - 0.364019268676587007586945274530e0 * t7 - 0.327332923993869911100496892757e2 * t9 + 0.611818688731327600711168700469e2 * t2;
                    t272 = t19 * t271;
                    a2ex = t33 * (t23 * (t13 * t265 + t272) + t13 * t25 * t271 - t25 * t19 * t265) + t23 * (-t13 * t271 * t36 + t19 * t265 * t36) + t44 * t25 * t265 + t46 * t272;
                    t295 = 0.489289861107465219758955163846e-27 * t4 - 0.811646584050938022357169465252e-25 * t7 + 0.201081780449649447982949592399e-24 * t9;
                    t301 = 0.995416508238033505861871474838e-3 * t4 - 0.364019268676587007586945341646e0 * t7 - 0.327332923993869911100496913944e2 * t9 + 0.611818688731327600711168700471e2 * t2;
                    t302 = t19 * t301;
                    a3ex = t33 * (t23 * (t13 * t295 + t302) + t13 * t25 * t301 - t25 * t19 * t295) + t23 * (-t13 * t301 * t36 + t19 * t295 * t36) + t44 * t25 * t295 + t46 * t302;

                    u(i1,i2,i3,u1c) = u1ex;
                    u(i1,i2,i3,u2c) = u2ex;
                    u(i1,i2,i3,u3c) = u3ex;
                    u(i1,i2,i3,pc ) = pex;
                }        
                if( computeVelocity )
                {
                    OV_ABORT("UDKS: finish me -- computeVelocity");
                }

            }
        }
        else if( iswCase==2 )
        {
            if( numberOfDimensions==2 )
            {
// ----------------------------------------------------------------------------------------------------
//     SurfWaveTractionTractionCase2 : Exact solution for the incompressible elasticity equations.
//    Periodic strip [0,1]x[0:Ly] :  x=0: traction BC, x=1 : traction BC.
//    File writtten by exactSolutionsIncompressibleElasticity.mw (maple)
//    Set D1 to scale the solution. Set mPeriod to set periods in y direction.
// ----------------------------------------------------------------------------------------------------
Real rho, c, mu, omega, Pi, k, Ly, LyPeriod1, A, B, C; 
rho = 1;
c = 1;
mu = 1;
omega = 5.7766060906888207001;
Pi=3.1415926535897932385;
k = 2*Pi;
Ly = mPeriod;
A = -.61982149405896696448e-2*D1;
B = 3.3190923798346935004*D1;
C = .84442232515900280636e-1*D1;
                int i1,i2,i3;
                FOR_3D(i1,i2,i3,I1,I2,I3)
                {
          // Reference coordinates:
                    Real x = vertex(i1,i2,i3,0);
                    Real y = vertex(i1,i2,i3,1);

                    u(i1,i2,i3,u1c) = D1*(cos(5.7766060906888207001*t)*cos(2*Pi*y)+sin(5.7766060906888207001*t)*sin(2*Pi*y))*(.19132481438887511774e-2*exp(2.4716876172717001887*x-6.283185307179586477)+.1580109802013600538e-3*exp(2.4716876172717001887*x+6.283185307179586477)-.17490174455709712477e-3*exp(2.4716876172717001887*x)+exp(-2.4716876172717001887*x)-.11958920445719399496e-2*exp(2*Pi*(x-1))-.11648477251471055807e-2*exp(2*Pi*x)-.11958920445719399496e-2*exp(-2*Pi*(x-1))+.15428081958872228858e-1*exp(-2*Pi*x));;
                    u(i1,i2,i3,u2c) = -.11958920445719399496e-2*D1*(cos(5.7766060906888207001*t)*sin(2*Pi*y)-1.*sin(5.7766060906888207001*t)*cos(2*Pi*y))*(.62935116044792865932*exp(2.4716876172717001887*x-6.283185307179586477)+.51976736039642007170e-1*exp(2.4716876172717001887*x+6.283185307179586477)-.57532848654772680303e-1*exp(2.4716876172717001887*x)-328.94382386215097900*exp(-2.4716876172717001887*x)-1.0000000000000000000*exp(2*Pi*(x-1))-.97404086801501682686*exp(2*Pi*x)+1.0000000000000000000*exp(-2*Pi*(x-1))-12.900898562624511937*exp(-2*Pi*x));;
                    u(i1,i2,i3,pc ) = D1*(cos(5.7766060906888207001*t)*cos(2*Pi*y)+sin(5.7766060906888207001*t)*sin(2*Pi*y))*(-.6351226721132286536e-2*exp(2*Pi*(x-1))-.61863543884118615911e-2*exp(2*Pi*x)+.6351226721132286536e-2*exp(-2*Pi*(x-1))-.81936531677557907269e-1*exp(-2*Pi*x));;
                }  

        // At the initial time we initialize some MOL variables
                    if( computeVelocity )
                    {
                        printF("\n $$$$$$$$$$ userDefinedKnownSolution: EVALUATE VELOCITY, ACCEL and PRESSURE AT PAST TIMES t=%9.3e dt=%16.6e $$$$$$$\n\n",t,dt);
                        assert( dt>0. );
                        assert( vgf!=NULL );
                        const int currentVelocity = dbase.get<int>("currentVelocity");
                        assert( currentVelocity>=0 && currentVelocity<numberOfVelocityFunctions );
            // fill in t0=0 and some past time levels
                        for( int level=0; level<numberOfVelocityFunctions; level++ )
                        {
                            real tc = 0. - level*dt;
                            const int prevVelocity = (currentVelocity - level + numberOfVelocityFunctions) % numberOfVelocityFunctions;
                            OV_GET_SERIAL_ARRAY(real,vgf[prevVelocity][grid],vLocal);
                            OV_GET_SERIAL_ARRAY(real,vtgf[prevVelocity][grid],vtLocal);
                            FOR_3D(i1,i2,i3,I1,I2,I3)
                            {
                                Real x = vertex(i1,i2,i3,0);
                                Real y = vertex(i1,i2,i3,1);
                                vLocal(i1,i2,i3,u1c)  = -.69081972684805749514e-2*D1*(cos(5.7766060906888207001*tc)*sin(2*Pi*y)-1.*sin(5.7766060906888207001*tc)*cos(2*Pi*y))*(-1.5998502143841780809*exp(2.4716876172717001887*x-6.283185307179586477)-.13212813056041260650*exp(2.4716876172717001887*x+6.283185307179586477)+.14625211811631526809*exp(2.4716876172717001887*x)-836.19587950176727728*exp(-2.4716876172717001887*x)+1.0000000000000000000*exp(2*Pi*(x-1))+.97404086801501682685*exp(2*Pi*x)+1.0000000000000000000*exp(-2*Pi*(x-1))-12.900898562624511937*exp(-2*Pi*x));;
                                vLocal(i1,i2,i3,u2c)  = D1*(cos(5.7766060906888207001*tc)*cos(2*Pi*y)+sin(5.7766060906888207001*tc)*sin(2*Pi*y))*(.43476819675214608236e-2*exp(2.4716876172717001887*x-6.283185307179586477)+.359065545933590771e-3*exp(2.4716876172717001887*x+6.283185307179586477)-.39744826792480695115e-3*exp(2.4716876172717001887*x)-2.2724088254880667638*exp(-2.4716876172717001887*x)-.69081972684805749512e-2*exp(2*Pi*(x-1))-.67288664638097874686e-2*exp(2*Pi*x)+.69081972684805749512e-2*exp(-2*Pi*(x-1))-.89121952211267628966e-1*exp(-2*Pi*x));;
                                vtLocal(i1,i2,i3,u1c) = .39905934416784763590e-1*D1*(cos(5.7766060906888207001*tc)*cos(2*Pi*y)+sin(5.7766060906888207001*tc)*sin(2*Pi*y))*(-1.5998502143841780809*exp(2.4716876172717001887*x-6.283185307179586477)-.13212813056041260650*exp(2.4716876172717001887*x+6.283185307179586477)+.14625211811631526809*exp(2.4716876172717001887*x)-836.19587950176727727*exp(-2.4716876172717001887*x)+1.0000000000000000000*exp(2*Pi*(x-1))+.97404086801501682683*exp(2*Pi*x)+1.0000000000000000000*exp(-2*Pi*(x-1))-12.900898562624511936*exp(-2*Pi*x));;
                                vtLocal(i1,i2,i3,u2c) = D1*(cos(5.7766060906888207001*tc)*sin(2*Pi*y)-sin(5.7766060906888207001*tc)*cos(2*Pi*y))*(.25114846133962426136e-1*exp(2.4716876172717001887*x-6.283185307179586477)+.20741802195964869641e-2*exp(2.4716876172717001887*x+6.283185307179586477)-.22959020852281620903e-2*exp(2.4716876172717001887*x)-13.126810661849395928*exp(-2.4716876172717001887*x)-.39905934416784763589e-1*exp(2*Pi*(x-1))-.38870010998275365402e-1*exp(2*Pi*x)+.39905934416784763589e-1*exp(-2*Pi*(x-1))-.51482241195768659763*exp(-2*Pi*x));;        
                            }         
                        }
            // --- save pressure for time extrapolation ---
                        const int currentPressure = dbase.get<int>("currentPressure");
                        assert( currentPressure>=0 && currentPressure<max(1,numberOfPressureFunctions) );
                        for( int level=0; level<numberOfPressureFunctions; level++ )
                        {
                            real tc = 0. - (level+3)*dt;  // -- save pressure starting at t-3*dt 
                            const int prevPressure = (currentPressure - level + numberOfPressureFunctions) % numberOfPressureFunctions;
                            assert( pgf !=NULL );
                            OV_GET_SERIAL_ARRAY(real,pgf[prevPressure][grid],pLocal);
                            FOR_3D(i1,i2,i3,I1,I2,I3)
                            {
                                Real x = vertex(i1,i2,i3,0);
                                Real y = vertex(i1,i2,i3,1);
                                pLocal(i1,i2,i3) = D1*(cos(5.7766060906888207001*tc)*cos(2*Pi*y)+sin(5.7766060906888207001*tc)*sin(2*Pi*y))*(-.6351226721132286536e-2*exp(2*Pi*(x-1))-.61863543884118615911e-2*exp(2*Pi*x)+.6351226721132286536e-2*exp(-2*Pi*(x-1))-.81936531677557907269e-1*exp(-2*Pi*x));; 
                            }
                        }
                    }
            }
            else
            {
        // 3D 
                assert( u1c>=0 && u2c>=0 && u3c>=0 );

                Real u1ex,u2ex,u3ex,pex, v1ex,v2ex,v3ex, a1ex,a2ex,a3ex;

// ----------------------------------------------------------------------------------------------------
//     caseName : Exact solution for the incompressible elasticity equations.
//  Period strip in y and z.
//    File written by ismStri3d.maple
// ----------------------------------------------------------------------------------------------------
Real rho, c, mu, omega; 
rho = 1;
c = 1;
mu = 1;
omega = 10.7649676392784327821062965544;
Real t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49
        ,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,t86,t87,t88,t89,t90,t91,t92,t93,t94,t95,t96,t97,t98,t99
        ,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t130,t131,t132,t133,t134,t135,t136,t137,t138,t139,t140,t141,t142,t143,t144,t145,t146,t147,t148,t149
        ,t150,t151,t152,t153,t154,t155,t156,t157,t158,t159,t160,t161,t162,t163,t164,t165,t166,t167,t168,t169,t170,t171,t172,t173,t174,t175,t176,t177,t178,t179,t180,t181,t182,t183,t184,t185,t186,t187,t188,t189,t190,t191,t192,t193,t194,t195,t196,t197,t198,t199
        ,t200,t201,t202,t203,t204,t205,t206,t207,t208,t209,t210,t211,t212,t213,t214,t215,t216,t217,t218,t219,t220,t221,t222,t223,t224,t225,t226,t227,t228,t229,t230,t231,t232,t233,t234,t235,t236,t237,t238,t239,t240,t241,t242,t243,t244,t245,t246,t247,t248,t249
        ,t250,t251,t252,t253,t254,t255,t256,t257,t258,t259,t260,t261,t262,t263,t264,t265,t266,t267,t268,t269,t270,t271,t272,t273,t274,t275,t276,t277,t278,t279,t280,t281,t282,t283,t284,t285,t286,t287,t288,t289,t290,t291,t292,t293,t294,t295,t296,t297,t298,t299
        ,t300,t301,t302,t303,t304,t305,t306,t307,t308,t309,t310,t311,t312,t313,t314,t315,t316,t317,t318,t319,t320,t321,t322,t323,t324,t325,t326,t327,t328,t329,t330,t331,t332,t333,t334,t335,t336,t337,t338,t339,t340,t341,t342,t343,t344,t345,t346,t347,t348,t349
        ,t350,t351,t352,t353,t354,t355,t356,t357,t358,t359,t360,t361,t362,t363,t364,t365,t366,t367,t368,t369,t370,t371,t372,t373,t374,t375,t376,t377,t378,t379,t380,t381,t382,t383,t384,t385,t386,t387,t388,t389,t390,t391,t392,t393,t394,t395,t396,t397,t398,t399
        ,t400,t401,t402,t403,t404,t405,t406,t407,t408,t409,t410,t411,t412,t413,t414,t415,t416,t417,t418,t419,t420,t421,t422,t423,t424,t425,t426,t427,t428,t429,t430,t431,t432,t433,t434,t435,t436,t437,t438,t439,t440,t441,t442,t443,t444,t445,t446,t447,t448,t449
        ,t450,t451,t452,t453,t454,t455,t456,t457,t458,t459,t460,t461,t462,t463,t464,t465,t466,t467,t468,t469,t470,t471,t472,t473,t474,t475,t476,t477,t478,t479,t480,t481,t482,t483,t484,t485,t486,t487,t488,t489,t490,t491,t492,t493,t494,t495,t496,t497,t498,t499
        ,t500;


                int i1,i2,i3;
                FOR_3D(i1,i2,i3,I1,I2,I3)
                {
          // Reference coordinates:
                    const Real x = vertex(i1,i2,i3,0);
                    const Real y = vertex(i1,i2,i3,1);
                    const Real z = vertex(i1,i2,i3,2);
                    t1 = 0.888576587631673249403176198012e1 * x;
                    t2 = exp(t1);
                    t4 = 0.607681603029061635492190309752e1 * x;
                    t5 = cos(t4);
                    t7 = sin(t4);
                    t9 = exp(-t1);
                    t11 = -0.753239723406067534184872382763e-4 * t2 + 0.204598093767530818350584185509e1 * t5 - 0.211866255833191900563999356400e0 * t7 - 0.544467928195366972426789429740e0 * t9;
                    t12 = 0.628318530717958647692528676656e1 * z;
                    t13 = cos(t12);
                    t19 = 0.2193926127766250949726433469e-1 * t7 - 0.211866255833191900563999359800e0 * t5 + 0.779997882698741460167705417491e-5 * t2 + 0.563809658456915764794509377163e-1 * t9;
                    t20 = sin(t12);
                    t21 = t20 * t19;
                    t23 = 0.628318530717958647692528676656e1 * y;
                    t24 = cos(t23);
                    t26 = sin(t23);
                    t33 = 0.107649676392784327821062965544e2 * t;
                    t34 = cos(t33);
                    t37 = sin(t33);
                    t45 = t13 * t37;
                    t47 = t37 * t26;
                    u1ex = t34 * (t24 * (t13 * t11 + t21) + t13 * t26 * t19 - t26 * t20 * t11) + t24 * (t37 * t20 * t11 - t13 * t37 * t19) + t45 * t26 * t11 + t47 * t21;
                    t51 = 0.102453787114446688167969810662e0 * t7;
                    t53 = 0.551541792167429340705776308962e-5 * t2 + 0.10609336515389131435859593393e-1 * t5 + t51 - 0.398673632793356429650421871652e-1 * t9;
                    t59 = 0.989390663484610868564140406607e0 * t7 + 0.102453787114446688167969810662e0 * t5 + 0.532620916279509775647851210991e-4 * t2 - 0.384996964165534223106925228171e0 * t9;
                    t60 = t20 * t59;
                    u2ex = t34 * (t24 * (t13 * t53 + t60) + t13 * t26 * t59 - t26 * t20 * t53) + t24 * (-t13 * t37 * t59 + t37 * t20 * t53) + t45 * t26 * t53 + t47 * t60;
                    t83 = 0.551541792167429340705776308960e-5 * t2 + 0.106093365153891314358595930604e-1 * t5 + t51 - 0.398673632793356429650421871650e-1 * t9;
                    t89 = 0.989390663484610868564140406600e0 * t7 + 0.102453787114446688167969807374e0 * t5 + 0.532620916279509775647851210990e-4 * t2 - 0.384996964165534223106925228170e0 * t9;
                    t90 = t20 * t89;
                    u3ex = t34 * (t24 * (t13 * t83 + t90) + t13 * t26 * t89 - t26 * t20 * t83) + t24 * (-t13 * t37 * t89 + t37 * t20 * t83) + t45 * t26 * t83 + t47 * t90;
                    t114 = t20 * t2;
                    t116 = t9 * t20;
                    t120 = t26 * t2;
                    t122 = t26 * t9;
                    t134 = t9 * t37;
                    t150 = t37 * t20;
                    pex = t34 * (t24 * (t13 * (-0.982344473618620248763243082755e-3 * t2 + 0.710072827687236792905290602940e1 * t9) + 0.101724137176215384777217667547e-3 * t114 - 0.735297522086965837907740907410e0 * t116) + t13 * (0.101724137176215384777217667547e-3 * t120 - 0.735297522086965837907740907410e0 * t122) + 0.982344473618620248763243082755e-3 * t20 * t120 - 0.710072827687236792905290602940e1 * t26 * t116) + t24 * (t13 * (-0.101724137176215384777217667547e-3 * t37 * t2 + 0.735297522086965837907740907410e0 * t134) - 0.982344473618620248763243082755e-3 * t37 * t114 + 0.710072827687236792905290602940e1 * t37 * t116) + t13 * (-0.982344473618620248763243082755e-3 * t37 * t120 + 0.710072827687236792905290602940e1 * t26 * t134) + 0.101724137176215384777217667547e-3 * t150 * t120 - 0.735297522086965837907740907410e0 * t150 * t122;
                    t159 = -0.839665196595764678524978767219e-4 * t2 + 0.228073338789939630269443895847e1 * t5 - 0.236175437683711317856477236823e0 * t7 - 0.606939272800132397554438883761e0 * t9;
                    t165 = -0.228073338789939630269443892187e1 * t7 + 0.220249185846552366494056465318e2 * t5 - 0.810860124708535449352960613025e-3 * t2 - 0.586117962764809884791888083821e1 * t9;
                    t166 = t20 * t165;
                    v1ex = t34 * (t24 * (t13 * t159 + t166) + t13 * t26 * t165 - t26 * t20 * t159) + t24 * (-t13 * t37 * t165 + t37 * t20 * t159) + t45 * t26 * t159 + t47 * t166;
                    t190 = -0.573364692775175013716097731084e-3 * t2 - 0.110291170280854028050252083644e1 * t5 - 0.106507584750160537696246687098e2 * t7 + 0.414447986046241432672814397105e1 * t9;
                    t192 = 0.110291170280854028050252083644e1 * t7;
                    t196 = t192 + 0.114209164262379012481627844618e0 * t5 + 0.593732954439200783815785753958e-4 * t2 - 0.429170875565405494811044982635e0 * t9;
                    t197 = t20 * t196;
                    v2ex = t34 * (t24 * (t13 * t190 + t197) + t13 * t26 * t196 - t26 * t20 * t190) + t24 * (-t13 * t37 * t196 + t37 * t20 * t190) + t45 * t26 * t190 + t47 * t197;
                    t221 = -0.573364692775175013716097731086e-3 * t2 - 0.110291170280854028050252080105e1 * t5 - 0.106507584750160537696246687097e2 * t7 + 0.414447986046241432672814397104e1 * t9;
                    t226 = t192 + 0.114209164262379012481627841038e0 * t5 + 0.593732954439200783815785753951e-4 * t2 - 0.429170875565405494811044982636e0 * t9;
                    t227 = t20 * t226;
                    v3ex = t34 * (t24 * (t13 * t221 + t227) + t13 * t26 * t226 - t26 * t20 * t221) + t24 * (-t13 * t37 * t226 + t37 * t20 * t221) + t45 * t26 * t221 + t47 * t227;
                    t251 = 0.872888300246865845979060058975e-2 * t2 - 0.237097535821555763859928806507e3 * t5 + 0.245520211145588663285774316485e2 * t7 + 0.630954090196297983278707421990e2 * t9;
                    t257 = -0.254242094385757243832246831452e1 * t7 + 0.245520211145588663285774320425e2 * t5 - 0.903896866918177004550662388841e-3 * t2 - 0.653368163070060996488433067136e1 * t9;
                    t258 = t20 * t257;
                    a1ex = t34 * (t24 * (t13 * t251 + t258) + t13 * t26 * t257 - t26 * t20 * t251) + t24 * (-t13 * t37 * t257 + t37 * t20 * t251) + t45 * t26 * t251 + t47 * t258;
                    t280 = 0.118728087897154083062681942597e2 * t7;
                    t282 = -0.639151604091117254913306627012e-3 * t2 - 0.122945795739354494986511003638e1 * t5 - t280 + 0.462001058718238122043635271806e1 * t9;
                    t288 = -0.114655070317318329149234456111e3 * t7 - 0.118728087897154083062681942597e2 * t5 - 0.617225236322957965178608717936e-2 * t2 + 0.446151915795190848609705083554e2 * t9;
                    t289 = t20 * t288;
                    a2ex = t34 * (t24 * (t13 * t282 + t289) + t13 * t26 * t288 - t26 * t20 * t282) + t24 * (-t13 * t37 * t288 + t37 * t20 * t282) + t45 * t26 * t282 + t47 * t289;
                    t312 = -0.639151604091117254913306627010e-3 * t2 - 0.122945795739354494986510999784e1 * t5 - t280 + 0.462001058718238122043635271804e1 * t9;
                    t318 = -0.114655070317318329149234456110e3 * t7 - 0.118728087897154083062681938787e2 * t5 - 0.617225236322957965178608717933e-2 * t2 + 0.446151915795190848609705083552e2 * t9;
                    t319 = t20 * t318;
                    a3ex = t34 * (t24 * (t13 * t312 + t319) + t13 * t26 * t318 - t26 * t20 * t312) + t24 * (-t13 * t37 * t318 + t37 * t20 * t312) + t45 * t26 * t312 + t47 * t319;

                    u(i1,i2,i3,u1c) = u1ex;
                    u(i1,i2,i3,u2c) = u2ex;
                    u(i1,i2,i3,u3c) = u3ex;
                    u(i1,i2,i3,pc ) = pex;
                }          
                if( computeVelocity )
                {
                    OV_ABORT("UDKS: finish me -- computeVelocity");
                }

            }              
        }
        else if( iswCase==3 )
        {
      // BC = DD
            if( numberOfDimensions==2 )
            {
                OV_ABORT("UDKS: finish me : strip BC=DD, 2D");
        // #Include "incompressible/ismStripDisplacementDisplacementCase3k2pi.h"
                int i1,i2,i3;
                FOR_3D(i1,i2,i3,I1,I2,I3)
                {
          // Reference coordinates:
                    Real x = vertex(i1,i2,i3,0);
                    Real y = vertex(i1,i2,i3,1);

                    u(i1,i2,i3,u1c) = D1*(cos(5.7766060906888207001*t)*cos(2*Pi*y)+sin(5.7766060906888207001*t)*sin(2*Pi*y))*(.19132481438887511774e-2*exp(2.4716876172717001887*x-6.283185307179586477)+.1580109802013600538e-3*exp(2.4716876172717001887*x+6.283185307179586477)-.17490174455709712477e-3*exp(2.4716876172717001887*x)+exp(-2.4716876172717001887*x)-.11958920445719399496e-2*exp(2*Pi*(x-1))-.11648477251471055807e-2*exp(2*Pi*x)-.11958920445719399496e-2*exp(-2*Pi*(x-1))+.15428081958872228858e-1*exp(-2*Pi*x));;
                    u(i1,i2,i3,u2c) = -.11958920445719399496e-2*D1*(cos(5.7766060906888207001*t)*sin(2*Pi*y)-1.*sin(5.7766060906888207001*t)*cos(2*Pi*y))*(.62935116044792865932*exp(2.4716876172717001887*x-6.283185307179586477)+.51976736039642007170e-1*exp(2.4716876172717001887*x+6.283185307179586477)-.57532848654772680303e-1*exp(2.4716876172717001887*x)-328.94382386215097900*exp(-2.4716876172717001887*x)-1.0000000000000000000*exp(2*Pi*(x-1))-.97404086801501682686*exp(2*Pi*x)+1.0000000000000000000*exp(-2*Pi*(x-1))-12.900898562624511937*exp(-2*Pi*x));;
                    u(i1,i2,i3,pc ) = D1*(cos(5.7766060906888207001*t)*cos(2*Pi*y)+sin(5.7766060906888207001*t)*sin(2*Pi*y))*(-.6351226721132286536e-2*exp(2*Pi*(x-1))-.61863543884118615911e-2*exp(2*Pi*x)+.6351226721132286536e-2*exp(-2*Pi*(x-1))-.81936531677557907269e-1*exp(-2*Pi*x));;
                }  

        // At the initial time we initialize some MOL variables
                    if( computeVelocity )
                    {
                        printF("\n $$$$$$$$$$ userDefinedKnownSolution: EVALUATE VELOCITY, ACCEL and PRESSURE AT PAST TIMES t=%9.3e dt=%16.6e $$$$$$$\n\n",t,dt);
                        assert( dt>0. );
                        assert( vgf!=NULL );
                        const int currentVelocity = dbase.get<int>("currentVelocity");
                        assert( currentVelocity>=0 && currentVelocity<numberOfVelocityFunctions );
            // fill in t0=0 and some past time levels
                        for( int level=0; level<numberOfVelocityFunctions; level++ )
                        {
                            real tc = 0. - level*dt;
                            const int prevVelocity = (currentVelocity - level + numberOfVelocityFunctions) % numberOfVelocityFunctions;
                            OV_GET_SERIAL_ARRAY(real,vgf[prevVelocity][grid],vLocal);
                            OV_GET_SERIAL_ARRAY(real,vtgf[prevVelocity][grid],vtLocal);
                            FOR_3D(i1,i2,i3,I1,I2,I3)
                            {
                                Real x = vertex(i1,i2,i3,0);
                                Real y = vertex(i1,i2,i3,1);
                                vLocal(i1,i2,i3,u1c)  = -.69081972684805749514e-2*D1*(cos(5.7766060906888207001*tc)*sin(2*Pi*y)-1.*sin(5.7766060906888207001*tc)*cos(2*Pi*y))*(-1.5998502143841780809*exp(2.4716876172717001887*x-6.283185307179586477)-.13212813056041260650*exp(2.4716876172717001887*x+6.283185307179586477)+.14625211811631526809*exp(2.4716876172717001887*x)-836.19587950176727728*exp(-2.4716876172717001887*x)+1.0000000000000000000*exp(2*Pi*(x-1))+.97404086801501682685*exp(2*Pi*x)+1.0000000000000000000*exp(-2*Pi*(x-1))-12.900898562624511937*exp(-2*Pi*x));;
                                vLocal(i1,i2,i3,u2c)  = D1*(cos(5.7766060906888207001*tc)*cos(2*Pi*y)+sin(5.7766060906888207001*tc)*sin(2*Pi*y))*(.43476819675214608236e-2*exp(2.4716876172717001887*x-6.283185307179586477)+.359065545933590771e-3*exp(2.4716876172717001887*x+6.283185307179586477)-.39744826792480695115e-3*exp(2.4716876172717001887*x)-2.2724088254880667638*exp(-2.4716876172717001887*x)-.69081972684805749512e-2*exp(2*Pi*(x-1))-.67288664638097874686e-2*exp(2*Pi*x)+.69081972684805749512e-2*exp(-2*Pi*(x-1))-.89121952211267628966e-1*exp(-2*Pi*x));;
                                vtLocal(i1,i2,i3,u1c) = .39905934416784763590e-1*D1*(cos(5.7766060906888207001*tc)*cos(2*Pi*y)+sin(5.7766060906888207001*tc)*sin(2*Pi*y))*(-1.5998502143841780809*exp(2.4716876172717001887*x-6.283185307179586477)-.13212813056041260650*exp(2.4716876172717001887*x+6.283185307179586477)+.14625211811631526809*exp(2.4716876172717001887*x)-836.19587950176727727*exp(-2.4716876172717001887*x)+1.0000000000000000000*exp(2*Pi*(x-1))+.97404086801501682683*exp(2*Pi*x)+1.0000000000000000000*exp(-2*Pi*(x-1))-12.900898562624511936*exp(-2*Pi*x));;
                                vtLocal(i1,i2,i3,u2c) = D1*(cos(5.7766060906888207001*tc)*sin(2*Pi*y)-sin(5.7766060906888207001*tc)*cos(2*Pi*y))*(.25114846133962426136e-1*exp(2.4716876172717001887*x-6.283185307179586477)+.20741802195964869641e-2*exp(2.4716876172717001887*x+6.283185307179586477)-.22959020852281620903e-2*exp(2.4716876172717001887*x)-13.126810661849395928*exp(-2.4716876172717001887*x)-.39905934416784763589e-1*exp(2*Pi*(x-1))-.38870010998275365402e-1*exp(2*Pi*x)+.39905934416784763589e-1*exp(-2*Pi*(x-1))-.51482241195768659763*exp(-2*Pi*x));;        
                            }         
                        }
            // --- save pressure for time extrapolation ---
                        const int currentPressure = dbase.get<int>("currentPressure");
                        assert( currentPressure>=0 && currentPressure<max(1,numberOfPressureFunctions) );
                        for( int level=0; level<numberOfPressureFunctions; level++ )
                        {
                            real tc = 0. - (level+3)*dt;  // -- save pressure starting at t-3*dt 
                            const int prevPressure = (currentPressure - level + numberOfPressureFunctions) % numberOfPressureFunctions;
                            assert( pgf !=NULL );
                            OV_GET_SERIAL_ARRAY(real,pgf[prevPressure][grid],pLocal);
                            FOR_3D(i1,i2,i3,I1,I2,I3)
                            {
                                Real x = vertex(i1,i2,i3,0);
                                Real y = vertex(i1,i2,i3,1);
                                pLocal(i1,i2,i3) = D1*(cos(5.7766060906888207001*tc)*cos(2*Pi*y)+sin(5.7766060906888207001*tc)*sin(2*Pi*y))*(-.6351226721132286536e-2*exp(2*Pi*(x-1))-.61863543884118615911e-2*exp(2*Pi*x)+.6351226721132286536e-2*exp(-2*Pi*(x-1))-.81936531677557907269e-1*exp(-2*Pi*x));; 
                            }
                        }
                    }
            }
            else
            {
                OV_ABORT("UDKS: finish me : strip BC=DD, 3D");

        // // 3D 
        // assert( u1c>=0 && u2c>=0 && u3c>=0 );

        // Real u1ex,u2ex,u3ex,pex, v1ex,v2ex,v3ex, a1ex,a2ex,a3ex;

        // #Include "incompressible/ismStrip3dTractionTractionRoot1.h"


        // int i1,i2,i3;
        // FOR_3D(i1,i2,i3,I1,I2,I3)
        // {
        //   // Reference coordinates:
        //   const Real x = vertex(i1,i2,i3,0);
        //   const Real y = vertex(i1,i2,i3,1);
        //   const Real z = vertex(i1,i2,i3,2);
        //   evalSolution(x,y,z,t, u1ex,u2ex,u3ex,pex, v1ex,v2ex,v3ex, a1ex,a2ex,a3ex);

        //   u(i1,i2,i3,u1c) = u1ex;
        //   u(i1,i2,i3,u2c) = u2ex;
        //   u(i1,i2,i3,u3c) = u3ex;
        //   u(i1,i2,i3,pc ) = pex;
        // }          
        // if( computeVelocity )
        // {
        //   OV_ABORT("UDKS: finish me -- computeVelocity");
        // }

            }              
        }
        else if( iswCase==5 )
        {
            
      // -- traction - displacement strip 
            if( numberOfDimensions==2 )
            {
// ----------------------------------------------------------------------------------------------------
//     SurfWaveDisplacementDisplacementCase5 : Exact solution for the incompressible elasticity equations.
//    Periodic strip [0,1]x[0:Ly] :  x=0: traction BC, x=1 : displacement BC.
//    File writtten by exactSolutionsIncompressibleElasticity.mw (maple)
//    Set D1 to scale the solution. Set mPeriod to set periods in y direction.
// ----------------------------------------------------------------------------------------------------
Real rho, c, mu, omega, Pi, k, Ly, LyPeriod1, A, B, C; 
rho = 1;
c = 1;
mu = 1;
omega = 6.1046733125517700372;
Pi=3.1415926535897932385;
k = 2*Pi;
Ly = mPeriod;
A = .15362742382885997805e-2*D1;
B = 2.8757257890665898734*D1;
C = -.82237422932641953477e-1*D1;
                int i1,i2,i3;
                FOR_3D(i1,i2,i3,I1,I2,I3)
                {
          // Reference coordinates:
                    Real x = vertex(i1,i2,i3,0);
                    Real y = vertex(i1,i2,i3,1);

                    u(i1,i2,i3,u1c) = -D1*(cos(6.1046733125517700372*t)*cos(2*Pi*y)+sin(6.1046733125517700372*t)*sin(2*Pi*y))*(.1552097244200338908e-3*exp(1.4870714009003177286*x+6.283185307179586477)+.18995939286205075122e-2*exp(1.4870714009003177286*x-6.283185307179586477)-.87963673345969899831e-3*exp(1.4870714009003177286*x)-exp(-1.4870714009003177286*x)-.25729430690529175941e-3*exp(2*Pi*x)-.92104789385792596055e-3*exp(2*Pi*(x-1))-.83689489742523677205e-2*exp(-2*Pi*x)+.92104789385792596055e-3*exp(-2*Pi*(x-1)));;
                    u(i1,i2,i3,u2c) = -.92104789385792596051e-3*D1*(cos(6.1046733125517700372*t)*sin(2*Pi*y)-1.*sin(6.1046733125517700372*t)*cos(2*Pi*y))*(-.39883077927999820340e-1*exp(1.4870714009003177286*x+6.283185307179586477)-.48812439407274663365*exp(1.4870714009003177286*x-6.283185307179586477)+.22603364911571255080*exp(1.4870714009003177286*x)-256.96249430909924217*exp(-1.4870714009003177286*x)+.27934954156138603277*exp(2*Pi*x)+1.0000000000000000000*exp(2*Pi*(x-1))-9.0863341961490874695*exp(-2*Pi*x)+1.0000000000000000000*exp(-2*Pi*(x-1)));;
                    u(i1,i2,i3,pc ) = D1*(cos(6.1046733125517700372*t)*cos(2*Pi*y)+sin(6.1046733125517700372*t)*sin(2*Pi*y))*(.54629496939894707374e-2*exp(2*Pi*(x-1))+.15260724925888727652e-2*exp(2*Pi*x)-.49638186616338720965e-1*exp(-2*Pi*x)+.54629496939894707369e-2*exp(-2*Pi*(x-1)));;
                } 

        // At the initial time we initialize some MOL variables
                    if( computeVelocity )
                    {
                        printF("\n $$$$$$$$$$ userDefinedKnownSolution: EVALUATE VELOCITY, ACCEL and PRESSURE AT PAST TIMES t=%9.3e dt=%16.6e $$$$$$$\n\n",t,dt);
                        assert( dt>0. );
                        assert( vgf!=NULL );
                        const int currentVelocity = dbase.get<int>("currentVelocity");
                        assert( currentVelocity>=0 && currentVelocity<numberOfVelocityFunctions );
            // fill in t0=0 and some past time levels
                        for( int level=0; level<numberOfVelocityFunctions; level++ )
                        {
                            real tc = 0. - level*dt;
                            const int prevVelocity = (currentVelocity - level + numberOfVelocityFunctions) % numberOfVelocityFunctions;
                            OV_GET_SERIAL_ARRAY(real,vgf[prevVelocity][grid],vLocal);
                            OV_GET_SERIAL_ARRAY(real,vtgf[prevVelocity][grid],vtLocal);
                            FOR_3D(i1,i2,i3,I1,I2,I3)
                            {
                                Real x = vertex(i1,i2,i3,0);
                                Real y = vertex(i1,i2,i3,1);
                                vLocal(i1,i2,i3,u1c)  = -.56226964972164959617e-2*D1*(cos(6.1046733125517700372*tc)*sin(2*Pi*y)-1.*sin(6.1046733125517700372*tc)*cos(2*Pi*y))*(.16851428189029157646*exp(1.4870714009003177286*x+6.283185307179586477)+2.0624268740942637988*exp(1.4870714009003177286*x-6.283185307179586477)-.95503908029714821575*exp(1.4870714009003177286*x)-1085.7198704525267692*exp(-1.4870714009003177286*x)-.27934954156138603278*exp(2*Pi*x)-1.0000000000000000000*exp(2*Pi*(x-1))-9.0863341961490874696*exp(-2*Pi*x)+1.0000000000000000000*exp(-2*Pi*(x-1)));;
                                vLocal(i1,i2,i3,u2c)  = .56226964972164959616e-2*D1*(cos(6.1046733125517700372*tc)*cos(2*Pi*y)+sin(6.1046733125517700372*tc)*sin(2*Pi*y))*(-.39883077927999820340e-1*exp(1.4870714009003177286*x+6.283185307179586477)-.48812439407274663365*exp(1.4870714009003177286*x-6.283185307179586477)+.22603364911571255080*exp(1.4870714009003177286*x)-256.96249430909924217*exp(-1.4870714009003177286*x)+.27934954156138603277*exp(2*Pi*x)+1.0000000000000000000*exp(2*Pi*(x-1))-9.0863341961490874695*exp(-2*Pi*x)+1.0000000000000000000*exp(-2*Pi*(x-1)));;
                                vtLocal(i1,i2,i3,u1c) = .34324725251135860639e-1*D1*(cos(6.1046733125517700372*tc)*cos(2*Pi*y)+sin(6.1046733125517700372*tc)*sin(2*Pi*y))*(.16851428189029157646*exp(1.4870714009003177286*x+6.283185307179586477)+2.0624268740942637988*exp(1.4870714009003177286*x-6.283185307179586477)-.95503908029714821575*exp(1.4870714009003177286*x)-1085.7198704525267692*exp(-1.4870714009003177286*x)-.27934954156138603278*exp(2*Pi*x)-1.0000000000000000000*exp(2*Pi*(x-1))-9.0863341961490874696*exp(-2*Pi*x)+1.0000000000000000000*exp(-2*Pi*(x-1)));;
                                vtLocal(i1,i2,i3,u2c) = .34324725251135860639e-1*D1*(cos(6.1046733125517700372*tc)*sin(2*Pi*y)-1.*sin(6.1046733125517700372*tc)*cos(2*Pi*y))*(-.39883077927999820340e-1*exp(1.4870714009003177286*x+6.283185307179586477)-.48812439407274663365*exp(1.4870714009003177286*x-6.283185307179586477)+.22603364911571255080*exp(1.4870714009003177286*x)-256.96249430909924217*exp(-1.4870714009003177286*x)+.27934954156138603277*exp(2*Pi*x)+1.0000000000000000000*exp(2*Pi*(x-1))-9.0863341961490874695*exp(-2*Pi*x)+1.0000000000000000000*exp(-2*Pi*(x-1)));;        
                            }         
                        }
            // --- save pressure for time extrapolation ---
                        const int currentPressure = dbase.get<int>("currentPressure");
                        assert( currentPressure>=0 && currentPressure<max(1,numberOfPressureFunctions) );
                        for( int level=0; level<numberOfPressureFunctions; level++ )
                        {
                            real tc = 0. - (level+3)*dt;  // -- save pressure starting at t-3*dt 
                            const int prevPressure = (currentPressure - level + numberOfPressureFunctions) % numberOfPressureFunctions;
                            assert( pgf !=NULL );
                            OV_GET_SERIAL_ARRAY(real,pgf[prevPressure][grid],pLocal);
                            FOR_3D(i1,i2,i3,I1,I2,I3)
                            {
                                Real x = vertex(i1,i2,i3,0);
                                Real y = vertex(i1,i2,i3,1);
                                pLocal(i1,i2,i3) = D1*(cos(6.1046733125517700372*tc)*cos(2*Pi*y)+sin(6.1046733125517700372*tc)*sin(2*Pi*y))*(.54629496939894707374e-2*exp(2*Pi*(x-1))+.15260724925888727652e-2*exp(2*Pi*x)-.49638186616338720965e-1*exp(-2*Pi*x)+.54629496939894707369e-2*exp(-2*Pi*(x-1)));; 
                            }
                        }
                    }
            }
            else
            {
    
                OV_ABORT("UDKS: finish me - case 5 ");
            }
        }        
        else
        {
            OV_ABORT("finish me");
        }


    }
    else if( userKnownSolution=="incompressibleAnnulusSolution" )
    {
    // --- incompressible annulus solution -----
        const int pc = dbase.get<int >("pc");
        assert( pc>=0 );

        MappedGrid & mg = cg[grid];
        mg.update(MappedGrid::THEvertex | MappedGrid::THEcenter);
        OV_GET_SERIAL_ARRAY_CONST(real,mg.vertex(),vertex);
        RealArray & u = ua;

        const int & iswCase = db.get<int>("iswCase");
        const int n = ipar[0];
        const int m = ipar[1];    

        if( t<3*dt )
            printF("UDKS: evaluate the incompressible annulus solution: iswCase=%i, n=%d, m=%d t=%9.3e ( Jn(lambda_nm r) \n",iswCase,n,m,t);

    // if( computeVelocity )
    // {
    //   OV_ABORT("UDKS: finish me -- computeVelocity");
    // }

    // Real C4 =1.; // defines the scale

    // annulusScaling : scale annulus solution to be size one
        if( !db.has_key("annulusScaling") )
        {
            db.put<Real>("annulusScaling") = -1.;  // this means scaling has not be set 
        }
        Real & annulusScaling = db.get<Real>("annulusScaling");
        Real uNorm; 
        if( iswCase==1 )
        {
      // ---- Annulus : displacement-displacement ----
      // printF("UDKS: EVAL NEW ANNULUS SOLUTION displacement-displacement\n");
      // int n=1, m=2; 
            Real c4=1./4.;  // old scaling 

// --- roots of the frequency equation -----
//   File written by diskAnnulusSolution.mw on 2021-10-18
//   Solution is of the form Jn(omega*r) with roots omega_{n,m} 
Real omegaArray[]={
    // n=1, m=1,2,..,4
        1.2510520647340855e+01, 1.7985120095545963e+01, 2.5104649597130992e+01, 3.0907710595982864e+01,
    // n=2, m=1,2,..,4
        1.2334522283695935e+01, 1.8035467567538374e+01, 2.5019752889578412e+01, 3.0937071838665323e+01,
    // n=3, m=1,2,..,4
        1.2176999775790313e+01, 1.8133089033841518e+01, 2.4949743731740190e+01, 3.0993818552125962e+01
};
#define omegaRoot(n,m) omegaArray[(m-1)+4*(n-1)]
// ------------------------------------------------------------------------------------------------ 
//    caseName: Exact solution for incompressible linear elasticity 
// Annulus of outer radius Rb=2, inner radius Ra = 1; r=1 Displacement BC; r=2 Displacement BC. 
// File Written by diskAnnulusSolution.mw (maple) on 2021-10-18 
// Set c3 to scale disk solution and c4 to scale the annulus solution 
// Set n and m to choose the solution, Jn(omega*r), omega=omega(n,m) 
// The eigenvalues are save in a separate include file. 
// ------------------------------------------------------------------------------------------------ 
Real rho,c,mu,omega; 
rho=1.;
mu=1.;
omega=omegaRoot(n,m);
Real t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49
        ,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,t86,t87,t88,t89,t90,t91,t92,t93,t94,t95,t96,t97,t98,t99
        ,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t130,t131,t132,t133,t134,t135,t136,t137,t138,t139,t140,t141,t142,t143,t144,t145,t146,t147,t148,t149
        ,t150,t151,t152,t153,t154,t155,t156,t157,t158,t159,t160,t161,t162,t163,t164,t165,t166,t167,t168,t169,t170,t171,t172,t173,t174,t175,t176,t177,t178,t179,t180,t181,t182,t183,t184,t185,t186,t187,t188,t189,t190,t191,t192,t193,t194,t195,t196,t197,t198,t199
        ,t200,t201,t202,t203,t204,t205,t206,t207,t208,t209,t210,t211,t212,t213,t214,t215,t216,t217,t218,t219,t220,t221,t222,t223,t224,t225,t226,t227,t228,t229,t230,t231,t232,t233,t234,t235,t236,t237,t238,t239,t240,t241,t242,t243,t244,t245,t246,t247,t248,t249
        ,t250,t251,t252,t253,t254,t255,t256,t257,t258,t259,t260,t261,t262,t263,t264,t265,t266,t267,t268,t269,t270,t271,t272,t273,t274,t275,t276,t277,t278,t279,t280,t281,t282,t283,t284,t285,t286,t287,t288,t289,t290,t291,t292,t293,t294,t295,t296,t297,t298,t299
        ,t300,t301,t302,t303,t304,t305,t306,t307,t308,t309,t310,t311,t312,t313,t314,t315,t316,t317,t318,t319,t320,t321,t322,t323,t324,t325,t326,t327,t328,t329,t330,t331,t332,t333,t334,t335,t336,t337,t338,t339,t340,t341,t342,t343,t344,t345,t346,t347,t348,t349
        ,t350,t351,t352,t353,t354,t355,t356,t357,t358,t359,t360,t361,t362,t363,t364,t365,t366,t367,t368,t369,t370,t371,t372,t373,t374,t375,t376,t377,t378,t379,t380,t381,t382,t383,t384,t385,t386,t387,t388,t389,t390,t391,t392,t393,t394,t395,t396,t397,t398,t399
        ,t400;



            if( annulusScaling <= 0. )
            { // -- first time: compute the norm of the solution so we can scale to be size one
                Real t0=0.; 
                    Real ur,vr,pr; 
                    Real u1Max=0., u2Max=0., u3Max=0., pMax=0.; 
                    uNorm=0; 
                    if( numberOfDimensions==2 )
                    {
                        FOR_3D(i1,i2,i3,I1,I2,I3)
                        {
                            Real x = vertex(i1,i2,i3,0);
                            Real y = vertex(i1,i2,i3,1);
                            Real r = sqrt( x*x + y*y );
                            Real theta = atan2(y,x); 
                            Real cosTheta = x/r;
                            Real sinTheta = y/r;
                            t1 = omega * c4;
                            t2 = jn(n, omega);
                            t5 = n + 0.1e1;
                            t6 = omega / 0.2e1;
                            t7 = jn(t5, t6);
                            t8 = omega * t7;
                            t9 = pow(0.1e1 / 0.2e1, n);
                            t12 = 0.1e1 / t9;
                            t15 = jn(n, t6);
                            t20 = 0.1e1 / (-0.8e1 * t12 * t15 * n + 0.8e1 * n * t2 + 0.2e1 * t12 * t8 - 0.2e1 * t9 * t8);
                            t21 = pow(r, n);
                            t22 = 0.1e1 / t21;
                            t23 = t22 * t20;
                            t24 = t23 * t1;
                            t25 = 0.1e1 / r;
                            t26 = yn(n, t6);
                            t27 = t26 * t25;
                            t28 = omega * t0;
                            t29 = cos(t28);
                            t30 = t29 * t7;
                            t31 = n * theta;
                            t32 = cos(t31);
                            t33 = t32 * t30;
                            t34 = t33 * t27;
                            t37 = c4 * t25;
                            t38 = omega * r;
                            t39 = jn(n, t38);
                            t41 = t39 * t20 * t37;
                            t42 = yn(t5, t6);
                            t43 = omega * t42;
                            t44 = sin(t28);
                            t45 = t44 * t9;
                            t46 = sin(t31);
                            t47 = t46 * t45;
                            t51 = t44 * t12;
                            t52 = t46 * t51;
                            t56 = t26 * n;
                            t60 = t21 * t20;
                            t61 = t60 * t1;
                            t62 = t42 * t25;
                            t63 = t44 * t15;
                            t64 = t46 * t63;
                            t65 = t64 * t62;
                            t68 = t44 * t7;
                            t69 = t46 * t68;
                            t70 = t69 * t27;
                            t77 = t29 * t9;
                            t78 = t32 * t77;
                            t82 = t29 * t12;
                            t83 = t32 * t82;
                            t90 = t29 * t15;
                            t91 = t32 * t90;
                            t92 = t91 * t62;
                            t100 = t25 * t60;
                            t101 = t100 * n * c4;
                            t102 = yn(n, omega);
                            t103 = t102 * t12;
                            t104 = t91 * t103;
                            t107 = 0.2e1 * t47 * t43 * t41 - 0.2e1 * t52 * t43 * t41 + 0.2e1 * t78 * t43 * t41 - 0.2e1 * t83 * t43 * t41 + 0.8e1 * t52 * t56 * t41 + 0.8e1 * t83 * t56 * t41 + 0.8e1 * t104 * t101 - 0.2e1 * t34 * t24 + 0.2e1 * t65 * t24 - 0.2e1 * t70 * t24 + 0.2e1 * t92 * t24 + 0.2e1 * t34 * t61 - 0.2e1 * t65 * t61 + 0.2e1 * t70 * t61 - 0.2e1 * t92 * t61;
                            t108 = t64 * t103;
                            t111 = t26 * t12;
                            t112 = t44 * t2;
                            t114 = t46 * t112 * t111;
                            t117 = t100 * t1;
                            t118 = t7 * t102;
                            t119 = t52 * t118;
                            t122 = t2 * t42;
                            t123 = t52 * t122;
                            t127 = t25 * t23 * t1;
                            t128 = t47 * t118;
                            t131 = t47 * t122;
                            t134 = t29 * t2;
                            t136 = t32 * t134 * t111;
                            t139 = t83 * t118;
                            t142 = t83 * t122;
                            t145 = t78 * t118;
                            t148 = t78 * t122;
                            t151 = yn(n, t38);
                            t152 = t29 * t151;
                            t155 = t44 * t151;
                            t158 = n * t102;
                            t167 = -0.8e1 * t32 * t29 * t158 * t41 - 0.8e1 * t46 * t44 * t158 * t41 + t32 * t152 * t37 + t46 * t155 * t37 + 0.8e1 * t108 * t101 - 0.8e1 * t114 * t101 - 0.8e1 * t136 * t101 - 0.2e1 * t119 * t117 + 0.2e1 * t123 * t117 - 0.2e1 * t139 * t117 + 0.2e1 * t142 * t117 + 0.2e1 * t128 * t127 - 0.2e1 * t131 * t127 + 0.2e1 * t145 * t127 - 0.2e1 * t148 * t127;
                            ur = t107 + t167;
                            t168 = jn(t5, t38);
                            t169 = t168 * omega;
                            t170 = c4 * t169;
                            t171 = t102 * t20;
                            t172 = t32 * t44;
                            t176 = t46 * t29;
                            t180 = t46 * t77;
                            t187 = t46 * t90;
                            t191 = t32 * t63;
                            t199 = t32 * t51;
                            t206 = t32 * t45;
                            t217 = t46 * t82;
                            t225 = t46 * t30 * t27;
                            t228 = t187 * t62;
                            t233 = omega * omega;
                            t234 = 0.1e1 / n;
                            t237 = c4 * t168 * t234 * t233;
                            t238 = t42 * t20;
                            t248 = -0.8e1 * t32 * t112 * t111 * t101 + 0.8e1 * t46 * t134 * t111 * t101 - 0.8e1 * t187 * t103 * t101 + 0.8e1 * t191 * t103 * t101 - 0.2e1 * t199 * t118 * t117 + 0.2e1 * t217 * t118 * t117 + 0.2e1 * t199 * t122 * t117 - 0.2e1 * t217 * t122 * t117 + 0.2e1 * t180 * t118 * t127 - 0.2e1 * t206 * t118 * t127 - 0.2e1 * t180 * t122 * t127 + 0.2e1 * t206 * t122 * t127 + 0.8e1 * t172 * t171 * t170 - 0.8e1 * t176 * t171 * t170 + 0.2e1 * t180 * t238 * t237 + 0.2e1 * t199 * t238 * t237 - 0.2e1 * t206 * t238 * t237 - 0.2e1 * t225 * t24 - 0.2e1 * t225 * t61 + 0.2e1 * t228 * t24;
                            t261 = t191 * t62;
                            t265 = t32 * t68 * t27;
                            t284 = yn(t5, t38);
                            t285 = t284 * t234 * omega;
                            t297 = t20 * c4 * t169;
                            t310 = t46 * t29 * c4 * t285 - t32 * t44 * c4 * t285 - 0.8e1 * t172 * t111 * t297 + 0.8e1 * t176 * t111 * t297 - t46 * t152 * t37 + t32 * t155 * t37 - 0.8e1 * t172 * t158 * t41 + 0.8e1 * t176 * t158 * t41 - 0.2e1 * t180 * t43 * t41 - 0.2e1 * t199 * t43 * t41 + 0.8e1 * t199 * t56 * t41 + 0.2e1 * t206 * t43 * t41 - 0.2e1 * t217 * t238 * t237 + 0.2e1 * t217 * t43 * t41 - 0.8e1 * t217 * t56 * t41 + 0.2e1 * t228 * t61 - 0.2e1 * t261 * t24 + 0.2e1 * t265 * t24 - 0.2e1 * t261 * t61 + 0.2e1 * t265 * t61;
                            vr = t248 + t310;
                            t312 = c4 * t233 * omega;
                            t313 = t234 * t20;
                            t315 = t22 * t313 * t312;
                            t319 = t21 * t313 * t312;
                            t335 = t60 * c4 * t233;
                            t342 = t313 * t312;
                            t343 = t42 * t22;
                            t347 = t42 * t21;
                            t351 = t26 * t21;
                            t361 = t26 * t22;
                            pr = 0.2e1 * t33 * t351 * t342 + 0.2e1 * t33 * t361 * t342 - 0.2e1 * t64 * t343 * t342 - 0.2e1 * t91 * t343 * t342 - 0.2e1 * t64 * t347 * t342 - 0.2e1 * t91 * t347 * t342 + 0.2e1 * t69 * t351 * t342 + 0.2e1 * t69 * t361 * t342 + 0.8e1 * t104 * t335 + 0.8e1 * t108 * t335 - 0.8e1 * t114 * t335 - 0.2e1 * t119 * t319 + 0.2e1 * t123 * t319 - 0.2e1 * t128 * t315 + 0.2e1 * t131 * t315 - 0.8e1 * t136 * t335 - 0.2e1 * t139 * t319 + 0.2e1 * t142 * t319 - 0.2e1 * t145 * t315 + 0.2e1 * t148 * t315;
              // (u1,u2) = ur*rHat + vr*thetaHat 
                            u(i1,i2,i3,u1c) = ur*cosTheta - vr*sinTheta; 
                            u(i1,i2,i3,u2c) = ur*sinTheta + vr*cosTheta;
                            u(i1,i2,i3,pc ) = pr;
                            u1Max=max(u1Max,u(i1,i2,i3,u1c));
                            u2Max=max(u2Max,u(i1,i2,i3,u2c));
                            pMax =max(pMax, u(i1,i2,i3,pc));
                        }
                        uNorm = max(u1Max,u2Max,u3Max,pMax);
            // if( scale==0. )
            // {
            //   printF("UDKS: ERROR: annulus solution is zero! scale=%e\n",scale);
            //   OV_ABORT("ERROR");
            // }
            // scale = 1./scale;
            // c4 *= scale;  // set this so the velocity and acceleration are scaled too
            // FOR_3(i1,i2,i3,I1,I2,I3)
            // {
            //   u(i1,i2,i3,u1c) *= scale;
            //   u(i1,i2,i3,u2c) *= scale;
            //   u(i1,i2,i3,pc ) *= scale;
            // }    
                    }
                    else
                    {
                        OV_ABORT("UDKS:ERROR: annulus solution : 3D called");
                    }
                if( uNorm==0. )
                {
                    printF("UDKS: ERROR: annulus solution is zero! uNorm=%e\n",uNorm);
                    OV_ABORT("ERROR");
                }        
                annulusScaling = c4/uNorm;
            }
            c4 = annulusScaling; 

                Real ur,vr,pr; 
                Real u1Max=0., u2Max=0., u3Max=0., pMax=0.; 
                uNorm=0; 
                if( numberOfDimensions==2 )
                {
                    FOR_3D(i1,i2,i3,I1,I2,I3)
                    {
                        Real x = vertex(i1,i2,i3,0);
                        Real y = vertex(i1,i2,i3,1);
                        Real r = sqrt( x*x + y*y );
                        Real theta = atan2(y,x); 
                        Real cosTheta = x/r;
                        Real sinTheta = y/r;
                        t1 = omega * c4;
                        t2 = jn(n, omega);
                        t5 = n + 0.1e1;
                        t6 = omega / 0.2e1;
                        t7 = jn(t5, t6);
                        t8 = omega * t7;
                        t9 = pow(0.1e1 / 0.2e1, n);
                        t12 = 0.1e1 / t9;
                        t15 = jn(n, t6);
                        t20 = 0.1e1 / (-0.8e1 * t12 * t15 * n + 0.8e1 * n * t2 + 0.2e1 * t12 * t8 - 0.2e1 * t9 * t8);
                        t21 = pow(r, n);
                        t22 = 0.1e1 / t21;
                        t23 = t22 * t20;
                        t24 = t23 * t1;
                        t25 = 0.1e1 / r;
                        t26 = yn(n, t6);
                        t27 = t26 * t25;
                        t28 = omega * t;
                        t29 = cos(t28);
                        t30 = t29 * t7;
                        t31 = n * theta;
                        t32 = cos(t31);
                        t33 = t32 * t30;
                        t34 = t33 * t27;
                        t37 = c4 * t25;
                        t38 = omega * r;
                        t39 = jn(n, t38);
                        t41 = t39 * t20 * t37;
                        t42 = yn(t5, t6);
                        t43 = omega * t42;
                        t44 = sin(t28);
                        t45 = t44 * t9;
                        t46 = sin(t31);
                        t47 = t46 * t45;
                        t51 = t44 * t12;
                        t52 = t46 * t51;
                        t56 = t26 * n;
                        t60 = t21 * t20;
                        t61 = t60 * t1;
                        t62 = t42 * t25;
                        t63 = t44 * t15;
                        t64 = t46 * t63;
                        t65 = t64 * t62;
                        t68 = t44 * t7;
                        t69 = t46 * t68;
                        t70 = t69 * t27;
                        t77 = t29 * t9;
                        t78 = t32 * t77;
                        t82 = t29 * t12;
                        t83 = t32 * t82;
                        t90 = t29 * t15;
                        t91 = t32 * t90;
                        t92 = t91 * t62;
                        t100 = t25 * t60;
                        t101 = t100 * n * c4;
                        t102 = yn(n, omega);
                        t103 = t102 * t12;
                        t104 = t91 * t103;
                        t107 = 0.2e1 * t47 * t43 * t41 - 0.2e1 * t52 * t43 * t41 + 0.2e1 * t78 * t43 * t41 - 0.2e1 * t83 * t43 * t41 + 0.8e1 * t52 * t56 * t41 + 0.8e1 * t83 * t56 * t41 + 0.8e1 * t104 * t101 - 0.2e1 * t34 * t24 + 0.2e1 * t65 * t24 - 0.2e1 * t70 * t24 + 0.2e1 * t92 * t24 + 0.2e1 * t34 * t61 - 0.2e1 * t65 * t61 + 0.2e1 * t70 * t61 - 0.2e1 * t92 * t61;
                        t108 = t64 * t103;
                        t111 = t26 * t12;
                        t112 = t44 * t2;
                        t114 = t46 * t112 * t111;
                        t117 = t100 * t1;
                        t118 = t7 * t102;
                        t119 = t52 * t118;
                        t122 = t2 * t42;
                        t123 = t52 * t122;
                        t127 = t25 * t23 * t1;
                        t128 = t47 * t118;
                        t131 = t47 * t122;
                        t134 = t29 * t2;
                        t136 = t32 * t134 * t111;
                        t139 = t83 * t118;
                        t142 = t83 * t122;
                        t145 = t78 * t118;
                        t148 = t78 * t122;
                        t151 = yn(n, t38);
                        t152 = t29 * t151;
                        t155 = t44 * t151;
                        t158 = n * t102;
                        t167 = -0.8e1 * t32 * t29 * t158 * t41 - 0.8e1 * t46 * t44 * t158 * t41 + t32 * t152 * t37 + t46 * t155 * t37 + 0.8e1 * t108 * t101 - 0.8e1 * t114 * t101 - 0.8e1 * t136 * t101 - 0.2e1 * t119 * t117 + 0.2e1 * t123 * t117 - 0.2e1 * t139 * t117 + 0.2e1 * t142 * t117 + 0.2e1 * t128 * t127 - 0.2e1 * t131 * t127 + 0.2e1 * t145 * t127 - 0.2e1 * t148 * t127;
                        ur = t107 + t167;
                        t168 = jn(t5, t38);
                        t169 = t168 * omega;
                        t170 = c4 * t169;
                        t171 = t102 * t20;
                        t172 = t32 * t44;
                        t176 = t46 * t29;
                        t180 = t46 * t77;
                        t187 = t46 * t90;
                        t191 = t32 * t63;
                        t199 = t32 * t51;
                        t206 = t32 * t45;
                        t217 = t46 * t82;
                        t225 = t46 * t30 * t27;
                        t228 = t187 * t62;
                        t233 = omega * omega;
                        t234 = 0.1e1 / n;
                        t237 = c4 * t168 * t234 * t233;
                        t238 = t42 * t20;
                        t248 = -0.8e1 * t32 * t112 * t111 * t101 + 0.8e1 * t46 * t134 * t111 * t101 - 0.8e1 * t187 * t103 * t101 + 0.8e1 * t191 * t103 * t101 - 0.2e1 * t199 * t118 * t117 + 0.2e1 * t217 * t118 * t117 + 0.2e1 * t199 * t122 * t117 - 0.2e1 * t217 * t122 * t117 + 0.2e1 * t180 * t118 * t127 - 0.2e1 * t206 * t118 * t127 - 0.2e1 * t180 * t122 * t127 + 0.2e1 * t206 * t122 * t127 + 0.8e1 * t172 * t171 * t170 - 0.8e1 * t176 * t171 * t170 + 0.2e1 * t180 * t238 * t237 + 0.2e1 * t199 * t238 * t237 - 0.2e1 * t206 * t238 * t237 - 0.2e1 * t225 * t24 - 0.2e1 * t225 * t61 + 0.2e1 * t228 * t24;
                        t261 = t191 * t62;
                        t265 = t32 * t68 * t27;
                        t284 = yn(t5, t38);
                        t285 = t284 * t234 * omega;
                        t297 = t20 * c4 * t169;
                        t310 = t46 * t29 * c4 * t285 - t32 * t44 * c4 * t285 - 0.8e1 * t172 * t111 * t297 + 0.8e1 * t176 * t111 * t297 - t46 * t152 * t37 + t32 * t155 * t37 - 0.8e1 * t172 * t158 * t41 + 0.8e1 * t176 * t158 * t41 - 0.2e1 * t180 * t43 * t41 - 0.2e1 * t199 * t43 * t41 + 0.8e1 * t199 * t56 * t41 + 0.2e1 * t206 * t43 * t41 - 0.2e1 * t217 * t238 * t237 + 0.2e1 * t217 * t43 * t41 - 0.8e1 * t217 * t56 * t41 + 0.2e1 * t228 * t61 - 0.2e1 * t261 * t24 + 0.2e1 * t265 * t24 - 0.2e1 * t261 * t61 + 0.2e1 * t265 * t61;
                        vr = t248 + t310;
                        t312 = c4 * t233 * omega;
                        t313 = t234 * t20;
                        t315 = t22 * t313 * t312;
                        t319 = t21 * t313 * t312;
                        t335 = t60 * c4 * t233;
                        t342 = t313 * t312;
                        t343 = t42 * t22;
                        t347 = t42 * t21;
                        t351 = t26 * t21;
                        t361 = t26 * t22;
                        pr = 0.2e1 * t33 * t351 * t342 + 0.2e1 * t33 * t361 * t342 - 0.2e1 * t64 * t343 * t342 - 0.2e1 * t91 * t343 * t342 - 0.2e1 * t64 * t347 * t342 - 0.2e1 * t91 * t347 * t342 + 0.2e1 * t69 * t351 * t342 + 0.2e1 * t69 * t361 * t342 + 0.8e1 * t104 * t335 + 0.8e1 * t108 * t335 - 0.8e1 * t114 * t335 - 0.2e1 * t119 * t319 + 0.2e1 * t123 * t319 - 0.2e1 * t128 * t315 + 0.2e1 * t131 * t315 - 0.8e1 * t136 * t335 - 0.2e1 * t139 * t319 + 0.2e1 * t142 * t319 - 0.2e1 * t145 * t315 + 0.2e1 * t148 * t315;
            // (u1,u2) = ur*rHat + vr*thetaHat 
                        u(i1,i2,i3,u1c) = ur*cosTheta - vr*sinTheta; 
                        u(i1,i2,i3,u2c) = ur*sinTheta + vr*cosTheta;
                        u(i1,i2,i3,pc ) = pr;
                        u1Max=max(u1Max,u(i1,i2,i3,u1c));
                        u2Max=max(u2Max,u(i1,i2,i3,u2c));
                        pMax =max(pMax, u(i1,i2,i3,pc));
                    }
                    uNorm = max(u1Max,u2Max,u3Max,pMax);
          // if( scale==0. )
          // {
          //   printF("UDKS: ERROR: annulus solution is zero! scale=%e\n",scale);
          //   OV_ABORT("ERROR");
          // }
          // scale = 1./scale;
          // c4 *= scale;  // set this so the velocity and acceleration are scaled too
          // FOR_3(i1,i2,i3,I1,I2,I3)
          // {
          //   u(i1,i2,i3,u1c) *= scale;
          //   u(i1,i2,i3,u2c) *= scale;
          //   u(i1,i2,i3,pc ) *= scale;
          // }    
                }
                else
                {
                    OV_ABORT("UDKS:ERROR: annulus solution : 3D called");
                }

                if( computeVelocity )
                {
                    printF("\n $$$$$$$$$$ userDefinedKnownSolution: EVALUATE VELOCITY, ACCEL and PRESSURE AT PAST TIMES t=%9.3e dt=%16.6e $$$$$$$\n\n",t,dt);
                    assert( dt>0. );
                    assert( vgf!=NULL );
                    const int currentVelocity = dbase.get<int>("currentVelocity");
                    assert( currentVelocity>=0 && currentVelocity<numberOfVelocityFunctions );
          // fill in t0=0 and some past time levels
                    for( int level=0; level<numberOfVelocityFunctions; level++ )
                    {
                        real tc = 0. - level*dt;
                        const int prevVelocity = (currentVelocity - level + numberOfVelocityFunctions) % numberOfVelocityFunctions;
                        OV_GET_SERIAL_ARRAY(real,vgf[prevVelocity][grid],vLocal);
                        OV_GET_SERIAL_ARRAY(real,vtgf[prevVelocity][grid],vtLocal);
                        FOR_3D(i1,i2,i3,I1,I2,I3)
                        {
                            Real x = vertex(i1,i2,i3,0);
                            Real y = vertex(i1,i2,i3,1);
                            Real r = sqrt( x*x + y*y );
                            Real theta = atan2(y,x); 
                            Real cosTheta = x/r;
                            Real sinTheta = y/r;
                            t1 = omega * tc;
                            t2 = sin(t1);
                            t3 = n * theta;
                            t4 = cos(t3);
                            t6 = cos(t1);
                            t7 = sin(t3);
                            t11 = pow(0.4e1, n);
                            t12 = pow(r, n);
                            t13 = pow(r, -n);
                            t15 = (t12 - t13) * t11;
                            t16 = omega / 0.2e1;
                            t17 = yn(n, t16);
                            t19 = pow(0.2e1, n);
                            t20 = pow(0.8e1, n);
                            t21 = t19 - t20;
                            t22 = omega * r;
                            t23 = yn(n, t22);
                            t24 = t23 * t21;
                            t25 = yn(n, omega);
                            t26 = t12 * t20;
                            t27 = t13 * t19;
                            t28 = t26 - t27;
                            t32 = n + 0.1e1;
                            t33 = jn(t32, t16);
                            t36 = jn(n, t16);
                            t38 = jn(n, t22);
                            t39 = t38 * t21;
                            t40 = jn(n, omega);
                            t44 = yn(t32, t16);
                            t49 = -t40 * t23 + t25 * t38;
                            t52 = pow(0.2e1, 0.2e1 * t32);
                            t60 = t36 * (-t25 * t12 + t23) + (t40 * t12 - t38) * t17;
                            t67 = omega / r;
                            t76 = 0.1e1 / (-t33 * t21 * omega + (-0.4e1 * t20 * t36 + t40 * t52) * n);
                            ur = 0.4e1 * t76 * t67 * (t33 * (-t17 * t15 + t28 * t25 + t24) * omega / 0.4e1 - t44 * omega * (-t36 * t15 + t28 * t40 + t39) / 0.4e1 + n * (t52 * t49 / 0.4e1 + t60 * t20)) * (t4 * t2 - t7 * t6) * c4;
                            t79 = yn(t32, t22);
                            t83 = (t12 + t13) * t11;
                            t85 = t26 + t27;
                            t93 = jn(t32, t22);
                            t104 = t93 * omega;
                            vr = -0.4e1 * t76 / n * t67 * (t7 * t2 + t6 * t4) * (t33 * omega * (-t79 * t21 * t22 + n * (-t17 * t83 + t85 * t25 + t24)) / 0.4e1 - t44 * (-t93 * t21 * t22 + n * (-t36 * t83 + t85 * t40 + t39)) * omega / 0.4e1 + n * (t52 * (-r * t25 * t104 / 0.4e1 + t22 * t40 * t79 / 0.4e1 + t49 * n / 0.4e1) + (-t79 * omega * t36 * r + t104 * t17 * r + n * t60) * t20)) * c4;
              // (u1,u2) = ur*rHat + vr*thetaHat 
                            vLocal(i1,i2,i3,u1c) = ur*cosTheta - vr*sinTheta;
                            vLocal(i1,i2,i3,u2c) = ur*sinTheta + vr*cosTheta;
                            t1 = pow(0.4e1, n);
                            t2 = pow(r, n);
                            t3 = pow(r, -n);
                            t5 = (t2 - t3) * t1;
                            t6 = omega / 0.2e1;
                            t7 = yn(n, t6);
                            t9 = pow(0.2e1, n);
                            t10 = pow(0.8e1, n);
                            t11 = t9 - t10;
                            t12 = omega * r;
                            t13 = yn(n, t12);
                            t14 = t13 * t11;
                            t15 = yn(n, omega);
                            t16 = t2 * t10;
                            t17 = t3 * t9;
                            t18 = t16 - t17;
                            t22 = n + 0.1e1;
                            t23 = jn(t22, t6);
                            t26 = jn(n, t6);
                            t28 = jn(n, t12);
                            t29 = t28 * t11;
                            t30 = jn(n, omega);
                            t34 = yn(t22, t6);
                            t39 = -t30 * t13 + t15 * t28;
                            t42 = pow(0.2e1, 0.2e1 * t22);
                            t50 = t26 * (-t15 * t2 + t13) + (t30 * t2 - t28) * t7;
                            t56 = omega * tc;
                            t57 = sin(t56);
                            t58 = n * theta;
                            t59 = sin(t58);
                            t61 = cos(t58);
                            t62 = cos(t56);
                            t66 = omega * omega;
                            t68 = 0.1e1 / r * t66;
                            t77 = 0.1e1 / (-t23 * t11 * omega + (-0.4e1 * t10 * t26 + t30 * t42) * n);
                            ur = 0.4e1 * t77 * t68 * (t59 * t57 + t62 * t61) * (t23 * (t18 * t15 - t7 * t5 + t14) * omega / 0.4e1 - t34 * omega * (t18 * t30 - t26 * t5 + t29) / 0.4e1 + n * (t42 * t39 / 0.4e1 + t50 * t10)) * c4;
                            t84 = yn(t22, t12);
                            t88 = (t2 + t3) * t1;
                            t90 = t16 + t17;
                            t98 = jn(t22, t12);
                            t109 = t98 * omega;
                            vr = 0.4e1 * t77 / n * t68 * (t23 * omega * (-t84 * t11 * t12 + n * (t90 * t15 - t7 * t88 + t14)) / 0.4e1 - t34 * (-t98 * t11 * t12 + n * (-t26 * t88 + t90 * t30 + t29)) * omega / 0.4e1 + n * (t42 * (-r * t15 * t109 / 0.4e1 + t12 * t30 * t84 / 0.4e1 + t39 * n / 0.4e1) + (-t84 * omega * t26 * r + t109 * t7 * r + n * t50) * t10)) * (t61 * t57 - t59 * t62) * c4;
              // (u1,u2) = ur*rHat + vr*thetaHat 
                            vtLocal(i1,i2,i3,u1c) = ur*cosTheta - vr*sinTheta;
                            vtLocal(i1,i2,i3,u2c) = ur*sinTheta + vr*cosTheta;        
                        }      
            // FOR_3D(i1,i2,i3,I1,I2,I3)
            // {
            //   Real x = vertex(i1,i2,i3,0);
            //   Real y = vertex(i1,i2,i3,1);
            //   vLocal(i1,i2,i3,u1c)  = v1e(x,y,tc);
            //   vLocal(i1,i2,i3,u2c)  = v2e(x,y,tc);
            //   vtLocal(i1,i2,i3,u1c) = v1et(x,y,tc);
            //   vtLocal(i1,i2,i3,u2c) = v2et(x,y,tc);        
            // }         
                    }
          // --- save pressure for time extrapolation ---
                    const int currentPressure = dbase.get<int>("currentPressure");
                    assert( currentPressure>=0 && currentPressure<max(1,numberOfPressureFunctions) );
                    for( int level=0; level<numberOfPressureFunctions; level++ )
                    {
                        real tc = 0. - (level+3)*dt;  // -- save pressure starting at t-3*dt 
                        const int prevPressure = (currentPressure - level + numberOfPressureFunctions) % numberOfPressureFunctions;
                        assert( pgf !=NULL );
                        OV_GET_SERIAL_ARRAY(real,pgf[prevPressure][grid],pLocal);
                        FOR_3D(i1,i2,i3,I1,I2,I3)
                        {
                            Real x = vertex(i1,i2,i3,0);
                            Real y = vertex(i1,i2,i3,1);
                            Real r = sqrt( x*x + y*y );
                            Real theta = atan2(y,x); 
                            Real cosTheta = x/r;
                            Real sinTheta = y/r;
                            t1 = omega * c4;
                            t2 = jn(n, omega);
                            t5 = n + 0.1e1;
                            t6 = omega / 0.2e1;
                            t7 = jn(t5, t6);
                            t8 = omega * t7;
                            t9 = pow(0.1e1 / 0.2e1, n);
                            t12 = 0.1e1 / t9;
                            t15 = jn(n, t6);
                            t20 = 0.1e1 / (-0.8e1 * t12 * t15 * n + 0.8e1 * n * t2 + 0.2e1 * t12 * t8 - 0.2e1 * t9 * t8);
                            t21 = pow(r, n);
                            t22 = 0.1e1 / t21;
                            t23 = t22 * t20;
                            t24 = t23 * t1;
                            t25 = 0.1e1 / r;
                            t26 = yn(n, t6);
                            t27 = t26 * t25;
                            t28 = omega * tc;
                            t29 = cos(t28);
                            t30 = t29 * t7;
                            t31 = n * theta;
                            t32 = cos(t31);
                            t33 = t32 * t30;
                            t34 = t33 * t27;
                            t37 = c4 * t25;
                            t38 = omega * r;
                            t39 = jn(n, t38);
                            t41 = t39 * t20 * t37;
                            t42 = yn(t5, t6);
                            t43 = omega * t42;
                            t44 = sin(t28);
                            t45 = t44 * t9;
                            t46 = sin(t31);
                            t47 = t46 * t45;
                            t51 = t44 * t12;
                            t52 = t46 * t51;
                            t56 = t26 * n;
                            t60 = t21 * t20;
                            t61 = t60 * t1;
                            t62 = t42 * t25;
                            t63 = t44 * t15;
                            t64 = t46 * t63;
                            t65 = t64 * t62;
                            t68 = t44 * t7;
                            t69 = t46 * t68;
                            t70 = t69 * t27;
                            t77 = t29 * t9;
                            t78 = t32 * t77;
                            t82 = t29 * t12;
                            t83 = t32 * t82;
                            t90 = t29 * t15;
                            t91 = t32 * t90;
                            t92 = t91 * t62;
                            t100 = t25 * t60;
                            t101 = t100 * n * c4;
                            t102 = yn(n, omega);
                            t103 = t102 * t12;
                            t104 = t91 * t103;
                            t107 = 0.2e1 * t47 * t43 * t41 - 0.2e1 * t52 * t43 * t41 + 0.2e1 * t78 * t43 * t41 - 0.2e1 * t83 * t43 * t41 + 0.8e1 * t52 * t56 * t41 + 0.8e1 * t83 * t56 * t41 + 0.8e1 * t104 * t101 - 0.2e1 * t34 * t24 + 0.2e1 * t65 * t24 - 0.2e1 * t70 * t24 + 0.2e1 * t92 * t24 + 0.2e1 * t34 * t61 - 0.2e1 * t65 * t61 + 0.2e1 * t70 * t61 - 0.2e1 * t92 * t61;
                            t108 = t64 * t103;
                            t111 = t26 * t12;
                            t112 = t44 * t2;
                            t114 = t46 * t112 * t111;
                            t117 = t100 * t1;
                            t118 = t7 * t102;
                            t119 = t52 * t118;
                            t122 = t2 * t42;
                            t123 = t52 * t122;
                            t127 = t25 * t23 * t1;
                            t128 = t47 * t118;
                            t131 = t47 * t122;
                            t134 = t29 * t2;
                            t136 = t32 * t134 * t111;
                            t139 = t83 * t118;
                            t142 = t83 * t122;
                            t145 = t78 * t118;
                            t148 = t78 * t122;
                            t151 = yn(n, t38);
                            t152 = t29 * t151;
                            t155 = t44 * t151;
                            t158 = n * t102;
                            t167 = -0.8e1 * t32 * t29 * t158 * t41 - 0.8e1 * t46 * t44 * t158 * t41 + t32 * t152 * t37 + t46 * t155 * t37 + 0.8e1 * t108 * t101 - 0.8e1 * t114 * t101 - 0.8e1 * t136 * t101 - 0.2e1 * t119 * t117 + 0.2e1 * t123 * t117 - 0.2e1 * t139 * t117 + 0.2e1 * t142 * t117 + 0.2e1 * t128 * t127 - 0.2e1 * t131 * t127 + 0.2e1 * t145 * t127 - 0.2e1 * t148 * t127;
                            ur = t107 + t167;
                            t168 = jn(t5, t38);
                            t169 = t168 * omega;
                            t170 = c4 * t169;
                            t171 = t102 * t20;
                            t172 = t32 * t44;
                            t176 = t46 * t29;
                            t180 = t46 * t77;
                            t187 = t46 * t90;
                            t191 = t32 * t63;
                            t199 = t32 * t51;
                            t206 = t32 * t45;
                            t217 = t46 * t82;
                            t225 = t46 * t30 * t27;
                            t228 = t187 * t62;
                            t233 = omega * omega;
                            t234 = 0.1e1 / n;
                            t237 = c4 * t168 * t234 * t233;
                            t238 = t42 * t20;
                            t248 = -0.8e1 * t32 * t112 * t111 * t101 + 0.8e1 * t46 * t134 * t111 * t101 - 0.8e1 * t187 * t103 * t101 + 0.8e1 * t191 * t103 * t101 - 0.2e1 * t199 * t118 * t117 + 0.2e1 * t217 * t118 * t117 + 0.2e1 * t199 * t122 * t117 - 0.2e1 * t217 * t122 * t117 + 0.2e1 * t180 * t118 * t127 - 0.2e1 * t206 * t118 * t127 - 0.2e1 * t180 * t122 * t127 + 0.2e1 * t206 * t122 * t127 + 0.8e1 * t172 * t171 * t170 - 0.8e1 * t176 * t171 * t170 + 0.2e1 * t180 * t238 * t237 + 0.2e1 * t199 * t238 * t237 - 0.2e1 * t206 * t238 * t237 - 0.2e1 * t225 * t24 - 0.2e1 * t225 * t61 + 0.2e1 * t228 * t24;
                            t261 = t191 * t62;
                            t265 = t32 * t68 * t27;
                            t284 = yn(t5, t38);
                            t285 = t284 * t234 * omega;
                            t297 = t20 * c4 * t169;
                            t310 = t46 * t29 * c4 * t285 - t32 * t44 * c4 * t285 - 0.8e1 * t172 * t111 * t297 + 0.8e1 * t176 * t111 * t297 - t46 * t152 * t37 + t32 * t155 * t37 - 0.8e1 * t172 * t158 * t41 + 0.8e1 * t176 * t158 * t41 - 0.2e1 * t180 * t43 * t41 - 0.2e1 * t199 * t43 * t41 + 0.8e1 * t199 * t56 * t41 + 0.2e1 * t206 * t43 * t41 - 0.2e1 * t217 * t238 * t237 + 0.2e1 * t217 * t43 * t41 - 0.8e1 * t217 * t56 * t41 + 0.2e1 * t228 * t61 - 0.2e1 * t261 * t24 + 0.2e1 * t265 * t24 - 0.2e1 * t261 * t61 + 0.2e1 * t265 * t61;
                            vr = t248 + t310;
                            t312 = c4 * t233 * omega;
                            t313 = t234 * t20;
                            t315 = t22 * t313 * t312;
                            t319 = t21 * t313 * t312;
                            t335 = t60 * c4 * t233;
                            t342 = t313 * t312;
                            t343 = t42 * t22;
                            t347 = t42 * t21;
                            t351 = t26 * t21;
                            t361 = t26 * t22;
                            pr = 0.2e1 * t33 * t351 * t342 + 0.2e1 * t33 * t361 * t342 - 0.2e1 * t64 * t343 * t342 - 0.2e1 * t91 * t343 * t342 - 0.2e1 * t64 * t347 * t342 - 0.2e1 * t91 * t347 * t342 + 0.2e1 * t69 * t351 * t342 + 0.2e1 * t69 * t361 * t342 + 0.8e1 * t104 * t335 + 0.8e1 * t108 * t335 - 0.8e1 * t114 * t335 - 0.2e1 * t119 * t319 + 0.2e1 * t123 * t319 - 0.2e1 * t128 * t315 + 0.2e1 * t131 * t315 - 0.8e1 * t136 * t335 - 0.2e1 * t139 * t319 + 0.2e1 * t142 * t319 - 0.2e1 * t145 * t315 + 0.2e1 * t148 * t315;
                            pLocal(i1,i2,i3) = pr;
                        }
                    }
                }

        }
        else if( iswCase==2 )
        {
      // ---- Annulus : traction-displacement ----
      // printF("UDKS: EVAL NEW ANNULUS SOLUTION\n");
      // int n=1, m=2; 
            Real c4=1./4.; 

// --- roots of the frequency equation -----
//   File written by diskAnnulusSolution.mw on 2021-10-18
//   Solution is of the form Jn(omega*r) with roots omega_{n,m} 
Real omegaArray[]={
    // n=1, m=1,2,..,4
        5.9301137018497625e+00, 1.0041228867524906e+01, 1.6222236245157285e+01, 2.2261799618829349e+01,
    // n=2, m=1,2,..,4
        7.3209270736757316e+00, 1.0788308757534492e+01, 1.6576115170026591e+01, 2.2312959733123941e+01,
    // n=3, m=1,2,..,4
        7.6148014507882346e+00, 1.2347886231096596e+01, 1.7321954724146530e+01, 2.2621546021972616e+01
};
#define omegaRoot(n,m) omegaArray[(m-1)+4*(n-1)]
// ------------------------------------------------------------------------------------------------ 
//    caseName: Exact solution for incompressible linear elasticity 
// Annulus of outer radius Rb=2, inner radius Ra = 1; r=1 Traction BC; r=2 Displacement BC. 
// File Written by diskAnnulusSolution.mw (maple) on 2021-10-18 
// Set c3 to scale disk solution and c4 to scale the annulus solution 
// Set n and m to choose the solution, Jn(omega*r), omega=omega(n,m) 
// The eigenvalues are save in a separate include file. 
// ------------------------------------------------------------------------------------------------ 
Real rho,c,mu,omega; 
rho=1.;
mu=1.;
omega=omegaRoot(n,m);
Real t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49
        ,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,t86,t87,t88,t89,t90,t91,t92,t93,t94,t95,t96,t97,t98,t99
        ,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t130,t131,t132,t133,t134,t135,t136,t137,t138,t139,t140,t141,t142,t143,t144,t145,t146,t147,t148,t149
        ,t150,t151,t152,t153,t154,t155,t156,t157,t158,t159,t160,t161,t162,t163,t164,t165,t166,t167,t168,t169,t170,t171,t172,t173,t174,t175,t176,t177,t178,t179,t180,t181,t182,t183,t184,t185,t186,t187,t188,t189,t190,t191,t192,t193,t194,t195,t196,t197,t198,t199
        ,t200,t201,t202,t203,t204,t205,t206,t207,t208,t209,t210,t211,t212,t213,t214,t215,t216,t217,t218,t219,t220,t221,t222,t223,t224,t225,t226,t227,t228,t229,t230,t231,t232,t233,t234,t235,t236,t237,t238,t239,t240,t241,t242,t243,t244,t245,t246,t247,t248,t249
        ,t250,t251,t252,t253,t254,t255,t256,t257,t258,t259,t260,t261,t262,t263,t264,t265,t266,t267,t268,t269,t270,t271,t272,t273,t274,t275,t276,t277,t278,t279,t280,t281,t282,t283,t284,t285,t286,t287,t288,t289,t290,t291,t292,t293,t294,t295,t296,t297,t298,t299
        ,t300,t301,t302,t303,t304,t305,t306,t307,t308,t309,t310,t311,t312,t313,t314,t315,t316,t317,t318,t319,t320,t321,t322,t323,t324,t325,t326,t327,t328,t329,t330,t331,t332,t333,t334,t335,t336,t337,t338,t339,t340,t341,t342,t343,t344,t345,t346,t347,t348,t349
        ,t350,t351,t352,t353,t354,t355,t356,t357,t358,t359,t360,t361,t362,t363,t364,t365,t366,t367,t368,t369,t370,t371,t372,t373,t374,t375,t376,t377,t378,t379,t380,t381,t382,t383,t384,t385,t386,t387,t388,t389,t390,t391,t392,t393,t394,t395,t396,t397,t398,t399
        ,t400,t401,t402,t403,t404,t405,t406,t407,t408,t409,t410,t411,t412,t413,t414,t415,t416,t417,t418,t419,t420,t421,t422,t423,t424,t425,t426,t427,t428,t429,t430,t431,t432,t433,t434,t435,t436,t437,t438,t439,t440,t441,t442,t443,t444,t445,t446,t447,t448,t449
        ,t450,t451,t452,t453,t454,t455,t456,t457,t458,t459,t460,t461,t462,t463,t464,t465,t466,t467,t468,t469,t470,t471,t472,t473,t474,t475,t476,t477,t478,t479,t480,t481,t482,t483,t484,t485,t486,t487,t488,t489,t490,t491,t492,t493,t494,t495,t496,t497,t498,t499
        ,t500,t501,t502,t503,t504,t505,t506,t507,t508,t509,t510,t511,t512,t513,t514,t515,t516,t517,t518,t519,t520,t521,t522,t523,t524,t525,t526,t527,t528,t529,t530,t531,t532,t533,t534,t535,t536,t537,t538,t539,t540,t541,t542,t543,t544,t545,t546,t547,t548,t549
        ,t550,t551,t552,t553,t554,t555,t556,t557,t558,t559,t560,t561,t562,t563,t564,t565,t566,t567,t568,t569,t570,t571,t572,t573,t574,t575,t576,t577,t578,t579,t580,t581,t582,t583,t584,t585,t586,t587,t588,t589,t590,t591,t592,t593,t594,t595,t596,t597,t598,t599
        ,t600,t601,t602,t603,t604,t605,t606,t607,t608,t609,t610,t611,t612,t613,t614,t615,t616,t617,t618,t619,t620,t621,t622,t623,t624,t625,t626,t627,t628,t629,t630,t631,t632,t633,t634,t635,t636,t637,t638,t639,t640,t641,t642,t643,t644,t645,t646,t647,t648,t649
        ,t650,t651,t652,t653,t654,t655,t656,t657,t658,t659,t660,t661,t662,t663,t664,t665,t666,t667,t668,t669,t670,t671,t672,t673,t674,t675,t676,t677,t678,t679,t680,t681,t682,t683,t684,t685,t686,t687,t688,t689,t690,t691,t692,t693,t694,t695,t696,t697,t698,t699
        ,t700,t701,t702,t703,t704,t705,t706,t707,t708,t709,t710,t711,t712,t713,t714,t715,t716,t717,t718,t719,t720,t721,t722,t723,t724,t725,t726,t727,t728,t729,t730,t731,t732,t733,t734,t735,t736,t737,t738,t739,t740,t741,t742,t743,t744,t745,t746,t747,t748,t749
        ,t750,t751,t752,t753,t754,t755,t756,t757,t758,t759,t760,t761,t762,t763,t764,t765,t766,t767,t768,t769,t770,t771,t772,t773,t774,t775,t776,t777,t778,t779,t780,t781,t782,t783,t784,t785,t786,t787,t788,t789,t790,t791,t792,t793,t794,t795,t796,t797,t798,t799
        ,t800,t801,t802,t803,t804,t805,t806,t807,t808,t809,t810,t811,t812,t813,t814,t815,t816,t817,t818,t819,t820,t821,t822,t823,t824,t825,t826,t827,t828,t829,t830,t831,t832,t833,t834,t835,t836,t837,t838,t839,t840,t841,t842,t843,t844,t845,t846,t847,t848,t849
        ,t850,t851,t852,t853,t854,t855,t856,t857,t858,t859,t860,t861,t862,t863,t864,t865,t866,t867,t868,t869,t870,t871,t872,t873,t874,t875,t876,t877,t878,t879,t880,t881,t882,t883,t884,t885,t886,t887,t888,t889,t890,t891,t892,t893,t894,t895,t896,t897,t898,t899
        ,t900,t901,t902,t903,t904,t905,t906,t907,t908,t909,t910,t911,t912,t913,t914,t915,t916,t917,t918,t919,t920,t921,t922,t923,t924,t925,t926,t927,t928,t929,t930,t931,t932,t933,t934,t935,t936,t937,t938,t939,t940,t941,t942,t943,t944,t945,t946,t947,t948,t949
        ,t950,t951,t952,t953,t954,t955,t956,t957,t958,t959,t960,t961,t962,t963,t964,t965,t966,t967,t968,t969,t970,t971,t972,t973,t974,t975,t976,t977,t978,t979,t980,t981,t982,t983,t984,t985,t986,t987,t988,t989,t990,t991,t992,t993,t994,t995,t996,t997,t998,t999
        ,t1000,t1001,t1002,t1003,t1004,t1005,t1006,t1007,t1008,t1009,t1010,t1011,t1012,t1013,t1014,t1015,t1016,t1017,t1018,t1019,t1020,t1021,t1022,t1023,t1024,t1025,t1026,t1027,t1028,t1029,t1030,t1031,t1032,t1033,t1034,t1035,t1036,t1037,t1038,t1039,t1040,t1041,t1042,t1043,t1044,t1045,t1046,t1047,t1048,t1049
        ,t1050,t1051,t1052,t1053,t1054,t1055,t1056,t1057,t1058,t1059,t1060,t1061,t1062,t1063,t1064,t1065,t1066,t1067,t1068,t1069,t1070,t1071,t1072,t1073,t1074,t1075,t1076,t1077,t1078,t1079,t1080,t1081,t1082,t1083,t1084,t1085,t1086,t1087,t1088,t1089,t1090,t1091,t1092,t1093,t1094,t1095,t1096,t1097,t1098,t1099
        ,t1100,t1101,t1102,t1103,t1104,t1105,t1106,t1107,t1108,t1109,t1110,t1111,t1112,t1113,t1114,t1115,t1116,t1117,t1118,t1119,t1120,t1121,t1122,t1123,t1124,t1125,t1126,t1127,t1128,t1129,t1130,t1131,t1132,t1133,t1134,t1135,t1136,t1137,t1138,t1139,t1140,t1141,t1142,t1143,t1144,t1145,t1146,t1147,t1148,t1149
        ,t1150,t1151,t1152,t1153,t1154,t1155,t1156,t1157,t1158,t1159,t1160,t1161,t1162,t1163,t1164,t1165,t1166,t1167,t1168,t1169,t1170,t1171,t1172,t1173,t1174,t1175,t1176,t1177,t1178,t1179,t1180,t1181,t1182,t1183,t1184,t1185,t1186,t1187,t1188,t1189,t1190,t1191,t1192,t1193,t1194,t1195,t1196,t1197,t1198,t1199
        ,t1200;



            if( annulusScaling <= 0. )
            { // -- first time: compute the norm of the solution so we can scale to be size one
                Real t0=0.; 
                    Real ur,vr,pr; 
                    Real u1Max=0., u2Max=0., u3Max=0., pMax=0.; 
                    uNorm=0; 
                    if( numberOfDimensions==2 )
                    {
                        FOR_3D(i1,i2,i3,I1,I2,I3)
                        {
                            Real x = vertex(i1,i2,i3,0);
                            Real y = vertex(i1,i2,i3,1);
                            Real r = sqrt( x*x + y*y );
                            Real theta = atan2(y,x); 
                            Real cosTheta = x/r;
                            Real sinTheta = y/r;
                            t1 = n * n;
                            t2 = omega * omega;
                            t3 = t2 * t1;
                            t4 = omega / 0.2e1;
                            t5 = jn(n, t4);
                            t6 = pow(0.1e1 / 0.2e1, n);
                            t7 = t6 * t5;
                            t10 = n + 0.1e1;
                            t11 = jn(t10, t4);
                            t12 = t6 * t11;
                            t13 = t1 * n;
                            t14 = omega * t13;
                            t17 = omega * t11;
                            t18 = t6 * n;
                            t21 = 0.1e1 / t6;
                            t25 = t2 * t2;
                            t27 = t2 * omega;
                            t30 = t5 * t21;
                            t31 = t1 * t1;
                            t35 = t21 * t11;
                            t44 = jn(n, omega);
                            t47 = t2 * n;
                            t50 = t1 * t44;
                            t54 = -0.32e2 * t21 * n * t17 + 0.128e3 * t1 * t30 - 0.32e2 * t14 * t12 + 0.4e1 * t27 * t12 + 0.32e2 * t14 * t35 + 0.32e2 * t18 * t17 - 0.16e2 * t2 * t50 - t25 * t30 - t25 * t7 + 0.4e1 * t27 * t35 + 0.16e2 * t3 * t30 + 0.16e2 * t3 * t7 - 0.128e3 * t31 * t30 + 0.128e3 * t31 * t44 - 0.16e2 * t47 * t7 - 0.128e3 * t50;
                            t55 = 0.1e1 / t54;
                            t56 = t55 * c4;
                            t57 = pow(r, n);
                            t58 = t57 * t56;
                            t59 = t58 * t3;
                            t60 = 0.1e1 / r;
                            t61 = t21 * t60;
                            t62 = t5 * t61;
                            t63 = yn(n, omega);
                            t64 = omega * t0;
                            t65 = sin(t64);
                            t66 = t65 * t63;
                            t67 = n * theta;
                            t68 = sin(t67);
                            t69 = t68 * t66;
                            t73 = yn(n, t4);
                            t74 = t73 * t61;
                            t75 = t65 * t44;
                            t76 = t68 * t75;
                            t80 = cos(t64);
                            t81 = t80 * t63;
                            t82 = cos(t67);
                            t83 = t82 * t81;
                            t87 = t80 * t44;
                            t88 = t82 * t87;
                            t92 = t58 * t14;
                            t93 = yn(t10, t4);
                            t94 = t93 * t61;
                            t95 = t88 * t94;
                            t98 = n * omega;
                            t99 = t58 * t98;
                            t100 = t11 * t60;
                            t101 = t63 * t100;
                            t102 = t80 * t21;
                            t103 = t82 * t102;
                            t104 = t103 * t101;
                            t109 = 0.1e1 / t57;
                            t110 = t109 * t56;
                            t111 = t110 * t3;
                            t112 = t5 * t60;
                            t113 = t63 * t112;
                            t114 = t65 * t6;
                            t115 = t68 * t114;
                            t116 = t115 * t113;
                            t119 = t110 * t14;
                            t120 = t115 * t101;
                            t124 = t44 * t73 * t60;
                            t125 = t115 * t124;
                            t129 = t44 * t93 * t60;
                            t130 = t115 * t129;
                            t133 = t11 * t61;
                            t134 = t69 * t133;
                            t137 = t76 * t94;
                            t140 = -0.16e2 * t69 * t62 * t59 - 0.16e2 * t83 * t62 * t59 + 0.16e2 * t76 * t74 * t59 + 0.16e2 * t88 * t74 * t59 + 0.32e2 * t104 * t99 - 0.16e2 * t116 * t111 + 0.16e2 * t125 * t111 + 0.32e2 * t120 * t119 - 0.32e2 * t130 * t119 - 0.32e2 * t134 * t92 + 0.32e2 * t137 * t92 + 0.32e2 * t95 * t92 - 0.32e2 * t95 * t99;
                            t145 = t80 * t6;
                            t146 = t82 * t145;
                            t147 = t146 * t113;
                            t150 = t146 * t101;
                            t153 = t146 * t124;
                            t156 = t146 * t129;
                            t159 = t110 * t47;
                            t164 = t110 * t98;
                            t171 = t65 * t93;
                            t172 = t68 * t171;
                            t173 = t172 * t112;
                            t176 = t65 * t73;
                            t177 = t68 * t176;
                            t178 = t177 * t100;
                            t181 = t55 * t60;
                            t182 = omega * r;
                            t183 = jn(n, t182);
                            t184 = t183 * c4;
                            t186 = t73 * t184 * t181;
                            t187 = t1 * t6;
                            t188 = t65 * t2;
                            t189 = t68 * t188;
                            t193 = -0.16e2 * t189 * t187 * t186 - 0.32e2 * t104 * t92 - 0.16e2 * t147 * t111 + 0.16e2 * t153 * t111 + 0.32e2 * t150 * t119 - 0.32e2 * t156 * t119 + 0.32e2 * t173 * t119 - 0.32e2 * t178 * t119 + 0.32e2 * t134 * t99 - 0.32e2 * t137 * t99 + 0.16e2 * t147 * t159 - 0.32e2 * t150 * t164 - 0.16e2 * t153 * t159 + 0.32e2 * t156 * t164;
                            t196 = t93 * t184 * t181;
                            t197 = t13 * t6;
                            t198 = t65 * omega;
                            t199 = t68 * t198;
                            t206 = t25 * c4;
                            t207 = t109 * t55;
                            t208 = t60 * t207;
                            t209 = t208 * t206;
                            t210 = t44 * t73;
                            t211 = t115 * t210;
                            t213 = t27 * c4;
                            t214 = t208 * t213;
                            t215 = t63 * t11;
                            t216 = t115 * t215;
                            t219 = t44 * t93;
                            t220 = t115 * t219;
                            t223 = t57 * t55;
                            t224 = t60 * t223;
                            t225 = t224 * t213;
                            t226 = t65 * t21;
                            t227 = t68 * t226;
                            t228 = t227 * t215;
                            t231 = t227 * t219;
                            t235 = t224 * c4 * t1;
                            t236 = t69 * t30;
                            t239 = t184 * t181;
                            t240 = t31 * t63;
                            t241 = t68 * t65;
                            t245 = t1 * t63;
                            t249 = t82 * t80;
                            t256 = c4 * t60;
                            t257 = yn(n, t182);
                            t258 = t65 * t257;
                            t261 = t80 * t257;
                            t264 = 0.16e2 * t189 * t18 * t186 + 0.32e2 * t199 * t197 * t196 - 0.128e3 * t241 * t240 * t239 - 0.128e3 * t249 * t240 * t239 + 0.128e3 * t241 * t245 * t239 + 0.128e3 * t249 * t245 * t239 + t68 * t258 * t256 + t82 * t261 * t256 - t211 * t209 - 0.4e1 * t216 * t214 + 0.4e1 * t220 * t214 - 0.4e1 * t228 * t225 + 0.4e1 * t231 * t225 - 0.128e3 * t236 * t235;
                            t265 = t83 * t30;
                            t268 = t73 * t21;
                            t269 = t88 * t268;
                            t272 = t224 * t206;
                            t276 = t224 * c4 * t31;
                            t282 = t21 * t184 * t181;
                            t283 = t13 * t93;
                            t284 = t80 * omega;
                            t285 = t82 * t284;
                            t295 = t80 * t73;
                            t296 = t82 * t295;
                            t297 = t296 * t100;
                            t300 = t80 * t93;
                            t301 = t82 * t300;
                            t302 = t301 * t112;
                            t305 = t76 * t268;
                            t311 = 0.32e2 * t103 * t98 * t196 - 0.32e2 * t146 * t98 * t196 - 0.32e2 * t285 * t283 * t282 - 0.128e3 * t265 * t235 + 0.128e3 * t269 * t235 + 0.128e3 * t305 * t235 + 0.128e3 * t236 * t276 + t265 * t272 + 0.128e3 * t265 * t276 - t269 * t272 - 0.128e3 * t269 * t276 - t305 * t272 - 0.32e2 * t297 * t99 + 0.32e2 * t302 * t99;
                            t322 = t27 * n;
                            t323 = t110 * t322;
                            t326 = t1 * t73;
                            t342 = t80 * t2;
                            t343 = t82 * t342;
                            t349 = -0.32e2 * t115 * t98 * t196 - 0.16e2 * t189 * t326 * t282 - 0.32e2 * t199 * t283 * t282 - 0.16e2 * t343 * t326 * t282 + 0.16e2 * t116 * t159 - 0.32e2 * t120 * t164 - 0.16e2 * t125 * t159 + 0.32e2 * t130 * t164 + 0.4e1 * t178 * t323 - 0.128e3 * t305 * t276 + 0.4e1 * t297 * t323 - 0.4e1 * t302 * t323 - 0.32e2 * t302 * t92;
                            t361 = t63 * t5;
                            t362 = t115 * t361;
                            t364 = t146 * t361;
                            t366 = t146 * t210;
                            t368 = t146 * t215;
                            t371 = t146 * t219;
                            t374 = t103 * t215;
                            t377 = t103 * t219;
                            t380 = t58 * t322;
                            t383 = -0.32e2 * t297 * t119 + 0.32e2 * t302 * t119 + 0.32e2 * t297 * t164 - 0.32e2 * t302 * t164 + t362 * t209 + t364 * t209 - t366 * t209 - 0.4e1 * t368 * t214 + 0.4e1 * t371 * t214 - 0.4e1 * t374 * t225 + 0.4e1 * t377 * t225 + t236 * t272 + 0.32e2 * t297 * t92 + 0.4e1 * t302 * t380;
                            t417 = 0.16e2 * t343 * t18 * t186 - 0.16e2 * t343 * t187 * t186 + 0.32e2 * t285 * t197 * t196 + 0.32e2 * t227 * t98 * t196 - 0.32e2 * t173 * t164 + 0.32e2 * t178 * t164 - 0.4e1 * t173 * t323 + 0.4e1 * t173 * t380 - 0.32e2 * t173 * t92 + 0.32e2 * t173 * t99 - 0.4e1 * t178 * t380 + 0.32e2 * t178 * t92 - 0.32e2 * t178 * t99 - 0.4e1 * t297 * t380;
                            t418 = t21 * t93;
                            t419 = t65 * t27;
                            t420 = t68 * t419;
                            t424 = t65 * t1;
                            t429 = t6 * t93;
                            t433 = t6 * t73;
                            t434 = t65 * t25;
                            t435 = t68 * t434;
                            t440 = t80 * t25;
                            t441 = t82 * t440;
                            t444 = t80 * t27;
                            t445 = t82 * t444;
                            t455 = t80 * t31;
                            t465 = t80 * t1;
                            t470 = t65 * t31;
                            t475 = -0.128e3 * t68 * t424 * t268 * t239 + 0.128e3 * t82 * t455 * t268 * t239 - 0.128e3 * t82 * t465 * t268 * t239 + 0.128e3 * t68 * t470 * t268 * t239 + 0.16e2 * t189 * t245 * t239 + 0.16e2 * t343 * t245 * t239 + t435 * t268 * t239 + t441 * t268 * t239 - 0.4e1 * t420 * t418 * t239 - 0.4e1 * t445 * t418 * t239 - 0.4e1 * t420 * t429 * t239 - 0.4e1 * t445 * t429 * t239 + t435 * t433 * t239 + t441 * t433 * t239;
                            ur = t140 + t193 + t264 + t311 + t349 + t383 + t417 + t475;
                            t478 = t2 * t60;
                            t481 = c4 * t57 * t1 * t478;
                            t482 = t21 * t55;
                            t483 = t5 * t482;
                            t484 = t82 * t66;
                            t488 = t68 * t81;
                            t492 = t73 * t482;
                            t493 = t68 * t87;
                            t497 = t82 * t114;
                            t498 = t497 * t113;
                            t501 = t497 * t101;
                            t504 = t497 * t124;
                            t507 = t497 * t129;
                            t518 = t484 * t133;
                            t521 = t82 * t75;
                            t522 = t521 * t94;
                            t525 = 0.1e1 / n;
                            t527 = yn(t10, t182);
                            t528 = t527 * t525 * omega;
                            t535 = t68 * t284;
                            t539 = t68 * t295;
                            t540 = t539 * t100;
                            t543 = t68 * t300;
                            t544 = t543 * t112;
                            t547 = -t82 * t65 * c4 * t528 + t68 * t80 * c4 * t528 - 0.32e2 * t535 * t197 * t196 - 0.16e2 * t484 * t483 * t481 + 0.16e2 * t488 * t483 * t481 - 0.16e2 * t493 * t492 * t481 + 0.16e2 * t498 * t111 - 0.16e2 * t504 * t111 - 0.32e2 * t501 * t119 + 0.32e2 * t507 * t119 - 0.16e2 * t498 * t159 + 0.16e2 * t504 * t159 + 0.32e2 * t501 * t164 - 0.32e2 * t507 * t164 + 0.32e2 * t540 * t164 - 0.32e2 * t544 * t164 - 0.32e2 * t518 * t92 + 0.32e2 * t522 * t92;
                            t552 = t68 * t342;
                            t556 = t68 * t145;
                            t567 = t68 * t102;
                            t574 = t488 * t30;
                            t577 = t493 * t268;
                            t584 = t27 * t60;
                            t587 = c4 * t57 * n * t584;
                            t588 = t11 * t55;
                            t589 = t539 * t588;
                            t592 = jn(t10, t182);
                            t593 = t55 * t592;
                            t594 = t593 * t322;
                            t595 = t63 * c4;
                            t596 = t68 * t80;
                            t597 = t596 * t595;
                            t600 = t82 * t65;
                            t601 = t600 * t595;
                            t605 = t56 * t592 * t2;
                            t609 = -0.16e2 * t552 * t18 * t186 - t556 * t210 * t209 + t556 * t361 * t209 - 0.4e1 * t556 * t215 * t214 + 0.4e1 * t556 * t219 * t214 + 0.4e1 * t567 * t215 * t225 - 0.4e1 * t567 * t219 * t225 - 0.32e2 * t596 * t429 * t605 - 0.32e2 * t540 * t119 + 0.32e2 * t544 * t119 + 0.128e3 * t574 * t235 - 0.128e3 * t577 * t235 - t574 * t272 + t577 * t272 - 0.128e3 * t574 * t276 + 0.4e1 * t589 * t587 + 0.16e2 * t597 * t594 - 0.16e2 * t601 * t594;
                            t612 = t56 * t592 * t27;
                            t631 = t593 * t14;
                            t634 = t593 * t98;
                            t654 = t488 * t133;
                            t657 = t493 * t94;
                            t662 = 0.128e3 * t596 * t240 * t239 - 0.128e3 * t600 * t240 * t239 - 0.128e3 * t596 * t245 * t239 + 0.128e3 * t600 * t245 * t239 + 0.32e2 * t596 * t418 * t605 - 0.32e2 * t600 * t418 * t605 + 0.32e2 * t600 * t429 * t605 + 0.16e2 * t596 * t433 * t612 - 0.16e2 * t600 * t433 * t612 + 0.32e2 * t518 * t99 - 0.32e2 * t522 * t99 - 0.128e3 * t597 * t631 + 0.128e3 * t597 * t634 + 0.128e3 * t601 * t631 - 0.128e3 * t601 * t634 + 0.32e2 * t654 * t92 - 0.32e2 * t654 * t99 - 0.32e2 * t657 * t92;
                            t665 = t556 * t113;
                            t672 = t556 * t101;
                            t675 = t556 * t124;
                            t678 = t556 * t129;
                            t692 = t183 * t1;
                            t694 = t55 * t692 * t478;
                            t695 = t21 * c4;
                            t696 = t539 * t695;
                            t701 = c4 * t109 * n * t584;
                            t702 = t5 * t55;
                            t703 = t543 * t702;
                            t715 = -t497 * t361 * t209 + t82 * t258 * t256 - t68 * t261 * t256 + 0.32e2 * t535 * t283 * t282 + 0.16e2 * t521 * t492 * t481 - 0.16e2 * t665 * t111 + 0.16e2 * t675 * t111 + 0.32e2 * t672 * t119 - 0.32e2 * t678 * t119 + 0.16e2 * t665 * t159 - 0.16e2 * t675 * t159 - 0.32e2 * t672 * t164 + 0.32e2 * t678 * t164 + 0.128e3 * t577 * t276 + 0.4e1 * t589 * t701 + 0.32e2 * t657 * t99 + 0.16e2 * t696 * t694 - 0.4e1 * t703 * t701;
                            t732 = t82 * t176;
                            t733 = t732 * t588;
                            t736 = t82 * t171;
                            t737 = t736 * t702;
                            t742 = t732 * t695;
                            t749 = t736 * t112;
                            t752 = t732 * t100;
                            t764 = 0.16e2 * t552 * t187 * t186 + 0.32e2 * t556 * t98 * t196 - 0.32e2 * t567 * t98 * t196 - 0.32e2 * t749 * t119 + 0.32e2 * t749 * t164 - 0.32e2 * t752 * t164 - 0.32e2 * t540 * t92 + 0.32e2 * t540 * t99 + 0.32e2 * t544 * t92 - 0.32e2 * t544 * t99 - 0.4e1 * t703 * t587 - 0.4e1 * t733 * t587 + 0.4e1 * t737 * t587 - 0.16e2 * t742 * t694 - 0.4e1 * t733 * t701 + 0.4e1 * t737 * t701 - 0.32e2 * t749 * t92 + 0.32e2 * t752 * t92;
                            t767 = t82 * t188;
                            t771 = t82 * t198;
                            t786 = t82 * t226;
                            t793 = t484 * t30;
                            t796 = t521 * t268;
                            t814 = 0.16e2 * t767 * t18 * t186 - 0.16e2 * t767 * t187 * t186 + 0.32e2 * t771 * t197 * t196 - 0.32e2 * t497 * t98 * t196 + 0.32e2 * t786 * t98 * t196 + t497 * t210 * t209 + 0.4e1 * t497 * t215 * t214 - 0.4e1 * t497 * t219 * t214 - 0.4e1 * t786 * t215 * t225 + 0.4e1 * t786 * t219 * t225 - 0.32e2 * t771 * t283 * t282 + 0.32e2 * t752 * t119 - 0.128e3 * t793 * t235 + 0.128e3 * t796 * t235 + t793 * t272 - t796 * t272 + 0.128e3 * t793 * t276 - 0.128e3 * t796 * t276;
                            t820 = t692 * t478;
                            t828 = t593 * t3;
                            t829 = t736 * t695;
                            t832 = t25 * omega;
                            t834 = t593 * t525 * t832;
                            t835 = t73 * c4;
                            t836 = t497 * t835;
                            t839 = t593 * t525 * t25;
                            t840 = t93 * c4;
                            t841 = t497 * t840;
                            t848 = t556 * t840;
                            t851 = t567 * t840;
                            t864 = -0.16e2 * t488 * t56 * t820 - 0.16e2 * t696 * t594 + 0.16e2 * t742 * t594 + 0.128e3 * t696 * t631 - 0.128e3 * t742 * t631 - 0.128e3 * t696 * t634 + 0.128e3 * t742 * t634 + t696 * t834 - t742 * t834 + 0.32e2 * t749 * t99 - 0.32e2 * t752 * t99 + 0.32e2 * t829 * t828 + 0.32e2 * t848 * t828 + 0.4e1 * t829 * t839 - t836 * t834 + 0.4e1 * t841 * t839 - 0.4e1 * t848 * t839 - 0.4e1 * t851 * t839;
                            t865 = t556 * t835;
                            t878 = t68 * t440;
                            t881 = t68 * t444;
                            t902 = t82 * t419;
                            t910 = t82 * t434;
                            t918 = -0.128e3 * t82 * t424 * t268 * t239 - 0.128e3 * t68 * t455 * t268 * t239 + 0.128e3 * t68 * t465 * t268 * t239 + 0.128e3 * t82 * t470 * t268 * t239 - t878 * t268 * t239 + t910 * t268 * t239 + 0.4e1 * t881 * t418 * t239 - 0.4e1 * t902 * t418 * t239 + 0.4e1 * t881 * t429 * t239 - 0.4e1 * t902 * t429 * t239 - t878 * t433 * t239 + t910 * t433 * t239 + 0.16e2 * t484 * t56 * t820 + 0.16e2 * t836 * t594 - 0.16e2 * t865 * t594 - 0.32e2 * t841 * t828 - 0.32e2 * t851 * t828 + t865 * t834;
                            vr = t547 + t609 + t662 + t715 + t764 + t814 + t864 + t918;
                            t921 = t832 * c4;
                            t922 = t525 * t55;
                            t923 = t57 * t922;
                            t924 = t923 * t921;
                            t927 = t1 * t55;
                            t929 = t57 * t927 * t213;
                            t932 = n * t55;
                            t933 = t57 * t932;
                            t934 = t933 * t206;
                            t938 = t109 * t927 * t213;
                            t948 = t109 * t932 * t206;
                            t958 = t25 * t2 * c4;
                            t959 = t109 * t922;
                            t960 = t959 * t958;
                            t964 = t2 * c4;
                            t965 = t933 * t964;
                            t968 = t923 * t958;
                            t971 = t959 * t921;
                            t976 = t223 * t921;
                            t977 = t73 * t11;
                            t978 = t249 * t977;
                            t981 = -0.32e2 * t228 * t929 + 0.4e1 * t231 * t924 + t236 * t968 - 0.16e2 * t265 * t934 - 0.128e3 * t265 * t965 + 0.16e2 * t269 * t934 + 0.16e2 * t305 * t934 + 0.16e2 * t364 * t948 - t364 * t960 - 0.16e2 * t366 * t948 + t366 * t960 - 0.32e2 * t368 * t938 + 0.4e1 * t368 * t971 + 0.32e2 * t371 * t938 - 0.4e1 * t371 * t971 - 0.32e2 * t374 * t929 + 0.4e1 * t377 * t924 + 0.32e2 * t377 * t929 - 0.4e1 * t978 * t976;
                            t982 = t207 * t213;
                            t985 = t93 * t5;
                            t986 = t249 * t985;
                            t989 = t207 * t921;
                            t994 = t223 * t213;
                            t999 = t241 * t985;
                            t1002 = t241 * t977;
                            t1025 = t927 * t213;
                            t1026 = t5 * t57;
                            t1030 = -0.32e2 * t301 * t1026 * t1025 - 0.4e1 * t1002 * t976 - 0.32e2 * t1002 * t982 - 0.4e1 * t1002 * t989 - 0.32e2 * t1002 * t994 - 0.32e2 * t371 * t982 + 0.32e2 * t374 * t994 - 0.32e2 * t377 * t994 + 0.4e1 * t986 * t976 + 0.4e1 * t999 * t976 - 0.32e2 * t978 * t982 - 0.4e1 * t978 * t989 - 0.32e2 * t978 * t994 + 0.32e2 * t986 * t982 + 0.32e2 * t999 * t982 + 0.4e1 * t986 * t989 + 0.32e2 * t986 * t994 + 0.4e1 * t999 * t989 + 0.32e2 * t999 * t994;
                            t1032 = t11 * t57;
                            t1036 = t207 * t206;
                            t1045 = t5 * t109;
                            t1049 = t11 * t109;
                            t1077 = -0.32e2 * t172 * t1026 * t1025 + 0.32e2 * t296 * t1032 * t1025 - 0.32e2 * t301 * t1045 * t1025 + 0.32e2 * t296 * t1049 * t1025 + 0.16e2 * t211 * t1036 - 0.16e2 * t362 * t1036 + 0.4e1 * t216 * t971 + 0.32e2 * t216 * t982 + 0.32e2 * t220 * t938 - 0.4e1 * t220 * t971 - 0.32e2 * t220 * t982 + 0.32e2 * t228 * t994 + 0.32e2 * t231 * t929 - 0.32e2 * t231 * t994 - 0.16e2 * t236 * t934 - 0.128e3 * t236 * t965 + 0.128e3 * t305 * t965 + 0.16e2 * t362 * t948 - t362 * t960;
                            t1081 = t57 * t13 * t55 * t964;
                            t1118 = 0.32e2 * t177 * t1032 * t1025 - 0.32e2 * t172 * t1045 * t1025 + 0.32e2 * t177 * t1049 * t1025 - 0.16e2 * t364 * t1036 + 0.16e2 * t366 * t1036 + 0.128e3 * t236 * t1081 + 0.128e3 * t265 * t1081 - 0.128e3 * t269 * t1081 - 0.128e3 * t305 * t1081 - 0.16e2 * t211 * t948 + t211 * t960 - 0.32e2 * t216 * t938 - 0.4e1 * t228 * t924 + t265 * t968 + 0.128e3 * t269 * t965 - t269 * t968 - t305 * t968 + 0.32e2 * t368 * t982 - 0.4e1 * t374 * t924;
                            pr = t981 + t1030 + t1077 + t1118;
              // (u1,u2) = ur*rHat + vr*thetaHat 
                            u(i1,i2,i3,u1c) = ur*cosTheta - vr*sinTheta; 
                            u(i1,i2,i3,u2c) = ur*sinTheta + vr*cosTheta;
                            u(i1,i2,i3,pc ) = pr;
                            u1Max=max(u1Max,u(i1,i2,i3,u1c));
                            u2Max=max(u2Max,u(i1,i2,i3,u2c));
                            pMax =max(pMax, u(i1,i2,i3,pc));
                        }
                        uNorm = max(u1Max,u2Max,u3Max,pMax);
            // if( scale==0. )
            // {
            //   printF("UDKS: ERROR: annulus solution is zero! scale=%e\n",scale);
            //   OV_ABORT("ERROR");
            // }
            // scale = 1./scale;
            // c4 *= scale;  // set this so the velocity and acceleration are scaled too
            // FOR_3(i1,i2,i3,I1,I2,I3)
            // {
            //   u(i1,i2,i3,u1c) *= scale;
            //   u(i1,i2,i3,u2c) *= scale;
            //   u(i1,i2,i3,pc ) *= scale;
            // }    
                    }
                    else
                    {
                        OV_ABORT("UDKS:ERROR: annulus solution : 3D called");
                    }
                if( uNorm==0. )
                {
                    printF("UDKS: ERROR: annulus solution is zero! uNorm=%e\n",uNorm);
                    OV_ABORT("ERROR");
                }        
                annulusScaling = c4/uNorm;
            }
            c4 = annulusScaling; 

                Real ur,vr,pr; 
                Real u1Max=0., u2Max=0., u3Max=0., pMax=0.; 
                uNorm=0; 
                if( numberOfDimensions==2 )
                {
                    FOR_3D(i1,i2,i3,I1,I2,I3)
                    {
                        Real x = vertex(i1,i2,i3,0);
                        Real y = vertex(i1,i2,i3,1);
                        Real r = sqrt( x*x + y*y );
                        Real theta = atan2(y,x); 
                        Real cosTheta = x/r;
                        Real sinTheta = y/r;
                        t1 = n * n;
                        t2 = omega * omega;
                        t3 = t2 * t1;
                        t4 = omega / 0.2e1;
                        t5 = jn(n, t4);
                        t6 = pow(0.1e1 / 0.2e1, n);
                        t7 = t6 * t5;
                        t10 = n + 0.1e1;
                        t11 = jn(t10, t4);
                        t12 = t6 * t11;
                        t13 = t1 * n;
                        t14 = omega * t13;
                        t17 = omega * t11;
                        t18 = t6 * n;
                        t21 = 0.1e1 / t6;
                        t25 = t2 * t2;
                        t27 = t2 * omega;
                        t30 = t5 * t21;
                        t31 = t1 * t1;
                        t35 = t21 * t11;
                        t44 = jn(n, omega);
                        t47 = t2 * n;
                        t50 = t1 * t44;
                        t54 = -0.32e2 * t21 * n * t17 + 0.128e3 * t1 * t30 - 0.32e2 * t14 * t12 + 0.4e1 * t27 * t12 + 0.32e2 * t14 * t35 + 0.32e2 * t18 * t17 - 0.16e2 * t2 * t50 - t25 * t30 - t25 * t7 + 0.4e1 * t27 * t35 + 0.16e2 * t3 * t30 + 0.16e2 * t3 * t7 - 0.128e3 * t31 * t30 + 0.128e3 * t31 * t44 - 0.16e2 * t47 * t7 - 0.128e3 * t50;
                        t55 = 0.1e1 / t54;
                        t56 = t55 * c4;
                        t57 = pow(r, n);
                        t58 = t57 * t56;
                        t59 = t58 * t3;
                        t60 = 0.1e1 / r;
                        t61 = t21 * t60;
                        t62 = t5 * t61;
                        t63 = yn(n, omega);
                        t64 = omega * t;
                        t65 = sin(t64);
                        t66 = t65 * t63;
                        t67 = n * theta;
                        t68 = sin(t67);
                        t69 = t68 * t66;
                        t73 = yn(n, t4);
                        t74 = t73 * t61;
                        t75 = t65 * t44;
                        t76 = t68 * t75;
                        t80 = cos(t64);
                        t81 = t80 * t63;
                        t82 = cos(t67);
                        t83 = t82 * t81;
                        t87 = t80 * t44;
                        t88 = t82 * t87;
                        t92 = t58 * t14;
                        t93 = yn(t10, t4);
                        t94 = t93 * t61;
                        t95 = t88 * t94;
                        t98 = n * omega;
                        t99 = t58 * t98;
                        t100 = t11 * t60;
                        t101 = t63 * t100;
                        t102 = t80 * t21;
                        t103 = t82 * t102;
                        t104 = t103 * t101;
                        t109 = 0.1e1 / t57;
                        t110 = t109 * t56;
                        t111 = t110 * t3;
                        t112 = t5 * t60;
                        t113 = t63 * t112;
                        t114 = t65 * t6;
                        t115 = t68 * t114;
                        t116 = t115 * t113;
                        t119 = t110 * t14;
                        t120 = t115 * t101;
                        t124 = t44 * t73 * t60;
                        t125 = t115 * t124;
                        t129 = t44 * t93 * t60;
                        t130 = t115 * t129;
                        t133 = t11 * t61;
                        t134 = t69 * t133;
                        t137 = t76 * t94;
                        t140 = -0.16e2 * t69 * t62 * t59 - 0.16e2 * t83 * t62 * t59 + 0.16e2 * t76 * t74 * t59 + 0.16e2 * t88 * t74 * t59 + 0.32e2 * t104 * t99 - 0.16e2 * t116 * t111 + 0.16e2 * t125 * t111 + 0.32e2 * t120 * t119 - 0.32e2 * t130 * t119 - 0.32e2 * t134 * t92 + 0.32e2 * t137 * t92 + 0.32e2 * t95 * t92 - 0.32e2 * t95 * t99;
                        t145 = t80 * t6;
                        t146 = t82 * t145;
                        t147 = t146 * t113;
                        t150 = t146 * t101;
                        t153 = t146 * t124;
                        t156 = t146 * t129;
                        t159 = t110 * t47;
                        t164 = t110 * t98;
                        t171 = t65 * t93;
                        t172 = t68 * t171;
                        t173 = t172 * t112;
                        t176 = t65 * t73;
                        t177 = t68 * t176;
                        t178 = t177 * t100;
                        t181 = t55 * t60;
                        t182 = omega * r;
                        t183 = jn(n, t182);
                        t184 = t183 * c4;
                        t186 = t73 * t184 * t181;
                        t187 = t1 * t6;
                        t188 = t65 * t2;
                        t189 = t68 * t188;
                        t193 = -0.16e2 * t189 * t187 * t186 - 0.32e2 * t104 * t92 - 0.16e2 * t147 * t111 + 0.16e2 * t153 * t111 + 0.32e2 * t150 * t119 - 0.32e2 * t156 * t119 + 0.32e2 * t173 * t119 - 0.32e2 * t178 * t119 + 0.32e2 * t134 * t99 - 0.32e2 * t137 * t99 + 0.16e2 * t147 * t159 - 0.32e2 * t150 * t164 - 0.16e2 * t153 * t159 + 0.32e2 * t156 * t164;
                        t196 = t93 * t184 * t181;
                        t197 = t13 * t6;
                        t198 = t65 * omega;
                        t199 = t68 * t198;
                        t206 = t25 * c4;
                        t207 = t109 * t55;
                        t208 = t60 * t207;
                        t209 = t208 * t206;
                        t210 = t44 * t73;
                        t211 = t115 * t210;
                        t213 = t27 * c4;
                        t214 = t208 * t213;
                        t215 = t63 * t11;
                        t216 = t115 * t215;
                        t219 = t44 * t93;
                        t220 = t115 * t219;
                        t223 = t57 * t55;
                        t224 = t60 * t223;
                        t225 = t224 * t213;
                        t226 = t65 * t21;
                        t227 = t68 * t226;
                        t228 = t227 * t215;
                        t231 = t227 * t219;
                        t235 = t224 * c4 * t1;
                        t236 = t69 * t30;
                        t239 = t184 * t181;
                        t240 = t31 * t63;
                        t241 = t68 * t65;
                        t245 = t1 * t63;
                        t249 = t82 * t80;
                        t256 = c4 * t60;
                        t257 = yn(n, t182);
                        t258 = t65 * t257;
                        t261 = t80 * t257;
                        t264 = 0.16e2 * t189 * t18 * t186 + 0.32e2 * t199 * t197 * t196 - 0.128e3 * t241 * t240 * t239 - 0.128e3 * t249 * t240 * t239 + 0.128e3 * t241 * t245 * t239 + 0.128e3 * t249 * t245 * t239 + t68 * t258 * t256 + t82 * t261 * t256 - t211 * t209 - 0.4e1 * t216 * t214 + 0.4e1 * t220 * t214 - 0.4e1 * t228 * t225 + 0.4e1 * t231 * t225 - 0.128e3 * t236 * t235;
                        t265 = t83 * t30;
                        t268 = t73 * t21;
                        t269 = t88 * t268;
                        t272 = t224 * t206;
                        t276 = t224 * c4 * t31;
                        t282 = t21 * t184 * t181;
                        t283 = t13 * t93;
                        t284 = t80 * omega;
                        t285 = t82 * t284;
                        t295 = t80 * t73;
                        t296 = t82 * t295;
                        t297 = t296 * t100;
                        t300 = t80 * t93;
                        t301 = t82 * t300;
                        t302 = t301 * t112;
                        t305 = t76 * t268;
                        t311 = 0.32e2 * t103 * t98 * t196 - 0.32e2 * t146 * t98 * t196 - 0.32e2 * t285 * t283 * t282 - 0.128e3 * t265 * t235 + 0.128e3 * t269 * t235 + 0.128e3 * t305 * t235 + 0.128e3 * t236 * t276 + t265 * t272 + 0.128e3 * t265 * t276 - t269 * t272 - 0.128e3 * t269 * t276 - t305 * t272 - 0.32e2 * t297 * t99 + 0.32e2 * t302 * t99;
                        t322 = t27 * n;
                        t323 = t110 * t322;
                        t326 = t1 * t73;
                        t342 = t80 * t2;
                        t343 = t82 * t342;
                        t349 = -0.32e2 * t115 * t98 * t196 - 0.16e2 * t189 * t326 * t282 - 0.32e2 * t199 * t283 * t282 - 0.16e2 * t343 * t326 * t282 + 0.16e2 * t116 * t159 - 0.32e2 * t120 * t164 - 0.16e2 * t125 * t159 + 0.32e2 * t130 * t164 + 0.4e1 * t178 * t323 - 0.128e3 * t305 * t276 + 0.4e1 * t297 * t323 - 0.4e1 * t302 * t323 - 0.32e2 * t302 * t92;
                        t361 = t63 * t5;
                        t362 = t115 * t361;
                        t364 = t146 * t361;
                        t366 = t146 * t210;
                        t368 = t146 * t215;
                        t371 = t146 * t219;
                        t374 = t103 * t215;
                        t377 = t103 * t219;
                        t380 = t58 * t322;
                        t383 = -0.32e2 * t297 * t119 + 0.32e2 * t302 * t119 + 0.32e2 * t297 * t164 - 0.32e2 * t302 * t164 + t362 * t209 + t364 * t209 - t366 * t209 - 0.4e1 * t368 * t214 + 0.4e1 * t371 * t214 - 0.4e1 * t374 * t225 + 0.4e1 * t377 * t225 + t236 * t272 + 0.32e2 * t297 * t92 + 0.4e1 * t302 * t380;
                        t417 = 0.16e2 * t343 * t18 * t186 - 0.16e2 * t343 * t187 * t186 + 0.32e2 * t285 * t197 * t196 + 0.32e2 * t227 * t98 * t196 - 0.32e2 * t173 * t164 + 0.32e2 * t178 * t164 - 0.4e1 * t173 * t323 + 0.4e1 * t173 * t380 - 0.32e2 * t173 * t92 + 0.32e2 * t173 * t99 - 0.4e1 * t178 * t380 + 0.32e2 * t178 * t92 - 0.32e2 * t178 * t99 - 0.4e1 * t297 * t380;
                        t418 = t21 * t93;
                        t419 = t65 * t27;
                        t420 = t68 * t419;
                        t424 = t65 * t1;
                        t429 = t6 * t93;
                        t433 = t6 * t73;
                        t434 = t65 * t25;
                        t435 = t68 * t434;
                        t440 = t80 * t25;
                        t441 = t82 * t440;
                        t444 = t80 * t27;
                        t445 = t82 * t444;
                        t455 = t80 * t31;
                        t465 = t80 * t1;
                        t470 = t65 * t31;
                        t475 = -0.128e3 * t68 * t424 * t268 * t239 + 0.128e3 * t82 * t455 * t268 * t239 - 0.128e3 * t82 * t465 * t268 * t239 + 0.128e3 * t68 * t470 * t268 * t239 + 0.16e2 * t189 * t245 * t239 + 0.16e2 * t343 * t245 * t239 + t435 * t268 * t239 + t441 * t268 * t239 - 0.4e1 * t420 * t418 * t239 - 0.4e1 * t445 * t418 * t239 - 0.4e1 * t420 * t429 * t239 - 0.4e1 * t445 * t429 * t239 + t435 * t433 * t239 + t441 * t433 * t239;
                        ur = t140 + t193 + t264 + t311 + t349 + t383 + t417 + t475;
                        t478 = t2 * t60;
                        t481 = c4 * t57 * t1 * t478;
                        t482 = t21 * t55;
                        t483 = t5 * t482;
                        t484 = t82 * t66;
                        t488 = t68 * t81;
                        t492 = t73 * t482;
                        t493 = t68 * t87;
                        t497 = t82 * t114;
                        t498 = t497 * t113;
                        t501 = t497 * t101;
                        t504 = t497 * t124;
                        t507 = t497 * t129;
                        t518 = t484 * t133;
                        t521 = t82 * t75;
                        t522 = t521 * t94;
                        t525 = 0.1e1 / n;
                        t527 = yn(t10, t182);
                        t528 = t527 * t525 * omega;
                        t535 = t68 * t284;
                        t539 = t68 * t295;
                        t540 = t539 * t100;
                        t543 = t68 * t300;
                        t544 = t543 * t112;
                        t547 = -t82 * t65 * c4 * t528 + t68 * t80 * c4 * t528 - 0.32e2 * t535 * t197 * t196 - 0.16e2 * t484 * t483 * t481 + 0.16e2 * t488 * t483 * t481 - 0.16e2 * t493 * t492 * t481 + 0.16e2 * t498 * t111 - 0.16e2 * t504 * t111 - 0.32e2 * t501 * t119 + 0.32e2 * t507 * t119 - 0.16e2 * t498 * t159 + 0.16e2 * t504 * t159 + 0.32e2 * t501 * t164 - 0.32e2 * t507 * t164 + 0.32e2 * t540 * t164 - 0.32e2 * t544 * t164 - 0.32e2 * t518 * t92 + 0.32e2 * t522 * t92;
                        t552 = t68 * t342;
                        t556 = t68 * t145;
                        t567 = t68 * t102;
                        t574 = t488 * t30;
                        t577 = t493 * t268;
                        t584 = t27 * t60;
                        t587 = c4 * t57 * n * t584;
                        t588 = t11 * t55;
                        t589 = t539 * t588;
                        t592 = jn(t10, t182);
                        t593 = t55 * t592;
                        t594 = t593 * t322;
                        t595 = t63 * c4;
                        t596 = t68 * t80;
                        t597 = t596 * t595;
                        t600 = t82 * t65;
                        t601 = t600 * t595;
                        t605 = t56 * t592 * t2;
                        t609 = -0.16e2 * t552 * t18 * t186 - t556 * t210 * t209 + t556 * t361 * t209 - 0.4e1 * t556 * t215 * t214 + 0.4e1 * t556 * t219 * t214 + 0.4e1 * t567 * t215 * t225 - 0.4e1 * t567 * t219 * t225 - 0.32e2 * t596 * t429 * t605 - 0.32e2 * t540 * t119 + 0.32e2 * t544 * t119 + 0.128e3 * t574 * t235 - 0.128e3 * t577 * t235 - t574 * t272 + t577 * t272 - 0.128e3 * t574 * t276 + 0.4e1 * t589 * t587 + 0.16e2 * t597 * t594 - 0.16e2 * t601 * t594;
                        t612 = t56 * t592 * t27;
                        t631 = t593 * t14;
                        t634 = t593 * t98;
                        t654 = t488 * t133;
                        t657 = t493 * t94;
                        t662 = 0.128e3 * t596 * t240 * t239 - 0.128e3 * t600 * t240 * t239 - 0.128e3 * t596 * t245 * t239 + 0.128e3 * t600 * t245 * t239 + 0.32e2 * t596 * t418 * t605 - 0.32e2 * t600 * t418 * t605 + 0.32e2 * t600 * t429 * t605 + 0.16e2 * t596 * t433 * t612 - 0.16e2 * t600 * t433 * t612 + 0.32e2 * t518 * t99 - 0.32e2 * t522 * t99 - 0.128e3 * t597 * t631 + 0.128e3 * t597 * t634 + 0.128e3 * t601 * t631 - 0.128e3 * t601 * t634 + 0.32e2 * t654 * t92 - 0.32e2 * t654 * t99 - 0.32e2 * t657 * t92;
                        t665 = t556 * t113;
                        t672 = t556 * t101;
                        t675 = t556 * t124;
                        t678 = t556 * t129;
                        t692 = t183 * t1;
                        t694 = t55 * t692 * t478;
                        t695 = t21 * c4;
                        t696 = t539 * t695;
                        t701 = c4 * t109 * n * t584;
                        t702 = t5 * t55;
                        t703 = t543 * t702;
                        t715 = -t497 * t361 * t209 + t82 * t258 * t256 - t68 * t261 * t256 + 0.32e2 * t535 * t283 * t282 + 0.16e2 * t521 * t492 * t481 - 0.16e2 * t665 * t111 + 0.16e2 * t675 * t111 + 0.32e2 * t672 * t119 - 0.32e2 * t678 * t119 + 0.16e2 * t665 * t159 - 0.16e2 * t675 * t159 - 0.32e2 * t672 * t164 + 0.32e2 * t678 * t164 + 0.128e3 * t577 * t276 + 0.4e1 * t589 * t701 + 0.32e2 * t657 * t99 + 0.16e2 * t696 * t694 - 0.4e1 * t703 * t701;
                        t732 = t82 * t176;
                        t733 = t732 * t588;
                        t736 = t82 * t171;
                        t737 = t736 * t702;
                        t742 = t732 * t695;
                        t749 = t736 * t112;
                        t752 = t732 * t100;
                        t764 = 0.16e2 * t552 * t187 * t186 + 0.32e2 * t556 * t98 * t196 - 0.32e2 * t567 * t98 * t196 - 0.32e2 * t749 * t119 + 0.32e2 * t749 * t164 - 0.32e2 * t752 * t164 - 0.32e2 * t540 * t92 + 0.32e2 * t540 * t99 + 0.32e2 * t544 * t92 - 0.32e2 * t544 * t99 - 0.4e1 * t703 * t587 - 0.4e1 * t733 * t587 + 0.4e1 * t737 * t587 - 0.16e2 * t742 * t694 - 0.4e1 * t733 * t701 + 0.4e1 * t737 * t701 - 0.32e2 * t749 * t92 + 0.32e2 * t752 * t92;
                        t767 = t82 * t188;
                        t771 = t82 * t198;
                        t786 = t82 * t226;
                        t793 = t484 * t30;
                        t796 = t521 * t268;
                        t814 = 0.16e2 * t767 * t18 * t186 - 0.16e2 * t767 * t187 * t186 + 0.32e2 * t771 * t197 * t196 - 0.32e2 * t497 * t98 * t196 + 0.32e2 * t786 * t98 * t196 + t497 * t210 * t209 + 0.4e1 * t497 * t215 * t214 - 0.4e1 * t497 * t219 * t214 - 0.4e1 * t786 * t215 * t225 + 0.4e1 * t786 * t219 * t225 - 0.32e2 * t771 * t283 * t282 + 0.32e2 * t752 * t119 - 0.128e3 * t793 * t235 + 0.128e3 * t796 * t235 + t793 * t272 - t796 * t272 + 0.128e3 * t793 * t276 - 0.128e3 * t796 * t276;
                        t820 = t692 * t478;
                        t828 = t593 * t3;
                        t829 = t736 * t695;
                        t832 = t25 * omega;
                        t834 = t593 * t525 * t832;
                        t835 = t73 * c4;
                        t836 = t497 * t835;
                        t839 = t593 * t525 * t25;
                        t840 = t93 * c4;
                        t841 = t497 * t840;
                        t848 = t556 * t840;
                        t851 = t567 * t840;
                        t864 = -0.16e2 * t488 * t56 * t820 - 0.16e2 * t696 * t594 + 0.16e2 * t742 * t594 + 0.128e3 * t696 * t631 - 0.128e3 * t742 * t631 - 0.128e3 * t696 * t634 + 0.128e3 * t742 * t634 + t696 * t834 - t742 * t834 + 0.32e2 * t749 * t99 - 0.32e2 * t752 * t99 + 0.32e2 * t829 * t828 + 0.32e2 * t848 * t828 + 0.4e1 * t829 * t839 - t836 * t834 + 0.4e1 * t841 * t839 - 0.4e1 * t848 * t839 - 0.4e1 * t851 * t839;
                        t865 = t556 * t835;
                        t878 = t68 * t440;
                        t881 = t68 * t444;
                        t902 = t82 * t419;
                        t910 = t82 * t434;
                        t918 = -0.128e3 * t82 * t424 * t268 * t239 - 0.128e3 * t68 * t455 * t268 * t239 + 0.128e3 * t68 * t465 * t268 * t239 + 0.128e3 * t82 * t470 * t268 * t239 - t878 * t268 * t239 + t910 * t268 * t239 + 0.4e1 * t881 * t418 * t239 - 0.4e1 * t902 * t418 * t239 + 0.4e1 * t881 * t429 * t239 - 0.4e1 * t902 * t429 * t239 - t878 * t433 * t239 + t910 * t433 * t239 + 0.16e2 * t484 * t56 * t820 + 0.16e2 * t836 * t594 - 0.16e2 * t865 * t594 - 0.32e2 * t841 * t828 - 0.32e2 * t851 * t828 + t865 * t834;
                        vr = t547 + t609 + t662 + t715 + t764 + t814 + t864 + t918;
                        t921 = t832 * c4;
                        t922 = t525 * t55;
                        t923 = t57 * t922;
                        t924 = t923 * t921;
                        t927 = t1 * t55;
                        t929 = t57 * t927 * t213;
                        t932 = n * t55;
                        t933 = t57 * t932;
                        t934 = t933 * t206;
                        t938 = t109 * t927 * t213;
                        t948 = t109 * t932 * t206;
                        t958 = t25 * t2 * c4;
                        t959 = t109 * t922;
                        t960 = t959 * t958;
                        t964 = t2 * c4;
                        t965 = t933 * t964;
                        t968 = t923 * t958;
                        t971 = t959 * t921;
                        t976 = t223 * t921;
                        t977 = t73 * t11;
                        t978 = t249 * t977;
                        t981 = -0.32e2 * t228 * t929 + 0.4e1 * t231 * t924 + t236 * t968 - 0.16e2 * t265 * t934 - 0.128e3 * t265 * t965 + 0.16e2 * t269 * t934 + 0.16e2 * t305 * t934 + 0.16e2 * t364 * t948 - t364 * t960 - 0.16e2 * t366 * t948 + t366 * t960 - 0.32e2 * t368 * t938 + 0.4e1 * t368 * t971 + 0.32e2 * t371 * t938 - 0.4e1 * t371 * t971 - 0.32e2 * t374 * t929 + 0.4e1 * t377 * t924 + 0.32e2 * t377 * t929 - 0.4e1 * t978 * t976;
                        t982 = t207 * t213;
                        t985 = t93 * t5;
                        t986 = t249 * t985;
                        t989 = t207 * t921;
                        t994 = t223 * t213;
                        t999 = t241 * t985;
                        t1002 = t241 * t977;
                        t1025 = t927 * t213;
                        t1026 = t5 * t57;
                        t1030 = -0.32e2 * t301 * t1026 * t1025 - 0.4e1 * t1002 * t976 - 0.32e2 * t1002 * t982 - 0.4e1 * t1002 * t989 - 0.32e2 * t1002 * t994 - 0.32e2 * t371 * t982 + 0.32e2 * t374 * t994 - 0.32e2 * t377 * t994 + 0.4e1 * t986 * t976 + 0.4e1 * t999 * t976 - 0.32e2 * t978 * t982 - 0.4e1 * t978 * t989 - 0.32e2 * t978 * t994 + 0.32e2 * t986 * t982 + 0.32e2 * t999 * t982 + 0.4e1 * t986 * t989 + 0.32e2 * t986 * t994 + 0.4e1 * t999 * t989 + 0.32e2 * t999 * t994;
                        t1032 = t11 * t57;
                        t1036 = t207 * t206;
                        t1045 = t5 * t109;
                        t1049 = t11 * t109;
                        t1077 = -0.32e2 * t172 * t1026 * t1025 + 0.32e2 * t296 * t1032 * t1025 - 0.32e2 * t301 * t1045 * t1025 + 0.32e2 * t296 * t1049 * t1025 + 0.16e2 * t211 * t1036 - 0.16e2 * t362 * t1036 + 0.4e1 * t216 * t971 + 0.32e2 * t216 * t982 + 0.32e2 * t220 * t938 - 0.4e1 * t220 * t971 - 0.32e2 * t220 * t982 + 0.32e2 * t228 * t994 + 0.32e2 * t231 * t929 - 0.32e2 * t231 * t994 - 0.16e2 * t236 * t934 - 0.128e3 * t236 * t965 + 0.128e3 * t305 * t965 + 0.16e2 * t362 * t948 - t362 * t960;
                        t1081 = t57 * t13 * t55 * t964;
                        t1118 = 0.32e2 * t177 * t1032 * t1025 - 0.32e2 * t172 * t1045 * t1025 + 0.32e2 * t177 * t1049 * t1025 - 0.16e2 * t364 * t1036 + 0.16e2 * t366 * t1036 + 0.128e3 * t236 * t1081 + 0.128e3 * t265 * t1081 - 0.128e3 * t269 * t1081 - 0.128e3 * t305 * t1081 - 0.16e2 * t211 * t948 + t211 * t960 - 0.32e2 * t216 * t938 - 0.4e1 * t228 * t924 + t265 * t968 + 0.128e3 * t269 * t965 - t269 * t968 - t305 * t968 + 0.32e2 * t368 * t982 - 0.4e1 * t374 * t924;
                        pr = t981 + t1030 + t1077 + t1118;
            // (u1,u2) = ur*rHat + vr*thetaHat 
                        u(i1,i2,i3,u1c) = ur*cosTheta - vr*sinTheta; 
                        u(i1,i2,i3,u2c) = ur*sinTheta + vr*cosTheta;
                        u(i1,i2,i3,pc ) = pr;
                        u1Max=max(u1Max,u(i1,i2,i3,u1c));
                        u2Max=max(u2Max,u(i1,i2,i3,u2c));
                        pMax =max(pMax, u(i1,i2,i3,pc));
                    }
                    uNorm = max(u1Max,u2Max,u3Max,pMax);
          // if( scale==0. )
          // {
          //   printF("UDKS: ERROR: annulus solution is zero! scale=%e\n",scale);
          //   OV_ABORT("ERROR");
          // }
          // scale = 1./scale;
          // c4 *= scale;  // set this so the velocity and acceleration are scaled too
          // FOR_3(i1,i2,i3,I1,I2,I3)
          // {
          //   u(i1,i2,i3,u1c) *= scale;
          //   u(i1,i2,i3,u2c) *= scale;
          //   u(i1,i2,i3,pc ) *= scale;
          // }    
                }
                else
                {
                    OV_ABORT("UDKS:ERROR: annulus solution : 3D called");
                }

                if( computeVelocity )
                {
                    printF("\n $$$$$$$$$$ userDefinedKnownSolution: EVALUATE VELOCITY, ACCEL and PRESSURE AT PAST TIMES t=%9.3e dt=%16.6e $$$$$$$\n\n",t,dt);
                    assert( dt>0. );
                    assert( vgf!=NULL );
                    const int currentVelocity = dbase.get<int>("currentVelocity");
                    assert( currentVelocity>=0 && currentVelocity<numberOfVelocityFunctions );
          // fill in t0=0 and some past time levels
                    for( int level=0; level<numberOfVelocityFunctions; level++ )
                    {
                        real tc = 0. - level*dt;
                        const int prevVelocity = (currentVelocity - level + numberOfVelocityFunctions) % numberOfVelocityFunctions;
                        OV_GET_SERIAL_ARRAY(real,vgf[prevVelocity][grid],vLocal);
                        OV_GET_SERIAL_ARRAY(real,vtgf[prevVelocity][grid],vtLocal);
                        FOR_3D(i1,i2,i3,I1,I2,I3)
                        {
                            Real x = vertex(i1,i2,i3,0);
                            Real y = vertex(i1,i2,i3,1);
                            Real r = sqrt( x*x + y*y );
                            Real theta = atan2(y,x); 
                            Real cosTheta = x/r;
                            Real sinTheta = y/r;
                            t1 = omega * tc;
                            t2 = cos(t1);
                            t3 = n * theta;
                            t4 = sin(t3);
                            t6 = sin(t1);
                            t7 = cos(t3);
                            t11 = omega * omega;
                            t12 = pow(r, -n);
                            t13 = yn(n, omega);
                            t14 = t13 * t12;
                            t15 = omega * r;
                            t16 = yn(n, t15);
                            t19 = n + 0.1e1;
                            t21 = pow(0.2e1, 0.2e1 * t19);
                            t24 = pow(0.8e1, n);
                            t25 = n * t24;
                            t26 = n * n;
                            t27 = t11 / 0.8e1;
                            t28 = t26 - t27 - 0.1e1;
                            t29 = pow(r, n);
                            t31 = (t12 - t29) * t28;
                            t32 = omega / 0.2e1;
                            t33 = yn(n, t32);
                            t36 = t26 * n;
                            t37 = -t36 + n;
                            t39 = pow(0.4e1, n);
                            t40 = t12 * t39;
                            t42 = -t36 - t27 + n;
                            t43 = pow(0.16e2, n);
                            t47 = -t39 * t37 + t43 * t42;
                            t48 = t16 * t47;
                            t49 = t29 * t43;
                            t50 = -t42;
                            t52 = t50 * t13 * t49;
                            t55 = jn(t19, t32);
                            t58 = jn(n, omega);
                            t59 = t58 * t12;
                            t60 = jn(n, t15);
                            t65 = jn(n, t32);
                            t70 = t60 * t47;
                            t72 = t50 * t58 * t49;
                            t75 = yn(t19, t32);
                            t84 = t65 * (-t13 * t29 + t16) + (t58 * t29 - t60) * t33;
                            t87 = pow(0.2e1, 0.4e1 * t19);
                            t88 = t87 * t11;
                            t92 = t26 - t11 / 0.16e2 - n;
                            t93 = t92 * t39;
                            t97 = t12 * t13 * t11 * t93 / 0.8e1;
                            t98 = t26 * t26;
                            t99 = t11 * t11;
                            t101 = t98 - t26 + t99 / 0.128e3;
                            t104 = t11 * t93 / 0.8e1;
                            t105 = t43 * t101 - t104;
                            t106 = t16 * t105;
                            t109 = t43 * t13 * t101 * t29;
                            t112 = t39 * t58;
                            t116 = t12 * t11 * t92 * t112 / 0.8e1;
                            t120 = t60 * (-t43 * t101 + t104);
                            t122 = t101 * t58 * t49;
                            t129 = t28 * t26;
                            t130 = t129 * (-t60 * t13 + t16 * t58) * t24;
                            t133 = 0.1e1 / r;
                            t137 = t43 * t50;
                            t138 = t39 * t37;
                            t147 = t65 * t105;
                            t149 = t129 * t58 * t24;
                            t151 = 0.1e1 / (-t55 * omega * (t21 * t11 / 0.32e2 + t137 + t138) / 0.4e1 - t11 * t26 * t65 * t87 / 0.128e3 + t147 - t149);
                            ur = -t151 * t133 * omega * (t55 * omega * (t21 * (t14 - t16) * t11 / 0.32e2 + t33 * t31 * t25 + t40 * t13 * t37 + t48 + t52) / 0.4e1 - t75 * (t21 * (t59 - t60) * t11 / 0.32e2 + t65 * t31 * t25 + t40 * t58 * t37 + t70 + t72) * omega / 0.4e1 - t88 * t26 * t84 / 0.128e3 + t65 * (t97 + t106 - t109) + t33 * (-t116 + t120 + t122) - t130) * (-t4 * t2 + t7 * t6) * c4;
                            t154 = yn(t19, t15);
                            t155 = omega * t154;
                            t164 = r * (t137 + t138);
                            t167 = (t29 + t12) * t24;
                            t168 = t28 * n;
                            t175 = t12 * t19 * (n - 0.1e1);
                            t183 = jn(t19, t15);
                            t184 = omega * t183;
                            vr = 0.1e1 / n * t133 * t151 * (t7 * t2 + t4 * t6) * omega * (t55 * (-t21 * t11 * (-r * t155 + (t14 + t16) * n) / 0.32e2 + t155 * t164 + n * (t175 * n * t13 * t39 - t33 * t168 * t167 + t48 + t52)) * omega / 0.4e1 - t75 * omega * (-t21 * t11 * (-r * t184 + (t59 + t60) * n) / 0.32e2 + t184 * t164 + (t175 * n * t112 - t65 * t168 * t167 + t70 + t72) * n) / 0.4e1 - t88 * t26 * (-t15 * t154 * t65 + t15 * t183 * t33 + n * t84) / 0.128e3 + t184 * r * (-t13 * t28 * t26 * t24 + t33 * t105) - t155 * (t147 - t149) * r + (t65 * (-t97 + t106 - t109) + t33 * (t116 + t120 + t122) - t130) * n) * c4;
              // (u1,u2) = ur*rHat + vr*thetaHat 
                            vLocal(i1,i2,i3,u1c) = ur*cosTheta - vr*sinTheta;
                            vLocal(i1,i2,i3,u2c) = ur*sinTheta + vr*cosTheta;
                            t1 = omega * omega;
                            t2 = pow(r, -n);
                            t3 = yn(n, omega);
                            t4 = t3 * t2;
                            t5 = omega * r;
                            t6 = yn(n, t5);
                            t9 = n + 0.1e1;
                            t11 = pow(0.2e1, 0.2e1 * t9);
                            t14 = pow(0.8e1, n);
                            t15 = n * t14;
                            t16 = n * n;
                            t17 = t1 / 0.8e1;
                            t18 = t16 - t17 - 0.1e1;
                            t19 = pow(r, n);
                            t21 = (t2 - t19) * t18;
                            t22 = omega / 0.2e1;
                            t23 = yn(n, t22);
                            t26 = t16 * n;
                            t27 = -t26 + n;
                            t29 = pow(0.4e1, n);
                            t30 = t2 * t29;
                            t32 = -t26 - t17 + n;
                            t33 = pow(0.16e2, n);
                            t37 = -t29 * t27 + t33 * t32;
                            t38 = t6 * t37;
                            t39 = t19 * t33;
                            t40 = -t32;
                            t42 = t40 * t3 * t39;
                            t45 = jn(t9, t22);
                            t48 = jn(n, omega);
                            t49 = t48 * t2;
                            t50 = jn(n, t5);
                            t55 = jn(n, t22);
                            t60 = t50 * t37;
                            t62 = t40 * t48 * t39;
                            t65 = yn(t9, t22);
                            t74 = t55 * (-t3 * t19 + t6) + (t48 * t19 - t50) * t23;
                            t77 = pow(0.2e1, 0.4e1 * t9);
                            t78 = t77 * t1;
                            t82 = t16 - t1 / 0.16e2 - n;
                            t83 = t82 * t29;
                            t87 = t2 * t3 * t1 * t83 / 0.8e1;
                            t88 = t16 * t16;
                            t89 = t1 * t1;
                            t91 = t88 - t16 + t89 / 0.128e3;
                            t94 = t1 * t83 / 0.8e1;
                            t95 = t33 * t91 - t94;
                            t96 = t6 * t95;
                            t99 = t33 * t3 * t91 * t19;
                            t102 = t29 * t48;
                            t106 = t2 * t1 * t82 * t102 / 0.8e1;
                            t110 = t50 * (-t33 * t91 + t94);
                            t112 = t91 * t48 * t39;
                            t119 = t18 * t16;
                            t120 = t119 * (-t50 * t3 + t6 * t48) * t14;
                            t124 = omega * tc;
                            t125 = cos(t124);
                            t126 = n * theta;
                            t127 = cos(t126);
                            t129 = sin(t124);
                            t130 = sin(t126);
                            t133 = 0.1e1 / r;
                            t137 = t33 * t40;
                            t138 = t29 * t27;
                            t147 = t55 * t95;
                            t149 = t119 * t48 * t14;
                            t151 = 0.1e1 / (-t45 * omega * (t1 * t11 / 0.32e2 + t137 + t138) / 0.4e1 - t1 * t16 * t55 * t77 / 0.128e3 + t147 - t149);
                            ur = -t151 * t133 * (t127 * t125 + t130 * t129) * t1 * (t45 * omega * (t11 * (t4 - t6) * t1 / 0.32e2 + t23 * t21 * t15 + t30 * t3 * t27 + t38 + t42) / 0.4e1 - t65 * (t11 * (t49 - t50) * t1 / 0.32e2 + t55 * t21 * t15 + t30 * t48 * t27 + t60 + t62) * omega / 0.4e1 - t78 * t16 * t74 / 0.128e3 + t55 * (t87 + t96 - t99) + t23 * (-t106 + t110 + t112) - t120) * c4;
                            t158 = yn(t9, t5);
                            t159 = omega * t158;
                            t168 = r * (t137 + t138);
                            t171 = (t19 + t2) * t14;
                            t172 = t18 * n;
                            t179 = t2 * t9 * (n - 0.1e1);
                            t187 = jn(t9, t5);
                            t188 = omega * t187;
                            vr = -0.1e1 / n * t151 * t133 * t1 * (t45 * (-t11 * t1 * (-r * t159 + (t4 + t6) * n) / 0.32e2 + t159 * t168 + n * (t179 * n * t3 * t29 - t23 * t172 * t171 + t38 + t42)) * omega / 0.4e1 - t65 * omega * (-t11 * t1 * (-r * t188 + (t49 + t50) * n) / 0.32e2 + t188 * t168 + (t179 * n * t102 - t55 * t172 * t171 + t60 + t62) * n) / 0.4e1 - t78 * t16 * (-t5 * t158 * t55 + t5 * t187 * t23 + n * t74) / 0.128e3 + t188 * r * (-t3 * t18 * t16 * t14 + t23 * t95) - t159 * (t147 - t149) * r + (t55 * (-t87 + t96 - t99) + t23 * (t106 + t110 + t112) - t120) * n) * (-t130 * t125 + t127 * t129) * c4;
              // (u1,u2) = ur*rHat + vr*thetaHat 
                            vtLocal(i1,i2,i3,u1c) = ur*cosTheta - vr*sinTheta;
                            vtLocal(i1,i2,i3,u2c) = ur*sinTheta + vr*cosTheta;        
                        }      
            // FOR_3D(i1,i2,i3,I1,I2,I3)
            // {
            //   Real x = vertex(i1,i2,i3,0);
            //   Real y = vertex(i1,i2,i3,1);
            //   vLocal(i1,i2,i3,u1c)  = v1e(x,y,tc);
            //   vLocal(i1,i2,i3,u2c)  = v2e(x,y,tc);
            //   vtLocal(i1,i2,i3,u1c) = v1et(x,y,tc);
            //   vtLocal(i1,i2,i3,u2c) = v2et(x,y,tc);        
            // }         
                    }
          // --- save pressure for time extrapolation ---
                    const int currentPressure = dbase.get<int>("currentPressure");
                    assert( currentPressure>=0 && currentPressure<max(1,numberOfPressureFunctions) );
                    for( int level=0; level<numberOfPressureFunctions; level++ )
                    {
                        real tc = 0. - (level+3)*dt;  // -- save pressure starting at t-3*dt 
                        const int prevPressure = (currentPressure - level + numberOfPressureFunctions) % numberOfPressureFunctions;
                        assert( pgf !=NULL );
                        OV_GET_SERIAL_ARRAY(real,pgf[prevPressure][grid],pLocal);
                        FOR_3D(i1,i2,i3,I1,I2,I3)
                        {
                            Real x = vertex(i1,i2,i3,0);
                            Real y = vertex(i1,i2,i3,1);
                            Real r = sqrt( x*x + y*y );
                            Real theta = atan2(y,x); 
                            Real cosTheta = x/r;
                            Real sinTheta = y/r;
                            t1 = n * n;
                            t2 = omega * omega;
                            t3 = t2 * t1;
                            t4 = omega / 0.2e1;
                            t5 = jn(n, t4);
                            t6 = pow(0.1e1 / 0.2e1, n);
                            t7 = t6 * t5;
                            t10 = n + 0.1e1;
                            t11 = jn(t10, t4);
                            t12 = t6 * t11;
                            t13 = t1 * n;
                            t14 = omega * t13;
                            t17 = omega * t11;
                            t18 = t6 * n;
                            t21 = 0.1e1 / t6;
                            t25 = t2 * t2;
                            t27 = t2 * omega;
                            t30 = t5 * t21;
                            t31 = t1 * t1;
                            t35 = t21 * t11;
                            t44 = jn(n, omega);
                            t47 = t2 * n;
                            t50 = t1 * t44;
                            t54 = -0.32e2 * t21 * n * t17 + 0.128e3 * t1 * t30 - 0.32e2 * t14 * t12 + 0.4e1 * t27 * t12 + 0.32e2 * t14 * t35 + 0.32e2 * t18 * t17 - 0.16e2 * t2 * t50 - t25 * t30 - t25 * t7 + 0.4e1 * t27 * t35 + 0.16e2 * t3 * t30 + 0.16e2 * t3 * t7 - 0.128e3 * t31 * t30 + 0.128e3 * t31 * t44 - 0.16e2 * t47 * t7 - 0.128e3 * t50;
                            t55 = 0.1e1 / t54;
                            t56 = t55 * c4;
                            t57 = pow(r, n);
                            t58 = t57 * t56;
                            t59 = t58 * t3;
                            t60 = 0.1e1 / r;
                            t61 = t21 * t60;
                            t62 = t5 * t61;
                            t63 = yn(n, omega);
                            t64 = omega * tc;
                            t65 = sin(t64);
                            t66 = t65 * t63;
                            t67 = n * theta;
                            t68 = sin(t67);
                            t69 = t68 * t66;
                            t73 = yn(n, t4);
                            t74 = t73 * t61;
                            t75 = t65 * t44;
                            t76 = t68 * t75;
                            t80 = cos(t64);
                            t81 = t80 * t63;
                            t82 = cos(t67);
                            t83 = t82 * t81;
                            t87 = t80 * t44;
                            t88 = t82 * t87;
                            t92 = t58 * t14;
                            t93 = yn(t10, t4);
                            t94 = t93 * t61;
                            t95 = t88 * t94;
                            t98 = n * omega;
                            t99 = t58 * t98;
                            t100 = t11 * t60;
                            t101 = t63 * t100;
                            t102 = t80 * t21;
                            t103 = t82 * t102;
                            t104 = t103 * t101;
                            t109 = 0.1e1 / t57;
                            t110 = t109 * t56;
                            t111 = t110 * t3;
                            t112 = t5 * t60;
                            t113 = t63 * t112;
                            t114 = t65 * t6;
                            t115 = t68 * t114;
                            t116 = t115 * t113;
                            t119 = t110 * t14;
                            t120 = t115 * t101;
                            t124 = t44 * t73 * t60;
                            t125 = t115 * t124;
                            t129 = t44 * t93 * t60;
                            t130 = t115 * t129;
                            t133 = t11 * t61;
                            t134 = t69 * t133;
                            t137 = t76 * t94;
                            t140 = -0.16e2 * t69 * t62 * t59 - 0.16e2 * t83 * t62 * t59 + 0.16e2 * t76 * t74 * t59 + 0.16e2 * t88 * t74 * t59 + 0.32e2 * t104 * t99 - 0.16e2 * t116 * t111 + 0.16e2 * t125 * t111 + 0.32e2 * t120 * t119 - 0.32e2 * t130 * t119 - 0.32e2 * t134 * t92 + 0.32e2 * t137 * t92 + 0.32e2 * t95 * t92 - 0.32e2 * t95 * t99;
                            t145 = t80 * t6;
                            t146 = t82 * t145;
                            t147 = t146 * t113;
                            t150 = t146 * t101;
                            t153 = t146 * t124;
                            t156 = t146 * t129;
                            t159 = t110 * t47;
                            t164 = t110 * t98;
                            t171 = t65 * t93;
                            t172 = t68 * t171;
                            t173 = t172 * t112;
                            t176 = t65 * t73;
                            t177 = t68 * t176;
                            t178 = t177 * t100;
                            t181 = t55 * t60;
                            t182 = omega * r;
                            t183 = jn(n, t182);
                            t184 = t183 * c4;
                            t186 = t73 * t184 * t181;
                            t187 = t1 * t6;
                            t188 = t65 * t2;
                            t189 = t68 * t188;
                            t193 = -0.16e2 * t189 * t187 * t186 - 0.32e2 * t104 * t92 - 0.16e2 * t147 * t111 + 0.16e2 * t153 * t111 + 0.32e2 * t150 * t119 - 0.32e2 * t156 * t119 + 0.32e2 * t173 * t119 - 0.32e2 * t178 * t119 + 0.32e2 * t134 * t99 - 0.32e2 * t137 * t99 + 0.16e2 * t147 * t159 - 0.32e2 * t150 * t164 - 0.16e2 * t153 * t159 + 0.32e2 * t156 * t164;
                            t196 = t93 * t184 * t181;
                            t197 = t13 * t6;
                            t198 = t65 * omega;
                            t199 = t68 * t198;
                            t206 = t25 * c4;
                            t207 = t109 * t55;
                            t208 = t60 * t207;
                            t209 = t208 * t206;
                            t210 = t44 * t73;
                            t211 = t115 * t210;
                            t213 = t27 * c4;
                            t214 = t208 * t213;
                            t215 = t63 * t11;
                            t216 = t115 * t215;
                            t219 = t44 * t93;
                            t220 = t115 * t219;
                            t223 = t57 * t55;
                            t224 = t60 * t223;
                            t225 = t224 * t213;
                            t226 = t65 * t21;
                            t227 = t68 * t226;
                            t228 = t227 * t215;
                            t231 = t227 * t219;
                            t235 = t224 * c4 * t1;
                            t236 = t69 * t30;
                            t239 = t184 * t181;
                            t240 = t31 * t63;
                            t241 = t68 * t65;
                            t245 = t1 * t63;
                            t249 = t82 * t80;
                            t256 = c4 * t60;
                            t257 = yn(n, t182);
                            t258 = t65 * t257;
                            t261 = t80 * t257;
                            t264 = 0.16e2 * t189 * t18 * t186 + 0.32e2 * t199 * t197 * t196 - 0.128e3 * t241 * t240 * t239 - 0.128e3 * t249 * t240 * t239 + 0.128e3 * t241 * t245 * t239 + 0.128e3 * t249 * t245 * t239 + t68 * t258 * t256 + t82 * t261 * t256 - t211 * t209 - 0.4e1 * t216 * t214 + 0.4e1 * t220 * t214 - 0.4e1 * t228 * t225 + 0.4e1 * t231 * t225 - 0.128e3 * t236 * t235;
                            t265 = t83 * t30;
                            t268 = t73 * t21;
                            t269 = t88 * t268;
                            t272 = t224 * t206;
                            t276 = t224 * c4 * t31;
                            t282 = t21 * t184 * t181;
                            t283 = t13 * t93;
                            t284 = t80 * omega;
                            t285 = t82 * t284;
                            t295 = t80 * t73;
                            t296 = t82 * t295;
                            t297 = t296 * t100;
                            t300 = t80 * t93;
                            t301 = t82 * t300;
                            t302 = t301 * t112;
                            t305 = t76 * t268;
                            t311 = 0.32e2 * t103 * t98 * t196 - 0.32e2 * t146 * t98 * t196 - 0.32e2 * t285 * t283 * t282 - 0.128e3 * t265 * t235 + 0.128e3 * t269 * t235 + 0.128e3 * t305 * t235 + 0.128e3 * t236 * t276 + t265 * t272 + 0.128e3 * t265 * t276 - t269 * t272 - 0.128e3 * t269 * t276 - t305 * t272 - 0.32e2 * t297 * t99 + 0.32e2 * t302 * t99;
                            t322 = t27 * n;
                            t323 = t110 * t322;
                            t326 = t1 * t73;
                            t342 = t80 * t2;
                            t343 = t82 * t342;
                            t349 = -0.32e2 * t115 * t98 * t196 - 0.16e2 * t189 * t326 * t282 - 0.32e2 * t199 * t283 * t282 - 0.16e2 * t343 * t326 * t282 + 0.16e2 * t116 * t159 - 0.32e2 * t120 * t164 - 0.16e2 * t125 * t159 + 0.32e2 * t130 * t164 + 0.4e1 * t178 * t323 - 0.128e3 * t305 * t276 + 0.4e1 * t297 * t323 - 0.4e1 * t302 * t323 - 0.32e2 * t302 * t92;
                            t361 = t63 * t5;
                            t362 = t115 * t361;
                            t364 = t146 * t361;
                            t366 = t146 * t210;
                            t368 = t146 * t215;
                            t371 = t146 * t219;
                            t374 = t103 * t215;
                            t377 = t103 * t219;
                            t380 = t58 * t322;
                            t383 = -0.32e2 * t297 * t119 + 0.32e2 * t302 * t119 + 0.32e2 * t297 * t164 - 0.32e2 * t302 * t164 + t362 * t209 + t364 * t209 - t366 * t209 - 0.4e1 * t368 * t214 + 0.4e1 * t371 * t214 - 0.4e1 * t374 * t225 + 0.4e1 * t377 * t225 + t236 * t272 + 0.32e2 * t297 * t92 + 0.4e1 * t302 * t380;
                            t417 = 0.16e2 * t343 * t18 * t186 - 0.16e2 * t343 * t187 * t186 + 0.32e2 * t285 * t197 * t196 + 0.32e2 * t227 * t98 * t196 - 0.32e2 * t173 * t164 + 0.32e2 * t178 * t164 - 0.4e1 * t173 * t323 + 0.4e1 * t173 * t380 - 0.32e2 * t173 * t92 + 0.32e2 * t173 * t99 - 0.4e1 * t178 * t380 + 0.32e2 * t178 * t92 - 0.32e2 * t178 * t99 - 0.4e1 * t297 * t380;
                            t418 = t21 * t93;
                            t419 = t65 * t27;
                            t420 = t68 * t419;
                            t424 = t65 * t1;
                            t429 = t6 * t93;
                            t433 = t6 * t73;
                            t434 = t65 * t25;
                            t435 = t68 * t434;
                            t440 = t80 * t25;
                            t441 = t82 * t440;
                            t444 = t80 * t27;
                            t445 = t82 * t444;
                            t455 = t80 * t31;
                            t465 = t80 * t1;
                            t470 = t65 * t31;
                            t475 = -0.128e3 * t68 * t424 * t268 * t239 + 0.128e3 * t82 * t455 * t268 * t239 - 0.128e3 * t82 * t465 * t268 * t239 + 0.128e3 * t68 * t470 * t268 * t239 + 0.16e2 * t189 * t245 * t239 + 0.16e2 * t343 * t245 * t239 + t435 * t268 * t239 + t441 * t268 * t239 - 0.4e1 * t420 * t418 * t239 - 0.4e1 * t445 * t418 * t239 - 0.4e1 * t420 * t429 * t239 - 0.4e1 * t445 * t429 * t239 + t435 * t433 * t239 + t441 * t433 * t239;
                            ur = t140 + t193 + t264 + t311 + t349 + t383 + t417 + t475;
                            t478 = t2 * t60;
                            t481 = c4 * t57 * t1 * t478;
                            t482 = t21 * t55;
                            t483 = t5 * t482;
                            t484 = t82 * t66;
                            t488 = t68 * t81;
                            t492 = t73 * t482;
                            t493 = t68 * t87;
                            t497 = t82 * t114;
                            t498 = t497 * t113;
                            t501 = t497 * t101;
                            t504 = t497 * t124;
                            t507 = t497 * t129;
                            t518 = t484 * t133;
                            t521 = t82 * t75;
                            t522 = t521 * t94;
                            t525 = 0.1e1 / n;
                            t527 = yn(t10, t182);
                            t528 = t527 * t525 * omega;
                            t535 = t68 * t284;
                            t539 = t68 * t295;
                            t540 = t539 * t100;
                            t543 = t68 * t300;
                            t544 = t543 * t112;
                            t547 = -t82 * t65 * c4 * t528 + t68 * t80 * c4 * t528 - 0.32e2 * t535 * t197 * t196 - 0.16e2 * t484 * t483 * t481 + 0.16e2 * t488 * t483 * t481 - 0.16e2 * t493 * t492 * t481 + 0.16e2 * t498 * t111 - 0.16e2 * t504 * t111 - 0.32e2 * t501 * t119 + 0.32e2 * t507 * t119 - 0.16e2 * t498 * t159 + 0.16e2 * t504 * t159 + 0.32e2 * t501 * t164 - 0.32e2 * t507 * t164 + 0.32e2 * t540 * t164 - 0.32e2 * t544 * t164 - 0.32e2 * t518 * t92 + 0.32e2 * t522 * t92;
                            t552 = t68 * t342;
                            t556 = t68 * t145;
                            t567 = t68 * t102;
                            t574 = t488 * t30;
                            t577 = t493 * t268;
                            t584 = t27 * t60;
                            t587 = c4 * t57 * n * t584;
                            t588 = t11 * t55;
                            t589 = t539 * t588;
                            t592 = jn(t10, t182);
                            t593 = t55 * t592;
                            t594 = t593 * t322;
                            t595 = t63 * c4;
                            t596 = t68 * t80;
                            t597 = t596 * t595;
                            t600 = t82 * t65;
                            t601 = t600 * t595;
                            t605 = t56 * t592 * t2;
                            t609 = -0.16e2 * t552 * t18 * t186 - t556 * t210 * t209 + t556 * t361 * t209 - 0.4e1 * t556 * t215 * t214 + 0.4e1 * t556 * t219 * t214 + 0.4e1 * t567 * t215 * t225 - 0.4e1 * t567 * t219 * t225 - 0.32e2 * t596 * t429 * t605 - 0.32e2 * t540 * t119 + 0.32e2 * t544 * t119 + 0.128e3 * t574 * t235 - 0.128e3 * t577 * t235 - t574 * t272 + t577 * t272 - 0.128e3 * t574 * t276 + 0.4e1 * t589 * t587 + 0.16e2 * t597 * t594 - 0.16e2 * t601 * t594;
                            t612 = t56 * t592 * t27;
                            t631 = t593 * t14;
                            t634 = t593 * t98;
                            t654 = t488 * t133;
                            t657 = t493 * t94;
                            t662 = 0.128e3 * t596 * t240 * t239 - 0.128e3 * t600 * t240 * t239 - 0.128e3 * t596 * t245 * t239 + 0.128e3 * t600 * t245 * t239 + 0.32e2 * t596 * t418 * t605 - 0.32e2 * t600 * t418 * t605 + 0.32e2 * t600 * t429 * t605 + 0.16e2 * t596 * t433 * t612 - 0.16e2 * t600 * t433 * t612 + 0.32e2 * t518 * t99 - 0.32e2 * t522 * t99 - 0.128e3 * t597 * t631 + 0.128e3 * t597 * t634 + 0.128e3 * t601 * t631 - 0.128e3 * t601 * t634 + 0.32e2 * t654 * t92 - 0.32e2 * t654 * t99 - 0.32e2 * t657 * t92;
                            t665 = t556 * t113;
                            t672 = t556 * t101;
                            t675 = t556 * t124;
                            t678 = t556 * t129;
                            t692 = t183 * t1;
                            t694 = t55 * t692 * t478;
                            t695 = t21 * c4;
                            t696 = t539 * t695;
                            t701 = c4 * t109 * n * t584;
                            t702 = t5 * t55;
                            t703 = t543 * t702;
                            t715 = -t497 * t361 * t209 + t82 * t258 * t256 - t68 * t261 * t256 + 0.32e2 * t535 * t283 * t282 + 0.16e2 * t521 * t492 * t481 - 0.16e2 * t665 * t111 + 0.16e2 * t675 * t111 + 0.32e2 * t672 * t119 - 0.32e2 * t678 * t119 + 0.16e2 * t665 * t159 - 0.16e2 * t675 * t159 - 0.32e2 * t672 * t164 + 0.32e2 * t678 * t164 + 0.128e3 * t577 * t276 + 0.4e1 * t589 * t701 + 0.32e2 * t657 * t99 + 0.16e2 * t696 * t694 - 0.4e1 * t703 * t701;
                            t732 = t82 * t176;
                            t733 = t732 * t588;
                            t736 = t82 * t171;
                            t737 = t736 * t702;
                            t742 = t732 * t695;
                            t749 = t736 * t112;
                            t752 = t732 * t100;
                            t764 = 0.16e2 * t552 * t187 * t186 + 0.32e2 * t556 * t98 * t196 - 0.32e2 * t567 * t98 * t196 - 0.32e2 * t749 * t119 + 0.32e2 * t749 * t164 - 0.32e2 * t752 * t164 - 0.32e2 * t540 * t92 + 0.32e2 * t540 * t99 + 0.32e2 * t544 * t92 - 0.32e2 * t544 * t99 - 0.4e1 * t703 * t587 - 0.4e1 * t733 * t587 + 0.4e1 * t737 * t587 - 0.16e2 * t742 * t694 - 0.4e1 * t733 * t701 + 0.4e1 * t737 * t701 - 0.32e2 * t749 * t92 + 0.32e2 * t752 * t92;
                            t767 = t82 * t188;
                            t771 = t82 * t198;
                            t786 = t82 * t226;
                            t793 = t484 * t30;
                            t796 = t521 * t268;
                            t814 = 0.16e2 * t767 * t18 * t186 - 0.16e2 * t767 * t187 * t186 + 0.32e2 * t771 * t197 * t196 - 0.32e2 * t497 * t98 * t196 + 0.32e2 * t786 * t98 * t196 + t497 * t210 * t209 + 0.4e1 * t497 * t215 * t214 - 0.4e1 * t497 * t219 * t214 - 0.4e1 * t786 * t215 * t225 + 0.4e1 * t786 * t219 * t225 - 0.32e2 * t771 * t283 * t282 + 0.32e2 * t752 * t119 - 0.128e3 * t793 * t235 + 0.128e3 * t796 * t235 + t793 * t272 - t796 * t272 + 0.128e3 * t793 * t276 - 0.128e3 * t796 * t276;
                            t820 = t692 * t478;
                            t828 = t593 * t3;
                            t829 = t736 * t695;
                            t832 = t25 * omega;
                            t834 = t593 * t525 * t832;
                            t835 = t73 * c4;
                            t836 = t497 * t835;
                            t839 = t593 * t525 * t25;
                            t840 = t93 * c4;
                            t841 = t497 * t840;
                            t848 = t556 * t840;
                            t851 = t567 * t840;
                            t864 = -0.16e2 * t488 * t56 * t820 - 0.16e2 * t696 * t594 + 0.16e2 * t742 * t594 + 0.128e3 * t696 * t631 - 0.128e3 * t742 * t631 - 0.128e3 * t696 * t634 + 0.128e3 * t742 * t634 + t696 * t834 - t742 * t834 + 0.32e2 * t749 * t99 - 0.32e2 * t752 * t99 + 0.32e2 * t829 * t828 + 0.32e2 * t848 * t828 + 0.4e1 * t829 * t839 - t836 * t834 + 0.4e1 * t841 * t839 - 0.4e1 * t848 * t839 - 0.4e1 * t851 * t839;
                            t865 = t556 * t835;
                            t878 = t68 * t440;
                            t881 = t68 * t444;
                            t902 = t82 * t419;
                            t910 = t82 * t434;
                            t918 = -0.128e3 * t82 * t424 * t268 * t239 - 0.128e3 * t68 * t455 * t268 * t239 + 0.128e3 * t68 * t465 * t268 * t239 + 0.128e3 * t82 * t470 * t268 * t239 - t878 * t268 * t239 + t910 * t268 * t239 + 0.4e1 * t881 * t418 * t239 - 0.4e1 * t902 * t418 * t239 + 0.4e1 * t881 * t429 * t239 - 0.4e1 * t902 * t429 * t239 - t878 * t433 * t239 + t910 * t433 * t239 + 0.16e2 * t484 * t56 * t820 + 0.16e2 * t836 * t594 - 0.16e2 * t865 * t594 - 0.32e2 * t841 * t828 - 0.32e2 * t851 * t828 + t865 * t834;
                            vr = t547 + t609 + t662 + t715 + t764 + t814 + t864 + t918;
                            t921 = t832 * c4;
                            t922 = t525 * t55;
                            t923 = t57 * t922;
                            t924 = t923 * t921;
                            t927 = t1 * t55;
                            t929 = t57 * t927 * t213;
                            t932 = n * t55;
                            t933 = t57 * t932;
                            t934 = t933 * t206;
                            t938 = t109 * t927 * t213;
                            t948 = t109 * t932 * t206;
                            t958 = t25 * t2 * c4;
                            t959 = t109 * t922;
                            t960 = t959 * t958;
                            t964 = t2 * c4;
                            t965 = t933 * t964;
                            t968 = t923 * t958;
                            t971 = t959 * t921;
                            t976 = t223 * t921;
                            t977 = t73 * t11;
                            t978 = t249 * t977;
                            t981 = -0.32e2 * t228 * t929 + 0.4e1 * t231 * t924 + t236 * t968 - 0.16e2 * t265 * t934 - 0.128e3 * t265 * t965 + 0.16e2 * t269 * t934 + 0.16e2 * t305 * t934 + 0.16e2 * t364 * t948 - t364 * t960 - 0.16e2 * t366 * t948 + t366 * t960 - 0.32e2 * t368 * t938 + 0.4e1 * t368 * t971 + 0.32e2 * t371 * t938 - 0.4e1 * t371 * t971 - 0.32e2 * t374 * t929 + 0.4e1 * t377 * t924 + 0.32e2 * t377 * t929 - 0.4e1 * t978 * t976;
                            t982 = t207 * t213;
                            t985 = t93 * t5;
                            t986 = t249 * t985;
                            t989 = t207 * t921;
                            t994 = t223 * t213;
                            t999 = t241 * t985;
                            t1002 = t241 * t977;
                            t1025 = t927 * t213;
                            t1026 = t5 * t57;
                            t1030 = -0.32e2 * t301 * t1026 * t1025 - 0.4e1 * t1002 * t976 - 0.32e2 * t1002 * t982 - 0.4e1 * t1002 * t989 - 0.32e2 * t1002 * t994 - 0.32e2 * t371 * t982 + 0.32e2 * t374 * t994 - 0.32e2 * t377 * t994 + 0.4e1 * t986 * t976 + 0.4e1 * t999 * t976 - 0.32e2 * t978 * t982 - 0.4e1 * t978 * t989 - 0.32e2 * t978 * t994 + 0.32e2 * t986 * t982 + 0.32e2 * t999 * t982 + 0.4e1 * t986 * t989 + 0.32e2 * t986 * t994 + 0.4e1 * t999 * t989 + 0.32e2 * t999 * t994;
                            t1032 = t11 * t57;
                            t1036 = t207 * t206;
                            t1045 = t5 * t109;
                            t1049 = t11 * t109;
                            t1077 = -0.32e2 * t172 * t1026 * t1025 + 0.32e2 * t296 * t1032 * t1025 - 0.32e2 * t301 * t1045 * t1025 + 0.32e2 * t296 * t1049 * t1025 + 0.16e2 * t211 * t1036 - 0.16e2 * t362 * t1036 + 0.4e1 * t216 * t971 + 0.32e2 * t216 * t982 + 0.32e2 * t220 * t938 - 0.4e1 * t220 * t971 - 0.32e2 * t220 * t982 + 0.32e2 * t228 * t994 + 0.32e2 * t231 * t929 - 0.32e2 * t231 * t994 - 0.16e2 * t236 * t934 - 0.128e3 * t236 * t965 + 0.128e3 * t305 * t965 + 0.16e2 * t362 * t948 - t362 * t960;
                            t1081 = t57 * t13 * t55 * t964;
                            t1118 = 0.32e2 * t177 * t1032 * t1025 - 0.32e2 * t172 * t1045 * t1025 + 0.32e2 * t177 * t1049 * t1025 - 0.16e2 * t364 * t1036 + 0.16e2 * t366 * t1036 + 0.128e3 * t236 * t1081 + 0.128e3 * t265 * t1081 - 0.128e3 * t269 * t1081 - 0.128e3 * t305 * t1081 - 0.16e2 * t211 * t948 + t211 * t960 - 0.32e2 * t216 * t938 - 0.4e1 * t228 * t924 + t265 * t968 + 0.128e3 * t269 * t965 - t269 * t968 - t305 * t968 + 0.32e2 * t368 * t982 - 0.4e1 * t374 * t924;
                            pr = t981 + t1030 + t1077 + t1118;
                            pLocal(i1,i2,i3) = pr;
                        }
                    }
                }

        }    
        else if( iswCase==3 )
        {
      // ---- Annulus : traction-traction ----
      // printF("UDKS: EVAL NEW ANNULUS SOLUTION\n");
      // int n=1, m=2; 
            Real c4=1./4.; 

// --- roots of the frequency equation -----
//   File written by diskAnnulusSolution.mw on 2021-10-18
//   Solution is of the form Jn(omega*r) with roots omega_{n,m} 
Real omegaArray[]={
    // n=1, m=1,2,..,4
        3.4346674494392022e+00, 7.7407996399021592e+00, 1.3009389736124731e+01, 1.9138855322118575e+01,
    // n=2, m=1,2,..,4
        4.7964119338586172e+00, 9.4022451914702559e+00, 1.3622383066594675e+01, 1.9446727547714004e+01,
    // n=3, m=1,2,..,4
        2.5439592523077787e+00, 6.2948319530425706e+00, 1.0674374290844092e+01, 1.4841467907331412e+01
};
#define omegaRoot(n,m) omegaArray[(m-1)+4*(n-1)]
// ------------------------------------------------------------------------------------------------ 
//    caseName: Exact solution for incompressible linear elasticity 
// Annulus of outer radius Rb=2, inner radius Ra = 1; r=1 Traction BC; r=2 Traction BC. 
// File Written by diskAnnulusSolution.mw (maple) on 2021-10-18 
// Set c3 to scale disk solution and c4 to scale the annulus solution 
// Set n and m to choose the solution, Jn(omega*r), omega=omega(n,m) 
// The eigenvalues are save in a separate include file. 
// ------------------------------------------------------------------------------------------------ 
Real rho,c,mu,omega; 
rho=1.;
mu=1.;
omega=omegaRoot(n,m);
Real t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49
        ,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,t86,t87,t88,t89,t90,t91,t92,t93,t94,t95,t96,t97,t98,t99
        ,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t130,t131,t132,t133,t134,t135,t136,t137,t138,t139,t140,t141,t142,t143,t144,t145,t146,t147,t148,t149
        ,t150,t151,t152,t153,t154,t155,t156,t157,t158,t159,t160,t161,t162,t163,t164,t165,t166,t167,t168,t169,t170,t171,t172,t173,t174,t175,t176,t177,t178,t179,t180,t181,t182,t183,t184,t185,t186,t187,t188,t189,t190,t191,t192,t193,t194,t195,t196,t197,t198,t199
        ,t200,t201,t202,t203,t204,t205,t206,t207,t208,t209,t210,t211,t212,t213,t214,t215,t216,t217,t218,t219,t220,t221,t222,t223,t224,t225,t226,t227,t228,t229,t230,t231,t232,t233,t234,t235,t236,t237,t238,t239,t240,t241,t242,t243,t244,t245,t246,t247,t248,t249
        ,t250,t251,t252,t253,t254,t255,t256,t257,t258,t259,t260,t261,t262,t263,t264,t265,t266,t267,t268,t269,t270,t271,t272,t273,t274,t275,t276,t277,t278,t279,t280,t281,t282,t283,t284,t285,t286,t287,t288,t289,t290,t291,t292,t293,t294,t295,t296,t297,t298,t299
        ,t300,t301,t302,t303,t304,t305,t306,t307,t308,t309,t310,t311,t312,t313,t314,t315,t316,t317,t318,t319,t320,t321,t322,t323,t324,t325,t326,t327,t328,t329,t330,t331,t332,t333,t334,t335,t336,t337,t338,t339,t340,t341,t342,t343,t344,t345,t346,t347,t348,t349
        ,t350,t351,t352,t353,t354,t355,t356,t357,t358,t359,t360,t361,t362,t363,t364,t365,t366,t367,t368,t369,t370,t371,t372,t373,t374,t375,t376,t377,t378,t379,t380,t381,t382,t383,t384,t385,t386,t387,t388,t389,t390,t391,t392,t393,t394,t395,t396,t397,t398,t399
        ,t400,t401,t402,t403,t404,t405,t406,t407,t408,t409,t410,t411,t412,t413,t414,t415,t416,t417,t418,t419,t420,t421,t422,t423,t424,t425,t426,t427,t428,t429,t430,t431,t432,t433,t434,t435,t436,t437,t438,t439,t440,t441,t442,t443,t444,t445,t446,t447,t448,t449
        ,t450,t451,t452,t453,t454,t455,t456,t457,t458,t459,t460,t461,t462,t463,t464,t465,t466,t467,t468,t469,t470,t471,t472,t473,t474,t475,t476,t477,t478,t479,t480,t481,t482,t483,t484,t485,t486,t487,t488,t489,t490,t491,t492,t493,t494,t495,t496,t497,t498,t499
        ,t500,t501,t502,t503,t504,t505,t506,t507,t508,t509,t510,t511,t512,t513,t514,t515,t516,t517,t518,t519,t520,t521,t522,t523,t524,t525,t526,t527,t528,t529,t530,t531,t532,t533,t534,t535,t536,t537,t538,t539,t540,t541,t542,t543,t544,t545,t546,t547,t548,t549
        ,t550,t551,t552,t553,t554,t555,t556,t557,t558,t559,t560,t561,t562,t563,t564,t565,t566,t567,t568,t569,t570,t571,t572,t573,t574,t575,t576,t577,t578,t579,t580,t581,t582,t583,t584,t585,t586,t587,t588,t589,t590,t591,t592,t593,t594,t595,t596,t597,t598,t599
        ,t600,t601,t602,t603,t604,t605,t606,t607,t608,t609,t610,t611,t612,t613,t614,t615,t616,t617,t618,t619,t620,t621,t622,t623,t624,t625,t626,t627,t628,t629,t630,t631,t632,t633,t634,t635,t636,t637,t638,t639,t640,t641,t642,t643,t644,t645,t646,t647,t648,t649
        ,t650,t651,t652,t653,t654,t655,t656,t657,t658,t659,t660,t661,t662,t663,t664,t665,t666,t667,t668,t669,t670,t671,t672,t673,t674,t675,t676,t677,t678,t679,t680,t681,t682,t683,t684,t685,t686,t687,t688,t689,t690,t691,t692,t693,t694,t695,t696,t697,t698,t699
        ,t700,t701,t702,t703,t704,t705,t706,t707,t708,t709,t710,t711,t712,t713,t714,t715,t716,t717,t718,t719,t720,t721,t722,t723,t724,t725,t726,t727,t728,t729,t730,t731,t732,t733,t734,t735,t736,t737,t738,t739,t740,t741,t742,t743,t744,t745,t746,t747,t748,t749
        ,t750,t751,t752,t753,t754,t755,t756,t757,t758,t759,t760,t761,t762,t763,t764,t765,t766,t767,t768,t769,t770,t771,t772,t773,t774,t775,t776,t777,t778,t779,t780,t781,t782,t783,t784,t785,t786,t787,t788,t789,t790,t791,t792,t793,t794,t795,t796,t797,t798,t799
        ,t800,t801,t802,t803,t804,t805,t806,t807,t808,t809,t810,t811,t812,t813,t814,t815,t816,t817,t818,t819,t820,t821,t822,t823,t824,t825,t826,t827,t828,t829,t830,t831,t832,t833,t834,t835,t836,t837,t838,t839,t840,t841,t842,t843,t844,t845,t846,t847,t848,t849
        ,t850,t851,t852,t853,t854,t855,t856,t857,t858,t859,t860,t861,t862,t863,t864,t865,t866,t867,t868,t869,t870,t871,t872,t873,t874,t875,t876,t877,t878,t879,t880,t881,t882,t883,t884,t885,t886,t887,t888,t889,t890,t891,t892,t893,t894,t895,t896,t897,t898,t899
        ,t900,t901,t902,t903,t904,t905,t906,t907,t908,t909,t910,t911,t912,t913,t914,t915,t916,t917,t918,t919,t920,t921,t922,t923,t924,t925,t926,t927,t928,t929,t930,t931,t932,t933,t934,t935,t936,t937,t938,t939,t940,t941,t942,t943,t944,t945,t946,t947,t948,t949
        ,t950,t951,t952,t953,t954,t955,t956,t957,t958,t959,t960,t961,t962,t963,t964,t965,t966,t967,t968,t969,t970,t971,t972,t973,t974,t975,t976,t977,t978,t979,t980,t981,t982,t983,t984,t985,t986,t987,t988,t989,t990,t991,t992,t993,t994,t995,t996,t997,t998,t999
        ,t1000,t1001,t1002,t1003,t1004,t1005,t1006,t1007,t1008,t1009,t1010,t1011,t1012,t1013,t1014,t1015,t1016,t1017,t1018,t1019,t1020,t1021,t1022,t1023,t1024,t1025,t1026,t1027,t1028,t1029,t1030,t1031,t1032,t1033,t1034,t1035,t1036,t1037,t1038,t1039,t1040,t1041,t1042,t1043,t1044,t1045,t1046,t1047,t1048,t1049
        ,t1050,t1051,t1052,t1053,t1054,t1055,t1056,t1057,t1058,t1059,t1060,t1061,t1062,t1063,t1064,t1065,t1066,t1067,t1068,t1069,t1070,t1071,t1072,t1073,t1074,t1075,t1076,t1077,t1078,t1079,t1080,t1081,t1082,t1083,t1084,t1085,t1086,t1087,t1088,t1089,t1090,t1091,t1092,t1093,t1094,t1095,t1096,t1097,t1098,t1099
        ,t1100,t1101,t1102,t1103,t1104,t1105,t1106,t1107,t1108,t1109,t1110,t1111,t1112,t1113,t1114,t1115,t1116,t1117,t1118,t1119,t1120,t1121,t1122,t1123,t1124,t1125,t1126,t1127,t1128,t1129,t1130,t1131,t1132,t1133,t1134,t1135,t1136,t1137,t1138,t1139,t1140,t1141,t1142,t1143,t1144,t1145,t1146,t1147,t1148,t1149
        ,t1150,t1151,t1152,t1153,t1154,t1155,t1156,t1157,t1158,t1159,t1160,t1161,t1162,t1163,t1164,t1165,t1166,t1167,t1168,t1169,t1170,t1171,t1172,t1173,t1174,t1175,t1176,t1177,t1178,t1179,t1180,t1181,t1182,t1183,t1184,t1185,t1186,t1187,t1188,t1189,t1190,t1191,t1192,t1193,t1194,t1195,t1196,t1197,t1198,t1199
        ,t1200,t1201,t1202,t1203,t1204,t1205,t1206,t1207,t1208,t1209,t1210,t1211,t1212,t1213,t1214,t1215,t1216,t1217,t1218,t1219,t1220,t1221,t1222,t1223,t1224,t1225,t1226,t1227,t1228,t1229,t1230,t1231,t1232,t1233,t1234,t1235,t1236,t1237,t1238,t1239,t1240,t1241,t1242,t1243,t1244,t1245,t1246,t1247,t1248,t1249
        ,t1250,t1251,t1252,t1253,t1254,t1255,t1256,t1257,t1258,t1259,t1260,t1261,t1262,t1263,t1264,t1265,t1266,t1267,t1268,t1269,t1270,t1271,t1272,t1273,t1274,t1275,t1276,t1277,t1278,t1279,t1280,t1281,t1282,t1283,t1284,t1285,t1286,t1287,t1288,t1289,t1290,t1291,t1292,t1293,t1294,t1295,t1296,t1297,t1298,t1299
        ,t1300,t1301,t1302,t1303,t1304,t1305,t1306,t1307,t1308,t1309,t1310,t1311,t1312,t1313,t1314,t1315,t1316,t1317,t1318,t1319,t1320,t1321,t1322,t1323,t1324,t1325,t1326,t1327,t1328,t1329,t1330,t1331,t1332,t1333,t1334,t1335,t1336,t1337,t1338,t1339,t1340,t1341,t1342,t1343,t1344,t1345,t1346,t1347,t1348,t1349
        ,t1350,t1351,t1352,t1353,t1354,t1355,t1356,t1357,t1358,t1359,t1360,t1361,t1362,t1363,t1364,t1365,t1366,t1367,t1368,t1369,t1370,t1371,t1372,t1373,t1374,t1375,t1376,t1377,t1378,t1379,t1380,t1381,t1382,t1383,t1384,t1385,t1386,t1387,t1388,t1389,t1390,t1391,t1392,t1393,t1394,t1395,t1396,t1397,t1398,t1399
        ,t1400,t1401,t1402,t1403,t1404,t1405,t1406,t1407,t1408,t1409,t1410,t1411,t1412,t1413,t1414,t1415,t1416,t1417,t1418,t1419,t1420,t1421,t1422,t1423,t1424,t1425,t1426,t1427,t1428,t1429,t1430,t1431,t1432,t1433,t1434,t1435,t1436,t1437,t1438,t1439,t1440,t1441,t1442,t1443,t1444,t1445,t1446,t1447,t1448,t1449
        ,t1450,t1451,t1452,t1453,t1454,t1455,t1456,t1457,t1458,t1459,t1460,t1461,t1462,t1463,t1464,t1465,t1466,t1467,t1468,t1469,t1470,t1471,t1472,t1473,t1474,t1475,t1476,t1477,t1478,t1479,t1480,t1481,t1482,t1483,t1484,t1485,t1486,t1487,t1488,t1489,t1490,t1491,t1492,t1493,t1494,t1495,t1496,t1497,t1498,t1499
        ,t1500,t1501,t1502,t1503,t1504,t1505,t1506,t1507,t1508,t1509,t1510,t1511,t1512,t1513,t1514,t1515,t1516,t1517,t1518,t1519,t1520,t1521,t1522,t1523,t1524,t1525,t1526,t1527,t1528,t1529,t1530,t1531,t1532,t1533,t1534,t1535,t1536,t1537,t1538,t1539,t1540,t1541,t1542,t1543,t1544,t1545,t1546,t1547,t1548,t1549
        ,t1550,t1551,t1552,t1553,t1554,t1555,t1556,t1557,t1558,t1559,t1560,t1561,t1562,t1563,t1564,t1565,t1566,t1567,t1568,t1569,t1570,t1571,t1572,t1573,t1574,t1575,t1576,t1577,t1578,t1579,t1580,t1581,t1582,t1583,t1584,t1585,t1586,t1587,t1588,t1589,t1590,t1591,t1592,t1593,t1594,t1595,t1596,t1597,t1598,t1599
        ,t1600,t1601,t1602,t1603,t1604,t1605,t1606,t1607,t1608,t1609,t1610,t1611,t1612,t1613,t1614,t1615,t1616,t1617,t1618,t1619,t1620,t1621,t1622,t1623,t1624,t1625,t1626,t1627,t1628,t1629,t1630,t1631,t1632,t1633,t1634,t1635,t1636,t1637,t1638,t1639,t1640,t1641,t1642,t1643,t1644,t1645,t1646,t1647,t1648,t1649
        ,t1650,t1651,t1652,t1653,t1654,t1655,t1656,t1657,t1658,t1659,t1660,t1661,t1662,t1663,t1664,t1665,t1666,t1667,t1668,t1669,t1670,t1671,t1672,t1673,t1674,t1675,t1676,t1677,t1678,t1679,t1680,t1681,t1682,t1683,t1684,t1685,t1686,t1687,t1688,t1689,t1690,t1691,t1692,t1693,t1694,t1695,t1696,t1697,t1698,t1699
        ,t1700,t1701,t1702,t1703,t1704,t1705,t1706,t1707,t1708,t1709,t1710,t1711,t1712,t1713,t1714,t1715,t1716,t1717,t1718,t1719,t1720,t1721,t1722,t1723,t1724,t1725,t1726,t1727,t1728,t1729,t1730,t1731,t1732,t1733,t1734,t1735,t1736,t1737,t1738,t1739,t1740,t1741,t1742,t1743,t1744,t1745,t1746,t1747,t1748,t1749
        ,t1750,t1751,t1752,t1753,t1754,t1755,t1756,t1757,t1758,t1759,t1760,t1761,t1762,t1763,t1764,t1765,t1766,t1767,t1768,t1769,t1770,t1771,t1772,t1773,t1774,t1775,t1776,t1777,t1778,t1779,t1780,t1781,t1782,t1783,t1784,t1785,t1786,t1787,t1788,t1789,t1790,t1791,t1792,t1793,t1794,t1795,t1796,t1797,t1798,t1799
        ,t1800,t1801,t1802,t1803,t1804,t1805,t1806,t1807,t1808,t1809,t1810,t1811,t1812,t1813,t1814,t1815,t1816,t1817,t1818,t1819,t1820,t1821,t1822,t1823,t1824,t1825,t1826,t1827,t1828,t1829,t1830,t1831,t1832,t1833,t1834,t1835,t1836,t1837,t1838,t1839,t1840,t1841,t1842,t1843,t1844,t1845,t1846,t1847,t1848,t1849
        ,t1850,t1851,t1852,t1853,t1854,t1855,t1856,t1857,t1858,t1859,t1860,t1861,t1862,t1863,t1864,t1865,t1866,t1867,t1868,t1869,t1870,t1871,t1872,t1873,t1874,t1875,t1876,t1877,t1878,t1879,t1880,t1881,t1882,t1883,t1884,t1885,t1886,t1887,t1888,t1889,t1890,t1891,t1892,t1893,t1894,t1895,t1896,t1897,t1898,t1899
        ,t1900,t1901,t1902,t1903,t1904,t1905,t1906,t1907,t1908,t1909,t1910,t1911,t1912,t1913,t1914,t1915,t1916,t1917,t1918,t1919,t1920,t1921,t1922,t1923,t1924,t1925,t1926,t1927,t1928,t1929,t1930,t1931,t1932,t1933,t1934,t1935,t1936,t1937,t1938,t1939,t1940,t1941,t1942,t1943,t1944,t1945,t1946,t1947,t1948,t1949
        ,t1950,t1951,t1952,t1953,t1954,t1955,t1956,t1957,t1958,t1959,t1960,t1961,t1962,t1963,t1964,t1965,t1966,t1967,t1968,t1969,t1970,t1971,t1972,t1973,t1974,t1975,t1976,t1977,t1978,t1979,t1980,t1981,t1982,t1983,t1984,t1985,t1986,t1987,t1988,t1989,t1990,t1991,t1992,t1993,t1994,t1995,t1996,t1997,t1998,t1999
        ,t2000,t2001,t2002,t2003,t2004,t2005,t2006,t2007,t2008,t2009,t2010,t2011,t2012,t2013,t2014,t2015,t2016,t2017,t2018,t2019,t2020,t2021,t2022,t2023,t2024,t2025,t2026,t2027,t2028,t2029,t2030,t2031,t2032,t2033,t2034,t2035,t2036,t2037,t2038,t2039,t2040,t2041,t2042,t2043,t2044,t2045,t2046,t2047,t2048,t2049
        ,t2050,t2051,t2052,t2053,t2054,t2055,t2056,t2057,t2058,t2059,t2060,t2061,t2062,t2063,t2064,t2065,t2066,t2067,t2068,t2069,t2070,t2071,t2072,t2073,t2074,t2075,t2076,t2077,t2078,t2079,t2080,t2081,t2082,t2083,t2084,t2085,t2086,t2087,t2088,t2089,t2090,t2091,t2092,t2093,t2094,t2095,t2096,t2097,t2098,t2099
        ,t2100,t2101,t2102,t2103,t2104,t2105,t2106,t2107,t2108,t2109,t2110,t2111,t2112,t2113,t2114,t2115,t2116,t2117,t2118,t2119,t2120,t2121,t2122,t2123,t2124,t2125,t2126,t2127,t2128,t2129,t2130,t2131,t2132,t2133,t2134,t2135,t2136,t2137,t2138,t2139,t2140,t2141,t2142,t2143,t2144,t2145,t2146,t2147,t2148,t2149
        ,t2150,t2151,t2152,t2153,t2154,t2155,t2156,t2157,t2158,t2159,t2160,t2161,t2162,t2163,t2164,t2165,t2166,t2167,t2168,t2169,t2170,t2171,t2172,t2173,t2174,t2175,t2176,t2177,t2178,t2179,t2180,t2181,t2182,t2183,t2184,t2185,t2186,t2187,t2188,t2189,t2190,t2191,t2192,t2193,t2194,t2195,t2196,t2197,t2198,t2199
        ,t2200,t2201,t2202,t2203,t2204,t2205,t2206,t2207,t2208,t2209,t2210,t2211,t2212,t2213,t2214,t2215,t2216,t2217,t2218,t2219,t2220,t2221,t2222,t2223,t2224,t2225,t2226,t2227,t2228,t2229,t2230,t2231,t2232,t2233,t2234,t2235,t2236,t2237,t2238,t2239,t2240,t2241,t2242,t2243,t2244,t2245,t2246,t2247,t2248,t2249
        ,t2250,t2251,t2252,t2253,t2254,t2255,t2256,t2257,t2258,t2259,t2260,t2261,t2262,t2263,t2264,t2265,t2266,t2267,t2268,t2269,t2270,t2271,t2272,t2273,t2274,t2275,t2276,t2277,t2278,t2279,t2280,t2281,t2282,t2283,t2284,t2285,t2286,t2287,t2288,t2289,t2290,t2291,t2292,t2293,t2294,t2295,t2296,t2297,t2298,t2299
        ,t2300,t2301,t2302,t2303,t2304,t2305,t2306,t2307,t2308,t2309,t2310,t2311,t2312,t2313,t2314,t2315,t2316,t2317,t2318,t2319,t2320,t2321,t2322,t2323,t2324,t2325,t2326,t2327,t2328,t2329,t2330,t2331,t2332,t2333,t2334,t2335,t2336,t2337,t2338,t2339,t2340,t2341,t2342,t2343,t2344,t2345,t2346,t2347,t2348,t2349
        ,t2350,t2351,t2352,t2353,t2354,t2355,t2356,t2357,t2358,t2359,t2360,t2361,t2362,t2363,t2364,t2365,t2366,t2367,t2368,t2369,t2370,t2371,t2372,t2373,t2374,t2375,t2376,t2377,t2378,t2379,t2380,t2381,t2382,t2383,t2384,t2385,t2386,t2387,t2388,t2389,t2390,t2391,t2392,t2393,t2394,t2395,t2396,t2397,t2398,t2399
        ,t2400,t2401,t2402,t2403,t2404,t2405,t2406,t2407,t2408,t2409,t2410,t2411,t2412,t2413,t2414,t2415,t2416,t2417,t2418,t2419,t2420,t2421,t2422,t2423,t2424,t2425,t2426,t2427,t2428,t2429,t2430,t2431,t2432,t2433,t2434,t2435,t2436,t2437,t2438,t2439,t2440,t2441,t2442,t2443,t2444,t2445,t2446,t2447,t2448,t2449
        ,t2450,t2451,t2452,t2453,t2454,t2455,t2456,t2457,t2458,t2459,t2460,t2461,t2462,t2463,t2464,t2465,t2466,t2467,t2468,t2469,t2470,t2471,t2472,t2473,t2474,t2475,t2476,t2477,t2478,t2479,t2480,t2481,t2482,t2483,t2484,t2485,t2486,t2487,t2488,t2489,t2490,t2491,t2492,t2493,t2494,t2495,t2496,t2497,t2498,t2499
        ,t2500;



            if( annulusScaling <= 0. )
            { // -- first time: compute the norm of the solution so we can scale to be size one
                Real t0=0.; 
                    Real ur,vr,pr; 
                    Real u1Max=0., u2Max=0., u3Max=0., pMax=0.; 
                    uNorm=0; 
                    if( numberOfDimensions==2 )
                    {
                        FOR_3D(i1,i2,i3,I1,I2,I3)
                        {
                            Real x = vertex(i1,i2,i3,0);
                            Real y = vertex(i1,i2,i3,1);
                            Real r = sqrt( x*x + y*y );
                            Real theta = atan2(y,x); 
                            Real cosTheta = x/r;
                            Real sinTheta = y/r;
                            t1 = 0.1e1 / r;
                            t2 = jn(n, omega);
                            t3 = n * n;
                            t4 = t3 * t3;
                            t5 = t4 * t3;
                            t8 = t4 * t2;
                            t10 = t4 * n;
                            t13 = t3 * n;
                            t14 = t13 * t2;
                            t16 = omega / 0.2e1;
                            t17 = jn(n, t16);
                            t18 = pow(0.1e1 / 0.2e1, n);
                            t19 = t18 * t17;
                            t20 = omega * omega;
                            t21 = t20 * t20;
                            t22 = t21 * n;
                            t25 = n + 0.1e1;
                            t26 = jn(t25, t16);
                            t27 = t18 * t26;
                            t28 = omega * t4;
                            t31 = t20 * omega;
                            t32 = t31 * t3;
                            t35 = t20 * t3;
                            t38 = omega * t13;
                            t41 = t31 * n;
                            t44 = 0.1e1 / t18;
                            t45 = t17 * t44;
                            t48 = t26 * t44;
                            t53 = t3 * t26;
                            t60 = t20 * t4;
                            t63 = t21 * t3;
                            t66 = omega * t10;
                            t69 = t31 * t13;
                            t74 = 0.64e2 * t18 * omega * t53 - 0.64e2 * t44 * omega * t53 + 0.256e3 * t10 * t2 + 0.14e2 * t22 * t19 - 0.32e2 * t35 * t19 - 0.256e3 * t5 * t2 - 0.2e1 * t22 * t45 - 0.64e2 * t28 * t27 + 0.8e1 * t32 * t27 + 0.64e2 * t38 * t27 - 0.24e2 * t41 * t27 + 0.128e3 * t35 * t45 + 0.64e2 * t38 * t48 - 0.24e2 * t41 * t48 - 0.160e3 * t60 * t45 + 0.18e2 * t63 * t45 - 0.64e2 * t66 * t48 + 0.32e2 * t69 * t48 - 0.256e3 * t14 + 0.256e3 * t8;
                            t91 = t21 * t20;
                            t93 = t21 * omega;
                            t98 = jn(t25, omega);
                            t99 = t13 * t98;
                            t122 = 0.32e2 * t44 * t20 * t13 * t17 + 0.256e3 * omega * t10 * t98 - 0.256e3 * t44 * t10 * t17 + 0.256e3 * t44 * t5 * t17 - 0.256e3 * omega * t99 + 0.256e3 * t13 * t45 - 0.32e2 * t20 * t14 + 0.32e2 * t60 * t19 - 0.18e2 * t63 * t19 + t91 * t19 + 0.32e2 * t20 * t8 - 0.64e2 * t66 * t27 + 0.32e2 * t69 * t27 - 0.4e1 * t93 * t27 + 0.64e2 * t28 * t48 - 0.32e2 * t31 * t99 - 0.8e1 * t32 * t48 - 0.256e3 * t4 * t45 - t91 * t45 + 0.4e1 * t93 * t48;
                            t124 = 0.1e1 / (t74 + t122);
                            t125 = t124 * t1;
                            t126 = pow(r, n);
                            t127 = 0.1e1 / t126;
                            t128 = t127 * t18;
                            t129 = omega * t0;
                            t130 = sin(t129);
                            t132 = t130 * t128 * t125;
                            t133 = n * theta;
                            t134 = sin(t133);
                            t135 = yn(n, t16);
                            t136 = t135 * t134;
                            t137 = t98 * t136;
                            t138 = t3 * c4;
                            t139 = t31 * t138;
                            t140 = t139 * t137;
                            t143 = t2 * t136;
                            t144 = t13 * c4;
                            t145 = t20 * t144;
                            t146 = t145 * t143;
                            t149 = n * c4;
                            t150 = t21 * t149;
                            t151 = t150 * t143;
                            t154 = yn(t25, t16);
                            t155 = t154 * t134;
                            t156 = t98 * t155;
                            t157 = t150 * t156;
                            t160 = t17 * t134;
                            t161 = yn(t25, omega);
                            t162 = t161 * t160;
                            t163 = t139 * t162;
                            t166 = t126 * t44;
                            t168 = t130 * t166 * t125;
                            t171 = cos(t129);
                            t173 = t171 * t128 * t125;
                            t174 = cos(t133);
                            t175 = t26 * t174;
                            t176 = yn(n, omega);
                            t177 = t176 * t175;
                            t178 = omega * t138;
                            t179 = t178 * t177;
                            t183 = t171 * t166 * t125;
                            t184 = t154 * t174;
                            t185 = t98 * t184;
                            t186 = t20 * t138;
                            t187 = t186 * t185;
                            t190 = t130 * t126;
                            t191 = t134 * t190;
                            t192 = t191 * t125;
                            t193 = t154 * t17;
                            t194 = t139 * t193;
                            t197 = t171 * t126;
                            t198 = t174 * t197;
                            t199 = t198 * t125;
                            t200 = t135 * t26;
                            t201 = t178 * t200;
                            t204 = t171 * t18;
                            t205 = omega * r;
                            t206 = jn(n, t205);
                            t208 = t206 * t204 * t125;
                            t209 = t135 * t174;
                            t210 = t186 * t209;
                            t213 = t171 * t44;
                            t215 = t206 * t213 * t125;
                            t216 = t4 * c4;
                            t217 = t20 * t216;
                            t218 = t217 * t209;
                            t221 = t21 * t138;
                            t222 = t221 * t209;
                            t225 = t10 * c4;
                            t226 = omega * t225;
                            t227 = t226 * t184;
                            t230 = t130 * t18;
                            t232 = t206 * t230 * t125;
                            t233 = t221 * t136;
                            t236 = t171 * t127;
                            t237 = t174 * t236;
                            t238 = t237 * t125;
                            t239 = omega * t144;
                            t240 = t239 * t200;
                            t243 = t31 * t149;
                            t244 = t243 * t200;
                            t247 = t145 * t209;
                            t250 = -0.32e2 * t140 * t132 + 0.64e2 * t146 * t132 - 0.2e1 * t151 * t132 + 0.8e1 * t157 * t132 + 0.32e2 * t163 * t132 - 0.2e1 * t151 * t168 - 0.64e2 * t179 * t173 - 0.64e2 * t187 * t183 + 0.8e1 * t194 * t192 - 0.64e2 * t201 * t199 + 0.32e2 * t210 * t208 + 0.18e2 * t222 * t208 + 0.64e2 * t227 * t208 + 0.160e3 * t218 * t215 - 0.32e2 * t247 * t215 + 0.18e2 * t233 * t232 - 0.64e2 * t240 * t238 + 0.32e2 * t244 * t238;
                            t251 = t150 * t209;
                            t254 = t130 * t127;
                            t255 = t134 * t254;
                            t256 = t255 * t125;
                            t257 = t178 * t193;
                            t262 = t226 * t193;
                            t265 = t31 * t144;
                            t266 = t265 * t193;
                            t269 = t130 * t44;
                            t271 = t206 * t269 * t125;
                            t272 = t217 * t136;
                            t281 = t239 * t184;
                            t284 = t243 * t184;
                            t289 = t239 * t193;
                            t292 = t243 * t193;
                            t295 = t265 * t184;
                            t298 = t17 * t174;
                            t299 = c4 * t176;
                            t300 = t5 * t299;
                            t301 = t300 * t298;
                            t304 = c4 * t2;
                            t305 = t5 * t304;
                            t306 = t305 * t209;
                            t309 = t10 * t299;
                            t310 = t309 * t298;
                            t313 = t10 * t304;
                            t314 = t313 * t209;
                            t317 = -0.256e3 * t301 * t183 + 0.256e3 * t306 * t183 + 0.256e3 * t310 * t183 - 0.256e3 * t314 * t183 + 0.64e2 * t201 * t238 + 0.64e2 * t201 * t256 - 0.32e2 * t218 * t208 - 0.32e2 * t295 * t208 - 0.128e3 * t210 * t215 + 0.2e1 * t251 * t215 - 0.64e2 * t281 * t215 + 0.24e2 * t284 * t215 - 0.18e2 * t233 * t271 + 0.64e2 * t289 * t238 - 0.32e2 * t292 * t238 - 0.64e2 * t257 * t256 - 0.64e2 * t262 * t256 + 0.40e2 * t266 * t256 + 0.160e3 * t272 * t271;
                            t319 = t4 * t299;
                            t320 = t319 * t298;
                            t323 = t4 * t304;
                            t324 = t323 * t209;
                            t327 = t13 * t299;
                            t328 = t327 * t298;
                            t331 = t13 * t304;
                            t332 = t331 * t209;
                            t339 = t176 * t160;
                            t340 = t145 * t339;
                            t343 = t150 * t339;
                            t346 = t26 * t134;
                            t347 = t161 * t346;
                            t348 = t150 * t347;
                            t351 = t161 * t298;
                            t352 = t265 * t351;
                            t355 = t93 * t149;
                            t356 = t355 * t351;
                            t359 = t176 * t298;
                            t360 = t217 * t359;
                            t363 = t221 * t359;
                            t366 = t161 * t175;
                            t367 = t217 * t366;
                            t370 = t2 * t155;
                            t371 = t239 * t370;
                            t374 = t243 * t370;
                            t377 = t176 * t346;
                            t378 = omega * t216;
                            t379 = t378 * t377;
                            t382 = t139 * t377;
                            t385 = -0.64e2 * t340 * t132 + 0.2e1 * t343 * t132 - 0.8e1 * t348 * t132 + 0.64e2 * t379 * t132 + 0.8e1 * t382 * t132 + 0.64e2 * t371 * t168 + 0.8e1 * t374 * t168 - 0.32e2 * t352 * t173 + 0.2e1 * t356 * t173 + 0.32e2 * t360 * t173 - 0.2e1 * t363 * t173 + 0.64e2 * t367 * t173 + 0.256e3 * t320 * t183 - 0.256e3 * t324 * t183 - 0.256e3 * t328 * t183 + 0.256e3 * t332 * t183 + 0.64e2 * t257 * t199 - 0.18e2 * t222 * t215;
                            t386 = t217 * t339;
                            t389 = t221 * t339;
                            t392 = t217 * t347;
                            t395 = t226 * t377;
                            t398 = t226 * t137;
                            t401 = t265 * t137;
                            t404 = t243 * t377;
                            t407 = t355 * t137;
                            t410 = t217 * t143;
                            t413 = t221 * t143;
                            t416 = t217 * t156;
                            t419 = t226 * t370;
                            t422 = t265 * t162;
                            t425 = t355 * t162;
                            t430 = t2 * t184;
                            t431 = t139 * t430;
                            t434 = t239 * t351;
                            t437 = t186 * t366;
                            t440 = t239 * t177;
                            t443 = 0.32e2 * t386 * t132 - 0.8e1 * t404 * t132 - 0.32e2 * t422 * t132 + 0.2e1 * t425 * t132 + 0.32e2 * t386 * t168 - 0.2e1 * t389 * t168 - 0.64e2 * t392 * t168 + 0.64e2 * t395 * t168 - 0.256e3 * t398 * t168 + 0.32e2 * t401 * t168 - 0.2e1 * t407 * t168 - 0.32e2 * t410 * t168 + 0.2e1 * t413 * t168 + 0.64e2 * t416 * t168 - 0.64e2 * t419 * t168 - 0.8e1 * t431 * t173 - 0.256e3 * t434 * t183 + 0.64e2 * t437 * t183 - 0.64e2 * t440 * t183;
                            t446 = t243 * t177;
                            t449 = t98 * t209;
                            t450 = t239 * t449;
                            t457 = t178 * t155;
                            t460 = t145 * t136;
                            t463 = t150 * t136;
                            t466 = t378 * t155;
                            t469 = t139 * t155;
                            t472 = t331 * t136;
                            t475 = t300 * t160;
                            t484 = t319 * t160;
                            t487 = t323 * t136;
                            t490 = t305 * t136;
                            t495 = 0.256e3 * t472 * t168 - 0.256e3 * t475 * t168 + 0.256e3 * t484 * t168 - 0.256e3 * t487 * t168 + 0.256e3 * t490 * t168 - 0.8e1 * t446 * t183 + 0.256e3 * t450 * t183 - 0.64e2 * t281 * t208 - 0.14e2 * t463 * t232 + 0.64e2 * t466 * t232 - 0.8e1 * t469 * t232 - 0.64e2 * t240 * t256 + 0.32e2 * t244 * t256 + 0.64e2 * t457 * t271 - 0.32e2 * t460 * t271 + 0.2e1 * t463 * t271 - 0.64e2 * t466 * t271 + 0.8e1 * t469 * t271;
                            t512 = t355 * t193;
                            t517 = t378 * t200;
                            t520 = t139 * t200;
                            t523 = t226 * t200;
                            t526 = t265 * t200;
                            t529 = t355 * t200;
                            t532 = t378 * t193;
                            t541 = 0.8e1 * t194 * t199 - 0.8e1 * t194 * t238 - 0.64e2 * t240 * t199 + 0.32e2 * t244 * t199 + 0.64e2 * t289 * t199 - 0.32e2 * t292 * t199 + 0.64e2 * t517 * t199 - 0.8e1 * t520 * t199 - 0.64e2 * t532 * t199 + 0.24e2 * t284 * t208 - 0.32e2 * t295 * t215 - 0.64e2 * t262 * t238 + 0.40e2 * t266 * t238 - 0.4e1 * t512 * t238 - 0.64e2 * t517 * t238 + 0.8e1 * t520 * t238 + 0.64e2 * t523 * t238 - 0.40e2 * t526 * t238 + 0.4e1 * t529 * t238;
                            t543 = t178 * t184;
                            t546 = t378 * t184;
                            t549 = t139 * t184;
                            t552 = t226 * t155;
                            t555 = t265 * t155;
                            t562 = t230 * t125;
                            t563 = t134 * t206;
                            t564 = c4 * t135;
                            t565 = t91 * t564;
                            t566 = t565 * t563;
                            t568 = t269 * t125;
                            t571 = t206 * t171 * t125;
                            t572 = t161 * t174;
                            t573 = t226 * t572;
                            t576 = t204 * t125;
                            t577 = t174 * t206;
                            t578 = c4 * t154;
                            t579 = t93 * t578;
                            t580 = t579 * t577;
                            t583 = t265 * t572;
                            t586 = t176 * t174;
                            t587 = t217 * t586;
                            t590 = t239 * t572;
                            t593 = t565 * t577;
                            t595 = t213 * t125;
                            t596 = t5 * t564;
                            t597 = t596 * t577;
                            t601 = t10 * t564;
                            t602 = t601 * t577;
                            t607 = -0.64e2 * t543 * t208 - 0.64e2 * t546 * t215 + 0.8e1 * t549 * t215 + 0.64e2 * t552 * t232 - 0.32e2 * t555 * t232 + 0.64e2 * t289 * t256 - 0.32e2 * t292 * t256 - t566 * t562 + t566 * t568 - 0.256e3 * t573 * t571 + 0.32e2 * t583 * t571 - 0.32e2 * t587 * t571 + 0.256e3 * t590 * t571 + 0.4e1 * t580 * t576 - t593 * t576 - 0.4e1 * t580 * t595 + t593 * t595 - 0.256e3 * t597 * t595 + 0.256e3 * t602 * t595;
                            t608 = t4 * t564;
                            t609 = t608 * t577;
                            t612 = t327 * t160;
                            t617 = t186 * t136;
                            t620 = t239 * t155;
                            t623 = t243 * t155;
                            t652 = -0.256e3 * t612 * t168 + 0.40e2 * t266 * t192 - 0.8e1 * t194 * t256 - 0.4e1 * t512 * t199 + 0.64e2 * t523 * t199 - 0.40e2 * t526 * t199 + 0.4e1 * t529 * t199 - 0.32e2 * t272 * t232 - 0.4e1 * t512 * t256 - 0.64e2 * t517 * t256 + 0.8e1 * t520 * t256 + 0.64e2 * t523 * t256 - 0.40e2 * t526 * t256 + 0.4e1 * t529 * t256 + 0.64e2 * t532 * t256 - 0.128e3 * t617 * t271 - 0.64e2 * t620 * t271 + 0.24e2 * t623 * t271 + 0.256e3 * t609 * t595;
                            t680 = t309 * t160;
                            t683 = t313 * t136;
                            t694 = 0.256e3 * t680 * t168 - 0.256e3 * t683 * t168 - 0.64e2 * t240 * t192 + 0.64e2 * t289 * t192 - 0.32e2 * t292 * t192 + 0.64e2 * t517 * t192 - 0.8e1 * t520 * t192 - 0.40e2 * t526 * t192 + 0.4e1 * t529 * t192 - 0.64e2 * t532 * t192 + 0.64e2 * t227 * t215 + 0.64e2 * t543 * t215 + 0.32e2 * t617 * t232 - 0.64e2 * t620 * t232 + 0.24e2 * t623 * t232 - 0.64e2 * t257 * t238 + 0.64e2 * t552 * t271 - 0.32e2 * t555 * t271;
                            t721 = t225 * t586;
                            t725 = t206 * t130 * t125;
                            t726 = t176 * t134;
                            t727 = t144 * t726;
                            t730 = t5 * c4;
                            t731 = t730 * t726;
                            t734 = t225 * t726;
                            t737 = t216 * t726;
                            t740 = t730 * t586;
                            t743 = -0.64e2 * t201 * t192 + 0.32e2 * t244 * t192 + 0.64e2 * t257 * t192 - 0.64e2 * t262 * t192 - 0.4e1 * t512 * t192 + 0.64e2 * t523 * t192 - 0.64e2 * t262 * t199 + 0.40e2 * t266 * t199 - 0.14e2 * t251 * t208 + 0.64e2 * t546 * t208 - 0.8e1 * t549 * t208 - 0.64e2 * t457 * t232 + 0.64e2 * t532 * t238 - 0.256e3 * t721 * t571 + 0.256e3 * t740 * t571 + 0.256e3 * t727 * t725 + 0.256e3 * t731 * t725 - 0.256e3 * t734 * t725 - 0.256e3 * t737 * t725;
                            t745 = t144 * t586;
                            t748 = t216 * t586;
                            t751 = t186 * t143;
                            t754 = t186 * t156;
                            t773 = t178 * t377;
                            t776 = t239 * t162;
                            t779 = t186 * t347;
                            t782 = t239 * t377;
                            t787 = t239 * t137;
                            t790 = -0.64e2 * t371 * t132 + 0.8e1 * t374 * t132 - 0.32e2 * t751 * t132 + 0.64e2 * t754 * t132 + 0.8e1 * t157 * t168 - 0.32e2 * t340 * t168 + 0.2e1 * t343 * t168 - 0.8e1 * t348 * t168 - 0.64e2 * t379 * t168 + 0.8e1 * t382 * t168 - 0.8e1 * t404 * t168 + 0.64e2 * t773 * t168 - 0.256e3 * t776 * t168 + 0.64e2 * t779 * t168 - 0.64e2 * t782 * t168 + 0.256e3 * t787 * t168 + 0.256e3 * t745 * t571 - 0.256e3 * t748 * t571;
                            t795 = t150 * t185;
                            t798 = t378 * t430;
                            t813 = t226 * t177;
                            t816 = t265 * t449;
                            t819 = t355 * t449;
                            t822 = t243 * t430;
                            t825 = t2 * t209;
                            t826 = t145 * t825;
                            t829 = t150 * t825;
                            t834 = t186 * t825;
                            t837 = t139 * t351;
                            t840 = -0.2e1 * t389 * t132 + 0.64e2 * t392 * t132 - 0.64e2 * t395 * t132 + 0.32e2 * t401 * t132 - 0.2e1 * t407 * t132 - 0.32e2 * t410 * t132 + 0.32e2 * t146 * t168 - 0.64e2 * t754 * t168 - 0.8e1 * t446 * t173 + 0.8e1 * t795 * t173 - 0.64e2 * t798 * t173 - 0.64e2 * t813 * t173 + 0.32e2 * t816 * t173 - 0.2e1 * t819 * t173 + 0.8e1 * t822 * t173 + 0.64e2 * t826 * t173 - 0.2e1 * t829 * t173 - 0.32e2 * t834 * t173 + 0.32e2 * t837 * t173;
                            t843 = t145 * t359;
                            t846 = t150 * t359;
                            t849 = t150 * t366;
                            t852 = t378 * t177;
                            t855 = t139 * t177;
                            t858 = t139 * t449;
                            t861 = t217 * t185;
                            t864 = t226 * t430;
                            t881 = t239 * t430;
                            t886 = t178 * t430;
                            t889 = -0.64e2 * t843 * t173 + 0.2e1 * t846 * t173 - 0.8e1 * t849 * t173 + 0.64e2 * t852 * t173 + 0.8e1 * t855 * t173 - 0.32e2 * t858 * t173 - 0.64e2 * t861 * t173 + 0.64e2 * t864 * t173 + 0.64e2 * t886 * t173 + 0.8e1 * t822 * t183 + 0.32e2 * t826 * t183 - 0.2e1 * t829 * t183 - 0.32e2 * t843 * t183 + 0.2e1 * t846 * t183 - 0.8e1 * t849 * t183 - 0.64e2 * t852 * t183 + 0.8e1 * t855 * t183 + 0.64e2 * t881 * t183;
                            t896 = t217 * t825;
                            t899 = t221 * t825;
                            t902 = t378 * t370;
                            t905 = t139 * t370;
                            t908 = t186 * t339;
                            t917 = t178 * t370;
                            t926 = t186 * t359;
                            t933 = t13 * t564;
                            t934 = t933 * t577;
                            t937 = 0.2e1 * t413 * t132 - 0.64e2 * t416 * t132 + 0.64e2 * t419 * t132 - 0.64e2 * t773 * t132 - 0.64e2 * t779 * t132 + 0.64e2 * t782 * t132 - 0.64e2 * t902 * t132 - 0.8e1 * t905 * t132 + 0.32e2 * t908 * t132 + 0.64e2 * t917 * t132 + 0.64e2 * t902 * t168 - 0.8e1 * t905 * t168 - 0.64e2 * t917 * t168 - 0.64e2 * t437 * t173 + 0.64e2 * t440 * t173 - 0.32e2 * t896 * t173 + 0.2e1 * t899 * t173 + 0.32e2 * t926 * t173 - 0.256e3 * t934 * t595;
                            t939 = t145 * t586;
                            t942 = t161 * t134;
                            t943 = t239 * t942;
                            t946 = t145 * t726;
                            t949 = t226 * t942;
                            t952 = t265 * t942;
                            t955 = t217 * t726;
                            t958 = t601 * t563;
                            t961 = t579 * t563;
                            t964 = t596 * t563;
                            t967 = t608 * t563;
                            t970 = t933 * t563;
                            t975 = t171 * t1;
                            t976 = yn(n, t205);
                            t978 = c4 * t174 * t976;
                            t980 = t130 * t1;
                            t982 = c4 * t134 * t976;
                            t994 = 0.64e2 * t187 * t173 - 0.64e2 * t881 * t173 - 0.8e1 * t431 * t183 + 0.8e1 * t795 * t183 + 0.64e2 * t798 * t183 + 0.4e1 * t961 * t562 + 0.256e3 * t958 * t568 - 0.4e1 * t961 * t568 - 0.256e3 * t964 * t568 + 0.256e3 * t967 * t568 - 0.256e3 * t970 * t568 + 0.32e2 * t939 * t571 + 0.256e3 * t943 * t725 + 0.32e2 * t946 * t725 - 0.256e3 * t949 * t725 + 0.32e2 * t952 * t725 - 0.32e2 * t955 * t725 + t978 * t975 + t982 * t980;
                            t995 = t226 * t162;
                            t1006 = t226 * t351;
                            t1021 = t226 * t449;
                            t1036 = 0.256e3 * t1006 * t183 - 0.256e3 * t1021 * t183 - 0.32e2 * t422 * t168 + 0.2e1 * t425 * t168 + 0.256e3 * t995 * t168 + 0.64e2 * t179 * t183 - 0.32e2 * t352 * t183 + 0.2e1 * t356 * t183 + 0.32e2 * t360 * t183 - 0.2e1 * t363 * t183 - 0.64e2 * t367 * t183 + 0.64e2 * t813 * t183 + 0.32e2 * t816 * t183 - 0.2e1 * t819 * t183 + 0.64e2 * t861 * t183 - 0.64e2 * t864 * t183 - 0.64e2 * t886 * t183 - 0.32e2 * t896 * t183 + 0.2e1 * t899 * t183;
                            ur = t250 + t317 + t385 + t443 + t495 + t541 + t607 + t652 + t694 + t743 + t790 + t840 + t889 + t937 + t994 + t1036;
                            t1053 = t134 * t197 * t125;
                            t1061 = t174 * t254 * t125;
                            t1066 = -0.8e1 * t194 * t1053 - 0.64e2 * t517 * t1053 + 0.64e2 * t532 * t1053 + 0.64e2 * t240 * t1061 + 0.32e2 * t292 * t1061 - 0.64e2 * t552 * t208 + 0.32e2 * t555 * t208 - 0.32e2 * t218 * t232 + 0.18e2 * t222 * t232 + 0.64e2 * t227 * t232 - 0.32e2 * t295 * t232;
                            t1074 = t134 * t236 * t125;
                            t1093 = -0.32e2 * t244 * t1061 - 0.64e2 * t257 * t1074 - 0.256e3 * t314 * t168 - 0.256e3 * t328 * t168 - 0.256e3 * t434 * t168 - 0.256e3 * t472 * t183 + 0.256e3 * t612 * t183 - 0.128e3 * t210 * t271 - 0.64e2 * t457 * t215 - 0.64e2 * t281 * t271 + 0.24e2 * t284 * t271 + 0.64e2 * t543 * t271;
                            t1116 = t174 * t190 * t125;
                            t1121 = -0.32e2 * t244 * t1053 - 0.64e2 * t262 * t1074 + 0.40e2 * t266 * t1074 - 0.4e1 * t512 * t1074 + 0.64e2 * t523 * t1074 - 0.40e2 * t526 * t1074 + 0.4e1 * t529 * t1074 + 0.64e2 * t289 * t1116 - 0.32e2 * t292 * t1116 + 0.256e3 * t310 * t168 + 0.256e3 * t332 * t168 - 0.32e2 * t295 * t271;
                            t1146 = -0.64e2 * t240 * t1116 + 0.32e2 * t244 * t1116 - 0.32e2 * t837 * t132 + 0.64e2 * t843 * t132 - 0.2e1 * t846 * t132 + 0.8e1 * t849 * t132 - 0.64e2 * t852 * t132 + 0.32e2 * t163 * t173 - 0.64e2 * t340 * t173 + 0.2e1 * t343 * t173 - 0.8e1 * t348 * t173 - 0.8e1 * t905 * t173;
                            t1173 = 0.64e2 * t179 * t132 - 0.64e2 * t886 * t132 + 0.64e2 * t179 * t168 - 0.8e1 * t849 * t168 - 0.64e2 * t852 * t168 - 0.64e2 * t886 * t168 - 0.32e2 * t896 * t168 + 0.2e1 * t899 * t168 - 0.64e2 * t416 * t173 - 0.64e2 * t773 * t173 + 0.64e2 * t917 * t173 + 0.64e2 * t917 * t183;
                            t1198 = 0.32e2 * t352 * t132 - 0.2e1 * t356 * t132 - 0.32e2 * t360 * t132 + 0.2e1 * t363 * t132 - 0.64e2 * t367 * t132 + 0.64e2 * t813 * t132 - 0.32e2 * t816 * t132 + 0.32e2 * t386 * t173 - 0.2e1 * t389 * t173 + 0.64e2 * t392 * t173 - 0.64e2 * t395 * t173 + 0.64e2 * t419 * t173;
                            t1224 = -0.64e2 * t187 * t132 + 0.8e1 * t446 * t132 - 0.8e1 * t822 * t132 + 0.32e2 * t834 * t132 + 0.64e2 * t881 * t132 + 0.32e2 * t360 * t168 - 0.2e1 * t363 * t168 - 0.64e2 * t367 * t168 + 0.32e2 * t401 * t173 - 0.2e1 * t407 * t173 + 0.32e2 * t340 * t183 - 0.64e2 * t773 * t183;
                            t1249 = -0.256e3 * t1021 * t168 + 0.2e1 * t819 * t132 + 0.64e2 * t861 * t132 - 0.64e2 * t864 * t132 + 0.32e2 * t896 * t132 - 0.2e1 * t899 * t132 + 0.64e2 * t813 * t168 + 0.8e1 * t855 * t168 - 0.32e2 * t386 * t183 + 0.2e1 * t389 * t183 + 0.64e2 * t392 * t183 - 0.64e2 * t395 * t183;
                            t1275 = 0.8e1 * t157 * t173 + 0.8e1 * t795 * t168 + 0.64e2 * t798 * t168 + 0.32e2 * t826 * t168 - 0.2e1 * t829 * t168 - 0.32e2 * t422 * t173 + 0.2e1 * t425 * t173 - 0.64e2 * t779 * t173 + 0.64e2 * t782 * t173 - 0.64e2 * t902 * t173 + 0.32e2 * t908 * t173;
                            t1300 = 0.32e2 * t816 * t168 - 0.2e1 * t819 * t168 + 0.64e2 * t861 * t168 - 0.64e2 * t864 * t168 + 0.256e3 * t398 * t183 + 0.8e1 * t404 * t183 + 0.32e2 * t422 * t183 - 0.2e1 * t425 * t183 + 0.256e3 * t776 * t183 - 0.64e2 * t779 * t183 + 0.64e2 * t782 * t183 - 0.256e3 * t995 * t183;
                            t1326 = 0.256e3 * t1006 * t168 - 0.32e2 * t352 * t168 + 0.2e1 * t356 * t168 - 0.64e2 * t371 * t183 - 0.32e2 * t401 * t183 + 0.2e1 * t407 * t183 + 0.32e2 * t410 * t183 - 0.2e1 * t413 * t183 - 0.64e2 * t416 * t183 + 0.64e2 * t419 * t183 + 0.64e2 * t754 * t183 - 0.256e3 * t787 * t183;
                            t1349 = -0.32e2 * t843 * t168 + 0.2e1 * t846 * t168 - 0.8e1 * t374 * t183 - 0.64e2 * t902 * t183 + 0.8e1 * t905 * t183 + 0.4e1 * t580 * t562 - t566 * t595 + t593 * t568 - 0.256e3 * t597 * t568 - 0.4e1 * t961 * t576 - 0.256e3 * t958 * t595 + 0.256e3 * t964 * t595;
                            t1354 = t44 * t124;
                            t1355 = t1354 * t38;
                            t1356 = jn(t25, t205);
                            t1357 = t130 * t1356;
                            t1359 = c4 * t209 * t1357;
                            t1362 = t3 * omega;
                            t1363 = t1354 * t1362;
                            t1364 = t171 * t1356;
                            t1366 = c4 * t136 * t1364;
                            t1369 = n * t20;
                            t1370 = t1354 * t1369;
                            t1372 = c4 * t155 * t1364;
                            t1378 = c4 * t184 * t1357;
                            t1381 = n * t93;
                            t1382 = t18 * t124;
                            t1383 = t1382 * t1381;
                            t1386 = t1382 * t60;
                            t1389 = t1382 * t63;
                            t1392 = t1354 * t66;
                            t1395 = t1354 * t69;
                            t1398 = t1354 * t1381;
                            t1401 = -0.256e3 * t1359 * t1355 + 0.256e3 * t1359 * t1363 - 0.256e3 * t1366 * t1363 + 0.18e2 * t1366 * t1383 - 0.256e3 * t1366 * t1392 + 0.160e3 * t1366 * t1395 - 0.18e2 * t1366 * t1398 + 0.64e2 * t1372 * t1370 - 0.64e2 * t1378 * t1370 + 0.64e2 * t1372 * t1386 - 0.32e2 * t1372 * t1389 + 0.4e1 * t961 * t595;
                            t1402 = t1382 * t1369;
                            t1405 = t1382 * t69;
                            t1428 = 0.8e1 * t431 * t132 - 0.8e1 * t795 * t132 + 0.64e2 * t798 * t132 - 0.64e2 * t826 * t132 + 0.2e1 * t829 * t132 - 0.8e1 * t855 * t132 + 0.32e2 * t858 * t132 - 0.18e2 * t1359 * t1383 + 0.32e2 * t1359 * t1405 - 0.64e2 * t1378 * t1386 + 0.64e2 * t1378 * t1402 - 0.2e1 * t151 * t173;
                            t1454 = 0.64e2 * t437 * t132 - 0.64e2 * t440 * t132 - 0.32e2 * t926 * t132 - 0.64e2 * t371 * t173 + 0.8e1 * t374 * t173 + 0.64e2 * t379 * t173 + 0.8e1 * t382 * t173 - 0.8e1 * t404 * t173 - 0.32e2 * t410 * t173 + 0.2e1 * t413 * t173 - 0.32e2 * t751 * t173 + 0.64e2 * t754 * t173;
                            t1479 = -0.32e2 * t140 * t173 + 0.64e2 * t146 * t173 + 0.256e3 * t602 * t568 - 0.256e3 * t943 * t571 - 0.32e2 * t946 * t571 + 0.256e3 * t949 * t571 - 0.32e2 * t952 * t571 + 0.32e2 * t955 * t571 - 0.256e3 * t573 * t725 + 0.32e2 * t583 * t725 - 0.32e2 * t587 * t725 + 0.256e3 * t590 * t725;
                            t1496 = 0.1e1 / n;
                            t1497 = t1496 * t91;
                            t1498 = t1354 * t1497;
                            t1503 = t1382 * t41;
                            t1506 = t1382 * t35;
                            t1509 = -0.32e2 * t1359 * t1503 - 0.4e1 * t1372 * t1498 + 0.4e1 * t1378 * t1498 + 0.64e2 * t1378 * t1506 - t593 * t562 + t566 * t576 - 0.4e1 * t580 * t568 + 0.256e3 * t609 * t568 - 0.256e3 * t934 * t568 - 0.256e3 * t967 * t595 + 0.256e3 * t970 * t595;
                            t1534 = 0.64e2 * t262 * t1053 - 0.40e2 * t266 * t1053 + 0.4e1 * t512 * t1053 - 0.64e2 * t523 * t1053 + 0.40e2 * t526 * t1053 - 0.4e1 * t529 * t1053 - 0.64e2 * t262 * t1116 + 0.40e2 * t266 * t1116 - 0.4e1 * t512 * t1116 - 0.64e2 * t1372 * t1402 - 0.256e3 * t301 * t168 + 0.256e3 * t306 * t168;
                            t1558 = t13 * t20;
                            t1559 = t1354 * t1558;
                            t1562 = 0.64e2 * t257 * t1061 + 0.64e2 * t201 * t1074 + 0.64e2 * t523 * t1116 - 0.40e2 * t526 * t1116 + 0.4e1 * t529 * t1116 + 0.32e2 * t1378 * t1389 + 0.64e2 * t1378 * t1559 - 0.18e2 * t233 * t208 + 0.32e2 * t272 * t208 + 0.18e2 * t233 * t215 - 0.160e3 * t272 * t215 - 0.18e2 * t222 * t271;
                            t1563 = t1354 * t22;
                            t1566 = t21 * t31;
                            t1567 = t1496 * t1566;
                            t1568 = t1354 * t1567;
                            t1570 = t1354 * t28;
                            t1573 = t1354 * t32;
                            t1586 = t1354 * t60;
                            t1589 = t1354 * t63;
                            t1592 = t1354 * t35;
                            t1595 = 0.18e2 * t1359 * t1398 - t1359 * t1568 - 0.256e3 * t1359 * t1570 + 0.32e2 * t1359 * t1573 + 0.256e3 * t1366 * t1570 - 0.32e2 * t1366 * t1573 - 0.64e2 * t1372 * t1559 + 0.8e1 * t1372 * t1563 - 0.8e1 * t1378 * t1563 - 0.64e2 * t1378 * t1586 + 0.32e2 * t1378 * t1589 + 0.64e2 * t1378 * t1592;
                            t1598 = t1382 * t1497;
                            t1605 = t1354 * t41;
                            t1610 = t1382 * t1558;
                            t1613 = t1382 * t22;
                            t1624 = 0.256e3 * t1366 * t1355 - 0.128e3 * t1366 * t1605 + 0.64e2 * t1372 * t1586 - 0.32e2 * t1372 * t1589 + 0.4e1 * t1372 * t1598 + 0.64e2 * t1372 * t1610 - 0.8e1 * t1372 * t1613 - 0.4e1 * t1378 * t1598 - 0.64e2 * t1378 * t1610 + 0.8e1 * t1378 * t1613 - t982 * t975 + t978 * t980;
                            t1637 = t1382 * t1567;
                            t1641 = t1356 * t124;
                            t1642 = t1641 * t63;
                            t1643 = t174 * t130;
                            t1644 = c4 * t161;
                            t1645 = t1644 * t1643;
                            t1648 = t1641 * t66;
                            t1649 = t299 * t1643;
                            t1652 = t1641 * t28;
                            t1655 = t1641 * t32;
                            t1658 = 0.256e3 * t1359 * t1392 - 0.160e3 * t1359 * t1395 + 0.128e3 * t1359 * t1605 - 0.32e2 * t1366 * t1405 + 0.32e2 * t1366 * t1503 - t1366 * t1637 - 0.64e2 * t1372 * t1506 - 0.32e2 * t1645 * t1642 - 0.256e3 * t1649 * t1648 + 0.256e3 * t1649 * t1652 - 0.32e2 * t1649 * t1655 + 0.32e2 * t939 * t725;
                            t1660 = t1641 * t69;
                            t1663 = t134 * t171;
                            t1664 = t299 * t1663;
                            t1669 = t1641 * t60;
                            t1670 = t1644 * t1663;
                            t1681 = t124 * t21;
                            t1682 = t1356 * t18;
                            t1683 = t1682 * t1681;
                            t1684 = t578 * t1643;
                            t1687 = t124 * t93;
                            t1688 = t1682 * t1687;
                            t1689 = t564 * t1663;
                            t1692 = t564 * t1643;
                            t1695 = t578 * t1663;
                            t1698 = 0.32e2 * t1670 * t1642 + 0.256e3 * t1645 * t1669 + 0.256e3 * t1664 * t1648 + 0.32e2 * t1649 * t1660 - 0.256e3 * t1664 * t1652 + 0.32e2 * t1664 * t1655 - 0.32e2 * t1664 * t1660 - 0.256e3 * t1670 * t1669 - 0.24e2 * t1684 * t1683 + 0.24e2 * t1695 * t1683 - 0.14e2 * t1689 * t1688 + 0.14e2 * t1692 * t1688;
                            t1723 = -0.64e2 * t201 * t1061 - 0.8e1 * t520 * t1061 + 0.14e2 * t463 * t208 - 0.64e2 * t466 * t208 - 0.64e2 * t552 * t215 + 0.32e2 * t555 * t215 + 0.160e3 * t218 * t271 + 0.256e3 * t734 * t571 + 0.256e3 * t737 * t571 - 0.256e3 * t721 * t725 + 0.256e3 * t740 * t725 - 0.256e3 * t748 * t725;
                            t1749 = 0.64e2 * t240 * t1053 - 0.64e2 * t289 * t1053 + 0.32e2 * t292 * t1053 - 0.64e2 * t240 * t1074 + 0.32e2 * t244 * t1074 + 0.64e2 * t289 * t1074 - 0.32e2 * t292 * t1074 + 0.8e1 * t469 * t208 - 0.14e2 * t251 * t232 + 0.64e2 * t546 * t232 - 0.8e1 * t549 * t232;
                            t1774 = -0.64e2 * t289 * t1061 - 0.64e2 * t187 * t168 - 0.8e1 * t431 * t168 + 0.64e2 * t437 * t168 - 0.64e2 * t440 * t168 - 0.8e1 * t446 * t168 + 0.256e3 * t450 * t168 + 0.8e1 * t822 * t168 + 0.64e2 * t881 * t168 - 0.2e1 * t343 * t183 + 0.8e1 * t348 * t183 + 0.64e2 * t227 * t271;
                            t1796 = t1356 * t44;
                            t1797 = t1796 * t1687;
                            t1800 = t1359 * t1637 + t1366 * t1568 - 0.64e2 * t1372 * t1592 - 0.32e2 * t146 * t183 + 0.2e1 * t151 * t183 - 0.8e1 * t157 * t183 + 0.2e1 * t1689 * t1797 + 0.64e2 * t379 * t183 - 0.8e1 * t382 * t183 - 0.256e3 * t727 * t571 - 0.256e3 * t731 * t571 + 0.256e3 * t745 * t725;
                            t1801 = t1796 * t1681;
                            t1808 = t1641 * t35;
                            t1811 = t1641 * t38;
                            t1818 = t1641 * t1362;
                            t1829 = 0.64e2 * t262 * t1061 - 0.40e2 * t266 * t1061 + 0.4e1 * t512 * t1061 - 0.256e3 * t1645 * t1808 + 0.256e3 * t1649 * t1811 - 0.256e3 * t1649 * t1818 - 0.256e3 * t1664 * t1811 + 0.256e3 * t1664 * t1818 + 0.256e3 * t1670 * t1808 - 0.24e2 * t1684 * t1801 - 0.2e1 * t1692 * t1797 + 0.24e2 * t1695 * t1801;
                            t1856 = 0.8e1 * t194 * t1061 + 0.64e2 * t517 * t1061 - 0.64e2 * t523 * t1061 + 0.40e2 * t526 * t1061 - 0.4e1 * t529 * t1061 - 0.64e2 * t532 * t1061 + 0.32e2 * t460 * t215 - 0.2e1 * t463 * t215 + 0.64e2 * t466 * t215 - 0.8e1 * t469 * t215 - 0.32e2 * t247 * t271 + 0.2e1 * t251 * t271;
                            t1881 = 0.256e3 * t320 * t168 - 0.256e3 * t324 * t168 - 0.256e3 * t484 * t183 + 0.256e3 * t487 * t183 - 0.256e3 * t680 * t183 + 0.32e2 * t210 * t232 + 0.128e3 * t617 * t215 + 0.64e2 * t620 * t215 - 0.24e2 * t623 * t215 - 0.64e2 * t281 * t232 + 0.24e2 * t284 * t232 - 0.64e2 * t546 * t271;
                            t1907 = 0.8e1 * t520 * t1053 - 0.8e1 * t194 * t1074 - 0.64e2 * t517 * t1074 + 0.8e1 * t520 * t1074 + 0.64e2 * t532 * t1074 + 0.8e1 * t194 * t1116 + 0.64e2 * t517 * t1116 - 0.8e1 * t520 * t1116 - 0.64e2 * t532 * t1116 + 0.256e3 * t475 * t183 + 0.256e3 * t683 * t183 + 0.8e1 * t549 * t271;
                            t1929 = yn(t25, t205);
                            t1930 = t1929 * t1496 * omega;
                            t1935 = -c4 * t1643 * t1930 + c4 * t1663 * t1930 + 0.64e2 * t201 * t1053 - 0.64e2 * t257 * t1053 - 0.64e2 * t201 * t1116 + 0.64e2 * t257 * t1116 - 0.256e3 * t490 * t183 + 0.64e2 * t457 * t208 - 0.32e2 * t617 * t208 + 0.64e2 * t620 * t208 - 0.24e2 * t623 * t208 - 0.64e2 * t543 * t232;
                            vr = t1935 + t1749 + t1173 + t1349 + t1249 + t1698 + t1198 + t1774 + t1479 + t1881 + t1146 + t1658 + t1723 + t1224 + t1907 + t1829 + t1800 + t1624 + t1401 + t1562 + t1428 + t1121 + t1300 + t1275 + t1509 + t1534 + t1093 + t1454 + t1066 + t1856 + t1595 + t1326;
                            t1940 = t198 * t1354;
                            t1941 = t161 * t26;
                            t1942 = t150 * t1941;
                            t1945 = t176 * t26;
                            t1946 = t139 * t1945;
                            t1949 = t98 * t135;
                            t1950 = t139 * t1949;
                            t1953 = t2 * t135;
                            t1954 = t145 * t1953;
                            t1957 = t98 * t154;
                            t1958 = t150 * t1957;
                            t1961 = t2 * t154;
                            t1962 = t139 * t1961;
                            t1965 = t237 * t1382;
                            t1966 = t176 * t17;
                            t1967 = t150 * t1966;
                            t1970 = t255 * t1382;
                            t1971 = t355 * t1961;
                            t1974 = t21 * t144;
                            t1975 = t1974 * t1966;
                            t1978 = t1974 * t1941;
                            t1981 = t31 * t216;
                            t1982 = t1981 * t1945;
                            t1985 = t1974 * t1953;
                            t1988 = t1974 * t1957;
                            t1991 = 0.64e2 * t1942 * t1940 - 0.64e2 * t1946 * t1940 + 0.256e3 * t1950 * t1940 - 0.256e3 * t1954 * t1940 - 0.64e2 * t1958 * t1940 + 0.64e2 * t1962 * t1940 - 0.32e2 * t1967 * t1965 - 0.32e2 * t1975 * t1965 - 0.64e2 * t1978 * t1965 + 0.64e2 * t1982 * t1965 + 0.32e2 * t1985 * t1965 + 0.64e2 * t1988 * t1965 + 0.8e1 * t1971 * t1970;
                            t1992 = t1981 * t1961;
                            t1995 = t217 * t1966;
                            t1998 = t217 * t1953;
                            t2001 = t191 * t1354;
                            t2006 = t221 * t1966;
                            t2009 = t355 * t1945;
                            t2012 = t221 * t1953;
                            t2017 = t161 * t17;
                            t2018 = t355 * t2017;
                            t2021 = t126 * t124;
                            t2022 = t134 * t130;
                            t2023 = t2022 * t2021;
                            t2024 = t93 * t138;
                            t2025 = t2024 * t193;
                            t2030 = t254 * t1382;
                            t2031 = t93 * t304;
                            t2032 = t2031 * t155;
                            t2035 = t1566 * t1644;
                            t2036 = t2035 * t160;
                            t2039 = -0.8e1 * t1971 * t1940 + 0.256e3 * t1995 * t1940 - 0.256e3 * t1998 * t1940 - 0.32e2 * t2006 * t1940 + 0.8e1 * t2009 * t1940 + 0.32e2 * t2012 * t1940 - 0.64e2 * t1992 * t1965 - 0.32e2 * t2018 * t1965 + 0.256e3 * t1995 * t2001 - 0.256e3 * t1998 * t2001 + 0.40e2 * t2025 * t2023 - 0.64e2 * t266 * t2023 - 0.8e1 * t2032 * t2030 - 0.2e1 * t2036 * t2030;
                            t2041 = c4 * t98;
                            t2042 = t1566 * t2041;
                            t2043 = t2042 * t136;
                            t2046 = t127 * t124;
                            t2047 = t2022 * t2046;
                            t2056 = t197 * t1354;
                            t2057 = t2035 * t298;
                            t2060 = t2042 * t209;
                            t2063 = t236 * t1382;
                            t2068 = t91 * t299;
                            t2069 = t2068 * t298;
                            t2072 = t91 * t1644;
                            t2073 = t2072 * t175;
                            t2076 = t91 * t304;
                            t2077 = t2076 * t209;
                            t2080 = t174 * t171;
                            t2081 = t2080 * t2021;
                            t2084 = 0.2e1 * t2043 * t2030 - 0.64e2 * t244 * t2047 + 0.64e2 * t292 * t2047 + 0.8e1 * t512 * t2047 - 0.8e1 * t529 * t2047 + 0.2e1 * t2057 * t2056 - 0.2e1 * t2060 * t2056 + 0.2e1 * t2069 * t2056 - 0.8e1 * t2073 * t2056 - 0.2e1 * t2077 * t2056 - 0.2e1 * t2057 * t2063 + 0.2e1 * t2060 * t2063 + 0.64e2 * t292 * t2081;
                            t2087 = t2080 * t2046;
                            t2106 = t150 * t1953;
                            t2109 = t2024 * t2017;
                            t2114 = t91 * t149;
                            t2115 = t2114 * t1966;
                            t2118 = 0.64e2 * t194 * t2081 - 0.64e2 * t194 * t2087 + 0.32e2 * t1975 * t1940 - 0.32e2 * t2109 * t1940 - 0.2e1 * t2115 * t1940 + 0.64e2 * t1942 * t1970 - 0.64e2 * t1946 * t1970 + 0.32e2 * t2106 * t1970 - 0.64e2 * t244 * t2081 + 0.8e1 * t512 * t2081 - 0.64e2 * t520 * t2081 + 0.8e1 * t512 * t2087 + 0.64e2 * t520 * t2087 - 0.8e1 * t529 * t2087;
                            t2131 = t186 * t1966;
                            t2134 = t243 * t1945;
                            t2137 = t186 * t1953;
                            t2140 = t243 * t1961;
                            t2147 = t2024 * t1949;
                            t2152 = -0.32e2 * t1985 * t1940 - 0.256e3 * t2131 * t1940 + 0.64e2 * t2134 * t1940 + 0.256e3 * t2137 * t1940 - 0.64e2 * t2140 * t1940 + 0.32e2 * t2147 * t1940 + 0.64e2 * t1942 * t1965 - 0.64e2 * t1946 * t1965 - 0.64e2 * t1958 * t1965 + 0.64e2 * t1962 * t1965 + 0.32e2 * t2106 * t1965 + 0.64e2 * t2134 * t1965 - 0.64e2 * t2140 * t1965;
                            t2153 = t2114 * t1953;
                            t2156 = t265 * t1945;
                            t2159 = t265 * t1961;
                            t2166 = t355 * t1949;
                            t2185 = 0.2e1 * t2153 * t1940 - 0.8e1 * t2009 * t1965 + 0.64e2 * t2159 * t1965 + 0.32e2 * t2166 * t1965 - 0.32e2 * t1975 * t1970 + 0.32e2 * t2109 * t1970 - 0.8e1 * t1971 * t2001 - 0.32e2 * t1985 * t2001 + 0.8e1 * t2009 * t2001 + 0.32e2 * t2012 * t2001 + 0.32e2 * t2147 * t2001 + 0.2e1 * t2153 * t2001 - 0.64e2 * t2156 * t2001 + 0.64e2 * t2159 * t2001;
                            t2203 = t139 * t2017;
                            t2206 = t145 * t1966;
                            t2215 = 0.64e2 * t1942 * t2001 - 0.64e2 * t1946 * t2001 + 0.256e3 * t1950 * t2001 + 0.64e2 * t2006 * t1965 - 0.64e2 * t2156 * t1965 + 0.2e1 * t2115 * t1970 - 0.32e2 * t2147 * t1970 - 0.2e1 * t2153 * t1970 + 0.32e2 * t1975 * t2001 - 0.32e2 * t2109 * t2001 - 0.2e1 * t2115 * t2001 - 0.256e3 * t2203 * t2001 + 0.256e3 * t2206 * t2001;
                            t2244 = -0.256e3 * t1954 * t2001 - 0.64e2 * t1958 * t2001 + 0.64e2 * t1962 * t2001 + 0.32e2 * t2109 * t1965 + 0.2e1 * t2115 * t1965 - 0.32e2 * t2147 * t1965 - 0.32e2 * t1967 * t1970 + 0.64e2 * t2006 * t1970 - 0.32e2 * t2018 * t1970 - 0.32e2 * t2006 * t2001 - 0.256e3 * t2131 * t2001 + 0.64e2 * t2134 * t2001 + 0.256e3 * t2137 * t2001 - 0.64e2 * t2140 * t2001;
                            t2262 = t1981 * t1949;
                            t2265 = t20 * t225;
                            t2266 = t2265 * t1953;
                            t2277 = 0.64e2 * t1988 * t1940 - 0.64e2 * t1992 * t1940 - 0.64e2 * t2156 * t1940 + 0.64e2 * t2159 * t1940 - 0.2e1 * t2153 * t1965 - 0.64e2 * t1978 * t1970 + 0.64e2 * t1982 * t1970 + 0.32e2 * t1985 * t1970 + 0.64e2 * t1982 * t2001 + 0.64e2 * t1988 * t2001 - 0.64e2 * t1992 * t2001 - 0.256e3 * t2262 * t2001 + 0.256e3 * t2266 * t2001;
                            t2306 = -0.256e3 * t2203 * t1940 + 0.256e3 * t2206 * t1940 - 0.64e2 * t1958 * t1970 + 0.64e2 * t1962 * t1970 + 0.8e1 * t1971 * t1965 - 0.64e2 * t2012 * t1965 + 0.64e2 * t1988 * t1970 - 0.64e2 * t1992 * t1970 - 0.8e1 * t2009 * t1970 - 0.64e2 * t2012 * t1970 + 0.64e2 * t2134 * t1970 - 0.64e2 * t2140 * t1970 - 0.64e2 * t2156 * t1970 + 0.32e2 * t2166 * t1970;
                            t2310 = t1981 * t2017;
                            t2313 = t2265 * t1966;
                            t2336 = -0.64e2 * t1978 * t1940 + 0.64e2 * t1982 * t1940 - 0.256e3 * t2262 * t1940 + 0.256e3 * t2266 * t1940 + 0.256e3 * t2310 * t1940 - 0.256e3 * t2313 * t1940 + 0.64e2 * t2159 * t1970 - 0.64e2 * t1978 * t2001 + 0.256e3 * t2310 * t2001 - 0.256e3 * t2313 * t2001 - 0.64e2 * t266 * t2047 + 0.64e2 * t526 * t2047 - 0.64e2 * t266 * t2087;
                            t2337 = t1981 * t193;
                            t2342 = t190 * t1354;
                            t2343 = t93 * t299;
                            t2344 = t2343 * t346;
                            t2357 = t1981 * t200;
                            t2360 = t2024 * t200;
                            t2371 = 0.64e2 * t194 * t2023 - 0.64e2 * t2337 * t2023 - 0.40e2 * t2025 * t2047 - 0.40e2 * t2025 * t2087 + 0.8e1 * t2344 * t2030 + 0.8e1 * t2032 * t2342 + 0.64e2 * t2337 * t2047 - 0.64e2 * t2357 * t2047 + 0.40e2 * t2360 * t2047 - 0.64e2 * t2337 * t2081 + 0.64e2 * t2337 * t2087 - 0.64e2 * t2357 * t2087 + 0.40e2 * t2360 * t2087 - 0.8e1 * t2344 * t2342;
                            t2384 = t2068 * t160;
                            t2401 = 0.64e2 * t2357 * t2023 - 0.40e2 * t2360 * t2023 + 0.8e1 * t512 * t2023 - 0.64e2 * t520 * t2023 + 0.64e2 * t526 * t2023 - 0.8e1 * t529 * t2023 + 0.40e2 * t2025 * t2081 + 0.64e2 * t2357 * t2081 - 0.40e2 * t2360 * t2081 - 0.64e2 * t266 * t2081 + 0.64e2 * t526 * t2081 + 0.64e2 * t526 * t2087 + 0.2e1 * t2384 * t2342;
                            t2402 = t93 * c4;
                            t2403 = t2402 * t193;
                            t2406 = t2402 * t200;
                            t2409 = t1566 * c4;
                            t2410 = t2409 * t193;
                            t2413 = t2409 * t200;
                            t2436 = -0.32e2 * t2403 * t2023 + 0.32e2 * t2406 * t2023 - 0.4e1 * t2410 * t2023 + 0.4e1 * t2413 * t2023 + 0.32e2 * t2403 * t2047 - 0.32e2 * t2406 * t2047 + 0.4e1 * t2410 * t2047 - 0.4e1 * t2413 * t2047 - 0.32e2 * t2403 * t2081 + 0.32e2 * t2406 * t2081 - 0.4e1 * t2410 * t2081 + 0.4e1 * t2413 * t2081 + 0.32e2 * t2403 * t2087 - 0.32e2 * t2406 * t2087;
                            t2448 = t2072 * t346;
                            t2451 = t2076 * t136;
                            t2454 = t91 * t2041;
                            t2455 = t2454 * t155;
                            t2468 = 0.64e2 * t292 * t2023 - 0.2e1 * t2384 * t2030 + 0.8e1 * t2448 * t2030 + 0.2e1 * t2451 * t2030 - 0.8e1 * t2455 * t2030 - 0.8e1 * t529 * t2081 + 0.4e1 * t2410 * t2087 - 0.4e1 * t2413 * t2087 - 0.64e2 * t244 * t2087 + 0.64e2 * t292 * t2087 - 0.8e1 * t2448 * t2342 - 0.2e1 * t2451 * t2342 + 0.8e1 * t2455 * t2342;
                            t2479 = t2454 * t184;
                            t2490 = t2343 * t175;
                            t2493 = t2031 * t184;
                            t2500 = -0.64e2 * t194 * t2047 - 0.64e2 * t244 * t2023 + 0.2e1 * t2036 * t2342 - 0.2e1 * t2043 * t2342 + 0.64e2 * t520 * t2047 + 0.8e1 * t2479 * t2056 - 0.8e1 * t2490 * t2056 + 0.8e1 * t2493 * t2056 - 0.2e1 * t2069 * t2063 + 0.8e1 * t2073 * t2063 + 0.2e1 * t2077 * t2063 - 0.8e1 * t2479 * t2063 + 0.8e1 * t2490 * t2063 - 0.8e1 * t2493 * t2063;
                            pr = t1991 + t2039 + t2084 + t2118 + t2152 + t2185 + t2215 + t2244 + t2277 + t2306 + t2336 + t2371 + t2401 + t2436 + t2468 + t2500;
              // (u1,u2) = ur*rHat + vr*thetaHat 
                            u(i1,i2,i3,u1c) = ur*cosTheta - vr*sinTheta; 
                            u(i1,i2,i3,u2c) = ur*sinTheta + vr*cosTheta;
                            u(i1,i2,i3,pc ) = pr;
                            u1Max=max(u1Max,u(i1,i2,i3,u1c));
                            u2Max=max(u2Max,u(i1,i2,i3,u2c));
                            pMax =max(pMax, u(i1,i2,i3,pc));
                        }
                        uNorm = max(u1Max,u2Max,u3Max,pMax);
            // if( scale==0. )
            // {
            //   printF("UDKS: ERROR: annulus solution is zero! scale=%e\n",scale);
            //   OV_ABORT("ERROR");
            // }
            // scale = 1./scale;
            // c4 *= scale;  // set this so the velocity and acceleration are scaled too
            // FOR_3(i1,i2,i3,I1,I2,I3)
            // {
            //   u(i1,i2,i3,u1c) *= scale;
            //   u(i1,i2,i3,u2c) *= scale;
            //   u(i1,i2,i3,pc ) *= scale;
            // }    
                    }
                    else
                    {
                        OV_ABORT("UDKS:ERROR: annulus solution : 3D called");
                    }
                if( uNorm==0. )
                {
                    printF("UDKS: ERROR: annulus solution is zero! uNorm=%e\n",uNorm);
                    OV_ABORT("ERROR");
                }        
                annulusScaling = c4/uNorm;
            }
            c4 = annulusScaling; 

                Real ur,vr,pr; 
                Real u1Max=0., u2Max=0., u3Max=0., pMax=0.; 
                uNorm=0; 
                if( numberOfDimensions==2 )
                {
                    FOR_3D(i1,i2,i3,I1,I2,I3)
                    {
                        Real x = vertex(i1,i2,i3,0);
                        Real y = vertex(i1,i2,i3,1);
                        Real r = sqrt( x*x + y*y );
                        Real theta = atan2(y,x); 
                        Real cosTheta = x/r;
                        Real sinTheta = y/r;
                        t1 = 0.1e1 / r;
                        t2 = jn(n, omega);
                        t3 = n * n;
                        t4 = t3 * t3;
                        t5 = t4 * t3;
                        t8 = t4 * t2;
                        t10 = t4 * n;
                        t13 = t3 * n;
                        t14 = t13 * t2;
                        t16 = omega / 0.2e1;
                        t17 = jn(n, t16);
                        t18 = pow(0.1e1 / 0.2e1, n);
                        t19 = t18 * t17;
                        t20 = omega * omega;
                        t21 = t20 * t20;
                        t22 = t21 * n;
                        t25 = n + 0.1e1;
                        t26 = jn(t25, t16);
                        t27 = t18 * t26;
                        t28 = omega * t4;
                        t31 = t20 * omega;
                        t32 = t31 * t3;
                        t35 = t20 * t3;
                        t38 = omega * t13;
                        t41 = t31 * n;
                        t44 = 0.1e1 / t18;
                        t45 = t17 * t44;
                        t48 = t26 * t44;
                        t53 = t3 * t26;
                        t60 = t20 * t4;
                        t63 = t21 * t3;
                        t66 = omega * t10;
                        t69 = t31 * t13;
                        t74 = 0.64e2 * t18 * omega * t53 - 0.64e2 * t44 * omega * t53 + 0.256e3 * t10 * t2 + 0.14e2 * t22 * t19 - 0.32e2 * t35 * t19 - 0.256e3 * t5 * t2 - 0.2e1 * t22 * t45 - 0.64e2 * t28 * t27 + 0.8e1 * t32 * t27 + 0.64e2 * t38 * t27 - 0.24e2 * t41 * t27 + 0.128e3 * t35 * t45 + 0.64e2 * t38 * t48 - 0.24e2 * t41 * t48 - 0.160e3 * t60 * t45 + 0.18e2 * t63 * t45 - 0.64e2 * t66 * t48 + 0.32e2 * t69 * t48 - 0.256e3 * t14 + 0.256e3 * t8;
                        t91 = t21 * t20;
                        t93 = t21 * omega;
                        t98 = jn(t25, omega);
                        t99 = t13 * t98;
                        t122 = 0.32e2 * t44 * t20 * t13 * t17 + 0.256e3 * omega * t10 * t98 - 0.256e3 * t44 * t10 * t17 + 0.256e3 * t44 * t5 * t17 - 0.256e3 * omega * t99 + 0.256e3 * t13 * t45 - 0.32e2 * t20 * t14 + 0.32e2 * t60 * t19 - 0.18e2 * t63 * t19 + t91 * t19 + 0.32e2 * t20 * t8 - 0.64e2 * t66 * t27 + 0.32e2 * t69 * t27 - 0.4e1 * t93 * t27 + 0.64e2 * t28 * t48 - 0.32e2 * t31 * t99 - 0.8e1 * t32 * t48 - 0.256e3 * t4 * t45 - t91 * t45 + 0.4e1 * t93 * t48;
                        t124 = 0.1e1 / (t74 + t122);
                        t125 = t124 * t1;
                        t126 = pow(r, n);
                        t127 = 0.1e1 / t126;
                        t128 = t127 * t18;
                        t129 = omega * t;
                        t130 = sin(t129);
                        t132 = t130 * t128 * t125;
                        t133 = n * theta;
                        t134 = sin(t133);
                        t135 = yn(n, t16);
                        t136 = t135 * t134;
                        t137 = t98 * t136;
                        t138 = t3 * c4;
                        t139 = t31 * t138;
                        t140 = t139 * t137;
                        t143 = t2 * t136;
                        t144 = t13 * c4;
                        t145 = t20 * t144;
                        t146 = t145 * t143;
                        t149 = n * c4;
                        t150 = t21 * t149;
                        t151 = t150 * t143;
                        t154 = yn(t25, t16);
                        t155 = t154 * t134;
                        t156 = t98 * t155;
                        t157 = t150 * t156;
                        t160 = t17 * t134;
                        t161 = yn(t25, omega);
                        t162 = t161 * t160;
                        t163 = t139 * t162;
                        t166 = t126 * t44;
                        t168 = t130 * t166 * t125;
                        t171 = cos(t129);
                        t173 = t171 * t128 * t125;
                        t174 = cos(t133);
                        t175 = t26 * t174;
                        t176 = yn(n, omega);
                        t177 = t176 * t175;
                        t178 = omega * t138;
                        t179 = t178 * t177;
                        t183 = t171 * t166 * t125;
                        t184 = t154 * t174;
                        t185 = t98 * t184;
                        t186 = t20 * t138;
                        t187 = t186 * t185;
                        t190 = t130 * t126;
                        t191 = t134 * t190;
                        t192 = t191 * t125;
                        t193 = t154 * t17;
                        t194 = t139 * t193;
                        t197 = t171 * t126;
                        t198 = t174 * t197;
                        t199 = t198 * t125;
                        t200 = t135 * t26;
                        t201 = t178 * t200;
                        t204 = t171 * t18;
                        t205 = omega * r;
                        t206 = jn(n, t205);
                        t208 = t206 * t204 * t125;
                        t209 = t135 * t174;
                        t210 = t186 * t209;
                        t213 = t171 * t44;
                        t215 = t206 * t213 * t125;
                        t216 = t4 * c4;
                        t217 = t20 * t216;
                        t218 = t217 * t209;
                        t221 = t21 * t138;
                        t222 = t221 * t209;
                        t225 = t10 * c4;
                        t226 = omega * t225;
                        t227 = t226 * t184;
                        t230 = t130 * t18;
                        t232 = t206 * t230 * t125;
                        t233 = t221 * t136;
                        t236 = t171 * t127;
                        t237 = t174 * t236;
                        t238 = t237 * t125;
                        t239 = omega * t144;
                        t240 = t239 * t200;
                        t243 = t31 * t149;
                        t244 = t243 * t200;
                        t247 = t145 * t209;
                        t250 = -0.32e2 * t140 * t132 + 0.64e2 * t146 * t132 - 0.2e1 * t151 * t132 + 0.8e1 * t157 * t132 + 0.32e2 * t163 * t132 - 0.2e1 * t151 * t168 - 0.64e2 * t179 * t173 - 0.64e2 * t187 * t183 + 0.8e1 * t194 * t192 - 0.64e2 * t201 * t199 + 0.32e2 * t210 * t208 + 0.18e2 * t222 * t208 + 0.64e2 * t227 * t208 + 0.160e3 * t218 * t215 - 0.32e2 * t247 * t215 + 0.18e2 * t233 * t232 - 0.64e2 * t240 * t238 + 0.32e2 * t244 * t238;
                        t251 = t150 * t209;
                        t254 = t130 * t127;
                        t255 = t134 * t254;
                        t256 = t255 * t125;
                        t257 = t178 * t193;
                        t262 = t226 * t193;
                        t265 = t31 * t144;
                        t266 = t265 * t193;
                        t269 = t130 * t44;
                        t271 = t206 * t269 * t125;
                        t272 = t217 * t136;
                        t281 = t239 * t184;
                        t284 = t243 * t184;
                        t289 = t239 * t193;
                        t292 = t243 * t193;
                        t295 = t265 * t184;
                        t298 = t17 * t174;
                        t299 = c4 * t176;
                        t300 = t5 * t299;
                        t301 = t300 * t298;
                        t304 = c4 * t2;
                        t305 = t5 * t304;
                        t306 = t305 * t209;
                        t309 = t10 * t299;
                        t310 = t309 * t298;
                        t313 = t10 * t304;
                        t314 = t313 * t209;
                        t317 = -0.256e3 * t301 * t183 + 0.256e3 * t306 * t183 + 0.256e3 * t310 * t183 - 0.256e3 * t314 * t183 + 0.64e2 * t201 * t238 + 0.64e2 * t201 * t256 - 0.32e2 * t218 * t208 - 0.32e2 * t295 * t208 - 0.128e3 * t210 * t215 + 0.2e1 * t251 * t215 - 0.64e2 * t281 * t215 + 0.24e2 * t284 * t215 - 0.18e2 * t233 * t271 + 0.64e2 * t289 * t238 - 0.32e2 * t292 * t238 - 0.64e2 * t257 * t256 - 0.64e2 * t262 * t256 + 0.40e2 * t266 * t256 + 0.160e3 * t272 * t271;
                        t319 = t4 * t299;
                        t320 = t319 * t298;
                        t323 = t4 * t304;
                        t324 = t323 * t209;
                        t327 = t13 * t299;
                        t328 = t327 * t298;
                        t331 = t13 * t304;
                        t332 = t331 * t209;
                        t339 = t176 * t160;
                        t340 = t145 * t339;
                        t343 = t150 * t339;
                        t346 = t26 * t134;
                        t347 = t161 * t346;
                        t348 = t150 * t347;
                        t351 = t161 * t298;
                        t352 = t265 * t351;
                        t355 = t93 * t149;
                        t356 = t355 * t351;
                        t359 = t176 * t298;
                        t360 = t217 * t359;
                        t363 = t221 * t359;
                        t366 = t161 * t175;
                        t367 = t217 * t366;
                        t370 = t2 * t155;
                        t371 = t239 * t370;
                        t374 = t243 * t370;
                        t377 = t176 * t346;
                        t378 = omega * t216;
                        t379 = t378 * t377;
                        t382 = t139 * t377;
                        t385 = -0.64e2 * t340 * t132 + 0.2e1 * t343 * t132 - 0.8e1 * t348 * t132 + 0.64e2 * t379 * t132 + 0.8e1 * t382 * t132 + 0.64e2 * t371 * t168 + 0.8e1 * t374 * t168 - 0.32e2 * t352 * t173 + 0.2e1 * t356 * t173 + 0.32e2 * t360 * t173 - 0.2e1 * t363 * t173 + 0.64e2 * t367 * t173 + 0.256e3 * t320 * t183 - 0.256e3 * t324 * t183 - 0.256e3 * t328 * t183 + 0.256e3 * t332 * t183 + 0.64e2 * t257 * t199 - 0.18e2 * t222 * t215;
                        t386 = t217 * t339;
                        t389 = t221 * t339;
                        t392 = t217 * t347;
                        t395 = t226 * t377;
                        t398 = t226 * t137;
                        t401 = t265 * t137;
                        t404 = t243 * t377;
                        t407 = t355 * t137;
                        t410 = t217 * t143;
                        t413 = t221 * t143;
                        t416 = t217 * t156;
                        t419 = t226 * t370;
                        t422 = t265 * t162;
                        t425 = t355 * t162;
                        t430 = t2 * t184;
                        t431 = t139 * t430;
                        t434 = t239 * t351;
                        t437 = t186 * t366;
                        t440 = t239 * t177;
                        t443 = 0.32e2 * t386 * t132 - 0.8e1 * t404 * t132 - 0.32e2 * t422 * t132 + 0.2e1 * t425 * t132 + 0.32e2 * t386 * t168 - 0.2e1 * t389 * t168 - 0.64e2 * t392 * t168 + 0.64e2 * t395 * t168 - 0.256e3 * t398 * t168 + 0.32e2 * t401 * t168 - 0.2e1 * t407 * t168 - 0.32e2 * t410 * t168 + 0.2e1 * t413 * t168 + 0.64e2 * t416 * t168 - 0.64e2 * t419 * t168 - 0.8e1 * t431 * t173 - 0.256e3 * t434 * t183 + 0.64e2 * t437 * t183 - 0.64e2 * t440 * t183;
                        t446 = t243 * t177;
                        t449 = t98 * t209;
                        t450 = t239 * t449;
                        t457 = t178 * t155;
                        t460 = t145 * t136;
                        t463 = t150 * t136;
                        t466 = t378 * t155;
                        t469 = t139 * t155;
                        t472 = t331 * t136;
                        t475 = t300 * t160;
                        t484 = t319 * t160;
                        t487 = t323 * t136;
                        t490 = t305 * t136;
                        t495 = 0.256e3 * t472 * t168 - 0.256e3 * t475 * t168 + 0.256e3 * t484 * t168 - 0.256e3 * t487 * t168 + 0.256e3 * t490 * t168 - 0.8e1 * t446 * t183 + 0.256e3 * t450 * t183 - 0.64e2 * t281 * t208 - 0.14e2 * t463 * t232 + 0.64e2 * t466 * t232 - 0.8e1 * t469 * t232 - 0.64e2 * t240 * t256 + 0.32e2 * t244 * t256 + 0.64e2 * t457 * t271 - 0.32e2 * t460 * t271 + 0.2e1 * t463 * t271 - 0.64e2 * t466 * t271 + 0.8e1 * t469 * t271;
                        t512 = t355 * t193;
                        t517 = t378 * t200;
                        t520 = t139 * t200;
                        t523 = t226 * t200;
                        t526 = t265 * t200;
                        t529 = t355 * t200;
                        t532 = t378 * t193;
                        t541 = 0.8e1 * t194 * t199 - 0.8e1 * t194 * t238 - 0.64e2 * t240 * t199 + 0.32e2 * t244 * t199 + 0.64e2 * t289 * t199 - 0.32e2 * t292 * t199 + 0.64e2 * t517 * t199 - 0.8e1 * t520 * t199 - 0.64e2 * t532 * t199 + 0.24e2 * t284 * t208 - 0.32e2 * t295 * t215 - 0.64e2 * t262 * t238 + 0.40e2 * t266 * t238 - 0.4e1 * t512 * t238 - 0.64e2 * t517 * t238 + 0.8e1 * t520 * t238 + 0.64e2 * t523 * t238 - 0.40e2 * t526 * t238 + 0.4e1 * t529 * t238;
                        t543 = t178 * t184;
                        t546 = t378 * t184;
                        t549 = t139 * t184;
                        t552 = t226 * t155;
                        t555 = t265 * t155;
                        t562 = t230 * t125;
                        t563 = t134 * t206;
                        t564 = c4 * t135;
                        t565 = t91 * t564;
                        t566 = t565 * t563;
                        t568 = t269 * t125;
                        t571 = t206 * t171 * t125;
                        t572 = t161 * t174;
                        t573 = t226 * t572;
                        t576 = t204 * t125;
                        t577 = t174 * t206;
                        t578 = c4 * t154;
                        t579 = t93 * t578;
                        t580 = t579 * t577;
                        t583 = t265 * t572;
                        t586 = t176 * t174;
                        t587 = t217 * t586;
                        t590 = t239 * t572;
                        t593 = t565 * t577;
                        t595 = t213 * t125;
                        t596 = t5 * t564;
                        t597 = t596 * t577;
                        t601 = t10 * t564;
                        t602 = t601 * t577;
                        t607 = -0.64e2 * t543 * t208 - 0.64e2 * t546 * t215 + 0.8e1 * t549 * t215 + 0.64e2 * t552 * t232 - 0.32e2 * t555 * t232 + 0.64e2 * t289 * t256 - 0.32e2 * t292 * t256 - t566 * t562 + t566 * t568 - 0.256e3 * t573 * t571 + 0.32e2 * t583 * t571 - 0.32e2 * t587 * t571 + 0.256e3 * t590 * t571 + 0.4e1 * t580 * t576 - t593 * t576 - 0.4e1 * t580 * t595 + t593 * t595 - 0.256e3 * t597 * t595 + 0.256e3 * t602 * t595;
                        t608 = t4 * t564;
                        t609 = t608 * t577;
                        t612 = t327 * t160;
                        t617 = t186 * t136;
                        t620 = t239 * t155;
                        t623 = t243 * t155;
                        t652 = -0.256e3 * t612 * t168 + 0.40e2 * t266 * t192 - 0.8e1 * t194 * t256 - 0.4e1 * t512 * t199 + 0.64e2 * t523 * t199 - 0.40e2 * t526 * t199 + 0.4e1 * t529 * t199 - 0.32e2 * t272 * t232 - 0.4e1 * t512 * t256 - 0.64e2 * t517 * t256 + 0.8e1 * t520 * t256 + 0.64e2 * t523 * t256 - 0.40e2 * t526 * t256 + 0.4e1 * t529 * t256 + 0.64e2 * t532 * t256 - 0.128e3 * t617 * t271 - 0.64e2 * t620 * t271 + 0.24e2 * t623 * t271 + 0.256e3 * t609 * t595;
                        t680 = t309 * t160;
                        t683 = t313 * t136;
                        t694 = 0.256e3 * t680 * t168 - 0.256e3 * t683 * t168 - 0.64e2 * t240 * t192 + 0.64e2 * t289 * t192 - 0.32e2 * t292 * t192 + 0.64e2 * t517 * t192 - 0.8e1 * t520 * t192 - 0.40e2 * t526 * t192 + 0.4e1 * t529 * t192 - 0.64e2 * t532 * t192 + 0.64e2 * t227 * t215 + 0.64e2 * t543 * t215 + 0.32e2 * t617 * t232 - 0.64e2 * t620 * t232 + 0.24e2 * t623 * t232 - 0.64e2 * t257 * t238 + 0.64e2 * t552 * t271 - 0.32e2 * t555 * t271;
                        t721 = t225 * t586;
                        t725 = t206 * t130 * t125;
                        t726 = t176 * t134;
                        t727 = t144 * t726;
                        t730 = t5 * c4;
                        t731 = t730 * t726;
                        t734 = t225 * t726;
                        t737 = t216 * t726;
                        t740 = t730 * t586;
                        t743 = -0.64e2 * t201 * t192 + 0.32e2 * t244 * t192 + 0.64e2 * t257 * t192 - 0.64e2 * t262 * t192 - 0.4e1 * t512 * t192 + 0.64e2 * t523 * t192 - 0.64e2 * t262 * t199 + 0.40e2 * t266 * t199 - 0.14e2 * t251 * t208 + 0.64e2 * t546 * t208 - 0.8e1 * t549 * t208 - 0.64e2 * t457 * t232 + 0.64e2 * t532 * t238 - 0.256e3 * t721 * t571 + 0.256e3 * t740 * t571 + 0.256e3 * t727 * t725 + 0.256e3 * t731 * t725 - 0.256e3 * t734 * t725 - 0.256e3 * t737 * t725;
                        t745 = t144 * t586;
                        t748 = t216 * t586;
                        t751 = t186 * t143;
                        t754 = t186 * t156;
                        t773 = t178 * t377;
                        t776 = t239 * t162;
                        t779 = t186 * t347;
                        t782 = t239 * t377;
                        t787 = t239 * t137;
                        t790 = -0.64e2 * t371 * t132 + 0.8e1 * t374 * t132 - 0.32e2 * t751 * t132 + 0.64e2 * t754 * t132 + 0.8e1 * t157 * t168 - 0.32e2 * t340 * t168 + 0.2e1 * t343 * t168 - 0.8e1 * t348 * t168 - 0.64e2 * t379 * t168 + 0.8e1 * t382 * t168 - 0.8e1 * t404 * t168 + 0.64e2 * t773 * t168 - 0.256e3 * t776 * t168 + 0.64e2 * t779 * t168 - 0.64e2 * t782 * t168 + 0.256e3 * t787 * t168 + 0.256e3 * t745 * t571 - 0.256e3 * t748 * t571;
                        t795 = t150 * t185;
                        t798 = t378 * t430;
                        t813 = t226 * t177;
                        t816 = t265 * t449;
                        t819 = t355 * t449;
                        t822 = t243 * t430;
                        t825 = t2 * t209;
                        t826 = t145 * t825;
                        t829 = t150 * t825;
                        t834 = t186 * t825;
                        t837 = t139 * t351;
                        t840 = -0.2e1 * t389 * t132 + 0.64e2 * t392 * t132 - 0.64e2 * t395 * t132 + 0.32e2 * t401 * t132 - 0.2e1 * t407 * t132 - 0.32e2 * t410 * t132 + 0.32e2 * t146 * t168 - 0.64e2 * t754 * t168 - 0.8e1 * t446 * t173 + 0.8e1 * t795 * t173 - 0.64e2 * t798 * t173 - 0.64e2 * t813 * t173 + 0.32e2 * t816 * t173 - 0.2e1 * t819 * t173 + 0.8e1 * t822 * t173 + 0.64e2 * t826 * t173 - 0.2e1 * t829 * t173 - 0.32e2 * t834 * t173 + 0.32e2 * t837 * t173;
                        t843 = t145 * t359;
                        t846 = t150 * t359;
                        t849 = t150 * t366;
                        t852 = t378 * t177;
                        t855 = t139 * t177;
                        t858 = t139 * t449;
                        t861 = t217 * t185;
                        t864 = t226 * t430;
                        t881 = t239 * t430;
                        t886 = t178 * t430;
                        t889 = -0.64e2 * t843 * t173 + 0.2e1 * t846 * t173 - 0.8e1 * t849 * t173 + 0.64e2 * t852 * t173 + 0.8e1 * t855 * t173 - 0.32e2 * t858 * t173 - 0.64e2 * t861 * t173 + 0.64e2 * t864 * t173 + 0.64e2 * t886 * t173 + 0.8e1 * t822 * t183 + 0.32e2 * t826 * t183 - 0.2e1 * t829 * t183 - 0.32e2 * t843 * t183 + 0.2e1 * t846 * t183 - 0.8e1 * t849 * t183 - 0.64e2 * t852 * t183 + 0.8e1 * t855 * t183 + 0.64e2 * t881 * t183;
                        t896 = t217 * t825;
                        t899 = t221 * t825;
                        t902 = t378 * t370;
                        t905 = t139 * t370;
                        t908 = t186 * t339;
                        t917 = t178 * t370;
                        t926 = t186 * t359;
                        t933 = t13 * t564;
                        t934 = t933 * t577;
                        t937 = 0.2e1 * t413 * t132 - 0.64e2 * t416 * t132 + 0.64e2 * t419 * t132 - 0.64e2 * t773 * t132 - 0.64e2 * t779 * t132 + 0.64e2 * t782 * t132 - 0.64e2 * t902 * t132 - 0.8e1 * t905 * t132 + 0.32e2 * t908 * t132 + 0.64e2 * t917 * t132 + 0.64e2 * t902 * t168 - 0.8e1 * t905 * t168 - 0.64e2 * t917 * t168 - 0.64e2 * t437 * t173 + 0.64e2 * t440 * t173 - 0.32e2 * t896 * t173 + 0.2e1 * t899 * t173 + 0.32e2 * t926 * t173 - 0.256e3 * t934 * t595;
                        t939 = t145 * t586;
                        t942 = t161 * t134;
                        t943 = t239 * t942;
                        t946 = t145 * t726;
                        t949 = t226 * t942;
                        t952 = t265 * t942;
                        t955 = t217 * t726;
                        t958 = t601 * t563;
                        t961 = t579 * t563;
                        t964 = t596 * t563;
                        t967 = t608 * t563;
                        t970 = t933 * t563;
                        t975 = t171 * t1;
                        t976 = yn(n, t205);
                        t978 = c4 * t174 * t976;
                        t980 = t130 * t1;
                        t982 = c4 * t134 * t976;
                        t994 = 0.64e2 * t187 * t173 - 0.64e2 * t881 * t173 - 0.8e1 * t431 * t183 + 0.8e1 * t795 * t183 + 0.64e2 * t798 * t183 + 0.4e1 * t961 * t562 + 0.256e3 * t958 * t568 - 0.4e1 * t961 * t568 - 0.256e3 * t964 * t568 + 0.256e3 * t967 * t568 - 0.256e3 * t970 * t568 + 0.32e2 * t939 * t571 + 0.256e3 * t943 * t725 + 0.32e2 * t946 * t725 - 0.256e3 * t949 * t725 + 0.32e2 * t952 * t725 - 0.32e2 * t955 * t725 + t978 * t975 + t982 * t980;
                        t995 = t226 * t162;
                        t1006 = t226 * t351;
                        t1021 = t226 * t449;
                        t1036 = 0.256e3 * t1006 * t183 - 0.256e3 * t1021 * t183 - 0.32e2 * t422 * t168 + 0.2e1 * t425 * t168 + 0.256e3 * t995 * t168 + 0.64e2 * t179 * t183 - 0.32e2 * t352 * t183 + 0.2e1 * t356 * t183 + 0.32e2 * t360 * t183 - 0.2e1 * t363 * t183 - 0.64e2 * t367 * t183 + 0.64e2 * t813 * t183 + 0.32e2 * t816 * t183 - 0.2e1 * t819 * t183 + 0.64e2 * t861 * t183 - 0.64e2 * t864 * t183 - 0.64e2 * t886 * t183 - 0.32e2 * t896 * t183 + 0.2e1 * t899 * t183;
                        ur = t250 + t317 + t385 + t443 + t495 + t541 + t607 + t652 + t694 + t743 + t790 + t840 + t889 + t937 + t994 + t1036;
                        t1053 = t134 * t197 * t125;
                        t1061 = t174 * t254 * t125;
                        t1066 = -0.8e1 * t194 * t1053 - 0.64e2 * t517 * t1053 + 0.64e2 * t532 * t1053 + 0.64e2 * t240 * t1061 + 0.32e2 * t292 * t1061 - 0.64e2 * t552 * t208 + 0.32e2 * t555 * t208 - 0.32e2 * t218 * t232 + 0.18e2 * t222 * t232 + 0.64e2 * t227 * t232 - 0.32e2 * t295 * t232;
                        t1074 = t134 * t236 * t125;
                        t1093 = -0.32e2 * t244 * t1061 - 0.64e2 * t257 * t1074 - 0.256e3 * t314 * t168 - 0.256e3 * t328 * t168 - 0.256e3 * t434 * t168 - 0.256e3 * t472 * t183 + 0.256e3 * t612 * t183 - 0.128e3 * t210 * t271 - 0.64e2 * t457 * t215 - 0.64e2 * t281 * t271 + 0.24e2 * t284 * t271 + 0.64e2 * t543 * t271;
                        t1116 = t174 * t190 * t125;
                        t1121 = -0.32e2 * t244 * t1053 - 0.64e2 * t262 * t1074 + 0.40e2 * t266 * t1074 - 0.4e1 * t512 * t1074 + 0.64e2 * t523 * t1074 - 0.40e2 * t526 * t1074 + 0.4e1 * t529 * t1074 + 0.64e2 * t289 * t1116 - 0.32e2 * t292 * t1116 + 0.256e3 * t310 * t168 + 0.256e3 * t332 * t168 - 0.32e2 * t295 * t271;
                        t1146 = -0.64e2 * t240 * t1116 + 0.32e2 * t244 * t1116 - 0.32e2 * t837 * t132 + 0.64e2 * t843 * t132 - 0.2e1 * t846 * t132 + 0.8e1 * t849 * t132 - 0.64e2 * t852 * t132 + 0.32e2 * t163 * t173 - 0.64e2 * t340 * t173 + 0.2e1 * t343 * t173 - 0.8e1 * t348 * t173 - 0.8e1 * t905 * t173;
                        t1173 = 0.64e2 * t179 * t132 - 0.64e2 * t886 * t132 + 0.64e2 * t179 * t168 - 0.8e1 * t849 * t168 - 0.64e2 * t852 * t168 - 0.64e2 * t886 * t168 - 0.32e2 * t896 * t168 + 0.2e1 * t899 * t168 - 0.64e2 * t416 * t173 - 0.64e2 * t773 * t173 + 0.64e2 * t917 * t173 + 0.64e2 * t917 * t183;
                        t1198 = 0.32e2 * t352 * t132 - 0.2e1 * t356 * t132 - 0.32e2 * t360 * t132 + 0.2e1 * t363 * t132 - 0.64e2 * t367 * t132 + 0.64e2 * t813 * t132 - 0.32e2 * t816 * t132 + 0.32e2 * t386 * t173 - 0.2e1 * t389 * t173 + 0.64e2 * t392 * t173 - 0.64e2 * t395 * t173 + 0.64e2 * t419 * t173;
                        t1224 = -0.64e2 * t187 * t132 + 0.8e1 * t446 * t132 - 0.8e1 * t822 * t132 + 0.32e2 * t834 * t132 + 0.64e2 * t881 * t132 + 0.32e2 * t360 * t168 - 0.2e1 * t363 * t168 - 0.64e2 * t367 * t168 + 0.32e2 * t401 * t173 - 0.2e1 * t407 * t173 + 0.32e2 * t340 * t183 - 0.64e2 * t773 * t183;
                        t1249 = -0.256e3 * t1021 * t168 + 0.2e1 * t819 * t132 + 0.64e2 * t861 * t132 - 0.64e2 * t864 * t132 + 0.32e2 * t896 * t132 - 0.2e1 * t899 * t132 + 0.64e2 * t813 * t168 + 0.8e1 * t855 * t168 - 0.32e2 * t386 * t183 + 0.2e1 * t389 * t183 + 0.64e2 * t392 * t183 - 0.64e2 * t395 * t183;
                        t1275 = 0.8e1 * t157 * t173 + 0.8e1 * t795 * t168 + 0.64e2 * t798 * t168 + 0.32e2 * t826 * t168 - 0.2e1 * t829 * t168 - 0.32e2 * t422 * t173 + 0.2e1 * t425 * t173 - 0.64e2 * t779 * t173 + 0.64e2 * t782 * t173 - 0.64e2 * t902 * t173 + 0.32e2 * t908 * t173;
                        t1300 = 0.32e2 * t816 * t168 - 0.2e1 * t819 * t168 + 0.64e2 * t861 * t168 - 0.64e2 * t864 * t168 + 0.256e3 * t398 * t183 + 0.8e1 * t404 * t183 + 0.32e2 * t422 * t183 - 0.2e1 * t425 * t183 + 0.256e3 * t776 * t183 - 0.64e2 * t779 * t183 + 0.64e2 * t782 * t183 - 0.256e3 * t995 * t183;
                        t1326 = 0.256e3 * t1006 * t168 - 0.32e2 * t352 * t168 + 0.2e1 * t356 * t168 - 0.64e2 * t371 * t183 - 0.32e2 * t401 * t183 + 0.2e1 * t407 * t183 + 0.32e2 * t410 * t183 - 0.2e1 * t413 * t183 - 0.64e2 * t416 * t183 + 0.64e2 * t419 * t183 + 0.64e2 * t754 * t183 - 0.256e3 * t787 * t183;
                        t1349 = -0.32e2 * t843 * t168 + 0.2e1 * t846 * t168 - 0.8e1 * t374 * t183 - 0.64e2 * t902 * t183 + 0.8e1 * t905 * t183 + 0.4e1 * t580 * t562 - t566 * t595 + t593 * t568 - 0.256e3 * t597 * t568 - 0.4e1 * t961 * t576 - 0.256e3 * t958 * t595 + 0.256e3 * t964 * t595;
                        t1354 = t44 * t124;
                        t1355 = t1354 * t38;
                        t1356 = jn(t25, t205);
                        t1357 = t130 * t1356;
                        t1359 = c4 * t209 * t1357;
                        t1362 = t3 * omega;
                        t1363 = t1354 * t1362;
                        t1364 = t171 * t1356;
                        t1366 = c4 * t136 * t1364;
                        t1369 = n * t20;
                        t1370 = t1354 * t1369;
                        t1372 = c4 * t155 * t1364;
                        t1378 = c4 * t184 * t1357;
                        t1381 = n * t93;
                        t1382 = t18 * t124;
                        t1383 = t1382 * t1381;
                        t1386 = t1382 * t60;
                        t1389 = t1382 * t63;
                        t1392 = t1354 * t66;
                        t1395 = t1354 * t69;
                        t1398 = t1354 * t1381;
                        t1401 = -0.256e3 * t1359 * t1355 + 0.256e3 * t1359 * t1363 - 0.256e3 * t1366 * t1363 + 0.18e2 * t1366 * t1383 - 0.256e3 * t1366 * t1392 + 0.160e3 * t1366 * t1395 - 0.18e2 * t1366 * t1398 + 0.64e2 * t1372 * t1370 - 0.64e2 * t1378 * t1370 + 0.64e2 * t1372 * t1386 - 0.32e2 * t1372 * t1389 + 0.4e1 * t961 * t595;
                        t1402 = t1382 * t1369;
                        t1405 = t1382 * t69;
                        t1428 = 0.8e1 * t431 * t132 - 0.8e1 * t795 * t132 + 0.64e2 * t798 * t132 - 0.64e2 * t826 * t132 + 0.2e1 * t829 * t132 - 0.8e1 * t855 * t132 + 0.32e2 * t858 * t132 - 0.18e2 * t1359 * t1383 + 0.32e2 * t1359 * t1405 - 0.64e2 * t1378 * t1386 + 0.64e2 * t1378 * t1402 - 0.2e1 * t151 * t173;
                        t1454 = 0.64e2 * t437 * t132 - 0.64e2 * t440 * t132 - 0.32e2 * t926 * t132 - 0.64e2 * t371 * t173 + 0.8e1 * t374 * t173 + 0.64e2 * t379 * t173 + 0.8e1 * t382 * t173 - 0.8e1 * t404 * t173 - 0.32e2 * t410 * t173 + 0.2e1 * t413 * t173 - 0.32e2 * t751 * t173 + 0.64e2 * t754 * t173;
                        t1479 = -0.32e2 * t140 * t173 + 0.64e2 * t146 * t173 + 0.256e3 * t602 * t568 - 0.256e3 * t943 * t571 - 0.32e2 * t946 * t571 + 0.256e3 * t949 * t571 - 0.32e2 * t952 * t571 + 0.32e2 * t955 * t571 - 0.256e3 * t573 * t725 + 0.32e2 * t583 * t725 - 0.32e2 * t587 * t725 + 0.256e3 * t590 * t725;
                        t1496 = 0.1e1 / n;
                        t1497 = t1496 * t91;
                        t1498 = t1354 * t1497;
                        t1503 = t1382 * t41;
                        t1506 = t1382 * t35;
                        t1509 = -0.32e2 * t1359 * t1503 - 0.4e1 * t1372 * t1498 + 0.4e1 * t1378 * t1498 + 0.64e2 * t1378 * t1506 - t593 * t562 + t566 * t576 - 0.4e1 * t580 * t568 + 0.256e3 * t609 * t568 - 0.256e3 * t934 * t568 - 0.256e3 * t967 * t595 + 0.256e3 * t970 * t595;
                        t1534 = 0.64e2 * t262 * t1053 - 0.40e2 * t266 * t1053 + 0.4e1 * t512 * t1053 - 0.64e2 * t523 * t1053 + 0.40e2 * t526 * t1053 - 0.4e1 * t529 * t1053 - 0.64e2 * t262 * t1116 + 0.40e2 * t266 * t1116 - 0.4e1 * t512 * t1116 - 0.64e2 * t1372 * t1402 - 0.256e3 * t301 * t168 + 0.256e3 * t306 * t168;
                        t1558 = t13 * t20;
                        t1559 = t1354 * t1558;
                        t1562 = 0.64e2 * t257 * t1061 + 0.64e2 * t201 * t1074 + 0.64e2 * t523 * t1116 - 0.40e2 * t526 * t1116 + 0.4e1 * t529 * t1116 + 0.32e2 * t1378 * t1389 + 0.64e2 * t1378 * t1559 - 0.18e2 * t233 * t208 + 0.32e2 * t272 * t208 + 0.18e2 * t233 * t215 - 0.160e3 * t272 * t215 - 0.18e2 * t222 * t271;
                        t1563 = t1354 * t22;
                        t1566 = t21 * t31;
                        t1567 = t1496 * t1566;
                        t1568 = t1354 * t1567;
                        t1570 = t1354 * t28;
                        t1573 = t1354 * t32;
                        t1586 = t1354 * t60;
                        t1589 = t1354 * t63;
                        t1592 = t1354 * t35;
                        t1595 = 0.18e2 * t1359 * t1398 - t1359 * t1568 - 0.256e3 * t1359 * t1570 + 0.32e2 * t1359 * t1573 + 0.256e3 * t1366 * t1570 - 0.32e2 * t1366 * t1573 - 0.64e2 * t1372 * t1559 + 0.8e1 * t1372 * t1563 - 0.8e1 * t1378 * t1563 - 0.64e2 * t1378 * t1586 + 0.32e2 * t1378 * t1589 + 0.64e2 * t1378 * t1592;
                        t1598 = t1382 * t1497;
                        t1605 = t1354 * t41;
                        t1610 = t1382 * t1558;
                        t1613 = t1382 * t22;
                        t1624 = 0.256e3 * t1366 * t1355 - 0.128e3 * t1366 * t1605 + 0.64e2 * t1372 * t1586 - 0.32e2 * t1372 * t1589 + 0.4e1 * t1372 * t1598 + 0.64e2 * t1372 * t1610 - 0.8e1 * t1372 * t1613 - 0.4e1 * t1378 * t1598 - 0.64e2 * t1378 * t1610 + 0.8e1 * t1378 * t1613 - t982 * t975 + t978 * t980;
                        t1637 = t1382 * t1567;
                        t1641 = t1356 * t124;
                        t1642 = t1641 * t63;
                        t1643 = t174 * t130;
                        t1644 = c4 * t161;
                        t1645 = t1644 * t1643;
                        t1648 = t1641 * t66;
                        t1649 = t299 * t1643;
                        t1652 = t1641 * t28;
                        t1655 = t1641 * t32;
                        t1658 = 0.256e3 * t1359 * t1392 - 0.160e3 * t1359 * t1395 + 0.128e3 * t1359 * t1605 - 0.32e2 * t1366 * t1405 + 0.32e2 * t1366 * t1503 - t1366 * t1637 - 0.64e2 * t1372 * t1506 - 0.32e2 * t1645 * t1642 - 0.256e3 * t1649 * t1648 + 0.256e3 * t1649 * t1652 - 0.32e2 * t1649 * t1655 + 0.32e2 * t939 * t725;
                        t1660 = t1641 * t69;
                        t1663 = t134 * t171;
                        t1664 = t299 * t1663;
                        t1669 = t1641 * t60;
                        t1670 = t1644 * t1663;
                        t1681 = t124 * t21;
                        t1682 = t1356 * t18;
                        t1683 = t1682 * t1681;
                        t1684 = t578 * t1643;
                        t1687 = t124 * t93;
                        t1688 = t1682 * t1687;
                        t1689 = t564 * t1663;
                        t1692 = t564 * t1643;
                        t1695 = t578 * t1663;
                        t1698 = 0.32e2 * t1670 * t1642 + 0.256e3 * t1645 * t1669 + 0.256e3 * t1664 * t1648 + 0.32e2 * t1649 * t1660 - 0.256e3 * t1664 * t1652 + 0.32e2 * t1664 * t1655 - 0.32e2 * t1664 * t1660 - 0.256e3 * t1670 * t1669 - 0.24e2 * t1684 * t1683 + 0.24e2 * t1695 * t1683 - 0.14e2 * t1689 * t1688 + 0.14e2 * t1692 * t1688;
                        t1723 = -0.64e2 * t201 * t1061 - 0.8e1 * t520 * t1061 + 0.14e2 * t463 * t208 - 0.64e2 * t466 * t208 - 0.64e2 * t552 * t215 + 0.32e2 * t555 * t215 + 0.160e3 * t218 * t271 + 0.256e3 * t734 * t571 + 0.256e3 * t737 * t571 - 0.256e3 * t721 * t725 + 0.256e3 * t740 * t725 - 0.256e3 * t748 * t725;
                        t1749 = 0.64e2 * t240 * t1053 - 0.64e2 * t289 * t1053 + 0.32e2 * t292 * t1053 - 0.64e2 * t240 * t1074 + 0.32e2 * t244 * t1074 + 0.64e2 * t289 * t1074 - 0.32e2 * t292 * t1074 + 0.8e1 * t469 * t208 - 0.14e2 * t251 * t232 + 0.64e2 * t546 * t232 - 0.8e1 * t549 * t232;
                        t1774 = -0.64e2 * t289 * t1061 - 0.64e2 * t187 * t168 - 0.8e1 * t431 * t168 + 0.64e2 * t437 * t168 - 0.64e2 * t440 * t168 - 0.8e1 * t446 * t168 + 0.256e3 * t450 * t168 + 0.8e1 * t822 * t168 + 0.64e2 * t881 * t168 - 0.2e1 * t343 * t183 + 0.8e1 * t348 * t183 + 0.64e2 * t227 * t271;
                        t1796 = t1356 * t44;
                        t1797 = t1796 * t1687;
                        t1800 = t1359 * t1637 + t1366 * t1568 - 0.64e2 * t1372 * t1592 - 0.32e2 * t146 * t183 + 0.2e1 * t151 * t183 - 0.8e1 * t157 * t183 + 0.2e1 * t1689 * t1797 + 0.64e2 * t379 * t183 - 0.8e1 * t382 * t183 - 0.256e3 * t727 * t571 - 0.256e3 * t731 * t571 + 0.256e3 * t745 * t725;
                        t1801 = t1796 * t1681;
                        t1808 = t1641 * t35;
                        t1811 = t1641 * t38;
                        t1818 = t1641 * t1362;
                        t1829 = 0.64e2 * t262 * t1061 - 0.40e2 * t266 * t1061 + 0.4e1 * t512 * t1061 - 0.256e3 * t1645 * t1808 + 0.256e3 * t1649 * t1811 - 0.256e3 * t1649 * t1818 - 0.256e3 * t1664 * t1811 + 0.256e3 * t1664 * t1818 + 0.256e3 * t1670 * t1808 - 0.24e2 * t1684 * t1801 - 0.2e1 * t1692 * t1797 + 0.24e2 * t1695 * t1801;
                        t1856 = 0.8e1 * t194 * t1061 + 0.64e2 * t517 * t1061 - 0.64e2 * t523 * t1061 + 0.40e2 * t526 * t1061 - 0.4e1 * t529 * t1061 - 0.64e2 * t532 * t1061 + 0.32e2 * t460 * t215 - 0.2e1 * t463 * t215 + 0.64e2 * t466 * t215 - 0.8e1 * t469 * t215 - 0.32e2 * t247 * t271 + 0.2e1 * t251 * t271;
                        t1881 = 0.256e3 * t320 * t168 - 0.256e3 * t324 * t168 - 0.256e3 * t484 * t183 + 0.256e3 * t487 * t183 - 0.256e3 * t680 * t183 + 0.32e2 * t210 * t232 + 0.128e3 * t617 * t215 + 0.64e2 * t620 * t215 - 0.24e2 * t623 * t215 - 0.64e2 * t281 * t232 + 0.24e2 * t284 * t232 - 0.64e2 * t546 * t271;
                        t1907 = 0.8e1 * t520 * t1053 - 0.8e1 * t194 * t1074 - 0.64e2 * t517 * t1074 + 0.8e1 * t520 * t1074 + 0.64e2 * t532 * t1074 + 0.8e1 * t194 * t1116 + 0.64e2 * t517 * t1116 - 0.8e1 * t520 * t1116 - 0.64e2 * t532 * t1116 + 0.256e3 * t475 * t183 + 0.256e3 * t683 * t183 + 0.8e1 * t549 * t271;
                        t1929 = yn(t25, t205);
                        t1930 = t1929 * t1496 * omega;
                        t1935 = -c4 * t1643 * t1930 + c4 * t1663 * t1930 + 0.64e2 * t201 * t1053 - 0.64e2 * t257 * t1053 - 0.64e2 * t201 * t1116 + 0.64e2 * t257 * t1116 - 0.256e3 * t490 * t183 + 0.64e2 * t457 * t208 - 0.32e2 * t617 * t208 + 0.64e2 * t620 * t208 - 0.24e2 * t623 * t208 - 0.64e2 * t543 * t232;
                        vr = t1935 + t1749 + t1173 + t1349 + t1249 + t1698 + t1198 + t1774 + t1479 + t1881 + t1146 + t1658 + t1723 + t1224 + t1907 + t1829 + t1800 + t1624 + t1401 + t1562 + t1428 + t1121 + t1300 + t1275 + t1509 + t1534 + t1093 + t1454 + t1066 + t1856 + t1595 + t1326;
                        t1940 = t198 * t1354;
                        t1941 = t161 * t26;
                        t1942 = t150 * t1941;
                        t1945 = t176 * t26;
                        t1946 = t139 * t1945;
                        t1949 = t98 * t135;
                        t1950 = t139 * t1949;
                        t1953 = t2 * t135;
                        t1954 = t145 * t1953;
                        t1957 = t98 * t154;
                        t1958 = t150 * t1957;
                        t1961 = t2 * t154;
                        t1962 = t139 * t1961;
                        t1965 = t237 * t1382;
                        t1966 = t176 * t17;
                        t1967 = t150 * t1966;
                        t1970 = t255 * t1382;
                        t1971 = t355 * t1961;
                        t1974 = t21 * t144;
                        t1975 = t1974 * t1966;
                        t1978 = t1974 * t1941;
                        t1981 = t31 * t216;
                        t1982 = t1981 * t1945;
                        t1985 = t1974 * t1953;
                        t1988 = t1974 * t1957;
                        t1991 = 0.64e2 * t1942 * t1940 - 0.64e2 * t1946 * t1940 + 0.256e3 * t1950 * t1940 - 0.256e3 * t1954 * t1940 - 0.64e2 * t1958 * t1940 + 0.64e2 * t1962 * t1940 - 0.32e2 * t1967 * t1965 - 0.32e2 * t1975 * t1965 - 0.64e2 * t1978 * t1965 + 0.64e2 * t1982 * t1965 + 0.32e2 * t1985 * t1965 + 0.64e2 * t1988 * t1965 + 0.8e1 * t1971 * t1970;
                        t1992 = t1981 * t1961;
                        t1995 = t217 * t1966;
                        t1998 = t217 * t1953;
                        t2001 = t191 * t1354;
                        t2006 = t221 * t1966;
                        t2009 = t355 * t1945;
                        t2012 = t221 * t1953;
                        t2017 = t161 * t17;
                        t2018 = t355 * t2017;
                        t2021 = t126 * t124;
                        t2022 = t134 * t130;
                        t2023 = t2022 * t2021;
                        t2024 = t93 * t138;
                        t2025 = t2024 * t193;
                        t2030 = t254 * t1382;
                        t2031 = t93 * t304;
                        t2032 = t2031 * t155;
                        t2035 = t1566 * t1644;
                        t2036 = t2035 * t160;
                        t2039 = -0.8e1 * t1971 * t1940 + 0.256e3 * t1995 * t1940 - 0.256e3 * t1998 * t1940 - 0.32e2 * t2006 * t1940 + 0.8e1 * t2009 * t1940 + 0.32e2 * t2012 * t1940 - 0.64e2 * t1992 * t1965 - 0.32e2 * t2018 * t1965 + 0.256e3 * t1995 * t2001 - 0.256e3 * t1998 * t2001 + 0.40e2 * t2025 * t2023 - 0.64e2 * t266 * t2023 - 0.8e1 * t2032 * t2030 - 0.2e1 * t2036 * t2030;
                        t2041 = c4 * t98;
                        t2042 = t1566 * t2041;
                        t2043 = t2042 * t136;
                        t2046 = t127 * t124;
                        t2047 = t2022 * t2046;
                        t2056 = t197 * t1354;
                        t2057 = t2035 * t298;
                        t2060 = t2042 * t209;
                        t2063 = t236 * t1382;
                        t2068 = t91 * t299;
                        t2069 = t2068 * t298;
                        t2072 = t91 * t1644;
                        t2073 = t2072 * t175;
                        t2076 = t91 * t304;
                        t2077 = t2076 * t209;
                        t2080 = t174 * t171;
                        t2081 = t2080 * t2021;
                        t2084 = 0.2e1 * t2043 * t2030 - 0.64e2 * t244 * t2047 + 0.64e2 * t292 * t2047 + 0.8e1 * t512 * t2047 - 0.8e1 * t529 * t2047 + 0.2e1 * t2057 * t2056 - 0.2e1 * t2060 * t2056 + 0.2e1 * t2069 * t2056 - 0.8e1 * t2073 * t2056 - 0.2e1 * t2077 * t2056 - 0.2e1 * t2057 * t2063 + 0.2e1 * t2060 * t2063 + 0.64e2 * t292 * t2081;
                        t2087 = t2080 * t2046;
                        t2106 = t150 * t1953;
                        t2109 = t2024 * t2017;
                        t2114 = t91 * t149;
                        t2115 = t2114 * t1966;
                        t2118 = 0.64e2 * t194 * t2081 - 0.64e2 * t194 * t2087 + 0.32e2 * t1975 * t1940 - 0.32e2 * t2109 * t1940 - 0.2e1 * t2115 * t1940 + 0.64e2 * t1942 * t1970 - 0.64e2 * t1946 * t1970 + 0.32e2 * t2106 * t1970 - 0.64e2 * t244 * t2081 + 0.8e1 * t512 * t2081 - 0.64e2 * t520 * t2081 + 0.8e1 * t512 * t2087 + 0.64e2 * t520 * t2087 - 0.8e1 * t529 * t2087;
                        t2131 = t186 * t1966;
                        t2134 = t243 * t1945;
                        t2137 = t186 * t1953;
                        t2140 = t243 * t1961;
                        t2147 = t2024 * t1949;
                        t2152 = -0.32e2 * t1985 * t1940 - 0.256e3 * t2131 * t1940 + 0.64e2 * t2134 * t1940 + 0.256e3 * t2137 * t1940 - 0.64e2 * t2140 * t1940 + 0.32e2 * t2147 * t1940 + 0.64e2 * t1942 * t1965 - 0.64e2 * t1946 * t1965 - 0.64e2 * t1958 * t1965 + 0.64e2 * t1962 * t1965 + 0.32e2 * t2106 * t1965 + 0.64e2 * t2134 * t1965 - 0.64e2 * t2140 * t1965;
                        t2153 = t2114 * t1953;
                        t2156 = t265 * t1945;
                        t2159 = t265 * t1961;
                        t2166 = t355 * t1949;
                        t2185 = 0.2e1 * t2153 * t1940 - 0.8e1 * t2009 * t1965 + 0.64e2 * t2159 * t1965 + 0.32e2 * t2166 * t1965 - 0.32e2 * t1975 * t1970 + 0.32e2 * t2109 * t1970 - 0.8e1 * t1971 * t2001 - 0.32e2 * t1985 * t2001 + 0.8e1 * t2009 * t2001 + 0.32e2 * t2012 * t2001 + 0.32e2 * t2147 * t2001 + 0.2e1 * t2153 * t2001 - 0.64e2 * t2156 * t2001 + 0.64e2 * t2159 * t2001;
                        t2203 = t139 * t2017;
                        t2206 = t145 * t1966;
                        t2215 = 0.64e2 * t1942 * t2001 - 0.64e2 * t1946 * t2001 + 0.256e3 * t1950 * t2001 + 0.64e2 * t2006 * t1965 - 0.64e2 * t2156 * t1965 + 0.2e1 * t2115 * t1970 - 0.32e2 * t2147 * t1970 - 0.2e1 * t2153 * t1970 + 0.32e2 * t1975 * t2001 - 0.32e2 * t2109 * t2001 - 0.2e1 * t2115 * t2001 - 0.256e3 * t2203 * t2001 + 0.256e3 * t2206 * t2001;
                        t2244 = -0.256e3 * t1954 * t2001 - 0.64e2 * t1958 * t2001 + 0.64e2 * t1962 * t2001 + 0.32e2 * t2109 * t1965 + 0.2e1 * t2115 * t1965 - 0.32e2 * t2147 * t1965 - 0.32e2 * t1967 * t1970 + 0.64e2 * t2006 * t1970 - 0.32e2 * t2018 * t1970 - 0.32e2 * t2006 * t2001 - 0.256e3 * t2131 * t2001 + 0.64e2 * t2134 * t2001 + 0.256e3 * t2137 * t2001 - 0.64e2 * t2140 * t2001;
                        t2262 = t1981 * t1949;
                        t2265 = t20 * t225;
                        t2266 = t2265 * t1953;
                        t2277 = 0.64e2 * t1988 * t1940 - 0.64e2 * t1992 * t1940 - 0.64e2 * t2156 * t1940 + 0.64e2 * t2159 * t1940 - 0.2e1 * t2153 * t1965 - 0.64e2 * t1978 * t1970 + 0.64e2 * t1982 * t1970 + 0.32e2 * t1985 * t1970 + 0.64e2 * t1982 * t2001 + 0.64e2 * t1988 * t2001 - 0.64e2 * t1992 * t2001 - 0.256e3 * t2262 * t2001 + 0.256e3 * t2266 * t2001;
                        t2306 = -0.256e3 * t2203 * t1940 + 0.256e3 * t2206 * t1940 - 0.64e2 * t1958 * t1970 + 0.64e2 * t1962 * t1970 + 0.8e1 * t1971 * t1965 - 0.64e2 * t2012 * t1965 + 0.64e2 * t1988 * t1970 - 0.64e2 * t1992 * t1970 - 0.8e1 * t2009 * t1970 - 0.64e2 * t2012 * t1970 + 0.64e2 * t2134 * t1970 - 0.64e2 * t2140 * t1970 - 0.64e2 * t2156 * t1970 + 0.32e2 * t2166 * t1970;
                        t2310 = t1981 * t2017;
                        t2313 = t2265 * t1966;
                        t2336 = -0.64e2 * t1978 * t1940 + 0.64e2 * t1982 * t1940 - 0.256e3 * t2262 * t1940 + 0.256e3 * t2266 * t1940 + 0.256e3 * t2310 * t1940 - 0.256e3 * t2313 * t1940 + 0.64e2 * t2159 * t1970 - 0.64e2 * t1978 * t2001 + 0.256e3 * t2310 * t2001 - 0.256e3 * t2313 * t2001 - 0.64e2 * t266 * t2047 + 0.64e2 * t526 * t2047 - 0.64e2 * t266 * t2087;
                        t2337 = t1981 * t193;
                        t2342 = t190 * t1354;
                        t2343 = t93 * t299;
                        t2344 = t2343 * t346;
                        t2357 = t1981 * t200;
                        t2360 = t2024 * t200;
                        t2371 = 0.64e2 * t194 * t2023 - 0.64e2 * t2337 * t2023 - 0.40e2 * t2025 * t2047 - 0.40e2 * t2025 * t2087 + 0.8e1 * t2344 * t2030 + 0.8e1 * t2032 * t2342 + 0.64e2 * t2337 * t2047 - 0.64e2 * t2357 * t2047 + 0.40e2 * t2360 * t2047 - 0.64e2 * t2337 * t2081 + 0.64e2 * t2337 * t2087 - 0.64e2 * t2357 * t2087 + 0.40e2 * t2360 * t2087 - 0.8e1 * t2344 * t2342;
                        t2384 = t2068 * t160;
                        t2401 = 0.64e2 * t2357 * t2023 - 0.40e2 * t2360 * t2023 + 0.8e1 * t512 * t2023 - 0.64e2 * t520 * t2023 + 0.64e2 * t526 * t2023 - 0.8e1 * t529 * t2023 + 0.40e2 * t2025 * t2081 + 0.64e2 * t2357 * t2081 - 0.40e2 * t2360 * t2081 - 0.64e2 * t266 * t2081 + 0.64e2 * t526 * t2081 + 0.64e2 * t526 * t2087 + 0.2e1 * t2384 * t2342;
                        t2402 = t93 * c4;
                        t2403 = t2402 * t193;
                        t2406 = t2402 * t200;
                        t2409 = t1566 * c4;
                        t2410 = t2409 * t193;
                        t2413 = t2409 * t200;
                        t2436 = -0.32e2 * t2403 * t2023 + 0.32e2 * t2406 * t2023 - 0.4e1 * t2410 * t2023 + 0.4e1 * t2413 * t2023 + 0.32e2 * t2403 * t2047 - 0.32e2 * t2406 * t2047 + 0.4e1 * t2410 * t2047 - 0.4e1 * t2413 * t2047 - 0.32e2 * t2403 * t2081 + 0.32e2 * t2406 * t2081 - 0.4e1 * t2410 * t2081 + 0.4e1 * t2413 * t2081 + 0.32e2 * t2403 * t2087 - 0.32e2 * t2406 * t2087;
                        t2448 = t2072 * t346;
                        t2451 = t2076 * t136;
                        t2454 = t91 * t2041;
                        t2455 = t2454 * t155;
                        t2468 = 0.64e2 * t292 * t2023 - 0.2e1 * t2384 * t2030 + 0.8e1 * t2448 * t2030 + 0.2e1 * t2451 * t2030 - 0.8e1 * t2455 * t2030 - 0.8e1 * t529 * t2081 + 0.4e1 * t2410 * t2087 - 0.4e1 * t2413 * t2087 - 0.64e2 * t244 * t2087 + 0.64e2 * t292 * t2087 - 0.8e1 * t2448 * t2342 - 0.2e1 * t2451 * t2342 + 0.8e1 * t2455 * t2342;
                        t2479 = t2454 * t184;
                        t2490 = t2343 * t175;
                        t2493 = t2031 * t184;
                        t2500 = -0.64e2 * t194 * t2047 - 0.64e2 * t244 * t2023 + 0.2e1 * t2036 * t2342 - 0.2e1 * t2043 * t2342 + 0.64e2 * t520 * t2047 + 0.8e1 * t2479 * t2056 - 0.8e1 * t2490 * t2056 + 0.8e1 * t2493 * t2056 - 0.2e1 * t2069 * t2063 + 0.8e1 * t2073 * t2063 + 0.2e1 * t2077 * t2063 - 0.8e1 * t2479 * t2063 + 0.8e1 * t2490 * t2063 - 0.8e1 * t2493 * t2063;
                        pr = t1991 + t2039 + t2084 + t2118 + t2152 + t2185 + t2215 + t2244 + t2277 + t2306 + t2336 + t2371 + t2401 + t2436 + t2468 + t2500;
            // (u1,u2) = ur*rHat + vr*thetaHat 
                        u(i1,i2,i3,u1c) = ur*cosTheta - vr*sinTheta; 
                        u(i1,i2,i3,u2c) = ur*sinTheta + vr*cosTheta;
                        u(i1,i2,i3,pc ) = pr;
                        u1Max=max(u1Max,u(i1,i2,i3,u1c));
                        u2Max=max(u2Max,u(i1,i2,i3,u2c));
                        pMax =max(pMax, u(i1,i2,i3,pc));
                    }
                    uNorm = max(u1Max,u2Max,u3Max,pMax);
          // if( scale==0. )
          // {
          //   printF("UDKS: ERROR: annulus solution is zero! scale=%e\n",scale);
          //   OV_ABORT("ERROR");
          // }
          // scale = 1./scale;
          // c4 *= scale;  // set this so the velocity and acceleration are scaled too
          // FOR_3(i1,i2,i3,I1,I2,I3)
          // {
          //   u(i1,i2,i3,u1c) *= scale;
          //   u(i1,i2,i3,u2c) *= scale;
          //   u(i1,i2,i3,pc ) *= scale;
          // }    
                }
                else
                {
                    OV_ABORT("UDKS:ERROR: annulus solution : 3D called");
                }

                if( computeVelocity )
                {
                    printF("\n $$$$$$$$$$ userDefinedKnownSolution: EVALUATE VELOCITY, ACCEL and PRESSURE AT PAST TIMES t=%9.3e dt=%16.6e $$$$$$$\n\n",t,dt);
                    assert( dt>0. );
                    assert( vgf!=NULL );
                    const int currentVelocity = dbase.get<int>("currentVelocity");
                    assert( currentVelocity>=0 && currentVelocity<numberOfVelocityFunctions );
          // fill in t0=0 and some past time levels
                    for( int level=0; level<numberOfVelocityFunctions; level++ )
                    {
                        real tc = 0. - level*dt;
                        const int prevVelocity = (currentVelocity - level + numberOfVelocityFunctions) % numberOfVelocityFunctions;
                        OV_GET_SERIAL_ARRAY(real,vgf[prevVelocity][grid],vLocal);
                        OV_GET_SERIAL_ARRAY(real,vtgf[prevVelocity][grid],vtLocal);
                        FOR_3D(i1,i2,i3,I1,I2,I3)
                        {
                            Real x = vertex(i1,i2,i3,0);
                            Real y = vertex(i1,i2,i3,1);
                            Real r = sqrt( x*x + y*y );
                            Real theta = atan2(y,x); 
                            Real cosTheta = x/r;
                            Real sinTheta = y/r;
                            t1 = omega * tc;
                            t2 = cos(t1);
                            t3 = n * theta;
                            t4 = sin(t3);
                            t6 = sin(t1);
                            t7 = cos(t3);
                            t11 = omega / 0.2e1;
                            t12 = yn(n, t11);
                            t13 = n * n;
                            t15 = omega * omega;
                            t16 = pow(r, -n);
                            t17 = pow(r, n);
                            t20 = n + 0.1e1;
                            t22 = pow(0.2e1, 0.3e1 * t20);
                            t23 = t22 * (t16 - t17) * t15;
                            t27 = pow(0.2e1, 0.2e1 * t20);
                            t28 = omega * r;
                            t29 = yn(n, t28);
                            t31 = t15 * t15;
                            t34 = pow(0.8e1, n);
                            t35 = n * t34;
                            t36 = t13 * t13;
                            t37 = t13 * n;
                            t40 = t13 * (-0.5e1 / 0.8e1 * t15 - 0.1e1);
                            t41 = t31 / 0.16e2;
                            t42 = t15 / 0.2e1;
                            t44 = t16 * (t36 - t37 + t40 + n + t41 + t42);
                            t46 = t17 * (t36 + t37 + t40 - n + t41 + t42);
                            t47 = t44 + t46;
                            t50 = yn(t20, omega);
                            t52 = yn(n, omega);
                            t53 = n - 0.1e1;
                            t55 = -t50 * omega + t53 * t52;
                            t56 = pow(0.4e1, n);
                            t57 = t56 * t55;
                            t58 = t15 / 0.8e1;
                            t59 = t37 - t58 - n;
                            t62 = t16 * t59 * n * t57;
                            t63 = t13 - t42 - n;
                            t64 = t37 + t58 - n;
                            t66 = pow(0.16e2, n);
                            t67 = t66 * t64 * t63;
                            t68 = n * t15;
                            t73 = n * t56;
                            t74 = t73 * t20 * (t37 - t68 / 0.2e1 + 0.3e1 / 0.8e1 * t15 - n);
                            t75 = -t67 - t74;
                            t76 = t29 * t75;
                            t77 = t17 * t55;
                            t80 = t66 * t64 * n * t77;
                            t83 = jn(t20, t11);
                            t86 = jn(n, t11);
                            t90 = jn(n, t28);
                            t96 = jn(t20, omega);
                            t98 = jn(n, omega);
                            t100 = -t96 * omega + t53 * t98;
                            t103 = t16 * t100 * t59 * t73;
                            t104 = t90 * t75;
                            t108 = t100 * t66 * t64 * n * t17;
                            t111 = yn(t20, t11);
                            t115 = t13 - t15 / 0.16e2 - n;
                            t119 = t16 * t68 * t115 * t57 / 0.8e1;
                            t123 = t36 + t13 * (-t58 - 0.1e1) + t31 / 0.128e3;
                            t126 = t115 * t56;
                            t131 = t66 * t123 * t63 + (t13 - t42 + n) * t15 * t126 / 0.8e1;
                            t132 = t29 * t131;
                            t135 = t66 * n * t123 * t77;
                            t142 = t16 * t100 * t15 * n * t126 / 0.8e1;
                            t144 = -t90 * t131;
                            t148 = t100 * t66 * n * t123 * t17;
                            t156 = t13 - t58 - 0.1e1;
                            t157 = t37 * t156;
                            t158 = t157 * (t100 * t29 - t90 * t55) * t34;
                            t161 = 0.1e1 / r;
                            t169 = t86 * t131;
                            t172 = t100 * t37 * t156 * t34;
                            t174 = 0.1e1 / (-t83 * omega * (t31 * t27 / 0.64e2 + t67 + t74) / 0.4e1 + t169 - t172);
                            ur = -t174 * t161 * omega * (t83 * omega * (t23 * t13 * t12 / 0.64e2 - t31 * t29 * t27 / 0.64e2 + t12 * t47 * t35 - t62 + t76 + t80) / 0.4e1 - t111 * omega * (t23 * t13 * t86 / 0.64e2 - t31 * t90 * t27 / 0.64e2 + t86 * t47 * t35 - t103 + t104 + t108) / 0.4e1 + t86 * (t119 + t132 - t135) + t12 * (-t142 + t144 + t148) - t158) * (-t4 * t2 + t7 * t6) * c4;
                            t177 = yn(t20, t28);
                            t187 = t22 * (t17 + t16) * t15;
                            t191 = -r * t75;
                            t192 = t177 * omega;
                            t195 = (t44 - t46) * t34;
                            t204 = jn(t20, t28);
                            t214 = t204 * omega;
                            vr = t174 / n * t161 * (t7 * t2 + t4 * t6) * omega * (t83 * omega * (-t27 * (t29 * n - t177 * t28) * t31 / 0.64e2 - t187 * t37 * t12 / 0.64e2 + t192 * t191 + (-t12 * n * t195 + t62 + t76 + t80) * n) / 0.4e1 - t111 * (-t27 * (t90 * n - t204 * t28) * t31 / 0.64e2 - t187 * t37 * t86 / 0.64e2 + t214 * t191 + (-t86 * n * t195 + t103 + t104 + t108) * n) * omega / 0.4e1 + t214 * (-t157 * t55 * t34 + t12 * t131) * r - t192 * r * (t169 - t172) + (t86 * (-t119 + t132 - t135) + t12 * (t142 + t144 + t148) - t158) * n) * c4;
              // (u1,u2) = ur*rHat + vr*thetaHat 
                            vLocal(i1,i2,i3,u1c) = ur*cosTheta - vr*sinTheta;
                            vLocal(i1,i2,i3,u2c) = ur*sinTheta + vr*cosTheta;
                            t1 = omega / 0.2e1;
                            t2 = yn(n, t1);
                            t3 = n * n;
                            t5 = omega * omega;
                            t6 = pow(r, -n);
                            t7 = pow(r, n);
                            t10 = n + 0.1e1;
                            t12 = pow(0.2e1, 0.3e1 * t10);
                            t13 = t12 * (t6 - t7) * t5;
                            t17 = pow(0.2e1, 0.2e1 * t10);
                            t18 = omega * r;
                            t19 = yn(n, t18);
                            t21 = t5 * t5;
                            t24 = pow(0.8e1, n);
                            t25 = n * t24;
                            t26 = t3 * t3;
                            t27 = t3 * n;
                            t30 = t3 * (-0.5e1 / 0.8e1 * t5 - 0.1e1);
                            t31 = t21 / 0.16e2;
                            t32 = t5 / 0.2e1;
                            t34 = t6 * (t26 - t27 + t30 + n + t31 + t32);
                            t36 = t7 * (t26 + t27 + t30 - n + t31 + t32);
                            t37 = t34 + t36;
                            t40 = yn(t10, omega);
                            t42 = yn(n, omega);
                            t43 = n - 0.1e1;
                            t45 = -t40 * omega + t43 * t42;
                            t46 = pow(0.4e1, n);
                            t47 = t46 * t45;
                            t48 = t5 / 0.8e1;
                            t49 = t27 - t48 - n;
                            t52 = t6 * t49 * n * t47;
                            t53 = t3 - t32 - n;
                            t54 = t27 + t48 - n;
                            t56 = pow(0.16e2, n);
                            t57 = t56 * t54 * t53;
                            t58 = n * t5;
                            t63 = n * t46;
                            t64 = t63 * t10 * (t27 - t58 / 0.2e1 + 0.3e1 / 0.8e1 * t5 - n);
                            t65 = -t57 - t64;
                            t66 = t19 * t65;
                            t67 = t7 * t45;
                            t70 = t56 * t54 * n * t67;
                            t73 = jn(t10, t1);
                            t76 = jn(n, t1);
                            t80 = jn(n, t18);
                            t86 = jn(t10, omega);
                            t88 = jn(n, omega);
                            t90 = -t86 * omega + t43 * t88;
                            t93 = t6 * t90 * t49 * t63;
                            t94 = t80 * t65;
                            t98 = t90 * t56 * t54 * n * t7;
                            t101 = yn(t10, t1);
                            t105 = t3 - t5 / 0.16e2 - n;
                            t109 = t6 * t58 * t105 * t47 / 0.8e1;
                            t113 = t26 + t3 * (-t48 - 0.1e1) + t21 / 0.128e3;
                            t116 = t105 * t46;
                            t121 = t56 * t113 * t53 + (t3 - t32 + n) * t5 * t116 / 0.8e1;
                            t122 = t19 * t121;
                            t125 = t56 * n * t113 * t67;
                            t132 = t6 * t90 * t5 * n * t116 / 0.8e1;
                            t134 = -t80 * t121;
                            t138 = t90 * t56 * n * t113 * t7;
                            t146 = t3 - t48 - 0.1e1;
                            t147 = t27 * t146;
                            t148 = t147 * (t90 * t19 - t80 * t45) * t24;
                            t152 = omega * tc;
                            t153 = cos(t152);
                            t154 = n * theta;
                            t155 = cos(t154);
                            t157 = sin(t152);
                            t158 = sin(t154);
                            t161 = 0.1e1 / r;
                            t169 = t76 * t121;
                            t172 = t90 * t27 * t146 * t24;
                            t174 = 0.1e1 / (-t73 * omega * (t21 * t17 / 0.64e2 + t57 + t64) / 0.4e1 + t169 - t172);
                            ur = -t174 * t161 * (t155 * t153 + t158 * t157) * t5 * (t73 * omega * (t13 * t3 * t2 / 0.64e2 - t21 * t19 * t17 / 0.64e2 + t2 * t37 * t25 - t52 + t66 + t70) / 0.4e1 - t101 * omega * (t13 * t3 * t76 / 0.64e2 - t21 * t80 * t17 / 0.64e2 + t76 * t37 * t25 - t93 + t94 + t98) / 0.4e1 + t76 * (t109 + t122 - t125) + t2 * (-t132 + t134 + t138) - t148) * c4;
                            t177 = yn(t10, t18);
                            t187 = t12 * (t7 + t6) * t5;
                            t191 = -r * t65;
                            t192 = t177 * omega;
                            t195 = (t34 - t36) * t24;
                            t204 = jn(t10, t18);
                            t214 = t204 * omega;
                            vr = -0.1e1 / n * t174 * t161 * t5 * (-t158 * t153 + t155 * t157) * c4 * (t73 * omega * (-t17 * (t19 * n - t177 * t18) * t21 / 0.64e2 - t187 * t27 * t2 / 0.64e2 + t192 * t191 + (-t2 * n * t195 + t52 + t66 + t70) * n) / 0.4e1 - t101 * (-t17 * (t80 * n - t204 * t18) * t21 / 0.64e2 - t187 * t27 * t76 / 0.64e2 + t214 * t191 + (-t76 * n * t195 + t93 + t94 + t98) * n) * omega / 0.4e1 + t214 * (-t147 * t45 * t24 + t2 * t121) * r - t192 * r * (t169 - t172) + (t76 * (-t109 + t122 - t125) + t2 * (t132 + t134 + t138) - t148) * n);
              // (u1,u2) = ur*rHat + vr*thetaHat 
                            vtLocal(i1,i2,i3,u1c) = ur*cosTheta - vr*sinTheta;
                            vtLocal(i1,i2,i3,u2c) = ur*sinTheta + vr*cosTheta;        
                        }      
            // FOR_3D(i1,i2,i3,I1,I2,I3)
            // {
            //   Real x = vertex(i1,i2,i3,0);
            //   Real y = vertex(i1,i2,i3,1);
            //   vLocal(i1,i2,i3,u1c)  = v1e(x,y,tc);
            //   vLocal(i1,i2,i3,u2c)  = v2e(x,y,tc);
            //   vtLocal(i1,i2,i3,u1c) = v1et(x,y,tc);
            //   vtLocal(i1,i2,i3,u2c) = v2et(x,y,tc);        
            // }         
                    }
          // --- save pressure for time extrapolation ---
                    const int currentPressure = dbase.get<int>("currentPressure");
                    assert( currentPressure>=0 && currentPressure<max(1,numberOfPressureFunctions) );
                    for( int level=0; level<numberOfPressureFunctions; level++ )
                    {
                        real tc = 0. - (level+3)*dt;  // -- save pressure starting at t-3*dt 
                        const int prevPressure = (currentPressure - level + numberOfPressureFunctions) % numberOfPressureFunctions;
                        assert( pgf !=NULL );
                        OV_GET_SERIAL_ARRAY(real,pgf[prevPressure][grid],pLocal);
                        FOR_3D(i1,i2,i3,I1,I2,I3)
                        {
                            Real x = vertex(i1,i2,i3,0);
                            Real y = vertex(i1,i2,i3,1);
                            Real r = sqrt( x*x + y*y );
                            Real theta = atan2(y,x); 
                            Real cosTheta = x/r;
                            Real sinTheta = y/r;
                            t1 = 0.1e1 / r;
                            t2 = jn(n, omega);
                            t3 = n * n;
                            t4 = t3 * t3;
                            t5 = t4 * t3;
                            t8 = t4 * t2;
                            t10 = t4 * n;
                            t13 = t3 * n;
                            t14 = t13 * t2;
                            t16 = omega / 0.2e1;
                            t17 = jn(n, t16);
                            t18 = pow(0.1e1 / 0.2e1, n);
                            t19 = t18 * t17;
                            t20 = omega * omega;
                            t21 = t20 * t20;
                            t22 = t21 * n;
                            t25 = n + 0.1e1;
                            t26 = jn(t25, t16);
                            t27 = t18 * t26;
                            t28 = omega * t4;
                            t31 = t20 * omega;
                            t32 = t31 * t3;
                            t35 = t20 * t3;
                            t38 = omega * t13;
                            t41 = t31 * n;
                            t44 = 0.1e1 / t18;
                            t45 = t17 * t44;
                            t48 = t26 * t44;
                            t53 = t3 * t26;
                            t60 = t20 * t4;
                            t63 = t21 * t3;
                            t66 = omega * t10;
                            t69 = t31 * t13;
                            t74 = 0.64e2 * t18 * omega * t53 - 0.64e2 * t44 * omega * t53 + 0.256e3 * t10 * t2 + 0.14e2 * t22 * t19 - 0.32e2 * t35 * t19 - 0.256e3 * t5 * t2 - 0.2e1 * t22 * t45 - 0.64e2 * t28 * t27 + 0.8e1 * t32 * t27 + 0.64e2 * t38 * t27 - 0.24e2 * t41 * t27 + 0.128e3 * t35 * t45 + 0.64e2 * t38 * t48 - 0.24e2 * t41 * t48 - 0.160e3 * t60 * t45 + 0.18e2 * t63 * t45 - 0.64e2 * t66 * t48 + 0.32e2 * t69 * t48 - 0.256e3 * t14 + 0.256e3 * t8;
                            t91 = t21 * t20;
                            t93 = t21 * omega;
                            t98 = jn(t25, omega);
                            t99 = t13 * t98;
                            t122 = 0.32e2 * t44 * t20 * t13 * t17 + 0.256e3 * omega * t10 * t98 - 0.256e3 * t44 * t10 * t17 + 0.256e3 * t44 * t5 * t17 - 0.256e3 * omega * t99 + 0.256e3 * t13 * t45 - 0.32e2 * t20 * t14 + 0.32e2 * t60 * t19 - 0.18e2 * t63 * t19 + t91 * t19 + 0.32e2 * t20 * t8 - 0.64e2 * t66 * t27 + 0.32e2 * t69 * t27 - 0.4e1 * t93 * t27 + 0.64e2 * t28 * t48 - 0.32e2 * t31 * t99 - 0.8e1 * t32 * t48 - 0.256e3 * t4 * t45 - t91 * t45 + 0.4e1 * t93 * t48;
                            t124 = 0.1e1 / (t74 + t122);
                            t125 = t124 * t1;
                            t126 = pow(r, n);
                            t127 = 0.1e1 / t126;
                            t128 = t127 * t18;
                            t129 = omega * tc;
                            t130 = sin(t129);
                            t132 = t130 * t128 * t125;
                            t133 = n * theta;
                            t134 = sin(t133);
                            t135 = yn(n, t16);
                            t136 = t135 * t134;
                            t137 = t98 * t136;
                            t138 = t3 * c4;
                            t139 = t31 * t138;
                            t140 = t139 * t137;
                            t143 = t2 * t136;
                            t144 = t13 * c4;
                            t145 = t20 * t144;
                            t146 = t145 * t143;
                            t149 = n * c4;
                            t150 = t21 * t149;
                            t151 = t150 * t143;
                            t154 = yn(t25, t16);
                            t155 = t154 * t134;
                            t156 = t98 * t155;
                            t157 = t150 * t156;
                            t160 = t17 * t134;
                            t161 = yn(t25, omega);
                            t162 = t161 * t160;
                            t163 = t139 * t162;
                            t166 = t126 * t44;
                            t168 = t130 * t166 * t125;
                            t171 = cos(t129);
                            t173 = t171 * t128 * t125;
                            t174 = cos(t133);
                            t175 = t26 * t174;
                            t176 = yn(n, omega);
                            t177 = t176 * t175;
                            t178 = omega * t138;
                            t179 = t178 * t177;
                            t183 = t171 * t166 * t125;
                            t184 = t154 * t174;
                            t185 = t98 * t184;
                            t186 = t20 * t138;
                            t187 = t186 * t185;
                            t190 = t130 * t126;
                            t191 = t134 * t190;
                            t192 = t191 * t125;
                            t193 = t154 * t17;
                            t194 = t139 * t193;
                            t197 = t171 * t126;
                            t198 = t174 * t197;
                            t199 = t198 * t125;
                            t200 = t135 * t26;
                            t201 = t178 * t200;
                            t204 = t171 * t18;
                            t205 = omega * r;
                            t206 = jn(n, t205);
                            t208 = t206 * t204 * t125;
                            t209 = t135 * t174;
                            t210 = t186 * t209;
                            t213 = t171 * t44;
                            t215 = t206 * t213 * t125;
                            t216 = t4 * c4;
                            t217 = t20 * t216;
                            t218 = t217 * t209;
                            t221 = t21 * t138;
                            t222 = t221 * t209;
                            t225 = t10 * c4;
                            t226 = omega * t225;
                            t227 = t226 * t184;
                            t230 = t130 * t18;
                            t232 = t206 * t230 * t125;
                            t233 = t221 * t136;
                            t236 = t171 * t127;
                            t237 = t174 * t236;
                            t238 = t237 * t125;
                            t239 = omega * t144;
                            t240 = t239 * t200;
                            t243 = t31 * t149;
                            t244 = t243 * t200;
                            t247 = t145 * t209;
                            t250 = -0.32e2 * t140 * t132 + 0.64e2 * t146 * t132 - 0.2e1 * t151 * t132 + 0.8e1 * t157 * t132 + 0.32e2 * t163 * t132 - 0.2e1 * t151 * t168 - 0.64e2 * t179 * t173 - 0.64e2 * t187 * t183 + 0.8e1 * t194 * t192 - 0.64e2 * t201 * t199 + 0.32e2 * t210 * t208 + 0.18e2 * t222 * t208 + 0.64e2 * t227 * t208 + 0.160e3 * t218 * t215 - 0.32e2 * t247 * t215 + 0.18e2 * t233 * t232 - 0.64e2 * t240 * t238 + 0.32e2 * t244 * t238;
                            t251 = t150 * t209;
                            t254 = t130 * t127;
                            t255 = t134 * t254;
                            t256 = t255 * t125;
                            t257 = t178 * t193;
                            t262 = t226 * t193;
                            t265 = t31 * t144;
                            t266 = t265 * t193;
                            t269 = t130 * t44;
                            t271 = t206 * t269 * t125;
                            t272 = t217 * t136;
                            t281 = t239 * t184;
                            t284 = t243 * t184;
                            t289 = t239 * t193;
                            t292 = t243 * t193;
                            t295 = t265 * t184;
                            t298 = t17 * t174;
                            t299 = c4 * t176;
                            t300 = t5 * t299;
                            t301 = t300 * t298;
                            t304 = c4 * t2;
                            t305 = t5 * t304;
                            t306 = t305 * t209;
                            t309 = t10 * t299;
                            t310 = t309 * t298;
                            t313 = t10 * t304;
                            t314 = t313 * t209;
                            t317 = -0.256e3 * t301 * t183 + 0.256e3 * t306 * t183 + 0.256e3 * t310 * t183 - 0.256e3 * t314 * t183 + 0.64e2 * t201 * t238 + 0.64e2 * t201 * t256 - 0.32e2 * t218 * t208 - 0.32e2 * t295 * t208 - 0.128e3 * t210 * t215 + 0.2e1 * t251 * t215 - 0.64e2 * t281 * t215 + 0.24e2 * t284 * t215 - 0.18e2 * t233 * t271 + 0.64e2 * t289 * t238 - 0.32e2 * t292 * t238 - 0.64e2 * t257 * t256 - 0.64e2 * t262 * t256 + 0.40e2 * t266 * t256 + 0.160e3 * t272 * t271;
                            t319 = t4 * t299;
                            t320 = t319 * t298;
                            t323 = t4 * t304;
                            t324 = t323 * t209;
                            t327 = t13 * t299;
                            t328 = t327 * t298;
                            t331 = t13 * t304;
                            t332 = t331 * t209;
                            t339 = t176 * t160;
                            t340 = t145 * t339;
                            t343 = t150 * t339;
                            t346 = t26 * t134;
                            t347 = t161 * t346;
                            t348 = t150 * t347;
                            t351 = t161 * t298;
                            t352 = t265 * t351;
                            t355 = t93 * t149;
                            t356 = t355 * t351;
                            t359 = t176 * t298;
                            t360 = t217 * t359;
                            t363 = t221 * t359;
                            t366 = t161 * t175;
                            t367 = t217 * t366;
                            t370 = t2 * t155;
                            t371 = t239 * t370;
                            t374 = t243 * t370;
                            t377 = t176 * t346;
                            t378 = omega * t216;
                            t379 = t378 * t377;
                            t382 = t139 * t377;
                            t385 = -0.64e2 * t340 * t132 + 0.2e1 * t343 * t132 - 0.8e1 * t348 * t132 + 0.64e2 * t379 * t132 + 0.8e1 * t382 * t132 + 0.64e2 * t371 * t168 + 0.8e1 * t374 * t168 - 0.32e2 * t352 * t173 + 0.2e1 * t356 * t173 + 0.32e2 * t360 * t173 - 0.2e1 * t363 * t173 + 0.64e2 * t367 * t173 + 0.256e3 * t320 * t183 - 0.256e3 * t324 * t183 - 0.256e3 * t328 * t183 + 0.256e3 * t332 * t183 + 0.64e2 * t257 * t199 - 0.18e2 * t222 * t215;
                            t386 = t217 * t339;
                            t389 = t221 * t339;
                            t392 = t217 * t347;
                            t395 = t226 * t377;
                            t398 = t226 * t137;
                            t401 = t265 * t137;
                            t404 = t243 * t377;
                            t407 = t355 * t137;
                            t410 = t217 * t143;
                            t413 = t221 * t143;
                            t416 = t217 * t156;
                            t419 = t226 * t370;
                            t422 = t265 * t162;
                            t425 = t355 * t162;
                            t430 = t2 * t184;
                            t431 = t139 * t430;
                            t434 = t239 * t351;
                            t437 = t186 * t366;
                            t440 = t239 * t177;
                            t443 = 0.32e2 * t386 * t132 - 0.8e1 * t404 * t132 - 0.32e2 * t422 * t132 + 0.2e1 * t425 * t132 + 0.32e2 * t386 * t168 - 0.2e1 * t389 * t168 - 0.64e2 * t392 * t168 + 0.64e2 * t395 * t168 - 0.256e3 * t398 * t168 + 0.32e2 * t401 * t168 - 0.2e1 * t407 * t168 - 0.32e2 * t410 * t168 + 0.2e1 * t413 * t168 + 0.64e2 * t416 * t168 - 0.64e2 * t419 * t168 - 0.8e1 * t431 * t173 - 0.256e3 * t434 * t183 + 0.64e2 * t437 * t183 - 0.64e2 * t440 * t183;
                            t446 = t243 * t177;
                            t449 = t98 * t209;
                            t450 = t239 * t449;
                            t457 = t178 * t155;
                            t460 = t145 * t136;
                            t463 = t150 * t136;
                            t466 = t378 * t155;
                            t469 = t139 * t155;
                            t472 = t331 * t136;
                            t475 = t300 * t160;
                            t484 = t319 * t160;
                            t487 = t323 * t136;
                            t490 = t305 * t136;
                            t495 = 0.256e3 * t472 * t168 - 0.256e3 * t475 * t168 + 0.256e3 * t484 * t168 - 0.256e3 * t487 * t168 + 0.256e3 * t490 * t168 - 0.8e1 * t446 * t183 + 0.256e3 * t450 * t183 - 0.64e2 * t281 * t208 - 0.14e2 * t463 * t232 + 0.64e2 * t466 * t232 - 0.8e1 * t469 * t232 - 0.64e2 * t240 * t256 + 0.32e2 * t244 * t256 + 0.64e2 * t457 * t271 - 0.32e2 * t460 * t271 + 0.2e1 * t463 * t271 - 0.64e2 * t466 * t271 + 0.8e1 * t469 * t271;
                            t512 = t355 * t193;
                            t517 = t378 * t200;
                            t520 = t139 * t200;
                            t523 = t226 * t200;
                            t526 = t265 * t200;
                            t529 = t355 * t200;
                            t532 = t378 * t193;
                            t541 = 0.8e1 * t194 * t199 - 0.8e1 * t194 * t238 - 0.64e2 * t240 * t199 + 0.32e2 * t244 * t199 + 0.64e2 * t289 * t199 - 0.32e2 * t292 * t199 + 0.64e2 * t517 * t199 - 0.8e1 * t520 * t199 - 0.64e2 * t532 * t199 + 0.24e2 * t284 * t208 - 0.32e2 * t295 * t215 - 0.64e2 * t262 * t238 + 0.40e2 * t266 * t238 - 0.4e1 * t512 * t238 - 0.64e2 * t517 * t238 + 0.8e1 * t520 * t238 + 0.64e2 * t523 * t238 - 0.40e2 * t526 * t238 + 0.4e1 * t529 * t238;
                            t543 = t178 * t184;
                            t546 = t378 * t184;
                            t549 = t139 * t184;
                            t552 = t226 * t155;
                            t555 = t265 * t155;
                            t562 = t230 * t125;
                            t563 = t134 * t206;
                            t564 = c4 * t135;
                            t565 = t91 * t564;
                            t566 = t565 * t563;
                            t568 = t269 * t125;
                            t571 = t206 * t171 * t125;
                            t572 = t161 * t174;
                            t573 = t226 * t572;
                            t576 = t204 * t125;
                            t577 = t174 * t206;
                            t578 = c4 * t154;
                            t579 = t93 * t578;
                            t580 = t579 * t577;
                            t583 = t265 * t572;
                            t586 = t176 * t174;
                            t587 = t217 * t586;
                            t590 = t239 * t572;
                            t593 = t565 * t577;
                            t595 = t213 * t125;
                            t596 = t5 * t564;
                            t597 = t596 * t577;
                            t601 = t10 * t564;
                            t602 = t601 * t577;
                            t607 = -0.64e2 * t543 * t208 - 0.64e2 * t546 * t215 + 0.8e1 * t549 * t215 + 0.64e2 * t552 * t232 - 0.32e2 * t555 * t232 + 0.64e2 * t289 * t256 - 0.32e2 * t292 * t256 - t566 * t562 + t566 * t568 - 0.256e3 * t573 * t571 + 0.32e2 * t583 * t571 - 0.32e2 * t587 * t571 + 0.256e3 * t590 * t571 + 0.4e1 * t580 * t576 - t593 * t576 - 0.4e1 * t580 * t595 + t593 * t595 - 0.256e3 * t597 * t595 + 0.256e3 * t602 * t595;
                            t608 = t4 * t564;
                            t609 = t608 * t577;
                            t612 = t327 * t160;
                            t617 = t186 * t136;
                            t620 = t239 * t155;
                            t623 = t243 * t155;
                            t652 = -0.256e3 * t612 * t168 + 0.40e2 * t266 * t192 - 0.8e1 * t194 * t256 - 0.4e1 * t512 * t199 + 0.64e2 * t523 * t199 - 0.40e2 * t526 * t199 + 0.4e1 * t529 * t199 - 0.32e2 * t272 * t232 - 0.4e1 * t512 * t256 - 0.64e2 * t517 * t256 + 0.8e1 * t520 * t256 + 0.64e2 * t523 * t256 - 0.40e2 * t526 * t256 + 0.4e1 * t529 * t256 + 0.64e2 * t532 * t256 - 0.128e3 * t617 * t271 - 0.64e2 * t620 * t271 + 0.24e2 * t623 * t271 + 0.256e3 * t609 * t595;
                            t680 = t309 * t160;
                            t683 = t313 * t136;
                            t694 = 0.256e3 * t680 * t168 - 0.256e3 * t683 * t168 - 0.64e2 * t240 * t192 + 0.64e2 * t289 * t192 - 0.32e2 * t292 * t192 + 0.64e2 * t517 * t192 - 0.8e1 * t520 * t192 - 0.40e2 * t526 * t192 + 0.4e1 * t529 * t192 - 0.64e2 * t532 * t192 + 0.64e2 * t227 * t215 + 0.64e2 * t543 * t215 + 0.32e2 * t617 * t232 - 0.64e2 * t620 * t232 + 0.24e2 * t623 * t232 - 0.64e2 * t257 * t238 + 0.64e2 * t552 * t271 - 0.32e2 * t555 * t271;
                            t721 = t225 * t586;
                            t725 = t206 * t130 * t125;
                            t726 = t176 * t134;
                            t727 = t144 * t726;
                            t730 = t5 * c4;
                            t731 = t730 * t726;
                            t734 = t225 * t726;
                            t737 = t216 * t726;
                            t740 = t730 * t586;
                            t743 = -0.64e2 * t201 * t192 + 0.32e2 * t244 * t192 + 0.64e2 * t257 * t192 - 0.64e2 * t262 * t192 - 0.4e1 * t512 * t192 + 0.64e2 * t523 * t192 - 0.64e2 * t262 * t199 + 0.40e2 * t266 * t199 - 0.14e2 * t251 * t208 + 0.64e2 * t546 * t208 - 0.8e1 * t549 * t208 - 0.64e2 * t457 * t232 + 0.64e2 * t532 * t238 - 0.256e3 * t721 * t571 + 0.256e3 * t740 * t571 + 0.256e3 * t727 * t725 + 0.256e3 * t731 * t725 - 0.256e3 * t734 * t725 - 0.256e3 * t737 * t725;
                            t745 = t144 * t586;
                            t748 = t216 * t586;
                            t751 = t186 * t143;
                            t754 = t186 * t156;
                            t773 = t178 * t377;
                            t776 = t239 * t162;
                            t779 = t186 * t347;
                            t782 = t239 * t377;
                            t787 = t239 * t137;
                            t790 = -0.64e2 * t371 * t132 + 0.8e1 * t374 * t132 - 0.32e2 * t751 * t132 + 0.64e2 * t754 * t132 + 0.8e1 * t157 * t168 - 0.32e2 * t340 * t168 + 0.2e1 * t343 * t168 - 0.8e1 * t348 * t168 - 0.64e2 * t379 * t168 + 0.8e1 * t382 * t168 - 0.8e1 * t404 * t168 + 0.64e2 * t773 * t168 - 0.256e3 * t776 * t168 + 0.64e2 * t779 * t168 - 0.64e2 * t782 * t168 + 0.256e3 * t787 * t168 + 0.256e3 * t745 * t571 - 0.256e3 * t748 * t571;
                            t795 = t150 * t185;
                            t798 = t378 * t430;
                            t813 = t226 * t177;
                            t816 = t265 * t449;
                            t819 = t355 * t449;
                            t822 = t243 * t430;
                            t825 = t2 * t209;
                            t826 = t145 * t825;
                            t829 = t150 * t825;
                            t834 = t186 * t825;
                            t837 = t139 * t351;
                            t840 = -0.2e1 * t389 * t132 + 0.64e2 * t392 * t132 - 0.64e2 * t395 * t132 + 0.32e2 * t401 * t132 - 0.2e1 * t407 * t132 - 0.32e2 * t410 * t132 + 0.32e2 * t146 * t168 - 0.64e2 * t754 * t168 - 0.8e1 * t446 * t173 + 0.8e1 * t795 * t173 - 0.64e2 * t798 * t173 - 0.64e2 * t813 * t173 + 0.32e2 * t816 * t173 - 0.2e1 * t819 * t173 + 0.8e1 * t822 * t173 + 0.64e2 * t826 * t173 - 0.2e1 * t829 * t173 - 0.32e2 * t834 * t173 + 0.32e2 * t837 * t173;
                            t843 = t145 * t359;
                            t846 = t150 * t359;
                            t849 = t150 * t366;
                            t852 = t378 * t177;
                            t855 = t139 * t177;
                            t858 = t139 * t449;
                            t861 = t217 * t185;
                            t864 = t226 * t430;
                            t881 = t239 * t430;
                            t886 = t178 * t430;
                            t889 = -0.64e2 * t843 * t173 + 0.2e1 * t846 * t173 - 0.8e1 * t849 * t173 + 0.64e2 * t852 * t173 + 0.8e1 * t855 * t173 - 0.32e2 * t858 * t173 - 0.64e2 * t861 * t173 + 0.64e2 * t864 * t173 + 0.64e2 * t886 * t173 + 0.8e1 * t822 * t183 + 0.32e2 * t826 * t183 - 0.2e1 * t829 * t183 - 0.32e2 * t843 * t183 + 0.2e1 * t846 * t183 - 0.8e1 * t849 * t183 - 0.64e2 * t852 * t183 + 0.8e1 * t855 * t183 + 0.64e2 * t881 * t183;
                            t896 = t217 * t825;
                            t899 = t221 * t825;
                            t902 = t378 * t370;
                            t905 = t139 * t370;
                            t908 = t186 * t339;
                            t917 = t178 * t370;
                            t926 = t186 * t359;
                            t933 = t13 * t564;
                            t934 = t933 * t577;
                            t937 = 0.2e1 * t413 * t132 - 0.64e2 * t416 * t132 + 0.64e2 * t419 * t132 - 0.64e2 * t773 * t132 - 0.64e2 * t779 * t132 + 0.64e2 * t782 * t132 - 0.64e2 * t902 * t132 - 0.8e1 * t905 * t132 + 0.32e2 * t908 * t132 + 0.64e2 * t917 * t132 + 0.64e2 * t902 * t168 - 0.8e1 * t905 * t168 - 0.64e2 * t917 * t168 - 0.64e2 * t437 * t173 + 0.64e2 * t440 * t173 - 0.32e2 * t896 * t173 + 0.2e1 * t899 * t173 + 0.32e2 * t926 * t173 - 0.256e3 * t934 * t595;
                            t939 = t145 * t586;
                            t942 = t161 * t134;
                            t943 = t239 * t942;
                            t946 = t145 * t726;
                            t949 = t226 * t942;
                            t952 = t265 * t942;
                            t955 = t217 * t726;
                            t958 = t601 * t563;
                            t961 = t579 * t563;
                            t964 = t596 * t563;
                            t967 = t608 * t563;
                            t970 = t933 * t563;
                            t975 = t171 * t1;
                            t976 = yn(n, t205);
                            t978 = c4 * t174 * t976;
                            t980 = t130 * t1;
                            t982 = c4 * t134 * t976;
                            t994 = 0.64e2 * t187 * t173 - 0.64e2 * t881 * t173 - 0.8e1 * t431 * t183 + 0.8e1 * t795 * t183 + 0.64e2 * t798 * t183 + 0.4e1 * t961 * t562 + 0.256e3 * t958 * t568 - 0.4e1 * t961 * t568 - 0.256e3 * t964 * t568 + 0.256e3 * t967 * t568 - 0.256e3 * t970 * t568 + 0.32e2 * t939 * t571 + 0.256e3 * t943 * t725 + 0.32e2 * t946 * t725 - 0.256e3 * t949 * t725 + 0.32e2 * t952 * t725 - 0.32e2 * t955 * t725 + t978 * t975 + t982 * t980;
                            t995 = t226 * t162;
                            t1006 = t226 * t351;
                            t1021 = t226 * t449;
                            t1036 = 0.256e3 * t1006 * t183 - 0.256e3 * t1021 * t183 - 0.32e2 * t422 * t168 + 0.2e1 * t425 * t168 + 0.256e3 * t995 * t168 + 0.64e2 * t179 * t183 - 0.32e2 * t352 * t183 + 0.2e1 * t356 * t183 + 0.32e2 * t360 * t183 - 0.2e1 * t363 * t183 - 0.64e2 * t367 * t183 + 0.64e2 * t813 * t183 + 0.32e2 * t816 * t183 - 0.2e1 * t819 * t183 + 0.64e2 * t861 * t183 - 0.64e2 * t864 * t183 - 0.64e2 * t886 * t183 - 0.32e2 * t896 * t183 + 0.2e1 * t899 * t183;
                            ur = t250 + t317 + t385 + t443 + t495 + t541 + t607 + t652 + t694 + t743 + t790 + t840 + t889 + t937 + t994 + t1036;
                            t1053 = t134 * t197 * t125;
                            t1061 = t174 * t254 * t125;
                            t1066 = -0.8e1 * t194 * t1053 - 0.64e2 * t517 * t1053 + 0.64e2 * t532 * t1053 + 0.64e2 * t240 * t1061 + 0.32e2 * t292 * t1061 - 0.64e2 * t552 * t208 + 0.32e2 * t555 * t208 - 0.32e2 * t218 * t232 + 0.18e2 * t222 * t232 + 0.64e2 * t227 * t232 - 0.32e2 * t295 * t232;
                            t1074 = t134 * t236 * t125;
                            t1093 = -0.32e2 * t244 * t1061 - 0.64e2 * t257 * t1074 - 0.256e3 * t314 * t168 - 0.256e3 * t328 * t168 - 0.256e3 * t434 * t168 - 0.256e3 * t472 * t183 + 0.256e3 * t612 * t183 - 0.128e3 * t210 * t271 - 0.64e2 * t457 * t215 - 0.64e2 * t281 * t271 + 0.24e2 * t284 * t271 + 0.64e2 * t543 * t271;
                            t1116 = t174 * t190 * t125;
                            t1121 = -0.32e2 * t244 * t1053 - 0.64e2 * t262 * t1074 + 0.40e2 * t266 * t1074 - 0.4e1 * t512 * t1074 + 0.64e2 * t523 * t1074 - 0.40e2 * t526 * t1074 + 0.4e1 * t529 * t1074 + 0.64e2 * t289 * t1116 - 0.32e2 * t292 * t1116 + 0.256e3 * t310 * t168 + 0.256e3 * t332 * t168 - 0.32e2 * t295 * t271;
                            t1146 = -0.64e2 * t240 * t1116 + 0.32e2 * t244 * t1116 - 0.32e2 * t837 * t132 + 0.64e2 * t843 * t132 - 0.2e1 * t846 * t132 + 0.8e1 * t849 * t132 - 0.64e2 * t852 * t132 + 0.32e2 * t163 * t173 - 0.64e2 * t340 * t173 + 0.2e1 * t343 * t173 - 0.8e1 * t348 * t173 - 0.8e1 * t905 * t173;
                            t1173 = 0.64e2 * t179 * t132 - 0.64e2 * t886 * t132 + 0.64e2 * t179 * t168 - 0.8e1 * t849 * t168 - 0.64e2 * t852 * t168 - 0.64e2 * t886 * t168 - 0.32e2 * t896 * t168 + 0.2e1 * t899 * t168 - 0.64e2 * t416 * t173 - 0.64e2 * t773 * t173 + 0.64e2 * t917 * t173 + 0.64e2 * t917 * t183;
                            t1198 = 0.32e2 * t352 * t132 - 0.2e1 * t356 * t132 - 0.32e2 * t360 * t132 + 0.2e1 * t363 * t132 - 0.64e2 * t367 * t132 + 0.64e2 * t813 * t132 - 0.32e2 * t816 * t132 + 0.32e2 * t386 * t173 - 0.2e1 * t389 * t173 + 0.64e2 * t392 * t173 - 0.64e2 * t395 * t173 + 0.64e2 * t419 * t173;
                            t1224 = -0.64e2 * t187 * t132 + 0.8e1 * t446 * t132 - 0.8e1 * t822 * t132 + 0.32e2 * t834 * t132 + 0.64e2 * t881 * t132 + 0.32e2 * t360 * t168 - 0.2e1 * t363 * t168 - 0.64e2 * t367 * t168 + 0.32e2 * t401 * t173 - 0.2e1 * t407 * t173 + 0.32e2 * t340 * t183 - 0.64e2 * t773 * t183;
                            t1249 = -0.256e3 * t1021 * t168 + 0.2e1 * t819 * t132 + 0.64e2 * t861 * t132 - 0.64e2 * t864 * t132 + 0.32e2 * t896 * t132 - 0.2e1 * t899 * t132 + 0.64e2 * t813 * t168 + 0.8e1 * t855 * t168 - 0.32e2 * t386 * t183 + 0.2e1 * t389 * t183 + 0.64e2 * t392 * t183 - 0.64e2 * t395 * t183;
                            t1275 = 0.8e1 * t157 * t173 + 0.8e1 * t795 * t168 + 0.64e2 * t798 * t168 + 0.32e2 * t826 * t168 - 0.2e1 * t829 * t168 - 0.32e2 * t422 * t173 + 0.2e1 * t425 * t173 - 0.64e2 * t779 * t173 + 0.64e2 * t782 * t173 - 0.64e2 * t902 * t173 + 0.32e2 * t908 * t173;
                            t1300 = 0.32e2 * t816 * t168 - 0.2e1 * t819 * t168 + 0.64e2 * t861 * t168 - 0.64e2 * t864 * t168 + 0.256e3 * t398 * t183 + 0.8e1 * t404 * t183 + 0.32e2 * t422 * t183 - 0.2e1 * t425 * t183 + 0.256e3 * t776 * t183 - 0.64e2 * t779 * t183 + 0.64e2 * t782 * t183 - 0.256e3 * t995 * t183;
                            t1326 = 0.256e3 * t1006 * t168 - 0.32e2 * t352 * t168 + 0.2e1 * t356 * t168 - 0.64e2 * t371 * t183 - 0.32e2 * t401 * t183 + 0.2e1 * t407 * t183 + 0.32e2 * t410 * t183 - 0.2e1 * t413 * t183 - 0.64e2 * t416 * t183 + 0.64e2 * t419 * t183 + 0.64e2 * t754 * t183 - 0.256e3 * t787 * t183;
                            t1349 = -0.32e2 * t843 * t168 + 0.2e1 * t846 * t168 - 0.8e1 * t374 * t183 - 0.64e2 * t902 * t183 + 0.8e1 * t905 * t183 + 0.4e1 * t580 * t562 - t566 * t595 + t593 * t568 - 0.256e3 * t597 * t568 - 0.4e1 * t961 * t576 - 0.256e3 * t958 * t595 + 0.256e3 * t964 * t595;
                            t1354 = t44 * t124;
                            t1355 = t1354 * t38;
                            t1356 = jn(t25, t205);
                            t1357 = t130 * t1356;
                            t1359 = c4 * t209 * t1357;
                            t1362 = t3 * omega;
                            t1363 = t1354 * t1362;
                            t1364 = t171 * t1356;
                            t1366 = c4 * t136 * t1364;
                            t1369 = n * t20;
                            t1370 = t1354 * t1369;
                            t1372 = c4 * t155 * t1364;
                            t1378 = c4 * t184 * t1357;
                            t1381 = n * t93;
                            t1382 = t18 * t124;
                            t1383 = t1382 * t1381;
                            t1386 = t1382 * t60;
                            t1389 = t1382 * t63;
                            t1392 = t1354 * t66;
                            t1395 = t1354 * t69;
                            t1398 = t1354 * t1381;
                            t1401 = -0.256e3 * t1359 * t1355 + 0.256e3 * t1359 * t1363 - 0.256e3 * t1366 * t1363 + 0.18e2 * t1366 * t1383 - 0.256e3 * t1366 * t1392 + 0.160e3 * t1366 * t1395 - 0.18e2 * t1366 * t1398 + 0.64e2 * t1372 * t1370 - 0.64e2 * t1378 * t1370 + 0.64e2 * t1372 * t1386 - 0.32e2 * t1372 * t1389 + 0.4e1 * t961 * t595;
                            t1402 = t1382 * t1369;
                            t1405 = t1382 * t69;
                            t1428 = 0.8e1 * t431 * t132 - 0.8e1 * t795 * t132 + 0.64e2 * t798 * t132 - 0.64e2 * t826 * t132 + 0.2e1 * t829 * t132 - 0.8e1 * t855 * t132 + 0.32e2 * t858 * t132 - 0.18e2 * t1359 * t1383 + 0.32e2 * t1359 * t1405 - 0.64e2 * t1378 * t1386 + 0.64e2 * t1378 * t1402 - 0.2e1 * t151 * t173;
                            t1454 = 0.64e2 * t437 * t132 - 0.64e2 * t440 * t132 - 0.32e2 * t926 * t132 - 0.64e2 * t371 * t173 + 0.8e1 * t374 * t173 + 0.64e2 * t379 * t173 + 0.8e1 * t382 * t173 - 0.8e1 * t404 * t173 - 0.32e2 * t410 * t173 + 0.2e1 * t413 * t173 - 0.32e2 * t751 * t173 + 0.64e2 * t754 * t173;
                            t1479 = -0.32e2 * t140 * t173 + 0.64e2 * t146 * t173 + 0.256e3 * t602 * t568 - 0.256e3 * t943 * t571 - 0.32e2 * t946 * t571 + 0.256e3 * t949 * t571 - 0.32e2 * t952 * t571 + 0.32e2 * t955 * t571 - 0.256e3 * t573 * t725 + 0.32e2 * t583 * t725 - 0.32e2 * t587 * t725 + 0.256e3 * t590 * t725;
                            t1496 = 0.1e1 / n;
                            t1497 = t1496 * t91;
                            t1498 = t1354 * t1497;
                            t1503 = t1382 * t41;
                            t1506 = t1382 * t35;
                            t1509 = -0.32e2 * t1359 * t1503 - 0.4e1 * t1372 * t1498 + 0.4e1 * t1378 * t1498 + 0.64e2 * t1378 * t1506 - t593 * t562 + t566 * t576 - 0.4e1 * t580 * t568 + 0.256e3 * t609 * t568 - 0.256e3 * t934 * t568 - 0.256e3 * t967 * t595 + 0.256e3 * t970 * t595;
                            t1534 = 0.64e2 * t262 * t1053 - 0.40e2 * t266 * t1053 + 0.4e1 * t512 * t1053 - 0.64e2 * t523 * t1053 + 0.40e2 * t526 * t1053 - 0.4e1 * t529 * t1053 - 0.64e2 * t262 * t1116 + 0.40e2 * t266 * t1116 - 0.4e1 * t512 * t1116 - 0.64e2 * t1372 * t1402 - 0.256e3 * t301 * t168 + 0.256e3 * t306 * t168;
                            t1558 = t13 * t20;
                            t1559 = t1354 * t1558;
                            t1562 = 0.64e2 * t257 * t1061 + 0.64e2 * t201 * t1074 + 0.64e2 * t523 * t1116 - 0.40e2 * t526 * t1116 + 0.4e1 * t529 * t1116 + 0.32e2 * t1378 * t1389 + 0.64e2 * t1378 * t1559 - 0.18e2 * t233 * t208 + 0.32e2 * t272 * t208 + 0.18e2 * t233 * t215 - 0.160e3 * t272 * t215 - 0.18e2 * t222 * t271;
                            t1563 = t1354 * t22;
                            t1566 = t21 * t31;
                            t1567 = t1496 * t1566;
                            t1568 = t1354 * t1567;
                            t1570 = t1354 * t28;
                            t1573 = t1354 * t32;
                            t1586 = t1354 * t60;
                            t1589 = t1354 * t63;
                            t1592 = t1354 * t35;
                            t1595 = 0.18e2 * t1359 * t1398 - t1359 * t1568 - 0.256e3 * t1359 * t1570 + 0.32e2 * t1359 * t1573 + 0.256e3 * t1366 * t1570 - 0.32e2 * t1366 * t1573 - 0.64e2 * t1372 * t1559 + 0.8e1 * t1372 * t1563 - 0.8e1 * t1378 * t1563 - 0.64e2 * t1378 * t1586 + 0.32e2 * t1378 * t1589 + 0.64e2 * t1378 * t1592;
                            t1598 = t1382 * t1497;
                            t1605 = t1354 * t41;
                            t1610 = t1382 * t1558;
                            t1613 = t1382 * t22;
                            t1624 = 0.256e3 * t1366 * t1355 - 0.128e3 * t1366 * t1605 + 0.64e2 * t1372 * t1586 - 0.32e2 * t1372 * t1589 + 0.4e1 * t1372 * t1598 + 0.64e2 * t1372 * t1610 - 0.8e1 * t1372 * t1613 - 0.4e1 * t1378 * t1598 - 0.64e2 * t1378 * t1610 + 0.8e1 * t1378 * t1613 - t982 * t975 + t978 * t980;
                            t1637 = t1382 * t1567;
                            t1641 = t1356 * t124;
                            t1642 = t1641 * t63;
                            t1643 = t174 * t130;
                            t1644 = c4 * t161;
                            t1645 = t1644 * t1643;
                            t1648 = t1641 * t66;
                            t1649 = t299 * t1643;
                            t1652 = t1641 * t28;
                            t1655 = t1641 * t32;
                            t1658 = 0.256e3 * t1359 * t1392 - 0.160e3 * t1359 * t1395 + 0.128e3 * t1359 * t1605 - 0.32e2 * t1366 * t1405 + 0.32e2 * t1366 * t1503 - t1366 * t1637 - 0.64e2 * t1372 * t1506 - 0.32e2 * t1645 * t1642 - 0.256e3 * t1649 * t1648 + 0.256e3 * t1649 * t1652 - 0.32e2 * t1649 * t1655 + 0.32e2 * t939 * t725;
                            t1660 = t1641 * t69;
                            t1663 = t134 * t171;
                            t1664 = t299 * t1663;
                            t1669 = t1641 * t60;
                            t1670 = t1644 * t1663;
                            t1681 = t124 * t21;
                            t1682 = t1356 * t18;
                            t1683 = t1682 * t1681;
                            t1684 = t578 * t1643;
                            t1687 = t124 * t93;
                            t1688 = t1682 * t1687;
                            t1689 = t564 * t1663;
                            t1692 = t564 * t1643;
                            t1695 = t578 * t1663;
                            t1698 = 0.32e2 * t1670 * t1642 + 0.256e3 * t1645 * t1669 + 0.256e3 * t1664 * t1648 + 0.32e2 * t1649 * t1660 - 0.256e3 * t1664 * t1652 + 0.32e2 * t1664 * t1655 - 0.32e2 * t1664 * t1660 - 0.256e3 * t1670 * t1669 - 0.24e2 * t1684 * t1683 + 0.24e2 * t1695 * t1683 - 0.14e2 * t1689 * t1688 + 0.14e2 * t1692 * t1688;
                            t1723 = -0.64e2 * t201 * t1061 - 0.8e1 * t520 * t1061 + 0.14e2 * t463 * t208 - 0.64e2 * t466 * t208 - 0.64e2 * t552 * t215 + 0.32e2 * t555 * t215 + 0.160e3 * t218 * t271 + 0.256e3 * t734 * t571 + 0.256e3 * t737 * t571 - 0.256e3 * t721 * t725 + 0.256e3 * t740 * t725 - 0.256e3 * t748 * t725;
                            t1749 = 0.64e2 * t240 * t1053 - 0.64e2 * t289 * t1053 + 0.32e2 * t292 * t1053 - 0.64e2 * t240 * t1074 + 0.32e2 * t244 * t1074 + 0.64e2 * t289 * t1074 - 0.32e2 * t292 * t1074 + 0.8e1 * t469 * t208 - 0.14e2 * t251 * t232 + 0.64e2 * t546 * t232 - 0.8e1 * t549 * t232;
                            t1774 = -0.64e2 * t289 * t1061 - 0.64e2 * t187 * t168 - 0.8e1 * t431 * t168 + 0.64e2 * t437 * t168 - 0.64e2 * t440 * t168 - 0.8e1 * t446 * t168 + 0.256e3 * t450 * t168 + 0.8e1 * t822 * t168 + 0.64e2 * t881 * t168 - 0.2e1 * t343 * t183 + 0.8e1 * t348 * t183 + 0.64e2 * t227 * t271;
                            t1796 = t1356 * t44;
                            t1797 = t1796 * t1687;
                            t1800 = t1359 * t1637 + t1366 * t1568 - 0.64e2 * t1372 * t1592 - 0.32e2 * t146 * t183 + 0.2e1 * t151 * t183 - 0.8e1 * t157 * t183 + 0.2e1 * t1689 * t1797 + 0.64e2 * t379 * t183 - 0.8e1 * t382 * t183 - 0.256e3 * t727 * t571 - 0.256e3 * t731 * t571 + 0.256e3 * t745 * t725;
                            t1801 = t1796 * t1681;
                            t1808 = t1641 * t35;
                            t1811 = t1641 * t38;
                            t1818 = t1641 * t1362;
                            t1829 = 0.64e2 * t262 * t1061 - 0.40e2 * t266 * t1061 + 0.4e1 * t512 * t1061 - 0.256e3 * t1645 * t1808 + 0.256e3 * t1649 * t1811 - 0.256e3 * t1649 * t1818 - 0.256e3 * t1664 * t1811 + 0.256e3 * t1664 * t1818 + 0.256e3 * t1670 * t1808 - 0.24e2 * t1684 * t1801 - 0.2e1 * t1692 * t1797 + 0.24e2 * t1695 * t1801;
                            t1856 = 0.8e1 * t194 * t1061 + 0.64e2 * t517 * t1061 - 0.64e2 * t523 * t1061 + 0.40e2 * t526 * t1061 - 0.4e1 * t529 * t1061 - 0.64e2 * t532 * t1061 + 0.32e2 * t460 * t215 - 0.2e1 * t463 * t215 + 0.64e2 * t466 * t215 - 0.8e1 * t469 * t215 - 0.32e2 * t247 * t271 + 0.2e1 * t251 * t271;
                            t1881 = 0.256e3 * t320 * t168 - 0.256e3 * t324 * t168 - 0.256e3 * t484 * t183 + 0.256e3 * t487 * t183 - 0.256e3 * t680 * t183 + 0.32e2 * t210 * t232 + 0.128e3 * t617 * t215 + 0.64e2 * t620 * t215 - 0.24e2 * t623 * t215 - 0.64e2 * t281 * t232 + 0.24e2 * t284 * t232 - 0.64e2 * t546 * t271;
                            t1907 = 0.8e1 * t520 * t1053 - 0.8e1 * t194 * t1074 - 0.64e2 * t517 * t1074 + 0.8e1 * t520 * t1074 + 0.64e2 * t532 * t1074 + 0.8e1 * t194 * t1116 + 0.64e2 * t517 * t1116 - 0.8e1 * t520 * t1116 - 0.64e2 * t532 * t1116 + 0.256e3 * t475 * t183 + 0.256e3 * t683 * t183 + 0.8e1 * t549 * t271;
                            t1929 = yn(t25, t205);
                            t1930 = t1929 * t1496 * omega;
                            t1935 = -c4 * t1643 * t1930 + c4 * t1663 * t1930 + 0.64e2 * t201 * t1053 - 0.64e2 * t257 * t1053 - 0.64e2 * t201 * t1116 + 0.64e2 * t257 * t1116 - 0.256e3 * t490 * t183 + 0.64e2 * t457 * t208 - 0.32e2 * t617 * t208 + 0.64e2 * t620 * t208 - 0.24e2 * t623 * t208 - 0.64e2 * t543 * t232;
                            vr = t1935 + t1749 + t1173 + t1349 + t1249 + t1698 + t1198 + t1774 + t1479 + t1881 + t1146 + t1658 + t1723 + t1224 + t1907 + t1829 + t1800 + t1624 + t1401 + t1562 + t1428 + t1121 + t1300 + t1275 + t1509 + t1534 + t1093 + t1454 + t1066 + t1856 + t1595 + t1326;
                            t1940 = t198 * t1354;
                            t1941 = t161 * t26;
                            t1942 = t150 * t1941;
                            t1945 = t176 * t26;
                            t1946 = t139 * t1945;
                            t1949 = t98 * t135;
                            t1950 = t139 * t1949;
                            t1953 = t2 * t135;
                            t1954 = t145 * t1953;
                            t1957 = t98 * t154;
                            t1958 = t150 * t1957;
                            t1961 = t2 * t154;
                            t1962 = t139 * t1961;
                            t1965 = t237 * t1382;
                            t1966 = t176 * t17;
                            t1967 = t150 * t1966;
                            t1970 = t255 * t1382;
                            t1971 = t355 * t1961;
                            t1974 = t21 * t144;
                            t1975 = t1974 * t1966;
                            t1978 = t1974 * t1941;
                            t1981 = t31 * t216;
                            t1982 = t1981 * t1945;
                            t1985 = t1974 * t1953;
                            t1988 = t1974 * t1957;
                            t1991 = 0.64e2 * t1942 * t1940 - 0.64e2 * t1946 * t1940 + 0.256e3 * t1950 * t1940 - 0.256e3 * t1954 * t1940 - 0.64e2 * t1958 * t1940 + 0.64e2 * t1962 * t1940 - 0.32e2 * t1967 * t1965 - 0.32e2 * t1975 * t1965 - 0.64e2 * t1978 * t1965 + 0.64e2 * t1982 * t1965 + 0.32e2 * t1985 * t1965 + 0.64e2 * t1988 * t1965 + 0.8e1 * t1971 * t1970;
                            t1992 = t1981 * t1961;
                            t1995 = t217 * t1966;
                            t1998 = t217 * t1953;
                            t2001 = t191 * t1354;
                            t2006 = t221 * t1966;
                            t2009 = t355 * t1945;
                            t2012 = t221 * t1953;
                            t2017 = t161 * t17;
                            t2018 = t355 * t2017;
                            t2021 = t126 * t124;
                            t2022 = t134 * t130;
                            t2023 = t2022 * t2021;
                            t2024 = t93 * t138;
                            t2025 = t2024 * t193;
                            t2030 = t254 * t1382;
                            t2031 = t93 * t304;
                            t2032 = t2031 * t155;
                            t2035 = t1566 * t1644;
                            t2036 = t2035 * t160;
                            t2039 = -0.8e1 * t1971 * t1940 + 0.256e3 * t1995 * t1940 - 0.256e3 * t1998 * t1940 - 0.32e2 * t2006 * t1940 + 0.8e1 * t2009 * t1940 + 0.32e2 * t2012 * t1940 - 0.64e2 * t1992 * t1965 - 0.32e2 * t2018 * t1965 + 0.256e3 * t1995 * t2001 - 0.256e3 * t1998 * t2001 + 0.40e2 * t2025 * t2023 - 0.64e2 * t266 * t2023 - 0.8e1 * t2032 * t2030 - 0.2e1 * t2036 * t2030;
                            t2041 = c4 * t98;
                            t2042 = t1566 * t2041;
                            t2043 = t2042 * t136;
                            t2046 = t127 * t124;
                            t2047 = t2022 * t2046;
                            t2056 = t197 * t1354;
                            t2057 = t2035 * t298;
                            t2060 = t2042 * t209;
                            t2063 = t236 * t1382;
                            t2068 = t91 * t299;
                            t2069 = t2068 * t298;
                            t2072 = t91 * t1644;
                            t2073 = t2072 * t175;
                            t2076 = t91 * t304;
                            t2077 = t2076 * t209;
                            t2080 = t174 * t171;
                            t2081 = t2080 * t2021;
                            t2084 = 0.2e1 * t2043 * t2030 - 0.64e2 * t244 * t2047 + 0.64e2 * t292 * t2047 + 0.8e1 * t512 * t2047 - 0.8e1 * t529 * t2047 + 0.2e1 * t2057 * t2056 - 0.2e1 * t2060 * t2056 + 0.2e1 * t2069 * t2056 - 0.8e1 * t2073 * t2056 - 0.2e1 * t2077 * t2056 - 0.2e1 * t2057 * t2063 + 0.2e1 * t2060 * t2063 + 0.64e2 * t292 * t2081;
                            t2087 = t2080 * t2046;
                            t2106 = t150 * t1953;
                            t2109 = t2024 * t2017;
                            t2114 = t91 * t149;
                            t2115 = t2114 * t1966;
                            t2118 = 0.64e2 * t194 * t2081 - 0.64e2 * t194 * t2087 + 0.32e2 * t1975 * t1940 - 0.32e2 * t2109 * t1940 - 0.2e1 * t2115 * t1940 + 0.64e2 * t1942 * t1970 - 0.64e2 * t1946 * t1970 + 0.32e2 * t2106 * t1970 - 0.64e2 * t244 * t2081 + 0.8e1 * t512 * t2081 - 0.64e2 * t520 * t2081 + 0.8e1 * t512 * t2087 + 0.64e2 * t520 * t2087 - 0.8e1 * t529 * t2087;
                            t2131 = t186 * t1966;
                            t2134 = t243 * t1945;
                            t2137 = t186 * t1953;
                            t2140 = t243 * t1961;
                            t2147 = t2024 * t1949;
                            t2152 = -0.32e2 * t1985 * t1940 - 0.256e3 * t2131 * t1940 + 0.64e2 * t2134 * t1940 + 0.256e3 * t2137 * t1940 - 0.64e2 * t2140 * t1940 + 0.32e2 * t2147 * t1940 + 0.64e2 * t1942 * t1965 - 0.64e2 * t1946 * t1965 - 0.64e2 * t1958 * t1965 + 0.64e2 * t1962 * t1965 + 0.32e2 * t2106 * t1965 + 0.64e2 * t2134 * t1965 - 0.64e2 * t2140 * t1965;
                            t2153 = t2114 * t1953;
                            t2156 = t265 * t1945;
                            t2159 = t265 * t1961;
                            t2166 = t355 * t1949;
                            t2185 = 0.2e1 * t2153 * t1940 - 0.8e1 * t2009 * t1965 + 0.64e2 * t2159 * t1965 + 0.32e2 * t2166 * t1965 - 0.32e2 * t1975 * t1970 + 0.32e2 * t2109 * t1970 - 0.8e1 * t1971 * t2001 - 0.32e2 * t1985 * t2001 + 0.8e1 * t2009 * t2001 + 0.32e2 * t2012 * t2001 + 0.32e2 * t2147 * t2001 + 0.2e1 * t2153 * t2001 - 0.64e2 * t2156 * t2001 + 0.64e2 * t2159 * t2001;
                            t2203 = t139 * t2017;
                            t2206 = t145 * t1966;
                            t2215 = 0.64e2 * t1942 * t2001 - 0.64e2 * t1946 * t2001 + 0.256e3 * t1950 * t2001 + 0.64e2 * t2006 * t1965 - 0.64e2 * t2156 * t1965 + 0.2e1 * t2115 * t1970 - 0.32e2 * t2147 * t1970 - 0.2e1 * t2153 * t1970 + 0.32e2 * t1975 * t2001 - 0.32e2 * t2109 * t2001 - 0.2e1 * t2115 * t2001 - 0.256e3 * t2203 * t2001 + 0.256e3 * t2206 * t2001;
                            t2244 = -0.256e3 * t1954 * t2001 - 0.64e2 * t1958 * t2001 + 0.64e2 * t1962 * t2001 + 0.32e2 * t2109 * t1965 + 0.2e1 * t2115 * t1965 - 0.32e2 * t2147 * t1965 - 0.32e2 * t1967 * t1970 + 0.64e2 * t2006 * t1970 - 0.32e2 * t2018 * t1970 - 0.32e2 * t2006 * t2001 - 0.256e3 * t2131 * t2001 + 0.64e2 * t2134 * t2001 + 0.256e3 * t2137 * t2001 - 0.64e2 * t2140 * t2001;
                            t2262 = t1981 * t1949;
                            t2265 = t20 * t225;
                            t2266 = t2265 * t1953;
                            t2277 = 0.64e2 * t1988 * t1940 - 0.64e2 * t1992 * t1940 - 0.64e2 * t2156 * t1940 + 0.64e2 * t2159 * t1940 - 0.2e1 * t2153 * t1965 - 0.64e2 * t1978 * t1970 + 0.64e2 * t1982 * t1970 + 0.32e2 * t1985 * t1970 + 0.64e2 * t1982 * t2001 + 0.64e2 * t1988 * t2001 - 0.64e2 * t1992 * t2001 - 0.256e3 * t2262 * t2001 + 0.256e3 * t2266 * t2001;
                            t2306 = -0.256e3 * t2203 * t1940 + 0.256e3 * t2206 * t1940 - 0.64e2 * t1958 * t1970 + 0.64e2 * t1962 * t1970 + 0.8e1 * t1971 * t1965 - 0.64e2 * t2012 * t1965 + 0.64e2 * t1988 * t1970 - 0.64e2 * t1992 * t1970 - 0.8e1 * t2009 * t1970 - 0.64e2 * t2012 * t1970 + 0.64e2 * t2134 * t1970 - 0.64e2 * t2140 * t1970 - 0.64e2 * t2156 * t1970 + 0.32e2 * t2166 * t1970;
                            t2310 = t1981 * t2017;
                            t2313 = t2265 * t1966;
                            t2336 = -0.64e2 * t1978 * t1940 + 0.64e2 * t1982 * t1940 - 0.256e3 * t2262 * t1940 + 0.256e3 * t2266 * t1940 + 0.256e3 * t2310 * t1940 - 0.256e3 * t2313 * t1940 + 0.64e2 * t2159 * t1970 - 0.64e2 * t1978 * t2001 + 0.256e3 * t2310 * t2001 - 0.256e3 * t2313 * t2001 - 0.64e2 * t266 * t2047 + 0.64e2 * t526 * t2047 - 0.64e2 * t266 * t2087;
                            t2337 = t1981 * t193;
                            t2342 = t190 * t1354;
                            t2343 = t93 * t299;
                            t2344 = t2343 * t346;
                            t2357 = t1981 * t200;
                            t2360 = t2024 * t200;
                            t2371 = 0.64e2 * t194 * t2023 - 0.64e2 * t2337 * t2023 - 0.40e2 * t2025 * t2047 - 0.40e2 * t2025 * t2087 + 0.8e1 * t2344 * t2030 + 0.8e1 * t2032 * t2342 + 0.64e2 * t2337 * t2047 - 0.64e2 * t2357 * t2047 + 0.40e2 * t2360 * t2047 - 0.64e2 * t2337 * t2081 + 0.64e2 * t2337 * t2087 - 0.64e2 * t2357 * t2087 + 0.40e2 * t2360 * t2087 - 0.8e1 * t2344 * t2342;
                            t2384 = t2068 * t160;
                            t2401 = 0.64e2 * t2357 * t2023 - 0.40e2 * t2360 * t2023 + 0.8e1 * t512 * t2023 - 0.64e2 * t520 * t2023 + 0.64e2 * t526 * t2023 - 0.8e1 * t529 * t2023 + 0.40e2 * t2025 * t2081 + 0.64e2 * t2357 * t2081 - 0.40e2 * t2360 * t2081 - 0.64e2 * t266 * t2081 + 0.64e2 * t526 * t2081 + 0.64e2 * t526 * t2087 + 0.2e1 * t2384 * t2342;
                            t2402 = t93 * c4;
                            t2403 = t2402 * t193;
                            t2406 = t2402 * t200;
                            t2409 = t1566 * c4;
                            t2410 = t2409 * t193;
                            t2413 = t2409 * t200;
                            t2436 = -0.32e2 * t2403 * t2023 + 0.32e2 * t2406 * t2023 - 0.4e1 * t2410 * t2023 + 0.4e1 * t2413 * t2023 + 0.32e2 * t2403 * t2047 - 0.32e2 * t2406 * t2047 + 0.4e1 * t2410 * t2047 - 0.4e1 * t2413 * t2047 - 0.32e2 * t2403 * t2081 + 0.32e2 * t2406 * t2081 - 0.4e1 * t2410 * t2081 + 0.4e1 * t2413 * t2081 + 0.32e2 * t2403 * t2087 - 0.32e2 * t2406 * t2087;
                            t2448 = t2072 * t346;
                            t2451 = t2076 * t136;
                            t2454 = t91 * t2041;
                            t2455 = t2454 * t155;
                            t2468 = 0.64e2 * t292 * t2023 - 0.2e1 * t2384 * t2030 + 0.8e1 * t2448 * t2030 + 0.2e1 * t2451 * t2030 - 0.8e1 * t2455 * t2030 - 0.8e1 * t529 * t2081 + 0.4e1 * t2410 * t2087 - 0.4e1 * t2413 * t2087 - 0.64e2 * t244 * t2087 + 0.64e2 * t292 * t2087 - 0.8e1 * t2448 * t2342 - 0.2e1 * t2451 * t2342 + 0.8e1 * t2455 * t2342;
                            t2479 = t2454 * t184;
                            t2490 = t2343 * t175;
                            t2493 = t2031 * t184;
                            t2500 = -0.64e2 * t194 * t2047 - 0.64e2 * t244 * t2023 + 0.2e1 * t2036 * t2342 - 0.2e1 * t2043 * t2342 + 0.64e2 * t520 * t2047 + 0.8e1 * t2479 * t2056 - 0.8e1 * t2490 * t2056 + 0.8e1 * t2493 * t2056 - 0.2e1 * t2069 * t2063 + 0.8e1 * t2073 * t2063 + 0.2e1 * t2077 * t2063 - 0.8e1 * t2479 * t2063 + 0.8e1 * t2490 * t2063 - 0.8e1 * t2493 * t2063;
                            pr = t1991 + t2039 + t2084 + t2118 + t2152 + t2185 + t2215 + t2244 + t2277 + t2306 + t2336 + t2371 + t2401 + t2436 + t2468 + t2500;
                            pLocal(i1,i2,i3) = pr;
                        }
                    }
                }

        }
        else
        {
            OV_ABORT("error");
        }

    // else if( iswCase==1 ) // ** OLD ***
    // {
    //   // files generated from AMP/ism/maple/DiskAnnulusSolutionWDH.mw
    //   C4 = 1./10.;  // Make |u|=.5, |p|=5
    //   #Include "incompressible/ismAnnulusDisplacementDisplacementWDH.h"
    //   assignIncompressibleAnnulusSolution();
    // }
    // else if( iswCase==2 )
    // {
    //   #Include "incompressible/ismAnnulusTractionDisplacementWDH.h"
    //   assignIncompressibleAnnulusSolution();
    // }
    // else
    // {
    //   C4 = 1./4.; // make u size 1
    //   #Include "incompressible/ismAnnulusTractionTractionWDH.h"
    //   assignIncompressibleAnnulusSolution();
    // }

    }
    else if( userKnownSolution=="incompressibleDiskSolution" )
    {
    // --- incompressible disk solution -----
        const int pc = dbase.get<int >("pc");
        assert( pc>=0 );

        MappedGrid & mg = cg[grid];
        mg.update(MappedGrid::THEvertex | MappedGrid::THEcenter | MappedGrid::THEmask );
        OV_GET_SERIAL_ARRAY_CONST(int,mg.mask(),maskLocal);
        OV_GET_SERIAL_ARRAY_CONST(real,mg.vertex(),vertex);
        RealArray & u = ua;

        const int & iswCase = db.get<int>("iswCase");
        const int n = ipar[0];
        const int m = ipar[1];     

        if( t<3*dt )
            printF("UDKS: evaluate the incompressible disk solution: iswCase=%i, n=%d, m=%d, t=%9.3e ( Jn(lambda_nm r) )\n",iswCase,n,m,t);

        const Real eps=1.e-9;   // ********************************** FIX ME *********

    // diskScaling : scale annulus solution to be size one
        if( !db.has_key("diskScaling") )
        {
            db.put<Real>("diskScaling") = -1.;  // this means scaling has not be set 
        }
        Real & diskScaling = db.get<Real>("diskScaling");
        Real uNorm; 

    // Real C4 =1.; // defines the scale
        if( iswCase==1 )
        {
      // ---- Disks : displacement ----
      // printF("UDKS: EVAL NEW DISK SOLUTION displacement-displacement\n");

            Real c3=1./4.;   // old scaling 

// --- roots of the frequency equation -----
//   File written by diskAnnulusSolution.mw on 2021-10-18
//   Solution is of the form Jn(omega*r) with roots omega_{n,m} 
Real omegaArray[]={
    // n=1, m=1,2,..,4
        5.1356223018406826e+00, 8.4172441403998649e+00, 1.1619841172149059e+01, 1.4795951782351261e+01,
    // n=2, m=1,2,..,4
        6.3801618959239835e+00, 9.7610231299816697e+00, 1.3015200721698434e+01, 1.6223466160318768e+01,
    // n=3, m=1,2,..,4
        7.5883424345038044e+00, 1.1064709488501185e+01, 1.4372536671617590e+01, 1.7615966049804833e+01
};
#define omegaRoot(n,m) omegaArray[(m-1)+4*(n-1)]
// ------------------------------------------------------------------------------------------------ 
//    caseName: Exact solution for incompressible linear elasticity 
// Solid disk of radius R=1; r=0 bounded solution; r=1 Displacement BC. 
// File Written by diskAnnulusSolution.mw (maple) on 2021-10-18 
// Set c3 to scale disk solution and c4 to scale the annulus solution 
// Set n and m to choose the solution, Jn(omega*r), omega=omega(n,m) 
// The eigenvalues are save in a separate include file. 
// ------------------------------------------------------------------------------------------------ 
Real rho,c,mu,omega; 
rho=1.;
mu=1.;
omega=omegaRoot(n,m);
Real t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49
        ,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,t86,t87,t88,t89,t90,t91,t92,t93,t94,t95,t96,t97,t98,t99
        ,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t130,t131,t132,t133,t134,t135,t136,t137,t138,t139,t140,t141,t142,t143,t144,t145,t146,t147,t148,t149
        ,t150,t151,t152,t153,t154,t155,t156,t157,t158,t159,t160,t161,t162,t163,t164,t165,t166,t167,t168,t169,t170,t171,t172,t173,t174,t175,t176,t177,t178,t179,t180,t181,t182,t183,t184,t185,t186,t187,t188,t189,t190,t191,t192,t193,t194,t195,t196,t197,t198,t199
        ,t200;



            if( diskScaling <= 0. )
            { // -- first time: compute the norm of the solution so we can scale to be size one
                Real t0=0.; 
                    Real ur,vr,pr; 
                    int i1,i2,i3;
                    Real theta; 
                    const Real eps=1.e-9;   // ********************************** FIX ME *********
                    if( t0<= dt )
                        printF("\n *** WARNING *** need a better fix for disk solution at r=0 ******\n\n");
                    Real u1Max=0., u2Max=0., u3Max=0., pMax=0.; 
                    if( numberOfDimensions==2 )
                    {
                        FOR_3D(i1,i2,i3,I1,I2,I3)
                        {
              // Reference coordinates:
                            Real x = vertex(i1,i2,i3,0);
                            Real y = vertex(i1,i2,i3,1);
                            Real r = sqrt( x*x + y*y );
                            if( r < eps ) //  avoid division by zero
                            {
                                x=eps;           // avoid atan(0,0)
                                y=eps;
                                r=sqrt(2.)*eps;  // r=sqrt(x**2+y**2)
                            }   
                            theta = atan2(y,x); 
                            Real cosTheta = cos(theta);
                            Real sinTheta = sin(theta);
                            t1 = omega * r;
                            t2 = jn(n, t1);
                            t3 = t2 * c3;
                            t4 = 0.1e1 / r;
                            t5 = omega * t0;
                            t6 = cos(t5);
                            t7 = t6 * t4;
                            t8 = n * theta;
                            t9 = cos(t8);
                            t10 = t9 * t7;
                            t12 = sin(t5);
                            t13 = t12 * t4;
                            t14 = sin(t8);
                            t15 = t14 * t13;
                            t17 = jn(n, omega);
                            t18 = t17 * c3;
                            t19 = pow(r, n);
                            t20 = t19 * t18;
                            ur = -t10 * t20 + t10 * t3 - t15 * t20 + t15 * t3;
                            t23 = 0.1e1 / n;
                            t26 = jn(n + 0.1e1, t1);
                            t27 = t26 * t23 * omega;
                            t34 = t14 * t7;
                            t36 = t9 * t13;
                            vr = -t9 * t12 * c3 * t27 + t14 * t6 * c3 * t27 + t34 * t20 - t36 * t20 - t34 * t3 + t36 * t3;
                            t40 = t23 * t18;
                            t41 = omega * omega;
                            t42 = t19 * t41;
                            pr = -t14 * t12 * t42 * t40 - t9 * t6 * t42 * t40;
              // (u1,u2) = ur*rHat + vr*thetaHat 
                            u(i1,i2,i3,u1c) = ur*cosTheta - vr*sinTheta;
                            u(i1,i2,i3,u2c) = ur*sinTheta + vr*cosTheta;
                            u(i1,i2,i3,pc ) = pr;
                            if( maskLocal(i1,i2,i3) > 0 )
                            {
                                u1Max=max(u1Max,u(i1,i2,i3,u1c));
                                u2Max=max(u2Max,u(i1,i2,i3,u2c));
                                pMax =max(pMax, u(i1,i2,i3,pc));
                            }
                        }
                        uNorm = max(u1Max,u2Max,u3Max,pMax);
                    }
                    else
                    {
            // ** FOR NOW USE THE 2D SOLUTION ***
                        OV_ABORT("UDKS:ERROR: disk solution : 3D called");
                        FOR_3D(i1,i2,i3,I1,I2,I3)
                        {
                            Real x = vertex(i1,i2,i3,0);
                            Real y = vertex(i1,i2,i3,1);
                            Real z = vertex(i1,i2,i3,2);
                            Real r = sqrt( x*x + y*y );
                            if( r < eps ) //  avoid division by zero
                            {
                                x=eps;           // avoid atan(0,0)
                                y=eps;
                                r=sqrt(2.)*eps;  // r=sqrt(x**2+y**2)
                            }   
                            theta = atan2(y,x); 
                            Real cosTheta = cos(theta);
                            Real sinTheta = sin(theta);
                            t1 = omega * r;
                            t2 = jn(n, t1);
                            t3 = t2 * c3;
                            t4 = 0.1e1 / r;
                            t5 = omega * t0;
                            t6 = cos(t5);
                            t7 = t6 * t4;
                            t8 = n * theta;
                            t9 = cos(t8);
                            t10 = t9 * t7;
                            t12 = sin(t5);
                            t13 = t12 * t4;
                            t14 = sin(t8);
                            t15 = t14 * t13;
                            t17 = jn(n, omega);
                            t18 = t17 * c3;
                            t19 = pow(r, n);
                            t20 = t19 * t18;
                            ur = -t10 * t20 + t10 * t3 - t15 * t20 + t15 * t3;
                            t23 = 0.1e1 / n;
                            t26 = jn(n + 0.1e1, t1);
                            t27 = t26 * t23 * omega;
                            t34 = t14 * t7;
                            t36 = t9 * t13;
                            vr = -t9 * t12 * c3 * t27 + t14 * t6 * c3 * t27 + t34 * t20 - t36 * t20 - t34 * t3 + t36 * t3;
                            t40 = t23 * t18;
                            t41 = omega * omega;
                            t42 = t19 * t41;
                            pr = -t14 * t12 * t42 * t40 - t9 * t6 * t42 * t40;
              // (u1,u2) = ur*rHat + vr*thetaHat 
                            u(i1,i2,i3,u1c) = ur*cosTheta - vr*sinTheta;
                            u(i1,i2,i3,u2c) = ur*sinTheta + vr*cosTheta;
                            u(i1,i2,i3,u3c) = 0.;
                            u(i1,i2,i3,pc ) = pr; 
                        }   
                    }
                if( uNorm==0. )
                {
                    printF("UDKS: ERROR: disk solution is zero! uNorm=%e\n",uNorm);
                    OV_ABORT("ERROR");
                }        
                diskScaling = c3/uNorm;
                printF("\n +++++ UDKS: uNorm=%9.2e +++++",uNorm);
            }
            c3 = diskScaling;       

                Real ur,vr,pr; 
                int i1,i2,i3;
                Real theta; 
                const Real eps=1.e-9;   // ********************************** FIX ME *********
                if( t<= dt )
                    printF("\n *** WARNING *** need a better fix for disk solution at r=0 ******\n\n");
                Real u1Max=0., u2Max=0., u3Max=0., pMax=0.; 
                if( numberOfDimensions==2 )
                {
                    FOR_3D(i1,i2,i3,I1,I2,I3)
                    {
            // Reference coordinates:
                        Real x = vertex(i1,i2,i3,0);
                        Real y = vertex(i1,i2,i3,1);
                        Real r = sqrt( x*x + y*y );
                        if( r < eps ) //  avoid division by zero
                        {
                            x=eps;           // avoid atan(0,0)
                            y=eps;
                            r=sqrt(2.)*eps;  // r=sqrt(x**2+y**2)
                        }   
                        theta = atan2(y,x); 
                        Real cosTheta = cos(theta);
                        Real sinTheta = sin(theta);
                        t1 = omega * r;
                        t2 = jn(n, t1);
                        t3 = t2 * c3;
                        t4 = 0.1e1 / r;
                        t5 = omega * t;
                        t6 = cos(t5);
                        t7 = t6 * t4;
                        t8 = n * theta;
                        t9 = cos(t8);
                        t10 = t9 * t7;
                        t12 = sin(t5);
                        t13 = t12 * t4;
                        t14 = sin(t8);
                        t15 = t14 * t13;
                        t17 = jn(n, omega);
                        t18 = t17 * c3;
                        t19 = pow(r, n);
                        t20 = t19 * t18;
                        ur = -t10 * t20 + t10 * t3 - t15 * t20 + t15 * t3;
                        t23 = 0.1e1 / n;
                        t26 = jn(n + 0.1e1, t1);
                        t27 = t26 * t23 * omega;
                        t34 = t14 * t7;
                        t36 = t9 * t13;
                        vr = -t9 * t12 * c3 * t27 + t14 * t6 * c3 * t27 + t34 * t20 - t36 * t20 - t34 * t3 + t36 * t3;
                        t40 = t23 * t18;
                        t41 = omega * omega;
                        t42 = t19 * t41;
                        pr = -t14 * t12 * t42 * t40 - t9 * t6 * t42 * t40;
            // (u1,u2) = ur*rHat + vr*thetaHat 
                        u(i1,i2,i3,u1c) = ur*cosTheta - vr*sinTheta;
                        u(i1,i2,i3,u2c) = ur*sinTheta + vr*cosTheta;
                        u(i1,i2,i3,pc ) = pr;
                        if( maskLocal(i1,i2,i3) > 0 )
                        {
                            u1Max=max(u1Max,u(i1,i2,i3,u1c));
                            u2Max=max(u2Max,u(i1,i2,i3,u2c));
                            pMax =max(pMax, u(i1,i2,i3,pc));
                        }
                    }
                    uNorm = max(u1Max,u2Max,u3Max,pMax);
                }
                else
                {
          // ** FOR NOW USE THE 2D SOLUTION ***
                    OV_ABORT("UDKS:ERROR: disk solution : 3D called");
                    FOR_3D(i1,i2,i3,I1,I2,I3)
                    {
                        Real x = vertex(i1,i2,i3,0);
                        Real y = vertex(i1,i2,i3,1);
                        Real z = vertex(i1,i2,i3,2);
                        Real r = sqrt( x*x + y*y );
                        if( r < eps ) //  avoid division by zero
                        {
                            x=eps;           // avoid atan(0,0)
                            y=eps;
                            r=sqrt(2.)*eps;  // r=sqrt(x**2+y**2)
                        }   
                        theta = atan2(y,x); 
                        Real cosTheta = cos(theta);
                        Real sinTheta = sin(theta);
                        t1 = omega * r;
                        t2 = jn(n, t1);
                        t3 = t2 * c3;
                        t4 = 0.1e1 / r;
                        t5 = omega * t;
                        t6 = cos(t5);
                        t7 = t6 * t4;
                        t8 = n * theta;
                        t9 = cos(t8);
                        t10 = t9 * t7;
                        t12 = sin(t5);
                        t13 = t12 * t4;
                        t14 = sin(t8);
                        t15 = t14 * t13;
                        t17 = jn(n, omega);
                        t18 = t17 * c3;
                        t19 = pow(r, n);
                        t20 = t19 * t18;
                        ur = -t10 * t20 + t10 * t3 - t15 * t20 + t15 * t3;
                        t23 = 0.1e1 / n;
                        t26 = jn(n + 0.1e1, t1);
                        t27 = t26 * t23 * omega;
                        t34 = t14 * t7;
                        t36 = t9 * t13;
                        vr = -t9 * t12 * c3 * t27 + t14 * t6 * c3 * t27 + t34 * t20 - t36 * t20 - t34 * t3 + t36 * t3;
                        t40 = t23 * t18;
                        t41 = omega * omega;
                        t42 = t19 * t41;
                        pr = -t14 * t12 * t42 * t40 - t9 * t6 * t42 * t40;
            // (u1,u2) = ur*rHat + vr*thetaHat 
                        u(i1,i2,i3,u1c) = ur*cosTheta - vr*sinTheta;
                        u(i1,i2,i3,u2c) = ur*sinTheta + vr*cosTheta;
                        u(i1,i2,i3,u3c) = 0.;
                        u(i1,i2,i3,pc ) = pr; 
                    }   
                }

                if( computeVelocity )
                {
                    printF("\n $$$$$$$$$$ userDefinedKnownSolution: EVALUATE VELOCITY, ACCEL and PRESSURE AT PAST TIMES t=%9.3e dt=%16.6e $$$$$$$\n\n",t,dt);
                    assert( dt>0. );
                    assert( vgf!=NULL );
                    const int currentVelocity = dbase.get<int>("currentVelocity");
                    assert( currentVelocity>=0 && currentVelocity<numberOfVelocityFunctions );
          // fill in t0=0 and some past time levels
                    for( int level=0; level<numberOfVelocityFunctions; level++ )
                    {
                        real tc = 0. - level*dt;
                        const int prevVelocity = (currentVelocity - level + numberOfVelocityFunctions) % numberOfVelocityFunctions;
                        OV_GET_SERIAL_ARRAY(real,vgf[prevVelocity][grid],vLocal);
                        OV_GET_SERIAL_ARRAY(real,vtgf[prevVelocity][grid],vtLocal);
                        FOR_3D(i1,i2,i3,I1,I2,I3)
                        {
                                Real x = vertex(i1,i2,i3,0);
                                Real y = vertex(i1,i2,i3,1);
                                Real r = sqrt( x*x + y*y );
                                if( r < eps ) //  avoid division by zero
                                {
                                    x=eps;           // avoid atan(0,0)
                                    y=eps;
                                    r=sqrt(2.)*eps;  // r=sqrt(x**2+y**2)
                                }   
                                Real theta = atan2(y,x); 
                                Real cosTheta = cos(theta);
                                Real sinTheta = sin(theta);
              // Real x = vertex(i1,i2,i3,0);
              // Real y = vertex(i1,i2,i3,1);
              // Real r = sqrt( x*x + y*y );
              // Real theta = atan2(y,x); 
              // Real cosTheta = x/r;
              // Real sinTheta = y/r;
                            t2 = jn(n, omega);
                            t3 = pow(r, n);
                            t5 = omega * r;
                            t6 = jn(n, t5);
                            t7 = -t3 * t2 + t6;
                            t8 = omega * tc;
                            t9 = cos(t8);
                            t10 = n * theta;
                            t11 = sin(t10);
                            t13 = sin(t8);
                            t14 = cos(t10);
                            t18 = 0.1e1 / r;
                            ur = t18 * (t11 * t9 - t14 * t13) * t7 * c3 * omega;
                            t26 = jn(n + 0.1e1, t5);
                            vr = t18 / n * (-r * omega * t26 + t7 * n) * omega * c3 * (t11 * t13 + t14 * t9);
              // (u1,u2) = ur*rHat + vr*thetaHat 
                            vLocal(i1,i2,i3,u1c) = ur*cosTheta - vr*sinTheta;
                            vLocal(i1,i2,i3,u2c) = ur*sinTheta + vr*cosTheta;
                            t1 = omega * omega;
                            t3 = jn(n, omega);
                            t4 = pow(r, n);
                            t6 = omega * r;
                            t7 = jn(n, t6);
                            t8 = -t4 * t3 + t7;
                            t9 = omega * tc;
                            t10 = sin(t9);
                            t11 = n * theta;
                            t12 = sin(t11);
                            t14 = cos(t9);
                            t15 = cos(t11);
                            t19 = 0.1e1 / r;
                            ur = -t19 * (t12 * t10 + t15 * t14) * t8 * t1 * c3;
                            t28 = jn(n + 0.1e1, t6);
                            vr = t19 / n * (-r * omega * t28 + t8 * n) * t1 * (-t15 * t10 + t12 * t14) * c3;
              // (u1,u2) = ur*rHat + vr*thetaHat 
                            vtLocal(i1,i2,i3,u1c) = ur*cosTheta - vr*sinTheta;
                            vtLocal(i1,i2,i3,u2c) = ur*sinTheta + vr*cosTheta;        
                        }      
                    }
          // --- save pressure for time extrapolation ---
                    const int currentPressure = dbase.get<int>("currentPressure");
                    assert( currentPressure>=0 && currentPressure<max(1,numberOfPressureFunctions) );
                    for( int level=0; level<numberOfPressureFunctions; level++ )
                    {
                        real tc = 0. - (level+3)*dt;  // -- save pressure starting at t-3*dt 
                        const int prevPressure = (currentPressure - level + numberOfPressureFunctions) % numberOfPressureFunctions;
                        assert( pgf !=NULL );
                        OV_GET_SERIAL_ARRAY(real,pgf[prevPressure][grid],pLocal);
                        FOR_3D(i1,i2,i3,I1,I2,I3)
                        {
                                Real x = vertex(i1,i2,i3,0);
                                Real y = vertex(i1,i2,i3,1);
                                Real r = sqrt( x*x + y*y );
                                if( r < eps ) //  avoid division by zero
                                {
                                    x=eps;           // avoid atan(0,0)
                                    y=eps;
                                    r=sqrt(2.)*eps;  // r=sqrt(x**2+y**2)
                                }   
                                Real theta = atan2(y,x); 
                                Real cosTheta = cos(theta);
                                Real sinTheta = sin(theta);
              // Real x = vertex(i1,i2,i3,0);
              // Real y = vertex(i1,i2,i3,1);
              // Real r = sqrt( x*x + y*y );
              // Real theta = atan2(y,x); 
              // Real cosTheta = x/r;
              // Real sinTheta = y/r;
                            t1 = omega * r;
                            t2 = jn(n, t1);
                            t3 = t2 * c3;
                            t4 = 0.1e1 / r;
                            t5 = omega * tc;
                            t6 = cos(t5);
                            t7 = t6 * t4;
                            t8 = n * theta;
                            t9 = cos(t8);
                            t10 = t9 * t7;
                            t12 = sin(t5);
                            t13 = t12 * t4;
                            t14 = sin(t8);
                            t15 = t14 * t13;
                            t17 = jn(n, omega);
                            t18 = t17 * c3;
                            t19 = pow(r, n);
                            t20 = t19 * t18;
                            ur = -t10 * t20 + t10 * t3 - t15 * t20 + t15 * t3;
                            t23 = 0.1e1 / n;
                            t26 = jn(n + 0.1e1, t1);
                            t27 = t26 * t23 * omega;
                            t34 = t14 * t7;
                            t36 = t9 * t13;
                            vr = -t9 * t12 * c3 * t27 + t14 * t6 * c3 * t27 + t34 * t20 - t36 * t20 - t34 * t3 + t36 * t3;
                            t40 = t23 * t18;
                            t41 = omega * omega;
                            t42 = t19 * t41;
                            pr = -t14 * t12 * t42 * t40 - t9 * t6 * t42 * t40;
                            pLocal(i1,i2,i3) = pr;
                        }
                    }
                }

        }
        else if( iswCase==2 )
        {
      // ---- Disk : traction ----
      // printF("UDKS: EVAL NEW DISK SOLUTION\n");

            Real c3=1./4.; 
            if( diskScaling>0 )
                c3 = diskScaling;

// --- roots of the frequency equation -----
//   File written by diskAnnulusSolution.mw on 2022-06-29
//   Solution is of the form Jn(omega*r) with roots omega_{n,m} 
Real omegaArray[]={
    // n=1, m=1,2,..,4
        3.0542369282271403e+00, 6.7061331941584591e+00, 9.9694678230875958e+00, 1.3170370856016123e+01,
    // n=2, m=1,2,..,4
        2.3538876347708396e+00, 4.7844442796765740e+00, 8.1766702433977513e+00, 1.1446029444951778e+01,
    // n=3, m=1,2,..,4
        3.6471263298525899e+00, 6.4781057915316590e+00, 9.6197905984024880e+00, 1.2882258435475911e+01,
    // n=4, m=1,2,..,4
        4.7880010353864664e+00, 8.0927304341395894e+00, 1.1059090301706138e+01, 1.4295359989629739e+01,
    // n=5, m=1,2,..,4
        5.8638944374942023e+00, 9.6143023902817394e+00, 1.2500478612658155e+01, 1.5695307037022381e+01,
    // n=6, m=1,2,..,4
        6.9057986668177223e+00, 1.1048755843098352e+01, 1.3939364538190292e+01, 1.7087954497231869e+01,
    // n=7, m=1,2,..,4
        7.9274893997322877e+00, 1.2410259297341431e+01, 1.5366324082026815e+01, 1.8476125448055522e+01,
    // n=8, m=1,2,..,4
        8.9360328180040246e+00, 1.3714343502448062e+01, 1.6772215029002498e+01, 1.9860168825158898e+01,
    // n=9, m=1,2,..,4
        9.9354478214627102e+00, 1.4974496337529492e+01, 1.8151048232718931e+01, 2.1238551161238865e+01,
    // n=10, m=1,2,..,4
        1.0928204807305053e+01, 1.6201207065273404e+01, 1.9500486215674639e+01, 2.2608628876597680e+01,
    // n=11, m=1,2,..,4
        1.1915916481800491e+01, 1.7402229137402822e+01, 2.0821022816609477e+01, 2.3967480401389517e+01,
    // n=12, m=1,2,..,4
        1.2899687104309149e+01, 1.8583196222563467e+01, 2.2114840204648401e+01, 2.5312585848839412e+01,
    // n=13, m=1,2,..,4
        1.3880302573704197e+01, 1.9748207857050163e+01, 2.3384859803037886e+01, 2.6642217549474359e+01,
    // n=14, m=1,2,..,4
        1.4858340299800952e+01, 2.0900280920341116e+01, 2.4634126054088785e+01, 2.7955539069249390e+01,
    // n=15, m=1,2,..,4
        1.5834235953539674e+01, 2.2041670960018678e+01, 2.5865477904105144e+01, 2.9252501050943966e+01,
    // n=16, m=1,2,..,4
        1.6808325753901503e+01, 2.3174094105061112e+01, 2.7081412738174955e+01, 3.0533642030186897e+01,
    // n=17, m=1,2,..,4
        1.7780874229213456e+01, 2.4298879285620925e+01, 2.8284058998729757e+01, 3.1799876859002021e+01,
    // n=18, m=1,2,..,4
        1.8752093010143026e+01, 2.5417073161571597e+01, 2.9475201161328325e+01, 3.3052317881946933e+01,
    // n=19, m=1,2,..,4
        1.9722153895251924e+01, 2.6529513241559757e+01, 3.0656324658780802e+01, 3.4292143950833705e+01,
    // n=20, m=1,2,..,4
        2.0691198149296995e+01, 2.7636879539908085e+01, 3.1828664268005112e+01, 3.5520514895283197e+01,
    // n=21, m=1,2,..,4
        2.1659343258396381e+01, 2.8739731619062381e+01, 3.2993248652456406e+01, 3.6738521880361787e+01,
    // n=22, m=1,2,..,4
        2.2626687928423105e+01, 2.9838535557474345e+01, 3.4150938528933265e+01, 3.7947163006349841e+01,
    // n=23, m=1,2,..,4
        2.3593315844640883e+01, 3.0933683876694453e+01, 3.5302458182458539e+01, 3.9147335155124392e+01,
    // n=24, m=1,2,..,4
        2.4559298541591919e+01, 3.2025510477614650e+01, 3.6448420984307281e+01, 4.0339835453424609e+01,
    // n=25, m=1,2,..,4
        2.5524697623193163e+01, 3.3114301988674151e+01, 3.7589349854618657e+01, 4.1525367878948303e+01,
    // n=26, m=1,2,..,4
        2.6489566501071567e+01, 3.4200306498858498e+01, 3.8725693609341911e+01, 4.2704552187144726e+01,
    // n=27, m=1,2,..,4
        2.7453951770778919e+01, 3.5283740359236365e+01, 3.9857840018530181e+01, 4.3877933485616522e+01,
    // n=28, m=1,2,..,4
        2.8417894312377217e+01, 3.6364793539971769e+01, 4.0986126261663405e+01, 4.5045991531969354e+01,
    // n=29, m=1,2,..,4
        2.9381430178797795e+01, 3.7443633894047535e+01, 4.2110847330860441e+01, 4.6209149294950587e+01,
    // n=30, m=1,2,..,4
        3.0344591319051433e+01, 3.8520410584177811e+01, 4.3232262817018217e+01, 4.7367780593101966e+01,
    // n=31, m=1,2,..,4
        3.1307406171659786e+01, 3.9595256862398819e+01, 4.4350602419409025e+01, 4.8522216779181885e+01,
    // n=32, m=1,2,..,4
        3.2269900155174526e+01, 4.0668292343901846e+01, 4.5466070444268971e+01, 4.9672752519096614e+01,
    // n=33, m=1,2,..,4
        3.3232096076399077e+01, 4.1739624881991048e+01, 4.6578849499299755e+01, 5.0819650751183382e+01,
    // n=34, m=1,2,..,4
        3.4194014472280698e+01, 4.2809352125674733e+01, 4.7689103545569914e+01, 5.1963146924194644e+01,
    // n=35, m=1,2,..,4
        3.5155673897950254e+01, 4.3877562822639638e+01, 4.8796980433204053e+01, 5.3103452611609671e+01,
    // n=36, m=1,2,..,4
        3.6117091170739943e+01, 4.4944337916350388e+01, 4.9902614020154451e+01, 5.4240758592641797e+01,
    // n=37, m=1,2,..,4
        3.7078281577983518e+01, 4.6009751475457672e+01, 5.1006125952406826e+01, 5.5375237480341195e+01,
    // n=38, m=1,2,..,4
        3.8039259054840106e+01, 4.7073871485667912e+01, 5.2107627167742237e+01, 5.6507045966612106e+01,
    // n=39, m=1,2,..,4
        3.9000036337166438e+01, 4.8136760528066574e+01, 5.3207219172556028e+01, 5.7636326743851346e+01,
    // n=40, m=1,2,..,4
        3.9960625093508963e+01, 4.9198476363122920e+01, 5.4304995131380034e+01, 5.8763210153769243e+01,
    // n=41, m=1,2,..,4
        4.0921036039534796e+01, 5.0259072435890889e+01, 5.5401040801025409e+01, 5.9887815605947101e+01,
    // n=42, m=1,2,..,4
        4.1881279037622434e+01, 5.1318598315005633e+01, 5.6495435335173524e+01, 6.1010252801816337e+01,
    // n=43, m=1,2,..,4
        4.2841363183854979e+01, 5.2377100075770668e+01, 5.7588251980420554e+01, 6.2130622793928676e+01,
    // n=44, m=1,2,..,4
        4.3801296884273899e+01, 5.3434620635796429e+01, 5.8679558680945158e+01, 6.3249018905504216e+01,
    // n=45, m=1,2,..,4
        4.4761087921940112e+01, 5.4491200050182180e+01, 5.9769418605901456e+01, 6.4365527531167069e+01,
    // n=46, m=1,2,..,4
        4.5720743516096001e+01, 5.5546875772049683e+01, 6.0857890611175446e+01, 6.5480228836384551e+01,
    // n=47, m=1,2,..,4
        4.6680270374515016e+01, 5.6601682883277986e+01, 6.1945029645153873e+01, 6.6593197370305636e+01,
    // n=48, m=1,2,..,4
        4.7639674739955421e+01, 5.7655654299507294e+01, 6.3030887106541324e+01, 6.7704502604351959e+01,
    // n=49, m=1,2,..,4
        4.8598962431494422e+01, 5.8708820952839823e+01, 6.4115511160946924e+01, 6.8814209406968355e+01,
    // n=50, m=1,2,..,4
        4.9558138881402518e+01, 5.9761211955138743e+01, 6.5198947021886237e+01, 6.9922378463321333e+01,
    // n=51, m=1,2,..,4
        5.0517209168121038e+01, 6.0812854744390541e+01, 6.6281237200959694e+01, 7.1029066647385868e+01
};
#define omegaRoot(n,m) omegaArray[(m-1)+4*(n-1)]
// ------------------------------------------------------------------------------------------------ 
//    caseName: Exact solution for incompressible linear elasticity 
// Solid disk of radius R=1; r=0 bounded solution; r=1 Traction BC. 
// File Written by diskAnnulusSolution.mw (maple) on 2021-10-18 
// Set c3 to scale disk solution and c4 to scale the annulus solution 
// Set n and m to choose the solution, Jn(omega*r), omega=omega(n,m) 
// The eigenvalues are save in a separate include file. 
// ------------------------------------------------------------------------------------------------ 
Real rho,c,mu,omega; 
rho=1.;
mu=1.;
omega=omegaRoot(n,m);
Real t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49
        ,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,t86,t87,t88,t89,t90,t91,t92,t93,t94,t95,t96,t97,t98,t99
        ,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t130,t131,t132,t133,t134,t135,t136,t137,t138,t139,t140,t141,t142,t143,t144,t145,t146,t147,t148,t149
        ,t150,t151,t152,t153,t154,t155,t156,t157,t158,t159,t160,t161,t162,t163,t164,t165,t166,t167,t168,t169,t170,t171,t172,t173,t174,t175,t176,t177,t178,t179,t180,t181,t182,t183,t184,t185,t186,t187,t188,t189,t190,t191,t192,t193,t194,t195,t196,t197,t198,t199
        ,t200;



            if( diskScaling <= 0. )
            { 
        // -- first time: compute the norm of the solution so we can scale to be size one
        // ** WE NEED TO LOOP OVER ALL GRIDS **
                Real t0=0.;
                uNorm=0.; 
                for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                {
                    MappedGrid & mg = cg[grid];
                    mg.update(MappedGrid::THEvertex | MappedGrid::THEcenter);
                    OV_GET_SERIAL_ARRAY_CONST(real,mg.vertex(),vertex);
                    Index I1,I2,I3;
                    getIndex(mg.gridIndexRange(),I1,I2,I3);
                    const int numberOfComponents=3; // u1,u2,p
                    RealArray u(I1,I2,I3,numberOfComponents);    // temp array

                    Real gridNorm=0.;
                        Real ur,vr,pr; 
                        int i1,i2,i3;
                        Real theta; 
                        const Real eps=1.e-9;   // ********************************** FIX ME *********
                        if( t0<= dt )
                            printF("\n *** WARNING *** need a better fix for disk solution at r=0 ******\n\n");
                        Real u1Max=0., u2Max=0., u3Max=0., pMax=0.; 
                        if( numberOfDimensions==2 )
                        {
                            FOR_3D(i1,i2,i3,I1,I2,I3)
                            {
                // Reference coordinates:
                                Real x = vertex(i1,i2,i3,0);
                                Real y = vertex(i1,i2,i3,1);
                                Real r = sqrt( x*x + y*y );
                                if( r < eps ) //  avoid division by zero
                                {
                                    x=eps;           // avoid atan(0,0)
                                    y=eps;
                                    r=sqrt(2.)*eps;  // r=sqrt(x**2+y**2)
                                }   
                                theta = atan2(y,x); 
                                Real cosTheta = cos(theta);
                                Real sinTheta = sin(theta);
                                t1 = omega * r;
                                t2 = jn(n, t1);
                                t3 = t2 * c3;
                                t4 = 0.1e1 / r;
                                t5 = omega * t0;
                                t6 = cos(t5);
                                t7 = t6 * t4;
                                t8 = n * theta;
                                t9 = cos(t8);
                                t12 = sin(t5);
                                t13 = t12 * t4;
                                t14 = sin(t8);
                                t17 = n * n;
                                t19 = omega * omega;
                                t23 = 0.1e1 / (0.2e1 * t17 - t19 - 0.2e1 * n) * c3;
                                t24 = pow(r, n);
                                t26 = t4 * t24 * t23;
                                t27 = n + 0.1e1;
                                t28 = jn(t27, omega);
                                t29 = t28 * n;
                                t30 = t6 * omega;
                                t35 = t12 * omega;
                                t40 = jn(n, omega);
                                t41 = t40 * t17;
                                t42 = t9 * t6;
                                t46 = t14 * t12;
                                t50 = t40 * n;
                                t51 = t42 * t50;
                                t54 = t46 * t50;
                                ur = 0.2e1 * t14 * t35 * t29 * t26 + 0.2e1 * t9 * t30 * t29 * t26 + t14 * t13 * t3 - 0.2e1 * t42 * t41 * t26 - 0.2e1 * t46 * t41 * t26 + t9 * t7 * t3 + 0.2e1 * t51 * t26 + 0.2e1 * t54 * t26;
                                t57 = t14 * t6;
                                t63 = jn(t27, t1);
                                t64 = t63 / n * omega;
                                t75 = t9 * t12;
                                vr = -t9 * t12 * c3 * t64 + t14 * t6 * c3 * t64 - 0.2e1 * t14 * t30 * t29 * t26 + 0.2e1 * t9 * t35 * t29 * t26 + t9 * t13 * t3 - t14 * t7 * t3 + 0.2e1 * t57 * t41 * t26 - 0.2e1 * t75 * t41 * t26 - 0.2e1 * t57 * t50 * t26 + 0.2e1 * t75 * t50 * t26;
                                t94 = t19 * omega * t23;
                                t95 = t28 * t24;
                                t101 = t24 * t19 * t23;
                                t104 = t19 * t23;
                                t105 = t40 * t24;
                                pr = 0.2e1 * t42 * t105 * t104 + 0.2e1 * t46 * t105 * t104 + 0.2e1 * t42 * t95 * t94 + 0.2e1 * t46 * t95 * t94 - 0.2e1 * t51 * t101 - 0.2e1 * t54 * t101;
                // (u1,u2) = ur*rHat + vr*thetaHat 
                                u(i1,i2,i3,u1c) = ur*cosTheta - vr*sinTheta;
                                u(i1,i2,i3,u2c) = ur*sinTheta + vr*cosTheta;
                                u(i1,i2,i3,pc ) = pr;
                                if( maskLocal(i1,i2,i3) > 0 )
                                {
                                    u1Max=max(u1Max,u(i1,i2,i3,u1c));
                                    u2Max=max(u2Max,u(i1,i2,i3,u2c));
                                    pMax =max(pMax, u(i1,i2,i3,pc));
                                }
                            }
                            gridNorm = max(u1Max,u2Max,u3Max,pMax);
                        }
                        else
                        {
              // ** FOR NOW USE THE 2D SOLUTION ***
                            OV_ABORT("UDKS:ERROR: disk solution : 3D called");
                            FOR_3D(i1,i2,i3,I1,I2,I3)
                            {
                                Real x = vertex(i1,i2,i3,0);
                                Real y = vertex(i1,i2,i3,1);
                                Real z = vertex(i1,i2,i3,2);
                                Real r = sqrt( x*x + y*y );
                                if( r < eps ) //  avoid division by zero
                                {
                                    x=eps;           // avoid atan(0,0)
                                    y=eps;
                                    r=sqrt(2.)*eps;  // r=sqrt(x**2+y**2)
                                }   
                                theta = atan2(y,x); 
                                Real cosTheta = cos(theta);
                                Real sinTheta = sin(theta);
                                t1 = omega * r;
                                t2 = jn(n, t1);
                                t3 = t2 * c3;
                                t4 = 0.1e1 / r;
                                t5 = omega * t0;
                                t6 = cos(t5);
                                t7 = t6 * t4;
                                t8 = n * theta;
                                t9 = cos(t8);
                                t12 = sin(t5);
                                t13 = t12 * t4;
                                t14 = sin(t8);
                                t17 = n * n;
                                t19 = omega * omega;
                                t23 = 0.1e1 / (0.2e1 * t17 - t19 - 0.2e1 * n) * c3;
                                t24 = pow(r, n);
                                t26 = t4 * t24 * t23;
                                t27 = n + 0.1e1;
                                t28 = jn(t27, omega);
                                t29 = t28 * n;
                                t30 = t6 * omega;
                                t35 = t12 * omega;
                                t40 = jn(n, omega);
                                t41 = t40 * t17;
                                t42 = t9 * t6;
                                t46 = t14 * t12;
                                t50 = t40 * n;
                                t51 = t42 * t50;
                                t54 = t46 * t50;
                                ur = 0.2e1 * t14 * t35 * t29 * t26 + 0.2e1 * t9 * t30 * t29 * t26 + t14 * t13 * t3 - 0.2e1 * t42 * t41 * t26 - 0.2e1 * t46 * t41 * t26 + t9 * t7 * t3 + 0.2e1 * t51 * t26 + 0.2e1 * t54 * t26;
                                t57 = t14 * t6;
                                t63 = jn(t27, t1);
                                t64 = t63 / n * omega;
                                t75 = t9 * t12;
                                vr = -t9 * t12 * c3 * t64 + t14 * t6 * c3 * t64 - 0.2e1 * t14 * t30 * t29 * t26 + 0.2e1 * t9 * t35 * t29 * t26 + t9 * t13 * t3 - t14 * t7 * t3 + 0.2e1 * t57 * t41 * t26 - 0.2e1 * t75 * t41 * t26 - 0.2e1 * t57 * t50 * t26 + 0.2e1 * t75 * t50 * t26;
                                t94 = t19 * omega * t23;
                                t95 = t28 * t24;
                                t101 = t24 * t19 * t23;
                                t104 = t19 * t23;
                                t105 = t40 * t24;
                                pr = 0.2e1 * t42 * t105 * t104 + 0.2e1 * t46 * t105 * t104 + 0.2e1 * t42 * t95 * t94 + 0.2e1 * t46 * t95 * t94 - 0.2e1 * t51 * t101 - 0.2e1 * t54 * t101;
                // (u1,u2) = ur*rHat + vr*thetaHat 
                                u(i1,i2,i3,u1c) = ur*cosTheta - vr*sinTheta;
                                u(i1,i2,i3,u2c) = ur*sinTheta + vr*cosTheta;
                                u(i1,i2,i3,u3c) = 0.;
                                u(i1,i2,i3,pc ) = pr; 
                            }   
                        }
                    printF(" +++++ UDKS: grid=%d : gridNorm=%9.2e numberOfTimeDerivatives=%d +++++",
                              grid,gridNorm,numberOfTimeDerivatives);
                    uNorm = max(uNorm,gridNorm);
                }
                if( uNorm==0. )
                {
                    printF("UDKS: ERROR: disk solution is zero! uNorm=%e\n",uNorm);
                    OV_ABORT("ERROR");
                }        
                diskScaling = c3/uNorm;
                printF("\n +++++ UDKS: uNorm=%9.2e +++++",uNorm);
            }
            printF("\n +++++ UDKS: t=%9.3e, diskScaling=%9.2e +++++",t, diskScaling);
            c3 = diskScaling;       

                Real ur,vr,pr; 
                int i1,i2,i3;
                Real theta; 
                const Real eps=1.e-9;   // ********************************** FIX ME *********
                if( t<= dt )
                    printF("\n *** WARNING *** need a better fix for disk solution at r=0 ******\n\n");
                Real u1Max=0., u2Max=0., u3Max=0., pMax=0.; 
                if( numberOfDimensions==2 )
                {
                    FOR_3D(i1,i2,i3,I1,I2,I3)
                    {
            // Reference coordinates:
                        Real x = vertex(i1,i2,i3,0);
                        Real y = vertex(i1,i2,i3,1);
                        Real r = sqrt( x*x + y*y );
                        if( r < eps ) //  avoid division by zero
                        {
                            x=eps;           // avoid atan(0,0)
                            y=eps;
                            r=sqrt(2.)*eps;  // r=sqrt(x**2+y**2)
                        }   
                        theta = atan2(y,x); 
                        Real cosTheta = cos(theta);
                        Real sinTheta = sin(theta);
                        t1 = omega * r;
                        t2 = jn(n, t1);
                        t3 = t2 * c3;
                        t4 = 0.1e1 / r;
                        t5 = omega * t;
                        t6 = cos(t5);
                        t7 = t6 * t4;
                        t8 = n * theta;
                        t9 = cos(t8);
                        t12 = sin(t5);
                        t13 = t12 * t4;
                        t14 = sin(t8);
                        t17 = n * n;
                        t19 = omega * omega;
                        t23 = 0.1e1 / (0.2e1 * t17 - t19 - 0.2e1 * n) * c3;
                        t24 = pow(r, n);
                        t26 = t4 * t24 * t23;
                        t27 = n + 0.1e1;
                        t28 = jn(t27, omega);
                        t29 = t28 * n;
                        t30 = t6 * omega;
                        t35 = t12 * omega;
                        t40 = jn(n, omega);
                        t41 = t40 * t17;
                        t42 = t9 * t6;
                        t46 = t14 * t12;
                        t50 = t40 * n;
                        t51 = t42 * t50;
                        t54 = t46 * t50;
                        ur = 0.2e1 * t14 * t35 * t29 * t26 + 0.2e1 * t9 * t30 * t29 * t26 + t14 * t13 * t3 - 0.2e1 * t42 * t41 * t26 - 0.2e1 * t46 * t41 * t26 + t9 * t7 * t3 + 0.2e1 * t51 * t26 + 0.2e1 * t54 * t26;
                        t57 = t14 * t6;
                        t63 = jn(t27, t1);
                        t64 = t63 / n * omega;
                        t75 = t9 * t12;
                        vr = -t9 * t12 * c3 * t64 + t14 * t6 * c3 * t64 - 0.2e1 * t14 * t30 * t29 * t26 + 0.2e1 * t9 * t35 * t29 * t26 + t9 * t13 * t3 - t14 * t7 * t3 + 0.2e1 * t57 * t41 * t26 - 0.2e1 * t75 * t41 * t26 - 0.2e1 * t57 * t50 * t26 + 0.2e1 * t75 * t50 * t26;
                        t94 = t19 * omega * t23;
                        t95 = t28 * t24;
                        t101 = t24 * t19 * t23;
                        t104 = t19 * t23;
                        t105 = t40 * t24;
                        pr = 0.2e1 * t42 * t105 * t104 + 0.2e1 * t46 * t105 * t104 + 0.2e1 * t42 * t95 * t94 + 0.2e1 * t46 * t95 * t94 - 0.2e1 * t51 * t101 - 0.2e1 * t54 * t101;
            // (u1,u2) = ur*rHat + vr*thetaHat 
                        u(i1,i2,i3,u1c) = ur*cosTheta - vr*sinTheta;
                        u(i1,i2,i3,u2c) = ur*sinTheta + vr*cosTheta;
                        u(i1,i2,i3,pc ) = pr;
                        if( maskLocal(i1,i2,i3) > 0 )
                        {
                            u1Max=max(u1Max,u(i1,i2,i3,u1c));
                            u2Max=max(u2Max,u(i1,i2,i3,u2c));
                            pMax =max(pMax, u(i1,i2,i3,pc));
                        }
                    }
                    uNorm = max(u1Max,u2Max,u3Max,pMax);
                }
                else
                {
          // ** FOR NOW USE THE 2D SOLUTION ***
                    OV_ABORT("UDKS:ERROR: disk solution : 3D called");
                    FOR_3D(i1,i2,i3,I1,I2,I3)
                    {
                        Real x = vertex(i1,i2,i3,0);
                        Real y = vertex(i1,i2,i3,1);
                        Real z = vertex(i1,i2,i3,2);
                        Real r = sqrt( x*x + y*y );
                        if( r < eps ) //  avoid division by zero
                        {
                            x=eps;           // avoid atan(0,0)
                            y=eps;
                            r=sqrt(2.)*eps;  // r=sqrt(x**2+y**2)
                        }   
                        theta = atan2(y,x); 
                        Real cosTheta = cos(theta);
                        Real sinTheta = sin(theta);
                        t1 = omega * r;
                        t2 = jn(n, t1);
                        t3 = t2 * c3;
                        t4 = 0.1e1 / r;
                        t5 = omega * t;
                        t6 = cos(t5);
                        t7 = t6 * t4;
                        t8 = n * theta;
                        t9 = cos(t8);
                        t12 = sin(t5);
                        t13 = t12 * t4;
                        t14 = sin(t8);
                        t17 = n * n;
                        t19 = omega * omega;
                        t23 = 0.1e1 / (0.2e1 * t17 - t19 - 0.2e1 * n) * c3;
                        t24 = pow(r, n);
                        t26 = t4 * t24 * t23;
                        t27 = n + 0.1e1;
                        t28 = jn(t27, omega);
                        t29 = t28 * n;
                        t30 = t6 * omega;
                        t35 = t12 * omega;
                        t40 = jn(n, omega);
                        t41 = t40 * t17;
                        t42 = t9 * t6;
                        t46 = t14 * t12;
                        t50 = t40 * n;
                        t51 = t42 * t50;
                        t54 = t46 * t50;
                        ur = 0.2e1 * t14 * t35 * t29 * t26 + 0.2e1 * t9 * t30 * t29 * t26 + t14 * t13 * t3 - 0.2e1 * t42 * t41 * t26 - 0.2e1 * t46 * t41 * t26 + t9 * t7 * t3 + 0.2e1 * t51 * t26 + 0.2e1 * t54 * t26;
                        t57 = t14 * t6;
                        t63 = jn(t27, t1);
                        t64 = t63 / n * omega;
                        t75 = t9 * t12;
                        vr = -t9 * t12 * c3 * t64 + t14 * t6 * c3 * t64 - 0.2e1 * t14 * t30 * t29 * t26 + 0.2e1 * t9 * t35 * t29 * t26 + t9 * t13 * t3 - t14 * t7 * t3 + 0.2e1 * t57 * t41 * t26 - 0.2e1 * t75 * t41 * t26 - 0.2e1 * t57 * t50 * t26 + 0.2e1 * t75 * t50 * t26;
                        t94 = t19 * omega * t23;
                        t95 = t28 * t24;
                        t101 = t24 * t19 * t23;
                        t104 = t19 * t23;
                        t105 = t40 * t24;
                        pr = 0.2e1 * t42 * t105 * t104 + 0.2e1 * t46 * t105 * t104 + 0.2e1 * t42 * t95 * t94 + 0.2e1 * t46 * t95 * t94 - 0.2e1 * t51 * t101 - 0.2e1 * t54 * t101;
            // (u1,u2) = ur*rHat + vr*thetaHat 
                        u(i1,i2,i3,u1c) = ur*cosTheta - vr*sinTheta;
                        u(i1,i2,i3,u2c) = ur*sinTheta + vr*cosTheta;
                        u(i1,i2,i3,u3c) = 0.;
                        u(i1,i2,i3,pc ) = pr; 
                    }   
                }

                if( computeVelocity )
                {
                    printF("\n $$$$$$$$$$ userDefinedKnownSolution: EVALUATE VELOCITY, ACCEL and PRESSURE AT PAST TIMES t=%9.3e dt=%16.6e $$$$$$$\n\n",t,dt);
                    assert( dt>0. );
                    assert( vgf!=NULL );
                    const int currentVelocity = dbase.get<int>("currentVelocity");
                    assert( currentVelocity>=0 && currentVelocity<numberOfVelocityFunctions );
          // fill in t0=0 and some past time levels
                    for( int level=0; level<numberOfVelocityFunctions; level++ )
                    {
                        real tc = 0. - level*dt;
                        const int prevVelocity = (currentVelocity - level + numberOfVelocityFunctions) % numberOfVelocityFunctions;
                        OV_GET_SERIAL_ARRAY(real,vgf[prevVelocity][grid],vLocal);
                        OV_GET_SERIAL_ARRAY(real,vtgf[prevVelocity][grid],vtLocal);
                        FOR_3D(i1,i2,i3,I1,I2,I3)
                        {
                                Real x = vertex(i1,i2,i3,0);
                                Real y = vertex(i1,i2,i3,1);
                                Real r = sqrt( x*x + y*y );
                                if( r < eps ) //  avoid division by zero
                                {
                                    x=eps;           // avoid atan(0,0)
                                    y=eps;
                                    r=sqrt(2.)*eps;  // r=sqrt(x**2+y**2)
                                }   
                                Real theta = atan2(y,x); 
                                Real cosTheta = cos(theta);
                                Real sinTheta = sin(theta);
              // Real x = vertex(i1,i2,i3,0);
              // Real y = vertex(i1,i2,i3,1);
              // Real r = sqrt( x*x + y*y );
              // Real theta = atan2(y,x); 
              // Real cosTheta = x/r;
              // Real sinTheta = y/r;
                            t1 = n * n;
                            t2 = omega * omega;
                            t4 = t1 - t2 / 0.2e1 - n;
                            t5 = omega * r;
                            t6 = jn(n, t5);
                            t8 = pow(r, n);
                            t9 = n + 0.1e1;
                            t10 = jn(t9, omega);
                            t12 = jn(n, omega);
                            t18 = t6 * t4 - n * (-t10 * omega + (n - 0.1e1) * t12) * t8;
                            t19 = omega * tc;
                            t20 = cos(t19);
                            t21 = n * theta;
                            t22 = sin(t21);
                            t24 = sin(t19);
                            t25 = cos(t21);
                            t34 = 0.1e1 / (0.2e1 * t1 - t2 - 0.2e1 * n) * omega;
                            t35 = 0.1e1 / r;
                            ur = 0.2e1 * t35 * t34 * c3 * (t22 * t20 - t25 * t24) * t18;
                            t43 = jn(t9, t5);
                            vr = 0.2e1 / n * t35 * t34 * (-t43 * omega * t4 * r + n * t18) * c3 * (t25 * t20 + t22 * t24);
              // (u1,u2) = ur*rHat + vr*thetaHat 
                            vLocal(i1,i2,i3,u1c) = ur*cosTheta - vr*sinTheta;
                            vLocal(i1,i2,i3,u2c) = ur*sinTheta + vr*cosTheta;
                            t1 = n * n;
                            t2 = omega * omega;
                            t4 = t1 - t2 / 0.2e1 - n;
                            t5 = omega * r;
                            t6 = jn(n, t5);
                            t8 = pow(r, n);
                            t9 = n + 0.1e1;
                            t10 = jn(t9, omega);
                            t12 = jn(n, omega);
                            t18 = t6 * t4 - n * (-t10 * omega + (n - 0.1e1) * t12) * t8;
                            t19 = omega * tc;
                            t20 = cos(t19);
                            t21 = n * theta;
                            t22 = cos(t21);
                            t24 = sin(t19);
                            t25 = sin(t21);
                            t34 = 0.1e1 / (0.2e1 * t1 - t2 - 0.2e1 * n) * t2;
                            t35 = 0.1e1 / r;
                            ur = -0.2e1 * t35 * t34 * c3 * (t22 * t20 + t25 * t24) * t18;
                            t44 = jn(t9, t5);
                            vr = 0.2e1 / n * t35 * t34 * (-t44 * omega * t4 * r + n * t18) * c3 * (t25 * t20 - t22 * t24);
              // (u1,u2) = ur*rHat + vr*thetaHat 
                            vtLocal(i1,i2,i3,u1c) = ur*cosTheta - vr*sinTheta;
                            vtLocal(i1,i2,i3,u2c) = ur*sinTheta + vr*cosTheta;        
                        }      
                    }
          // --- save pressure for time extrapolation ---
                    const int currentPressure = dbase.get<int>("currentPressure");
                    assert( currentPressure>=0 && currentPressure<max(1,numberOfPressureFunctions) );
                    for( int level=0; level<numberOfPressureFunctions; level++ )
                    {
                        real tc = 0. - (level+3)*dt;  // -- save pressure starting at t-3*dt 
                        const int prevPressure = (currentPressure - level + numberOfPressureFunctions) % numberOfPressureFunctions;
                        assert( pgf !=NULL );
                        OV_GET_SERIAL_ARRAY(real,pgf[prevPressure][grid],pLocal);
                        FOR_3D(i1,i2,i3,I1,I2,I3)
                        {
                                Real x = vertex(i1,i2,i3,0);
                                Real y = vertex(i1,i2,i3,1);
                                Real r = sqrt( x*x + y*y );
                                if( r < eps ) //  avoid division by zero
                                {
                                    x=eps;           // avoid atan(0,0)
                                    y=eps;
                                    r=sqrt(2.)*eps;  // r=sqrt(x**2+y**2)
                                }   
                                Real theta = atan2(y,x); 
                                Real cosTheta = cos(theta);
                                Real sinTheta = sin(theta);
              // Real x = vertex(i1,i2,i3,0);
              // Real y = vertex(i1,i2,i3,1);
              // Real r = sqrt( x*x + y*y );
              // Real theta = atan2(y,x); 
              // Real cosTheta = x/r;
              // Real sinTheta = y/r;
                            t1 = omega * r;
                            t2 = jn(n, t1);
                            t3 = t2 * c3;
                            t4 = 0.1e1 / r;
                            t5 = omega * tc;
                            t6 = cos(t5);
                            t7 = t6 * t4;
                            t8 = n * theta;
                            t9 = cos(t8);
                            t12 = sin(t5);
                            t13 = t12 * t4;
                            t14 = sin(t8);
                            t17 = n * n;
                            t19 = omega * omega;
                            t23 = 0.1e1 / (0.2e1 * t17 - t19 - 0.2e1 * n) * c3;
                            t24 = pow(r, n);
                            t26 = t4 * t24 * t23;
                            t27 = n + 0.1e1;
                            t28 = jn(t27, omega);
                            t29 = t28 * n;
                            t30 = t6 * omega;
                            t35 = t12 * omega;
                            t40 = jn(n, omega);
                            t41 = t40 * t17;
                            t42 = t9 * t6;
                            t46 = t14 * t12;
                            t50 = t40 * n;
                            t51 = t42 * t50;
                            t54 = t46 * t50;
                            ur = 0.2e1 * t14 * t35 * t29 * t26 + 0.2e1 * t9 * t30 * t29 * t26 + t14 * t13 * t3 - 0.2e1 * t42 * t41 * t26 - 0.2e1 * t46 * t41 * t26 + t9 * t7 * t3 + 0.2e1 * t51 * t26 + 0.2e1 * t54 * t26;
                            t57 = t14 * t6;
                            t63 = jn(t27, t1);
                            t64 = t63 / n * omega;
                            t75 = t9 * t12;
                            vr = -t9 * t12 * c3 * t64 + t14 * t6 * c3 * t64 - 0.2e1 * t14 * t30 * t29 * t26 + 0.2e1 * t9 * t35 * t29 * t26 + t9 * t13 * t3 - t14 * t7 * t3 + 0.2e1 * t57 * t41 * t26 - 0.2e1 * t75 * t41 * t26 - 0.2e1 * t57 * t50 * t26 + 0.2e1 * t75 * t50 * t26;
                            t94 = t19 * omega * t23;
                            t95 = t28 * t24;
                            t101 = t24 * t19 * t23;
                            t104 = t19 * t23;
                            t105 = t40 * t24;
                            pr = 0.2e1 * t42 * t105 * t104 + 0.2e1 * t46 * t105 * t104 + 0.2e1 * t42 * t95 * t94 + 0.2e1 * t46 * t95 * t94 - 0.2e1 * t51 * t101 - 0.2e1 * t54 * t101;
                            pLocal(i1,i2,i3) = pr;
                        }
                    }
                }

        }    
        else
        {
            OV_ABORT("error");
        }

    }

    else if( userKnownSolution=="gaussianPlaneWave" )
    {
    // --- Gaussian plane wave : incompressible elasticity -----
        const int pc = dbase.get<int >("pc");
        assert( pc>=0 );

        MappedGrid & mg = cg[grid];
        mg.update(MappedGrid::THEvertex | MappedGrid::THEcenter);
        OV_GET_SERIAL_ARRAY_CONST(real,mg.vertex(),vertex);
        RealArray & u = ua;

                      
        const Real kx   = rpar[ 0];
        const Real ky   = rpar[ 1];      
        const Real kz   = rpar[ 2];      
        const Real k0   = rpar[ 3];      
        const Real beta = rpar[ 4];      
        const Real ax   = rpar[ 5];      
        const Real ay   = rpar[ 6];      
        const Real az   = rpar[ 7]; 
        const Real x0   = rpar[ 8];      
        const Real y0   = rpar[ 9];      
        const Real z0   = rpar[10];   

        const Real kNorm = sqrt( kx*kx + ky*ky + kz*kz );
        const Real omega = cs*kNorm;  
                
        if( t<3*dt )
            printF("UDKS: evaluate the Gaussan Plane Wave kx=%g, ky=%g, kz=%g, k0=%g, beta=%g, (ax,ay,az)=(%g,%g,%g) cs=%g, t=%9.3e\n",kx,ky,kz,k0,beta,ax,ay,az,cs,t);

        if( numberOfDimensions==2 )
        {
            FOR_3D(i1,i2,i3,I1,I2,I3)
            {
                const Real x = vertex(i1,i2,i3,0);
                const Real y = vertex(i1,i2,i3,1);

                const Real xi = kx*(x-x0) + ky*(y-y0) - omega*t; 
                const Real expxi = exp( -beta*xi*xi )*cos(k0*xi); 

                u(i1,i2,i3,u1c) = ax*expxi;
                u(i1,i2,i3,u2c) = ay*expxi;
                u(i1,i2,i3,pc ) = 0.;
            }
        }
        else
        {
            FOR_3D(i1,i2,i3,I1,I2,I3)
            {
                const Real x = vertex(i1,i2,i3,0);
                const Real y = vertex(i1,i2,i3,1);
                const Real z = vertex(i1,i2,i3,2);

                const Real xi = kx*(x-x0) + ky*(y-y0) + kz*(z-z0) - omega*t; 
                const Real expxi = exp( -beta*xi*xi )*cos(k0*xi); 

                u(i1,i2,i3,u1c) = ax*expxi;
                u(i1,i2,i3,u2c) = ay*expxi;
                u(i1,i2,i3,u3c) = az*expxi;
                u(i1,i2,i3,pc ) = 0.;
            }   

        }       

    }
    else if( userKnownSolution=="linearBeamExactSolution" )
    {
  
        MappedGrid & mg = cg[grid];
        mg.update(MappedGrid::THEvertex | MappedGrid::THEcenter);

        OV_GET_SERIAL_ARRAY_CONST(real,mg.vertex(),vertex);

        RealArray & u = ua;

        const int & uc = dbase.get<int >("uc");   //  u velocity component =u(all,all,all,uc)
        const int & vc = dbase.get<int >("vc");  
        const int & pc = dbase.get<int >("pc");  
        double E=1.4e6;
        double rhos=10000.0;
        double h=0.02;
        double Ioverb=6.6667e-7;
        double rhof=1000;
        double nu=0.001;
        double L=0.3;
        double H=0.3;
        double k=2.0*3.141592653589/L;
        double omega0=sqrt(E*Ioverb*k*k*k*k/(rhos*h));
        double what = 0.00001;
    //double beta=1.0/nu*sqrt(E*Ioverb/(rhos*h));
    //std::complex<double> omegatilde(1.065048891,-5.642079778e-4);
        std::cout << "t = " << t << std::endl;
        double omegar = 0.8907148069, omegai = -0.9135887123e-2;
        for ( int i3=I3.getBase(); i3<=I3.getBound(); i3++ ) {
            for ( int i2=I2.getBase(); i2<=I2.getBound(); i2++ ) {
                for ( int i1=I1.getBase(); i1<=I1.getBound(); i1++ ) {
                    
                    double y = vertex(i1,i2,i3,1);
                    double x = vertex(i1,i2,i3,0);
                    
                    BeamModel::exactSolutionVelocity(x,y,t,k,H,
                                                                                      omegar,omegai, 
                                                                                      omega0,nu,
                                                                                      what,u(i1,i2,i3,uc),
                                                                                      u(i1,i2,i3,vc));

                    BeamModel::exactSolutionPressure(x,y,t,k,H,
                                                                                      omegar,omegai, 
                                                                                      omega0,nu,
                                                                                      what,u(i1,i2,i3,pc));

          //std::cout << x << " " << y << " " << u(i1,i2,i3,uc) << " " << u(i1,i2,i3,vc) << std::endl;
                }
            }
        }   
    }

    else if( userKnownSolution=="bulkSolidPiston" )
    {
    // ---- return the exact solution for the FSI INS+elastic piston ---
    //     y_I(t) = F(t + Hbar/cp) - F(t - Hbar/cp)
    //        F(z) = amp * R(z) * sin( 2*Pi*k(t-t0) )
    //        R(z) = ramp function that smoothly transitions from 0 to 1 

    // assert( v1c>=0 && u1c>=0 && s22c >=0 );

        const real & H        = rpar[0];
        const real & Hbar     = rpar[1];
        const real & rhoFluid = rpar[2];
        const real & rhoBar   = rpar[3];
        const real & lambdaBar= rpar[4];
        const real & muBar    = rpar[5];
        const real & thetaR   = rpar[6]; // rotation of domain (radians)

        const real cp2 = sqrt((lambdaBar+2.*muBar)/rhoBar);
        assert( cp==cp2 ); // sanity check
        
        if( t<= 2.*dt )
        {
            printF("--SM-- userDefinedKnownSolution: bulkSolidPiston, Hbar=%g, t=%9.3e\n",Hbar,t);
            printF("--SM-- lambda   =%g mu   =%g rho   =%g\n",lambda,mu,rho);
            printF("--SM-- lambdaBar=%g muBar=%g rhoBar=%g\n",lambdaBar,muBar,rhoBar);
            
        }
        assert( lambda==lambdaBar && mu==muBar && rho==rhoBar );
    // assert( numberOfTimeDerivatives==0 );

        TimeFunction & bsp = db.get<TimeFunction>("timeFunctionBSP");

        MappedGrid & mg = cg[grid];
        mg.update(MappedGrid::THEvertex | MappedGrid::THEcenter);

        const realArray & center = mg.center();
        RealArray & u = ua;

        const real ct = cos(thetaR);
        const real st = sin(thetaR);

        int i1,i2,i3;
        FOR_3D(i1,i2,i3,I1,I2,I3)
        {
      // Reference coordinates:
            const real x= center(i1,i2,i3,0);
            const real y= center(i1,i2,i3,1);
            const real yRef = -st*x+ct*y;

      // printF("thetaR=%f\n",thetaR);
            real xim,xip, fm,fp, fmd,fpd;
            xim=t-(yRef+Hbar)/cp;
            xip=t+(yRef+Hbar)/cp;  
                
            bsp.eval(xim, fm,fmd );  // fmd = d(fm(xi))/d(xi)
            bsp.eval(xip, fp,fpd );

            u(i1,i2,i3,u1c)= -(fp - fm)*st;
            u(i1,i2,i3,u2c)=  (fp - fm)*ct;

      // velocities
            if( assignVelocity )
            {
                if( numberOfTimeDerivatives==0 )
                {
                    u(i1,i2,i3,v1c)= -(fpd - fmd)*st;
                    u(i1,i2,i3,v2c)=  (fpd - fmd)*ct;
                }
                else if( numberOfTimeDerivatives==1 )
                {
          // ---- return the acceleration ----
          // eval F''
                    real fmdd,fpdd;
                    bsp.evalDerivative(xim, fmdd, 2 ); // 2 derivatives 
                    bsp.evalDerivative(xip, fpdd, 2 );        
                    u(i1,i2,i3,v1c)= -(fpdd - fmdd)*st;
                    u(i1,i2,i3,v2c)=  (fpdd - fmdd)*ct;
                }
                else
                {
                    OV_ABORT("error");
                }
                
            }
            
      // stresses
            real u2y;
            if( assignStress )
            {
                u2y = (fpd + fmd)/cp; // note "+" sign
                
                const real s11 = lambda*u2y;
                const real s22 = (lambda+2.*mu)*u2y;

        // e1 = [ct,st], e2 = [-st,ct]
        // sigma = s11 e1 e1^T + s22 e2 e2^T

                u(i1,i2,i3,s11c)= s11*ct*ct+s22*st*st;
                u(i1,i2,i3,s12c)= s11*st*ct-s22*st*ct;
                u(i1,i2,i3,s21c)= s11*st*ct-s22*st*ct;
                u(i1,i2,i3,s22c)= s11*st*st+s22*ct*ct;;
            }
            
        }

    }

    else if( userKnownSolution=="radialElasticPiston" )
    {
    // ---- return the exact solution for radial elastic piston ----

        assert( numberOfTimeDerivatives==0 );

    // -- we could avoid building the vertex array on Cartesian grids ---
        GET_VERTEX_ARRAY(xLocal);
        const real & R        = rpar[0];
        const real & Rbar     = rpar[1];
        const real & rho      = rpar[2];
        const real & rhoBar   = rpar[3];
        const real & lambdaBar= rpar[4];
        const real & muBar    = rpar[5];
        const real & k        = rpar[6];

        const real cp = sqrt((lambdaBar+2.*muBar)/rhoBar);
        
    // uI = uI(t) =  interface displacement in the radial direction 
    // eval uI and vI = uI_t 
        real uI,vI;
        TimeFunction & bsp = db.get<TimeFunction>("timeFunctionREP");
        bsp.eval(t, uI,vI );  

        if( t <= 2.*dt )
        {
            printF("--SM-- getUserDefinedDeformingBodyKnownSolution: radialElasticPiston, t=%9.3e uI=%9.3e vI=%9.3e, Rbar=%6.3f\n",
                          t,uI,vI,Rbar );
        }

        RealArray & u = ua;
        const real eps = 10.*REAL_EPSILON;
        const real sqrtEps = sqrt(REAL_EPSILON);
        
        int i1,i2,i3;
        FOR_3D(i1,i2,i3,I1,I2,I3)
        {
      // Reference coordinates:
            real x= xLocal(i1,i2,i3,0);
            real y= xLocal(i1,i2,i3,1);
            real r = sqrt( SQR(x) + SQR(y) );
            real ct,st;
            if( r>eps )
            {
                ct=x/r; st=y/r;
            }
            else
            {
                ct=1.; st=0.;  // at the origin we just pick an angle, should not matter
            }
            
            
            real kr=k*r;
            real jnkr = jn(1,kr);
            real ur = uI*jnkr;    // radial displacement 
            real vr = vI*jnkr;    // radial velocity

            u(i1,i2,i3,u1c)=ur*ct;
            u(i1,i2,i3,u2c)=ur*st;

      // velocities
            if( assignVelocity )
            {
                u(i1,i2,i3,v1c)=vr*ct;
                u(i1,i2,i3,v2c)=vr*st;
            }
            
      // stresses
            if( assignStress )
            {
                real jnkrp = .5*k*(jn(0,kr)-jn(2,kr));  // Jn' = .5*( J(n-1) - J(n+1) ) check me 
        // real jnkrp = (jn(0,kr)-jnkr)/r;  // J1'(z) = (J0(z)-J1(z))/z
                
        // ur = amp*J_1(k*r)* ...
                real urr=uI*jnkrp;      // r-derivative of the radial displacement
                
                real urByr;
                if( fabs(r)>sqrtEps )
                    urByr=ur/r;
                else
                {
                    urByr=urr;  // for r small, ur/r = (ur(r)-ur(0))/r =  urr(r) + ...
          // printF(" --UDKDBS: i=(%i,%i) r=%8.2e urr=%9.3e ur=%9.3e ct=%9.3e st=%9.3e\n",i1,i2,r,urr,ur,ct,st);
                }
                
                real sigmarr = (lambdaBar+2.*muBar)*urr + lambdaBar*urByr;
                real sigmart=0.;
                real sigmatt = lambdaBar*urr +  (lambdaBar+2.*muBar)*urByr;

        // Cartesian components of the stress tensor:
        //  [ s11 s12 ] = srr rHat rHat^T + srt rHat thetaHat^T + str thetaHat^t rHat + stt thetaHat thetaHat^T
        //  [ s21 s22 ]
        // where
        //   rHat=[cos,sin], thetaHat=[-sin,cos]
        /// **CHECK ME**
                u(i1,i2,i3,s11c)= sigmarr*ct*ct - 2.*sigmart*ct*st + sigmatt*st*st;
                u(i1,i2,i3,s12c)= sigmarr*ct*st + sigmart*(ct*ct-st*st) - sigmatt*ct*st ;
                u(i1,i2,i3,s21c)= u(i1,i2,i3,s12c);
                u(i1,i2,i3,s22c)= sigmarr*st*st + 2.*sigmart*ct*st + sigmatt*ct*ct;
            }
            
        }
        
    }
    else if ( userKnownSolution == "fibShear" ) {
    // --------------------------------------------------------------------------------
    // ------ Exact solution for a parallel flow shearing a bulk elastic solid --------
    //  \bar{u}_1(y,t) = amp         exp(i omega t) (A cos(ks y) + B sin( ks y))
    //      {v}_1(y,t) = amp i omega exp(i omega t) (C exp(kf y) + D exp(-kf y))
    //              ks = omega / cs
    //              kf = sqrt(i rho omega / mu)
    //  Parameters:
    //  amp    : maximum amplitude of the displacement 
    //  omega  : time frequency of solution 
    //  H,Hbar : Height of fluid and solid domains
    //  rhoBar,lambaBar,muBar : solid density and Lame parameters
    // --------------------------------------------------------------------------------

        const real & omegar = rpar[0];
        const real & omegai = rpar[1];
        const real & ar     = rpar[2];
        const real & ai     = rpar[3];
        const real & br     = rpar[4];
        const real & bi     = rpar[5];
        const real & cr     = rpar[6];
        const real & ci     = rpar[7];
        const real & dr     = rpar[8];
        const real & di     = rpar[9];
        const real & ksr    = rpar[10];
        const real & ksi    = rpar[11];
        const real & kfr    = rpar[12];
        const real & kfi    = rpar[13];
        const real & amp    = rpar[14];
        const real & mu     = rpar[15];
        const real & thetaR = rpar[16];
        const real & muBar  = rpar[17];

        printF("--SM-- userDefinedKnownSolution: fibShear, t=%9.3e, "
                      "rhoBar=%9.3e, muBar=%9.3e\n",t,rho,muBar);

        if( numberOfTimeDerivatives!=0 )
        {
            printF("\n --SM--UDKS:fibShear: ERROR: numberOfTimeDerivatives=%i NOT IMPLEMENTED! *****FIX ME**** \n\n",
                  numberOfTimeDerivatives);
        }

    // const real cs2 = sqrt((muBar)/rhoBar);

    // sanity checks
    // assert( cs==cs2 ); 
    // assert( mu==muBar );


    // fill in the solution
        MappedGrid & mg = cg[grid];
        mg.update(MappedGrid::THEvertex | MappedGrid::THEcenter);

        const realArray & center = mg.center();
        RealArray & u = ua;

        real u0_r = amp*exp(-omegai*t)*cos(omegar*t);
        real u0_i = amp*exp(-omegai*t)*sin(omegar*t);
        real v0_r =  omegai*u0_r + omegar*u0_i;
        real v0_i = -omegar*u0_r + omegai*u0_i;

        const real ct = cos(thetaR);
        const real st = sin(thetaR);

        int i1,i2,i3;
        FOR_3D(i1,i2,i3,I1,I2,I3)
        {
      // Reference coordinates:
            const real x= center(i1,i2,i3,0);
            const real y= center(i1,i2,i3,1);
            real yRef = -st*x+ct*y;      
      // Evaluate the solution for shear flow (FSI)
      //
      // u1  = amp    ( A cos(ks y) + B sin(ks y)) exp(-i omega t)
      // u1y = amp ks (-A sin(ks y) + B cos(ks y)) exp(-i omega t)
      //
      // Return:
      //  ur  = real( A cos(ks y) + B sin(ks y) )
      //  ui  = imag( A cos(ks y) + B sin(ks y) )
      //  uyr = real( ks (-A sin(ks y) + B cos(ks y)) )
      //  uyi = imag( ks (-A sin(ks y) + B cos(ks y)) )

            real ur,ui,vr,vi,uyr,uyi;
            evalFibShearSolidFull(ksr,ksi,ar,ai,br,bi,yRef,t,ur,ui,vr,vi,uyr,uyi,omegar,omegai);

      // u(i1,i2,i3,u1c)=u0_r*ur-u0_i*ui;
            u(i1,i2,i3,u1c)=amp*ur*ct;
            u(i1,i2,i3,u2c)=amp*ur*st;


      // velocities
            if( assignVelocity )
            {
        // u(i1,i2,i3,v1c)=v0_r*ur-v0_i*ui;
                if( numberOfTimeDerivatives==0 )
                {
                    u(i1,i2,i3,v1c)=amp*vr*ct;
                    u(i1,i2,i3,v2c)=amp*vr*st;
                }
                else
                {
                    printF("--SM--UDKS: WARNING numberOfTimeDerivatives=%i ! Setting acceleration to zero\n");
                    u(i1,i2,i3,v1c)=0.;
                    u(i1,i2,i3,v2c)=0.;
                    
                }
                
            }
            
      // stresses
            if( assignStress )
            {
                if( numberOfTimeDerivatives==0 )
                {
                    const real srt = amp*muBar*uyr;
                    u(i1,i2,i3,s11c)=-2.*ct*st*srt;
                    u(i1,i2,i3,s12c)=(SQR(ct)-SQR(st))*srt;
                    u(i1,i2,i3,s21c)=(SQR(ct)-SQR(st))*srt;
                    u(i1,i2,i3,s22c)= 2.*ct*st*srt;
                }
                else
                {
                    printF("--SM--UDKS: WARNING numberOfTimeDerivatives=%i ! Setting stress-rate to zero\n");
                    u(i1,i2,i3,s11c)=0.;
                    u(i1,i2,i3,s12c)=0.;
                    u(i1,i2,i3,s21c)=0.;
                    u(i1,i2,i3,s22c)=0.;
                }
            }
            
        }


    }
    else if ( userKnownSolution == "radialFibShear" ) 
    {
    // --------------------------------------------------------------------------------
    // ------ Exact solution for a radial flow shearing a bulk elastic solid ----------
    // \bar{u}_theta(y,t) = amp exp(i omega t) (A J1(ks y) + B Y1(ks y))
    //     {v}_theta(y,t) = amp exp(i omega t) (C J1(kf y) + D Y1(kf y))
    //              ks = omega / cs
    //              kf = sqrt(- i rho omega / mu)
    //  Parameters:
    //  amp    : maximum amplitude of the displacement 
    //  omega  : time frequency of solution 
    //  H,Hbar : Height of fluid and solid domains
    //  rhoBar,lambaBar,muBar : solid density and Lame parameters
    // --------------------------------------------------------------------------------

        const real & omegar = rpar[0];
        const real & omegai = rpar[1];
        const real & ar     = rpar[2];
        const real & ai     = rpar[3];
        const real & br     = rpar[4];
        const real & bi     = rpar[5];
        const real & cr     = rpar[6];
        const real & ci     = rpar[7];
        const real & R      = rpar[8];
        const real & Rbar   = rpar[9];
        const real & ksr    = rpar[10];
        const real & ksi    = rpar[11];
        const real & kfr    = rpar[12];
        const real & kfi    = rpar[13];
        const real & amp    = rpar[14];
        const real & muFluid= rpar[15];
        const real & muBar  = rpar[16];

        printF("--SM-- userDefinedKnownSolution: radialFibShear, t=%9.3e, "
                      "rhoBar=%9.3e, muBar=%9.3e\n",t,rho,mu);

        assert( mu==muBar );
        
        assert( numberOfTimeDerivatives==0 );


    // fill in the solution
        MappedGrid & mg = cg[grid];
        mg.update(MappedGrid::THEvertex | MappedGrid::THEcenter);
        const real eps = 10.*REAL_EPSILON;

        const realArray & center = mg.center();
        RealArray & u = ua;

        int i1,i2,i3;
        FOR_3D(i1,i2,i3,I1,I2,I3)
        {
      // Reference coordinates:
            const real x= center(i1,i2,i3,0);
            const real y= center(i1,i2,i3,1);
            const real r= sqrt(SQR(x) + SQR(y));

      // compute trig functions
            real cosTheta, sinTheta;
            if( fabs(r)>eps )
                {
                    cosTheta=x/r; sinTheta=y/r;
                }
            else
                {
                    cosTheta=1.; sinTheta=0.;
                }
            
      // Evaluate the solution for shear flow (FSI)
      //
      // u_theta = 
      // v_theta = 
      // a_theta = 
      // d/dr u_theta = 
      // d/dr v_theta = 
      //
      // srt is the stress without mubar
            real ut,vt,at,srt,vtr;
            evalRadialFibShearSolid(ksr,ksi,omegar,omegai,cr,ci,r,t,
                                                            ut,vt,at,srt,vtr);

      //  printF("--DS-- UDKS:RadialFibShear: ut=%e vt=%e srt=%e at=%e vtr=%e\n",ut,vt,srt,at,srt,vtr);

            if( fabs(ut)>1e10 || fabs(vt)>1e10 || fabs(vtr)>1e10 || fabs(srt)>1e10 || ut!=ut || vt!=vt || srt!=srt )
            {
                printF("--SM-- UDKS:RadialFibShear: ERROR ut=%e vt=%e srt=%e\n",ut,vt,srt);
                OV_ABORT("error");
            }

            real u1 = -sinTheta*ut;
            real u2 =  cosTheta*ut;

            u(i1,i2,i3,u1c)=amp*u1;
            u(i1,i2,i3,u2c)=amp*u2;

      // velocities
            if( assignVelocity )
            {
                real v1 = -sinTheta*vt;
                real v2 =  cosTheta*vt;

                u(i1,i2,i3,v1c)=amp*v1;
                u(i1,i2,i3,v2c)=amp*v2;
            }
            
      // stresses
            if( assignStress )
            {
                real sigmart;
                sigmart = amp*mu*srt;

                u(i1,i2,i3,s11c)= -2.*cosTheta*sinTheta*sigmart;
                u(i1,i2,i3,s21c)=  (SQR(cosTheta)-SQR(sinTheta))*sigmart;
                u(i1,i2,i3,s12c)=  (SQR(cosTheta)-SQR(sinTheta))*sigmart;
                u(i1,i2,i3,s22c)=  2.*cosTheta*sinTheta*sigmart;
            }
            
        }
    // if( grid==1 )
    //   ::display(u,sPrintF("--SM-- UDKS u t=%.2e",t),"%8.2e ");

    }
    else if (userKnownSolution=="fibCartWave") 
    {
    // -- traveling wave solution for elastic solid and linearized fluid --
    // 
    // linearized fluid: 0 < x < L,      0 < y < H
    // solid reference:  0 < x < L,  -Hbar < y < 0

        const real & omegar = rpar[0];
        const real & omegai = rpar[1];
        const real & k      = rpar[2];
        const real & k1r    = rpar[3];
        const real & k1i    = rpar[4];
        const real & k2r    = rpar[5];
        const real & k2i    = rpar[6];
        const real & Ar     = rpar[7];
        const real & Ai     = rpar[8];
        const real & Br     = rpar[9];
        const real & Bi     = rpar[10];
        const real & amp    = rpar[11];
        const real & mu     = rpar[12];
        const real & rho    = rpar[13];
        const real & muBar  = rpar[14];
        const real & lambdaBar = rpar[15];
        const real & rhoBar = rpar[16];
        const real & H      = rpar[17];
        const real & HBar   = rpar[18];

        printF("--SM-- userDefinedKnownSolution: fibCartWave, t=%9.3e, "
                      "rhoBar=%9.3e, muBar=%9.3e\n",t,rho,mu);

        assert( numberOfTimeDerivatives==0 );

    // fill in the solution
        MappedGrid & mg = cg[grid];
        mg.update(MappedGrid::THEvertex | MappedGrid::THEcenter);

        const realArray & center = mg.center();
        RealArray & u = ua;

        int i1,i2,i3;
        FOR_3D(i1,i2,i3,I1,I2,I3)
            {
        // Reference coordinates:
                const real x= center(i1,i2,i3,0);
                const real y= center(i1,i2,i3,1);

                real u1Barr, u2Barr, v1Barr, v2Barr, s11Barr, s12Barr, s22Barr;

                evalFibCartWaveSolid(omegar, omegai, k, muBar,
                                                          rhoBar, lambdaBar, k1r, k1i,
                                                          k2r, k2i, amp, x, y, t, HBar, 
                                                          u1Barr, u2Barr, v1Barr, 
                                                          v2Barr, s11Barr, 
                                                          s12Barr, s22Barr);

                u(i1,i2,i3,u1c)=u1Barr;
                u(i1,i2,i3,u2c)=u2Barr;


        // velocities
                if( assignVelocity )
                    {
                        u(i1,i2,i3,v1c)=v1Barr;
                        u(i1,i2,i3,v2c)=v2Barr;
                    }
            
        // stresses
                if( assignStress )
                    {
                        u(i1,i2,i3,s11c)=s11Barr;
                        u(i1,i2,i3,s12c)=s12Barr;
                        u(i1,i2,i3,s21c)=s12Barr;
                        u(i1,i2,i3,s22c)=s22Barr;
                    }
            
            }



    }
    else if (userKnownSolution=="fibRadialWave") 
    {
      // -- traveling wave solution for elastic solid and linearized fluid --
      // 
      // linearized fluid:     0 < r < Rbar
      // solid reference:   Rbar < r < R
      //
      // Fluid solution: (hat vars)
      //  p       = b (r/R)^n 
      //  v_r     = a J_n(lambda r) / r + b n r^(n-1)/(mu lambda^2 R^n)
      //  v_theta = i a (J_n(lambda r) / r - lambda/n J_{n+1}(lambda r)) 
      //               + b i n r^(n-1) / (mu lambda^2 R^n)
      //
      // Solid solution: (hat vars)
      //  u_r     = d1 (n J_n(kp r) / r - kp J_{n+1}(kp r))
      //           +d2 (n Y_n(kp r) / r - kp Y_{n+1}(kp r))
      //           +in/r( d3 J_n(ks r) + d4 Y_n(ks r))
      //  u_theta =-d3 (n J_n(ks r) / r - ks J_{n+1}(ks r))
      //           -d4 (n Y_n(ks r) / r - ks Y_{n+1}(ks r))
      //           +in/r( d1 J_n(kp r) + d2 Y_n(kp r))
      //
      // kp = omega/cp, ks = omega/cs
      // lambda = sqrt(i omega rho / mu)
      //
      // All scaled by: amp exp(i(n theta - omega t))
      //

            const real & omegar = rpar[0];
            const real & omegai = rpar[1];
            const real & n      = rpar[2];
            const real & d1r    = rpar[3];
            const real & d1i    = rpar[4];
            const real & d2r    = rpar[5];
            const real & d2i    = rpar[6];
            const real & d3r    = rpar[7];
            const real & d3i    = rpar[8];
            const real & d4r    = rpar[9];
            const real & d4i    = rpar[10];
            const real & ar     = rpar[11];
            const real & ai     = rpar[12];
            const real & br     = rpar[13];
            const real & bi     = rpar[14];
            const real & nu     = rpar[15];
            const real & cp     = rpar[16];
            const real & cs     = rpar[17];
            const real & amp    = rpar[18];

        printF("--SM-- userDefinedKnownSolution: fibRadialWave, t=%9.3e, "
                      "rhoBar=%9.3e, muBar=%9.3e\n",t,rho,mu);

        assert( numberOfTimeDerivatives==0 );

    // fill in the solution
        MappedGrid & mg = cg[grid];
        mg.update(MappedGrid::THEvertex | MappedGrid::THEcenter);

        const realArray & center = mg.center();
        RealArray & u = ua;

        const real eps = 10.*REAL_EPSILON;
        int i1,i2,i3;
        FOR_3D(i1,i2,i3,I1,I2,I3)
            {
        // Reference coordinates:
                const real x= center(i1,i2,i3,0);
                const real y= center(i1,i2,i3,1);
                const real r= sqrt(SQR(x)+SQR(y));

                real ct,st;
        // compute trig functions
                real theta = atan2(y,x);
                if( r>eps ) {
                    ct=x/r; st=y/r;
                } else {
                    ct=1.; st=0.;  // at the origin we just pick an angle, should not matter
                }

                real ur,ut,vr,vt,ar,at,srr,srt,stt,sdrr,sdrt,sdtt;
                evalFibRadialWaveSolid(omegar,omegai,d1r,d1i,d2r,d2i,d3r,d3i,d4r,d4i,
                                                              n,rho,mu,lambda,r,theta,t,
                                                              ur,ut,vr,vt,ar,at,srr,srt,stt,sdrr,sdrt,sdtt);

                real u1 = amp*(ct*ur-st*ut);
                real u2 = amp*(st*ur+ct*ut);
                u(i1,i2,i3,u1c)=u1;
                u(i1,i2,i3,u2c)=u2;

        // velocities
                if( assignVelocity )
                    {
                        real v1 = amp*(ct*vr-st*vt);
                        real v2 = amp*(st*vr+ct*vt);
                        u(i1,i2,i3,v1c) = v1;
                        u(i1,i2,i3,v2c) = v2;
            // printF("(%d,%d,%d), amp=%f, t=%f\n",i1,i2,i3,amp,t);
            // printF("  r=%f, theta=%f\n",r,theta);
            // printF("  ur=%e, ut=%e\n",ur,ut);
            // printF("  u1=%e, u2=%e\n",u1,u2);
            // printF("  vr=%e, vt=%e\n",vr,vt);
            // printF("  v1=%e, v2=%e\n",v1,v2);
            // printF("  omega = %e + i %e\n",omegar,omegai);

                    }
            
        // stresses
                if( assignStress )
                    {
                        u(i1,i2,i3,s11c)=amp*(srr*SQR(ct)
                                                                    -2.0*srt*ct*st
                                                                    +stt*SQR(st));
                        u(i1,i2,i3,s12c)=amp*(srr*ct*st
                                                                    +srt*(2.0*SQR(ct)-1)
                                                                    -stt*ct*st);
                        u(i1,i2,i3,s21c)=u(i1,i2,i3,s12c);
                        u(i1,i2,i3,s22c)=amp*(srr*SQR(st)
                                                                    +2.0*srt*ct*st
                                                                    +stt*SQR(ct));
            // printF("(%d,%d,%d), amp=%f, t=%f\n",i1,i2,i3,amp,t);
            // printF("  r=%f, theta=%f\n",r,theta);
            // printF("  srr=%e, srt=%e, stt=%e\n",srr,srt,stt);

                    }
            }

    }

    else 
    {
    // look for a solution in the base class
        Parameters::getUserDefinedKnownSolution( t, cg, grid, ua, I1,I2,I3 );



    }
    
  // else
  // {
  //   printF("getUserDefinedKnownSolution:ERROR: unknown value for userDefinedKnownSolution=%s\n",
  //       (const char*)userKnownSolution);
  //   OV_ABORT("ERROR");
  // }
    
    return 0;
}



int SmParameters::
updateUserDefinedKnownSolution(GenericGraphicsInterface & gi, CompositeGrid & cg)
// ==========================================================================================
/// \brief This function is called to set the user defined know solution.
/// 
/// \return   >0 : known solution was chosen, 0 : no known solution was chosen
///
// ==========================================================================================
{
  // Make  dbase.get<real >("a") sub-directory in the data-base to store variables used here
    if( ! dbase.get<DataBase >("modelData").has_key("userDefinedKnownSolutionData") )
        dbase.get<DataBase >("modelData").put<DataBase>("userDefinedKnownSolutionData");
    DataBase & db =  dbase.get<DataBase >("modelData").get<DataBase>("userDefinedKnownSolutionData");

    if( !db.has_key("userKnownSolution") )
    {
        db.put<aString>("userKnownSolution");
        db.get<aString>("userKnownSolution")="unknownSolution";
        
        db.put<real[20]>("rpar");
        db.put<int[20]>("ipar");
    }
    aString & userKnownSolution = db.get<aString>("userKnownSolution");
    real *rpar = db.get<real[20]>("rpar");
    int *ipar = db.get<int[20]>("ipar");


    const aString menu[]=
        {
            "no known solution",
            "rotating disk",  // for cgsm SVK model
            "linear beam exact solution",
            "incompressible rectangle eigenmodes",
            "incompressible surface wave",
            "incompressible annulus solution",
            "incompressible disk solution",
            "gaussian plane wave",
            "incompressible hollow cylinder",
            "incompressible hollow sphere",
            "incompressible solid sphere",
            "choose a common known solution",
            "done",
            ""
        }; 

    gi.appendToTheDefaultPrompt("userDefinedKnownSolution>");
    aString answer;
    for( ;; ) 
    {

        int response=gi.getMenuItem(menu,answer,"Choose a known solution");
        
        if( answer=="done" || answer=="exit" )
        {
            break;
        }
        else if( answer=="choose a common known solution" )
        {
      // Look for a known solution from the base class (in common/src)
            Parameters::updateUserDefinedKnownSolution(gi,cg);
        }
        else if( answer=="no known solution" )
        {
            userKnownSolution="unknownSolution";
        }
        else if( answer=="exact solution from a file" )
        {
            userKnownSolution="exactSolutionFromAFile";
            dbase.get<bool>("knownSolutionIsTimeDependent")=false;  // known solution does NOT depend on time

            printF("The exact solution can be defined by a solution in a show file (e.g. from a fine grid solution)\n");
            
            gi.inputString(answer,"Enter the the name of the file holding the exact solution");

      // sScanF(answer,"%e %e",&rpar[0],&rpar[1]);
      // printF("forced piston: mass=%e, height=%e\n",rpar[0],rpar[1]);

        }
        else if( answer=="rotating disk" )
        {
            userKnownSolution="rotatingDisk";
            dbase.get<bool>("knownSolutionIsTimeDependent")=true;  // known solution IS time dependent

            if( !db.has_key("tDisk") )
            { // Create parameters for the rotating disk solution
                db.put<real>("tDisk");
                db.put<int>("numberOfGridPointsDisk");
                db.put<real>("omegaDisk");
                db.put<real>("innerRadiusDisk");
                db.put<real>("outerRadiusDisk");
                db.put<RealArray>("uDisk"); // exact solution is stored here
            }
            

            real & tDisk = db.get<real>("tDisk");
            int & numberOfGridPoints = db.get<int>("numberOfGridPointsDisk");
            real & omega = db.get<real>("omegaDisk");
            real & innerRadius = db.get<real>("innerRadiusDisk");
            real & outerRadius = db.get<real>("outerRadiusDisk");
            RealArray & uDisk = db.get<RealArray>("uDisk"); // exact solution is stored here

      // Defaults:
            tDisk=-1.;  // this will cause the solution to be computed 
            numberOfGridPoints=101;
            omega=.5;
            innerRadius=0.;
            outerRadius=1.;

      // Prompt for input:
            printF("--- The rotating disk exact solution requires: ---\n"
                          " n : number of points to use when computing the exact solution\n"
                          " omega : rotation rate\n"
                          " ra,rb : radial bounds\n");
            gi.inputString(answer,"Enter n,omega,ra,rb for the exact solution");
            sScanF(answer,"%i %e %e %e",&numberOfGridPoints,&omega,&innerRadius,&outerRadius);

            printF("rotatingDisk: setting n=%i, omega=%9.3e, ra=%9.3e, rb=%9.3e\n",
                          numberOfGridPoints,omega,innerRadius,outerRadius);

//       // We need to keep a FlowSolutions object around
//       db.put<FlowSolutions*>("flowSolutions",NULL);

//       db.get<FlowSolutions*>("flowSolutions")=new FlowSolutions;

//       FlowSolutions & flowSolutions = *db.get<FlowSolutions*>("flowSolutions");

        }
        
        else if( answer=="linear beam exact solution" ) 
        {

            userKnownSolution="linearBeamExactSolution";
            dbase.get<bool>("knownSolutionIsTimeDependent")=true;  // known solution IS time dependent 
            double omega;
            gi.inputString(answer,"Enter omega");
            sScanF(answer,"%e",&omega);
        }

        else if( answer=="incompressible surface wave" )
        {
            userKnownSolution="incompressibleSurfaceWave";
            dbase.get<bool>("knownSolutionIsTimeDependent")=true;  // known solution IS time dependent

            printF("--- Incompressible surface wave ---\n"
                          " This is an exact solution to the incompressible linear elastic equations.\n"
                          " The solution is defined on a rectange and is periodic in y\n"
                          " The are different cases defined:\n"
                          "     Case 1: traction-displacement BC: surface wave as speed cs=omega/k, region=[0,1]x[0,Ly] Ly is given.\n"
                          "     Case 2: traction-traction BC\n"
                          "     Case 3: displacement-displacement BC.\n"
                          );      

            if( !db.has_key("iswCase") )
            { // Create parameters for the incompressible surface wave
                db.put<int>("iswCase")=1;
                db.put<int>("mPeriod")=1;
            }
            int & iswCase = db.get<int>("iswCase");
            int & mPeriod = db.get<int>("mPeriod");
            gi.inputString(answer,"Enter iswCase and mPeriod");
            sScanF(answer,"%i %i",&iswCase,&mPeriod);

        }

        else if( answer=="incompressible annulus solution" )
        {
            userKnownSolution="incompressibleAnnulusSolution";
            dbase.get<bool>("knownSolutionIsTimeDependent")=true;  // known solution IS time dependent

            printF("--- Incompressible Annulus Solution ---\n"
                          " This is an exact solution to the incompressible linear elastic equations.\n"
                          " The solution is defined on an annulus and involves:\n"
                          "         Jn(lambda_mn r) : n = order of Bessel, m=1,2,3, root number\n"
                          " The are different cases defined:\n"
                          "     Case 1: displacement + displacement (BC on inner and outer boundaries)\n"
                          "     Case 2: traction + displacement\n"
                          "     Case 3: traction + traction\n"
                          );      

            if( !db.has_key("iswCase") )
            { // Create parameters for the incompressible annulus solution
                db.put<int>("iswCase")=1;
            }
            int & iswCase = db.get<int>("iswCase");
            int n=1, m=1;
            gi.inputString(answer,"Enter iswCase, n,m");
            sScanF(answer,"%i %i %i",&iswCase,&n,&m);
            printF("Setting iswCase=%d, n=%d, m=%d ( J_n(lamnda_mn r)\n",iswCase,n,m);
            ipar[0]=n;
            ipar[1]=m;




        }

        else if( answer=="incompressible disk solution" )
        {
            userKnownSolution="incompressibleDiskSolution";
            dbase.get<bool>("knownSolutionIsTimeDependent")=true;  // known solution IS time dependent

            printF("--- Incompressible Disk Solution ---\n"
                          " This is an exact solution to the incompressible linear elastic equations.\n"
                          " The solution is defined on the unit disk and involves:\n"
                          "         Jn(lambda_mn r) : n = order of Bessel, m=1,2,3, root number\n"             
                          " The are different cases defined:\n"
                          "     Case 1: displacement BC\n"
                          "     Case 2: traction BC\n"
                          );      

            if( !db.has_key("iswCase") )
            { // Create parameters for the incompressible annulus solution
                db.put<int>("iswCase")=1;
            }
            int & iswCase = db.get<int>("iswCase");
            int n=1, m=1;
            gi.inputString(answer,"Enter iswCase, n,m");
            sScanF(answer,"%i %i %i",&iswCase,&n,&m);
            printF("Setting iswCase=%d, n=%d, m=%d ( J_n(lamnda_mn r)\n",iswCase,n,m);
            ipar[0]=n;
            ipar[1]=m;      


        }    
        else if( answer=="gaussian plane wave" )
        {
            userKnownSolution="gaussianPlaneWave";
            dbase.get<bool>("knownSolutionIsTimeDependent")=true;  // known solution IS time dependent

            printF("--- Gaussian plane wave ---\n"
                          " This is an exact solution to the incompressible linear elastic equations.\n"
                          " The solution is \n"
                          "       uv = av exp( -beta * xi^2 ) cos( 2*pi*k0*xi )\n"
                          "    xi = kv.(xv-xv0) - omega*t \n"
                          "    kv = ( kx,ky,kz )\n"
                          "    xv0 = (x0,y0,z0)\n"
                          "    k0 = modulation frequency\n"
                          "    av = (ax,ay,az) if ax=ay=az=0 then av is computed automatically so that av.kv =0\n"
                          "    c = shear wave speed\n");

            Real kx=1., ky=0, kz=0., k0=0., beta=10., ax=0., ay=1., az=0., x0=0., y0=0., z0=0.; 
            gi.inputString(answer,"Enter kx,ky,kz,k0, beta, x0,y0,z0,ax,ay,az");
            sScanF(answer,"%e %e %e %e %e %e %e %e %e %e %e",&kx, &ky, &kz, &k0, &beta, &x0,&y0,&z0, &ax, &ay, &az);

            k0*=twoPi;

            if( fabs(ax)+fabs(ay)+fabs(az) == 0. )
            {
                const real kNorm = sqrt(kx*kx+ky*ky+kz*kz);
        // const real cs = sqrt(mu/rho);
        // const real cc = cs*kNorm;
                if( fabs(kx)+fabs(ky)>0. )
                {
                    ax = -ky/kNorm;
                    ay =  kx/kNorm;
                    az =0.;
                }
                else
                {
                    ay = -kz/kNorm;
                    az =  ky/kNorm;
                    ax =0.;
                }
            }      

            printF("Setting (kx,ky,kz)=(%g,%g,%g), k0=%g, beta=%g, (x0,y0,z0)=(%g,%g,%g) (ax,ay,az)=(%g,%g,%g) \n",kx/twoPi,ky/twoPi,kz/twoPi,k0,beta,x0,y0,z0,ax,ay,az);
            rpar[ 0]=kx;
            rpar[ 1]=ky;      
            rpar[ 2]=kz;      
            rpar[ 3]=k0;      
            rpar[ 4]=beta;      
            rpar[ 5]=ax;      
            rpar[ 6]=ay;      
            rpar[ 7]=az;  
            rpar[ 8]=x0;      
            rpar[ 9]=y0;      
            rpar[10]=z0;              


        }        
        else if( answer=="incompressible hollow cylinder" ||
                          answer=="incompressible solid cylinder"  )
        {
            if( answer=="incompressible hollow cylinder" )
                userKnownSolution="incompressibleHollowCylinder";
            else
                userKnownSolution="incompressibleSolidCylinder";
            dbase.get<bool>("knownSolutionIsTimeDependent")=true;  // known solution IS time dependent

            printF("--- Incompressible hollow or solid cylinder ---\n"
                          " This is an exact solution to the incompressible linear elastic equations.\n");

            aString caseName="hollowCylinderDD"; 
            gi.inputString(caseName,"Enter the caseName (hollowCylinderDD, hollowCylinderTT, solidCylinderD, solidCylinderT )");
            if( !dbase.has_key("SmCylinderExactSolution") )
            {
                dbase.put<SmCylinderExactSolution>("cylExact");
            }

      // // int & iswCase = db.get<int>("iswCase");
      // int iswCase=0; // FINISH ME 
      // int n=1, m=1;
      // gi.inputString(answer,"Enter iswCase, n,m");
      // sScanF(answer,"%i %i %i",&iswCase,&n,&m);
      // printF("Setting iswCase=%d, n=%d, m=%d\n",iswCase,n,m);
      // ipar[0]=n;
      // ipar[1]=m;


            SmCylinderExactSolution & cylExact = dbase.get<SmCylinderExactSolution>("cylExact");
            cylExact.initialize( cg, caseName );

        }
        else if( answer=="incompressible hollow sphere" ||
                          answer=="incompressible solid sphere"  )
        {
            if( answer=="incompressible hollow sphere" )
                userKnownSolution="incompressibleHollowSphere";
            else
                userKnownSolution="incompressibleSolidSphere";
            
            dbase.get<bool>("knownSolutionIsTimeDependent")=true;  // known solution IS time dependent

            printF("--- Incompressible hollow or solid sphere ---\n"
                          " This is an exact solution to the incompressible linear elastic equations.\n");

            aString caseName="solidSphereD"; 
            gi.inputString(caseName,"Enter the caseName (solidSphereD, solidSphereT )");
            if( !dbase.has_key("SmSphereExactSolution") )
            {
                dbase.put<SmSphereExactSolution>("sphereExact");
            }

      // // int & iswCase = db.get<int>("iswCase");
      // int iswCase=0; // FINISH ME 
      // int n=1, m=1;
      // gi.inputString(answer,"Enter iswCase, n,m");
      // sScanF(answer,"%i %i %i",&iswCase,&n,&m);
      // printF("Setting iswCase=%d, n=%d, m=%d\n",iswCase,n,m);
      // ipar[0]=n;
      // ipar[1]=m;


            SmSphereExactSolution & sphereExact = dbase.get<SmSphereExactSolution>("sphereExact");
            sphereExact.initialize( cg, caseName );


        } 
        else if( answer=="incompressible rectangle eigenmodes" )
        {
            userKnownSolution="incompressibleRectangleEigenmode";

            
            dbase.get<bool>("knownSolutionIsTimeDependent")=true;  // known solution IS time dependent

            printF("--- Eigenmodes of an incompressible solid on a rectangular domain ---\n"
                          " This is an exact solution to the incompressible linear elastic equations.\n"
                          "       u1 = uHat(x) * exp( 2*pi*my*y ) * exp( i omega t).\n"
                          "       omega = omega(mx,my), mx=1,2,...,  my=1,2,...\n"
                          " caseName:\n"
                          "   rectangleDD : periodic in y, BC=displacement on left and right.\n"
                          "   rectangleDD : periodic in y, BC=traction on left and BC=dispplacement on right.\n"
                          "   rectangleTT : periodic in y, BC=traction on left and right.\n"
                          );

            aString caseName="rectangleTT"; 
            gi.inputString(caseName,"Enter the caseName (rectangleDD, rectangleTD, rectangleTT )");
            printF("Setting caseName=[%s]\n",(const char*)caseName);
            if( !dbase.has_key("SmRectangleExactSolution") )
            {
                dbase.put<SmRectangleExactSolution>("rectangleExact");
            }

            int mx=2, my=1, mz=1;
            gi.inputString(answer,"Enter mode numbers: mx,my,mz");
            sScanF(answer,"%i %i %i",&mx,&my,&mz);

            SmRectangleExactSolution & rectangleExact = dbase.get<SmRectangleExactSolution>("rectangleExact");
            rectangleExact.initialize( cg, caseName );
            rectangleExact.setParameter( "mx",mx );
            rectangleExact.setParameter( "my",my );
            rectangleExact.setParameter( "mz",mz );

        }    

        else
        {
            printF("unknown response=[%s]\n",(const char*)answer);
            gi.stopReadingCommandFile();
        }
        
    }

    gi.unAppendTheDefaultPrompt();
    bool knownSolutionChosen = userKnownSolution!="unknownSolution";
    return knownSolutionChosen;
}
