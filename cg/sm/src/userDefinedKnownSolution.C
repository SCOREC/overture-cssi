// This file automatically generated from userDefinedKnownSolution.bC with bpp.
#include "SmParameters.h"
#include "FlowSolutions.h"
#include "GenericGraphicsInterface.h"
#include "FluidPiston.h"
#include "PistonMotion.h"
#include "ParallelUtility.h"
#include "TimeFunction.h"
#include "BeamModel.h"

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

#define FOR_3D(i1,i2,i3,I1,I2,I3) int I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase();  int I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(); for(i3=I3Base; i3<=I3Bound; i3++) for(i2=I2Base; i2<=I2Bound; i2++) for(i1=I1Base; i1<=I1Bound; i1++)

// Macro to get the vertex array
#define GET_VERTEX_ARRAY(x)                                     mg.update(MappedGrid::THEvertex | MappedGrid::THEcenter);       OV_GET_SERIAL_ARRAY_CONST(real,mg.vertex(),x);


// --------------------------------------------------------------------------------
// Assign the incompressible annulus solution
// --------------------------------------------------------------------------------


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


int SmParameters::
getUserDefinedKnownSolution(real t, CompositeGrid & cg, int grid, RealArray & ua, 
                                                        const Index & I1, const Index &I2, const Index &I3, 
                                                        int numberOfTimeDerivatives /* = 0 */  )
// ==========================================================================================
///  \brief Evaluate a user defined known solution.
// ==========================================================================================
{
    MappedGrid & mg = cg[grid];

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

    const int u1c = dbase.get<int >("uc");
    const int u2c = dbase.get<int >("vc");

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
        else if( iswCase==2 )
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

        int & iswCase = db.get<int>("iswCase");

        if( t<3*dt )
            printF("UDKS: evaluate the incompressible annulus solution: iswCase=%i, t=%9.3e\n",iswCase,t);


        Real C4 =1.; // defines the scale

        if( iswCase==1 )
        {
      // files generated from AMP/ism/maple/DiskAnnulusSolutionWDH.mw
            C4 = 1./10.;  // Make |u|=.5, |p|=5
// ------------------------------------------------------------------------------------------------ 
//    AnnulusDispDispBCs: Exact solution for incompressible linear elasticity 
// Annulus of outer radius Rb=2, inner radius Ra = 1; r=1 Displacement BC; r=2 Displacement BC. 
// File Written by Disk_Annulus_ExactSoln.mw (maple) 
// Set C4 to scale the solution 
// ------------------------------------------------------------------------------------------------ 
Real rho,c,mu,omegaRoot,betaRoot,Pi,k,C1,C2,C3,ra,rb; 
ra=.5;
rb=1;
rho=1;
mu=1;
betaRoot=12.5105206473408546981549;
k=1;
C1 = 0.192992035693309510535441e2 * C4;
C2 = -0.304935551313522923733292e2 * C4;
C3 = 0.992579994490862297289017e0 * C4;
omegaRoot = betaRoot * sqrt(mu / rho);
Pi=3.14159265358979323846264;
                Real ur,vr,pr; 
                int i1,i2,i3;
                FOR_3D(i1,i2,i3,I1,I2,I3)
                {
          // Reference coordinates:
                    Real x = vertex(i1,i2,i3,0);
                    Real y = vertex(i1,i2,i3,1);
                    Real r = sqrt( x*x + y*y );
                    Real theta = atan2(y,x); 
                    Real cosTheta = x/r;
                    Real sinTheta = y/r;
                        ur = 0.1e1 / r * jn(1, 0.125105206473408546981549e2 * r) * C4 * cos(t * omegaRoot) * cos(theta) + 0.1e1 / r * jn(1, 0.125105206473408546981549e2 * r) * C4 * sin(t * omegaRoot) * sin(theta) + 0.992579994490862297289017e0 / r * yn(1, 0.125105206473408546981549e2 * r) * C4 * cos(t * omegaRoot) * cos(theta) + 0.992579994490862297289017e0 / r * yn(1, 0.125105206473408546981549e2 * r) * C4 * sin(t * omegaRoot) * sin(theta) + 0.123307251957620161779820e0 * C4 * cos(t * omegaRoot) * cos(theta) + 0.123307251957620161779820e0 * C4 * sin(t * omegaRoot) * sin(theta) + 0.194830655687808236749749e0 * pow(r, -0.2e1) * C4 * cos(t * omegaRoot) * cos(theta) + 0.194830655687808236749749e0 * pow(r, -0.2e1) * C4 * sin(t * omegaRoot) * sin(theta);
                        vr = -0.100000000000000000000000e1 / r * jn(1, 0.125105206473408546981549e2 * r) * C4 * sin(t * omegaRoot) * cos(theta) + 0.100000000000000000000000e1 / r * jn(1, 0.125105206473408546981549e2 * r) * C4 * cos(t * omegaRoot) * sin(theta) + 0.124176925152154045774849e2 * yn(0, 0.125105206473408546981549e2 * r) * C4 * sin(t * omegaRoot) * cos(theta) - 0.124176925152154045774849e2 * yn(0, 0.125105206473408546981549e2 * r) * C4 * cos(t * omegaRoot) * sin(theta) + 0.125105206473408546981549e2 * jn(0, 0.125105206473408546981549e2 * r) * C4 * sin(t * omegaRoot) * cos(theta) - 0.125105206473408546981549e2 * jn(0, 0.125105206473408546981549e2 * r) * C4 * cos(t * omegaRoot) * sin(theta) - 0.992579994490862297289018e0 / r * yn(1, 0.125105206473408546981549e2 * r) * C4 * sin(t * omegaRoot) * cos(theta) + 0.992579994490862297289018e0 / r * yn(1, 0.125105206473408546981549e2 * r) * C4 * cos(t * omegaRoot) * sin(theta) + 0.123307251957620161779820e0 * C4 * sin(t * omegaRoot) * cos(theta) - 0.123307251957620161779820e0 * C4 * cos(t * omegaRoot) * sin(theta) - 0.194830655687808236749749e0 * pow(r, -0.2e1) * C4 * sin(t * omegaRoot) * cos(theta) + 0.194830655687808236749749e0 * pow(r, -0.2e1) * C4 * cos(t * omegaRoot) * sin(theta);
                        pr = 0.192992035693309510535441e2 * C4 * mu * pow(r, k) * cos(k * theta) * cos(t * omegaRoot) + 0.192992035693309510535441e2 * C4 * mu * pow(r, k) * sin(k * theta) * sin(t * omegaRoot) - 0.304935551313522923733292e2 * C4 * mu / pow(r, k) * cos(k * theta) * cos(t * omegaRoot) - 0.304935551313522923733292e2 * C4 * mu / pow(r, k) * sin(k * theta) * sin(t * omegaRoot);
          // (u1,u2) = ur*rHat + vr*thetaHat 
                    u(i1,i2,i3,u1c) = ur*cosTheta - vr*sinTheta;
                    u(i1,i2,i3,u2c) = ur*sinTheta + vr*cosTheta;
                    u(i1,i2,i3,pc ) = pr;
                }
        }
        else if( iswCase==2 )
        {
// ------------------------------------------------------------------------------------------------ 
//    AnnulusTracDispBCs: Exact solution for incompressible linear elasticity 
// Annulus of outer radius Rb=2, inner radius Ra = 1; r=1 Traction BC; r=2 Displacement BC. 
// File Written by Disk_Annulus_ExactSoln.mw (maple) 
// Set C4 to scale the solution 
// ------------------------------------------------------------------------------------------------ 
Real rho,c,mu,omegaRoot,betaRoot,Pi,k,C1,C2,C3,ra,rb; 
ra=.5;
rb=1.;
rho=1;
mu=1;
betaRoot=2.70180213338953294069697;
k=1;
C1 = -0.421761043357083014501791e1 * C4;
C2 = 0.187798693127866841300838e1 * C4;
C3 = -0.527973122762447522209045e0 * C4;
omegaRoot = betaRoot * sqrt(mu / rho);
Pi=3.14159265358979323846264;
                Real ur,vr,pr; 
                int i1,i2,i3;
                FOR_3D(i1,i2,i3,I1,I2,I3)
                {
          // Reference coordinates:
                    Real x = vertex(i1,i2,i3,0);
                    Real y = vertex(i1,i2,i3,1);
                    Real r = sqrt( x*x + y*y );
                    Real theta = atan2(y,x); 
                    Real cosTheta = x/r;
                    Real sinTheta = y/r;
                        ur = 0.1e1 / r * jn(1, 0.270180213338953294069697e1 * r) * C4 * cos(t * omegaRoot) * cos(theta) + 0.1e1 / r * jn(1, 0.270180213338953294069697e1 * r) * C4 * sin(t * omegaRoot) * sin(theta) + 0.185886219453615487107146e2 * C4 * mu / (-0.458941557224063735118251e2 * mu - 0.98912835974785841653545e0) * cos(t * omegaRoot) * cos(theta) + 0.185886219453615487107146e2 * C4 * mu / (-0.458941557224063735118251e2 * mu - 0.98912835974785841653545e0) * sin(t * omegaRoot) * sin(theta) - 0.652706113525886387519789e1 * pow(r, -0.2e1) * C4 / (-0.458941557224063735118251e2 * mu - 0.98912835974785841653545e0) * cos(t * omegaRoot) * cos(theta) - 0.652706113525886387519789e1 * pow(r, -0.2e1) * C4 / (-0.458941557224063735118251e2 * mu - 0.98912835974785841653545e0) * sin(t * omegaRoot) * sin(theta) + 0.429478687316422375462374e2 / r * yn(1, 0.270180213338953294069697e1 * r) * C4 / (-0.458941557224063735118251e2 * mu - 0.98912835974785841653545e0) * mu * cos(t * omegaRoot) * cos(theta) + 0.429478687316422375462374e2 / r * yn(1, 0.270180213338953294069697e1 * r) * C4 / (-0.458941557224063735118251e2 * mu - 0.98912835974785841653545e0) * mu * sin(t * omegaRoot) * sin(theta) - 0.181947548294283194480474e2 / r * yn(1, 0.270180213338953294069697e1 * r) * C4 / (-0.458941557224063735118251e2 * mu - 0.98912835974785841653545e0) * cos(t * omegaRoot) * cos(theta) - 0.181947548294283194480474e2 / r * yn(1, 0.270180213338953294069697e1 * r) * C4 / (-0.458941557224063735118251e2 * mu - 0.98912835974785841653545e0) * sin(t * omegaRoot) * sin(theta) + 0.185886219453615487107146e2 * mu * pow(r, -0.2e1) * C4 / (-0.458941557224063735118251e2 * mu - 0.98912835974785841653545e0) * cos(t * omegaRoot) * cos(theta) + 0.185886219453615487107146e2 * mu * pow(r, -0.2e1) * C4 / (-0.458941557224063735118251e2 * mu - 0.98912835974785841653545e0) * sin(t * omegaRoot) * sin(theta) + 0.849940719401138113113085e1 * C4 / (-0.458941557224063735118251e2 * mu - 0.98912835974785841653545e0) * cos(t * omegaRoot) * cos(theta) + 0.849940719401138113113085e1 * C4 / (-0.458941557224063735118251e2 * mu - 0.98912835974785841653545e0) * sin(t * omegaRoot) * sin(theta);
                        vr = 0.116036643363684611599021e3 * mu * C4 / (-0.458941557224063735118251e2 * mu - 0.98912835974785841653545e0) * yn(0, 0.270180213338953294069697e1 * r) * sin(t * omegaRoot) * cos(theta) - 0.491586274146489410093073e2 * C4 / (-0.458941557224063735118251e2 * mu - 0.98912835974785841653545e0) * yn(0, 0.270180213338953294069697e1 * r) * sin(t * omegaRoot) * cos(theta) + 0.491586274146489410093073e2 * C4 / (-0.458941557224063735118251e2 * mu - 0.98912835974785841653545e0) * yn(0, 0.270180213338953294069697e1 * r) * cos(t * omegaRoot) * sin(theta) - 0.116036643363684611599021e3 * mu * C4 / (-0.458941557224063735118251e2 * mu - 0.98912835974785841653545e0) * yn(0, 0.270180213338953294069697e1 * r) * cos(t * omegaRoot) * sin(theta) - 0.999999999999999999999999e0 / r * jn(1, 0.270180213338953294069697e1 * r) * C4 * sin(t * omegaRoot) * cos(theta) + 0.999999999999999999999999e0 / r * jn(1, 0.270180213338953294069697e1 * r) * C4 * cos(t * omegaRoot) * sin(theta) + 0.185886219453615487107146e2 * C4 * mu / (-0.458941557224063735118251e2 * mu - 0.98912835974785841653545e0) * sin(t * omegaRoot) * cos(theta) - 0.185886219453615487107146e2 * C4 * mu / (-0.458941557224063735118251e2 * mu - 0.98912835974785841653545e0) * cos(t * omegaRoot) * sin(theta) + 0.652706113525886387519789e1 * pow(r, -0.2e1) * C4 / (-0.458941557224063735118251e2 * mu - 0.98912835974785841653545e0) * sin(t * omegaRoot) * cos(theta) - 0.652706113525886387519789e1 * pow(r, -0.2e1) * C4 / (-0.458941557224063735118251e2 * mu - 0.98912835974785841653545e0) * cos(t * omegaRoot) * sin(theta) - 0.429478687316422375462373e2 / r * yn(1, 0.270180213338953294069697e1 * r) * C4 / (-0.458941557224063735118251e2 * mu - 0.98912835974785841653545e0) * mu * sin(t * omegaRoot) * cos(theta) + 0.429478687316422375462373e2 / r * yn(1, 0.270180213338953294069697e1 * r) * C4 / (-0.458941557224063735118251e2 * mu - 0.98912835974785841653545e0) * mu * cos(t * omegaRoot) * sin(theta) + 0.270180213338953294069696e1 * C4 * jn(0, 0.270180213338953294069697e1 * r) * sin(t * omegaRoot) * cos(theta) - 0.270180213338953294069696e1 * C4 * jn(0, 0.270180213338953294069697e1 * r) * cos(t * omegaRoot) * sin(theta) + 0.181947548294283194480474e2 / r * yn(1, 0.270180213338953294069697e1 * r) * C4 / (-0.458941557224063735118251e2 * mu - 0.98912835974785841653545e0) * sin(t * omegaRoot) * cos(theta) - 0.181947548294283194480474e2 / r * yn(1, 0.270180213338953294069697e1 * r) * C4 / (-0.458941557224063735118251e2 * mu - 0.98912835974785841653545e0) * cos(t * omegaRoot) * sin(theta) - 0.185886219453615487107146e2 * mu * pow(r, -0.2e1) * C4 / (-0.458941557224063735118251e2 * mu - 0.98912835974785841653545e0) * sin(t * omegaRoot) * cos(theta) + 0.185886219453615487107146e2 * mu * pow(r, -0.2e1) * C4 / (-0.458941557224063735118251e2 * mu - 0.98912835974785841653545e0) * cos(t * omegaRoot) * sin(theta) + 0.849940719401138113113085e1 * C4 / (-0.458941557224063735118251e2 * mu - 0.98912835974785841653545e0) * sin(t * omegaRoot) * cos(theta) - 0.849940719401138113113085e1 * C4 / (-0.458941557224063735118251e2 * mu - 0.98912835974785841653545e0) * cos(t * omegaRoot) * sin(theta);
                        pr = 0.135692009903544734172044e3 * C4 * mu * mu / (-0.458941557224063735118251e2 * mu - 0.98912835974785841653545e0) * pow(r, k) * cos(k * theta) * cos(t * omegaRoot) + 0.135692009903544734172044e3 * C4 * mu * mu / (-0.458941557224063735118251e2 * mu - 0.98912835974785841653545e0) * pow(r, k) * sin(k * theta) * sin(t * omegaRoot) + 0.620434182014141753756432e2 * C4 * mu / (-0.458941557224063735118251e2 * mu - 0.98912835974785841653545e0) * pow(r, k) * cos(k * theta) * cos(t * omegaRoot) + 0.620434182014141753756432e2 * C4 * mu / (-0.458941557224063735118251e2 * mu - 0.98912835974785841653545e0) * pow(r, k) * sin(k * theta) * sin(t * omegaRoot) - 0.135692009903544734172044e3 * C4 * mu * mu / (-0.458941557224063735118251e2 * mu - 0.98912835974785841653545e0) / pow(r, k) * cos(k * theta) * cos(t * omegaRoot) - 0.135692009903544734172044e3 * C4 * mu * mu / (-0.458941557224063735118251e2 * mu - 0.98912835974785841653545e0) / pow(r, k) * sin(k * theta) * sin(t * omegaRoot) + 0.476458151018338659124449e2 * C4 * mu / (-0.458941557224063735118251e2 * mu - 0.98912835974785841653545e0) / pow(r, k) * cos(k * theta) * cos(t * omegaRoot) + 0.476458151018338659124449e2 * C4 * mu / (-0.458941557224063735118251e2 * mu - 0.98912835974785841653545e0) / pow(r, k) * sin(k * theta) * sin(t * omegaRoot);
          // (u1,u2) = ur*rHat + vr*thetaHat 
                    u(i1,i2,i3,u1c) = ur*cosTheta - vr*sinTheta;
                    u(i1,i2,i3,u2c) = ur*sinTheta + vr*cosTheta;
                    u(i1,i2,i3,pc ) = pr;
                }
        }
        else
        {
            C4 = 1./4.; // make u size 1
// ------------------------------------------------------------------------------------------------ 
//    AnnulusTracTracBCs: Exact solution for incompressible linear elasticity 
// Annulus of outer radius Rb=2, inner radius Ra = 1; r=1 Traction BC; r=2 Traction BC. 
// File Written by Disk_Annulus_ExactSoln.mw (maple) 
// Set C4 to scale the solution 
// ------------------------------------------------------------------------------------------------ 
Real rho,c,mu,omegaRoot,betaRoot,Pi,k,C1,C2,C3,ra,rb; 
ra=.5;
rb=1;
rho=1;
mu=1;
betaRoot=13.0093897361247308647679;
k=1;
C1 = 0.104222070006098753381362e2 * C4;
C2 = -0.432673264599913402038813e1 * C4;
C3 = -0.153360524625649062085361e0 * C4;
omegaRoot = betaRoot * sqrt(mu / rho);
Pi=3.14159265358979323846264;
                Real ur,vr,pr; 
                int i1,i2,i3;
                FOR_3D(i1,i2,i3,I1,I2,I3)
                {
          // Reference coordinates:
                    Real x = vertex(i1,i2,i3,0);
                    Real y = vertex(i1,i2,i3,1);
                    Real r = sqrt( x*x + y*y );
                    Real theta = atan2(y,x); 
                    Real cosTheta = x/r;
                    Real sinTheta = y/r;
                        ur = 0.1e1 / r * jn(1, 0.130093897361247308647679e2 * r) * C4 * cos(t * omegaRoot) * cos(theta) + 0.1e1 / r * jn(1, 0.130093897361247308647679e2 * r) * C4 * sin(t * omegaRoot) * sin(theta) + 0.497519350531880418705242e4 / r * yn(1, 0.130093897361247308647679e2 * r) * C4 / (-0.259171103056823175687526e5 * mu + 0.280781480658330616286753e4) * mu * cos(t * omegaRoot) * cos(theta) + 0.497519350531880418705242e4 / r * yn(1, 0.130093897361247308647679e2 * r) * C4 / (-0.259171103056823175687526e5 * mu + 0.280781480658330616286753e4) * mu * sin(t * omegaRoot) * sin(theta) - 0.143113982384782921549911e4 / r * yn(1, 0.130093897361247308647679e2 * r) * C4 / (-0.259171103056823175687526e5 * mu + 0.280781480658330616286753e4) * cos(t * omegaRoot) * cos(theta) - 0.143113982384782921549911e4 / r * yn(1, 0.130093897361247308647679e2 * r) * C4 / (-0.259171103056823175687526e5 * mu + 0.280781480658330616286753e4) * sin(t * omegaRoot) * sin(theta) + 0.520639974058253950243800e3 * C4 / (-0.259171103056823175687526e5 * mu + 0.280781480658330616286753e4) * cos(t * omegaRoot) * cos(theta) + 0.520639974058253950243800e3 * C4 / (-0.259171103056823175687526e5 * mu + 0.280781480658330616286753e4) * sin(t * omegaRoot) * sin(theta) + 0.701497290712961400003724e2 / mu * C4 / (-0.259171103056823175687526e5 * mu + 0.280781480658330616286753e4) * cos(t * omegaRoot) * cos(theta) + 0.701497290712961400003724e2 / mu * C4 / (-0.259171103056823175687526e5 * mu + 0.280781480658330616286753e4) * sin(t * omegaRoot) * sin(theta) + 0.142309060522577509283130e4 * pow(r, -0.2e1) * C4 / (-0.259171103056823175687526e5 * mu + 0.280781480658330616286753e4) * cos(t * omegaRoot) * cos(theta) + 0.142309060522577509283130e4 * pow(r, -0.2e1) * C4 / (-0.259171103056823175687526e5 * mu + 0.280781480658330616286753e4) * sin(t * omegaRoot) * sin(theta);
                        vr = 0.143113982384782921549911e4 / r * yn(1, 0.130093897361247308647679e2 * r) * C4 / (-0.259171103056823175687526e5 * mu + 0.280781480658330616286753e4) * sin(t * omegaRoot) * cos(theta) - 0.143113982384782921549911e4 / r * yn(1, 0.130093897361247308647679e2 * r) * C4 / (-0.259171103056823175687526e5 * mu + 0.280781480658330616286753e4) * cos(t * omegaRoot) * sin(theta) + 0.130093897361247308647679e2 * jn(0, 0.130093897361247308647679e2 * r) * C4 * sin(t * omegaRoot) * cos(theta) - 0.130093897361247308647679e2 * jn(0, 0.130093897361247308647679e2 * r) * C4 * cos(t * omegaRoot) * sin(theta) + 0.647242313233288727871400e5 * mu * yn(0, 0.130093897361247308647679e2 * r) * C4 / (-0.259171103056823175687526e5 * mu + 0.280781480658330616286753e4) * sin(t * omegaRoot) * cos(theta) - 0.647242313233288727871400e5 * mu * yn(0, 0.130093897361247308647679e2 * r) * C4 / (-0.259171103056823175687526e5 * mu + 0.280781480658330616286753e4) * cos(t * omegaRoot) * sin(theta) - 0.186182557353253047298274e5 * yn(0, 0.130093897361247308647679e2 * r) * C4 / (-0.259171103056823175687526e5 * mu + 0.280781480658330616286753e4) * sin(t * omegaRoot) * cos(theta) + 0.186182557353253047298274e5 * yn(0, 0.130093897361247308647679e2 * r) * C4 / (-0.259171103056823175687526e5 * mu + 0.280781480658330616286753e4) * cos(t * omegaRoot) * sin(theta) - 0.999999999999999999999997e0 / r * jn(1, 0.130093897361247308647679e2 * r) * C4 * sin(t * omegaRoot) * cos(theta) + 0.999999999999999999999997e0 / r * jn(1, 0.130093897361247308647679e2 * r) * C4 * cos(t * omegaRoot) * sin(theta) + 0.701497290712961400003724e2 / mu * C4 / (-0.259171103056823175687526e5 * mu + 0.280781480658330616286753e4) * sin(t * omegaRoot) * cos(theta) - 0.701497290712961400003724e2 / mu * C4 / (-0.259171103056823175687526e5 * mu + 0.280781480658330616286753e4) * cos(t * omegaRoot) * sin(theta) - 0.142309060522577509283130e4 * pow(r, -0.2e1) * C4 / (-0.259171103056823175687526e5 * mu + 0.280781480658330616286753e4) * sin(t * omegaRoot) * cos(theta) + 0.142309060522577509283130e4 * pow(r, -0.2e1) * C4 / (-0.259171103056823175687526e5 * mu + 0.280781480658330616286753e4) * cos(t * omegaRoot) * sin(theta) + 0.520639974058253950243800e3 * C4 / (-0.259171103056823175687526e5 * mu + 0.280781480658330616286753e4) * sin(t * omegaRoot) * cos(theta) - 0.520639974058253950243800e3 * C4 / (-0.259171103056823175687526e5 * mu + 0.280781480658330616286753e4) * cos(t * omegaRoot) * sin(theta) - 0.497519350531880418705240e4 / r * yn(1, 0.130093897361247308647679e2 * r) * C4 / (-0.259171103056823175687526e5 * mu + 0.280781480658330616286753e4) * mu * sin(t * omegaRoot) * cos(theta) + 0.497519350531880418705240e4 / r * yn(1, 0.130093897361247308647679e2 * r) * C4 / (-0.259171103056823175687526e5 * mu + 0.280781480658330616286753e4) * mu * cos(t * omegaRoot) * sin(theta);
                        pr = -0.240849861329871999774099e6 * C4 * mu / (-0.259171103056823175687526e5 * mu + 0.280781480658330616286753e4) / pow(r, k) * cos(k * theta) * cos(t * omegaRoot) - 0.240849861329871999774099e6 * C4 * mu / (-0.259171103056823175687526e5 * mu + 0.280781480658330616286753e4) / pow(r, k) * sin(k * theta) * sin(t * omegaRoot) + 0.881153069904669756485534e5 * C4 / (-0.259171103056823175687526e5 * mu + 0.280781480658330616286753e4) * pow(r, k) * mu * cos(k * theta) * cos(t * omegaRoot) + 0.881153069904669756485534e5 * C4 / (-0.259171103056823175687526e5 * mu + 0.280781480658330616286753e4) * pow(r, k) * mu * sin(k * theta) * sin(t * omegaRoot) + 0.118724362715255684088706e5 * C4 / (-0.259171103056823175687526e5 * mu + 0.280781480658330616286753e4) * pow(r, k) * cos(k * theta) * cos(t * omegaRoot) + 0.118724362715255684088706e5 * C4 / (-0.259171103056823175687526e5 * mu + 0.280781480658330616286753e4) * pow(r, k) * sin(k * theta) * sin(t * omegaRoot);
          // (u1,u2) = ur*rHat + vr*thetaHat 
                    u(i1,i2,i3,u1c) = ur*cosTheta - vr*sinTheta;
                    u(i1,i2,i3,u2c) = ur*sinTheta + vr*cosTheta;
                    u(i1,i2,i3,pc ) = pr;
                }
        }

    }
      else if( userKnownSolution=="incompressibleDiskSolution" )
    {
    // --- incompressible disk solution -----
        const int pc = dbase.get<int >("pc");
        assert( pc>=0 );

        MappedGrid & mg = cg[grid];
        mg.update(MappedGrid::THEvertex | MappedGrid::THEcenter);
        OV_GET_SERIAL_ARRAY_CONST(real,mg.vertex(),vertex);
        RealArray & u = ua;

        int & iswCase = db.get<int>("iswCase");

        if( t<3*dt )
            printF("UDKS: evaluate the incompressible disk solution: iswCase=%i, t=%9.3e\n",iswCase,t);


        Real C4 =1.; // defines the scale

        if( iswCase==1 )
        {
      // files generated from AMP/ism/maple/DiskAnnulusSolutionWDH.mw
      // C4 = 1./50.;  // pressure was size 100 with C4=1
// ------------------------------------------------------------------------------------------------ 
//    DiskDispBC: Exact solution for incompressible linear elasticity 
// Solid disk of radius R=1; r=0 bounded solution; r=1 Displacement BC. 
// File Written by Disk_Annulus_ExactSoln.mw (maple) 
// Set C4 to scale the solution 
// ------------------------------------------------------------------------------------------------ 
Real rho,c,mu,omegaRoot,betaRoot,Pi,k,C1,C2,C3,ra,rb; 

rho=1;
mu=1;
betaRoot=6.38016189592398350623661;
k=2;
C1 = 0.0e0;
C2 = 0.607082650035880738991908e1 * C4;
omegaRoot = betaRoot * sqrt(mu / rho);
Pi=3.14159265358979323846264;
                Real ur,vr,pr; 
                int i1,i2,i3;
                Real theta; 
                const Real eps=1.e-9;   // ********************************** FIX ME *********
                if( t<= dt )
                    printF("\n *** WARNING *** need a better fix for disk solution at r=0 ******\n\n");
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
                        ur = 0.1e1 / r * jn(k, 0.638016189592398350623661e1 * r) * C4 * cos(k * theta) * cos(t * omegaRoot) + 0.1e1 / r * jn(k, 0.638016189592398350623661e1 * r) * C4 * sin(k * theta) * sin(t * omegaRoot) - pow(r, k) / r * jn(k, 0.638016189592398350623661e1) * C4 * cos(k * theta) * cos(t * omegaRoot) - pow(r, k) / r * jn(k, 0.638016189592398350623661e1) * C4 * sin(k * theta) * sin(t * omegaRoot);
                        vr = 0.638016189592398350623661e1 / (double) k * jn(k + 1, 0.638016189592398350623661e1 * r) * C4 * sin((double) (k * theta)) * cos(t * omegaRoot) - 0.638016189592398350623661e1 / (double) k * jn(k + 1, 0.638016189592398350623661e1 * r) * C4 * cos((double) (k * theta)) * sin(t * omegaRoot) - 0.100000000000000000000000e1 / r * jn(k, 0.638016189592398350623661e1 * r) * C4 * sin((double) (k * theta)) * cos(t * omegaRoot) + 0.100000000000000000000000e1 / r * jn(k, 0.638016189592398350623661e1 * r) * C4 * cos((double) (k * theta)) * sin(t * omegaRoot) + 0.100000000000000000000000e1 * pow(r, (double) k) / r * jn(k, 0.638016189592398350623661e1) * C4 * sin((double) (k * theta)) * cos(t * omegaRoot) - 0.100000000000000000000000e1 * pow(r, (double) k) / r * jn(k, 0.638016189592398350623661e1) * C4 * cos((double) (k * theta)) * sin(t * omegaRoot);
                        pr = -0.407064658182003197420524e2 * jn(k, 0.638016189592398350623661e1) * C4 / k * mu * pow(r, k) * cos(k * theta) * cos(t * omegaRoot) - 0.407064658182003197420524e2 * jn(k, 0.638016189592398350623661e1) * C4 / k * mu * pow(r, k) * sin(k * theta) * sin(t * omegaRoot);
          // (u1,u2) = ur*rHat + vr*thetaHat 
                    u(i1,i2,i3,u1c) = ur*cosTheta - vr*sinTheta;
                    u(i1,i2,i3,u2c) = ur*sinTheta + vr*cosTheta;
                    u(i1,i2,i3,pc ) = pr;
                }
        }
        else if( iswCase==2 )
        {
// ------------------------------------------------------------------------------------------------ 
//    DiskTracBC: Exact solution for incompressible linear elasticity 
// Solid disk of radius R=1; r=0 bounded solution; r=1 Traction BC. 
// File Written by Disk_Annulus_ExactSoln.mw (maple) 
// Set C4 to scale the solution 
// ------------------------------------------------------------------------------------------------ 
Real rho,c,mu,omegaRoot,betaRoot,Pi,k,C1,C2,C3,ra,rb; 
ra=0.;
rb=1;
rho=1;
mu=1;
betaRoot=8.17667024339775129040010;
k=2;
C1 = 0.0e0;
C2 = 0.486843886152537642280292e1 * C4;
omegaRoot = betaRoot * sqrt(mu / rho);
Pi=3.14159265358979323846264;
                Real ur,vr,pr; 
                int i1,i2,i3;
                Real theta; 
                const Real eps=1.e-9;   // ********************************** FIX ME *********
                if( t<= dt )
                    printF("\n *** WARNING *** need a better fix for disk solution at r=0 ******\n\n");
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
                        ur = 0.1e1 / r * jn(k, 0.817667024339775129040010e1 * r) * C4 * cos(k * theta) * cos(t * omegaRoot) + 0.1e1 / r * jn(k, 0.817667024339775129040010e1 * r) * C4 * sin(k * theta) * sin(t * omegaRoot) - 0.163533404867955025808002e2 * pow(r, k) / r * C4 / (-0.2e1 * k * k + 0.2e1 * k + 0.668579362692662413298169e2 * mu) * k * jn(k + 0.1e1, 0.817667024339775129040010e1) * cos(k * theta) * cos(t * omegaRoot) - 0.163533404867955025808002e2 * pow(r, k) / r * C4 / (-0.2e1 * k * k + 0.2e1 * k + 0.668579362692662413298169e2 * mu) * k * jn(k + 0.1e1, 0.817667024339775129040010e1) * sin(k * theta) * sin(t * omegaRoot) + 0.2e1 * pow(r, k) / r * C4 / (-0.2e1 * k * k + 0.2e1 * k + 0.668579362692662413298169e2 * mu) * k * k * jn(k, 0.817667024339775129040010e1) * cos(k * theta) * cos(t * omegaRoot) + 0.2e1 * pow(r, k) / r * C4 / (-0.2e1 * k * k + 0.2e1 * k + 0.668579362692662413298169e2 * mu) * k * k * jn(k, 0.817667024339775129040010e1) * sin(k * theta) * sin(t * omegaRoot) - 0.2e1 * pow(r, k) / r * C4 / (-0.2e1 * k * k + 0.2e1 * k + 0.668579362692662413298169e2 * mu) * k * jn(k, 0.817667024339775129040010e1) * cos(k * theta) * cos(t * omegaRoot) - 0.2e1 * pow(r, k) / r * C4 / (-0.2e1 * k * k + 0.2e1 * k + 0.668579362692662413298169e2 * mu) * k * jn(k, 0.817667024339775129040010e1) * sin(k * theta) * sin(t * omegaRoot);
                        vr = 0.200000000000000000000000e1 * pow(r, k) / r * C4 / (-0.2e1 * k * k + 0.2e1 * k + 0.668579362692662413298169e2 * mu) * k * k * jn(k, 0.817667024339775129040010e1) * cos(k * theta) * sin(t * omegaRoot) - 0.817667024339775129040010e1 / k * jn(k + 0.1e1, 0.817667024339775129040010e1 * r) * C4 * cos(k * theta) * sin(t * omegaRoot) - 0.100000000000000000000000e1 / r * jn(k, 0.817667024339775129040010e1 * r) * C4 * sin(k * theta) * cos(t * omegaRoot) - 0.200000000000000000000000e1 * pow(r, k) / r * C4 / (-0.2e1 * k * k + 0.2e1 * k + 0.668579362692662413298169e2 * mu) * k * k * jn(k, 0.817667024339775129040010e1) * sin(k * theta) * cos(t * omegaRoot) + 0.100000000000000000000000e1 / r * jn(k, 0.817667024339775129040010e1 * r) * C4 * cos(k * theta) * sin(t * omegaRoot) + 0.817667024339775129040010e1 / k * jn(k + 0.1e1, 0.817667024339775129040010e1 * r) * C4 * sin(k * theta) * cos(t * omegaRoot) - 0.200000000000000000000000e1 * pow(r, k) / r * C4 / (-0.2e1 * k * k + 0.2e1 * k + 0.668579362692662413298169e2 * mu) * k * jn(k, 0.817667024339775129040010e1) * cos(k * theta) * sin(t * omegaRoot) - 0.163533404867955025808002e2 * pow(r, k) / r * C4 / (-0.2e1 * k * k + 0.2e1 * k + 0.668579362692662413298169e2 * mu) * k * jn(k + 0.1e1, 0.817667024339775129040010e1) * cos(k * theta) * sin(t * omegaRoot) + 0.200000000000000000000000e1 * pow(r, k) / r * C4 / (-0.2e1 * k * k + 0.2e1 * k + 0.668579362692662413298169e2 * mu) * k * jn(k, 0.817667024339775129040010e1) * sin(k * theta) * cos(t * omegaRoot) + 0.163533404867955025808002e2 * pow(r, k) / r * C4 / (-0.2e1 * k * k + 0.2e1 * k + 0.668579362692662413298169e2 * mu) * k * jn(k + 0.1e1, 0.817667024339775129040010e1) * sin(k * theta) * cos(t * omegaRoot);
                        pr = -0.109335059605578508270112e4 * C4 / (-(double) (2 * k * k) + (double) (2 * k) + 0.668579362692662413298169e2 * mu) * mu * pow(r, (double) k) * jn(k + 1, 0.817667024339775129040010e1) * cos((double) (k * theta)) * cos(t * omegaRoot) - 0.109335059605578508270112e4 * C4 / (-(double) (2 * k * k) + (double) (2 * k) + 0.668579362692662413298169e2 * mu) * mu * pow(r, (double) k) * jn(k + 1, 0.817667024339775129040010e1) * sin((double) (k * theta)) * sin(t * omegaRoot) + 0.133715872538532482659634e3 * C4 / (-(double) (2 * k * k) + (double) (2 * k) + 0.668579362692662413298169e2 * mu) * mu * pow(r, (double) k) * (double) k * jn(k, 0.817667024339775129040010e1) * cos((double) (k * theta)) * cos(t * omegaRoot) + 0.133715872538532482659634e3 * C4 / (-(double) (2 * k * k) + (double) (2 * k) + 0.668579362692662413298169e2 * mu) * mu * pow(r, (double) k) * (double) k * jn(k, 0.817667024339775129040010e1) * sin((double) (k * theta)) * sin(t * omegaRoot) - 0.133715872538532482659634e3 * C4 / (-(double) (2 * k * k) + (double) (2 * k) + 0.668579362692662413298169e2 * mu) * mu * pow(r, (double) k) * jn(k, 0.817667024339775129040010e1) * cos((double) (k * theta)) * cos(t * omegaRoot) - 0.133715872538532482659634e3 * C4 / (-(double) (2 * k * k) + (double) (2 * k) + 0.668579362692662413298169e2 * mu) * mu * pow(r, (double) k) * jn(k, 0.817667024339775129040010e1) * sin((double) (k * theta)) * sin(t * omegaRoot);
          // (u1,u2) = ur*rHat + vr*thetaHat 
                    u(i1,i2,i3,u1c) = ur*cosTheta - vr*sinTheta;
                    u(i1,i2,i3,u2c) = ur*sinTheta + vr*cosTheta;
                    u(i1,i2,i3,pc ) = pr;
                }
        }
        else
        {
            OV_ABORT("UDKS: unknown icase=");
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
            "incompressible surface wave",
            "incompressible annulus solution",
            "incompressible disk solution",
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
                          "     Case 1: traction-dispalcement BC: surface wave as speed cs=omega/k, region=[0,1]x[0,Ly] Ly is given.\n"
                          "     Case 2: traction-traction BC\n"
          //    "     Case 2: surface wave as speed different from cs, region=[0,1]x[0,1].\n"
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
                          " The solution is defined on an Annulus\n"
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
            gi.inputString(answer,"Enter iswCase");
            sScanF(answer,"%i",&iswCase);

        }

        else if( answer=="incompressible disk solution" )
        {
            userKnownSolution="incompressibleDiskSolution";
            dbase.get<bool>("knownSolutionIsTimeDependent")=true;  // known solution IS time dependent

            printF("--- Incompressible Disk Solution ---\n"
                          " This is an exact solution to the incompressible linear elastic equations.\n"
                          " The solution is defined on the unit disk\n"
                          " The are different cases defined:\n"
                          "     Case 1: displacement BC\n"
                          "     Case 2: traction BC\n"
                          );      

            if( !db.has_key("iswCase") )
            { // Create parameters for the incompressible annulus solution
                db.put<int>("iswCase")=1;
            }
            int & iswCase = db.get<int>("iswCase");
            gi.inputString(answer,"Enter iswCase");
            sScanF(answer,"%i",&iswCase);

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
