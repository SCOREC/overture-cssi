#include "AdParameters.h"
#include "GenericGraphicsInterface.h"
#include "ParallelUtility.h"
#include "Interface.h"  

#include "AdExactSolutions.h"

// WARNING: Overture type "real" conflicts with complex "real" 
// Use Overture type "Real" Instead
#include <complex>

// typedef ::real Real;
// typedef ::real Real;    

#define ForBoundary(side,axis)   for( int axis=0; axis<mg.numberOfDimensions(); axis++ ) \
                                 for( int side=0; side<=1; side++ )

// Forward declaration: 
Parameters & getInterfaceParameters( int grid, int side, int axis, Parameters & parameters);                                  

// Bessel functions taking complex arguments: 
#define cBesselJ EXTERN_C_NAME(cbesselj)
#define cBesselY EXTERN_C_NAME(cbessely)

extern "C"
{

  void cBesselJ( const Real & nu, const Real &zr, const Real &zi, Real &jr,Real &ji );
  void cBesselY( const Real & nu, const Real &zr, const Real &zi, Real &yr,Real &yi );

}

// namespace
// {
// // These are for passing to the rotating disk function
// Real innerRadius, outerRadius, inertiaTerm, massTerm;
// }



#define FOR_3D(i1,i2,i3,I1,I2,I3)                                       \
int I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase(); \
int I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(); \
for(int i3=I3Base; i3<=I3Bound; i3++)                                       \
  for(int i2=I2Base; i2<=I2Bound; i2++)                                     \
    for(int i1=I1Base; i1<=I1Bound; i1++)


// Macro to get the vertex array
#define GET_VERTEX_ARRAY(x)                                     \
mg.update(MappedGrid::THEvertex | MappedGrid::THEcenter);       \
OV_GET_SERIAL_ARRAY_CONST(real,mg.vertex(),x);                  \
if( !thisProcessorHasPoints )                                   \
  return 0; // no points on this processor


// // =========================================================================================
// /// \brief Class to evaluate the exact solution pressure for the time-dependent eigenfunction for
// ///    circular Couette flow. 
// //
// //     p(r)  = int_a^r v(s)^2/s  ds
// //     v(r) = BesselJ(1,lambda*r)*Y1a - J1a*BesselY(1,lambda*r)
// // 
// // See file Dropbox/CG/fins/codes/besselZeros.maple
// // =========================================================================================
// class CircularCouettePressure
// {
//   public:

//   CircularCouettePressure(){ pData=NULL;  a=.5; b=1.; ra=a-.125; rb=b+.125; numData=0; };
//   ~CircularCouettePressure(){ delete [] pData; };

//   int initialize( Real amp1, int nFreq1, Real rInner, Real rOuter, Real omegaInner, Real omegaOuter, Real rho_, Real nu, Real t );
//   Real evalPressure( Real r );
//   Real evalF(  Real r );
//   Real getLambda();
//   Real getIntegral( Real x0, Real x1 );
  
//   int numData;
//   Real a,b;    // inner and outer radii
//   Real ra,rb;  // extended region overwhich data is computed (to allow for ghost)
//   Real *pData;

//   Real rho,lambda,A,B,amp,j1a,y1a;
// };

// Real CircularCouettePressure::
// getLambda()
// {
//   return lambda;
// }


// // =========================================================================================
// /// \brief Evaluate the integrand for the pressure 
// // =========================================================================================
// Real CircularCouettePressure::
// evalF( Real r )
// {
  
//   Real lamr = lambda*r;
//   Real uTheta;
//   if( a!=0. )
//    uTheta=A/r + B*r + amp* ( jn(1,lamr)*y1a - yn(1,lamr)*j1a );
//   else
//     uTheta= B*r + amp* ( jn(1,lamr) );

//   Real f;
//   if( fabs(r) > REAL_MIN*10. )
//   {
//    f = rho*uTheta*uTheta/r;
//   }
//   else
//   {
//     // Note that J1(r) = r/2  as r -> 0 
//     f = rho*uTheta*( B + amp*lamr/2. );
//   }
  

//   return f;
  
// }

// // =========================================================================================
// /// \brief Evaluate the integral from x0 to x1 
// // =========================================================================================
// Real CircularCouettePressure::
// getIntegral( Real x0, Real x1 )
// {

//   Real xm=.5*(x0+x1);
//   Real f0 = evalF(x0);
//   Real fm = evalF(xm);
//   Real f1 = evalF(x1);
//   Real h = .5*(x1-x0);
//   Real fint = (f0 +4.*fm + f1)*h/3.; // Simpson's rule 
// //  Real fint = (f0 + f1)*h; // Trap
//   return fint;
// }





int AdParameters::
getUserDefinedKnownSolution(Real t, CompositeGrid & cg, int grid, RealArray & ua, 
                            const Index & I1_, const Index &I2_, const Index &I3_, 
                            int numberOfTimeDerivatives /* = 0 */ )
// ==========================================================================================
///  \brief Evaluate a user defined known solution.
/// \param numberOfTimeDerivatives (input) : number of time derivatives to evaluate (0 means evaluate the solution).
// ==========================================================================================
{
  MappedGrid & mg = cg[grid];
  const int numberOfDimensions = mg.numberOfDimensions();
    
  // Adjust index bounds for parallel *wdh* 2017/05/31 
  Index I1=I1_, I2=I2_, I3=I3_;
  
  mg.update(MappedGrid::THEmask); // *wdh* Dec 17, 2018 -- needed in paralle 
  
  OV_GET_SERIAL_ARRAY_CONST(int,mg.mask(),mask);
  int includeGhost=1;
  bool thisProcessorHasPoints=ParallelUtility::getLocalArrayBounds(mg.mask(),mask,I1,I2,I3,includeGhost);
  
  // printF("--AdParameters-- getUserDefinedKnownSolution: START---\n");

  if( ! dbase.get<DataBase >("modelData").has_key("userDefinedKnownSolutionData") )
  {
    printf("getUserDefinedKnownSolution:ERROR: sub-directory `userDefinedKnownSolutionData' not found!\n");
    OV_ABORT("error");
  }
  DataBase & db =  dbase.get<DataBase >("modelData").get<DataBase>("userDefinedKnownSolutionData");

  const aString & userKnownSolution = db.get<aString>("userKnownSolution");

  Real *rpar = db.get<Real[20]>("rpar");
  int *ipar = db.get<int[20]>("ipar");
  
  const int numberOfComponents =dbase.get<int >("numberOfComponents") - dbase.get<int >("numberOfExtraVariables");
  std::vector<Real> & kappa = dbase.get<std::vector<Real> >("kappa");
  const Real & thermalConductivity = dbase.get<Real>("thermalConductivity");

  const Real dt = dbase.get<Real>("dt");

  if( numberOfTimeDerivatives !=0 )
  {
    OV_ABORT("CgAd:userDefinedKnownSolution: numberOfTimeDerivatives !=0");
  }  
  
  if( userKnownSolution=="squareEigenfunction" )
  {
    // --- Eigenfunction for the heat equation in a square ---
 
    const Real & kx = rpar[0];
    const Real & ky = rpar[1];
    const Real & kz = rpar[2];

    const int bcOpt = ipar[0];

    if( t<=3.*dt )
      printF("heat equation eigenfunction: kx=%g, ky=%g, kz=%g, bcOpt=%d, t=%9.3e\n",kx,ky,kz,bcOpt,t);

    // --- we could avoid building the vertex array on Cartesian grids ---
    GET_VERTEX_ARRAY(x);

    if( mg.numberOfDimensions()==2 )
    {
      for( int n=0; n<numberOfComponents; n++ )
      {
        const Real b = kappa[n]*SQR(Pi)*( SQR(kx) + SQR(ky) );
        if( bcOpt==0 )
        { // Dirichlet BCs
          ua(I1,I2,I3,n) = exp( -b*t )*sin( (Pi*kx)*x(I1,I2,I3,0) )*sin( (Pi*ky)*x(I1,I2,I3,1) );
        }
        else
        { // Neumann BCs
          ua(I1,I2,I3,n) = exp( -b*t )*cos( (Pi*kx)*x(I1,I2,I3,0) )*cos( (Pi*ky)*x(I1,I2,I3,1) );
        }
      }
    }
    else
    { // ----- 3D ----
      for( int n=0; n<numberOfComponents; n++ )
      {
        const Real b = kappa[n]*SQR(Pi)*( SQR(kx) + SQR(ky) + SQR(kz) );
        if( bcOpt==0 )
        { // Dirichlet BCs
          ua(I1,I2,I3,n) = exp( -b*t )*sin( (Pi*kx)*x(I1,I2,I3,0) ) *sin( (Pi*ky)*x(I1,I2,I3,1) ) *sin( (Pi*kz)*x(I1,I2,I3,2) );
        }
        else
        { // Neumann BCs
          ua(I1,I2,I3,n) = exp( -b*t )*cos( (Pi*kx)*x(I1,I2,I3,0) ) *cos( (Pi*ky)*x(I1,I2,I3,1) ) *cos( (Pi*kz)*x(I1,I2,I3,2) );
        }        
      }
    }


  }

  else if( userKnownSolution=="diskEigenfunction" )
  {
    // --- Eigenfunction for the heat equation in a disk ---

    // --- we could avoid building the vertex array on Cartesian grids ---
    GET_VERTEX_ARRAY(x);

   
 
    const int n     = ipar[0];  // angular number, n=0,1,... --> Jn(omega*r)
    const int m     = ipar[1];  // radial number m=0,... 
    const int bcOpt = ipar[2];

    const Real a   = rpar[0];  // radius of disk
    const Real amp = rpar[1];  // amplitude 

    Real lambda; 
    if( bcOpt==0 )
    {
      // DIRICHLET BC

      #include "../../mx/src/besselZeros.h"    
      assert( m<mdbz && n<ndbz );
   
      const Real jzmn = besselZeros[n][m];  // m'th zero of Jn

      lambda=jzmn/a;
    }
    else
    {
      // Neumann BC
      #include "../../mx/src/besselPrimeZeros.h"    
      assert( m<mdbpz && n<ndbpz );
   
      const Real jzmn = besselPrimeZeros[n][m];  // m'th zero of Jn' (excluding r=0 for J0)

      lambda=jzmn/a;

    }
     

    if( t<=3.*dt )
      printF("Disk: Bessel function solution: a=%g, n=%i, m=%i, lambda=%e, bcOpt=%d \n",a,n,m,lambda,bcOpt);

    for( int nn=0; nn<numberOfComponents; nn++ )
    {
      const Real expt = amp*exp(-kappa[nn]*lambda*lambda*t); 

      FOR_3D(i1,i2,i3,I1,I2,I3)
      { 
        const Real xd = x(i1,i2,i3,0), yd = x(i1,i2,i3,1);
        const Real theta = atan2(yd,xd);
        const Real r = sqrt( xd*xd + yd*yd );

        ua(i1,i2,i3,nn) = jn(n,lambda*r)*cos(n*theta)*expt;
      }
    }

  }

  else if( userKnownSolution=="annulusEigenfunction" )
  {
    // --- Eigenfunction for the heat equation in an annulus ---

    // We can NOT check the boundaryCondition array since it may not be set yet !
   //  bool dirichlet = (mg.boundaryCondition(0,1)==AdParameters::dirichletBoundaryCondition &&
   //                    mg.boundaryCondition(1,1)==AdParameters::dirichletBoundaryCondition );
   //  bool neumann   = (mg.boundaryCondition(0,1)==AdParameters::neumannBoundaryCondition &&
   //                    mg.boundaryCondition(1,1)==AdParameters::neumannBoundaryCondition );

   // dirichlet=true;

   // assert( dirichlet || neumann );

    // --- we could avoid building the vertex array on Cartesian grids ---
    GET_VERTEX_ARRAY(x);

    const int n     = ipar[0];  // angular number, n=0,1,... --> Jn(lambda*r), Yn(lambda*r)
    const int m     = ipar[1];  // radial number m=0,... 
    const int bcOpt = ipar[2];

    const Real amp  = rpar[0];  // amplitude 

    Real lambda, cJ,cY;
    if( bcOpt==0 )
    {
      // Dirichlet boundary conditions:
      //     det(lambda) =  Jn(lambda*a)*Yn(lambda*b) - Jn(lambda*b)*Yn(lambda*a) = 0 
      // This solution assumes ra=.5 and rb=1 (set in the next file)
      #include "../codes/annulusEigenvaluesHeatEquationDirichlet.h" 
      lambda = annulusEigs[n][m];  // m'th zero of the determinant condition)  
      assert( m<numRoot && n<numBesselOrder );
      // Here is the eignevctor [cJ,cY] 
      cJ = yn(n,lambda*ra);
      cY =-jn(n,lambda*ra);
    } 
    else 
    {
      // Neumann boundary conditions:
      //     det(lambda) =  Jn'(lambda*a)*Yn'(lambda*b) - Jn'(lambda*b)*Yn'(lambda*a) = 0 
      // This solution assumes ra=.5 and rb=1 (set in the next file)
      #include "../codes/annulusEigenvaluesHeatEquationNeumann.h"   
      lambda = annulusEigs[n][m];  // m'th zero of the detreminant condition)
      assert( m<numRoot && n<numBesselOrder );

      // Here is the eignevctor [cJ,cY] 
      // Use Jn'(z) = (n/z)Jn - J_{n+1}
      const Real z = lambda*ra;
      cJ = ( (n/z)*yn(n,z) - yn(n+1,z) );    //  Yn'(lambda*a)
      cY =-( (n/z)*jn(n,z) - jn(n+1,z) );    // -Jn'(lambda*a)
    }

    // scale eigenvector so solution is of size amp
    Real cNorm = sqrt( cJ*cJ + cY*cY );
    cJ *= amp/cNorm; 
    cY *= amp/cNorm;

    if( t<=3.*dt )
      printF("Annulus: Bessel function solution: amp=%g, n=%i, m=%i, lambda=%e, bcOpt=%d\n", amp,n,m,lambda,bcOpt);

    for( int nn=0; nn<numberOfComponents; nn++ )
    {
      const Real expt = amp*exp(-kappa[nn]*lambda*lambda*t); 

      FOR_3D(i1,i2,i3,I1,I2,I3)
      { 
        const Real xd = x(i1,i2,i3,0), yd = x(i1,i2,i3,1);
        const Real theta = atan2(yd,xd);
        const Real r = sqrt( xd*xd + yd*yd );

        ua(i1,i2,i3,nn) = ( cJ*jn(n,lambda*r) + cY*yn(n,lambda*r) )*cos(n*theta)*expt;
      }
    }



  }

  else if( userKnownSolution=="concentricCylinders" )
  {
    // ---------- CONCENTRIC CYLINDERS (ANNULII) EIGENFUNCTIONS : TWO DOMAIN  CHT ----------------
    // See Cgmp reference guide for derivation of the solution

    // --- we could avoid building the vertex array on Cartesian grids ---
    GET_VERTEX_ARRAY(x);

    const int n      = ipar[0];
    const int m      = ipar[1];
    const int option = ipar[2];

    const Real amp  = rpar[0]; 
    const Real a    = rpar[1];
    const Real b    = rpar[2]; 
    const Real c    = rpar[3]; 
    const Real D1   = rpar[4];
    const Real K1   = rpar[5]; 
    const Real D2   = rpar[6]; 
    const Real K2   = rpar[7]; 
     
    // Solution is u = amp * [ c1J*Jn(lam*r) + c1Y*Yn(lam*r) ] cos(m*theta) exp(-Dl*lam^2 t),   a < r < b\n");
    //             u = amp * [ c2J*Jn(lam*r) + c2Y*Yn(lam*r) ] cos(m*theta) exp(-Dr*lam^2 t),   b < r < c\n"); 


    Real s;  // eigenvalue 
    Real D1a = D1, D2a=D2; // for checking value in the include files below

    const Real eps = REAL_EPSILON*100.;
    if( fabs(D1-.1)<eps && fabs(D2-.05)<eps && fabs(K2/K1-2.)<eps )
    {
      // --- Case 1 ----

      // This solution assumes radii : ra=.5, rb=1, rc=1.5 (set in the next file)
      // const Real ra=5.00000000000000e-01, rb=1.00000000000000e+00, rc=1.50000000000000e+00, Kr=2.00000000000000e+00; 
      // const Real D1=1.00000000000000e-01, D2=5.00000000000000e-02; 
      // This next include file was created by ad/codes/concentricAnnulusEigenvalues.maple
      #include "../codes/concentricAnnulusEigenvaluesHeatEquationDirichletCase1.h"  

      if( ra!=a || rb!=b || rc!=c || Kr!=K2/K1 || D1!=D1a || D2!=D2a )
      {
        printF("userDefinedKnownSolution:ERROR: concentricCylinders requires a=%g, b=%g, c=%g, K2/K1=%g, D1=%g, D2=%g\n",ra,rb,rc,Kr,D1,D2);
        printF(" a=%g, b=%g, c=%g, K1=%g, K2=%g => K2/K1=%g D1=%g, D2=%g\n",a,b,c,K1,K2,K2/K1,D1a,D2a);
        OV_ABORT("ERROR");
      }  

      assert( m<numRoot && n<numBesselOrder );
      s = concentricAnnulusEigs[n][m];  // m'th zero of the determinant condition)
    }
    else
    {
      // --- Case 2 ----

      // const Real ra=5.00000000000000e-01, rb=1.00000000000000e+00, rc=1.50000000000000e+00, Kr=5.00000000000000e-01; 
      // const Real D1=9.00000000000000e-01, D2=8.00000000000000e-01; 

      #include "../codes/concentricAnnulusEigenvaluesHeatEquationDirichletCase2.h" 

      if( ra!=a || rb!=b || rc!=c || Kr!=K2/K1 || D1!=D1a || D2!=D2a )
      {
        printF("userDefinedKnownSolution:ERROR: concentricCylinders requires a=%g, b=%g, c=%g, K2/K1=%g, D1=%g, D2=%g\n",ra,rb,rc,Kr,D1,D2);
        printF(" a=%g, b=%g, c=%g, K1=%g, K2=%g => K2/K1=%g D1=%g, D2=%g\n",a,b,c,K1,K2,K2/K1,D1a,D2a);
        OV_ABORT("ERROR");
      } 

      assert( m<numRoot && n<numBesselOrder );
      s = concentricAnnulusEigs[n][m];  // m'th zero of the determinant condition)

    }

    //  c1J = Jn(lam*b) - (Jn(lam*c)/(Yn(lam*c)*Yn(lam*b) 
    //  c2J = Jn(lam*b) - (Jn(lam*a)/(Yn(lam*a)*Yn(lam*b) 
    //  c1Y = -Jn(lamba*a)/Yn(lambda*a) * c1J 
    //  c2Y = -Jn(lamba*c)/Yn(lambda*c) * c2J

    const Real alpha1=1./sqrt(D1), alpha2=1./sqrt(D2);
    const Real lambda1=alpha1*s, lambda2=alpha2*s;

    Real c1J = jn(n,lambda2*b) - (jn(n,lambda2*c)/yn(n,lambda2*c)) *yn(n,lambda2*b);
    Real c2J = jn(n,lambda1*b) - (jn(n,lambda1*a)/yn(n,lambda1*a)) *yn(n,lambda1*b);

    Real c1Y = -(jn(n,lambda1*a)/yn(n,lambda1*a)) * c1J;
    Real c2Y = -(jn(n,lambda2*c)/yn(n,lambda2*c)) * c2J;

    // Normalize solution to be approximately of size "amp" 
    const Real cNorm = sqrt( SQR(c1J) + SQR(c1Y) );
    const Real scale = amp/cNorm;
    c1J *=scale; c2J*=scale; c1Y*=scale; c2Y*=scale;


    if( true && t==0. )
    {
      // CHECK THE SOLUTION : jump conditions
      printF("concentric cylinders: option=%d (0=inner, 1=outer): gridName=[%s]\n",option,(const char*)mg.getName());

      Real r = b; 
      Real z1 = lambda1*r; 
      Real z2 = lambda2*r; 

      Real jn1 = jn(n,z1), yn1 = yn(n,z1);
      Real jn2 = jn(n,z2), yn2 = yn(n,z2);

      Real u1 = ( c1J*jn1 + c1Y*yn1 ); // *cos(m*theta); // no need to include theta, 
      Real u2 = ( c2J*jn2 + c2Y*yn2 ); // *cos(m*theta);

      // Use Jn'(z) = (n/z)Jn - J_{n+1}
      Real jnp1 = lambda1*( (n/z1)*jn1 - jn(n+1,z1) );
      Real ynp1 = lambda1*( (n/z1)*yn1 - yn(n+1,z1) );

      Real jnp2 = lambda2*( (n/z2)*jn2 - jn(n+1,z2) );
      Real ynp2 = lambda2*( (n/z2)*yn2 - yn(n+1,z2) );

      Real u1r = ( c1J*jnp1 + c1Y*ynp1 ); // *cos(m*theta);
      Real u2r = ( c2J*jnp2 + c2Y*ynp2 ); // *cos(m*theta);

      Real jumpFlux = K1*u1r - K2*u2r; 

      if( true )
        printF(" concentricCylinders:ExactSolution:  **CHECK JUMPS** [T]=%8.2e, [K T_r]=%8.2e\n", abs(u1-u2), abs(jumpFlux));

      // --- check that grids have the correct inner and outer radii ---
      int axisI=0;  // radial direction is along this axis
      Index Ib1,Ib2,Ib3;
      for( int sideI=0; sideI<=1; sideI++ )
      {
        getBoundaryIndex(mg.gridIndexRange(),sideI,axisI,Ib1,Ib2,Ib3);
        RealArray rad(Ib1,Ib2,Ib3);
        rad = sqrt( SQR(x(Ib1,Ib2,Ib3,0)) + SQR(x(Ib1,Ib2,Ib3,1)) );

        Real re = option==0 ? ( sideI==0 ? a : b ) : ( sideI==0 ? b : c); // expected radius 
        Real err = max(fabs(rad-re))/re; 
        if( err > REAL_EPSILON*100. )
        {
          Real radEst = max(rad);
          printF("conCyls:ERROR: grid=%d (%s): radius on side %d is not =%g (est. rad=%g)\n",grid,(const char*)mg.getName(),sideI,re,radEst);
          printF(" min(rad)=%9.3e, max(rad)=%9.3e\n",min(rad),max(rad));
          printF(" max(fabs(rad-re))/re = %9.3e\n",err);
          OV_ABORT("error: fix grids, domain order may be wrong.");
        }
        else
        {
          printF("conCyls: grid=%d (%s): radius on side %d is =%g (as expected)\n",grid,(const char*)mg.getName(),sideI,re);

        }
      }

      // OV_ABORT("stop here for now");


    }

    Real D, cJ, cY, lambda;
    if( option==0 )
    {
       // ---- INNER annulus ----
      D = D1; cJ = c1J; cY = c1Y; lambda=lambda1;
      if( kappa[0]!=D1 || thermalConductivity != K1 )
      {
        printF("CGAD:UDKS: error kappa[0]=%g, and thermalConductivity=%g should match D1=%g, and K1=%g\n",kappa[0],thermalConductivity,D1,K1);
        OV_ABORT("ERROR"); 
      }
    }
    else if( option==1 )
    {
       // ---- OUTER annulus ----
      D = D2; cJ = c2J; cY = c2Y; lambda=lambda2; 
      if( kappa[0]!=D2 || thermalConductivity != K2 )
      {
        printF("CGAD:UDKS: error kappa[0]=%g, and thermalConductivity=%g should match D2=%g, and K2=%g\n",kappa[0],thermalConductivity,D2,K2);
        OV_ABORT("ERROR"); 
      }      
    }
    else
    {
      OV_ABORT("ERROR- unexpected value for option"); 
    }

    if( t<=3.*dt )
      printF("CGAD:UDKS: concentric cylinders: n=%d, m=%d, option=%d (0=inner,1=outer), a=%g, b=%g, c=%g,  t=%9.3e\n",n,m,option,a,b,c, t );

  


    // ---- evaluate the solution on the grid ----
    for( int nn=0; nn<numberOfComponents; nn++ )
    {
      const Real expt = amp*exp( - s*s*t ); 

      FOR_3D(i1,i2,i3,I1,I2,I3)
      { 
        const Real xd = x(i1,i2,i3,0), yd = x(i1,i2,i3,1);
        const Real theta = atan2(yd,xd);
        const Real r = sqrt( xd*xd + yd*yd );

        ua(i1,i2,i3,nn) = ( cJ*jn(n,lambda*r) + cY*yn(n,lambda*r) )*cos(n*theta)*expt;

      }
    }


  }

  else if( userKnownSolution=="doubleAnnulus"   ||
           userKnownSolution=="doubleRectangle"    )
  {
    if( userKnownSolution=="doubleAnnulus" )
      printF("Call AdExactSolutions to evaluate the double annulus solution at t=%9.3e\n",t);
    else
      printF("Call AdExactSolutions to evaluate the double rectangle solution at t=%9.3e\n",t);

    AdExactSolutions & adExactSolution = dbase.get<AdExactSolutions>("adExactSolution");
    adExactSolution.evalSolution( t, cg, grid, ua, I1,I2,I3, numberOfTimeDerivatives );

  }
  else if( userKnownSolution=="planeInterfaceMode" )
  {
    // -- Plane interface mode for the Heat equation on two adjacent squares---\n"
    //    u = amp*e^{s*t} * uHat(x) * e^{i*pi*k*y} (form of solution on each side),\n"
    //    uHat(x) = c1*exp( alpha*x) + c2*exp( alpha*x ),\n"
    //    s = sr + I*si (Laplace transform parameter, input)\n"
    //    k = Fourier mode in Y (input)\n"
    //    c1,c2 : differ on two sides of the interface. \n"
    //    option : 0=left-side, 1=right-side.\n"


    const int & option= ipar[0];

    const Real & amp = rpar[0];
    const Real & sr  = rpar[1];
    const Real & si  = rpar[2];
    const Real & ky  = rpar[3];
    const Real & c1  = rpar[4];
    const Real & c2  = rpar[5];
    const Real & Dl  = rpar[6];
    const Real & Kl  = rpar[7];
    const Real & Dr  = rpar[8];
    const Real & Kr  = rpar[9];
    // Advection velocities: 
    const Real & uL  = rpar[10];
    const Real & vL  = rpar[11];
    const Real & wL  = rpar[12];
    const Real & uR  = rpar[13];
    const Real & vR  = rpar[14];
    const Real & wR  = rpar[15];       

    // --- we could avoid building the vertex array on Cartesian grids ---
    GET_VERTEX_ARRAY(x);

    int nn=0;  // component number 
    assert( numberOfComponents==1 );  // do this case for now

    const Real kypi = ky*Pi; 
    std::complex<Real> I(0.0,1.0);             // sqrt(-1)
    std::complex<Real> ss(sr,si);              // complex "s"
    std::complex<Real> alphap, alpham;         // alpha+ and alpha- for the appropriate side
    std::complex<Real> expst = amp*exp(ss*t);


    const int & multiDomainProblem = dbase.get<int>("multiDomainProblem"); 
    std::complex<Real> c1a, c2a;

    // Define velocities scaled by diffusion coeff
    const Real uDL = uL/Dl, vDL = vL/Dl;
    const Real uDR = uR/Dr, vDR = vR/Dr;

    // ****  see formula in CgmpRef/tex/planeInterface.tex ******

    std::complex<Real> alphaLp = .5*uDL + sqrt( SQR(.5*uDL) + ss/Dl + kypi*kypi + I*kypi*vDL ); 
    std::complex<Real> alphaLm = .5*uDL - sqrt( SQR(.5*uDL) + ss/Dl + kypi*kypi + I*kypi*vDL ); 

    std::complex<Real> alphaRp = .5*uDR + sqrt( SQR(.5*uDR) + ss/Dr + kypi*kypi + I*kypi*vDR ); 
    std::complex<Real> alphaRm = .5*uDR - sqrt( SQR(.5*uDR) + ss/Dr + kypi*kypi + I*kypi*vDR ); 

    // std::complex<Real> alphaR = sqrt( ss/Dr + kypi*kypi );     

    if( option==0 )
    {
      // left domain coefficients
      if( Dl != kappa[0] || Kl != thermalConductivity )
      {
        printF(">>>> CGAD:UDKS: Plane interface mode:ERROR: left domain: [Dl,Kl]=[%g,%g] != [kappa,thermalConductivity]=[%g,%g] \n",Dl,Kl,kappa[0],thermalConductivity);
        OV_ABORT("ERROR");
      }
      // alpha = sqrt( ss/Dl + kypi*kypi ); 
      c1a = c1; c2a= c2;
      alphap = alphaLp;
      alpham = alphaLm;
    }
    else if( option==1 )
    {
      // Right domain -- choose coefficients to match the jump conditions:
      //    Tl(0) = Tr(0)
      //   Kl*Tl.x(0)  = Kr*Tr.x(0)
      if( Dr != kappa[0] || Kr != thermalConductivity )
      {
        printF(">>>> CGAD:UDKS: Plane interface mode:ERROR: right domain: [Dr,Kr]=[%g,%g] != [kappa,thermalConductivity]=[%g,%g] \n",Dr,Kr,kappa[0],thermalConductivity);        
        OV_ABORT("ERROR");
      }
      // ****  see formula in CgmpRef/tex/planeInterface.tex ******

      // we need c1a and c2a to satisfy:
      //    c1a + c2a = c1 + c2
      //  Kr*alphaR(c1a-c2a) = Kl*alphaL( c1-c2) 
      // std::complex<Real> alphaL = sqrt( ss/Dl + kypi*kypi ); 
      // std::complex<Real> alphaR = sqrt( ss/Dr + kypi*kypi ); 
      // std::complex<Real> beta   = (Kl*alphaL)/(Kr*alphaR); 

      // c1a = .5*( (1.+beta)*c1 + (1.-beta)*c2 );
      // c2a = .5*( (1.-beta)*c1 + (1.+beta)*c2 );

      // alpha = alphaR;
      std::complex<Real> gamma1 = (Kl/Kr)*(alphaLp-uDL);
      std::complex<Real> gamma2 = (Kl/Kr)*(alphaLm-uDL);
      std::complex<Real> a2 = alphaRm-uDR;

      c1a = ( (gamma1-a2)*c1 + (gamma2-a2)*c2 )/( alphaRp-alphaRm );
      c2a = c1 + c2 - c1a;
      alphap = alphaRp;
      alpham = alphaRm;

      // --- check solution: -----
      Real xd=0., yd=.3;

      Real uL  = std::real( ( c1 *exp(alphaLp*xd) + c2 *exp(alphaLm*xd) )* exp( I*kypi*yd ) * expst );
      Real uR  = std::real( ( c1a*exp(alphaRp*xd) + c2a*exp(alphaRm*xd) )* exp( I*kypi*yd ) * expst );

      Real uxL = std::real( ( c1 *alphaLp*exp(alphaLp*xd) + c2 *alphaLm*exp(alphaLm*xd) )* exp( I*kypi*yd ) * expst );
      Real uxR = std::real( ( c1a*alphaRp*exp(alphaRp*xd) + c2a*alphaRm*exp(alphaRm*xd) )* exp( I*kypi*yd ) * expst );      
      
      Real jumpT = uL - uR;
      Real jumpFlux = Kl*( uxL - uDL*uL ) -
                      Kr*( uxR - uDR*uR );

      // T.t + u*T.x + v*T.y = D*( T.xx + T.yy )
      // s/D + uD*z + i*k*vD = z^2 - k^2 
      //   z^2 - uD*z -ik*vD - s/D - ky^2  = 0 
      Real resLp = abs(  SQR(alphaLp) - uDL*alphaLp - I*kypi*vDL - SQR(kypi) - ss/Dl );
      Real resLm = abs(  SQR(alphaLm) - uDL*alphaLm - I*kypi*vDL - SQR(kypi) - ss/Dl );

      Real resRp = abs(  SQR(alphaRp) - uDR*alphaRp - I*kypi*vDR - SQR(kypi) - ss/Dr );
      Real resRm = abs(  SQR(alphaRm) - uDR*alphaRm - I*kypi*vDR - SQR(kypi) - ss/Dr );


      printf("CHECK: t=%9.3e:  [T]=%8.2e, [K(T.x-uT)]=%8.2e, resL=%8.2e, resR=%8.2e\n",t,jumpT,jumpFlux,max(resLp,resLm),max(resRp,resRm));

      // OV_ABORT("UDKS:PIM: stop here for now");


    }

    if( t<=dt )
    {
      if( option==0 )
      {
        printF(">>>> CGAD:UDKS: Plane interface mode:LEFT: amp=%g, [Dl,Kl]=[%g,%g], [sr,si]=[%g,%g], ky=%g, c1a=[%g,%g], c2a=[%g,%g]\n",
                amp,Dl,Kl,sr,si,ky, std::real(c1a),std::imag(c1a), std::real(c2a),std::imag(c2a) );
      }
      else
      {
        printF(">>>> CGAD:UDKS: Plane interface mode:RIGHT: amp=%g, [Dr,Kr]=[%g,%g], [sr,si]=[%g,%g], ky=%g, c1a=[%g,%g], c2a=[%g,%g]\n",
                  amp,Dr,Kr,sr,si,ky, std::real(c1a),std::imag(c1a), std::real(c2a),std::imag(c2a) );
      }
      printF(">>>  alphap=[%g,%g] alpham=[%g,%g]\n",std::real(alphap),std::imag(alphap),std::real(alpham),std::imag(alpham));
      printF(">>>> CGAD:UDKS: [uL,vL]=[%g,%g], [uR,vR]=[%g,%g] (advection velocities)\n",uL,vL,uR,vR);


      // if( (uL !=0. || uR!=0.) )
      // {
      //   OV_ABORT("stop here for now");
      // }
    }

    if( numberOfDimensions==2 )
    {

      FOR_3D(i1,i2,i3,I1,I2,I3)
      { 
        const Real xd = x(i1,i2,i3,0), yd = x(i1,i2,i3,1);
        ua(i1,i2,i3,nn) = std::real( ( c1a*exp(alphap*xd) + c2a*exp(alpham*xd) )* exp( I*kypi*yd ) * expst );
      }
    }
    else
    { // --- 3D ----
      OV_ABORT("finish me");
    }


  }

  else 
  {
    // look for a solution in the base class
    Parameters::getUserDefinedKnownSolution( t, cg, grid, ua, I1,I2,I3 );
  }
  

  
  return 0;
}




int AdParameters::
updateUserDefinedKnownSolution( GenericGraphicsInterface & gi, CompositeGrid & cg)
// ==========================================================================================
/// \brief This function is called to set the user defined known solution.
/// 
/// \return   >0 : known solution was chosen, 0 : no known solution was chosen
///
// ==========================================================================================
{

  // KnownSolutionsEnum & knownSolution = dbase.get<Parameters::KnownSolutionsEnum>("knownSolution");

  // Make  dbase.get<Real >("a") sub-directory in the data-base to store variables used here
  if( ! dbase.get<DataBase >("modelData").has_key("userDefinedKnownSolutionData") )
     dbase.get<DataBase >("modelData").put<DataBase>("userDefinedKnownSolutionData");
  DataBase & db =  dbase.get<DataBase >("modelData").get<DataBase>("userDefinedKnownSolutionData");

  if( !db.has_key("userKnownSolution") )
  {
    db.put<aString>("userKnownSolution");
    db.get<aString>("userKnownSolution")="unknownSolution";
    
    db.put<Real[20]>("rpar");
    db.put<int[20]>("ipar");
  }
  aString & userKnownSolution = db.get<aString>("userKnownSolution");
  Real *rpar = db.get<Real[20]>("rpar");
  int *ipar = db.get<int[20]>("ipar");


  const aString menu[]=
    {
      "no known solution",
      "choose a common known solution",
      "square eigenfunction",
      "disk eigenfunction",
      "annulus eigenfunction",
      "plane interface mode",
      "concentric cylinders",
      "double annulus",
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
    else if( answer=="no known solution" )
    {
      userKnownSolution="unknownSolution";
    }

    else if( answer=="choose a common known solution" )
    {
      // Look for a known solution from the base class (in common/src)
      Parameters::updateUserDefinedKnownSolution(gi,cg);
    }

    else if( answer=="square eigenfunction" )
    {
      userKnownSolution="squareEigenfunction";

      int & bcOpt = ipar[0]; 

      Real & kx = rpar[0];
      Real & ky = rpar[1];
      Real & kz = rpar[2];

      kx=1., ky=1., kz=1.; 
      
      printF("--- Eigenfunction for the Heat Equation in a Square Domain ---\n"
             "    u = e^{- b*t } * sin( pi*kx*x ) * sin( pi*ky*y )  [2D, Dirichlet BCs]\n"
             "    u = e^{- b*t } * cos( pi*kx*x ) * cos( pi*ky*y )  [2D, Neumann BCs]\n"
             "       b = D*(pi)^2*( kx^2 + ky^2 )\n"
             "    u = e^{- b*t } * sin( pi*kx*x ) * sin( pi*ky*y ) *sin( pi*kz*z ) [3D, Dirichlet BCs]\n"
            "         b = D*(pi)^2*( kx^2 + ky^2 + kz^2 )\n"
             "    bcOpt : 0=Dirichlet, 1=Neumann BCs on the annulus.\n"               
        );
      
      gi.inputString(answer,"Enter kx,ky,kz,bcOpt for the exact solution");
      sScanF(answer,"%e %e %e %i",&kx,&ky,&kz,&bcOpt);

      printF("Square eigenfuncton: kx=%g, ky=%g, kz=%g, bcOpt=%d\n",kx,ky,kz,bcOpt);

      dbase.get<bool>("knownSolutionIsTimeDependent")= true;  // known solution is time dependent

    }

    else if( answer=="disk eigenfunction" )
    {
      userKnownSolution="diskEigenfunction";

      int & n     = ipar[0];
      int & m     = ipar[1];
      int & bcOpt = ipar[2];      

      Real & a   = rpar[0];
      Real & amp = rpar[1];

      m=0; n=0; amp=1.; a=1.;
      
      printF("--- Eigenfunction for the Heat Equation in the unit disk ---\n"
             "    u = amp*e^{- lambda^2*t } * Jn(lambda_m*r) * cos(n*theta)  [2D, Dirichlet BCs]\n"
             "    n = angular number, n=0,1,2,\n"
             "    m = radial number (m'th zero of Bessel Jn(lambda*a)=0 \n"
             "    a = radius of the disk\n"
             "    amp = amplitude\n"
             "    bcOpt : 0=Dirichlet, 1=Neumann BCs on the annulus.\n"             
        );
      
      gi.inputString(answer,"Enter n,m,a,amp,bcOpt for the exact solution");
      sScanF(answer,"%i %i %e %e %i",&n,&m,&a,&amp,&bcOpt);

      printF("Disk eigenfuncton: n=%i, m=%i, a=%g, amp=%g, bcOpt=%d\n",n,m,a,amp,bcOpt);

      dbase.get<bool>("knownSolutionIsTimeDependent")= true;  // known solution is time dependent

    }  

    else if( answer=="annulus eigenfunction" )
    {
      userKnownSolution="annulusEigenfunction";

      int & n     = ipar[0];
      int & m     = ipar[1];
      int & bcOpt = ipar[2];

      Real & amp = rpar[0];

      m=0; n=0; amp=1.; bcOpt=0;
      
      printF("--- Eigenfunction for the Heat Equation in an annulus---\n"
             "    u = amp*e^{- lambda^2*t } * uHat(r) * cos(n*theta),\n"
             "    uHat(r) = c1*Jn(lambda*r) + c2*Yn(lambda*r),\n"
             "    Annulus: ra=.5, rb=1.\n"
             "    n = angular number, n=0,1,2,\n"
             "    m = radial number (m'th zero of determinant condition). \n"
             "    amp = amplitude.\n"
             "    bcOpt : 0=Dirichlet, 1=Neumann BCs on the annulus.\n"
        );
      
      gi.inputString(answer,"Enter n,m,amp,bcOpt for the exact solution");
      sScanF(answer,"%i %i %e %i",&n,&m,&amp,&bcOpt);

      printF("Annulus eigenfuncton: n=%i, m=%i, amp=%g, bcOpt=%d\n",n,m,amp,bcOpt);

      dbase.get<bool>("knownSolutionIsTimeDependent")= true;  // known solution is time dependent

    } 

    else if( answer=="concentric cylinders" )
    {
      userKnownSolution="concentricCylinders";

      int & n     = ipar[0];
      int & m     = ipar[1];
      int & option= ipar[2];

      Real & amp = rpar[0];
      Real & a   = rpar[1];
      Real & b   = rpar[2];
      Real & c   = rpar[3];
      Real & D1  = rpar[4];
      Real & K1  = rpar[5];
      Real & D2  = rpar[6];
      Real & K2  = rpar[7];      


      printF("--- Concentric cylinders (or annulii in 2D) eigenfunction for the Heat Equation---\n"
             "  Solution is u = amp * [ c1J*Jn(lam1*r) + c1Y*Yn(lam1*r) ] cos(n*theta) exp(- s*t ),   a < r < b\n"
             "              u = amp * [ c2J*Jn(lam2*r) + c2Y*Yn(lam2*r) ] cos(n*theta) exp(- s*t ) ,   b < r < c\n"       
             "    n = angular number, n=0,1,2,\n"
             "    m = radial number, m=,1,2,3  (m'th zero of determinant condition for given n). \n"
             "    amp = amplitude.\n"
             "    a,b,c : inner, middle and outer radii\n"
             "    D1, K2 : diffusivity and thermal thermalConductivity (inner)\n"
             "    D2, K2 : diffusivity and thermal thermalConductivity (outer)\n"
             "    option : 0=left-side, 1=right-side.\n"   
             " Note: eigenvalues have only been computed for certain parameter values.\n"          
        );
      
      n=0; m=1; amp=1; a=.5; b=1.; c=1.5; D1=.1; K1=1.; D2=.05; K2=2.; option=0;

      gi.inputString(answer,"Enter n, m, amp, a, b, c, D1,K1, D2, K2, option, for the exact solution");
      sScanF(answer,"%i %i %e %e %e %e %e %e %e %e %i",&n, &m, &amp, &a, &b, &c, &D1,&K1, &D2,&K2,  &option);

      printF("Concentric cylinders eigenfunction: n=%i, m=%i, amp=%g, [D1,K1]=[%g,%g], [D2,K2]=[%g,%g] [a,b,c]=[%g,%g,%g] option=%d\n",n,m,amp,D1,K1,D2,K2,a,b,c,option);

      dbase.get<bool>("knownSolutionIsTimeDependent")= true;  // known solution is time dependent

    } 

    else if( answer=="double annulus" )
    {
      userKnownSolution="doubleAnnulus";

      printF("--- Double Annulus CHT solution - supports advection diffusion ---\n");
      printF("  caseName=doubleAnnulusCase5 : no advection.\n");

      aString caseName="doubleAnnulusCase5";    
      gi.inputString(answer,"Enter caseName");
      caseName = answer;
      printF("Setting caseName=[%s]\n",(const char*)caseName); 

      if( !dbase.has_key("adExactSolution") )
      {
        AdExactSolutions & adExactSolution = dbase.put<AdExactSolutions>("adExactSolution");

        adExactSolution.initialize( cg, "doubleAnnulus", caseName );
      }

      dbase.get<bool>("knownSolutionIsTimeDependent") = true;  // known solution is time dependent

    } 

    else if( answer=="double rectangle" )
    {
      userKnownSolution="doubleRectangle";

      printF("--- Double Rectangle solution - supports advection diffusion ---\n");
      printF("  caseName=doubleRectangleDL1KL1DR1KR1 : no advection.\n");

      aString caseName="doubleRectangleDL1KL1DR1KR1";    
      gi.inputString(answer,"Enter caseName");
      caseName = answer;
      printF("Setting caseName=[%s]\n",(const char*)caseName); 

      if( !dbase.has_key("adExactSolution") )
      {
        AdExactSolutions & adExactSolution = dbase.put<AdExactSolutions>("adExactSolution");

        adExactSolution.initialize( cg, "doubleRectangle", caseName );
      }

      dbase.get<bool>("knownSolutionIsTimeDependent") = true;  // known solution is time dependent

    } 


    else if( answer=="plane interface mode" )
    {
      userKnownSolution="planeInterfaceMode";

      int & option= ipar[0];

      Real & amp = rpar[0];
      Real & sr  = rpar[1];
      Real & si  = rpar[2];
      Real & ky  = rpar[3];
      Real & c1  = rpar[4];
      Real & c2  = rpar[5];
      Real & Dl  = rpar[6];
      Real & Kl  = rpar[7];
      Real & Dr  = rpar[8];
      Real & Kr  = rpar[9];

      // Advection velocities: 
      Real & uL  = rpar[10];
      Real & vL  = rpar[11];
      Real & wL  = rpar[12];

      Real & uR  = rpar[13];
      Real & vR  = rpar[14];
      Real & wR  = rpar[15];     

      option=0; amp=1; sr=-.1; si=1.; ky=2; c1=1; c2=1; Dl=1; Kl=1.; Dr=1.; Kr=1.;
      uL=0.; vL=0.; wL=0.; uR=0.; vR=0.; wR=0; 

      
      printF("--- Plane interface mode for the Heat equation on two adjacent squares---\n"
             "    u = amp*e^{s*t} * uHat(x) * e^{i*pi*ky*y} (form of solution on each side),\n"
             "    uHat(x) = c1*exp( alpha*x) + c2*exp( alpha*x ),\n"
             "    s = sr + I*si (Laplace transform parameter, input)\n"
             "    ky = Fourier mode in Y (input)\n"
             "    c1,c2 : coefficients in solution, differ on two sides of the interface. \n"
             "    Dl, Kl : diffusivity and thermal thermalConductivity on left\n"
             "    Dr, Kr : diffusivity and thermal thermalConductivity on right\n"
             "    uL,vL,wL : advection velocity on left\n"
             "    uR,vR,wR : advection velocity on right\n"
             "    option : 0=left-side, 1=right-side.\n"
        );
      
      gi.inputString(answer,"Enter option, amp,sr,si,ky,c1,c2, Dl,Kl, Dr, Kr, uL,vL,wL, uR,vR, wR  for the exact solution");
      sScanF(answer,"%i %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e ",&option,&amp,&sr,&si,&ky,&c1,&c2, &Dl,&Kl, &Dr,&Kr,
                     &uL,&vL,&wL, &uR,&vR, &wR);

      printF("Plane interface mode: Setting amp=%g, [sr,si]=[%g,%g] ky=%g, [c1,c2]=[%g,%g], [Dl,Kl]=[%g,%g], [Dr,Kr]=[%g,%g] option=%d\n",
               amp,sr,si,ky,c1,c2,Dl,Kl,Dr,Kr,option);
      printF("Plane interface mode: [uL,vL,wL]=[%g,%g,%g], [uR,vR,wR]=[%g,%g,%g] (advection velocities)\n",uL,vL,wL,uR,vR,wR);

      dbase.get<bool>("knownSolutionIsTimeDependent")= true;  // known solution is time dependent

      // std::vector<Real> & a = dbase.get<std::vector<Real> >("a");
      // std::vector<Real> & b = dbase.get<std::vector<Real> >("b");
      // std::vector<Real> & c = dbase.get<std::vector<Real> >("c"); 
      // if( option==0 )
      // {
      //   uL = a[0]; vL = b[0]; wL = c[0];  // advection velocity for the left side
      // }
      // else
      // {
      //   uR = a[0]; vR = b[0]; wR = c[0];  // advection velocity for the right side
      // }

      // // --- Find some parameters from the opposite side ---
      // const IntegerArray & interfaceType = dbase.get<IntegerArray >("interfaceType");
      // for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
      // {
      //   MappedGrid & mg = cg[grid];
      //   ForBoundary(side,axis)
      //   {
      //     if( interfaceType(side,axis,grid) == Parameters::heatFluxInterface && mg.boundaryCondition(side,axis)>0 )
      //     {
      //        Parameters & params2 = getInterfaceParameters( grid,side,axis,*this );
      //        std::vector<Real> & a2 = params2.dbase.get<std::vector<Real> >("a");
      //        std::vector<Real> & b2 = params2.dbase.get<std::vector<Real> >("b");
      //        std::vector<Real> & c2 = params2.dbase.get<std::vector<Real> >("c"); 

      //        if( option==0 )
      //        {
      //          uR = a2[0]; vR = b2[0]; wR = c2[0];  // advection velocity for the right side
      //        }
      //        else
      //        {
      //          uL = a2[0]; vL = b2[0]; wL = c2[0];  // advection velocity for the left side
      //        }
      //     }
      //   }
      // }

    }    

    else
    {
      printF("unknown response=[%s]\n",(const char*)answer);
      gi.stopReadingCommandFile();
    }
    
  }

  gi.unAppendTheDefaultPrompt();
  bool knownSolutionChosen = userKnownSolution!="unknownSolution";

  // if( knownSolutionChosen )
  //   knownSolution = Parameters::userDefinedKnownSolution;  // added Nov 27, 2021 so cgmp will plot errors

  return knownSolutionChosen;
}
