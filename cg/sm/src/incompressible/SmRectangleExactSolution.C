// This file automatically generated from SmRectangleExactSolution.bC with bpp.
// ==========================================================================================================
// Class to define exact solutions for CgSm
// 
//   VIBRATIONAL MODES OF A RECTANGLE IN 2D OR 3
//      (1) PERIODIC STRIP IN 2D
//      (1) PERIODIC STRIP IN 3D
// 
// ==========================================================================================================


#include "SmRectangleExactSolution.h"

#include "NurbsMapping.h"

#include "PlotStuff.h"
#include "ParallelUtility.h"

#include <complex>


typedef ::real LocalReal;
// typedef ::real OV_real;

typedef std::vector<std::complex<Real> > ComplexVector;

#define FOR_3D(i1,i2,i3,I1,I2,I3)                                       int I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase(); int I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(); for(i3=I3Base; i3<=I3Bound; i3++)                                       for(i2=I2Base; i2<=I2Bound; i2++)                                     for(i1=I1Base; i1<=I1Bound; i1++)



// ===============================================================================
/// \brief  Constructor for the class that defines exact solutions for a Rectangle
// ===============================================================================
SmRectangleExactSolution::
SmRectangleExactSolution()
{

  // dbase.put<int>("numberOfDimensions");
    dbase.put<bool>("initialized")=false;

    dbase.put<Real>("omega") = -1.;
    dbase.put<aString>("caseName") = "unknown";

  // mode numbers 
    dbase.put<int>("mx")=1;
    dbase.put<int>("my")=1;
    dbase.put<int>("mz")=1;

}


// ===============================================================================
/// \brief destructor 
// ===============================================================================
SmRectangleExactSolution::
~SmRectangleExactSolution()
{
    if( dbase.has_key("nurbs") )
    {
        NurbsMapping *nurbs = dbase.get<NurbsMapping*>("nurbs");
        delete [] nurbs;
    }
    

}


// ===============================================================================
/// \brief Return a solution parameter by name
// ===============================================================================
int SmRectangleExactSolution::
getParameter( const aString & name, Real & value )
{
    if( dbase.has_key(name) )
    {
        value = dbase.get<Real>(name);
    // printF("\n ***** SmRectangleExactSolution:: omega=%g ****\n",value);
    }
    else
    {
        printF("SmRectangleExactSolution::getParameter: ERROR: unknown parameter=[%s]\n",(const char*)name);
    }
    return 0;
}

// ===============================================================================
/// \brief Set an integer parameter value by name
///   Common parameters: mx, my, mz
// ===============================================================================
int SmRectangleExactSolution::
setParameter( const aString & name, const int value )
{
    if( dbase.has_key(name) )
    {
        dbase.get<int>(name) = value;
    // printF("\n ***** SmRectangleExactSolution:: omega=%g ****\n",value);
    }
    else
    {
        printF("SmRectangleExactSolution::setParameter: ERROR: unknown parameter=[%s]\n",(const char*)name);
    }
    return 0;
}



// ===============================================================================
/// \brief Initialize the exact solution.
// ===============================================================================
int SmRectangleExactSolution::
initialize( CompositeGrid & cg, const aString & caseName )
{
    bool & initialized=dbase.get<bool>("initialized");
    initialized=true;

    dbase.get<aString>("caseName")=caseName;  // save the case name


    return 0;
}



//========================================================================================================
/// \brief Evaluate the exact solution
/// \param numberOfTimeDerivatives (input) : 
///      NOTE: numberOfTimeDerivatives=-1 means eval the pressure only.
//========================================================================================================
int SmRectangleExactSolution::
evalSolution(Real t, CompositeGrid & cg, int grid, RealArray & ua, 
                                                    const Index & I1, const Index &I2, const Index &I3, 
                                                    int numberOfTimeDerivatives /* = 0 */  )
{

    const aString & caseName = dbase.get<aString>("caseName");

 // The eigenvalue defined by (mx,my), in the range 1,2,3,...
    const int & mx = dbase.get<int>("mx");
    const int & my = dbase.get<int>("my");
    const int & mz = dbase.get<int>("mz");

    if( true )
        printF("\n +++++ SmRectangleExactSolution:: evalSolution caseName=%s, (mx,my,mz)=(%d,%d,%d) numberOfTimeDerivatives=%d, t=%9.3e ++++++\n",
                      (const char*)caseName,mx,my,mz,numberOfTimeDerivatives,t);

  

    Real rho =1., mu=1., c=1.; // FIX ME ****************

    MappedGrid & mg = cg[grid];
    mg.update(MappedGrid::THEcenter | MappedGrid::THEvertex );
    const int numberOfDimensions = mg.numberOfDimensions();

    const int u1c = 0;
    const int u2c = 1;
    const int u3c = 2;
    int pc        = numberOfDimensions;       // may change below if we only eval the pressure

    OV_GET_SERIAL_ARRAY(real,mg.center(),xLocal);
    OV_GET_SERIAL_ARRAY(real,ua,uLocal);


    Real A1,B1, A2,B2, A3, B3, P1, P2, P1p, P1m, P2p, P2m, P3p, P3m;
    Real omega, beta, alpha, kx, ky, kz, k;

    Real scaleFactor=1.; // sets the scale 
    bool isSurfaceWave=false;

  // ---> File written by AMP/ism/maple/ismStrip2d.mpl
  //
    if( caseName=="periodicStripTT" || caseName=="rectangleTT" )
    {
    // periodic strip BC=TT

    // --- Eigenvalues of a 2D periodic strip for incompressible elasticity -----
    //   BC = TT
    //   File written by ismStrip2d.mw on 2022-06-20
    //   Solution is of the form  : u(x,y,t) = uHat(x) *exp( I*k*y )*exp( I*omega t) , k=2*Pi*m 
    //   uHat1 =  A1*cos(beta*x) + B1*sin(beta*x) + pHat'(x)/(mu*(beta^2+k^2) 
    //   uHat2 =  from divergence
    //   pHat  =  A3*exp(k*x) + B3*exp(-k*x) 
    //   beta^2 = (omega/c)^2 - k^2 
    //   Note that beta is generally real except for surface waves where it is pure imaginary
        if( numberOfDimensions==2 )
        {         
// --- Eigenvalues of a 2D periodic strip for incompressible elasticity -----
//   BC = TT
//   File written by ismStrip2d.mw on 2022-06-21
//   Solution is of the form  : u(x,y,t) = uHat(x) *exp( I*k*y )*exp( I*omega t) , k=2*Pi*m 
//   uHat1 =  A1*cos(beta*x) + B1*sin(beta*x) + pHat'(x)/(mu*(beta^2+k^2) 
//   uHat2 =  from divergence
//   pHat  =  A3*exp(k*x) + B3*exp(-k*x) 
//   beta^2 = (omega/c)^2 - k^2 
//   Note that beta is generally real except for surface waves where it is pure imaginary 
const int mxd=4, myd=3; // max number of roots in the table
const Real omegaArray[]={
    // m=1, n=1,2,..,4
        5.7766060906888207e+00, 6.6005868537057281e+00, 8.8857658763167325e+00, 1.1003069256654355e+01,
    // m=2, n=1,2,..,4
        1.2106218072213127e+01, 1.3753223189861312e+01, 1.5624812844669151e+01, 1.7771531752633465e+01,
    // m=3, n=1,2,..,4
        1.8027501109917280e+01, 1.9472778749087922e+01, 2.0798488827365205e+01, 2.2516638126074782e+01
};
#define omegaRoot(mx,my) omegaArray[(mx-1)+mxd*(my-1)]
omega = omegaRoot(mx,my);
ky = twoPi*my; kx=ky; k=kx;
B3 = scaleFactor;  // scaleFactor should be set before this point 
Real betaSq = SQR(omega/c) - SQR(kx);
if( betaSq >0 )
{
    isSurfaceWave=false;
    beta = sqrt( betaSq );
    A1 = -0.4e1 / (0.4e1 * exp(k) * pow(beta, 0.5e1) * pow(k, 0.3e1) - 0.4e1 * exp(k) * pow(k, 0.7e1) * beta - sin(beta) * pow(beta, 0.8e1) + 0.2e1 * sin(beta) * pow(beta, 0.6e1) * k * k - 0.2e1 * sin(beta) * beta * beta * pow(k, 0.6e1) + sin(beta) * pow(k, 0.8e1) - 0.4e1 * cos(beta) * pow(beta, 0.5e1) * pow(k, 0.3e1) + 0.4e1 * cos(beta) * beta * pow(k, 0.7e1)) * pow(k, 0.3e1) * (0.2e1 * pow(k, 0.3e1) * exp(k) * beta - sin(beta) * pow(beta, 0.4e1) + 0.2e1 * sin(beta) * beta * beta * k * k - sin(beta) * pow(k, 0.4e1) - 0.2e1 * pow(k, 0.3e1) * exp(-k) * beta) * B3;
    A3 = (sin(beta) * pow(beta, 0.4e1) - 0.2e1 * sin(beta) * beta * beta * k * k + sin(beta) * pow(k, 0.4e1) - 0.4e1 * cos(beta) * beta * pow(k, 0.3e1) + 0.4e1 * pow(k, 0.3e1) * exp(-k) * beta) * B3 / (0.4e1 * pow(k, 0.3e1) * exp(k) * beta - sin(beta) * pow(beta, 0.4e1) + 0.2e1 * sin(beta) * beta * beta * k * k - sin(beta) * pow(k, 0.4e1) - 0.4e1 * cos(beta) * beta * pow(k, 0.3e1));
    B1 = 0.2e1 * (exp(k) * beta * beta - exp(k) * k * k - 0.2e1 * cos(beta) * beta * beta + 0.2e1 * cos(beta) * k * k + exp(-k) * beta * beta - exp(-k) * k * k) * pow(k, 0.3e1) * B3 / (beta * beta + k * k) / (0.4e1 * pow(k, 0.3e1) * exp(k) * beta - sin(beta) * pow(beta, 0.4e1) + 0.2e1 * sin(beta) * beta * beta * k * k - sin(beta) * pow(k, 0.4e1) - 0.4e1 * cos(beta) * beta * pow(k, 0.3e1));
  // uHat2 = -uHat1.x/ky;
    A2 = -beta*B1/ky;
    B2 =  beta*A1/ky;
}
else
{
    isSurfaceWave=true;
    alpha = sqrt( -betaSq );
    A1 = 0.4e1 * B3 * (0.2e1 * pow(k, 0.3e1) * exp(-k) * alpha + pow(alpha * alpha + k * k, 0.2e1) * sinh(alpha) - 0.2e1 * pow(k, 0.3e1) * exp(k) * alpha) * pow(k, 0.3e1) / (alpha * alpha + k * k) / (alpha + k) / (pow(alpha * alpha + k * k, 0.2e1) * sinh(alpha) - 0.4e1 * pow(k, 0.3e1) * alpha * (exp(k) - cosh(alpha))) / (k - alpha);
    B1 = 0.2e1 * (alpha * alpha + k * k) * (exp(k) - 0.2e1 * cosh(alpha) + exp(-k)) * pow(k, 0.3e1) * B3 / (alpha + k) / (pow(alpha * alpha + k * k, 0.2e1) * sinh(alpha) - 0.4e1 * pow(k, 0.3e1) * alpha * (exp(k) - cosh(alpha))) / (k - alpha);
    A3 = -B3 * (0.4e1 * pow(k, 0.3e1) * exp(-k) * alpha + pow(alpha * alpha + k * k, 0.2e1) * sinh(alpha) - 0.4e1 * cosh(alpha) * alpha * pow(k, 0.3e1)) / (pow(alpha * alpha + k * k, 0.2e1) * sinh(alpha) - 0.4e1 * pow(k, 0.3e1) * alpha * (exp(k) - cosh(alpha)));
  // uHat2 = -uHat1.x/ky;
    A2 = -alpha*B1/ky;
    B2 = -alpha*A1/ky;
}
P1p=  A3*kx/(rho*SQR(omega));
P1m= -B3*kx/(rho*SQR(omega));
P2p= -kx*P1p/ky;
P2m= +kx*P1m/ky;
        }
        else
        {
            OV_ABORT("Finish me");
        }      

    }
    else if( caseName=="periodicStripTD" || caseName=="rectangleTD" )
    {
    // periodic strip BC=TRACTION - DISPLACEMENT 

        printF(" \n @@@@@@@ case TD @@@@ \n");  
        if( numberOfDimensions==2 )
        {
// --- Eigenvalues of a 2D periodic strip for incompressible elasticity -----
//   BC = TD
//   File written by ismStrip2d.mw on 2022-06-21
//   Solution is of the form  : u(x,y,t) = uHat(x) *exp( I*k*y )*exp( I*omega t) , k=2*Pi*m 
//   uHat1 =  A1*cos(beta*x) + B1*sin(beta*x) + pHat'(x)/(mu*(beta^2+k^2) 
//   uHat2 =  from divergence
//   pHat  =  A3*exp(k*x) + B3*exp(-k*x) 
//   beta^2 = (omega/c)^2 - k^2 
//   Note that beta is generally real except for surface waves where it is pure imaginary 
const int mxd=4, myd=3; // max number of roots in the table
const Real omegaArray[]={
    // m=1, n=1,2,..,4
        6.1046733125517700e+00, 8.2872475749655039e+00, 1.0694170271326492e+01, 1.3041586135065603e+01,
    // m=2, n=1,2,..,4
        1.2008573987474925e+01, 1.3377455285429044e+01, 1.5071190101197233e+01, 1.7199046325461177e+01,
    // m=3, n=1,2,..,4
        1.8007351701186442e+01, 1.9295616953990335e+01, 2.0422956693880925e+01, 2.2014459535174920e+01
};
#define omegaRoot(mx,my) omegaArray[(mx-1)+mxd*(my-1)]
omega = omegaRoot(mx,my);
ky = twoPi*my; kx=ky; k=kx;
B3 = scaleFactor;  // scaleFactor should be set before this point 
Real betaSq = SQR(omega/c) - SQR(kx);
if( betaSq >0 )
{
    isSurfaceWave=false;
    beta = sqrt( betaSq );
    A1 = -0.4e1 / (sin(beta) * pow(beta, 0.6e1) - sin(beta) * pow(beta, 0.4e1) * k * k - sin(beta) * beta * beta * pow(k, 0.4e1) + sin(beta) * pow(k, 0.6e1) + 0.4e1 * cos(beta) * pow(beta, 0.3e1) * pow(k, 0.3e1) + 0.4e1 * cos(beta) * beta * pow(k, 0.5e1) + 0.2e1 * exp(k) * pow(beta, 0.5e1) * k - 0.2e1 * exp(k) * beta * pow(k, 0.5e1)) * pow(k, 0.3e1) * (sin(beta) * beta * beta - sin(beta) * k * k + k * exp(k) * beta - k * exp(-k) * beta) * B3;
    A3 = -(sin(beta) * pow(beta, 0.4e1) - 0.2e1 * sin(beta) * beta * beta * k * k + sin(beta) * pow(k, 0.4e1) - 0.4e1 * cos(beta) * beta * pow(k, 0.3e1) - 0.2e1 * exp(-k) * pow(beta, 0.3e1) * k + 0.2e1 * exp(-k) * beta * pow(k, 0.3e1)) * B3 / (sin(beta) * pow(beta, 0.4e1) - 0.2e1 * sin(beta) * beta * beta * k * k + sin(beta) * pow(k, 0.4e1) + 0.4e1 * cos(beta) * beta * pow(k, 0.3e1) + 0.2e1 * exp(k) * pow(beta, 0.3e1) * k - 0.2e1 * exp(k) * beta * pow(k, 0.3e1));
    B1 = (0.4e1 * cos(beta) * k * k + exp(k) * beta * beta - k * k * exp(k) + exp(-k) * beta * beta - k * k * exp(-k)) * k * B3 * (beta * beta - k * k) / (beta * beta + k * k) / (sin(beta) * pow(beta, 0.4e1) - 0.2e1 * sin(beta) * beta * beta * k * k + sin(beta) * pow(k, 0.4e1) + 0.4e1 * cos(beta) * beta * pow(k, 0.3e1) + 0.2e1 * exp(k) * pow(beta, 0.3e1) * k - 0.2e1 * exp(k) * beta * pow(k, 0.3e1));
  // uHat2 = -uHat1.x/ky;
    A2 = -beta*B1/ky;
    B2 =  beta*A1/ky;
}
else
{
    isSurfaceWave=true;
    alpha = sqrt( -betaSq );
    A1 = 0.4e1 * pow(k, 0.3e1) * B3 * (k * exp(-k) * alpha + (alpha * alpha + k * k) * sinh(alpha) - k * exp(k) * alpha) / (pow(alpha * alpha + k * k, 0.2e1) * sinh(alpha) - 0.2e1 * alpha * k * ((alpha * alpha + k * k) * exp(k) - 0.2e1 * cosh(alpha) * k * k)) / (k - alpha) / (alpha + k);
    B1 = k * ((alpha * alpha + k * k) * exp(-k) + (alpha * alpha + k * k) * exp(k) - 0.4e1 * cosh(alpha) * k * k) * B3 * (alpha * alpha + k * k) / ((-0.2e1 * pow(alpha, 0.3e1) * k - 0.2e1 * alpha * pow(k, 0.3e1)) * exp(k) + pow(alpha * alpha + k * k, 0.2e1) * sinh(alpha) + 0.4e1 * cosh(alpha) * alpha * pow(k, 0.3e1)) / (k - alpha) / (alpha + k);
    A3 = -((0.2e1 * pow(alpha, 0.3e1) * k + 0.2e1 * alpha * pow(k, 0.3e1)) * exp(-k) + pow(alpha * alpha + k * k, 0.2e1) * sinh(alpha) - 0.4e1 * cosh(alpha) * alpha * pow(k, 0.3e1)) * B3 / (pow(alpha * alpha + k * k, 0.2e1) * sinh(alpha) - 0.2e1 * alpha * k * ((alpha * alpha + k * k) * exp(k) - 0.2e1 * cosh(alpha) * k * k));
  // uHat2 = -uHat1.x/ky;
    A2 = -alpha*B1/ky;
    B2 = -alpha*A1/ky;
}
P1p=  A3*kx/(rho*SQR(omega));
P1m= -B3*kx/(rho*SQR(omega));
P2p= -kx*P1p/ky;
P2m= +kx*P1m/ky;
        }
        else
        {
            OV_ABORT("Finish me");
        }

    }  
    else if( caseName=="periodicStripDD" || caseName=="rectangleDD" )
    {
    // periodic strip BC=DISPLACEMENT - DISPLACEMENT 
        printF(" \n @@@@@@@ case DD @@@@ \n");  
        if( numberOfDimensions==2 )
        {    
// --- Eigenvalues of a 2D periodic strip for incompressible elasticity -----
//   BC = DD
//   File written by ismStrip2d.mw on 2022-06-21
//   Solution is of the form  : u(x,y,t) = uHat(x) *exp( I*k*y )*exp( I*omega t) , k=2*Pi*m 
//   uHat1 =  A1*cos(beta*x) + B1*sin(beta*x) + pHat'(x)/(mu*(beta^2+k^2) 
//   uHat2 =  from divergence
//   pHat  =  A3*exp(k*x) + B3*exp(-k*x) 
//   beta^2 = (omega/c)^2 - k^2 
//   Note that beta is generally real except for surface waves where it is pure imaginary 
const int mxd=4, myd=3; // max number of roots in the table
const Real omegaArray[]={
    // m=1, n=1,2,..,4
        7.6465512957455341e+00, 1.0252755871737867e+01, 1.3170257593343554e+01, 1.6177532569133238e+01,
    // m=2, n=1,2,..,4
        1.3104495532075780e+01, 1.4553250335532520e+01, 1.6601752731685955e+01, 1.9009502048418872e+01,
    // m=3, n=1,2,..,4
        1.9173530688079666e+01, 2.0105183319631333e+01, 2.1545698437568630e+01, 2.3380429851173121e+01
};
#define omegaRoot(mx,my) omegaArray[(mx-1)+mxd*(my-1)]
omega = omegaRoot(mx,my);
ky = twoPi*my; kx=ky; k=kx;
B3 = scaleFactor;  // scaleFactor should be set before this point 
Real betaSq = SQR(omega/c) - SQR(kx);
if( betaSq >0 )
{
    isSurfaceWave=false;
    beta = sqrt( betaSq );
    A1 = 0.1e1 / (cos(beta) * pow(beta, 0.3e1) + cos(beta) * k * k * beta + beta * beta * k * sin(beta) + sin(beta) * pow(k, 0.3e1) - exp(k) * pow(beta, 0.3e1) - k * k * exp(k) * beta) * (0.2e1 * sin(beta) * k - beta * exp(k) + beta * exp(-k)) * k * B3;
    A3 = (beta * cos(beta) - sin(beta) * k - beta * exp(-k)) * B3 / (beta * cos(beta) + sin(beta) * k - beta * exp(k));
    B1 = -(-exp(k) - exp(-k) + 0.2e1 * cos(beta)) * k * k * B3 / (beta * beta + k * k) / (beta * cos(beta) + sin(beta) * k - beta * exp(k));
  // uHat2 = -uHat1.x/ky;
    A2 = -beta*B1/ky;
    B2 =  beta*A1/ky;
}
else
{
    isSurfaceWave=true;
    alpha = sqrt( -betaSq );
    A1 = (-alpha * exp(k) + alpha * exp(-k) + 0.2e1 * sinh(alpha) * k) * k * B3 / (-alpha * exp(k) + alpha * cosh(alpha) + sinh(alpha) * k) / (-alpha + k) / (alpha + k);
    B1 = -(exp(k) + exp(-k) - 0.2e1 * cosh(alpha)) * k * k * B3 / (-alpha * exp(k) + alpha * cosh(alpha) + sinh(alpha) * k) / (-alpha * alpha + k * k);
    A3 = (alpha * exp(-k) - alpha * cosh(alpha) + sinh(alpha) * k) / (alpha * exp(k) - alpha * cosh(alpha) - sinh(alpha) * k) * B3;
  // uHat2 = -uHat1.x/ky;
    A2 = -alpha*B1/ky;
    B2 = -alpha*A1/ky;
}
P1p=  A3*kx/(rho*SQR(omega));
P1m= -B3*kx/(rho*SQR(omega));
P2p= -kx*P1p/ky;
P2m= +kx*P1m/ky;
        }
        else
        {
            OV_ABORT("Finish me");
        }
    }    
    else
    {
        printF("SmRectangleExactSolution:ERROR: unknown caseName=[%s]\n",(const char*)caseName);
        printF(" Expected caseName=rectangleTT or caseName=periodicStripTT\n");
        OV_ABORT("ERROR");
    }


    Real & omegadb = dbase.get<Real>("omega");
    omegadb = omega;  // save 

    Real coswt, sinwt;
  // By default we eval the displacements and pressure
    bool evalDisplacements=true;
    bool evalPressure=true;

    if( numberOfTimeDerivatives==-1 )
    { // Only eval the pressure 
        evalDisplacements=false;
        numberOfTimeDerivatives=0;
        pc=0; // NOTE - save pressure only
    }

    if( numberOfTimeDerivatives==0 )
    {
        coswt = cos(omega*t);  
        sinwt = sin(omega*t);
    }
    else if( numberOfTimeDerivatives==1 )
    { // do not eval pressure in this case
        evalPressure=false;
        coswt = -omega*sin(omega*t); 
        sinwt =  omega*cos(omega*t);
    }
    else if( numberOfTimeDerivatives==2 )
    { // do not eval pressure in this case
        evalPressure=false;
        coswt = -omega*omega*cos(omega*t); 
        sinwt = -omega*omega*sin(omega*t);
    }  
    else
    {
        OV_ABORT("error");
    }

    Real cosbx, sinbx; 
    int i1,i2,i3;
    if( numberOfDimensions==2 )
    {
        FOR_3D(i1,i2,i3,I1,I2,I3)
        {
            Real x = xLocal(i1,i2,i3,0);
            Real y = xLocal(i1,i2,i3,1);

            if( isSurfaceWave )
            {
                cosbx = cosh(alpha*x); sinbx=sinh(alpha*x);
            }
            else
            {
                cosbx = cos(beta*x), sinbx=sin(beta*x);
            }

            const Real exppkx=exp(kx*x), expmkx=1./exppkx; 
            const Real cosky = cos(ky*y), sinky=sin(ky*y); 

      // u1 = Re( u1Hat * exp( i*ky*y) * exp(i*omega t) )
      // u2 = Im( u1Hat * exp( i*ky*y) * exp(i*omega t) )
      // p  = Re( pHat  * exp( i*ky*y) * exp(i*omega t) )

            if( evalDisplacements ) // eval u or its time derivative 
            {
                uLocal(i1,i2,i3,u1c) =  ( A1*cosbx + B1*sinbx + P1p*exppkx + P1m*expmkx )*(cosky*coswt - sinky*sinwt);
                uLocal(i1,i2,i3,u2c) =  ( A2*cosbx + B2*sinbx + P2p*exppkx + P2m*expmkx )*(cosky*sinwt + sinky*coswt);
            }
            if( evalPressure )
            {
                uLocal(i1,i2,i3,pc ) =  ( A3*exppkx + B3*expmkx                         )*(cosky*coswt - sinky*sinwt);
            }

        }
    }
    else
    {
        OV_ABORT("finish me");
        FOR_3D(i1,i2,i3,I1,I2,I3)
        {
            Real x = xLocal(i1,i2,i3,0);
            Real y = xLocal(i1,i2,i3,1);
            Real z = xLocal(i1,i2,i3,2);

            if( isSurfaceWave )
            {
                cosbx = cosh(alpha*x); sinbx=sinh(alpha*x);
            }
            else
            {
                cosbx = cos(beta*x), sinbx=sin(beta*x);
            }

            const Real exppkx=exp(kx*x), expmkx=1./exppkx; 
            const Real cosky = cos(ky*y +kz*z), sinky=sin(ky*y + kz*z); 

      // u1 = Re( u1Hat * exp( i*ky*y + i*kz*z ) * exp(i*omega t) )
      // u2 = Im( u1Hat * exp( i*ky*y + i*kz*z ) * exp(i*omega t) )
      // p  = Re( pHat  * exp( i*ky*y + i*kz*z ) * exp(i*omega t) )

            if( evalDisplacements ) // eval u or its time derivative 
            {
                uLocal(i1,i2,i3,u1c) =  ( A1*cosbx + B1*sinbx + P1p*exppkx + P1m*expmkx )*(cosky*coswt - sinky*sinwt);
                uLocal(i1,i2,i3,u2c) =  ( A2*cosbx + B2*sinbx + P2p*exppkx + P2m*expmkx )*(cosky*sinwt + sinky*coswt);
                uLocal(i1,i2,i3,u3c) =  ( A3*cosbx + B3*sinbx + P3p*exppkx + P3m*expmkx )*(cosky*sinwt + sinky*coswt);
            }
            if( evalPressure )
            {
                uLocal(i1,i2,i3,pc ) =  ( P1*exppkx + P2*expmkx                         )*(cosky*coswt - sinky*sinwt);
            }

        }    
    }


    return 0;
}


