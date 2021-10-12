// This file automatically generated from NonlinearExactSolution.bC with bpp.
#include "NonlinearExactSolution.h"

#include "DispersiveMaterialParameters.h"

#include "ParallelUtility.h"

// ===============================================================================
// Class to define exact solutions to Maxwell's equations for
//                 NONLINEAR MULTILEVEL ATOMIC MODELS
// 
// ===============================================================================

#define mbe1d EXTERN_C_NAME(mbe1d)

extern "C"
{

// void mbe1d(real & epsilon, real & dt, real & xa, real &xb, int &n, real &tfinal, real &un);
    void mbe1d(real & rpar, int & ipar, int & n, real &un);

}


#define FOR_3D(i1,i2,i3,I1,I2,I3)                                       int I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase(); int I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(); for(i3=I3Base; i3<=I3Bound; i3++)                                       for(i2=I2Base; i2<=I2Bound; i2++)                                     for(i1=I1Base; i1<=I1Bound; i1++)


// ===============================================================================
/// \brief  Constructor 
// ===============================================================================
NonlinearExactSolution::
NonlinearExactSolution( )
{

    dbase.put<int>("initialized")=0;
    dbase.put<aString>("caseName");
  // caseName:
  //     soliton : 1D soliton-like defined by the solution to a 1D fortran code
  //     asymptoticSoliton : asymptotic formula for a soliton
    aString & caseName = dbase.get<aString>("caseName");
    caseName = "soliton";

    
}



// ===============================================================================
/// \brief destructor for the class that defines exact solutions to Maxwell's equations for a sphere
// ===============================================================================
NonlinearExactSolution::
~NonlinearExactSolution()
{
  // std::complex<LocalReal> *& Vvc = dbase.get<std::complex<LocalReal>* >("Vvc"); 
  // delete [] Vvc;

  // std::complex<LocalReal> *& kvc = dbase.get<std::complex<LocalReal>* >("kvc"); 
  // delete [] kvc;

  // std::complex<LocalReal> *& rtc = dbase.get<std::complex<LocalReal>* >("rtc"); 
  // delete [] rtc;

}

//=========================================================================
/// \brief set the case name which defines the solution
/// caseName:
///     soliton : 1D soliton-like defined by the solution to a 1D fortran code
///     asymptoticSoliton : asymptotic formula for a soliton
//=========================================================================
int NonlinearExactSolution::setCase( const aString & caseName )
{ 
    dbase.get<aString>("caseName") = caseName;
    
    return 0;
}

// ------------------ Here is how to call Fortran routines -------
//- 
//- // lapack routines
//- #ifdef OV_USE_DOUBLE
//-   #define GESV EXTERN_C_NAME(zgesv)
//-   #define GEEV EXTERN_C_NAME(zgeev)
//- #else
//-   #define GESV EXTERN_C_NAME(cgesv)
//-   #define GEEV EXTERN_C_NAME(cgeev)
//- #endif
//- 
//- extern "C"
//- {
//-   /* Solve  A*X = B (complex) */
//-   void GESV( int & N, int & NRHS, std::complex<LocalReal>  & a, const int & lda, int & ipvt, std::complex<LocalReal> & b, int & LDB, int & info );
//- 
//-   void GEEV( char *jobvl,
//-          char* jobvr,
//-          int & n,
//-          std::complex<LocalReal> & a,
//-          const int & lda,
//-          std::complex<LocalReal> & w,
//-          std::complex<LocalReal> &vl,
//-          int & ldvl,
//-          std::complex<LocalReal> &vr,
//-          int & ldvr,
//-          std::complex<LocalReal> & work,
//-          int & lwork,
//-          LocalReal & rwork,
//-          int & info );
//- }
//- 

//==============================================================================
//==================================  INITIALIZE  ==============================
///
/// \brief Init the nonlinear solution.
///
//==============================================================================
int NonlinearExactSolution::
initialize( CompositeGrid & cg, int numberOfDomains,
                std::vector<DispersiveMaterialParameters> & dispersiveMaterialParameters,
                const real & omega, const RealArray & kvI, const RealArray & asymParams, const int solveForAllFields )
{

    int & initialized = dbase.get<int>("initialized");
    if( initialized )
    {
        printF("NonlinearExactSolution::initialize - WARNING: already initialized! Nothing to be done\n");
        return 0;
    }

    const aString & caseName = dbase.get<aString>("caseName");
    if( caseName != "soliton" && caseName !="asymptoticSoliton" )
    {
        printF("NonlinearExactSolution::initialize:ERROR: Unknown caseName=[%s]\n",(const char*)caseName);
        printF(" Valid caseNames are `soliton' and 'asymptoticSoliton'\n");
        OV_ABORT("error");
    }

  // initialize
    if( numberOfDomains==1 )
    {
    // ---- 1 domain problem
    // parameters for asymptotic solutions
        const real eps = 0.1;
        const real U = 1./2;
        const real eta = 1.;
        const real x0 = 0.;
        real epsHat = eps/2.*sqrt(eta/(U-U*U));

    // parameters for asymptotic solutions
        asymParams(0) = U;
        asymParams(1) = x0;
        asymParams(2) = epsHat;
        asymParams(3) = 2.*sqrt(eta*U/(1.-U)); // coefficients for the asymptotic solution of E

        printF("Nonlinear MLA Exact Solution: initialize for one domain. caseName=%s\n",(const char*)caseName);
    }
    else
    {
    // ---- 2 domain scattering problem (TODO)
        printF("Nonlinear MLA Exact Solution: initialize for two domains (TODO)\n");
    
    }
    
    initialized=1;

  // Keep a pointer to the vector of DispersiveMaterialParameters
    dbase.put<std::vector<DispersiveMaterialParameters>* >("pDispersiveMaterialParameters") = &dispersiveMaterialParameters;
    

  // OV_ABORT("NonlinearExactSolution::initialize *new*  stop here for now");
      

    return 0;
}





// ==========================================================================================
/// \brief  Evaluate the solution and save in an array.
///
/// \param numberOfTimeDerivatives (input) : evaluate this many time-derivatives of the solution.
/// \param computeMagneticField (input): if true return the magnetic field in 3D (in 2D the magnetic field is always computed). 
// ==========================================================================================
int NonlinearExactSolution::
eval(real dt, real t, CompositeGrid & cg, int grid, 
          realArray & ua, realArray & pv, realArray & qv,
          const Index & I1a, const Index &I2a, const Index &I3a, 
          int numberOfTimeDerivatives /* = 0 */,
          bool computeMagneticField /* = false */ )
{

  // domain number for this grid: 
    const int & numberOfDomains = cg.numberOfDomains();
    const int myDomain = cg.domainNumber(grid);

    const aString & caseName = dbase.get<aString>("caseName");
    if( t <= 0. )
    {
        printF("--NonlinearExactSolution--  eval %s on grid=%i, domain=%i, at t=%9.3e (numberOfDomains=%d) \n",
                  (const char*)caseName,grid,myDomain,t,numberOfDomains);
    }

  // const real & dt= deltaT;

    MappedGrid & mg = cg[grid];
    const int numberOfDimensions = cg.numberOfDimensions();
    const IntegerArray & dw = mg.discretizationWidth();
    
  // dw.display("dw");
    Range R = numberOfDimensions;
    const int minDiscretizationWidth=min(dw(R));

    const int orderOfAccuracyInSpace = minDiscretizationWidth-1;

    std::vector<DispersiveMaterialParameters> & dmpVector = *dbase.get<std::vector<DispersiveMaterialParameters>* >("pDispersiveMaterialParameters");

  //    parameters.dbase.put<std::vector<DispersiveMaterialParameters> >("materialRegionParameters");

    DispersiveMaterialParameters & dmp    = dmpVector[myDomain];
    const int numberOfPolarizationVectors = dmp.numberOfPolarizationVectors;
    RealArray & gdmPar                    = dmp.modelParameters;

  // --- Here are the parameters in the nonlinear model ---
    RealArray multilevelAtomicParams;
    int numberOfAtomicLevels=-1;
    dmp.getNonlinearParameters( numberOfAtomicLevels, multilevelAtomicParams );  

    if( true && t<=0. && dmp.isDispersiveMaterial() )
    {
    //   OV_ABORT("NonlinearExactSolution::eval -- finish me for dispersive materials");
        printF("NES:eval: grid=%d: numberOfPolarizationVectors=%d, numberOfAtomicLevels=%d\n",
                    grid,numberOfPolarizationVectors,numberOfAtomicLevels);
        ::display(multilevelAtomicParams,"multilevelAtomicParams");
    }
    
    
    OV_GET_SERIAL_ARRAY(real,ua,uLocal);

    Index I1=I1a, I2=I2a, I3=I3a;
    bool ok = ParallelUtility::getLocalArrayBounds(ua,uLocal,I1,I2,I3,1);   
    if( !ok ) return 0;  // no points on this processor (NOTE: no communication should be done after this point)

  // -- we optimize for Cartesian grids (we can avoid creating the vertex array)
    const bool isRectangular=mg.isRectangular();
    if( !isRectangular )
        mg.update(MappedGrid::THEvertex | MappedGrid::THEcenter);
    OV_GET_SERIAL_ARRAY(real,mg.center(),xLocal);

    real dvx[3]={1.,1.,1.}, xab[2][3]={{0.,0.,0.},{0.,0.,0.}};
    int iv0[3]={0,0,0}; //
    int iv[3], &i1=iv[0], &i2=iv[1], &i3=iv[2];  // NOTE: iv[0]==i1, iv[1]==i2, iv[2]==i3
    real xv[3]={0.,0.,0.};
    if( isRectangular )
    {
        mg.getRectangularGridParameters( dvx, xab );
        for( int dir=0; dir<mg.numberOfDimensions(); dir++ )
        {
            iv0[dir]=mg.gridIndexRange(0,dir);
            if( mg.isAllCellCentered() )
                xab[0][dir]+=.5*dvx[dir];  // offset for cell centered
        }
    }
  // This macro defines the grid points for rectangular grids:
    #undef XC
    #define XC(iv,axis) (xab[0][axis]+dvx[axis]*(iv[axis]-iv0[axis]))

    
  // -- Store components here: 
    const int ex=0, ey=1, ez=2;
    const int hx=3, hy=4, hz=numberOfDimensions==2 ? 2 :  5;

    if( computeMagneticField && numberOfDimensions==3 && ua.getLength(3)<6 )
    {
        printF(" NonlinearExactSolution::ERROR: Not enough space in ua to hold the H field\n");
        OV_ABORT("error");
    }
    

  // --- Get Arrays for the dispersive model ----

    RealArray pLocal, qLocal;
    if( numberOfPolarizationVectors>0 )
    {
        OV_GET_SERIAL_ARRAY(real, pv,pLoc);
        pLocal.reference(pLoc);
    }
    if( numberOfAtomicLevels>0 )
    {
        OV_GET_SERIAL_ARRAY(real, qv,qLoc);
        qLocal.reference(qLoc);
    }
    
    if( caseName == "asymptoticSoliton" )
    {
    // ---- Asymptotic Solution ----

    // parameters for asymptotic solutions -- **FIX ME**
    // const real eps = 0.1;
        const real eps = sqrt(multilevelAtomicParams(0,0,0));
        const real U = 0.5;
        const real eta = 1.;
        const real x0 = 0.;
        const real epsHat = eps/2.*sqrt(eta/(U-U*U));
    // const real epsHat =.1;  // **** FIX ME ***
    // const real epsHat = sqrt(multilevelAtomicParams(0,0,0));
        real x,y,z=0.;
        if( t<1.5*dt )
        {
            printF("NES: asymptoticSoliton: U=%g (envelope speed), eta=%g, x0=%g, eps=%g, epsHat=%g\n",U,eta,x0,eps,epsHat);
        }
    // real epsilon=sqrt(multilevelAtomicParams(0,0,0));
    // printF("asymptoticSoliton: epsHat=%g, epsilon=%g\n ",epsHat,epsilon);
    // OV_ABORT("NonlinearExactSolution::finish me");

        if( numberOfTimeDerivatives==0 )
        {
            FOR_3D(i1,i2,i3,I1,I2,I3)
            {
                if( !isRectangular )
                {
                    x= xLocal(i1,i2,i3,0);  
                    y= xLocal(i1,i2,i3,1);
                    if( numberOfDimensions==3 ) z= xLocal(i1,i2,i3,2);
                }
                else
                {
                    x=XC(iv,0);
                    y=XC(iv,1);
                    if( numberOfDimensions==3 ) z=XC(iv,2);
                }

        // x0=0; U=0.5; eta=1.; 
        // epsHat=.1; % check -- sqrt(polarizationNECoefficients=[0.01])
                
        // E = @(x,t) ampE*sech( epsHat*(x-x0-U*t)).*sin(x-t);
        // P = @(x,t) 2*epsHat*tanh( epsHat*(x-x0-U*t) ).*sech( epsHat*(x-x0-U*t) ).*cos(x-t);
        // D = @(x,t) 1 - 2*sech( epsHat*(x-x0-U*t) ).^2;

                const real ampE = 2.*sqrt( eta*U/(1-U) );
                const real seche = 1./cosh( epsHat*(x-x0-U*t) );

                const real efield = 2.*sqrt(eta*U/(1.-U))/cosh(epsHat*(x-x0-U*t))*sin(x-t);
        // printF("efield %g\n",efield);
                const real et = - 2*sqrt(eta*U/(1. - U))*sin(-x + t)*epsHat*U*sinh(epsHat*(-U*t + x - x0))/pow(cosh(epsHat*(-U*t + x - x0)),2.) - 2*sqrt(eta*U/(1. - U))*cos(-x + t)/cosh(epsHat*(-U*t + x - x0));
        // printF("et %g\n",et);
                const real exxt = -12.*sqrt(eta*U/(1. - U))*sin(-x + t)*pow(epsHat,3)*pow(sinh(epsHat*(-U*t + x - x0)),3)*U/pow(cosh(epsHat*(-U*t + x - x0)),4) - 4.*sqrt(eta*U/(1. - U))*cos(-x + t)*pow(epsHat,2)*pow(sinh(epsHat*(-U*t + x - x0)),2)/pow(cosh(epsHat*(-U*t + x - x0)),3) + 10.*sqrt(eta*U/(1. - U))*sin(-x + t)*pow(epsHat,3)*sinh(epsHat*(-U*t + x - x0))*U/pow(cosh(epsHat*(-U*t + x - x0)),2) - 8.*sqrt(eta*U/(1. - U))*cos(-x + t)*pow(epsHat,2)*pow(sinh(epsHat*(-U*t + x - x0)),2)*U/pow(cosh(epsHat*(-U*t + x - x0)),3) + 4.*sqrt(eta*U/(1. - U))*sin(-x + t)*epsHat*sinh(epsHat*(-U*t + x - x0))/pow(cosh(epsHat*(-U*t + x - x0)),2) + 4.*sqrt(eta*U/(1. - U))*cos(-x + t)*pow(epsHat,2)*U/cosh(epsHat*(-U*t + x - x0)) + 2.*sqrt(eta*U/(1. - U))*cos(-x + t)*pow(epsHat,2)/cosh(epsHat*(-U*t + x - x0)) + 2.*sqrt(eta*U/(1. - U))*sin(-x + t)*epsHat*U*sinh(epsHat*(-U*t + x - x0))/pow(cosh(epsHat*(-U*t + x - x0)),2) + 2.*sqrt(eta*U/(1. - U))*cos(-x + t)/cosh(epsHat*(-U*t + x - x0));
        // printF("exxt %g\n",exxt);
                const real e1x = 2.*sqrt(eta*U/(1. - U))*sin(-x + t)*epsHat*sinh(epsHat*(-U*t + x - x0))/pow(cosh(epsHat*(-U*t + x - x0)),2) + 2.*sqrt(eta*U/(1. - U))*cos(-x + t)/cosh(epsHat*(-U*t + x - x0));
        // printF("e1x %g\n",ex1);
                const real exx = - 4.*sqrt(eta*U/(1. - U))*sin(-x + t)*pow(epsHat,2)*pow(sinh(epsHat*(-U*t + x - x0)),2)/pow(cosh(epsHat*(-U*t + x - x0)),3) - 4.*sqrt(eta*U/(1. - U))*cos(-x + t)*epsHat*sinh(epsHat*(-U*t + x - x0))/pow(cosh(epsHat*(-U*t + x - x0)),2) + 2.*sqrt(eta*U/(1. - U))*sin(-x + t)*pow(epsHat,2)/cosh(epsHat*(-U*t + x - x0)) + 2.*sqrt(eta*U/(1. - U))*sin(-x + t)/cosh(epsHat*(-U*t + x - x0));
        // printF("exx %g\n",exx);
                const real exxxx = - 48.*sqrt(eta*U/(1. - U))*sin(-x + t)*pow(eps,4)*pow(sinh(eps*(-U*t + x - x0)),4)/pow(cosh(eps*(-U*t + x - x0)),5) - 48.*sqrt(eta*U/(1. - U))*cos(-x + t)*pow(eps,3)*pow(sinh(eps*(-U*t + x - x0)),3)/pow(cosh(eps*(-U*t + x - x0)),4) + 56.*sqrt(eta*U/(1. - U))*sin(-x + t)*pow(eps,4)*pow(sinh(eps*(-U*t + x - x0)),2)/pow(cosh(eps*(-U*t + x - x0)),3) + 24.*sqrt(eta*U/(1. - U))*sin(-x + t)*pow(eps,2)*pow(sinh(eps*(-U*t + x - x0)),2)/pow(cosh(eps*(-U*t + x - x0)),3) + 40.*sqrt(eta*U/(1. - U))*cos(-x + t)*pow(eps,3)*sinh(eps*(-U*t + x - x0))/pow(cosh(eps*(-U*t + x - x0)),2) - 10.*sqrt(eta*U/(1. - U))*sin(-x + t)*pow(eps,4)/cosh(eps*(-U*t + x - x0)) + 8.*sqrt(eta*U/(1. - U))*cos(-x + t)*eps*sinh(eps*(-U*t + x - x0))/pow(cosh(eps*(-U*t + x - x0)),2) - 12.*sqrt(eta*U/(1. - U))*sin(-x + t)*pow(eps,2)/cosh(eps*(-U*t + x - x0)) - 2.*sqrt(eta*U/(1. - U))*sin(-x + t)/cosh(eps*(-U*t + x - x0));
        // printF("exxxx %g\n",exxxx);


                const real p = 2.*tanh(epsHat*(x-x0-U*t))/cosh(epsHat*(x-x0-U*t))*cos(x-t);
                const real pt = - 2*epsHat*U*(1 - pow(tanh(epsHat*(-U*t + x - x0)),2))*cos(-x + t)/cosh(epsHat*(-U*t + x - x0)) + 2*tanh(epsHat*(-U*t + x - x0))*cos(-x + t)*epsHat*U*sinh(epsHat*(-U*t + x - x0))/pow(cosh(epsHat*(-U*t + x - x0)),2) - 2*tanh(epsHat*(-U*t + x - x0))*sin(-x + t)/cosh(epsHat*(-U*t + x - x0));
                const real pxx = -4.*pow(epsHat,2)*tanh(epsHat*(-U*t + x - x0))*(1 - pow(tanh(epsHat*(-U*t + x - x0)),2))*cos(-x + t)/cosh(epsHat*(-U*t + x - x0)) - 4.*pow(epsHat,2)*(1 - pow(tanh(epsHat*(-U*t + x - x0)),2))*cos(-x + t)*sinh(epsHat*(-U*t + x - x0))/pow(cosh(epsHat*(-U*t + x - x0)),2) + 4.*epsHat*(1 - pow(tanh(epsHat*(-U*t + x - x0)),2))*sin(-x + t)/cosh(epsHat*(-U*t + x - x0)) + 4.*tanh(epsHat*(-U*t + x - x0))*cos(-x + t)*pow(epsHat,2)*pow(sinh(epsHat*(-U*t + x - x0)),2)/pow(cosh(epsHat*(-U*t + x - x0)),3) - 4.*tanh(epsHat*(-U*t + x - x0))*sin(-x + t)*epsHat*sinh(epsHat*(-U*t + x - x0))/pow(cosh(epsHat*(-U*t + x - x0)),2) - 2.*tanh(epsHat*(-U*t + x - x0))*cos(-x + t)*pow(epsHat,2)/cosh(epsHat*(-U*t + x - x0)) - 2.*tanh(epsHat*(-U*t + x - x0))*cos(-x + t)/cosh(epsHat*(-U*t + x - x0));
                const real ptxx = 2.*tanh(epsHat*(-U*t + x - x0))*sin(-x + t)*pow(epsHat,2)/cosh(epsHat*(-U*t + x - x0)) - 12.*pow(epsHat,3)*(1 - pow(tanh(epsHat*(-U*t + x - x0)),2))*cos(-x + t)*pow(sinh(epsHat*(-U*t + x - x0)),2)*U/pow(cosh(epsHat*(-U*t + x - x0)),3) + 8.*pow(epsHat,2)*tanh(epsHat*(-U*t + x - x0))*U*(1 - pow(tanh(epsHat*(-U*t + x - x0)),2))*sin(-x + t)/cosh(epsHat*(-U*t + x - x0)) + 8.*pow(epsHat,2)*(1 - pow(tanh(epsHat*(-U*t + x - x0)),2))*sin(-x + t)*U*sinh(epsHat*(-U*t + x - x0))/pow(cosh(epsHat*(-U*t + x - x0)),2) + 12.*tanh(epsHat*(-U*t + x - x0))*cos(-x + t)*pow(epsHat,3)*pow(sinh(epsHat*(-U*t + x - x0)),3)*U/pow(cosh(epsHat*(-U*t + x - x0)),4) - 10.*tanh(epsHat*(-U*t + x - x0))*cos(-x + t)*pow(epsHat,3)*sinh(epsHat*(-U*t + x - x0))*U/pow(cosh(epsHat*(-U*t + x - x0)),2) - 8.*tanh(epsHat*(-U*t + x - x0))*sin(-x + t)*pow(epsHat,2)*pow(sinh(epsHat*(-U*t + x - x0)),2)*U/pow(cosh(epsHat*(-U*t + x - x0)),3) - 2.*tanh(epsHat*(-U*t + x - x0))*cos(-x + t)*epsHat*U*sinh(epsHat*(-U*t + x - x0))/pow(cosh(epsHat*(-U*t + x - x0)),2) - 8.*pow(epsHat,3)*pow(tanh(epsHat*(-U*t + x - x0)),2)*U*(1 - pow(tanh(epsHat*(-U*t + x - x0)),2))*cos(-x + t)/cosh(epsHat*(-U*t + x - x0)) + 2.*tanh(epsHat*(-U*t + x - x0))*sin(-x + t)/cosh(epsHat*(-U*t + x - x0)) + 2.*epsHat*U*(1 - pow(tanh(epsHat*(-U*t + x - x0)),2))*cos(-x + t)/cosh(epsHat*(-U*t + x - x0)) + 4.*pow(epsHat,2)*tanh(epsHat*(-U*t + x - x0))*(1 - pow(tanh(epsHat*(-U*t + x - x0)),2))*sin(-x + t)/cosh(epsHat*(-U*t + x - x0)) + 4.*pow(epsHat,2)*(1 - pow(tanh(epsHat*(-U*t + x - x0)),2))*sin(-x + t)*sinh(epsHat*(-U*t + x - x0))/pow(cosh(epsHat*(-U*t + x - x0)),2) + 6.*pow(epsHat,3)*(1 - pow(tanh(epsHat*(-U*t + x - x0)),2))*cos(-x + t)*U/cosh(epsHat*(-U*t + x - x0)) - 4.*tanh(epsHat*(-U*t + x - x0))*sin(-x + t)*pow(epsHat,2)*pow(sinh(epsHat*(-U*t + x - x0)),2)/pow(cosh(epsHat*(-U*t + x - x0)),3) + 4.*tanh(epsHat*(-U*t + x - x0))*sin(-x + t)*pow(epsHat,2)*U/cosh(epsHat*(-U*t + x - x0)) - 4.*tanh(epsHat*(-U*t + x - x0))*cos(-x + t)*epsHat*sinh(epsHat*(-U*t + x - x0))/pow(cosh(epsHat*(-U*t + x - x0)),2) + 4.*pow(epsHat,3)*U*pow((1 - pow(tanh(epsHat*(-U*t + x - x0)),2)),2)*cos(-x + t)/cosh(epsHat*(-U*t + x - x0)) - 12.*pow(epsHat,3)*tanh(epsHat*(-U*t + x - x0))*(1 - pow(tanh(epsHat*(-U*t + x - x0)),2))*cos(-x + t)*U*sinh(epsHat*(-U*t + x - x0))/pow(cosh(epsHat*(-U*t + x - x0)),2) + 4.*epsHat*(1 - pow(tanh(epsHat*(-U*t + x - x0)),2))*cos(-x + t)/cosh(epsHat*(-U*t + x - x0));

                const real d = 1-2./pow((cosh(epsHat*(x-x0-U*t))),2);
                const real dx = 4.*epsHat*sinh(epsHat*(-U*t + x - x0))/pow(cosh(epsHat*(-U*t + x - x0)),3);
                const real dxx = -12.*pow(epsHat,2)*pow(sinh(epsHat*(-U*t + x - x0)),2)/pow(cosh(epsHat*(-U*t + x - x0)),4) + 4.*pow(epsHat,2)/pow(cosh(epsHat*(-U*t + x - x0)),2);


                const real pnec = eps;
                const real b1 = 0.;
                const real b0 = 1.;
                const real alphaP = eps;
                const real peptc = -eps;
                const real c = 1.;

                const real ptt = -b0*p-b1*pt+pnec*d*efield;
                const real ett = c*c*exx-alphaP*ptt;
                const real dt = peptc*efield*pt;
                const real pttxx = -b0*pxx-b1*ptxx+pnec*(dxx*efield+2.*dx*e1x+d*exx);

                const real pttt = -b0*pt-b1*ptt+pnec*dt*efield+pnec*d*et;
                const real ettt = c*c*exxt-alphaP*pttt;
                const real dtt = peptc*et*pt+peptc*efield*ptt;

                const real ptttt = -b0*ptt-b1*pttt+pnec*dtt*efield+2.*pnec*dt*et+pnec*d*ett;
                const real exxtt = c*c*exxxx-alphaP*pttxx;
                const real etttt = c*c*exxtt-alphaP*ptttt;

                if( numberOfDimensions==2 )
                {
                    uLocal(i1,i2,i3,ex) = 0.;
                    uLocal(i1,i2,i3,ey) = ampE*seche*sin(x-t);
                    if (t<0.){
                        uLocal(i1,i2,i3,ey) = ampE*seche*sin(x-0.) - dt*et + pow(dt,2)/2.*ett - pow(dt,3)/6.*ettt + pow(dt,4)/24.*etttt;
                    }
                    uLocal(i1,i2,i3,hz) = 0.;
                }
                else if( numberOfDimensions==3 )
                {
                    uLocal(i1,i2,i3,ex) = 0.;
                    uLocal(i1,i2,i3,ey) = ampE*seche*sin(x-t);
                    if (t<0.){
                        uLocal(i1,i2,i3,ey) = ampE*seche*sin(x-0.) - dt*et + pow(dt,2)/2.*ett - pow(dt,3)/6.*ettt + pow(dt,4)/24.*etttt;
                    }
                    uLocal(i1,i2,i3,ez) = 0.;
                }



        // --- assign polarization vectors ---
                for( int ip=0; ip<numberOfPolarizationVectors; ip++ )
                {
                    const int pc= ip*numberOfDimensions;
                    pLocal(i1,i2,i3,pc  ) = 0.;
                    pLocal(i1,i2,i3,pc+1) = 2*epsHat*tanh( epsHat*(x-x0-U*t) )*seche*cos(x-t);
                    if (t<0.){
                        pLocal(i1,i2,i3,pc+1) = 2*epsHat*tanh( epsHat*(x-x0-U*0.) )*seche*cos(x-0.) - epsHat*dt*pt + epsHat*pow(dt,2)/2.*ptt - epsHat*pow(dt,3)/6.*pttt + epsHat*pow(dt,4)/24.*ptttt;
                    }

                    if( numberOfDimensions==3 )
                    {
                        pLocal(i1,i2,i3,pc+2) = 0.;
                    }
                
                }
        // -- assign population densities ----
                for( int na=0; na<numberOfAtomicLevels; na++ )
                {
                    qLocal(i1,i2,i3,na) = 1. - 2.*SQR(seche);
                }


            }  // end for i1,i2,i3


        }
        else
        {
            OV_ABORT("NonlinearExactSolution::ERROR: numberOfTimeDerivatives != 0 ");
            
        } // end if number of time derivatives 
    }
    else if( caseName == "soliton")
    {
    //  =============== Soliton-like solution defined by the numerical solution of a 1D problem ====
    


        const IntegerArray &gid=mg.gridIndexRange();
        int n=gid(1,0)-gid(0,0);

    // ::display(gid,"gid");
    // get 
    // realArray un;
    // const int size = n+3; // index goes from -1 to n+1
    // const int base = -1;
    // const int bound = n+1;

        const int size = n+5; // index goes from -2 to n+2
        const int base = -2;
        const int bound = n+2;

        real *pun = new real[3*size]; // hardcoded 3 variables [E,P,D]; see mbe1d.f90 for details.
        #define un(i,j) pun[(i-base)+(size)*(j)] //

    // real dvx[3]={1.,1.,1.}, xab[2][3]={{0.,0.,0.},{0.,0.,0.}};
    // mg.getRectangularGridParameters( dvx, xab );

        real xa = xab[0][0];
        real xb = xab[1][0];

        real epsilon=sqrt(multilevelAtomicParams(0,0,0));

        real rpar[5];
        rpar[0] = t;
        rpar[1] = dt;
        rpar[2] = xa;
        rpar[3] = xb;
        rpar[4] = epsilon;


        int ipar[5];
    // ipar[0] = n;
        ipar[1] = orderOfAccuracyInSpace;
      
    // Call a Fortran subroutine
    // mbe1d(epsilon,dt,xa,xb,n,t,un(base,0));
        mbe1d(rpar[0],ipar[0],n,un(base,0));
    // ::display(multilevelAtomicParams,"params");
    // cout << "At time (NES): " << t << endl;
    // cout << "Time step (NES): " << dt << endl;
    // cout << "Bounds in x are (NES): " << I1.getBase() << " "<<  I1.getBound() << " " << n <<" " << gid(1,0) << " " << gid(0,0) << endl;
    // cout << "Bounds in (x,y,z) are (NES): " << I1a.getBase() << " "<<  I1a.getBound() << " " << I2a.getBase() << " "<<  I2a.getBound() << " " << I3a.getBase() << " "<<  I3a.getBound() << endl;


        realArray v(I1,3);
        v = 0.;
        for (int i = I1.getBase(); i<=I1.getBound(); i++){
            for (int j=0;j<3;j++){
                    v(i,j) = un(i,j);
            }
        }

    // ::display(v,"v");

        real x0[3]={0.,0.,0.};   //     
        real x,y,z=0.;

    // ::display(pLocal,"pLocal");
        if (numberOfPolarizationVectors*numberOfDimensions > pLocal.getLength(3)){
            printF("Error: too many polarization vectors");
            OV_ABORT("Error");
        }

        if( numberOfTimeDerivatives==0 )
        {
            FOR_3D(i1,i2,i3,I1,I2,I3)
            {
                if( !isRectangular )
                {
                    x= xLocal(i1,i2,i3,0)-x0[0];   // shift point to reference coordinates 
                    y= xLocal(i1,i2,i3,1)-x0[1];
                    if( numberOfDimensions==3 ) z= xLocal(i1,i2,i3,2)-x0[2];
                }
                else
                {
                    x=XC(iv,0)-x0[0];
                    y=XC(iv,1)-x0[1];
                    if( numberOfDimensions==3 ) z=XC(iv,2)-x0[2];
                }

      
                if( numberOfDimensions==2 )
                {
                    uLocal(i1,i2,i3,ex) = 0.;
                    uLocal(i1,i2,i3,ey) = v(i1,0);
                    uLocal(i1,i2,i3,hz) = 0.;
                }
                else if( numberOfDimensions==3 )
                {
                    uLocal(i1,i2,i3,ex) = 0.;
                    uLocal(i1,i2,i3,ey) = v(i1,0);
                    uLocal(i1,i2,i3,ez) = 0.;
                }



        // --- assign polarization vectors ---
                for( int ip=0; ip<numberOfPolarizationVectors; ip++ )
                {
                    const int pc= ip*numberOfDimensions;
                    pLocal(i1,i2,i3,pc  ) = 0.;
                    pLocal(i1,i2,i3,pc+1) = v(i1,1);

                    if( numberOfDimensions==3 )
                    {
                        pLocal(i1,i2,i3,pc+2) = 0.;
                    }
                
                }
        // -- assign population densities ----
                for( int na=0; na<numberOfAtomicLevels; na++ )
                {
                    qLocal(i1,i2,i3,na) = v(i1,2);
                }


            }  // end for i1,i2,i3


        }
        else
        {
            OV_ABORT("NonlinearExactSolution::ERROR: numberOfTimeDerivatives != 0 ");
            
        } // end if number of time derivatives 
            
    }
    else
    { 
        OV_ABORT("ERROR: unknown caseName");

    }
    return 0;
}







// ===============================================================================
/// \brief Check the solution.
// ===============================================================================
int NonlinearExactSolution::
NonlinearExactSolution::check()
{

  // ------------- CHECK THAT THE EQUATIONS ARE SATISFIED AT POINTS INSIDE AND OUTSIDE ------

    printF("------------ NonlinearExactSolution::check: CHECK THE EQUATIONS ------------\n\n");

    real maxErr=0.;

    real errE = 0.;
    real errP = 0.;
    real errD = 0.;

  // E_tt - E_xx = epsilon*alphaP*P_tt
  // P_tt + P = epsilon*(N0-N1)*E
  // D_t = -epsilon*E*P_t (D=N0-N1)



    printF("\n Max-error in tests = %9.2e. TESTS %s.\n",maxErr,(maxErr<1.e-5 ? "PASSED" : "***FAILED***"));

    printF("\n------------ FINSHED CHECK EQUATIONS ------------\n\n");


    
    return 0;
}


