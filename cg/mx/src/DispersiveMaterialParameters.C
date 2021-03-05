// Generalized Dispersion Model:
//       E_tt - c^2 Delta(E) = -alphaP P_tt
//       P_tt + b1 P_1 + b0 = a0*E + a1*E_t 
// 
#include "DispersiveMaterialParameters.h"
#include "PlotStuff.h"
#include "display.h"
#include "PlaneInterfaceExactSolution.h"

int 
getLineFromFile( FILE *file, char s[], int lim);

// define evalDispersionRelation EXTERN_C_NAME(evaldispersionrelation)
// define evalGeneralizedDispersionRelation EXTERN_C_NAME(evalgeneralizeddispersionrelation)
#define evalEigGDM EXTERN_C_NAME(evaleiggdm)
#define evalInverseGDM EXTERN_C_NAME(evalinversegdm)
#define evalComplexIndexOfRefraction EXTERN_C_NAME(evalcomplexindexofrefraction)


extern "C"
{
  // void evalDispersionRelation( const real& cc, const real& eps, const real& gam, const real& omegap, const real& k, 
  //                              real& reS, real& imS);

  // void evalGeneralizedDispersionRelation( const real& c, const real& k, const real& a0, const real& a1, 
  //                                         const real& b0, const real& b1,const real& alphaP, 
  //                                         real& reS, real& imS, real & psir, real & psi  );

  // compute GGM eigenvalues for multiple polarization vectors 
  void evalEigGDM( const int & mode, const int & Np, const real& c, const real& k, const real& a0, const real& a1, 
                   const real& b0, const real& b1, 
                   real& reS, real& imS, real & srm, real & sim, real & chir, real & chi, real & chiSumr, real & chiSumi  );


  // Evaluate the "INVERSE" dispersion relation (compute k=*(kr,ki) given s=(sr,si) 
  // for the generalized dispersion model (GDM) With multiple polarization vectors 
  void evalInverseGDM( const real&c, const real&sr,const real&si, const int&Np,
                       const real&a0,const real&a1,const real&b0,const real&b1,
                       real&kr,real&ki,real&chir,real&chii, real & chiSumr, real & chiSumi );

  void evalComplexIndexOfRefraction(const real&mu1,const real&eps1,const real&chi1r,const real&chi1i, 
                                    const real&mu2,const real&eps2,const real&chi2r,const real&chi2i,
                                    real & mr, real & mi);


}

// lapack routines
#ifdef OV_USE_DOUBLE
  #define GETRF EXTERN_C_NAME(dgetrf)
  #define GETRI EXTERN_C_NAME(dgetri)
  #define GECON EXTERN_C_NAME(dgecon)
  #define LANGE EXTERN_C_NAME(dlange)
  #define GEEV  EXTERN_C_NAME(dgeev)
#else
  #define GETRF EXTERN_C_NAME(sgetrf)
  #define GETRI EXTERN_C_NAME(sgetri)
  #define GECON EXTERN_C_NAME(sgecon)
  #define LANGE EXTERN_C_NAME(slange)
  #define GEEV  EXTERN_C_NAME(sgeev)
#endif

extern "C"
{
  void GETRF( int & m, int & n, real & a, const int & lda, int & ipvt, int & info );
  void GETRI( int & n, real & a, const int & lda, const int & ipvt, real & work, const int & iwork, int & info );

  void GECON( char *norm, int & n, real & a, const int & lda, real & anorm, real & rcond, real & work, int & iwork, int & info );
  real LANGE( char *norm, int & m, int & n, real & a, const int & lda, real & work );

   void GEEV( char *jobvl, char* jobvr, int & n, real & a, const int & lda,
              real & wr, real & wi, real &vl, int & ldvl, real & vr, int & ldvr, real & work, int & lwork, int & info );

}


// ============================================================================
/// \brief Class to define parameters of a dispersive material.
///
/// Generalized Dispersion Model:
///       E_tt - c^2 Delta(E) = -(1/eps) P_tt
///       P_tt + b1 P_1 + b0 = eps*( a0*E + a1*E_t )
///
///  Bianisotropic Material tensor: (*new* May 5, 2019) 
///        [ D ] = K  [ E ]
///        [ B ]      [ H ]
/// 
///    K(6,6) : material tensor 
///
///     K0(6,6)  : constant part of material tensor 
///     bianisotropicParameters(4,Np,6,6)        : GDM 
///     Np(6,6) : number of polarization vectors 
// ============================================================================
DispersiveMaterialParameters::
DispersiveMaterialParameters()
{

  // general dispersive model parameters:
  //   modelParameters(i,k)  : i=0,1,2,3,4 are the parmeters in the equation 
  //                           for P_k , k=1,2,...,numberOfPolarizationVectors
  // Polarization equation for vector P_k
  //   (P_k)_tt + b1_k (P_k)_t + b0_k P_k = a0_k E + a1_k E_t
  // modelParameters(0:3,k) = [a0,a1,b0,b1] 

  alphaP=-1.;  // this means set to default, alphaP=1/eps
  epsInf=1.;    // eps(s=infinity)
  muInf =1.;

  dbase.put<bool>("isDispersive")= false;  // by default a domain is non-dispersive

  numberOfPolarizationVectors=0; // by default a domain is non-dispersive
  numberOfModelParameters=4;     // [a0,a1,b0,b1] 
  modelParameters.redim(numberOfModelParameters,1); // fill in defaults of zero
  modelParameters=0.;

  nonlinearModel = noNonlinearModel;
  dbase.put<aString>("nonlinearModelName")="none"; 
  dbase.put<bool>("isNonlinearMaterial")=false;
  dbase.put<RealArray>("nonlinearModelParameters");  // holds parameters for nonlinear models

  dbase.put<int>("numberOfAtomicLevels")=0;  // number of atomic levels ( N's ) for the multilevelAtomic (Maxwell-Bloch) model

  // **OLD WAY: 
  // Drude-Lorentz model:  
  //     P_tt + gamma P_t = omegap^2 E 
  gamma=1.;  // damping 
  omegap=1.; // plasma frequency

  // We save the root after it has been computed to save computations
  mode=-1;
  rootComputed=false;
  ck0=0; sr0=0; si0=0; 
  chiSumr0=0.; chiSumi0=0.;
  
  dbase.put<int>("debug")=-1;  // -1 : use default 

  // new way to save parameters
  dbase.put<aString>("name")="generic";

  dbase.put<MaterialTypeEnum>("materialType")=isotropic;

  const real nm = 1e-9;         // nanometers  (meter-per-nm)
  const real um = 1e-6;         // micrometers
  const real c0 = 299792458;    // the speed of light, [m/c]
  real L0 = 100*nm;             // length scale

  dbase.put<real>("velocityScale")= c0;  // velocity scale 
  dbase.put<real>("lengthScale")  = L0;  // length scale 

  // Scale for omega: 
  dbase.put<real>("omegaScale")=1.;

  // Range in omega for which the fit was made
  dbase.put<real>("omegaMin")=.2;
  dbase.put<real>("omegaMax")=1.;

  dbase.put<DispersionRelationOptionEnum>("dispersionRelationComputeOption")=computeComplexFrequency;

  dbase.put<bool>("normalizedUnits")=true;      // if true, do not scale GDM parameters by omegaScale


}

// ============================================================================
/// \brief Copy constructor
// ============================================================================
DispersiveMaterialParameters::
DispersiveMaterialParameters(const DispersiveMaterialParameters& x)
{
  *this=x;  

  // why do we have all these following that should be handled by operator= ??? ********
  if( !dbase.has_key("name" ) ) dbase.put<aString>("name");
  dbase.get<aString>("name") = x.dbase.get<aString>("name");

  if( !dbase.has_key("debug") ) dbase.put<int>("debug");
  dbase.get<int>("debug")= x.dbase.get<int>("debug");

  if( !dbase.has_key("materialType") ) dbase.put<MaterialTypeEnum>("materialType");
  dbase.get<MaterialTypeEnum>("materialType")= x.dbase.get<MaterialTypeEnum>("materialType");

  if( !dbase.has_key("velocityScale") ) dbase.put<real>("velocityScale");
  dbase.get<real>("velocityScale")  = x.dbase.get<real>("velocityScale");

  if( !dbase.has_key("lengthScale") ) dbase.put<real>("lengthScale");
  dbase.get<real>("lengthScale")  = x.dbase.get<real>("lengthScale");

  if( !dbase.has_key("omegaScale") ) dbase.put<real>("omegaScale");
  dbase.get<real>("omegaScale")  = x.dbase.get<real>("omegaScale");

  if( !dbase.has_key("omegaMin") ) dbase.put<real>("omegaMin");
  dbase.get<real>("omegaMin")  = x.dbase.get<real>("omegaMin");

  if( !dbase.has_key("omegaMax") ) dbase.put<real>("omegaMax");
  dbase.get<real>("omegaMax")  = x.dbase.get<real>("omegaMax");

  if( x.dbase.has_key("NpBA") )
  {
    dbase.put<IntegerArray>("NpBA") = x.dbase.get<IntegerArray>("NpBA");
    dbase.put<RealArray>("K0")      = x.dbase.get<RealArray>("K0");
    dbase.put<RealArray>("bianisotropicParameters") = x.dbase.get<RealArray>("bianisotropicParameters");
  }
  
  if( !dbase.has_key("normalizedUnits") ) dbase.put<bool>("normalizedUnits");
  dbase.get<bool>("normalizedUnits")  = x.dbase.get<bool>("normalizedUnits");


}

// ============================================================================
/// \brief Destructor.
// ============================================================================
DispersiveMaterialParameters::
~DispersiveMaterialParameters()
{
}


DispersiveMaterialParameters & DispersiveMaterialParameters::
operator =( const DispersiveMaterialParameters & x )
{
  numberOfPolarizationVectors=x.numberOfPolarizationVectors;
  numberOfModelParameters=x.numberOfModelParameters;
  modelParameters.redim(0);
  modelParameters=x.modelParameters;
  chir0.redim(0);
  chir0 = x.chir0;
  chii0.redim(0);
  chii0 = x.chii0;
  chiSumr0 = x.chiSumr0;
  chiSumi0 = x.chiSumi0;

  epsInf=x.epsInf;
  muInf =x.muInf;
  
  alphaP=x.alphaP;
  mode  =x.mode;
  rootComputed=x.rootComputed;
  ck0=x.ck0;
  
  gamma=x.gamma;
  omegap=x.omegap;

  if( !dbase.has_key("isDispersive" ) ) dbase.put<bool>("isDispersive");
  dbase.get<bool>("isDispersive")= x.dbase.get<bool>("isDispersive");
  

  nonlinearModel = x.nonlinearModel;
  if( !dbase.has_key("nonlinearModelName") ) dbase.put<aString>("nonlinearModelName");
  dbase.get<aString>("nonlinearModelName")= x.dbase.get<aString>("nonlinearModelName");

  if( !dbase.has_key("isNonlinearMaterial") ) dbase.put<bool>("isNonlinearMaterial");
  dbase.get<bool>("isNonlinearMaterial")= x.dbase.get<bool>("isNonlinearMaterial");

  if( !dbase.has_key("numberOfAtomicLevels") ) dbase.put<int>("numberOfAtomicLevels");
  dbase.get<int>("numberOfAtomicLevels")= x.dbase.get<int>("numberOfAtomicLevels");

  if( !dbase.has_key("nonlinearModelParameters") ) dbase.put<RealArray>("nonlinearModelParameters");
  RealArray & nonlinearModelParameters = dbase.get<RealArray>("nonlinearModelParameters");
  nonlinearModelParameters.redim(0);
  nonlinearModelParameters = x.dbase.get<RealArray>("nonlinearModelParameters"); 

  if( !dbase.has_key("debug" ) ) dbase.put<int>("debug");
  dbase.get<int>("debug")= x.dbase.get<int>("debug");

  if( !dbase.has_key("name" ) ) dbase.put<aString>("name");
  dbase.get<aString>("name")= x.dbase.get<aString>("name");
  
  if( !dbase.has_key("materialType") ) dbase.put<MaterialTypeEnum>("materialType");
  dbase.get<MaterialTypeEnum>("materialType")= x.dbase.get<MaterialTypeEnum>("materialType");

  if( !dbase.has_key("velocityScale") ) dbase.put<real>("velocityScale");
  dbase.get<real>("velocityScale")  = x.dbase.get<real>("velocityScale");

  if( !dbase.has_key("lengthScale") ) dbase.put<real>("lengthScale");
  dbase.get<real>("lengthScale")    = x.dbase.get<real>("lengthScale");

  if( !dbase.has_key("omegaScale") ) dbase.put<real>("omegaScale");
  dbase.get<real>("omegaScale")    = x.dbase.get<real>("omegaScale");

  if( !dbase.has_key("omegaMin") ) dbase.put<real>("omegaMin");
  dbase.get<real>("omegaMin")    = x.dbase.get<real>("omegaMin");

  if( !dbase.has_key("omegaMax") ) dbase.put<real>("omegaMax");
  dbase.get<real>("omegaMax")    = x.dbase.get<real>("omegaMax");

  if( !dbase.has_key("dispersionRelationComputeOption") )
    dbase.put<DispersionRelationOptionEnum>("dispersionRelationComputeOption");
  dbase.get<DispersionRelationOptionEnum>("dispersionRelationComputeOption") =
             x.dbase.get<DispersionRelationOptionEnum>("dispersionRelationComputeOption");

  if( !dbase.has_key("normalizedUnits") ) dbase.put<bool>("normalizedUnits");
  dbase.get<bool>("normalizedUnits")  = x.dbase.get<bool>("normalizedUnits");


  return *this;
}

// ==========================================================================================
/// \brief Set the dispersionRelationComputeOption
// ==========================================================================================
int DispersiveMaterialParameters::setDispersionRelationComputeOption( DispersionRelationOptionEnum  computeOption )
{
  dbase.get<DispersionRelationOptionEnum>("dispersionRelationComputeOption")=computeOption;
  return 0;
}


  
// ==========================================================================================
/// \brief Return true if the material is dispersive.
// ==========================================================================================
bool DispersiveMaterialParameters::isDispersiveMaterial() const
{
  return dbase.get<bool>("isDispersive");
}


// ==========================================================================================
/// \brief Return true if the material is nonlinear
// ==========================================================================================
bool DispersiveMaterialParameters::isNonlinearMaterial() const
{
  return dbase.get<bool>("isNonlinearMaterial");
}

// ==========================================================================================
/// \brief Return the number of atomic levels (N's) in the multilevelAtomic model
// ==========================================================================================
int DispersiveMaterialParameters::getNumberOfAtomicLevels() const
{
  return dbase.get<int>("numberOfAtomicLevels");
}

// ==========================================================================================
/// \brief Return the number of polarization vectors for isotropic materials
// ==========================================================================================
int DispersiveMaterialParameters::getNumberOfPolarizationVectors() const
{
  return numberOfPolarizationVectors;
}



// ==========================================================================================
/// \brief Return the material name 
// ==========================================================================================
aString DispersiveMaterialParameters::getMaterialName() const
{
  return dbase.get<aString>("name");
}

// ==========================================================================================
/// \brief Set the debug (bit) flag
// ==========================================================================================
int DispersiveMaterialParameters::setDebug( int debugFlag )
{
  dbase.get<int>("debug")= debugFlag;
}


// ==========================================================================================
/// \brief Display parameters
// ==========================================================================================
int DispersiveMaterialParameters::
display( FILE *file /* = stdout */, const aString & label /* = nullString */   ) const
{
  const aString & name = dbase.get<aString>("name");
  const bool & isDispersive = dbase.get<bool>("isDispersive");
  const DispersionRelationOptionEnum & dispersionRelationComputeOption =
                dbase.get<DispersionRelationOptionEnum>("dispersionRelationComputeOption");
  
  MaterialTypeEnum & materialType = dbase.get<MaterialTypeEnum>("materialType");
  const real & omegaMin = dbase.get<real>("omegaMin");   
  const real & omegaMax = dbase.get<real>("omegaMax");   
  
  const real V0 = dbase.get<real>("velocityScale");
  const real L0 = dbase.get<real>("lengthScale"); 
  const real & omegaScale = dbase.get<real>("omegaScale");

  const real T0 = L0/V0;
  const real Omega0 = 1./T0;
  const real wScale = omegaScale/Omega0;
  
  // normalizedUnits : if true, do not scale GDM parameters by wScale = omegaScale/Omega0
  const bool & normalizedUnits = dbase.get<bool>("normalizedUnits");  


  fPrintF(file,"Label=%s\n",(const char*)label);
  fPrintF(file,"----------------------------------------------------------------------------\n");
  fPrintF(file,"-------------------- Dispersive Material Parameters ------------------------\n");
  
  fPrintF(file," name=[%s]\n",(const char*)name);
  fPrintF(file," materialType=%s\n",(materialType==isotropic ? "isotropic" :
                               materialType==bianisotropic ? "bianisotropic" : "unknown"));
  fPrintF(file," isDispersive=%s.\n",(isDispersive ? "yes" : "no"));
  fPrintF(file," dispersionRelationComputeOption=%s.\n",
	     (dispersionRelationComputeOption==computeComplexFrequency ? "computeComplexFrequency" : "computeComplexWaveNumber"));
  fPrintF(file," Length-scale L0=%9.3e, velocity-scale V0=%9.3e, Omega0=V0/L0=%9.3e, omegaScale=%9.3e \n",
                L0,V0,Omega0,omegaScale);
  if( isDispersive )
  {
    if( normalizedUnits )
      fPrintF(file," normalizedUnits=true : GDM parameters are NOT scaled by omegaScale/Omega0=%g\n",wScale);
    else
      fPrintF(file," normalizedUnits=false : GDM parameters ARE scaled by omegaScale/Omega0=%g.\n",wScale);
  }
  
  if( materialType==isotropic )
  {
    fPrintF(file," Generalized Dispersion Material Model (GDM) parameters:\n");
    fPrintF(file," numberOfPolarizationVectors=%i\n",numberOfPolarizationVectors);
    fPrintF(file," epsInf=%10.4e, alphaP*epsInf=%g\n",epsInf,alphaP*epsInf);
    fPrintF(file,"    Susceptibilities:  \n"
                 "  Chi_j(s) = ( a_{0,j} + a_{1,j} s )/( b_{0,j} + b_{1,j}s + s^2 ) , j=0,1,2,...,Np-1\n");
    
    for( int j=0; j<numberOfPolarizationVectors; j++ )
    {
      fPrintF(file,"  j=%d a0=%10.4e a1=%10.4e b0=%10.4e b1=%10.4e\n",
             j,modelParameters(0,j),modelParameters(1,j),modelParameters(2,j),modelParameters(3,j));
    }
  }
  else if( materialType==bianisotropic )
  {
    RealArray & K0 = dbase.get<RealArray>("K0");
    // ::display(K0,"K0",file);
    for( int k1=0; k1<6; k1++ )
    {
      if( k1==2 )
        fPrintF(file,"K0 = [");
      else
        fPrintF(file,"     [");
      for( int k2=0; k2<6; k2++ )
      {
        fPrintF(file,"%9.3e ",K0(k1,k2));
      }
      fPrintF(file,"]\n");
    }
    
    // bianisotropicParameters(4,Np,6,6) 
    RealArray & bianisotropicParameters = dbase.get<RealArray>("bianisotropicParameters");
    IntegerArray & NpBA = dbase.get<IntegerArray>("NpBA");

    for( int k1=0; k1<6; k1++ )
    {
      for( int k2=0; k2<6; k2++ )
      {
	if( NpBA(k1,k2)>0 )
	{
	  fPrintF(file,"K(%d,%d): Np=%d:\n",k1,k2,NpBA(k1,k2));
	  for( int n=0; n<NpBA(k1,k2); n++ )
	  {
	    fPrintF(file,"   n=%d: [a0,a1,b0,b1]=[%9.3e,%9.3e,%9.3e,%9.3e]\n",n,
		    bianisotropicParameters(0,n,k1,k2),
		    bianisotropicParameters(1,n,k1,k2),
		    bianisotropicParameters(2,n,k1,k2),
		    bianisotropicParameters(3,n,k1,k2));
	  }
	}
      }
    }
    


  }
  
  fPrintF(file,"----------------------------------------------------------------------------\n");
  



  return 0;
}



// ==========================================================================================
/// \brief Compute the real and imaginary parts of the disperion relation parameter "s"
///
/// \param c,k (input) 
/// \param  sr,si (ouptut) : real and imaginary parts of s omega in exp(s*t)*exp(i k*x )
/// \param  chir,chii  (ouptut) : real and imaginary parts of the electric susceptibility chi, P=eps*chi*E
/// \param chiSumi,pshiSumr (output) real and imaginary parts of chi 
// ==========================================================================================
int DispersiveMaterialParameters::
evaluateDispersionRelation( const real c, const real k, real & sr, real & si, real chir[], real chii[], 
                            real & chiSumr, real & chiSumi  )
{
  // assert( numberOfPolarizationVectors==1 );
  assert( numberOfModelParameters==4 );
  
  const real a0=modelParameters(0,0);
  const real a1=modelParameters(1,0);
  const real b0=modelParameters(2,0);
  const real b1=modelParameters(3,0);
  
  // We save the root after it has been computed to save computations

  if( !rootComputed || fabs(c*k-ck0) > REAL_EPSILON*10*abs(ck0+1.) )
  { 
    printF("--DMP-- GDM: RECOMPUTE root: rootComputed=%i c*k=%e ck0=%e (Np=%i, mode=%i)\n",
           (int)rootComputed,c*k,ck0,numberOfPolarizationVectors,mode);
    rootComputed=true;
    ck0=c*k;  

    if( chir0.getLength(0)!=numberOfPolarizationVectors )
    {
      chir0.redim(numberOfPolarizationVectors); chir0=0.;
      chii0.redim(numberOfPolarizationVectors); chii0=0.;
    }

    int Np=numberOfPolarizationVectors;
    RealArray a0v(Np), a1v(Np), b0v(Np), b1v(Np);
    for( int j=0; j<Np; j++ )
    {
      a0v(j)=modelParameters(0,j); 
      a1v(j)=modelParameters(1,j);
      b0v(j)=modelParameters(2,j);
      b1v(j)=modelParameters(3,j);
    }
      
    int neig = 2*Np+2; // total number of eigenvalues "s"
    RealArray srv(neig), siv(neig);
    // (sr,si) = eigenvalue with largest imaginary part
    evalEigGDM( mode, Np, c, k, a0v(0),a1v(0),b0v(0),b1v(0), srv(0), siv(0), sr0,si0, chir0(0),chii0(0),chiSumr0,chiSumi0 );
    if( false )
    {
      for( int i=0; i<neig; i++ )
      {
        printF("--DMP-- GDM: c=%g k=%g i=%d: real(s)=%g, Im(s)=%g *NEW*\n",
               c,k, i, srv(i),siv(i));
      }
    }

    bool printResults=true;
    if( printResults )
      printF("--DMP-- GDM: s=(%20.12e,%20.12e) :\n",sr0,si0);
      
    for( int j=0; j<numberOfPolarizationVectors; j++ )
    {
      if( printResults )
        printF("  j=%d a0=%9.3e a1=%9.3e b0=%9.3e b1=%9.3e chir=%20.12e chii=%20.12e\n",
               j,a0v(j),a1v(j), b0v(j), b1v(j), chir0(j),chii0(j));
    }
    

  }

  sr=sr0; si=si0;
  chiSumr = chiSumr0;
  chiSumi = chiSumi0;
  for( int j=0; j<numberOfPolarizationVectors; j++ )
  {
    chir[j]=chir0(j); 
    chii[j]=chii0(j);  // 
  }

  return 0;
}

// ==========================================================================================
/// \brief Evaluate the "inverse" dispersion relation and 
///       return the complex wave number k=(kr,ki) given s=(sr,si). Also return chi(j) 
///
/// \param (sr,si) (input) : "s"
/// \param kr,ki (output) :
/// \param  chir,chii  (ouptut) : real and imaginary parts of the electric susceptibility chi, P=eps*chi*E
/// \param chiSumi,pshiSumr (output) real and imaginary parts of chi 
// ==========================================================================================
int DispersiveMaterialParameters::
evaluateComplexWaveNumber( const real c, const real & sr, const real & si, 
                           real & kr, real &ki, real chir[], real chii[], real & chiSumr, real & chiSumi  )
{

  if( chir0.getLength(0)!=numberOfPolarizationVectors )
  {
    chir0.redim(numberOfPolarizationVectors); chir0=0.;
    chii0.redim(numberOfPolarizationVectors); chii0=0.;
  }

  const int Np=numberOfPolarizationVectors;
  RealArray a0v(Np), a1v(Np), b0v(Np), b1v(Np);
  for( int j=0; j<Np; j++ )
  {
    a0v(j)=modelParameters(0,j); 
    a1v(j)=modelParameters(1,j);
    b0v(j)=modelParameters(2,j);
    b1v(j)=modelParameters(3,j);
  }

  evalInverseGDM( c, sr,si, numberOfPolarizationVectors,a0v(0),a1v(0),b0v(0),b1v(0), 
                  kr,ki,chir0(0),chii0(0),chiSumr0,chiSumi0 );

  if( false )
    printF("--DMP-- evaluateComplexWaveNumber: numberOfPolarizationVectors=%i s=(%9.3e,%9.3e) --> k=(%9.3e,%9.3e)\n",numberOfPolarizationVectors,sr,si,kr,ki);

  
  chiSumr = chiSumr0;
  chiSumi = chiSumi0;
  for( int j=0; j<numberOfPolarizationVectors; j++ )
  {
    chir[j]=chir0(j);
    chii[j]=chii0(j);
        
    if( false )
      printF("evaluateComplexWaveNumber:  j=%d a0=%9.3e a1=%9.3e b0=%9.3e b1=%9.3e chir=%20.12e chii=%20.12e\n",
             j,a0v(j),a1v(j), b0v(j), b1v(j), chir0(j),chii0(j));
  }

  return 0;
}



// // ==========================================================================================
// /// \brief Compute the real and imaginary parts of the disperion relation parameter "s"
// ///
// /// \param c,eps,mu,k (input) 
// /// \param  reS,imS (ouptut) : real and imaginary parts of s omega in exp(s*t)*exp(i k*x )
// // ==========================================================================================
// int DispersiveMaterialParameters::
// computeDispersionRelation( const real c, const real eps, const real mu, const real k, 
//                            real & reS, real & imS )
// {
//   // ****** OLD WAY ********

//   evalDispersionRelation( c, eps, gamma, omegap, k,  reS, imS );
  
//   printF("--DispersiveMaterialParameters-- dispersion-relation: c=%g eps=%g mu=%g gamma=%g omegap=%g"
//          " -> real(s)=%g, Im(s)=%g\n",
// 	 c,eps,mu,gamma,omegap, reS,imS );

//   return 0;
// }


// // ==========================================================================================
// /// \brief Compute the real and imaginary parts of the dispersive plane wave "omega"
// ///
// /// \param c,eps,mu,k (input) 
// /// \param  omegar,omegai (ouptut) : real and imaginary parts of omega in exp(i(k*x-omega*t))
// // ==========================================================================================
// int DispersiveMaterialParameters::
// computeDispersivePlaneWaveParameters( const real c, const real eps, const real mu, const real k, 
//                                       real & omegar, real & omegai )
// {
//   // ****** OLD WAY ********

//   real reS, imS;
//   evalDispersionRelation( c, eps, gamma, omegap, k,  reS, imS );
//   omegar=imS;
//   omegai=reS;
  
//   printF("--DispersiveMaterialParameters-- dispersion-relation: c=%g eps=%g mu=%g gamma=%g omegap=%g -> omegar=%g, omegai=%g\n",
// 	 c,eps,mu,gamma,omegap, omegar,omegai);

//   return 0;
// }


// ==========================================================================================
/// \brief  Evaluate the complex index of refraction m=(mr,mi) for scattering off dispersive dielectrics
/// \params mu1,eps1 : region 1
/// \params chi1r,chi1i : chi=chi1r + I*chi1i , real and imaginarty parts of the electric susceptibility
///
///    m = mr + i*mi = sqrt{ mu_2 eps_2 (1+ Chi_2(s) ) / mu_1 eps_1 (1+Chi_1(s) ) }
///
// ==========================================================================================
int DispersiveMaterialParameters::
evaluateComplexIndexOfRefraction( const real mu1, const real eps1, const real chi1r, const real chi1i, 
                                  const real mu2, const real eps2, const real chi2r, const real chi2i, 
                                  real  & mr, real & mi ) const
{

  evalComplexIndexOfRefraction(mu1,eps1,chi1r,chi1i, mu2,eps2,chi2r,chi2i, mr,mi);
  
  return 0;
}

// ==========================================================================================
/// \brief return the parameter epsInf (large s limit of epsHat)
// ==========================================================================================
real DispersiveMaterialParameters::
getEpsInf() const
{
  return epsInf;
}

// ==========================================================================================
/// \brief return the parameter muInf (large s limit of muHat)
// ==========================================================================================
real DispersiveMaterialParameters::
getMuInf() const
{
  return muInf;
}


// ==========================================================================================
/// \brief return the parameter alphaP -- normally = 1/eps if set, -1 if not set.
// ==========================================================================================
real DispersiveMaterialParameters::
getAlphaP() const
{
  return alphaP;
}



// ==========================================================================================
/// \brief Return the isotropic parameters 
/// \param Np (output) : number Of polarization vectors
/// \param epsInf (output) : constant part of epsHat(s)
/// \param modelParams(0:3,0:Np-1) (output) : GDM parameters a0(j),a1(j),b0(j),b1(j), j=0,..,Np-1
// ==========================================================================================
int DispersiveMaterialParameters::
getIsotropicParameters( int & Np,  real & epsInf_, RealArray & modelParams ) const
{
  Np = numberOfPolarizationVectors;
  epsInf_=epsInf;
  modelParams.redim(modelParameters);
  modelParams=modelParameters;

  return 0;
}

// ==========================================================================================
/// \brief Return the parameters from nonlinear models 
/// \param numberOfNonlinearVariables (output) : number Of nonlinear unknowns
/// \param nlPar (output) : nonlinear parameters:
///    nlPar(0:numPar,0:numberOfNonlinearVariables-1) : 
/// 
// ==========================================================================================
int DispersiveMaterialParameters::
getNonlinearParameters( int & numberOfNonlinearVariables, RealArray & nlPar ) const 
{
  int & numberOfAtomicLevels = dbase.get<int>("numberOfAtomicLevels"); // for multilevelAtomic model 
  RealArray & nonlinearModelParameters = dbase.get<RealArray>("nonlinearModelParameters");


  numberOfNonlinearVariables = numberOfAtomicLevels;
  nlPar.redim(0);
  nlPar = nonlinearModelParameters;
  
  return 0;
  
}




// ===================================================================================================
/// \brief Return the bi-anisotropic Np array : number of polarization vectors for each entry of K
/// \param Np(0:5,0:5) (output) : number of polarization vectors for each element of K
// ===================================================================================================
IntegerArray& 
DispersiveMaterialParameters::getBianisotropicNp() const
{
  MaterialTypeEnum & materialType = dbase.get<MaterialTypeEnum>("materialType");
  if( materialType == bianisotropic )
  {
    if( ! dbase.has_key("bianisotropicParameters" ) )
    {
      printF("DispersiveMaterialParameters::getBianisotropicNp: bianisotropicParameters have not been assigned\n");
      return Overture::nullIntArray();
    }
   IntegerArray & NpBA                  = dbase.get<IntegerArray>("NpBA");
   return NpBA;
  }
  else
  {
    printF("DispersiveMaterialParameters::getBianisotropicNp:ERROR: materialType is not equal to bianisotropic\n");
    return Overture::nullIntArray();
  }
  
}

// ==========================================================================================
/// \brief Return the bi-anisotropic GDM parameters 
/// 
/// \param bianisotropicParameters(0:3,NpMax,0:5,0;5) : GDM parameters a0(j),a1(j),b0(j),b1(j), j=0,..,Np(k1,k2)-1
// ==========================================================================================
RealArray & DispersiveMaterialParameters::
getBianisotropicGDMParameters() const
{
  MaterialTypeEnum & materialType = dbase.get<MaterialTypeEnum>("materialType");
  if( materialType == bianisotropic )
  {
    if( ! dbase.has_key("bianisotropicParameters" ) )
    {
      printF("DispersiveMaterialParameters::getBianisotropicGDMParameters: bianisotropicParameters have not been assigned\n");
      return Overture::nullRealArray();
    }
    
   return dbase.get<RealArray>("bianisotropicParameters");

  }
  else
  {
    printF("DispersiveMaterialParameters::getBianisotropicGDMParameters:ERROR: materialType is not equal to bianisotropic\n");
    return Overture::nullRealArray();
  }

}



// ==========================================================================================
/// \brief Return the bi-anisotropic parameters 
/// \param K0(0:5,0:5) (output) : constant part of the BA material matrix
/// \param Np(0:5,0:5) (output) : number of polarization vectors for each element of K
/// \param bianisotropicParameters(0:3,NpMax,0:5,0;5) : GDM parameters a0(j),a1(j),b0(j),b1(j), j=0,..,Np(k1,k2)-1
// ==========================================================================================
int DispersiveMaterialParameters::
getBianisotropicParameters( RealArray & K0, RealArray & bianisotropicParameters, IntegerArray & Np )
{
  MaterialTypeEnum & materialType = dbase.get<MaterialTypeEnum>("materialType");
  if( materialType == bianisotropic )
  {
    if( ! dbase.has_key("bianisotropicParameters" ) )
    {
      printF("DispersiveMaterialParameters::getBianisotropicParametersERROR: bianisotropicParameters have not been assigned\n");
      return 1;
    }
    
   IntegerArray & NpBA                  = dbase.get<IntegerArray>("NpBA");
   RealArray & K0_                      = dbase.get<RealArray>("K0");
   RealArray & bianisotropicParameters_ = dbase.get<RealArray>("bianisotropicParameters");

   Np.redim(NpBA); Np=NpBA;
    
   K0.redim(K0_);    K0=K0_;

   bianisotropicParameters.redim(bianisotropicParameters_);
   bianisotropicParameters=bianisotropicParameters_;
   


  }
  else
  {
    printF("DispersiveMaterialParameters::getBianisotropicParameters:ERROR: materialType is not equal to bianisotropic\n");
    return 1;
  }

  return 0;
}


// ==========================================================================================
/// \brief Return the bi-anisotropic material matrix 
/// \param K0(0:5,0:5) (output) : constant part of the inverse BA material matrix
// ==========================================================================================
int DispersiveMaterialParameters::
getBianisotropicMaterialMatrix( RealArray & K0 ) const 
{
  MaterialTypeEnum & materialType = dbase.get<MaterialTypeEnum>("materialType");
  if( materialType == bianisotropic )
  {
    if( ! dbase.has_key("bianisotropicParameters" ) )
    {
      printF("DispersiveMaterialParameters::getBianisotropicParametersERROR: bianisotropicParameters have not been assigned\n");
      return 1;
    }
    K0.redim(6,6);
    K0=dbase.get<RealArray>("K0");
  }
  else
  {
    printF("DispersiveMaterialParameters::getBianisotropicMaterialMatrix:WARNING: materialType is not equal to bianisotropic.\n");
    printF("DMP: Returning K0 for an isotropic material with epsInf=%9.3e, muInf=%9.3e\n",epsInf,muInf);
    K0.redim(6,6);
    K0=0.;
    for( int i=0; i<3; i++ )
    {
      K0(i  ,i  )=epsInf;
      K0(i+3,i+3)=muInf;
    }
    


    return 1;
  }

  return 0;
}


  
// ==========================================================================================
/// \brief Return the inverse of the bi-anisotropic material matrix 
/// \param K0i(0:5,0:5) (output) : constant part of the inverse BA material matrix
// ==========================================================================================
int DispersiveMaterialParameters::
getBianisotropicMaterialMatrixInverse( RealArray & K0i )
{
  MaterialTypeEnum & materialType = dbase.get<MaterialTypeEnum>("materialType");
  if( materialType == bianisotropic )
  {
    if( ! dbase.has_key("bianisotropicParameters" ) )
    {
      printF("DispersiveMaterialParameters::getBianisotropicMaterialMatrixInverse: bianisotropicParameters have not been assigned\n");
      return 1;
    }

    if( !dbase.has_key("K0Inverse") )
    {
      // First time: compute the inverse of K0
      RealArray & K0Inverse = dbase.put<RealArray>("K0Inverse");
      
      K0Inverse.redim(6,6);

      RealArray & K0 = dbase.get<RealArray>("K0");

      // K0(3,0)=1.;
      // K0(0,1)=1.;
      
      K0Inverse=K0;

      int n=6, info=0;
      int ipvt[6];
      


      const int lwork=36;
      real work[lwork], anorm,rcond;
      int iwork[6];

      char norm = '1'; // 'I'; // '1'
      anorm = LANGE( "1", n, n, K0Inverse(0,0), n, work[0] ); // '1' = 1-norm 

      GETRF( n, n, K0Inverse(0,0), n, ipvt[0],  info );
      if( info!=0 )
      {
        printF("ERROR return from laplack routine GETRF, LU factorization: info=%d\n",info);
        OV_ABORT("GETRF: info !=0 ");
      }
      
      GECON( "1", n, K0Inverse(0,0), n, anorm, rcond, work[0], iwork[0], info ); // '1' = 1-norm 
      if( info!=0 )
        OV_ABORT("GECON: info !=0 ");

      printF(" 1-norm inverse condition number of K0 = %9.2e, anorm=%9.2e\n",rcond,anorm);


      GETRI( n, K0Inverse(0,0), n, ipvt[0], work[0], lwork, info );
      if( info!=0 )
        OV_ABORT("GETRI: info !=0 ");

      if( false )
      {
        ::display(K0,"K0");
        ::display(K0Inverse,"K0Inverse");
      }
      

    }
    
    K0i.redim(6,6);
    K0i = dbase.get<RealArray>("K0Inverse");

  }
  else
  {
    printF("DispersiveMaterialParameters::getBianisotropicMaterialMatrixInverse: ERROR: this is not a BA material!\n");
    OV_ABORT("Error");
  }
  
  return 0;
}



// ==========================================================================================
/// \brief Write current parameters to a material database file,
///
/// \param fileName (input) : name of material file  // write parameters to a material database file
// ==========================================================================================
int DispersiveMaterialParameters::
writeToFile( const aString & fileName )
{

  FILE *datFile = fopen((const char*)fileName, "w");
  if( datFile==NULL )
  {
    printF("DispersiveMaterialParameters::writeToFile:ERROR: unable to open fileName=[%s]\n",(const char*)fileName);
    OV_ABORT("error");
  }

  // --- Define the frequency scale --
  real V0 = dbase.get<real>("velocityScale");
  real L0 = dbase.get<real>("lengthScale"); 
  real T0 = L0/V0;  // time scale 
  
  // real Omega0 = 2.*Pi/T0;       // angular frequency scale : omega = c0*k = c0*(2*pi)/lambda = 2*pi * c0/L0 
  real Omega0 = 1./T0;             // angular frequency scales by this if we scale x/L0 and t/T0
  
  const MaterialTypeEnum & materialType = dbase.get<MaterialTypeEnum>("materialType");
  const aString & name = dbase.get<aString>("name");
  const real & omegaMin = dbase.get<real>("omegaMin");   
  const real & omegaMax = dbase.get<real>("omegaMax");   
  const real & omegaScale = dbase.get<real>("omegaScale");

  aString gdmDataFile="none";

  // Get the current date
  time_t *tp= new time_t;
  time(tp);
  // tm *ptm=localtime(tp);
  const char *dateString = ctime(tp);

  fPrintF(datFile,"# Dispersive material parameters file\n");
  fPrintF(datFile,"# File written by CgMx on %s",dateString);  // Note: dateString includes newline
  fPrintF(datFile,"units=\"MKS\";\n");
  fPrintF(datFile,"name=\"%s\";\n",(const char*)name);
  fPrintF(datFile,"gdmDataFile=\"%s\";\n",(const char*)gdmDataFile);
  fPrintF(datFile,"velocityScale=%20.14e;\n",V0);
  fPrintF(datFile,"lengthScale=%20.14e;\n",L0);
  fPrintF(datFile,"timeScale=%20.14e;\n",T0);
  fPrintF(datFile,"omegaScale=%20.14e;\n",omegaScale*Omega0);  // note scale by Omega0 -- is this right?
  //  min and max values for data 
  fPrintF(datFile,"omegaDataMin=%20.14e;\n",omegaMin);
  fPrintF(datFile,"omegaDataMax=%20.14e;\n",omegaMax);
  

  if( materialType==isotropic )
  {
    fPrintF(datFile,"# epsGDMPars%d=[epsInf,a0,a1,b0,b1, a0,a1,b0,b1, ...] (%d polarization vectors)\n",
	       numberOfPolarizationVectors,numberOfPolarizationVectors);
    fPrintF(datFile,"epsGDMPars%d=[%20.14e",numberOfPolarizationVectors,epsInf);
    for( int j=0; j<numberOfPolarizationVectors; j++ )
    {
      fPrintF(datFile,",%20.14e, %20.14e, %20.14e, %20.14e",
	      modelParameters(0,j),modelParameters(1,j),modelParameters(2,j),modelParameters(3,j));
    }
    fPrintF(datFile,"];\n");
  }
  else
  {
    OV_ABORT("DispersiveMaterialParameters::writeToFile: finish me");
  }
  
  fclose(datFile);
  printF("Wrote material data file=[%s]\n",(const char*)fileName);



  // Files were first written by the Matlab file plotGDMFit.m 

/* ----
  % ---------------------------------------------------------------------------
  % ---------------- Output a data file for C++ --------------------------------
  % ---------------------------------------------------------------------------
  % datFileName=sprintf('%s.gdm',name); 
  datFileName=sprintf('%sDispersionFits.txt',name); 
  datFile = fopen(datFileName,'w');

  fprintf(datFile,'# File written by plotGDMFit.m on %s\n',datestr(now));
  fprintf(datFile,'units="MKS";\n');
  fprintf(datFile,'name="%s";\n',name);
  fprintf(datFile,'gdmDataFile="%s";\n',gdmDataFile);
  fprintf(datFile,'omegaScale=%20.14e;\n',Omega0*dataScale);
  % min and max values for data 
  fprintf(datFile,'omegaDataMin=%20.14e;\n',min(omega));
  fprintf(datFile,'omegaDataMax=%20.14e;\n',max(omega));
  
 for( np=1:numFits )
    gdmPars = gdmFitData.gdmPar{np};
 
    fprintf(datFile,'epsGDM%dErrorNorm=%e;\n',np,gdmFitData.errorNorm{np});
    fprintf(datFile,'epsGDMPars%d=[%20.14e',np,gdmPars(1)*epsHatNorm);
    for ip=1:np
        i0 = 4*(ip-1) + 1;
        a0 = gdmPars(i0+1)*epsHatNorm;
        a1 = gdmPars(i0+2)*epsHatNorm;
        b0 = gdmPars(i0+3);
        b1 = gdmPars(i0+4);
	fprintf(datFile,',%20.14e, %20.14e, %20.14e, %20.14e',a0,a1,b0,b1);
    end;
    fprintf(datFile,'];\n');
     

  end % end for np 
  
  fclose(datFile);
  fprintf('Wrote text data file=[%s]\n',datFileName);

  ------ */

  return 0;
}


// ===============================================================================================================
/// \brief Fill in a matrix of values from an input string 
/// \param aline (input) : string containing comma separated values of the form 
///        aline(len,end)  = "[ val, val, val, ..., val ];" 
///  \param par (output) : fill in values par(n1a:n1b,n2a:n2b,n3a:n3b) 
/// \return value: 1=success, 0=failed to read all values. 
// ===============================================================================================================
int getMatrixValuesFromString( const aString & aline, int len,
			       int n1a, int n1b,
			       int n2a, int n2b,
			       int n3a, int n3b, 
			       RealArray & par )
{
  int debug=0; // 1 
  if( debug )
  {
    printF("getMatrixValuesFromString: aline=%s\n",(const char*)aline); 
    printF("getMatrixValuesFromString: [n1a,n1b][n2a,n2b][n3a,n3b]=[%d,%d][%d,%d][%d,%d]\n",n1a,n1b,n2a,n2b,n3a,n3b);
  }


  int length=aline.length();
  int ia=len+1;  // parameters start here
  for( int m3=n3a; m3<=n3b; m3++ )
  {
    for( int m1=n1a; m1<=n1b; m1++ )  // *NOTE* since it is a MATRIX and not an array, reverse order of m1 and m2 when reading
    {
      for( int m2=n2a; m2<=n2b; m2++ )
      {
        //  Read one value at a time ...

	int ie=ia+1;  // find ie = end of next number 
	while( ie<length && aline[ie] != ',' && aline[ie] != ']' )
	{
	  if( debug ) { printF(" aline[%d]=%c\n",ie,aline[ie]);  } 
	  
	  ie++;
	}
        if( debug ) { printF("read aline(ia,ie-1)=[%s]\n",(const char*)aline(ia,ie-1));;  } // 
	
	sScanF(aline(ia,ie-1),"%e",&par(m1,m2,m3));

	if( debug ) {  printF("par(%i,%i,%i)=%g\n",m1,m2,m3,par(m1,m2,m3));  } // 
	ia=ie+1;	
	if( ia>=length )
	{
	  printF("getMatrixValuesFromString:ERROR not enough values found from string=[%s]\n",(const char*)aline);
	  return 0;
	}
    
      }
      if( ia>=length )
	break;
    }
  }
  return 1;
}


// ==========================================================================================
/// \brief Read dispersive material parameters from a file.
///
/// \param fileName (input) : name of material file
/// \param numberOfPolarizationVectorsRequested (input): optionally specify the number of polarization vectors
 //   if more than one fit is available from the file (-1 = choose the largest number).
// ==========================================================================================
int DispersiveMaterialParameters::
readFromFile( const aString & fileName, int numberOfPolarizationVectorsRequested /* = -1 */ )
{
  const int myid = max(0,Communication_Manager::My_Process_Number);

  if( numberOfPolarizationVectorsRequested <=0 )
    numberOfPolarizationVectorsRequested=INT_MAX;  // choose maximum number available

  FILE *file = fopen((const char*)fileName, "r");
  if( file==NULL )
  {
    printF("DispersiveMaterialParameters::readFromFile:ERROR: unable to open fileName=[%s]\n",(const char*)fileName);
    OV_ABORT("error");
  }

  const int NpMax=10;  // **FIX ME**


  MaterialTypeEnum & materialType = dbase.get<MaterialTypeEnum>("materialType");
  if( materialType == bianisotropic )
  {
    //    K(6,6) : material tensor 
    //
    //     K0(6,6)  : constant part of material tensor 
    //     bianisotropicParameters(4,Np,6,6)        : GDM 
    //     Np(6,6) : number of polarization vectors 

    printF("DispersiveMaterialParameters::readFromFile: This is an bianisotropic material\n");

    if( ! dbase.has_key("bianisotropicParameters" ) )
    {
      dbase.put<IntegerArray>("NpBA");
      dbase.put<RealArray>("K0");
      dbase.put<RealArray>("bianisotropicParameters");
     
      RealArray & bianisotropicParameters = dbase.get<RealArray>("bianisotropicParameters");
      bianisotropicParameters.redim(4,NpMax,6,6);

      bianisotropicParameters=0.;
	
      IntegerArray & NpBA = dbase.get<IntegerArray>("NpBA");
      NpBA.redim(6,6); NpBA=0;
      
    }
    
  }
  

  aString & name = dbase.get<aString>("name");
  bool & isDispersive = dbase.get<bool>("isDispersive");
  isDispersive=false;
  
  aString & nonlinearModelName = dbase.get<aString>("nonlinearModelName");
  int & numberOfAtomicLevels = dbase.get<int>("numberOfAtomicLevels"); // for multilevelAtomic model 


  real & omegaMin = dbase.get<real>("omegaMin");   
  real & omegaMax = dbase.get<real>("omegaMax");   
  
  aString units;
  
  real & omegaScale = dbase.get<real>("omegaScale");
  RealArray gdmPars;
  real epsInf0;
  RealArray a0v(NpMax),a1v(NpMax),b0v(NpMax),b1v(NpMax);


  const real nm = 1e-9;         // nanometers  (meter-per-nm)
  const real um = 1e-6;         // micrometers
  const real c0 = 299792458;    // the speed of light, [m/c]
  // real L0 = 100*nm;              // length scale
  // real T0 = L0/c0;

  // --- Define the frequency scale --
  real V0 = dbase.get<real>("velocityScale");
  real L0 = dbase.get<real>("lengthScale"); 
  real T0 = L0/V0;  // time scale 
  
  // real Omega0 = 2.*Pi/T0;       // angular frequency scale : omega = c0*k = c0*(2*pi)/lambda = 2*pi * c0/L0 
  real Omega0 = 1./T0;             // angular frequency scales by this if we scale x/L0 and t/T0 

  real wScale=1.; // Frequency scale -- set below 

  // if true, do not scale GDM parameters by wScale = omegaScale/Omega0
  bool & normalizedUnits = dbase.get<bool>("normalizedUnits");  
  normalizedUnits=false;


  int currentBestNp=-1;
  
  int len;
  const int buffLength=2000;  // allow for continuation lines
  char line[buffLength];
  int numRead=getLineFromFile(file,line,buffLength);  // read a line from the file.

  bool parametersSet=false;
  GenericGraphicsInterface & gi = *Overture::getGraphicsInterface(); // needed for parser

  // Allow for empty lines in the file: 
  // while( numRead>0 )  
  bool parseForPerl=false;      // this is changed by the line "parser="perl" in the file. 
  bool semiColonEndsLine=true;  // Format V1, may be changed by file wit a line "format="V2.0".
  int shiftForSemiColon=1;
  
  while( !feof(file) )
  { 
    if( numRead>0 )
    {
      aString aline=line;
      printF("DispersiveMP: readFromFile: line=[%s]\n",(const char*)line);

      // Parse for perl variables in the file:
      if( parseForPerl )
      {
        int returnValue=gi.parseLine(aline);
        printF("DispersiveMP: after parsing: aline=[%s]\n",(const char*)aline);
      }
      
      if( line[0]=='#' )
	printF("Comment: %s\n",(const char*)line);
      else
      {
	if( len=aline.matches("format=") )
        {
          // The format= command defines the file format 
	  aString format=aline(len,aline.length()-1); 
	  int nl=format.length()-1;
	  if( format[0]=='"' && format[nl]=='"' )
	    format = format(1,nl-1);  // remove double quotes 
	  printF("format=[%s]\n",(const char*)format);
          if( format == "V2.0" )
          { // In the new format we do no end lines with semi-colon so we can parse lines with perl.
            printF("File format V2.0: do no end lines with semi-colon so we can parse lines with perl.\n");
            semiColonEndsLine=false;
            shiftForSemiColon=0;
          }
        }
        
	else if( len=aline.matches("parser=") )
        {
          // The parser="perl" will turn on parsing by perl.
          //     parser="none" will turn off parsing 
	  aString parser=aline(len,aline.length()-1); 
	  int nl=parser.length()-1;
	  if( parser[0]=='"' && parser[nl]=='"' )
	    parser = parser(1,nl-1);  // remove double quotes 
	  printF("parser=[%s]\n",(const char*)parser);
          if( parser == "perl" )
          { 
            printF("parser=perl: parse lines with a `;' for perl commands and substitute perl dollar variables.\n");
            parseForPerl=true;
          }
          else if( parser == "none" )
          {
            printF("Parser=none: lines will NOT be parsed for perl.\n");
            parseForPerl=false;
          }
          else
          {
            printF("Unknown parser=[%s]\n",(const char*)parser);
          }
          
        }
        
	else if( len=aline.matches("omegaScale=") )
	{
	  sScanF(aline(len,aline.length()-1),"%e",&omegaScale);
	  printF("omegaScale=%22.16e, Omega0=%22.16e\n",omegaScale,Omega0);

	  if( !normalizedUnits )
	  {
	    
	    wScale=omegaScale/Omega0; // Frequency scale 
	    printF("DispersiveMaterialParameters::readFromFile: non-dimensionalize: L0=%g(nm) = %g(m) T0=%g(s) Omega0=%g omegaScale=%g,  "
		   "wscale= omegaScale/Omega0 = omegaScale*T0/(2*pi)=%g\n", L0/nm,L0,T0,Omega0,omegaScale,wScale);
	  }
	  else
	  {
            printF("DMP: Ignoring omegaScale=%g since units=`normalizedUnits'\n",omegaScale);
	  }
	  
	}
	else if( len=aline.matches("omegaDataMin=") )
	{
	  sScanF(aline(len,aline.length()-1),"%e",&omegaMin);
	  printF("omegaMin=%22.16e\n",omegaMin);
	}
	else if( len=aline.matches("omegaDataMax=") )
	{
	  sScanF(aline(len,aline.length()-1),"%e",&omegaMax);
	  printF("omegaMax=%22.16e\n",omegaMax);
	}
      

	else if( len=aline.matches("name=") )
	{
	  name=aline(len,aline.length()-1-shiftForSemiColon); // skip final ";"
	  printF("DMP: name=%s\n",(const char*)name);
	  int nl=name.length()-1;
	  if( name[0]=='"' && name[nl]=='"' )
	    name = name(1,nl-1);  // remove double quotes 
	}
	else if( len=aline.matches("units=") )
	{
	  units=aline(len,aline.length()-1-shiftForSemiColon); // skip final ";"
          int nl=units.length()-1;
	  if( units[0]=='"' && units[nl]=='"' )
	    units = units(1,nl-1);  // remove double quotes 

          if( units == "normalized" )  
	  {
	    printF("DMP::INFO: units=%s : GDM parameters will NOT be scaled by a reference frequency.\n",(const char*)units);
	    normalizedUnits=true;   // do not scale GDM parameters
	    wScale=1.;
	  }
	  else
	    printF("DMP: units=[%s].\n",(const char*)units);

	}
        else if( len=aline.matches("nonlinearModel=") )
	{
	  nonlinearModelName=aline(len,aline.length()-1-shiftForSemiColon); // skip final ";"
          int nl=nonlinearModelName.length()-1;
	  if( nonlinearModelName[0]=='"' && nonlinearModelName[nl]=='"' )
	    nonlinearModelName = nonlinearModelName(1,nl-1);  // remove double quotes 

	  printF("Setting nonlinearModel=[%s]\n",(const char*)nonlinearModelName);
	  if( nonlinearModelName == "multilevelAtomic" )
	  {
	    nonlinearModel=multilevelAtomic;
	    dbase.get<bool>("isNonlinearMaterial")=true;
	  }
	  else
	  {
            printF("DispersiveMaterialParameters:WARNING: IGNORING UNKNOWN NONLINEAR MDOEL model=[%s]]\\n",(const char*)nonlinearModelName);
	  }

	}
        else if( len=aline.matches("numberOfAtomicLevels=") )
	{
           sScanF(aline(len,aline.length()-1),"%i",&numberOfAtomicLevels);
	   printF("DMP: Setting numberOfAtomicLevels=%d\n",numberOfAtomicLevels);
	}
        else if( len=aline.matches("numberOfPolarizationVectors=") )
	{
	  // **do this for now for the MLA model ***
           sScanF(aline(len,aline.length()-1),"%i",&numberOfPolarizationVectors);
	   printF("DMP: Setting numberOfPolarizationVectors=%d\n",numberOfPolarizationVectors);

	   if( numberOfPolarizationVectors > numberOfPolarizationVectorsRequested )
	   {
	     printF("\n ++++++ ERROR: numberOfPolarizationVectors=%d (from file) is greater than numberOfPolarizationVectorsRequested=%d.\n"
		    " You should change numberOfPolarizationVectorsRequested to %d.\n",
		    numberOfPolarizationVectors,numberOfPolarizationVectorsRequested,numberOfPolarizationVectors);
	     OV_ABORT("DMP:ERROR");
	   }
	   

	}

        else if( (len=aline.matches("polarizationNECoefficients=")) || 
	         (len=aline.matches("populationRelaxationCoefficients="))|| 
	         (len=aline.matches("populationEPtCoefficients="))  )
	{
	  int & numberOfAtomicLevels = dbase.get<int>("numberOfAtomicLevels"); // for multilevelAtomic model 
	  RealArray & nonlinearModelParameters = dbase.get<RealArray>("nonlinearModelParameters");

          // We store 3 arrays of coefficients  ** FIX ME ***
	  //    P[ij]_tt + ...   =    alpha(i,j) * P[ij]_t E_j 
	  //   alpha(i,j), beta(i,j), gamma(i,j)   i,j = 0,1,...,numberOfAtomicLevels

          const int numberOfCoefficientTensors = 3; // a, alpha, beta 
          const int maxCoeff = max(numberOfAtomicLevels,numberOfPolarizationVectors);
          if( nonlinearModelParameters.getLength(0) != maxCoeff )
	  {
	    nonlinearModelParameters.redim(maxCoeff,maxCoeff,numberOfCoefficientTensors);
            nonlinearModelParameters=0.;
	  }
	  
          if( len=aline.matches("polarizationNECoefficients=") )
	  {
            printF("Reading polarizationNECoefficients matrix: numberOfPolarizationVectors=%d, numberOfAtomicLevels=%d\n",
		   numberOfPolarizationVectors,numberOfAtomicLevels );
	
            // read matrix values:
	    //     nonlinearModelParameters(0:np-1,0:na-1,0:0) 
            int ok= getMatrixValuesFromString( aline,len,
					       0,numberOfPolarizationVectors-1,  // range for first component
					       0,numberOfAtomicLevels-1,         // range for second component
					       0,0,                              // range for third component 
					       nonlinearModelParameters );
            Range R1=numberOfPolarizationVectors, R2=numberOfAtomicLevels;
	    ::display(nonlinearModelParameters(R1,R2,0),"polarizationNECoefficients (TRANSPOSE)","%6.2f ");
	    
            if( !ok )
	    {
	      printf("DMP:readFromFile:ERROR reading matrix values for polarizationNECoefficients!\n");
	    }
	    // OV_ABORT("polarizationNECoefficients: stop here for now");
	    


	  }
	  else if( len=aline.matches("populationRelaxationCoefficients=") )
	  {
            
            printF("Reading populationRelaxationCoefficients matrix: numberOfAtomicLevels=%d\n",
		   numberOfAtomicLevels );
	
            // read matrix values:
	    //     nonlinearModelParameters(0:np-1,0:na-1,0:0) 
            int ok= getMatrixValuesFromString( aline,len,
					       0,numberOfAtomicLevels-1,         // range for first component
					       0,numberOfAtomicLevels-1,         // range for second component
					       1,1,                              // range for third component *NOTE*
					       nonlinearModelParameters );
            Range R1=numberOfAtomicLevels, R2=numberOfAtomicLevels;
	    ::display(nonlinearModelParameters(R1,R2,1),"populationRelaxationCoefficients (TRANSPOSE)","%6.2f ");
	    
            if( !ok )
	    {
	      printf("DMP:readFromFile:ERROR reading matrix values for populationRelaxationCoefficients!\n");
	    }
	    // OV_ABORT("populationRelaxationCoefficients: stop here for now");

	  }
	  else if( len=aline.matches("populationEPtCoefficients=") )
	  {
            
            printF("Reading populationEPtCoefficients matrix: numberOfAtomicLevels=%d, numberOfPolarizationVectors=%d\n",
		   numberOfAtomicLevels,numberOfPolarizationVectors );
	
            // read matrix values:
	    //     nonlinearModelParameters(0:na-1,0:np-1,0:0) 
            int ok= getMatrixValuesFromString( aline,len,
					       0,numberOfAtomicLevels-1,            // range for first component
					       0,numberOfPolarizationVectors-1,     // range for second component
					       2,2,                                 // range for third component *NOTE*
					       nonlinearModelParameters );
            Range R1=numberOfAtomicLevels, R2=numberOfPolarizationVectors;
	    ::display(nonlinearModelParameters(R1,R2,2),"populationEPtCoefficients (TRANSPOSE)","%6.2f ");
	    
            if( !ok )
	    {
	      printf("DMP:readFromFile:ERROR reading matrix values for populationEPtCoefficients!\n");
	    }
	    // OV_ABORT("populationEPtCoefficients: stop here for now");

	  }

	  else
	  {
	    OV_ABORT("ERROR -- this should not happen!");
	  }
	  

	  
	}
	

	else if( len=aline.matches("epsGDMPars") )	
	{
	  // "epsGDMPars1=["
	  // "epsGDMPars2=["
	  // "epsGDMPars3=["

	  // GDM parameters
          parametersSet=true;  // parameters have been read
	  isDispersive=true;
	  
	  int Np=1;
	  sScanF(aline(len,len+1),"%i",&Np);
          if( Np==0 )
          {
            // non-dispersive material parameters
            sScanF(aline(len+3,aline.length()-1),"%e",&epsInf0);
          }
          else if( Np==1 )
	    sScanF(aline(len+3,aline.length()-1),"%e %e %e %e %e",&epsInf0,&a0v(0),&a1v(0),&b0v(0),&b1v(0));
	  else if( Np==2 )
	    sScanF(aline(len+3,aline.length()-1),"%e %e %e %e %e %e %e %e %e",&epsInf0,
		   &a0v(0),&a1v(0),&b0v(0),&b1v(0),
		   &a0v(1),&a1v(1),&b0v(1),&b1v(1));
	  else if( Np==3 )
	    sScanF(aline(len+3,aline.length()-1),"%e %e %e %e %e %e %e %e %e %e %e %e %e",&epsInf0,
		   &a0v(0),&a1v(0),&b0v(0),&b1v(0),
		   &a0v(1),&a1v(1),&b0v(1),&b1v(1),
		   &a0v(2),&a1v(2),&b0v(2),&b1v(2));
	  else if( Np==4 )
	    sScanF(aline(len+3,aline.length()-1),"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e",&epsInf0,
		   &a0v(0),&a1v(0),&b0v(0),&b1v(0),
		   &a0v(1),&a1v(1),&b0v(1),&b1v(1),
		   &a0v(2),&a1v(2),&b0v(2),&b1v(2),
		   &a0v(3),&a1v(3),&b0v(3),&b1v(3));
	  else if( Np==5 )
	    sScanF(aline(len+3,aline.length()-1),"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e",&epsInf0,
		   &a0v(0),&a1v(0),&b0v(0),&b1v(0),
		   &a0v(1),&a1v(1),&b0v(1),&b1v(1),
		   &a0v(2),&a1v(2),&b0v(2),&b1v(2),
		   &a0v(3),&a1v(3),&b0v(3),&b1v(3),
		   &a0v(4),&a1v(4),&b0v(4),&b1v(4));
	  else
	  {
	    printF("DispersiveMaterialParameters:readFromFile:WARNING: finish me for -- Np=%d",Np);
	  }
	  printF("Np=%d, epsInf0=%g\n",Np,epsInf0);
	  for( int j=0; j<Np; j++ )
	    printF("j=%d: [a0,a1,b0,b1]=[%g,%g,%g,%g]\n",j,a0v(j),a1v(j),b0v(j),b1v(j));


	  // -- See bamx.pdf notes for scaling ---
	  // MKS (Hertz) scaling of parameters would be
	  //    a0*omegaScale^2
	  //    a1*omegaScale
	  //    b0*omegaScale^2
	  //    b1*omegaScale
	  //
	  // We define a scaled s
	  //      sTilde = s/Omega0
	  //   Omega0 = 2*pi*c0/L0 = 2*pi/T0
	  // and then the GDM parameters become
	  //    a0*omegaScale^2/Omega0^2  = a0*wScale^2 
	  //    a1*omegaScale/Omega0      = a1*wScale
	  //    b0*omegaScale^2/Omega0^2  = b0*wScale^2 
	  //    b1*omegaScale/Omega0      = b1*wScale
	  // where
	  //    wScale = omegaScale/Omega0 = omegaScale*T0/(2*pi)

	  // printF("\n DMP ***** TEMP******* set wScale=1 (was %9.2e)\n",wScale);
	  // wScale=1; // ********** TEMP ***********
        
	  if( normalizedUnits )
	     wScale=1.;
	   
          if( Np <= numberOfPolarizationVectorsRequested && Np>currentBestNp )
	  {
	    currentBestNp=Np;
	  
	    epsInf=epsInf0;
	    alphaP=1./epsInf;

	    numberOfPolarizationVectors=Np;
	    printF("\n DispersiveMaterialParameters: SETTING numberOfPolarizationVectors=%d\n\n",numberOfPolarizationVectors);

            if( Np==0 )
            {
              isDispersive=false;
            }
            else if( Np>0 )
            {
              modelParameters.redim(4,Np);
              for( int j=0; j<Np; j++ )
              {
                // Check the scaling
                // modelParameters(0,j)=a0v(j)*(wScale*wScale)/epsInf; 
                // modelParameters(1,j)=a1v(j)*(       wScale)/epsInf;
                // modelParameters(2,j)=b0v(j)*(wScale*wScale)/epsInf;
                // modelParameters(3,j)=b1v(j)*(       wScale)/epsInf;

                modelParameters(0,j)=a0v(j)*(wScale*wScale); 
                modelParameters(1,j)=a1v(j)*(       wScale);
                modelParameters(2,j)=b0v(j)*(wScale*wScale);
                modelParameters(3,j)=b1v(j)*(       wScale);
              }
            }
	  }
	

	}
	else if( len=aline.matches("dispersiveModel=") )
	{
	  aString dm = aline(len,aline.length()-1);
	  if( dm == "\"GDM\";" || dm == "\"GDM\"" || dm == "GDM" )
	  {
	    printF("dispersiveModel = GDM. Setting isDispersive=true\n");
	    // isDispersive=true;
	  }
          printF("\n DMP:readFromFile: aline=[%s] dm=[%s]\n",(const char*)aline,(const char*)dm);
          if( dm=="\"none\"" || dm=="\"none\";" )
          {
            isDispersive=false;
          }
          
	  
	}


	else if( len=aline.matches("bianisotropicK0")   )	
	{
	  // ---- BI-ANISOTROPIC MATERIAL MATRIX K0 ---
	  //   bianisotropicK0=[...]; 
	  // K0 matrix

          parametersSet=true;
          
          if( ! dbase.has_key("K0") )
	    dbase.put<RealArray>("K0");

	  RealArray & K0 = dbase.get<RealArray>("K0");
	  K0.redim(6,6);
	
	  printF("Reading the bianisotropic K0 matrix\n");
	
	  // --- NOTE: there are too many input parameters for sScanF to read 36 values at once
	  //     so read one value at a time ...
	  int length=aline.length();
	  int ia=len+2;  // parameters start here
	  for( int k2=0; k2<6; k2++ )
	  {
	    for( int k1=0; k1<6; k1++ )
	    {
	      int ie=ia+1;  // ie = end of number 
	      while( ie<length && aline[ie] != ',' && aline[ie] != ']' )
	      {
		// printF(" aline[%d]=%c\n",ie,aline[ie]);
		ie++;
	      }
	      sScanF(aline(ia,ie-1),"%e",&K0(k1,k2));

	      // printF("K0(%i,%i)=%g\n",k1,k2,K0(k1,k2));
	      ia=ie+1;	
	      if( ia>=length )
	      {
		printF("ERROR reading K0 from string=[%s]\n",(const char*)aline);
		break;
	      }
    
	    }
	    if( ia>=length )
	      break;
	  }
	
          if( materialType == isotropic )
	  {
            // Isotropic material: 
	    epsInf= K0(0,0);
            muInf = K0(3,3);
	    alphaP=1./epsInf;

            printF("DMP:readFromFile: Setting isotropic epsInf = K0(0,0) = %g, muInf = K0(3,3) = %g \n",epsInf,muInf);
   
	  }
	  
          if( myid==0 && false  )
            ::display(K0,"K0");
	}
      
	else if( len=aline.matches("bianistropicK0")  )	
	{
          // **** OLD VERSION FOR BACKWARD COMPAT -- bianistropic is SPELLED WRONG !! ****


	  // ---- BI-ANISTROPIC MATERIAL MATRIX K0 ---
	  //   bianistropicK0=[...]; 
	  // K0 matrix

          parametersSet=true;
          
          if( ! dbase.has_key("K0") )
	    dbase.put<RealArray>("K0");

	  RealArray & K0 = dbase.get<RealArray>("K0");
	  K0.redim(6,6);
	
	  printF("Reading the bianistropic K0 matrix\n");
	
	  // --- NOTE: there are too many input parameters for sScanF to read 36 values at once
	  //     so read one value at a time ...
	  int length=aline.length();
	  int ia=len+2;  // parameters start here
	  for( int k2=0; k2<6; k2++ )
	  {
	    for( int k1=0; k1<6; k1++ )
	    {
	      int ie=ia+1;  // ie = end of number 
	      while( ie<length && aline[ie] != ',' && aline[ie] != ']' )
	      {
		// printF(" aline[%d]=%c\n",ie,aline[ie]);
		ie++;
	      }
	      sScanF(aline(ia,ie-1),"%e",&K0(k1,k2));

	      // printF("K0(%i,%i)=%g\n",k1,k2,K0(k1,k2));
	      ia=ie+1;	
	      if( ia>=length )
	      {
		printF("ERROR reading K0 from string=[%s]\n",(const char*)aline);
		break;
	      }
    
	    }
	    if( ia>=length )
	      break;
	  }
	
          if( materialType == isotropic )
	  {
            // Isotropic material: 
	    epsInf= K0(0,0);
            muInf = K0(3,3);
	    alphaP=1./epsInf;

            printF("DMP:readFromFile: Setting isotropic epsInf = K0(0,0) = %g, muInf = K0(3,3) = %g \n",epsInf,muInf);
   
	  }
	  
          if( myid==0 && false  )
            ::display(K0,"K0");
	}
      

	else if( len=aline.matches("bianisotropicPars") )	
	{
	  // ---- BI-ANISOTROPIC MATERIAL GDM COEFFICIENTS ---

	  // Entry form: 
	  // "bianisotropicPars[1-6][1-6]GDM[1-]=["

	  isDispersive=true;
	  
	  // bianisotropicParameters(4,Np,6,6) 
          if( ! dbase.has_key("bianisotropicParameters") )
	  {
	    dbase.put<RealArray>("bianisotropicParameters");
            RealArray & bianisotropicParameters = dbase.get<RealArray>("bianisotropicParameters");
	    bianisotropicParameters.redim(4,NpMax,6,6);
            bianisotropicParameters=0.;

	    dbase.put<IntegerArray>("NpBA");
  	    IntegerArray & NpBA = dbase.get<IntegerArray>("NpBA");
            NpBA.redim(6,6); NpBA=0;

	  }

	  RealArray & bianisotropicParameters = dbase.get<RealArray>("bianisotropicParameters");
	  IntegerArray & NpBA = dbase.get<IntegerArray>("NpBA");

	  // GDM parameters
	  if( normalizedUnits )
	     wScale=1.;
	   
	  int k1=1,k2=1;
	  sScanF(aline(len  ,len  ),"%i1",&k1);
	  sScanF(aline(len+1,len+1),"%i1",&k2);

	  assert( k1>=1 && k1<=6 );
	  assert( k2>=1 && k2<=6 );
	  int Np; 
	  // printf("aline(len+5,len+5)=[%s]\n",(const char*)aline(len+5,len+5));
	
	  sScanF(aline(len+5,len+5),"%i1",&Np);

	  if( true )
	    printF("Reading GDM parameters for material matrix K%d%d, Np=%d, numberOfPolarizationVectorsRequested=%d.\n",
		   k1,k2,Np,numberOfPolarizationVectorsRequested);

	  // --- Read the GDM data ----
	  if( Np <= numberOfPolarizationVectorsRequested )
	  {
	  
	    printF("Reading GDM parameters for material matrix K%d%d, Np=%d (scaling factor wScale=%g).\n",k1,k2,Np,wScale);
	
	    int ia=len+8;  // parameters start here
	    if( Np==1 )
	      sScanF(aline(ia,aline.length()-1),"%e %e %e %e",&a0v(0),&a1v(0),&b0v(0),&b1v(0));
	    else if( Np==2 )
	      sScanF(aline(ia,aline.length()-1),"%e %e %e %e %e %e %e %e",
		     &a0v(0),&a1v(0),&b0v(0),&b1v(0),
		     &a0v(1),&a1v(1),&b0v(1),&b1v(1));
	    else if( Np==3 )
	      sScanF(aline(ia,aline.length()-1),"%e %e %e %e %e %e %e %e %e %e %e %e",
		     &a0v(0),&a1v(0),&b0v(0),&b1v(0),
		     &a0v(1),&a1v(1),&b0v(1),&b1v(1),
		     &a0v(2),&a1v(2),&b0v(2),&b1v(2));
	    else if( Np==4 )
	      sScanF(aline(ia,aline.length()-1),"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e",
		     &a0v(0),&a1v(0),&b0v(0),&b1v(0),
		     &a0v(1),&a1v(1),&b0v(1),&b1v(1),
		     &a0v(2),&a1v(2),&b0v(2),&b1v(2),
		     &a0v(3),&a1v(3),&b0v(3),&b1v(3));
	    else if( Np==5 )
	      sScanF(aline(ia,aline.length()-1),"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e",
		     &a0v(0),&a1v(0),&b0v(0),&b1v(0),
		     &a0v(1),&a1v(1),&b0v(1),&b1v(1),
		     &a0v(2),&a1v(2),&b0v(2),&b1v(2),
		     &a0v(3),&a1v(3),&b0v(3),&b1v(3),
		     &a0v(4),&a1v(4),&b0v(4),&b1v(4));
	    else
	    {
	      printF("DispersiveMaterialParameters:readFromFile:WARNING: finish me for -- Np=%d",Np);
	    }
	    NpBA(k1-1,k2-1)=Np;
	
	    // printF("k1=%d, k2=%d, Np=%d\n",k1,k2,Np);
	    for( int j=0; j<Np; j++ )
	    {
	      assert( j<bianisotropicParameters.getLength(1) );
	  
	      bianisotropicParameters(0,j,k1-1,k2-1) = a0v(j)*(wScale*wScale);
	      bianisotropicParameters(1,j,k1-1,k2-1) = a1v(j)*(wScale       );
	      bianisotropicParameters(2,j,k1-1,k2-1) = b0v(j)*(wScale*wScale);
	      bianisotropicParameters(3,j,k1-1,k2-1) = b1v(j)*(wScale       );
	  
	      printF("j=%d: [a0,a1,b0,b1]=[%g,%g,%g,%g] -> scaled =[%g,%g,%g,%g]\n",
		     j,a0v(j),a1v(j),b0v(j),b1v(j),a0v(j)*(wScale*wScale),a1v(j)*wScale,b0v(j)*(wScale*wScale),b1v(j)*wScale);
	    }
	

	    
	  } // end if Np <= 
	

	}  // end if( len=aline.matches("bianisotropicPars") )

	else if( len=aline.matches("bianistropicPars") )	
	{
          // **** OLD VERSION FOR BACKWARD COMPAT -- bianistropic is SPELLED WRONG !! ****

	  // ---- BI-ANISTROPIC MATERIAL GDM COEFFICIENTS ---

	  // Entry form: 
	  // "bianistropicPars[1-6][1-6]GDM[1-]=["

	  isDispersive=true;
	  
	  // bianisotropicParameters(4,Np,6,6) 
          if( ! dbase.has_key("bianisotropicParameters") )
	  {
	    dbase.put<RealArray>("bianisotropicParameters");
            RealArray & bianisotropicParameters = dbase.get<RealArray>("bianisotropicParameters");
	    bianisotropicParameters.redim(4,NpMax,6,6);
            bianisotropicParameters=0.;

	    dbase.put<IntegerArray>("NpBA");
  	    IntegerArray & NpBA = dbase.get<IntegerArray>("NpBA");
            NpBA.redim(6,6); NpBA=0;

	  }

	  RealArray & bianisotropicParameters = dbase.get<RealArray>("bianisotropicParameters");
	  IntegerArray & NpBA = dbase.get<IntegerArray>("NpBA");

	  // GDM parameters
	  if( normalizedUnits )
	     wScale=1.;
	   
	  int k1=1,k2=1;
	  sScanF(aline(len  ,len  ),"%i1",&k1);
	  sScanF(aline(len+1,len+1),"%i1",&k2);

	  assert( k1>=1 && k1<=6 );
	  assert( k2>=1 && k2<=6 );
	  int Np; 
	  // printf("aline(len+5,len+5)=[%s]\n",(const char*)aline(len+5,len+5));
	
	  sScanF(aline(len+5,len+5),"%i1",&Np);

	  if( true )
	    printF("Reading GDM parameters for material matrix K%d%d, Np=%d, numberOfPolarizationVectorsRequested=%d\n",
		   k1,k2,Np,numberOfPolarizationVectorsRequested);

	  // --- Read the GDM data ----
	  if( Np <= numberOfPolarizationVectorsRequested )
	  {
	  
	    printF("Reading GDM parameters for material matrix K%d%d, Np=%d\n",k1,k2,Np);
	
	    int ia=len+8;  // parameters start here
	    if( Np==1 )
	      sScanF(aline(ia,aline.length()-1),"%e %e %e %e",&a0v(0),&a1v(0),&b0v(0),&b1v(0));
	    else if( Np==2 )
	      sScanF(aline(ia,aline.length()-1),"%e %e %e %e %e %e %e %e",
		     &a0v(0),&a1v(0),&b0v(0),&b1v(0),
		     &a0v(1),&a1v(1),&b0v(1),&b1v(1));
	    else if( Np==3 )
	      sScanF(aline(ia,aline.length()-1),"%e %e %e %e %e %e %e %e %e %e %e %e",
		     &a0v(0),&a1v(0),&b0v(0),&b1v(0),
		     &a0v(1),&a1v(1),&b0v(1),&b1v(1),
		     &a0v(2),&a1v(2),&b0v(2),&b1v(2));
	    else if( Np==4 )
	      sScanF(aline(ia,aline.length()-1),"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e",
		     &a0v(0),&a1v(0),&b0v(0),&b1v(0),
		     &a0v(1),&a1v(1),&b0v(1),&b1v(1),
		     &a0v(2),&a1v(2),&b0v(2),&b1v(2),
		     &a0v(3),&a1v(3),&b0v(3),&b1v(3));
	    else if( Np==5 )
	      sScanF(aline(ia,aline.length()-1),"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e",
		     &a0v(0),&a1v(0),&b0v(0),&b1v(0),
		     &a0v(1),&a1v(1),&b0v(1),&b1v(1),
		     &a0v(2),&a1v(2),&b0v(2),&b1v(2),
		     &a0v(3),&a1v(3),&b0v(3),&b1v(3),
		     &a0v(4),&a1v(4),&b0v(4),&b1v(4));
	    else
	    {
	      printF("DispersiveMaterialParameters:readFromFile:WARNING: finish me for -- Np=%d",Np);
	    }
	    NpBA(k1-1,k2-1)=Np;
	
	    // printF("k1=%d, k2=%d, Np=%d\n",k1,k2,Np);
	    for( int j=0; j<Np; j++ )
	    {
	      assert( j<bianisotropicParameters.getLength(1) );
	  
	      bianisotropicParameters(0,j,k1-1,k2-1) = a0v(j)*(wScale*wScale);
	      bianisotropicParameters(1,j,k1-1,k2-1) = a1v(j)*(wScale       );
	      bianisotropicParameters(2,j,k1-1,k2-1) = b0v(j)*(wScale*wScale);
	      bianisotropicParameters(3,j,k1-1,k2-1) = b1v(j)*(wScale       );
	  
	      printF("j=%d: [a0,a1,b0,b1]=[%g,%g,%g,%g] -> scaled =[%g,%g,%g,%g]\n",
		     j,a0v(j),a1v(j),b0v(j),b1v(j),a0v(j)*(wScale*wScale),a1v(j)*wScale,b0v(j)*(wScale*wScale),b1v(j)*wScale);
	    }
	

	    
	  } // end if Np <= 
	

	}  // end if( len=aline.matches("bianistropicPars") )



	else
	{
	  printF("Ignored: line=[%s]\n",(const char*)line);      
	}
      }
    } // end if numRead>0 

    numRead=getLineFromFile(file,line,buffLength);

  } // end if !feof(file) 
     

  if( nonlinearModel == multilevelAtomic )
  {

    IntegerArray & NpBA = dbase.get<IntegerArray>("NpBA");

    printF("DispersiveMaterialParameters: nonlinearModel=multilevelAtomic  NpBA(0,0)=%d, numberOfPolarizationVectors=%d\n",
	   NpBA(0,0),numberOfPolarizationVectors);

    if( NpBA(0,0) != numberOfPolarizationVectors )
    {
      printF("\n +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
      printF("DispersiveMaterialParameters:ERROR reading file=[%s]\n",(const char*)fileName);
      printF("DMP: nonlinearModel = multilevelAtomic, but numberOfPolarizationVectors=%d is not equal to NpBA(0,0)=%d\n",
	     numberOfPolarizationVectors, NpBA(0,0));
      printF("DMP: Something is wrong here. Maybe numberOfPolarizationVectorsRequested=%d is wrong.\n",
	     numberOfPolarizationVectorsRequested);
      OV_ABORT("DMP::error");
    }
    
  }
  

  if( materialType==isotropic && dbase.has_key("bianisotropicParameters") )
  {
    // Fill in isotropic parameters
    RealArray & bianisotropicParameters = dbase.get<RealArray>("bianisotropicParameters");
    IntegerArray & NpBA = dbase.get<IntegerArray>("NpBA");

    printf("DMP::readFromFile:Setting [a0,a1,b0,b1] from BA parameters for an non-BA material.\n");
    printf("DMP::INFO scaling GDM parameters a0,a1 by epsInf=%g since \n"
	   "        K(s)   = K0 + K1(s)\n"
	   "        eps(s) = epsInf*( 1 + chi(s) ) = epsInf + epsInf*chi(s) \n"
	   "   =>   chi(s) = K1(s)/epsInf \n",epsInf);
    int k1=1, k2=1;
    int Np = NpBA(k1-1,k2-1);
    numberOfPolarizationVectors = Np;

    modelParameters.redim(4,Np);

    for( int j=0; j<Np; j++ )
    {
      assert( j<bianisotropicParameters.getLength(1) );
	  
      a0v(j) = bianisotropicParameters(0,j,k1-1,k2-1)/epsInf;
      a1v(j) = bianisotropicParameters(1,j,k1-1,k2-1)/epsInf;
      b0v(j) = bianisotropicParameters(2,j,k1-1,k2-1);
      b1v(j) = bianisotropicParameters(3,j,k1-1,k2-1);
	  
      modelParameters(0,j)=a0v(j);
      modelParameters(1,j)=a1v(j);
      modelParameters(2,j)=b0v(j);
      modelParameters(3,j)=b1v(j);


      printF("j=%d: [a0,a1,b0,b1]=[%g,%g,%g,%g]\n",j,a0v(j),a1v(j),b0v(j),b1v(j));
    }
  }


  if( !parametersSet )
  {
    printF("DispersiveMaterialParameters::readFromFile:ERROR: NO GDM parameters found in fileName=[%s]\n",(const char*)fileName);
    OV_ABORT("error");
  }
  

  fclose(file);

  // ----- sanity check ---
  if( numberOfPolarizationVectors>0 && !isDispersiveMaterial() )
  {
    printF("DispersiveMaterialParameters::readFromFile:ERROR: numberOfPolarizationVectors=%d >0 BUT isDispersiveMaterial()=%d\n",
           numberOfPolarizationVectors,(int)isDispersiveMaterial());
    OV_ABORT("error");
  }   

  return 0;
  
}

// ==========================================================================================
/// \brief  Set a velocity and length scale for non-dimensionalization of the parameters .
/// \param velocityScale (input): scale for the velocity, e.g. speed of light in a vacuum.
/// \param lengthScale (input): length scale, e.g. 100 nm = 10^(-7) 
// ==========================================================================================
int DispersiveMaterialParameters::
setScales( real velocityScale, real lengthScale )
{
  dbase.get<real>("velocityScale")=velocityScale;
  dbase.get<real>("lengthScale")  =lengthScale;
  
  return 0;
}



// ==========================================================================================
/// \brief Set the material type.
/// \param matType (input): isotropic, bianisotropic
// ==========================================================================================
int DispersiveMaterialParameters::
setMaterialType( MaterialTypeEnum matType ) 
{
  dbase.get<MaterialTypeEnum>("materialType")=matType;
  
  return 0;
}

// ==========================================================================================
/// \brief return the material type.
// ==========================================================================================
DispersiveMaterialParameters::MaterialTypeEnum DispersiveMaterialParameters::
getMaterialType() const 
{
  return dbase.get<MaterialTypeEnum>("materialType");
}



// ==========================================================================================
/// \brief Specify the number of polarization vectors (number of GDM equation)
// ==========================================================================================
int DispersiveMaterialParameters::
setNumberOfPolarizationVectors(  const int numPolarizationVectors )
{
  if( numberOfPolarizationVectors!=numPolarizationVectors )
  {
    numberOfPolarizationVectors=numPolarizationVectors;
    numberOfModelParameters=4;     
    modelParameters.redim(numberOfModelParameters,numberOfPolarizationVectors); 
    modelParameters=0.;

    chir0.redim(numberOfPolarizationVectors); chir0=0.;
    chii0.redim(numberOfPolarizationVectors); chii0=0.;
    
    dbase.get<bool>("isDispersive")=true;

  }
  return 0;
}

// ==========================================================================================
/// \brief Specify the mode (i.e. the root of the dispersion relation to use for exact solutions.
/// \parammodeToCHoose (input) : mode number=0,1,2,... . Choose modeToChoose=-1 for default.
///     The default root is the one with largest imaginary part. 
// ==========================================================================================
int DispersiveMaterialParameters::
setMode( const int modeToChoose )
{
  mode=modeToChoose;
  
  return 0;
}

// ==========================================================================================
/// \brief Set the parameters in the GDM model for equation "eqn"
/// \param a0,a1,b0,b1 (input)
/// 
/// Generalized Dispersion Model:
///       E_tt - c^2 Delta(E) = -alphaP P_tt
///       Pi_tt + b1i Pi_1 + b0i = a0i*E + a1i*E_t    i=0,1,2,...,numPolarVectors-1
// ==========================================================================================
/// \brief Specify the number of polarization vectors (number of GDM equation)
// ==========================================================================================
int DispersiveMaterialParameters::
setParameters( const int eqn, const real a0, const real a1, const real b0, const real b1 )
{
  if( eqn<0 || eqn>=numberOfPolarizationVectors )
  {
    printF("DispersiveMaterialParameters::setParameters:ERROR: Trying to GDM eqn=%i, "
           "but numberOfPolarizationVectors=%i\n",eqn,numberOfPolarizationVectors);
    return 1;
    
  }

  printF("DispersiveMaterialParameters::setParameters: Setting GDM parameters eqn=%i: "
         "a0=%9.3e, a1=%9.3e, b0=%9.3e, b1=%9.3e\n",eqn,a0,a1,b0,b1);

  modelParameters(0,eqn)=a0;
  modelParameters(1,eqn)=a1;
  modelParameters(2,eqn)=b0;
  modelParameters(3,eqn)=b1;

  return 0;
}

// ==========================================================================================
/// \brief Set the parameter alphapP in the GDM model 
///
/// Generalized Dispersion Model:
///       E_tt - c^2 Delta(E) = -alphaP P_tt
///       Pi_tt + b1i Pi_1 + b0i = a0i*E + a1i*E_t    i=0,1,2,...,numPolarVectors-1
// ==========================================================================================
int DispersiveMaterialParameters::
setParameter( const real alphaP_ )
{
  alphaP = alphaP_;
  
  return 0;
}


// ==========================================================================================
/// \brief Set the parameter epsInf in the GDM model 
///
// ==========================================================================================
int DispersiveMaterialParameters::
setEpsInf( const real epsInf_ )
{
  epsInf=epsInf_;
  
  return 0;
}

// ==========================================================================================
/// \brief Set the parameter muInf in the GDM model 
///
// ==========================================================================================
int DispersiveMaterialParameters::
setMuInf( const real muInf_ )
{
  muInf=muInf_;
  
  return 0;
}



// ==========================================================================================
/// \brief Set the parameters in the GDM model for 1 polarization vector
/// \param a0,a1,b0,b1 (input)
/// 
/// Generalized Dispersion Model:
///       E_tt - c^2 Delta(E) = -alphaP P_tt
///       P_tt + b1 P_1 + b0 = a0*E + a1*E_t  
// ==========================================================================================
int DispersiveMaterialParameters::
setParameters( const real a0, const real a1, const real b0, const real b1 )
{
  numberOfPolarizationVectors=1;
  numberOfModelParameters=4;     
  modelParameters.redim(numberOfModelParameters,numberOfPolarizationVectors); 

  modelParameters(0,0)=a0;
  modelParameters(1,0)=a1;
  modelParameters(2,0)=b0;
  modelParameters(3,0)=b1;
  
  dbase.get<bool>("isDispersive")=true;

  return 0;
}

//=======================================================================================
// Multiply two matrices (2D A++ arrays) together
//   c <- a*b
//   
//======================================================================================
//=======================================================================================
// Multiply two matrices (2D A++ arrays) together
//   c <- a*b
//   
//======================================================================================
void DispersiveMaterialParameters::
myMatrixMultiply( RealArray & a, RealArray & b, RealArray & c ) 
{
  // ::display(a,"a");
  // ::display(b,"b");

  for( int i=a.getBase(0); i<=a.getBound(0); i++ )
    for( int j=a.getBase(1); j<=a.getBound(1); j++ )
    {
      real t=0;
      for( int k=a.getBase(1); k<=a.getBound(1); k++ )
        t=t+a(i,k)*b(k,j);
      c(i,j)=t;
    }
  
  // ::display(c,"c= a* b");

}

// ================================================================================
/// \brief Evaluate the eigenvalues of the BA Maxwell equation matrix
///
/// numberOfDimensions (input) number of space dimensions
/// reLambda (output) : maximum absolute value of the real part of all eigenvalues (max. wave speed)
/// imLambda (output) :  maximum absolute value of the imgainary part of all eigenvalues (non-zero implies
///                the problem is not hyperbolic (i.e. ill-posed)  
///
/// The Eigenvalues are found from the matrix (for all unit vectors (kx,ky,kz) )
///       A = K0^{-1} [ kx*Cx + ky*Cy + kz*Cz ]
///
///         = K0^{-1} [  0   0   0   0 -kz  ky]
///                   [  0   0   0  kz   0 -kx] 
///                   [  0   0   0 -ky  kx  0 ] 
///                   [  0  kz -ky   0   0  0 ] 
///                   [-kz   0  kx   0   0  0 ] 
///                   [ ky -kx   0   0   0  0 ] 
///        
// ================================================================================
int DispersiveMaterialParameters::
evalBianisotropicEigenValues( int numberOfDimensions, real & reLambda, real & imLambda ) 
{
  int debug=0;

  if( !dbase.has_key("eigenvaluesComputed") )
  {
    dbase.put<bool>("eigenvaluesComputed")=false;
    dbase.put<real>("reLambdaBA")=0.;
    dbase.put<real>("imLambdaBA")=0.;
  }


  bool & eigenvaluesComputed = dbase.get<bool>("eigenvaluesComputed");
  real & reLambdaBA = dbase.get<real>("reLambdaBA");
  real & imLambdaBA = dbase.get<real>("imLambdaBA");
  if( eigenvaluesComputed )
  {
    reLambda = reLambdaBA;
    imLambda = imLambdaBA;
    if( debug & 1 )
      printF("evalBianisotropicEigenValues: reLambda=%9.3e, imLambda=%9.3e (using saved values).\n",reLambda,imLambda);

    return 0;
  }
  
  eigenvaluesComputed=true;


  const int ex=0, ey=1, ez=2, hx=3, hy=4, hz=5;

  RealArray K0i;
  getBianisotropicMaterialMatrixInverse( K0i);
  
  RealArray A(6,6);
  RealArray C(6,6); C=0.;

  int n=6, lwork=4*n, info;
  RealArray wr(n), wi(n), work(lwork), vl(6,6), vr(6,6);
  real dummy;

  real lambdaMin=REAL_MAX, lambdaMax=-lambdaMin;
  real imPart=0.;
  // --- Loop over different directions to compute the eigenvalues ---
  int numx=5;
  int numy= numberOfDimensions<2 ? 1 : numx;
  int numz= numberOfDimensions<3 ? 1 : numx;

  for( int iz=0; iz<numz; iz++ )
  {
    for( int iy=0; iy<numy; iy++ )
    {
      for( int ix=0; ix<numx; ix++ )
      {
        if( ix==0 && iy==0 ) continue;
      
        real kx = real(ix)/max(1.,numx-1.);
        real ky = real(iy)/max(1.,numy-1.);
        real kz = real(iz)/max(1.,numz-1.);
        real kNorm= sqrt(kx*kx + ky*ky + kz*kz);
        kx/=kNorm; ky/=kNorm; kz/=kNorm;

        C(0,hz) = ky; C(0,hy)=-kz;    // (Hz)_y - (Hy)_z
        C(1,hx) = kz; C(1,hz)=-kx;    // (Hx)_z - (Hz)_x
        C(2,hy) = kx; C(2,hx)=-ky;    // (Hy)_x - (Hx)_y
  
        C(3,ez) =-ky; C(3,ey)= kz;    // -[ (Ez)_y - (Ey)_z ]
        C(4,ex) =-kz; C(4,ez)= kx;    // -[ (Ex)_z - (Ez)_x ]
        C(5,ey) =-kx; C(5,ex)= ky;    // -[ (Ey)_x - (Ex)_y ]
  

        myMatrixMultiply(K0i,C,A);

        // ::display(C,"C");
        // ::display(A,"A");

        // "N", "N"  = do not compute left and right eigenvectors 
        GEEV("N", "N", n, A(0,0), n, wr(0), wi(0), vl(0,0), n, vr(0,0), n, work(0), lwork, info);
        if( info !=0 )
        {
          printF("ERROR return info=%d from eigenvalue routine DGEEV\n",info);
        }
  
        real lamMin=min(wr), lamMax=max(wr);
        lambdaMin=min(lambdaMin,lamMin);
        lambdaMax=max(lambdaMax,lamMax);

        imPart = max(imPart,max(fabs(wi)));
      
        if( debug & 2 )
          for (int i=0; i<n; i++ )
            printF("kx=%5.3f ky=%5.3f kz=%5.3f lambda(%d)=(%9.2e,%9.2e)\n",kx,ky,kz, i,wr(i),wi(i));

        if( debug & 1 )
          printF("kx=%5.3f ky=%5.3f kz=%5.3f lambda: min=%9.3e max=%9.3e, imPart=%9.3e\n",kx,ky,kz,lamMin,lamMax,imPart);
        
      } // end for ix 
    } // end for iy 
  } // end for iz
  
  reLambda=max(fabs(lambdaMin),fabs(lambdaMax));
  imLambda=imPart;

  // Save the values so we do not need to recompute
  reLambdaBA = reLambda;
  imLambdaBA = imLambda;

  printF("\nlambda: min=%9.3e max=%9.3e (max-Im part=%8.2e -- zero=Hyperbolic)\n",lambdaMin,lambdaMax,imPart);
  if( fabs(imPart) > 100.*REAL_EPSILON*reLambda )
  {
    printF("DispersiveMaterialParameters: WARNING: Problem is ill-posed!\n");
    OV_ABORT("error");
  }
  
  
  
  return 0;
}

// ================================================================================
/// \brief Compute the time-stepping eigenvalue (reLambda,imLambda) for the BA Maxwell equations
///
/// \returnValue : dtMax 
///
/// The Eigenvalues are found from the matrix 
///       A = K0^{-1} [ Dx*Cx + Dy*Cy + Dz*Cz ]
///
///         = K0^{-1} [  0   0   0   0 -Dz  Dy]
///                   [  0   0   0  Dz   0 -Dx] 
///                   [  0   0   0 -Dy  Dx  0 ] 
///                   [  0  Dz -Dy   0   0  0 ] 
///                   [-Dz   0  Dx   0   0  0 ] 
///                   [ Dy -Dx   0   0   0  0 ] 
///        
/// where  Dx is the symbol of the first derivative operator
/// 
// ================================================================================
int DispersiveMaterialParameters::
evalBianisotropicTimeSteppingEigenvalue( MappedGrid & mg, const int orderOfAccuracy, real & reLambda, real & imLambda )
{
  int debug=0; // 3;

  if( !dbase.has_key("timeSteppingEigenvalueComputed") )
  {
    dbase.put<bool>("timeSteppingEigenvalueComputed")=false;
    dbase.put<real>("reLambdaTSE")=0.;
    dbase.put<real>("imLambdaTSE")=0.;
  }


  bool & timeSteppingEigenvalueComputed = dbase.get<bool>("timeSteppingEigenvalueComputed");
  real & reLambdaTSE = dbase.get<real>("reLambdaTSE");
  real & imLambdaTSE = dbase.get<real>("imLambdaTSE");

  if( timeSteppingEigenvalueComputed )
  {
    reLambda = reLambdaTSE;
    imLambda = imLambdaTSE;
    if( debug & 1 )
      printF("evalBianisotropicTimeSteppingEigenvalue: reLambda=%9.3e, imLambda=%9.3e (using saved values).\n",reLambda,imLambda);

    return 0;
  }
  
  timeSteppingEigenvalueComputed=true;

  const int numberOfDimensions = mg.numberOfDimensions();

  const bool isRectangular=mg.isRectangular();
  if( !isRectangular )
  {
    printF("evalBianisotropicTimeSteppingEigenvalue:ERROR: not implemented for curvilinear grids\n");
    OV_ABORT("fix me");
  }
 
 
  if( !( orderOfAccuracy==2 || orderOfAccuracy==4 ) )
  {
    printF("evalBianisotropicTimeSteppingEigenvalue:ERROR: not implemented for orderOfAccuracy=%d\n",orderOfAccuracy);
    OV_ABORT("fix me");
  }
 
  // Cartesian grid spacing: 
  real dx[3]={1.,1.,1.};
  mg. getDeltaX( dx );

  const int ex=0, ey=1, ez=2, hx=3, hy=4, hz=5;

  RealArray K0i;
  getBianisotropicMaterialMatrixInverse( K0i);
  
  RealArray A(6,6);
  RealArray C(6,6); C=0.;

  int n=6, lwork=4*n, info;
  RealArray wr(n), wi(n), work(lwork), vl(6,6), vr(6,6);
  real dummy;

  real reLambdaMin=REAL_MAX, reLambdaMax=-reLambdaMin;
  reLambda=0.; imLambda=0.;
  
  // --- Loop over different values of Dx, Dy and Dz to compute the time-stepping eigenvalues ---
  int numx=21;  //  5; // 21;  // number of wave numbers to check
  int numy= numberOfDimensions<2 ? 1 : numx;
  int numz= numberOfDimensions<3 ? 1 : numx;

  real xiMax[3]={0.,0.,0.}; // 
  real kMax[3]={0.,0.,0.}; // 
  
  for( int iz=0; iz<numz; iz++ )
  {
    for( int iy=0; iy<numy; iy++ )
    {
      for( int ix=0; ix<numx; ix++ )
      {
        real xix = Pi*ix/max(1.,numx-1.);  // 0 <= xix <= pi
        real xiy = Pi*iy/max(1.,numy-1.);
        real xiz = Pi*iz/max(1.,numz-1.);

	real kx,ky,kz;
        if( orderOfAccuracy==2 )
	{
	  // Symbols of D0
	  kx = sin(xix)/dx[0];
	  ky = sin(xiy)/dx[1];
	  kz = sin(xiz)/dx[2];
	}
	else if( orderOfAccuracy==4 )
	{
          real s2x= sin(xix/2.);
          real s2y= sin(xiy/2.);
          real s2z= sin(xiz/2.);
	  
          // Symbols of D0*( 1 - (1/6) D+D- )
	  kx = sin(xix)*( 1. + (2./3.)*SQR(s2x) )/dx[0];
	  ky = sin(xiy)*( 1. + (2./3.)*SQR(s2y) )/dx[1];
	  kz = sin(xiz)*( 1. + (2./3.)*SQR(s2z) )/dx[2];
	}
	
        C(0,hz) = ky; C(0,hy)=-kz;    // (Hz)_y - (Hy)_z
        C(1,hx) = kz; C(1,hz)=-kx;    // (Hx)_z - (Hz)_x
        C(2,hy) = kx; C(2,hx)=-ky;    // (Hy)_x - (Hx)_y
  
        C(3,ez) =-ky; C(3,ey)= kz;    // -[ (Ez)_y - (Ey)_z ]
        C(4,ex) =-kz; C(4,ez)= kx;    // -[ (Ex)_z - (Ez)_x ]
        C(5,ey) =-kx; C(5,ex)= ky;    // -[ (Ey)_x - (Ex)_y ]
  

        myMatrixMultiply(K0i,C,A);

        // ::display(C,"C");
        // ::display(A,"A");

        // "N", "N"  = do not compute left and right eigenvectors 
        GEEV("N", "N", n, A(0,0), n, wr(0), wi(0), vl(0,0), n, vr(0,0), n, work(0), lwork, info);
        if( info !=0 )
        {
          printF("ERROR return info=%d from eigenvalue routine DGEEV\n",info);
        }
  
	
	// we have left off an factor of "i" in the matrix so real-part of eigs is really imaginary part 
        real wrMax = max(fabs(wr));
	if( wrMax>imLambda)
	{ 
	  imLambda = wrMax;
	  xiMax[0]=xix; xiMax[1]=xiy; xiMax[2]=xiz;
	  kMax[0]=kx; kMax[1]=ky; kMax[2]=kz;
	}
	
        real wiMin = min(wi), wiMax=max(wi);
        reLambdaMin=min(reLambdaMin,wiMin);
        reLambdaMax=max(reLambdaMax,wiMax);

        reLambda = max(fabs(wiMin),fabs(wiMax));

        if( debug & 2 )
          for (int i=0; i<n; i++ )
            printF("BA-TSE: xix=%5.3f, xiy=%5.3f, xiz=%5.3f, kx=%7.2f ky=%7.2f kz=%7.2f lambda(%d)=(%9.2e,%9.2e)\n",
		   xix,xiy,xiz,kx,ky,kz, i,wi(i),wr(i));

        
      } // end for ix 
    } // end for iy 
  } // end for iz
  
  // Save the values so we do not need to recompute
  reLambdaTSE = reLambda;
  imLambdaTSE = imLambda;

  printF("\n evalBianisotropicTimeSteppingEigenvalue real-lambda: min=%9.3e max=%9.3e (max-Re part=%8.2e => zero=Hyperbolic)\n",
	 reLambdaMin,reLambdaMax,reLambda);

  if( fabs(reLambda) > 100.*REAL_EPSILON*imLambda )
  {
    printF("DispersiveMaterialParameters:evalBianisotropicTimeSteppingEigenvalue: WARNING: Problem is ill-posed!\n");
    OV_ABORT("error");
  }
  
  printF(" imLambda=%9.3e occured at xi=[%g,%g,%g] DHat=[%g,%g,%g]\n\n",imLambda,xiMax[0],xiMax[1],xiMax[2],kMax[0],kMax[1],kMax[2]);
  
  
  return 0;
}


// ==========================================================================================
/// \brief Solve the dispersion relation for k given s, but adjust s so that the resulting
///   solution gives a value for k that is two-pi periodic in space. This only works if
///   we have a root with Im(k)=0 
///
/// \param kv[3] (input) : wave vector direction, entries assumed to be multiples of 2*pi
/// \param sr,si (input/output) : On input, a guess for s, on output the value for s that gives
///             a two-pi periodic solution in space
/// \param kr,ki (output) : complex wave vector coefficient
/// \param evr[6], evi[6], chi (output) : eigenvector and susceptibilities.
// ==========================================================================================
int DispersiveMaterialParameters::
findPeriodicSolution( real kv[3], real & sr, real & si, real & kr, real & ki, real evr[6], real evi[6], RealArray & chi )
{
  DispersionRelationOptionEnum computeOption =computeComplexWaveNumber;
  setDispersionRelationComputeOption( computeOption );

  // Guess:
  // sr=0.;
  // si = -sqrt( SQR(kv[0]) + SQR(kv[1]) + SQR(kv[2]) );
  
  int & debug = dbase.get<int>("debug");
  debug=1;
  
  printF("+++ findPeriodicSolution: kv=[%e,%e,%e]\n",kv[0],kv[1],kv[2]);

  // --- find the target wave number ---

  real signS = si >=0. ? 1. : -1.;

  // (1) get an initial guess for kr 
  real skr=sr, ski=si;
  getBianisotropicDispersivePlaneWaveSolution( kv, skr,ski,evr,evi,chi,noPolarization );
  kr=skr; ki=ski;
  printF("START: s= %e + %e I --> k = %e + %e I\n",sr,si,kr,ki );
  
  debug=0;

  // Target kr has an integer value 
  real krTarget = floor( kr + .5 );
  if( krTarget>0. ){ krTarget=max(1.,krTarget); } else{  krTarget=min(krTarget,-1.);}  // 
  assert( fabs(krTarget) > REAL_EPSILON*10. );

  printF("+++ findPeriodicSolution: target kr=%e\n",krTarget);

  // Find a interval [sia,sib] or [kra,krb] which holds the target
  real sia, sib, kra, krb, ds;
  if( kr >= krTarget )
  { // we have an upper bound, search backward 
    sib = si; krb = kr;  ds =-.2;
  }
  else
  { // we have a lower bound, search forward 
    sia = si; kra = kr; ds =+.2;
  }

  int maxCheck=101;
  bool found=false;
  for( int it=0; it<maxCheck; it++ )
  {
    si += ds*signS;

    skr=sr, ski=si;
    getBianisotropicDispersivePlaneWaveSolution( kv, skr,ski,evr,evi,chi,noPolarization );
    kr=skr; ki=ski;

    printF("it=%i: s= %e + %e I --> k = %e + %e I\n",it,sr,si,kr,ki );

    if( ds>0. )
    {
      if( kr >= krTarget )
      {
	sib=si;  krb=kr; found=true; // upper bound found 
	break;
      }
      else
      {
        sia=si;  kra = kr; // new lower bound 
      }
    }

    if( ds<0. )
    {
      if( kr <= krTarget )
      {
	sia=si; kra=kr; found=true;  // lower bound found 
        break;
      }
      else
      {
	sib = si; krb = kr;  // new upper bound 
      }
    }
  }
  assert( found );

  printF("+++ findPeriodicSolution: Range [sia,sib]=[%e,%e] [kra,krb]=[%e,%e]\n",sia,sib,kra,krb);

  // -- Now search for the root ---
  const real tol = REAL_EPSILON*100;
  int maxIt=100; found=false;
  for( int it=0; it<maxIt; it++ )
  {
    real sim = .5*( sia+sib );
    skr=sr, ski=sim;
    getBianisotropicDispersivePlaneWaveSolution( kv, skr,ski,evr,evi,chi,noPolarization );
    kr=skr; ki=ski;

    if( kr >= krTarget )
    {
      sib = sim; krb = kr;  // new upper bound 
    }
    else 
    {
      sia = sim;  kra = kr; // new lower bound 
    }
    printF("bisect: it=%d: Range [sia,sib]=[%20.16e,%20.16e] [kra,krb]=[%20.16e,%20.16e]\n",it,sia,sib,kra,krb);
    
    if( fabs( sib-sia ) <= tol )
    {
      si=sim; found=true;
      break;
    }
  }
  if( !found )
  {
    printF("+++ findPeriodicSolution: WARNING: solution not found to required tolerance\n");
  }
  
  printF("+++ Solution: s= %e + %e I --> k = %20.16e + %e I, |kr-krTarget|/|krTarget| =%9.2e\n",sr,si,kr,ki,
	 fabs(kr-krTarget)/fabs(krTarget) );
  
  const real tolki = tol*1000.;
  if( fabs(ki)>tolki )
  {
    printF("+++ findPeriodicSolution: WARNING: This solution is NOT periodic since Im(k) =%e != 0 (tol=%8.2e)\n",ki,tolki);
    printF("You could try using a different solution mode.\n");
    OV_ABORT("error");
  }
  

  return 0;
  
}

#include <complex>

// ==========================================================================================
/// \brief Plot dispersive properties versus frequency or wavelength
// ==========================================================================================
int DispersiveMaterialParameters::
update( GenericGraphicsInterface & gi )
{

  DispersiveMaterialParameters & dmp = *this;
  DispersionRelationOptionEnum & computeOption = dbase.get<DispersionRelationOptionEnum>("dispersionRelationComputeOption");

  // ------- Plot some quantities --------
  const Real PHz = 1e15;
  const Real nm = 1.e-9;

  const Real & V0 = dbase.get<Real>("velocityScale");  // velocity scale 
  const Real & L0 = dbase.get<Real>("lengthScale");    // length scale 

  Real omegaMin = dbase.get<Real>("omegaMin");   
  Real omegaMax = dbase.get<Real>("omegaMax");   

  const MaterialTypeEnum & materialType = dbase.get<MaterialTypeEnum>("materialType");

  int iBA=1, jBA=1; // Entry in BA tensor K to plot 

  Real kv[3] = { 1./sqrt(2.),1./sqrt(2.),0.}; //  wave vector for plane wave solution
  Real sr =0 ,si=-1.;

  int Np;
  Real epsInf;
  RealArray modelParams;
  
  dmp.getIsotropicParameters( Np,  epsInf, modelParams );
  printF("Np=%d, epsInf=%g\n",Np,epsInf);
  for( int j=0; j<Np; j++ )
    printF("j=%d: [a0,a1,b0,b1]=[%g,%g,%g,%g]\n",j,modelParams(0,j),modelParams(1,j),modelParams(2,j),modelParams(3,j));



  int nw=101;

  const Real & omegaScale = dbase.get<Real>("omegaScale");
  // Real OmegaScale = 9.4182578365442600e+15; // From Gold approx 10 PHz  *** FIX ME ***

  // Real omega0 = 2.*Pi*V0/L0;
  Real omega0 = V0/L0;
  // printF(" L0=%9.2e (nm), Omega0 = 2.*Pi*V0/L = %9.3e, OmegaScale=%9.3e (from file)\n",L0/nm, omega0,omegaScale);
  printF(" L0=%9.2e (nm), Omega0 = V0/L = %9.3e, OmegaScale=%9.3e (from file)\n",L0/nm, omega0,omegaScale);

  // By default plot over the range of omega for which the fit was made:
  // Real wMin=.15, wMax=1.1;
  RealArray w(nw), epsHatr(nw), epsHati(nw), nHatr(nw), nHati(nw), om(nw), lam(nw);

  // -- Array to hold results to plot --
  int numFields = materialType==isotropic ? 4 : 2;
  RealArray fields(nw,numFields);


  PlotStuffParameters psp;
  GUIState dialog;

  dialog.setWindowTitle("DispersiveMaterialParameters");
  dialog.setExitCommand("exit", "exit");

  aString cmds[] = {"plot versus omega",
                    "plot versus lambda",
                    "display parameters",
                    "bianisotropic eigenvalues",
                    "bianisotropic plane wave s",
                    "bianisotropic plane wave k",
                    "find periodic solution",
                    // "find BA plane wave solution",
		    ""};

  int numberOfPushButtons=0;  // number of entries in cmds
  while( cmds[numberOfPushButtons]!="" ){numberOfPushButtons++;}; // 
  // int numRows=(numberOfPushButtons+1)/2;
  int numRows=numberOfPushButtons;
  dialog.setPushButtons( cmds, cmds, numRows ); 

  aString opCommand1[] = {"compute s given k",
         		  "compute k given s",
        		  ""};

  dialog.setOptionMenuColumns(1);
  dialog.addOptionMenu( "Dispersion relation:", opCommand1, opCommand1, (int) computeOption );

  const int numberOfTextStrings=15;  // max number allowed
  aString textLabels[numberOfTextStrings];
  aString textStrings[numberOfTextStrings];

  int nt=0;

  textLabels[nt] = "omegaMin:"; 
  sPrintF(textStrings[nt],"%g",omegaMin);  nt++; 

  textLabels[nt] = "omegaMax:"; 
  sPrintF(textStrings[nt],"%g",omegaMax);  nt++; 

  textLabels[nt] = "BA tensor entry:"; 
  sPrintF(textStrings[nt],"%i, %i",iBA,jBA);  nt++; 

  textLabels[nt] = "kv:"; 
  sPrintF(textStrings[nt],"%g, %g, %g (for plane wave solution)",kv[0],kv[1],kv[2]);  nt++; 

  textLabels[nt] = "s:"; 
  sPrintF(textStrings[nt],"%g, %g (for plane wave solution)",sr,si);  nt++; 

  //  textLabels[nt] = "output file:"; 
  // sPrintF(textStrings[nt],"%s",(const char*)outputFileName);  nt++; 

  // null strings terminal list
  textLabels[nt]="";   textStrings[nt]="";  assert( nt<numberOfTextStrings );
  dialog.setTextBoxes(textLabels, textLabels, textStrings);


  gi.pushGUI(dialog);

  bool recompute=true;
  aString answer;
  int len;
  
  for(;;)
  {
    
    gi.getAnswer(answer,"");  

    if( answer=="exit" )
    {
      break;
    }
    else if( dialog.getTextValue(answer,"omegaMin:","%e",omegaMin) ){recompute=true;} //
    else if( dialog.getTextValue(answer,"omegaMax:","%e",omegaMax) ){recompute=true;} //

    else if( len=answer.matches("BA tensor entry:") )
    {
      sScanF(answer(len,answer.length()-1),"%i %i",&iBA,&jBA);
      if( iBA<1 || iBA>3 )
      {
	printF("ERROR: iBA =%i must be between 1 and 3\n",iBA); iBA=max(1,min(iBA,3));
      }
      if( jBA<1 || jBA>3 )
      {
	printF("ERROR: jBA =%i must be between 1 and 3\n",jBA); jBA=max(1,min(jBA,3));
      }
      
      printF("Setting tensor entry to K(%i,%i)\n",iBA,jBA);
    }
    

    else if( answer=="compute s given k" ||
             answer=="compute k given s" )
    {
      if( answer=="compute s given k" )
      {
	computeOption=computeComplexFrequency;
	printF("Compute s given k when solving the dispersion relation.\n");
      }
      else
      {
	computeOption=computeComplexWaveNumber;
	printF("Compute k given s when solving the dispersion relation.\n");
      }
    }

    else if( answer=="display parameters" )
    {
      display();
    }
    
    else if( answer=="plot versus omega" || answer=="plot versus lambda" )
    {

      if( recompute )
      {
        Real wMin=omegaMin, wMax=omegaMax;  
        for( int i=0; i<nw; i++ )
        {
          w(i) = wMin + (wMax-wMin)*(i)/Real(nw-1);

          Real omega = w(i)*omegaScale;          // omega in MKS units (Hz)
          Real lambda = (2.*Pi*V0/(omega)) / nm; // lambda in nm 

          w(i) *= omegaScale/omega0;  // non-dimensional omega

          om(i)=omega/PHz;
          lam(i)=lambda;

	  if( materialType==isotropic )
	  {
	    dmp.evalEpsAndN( w(i), epsHatr(i), epsHati(i), nHatr(i), nHati(i) );

	    printF(" w=%7.3g (= %8.2e PHz) lambda=%8.2e (nm) epsHat=(%9.2e,%9.2e) n=(%9.2e,%9.2e)\n",
                   w(i),omega/PHz,lambda,epsHatr(i), epsHati(i), nHatr(i), nHati(i) );
	  }
	  else
	  {
	    dmp.evalMaterialTensor( w(i), epsHatr(i), epsHati(i), iBA,jBA );

	    printF(" w=%7.3g (= %8.2e PHz) lambda=%8.2e (nm) epsHat=(%9.2e,%9.2e)\n",
                   w(i),omega/PHz,lambda,epsHatr(i), epsHati(i) );
	  }
	  

    

        }

        Range all;
        fields(all,0) = epsHatr;
        fields(all,1) = epsHati;
	if( materialType==isotropic )
	{
	  fields(all,2) = nHatr;
	  fields(all,3) = nHati;
	}
	
        recompute=false;
      }
      
      // psp.set(GI_TOP_LABEL,sPrintF("Test Read Material File"));
      aString title;
      if( materialType==isotropic )
        sPrintF(title,"DMP: %s, Np=%i",(const char*)dbase.get<aString>("name"),numberOfPolarizationVectors );
      else
      {
        const IntegerArray & NpBA = dbase.get<IntegerArray>("NpBA");

        sPrintF(title,"DMP: %s, bianisotropic K(%i,%i) Np=%i",(const char*)dbase.get<aString>("name"),iBA,jBA,NpBA(iBA-1,jBA-1) );
      }
      
      const aString names[]={"epsr","epsi","nr","ni"};

      #ifndef USE_PPP
      if( answer=="plot versus omega"  )
        PlotIt::plot(gi, om, fields, title,"omega (PHz)", names,psp );
      else
        PlotIt::plot(gi, lam, fields, title,"lambda (nm)", names,psp );
      #else
      printF("DMP: FINISH PLOTTING IN PARALLEL\n");
      #endif

    }
    else if( answer=="bianisotropic eigenvalues" )
    {
      Real reLambda, imLambda;
      int nd=2;  // 2D 
      evalBianisotropicEigenValues( nd, reLambda, imLambda);
    }
    else if( len=answer.matches("kv:") )
    {
      sScanF(answer(len,answer.length()-1),"%e %e %e",&kv[0],&kv[1],&kv[2]);
      printF("Setting kv=[%e,%e,%e]\n",kv[0],kv[1],kv[2]);
    }
    else if( len=answer.matches("s:") )
    {
      sScanF(answer(len,answer.length()-1),"%e %e",&sr,&si);
      printF("Setting s=%e + %e I]\n",sr,si);
    }
    else if( answer=="bianisotropic plane wave s" ||
             answer=="bianisotropic plane wave k" )
    {
      Real evr[6], evi[6];  // eigenvector 
      RealArray chi;
      DispersionRelationOptionEnum computeOption =
          answer=="bianisotropic plane wave s" ? computeComplexFrequency : computeComplexWaveNumber;

      setDispersionRelationComputeOption( computeOption );
      getBianisotropicDispersivePlaneWaveSolution( kv, sr,si,evr,evi,chi,noPolarization );

    }
    else if( answer=="find periodic solution" )
    {
      //  *** Put this here for now for testing *****

      Real evr[6], evi[6];  // eigenvector 
      RealArray chi;

      Real kx =1, ky=1, kz=0;
      Real kvp[3]= { kx*twoPi, ky*twoPi, kz*twoPi}; // 

      Real sr=0.;
      Real si = -sqrt( SQR(kvp[0]) + SQR(kvp[1]) + SQR(kvp[2]) );

      Real kr,ki;
      findPeriodicSolution( kvp, sr,si,kr,ki,evr,evi,chi );

    }
    
/* -- Moved to texact.C 
    else if( answer=="find BA plane wave solution" )
    {
      Real evr[6], evi[6];  // eigenvector 
      RealArray chi, qv;

      Real kx =1, ky=1, kz=0;
      Real kvp[3]= { kx*twoPi, ky*twoPi, kz*twoPi}; // 

      Real sr=0.;
      Real si = -sqrt( SQR(kvp[0]) + SQR(kvp[1]) + SQR(kvp[2]) );

      // DispersionRelationOptionEnum computeOption =computeComplexWaveNumber;
      DispersionRelationOptionEnum computeOption =computeComplexFrequency;
      setDispersionRelationComputeOption( computeOption );

      DispersiveMaterialParameters & dmp2 = *this;  // do this for now 

      Real skr=sr, ski=si;
      
      // typedef std::vector<std::complex<Real> > ComplexVector;
      // ComplexVector qvi;
      // getBianisotropicPlaneInterfaceSolution(  dmp2, kv, skr, ski, qvi,qv,  chi, noPolarization );
      PlaneInterfaceExactSolution pies;
      
      pies. getBianisotropicPlaneInterfaceSolution( *this, dmp2, kv, skr, ski, qv,  chi );
    }
    --- */
    
  }
  

  gi.popGUI(); // restore the previous GUI  


  return 0;
}





// ------ routines use COMPLEX arithmetic go at the bottom to avoid name conflicts with "real" ---------------

#include <complex>
typedef ::real LocalReal;
typedef ::real OvReal;

#include "ComplexArray.h"


  

// ===============================================================================================================================
/// \brief Return the frequency space BA material matrix evaluated at "s"
// ===============================================================================================================================
int DispersiveMaterialParameters::
getK( ComplexArray & K, Complex & s ) const
{
  printF("DMP:: ENTERING getK\n");
  if( !dbase.has_key("K0") )
  {
    // --- isotropic case ---
    printF("DMP::getK: there is no K0 matrix!\n");
    for( int i1=0; i1<6; i1++ )
    {
      for( int i2=0; i2<6; i2++ )
      {
	K(i1,i2) = 0.;
      }
    }
    for( int i=0; i<3; i++ ){ K(i,i)=epsInf; } // 
    for( int i=3; i<6; i++ ){ K(i,i)=muInf; } // 

    if( isDispersiveMaterial() )
    {
      OV_ABORT("DMP::getK: finish me for dispersive");
    }
    
    return 0;

  }
  
  RealArray & K0 = dbase.get<RealArray>("K0");

  for( int i1=0; i1<6; i1++ )
  {
    for( int i2=0; i2<6; i2++ )
    {
      K(i1,i2) = K0(i1,i2);
    }
  }
  if( !isDispersiveMaterial() )
  {
    return 0;
  }
  

  // bianisotropicParameters(4,Np,6,6) 
  if(!dbase.has_key("bianisotropicParameters") )
  {
   printF("DMP::getK: there is no array bianisotropicParameters!\n");
  }
  

  RealArray & bianisotropicParameters = dbase.get<RealArray>("bianisotropicParameters");
  IntegerArray & NpBA = dbase.get<IntegerArray>("NpBA");

  for( int i1=0; i1<6; i1++ )
  {
    for( int i2=0; i2<6; i2++ )
    {
      // K(i1,i2) = K0(i1,i2);
      for( int n=0; n<NpBA(i1,i2); n++ )
      {
	OvReal a0=bianisotropicParameters(0,n,i1,i2);
	OvReal a1=bianisotropicParameters(1,n,i1,i2);
	OvReal b0=bianisotropicParameters(2,n,i1,i2);
	OvReal b1=bianisotropicParameters(3,n,i1,i2);
	  
	K(i1,i2) = K(i1,i2) + (a0+a1*s)/( b0+b1*s+s*s);
      }
    }
  }


  return 0;
}

