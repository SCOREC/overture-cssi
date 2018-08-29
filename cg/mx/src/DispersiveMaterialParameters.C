// Generalized Dispersion Model:
//       E_tt - c^2 Delta(E) = -alphaP P_tt
//       P_tt + b1 P_1 + b0 = a0*E + a1*E_t 
// 
#include "DispersiveMaterialParameters.h"

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



// ============================================================================
/// \brief Class to define parameters of a dispersive material.
///
/// Generalized Dispersion Model:
///       E_tt - c^2 Delta(E) = -(1/eps) P_tt
///       P_tt + b1 P_1 + b0 = eps*( a0*E + a1*E_t )
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
  
  numberOfPolarizationVectors=0; // by default a domain is non-dispersive
  numberOfModelParameters=4;     // [a_k,b_k,c_k,d_k] 
  modelParameters.redim(numberOfModelParameters,1); // fill in defaults of zero
  modelParameters=0.;

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
  
}

// ============================================================================
/// \brief Copy constructor
// ============================================================================
DispersiveMaterialParameters::
DispersiveMaterialParameters(const DispersiveMaterialParameters& x)
{
  *this=x;  
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

  alphaP=x.alphaP;
  mode  =x.mode;
  rootComputed=x.rootComputed;
  ck0=x.ck0;
  
  gamma=x.gamma;
  omegap=x.omegap;

  return *this;
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
    printF("--DMP-- GDM: RECOMPUTE root: rootComputed=%i c*k=%e ck0=%e\n",(int)rootComputed,c*k,ck0);
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
/// \brief return the parameter alphaP -- normally = 1/eps if set, -1 if not set.
// ==========================================================================================
real DispersiveMaterialParameters::
getAlphaP() const
{
  return alphaP;
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
/// \brief Specify the number of polarization vectors (number of GDM equation)
// ==========================================================================================
int DispersiveMaterialParameters::
setParameter( const real alphaP_ )
{
  alphaP = alphaP_;
  
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
  
  return 0;
}

