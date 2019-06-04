// Generalized Dispersion Model:
//       E_tt - c^2 Delta(E) = -alphaP P_tt
//       P_tt + b1 P_1 + b0 = a0*E + a1*E_t 
// 
#include "DispersiveMaterialParameters.h"
#include "PlotStuff.h"
#include "display.h"

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
  epsInf=1;    // eps(s=infinity)
  
  numberOfPolarizationVectors=0; // by default a domain is non-dispersive
  numberOfModelParameters=4;     // [a0,a1,b0,b1] 
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

}

// ============================================================================
/// \brief Copy constructor
// ============================================================================
DispersiveMaterialParameters::
DispersiveMaterialParameters(const DispersiveMaterialParameters& x)
{
  *this=x;  
  if( !dbase.has_key("name" ) ) dbase.put<aString>("name");
  dbase.get<aString>("name") = x.dbase.get<aString>("name");

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
  alphaP=x.alphaP;
  mode  =x.mode;
  rootComputed=x.rootComputed;
  ck0=x.ck0;
  
  gamma=x.gamma;
  omegap=x.omegap;

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

  return *this;
}


// ==========================================================================================
/// \brief Return the material name 
// ==========================================================================================
aString DispersiveMaterialParameters::getMaterialName() const
{
  return dbase.get<aString>("name");
}


// ==========================================================================================
/// \brief Display parameters
// ==========================================================================================
int DispersiveMaterialParameters::
display( FILE *file /* = stdout */ ) const
{
  const aString & name = dbase.get<aString>("name");

  MaterialTypeEnum & materialType = dbase.get<MaterialTypeEnum>("materialType");
  const real & omegaMin = dbase.get<real>("omegaMin");   
  const real & omegaMax = dbase.get<real>("omegaMax");   
  
  const real V0 = dbase.get<real>("velocityScale");
  const real L0 = dbase.get<real>("lengthScale"); 
  const real & omegaScale = dbase.get<real>("omegaScale");

  fPrintF(file,"----------------------------------------------------------------------------\n");
  fPrintF(file,"-------------------- Dispersive Material Parameters ------------------------\n");
  
  fPrintF(file," name=[%s]\n",(const char*)name);
  fPrintF(file," materialType=%s\n",(materialType==isotropic ? "isotropic" :
                               materialType==bianisotropic ? "bianisotropic" : "unknown"));
  fPrintF(file," Length-scale L0=%9.3e, velocity-scale V0=%9.3e, omegaScale=%9.3e \n",L0,V0,omegaScale);
  
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
/// \brief return the parameter alphaP -- normally = 1/eps if set, -1 if not set.
// ==========================================================================================
real DispersiveMaterialParameters::
getAlphaP() const
{
  return alphaP;
}

/* -----
// ==========================================================================================
/// \brief interactive update: define a dispersive material or read from a data file.
// ==========================================================================================
int DispersiveMaterialParameters::
update()
{
  GenericGraphicsInterface & gi = *gip;

  aString fileName="myMaterial.txt";

  GUIState gui;

  DialogData & dialog=gui;

  dialog.setWindowTitle("Dispersive Material Parameters");
  dialog.setExitCommand("exit", "exit");

  dialog.setOptionMenuColumns(1);

  // ************** PUSH BUTTONS *****************
  aString pushButtonCommands[] = {"read file",
				  ""};
  int numRows=1;
  dialog.setPushButtons(pushButtonCommands,  pushButtonCommands, numRows ); 


  // ----- Text strings ------
  const int numberOfTextStrings=30;
  aString textCommands[numberOfTextStrings];
  aString textLabels[numberOfTextStrings];
  aString textStrings[numberOfTextStrings];

  int nt=0;

  textCommands[nt] = "file name:";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%s",fileName);  nt++; 


  textCommands[nt] = "numberOfPolarizationVectors:";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%i",numberOfPolarizationVectors);  nt++; 


  // null strings terminal list
  assert( nt<numberOfTextStrings );
  textCommands[nt]="";   textLabels[nt]="";   textStrings[nt]="";  
  dialog.setTextBoxes(textCommands, textLabels, textStrings);

  gi.pushGUI(gui);
  aString answer,line;
  int len=0;
  for(;;) 
  {
    gi.getAnswer(answer,"");      
    // printF("Start: answer=[%s]\n",(const char*) answer);
    
    if( answer=="continue" || answer=="exit" )
    {
      break;
    }
    else if( answer=="time stepping options..." )
    {
      timeSteppingOptionsDialog.showSibling();
    }


    else
    {
      printF("Unknown command = [%s]\n",(const char*)answer);
      gi.stopReadingCommandFile();
       
    }
  }

  gi.popGUI();  // pop dialog
  
  return 0;
}

---- */

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
/// \brief Read dispersive material parameters from a file.
///
/// \param fileName (input) : name of material file
/// \param numberOfPolarizationVectorsRequested (input): optionally specify the number of polarization vectors
 //   if more than one fit is available from the file (-1 = choose the largest number).
// ==========================================================================================
int DispersiveMaterialParameters::
readFromFile( const aString & fileName, int numberOfPolarizationVectorsRequested /* = -1 */ )
{
  if( numberOfPolarizationVectorsRequested <=0 )
    numberOfPolarizationVectorsRequested=INT_MAX;  // choose maximum number available

  FILE *file = fopen((const char*)fileName, "r");
  if( file==NULL )
  {
    printF("DispersiveMaterialParameters::readFromFile:ERROR: unable to open fileName=[%s]\n",(const char*)fileName);
    OV_ABORT("error");
  }

  MaterialTypeEnum & materialType = dbase.get<MaterialTypeEnum>("materialType");
  if( materialType == bianisotropic )
  {
     //    K(6,6) : material tensor 
     //
     //     K0(6,6)  : constant part of material tensor 
     //     bianisotropicParameters(4,Np,6,6)        : GDM 
     //     Np(6,6) : number of polarization vectors 

    if( ! dbase.has_key("bianisotropicParameters" ) )
    {
      dbase.put<IntegerArray>("NpBA");
      dbase.put<RealArray>("K0");
      dbase.put<RealArray>("bianisotropicParameters");
     
      RealArray & bianisotropicParameters = dbase.get<RealArray>("bianisotropicParameters");
      int NpMax=10;   // max Np **FIX ME**
      bianisotropicParameters.redim(4,NpMax,6,6);

      bianisotropicParameters=0.;
	
      IntegerArray & NpBA = dbase.get<IntegerArray>("NpBA");
      NpBA.redim(6,6); NpBA=0;
      
    }
    
  }
  

  aString & name = dbase.get<aString>("name");
  real & omegaMin = dbase.get<real>("omegaMin");   
  real & omegaMax = dbase.get<real>("omegaMax");   
  
  aString units;
  
  real & omegaScale = dbase.get<real>("omegaScale");
  RealArray gdmPars;
  real epsInf0;
  int npMax=10;
  RealArray a0v(npMax),a1v(npMax),b0v(npMax),b1v(npMax);


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
  real Omega0 = 1./T0;       // angular frequency scales by this if we scale x/LO0 and t/T0 

  real wScale=1.; // Frequency scale -- set below 

  int currentBestNp=0;
  
  int len;
  const int buffLength=2000;
  char line[buffLength];
  int numRead=getLineFromFile(file,line,buffLength);  // read a line from the file.
  while( numRead>0 )
  { 
    aString aline=line;
    // printF("line=[%s]\n",(const char*)line);

    if( line[0]=='#' )
      printF("Comment: %s\n",(const char*)line);
    else
    {
      if( len=aline.matches("omegaScale=") )
      {
	sScanF(aline(len,aline.length()-1),"%e",&omegaScale);
        printF("omegaScale=%22.16e\n",omegaScale);

	wScale=omegaScale/Omega0; // Frequency scale 
	printF("DispersiveMaterialParameters::readFromFile: non-dimensionalize: L0=%g(nm) = %g(m) T0=%g(s) Omega0=%g omegaScale=%g,  "
	       "wscale= omegaScale*T0/(2*pi)=%g\n", L0/nm,L0,T0,Omega0,omegaScale,wScale);

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
	name=aline(len,aline.length()-2); // skip final ";"
        printF("name=%s\n",(const char*)name);
        int nl=name.length()-1;
        if( name[0]=='"' && name[nl]=='"' )
	  name = name(1,nl-1);  // remove double quotes 
      }
      else if( len=aline.matches("units=") )
      {
	units=aline(len,aline.length()-2); // skip final ";"
        printF("units=%s\n",(const char*)units);
      }
      else if( len=aline.matches("epsGDMPars") )	
      {
        // "epsGDMPars1=["
        // "epsGDMPars2=["
        // "epsGDMPars3=["

        // GDM parameters
	int Np=1;
        sScanF(aline(len,len+1),"%i",&Np);
	if( Np==1 )
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
        
	if( Np <= numberOfPolarizationVectorsRequested && Np>currentBestNp )
	{
	  currentBestNp=Np;
	  
	  epsInf=epsInf0;
	  alphaP=1./epsInf;

	  numberOfPolarizationVectors=Np;
	  printF("\n DispersiveMaterialParameters: SETTING numberOfPolarizationVectors=%d\n\n",numberOfPolarizationVectors);

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


      else if( materialType == bianisotropic && (len=aline.matches("bianistropicK0"))   )	
      {
	// ---- BI-ANISTROPIC MATERIAL MATRIX K0 ---
        //   bianistropicK0=[...]; 
        // K0 matrix
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
	
	::display(K0,"K0");
      }
      
      else if( materialType == bianisotropic && (len=aline.matches("bianistropicPars")) )	
      {
	// ---- BI-ANISTROPIC MATERIAL GDM COEFFICIENTS ---

        // Entry form: 
        // "bianistropicPars[1-6][1-6]GDM[1-]=["

        // bianisotropicParameters(4,Np,6,6) 
        RealArray & bianisotropicParameters = dbase.get<RealArray>("bianisotropicParameters");
        IntegerArray & NpBA = dbase.get<IntegerArray>("NpBA");

        // GDM parameters
	int k1=1,k2=1;
        sScanF(aline(len  ,len  ),"%i1",&k1);
        sScanF(aline(len+1,len+1),"%i1",&k2);
        // printF("Reading GDM parameters for material matrix K%d%d\n",k1,k2);

        assert( k1>=1 && k1<=6 );
        assert( k2>=1 && k2<=6 );
	int Np; 
        // printf("aline(len+5,len+5)=[%s]\n",(const char*)aline(len+5,len+5));
	
	sScanF(aline(len+5,len+5),"%i1",&Np);
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
	

      }


      else
      {
        printF("Ignored: line=[%s]\n",(const char*)line);      
      }
      

    }
    



    numRead=getLineFromFile(file,line,buffLength);
  }
     

  fclose(file);

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


// ==========================================================================================
/// \brief Plot dispersive properties versus frequency or wavelength
// ==========================================================================================
int DispersiveMaterialParameters::
update( GenericGraphicsInterface & gi )
{

  DispersiveMaterialParameters & dmp = *this;

  // ------- Plot some quantities --------
  const real PHz = 1e15;
  const real nm = 1.e-9;

  const real & V0 = dbase.get<real>("velocityScale");  // velocity scale 
  const real & L0 = dbase.get<real>("lengthScale");    // length scale 

  real omegaMin = dbase.get<real>("omegaMin");   
  real omegaMax = dbase.get<real>("omegaMax");   

  const MaterialTypeEnum & materialType = dbase.get<MaterialTypeEnum>("materialType");

  int iBA=1, jBA=1; // Entry in BA tensor K to plot 

   

  int Np;
  real epsInf;
  RealArray modelParams;
  
  dmp.getIsotropicParameters( Np,  epsInf, modelParams );
  printF("Np=%d, epsInf=%g\n",Np,epsInf);
  for( int j=0; j<Np; j++ )
    printF("j=%d: [a0,a1,b0,b1]=[%g,%g,%g,%g]\n",j,modelParams(0,j),modelParams(1,j),modelParams(2,j),modelParams(3,j));



  int nw=101;

  const real & omegaScale = dbase.get<real>("omegaScale");
  // real OmegaScale = 9.4182578365442600e+15; // From Gold approx 10 PHz  *** FIX ME ***

  // real omega0 = 2.*Pi*V0/L0;
  real omega0 = V0/L0;
  // printF(" L0=%9.2e (nm), Omega0 = 2.*Pi*V0/L = %9.3e, OmegaScale=%9.3e (from file)\n",L0/nm, omega0,omegaScale);
  printF(" L0=%9.2e (nm), Omega0 = V0/L = %9.3e, OmegaScale=%9.3e (from file)\n",L0/nm, omega0,omegaScale);

  // By default plot over the range of omega for which the fit was made:
  // real wMin=.15, wMax=1.1;
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
		    ""};

  int numberOfPushButtons=0;  // number of entries in cmds
  while( cmds[numberOfPushButtons]!="" ){numberOfPushButtons++;}; // 
  // int numRows=(numberOfPushButtons+1)/2;
  int numRows=numberOfPushButtons;
  dialog.setPushButtons( cmds, cmds, numRows ); 

  // aString opCommand1[] = {"omega",
  //       		  "lambda",
  //       		  ""};

  // dialog.setOptionMenuColumns(1);
  // dialog.addOptionMenu( "x-axis:", opCommand1, opCommand1, xAxisType );

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
    

    else if( answer=="display parameters" )
    {
      display();
    }
    
    else if( answer=="plot versus omega" || answer=="plot versus lambda" )
    {

      if( recompute )
      {
        real wMin=omegaMin, wMax=omegaMax;  
        for( int i=0; i<nw; i++ )
        {
          w(i) = wMin + (wMax-wMin)*(i)/real(nw-1);

          real omega = w(i)*omegaScale;          // omega in MKS units (Hz)
          real lambda = (2.*Pi*V0/(omega)) / nm; // lambda in nm 

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

      if( answer=="plot versus omega"  )
        PlotIt::plot(gi, om, fields, title,"omega (PHz)", names,psp );
      else
        PlotIt::plot(gi, lam, fields, title,"lambda (nm)", names,psp );

    }
    
  }
  

  gi.popGUI(); // restore the previous GUI  


  return 0;
}










// Include complex down here to minimize name conflicts
#include <complex>

typedef ::real LocalReal;
typedef ::real OV_real;

// ==========================================================================================
/// \brief Evaluate the real and imaginary parts of epsHat(omega) and 
///       the index of refraction n(omega) = sqrt(epsHat) 
///
/// epsHat(s) = epsInf + SUM_j ( a0(j) + a1(j)*s )/( b0(j) + b1(j)*s + s*s )
/// where s = -i omega
// ==========================================================================================
int DispersiveMaterialParameters::
evalEpsAndN( const LocalReal omega, LocalReal & epsHatr, LocalReal & epsHati, LocalReal & nHatr, LocalReal & nHati  ) const
{

  // std::complex<LocalReal> I(0.0,1.0); 
  std::complex<LocalReal> s(0,-omega); // s = -i omega
  std::complex<LocalReal> epsHat, nHat, a0, a1, b0, b1;


  // LocalReal a0,a1,b0,b1;


  epsHat = epsInf;
  for( int j=0; j<numberOfPolarizationVectors; j++ )
  {
    a0=modelParameters(0,j);
    a1=modelParameters(1,j);
    b0=modelParameters(2,j);
    b1=modelParameters(3,j);

    epsHat += (a0 + a1*s)/( b0 +b1*s +s*s );
  }
  
  epsHatr = std::real(epsHat);
  epsHati = std::imag(epsHat);
  
  nHat = std::sqrt(epsHat);
  nHatr = std::real(nHat);
  nHati = std::imag(nHat);
  

  return 0;
}



// ==========================================================================================
/// \brief Evaluate the real and imaginary parts of the BA material tensor, K(i,j)
///
/// \param i,j (input) : index into the material matrix (i=1,2, or 3 and j=1,2,or 3)
/// K_ij(s) = K0_ij + SUM_m ( a0(m) + a1(m)*s )/( b0(m) + b1(m)*s + s*s )
/// where s = -i omega
// ==========================================================================================
int DispersiveMaterialParameters::
evalMaterialTensor( const LocalReal omega, LocalReal & kHatr, LocalReal & kHati, const int iBA, const int jBA ) const
{

  const MaterialTypeEnum & materialType = dbase.get<MaterialTypeEnum>("materialType");

  if( materialType != bianisotropic )
  {
    printF("DispersiveMaterialParameters::evalMaterialTensor:ERROR: materialType != bianisotropic\n");
    kHatr=0;
    kHati=0;
    return 0;
  }
  if( iBA<1 || iBA>3 | jBA<1 || jBA>3 )
  {
    printF("DispersiveMaterialParameters::evalMaterialTensor:ERROR: iBA=%i or jBA=%i is not between 1 and 3\n",
	   iBA,jBA);
    kHatr=0;
    kHati=0;
    return 0;
  }
  
  
  const RealArray & K0 = dbase.get<RealArray>("K0");
  const RealArray & bianisotropicParameters = dbase.get<RealArray>("bianisotropicParameters");
  const IntegerArray & NpBA = dbase.get<IntegerArray>("NpBA");
	
  // std::complex<LocalReal> I(0.0,1.0); 
  std::complex<LocalReal> s(0,-omega); // s = -i omega
  std::complex<LocalReal> kHat, a0, a1, b0, b1;

  const int ik=iBA-1, jk=jBA-1;  // switch to base 0 
  kHat = K0(ik,jk);
  for( int j=0; j<NpBA(ik,jk); j++ )
  {
    a0=bianisotropicParameters(0,j,ik,jk);
    a1=bianisotropicParameters(1,j,ik,jk);
    b0=bianisotropicParameters(2,j,ik,jk);
    b1=bianisotropicParameters(3,j,ik,jk);

    kHat += (a0 + a1*s)/( b0 +b1*s +s*s );
  }
  
  kHatr = std::real(kHat);
  kHati = std::imag(kHat);
  

  return 0;
}
