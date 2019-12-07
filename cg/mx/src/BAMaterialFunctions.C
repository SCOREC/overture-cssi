// 
// BA related functions from the DispersiveMaterialParameters class
// 
#include "DispersiveMaterialParameters.h"
#include "PlotStuff.h"
#include "display.h"

/* -----------
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
--- */

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


// ==========================================================================================
/// \brief Compute a plane wave solution to Bianisotropic Maxwell Equations (non-dispersive)
///    The solution is
///         [ E ]  = e^{ s*t} e^{i kv.xv } [ ev ]
///         [ H ]                          [ ev ]
///
/// Note: We allow for complex s in the non-dispersive case for the future addition of damping. 
/// 
/// \param kv[3] (input) : wave vector 
/// \param skr,ski  (output) : either real and imaginary parts of s (if computeOption==computeComplexFrequency)
///                 or real and imag parts of k (if  computeOption==computeComplexWaveNumber)
/// \param evr[6], evi[6] (output) : eigenvector ev= evr + I*evi (real and imaginary parts)
/// \param chi(;,2) : real and imaginary parts of the component susceptibilities
/// \param polarization : look for an eigenvector with this polarization if possible 
/// \param computeOption : computeComplexFrequency (i.e. s) or computeComplexWaveNumber (compute wave number)
///
/// The Eigenvalues are found from the matrix (for all unit vectors (kx,ky,kz) )
///      A= - K0^{-1} [ kx*Cx + ky*Cy + kz*Cz ]
///
///       = - K0^{-1} [  0   0   0   0 -kz  ky]
///                   [  0   0   0  kz   0 -kx] 
///                   [  0   0   0 -ky  kx  0 ] 
///                   [  0  kz -ky   0   0  0 ] 
///                   [-kz   0  kx   0   0  0 ] 
///                   [ ky -kx   0   0   0  0 ] 
///        
// ========================================================================================
int DispersiveMaterialParameters::
getBianisotropicPlaneWaveSolution( const real kv[3],
                                   real & skr, real & ski,
                                   real evr[6], real evi[6],
                                   RealArray & chi,
                                   PolarizationEnum polarization /* = noPolarization */ )
{

  const bool & isDispersive=  dbase.get<bool>("isDispersive");
  const DispersionRelationOptionEnum & computeOption= dbase.get<DispersionRelationOptionEnum>("dispersionRelationComputeOption");

  if( isDispersive || computeOption==computeComplexWaveNumber )
  {
    // This version can do either dispersive or non-dispersive
    // Maybe do this always ?? .. but does not find polarizaed state ?
    return getBianisotropicDispersivePlaneWaveSolution( kv, skr,ski,evr,evi,chi,polarization );
  }
  
  if( computeOption==computeComplexWaveNumber )
  {
    printF("getBianisotropicPlaneWaveSolution:ERROR: finish computeComplexWaveNumber for non-dispersive material\n");
    OV_ABORT("ERROR");
  }
  
  real & sr = skr;
  real & si = ski;

  int debug = dbase.get<int>("debug");
  if( debug < 0 ) debug=0;   // current default debug flag 

  if( !dbase.has_key("baPlaneWaveComputed") )
  {
    dbase.put<bool>("baPlaneWaveComputed")=false;
    dbase.put<real>("skrBA")=0.;
    dbase.put<real>("skiBA")=0.;
    dbase.put<real[6]>("evBA");
    dbase.put<real[3]>("kvBA");
  }


  bool & baPlaneWaveComputed = dbase.get<bool>("baPlaneWaveComputed");
  real & skrBA = dbase.get<real>("skrBA");
  real & skiBA = dbase.get<real>("skiBA");
  real * evBA = dbase.get<real[6]>("evBA");
  real * kvBA = dbase.get<real[3]>("kvBA");
  if( baPlaneWaveComputed )
  {
    // The plane wave solution was previously computed -- we can reuse the result if kv is the same as before:
    real kDiff = fabs(kv[0]-kvBA[0]) + fabs(kv[1]-kvBA[1]) + fabs(kv[2]-kvBA[2]);
    if( kDiff < 100.*REAL_EPSILON )
    {
      skr = skrBA;
      ski = skiBA; 
      // si = -omegaBA; // note minus
      for( int m=0; m<6; m++ ){ evr[m]=evBA[m]; evi[m]=0.; }  // 
      if( debug & 1 )
	printF("getBianisotropicPlaneWaveSolution: s=%9.3e + %9.3e I, ev[%6.3g,%6.3g,%6.3g,%6.3g,%6.3g,%6.3g] (using saved values).\n",
	       sr,si,evr[0],evr[1],evr[2],evr[3],evr[4],evr[5]);

      return 0;
    }
  }
  
  baPlaneWaveComputed=true;


  const int ex=0, ey=1, ez=2, hx=3, hy=4, hz=5;

  real kx = kv[0];
  real ky = kv[1];
  real kz = kv[2];
  real kNorm= sqrt(kx*kx + ky*ky + kz*kz);
  kx/=kNorm; ky/=kNorm; kz/=kNorm;

  RealArray K0i;
  getBianisotropicMaterialMatrixInverse( K0i);
  
  RealArray A(6,6);
  RealArray C(6,6); C=0.;

  int n=6, lwork=4*n, info;
  RealArray wr(n), wi(n), work(lwork), vl(6,6), vr(6,6);
  real dummy;


  C(0,hz) = ky; C(0,hy)=-kz;    // (Hz)_y - (Hy)_z
  C(1,hx) = kz; C(1,hz)=-kx;    // (Hx)_z - (Hz)_x
  C(2,hy) = kx; C(2,hx)=-ky;    // (Hy)_x - (Hx)_y
  
  C(3,ez) =-ky; C(3,ey)= kz;    // -[ (Ez)_y - (Ey)_z ]
  C(4,ex) =-kz; C(4,ez)= kx;    // -[ (Ex)_z - (Ez)_x ]
  C(5,ey) =-kx; C(5,ex)= ky;    // -[ (Ey)_x - (Ex)_y ]
  
  C = - C ;  // flip sign so our solution is k*x - omega t 

  myMatrixMultiply(K0i,C,A);

  // ::display(C,"C");
  // ::display(A,"A");

  // "N"  = do not compute left  eigenvectors 
  // "V"  = DO compute right eigenvectors 
  GEEV("N", "V", n, A(0,0), n, wr(0), wi(0), vl(0,0), n, vr(0,0), n, work(0), lwork, info);
  if( info !=0 )
  {
    printF("getBianisotropicPlaneWaveSolution:ERROR return info=%d from eigenvalue routine DGEEV\n",info);
  }
  
  aString fieldName[6]={"Ex","Ey","Ez","Hx","Hy","Hz"};  // 
  if( debug & 2 )
  {
    printF("Normalized kv: kx=%6.3f, ky=%6.3f, kz=%6.3f\n",kx,ky,kz);
    for (int i=0; i<n; i++ )
    {
      printF(" s(%d)=(%9.2e,%9.2e) ev=[%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f]",
	     i,wi(i),wr(i),vr(0,i),vr(1,i),vr(2,i),vr(3,i),vr(4,i),vr(5,i));
      printF(", Polarization= ");
      for( int j=0; j<n; j++ )
      {
	if( fabs(vr(j,i))> 100.*REAL_EPSILON )
	  printF("%s ",(const char*) fieldName[j]);
      }
      printF("\n");
      
    }
  }
  // Check for the requested polarization
  aString polarizationName[7]={"none","ExEyHz", "EyEzHx", "ExEzHy", "HxHyEz", "HyHzEx", "HxHzEy" };   // 
  IntegerArray p(6);
  PolarizationEnum polar[6];
  bool polarizationFound=false;
  for (int i=0; i<n; i++ )
  {
    p=0;
    for( int j=0; j<n; j++ )
    {
      if( fabs(vr(j,i))> 100.*REAL_EPSILON && ( fabs(wr(i))+ fabs(wi(i)) >100.*REAL_EPSILON)  )
        p(j)=1;
    }
    if( sum(p)==3 )
    {
      if(      (p(ex)+p(ey)+p(hz)) == 3 )
        polar[i]=ExEyHzPolarization;

      else if( (p(ey)+p(ez)+p(hx)) == 3 )
        polar[i]=EyEzHxPolarization;

      else if( (p(ex)+p(ez)+p(hy)) == 3 )
        polar[i]=ExEzHyPolarization;

      else if( (p(hx)+p(hy)+p(ez)) == 3 )
        polar[i]=HxHyEzPolarization;

      else if( (p(hy)+p(hz)+p(ex)) == 3 )
        polar[i]=HyHzExPolarization;

      else if( (p(hx)+p(hz)+p(ey)) == 3 )
        polar[i]=HxHzEyPolarization;

      else
      {
	printF("Error -- unexpected polarization!");
	OV_ABORT("getBianisotropicPlaneWaveSolution:ERROR: this case should not happen");
      }
      printF(" i=%d : polarization=%s\n",i,(const char*)polarizationName[polar[i]]);

      if( polar[i]==polarization )
      {
	polarizationFound=true;
      }
    }
    else if( sum(p)==2 )
    {
      if(      (p(ex)+p(ey)+p(hz)) == 2 )
        polar[i]=ExEyHzPolarization;

      else if( (p(ey)+p(ez)+p(hx)) == 2 )
        polar[i]=EyEzHxPolarization;

      else if( (p(ex)+p(ez)+p(hy)) == 2 )
        polar[i]=ExEzHyPolarization;

      else if( (p(hx)+p(hy)+p(ez)) == 2 )
        polar[i]=HxHyEzPolarization;

      else if( (p(hy)+p(hz)+p(ex)) == 2 )
        polar[i]=HyHzExPolarization;

      else if( (p(hx)+p(hz)+p(ey)) == 2 )
        polar[i]=HxHzEyPolarization;

      else
      {
	OV_ABORT("getBianisotropicPlaneWaveSolution:ERROR: this case should not happen");
      }
      printF(" i=%d : polarization=%s\n",i,(const char*)polarizationName[polar[i]]);

      if( polar[i]==polarization )
      {
	polarizationFound=true;
      }

    }
    else
    {
      polar[i]=noPolarization;
    }
  }


  // Choose solution with maximum omega and given polarization (if it exists)
  int iMax=-1;
  real lamMax=-REAL_MAX;
  for( int i=0; i<n; i++ )
  {
    if( wr(i) > lamMax && ( polar[i]==polarization || !polarizationFound ) )
    {
      lamMax=wr(i); iMax=i;
    }
  }
  real omega = wr(iMax)*kNorm;
  sr=0.;
  si=-omega; // note minus 
  
  for( int i=0; i<n; i++ ){ evr[i]=vr(i,iMax); evi[i]=0.; }  // 
  // scale so max-norm is 1 
  real evMax=0.;
  for( int i=0; i<n; i++ ){ if( fabs(evr[i])>evMax ){ evMax=fabs(evr[i]); } }; //
  for( int i=0; i<n; i++ ){ evr[i] /= evMax; }  // 

  if( polarizationFound )
    printF("Requested polarization=[%s] found.\n",(const char*)polarizationName[polarization]);
  else
    printF("Requested polarization=[%s] NOT found.\n",(const char*)polarizationName[polarization]);
  

  printF("Choosing mode i=%d: omega=lam*kNorm=%6.3f (lam=%6.3f,kNorm=%6.3f) evr=[%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f]\n",
	 iMax,omega,wr(iMax),kNorm,evr[0],evr[1],evr[2],evr[3],evr[4],evr[5]);

  // save values 
  // omegaBA = omega;
  skrBA=0.;
  skiBA = -omega;  // note sign 
  for( int m=0; m<6; m++ ){ evBA[m]=evr[m]; }  // 
  for( int m=0; m<3; m++ ){ kvBA[m]=kv[m]; }  // 

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



// lapack routines
#ifdef OV_USE_DOUBLE
  #define ZGESV  EXTERN_C_NAME(zgesv)
  #define ZGEEV  EXTERN_C_NAME(zgeev)
#else
  #define ZGESV  EXTERN_C_NAME(cgesv)
  #define ZGEEV  EXTERN_C_NAME(cgeev)
#endif

extern "C"
{
  void ZGESV( int & n, int & nrhs, std::complex<LocalReal> & a, const int & lda, int & ipvt, std::complex<LocalReal> & b, const int & ldb, int & info );

  // complex eigenvalues
  void ZGEEV( char *jobvl, char* jobvr, int & n, std::complex<LocalReal>  & a, const int & lda,
              std::complex<LocalReal> & wc, std::complex<LocalReal> &vl, int & ldvl, std::complex<LocalReal> & vr, int & ldvr,
	      std::complex<LocalReal> & cwork, int & lwork, LocalReal & rwork,  int & info );
}



// ==========================================================================================
/// \brief Compute a plane wave solution to Bianisotropic Maxwell Equations (DISPERSIVE)
///    The solution is
///         [ E ]  = e^{ s*t} e^{i kv.xv } [ ev ]
///         [ H ]                          [ ev ]
///
/// \param kv[3] (input or output) : wave vector 
/// \param skr,ski  (output) : either real and imaginary parts of s (if computeOption==computeComplexFrequency)
///                 or real and imag parts of k (if  computeOption==computeComplexWaveNumber)
/// \param evr[6], evi[6] (output) : eigenvector ev= evr + I*evi (real and imaginary parts)
/// \param chi(;,2) : real and imaginary parts of the component susceptibilities
/// \param polarization : look for an eigenvector with this polarization if possible 
/// \param computeOption : computeComplexFrequency (i.e. s) or computeComplexWaveNumber (compute wave number)
/// 
/// The Eigenvalues are found from the matrix (for all unit vectors (kx,ky,kz) )
///      A (s) = i*s K(s) + L(k)
///
///       L(k) = [  0   0   0   0 -kz  ky]
///              [  0   0   0  kz   0 -kx] 
///              [  0   0   0 -ky  kx  0 ] 
///              [  0  kz -ky   0   0  0 ] 
///              [-kz   0  kx   0   0  0 ] 
///              [ ky -kx   0   0   0  0 ] 
///        
///
/// Notes: The polarization terms will be defined from
///          P_{i,j,n} = Chi_{i,j,m} U_i
// ========================================================================================
int DispersiveMaterialParameters::
getBianisotropicDispersivePlaneWaveSolution( const LocalReal kv[3],
                                             LocalReal & skr, LocalReal & ski,
                                             LocalReal evr[6], LocalReal evi[6],
                                             RealArray & chi,
                                             PolarizationEnum polarization /* = noPolarization */ ) 
{
  // *wdh* Dec 3, 2019 -- this version can also do non-dispersive case too --
  
  // if( !dbase.has_key("bianisotropicParameters" ) )
  // {
  //   printF("DispMatParam:getBianisotropicDispersivePlaneWaveSolution: this is not a dispersive BA material\n");
  //   OV_ABORT("Error");
  // }

  const bool & isDispersive=  dbase.get<bool>("isDispersive");

  const DispersionRelationOptionEnum & computeOption= dbase.get<DispersionRelationOptionEnum>("dispersionRelationComputeOption");

  int debug = dbase.get<int>("debug");
  if( debug < 0 ) debug=0;   // current default debug flag 

  if( debug & 1 )
  {
    printF("***DispMatParam: entering getBianisotropicDispersivePlaneWaveSolution");
    if( computeOption==computeComplexFrequency )
      printF(", computeComplexFrequency\n ****\n");
    else 
      printF(", computeComplexWaveNumber\n ****\n");
  }
  
  LocalReal & sr = skr;
  LocalReal & si = ski;

  LocalReal kr=1., ki=0;  // If we compute k = kr + ki*I
  

  std::complex<LocalReal> I(0,1.); 
  // printF(" I= %f + %f \n",std::real(I),imag(I));

  // --- we save any previously computed values ----
  if( !dbase.has_key("baPlaneWaveFrequencyComputed") )
  {
    dbase.put<bool>("baPlaneWaveFrequencyComputed")=false;
    dbase.put<bool>("baPlaneWaveNumberComputed")=false;
    dbase.put<LocalReal>("srBA")=0.;
    dbase.put<LocalReal>("siBA")=0.;
    dbase.put<LocalReal>("krBA")=0.;
    dbase.put<LocalReal>("kiBA")=0.;
    dbase.put<LocalReal[6]>("evrBA");
    dbase.put<LocalReal[6]>("eviBA");
    dbase.put<LocalReal[3]>("kvBA");
  }

  bool & baPlaneWaveFrequencyComputed = dbase.get<bool>("baPlaneWaveFrequencyComputed");
  bool & baPlaneWaveNumberComputed = dbase.get<bool>("baPlaneWaveNumberComputed");
  LocalReal & srBA  = dbase.get<LocalReal>("srBA");
  LocalReal & siBA  = dbase.get<LocalReal>("siBA");
  LocalReal & krBA  = dbase.get<LocalReal>("krBA");
  LocalReal & kiBA  = dbase.get<LocalReal>("kiBA");
  LocalReal * evrBA = dbase.get<LocalReal[6]>("evrBA");
  LocalReal * eviBA = dbase.get<LocalReal[6]>("eviBA");
  LocalReal * kvBA  = dbase.get<LocalReal[3]>("kvBA");

  if( (computeOption==computeComplexFrequency  && baPlaneWaveFrequencyComputed) ||
      (computeOption==computeComplexWaveNumber && baPlaneWaveNumberComputed ) )
  {
    // The plane wave solution was previously computed -- we can reuse the result if kv is the same as before:

    LocalReal kvDiff = fabs(kv[0]-kvBA[0]) + fabs(kv[1]-kvBA[1]) + fabs(kv[2]-kvBA[2]);

    // For comuting the waveNumber also check that sr,si are the same
    LocalReal sDiff=0.;
    if( baPlaneWaveNumberComputed )
      sDiff = fabs(sr-srBA) + fabs(si-siBA);
        
    if( max(kvDiff,sDiff) < 100.*REAL_EPSILON )
    {
      if( computeOption==computeComplexFrequency && baPlaneWaveFrequencyComputed )
      {
        skr = srBA;
        ski = siBA;
        if( debug & 1 )
          printF("getDispersiveBianisotropicPlaneWaveSolution: s=%9.3e + %9.3e I (using saved values).\n",skr,ski);
      }
      else
      {
        skr = krBA;
        ski = kiBA;
        if( debug & 1 )
          printF("getDispersiveBianisotropicPlaneWaveSolution: k=%9.3e + %9.3e I (using saved values).\n",skr,ski);
      }
  
      for( int m=0; m<6; m++ ){ evr[m]=evrBA[m]; evi[m]=eviBA[m];  }  // 

      if( dbase.has_key("baChi") )
      {
        RealArray & baChi = dbase.get<RealArray>("baChi");
        chi.redim(0);
        chi=baChi;
      }
      

      return 0;
    }
  }
  
  if( debug <= 0 ) debug=3;   // current default debug flag 


  const int ex=0, ey=1, ez=2, hx=3, hy=4, hz=5;

  LocalReal kx = kv[0];
  LocalReal ky = kv[1];
  LocalReal kz = kv[2];
  LocalReal kNorm= 1.;
  // LocalReal kNorm= sqrt(kx*kx + ky*ky + kz*kz);
  kx/=kNorm; ky/=kNorm; kz/=kNorm;

  LocalReal kSize= sqrt(kx*kx + ky*ky + kz*kz);
  const LocalReal tolk  =        kSize*REAL_EPSILON*100.;  // tolerance for quantities proportional to k 
  const LocalReal tolk3 = pow(kSize,3)*REAL_EPSILON*100.;  // tolerance for quantities proportional to k^3 

  // For non-dipserive case: 
  IntegerArray Np0(6,6); Np0=0;
  RealArray biPar0;

  const  RealArray & K0                       = dbase.get<RealArray>("K0")         ;
  const IntegerArray & NpBA                   = isDispersive ? dbase.get<IntegerArray>("NpBA") : Np0;
  const  RealArray & bianisotropicParameters  = isDispersive ? dbase.get<RealArray>("bianisotropicParameters") : biPar0;

  int numberOfPolarizationTerms=0;
  std::complex<LocalReal> s, sp2;

  std::complex<LocalReal> Ac[36], A2c[36], bc[6];
#define A(i,j) Ac[(i)+6*(j)]
#define A2(i,j) A2c[(i)+6*(j)]
#define b(i) bc[(i)]


  if( computeOption==computeComplexFrequency )
  {
    // -------------- Compute s = sr + si I  -------------------

    baPlaneWaveFrequencyComputed=true;
    // --- Define some variables used in the include file below ----
    int NpMax=3; // max number of polarization terms expected in include files 
    // Range I3(1,3), I6(1,6);
    RealArray a0v(NpMax,6,6), a1v(NpMax,6,6), b0v(NpMax,6,6), b1v(NpMax,6,6);
    // RealArray eps(I3,I3);
    // LocalReal mu1,mu2,mu3;
  
    const LocalReal kxp2=kx*kx,     kyp2=ky*ky,     kzp2=kz*kz;
    const LocalReal kxp3=kxp2*kx,   kyp3=kyp2*ky,   kzp3=kzp2*kz;
    const LocalReal kxp4=kxp2*kxp2, kyp4=kyp2*kyp2, kzp4=kzp2*kzp2;

    a0v=0.; a1v=0.;  b0v=0.; b1v=0.;
    // eps=1.; mu1=1.; mu2=1.; mu3=1.;
  
    // do this for now
    // eps(1,1) = K0(1,1);
    // eps(2,2) = K0(2,2);
    // eps(3,3) = K0(3,3);
  
    for( int k1=0; k1<6; k1++ )
    {
      for( int k2=0; k2<6; k2++ )
      {
        for( int n=0; n<NpBA(k1,k2); n++ )
        {
          numberOfPolarizationTerms++;

          a0v(n,k1,k2) = bianisotropicParameters(0,n,k1,k2);
          a1v(n,k1,k2) = bianisotropicParameters(1,n,k1,k2);
          b0v(n,k1,k2) = bianisotropicParameters(2,n,k1,k2);
          b1v(n,k1,k2) = bianisotropicParameters(3,n,k1,k2);
	  if( debug & 1 )
	    printF("Setting: K(%d,%d) : n=%d: [a0,a1,b0,b1]=[%e,%e,%e,%e]\n",k1,k2,n,
		   a0v(n,k1,k2),a1v(n,k1,k2),b0v(n,k1,k2),n,b1v(k1,k2));
        }

      }
    }
    
    LocalReal * pa0 = a0v.getDataPointer();
#define a0(i,j,k) pa0[(i)+NpMax*((j)+6*(k))]

    LocalReal * pa1 = a1v.getDataPointer();
#define a1(i,j,k) pa1[(i)+NpMax*((j)+6*(k))]

    LocalReal * pb0 = b0v.getDataPointer();
#define b0(i,j,k) pb0[(i)+NpMax*((j)+6*(k))]

    LocalReal * pb1 = b1v.getDataPointer();
#define b1(i,j,k) pb1[(i)+NpMax*((j)+6*(k))]

    LocalReal * pK0 = K0.getDataPointer();
#define K0a(i,j) pK0[(i)+6*(j)]


    std::complex<LocalReal> *pcmc, *pwc, *pvlc, *pvrc, *pworkc;
    // pcmc = new std::complex<LocalReal> [Nd*Nd];  
#define cmc(i,j) pcmc[(i)+Nd*(j)]
#define cm(i,j) pcmc[(i)+Nd*(j)]

    // int n=0; // fix me 
    std::complex<LocalReal> K15, K51, K24, K42;  // ** fix me **
  
#define BACASE 0 

#if BACASE == 0
  #include "baDispersionRelation0CompanionMatrixOpt.h"
#elif BACASE == 1
  #include "baDispersionRelation1CompanionMatrixOpt.h"
#elif BACASE == 2
 #include "baDispersionRelation2CompanionMatrixOpt.h"
#elif BACASE == 3
  //  Derek's BA form 
 #include "baDispersionRelation3CompanionMatrixOpt.h"
#endif
    
    // ::display(eps,"eps","%6.3f ");
    // ::display(cm,"Companion matrix","%6.3f ");


    int n=Nd, info;   // ** FIX ME ***
    pvrc = new std::complex<LocalReal> [Nd*Nd];
    pvlc = pvrc; // not used 
    pwc = new std::complex<LocalReal> [Nd];

    int lworkc = 2*Nd;
    RealArray rwork(2*Nd);
    pworkc = new std::complex<LocalReal> [lworkc];
#define wc(i) pwc[(i)]
#define vlc(i,j) pvlc[(i)+Nd*(j)]
#define vrc(i,j) pvrc[(i)+Nd*(j)]
#define workc(i) pworkc[(i)]

  

    for( int i=0; i<n; i++ )
      for( int j=0; j<n; j++ )
      {
        // cmc(i,j)=cm(i,j);
        // printF(" cm(%d,%d) = %e + %e I\n",i,j,std::real(cm(i,j)),imag(cm(i,j)));
      }
  
    // Compute eigenvalues of a complex matrix -- do we need eigenvectors ????
    ZGEEV( "N", "V", n, cmc(0,0), n, wc(0), vlc(0,0), n, vrc(0,0), n, workc(0), lworkc, rwork(0), info); 
    if( info !=0 )
    {
      printF("getBianisotropicPlaneWaveSolution:ERROR return info=%d from eigenvalue routine ZGEEV\n",info);
    }
  

    // aString fieldName[6]={"Ex","Ey","Ez","Hx","Hy","Hz"};  // 
    if( debug & 2 )
    {
      printF("Normalized kv: kx=%6.3f, ky=%6.3f, kz=%6.3f\n",kx,ky,kz);
      for (int i=0; i<n; i++ )
      {
        //printF("old: s(%3d)=(%9.2e,%9.2e) v=[",i,wr(i),wi(i));
        //for( int j=0; j<n; j++ )
//	printF("%6.3f,",vr(j,i));
        //printF("]\n");

        if( true )
        {
          LocalReal sr = std::real(wc(i)), si=imag(wc(i));
          printF(" s(%3d)=(%9.2e,%9.2e)  ",i,sr,si);
	  if( sr >  tolk ) printF(" *growth*  sr>0 ");
	  if( sr < -tolk ) printF(" *decay *  sr<0");
	  printF("\n");
        }
        else
        {
          printF(" s(%3d)=(%9.2e,%9.2e) v=[",i,std::real(wc(i)),imag(wc(i)));
          for( int j=0; j<n; j++ )
            printF("%6.3f+%6.3f I,",std::real(vrc(j,i)),imag(vrc(i,j)));
          printF("]\n");
        }

      }
    }

    RealArray wr(n), wi(n);
    for (int i=0; i<n; i++ )
    {
      wr(i) = std::real(wc(i));
      wi(i) = imag(wc(i));
    }
  

    // choose the eigenvalue with largest negative imaginary part
  
    LocalReal wiMin=REAL_MAX;
    int iMax=-1;
    for (int i=0; i<n; i++ )
    {
      if( wi(i) < wiMin )
      {
        wiMin= wi(i);
        iMax=i;
      }
    
    }
    int eigIndex = iMax;

    // eigIndex=0;   // ** TEMP 
  

    s = wr(eigIndex) + I*wi(eigIndex);

 
    s = s*kNorm;  // unscale s 
 
    sr = std::real(s);
    si =      imag(s);

 
    // -- clean up ---
    delete [] pcmc;
    delete [] pvrc;
    delete [] pwc;
    delete [] pworkc;


  }
  // end computeOption==computeComplexFrequency 
 


  if( computeOption==computeComplexWaveNumber )
  {
    // --- Compute k given s and wave vector direction -----
    baPlaneWaveNumberComputed=true;
   
    LocalReal kx = kv[0];
    LocalReal ky = kv[1];
    LocalReal kz = kv[2];
    LocalReal kNorm= sqrt(kx*kx + ky*ky + kz*kz);
    kx/=kNorm; ky/=kNorm; kz/=kNorm;

    if( debug & 1 )
      printF("Compute complex wave number for s=[%e,%e] and kv=[%e,%e,%e]\n",sr,si,kx,ky,kz);
  
    const LocalReal kxp2=kx*kx,     kyp2=ky*ky,     kzp2=kz*kz;
    const LocalReal kxp3=kxp2*kx,   kyp3=kyp2*ky,   kzp3=kzp2*kz;
    const LocalReal kxp4=kxp2*kxp2, kyp4=kyp2*kyp2, kzp4=kzp2*kzp2;

    s = sr + si*I; //  = imag(s)*I;  // remove real part 

    sp2 = s*s;

  
    std::complex<LocalReal> pKs[36];
#define Ks(i,j) pKs[(i)+6*(j)]
    numberOfPolarizationTerms=0;
    for( int k1=0; k1<6; k1++ )
    {
      for( int k2=0; k2<6; k2++ )
      {
        Ks(k1,k2)=K0(k1,k2);
        for( int n=0; n<NpBA(k1,k2); n++ )
        {
          numberOfPolarizationTerms++;
          LocalReal a0 = bianisotropicParameters(0,n,k1,k2);
          LocalReal a1 = bianisotropicParameters(1,n,k1,k2);
          LocalReal b0 = bianisotropicParameters(2,n,k1,k2);
          LocalReal b1 = bianisotropicParameters(3,n,k1,k2);

          Ks(k1,k2) += (a0+a1*s)/(b0+b1*s+s*s);

        }
      }
    }

    const int Nd=4;  // dimension of companion matrix (4 roots)

    std::complex<LocalReal> *pcmc, *pwc, *pvlc, *pvrc, *pworkc;
    // pcmc = new std::complex<LocalReal> [Nd*Nd];  
#define cmc(i,j) pcmc[(i)+Nd*(j)]
#define cm(i,j) pcmc[(i)+Nd*(j)]

    int n=Nd, info;   // ** FIX ME ***
    pvrc = new std::complex<LocalReal> [Nd*Nd];
    pvlc = pvrc; // not used 
    pwc = new std::complex<LocalReal> [Nd];

    int lworkc = 2*Nd;
    RealArray rwork(2*Nd);
    pworkc = new std::complex<LocalReal> [lworkc];
#define wc(i) pwc[(i)]
#define vlc(i,j) pvlc[(i)+Nd*(j)]
#define vrc(i,j) pvrc[(i)+Nd*(j)]
#define workc(i) pworkc[(i)]

    pcmc = new std::complex<LocalReal> [Nd*Nd];
 
#include "baDispersionRelationWaveNumberPolynomialMatrix.h"

    // Compute eigenvalues of a complex matrix 
    ZGEEV( "N", "N", n, cmc(0,0), n, wc(0), vlc(0,0), n, vrc(0,0), n, workc(0), lworkc, rwork(0), info); 
    if( info !=0 )
    {
      printF("getBianisotropicPlaneWaveSolution:ERROR return info=%d from eigenvalue routine ZGEEV\n",info);
    }

    if( debug & 1 ) 
      printF("\n --- Compute wave number k given s=%e + %e I, and direction kHat.\n",std::real(s),imag(s));

    if( debug & 2 )
    {
      printF("Normalized kv: kx=%6.3f, ky=%6.3f, kz=%6.3f\n",kx,ky,kz);
      for (int i=0; i<n; i++ )
      {
	LocalReal kr=std::real(wc(i)), ki=imag(wc(i));
        printF(" k(%3d)=(%9.2e,%9.2e)",i,kr,ki);
        // For s = -i omega 
	//    no-growth:  omega *Im(k)/Re(k) >=0 
        // or             omega *Im(k) *Re(k) >=0 
	// More general: (set bamxp/notes)
	//    no growth: - Im(s) * Im(k) * Re(k) - Re(k)^2 Re(s)  >= 0 
	if( -si*ki*kr - kr*kr*sr < -tolk3 ) printF("  **growth**  -si*ki*kr - kr*kr*sr < 0 ");
	if( -si*ki*kr - kr*kr*sr >  tolk3 ) printF("  **decay **  -si*ki*kr - kr*kr*sr > 0 ");
	printF("\n");

      }
    }

    // ---- Pick a k ----
    int kIndex=0; 
    if( mode>=0 && mode<=4 )
    {
      kIndex=mode;  // user specified mode
    }
    else
    {
      LocalReal kiMax=0.;      
      for( int i=0; i<n; i++ )
      {
        LocalReal kr0 = std::real(wc(i));
        LocalReal ki0 = std::imag(wc(i));
        if( kr0 >0. && fabs(ki0) >= kiMax ) // chose a mode with kr >0 (right moving) and dispersive, ki !=0 
        {
          kIndex=i; kiMax=fabs(ki0);
        }
      }
    }
    
    kr = std::real(wc(kIndex));
    ki = std::imag(wc(kIndex));

    if( debug & 2 )
      printF("Choosing k=%e + %e I  (kIndex=%d) (mode=%d)\n",kr,ki,kIndex,mode);

    // -- clean up ---
    delete [] pcmc;
    delete [] pvrc;
    delete [] pwc;
    delete [] pworkc;

  }
 
  if( true )
  {
    // ---- Compute the EM EIGENVECTOR for a given s or k -----
    // ** new way **

    // --- First fill in the material matrix K(s) ---

    std::complex<LocalReal> pKs[36];
#define Ks(i,j) pKs[(i)+6*(j)]
    numberOfPolarizationTerms=0;
    for( int k1=0; k1<6; k1++ )
    {
      for( int k2=0; k2<6; k2++ )
      {
        Ks(k1,k2)=K0(k1,k2);
        for( int n=0; n<NpBA(k1,k2); n++ )
        {
          numberOfPolarizationTerms++;
          LocalReal a0 = bianisotropicParameters(0,n,k1,k2);
          LocalReal a1 = bianisotropicParameters(1,n,k1,k2);
          LocalReal b0 = bianisotropicParameters(2,n,k1,k2);
          LocalReal b1 = bianisotropicParameters(3,n,k1,k2);

          Ks(k1,k2) += (a0+a1*s)/(b0+b1*s+s*s);

        }
      }
    }

    // kx,ky,kz must be non-normalized here 
    std::complex<LocalReal> kc;
   if( computeOption==computeComplexFrequency )
     kc=1.;
   else
   {
     // if we computed the complex wave number then scale by this
     LocalReal kNorm = sqrt( SQR(kv[0])+SQR(kv[1])+SQR(kv[2]) );
     kc = (kr + ki*I)/kNorm;
   }
   
    std::complex<LocalReal> kx = kv[0]*kc;
    std::complex<LocalReal> ky = kv[1]*kc;
    std::complex<LocalReal> kz = kv[2]*kc;

    // Evaluate A = I*s*Ks + L(kv) 
#include "baEigMatrix.h"


    int n=6;
    for (int i=0; i<n; i++ )
    {
      for( int j=0; j<n; j++ )
      {
        A2(i,j)=A(i,j); // save for checking 
      }
    }
    if( debug & 2 )
    {
      printF("Here is A(s) for eigenvalue s=[%8.5f,%8.5f]\n",std::real(s),imag(s));
      for (int i=0; i<n; i++ )
      {
	for( int j=0; j<n; j++ )
	{
	  printF("A%d%d=[%5.2f,%5.2f], ",i,j,std::real(A(i,j)),imag(A(i,j)));
	}
	printF("\n");
      }
    }
    

    // Solve A x = b 
    int ipvt[6], info;
    int nrhs=1;
    for (int i=0; i<n; i++ )
    {
      b(i)=1.;  // guess at eigenvector 
    }
 
    ZGESV( n, nrhs, A(0,0), n, ipvt[0], b(0), n, info );
    if( info!=0 )
    {
      printF("ERROR return from laplack routine ZGESV, LU factorization: info=%d\n",info);
      OV_ABORT("ZGESV: info !=0 ");
    }
    // normalize the eigenvector: 
    LocalReal bNorm=0., br=0., bi=0., bMax=-1.;
    int iMax=-1;
    for (int i=0; i<n; i++ )
    {
      LocalReal bir=std::real(b(i)), bii=std::imag(b(i));
      br += SQR(bir);
      bi += SQR(bii);
      if( fabs(b(i)) > bMax )
      {
        bMax=abs(b(i));
        iMax=i;
      }
      // bNorm += abs(b(i))*abs(b(i));
      // printF(" b(%d) =[%5.2f,%5.2f] |b(i)|=%9.2e iMax=%d\n",i,std::real(b(i)),imag(b(i)),abs(b(i)),iMax);

    }
 
    bNorm=sqrt(br+bi);
    std::complex<LocalReal> bScale=b(iMax);
    for (int i=0; i<n; i++ )
    {
      // b(i) /= bNorm* b(iMax)/abs(b(iMax));  // scale to largest entry is 1
      b(i) /= bScale;  // scale to largest entry is 1
    }
    if( debug & 2 )
    {
      printF("Eigenvector:\n");
      for (int i=0; i<n; i++ )
      {
	printF(" ev(%d) =[%9.2e,%9.2e], ",i,std::real(b(i)),imag(b(i)));
      }
      printF("\n");
    }
    

    // CHECK that A*b = 0 
    std::complex<LocalReal> resc[6];
#define res(i) resc[(i)]
    LocalReal resMax=0.;
    for (int i=0; i<n; i++ )
    {
      res(i)=0;
      for (int j=0; j<n; j++ )
      {
        res(i) += A2(i,j)*b(j);
      }
      resMax=max(resMax,abs(res(i)));
    }

    if( debug & 2 )
      printF("check... norm(A*b)=%8.2e (should be zero)\n",resMax);

    const LocalReal tol = 100.*sqrt(REAL_EPSILON);  // eigenvalues are not so accurate when there is a double root
    if( resMax > tol )
    {
      printF("DispMatPar: ERROR computing eigenvalue for dispersive plane wave: A(s)*b =%9.2e is too big"
             "tol=%8.2e\n",resMax,tol);
      OV_ABORT("error");
    }
 
    for (int i=0; i<n; i++ )
    {
      evr[i]=std::real(b(i));
      evi[i]=std::imag(b(i));
    }

  } // end compute eigenvector 

 

  // --- save values so we don't need to recompute ---
  if( computeOption==computeComplexFrequency )
  {
    // save complex s 
    srBA = sr;
    siBA = si;
  }
  else
  {
    // save complex s 
    srBA = sr;
    siBA = si;

    // save complex k 
    LocalReal kNorm = sqrt( SQR(kv[0])+SQR(kv[1])+SQR(kv[2]) );
    krBA = kr/kNorm;
    kiBA = ki/kNorm;

    // We return kr,ki
    skr=krBA;
    ski=kiBA;
  }
 
  for( int m=0; m<6; m++ ){ evrBA[m]=evr[m]; eviBA[m]=evi[m];  }  // 
  for( int m=0; m<3; m++ ){ kvBA[m]=kv[m]; }  // 

  if( true || computeOption==computeComplexFrequency )
  {
   
    // ----- compute susceptibilities ----
    if( numberOfPolarizationTerms>0 )
    {
      chi.redim(numberOfPolarizationTerms*2,2);
      if( debug & 1 )
	printF("susceptibilities: chi = (a0+a1*s)/(b0 + b1*s + s*2) * U_i \n");
      std::complex<LocalReal> chin;
      int m=0;
      for( int k1=0; k1<6; k1++ )
      {
        for( int k2=0; k2<6; k2++ )
        {
          for( int n=0; n<NpBA(k1,k2); n++ )
          {
            LocalReal a0 = bianisotropicParameters(0,n,k1,k2);
            LocalReal a1 = bianisotropicParameters(1,n,k1,k2);
            LocalReal b0 = bianisotropicParameters(2,n,k1,k2);
            LocalReal b1 = bianisotropicParameters(3,n,k1,k2);

            // P_{ijn}
            int ec = k2; // This GDM term multiplies this E or H component

            chin = (a0+a1*s)/(b0+b1*s+s*s) * b(ec);
            chi(m,0) = std::real(chin);
            chi(m,1) = std::imag(chin);
            m++;
	    if( debug & 2 )
	      printF("P: Chi(n=%d,k1=%d,k2=%d) = %10.3e + %10.3e I  (ec=%d)\n",n,k1,k2,std::real(chin),imag(chin),ec);

            // Q_{ijn} = time-derivative of P
            chin = s*chin;  
            chi(m,0) = std::real(chin);
            chi(m,1) = std::imag(chin);
            m++;
	    if( debug & 2 )
              printF("Q: Chi(n=%d,k1=%d,k2=%d) = %10.3e + %10.3e I  (ec=%d)\n",n,k1,k2,std::real(chin),imag(chin),ec);
         
          }
        }
      }

      // -- save chi ---
      if( !dbase.has_key("baChi") )
        dbase.put<RealArray>("baChi");
      RealArray & baChi = dbase.get<RealArray>("baChi");
      baChi.redim(0);
      baChi=chi;

    }
  }



  return 0;
}

