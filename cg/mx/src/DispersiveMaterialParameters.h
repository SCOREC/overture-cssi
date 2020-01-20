// Class to define parameters of a dispersive material
//
// Generalized Dispersion Model:
//       E_tt - c^2 Delta(E) = - (1/eps) P_tt
//       P_tt + b1 P_1 + b0 = eps*( a0*E + a1*E_t )
// 
#include "Overture.h"

#define KK_DEBUG
#include "DBase.hh"
using namespace DBase;

class DispersiveMaterialParameters
{
public:


  enum MaterialTypeEnum
  {
    isotropic,
    bianisotropic
  };
  
  enum PolarizationEnum
  {
    noPolarization,
    ExEyHzPolarization,        // TEz 
    EyEzHxPolarization,
    ExEzHyPolarization,
    HyHzExPolarization,
    HxHzEyPolarization,
    HxHyEzPolarization
  };
  
  // We have choices when solving the dispersion relation
  enum DispersionRelationOptionEnum
  {
    computeComplexFrequency=0,
    computeComplexWaveNumber
  };
  

  DispersiveMaterialParameters();
  DispersiveMaterialParameters(const DispersiveMaterialParameters& x);
  ~DispersiveMaterialParameters();

  DispersiveMaterialParameters& operator=( const DispersiveMaterialParameters& x);

  // diplay parameters
  int display( FILE *file = stdout, const aString & label = nullString ) const;

  // Evaluate the real and imaginary parts of epsHat(omega) and nHat(omega) 
  int evalEpsAndN( const real omega, real & epsHatr, real & epsHati, real & nHatr, real & nHati ) const;

  // Evaluate the eigenvalues for a BA material 
  int evalBianisotropicEigenValues( int numberOfDimensions, real & reLambda, real & imLambda );

  // compute the time-stepping eigenvalue for the BA equations 
  int evalBianisotropicTimeSteppingEigenvalue( MappedGrid & mg, const int orderOfAccuracy, real & reLambda, real & imLambda );

  // Evaluate the real and imaginary parts of the BA material tensor, K(iBA,jBA) iBA=1,2,3,  jBA=1,2,3
  int evalMaterialTensor( const real omega, real & kHatr, real & kHati, const int iBA, const int jBA ) const;

  // Return the complex wave number given s=(sr,si)
  int evaluateComplexWaveNumber( const real c, const real & sr, const real & si, 
                                 real & kr, real &ki, real chir[], real chii[], real & chiSumr, real & chiSumi );
  // return the and imaginary parts of "s" in the dispersion relation
  int evaluateDispersionRelation( const real c, const real k, real & sr, real & si, real chir[], real chii[], 
                                  real & chiSumr, real & chiSumi  );

  // Evaluate the complex index of refraction for scattering off dispersive dielectrics
  int evaluateComplexIndexOfRefraction( const real mu1, const real eps1, const real chi1r, const real chi1i, 
                                        const real mu2, const real eps2, const real chi2r, const real chi2i, 
                                        real & mr, real & mi ) const;


  // // *old* return the and imaginary parts of "s" in the dispersion relation
  // int 
  // computeDispersionRelation( const real cc, const real eps, const real mu, const real k, 
  //                            real & reS, real & imS );

  int findPeriodicSolution( real kv[3], real & sr, real & si, real & kr, real & ki, real evr[6], real evi[6], RealArray & chi );

  int getBianisotropicMaterialMatrix( RealArray & K0 ) const;
  int getBianisotropicMaterialMatrixInverse( RealArray & K0i );

  IntegerArray& getBianisotropicNp() const;

  RealArray & getBianisotropicGDMParameters() const;


  // Return the bi-anisotropic parameters 
  int getBianisotropicParameters( RealArray & K0, RealArray & bianisotropicParameters, IntegerArray & Np );

  int getBianisotropicPlaneWaveSolution( const real kv[3],
					 real & skr, real & ski,
					 real evr[6], real evi[6],
                                         RealArray & chi,
					 PolarizationEnum polarization = noPolarization );
  
  int getBianisotropicDispersivePlaneWaveSolution( const real kv[3],
                                                   real & skr, real & ski,
                                                   real evr[6], real evi[6],
                                                   RealArray & chi,
                                                   PolarizationEnum polarization = noPolarization );

  // Return the isotropic parameters
  int getIsotropicParameters( int & Np,  real & epsInf, RealArray & modelParams ) const;

  real getEpsInf() const;
  real getMuInf() const;

  real getAlphaP() const;

  // Return the material name 
  aString getMaterialName() const;

  // Return the material type
  MaterialTypeEnum getMaterialType() const;

  bool isDispersiveMaterial() const;

  // Plot dispersive properties versus frequency or wavelength
  int update( GenericGraphicsInterface & gi );

  
  // read dispersive material parameters from a file 
  int readFromFile( const aString & fileName, int numberOfPolarizationVectorsRequested = -1 );

  int setDebug( int debugFlag );

  int setDispersionRelationComputeOption( DispersionRelationOptionEnum computeOption );

  int setEpsInf( const real epsInf );
  int setMuInf( const real muInf );

  int setMaterialType( MaterialTypeEnum matType );

  int setMode( const int modeToChoose );

  int setNumberOfPolarizationVectors( const int numPolarizationVectors );

  int setParameter( const real alphaP );

  int setParameters( const int eqn, const real a0, const real a1, const real b0, const real b1 );

  // Set a velocity and length scale for non-dimensionalization of the parameters 
  int setScales( real velocityScale, real lengthScale );

  // //  *OLD*
  // int 
  // computeDispersivePlaneWaveParameters( const real cc, const real eps, const real mu, const real k, 
  //                                       real & omegar, real & omegai );
  int setParameters( const real a0, const real a1, const real b0, const real b1 );

  // write parameters to a material database file
  int writeToFile( const aString & fileName );


  static void myMatrixMultiply( RealArray & a, RealArray & b, RealArray & c );

  // Data members -- make public for now
  public:

  real epsInf, muInf;
  real alphaP;
  real gamma, omegap;  // Drude-Lorentz model

  // general dispersive model parameters:
  //   modelParameters(i,k)  : i=0,1,2,3,numberOfModelParameters-1 are the parmeters in the equation 
  //                           for P_k , k=1,2,...,numberOfPolarizationVectors
  //   modelParameter(0,k) : a0
  //   modelParameter(1,k) : a1
  //   modelParameter(2,k) : b0
  //   modelParameter(3,k) : b1
  // 
  int numberOfPolarizationVectors;
  int numberOfModelParameters;
  RealArray modelParameters;

  // Save the root after it has been computed to save computations
  int mode;  // which root to choose.
  bool rootComputed;
  real ck0,sr0,si0, chiSumr0,chiSumi0;
  RealArray chir0,chii0;

  // The database is a place to store parameters
  mutable DataBase dbase;

};
  
