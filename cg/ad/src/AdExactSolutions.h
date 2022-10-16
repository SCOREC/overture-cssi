#ifndef SM_CYLINDER_EXACT_SOLUTION
#define SM_CYLINDER_EXACT_SOLUTION

// ==========================================================================================================
// Class to define exact solutions for CgAd
// 
//   (1) DOUBLE ANNULUS
// 
// ==========================================================================================================

#include "Cgad.h"

#define KK_DEBUG
#include "DBase.hh"
using namespace DBase;

typedef ::real LocalReal;

class AdExactSolutions 
{

  public:
  
  AdExactSolutions();
  ~AdExactSolutions();
  
  int evalSolution(Real t, CompositeGrid & cg, int grid, RealArray & ua, 
                   const Index & I1, const Index &I2, const Index &I3, 
                   int numberOfTimeDerivatives /* = 0 */  );
  
  int getParameter( const aString & name, Real & value );
  
  int initialize( CompositeGrid & cg, const aString & solution, const aString & caseName  );

  // int setDomain( int domain );
  
  private:

  int setup();    
  
  // The database is a place to store parameters
  mutable DataBase dbase;

};


#endif
