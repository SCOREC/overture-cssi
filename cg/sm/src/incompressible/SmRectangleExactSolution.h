#ifndef SM_RECTANGLE_EXACT_SOLUTION
#define SM_RECTANGLE_EXACT_SOLUTION

// ==========================================================================================================
// Class to define exact solutions for CgSm
// 
//   VIBRATIONAL MODES OF A RECTANGLE IN 2D OR 3
//      (1) PERIODIC STRIP IN 2D
//      (1) PERIODIC STRIP IN 3D
// 
// ==========================================================================================================

#include "Cgsm.h"

#define KK_DEBUG
#include "DBase.hh"
using namespace DBase;

typedef ::real LocalReal;

class SmRectangleExactSolution 
{

public:

SmRectangleExactSolution();

~SmRectangleExactSolution();

int evalSolution(Real t, CompositeGrid & cg, int grid, RealArray & ua, 
                 const Index & I1, const Index &I2, const Index &I3, 
                 int numberOfTimeDerivatives /* = 0 */  );

int getParameter( const aString & name, Real & value );
int setParameter( const aString & name, const int value );

int initialize( CompositeGrid & cg, const aString & caseName  );

private:

// The database is a place to store parameters
mutable DataBase dbase;

};


#endif
