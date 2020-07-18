#ifndef SLABS_EXACT_SOLUTIONS
#define SLABS_EXACT_SOLUTIONS


// ===============================================================================
// Class to define exact solutions to Maxwell's equations
// 
//     Scattering from one or more SLABS
//
//     -------------------------------------------------------------------
//     |                 |     |       |      |                          |
//     |         0       |  1  |   2   |   3  |           4              |
//     |                 |     |       |      |                          |
//     -------------------------------------------------------------------
//
// ===============================================================================

#include "Maxwell.h"

#define KK_DEBUG
#include "DBase.hh"
using namespace DBase;

class MxParameters;

typedef ::real LocalReal;

class SlabsExactSolution 
{

public:

SlabsExactSolution( );

~SlabsExactSolution();

// int evalTest();

// Check that the solution satisfies Maxwell equations and the interface conditions **finish me ***
int check();


// evaluate the solution
int eval(real t, CompositeGrid & cg, int grid, 
	 realArray & ua, realArray & pv,
	 const Index & I1a, const Index &I2a, const Index &I3a, 
	 int numberOfTimeDerivatives = 0,
         bool computeMagneticField = false );

// evaluate the solution in frequency space at a point x
// int eval( real x[3], real *Ev, real *Hv = NULL );


int initialize( CompositeGrid & cg, int numberOfDomains,
		std::vector<DispersiveMaterialParameters> & dispersiveMaterialParameters,
		const real & omega, const RealArray & kvI, const int solveForAllFields );

// There are four scattering cases, 2 polarizations, forward/backward 
int setScatteringCase( int scatCase );



private:

  // The database is a place to store parameters
  mutable DataBase dbase;

};


#endif
