#ifndef NONLINEAR_EXACT_SOLUTIONS
#define NONLINEAR_EXACT_SOLUTIONS


// ===============================================================================
// Class to define exact solutions for nonlinear models
// 
//   (1) Solution to 1D Maxwell-Bloch equations
//
// ===============================================================================

#include "Maxwell.h"

#define KK_DEBUG
#include "DBase.hh"
using namespace DBase;

// class MxParameters;

// typedef ::real LocalReal;

class NonlinearExactSolution 
{

public:

NonlinearExactSolution( );

~NonlinearExactSolution();


int check();


// evaluate the solution
int eval(real dt, real t, CompositeGrid & cg, int grid, 
	 realArray & ua, realArray & pv, realArray & qv,
	 const Index & I1a, const Index &I2a, const Index &I3a, 
	 int numberOfTimeDerivatives = 0,
         bool computeMagneticField = false );

// evaluate the solution in frequency space at a point x
// int eval( real x[3], real *Ev, real *Hv = NULL );


int initialize( CompositeGrid & cg, int numberOfDomains,
		std::vector<DispersiveMaterialParameters> & dispersiveMaterialParameters,
		const real & omega, const RealArray & kvI, const RealArray & asymParams, const int solveForAllFields );

private:

  // The database is a place to store parameters
  mutable DataBase dbase;

};


#endif
