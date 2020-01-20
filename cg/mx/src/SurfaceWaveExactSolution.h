#ifndef SURFACE_WAVE_EXACT_SOLUTION
#define SURFACE_WAVE_EXACT_SOLUTION


// ==========================================================================================================
// Class to define exact solutions to Maxwell's equations
// 
//    SURFACE WAVE BETWEEN TWO MATERIALS
// 
//       - Surface Plasmon Polariton : surface wave between a metal (plasmon wave) and a dieletric (polariton wave)
// ==========================================================================================================

#include "Maxwell.h"

#define KK_DEBUG
#include "DBase.hh"
using namespace DBase;

typedef ::real LocalReal;

class SurfaceWaveExactSolution 
{

public:

SurfaceWaveExactSolution();

~SurfaceWaveExactSolution();

// evaluate the BA solution
int evalBA(DispersiveMaterialParameters & dmp1,
	   DispersiveMaterialParameters & dmp2,
	   real t, CompositeGrid & cg, int grid,
	   IntegerArray & matMask,
	   realArray & ua, realArray & pv,
	   const Index & I1a, const Index &I2a, const Index &I3a, 
	   int numberOfTimeDerivatives = 0,
	   int solveForAllFields = 1 );


int initialize( CompositeGrid & cg, DispersiveMaterialParameters & dmp1, DispersiveMaterialParameters & dmp2,
                const aString & caseName  );

private:

// The database is a place to store parameters
mutable DataBase dbase;

};


#endif
