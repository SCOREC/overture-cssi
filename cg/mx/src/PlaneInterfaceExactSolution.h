#ifndef PLANE_INTERFACE_EXACT_SOLUTIONS
#define PLANE_INTERFACE_EXACT_SOLUTIONS


// ===============================================================================
// Class to define exact solutions to Maxwell's equations
// 
//     Scattering from a plane material interface in 2D and 3D , dispersive or not
// ===============================================================================

#include "Maxwell.h"

#define KK_DEBUG
#include "DBase.hh"
using namespace DBase;

typedef ::real LocalReal;

class PlaneInterfaceExactSolution 
{

public:

PlaneInterfaceExactSolution();

~PlaneInterfaceExactSolution();

// Check that the solution satisfies Maxwell equations and the interface conditions
int check();


// evaluate the solution
int eval(real t, CompositeGrid & cg, int grid, 
	 realArray & ua, realArray & pv,
	 const Index & I1a, const Index &I2a, const Index &I3a, 
	 int numberOfTimeDerivatives = 0,
         bool computeMagneticField = false );

// evaluate the BA solution
int evalBA(DispersiveMaterialParameters & dmp1,
	   DispersiveMaterialParameters & dmp2,
	   real t, CompositeGrid & cg, int grid,
	   IntegerArray & matMask,
	   realArray & ua, realArray & pv,
	   const Index & I1a, const Index &I2a, const Index &I3a, 
	   int numberOfTimeDerivatives = 0,
	   int solveForAllFields = 1 );


// evaluate the solution in frequency space at a point x
int eval( real x[3], real *Ev, real *Hv = NULL );


// Evaluate some derived quantities
// int getDispersiveParameters( real beta0v[2], real beta1v[2] );

int initialize( CompositeGrid & cg, DispersiveMaterialParameters & dmp1, DispersiveMaterialParameters & dmp2,
		real *av, real *kvr, real *kvi );


  // could make this private:
int initializeBAPlaneInterfaceSolution( DispersiveMaterialParameters & dmp1,
					DispersiveMaterialParameters & dmp2, 
					const real kv[3],
					real & skr, real & ski );

private:

  // utility routine (using complex numbers) defined at the bottom of the file.
  void
  getTransmisionWaveNumber( const LocalReal & kr,  const LocalReal & ki, 
                            const LocalReal & kxr, const LocalReal & kxi, 
                            const LocalReal & kyr, const LocalReal & kyi, 
                            LocalReal & kxpr, LocalReal & kxpi, 
                            LocalReal & kypr, LocalReal & kypi );

  void
  checkPlaneMaterialInterfaceJumps( 
    const LocalReal & c1, const LocalReal & c2,
    const LocalReal & eps1, const LocalReal & eps2,
    const LocalReal & mu1, const LocalReal & mu2,

    const LocalReal & sr, const LocalReal & si,
    const LocalReal & rr, const LocalReal & ri, 
    const LocalReal & taur, const LocalReal & taui, 

    const LocalReal & eps1Hatr, const LocalReal & eps1Hati,
    const LocalReal & eps2Hatr, const LocalReal & eps2Hati,

    const LocalReal & psiSum1r, const LocalReal & psiSum1i,
    const LocalReal & psiSum2r, const LocalReal & psiSum12i,
    const LocalReal & kxr, const LocalReal & kxi,
    const LocalReal & kyr, const LocalReal & kyi,
    const LocalReal & kxpr, const LocalReal & kxpi,
    const LocalReal & kypr, const LocalReal & kypi
    );




  // The database is a place to store parameters
  mutable DataBase dbase;

};


#endif
