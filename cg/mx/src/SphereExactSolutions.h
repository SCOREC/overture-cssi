#ifndef SPHERE_EXACT_SOLUTIONS
#define SPHERE_EXACT_SOLUTIONS


// ===============================================================================
// Class to define exact solutions to Maxwell's equations for a sphere
//     Scattering from a PEC sphere
//     Scattering from a dieletric sphere
// ===============================================================================

#include "Maxwell.h"

#define KK_DEBUG
#include "DBase.hh"
using namespace DBase;


class SphereExactSolutions 
{

public:

SphereExactSolutions();

~SphereExactSolutions();

// Check that the solution satisfiesd Maxwell equations: 
int check();

int initialize( int numberOfDomains, real a, real *sc, real *kc, 
                real *epsOut, real *muOut, real *epsIn=NULL, real *muIn=NULL );


// evaluate the solution in frequency space at a point x
int eval( real x[3], real *Ev, real *Hv=NULL );

// evaluate the solution in frequency space on a grid 
int eval( realMappedGridFunction & u , int domain, int grid, bool computeMagneticField=false );

// Evaluate some derived quantities
int getDispersiveParameters( real beta0v[2], real beta1v[2] );


private:

 int setPars();

  // The database is a place to store parameters
  mutable DataBase dbase;

};


#endif
