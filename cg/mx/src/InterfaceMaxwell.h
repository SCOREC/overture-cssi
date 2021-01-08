#ifndef INTERFACE_MAXWELL_H
#define INTERFACE_MAXWELL_H

// ======================================================================================
// Class to assign interface conditions for Maxwell's equations
// ======================================================================================


#include "Maxwell.h"

#define KK_DEBUG
#include "DBase.hh"
using namespace DBase;


class InterfaceMaxwell
{

public:

InterfaceMaxwell( Maxwell & cgmx );
~InterfaceMaxwell();


// void initializeInterfaces( Maxwell & mx );

// -- compute the interface stencil coefficients ---
int computeStencilCoefficients( int current, real t, real dt );


private:

  Maxwell & mx;  // Here is a reference to the Maxwell solver 

  // The database is a place to store parameters
  mutable DataBase dbase;


};


#endif
