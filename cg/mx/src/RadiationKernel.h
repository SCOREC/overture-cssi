#ifndef RADIATION_KERNEL
#define RADIATION_KERNEL


#include "Overture.h"

#define KK_DEBUG
#include "DBase.hh"
using namespace DBase;

class RadiationKernel
{
public:

enum KernelTypeEnum
{
  planar,
  slab,
  cylindrical,
  spherical
};


RadiationKernel();
~RadiationKernel();

int setKernelType( KernelTypeEnum type );
KernelTypeEnum getKernelType() const;

// New initialize routine for parallel, u provides the grid and parallel distribution
int initialize( realMappedGridFunction & u,
		int side, int axis,
		int numberOfFields_, 
		int numberOfModes1_,   int numberOfModes2_,
		real c_, 
		int orderOfTimeStepping_, int numberOfPoles_,
		real radius =1. );

  
// initialize for 3D domains
int initialize( int numberOfDimensions,
		int numberOfGridPoints1_, int numberOfGridPoints2_,
		int numberOfFields_, 
		int numberOfModes1_, int numberOfModes2_,
		real period1_, real period2_,
		real c_, 
		int orderOfTimeStepping_, int numberOfPoles_, 
                real radius=1. );

// initialize for 2D domains 
int initialize( int numberOfGridPoints1_, 
		int numberOfFields_, 
		int numberOfModes1_, 
		real period1_,  
		real c_, 
		int orderOfTimeStepping_, int numberOfPoles_,
		real radius =1. );
  
int evaluateKernel( double dt, RealArray & u, RealArray & Hu );

// a new serial version 
int evaluateKernelNew( double dt, RealArray & u, RealArray & Hu );

// latest parallel version 
int evaluateKernelParallel( double dt, RealArray & up, RealArray & hup );

// parallel version 3 
int evaluateKernelParallelV3( double dt, RealArray & up, RealArray & hup );

// parallel version 2 
int evaluateKernelParallelV2( double dt, RealArray & up, RealArray & hup );

// parallel version 1
int evaluateKernelParallelV1( double dt, RealArray & u, RealArray & Hu, RealArray & up, RealArray & hup );

// Specify whether to use the new FourierTransform class (Supports parallel FFTs)
int setUseFourierTransformClass( bool trueOrFalse = true );

// set the debug flag, non-zero to output info  
int setDebug( int value );

// set the name of the debug file
int setDebugFileName( const aString & fileName );

// CPU time used by all calls
static real cpuTime;

// -----------------------------------------------------------------------------
protected:

int initializePoles();

KernelTypeEnum kernelType;

int numberOfGridPoints1, numberOfGridPoints2;
int numberOfModes1, numberOfModes2;
int ns1, ns2;
int numberOfFields,numberOfPoles,orderOfTimeStepping,bcinit;
double c,period1,period2,radius;
double *ploc,*fold,*phi,*amc, *fftsave1, *fftsave2;

double *zl, *pl1, *zl2;  

double *alpha, *beta;
int *npoles;


// The database is the new place to store parameters
mutable DataBase dbase;

};

#endif
