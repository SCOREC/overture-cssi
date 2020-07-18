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

int initialize( int numberOfDimensions,
		int numberOfGridPoints1_, int numberOfGridPoints2_,
		int numberOfFields_, 
		int numberOfModes1_, int numberOfModes2_,
		real period1_, real period2_,
		real c_, 
		int orderOfTimeStepping_, int numberOfPoles_, 
                real radius=1. );

int evaluateKernel( double dt, RealArray & u, RealArray & Hu );

static real cpuTime;

protected:

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
