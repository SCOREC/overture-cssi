#ifndef RADIATION_BOUNDARY_CONDITION
#define RADIATION_BOUNDARY_CONDITION


#include "Overture.h"

class RadiationKernel;
class OGFunction;

#define KK_DEBUG
#include "DBase.hh"
using namespace DBase;

class RadiationBoundaryCondition
{
public:

RadiationBoundaryCondition(int orderOfAccuracy=4);
~RadiationBoundaryCondition();


int initialize( realMappedGridFunction & u, 
		int side, int axis,
                int nc1=0, int nc2=0, 
		real c_=1., real period_=-1.,  
		int numberOfModes_=-1, 
		int orderOfTimeStepping_=-1, int numberOfPoles_=-1 );

int assignBoundaryConditions( realMappedGridFunction & u, real t, real dt,
			   realMappedGridFunction & u2 );

int setDebug( int debugFlag, FILE *pDebugFile=stdout );

int setOrderOfAccuracy( int orderOfAccuracy );
int setNumberOfPoles( int numPoles  );

// use new parallel version
int useParallelVersion( bool trueOrFalse );

static int debug;
static real cpuTime;

OGFunction *tz;

protected:

RadiationKernel *radiationKernel; 

int nc1,nc2,numberOfGridPoints1,numberOfGridPoints2, numberOfModes1,numberOfModes2, orderOfTimeStepping, numberOfPoles;
double period1, period2, c, radius; 

int rside,raxis;   // apply BC on this face

int orderOfAccuracy,numberOfDerivatives,numberOfTimeLevels,currentTimeLevel;
real currentTime;
RealArray uSave, uxSave;

// The database is the new place to store parameters
mutable DataBase dbase;

};

#endif
