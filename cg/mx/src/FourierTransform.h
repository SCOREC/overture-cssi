#ifndef FOURIER_TRANSFORM
#define FOURIER_TRANSFORM

#include "Overture.h"

#define KK_DEBUG
#include "DBase.hh"
using namespace DBase;

// forward declaration:
class IndexBox;

class FourierTransform
{
public:

FourierTransform();
~FourierTransform();

int initialize( realMappedGridFunction & u, const int ndfft=1, const int sidefft=-1, const int axisfft=-1 );

int initialize( realArray & u, const IntegerArray & indexRange, const IntegerArray & gridIndexRange,
		const int ndfft=1, const int sidefft=-1, const int axisfft=-1 );

// int initialize( realArray & x, const int numberOfDimensions=1, const int dir1=0, const int dir2=1, const int dir3=2 );

// forward transform 
int forwardTransform( const RealArray & u, const int mc, void *zl  );

// backward transform
int backwardTransform( void *zl, RealArray & u, const int mc  );

int getLocalIndexBox( IndexBox & fftBox ) const;

// Update periodic boundaries
int periodicUpdate( const IntegerArray & gridIndexRange, const IntegerArray & dimension,
		    const IntegerArray & indexRange, const IntegerArray & isPeriodic,
		    realArray & x, RealArray & u, const int mc  );

// Update periodic boundaries
int periodicUpdate( realMappedGridFunction & u, const Range & C = nullRange );

// Update periodic boundaries - this version takes a local array
int periodicUpdate( RealArray & u, const Range & C = nullRange );

// set the debug flag, non-zero to output info  
int setDebug( int value );

// set the name of the debug file
int setDebugFileName( const aString & fileName );

// Physical transforms are the reverse of the default DFTs
int usePhysicalTransforms( bool usePhysicalTransforms=true );

// CPU time used by all calls
static real cpuTime;

protected:

 // The database is used to hold parameters
 mutable DataBase dbase;

};





#endif
