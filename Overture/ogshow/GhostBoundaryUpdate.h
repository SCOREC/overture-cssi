#ifndef GHOST_BOUNDARY_UPDATE_H
#define GHOST_BOUNDARY_UPDATE_H

// ==========================================================================
// GhostBoundaryUpdate: update for parallel ghost points and periodic points
// ==========================================================================

#include "Overture.h"

// Kyle Chand's nice data base class:
#define KK_DEBUG
#include "DBase.hh"
using namespace DBase;



class GhostBoundaryUpdate
{
public:

GhostBoundaryUpdate();
GhostBoundaryUpdate( const aString & debugFileName );
~GhostBoundaryUpdate();

// initialize communication schedules for parallel ghost updates
int initialize( realArray & u, const int numberOfGridDimensions );

// initialize communication schedules for parallel ghost updates and periodic updates (possibly on a face only)
int initialize( realMappedGridFunction & u,
 		int fside=-1, int faxis=-1 );

// initialize communication schedules for parallel ghost updates and periodic updates on a face
int initialize( realArray & u,
		const int numberOfGridDimensions,
		const IntegerArray & gid, const IntegerArray & dim, const IntegerArray & indexRange, const IntegerArray & isPeriodic,
		int fside=-1, int faxis=-1 );

// set the debug flag 
int setDebug( int debug );

// set the name of the debug file
int setDebugFileName( const aString & fileName );

// update parallel ghost and optionally periodic ghost 
#ifdef USE_PPP
int updateGhostBoundaries( realArray & u, const Range & C = nullRange );
#endif

// version to assign uLocal that may be different from the local array of u
int updateGhostBoundaries( realSerialArray & uLocal, const Range & C = nullRange );

protected:

void initializeClass(const aString & debugFileName);

// The database is a place to store parameters
mutable DataBase dbase;

};
#endif
