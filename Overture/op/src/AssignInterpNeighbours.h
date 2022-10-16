#ifndef ASSIGN_INTERP_NEIGHBOURS
#define ASSIGN_INTERP_NEIGHBOURS "AssignInterpNeighbours"

// ==========================================================================
//  This class is used to assign the unused points next to interpolation
// points so that a wide (e.g. 5-pt) stencil can be used with fewer layers
// of interpolation points. This assignment is currently done with extrapolation
// although in the future we could consider interpolating these points instead.
// ==========================================================================

#include "Overture.h"

// forward declarations
class InterpolatePointsOnAGrid;
class OGFunction;

#define KK_DEBUG
#include "DBase.hh"
using namespace DBase;

class AssignInterpNeighbours
{
public:

  enum AssignmentTypeEnum
  {
    extrapolateInterpolationNeighbours=0,
    interpolateInterpolationNeighbours
  };

AssignInterpNeighbours( AssignmentTypeEnum assignmentType = extrapolateInterpolationNeighbours );

// copy constructor
AssignInterpNeighbours( const AssignInterpNeighbours & x );

~AssignInterpNeighbours();

// Assign values to the unused points next to interpolation points
int 
assign( realMappedGridFunction & uA, Range & C, const BoundaryConditionParameters & bcParams );

// Call this routine when the grid has changed and we need to re-initialize
int gridHasChanged();

// Assign values to the unused points next to interpolation points
int assignInterpolationNeighbours( realCompositeGridFunction & u, const Range & C = nullRange, OGFunction *TZFlow=NULL, real t=0. );

AssignInterpNeighbours & operator= ( const AssignInterpNeighbours & x );

void setAssignmentType( AssignmentTypeEnum assignmentType );

// Provide the interpolation point array (used in serial only)
void setInterpolationPoint( intArray & interpolationPoint );

// For interpolating interp neighbours:
int setInterpolationWidth( int width );

// For interpolating interp neighbours:
int setNumberOfValidGhostPoints( int numValidGhost );

// return size of this object  
real sizeOf(FILE *file = NULL ) const;


static int debug;

protected:

  // This next enum is used for setting an error status
  enum ErrorStatusEnum
  {
    noErrors=0,
    errorInFindInterpolationNeighbours
  } errorStatus;


// setup routine
int 
setup();

// routine for setting up arrays for assigning the neighbours to interpolation points
int 
findInterpolationNeighbours( MappedGrid & mg );

// For interpolating interp neighbours:
int setupInterpolation( CompositeGrid& cg );


int isInitialized;
AssignmentTypeEnum assignmentType;

int numberOfInterpolationNeighbours;

int maximumWidthToExtrapolationInterpolationNeighbours; 
IntegerArray *extrapolateInterpolationNeighbourPoints;
IntegerArray *extrapolateInterpolationNeighboursDirection;
IntegerArray *extrapolateInterpolationNeighboursVariableWidth;  

intArray *interpolationPoint;

// For interpolating interp-neighbours June 2022
int ipogIsInitialized;
InterpolatePointsOnAGrid *ipog;  // new way for parallel 
IntegerArray periodicUpdateNeeded;


static FILE *debugFile;  // make one debug file for all instances (we use the same name)

// for communicating values:

int npr, nps;     // number of proc. that we receieve or send data to
int *ppr,*pps;    // list of processors for rec. and sending to
int *nar, **iar;  // list of points to recieve from other processors
int *nas, **ias;  // list of points to send to other processor


  #ifdef USE_PPP
    MPI_Comm AIN_COMM;  // Communicator for the parallel interpolator
  #else
    int AIN_COMM;
  #endif

 // This database contains parameters and data *new way* June 2022.
 DataBase dbase;    

};


#endif
