//                                   -*- c++ -*-
#ifndef DEFORMING_BODY_MOTION_H
#define DEFORMING_BODY_MOTION_H
 
//
// Master class for keeping track of deforming bodies in flow:
//   The particular physics of the object is determined by
//   a separate class that handles the evolution of the body
//   under the surface stresses arising from the fluid:
//

#include "Overture.h"
#include "HyperbolicMapping.h"
#include "GenericGraphicsInterface.h"
#include "GridEvolution.h"
#include "DBase.hh"
using namespace DBase;

// Forward defs:
class GenericGraphicsInterface;
class MappingInformation;
class DeformingGrid;
class GridFunction; 
class Parameters;

// Forward defs for the physicsObjects
class ElasticFilament;
class BeamModel;
class NonlinearBeamModel;
class BodyForce;

//................................
class DeformingBodyMotion
{
 
public:

friend class Parameters;
friend class AdParameters;
friend class AsfParameters;
friend class CnsParameters;
friend class InsParameters;
friend class SmParameters;
friend class Cgins;

enum DeformingBodyType
{
  elasticFilament,
  //elasticRod,
  // elasticShell,
  elasticBody,
  userDefinedDeformingBody,
  unknownBody
};

// ---- "user" defined motions ----
enum UserDefinedDeformingBodyMotionEnum
{
  iceDeform,
  ellipseDeform,
  sphereDeform,
  freeSurface,    // previously advectBody
  elasticShell,
  elasticBeam,
  nonlinearBeam,
  interfaceDeform,
  userDeformingSurface,  // deforming surface is defined by the function: userDefinedDeformingSurface
  linearDeform, // for testing the elastic piston problem 
  cylDeform // beam cylinder
}; 
 
typedef ArraySimpleFixed<int,2,2,1,1> BcArray;
enum DeformingBoundaryBoundaryConditionEnum
{
  periodicBoundaryCondition=0,
  dirichletBoundaryCondition,
  neumannBoundaryCondition,
  slideBoundaryCondition,
  numberOfBoundaryConditions // counts number of extries
};

  

enum InitialStateOptionEnum
{
  initialPosition,
  initialVelocity,
  initialAcceleration,
  initialBoundaryPosition,
  initialBoundaryVelocity,
  initialBoundaryAcceleration 
};

// Starting curve (or surface) for deforming body can be one of 
enum StartCurveTypesEnum
{
  nurbsStartCurve=0,
  splineStartCurve
};

DeformingBodyMotion(Parameters & params, 
		    int numberOfTimeLevels               = 3,
		    GenericGraphicsInterface *pGIDebug00 =NULL,
		    int debug00                          =0);

~DeformingBodyMotion();

int advanceInterfaceDeform( real t1, real t2, real t3, 
			    GridFunction & cgf1,
			    GridFunction & cgf2,
			    GridFunction & cgf3,
			    int option  );

int buildBeamFluidInterfaceData( CompositeGrid & cg );
// Longfei 20170119: 
int buildBeamFluidInterfaceData3D( CompositeGrid & cg );

int buildElasticShellOptionsDialog(DialogData & dialog );

int buildElasticBeamOptionsDialog(DialogData & dialog );

int buildFreeSurfaceOptionsDialog(DialogData & dialog );

// apply correction at time t using new values of the forces at time t.
int correct( real t1, real t2, 
	     GridFunction & cgf1,GridFunction & cgf2 );

// define faces and grids that form the deforming body 
int defineBody( int numberOfFaces, IntegerArray & boundaryFaces );

const IntegerArray & getBoundaryFaces() const;

// Construct a grid from the past time, needed to start some PC schemes.
int getPastTimeGrid(  real pastTime , CompositeGrid & cg );

int getAccelerationBC( const real time0, const int grid, const int side, const int axis,
                       MappedGrid & mg, const Index &I1, const Index &I2, const Index &I3, 
		       realSerialArray & bcAcceleration);

// return the order of accuracy used to compute the acceleration 
int getAccelerationOrderOfAccuracy() const;

// return body as a BodyForce object for plotting 
int getBody( BodyForce & body );

int getBodyVolumeAndSurfaceArea( CompositeGrid & cg, real & volume, real & area );

// Return the beamModel (if it exists)
BeamModel& getBeamModel();

int getBulkSolidParameters( real & impedance );

int getElasticShellOption(const aString & answer, DialogData & dialog );

int getElasticBeamOption(const aString & answer, DialogData & dialog );

int getFreeSurfaceOption(const aString & answer, DialogData & dialog );

// return the initial state (position, velocity, acceleration)
int getInitialState( InitialStateOptionEnum stateOption, 
		     const real time, 
                     const int grid, MappedGrid & mg, const Index &I1, const Index &I2, const Index &I3, 
		     realSerialArray & state );

real getMaximumRelativeCorrection() const;

int getNumberOfGrids(); 

real getTimeStep() const;

int getVelocity( const real time0, 
		 const int grid, 
		 CompositeGrid & cg,
		 realArray & gridVelocity);


int getVelocityBC( const real time0, const int side, const int axis, const int grid, MappedGrid & mg, 
                   const Index &I1, const Index &I2, const Index &I3, 
		   realSerialArray & bcVelocity);

// return the order of accuracy used to compute the velocity 
int getVelocityOrderOfAccuracy() const;



//..Grid position, velocity & boundary acceleration

int initialize(CompositeGrid & cg, real t=0. );

int initializeGrid(CompositeGrid & cg, real t=0. ); 

int initializePast( real time00, real dt00, CompositeGrid & cg);

// return true if the deforming body is a beam model
bool isBeamModel() const;

// return true if the deforming body is a bulk solid model
bool isBulkSolidModel() const;

// return true if this is a beam model with fluid on two sides
bool beamModelHasFluidOnTwoSides() const;

// integrate the BODY to a new time
int integrate( real t1, real t2, real t3, 
	       GridFunction & cgf1,GridFunction & cgf2,GridFunction & cgf3,
	       realCompositeGridFunction & stress );

// output probe info
int outputProbes( GridFunction & gf0, int stepNumber );

// --- plot things related to moving grids (e.g. the center lines of beams or shells)
int plot(GenericGraphicsInterface & gi, GridFunction & cgf, GraphicsParameters & psp );

void printFilamentHyperbolicDimensions(CompositeGrid & cg00, int gridToMove00);

// print time step info
void printTimeStepInfo( FILE *file=stdout );

// Project the interface velocity (for added mass schemes)
int projectInterfaceVelocity( GridFunction & cgf );

int regenerateComponentGrids( const real newT, CompositeGrid & cg);

void registerDeformingComponentGrid( const int grid, CompositeGrid & cg);

// Longfei 20170220:
int saveShow();

// set the order of accuracy used to compute the acceleration
int setAccelerationOrderOfAccuracy( int order );

// set the order of accuracy used to compute the velocity
int setVelocityOrderOfAccuracy( int order );

int setType( const DeformingBodyType bodyType );

// interactive update
int update(CompositeGrid & cg, GenericGraphicsInterface & gi );
int update(  GenericGraphicsInterface & gi );


// user defined deforming surface: 
int userDefinedDeformingSurface(real t1, real t2, real t3, 
				GridFunction & cgf1,
				GridFunction & cgf2,
				GridFunction & cgf3,
				int option );

bool hasCorrectionConverged() const;

int get( const GenericDataBase & dir, const aString & name);

int put( GenericDataBase & dir, const aString & name) const;

// Write information to the `check file' 
int writeCheckFile( real t, FILE *file );

// Write information about the moving grids
void writeParameterSummary( FILE *file= stdout );

protected: 

int advanceElasticShell(real t1, real t2, real t3, 
			GridFunction & cgf1,
			GridFunction & cgf2,
			GridFunction & cgf3,
			realCompositeGridFunction & stress,
			int advanceOption );

int advanceElasticBeam(real t1, real t2, real t3, 
		       GridFunction & cgf1,
		       GridFunction & cgf2,
		       GridFunction & cgf3,
		       realCompositeGridFunction & stress,
		       int advanceOption );

int advanceFreeSurface(real t1, real t2, real t3, 
                       GridFunction & cgf1,
                       GridFunction & cgf2,
                       GridFunction & cgf3,
                       int advanceOption );

int advanceNonlinearBeam(real t1, real t2, real t3, 
			 GridFunction & cgf3,
			 realCompositeGridFunction & stress,
			 int advanceOption );

  int getFace( int grid ) const;

  int getPastLevelGrid( const int level, 
		  const int grid, 
		  CompositeGrid & cg,
		  realArray & gridVelocity);

  void simpleGetVelocity( const real vTime, 
			  const int grid00, 
			  CompositeGrid & cg,
			  realArray & gridVelocity);

int smoothInterface( RealArray & x0, CompositeGrid & cg, 
                     const int sideToMove, const int axisToMove, const int gridToMove, const real *vScale );

  ElasticFilament *pElasticFilament;  // the Mapping for an elastic filament 
  DeformingGrid *pDeformingGrid;      // The "grid" associated with the elastic filament

  BeamModel* pBeamModel;

  NonlinearBeamModel* pNonlinearBeamModel;

  int debug;
  GenericGraphicsInterface *pGIDebug;
  MappingInformation       *pMapInfoDebug;

  Parameters & parameters;  // parameters from the DomainSolver

  mutable DataBase deformingBodyDataBase;  // save DeformingBodyMotion parameters in here 


};

#endif

