#ifndef RIGID_BODY_MOTION_H
#define RIGID_BODY_MOTION_H

// Class for keeping track of the position and orientation of a rigid body moving
// under the influence of forces and torques.

#include "Overture.h"
#define KK_DEBUG
#include "DBase.hh"
using namespace DBase;


class GenericGraphicsInterface;
class Parameters;

class RigidBodyMotion
{
  
public:

enum PositionConstraintEnum
{
  positionHasNoConstraint=0,
  positionConstrainedToAPlane,
  positionConstrainedToALine,
  positionIsFixed
};

  enum RotationConstraintEnum
  {
    rotationHasNoConstraint=0,
    rotationConstrainedToAPlane,
    rotationConstrainedToALine,
    rotationIsFixed
  };

  enum TimeSteppingMethodEnum
  {
    leapFrogTrapezoidal=0,
    improvedEuler,
    implicitRungeKutta
  } timeSteppingMethod;

  enum TwilightZoneTypeEnum
  {
    polynomialTwilightZone=0,
    trigonometricTwilightZone
  };
  

  enum BodyForceTypeEnum
  {
    timePolynomialBodyForce=0,
    timeFunctionBodyForce,
    restrictAngle
  };

  enum
  {
    defaultOrderOfAccuracy=-1234
  };


  RigidBodyMotion(int numberOfDimensions = 3);
  ~RigidBodyMotion();
  
  bool axesOfInertiaHaveBeenInitialized() const;

  bool centerOfMassHasBeenInitialized() const;

  // correct solution at time t using new values of the forces at time t.
  int correct( real t, 
	       const RealArray & force, 
	       const RealArray & torque,
	       RealArray & xCM, 
	       RealArray & rotation );

  // Version of correct that takes added mass matrices.
  int correct( real t, 
	       const RealArray & force, 
	       const RealArray & torque,
               const RealArray & A11, const RealArray & A12, const RealArray & A21, const RealArray & A22, 
	       RealArray & xCM, 
	       RealArray & rotation );

  int displayAddedDampingTensors( const aString & label, const real t = 0., FILE *file=stdout ) const;

  // get from a data base file
  int get( const GenericDataBase & dir, const aString & name);

  int getAcceleration( real t, RealArray & aCM  ) const;
 
  // return the added-daping scale factor 
  real getAddedDampingScaleFactor() const;

  // Get the added damping tensors: 
  int getAddedDampingTensors( RealArray & addedDampingTensors, const real t ) const;

  // evaluate the added mass matrices at time t
  int getAddedMassMatrices( const real t, RealArray & A11 , RealArray & A12 , RealArray & A21, RealArray & A22 ) const;

  int getAngularAcceleration( real t, RealArray & omegaDot  ) const;

  int getAngularVelocities( real t, RealArray & omega  ) const;

  int getAxesOfInertia( real t, RealArray & axesOfInertia  ) const;

  int getBodyForces( const real t, RealArray & bodyForce, RealArray & bodyTorque ) const;

  bool getCorrectionHasConverged();

  // general purpose routine:
  int getCoordinates( real t, 
		      RealArray & xCM      = Overture::nullRealArray(), 
		      RealArray & vCM      = Overture::nullRealArray(),
		      RealArray & rotation = Overture::nullRealArray(), 
		      RealArray & omega    = Overture::nullRealArray(), 
		      RealArray & omegaDot = Overture::nullRealArray(), 
		      RealArray & aCM      = Overture::nullRealArray(),
                      RealArray & axesOfInertia = Overture::nullRealArray(),
                      RealArray & momentOfInertiaTensor = Overture::nullRealArray()
                     ) const;

  real getDensity() const; // return the density, if known. If not known return a negative value.

  // Return the value of the directProjectionAddedMass option 
  bool getDirectProjectionAddedMass() const;

  // get exact solution (e.g. for twilightzone forcing)
  int getExactSolution( const int deriv, const real t, RealArray & xe, RealArray & ve , RealArray & we ) const;

  // return the force and torque at a given time t
  int getForce(const real t, RealArray & fv, RealArray & gv ) const;

  int getInitialConditions(RealArray & x0 = Overture::nullRealArray() , 
			   RealArray & v0 = Overture::nullRealArray() , 
			   RealArray & w0 = Overture::nullRealArray() ,
			   RealArray & axesOfInertia = Overture::nullRealArray()  ) const ;

  // return the known solution (if available)
  int getKnownSolution( real t, 
			RealArray & xCM      = Overture::nullRealArray(), 
			RealArray & vCM      = Overture::nullRealArray(),
			RealArray & aCM      = Overture::nullRealArray(),
			RealArray & omega    = Overture::nullRealArray(), 
			RealArray & omegaDot = Overture::nullRealArray() ) const;

  real getMass() const;

  // return the "total force" on the body 
  int getMassTimesAcceleration( real t,
       		                RealArray & mvDot,
			        RealArray & mOmegaDot );

  real getMaximumRelativeCorrection();

  // return principal axes of inertia:
  RealArray getMomentsOfInertia() const;

  int getMomentOfInertiaTensor( real t, RealArray & momentOfInertiaTensor  ) const;

  int getPosition( real t, RealArray & xCM ) const;

  // determine the velocity or acceleration of a point p(t) belonging to the body.
  int getPointAcceleration( real t, const RealArray & p, RealArray & ap) const;
  int getPointVelocity( real t, const RealArray & p, RealArray & vp) const;


  // use this next routine to compute the point velocity or acceleration for many points.
  int getPointTransformationMatrix( real t, 
				    RealArray & rDotRt  = Overture::nullRealArray(),
				    RealArray & rDotDotRt = Overture::nullRealArray() ) const ;

  int getRotationMatrix(real t, 
                        RealArray & r,                     
                        RealArray & rDot = Overture::nullRealArray(),
                        RealArray & rDotDot = Overture::nullRealArray() ) const;
  
  real getTimeStepEstimate() const;

  // Return the name of the time stepping method as a string.
  aString getTimeSteppingMethodName() const;


  int getVelocity( real t, RealArray & vCM  ) const;

  // return the volume of the body     
  real getVolume() const;

  // Indicate whether added mass matrices will be used (and provided by the user)
  int includeAddedMass( bool trueOrFalse = true );

  // integrate to a new time
  int integrate( real tf, 
		 const RealArray & force, 
		 const RealArray & torque,
                 real t, 
                 RealArray & xCM, 
		 RealArray & rotation );

  // Version of integrate that takes added mass matrices.
  int integrate( real tf, 
		 const RealArray & force, 
		 const RealArray & torque,
                 real t, 
                 const RealArray & A11, const RealArray & A12, const RealArray & A21, const RealArray & A22, 
                 RealArray & xCM, 
		 RealArray & rotation );  

  bool massHasBeenInitialized() const;

  bool momentsOfInertiaHaveBeenInitialized() const;

  int momentumTransfer( real t, RealArray & v );

  // print time step info
  void printTimeStepInfo( FILE *file=stdout );

  // put to a data base file
  int put( GenericDataBase & dir, const aString & name) const;

  // Reset the rigid body to it's initial state. (e.g. remove saved solutions etc.)
  int reset();

  // Supply the acceleration to over-ride the default values (e.g. for AMP schemes)
  int setAcceleration( real t, RealArray & vDot, RealArray & wDot );

  // Specify the added damping tensors: 
  int setAddedDampingTensors( const RealArray & addedDampingTensors, const real t, const real scaleFactor );

  int setAxesOfInertial( const RealArray & axesOfInertia );

  // set the body number corresponding to the MovingGrids class
  int setBodyNumber( int bodyNumber );

  // set the name of the check file
  int setCheckFileName( const aString & checkFileName );

  void setDensity( const real bodyDensity );

  int setInitialCentreOfMass( const RealArray & x0 );

  int setInitialConditions( real t0=0., 
			    const RealArray & x0 = Overture::nullRealArray() , 
			    const RealArray & v0 = Overture::nullRealArray() , 
			    const RealArray & w0 = Overture::nullRealArray() ,
			    const RealArray & axesOfInertia = Overture::nullRealArray()  );

  void setMass( const real totalMass );

  int setMomentsOfInertia( const real mI1, const real mI2, const real mI3 );

  // Set the tolerance for the Newton iteration for implicit schemes.
  int setNewtonTolerance( real tol );

  // set the Parameters
  int setParameters( Parameters & parameters );

  // Supply forcing at "negative" times for startup (the fourth-order scheme requires values at t=-dt).
  int setPastTimeForcing( real t, const RealArray & force, const RealArray & torque );
  int setPastTimeForcing( real t, const RealArray & force, const RealArray & torque,
			const RealArray & A11, const RealArray & A12, const RealArray & A21, const RealArray & A22  );

  int setPositionConstraint( PositionConstraintEnum constraint, const RealArray & values );

  int setProperties(real totalMass, 
		    const RealArray & momentsOfInertia, 
		    const int numberOfDimensions = 3 );

  int setRotationConstraint( RotationConstraintEnum constraint, const RealArray & values );

  int setTwilightZone( bool trueOrFalse, TwilightZoneTypeEnum type = trigonometricTwilightZone, 
                       int degreeOfPoly = 2 );

  // choose the time stepping method
  int setTimeSteppingMethod( const TimeSteppingMethodEnum method, int orderOfAccuracy=defaultOrderOfAccuracy );

  void setVolume( real volume );

  int update( GenericGraphicsInterface & gi );

  // Return true if the added mass matrices are being used.
  bool useAddedMass() const;

  // Write the current solution (and errors) in the check file (for regression tests)
  int writeCheckFile( real t );

  // Write information about the moving grids
  void writeParameterSummary( FILE *file= stdout );


  // --- Martix utility routines (these should be put elsewhere) ---
  static RealArray mult( const RealArray & a, const RealArray & b );
  static RealArray trans( const RealArray &a );
  static real dot( const RealArray &a, const RealArray &b );
  static RealArray getCrossProductMatrix( const RealArray & w );
  static RealArray solve( const RealArray & a, const RealArray & b );

  static int debug;               // debug flag



 protected:


  enum ConstraintApplicationEnum
  {
    applyConstraintsToPosition=0,
    applyConstraintsToVelocity,
    applyConstraintsToForces
  };


  int applyConstraints( const int step, const ConstraintApplicationEnum option );

  int buildBodyForceOptionsDialog(DialogData & dialog );

  // Solve the DIRK equations
  int dirkImplicitSolve( const real dt, const real aii, const real tc, const RealArray & yv, const RealArray &yv0, 
			 RealArray & kv );

  int getBodyForceOption(const aString & answer, DialogData & dialog );

  // get forcing (protected routine) 
  int getForceInternal(const real t, RealArray & fv, RealArray & gv, RealArray *pA=NULL ) const;
 
  // Make a matrix orthonormal
  int orthoNormalizeMatrix( RealArray & e, int k =0 ) const;

  // Under-relax the forces on the body.
  int relaxForce( real t, int next, 
  		  const RealArray & f, const RealArray & force, const RealArray & bodyForce, 
		  const RealArray & g, const RealArray & torque, const RealArray & bodyTorque );

  // Take a step with implicit Runge-Kutta.
  int takeStepImplicitRungeKutta( const real t0, const real dt );

  // Take a leap frog step: (predictor)
  int takeStepLeapFrog( const real t0, const real dt );

  // Take a step with the trapezoidal rule (corrector)
  int takeStepTrapezoid( const real t0, const real dt );

  // Take a step of the trapezoidal rule (corrector)
  // int takeStepTrapezoid( const real t, const real dt );
 
  // Take a step with extrapolation (predictor)
  int takeStepExtrapolation( const real t0, const real dt );

  real mass;               // total mass
  int numberOfDimensions;  // 2 or 3 dimensional problem
  int current;             // position of current solution in arrays.
  int numberOfSteps;       // number of time steps taken
  int numberSaved;         // number of different times we have saved force and moment data.
  int maximumNumberToSave; // save at most this many previous time values of various quantities
  int initialConditionsGiven; // keeps track of which initial conditions have been given
  real density;            // In cases where the volume can be determined we can use the density to get the mass.

  real damping;  // for damping oscillations.

  RealArray mI;     // 3 moment of inertia

  RealArray e, e0;    // e(0:2,0:2) : 3 axes of inertia, e0=e at t=0

  RealArray x, v;  // centre of mass: position, velocity
  RealArray x0, v0, w0;  // initial centre of mass,  velocity and angular velocity
  
  RealArray f, g;  // force and torque

  RealArray w;      // angular velocities about axes of inertia

  RealArray time;   // arrays of times for which we know data.

  RealArray bodyForceCoeff, bodyTorqueCoeff; // coefficients in the body force

  const Range R;          // R=Range(0,2) 

  PositionConstraintEnum positionConstraint;
  RotationConstraintEnum rotationConstraint;
  RealArray constraintValues;
  
  // keep track of the convergence of the correction steps (used for light bodies)
  bool relaxCorrectionSteps;
  bool correctionHasConverged;
  real maximumRelativeCorrection;

  real correctionAbsoluteToleranceForce;
  real correctionRelativeToleranceForce;
  real correctionRelaxationParameterForce;

  real correctionAbsoluteToleranceTorque;
  real correctionRelativeToleranceTorque;
  real correctionRelaxationParameterTorque;

  bool twilightZone;
  TwilightZoneTypeEnum twilightZoneType;

  FILE *logFile;

  BodyForceTypeEnum bodyForceType;

  // The database is the new way to hold parameters.
  mutable DataBase dbase; 
};


#endif
