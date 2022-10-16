#include "DomainSolver.h"

#ifndef CGAD_H
#define CGAD_H


// CG equation-domain solver for the advection-diffusion (AD) equations

class Cgad : public DomainSolver
{
public:

Cgad(CompositeGrid & cg, GenericGraphicsInterface *ps=NULL, Ogshow *show=NULL, const int & plotOption=1 );


virtual ~Cgad();


virtual
void addForcing(realMappedGridFunction & dvdt, const realMappedGridFunction & u, int iparam[], real rparam[],
		realMappedGridFunction & dvdtImplicit = Overture::nullRealMappedGridFunction(),
                realMappedGridFunction *referenceFrameVelocity=NULL);

virtual void
advanceADI( real & t, real & dt, int & numberOfSubSteps, int & init, int initialStep  );

virtual int 
applyBoundaryConditions(const real & t, realMappedGridFunction & u, 
			realMappedGridFunction & gridVelocity,
			const int & grid,
			const int & option=-1,
			realMappedGridFunction *puOld=NULL, 
			realMappedGridFunction *pGridVelocityOld=NULL,
			const real & dt=-1.);


virtual int
applyBoundaryConditionsForImplicitTimeStepping(realMappedGridFunction & rhs, 
					       realMappedGridFunction & uL,
                                               realMappedGridFunction &uOld, // *wdh* Dec 20, 2017 -- added uOld
					       realMappedGridFunction & gridVelocity,
					       real t,
					       int scalarSystem,
					       int grid );

virtual void 
buildImplicitSolvers(CompositeGrid & cg);

virtual int 
buildRunTimeDialog();

virtual int 
buildTimeSteppingDialog(DialogData & dialog );

virtual void 
formMatrixForImplicitSolve(const real & dt0,
			   GridFunction & cgf1,
			   GridFunction & cgf0 );

// virtual int
// iterativeInterfaceRightHandSide( IterativeInterfaceOptionsEnum option, GridFaceDescriptor & info, 
//                                 int gfIndex, real t );

// virtual int 
// formImplicitTimeSteppingMatrix(realMappedGridFunction & coeff,
// 			       const real & dt0, 
// 			       int scalarSystem, 
// 			       realMappedGridFunction & uL,
// 			       const int & grid );

  
// Evaluate variable advection coefficients:
virtual void 
getAdvectionCoefficients( GridFunction & cgf );

// setup champ boundary conditions
int champBoundaryConditions( realCompositeGridFunction & coeffcg, Parameters & parameters, Real dt );

// Evaluate variable diffusion coefficients:
virtual void 
getDiffusionCoefficients( GridFunction & cgf );


// Return the list of interface data needed by a given interface:
virtual int
getInterfaceDataOptions( GridFaceDescriptor & info, int & interfaceDataOptions ) const;

// Optionally evaluate the max residual in any interface conditions (only used for some interface conditions such as CHAMP)
virtual int
getInterfaceResidual( real t, real dt, GridFunction & cgf, Real & residual );

virtual int
getUt(const realMappedGridFunction & v, 
	  const realMappedGridFunction & gridVelocity, 
	  realMappedGridFunction & dvdt, 
	  int iparam[], real rparam[],
	  realMappedGridFunction & dvdtImplicit = Overture::nullRealMappedGridFunction(),
	  MappedGrid *pmg2=NULL,
	  const realMappedGridFunction *pGridVelocity2= NULL);


virtual void
getTimeSteppingEigenvalue(MappedGrid & mg, 
 			       realMappedGridFunction & u, 
 			       realMappedGridFunction & gridVelocity,  
 			       real & reLambda,
 			       real & imLambda, 
 			       const int & grid);

virtual int
getTimeSteppingOption(const aString & answer,
		      DialogData & dialog );

virtual void 
implicitSolve(const real & dt0,
	      GridFunction & cgf1,
              GridFunction & cgf0);

// routine called by takeTimeStepBDF
virtual int 
implicitTimeStep( real & t0, real & dt0, int correction, AdvanceOptions & advanceOptions );

virtual int
interfaceRightHandSide( InterfaceOptionsEnum option, 
                        int interfaceDataOptions,
                        GridFaceDescriptor & info, 
                        GridFaceDescriptor & gfd, 
			int gfIndex, real t, bool saveTimeHistory = false );

virtual int 
plot(const real & t, 
     const int & optionIn,
     real & tFinal,
     int solutionToPlot = -1 );

// -- project interface values such as Temperature for conjugate heat transfer  ---
virtual int
projectInterface(real t, real dt, GridFunction & cgf );


virtual void
saveShowFileComments( Ogshow &show );


virtual int
setOgesBoundaryConditions( GridFunction &cgf, IntegerArray & boundaryConditions, RealArray &boundaryConditionData,
                           const int imp );

virtual int
setPlotTitle(const real &t, const real &dt);

virtual int
setupGridFunctions();

virtual int 
setupPde(aString & reactionName,bool restartChosen, IntegerArray & originalBoundaryCondition);

virtual int 
setupUserDefinedInitialConditions();

int 
thinFilmSolver(  real & t0, real & dt0, int correction, AdvanceOptions & advanceOptions );

virtual int
updateGeometryArrays(GridFunction & cgf);

virtual int
updateToMatchGrid(CompositeGrid & cg);

virtual int 
updateStateVariables(GridFunction & cgf, int stage=-1);

virtual int
userDefinedBoundaryValues(const real & t, 
                          GridFunction & gf0,
			  const int & grid,
			  int side0 = -1,
			  int axis0 = -1,
			  ForcingTypeEnum forcingType =computeForcing );

virtual int 
userDefinedInitialConditions(CompositeGrid & cg, realCompositeGridFunction & u );

virtual void 
userDefinedInitialConditionsCleanup();

virtual
void writeParameterSummary( FILE *file );

protected:


private:


};

#endif
