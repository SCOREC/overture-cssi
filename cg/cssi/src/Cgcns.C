// ========================================================================================================
/// \class Cgcssi
/// \brief Solve the Compressible Navier-Stokes and Euler equations.
/// \details Cgcssi can be used to solve the compressible Navier-Stokes and reactive Euler equations on moving grids
///    with adaptive mesh refinement. 
///
// ========================================================================================================

#include "Cgcssi.h"
#include "CssiParameters.h"
#include "Ogshow.h"

// ===================================================================================================================
/// \brief Constructor for the Cgcssi class.
///
/// \param cg_ (input) : use this CompositeGrid.
/// \param ps (input) : pointer to a graphics object to use.
/// \param show (input) : pointer to a show file to use.
/// \param plotOption_ (input) : plot option 
/// 
///  \note CssiParameters (passed to the DomainSolver constructor above) replaces the base class Parameters
// ==================================================================================================================
Cgcssi::
Cgcssi(CompositeGrid & cg_, 
      GenericGraphicsInterface *ps /* =NULL */, 
      Ogshow *show /* =NULL */ , 
      const int & plotOption_ /* =1 */) 
   : DomainSolver(*(new CssiParameters),cg_,ps,show,plotOption_)
{
  className="Cgcssi";
  name="cssi";

  // should this be somewhere else? setup?
  if( realPartOfEigenvalue.size() != cg.numberOfComponentGrids() )
    realPartOfEigenvalue.resize(cg.numberOfComponentGrids(),-1.);
  if( imaginaryPartOfEigenvalue.size() != cg.numberOfComponentGrids() )
    imaginaryPartOfEigenvalue.resize(cg.numberOfComponentGrids(),-1.);
}


// ===================================================================================================================
/// \brief Destructor.
// ==================================================================================================================
Cgcssi::
~Cgcssi()
{
  delete & parameters;
}



// ===================================================================================================================
/// \brief Update geometry arrays when the grid has changed (called by adaptGrids for example).
/// \param cgf (input) : 
// ===================================================================================================================
int Cgcssi::
updateGeometryArrays(GridFunction & cgf)
{
  if( debug() & 8 ) printF(" --- Cgcssi::updateGeometryArrays ---\n");
  
  CompositeGrid & cg = cgf.cg;

  if( realPartOfEigenvalue.size() != cg.numberOfComponentGrids() )
  {
    realPartOfEigenvalue.resize(cg.numberOfComponentGrids(),-1.);
    //for( int grid=0; grid<realPartOfEigenvalue.size(); grid++ )
    //  printF("Cgcssi::updateGeometryArrays: realPartOfEigenvalue[%i]=%e\n",grid,realPartOfEigenvalue[grid]);
  }
  
  if( imaginaryPartOfEigenvalue.size() != cg.numberOfComponentGrids() )
    imaginaryPartOfEigenvalue.resize(cg.numberOfComponentGrids(),-1.);
  
  
  const CssiParameters::PDE & pde = parameters.dbase.get<CssiParameters::PDE >("pde");
  const CssiParameters::GodunovVariation & conservativeGodunovMethod = 
                           parameters.dbase.get<CssiParameters::GodunovVariation >("conservativeGodunovMethod");
  const CssiParameters::PDEVariation & pdeVariation = 
                 parameters.dbase.get<CssiParameters::PDEVariation >("pdeVariation");
  

  if( pde==CssiParameters::compressibleNavierStokes ||
      pde==CssiParameters::compressibleMultiphase )
  {
    if( pdeVariation==CssiParameters::conservativeWithArtificialDissipation )
    {

      cgf.cg.update(MappedGrid::THEinverseVertexDerivative | MappedGrid::THEinverseCenterDerivative | 
		    MappedGrid::THEcenterJacobian |
		    MappedGrid::THEvertex | MappedGrid::THEcenter );
    }
    else
    {
      int grid;
      for( grid=0; grid<cgf.cg.numberOfComponentGrids(); grid++ )
      {
	MappedGrid & mg = cgf.cg[grid];
	if( mg.isRectangular() )
	{
	  mg.update(MappedGrid::THEmask);
	}
	else
	{
	  mg.update(MappedGrid::THEinverseVertexDerivative | MappedGrid::THEinverseCenterDerivative | 
		    MappedGrid::THEcenterJacobian );
	}
      }
    }

  }

  return DomainSolver::updateGeometryArrays(cgf);
}

//===============================================================================================
/// \brief Update the solver to match a new grid.
/// \param cg (input): composite grid.
//===============================================================================================
int Cgcssi::
updateToMatchGrid(CompositeGrid & cg)
{
  // printF("\n $$$$$$$$$$$$$$$ Cgcssi: updateToMatchGrid(CompositeGrid & cg) $$$$$$$$$$$$\n\n");
  
  int returnValue =DomainSolver::updateToMatchGrid(cg);

  if( realPartOfEigenvalue.size() != cg.numberOfComponentGrids() )
    realPartOfEigenvalue.resize(cg.numberOfComponentGrids(),-1.);
  if( imaginaryPartOfEigenvalue.size() != cg.numberOfComponentGrids() )
    imaginaryPartOfEigenvalue.resize(cg.numberOfComponentGrids(),-1.);

  return returnValue;
  
}

//===============================================================================================
/// \brief Project the initial conditions of a steady state newton solver on the linear solution.
/// \param cgf (input) : project this grid function.
/// \author kkc. 
//===============================================================================================
int Cgcssi::
project(GridFunction & cgf)
{ // project the initial conditions of a steady state newton solver on the linear solution
  real time0=getCPU();
  const int myid=Communication_Manager::My_Process_Number;

  if ( parameters.dbase.get<Parameters::TimeSteppingMethod >("timeSteppingMethod")!=Parameters::steadyStateNewton)
    return 0;

  if( parameters.dbase.get<bool >("projectInitialConditions") )
  {

    updateToMatchGrid(gf[current].cg); 
      
    int nss = 1,s=0,init=0;
    real dt=0,t=0;
    real impFacSave = parameters.dbase.get<real >("implicitFactor");
    parameters.dbase.get<real >("implicitFactor") = 1.;
    DomainSolver::advanceNewton(t,dt,nss,init,s);
    parameters.dbase.get<int >("globalStepNumber")--;
    parameters.dbase.get<real >("implicitFactor") = impFacSave;
  }
  else
    updateToMatchGrid(gf[current].cg); 

  if( debug() & 16 ) 
  {
    fprintf(parameters.dbase.get<FILE* >("debugFile")," \n ****Solution after projectVelocity and BC's**** \n");
    outputSolution( cgf.u,0. );
  }

  real time=getCPU()-time0;
  printF(">>>>>Time to project = %8.2e s <<<<\n",time);

  return 0;
}

//===============================================================================================
/// \brief Save comments to the show file that will appear as the top label when viewed with plotStuff.
/// \param show (input) : show file to use.
//===============================================================================================
void Cgcssi::
saveShowFileComments( Ogshow &show )
{
  char buffer[80]; 
  aString timeLine="";
  if(  parameters.dbase.has_key("timeLine") )
    timeLine=parameters.dbase.get<aString>("timeLine");

  const CssiParameters::PDE & pde = parameters.dbase.get<CssiParameters::PDE >("pde");
  const CssiParameters::GodunovVariation & conservativeGodunovMethod = 
                           parameters.dbase.get<CssiParameters::GodunovVariation >("conservativeGodunovMethod");
  const CssiParameters::PDEVariation & pdeVariation = 
                 parameters.dbase.get<CssiParameters::PDEVariation >("pdeVariation");
  const real mu = parameters.dbase.get<real >("mu");
  const real kThermal = parameters.dbase.get<real >("kThermal");

  if( pde==CssiParameters::compressibleNavierStokes )
  {
    aString showFileTitle[5];

    if( pdeVariation==CssiParameters::conservativeGodunov )
    {
      if( parameters.dbase.get<int >("numberOfSpecies")==0 )
      {
	if( mu>0 || kThermal>0 )
          showFileTitle[0]=sPrintF(buffer,"N-S mu=%7.1e k=%7.1e",mu,kThermal);
	else
	  showFileTitle[0]="Euler";
      }
      else 
	showFileTitle[0]="Reactive Euler (" + parameters.dbase.get<aString >("reactionName") + ")"; 
    }
    else if( parameters.dbase.get<real >("reynoldsNumber")<1.e20 )
    {
      showFileTitle[0]=sPrintF(buffer,"Compressible NS, Reynolds=%6.2e",parameters.dbase.get<real >("reynoldsNumber"));
    }
    else if( parameters.dbase.get<real >("mu")>0. || parameters.dbase.get<real >("kThermal")>0. )
    {
      showFileTitle[0]=sPrintF(buffer,"Compressible NS, mu=%3.1f, k=%3.1f",parameters.dbase.get<real >("mu"),
			       parameters.dbase.get<real >("kThermal"));
    }
    else
    {
      showFileTitle[0]=sPrintF(buffer,"Euler");
    }
    
    
    showFileTitle[1]=timeLine;
    showFileTitle[2]="";  // marks end of titles
    
    for( int i=0; showFileTitle[i]!=""; i++ )
    {
      // printf("** saveShow: showFileTitle[i]=[%s]\n",(const char*)showFileTitle[i]);
      show.saveComment(i,showFileTitle[i]);
    }
  }
  else if( pde==CssiParameters::compressibleMultiphase )
  {
    // save comments that go at the top of each plot

    aString showFileTitle[5];
    if( pdeVariation==CssiParameters::conservativeGodunov )
    {
      if( parameters.dbase.get<int >("numberOfSpecies")==0 )
      {
	showFileTitle[0]="Compressible Multiphase";
      }
      else 
	showFileTitle[0]="Reactive Compressible Multiphase (" + parameters.dbase.get<aString >("reactionName") + ")"; 
    }
    else if( parameters.dbase.get<real >("reynoldsNumber")<1.e20 )
    {
      showFileTitle[0]=sPrintF(buffer,"Compressible Multiphase, Reynolds=%6.2e",parameters.dbase.get<real >("reynoldsNumber"));
    }
    else if( parameters.dbase.get<real >("mu")>0. || parameters.dbase.get<real >("kThermal")>0. )
    {
      showFileTitle[0]=sPrintF(buffer,"Compressible Multiphase, mu=%3.1f, k=%3.1f",parameters.dbase.get<real >("mu"),
			       parameters.dbase.get<real >("kThermal"));
    }
    else
    {
      showFileTitle[0]=sPrintF(buffer,"Compressible Multiphase");
    }
    
    showFileTitle[1]=timeLine;
    showFileTitle[2]="";  // marks end of titles

    for( int i=0; showFileTitle[i]!=""; i++ )
      show.saveComment(i,showFileTitle[i]);

  }
}

//===============================================================================================
/// \brief Output parameter values to the header information that is printed near the start of the computation.
/// \param file (input) : write infomation to this file.
//===============================================================================================
void Cgcssi::
writeParameterSummary( FILE * file )
{
  DomainSolver::writeParameterSummary( file );

  FILE * checkFile = parameters.dbase.get<FILE* >("checkFile");
  
  const CssiParameters::PDE & pde = parameters.dbase.get<CssiParameters::PDE >("pde");
  const CssiParameters::GodunovVariation & conservativeGodunovMethod = 
                           parameters.dbase.get<CssiParameters::GodunovVariation >("conservativeGodunovMethod");
  const CssiParameters::PDEVariation & pdeVariation = 
                 parameters.dbase.get<CssiParameters::PDEVariation >("pdeVariation");

  if ( file==checkFile )
  {
    if( pde==CssiParameters::compressibleNavierStokes )
    {
      fPrintF(checkFile,"\\caption{Compressible Navier Stokes, gridName, $\\mu=%8.1e$, $t=%2.1f$, ",
	      parameters.dbase.get<real >("mu"),parameters.dbase.get<real >("tFinal"));
    }
    else if( pde==CssiParameters::compressibleMultiphase )
    {
      fPrintF(checkFile,"\\caption{Compressible Multiphase, gridName, $t=%2.1f$, ",
	      parameters.dbase.get<real >("tFinal"));
    }
    return;
  }

  fPrintF(file," machNumber=%e, reynoldsNumber=%e, prandtlNumber=%e \n",parameters.dbase.get<real >("machNumber"),
	  parameters.dbase.get<real >("reynoldsNumber"),parameters.dbase.get<real >("prandtlNumber"));
  fPrintF(file," mu=%e, kThermal=%e, thermalConductivity=%e Rg=%e, gamma=%e \n",parameters.dbase.get<real >("mu"),
          parameters.dbase.get<real >("kThermal"),parameters.dbase.get<real >("thermalConductivity"),
          parameters.dbase.get<real >("Rg"),
	  parameters.dbase.get<real >("gamma"));

  const ArraySimpleFixed<real,3,1,1,1> & gravity = parameters.dbase.get<ArraySimpleFixed<real,3,1,1,1> >("gravity");

  if( gravity[0]!=0. || 
      gravity[1]!=0. || 
      gravity[2]!=0. )
    fPrintF(file," gravity is on, acceleration due to gravity = (%8.2e,%8.2e,%8.2e) \n",
	    gravity[0],
            gravity[1],
            gravity[2]);
      
  if( pdeVariation==CssiParameters::conservativeWithArtificialDissipation )
  {
    fPrintF(file," conservative with artificial diffusion: av2=%7.3e, aw2=%7.3e, av4=%7.3e, aw4=%7.3e\n",
	    parameters.dbase.get<real >("av2"),parameters.dbase.get<real >("aw2"),parameters.dbase.get<real >("av4"),
            parameters.dbase.get<real >("aw4"));
  }
  else if( pdeVariation==CssiParameters::conservativeGodunov )
  {
    fPrintF(file," conservative Godunov method (variation=%s)\n",
	    conservativeGodunovMethod==CssiParameters::fortranVersion ? "fortran-version" :
	    conservativeGodunovMethod==CssiParameters::multiComponentVersion ? "multi-component" :
	    conservativeGodunovMethod==CssiParameters::multiFluidVersion ? "multi-fluid" :
	    "unknown" );
    fPrintF(file," Riemann solver = %s\n",parameters.dbase.get<CssiParameters::RiemannSolverEnum >("riemannSolver")==CssiParameters::exactRiemannSolver ? "exact" :
	    parameters.dbase.get<CssiParameters::RiemannSolverEnum >("riemannSolver")==CssiParameters::roeRiemannSolver ? "Roe" :
	    parameters.dbase.get<CssiParameters::RiemannSolverEnum >("riemannSolver")==CssiParameters::hLLRiemannSolver ? "HLL" : "unknown");
    
    fPrintF(file," order of accuracy of the Godunov method is %i.\n",
            parameters.dbase.get<int >("orderOfAccuracyForGodunovMethod"));
    int slopeLimiter=1; // default : use slope limiter
    parameters.dbase.get<ListOfShowFileParameters >("pdeParameters").getParameter("SlopeLimiter",slopeLimiter);
    fPrintF(file," slope limiter is %s.\n",(slopeLimiter==0 ? "off" : slopeLimiter==1 ? "on" : "invalid!"));
    fPrintF(file," Godunov artificial viscosity=%8.2e ( multiplies max(0,-div) ).\n",
            parameters.dbase.get<real >("godunovArtificialViscosity"));

  }
  else if( pdeVariation==CssiParameters::nonConservative )
  {
    fPrintF(file," nonconservative method: nuRho = %7.3e, anu=%e\n",parameters.dbase.get<real >("nuRho"),parameters.dbase.get<real >("anu"));
  }

  CssiParameters::EquationOfStateEnum & equationOfState = 
    parameters.dbase.get<CssiParameters::EquationOfStateEnum >("equationOfState");
  aString eosName;
  eosName=( equationOfState==CssiParameters::idealGasEOS ? "ideal gas" :
	    equationOfState==CssiParameters::jwlEOS ? "JWL" : 
	    equationOfState==CssiParameters::mieGruneisenEOS ? "Mie-Gruneisen" : 
	    equationOfState==CssiParameters::userDefinedEOS ? "user defined" : 
	    equationOfState==CssiParameters::stiffenedGasEOS ? "stiffened gas" : 
	    equationOfState==CssiParameters::taitEOS ? "Tait" : 
	    "unknown");

  if( equationOfState==CssiParameters::userDefinedEOS )
  {
    const aString & userDefinedEquationOfStateName = parameters.dbase.get<aString>("userDefinedEquationOfStateName");
    eosName = "user defined: " + userDefinedEquationOfStateName;
  }
  
  fPrintF(file," equation of state is: %s.\n",(const char*)eosName);
  
  fPrintF(file," slip-wall BC option=%i. 0=default, 1=slipWallPressureEntropySymmetry, 2=slipWallTaylor, 3=slipWallCharacteristic, 4=slipWallDerivative.\n",
          parameters.dbase.get<int >("slipWallBoundaryConditionOption"));
  
  fPrintF(file," densityLowerBound=%8.2e, pressureLowerBound=%8.2e, velocityLimiterEpsilon=%8.2e\n",
	  parameters.dbase.get<real >("densityLowerBound"),parameters.dbase.get<real >("pressureLowerBound"),parameters.dbase.get<real >("velocityLimiterEps"));

}

