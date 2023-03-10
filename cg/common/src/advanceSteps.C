// This file automatically generated from advanceSteps.bC with bpp.
// ==========================================================================================
//   This file contains functions that implement separate steps in an advance routine
//   These separate steps can be combined to form a time stepping algorithm such as 
//   a predictor corrector method.
//
// These functions should probably be virtual members of an Advance class so they can be 
// over-loaded? They now implement a PC method. 
// 
//      initializeTimeStepping( t,dt,init );
//      startTimeStep( t0,dt0,advanceOptions );
//      takeTimeStep( t0,dt0,correction,advanceOptions );
//      endTimeStep( t0,dt0,advanceOptions );
// 
// Here is the anticipated usage: 
//
//   initializeTimeStepping( t,dt,init )
//   for( int subStep=0; subStep<numberOfSubSteps; subStep++ )
//   {
//     startTimeStep( t0,dt0,advanceOptions );
//     for( int correction=0; correction<numberOfCorrections; correction++ )  // these could also be stages of a RK 
//     {    
//       takeTimeStep( t0,dt0,correction,advanceOptions );
//     }
//     endTimeStep( t0,dt0,advanceOptions );
// 
//   } // end  substeps
//
//
// ==========================================================================================


#include "DomainSolver.h"
#include "AdvanceOptions.h"

// ===================================================================================================================
/// \brief Initialize the time stepping (a time sub-step function). 
/// \details 
/// \param t0 (input) : current time
/// \param dt0 (input) : current time step
///
// ===================================================================================================================
int DomainSolver::
initializeTimeStepping( real & t0, real & dt0 )
{
    real cpu0=getCPU();
    int returnValue=0;
    const Parameters::TimeSteppingMethod & timeSteppingMethod = 
                                        parameters.dbase.get<Parameters::TimeSteppingMethod >("timeSteppingMethod");

    if( timeSteppingMethod==Parameters::forwardEuler )
    {
        returnValue=initializeTimeSteppingFE( t0,dt0 );
    }
    else if( timeSteppingMethod==Parameters::implicit )
    {
        const Parameters::ImplicitMethod & method = parameters.dbase.get<Parameters::ImplicitMethod>("implicitMethod");
        if ( method == Parameters::approximateFactorization )
            returnValue=initializeTimeSteppingAF( t0,dt0 );
        else if( method == Parameters::backwardDifferentiationFormula )
            returnValue=initializeTimeSteppingBDF( t0,dt0 );
        else 
            returnValue=initializeTimeSteppingIM( t0,dt0 );
    }
    else if( timeSteppingMethod==Parameters::adamsBashforth2 ||
                      timeSteppingMethod==Parameters::adamsPredictorCorrector2 ||
                      timeSteppingMethod==Parameters::adamsPredictorCorrector4 )
    {
        returnValue=initializeTimeSteppingPC( t0,dt0 );
    }
    else
    {
        printF("DomainSolver::initializeTimeStepping:ERROR: un-implemented time stepping=%i\n",
         	   (int)timeSteppingMethod);
        Overture::abort("error");
    }

  // Moved to multiStageAdvance.bC
  // // Initialize interface data such as time histories of interface heat fluxes -- *wdh* Dec 12, 2021
  // initializeInterfaceData(  t0, dt0 );

    if( parameters.dbase.get<int>("multiDomainProblem") )
    { // for multi-domain problems we need to increment the timings
        RealArray & timing = parameters.dbase.get<RealArray>("timing");
        timing(parameters.dbase.get<int>("timeForAdvance"))+=getCPU()-cpu0;
        timing(parameters.dbase.get<int>("totalTime"))+=getCPU()-cpu0;
    }

    if( false ) // *wdh* 2013/01/03 -- output is already called for globalStepNumber=-1 in advance
        output( gf[current],parameters.dbase.get<int >("globalStepNumber")+1 ); 

    return returnValue;
}


// ===================================================================================================================
/// \brief Start an individual time step (a time sub-step function).
/// \details 
/// \param t0 (input) : current time
/// \param dt0 (input) : current time step
/// \param correction (input) : for predictor corrector methods this indicates the correction step number.
/// \param currentGF (output) : points to the grid-function holding the current solution (time t0)
/// \param nextGF (output) : points to the grid-function holding the new solution (time t0+dt0)
/// \param numberOfCorrectorSteps (output) : return the number of corrector steps that will be used.
///
// ===================================================================================================================
int DomainSolver::
startTimeStep( real & t0, real & dt0, int & currentGF, int & nextGF, AdvanceOptions & advanceOptions )
{
    real cpu0=getCPU();
    int returnValue=0;
    const Parameters::TimeSteppingMethod & timeSteppingMethod = 
                                        parameters.dbase.get<Parameters::TimeSteppingMethod >("timeSteppingMethod");

    if( timeSteppingMethod==Parameters::forwardEuler )
    {
        returnValue=startTimeStepFE( t0,dt0,currentGF,nextGF,advanceOptions );
    }
    else if( timeSteppingMethod==Parameters::implicit )
    {
        const Parameters::ImplicitMethod &method = parameters.dbase.get<Parameters::ImplicitMethod>("implicitMethod");
        if ( method == Parameters::approximateFactorization )
            returnValue=startTimeStepAF( t0,dt0,currentGF,nextGF,advanceOptions );
        else if( method == Parameters::backwardDifferentiationFormula )
            returnValue=startTimeStepBDF( t0,dt0,currentGF,nextGF,advanceOptions );
        else
            returnValue=startTimeStepIM( t0,dt0,currentGF,nextGF,advanceOptions );
    }
    else if( timeSteppingMethod==Parameters::adamsBashforth2 ||
                      timeSteppingMethod==Parameters::adamsPredictorCorrector2 ||
                      timeSteppingMethod==Parameters::adamsPredictorCorrector4 )
    {
        returnValue=startTimeStepPC( t0,dt0,currentGF,nextGF,advanceOptions );
    }
    else
    {
        printF("DomainSolver::startTimeStep:ERROR: un-implemented time stepping=%i\n",
         	   (int)timeSteppingMethod);
        Overture::abort("error");
    }

  // save grid-function current and next: 
    parameters.dbase.get<int>("currentGF")=currentGF;
    parameters.dbase.get<int>("nextGF")   =nextGF;
    

    if( parameters.dbase.get<int>("multiDomainProblem") )
    { // for multi-domain problems we need to increment the timings
        RealArray & timing = parameters.dbase.get<RealArray>("timing");
        timing(parameters.dbase.get<int>("timeForAdvance"))+=getCPU()-cpu0;
        timing(parameters.dbase.get<int>("totalTime"))+=getCPU()-cpu0;
    }

    return returnValue;
}

// ===================================================================================================================
/// \brief Take a single time step (a time sub-step function).
/// \details 
/// \param t0 (input) : current time
/// \param dt0 (input) : current time step
/// \param correction (input) : for predictor corrector methods this indicates the correction step number.
/// \param advanceOptions (input) : additional options that adjust the behaviour of this function.
///       advanceOptions.takeTimeStepOption can be used to not apply or only apply the boundary conditions.
// ===================================================================================================================
int DomainSolver::
takeTimeStep( real & t0, real & dt0, int correction, AdvanceOptions & advanceOptions )
{
    real cpu0=getCPU();
    int returnValue=0;
    const Parameters::TimeSteppingMethod & timeSteppingMethod = 
                                        parameters.dbase.get<Parameters::TimeSteppingMethod >("timeSteppingMethod");

  // kkc 101006 added filter here, is this the right place? maybe it should be in startTimeStep instead?
    applyFilter(current); // note that this returns immediately if filtering is turned off

    if( timeSteppingMethod==Parameters::forwardEuler )
    {
        returnValue=takeTimeStepFE( t0,dt0,correction,advanceOptions );
    }
    else if( timeSteppingMethod==Parameters::implicit )
    {
        const Parameters::ImplicitMethod &method = parameters.dbase.get<Parameters::ImplicitMethod>("implicitMethod");
        if ( method == Parameters::approximateFactorization )
            returnValue=takeTimeStepAF( t0,dt0,correction,advanceOptions );
        else if ( method == Parameters::backwardDifferentiationFormula )
            returnValue=takeTimeStepBDF( t0,dt0,correction,advanceOptions );
        else
            returnValue=takeTimeStepIM( t0,dt0,correction,advanceOptions );
    }
    else if( timeSteppingMethod==Parameters::adamsBashforth2 ||
                      timeSteppingMethod==Parameters::adamsPredictorCorrector2 ||
                      timeSteppingMethod==Parameters::adamsPredictorCorrector4 )
    {
        returnValue=takeTimeStepPC( t0,dt0,correction,advanceOptions );
    }
    else
    {
        printF("DomainSolver::takeTimeStep:ERROR: un-implemented time stepping=%i\n",
         	   (int)timeSteppingMethod);
        Overture::abort("error");
    }

    if( parameters.dbase.get<int>("multiDomainProblem") )
    { // for multi-domain problems we need to increment the timings
        RealArray & timing = parameters.dbase.get<RealArray>("timing");
        timing(parameters.dbase.get<int>("timeForAdvance"))+=getCPU()-cpu0;
        timing(parameters.dbase.get<int>("totalTime"))+=getCPU()-cpu0;
    }
    return returnValue;
}

// ===================================================================================================================
/// \brief End an individual time step (a time sub-step function).
/// \details 
/// \param t0 (input) : current time
/// \param dt0 (input) : current time step
/// \param correction (input) : for predictor corrector methods this indicates the correction step number.
///
// ===================================================================================================================
int DomainSolver::
endTimeStep( real & t0, real & dt0, AdvanceOptions & advanceOptions )
{
    real cpu0=getCPU();
    int returnValue=0;
    const Parameters::TimeSteppingMethod & timeSteppingMethod = 
                                        parameters.dbase.get<Parameters::TimeSteppingMethod >("timeSteppingMethod");

    if( timeSteppingMethod==Parameters::forwardEuler )
    {
        returnValue=endTimeStepFE( t0,dt0,advanceOptions );
    }
    else if( timeSteppingMethod==Parameters::implicit )
    {
        const Parameters::ImplicitMethod &method = parameters.dbase.get<Parameters::ImplicitMethod>("implicitMethod");
        if ( method == Parameters::approximateFactorization )
            returnValue=endTimeStepAF( t0,dt0,advanceOptions );
        else if ( method == Parameters::backwardDifferentiationFormula )
            returnValue=endTimeStepBDF( t0,dt0,advanceOptions );
        else
            returnValue=endTimeStepIM( t0,dt0,advanceOptions );
    }
    else if( timeSteppingMethod==Parameters::adamsBashforth2 ||
                      timeSteppingMethod==Parameters::adamsPredictorCorrector2 ||
                      timeSteppingMethod==Parameters::adamsPredictorCorrector4 )
    {
        returnValue=endTimeStepPC( t0,dt0,advanceOptions );
    }
    else
    {
        printF("DomainSolver::endTimeStep:ERROR: un-implemented time stepping=%i\n",
         	   (int)timeSteppingMethod);
        Overture::abort("error");
    }

    if( parameters.dbase.get<int>("multiDomainProblem") )
    { // for multi-domain problems we need to increment the timings
        RealArray & timing = parameters.dbase.get<RealArray>("timing");
        timing(parameters.dbase.get<int>("timeForAdvance"))+=getCPU()-cpu0;
        timing(parameters.dbase.get<int>("totalTime"))+=getCPU()-cpu0;
    }
    return returnValue;
}


