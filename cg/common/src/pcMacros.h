// This file contains some macros that are shared amongst the different predictor-corrector methods


// ==================================================================================================
// MACRO: This macro saves past values of the pressure and values of the velocity on the ghost lines
// For use with the fourth-order accurate INS solver. 
//
// tp (input) : past value of time
// nab (input) : save results in fn[nab] (NOTE: use the grid from gf[mOld] not the one with fn[nab] !)
// ==================================================================================================
#beginMacro savePressureAndGhostVelocity(tp,nab)
if( orderOfAccuracy==4 )
{
  const int uc = parameters.dbase.get<int >("uc");
  const int pc = parameters.dbase.get<int >("pc");
  if( uc>=0 && pc>=0 ) // *wdh* April 16 2021 only do this for INS
  {
    OGFunction & e = *(parameters.dbase.get<OGFunction* >("exactSolution"));
         
    const int numberOfDimensions=cg.numberOfDimensions();
    const int numberOfGhostLines=2;
    Range V(uc,uc+numberOfDimensions-1);

    for( int grid=0; grid<gf[mOld].cg.numberOfComponentGrids(); grid++ )
    {
      MappedGrid & c = gf[mOld].cg[grid];

      realArray & fng = fn[nab][grid];
      realArray & uOld = gf[mOld].u[grid];
  #ifdef USE_PPP
      realSerialArray fnLocal; getLocalArrayWithGhostBoundaries(fng,fnLocal);
      realSerialArray uOldLocal; getLocalArrayWithGhostBoundaries(uOld,uOldLocal);
  #else
      realSerialArray & fnLocal = fng;
      realSerialArray & uOldLocal = uOld;
  #endif
      OV_GET_SERIAL_ARRAY_CONST(real,c.vertex(),xLocal);
      const int isRectangular=false; // for e.gd(..)

      const IntegerArray & gridIndexRange = c.gridIndexRange();
      getIndex(c.dimension(),I1,I2,I3);


      // save p for use when extrapolating in time
      //    ua(.,.,.,pc)= p(t-2*dt)  (for 2nd/4th order)
      //    ub(.,.,.,pc)= p(t-3*dt)  (for 4th order)
      //    uc(.,.,.,pc)= p(t-4*dt)  (for 4th order)
      if( parameters.dbase.get<bool >("twilightZoneFlow") )
      {
        // *wdh* 050416 fn[nab][grid](I1,I2,I3,pc)=e(c,I1,I2,I3,pc,tp);  
        //  fn[nab][grid](I1,I2,I3,pc)=e(c,I1,I2,I3,pc,tp);
        // e.gd(fn[nab][grid],0,0,0,0,I1,I2,I3,pc,tp);
        e.gd(fnLocal,xLocal,numberOfDimensions,isRectangular,0,0,0,0,I1,I2,I3,pc,tp);
        //  display(fn[nab][grid],"fn[nab][grid] after assigning for fourth order",debugFile,"%5.2f ");
        fprintf(debugFile,"savePressureAndGhostVelocity: Set p at old time for fourth-order: nab=%i, t=%9.3e\n",nab,tp);

        if( debug() & 4 )
        {
          display(xLocal,"savePressureAndGhostVelocity: xLocal from gf[mOld] ",debugFile,"%6.3f ");
          display(fn[nab][grid],"savePressureAndGhostVelocity: fn[nab][grid] after assigning p for fourth order",debugFile,"%6.3f ");
        }
            
          
      }
      else
      {
        printF("pcMacros: savePressureAndGhostVel: save past: TP=%9.3e, NAB=%i *** FIX ME ****\n",tp,nab);
        

        bool ok = ParallelUtility::getLocalArrayBounds(fng,fnLocal,I1,I2,I3);
        if( ok )
          fnLocal(I1,I2,I3,pc)=uOldLocal(I1,I2,I3,pc); // *** fix this ****
      }
        
      // We also extrapolate, in time, the ghost values of u -- used in the BC's
      getIndex(gridIndexRange,I1,I2,I3,numberOfGhostLines);
      for( int axis=0; axis<c.numberOfDimensions(); axis++ )
      {
        for( int side=0; side<=1; side++ )
        {
          const int is=1-2*side;
          if( c.boundaryCondition(side,axis)>0 )
          {
            // set values on the two ghost lines
            if( side==0 )
              Iv[axis]=Range(gridIndexRange(side,axis)-2,gridIndexRange(side,axis)-1);
            else
              Iv[axis]=Range(gridIndexRange(side,axis)+1,gridIndexRange(side,axis)+2);
   
            if( parameters.dbase.get<bool >("twilightZoneFlow") )
            {
              // *wdh* 050416 fn[nab][grid](I1,I2,I3,V)=e(c,I1,I2,I3,V,tp);
              // fn[nab][grid](I1,I2,I3,V)=e(c,I1,I2,I3,V,tp);
              // display(fn[nab][grid],"fn[nab][grid] before assign V on ghost",debugFile,"%5.2f ");
              // e.gd(fn[nab][grid],0,0,0,0,I1,I2,I3,V,tp);
              e.gd(fnLocal,xLocal,numberOfDimensions,isRectangular,0,0,0,0,I1,I2,I3,V,tp);
              // display(fn[nab][grid],"fn[nab][grid] after assign V on ghost",debugFile,"%5.2f ");

            }
            else
            {
              bool ok = ParallelUtility::getLocalArrayBounds(fng,fnLocal,I1,I2,I3);
              if( ok )
                fnLocal(I1,I2,I3,V)=uOldLocal(I1,I2,I3,V); // ***** fix this ****
            }
          }
        }
        // set back to gridIndexRange to avoid re-doing corners: *** is this ok for 3D ???
        Iv[axis]=Range(gridIndexRange(0,axis),gridIndexRange(1,axis));
      }
        
    }  // end for grid
  } 
} // end if orderOfAccuracy==4 
#endMacro



// ===============================================================================
///  MACRO:  Perform the initialization step for the PC method
///
///  \METHOD (input) : name of the method: e.g. adamsPC or implicitPC
///  Parameters:
///   numberOfPastTimes (input) : method needs u(t-dt), ... u(t-n*dt), n=numberOfPastTimes  
///   numberOfPastTimeDerivatives (input) : method needs u_t(t-dt), ..., u_t(t-m*dt) m=numberOfPastTimeDerivatives
///
// ===============================================================================
#beginMacro initializePredictorCorrector(METHOD,utImplicit)

const int orderOfPredictorCorrector = parameters.dbase.get<int >("orderOfPredictorCorrector");
const int orderOfTimeExtrapolationForPressure = parameters.dbase.get<int >("orderOfTimeExtrapolationForPressure");

printP("--METHOD-- initializePredictorCorrector: numberOfPastTimes=%i, numberOfPastTimeDerivatives=%i\n",numberOfPastTimes,numberOfPastTimeDerivatives);
printP("--METHOD-- initializePredictorCorrector: mCur=%i, mOld=%i gf[mCur].t=%9.2e\n",mCur,mOld,gf[mCur].t);
fPrintF(debugFile,"--METHOD-- initializePredictorCorrector: mCur=%i, mOld=%i gf[mCur].t=%9.2e\n",mCur,mOld,gf[mCur].t);

if( movingGridProblem() )
{ 
  getGridVelocity( gf[mCur],t0 );
}
 

if( orderOfTimeExtrapolationForPressure!=-1 )
{

  // if( orderOfPredictorCorrector==2 && orderOfTimeExtrapolationForPressure>1 &&
  //     poisson!=NULL && poisson->isSolverIterative()  )

  // *wdh* 2015/01/26: we may need past time pressure for other reasons:
  const bool & predictedPressureNeeded = parameters.dbase.get<bool>("predictedPressureNeeded");
  const bool predictPressure = predictedPressureNeeded || (poisson!=NULL && poisson->isSolverIterative());
  if( orderOfPredictorCorrector==2 && orderOfTimeExtrapolationForPressure>1 && predictPressure )
  {
    // orderOfTimeExtrapolationForPressure==1 :  p(t+dt) = 2*p(t) - p(t-dt)
    //                                      2 :  p(t+dt) = 3*p(t) - 3*p(t-dt) + p(t-2*dt)
    assert( previousPressure==NULL );
    assert( !parameters.isMovingGridProblem() );  // fix for this case
    
    numberOfExtraPressureTimeLevels = orderOfTimeExtrapolationForPressure - 1;
    printF("--DS-- ***initPC: allocate %i extra grid functions to store the pressure at previous times ****\n",
           numberOfExtraPressureTimeLevels);
    
    previousPressure = new realCompositeGridFunction [numberOfExtraPressureTimeLevels];
    for( int i=0; i<numberOfExtraPressureTimeLevels; i++ )
    {
      previousPressure[i].updateToMatchGrid(gf[mCur].cg);
    }
  }
  
  printF("--METHOD-- orderOfPredictorCorrector=%i, orderOfTimeExtrapolationForPressure=%i, predictPressure=%i\n",
         "           numberOfExtraPressureTimeLevels=%i\n",
         orderOfPredictorCorrector,orderOfTimeExtrapolationForPressure,(int)predictPressure,numberOfExtraPressureTimeLevels);
  

}

fn[nab0]=0.; 
if( numberOfPastTimeDerivatives>0 )
  fn[nab1]=0.; 
 
 
if( parameters.dbase.get<bool >("twilightZoneFlow") )
{
  OGFunction & e = *(parameters.dbase.get<OGFunction* >("exactSolution"));
 
  if( orderOfAccuracy==4 )
  {
    // For the fourth-order PC method, first compute u.t(t-2*dt) and u.t(t-3*dt)
    // Even for 2nd-order in time methods -- save p and u-ghost at t-2*dt
    const int numberOfPreviousValuesOfPressureToSave= orderOfPredictorCorrector==2 ? 1 : 2;

    int grid;
    for( int m=0; m<numberOfPreviousValuesOfPressureToSave; m++ )
    {
      // *** FIX ME: 4 -> numberOfExtraFunctionsToUse
      const int nab=(nab2+m) % 4; // save du/dt in fn[nab] 
      real tp=t0-(m+2)*dt0;       // move grid to this previous time

      if( movingGridProblem() )
      {
        // move gf[mOld] to t-(m+2)*dt
        moveGrids( t0,t0,tp,dt0,gf[mCur],gf[mCur],gf[mOld] );   // Is this correct? dt0?   
 
        gf[mOld].u.updateToMatchGrid(gf[mOld].cg); // *wdh* 040826

        // *wdh* 111125: the vertex is used below for error checking
        fn[nab].updateToMatchGrid(gf[mOld].cg);    


        // *wdh* 090806
        real cpu0=getCPU();
        gf[mOld].u.getOperators()->updateToMatchGrid(gf[mOld].cg); 
        parameters.dbase.get<RealArray>("timing")(parameters.dbase.get<int>("timeForUpdateOperators"))+=getCPU()-cpu0;


      }
      gf[mOld].t=tp;
 
      e.assignGridFunction( gf[mOld].u,tp );

      updateStateVariables(gf[mOld]); // *wdh* 080204 
      
      if( parameters.useConservativeVariables() )
        gf[mOld].primitiveToConservative();
 
      if( orderOfPredictorCorrector==4 ) 
      { // we only need du/dt at old times for pc4
        for( int grid=0; grid<gf[mCur].cg.numberOfComponentGrids(); grid++ )
        {
          rparam[0]=gf[mOld].t;
          rparam[1]=gf[mOld].t; // tforce
          rparam[2]=gf[mOld].t+dt0; //   ** check me **
          iparam[0]=grid;
          iparam[1]=gf[mOld].cg.refinementLevelNumber(grid);
          iparam[2]=numberOfStepsTaken;

          getUt(gf[mOld].u[grid],gf[mOld].getGridVelocity(grid),fn[nab][grid],iparam,rparam,
                utImplicit[grid],&gf[mOld].cg[grid]);
        }
      }
      
      // save past time values of p and ghost u for the 4th order method
      // NOTE: PAST time values are saved in a funny place:
      // save p for use when extrapolating in time
      //   ua(.,.,.,pc)= p(t-2*dt) : needed for 3rd order extrapolation: uCur(t), uOld(t-dt), ua(t-2*dt)
      //   ub(.,.,.,pc)= p(t-3*dt) : needed for 4th-order extrapolation: uCur(t), uOld(t-dt), ua(t-2*dt), ub(t-3*dt)
      //   uc(.,.,.,pc)= p(t-4*dt) : needed for 5th-order extrapolation: uCur(t), uOld(t-dt), ua(t-2*dt), ub(t-3*dt), ub(t-4*dt)
      assert( nab0==0 );
      const int nabPastTime=(nab0+m);
      savePressureAndGhostVelocity(tp,nabPastTime);
      
      if( debug() & 32 )
      {
        // determineErrors( gf[mOld].u,gf[mOld].gridVelocity, tp, 0, error,
        //               sPrintF(" adams:startup: errors in u at t=%e \n",nab,tp) );

        if( movingGridProblem() && debug() & 64 )
        {
          CompositeGrid & cg = *fn[nab].getCompositeGrid();
          
          for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
          {
            if( parameters.gridIsMoving(grid) )
            {
              display(cg[grid].vertex()(I1,I2,I3,0),sPrintF("\n *** PC: AFTER moveGrids:  fn[nab] "
                      "grid=%i vertex after move back t=%e",grid,gf[mOld].t),debugFile,"%8.5f ");
            }
          }
      
        }
        determineErrors( fn[nab],gf[mOld].gridVelocity, tp, 1, error,
                         sPrintF(" adams:startup: errors in ut (nab=%i) at t=%e \n",nab,tp) );
      }

    }
  }
 
       
  // get solution and derivative at t-dt
  if( movingGridProblem() )
  {
    // move gf[mOld] to t-dt
    if( debug() & 2 )
      fPrintF(debugFile,"--METHOD-- take an initial step backwards\n");
        
         // display(gf[mOld].cg[0].vertex()(I1,I2,I3,0),sPrintF(" gf[mOld] vertex before move back at t=%e",gf[mOld].t),
         //                  debugFile,"%5.2f ");
 
    moveGrids( t0,t0,t0-dt0,dt0,gf[mCur],gf[mCur],gf[mOld] );          // this will set gf[mOld].t=t-dt
 
    // *wdh* 090806
    if( parameters.isAdaptiveGridProblem() )
    { // both moving and AMR 
      parameters.dbase.get<Ogen* >("gridGenerator")->updateRefinement(gf[mOld].cg);
    }
    // *wdh* 090806
    real cpu0=getCPU();
    gf[mOld].u.getOperators()->updateToMatchGrid(gf[mOld].cg); 
    parameters.dbase.get<RealArray>("timing")(parameters.dbase.get<int>("timeForUpdateOperators"))+=getCPU()-cpu0;

    if( debug() & 64 )
    {
      for( int grid=0; grid<gf[mOld].cg.numberOfComponentGrids(); grid++ )
      {
        if( parameters.gridIsMoving(grid) )
        {
          display(gf[mOld].cg[grid].vertex()(I1,I2,I3,0),sPrintF("\n --METHOD-- AFTER moveGrids:  gf[mOld] grid=%i vertex after move back t=%e",grid,gf[mOld].t),
                  debugFile,"%10.7f ");
        }
      }
      
    }
    
 
    gf[mOld].u.updateToMatchGrid(gf[mOld].cg); // make sure the grid is correct, vertex used in TZ  *wdh* 040826

    // *wdh* 111125: the vertex is used below for error checking and computing ghost values of u
    if( numberOfPastTimeDerivatives>0 )
      fn[nab1].updateToMatchGrid(gf[mOld].cg);

  }
  else
    gf[mOld].t=t0-dt0; 
 
  // assign u(t-dt) with the TZ solution: 
  if( debug() & 1 )
    printP("--METHOD-- initializePredictorCorrector: get past solution mOld=%d, t0-dt0=%9.3e\n",mOld,t0-dt0);
  e.assignGridFunction( gf[mOld].u,t0-dt0 );
  updateStateVariables(gf[mOld]); // *wdh* 080204 
  if( parameters.useConservativeVariables() )
    gf[mOld].primitiveToConservative();

  // For BDF or IMEX-BDF schemes we need more past solutions
  for( int kgf=2; kgf<=numberOfPastTimes; kgf++ )
  {
     #If #METHOD == "BDF"
       // BDF scheme grid index counts backward for past times
       const int mgf = (mCur - kgf + numberOfGridFunctions) % numberOfGridFunctions;
     #Else
       // PC and IMEX-BDF scheme grid index counts forward for past time 
       const int mgf = (mCur + kgf + numberOfGridFunctions) % numberOfGridFunctions;
     #End
     const real tgf = t0-dt0*kgf;
     if( true )
       printF("\n --METHOD-- init past time solution gf[mgf=%i] at t=%9.3e numberOfGridFunctions=%i " 
              "numberOfPastTimes=%i orderOfTimeAccuracy=%i mCur=%i, kgf=%i\n",
              mgf,tgf,numberOfGridFunctions,numberOfPastTimes,orderOfTimeAccuracy,mCur,kgf);
     
     if( movingGridProblem() )
     {
       // **CHECK ME: dt0*kgf ? or -dt0*kgf
       // Note: on input gf[mgf].t=0 indicates the initial grid in gf[mgf] is located at t=0
       moveGrids( t0,t0,tgf,dt0*kgf,gf[mCur],gf[mCur],gf[mgf] );// this will set gf[mgf].t=tgf
     }
     else
     {
       gf[mgf].t=tgf;
     }
     
     gf[mgf].u.updateToMatchGrid(gf[mgf].cg); 
     e.assignGridFunction( gf[mgf].u,tgf );
     updateStateVariables(gf[mgf]); 
     if( parameters.useConservativeVariables() )
       gf[mgf].primitiveToConservative();
     if( false )
     {
       ::display(gf[mgf].u[0],sPrintF("--METHOD-- past time solution gf[mgf=%i].u t=%9.3e",mgf,tgf),"%6.3f ");
     }
     

  }
  
  // For IMEX-BDF schemes we need more past time-derivatives
  if( true && implicitMethod==Parameters::implicitExplicitMultistep  ) // *wdh* Feb. 3, 2017
  {
    for( int kgf=1; kgf<=numberOfPastTimeDerivatives; kgf++ )
    {
      const int mgf = (mCur + kgf + numberOfGridFunctions) % numberOfGridFunctions;
      const int ngf = (nab0 + kgf + numberOfTimeDerivativeLevels) % numberOfTimeDerivativeLevels;
      const real tgf = t0-dt0*kgf;
      gf[mgf].t=tgf;
      
      if( true )
        printF("--METHOD-- init past time du/dt at t=%9.3e (gf[mgf=%i].t=%9.3e) fn[ngf=%i]\n",
           tgf,mgf,gf[mgf].t,ngf);

      // -- evaluate du/dt(t-dt) --
      for( int grid=0; grid<gf[mCur].cg.numberOfComponentGrids(); grid++ )
      {
         
        rparam[0]=gf[mgf].t;
        rparam[1]=gf[mgf].t;     // tforce
        rparam[2]=gf[mgf].t+dt0; // tImplicit  ******** CHECK ME *************
        iparam[0]=grid;
        iparam[1]=gf[mgf].cg.refinementLevelNumber(grid);
        iparam[2]=numberOfStepsTaken;
        getUt(gf[mgf].u[grid],gf[mgf].getGridVelocity(grid),fn[ngf][grid],iparam,rparam,
              utImplicit[grid],&gf[mgf].cg[grid]);

        if( false )
        {
          ::display(gf[mgf].u[grid],sPrintF("--METHOD-- past time u gf[mgf=%i] t=%9.3e",mgf,tgf),"%6.3f ");
          ::display(fn[ngf][grid],sPrintF("--METHOD-- past time du/dt fn[ngf=%i] t=%9.3e",ngf,tgf),"%6.3f ");
        }
      }
      // *wdh* *new* June 7, 2017 **CHECK ME**
      // save past time values of p and ghost u for the 4th order method
      // NOTE: PAST time values are saved in a funny place:
      // save p for use when extrapolating in time
      //    ua(.,.,.,pc)= p(t-2*dt)  (for 2nd/4th order)
      //    ub(.,.,.,pc)= p(t-3*dt)  (for 4th order)
      //    uc(.,.,.,pc)= p(t-4*dt)  (for 4th order)
      // *** savePressureAndGhostVelocity(tgf,ngf);

    }
  }
  else
  {
    // *old* 
    if( numberOfPastTimeDerivatives>0 )
    {
      // -- evaluate du/dt(t-dt) --
      for( int grid=0; grid<gf[mCur].cg.numberOfComponentGrids(); grid++ )
      {
        rparam[0]=gf[mOld].t;
        rparam[1]=gf[mOld].t; // tforce
        rparam[2]=gf[mOld].t+dt0; // tImplicit  **** check me ***
        iparam[0]=grid;
        iparam[1]=gf[mOld].cg.refinementLevelNumber(grid);
        iparam[2]=numberOfStepsTaken;
        getUt(gf[mOld].u[grid],gf[mOld].getGridVelocity(grid),fn[nab1][grid],iparam,rparam,
              utImplicit[grid],&gf[mOld].cg[grid]);
      }
    }
  }
  
  // display(fn[nab1][0],sPrintF("ut(t-dt) from getUt at t=%e\n",gf[mOld].t),debugFile,"%5.2f ");
 
  if( false ) // for testing assign du/dt(t-dt) from TZ directly
  {
    for( int grid=0; grid<gf[mOld].cg.numberOfComponentGrids(); grid++ )
    {
      MappedGrid & c = gf[mOld].cg[grid];
      getIndex(c.dimension(),I1,I2,I3);
      Range Na(0,parameters.dbase.get<int >("numberOfComponents")-1);
      fn[nab1][grid](I1,I2,I3,Na)=e.t(c,I1,I2,I3,Na,t0-dt0); 
 
      if( parameters.gridIsMoving(grid) )
      { // add on gDot.grad(u)
        const realArray & gridVelocity = gf[mOld].getGridVelocity(grid);
        const int na=parameters.dbase.get<int >("uc"), nb=na+c.numberOfDimensions()-1;   // ***** watch out ***
        for( int n=na; n<=nb; n++ )
        {
          fn[nab1][grid](I1,I2,I3,n)+=gridVelocity(I1,I2,I3,0)*e.x(c,I1,I2,I3,n,t0-dt0)+
            gridVelocity(I1,I2,I3,1)*e.y(c,I1,I2,I3,n,t0-dt0);
          if( c.numberOfDimensions()>2 )
            fn[nab1][grid](I1,I2,I3,n)+=gridVelocity(I1,I2,I3,2)*e.z(c,I1,I2,I3,n,t0-dt0);
        }
            
        display(fn[nab1][grid],sPrintF("METHOD:init: ut(t-dt) grid=%i from TZ at t=%e\n",grid,gf[mOld].t),debugFile,"%5.2f ");

      }
    }
  }
       
 
  if( numberOfExtraPressureTimeLevels>0 )
  {
    // get extra time levels for extrapolating the pressure in time (as an initial guess for iterative solvers)
    for( int i=0; i<numberOfExtraPressureTimeLevels; i++ )
    {
      // if orderOfPC==2 : we need p(t-2*dt), p(t-3*dt) ...
      //             ==4 : we need p(t-5*dt), ... **check**
      const real tp = t0-dt0*(i+orderOfPredictorCorrector);
      realCompositeGridFunction & pp = previousPressure[i];
      Range all;
      for( int grid=0; grid<gf[mCur].cg.numberOfComponentGrids(); grid++ )
      {
        e.gd( pp[grid],0,0,0,0,all,all,all,parameters.dbase.get<int >("pc"),tp);
      }
    }
  }
  

  if( debug() & 4 || debug() & 64 )
  {
    for( int grid=0; grid<gf[mCur].cg.numberOfComponentGrids(); grid++ )
    {
      aString buff;
      display(gf[mOld].u[grid],sPrintF(buff,"\n--METHOD-- Init:gf[mOld].u grid=%i : du/dt(t) t=%9.3e",grid,gf[mOld].t),
              debugFile,"%9.3e ");
      if( numberOfPastTimeDerivatives>0 )
        display(fn[nab1][grid],sPrintF(buff,"\n--METHOD-- Init:fn[nab1] grid=%i : du/dt(t) t=%9.3e",grid,gf[mOld].t),
               debugFile,"%9.3e ");
      if( parameters.isMovingGridProblem() )
      {
        display(gf[mOld].getGridVelocity(grid),sPrintF("--METHOD-- t=-dt: gridVelocity[%i] at t=%9.3e\n",grid,gf[mOld].t),debugFile,"%5.2f ");
        display(gf[mCur].getGridVelocity(grid),sPrintF("--METHOD-- t=0 : gridVelocity[%i] at t=%9.3e\n",grid,gf[mCur].t),debugFile,"%5.2f ");
      }
      if( debug() & 64 && parameters.isMovingGridProblem() )
      {
        display(gf[mOld].cg[grid].vertex(),sPrintF("--METHOD-- gf[mOld].cg[%i].vertex at t=%9.3e\n",grid,gf[mOld].t),debugFile,"%7.4f ");
        display(gf[mCur].cg[grid].vertex(),sPrintF("--METHOD-- gf[mCur].cg[%i].vertex at t=%9.3e\n",grid,gf[mCur].t),debugFile,"%7.4f ");
      }
      
    }
    
  }
  if( debug() & 4 )
  {
    if( parameters.isMovingGridProblem() )
    {
      determineErrors( gf[mOld].u,gf[mOld].gridVelocity, gf[mOld].t, 0, error,
                       sPrintF("--METHOD-- errors in u at t=%9.3e (t0-dt0=%9.3e)\n",gf[mOld].t,t0-dt0) );
      #If #METHOD == "BDF"
      #Else
        if( numberOfPastTimeDerivatives>0 )
        {
          fn[nab1].updateToMatchGrid(gf[mOld].cg);  // for moving grid TZ to get errors correct
          determineErrors( fn[nab1],gf[mOld].gridVelocity, gf[mOld].t, 1, error,
                         sPrintF("--METHOD-- errors in ut (fn[nab1]) at t=%9.3e (t0-dt0=%9.3e)\n",gf[mOld].t,t0-dt0) );
        }
      #End
    }
    
  }
 
       
}
else  
{

  // **************************************************************************
  // ************************ REAL RUN ****************************************
  // ****************** Initialize for NOT twilightZoneFlow *******************
  // **************************************************************************

  // printF(" **************** METHOD: still need correct initial values for du/dt(t-dt)  ****** \n");
  // printF(" **************** use values from du/dt(t)                                  ****** \n");
  printF("\n--METHOD-- Initialize past time values for scheme, numberOfPastTimes=%i"
         " numberOfPastTimeDerivatives=%i ---\n",numberOfPastTimes,numberOfPastTimeDerivatives);
  if( debug() & 2 )
    fPrintF(debugFile,"--METHOD-- Initialize past time values for scheme, numberOfPastTimes=%i"
            " numberOfPastTimeDerivatives=%i ---\n",numberOfPastTimes,numberOfPastTimeDerivatives);

  if( parameters.useConservativeVariables() )
    gf[mCur].primitiveToConservative();
 
  // if( parameters.isAdaptiveGridProblem() )
  //   gf[mOld].u.updateToMatchGrid(gf[mOld].cg);  // 040928 -- why is this needed here ?
  
  if( debug() & 8 )
  {
    printF(" PC: init: gf[mOld].u.numberOfGrids=%i \n",gf[mOld].u.numberOfGrids());
    printF(" PC: init: gf[mOld].cg.numberOfComponentGrids=%i \n",gf[mOld].cg.numberOfComponentGrids());
  
    printF(" PC: init: gf[mCur].u.numberOfGrids=%i \n",gf[mCur].u.numberOfGrids());
    printF(" PC: init: gf[mCur].cg.numberOfComponentGrids=%i \n",gf[mCur].cg.numberOfComponentGrids());
  }
  
  if( !parameters.dbase.get<bool>("useNewTimeSteppingStartup") )
  {
    printF(" -- METHOD-- USE OLD STARTUP, Set past time solutions to t=0 solution \n");
    if( numberOfPastTimes==1 )
    {
      // uOld=uCur 
      assign(gf[mOld].u,gf[mCur].u); 
      gf[mOld].form=gf[mCur].form;
    }
    else
    { // June 8, 2017 *wdh*
      for( int kgf=1; kgf<=numberOfPastTimes; kgf++ )
      {
        const int mgf = (mCur + kgf + numberOfGridFunctions) % numberOfGridFunctions;
        assign(gf[mgf].u,gf[mCur].u); 
        gf[mgf].t=t0-dt0*kgf;
        gf[mgf].form=gf[mCur].form;
      }
    }
    
  }
  else
  {
    // *new* way to initialize past time solution  // *wdh* 2014/06/28 
    printF(" -- METHOD-- USE NEW STARTUP numberOfPastTimes=%i mCur=%i mOld=%i\n",numberOfPastTimes,mCur,mOld);
    if( numberOfPastTimes==1 )
    {
      gf[mOld].t=t0-dt0;
      int previous[1]={mOld};  // 
      getPastTimeSolutions( mCur, numberOfPastTimes, previous  ); 

    }
    else
    {
      // For BDF schemes we need more past solutions (NOTE: this does not work for PC since previous[0]!=mOld)
      int *previous = new int[numberOfPastTimes];
      for( int kgf=1; kgf<=numberOfPastTimes; kgf++ )
      {
        #If #METHOD == "BDF"
          // BDF scheme grid index counts backward for past times        
          const int mgf = (mCur - kgf + numberOfGridFunctions) % numberOfGridFunctions; // *wdh* April 21, 2021
        #Else
          const int mgf = (mCur + kgf + numberOfGridFunctions) % numberOfGridFunctions;
        #End
        gf[mgf].t=t0-dt0*kgf;
        gf[mgf].form=gf[mCur].form;
        previous[kgf-1]=mgf;
      }
      getPastTimeSolutions( mCur, numberOfPastTimes, previous  );
      delete [] previous;
    }
    
  }
  
  // gf[mOld].form=gf[mCur].form;
 
  // For PC and IMEX-BDF schemes we need more past time-derivatives
  if( parameters.dbase.get<Parameters::TimeSteppingMethod >("timeSteppingMethod")==Parameters::adamsPredictorCorrector4 || 
      implicitMethod==Parameters::implicitExplicitMultistep )
  {
    // ** THIS SECTION IS REPEATED FROM ABOVE -- *FIX ME*

    for( int kgf=1; kgf<=numberOfPastTimeDerivatives; kgf++ )
    {
      const int mgf = (mCur + kgf + numberOfGridFunctions) % numberOfGridFunctions;
      const int ngf = (nab0 + kgf + numberOfTimeDerivativeLevels) % numberOfTimeDerivativeLevels;
      const real tgf = t0-dt0*kgf;
      gf[mgf].t=tgf;
      
      if( true )
        printF("--METHOD-- init past time du/dt at t=%9.3e (gf[mgf=%i].t=%9.3e) fn[ngf=%i]\n",
           tgf,mgf,gf[mgf].t,ngf);

      // -- evaluate du/dt(t-dt) --
      for( int grid=0; grid<gf[mCur].cg.numberOfComponentGrids(); grid++ )
      {
         
        rparam[0]=gf[mgf].t;
        rparam[1]=gf[mgf].t; // tforce
        rparam[2]=gf[mCur].t-gf[mgf].t; // tImplicit  *************** check me 090806 **********************
        iparam[0]=grid;
        iparam[1]=gf[mgf].cg.refinementLevelNumber(grid);
        iparam[2]=numberOfStepsTaken;
        getUt(gf[mgf].u[grid],gf[mgf].getGridVelocity(grid),fn[ngf][grid],iparam,rparam,
              utImplicit[grid],&gf[mgf].cg[grid]);

        if( false )
        {
          ::display(fn[ngf][grid],sPrintF("--METHOD-- past time du/dt fn[ngf=%i] t=%9.3e",ngf,tgf),"%6.3f ");
        }
      }
    }



    // *wdh* *new* June 7, 2017 **CHECK ME**
    if( orderOfPredictorCorrector==4 )
    {
      for( int m=0; m<=1; m++ )
      {
        // save past time values of p and ghost u for the 4th order method
        // NOTE: PAST time values are saved in a funny place:
        // save p for use when extrapolating in time
        //    ua(.,.,.,pc)= p(t-2*dt)  (for 2nd/4th order)
        //    ub(.,.,.,pc)= p(t-3*dt)  (for 4th order)
        //    uc(.,.,.,pc)= p(t-4*dt)  (for 4th order)
        assert( nab0==0 );
        const int nabPastTime=(nab0+m);
        real tp=t0-(m+2)*dt0;     
        savePressureAndGhostVelocity(tp,nabPastTime);
      }
    }

  }
  else
  {
    // *********** OLD WAY *******
            
    if( numberOfPastTimeDerivatives>0 )
    {
      if( debug() & 2 )
        fPrintF(debugFile,"--METHOD-- get past time du/dt at t=%9.3e...\n",gf[mOld].t);

      for( int grid=0; grid<gf[mCur].cg.numberOfComponentGrids(); grid++ )
      {
        rparam[0]=gf[mOld].t;
        rparam[1]=gf[mOld].t; // tforce
        // *wdh* 090806 : what was this? rparam[2]=gf[mCur].t-gf[mOld].t; // tImplicit
        rparam[2]=gf[mCur].t; // tImplicit = apply forcing for implicit time stepping at this time
        iparam[0]=grid;
        iparam[1]=gf[mOld].cg.refinementLevelNumber(grid);
        iparam[2]=numberOfStepsTaken;
        getUt(gf[mOld].u[grid],gf[mOld].getGridVelocity(grid),fn[nab1][grid],iparam,rparam,
              utImplicit[grid],&gf[mOld].cg[grid]);
      }
    }
  
    if( debug() & 4 )
    {
      #If #METHOD == "BDF"
      #Else      
      determineErrors( fn[nab1],gf[mOld].gridVelocity, gf[mOld].t, 1, error,
                       sPrintF(" PC:init: du/dt at past time t=%e \n",gf[mOld].t) );
      #End
    }
  
    for( int grid=0; grid<gf[mOld].cg.numberOfComponentGrids(); grid++ )
    {
      MappedGrid & c = gf[mOld].cg[grid];
      getIndex(c.dimension(),I1,I2,I3);
      // fn[nab1][grid](I1,I2,I3,N)=fn[nab0][grid](I1,I2,I3,N);
      if( orderOfPredictorCorrector==4 )
      {
        printF("INFO -- PC order 4 : Setting 2 past time derivatives to values at t=0. **FIX ME**\n");
        
        for( int m=0; m<=1; m++ )
        {
          const int nab=(mOld+m+1) % 4;

          // *** WE COULD DO BETTER HERE ***
          assign(fn[nab][grid],fn[nab0][grid],I1,I2,I3,N);

        }
      }
    }
    // *wdh* *new* June 7, 2017 **CHECK ME**
    if( orderOfPredictorCorrector==4 )
    {
      for( int m=0; m<=1; m++ )
      {
        // save past time values of p and ghost u for the 4th order method
        // NOTE: PAST time values are saved in a funny place:
        // save p for use when extrapolating in time
        //    ua(.,.,.,pc)= p(t-2*dt)  (for 2nd/4th order)
        //    ub(.,.,.,pc)= p(t-3*dt)  (for 4th order)
        //    uc(.,.,.,pc)= p(t-4*dt)  (for 4th order)
        assert( nab0==0 );
        const int nabPastTime=(nab0+m);
        real tp=t0-(m+2)*dt0;     
        savePressureAndGhostVelocity(tp,nabPastTime);
      }
    }
    

  } // end OLD WAY
  
 
}
     
 
dtb=dt0;    // delta t to go from ub to ua
 
dtp[0]=dt0;
dtp[1]=dt0;
dtp[2]=dt0;
dtp[3]=dt0;
dtp[4]=dt0;
 
//       if( debug() & 8 )
//       {
//         fprintf(debugFile," advance Adams PC: ut at t0=%e\n",t0);
//         outputSolution( fn[nab0],t0 );
//         fprintf(debugFile," advance Adams PC: Errors in ut at t0=%e\n",t0);
//         determineErrors( fn[1],t0-dt0 );
//       }
 
init=false;


#endMacro


// =======================================================================================================
//    Macro to move the grids at the start of a PC time step.
// Arguments:
//    METHOD : name of the calling function (for debug output)
//    utImplicit : name of the grid function that holds the explicit part of the implicit operator.
//    
//    predictorOrder : order of the predictor corrector
//    ub,uc,ud : grid functions that hold du/dt at times tb, tc, td
//               If predictorOrder==2 then explosed points are filled in on ub.
//               If predictorOrder==3 then explosed points are filled in on ub and uc.
//               If predictorOrder==4 then explosed points are filled in on ub, uc and ud.
// =======================================================================================================
#beginMacro moveTheGridsMacro(METHOD,utImplicit,predictorOrder,tb,ub,tc,uc,td,ud)

if( movingGridProblem() )
{
  checkArrays(" METHOD : before move grids"); 

  if( debug() & 8 )
    printF(" METHOD: before moveTheGridsMacro: t0=%9.3e, gf[mNew].t=%9.3e, gf[mNew].gridVelocityTime=%9.3e\n",
           t0,gf[mNew].t,gf[mNew].gridVelocityTime);
  if( debug() & 4 )
    fPrintF(debugFile," METHOD: before moveTheGridsMacro: t0=%9.3e, gf[mNew].t=%9.3e, gf[mNew].gridVelocityTime=%9.3e\n",
           t0,gf[mNew].t,gf[mNew].gridVelocityTime);

  // generate gf[mNew] from gf[mCur] (compute grid velocity on gf[mCur] and gf[mNew]
  moveGrids( t0,t0,t0+dt0,dt0,gf[mCur],gf[mCur],gf[mNew] ); 

  checkArrayIDs(sPrintF(" METHOD : after move grids t=%9.3e",gf[mCur].t));

  if( parameters.isAdaptiveGridProblem() )
  {
    // both moving and AMR 
    parameters.dbase.get<Ogen* >("gridGenerator")->updateRefinement(gf[mNew].cg);
  }
      
  if( debug() & 16 )
  {
    if( twilightZoneFlow() )
    {
      fprintf(debugFile,"\n ---> METHOD : Errors in u after moveGrids t=%e  \n",gf[mCur].t);
      determineErrors( gf[mCur] );
    }
  }


  real cpu0=getCPU();
  gf[mNew].cg.rcData->interpolant->updateToMatchGrid( gf[mNew].cg );  
  parameters.dbase.get<RealArray>("timing")(parameters.dbase.get<int>("timeForUpdateInterpolant"))+=getCPU()-cpu0;

  cpu0=getCPU();
  gf[mNew].u.getOperators()->updateToMatchGrid(gf[mNew].cg); 
  parameters.dbase.get<RealArray>("timing")(parameters.dbase.get<int>("timeForUpdateOperators"))+=getCPU()-cpu0;

  if( debug() & 4 ) printf("METHOD : step: update gf[mNew] for moving grids, gf[mNew].t=%9.3e,...\n",gf[mNew].t);

  if( debug() & 16 )
  {
    if( twilightZoneFlow() )
    {
      fprintf(debugFile,"\n ---> METHOD: Errors in u before updateForMovingGrids t=%e  \n",gf[mCur].t);
      fprintf(debugFile,"*** mCur=%i mNew=%i numberOfGridFunctions=%i *** \n",
              mCur,mNew,numberOfGridFunctions);
          
      determineErrors( gf[mCur] );
    }
  }

  updateForMovingGrids(gf[mNew]);
  // ****      gf[mNew].u.updateToMatchGrid( gf[mNew].cg );  

  checkArrayIDs(sPrintF(" METHOD : after updateForMovingGrids t=%9.3e",gf[mCur].t));

  if( debug() & 16 )
  {
    if( twilightZoneFlow() )
    {
      fprintf(debugFile,"\n ---> METHOD: Errors in u after updateForMovingGrids t=%e  \n",gf[mCur].t);
      determineErrors( gf[mCur] );
    }
  }

  // get values for exposed points on gf[mCur]
  if( parameters.useConservativeVariables() )
    gf[mCur].primitiveToConservative();  // *wdh* 010318

  bool useNewIMEX= parameters.dbase.get<Parameters::TimeSteppingMethod >("timeSteppingMethod")==Parameters::adamsPredictorCorrector4 
                     || implicitMethod==Parameters::implicitExplicitMultistep;

  cpu0=getCPU();
  if( useNewExposedPoints && parameters.dbase.get<int>("simulateGridMotion")==0 )
  {
    if( Parameters::checkForFloatingPointErrors )
      checkSolution(gf[mCur].u,"Before interp exposed points",true);

    if( debug() & 16 )
    {
      if( twilightZoneFlow() )
      {
        fprintf(debugFile,"\n ---> METHOD: Errors in u BEFORE interp exposed t=%e  \n",gf[mCur].t);
        determineErrors( gf[mCur] );
      }
    }

    // Save an array of the number of exposed points on each grid and each time level so that
    // we can avoid recomputing du/dt where it is not needed *wdh* Jan 7, 2019
    IntegerArray numExposed;
    if( numberOfPastTimeDerivatives>0 )
    {
      if( useNewIMEX && numberOfPastTimes < numberOfPastTimeDerivatives )
      {
        printF("ims:ERROR: numberOfPastTimes=%i must be at least numberOfPastTimeDerivatives=%i for moving grids\n"
               "           so that exposed points on past du/dt can be computed\n",numberOfPastTimes,numberOfPastTimeDerivatives);
        OV_ABORT("Error");
      }
      
      
      numExposed.redim(gf[mCur].cg.numberOfComponentGrids(),max(numberOfPastTimes,numberOfPastTimeDerivatives)+1);
    }
    
     
    ExposedPoints exposedPoints;
    exposedPoints.setAssumeInterpolationNeighboursAreAssigned(parameters.dbase.get<int >("extrapolateInterpolationNeighbours"));
    exposedPoints.initialize(gf[mCur].cg,gf[mNew].cg,parameters.dbase.get<int >("stencilWidthForExposedPoints"));
    exposedPoints.interpolate(gf[mCur].u,(twilightZoneFlow() ? parameters.dbase.get<OGFunction* >("exactSolution") : NULL),t0);

    if( numberOfPastTimeDerivatives>0 )
    { // save the number of exposed points
      for( int grid=0; grid<gf[mCur].cg.numberOfComponentGrids(); grid++ )
        numExposed(grid,0)= exposedPoints.getNumberOfExposedPoints(grid);
    }
    
    if( true )
    {
      printF("gf[%i]: t=%9.3e, exposed points=[",mCur,gf[mCur].t);
      for( int grid=0; grid<gf[mCur].cg.numberOfComponentGrids(); grid++ )
        printF("%i,",exposedPoints.getNumberOfExposedPoints(grid));
      printF("]\n");
    }
    
    // exposedPoints.getNumberOfExposedPoints(grid)
    
    // Added for BDF: *wdh* 2015/04/05
    for( int kp=1; kp<=numberOfPastTimes; kp++ )
    {
      #ifdef USE_NEW_OLD_TIME_INDEX
        const int mPast = gfList[kp];
        const real tgf = tList[kp];
      #else
        // **** CHECK ME SIGN of kp ****
        // BDF scheme may use -kp
        // const int mPast = (mCur -kp + numberOfGridFunctions) % numberOfGridFunctions; 
        const int mPast = (mCur + kp + numberOfGridFunctions) % numberOfGridFunctions; // *wdh* Jan9, 2019
        const real tgf = t0-dt0*kp;
      #endif

      printF("exposedPoints: mCur=%i, mPast=%i, t0=%9.3e, t=%9.3e\n",mCur,mPast,t0,gf[mPast].t);
      exposedPoints.initialize(gf[mPast].cg,gf[mNew].cg,parameters.dbase.get<int >("stencilWidthForExposedPoints"));
      exposedPoints.interpolate(gf[mPast].u,(twilightZoneFlow() ? parameters.dbase.get<OGFunction* >("exactSolution") : NULL),t0);

      if( numberOfPastTimeDerivatives>0 )
      { // save the number of exposed points
        assert( kp<=numExposed.getBound(1) );
        
        for( int grid=0; grid<gf[mCur].cg.numberOfComponentGrids(); grid++ )
          numExposed(grid,kp)= exposedPoints.getNumberOfExposedPoints(grid);
      }

      if( true )
      {
        printF("gf[%i]: t=%9.3e, exposed points=[",mPast,gf[mPast].t);
        for( int grid=0; grid<gf[mPast].cg.numberOfComponentGrids(); grid++ )
          printF("%i,",exposedPoints.getNumberOfExposedPoints(grid));
        printF("]\n");
      }
    }

    // For IMEX-BDF schemes we need more past time-derivatives
    if( useNewIMEX ) 
    {
      // for( int kgf=1; kgf<numberOfPastTimeDerivatives; kgf++ )   //  DO NOT FIXUP current time *wdh* Jan 7, 2019
      for( int kgf=1; kgf<=numberOfPastTimeDerivatives; kgf++ )   //  DO NOT FIXUP current time *wdh* Jan 7, 2019
      {
        const int & numberOfTimeDerivativeLevels = parameters.dbase.get<int>("numberOfTimeDerivativeLevels"); // *CHECK ME**
        #ifdef USE_NEW_OLD_TIME_INDEX
          const int mgf = gfList[kgf];
          const int ngf = fnList[kgf];
          const real tgf = tList[kgf];
        #else
          // **** CHECK ME SIGNS ****
          const int mgf = (mCur + kgf + numberOfGridFunctions) % numberOfGridFunctions;
          const int ngf = (nab0 + kgf + numberOfTimeDerivativeLevels) % numberOfTimeDerivativeLevels;
          const real tgf = t0-dt0*kgf;
        #endif
        if( fabs(gf[mgf].t -tgf) > dt0*1.e-5  )
        {
          printF("moveGrids: ERROR: gf[mgf].t =%10.3e IS NOT EQUAL TO tgf=%10.3e, t0=%10.3e, mCur=%i, mgf=%i, kgf=%i\n",
                 gf[mgf].t,tgf,t0,mCur,mgf,kgf);
          printF("numberOfGridFunctions=%i, numberOfPastTimeDerivatives=%i\n",numberOfGridFunctions,numberOfPastTimeDerivatives);
          #ifdef USE_NEW_OLD_TIME_INDEX
            printF("gfList=%i,%i,%i\n",gfList[0],gfList[1],gfList[2]);
            printF("fnList=%i,%i,%i\n",fnList[0],fnList[1],fnList[2]);
            printF("tList=%g,%g,%g\n",tList[0],tList[1],tList[2]);
          
          #endif
          printF("gf[0].t=%9.3e., gf[1].t=%9.3e, gf[2].t=%9.3e\n",gf[0].t,gf[1].t,gf[2].t);
          
          OV_ABORT("error");
        }
        
        if( true )
        {
          printF("--METHOD-- moveGrids: numberOfGridFunctions=%i, numberOfTimeDerivativeLevels=%i"
                 " numberOfPastTimes=%i, numberOfPastTimeDerivatives=%i\n",
                 numberOfGridFunctions,numberOfTimeDerivativeLevels,numberOfPastTimes,numberOfPastTimeDerivatives);
        }
        
        if( true )
          printF("--METHOD-- Evaluate previous time du/dt at t=%9.3e (gf[mgf=%i].t=%9.3e) fn[ngf=%i]\n",
                 tgf,mgf,gf[mgf].t,ngf);
  
        // -- evaluate du/dt(t-dt) --
        for( int grid=0; grid<gf[mCur].cg.numberOfComponentGrids(); grid++ )
        {
          assert( kgf <= numExposed.getBound(1) );
          if( numExposed(grid,kgf)==0 )
          {
            if( true )
              printF(" ---> Skip computing du/dt on past time level kgf=%i grid=%i since numExposed=%i\n",
                     kgf,grid,numExposed(grid,kgf));
            continue;
          }
          else
          {
            if( true )
              printF(" ---> Compute du/dt on past time level kgf=%i grid=%i since numExposed=%i\n",
                     kgf,grid,numExposed(grid,kgf));
        
          }
          
           
          rparam[0]=gf[mgf].t;
          rparam[1]=gf[mgf].t; // tforce
          rparam[2]=gf[mgf].t+dt0; // tImplicit:  **CHECK ME Jan 7, 2019
          iparam[0]=grid;
          iparam[1]=gf[mgf].cg.refinementLevelNumber(grid);
          iparam[2]=numberOfStepsTaken-kgf;
          // RealArray f0;
          // f0=fn[ngf][grid];
          
          getUt(gf[mgf].u[grid],gf[mgf].getGridVelocity(grid),fn[ngf][grid],iparam,rparam,utImplicit[grid],
                &gf[mNew].cg[grid]); // Use mask from new time 
  
          /* ----
          real maxDiff = max(fabs(f0-fn[ngf][grid]));
          if( true && maxDiff > 1.e-10 )
          {
            printF("ERROR: computing du/dt for  grid=%i, answer has changed, maxDiff=%9.3e\n",grid,maxDiff);
            if( false )
            {
              ::display(gf[mgf].u[grid],sPrintF("--METHOD-- u, mgf=%i, at t =%9.3e",mgf,tgf),"%6.3f ");
              ::display(f0,sPrintF("--METHOD-- INITIAL past time du/dt fn[ngf=%i] t=%9.3e",ngf,tgf),"%6.3f ");
              // ::display(fn[ngf][grid],sPrintF("--METHOD-- past time du/dt fn[ngf=%i] t=%9.3e",ngf,tgf),"%6.3f ");
              ::display(fn[ngf][grid]-f0,sPrintF("--METHOD-- DIFFERENCE t=%9.3e",tgf),"%9.2e ");
              // OV_ABORT("stop here for now");
            }
          }
          --- */

          if( false )
          {
            ::display(fn[ngf][grid],sPrintF("--METHOD-- past time du/dt fn[ngf=%i] t=%9.3e",ngf,tgf),"%6.3f ");
          }
        }
      }
    }
    

    if( debug() & 16 )
    {
      if( twilightZoneFlow() )
      {
        fprintf(debugFile,"\n ---> METHOD: Errors in u AFTER interp exposed t=%e  \n",gf[mCur].t);
        determineErrors( gf[mCur] );
      }
    }

    if( predictorOrder>=2 && !useNewIMEX  )
    {
      // -------------------------
      // --- fixup du/dt(t-dt) ---
      // -------------------------

      // NOTE: we CANNOT directly interpolate points on du/dt since for moving grids
      // du/dt includes the -gDot.grad(u) term which differs from grid to grid 

      // Current procedure: 
      //   1. Interpolate exposed points on u(t-dt)
      //   2. Recompute du/dt(t-dt) 

      // Optimizations: 
      //   - only recompute du/dt(t-dt) on grids with exposed points
      //   - could only compute du/dt(t-dt) on those points where is not already known.

      if( gf[mCur].t<=0. || debug() & 4 )
      {
        printF(" --- INFO: Fixup exposed points of u(t-dt) and recompute du/dt(t-dt) t=%9.3e, tb=%9.3e ----- \n"
               "     The extra work involved in recomputing du/dt(t-dt) can be avoided by using "
               "the option 'first order predictor'.\n",gf[mCur].t,tb);
      }
      
      ExposedPoints exposedPoints;

      //            
      // exposedPoints.setExposedPointType(ExposedPoints::exposedDiscretization);

      const int extrapolateInterpolationNeighbours=parameters.dbase.get<int >("extrapolateInterpolationNeighbours");
      const int stencilWidthForExposedPoints=parameters.dbase.get<int >("stencilWidthForExposedPoints");
      
      if( debug() & 4 )
      {
        fPrintF(debugFile," ---- compute exposed for du/dt(t-dt), extrapolateInterpolationNeighbours=%i, "
                "stencilWidthForExposedPoints=%i\n",extrapolateInterpolationNeighbours,stencilWidthForExposedPoints);
      }
      

      exposedPoints.setAssumeInterpolationNeighboursAreAssigned(extrapolateInterpolationNeighbours);

      exposedPoints.initialize(gf[mOld].cg,gf[mNew].cg,stencilWidthForExposedPoints);
      exposedPoints.interpolate(gf[mOld].u,(twilightZoneFlow() ? 
                                parameters.dbase.get<OGFunction* >("exactSolution") : NULL),gf[mOld].t);
          
      // For now recompute du/dt(t-dt) using the mask values from cg(t+dt)
      for( int grid=0; grid<gf[mOld].cg.numberOfComponentGrids(); grid++ )
      {
             
        if( gridWasAdapted || exposedPoints.getNumberOfExposedPoints(grid)>0 )
        {
          if( debug() & 2 )
          {
            printF(" ---- METHOD: recompute du/dt(t-dt) for grid=%i t-dt = %9.3e  (%i exposed)-----\n",grid,gf[mOld].t,
                   exposedPoints.getNumberOfExposedPoints(grid));
            fPrintF(debugFile," ---- METHOD: recompute du/dt(t-dt) for grid=%i t-dt = %9.3e  (%i exposed)-----\n",
                   grid,gf[mOld].t,exposedPoints.getNumberOfExposedPoints(grid));
          }
          
          // This is only necesssary if there are exposed points on this grid
          rparam[0]=gf[mOld].t;
          rparam[1]=gf[mOld].t;
          rparam[2]=gf[mCur].t; // tImplicit
          iparam[0]=grid;
          iparam[1]=gf[mOld].cg.refinementLevelNumber(grid);
          iparam[2]=numberOfStepsTaken;

          getUt(gf[mOld].u[grid],gf[mOld].getGridVelocity(grid),ub[grid],iparam,rparam,
                utImplicit[grid],&gf[mNew].cg[grid]);
              
        }
        else
        {
          if( debug() & 2 )
          {
            printF(" ---- METHOD: fixp du/dt(t-dt) for grid=%i t-dt = %9.3e  (%i exposed) ...ok -----\n",
                   grid,gf[mOld].t,exposedPoints.getNumberOfExposedPoints(grid));
            fPrintF(debugFile," ---- METHOD: fixp du/dt(t-dt) for grid=%i t-dt = %9.3e  (%i exposed) ...ok -----\n",
                   grid,gf[mOld].t,exposedPoints.getNumberOfExposedPoints(grid));
          }
        }
        
      }
      if( debug() & 4 )
      { 
        if( twilightZoneFlow() )
        {
          fprintf(debugFile," ***METHOD: gf[mOld] after interp exposed, gf[mOld].t=%e",gf[mOld].t);
          for( int grid=0; grid<gf[mOld].cg.numberOfComponentGrids(); grid++ )
          {
            display(gf[mOld].u[grid],sPrintF("\n ****gf[mOld].u[grid=%i]",grid),debugFile,"%7.1e ");
          }
          determineErrors( gf[mOld] );

          fprintf(debugFile," ***METHOD: du/dt(t-dt)  after interp exposed, gf[mOld].t=%e",gf[mOld].t);
          for( int grid=0; grid<gf[mOld].cg.numberOfComponentGrids(); grid++ )
          {
            display(ub[grid],sPrintF("\n ****ub[grid=%i]: du/dt(t-dt)",grid),debugFile,"%7.1e ");
          }
          determineErrors( ub,gf[mOld].gridVelocity, gf[mOld].t, 1, error );
        }
      }
      
            
      if( predictorOrder>=3 )
      {
        printF("\n XXXXX METHOD: moveTheGridsMacro:Error: finish me for predictorOrder>=3 XXXXX \n\n");
        // OV_ABORT("METHOD: moveTheGridsMacro:Error: finish me for predictorOrder>=3");
      }
      
    } // end if predictorOrder>=2
        
    if( Parameters::checkForFloatingPointErrors )
      checkSolution(gf[mCur].u,"METHOD: After interp exposed points",true);

  }
  else if( parameters.dbase.get<int>("simulateGridMotion")==0  )
  {
    // *old way* 
    interpolateExposedPoints(gf[mCur].cg,gf[mNew].cg,gf[mCur].u, 
                             (twilightZoneFlow() ? parameters.dbase.get<OGFunction* >("exactSolution") : NULL),t0,
                             false,Overture::nullIntArray(),Overture::nullIntegerDistributedArray(),
                             parameters.dbase.get<int >("stencilWidthForExposedPoints") ); 
  }

    
  if( twilightZoneFlow() && false ) // **** wdh **** 
  {
    // for testing ***
    OGFunction & e = *(parameters.dbase.get<OGFunction* >("exactSolution"));
        
    int grid=0; 
    MappedGrid & c = gf[mNew].cg[grid];
    getIndex(c.dimension(),I1,I2,I3);

    ub[grid](I1,I2,I3,N)=e.t(c,I1,I2,I3,N,t0-dt0); 
        
  }



  parameters.dbase.get<RealArray>("timing")(parameters.dbase.get<int>("timeForInterpolateExposedPoints"))+=getCPU()-cpu0;
  // compute dudt now -- after exposed points have been computed!

  checkArrayIDs(sPrintF(" METHOD : after moving grids update t=%9.3e",gf[mCur].t));


  if( debug() & 16 )
  {
    if( twilightZoneFlow() )
    {
      fprintf(debugFile,"\n ---> METHOD: Errors in u after move grids t=%e  \n",gf[mCur].t);
      determineErrors( gf[mCur] );
    }
  }

  if( debug() & 16 )
    printf(" METHOD: AFTER moveTheGridsMacro: t0=%9.3e, gf[mNew].t=%9.3e, gf[mNew].gridVelocityTime=%9.3e\n",
           t0,gf[mNew].t,gf[mNew].gridVelocityTime);
}


#endMacro



// =======================================================================================================
//    Macro to correct for moving grids.
// Arguments:
//    METHOD : name of the calling function (for debug output)
// Retrun:
//   movingGridCorrectionsHaveConverged = true if this is a moving grid problem and the sub-iteration
//             corrections have converged.
// =======================================================================================================
#beginMacro correctForMovingGridsMacro(METHOD)

  // Correct for forces on moving bodies if we have more corrections.
  bool movingGridCorrectionsHaveConverged = false;
  real delta =0.; // holds relative correction when we are sub-cycling 
  const bool useMovingGridSubIterations= parameters.dbase.get<bool>("useMovingGridSubIterations");
  const int & multiDomainProblem = parameters.dbase.get<int>("multiDomainProblem");

  // *wdh* 2015/12/16 -- explicitly check for useMovingGridSubIterations, otherwise we can do multiple
  //                     corrections always if requested,

  if( false && movingGridProblem() )
    printF("METHOD: moving grid correction step : numberOfCorrections=%i useMovingGridSubIterations=%i\n",
       numberOfCorrections,(int)useMovingGridSubIterations);


  if( movingGridProblem() && ( (numberOfCorrections==1 && !multiDomainProblem ) // *wdh* Sept 15, 2018 added !multiDomainProblem
                              || !useMovingGridSubIterations)  ) // *wdh* 2015/12/16 
  {
    if( numberOfCorrections>10 )
    {
      printF("WARNING: movingGrid problem, useMovingGridSubIterations=false but numberOfCorrections>10\n");
      OV_ABORT("ERROR: this is an error for now");
    }
    
    correctMovingGrids( t0,t0+dt0,gf[mCur],gf[mNew] ); 

  }
// else if( movingGridProblem() && (correction+1)<numberOfCorrections)
  else if( movingGridProblem() )
  {
    // --- we may be iterating on the moving body motion (e.g.for light bodies) ---
    //     After correcting for the motion, check for convergence
    correctMovingGrids( t0,t0+dt0,gf[mCur],gf[mNew] ); 

    // Check if the correction step has converged
    bool isConverged = getMovingGridCorrectionHasConverged();
    delta = getMovingGridMaximumRelativeCorrection();
    if( true || debug() & 2 )
      printF("METHOD: moving grid correction step : delta =%8.2e (correction=%i, isConverged=%i)\n",
             delta,correction+1,(int)isConverged);

    if( isConverged && (correction+1) >=minimumNumberOfPCcorrections )  // note correction+1 
    {
      movingGridCorrectionsHaveConverged=true;  // we have converged -- we can break from correction steps
      if( delta!=0. && debug() & 1 )
        printF("METHOD: moving grid correction step : sub-iterations converged after %i corrections, rel-err =%8.2e\n",
               correction+1,delta);
      // break;  // we have converged -- break from correction steps
    }
    if( !isConverged && (correction+1)>=numberOfCorrections )
    {
      printF("METHOD:ERROR: moving grid corrections have not converged! numberOfCorrections=%i, rel-err =%8.2e\n",
             correction+1,delta);
    }
    

  }
  else 
  {
  }

#endMacro
