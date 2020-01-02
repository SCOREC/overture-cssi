// This file automatically generated from advanceBAMX.bC with bpp.
//
//           BIANISTROPIC MAXWELL'S EQUATIONS
//
#include "Maxwell.h"
#include "display.h"
#include "CompositeGridOperators.h"
#include "ParallelUtility.h"
#include "ParallelGridUtility.h"
#include "DispersiveMaterialParameters.h"

#include "updateOpt.h"

#define advBA EXTERN_C_NAME(advba)
#define mxFilter EXTERN_C_NAME(mxfilter)
#define getGDMParameters EXTERN_C_NAME(getgdmparameters)

extern "C"
{
            void advBA(const int&nd,
            const int&n1a,const int&n1b,const int&n2a,const int&n2b,const int&n3a,const int&n3b,
            const int&nd1a,const int&nd1b,const int&nd2a,const int&nd2b,const int&nd3a,const int&nd3b,
            const int&nd4a,const int&nd4b,
            const int&mask,const real&rx,  
            const real&um, const real&u, real&un, const real&f, const real&fa, real & K0i, int & matMask, 
            const real& pm,const real&p,const real& pt, const real&xy,
            const real& etax,const real & etay, const real & etaz,
            const int&bc, const real &dis, const real &varDis, const int&ipar, const real&rpar, int&ierr );

  void mxFilter(const int&nd,
            const int&nd1a,const int&nd1b,const int&nd2a,const int&nd2b,const int&nd3a,const int&nd3b,
            const int & gridIndexRange, const real & u, const real & d, 
            const int&mask, const int&boundaryCondition, const int&ipar, const real&rpar, int&ierr );
}

static Maxwell *cgmxPointer=NULL; // for getBAGDMParameters

#define getBAGDMParameters EXTERN_C_NAME(getbagdmparameters)
extern "C"
{

// ================================================================================
///  \brief Return the BA gdm parameters  
/// \param grid (input) : return parameters for this grid 
/// 
/// \param gdmPar(0:3,0:NpMax-1,0:5,0:5,0:maxRegions-1) (input/output) : must be allocated on input. Values are returned here.
///
/// \param Npt(1:6,1:6,1:maxRegions) (input/output) : number of polarization terms for K(i,j) i=1,..,6, j=1,...,6
// ================================================================================
void getBAGDMParameters( int & grid, real *gdmPar, int *Npt, 
                   			 const int & NpMax,
                   			 const int & maxRegions )
{
#define Npta(k1,k2,mr) Npt[(k1) + 6*((k2) + 6*( mr ) ) ]
#define gdmVar(m,n,k1,k2,mr) gdmPar[(m) + 4*( (n) + NpMax*( (k1) + 6*( (k2) + 6*( (mr) ) ) ) ) ]

    assert( cgmxPointer != NULL );
    assert( grid==0 );
    
//   CompositeGrid *&cgp = cgmxPointer->cgp;
//   assert( cgp!=NULL );
//   CompositeGrid & cg= *cgp;
// 
//   const int domain = cg.domainNumber(grid);
//   const DispersiveMaterialParameters & dmp =  cgmxPointer->getDomainDispersiveMaterialParameters(domain);
    
    const int numberOfMaterialRegions = cgmxPointer->numberOfMaterialRegions;
  // printF("getBAGDMParameters: number of material regions = %i.\n",numberOfMaterialRegions);
    assert( numberOfMaterialRegions<maxRegions );
    
    std::vector<DispersiveMaterialParameters> & dmpVector = 
        cgmxPointer->dbase.get<std::vector<DispersiveMaterialParameters> >("materialRegionParameters");

    for( int mr=0; mr<numberOfMaterialRegions; mr++ )
    {
        DispersiveMaterialParameters & dmp = dmpVector[mr]; 

        const IntegerArray & Np = dmp.getBianisotropicNp();
        const RealArray & biPar = dmp.getBianisotropicGDMParameters();

        for( int k1=0; k1<6; k1++ )
        {
            for( int k2=0; k2<6; k2++ )
            {
      	assert( Np(k1,k2)>=0 && Np(k1,k2)<NpMax );

                Npta(k1,k2,mr)= Np(k1,k2);
      	
      	for( int n=0; n<Np(k1,k2); n++ )
      	{
        	  for( int m=0; m<4; m++ )
        	  {
          	    gdmVar(m,n,k1,k2,mr)=biPar(m,n,k1,k2);
        	  }
      	}
            }
        }
        
    }
    
#undef gdmVar  
}

}

#define evalUserDefinedKnownSolution EXTERN_C_NAME(evaluserdefinedknownsolution)


// ==================================================================================
//   Optimized update for RK stages
// 
// Perform the update: 
//    uNew(I1,I2,I3,C) = uOld(I1,I2,I3,C) + SUM_k=1^{numTerms-1} ctk * utk(I1,I2,I3,Ct) 
//
//   maskOption=0 : assign points where  mask>0 otherwise set u2=u1
//             =1 : assign all points
// ==================================================================================




// fourth order dissipation 2D: ***** NOTE: this is minus of the 4th difference:  -(D+D-)^2 *********
#define FD4_2D(u,i1,i2,i3,c) (    -( u(i1-2,i2,i3,c)+u(i1+2,i2,i3,c)+u(i1,i2-2,i3,c)+u(i1,i2+2,i3,c) )   +4.*( u(i1-1,i2,i3,c)+u(i1+1,i2,i3,c)+u(i1,i2-1,i3,c)+u(i1,i2+1,i3,c) ) -12.*u(i1,i2,i3,c) )

// fourth order dissipation 3D:
#define FD4_3D(u,i1,i2,i3,c) (    -( u(i1-2,i2,i3,c)+u(i1+2,i2,i3,c)+u(i1,i2-2,i3,c)+u(i1,i2+2,i3,c)+u(i1,i2,i3-2,c)+u(i1,i2,i3+2,c) )   +4.*( u(i1-1,i2,i3,c)+u(i1+1,i2,i3,c)+u(i1,i2-1,i3,c)+u(i1,i2+1,i3,c)+u(i1,i2,i3-1,c)+u(i1,i2,i3+1,c) ) -18.*u(i1,i2,i3,c) )

#define FOR_3(i1,i2,i3,I1,I2,I3) for( i3=I3.getBase(); i3<=I3.getBound(); i3++ )  for( i2=I2.getBase(); i2<=I2.getBound(); i2++ )  for( i1=I1.getBase(); i1<=I1.getBound(); i1++ )  

#define FN(m) fn[m+numberOfFunctions*(grid)]

// =======================================================================================================
// Macro: Compute the RHS forcing and/or curvilinear operators (e.g. conservative-form of the operators)
//     Note: Some curvilinear opertaors are now computed in advOpt if useCurvilinearOptNew==true
// =======================================================================================================


// ==============================================================================================
// Macro: call the optimized advance routine 
// ==============================================================================================


// =================================================================================
/// \brief Advance a time-step : Bianistropic Maxwell's Equations
// =================================================================================
void Maxwell::
advanceBAMX(  int numberOfStepsTaken, int current, real t, real dt )
{
    checkArrays("advanceBAMX:start");

//    printF("advanceBAMX: t=%e current=%i, numberOfFunctions=%i, numberOfTimeLevels=%i\n",t,
//        current,numberOfFunctions,numberOfTimeLevels);
    
    assert( cgp!=NULL );
    CompositeGrid & cg= *cgp;
    const int numberOfDimensions = cg.numberOfDimensions();
    const int numberOfComponentGrids = cg.numberOfComponentGrids();
    const int & solveForAllFields = dbase.get<int>("solveForAllFields");
    const int & maxNumberOfPolarizationComponents = parameters.dbase.get<int>("maxNumberOfPolarizationComponents");

    cgmxPointer = this; // for getGMDParameters 
    
    Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
    Range C(ex,hz);
    const int numberOfComponents=cgfields[0][0].getLength(3);
    Range Ca = numberOfComponents; // includes dispersion variables 
    const int prev = (current-1+numberOfTimeLevels) % numberOfTimeLevels;
    const int next = (current+1) % numberOfTimeLevels;

    const real cMax=max(cGrid);
    const BoundaryForcingEnum & boundaryForcingOption =dbase.get<BoundaryForcingEnum>("boundaryForcingOption");
    const int & useSosupDissipation = parameters.dbase.get<int>("useSosupDissipation");
    const real & sosupParameter = parameters.dbase.get<real>("sosupParameter");    // scaling of sosup dissipation
    const int & sosupDissipationOption = parameters.dbase.get<int>("sosupDissipationOption"); 
    const int & sosupDissipationFrequency = parameters.dbase.get<int>("sosupDissipationFrequency"); 

    const IntegerArray & totalNumberOfPolarizationComponents =
        parameters.dbase.get<IntegerArray>("totalNumberOfPolarizationComponents");

  // currently there is no reason to update the parallel ghost for P
    const bool updatePolarizationParallelGhost=false;
    
  // We may be able to avoid some parallel ghost updates if there is only one grid 
  // Maybe not if there are far field BC's 
    const bool updateFieldParallelGhost = true; // numberOfComponentGrids>1;
    

    const int dw = max(cg[0].discretizationWidth()); // discretization width

    const int & extrapolateInterpolationNeighbours = dbase.get<int>("extrapolateInterpolationNeighbours");

    const int & orderOfRungeKutta = dbase.get<int>("orderOfRungeKutta");

    const int & useSuperGrid = parameters.dbase.get<int>("useSuperGrid");
    RealArray *etaxSuperGrid = useSuperGrid ? parameters.dbase.get<RealArray*>("etaxSuperGrid" ) : NULL;
    RealArray *etaySuperGrid = useSuperGrid ? parameters.dbase.get<RealArray*>("etaySuperGrid" ) : NULL;
    RealArray *etazSuperGrid = useSuperGrid ? parameters.dbase.get<RealArray*>("etazSuperGrid" ) : NULL;

  // -- now done in setupGridFunctions
  // // -- We normally extrapolate interpolation neighbours for Sosup dissipation which uses a wider stencil
  // //    the stencil is wider. This is not necessary if the grid was generated with more
  // //    layers of interpolation points. 

  // if( dw > orderOfAccuracyInSpace+1 )
  //   extrapolateInterpolationNeighbours=false; // *wdh* added this check, June 15, 2016
  // else if( useSosupDissipation )
  //   extrapolateInterpolationNeighbours=true;

    
    if( useSosupDissipation && sosupDissipationOption==1 )
    {
    // 2 stage sosup dissipation assumes 3 levels 
        assert( numberOfTimeLevels==3 );
    }
    

 // *** add higher order dissipation that requires interpolation:
//  bool useComputeArtificialDissipation=(artificialDissipation>0. || artificialDissipationCurvilinear>0. ) && 
//    orderOfArtificialDissipation > orderOfAccuracyInSpace && orderOfArtificialDissipation>6;

//   if( useComputeArtificialDissipation )
//   {
//     computeDissipation( current,t,dt );
//   }

    if( applyFilter &&  (numberOfStepsTaken % filterFrequency) ==0  )
    {
    // apply high order filter to u[current]
        addFilter( current,t,dt );
    }
    
  // --- arrays holding external forcings : 
    const bool & useNewForcingMethod= dbase.get<bool>("useNewForcingMethod");
    const int & numberOfForcingFunctions= dbase.get<int>("numberOfForcingFunctions"); 
    int & fCurrent = dbase.get<int>("fCurrent");         // forcingArray[fCurrent] : current forcing
    realArray *& forcingArray = dbase.get<realArray*>("forcingArray");  

    sizeOfLocalArraysForAdvance=0.;

  // --------------------------------------------------------------
  // -------------------- ADVANCE GRIDS ---------------------------
  // --------------------------------------------------------------

  // We may option turn off dissipation on some grids:
    RealArray & useDissipation = parameters.dbase.get<RealArray>("useDissipation");
    const bool selectiveDissipation = useDissipation.getLength(0)==cg.numberOfComponentGrids();

  // *** FIX ME ***
  // Runge Kutta Weights
  //     k_i = f( t+ c_i dt , y^n + dt*SUM_j a_ij k j
  //     y^{n+1} = y^n + dt* SUM_i b_i k_i 
  // 
    int numberOfRungeKuttaStages=orderOfRungeKutta; 
  // if( timeSteppingMethod==rungeKutta )
  //   printF(" @@@@@@@@@@@ numberOfRungeKuttaStages=%d @@@@@@@@@@@@@\n",numberOfRungeKuttaStages);
    
    const int maxRKStages=4;
    RealArray aRK(maxRKStages,maxRKStages), bRK(maxRKStages), cRK(maxRKStages+1);
    aRK=0.; bRK=0.; cRK=0.;
    

    if( numberOfRungeKuttaStages==1 )
    {
    // forward Euler
        cRK(0)=0.; cRK(1)=1.;
        bRK(0)=1.;
    }
    else if( numberOfRungeKuttaStages==2 )
    {
    // RK2: explicit trapezoidal rule 
    //   k1 = f(y^n,t)
    //   k2 = f(y^n+dt*k1,t+dt)
    //   y^{n+1} = y^n + (dt/2)*( k1 + k2 )
        
        aRK(1,0)=1.;
        cRK(0)=0.;  cRK(1)=1.;  cRK(2)=1.;
        bRK(0)=.5;  bRK(1)=.5;
    }
    else if( numberOfRungeKuttaStages==3 )
    {
    // SSP-RK3 
    //   k1 = f(y^n,t)
    //   k2 = f(y^n+dt*k1,t+dt)
    //   k3 = f(y^n+ (dt/4)*(k1+k2),t+dt/2)
    //   y^{n+1} = y^n + (dt/6)*( k1 + k2 + 4*k3 )
      
        aRK(1,0)=1.;
        aRK(2,0)=.25; aRK(2,1)=.25;
        
        cRK(0)=0.;     cRK(1)=1.;    cRK(2)=.5;     cRK(3)=1.;
        bRK(0)=1./6.;  bRK(1)=1./6.; bRK(2)=4./6.;


    }
    else if( numberOfRungeKuttaStages==4 )
    {
    // RK4:
    //   k1 = f(y^n,t)
    //   k2 = f(y^n+(dt/2)*k1,t+dt/2)
    //   k3 = f(y^n+(dt/2)*k2,t+dt/2)
    //   k4 = f(y^n+    dt*k3,t+dt)
    //   y^{n+1} = y^n + (dt/6)*( k1 + 2*k2 + 2*k3 + k4 )

        aRK(1,0)=.5;
        aRK(2,0)=.0; aRK(2,1)=.5;
        aRK(3,0)=.0; aRK(3,1)=0.; aRK(3,2)=1.;
        
        cRK(0)=0.;    cRK(1)=.5;    cRK(2)=.5;     cRK(3)=1.;     cRK(4)=1.;
        bRK(0)=1./6.; bRK(1)=2./6.; bRK(2)=2./6.;  bRK(3)=1./6.;

    }
    else
    {
        OV_ABORT("Unexpected numberOfRungeKuttaStages -- finish me");
    }
    
    

  // Here is a list of the stages in the multi-stage FD algorithm 
    if( !dbase.has_key("stageInfoList") )
    {
        std::vector<StageOptionEnum> & stageInfoList = dbase.put<std::vector<StageOptionEnum> >("stageInfoList");
    // -- defaults:
        if( useSosupDissipation && sosupDissipationOption==1 )
        {
      // default two stage algorithm 
            int numberOfStages=2;
            stageInfoList.resize(numberOfStages);
            stageInfoList[0]=StageOptionEnum( addDissipationInStage );
            stageInfoList[1]=StageOptionEnum( updateInteriorInStage | applyBCInStage );
        }
        else if( timeSteppingMethod==rungeKutta )
        {
      // printF(">>> Set stages for RK%d\n",orderOfRungeKutta);

            int numStages=numberOfRungeKuttaStages;
            const real adc = artificialDissipation;   // FIX ME 
            if( adc>0. ) numStages++;

            stageInfoList.resize(numStages);
            for( int k=0; k<numStages; k++ )
      	stageInfoList[k]=StageOptionEnum( updateInteriorInStage | applyBCInStage );

      // Final stage -- add dissipation
            if( adc>0. )
                stageInfoList[numberOfRungeKuttaStages]=StageOptionEnum( addDissipationInStage | applyBCInStage );


      // finish me ... add dissipation on last stage


        }
        else
        {
      // default single stage algorithm 
            int numberOfStages=1;
            stageInfoList.resize(numberOfStages);
      // Oct 22 -- turn off dissipation in ME scheme for now
      // stageInfoList[0]=StageOptionEnum(updateInteriorInStage | addDissipationInStage | applyBCInStage );
            stageInfoList[0]=StageOptionEnum(updateInteriorInStage | applyBCInStage );
        }
        
        
    }
    std::vector<StageOptionEnum> & stageInfoList = dbase.get<std::vector<StageOptionEnum> >("stageInfoList");



  // -- we may advance in two stages for sosup dissipation --
  // const int numberOfStages= (useSosupDissipation && sosupDissipationOption==1 )? 2 : 1;
    const int numberOfStages= stageInfoList.size();

  // +++++++++++++++++++++++++ START STAGES +++++++++++++++++++++++++++++-
    for( int stage=0; stage<numberOfStages; stage ++ )
    { 
  
        StageOptionEnum & stageInfo = stageInfoList[stage];
        const bool computeUt      = stageInfo & computeUtInStage;
        const bool updateInterior = stageInfo & updateInteriorInStage;
        bool addDissipation = stageInfo & addDissipationInStage;
        const bool applyBC        = stageInfo & applyBCInStage;
        
        const bool lastStage = stage==numberOfStages-1;

        real tForce = t;       // get forcing at this time 
        real tBC    = t + dt;  // apply BC's at this time  
        if( timeSteppingMethod==rungeKutta && stage<numberOfRungeKuttaStages )
        {
            tForce = t + cRK(stage)*dt;
            tBC    = t + cRK(stage+1)*dt;
        }
        
    // check me 
        if( numberOfStages>1 && addDissipation && ( numberOfStepsTaken % sosupDissipationFrequency != 0 ) )
            continue;  // skip sosup dissipation 

        if( numberOfStepsTaken<=1 )
        {
            printF("\n +++++ FDBA t=%9.3e STAGE %i: computeUt=%i updateInterior=%i addDissipation=%i applyBCs=%i"
           	     " tBC=%9.3e tForce=%9.3e\n",
                          t,stage,(stageInfo & computeUtInStage ? 1:0),(stageInfo & updateInteriorInStage ? 1 : 0),
                          (stageInfo & addDissipationInStage ? 1:0),(stageInfo & applyBCInStage ? 1:0),tBC,tForce);

            fprintf(pDebugFile,"\n +++++ FDBA t=%9.3e STAGE %i: computeUt=%i updateInterior=%i addDissipation=%i applyBCs=%i"
           	     " tBC=%9.3e tForce=%9.3e\n",
                          t,stage,(stageInfo & computeUtInStage ? 1:0),(stageInfo & updateInteriorInStage ? 1 : 0),
                          (stageInfo & addDissipationInStage ? 1:0),(stageInfo & applyBCInStage ? 1:0),tBC,tForce);

        }
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            real time0=getCPU();

            bool turnOffDissipation = selectiveDissipation && !useDissipation(grid);
            if( turnOffDissipation )
            {
                if( addDissipation && !updateInterior )
                {
                    if( t<= 3.*dt && addDissipation )
                        printF("--ADV-- Skip this dissipation stage for grid=%i (%s)\n",grid,(const char*)cg[grid].getName());

          // skip this step 
                    continue;
                }
                else if( addDissipation && updateInterior )
                {
                    if( t<= 3.*dt && addDissipation )
                        printF("--ADV-- Turn off dissipation on grid=%i (%s)\n",grid,(const char*)cg[grid].getName());

          // turn off dissipation
                    addDissipation=false;
                }
                
            }
            
      // --- Lookup info for the dispersion model ---
            const int domain = cg.domainNumber(grid);
            DispersiveMaterialParameters & dmp = getDomainDispersiveMaterialParameters(domain);
            const int numberOfPolarizationVectors = dmp.numberOfPolarizationVectors;      

      // Each grid may or may not have dispersion model: 
      // const DispersionModelEnum localDispersionModel = numberOfPolarizationVectors>0 ? dispersionModel : noDispersion;
            const DispersionModelEnum localDispersionModel = totalNumberOfPolarizationComponents(grid)>0
                                          ? dispersionModel : noDispersion;
            
      // printF("***** localDispersionModel=%d\n",(int)localDispersionModel);

            MappedGrid & mg = cg[grid];
            assert( mgp==NULL || op!=NULL );
            MappedGridOperators & mgop = mgp!=NULL ? *op : (*cgop)[grid];

            getIndex(mg.gridIndexRange(),I1,I2,I3);

    
            getBoundsForPML( mg,Iv );

            realMappedGridFunction & fieldPrev    =mgp!=NULL ? fields[prev]    : cgfields[prev][grid];
            realMappedGridFunction & fieldCurrent =mgp!=NULL ? fields[current] : cgfields[current][grid];
            realMappedGridFunction & fieldNext    =mgp!=NULL ? fields[next]    : cgfields[next][grid];


            realArray & um = fieldPrev;
            realArray & u  = fieldCurrent;
            realArray & un = fieldNext;

            if( debug & 8 )
            {
                Communication_Manager::Sync();
                fPrintF(debugFile," **** start of advanceBAMX stage=%d, t=%8.2e\n",stage,t);
        // display(um,sPrintF("um start of advanceBAMX, t=%8.2e",t),debugFile,"%8.2e ");
                display(u,sPrintF("u start of advanceBAMX, stage=%d, t=%8.2e",stage,t),debugFile,"%8.2e ");
                Communication_Manager::Sync();
            }

      // 2D TEz mode:
      //   (Ex).t = (1/eps)*[  (Hz).y ]
      //   (Ey).t = (1/eps)*[ -(Hz).x ]
      //   (Hz).t = (1/mu) *[ (Ex).y - (Ey).x ]

            c = cGrid(grid);
            eps = epsGrid(grid);
            mu = muGrid(grid);
        
            if( numberOfStepsTaken<1 ) 
                printF(" advanceBAMX:INFO eps,mu,c=%8.2e %8.2e %8.2e for grid=%i (%s) \n",eps,mu,c,grid,
                              (const char*)cg[grid].getName());
        
      // --- Sanity check: We should now have alphaP=1/eps (or alphaP=0 for testing) ----
            if( localDispersionModel!=noDispersion )
            {
                if( dmp.alphaP !=0. && fabs(dmp.alphaP*eps -1.) > REAL_EPSILON*100. )
                {
                    printF("advance: ERROR: grid=%i: alphaP=%.3e is NOT equal to 1/eps=%.3e -- correct this! Set alphaP=-1 to get default.\n",grid,dmp.alphaP,1./eps);
                    OV_ABORT("error");
                }
                
            }
            

            const bool isRectangular=mg.isRectangular();

      // -------------------------------------------------------------------
      // ------- Compute the RHS forcing and/or curvilinear operators ------
      // -------------------------------------------------------------------
      //     Note: Some curvilinear operators are now computed in advOpt if useCurvilinearOptNew==true
      // const real dtsq=dt*dt; 
            const real csq=c*c;
      // const real cdtsq=c*c*dt*dt;
      // real adc = isRectangular ? c*artificialDissipation : c*artificialDissipationCurvilinear;  
            real adc = artificialDissipation; // ** FIX ME **
            const bool addForcing = forcingOption!=noForcing; 
            bool useCurvilinearOpt=true && !isRectangular;   // use advOpt to advance curvilinear grids given the RHS
      // useCurvilinearOptNew: if true, use advOpt to advance full equations in curvilinear case
      //                       if false, evaluate RHS for curvilinear below and then use advOpt to update solution
      // bool useCurvilinearOptNew = !isRectangular && !useConservative && numberOfDimensions==2 ;
      // *wdh* Feb 25, 2018 -- dispersion model implements optimized curvilinear non-conservative
            bool useCurvilinearOptNew = (
                ( localDispersionModel!=noDispersion && !useConservative ) || 
                ( (!isRectangular  && 
           !useConservative && 
                      (orderOfAccuracyInSpace==4 || orderOfAccuracyInSpace==6) )
                    ) );
            realArray f; // *** SAVE FORCING HERE ***
            if( updateInterior )
            {
                if( addForcing || !useCurvilinearOptNew || (useNewForcingMethod && addForcing) )
                {
          // --- allocate temp space for the forcing ---
                    Index D1,D2,D3;
                    getIndex(mg.dimension(),D1,D2,D3);
                    f.partition(mg.getPartition());
                    f.redim(D1,D2,D3,Ca);  // could use some other array for work space ??
                }
                if( useNewForcingMethod && addForcing )
                { // *new way* 2015/05/18 
                    real timef = getCPU();
                    assert( forcingArray !=NULL );
                    const int fNext = (fCurrent+1) % numberOfForcingFunctions;
                    printF("--MX-ADVS-- evaluate external forcing: t=%9.3e, fCurrent=%i, fNext=%i, (%i)\n",t,
                     	   fCurrent,fNext,numberOfForcingFunctions);
                    realArray & fa = forcingArray[grid];
                    realArray & fb = f;  // we re-use f here for work-space 
                    OV_GET_SERIAL_ARRAY(real,fa,faLocal);
                    OV_GET_SERIAL_ARRAY(real,fb,fbLocal);
                    int includeGhost=1;
                    bool ok = ParallelUtility::getLocalArrayBounds(fb,fbLocal,I1,I2,I3,includeGhost);
                    const int option=1;  // do not append forcing to the "f" array 
                    getForcing( next, grid,fb,t+dt,dt,option );  // **NOTE: get forcing at t+dt 
                    if( ok )
                        faLocal(I1,I2,I3,C,fNext)=fbLocal(I1,I2,I3,C);  // save in fa array
          // faLocal(I1,I2,I3,C,fCurrent)=f;  // *** TEST
          // printF("--MX-ADVR-- max( faLocal(fCurrent)-f)=%8.2e\n",max(fabs(faLocal(I1,I2,I3,C,fCurrent)-f(I1,I2,I3,C))));
                    timing(timeForForcing) += getCPU()-timef;
                }
                if( addForcing || !useCurvilinearOptNew )
                {
                    int option=1;  
                    if( !isRectangular && !useCurvilinearOptNew )
                    {
                        if( timeSteppingMethod == modifiedEquationTimeStepping && orderOfAccuracyInTime>=4 )
                        {
      	// Compute the square of the spatial operator
                  	assert( numberOfFunctions>=1 && fn!=NULL );
                  	const int m0=currentFn;
                  	realArray & lapSq = FN(m0);  
                  	Index J1,J2,J3;
                  	const int extra=1; // orderOfAccuracyInSpace/2-1;
                  	getIndex(mg.gridIndexRange(),J1,J2,J3,extra);
                  	mgop.setOrderOfAccuracy(orderOfAccuracyInSpace-2);
                  	mgop.derivative(MappedGridOperators::laplacianOperator,u,f,J1,J2,J3,C);  // *** use f as a temporary
                            #ifdef USE_PPP
                        	  f.updateGhostBoundaries();
                            #endif
      	// display(f,sPrintF("f=lap(order=2) t=%e processor=%i",t,myid),debugFile,"%6.2f ");
                  	mgop.derivative(MappedGridOperators::laplacianOperator,f,lapSq,I1,I2,I3,C);
                  	mgop.setOrderOfAccuracy(orderOfAccuracyInSpace);
                  	lapSq(I1,I2,I3,C)*=csq*csq;
      	// display(lapSq,sPrintF("lapSq t=%e processor=%i",t,myid),debugFile,"%6.2f ");
      	// printF(" max(fabs(lapSq))=%8.2e min=%8.2e\n",max(fabs(lapSq(I1,I2,I3,C))),min(fabs(lapSq(I1,I2,I3,C))));
                        }
            // compute laplacian for curvilinear grids
                        if( t<3.*dt )
                        {
                  	printF("--MX-- advStr: compute laplacian for curvilinear grids useConservative=%i\n",(int)useConservative);
                        }
                        mgop.derivative(MappedGridOperators::laplacianOperator,u,f,I1,I2,I3,C);
            // * mgop.derivative(MappedGridOperators::laplacianOperator,u,f,I1,I2,I3);
                        f(I1,I2,I3,C)*=csq;
            // if( localDispersionModel != noDispersion )
            // { // set f for dispersive modes to zero 
            //   Range P(pxc,pxc+numberOfDimensions-1);
            //   f(I1,I2,I3,P)=0.;
            // }
            // display(f,sPrintF("lap*csq t=%e processor=%i",t,myid),debugFile,"%6.2f ");
            //f = csq*f + (csq*cdt*cdt/12.)*lapSq;  // put all into f
                        option=0;  // append any forcing below to the "f" array
                    }
                    else
                    {
                        option=1;  // do not append forcing to the "f" array 
                    }
                    if( !useNewForcingMethod && addForcing )
                    {
            //kkc getForcing is called from advance but is also called from assignIC and getErrors.
            //    we have to add the timing in external to getForcing to avoid double counting the time
            //    in assignIC and getErrors
                        real timef = getCPU();
                        getForcing( next, grid,f,tForce,dt,option );
                        timing(timeForForcing) += getCPU()-timef;
                    }
                }
            } // end if updateInterior
        
      // -----------------------------------------------
      // ----- call the optimized advance routine ------
      // -----------------------------------------------
                OV_GET_SERIAL_ARRAY(real,um,umLocal);
                OV_GET_SERIAL_ARRAY(real,u,uLocal);
                OV_GET_SERIAL_ARRAY(real,un,unLocal);
                OV_GET_SERIAL_ARRAY(real,f,fLocal);
                #ifdef USE_PPP
                  realSerialArray varDis; 
                  if( useVariableDissipation ) getLocalArrayWithGhostBoundaries((*variableDissipation)[grid],varDis);
                #else
                  const realSerialArray & varDis = useVariableDissipation ? (*variableDissipation)[grid] : uLocal;
                #endif
        // --- Get pointers to arrays for the dispersive model ----
        // We only need polarization arrays if numberOfPolarizationVectors>0
        // const bool getP = numberOfPolarizationVectors>0;
                const bool getP = totalNumberOfPolarizationComponents(grid)>0;
                realMappedGridFunction & pNext= getP ? getDispersionModelMappedGridFunction( grid,next    ) : fieldNext;
                realMappedGridFunction & pCur = getP ? getDispersionModelMappedGridFunction( grid,current ) : fieldCurrent;
                realMappedGridFunction & pPrev= getP ? getDispersionModelMappedGridFunction( grid,prev )    : fieldPrev;
                OV_GET_SERIAL_ARRAY(real,pNext,pnLocal);
                OV_GET_SERIAL_ARRAY(real, pCur,pLocal);
                OV_GET_SERIAL_ARRAY(real,pPrev,pmLocal);
                real *pnptr=pnLocal.getDataPointer();
                real *pptr = pLocal.getDataPointer();
                real *pmptr=pmLocal.getDataPointer();
        // ::display(pLocal,"pLocal","%5.2f ");
            /* -- OLD 
                real *pnptr=unLocal.getDataPointer();
                real *pptr = uLocal.getDataPointer();
                real *pmptr=umLocal.getDataPointer(); // set defaults if not used  -- just point to u 
                if( numberOfPolarizationVectors>0 )
                {
          // --- Get grid functions for the dispersive model ----
          //  *** FIX ME *** only need 2 levels I think 
                    realMappedGridFunction & pNext= getDispersionModelMappedGridFunction( grid,next );
                    realMappedGridFunction & pCur = getDispersionModelMappedGridFunction( grid,current );
                    realMappedGridFunction & pPrev= getDispersionModelMappedGridFunction( grid,prev );
                    OV_GET_SERIAL_ARRAY(real,pNext,pnLocal);
                    OV_GET_SERIAL_ARRAY(real, pCur,pLocal);
                    OV_GET_SERIAL_ARRAY(real,pPrev,pmLocal);
                    pnptr=pnLocal.getDataPointer();
                    pptr = pLocal.getDataPointer();
                    pmptr=pmLocal.getDataPointer();
            // ::display(pLocal,"pLocal");
                }
                --- */
                bool ok = ParallelUtility::getLocalArrayBounds(u,uLocal,I1,I2,I3);
                real timeAdv=getCPU();
        // In some cases we combine the artificial dissipation loop with the main loop
                int combineDissipationWithAdvance = adc>0. && isRectangular && 
          !useVariableDissipation &&
                    timeSteppingMethod==modifiedEquationTimeStepping &&
                    orderOfAccuracyInSpace==4 && orderOfAccuracyInTime==4;
        // combineDissipationWithAdvance=0;
                const int useWhereMask = numberOfComponentGrids>1;
       // const bool updateSolution = updateInterior;
       // const bool updateDissipation = addDissipation;
                int gridType = isRectangular? 0 : 1;
                int option=(isRectangular || useCurvilinearOpt) ? 0 : 1;   // 0=Maxwell+AD, 1=AD
                int ipar[]={option,
                        	      gridType,
                        	      orderOfAccuracyInSpace,
                        	      orderOfAccuracyInTime,
                        	      (int)addForcing,
                        	      orderOfArtificialDissipation,
                        	      ex,ey,ez,hx,hy,hz,
                        	      int(solveForElectricField),
                        	      int(solveForMagneticField),
                        	      useWhereMask,
                        	      (int)timeSteppingMethod,
                        	      (int)useVariableDissipation,
                        	      (int)useCurvilinearOptNew,
                        	      (int)useConservative,
                        	      combineDissipationWithAdvance,
                        	      (int)useDivergenceCleaning, 
                        	      (int)useNewForcingMethod,
                        	      numberOfForcingFunctions,
                        	      fCurrent,
                        	      localDispersionModel,         // ipar(24) 
                        	      pxc,pyc,pzc, 
                                        totalNumberOfPolarizationComponents(grid),  // ipar(28)
                        	      grid,                         // ipar(29)
                                        0, 0,0,0,  // for future use
                    // qxc,qyc,qzc, 
                    // rxc,ryc,rzc,
                        	      useSosupDissipation,
                        	      sosupDissipationOption,
                        	      updateInterior,
                        	      addDissipation,
                                        computeUt,                  // ipar(38)
                                        (int)forcingOption,
                                        (int)solveForAllFields,     // ipar(40)
                                        (int)dmp.getMaterialType(), // ipar(41)
                        	      numberOfMaterialRegions,    // ipar(42) 
                                        useSuperGrid                // ipar(43)
                                      };  //
                real dx[3]={1.,1.,1.};
                if( isRectangular )
                    mg.getDeltaX(dx);
                real rpar[35];
                rpar[ 0]=c;
                rpar[ 1]=dt;
                rpar[ 2]=dx[0];
                rpar[ 3]=dx[1];
                rpar[ 4]=dx[2];
                rpar[ 5]=adc;
                rpar[ 6]=divergenceDamping;
                rpar[ 7]=mg.gridSpacing(0);
                rpar[ 8]=mg.gridSpacing(1);
                rpar[ 9]=mg.gridSpacing(2);
                rpar[10]=eps;
                rpar[11]=mu;
                rpar[12]=kx; // for plane wave scattering
                rpar[13]=ky;
                rpar[14]=kz;
                rpar[15]=sigmaEGrid(grid);
                rpar[16]=sigmaHGrid(grid);
                rpar[17]=divergenceCleaningCoefficient;
                rpar[18]=tForce;     // **check me : Nov 12, 2019 ***
                rpar[19]= (real &)tz;  // twilight zone pointer
                rpar[20]=0.;  // return cpu for dissipation
        // Dispersive material parameters
                rpar[21]=dmp.gamma;
                rpar[22]=dmp.omegap;
        // ADD THIS AS AN OPTION
                if( false && isRectangular )
                    rpar[23]=0.;      // TEST-- turn off sosup diss on rectangular grids
                else
                    rpar[23]=sosupParameter;
        // New way: GDM coefficients
                const RealArray & mp = dmp.modelParameters;
                rpar[24]=dmp.alphaP;
                rpar[25]=mp(0,0); // a0 
                rpar[26]=mp(1,0); // a1 
                rpar[27]=mp(2,0); // b0 
                rpar[28]=mp(3,0); // b1 
                int ierr=0;
                real *umptr=umLocal.getDataPointer();
                real *uptr =uLocal.getDataPointer();
                real *unptr=unLocal.getDataPointer();
                real *ut1ptr = uptr; 
                real *ut2ptr = uptr; 
                real *ut3ptr = uptr; 
                real *ut4ptr = uptr; 
                real *ut5ptr = uptr; 
                real *ut6ptr = uptr; 
                real *ut7ptr = uptr; 
                real *fptr   = (addForcing || useCurvilinearOpt) ? fLocal.getDataPointer() : uptr;
                real *kvp[10];  // pointers to RK stage slope functions k1,k2,k3,k4 
        // external forcings at different time levels are stored here: 
                real *faptr = fptr;
                if( useNewForcingMethod && addForcing )
                {
                    OV_GET_SERIAL_ARRAY(real,forcingArray[grid],faLocal);
                    faptr   = faLocal.getDataPointer();
                }
                assert( !useVariableDissipation || variableDissipation!=NULL );
                real *pVarDis = useVariableDissipation ? varDis.getDataPointer() : uptr;
        // Dispersion model: polarization vectors 
                ut1ptr = pmptr;
                ut2ptr = pptr;
                ut3ptr = pnptr;
                if( timeSteppingMethod==modifiedEquationTimeStepping )
                {
          // --- WORK-SPACE FOR Modified Equation Time Stepping ---
                    if( useConservative )
                    {
            // one work space array needed
                        assert( numberOfFunctions>=1 && fn!=NULL );
                        const int m0=currentFn;
                        OV_GET_SERIAL_ARRAY(real,FN(m0),f0Local); ut1ptr=f0Local.getDataPointer();
                    }
                    else
                    {
            // one work space array needed
                        assert( numberOfFunctions>=1 && fn!=NULL );
                        const int m0=currentFn;
                        OV_GET_SERIAL_ARRAY(real,FN(m0),f0Local); ut1ptr=f0Local.getDataPointer();
                    }
                }
                else if( timeSteppingMethod==rungeKutta )
                {
                    assert( numberOfFunctions>=4 );
                    for( int k=0; k<numberOfRungeKuttaStages; k++ )
                    {
                        OV_GET_SERIAL_ARRAY(real,FN(k),f0Local); kvp[k]=f0Local.getDataPointer();
            // cout << "kvp[" << k << "]=" << kvp[k] << endl;
                    }
                }
                else  // MOL time-stepping 
                {
          //   ****** FIX ME *****
                    OV_ABORT("error");
            /* ---- Not needed: 
                    if( orderOfAccuracyInTime>=4 )
                    {
                        assert( numberOfFunctions>=3 && fn!=NULL );
                        const int m0=currentFn, m1=(m0+1)%numberOfFunctions, m2=(m1+1)%numberOfFunctions;
                        OV_GET_SERIAL_ARRAY(real,FN(m0),f0Local); ut1ptr=f0Local.getDataPointer();
                        OV_GET_SERIAL_ARRAY(real,FN(m1),f1Local); ut2ptr=f1Local.getDataPointer();
                        OV_GET_SERIAL_ARRAY(real,FN(m2),f2Local); ut3ptr=f2Local.getDataPointer();
                        if( orderOfAccuracyInTime>=6 && timeSteppingMethod==stoermerTimeStepping )
                        {
                  	assert( numberOfFunctions>=5 );
                  	const int m3=(m2+1)%numberOfFunctions, m4=(m3+1)%numberOfFunctions;
                  	OV_GET_SERIAL_ARRAY(real,FN(m3),f3Local); ut4ptr=f3Local.getDataPointer();
                  	OV_GET_SERIAL_ARRAY(real,FN(m4),f4Local); ut5ptr=f4Local.getDataPointer();
                  	if( orderOfAccuracyInTime>=8 && timeSteppingMethod==stoermerTimeStepping )
                  	{
                    	  assert( numberOfFunctions>=7 );
                    	  const int m5=(m4+1)%numberOfFunctions, m6=(m5+1)%numberOfFunctions;
                    	  OV_GET_SERIAL_ARRAY(real,FN(m5),f5Local); ut6ptr=f5Local.getDataPointer();
                    	  OV_GET_SERIAL_ARRAY(real,FN(m6),f6Local); ut7ptr=f6Local.getDataPointer();
                  	}
                        }
                    }
                    --- */
                }
                intArray & mask = mg.mask();
                OV_GET_SERIAL_ARRAY(int,mask,maskLocal);
                real *rxptr;
                if( isRectangular )
                {
                    rxptr=uptr;
                }
                else
                {
                    OV_GET_SERIAL_ARRAY(real,mg.inverseVertexDerivative(),rxLocal);
                    rxptr = rxLocal.getDataPointer();
                }
                const bool centerNeeded = forcingOption == twilightZoneForcing;
                real *xyptr = uptr;
                if( centerNeeded )
                {
          // Pass the xy array for twilight-zone -- needed for TZ for dispersive models 
                    OV_GET_SERIAL_ARRAY(real,mg.center(),xLocal);
                    xyptr = xLocal.getDataPointer();
                    assert( xyptr!=NULL );
                    ut4ptr = xyptr;   // note: this is ut4ptr for ba solver
                }
                if( numberOfMaterialRegions==0 )
                    numberOfMaterialRegions=1;   // ** FIX ME **
        // --- Get the inverse of the BA material matrix ---
                if( !dbase.has_key("K0Inverse") )
                {
          // This must be a single material problem since the the K0Inverse has not been created by defineMaterialRegions
                    RealArray & K0i = dbase.put<RealArray>("K0Inverse");
                    K0i.redim(6,6,1);
                    RealArray K0iTemp;
                    dmp.getBianisotropicMaterialMatrixInverse( K0iTemp );
                    Range all;
                    K0i(all,all,0)=K0iTemp;
                }
                RealArray & K0i = dbase.get<RealArray>("K0Inverse");
        // RealArray K0i(6,6,numberOfMaterialRegions);
        // RealArray K0iTemp;
        // dmp.getBianisotropicMaterialMatrixInverse( K0iTemp );
        // Range all;
        // K0i(all,all,0)=K0iTemp;
        // if( numberOfMaterialRegions>1 )
        // {
        //   // --- Fill in the material matrix for the embedded material regions ----
        //   std::vector<DispersiveMaterialParameters> & dmpVector = 
        //       dbase.get<std::vector<DispersiveMaterialParameters> >("materialRegionParameters");
        //   printF("materialRegionParameters: dmpVector.size()=%i\n",dmpVector.size());
        //   for( int mr=1; mr<numberOfMaterialRegions; mr++ )
        //   {
        //     // DispersiveMaterialParameters & dmp2 = dmpVector[mr-1];
        //     DispersiveMaterialParameters & dmp2 = dmpVector[mr];  // ** FIX ME ***
        //     dmp2.display();
        //     K0iTemp=0;
        //     dmp2.getBianisotropicMaterialMatrixInverse( K0iTemp );
        //     // ::display(K0iTemp,sPrintF("Material matrix: region %d",mr));
        //     K0i(all,all,mr)=K0iTemp;  
        //   }
        // }
        // ***NOTE*** pmask points to the bodyMask
                if( numberOfMaterialRegions>1 )
                {
                    assert( pBodyMask!=NULL );
                }
                int *pmask = numberOfMaterialRegions>1 ? pBodyMask->getDataPointer() : maskLocal.getDataPointer();
                if( debug & 8 )
                {
                    fprintf(debugFile,"addForcing=%i, useNewForcingMethod=%i\n",(int)addForcing,(int)useNewForcingMethod);
                    display(fLocal,sPrintF("fLocal BEFORE advBA t=%8.2e",t),debugFile,"%8.2e ");
          // display(un,sPrintF("un BEFORE advBA t=%8.2e",t),debugFile,"%8.2e ");
          // display(un,sPrintF("un BEFORE advBA t=%8.2e",t),debugFile,"%16.10e ");
                    display(u,sPrintF("u BEFORE advBA t=%8.2e",t),debugFile,"%16.10e ");
                }
        // int *maskptr = useWhereMask ? maskLocal.getDataPointer() : ipar;
                int *maskptr = maskLocal.getDataPointer(); // *wdh* Jan 5, 2017 -- do this always
                realSerialArray *dis = NULL;
                real *pdis=uptr;
                if( adc>0. && !combineDissipationWithAdvance && ok )
                {
          // create a temp array to hold the artificial dissipation
                    dis = new realSerialArray(uLocal.dimension(0),uLocal.dimension(1),uLocal.dimension(2),uLocal.dimension(3));
                    pdis = dis->getDataPointer();
                    assert( pdis!=NULL );
                    sizeOfLocalArraysForAdvance=max(sizeOfLocalArraysForAdvance,(double)(dis->elementCount()*sizeof(real)));
                }
                real *pEtaxSuperGrid=uptr, *pEtaySuperGrid=uptr, *pEtazSuperGrid=uptr;
                if( useSuperGrid )
                {
                    pEtaxSuperGrid = etaxSuperGrid[grid].getDataPointer();
                    pEtaySuperGrid = etaySuperGrid[grid].getDataPointer();
                    pEtazSuperGrid = etazSuperGrid[grid].getDataPointer();
                }
                if( timeSteppingMethod==rungeKutta )
                {
                    if( stage<numberOfRungeKuttaStages )
                    {
                        if( stage>0 )
                        {
                  	uptr =unptr;         // evaluate time-derivative of this variable 
                  	pptr =pnptr;         // evaluate time-derivative of this variable 
                        }
                        unptr = kvp[stage];    // store  du/dt here
                    }
                    else
                    {
            // dissipation stage
            // uptr =unptr;         // un <- un + diss(un) 
                        uptr =unptr;         // k1 <- un + diss(un) 
                        unptr = kvp[0];
                    }
                    assert( unptr != NULL );
                }
                if( ok )
                {
                    if( combineDissipationWithAdvance )
                    {
                        assert( umptr!=unptr );
                    }
                    advBA(mg.numberOfDimensions(),
                         	       I1.getBase(),I1.getBound(),I2.getBase(),I2.getBound(),I3.getBase(),I3.getBound(),
                         	       uLocal.getBase(0),uLocal.getBound(0),uLocal.getBase(1),uLocal.getBound(1),
                         	       uLocal.getBase(2),uLocal.getBound(2),
                         	       uLocal.getBase(3),uLocal.getBound(3),
                         	       *maskptr,*rxptr,  
                         	       *umptr,*uptr,*unptr, *fptr,
                         	       *faptr,  // forcing at multiple time levels 
                                 	       K0i(0,0),*pmask,
      //	       *ut1ptr,*ut2ptr,*ut3ptr,*ut4ptr,*ut5ptr,*ut6ptr,*ut7ptr,   
                         	       *pmptr,*pptr,*pnptr,*xyptr,
                         	       *pEtaxSuperGrid,*pEtaySuperGrid,*pEtazSuperGrid,   
                         	       mg.boundaryCondition(0,0), *pdis, *pVarDis, ipar[0], rpar[0], ierr );
                    const bool useOptUpdate=true;
                    if( timeSteppingMethod==rungeKutta && stage<numberOfRungeKuttaStages )
                    {
                        OV_GET_SERIAL_ARRAY(real,FN(0),k1Local);
                        OV_GET_SERIAL_ARRAY(real,FN(1),k2Local);
                        OV_GET_SERIAL_ARRAY(real,FN(2),k3Local);
                        OV_GET_SERIAL_ARRAY(real,FN(3),k4Local);
                        const int numPolarizationTerms = maxNumberOfPolarizationComponents*2;
                        Range P=  max(1,numPolarizationTerms);  // range of all p, p.t terms 
                        Range Pt=P+C.length();           // time derivatives are here in the k1,k2,.. arrays 
            // ****** OPTIMIZE THIS *****
                        if( stage==numberOfRungeKuttaStages-1 )
                        {
              // Last RK stage -- perform the quadrature 
                  	if( useOptUpdate )
                  	{
      	  // **new** optimized update 
                    	  int numTerms=1+numberOfRungeKuttaStages;
                    	  real ct1 = dt*bRK(0), ct2=dt*bRK(1), ct3=dt*bRK(2), ct4=dt*bRK(3);
                            {
                                int ierr=0;
                                const int maskOption=0; // assign pts where mask>0
                                int is4 = C.getBase()-C.getBase();  // index offset component values in k1Local, k2Local, etc. 
                                int ipar[]={numTerms,maskOption,I1.getBase(),I1.getBound(),I2.getBase(),I2.getBound(),I3.getBase(),I3.getBound(),
                                                        C.getBase(),C.getBound(), is4}; //
                                real cOld=1.;
                                real rpar[10]={cOld,ct1,ct2,ct3,ct4,0.,0.,0.,0.,0.};
                                updateOptNew(unLocal.getBase(0),unLocal.getBound(0),unLocal.getBase(1),unLocal.getBound(1),
                                                          unLocal.getBase(2),unLocal.getBound(2),unLocal.getBase(3),unLocal.getBound(3),
                                                          *mask.getDataPointer(),  
                                                          *unLocal.getDataPointer(),  *uLocal.getDataPointer(), 
                                                          *k1Local.getDataPointer(),*k2Local.getDataPointer(),*k3Local.getDataPointer(),*k4Local.getDataPointer(),
                                                          *k1Local.getDataPointer(),*k1Local.getDataPointer(),*k1Local.getDataPointer(),*k1Local.getDataPointer(), // not used 
                                                          *k1Local.getDataPointer(),                                                                   // not used 
                                                          ipar[0], rpar[0], ierr );
                            }
                    	  if( dispersionModel != noDispersion )
                                {
                                    int ierr=0;
                                    const int maskOption=0; // assign pts where mask>0
                                    int is4 = Pt.getBase()-P.getBase();  // index offset component values in k1Local, k2Local, etc. 
                                    int ipar[]={numTerms,maskOption,I1.getBase(),I1.getBound(),I2.getBase(),I2.getBound(),I3.getBase(),I3.getBound(),
                                                            P.getBase(),P.getBound(), is4}; //
                                    real cOld=1.;
                                    real rpar[10]={cOld,ct1,ct2,ct3,ct4,0.,0.,0.,0.,0.};
                                    updateOptNew(pnLocal.getBase(0),pnLocal.getBound(0),pnLocal.getBase(1),pnLocal.getBound(1),
                                                              pnLocal.getBase(2),pnLocal.getBound(2),pnLocal.getBase(3),pnLocal.getBound(3),
                                                              *mask.getDataPointer(),  
                                                              *pnLocal.getDataPointer(),  *pLocal.getDataPointer(), 
                                                              *k1Local.getDataPointer(),*k2Local.getDataPointer(),*k3Local.getDataPointer(),*k4Local.getDataPointer(),
                                                              *k1Local.getDataPointer(),*k1Local.getDataPointer(),*k1Local.getDataPointer(),*k1Local.getDataPointer(), // not used 
                                                              *k1Local.getDataPointer(),                                                                   // not used 
                                                              ipar[0], rpar[0], ierr );
                                }
                  	}
      	// **** OLD *****
                            else if( numberOfRungeKuttaStages==1 )
                  	{ // RK1
                    	  unLocal(I1,I2,I3,C) = uLocal(I1,I2,I3,C) + (dt*bRK(0))*k1Local(I1,I2,I3,C);
                                if( dispersionModel != noDispersion )
                                    pnLocal(I1,I2,I3,P) = pLocal(I1,I2,I3,P) + (dt*bRK(0))*k1Local(I1,I2,I3,Pt);
                  	}
                            else if( numberOfRungeKuttaStages==2 )
                  	{ // RK2
                    	  unLocal(I1,I2,I3,C) = uLocal(I1,I2,I3,C) + (dt*bRK(0))*k1Local(I1,I2,I3,C)
                                                             	                                           + (dt*bRK(1))*k2Local(I1,I2,I3,C);
                                if( dispersionModel != noDispersion )
                                    pnLocal(I1,I2,I3,P) = pLocal(I1,I2,I3,P) + (dt*bRK(0))*k1Local(I1,I2,I3,Pt)
                                                                                                                      + (dt*bRK(1))*k2Local(I1,I2,I3,Pt);
                  	}
                            else if( numberOfRungeKuttaStages==3 )
                  	{ // RK3 
                    	  unLocal(I1,I2,I3,C) = uLocal(I1,I2,I3,C) + (dt*bRK(0))*k1Local(I1,I2,I3,C)
                                                             	                                           + (dt*bRK(1))*k2Local(I1,I2,I3,C)
                                                             	                                           + (dt*bRK(2))*k3Local(I1,I2,I3,C);
                                if( dispersionModel != noDispersion )
                                    pnLocal(I1,I2,I3,P) = pLocal(I1,I2,I3,P) + (dt*bRK(0))*k1Local(I1,I2,I3,Pt)
                                                                                                                      + (dt*bRK(1))*k2Local(I1,I2,I3,Pt)
                                                               	                                             + (dt*bRK(2))*k3Local(I1,I2,I3,Pt);
                  	}
                            else if( numberOfRungeKuttaStages==4 )
                  	{ // RK4 
                    	  unLocal(I1,I2,I3,C) = uLocal(I1,I2,I3,C) + (dt*bRK(0))*k1Local(I1,I2,I3,C)
                                                             	                                           + (dt*bRK(1))*k2Local(I1,I2,I3,C)
                                                             	                                           + (dt*bRK(2))*k3Local(I1,I2,I3,C)
                                                                     	                                   + (dt*bRK(3))*k4Local(I1,I2,I3,C);
                                if( dispersionModel != noDispersion )
                                    pnLocal(I1,I2,I3,P) = pLocal(I1,I2,I3,P) + (dt*bRK(0))*k1Local(I1,I2,I3,Pt)
                                                                                                                      + (dt*bRK(1))*k2Local(I1,I2,I3,Pt)
                                                                                                                      + (dt*bRK(2))*k3Local(I1,I2,I3,Pt)
                                                                                                                      + (dt*bRK(3))*k4Local(I1,I2,I3,Pt);
                  	}
                        }
                        else if( useOptUpdate && stage>=0 && stage<=2 )
                        {
              // **new** optimized update 
                  	int numTerms=2+stage;
                  	real ct1 = dt*aRK(stage+1,0), ct2=dt*aRK(stage+1,1), ct3=dt*aRK(stage+1,2), ct4=0.;
                        {
                            int ierr=0;
                            const int maskOption=0; // assign pts where mask>0
                            int is4 = C.getBase()-C.getBase();  // index offset component values in k1Local, k2Local, etc. 
                            int ipar[]={numTerms,maskOption,I1.getBase(),I1.getBound(),I2.getBase(),I2.getBound(),I3.getBase(),I3.getBound(),
                                                    C.getBase(),C.getBound(), is4}; //
                            real cOld=1.;
                            real rpar[10]={cOld,ct1,ct2,ct3,ct4,0.,0.,0.,0.,0.};
                            updateOptNew(unLocal.getBase(0),unLocal.getBound(0),unLocal.getBase(1),unLocal.getBound(1),
                                                      unLocal.getBase(2),unLocal.getBound(2),unLocal.getBase(3),unLocal.getBound(3),
                                                      *mask.getDataPointer(),  
                                                      *unLocal.getDataPointer(),  *uLocal.getDataPointer(), 
                                                      *k1Local.getDataPointer(),*k2Local.getDataPointer(),*k3Local.getDataPointer(),*k4Local.getDataPointer(),
                                                      *k1Local.getDataPointer(),*k1Local.getDataPointer(),*k1Local.getDataPointer(),*k1Local.getDataPointer(), // not used 
                                                      *k1Local.getDataPointer(),                                                                   // not used 
                                                      ipar[0], rpar[0], ierr );
                        }
                  	if( dispersionModel != noDispersion )
                            {
                                int ierr=0;
                                const int maskOption=0; // assign pts where mask>0
                                int is4 = Pt.getBase()-P.getBase();  // index offset component values in k1Local, k2Local, etc. 
                                int ipar[]={numTerms,maskOption,I1.getBase(),I1.getBound(),I2.getBase(),I2.getBound(),I3.getBase(),I3.getBound(),
                                                        P.getBase(),P.getBound(), is4}; //
                                real cOld=1.;
                                real rpar[10]={cOld,ct1,ct2,ct3,ct4,0.,0.,0.,0.,0.};
                                updateOptNew(pnLocal.getBase(0),pnLocal.getBound(0),pnLocal.getBase(1),pnLocal.getBound(1),
                                                          pnLocal.getBase(2),pnLocal.getBound(2),pnLocal.getBase(3),pnLocal.getBound(3),
                                                          *mask.getDataPointer(),  
                                                          *pnLocal.getDataPointer(),  *pLocal.getDataPointer(), 
                                                          *k1Local.getDataPointer(),*k2Local.getDataPointer(),*k3Local.getDataPointer(),*k4Local.getDataPointer(),
                                                          *k1Local.getDataPointer(),*k1Local.getDataPointer(),*k1Local.getDataPointer(),*k1Local.getDataPointer(), // not used 
                                                          *k1Local.getDataPointer(),                                                                   // not used 
                                                          ipar[0], rpar[0], ierr );
                            }
      	// ::display(k1Local,"k1","%6.2f ");
      	// ::display(unLocal,"un","%6.2f ");
                        }
                        else if( stage==0 )
                        { // form the temporary solution use to compute the next slope function k[stage+1] 
              // if( useOptUpdate )
              // {
              //   int numTerms=2;
              //   real ct1 = dt*aRK(stage+1,0), ct2=0., ct3=0., ct4=0.;
              //   updateRK( unLocal, uLocal,k1Local,k2Local,k3Local,k4Local, ct1,ct2,ct3,ct4,numTerms, I1,I2,I3,C,C );
              //   if( dispersionModel != noDispersion )
              //     updateRK( pnLocal, pLocal,k1Local,k2Local,k3Local,k4Local, ct1,ct2,ct3,ct4,numTerms, I1,I2,I3,P,Pt );
              // }
                  	unLocal(I1,I2,I3,C) = uLocal(I1,I2,I3,C) + (dt*aRK(stage+1,0))*k1Local(I1,I2,I3,C);
                  	if( dispersionModel != noDispersion )
                    	  pnLocal(I1,I2,I3,P) = pLocal(I1,I2,I3,P) + (dt*aRK(stage+1,0))*k1Local(I1,I2,I3,Pt);
        	// printF("+++++ STAGE 0 +++++\n");
      	// ::display(k1Local(I1,I2,I3,C),"k1Local(I1,I2,I3,C)");
      	// ::display(unLocal(I1,I2,I3,C),"unLocal(I1,I2,I3,C)");
                        }
                        else if( stage==1 )
                        {
                  	unLocal(I1,I2,I3,C) = uLocal(I1,I2,I3,C) + (dt*aRK(stage+1,0))*k1Local(I1,I2,I3,C)
                                                           	                                         + (dt*aRK(stage+1,1))*k2Local(I1,I2,I3,C);
                            if( dispersionModel != noDispersion )
                                pnLocal(I1,I2,I3,P) = pLocal(I1,I2,I3,P) + (dt*aRK(stage+1,0))*k1Local(I1,I2,I3,Pt)
                                                             	                                           + (dt*aRK(stage+1,1))*k2Local(I1,I2,I3,Pt);
                        }
                        else if( stage==2 )
                        {
                  	unLocal(I1,I2,I3,C) = uLocal(I1,I2,I3,C) + (dt*aRK(stage+1,0))*k1Local(I1,I2,I3,C)
                                                           	                                         + (dt*aRK(stage+1,1))*k2Local(I1,I2,I3,C)
                                                           	                                         + (dt*aRK(stage+1,2))*k3Local(I1,I2,I3,C);
                            if( dispersionModel != noDispersion )
                                pnLocal(I1,I2,I3,P) = pLocal(I1,I2,I3,P) + (dt*aRK(stage+1,0))*k1Local(I1,I2,I3,Pt)
                                                             	                                           + (dt*aRK(stage+1,1))*k2Local(I1,I2,I3,Pt)
                                                             	                                           + (dt*aRK(stage+1,2))*k3Local(I1,I2,I3,Pt);
                        }
                        else 	
                        {
                            OV_ABORT("Unexpected RK stage ?! ");
                        }
                    }
                    else if( timeSteppingMethod==rungeKutta && stage==numberOfRungeKuttaStages  )
                    {
            // dissipation stage
            //   Filtered solution was returned in k1 
            // -- copy un = k1 
                        assert( addDissipation );
                        OV_GET_SERIAL_ARRAY(real,FN(0),k1Local);
                        if( useOptUpdate )
                        {
      	// **new** optimized update 
                  	int numTerms=1;  
                  	real ct1 = 0, ct2=0, ct3=0, ct4=0; // Not used 
                        {
                            int ierr=0;
                            const int maskOption=0; // assign pts where mask>0
                            int is4 = C.getBase()-C.getBase();  // index offset component values in k1Local, k1Local, etc. 
                            int ipar[]={numTerms,maskOption,I1.getBase(),I1.getBound(),I2.getBase(),I2.getBound(),I3.getBase(),I3.getBound(),
                                                    C.getBase(),C.getBound(), is4}; //
                            real cOld=1.;
                            real rpar[10]={cOld,ct1,ct2,ct3,ct4,0.,0.,0.,0.,0.};
                            updateOptNew(unLocal.getBase(0),unLocal.getBound(0),unLocal.getBase(1),unLocal.getBound(1),
                                                      unLocal.getBase(2),unLocal.getBound(2),unLocal.getBase(3),unLocal.getBound(3),
                                                      *mask.getDataPointer(),  
                                                      *unLocal.getDataPointer(),  *k1Local.getDataPointer(), 
                                                      *k1Local.getDataPointer(),*k1Local.getDataPointer(),*k1Local.getDataPointer(),*k1Local.getDataPointer(),
                                                      *k1Local.getDataPointer(),*k1Local.getDataPointer(),*k1Local.getDataPointer(),*k1Local.getDataPointer(), // not used 
                                                      *k1Local.getDataPointer(),                                                                   // not used 
                                                      ipar[0], rpar[0], ierr );
                        }
                        }
                        else
                        {
                            unLocal(I1,I2,I3,C) = k1Local(I1,I2,I3,C);
                        }
                    }
                }
                timeAdv=getCPU()-timeAdv;
                timing(timeForAdvOpt)+=timeAdv;
                if( debug & 8 )
                {
          // OV_GET_SERIAL_ARRAY(real,FN(0),k1Local);
          // display(k1Local,sPrintF("k1Local after advBA, processor=%i before BC's stage=%d, t=%8.2e",
          //		    myid,stage,t),pDebugFile,"%8.2e ");
                    display(unLocal,sPrintF("unLocal after advBA, processor=%i before BC's stage=%d, t=%8.2e",
                                  			    myid,stage,t),pDebugFile,"%8.2e ");
          // display(unLocal,sPrintF("unLocal after advBA, processor=%i before BC's stage=%d, t=%8.2e",
          //     		    Communication_Manager::My_Process_Number,stage,t),pDebugFile,"%8.2e ");
                    display(un,sPrintF("un after advBA, before BC's t=%8.2e",t),debugFile,"%8.2e ");
                }
      //        printF(" p=%i time for advBA=%e I1,I2,I3=[%i,%i][%i,%i][%i,%i]\n",Communication_Manager::My_Process_Number,timeAdv,
      //                   I1.getBase(),I1.getBound(),I2.getBase(),I2.getBound(),I3.getBase(),I3.getBound());
                if( dis!=NULL )
                {
                    delete dis;
                }
                timing(timeForDissipation)+=rpar[20];
                if( isRectangular )   
                    timing(timeForAdvanceRectangularGrids)+=getCPU()-time0;
                else
                    timing(timeForAdvanceCurvilinearGrids)+=getCPU()-time0;
        
            if( applyBC )
            {
        // Is this the right place?
                real timea=getCPU();
                if( debug & 8 )
                {
                    Communication_Manager::Sync();
                    display(unLocal,sPrintF(" Before updateGhostBoundaries: t=%e",t),pDebugFile,"%8.2e ");
                }
            
        // **** at this point we really only need to update interior-ghost points needed for
        //      interpolation or boundary conditions
      	if( updateFieldParallelGhost )
        	  un.updateGhostBoundaries();

                if( updatePolarizationParallelGhost )
                    pNext.updateGhostBoundaries();

                if( debug & 8 )
                {
                    display(unLocal,sPrintF(" After updateGhostBoundaries: t=%e",t),pDebugFile,"%8.2e ");
                }
            
                timing(timeForUpdateGhostBoundaries)+=getCPU()-timea;

                if( debug & 8 )
                {
                    display(u,sPrintF("u after advOpt and updateGhost, stage=%d, t=%8.2e",stage,t),debugFile,"%8.2e ");
                }


                if( debug & 16 )
                {
                    if( addForcing ) 
                        ::display(f,sPrintF("  *** advanceBAMX: Here is the forcing f grid=%i t=%9.3e ********\n",grid,t),
                                            pDebugFile,"%7.4f ");
                    fprintf(pDebugFile," *** advanceBAMX: After advance, before BC, grid=%i t=%9.3e ********\n",grid,t);
                    getErrors( next,tBC,dt );
                }
            } // end if applyBC
        

        } // end grid

        if( applyBC )
        {
      // ------ APPLY BOUNDARY CONDITIONS ------
            if( false && numberOfStages>1 && numberOfStepsTaken<=4 )
            {
                printF(" +++ apply boundary conditions...\n");
            }
            

      //    Communication_Manager::Sync();

            if( orderOfAccuracyInTime>=4 )
            {
                currentFn=(currentFn+orderOfAccuracyInTime-2)%numberOfFunctions;
            }
    
            if( cg.numberOfComponentGrids()>1 )
            {
                real timei=getCPU();
        
                if( debug & 8 )
                    cgfields[next].display(sPrintF("cgfields[next] before interpolate, t=%8.2e",t),debugFile,"%8.2e ");

                cgfields[next].interpolate();
  
                if( dispersionModel !=noDispersion )
                {
                    if( !parameters.dbase.has_key("domainInterpolant") )
                    {
            // Create Interpolants for the separate domains so we can interpolate the Polarization Vectors 

                        Interpolant **& domainInterpolant = parameters.dbase.put<Interpolant**>("domainInterpolant");
                        domainInterpolant = new Interpolant* [cg.numberOfDomains()];
                    
                        for( int domain=0; domain<cg.numberOfDomains(); domain++ )
                        {
                            domainInterpolant[domain]=NULL;
                            const DispersiveMaterialParameters & dmp = getDomainDispersiveMaterialParameters(domain);
                            if( dmp.numberOfPolarizationVectors>0 )
                            {
                                domainInterpolant[domain] = new Interpolant(cg.domain[domain]);
                                domainInterpolant[domain]-> incrementReferenceCount();
                            }
                        }
                        
                    }
                    

                    for( int domain=0; domain<cg.numberOfDomains(); domain++ )
                    {
                        realCompositeGridFunction *ppNext= getDispersionModelCompositeGridFunction( domain,next );
                        if( ppNext!=NULL )
                        {
                            ppNext->interpolate();
                        }
                        
                    }
                    
                }
                

                if( debug & 8 )
                    cgfields[next].display(sPrintF("cgfields[next] after interpolate, t=%8.2e",t),debugFile,"%8.2e ");


                if( projectInterpolation )
                    projectInterpolationPoints( numberOfStepsTaken, next, tBC, dt );

                timing(timeForInterpolate)+=getCPU()-timei;

            }


      // ================= Project Fields =================================

      // Here we project the field to satisfy Gauss's Law  div(E)=0 or div(E) = stuff
      // It seems important to project two steps in a row

      // We project fields here before the BC's -- to do this we use:
      //    --> Boundary values of the field should be accurate
      //    --> the divergence on the boundary is not needed for the projection
      //    --> for 4th-order we do need values on the first ghost line

            const bool projectBeforeBCs=false;
            const bool projectThisStep = projectFields && (
                ( (numberOfStepsTaken % frequencyToProjectFields) < numberOfConsecutiveStepsToProject ) ||
                numberOfStepsTaken < numberOfInitialProjectionSteps );

            if( lastStage && projectBeforeBCs && projectThisStep )
            {

                if( cg.numberOfComponentGrids()==1 )
                    cgfields[next].periodicUpdate();   // this seems to be needed

                project( numberOfStepsTaken, next, tBC , dt );

                if( false )
                {
          // re-apply the BC;s
                    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                    {
                        realMappedGridFunction & fieldNext    =mgp!=NULL ? fields[next]    : cgfields[next][grid];
                        realMappedGridFunction & fieldCurrent =mgp!=NULL ? fields[current] : cgfields[current][grid];

                        int option=0;
                        assignBoundaryConditions( option, grid, tBC, dt, fieldNext, fieldCurrent,current );
                    }
                }
        
            }

      // ================================================================================
      // ============== MATERIAL INTERFACES : STAGE I - BOUNDARY VALUES =================
      // ================================================================================

      // ---- Assign values on the material interfaces BOUNDARY (but not ghost)------
            bool assignInterfaceValues=true;
            bool assignInterfaceGhostValues=false;
            assignInterfaceBoundaryConditions( current, tBC, dt,assignInterfaceValues,assignInterfaceGhostValues );


      // ======================================================================
      // ====================== Boundary Conditions ===========================
      // ======================================================================
            BoundaryConditionParameters extrapParams;

            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
            {
                realMappedGridFunction & fieldNext    =mgp!=NULL ? fields[next]    : cgfields[next][grid];
                realMappedGridFunction & fieldCurrent =mgp!=NULL ? fields[current] : cgfields[current][grid];


                if( true )
                {
          // ** 050402 wdh: This next line seems to be needed for Hz (annulus.rbc) -- could fix in PEC
                    fieldNext.periodicUpdate(); 
                }
        
                int option=0; // not used.
                assignBoundaryConditions( option, grid, tBC, dt, fieldNext, fieldCurrent,current );

        // if( true )
        // {
        //   printF("CGMX: TESTING: re-apply BC's\n");
                    
        //   fieldNext.periodicUpdate(); 
        //   assignBoundaryConditions( option, grid, tBC, dt, fieldNext, fieldCurrent,current );
        // }
                


        // Extrapolate neighbours of interpolation points for the wider upwind stencil
                if( extrapolateInterpolationNeighbours && cg.numberOfComponentGrids()>1 ) 
                {
                    extrapParams.orderOfExtrapolation=dbase.get<int>("orderOfExtrapolationForInterpolationNeighbours");
                    if( debug & 4 )
                        printf("***advanceStructured: orderOfExtrapolationForInterpolationNeighbours=%i\n",
                                      dbase.get<int>("orderOfExtrapolationForInterpolationNeighbours"));
        
          // MappedGridOperators & mgop = mgp!=NULL ? *op : (*cgop)[grid];
          // fieldNext.setOperators(mgop);
                    fieldNext.applyBoundaryCondition(Ca,BCTypes::extrapolateInterpolationNeighbours,BCTypes::allBoundaries,0.,tBC,
                                                                                      extrapParams,grid);
                }

        // --------------- Update parallel ghost and periodic ----------------
                real timeBC=getCPU();
#ifdef USE_PPP
                if( orderOfAccuracyInSpace>2 )  // this doesn't seem to be needed for 2nd order ?
                {
                    real timea=getCPU();

            	  fieldNext.updateGhostBoundaries(); // this update is needed *wdh* Dec 7, 2019

        	  if( updatePolarizationParallelGhost )
        	  {
          	    realMappedGridFunction & pNext= getDispersionModelMappedGridFunction( grid,next    );
          	    pNext.updateGhostBoundaries();
        	  }
        	  
                    timing(timeForUpdateGhostBoundaries)+=getCPU()-timea;
                }
#endif

                fieldNext.periodicUpdate();
            
                timing(timeForBoundaryConditions)+=getCPU()-timeBC;
        // display(fieldNext,"fieldNext after finishBoundaryConditions","%7.4f ");
        

            }
            
            if( debug & 4 )
            {
                if( mgp!=NULL )
                    display(fields[next],sPrintF("fields[next] after advanceBAMX, stage=%d, t=%8.2e",stage,tBC),debugFile,"%8.2e ");
                else
                    cgfields[next].display(sPrintF("cgfields[next] after advanceBAMX, stage=%d, t=%8.2e",stage,tBC),debugFile,"%8.2e ");
            }
    
            if( debug & 4 )
            {
                fPrintF(debugFile,"\n ***************** advanceStructured Errors AFTER BC's and BEFORE assignInterface t=%9.3e ********\n",tBC);
                fprintf(pDebugFile,"\n ***************** advanceStructured Errors AFTER BC's and BEFORE assignInterface t=%9.3e ********\n",tBC);
                getErrors( next,tBC,dt );
            }

      // ================================================================================
      // ================ MATERIAL INTERFACES : STAGE II - GHOST VALUES =================
      // ================================================================================
            assignInterfaceValues=false;      // do not project values on the interface
            assignInterfaceGhostValues=true;  // assign ghost 
            assignInterfaceBoundaryConditions( current, tBC, dt,assignInterfaceValues,assignInterfaceGhostValues );  

            if( debug & 4 && gridHasMaterialInterfaces )
            {
                fPrintF(debugFile,"\n ***************** advanceStructured Errors after assignInterface t=%9.3e ********\n",tBC);
                fprintf(pDebugFile,"\n ***************** advanceStructured Errors after assignInterface t=%9.3e ********\n",tBC);
                getErrors( next,tBC,dt );
            }
    

      // Here we project the field to satsify Gauss's Law  div(E)=0 or div(E) = stuff
      // It seems important to project two steps in a row
            if( !projectBeforeBCs && projectThisStep )
            {
        
                project( numberOfStepsTaken, next, tBC , dt );

                if( true )
                {
          // re-apply the BC;s
                    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                    {
                        realMappedGridFunction & fieldNext    =mgp!=NULL ? fields[next]    : cgfields[next][grid];
                        realMappedGridFunction & fieldCurrent =mgp!=NULL ? fields[current] : cgfields[current][grid];

                        int option=0;
                        assignBoundaryConditions( option, grid, tBC, dt, fieldNext, fieldCurrent,current );
                    }
                }

                const bool computeMaxNorms = numberOfStepsTaken<10  || (debug & 1);
                if( computeMaxNorms )
                {
                    getMaxDivergence(next,tBC);

                    printF("===>> project: After project and BC's: |div(E)-rho|=%8.2e, |div(E)-rho|/|grad(E)|=%8.2e, step=%i, t=%9.3e\n",
                                  divEMax,divEMax/gradEMax,numberOfStepsTaken,tBC);
        
                }    

            }

        } // end if applyBC
        
    } // end for stage 
  //  ----------------------- END STAGE -----------------------------------  

    if( useNewForcingMethod && numberOfForcingFunctions>0 )
        fCurrent = (fCurrent +1) % numberOfForcingFunctions;  // increment current forcing index


    if( debug & 16 )
    {
        fPrintF(debugFile," ******************* advanceStructured: Errors at end t=%9.3e ********\n",t+dt);
        fprintf(pDebugFile," ******************* advanceStructured: Errors at end t=%9.3e ********\n",t+dt);
        getErrors( next,t+dt,dt );
    }


    checkArrays("advanceStructured:end");
    
}

#undef FN



