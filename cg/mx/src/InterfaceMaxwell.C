// This file automatically generated from InterfaceMaxwell.bC with bpp.
// ======================================================================================
// Class to assign interface conditions for Maxwell's equations
// 
// Dec 2020 - initial version wdh
// ======================================================================================


#include "InterfaceMaxwell.h"

#include "ParallelUtility.h"
#include "ParallelGridUtility.h"
#include "DispersiveMaterialParameters.h"
#include "CheckParallel.h"

#define interfaceMaxwell EXTERN_C_NAME(interfacemaxwell)
#define newInterfaceMaxwell EXTERN_C_NAME(newinterfacemaxwell)
#define interface3dMaxwell EXTERN_C_NAME(interface3dmaxwell)
#define interfaceOpt EXTERN_C_NAME(interfaceopt)
extern "C"
{
void interfaceMaxwell( const int&nd, 
                                              const int&nd1a,const int&nd1b,const int&nd2a,const int&nd2b,const int&nd3a,const int&nd3b,
                   		       const int&gridIndexRange1, real&u1, const int&mask1,const real&rsxy1, const real&xy1, 
                                              const int&boundaryCondition1, 
                   		       const int&md1a,const int&md1b,const int&md2a,const int&md2b,const int&md3a,const int&md3b,
                   		       const int&gridIndexRange2, real&u2, const int&mask2,const real&rsxy2, const real&xy2, 
                                              const int&boundaryCondition2,
                   		       const int&ipar, const real&rpar, 
                                              real&aa2, real&aa4, real&aa8, 
                                              int&ipvt2, int&ipvt4, int&ipvt8,
                                              int&ierr );

void newInterfaceMaxwell( const int&nd, 
                                              const int&nd1a,const int&nd1b,const int&nd2a,const int&nd2b,const int&nd3a,const int&nd3b,
                   		       const int&gridIndexRange1, real&u1, const int&mask1,const real&rsxy1, const real&xy1, 
                                              const int&boundaryCondition1, 
                   		       const int&md1a,const int&md1b,const int&md2a,const int&md2b,const int&md3a,const int&md3b,
                   		       const int&gridIndexRange2, real&u2, const int&mask2,const real&rsxy2, const real&xy2, 
                                              const int&boundaryCondition2,
                   		       const int&ipar, const real&rpar, int&ierr );

void interface3dMaxwell( const int&nd, 
                                              const int&nd1a,const int&nd1b,const int&nd2a,const int&nd2b,const int&nd3a,const int&nd3b,
                   		       const int&gridIndexRange1, 
                                              real&u1, const real&u1n, const real&u1m, real& v1, 
                                              const int&mask1,const real&rsxy1, const real&xy1, 
                                              real&p1, const real&p1n, const real&p1m, 
                                              real&q1, const real&q1n, const real&q1m, 
                                              const int&boundaryCondition1, 
                   		       const int&md1a,const int&md1b,const int&md2a,const int&md2b,const int&md3a,const int&md3b,
                                              const int&gridIndexRange2, 
                                              real&u2, const real&u2n, const real&u2m, real & v2,
                                              const int&mask2,const real&rsxy2, const real&xy2, 
                                              real&p2, const real&p2n, const real&p2m,  
                                              real&q2, const real&q2n, const real&q2m,  
                                              const int&boundaryCondition2,
                   		       const int&ipar, const real&rpar, 
                                              real&aa2, real&aa4, real&aa8, 
                                              int&ipvt2, int&ipvt4, int&ipvt8,
                                              int&ierr );

void interfaceOpt( const int&nd, 
                                              const int&nd1a,const int&nd1b,const int&nd2a,const int&nd2b,const int&nd3a,const int&nd3b,
                   		       const int&gridIndexRange1, 
                                              real&u1, const real&u1n, const real&u1m, real& v1, 
                                              const int&mask1,const real&rsxy1, const real&xy1, 
                                              real&p1, const real&p1n, const real&p1m, 
                                              real&q1, const real&q1n, const real&q1m, 
                                              const int&boundaryCondition1, 
                   		       const int&md1a,const int&md1b,const int&md2a,const int&md2b,const int&md3a,const int&md3b,
                                              const int&gridIndexRange2, 
                                              real&u2, const real&u2n, const real&u2m, real & v2,
                                              const int&mask2,const real&rsxy2, const real&xy2, 
                                              real&p2, const real&p2n, const real&p2m,  
                                              real&q2, const real&q2n, const real&q2m,  
                                              const int&boundaryCondition2,
                   		       const int&ipar, const real&rpar, 
                                              real&aa2, real&aa4, real&aa8, 
                                              int&ipvt2, int&ipvt4, int&ipvt8,
                                              int&ierr );

}


static CheckParallel *pCheckParallel=NULL;


// ==================================================================================================
// ------ This macro next assigns the local arrays that hold the interface values in parallel ------
// ==================================================================================================

// ===========================================================================================
//  Macro:  Parallel copy near the interface
//
//     Copy data from array u1 into u1nCopy (which is disptributed like u2n) 
//     v1Copy(Jv2)=v1(Jv1) 
// ===========================================================================================

// ======================================================================================
// Macro: Check that the mask values agree across the interface 
// ======================================================================================
        


// ===========================================================================================
//  This macro determines the local arrays that hold in the interface values in PARALLEL.
//
//  In parallel we build new local arrays with a copy of the values from the
//  other side. 
//    
// NOTES:
//   - called by assignInterfaceBoundaryConditions
// 
//   - Since the arrays on either side the of interface may be distributed differently
//   we copy values from one side to the other so we can solve the interface equations.
// 
//   - We could copy values to one side, solve, and copy the results back -- instead we copy
//   values to both sides and solve on both sides (and do not copy back)
//               
// ===========================================================================================



// ============================================================================
//   Assign the parameters for the optimized interface routines.
// 
//  OPTION: serial, parallel1, parallel2 
//    serial : solve on both grids simultaneously in serial
//    parallel1 : solve on grid 1 
//    parallel2 : solve on grid 2
// ============================================================================
// =============== end assignOptParameters() ==============


// ===============================================================================================
//   Call the optimized interface routines.
//
//  OPTION: serial, parallel1, parallel2 
//    serial : solve on both grids simultaneously in serial
//    parallel1 : solve on grid 1 
//    parallel2 : solve on grid 2
// ===============================================================================================
// ======================= end  assignOptInterface() ====================================







// ======================================================================================
/// \brief Constructor for class to assign interface conditions for Maxwell's equations
// ======================================================================================
InterfaceMaxwell::InterfaceMaxwell( Maxwell & cgmx ) : mx(cgmx)
{
}



// ======================================================================================
/// \brief Destructor for class to assign interface conditions for Maxwell's equations
// ======================================================================================
InterfaceMaxwell::~InterfaceMaxwell()
{
}



// ======================================================================================
/// \brief Pre-compute the coefficients in the stencil for the interface equations
// ======================================================================================
int InterfaceMaxwell::computeStencilCoefficients(int current, real t, real dt )
{
    printF("\n >>>>>>>>>> InterfaceMaxwell:: Entering computeStencilCoefficients\n");
    
    CompositeGrid *& cgp = mx.cgp;
    int & debug = mx.debug;
    int & orderOfAccuracyInSpace = mx.orderOfAccuracyInSpace;

    std::vector<InterfaceInfo> & interfaceInfo = mx.interfaceInfo;  // holds info and work-space for interfaces

    int & numberOfTimeLevels = mx.numberOfTimeLevels;
    bool & useNewInterfaceRoutines = mx.useNewInterfaceRoutines;
    int & interfaceEquationsOption = mx.interfaceEquationsOption;
    int & numberOfIterationsForInterfaceBC = mx.numberOfIterationsForInterfaceBC;
    int & materialInterfaceOption = mx.materialInterfaceOption;
    
    bool & solveForElectricField = mx.solveForElectricField;
    bool & solveForMagneticField = mx.solveForMagneticField;
    
    real & omegaForInterfaceIteration = mx.omegaForInterfaceIteration;
  // real & atolForInterfaceIterations = mx.atolForInterfaceIterations;

    RealArray & epsGrid = mx.epsGrid;
    RealArray & muGrid = mx.muGrid;
    RealArray & cGrid = mx.cGrid;
    

    int & myid = mx.myid;

    int & ex = mx.ex;
    int & ey = mx.ey;
    int & ez = mx.ez;
    int & hx = mx.hx;
    int & hy = mx.hy;
    int & hz = mx.hz;

    int & ext = mx.ext;
    int & eyt = mx.eyt;
    int & ezt = mx.ezt;
    int & hxt = mx.hxt;
    int & hyt = mx.hyt;
    int & hzt = mx.hzt;

    int & pxc = mx.pxc;

    Maxwell::MethodEnum & method = mx.method;
    Maxwell::ForcingEnum & forcingOption  = mx.forcingOption;
    Maxwell::DispersionModelEnum & dispersionModel = mx.dispersionModel;
    Maxwell::KnownSolutionEnum & knownSolutionOption = mx.knownSolutionOption;
    
    OGFunction *& tz = mx.tz;
    

    FILE *&debugFile = mx.debugFile;
    FILE *& pDebugFile = mx.pDebugFile; 

    realCompositeGridFunction *& cgfields = mx.cgfields;

    int assignInterfaceValues=1;
    int assignInterfaceGhostValues=1;


  // --- COPIED FROM assignInterfaceBoundaryConditions.bC ------


    assert( cgp!=NULL );
    CompositeGrid & cg= *cgp;
    const int numberOfDimensions = cg.numberOfDimensions();
    const int numberOfComponentGrids = cg.numberOfComponentGrids();

    Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
    Index Jv[3], &J1=Jv[0], &J2=Jv[1], &J3=Jv[2];

    const int prev = (current-1+numberOfTimeLevels) % numberOfTimeLevels;
    const int next = (current+1) % numberOfTimeLevels;

  // *** could iterate to solve the interface conditions for 4th order

    bool interpolateThisDomain=false;  // set to true if we need to interpolate again after setting the interface values
    
    bool reduceOrderOfAccuracyForSosup=true;

    CheckParallel & checkParallel = mx.dbase.get<CheckParallel>("checkParallel");  // use to check parallel computations  
    pCheckParallel = &checkParallel;  // set pointer 

    int bcOrderOfAccuracy=orderOfAccuracyInSpace;
    if( reduceOrderOfAccuracyForSosup && method==Maxwell::sosup && orderOfAccuracyInSpace==6 )
    {
    // NOTE: for now apply 4th order BC's for sosup order 6
        bcOrderOfAccuracy=4;
    }

  //  ----------------------------------------------------------------------------
  //  ------------------------- loop over interfaces -----------------------------
  //  ----------------------------------------------------------------------------
    for( int inter=0; inter < interfaceInfo.size(); inter++ )
    {
        InterfaceInfo & interface = interfaceInfo[inter]; 

        const int grid1=interface.grid1, side1=interface.side1, dir1=interface.dir1;
        const int grid2=interface.grid2, side2=interface.side2, dir2=interface.dir2;
        
        MappedGrid & mg1 = cg[grid1];
        const IntegerArray & bc1 = mg1.boundaryCondition();
        const IntegerArray & share1 = mg1.sharedBoundaryFlag();

        MappedGrid & mg2 = cg[grid2];
        const IntegerArray & bc2 = mg2.boundaryCondition();
        const IntegerArray & share2 = mg2.sharedBoundaryFlag();


        printF("--InterfaceMaxwell-- Interface: %s=(grid1,side,dir1)=(%i,%i,%i) "
                      "AND %s=(grid2,side2,dir2)=(%i,%i,%i) share=%i\n",
                      (const char*)mg1.getName(),grid1,side1,dir1,
                      (const char*)mg2.getName(),grid2,side2,dir2,share1(side1,dir1));


        realArray & u1 = cgfields[next][grid1];
        realArray & u2 = cgfields[next][grid2];

        realArray & u1n = cgfields[current][grid1];
        realArray & u2n = cgfields[current][grid2];

        realArray & u1m = cgfields[prev][grid1];
        realArray & u2m = cgfields[prev][grid2];


        const int numParallelGhost = u1.getGhostBoundaryWidth(0);

        IntegerArray bc1Local(2,3), bc2Local(2,3);

        const int internalGhostBC=-9876; // set bc at internal parallel ghost to this value 
          
        ParallelGridUtility::getLocalBoundaryConditions( cgfields[next][grid1],bc1Local,internalGhostBC );
        ParallelGridUtility::getLocalBoundaryConditions( cgfields[next][grid2],bc2Local,internalGhostBC  );
    // ::display(bc1Local,"bc1Local");
    // ::display(bc2Local,"bc2Local");


        const int extra=0; // orderOfAccuracyInSpace/2;
        getBoundaryIndex(mg1.gridIndexRange(),side1,dir1,I1,I2,I3,extra);
        getBoundaryIndex(mg2.gridIndexRange(),side2,dir2,J1,J2,J3,extra);
    
    // check that the number of points in the tangential directions match -- eventually we will fix this
        for( int dir=1; dir<mg1.numberOfDimensions(); dir++ )
        {
            int dir1p = (dir1+dir) % mg1.numberOfDimensions();
            int dir2p = (dir2+dir) % mg2.numberOfDimensions();
            if( Iv[dir1p].getLength()!=Jv[dir2p].getLength() )
            {
      	printF("Cgmx::applyInterfaceBC:ERROR: The number of grid points on the two interfaces do not match\n"
             	       " (grid1,side1,dir1,bc1)=(%i,%i,%i,%i) Iv=[%i,%i][%i,%i][%i,%i]\n"
             	       " (grid2,side2,dir2,bc2)=(%i,%i,%i,%i) Jv=[%i,%i][%i,%i][%i,%i]\n",
             	       grid1,side1,dir1,mg1.boundaryCondition(side1,dir1),
             	       I1.getBase(),I1.getBound(),I2.getBase(),I2.getBound(),I3.getBase(),I3.getBound(),  
             	       grid2,side2,dir2,mg2.boundaryCondition(side2,dir2),
             	       J1.getBase(),J1.getBound(),J2.getBase(),J2.getBound(),J3.getBase(),J3.getBound());
      	printF("grid names are %s and %s\n",(const char*)mg1.getName(),(const char*)mg2.getName());
      	OV_ABORT("error"); // this prints the file and line number and calls Overture::abort. 
            }

      // We need to interpolate the grid function if an interface has interpolation points on it:
            if( bc1(0,dir1p)==0 || bc1(1,dir1p)==0 )
            {
      	interpolateThisDomain=true;
            }
            if( bc2(0,dir2p)==0 || bc2(1,dir2p)==0 )
            {
      	interpolateThisDomain=true;
            }
        } // end for dir 
    

        intArray & mask1 = mg1.mask();
        intArray & mask2 = mg2.mask();
        OV_GET_SERIAL_ARRAY(int,mask1,mask1Local);
        OV_GET_SERIAL_ARRAY(int,mask2,mask2Local);

    // #ifdef USE_PPP
    //  intSerialArray mask1Local; getLocalArrayWithGhostBoundaries(mask1,mask1Local);
    //  intSerialArray mask2Local; getLocalArrayWithGhostBoundaries(mask2,mask2Local);
    // #else
    //  intSerialArray & mask1Local = mask1;
    //  intSerialArray & mask2Local = mask2;
    // #endif

        bool isRectangular1= mg1.isRectangular();
        real dx1[3]={0.,0.,0.}; //
        if( isRectangular1 )
            mg1.getDeltaX(dx1);

        bool isRectangular2= mg2.isRectangular();
        real dx2[3]={0.,0.,0.}; //
        if( isRectangular2 )
            mg2.getDeltaX(dx2);


        if( orderOfAccuracyInSpace >=4 && useNewInterfaceRoutines )
        { 
      // WARNING: Do not use the fourth-order rectangular grid interface conditions in interface3d.bf 
      //   since these require iterations which have not been implemented

            if( numberOfDimensions==3 && interfaceEquationsOption==0 && ( !isRectangular1 || !isRectangular2) )
            {
      	printF("--MX--assignInterface:WARNING: interfaceEquationsOption=0 is NOT recommended"
                              " for order=4 in 3D (curvilinear)!! Use new routines: interfaceEquationsOption=1.\n");
                if( numberOfIterationsForInterfaceBC !=1 ) 
      	{
        	  printF("--MX--assignInterface: ERROR numberOfIterationsForInterfaceBC=%i should be 1"
             		 " for interfaceEquationsOption=0 (3D, order=4)\n",numberOfIterationsForInterfaceBC);
                    printf("This is a fatal error to make sure that you really wanted to use interfaceEquationsOption=0 \n");
        	  OV_ABORT("error");
      	}
      	
            }

            if( (isRectangular1 && isRectangular2) && t<= 1.5*dt  )
            {
	// printF("--MX-- assignInterfaceBC: INFO - using curvilinear version since iterations are required.\n");
      	printF("--MX-- assignInterfaceBC: INFO - grids are rectangular but using curvilinear interface version ...\n"
                              "      ... for 4th-order since latest version not implemented for rectangular grids.\n");
            }
            
            isRectangular1=false;
            isRectangular2=false;

        }

    // Extra check added *wdh* May 11, 2020
        if( numberOfDimensions==3 && (isRectangular1 || isRectangular2) )
        {
            if( t<= 1.5*dt  )
                printf("CGMX:assignInterfaceBCs: WARNING: 3D rectangular grids not implemented. Use curvilinear version instead\n");

            isRectangular1=false;
            isRectangular2=false;
            
        }
        

        if( !isRectangular1 )
        {
            mg1.update(MappedGrid::THEinverseVertexDerivative);
            mg2.update(MappedGrid::THEinverseVertexDerivative);
        }
              		  
        assert(isRectangular1==isRectangular2);

        const bool centerNeeded = forcingOption==Maxwell::twilightZoneForcing;  // we need the grid points 

    // For dispersion materials, look up arrays of polarization vectors
    // The number of Pv vectors may be different on either side
    // --- Lookup info for the dispersion model ---
        const DispersiveMaterialParameters & dmp1 = mx.getDispersiveMaterialParameters(grid1);
        const int numberOfPolarizationVectors1 = dmp1.numberOfPolarizationVectors; 

        const DispersiveMaterialParameters & dmp2 = mx.getDispersiveMaterialParameters(grid2);
        const int numberOfPolarizationVectors2 = dmp2.numberOfPolarizationVectors; 

        const bool isDispersive = dmp1.isDispersiveMaterial() || dmp2.isDispersiveMaterial();
        if(  isDispersive && dispersionModel==Maxwell::noDispersion )
        {
      // Sanity check 
            printF("CgMx::assignInterfaceBC:ERROR: Some materials are dispersive but dispersionModel==noDispersion\n"
                          "  Something is wrong here\n");
            OV_ABORT("ERROR");
        }
        


    // --- Get pointers to arrays for the dispersive model ----
        real tempValue;
        real *p1ptr = &tempValue;  // set default when not used
        real *p2ptr = &tempValue;  // set default when not used

        real *p1nptr = &tempValue;  // set default when not used
        real *p2nptr = &tempValue;  // set default when not used

        real *p1mptr = &tempValue;  // set default when not used
        real *p2mptr = &tempValue;  // set default when not used

    // real *p1ptr = u1Local.getDataPointer();  // set default when not used
    // real *p2ptr = u2Local.getDataPointer();  // set default when not used

    // real *p1nptr = u1Local.getDataPointer();  // set default when not used
    // real *p2nptr = u2Local.getDataPointer();  // set default when not used

    // real *p1mptr = u1Local.getDataPointer();  // set default when not used
    // real *p2mptr = u2Local.getDataPointer();  // set default when not used

        if( numberOfPolarizationVectors1>0 )
        {
            realMappedGridFunction & p1 = mx.getDispersionModelMappedGridFunction( grid1,next );
            OV_GET_SERIAL_ARRAY(real,p1,p1Local);
            p1ptr=p1Local.getDataPointer();

            realMappedGridFunction & p1n = mx.getDispersionModelMappedGridFunction( grid1,current );
            OV_GET_SERIAL_ARRAY(real,p1n,p1nLocal);
            p1nptr=p1nLocal.getDataPointer();

            realMappedGridFunction & p1m = mx.getDispersionModelMappedGridFunction( grid1,prev );
            OV_GET_SERIAL_ARRAY(real,p1m,p1mLocal);
            p1mptr=p1mLocal.getDataPointer();
        }
        if( numberOfPolarizationVectors2>0 )
        {
            realMappedGridFunction & p2 = mx.getDispersionModelMappedGridFunction( grid2,next );
            OV_GET_SERIAL_ARRAY(real,p2,p2Local);
            p2ptr=p2Local.getDataPointer();

            realMappedGridFunction & p2n = mx.getDispersionModelMappedGridFunction( grid2,current );
            OV_GET_SERIAL_ARRAY(real,p2n,p2nLocal);
            p2nptr=p2nLocal.getDataPointer();

            realMappedGridFunction & p2m = mx.getDispersionModelMappedGridFunction( grid2,prev );
            OV_GET_SERIAL_ARRAY(real,p2m,p2mLocal);
            p2mptr=p2mLocal.getDataPointer();
        }
        realArray & p1  = numberOfPolarizationVectors1>0 ? mx.getDispersionModelMappedGridFunction( grid1,next ) : u1;
        realArray & p1n = numberOfPolarizationVectors1>0 ? mx.getDispersionModelMappedGridFunction( grid1,current ) : u1n;
        
        realArray & p2  = numberOfPolarizationVectors2>0 ? mx.getDispersionModelMappedGridFunction( grid2,next ) : u2;
        realArray & p2n = numberOfPolarizationVectors2>0 ? mx.getDispersionModelMappedGridFunction( grid2,current ) : u2n;
        

    // --- Get pointers to arrays for the nonlinear model ----
        real *q1ptr = &tempValue;   // set default when not used
        real *q2ptr = &tempValue;   // set default when not used

        real *q1nptr = &tempValue;  // set default when not used
        real *q2nptr = &tempValue;  // set default when not used

        real *q1mptr = &tempValue;  // set default when not used
        real *q2mptr = &tempValue;  // set default when not used

        if( dmp1.isNonlinearMaterial() )
        {
            realMappedGridFunction & q1 = mx.getNonlinearModelMappedGridFunction( grid1,next );
            OV_GET_SERIAL_ARRAY(real,q1,q1Local);
            q1ptr=q1Local.getDataPointer();

            realMappedGridFunction & q1n = mx.getNonlinearModelMappedGridFunction( grid1,current );
            OV_GET_SERIAL_ARRAY(real,q1n,q1nLocal);
            q1nptr=q1nLocal.getDataPointer();

            realMappedGridFunction & q1m = mx.getNonlinearModelMappedGridFunction( grid1,prev );
            OV_GET_SERIAL_ARRAY(real,q1m,q1mLocal);
            q1mptr=q1mLocal.getDataPointer();
        }
        if( dmp2.isNonlinearMaterial() )
        {
            realMappedGridFunction & q2 = mx.getNonlinearModelMappedGridFunction( grid2,next );
            OV_GET_SERIAL_ARRAY(real,q2,q2Local);
            q2ptr=q2Local.getDataPointer();

            realMappedGridFunction & q2n = mx.getNonlinearModelMappedGridFunction( grid2,current );
            OV_GET_SERIAL_ARRAY(real,q2n,q2nLocal);
            q2nptr=q2nLocal.getDataPointer();

            realMappedGridFunction & q2m = mx.getNonlinearModelMappedGridFunction( grid2,prev );
            OV_GET_SERIAL_ARRAY(real,q2m,q2mLocal);
            q2mptr=q2mLocal.getDataPointer();
        }


    // ----------- PARALLEL COPY ---------
    //  Determine the local arrays that hold in the interface values.
    //  In parallel we build new local arrays with a copy of the values from the
    //  other side. 
        bool ok1=true, ok2=true;
        #ifdef USE_PPP
      // ------------- PARALLEL VERSION -------------
            realSerialArray u1Local; getLocalArrayWithGhostBoundaries(u1,u1Local);
            realSerialArray u2Local; getLocalArrayWithGhostBoundaries(u2,u2Local);
            realSerialArray u1nLocal; getLocalArrayWithGhostBoundaries(u1n,u1nLocal);
            realSerialArray u2nLocal; getLocalArrayWithGhostBoundaries(u2n,u2nLocal);
            realSerialArray u1mLocal; getLocalArrayWithGhostBoundaries(u1m,u1mLocal);
            realSerialArray u2mLocal; getLocalArrayWithGhostBoundaries(u2m,u2mLocal);
      // Parallel support for dispersion added *wdh* Nov 18, 2020
      // if( numberOfPolarizationVectors1 >0 || numberOfPolarizationVectors2>0 )
      // {
      //   printF("--MX-- INTERFACE: numberOfPolarizationVectors1=%d, numberOfPolarizationVectors2=%d\n",
      //          numberOfPolarizationVectors1,numberOfPolarizationVectors2);
      //   printF(" dmp1.isDispersiveMaterial=%d, dmp2.isDispersiveMaterial=%d\n",
      //          (int)dmp1.isDispersiveMaterial(), (int)dmp2.isDispersiveMaterial());
      //   OV_ABORT("--MX-- INTERFACE - finish me for dispersive and parallel");
      // }
      // *wdh* 081122 -- We need to check that the points are not traversed in the reverse order from one side to the other **** TODO ****
            if( dir1!=dir2 )
            {
                printF("Cgmx:assignInterfaceBC: Error: in parallel we assume that the interface satisfies\n"
                              " dir1==dir2 : you may have to remake the grid to satisfy this\n");
                OV_ABORT("Error");
        //   printF("ERROR in file %s line %d.\n",__FILE__,__LINE__);
        //   Overture::abort("error");
            }
      // stencil half width: (we copy points (-halfWidth,+halfWidth) 
            const int halfWidth= bcOrderOfAccuracy/2;
      // The total maximum extrapolation width is extrapWidth+1, e.g. extrapWidth=2 -> 1 3 3 1 extrapolation 
      // Here is the desired width for extrapolation: (we may not have enough room for this) (add 1 to get the order of extrapolation)
            const int extrapWidth=bcOrderOfAccuracy;   
            int extrapolationWidth1=extrapWidth;  // actual extrapolation width allowed for grid1 (changed below)
            int extrapolationWidth2=extrapWidth;  // actual extrapolation width allowed for grid2 (changed below)
            int width[2];
            const int nd=4;
            Index Iv1[4], Iv2[4];
            getIndex(mg1.dimension(),Iv1[0],Iv1[1],Iv1[2]);
            getIndex(mg2.dimension(),Iv2[0],Iv2[1],Iv2[2]);
            if( interface.pmask1==NULL )
            {
        // ********** Initialization stage : copy geometry data from near the interface *********
                printF("***** getLocalInterfaceArray: parallel copy of mask1, rsxy1, xy1 etc. (DONE ONCE AT START) ****\n");
                width[0]  =halfWidth;     // copy this many ghost pts
                width[1]  =halfWidth;   
                Iv1[dir1]=Range(mg1.gridIndexRange(side1,dir1)-width[0],mg1.gridIndexRange(side1,dir1)+width[1]);
                Iv2[dir2]=Range(mg2.gridIndexRange(side2,dir2)-width[0],mg2.gridIndexRange(side2,dir2)+width[1]);
        // ==== copy the mask ====
                assert( interface.pmask1==NULL && interface.pmask2==NULL );
                interface.pmask1 = new intSerialArray;
                interface.pmask2 = new intSerialArray;
                intSerialArray & mask1i = *interface.pmask1;
                intSerialArray & mask2i = *interface.pmask2;
        // --- copy values from mask2 into an array mask2b that is distributed in the same way as mask1 ---
                intArray mask2b; mask2b.partition(mask1.getPartition());
                mask2b.redim(mask1.dimension(0),mask1.dimension(1),mask1.dimension(2));
                intSerialArray mask2bLocal; getLocalArrayWithGhostBoundaries(mask2b,mask2bLocal);
                mask2bLocal=0; 
                Iv1[3]=0; Iv2[3]=0;
                ParallelUtility::copy(mask2b,Iv1,mask2,Iv2,nd); 
                mask2b.updateGhostBoundaries();  // I think this IS needed
                mask2i.redim(mask2bLocal.dimension(0),mask2bLocal.dimension(1),mask2bLocal.dimension(2));
                mask2i = mask2bLocal;  // save here -- we only really need to save the pts near the interface
                if( debug & 8 )
                {
                    fprintf(debugFile," [grid1,side1,dir1]=[%d,%d,%d][grid2,side2,dir2]=[%d,%d,%d] \n",grid1,side1,dir1,grid2,side2,dir2);
                    displayMask(mask2,"getLocalInterfaceArray: parallel copy: mask2",debugFile);
                    displayMask(mask2i,"getLocalInterfaceArray: parallel copy: mask2i",pDebugFile);
                }
        // ===== check that the mask arrays match on the interface ======
        // Added by *wdh* Nov 26, 2020
                if( debug & 8 || debug & 8 )
                {
                    fprintf(debugFile," [grid1,side1,dir1]=[%d,%d,%d][grid2,side2,dir2]=[%d,%d,%d] \n",grid1,side1,dir1,grid2,side2,dir2);
                    displayMask(mask1,"mask1",debugFile);
                    displayMask(mask2,"mask2",debugFile);
                    displayMask(mask1Local,"mask1Local",pDebugFile);
                    displayMask(mask2Local,"mask2Local",pDebugFile);
                    fflush(debugFile);
                }
        // ===== check that the mask arrays match on the interface ======  Added by *wdh* Nov 26, 2020
                Index Ib1,Ib2,Ib3;
                getBoundaryIndex(mg1.gridIndexRange(),side1,dir1,Ib1,Ib2,Ib3);
                Index Jb1=Ib1, Jb2=Ib2, Jb3=Ib3;
                {
                    int includeGhost=1;
                    OV_GET_SERIAL_ARRAY(int,mask1,mask1Local);
                    OV_GET_SERIAL_ARRAY(int,mask2b,mask2Local);
                    bool ok1 = ParallelUtility::getLocalArrayBounds(mask1,mask1Local,Ib1,Ib2,Ib3,includeGhost);
                    bool ok2 = ParallelUtility::getLocalArrayBounds(mask2b,mask2Local,Jb1,Jb2,Jb3,includeGhost);
                    if( ok1 && ok2 )
                    {
            // int maskDiff = max(abs(mask1Local(Ib1,Ib2,Ib3)-mask2Local(Jb1,Jb2,Jb3)));
            // Check that both masks are discretization points (ignore isNeeded etc.) *wdh* Nov 29, 2020
                        int maskDiff = max( (mask1Local(Ib1,Ib2,Ib3) & MappedGrid::ISdiscretizationPoint) -
                                                                (mask2Local(Jb1,Jb2,Jb3) & MappedGrid::ISdiscretizationPoint) );
                        if( maskDiff > 0 )
                        {
                            printf("MX:assignInterfaceBC:ERROR: myid=%d: The mask arrays mask1 and mask2b do not match on the interface. \n"
                                          "  This is currently required.\n"
                                          "  Try re-generating the grid with more lines in the normal direction, this sometimes fixes this problem.\n",myid);
                            printf("debug=%d\n",debug);
                            if( debug & 2 )
                            {
                                fprintf(pDebugFile,"MX:assignInterfaceBC:ERROR: The mask arrays mask1 and mask2b do not match on the interface. \n"
                                                " maskDiff=%d\n",maskDiff);
                                fprintf(pDebugFile," [grid1,side1,dir1]=[%d,%d,%d][grid2,side2,dir2]=[%d,%d,%d] \n",grid1,side1,dir1,grid2,side2,dir2);
                                ::display(mask1Local(Ib1,Ib2,Ib3),"mask1Local on the boundary",pDebugFile);
                                ::display(mask2Local(Jb1,Jb2,Jb3),"mask2Local on the boundary",pDebugFile);
                                intSerialArray diff(Ib1,Ib2,Ib3);
                                diff = abs(mask1Local(Ib1,Ib2,Ib3)-mask2Local(Jb1,Jb2,Jb3));
                                ::display(diff,"difference |mask1Local-mask2Local| Local on the boundary",pDebugFile);
                                ::displayMask(mask1Local(Ib1,Ib2,Ib3),"mask1Local on the boundary (displayMask)",pDebugFile);
                                ::displayMask(mask2Local(Jb1,Jb2,Jb3),"mask2Local on the boundary (displayMask)",pDebugFile);
                                fflush(pDebugFile);
                                fclose(pDebugFile);
                            }
                            OV_ABORT("ERROR");
                        }
                        else
                        {
                            if( debug & 2 )
                            {
                                fprintf(pDebugFile,"\n >>> assignInterfaceBC: INFO: interface grid points match: mask arrays agree on the interface: \n");
                                fprintf(pDebugFile,"    [grid1,side1,dir1]=[%d,%d,%d][grid2,side2,dir2]=[%d,%d,%d] \n\n",
                                                            grid1,side1,dir1,grid2,side2,dir2);
                            }
                        }
                    }
                }
        // --- copy values from mask1 into an array mask1b that is distributed in the same way as mask2 ---
                intArray mask1b; mask1b.partition(mask2.getPartition());
                mask1b.redim(mask2.dimension(0),mask2.dimension(1),mask2.dimension(2));
                intSerialArray mask1bLocal; getLocalArrayWithGhostBoundaries(mask1b,mask1bLocal);
                mask1bLocal=0; 
                ParallelUtility::copy(mask1b,Iv2,mask1,Iv1,nd);  
                mask1b.updateGhostBoundaries();  // I think this IS needed
                mask1i.redim(mask1bLocal.dimension(0),mask1bLocal.dimension(1),mask1bLocal.dimension(2));
                mask1i = mask1bLocal;
                if( debug & 8 || debug & 8 )
                {
                    displayMask(mask1bLocal,"getLocalInterfaceArray: parallel copy: mask1bLocal",pDebugFile);
                    displayMask(mask1i,"getLocalInterfaceArray: parallel copy: mask1i",pDebugFile);
                }
        // ===== check that the mask arrays match on the interface ======  Added by *wdh* Nov 26, 2020
                getBoundaryIndex(mg2.gridIndexRange(),side2,dir2,Jb1,Jb2,Jb3);
                Ib1=Jb1, Ib2=Jb2, Ib3=Jb3;
                {
                    int includeGhost=1;
                    OV_GET_SERIAL_ARRAY(int,mask1b,mask1Local);
                    OV_GET_SERIAL_ARRAY(int,mask2,mask2Local);
                    bool ok1 = ParallelUtility::getLocalArrayBounds(mask1b,mask1Local,Ib1,Ib2,Ib3,includeGhost);
                    bool ok2 = ParallelUtility::getLocalArrayBounds(mask2,mask2Local,Jb1,Jb2,Jb3,includeGhost);
                    if( ok1 && ok2 )
                    {
            // int maskDiff = max(abs(mask1Local(Ib1,Ib2,Ib3)-mask2Local(Jb1,Jb2,Jb3)));
            // Check that both masks are discretization points (ignore isNeeded etc.) *wdh* Nov 29, 2020
                        int maskDiff = max( (mask1Local(Ib1,Ib2,Ib3) & MappedGrid::ISdiscretizationPoint) -
                                                                (mask2Local(Jb1,Jb2,Jb3) & MappedGrid::ISdiscretizationPoint) );
                        if( maskDiff > 0 )
                        {
                            printf("MX:assignInterfaceBC:ERROR: myid=%d: The mask arrays mask1b and mask2 do not match on the interface. \n"
                                          "  This is currently required.\n"
                                          "  Try re-generating the grid with more lines in the normal direction, this sometimes fixes this problem.\n",myid);
                            printf("debug=%d\n",debug);
                            if( debug & 2 )
                            {
                                fprintf(pDebugFile,"MX:assignInterfaceBC:ERROR: The mask arrays mask1b and mask2 do not match on the interface. \n"
                                                " maskDiff=%d\n",maskDiff);
                                fprintf(pDebugFile," [grid1,side1,dir1]=[%d,%d,%d][grid2,side2,dir2]=[%d,%d,%d] \n",grid1,side1,dir1,grid2,side2,dir2);
                                ::display(mask1Local(Ib1,Ib2,Ib3),"mask1Local on the boundary",pDebugFile);
                                ::display(mask2Local(Jb1,Jb2,Jb3),"mask2Local on the boundary",pDebugFile);
                                intSerialArray diff(Ib1,Ib2,Ib3);
                                diff = abs(mask1Local(Ib1,Ib2,Ib3)-mask2Local(Jb1,Jb2,Jb3));
                                ::display(diff,"difference |mask1Local-mask2Local| Local on the boundary",pDebugFile);
                                ::displayMask(mask1Local(Ib1,Ib2,Ib3),"mask1Local on the boundary (displayMask)",pDebugFile);
                                ::displayMask(mask2Local(Jb1,Jb2,Jb3),"mask2Local on the boundary (displayMask)",pDebugFile);
                                fflush(pDebugFile);
                                fclose(pDebugFile);
                            }
                            OV_ABORT("ERROR");
                        }
                        else
                        {
                            if( debug & 2 )
                            {
                                fprintf(pDebugFile,"\n >>> assignInterfaceBC: INFO: interface grid points match: mask arrays agree on the interface: \n");
                                fprintf(pDebugFile,"    [grid1,side1,dir1]=[%d,%d,%d][grid2,side2,dir2]=[%d,%d,%d] \n\n",
                                                            grid1,side1,dir1,grid2,side2,dir2);
                            }
                        }
                    }
                }
        // === copy rsxy ===
                if( !isRectangular1 )
                {
                    assert( isRectangular1==isRectangular2 );
                    realArray & rsxy1 = mg1.inverseVertexDerivative();
                    realArray & rsxy2 = mg2.inverseVertexDerivative();
                    realSerialArray rsxy1Local; getLocalArrayWithGhostBoundaries(rsxy1,rsxy1Local);
                    realSerialArray rsxy2Local; getLocalArrayWithGhostBoundaries(rsxy2,rsxy2Local);
                    assert( interface.prsxy1==NULL && interface.prsxy2==NULL );
                    interface.prsxy1 = new realSerialArray;
                    interface.prsxy2 = new realSerialArray;
                    realSerialArray & rsxy1i = *interface.prsxy1;
                    realSerialArray & rsxy2i = *interface.prsxy2;
                    realArray rsxy2b; rsxy2b.partition(rsxy1.getPartition());
                    rsxy2b.redim(rsxy1.dimension(0),rsxy1.dimension(1),rsxy1.dimension(2),rsxy1.dimension(3));
                    realSerialArray rsxy2bLocal; getLocalArrayWithGhostBoundaries(rsxy2b,rsxy2bLocal);
                    rsxy2bLocal=0.; 
                    Iv1[3]=rsxy1.dimension(3); Iv2[3]=Iv1[3];
                    ParallelUtility::copy(rsxy2b,Iv1,rsxy2,Iv2,nd); 
                    rsxy2b.updateGhostBoundaries();  // I think this IS needed
                    rsxy2i.redim(rsxy2bLocal.dimension(0),rsxy2bLocal.dimension(1),rsxy2bLocal.dimension(2),rsxy2bLocal.dimension(3));
                    rsxy2i = rsxy2bLocal;  // save here -- we only really need to save the pts near the interface
                    if( debug & 8 )
                    {
                        display(rsxy2Local,"getLocalInterfaceArray: parallel copy: rsxy2Local",pDebugFile,"%5.2f");
                        display(rsxy2i,"getLocalInterfaceArray: parallel copy: rsxy2i",pDebugFile,"%5.2f");
                    }
                    realArray rsxy1b; rsxy1b.partition(rsxy2.getPartition());
                    rsxy1b.redim(rsxy2.dimension(0),rsxy2.dimension(1),rsxy2.dimension(2),rsxy2.dimension(3));
                    realSerialArray rsxy1bLocal; getLocalArrayWithGhostBoundaries(rsxy1b,rsxy1bLocal);
                    rsxy1bLocal=0.; 
                    ParallelUtility::copy(rsxy1b,Iv2,rsxy1,Iv1,nd);  
                    rsxy1b.updateGhostBoundaries();  // I think this IS needed
                    rsxy1i.redim(rsxy1bLocal.dimension(0),rsxy1bLocal.dimension(1),rsxy1bLocal.dimension(2),rsxy1bLocal.dimension(3));
                    rsxy1i = rsxy1bLocal;
                }
        // === copy xy ===
                if( centerNeeded )
                {
                    realArray & xy1 = mg1.center();
                    realArray & xy2 = mg2.center();
                    realSerialArray xy1Local; getLocalArrayWithGhostBoundaries(xy1,xy1Local);
                    realSerialArray xy2Local; getLocalArrayWithGhostBoundaries(xy2,xy2Local);
                    assert( interface.pxy1==NULL && interface.pxy2==NULL );
                    interface.pxy1 = new realSerialArray;
                    interface.pxy2 = new realSerialArray;
                    realSerialArray & xy1i = *interface.pxy1;
                    realSerialArray & xy2i = *interface.pxy2;
                    realArray xy2b; xy2b.partition(xy1.getPartition());
                    xy2b.redim(xy1.dimension(0),xy1.dimension(1),xy1.dimension(2),xy1.dimension(3));
                    realSerialArray xy2bLocal; getLocalArrayWithGhostBoundaries(xy2b,xy2bLocal);
                    xy2bLocal=0.; 
                    Iv1[3]=xy1.dimension(3); Iv2[3]=Iv1[3];
                    ParallelUtility::copy(xy2b,Iv1,xy2,Iv2,nd); 
                    xy2b.updateGhostBoundaries();  // I think this IS needed
                    xy2i.redim(xy2bLocal.dimension(0),xy2bLocal.dimension(1),xy2bLocal.dimension(2),xy2bLocal.dimension(3));
                    xy2i = xy2bLocal;  // save here -- we only really need to save the pts near the interface
                    if( debug & 8 )
                    {
                        display(xy2Local,"getLocalInterfaceArray: parallel copy: xy2Local",pDebugFile,"%6.2f");
                        display(xy2i,"getLocalInterfaceArray: parallel copy: xy2i",pDebugFile,"%6.2f");
                    }
                    realArray xy1b; xy1b.partition(xy2.getPartition());
                    xy1b.redim(xy2.dimension(0),xy2.dimension(1),xy2.dimension(2),xy2.dimension(3));
                    realSerialArray xy1bLocal; getLocalArrayWithGhostBoundaries(xy1b,xy1bLocal);
                    xy1bLocal=0.; 
                    ParallelUtility::copy(xy1b,Iv2,xy1,Iv1,nd);  
                    xy1b.updateGhostBoundaries();  // I think this IS needed
                    xy1i.redim(xy1bLocal.dimension(0),xy1bLocal.dimension(1),xy1bLocal.dimension(2),xy1bLocal.dimension(3));
                    xy1i = xy1bLocal;
                    if( debug & 8 )
                    {
                        display(xy1Local,"getLocalInterfaceArray: parallel copy: xy1Local",pDebugFile,"%6.2f");
                        display(xy1i,"getLocalInterfaceArray: parallel copy: xy1i",pDebugFile,"%6.2f");
                    }
                }
                printF("***** getLocalInterfaceArray: FINISHED parallel copy of mask1, rsxy1, xy1 etc. ****\n");
            }
      // We are copying values from (side2,dir2) of u2 into (side1,dir1) of u2Copy (an array distributed as u1)
      //
      //            u2 side2=0                        side2=1
      //         X--X--X--X--X--X--X--X--X-- ...  X--X--X--X--X--X
      //             w[0] 0  1    w[1]           w[0]   N w[1]
      // 
      //            u1 side1=0                        side1=1
      //         X--X--X--X--X--X--X--X--X-- ...   --X--X--X--X--X
      //                  0  1    hw                    N 
            Iv1[3]=u1.dimension(3); Iv2[3]=u2.dimension(3);
            width[side2]  =halfWidth;     // copy this many ghost pts
            width[1-side2]=extrapWidth;   //  copy extra interior values 
      // we can only copy as many values as are available in the u1 distribution: 
      // NOTE: parallel ghost points ARE added to the ends of a distributed array so we should have both
      // the normal ghost points and the parallel ghostpoints
            if( side1==0 )
                width[0] = min(width[0], mg1.gridIndexRange(side1,dir1)-(u1.getBase(dir1)-u1.getGhostBoundaryWidth(dir1)));
            else
                width[1] = min(width[1], (u1.getBound(dir1)+u1.getGhostBoundaryWidth(dir1))-mg1.gridIndexRange(side1,dir1));
            Iv1[dir1]=Range(mg1.gridIndexRange(side1,dir1)-width[0],mg1.gridIndexRange(side1,dir1)+width[1]);
            Iv2[dir2]=Range(mg2.gridIndexRange(side2,dir2)-width[0],mg2.gridIndexRange(side2,dir2)+width[1]);
            extrapolationWidth1=width[1-side2];  // here is the actual maximum extrapolation width allowed. 
            if( debug & 8 )
            {
                fprintf(pDebugFile,
                                "interfaceBC:copy u2 on u1-distribution: side1=%i side2=%i extrapWidth=%i \n"
                                "                 u1Local.getBase=%i u1Local.getBound=%i,\n"
                                "                 Iv1=[%i,%i][%i,%i][%i,%i][%i,%i]  Iv2=[%i,%i][%i,%i][%i,%i][%i,%i]\n",
                                side1,side2,
                                extrapolationWidth1,
                                u1Local.getBase(dir1),u1Local.getBound(dir1),
                                Iv1[0].getBase(),Iv1[0].getBound(),Iv1[1].getBase(),Iv1[1].getBound(),
                                Iv1[2].getBase(),Iv1[2].getBound(),Iv1[3].getBase(),Iv1[3].getBound(),
                                Iv2[0].getBase(),Iv2[0].getBound(),Iv2[1].getBase(),Iv2[1].getBound(),
                                Iv2[2].getBase(),Iv2[2].getBound(),Iv2[3].getBase(),Iv2[3].getBound());
            }
      // --- copy values from u2 into an array u2Copy that is distributed in the same way as u1 ---
            realArray u2Copy;
            realSerialArray u2CopyLocal;
                u2Copy.partition(u1.getPartition());
                u2Copy.redim(u1.dimension(0),u1.dimension(1),u1.dimension(2),u2.dimension(3));
                getLocalArrayWithGhostBoundaries(u2Copy,u2CopyLocal);
                u2CopyLocal=0.; 
        // -- components to copy : 
                Iv2[3]=u2.dimension(3);
                Iv1[3]=Iv2[3];
                ParallelUtility::copy(u2Copy,Iv1,u2,Iv2,nd);  // u2Copy(Iv1)=u2(Iv2)
                u2Copy.updateGhostBoundaries(); // *********** these are currently needed ********************   **FIX ME**
      // realArray u2Copy; u2Copy.partition(u1.getPartition());
      // // the next line causes a bug in u2Copy.updateGhostBoundaries(); below -- doesn't like non-zero base
      // // u2Copy.redim(u1.dimension(0),u1.dimension(1),u1.dimension(2),Range(tc2,tc2)); // note last arg
      // u2Copy.redim(u1.dimension(0),u1.dimension(1),u1.dimension(2),u2.dimension(3)); // note last arg
      // realSerialArray u2CopyLocal; getLocalArrayWithGhostBoundaries(u2Copy,u2CopyLocal);
      // u2CopyLocal=0.; 
      // ParallelUtility::copy(u2Copy,Iv1,u2,Iv2,nd);  // u2Copy(Iv1)=u2(Iv2)
      // u2Copy.updateGhostBoundaries(); // *********** these are currently needed ********************
      // // u2Copy(Iv1[0],Iv1[1],Iv1[2],Iv1[3])=u2(Iv2[0],Iv2[1],Iv2[2],Iv2[3]);
      // *** FIX ME *** avoid creating all these arrays (?) 
            realArray u2nCopy; realSerialArray u2nCopyLocal;  
            realArray p2Copy;  realSerialArray p2CopyLocal;  
            realArray p2nCopy; realSerialArray p2nCopyLocal;  
            if( isDispersive )
            {
                    u2nCopy.partition(u1n.getPartition());
                    u2nCopy.redim(u1n.dimension(0),u1n.dimension(1),u1n.dimension(2),u2n.dimension(3));
                    getLocalArrayWithGhostBoundaries(u2nCopy,u2nCopyLocal);
                    u2nCopyLocal=0.; 
          // -- components to copy : 
                    Iv2[3]=u2n.dimension(3);
                    Iv1[3]=Iv2[3];
                    ParallelUtility::copy(u2nCopy,Iv1,u2n,Iv2,nd);  // u2nCopy(Iv1)=u2n(Iv2)
                    u2nCopy.updateGhostBoundaries(); // *********** these are currently needed ********************   **FIX ME**
            }
            if( isDispersive && numberOfPolarizationVectors2>0 )
            {
        // --- dispersive material ---
        //   The interface conditions rely on the solutions at t+dt and t
        // Copy:
        //   u1n, u2n : solution at current time (we are assign interface values at the "next" time
        //   p1,p1n, p2,p2n
                    p2Copy.partition(p1.getPartition());
                    p2Copy.redim(p1.dimension(0),p1.dimension(1),p1.dimension(2),p2.dimension(3));
                    getLocalArrayWithGhostBoundaries(p2Copy,p2CopyLocal);
                    p2CopyLocal=0.; 
          // -- components to copy : 
                    Iv2[3]=p2.dimension(3);
                    Iv1[3]=Iv2[3];
                    ParallelUtility::copy(p2Copy,Iv1,p2,Iv2,nd);  // p2Copy(Iv1)=p2(Iv2)
                    p2Copy.updateGhostBoundaries(); // *********** these are currently needed ********************   **FIX ME**
                    p2nCopy.partition(p1n.getPartition());
                    p2nCopy.redim(p1n.dimension(0),p1n.dimension(1),p1n.dimension(2),p2n.dimension(3));
                    getLocalArrayWithGhostBoundaries(p2nCopy,p2nCopyLocal);
                    p2nCopyLocal=0.; 
          // -- components to copy : 
                    Iv2[3]=p2n.dimension(3);
                    Iv1[3]=Iv2[3];
                    ParallelUtility::copy(p2nCopy,Iv1,p2n,Iv2,nd);  // p2nCopy(Iv1)=p2n(Iv2)
                    p2nCopy.updateGhostBoundaries(); // *********** these are currently needed ********************   **FIX ME**
            }
      // ----
            width[side1]  =halfWidth;     // copy this many ghost pts
            width[1-side1]=extrapWidth;   //  copy extra interior values 
            if( side2==0 )
                width[0] = min(width[0], mg2.gridIndexRange(side2,dir2)-(u2.getBase(dir2)-u2.getGhostBoundaryWidth(dir2))); 
            else
                width[1] = min(width[1], (u2.getBound(dir2)+u2.getGhostBoundaryWidth(dir2))-mg2.gridIndexRange(side2,dir2));
            Iv1[dir1]=Range(mg1.gridIndexRange(side1,dir1)-width[0],mg1.gridIndexRange(side1,dir1)+width[1]);
            Iv2[dir2]=Range(mg2.gridIndexRange(side2,dir2)-width[0],mg2.gridIndexRange(side2,dir2)+width[1]);
            extrapolationWidth2=width[1-side1];  // here is the actual maximum extrapolation width allowed. 
            if( debug & 8 )
            {
                fprintf(pDebugFile,
                                "interfaceBC:copy u1 on u2-distribution: extrapWidth=%i\n"
                                "                 Iv1=[%i,%i][%i,%i][%i,%i][%i,%i]  Iv2=[%i,%i][%i,%i][%i,%i][%i,%i]\n",extrapolationWidth2,
                                Iv1[0].getBase(),Iv1[0].getBound(),Iv1[1].getBase(),Iv1[1].getBound(),
                                Iv1[2].getBase(),Iv1[2].getBound(),Iv1[3].getBase(),Iv1[3].getBound(),
                                Iv2[0].getBase(),Iv2[0].getBound(),Iv2[1].getBase(),Iv2[1].getBound(),
                                Iv2[2].getBase(),Iv2[2].getBound(),Iv2[3].getBase(),Iv2[3].getBound());
            }
            realArray u1Copy; realSerialArray u1CopyLocal; 
                u1Copy.partition(u2.getPartition());
                u1Copy.redim(u2.dimension(0),u2.dimension(1),u2.dimension(2),u1.dimension(3));
                getLocalArrayWithGhostBoundaries(u1Copy,u1CopyLocal);
                u1CopyLocal=0.; 
        // -- components to copy : 
                Iv1[3]=u1.dimension(3);
                Iv2[3]=Iv1[3];
                ParallelUtility::copy(u1Copy,Iv2,u1,Iv1,nd);  // u1Copy(Iv2)=u1(Iv1)
                u1Copy.updateGhostBoundaries(); // *********** these are currently needed ********************   **FIX ME**
      // // copy values from u1 into an array u1Copy that is distributed in the same way as u2
      // realArray u1Copy; u1Copy.partition(u2.getPartition());
      // // u1Copy.redim(u2.dimension(0),u2.dimension(1),u2.dimension(2),Range(tc1,tc1));
      // u1Copy.redim(u2.dimension(0),u2.dimension(1),u2.dimension(2),u1.dimension(3));
      // realSerialArray u1CopyLocal; getLocalArrayWithGhostBoundaries(u1Copy,u1CopyLocal);
      // u1CopyLocal=0.; 
      // ParallelUtility::copy(u1Copy,Iv2,u1,Iv1,nd);  // u1Copy(Iv2)=u1(Iv1)
      // u1Copy.updateGhostBoundaries(); // *********** these are currently needed ********************   **FIX ME**
      // // u1Copy(Iv2[0],Iv2[1],Iv2[2],Iv2[3])=u1(Iv1[0],Iv1[1],Iv1[2],Iv1[3]);
      // *** FIX ME *** avoid creating all these arrays (?) 
            realArray u1nCopy; realSerialArray u1nCopyLocal;  
            realArray p1Copy;  realSerialArray p1CopyLocal;  
            realArray p1nCopy; realSerialArray p1nCopyLocal;  
            if( isDispersive )
            {
                    u1nCopy.partition(u2n.getPartition());
                    u1nCopy.redim(u2n.dimension(0),u2n.dimension(1),u2n.dimension(2),u1n.dimension(3));
                    getLocalArrayWithGhostBoundaries(u1nCopy,u1nCopyLocal);
                    u1nCopyLocal=0.; 
          // -- components to copy : 
                    Iv1[3]=u1n.dimension(3);
                    Iv2[3]=Iv1[3];
                    ParallelUtility::copy(u1nCopy,Iv2,u1n,Iv1,nd);  // u1nCopy(Iv2)=u1n(Iv1)
                    u1nCopy.updateGhostBoundaries(); // *********** these are currently needed ********************   **FIX ME**
            }
            if( isDispersive && numberOfPolarizationVectors1>0 )
            {
        // --- dispersive material ---
        //   The interface conditions rely on the solutions at t+dt and t
        // Copy:
        //   u1n, u2n : solution at current time (we are assign interface values at the "next" time
        //   p1,p1n, p2,p2n
                    p1Copy.partition(p2.getPartition());
                    p1Copy.redim(p2.dimension(0),p2.dimension(1),p2.dimension(2),p1.dimension(3));
                    getLocalArrayWithGhostBoundaries(p1Copy,p1CopyLocal);
                    p1CopyLocal=0.; 
          // -- components to copy : 
                    Iv1[3]=p1.dimension(3);
                    Iv2[3]=Iv1[3];
                    ParallelUtility::copy(p1Copy,Iv2,p1,Iv1,nd);  // p1Copy(Iv2)=p1(Iv1)
                    p1Copy.updateGhostBoundaries(); // *********** these are currently needed ********************   **FIX ME**
                    p1nCopy.partition(p2n.getPartition());
                    p1nCopy.redim(p2n.dimension(0),p2n.dimension(1),p2n.dimension(2),p1n.dimension(3));
                    getLocalArrayWithGhostBoundaries(p1nCopy,p1nCopyLocal);
                    p1nCopyLocal=0.; 
          // -- components to copy : 
                    Iv1[3]=p1n.dimension(3);
                    Iv2[3]=Iv1[3];
                    ParallelUtility::copy(p1nCopy,Iv2,p1n,Iv1,nd);  // p1nCopy(Iv2)=p1n(Iv1)
                    p1nCopy.updateGhostBoundaries(); // *********** these are currently needed ********************   **FIX ME**
            }
            int includeGhost=0;  // do NOT include parallel ghost since we can't apply the stencil there
            ok1 = ParallelUtility::getLocalArrayBounds(u1,u1Local,I1,I2,I3,includeGhost);
            ok2 = ParallelUtility::getLocalArrayBounds(u2,u2Local,J1,J2,J3,includeGhost);
      // ::display(u1Local,"interfaceBC: u1Local","%5.2f ");
      // ::display(u1CopyLocal,"interfaceBC:u1CopyLocal after copy","%5.2f ");
      // ::display(u2Local,"interfaceBC: u2Local","%5.2f ");
      // ::display(u2CopyLocal,"interfaceBC:u2CopyLocal after copy","%5.2f ");
        #else
      // ------------- SERIAL VERSION -------------
            realSerialArray & u1Local = u1;
            realSerialArray & u2Local = u2;
            realSerialArray & u1nLocal = u1n;
            realSerialArray & u2nLocal = u2n;
            realSerialArray & u1mLocal = u1m;
            realSerialArray & u2mLocal = u2m;
        #endif

        if( debug & 4 ) 
        {
            ::display(u1Local,sPrintF("u1Local before assignOptInterface grid1=%i, t=%8.2e",grid1,t),pDebugFile,"%9.2e ");
            ::display(u2Local,sPrintF("u2Local before assignOptInterface grid2=%i, t=%8.2e",grid2,t),pDebugFile,"%9.2e ");
        }
        if( debug & 4 )
        {
            checkParallel.checkDiff(u1Local,sPrintF("u1Local before assignOptInterface grid1=%i, t=%8.2e",grid1,t),pDebugFile);
            checkParallel.checkDiff(u2Local,sPrintF("u2Local before assignOptInterface grid2=%i, t=%8.2e",grid2,t),pDebugFile);
        }


        int n1a=I1.getBase(),n1b=I1.getBound(),
            n2a=I2.getBase(),n2b=I2.getBound(),
            n3a=I3.getBase(),n3b=I3.getBound();

        int m1a=J1.getBase(),m1b=J1.getBound(),
            m2a=J2.getBase(),m2b=J2.getBound(),
            m3a=J3.getBase(),m3b=J3.getBound();


        real temp;
        real *ptemp=&temp;

        int gridType = isRectangular1 ? 0 : 1;
        int orderOfExtrapolation=orderOfAccuracyInSpace+1;  // not used
        int useForcing = forcingOption==Maxwell::twilightZoneForcing;
        int useWhereMask=true;
        int parallel=0;
        #ifdef USE_PPP
            parallel=1;
        #endif

        real *u1p=u1Local.getDataPointer();
        real *u1np=u1nLocal.getDataPointer();
        real *u1mp=u1mLocal.getDataPointer();
        real *prsxy1=isRectangular1 ? ptemp : mg1.inverseVertexDerivative().getLocalArray().getDataPointer();
        real *pxy1= centerNeeded ? mg1.center().getLocalArray().getDataPointer() : ptemp; 
        int *mask1p=mask1Local.getDataPointer();

        real *u2p=u2Local.getDataPointer();
        real *u2np=u2nLocal.getDataPointer();
        real *u2mp=u2mLocal.getDataPointer();
        real *prsxy2=isRectangular2 ? ptemp : mg2.inverseVertexDerivative().getLocalArray().getDataPointer();
        real *pxy2= centerNeeded ? mg2.center().getLocalArray().getDataPointer() : ptemp; 
        int *mask2p=mask2Local.getDataPointer();

        #ifdef USE_PPP
     // pointers to copies of the interface geometry and mask data 
          int *pmask1b = interface.pmask1->getDataPointer();
          int *pmask2b = interface.pmask2->getDataPointer();
        
          real *prsxy1b = isRectangular1 ? ptemp : interface.prsxy1->getDataPointer();
          real *prsxy2b = isRectangular2 ? ptemp : interface.prsxy2->getDataPointer();
          
          real *pxy1b = !centerNeeded ? ptemp : interface.pxy1->getDataPointer();
          real *pxy2b = !centerNeeded ? ptemp : interface.pxy2->getDataPointer();
        #endif

     // We need some temp space for jacobi updates
        const int & useJacobiInterfaceUpdate = mx.dbase.get<int>("useJacobiInterfaceUpdate");

        bool useJacobiUpdate=false;
        if( orderOfAccuracyInSpace==4 )
        {
            useJacobiUpdate = useJacobiInterfaceUpdate; 
        }
        RealArray v1,v2, v1b, v2b;
        real *v1p =u1p, *v2p =u2p; // default when not used
        real *v1bp=u1p, *v2bp=u2p; // default when not used
        if( useJacobiUpdate )
        {
      // --- allocate temp-space for Jacobi update -----

      // Do this for now -- we can optimize ----

            v1.redim(u1Local);         v2.redim(u2Local);
            v1p=v1.getDataPointer();   v2p=v2.getDataPointer();
      // Fix for parallel -- use u1Copy, u2Copy
            #ifdef USE_PPP
                v1b.redim(u1CopyLocal);         v2b.redim(u2CopyLocal);
                v1bp=v1b.getDataPointer();      v2bp=v2b.getDataPointer();
            #endif
        }



        int ierr=0;

        if( bcOrderOfAccuracy<6 )
        { 
      // ------------------------------------
      // ----- order of accuracy <= 4 -------
      // ------------------------------------

      // macro: 

            #ifdef USE_PPP
                assert( dir1==dir2 );
        // In parallel we solve the equations in serial on both sides of the interface
                if( ok1 )
      	{
                    orderOfExtrapolation=extrapolationWidth1+1;
          // Each grid may or may not have dispersion model: 
                    const Maxwell::DispersionModelEnum dispersionModel1 = dmp1.numberOfPolarizationVectors>0 ? dispersionModel : Maxwell::noDispersion;
                    const Maxwell::DispersionModelEnum dispersionModel2 = dmp2.numberOfPolarizationVectors>0 ? dispersionModel : Maxwell::noDispersion;
                    const int nonlinearModel1 = dmp1.getNonlinearModel();
                    const int nonlinearModel2 = dmp2.getNonlinearModel();
                    int ipar[]={ //
                        side1, dir1, grid1,         // keep side1,dir1 since we don't reverse the points.
                            n1a,n1b,n2a,n2b,n3a,n3b,
                        side2, dir2, grid2,         //  keep side2,dir2 since we don't reverse the points.
                            n1a,n1b,n2a,n2b,n3a,n3b,  // use grid1 dimensions for grid2 when we solve on grid1
                        gridType,            
                        bcOrderOfAccuracy,
                        orderOfExtrapolation,
                        useForcing,          
                        ex,                  
                        ey,                  
                        ez,                  
                        hx ,                 
                        hy,                  
                        hz,                  
                        (int)solveForElectricField,          
                        (int)solveForMagneticField,          
                        useWhereMask,       
                        debug,
                        numberOfIterationsForInterfaceBC,
                        materialInterfaceOption,
                        (int)interface.initialized,
                        myid,
                        parallel,
                        (int)forcingOption,
                        interfaceEquationsOption,
                        (int)assignInterfaceValues,
                        (int)assignInterfaceGhostValues,
                        mx.dbase.get<int>("setDivergenceAtInterfaces"),
                        mx.dbase.get<int>("useImpedanceInterfaceProjection"),
                        0,   // numberOfInterfaceIterationsUsed : returned value ipar[43]
                        dispersionModel1,   // ipar[44]
                        dispersionModel2,   // ipar[45]
                        pxc,                // ipar[46]
                        knownSolutionOption,
                        useJacobiUpdate,
                        nonlinearModel1,      // ipar[49]
                        nonlinearModel2,      // ipar[50]
                        numParallelGhost,     // ipar[51]
                        internalGhostBC       // ipar[52]
                    };
                const real & rtolForInterfaceIterations = mx.dbase.get<real>("rtolForInterfaceIterations");
                const real & atolForInterfaceIterations = mx.dbase.get<real>("atolForInterfaceIterations");
                real rpar[]={ //
                    dx1[0],
                    dx1[1],
                    dx1[2],
                    mg1.gridSpacing(0),
                    mg1.gridSpacing(1),
                    mg1.gridSpacing(2),
                    dx2[0],
                    dx2[1],
                    dx2[2],
                    mg2.gridSpacing(0),
                    mg2.gridSpacing(1),
                    mg2.gridSpacing(2),
                    t,    
                    (real &)tz,  // twilight zone pointer
                    dt,    
                    epsGrid(grid1),
                    muGrid(grid1),   
                    cGrid(grid1),    
                    epsGrid(grid2),  
                    muGrid(grid2),   
                    cGrid(grid2),
                    omegaForInterfaceIteration,
                    0., // return value averageInterfaceConvergenceRate
                    0.,  // return value maxFinalResidual
                    rtolForInterfaceIterations,
                    atolForInterfaceIterations
                };
        // work space: 
                real *rwk=interface.rwk;
                int *iwk=interface.iwk;
                assert( rwk!=NULL && iwk!=NULL );
                const int ndf = max(interface.ndf1,interface.ndf2); 
        // assign pointers into the work spaces
                int pa2=0,pa4=0,pa8=0, pipvt2=0,pipvt4=0,pipvt8=0;
                if( bcOrderOfAccuracy==2 )
                {
                    if( !( mg1.numberOfDimensions()==3 ) )  // new 3D doesn't use work-space yet
                    {
                        pa2=0; 
                        pa4=pa2 + 2*2*2*ndf;
                        pa8=0;  // not used
                        pipvt2=0;
                        pipvt4=pipvt2 + 2*ndf; 
                        pipvt8=0;
                    }
                }
                else if( bcOrderOfAccuracy==4 )
                {
                    pa2=0; // not used
                    pa4=0;
                    pa8=pa4+4*4*2*ndf;
                    pipvt2=0;
                    pipvt4=0;
                    pipvt8=pipvt4+4*ndf;
                }
                if( mg1.numberOfDimensions()==2 && !useNewInterfaceRoutines )
                {
          // OLD interface routines
                    interfaceMaxwell( mg1.numberOfDimensions(), 
                                  		      u1Local.getBase(0),u1Local.getBound(0),
                                  		      u1Local.getBase(1),u1Local.getBound(1),
                                  		      u1Local.getBase(2),u1Local.getBound(2),
                                  		      mg1.gridIndexRange(0,0), *u1p, *mask1p,*prsxy1, *pxy1, bc1Local(0,0), 
                                  		      u2CopyLocal.getBase(0),u2CopyLocal.getBound(0),
                                  		      u2CopyLocal.getBase(1),u2CopyLocal.getBound(1),
                                  		      u2CopyLocal.getBase(2),u2CopyLocal.getBound(2),
        		      // note: use grid1 mesh data here ASSUMES GRIDS MATCH *FIX ME*
        		      // mg1.gridIndexRange(0,0), *u2CopyLocal.getDataPointer(), *mask1p,*prsxy1, *pxy1, bc1Local(0,0), 
                              // fixed version: but note bc1Local 
                                                            mg1.gridIndexRange(0,0), *u2CopyLocal.getDataPointer(), *pmask2b,*prsxy2b,*pxy2b, bc1Local(0,0),
                                		    ipar[0], rpar[0], 
                                		    rwk[pa2],rwk[pa4],rwk[pa8], iwk[pipvt2],iwk[pipvt4],iwk[pipvt8],
                                		    ierr );
                }
                else
                {
          // new interface routines -- 2D versions are done here too
          // interface3dMaxwell( mg1.numberOfDimensions(), 
                    interfaceOpt( mg1.numberOfDimensions(), 
                                // --- serial case ---
                                    		        u1Local.getBase(0),u1Local.getBound(0),
                                    		        u1Local.getBase(1),u1Local.getBound(1),
                                    		        u1Local.getBase(2),u1Local.getBound(2),
                                                                mg1.gridIndexRange(0,0), 
                                                                *u1p, *u1np, *u1mp, *v1p,
                                                                *mask1p,*prsxy1, *pxy1, 
                                                                *p1ptr,*p1nptr,*p1mptr,  
                                                                *q1ptr,*q1nptr,*q1mptr,   
                                                                bc1Local(0,0), 
                                    		        u2CopyLocal.getBase(0),u2CopyLocal.getBound(0),
                                    		        u2CopyLocal.getBase(1),u2CopyLocal.getBound(1),
                                    		        u2CopyLocal.getBase(2),u2CopyLocal.getBound(2),
        		        // note: use grid1 mesh data here ASSUMES GRIDS MATCH *FIX ME*
        		        // mg1.gridIndexRange(0,0), *u2CopyLocal.getDataPointer(), *mask1p,*prsxy1, *pxy1, bc1Local(0,0),
                                // fixed version: but note bc1Local 
                                                                mg1.gridIndexRange(0,0), 
                                                                *u2CopyLocal.getDataPointer(),
                                // *u2np,
                                                                *u2nCopyLocal.getDataPointer(),
                                                                *u2mp, *v2bp,
                                                                *pmask2b,*prsxy2b,*pxy2b, 
                                // *p2ptr,*p2nptr,
                                                                *p2CopyLocal.getDataPointer(),
                                                                *p2nCopyLocal.getDataPointer(),
                                                                *p2mptr, 
                                                                *q2ptr,*q2nptr,*q2mptr, 
                                                                bc1Local(0,0),
                                  		      ipar[0], rpar[0], 
                                  		      rwk[pa2],rwk[pa4],rwk[pa8], iwk[pipvt2],iwk[pipvt4],iwk[pipvt8],
                                  		      ierr );
                }
        // sosup solves for E and E.t so we need to divide some counts by 2: 
                const real numSolvesPerStep = method==Maxwell::sosup ? 2. : 1.;
                interface.totalInterfaceIterations+=ipar[43]/numSolvesPerStep; // counts total number of interface iterations
                interface.averageInterfaceConvergenceRate+=rpar[22]/numSolvesPerStep;
                interface.maxFinalResidual=max(rpar[23],interface.maxFinalResidual);
                interface.averageFinalResidual+=rpar[23]/numSolvesPerStep;  // keeps sums of residuals 
        	  if( method==Maxwell::sosup )
        	  { // for sosup we assign E.t
            // Each grid may or may not have dispersion model: 
                        const Maxwell::DispersionModelEnum dispersionModel1 = dmp1.numberOfPolarizationVectors>0 ? dispersionModel : Maxwell::noDispersion;
                        const Maxwell::DispersionModelEnum dispersionModel2 = dmp2.numberOfPolarizationVectors>0 ? dispersionModel : Maxwell::noDispersion;
                        const int nonlinearModel1 = dmp1.getNonlinearModel();
                        const int nonlinearModel2 = dmp2.getNonlinearModel();
                        int ipar[]={ //
                            side1, dir1, grid1,         // keep side1,dir1 since we don't reverse the points.
                                n1a,n1b,n2a,n2b,n3a,n3b,
                            side2, dir2, grid2,         //  keep side2,dir2 since we don't reverse the points.
                                n1a,n1b,n2a,n2b,n3a,n3b,  // use grid1 dimensions for grid2 when we solve on grid1
                            gridType,            
                            bcOrderOfAccuracy,
                            orderOfExtrapolation,
                            useForcing,          
                            ext,                  
                            eyt,                  
                            ezt,                  
                            hxt ,                 
                            hyt,                  
                            hzt,                  
                            (int)solveForElectricField,          
                            (int)solveForMagneticField,          
                            useWhereMask,       
                            debug,
                            numberOfIterationsForInterfaceBC,
                            materialInterfaceOption,
                            (int)interface.initialized,
                            myid,
                            parallel,
                            (int)forcingOption,
                            interfaceEquationsOption,
                            (int)assignInterfaceValues,
                            (int)assignInterfaceGhostValues,
                            mx.dbase.get<int>("setDivergenceAtInterfaces"),
                            mx.dbase.get<int>("useImpedanceInterfaceProjection"),
                            0,   // numberOfInterfaceIterationsUsed : returned value ipar[43]
                            dispersionModel1,   // ipar[44]
                            dispersionModel2,   // ipar[45]
                            pxc,                // ipar[46]
                            knownSolutionOption,
                            useJacobiUpdate,
                            nonlinearModel1,      // ipar[49]
                            nonlinearModel2,      // ipar[50]
                            numParallelGhost,     // ipar[51]
                            internalGhostBC       // ipar[52]
                        };
                    const real & rtolForInterfaceIterations = mx.dbase.get<real>("rtolForInterfaceIterations");
                    const real & atolForInterfaceIterations = mx.dbase.get<real>("atolForInterfaceIterations");
                    real rpar[]={ //
                        dx1[0],
                        dx1[1],
                        dx1[2],
                        mg1.gridSpacing(0),
                        mg1.gridSpacing(1),
                        mg1.gridSpacing(2),
                        dx2[0],
                        dx2[1],
                        dx2[2],
                        mg2.gridSpacing(0),
                        mg2.gridSpacing(1),
                        mg2.gridSpacing(2),
                        t,    
                        (real &)tz,  // twilight zone pointer
                        dt,    
                        epsGrid(grid1),
                        muGrid(grid1),   
                        cGrid(grid1),    
                        epsGrid(grid2),  
                        muGrid(grid2),   
                        cGrid(grid2),
                        omegaForInterfaceIteration,
                        0., // return value averageInterfaceConvergenceRate
                        0.,  // return value maxFinalResidual
                        rtolForInterfaceIterations,
                        atolForInterfaceIterations
                    };
          // work space: 
                    real *rwk=interface.rwk;
                    int *iwk=interface.iwk;
                    assert( rwk!=NULL && iwk!=NULL );
                    const int ndf = max(interface.ndf1,interface.ndf2); 
          // assign pointers into the work spaces
                    int pa2=0,pa4=0,pa8=0, pipvt2=0,pipvt4=0,pipvt8=0;
                    if( bcOrderOfAccuracy==2 )
                    {
                        if( !( mg1.numberOfDimensions()==3 ) )  // new 3D doesn't use work-space yet
                        {
                            pa2=0; 
                            pa4=pa2 + 2*2*2*ndf;
                            pa8=0;  // not used
                            pipvt2=0;
                            pipvt4=pipvt2 + 2*ndf; 
                            pipvt8=0;
                        }
                    }
                    else if( bcOrderOfAccuracy==4 )
                    {
                        pa2=0; // not used
                        pa4=0;
                        pa8=pa4+4*4*2*ndf;
                        pipvt2=0;
                        pipvt4=0;
                        pipvt8=pipvt4+4*ndf;
                    }
                    if( mg1.numberOfDimensions()==2 && !useNewInterfaceRoutines )
                    {
            // OLD interface routines
                        interfaceMaxwell( mg1.numberOfDimensions(), 
                                      		      u1Local.getBase(0),u1Local.getBound(0),
                                      		      u1Local.getBase(1),u1Local.getBound(1),
                                      		      u1Local.getBase(2),u1Local.getBound(2),
                                      		      mg1.gridIndexRange(0,0), *u1p, *mask1p,*prsxy1, *pxy1, bc1Local(0,0), 
                                      		      u2CopyLocal.getBase(0),u2CopyLocal.getBound(0),
                                      		      u2CopyLocal.getBase(1),u2CopyLocal.getBound(1),
                                      		      u2CopyLocal.getBase(2),u2CopyLocal.getBound(2),
          		      // note: use grid1 mesh data here ASSUMES GRIDS MATCH *FIX ME*
          		      // mg1.gridIndexRange(0,0), *u2CopyLocal.getDataPointer(), *mask1p,*prsxy1, *pxy1, bc1Local(0,0), 
                                // fixed version: but note bc1Local 
                                                                mg1.gridIndexRange(0,0), *u2CopyLocal.getDataPointer(), *pmask2b,*prsxy2b,*pxy2b, bc1Local(0,0),
                                    		    ipar[0], rpar[0], 
                                    		    rwk[pa2],rwk[pa4],rwk[pa8], iwk[pipvt2],iwk[pipvt4],iwk[pipvt8],
                                    		    ierr );
                    }
                    else
                    {
            // new interface routines -- 2D versions are done here too
            // interface3dMaxwell( mg1.numberOfDimensions(), 
                        interfaceOpt( mg1.numberOfDimensions(), 
                                  // --- serial case ---
                                        		        u1Local.getBase(0),u1Local.getBound(0),
                                        		        u1Local.getBase(1),u1Local.getBound(1),
                                        		        u1Local.getBase(2),u1Local.getBound(2),
                                                                    mg1.gridIndexRange(0,0), 
                                                                    *u1p, *u1np, *u1mp, *v1p,
                                                                    *mask1p,*prsxy1, *pxy1, 
                                                                    *p1ptr,*p1nptr,*p1mptr,  
                                                                    *q1ptr,*q1nptr,*q1mptr,   
                                                                    bc1Local(0,0), 
                                        		        u2CopyLocal.getBase(0),u2CopyLocal.getBound(0),
                                        		        u2CopyLocal.getBase(1),u2CopyLocal.getBound(1),
                                        		        u2CopyLocal.getBase(2),u2CopyLocal.getBound(2),
          		        // note: use grid1 mesh data here ASSUMES GRIDS MATCH *FIX ME*
          		        // mg1.gridIndexRange(0,0), *u2CopyLocal.getDataPointer(), *mask1p,*prsxy1, *pxy1, bc1Local(0,0),
                                  // fixed version: but note bc1Local 
                                                                    mg1.gridIndexRange(0,0), 
                                                                    *u2CopyLocal.getDataPointer(),
                                  // *u2np,
                                                                    *u2nCopyLocal.getDataPointer(),
                                                                    *u2mp, *v2bp,
                                                                    *pmask2b,*prsxy2b,*pxy2b, 
                                  // *p2ptr,*p2nptr,
                                                                    *p2CopyLocal.getDataPointer(),
                                                                    *p2nCopyLocal.getDataPointer(),
                                                                    *p2mptr, 
                                                                    *q2ptr,*q2nptr,*q2mptr, 
                                                                    bc1Local(0,0),
                                      		      ipar[0], rpar[0], 
                                      		      rwk[pa2],rwk[pa4],rwk[pa8], iwk[pipvt2],iwk[pipvt4],iwk[pipvt8],
                                      		      ierr );
                    }
          // sosup solves for E and E.t so we need to divide some counts by 2: 
                    const real numSolvesPerStep = method==Maxwell::sosup ? 2. : 1.;
                    interface.totalInterfaceIterations+=ipar[43]/numSolvesPerStep; // counts total number of interface iterations
                    interface.averageInterfaceConvergenceRate+=rpar[22]/numSolvesPerStep;
                    interface.maxFinalResidual=max(rpar[23],interface.maxFinalResidual);
                    interface.averageFinalResidual+=rpar[23]/numSolvesPerStep;  // keeps sums of residuals 
        	  }
        	  if( debug & 16 )
                        ::display(u1Local,sPrintF("u1Local after assignOptInterface t=%8.2e",t),pDebugFile,"%5.2f ");
      	}
      	if( ok2 )
      	{
                    orderOfExtrapolation=extrapolationWidth2+1;
          // Each grid may or may not have dispersion model: 
                    const Maxwell::DispersionModelEnum dispersionModel1 = dmp1.numberOfPolarizationVectors>0 ? dispersionModel : Maxwell::noDispersion;
                    const Maxwell::DispersionModelEnum dispersionModel2 = dmp2.numberOfPolarizationVectors>0 ? dispersionModel : Maxwell::noDispersion;
                    const int nonlinearModel1 = dmp1.getNonlinearModel();
                    const int nonlinearModel2 = dmp2.getNonlinearModel();
                    int ipar[]={ //
                        side1, dir1, grid1,         // keep side1,dir1 since we don't reverse the points.
                            m1a,m1b,m2a,m2b,m3a,m3b,  // use grid2 dimensions for grid1 when we solve on grid2
                        side2, dir2, grid2,         //  keep side2,dir2 since we don't reverse the points.
                            m1a,m1b,m2a,m2b,m3a,m3b,
                        gridType,            
                        bcOrderOfAccuracy,
                        orderOfExtrapolation,
                        useForcing,          
                        ex,                  
                        ey,                  
                        ez,                  
                        hx ,                 
                        hy,                  
                        hz,                  
                        (int)solveForElectricField,          
                        (int)solveForMagneticField,          
                        useWhereMask,       
                        debug,
                        numberOfIterationsForInterfaceBC,
                        materialInterfaceOption,
                        (int)interface.initialized,
                        myid,
                        parallel,
                        (int)forcingOption,
                        interfaceEquationsOption,
                        (int)assignInterfaceValues,
                        (int)assignInterfaceGhostValues,
                        mx.dbase.get<int>("setDivergenceAtInterfaces"),
                        mx.dbase.get<int>("useImpedanceInterfaceProjection"),
                        0,   // numberOfInterfaceIterationsUsed : returned value ipar[43]
                        dispersionModel1,   // ipar[44]
                        dispersionModel2,   // ipar[45]
                        pxc,                // ipar[46]
                        knownSolutionOption,
                        useJacobiUpdate,
                        nonlinearModel1,      // ipar[49]
                        nonlinearModel2,      // ipar[50]
                        numParallelGhost,     // ipar[51]
                        internalGhostBC       // ipar[52]
                    };
                const real & rtolForInterfaceIterations = mx.dbase.get<real>("rtolForInterfaceIterations");
                const real & atolForInterfaceIterations = mx.dbase.get<real>("atolForInterfaceIterations");
                real rpar[]={ //
                    dx1[0],
                    dx1[1],
                    dx1[2],
                    mg1.gridSpacing(0),
                    mg1.gridSpacing(1),
                    mg1.gridSpacing(2),
                    dx2[0],
                    dx2[1],
                    dx2[2],
                    mg2.gridSpacing(0),
                    mg2.gridSpacing(1),
                    mg2.gridSpacing(2),
                    t,    
                    (real &)tz,  // twilight zone pointer
                    dt,    
                    epsGrid(grid1),
                    muGrid(grid1),   
                    cGrid(grid1),    
                    epsGrid(grid2),  
                    muGrid(grid2),   
                    cGrid(grid2),
                    omegaForInterfaceIteration,
                    0., // return value averageInterfaceConvergenceRate
                    0.,  // return value maxFinalResidual
                    rtolForInterfaceIterations,
                    atolForInterfaceIterations
                };
        // work space: 
                real *rwk=interface.rwk;
                int *iwk=interface.iwk;
                assert( rwk!=NULL && iwk!=NULL );
                const int ndf = max(interface.ndf1,interface.ndf2); 
        // assign pointers into the work spaces
                int pa2=0,pa4=0,pa8=0, pipvt2=0,pipvt4=0,pipvt8=0;
                if( bcOrderOfAccuracy==2 )
                {
                    if( !( mg1.numberOfDimensions()==3 ) )  // new 3D doesn't use work-space yet
                    {
                        pa2=0; 
                        pa4=pa2 + 2*2*2*ndf;
                        pa8=0;  // not used
                        pipvt2=0;
                        pipvt4=pipvt2 + 2*ndf; 
                        pipvt8=0;
                    }
                }
                else if( bcOrderOfAccuracy==4 )
                {
                    pa2=0; // not used
                    pa4=0;
                    pa8=pa4+4*4*2*ndf;
                    pipvt2=0;
                    pipvt4=0;
                    pipvt8=pipvt4+4*ndf;
                }
                if( mg1.numberOfDimensions()==2 && !useNewInterfaceRoutines )
                {
          // OLD interface routines
                    interfaceMaxwell( mg1.numberOfDimensions(), 
                                  		      u1CopyLocal.getBase(0),u1CopyLocal.getBound(0),
                                  		      u1CopyLocal.getBase(1),u1CopyLocal.getBound(1),
                                  		      u1CopyLocal.getBase(2),u1CopyLocal.getBound(2),
        		      // note: use grid2 mesh data here  ASSUMES GRIDS MATCH *FIX ME*
                              // mg2.gridIndexRange(0,0), *u1CopyLocal.getDataPointer(), *mask2p,*prsxy2, *pxy2, bc2Local(0,0),
                              // fixed version: but note bc2Local 
                                                            mg2.gridIndexRange(0,0), *u1CopyLocal.getDataPointer(), *pmask1b,*prsxy1b, *pxy1b, bc2Local(0,0),
                                  		      u2Local.getBase(0),u2Local.getBound(0),
                                  		      u2Local.getBase(1),u2Local.getBound(1),
                                  		      u2Local.getBase(2),u2Local.getBound(2),
                                  		      mg2.gridIndexRange(0,0), *u2p, *mask2p,*prsxy2, *pxy2, bc2Local(0,0), 
                                		    ipar[0], rpar[0], 
                                		    rwk[pa2],rwk[pa4],rwk[pa8], iwk[pipvt2],iwk[pipvt4],iwk[pipvt8],
                                		    ierr );
                }
                else
                {
          // new interface routines -- 2D versions are done here too
          // interface3dMaxwell( mg1.numberOfDimensions(), 
                    interfaceOpt( mg1.numberOfDimensions(), 
                                    		        u1CopyLocal.getBase(0),u1CopyLocal.getBound(0),
                                    		        u1CopyLocal.getBase(1),u1CopyLocal.getBound(1),
                                    		        u1CopyLocal.getBase(2),u1CopyLocal.getBound(2),
        		        // note: use grid2 mesh data here  ASSUMES GRIDS MATCH *FIX ME*
                                // mg2.gridIndexRange(0,0), *u1CopyLocal.getDataPointer(), *mask2p,*prsxy2, *pxy2, bc2Local(0,0),
                                // fixed version: but note bc2Local 
                                                                mg2.gridIndexRange(0,0), 
                                                                *u1CopyLocal.getDataPointer(),   
                                // *u1np,
                                                                *u1nCopyLocal.getDataPointer(),
                                                                *u1mp, *v1bp, 
                                                                *pmask1b,*prsxy1b,*pxy1b, 
                                // *p1ptr,*p1nptr,
                                                                *p1CopyLocal.getDataPointer(),
                                                                *p1nCopyLocal.getDataPointer(),
                                                                *p1mptr,   
                                                                *q1ptr,*q1nptr,*q1mptr,   
                                                                bc2Local(0,0),
                                // --- serial case ---
                                    		        u2Local.getBase(0),u2Local.getBound(0),
                                    		        u2Local.getBase(1),u2Local.getBound(1),
                                    		        u2Local.getBase(2),u2Local.getBound(2),
                                                                mg2.gridIndexRange(0,0), 
                                                                *u2p, *u2np, *u2mp, *v2p,
                                                                *mask2p,*prsxy2, *pxy2, 
                                                                *p2ptr,*p2nptr,*p2mptr, 
                                                                *q2ptr,*q2nptr,*q2mptr, 
                                                                bc2Local(0,0), 
                                  		      ipar[0], rpar[0], 
                                  		      rwk[pa2],rwk[pa4],rwk[pa8], iwk[pipvt2],iwk[pipvt4],iwk[pipvt8],
                                  		      ierr );
                }
        // sosup solves for E and E.t so we need to divide some counts by 2: 
                const real numSolvesPerStep = method==Maxwell::sosup ? 2. : 1.;
                interface.totalInterfaceIterations+=ipar[43]/numSolvesPerStep; // counts total number of interface iterations
                interface.averageInterfaceConvergenceRate+=rpar[22]/numSolvesPerStep;
                interface.maxFinalResidual=max(rpar[23],interface.maxFinalResidual);
                interface.averageFinalResidual+=rpar[23]/numSolvesPerStep;  // keeps sums of residuals 
        	  if( method==Maxwell::sosup )
        	  { // for sosup we assign E.t
            // Each grid may or may not have dispersion model: 
                        const Maxwell::DispersionModelEnum dispersionModel1 = dmp1.numberOfPolarizationVectors>0 ? dispersionModel : Maxwell::noDispersion;
                        const Maxwell::DispersionModelEnum dispersionModel2 = dmp2.numberOfPolarizationVectors>0 ? dispersionModel : Maxwell::noDispersion;
                        const int nonlinearModel1 = dmp1.getNonlinearModel();
                        const int nonlinearModel2 = dmp2.getNonlinearModel();
                        int ipar[]={ //
                            side1, dir1, grid1,         // keep side1,dir1 since we don't reverse the points.
                                m1a,m1b,m2a,m2b,m3a,m3b,  // use grid2 dimensions for grid1 when we solve on grid2
                            side2, dir2, grid2,         //  keep side2,dir2 since we don't reverse the points.
                                m1a,m1b,m2a,m2b,m3a,m3b,
                            gridType,            
                            bcOrderOfAccuracy,
                            orderOfExtrapolation,
                            useForcing,          
                            ext,                  
                            eyt,                  
                            ezt,                  
                            hxt ,                 
                            hyt,                  
                            hzt,                  
                            (int)solveForElectricField,          
                            (int)solveForMagneticField,          
                            useWhereMask,       
                            debug,
                            numberOfIterationsForInterfaceBC,
                            materialInterfaceOption,
                            (int)interface.initialized,
                            myid,
                            parallel,
                            (int)forcingOption,
                            interfaceEquationsOption,
                            (int)assignInterfaceValues,
                            (int)assignInterfaceGhostValues,
                            mx.dbase.get<int>("setDivergenceAtInterfaces"),
                            mx.dbase.get<int>("useImpedanceInterfaceProjection"),
                            0,   // numberOfInterfaceIterationsUsed : returned value ipar[43]
                            dispersionModel1,   // ipar[44]
                            dispersionModel2,   // ipar[45]
                            pxc,                // ipar[46]
                            knownSolutionOption,
                            useJacobiUpdate,
                            nonlinearModel1,      // ipar[49]
                            nonlinearModel2,      // ipar[50]
                            numParallelGhost,     // ipar[51]
                            internalGhostBC       // ipar[52]
                        };
                    const real & rtolForInterfaceIterations = mx.dbase.get<real>("rtolForInterfaceIterations");
                    const real & atolForInterfaceIterations = mx.dbase.get<real>("atolForInterfaceIterations");
                    real rpar[]={ //
                        dx1[0],
                        dx1[1],
                        dx1[2],
                        mg1.gridSpacing(0),
                        mg1.gridSpacing(1),
                        mg1.gridSpacing(2),
                        dx2[0],
                        dx2[1],
                        dx2[2],
                        mg2.gridSpacing(0),
                        mg2.gridSpacing(1),
                        mg2.gridSpacing(2),
                        t,    
                        (real &)tz,  // twilight zone pointer
                        dt,    
                        epsGrid(grid1),
                        muGrid(grid1),   
                        cGrid(grid1),    
                        epsGrid(grid2),  
                        muGrid(grid2),   
                        cGrid(grid2),
                        omegaForInterfaceIteration,
                        0., // return value averageInterfaceConvergenceRate
                        0.,  // return value maxFinalResidual
                        rtolForInterfaceIterations,
                        atolForInterfaceIterations
                    };
          // work space: 
                    real *rwk=interface.rwk;
                    int *iwk=interface.iwk;
                    assert( rwk!=NULL && iwk!=NULL );
                    const int ndf = max(interface.ndf1,interface.ndf2); 
          // assign pointers into the work spaces
                    int pa2=0,pa4=0,pa8=0, pipvt2=0,pipvt4=0,pipvt8=0;
                    if( bcOrderOfAccuracy==2 )
                    {
                        if( !( mg1.numberOfDimensions()==3 ) )  // new 3D doesn't use work-space yet
                        {
                            pa2=0; 
                            pa4=pa2 + 2*2*2*ndf;
                            pa8=0;  // not used
                            pipvt2=0;
                            pipvt4=pipvt2 + 2*ndf; 
                            pipvt8=0;
                        }
                    }
                    else if( bcOrderOfAccuracy==4 )
                    {
                        pa2=0; // not used
                        pa4=0;
                        pa8=pa4+4*4*2*ndf;
                        pipvt2=0;
                        pipvt4=0;
                        pipvt8=pipvt4+4*ndf;
                    }
                    if( mg1.numberOfDimensions()==2 && !useNewInterfaceRoutines )
                    {
            // OLD interface routines
                        interfaceMaxwell( mg1.numberOfDimensions(), 
                                      		      u1CopyLocal.getBase(0),u1CopyLocal.getBound(0),
                                      		      u1CopyLocal.getBase(1),u1CopyLocal.getBound(1),
                                      		      u1CopyLocal.getBase(2),u1CopyLocal.getBound(2),
          		      // note: use grid2 mesh data here  ASSUMES GRIDS MATCH *FIX ME*
                                // mg2.gridIndexRange(0,0), *u1CopyLocal.getDataPointer(), *mask2p,*prsxy2, *pxy2, bc2Local(0,0),
                                // fixed version: but note bc2Local 
                                                                mg2.gridIndexRange(0,0), *u1CopyLocal.getDataPointer(), *pmask1b,*prsxy1b, *pxy1b, bc2Local(0,0),
                                      		      u2Local.getBase(0),u2Local.getBound(0),
                                      		      u2Local.getBase(1),u2Local.getBound(1),
                                      		      u2Local.getBase(2),u2Local.getBound(2),
                                      		      mg2.gridIndexRange(0,0), *u2p, *mask2p,*prsxy2, *pxy2, bc2Local(0,0), 
                                    		    ipar[0], rpar[0], 
                                    		    rwk[pa2],rwk[pa4],rwk[pa8], iwk[pipvt2],iwk[pipvt4],iwk[pipvt8],
                                    		    ierr );
                    }
                    else
                    {
            // new interface routines -- 2D versions are done here too
            // interface3dMaxwell( mg1.numberOfDimensions(), 
                        interfaceOpt( mg1.numberOfDimensions(), 
                                        		        u1CopyLocal.getBase(0),u1CopyLocal.getBound(0),
                                        		        u1CopyLocal.getBase(1),u1CopyLocal.getBound(1),
                                        		        u1CopyLocal.getBase(2),u1CopyLocal.getBound(2),
          		        // note: use grid2 mesh data here  ASSUMES GRIDS MATCH *FIX ME*
                                  // mg2.gridIndexRange(0,0), *u1CopyLocal.getDataPointer(), *mask2p,*prsxy2, *pxy2, bc2Local(0,0),
                                  // fixed version: but note bc2Local 
                                                                    mg2.gridIndexRange(0,0), 
                                                                    *u1CopyLocal.getDataPointer(),   
                                  // *u1np,
                                                                    *u1nCopyLocal.getDataPointer(),
                                                                    *u1mp, *v1bp, 
                                                                    *pmask1b,*prsxy1b,*pxy1b, 
                                  // *p1ptr,*p1nptr,
                                                                    *p1CopyLocal.getDataPointer(),
                                                                    *p1nCopyLocal.getDataPointer(),
                                                                    *p1mptr,   
                                                                    *q1ptr,*q1nptr,*q1mptr,   
                                                                    bc2Local(0,0),
                                  // --- serial case ---
                                        		        u2Local.getBase(0),u2Local.getBound(0),
                                        		        u2Local.getBase(1),u2Local.getBound(1),
                                        		        u2Local.getBase(2),u2Local.getBound(2),
                                                                    mg2.gridIndexRange(0,0), 
                                                                    *u2p, *u2np, *u2mp, *v2p,
                                                                    *mask2p,*prsxy2, *pxy2, 
                                                                    *p2ptr,*p2nptr,*p2mptr, 
                                                                    *q2ptr,*q2nptr,*q2mptr, 
                                                                    bc2Local(0,0), 
                                      		      ipar[0], rpar[0], 
                                      		      rwk[pa2],rwk[pa4],rwk[pa8], iwk[pipvt2],iwk[pipvt4],iwk[pipvt8],
                                      		      ierr );
                    }
          // sosup solves for E and E.t so we need to divide some counts by 2: 
                    const real numSolvesPerStep = method==Maxwell::sosup ? 2. : 1.;
                    interface.totalInterfaceIterations+=ipar[43]/numSolvesPerStep; // counts total number of interface iterations
                    interface.averageInterfaceConvergenceRate+=rpar[22]/numSolvesPerStep;
                    interface.maxFinalResidual=max(rpar[23],interface.maxFinalResidual);
                    interface.averageFinalResidual+=rpar[23]/numSolvesPerStep;  // keeps sums of residuals 
        	  }
        	  if( debug & 16 )
                        ::display(u1Local,sPrintF("u2Local after assignOptInterface t=%8.2e",t),pDebugFile,"%5.2f ");
      	}
            #else
        // serial
          // Each grid may or may not have dispersion model: 
                    const Maxwell::DispersionModelEnum dispersionModel1 = dmp1.numberOfPolarizationVectors>0 ? dispersionModel : Maxwell::noDispersion;
                    const Maxwell::DispersionModelEnum dispersionModel2 = dmp2.numberOfPolarizationVectors>0 ? dispersionModel : Maxwell::noDispersion;
                    const int nonlinearModel1 = dmp1.getNonlinearModel();
                    const int nonlinearModel2 = dmp2.getNonlinearModel();
                    int ipar[]={ //
                        side1, dir1, grid1,         // keep side1,dir1 since we don't reverse the points.
                            n1a,n1b,n2a,n2b,n3a,n3b,
                        side2, dir2, grid2,         //  keep side2,dir2 since we don't reverse the points.
                            m1a,m1b,m2a,m2b,m3a,m3b,
                        gridType,            
                        bcOrderOfAccuracy,
                        orderOfExtrapolation,
                        useForcing,          
                        ex,                  
                        ey,                  
                        ez,                  
                        hx ,                 
                        hy,                  
                        hz,                  
                        (int)solveForElectricField,          
                        (int)solveForMagneticField,          
                        useWhereMask,       
                        debug,
                        numberOfIterationsForInterfaceBC,
                        materialInterfaceOption,
                        (int)interface.initialized,
                        myid,
                        parallel,
                        (int)forcingOption,
                        interfaceEquationsOption,
                        (int)assignInterfaceValues,
                        (int)assignInterfaceGhostValues,
                        mx.dbase.get<int>("setDivergenceAtInterfaces"),
                        mx.dbase.get<int>("useImpedanceInterfaceProjection"),
                        0,   // numberOfInterfaceIterationsUsed : returned value ipar[43]
                        dispersionModel1,   // ipar[44]
                        dispersionModel2,   // ipar[45]
                        pxc,                // ipar[46]
                        knownSolutionOption,
                        useJacobiUpdate,
                        nonlinearModel1,      // ipar[49]
                        nonlinearModel2,      // ipar[50]
                        numParallelGhost,     // ipar[51]
                        internalGhostBC       // ipar[52]
                    };
                const real & rtolForInterfaceIterations = mx.dbase.get<real>("rtolForInterfaceIterations");
                const real & atolForInterfaceIterations = mx.dbase.get<real>("atolForInterfaceIterations");
                real rpar[]={ //
                    dx1[0],
                    dx1[1],
                    dx1[2],
                    mg1.gridSpacing(0),
                    mg1.gridSpacing(1),
                    mg1.gridSpacing(2),
                    dx2[0],
                    dx2[1],
                    dx2[2],
                    mg2.gridSpacing(0),
                    mg2.gridSpacing(1),
                    mg2.gridSpacing(2),
                    t,    
                    (real &)tz,  // twilight zone pointer
                    dt,    
                    epsGrid(grid1),
                    muGrid(grid1),   
                    cGrid(grid1),    
                    epsGrid(grid2),  
                    muGrid(grid2),   
                    cGrid(grid2),
                    omegaForInterfaceIteration,
                    0., // return value averageInterfaceConvergenceRate
                    0.,  // return value maxFinalResidual
                    rtolForInterfaceIterations,
                    atolForInterfaceIterations
                };
        // work space: 
                real *rwk=interface.rwk;
                int *iwk=interface.iwk;
                assert( rwk!=NULL && iwk!=NULL );
                const int ndf = max(interface.ndf1,interface.ndf2); 
        // assign pointers into the work spaces
                int pa2=0,pa4=0,pa8=0, pipvt2=0,pipvt4=0,pipvt8=0;
                if( bcOrderOfAccuracy==2 )
                {
                    if( !( mg1.numberOfDimensions()==3 ) )  // new 3D doesn't use work-space yet
                    {
                        pa2=0; 
                        pa4=pa2 + 2*2*2*ndf;
                        pa8=0;  // not used
                        pipvt2=0;
                        pipvt4=pipvt2 + 2*ndf; 
                        pipvt8=0;
                    }
                }
                else if( bcOrderOfAccuracy==4 )
                {
                    pa2=0; // not used
                    pa4=0;
                    pa8=pa4+4*4*2*ndf;
                    pipvt2=0;
                    pipvt4=0;
                    pipvt8=pipvt4+4*ndf;
                }
                if( mg1.numberOfDimensions()==2 && !useNewInterfaceRoutines )
                {
          // OLD interface routines
                    interfaceMaxwell( mg1.numberOfDimensions(), 
                                  		      u1Local.getBase(0),u1Local.getBound(0),
                                  		      u1Local.getBase(1),u1Local.getBound(1),
                                  		      u1Local.getBase(2),u1Local.getBound(2),
                                  		      mg1.gridIndexRange(0,0), *u1p, *mask1p,*prsxy1, *pxy1, bc1Local(0,0), 
                                  		      u2Local.getBase(0),u2Local.getBound(0),
                                  		      u2Local.getBase(1),u2Local.getBound(1),
                                  		      u2Local.getBase(2),u2Local.getBound(2),
                                  		      mg2.gridIndexRange(0,0), *u2p, *mask2p,*prsxy2, *pxy2, bc2Local(0,0), 
                                		    ipar[0], rpar[0], 
                                		    rwk[pa2],rwk[pa4],rwk[pa8], iwk[pipvt2],iwk[pipvt4],iwk[pipvt8],
                                		    ierr );
                }
                else
                {
          // new interface routines -- 2D versions are done here too
          // interface3dMaxwell( mg1.numberOfDimensions(), 
                    interfaceOpt( mg1.numberOfDimensions(), 
                                // --- serial case ---
                                    		        u1Local.getBase(0),u1Local.getBound(0),
                                    		        u1Local.getBase(1),u1Local.getBound(1),
                                    		        u1Local.getBase(2),u1Local.getBound(2),
                                                                mg1.gridIndexRange(0,0), 
                                                                *u1p, *u1np, *u1mp, *v1p,
                                                                *mask1p,*prsxy1, *pxy1, 
                                                                *p1ptr,*p1nptr,*p1mptr,  
                                                                *q1ptr,*q1nptr,*q1mptr,   
                                                                bc1Local(0,0), 
                                // --- serial case ---
                                    		        u2Local.getBase(0),u2Local.getBound(0),
                                    		        u2Local.getBase(1),u2Local.getBound(1),
                                    		        u2Local.getBase(2),u2Local.getBound(2),
                                                                mg2.gridIndexRange(0,0), 
                                                                *u2p, *u2np, *u2mp, *v2p,
                                                                *mask2p,*prsxy2, *pxy2, 
                                                                *p2ptr,*p2nptr,*p2mptr, 
                                                                *q2ptr,*q2nptr,*q2mptr, 
                                                                bc2Local(0,0), 
                                  		      ipar[0], rpar[0], 
                                  		      rwk[pa2],rwk[pa4],rwk[pa8], iwk[pipvt2],iwk[pipvt4],iwk[pipvt8],
                                  		      ierr );
                }
        // sosup solves for E and E.t so we need to divide some counts by 2: 
                const real numSolvesPerStep = method==Maxwell::sosup ? 2. : 1.;
                interface.totalInterfaceIterations+=ipar[43]/numSolvesPerStep; // counts total number of interface iterations
                interface.averageInterfaceConvergenceRate+=rpar[22]/numSolvesPerStep;
                interface.maxFinalResidual=max(rpar[23],interface.maxFinalResidual);
                interface.averageFinalResidual+=rpar[23]/numSolvesPerStep;  // keeps sums of residuals 
                if( method==Maxwell::sosup )
      	{ // for sosup we assign E.t
          // printF("Assign interface values for sosup: E.t at t=%9.3e\n",t);
            // Each grid may or may not have dispersion model: 
                        const Maxwell::DispersionModelEnum dispersionModel1 = dmp1.numberOfPolarizationVectors>0 ? dispersionModel : Maxwell::noDispersion;
                        const Maxwell::DispersionModelEnum dispersionModel2 = dmp2.numberOfPolarizationVectors>0 ? dispersionModel : Maxwell::noDispersion;
                        const int nonlinearModel1 = dmp1.getNonlinearModel();
                        const int nonlinearModel2 = dmp2.getNonlinearModel();
                        int ipar[]={ //
                            side1, dir1, grid1,         // keep side1,dir1 since we don't reverse the points.
                                n1a,n1b,n2a,n2b,n3a,n3b,
                            side2, dir2, grid2,         //  keep side2,dir2 since we don't reverse the points.
                                m1a,m1b,m2a,m2b,m3a,m3b,
                            gridType,            
                            bcOrderOfAccuracy,
                            orderOfExtrapolation,
                            useForcing,          
                            ext,                  
                            eyt,                  
                            ezt,                  
                            hxt ,                 
                            hyt,                  
                            hzt,                  
                            (int)solveForElectricField,          
                            (int)solveForMagneticField,          
                            useWhereMask,       
                            debug,
                            numberOfIterationsForInterfaceBC,
                            materialInterfaceOption,
                            (int)interface.initialized,
                            myid,
                            parallel,
                            (int)forcingOption,
                            interfaceEquationsOption,
                            (int)assignInterfaceValues,
                            (int)assignInterfaceGhostValues,
                            mx.dbase.get<int>("setDivergenceAtInterfaces"),
                            mx.dbase.get<int>("useImpedanceInterfaceProjection"),
                            0,   // numberOfInterfaceIterationsUsed : returned value ipar[43]
                            dispersionModel1,   // ipar[44]
                            dispersionModel2,   // ipar[45]
                            pxc,                // ipar[46]
                            knownSolutionOption,
                            useJacobiUpdate,
                            nonlinearModel1,      // ipar[49]
                            nonlinearModel2,      // ipar[50]
                            numParallelGhost,     // ipar[51]
                            internalGhostBC       // ipar[52]
                        };
                    const real & rtolForInterfaceIterations = mx.dbase.get<real>("rtolForInterfaceIterations");
                    const real & atolForInterfaceIterations = mx.dbase.get<real>("atolForInterfaceIterations");
                    real rpar[]={ //
                        dx1[0],
                        dx1[1],
                        dx1[2],
                        mg1.gridSpacing(0),
                        mg1.gridSpacing(1),
                        mg1.gridSpacing(2),
                        dx2[0],
                        dx2[1],
                        dx2[2],
                        mg2.gridSpacing(0),
                        mg2.gridSpacing(1),
                        mg2.gridSpacing(2),
                        t,    
                        (real &)tz,  // twilight zone pointer
                        dt,    
                        epsGrid(grid1),
                        muGrid(grid1),   
                        cGrid(grid1),    
                        epsGrid(grid2),  
                        muGrid(grid2),   
                        cGrid(grid2),
                        omegaForInterfaceIteration,
                        0., // return value averageInterfaceConvergenceRate
                        0.,  // return value maxFinalResidual
                        rtolForInterfaceIterations,
                        atolForInterfaceIterations
                    };
          // work space: 
                    real *rwk=interface.rwk;
                    int *iwk=interface.iwk;
                    assert( rwk!=NULL && iwk!=NULL );
                    const int ndf = max(interface.ndf1,interface.ndf2); 
          // assign pointers into the work spaces
                    int pa2=0,pa4=0,pa8=0, pipvt2=0,pipvt4=0,pipvt8=0;
                    if( bcOrderOfAccuracy==2 )
                    {
                        if( !( mg1.numberOfDimensions()==3 ) )  // new 3D doesn't use work-space yet
                        {
                            pa2=0; 
                            pa4=pa2 + 2*2*2*ndf;
                            pa8=0;  // not used
                            pipvt2=0;
                            pipvt4=pipvt2 + 2*ndf; 
                            pipvt8=0;
                        }
                    }
                    else if( bcOrderOfAccuracy==4 )
                    {
                        pa2=0; // not used
                        pa4=0;
                        pa8=pa4+4*4*2*ndf;
                        pipvt2=0;
                        pipvt4=0;
                        pipvt8=pipvt4+4*ndf;
                    }
                    if( mg1.numberOfDimensions()==2 && !useNewInterfaceRoutines )
                    {
            // OLD interface routines
                        interfaceMaxwell( mg1.numberOfDimensions(), 
                                      		      u1Local.getBase(0),u1Local.getBound(0),
                                      		      u1Local.getBase(1),u1Local.getBound(1),
                                      		      u1Local.getBase(2),u1Local.getBound(2),
                                      		      mg1.gridIndexRange(0,0), *u1p, *mask1p,*prsxy1, *pxy1, bc1Local(0,0), 
                                      		      u2Local.getBase(0),u2Local.getBound(0),
                                      		      u2Local.getBase(1),u2Local.getBound(1),
                                      		      u2Local.getBase(2),u2Local.getBound(2),
                                      		      mg2.gridIndexRange(0,0), *u2p, *mask2p,*prsxy2, *pxy2, bc2Local(0,0), 
                                    		    ipar[0], rpar[0], 
                                    		    rwk[pa2],rwk[pa4],rwk[pa8], iwk[pipvt2],iwk[pipvt4],iwk[pipvt8],
                                    		    ierr );
                    }
                    else
                    {
            // new interface routines -- 2D versions are done here too
            // interface3dMaxwell( mg1.numberOfDimensions(), 
                        interfaceOpt( mg1.numberOfDimensions(), 
                                  // --- serial case ---
                                        		        u1Local.getBase(0),u1Local.getBound(0),
                                        		        u1Local.getBase(1),u1Local.getBound(1),
                                        		        u1Local.getBase(2),u1Local.getBound(2),
                                                                    mg1.gridIndexRange(0,0), 
                                                                    *u1p, *u1np, *u1mp, *v1p,
                                                                    *mask1p,*prsxy1, *pxy1, 
                                                                    *p1ptr,*p1nptr,*p1mptr,  
                                                                    *q1ptr,*q1nptr,*q1mptr,   
                                                                    bc1Local(0,0), 
                                  // --- serial case ---
                                        		        u2Local.getBase(0),u2Local.getBound(0),
                                        		        u2Local.getBase(1),u2Local.getBound(1),
                                        		        u2Local.getBase(2),u2Local.getBound(2),
                                                                    mg2.gridIndexRange(0,0), 
                                                                    *u2p, *u2np, *u2mp, *v2p,
                                                                    *mask2p,*prsxy2, *pxy2, 
                                                                    *p2ptr,*p2nptr,*p2mptr, 
                                                                    *q2ptr,*q2nptr,*q2mptr, 
                                                                    bc2Local(0,0), 
                                      		      ipar[0], rpar[0], 
                                      		      rwk[pa2],rwk[pa4],rwk[pa8], iwk[pipvt2],iwk[pipvt4],iwk[pipvt8],
                                      		      ierr );
                    }
          // sosup solves for E and E.t so we need to divide some counts by 2: 
                    const real numSolvesPerStep = method==Maxwell::sosup ? 2. : 1.;
                    interface.totalInterfaceIterations+=ipar[43]/numSolvesPerStep; // counts total number of interface iterations
                    interface.averageInterfaceConvergenceRate+=rpar[22]/numSolvesPerStep;
                    interface.maxFinalResidual=max(rpar[23],interface.maxFinalResidual);
                    interface.averageFinalResidual+=rpar[23]/numSolvesPerStep;  // keeps sums of residuals 
      	}
      	
            #endif

            
        }
        else
        {
      // *** test the new interface routines for order of accuracy >= 6 

            printF("Call new interface routines for order=%i at t=%9.3e\n",orderOfAccuracyInSpace,t);

            assert( mg1.numberOfDimensions()==2 );

      // macro: 
        // Each grid may or may not have dispersion model: 
                const Maxwell::DispersionModelEnum dispersionModel1 = dmp1.numberOfPolarizationVectors>0 ? dispersionModel : Maxwell::noDispersion;
                const Maxwell::DispersionModelEnum dispersionModel2 = dmp2.numberOfPolarizationVectors>0 ? dispersionModel : Maxwell::noDispersion;
                const int nonlinearModel1 = dmp1.getNonlinearModel();
                const int nonlinearModel2 = dmp2.getNonlinearModel();
                int ipar[]={ //
                    side1, dir1, grid1,         // keep side1,dir1 since we don't reverse the points.
                        n1a,n1b,n2a,n2b,n3a,n3b,
                    side2, dir2, grid2,         //  keep side2,dir2 since we don't reverse the points.
                        m1a,m1b,m2a,m2b,m3a,m3b,
                    gridType,            
                    bcOrderOfAccuracy,
                    orderOfExtrapolation,
                    useForcing,          
                    ex,                  
                    ey,                  
                    ez,                  
                    hx ,                 
                    hy,                  
                    hz,                  
                    (int)solveForElectricField,          
                    (int)solveForMagneticField,          
                    useWhereMask,       
                    debug,
                    numberOfIterationsForInterfaceBC,
                    materialInterfaceOption,
                    (int)interface.initialized,
                    myid,
                    parallel,
                    (int)forcingOption,
                    interfaceEquationsOption,
                    (int)assignInterfaceValues,
                    (int)assignInterfaceGhostValues,
                    mx.dbase.get<int>("setDivergenceAtInterfaces"),
                    mx.dbase.get<int>("useImpedanceInterfaceProjection"),
                    0,   // numberOfInterfaceIterationsUsed : returned value ipar[43]
                    dispersionModel1,   // ipar[44]
                    dispersionModel2,   // ipar[45]
                    pxc,                // ipar[46]
                    knownSolutionOption,
                    useJacobiUpdate,
                    nonlinearModel1,      // ipar[49]
                    nonlinearModel2,      // ipar[50]
                    numParallelGhost,     // ipar[51]
                    internalGhostBC       // ipar[52]
                };
            const real & rtolForInterfaceIterations = mx.dbase.get<real>("rtolForInterfaceIterations");
            const real & atolForInterfaceIterations = mx.dbase.get<real>("atolForInterfaceIterations");
            real rpar[]={ //
                dx1[0],
                dx1[1],
                dx1[2],
                mg1.gridSpacing(0),
                mg1.gridSpacing(1),
                mg1.gridSpacing(2),
                dx2[0],
                dx2[1],
                dx2[2],
                mg2.gridSpacing(0),
                mg2.gridSpacing(1),
                mg2.gridSpacing(2),
                t,    
                (real &)tz,  // twilight zone pointer
                dt,    
                epsGrid(grid1),
                muGrid(grid1),   
                cGrid(grid1),    
                epsGrid(grid2),  
                muGrid(grid2),   
                cGrid(grid2),
                omegaForInterfaceIteration,
                0., // return value averageInterfaceConvergenceRate
                0.,  // return value maxFinalResidual
                rtolForInterfaceIterations,
                atolForInterfaceIterations
            };

            newInterfaceMaxwell( mg1.numberOfDimensions(), 
                     			   u1Local.getBase(0),u1Local.getBound(0),
                     			   u1Local.getBase(1),u1Local.getBound(1),
                     			   u1Local.getBase(2),u1Local.getBound(2),
                     			   mg1.gridIndexRange(0,0), *u1p, *mask1p,*prsxy1, *pxy1, bc1(0,0), 
                     			   u2Local.getBase(0),u2Local.getBound(0),
                     			   u2Local.getBase(1),u2Local.getBound(1),
                     			   u2Local.getBase(2),u2Local.getBound(2),
                     			   mg2.gridIndexRange(0,0), *u2p, *mask2p,*prsxy2, *pxy2, bc2(0,0), 
                     			   ipar[0], rpar[0], ierr );

        }
              		  
    // wait to set initialized=true until ghost values have been assigned the first time : 
        if( assignInterfaceGhostValues )
            interface.initialized=true;
        
    // In some cases we have an optimized periodic update implemented
        bool updatePeriodic = mg1.numberOfDimensions()==3; 
        #ifdef USE_PPP
      // in parallel we cannot use the optimized periodic update in the interface routines.
            updatePeriodic=true;
        #endif
        if(  updatePeriodic )
        {
      // printF("After assign interfaces: periodic update...\n");
            
            cgfields[next][grid1].periodicUpdate(); 
            cgfields[next][grid2].periodicUpdate();
            u1.updateGhostBoundaries(); // *wdh* 081127 
            u2.updateGhostBoundaries();
        }
        
        if( debug & 4 )
            fprintf(pDebugFile," **** After assigning interfaces t=%8.2e, tz=%i\n",t,(tz==NULL ? 0 : 1));
        
        if( debug & 4 ) 
        {
            if( tz!=NULL && pDebugFile!=NULL )
            {
                realArray & x1 = mg1.center();
                realArray & x2 = mg2.center();
                #ifdef USE_PPP
                    realSerialArray x1Local; getLocalArrayWithGhostBoundaries(x1,x1Local);
                    realSerialArray x2Local; getLocalArrayWithGhostBoundaries(x2,x2Local);
                #else
                    realSerialArray & x1Local = x1;
                    realSerialArray & x2Local = x2;
                #endif
                		    
      	OGFunction & e = *tz;
                Range E(ex,ex+numberOfDimensions-1);
      	int isRectangular=0;
      	getGhostIndex(mg1.gridIndexRange(),side1,dir1,I1,I2,I3);

      	int includeGhost=0; 
      	bool ok = ParallelUtility::getLocalArrayBounds(u1,u1Local,I1,I2,I3,includeGhost);
      	if( ok )
      	{
        	  realSerialArray err(I1,I2,I3,E);
        	  realSerialArray ue(I1,I2,I3,E);
              		  
	  // err=u1(I1,I2,I3,ex)-e(mg1,I1,I2,I3,ex,t);


        	  e.gd( ue  ,x1Local,numberOfDimensions,isRectangular,0,0,0,0,I1,I2,I3,E,t);
        	  err=u1Local(I1,I2,I3,E)-ue(I1,I2,I3,E);
        	  ::display(u1Local(I1,I2,I3,E),sPrintF("u1 (ghost) (grid1=%i) after interface, t=%e",grid1,t),pDebugFile,"%8.1e ");
        	  ::display(err,sPrintF("err in u1 (ghost) (grid1=%i) after interface, t=%e",grid1,t),pDebugFile,"%8.1e ");

	  // getGhostIndex(mg1.gridIndexRange(),side1,dir1,I1,I2,I3,-1);
	  // err=u1(I1,I2,I3,ex)-e(mg1,I1,I2,I3,ex,t);
	  // ::display(err,sPrintF("err in u1 (ex,line 1) after interface, t=%e",t),pDebugFile,"%8.1e ");
      	}
      	

      	getGhostIndex(mg2.gridIndexRange(),side2,dir2,I1,I2,I3);
      	ok = ParallelUtility::getLocalArrayBounds(u2,u2Local,I1,I2,I3,includeGhost);
      	if( ok )
      	{
        	  realSerialArray err(I1,I2,I3,E);
        	  realSerialArray ue(I1,I2,I3,E);
                  
        	  e.gd( ue  ,x2Local,numberOfDimensions,isRectangular,0,0,0,0,I1,I2,I3,E,t);
        	  err=u2Local(I1,I2,I3,E)-ue(I1,I2,I3,E);
        	  ::display(u2Local(I1,I2,I3,E),sPrintF("u2 (ghost) (grid2=%i) after interface, t=%e",grid2,t),pDebugFile,"%8.1e ");
        	  ::display(err,sPrintF("err in u2 (ghost) (grid2=%i) after interface, t=%e",grid2,t),pDebugFile,"%8.1e ");

	  // err=u2(I1,I2,I3,ex)-e(mg2,I1,I2,I3,ex,t);
	  // ::display(err,sPrintF("err in u2 (ex,ghost) after interface, t=%e",t),pDebugFile,"%8.1e ");
	  // getGhostIndex(mg2.gridIndexRange(),side2,dir2,I1,I2,I3,-1);
	  // err=u2(I1,I2,I3,ex)-e(mg2,I1,I2,I3,ex,t);
	  // ::display(err,sPrintF("err in u2 (ex,line 1) after interface, t=%e",t),pDebugFile,"%8.1e ");
      	}
      	
            }
        }
        
        if( debug & 4 ) 
        {
            ::display(u1Local,sPrintF("u1Local after assignOptInterface grid1=%i, t=%8.2e",grid1,t),pDebugFile,"%9.2e ");
            ::display(u2Local,sPrintF("u2Local after assignOptInterface grid2=%i, t=%8.2e",grid2,t),pDebugFile,"%9.2e ");
        }
        if( debug & 4 )
        {
            checkParallel.checkDiff(u1Local,sPrintF("u1Local after assignOptInterface grid1=%i, t=%8.2e",grid1,t),pDebugFile);
            checkParallel.checkDiff(u2Local,sPrintF("u2Local after assignOptInterface grid2=%i, t=%8.2e",grid2,t),pDebugFile);
        }
        

    }
    

    printF("\n <<<<<<<<<<< InterfaceMaxwell:: done computeStencilCoefficients.\n\n");

    return 0;
}


/* ---
void InterfaceMaxwell::
initializeInterfaces( Maxwell & mx )
// =====================================================================================================
//   /Description:
//      Find the interfaces and initialize the work-space. 
// =====================================================================================================
{
    CompositeGrid *& cgp = mx.cgp;
    int & debug = mx.debug;
    int & orderOfAccuracyInSpace = mx.orderOfAccuracyInSpace;
    
}
--- */


