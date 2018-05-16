// ================================================================================================
/// Macro: assignTractionBC
///  Assign a free-surface or traction free BC
// ================================================================================================
#beginMacro assignTractionBC()
{
 // twilight always needs the vertex: 
  bool vertexNeeded = !isRectangular || parameters.dbase.get<bool >("twilightZoneFlow");

  OV_GET_SERIAL_ARRAY(real,u,uLocal);
  real *pu = uLocal.getDataPointer();

  OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
  const int *pmask = maskLocal.getDataPointer();

  real temp;
  real *pxy=&temp, *prsxy=&temp;
  if( !isRectangular )
  {
    OV_GET_SERIAL_ARRAY(real,mg.inverseVertexDerivative(),rxLocal);
    prsxy=rxLocal.getDataPointer();
  }
  if( vertexNeeded )
  {
    OV_GET_SERIAL_ARRAY(real,mg.vertex(),xLocal);
    pxy=xLocal.getDataPointer();
  }
  
  // Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
  // getIndex(mg.gridIndexRange(),I1,I2,I3);
  // int includeGhost=1;
  // bool ok = ParallelUtility::getLocalArrayBounds(u0,uLocal,I1,I2,I3,includeGhost);


  OV_GET_SERIAL_ARRAY_CONDITIONAL(real,gridVelocity,gridVelocityLocal,parameters.gridIsMoving(grid));
  realSerialArray gtt;  // holds the boundary acceleration
  // if( parameters.gridIsMoving(grid) )
  // {
  //   // Moving grid problem: compute the acceleration on the boundary

  //   MovingGrids & movingGrids = parameters.dbase.get<MovingGrids >("movingGrids");
  //   int boundaryAccelerationOption=4; // return gtt
  //   if( ok )
  //   {
  //     // note: insbc4 assumes gtt has the same first 3 dimensions as uLocal
  //     gtt.redim(uLocal.dimension(0),uLocal.dimension(1),uLocal.dimension(2),numberOfDimensions);
  //     gtt=0.;
  //     movingGrids.getBoundaryAcceleration( mg, gtt, grid, t, boundaryAccelerationOption );
  //   }
    
  // }
  real *pgt  = gridVelocityLocal.getDataPointer(); // pointer to the grid velocity, g'
  real *pgtt = gtt.getDataPointer();               // pointer to the grid acceleration, g''
  


  // check -- is it ok to use gid instead of ir?
  IntegerArray indexRangeLocal(2,3), dimLocal(2,3), bcLocal(2,3);
  ParallelGridUtility::getLocalIndexBoundsAndBoundaryConditions( u,indexRangeLocal,dimLocal,bcLocal ); 


  // *wdh* 110311 - Add Boussinesq terms and boundary conditions for T
  const InsParameters::PDEModel & pdeModel = parameters.dbase.get<InsParameters::PDEModel >("pdeModel");


  const bool assignTemperature = pdeModel==InsParameters::BoussinesqModel ||
                                 pdeModel==InsParameters::viscoPlasticModel;

  // -- We may adjust the artificial dissipation since the fourth-order dissipation doesn't work too well ---
  bool useSecondOrderArtificialDiffusion = parameters.dbase.get<bool >("useSecondOrderArtificialDiffusion");
  bool useFourthOrderArtificialDiffusion = parameters.dbase.get<bool >("useFourthOrderArtificialDiffusion");
  
  real ad21 = parameters.dbase.get<real >("ad21");
  real ad22 = parameters.dbase.get<real >("ad22");
  real ad41 = parameters.dbase.get<real >("ad41");
  real ad42 = parameters.dbase.get<real >("ad42");

  const int np = max(Communication_Manager::Number_Of_Processors,1);

  int useWhereMask=false; // **NOTE** for  moving grids we may need to evaluate at more points than just mask >0 

  real ajs=1.;  // is this used?
  //   real ajs=getSignForJacobian(mg);
  real thermalExpansivity=1.;
  // parameters.dbase.get<ListOfShowFileParameters >("pdeParameters").getParameter("thermalExpansivity",thermalExpansivity)

  real dx[3]={1.,1.,1.};
  bool isRectangular=mg.isRectangular();
  if( isRectangular )
  {
    mg.getDeltaX(dx);
  }

  int gridType = isRectangular ? 0 : 1;
  int iparam[] ={parameters.dbase.get<int >("pc"),
		 parameters.dbase.get<int >("uc"),
		 parameters.dbase.get<int >("vc"),
		 parameters.dbase.get<int >("wc"),
		 parameters.dbase.get<int >("sc"),
		 grid,
		 gridType,
		 min(4,parameters.dbase.get<int >("orderOfAccuracy")),
		 (int)parameters.gridIsMoving(grid),
		 useWhereMask,
		 parameters.getGridIsImplicit(grid),
		 (int)parameters.dbase.get<Parameters::ImplicitMethod >("implicitMethod"),
		 (int)parameters.dbase.get<Parameters::ImplicitOption >("implicitOption"),
		 (int)parameters.isAxisymmetric(),
		 useSecondOrderArtificialDiffusion,
		 useFourthOrderArtificialDiffusion,
		 (int)parameters.dbase.get<bool >("twilightZoneFlow"),
		 np,
		 parameters.dbase.get<int>("outflowOption"),
                 parameters.dbase.get<int >("orderOfExtrapolationForOutflow"),
                 parameters.dbase.get<int >("debug"),
		 parameters.dbase.get<int >("myid"),
                 (int)assignTemperature,
                 parameters.dbase.get<int >("tc"),
                 parameters.dbase.get<int >("numberOfComponents")
  };

  const real adcPassiveScalar=1.; // coeff or linear artificial diffusion for the passive scalar ** add to params
  real gravity[3];

  real rparam[]={dx[0],dx[1],dx[2],
                 mg.gridSpacing(0),
                 mg.gridSpacing(1),
                 mg.gridSpacing(2),
                 parameters.dbase.get<real >("nu"),
                 t,
                 ad21, //kkc 101029 parameters.dbase.get<real >("ad21"),
                 ad22, //kkc 101029 parameters.dbase.get<real >("ad22"),
                 ad41, //kkc 101029 parameters.dbase.get<real >("ad41"),
                 ad42, //kkc 101029 parameters.dbase.get<real >("ad42"),
                 parameters.dbase.get<real >("nuPassiveScalar"),
                 adcPassiveScalar,
                 ajs,
                 gravity[0],
                 gravity[1],
                 gravity[2],
                 thermalExpansivity,
                 (real &)(parameters.dbase.get<OGFunction* >("exactSolution")), // pointer to TZ
                 REAL_MIN
  };
    
  int ierr=0;
  int bcOption=0;
  insTractionBC(bcOption,
	 mg.numberOfDimensions(),
	 uLocal.getBase(0),uLocal.getBound(0),uLocal.getBase(1),uLocal.getBound(1),
	 uLocal.getBase(2),uLocal.getBound(2),uLocal.getBase(3),uLocal.getBound(3),
	 iparam[0],rparam[0], *pu, *pmask, *pxy, *prsxy, *pgt, *pgtt,
	 bcLocal(0,0), indexRangeLocal(0,0), ierr ) ;
}
#endMacro
