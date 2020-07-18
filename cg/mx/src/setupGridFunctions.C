#include "Maxwell.h"
#include "SquareMapping.h"
#include "BoxMapping.h"
#include "AnnulusMapping.h"
#include "MatrixTransform.h"
#include "DataPointMapping.h"
#include "CompositeGridOperators.h"
#include "display.h"
#include "UnstructuredMapping.h"
#include "ParallelUtility.h"
#include "GridStatistics.h"
#include "DispersiveMaterialParameters.h"

#include "ULink.h"

extern bool verifyUnstructuredConnectivity( UnstructuredMapping &umap, bool verbose );

#define FOR_3D(i1,i2,i3,I1,I2,I3) \
int I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase();  \
int I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(); \
for(i3=I3Base; i3<=I3Bound; i3++) \
for(i2=I2Base; i2<=I2Bound; i2++) \
for(i1=I1Base; i1<=I1Bound; i1++)

#define FOR_3(i1,i2,i3,I1,I2,I3) \
I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase();  \
I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(); \
for(i3=I3Base; i3<=I3Bound; i3++) \
for(i2=I2Base; i2<=I2Bound; i2++) \
for(i1=I1Base; i1<=I1Bound; i1++)


//! Setup and initialization. Build the solution fields. 
int Maxwell::
setupGridFunctions()
// ===================================================================================
//    Build grid functions
// ===================================================================================
{
  real time0=getCPU();

  printF("\n >>>>>>>>>>>>>>>> Maxwell: ENTERING setupGridFunctions\n");

  const int & solveForAllFields = dbase.get<int>("solveForAllFields");
  printF(" solveForAllFields=%i, solveForMagneticField=%i\n",(int)solveForAllFields,(int)solveForMagneticField);
  
  const int & np = dbase.get<int>("np"); // number of processors 

  if( np > 1  )
  {
    // --- make sure there are enough parallel ghost points ---
    
    const int & numberOfParallelGhost = dbase.get<int>("numberOfParallelGhost");  // number of parallel ghost points 
    int numParGhostNeeded = orderOfAccuracyInSpace/2;

    if( method==bamx && artificialDissipation>0. ) numParGhostNeeded++;

    const int useSosupDissipation = parameters.dbase.get<int>("useSosupDissipation");
    if( method==nfdtd && useSosupDissipation ) numParGhostNeeded++;

    if( numParGhostNeeded > numberOfParallelGhost )
    {
      printF("CgMx:ERROR: Not enough parallel ghost points. Current=%d\n",numberOfParallelGhost);
      printF("     orderOfAccuracyInSpace=%d requires %d ghost.\n",orderOfAccuracyInSpace,orderOfAccuracyInSpace/2);
      if( method==bamx && artificialDissipation>0. )
	printF(" Method BAMX + dissipation requires 1 additional parallel ghost.\n");
      if( method==nfdtd && useSosupDissipation )
	printF(" Method NFDTD + sosup dissipation requires 1 additional parallel ghost.\n");

      printF("Add the cgmx command line option \"-numberOfParallelGhost=%d\" \n",numParGhostNeeded);
      OV_ABORT("ERROR");
    }
  }
  


  assert( cgp!=NULL );
  CompositeGrid & cg= *cgp;
  const int numberOfDimensions = cg.numberOfDimensions();

  parameters.dbase.get<int>("numberOfDimensions")=numberOfDimensions;
  
  int grid,side,axis;
  for( grid=0; grid<cg.numberOfComponentGrids(); grid++ )
  {
    MappedGrid & mg = cg[grid];

    if( mg.getGridType()==MappedGrid::structuredGrid )
    {
      for( axis=0; axis<mg.numberOfDimensions(); axis++ )
	for( side=0; side<=1; side++ )
	{
	  int bc0 = mg.boundaryCondition(side,axis);
	  if( bc0>0 )
	  {
	    if( bcOption==useAllDirichletBoundaryConditions )
	    {
	      mg.setBoundaryCondition(side,axis,dirichlet);
	    }
	    else if( bcOption==useAllPerfectElectricalConductorBoundaryConditions )
	    {
	      mg.setBoundaryCondition(side,axis,perfectElectricalConductor);
	    }
	  }
      
	}
    }
    else
    {
      UnstructuredMapping & umap = (UnstructuredMapping &) mg.mapping().getMapping();
      int dDim = umap.getDomainDimension();
      int rDim = umap.getRangeDimension();

	
      UnstructuredMapping::EntityTypeEnum etype,etypeBdy;
      etype = UnstructuredMapping::EntityTypeEnum(dDim);
      etypeBdy = UnstructuredMapping::Edge;//UnstructuredMapping::EntityTypeEnum(dDim-1);

      string ghostCellTag = string("Ghost ")+UnstructuredMapping::EntityTypeStrings[etype].c_str();
      string ghostBdyTag = string("Ghost ")+UnstructuredMapping::EntityTypeStrings[etypeBdy].c_str();
      string bctag = string("__bcnum ")+UnstructuredMapping::EntityTypeStrings[etypeBdy].c_str();
      string perTag = string("periodic ")+UnstructuredMapping::EntityTypeStrings[etype].c_str();
      // set a tag on all the bounding entities to specify the pec bc
      std::string pecTag("PEC BC");
      if ( bcOption )
      {
	int bcv = dirichlet;
	    
	if ( bcOption==useAllPerfectElectricalConductorBoundaryConditions )
	  bcv = perfectElectricalConductor;

	UnstructuredMapping::tag_entity_iterator git = umap.tag_entity_begin(ghostCellTag);
	UnstructuredMapping::tag_entity_iterator git_end = umap.tag_entity_end(ghostCellTag);

	UnstructuredMappingAdjacencyIterator ebdy,ebdy_end;

	int nPec=0;
	for ( ; git!=git_end; git++ )
	{
	  string bctagc = string("__bcnum ")+UnstructuredMapping::EntityTypeStrings[etype].c_str();
	  //		umap.addTag(git->et,git->e,bctagc,(void *)bcv);
	  umap.setBC(git->et,git->e,bcv);

	  etypeBdy = UnstructuredMapping::Edge;//UnstructuredMapping::EntityTypeEnum(dDim-1);
	  string ghostBdyTag = string("Ghost ")+UnstructuredMapping::EntityTypeStrings[etypeBdy].c_str();
	  string bctag = string("__bcnum ")+UnstructuredMapping::EntityTypeStrings[etypeBdy].c_str();
	  ebdy_end = umap.adjacency_end(*git,etypeBdy);
	  for ( ebdy=umap.adjacency_begin(*git, etypeBdy); ebdy!=ebdy_end; ebdy++ )
	  {
	    if ( !umap.hasTag(etypeBdy,*ebdy,bctag) )
	    {
			
	      if ( bcOption==useAllPerfectElectricalConductorBoundaryConditions )
		umap.addTag(etypeBdy,*ebdy,pecTag,0);

	      //			umap.addTag(etypeBdy,*ebdy,bctag,(void *)bcv);
	      umap.setBC(etypeBdy,*ebdy,bcv);
	      nPec++;
	    }

	    if ( umap.hasTag(git->et, git->e, perTag ) )
	      umap.deleteTag(git->et, git->e, perTag );
	  }

	  if ( rDim==3 )
	  {
	    etypeBdy = UnstructuredMapping::Face;//UnstructuredMapping::EntityTypeEnum(dDim-1);
	    ghostBdyTag = string("Ghost ")+UnstructuredMapping::EntityTypeStrings[etypeBdy].c_str();
	    bctag = string("__bcnum ")+UnstructuredMapping::EntityTypeStrings[etypeBdy].c_str();
	    ebdy_end = umap.adjacency_end(*git,etypeBdy);
	    for ( ebdy=umap.adjacency_begin(*git, etypeBdy); ebdy!=ebdy_end; ebdy++ )
	    {
	      if (  !umap.hasTag(etypeBdy,*ebdy,bctag) )
	      {
			    
		if ( bcOption==useAllPerfectElectricalConductorBoundaryConditions )
		  umap.addTag(etypeBdy,*ebdy,pecTag,0);
			    
		//			    umap.addTag(etypeBdy,*ebdy,bctag,(void *)bcv);
		umap.setBC(etypeBdy,*ebdy,bcv);
		nPec++;
	      }
			
	      if ( umap.hasTag(git->et, git->e, perTag ) )
		umap.deleteTag(git->et, git->e, perTag );
	    }
	  }
	}

	if ( !nPec )
	{
	  string bdyEdgeTag = string("boundary ")+UnstructuredMapping::EntityTypeStrings[int(UnstructuredMapping::Edge)].c_str();
	  string bdyFaceTag = string("boundary ")+UnstructuredMapping::EntityTypeStrings[int(UnstructuredMapping::Face)].c_str();
	  string bdyVertTag = string("boundary ")+UnstructuredMapping::EntityTypeStrings[int(UnstructuredMapping::Vertex)].c_str();
		
	  UnstructuredMapping::tag_entity_iterator bit = umap.tag_entity_begin(bdyEdgeTag);
	  UnstructuredMapping::tag_entity_iterator bit_end = umap.tag_entity_end(bdyEdgeTag);
		
	  UnstructuredMappingAdjacencyIterator ebdy,ebdy_end;

	  string bctagc = string("__bcnum ")+UnstructuredMapping::EntityTypeStrings[int(UnstructuredMapping::Edge)].c_str();
	  for ( ; bit!=bit_end; bit++ )
	  {
	    if ( !umap.hasTag(bit->et,bit->e,bctagc) )
	    {
	      //			umap.addTag(bit->et,bit->e,bctagc,(void *)bcv);
	      umap.setBC(bit->et,bit->e,bcv);
	      if ( bcOption==useAllPerfectElectricalConductorBoundaryConditions )
		umap.addTag(bit->et,bit->e,pecTag,0);
	      nPec++;
	    }
		    
	  }

	  bit = umap.tag_entity_begin(bdyVertTag);
	  bit_end = umap.tag_entity_end(bdyVertTag);
	  bctagc = string("__bcnum ")+UnstructuredMapping::EntityTypeStrings[int(UnstructuredMapping::Vertex)].c_str();
	  for ( ; bit!=bit_end; bit++ )
	  {
	    if ( !umap.hasTag(bit->et,bit->e,bctagc) )
	    {
	      //			umap.addTag(bit->et,bit->e,bctagc,(void *)bcv);
	      umap.setBC(bit->et,bit->e,bcv);
	      if ( bcOption==useAllPerfectElectricalConductorBoundaryConditions )
		umap.addTag(bit->et,bit->e,pecTag,0);
	      nPec++;
	    }
		    
	  }

	  if ( rDim==3 )
	  {
	    bit = umap.tag_entity_begin(bdyFaceTag);
	    bit_end = umap.tag_entity_end(bdyFaceTag);
	    bctagc = string("__bcnum ")+UnstructuredMapping::EntityTypeStrings[int(UnstructuredMapping::Face)].c_str();
	    for ( ; bit!=bit_end; bit++ )
	    {
	      if ( !umap.hasTag(bit->et,bit->e,bctagc) )
	      {
		//			    umap.addTag(bit->et,bit->e,bctagc,(void *)bcv);
		umap.setBC(bit->et,bit->e,bcv);
		if ( bcOption==useAllPerfectElectricalConductorBoundaryConditions )
		  umap.addTag(bit->et,bit->e,pecTag,0);
		nPec++;
	      }
	    }
	  }
		
	}

	cout<<"unstructured grid "<<grid<<" has "<<nPec<<" non-periodic boundary entities"<<endl;
      }
      cout<<"unstructured grid "<<grid<<" has "<<std::distance(umap.tag_entity_begin(perTag),umap.tag_entity_end(perTag))<<" periodic boundary entities"<<endl;
	
    }
  
  }

  numberOfTimeLevels=2;  // keep this many levels of u
  if( method==nfdtd || method==yee || method==sosup || method==bamx )
  {
    if( usingPMLBoundaryConditions() ||
        true )  // for combined dissipation with advance
      numberOfTimeLevels=3;  // needed for PML and ABC's 
  }

  int & numberOfComponents      = dbase.get<int>("numberOfComponents");
  int & numberOfComponentsForTZ = dbase.get<int>("numberOfComponentsForTZ");
  
  numberOfComponents=-1;  // defines how many components are stored in the grid functions

  if( method!=nfdtd && method!=sosup )
  {// both fields are needed by the staggered schemes
    solveForElectricField = solveForMagneticField = true;
  }

  if( gridType!=compositeGrid )
  {
    // ****************************************
    // *********** single grid  ***************
    // ****************************************

    MappedGrid & mg = cg[0];
    Mapping & map = mg.mapping().getMapping();
  
    bool unstructured=map.getClassName()=="UnstructuredMapping";

    if( (method==nfdtd || method==sosup|| method==bamx ) && !mg.isRectangular() )
    {
      // The energy estimate may not be correct unless we ensure that the
      // jacobian derivatives are periodic
      mg.update( MappedGrid::THEinverseVertexDerivative );
      realMappedGridFunction & rx = mg.inverseVertexDerivative();
      rx.periodicUpdate();
    }

    if( method==defaultMethod )
    {
      if( mg.isRectangular() )
	method=yee;
      else
	method=dsi;
    }
  

    Range all;
    if( numberOfDimensions==2 )
    {
      numberOfComponents= method==yee || method==nfdtd || method==sosup || method==bamx ? 
	(int)numberOfComponentsRectangularGrid : 
	(int)numberOfComponentsCurvilinearGrid; // number of components
    }
    else
    {
      numberOfComponents=(int(solveForElectricField)+int(solveForMagneticField))*3;
    }

    if( method==sosup )
    {
      // With sosup we also advance/store the time-derivative of the fields
      numberOfComponents*= 2;
    }
    numberOfComponentsForTZ=numberOfComponents;
    

    //kkc 040304 ops apparently needed now for all schemes
    MappedGridOperators *ops = new MappedGridOperators(mg);
    this->op = ops; //kkc 040304 I guess we need this too?

    if( method==nfdtd || method==bamx )
    {      
      assert( method!=bamx );  // finish me 
      

      fields = new realMappedGridFunction [numberOfTimeLevels];
      for( int n=0; n<numberOfTimeLevels; n++ )
      {
	fields[n].updateToMatchGrid(mg,all,all,all,numberOfComponents);
	fields[n]=0.;
	      
	if( numberOfDimensions==2 )
	{
	  if( numberOfComponents==3 )
	  {
	    fields[n].setName("Ex",ex);
	    fields[n].setName("Ey",ey);
	    fields[n].setName("Hz",hz);
	  }
	  else
	  {
	    fields[n].setName("Ex10",ex10);
	    fields[n].setName("Ey10",ey10);
	    fields[n].setName("Ex01",ex01);
	    fields[n].setName("Ey01",ey01);
	    fields[n].setName("Hz",hz11);
	  }
	  fields[n].setOperators(*ops);
	}
	else
	{
	  //kkc 040304 XXX  BILL, I assumed the operators were set for these guys since you
	  //            use them elsewhere 
	  int c=0;
	  if( solveForElectricField )
	  {
	    ex=c; fields[n].setName("Ex",ex); c++;
	    ey=c; fields[n].setName("Ey",ey); c++;
	    ez=c; fields[n].setName("Ez",ez); c++;
	  }
	  else
	  {
	    ex=c;  //  ex always points to first component
	  }
	  if( solveForMagneticField )
	  {
	    hx=c; fields[n].setName("Hx",hx); c++;
	    hy=c; fields[n].setName("Hy",hy); c++;
	    hz=c; fields[n].setName("Hz",hz); c++;
	  }
	  else
	  {
	    hz=c-1;   // hz always points to the last component
	  }
	}
      }

    }
    else if( method==sosup )
    {      
      // SOS Upwind: we store the fields and the time derivatives 
      fields = new realMappedGridFunction [numberOfTimeLevels];
      for( int n=0; n<numberOfTimeLevels; n++ )
      {
	fields[n].updateToMatchGrid(mg,all,all,all,numberOfComponents);
	fields[n]=0.;
	      
	if( numberOfDimensions==2 )
	{
          int c=0; 
	  ex=c;  fields[n].setName("Ex",ex); c++;
	  ey=c;  fields[n].setName("Ey",ey); c++;
	  hz=c;  fields[n].setName("Hz",hz); c++;
          // time-derivatives:
	  ext=c; fields[n].setName("Ext",ext); c++;
	  eyt=c; fields[n].setName("Eyt",eyt); c++;
	  hzt=c; fields[n].setName("Hzt",hzt); c++;

	  fields[n].setOperators(*ops);
	}
	else
	{
	  int c=0;
	  if( solveForElectricField )
	  {
	    ex=c;  fields[n].setName("Ex",ex);   c++;
	    ey=c;  fields[n].setName("Ey",ey);   c++;
	    ez=c;  fields[n].setName("Ez",ez);   c++;
            // time-derivatives:
	    ext=c; fields[n].setName("Ext",ext); c++;
	    eyt=c; fields[n].setName("Eyt",eyt); c++;
	    ezt=c; fields[n].setName("Ezt",ezt); c++;
	  }
	  else
	  {
	    ex=c;  //  ex always points to first component
	  }
	  if( solveForMagneticField )
	  {
	    hx=c;  fields[n].setName("Hx",hx);   c++;
	    hy=c;  fields[n].setName("Hy",hy);   c++;
	    hz=c;  fields[n].setName("Hz",hz);   c++;
            // time-derivatives:
	    hxt=c; fields[n].setName("Hxt",hxt); c++;
	    hyt=c; fields[n].setName("Hyt",hyt); c++;
	    hzt=c; fields[n].setName("Hzt",hzt); c++;
	  }
	  else
	  {
	    hz=c-1;   // hz always points to the last component
	  }
	}
      }

    }
    else
    {
      // yee and dsi gridfunction setup
      //	  fields = new realMappedGridFunction [3*numberOfTimeLevels];
      //kkc 040304 dsi and yee schemes now use faceCenteredAll gridfunctions
      fields = new realMappedGridFunction [2*numberOfTimeLevels];
      if ( numberOfDimensions==2 )
      {
	ex = 0;
	ey = ez = 1;
	hz = 0;
      }
      else
      {
	ex=hx=0;
	ey=hy=1;
	ez=hz=2;

      }

      MappedGridOperators *ops = new MappedGridOperators(mg);
      for( int n=0; n<numberOfTimeLevels; n++ )
      {
	if ( numberOfDimensions==2 )
	{
	  fields[n].updateToMatchGrid(mg,GridFunctionParameters::cellCentered);
	  if ( method==dsiMatVec )
	  {
	    //		      fields[n+numberOfTimeLevels].updateToMatchGrid(mg,GridFunctionParameters::faceCenteredAll);
	    //		      fields[n+numberOfTimeLevels].updateToMatchGrid(mg,GridFunctionParameters::edgeCentered,2);
	    fields[n+numberOfTimeLevels].updateToMatchGrid(mg,GridFunctionParameters::edgeCentered);
	    fields[n+numberOfTimeLevels].setName("E.n",ex);
	  }
	  else
	  {
	    fields[n+numberOfTimeLevels].updateToMatchGrid(mg,GridFunctionParameters::faceCenteredAll,mg.numberOfDimensions()); //XXX the last argument should go away when the dsi schemes are fully implemented
	    fields[n+numberOfTimeLevels].setName("Ex",ex);
	    fields[n+numberOfTimeLevels].setName("Ey",ey);
	  }
	  fields[n].setName("Hz",hz);
	}
	else
	{
	  if ( method!=dsiMatVec )
	  {
	    fields[n].updateToMatchGrid(mg,GridFunctionParameters::faceCenteredAll,numberOfDimensions);
	    fields[n+numberOfTimeLevels].updateToMatchGrid(mg,GridFunctionParameters::edgeCentered,numberOfDimensions); 		   
	    fields[n].setName("Hx",hx);
	    fields[n].setName("Hy",hy);
	    fields[n].setName("Hz",hz);
	    fields[n+numberOfTimeLevels].setName("Ex",ex);
	    fields[n+numberOfTimeLevels].setName("Ey",ey);
	    fields[n+numberOfTimeLevels].setName("Ez",ez);

	  }
	  else
	  {
	    fields[n].updateToMatchGrid(mg,GridFunctionParameters::faceCenteredAll);
	    fields[n+numberOfTimeLevels].updateToMatchGrid(mg,GridFunctionParameters::edgeCentered); 
	    fields[n].setName("H.n",0);
	    fields[n+numberOfTimeLevels].setName("E.n",0);
	  }

	}
	//	      fields[n+numberOfTimeLevels].updateToMatchGrid(mg,GridFunctionParameters::faceCenteredAxis1,all,all,all,2);
	//	      fields[n+2*numberOfTimeLevels].updateToMatchGrid(mg,GridFunctionParameters::faceCenteredAxis2,all,all,all,2);

	//	      fields[n].setName("Hz",hz);
	//	      fields[n+numberOfTimeLevels].setName("E{.n,x}",ex);
	//	      fields[n+numberOfTimeLevels].setName("Ey",ey);
// 	      fields[n+numberOfTimeLevels].setName("Ex{10}",ey);
// 	      fields[n+2*numberOfTimeLevels].setName("Ex01",ex);
// 	      fields[n+2*numberOfTimeLevels].setName("Ex01",ey);
	fields[n].setOperators(*ops);
	fields[n+numberOfTimeLevels].setOperators(*ops);
	//	      fields[n+2*numberOfTimeLevels].setOperators(*ops);
      }

    }
  
    if( artificialDissipation>0. || artificialDissipationCurvilinear>0. )
    {
      dissipation=new realMappedGridFunction;
      if ( method==nfdtd || method==bamx || numberOfDimensions==2 )
	dissipation->updateToMatchGrid(mg,GridFunctionParameters::cellCentered);
      else
      {
	if ( method==dsiMatVec )
	  dissipation->updateToMatchGrid(mg,GridFunctionParameters::faceCenteredAll); 
	else
	  dissipation->updateToMatchGrid(mg,GridFunctionParameters::faceCenteredAll,numberOfDimensions); 
      }
	
      *dissipation = 0;
    }
  
    if( method==nfdtd || method==sosup || method==bamx )
    {
      op = new MappedGridOperators(mg);
      op->useConservativeApproximations(useConservative);
      op->setOrderOfAccuracy(orderOfAccuracyInSpace);
    }

  }
  else
  {
    // *********************************************
    // *********** CompositeGrid *******************
    // *********************************************

    if ( method==nfdtd || method==sosup || method==yee || method==bamx  ) 
    {
      Range all;
      if( numberOfDimensions==2 )
      {
	numberOfComponents= method==yee || method==nfdtd || method==sosup || method==bamx  ? 
	  (int)numberOfComponentsRectangularGrid : 
	  (int)numberOfComponentsCurvilinearGrid; 
	    
      }
      else
      {
	numberOfComponents=(int(solveForElectricField)+int(solveForMagneticField))*3;
      }

      if( method==bamx && solveForAllFields )
      {
        // Solve for all components of E and H 
        numberOfComponents=6;
      }
      


      if( method==sosup )
      {
	// With sosup we also advance/store the time-derivative of the fields
	numberOfComponents*= 2;
      }


      const int & useSosupDissipation = parameters.dbase.get<int>("useSosupDissipation");
      const int dw = max(cg[0].discretizationWidth()); // discertization width

      int & extrapolateInterpolationNeighbours = dbase.get<int>("extrapolateInterpolationNeighbours");
      if( dw > orderOfAccuracyInSpace+1 )
        extrapolateInterpolationNeighbours=false; // *wdh* added this check, June 15, 2016
      else if( useSosupDissipation ||
	      ( method==bamx && artificialDissipation>0. ) )
      {
        extrapolateInterpolationNeighbours=true;
        // Sosup dissipation requires extra ghost points *wdh* June 18, 2018

        // what order should this be? 
        dbase.get<int>("orderOfExtrapolationForInterpolationNeighbours")=orderOfAccuracyInSpace+1; 
        // dbase.get<int>("orderOfExtrapolationForInterpolationNeighbours")=orderOfAccuracyInSpace+2; // try this 

        // We need to increase the maximum allowable width to extrap interp neighbours
        GenericMappedGridOperators::setDefaultMaximumWidthForExtrapolateInterpolationNeighbours(
          dbase.get<int>("orderOfExtrapolationForInterpolationNeighbours")+1);
      }


      numberOfComponentsForTZ=numberOfComponents;
	
      // dispersionModelGridFunction[domain][numTimeLevels] : 
      realCompositeGridFunction **& dmgf = 
                   parameters.dbase.get<realCompositeGridFunction**>("dispersionModelGridFunction");

      // nonlinearModelGridFunction[domain][numTimeLevels] : 
      realCompositeGridFunction **& nmgf = parameters.dbase.get<realCompositeGridFunction**>("nonlinearModelGridFunction");


      const int & numberOfDomains = cg.numberOfDomains();
      const int numberOfComponentGrids = cg.numberOfComponentGrids();
      
      int & maxNumberOfPolarizationVectors = parameters.dbase.get<int>("maxNumberOfPolarizationVectors");
      maxNumberOfPolarizationVectors=0;
      
      // max number of extries on a grid (for TZ) for bamx
      int & maxNumberOfPolarizationComponents = parameters.dbase.get<int>("maxNumberOfPolarizationComponents");
      maxNumberOfPolarizationComponents=0;


      // total number of polarization components per grid 
      IntegerArray & totalNumberOfPolarizationComponents =
	parameters.dbase.get<IntegerArray>("totalNumberOfPolarizationComponents");
      totalNumberOfPolarizationComponents.redim(cg.numberOfComponentGrids());
      totalNumberOfPolarizationComponents=0;

      IntegerArray & numberOfMaterials = parameters.dbase.get<IntegerArray>("numberOfMaterials");
      std::vector<IntegerArray> & materialList = parameters.dbase.get<std::vector<IntegerArray> >("materialList");  

      if( method==bamx && numberOfMaterials.getLength(0) != cg.numberOfComponentGrids() )
      {
        // There must be a single material for BAMX (these arrays are normally set in defineMaterialRegions) 
        assert( numberOfMaterialRegions==1 );
        numberOfMaterials.redim(cg.numberOfComponentGrids());
        numberOfMaterials=0;
      }

      if( dispersionModel!=noDispersion )
      {
	// -- add Grid functions to hold the data for the dispersion models ---
        //    We only create grid functions on domains that are dispersive
        assert( dmgf==NULL );
        assert( numberOfDomains>=1 );
        
        // dispersionModelGridFunction[domain][numTimeLevels] : 
        dmgf = new realCompositeGridFunction* [numberOfDomains];

        cg.update(GridCollection::THEdomain); // create CG's for each domain
        
        maxNumberOfPolarizationVectors=0; // keep track of the max number of polarization vectors for TZ

        int & maxNumberOfNonlinearVectors = parameters.dbase.get<int>("maxNumberOfNonlinearVectors");
        maxNumberOfNonlinearVectors=0;    // keep track for TZ 

	if(  !( nonlinearModel == noNonlinearModel)  )
	{
          // --- nonlinear model ----
	  // Allocate pointers to the composite grid functions on different domains 
	  assert( nmgf==NULL );
	  nmgf = new realCompositeGridFunction* [numberOfDomains];

	}
	
        // ------------------------- START DOMAIN --------------------------------
        for( int domain=0; domain<cg.numberOfDomains(); domain++ ) 
        {
          dmgf[domain]=NULL;  // default 
          
	  if( method==nfdtd )
	  {
	    const DispersiveMaterialParameters & dmp = getDomainDispersiveMaterialParameters(domain);
	    const int numberOfPolarizationVectors = dmp.numberOfPolarizationVectors;
	    maxNumberOfPolarizationVectors=max(maxNumberOfPolarizationVectors,numberOfPolarizationVectors);

	    if( true )
              printF("setupGF: domain=%d: numberOfPolarizationVectors=%d\n",domain,numberOfPolarizationVectors);
	    

	    if( numberOfPolarizationVectors>0 )
	    {
	      dmgf[domain] = new realCompositeGridFunction [numberOfTimeLevels];

	      CompositeGrid & cgd = cg.domain[domain]; // Here is the CompositeGrid for just this domain
	      for( int n=0; n<numberOfTimeLevels; n++ )
	      {
              
		realCompositeGridFunction & u = dmgf[domain][n];
		const int numComp = numberOfPolarizationVectors*numberOfDimensions;
		u.updateToMatchGrid(cgd,all,all,all,numComp);
		u=0.;
		for( int iv=0; iv<numberOfPolarizationVectors; iv++ )
		{
		  int pc = iv*numberOfDimensions;
		  u.setName(sPrintF("Px%i",iv),pc);   pc++;
		  u.setName(sPrintF("Py%i",iv),pc);   pc++;
		  if( numberOfDimensions==3 )
		  {
		    u.setName(sPrintF("Pz%i",iv),pc);  pc++;
		  }
		}
              

	      }
            
	    }

	    pxc=numberOfComponents;
	    pyc=pxc+1;
	    if( numberOfDimensions==3 ) pzc=pyc+1;
        
            if( dmp.isNonlinearMaterial() )
	    {
	      // -------- nonlinear material -----
              // allocate grid functions to hold N's 
	      assert( nmgf != NULL );

              const int numberOfAtomicLevels = dmp.getNumberOfAtomicLevels();
	      printF("MX:setupGridFunctions: numberOfAtomicLevels=%d\n",numberOfAtomicLevels);

	      maxNumberOfNonlinearVectors = max(maxNumberOfNonlinearVectors,numberOfAtomicLevels);
	       
	      if( numberOfAtomicLevels>0 )
	      {
		nmgf[domain] = new realCompositeGridFunction [numberOfTimeLevels];

		CompositeGrid & cgd = cg.domain[domain]; // Here is the CompositeGrid for just this domain
		for( int n=0; n<numberOfTimeLevels; n++ )
		{
              
		  realCompositeGridFunction & unl = nmgf[domain][n];
		  const int numComp = numberOfAtomicLevels;
		  unl.updateToMatchGrid(cgd,all,all,all,numComp);
		  unl=0.;
		  for( int inl=0; inl<numComp; inl++ )
		  {
		    unl.setName(sPrintF("N%i",inl),inl);  
		  }
		}
	      }
	    }

	  }
	  else if( method==bamx )
	  {
            // BA GDM equations
            assert( domain==0 );


	    if( !dbase.has_key("materialRegionParameters") )
	    {
	      std::vector<DispersiveMaterialParameters> & dmpVector =
		dbase.put<std::vector<DispersiveMaterialParameters> >("materialRegionParameters");

	      // Material "0" is the background material: 
              int grid=0;
	      const int domain = cg.domainNumber(grid);
	      const DispersiveMaterialParameters & dmp0 = getDomainDispersiveMaterialParameters(domain);
	      dmpVector.push_back(dmp0);
	      // aString label="dmpVector[0]";
	      // dmpVector[0].display(stdout,label);

     
	    }

	    std::vector<DispersiveMaterialParameters> & dmpVector =
              dbase.get<std::vector<DispersiveMaterialParameters> >("materialRegionParameters");

	    // const DispersiveMaterialParameters & dmp = getDomainDispersiveMaterialParameters(domain);
	    int numberOfGDMTerms=0;
	    for( int mr=0; mr<numberOfMaterialRegions; mr++ )
	    {
	      const DispersiveMaterialParameters & dmp = dmpVector[mr];
	      const IntegerArray & Np = dmp.getBianisotropicNp();
	      numberOfGDMTerms += sum(Np);
	    }
	    
            // *****************fix me for more that one grid 
	    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
	      totalNumberOfPolarizationComponents(grid)=numberOfGDMTerms;

	    maxNumberOfPolarizationComponents=numberOfGDMTerms;

	    printF("setupGF: BA-GDM: numberOfMaterialRegions=%d: numberOfGDMTerms=%d (total)\n",
		   numberOfMaterialRegions,numberOfGDMTerms);
	    

	    if( numberOfGDMTerms>0 )
	    {
	      dmgf[domain] = new realCompositeGridFunction [numberOfTimeLevels];

	      CompositeGrid & cgd = cg.domain[domain]; // Here is the CompositeGrid for just this domain
	      for( int n=0; n<numberOfTimeLevels; n++ )
	      {
              
		realCompositeGridFunction & u = dmgf[domain][n];
		const int numComp = numberOfGDMTerms*2;  // we store P and Q=Pt 
		u.updateToMatchGrid(cgd,all,all,all,numComp);
		u=0.;
		for( int iv=0; iv<numberOfGDMTerms; iv++ )
		{
		  int pc = iv*2;
		  u.setName(sPrintF("P%i",iv),pc);   pc++;
		  u.setName(sPrintF("Q%i",iv),pc);   pc++;
		}
              

	      }
            
	    }

            // numberOfComponentsForTZ += maxNumberOfPolarizationComponents*2;  // July 17, 2020 -- move below 

	  }
	  
        } // ----- end for domain ------
	

	if(  method!=bamx )
	{
	  numberOfComponentsForTZ += numberOfDimensions*maxNumberOfPolarizationVectors; // *** check me *** why inside the loop ?
	}
	else
	{
	  numberOfComponentsForTZ += maxNumberOfPolarizationComponents*2;
	}
	

	numberOfComponentsForTZ += maxNumberOfNonlinearVectors;
	

	// *** CHECK numberOfComponentsForTZ  July 15, 2020 **FINISH ME***
	printf("\n +++++++ setupGridFunctions: numberOfDomains=%d, numberOfComponents=%d maxNumberOfPolarizationVectors=%d, "
	       " maxNumberOfNonlinearVectors=%d numberOfComponentsForTZ=%d\n",
	       numberOfDomains,numberOfComponents,maxNumberOfPolarizationVectors,maxNumberOfNonlinearVectors,numberOfComponentsForTZ);
	

        
        // We need the mapping from global "grid" number to the domain grid number
        //       cg[grid] <-> cg.domain[domainGridNumber(grid)]
        IntegerArray & domainGridNumber = parameters.dbase.put<IntegerArray>("domainGridNumber");
        domainGridNumber.redim(numberOfComponentGrids);
        IntegerArray nG(numberOfDomains); nG=0;  // nG = counts grids per domain 
        for( int grid=0; grid<numberOfComponentGrids; grid++ ) 
        {
          int d = cg.domainNumber(grid);      // domain number 
          domainGridNumber(grid) = nG(d);     //  gc[grid] <-> gc.domain[domainGridNumber(grid)]
          nG(d)++;
        }
      }
      


      printF("\n >>>>>>>>>>> setupGridFunctions: numberOfDomains=%i\n",numberOfDomains);
      if( numberOfDomains>1 )
      {
        for( int domain=0; domain<cg.numberOfDomains(); domain++ )
	{
  	  const DispersiveMaterialParameters & dmp = getDomainDispersiveMaterialParameters(domain);
          const int numberOfPolarizationVectors = dmp.numberOfPolarizationVectors;

          printF(" Domain %i (%s) number of polarization vectors=%i\n",
                  domain,(const char*)cg.getDomainName(domain),numberOfPolarizationVectors);
	}

      }
      

      if( true || method==yee || method==sosup )
	printF("\n *********** setupGridFunctions: numberOfComponents=%i numberOfTimeLevels=%i dispersionModel=%i\n\n",
               numberOfComponents, numberOfTimeLevels,(int)dispersionModel);

      numberOfFields=2;
	
      // *new* Store grid functions in array gf[] -- to transition to a DomainSolver 
      gf = new GridFunction [numberOfTimeLevels]; //  *new* way 2017/05/18

      cgfields = new realCompositeGridFunction [numberOfTimeLevels]; // old way
	
      for( int n=0; n<numberOfTimeLevels; n++ )
      {
        cgfields[n].updateToMatchGrid(cg,all,all,all,numberOfComponents);

        // do this for now: (keep cgfields array for now)
        gf[n].cg.reference(cg);
        gf[n].u.reference(cgfields[n]);
        gf[n].setParameters(parameters);

        // else
        // {
        //   // ---- Allocate grid functions for the dispersive model ----
        //   // Different domains will have possibly different numbers of components since 
        //   //    each domain can have a different number of polarization vectors 

        //   // first allocate with minimal number of components
        //   cgfields[n].updateToMatchGrid(cg,all,all,all,numberOfComponents);

 	//   for( int domain=0; domain<cg.numberOfDomains(); domain++ )
        //   {
        //     const DispersiveMaterialParameters & dmp = getDomainDispersiveMaterialParameters(domain);
        //     const int numberOfPolarizationVectors = dmp.numberOfPolarizationVectors;

        //     if( numberOfPolarizationVectors==0 )
        //     {
        //       printF(" Domain %i (%s) number of polarization vectors=%i, numberOfComponents=%i\n",
        //              domain,(const char*)cg.getDomainName(domain),numberOfPolarizationVectors,numberOfComponents);

        //     }
        //     else if( numberOfPolarizationVectors>0 )
        //     {
        //       const int numComp = numberOfComponents + numberOfDimensions*numberOfPolarizationVectors;
        //       printF(" Domain %i (%s) number of polarization vectors=%i, --> numComp=%i\n",
        //              domain,(const char*)cg.getDomainName(domain),numberOfPolarizationVectors,numComp);

        //       for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        //       {
        //         if( cg.domainNumber(grid)==domain )
        //         {
        //           realMappedGridFunction & ud = cgfields[n][grid];
        //           ud.updateToMatchGrid(cg[grid],all,all,all,numComp);
        //         }
        //       }
        //     }
            
        //     // CompositeGrid & cgd = cg.domain[domain];
        //     // realCompositeGridFunction & ud = cgfields[n].domain[domain];
        //     // ud.updateToMatchGrid(cgd,all,all,all,numberOfComponents);
        //   }
          
        // }
        
	cgfields[n]=0.;
	    
        int c=0; // counts components
	if( method!=sosup )
	{
          //  --- set field component names for nfdtd or Yee ---
          assert( method==nfdtd || method==yee || method==bamx );
	  
          if( numberOfComponents==6 )
          {
            ex=c; cgfields[n].setName("Ex",ex); c++;
            ey=c; cgfields[n].setName("Ey",ey); c++;
            ez=c; cgfields[n].setName("Ez",ez); c++;
            hx=c; cgfields[n].setName("Hx",hx); c++;
            hy=c; cgfields[n].setName("Hy",hy); c++;
            hz=c; cgfields[n].setName("Hz",hz); c++;
          }
          else if( numberOfDimensions==2 )
	  {
	    if(  method!=dsiMatVec )
	    {
	      cgfields[n].setName("Ex",ex);
	      cgfields[n].setName("Ey",ey);
	      cgfields[n].setName("Hz",hz);
	      c+=3;
	    }
	    else
	    {
	      cgfields[n].setName("Ex10",ex10);
	      cgfields[n].setName("Ey10",ey10);
	      cgfields[n].setName("Ex01",ex01);
	      cgfields[n].setName("Ey01",ey01);
	      cgfields[n].setName("Hz",hz11);
	      c+=5;
	    }
	  }
	  else
	  {
	    if( solveForElectricField )
	    {
	      ex=c; cgfields[n].setName("Ex",ex); c++;
	      ey=c; cgfields[n].setName("Ey",ey); c++;
	      ez=c; cgfields[n].setName("Ez",ez); c++;
	    }
	    if( solveForMagneticField )
	    {
	      hx=c; cgfields[n].setName("Hx",hx); c++;
	      hy=c; cgfields[n].setName("Hy",hy); c++;
	      hz=c; cgfields[n].setName("Hz",hz); c++;
	    }
	    
	  }
	}
	else 
	{
	  // --- set field component names for SOSUP ---
	  assert( method==sosup );
	  if( numberOfDimensions==2 )
	  {
	    ex=c;  cgfields[n].setName("Ex",ex); c++;
	    ey=c;  cgfields[n].setName("Ey",ey); c++;
	    hz=c;  cgfields[n].setName("Hz",hz); c++;
	    // time-derivatives:
	    ext=c; cgfields[n].setName("Ext",ext); c++;
	    eyt=c; cgfields[n].setName("Eyt",eyt); c++;
	    hzt=c; cgfields[n].setName("Hzt",hzt); c++;
	  }
	  else
	  {
	    if( solveForElectricField )
	    {
	      ex=c;  cgfields[n].setName("Ex",ex);   c++;
	      ey=c;  cgfields[n].setName("Ey",ey);   c++;
	      ez=c;  cgfields[n].setName("Ez",ez);   c++;
	      // time-derivatives:
	      ext=c; cgfields[n].setName("Ext",ext); c++;
	      eyt=c; cgfields[n].setName("Eyt",eyt); c++;
	      ezt=c; cgfields[n].setName("Ezt",ezt); c++;
	    }
	    else
	    {
	      ex=c;  //  ex always points to first component
	    }
	    if( solveForMagneticField )
	    {
	      hx=c;  cgfields[n].setName("Hx",hx);   c++;
	      hy=c;  cgfields[n].setName("Hy",hy);   c++;
	      hz=c;  cgfields[n].setName("Hz",hz);   c++;
	      // time-derivatives:
	      hxt=c; cgfields[n].setName("Hxt",hxt); c++;
	      hyt=c; cgfields[n].setName("Hyt",hyt); c++;
	      hzt=c; cgfields[n].setName("Hzt",hzt); c++;
	    }
	    else
	    {
	      hz=c-1;   // hz always points to the last component
	    }
	  }
	}

	
      } // end for n
      
	
      // *wdh* Do this as needed in computeDissipation
//  	if( artificialDissipation>0. )
//  	{
//  	  cgdissipation=new realCompositeGridFunction;
//  	  cgdissipation->updateToMatchGrid(cg,all,all,all,numberOfComponents);
//  	  if( numberOfComponents==3 )
//  	  {
//  	    cgdissipation->setName("Ex dissipation",ex);
//  	    cgdissipation->setName("Ey dissipation",ey);
//  	    cgdissipation->setName("Hz dissipation",hz);
//  	  }
//  	}
	
      if( method==nfdtd || method==sosup || method==bamx )
      {
	cgop = new CompositeGridOperators(cg);
	cgop->useConservativeApproximations(useConservative);
	cgop->setOrderOfAccuracy(orderOfAccuracyInSpace);
      }
    }
    else // dsi scheme
    {
      Range all;
      if ( numberOfDimensions==2 )
      {
	ex = 0;
	ey = ez = 1;
	hz = 0;
      }
      else
      {
	ex=hx=0;
	ey=hy=1;
	ez=hz=2;

      }

      CompositeGridOperators *ops = new CompositeGridOperators(cg);
      dsi_cgfieldsE0 = new realCompositeGridFunction[numberOfTimeLevels];
      //		  dsi_cgfieldsE1 = new realCompositeGridFunction[numberOfTimeLevels];
      dsi_cgfieldsH  = new realCompositeGridFunction[numberOfTimeLevels];
      for( int n=0; n<numberOfTimeLevels; n++ )
      {
	if ( numberOfDimensions==2 )
	{
	  dsi_cgfieldsH[n].updateToMatchGrid(cg,GridFunctionParameters::cellCentered);
	  if ( method==dsiMatVec )
	  {
	    //		      dsi_cgfieldsE0[n].updateToMatchGrid(cg,GridFunctionParameters::edgeCentered,2);
	    dsi_cgfieldsE0[n].updateToMatchGrid(cg,GridFunctionParameters::edgeCentered);
	    dsi_cgfieldsE0[n].setName("E.n",ex);
	    if( artificialDissipation>0. )
	    {
	      cgdissipation=new realCompositeGridFunction;
	      cgdissipation->updateToMatchGrid(cg,GridFunctionParameters::cellCentered);
	      cgdissipation->setName("H dissp",0);
	      e_cgdissipation=new realCompositeGridFunction;
	      e_cgdissipation->updateToMatchGrid(cg,GridFunctionParameters::edgeCentered);
	      e_cgdissipation->setName("E.n dissp",0);
	    }
	  }
	  else
	  {
	    dsi_cgfieldsE0[n].updateToMatchGrid(cg,GridFunctionParameters::faceCenteredAll,cg.numberOfDimensions()); //XXX the last argument should go away when the dsi schemes are fully implemented
	    dsi_cgfieldsE0[n].setName("Ex",ex);
	    dsi_cgfieldsE0[n].setName("Ey",ey);
	  }
	  dsi_cgfieldsH[n].setName("Hz",hz);
	}
	else
	{

	  if ( method!=dsiMatVec )
	  {
	    dsi_cgfieldsH[n].updateToMatchGrid(cg,GridFunctionParameters::faceCenteredAll,cg.numberOfDimensions());
	    dsi_cgfieldsE0[n].updateToMatchGrid(cg,GridFunctionParameters::edgeCentered,cg.numberOfDimensions()); 		   
	    dsi_cgfieldsH[n].setName("Hx",hx);
	    dsi_cgfieldsH[n].setName("Hy",hy);
	    dsi_cgfieldsH[n].setName("Hz",hz);
	    dsi_cgfieldsE0[n].setName("Ex",ex);
	    dsi_cgfieldsE0[n].setName("Ey",ey);
	    dsi_cgfieldsE0[n].setName("Ez",ez);

	  }
	  else
	  {
	    dsi_cgfieldsH[n].updateToMatchGrid(cg,GridFunctionParameters::faceCenteredAll);
	    dsi_cgfieldsE0[n].updateToMatchGrid(cg,GridFunctionParameters::edgeCentered); 
	    dsi_cgfieldsH[n].setName("H.n",0);
	    dsi_cgfieldsE0[n].setName("E.n",0);
	    if( artificialDissipation>0. && !cgdissipation )
	    {
	      cgdissipation=new realCompositeGridFunction;
	      cgdissipation->updateToMatchGrid(cg,GridFunctionParameters::faceCenteredAll);
	      cgdissipation->setName("H dissp",0);

	      e_cgdissipation=new realCompositeGridFunction;
	      e_cgdissipation->updateToMatchGrid(cg,GridFunctionParameters::edgeCentered);
	      e_cgdissipation->setName("E dissp",0);
	    }

	  }

	}
	//	      fields[n+numberOfTimeLevels].updateToMatchGrid(mg,GridFunctionParameters::faceCenteredAxis1,all,all,all,2);
	//	      fields[n+2*numberOfTimeLevels].updateToMatchGrid(mg,GridFunctionParameters::faceCenteredAxis2,all,all,all,2);

	//	      fields[n].setName("Hz",hz);
	//	      fields[n+numberOfTimeLevels].setName("E{.n,x}",ex);
	//	      fields[n+numberOfTimeLevels].setName("Ey",ey);
// 	      fields[n+numberOfTimeLevels].setName("Ex{10}",ey);
// 	      fields[n+2*numberOfTimeLevels].setName("Ex01",ex);
// 	      fields[n+2*numberOfTimeLevels].setName("Ex01",ey);
	dsi_cgfieldsH[n].setOperators(*ops);
	dsi_cgfieldsE0[n].setOperators(*ops);
	//	      fields[n+2*numberOfTimeLevels].setOperators(*ops);
      }

#if 0
      ex = 0;
      ey = 1;
      hz = 0;

      for( int n=0; n<numberOfTimeLevels; n++ )
      {
	dsi_cgfieldsE0[n].updateToMatchGrid(cg,GridFunctionParameters::faceCenteredAll,all,all,all,Range(2));
	dsi_cgfieldsE1[n].updateToMatchGrid(cg,GridFunctionParameters::faceCenteredAll,all,all,all,Range(2));
	dsi_cgfieldsH[n].updateToMatchGrid(cg,GridFunctionParameters::cellCentered,all,all,all);

	dsi_cgfieldsE0[n].setName("Ex{10}",ex);
	dsi_cgfieldsE0[n].setName("Ey{10}",ey);
	dsi_cgfieldsE1[n].setName("Ex01",ex);
	dsi_cgfieldsE1[n].setName("Ey01",ey);
	dsi_cgfieldsH[n].setName("Hz",hz);
      }
#endif
    }
  }
  
  if( method==yee )
  {
    methodName="Yee";
  }
  else if( method==dsi )
  {
    methodName="DSI";
  }
  else if ( method==dsiMatVec)
  {
    methodName="DSI-MatVec";
  }
  else if( method==nfdtd )
  {
    methodName="FD";
  }
  else if( method==dsiNew )
  {
    methodName="DSI-new";
  }
  else if( method==sosup )
  {
    methodName="SOSUP";
  }
  else if( method==bamx )
  {
    methodName="BAMX";
  }
  else
  {
    printF("Maxwell:setup:ERROR: unknown method=%i\n",(int)method);
    OV_ABORT("ERROR");
  }

  dxMinMax.redim(cgp->numberOfComponentGrids(),2);
  dxMinMax=0.;
  
  if( useVariableDissipation && method==nfdtd )
    buildVariableDissipation();

  if ( method==dsiMatVec )
  {
    real tm0 = getCPU();
    setupDSICoefficients();
    timing(timeForDSIMatrix) += getCPU()-tm0;
  }

  if( initialConditionOption==annulusEigenfunctionInitialCondition )
  {
    // try to guess the cylinder radius and length

    if( cg[0].isRectangular() )
    {
      // assume grid 0 is the inner core 
      MappedGrid & mg = cg[0];
      real dx[3]={1.,1.,1.}, xab[2][3]={0.,0.,0.,0.,0.,0.};
      mg.getRectangularGridParameters( dx, xab );

      cylinderAxisStart=min(xab[0][2],xab[1][2]);
      cylinderAxisEnd  =max(xab[0][2],xab[1][2]);
 
      printF(" setupGridFunctions: I guess that the cylinder extends from [%8.2e,%8.2e] in the axial direction\n",
             cylinderAxisStart,cylinderAxisEnd    );
	
    }
      
  }

  if( useChargeDensity )
  {
    
    rc=numberOfComponentsForTZ;  // position of density in TZ functions

    if( useChargeDensity && pRho==NULL )
    {
      // allocate the grid function that holds the charge density
      pRho=new realCompositeGridFunction(cg);
    }
    numberOfComponentsForTZ += 1;
    
  }
  if( useTwilightZoneMaterials ) 
  {
    // for variable coefficients we define eps, mu, sigmaE and sigmaH with TZ functions:
    // ** FIX ME -- only do if we have variable coefficients ?? ***
    epsc   =rc+1;
    muc    =epsc+1;
    sigmaEc=muc+1;
    sigmaHc=sigmaEc+1;

    numberOfComponentsForTZ += 4;
  }
  
  numberOfSequences=numberOfComponents;
  if( computeEnergy )
    numberOfSequences+=2; // save the energy and delta(energy)

  computeTimeStep();

  // --- check for negative volumes : this is usually bad news --- *wdh* 2013/09/26
  const int numberOfGhost = orderOfAccuracyInSpace/2;
  int numberOfNegativeVolumes= GridStatistics::checkForNegativeVolumes( cg,numberOfGhost,stdout ); 
  if( numberOfNegativeVolumes>0 )
  {
    printF("Cgmx::FATAL Error: this grid has negative volumes (maybe only in ghost points).\n"
           "  This will normally cause severe or subtle errors. Please remake the grid.\n");
    OV_ABORT("ERROR");
  }
  else
  {
    printF("Cgmx:: No negative volumes were found.\n");
  }

  const int numberOfComponentGrids = cg.numberOfComponentGrids();

  const bool & useNewForcingMethod= dbase.get<bool>("useNewForcingMethod");
  if( useNewForcingMethod )
  {
    // --- Allocate space for forcing arrays ----
    //           (*new* wdh 2015/05/18)
    //   (1) Used to store external forcing at different time levels (e.g. for modified equation or SOSUP)
    //   (2) Used to store the RHS for method-of-lines schemes (.e.g Stoermer)


    int & numberOfForcingFunctions= dbase.get<int>("numberOfForcingFunctions"); // number of elements in forcingArray
    int & fCurrent = dbase.get<int>("fCurrent");          // forcingArray[fCurrent] : current forcing
    realArray *& forcingArray = dbase.get<realArray*>("forcingArray");
  
    numberOfForcingFunctions=0;
    if( (method==nfdtd || method==bamx ) && timeSteppingMethod==modifiedEquationTimeStepping )
    {
      // -- standard finite difference scheme + modified equation time-stepping

      if( forcingIsOn() )
      {
	// We need multiple time levels to compute time-differences in the modified equations:
	//    order=2 : 1 level needed
	//    order=4 : 3 levels for 2nd-order approx. to f_tt in ME forcing
	//    order=6 : 5 levels for 4th-order approx. to f_tt and 2nd-order to f_tttt
	numberOfForcingFunctions=orderOfAccuracyInTime-1;
      }
    }
    else if(  false &&  // finish me 
	      method==nfdtd && timeSteppingMethod== stoermerTimeStepping )
    {
      // -- standard finite difference scheme + Stoermer time-stepping

      // We need to store past time levels of the RHS
      numberOfForcingFunctions=orderOfAccuracyInTime-1; // *check me*
    }

    delete [] forcingArray;
    if( numberOfForcingFunctions>0 )
      forcingArray = new realArray [numberOfComponentGrids];
    fCurrent=0;

    // --- allocate space ---
    if( numberOfForcingFunctions>0 )
    {
      for( int grid=0; grid<numberOfComponentGrids; grid++ )
      {
	MappedGrid & mg = cg[grid];
	Index I1,I2,I3;
	getIndex(mg.dimension(),I1,I2,I3);
	Range C(ex,hz);
	
	realArray & fa = forcingArray[grid];
	fa.partition(mg.getPartition());
	fa.redim(I1,I2,I3,C,numberOfForcingFunctions);
      }
    }
    
  }
  
  // ----------------------------------------------------------
  // ---- Allocate work space for right-hand-sides etc. -------
  // ----------------------------------------------------------

  delete [] fn;

  bool allocateWorkSpaceForRectangularGrids=false;
  
  if( method==nfdtd || method==bamx )
  {
    // WORK-SPACE REQUIREMENTS

    if( timeSteppingMethod==modifiedEquationTimeStepping )
    {
      if( useConservative )
      {  // fn[0] : store uLapSq in advanceStructured
	numberOfFunctions=1;
      }
      else
      {
	// fn[0] : for ut1ptr -> advOpt 
         numberOfFunctions=1;
      }
      
    }
    else if( timeSteppingMethod==rungeKutta )
    {
      numberOfFunctions=4;  // IS THIS RIGHT ?  Forcing + k1,k2,k3,k4
      allocateWorkSpaceForRectangularGrids=true;
      
    }
    else
    {
      numberOfFunctions=3;  // *********************** fix *************** how many are needed?
    }
    
  }
  else if( method==sosup )
  {
    numberOfFunctions=3;  // *********************** fix *************** how many are needed?
  }
  
  currentFn=0; 
  if( numberOfFunctions>0 )
  {
    fn = new realArray [numberOfFunctions*numberOfComponentGrids];

    // *wdh* Jan 6, 2017 Range C(ex,hz);
    int numberOfComponents=cgfields[0][0].getLength(3);
    if( dispersionModel != noDispersion && method==bamx && timeSteppingMethod==rungeKutta )
    {
      // FOR BAMX + RK we store time derivatives of p,q in fn 
      const int & maxNumberOfPolarizationComponents = parameters.dbase.get<int>("maxNumberOfPolarizationComponents");
	
     // ** FIX ME FOR more than one grid --- NOT ALL GRIDS HAVE ALL POLARIZATION VECTORS 
      numberOfComponents += maxNumberOfPolarizationComponents*2; 
    }

    for( int grid=0; grid<numberOfComponentGrids; grid++ )
    {
      MappedGrid & mg = cg[grid];
      const bool isRectangular=mg.isRectangular();

      // **FIX ME** work space may not be needed on Cartesian grids *****************

      Index I1,I2,I3;
      getIndex(mg.dimension(),I1,I2,I3);

      Range C = numberOfComponents;   

      #define FN(m) fn[m+numberOfFunctions*(grid)]
      for( int m=0; m<numberOfFunctions; m++ )
      {
        // -- no need to allocate work space for rectangular grids and these methods
	if( allocateWorkSpaceForRectangularGrids || !isRectangular )
	{
	  FN(m).partition(mg.getPartition());
	  FN(m).redim(I1,I2,I3,C);
	}
      }
      
    }
    
  }
  

  if( method==nfdtd || method==bamx )
  {
    // -- generate geometry arrays ---
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
      MappedGrid & mg=cg[grid];
      const bool isRectangular=mg.isRectangular();
      if( !isRectangular )
      {
        if( useConservative )
        {
  	  // The conservative operators need the jacobian
	   mg.update( MappedGrid::THEinverseVertexDerivative | 
                     MappedGrid::THEinverseCenterDerivative | MappedGrid::THEcenterJacobian );
        }
        else
        {
    	  mg.update( MappedGrid::THEinverseVertexDerivative | 
                   MappedGrid::THEinverseCenterDerivative );
        }
      }
    }
  }

  const int & useSosupDissipation = parameters.dbase.get<int>("useSosupDissipation");
  bool addBAMXDissipation = method==bamx && artificialDissipation>0.;
  if( method==sosup || addBAMXDissipation  || ( (method==nfdtd ) && useSosupDissipation) )
  {
    // Sosup requires an extra ghost line  -- check there are enough

    const int minGhostNeeded = orderOfAccuracyInSpace/2 +1;
    Range Rx=cg.numberOfDimensions();
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
      MappedGrid & mg=cg[grid];
      const IntegerArray & numberOfGhostPoints = mg.numberOfGhostPoints();
      int numGhost = min(numberOfGhostPoints(Range(0,1),Rx));
      if( numGhost < minGhostNeeded )
      {
        printF("--MX-- setupGridFunctions: ERROR: the grid does not have enough ghost points for sosup dissipation.\n"
	"   orderOfAccuracy=%i requires at least %i ghost points.\n"
	"   You could remake the grid with more ghost points to fix this error.\n",
	orderOfAccuracyInSpace,minGhostNeeded);
	OV_ABORT("ERROR");
	
      } 
    }
  }
  
   
  // ---------------------- SETUP ERROR GRID FUNCTION ----------------
  Range all;
  if( mgp!=NULL )
  {

    MappedGrid & mg = *mgp;

    if( errp==NULL && checkErrors ) 
    {
      // Create a grid function to hold the errors for plotting
      if ( method==nfdtd || method==bamx )
      {
	int numberOfComponents = fields[0].getLength(3);
	errp = new realMappedGridFunction[1];
	errp->updateToMatchGrid(mg,all,all,all,numberOfComponents);
	errp->setName("Ex error",ex);
	errp->setName("Ey error",ey);
	errp->setName("Hz error",hz);
      }
      else
      {
	errp = new realMappedGridFunction[2];
	if ( cg.numberOfDimensions()==2 )
	{
	  errp[0].updateToMatchGrid(mg,GridFunctionParameters::edgeCentered,2);
	  errp[1].updateToMatchGrid(mg,GridFunctionParameters::cellCentered);
	  errp[0].setName("Ex error",ex);
	  errp[0].setName("Ey error",ey);
	  errp[1].setName("Hz error",hz);
	}
	else
	{
	  errp[0].updateToMatchGrid(mg,GridFunctionParameters::edgeCentered,3);
	  errp[0].setName("Ex error",ex);
	  errp[0].setName("Ey error",ey);
	  errp[0].setName("Ez error",ez);
	  errp[1].updateToMatchGrid(mg,GridFunctionParameters::faceCenteredAll,3);
	  errp[1].setName("Hx error",hx);
	  errp[1].setName("Hy error",hy);
	  errp[1].setName("Hz error",hz);
	}
      }
    }
  }
  else
  {
    if( checkErrors )//&& cg[0].getGridType()==MappedGrid::structuredGrid )
    {
      if( method==nfdtd || method==yee || method==bamx )
      {
	cgerrp = new realCompositeGridFunction [1];
	int numberOfComponents = cgfields[0][0].getLength(3);
	cgerrp->updateToMatchGrid(*cgp,all,all,all,numberOfComponents);
	
        if( solveForAllFields )
        {
          cgerrp->setName("Ex error",ex);
          cgerrp->setName("Ey error",ey);
          cgerrp->setName("Ez error",ez);
          cgerrp->setName("Hx error",hx);
          cgerrp->setName("Hy error",hy);
          cgerrp->setName("Hz error",hz);
        }
	else if( cg.numberOfDimensions()==2 )
	{
	  cgerrp->setName("Ex error",ex);
	  cgerrp->setName("Ey error",ey);
	  cgerrp->setName("Hz error",hz);
	  
	}
	else
	{
	  if( solveForElectricField )
	  {
	    cgerrp->setName("Ex error",ex);
	    cgerrp->setName("Ey error",ey);
	    cgerrp->setName("Ez error",ez);
	  }
	  if( solveForMagneticField )
	  {
	    cgerrp->setName("Hx error",hx);
	    cgerrp->setName("Hy error",hy);
	    cgerrp->setName("Hz error",hz);
	  }
	      
	}
      }
      else if( method==sosup )
      {
	cgerrp = new realCompositeGridFunction [1];
	int numberOfComponents = cgfields[0][0].getLength(3);
	cgerrp->updateToMatchGrid(*cgp,all,all,all,numberOfComponents);
	
	if( cg.numberOfDimensions()==2 )
	{
	  assert( numberOfComponents==6 );
	  cgerrp->setName("Ex error",ex);
	  cgerrp->setName("Ey error",ey);
	  cgerrp->setName("Hz  error",hz );
	  cgerrp->setName("Ext error",ext);
	  cgerrp->setName("Eyt error",eyt);
	  cgerrp->setName("Hzt error",hzt);
	}
	else
	{
	  if( solveForElectricField )
	  {
	    cgerrp->setName("Ex error",ex);
	    cgerrp->setName("Ey error",ey);
	    cgerrp->setName("Ez error",ez);
	    cgerrp->setName("Ext error",ext);
	    cgerrp->setName("Eyt error",eyt);
	    cgerrp->setName("Ezt error",ezt);
	  }
	  if( solveForMagneticField )
	  {
	    cgerrp->setName("Hx error",hx);
	    cgerrp->setName("Hy error",hy);
	    cgerrp->setName("Hz error",hz);
	    cgerrp->setName("Hxt error",hxt);
	    cgerrp->setName("Hyt error",hyt);
	    cgerrp->setName("Hzt error",hzt);
	  }
	      
	}
      }
      else
      {
	cgerrp = new realCompositeGridFunction [2];
	if ( cg.numberOfDimensions()==2 )
	{
	  cgerrp[0].updateToMatchGrid(*cgp,GridFunctionParameters::edgeCentered,2);
	  cgerrp[1].updateToMatchGrid(*cgp,GridFunctionParameters::cellCentered);
	  cgerrp[0].setName("Ex error",ex);
	  cgerrp[0].setName("Ey error",ey);
	  cgerrp[1].setName("Hz error",hz);	  
	}
	else
	{
	  cgerrp[0].updateToMatchGrid(*cgp,GridFunctionParameters::edgeCentered,3);
	  cgerrp[1].updateToMatchGrid(*cgp,GridFunctionParameters::faceCentered,3);

	  cgerrp[0].setName("Ex error",ex);
	  cgerrp[0].setName("Ey error",ey);
	  cgerrp[0].setName("Ez error",ez);
	      
	  cgerrp[1].setName("Hx error",hx);
	  cgerrp[1].setName("Hy error",hy);
	  cgerrp[1].setName("Hz error",hz);

	}
      }
    }
    
  }

  // Initialize super-grid absorbing layer functions if needed.
  buildSuperGrid();

  timing(timeForInitialize)+=getCPU()-time0;

  printF("\n <<<<<<<<<<<<<<<<<<<<<<<<< Maxwell: LEAVING setupGridFunctions\n");
  
  return 0;
}

