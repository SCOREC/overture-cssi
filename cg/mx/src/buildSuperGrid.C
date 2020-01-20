#include "Maxwell.h"
#include "ParallelUtility.h"
#include "PlotStuff.h"
#include "display.h"

// ================================================================================================
///  \brief  Adjust index bounds to account for the an absorbing BC (.e.g. supergrid)
//   when computing errors or plotting 
//
/// /param extra (input) : an additional offset 
/// /Iv (input/output) :
/// /Return value: true if the Index Iv was changed.
// ================================================================================================
bool Maxwell::
adjustBoundsForAbsorbingLayer( MappedGrid & mg, Index Iv[3], int extra /* =0 */ )
{
  bool useAbsorbing = (mg.boundaryCondition(0,0)==absorbing || mg.boundaryCondition(1,0)==absorbing ||
		       mg.boundaryCondition(0,1)==absorbing || mg.boundaryCondition(1,1)==absorbing ||
		       mg.boundaryCondition(0,2)==absorbing || mg.boundaryCondition(1,2)==absorbing);
  
  if( !useAbsorbing ) return false;

  const bool isRectangular=mg.isRectangular();
  assert( isRectangular );  // 
  real dx[3]={1.,1.,1.};
  mg.getDeltaX( dx );
  const IntegerArray & gid = mg.gridIndexRange();
      
  const real & superGridWidth = parameters.dbase.get<real>("superGridWidth");

  //printF(" +++ adjustBoundsForAbsorbingLayer: superGridWidth=%g, extra=%d\n",superGridWidth,extra);
  //printF(" Input: Iv=[%d,%d][%d,%d]\n",Iv[0].getBase(),Iv[0].getBound(),Iv[1].getBase(),Iv[1].getBound());
  for( int axis=0; axis<mg.numberOfDimensions(); axis++ )
  {
    int offset =  int( superGridWidth/dx[axis] + .5 ) + extra;  // what should this be ? 
    //printF("[axis=%d dx=%e offset= %d]",axis,dx[axis],offset);

    int na=Iv[axis].getBase();
    if( mg.boundaryCondition(0,axis)==absorbing )
    {
      na = max(na,gid(0,axis) + offset);
      na = min(na, int( gid(0,axis)+(gid(0,axis)+gid(1,axis)-1)/2  ));
    }
    
    int nb=Iv[axis].getBound();
    if( mg.boundaryCondition(1,axis)==absorbing )
    {
      nb = min(nb,gid(1,axis) - offset);
      nb = max(nb, int( gid(0,axis)+ (gid(0,axis)+gid(1,axis)+1)/2  ));
    }
    

    Iv[axis]=Range(na,nb);
  }
  //printF("\n");
  //printF(" Output: Iv=[%d,%d][%d,%d]\n",Iv[0].getBase(),Iv[0].getBound(),Iv[1].getBase(),Iv[1].getBound());

  
  return useAbsorbing;
}




// =================================================================================================
/// \brief Build the super-grid absorbing layer functions 
// =================================================================================================
int Maxwell::
buildSuperGrid( )
{
  int & useSuperGrid = parameters.dbase.get<int>("useSuperGrid");

  // useSuperGrid=true;  // *** TEST ****
  
  if( useSuperGrid )
  {
    printF("\n ++++++Build Supergrid absrobing layer functions +++++++\n");

    // Here is the CompositeGrid: 
    assert( cgp!=NULL );
    CompositeGrid & cg = *cgp;
    const int numberOfComponentGrids = cg.numberOfComponentGrids();
    const int numberOfDimensions = cg.numberOfDimensions();

    // --- Build the Super-grid layer functions ----

    if( !parameters.dbase.has_key("etaxSuperGrid") )
    {
      parameters.dbase.put<RealArray*>("etaxSuperGrid" )=NULL;  // ** delete me ** 
      parameters.dbase.put<RealArray*>("etaySuperGrid" )=NULL;
      parameters.dbase.put<RealArray*>("etazSuperGrid" )=NULL;
      
    }
    RealArray *& etaxSuperGrid = parameters.dbase.get<RealArray*>("etaxSuperGrid" );
    RealArray *& etaySuperGrid = parameters.dbase.get<RealArray*>("etaySuperGrid" );
    RealArray *& etazSuperGrid = parameters.dbase.get<RealArray*>("etazSuperGrid" );
    
    etaxSuperGrid = new RealArray[numberOfComponentGrids];
    etaySuperGrid = new RealArray[numberOfComponentGrids];
    etazSuperGrid = new RealArray[numberOfComponentGrids];

    real epsL,pSG,qSG;
    real ax,bx,ay,by;
    
#define sigmaSG(z) ( (1.-epsL)*pow( 1 - pow( (1 -(z)/superGridWidth),pSG),qSG) )
#define etaSG(z)  (1. -sigmaSG(z)) 

    epsL=1e-4;
    pSG=4;
    qSG=4;
    // superGridWidth=.2;
    const real & superGridWidth = parameters.dbase.get<real>("superGridWidth");
   
    // useAbsorbingLayer(axis,grid) = true or false if we use an absorbing layer for this axis and grid 
    IntegerArray & useAbsorbingLayer = parameters.dbase.get<IntegerArray>("useAbsorbingLayer");
    useAbsorbingLayer.redim(3,numberOfComponentGrids);
    useAbsorbingLayer=false;
    
    for( int grid=0; grid<numberOfComponentGrids; grid++ )
    {
      MappedGrid & mg = cg[grid];
      const bool isRectangular=mg.isRectangular();
      const IntegerArray & dimension = mg.dimension();
      const IntegerArray & gid = mg.gridIndexRange();
      const IntegerArray & bc = mg.boundaryCondition();

      bool buildLayersThisGrid=false;
      // bool buildLayersThisAxis[3] ={false,false,false }; // 
      for( int axis=0; axis<numberOfDimensions; axis++)
      {
        for( int side=0; side<=1; side++ )
        {
          if( bc(side,axis)==absorbing )   
          {
            buildLayersThisGrid=true;
            // buildLayersThisAxis[axis]=true;
	    useAbsorbingLayer(axis,grid)=true;
          }
        }
      }
      if( !buildLayersThisGrid ) continue;

      // // *** TEMP: For now build layers in all directions if any direction is needed -- need to optimize implementation
      // //  to allow only some directions 
      // 	for( int axis=0; axis<numberOfDimensions; axis++)
      // 	{
      // 	  buildLayersThisAxis[axis]=buildLayersThisGrid;
      // 	}
      

      assert( isRectangular );  // assume this for now 

      RealArray & etax = etaxSuperGrid[grid];
      RealArray & etay = etaySuperGrid[grid];
      RealArray & etaz = etazSuperGrid[grid];
      
      Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
      getIndex( mg.dimension(),I1,I2,I3 );          // all points including ghost points.
      OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
      bool ok = ParallelUtility::getLocalArrayBounds(mg.mask(),maskLocal,I1,I2,I3,1);   
      if( !ok ) return 0;  // no points on this processor (NOTE: no communication should be done after this point)

      if( useAbsorbingLayer(0,grid) )
      {
	Range D1=maskLocal.dimension(0);
        etax.redim(D1); etax=1.;         // local dimensions of etax should match dimensions of uLocal (for call to advBA)
      }
      if( useAbsorbingLayer(1,grid) )
      {
	Range D2=maskLocal.dimension(1);
        etay.redim(D2);
        etay=1.;
      }
      if( useAbsorbingLayer(2,grid) )
      {
	Range D3=maskLocal.dimension(2);
        etaz.redim(D3);
        etaz=1.;
      }
      

      // -- we optimize for Cartesian grids (we can avoid creating the vertex array)
      if( !isRectangular )
      {
        mg.update(MappedGrid::THEvertex | MappedGrid::THEcenter);
        OV_GET_SERIAL_ARRAY(real,mg.center(),xLocal);

        OV_ABORT("buildSuperGrid: FINISH ME FOR CURVILINEAR GRIDS");
      }
      
      real dvx[3]={1.,1.,1.}, xab[2][3]={{0.,0.,0.},{0.,0.,0.}};
      int iv0[3]={0,0,0}; //
      int iv[3], &i1=iv[0], &i2=iv[1], &i3=iv[2];  // NOTE: iv[0]==i1, iv[1]==i2, iv[2]==i3
      real xv[3]={0.,0.,0.};
      if( isRectangular )
      {
        mg.getRectangularGridParameters( dvx, xab );
        for( int dir=0; dir<mg.numberOfDimensions(); dir++ )
        {
          iv0[dir]=mg.gridIndexRange(0,dir);
          if( mg.isAllCellCentered() )
            xab[0][dir]+=.5*dvx[dir];  // offset for cell centered
        }
      }
      // This macro defines the grid points for rectangular grids:
#undef XC
#define XC(iv,axis) (xab[0][axis]+dvx[axis]*(iv[axis]-iv0[axis]))
   

      for( int axis=0; axis<numberOfDimensions; axis++)
      {
        if( useAbsorbingLayer(axis,grid) )
        {
          // ---- Build layer function for this axis ----
          RealArray & eta = axis==0 ? etax : axis==1 ? etay : etaz;
          for( int dir=0; dir<3; dir++ ){ iv[dir]=gid(0,dir); } // maybe not needed 

          for( int i=Iv[axis].getBase(); i<=Iv[axis].getBound(); i++ )
          {
            assert( i>=eta.getBase(0) && i<=eta.getBound(0) );

            iv[axis]=i;
            real x = XC(iv,axis);
            if( bc(0,axis) == absorbing ) 
            {
              real z = xab[0][axis] +superGridWidth - x;
              if( z>superGridWidth )
                eta(i)=epsL;
              else if( z>0 )
                eta(i)=etaSG(z);

              // printF(" i=%d, x=%e, z=%e, eta=%e\n",i,x,z,eta(i));
            }
            if( bc(1,axis) ==absorbing )
            {

              real z = x - ( xab[1][axis] -superGridWidth); 
              if( z>superGridWidth )
                eta(i)=epsL;
              else if( z>0 )
                eta(i)=etaSG(z);
            }
          }
        } // end for buildLayers
      } // end for axis
    

      bool plotLayers=false;
      if( plotLayers )
      {
        assert( useAbsorbingLayer(0,grid) );
        
        assert( gip !=NULL );
        GenericGraphicsInterface & gi = *gip;

        PlotStuffParameters psp;
        aString title=sPrintF("Super-grid functions: width=%g p=%g q=%g",superGridWidth,pSG,qSG);
        // int numFields = 1;
        // int numPoints = dimension(1,0)-dimension(0,0)+1;
        // RealArray fields(numPoints,numFields);

        const aString names[]={"etax"};

        RealArray xv(I1);
        int axis=0;
        for( int i=dimension(0,axis); i<=dimension(1,axis); i++ )
        {
          iv[axis]=i;
          xv(i)=XC(iv,axis);
        }
         
        ::display(xv,"xv","%8.2e ");
        ::display(etax,"etax","%8.2e ");

        psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);

        #ifndef USE_PPP
           PlotIt::plot(gi, xv, etax, title, "super-grid", names,psp );
        // PlotIt::plot(gi, xv, etay, title, "super-grid", names,psp );
        #else
         printF("DMP: FINISH PLOTTING IN PARALLEL\n");
        #endif

      }
      


    }  // end for grid

  }
  
  return 0;
}
