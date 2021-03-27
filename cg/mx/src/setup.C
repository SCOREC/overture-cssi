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



//! Setup and initialization. Build the grid.
int Maxwell::
setupGrids()
// ===================================================================================
// Build a grid
// ===================================================================================
{
  real time0=getCPU();

  // real xa=0, xb=1., ya=0., yb=1.;
  real xa=xab[0][0], xb=xab[1][0], ya=xab[0][1], yb=xab[1][1], za=xab[0][2], zb=xab[1][2];
  
  real dx[3];
  dx[0]=(xb-xa)/(nx[0]-1);
  dx[1]=(yb-ya)/(nx[1]-1);
  dx[2]=(zb-za)/(nx[2]-1);

//   xa=0, xb=1., ya=0., yb=1.25;
//   SquareMapping rectangle(xa,xb,ya,yb);
//square.setMappingCoordinateSystem( Mapping::general ); 

  Mapping *mapPointer;

  if( gridType==square && elementType==structuredElements )
  {
    mapPointer= new SquareMapping(xa,xb,ya,yb);  mapPointer->incrementReferenceCount();
  }
  else if ( gridType==box )
  {
    mapPointer = new BoxMapping(xa,xb,ya,yb,za,zb); mapPointer->incrementReferenceCount();

    for( int axis=0; axis<mapPointer->getDomainDimension(); axis++ )
    {
      mapPointer->setGridDimensions(axis,nx[axis]);
      if( bcOption==useAllPeriodicBoundaryConditions )
	mapPointer->setIsPeriodic(axis,Mapping::derivativePeriodic);
    }

    if ( elementType!=structuredElements )
    {
      UnstructuredMapping & uns = * new UnstructuredMapping;
      BoxMapping &bmap = (BoxMapping &)*mapPointer;
      mapPointer->decrementReferenceCount();
      mapPointer= &uns;  mapPointer->incrementReferenceCount();
      if ( bcOption==useAllPeriodicBoundaryConditions )
	uns.addGhostElements(true);
      else
	uns.addGhostElements(false);

      printF(" **** elementType=%i \n",elementType);
	  
      //	  UnstructuredMapping etype;
      uns.addGhostElements(true);
      uns.buildFromAMapping(bmap);
      verifyUnstructuredConnectivity(uns,true);
      for( int axis=0; axis<bmap.getDomainDimension() && bcOption==useAllPeriodicBoundaryConditions; axis++ )
      {
	uns.setIsPeriodic(axis,Mapping::derivativePeriodic);
	uns.setBoundaryCondition(Start,axis,-1);
	uns.setBoundaryCondition(End  ,axis,-1);
      }
      //	  if ( !bcOption==useAllPeriodicBoundaryConditions )
      //	    uns.expandGhostBoundary();

    }
  }
  else if( gridType==annulus )
  {
    mapPointer= new AnnulusMapping;  mapPointer->incrementReferenceCount();
  }
  else if( gridType==rotatedSquare )
  {
    Mapping *sq=new SquareMapping(xa,xb,ya,yb); 
    sq->incrementReferenceCount();
    MatrixTransform & mat = * new MatrixTransform(*sq);
    mapPointer= &mat;   mapPointer->incrementReferenceCount();
    if( sq->decrementReferenceCount()==0 ) delete sq;
    
    // mat.rotate(axis3,90.*twoPi/360.);
    // mat.rotate(axis3,45.*twoPi/360.);
    mat.rotate(axis3,30.*twoPi/360.);
    for( int axis=0; axis<mapPointer->getDomainDimension(); axis++ )
    {
      mapPointer->setGridDimensions(axis,nx[axis]);
    }

    UnstructuredMapping & uns = * new UnstructuredMapping;
    mapPointer->decrementReferenceCount();
    mapPointer= &uns;  mapPointer->incrementReferenceCount();
    uns.addGhostElements(true);

    printF(" **** elementType=%i \n",elementType);
      
    int domainDimension = 2;
    if ( elementType==defaultUnstructured && domainDimension==2)
      uns.buildFromARegularMapping(mat);
    else if ( domainDimension==2 )
    {
      uns.buildFromARegularMapping(mat,elementType==triangles ? UnstructuredMapping::triangle : elementType==quadrilaterals ? UnstructuredMapping::quadrilateral : UnstructuredMapping::hexahedron);
    }
    else
    {
      uns.buildFromAMapping(mat);
      for( int axis=0; axis<uns.getDomainDimension() && 
	     bcOption==useAllPeriodicBoundaryConditions; axis++ )
      {
	uns.setIsPeriodic(axis,Mapping::derivativePeriodic);
	uns.setBoundaryCondition(Start,axis,-1);
	uns.setBoundaryCondition(End  ,axis,-1);
      }
    }
    

    //            verifyUnstructuredConnectivity(uns,true);
    
    if( mat.decrementReferenceCount() == 0 )
      delete &mat;

  }
  else if( gridType==skewedSquare || gridType==sineSquare ||
           gridType==chevron || gridType==chevbox || gridType==sineByTriangles )
  {
    // build a DataPointMapping
    DataPointMapping & dpm = *new DataPointMapping;
    mapPointer= &dpm;   mapPointer->incrementReferenceCount();
    Index I1,I2,I3;

    // *wdh* 050818 -- use 2nd-order interpolation so the metrics are periodic when the grid is periodic
//    dpm.setOrderOfInterpolation(4);
    dpm.setOrderOfInterpolation(2);

//    int numberOfGhostLines=0; // add more for 4th-order and periodic derivatives of mapping at ghost pts (?)
    // int numberOfGhostLines=4; // add more for 4th-order and periodic derivatives of mapping at ghost pts (?)
    int numberOfGhostLines=max(4,(orderOfAccuracyInSpace+2)/2);
//   int numberOfGhostLines=8; // add more for 4th-order and periodic derivatives of mapping at ghost pts (?)

    I1=Range(-numberOfGhostLines,nx[0]+numberOfGhostLines-1);
    I2=Range(-numberOfGhostLines,nx[1]+numberOfGhostLines-1);
    I3= gridType!=chevbox ? Range(0,0) : Range(-numberOfGhostLines,nx[2]+numberOfGhostLines-1);
    
    int rDim = gridType==chevbox ? 3 : 2;
    realArray x(I1,I2,I3,rDim);
    realArray r1(I1),r2(I2),r3(I3);

    r1.seqAdd(-numberOfGhostLines*dx[0],dx[0]);
    r2.seqAdd(-numberOfGhostLines*dx[1],dx[1]);
    if ( gridType==chevbox )
      r3.seqAdd(-numberOfGhostLines*dx[2],dx[2]);

    if( gridType==sineSquare || gridType==sineByTriangles )
    {
      realArray bottom(I1);
      real amplitude=.1;    // amplitude of the sine wave
      bottom=(cos(twoPi*r1)-1.)*amplitude;
    
      int i2=0, i3=0;
      for( i2=I2.getBase(); i2<=I2.getBound(); i2++ )
      {
	x(I1,i2,i3,0)=r1(I1);
	x(I1,i2,i3,1)=bottom(I1)+r2(i2);
      }
    }
    else if( gridType==chevron || gridType==chevbox)
    {
      real amplitudeY=dx[1]*chevronAmplitude;    // amplitude of the chevron oscillation
      realArray yShift(I1);
      // .5 seemed stable, .75 was not stable
      real freqX=(Pi/dx[0])*chevronFrequency; 
      // real freq=(Pi/dx[0])*.75;  // this is a lower frequency perturbation

      bool useMoreFrequencies=true;
      if( !useMoreFrequencies )
      {
        yShift=cos(freqX*r1)*amplitudeY;
      }
      else
      {
        yShift=cos(freqX*r1)*amplitudeY + cos(.5*freqX*r1)*amplitudeY*1.5;
      }
      
      // display(yShift,"yShift for the chevron grid");
      

      // remove perturbation from the ends
      bool removePerturbationAtEdges=false; // false; // true;
      int extra=numberOfGhostLines+5;
      if( removePerturbationAtEdges )
      {
	yShift(Range(I1.getBase(),I1.getBase()+extra))=0.;
	yShift(Range(I1.getBound()-extra,I1.getBound()))=0.;
      }
      
      // here is a perturbation in the x-direction
      real amplitudeX=dx[0]*chevronAmplitude;  
      realArray xShift(I2);
      real freqY=(Pi/dx[1])*chevronFrequency; 
      xShift=cos(freqY*r2)*amplitudeX; 
      if( useMoreFrequencies )
      {
        xShift=cos(freqY*r2)*amplitudeX+ cos(.5*freqY*r2)*amplitudeX*1.5;
      }
      
      if( removePerturbationAtEdges )
      {
	xShift(Range(I2.getBase(),I2.getBase()+extra))=0.;
	xShift(Range(I2.getBound()-extra,I2.getBound()))=0.;
      }

      int i2=0, i3=0;
      if ( gridType==chevron )
      {
        if( true )
	{
	  // box with random perturbations
          int seed=184273654;
          srand(seed);
	  int i1;
	  for( i2=I2.getBase(); i2<=I2.getBound(); i2++ )
	    for( i1=I1.getBase(); i1<=I1.getBound(); i1++ )
	    {
	      real d1=(-1.+2.*rand()/RAND_MAX)*.25*dx[0];  // random number between [-.25,.25]*dx[0]
	      real d2=(-1.+2.*rand()/RAND_MAX)*.25*dx[1];
		
	      x(i1,i2,i3,0)= i1*dx[0]+d1;
	      x(i1,i2,i3,1)= i2*dx[1]+d2;
	    }

          // make periodic
          real xba=1., yba=1., zba=1.;
          Range Rx=2;
          for( i1=I1.getBase(); i1<0; i1++ )
	  {
            x(i1,I2,I3,Rx)=x(i1+nx[0]-1,I2,I3,Rx);
            x(i1,I2,I3,0)-=xba; 
	  }
          for( i1=nx[0]-1; i1<=I1.getBound(); i1++ )
	  {
            x(i1,I2,I3,Rx)=x(i1-nx[0]+1,I2,I3,Rx); 
            x(i1,I2,I3,0)+=xba; 
	  }
	  
          for( i2=I2.getBase(); i2<0; i2++ )
	  {
            x(I1,i2,I3,Rx)=x(I1,i2+nx[1]-1,I3,Rx);
            x(I1,i2,I3,1)-=yba; 
	  }
          for( i2=nx[1]-1; i2<=I2.getBound(); i2++ )
	  {
            x(I1,i2,I3,Rx)=x(I1,i2-nx[1]+1,I3,Rx); 
            x(I1,i2,I3,1)+=yba; 
	  }
	  
	}
	else
	{
	  for( i2=I2.getBase(); i2<=I2.getBound(); i2++ )
	  {
	    x(I1,i2,i3,0)=r1(I1);
	    if( !removePerturbationAtEdges ||
		(i2>=I2.getBase() + extra &&
		 i2<=I2.getBound()- extra) )
	    {
	      x(I1,i2,i3,1)=yShift(I1)+r2(i2);
	    }
	    else
	    {
	      x(I1,i2,i3,1)=r2(i2);
	    }
	  
	  }
	  if( amplitudeX != 0. )
	  {
	    for( i2=I2.getBase(); i2<=I2.getBound(); i2++ )
	      x(I1,i2,i3,0)+=xShift(i2);
	  }
	}
	
      }
      else 
      {
        // chevBox
        if( true )
	{
	  // box with random perturbations
          int seed=184273654;
          srand(seed);
	  int i1;
	  for( i3=I3.getBase(); i3<=I3.getBound(); i3++ )
	    for( i2=I2.getBase(); i2<=I2.getBound(); i2++ )
	      for( i1=I1.getBase(); i1<=I1.getBound(); i1++ )
	      {
                real d1=(-1.+2.*rand()/RAND_MAX)*.25*dx[0];  // random number between [-.25,.25]*dx[0]
                real d2=(-1.+2.*rand()/RAND_MAX)*.25*dx[1];
                real d3=(-1.+2.*rand()/RAND_MAX)*.25*dx[2];
		
                // if( i1<=0 || i1>=nx[0]-1 ) d1=0.;
                // if( i2<=0 || i2>=nx[1]-1 ) d2=0.;
                // if( i3<=0 || i3>=nx[2]-1 ) d3=0.;

		x(i1,i2,i3,0)= i1*dx[0]+d1;
		x(i1,i2,i3,1)= i2*dx[1]+d2;
		x(i1,i2,i3,2)= i3*dx[2]+d3;
	      }

          // make periodic
          real xba=1., yba=1., zba=1.;
          Range Rx=3;
          for( int mm=0; mm<=1; mm++ )
	  {
	    for( i1=I1.getBase(); i1<0; i1++ )
	    {
	      x(i1,I2,I3,Rx)=x(i1+nx[0]-1,I2,I3,Rx);
	      x(i1,I2,I3,0)-=xba; 
	    }
	    for( i1=nx[0]-1; i1<=I1.getBound(); i1++ )
	    {
	      x(i1,I2,I3,Rx)=x(i1-nx[0]+1,I2,I3,Rx); 
	      x(i1,I2,I3,0)+=xba; 
	    }
	  
	    for( i2=I2.getBase(); i2<0; i2++ )
	    {
	      x(I1,i2,I3,Rx)=x(I1,i2+nx[1]-1,I3,Rx);
	      x(I1,i2,I3,1)-=yba; 
	    }
	    for( i2=nx[1]-1; i2<=I2.getBound(); i2++ )
	    {
	      x(I1,i2,I3,Rx)=x(I1,i2-nx[1]+1,I3,Rx); 
	      x(I1,i2,I3,1)+=yba; 
	    }
	  
	    for( i3=I3.getBase(); i3<0; i3++ )
	    {
	      x(I1,I2,i3,Rx)=x(I1,I2,i3+nx[2]-1,Rx);
	      x(I1,I2,i3,2)-=zba; 
	    }
	    for( i3=nx[2]-3; i3<=I3.getBound(); i3++ )
	    {
	      x(I1,I2,i3,Rx)=x(I1,I2,i3-nx[2]+1,Rx);
	      x(I1,I2,i3,2)+=zba; 
	    }
	  }
	  

	}
	else
	{
	  for( i3=I3.getBase(); i3<=I3.getBound(); i3++ )
	    for( i2=I2.getBase(); i2<=I2.getBound(); i2++ )
	    {
	      x(I1,i2,i3,0)=r1(I1);
	      x(I1,i2,i3,1)=yShift(I1)+r2(i2);
	    }

	  for( i3=I3.getBase(); i3<=I3.getBound(); i3++ )
	    x(I1,I2,i3,2)=r3(i3);

	  if( amplitudeX != 0. )
	  {
	    for( i2=I2.getBase(); i2<=I2.getBound(); i2++ )
	      x(I1,i2,I3,0)+=xShift(i2);
	  }
	  real amplitudeZ=dx[2]*chevronAmplitude;  
	  if( amplitudeZ != 0. )
	  {
	
	    realArray zShift(I3);
	    real freqY=(Pi/dx[1])*chevronFrequency; 
	    zShift=cos(freqY*r3)*amplitudeZ; 
	    if( useMoreFrequencies )
	    {
	      zShift=cos(freqY*r3)*amplitudeZ+ cos(.5*freqY*r3)*amplitudeZ*1.5;
	    }
	    for( i3=I3.getBase(); i3<=I3.getBound(); i3++ )
	    {
	      x(I1,I2,i3,2)+=zShift(i3);
	    }
	  }
	}
	
      }

    }
    else
    {
      OV_ABORT("error");
    }
    // display(x,"Data points for DPM");

    const int domainDimension= gridType==chevbox ? 3 : 2;
    int ng=max(2,numberOfGhostLines);
    // int numberOfGhostLinesNew[2][3]={ng,ng,ng,ng,ng,ng}; //
    Mapping::IndexRangeType numberOfGhostLinesNew;
    for( int axis=0; axis<3; axis++ )for( int side=0; side<=1; side++ ) numberOfGhostLinesNew(side,axis)=ng;
    
    // tell the dpm to use the extra ghost lines (*wdh* 050818)
    dpm.setDomainDimension(domainDimension); 
    dpm.setRangeDimension(domainDimension); 
    dpm.setNumberOfGhostLines(numberOfGhostLinesNew);
    
    // Set periodic BC's before supplying data points so the ghost points are computed properly
    // *************** wdh* 050818 **
    if( bcOption==useAllPeriodicBoundaryConditions )
    {
      for( int axis=0; axis<domainDimension; axis++ )
      {
 	dpm.setIsPeriodic(axis,Mapping::derivativePeriodic);
 	dpm.setBoundaryCondition(Start,axis,-1);
 	dpm.setBoundaryCondition(End  ,axis,-1);
      }
    }

    dpm.setDataPoints(x,3,domainDimension,numberOfGhostLines);

    if( gridType==sineByTriangles || elementType!=structuredElements )
    {
      int axis;
      for( axis=0; axis<domainDimension; axis++ )
      {
	dpm.setIsPeriodic(axis,Mapping::derivativePeriodic);
	dpm.setBoundaryCondition(Start,axis,-1);
	dpm.setBoundaryCondition(End  ,axis,-1);
      }

      UnstructuredMapping & uns = * new UnstructuredMapping;
      mapPointer->decrementReferenceCount();
      mapPointer= &uns;  mapPointer->incrementReferenceCount();
      if ( bcOption==useAllPeriodicBoundaryConditions )
	uns.addGhostElements(true);
      else
	uns.addGhostElements(false);

      printF(" **** elementType=%i \n",elementType);
      

      if ( elementType==defaultUnstructured && domainDimension==2)
	uns.buildFromARegularMapping(dpm);
      else if ( domainDimension==2 )
      {
	uns.buildFromARegularMapping(dpm,elementType==triangles ? UnstructuredMapping::triangle : elementType==quadrilaterals ? UnstructuredMapping::quadrilateral : UnstructuredMapping::hexahedron);
      }
      else
      {
// 	  uns.buildFromAMapping(dpm);
// 	  verifyUnstructuredConnectivity(uns,true);

	uns.addGhostElements(true);
	uns.buildFromAMapping(dpm);
	for( int axis=0; axis<uns.getDomainDimension() && 
	       bcOption==useAllPeriodicBoundaryConditions; axis++ )
	{
	  uns.setIsPeriodic(axis,Mapping::derivativePeriodic);
	  uns.setBoundaryCondition(Start,axis,-1);
	  uns.setBoundaryCondition(End  ,axis,-1);
	}
      }

      //            verifyUnstructuredConnectivity(uns,true);

      if( dpm.decrementReferenceCount() == 0 )
        delete &dpm;
      
    }
  }
  else if( gridType==perturbedSquare || gridType==perturbedBox )
  {
    // ************ perturbed square or box **********

    const int rangeDimension = gridType==perturbedSquare ? 2 : 3; 
    const int domainDimension= rangeDimension; 

    // build a DataPointMapping
    DataPointMapping & dpm = *new DataPointMapping;
    mapPointer= &dpm;   mapPointer->incrementReferenceCount();
    Index I1,I2,I3;

    dpm.setOrderOfInterpolation(2);

    int numberOfGhostLines=max(4,(orderOfAccuracyInSpace+2)/2);

    I1=Range(-numberOfGhostLines,nx[0]+numberOfGhostLines-1);
    I2=Range(-numberOfGhostLines,nx[1]+numberOfGhostLines-1);
    I3= domainDimension==2 ? Range(0,0) : Range(-numberOfGhostLines,nx[2]+numberOfGhostLines-1);
    
    realArray x(I1,I2,I3,rangeDimension);

    int seed=184273654;
    srand(seed);          // supply seed to the random number generator

    int i1=0, i2=0, i3=0;
    if ( gridType==perturbedSquare )
    {
      //  === square with random perturbations ===

      for( i2=I2.getBase(); i2<=I2.getBound(); i2++ )
	for( i1=I1.getBase(); i1<=I1.getBound(); i1++ )
	{
	  real d1=(-1.+2.*rand()/RAND_MAX)*.25*dx[0];  // random number between [-.25,.25]*dx[0]
	  real d2=(-1.+2.*rand()/RAND_MAX)*.25*dx[1];
		
	  x(i1,i2,i3,0)= i1*dx[0]+d1;
	  x(i1,i2,i3,1)= i2*dx[1]+d2;
	}

      // make periodic
      real xba=1., yba=1., zba=1.;
      Range Rx=2;
      for( i1=I1.getBase(); i1<0; i1++ )
      {
	x(i1,I2,I3,Rx)=x(i1+nx[0]-1,I2,I3,Rx);
	x(i1,I2,I3,0)-=xba; 
      }
      for( i1=nx[0]-1; i1<=I1.getBound(); i1++ )
      {
	x(i1,I2,I3,Rx)=x(i1-nx[0]+1,I2,I3,Rx); 
	x(i1,I2,I3,0)+=xba; 
      }
	  
      for( i2=I2.getBase(); i2<0; i2++ )
      {
	x(I1,i2,I3,Rx)=x(I1,i2+nx[1]-1,I3,Rx);
	x(I1,i2,I3,1)-=yba; 
      }
      for( i2=nx[1]-1; i2<=I2.getBound(); i2++ )
      {
	x(I1,i2,I3,Rx)=x(I1,i2-nx[1]+1,I3,Rx); 
	x(I1,i2,I3,1)+=yba; 
      }
	
    }
    else 
    {
      //  === Box with random perturbations ===
      // box with random perturbations
      for( i3=I3.getBase(); i3<=I3.getBound(); i3++ )
	for( i2=I2.getBase(); i2<=I2.getBound(); i2++ )
	  for( i1=I1.getBase(); i1<=I1.getBound(); i1++ )
	  {
	    real d1=(-1.+2.*rand()/RAND_MAX)*.25*dx[0];  // random number between [-.25,.25]*dx[0]
	    real d2=(-1.+2.*rand()/RAND_MAX)*.25*dx[1];
	    real d3=(-1.+2.*rand()/RAND_MAX)*.25*dx[2];
		
	    // if( i1<=0 || i1>=nx[0]-1 ) d1=0.;
	    // if( i2<=0 || i2>=nx[1]-1 ) d2=0.;
	    // if( i3<=0 || i3>=nx[2]-1 ) d3=0.;

	    x(i1,i2,i3,0)= i1*dx[0]+d1;
	    x(i1,i2,i3,1)= i2*dx[1]+d2;
	    x(i1,i2,i3,2)= i3*dx[2]+d3;
	  }

      // make periodic
      real xba=1., yba=1., zba=1.;
      Range Rx=3;
      for( int mm=0; mm<=1; mm++ )  // loop twice to make sure all corners points are correct
      {
	for( i1=I1.getBase(); i1<0; i1++ )
	{
	  x(i1,I2,I3,Rx)=x(i1+nx[0]-1,I2,I3,Rx);
	  x(i1,I2,I3,0)-=xba; 
	}
	for( i1=nx[0]-1; i1<=I1.getBound(); i1++ )
	{
	  x(i1,I2,I3,Rx)=x(i1-nx[0]+1,I2,I3,Rx); 
	  x(i1,I2,I3,0)+=xba; 
	}
	  
	for( i2=I2.getBase(); i2<0; i2++ )
	{
	  x(I1,i2,I3,Rx)=x(I1,i2+nx[1]-1,I3,Rx);
	  x(I1,i2,I3,1)-=yba; 
	}
	for( i2=nx[1]-1; i2<=I2.getBound(); i2++ )
	{
	  x(I1,i2,I3,Rx)=x(I1,i2-nx[1]+1,I3,Rx); 
	  x(I1,i2,I3,1)+=yba; 
	}
	  
	for( i3=I3.getBase(); i3<0; i3++ )
	{
	  x(I1,I2,i3,Rx)=x(I1,I2,i3+nx[2]-1,Rx);
	  x(I1,I2,i3,2)-=zba; 
	}
	for( i3=nx[2]-3; i3<=I3.getBound(); i3++ )
	{
	  x(I1,I2,i3,Rx)=x(I1,I2,i3-nx[2]+1,Rx);
	  x(I1,I2,i3,2)+=zba; 
	}
      }
      
    }
    
    int ng=max(2,numberOfGhostLines);
    // int numberOfGhostLinesNew[2][3]={ng,ng,ng,ng,ng,ng}; //
    Mapping::IndexRangeType numberOfGhostLinesNew;
    for( int axis=0; axis<3; axis++ )for( int side=0; side<=1; side++ ) numberOfGhostLinesNew(side,axis)=ng;
    
    dpm.setDomainDimension(domainDimension); 
    dpm.setRangeDimension(rangeDimension); 
    dpm.setNumberOfGhostLines(numberOfGhostLinesNew);
    
    // Set periodic BC's before supplying data points so the ghost points are computed properly
    if( bcOption==useAllPeriodicBoundaryConditions )
    {
      for( int axis=0; axis<domainDimension; axis++ )
      {
 	dpm.setIsPeriodic(axis,Mapping::derivativePeriodic);
 	dpm.setBoundaryCondition(Start,axis,-1);
 	dpm.setBoundaryCondition(End  ,axis,-1);
      }
    }

    dpm.setDataPoints(x,3,domainDimension,numberOfGhostLines);

    
  }
  else if( gridType==squareByTriangles || gridType==squareByQuads ||
	   (gridType==square && elementType==Maxwell::defaultUnstructured) )
  {
    Mapping & sq= *new SquareMapping(xa,xb,ya,yb); 
    sq.incrementReferenceCount();

    int side,axis;
    for( axis=0; axis<sq.getDomainDimension(); axis++ )
    {
      sq.setGridDimensions(axis,nx[axis]);
      sq.setIsPeriodic(axis,Mapping::derivativePeriodic);
      sq.setBoundaryCondition(Start,axis,-1);
      sq.setBoundaryCondition(End  ,axis,-1);
    }

    UnstructuredMapping & uns = * new UnstructuredMapping;
    mapPointer= &uns;   mapPointer->incrementReferenceCount();
    uns.addGhostElements(true);
    if( gridType==squareByTriangles )
    {
      uns.buildFromARegularMapping(sq,UnstructuredMapping::triangle);
    }
    else
    {

      uns.buildFromARegularMapping(sq,UnstructuredMapping::quadrilateral);

//        // For debugging:
//        UnstructuredMappingIterator iter;
//        for( iter=uns.begin(UnstructuredMapping::Region); iter!=uns.end(UnstructuredMapping::Region); iter++ )
//        {
//  	printf(" Element %i is valid \n",*iter);
//        }

    }

    //    uns.getEntities(UnstructuredMapping::Face).display("Faces");
    //    uns.getEntities(UnstructuredMapping::Edge).display("Edges");
    //        verifyUnstructuredConnectivity(uns,true);

//      GenericGraphicsInterface & ps = *gip;
//      GraphicsParameters params;
//      params.set(GI_PLOT_UNS_EDGES,true);
//      params.set(GI_PLOT_UNS_FACES,true);

//      PlotIt::plot(ps,uns,params);
    
    if( sq.decrementReferenceCount()==0 ) delete &sq;
  }
  else if( gridType==compositeGrid )
  {
    // *wdh* 090427 The CompositeGrid is now created in the main program

    //  // In this case we grab the Mapping from the first component grid.
    //  // create and read in a CompositeGrid
    //  delete cgp;
    //  cgp=new CompositeGrid;
    //  CompositeGrid & cg = *cgp;
    //  getFromADataBase(cg,nameOfGridFile);

    assert( cgp!=NULL);
    CompositeGrid & cg = *cgp;
    mapPointer=&(cg[0].mapping().getMapping());
  }
  else
  {
    printF("Cgmx:ERROR: unknown gridType=",(int)gridType);
    OV_ABORT("ERROR");
  }
  
  if( cgp==NULL )
  {
    // ****************************************
    // *********** single grid  ***************
    // ****************************************

    Mapping & map = *mapPointer;
  
    bool unstructured=!(elementType==structuredElements);//map.getClassName()=="UnstructuredMapping";

    int side,axis;
    if( !unstructured )
    {
      for( axis=0; axis<map.getDomainDimension(); axis++ )
      {
	map.setGridDimensions(axis,nx[axis]);
	if( map.getIsPeriodic(axis)==Mapping::notPeriodic )
	{
	  if( bcOption==useAllPeriodicBoundaryConditions )
	  {
	    map.setIsPeriodic(axis,Mapping::derivativePeriodic);
	    map.setBoundaryCondition(Start,axis,-1);
	    map.setBoundaryCondition(End  ,axis,-1);
	  }
	  else if( bcOption==useAllDirichletBoundaryConditions )
	  {
	    map.setBoundaryCondition(Start,axis,dirichlet);
	    map.setBoundaryCondition(End  ,axis,dirichlet);
	  }
	  else if( bcOption==useAllPerfectElectricalConductorBoundaryConditions )
	  {
	    map.setBoundaryCondition(Start,axis,perfectElectricalConductor);
	    map.setBoundaryCondition(End  ,axis,perfectElectricalConductor);
	  }
	}
      
      }
    }
  
    // Build a MappedGrid
    mgp = new MappedGrid(map);
    MappedGrid & mg = *mgp;
  
    if( !unstructured )
    {
      const int numberOfGhostPoints=max(2,orderOfAccuracyInSpace/2); 
      printF(">>>>Create the MappedGrid with %i ghost points, orderOfAccuracyInSpace=%i\n",numberOfGhostPoints,
	     orderOfAccuracyInSpace  );
      
      for( axis=0; axis<map.getDomainDimension(); axis++ )
	for( side=Start; side<=End; side++ )
	  mg.setNumberOfGhostPoints(side,axis,numberOfGhostPoints);
  
    }

    cgp= new CompositeGrid(mg.numberOfDimensions(),1);
    CompositeGrid & cg = *cgp;
    cg[0].reference(mg);
    cg.updateReferences();
    
//    mg.update(MappedGrid::THEcenter | MappedGrid::THEvertex );  

//    cg.update(MappedGrid::THEmask | MappedGrid::THEcenter | MappedGrid::THEvertex );  
    cg.update( MappedGrid::THEcenter | MappedGrid::THEvertex | MappedGrid::THEmask );  

    if ( unstructured )
    {
      cg.update(MappedGrid::THEcorner| MappedGrid::THEfaceArea | MappedGrid::THEfaceNormal | MappedGrid::THEcellVolume | MappedGrid::THEcenterNormal | MappedGrid::THEcenterArea );

      for ( int g=0; false && bcOption==useAllPeriodicBoundaryConditions && g<cg.numberOfGrids(); g++ )
      {
	// we do this trick for convenience when working with periodic boundaries and dsi schemes
	bool vCent = mg.isAllVertexCentered();
	realArray &cFArea = vCent ? mg.centerArea() : mg.faceArea();
	realArray &cFNorm = vCent ? mg.centerNormal() : mg.faceNormal();
	realArray &cEArea = vCent ? mg.faceArea() : mg.centerArea();
	realArray &cENorm = vCent ? mg.faceNormal() : mg.centerNormal();

	MappedGrid &mg = cg[g];
	const IntegerArray &perEImages = *mg.getUnstructuredPeriodicBC(UnstructuredMapping::Edge);
	    
	for ( int e=perEImages.getBase(0); e<=perEImages.getBound(0); e++ )
	{
	  cEArea(perEImages(e,0),0,0) = cEArea(perEImages(e,1),0,0);
	  for ( int r=0; r<mg.numberOfDimensions(); r++ )
	    cENorm(perEImages(e,0),0,0,r) = cENorm(perEImages(e,1),0,0,r);
	}

	const IntegerArray &perHImages = *mg.getUnstructuredPeriodicBC(UnstructuredMapping::Face);
	for ( int h=perHImages.getBase(0); mg.numberOfDimensions()>2 &&h<=perHImages.getBound(0); h++ )
	{
	  cFArea(perHImages(h,0),0,0) = cFArea(perHImages(h,1),0,0);
	  for ( int r=0; r<mg.numberOfDimensions(); r++ )
	    cFNorm(perHImages(h,0),0,0,r) = cFNorm(perHImages(h,1),0,0,r);
	}

      }
    }

//  display(mg.vertex(),"mg.vertex()");
    // display(mg.center(),"mg.center()");

//   display(mg.dimension(),"mg.dimension()");
//   display(mg.gridIndexRange(),"mg.gridIndexRange()");
  
//    GenericGraphicsInterface & ps = *gip;
//    PlotIt::plot(ps,mg);
  
//  mg.updateReferences();

//   mg.update(MappedGrid::THEcenterJacobian | MappedGrid::THEinverseCenterDerivative );
//   display(mg.center(),"mg.center()");
//   display(mg.centerJacobian(),"mg.centerJacobian()");
//   display(mg.inverseCenterDerivative(),"inverseCenterDerivative");
  

    if( method==defaultMethod )
    {
      if( mg.isRectangular() )
	method=yee;
      else
	method=dsi;
    }
  
    if( mapPointer->decrementReferenceCount()==0 ) delete mapPointer;

  }
  else
  {
    // *********************************************
    // *********** CompositeGrid *******************
    // *********************************************

    CompositeGrid & cg = *cgp;

    //    cg.update( MappedGrid::THEcenter | MappedGrid::THEvertex | MappedGrid::THEmask );  

    if ( cg[0].getGridType()==MappedGrid::unstructuredGrid )
    {
      UnstructuredMapping &umap = (UnstructuredMapping &) cg[0].mapping().getMapping();
      umap.expandGhostBoundary();
      verifyUnstructuredConnectivity(umap,true);
      //	umap.expandGhostBoundary();
      //	verifyUnstructuredConnectivity(umap,true);

      cg.destroy( MappedGrid::THEcenter | MappedGrid::THEvertex | MappedGrid::THEmask |
		  MappedGrid::THEcorner | MappedGrid::THEcellVolume | MappedGrid::THEcenterNormal |
		  MappedGrid::THEfaceArea | MappedGrid::THEfaceNormal | 
		  MappedGrid::THEcellVolume  | MappedGrid::THEcenterArea );	

      cg.update( MappedGrid::THEcenter | MappedGrid::THEvertex | MappedGrid::THEmask |
		 MappedGrid::THEcorner | MappedGrid::THEcellVolume | MappedGrid::THEcenterNormal |
		 MappedGrid::THEfaceArea | MappedGrid::THEfaceNormal | 
		 MappedGrid::THEcellVolume  | MappedGrid::THEcenterArea );	

    }
    else
    {
      // cg.update(MappedGrid::THEmask );
      
      // *wdh* 031202 cg.update(MappedGrid::THEcenter | MappedGrid::THEvertex );  
    }

    pinterpolant = new Interpolant(cg);
    pinterpolant->incrementReferenceCount();

    // we allow implicit interpolation in parallel now *wdh* March 20, 2021
// #ifdef USE_PPP
//     if( !pinterpolant->interpolationIsExplicit() )
//     {
//       printF("*** ERROR: The parallel composite grid interpolator needs explicit interpolation ****\n");
//       Overture::abort();
//     }
// #endif

    //kkc can now use dsi, this would override the command line spec    method=nfdtd;
    
    // Set the default order of accuracy from the grid parameters
    int minDiscretizationWidth=INT_MAX;
    int minInterpolationWidth=INT_MAX;
    Range R=cg.numberOfDimensions();
    const IntegerArray & iw = cg.interpolationWidth;
    // iw.display("iw");
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
      MappedGrid & mg = cg[grid];
      const IntegerArray & dw = mg.discretizationWidth();
      
      // dw.display("dw");
      
      minDiscretizationWidth=min(minDiscretizationWidth,min(dw(R)));
      
      for( int grid2=0; grid2<cg.numberOfComponentGrids(); grid2++ )
      {
        if( grid!=grid2 )
  	  minInterpolationWidth=min( minInterpolationWidth,min(iw(R,grid,grid2)));
      }
    }
    if( minInterpolationWidth==INT_MAX ) minInterpolationWidth=minDiscretizationWidth;
    printF(" *** minDiscretizationWidth=%i, minInterpolationWidth=%i ****\n",minDiscretizationWidth,
	   minInterpolationWidth);

    const int maxOrderOfAccuracy=8;  // *************
    
    orderOfAccuracyInSpace=min(maxOrderOfAccuracy,minDiscretizationWidth-1,minInterpolationWidth-1);
    if( orderOfAccuracyInSpace%2 ==1 )
      orderOfAccuracyInSpace--;   // must be even
    
    orderOfAccuracyInTime =orderOfAccuracyInSpace;
    orderOfArtificialDissipation=orderOfAccuracyInSpace;
    
    printF("***Setting orderOfAccuracyInSpace=%i, orderOfAccuracyInTime=%i, orderOfArtificialDissipation=%i\n",
	   orderOfAccuracyInSpace,orderOfAccuracyInTime,orderOfArtificialDissipation);

    if( orderOfAccuracyInSpace>4 )
    {
      printF("***Setting useConservative=false by default for order of accuracy >4.\n");
      useConservative=false;
    }
    
  } // end compositeGrid

  timing(timeForInitialize)+=getCPU()-time0;

  // ********** By default do not solve for H in 3D ************
  if( cgp!=NULL && cgp->numberOfDimensions()==3 && method!=yee )
  {
    solveForMagneticField=false;
  }
  else if( cgp!=NULL && cgp->numberOfDimensions()==2 )
  {
    kz=0; // *wdh* 040626 
  }
  
  // These next arrays hold (eps,mu,c) in the case when they are constant on each grid but
  // may vary from grid to grid
  const int numberOfComponentGrids = mgp!=NULL ? 1 : cgp->numberOfComponentGrids();
  epsGrid.redim(numberOfComponentGrids); epsGrid=eps;
  muGrid.redim(numberOfComponentGrids);  muGrid=mu;
  cGrid.redim(numberOfComponentGrids);   cGrid=c;
  sigmaEGrid.redim(numberOfComponentGrids);  sigmaEGrid=0.;
  sigmaHGrid.redim(numberOfComponentGrids);  sigmaHGrid=0.;

  // Dispersive material parameters may vary from domain to domain
  std::vector<DispersiveMaterialParameters> & dmpVector = 
    dbase.get<std::vector<DispersiveMaterialParameters> >("dispersiveMaterialParameters");

  assert( cgp!=NULL );
  CompositeGrid & cg = *cgp;

  if( dmpVector.size()<cg.numberOfDomains() )
    dmpVector.resize(cg.numberOfDomains());

 // subtract out the incident field before apply NRBC's
  adjustFarFieldBoundariesForIncidentField.redim(numberOfComponentGrids);
  adjustFarFieldBoundariesForIncidentField=0;

  return 0;
}
