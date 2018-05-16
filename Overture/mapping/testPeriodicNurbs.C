//
// Test periodic NurbsMapping *new* *wdh* May 9, 2018
//

#include "Mapping.h"
#include "NurbsMapping.h"
#include "PlotStuff.h"
#include "MappingInformation.h"

int
main()
{
  ios::sync_with_stdio(); // Synchronize C++ and C I/O subsystems
  Index::setBoundsCheck(on);  //  Turn on A++ array bounds checking


  cout << "====== Test Periodic Nurbs =====" << endl;

  int option=0;
  // cout <<
  //   "Enter option: 0 = define a nurbs passing through given points (interpolate) \n"
  //   "              1 = test reparameterization \n"
  //   "              2 = test insertKnot, removeKnot \n"
  //   "              3 = test merge \n";
  // cin >> option;

  NurbsMapping nurbs;
  int rangeDimension=2;
  
/* ---
  cout << "nurbs.getRangeBound(0,0) =" << nurbs.getRangeBound(0,0) << endl;
  cout << "nurbs.getRangeBound(1,0) =" << nurbs.getRangeBound(1,0) << endl;
  cout << "nurbs.getRangeBound(0,1) =" << nurbs.getRangeBound(0,1) << endl;
  cout << "nurbs.getRangeBound(1,1) =" << nurbs.getRangeBound(1,1) << endl;
--- */  

  PlotStuff ps;
  GraphicsParameters params;

  int degree=3;

  if( option==0 )
  {
    // int numPts=11;
    int numPts=51;
    RealArray x0(numPts,2);
    // x0(0,0)=0.; x0(0,1)=0.;
    // x0(1,0)=1.; x0(1,1)=0.;
    // x0(2,0)=1.; x0(2,1)=1.;
    // x0(3,0)=0.; x0(3,1)=1.;
    // x0(4,0)=0.; x0(4,1)=0.;

    // put points on an ellipse
    real a=.75, b=.5;
    for( int i=0; i<numPts; i++ )
    {
      real theta=twoPi*i/(numPts-1.);
      x0(i,0) = a*cos(theta);
      x0(i,1) = b*sin(theta);
    }
    
    // nurbs.interpolate(x0);
    nurbs.interpolate(x0,0,Overture::nullRealArray(),degree,NurbsMapping::parameterizeByIndex);

    int numGridPoints=101;
    nurbs.setGridDimensions(0,numGridPoints);
    
    MappingInformation mapInfo;
    mapInfo.graphXInterface=&ps;
    
    int num=3;
    RealArray r(num,1), x(num,2), xr(num,2,1);
    r=0.; 
    r(1,0)=.5;
    r(2,0)=1.;
    nurbs.map( r,x,xr );
    for( int i=0; i<num; i++ )
    {
      printF(" nurbs (not periodic) : r=%.6e x=(%16.8e,%16.8e) xr=(%16.8e,%16.8e)\n",r(i,0),x(i,0),x(i,1),xr(i,0,0),xr(i,1,0));
    }

    if( true )
      nurbs.update(mapInfo);
    
    printF("testPeriodicNurbs: make a function periodic NURBS curve\n");
    NurbsMapping nurbs2;
    nurbs2.setIsPeriodic(0,Mapping::functionPeriodic);
  // void interpolate(const RealArray & x, 
  //       	   const int & option     = 0 ,
  //       	   RealArray & parameterization  =Overture::nullRealArray(),
  //                  int degree = 3,
  //                  ParameterizationTypeEnum parameterizationType=parameterizeByChordLength,
  //                  int numberOfGhostPoints=0 );

    // TWO KEY THINGS to get a good periodic curve: 
    //    (1) Use Equally spaced knots when interpolating
    //    (2) parameterizeByIndex
    nurbs2.setInterpolateKnotsOption( NurbsMapping::interpolateWithEquallySpacedKnots );

    // nurbs2.interpolate(x0);
    nurbs2.interpolate(x0,0,Overture::nullRealArray(),degree,NurbsMapping::parameterizeByIndex);
    nurbs2.setGridDimensions(0,numGridPoints);
    nurbs2.map( r,x,xr );
    for( int i=0; i<num; i++ )
    {
      printF(" nurbs2 (periodic)   : r=%.6e x=(%16.8e,%16.8e) xr=(%16.8e,%16.8e)\n",r(i,0),x(i,0),x(i,1),xr(i,0,0),xr(i,1,0));
    }

    nurbs2.update(mapInfo);
    
    if( false )
    {
      // Did not get this to work very well:
      printF("testPeriodicNurbs: Wrap the points to make periodic\n");
      int extra=6; // wrap this many extra points
      RealArray xp(numPts+2*extra,2);

      // put points on an ellipse
      for( int i=-extra; i<numPts+extra; i++ )
      {
        real theta=twoPi*i/(numPts-1.);
        xp(i+extra,0) = a*cos(theta);
        xp(i+extra,1) = b*sin(theta);
      }

      NurbsMapping nurbs3; 
      // nurbs3.setIsPeriodic(0,Mapping::functionPeriodic);
      // nurbs3.interpolate(xp);
      nurbs3.interpolate(xp,0,Overture::nullRealArray(),degree,NurbsMapping::parameterizeByIndex);
      nurbs3.setGridDimensions(0,numGridPoints);

      RealArray x2(1,2),r2(1,1);
      Range all;
      x2(0,all)=xp(extra,all); r2=0.;
      nurbs3.inverseMap(x2,r2);
      printF(" r2=%.6e\n",r2(0,0));

      // real ra = extra/(numPts-1.);
      if( r2(0,0)>.5 )
        r2(0,0)=1.-r2(0,0);
    
      real ra = r2(0,0);
      real rb=1.-ra;
      nurbs3.reparameterize(ra,rb);
    
      nurbs3.map( r,x,xr );
      for( int i=0; i<num; i++ )
      {
        printF(" nurbs3 (wrapped)   : r=%.6e x=(%16.8e,%16.8e) xr=(%16.8e,%16.8e)\n",r(i,0),x(i,0),x(i,1),xr(i,0,0),xr(i,1,0));
      }

      nurbs3.update(mapInfo);
    }
    
    // num=101;
    // r.redim(num,1), x.redim(num,2), xr.redim(num,2,1);
    // for( int i=0; i<num; i++ )
    // {
    //   r(i,0) = i*1./(num-1);
    // }
    // nurbs.map(r,x,xr);
    // RealArray x2(num,2), xr2(num,2,1);
    // nurbs2.map(r,x2,xr2);
    
    // RealArray xe(num,2);
    // for( int i=0; i<numPts; i++ )
    // {
    //   real theta=twoPi*i/(numPts-1.);
    //   xe(i,0) = a*cos(theta);
    //   xe(i,1) = b*sin(theta);
    // }

    // PlotIt::plot(ps,nurbs);
    // GraphicsParameters params;
    // bool plotControlPoints=true;
    // nurbs.plot(ps,params,plotControlPoints);
  }
//   else if( option==1 )
//   {
  
//     PlotIt::plot(ps,nurbs);
//     nurbs.reparameterize(.1,.9, .5,1.);
//     PlotIt::plot(ps,nurbs);
//     nurbs.checkMapping();
//     nurbs.reparameterize(1.,0.);
//     PlotIt::plot(ps,nurbs);
//     nurbs.checkMapping();
//   }
//   else if( option==2 )
//   {
//     // test insertKnot, removeKnot    
 
//     int p1=3;
//     int n1=p1+1;
//     int m1=n1+p1+1;
//     rangeDimension=2;

//     // knots are clamped
//     realArray uKnot,cPoint;
  
//     uKnot.redim(m1+1);
//     uKnot(0)=0.; uKnot(1)=0.; uKnot(2)=0.; uKnot(3)=0.;
//     uKnot(4)=.5;
//     uKnot(m1-3)=1.; uKnot(m1-2)=1.; uKnot(m1-1)=1.; uKnot(m1)=1.; 

//     // control points (holds weight in last position)
//     cPoint.redim(n1+1,3);   
//     cPoint(0,0)=0.;  cPoint(0,1)=0.; cPoint(0,2)=1.;
//     cPoint(1,0)=.25; cPoint(1,1)=.7; cPoint(1,2)=1.;

//     cPoint(2,0)=.75; cPoint(2,1)=.7; cPoint(2,2)=1.;

//     cPoint(n1-1,0)=.75; cPoint(n1-1,1)=1.; cPoint(n1-1,2)=1.;
//     cPoint(n1,0)=1.;  cPoint(n1,1)=0.; cPoint(n1,2)=1.;

//     nurbs.specify( m1,n1,p1,uKnot,cPoint,rangeDimension);
    
//     params.set(GI_TOP_LABEL,"original nurb");  // set title
//     PlotIt::plot(ps,nurbs,params);

//     nurbs.insertKnot(uKnot(4));
//     params.set(GI_TOP_LABEL,"knot inserted"); 
//     PlotIt::plot(ps,nurbs,params);
  
//     int numberRemoved;
//     nurbs.removeKnot(5,1,numberRemoved);
//     params.set(GI_TOP_LABEL,"knot removed"); 
//     PlotIt::plot(ps,nurbs,params);
//   }
//   else if( option==3 )
//   {
//     // test merge
  
//     realArray uKnot,cPoint;

//     NurbsMapping nurbs2;
//     int p1=2;
//     int n1=2;
//     int m1=n1+p1+1;

//     uKnot.redim(m1+1);
//     uKnot(0)=0.; uKnot(1)=0.; uKnot(2)=0.; 
//     uKnot(m1-2)=1.; uKnot(m1-1)=1.; uKnot(m1)=1.; 

//     // control points (holds weight in last position)
//     cPoint.redim(n1+1,3);   
//     cPoint(0,0)=1.; cPoint(0,1)=0.;  cPoint(0,2)=1.;
//     cPoint(1,0)=1.; cPoint(1,1)=1.;  cPoint(1,2)=SQRT(2.)/2.;
//     cPoint(2,0)=0.; cPoint(2,1)=1.;  cPoint(2,2)=1.;

//     nurbs2.specify( m1,n1,p1,uKnot,cPoint,rangeDimension);
//     params.set(GI_TOP_LABEL,"nurbs 2, 90 degree arc");
//     PlotIt::plot(ps,nurbs2,params);

//     // another 90 degree arc
//     NurbsMapping nurbs3;
//     p1=2;
//     n1=2;
//     m1=n1+p1+1;

//     uKnot.redim(m1+1);
//     uKnot(0)=0.; uKnot(1)=0.; uKnot(2)=0.; 
//     uKnot(m1-2)=1.; uKnot(m1-1)=1.; uKnot(m1)=1.; 

//     // control points (holds weight in last position)
//     cPoint.redim(n1+1,3);   
//     cPoint(0,0)=0.;  cPoint(0,1)=1.;  cPoint(0,2)=1.;
//     cPoint(1,0)=-1.; cPoint(1,1)=1.;  cPoint(1,2)=SQRT(2.)/2.;
//     cPoint(2,0)=-1.; cPoint(2,1)=0.;  cPoint(2,2)=1.;

//     nurbs3.specify( m1,n1,p1,uKnot,cPoint,rangeDimension);
//     params.set(GI_TOP_LABEL,"nurbs 3, 90 degree arc");
//     PlotIt::plot(ps,nurbs3,params);

//     printf("check curve \n");
//     nurbs3.checkMapping();

//     nurbs3.merge(nurbs2);
//     params.set(GI_TOP_LABEL,"merged nurb, 180 degree arc");

// /* ---
//   Mapping::debug=63;
//   realArray r(1,1),x(1,2),t(1,1);
//   r=.5;
//   nurbs3.map(r,x);
//   nurbs3.inverseMap(x,t);
//   r.display("Here is r");
//   x.display("here is x");
//   t.display(" inverseMap(map(r)) ");
//   Mapping::debug=0;
// ---- */

//     PlotIt::plot(ps,nurbs3,params);
// //  printf("check merged curve \n");
// //  nurbs3.checkMapping();
  

//  }
  
  return 0;
}
