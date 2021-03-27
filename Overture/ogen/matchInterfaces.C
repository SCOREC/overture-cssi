// ====================================================================================
//     Locate matching interfaces
// ====================================================================================

// #define BOUNDS_CHECK

#include "Overture.h"
#include "display.h"
#include "ParallelUtility.h"
#include "InterfaceInfo.h"
//#include "ParallelGridUtility.h"
#include "Ogen.h"

static const int ISneededPoint = CompositeGrid::ISreservedBit2;  // from Cgsh.h


#define FOR_3D(i1,i2,i3,I1,I2,I3) \
int I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase();  \
int I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(); \
for(int i3=I3Base; i3<=I3Bound; i3++) \
for(int i2=I2Base; i2<=I2Bound; i2++) \
for(int i1=I1Base; i1<=I1Bound; i1++)

#define FOR_3(i1,i2,i3,I1,I2,I3) \
I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase();  \
I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(); \
for(int i3=I3Base; i3<=I3Bound; i3++) \
for(int i2=I2Base; i2<=I2Bound; i2++) \
for(int i1=I1Base; i1<=I1Bound; i1++)

#define FOR_3IJD(i1,i2,i3,I1,I2,I3,j1,j2,j3,J1,J2,J3) \
int I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase();  \
int I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(); \
int J1Base =J1.getBase(),   J2Base =J2.getBase(),  J3Base =J3.getBase();  \
for(int i3=I3Base,j3=J3Base; i3<=I3Bound; i3++,j3++) \
for(int i2=I2Base,j2=J2Base; i2<=I2Bound; i2++,j2++) \
for(int i1=I1Base,j1=J1Base; i1<=I1Bound; i1++,j1++)


#define FOR_4IJD(i1,i2,i3,i4,I1,I2,I3,I4,j1,j2,j3,j4,J1,J2,J3,J4) \
int I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase(),  I4Base =I4.getBase();  \
int I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(), I4Bound=I4.getBound(); \
int J1Base =J1.getBase(),   J2Base =J2.getBase(),  J3Base =J3.getBase(),  J4Base =J4.getBase();  \
for(i4=I4Base,j4=J4Base; i4<=I4Bound; i4++,j4++) \
for(i3=I3Base,j3=J3Base; i3<=I3Bound; i3++,j3++) \
for(i2=I2Base,j2=J2Base; i2<=I2Bound; i2++,j2++) \
for(i1=I1Base,j1=J1Base; i1<=I1Bound; i1++,j1++)

#define FOR_4IJ(i1,i2,i3,i4,I1,I2,I3,I4,j1,j2,j3,j4,J1,J2,J3,J4) \
I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase(),  I4Base =I4.getBase();  \
I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(), I4Bound=I4.getBound(); \
J1Base =J1.getBase(),   J2Base =J2.getBase(),  J3Base =J3.getBase(),  J4Base =J4.getBase();  \
for(i4=I4Base,j4=J4Base; i4<=I4Bound; i4++,j4++) \
for(i3=I3Base,j3=J3Base; i3<=I3Bound; i3++,j3++) \
for(i2=I2Base,j2=J2Base; i2<=I2Bound; i2++,j2++) \
for(i1=I1Base,j1=J1Base; i1<=I1Bound; i1++,j1++)

#define FOR_4IJD_WITH_STRIDE(i0,i1,i2,i3,indexi,j0,j1,j2,j3,indexj)\
      for( int i3=indexi[3][0],j3=indexj[3][0]; i3<=indexi[3][1]; i3+=indexi[3][2], j3+=indexj[3][2] )\
      for( int i2=indexi[2][0],j2=indexj[2][0]; i2<=indexi[2][1]; i2+=indexi[2][2], j2+=indexj[2][2] )\
      for( int i1=indexi[1][0],j1=indexj[1][0]; i1<=indexi[1][1]; i1+=indexi[1][2], j1+=indexj[1][2] )\
      for( int i0=indexi[0][0],j0=indexj[0][0]; i0<=indexi[0][1]; i0+=indexi[0][2], j0+=indexj[0][2] )



// =====================================================================================================
/// \brief  Find the matching interfaces between domains.
// =====================================================================================================
void 
matchInterfaces( CompositeGrid & cg, std::vector<InterfaceInfo> & interfaceInfo )
{
  // real time0=getCPU();
  const int interfaceShareFlag=100;   // INTERFACE SHARE FLAGS ARE GREATER THAN OR EQUAL TO THIS 

  int debug=1;
  FILE *debugFile  = stdout;   // do this for now 
  FILE *pDebugFile = stdout;   // do this for now 

  if( true || debug & 4 )
    printF("matchInterfaces....\n");
  

  // ************** THIS ROUTINE SHOULD USE THE ONE IN cg/mp/src/assignInterfaceConditions.C ****************

  // assert( cgp!=NULL );
  // CompositeGrid & cg= *cgp;
  const int numberOfDimensions = cg.numberOfDimensions();
  const int numberOfComponentGrids = cg.numberOfComponentGrids();

  // int bcOrderOfAccuracy=orderOfAccuracyInSpace;
  // if( method==sosup && orderOfAccuracyInSpace==6 )
  // {
  //   // NOTE: for now apply 4th order BC's for sosup order 6
  //   bcOrderOfAccuracy=4;
  // }

  Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
  Index Jv[3], &J1=Jv[0], &J2=Jv[1], &J3=Jv[2];

  // --- The matched array is used to avoid checking faces that have already been matched ---
  //       matched(side,axis,grid)=1 if this face has been natched to another.
  IntegerArray matched(2,cg.numberOfDimensions(),cg.numberOfComponentGrids());
  matched=0;  
  
  int grid1;
  for( grid1=0; grid1<cg.numberOfComponentGrids(); grid1++ )
  {

    MappedGrid & mg1 = cg[grid1];
    const IntegerArray & bc1 = mg1.boundaryCondition();
    const IntegerArray & share1 = mg1.sharedBoundaryFlag();

    // check for interface boundary conditions
    for( int dir1=0; dir1<mg1.numberOfDimensions(); dir1++ )
    {
      for( int side1=0; side1<=1; side1++ )
      {
        if( share1(side1,dir1) >= interfaceShareFlag && matched(side1,dir1,grid1)==0 )
        {
          bool interfaceFound=false;
          for( int grid2=grid1+1; grid2<cg.numberOfComponentGrids(); grid2++ ) // only check higher numbered grids
          {
            // When there are multiple domains, matching interfaces grids should not be in the same domain:
            if( cg.numberOfDomains()>1 && cg.domainNumber(grid1)==cg.domainNumber(grid2) )
              continue;

            MappedGrid & mg2 = cg[grid2];
            const IntegerArray & bc2 = mg2.boundaryCondition();
            const IntegerArray & share2 = mg2.sharedBoundaryFlag();

            for( int dir2=0; dir2<mg2.numberOfDimensions(); dir2++ )
            {
              for( int side2=0; side2<=1; side2++ )
              {
                if( share2(side2,dir2) >= interfaceShareFlag &&
                    share1(side1,dir1)==share2(side2,dir2) && 
                    matched(side2,dir2,grid2)==0 && 
                    matched(side1,dir1,grid1)==0 )  
                {
                  // ********************************************
                  // **** potential interface found *************
                  // ********************************************

                  // ---- Check that the interfaces are close in space ---- *wdh* 2017/05/06 
                  Mapping & map1 = mg1.mapping().getMapping();
                  Mapping & map2 = mg2.mapping().getMapping();
                  
                  RealArray r1(1,3), x1(1,3), r2(1,3), x2(1,3); x1=0.; x2=0.;
                  r1=.5; r1(0,dir1)=side1;  // r-coords for a point on the middle of the face of grid1
                  map1.mapS(r1,x1);           // x coords on grid 1
                  r2=-1.;
                  map2.inverseMapS(x1,r2);    // r-coords on grid 2
                  map2.mapS(r2,x2);           // x-coords on grid 2
                  real rDist=fabs( r2(0,dir2)-side2 );
                  real xDist=fabs( x1(0,0)-x2(0,0) ) + fabs( x1(0,1)-x2(0,1) ) + fabs( x1(0,2)-x2(0,2) );
                  const RealArray & bb = map1.getBoundingBox();
                  real xNorm = fabs(bb(1,0)-bb(0,0)) + fabs(bb(1,1)-bb(0,1));
                  if( numberOfDimensions==3 ) xNorm +=  fabs(bb(1,2)-bb(0,2));

                  xDist /=max(REAL_MIN*100,xNorm);  // normalize x-distance 
                  
                  if( debug & 2 )
                  {
                    printF("--matchInterfaces-- POTENTIAL Interface: %s=(grid1,side,dir1)=(%i,%i,%i) "
                           "AND %s=(grid2,side2,dir2)=(%i,%i,%i) share=%i\n",
                           (const char*)mg1.getName(),grid1,side1,dir1,
                           (const char*)mg2.getName(),grid2,side2,dir2,share1(side1,dir1));
                    printF(" grid1: r1=[%8.2e,%8.2e,%8.2e] x1=[%8.2e,%8.2e,%8.2e]\n"
                           " grid2: r2=[%8.2e,%8.2e,%8.2e] x2=[%8.2e,%8.2e,%8.2e] .. rDist=%8.2e xDist=%8.2e, xNorm=%8.2e\n",
                           r1(0,0),r1(0,1),r1(0,2), x1(0,0),x1(0,1),x1(0,2), 
                           r2(0,0),r2(0,1),r2(0,2), x2(0,0),x2(0,1),x2(0,2), rDist, xDist, xNorm);
                  }
                  // --- tolerances for surfaces to match ---
                  const real rTol=1.e-2;
                  const real xTol=1.e-5;
                   
                  if( rDist<rTol && xDist < xTol )
                  {
                    if( interfaceFound )
                    {
                      // -- ERROR a previous matching interface has already been found! ---

                      printF("matchInterfaces:ERROR: A second valid interface match has been found:\n"
                             " (grid1,side1,dir1,share1)=(%i,%i,%i,%i) also matches"
                             " (grid2,side2,dir2,share2)=(%i,%i,%i,%i).\n"
                             " Matching interface grids should be marked with distinct share values.\n",
                             grid1,side1,dir1,share1(side1,dir1), grid2,side2,dir2,share2(side2,dir2));
                      OV_ABORT("error");
                    }

                    // ---- ADD an interface to the list ----
                    interfaceFound=true;
                    matched(side1,dir1,grid1)=1; // indicates these faces have been matched
                    matched(side2,dir2,grid2)=1;
                    if( debug & 1 )
                    {
                      printF("--matchInterfaces-- Interface found: %s=(grid1,side,dir1)=(%i,%i,%i) matches "
                             "%s=(grid2,side2,dir2)=(%i,%i,%i) share=%i\n",
                             (const char*)mg1.getName(),grid1,side1,dir1,
                             (const char*)mg2.getName(),grid2,side2,dir2,share1(side1,dir1));
                    }
                    // printF("interfaceInfo.size() = %d\n",interfaceInfo.size());
                    // int interfaceInfoSize=interfaceInfo.size();
                    // interfaceInfo.resize(interfaceInfoSize+1);
                    // interfaceInfo[interfaceInfoSize]=InterfaceInfo(grid1,side1,dir1, grid2,side2,dir2);
                    
                    interfaceInfo.push_back(InterfaceInfo(grid1,side1,dir1, grid2,side2,dir2));   

                    // InterfaceInfo & info12 = *(new InterfaceInfo(grid1,side1,dir1, grid2,side2,dir2));
                    // interfaceInfo.push_back(info12);
                  
                    // ***** Check that the valid points match on the interface ******
                    //   *wdh* 2015/08/11
                    const intArray & mask1 = mg1.mask();
                    const intArray & mask2 = mg2.mask();
                    OV_GET_SERIAL_ARRAY_CONST(int,mask1,mask1Local);
                    OV_GET_SERIAL_ARRAY_CONST(int,mask2,mask2Local);
                    // *wdh* Nov 29, 2020
                    // getBoundaryIndex(mg1.dimension(),side1,dir1,I1,I2,I3);
                    // getBoundaryIndex(mg2.dimension(),side2,dir2,J1,J2,J3);
                    // Fixed: check boundary points *wdh* Nov 29, 2020
                    getBoundaryIndex(mg1.gridIndexRange(),side1,dir1,I1,I2,I3);
                    getBoundaryIndex(mg2.gridIndexRange(),side2,dir2,J1,J2,J3);

                    if( debug & 8 )
                    {
                      fprintf(debugFile,"--matchInterfaces-- Interface found: %s=(grid1,side,dir1)=(%i,%i,%i) matches "
                           "%s=(grid2,side2,dir2)=(%i,%i,%i) share=%i\n",
                           (const char*)mg1.getName(),grid1,side1,dir1,
                           (const char*)mg2.getName(),grid2,side2,dir2,share1(side1,dir1));

                      displayMask(mask1,"mask1",debugFile);
                      displayMask(mask2,"mask2",debugFile);

                      fflush(debugFile);
                    }

                    
                    // -- Check that interfaces have the same number of grid points -- 
                    const int d1=cg.domainNumber(grid1), d2=cg.domainNumber(grid2);
                    
                    for( int dir=1; dir<mg1.numberOfDimensions(); dir++ )
                    {
                      int dir1p = (dir1+dir) % mg1.numberOfDimensions();
                      int dir2p = (dir2+dir) % mg2.numberOfDimensions();
                      if( Iv[dir1p].getLength()!=Jv[dir2p].getLength() )
                      {
                        printF("matchInterfaces::ERROR: The number of grid points on the two interfaces do not match\n"
                               " (d1,grid1,side1,dir1,bc1)=(%i,%i,%i,%i,%i) Iv=[%i,%i][%i,%i][%i,%i]\n"
                               " (d2,grid2,side2,dir2,bc2)=(%i,%i,%i,%i,%i) Jv=[%i,%i][%i,%i][%i,%i]\n",
                               d1,grid1,side1,dir1,mg1.boundaryCondition(side1,dir1),
                               I1.getBase(),I1.getBound(),I2.getBase(),I2.getBound(),I3.getBase(),I3.getBound(),  
                               d2,grid2,side2,dir2,mg2.boundaryCondition(side2,dir2),
                               J1.getBase(),J1.getBound(),J2.getBase(),J2.getBound(),J3.getBase(),J3.getBound());
                        printF("grid names are [%s] and [%s]\n",(const char*)mg1.getName(),(const char*)mg2.getName());
                        OV_ABORT("error");
                      }
                    }
                                            
                      

                    
                    // ------ check masks moved -----
                    // int includeGhost=1;
                    // bool ok1 = ParallelUtility::getLocalArrayBounds(mask1,mask1Local,I1,I2,I3,includeGhost);
                    // bool ok2 = ParallelUtility::getLocalArrayBounds(mask2,mask2Local,J1,J2,J3,includeGhost);

                    // // This next check does not always work in parallel -- fix me -- could add check when we create 
                    // // copies of the mask array *fix me*
                    // if( ok1 && ok2 &&
                    //     I1.getLength()==J1.getLength() &&  
                    //     I2.getLength()==J2.getLength() &&
                    //     I3.getLength()==J3.getLength() 
                    //   )
                    // {
                    //   // ---- The check for parallel is now done in getLocalInterfaceArrays macro ---
                      
                    //   // int maskDiff = max(abs(mask1Local(I1,I2,I3)-mask2Local(J1,J2,J3)));

                    //   // Check that both masks are discretization points (ignore isNeeded etc.) *wdh* Nov 29, 2020
                    //   int maskDiff = max( (mask1Local(I1,I2,I3) & MappedGrid::ISdiscretizationPoint) -
                    //                       (mask2Local(J1,J2,J3) & MappedGrid::ISdiscretizationPoint) );

                    //   if( maskDiff > 0 )
                    //   {
                    //     printF("       ERROR: the mask arrays do not match on the interface. This is currently required.\n");
                    //     printF("       Try re-generating the grid with more lines in the normal direction, this sometimes fixes this problem.\n");
                    //     printF("       See debug file for more info if debug &2 \n");
                    //     if( debug & 2 )
                    //     {
                    //       displayMask(mask1Local(I1,I2,I3),"mask1Local(I1,I2,I3)",pDebugFile);
                    //       displayMask(mask2Local(J1,J2,J3),"mask2Local(J1,J2,J3)",pDebugFile);
                    //       IntegerArray md(I1,I2,I3);
                    //       md = ( (mask1Local(I1,I2,I3) & MappedGrid::ISdiscretizationPoint) -
                    //              (mask2Local(J1,J2,J3) & MappedGrid::ISdiscretizationPoint) );
                    //       where( md>0 )
                    //       {
                    //         md=1;
                    //       }
                          
                    //      ::display(md,"difference)",pDebugFile,"%2i ");

                    //       fprintf(pDebugFile,"       INFO: interface grid points match: mask arrays agree on the interface.\n\n");
                    //     }
                    //     // OV_ABORT("ERROR");
                    //   }
                    //   else
                    //   {
                    //     printF("       INFO: interface grid points match: mask arrays agree on the interface.\n\n");
                    //     if( debug & 2 )
                    //     {
                    //       displayMask(mask1Local(I1,I2,I3),"mask1Local(I1,I2,I3)",pDebugFile);
                    //       displayMask(mask2Local(J1,J2,J3),"mask2Local(J1,J2,J3)",pDebugFile);

                    //       fprintf(pDebugFile,"       INFO: interface grid points match: mask arrays agree on the interface.\n\n");
                    //     }
                        
                    //   }
                    // }

                  } // end if rDist < rTol and xDist < xTol 
                  
                }
              }
            }
          } // end for grid 2
          if( !interfaceFound )
          { // if no interface was found, check that we haven't found a match already: 
            for( int inter=0; inter < interfaceInfo.size(); inter++ )
            {
              InterfaceInfo & interface = interfaceInfo[inter]; 
              if( (grid1==interface.grid1 && side1==interface.side1 && dir1==interface.dir1) ||
                  (grid1==interface.grid2 && side1==interface.side2 && dir1==interface.dir2) )
              {
                interfaceFound=true;
                break;
              }
            }
            if( !interfaceFound )
            {
              printF("initializeInterfaces:ERROR: No matching interface found for (grid1,side,dir1)=(%i,%i,%i).\n",
                     grid1,side1,dir1);
              OV_ABORT("error");
            }
          }
        }
      }
    }
  }



  // timing(timeForInterfaceBC)+=getCPU()-time0;
}

// =====================================================================================================
/// \brief Check that the mask values match on interfaces
// =====================================================================================================
int 
checkMatchingMasksOnInterfaces( CompositeGrid & cg, std::vector<InterfaceInfo> & interfaceInfo )
{
  // real time0=getCPU();

  int debug=1;
  FILE *debugFile  = stdout;   // do this for now 
  FILE *pDebugFile = stdout;   // do this for now 

  if( true || debug & 4 )
    printF("checkMatchingMasksOnInterfaces....\n");
  

  const int numberOfDimensions = cg.numberOfDimensions();
  const int numberOfComponentGrids = cg.numberOfComponentGrids();

  
  Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
  Index Jv[3], &J1=Jv[0], &J2=Jv[1], &J3=Jv[2];

  //  ----------------------------------------------------------------------------
  //  ------------------------- loop over interfaces -----------------------------
  //  ----------------------------------------------------------------------------
  bool masksMatch=true; // set to false if some interface masks do not match
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

    // ***** Check that the valid points match on the interface ******
    
    const intArray & mask1 = mg1.mask();
    const intArray & mask2 = mg2.mask();
    OV_GET_SERIAL_ARRAY_CONST(int,mask1,mask1Local);
    OV_GET_SERIAL_ARRAY_CONST(int,mask2,mask2Local);
    
    getBoundaryIndex(mg1.gridIndexRange(),side1,dir1,I1,I2,I3);
    getBoundaryIndex(mg2.gridIndexRange(),side2,dir2,J1,J2,J3);

    if( debug & 2 )
    {
      fprintf(debugFile,"--checkMatchingMasksOnInterfaces-- Interface : %s=(grid1,side,dir1)=(%i,%i,%i) matches "
           "%s=(grid2,side2,dir2)=(%i,%i,%i) share=%i\n",
           (const char*)mg1.getName(),grid1,side1,dir1,
           (const char*)mg2.getName(),grid2,side2,dir2,share1(side1,dir1));
    }

    if( debug & 8 )
    {

      displayMask(mask1,"mask1",debugFile);
      displayMask(mask2,"mask2",debugFile);

      fflush(debugFile);
    }

    if( true)
    {
      // -- Check that interfaces have the same number of grid points -- 
      const int d1=cg.domainNumber(grid1), d2=cg.domainNumber(grid2);
      
      for( int dir=1; dir<mg1.numberOfDimensions(); dir++ )
      {
        int dir1p = (dir1+dir) % mg1.numberOfDimensions();
        int dir2p = (dir2+dir) % mg2.numberOfDimensions();
        if( Iv[dir1p].getLength()!=Jv[dir2p].getLength() )
        {
          printF("checkMatchingMasksOnInterfaces::ERROR: The number of grid points on the two interfaces do not match\n"
                 " (d1,grid1,side1,dir1,bc1)=(%i,%i,%i,%i,%i) Iv=[%i,%i][%i,%i][%i,%i]\n"
                 " (d2,grid2,side2,dir2,bc2)=(%i,%i,%i,%i,%i) Jv=[%i,%i][%i,%i][%i,%i]\n",
                 d1,grid1,side1,dir1,mg1.boundaryCondition(side1,dir1),
                 I1.getBase(),I1.getBound(),I2.getBase(),I2.getBound(),I3.getBase(),I3.getBound(),  
                 d2,grid2,side2,dir2,mg2.boundaryCondition(side2,dir2),
                 J1.getBase(),J1.getBound(),J2.getBase(),J2.getBound(),J3.getBase(),J3.getBound());
          printF("grid names are [%s] and [%s]\n",(const char*)mg1.getName(),(const char*)mg2.getName());
          OV_ABORT("error");
        }
      }

           
    }

    int includeGhost=1;
    bool ok1 = ParallelUtility::getLocalArrayBounds(mask1,mask1Local,I1,I2,I3,includeGhost);
    bool ok2 = ParallelUtility::getLocalArrayBounds(mask2,mask2Local,J1,J2,J3,includeGhost);

    // This next check does not always work in parallel -- fix me -- could add check when we create 
    // copies of the mask array *fix me*
    if( ok1 && ok2 &&
        I1.getLength()==J1.getLength() &&  
        I2.getLength()==J2.getLength() &&
        I3.getLength()==J3.getLength() 
       )
     {
       // ****** FIX ME ******

      // ---- The check for parallel is now done in (MAXWELL) getLocalInterfaceArrays macro ---
      
      // int maskDiff = max(abs(mask1Local(I1,I2,I3)-mask2Local(J1,J2,J3)));

      // Check that both masks are discretization points (ignore isNeeded etc.) *wdh* Nov 29, 2020
      int maskDiff = max( (mask1Local(I1,I2,I3) & MappedGrid::ISdiscretizationPoint) -
                          (mask2Local(J1,J2,J3) & MappedGrid::ISdiscretizationPoint) );

      if( maskDiff > 0 )
      {
        masksMatch=false;
        printF("       ERROR: the mask arrays do not match on the interface. This is currently required.\n");
        printF("       Try re-generating the grid with more lines in the normal direction, this sometimes fixes this problem.\n");
        printF("       See debug file for more info if debug &2 \n");
        if( debug & 2 )
        {
          displayMask(mask1Local(I1,I2,I3),"mask1Local(I1,I2,I3)",pDebugFile);
          displayMask(mask2Local(J1,J2,J3),"mask2Local(J1,J2,J3)",pDebugFile);
          IntegerArray md(I1,I2,I3);
          md = ( (mask1Local(I1,I2,I3) & MappedGrid::ISdiscretizationPoint) -
                 (mask2Local(J1,J2,J3) & MappedGrid::ISdiscretizationPoint) );
          where( md>0 )
          {
            md=1;
          }
          
         ::display(md,"difference)",pDebugFile,"%2i ");

          fprintf(pDebugFile,"       INFO: interface grid points match: mask arrays agree on the interface.\n\n");
        }
        // OV_ABORT("ERROR");
      }
      else
      {
         printF("       INFO: interface grid points match: mask arrays agree on the interface.\n\n");
         if( debug & 8 )
         {
           displayMask(mask1Local(I1,I2,I3),"mask1Local(I1,I2,I3)",pDebugFile);
           displayMask(mask2Local(J1,J2,J3),"mask2Local(J1,J2,J3)",pDebugFile);

           fprintf(pDebugFile,"       INFO: interface grid points match: mask arrays agree on the interface.\n\n");
         }
       }
     } // end if ok 

  }
  if( masksMatch )
  {
    printF("\n =============  SUCCESS : ALL INTERFACES HAVE MATCHING MASKS =============\n\n");
  }
  else
  {
  printF("\n =============  **FAILURE**: SOME INTERFACES DO NOT HAVE MATCHING MASKS =============\n\n");
  }
  // timing(timeForInterfaceBC)+=getCPU()-time0;
  return 0;
}


// ** MOVED TO GRID-STATISTICS ***
// // ========================================================================================================
// /// \brief Check for "tall" cells at interfaces, cells that are taller in the normal direction. This may
// ///   cause instabilities for some PDE solvers.
// /// \param interfaceInfo (input) : list of interfaces to check
// /// \param tallCellRatioBound (input) : print an error message if the tall-cell ratio exceeds this value (default =1.3)
// // ========================================================================================================
// real
// checkForTallCells( CompositeGrid & cg, std::vector<InterfaceInfo> & interfaceInfo, real tallCellRatioBound=-1 )
// {
//   if( interfaceInfo.size()==0 ) return 1.;

//   if( tallCellRatioBound<0 )
//       tallCellRatioBound=1.3;   // default worst case allowed

//   // ----- check for tall cells at interface -- these may cause instabilities for order=4 ----
//   // *wdh* March 3, 2021
//   //              |
//   //              +---------+---
//   //            /\|         |
//   //     t-dist   | n-dist  |        tall-cell-ratio = n-dist/t-dist 
//   //            \/+<------->+---
//   //              |         |
//   //              |         |
//   //              +---------+---
//   //              |         |
//   //              |         |
//   //              +---------+---
//   //              |
//   //              ^
//   //           interface
//   // 
  
//   printF("--- checkForTallCells: CHECK FOR TALL CELLS AT THE INTERFACE--- \n");

//   // ---- check for cells with normal-dist/tangential-dist > tallCellRatioBound ----
      
//   real maxTallCellRatio=0.;     // worst tall-cell ratio
//   int numTallTotal     =0;      // counts total number of tall cells exceeding bound
//   int numPointsTotal   =0;      // total interface cells checked
  
//   const int numberOfDimensions = cg.numberOfDimensions();
//   Index I1,I2,I3;
//   for( int inter=0; inter < interfaceInfo.size(); inter++ )
//   {
//     InterfaceInfo & interface = interfaceInfo[inter]; 
//     for( int iside=0; iside<=1; iside++ )// two sides of the interface
//     {
//       int grid, side, dir;
//       if( iside==0 )
//       {
//         grid=interface.grid1, side=interface.side1, dir=interface.dir1;
//       }
//       else
//       {
//         grid=interface.grid2, side=interface.side2, dir=interface.dir2;
//       }
//       MappedGrid & mg = cg[grid];
  
//       intArray & mask = mg.mask();
//       OV_GET_SERIAL_ARRAY(int,mask,maskLocal);

//       mg.update( MappedGrid::THEvertex );  // *** FIX ME **********************************

//       OV_GET_SERIAL_ARRAY(real,mg.vertex(),xLocal);
    
//       const int extra=0; 
//       getBoundaryIndex(mg.gridIndexRange(),side,dir,I1,I2,I3,extra);

//       int includeGhost=0;
//       bool ok1 = ParallelUtility::getLocalArrayBounds(mask,maskLocal,I1,I2,I3,includeGhost);
//       if( ok1  )
//       {
//         int isv[3], &is1 = isv[0], &is2=isv[1], &is3=isv[2]; // normal index shift
//         isv[0]=0; isv[1]=0; isv[2]=0;
//         isv[dir]= 1-2*side;

//         int dirp1 = (dir+1) % numberOfDimensions;
//         int itv[3], &it1 = itv[0], &it2=itv[1], &it3=itv[2]; // tangential index shift
//         itv[0]=0; itv[1]=0; itv[2]=0;
//         itv[dirp1]= 1;

      
//         real tallCellRatio=0.;  // worst TCR on this grid 
//         int numTall=0;          // counts tall cells on this grid exceeding bound
//         int numPoints=0;        // counts total points checked on this grid 
//         if( numberOfDimensions==2 )
//         {
//           // tallCellRatio = normal-dist / tangential-dist 
//           FOR_3D(i1,i2,i3,I1,I2,I3) 
//           {
//             if( mask(i1,i2,i3)>0 )
//             {
//               // distance "normal" to the boundary for this point
//               real nDist = sqrt( SQR( xLocal(i1+is1,i2+is2,i3+is3,0) - xLocal(i1,i2,i3,0) ) +
//                                  SQR( xLocal(i1+is1,i2+is2,i3+is3,1) - xLocal(i1,i2,i3,1) )  );
//               // tangential distance 
//               real tDist = sqrt( SQR( xLocal(i1+it1,i2+it2,i3+it3,0) - xLocal(i1,i2,i3,0) ) +
//                                  SQR( xLocal(i1+it1,i2+it2,i3+it3,1) - xLocal(i1,i2,i3,1) )  );
//               tDist = max( tDist, REAL_MIN*1000.); // avoid division by zero

//               numPoints++;
//               real cellAspectRatio = nDist/tDist;
//               if( cellAspectRatio > tallCellRatioBound )
//                 numTall++;
//               tallCellRatio = max( tallCellRatio, cellAspectRatio );
//             }
//           }
        
//         }
//         else if( numberOfDimensions==3 )
//         {
//           int dirp2 = (dir+2) % numberOfDimensions;
//           int irv[3], &ir1 = irv[0], &ir2=irv[1], &ir3=irv[2]; // tangential index shift
//           irv[0]=0; irv[1]=0; irv[2]=0;
//           irv[dirp2]= 1;

//           FOR_3D(i1,i2,i3,I1,I2,I3) 
//           {
//             if( mask(i1,i2,i3)>0 )
//             {
//               // distance "normal" to the boundary for this point
//               real nDist = sqrt( SQR( xLocal(i1+is1,i2+is2,i3+is3,0) - xLocal(i1,i2,i3,0) ) +
//                                  SQR( xLocal(i1+is1,i2+is2,i3+is3,1) - xLocal(i1,i2,i3,1) ) +
//                                  SQR( xLocal(i1+is1,i2+is2,i3+is3,2) - xLocal(i1,i2,i3,2) )  );

//               // tangential distances
//               real tDist1 = sqrt( SQR( xLocal(i1+it1,i2+it2,i3+it3,0) - xLocal(i1,i2,i3,0) ) +
//                                   SQR( xLocal(i1+it1,i2+it2,i3+it3,1) - xLocal(i1,i2,i3,1) ) +
//                                   SQR( xLocal(i1+it1,i2+it2,i3+it3,2) - xLocal(i1,i2,i3,2) )  );

//               real tDist2 = sqrt( SQR( xLocal(i1+ir1,i2+ir2,i3+ir3,0) - xLocal(i1,i2,i3,0) ) +
//                                   SQR( xLocal(i1+ir1,i2+ir2,i3+ir3,1) - xLocal(i1,i2,i3,1) ) +
//                                   SQR( xLocal(i1+ir1,i2+ir2,i3+ir3,2) - xLocal(i1,i2,i3,2) )  );
//               real tDist = max( min(tDist1,tDist2), REAL_MIN*1000.);

//               numPoints++;                
//               real cellAspectRatio = nDist/tDist;
//               if( cellAspectRatio > tallCellRatioBound )
//                 numTall++;

//               tallCellRatio = max( tallCellRatio, cellAspectRatio );
//             }
//           }

//         }
//         else
//         {
//           OV_ABORT("numberOfDimensions!?");
//         }
//         numTall       = ParallelUtility::getSum( numTall );
//         numPoints     = ParallelUtility::getSum( numPoints );
//         tallCellRatio = ParallelUtility::getMaxValue( tallCellRatio );

//         numTallTotal   += numTall;
//         numPointsTotal += numPoints;
//         maxTallCellRatio = max( maxTallCellRatio, tallCellRatio );

        
//         real fractionBad = numTall/max(1.,numPoints);
//         printF(" grid=%3d (side,dir)=(%d,%d) worst tallCellRatio=%9.3e, num bad=%5d (%5.1f %% of cells exceed ratio %4.2f) (name=%s)\n",
//                grid,side,dir,tallCellRatio,numTall,fractionBad*100,tallCellRatioBound, (const char*)mg.getName());
//       }
      
//     } // end for iside 
//   } // for end interface

//   real fractionBad = numTallTotal/max(1.,numPointsTotal);
//   printF(" Maximum tall cell ratio=%8.2e (tallCellRatioBound=%g), total bad cells=%6d (%5.1f %% of %d points).\n",
//          maxTallCellRatio,tallCellRatioBound,numTallTotal, 100.*fractionBad,numPointsTotal);

//   if( maxTallCellRatio > tallCellRatioBound )
//   {
//     printF("\n");
//     printF("ERROR: There are some cells at the interface that are too `tall' in the normal direction.\n");
//     printF("This may cause instablilities for some PDE solvers (e.g. The current fourth-order accurate scheme in CgMx).\n");
//     printF("You should add more grid lines in the normal direction to the offending grids.\n");
//     // printF("Otherwise, increase the tallCellRatioBound if you want to run anyway at your own risk!\n\n");
//     // OV_ABORT("checkForTallCells:ERROR");
//   }
//   else
//   {
//     printF("\n >>> No tall cells exceeding the ratio %4.2f were found. <<<< \n\n",tallCellRatioBound);
//   }
  


//   return maxTallCellRatio;
// }



// =====================================================================================================
/// \brief  Make sure masks match at interfaces.
/// \details Routines that apply interface conditions make require that both sides
///          of the interface have discretization points at the same locations. This function
///          adjusts the mask to make this so. This function should be called after
///          needed points have been marked. 
/// \return number of changes made.       
// =====================================================================================================
int 
matchInterfaceMasks( CompositeGrid & cg, std::vector<InterfaceInfo> & interfaceInfo, Ogen & ogen )
{
  // ---- Make a list of the grids that match at interfaces ----
  // matchInterfaces( cg,interfaceInfo );


  int debug=1;
  FILE *debugFile  = stdout;   // do this for now 
  FILE *pDebugFile = stdout;   // do this for now 

  if( true || debug & 4 )
    printF(" ++++++++++  matchInterfaceMasks : numberOfInterfaces to check = %d....\n",interfaceInfo.size());
  

  const int numberOfDimensions = cg.numberOfDimensions();
  const int numberOfComponentGrids = cg.numberOfComponentGrids();

  
  Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
  Index Jv[3], &J1=Jv[0], &J2=Jv[1], &J3=Jv[2];

  //  ----------------------------------------------------------------------------
  //  ------------------------- loop over interfaces -----------------------------
  //  ----------------------------------------------------------------------------
  int numChanges=0; // counts changes 
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

    // ***** Check that the valid points match on the interface ******
    
    const intArray & mask1 = mg1.mask();
    const intArray & mask2 = mg2.mask();
    OV_GET_SERIAL_ARRAY_CONST(int,mask1,mask1Local);
    OV_GET_SERIAL_ARRAY_CONST(int,mask2,mask2Local);
    
    getBoundaryIndex(mg1.gridIndexRange(),side1,dir1,I1,I2,I3);
    getBoundaryIndex(mg2.gridIndexRange(),side2,dir2,J1,J2,J3);

    if( debug & 2 )
    {
      fprintf(debugFile,"\n --matchInterfaceMasks: Interface : %s=(grid1,side,dir1)=(%i,%i,%i) matches "
           "%s=(grid2,side2,dir2)=(%i,%i,%i) share=%i\n",
           (const char*)mg1.getName(),grid1,side1,dir1,
           (const char*)mg2.getName(),grid2,side2,dir2,share1(side1,dir1));
    }

    if( debug & 8 )
    {

      displayMask(mask1,"mask1",debugFile);
      displayMask(mask2,"mask2",debugFile);

      fflush(debugFile);
    }

    

    int includeGhost=1;
    bool ok1 = ParallelUtility::getLocalArrayBounds(mask1,mask1Local,I1,I2,I3,includeGhost);
    bool ok2 = ParallelUtility::getLocalArrayBounds(mask2,mask2Local,J1,J2,J3,includeGhost);

    // This next check does not always work in parallel -- fix me -- could add check when we create 
    // copies of the mask array *fix me*
    if( ok1 && ok2 &&
        I1.getLength()==J1.getLength() &&  
        I2.getLength()==J2.getLength() &&
        I3.getLength()==J3.getLength() 
       )
     {
       // ****** FIX ME ******
       if( (numberOfDimensions==2  && debug & 2 ) || debug & 4 )
       { 

         displayMask(mask1Local(I1,I2,I3),"mask1Local(I1,I2,I3)",pDebugFile);
         displayMask(mask2Local(J1,J2,J3),"mask2Local(J1,J2,J3)",pDebugFile);
       }

       FOR_3IJD(i1,i2,i3,I1,I2,I3,j1,j2,j3,J1,J2,J3) 
       {
         if(  (mask1Local(i1,i2,i3) & MappedGrid::ISdiscretizationPoint) &&
              (mask2Local(j1,j2,j3) & MappedGrid::ISinterpolationPoint) )
         {
           int jv[3]={j1,j2,j3};
           if( ogen.canDiscretize(mg2, jv) )
           {
             if( debug & 1 )
               printF(">>> mask1(%2d,%2d,%2d)=ISdiscretizationPoint but mask2(%2d,%2d,%2d)=ISinterpolationPoint "
                      "Resetting interp. pt to be discretization.\n", i1,i2,i3,j1,j2,j3);
             mask2Local(j1,j2,j3) = mask1Local(i1,i2,i3);
             numChanges++;
           }
           else
           {
             if( debug & 1 )
               printF(">>> mask1(%2d,%2d,%2d)=ISdiscretizationPoint but mask2(%2d,%2d,%2d)=ISinterpolationPoint "
                      "Resetting discretization pt to interp. pt.\n", i1,i2,i3,j1,j2,j3);
              // ** FIX ME -- CHECK IF THIS IS VALID
              mask1Local(i1,i2,i3)=mask2Local(j1,j2,j3);
             numChanges++;
           }

         }
         else if(  (mask1Local(i1,i2,i3) & MappedGrid::ISinterpolationPoint) &&
                   (mask2Local(j1,j2,j3) & MappedGrid::ISdiscretizationPoint) )
         {
           int iv[3]={i1,i2,i3};
           if( ogen.canDiscretize(mg1,iv) )
           {
             if( debug & 1 )
               printF(">>> mask2(%2d,%2d,%2d)=ISdiscretizationPoint but mask1(%2d,%2d,%2d)=ISinterpolationPoint"
                      " Resetting interp. pt to be discretization.  **CHECK ME**\n",j1,j2,j3,i1,i2,i3);
             // ** FIX ME -- CHECK IF THIS IS VALID
             mask1Local(i1,i2,i3) = mask2Local(j1,j2,j3);
             numChanges++;
           }
           else
           {
             if( debug & 1 )
               printF(">>> mask2(%2d,%2d,%2d)=ISdiscretizationPoint but mask1(%2d,%2d,%2d)=ISinterpolationPoint"
                      " Resetting discr. pt to be an interp. pt. **CHECK ME** \n",j1,j2,j3,i1,i2,i3);
             mask2Local(j1,j2,j3) = mask1Local(i1,i2,i3);
             numChanges++;
           }
         }
         else if(  !(mask1Local(i1,i2,i3) & ISneededPoint) &&
                    (mask2Local(j1,j2,j3) & ISneededPoint) )
         {
           if( mask1Local(i1,i2,i3)==0 )
           {
             if( debug & 1 )
               printF("ERROR: >>> mask2(%2d,%2d,%2d)=ISneededPoint but mask1(%2d,%2d,%2d)=%d"
                     " This should not happen! Resetting maks1=mask2.\n",j1,j2,j3,i1,i2,i3,mask1(i1,i2,i3));
             mask1Local(i1,i2,i3) = mask2Local(j1,j2,j3); 
             numChanges++;
           }
           else
           {
             if( debug & 1 )
               printF(">>> mask2(%2d,%2d,%2d)=ISneededPoint but mask1(%2d,%2d,%2d)=%d"
                    " Setting mask1 to be ISneededPoint.\n",j1,j2,j3,i1,i2,i3,mask1(i1,i2,i3));
             mask1Local(i1,i2,i3) |= ISneededPoint;
             numChanges++;
           }
         }
         else if(  (mask1Local(i1,i2,i3) & ISneededPoint) &&
                  !(mask2Local(j1,j2,j3) & ISneededPoint) )
         {
           if( mask2Local(j1,j2,j3)==0 )
           {
             if( debug & 1 )
               printF("ERROR >>> mask1(%2d,%2d,%2d)=ISneededPoint but mask2(%2d,%2d,%2d)=%d"
                     " THIS SHOULD NOT HAPPEN. Resetting mask2=mask1\n",i1,i2,i3,j1,j2,j3,mask2(j1,j2,j3));
             mask2Local(j1,j2,j3)=mask1Local(i1,i2,i3);
             numChanges++;
           }
           else
           {
             if( debug & 1 )
               printF(">>> mask1(%2d,%2d,%2d)=ISneededPoint but mask2(%2d,%2d,%2d)=%d"
                    " Setting mask2 to be ISneededPoint.\n",i1,i2,i3,j1,j2,j3,mask2(j1,j2,j3));
             mask2Local(j1,j2,j3) |= ISneededPoint;
             numChanges++;
           }
         }

       }


      // // ---- The check for parallel is now done in (MAXWELL) getLocalInterfaceArrays macro ---
      
      // // int maskDiff = max(abs(mask1Local(I1,I2,I3)-mask2Local(J1,J2,J3)));

      // // Check that both masks are discretization points (ignore isNeeded etc.) *wdh* Nov 29, 2020
      // int maskDiff = max( (mask1Local(I1,I2,I3) & MappedGrid::ISdiscretizationPoint) -
      //                     (mask2Local(J1,J2,J3) & MappedGrid::ISdiscretizationPoint) );

      // if( maskDiff > 0 )
      // {
      //   printF("       ERROR: the mask arrays do not match on the interface. This is currently required.\n");
      //   printF("       Try re-generating the grid with more lines in the normal direction, this sometimes fixes this problem.\n");
      //   printF("       See debug file for more info if debug &2 \n");
      //   if( debug & 2 )
      //   {
      //     displayMask(mask1Local(I1,I2,I3),"mask1Local(I1,I2,I3)",pDebugFile);
      //     displayMask(mask2Local(J1,J2,J3),"mask2Local(J1,J2,J3)",pDebugFile);
      //     IntegerArray md(I1,I2,I3);
      //     md = ( (mask1Local(I1,I2,I3) & MappedGrid::ISdiscretizationPoint) -
      //            (mask2Local(J1,J2,J3) & MappedGrid::ISdiscretizationPoint) );
      //     where( md>0 )
      //     {
      //       md=1;
      //     }
          
      //    ::display(md,"difference)",pDebugFile,"%2i ");

      //     fprintf(pDebugFile,"       INFO: interface grid points match: mask arrays agree on the interface.\n\n");
      //   }
      //   // OV_ABORT("ERROR");
      // }
      // else
      // {
      //    printF("       INFO: interface grid points match: mask arrays agree on the interface.\n\n");
      //    if( debug & 2 )
      //    {
      //      displayMask(mask1Local(I1,I2,I3),"mask1Local(I1,I2,I3)",pDebugFile);
      //      displayMask(mask2Local(J1,J2,J3),"mask2Local(J1,J2,J3)",pDebugFile);

      //      fprintf(pDebugFile,"       INFO: interface grid points match: mask arrays agree on the interface.\n\n");
      //    }
      //  }
     } // end if ok 

  }
  // timing(timeForInterfaceBC)+=getCPU()-time0;



  return numChanges;
}
