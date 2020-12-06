
// ===========================================================================================
//  Macro:  Parallel copy near the interface
//
//     Copy data from array u1 into u1nCopy (which is disptributed like u2n) 
//     v1Copy(Jv2)=v1(Jv1) 
// ===========================================================================================
#beginMacro copyLocalInterfaceArrayMacro(v1Copy,v1CopyLocal,v1,v2,Jv1,Jv2)
  v1Copy.partition(v2.getPartition());
  v1Copy.redim(v2.dimension(0),v2.dimension(1),v2.dimension(2),v1.dimension(3));
  getLocalArrayWithGhostBoundaries(v1Copy,v1CopyLocal);
    
  v1CopyLocal=0.; 
  // -- components to copy : 
  Jv1[3]=v1.dimension(3);
  Jv2[3]=Jv1[3];

  ParallelUtility::copy(v1Copy,Jv2,v1,Jv1,nd);  // v1Copy(Jv2)=v1(Jv1)
  v1Copy.updateGhostBoundaries(); // *********** these are currently needed ********************   **FIX ME**
#endMacro 

// ======================================================================================
// Macro: Check that the mask values agree across the interface 
// ======================================================================================
#beginMacro checkMaskOnInterfaceMacro(mask1,I1,I2,I3,mask2,J1,J2,J3)
{
  int includeGhost=1;

  OV_GET_SERIAL_ARRAY(int,mask1,mask1Local);
  OV_GET_SERIAL_ARRAY(int,mask2,mask2Local);
  
  bool ok1 = ParallelUtility::getLocalArrayBounds(mask1,mask1Local,I1,I2,I3,includeGhost);
  bool ok2 = ParallelUtility::getLocalArrayBounds(mask2,mask2Local,J1,J2,J3,includeGhost);

  if( ok1 && ok2 )
  {
    // int maskDiff = max(abs(mask1Local(I1,I2,I3)-mask2Local(J1,J2,J3)));
    // Check that both masks are discretization points (ignore isNeeded etc.) *wdh* Nov 29, 2020
    int maskDiff = max( (mask1Local(I1,I2,I3) & MappedGrid::ISdiscretizationPoint) -
                        (mask2Local(J1,J2,J3) & MappedGrid::ISdiscretizationPoint) );

    if( maskDiff > 0 )
    {
      printf("MX:assignInterfaceBC:ERROR: myid=%d: The mask arrays mask1 and mask2 do not match on the interface. \n"
             "  This is currently required.\n"
             "  Try re-generating the grid with more lines in the normal direction, this sometimes fixes this problem.\n",myid);
      printf("debug=%d\n",debug);
      
      if( debug & 2 )
      {
        fprintf(pDebugFile,"MX:assignInterfaceBC:ERROR: The mask arrays mask1 and mask2 do not match on the interface. \n"
                " maskDiff=%d\n",maskDiff);
        
        fprintf(pDebugFile," [grid1,side1,dir1]=[%d,%d,%d][grid2,side2,dir2]=[%d,%d,%d] \n",grid1,side1,dir1,grid2,side2,dir2);
        
        ::display(mask1Local(I1,I2,I3),"mask1Local on the boundary",pDebugFile);
        ::display(mask2Local(J1,J2,J3),"mask2Local on the boundary",pDebugFile);
        intSerialArray diff(I1,I2,I3);
        diff = abs(mask1Local(I1,I2,I3)-mask2Local(J1,J2,J3));
        ::display(diff,"difference |mask1Local-mask2Local| Local on the boundary",pDebugFile);

        ::displayMask(mask1Local(I1,I2,I3),"mask1Local on the boundary (displayMask)",pDebugFile);
        ::displayMask(mask2Local(J1,J2,J3),"mask2Local on the boundary (displayMask)",pDebugFile);

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
#endMacro
    


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
#beginMacro getLocalInterfaceArraysMacro()

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
    checkMaskOnInterfaceMacro(mask1,Ib1,Ib2,Ib3, mask2b,Jb1,Jb2,Jb3 );
      
      
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
    checkMaskOnInterfaceMacro(mask1b,Ib1,Ib2,Ib3, mask2,Jb1,Jb2,Jb3 );

    

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
  copyLocalInterfaceArrayMacro(u2Copy,u2CopyLocal,u2,u1,Iv2,Iv1);



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
    copyLocalInterfaceArrayMacro(u2nCopy,u2nCopyLocal,u2n,u1n,Iv2,Iv1);
  }

  if( isDispersive && numberOfPolarizationVectors2>0 )
  {
    // --- dispersive material ---
    //   The interface conditions rely on the solutions at t+dt and t
    // Copy:
    //   u1n, u2n : solution at current time (we are assign interface values at the "next" time
    //   p1,p1n, p2,p2n

    copyLocalInterfaceArrayMacro(p2Copy,p2CopyLocal,p2,p1,Iv2,Iv1);
    copyLocalInterfaceArrayMacro(p2nCopy,p2nCopyLocal,p2n,p1n,Iv2,Iv1);
      
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
  copyLocalInterfaceArrayMacro(u1Copy,u1CopyLocal,u1,u2,Iv1,Iv2);


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
    copyLocalInterfaceArrayMacro(u1nCopy,u1nCopyLocal,u1n,u2n,Iv1,Iv2);
  }
  if( isDispersive && numberOfPolarizationVectors1>0 )
  {
    // --- dispersive material ---
    //   The interface conditions rely on the solutions at t+dt and t
    // Copy:
    //   u1n, u2n : solution at current time (we are assign interface values at the "next" time
    //   p1,p1n, p2,p2n

    copyLocalInterfaceArrayMacro(p1Copy,p1CopyLocal,p1,p2,Iv1,Iv2);
    copyLocalInterfaceArrayMacro(p1nCopy,p1nCopyLocal,p1n,p2n,Iv1,Iv2);
      
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

#endMacro

