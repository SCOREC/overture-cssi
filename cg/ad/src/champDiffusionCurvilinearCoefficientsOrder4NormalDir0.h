// Coefficients for CHAMP curvilinear order 4, NORMAL DIRECTION=0.
// File written by champGenCode.maple
const Real KLR = theta;
const Real DLR = beta;
printF("WDH: KLR=%9.3e, DLR=%9.3e\n", KLR,DLR);
#defineMacro diffr2(u,i1,i2,i3,dr)  ((u(i1+1,i2,i3)-u(i1-1,i2,i3))/(2.*dr[0]))
#defineMacro diffs2(u,i1,i2,i3,dr)  ((u(i1,i2+1,i3)-u(i1,i2-1,i3))/(2.*dr[1]))
#defineMacro diffrr2(u,i1,i2,i3,dr) ((u(i1+1,i2,i3)-2.*u(i1,i2,i3)+u(i1-1,i2,i3))/(SQR(dr[0])))
#defineMacro diffss2(u,i1,i2,i3,dr) ((u(i1,i2+1,i3)-2.*u(i1,i2,i3)+u(i1,i2-1,i3))/(SQR(dr[1])))
#defineMacro diffr4(u,i1,i2,i3,dr)  (((1./12.)*u(i1-2,i2,i3)-(2./3.)*u(i1-1,i2,i3)+(2./3.)*u(i1+1,i2,i3)-(1./12.)*u(i1+2,i2,i3))/(dr[0]))
#defineMacro diffs4(u,i1,i2,i3,dr)  (((1./12.)*u(i1,i2-2,i3)-(2./3.)*u(i1,i2-1,i3)+(2./3.)*u(i1,i2+1,i3)-(1./12.)*u(i1,i2+2,i3))/(dr[1]))
#defineMacro diffrr4(u,i1,i2,i3,dr)  ((-1.*u(i1-2,i2,i3)+16.*u(i1-1,i2,i3)-30.*u(i1,i2,i3)+16.*u(i1+1,i2,i3)-1.*u(i1+2,i2,i3))/(12.*SQR(dr[0])))
#defineMacro diffss4(u,i1,i2,i3,dr)  ((-1.*u(i1,i2-2,i3)+16.*u(i1,i2-1,i3)-30.*u(i1,i2,i3)+16.*u(i1,i2+1,i3)-1.*u(i1,i2+2,i3))/(12.*SQR(dr[1])))
#defineMacro diffrrr2(u,i1,i2,i3,dr)  (((-1./2.)*u(i1-2,i2,i3)+1.*u(i1-1,i2,i3)-1.*u(i1+1,i2,i3)+(1./2.)*u(i1+2,i2,i3))/(pow(dr[0],3)))
#defineMacro diffsss2(u,i1,i2,i3,dr)  (((-1./2.)*u(i1,i2-2,i3)+1.*u(i1,i2-1,i3)-1.*u(i1,i2+1,i3)+(1./2.)*u(i1,i2+2,i3))/(pow(dr[1],3)))
#defineMacro diffrrrr2(u,i1,i2,i3,dr)  ((u(i1-2,i2,i3)-4.*u(i1-1,i2,i3)+6.*u(i1,i2,i3)-4.*u(i1+1,i2,i3)+u(i1+2,i2,i3))/(pow(dr[0],4)))
#defineMacro diffssss2(u,i1,i2,i3,dr)  ((u(i1,i2-2,i3)-4.*u(i1,i2-1,i3)+6.*u(i1,i2,i3)-4.*u(i1,i2+1,i3)+u(i1,i2+2,i3))/(pow(dr[1],4)))
#defineMacro diffrs4(u,i1,i2,i3,dr) ((64.*u(i1+1,i2+1,i3)-64.*u(i1+1,i2-1,i3)-8.*u(i1+1,i2+2,i3)+8.*u(i1+1,i2-2,i3)-64.*u(i1-1,i2+1,i3)+64.*u(i1-1,i2-1,i3)+8.*u(i1-1,i2+2,i3)-8.*u(i1-1,i2-2,i3)-8.*u(i1+2,i2+1,i3)+8.*u(i1+2,i2-1,i3)+u(i1+2,i2+2,i3)-u(i1+2,i2-2,i3)+8.*u(i1-2,i2+1,i3)-8.* u(i1-2,i2-1,i3)-u(i1-2,i2+2,i3)+u(i1-2,i2-2,i3))/(dr[0]*dr[1]*144.))
#defineMacro diffrs2(u,i1,i2,i3,dr) ((diffs2(u,i1+1,i2,i3,dr)-diffs2(u,i1-1,i2,i3,dr))/(2.*dr[0]))
// Derivatives of the metrics: order 4
#defineMacro diffMr4(rx,i1,i2,i3,m1,m2,dr)  (((1./12.)*rx(i1-2,i2,i3,m1,m2)-(2./3.)*rx(i1-1,i2,i3,m1,m2)+(2./3.)*rx(i1+1,i2,i3,m1,m2)-(1./12.)*rx(i1+2,i2,i3,m1,m2))/(dr[0]))
#defineMacro diffMs4(rx,i1,i2,i3,m1,m2,dr)  (((1./12.)*rx(i1,i2-2,i3,m1,m2)-(2./3.)*rx(i1,i2-1,i3,m1,m2)+(2./3.)*rx(i1,i2+1,i3,m1,m2)-(1./12.)*rx(i1,i2+2,i3,m1,m2))/(dr[1]))
#defineMacro diffMx4(rx,i1,i2,i3,m1,m2,dr)  rx(i1,i2,i3,0,0)*diffMr4(rx,i1,i2,i3,m1,m2,dr) + rx(i1,i2,i3,1,0)*diffMs4(rx,i1,i2,i3,m1,m2,dr)
#defineMacro diffMy4(rx,i1,i2,i3,m1,m2,dr)  rx(i1,i2,i3,0,1)*diffMr4(rx,i1,i2,i3,m1,m2,dr) + rx(i1,i2,i3,1,1)*diffMs4(rx,i1,i2,i3,m1,m2,dr)
// Derivatives of the metrics: order 2
#defineMacro diffMr2(rx,i1,i2,i3,m1,m2,dr)  (-rx(i1-1,i2,i3,m1,m2)+rx(i1+1,i2,i3,m1,m2))/(2.*dr[0])
#defineMacro diffMs2(rx,i1,i2,i3,m1,m2,dr)  (-rx(i1,i2-1,i3,m1,m2)+rx(i1,i2+1,i3,m1,m2))/(2.*dr[1])
#defineMacro diffMx2(rx,i1,i2,i3,m1,m2,dr)  rx(i1,i2,i3,0,0)*diffMr2(rx,i1,i2,i3,m1,m2,dr) + rx(i1,i2,i3,1,0)*diffMs2(rx,i1,i2,i3,m1,m2,dr)
#defineMacro diffMy2(rx,i1,i2,i3,m1,m2,dr)  rx(i1,i2,i3,0,1)*diffMr2(rx,i1,i2,i3,m1,m2,dr) + rx(i1,i2,i3,1,1)*diffMs2(rx,i1,i2,i3,m1,m2,dr)
// Fourth-order one-sided approximations
// Here are fully one-sided: *check me* *wdh* April 15, 2022
#defineMacro diffMr4OneSidedPlus(rx,i1,i2,i3,m1,m2,dr)  (( (-25./12.)*rx(i1,i2,i3,m1,m2)+(4.)*rx(i1+1,i2,i3,m1,m2)-(3.)*rx(i1+2,i2,i3,m1,m2)+(4./3.)*rx(i1+3,i2,i3,m1,m2)-(1./4.)*rx(i1+4,i2,i3,m1,m2))/(dr[0])) 
#defineMacro diffMr4OneSidedMinus(rx,i1,i2,i3,m1,m2,dr) (( (+25./12.)*rx(i1,i2,i3,m1,m2)-(4.)*rx(i1-1,i2,i3,m1,m2)+(3.)*rx(i1-2,i2,i3,m1,m2)-(4./3.)*rx(i1-3,i2,i3,m1,m2)+(1./4.)*rx(i1-4,i2,i3,m1,m2))/(dr[0])) 
#defineMacro diffMs4OneSidedPlus(rx,i1,i2,i3,m1,m2,dr)  (( (-25./12.)*rx(i1,i2,i3,m1,m2)+(4.)*rx(i1,i2+1,i3,m1,m2)-(3.)*rx(i1,i2+2,i3,m1,m2)+(4./3.)*rx(i1,i2+3,i3,m1,m2)-(1./4.)*rx(i1,i2+4,i3,m1,m2))/(dr[1])) 
#defineMacro diffMs4OneSidedMinus(rx,i1,i2,i3,m1,m2,dr) (( (+25./12.)*rx(i1,i2,i3,m1,m2)-(4.)*rx(i1,i2-1,i3,m1,m2)+(3.)*rx(i1,i2-2,i3,m1,m2)-(4./3.)*rx(i1,i2-3,i3,m1,m2)+(1./4.)*rx(i1,i2-4,i3,m1,m2))/(dr[1])) 
#defineMacro diffa2(u,i1,i2,i3,dr) diffs2(u,i1,i2,i3,dr)
#defineMacro diffaa2(u,i1,i2,i3,dr) diffss2(u,i1,i2,i3,dr)
#defineMacro diffa(u,i1,i2,i3,dr) diffs4(u,i1,i2,i3,dr)
#defineMacro diffaa(u,i1,i2,i3,dr) diffss4(u,i1,i2,i3,dr)
#defineMacro diffaaa(u,i1,i2,i3,dr) diffsss2(u,i1,i2,i3,dr)
#defineMacro diffaaaa(u,i1,i2,i3,dr) diffssss2(u,i1,i2,i3,dr)
#defineMacro extrap3(uu,k1,k2,k3,ks1,ks2,ks3) \
            (  3.*uu(k1,k2,k3)\
              -3.*uu(k1+  ks1,k2+  ks2,k3+  ks3)\
              +1.*uu(k1+2*ks1,k2+2*ks2,k3+2*ks3))
#defineMacro extrap5(uu,k1,k2,k3,ks1,ks2,ks3) \
            (   5.*uu(k1,k2,k3)\
              -10.*uu(k1+  ks1,k2+  ks2,k3+  ks3)\
              +10.*uu(k1+2*ks1,k2+2*ks2,k3+2*ks3)\
               -5.*uu(k1+3*ks1,k2+3*ks2,k3+3*ks3)\
                  +uu(k1+4*ks1,k2+4*ks2,k3+4*ks3))
// Macro : fill end ghost points on b1, b2, by periodicity or extrapolation
#beginMacro fillEndGhost2(side1,axis1,mg,gid,b1,b2,startGhost,endGhost)
  axist = (axis1 + 1) % numberOfDimensions; // tangential direction
  for( int sidet=0; sidet<=1; sidet++ ) // top and bottom of face
  { 
    ks1=ks2=ks3=0;  ksv[axist]=1-2*sidet;
    I1=Range(gid(0,0),gid(0,0)); 
    I2=Range(gid(0,1),gid(0,1)); 
    I3=Range(gid(0,2),gid(0,2)); 
    Iv[axis1]=Range(gid(side1,axis1),gid(side1,axis1));
    Iv[axist]=gid(sidet,axist);
    for( int ghost=startGhost; ghost<=endGhost; ghost++ )
    { 
      Ig1=I1-ks1*ghost; Ig2=I2-ks2*ghost; Ig3=I3-ks3*ghost;
      if( mg.boundaryCondition(sidet,axist)<0 )
      { // periodic
        Index Ip1 = Ig1 + (gid(1,0)-gid(0,0))*ks1;
        Index Ip2 = Ig2 + (gid(1,1)-gid(0,1))*ks2;
        Index Ip3 = Ig3 + (gid(1,2)-gid(0,2))*ks3;
        b1(Ig1,Ig2,Ig3) = b1(Ip1,Ip2,Ip3); 
        b2(Ig1,Ig2,Ig3) = b2(Ip1,Ip2,Ip3); 
      }
      else;
      {
        b1(Ig1,Ig2,Ig3) = extrap5(b1,Ig1+ks1,Ig2+ks2,Ig3+ks3,ks1,ks2,ks3); 
        b2(Ig1,Ig2,Ig3) = extrap5(b2,Ig1+ks1,Ig2+ks2,Ig3+ks3,ks1,ks2,ks3); 
      }                  
    } // end for ghost 
  } // end for side t
#endMacro
         
// Macro : fill end ghost points on crr, crs, ... by periodicity or extrapolation
#beginMacro fillEndGhost(side1,axis1,mg,gid,crr,crs,css,cr,cs)
  axist = (axis1 + 1) % numberOfDimensions; // tangential direction
  for( int sidet=0; sidet<=1; sidet++ ) // top and bottom of face
  { 
    ks1=ks2=ks3=0;  ksv[axist]=1-2*sidet;
    I1=Range(gid(0,0),gid(0,0)); 
    I2=Range(gid(0,1),gid(0,1)); 
    I3=Range(gid(0,2),gid(0,2)); 
    Iv[axis1]=Range(gid(side1,axis1)-1,gid(side1,axis1)+1);
    Iv[axist]=gid(sidet,axist);
    for( int ghost=1; ghost<=numGhost+1; ghost++ )
    { 
      Ig1=I1-ks1*ghost; Ig2=I2-ks2*ghost; Ig3=I3-ks3*ghost;
      if( mg.boundaryCondition(sidet,axist)<0 )
      { // periodic
        Index Ip1 = Ig1 + (gid(1,0)-gid(0,0))*ks1;
        Index Ip2 = Ig2 + (gid(1,1)-gid(0,1))*ks2;
        Index Ip3 = Ig3 + (gid(1,2)-gid(0,2))*ks3;
        crr(Ig1,Ig2,Ig3) = crr(Ip1,Ip2,Ip3); 
        crs(Ig1,Ig2,Ig3) = crs(Ip1,Ip2,Ip3); 
        css(Ig1,Ig2,Ig3) = css(Ip1,Ip2,Ip3); 
        cr(Ig1,Ig2,Ig3)  =  cr(Ip1,Ip2,Ip3); 
        cs(Ig1,Ig2,Ig3)  =  cs(Ip1,Ip2,Ip3); 
      }
      else;
      {
        if( ghost==3 )
        { 
          crr(Ig1,Ig2,Ig3) = extrap5(crr,Ig1+ks1,Ig2+ks2,Ig3+ks3,ks1,ks2,ks3); 
          crs(Ig1,Ig2,Ig3) = extrap5(crs,Ig1+ks1,Ig2+ks2,Ig3+ks3,ks1,ks2,ks3); 
          css(Ig1,Ig2,Ig3) = extrap5(css,Ig1+ks1,Ig2+ks2,Ig3+ks3,ks1,ks2,ks3); 
        } 
        cr(Ig1,Ig2,Ig3) = extrap5(cr,Ig1+ks1,Ig2+ks2,Ig3+ks3,ks1,ks2,ks3); 
        cs(Ig1,Ig2,Ig3) = extrap5(cs,Ig1+ks1,Ig2+ks2,Ig3+ks3,ks1,ks2,ks3); 
      }                  
    } // end for ghost 
  } // end for side t
#endMacro
         
const int axis1L=axis,  axis2L=(axis1L+1) % numberOfDimensions;
const int axis1R=axis2, axis2R=(axis1R+1) % numberOfDimensions;
const IntegerArray & gir1 = mg.indexRange();
const IntegerArray & gir2 = mg2.indexRange();
int ksv[3], &ks1=ksv[0], &ks2=ksv[1], &ks3=ksv[2];
int axist;
//  ----- We need b1,b2 at 4 extra tangential ghosts (order=4)---- 
//  since L1 is evaluated at 4 extra ghost 
int extra=2*numGhost;
getBoundaryIndex(mg.indexRange(), side, axis, Ib1,Ib2,Ib3,extra);
getBoundaryIndex(mg2.indexRange(),side2,axis2,Jb1,Jb2,Jb3,extra);
RealArray n1L(Ib1,Ib2,Ib3), n2L(Ib1,Ib2,Ib3)  ; // save the normal, it needs to be extended too;
RealArray b1L(Ib1,Ib2,Ib3), b2L(Ib1,Ib2,Ib3);
RealArray b1R(Jb1,Jb2,Jb3), b2R(Jb1,Jb2,Jb3);
// Evaluate the normal derivative coefficients at points on the boundary. Use right-normal = -left-normal.
extra=numGhost;  // we can only directly evaluate 2 ghost using metrics 
getBoundaryIndex(mg.indexRange(), side, axis, Ib1,Ib2,Ib3,extra);
getBoundaryIndex(mg2.indexRange(),side2,axis2,Jb1,Jb2,Jb3,extra);
FOR_3IJD(i1,i2,i3,Ib1,Ib2,Ib3,j1,j2,j3,Jb1,Jb2,Jb3)
{
  n1L(i1,i2,i3) = normal(i1,i2,i3,0);
  n2L(i1,i2,i3) = normal(i1,i2,i3,1);
  b1L(i1,i2,i3) = -( normal(i1,i2,i3,0)*rxL(i1,i2,i3,0,0) + normal(i1,i2,i3,1)*rxL(i1,i2,i3,0,1)) ;
  b2L(i1,i2,i3) = -( normal(i1,i2,i3,0)*rxL(i1,i2,i3,1,0) + normal(i1,i2,i3,1)*rxL(i1,i2,i3,1,1)) ;
  b1R(j1,j2,j3) = -( normal(i1,i2,i3,0)*rxR(j1,j2,j3,0,0) + normal(i1,i2,i3,1)*rxR(j1,j2,j3,0,1)) ;
  b2R(j1,j2,j3) = -( normal(i1,i2,i3,0)*rxR(j1,j2,j3,1,0) + normal(i1,i2,i3,1)*rxR(j1,j2,j3,1,1)) ;
} // end FOR_3IJD

// Fill 2nd ghost in tangential directions for b1, b2, and nL 
const int ghostToFillStart = numGhost+1;
const int ghostToFillEnd   = 2*numGhost;
fillEndGhost2(side ,axis ,mg, gir1,n1L,n2L,ghostToFillStart,ghostToFillEnd);
fillEndGhost2(side ,axis ,mg, gir1,b1L,b2L,ghostToFillStart,ghostToFillEnd);
fillEndGhost2(side2,axis2,mg2,gir2,b1R,b2R,ghostToFillStart,ghostToFillEnd);

// --- Set the optimized Schwarz parameter: Sl = pl/dx (now that we know b1,b2)  ---- 
//  For now set a constant value on the face 
if( Sl<0 )
{
  j1=Jb1.getBase(); j2=Jb2.getBase(); j3=Jb3.getBase(); // take b1R from this point 
  const Real dnR = fabs( drR[axis2]/b1R(j1,j2,j3) );    // approx. normal grid spacing on right
  Sl = pl/dnR; 
  champParameters(4,side,axis,grid)=Sl; // save Sl
}
// *NOTE* WE NEED AN EXTRA GHOST IN TANGENTIAL DIRECTION FOR CURVILINEAR 
// Index bounds for points near the boundary: add 3 extra points in tangential directions and 1 in normal direction.
extra=numGhost+1;
getBoundaryIndex(mg.indexRange(),  side,axis, I1,I2,I3,extra);
getBoundaryIndex(mg2.indexRange(),side2,axis2,J1,J2,J3,extra);
int extraNormal=numGhost-1;
Iv[axis1L]=Range( Iv[axis1L].getBase()-extraNormal,Iv[axis1L].getBound()+extraNormal );
Jv[axis1R]=Range( Jv[axis1R].getBase()-extraNormal,Jv[axis1R].getBound()+extraNormal );
RealArray crrL(I1,I2,I3), crsL(I1,I2,I3), cssL(I1,I2,I3), crL(I1,I2,I3), csL(I1,I2,I3);
RealArray crrR(J1,J2,J3), crsR(J1,J2,J3), cssR(J1,J2,J3), crR(J1,J2,J3), csR(J1,J2,J3);
RealArray rxxL(3,3), rxxR(3,3);
RealArray rxyL(3,3), rxyR(3,3);
// We can fill in crr, crs, and ss to 2 ghost since they do not involve terms using rxx, sxx, etc. 
extra=numGhost;
getBoundaryIndex(mg.indexRange(),  side,axis, I1,I2,I3,extra);
getBoundaryIndex(mg2.indexRange(),side2,axis2,J1,J2,J3,extra);
Iv[axis1L]=Range( Iv[axis1L].getBase()-extraNormal,Iv[axis1L].getBound()+extraNormal );
Jv[axis1R]=Range( Jv[axis1R].getBase()-extraNormal,Jv[axis1R].getBound()+extraNormal );
FOR_3IJ(i1,i2,i3,I1,I2,I3,j1,j2,j3,J1,J2,J3)
{
  // Evaluate coefficients of the Laplacian at points near the boundary
  // uxx = (rx*Dr + sx*Ds)*[ (rx)*ur + sx*us ]
  // uxx = (rx)^2 urr + 2*(rx*sx)*urs + (sx)^2 *uss + (rxx)*ur + (sxx)*us
  // uyy = (ry)^2 urr + 2*(ry*sy)*urs + (sy)^2 *uss + (ryy)*ur + (syy)*us
  crrL(i1,i2,i3) = SQR(rxL(i1,i2,i3,0,0)) + SQR(rxL(i1,i2,i3,0,1));                                // (rx)^2 + (ry)^2 
  crsL(i1,i2,i3) = 2.*(rxL(i1,i2,i3,0,0)*rxL(i1,i2,i3,1,0) + rxL(i1,i2,i3,0,1)*rxL(i1,i2,i3,1,1)); // 2( rx*sx + ry*sy ) 
  cssL(i1,i2,i3) = SQR(rxL(i1,i2,i3,1,0)) + SQR(rxL(i1,i2,i3,1,1));                                // (sx)^2 + (sy)^2
  // These are done below
  //crL(i1,i2,i3)  = rxxL(0,0) + rxyL(0,1);  // rxx + ryy
  //csL(i1,i2,i3)  = rxxL(1,0) + rxyL(1,1);  // sxx + syy
  // Evaluate coefficients of the Laplacian at points near the boundary
  crrR(j1,j2,j3) = SQR(rxR(j1,j2,j3,0,0)) + SQR(rxR(j1,j2,j3,0,1));
  crsR(j1,j2,j3) = 2.*( rxR(j1,j2,j3,0,0)*rxR(j1,j2,j3,1,0) + rxR(j1,j2,j3,0,1)*rxR(j1,j2,j3,1,1) );
  cssR(j1,j2,j3) = SQR(rxR(j1,j2,j3,1,0)) + SQR(rxR(j1,j2,j3,1,1));
  // These are done below
  //crR(j1,j2,j3)  = rxxR(0,0) + rxyR(0,1);  // rxx + ryy
  //csR(j1,j2,j3)  = rxxR(1,0) + rxyR(1,1);  // sxx + syy
  // printF("WDH: crrL=%9.3e, crsL=%9.3e, cssL=%9.3e\n",crrL(i1,i2,i3),crsL(i1,i2,i3),cssL(i1,i2,i3) );
  // printF("WDH: crrR=%9.3e, crsR=%9.3e, cssR=%9.3e\n",crrR(j1,j2,j3),crsR(j1,j2,j3),cssR(j1,j2,j3) );
} // end FOR_3IJ

// We can only compute cr and cs up to the boundary since we need two ghost points to evaluate rxx, rxy, sxx, sxy
getBoundaryIndex(mg.indexRange(),side,axis,I1,I2,I3);
getBoundaryIndex(mg2.indexRange(),side2,axis2,J1,J2,J3);
Iv[axis1L]= Range( Iv[axis1L].getBase()-extraNormal,Iv[axis1L].getBound()+extraNormal );
Jv[axis1R]= Range( Jv[axis1R].getBase()-extraNormal,Jv[axis1R].getBound()+extraNormal );
Real rxr,rxs;
FOR_3IJ(i1,i2,i3,I1,I2,I3,j1,j2,j3,J1,J2,J3)
{
  for( int m2=0; m2<numberOfDimensions; m2++ )
  {
    for( int m1=0; m1<numberOfDimensions; m1++ )
    {
        // Evaluate derivatives of metrics to order 4 .. we are computing more than we need here *fix me*
        if( i1-2 >= mg.dimension(0,0) && i1+2 <= mg.dimension(1,0) )
          rxr = diffMr4(rxL,i1,i2,i3,m1,m2,drL);              // centered 
        else if( i1-2 < mg.dimension(0,0) )
          rxr = diffMr4OneSidedPlus(rxL,i1,i2,i3,m1,m2,drL);   // one-sided to right
        else 
          rxr = diffMr4OneSidedMinus(rxL,i1,i2,i3,m1,m2,drL);  // one-sided to left
        if( i2-2 >= mg.dimension(0,1) && i2+2 <= mg.dimension(1,1) )
          rxs = diffMs4(rxL,i1,i2,i3,m1,m2,drL);              // centered 
        else if( i2-2 < mg.dimension(0,1) )
          rxs = diffMs4OneSidedPlus(rxL,i1,i2,i3,m1,m2,drL);   // one-sided to right
        else 
          rxs = diffMs4OneSidedMinus(rxL,i1,i2,i3,m1,m2,drL);  // one-sided to left
        rxxL(m1,m2) = rxL(i1,i2,i3,0,0)*rxr + rxL(i1,i2,i3,1,0)*rxs;
        rxyL(m1,m2) = rxL(i1,i2,i3,0,1)*rxr + rxL(i1,i2,i3,1,1)*rxs;
                                                      
        if( j1-2 >= mg2.dimension(0,0) && j1+2 <= mg2.dimension(1,0) )
          rxr = diffMr4(rxR,j1,j2,j3,m1,m2,drR);  // centered 
        else if( j1-2 < mg2.dimension(0,0) )
          rxr = diffMr4OneSidedPlus(rxR,j1,j2,j3,m1,m2,drR);  // one-sided to right
        else 
          rxr = diffMr4OneSidedMinus(rxR,j1,j2,j3,m1,m2,drR);  // one-sided to left
        if( j2-2 >= mg2.dimension(0,1) && j2+2 <= mg2.dimension(1,1) )
          rxs = diffMs4(rxR,j1,j2,j3,m1,m2,drR);  // centered 
        else if( j2-2 < mg2.dimension(0,1) )
          rxs = diffMs4OneSidedPlus(rxR,j1,j2,j3,m1,m2,drR);  // one-sided to right
        else 
          rxs = diffMs4OneSidedMinus(rxR,j1,j2,j3,m1,m2,drR);  // one-sided to left
        rxxR(m1,m2) = rxR(j1,j2,j3,0,0)*rxr + rxR(j1,j2,j3,1,0)*rxs;
        rxyR(m1,m2) = rxR(j1,j2,j3,0,1)*rxr + rxR(j1,j2,j3,1,1)*rxs;
                                                                    
    }
  }
  // Evaluate coefficients of the Laplacian at points near the boundary
  // uxx = (rx*Dr + sx*Ds)*[ (rx)*ur + sx*us ]
  // uxx = (rx)^2 urr + 2*(rx*sx)*urs + (sx)^2 *uss + (rxx)*ur + (sxx)*us
  // uyy = (ry)^2 urr + 2*(ry*sy)*urs + (sy)^2 *uss + (ryy)*ur + (syy)*us
  crL(i1,i2,i3)  = rxxL(0,0) + rxyL(0,1);  // rxx + ryy
  csL(i1,i2,i3)  = rxxL(1,0) + rxyL(1,1);  // sxx + syy
  // Evaluate coefficients of the Laplacian at points near the boundary
  crR(j1,j2,j3)  = rxxR(0,0) + rxyR(0,1);  // rxx + ryy
  csR(j1,j2,j3)  = rxxR(1,0) + rxyR(1,1);  // sxx + syy
  if( advectionIsOn )  // include advection terms scaled by 1/D : (u/D).grad 
  {
    getAdvectionCoefficients(i1,i2,i3,u1DL,u2DL, j1,j2,j3,u1DR,u2DR );
    crL(i1,i2,i3) = crL(i1,i2,i3) - ( u1DL*rxL(i1,i2,i3,0,0) + u2DL*rxL(i1,i2,i3,0,1) ); 
    csL(i1,i2,i3) = csL(i1,i2,i3) - ( u1DL*rxL(i1,i2,i3,1,0) + u2DL*rxL(i1,i2,i3,1,1) ); 
    crR(j1,j2,j3) = crR(j1,j2,j3) - ( u1DR*rxR(j1,j2,j3,0,0) + u2DR*rxR(j1,j2,j3,0,1) ); 
    csR(j1,j2,j3) = csR(j1,j2,j3) - ( u1DR*rxR(j1,j2,j3,1,0) + u2DR*rxR(j1,j2,j3,1,1) ); 
  }
  // printF("WDH: (i1,i2)=(%4d,%4d) crL=%9.3e csL=%9.3e\n",i1,i2,crL(i1,i2,i3),csL(i1,i2,i3) );
  // printF("WDH: (j1,j2)=(%4d,%4d) crR=%9.3e csR=%9.3e\n",j1,j2,crR(j1,j2,j3),csR(j1,j2,j3) );
} // end FOR_3IJ

// --- Now fill in one extra ghost for crr,crs,css and two extra ghost for cr,cs --- 
// --- Fill ghost points on ends through extrapolation or periodicity --- 
// Fill ghost in tangential directions for crr, crs, etc.
fillEndGhost(side ,axis ,mg, gir1,crrL,crsL,cssL,crL,csL);
fillEndGhost(side2,axis2,mg2,gir2,crrR,crsR,cssR,crR,csR);
         
// ---- Evaluate derivatives of coefficients : (crr).rr, (csr).rs etc. 
//  L3 : needs (crr).r (crr).s etc at one tangential ghost. 
extra=numGhost-1;
getBoundaryIndex(mg.indexRange(),  side,axis, Ib1,Ib2,Ib3,extra);
getBoundaryIndex(mg2.indexRange(),side2,axis2,Jb1,Jb2,Jb3,extra);
RealArray b1Lrr(Ib1,Ib2,Ib3), b1Lrs(Ib1,Ib2,Ib3), b1Lss(Ib1,Ib2,Ib3), b1Lr(Ib1,Ib2,Ib3), b1Ls(Ib1,Ib2,Ib3);
RealArray b2Lrr(Ib1,Ib2,Ib3), b2Lrs(Ib1,Ib2,Ib3), b2Lss(Ib1,Ib2,Ib3), b2Lr(Ib1,Ib2,Ib3), b2Ls(Ib1,Ib2,Ib3);
RealArray crrLrr(Ib1,Ib2,Ib3), crrLrs(Ib1,Ib2,Ib3), crrLss(Ib1,Ib2,Ib3), crrLr(Ib1,Ib2,Ib3), crrLs(Ib1,Ib2,Ib3);
RealArray crsLrr(Ib1,Ib2,Ib3), crsLrs(Ib1,Ib2,Ib3), crsLss(Ib1,Ib2,Ib3), crsLr(Ib1,Ib2,Ib3), crsLs(Ib1,Ib2,Ib3);
RealArray cssLrr(Ib1,Ib2,Ib3), cssLrs(Ib1,Ib2,Ib3), cssLss(Ib1,Ib2,Ib3), cssLr(Ib1,Ib2,Ib3), cssLs(Ib1,Ib2,Ib3);
RealArray crLrr(Ib1,Ib2,Ib3), crLrs(Ib1,Ib2,Ib3), crLss(Ib1,Ib2,Ib3), crLr(Ib1,Ib2,Ib3), crLs(Ib1,Ib2,Ib3);
RealArray csLrr(Ib1,Ib2,Ib3), csLrs(Ib1,Ib2,Ib3), csLss(Ib1,Ib2,Ib3), csLr(Ib1,Ib2,Ib3), csLs(Ib1,Ib2,Ib3);
// Evaluate derivatives of coefficients at points on the boundary (left)
FOR_3IJ(i1,i2,i3,Ib1,Ib2,Ib3,j1,j2,j3,Jb1,Jb2,Jb3)
{
  // tangential directives of b1,b2 only 
  b1Lss(i1,i2,i3) = diffss2(b1L,i1,i2,i3,drL);
   b1Ls(i1,i2,i3) =  diffs2(b1L,i1,i2,i3,drL);
  b2Lss(i1,i2,i3) = diffss2(b2L,i1,i2,i3,drL);
   b2Ls(i1,i2,i3) =  diffs2(b2L,i1,i2,i3,drL);
  // normal and tangential directives of c 
  crrLrr(i1,i2,i3) = diffrr2(crrL,i1,i2,i3,drL);
  crrLss(i1,i2,i3) = diffss2(crrL,i1,i2,i3,drL);
  crrLrs(i1,i2,i3) = diffrs2(crrL,i1,i2,i3,drL);
   crrLr(i1,i2,i3) =  diffr2(crrL,i1,i2,i3,drL);
   crrLs(i1,i2,i3) =  diffs2(crrL,i1,i2,i3,drL);
  crsLrr(i1,i2,i3) = diffrr2(crsL,i1,i2,i3,drL);
  crsLss(i1,i2,i3) = diffss2(crsL,i1,i2,i3,drL);
  crsLrs(i1,i2,i3) = diffrs2(crsL,i1,i2,i3,drL);
   crsLr(i1,i2,i3) =  diffr2(crsL,i1,i2,i3,drL);
   crsLs(i1,i2,i3) =  diffs2(crsL,i1,i2,i3,drL);
  cssLrr(i1,i2,i3) = diffrr2(cssL,i1,i2,i3,drL);
  cssLss(i1,i2,i3) = diffss2(cssL,i1,i2,i3,drL);
  cssLrs(i1,i2,i3) = diffrs2(cssL,i1,i2,i3,drL);
   cssLr(i1,i2,i3) =  diffr2(cssL,i1,i2,i3,drL);
   cssLs(i1,i2,i3) =  diffs2(cssL,i1,i2,i3,drL);
  crLrr(i1,i2,i3) = diffrr2(crL,i1,i2,i3,drL);
  crLss(i1,i2,i3) = diffss2(crL,i1,i2,i3,drL);
  crLrs(i1,i2,i3) = diffrs2(crL,i1,i2,i3,drL);
   crLr(i1,i2,i3) =  diffr2(crL,i1,i2,i3,drL);
   crLs(i1,i2,i3) =  diffs2(crL,i1,i2,i3,drL);
  csLrr(i1,i2,i3) = diffrr2(csL,i1,i2,i3,drL);
  csLss(i1,i2,i3) = diffss2(csL,i1,i2,i3,drL);
  csLrs(i1,i2,i3) = diffrs2(csL,i1,i2,i3,drL);
   csLr(i1,i2,i3) =  diffr2(csL,i1,i2,i3,drL);
   csLs(i1,i2,i3) =  diffs2(csL,i1,i2,i3,drL);
} // end FOR_3D
RealArray b1Rrr(Jb1,Jb2,Jb3), b1Rrs(Jb1,Jb2,Jb3), b1Rss(Jb1,Jb2,Jb3), b1Rr(Jb1,Jb2,Jb3), b1Rs(Jb1,Jb2,Jb3);
RealArray b2Rrr(Jb1,Jb2,Jb3), b2Rrs(Jb1,Jb2,Jb3), b2Rss(Jb1,Jb2,Jb3), b2Rr(Jb1,Jb2,Jb3), b2Rs(Jb1,Jb2,Jb3);
RealArray crrRrr(Jb1,Jb2,Jb3), crrRrs(Jb1,Jb2,Jb3), crrRss(Jb1,Jb2,Jb3), crrRr(Jb1,Jb2,Jb3), crrRs(Jb1,Jb2,Jb3);
RealArray crsRrr(Jb1,Jb2,Jb3), crsRrs(Jb1,Jb2,Jb3), crsRss(Jb1,Jb2,Jb3), crsRr(Jb1,Jb2,Jb3), crsRs(Jb1,Jb2,Jb3);
RealArray cssRrr(Jb1,Jb2,Jb3), cssRrs(Jb1,Jb2,Jb3), cssRss(Jb1,Jb2,Jb3), cssRr(Jb1,Jb2,Jb3), cssRs(Jb1,Jb2,Jb3);
RealArray crRrr(Jb1,Jb2,Jb3), crRrs(Jb1,Jb2,Jb3), crRss(Jb1,Jb2,Jb3), crRr(Jb1,Jb2,Jb3), crRs(Jb1,Jb2,Jb3);
RealArray csRrr(Jb1,Jb2,Jb3), csRrs(Jb1,Jb2,Jb3), csRss(Jb1,Jb2,Jb3), csRr(Jb1,Jb2,Jb3), csRs(Jb1,Jb2,Jb3);
// Evaluate derivatives of coefficients at points on the boundary (right)
FOR_3IJ(i1,i2,i3,Ib1,Ib2,Ib3,j1,j2,j3,Jb1,Jb2,Jb3)
{
 // tangential directives of b1,b2 only 
  b1Rss(j1,j2,j3) = diffss2(b1R,j1,j2,j3,drR);
   b1Rs(j1,j2,j3) =  diffs2(b1R,j1,j2,j3,drR);
  b2Rss(j1,j2,j3) = diffss2(b2R,j1,j2,j3,drR);
   b2Rs(j1,j2,j3) =  diffs2(b2R,j1,j2,j3,drR);
  // normal and tangential directives of c 
  crrRrr(j1,j2,j3) = diffrr2(crrR,j1,j2,j3,drR);
  crrRss(j1,j2,j3) = diffss2(crrR,j1,j2,j3,drR);
  crrRrs(j1,j2,j3) = diffrs2(crrR,j1,j2,j3,drR);
   crrRr(j1,j2,j3) =  diffr2(crrR,j1,j2,j3,drR);
   crrRs(j1,j2,j3) =  diffs2(crrR,j1,j2,j3,drR);
  crsRrr(j1,j2,j3) = diffrr2(crsR,j1,j2,j3,drR);
  crsRss(j1,j2,j3) = diffss2(crsR,j1,j2,j3,drR);
  crsRrs(j1,j2,j3) = diffrs2(crsR,j1,j2,j3,drR);
   crsRr(j1,j2,j3) =  diffr2(crsR,j1,j2,j3,drR);
   crsRs(j1,j2,j3) =  diffs2(crsR,j1,j2,j3,drR);
  cssRrr(j1,j2,j3) = diffrr2(cssR,j1,j2,j3,drR);
  cssRss(j1,j2,j3) = diffss2(cssR,j1,j2,j3,drR);
  cssRrs(j1,j2,j3) = diffrs2(cssR,j1,j2,j3,drR);
   cssRr(j1,j2,j3) =  diffr2(cssR,j1,j2,j3,drR);
   cssRs(j1,j2,j3) =  diffs2(cssR,j1,j2,j3,drR);
  crRrr(j1,j2,j3) = diffrr2(crR,j1,j2,j3,drR);
  crRss(j1,j2,j3) = diffss2(crR,j1,j2,j3,drR);
  crRrs(j1,j2,j3) = diffrs2(crR,j1,j2,j3,drR);
   crRr(j1,j2,j3) =  diffr2(crR,j1,j2,j3,drR);
   crRs(j1,j2,j3) =  diffs2(crR,j1,j2,j3,drR);
  csRrr(j1,j2,j3) = diffrr2(csR,j1,j2,j3,drR);
  csRss(j1,j2,j3) = diffss2(csR,j1,j2,j3,drR);
  csRrs(j1,j2,j3) = diffrs2(csR,j1,j2,j3,drR);
   csRr(j1,j2,j3) =  diffr2(csR,j1,j2,j3,drR);
   csRs(j1,j2,j3) =  diffs2(csR,j1,j2,j3,drR);
} // end FOR_3IJ

// ----- Evaluate coefficients of L1 on the boundary -----

// 2*numGhost = 4 extra ghost needed for fourth order 
extra=2*numGhost;
getBoundaryIndex(mg.indexRange(), side ,axis ,Ib1,Ib2,Ib3,extra);
getBoundaryIndex(mg2.indexRange(),side2,axis2,Jb1,Jb2,Jb3,extra);
RealArray L1ICoeff(Ib1,Ib2,Ib3), L1ICoeffs(Ib1,Ib2,Ib3), L1ICoeffss(Ib1,Ib2,Ib3), L1ICoeffsss(Ib1,Ib2,Ib3);
RealArray L1rCoeff(Ib1,Ib2,Ib3), L1rCoeffs(Ib1,Ib2,Ib3), L1rCoeffss(Ib1,Ib2,Ib3), L1rCoeffsss(Ib1,Ib2,Ib3);
RealArray L1sCoeff(Ib1,Ib2,Ib3), L1sCoeffs(Ib1,Ib2,Ib3), L1sCoeffss(Ib1,Ib2,Ib3), L1sCoeffsss(Ib1,Ib2,Ib3);
FOR_3IJ(i1,i2,i3,Ib1,Ib2,Ib3,j1,j2,j3,Jb1,Jb2,Jb3)
{
  const Real n1R = -n1L(i1,i2,i3), n2R = -n2L(i1,i2,i3);    // for advection, use right normal 
  getAdvectionCoefficients(i1,i2,i3,u1DL,u2DL, j1,j2,j3,u1DR,u2DR );
  L1ICoeff(i1,i2,i3) = ((-1.*KLR*u2DL+1.*u2DR)*n2R+(-1.*KLR*u1DL+1.*u1DR)*n1R)/b1R(j1,j2,j3);
  L1rCoeff(i1,i2,i3) = 1.*KLR*b1L(i1,i2,i3)/b1R(j1,j2,j3);
  L1sCoeff(i1,i2,i3) = (1.*KLR*b2L(i1,i2,i3)-1.*b2R(j1,j2,j3))/b1R(j1,j2,j3);
} // end FOR_3D

// ----- Evaluate coefficients of L2 on the boundary -----
extra=numGhost;
getBoundaryIndex(mg.indexRange(), side,  axis,Ib1,Ib2,Ib3,extra);
getBoundaryIndex(mg2.indexRange(),side2,axis2,Jb1,Jb2,Jb3,extra);
RealArray L2rrCoeff(Ib1,Ib2,Ib3), L2rrCoeffs(Ib1,Ib2,Ib3), L2rrCoeffss(Ib1,Ib2,Ib3);
RealArray L2rsCoeff(Ib1,Ib2,Ib3), L2rsCoeffs(Ib1,Ib2,Ib3), L2rsCoeffss(Ib1,Ib2,Ib3);
RealArray L2ssCoeff(Ib1,Ib2,Ib3), L2ssCoeffs(Ib1,Ib2,Ib3), L2ssCoeffss(Ib1,Ib2,Ib3);
RealArray L2rCoeff(Ib1,Ib2,Ib3), L2rCoeffs(Ib1,Ib2,Ib3), L2rCoeffss(Ib1,Ib2,Ib3);
RealArray L2sCoeff(Ib1,Ib2,Ib3), L2sCoeffs(Ib1,Ib2,Ib3), L2sCoeffss(Ib1,Ib2,Ib3);
RealArray L2ICoeff(Ib1,Ib2,Ib3), L2ICoeffs(Ib1,Ib2,Ib3), L2ICoeffss(Ib1,Ib2,Ib3);
FOR_3IJ(i1,i2,i3,Ib1,Ib2,Ib3,j1,j2,j3,Jb1,Jb2,Jb3)
{
  // Compute tangential derivatives of L1 - label tangential directions a and b (for 3D)
  L1rCoeffs(i1,i2,i3) = diffa(L1rCoeff,i1,i2,i3,drL);
  L1sCoeffs(i1,i2,i3) = diffa(L1sCoeff,i1,i2,i3,drL);
  L1ICoeffs(i1,i2,i3) = diffa(L1ICoeff,i1,i2,i3,drL);
  // Compute L2 coefficients
  L2rrCoeff(i1,i2,i3) = 1.*DLR*crrL(i1,i2,i3)/crrR(j1,j2,j3);
  L2rsCoeff(i1,i2,i3) = (1.*DLR*crsL(i1,i2,i3)-1.*crsR(j1,j2,j3)*L1rCoeff(i1,i2,i3))/crrR(j1,j2,j3);
  L2ssCoeff(i1,i2,i3) = 1.*(DLR*cssL(i1,i2,i3)-crsR(j1,j2,j3)*L1sCoeff(i1,i2,i3)-cssR(j1,j2,j3))/crrR(j1,j2,j3);
  L2rCoeff(i1,i2,i3) = 1.*(DLR*crL(i1,i2,i3)-crR(j1,j2,j3)*L1rCoeff(i1,i2,i3)-crsR(j1,j2,j3)*L1rCoeffs(i1,i2,i3))/crrR(j1,j2,j3);
  L2sCoeff(i1,i2,i3) = ((-1.*L1ICoeff(i1,i2,i3)-1.*L1sCoeffs(i1,i2,i3))*crsR(j1,j2,j3)+1.*DLR*csL(i1,i2,i3)-1.*crR(j1,j2,j3)*L1sCoeff(i1,i2,i3)-1.*csR(j1,j2,j3))/crrR(j1,j2,j3);
  L2ICoeff(i1,i2,i3) = (-1.*crR(j1,j2,j3)*L1ICoeff(i1,i2,i3)-1.*crsR(j1,j2,j3)*L1ICoeffs(i1,i2,i3))/crrR(j1,j2,j3);
} // end FOR_3D

// ----- Evaluate coefficients of L3 on the boundary -----
extra=numGhost-1;
getBoundaryIndex(mg.indexRange(),side,axis,Ib1,Ib2,Ib3,extra);
getBoundaryIndex(mg2.indexRange(),side2,axis2,Jb1,Jb2,Jb3,extra);
RealArray L3rrrCoeff(Ib1,Ib2,Ib3), L3rrrCoeffs(Ib1,Ib2,Ib3);
RealArray L3rrsCoeff(Ib1,Ib2,Ib3), L3rrsCoeffs(Ib1,Ib2,Ib3);
RealArray L3rssCoeff(Ib1,Ib2,Ib3), L3rssCoeffs(Ib1,Ib2,Ib3);
RealArray L3sssCoeff(Ib1,Ib2,Ib3), L3sssCoeffs(Ib1,Ib2,Ib3);
RealArray L3rrCoeff(Ib1,Ib2,Ib3), L3rrCoeffs(Ib1,Ib2,Ib3);
RealArray L3rsCoeff(Ib1,Ib2,Ib3), L3rsCoeffs(Ib1,Ib2,Ib3);
RealArray L3ssCoeff(Ib1,Ib2,Ib3), L3ssCoeffs(Ib1,Ib2,Ib3);
RealArray L3rCoeff(Ib1,Ib2,Ib3), L3rCoeffs(Ib1,Ib2,Ib3);
RealArray L3sCoeff(Ib1,Ib2,Ib3), L3sCoeffs(Ib1,Ib2,Ib3);
RealArray L3ICoeff(Ib1,Ib2,Ib3), L3ICoeffs(Ib1,Ib2,Ib3);
FOR_3IJ(i1,i2,i3,Ib1,Ib2,Ib3,j1,j2,j3,Jb1,Jb2,Jb3)
{
  const Real n1R = -n1L(i1,i2,i3), n2R = -n2L(i1,i2,i3);    // for advection, use right normal 
  getAdvectionCoefficients(i1,i2,i3,u1DL,u2DL, j1,j2,j3,u1DR,u2DR );
  // Compute tangential derivatives of L1 
  L1rCoeffss(i1,i2,i3) = diffaa2(L1rCoeff,i1,i2,i3,drL);
  L1sCoeffss(i1,i2,i3) = diffaa2(L1sCoeff,i1,i2,i3,drL);
  L1ICoeffss(i1,i2,i3) = diffaa2(L1ICoeff,i1,i2,i3,drL);
  // Compute tangential derivatives of L2
  L2rrCoeffs(i1,i2,i3) = diffa2(L2rrCoeff,i1,i2,i3,drL);
  L2rsCoeffs(i1,i2,i3) = diffa2(L2rsCoeff,i1,i2,i3,drL);
  L2ssCoeffs(i1,i2,i3) = diffa2(L2ssCoeff,i1,i2,i3,drL);
  L2rCoeffs(i1,i2,i3)  = diffa2(L2rCoeff, i1,i2,i3,drL);
  L2sCoeffs(i1,i2,i3)  = diffa2(L2sCoeff, i1,i2,i3,drL);
  L2ICoeffs(i1,i2,i3)  = diffa2(L2ICoeff, i1,i2,i3,drL);
  // Compute L3 coefficients
  L3rrrCoeff(i1,i2,i3) = 1.*DLR*KLR*b1L(i1,i2,i3)*crrL(i1,i2,i3)/b1R(j1,j2,j3)/crrR(j1,j2,j3);
  L3rrsCoeff(i1,i2,i3) = ((-1.*crsR(j1,j2,j3)*b1R(j1,j2,j3)-1.*crrR(j1,j2,j3)*b2R(j1,j2,j3))*L2rrCoeff(i1,i2,i3)+(1.*crsL(i1,i2,i3)*b1L(i1,i2,i3)+1.*b2L(i1,i2,i3)*crrL(i1,i2,i3))*KLR*DLR)/b1R(j1,j2,j3)/crrR(j1,j2,j3);
  L3rssCoeff(i1,i2,i3) = 1.*(DLR*KLR*b1L(i1,i2,i3)*cssL(i1,i2,i3)+DLR*KLR*b2L(i1,i2,i3)*crsL(i1,i2,i3)-b1R(j1,j2,j3)*crsR(j1,j2,j3)*L2rsCoeff(i1,i2,i3)-b1R(j1,j2,j3)*cssR(j1,j2,j3)*L1rCoeff(i1,i2,i3)-b2R(j1,j2,j3)*crrR(j1,j2,j3)*L2rsCoeff(i1,i2,i3)-b2R(j1,j2,j3)*crsR(j1,j2,j3)*L1rCoeff(i1,i2,i3))/b1R(j1,j2,j3)/crrR(j1,j2,j3);
  L3sssCoeff(i1,i2,i3) = ((-1.*L2ssCoeff(i1,i2,i3)*crsR(j1,j2,j3)-1.*L1sCoeff(i1,i2,i3)*cssR(j1,j2,j3))*b1R(j1,j2,j3)+(-1.*crrR(j1,j2,j3)*L2ssCoeff(i1,i2,i3)-1.*crsR(j1,j2,j3)*L1sCoeff(i1,i2,i3)-1.*cssR(j1,j2,j3))*b2R(j1,j2,j3)+1.*DLR*KLR*b2L(i1,i2,i3)*cssL(i1,i2,i3))/b1R(j1,j2,j3)/crrR(j1,j2,j3);
  L3rrCoeff(i1,i2,i3) = (((1.*n1R*u1DR+1.*n2R*u2DR)*crrR(j1,j2,j3)+(-1.*crR(j1,j2,j3)-1.*crrRr(j1,j2,j3))*b1R(j1,j2,j3)-1.*crrRs(j1,j2,j3)*b2R(j1,j2,j3))*L2rrCoeff(i1,i2,i3)+((-1.*n1R*u1DL-1.*n2R*u2DL)*crrL(i1,i2,i3)+(1.*crL(i1,i2,i3)+1.*crrLr(i1,i2,i3))*b1L(i1,i2,i3)+1.*b2L(i1,i2,i3)*crrLs(i1,i2,i3))*KLR*DLR-1.*b1R(j1,j2,j3)*crsR(j1,j2,j3)*L2rrCoeffs(i1,i2,i3)-1.*b2R(j1,j2,j3)*crrR(j1,j2,j3)*L2rrCoeffs(i1,i2,i3))/b1R(j1,j2,j3)/crrR(j1,j2,j3);
  L3rsCoeff(i1,i2,i3) = (((-1.*crsRr(j1,j2,j3)-1.*csR(j1,j2,j3))*L1rCoeff(i1,i2,i3)+(-1.*crR(j1,j2,j3)-1.*crrRr(j1,j2,j3))*L2rsCoeff(i1,i2,i3)+(-1.*L2rsCoeffs(i1,i2,i3)-1.*L2rCoeff(i1,i2,i3))*crsR(j1,j2,j3)-2.*cssR(j1,j2,j3)*L1rCoeffs(i1,i2,i3))*b1R(j1,j2,j3)+((-1.*crR(j1,j2,j3)-1.*crsRs(j1,j2,j3))*b2R(j1,j2,j3)+crsR(j1,j2,j3)*(n1R*u1DR+n2R*u2DR))*L1rCoeff(i1,i2,i3)+(-1.*crsL(i1,i2,i3)*n1R*u1DL-1.*crsL(i1,i2,i3)*n2R*u2DL+(crsLr(i1,i2,i3)+csL(i1,i2,i3))*b1L(i1,i2,i3)+b2L(i1,i2,i3)*(crsLs(i1,i2,i3)+crL(i1,i2,i3)))*KLR*DLR+(-1.*crrRs(j1,j2,j3)*L2rsCoeff(i1,i2,i3)+(-1.*L2rsCoeffs(i1,i2,i3)-1.*L2rCoeff(i1,i2,i3))*crrR(j1,j2,j3)-2.*crsR(j1,j2,j3)*L1rCoeffs(i1,i2,i3))*b2R(j1,j2,j3)+crrR(j1,j2,j3)*L2rsCoeff(i1,i2,i3)*(n1R*u1DR+n2R*u2DR))/b1R(j1,j2,j3)/crrR(j1,j2,j3);
  L3ssCoeff(i1,i2,i3) = (((-1.*crsRr(j1,j2,j3)-1.*csR(j1,j2,j3))*L1sCoeff(i1,i2,i3)+(-1.*L2ssCoeffs(i1,i2,i3)-1.*L2sCoeff(i1,i2,i3))*crsR(j1,j2,j3)+(-1.*crR(j1,j2,j3)-1.*crrRr(j1,j2,j3))*L2ssCoeff(i1,i2,i3)+(-2.*L1sCoeffs(i1,i2,i3)-1.*L1ICoeff(i1,i2,i3))*cssR(j1,j2,j3)-1.*cssRr(j1,j2,j3))*b1R(j1,j2,j3)+((-1.*crR(j1,j2,j3)-1.*crsRs(j1,j2,j3))*L1sCoeff(i1,i2,i3)+(-2.*L1sCoeffs(i1,i2,i3)-1.*L1ICoeff(i1,i2,i3))*crsR(j1,j2,j3)-1.*crrRs(j1,j2,j3)*L2ssCoeff(i1,i2,i3)+(-1.*L2ssCoeffs(i1,i2,i3)-1.*L2sCoeff(i1,i2,i3))*crrR(j1,j2,j3)-1.*cssRs(j1,j2,j3)-1.*csR(j1,j2,j3))*b2R(j1,j2,j3)+crsR(j1,j2,j3)*(n1R*u1DR+n2R*u2DR)*L1sCoeff(i1,i2,i3)+crrR(j1,j2,j3)*(n1R*u1DR+n2R*u2DR)*L2ssCoeff(i1,i2,i3)+(-1.*cssL(i1,i2,i3)*n1R*u1DL-1.*cssL(i1,i2,i3)*n2R*u2DL+(cssLs(i1,i2,i3)+csL(i1,i2,i3))*b2L(i1,i2,i3)+cssLr(i1,i2,i3)*b1L(i1,i2,i3))*KLR*DLR+cssR(j1,j2,j3)*(n1R*u1DR+n2R*u2DR))/b1R(j1,j2,j3)/crrR(j1,j2,j3);
  L3rCoeff(i1,i2,i3) = (((-1.*crsRr(j1,j2,j3)-1.*csR(j1,j2,j3))*L1rCoeffs(i1,i2,i3)+(-1.*crR(j1,j2,j3)-1.*crrRr(j1,j2,j3))*L2rCoeff(i1,i2,i3)-1.*crsR(j1,j2,j3)*L2rCoeffs(i1,i2,i3)-1.*cssR(j1,j2,j3)*L1rCoeffss(i1,i2,i3)-1.*crRr(j1,j2,j3)*L1rCoeff(i1,i2,i3))*b1R(j1,j2,j3)+((-1.*crR(j1,j2,j3)-1.*crsRs(j1,j2,j3))*b2R(j1,j2,j3)+(1.*n1R*u1DR+1.*n2R*u2DR)*crsR(j1,j2,j3))*L1rCoeffs(i1,i2,i3)+(-1.*crsR(j1,j2,j3)*L1rCoeffss(i1,i2,i3)-1.*crRs(j1,j2,j3)*L1rCoeff(i1,i2,i3)-1.*crrRs(j1,j2,j3)*L2rCoeff(i1,i2,i3)-1.*L2rCoeffs(i1,i2,i3)*crrR(j1,j2,j3))*b2R(j1,j2,j3)+(1.*n1R*u1DR+1.*n2R*u2DR)*crrR(j1,j2,j3)*L2rCoeff(i1,i2,i3)+(1.*n1R*u1DR+1.*n2R*u2DR)*crR(j1,j2,j3)*L1rCoeff(i1,i2,i3)+(-1.*crL(i1,i2,i3)*n1R*u1DL-1.*crL(i1,i2,i3)*n2R*u2DL+1.*crLr(i1,i2,i3)*b1L(i1,i2,i3)+1.*crLs(i1,i2,i3)*b2L(i1,i2,i3))*KLR*DLR)/b1R(j1,j2,j3)/crrR(j1,j2,j3);
  L3sCoeff(i1,i2,i3) = (((-1.*L2ICoeff(i1,i2,i3)-1.*L2sCoeffs(i1,i2,i3))*crsR(j1,j2,j3)+(-1.*crsRr(j1,j2,j3)-1.*csR(j1,j2,j3))*L1sCoeffs(i1,i2,i3)+(-1.*crsRr(j1,j2,j3)-1.*csR(j1,j2,j3))*L1ICoeff(i1,i2,i3)+(-1.*crR(j1,j2,j3)-1.*crrRr(j1,j2,j3))*L2sCoeff(i1,i2,i3)-1.*csRr(j1,j2,j3)-1.*cssR(j1,j2,j3)*L1sCoeffss(i1,i2,i3)-2.*cssR(j1,j2,j3)*L1ICoeffs(i1,i2,i3)-1.*crRr(j1,j2,j3)*L1sCoeff(i1,i2,i3))*b1R(j1,j2,j3)+((-1.*L1sCoeffss(i1,i2,i3)-2.*L1ICoeffs(i1,i2,i3))*crsR(j1,j2,j3)+(-1.*crR(j1,j2,j3)-1.*crsRs(j1,j2,j3))*L1sCoeffs(i1,i2,i3)+(-1.*crR(j1,j2,j3)-1.*crsRs(j1,j2,j3))*L1ICoeff(i1,i2,i3)-1.*crrRs(j1,j2,j3)*L2sCoeff(i1,i2,i3)+(-1.*L2ICoeff(i1,i2,i3)-1.*L2sCoeffs(i1,i2,i3))*crrR(j1,j2,j3)-1.*csRs(j1,j2,j3)-1.*crRs(j1,j2,j3)*L1sCoeff(i1,i2,i3))*b2R(j1,j2,j3)+(L1ICoeff(i1,i2,i3)+L1sCoeffs(i1,i2,i3))*(n1R*u1DR+n2R*u2DR)*crsR(j1,j2,j3)+(u1DR*L2sCoeff(i1,i2,i3)*crrR(j1,j2,j3)+(crR(j1,j2,j3)*L1sCoeff(i1,i2,i3)+csR(j1,j2,j3))*u1DR-1.*csL(i1,i2,i3)*u1DL*KLR*DLR)*n1R+(u2DR*L2sCoeff(i1,i2,i3)*crrR(j1,j2,j3)+(crR(j1,j2,j3)*L1sCoeff(i1,i2,i3)+csR(j1,j2,j3))*u2DR-1.*csL(i1,i2,i3)*u2DL*KLR*DLR)*n2R+DLR*KLR*(b1L(i1,i2,i3)*csLr(i1,i2,i3)+b2L(i1,i2,i3)*csLs(i1,i2,i3)))/b1R(j1,j2,j3)/crrR(j1,j2,j3);
  L3ICoeff(i1,i2,i3) = (((-1.*crsRr(j1,j2,j3)-1.*csR(j1,j2,j3))*L1ICoeffs(i1,i2,i3)+(-1.*crR(j1,j2,j3)-1.*crrRr(j1,j2,j3))*L2ICoeff(i1,i2,i3)-1.*crRr(j1,j2,j3)*L1ICoeff(i1,i2,i3)-1.*crsR(j1,j2,j3)*L2ICoeffs(i1,i2,i3)-1.*cssR(j1,j2,j3)*L1ICoeffss(i1,i2,i3))*b1R(j1,j2,j3)+((-1.*crR(j1,j2,j3)-1.*crsRs(j1,j2,j3))*b2R(j1,j2,j3)+(1.*n1R*u1DR+1.*n2R*u2DR)*crsR(j1,j2,j3))*L1ICoeffs(i1,i2,i3)+(-1.*L1ICoeff(i1,i2,i3)*crRs(j1,j2,j3)-1.*L2ICoeffs(i1,i2,i3)*crrR(j1,j2,j3)-1.*crrRs(j1,j2,j3)*L2ICoeff(i1,i2,i3)-1.*crsR(j1,j2,j3)*L1ICoeffss(i1,i2,i3))*b2R(j1,j2,j3)+(1.*n1R*u1DR+1.*n2R*u2DR)*crrR(j1,j2,j3)*L2ICoeff(i1,i2,i3)+(1.*n1R*u1DR+1.*n2R*u2DR)*crR(j1,j2,j3)*L1ICoeff(i1,i2,i3))/b1R(j1,j2,j3)/crrR(j1,j2,j3);
} // end FOR_3D

// ----- Evaluate coefficients of L4 on the boundary -----
extra=0;
getBoundaryIndex(mg.indexRange(), side,axis  ,Ib1,Ib2,Ib3,extra);
getBoundaryIndex(mg2.indexRange(),side2,axis2,Jb1,Jb2,Jb3,extra);
RealArray L4rrrrCoeff(Ib1,Ib2,Ib3);
RealArray L4rrrsCoeff(Ib1,Ib2,Ib3);
RealArray L4rrssCoeff(Ib1,Ib2,Ib3);
RealArray L4rsssCoeff(Ib1,Ib2,Ib3);
RealArray L4ssssCoeff(Ib1,Ib2,Ib3);
RealArray L4rrrCoeff(Ib1,Ib2,Ib3);
RealArray L4rrsCoeff(Ib1,Ib2,Ib3);
RealArray L4rssCoeff(Ib1,Ib2,Ib3);
RealArray L4sssCoeff(Ib1,Ib2,Ib3);
RealArray L4rrCoeff(Ib1,Ib2,Ib3);
RealArray L4rsCoeff(Ib1,Ib2,Ib3);
RealArray L4ssCoeff(Ib1,Ib2,Ib3);
RealArray L4rCoeff(Ib1,Ib2,Ib3);
RealArray L4sCoeff(Ib1,Ib2,Ib3);
RealArray L4ICoeff(Ib1,Ib2,Ib3);
FOR_3IJ(i1,i2,i3,Ib1,Ib2,Ib3,j1,j2,j3,Jb1,Jb2,Jb3)
{
  // Compute tangential derivatives of L1
  L1rCoeffsss(i1,i2,i3) = diffaaa(L1rCoeff,i1,i2,i3,drL);
  L1sCoeffsss(i1,i2,i3) = diffaaa(L1sCoeff,i1,i2,i3,drL);
  L1ICoeffsss(i1,i2,i3) = diffaaa(L1ICoeff,i1,i2,i3,drL);
  // Compute tangential derivatives of L2
  L2rrCoeffss(i1,i2,i3) = diffaa2(L2rrCoeff,i1,i2,i3,drL);
  L2rsCoeffss(i1,i2,i3) = diffaa2(L2rsCoeff,i1,i2,i3,drL);
  L2ssCoeffss(i1,i2,i3) = diffaa2(L2ssCoeff,i1,i2,i3,drL);
  L2rCoeffss(i1,i2,i3)  = diffaa2(L2rCoeff, i1,i2,i3,drL);
  L2sCoeffss(i1,i2,i3)  = diffaa2(L2sCoeff, i1,i2,i3,drL);
  L2ICoeffss(i1,i2,i3)  = diffaa2(L2ICoeff, i1,i2,i3,drL);
  // Compute tangential derivatives of L3
  L3rrrCoeffs(i1,i2,i3) = diffa2(L3rrrCoeff,i1,i2,i3,drL);
  L3rrsCoeffs(i1,i2,i3) = diffa2(L3rrsCoeff,i1,i2,i3,drL);
  L3rssCoeffs(i1,i2,i3) = diffa2(L3rssCoeff,i1,i2,i3,drL);
  L3sssCoeffs(i1,i2,i3) = diffa2(L3sssCoeff,i1,i2,i3,drL);
  L3rrCoeffs(i1,i2,i3)  = diffa2(L3rrCoeff, i1,i2,i3,drL);
  L3rsCoeffs(i1,i2,i3)  = diffa2(L3rsCoeff, i1,i2,i3,drL);
  L3ssCoeffs(i1,i2,i3)  = diffa2(L3ssCoeff, i1,i2,i3,drL);
  L3rCoeffs(i1,i2,i3)   = diffa2(L3rCoeff,  i1,i2,i3,drL);
  L3sCoeffs(i1,i2,i3)   = diffa2(L3sCoeff,  i1,i2,i3,drL);
  L3ICoeffs(i1,i2,i3)   = diffa2(L3ICoeff,  i1,i2,i3,drL);
  // Compute L4 coefficients
  L4rrrrCoeff(i1,i2,i3) = 1.*pow(DLR,2.)*pow(crrL(i1,i2,i3),2.)/pow(crrR(j1,j2,j3),2.);
  L4rrrsCoeff(i1,i2,i3) = (2.*pow(DLR,2.)*crrL(i1,i2,i3)*crsL(i1,i2,i3)-2.*crrR(j1,j2,j3)*crsR(j1,j2,j3)*L3rrrCoeff(i1,i2,i3))/pow(crrR(j1,j2,j3),2.);
  L4rrssCoeff(i1,i2,i3) = (2.*pow(DLR,2.)*crrL(i1,i2,i3)*cssL(i1,i2,i3)+pow(DLR,2.)*pow(crsL(i1,i2,i3),2.)-2.*crrR(j1,j2,j3)*crsR(j1,j2,j3)*L3rrsCoeff(i1,i2,i3)-2.*crrR(j1,j2,j3)*cssR(j1,j2,j3)*L2rrCoeff(i1,i2,i3)-pow(crsR(j1,j2,j3),2.)*L2rrCoeff(i1,i2,i3))/pow(crrR(j1,j2,j3),2.);
  L4rsssCoeff(i1,i2,i3) = (2.*pow(DLR,2.)*crsL(i1,i2,i3)*cssL(i1,i2,i3)-2.*crrR(j1,j2,j3)*crsR(j1,j2,j3)*L3rssCoeff(i1,i2,i3)-2.*crrR(j1,j2,j3)*cssR(j1,j2,j3)*L2rsCoeff(i1,i2,i3)-pow(crsR(j1,j2,j3),2.)*L2rsCoeff(i1,i2,i3)-2.*crsR(j1,j2,j3)*cssR(j1,j2,j3)*L1rCoeff(i1,i2,i3))/pow(crrR(j1,j2,j3),2.);
  L4ssssCoeff(i1,i2,i3) = (pow(DLR,2.)*pow(cssL(i1,i2,i3),2.)-2.*crrR(j1,j2,j3)*crsR(j1,j2,j3)*L3sssCoeff(i1,i2,i3)-2.*crrR(j1,j2,j3)*cssR(j1,j2,j3)*L2ssCoeff(i1,i2,i3)-pow(crsR(j1,j2,j3),2.)*L2ssCoeff(i1,i2,i3)-2.*crsR(j1,j2,j3)*cssR(j1,j2,j3)*L1sCoeff(i1,i2,i3)-pow(cssR(j1,j2,j3),2.))/pow(crrR(j1,j2,j3),2.);
  L4rrrCoeff(i1,i2,i3) = (((-2.*crR(j1,j2,j3)-2.*crrRr(j1,j2,j3))*L3rrrCoeff(i1,i2,i3)-2.*crsR(j1,j2,j3)*L3rrrCoeffs(i1,i2,i3))*crrR(j1,j2,j3)-1.*crrRs(j1,j2,j3)*crsR(j1,j2,j3)*L3rrrCoeff(i1,i2,i3)+((2.*crL(i1,i2,i3)+2.*crrLr(i1,i2,i3))*crrL(i1,i2,i3)+crrLs(i1,i2,i3)*crsL(i1,i2,i3))*pow(DLR,2.))/pow(crrR(j1,j2,j3),2.);
  L4rrsCoeff(i1,i2,i3) = (((-2.*L3rrCoeff(i1,i2,i3)-2.*L3rrsCoeffs(i1,i2,i3))*crsR(j1,j2,j3)+(-2.*csR(j1,j2,j3)-2.*crsRr(j1,j2,j3))*L2rrCoeff(i1,i2,i3)+(-2.*crR(j1,j2,j3)-2.*crrRr(j1,j2,j3))*L3rrsCoeff(i1,i2,i3)-4.*cssR(j1,j2,j3)*L2rrCoeffs(i1,i2,i3))*crrR(j1,j2,j3)-2.*pow(crsR(j1,j2,j3),2.)*L2rrCoeffs(i1,i2,i3)+((-2.*crR(j1,j2,j3)-1.*crrRr(j1,j2,j3)-1.*crsRs(j1,j2,j3))*L2rrCoeff(i1,i2,i3)-1.*crrRs(j1,j2,j3)*L3rrsCoeff(i1,i2,i3))*crsR(j1,j2,j3)-2.*crrRs(j1,j2,j3)*cssR(j1,j2,j3)*L2rrCoeff(i1,i2,i3)+((2.*crL(i1,i2,i3)+crrLr(i1,i2,i3)+crsLs(i1,i2,i3))*crsL(i1,i2,i3)+(2.*crsLr(i1,i2,i3)+2.*csL(i1,i2,i3))*crrL(i1,i2,i3)+2.*crrLs(i1,i2,i3)*cssL(i1,i2,i3))*pow(DLR,2.))/pow(crrR(j1,j2,j3),2.);
  L4rssCoeff(i1,i2,i3) = ((-1.*L2rCoeff(i1,i2,i3)-2.*L2rsCoeffs(i1,i2,i3))*pow(crsR(j1,j2,j3),2.)+((-2.*L3rsCoeff(i1,i2,i3)-2.*L3rssCoeffs(i1,i2,i3))*crrR(j1,j2,j3)+(-2.*crR(j1,j2,j3)-1.*crrRr(j1,j2,j3)-1.*crsRs(j1,j2,j3))*L2rsCoeff(i1,i2,i3)+(-1.*cssRs(j1,j2,j3)-1.*crsRr(j1,j2,j3)-2.*csR(j1,j2,j3))*L1rCoeff(i1,i2,i3)-1.*L3rssCoeff(i1,i2,i3)*crrRs(j1,j2,j3)-6.*cssR(j1,j2,j3)*L1rCoeffs(i1,i2,i3))*crsR(j1,j2,j3)+((-2.*csR(j1,j2,j3)-2.*crsRr(j1,j2,j3))*L2rsCoeff(i1,i2,i3)-2.*cssRr(j1,j2,j3)*L1rCoeff(i1,i2,i3)+(-2.*L2rCoeff(i1,i2,i3)-4.*L2rsCoeffs(i1,i2,i3))*cssR(j1,j2,j3)+(-2.*crR(j1,j2,j3)-2.*crrRr(j1,j2,j3))*L3rssCoeff(i1,i2,i3))*crrR(j1,j2,j3)-2.*crrRs(j1,j2,j3)*cssR(j1,j2,j3)*L2rsCoeff(i1,i2,i3)+(-2.*crR(j1,j2,j3)-2.*crsRs(j1,j2,j3))*cssR(j1,j2,j3)*L1rCoeff(i1,i2,i3)+((crsLr(i1,i2,i3)+cssLs(i1,i2,i3)+2.*csL(i1,i2,i3))*crsL(i1,i2,i3)+(2.*crL(i1,i2,i3)+2.*crsLs(i1,i2,i3))*cssL(i1,i2,i3)+2.*crrL(i1,i2,i3)*cssLr(i1,i2,i3))*pow(DLR,2.))/pow(crrR(j1,j2,j3),2.);
  L4sssCoeff(i1,i2,i3) = ((-1.*L2sCoeff(i1,i2,i3)-2.*L2ssCoeffs(i1,i2,i3))*pow(crsR(j1,j2,j3),2.)+((-2.*L3ssCoeff(i1,i2,i3)-2.*L3sssCoeffs(i1,i2,i3))*crrR(j1,j2,j3)+(-6.*L1sCoeffs(i1,i2,i3)-2.*L1ICoeff(i1,i2,i3))*cssR(j1,j2,j3)+(-2.*crR(j1,j2,j3)-1.*crrRr(j1,j2,j3)-1.*crsRs(j1,j2,j3))*L2ssCoeff(i1,i2,i3)+(-1.*cssRs(j1,j2,j3)-1.*crsRr(j1,j2,j3)-2.*csR(j1,j2,j3))*L1sCoeff(i1,i2,i3)-1.*cssRr(j1,j2,j3)-1.*crrRs(j1,j2,j3)*L3sssCoeff(i1,i2,i3))*crsR(j1,j2,j3)+((-2.*L2sCoeff(i1,i2,i3)-4.*L2ssCoeffs(i1,i2,i3))*cssR(j1,j2,j3)+(-2.*csR(j1,j2,j3)-2.*crsRr(j1,j2,j3))*L2ssCoeff(i1,i2,i3)-2.*cssRr(j1,j2,j3)*L1sCoeff(i1,i2,i3)+(-2.*crR(j1,j2,j3)-2.*crrRr(j1,j2,j3))*L3sssCoeff(i1,i2,i3))*crrR(j1,j2,j3)+(-2.*crrRs(j1,j2,j3)*L2ssCoeff(i1,i2,i3)+(-2.*crR(j1,j2,j3)-2.*crsRs(j1,j2,j3))*L1sCoeff(i1,i2,i3)-2.*cssRs(j1,j2,j3)-2.*csR(j1,j2,j3))*cssR(j1,j2,j3)+((2.*csL(i1,i2,i3)+2.*cssLs(i1,i2,i3))*cssL(i1,i2,i3)+crsL(i1,i2,i3)*cssLr(i1,i2,i3))*pow(DLR,2.))/pow(crrR(j1,j2,j3),2.);
  L4rrCoeff(i1,i2,i3) = (((-2.*crRr(j1,j2,j3)-1.*crrRrr(j1,j2,j3))*L2rrCoeff(i1,i2,i3)-2.*L3rrCoeffs(i1,i2,i3)*crsR(j1,j2,j3)+(-2.*csR(j1,j2,j3)-2.*crsRr(j1,j2,j3))*L2rrCoeffs(i1,i2,i3)-2.*crR(j1,j2,j3)*L3rrCoeff(i1,i2,i3)-2.*crrRr(j1,j2,j3)*L3rrCoeff(i1,i2,i3)-2.*cssR(j1,j2,j3)*L2rrCoeffss(i1,i2,i3))*crrR(j1,j2,j3)+((-1.*crrRrs(j1,j2,j3)-1.*crRs(j1,j2,j3))*crsR(j1,j2,j3)-1.*crrRr(j1,j2,j3)*crR(j1,j2,j3)-1.*cssR(j1,j2,j3)*crrRss(j1,j2,j3)-1.*csR(j1,j2,j3)*crrRs(j1,j2,j3)-1.*pow(crR(j1,j2,j3),2.))*L2rrCoeff(i1,i2,i3)-1.*pow(crsR(j1,j2,j3),2.)*L2rrCoeffss(i1,i2,i3)+((-2.*crR(j1,j2,j3)-1.*crrRr(j1,j2,j3)-1.*crsRs(j1,j2,j3))*L2rrCoeffs(i1,i2,i3)-1.*crrRs(j1,j2,j3)*L3rrCoeff(i1,i2,i3))*crsR(j1,j2,j3)+((crrLrr(i1,i2,i3)+2.*crLr(i1,i2,i3))*crrL(i1,i2,i3)+(crrLrs(i1,i2,i3)+crLs(i1,i2,i3))*crsL(i1,i2,i3)+cssL(i1,i2,i3)*crrLss(i1,i2,i3)+crrLr(i1,i2,i3)*crL(i1,i2,i3)+csL(i1,i2,i3)*crrLs(i1,i2,i3)+pow(crL(i1,i2,i3),2.))*pow(DLR,2.)-2.*cssR(j1,j2,j3)*crrRs(j1,j2,j3)*L2rrCoeffs(i1,i2,i3))/pow(crrR(j1,j2,j3),2.);
  L4rsCoeff(i1,i2,i3) = ((-2.*L2rCoeffs(i1,i2,i3)-1.*L2rsCoeffss(i1,i2,i3))*pow(crsR(j1,j2,j3),2.)+((-2.*L3rsCoeffs(i1,i2,i3)-2.*L3rCoeff(i1,i2,i3))*crrR(j1,j2,j3)+(-1.*crsRrs(j1,j2,j3)-1.*crRr(j1,j2,j3)-1.*csRs(j1,j2,j3))*L1rCoeff(i1,i2,i3)-6.*cssR(j1,j2,j3)*L1rCoeffss(i1,i2,i3)+(-1.*crrRrs(j1,j2,j3)-1.*crRs(j1,j2,j3))*L2rsCoeff(i1,i2,i3)+(-2.*L2rsCoeffs(i1,i2,i3)-2.*L2rCoeff(i1,i2,i3))*crR(j1,j2,j3)+(-1.*crrRr(j1,j2,j3)-1.*crsRs(j1,j2,j3))*L2rsCoeffs(i1,i2,i3)+(-2.*crsRr(j1,j2,j3)-2.*cssRs(j1,j2,j3)-4.*csR(j1,j2,j3))*L1rCoeffs(i1,i2,i3)+(-1.*crrRr(j1,j2,j3)-1.*crsRs(j1,j2,j3))*L2rCoeff(i1,i2,i3)-1.*crrRs(j1,j2,j3)*L3rsCoeff(i1,i2,i3))*crsR(j1,j2,j3)+((-1.*crsRrr(j1,j2,j3)-2.*csRr(j1,j2,j3))*L1rCoeff(i1,i2,i3)+(-2.*L2rsCoeffss(i1,i2,i3)-4.*L2rCoeffs(i1,i2,i3))*cssR(j1,j2,j3)+(-2.*crRr(j1,j2,j3)-1.*crrRrr(j1,j2,j3))*L2rsCoeff(i1,i2,i3)-2.*crR(j1,j2,j3)*L3rsCoeff(i1,i2,i3)+(-2.*csR(j1,j2,j3)-2.*crsRr(j1,j2,j3))*L2rsCoeffs(i1,i2,i3)-4.*cssRr(j1,j2,j3)*L1rCoeffs(i1,i2,i3)+(-2.*csR(j1,j2,j3)-2.*crsRr(j1,j2,j3))*L2rCoeff(i1,i2,i3)-2.*crrRr(j1,j2,j3)*L3rsCoeff(i1,i2,i3))*crrR(j1,j2,j3)+((-1.*crsRss(j1,j2,j3)-2.*crRs(j1,j2,j3))*cssR(j1,j2,j3)+(-1.*crsRr(j1,j2,j3)-2.*csR(j1,j2,j3))*crR(j1,j2,j3)-1.*crsRs(j1,j2,j3)*csR(j1,j2,j3))*L1rCoeff(i1,i2,i3)+(-4.*crR(j1,j2,j3)*L1rCoeffs(i1,i2,i3)-1.*crrRss(j1,j2,j3)*L2rsCoeff(i1,i2,i3)-2.*crrRs(j1,j2,j3)*L2rCoeff(i1,i2,i3)-2.*crrRs(j1,j2,j3)*L2rsCoeffs(i1,i2,i3)-4.*crsRs(j1,j2,j3)*L1rCoeffs(i1,i2,i3))*cssR(j1,j2,j3)+((crsLrs(i1,i2,i3)+crLr(i1,i2,i3)+csLs(i1,i2,i3))*crsL(i1,i2,i3)+(crsLrr(i1,i2,i3)+2.*csLr(i1,i2,i3))*crrL(i1,i2,i3)+(crsLss(i1,i2,i3)+2.*crLs(i1,i2,i3))*cssL(i1,i2,i3)+(crsLr(i1,i2,i3)+2.*csL(i1,i2,i3))*crL(i1,i2,i3)+csL(i1,i2,i3)*crsLs(i1,i2,i3))*pow(DLR,2.)+(-1.*crrRr(j1,j2,j3)*crR(j1,j2,j3)-1.*csR(j1,j2,j3)*crrRs(j1,j2,j3)-1.*pow(crR(j1,j2,j3),2.))*L2rsCoeff(i1,i2,i3))/pow(crrR(j1,j2,j3),2.);
  L4ssCoeff(i1,i2,i3) = ((-1.*L2ssCoeffss(i1,i2,i3)-2.*L2sCoeffs(i1,i2,i3)-1.*L2ICoeff(i1,i2,i3))*pow(crsR(j1,j2,j3),2.)+((-2.*L3ssCoeffs(i1,i2,i3)-2.*L3sCoeff(i1,i2,i3))*crrR(j1,j2,j3)+(-6.*L1ICoeffs(i1,i2,i3)-6.*L1sCoeffss(i1,i2,i3))*cssR(j1,j2,j3)+(-1.*crsRrs(j1,j2,j3)-1.*crRr(j1,j2,j3)-1.*csRs(j1,j2,j3))*L1sCoeff(i1,i2,i3)+(-2.*L2sCoeff(i1,i2,i3)-2.*L2ssCoeffs(i1,i2,i3))*crR(j1,j2,j3)+(-4.*L1sCoeffs(i1,i2,i3)-2.*L1ICoeff(i1,i2,i3))*csR(j1,j2,j3)+(-1.*crrRrs(j1,j2,j3)-1.*crRs(j1,j2,j3))*L2ssCoeff(i1,i2,i3)+(-1.*crrRr(j1,j2,j3)-1.*crsRs(j1,j2,j3))*L2ssCoeffs(i1,i2,i3)+(-2.*cssRs(j1,j2,j3)-2.*crsRr(j1,j2,j3))*L1sCoeffs(i1,i2,i3)+(-1.*cssRs(j1,j2,j3)-1.*crsRr(j1,j2,j3))*L1ICoeff(i1,i2,i3)+(-1.*crrRr(j1,j2,j3)-1.*crsRs(j1,j2,j3))*L2sCoeff(i1,i2,i3)-1.*crrRs(j1,j2,j3)*L3ssCoeff(i1,i2,i3)-1.*cssRrs(j1,j2,j3)-1.*csRr(j1,j2,j3))*crsR(j1,j2,j3)+((-2.*L2ICoeff(i1,i2,i3)-2.*L2ssCoeffss(i1,i2,i3)-4.*L2sCoeffs(i1,i2,i3))*cssR(j1,j2,j3)+(-1.*crsRrr(j1,j2,j3)-2.*csRr(j1,j2,j3))*L1sCoeff(i1,i2,i3)-2.*crR(j1,j2,j3)*L3ssCoeff(i1,i2,i3)+(-2.*L2sCoeff(i1,i2,i3)-2.*L2ssCoeffs(i1,i2,i3))*csR(j1,j2,j3)+(-2.*crRr(j1,j2,j3)-1.*crrRrr(j1,j2,j3))*L2ssCoeff(i1,i2,i3)-2.*cssRr(j1,j2,j3)*L1ICoeff(i1,i2,i3)-4.*cssRr(j1,j2,j3)*L1sCoeffs(i1,i2,i3)-1.*cssRrr(j1,j2,j3)-2.*crrRr(j1,j2,j3)*L3ssCoeff(i1,i2,i3)-2.*crsRr(j1,j2,j3)*L2ssCoeffs(i1,i2,i3)-2.*crsRr(j1,j2,j3)*L2sCoeff(i1,i2,i3))*crrR(j1,j2,j3)+((-1.*crsRss(j1,j2,j3)-2.*crRs(j1,j2,j3))*L1sCoeff(i1,i2,i3)+(-4.*L1sCoeffs(i1,i2,i3)-2.*L1ICoeff(i1,i2,i3))*crR(j1,j2,j3)-2.*crsRs(j1,j2,j3)*L1ICoeff(i1,i2,i3)-2.*crrRs(j1,j2,j3)*L2ssCoeffs(i1,i2,i3)-2.*csRs(j1,j2,j3)-4.*crsRs(j1,j2,j3)*L1sCoeffs(i1,i2,i3)-1.*crrRss(j1,j2,j3)*L2ssCoeff(i1,i2,i3)-2.*crrRs(j1,j2,j3)*L2sCoeff(i1,i2,i3)-1.*cssRss(j1,j2,j3))*cssR(j1,j2,j3)+((-1.*crsRr(j1,j2,j3)-2.*csR(j1,j2,j3))*crR(j1,j2,j3)-1.*crsRs(j1,j2,j3)*csR(j1,j2,j3))*L1sCoeff(i1,i2,i3)-1.*pow(crR(j1,j2,j3),2.)*L2ssCoeff(i1,i2,i3)+(-1.*crrRr(j1,j2,j3)*L2ssCoeff(i1,i2,i3)-1.*cssRr(j1,j2,j3))*crR(j1,j2,j3)-1.*pow(csR(j1,j2,j3),2.)+(-1.*crrRs(j1,j2,j3)*L2ssCoeff(i1,i2,i3)-1.*cssRs(j1,j2,j3))*csR(j1,j2,j3)+((cssLrs(i1,i2,i3)+csLr(i1,i2,i3))*crsL(i1,i2,i3)+(cssLss(i1,i2,i3)+2.*csLs(i1,i2,i3))*cssL(i1,i2,i3)+cssLr(i1,i2,i3)*crL(i1,i2,i3)+csL(i1,i2,i3)*cssLs(i1,i2,i3)+crrL(i1,i2,i3)*cssLrr(i1,i2,i3)+pow(csL(i1,i2,i3),2.))*pow(DLR,2.))/pow(crrR(j1,j2,j3),2.);
  L4rCoeff(i1,i2,i3) = (-1.*pow(crsR(j1,j2,j3),2.)*L2rCoeffss(i1,i2,i3)+(-2.*L3rCoeffs(i1,i2,i3)*crrR(j1,j2,j3)+(-1.*crsRrs(j1,j2,j3)-1.*crRr(j1,j2,j3)-1.*csRs(j1,j2,j3))*L1rCoeffs(i1,i2,i3)-2.*cssR(j1,j2,j3)*L1rCoeffsss(i1,i2,i3)+(-1.*crrRrs(j1,j2,j3)-1.*crRs(j1,j2,j3))*L2rCoeff(i1,i2,i3)-2.*crR(j1,j2,j3)*L2rCoeffs(i1,i2,i3)+(-1.*cssRs(j1,j2,j3)-1.*crsRr(j1,j2,j3)-2.*csR(j1,j2,j3))*L1rCoeffss(i1,i2,i3)+(-1.*crrRr(j1,j2,j3)-1.*crsRs(j1,j2,j3))*L2rCoeffs(i1,i2,i3)-1.*crrRs(j1,j2,j3)*L3rCoeff(i1,i2,i3)-1.*crRrs(j1,j2,j3)*L1rCoeff(i1,i2,i3))*crsR(j1,j2,j3)+((-1.*crsRrr(j1,j2,j3)-2.*csRr(j1,j2,j3))*L1rCoeffs(i1,i2,i3)-2.*cssR(j1,j2,j3)*L2rCoeffss(i1,i2,i3)+(-2.*crRr(j1,j2,j3)-1.*crrRrr(j1,j2,j3))*L2rCoeff(i1,i2,i3)-2.*crR(j1,j2,j3)*L3rCoeff(i1,i2,i3)-2.*cssRr(j1,j2,j3)*L1rCoeffss(i1,i2,i3)+(-2.*csR(j1,j2,j3)-2.*crsRr(j1,j2,j3))*L2rCoeffs(i1,i2,i3)-2.*crrRr(j1,j2,j3)*L3rCoeff(i1,i2,i3)-1.*crRrr(j1,j2,j3)*L1rCoeff(i1,i2,i3))*crrR(j1,j2,j3)+((-1.*crsRss(j1,j2,j3)-2.*crRs(j1,j2,j3))*cssR(j1,j2,j3)+(-1.*crsRr(j1,j2,j3)-2.*csR(j1,j2,j3))*crR(j1,j2,j3)-1.*crsRs(j1,j2,j3)*csR(j1,j2,j3))*L1rCoeffs(i1,i2,i3)+(-2.*crR(j1,j2,j3)*L1rCoeffss(i1,i2,i3)-1.*crrRss(j1,j2,j3)*L2rCoeff(i1,i2,i3)-2.*crrRs(j1,j2,j3)*L2rCoeffs(i1,i2,i3)-2.*crsRs(j1,j2,j3)*L1rCoeffss(i1,i2,i3)-1.*crRss(j1,j2,j3)*L1rCoeff(i1,i2,i3))*cssR(j1,j2,j3)+(-1.*crrRr(j1,j2,j3)*crR(j1,j2,j3)-1.*csR(j1,j2,j3)*crrRs(j1,j2,j3)-1.*pow(crR(j1,j2,j3),2.))*L2rCoeff(i1,i2,i3)-1.*crRr(j1,j2,j3)*crR(j1,j2,j3)*L1rCoeff(i1,i2,i3)-1.*csR(j1,j2,j3)*crRs(j1,j2,j3)*L1rCoeff(i1,i2,i3)+pow(DLR,2.)*(crL(i1,i2,i3)*crLr(i1,i2,i3)+crLrr(i1,i2,i3)*crrL(i1,i2,i3)+crLrs(i1,i2,i3)*crsL(i1,i2,i3)+crLs(i1,i2,i3)*csL(i1,i2,i3)+crLss(i1,i2,i3)*cssL(i1,i2,i3)))/pow(crrR(j1,j2,j3),2.);
  L4sCoeff(i1,i2,i3) = ((-1.*L2sCoeffss(i1,i2,i3)-2.*L2ICoeffs(i1,i2,i3))*pow(crsR(j1,j2,j3),2.)+((-2.*L3sCoeffs(i1,i2,i3)-2.*L3ICoeff(i1,i2,i3))*crrR(j1,j2,j3)+(-6.*L1ICoeffss(i1,i2,i3)-2.*L1sCoeffsss(i1,i2,i3))*cssR(j1,j2,j3)+(-2.*L2ICoeff(i1,i2,i3)-2.*L2sCoeffs(i1,i2,i3))*crR(j1,j2,j3)+(-2.*L1sCoeffss(i1,i2,i3)-4.*L1ICoeffs(i1,i2,i3))*csR(j1,j2,j3)+(-1.*crsRrs(j1,j2,j3)-1.*crRr(j1,j2,j3)-1.*csRs(j1,j2,j3))*L1sCoeffs(i1,i2,i3)+(-1.*crsRrs(j1,j2,j3)-1.*crRr(j1,j2,j3)-1.*csRs(j1,j2,j3))*L1ICoeff(i1,i2,i3)+(-1.*crrRrs(j1,j2,j3)-1.*crRs(j1,j2,j3))*L2sCoeff(i1,i2,i3)+(-1.*cssRs(j1,j2,j3)-1.*crsRr(j1,j2,j3))*L1sCoeffss(i1,i2,i3)+(-2.*cssRs(j1,j2,j3)-2.*crsRr(j1,j2,j3))*L1ICoeffs(i1,i2,i3)+(-1.*crrRr(j1,j2,j3)-1.*crsRs(j1,j2,j3))*L2sCoeffs(i1,i2,i3)+(-1.*crrRr(j1,j2,j3)-1.*crsRs(j1,j2,j3))*L2ICoeff(i1,i2,i3)-1.*csRrs(j1,j2,j3)-1.*crRrs(j1,j2,j3)*L1sCoeff(i1,i2,i3)-1.*crrRs(j1,j2,j3)*L3sCoeff(i1,i2,i3))*crsR(j1,j2,j3)+((-2.*L2sCoeffss(i1,i2,i3)-4.*L2ICoeffs(i1,i2,i3))*cssR(j1,j2,j3)-2.*crR(j1,j2,j3)*L3sCoeff(i1,i2,i3)+(-2.*L2ICoeff(i1,i2,i3)-2.*L2sCoeffs(i1,i2,i3))*csR(j1,j2,j3)+(-1.*crsRrr(j1,j2,j3)-2.*csRr(j1,j2,j3))*L1sCoeffs(i1,i2,i3)+(-1.*crsRrr(j1,j2,j3)-2.*csRr(j1,j2,j3))*L1ICoeff(i1,i2,i3)+(-2.*crRr(j1,j2,j3)-1.*crrRrr(j1,j2,j3))*L2sCoeff(i1,i2,i3)-1.*csRrr(j1,j2,j3)-2.*cssRr(j1,j2,j3)*L1sCoeffss(i1,i2,i3)-4.*cssRr(j1,j2,j3)*L1ICoeffs(i1,i2,i3)-1.*crRrr(j1,j2,j3)*L1sCoeff(i1,i2,i3)-2.*crsRr(j1,j2,j3)*L2ICoeff(i1,i2,i3)-2.*crrRr(j1,j2,j3)*L3sCoeff(i1,i2,i3)-2.*crsRr(j1,j2,j3)*L2sCoeffs(i1,i2,i3))*crrR(j1,j2,j3)+((-2.*L1sCoeffss(i1,i2,i3)-4.*L1ICoeffs(i1,i2,i3))*crR(j1,j2,j3)+(-1.*crsRss(j1,j2,j3)-2.*crRs(j1,j2,j3))*L1sCoeffs(i1,i2,i3)+(-1.*crsRss(j1,j2,j3)-2.*crRs(j1,j2,j3))*L1ICoeff(i1,i2,i3)-1.*csRss(j1,j2,j3)-1.*crrRss(j1,j2,j3)*L2sCoeff(i1,i2,i3)-1.*crRss(j1,j2,j3)*L1sCoeff(i1,i2,i3)-2.*crrRs(j1,j2,j3)*L2sCoeffs(i1,i2,i3)-2.*crrRs(j1,j2,j3)*L2ICoeff(i1,i2,i3)-2.*crsRs(j1,j2,j3)*L1sCoeffss(i1,i2,i3)-4.*crsRs(j1,j2,j3)*L1ICoeffs(i1,i2,i3))*cssR(j1,j2,j3)-1.*pow(crR(j1,j2,j3),2.)*L2sCoeff(i1,i2,i3)+((-2.*L1sCoeffs(i1,i2,i3)-2.*L1ICoeff(i1,i2,i3))*csR(j1,j2,j3)-1.*csRr(j1,j2,j3)-1.*crsRr(j1,j2,j3)*L1ICoeff(i1,i2,i3)-1.*crsRr(j1,j2,j3)*L1sCoeffs(i1,i2,i3)-1.*crRr(j1,j2,j3)*L1sCoeff(i1,i2,i3)-1.*crrRr(j1,j2,j3)*L2sCoeff(i1,i2,i3))*crR(j1,j2,j3)+(-1.*csRs(j1,j2,j3)-1.*crrRs(j1,j2,j3)*L2sCoeff(i1,i2,i3)-1.*crRs(j1,j2,j3)*L1sCoeff(i1,i2,i3)-1.*crsRs(j1,j2,j3)*L1ICoeff(i1,i2,i3)-1.*crsRs(j1,j2,j3)*L1sCoeffs(i1,i2,i3))*csR(j1,j2,j3)+pow(DLR,2.)*(crL(i1,i2,i3)*csLr(i1,i2,i3)+crrL(i1,i2,i3)*csLrr(i1,i2,i3)+crsL(i1,i2,i3)*csLrs(i1,i2,i3)+csL(i1,i2,i3)*csLs(i1,i2,i3)+csLss(i1,i2,i3)*cssL(i1,i2,i3)))/pow(crrR(j1,j2,j3),2.);
  L4ICoeff(i1,i2,i3) = (-1.*pow(crsR(j1,j2,j3),2.)*L2ICoeffss(i1,i2,i3)+(-2.*L3ICoeffs(i1,i2,i3)*crrR(j1,j2,j3)+(-1.*crsRrs(j1,j2,j3)-1.*crRr(j1,j2,j3)-1.*csRs(j1,j2,j3))*L1ICoeffs(i1,i2,i3)-2.*cssR(j1,j2,j3)*L1ICoeffsss(i1,i2,i3)+(-1.*crrRrs(j1,j2,j3)-1.*crRs(j1,j2,j3))*L2ICoeff(i1,i2,i3)-2.*crR(j1,j2,j3)*L2ICoeffs(i1,i2,i3)+(-1.*cssRs(j1,j2,j3)-1.*crsRr(j1,j2,j3)-2.*csR(j1,j2,j3))*L1ICoeffss(i1,i2,i3)+(-1.*crrRr(j1,j2,j3)-1.*crsRs(j1,j2,j3))*L2ICoeffs(i1,i2,i3)-1.*crrRs(j1,j2,j3)*L3ICoeff(i1,i2,i3)-1.*crRrs(j1,j2,j3)*L1ICoeff(i1,i2,i3))*crsR(j1,j2,j3)+((-1.*crsRrr(j1,j2,j3)-2.*csRr(j1,j2,j3))*L1ICoeffs(i1,i2,i3)-2.*cssR(j1,j2,j3)*L2ICoeffss(i1,i2,i3)+(-2.*crRr(j1,j2,j3)-1.*crrRrr(j1,j2,j3))*L2ICoeff(i1,i2,i3)-2.*crR(j1,j2,j3)*L3ICoeff(i1,i2,i3)-2.*cssRr(j1,j2,j3)*L1ICoeffss(i1,i2,i3)+(-2.*csR(j1,j2,j3)-2.*crsRr(j1,j2,j3))*L2ICoeffs(i1,i2,i3)-2.*crrRr(j1,j2,j3)*L3ICoeff(i1,i2,i3)-1.*crRrr(j1,j2,j3)*L1ICoeff(i1,i2,i3))*crrR(j1,j2,j3)+((-1.*crsRss(j1,j2,j3)-2.*crRs(j1,j2,j3))*cssR(j1,j2,j3)+(-1.*crsRr(j1,j2,j3)-2.*csR(j1,j2,j3))*crR(j1,j2,j3)-1.*crsRs(j1,j2,j3)*csR(j1,j2,j3))*L1ICoeffs(i1,i2,i3)+(-1.*crrRss(j1,j2,j3)*L2ICoeff(i1,i2,i3)-2.*crsRs(j1,j2,j3)*L1ICoeffss(i1,i2,i3)-1.*crRss(j1,j2,j3)*L1ICoeff(i1,i2,i3)-2.*crrRs(j1,j2,j3)*L2ICoeffs(i1,i2,i3)-2.*crR(j1,j2,j3)*L1ICoeffss(i1,i2,i3))*cssR(j1,j2,j3)+(-1.*crrRr(j1,j2,j3)*crR(j1,j2,j3)-1.*csR(j1,j2,j3)*crrRs(j1,j2,j3)-1.*pow(crR(j1,j2,j3),2.))*L2ICoeff(i1,i2,i3)-1.*crRr(j1,j2,j3)*crR(j1,j2,j3)*L1ICoeff(i1,i2,i3)-1.*csR(j1,j2,j3)*crRs(j1,j2,j3)*L1ICoeff(i1,i2,i3))/pow(crrR(j1,j2,j3),2.);
} // end FOR_3D

// ----- Define coefficients in difference operators -----
Range R4(-2,2);
RealArray iCoeff(R4);
RealArray rCoeff(R4), rrCoeff(R4), rrrCoeff(R4), rrrrCoeff(R4);
RealArray sCoeff(R4), ssCoeff(R4), sssCoeff(R4), ssssCoeff(R4);
   iCoeff(-2)= 0.;    iCoeff(-1)=  0.;    iCoeff(0)=  1.;    iCoeff(1)= 0.;    iCoeff(2)= 0.;
   rCoeff(-2)= 1.;    rCoeff(-1)= -8.;    rCoeff(0)=  0.;    rCoeff(1)= 8.;    rCoeff(2)=-1.; rCoeff    /=(    12.*drL[0]);
  rrCoeff(-2)=-1.;   rrCoeff(-1)= 16.;   rrCoeff(0)=-30.;   rrCoeff(1)=16.;   rrCoeff(2)=-1.; rrCoeff   /=(12.*SQR(drL[0]));
 rrrCoeff(-2)=-1.;  rrrCoeff(-1)=  2.;  rrrCoeff(0)=  0.;  rrrCoeff(1)=-2.;  rrrCoeff(2)= 1.; rrrCoeff  /=( 2.*pow(drL[0],3));
rrrrCoeff(-2)= 1.; rrrrCoeff(-1)= -4.; rrrrCoeff(0)=  6.; rrrrCoeff(1)=-4.; rrrrCoeff(2)= 1.; rrrrCoeff /=(    pow(drL[0],4));
   sCoeff(-2)= 1.;    sCoeff(-1)= -8.;    sCoeff(0)=  0.;    sCoeff(1)= 8.;    sCoeff(2)=-1.; sCoeff    /=(    12.*drL[1]);
  ssCoeff(-2)=-1.;   ssCoeff(-1)= 16.;   ssCoeff(0)=-30.;   ssCoeff(1)=16.;   ssCoeff(2)=-1.; ssCoeff   /=(12.*SQR(drL[1]));
 sssCoeff(-2)=-1.;  sssCoeff(-1)=  2.;  sssCoeff(0)=  0.;  sssCoeff(1)=-2.;  sssCoeff(2)= 1.; sssCoeff  /=( 2.*pow(drL[1],3));
ssssCoeff(-2)= 1.; ssssCoeff(-1)= -4.; ssssCoeff(0)=  6.; ssssCoeff(1)=-4.; ssssCoeff(2)= 1.; ssssCoeff /=(    pow(drL[1],4));

// ----- Fill in the Matrix Coefficients for CHAMP -----
const Real dxs = (1-2*side2)*drR[axis2];
const Real h = dxs;
printF("WDH: grid=%d, (side,axis)=(%d,%d) (side2,axis2)=(%d,%d) dxs=%9.3e, Sl=%9.3e\n",grid, side,axis, side2,axis2, dxs,Sl);
RealArray coeff4(R4,R4);
const int e=0, c=0; // eqn number and component number
// ------- FILL CHAMP conditions into the matrix -----
// NOTE: skip top boundary if periodic (use indexRange) 
// NOTE: skip adjacent boundaries if Dirichlet BC *finish me**
extra=0; 
getBoundaryIndex(mg.indexRange(), side, axis, Ib1,Ib2,Ib3,extra);
getBoundaryIndex(mg2.indexRange(),side2,axis2,Jb1,Jb2,Jb3,extra);
int axisp = (axis+1) % numberOfDimensions;
FOR_3IJ(i1,i2,i3,Ib1,Ib2,Ib3,j1,j2,j3,Jb1,Jb2,Jb3)
{
  i1m=i1-is1, i2m=i2-is2, i3m=i3-is3; //  ghost point is (i1m,i2m,i3m)
  const Real bn = b1R(j1,j2,j3), bt = b2R(j1,j2,j3);  // b in normal and tangential directions
  // The next macro defines cI, cr, cs, crr, ...
  getChampCoeffOrder4NormalDirection0()
//  if( orderOfAccuracy==4 )
//  {
//    printF("domain2=%d: cI   =%12.4e, cr   =%12.4e, cs   =%12.4e, crr  =%12.4e, crs=%12.4e, css=%12.4e\n",domain2,cI,cr,cs,crr,crs,css);
//    printF("            crrr =%12.4e, crrs =%12.4e, crss =%12.4e, csss =%12.4e\n",crrr,crrs,crss,csss);
//    printF("            crrrr=%12.4e, crrrs=%12.4e, crrss=%12.4e, crsss=%12.4e, cssss=%12.4e\n",crrrr,crrrs,crrss,crsss,cssss);
//    printF("            L4ssssC=%12.4e, (Sl*h4By24-bn*h3By6)=%12.4e h4By24=%12.4e h3By6=%12.4e, Sl=%12.4e, bn=%12.4e, h=%12.4e, bt=%12.4e\n",L4ssssC,(Sl*h4By24-bn*h3By6),h4By24,h3By6,Sl,bn,h,bt);
//    if( domain2==0 )
//    {
//      OV_ABORT("stop here for now");
//    }
//  }
  ForStencil(m1,m2,m3)
  {
    int m  = M123(m1,m2,m3);        // the single-component coeff-index
    int mm = M123CE(m1,m2,m3,c,e);  // the system coeff-index
    coeff4(m1,m2) = 
             + cI   *   iCoeff(m1) * iCoeff(m2) 
             + cr   *   rCoeff(m1) * iCoeff(m2) 
             + cs   *   iCoeff(m1) * sCoeff(m2) 
             + crr  *  rrCoeff(m1) * iCoeff(m2) 
             + crs  *   rCoeff(m1) * sCoeff(m2) 
             + css  *   iCoeff(m1) *ssCoeff(m2) 
             + crrr   * rrrCoeff(m1)  *   iCoeff(m2) 
             + crrs   *  rrCoeff(m1)  *   sCoeff(m2) 
             + crss   *   rCoeff(m1)  *  ssCoeff(m2) 
             + csss   *   iCoeff(m1)  * sssCoeff(m2) 
             + crrrr  *rrrrCoeff(m1)  *   iCoeff(m2) 
             + crrrs  * rrrCoeff(m1)  *   sCoeff(m2) 
             + crrss  *  rrCoeff(m1)  *  ssCoeff(m2) 
             + crsss  *   rCoeff(m1)  * sssCoeff(m2) 
             + cssss  *   iCoeff(m1)  *ssssCoeff(m2) 
             ;
    if( fillMatrixWDH )
    {
      coeff(mm,i1m,i2m,i3m) = coeff4(m1,m2);
      // Specify that the above coeff value is the coefficient of component c at the grid point (j1,j2,j3).
      const int k1=i1+m1, k2=i2+m2, k3=i3+m3;   // the stencil is centred on the boundary pt (i1,i2,i3)    
      setEquationNumber(mm, e,i1m,i2m,i3m,  c,k1,k2,k3 );  // macro to set equationNumber                  
                                                                                                           
       // Fill Ghost 2 -- extrapolation for now                                                            
       const int ghost=2;                                                                                  
       fillGhostExtrapolation(ghost);                                                                      
    }
    if( twilightZoneFlow )
    {
      // --- For twilightZone we save some coefficients that go into the CHAMP matrix ---
      cc(i1,i2,i3, 0) = cI;   // coeff of I
      cc(i1,i2,i3, 1) = cr;   // coeff of ur
      cc(i1,i2,i3, 2) = cs;   // coeff of us
      cc(i1,i2,i3, 3) = crr;  // coeff of urr
      cc(i1,i2,i3, 4) = crs;  // coeff of urs
      cc(i1,i2,i3, 5) = css;  // coeff of uss
      cc(i1,i2,i3, 6) = crrr;    // coeff of urrr
      cc(i1,i2,i3, 7) = crrs;    // coeff of urrs
      cc(i1,i2,i3, 8) = crss;    // coeff of urss
      cc(i1,i2,i3, 9) = csss;    // coeff of usss
      cc(i1,i2,i3,10) = crrrr;   // coeff of urrrr
      cc(i1,i2,i3,11) = crrrs;   // coeff of urrrs
      cc(i1,i2,i3,12) = crrss;   // coeff of urrss
      cc(i1,i2,i3,13) = crsss;   // coeff of ursss
      cc(i1,i2,i3,14) = cssss;   // coeff of ussss
    } 
  } // end ForStencil
 if( debug & 4 ) 
 {
   printF("WDH:order=4: (i1,i2,i3)=(%3d,%3d,%3d)\n",i1,i2,i3);     
   printF("    coeff4=(%10.3e,%10.3e,%10.3e,%10.3e,%10.3e,\n"    
          "            %10.3e,%10.3e,%10.3e,%10.3e,%10.3e,\n"    
          "            %10.3e,%10.3e,%10.3e,%10.3e,%10.3e,\n"    
          "            %10.3e,%10.3e,%10.3e,%10.3e,%10.3e,\n"    
          "            %10.3e,%10.3e,%10.3e,%10.3e,%10.3e )\n",  
           coeff4(-2,-2),coeff4(-1,-2),coeff4(0,-2),coeff4(1,-2),coeff4(2,-2),  
           coeff4(-2,-1),coeff4(-1,-1),coeff4(0,-1),coeff4(1,-1),coeff4(2,-1),  
           coeff4(-2, 0),coeff4(-1, 0),coeff4(0, 0),coeff4(1, 0),coeff4(2, 0),  
           coeff4(-2, 1),coeff4(-1, 1),coeff4(0, 1),coeff4(1, 1),coeff4(2, 1),  
           coeff4(-2, 2),coeff4(-1, 2),coeff4(0, 2),coeff4(1, 2),coeff4(2, 2)); 
 }
} // end FOR_3D
