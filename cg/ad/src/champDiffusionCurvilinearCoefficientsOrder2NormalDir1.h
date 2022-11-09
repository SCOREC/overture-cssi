// Coefficients for CHAMP curvilinear order 2, NORMAL DIRECTION=1.
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
#defineMacro diffa2(u,i1,i2,i3,dr) diffr2(u,i1,i2,i3,dr)
#defineMacro diffaa2(u,i1,i2,i3,dr) diffrr2(u,i1,i2,i3,dr)
#defineMacro diffa(u,i1,i2,i3,dr) diffr2(u,i1,i2,i3,dr)
#defineMacro diffaa(u,i1,i2,i3,dr) diffrr2(u,i1,i2,i3,dr)
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

// Fill 3rd and 4th ghost in tangential directions for b1, b2 
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
  const Real dnR = fabs( drR[axis2]/b2R(j1,j2,j3) );    // approx. normal grid spacing on right
  Sl = pl/dnR; 
  champParameters(4,side,axis,grid)=Sl; // save Sl
}
extra=0;
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
extra=0;
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
        // Evaluate derivatives of metrics to order 2 .. we are computing more than we need here *fix me*
        rxxL(m1,m2) = diffMx2(rxL,i1,i2,i3,m1,m2,drL);
        rxyL(m1,m2) = diffMy2(rxL,i1,i2,i3,m1,m2,drL);
                                                      
        rxxR(m1,m2) = diffMx2(rxR,j1,j2,j3,m1,m2,drR);
        rxyR(m1,m2) = diffMy2(rxR,j1,j2,j3,m1,m2,drR);
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


// ----- Evaluate coefficients of L1 on the boundary -----

// 2*numGhost = 2 extra ghost needed for 2nd order, true ? 
extra=2*numGhost;
getBoundaryIndex(mg.indexRange(), side ,axis ,Ib1,Ib2,Ib3,extra);
getBoundaryIndex(mg2.indexRange(),side2,axis2,Jb1,Jb2,Jb3,extra);
RealArray L1ICoeff(Ib1,Ib2,Ib3), L1ICoeffr(Ib1,Ib2,Ib3);
RealArray L1rCoeff(Ib1,Ib2,Ib3), L1rCoeffr(Ib1,Ib2,Ib3);
RealArray L1sCoeff(Ib1,Ib2,Ib3), L1sCoeffr(Ib1,Ib2,Ib3);
FOR_3IJ(i1,i2,i3,Ib1,Ib2,Ib3,j1,j2,j3,Jb1,Jb2,Jb3)
{
  const Real n1R = -n1L(i1,i2,i3), n2R = -n2L(i1,i2,i3);    // for advection, use right normal 
  getAdvectionCoefficients(i1,i2,i3,u1DL,u2DL, j1,j2,j3,u1DR,u2DR );
  L1ICoeff(i1,i2,i3) = ((-1.*KLR*u2DL+1.*u2DR)*n2R+(-1.*KLR*u1DL+1.*u1DR)*n1R)/b2R(j1,j2,j3);
  L1rCoeff(i1,i2,i3) = (1.*KLR*b1L(i1,i2,i3)-1.*b1R(j1,j2,j3))/b2R(j1,j2,j3);
  L1sCoeff(i1,i2,i3) = 1.*b2L(i1,i2,i3)*KLR/b2R(j1,j2,j3);
} // end FOR_3D

// ----- Evaluate coefficients of L2 on the boundary -----
extra=0;
getBoundaryIndex(mg.indexRange(), side,  axis,Ib1,Ib2,Ib3,extra);
getBoundaryIndex(mg2.indexRange(),side2,axis2,Jb1,Jb2,Jb3,extra);
RealArray L2rrCoeff(Ib1,Ib2,Ib3);
RealArray L2rsCoeff(Ib1,Ib2,Ib3);
RealArray L2ssCoeff(Ib1,Ib2,Ib3);
RealArray L2rCoeff(Ib1,Ib2,Ib3);
RealArray L2sCoeff(Ib1,Ib2,Ib3);
RealArray L2ICoeff(Ib1,Ib2,Ib3);
FOR_3IJ(i1,i2,i3,Ib1,Ib2,Ib3,j1,j2,j3,Jb1,Jb2,Jb3)
{
  // Compute tangential derivatives of L1 - label tangential directions a and b (for 3D)
  L1rCoeffr(i1,i2,i3) = diffa(L1rCoeff,i1,i2,i3,drL);
  L1sCoeffr(i1,i2,i3) = diffa(L1sCoeff,i1,i2,i3,drL);
  L1ICoeffr(i1,i2,i3) = diffa(L1ICoeff,i1,i2,i3,drL);
  // Compute L2 coefficients
  L2rrCoeff(i1,i2,i3) = 1.*(DLR*crrL(i1,i2,i3)-crsR(j1,j2,j3)*L1rCoeff(i1,i2,i3)-crrR(j1,j2,j3))/cssR(j1,j2,j3);
  L2rsCoeff(i1,i2,i3) = (1.*DLR*crsL(i1,i2,i3)-1.*crsR(j1,j2,j3)*L1sCoeff(i1,i2,i3))/cssR(j1,j2,j3);
  L2ssCoeff(i1,i2,i3) = 1.*cssL(i1,i2,i3)*DLR/cssR(j1,j2,j3);
  L2rCoeff(i1,i2,i3) = ((-1.*L1ICoeff(i1,i2,i3)-1.*L1rCoeffr(i1,i2,i3))*crsR(j1,j2,j3)+1.*DLR*crL(i1,i2,i3)-1.*csR(j1,j2,j3)*L1rCoeff(i1,i2,i3)-1.*crR(j1,j2,j3))/cssR(j1,j2,j3);
  L2sCoeff(i1,i2,i3) = 1.*(DLR*csL(i1,i2,i3)-crsR(j1,j2,j3)*L1sCoeffr(i1,i2,i3)-csR(j1,j2,j3)*L1sCoeff(i1,i2,i3))/cssR(j1,j2,j3);
  L2ICoeff(i1,i2,i3) = (-1.*crsR(j1,j2,j3)*L1ICoeffr(i1,i2,i3)-1.*csR(j1,j2,j3)*L1ICoeff(i1,i2,i3))/cssR(j1,j2,j3);
} // end FOR_3D

// ----- Define coefficients in difference operators -----
Range R4(-1,1);
RealArray iCoeff(R4);
RealArray rCoeff(R4), rrCoeff(R4), rrrCoeff(R4), rrrrCoeff(R4);
RealArray sCoeff(R4), ssCoeff(R4), sssCoeff(R4), ssssCoeff(R4);
   iCoeff(-1)=  0.;    iCoeff(0)=  1.;    iCoeff(1)= 0.;   
   rCoeff(-1)= -1.;    rCoeff(0)=  0.;    rCoeff(1)= 1.;    rCoeff    /=(2.*drL[0])  ;
  rrCoeff(-1)=  1.;   rrCoeff(0)= -2.;   rrCoeff(1)= 1.;    rrCoeff   /=(SQR(drL[0]));
   sCoeff(-1)= -1.;    sCoeff(0)=  0.;    sCoeff(1)= 1.;    sCoeff    /=(2.*drL[1])  ;
  ssCoeff(-1)=  1.;   ssCoeff(0)= -2.;   ssCoeff(1)= 1.;    ssCoeff   /=(SQR(drL[1]));

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
  const Real bn = b2R(j1,j2,j3), bt = b1R(j1,j2,j3);  // b in normal and tangential directions
  // The next macro defines cI, cr, cs, crr, ...
  getChampCoeffOrder2NormalDirection1()
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
    } 
  } // end ForStencil
 if( debug & 4 ) 
 {
   printF("WDH:order=2: (i1,i2,i3)=(%3d,%3d,%3d)\n",i1,i2,i3);     
   printF("    coeff =(%10.3e,%10.3e,%10.3e,\n"    
          "            %10.3e,%10.3e,%10.3e,\n"    
          "            %10.3e,%10.3e,%10.3e )\n",  
           coeff4(-1,-1),coeff4(0,-1),coeff4(1,-1),  
           coeff4(-1, 0),coeff4(0, 0),coeff4(1, 0),  
           coeff4(-1, 1),coeff4(0, 1),coeff4(1, 1)); 
 }
} // end FOR_3D
