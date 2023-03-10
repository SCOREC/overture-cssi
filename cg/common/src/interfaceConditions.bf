c *******************************************************************************
c   Interface boundary conditions
c *******************************************************************************

c These next include files will define the macros that will define the difference approximations
c The actual macro is called below
#Include "defineDiffNewerOrder2f.h"
#Include "defineDiffNewerOrder4f.h"


#beginMacro beginLoops(n1a,n1b,n2a,n2b,n3a,n3b,na,nb)
do i3=n3a,n3b
do i2=n2a,n2b
do i1=n1a,n1b
do n=na,nb
  ! write(*,'(" periodic i1,i2,i3,n=",4i4)') i1,i2,i3,n
#endMacro

#beginMacro endLoops()
end do
end do
end do
end do
#endMacro


#beginMacro beginLoops2d()
 i3=n3a
 j3=m3a

 j2=m2a
 do i2=n2a,n2b
  j1=m1a
  do i1=n1a,n1b
#endMacro
#beginMacro endLoops2d()
   j1=j1+1
  end do
  j2=j2+1
 end do
#endMacro

#beginMacro beginGhostLoops2d()
 i3=n3a
 j3=m3a
 j2=mm2a
 do i2=nn2a,nn2b
  j1=mm1a
  do i1=nn1a,nn1b
#endMacro

#beginMacro beginLoops3d()
 j3=m3a
 do i3=n3a,n3b
 j2=m2a
 do i2=n2a,n2b
  j1=m1a
  do i1=n1a,n1b
#endMacro
#beginMacro endLoops3d()
   j1=j1+1
  end do
  j2=j2+1
 end do
  j3=j3+1
 end do
#endMacro

#beginMacro beginGhostLoops3d()
 j3=mm3a
 do i3=nn3a,nn3b
 j2=mm2a
 do i2=nn2a,nn2b
  j1=mm1a
  do i1=nn1a,nn1b
#endMacro

#defineMacro extrap2(uu,k1,k2,k3,kc,ks1,ks2,ks3) \
            (2.*uu(k1,k2,k3,kc)-uu(k1+ks1,k2+ks2,k3+ks3,kc))

#defineMacro extrap3(uu,k1,k2,k3,kc,ks1,ks2,ks3) \
            (3.*uu(k1,k2,k3,kc)-3.*uu(k1+ks1,k2+ks2,k3+ks3,kc)\
            +   uu(k1+2*ks1,k2+2*ks2,k3+2*ks3,kc))

#defineMacro extrap4(uu,k1,k2,k3,kc,ks1,ks2,ks3) \
            (4.*uu(k1,k2,k3,kc)-6.*uu(k1+ks1,k2+ks2,k3+ks3,kc)\
            +4.*uu(k1+2*ks1,k2+2*ks2,k3+2*ks3,kc)-uu(k1+3*ks1,k2+3*ks2,k3+3*ks3,kc))

#defineMacro extrap5(uu,k1,k2,k3,kc,ks1,ks2,ks3) \
            (5.*uu(k1,k2,k3,kc)-10.*uu(k1+ks1,k2+ks2,k3+ks3,kc)\
            +10.*uu(k1+2*ks1,k2+2*ks2,k3+2*ks3,kc)-5.*uu(k1+3*ks1,k2+3*ks2,k3+3*ks3,kc)\
            +uu(k1+4*ks1,k2+4*ks2,k3+4*ks3,kc))

c This macro will assign the jump conditions on the boundary
c DIM (input): number of dimensions (2 or 3)
c GRIDTYPE (input) : curvilinear or rectangular
#beginMacro boundaryJumpConditions(DIM,GRIDTYPE)
 #If #DIM eq "2"
  if( eps1.lt.eps2 )then
    epsRatio=eps1/eps2
    beginGhostLoops2d()
      ! eps2 n.u2 = eps1 n.u1
      !     tau.u2 = tau.u1

      #If #GRIDTYPE eq "curvilinear"
       an1=rsxy1(i1,i2,i3,axis1,0)   ! normal (an1,an2)
       an2=rsxy1(i1,i2,i3,axis1,1)
       aNorm=max(epsx,sqrt(an1**2+an2**2))
       an1=an1/aNorm
       an2=an2/aNorm
      #Elif #GRIDTYPE eq "rectangular"
       an1=an1Cartesian
       an2=an2Cartesian
      #Else
         stop 1111
      #End
      ua=u1(i1,i2,i3,ex)
      ub=u1(i1,i2,i3,ey)
      nDotU = an1*ua+an2*ub
      ! u2 equals u1 but with normal component = eps1/eps2*(n.u1)
      u2(j1,j2,j3,ex) = ua + (nDotU*epsRatio - nDotU)*an1
      u2(j1,j2,j3,ey) = ub + (nDotU*epsRatio - nDotU)*an2
      u2(j1,j2,j3,hz) = u1(i1,i2,i3,hz)
    endLoops2d()
  else
    epsRatio=eps2/eps1
    beginGhostLoops2d()
      ! eps2 n.u2 = eps1 n.u1
      !     tau.u2 = tau.u1

      #If #GRIDTYPE eq "curvilinear"
       an1=rsxy1(i1,i2,i3,axis1,0)
       an2=rsxy1(i1,i2,i3,axis1,1)
       aNorm=max(epsx,sqrt(an1**2+an2**2))
       an1=an1/aNorm
       an2=an2/aNorm
      #Elif #GRIDTYPE eq "rectangular"
       an1=an1Cartesian
       an2=an2Cartesian
      #Else
        stop 1112
      #End
      ua=u2(j1,j2,j3,ex)
      ub=u2(j1,j2,j3,ey)

      nDotU = an1*ua+an2*ub

      u1(i1,i2,i3,ex) = ua + (nDotU*epsRatio - nDotU)*an1
      u1(i1,i2,i3,ey) = ub + (nDotU*epsRatio - nDotU)*an2
      u1(i1,i2,i3,hz) = u2(j1,j2,j3,hz)
    endLoops2d()
  end if
 #Else
   stop 7742
 #End
#endMacro

c ** Precompute the derivatives of rsxy ***
c assign rvx(m) = (rx,sy)
c        rvxx(m) = (rxx,sxx)
#beginMacro computeRxDerivatives(rv,rsxy,i1,i2,i3)
do m=0,nd-1
 rv ## x(m)   =rsxy(i1,i2,i3,m,0)
 rv ## y(m)   =rsxy(i1,i2,i3,m,1)

 rv ## xx(m)  =rsxy ## x22(i1,i2,i3,m,0)
 rv ## xy(m)  =rsxy ## x22(i1,i2,i3,m,1)
 rv ## yy(m)  =rsxy ## y22(i1,i2,i3,m,1)

 rv ## xxx(m) =rsxy ## xx22(i1,i2,i3,m,0)
 rv ## xxy(m) =rsxy ## xx22(i1,i2,i3,m,1)
 rv ## xyy(m) =rsxy ## xy22(i1,i2,i3,m,1)
 rv ## yyy(m) =rsxy ## yy22(i1,i2,i3,m,1)

 rv ## xxxx(m)=rsxy ## xxx22(i1,i2,i3,m,0)
 rv ## xxyy(m)=rsxy ## xyy22(i1,i2,i3,m,0)
 rv ## yyyy(m)=rsxy ## yyy22(i1,i2,i3,m,1)
end do
#endMacro

c assign some temporary variables that are used in the evaluation of the operators
#beginMacro setJacobian(rv,axis1,axisp1)
 rx   =rv ## x(axis1)   
 ry   =rv ## y(axis1)   
                    
 rxx  =rv ## xx(axis1)  
 rxy  =rv ## xy(axis1)  
 ryy  =rv ## yy(axis1)  
                    
 rxxx =rv ## xxx(axis1) 
 rxxy =rv ## xxy(axis1) 
 rxyy =rv ## xyy(axis1) 
 ryyy =rv ## yyy(axis1) 
                    
 rxxxx=rv ## xxxx(axis1)
 rxxyy=rv ## xxyy(axis1)
 ryyyy=rv ## yyyy(axis1)

 sx   =rv ## x(axis1p1)   
 sy   =rv ## y(axis1p1)   
                    
 sxx  =rv ## xx(axis1p1)  
 sxy  =rv ## xy(axis1p1)  
 syy  =rv ## yy(axis1p1)  
                    
 sxxx =rv ## xxx(axis1p1) 
 sxxy =rv ## xxy(axis1p1) 
 sxyy =rv ## xyy(axis1p1) 
 syyy =rv ## yyy(axis1p1) 
                    
 sxxxx=rv ## xxxx(axis1p1)
 sxxyy=rv ## xxyy(axis1p1)
 syyyy=rv ## yyyy(axis1p1)

#endMacro


! update the periodic ghost points
#beginMacro periodicUpdate2d(u,bc,gid,side,axis)

axisp1=mod(axis+1,nd)
if( bc(0,axisp1).lt.0 )then
  ! direction axisp1 is periodic
  diff(axis)=0
  diff(axisp1)=gid(1,axisp1)-gid(0,axisp1)

  if( side.eq.0 )then
    ! assign 4 ghost points outside lower corner
    np1a=gid(0,0)-2
    np1b=gid(0,0)-1
    np2a=gid(0,1)-2
    np2b=gid(0,1)-1

    beginLoops(np1a,np1b,np2a,np2b,n3a,n3b,ex,hz)
     u(i1,i2,i3,n) = u(i1+diff(0),i2+diff(1),i3,n)
    endLoops()

    ! assign 4 ghost points outside upper corner
    if( axis.eq.0 )then
      np2a=gid(1,axisp1)+1
      np2b=gid(1,axisp1)+2
    else
      np1a=gid(1,axisp1)+1
      np1b=gid(1,axisp1)+2
    end if

    beginLoops(np1a,np1b,np2a,np2b,n3a,n3b,ex,hz)
     u(i1,i2,i3,n) = u(i1-diff(0),i2-diff(1),i3,n)
    endLoops()

  else

    ! assign 4 ghost points outside upper corner
    np1a=gid(1,0)+1
    np1b=gid(1,0)+2
    np2a=gid(1,1)+1
    np2b=gid(1,1)+2

    beginLoops(np1a,np1b,np2a,np2b,n3a,n3b,ex,hz)
     u(i1,i2,i3,n) = u(i1-diff(0),i2-diff(1),i3,n)
    endLoops()

    if( axis.eq.0 )then
      np2a=gid(0,axisp1)-2
      np2b=gid(0,axisp1)-1
    else
      np1a=gid(0,axisp1)-2
      np1b=gid(0,axisp1)-1
    end if

    beginLoops(np1a,np1b,np2a,np2b,n3a,n3b,ex,hz)
     u(i1,i2,i3,n) = u(i1+diff(0),i2+diff(1),i3,n)
    endLoops()
  end if

endif


#endMacro


#beginMacro getExact(ep,xy, i1,i2,i3,m, ue )
 ue=ogf(ep,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,kd3),m,t)
#endMacro

#beginMacro getDerivs(ep,xy, i1,i2,i3,m, uex,uey,uez, uexx,ueyy,uezz )
 call ogderiv(ep, 0,1,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,kd3),t,m,uex)
 call ogderiv(ep, 0,0,1,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,kd3),t,m,uey)
 call ogderiv(ep, 0,2,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,kd3),t,m,uexx)
 call ogderiv(ep, 0,0,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,kd3),t,m,ueyy)
if( nd.gt.2 )then
  call ogderiv(ep, 0,0,0,1, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,kd3),t,m,uez)
  call ogderiv(ep, 0,0,0,2, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,kd3),t,m,uezz)
end if
#endMacro

      subroutine interfaceCgCm( nd, nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,\
                               gridIndexRange1, u1, mask1,rsxy1, xy1, boundaryCondition1, \
                               md1a,md1b,md2a,md2b,md3a,md3b,\
                               gridIndexRange2, u2, mask2,rsxy2, xy2, boundaryCondition2, \
                               ipar, rpar, \
                               aa2,aa4,aa8, ipvt2,ipvt4,ipvt8, \
                               ierr )
c ===================================================================================
c  Interface boundary conditions for Cg
c
c  gridType : 0=rectangular, 1=curvilinear
c
c  u1: solution on the "left" of the interface
c  u2: solution on the "right" of the interface
c
c  aa2,aa4,aa8 : real work space arrays that must be saved from call to call
c  ipvt2,ipvt4,ipvt8: integer work space arrays that must be saved from call to call
c ===================================================================================

      implicit none

      integer nd, \
              nd1a,nd1b,nd2a,nd2b,nd3a,nd3b, \
              md1a,md1b,md2a,md2b,md3a,md3b, \
              n1a,n1b,n2a,n2b,n3a,n3b,  \
              m1a,m1b,m2a,m2b,m3a,m3b,  \
              ierr

      real u1(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
      integer mask1(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
      real rsxy1(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1,0:nd-1)
      real xy1(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1)
      integer gridIndexRange1(0:1,0:2),boundaryCondition1(0:1,0:2)

      real u2(md1a:md1b,md2a:md2b,md3a:md3b,0:*)
      integer mask2(md1a:md1b,md2a:md2b,md3a:md3b)
      real rsxy2(md1a:md1b,md2a:md2b,md3a:md3b,0:nd-1,0:nd-1)
      real xy2(md1a:md1b,md2a:md2b,md3a:md3b,0:nd-1)
      integer gridIndexRange2(0:1,0:2),boundaryCondition2(0:1,0:2)

      integer ipar(0:*)
      real rpar(0:*)

      ! work space arrays that must be saved from call to call:
      real aa2(0:1,0:1,0:1,0:*),aa4(0:3,0:3,0:1,0:*),aa8(0:7,0:7,0:1,0:*)
      integer ipvt2(0:1,0:*), ipvt4(0:3,0:*), ipvt8(0:7,0:*)

      real ogf
c     --- local variables ----
      
      integer side1,axis1,grid1,side2,axis2,grid2,gridType,orderOfAccuracy,orderOfExtrapolation,useForcing,\
        tc1,tc2,useWhereMask,debug,axis1p1,axis2p1,nn,n1,n2,np,myid,iofile,normalSign1,normalSign2
      real dx1(0:2),dr1(0:2),dx2(0:2),dr2(0:2)
      real dx(0:2),dr(0:2)
      real t,ep1,ep2,dt,ktc1,ktc2,kappa1,kappa2
      integer axisp1,axisp2,i1,i2,i3,is1,is2,is3,j1,j2,j3,js1,js2,js3,ks1,ks2,ks3,is,js,it,nit
      integer option,initialized
      integer id1(0:2),id2(0:2),id3(0:2), jd1(0:2),jd2(0:2),jd3(0:2)
      integer i1p,i2p,i3p, j1p,j2p,j3p, side, ia1,ia2,ia3, ja1,ja2,ja3

      real u1e,u1ex,u1ey,u1ez, u1exx,u1eyy,u1ezz
      real u2e,u2ex,u2ey,u2ez, u2exx,u2eyy,u2ezz

      integer numGhost
      integer nn1a,nn1b,nn2a,nn2b,nn3a,nn3b
      integer mm1a,mm1b,mm2a,mm2b,mm3a,mm3b

      real rx1,ry1,rx2,ry2

      real aLap0,aLap1,bLap0,bLap1,aLapX0,aLapX1,bLapY0,bLapY1,cLapX0,cLapX1,dLapY0,dLapY1,aLapSq0,aLapSq1,bLapSq0,bLapSq1
      real a11,a12,a21,a22,det,b0,b1,b2

      real a0,a1,cc0,cc1,d0,d1,dr0,ds0
      real aNormSq,divu,uAve,ktcAve

      real epsRatio,an1,an2,an3,aNorm,ua,ub,nDotU
      real epsx

      real tau1,tau2,tau3,clap1,clap2,ulap1,vlap1,wlap1,ulap2,vlap2,wlap2,an1Cartesian,an2Cartesian
      real ulapSq1,vlapSq1,ulapSq2,vlapSq2,wlapSq1,wlapSq2

      integer np1a,np1b,np2a,np2b,np3a,np3b,diff(0:2)

      real rx,ry,rxx,rxy,ryy,rxxx,rxxy,rxyy,ryyy,rxxxx,rxxyy,ryyyy
      real sx,sy,sxx,sxy,syy,sxxx,sxxy,sxyy,syyy,sxxxx,sxxyy,syyyy

      real rv1x(0:2),rv1y(0:2),rv1xx(0:2),rv1xy(0:2),rv1yy(0:2),rv1xxx(0:2),rv1xxy(0:2),rv1xyy(0:2),rv1yyy(0:2),\
           rv1xxxx(0:2),rv1xxyy(0:2),rv1yyyy(0:2)
      real sv1x(0:2),sv1y(0:2),sv1xx(0:2),sv1xy(0:2),sv1yy(0:2),sv1xxx(0:2),sv1xxy(0:2),sv1xyy(0:2),sv1yyy(0:2),\
           sv1xxxx(0:2),sv1xxyy(0:2),sv1yyyy(0:2)
      real rv2x(0:2),rv2y(0:2),rv2xx(0:2),rv2xy(0:2),rv2yy(0:2),rv2xxx(0:2),rv2xxy(0:2),rv2xyy(0:2),rv2yyy(0:2),\
           rv2xxxx(0:2),rv2xxyy(0:2),rv2yyyy(0:2)
      real sv2x(0:2),sv2y(0:2),sv2xx(0:2),sv2xy(0:2),sv2yy(0:2),sv2xxx(0:2),sv2xxy(0:2),sv2xyy(0:2),sv2yyy(0:2),\
           sv2xxxx(0:2),sv2xxyy(0:2),sv2yyyy(0:2)

      integer numberOfEquations,job
      real a2(0:1,0:1),a4(0:3,0:3),a8(0:7,0:7),q(0:11),f(0:11),rcond,work(0:11)
      integer ipvt(0:11)
      real aa(0:1,0:1), ff(0:11), rr(0:11)
      real err

c**   ! boundary conditions parameters
c**   #Include "bcDefineFortranInclude.h"
 
      integer rectangular,curvilinear
      parameter(\
        rectangular=0,\
        curvilinear=1)


      integer kd,m,n,kd3
      logical extrapTangential
c     real rx,ry,rz,sx,sy,sz,tx,ty,tz
      declareDifferenceNewOrder2(u1,rsxy1,dr1,dx1,RX)
      declareDifferenceNewOrder2(u2,rsxy2,dr2,dx2,RX)

      declareDifferenceNewOrder4(u1,rsxy1,dr1,dx1,RX)
      declareDifferenceNewOrder4(u2,rsxy2,dr2,dx2,RX)

c     --- start statement function ----
c**      extrapTangential(bc) = bc.lt.0
      extrapTangential(m) = .true.

c.......statement functions for jacobian
c     rx(i1,i2,i3)=rsxy1(i1,i2,i3,0,0)
c     ry(i1,i2,i3)=rsxy1(i1,i2,i3,0,1)
c     rz(i1,i2,i3)=rsxy1(i1,i2,i3,0,2)
c     sx(i1,i2,i3)=rsxy1(i1,i2,i3,1,0)
c     sy(i1,i2,i3)=rsxy1(i1,i2,i3,1,1)
c     sz(i1,i2,i3)=rsxy1(i1,i2,i3,1,2)
c     tx(i1,i2,i3)=rsxy1(i1,i2,i3,2,0)
c     ty(i1,i2,i3)=rsxy1(i1,i2,i3,2,1)
c     tz(i1,i2,i3)=rsxy1(i1,i2,i3,2,2) 


c     The next macro call will define the difference approximation statement functions
      defineDifferenceNewOrder2Components1(u1,rsxy1,dr1,dx1,RX)
      defineDifferenceNewOrder2Components1(u2,rsxy2,dr2,dx2,RX)

      defineDifferenceNewOrder4Components1(u1,rsxy1,dr1,dx1,RX)
      defineDifferenceNewOrder4Components1(u2,rsxy2,dr2,dx2,RX)

c............... end statement functions

      ierr=0

      side1                =ipar(0)
      axis1                =ipar(1)
      grid1                =ipar(2)
      n1a                  =ipar(3)
      n1b                  =ipar(4)
      n2a                  =ipar(5)
      n2b                  =ipar(6)
      n3a                  =ipar(7)
      n3b                  =ipar(8)

      side2                =ipar(9)
      axis2                =ipar(10)
      grid2                =ipar(11)
      m1a                  =ipar(12)
      m1b                  =ipar(13)
      m2a                  =ipar(14)
      m2b                  =ipar(15)
      m3a                  =ipar(16)
      m3b                  =ipar(17)

      gridType             =ipar(18)
      orderOfAccuracy      =ipar(19)
      orderOfExtrapolation =ipar(20)
      useForcing           =ipar(21)
      tc1                  =ipar(22)  ! T
      tc2                  =ipar(23) 
      np                   =ipar(24)
      myid                 =ipar(25)
      normalSign1          =ipar(26)  ! for parallel we may flip the sign of the normal
      normalSign2          =ipar(27)  ! for parallel we may flip the sign of the normal
c      solveForE            =ipar(28)
c      solveForH            =ipar(29)
      useWhereMask         =ipar(30)
      debug                =ipar(31)
      nit                  =ipar(32)
      option               =ipar(33)
      initialized          =ipar(34)
     
      dx1(0)                =rpar(0)
      dx1(1)                =rpar(1)
      dx1(2)                =rpar(2)
      dr1(0)                =rpar(3)
      dr1(1)                =rpar(4)
      dr1(2)                =rpar(5)

      dx2(0)                =rpar(6)
      dx2(1)                =rpar(7)
      dx2(2)                =rpar(8)
      dr2(0)                =rpar(9)
      dr2(1)                =rpar(10)
      dr2(2)                =rpar(11)

      t                    =rpar(12)
      ep1                  =rpar(13) ! pointer for exact solution
      ep2                  =rpar(14) ! pointer for exact solution
      dt                   =rpar(15)
      ktc1                 =rpar(16)
      ktc2                 =rpar(17)
      kappa1               =rpar(18)
      kappa2               =rpar(19)
     
      kd3=min(nd-1,2)  ! for indexing xy
      
      ! kkc 080516 
      ktcAve = ktc1+ktc2

      ! iofile = file for debug output
      if( np.eq.1 )then
        iofile=6
      else
        iofile=10+myid
      end if

      if( debug.gt.3 )then
        write(iofile,'(" interfaceCgCm: ktc1,ktc2=",2f10.5," kappa1,kappa2=",2e10.2," gridType,tc1,tc2=",3i2)') ktc1,ktc2,kappa1,kappa2,gridType,tc1,tc2
           ! '
      end if

      if( nit.lt.0 .or. nit.gt.100 )then
        write(*,'(" interfaceBC: ERROR: nit=",i9)') nit
        nit=max(1,min(100,nit))
      end if

      if( debug.gt.3 )then
        write(iofile,'(" interfaceCgCm: **START** grid1=",i4," side1,axis1=",2i2)') grid1,side1,axis1
           ! '
        write(iofile,'(" interfaceCgCm: **START** grid2=",i4," side2,axis2=",2i2)') grid2,side2,axis2
           ! '
        write(iofile,'("n1a,n1b,...=",6i5)') n1a,n1b,n2a,n2b,n3a,n3b
        write(iofile,'("m1a,m1b,...=",6i5)') m1a,m1b,m2a,m2b,m3a,m3b
        write(iofile,'("dr1,dr2=",6e9.2)') dr1(0),dr1(1),dr1(2),dr2(0),dr2(1),dr2(2)

      end if
      if( debug.gt.7 )then
      write(iofile,*) 'u1=',((((u1(i1,i2,i3,m),m=0,2),i1=n1a,n1b),i2=n2a,n2b),i3=n3a,n3b)
      write(iofile,*) 'u2=',((((u2(i1,i2,i3,m),m=0,2),i1=m1a,m1b),i2=m2a,m2b),i3=m3a,m3b)

      end if
     
      ! *** do this for now --- assume grids have equal spacing
      dx(0)=dx1(0)
      dx(1)=dx1(1)
      dx(2)=dx1(2)

      dr(0)=dr1(0)
      dr(1)=dr1(1)
      dr(2)=dr1(2)

      epsx=1.e-20  ! fix this 


      numGhost=orderOfAccuracy/2


      ! bounds for loops that include ghost points in the tangential directions:
      nn1a=n1a
      nn1b=n1b
      nn2a=n2a
      nn2b=n2b
      nn3a=n3a
      nn3b=n3b

      mm1a=m1a
      mm1b=m1b
      mm2a=m2a
      mm2b=m2b
      mm3a=m3a
      mm3b=m3b


      i3=n3a
      j3=m3a

      axis1p1=mod(axis1+1,nd)
      axis2p1=mod(axis2+1,nd)

      is1=0
      is2=0
      is3=0

      do m=0,2
        id1(m)=0
        id2(m)=0
        id3(m)=0
        jd1(m)=0
        jd2(m)=0
        jd3(m)=0
      end do

      if( axis1.eq.0 ) then
        is1=1-2*side1
        id1(axis1)=is1
        ! include ghost lines in tangential directions (for extrapolating) :
        if( extrapTangential(boundaryCondition1(0,1)) )then 
          nn2a=nn2a-numGhost
          nn2b=nn2b+numGhost
        end if
        if( nd.eq.3 .and. extrapTangential(boundaryCondition1(0,2)) )then
          nn3a=nn3a-numGhost
          nn3b=nn3b+numGhost
        end if
      else if( axis1.eq.1 )then
        is2=1-2*side1
        id2(axis1)=is2
        if( extrapTangential(boundaryCondition1(0,0)) )then
          nn1a=nn1a-numGhost
          nn1b=nn1b+numGhost
        end if
        if( nd.eq.3 .and. extrapTangential(boundaryCondition1(0,2)) )then
          nn3a=nn3a-numGhost
          nn3b=nn3b+numGhost
        end if
      else if( axis1.eq.2 )then
        is3=1-2*side1
        id3(axis1)=is3
        if( extrapTangential(boundaryCondition1(0,0)) )then
          nn1a=nn1a-numGhost
          nn1b=nn1b+numGhost
        end if
        if( extrapTangential(boundaryCondition1(0,1)) )then
          nn2a=nn2a-numGhost
          nn2b=nn2b+numGhost
        end if
      else
        ! invalid value for axis1
        stop 1143
      end if


      js1=0
      js2=0
      js3=0
      if( axis2.eq.0 ) then
        js1=1-2*side2
        jd1(axis2)=js1
        if( extrapTangential(boundaryCondition2(0,1)) )then
          mm2a=mm2a-numGhost
          mm2b=mm2b+numGhost
        end if
        if( nd.eq.3 .and. extrapTangential(boundaryCondition2(0,2)) )then
          mm3a=mm3a-numGhost
          mm3b=mm3b+numGhost
        end if
      else if( axis2.eq.1 ) then
        js2=1-2*side2
        jd2(axis2)=js2
        if( extrapTangential(boundaryCondition2(0,0)) )then
          mm1a=mm1a-numGhost
          mm1b=mm1b+numGhost
        end if
        if( nd.eq.3 .and. extrapTangential(boundaryCondition2(0,2)) )then
          mm3a=mm3a-numGhost
          mm3b=mm3b+numGhost
        end if
      else if( axis2.eq.2 ) then
        js3=1-2*side2
        jd3(axis2)=js3
        if( extrapTangential(boundaryCondition2(0,0)) )then
          mm1a=mm1a-numGhost
          mm1b=mm1b+numGhost
        end if
        if( extrapTangential(boundaryCondition2(0,1)) )then
          mm2a=mm2a-numGhost
          mm2b=mm2b+numGhost
        end if
      else
        ! invalid value for axis2
        stop 1144
      end if

      is=normalSign1  ! 1-2*side1
      js=normalSign2  ! 1-2*side2


      if( debug.gt.3 )then
        write(iofile,'("nn1a,nn1b,...=",6i5)') nn1a,nn1b,nn2a,nn2b,nn3a,nn3b
        write(iofile,'("mm1a,mm1b,...=",6i5)') mm1a,mm1b,mm2a,mm2b,mm3a,mm3b

      end if


      if( nd.eq.2 )then

        ! *********************************** 
        ! **************** 2D ***************
        ! *********************************** 

       if( orderOfAccuracy.eq.2 .and. gridType.eq.rectangular )then
  
         ! testing: just copy values to ghost points from first line in 
         !beginLoops2d()
         !  u1(i1-is1,i2-is2,i3,tc1)=u2(j1+js1,j2+js2,j3,tc2)
         !  u2(j1-js1,j2-js2,j3,tc2)=u1(i1+is1,i2+is2,i3,tc1)
         !endLoops2d()

        
         ! begin by satisfying the boundary jump conditions, [u]=0 
         !   --> Set the boundary values to the average the two values .
         beginGhostLoops2d()
!kkc 080624            uAve=.5*(u1(i1,i2,i3,tc1)+u2(j1,j2,j3,tc2))
            uAve=(ktc1*u1(i1,i2,i3,tc1) + ktc2*u2(j1,j2,j3,tc2))/ktcAve
            u1(i1,i2,i3,tc1)=uAve
            u2(j1,j2,j3,tc2)=uAve
         endLoops2d()
         if( useForcing.ne.0 )then
           ! Here is how we fix up the interface value when there is twilightzone: 
           ! [ u ] = g        (the twilightzone functions across the interface may be different)
           ! [ k u.x ] = 0    (this is may not be true for TZ but we use this condition to derived the TZ forcing)
           ! ->  u2i-u1i = g 
           ! ->  k2*( u2-u2i ) = k1*( u1-u1i ) 
           ! -> u1i = (k1*u1+k2*u2)/(k1+k2)  - k2*g/(k1+k2) 
           ! -> u2i = (k1*u1+k2*u2)/(k1+k2)  + k1*g/(k1+k2) 
           beginGhostLoops2d()
             getExact(ep1,xy1,i1,i2,i3,tc1, u1e)
             getExact(ep2,xy2,j1,j2,j3,tc2, u2e)
             u1(i1,i2,i3,tc1) = u1(i1,i2,i3,tc1) - ktc2*(u2e-u1e)/ktcAve 
             u2(j1,j2,j3,tc2) = u2(j1,j2,j3,tc2) + ktc1*(u2e-u1e)/ktcAve 
           endLoops2d()
         end if

         ! Solve the Temperature jump conditions:
         !    [ k*ux ] = 0
         !    [ kappa*(uxx+uyy) ] = 0


         a11= ktc1/(2.*dx1(axis1))
         a12= ktc2/(2.*dx2(axis2))
         a21= kappa1/(dx1(axis1)**2)
         a22=-kappa2/(dx2(axis2)**2)
         det=a11*a22-a12*a21
         f(0)=0.
         f(1)=0.
         an1=0.
         an2=0.
         if( axis1.eq.0 )then
           an1=1.
         else
           an2=1.
         end if
         beginLoops2d()
           ! eval Laplacian with wrong values on the ghost points:
           ulap1=u1Laplacian22r(i1,i2,i3,tc1)
           ulap2=u2Laplacian22r(j1,j2,j3,tc2)

           if( useForcing.ne.0 )then
             getDerivs(ep1,xy1, i1,i2,i3,tc1, u1ex,u1ey,u1ez, u1exx,u1eyy,u1ezz)
             getDerivs(ep2,xy2, j1,j2,j3,tc2, u2ex,u2ey,u2ez, u2exx,u2eyy,u2ezz)
             f(0) = (an1*u1ex+an2*u1ey)*ktc1 - (an1*u2ex+an2*u2ey)*ktc2
             f(1) = (u1exx+u1eyy)*kappa1 - (u2exx+u2eyy)*kappa2
           end if

           b1 =-is*ktc1*u1(i1+is1,i2+is2,i3,tc1)/(2.*dx1(axis1)) + \
                js*ktc2*u2(j1+js1,j2+js2,j3,tc2)/(2.*dx2(axis2)) + f(0)
           b2 =-kappa1*(ulap1-u1(i1-is1,i2-is2,i3,tc1)/dx1(axis1)**2) \
               +kappa2*(ulap2-u2(j1-js1,j2-js2,j3,tc2)/dx2(axis2)**2) + f(1)

           u1(i1-is1,i2-is2,i3,tc1)=( b1*a22-b2*a12)/det
           u2(j1-js1,j2-js2,j3,tc2)=(-b1*a21+b2*a11)/det
           if( debug.gt.3 .and. useForcing.ne.0 )then ! re-evaluate

             do n=0,1
               rr(n) = (aa(n,0)*f(0)+aa(n,1)*f(1)) - ff(n)
             end do

             ulap1=u1Laplacian22r(i1,i2,i3,tc1)
             ulap2=u2Laplacian22r(j1,j2,j3,tc2)

             f(0) = (an1*u1x22r(i1,i2,i3,tc1)+an2*u1y22r(i1,i2,i3,tc1))*ktc1 -\
                    (an1*u2x22r(j1,j2,j3,tc2)+an2*u2y22r(j1,j2,j3,tc2))*ktc2
             f(1) = ulap1*kappa1 - ulap2*kappa2

             getDerivs(ep1,xy1, i1,i2,i3,tc1, u1ex,u1ey,u1ez, u1exx,u1eyy,u1ezz)
             getDerivs(ep2,xy2, j1,j2,j3,tc2, u2ex,u2ey,u2ez, u2exx,u2eyy,u2ezz)
             f(0) = f(0) - (  (an1*u1ex+an2*u1ey)*ktc1 - (an1*u2ex+an2*u2ey)*ktc2 )
             f(1) = f(1) - ( (u1exx+u1eyy)*kappa1 - (u2exx+u2eyy)*kappa2 )

             getExact(ep1,xy1,i1,i2,i3,tc1, u1e)
             getExact(ep2,xy2,j1,j2,j3,tc2, u2e)

             write(iofile,'(" --> order2-rect: i1,i2=",2i4," u1,u2=",2f6.3," [k*grad],[kappa*Lap]=",4e10.2)') i1,i2,u1(i1,i2,i3,tc1),u2(j1,j2,j3,tc2),f(0),f(1)
             write(iofile,'("   t=",e10.3," u1e,u2e=",2f6.3)') t,u1e,u2e
             ! '
           end if


         endLoops2d()

       else if( orderOfAccuracy.eq.2 .and. gridType.eq.curvilinear )then
  
         ! testing: just copy values to ghost points from first line in 
         ! beginLoops2d()
         !   u1(i1-is1,i2-is2,i3,tc1)=u2(j1+js1,j2+js2,j3,tc2)
         !   u2(j1-js1,j2-js2,j3,tc2)=u1(i1+is1,i2+is2,i3,tc1)
         ! endLoops2d()
        

         ! begin by satisfying the boundary jump conditions, [u]=0 
         !   --> Set the boundary values to the average the two values .
         ! Also extrapolate ghost values -- this may be needed if the grid is not orthogonal
         beginGhostLoops2d()
!kkc 080624            uAve=.5*(u1(i1,i2,i3,tc1)+u2(j1,j2,j3,tc2))
            uAve=(ktc1*u1(i1,i2,i3,tc1) + ktc2*u2(j1,j2,j3,tc2))/ktcAve
            u1(i1,i2,i3,tc1)=uAve
            u2(j1,j2,j3,tc2)=uAve
            u1(i1-is1,i2-is2,i3,tc1)=extrap3(u1,i1,i2,i3,tc1,is1,is2,is3)
            u2(j1-js1,j2-js2,j3,tc2)=extrap3(u2,j1,j2,j3,tc2,js1,js2,js3)
         endLoops2d()
         if( useForcing.ne.0 )then
           beginGhostLoops2d()
             getExact(ep1,xy1,i1,i2,i3,tc1, u1e)
             getExact(ep2,xy2,j1,j2,j3,tc2, u2e)
             u1(i1,i2,i3,tc1) = u1(i1,i2,i3,tc1) - ktc2*(u2e-u1e)/ktcAve 
             u2(j1,j2,j3,tc2) = u2(j1,j2,j3,tc2) + ktc1*(u2e-u1e)/ktcAve 
           endLoops2d()
         end if

  
         ! assign the extended boundary in the tangential direction  ! ***** finish this ******
         if( .false. )then
         do side=0,1
           ! (i1,i2,i3) = ghost point in the tangential direction
           i1=gridIndexRange1(side,0)-id1(axis1p1)
           i2=gridIndexRange1(side,1)-id2(axis1p1)
           i3=gridIndexRange1(side,2)-id3(axis1p1)
           if( boundaryCondition1(side,axis1p1).lt.0 )then
             ! periodic
             i1p= i1 + id1(axis1)*(gridIndexRange1(1,0)-gridIndexRange1(0,0))
             i2p= i2 + id2(axis1)*(gridIndexRange1(1,1)-gridIndexRange1(0,1))
             i3p= i3 + id3(axis1)*(gridIndexRange1(1,2)-gridIndexRange1(0,2))
             u1(i1,i2,i3,tc1)=u1(i1p,i2p,i3p,tc1)
           else
             ! extrap
             ia1=id1(axis1)
             ia2=id2(axis1)
             ia3=id3(axis1)
             u1(i1,i2,i3,tc1)=extrap3(u1,i1,i2,i3,tc1,ia1,ia2,ia3)
           end if

         end do
         end if

         beginLoops2d()

           ! Solve the Temperature jump conditions:
           !    [ kappa*n.grad(u) ] = 0
           !    [ kappa*(uxx+uyy) ] = 0

           ! here is the normal (assumed to be the same on both sides)
           an1=rsxy1(i1,i2,i3,axis1,0)   ! normal (an1,an2)
           an2=rsxy1(i1,i2,i3,axis1,1)
           aNorm=max(epsx,sqrt(an1**2+an2**2))
           an1=an1/aNorm
           an2=an2/aNorm

           ulap1=u1Laplacian22(i1,i2,i3,tc1)
           ulap2=u2Laplacian22(j1,j2,j3,tc2)

           f(0) = (an1*u1x22(i1,i2,i3,tc1)+an2*u1y22(i1,i2,i3,tc1))*ktc1 -\
                  (an1*u2x22(j1,j2,j3,tc2)+an2*u2y22(j1,j2,j3,tc2))*ktc2
           f(1) = ulap1*kappa1 - ulap2*kappa2


           if( useForcing.ne.0 )then

             getDerivs(ep1,xy1, i1,i2,i3,tc1, u1ex,u1ey,u1ez, u1exx,u1eyy,u1ezz)
             getDerivs(ep2,xy2, j1,j2,j3,tc2, u2ex,u2ey,u2ez, u2exx,u2eyy,u2ezz)
             f(0) = f(0) - (  (an1*u1ex+an2*u1ey)*ktc1 - (an1*u2ex+an2*u2ey)*ktc2 )
             f(1) = f(1) - ( (u1exx+u1eyy)*kappa1 - (u2exx+u2eyy)*kappa2 )

           end if

           a2(0,0)=-is*(an1*rsxy1(i1,i2,i3,axis1,0)+an2*rsxy1(i1,i2,i3,axis1,1))*ktc1/(2.*dr1(axis1))
           a2(0,1)= js*(an1*rsxy2(j1,j2,j3,axis2,0)+an2*rsxy2(j1,j2,j3,axis2,1))*ktc2/(2.*dr2(axis2))

           ! coeff of u(-1) from lap = u.xx + u.yy
           clap1=(rsxy1(i1,i2,i3,axis1,0)**2+rsxy1(i1,i2,i3,axis1,1)**2)/(dr1(axis1)**2) \
                     -is*(rsxy1x22(i1,i2,i3,axis1,0)+rsxy1y22(i1,i2,i3,axis1,1))/(2.*dr1(axis1))
           clap2=(rsxy2(j1,j2,j3,axis2,0)**2+rsxy2(j1,j2,j3,axis2,1)**2)/(dr2(axis2)**2) \
                       -js*(rsxy2x22(j1,j2,j3,axis2,0)+rsxy2y22(j1,j2,j3,axis2,1))/(2.*dr2(axis2))

           a2(1,0)= clap1*kappa1
           a2(1,1)=-clap2*kappa2

           q(0) = u1(i1-is1,i2-is2,i3,tc1)
           q(1) = u2(j1-js1,j2-js2,j3,tc2)

           ! subtract off the contributions from the wrong values at the ghost points:
           do n=0,1
             f(n) = (a2(n,0)*q(0)+a2(n,1)*q(1)) - f(n)
           end do

           if( debug.gt.3 )then
             write(iofile,'("before solve t=",e10.3," ulap1,uLap2=",2e10.2)') t,ulap1,uLap2
             write(iofile,'("before solve t=",e10.3," an1,an2=",2e10.2)') t,an1,an2
             write(iofile,'("before solve t=",e10.3," f(0),f(1)=",2e10.2)') t,f(0),f(1)
             write(iofile,'("before solve a2=",4e10.2)') a2(0,0),a2(1,0),a2(0,1),a2(1,1)
           endif

           ! save for testing:
c           aa(0,0)=a2(0,0)
c           aa(0,1)=a2(0,1)
c           aa(1,0)=a2(1,0)
c           aa(1,1)=a2(1,1)
c           ff(0)=f(0)
c           ff(1)=f(1)

           call dgeco( a2(0,0), 2, 2, ipvt(0),rcond,work(0))
           job=0
           call dgesl( a2(0,0), 2, 2, ipvt(0), f(0), job)

           if( debug.gt.3 )then
             write(iofile,'("after solve rcond=",e10.3," f(0),f(1)=",2e10.2)') rcond,f(0),f(1)
             ! ' 
           endif

           u1(i1-is1,i2-is2,i3,tc1)=f(0)
           u2(j1-js1,j2-js2,j3,tc2)=f(1)

           if( debug.gt.3 )then ! re-evaluate

             do n=0,1
               rr(n) = (aa(n,0)*f(0)+aa(n,1)*f(1)) - ff(n)
             end do

             ulap1=u1Laplacian22(i1,i2,i3,tc1)
             ulap2=u2Laplacian22(j1,j2,j3,tc2)

             f(0) = (an1*u1x22(i1,i2,i3,tc1)+an2*u1y22(i1,i2,i3,tc1))*ktc1 -\
                    (an1*u2x22(j1,j2,j3,tc2)+an2*u2y22(j1,j2,j3,tc2))*ktc2
             f(1) = ulap1*kappa1 - ulap2*kappa2

             write(iofile,'("after solve   t=",e10.3," ulap1,uLap2=",2f6.3)') t,ulap1,uLap2
             ! '

             if( useForcing.ne.0 )then
               getDerivs(ep1,xy1, i1,i2,i3,tc1, u1ex,u1ey,u1ez, u1exx,u1eyy,u1ezz)
               getDerivs(ep2,xy2, j1,j2,j3,tc2, u2ex,u2ey,u2ez, u2exx,u2eyy,u2ezz)

               write(iofile,'("   t=",e10.3," u1ex,u1ey,u1ez=",3f6.3)') t,u1ex,u1ey,u1ez
               write(iofile,'("   t=",e10.3," u2ex,u2ey,u2ez=",3f6.3)') t,u2ex,u2ey,u2ez

               f(0) = f(0) - (  (an1*u1ex+an2*u1ey)*ktc1 - (an1*u2ex+an2*u2ey)*ktc2 )
               f(1) = f(1) - ( (u1exx+u1eyy)*kappa1 - (u2exx+u2eyy)*kappa2 )
             end if

             getExact(ep1,xy1,i1,i2,i3,tc1, u1e)
             getExact(ep2,xy2,j1,j2,j3,tc2, u2e)

             write(iofile,'(" --> order2-curv: i1,i2=",2i4," u1,u2=",2f6.3," [k*grad],[k*Lap]=",4e10.2)') i1,i2,u1(i1,i2,i3,tc1),u2(j1,j2,j3,tc2),f(0),f(1)
             write(iofile,'("   t=",e10.3," u1e,u2e=",2f6.3)') t,u1e,u2e
             write(iofile,'("   xy1=",2e10.3)') xy1(i1,i2,i3,0),xy1(i1,i2,i3,1)
             write(iofile,'("   xy2=",2e10.3)') xy2(j1,j2,j3,0),xy2(j1,j2,j3,1)
             write(iofile,'(" rsxy1=",4e10.3)') rsxy1(i1,i2,i3,0,0),rsxy1(i1,i2,i3,1,0),rsxy1(i1,i2,i3,0,1),rsxy1(i1,i2,i3,1,1)
             write(iofile,'(" rsxy2=",4e10.3)') rsxy2(j1,j2,j3,0,0),rsxy2(j1,j2,j3,1,0),rsxy2(j1,j2,j3,0,1),rsxy2(j1,j2,j3,1,1)
c             write(iofile,'("     rcond=",e9.2," resid=",2e8.1," clap1,clap2=",4e10.2)') rcond,rr(0),rr(1),clap1,clap2
c             write(iofile,'("     aa=",4e9.2," ff=",4e9.2)') aa(0,0),aa(0,1),aa(1,0),aa(1,1),ff(0),ff(1)
             ! '
           end if

         endLoops2d()

       else
         stop 6663
       end if

      else if( nd.eq.3 )then
       !    *************************
       !    ********** 3D ***********
       !    *************************

       if( orderOfAccuracy.eq.2 .and. gridType.eq.rectangular )then
  
         ! begin by satisfying the boundary jump conditions, [u]=0 
         !   --> Set the boundary values to the average the two values .
         beginGhostLoops3d()
!kkc 080624            uAve=.5*(u1(i1,i2,i3,tc1)+u2(j1,j2,j3,tc2))
            uAve=(ktc1*u1(i1,i2,i3,tc1) + ktc2*u2(j1,j2,j3,tc2))/ktcAve
            u1(i1,i2,i3,tc1)=uAve
            u2(j1,j2,j3,tc2)=uAve
         endLoops3d()
         if( useForcing.ne.0 )then
           beginGhostLoops3d()
             getExact(ep1,xy1,i1,i2,i3,tc1, u1e)
             getExact(ep2,xy2,j1,j2,j3,tc2, u2e)
             u1(i1,i2,i3,tc1) = u1(i1,i2,i3,tc1) - ktc2*(u2e-u1e)/ktcAve 
             u2(j1,j2,j3,tc2) = u2(j1,j2,j3,tc2) + ktc1*(u2e-u1e)/ktcAve 
           endLoops3d()
         end if

         ! Solve the Temperature jump conditions:
         !    [ kappa*ux ] = 0
         !    [ kappa*(uxx+uyy) ] = 0
         an1=0.
         an2=0.
         an3=0.
         if( axis1.eq.0 )then
           an1=1.
         else if( axis1.eq.1 )then
           an2=1.
         else
           an3=1.
         end if
         f(0)=0.
         f(1)=0.
         a11= ktc1/(2.*dx1(axis1))
         a12= ktc2/(2.*dx2(axis2))
         a21= kappa1/(dx1(axis1)**2)
         a22=-kappa2/(dx2(axis2)**2)
         det=a11*a22-a12*a21
         beginLoops3d()
           ! eval Laplacian with wrong values on the ghost points:
           ulap1=u1Laplacian23r(i1,i2,i3,tc1)
           ulap2=u2Laplacian23r(j1,j2,j3,tc2)

           if( useForcing.ne.0 )then
             getDerivs(ep1,xy1, i1,i2,i3,tc1, u1ex,u1ey,u1ez, u1exx,u1eyy,u1ezz)
             getDerivs(ep2,xy2, j1,j2,j3,tc2, u2ex,u2ey,u2ez, u2exx,u2eyy,u2ezz)
             f(0) = (an1*u1ex+an2*u1ey+an3*u1ez)*ktc1 - (an1*u2ex+an2*u2ey+an3*u2ez)*ktc2
             f(1) = (u1exx+u1eyy+u1ezz)*kappa1 - (u2exx+u2eyy+u2ezz)*kappa2
           end if

           b1 =-is*ktc1*u1(i1+is1,i2+is2,i3+is3,tc1)/(2.*dx1(axis1)) + \
                js*ktc2*u2(j1+js1,j2+js2,j3+js3,tc2)/(2.*dx2(axis2)) + f(0)
           b2 =-kappa1*(ulap1-u1(i1-is1,i2-is2,i3-is3,tc1)/dx1(axis1)**2) \
               +kappa2*(ulap2-u2(j1-js1,j2-js2,j3-js3,tc2)/dx2(axis2)**2) + f(1)

           u1(i1-is1,i2-is2,i3-is3,tc1)=( b1*a22-b2*a12)/det
           u2(j1-js1,j2-js2,j3-js3,tc2)=(-b1*a21+b2*a11)/det
         endLoops3d()

       else if( orderOfAccuracy.eq.2 .and. gridType.eq.curvilinear )then
  
         ! begin by satisfying the boundary jump conditions, [u]=0 
         !   --> Set the boundary values to the average the two values .
         ! Also extrapolate ghost values -- this may be need if the grid is not orthogonal
         beginGhostLoops3d()
!kkc 080624             uAve=.5*(u1(i1,i2,i3,tc1)+u2(j1,j2,j3,tc2))
            uAve=(ktc1*u1(i1,i2,i3,tc1) + ktc2*u2(j1,j2,j3,tc2))/ktcAve
            u1(i1,i2,i3,tc1)=uAve
            u2(j1,j2,j3,tc2)=uAve
            u1(i1-is1,i2-is2,i3-is3,tc1)=extrap3(u1,i1,i2,i3,tc1,is1,is2,is3)
            u2(j1-js1,j2-js2,j3-js3,tc2)=extrap3(u2,j1,j2,j3,tc2,js1,js2,js3)

         endLoops3d()
         if( useForcing.ne.0 )then
           beginGhostLoops3d()
             getExact(ep1,xy1,i1,i2,i3,tc1, u1e)
             getExact(ep2,xy2,j1,j2,j3,tc2, u2e)
             u1(i1,i2,i3,tc1) = u1(i1,i2,i3,tc1) - ktc2*(u2e-u1e)/ktcAve 
             u2(j1,j2,j3,tc2) = u2(j1,j2,j3,tc2) + ktc1*(u2e-u1e)/ktcAve 
c      write(iofile,'(" --> 3d-curv: i1,i2,i3=",3i3," j1,j2,j3=",3i3," u1,u2=",2f6.3," u1e,u2e=",2f6.3)') i1,i2,i3,j1,j2,j3,u1(i1,i2,i3,tc1),u2(j1,j2,j3,tc2),u1e,u2e
      ! '

           endLoops3d()
         end if

         beginLoops3d()

           ! Solve the Temperature jump conditions:
           !    [ kappa*n.grad(u) ] = 0
           !    [ kappa*(uxx+uyy+uzz) ] = 0

           ! here is the normal (assumed to be the same on both sides)
           an1=rsxy1(i1,i2,i3,axis1,0)   ! normal (an1,an2,an3)
           an2=rsxy1(i1,i2,i3,axis1,1)
           an3=rsxy1(i1,i2,i3,axis1,2)
           aNorm=max(epsx,sqrt(an1**2+an2**2+an3**2))
           an1=an1/aNorm
           an2=an2/aNorm
           an3=an3/aNorm

           ulap1=u1Laplacian23(i1,i2,i3,tc1)
           ulap2=u2Laplacian23(j1,j2,j3,tc2)

           f(0) = (an1*u1x23(i1,i2,i3,tc1)+an2*u1y23(i1,i2,i3,tc1)+an3*u1z23(i1,i2,i3,tc1))*ktc1 -\
                  (an1*u2x23(j1,j2,j3,tc2)+an2*u2y23(j1,j2,j3,tc2)+an3*u2z23(j1,j2,j3,tc2))*ktc2
           f(1) = ulap1*kappa1 - ulap2*kappa2

           if( useForcing.ne.0 )then
             getDerivs(ep1,xy1, i1,i2,i3,tc1, u1ex,u1ey,u1ez, u1exx,u1eyy,u1ezz)
             getDerivs(ep2,xy2, j1,j2,j3,tc2, u2ex,u2ey,u2ez, u2exx,u2eyy,u2ezz)
             f(0) = f(0) - (  (an1*u1ex+an2*u1ey+an3*u1ez)*ktc1 - (an1*u2ex+an2*u2ey+an3*u2ez)*ktc2 )
             f(1) = f(1) - ( (u1exx+u1eyy+u1ezz)*kappa1 - (u2exx+u2eyy+u2ezz)*kappa2 )
           end if


           a2(0,0)=-is*(an1*rsxy1(i1,i2,i3,axis1,0)+\
                        an2*rsxy1(i1,i2,i3,axis1,1)+\
                        an3*rsxy1(i1,i2,i3,axis1,2))*ktc1/(2.*dr1(axis1))
           a2(0,1)= js*(an1*rsxy2(j1,j2,j3,axis2,0)+\
                        an2*rsxy2(j1,j2,j3,axis2,1)+\
                        an3*rsxy2(j1,j2,j3,axis2,2))*ktc2/(2.*dr2(axis2))

           ! coeff of u(-1) from lap = u.xx + u.yy + u.zz
           clap1=(rsxy1(i1,i2,i3,axis1,0)**2+rsxy1(i1,i2,i3,axis1,1)**2+rsxy1(i1,i2,i3,axis1,2)**2)/(dr1(axis1)**2) \
               -is*(rsxy1x23(i1,i2,i3,axis1,0)+rsxy1y23(i1,i2,i3,axis1,1)+rsxy1z23(i1,i2,i3,axis1,2))/(2.*dr1(axis1))
           clap2=(rsxy2(j1,j2,j3,axis2,0)**2+rsxy2(j1,j2,j3,axis2,1)**2+rsxy2(j1,j2,j3,axis2,2)**2)/(dr2(axis2)**2) \
               -js*(rsxy2x23(j1,j2,j3,axis2,0)+rsxy2y23(j1,j2,j3,axis2,1)+rsxy2z23(j1,j2,j3,axis2,2))/(2.*dr2(axis2))

           a2(1,0)= clap1*kappa1
           a2(1,1)=-clap2*kappa2

           q(0) = u1(i1-is1,i2-is2,i3-is3,tc1)
           q(1) = u2(j1-js1,j2-js2,j3-js3,tc2)

           ! subtract off the contributions from the wrong values at the ghost points:
           do n=0,1
             f(n) = (a2(n,0)*q(0)+a2(n,1)*q(1)) - f(n)
           end do

           call dgeco( a2(0,0), 2, 2, ipvt(0),rcond,work(0))
           job=0
           call dgesl( a2(0,0), 2, 2, ipvt(0), f(0), job)

           u1(i1-is1,i2-is2,i3-is3,tc1)=f(0)
           u2(j1-js1,j2-js2,j3-js3,tc2)=f(1)

           if( debug.gt.3 .and. useForcing.ne.0 )then ! re-evaluate

             do n=0,1
               rr(n) = (aa(n,0)*f(0)+aa(n,1)*f(1)) - ff(n)
             end do

             ulap1=u1Laplacian23(i1,i2,i3,tc1)
             ulap2=u2Laplacian23(j1,j2,j3,tc2)

             f(0) = (an1*u1x23(i1,i2,i3,tc1)+an2*u1y23(i1,i2,i3,tc1)+an3*u1z23(i1,i2,i3,tc1))*ktc1 -\
                    (an1*u2x23(j1,j2,j3,tc2)+an2*u2y23(j1,j2,j3,tc2)+an3*u2z23(j1,j2,j3,tc2))*ktc2
             f(1) = ulap1*kappa1 - ulap2*kappa2

             getDerivs(ep1,xy1, i1,i2,i3,tc1, u1ex,u1ey,u1ez, u1exx,u1eyy,u1ezz)
             getDerivs(ep2,xy2, j1,j2,j3,tc2, u2ex,u2ey,u2ez, u2exx,u2eyy,u2ezz)
             f(0) = f(0) - (  (an1*u1ex+an2*u1ey+an3*u1ez)*ktc1 - (an1*u2ex+an2*u2ey+an3*u2ez)*ktc2 )
             f(1) = f(1) - ( (u1exx+u1eyy+u1ezz)*kappa1 - (u2exx+u2eyy+u2ezz)*kappa2 )

             getExact(ep1,xy1,i1,i2,i3,tc1, u1e)
             getExact(ep2,xy2,j1,j2,j3,tc2, u2e)

             write(iofile,'(" --> 3d-curv: i1,i2,i3=",3i4," u1,u2=",2f6.3," [k*grad],[kappa*Lap]=",4e10.2)') i1,i2,i3,u1(i1,i2,i3,tc1),u2(j1,j2,j3,tc2),f(0),f(1)
             write(iofile,'("   t=",e10.3," u1e,u2e=",2f6.3)') t,u1e,u2e
             ! '
           end if

         endLoops3d()


       else
         stop 6675
       end if

      else  
         ! 3D
        stop 6676
      end if


c -----------------------------
c**       if( orderOfAccuracy.eq.2 .and. gridType.eq.rectangular )then
c**  
c**        if( .true. )then
c**
c**         ! just copy values from ghost points for now
c**         beginLoops2d()
c**           u1(i1-is1,i2-is2,i3,ex)=u2(j1+js1,j2+js2,j3,ex)
c**           u1(i1-is1,i2-is2,i3,ey)=u2(j1+js1,j2+js2,j3,ey)
c**           u1(i1-is1,i2-is2,i3,hz)=u2(j1+js1,j2+js2,j3,hz) 
c**
c**           u2(j1-js1,j2-js2,j3,ex)=u1(i1+is1,i2+is2,i3,ex)
c**           u2(j1-js1,j2-js2,j3,ey)=u1(i1+is1,i2+is2,i3,ey)
c**           u2(j1-js1,j2-js2,j3,hz)=u1(i1+is1,i2+is2,i3,hz)
c**         endLoops2d()
c**
c**       else
c**
c**         ! ---- first satisfy the jump conditions on the boundary --------
c**         !    [ eps n.u ] = 0
c**         !    [ tau.u ] = 0
c**         boundaryJumpConditions(2,rectangular)
c**
c**         ! initialization step: assign first ghost line by extrapolation
c**         ! NOTE: assign ghost points outside the ends
c**         beginGhostLoops2d()
c**            u1(i1-is1,i2-is2,i3,ex)=extrap3(u1,i1,i2,i3,ex,is1,is2,is3)
c**            u1(i1-is1,i2-is2,i3,ey)=extrap3(u1,i1,i2,i3,ey,is1,is2,is3)
c**            u1(i1-is1,i2-is2,i3,hz)=extrap3(u1,i1,i2,i3,hz,is1,is2,is3)
c**c
c**            u2(j1-js1,j2-js2,j3,ex)=extrap3(u2,j1,j2,j3,ex,js1,js2,js3)
c**            u2(j1-js1,j2-js2,j3,ey)=extrap3(u2,j1,j2,j3,ey,js1,js2,js3)
c**            u2(j1-js1,j2-js2,j3,hz)=extrap3(u2,j1,j2,j3,hz,js1,js2,js3)
c**
c**         endLoops2d()
c**
c**         ! here are the real jump conditions
c**         !   [ u.x + v.y ] = 0
c**         !   [ u.xx + u.yy ] = 0
c**         !   [ v.x - u.y ] =0 
c**         !   [ (v.xx+v.yy)/eps ] = 0
c**         beginLoops2d()
c**           ! first evaluate the equations we want to solve with the wrong values at the ghost points:
c**           f(0)=(u1x22r(i1,i2,i3,ex)+u1y22r(i1,i2,i3,ey)) - \
c**                (u2x22r(j1,j2,j3,ex)+u2y22r(j1,j2,j3,ey))
c**           f(1)=(u1xx22r(i1,i2,i3,ex)+u1yy22r(i1,i2,i3,ex)) - \
c**                (u2xx22r(j1,j2,j3,ex)+u2yy22r(j1,j2,j3,ex))
c**
c**           f(2)=(u1x22r(i1,i2,i3,ey)-u1y22r(i1,i2,i3,ex)) - \
c**                (u2x22r(j1,j2,j3,ey)-u2y22r(j1,j2,j3,ex))
c**           
c**           f(3)=(u1xx22r(i1,i2,i3,ey)+u1yy22r(i1,i2,i3,ey))/eps1 - \
c**                (u2xx22r(j1,j2,j3,ey)+u2yy22r(j1,j2,j3,ey))/eps2
c**    
c**      ! write(*,'(" --> i1,i2=",2i4," f(start)=",4f8.3)') i1,i2,f(0),f(1),f(2),f(3)
c**
c**           ! here is the matrix of coefficients for the unknowns u1(-1),v1(-1),u2(-1),v2(-1)
c**           ! Solve:
c**           !     
c**           !       A [ U ] = A [ U(old) ] - [ f ]
c**           if( axis1.eq.0 )then
c**             a4(0,0) = -is1/(2.*dx1(0))    ! coeff of u1(-1) from [u.x+v.y] 
c**             a4(0,1) = 0.                  ! coeff of v1(-1) from [u.x+v.y] 
c**           
c**             a4(2,0) = 0.
c**             a4(2,1) = -is1/(2.*dx1(0))    ! coeff of v1(-1) from [v.x - u.y] 
c**           else 
c**             a4(0,0) = 0.                 
c**             a4(0,1) = -is2/(2.*dx1(1))    ! coeff of v1(-1) from [u.x+v.y] 
c**
c**             a4(2,0) =  is2/(2.*dx1(1))    ! coeff of u1(-1) from [v.x - u.y] 
c**             a4(2,1) = 0.
c**           end if
c**           if( axis2.eq.0 )then
c**             a4(0,2) = js1/(2.*dx2(0))    ! coeff of u2(-1) from [u.x+v.y] 
c**             a4(0,3) = 0. 
c**           
c**             a4(2,2) = 0.
c**             a4(2,3) = js1/(2.*dx2(0))    ! coeff of v2(-1) from [v.x - u.y]
c**           else
c**             a4(0,2) = 0. 
c**             a4(0,3) = js2/(2.*dx2(1))    ! coeff of v2(-1) from [u.x+v.y] 
c**
c**             a4(2,2) =-js2/(2.*dx2(1))    ! coeff of u2(-1) from [v.x - u.y] 
c**             a4(2,3) = 0.
c**           end if
c**
c**           a4(1,0) = 1./(dx1(axis1)**2)   ! coeff of u1(-1) from [u.xx + u.yy]
c**           a4(1,1) = 0. 
c**           a4(1,2) =-1./(dx2(axis2)**2)   ! coeff of u2(-1) from [u.xx + u.yy]
c**           a4(1,3) = 0. 
c**             
c**           a4(3,0) = 0.                      
c**           a4(3,1) = 1./(dx1(axis1)**2)/eps1 ! coeff of v1(-1) from [(v.xx+v.yy)/eps]
c**           a4(3,2) = 0. 
c**           a4(3,3) =-1./(dx2(axis2)**2)/eps2 ! coeff of v2(-1) from [(v.xx+v.yy)/eps]
c**             
c**
c**           q(0) = u1(i1-is1,i2-is2,i3,ex)
c**           q(1) = u1(i1-is1,i2-is2,i3,ey)
c**           q(2) = u2(j1-js1,j2-js2,j3,ex)
c**           q(3) = u2(j1-js1,j2-js2,j3,ey)
c**
c**           ! subtract off the contributions from the wrong values at the ghost points:
c**           do n=0,3
c**             f(n) = (a4(n,0)*q(0)+a4(n,1)*q(1)+a4(n,2)*q(2)+a4(n,3)*q(3)) - f(n)
c**           end do
c**      ! write(*,'(" --> i1,i2=",2i4," f(subtract)=",4f8.3)') i1,i2,f(0),f(1),f(2),f(3)
c**           ! solve A Q = F
c**           ! factor the matrix
c**           numberOfEquations=4
c**           call dgeco( a4(0,0), numberOfEquations, numberOfEquations, ipvt(0),rcond,work(0))
c**           ! solve
c**      ! write(*,'(" --> i1,i2=",2i4," rcond=",e10.2)') i1,i2,rcond
c**           job=0
c**           call dgesl( a4(0,0), numberOfEquations, numberOfEquations, ipvt(0), f(0), job)
c**      ! write(*,'(" --> i1,i2=",2i4," f(solve)=",4f8.3)') i1,i2,f(0),f(1),f(2),f(3)
c**
c**           u1(i1-is1,i2-is2,i3,ex)=f(0)
c**           u1(i1-is1,i2-is2,i3,ey)=f(1)
c**           u2(j1-js1,j2-js2,j3,ex)=f(2)
c**           u2(j1-js1,j2-js2,j3,ey)=f(3)
c**
c**      if( debug.gt.3 )then ! re-evaluate
c**           f(0)=(u1x22r(i1,i2,i3,ex)+u1y22r(i1,i2,i3,ey)) - \
c**                (u2x22r(j1,j2,j3,ex)+u2y22r(j1,j2,j3,ey))
c**           f(1)=(u1xx22r(i1,i2,i3,ex)+u1yy22r(i1,i2,i3,ex)) - \
c**                (u2xx22r(j1,j2,j3,ex)+u2yy22r(j1,j2,j3,ex))
c**
c**           f(2)=(u1x22r(i1,i2,i3,ey)-u1y22r(i1,i2,i3,ex)) - \
c**                (u2x22r(j1,j2,j3,ey)-u2y22r(j1,j2,j3,ex))
c**           
c**           f(3)=(u1xx22r(i1,i2,i3,ey)+u1yy22r(i1,i2,i3,ey))/eps1 - \
c**                (u2xx22r(j1,j2,j3,ey)+u2yy22r(j1,j2,j3,ey))/eps2
c**    
c**        write(*,'(" --> i1,i2=",2i4," f(re-eval)=",4e10.2)') i1,i2,f(0),f(1),f(2),f(3)
c**      end if
c**
c**           ! do this for now
c**           u1(i1-is1,i2-is2,i3,hz)=u2(j1+js1,j2+js2,j3,hz) 
c**           u2(j1-js1,j2-js2,j3,hz)=u1(i1+is1,i2+is2,i3,hz)
c**
c**
c**         endLoops2d()
c**
c**         ! periodic update
c**         periodicUpdate2d(u1,boundaryCondition1,gridIndexRange1,side1,axis1)
c**         periodicUpdate2d(u2,boundaryCondition2,gridIndexRange2,side2,axis2)
c**
c**       end if
c**
c**       else if( orderOfAccuracy.eq.2 .and. gridType.eq.curvilinear )then
c**         ! ***** curvilinear case *****
c**
c**         ! ---- first satisfy the jump conditions on the boundary --------
c**         !    [ eps n.u ] = 0
c**         !    [ tau.u ] = 0
c**         boundaryJumpConditions(2,curvilinear)
c**
c**         ! initialization step: assign first ghost line by extrapolation
c**         ! NOTE: assign ghost points outside the ends
c**         beginGhostLoops2d()
c**            u1(i1-is1,i2-is2,i3,ex)=extrap3(u1,i1,i2,i3,ex,is1,is2,is3)
c**            u1(i1-is1,i2-is2,i3,ey)=extrap3(u1,i1,i2,i3,ey,is1,is2,is3)
c**            u1(i1-is1,i2-is2,i3,hz)=extrap3(u1,i1,i2,i3,hz,is1,is2,is3)
c**c
c**            u2(j1-js1,j2-js2,j3,ex)=extrap3(u2,j1,j2,j3,ex,js1,js2,js3)
c**            u2(j1-js1,j2-js2,j3,ey)=extrap3(u2,j1,j2,j3,ey,js1,js2,js3)
c**            u2(j1-js1,j2-js2,j3,hz)=extrap3(u2,j1,j2,j3,hz,js1,js2,js3)
c**
c**         endLoops2d()
c**
c**         ! here are the real jump conditions for the ghost points
c**         !   [ u.x + v.y ] = 0 = [ rx*ur + ry*vr + sx*us + sy*vs ] 
c**         !   [ n.(uv.xx + uv.yy) ] = 0
c**         !   [ v.x - u.y ] =0 
c**         !   [ tau.(uv.xx+uv.yy)/eps ] = 0
c**
c**
c**         beginLoops2d()
c**
c**           ! here is the normal (assumed to be the same on both sides)
c**           an1=rsxy1(i1,i2,i3,axis1,0)   ! normal (an1,an2)
c**           an2=rsxy1(i1,i2,i3,axis1,1)
c**           aNorm=max(epsx,sqrt(an1**2+an2**2))
c**           an1=an1/aNorm
c**           an2=an2/aNorm
c**           tau1=-an2
c**           tau2= an1
c**
c**           ulap1=u1Laplacian22(i1,i2,i3,ex)
c**           vlap1=u1Laplacian22(i1,i2,i3,ey)
c**           ulap2=u2Laplacian22(j1,j2,j3,ex)
c**           vlap2=u2Laplacian22(j1,j2,j3,ey)
c**
c**           ! first evaluate the equations we want to solve with the wrong values at the ghost points:
c**           if( giveDiv.eq.0 )then
c**             f(0)=(u1x22(i1,i2,i3,ex)+u1y22(i1,i2,i3,ey)) - \
c**                  (u2x22(j1,j2,j3,ex)+u2y22(j1,j2,j3,ey))
c**             f(1)=( an1*ulap1 +an2*vlap1 )- \
c**                  ( an1*ulap2 +an2*vlap2 )
c**           else
c**             ! *** give div(u)=0 on both sides ***
c**             f(0)=u1x22(i1,i2,i3,ex)+u1y22(i1,i2,i3,ey)
c**             f(1)=u2x22(j1,j2,j3,ex)+u2y22(j1,j2,j3,ey)
c**           end if
c**
c**           f(2)=(u1x22(i1,i2,i3,ey)-u1y22(i1,i2,i3,ex)) - \
c**                (u2x22(j1,j2,j3,ey)-u2y22(j1,j2,j3,ex))
c**           
c**           f(3)=( tau1*ulap1 +tau2*vlap1 )/eps1 - \
c**                ( tau1*ulap2 +tau2*vlap2 )/eps2
c**    
c**      ! write(*,'(" --> order2-curv: i1,i2=",2i4," f(start)=",4f8.3)') i1,i2,f(0),f(1),f(2),f(3)
c**
c**           ! here is the matrix of coefficients for the unknowns u1(-1),v1(-1),u2(-1),v2(-1)
c**           ! Solve:
c**           !     
c**           !       A [ U ] = A [ U(old) ] - [ f ]
c**           if( giveDiv.eq.0 )then
c**             a4(0,0) = -is*rsxy1(i1,i2,i3,axis1,0)/(2.*dr1(axis1))    ! coeff of u1(-1) from [u.x+v.y] 
c**             a4(0,1) = -is*rsxy1(i1,i2,i3,axis1,1)/(2.*dr1(axis1))    ! coeff of v1(-1) from [u.x+v.y] 
c**             a4(0,2) =  js*rsxy2(j1,j2,j3,axis2,0)/(2.*dr2(axis2))    ! coeff of u2(-1) from [u.x+v.y] 
c**             a4(0,3) =  js*rsxy2(j1,j2,j3,axis2,1)/(2.*dr2(axis2))    ! coeff of v2(-1) from [u.x+v.y] 
c**           else
c**             a4(0,0) = -is*rsxy1(i1,i2,i3,axis1,0)/(2.*dr1(axis1))    ! coeff of u1(-1) from u.x+v.y=0
c**             a4(0,1) = -is*rsxy1(i1,i2,i3,axis1,1)/(2.*dr1(axis1))    ! coeff of v1(-1) from u.x+v.y=0
c**             a4(0,2) =  0.
c**             a4(0,3) =  0.
c**
c**             a4(1,0) = 0.
c**             a4(1,1) = 0.
c**             a4(1,2) = -js*rsxy2(j1,j2,j3,axis2,0)/(2.*dr2(axis2))    ! coeff of u2(-1) from u.x+v.y=0
c**             a4(1,3) = -js*rsxy2(j1,j2,j3,axis2,1)/(2.*dr2(axis2))    ! coeff of v2(-1) from u.x+v.y=0
c**           end if
c**
c**           a4(2,0) =  is*rsxy1(i1,i2,i3,axis1,1)/(2.*dr1(axis1))   ! coeff of u1(-1) from [v.x - u.y] 
c**           a4(2,1) = -is*rsxy1(i1,i2,i3,axis1,0)/(2.*dr1(axis1))   ! coeff of v1(-1) from [v.x - u.y] 
c**
c**           a4(2,2) = -js*rsxy2(j1,j2,j3,axis2,1)/(2.*dr2(axis2))   ! coeff of u2(-1) from [v.x - u.y] 
c**           a4(2,3) =  js*rsxy2(j1,j2,j3,axis2,0)/(2.*dr2(axis2))   ! coeff of v2(-1) from [v.x - u.y] 
c**
c**
c**           ! coeff of u(-1) from lap = u.xx + u.yy
c**           clap1=(rsxy1(i1,i2,i3,axis1,0)**2+rsxy1(i1,i2,i3,axis1,1)**2)/(dr1(axis1)**2) \
c**                     -is*(rsxy1x22(i1,i2,i3,axis1,0)+rsxy1y22(i1,i2,i3,axis1,1))/(2.*dr1(axis1))
c**           clap2=(rsxy2(j1,j2,j3,axis2,0)**2+rsxy2(j1,j2,j3,axis2,1)**2)/(dr2(axis2)**2) \
c**                       -js*(rsxy2x22(j1,j2,j3,axis2,0)+rsxy2y22(j1,j2,j3,axis2,1))/(2.*dr2(axis2)) 
c**
c**           !   [ n.(uv.xx + u.yy) ] = 0
c**           if( giveDiv.eq.0 )then
c**             a4(1,0) = an1*clap1
c**             a4(1,1) = an2*clap1
c**             a4(1,2) =-an1*clap2
c**             a4(1,3) =-an2*clap2
c**           end if 
c**           !   [ tau.(uv.xx+uv.yy)/eps ] = 0
c**           a4(3,0) = tau1*clap1/eps1
c**           a4(3,1) = tau2*clap1/eps1
c**           a4(3,2) =-tau1*clap2/eps2
c**           a4(3,3) =-tau2*clap2/eps2
c**             
c**
c**           q(0) = u1(i1-is1,i2-is2,i3,ex)
c**           q(1) = u1(i1-is1,i2-is2,i3,ey)
c**           q(2) = u2(j1-js1,j2-js2,j3,ex)
c**           q(3) = u2(j1-js1,j2-js2,j3,ey)
c**
c**           ! subtract off the contributions from the wrong values at the ghost points:
c**           do n=0,3
c**             f(n) = (a4(n,0)*q(0)+a4(n,1)*q(1)+a4(n,2)*q(2)+a4(n,3)*q(3)) - f(n)
c**           end do
c**      ! write(*,'(" --> order2-curv: i1,i2=",2i4," f(subtract)=",4f8.3)') i1,i2,f(0),f(1),f(2),f(3)
c**           ! solve A Q = F
c**           ! factor the matrix
c**           numberOfEquations=4
c**           call dgeco( a4(0,0), numberOfEquations, numberOfEquations, ipvt(0),rcond,work(0))
c**           ! solve
c**      !   write(*,'(" --> order2-curv: i1,i2=",2i4," rcond=",e10.2)') i1,i2,rcond
c**           job=0
c**           call dgesl( a4(0,0), numberOfEquations, numberOfEquations, ipvt(0), f(0), job)
c**      ! write(*,'(" --> order2-curv: i1,i2=",2i4," f(solve)=",4f8.3)') i1,i2,f(0),f(1),f(2),f(3)
c**
c**           u1(i1-is1,i2-is2,i3,ex)=f(0)
c**           u1(i1-is1,i2-is2,i3,ey)=f(1)
c**           u2(j1-js1,j2-js2,j3,ex)=f(2)
c**           u2(j1-js1,j2-js2,j3,ey)=f(3)
c**
c**           if( debug.gt.3 )then ! re-evaluate
c**             ulap1=u1Laplacian22(i1,i2,i3,ex)
c**             vlap1=u1Laplacian22(i1,i2,i3,ey)
c**             ulap2=u2Laplacian22(j1,j2,j3,ex)
c**             vlap2=u2Laplacian22(j1,j2,j3,ey)
c**  
c**             if( giveDiv.eq.0 )then
c**               f(0)=(u1x22(i1,i2,i3,ex)+u1y22(i1,i2,i3,ey)) - \
c**                    (u2x22(j1,j2,j3,ex)+u2y22(j1,j2,j3,ey))
c**               f(1)=( an1*ulap1 +an2*vlap1 )- \
c**                    ( an1*ulap2 +an2*vlap2 )
c**             else
c**               ! *** give div(u)=0 on both sides ***
c**               f(0)=u1x22(i1,i2,i3,ex)+u1y22(i1,i2,i3,ey)
c**               f(1)=u2x22(j1,j2,j3,ex)+u2y22(j1,j2,j3,ey)
c**             end if
c**             f(2)=(u1x22(i1,i2,i3,ey)-u1y22(i1,i2,i3,ex)) - \
c**                  (u2x22(j1,j2,j3,ey)-u2y22(j1,j2,j3,ex))
c**             f(3)=( tau1*ulap1 +tau2*vlap1 )/eps1 - \
c**                  ( tau1*ulap2 +tau2*vlap2 )/eps2
c**             write(*,'(" --> order2-curv: i1,i2=",2i4," f(re-eval)=",4e10.2)') i1,i2,f(0),f(1),f(2),f(3)
c**               ! '
c**           end if
c**
c**           ! solve for Hz
c**           !  [ w.n/eps] = 0
c**           !  [ Lap(w)/eps] = 0
c**
c**           wlap1=u1Laplacian22(i1,i2,i3,hz)
c**           wlap2=u2Laplacian22(j1,j2,j3,hz)
c**
c**           f(0) = (an1*u1x22(i1,i2,i3,hz)+an2*u1y22(i1,i2,i3,hz))/eps1 -\
c**                  (an1*u2x22(j1,j2,j3,hz)+an2*u2y22(j1,j2,j3,hz))/eps2
c**           f(1) = wlap1/eps1 - wlap2/eps2
c**
c**           a2(0,0)=-is*(an1*rsxy1(i1,i2,i3,axis1,0)+an2*rsxy1(i1,i2,i3,axis1,1))/(2.*dr1(axis1)*eps1)
c**           a2(0,1)= js*(an1*rsxy2(j1,j2,j3,axis2,0)+an2*rsxy2(j1,j2,j3,axis2,1))/(2.*dr2(axis2)*eps2)
c**
c**           a2(1,0)= clap1/eps1
c**           a2(1,1)=-clap2/eps2
c**
c**           q(0) = u1(i1-is1,i2-is2,i3,hz)
c**           q(1) = u2(j1-js1,j2-js2,j3,hz)
c**
c**           ! subtract off the contributions from the wrong values at the ghost points:
c**           do n=0,1
c**             f(n) = (a2(n,0)*q(0)+a2(n,1)*q(1)) - f(n)
c**           end do
c**
c**           call dgeco( a2(0,0), 2, 2, ipvt(0),rcond,work(0))
c**           job=0
c**           call dgesl( a2(0,0), 2, 2, ipvt(0), f(0), job)
c**
c**           u1(i1-is1,i2-is2,i3,hz)=f(0)
c**           u2(j1-js1,j2-js2,j3,hz)=f(1)
c**
c**           ! u1(i1-is1,i2-is2,i3,hz)=u2(j1+js1,j2+js2,j3,hz) 
c**           ! u2(j1-js1,j2-js2,j3,hz)=u1(i1+is1,i2+is2,i3,hz)
c**
c**           if( debug.gt.3 )then ! re-evaluate
c**
c**             wlap1=u1Laplacian22(i1,i2,i3,hz)
c**             wlap2=u2Laplacian22(j1,j2,j3,hz)
c**
c**             f(0) = (an1*u1x22(i1,i2,i3,hz)+an2*u1y22(i1,i2,i3,hz))/eps1 -\
c**                    (an1*u2x22(j1,j2,j3,hz)+an2*u2y22(j1,j2,j3,hz))/eps2
c**             f(1) = wlap1/eps1 - wlap2/eps2
c**
c**             write(*,'(" --> order2-curv: i1,i2=",2i4," hz-f(re-eval)=",4e10.2)') i1,i2,f(0),f(1)
c**               ! '
c**           end if
c**
c**         endLoops2d()
c**
c**         ! now make sure that div(u)=0 etc.
c**         if( .false. )then
c**         beginLoops2d() ! =============== start loops =======================
c**
c**           ! 0  [ u.x + v.y ] = 0
c**           ! first evaluate the equations we want to solve with the wrong values at the ghost points:
c**           divu=u1x22(i1,i2,i3,ex)+u1y22(i1,i2,i3,ey)
c**           a0=-is*rsxy1(i1,i2,i3,axis1,0)*dr112(axis1)
c**           a1=-is*rsxy1(i1,i2,i3,axis1,1)*dr112(axis1)
c**           aNormSq=a0**2+a1**2
c**           ! now project:  a.uNew = a.uOld - div  ->  (div-a.uOld)+a.uNew = div(uNew) = 0
c**           u1(i1-is1,i2-is2,i3,ex)=u1(i1-is1,i2-is2,i3,ex)-divu*a0/aNormSq
c**           u1(i1-is1,i2-is2,i3,ey)=u1(i1-is1,i2-is2,i3,ey)-divu*a1/aNormSq
c**
c**           divu=u2x22(j1,j2,j3,ex)+u2y22(j1,j2,j3,ey)
c**           a0=-js*rsxy2(j1,j2,j3,axis2,0)*dr212(axis2) 
c**           a1=-js*rsxy2(j1,j2,j3,axis2,1)*dr212(axis2) 
c**           aNormSq=a0**2+a1**2
c**
c**           u2(j1-js1,j2-js2,j3,ex)=u2(j1-js1,j2-js2,j3,ex)-divu*a0/aNormSq
c**           u2(j1-js1,j2-js2,j3,ey)=u2(j1-js1,j2-js2,j3,ey)-divu*a1/aNormSq
c**
c**           if( debug.gt.3 )then
c**             write(*,'(" --> 2cth: eval div1,div2=",2e10.2)') u1x22(i1,i2,i3,ex)+u1y22(i1,i2,i3,ey),u2x22(j1,j2,j3,ex)+u2y22(j1,j2,j3,ey)
c**           end if
c**         endLoops2d()
c**         end if
c**
c**         ! periodic update
c**         periodicUpdate2d(u1,boundaryCondition1,gridIndexRange1,side1,axis1)
c**         periodicUpdate2d(u2,boundaryCondition2,gridIndexRange2,side2,axis2)
c**
c**       else if( .false. .and. orderOfAccuracy.eq.4 )then
c**
c**         ! for testing -- just assign from the other ghost points
c**
c**         beginLoops2d()
c**           u1(i1-is1,i2-is2,i3,ex)=u2(j1+js1,j2+js2,j3,ex)
c**           u1(i1-is1,i2-is2,i3,ey)=u2(j1+js1,j2+js2,j3,ey)
c**           u1(i1-is1,i2-is2,i3,hz)=u2(j1+js1,j2+js2,j3,hz) 
c**
c**           u2(j1-js1,j2-js2,j3,ex)=u1(i1+is1,i2+is2,i3,ex)
c**           u2(j1-js1,j2-js2,j3,ey)=u1(i1+is1,i2+is2,i3,ey)
c**           u2(j1-js1,j2-js2,j3,hz)=u1(i1+is1,i2+is2,i3,hz)
c**
c**           u1(i1-2*is1,i2-2*is2,i3,ex)=u2(j1+2*js1,j2+2*js2,j3,ex)
c**           u1(i1-2*is1,i2-2*is2,i3,ey)=u2(j1+2*js1,j2+2*js2,j3,ey)
c**           u1(i1-2*is1,i2-2*is2,i3,hz)=u2(j1+2*js1,j2+2*js2,j3,hz) 
c**
c**           u2(j1-2*js1,j2-2*js2,j3,ex)=u1(i1+2*is1,i2+2*is2,i3,ex)
c**           u2(j1-2*js1,j2-2*js2,j3,ey)=u1(i1+2*is1,i2+2*is2,i3,ey)
c**           u2(j1-2*js1,j2-2*js2,j3,hz)=u1(i1+2*is1,i2+2*is2,i3,hz)
c**
c**         endLoops2d()
c**
c**       else if( orderOfAccuracy.eq.4 .and. gridType.eq.rectangular )then
c**  
c**         ! --------------- 4th Order Rectangular ---------------
c**         ! ---- first satisfy the jump conditions on the boundary --------
c**         !    [ eps n.u ] = 0
c**         !    [ tau.u ] = 0
c**         boundaryJumpConditions(2,rectangular)
c**
c**         ! here are the real jump conditions for the ghost points
c**         ! 0  [ u.x + v.y ] = 0
c**         ! 1  [ u.xx + u.yy ] = 0
c**         ! 2  [ v.x - u.y ] =0 
c**         ! 3  [ (v.xx+v.yy)/eps ] = 0
c**         ! 4  [ (u.xx+u.yy).x + (v.xx+v.yy).y ] = 0  OR [ (u.xx).x + (v.xx).y ] = 0 OR  [ (u.yy).x + (v.yy).y ] = 0 
c**         ! 5  [ {(Delta v).x - (Delta u).y}/eps ] =0  -> [ {(v.xxx+v.xyy)-(u.xxy+u.yyy)}/eps ] = 0
c**         ! 6  [ Delta^2 u/eps ] = 0
c**         ! 7  [ Delta^2 v/eps^2 ] = 0 
c**
c**
c**         ! initialization step: assign first ghost line by extrapolation
c**         ! NOTE: assign ghost points outside the ends
c**         beginGhostLoops2d()
c**           u1(i1-is1,i2-is2,i3,ex)=extrap4(u1,i1,i2,i3,ex,is1,is2,is3)
c**           u1(i1-is1,i2-is2,i3,ey)=extrap4(u1,i1,i2,i3,ey,is1,is2,is3)
c**           u1(i1-is1,i2-is2,i3,hz)=extrap4(u1,i1,i2,i3,hz,is1,is2,is3)
c**
c**           u2(j1-js1,j2-js2,j3,ex)=extrap4(u2,j1,j2,j3,ex,js1,js2,js3)
c**           u2(j1-js1,j2-js2,j3,ey)=extrap4(u2,j1,j2,j3,ey,js1,js2,js3)
c**           u2(j1-js1,j2-js2,j3,hz)=extrap4(u2,j1,j2,j3,hz,js1,js2,js3)
c**
c**           ! --- also extrap 2nd line for now
c**           ! u1(i1-2*is1,i2-2*is2,i3,ex)=extrap4(u1,i1-is1,i2-is2,i3,ex,is1,is2,is3)
c**           ! u1(i1-2*is1,i2-2*is2,i3,ey)=extrap4(u1,i1-is1,i2-is2,i3,ey,is1,is2,is3)
c**           ! u1(i1-2*is1,i2-2*is2,i3,hz)=extrap4(u1,i1-is1,i2-is2,i3,hz,is1,is2,is3)
c**
c**           ! u2(j1-2*js1,j2-2*js2,j3,ex)=extrap4(u2,j1-js1,j2-js2,j3,ex,js1,js2,js3)
c**           ! u2(j1-2*js1,j2-2*js2,j3,ey)=extrap4(u2,j1-js1,j2-js2,j3,ey,js1,js2,js3)
c**           ! u2(j1-2*js1,j2-2*js2,j3,hz)=extrap4(u2,j1-js1,j2-js2,j3,hz,js1,js2,js3)
c**         endLoops2d()
c**
c**         beginLoops2d() ! =============== start loops =======================
c**
c**           ! first evaluate the equations we want to solve with the wrong values at the ghost points:
c**           f(0)=(u1x42r(i1,i2,i3,ex)+u1y42r(i1,i2,i3,ey)) - \
c**                (u2x42r(j1,j2,j3,ex)+u2y42r(j1,j2,j3,ey))
c**
c**           f(1)=(u1xx42r(i1,i2,i3,ex)+u1yy42r(i1,i2,i3,ex)) - \
c**                (u2xx42r(j1,j2,j3,ex)+u2yy42r(j1,j2,j3,ex))
c**
c**           f(2)=(u1x42r(i1,i2,i3,ey)-u1y42r(i1,i2,i3,ex)) - \
c**                (u2x42r(j1,j2,j3,ey)-u2y42r(j1,j2,j3,ex))
c**           
c**           f(3)=(u1xx42r(i1,i2,i3,ey)+u1yy42r(i1,i2,i3,ey))/eps1 - \
c**                (u2xx42r(j1,j2,j3,ey)+u2yy42r(j1,j2,j3,ey))/eps2
c**    
c**           ! These next we can do to 2nd order -- these need a value on the first ghost line --
c**           f(4)=(u1xxx22r(i1,i2,i3,ex)+u1xyy22r(i1,i2,i3,ex)+u1xxy22r(i1,i2,i3,ey)+u1yyy22r(i1,i2,i3,ey)) - \
c**                (u2xxx22r(j1,j2,j3,ex)+u2xyy22r(j1,j2,j3,ex)+u2xxy22r(j1,j2,j3,ey)+u2yyy22r(j1,j2,j3,ey))
c**
c**           f(5)=((u1xxx22r(i1,i2,i3,ey)+u1xyy22r(i1,i2,i3,ey))-(u1xxy22r(i1,i2,i3,ex)+u1yyy22r(i1,i2,i3,ex)))/eps1 - \
c**                ((u2xxx22r(j1,j2,j3,ey)+u2xyy22r(j1,j2,j3,ey))-(u2xxy22r(j1,j2,j3,ex)+u2yyy22r(j1,j2,j3,ex)))/eps2
c**
c**           f(6)=(u1LapSq22r(i1,i2,i3,ex))/eps1 - \
c**                (u2LapSq22r(j1,j2,j3,ex))/eps2
c**
c**           f(7)=(u1LapSq22r(i1,i2,i3,ey))/eps1**2 - \
c**                (u2LapSq22r(j1,j2,j3,ey))/eps2**2
c**           
c**       write(*,'(" --> 4th: j1,j2=",2i4," u1xx,u1yy,u2xx,u2yy=",4e10.2)') j1,j2,u1xx42r(i1,i2,i3,ex),\
c**           u1yy42r(i1,i2,i3,ex),u2xx42r(j1,j2,j3,ex),u2yy42r(j1,j2,j3,ex)
c**       write(*,'(" --> 4th: i1,i2=",2i4," f(start)=",8e10.2)') i1,i2,f(0),f(1),f(2),f(3),f(4),f(5),f(6),f(7)
c**
c**           ! here is the matrix of coefficients for the unknowns u1(-1),v1(-1),u2(-1),v2(-1)
c**           ! Solve:
c**           !     
c**           !       A [ U ] = A [ U(old) ] - [ f ]
c**c      u1x43r(i1,i2,i3,kd)=(8.*(u1(i1+1,i2,i3,kd)-u1(i1-1,i2,i3,kd))-(
c**c     & u1(i1+2,i2,i3,kd)-u1(i1-2,i2,i3,kd)))*dx141(0)
c**
c**
c**           ! 0  [ u.x + v.y ] = 0
c**           a8(0,0) = -is*8.*rx1*dx141(axis1)     ! coeff of u1(-1) from [u.x+v.y] 
c**           a8(0,1) = -is*8.*ry1*dx141(axis1)     ! coeff of v1(-1) from [u.x+v.y] 
c**           a8(0,4) =  is*rx1*dx141(axis1)        ! u1(-2)
c**           a8(0,5) =  is*ry1*dx141(axis1)        ! v1(-2) 
c**
c**           a8(0,2) =  js*8.*rx2*dx241(axis2)     ! coeff of u2(-1) from [u.x+v.y] 
c**           a8(0,3) =  js*8.*ry2*dx241(axis2) 
c**           a8(0,6) = -js*   rx2*dx241(axis2) 
c**           a8(0,7) = -js*   ry2*dx241(axis2) 
c**
c**           ! 1  [ u.xx + u.yy ] = 0
c**c      u1xx43r(i1,i2,i3,kd)=( -30.*u1(i1,i2,i3,kd)+16.*(u1(i1+1,i2,i3,
c**c     & kd)+u1(i1-1,i2,i3,kd))-(u1(i1+2,i2,i3,kd)+u1(i1-2,i2,i3,kd)) )*
c**c     & dx142(0)
c**           
c**           a8(1,0) = 16.*dx142(axis1)         ! coeff of u1(-1) from [u.xx + u.yy]
c**           a8(1,1) = 0. 
c**           a8(1,4) =    -dx142(axis1)         ! coeff of u1(-2) from [u.xx + u.yy]
c**           a8(1,5) = 0. 
c**
c**           a8(1,2) =-16.*dx242(axis2)         ! coeff of u2(-1) from [u.xx + u.yy]
c**           a8(1,3) = 0. 
c**           a8(1,6) =     dx242(axis2)         ! coeff of u2(-2) from [u.xx + u.yy]
c**           a8(1,7) = 0. 
c**
c**
c**           ! 2  [ v.x - u.y ] =0 
c**           a8(2,0) =  is*8.*ry1*dx141(axis1)
c**           a8(2,1) = -is*8.*rx1*dx141(axis1)    ! coeff of v1(-1) from [v.x - u.y] 
c**           a8(2,4) = -is*   ry1*dx141(axis1)
c**           a8(2,5) =  is*   rx1*dx141(axis1)
c**
c**           a8(2,2) = -js*8.*ry2*dx241(axis2)
c**           a8(2,3) =  js*8.*rx2*dx241(axis2)
c**           a8(2,6) =  js*   ry2*dx241(axis2)
c**           a8(2,7) = -js*   rx2*dx241(axis2)
c**
c**           ! 3  [ (v.xx+v.yy)/eps ] = 0
c**           a8(3,0) = 0.                      
c**           a8(3,1) = 16.*dx142(axis1)/eps1 ! coeff of v1(-1) from [(v.xx+v.yy)/eps]
c**           a8(3,4) = 0.                      
c**           a8(3,5) =    -dx142(axis1)/eps1 ! coeff of v1(-2) from [(v.xx+v.yy)/eps]
c**
c**           a8(3,2) = 0. 
c**           a8(3,3) =-16.*dx242(axis2)/eps2 ! coeff of v2(-1) from [(v.xx+v.yy)/eps]
c**           a8(3,6) = 0. 
c**           a8(3,7) =     dx242(axis2)/eps2 ! coeff of v2(-2) from [(v.xx+v.yy)/eps]
c**
c**           ! 4  [ (u.xx+u.yy).x + (v.xx+v.yy).y ] = 0
c**c     u1xxx22r(i1,i2,i3,kd)=(-2.*(u1(i1+1,i2,i3,kd)-u1(i1-1,i2,i3,kd))+
c**c    & (u1(i1+2,i2,i3,kd)-u1(i1-2,i2,i3,kd)) )*dx122(0)*dx112(0)
c**c    u1xxy22r(i1,i2,i3,kd)=( u1xx22r(i1,i2+1,i3,kd)-u1xx22r(i1,i2-1,
c**c     & i3,kd))/(2.*dx1(1))
c**c      u1yy23r(i1,i2,i3,kd)=(-2.*u1(i1,i2,i3,kd)+(u1(i1,i2+1,i3,kd)+u1(
c**c     & i1,i2-1,i3,kd)) )*dx122(1)
c**c     u1xyy22r(i1,i2,i3,kd)=( u1yy22r(i1+1,i2,i3,kd)-u1yy22r(i1-1,i2,
c**c     & i3,kd))/(2.*dx1(0))
c**          a8(4,0)= ( is*rx1*2.*dx122(axis1)*dx112(axis1)+is*rx1*2.*dx122(1)/(2.*dx1(0)))
c**          a8(4,1)= ( is*ry1*2.*dx122(axis1)*dx112(axis1)+is*ry1*2.*dx122(0)/(2.*dx1(1)))
c**          a8(4,4)= (-is*rx1   *dx122(axis1)*dx112(axis1) )  
c**          a8(4,5)= (-is*ry1   *dx122(axis1)*dx112(axis1))
c**
c**          a8(4,2)=-( js*rx2*2.*dx222(axis2)*dx212(axis2)+js*rx2*2.*dx222(1)/(2.*dx2(0)))
c**          a8(4,3)=-( js*ry2*2.*dx222(axis2)*dx212(axis2)+js*ry2*2.*dx222(0)/(2.*dx2(1)))
c**          a8(4,6)=-(-js*rx2   *dx222(axis2)*dx212(axis2))   
c**          a8(4,7)=-(-js*ry2   *dx222(axis2)*dx212(axis2))
c**
c**          ! 5  [ {(Delta v).x - (Delta u).y}/eps ] =0  -> [ {(v.xxx+v.xyy)-(u.xxy+u.yyy)}/eps ] = 0
c**
c**          a8(5,0)=-( is*ry1*2.*dx122(axis1)*dx112(axis1)+is*ry1*2.*dx122(0)/(2.*dx1(1)))/eps1
c**          a8(5,1)= ( is*rx1*2.*dx122(axis1)*dx112(axis1)+is*rx1*2.*dx122(1)/(2.*dx1(0)))/eps1
c**          a8(5,4)=-(-is*ry1   *dx122(axis1)*dx112(axis1))/eps1
c**          a8(5,5)= (-is*rx1   *dx122(axis1)*dx112(axis1))/eps1   
c**
c**          a8(5,2)= ( js*ry2*2.*dx222(axis2)*dx212(axis2)+js*ry2*2.*dx222(0)/(2.*dx2(1)))/eps2
c**          a8(5,3)=-( js*rx2*2.*dx222(axis2)*dx212(axis2)+js*rx2*2.*dx222(1)/(2.*dx2(0)))/eps2
c**          a8(5,6)= (-js*ry2   *dx222(axis2)*dx212(axis2))/eps2
c**          a8(5,7)=-(-js*rx2   *dx222(axis2)*dx212(axis2))/eps2   
c**
c**           ! 6  [ Delta^2 u/eps ] = 0
c**c     u1LapSq22r(i1,i2,i3,kd)= ( 6.*u1(i1,i2,i3,kd)- 4.*(u1(i1+1,i2,i3,
c**c    & kd)+u1(i1-1,i2,i3,kd))+(u1(i1+2,i2,i3,kd)+u1(i1-2,i2,i3,kd)) )
c**c    & /(dx1(0)**4)+( 6.*u1(i1,i2,i3,kd)-4.*(u1(i1,i2+1,i3,kd)+u1(i1,
c**c    & i2-1,i3,kd)) +(u1(i1,i2+2,i3,kd)+u1(i1,i2-2,i3,kd)) )/(dx1(1)**
c**c    & 4)+( 8.*u1(i1,i2,i3,kd)-4.*(u1(i1+1,i2,i3,kd)+u1(i1-1,i2,i3,kd)
c**c    & +u1(i1,i2+1,i3,kd)+u1(i1,i2-1,i3,kd))+2.*(u1(i1+1,i2+1,i3,kd)+
c**c    & u1(i1-1,i2+1,i3,kd)+u1(i1+1,i2-1,i3,kd)+u1(i1-1,i2-1,i3,kd)) )
c**c    & /(dx1(0)**2*dx1(1)**2)
c**
c**           a8(6,0) = -(4./(dx1(axis1)**4) +4./(dx1(0)**2*dx1(1)**2) )/eps1
c**           a8(6,1) = 0.
c**           a8(6,4) =   1./(dx1(axis1)**4)/eps1
c**           a8(6,5) = 0.
c**
c**           a8(6,2) = (4./(dx2(axis2)**4) +4./(dx1(0)**2*dx1(1)**2) )/eps2
c**           a8(6,3) = 0.
c**           a8(6,6) =  -1./(dx2(axis2)**4)/eps2
c**           a8(6,7) = 0.
c**
c**           ! 7  [ Delta^2 v/eps^2 ] = 0 
c**           a8(7,0) = 0.
c**           a8(7,1) = -(4./(dx1(axis1)**4) +4./(dx2(0)**2*dx2(1)**2) )/eps1**2
c**           a8(7,4) = 0.
c**           a8(7,5) =   1./(dx1(axis1)**4)/eps1**2
c**
c**           a8(7,2) = 0.
c**           a8(7,3) =  (4./(dx2(axis2)**4) +4./(dx2(0)**2*dx2(1)**2) )/eps2**2
c**           a8(7,6) = 0.
c**           a8(7,7) =  -1./(dx2(axis2)**4)/eps2**2
c**
c**           q(0) = u1(i1-is1,i2-is2,i3,ex)
c**           q(1) = u1(i1-is1,i2-is2,i3,ey)
c**           q(2) = u2(j1-js1,j2-js2,j3,ex)
c**           q(3) = u2(j1-js1,j2-js2,j3,ey)
c**
c**           q(4) = u1(i1-2*is1,i2-2*is2,i3,ex)
c**           q(5) = u1(i1-2*is1,i2-2*is2,i3,ey)
c**           q(6) = u2(j1-2*js1,j2-2*js2,j3,ex)
c**           q(7) = u2(j1-2*js1,j2-2*js2,j3,ey)
c**
c**       write(*,'(" --> 4th: i1,i2=",2i4," q=",8e10.2)') i1,i2,q(0),q(1),q(2),q(3),q(4),q(5),q(6),q(7)
c**
c**           ! subtract off the contributions from the initial (wrong) values at the ghost points:
c**           do n=0,7
c**             f(n) = (a8(n,0)*q(0)+a8(n,1)*q(1)+a8(n,2)*q(2)+a8(n,3)*q(3)+\
c**                     a8(n,4)*q(4)+a8(n,5)*q(5)+a8(n,6)*q(6)+a8(n,7)*q(7)) - f(n)
c**           end do
c**
c**           ! solve A Q = F
c**           ! factor the matrix
c**           numberOfEquations=8
c**           call dgeco( a8(0,0), numberOfEquations, numberOfEquations, ipvt(0),rcond,work(0))
c**           ! solve
c**       write(*,'(" --> 4th: i1,i2=",2i4," rcond=",e10.2)') i1,i2,rcond
c**           job=0
c**           call dgesl( a8(0,0), numberOfEquations, numberOfEquations, ipvt(0), f(0), job)
c**
c**       write(*,'(" --> 4th: i1,i2=",2i4," f(solve)=",8e10.2)') i1,i2,f(0),f(1),f(2),f(3),f(4),f(5),f(6),f(7)
c**
c**           if( .true. )then
c**           u1(i1-is1,i2-is2,i3,ex)=f(0)
c**           u1(i1-is1,i2-is2,i3,ey)=f(1)
c**           u2(j1-js1,j2-js2,j3,ex)=f(2)
c**           u2(j1-js1,j2-js2,j3,ey)=f(3)
c**
c**           u1(i1-2*is1,i2-2*is2,i3,ex)=f(4)
c**           u1(i1-2*is1,i2-2*is2,i3,ey)=f(5)
c**           u2(j1-2*js1,j2-2*js2,j3,ex)=f(6)
c**           u2(j1-2*js1,j2-2*js2,j3,ey)=f(7)
c**           end if
c**
c**          if( debug.gt.3 )then ! re-evaluate
c**           f(0)=(u1x42r(i1,i2,i3,ex)+u1y42r(i1,i2,i3,ey)) - \
c**                (u2x42r(j1,j2,j3,ex)+u2y42r(j1,j2,j3,ey))
c**           f(1)=(u1xx42r(i1,i2,i3,ex)+u1yy42r(i1,i2,i3,ex)) - \
c**                (u2xx42r(j1,j2,j3,ex)+u2yy42r(j1,j2,j3,ex))
c**
c**           f(2)=(u1x42r(i1,i2,i3,ey)-u1y42r(i1,i2,i3,ex)) - \
c**                (u2x42r(j1,j2,j3,ey)-u2y42r(j1,j2,j3,ex))
c**           
c**           f(3)=(u1xx42r(i1,i2,i3,ey)+u1yy42r(i1,i2,i3,ey))/eps1 - \
c**                (u2xx42r(j1,j2,j3,ey)+u2yy42r(j1,j2,j3,ey))/eps2
c**    
c**           ! These next we can do to 2nd order -- these need a value on the first ghost line --
c**           f(4)=(u1xxx22r(i1,i2,i3,ex)+u1xyy22r(i1,i2,i3,ex)+u1xxy22r(i1,i2,i3,ey)+u1yyy22r(i1,i2,i3,ey)) - \
c**                (u2xxx22r(j1,j2,j3,ex)+u2xyy22r(j1,j2,j3,ex)+u2xxy22r(j1,j2,j3,ey)+u2yyy22r(j1,j2,j3,ey))
c**
c**           f(5)=((u1xxx22r(i1,i2,i3,ey)+u1xyy22r(i1,i2,i3,ey))-(u1xxy22r(i1,i2,i3,ex)+u1yyy22r(i1,i2,i3,ex)))/eps1 - \
c**                ((u2xxx22r(j1,j2,j3,ey)+u2xyy22r(j1,j2,j3,ey))-(u2xxy22r(j1,j2,j3,ex)+u2yyy22r(j1,j2,j3,ex)))/eps2
c**
c**           f(6)=(u1LapSq22r(i1,i2,i3,ex))/eps1 - \
c**                (u2LapSq22r(j1,j2,j3,ex))/eps2
c**
c**           f(7)=(u1LapSq22r(i1,i2,i3,ey))/eps1**2 - \
c**                (u2LapSq22r(j1,j2,j3,ey))/eps2**2
c**    
c**           write(*,'(" --> 4th: i1,i2=",2i4," f(re-eval)=",8e10.2)') i1,i2,f(0),f(1),f(2),f(3),f(4),f(5),f(6),f(7)
c**          end if
c**
c**           ! do this for now
c**           u1(i1-is1,i2-is2,i3,hz)=u2(j1+js1,j2+js2,j3,hz) 
c**           u2(j1-js1,j2-js2,j3,hz)=u1(i1+is1,i2+is2,i3,hz)
c**
c**           u1(i1-2*is1,i2-2*is2,i3,hz)=u2(j1+2*js1,j2+2*js2,j3,hz) 
c**           u2(j1-2*js1,j2-2*js2,j3,hz)=u1(i1+2*is1,i2+2*is2,i3,hz)
c**
c**         endLoops2d()
c**
c**         ! periodic update
c**         periodicUpdate2d(u1,boundaryCondition1,gridIndexRange1,side1,axis1)
c**         periodicUpdate2d(u2,boundaryCondition2,gridIndexRange2,side2,axis2)
c**
c**       else if( orderOfAccuracy.eq.4 .and. gridType.eq.curvilinear )then
c**  
c**         ! --------------- 4th Order Curvilinear ---------------
c**
c**         ! ---- first satisfy the jump conditions on the boundary --------
c**         !    [ eps n.u ] = 0
c**         !    [ tau.u ] = 0
c**         !    [ w ] = 0 
c**         boundaryJumpConditions(2,curvilinear)
c**
c**         ! here are the real jump conditions for the ghost points
c**         ! 0  [ u.x + v.y ] = 0
c**         ! 1  [ n.(uv.xx + uv.yy) ] = 0
c**         ! 2  [ v.x - u.y ] =0 
c**         ! 3  [ tau.(v.xx+v.yy)/eps ] = 0
c**         ! 4  [ (u.xx+u.yy).x + (v.xx+v.yy).y ] = 0  OR [ (u.xx).x + (v.xx).y ] = 0 OR  [ (u.yy).x + (v.yy).y ] = 0 
c**         ! 5  [ {(Delta v).x - (Delta u).y}/eps ] =0  -> [ {(v.xxx+v.xyy)-(u.xxy+u.yyy)}/eps ] = 0
c**         ! 6  [ n.Delta^2 uv/eps ] = 0
c**         ! 7  [ tau.Delta^2 uv/eps^2 ] = 0 
c**
c**
c**
c**         ! initialization step: assign first ghost line by extrapolation
c**         ! NOTE: assign ghost points outside the ends
c**
c**         
c**         beginGhostLoops2d()
c**
c**c          u1(i1-is1,i2-is2,i3,ex)=extrap2(u1,i1,i2,i3,ex,is1,is2,is3)
c**c          u1(i1-is1,i2-is2,i3,ey)=extrap2(u1,i1,i2,i3,ey,is1,is2,is3)
c**c          u1(i1-is1,i2-is2,i3,hz)=extrap2(u1,i1,i2,i3,hz,is1,is2,is3)
c**
c**c          u2(j1-js1,j2-js2,j3,ex)=extrap2(u2,j1,j2,j3,ex,js1,js2,js3)
c**c          u2(j1-js1,j2-js2,j3,ey)=extrap2(u2,j1,j2,j3,ey,js1,js2,js3)
c**c          u2(j1-js1,j2-js2,j3,hz)=extrap2(u2,j1,j2,j3,hz,js1,js2,js3)
c**c
c**            u1(i1-is1,i2-is2,i3,ex)=extrap4(u1,i1,i2,i3,ex,is1,is2,is3)
c**            u1(i1-is1,i2-is2,i3,ey)=extrap4(u1,i1,i2,i3,ey,is1,is2,is3)
c**            u1(i1-is1,i2-is2,i3,hz)=extrap4(u1,i1,i2,i3,hz,is1,is2,is3)
c**c
c**            u2(j1-js1,j2-js2,j3,ex)=extrap4(u2,j1,j2,j3,ex,js1,js2,js3)
c**            u2(j1-js1,j2-js2,j3,ey)=extrap4(u2,j1,j2,j3,ey,js1,js2,js3)
c**            u2(j1-js1,j2-js2,j3,hz)=extrap4(u2,j1,j2,j3,hz,js1,js2,js3)
c**
c**           ! --- also extrap 2nd line for now
c**           u1(i1-2*is1,i2-2*is2,i3,ex)=extrap4(u1,i1-is1,i2-is2,i3,ex,is1,is2,is3)
c**           u1(i1-2*is1,i2-2*is2,i3,ey)=extrap4(u1,i1-is1,i2-is2,i3,ey,is1,is2,is3)
c**           u1(i1-2*is1,i2-2*is2,i3,hz)=extrap4(u1,i1-is1,i2-is2,i3,hz,is1,is2,is3)
c**
c**           u2(j1-2*js1,j2-2*js2,j3,ex)=extrap4(u2,j1-js1,j2-js2,j3,ex,js1,js2,js3)
c**           u2(j1-2*js1,j2-2*js2,j3,ey)=extrap4(u2,j1-js1,j2-js2,j3,ey,js1,js2,js3)
c**           u2(j1-2*js1,j2-2*js2,j3,hz)=extrap4(u2,j1-js1,j2-js2,j3,hz,js1,js2,js3)
c**         endLoops2d()
c**
c**         ! write(*,'(">>> interface: order=4 initialized=",i4)') initialized
c**
c**         do it=1,nit ! *** begin iteration ****
c**
c**           err=0.
c**         ! =============== start loops ======================
c**         nn=-1 ! counts points on the interface
c**         beginLoops2d() 
c**
c**           nn=nn+1
c**
c**           ! here is the normal (assumed to be the same on both sides)
c**           an1=rsxy1(i1,i2,i3,axis1,0)   ! normal (an1,an2)
c**           an2=rsxy1(i1,i2,i3,axis1,1)
c**           aNorm=max(epsx,sqrt(an1**2+an2**2))
c**           an1=an1/aNorm
c**           an2=an2/aNorm
c**           tau1=-an2
c**           tau2= an1
c**
c**           ulap1=u1Laplacian42(i1,i2,i3,ex)
c**           vlap1=u1Laplacian42(i1,i2,i3,ey)
c**           ulap2=u2Laplacian42(j1,j2,j3,ex)
c**           vlap2=u2Laplacian42(j1,j2,j3,ey)
c**
c**           ulapSq1=u1LapSq22(i1,i2,i3,ex)
c**           vlapSq1=u1LapSq22(i1,i2,i3,ey)
c**           ulapSq2=u2LapSq22(j1,j2,j3,ex)
c**           vlapSq2=u2LapSq22(j1,j2,j3,ey)
c**
c**         
c**
c**           ! first evaluate the equations we want to solve with the wrong values at the ghost points:
c**           f(0)=(u1x42(i1,i2,i3,ex)+u1y42(i1,i2,i3,ey)) - \
c**                (u2x42(j1,j2,j3,ex)+u2y42(j1,j2,j3,ey))
c**
c**           f(1)=(an1*ulap1+an2*vlap1) - \
c**                (an1*ulap2+an2*vlap2)
c**
c**           f(2)=(u1x42(i1,i2,i3,ey)-u1y42(i1,i2,i3,ex)) - \
c**                (u2x42(j1,j2,j3,ey)-u2y42(j1,j2,j3,ex))
c**           
c**           f(3)=(tau1*ulap1+tau2*vlap1)/eps1 - \
c**                (tau1*ulap2+tau2*vlap2)/eps2
c**    
c**           ! These next we can do to 2nd order -- these need a value on the first ghost line --
c**           f(4)=(u1xxx22(i1,i2,i3,ex)+u1xyy22(i1,i2,i3,ex)+u1xxy22(i1,i2,i3,ey)+u1yyy22(i1,i2,i3,ey)) - \
c**                (u2xxx22(j1,j2,j3,ex)+u2xyy22(j1,j2,j3,ex)+u2xxy22(j1,j2,j3,ey)+u2yyy22(j1,j2,j3,ey))
c**
c**           f(5)=((u1xxx22(i1,i2,i3,ey)+u1xyy22(i1,i2,i3,ey))-(u1xxy22(i1,i2,i3,ex)+u1yyy22(i1,i2,i3,ex)))/eps1 - \
c**                ((u2xxx22(j1,j2,j3,ey)+u2xyy22(j1,j2,j3,ey))-(u2xxy22(j1,j2,j3,ex)+u2yyy22(j1,j2,j3,ex)))/eps2
c**
c**           f(6)=(an1*ulapSq1+an2*vlapSq1)/eps1 - \
c**                (an1*ulapSq2+an2*vlapSq2)/eps2
c**
c**           f(7)=(tau1*ulapSq1+tau2*vlapSq1)/eps1**2 - \
c**                (tau1*ulapSq2+tau2*vlapSq2)/eps2**2
c**           
c**       if( debug.gt.3 ) write(*,'(" --> 4cth: j1,j2=",2i4," u1xx,u1yy,u2xx,u2yy=",4e10.2)') j1,j2,u1xx42(i1,i2,i3,ex),\
c**           u1yy42(i1,i2,i3,ex),u2xx42(j1,j2,j3,ex),u2yy42(j1,j2,j3,ex)
c**       if( debug.gt.3 ) write(*,'(" --> 4cth: i1,i2=",2i4," f(start)=",8e10.2)') i1,i2,f(0),f(1),f(2),f(3),f(4),f(5),f(6),f(7)
c**
c**
c**c here are the macros from deriv.maple (file=derivMacros.h)
c**
c**#defineMacro lapCoeff4a(is,dr,ds) ( (-2/3.*rxx*is-2/3.*ryy*is)/dr+(4/3.*rx**2+4/3.*ry**2)/dr**2 )
c**
c**#defineMacro lapCoeff4b(is,dr,ds) ( (1/12.*rxx*is+1/12.*ryy*is)/dr+(-1/12.*rx**2-1/12.*ry**2)/dr**2 )
c**
c**#defineMacro xLapCoeff4a(is,dr,ds) ( (-1/2.*rxyy*is-1/2.*rxxx*is+(sy*(ry*sx*is+sy*rx*is)+3*rx*sx**2*is+ry*sy*sx*is)/ds**2)/dr+(2*ry*rxy+3*rx*rxx+ryy*rx)/dr**2+(ry**2*rx*is+rx**3*is)/dr**3 )
c**
c**#defineMacro xLapCoeff4b(is,dr,ds) ( (-1/2.*rx**3*is-1/2.*ry**2*rx*is)/dr**3 )
c**
c**#defineMacro yLapCoeff4a(is,dr,ds) ( (-1/2.*ryyy*is-1/2.*rxxy*is+(3*ry*sy**2*is+ry*sx**2*is+2*sy*rx*sx*is)/ds**2)/dr+(2*rxy*rx+ry*rxx+3*ry*ryy)/dr**2+(ry**3*is+ry*rx**2*is)/dr**3 )
c**
c**#defineMacro yLapCoeff4b(is,dr,ds) ( (-1/2.*ry*rx**2*is-1/2.*ry**3*is)/dr**3 )
c**
c**#defineMacro lapSqCoeff4a(is,dr,ds) ( (-1/2.*rxxxx*is-rxxyy*is-1/2.*ryyyy*is+(2*sy*(2*rxy*sx*is+2*rx*sxy*is)+2*ry*(2*sxy*sx*is+sy*sxx*is)+7*rx*sxx*sx*is+sy*(3*ry*syy*is+3*sy*ryy*is)+sx*(3*rx*sxx*is+3*rxx*sx*is)+sx*(2*rxx*sx*is+2*rx*sxx*is)+2*sy*(2*rx*sxy*is+ry*sxx*is+2*rxy*sx*is+sy*rxx*is)+7*ry*sy*syy*is+rxx*sx**2*is+4*ry*sxy*sx*is+4*syy*rx*sx*is+2*ryy*sx**2*is+ryy*sy**2*is+sy*(2*sy*ryy*is+2*ry*syy*is))/ds**2)/dr+(3*ryy**2+3*rxx**2+4*rxy**2+4*ry*rxxy+4*rx*rxxx+4*ry*ryyy+2*ryy*rxx+4*rx*rxyy+(2*ry*(-4*sy*rx*sx-2*ry*sx**2)-12*ry**2*sy**2+2*sy*(-2*sy*rx**2-4*ry*rx*sx)-12*rx**2*sx**2)/ds**2)/dr**2+(6*ry**2*ryy*is+4*ry*rxy*rx*is+2*ry*(ry*rxx*is+2*rxy*rx*is)+6*rxx*rx**2*is+2*ryy*rx**2*is)/dr**3+(-8*ry**2*rx**2-4*ry**4-4*rx**4)/dr**4 )
c**
c**#defineMacro lapSqCoeff4b(is,dr,ds) ( (-3*rxx*rx**2*is-ryy*rx**2*is-2*ry*rxy*rx*is-3*ry**2*ryy*is+2*ry*(-rxy*rx*is-1/2.*ry*rxx*is))/dr**3+(rx**4+2*ry**2*rx**2+ry**4)/dr**4 )
c**
c**
c**           ! here is the matrix of coefficients for the unknowns u1(-1),v1(-1),u2(-1),v2(-1)
c**           ! Solve:
c**           !     
c**           !       A [ U ] = A [ U(old) ] - [ f ]
c**c      u1r4(i1,i2,i3,kd)=(8.*(u1(i1+1,i2,i3,kd)-u1(i1-1,i2,i3,kd))-(u1(
c**c     & i1+2,i2,i3,kd)-u1(i1-2,i2,i3,kd)))*dr114(0)
c**c      u1x42(i1,i2,i3,kd)= rsxy1(i1,i2,i3,0,0)*u1r4(i1,i2,i3,kd)+rsxy1(
c**c     & i1,i2,i3,1,0)*u1s4(i1,i2,i3,kd)
c**c      u1y42(i1,i2,i3,kd)= rsxy1(i1,i2,i3,0,1)*u1r4(i1,i2,i3,kd)+rsxy1(
c**c     & i1,i2,i3,1,1)*u1s4(i1,i2,i3,kd)
c**c          a4(0,0) = -is*rsxy1(i1,i2,i3,axis1,0)/(2.*dr1(axis1))    ! coeff of u1(-1) from [u.x+v.y] 
c**c          a4(0,1) = -is*rsxy1(i1,i2,i3,axis1,1)/(2.*dr1(axis1))    ! coeff of v1(-1) from [u.x+v.y] 
c**c
c**c          a4(2,0) =  is*rsxy1(i1,i2,i3,axis1,1)/(2.*dr1(axis1))   ! coeff of u1(-1) from [v.x - u.y] 
c**c          a4(2,1) = -is*rsxy1(i1,i2,i3,axis1,0)/(2.*dr1(axis1))   ! coeff of v1(-1) from [v.x - u.y] 
c**c
c**c          a4(0,2) =  js*rsxy2(j1,j2,j3,axis2,0)/(2.*dr2(axis2))    ! coeff of u2(-1) from [u.x+v.y] 
c**c          a4(0,3) =  js*rsxy2(j1,j2,j3,axis2,1)/(2.*dr2(axis2))    ! coeff of v2(-1) from [u.x+v.y] 
c**c
c**c          a4(2,2) = -js*rsxy2(j1,j2,j3,axis2,1)/(2.*dr2(axis2))   ! coeff of u2(-1) from [v.x - u.y] 
c**c          a4(2,3) =  js*rsxy2(j1,j2,j3,axis2,0)/(2.*dr2(axis2))   ! coeff of v2(-1) from [v.x - u.y] 
c**
c**
c**           ! write(*,'(" interface:E: initialized,it=",2i4)') initialized,it
c**           if( .false. .or. (initialized.eq.0 .and. it.eq.1) )then
c**             ! form the matrix (and save factor for later use)
c**
c**             computeRxDerivatives(rv1,rsxy1,i1,i2,i3)
c**             computeRxDerivatives(rv2,rsxy2,j1,j2,j3)
c** 
c**             ! 0  [ u.x + v.y ] = 0
c**             aa8(0,0,0,nn) = -is*8.*rsxy1(i1,i2,i3,axis1,0)*dr114(axis1)     ! coeff of u1(-1) from [u.x+v.y] 
c**             aa8(0,1,0,nn) = -is*8.*rsxy1(i1,i2,i3,axis1,1)*dr114(axis1)     ! coeff of v1(-1) from [u.x+v.y] 
c**             aa8(0,4,0,nn) =  is*   rsxy1(i1,i2,i3,axis1,0)*dr114(axis1)     ! u1(-2)
c**             aa8(0,5,0,nn) =  is*   rsxy1(i1,i2,i3,axis1,1)*dr114(axis1)     ! v1(-2) 
c**  
c**             aa8(0,2,0,nn) =  js*8.*rsxy2(j1,j2,j3,axis2,0)*dr214(axis2)     ! coeff of u2(-1) from [u.x+v.y] 
c**             aa8(0,3,0,nn) =  js*8.*rsxy2(j1,j2,j3,axis2,1)*dr214(axis2)  
c**             aa8(0,6,0,nn) = -js*   rsxy2(j1,j2,j3,axis2,0)*dr214(axis2) 
c**             aa8(0,7,0,nn) = -js*   rsxy2(j1,j2,j3,axis2,1)*dr214(axis2)  
c**
c**           ! 1  [ u.xx + u.yy ] = 0
c**c this macro comes from deriv.maple
c**c return the coefficient of u(-1) in uxxx+uxyy
c**c#defineMacro lapCoeff4a(is,dr,ds) ((-1/3.*rxx*is-1/3.*ryy*is)/dr+(4/3.*rx**2+4/3.*ry**2)/dr**2)
c**
c**c return the coefficient of u(-2) in uxxx+uxyy
c**c#defineMacro lapCoeff4b(is,dr,ds) ((1/24.*rxx*is+1/24.*ryy*is)/dr+(-1/12.*rx**2-1/12.*ry**2)/dr**2 )
c**
c**             setJacobian(rv1,axis1,axis1p1)
c**             dr0=dr1(axis1)
c**             ds0=dr1(axis1p1)
c**             aLap0 = lapCoeff4a(is,dr0,ds0)
c**             aLap1 = lapCoeff4b(is,dr0,ds0)
c**  
c**             setJacobian(rv2,axis2,axis2p1)
c**             dr0=dr2(axis2)
c**             ds0=dr2(axis2p1)
c**             bLap0 = lapCoeff4a(js,dr0,ds0)
c**             bLap1 = lapCoeff4b(js,dr0,ds0)
c**  
c**            if( debug.gt.3 )then
c**             aa8(1,0,0,nn) = 16.*dx142(axis1)         ! coeff of u1(-1) from [u.xx + u.yy]
c**             aa8(1,4,0,nn) =    -dx142(axis1)         ! coeff of u1(-2) from [u.xx + u.yy]
c**              write(*,'(" 4th: lap4: aLap0: rect=",e12.4," curv=",e12.4)') aLap0,aa8(1,0,0,nn)
c**              write(*,'(" 4th: lap4: aLap1: rect=",e12.4," curv=",e12.4)') aLap1,aa8(1,4,0,nn)
c**            end if
c**  
c**             aa8(1,0,0,nn) = an1*aLap0       ! coeff of u1(-1) from [n.(u.xx + u.yy)]
c**             aa8(1,1,0,nn) = an2*aLap0 
c**             aa8(1,4,0,nn) = an1*aLap1       ! coeff of u1(-2) from [n.(u.xx + u.yy)]
c**             aa8(1,5,0,nn) = an2*aLap1  
c**             
c**             aa8(1,2,0,nn) =-an1*bLap0       ! coeff of u2(-1) from [n.(u.xx + u.yy)]
c**             aa8(1,3,0,nn) =-an2*bLap0
c**             aa8(1,6,0,nn) =-an1*bLap1       ! coeff of u2(-2) from [n.(u.xx + u.yy)]
c**             aa8(1,7,0,nn) =-an2*bLap1
c**  
c**           ! 2  [ v.x - u.y ] =0 
c**c          a8(2,0) =  is*8.*ry1*dx114(axis1)
c**c          a8(2,1) = -is*8.*rx1*dx114(axis1)    ! coeff of v1(-1) from [v.x - u.y] 
c**c          a8(2,4) = -is*   ry1*dx114(axis1)
c**c          a8(2,5) =  is*   rx1*dx114(axis1)
c**c          a8(2,2) = -js*8.*ry2*dx214(axis2)
c**c          a8(2,3) =  js*8.*rx2*dx214(axis2)
c**c          a8(2,6) =  js*   ry2*dx214(axis2)
c**c          a8(2,7) = -js*   rx2*dx214(axis2)
c**
c**             aa8(2,0,0,nn) =  is*8.*rsxy1(i1,i2,i3,axis1,1)*dr114(axis1)    
c**             aa8(2,1,0,nn) = -is*8.*rsxy1(i1,i2,i3,axis1,0)*dr114(axis1)    
c**             aa8(2,4,0,nn) = -is*   rsxy1(i1,i2,i3,axis1,1)*dr114(axis1)       
c**             aa8(2,5,0,nn) =  is*   rsxy1(i1,i2,i3,axis1,0)*dr114(axis1)       
c**  
c**             aa8(2,2,0,nn) = -js*8.*rsxy2(j1,j2,j3,axis2,1)*dr214(axis2)  
c**             aa8(2,3,0,nn) =  js*8.*rsxy2(j1,j2,j3,axis2,0)*dr214(axis2)    
c**             aa8(2,6,0,nn) =  js*   rsxy2(j1,j2,j3,axis2,1)*dr214(axis2)  
c**             aa8(2,7,0,nn) = -js*   rsxy2(j1,j2,j3,axis2,0)*dr214(axis2) 
c**  
c**             ! 3  [ tau.(uv.xx+uv.yy)/eps ] = 0
c**             aa8(3,0,0,nn) =tau1*aLap0/eps1
c**             aa8(3,1,0,nn) =tau2*aLap0/eps1
c**             aa8(3,4,0,nn) =tau1*aLap1/eps1
c**             aa8(3,5,0,nn) =tau2*aLap1/eps1
c**  
c**             aa8(3,2,0,nn) =-tau1*bLap0/eps2
c**             aa8(3,3,0,nn) =-tau2*bLap0/eps2
c**             aa8(3,6,0,nn) =-tau1*bLap1/eps2
c**             aa8(3,7,0,nn) =-tau2*bLap1/eps2
c**  
c**  
c**             ! 4  [ (u.xx+u.yy).x + (v.xx+v.yy).y ] = 0
c**  
c**            setJacobian(rv1,axis1,axis1p1)
c**  
c**  
c**            dr0=dr1(axis1)
c**            ds0=dr1(axis1p1)
c**            aLapX0 = xLapCoeff4a(is,dr0,ds0)
c**            aLapX1 = xLapCoeff4b(is,dr0,ds0)
c**  
c**            bLapY0 = yLapCoeff4a(is,dr0,ds0)
c**            bLapY1 = yLapCoeff4b(is,dr0,ds0)
c**  
c**            setJacobian(rv2,axis2,axis2p1)
c**  
c**            dr0=dr2(axis2)
c**            ds0=dr2(axis2p1)
c**            cLapX0 = xLapCoeff4a(js,dr0,ds0)
c**            cLapX1 = xLapCoeff4b(js,dr0,ds0)
c**  
c**            dLapY0 = yLapCoeff4a(js,dr0,ds0)
c**            dLapY1 = yLapCoeff4b(js,dr0,ds0)
c**  
c**  
c**            ! 4  [ (u.xx+u.yy).x + (v.xx+v.yy).y ] = 0
c**            if( debug.gt.3 )then
c**            aa8(4,0,0,nn)= ( is*rx1*2.*dx122(axis1)*dx112(axis1)+is*rx1*2.*dx122(1)/(2.*dx1(0)))
c**            aa8(4,1,0,nn)= ( is*ry1*2.*dx122(axis1)*dx112(axis1)+is*ry1*2.*dx122(0)/(2.*dx1(1)))
c**            aa8(4,4,0,nn)= (-is*rx1   *dx122(axis1)*dx112(axis1) )  
c**            aa8(4,5,0,nn)= (-is*ry1   *dx122(axis1)*dx112(axis1))
c**              write(*,'(" 4th: xlap4: aLapX0: rect=",e12.4," curv=",e12.4)') aLapX0,aa8(4,0,0,nn)
c**              write(*,'(" 4th: xlap4: aLapX1: rect=",e12.4," curv=",e12.4)') aLapX1,aa8(4,4,0,nn)
c**              write(*,'(" 4th: ylap4: bLapY0: rect=",e12.4," curv=",e12.4)') bLapY0,aa8(4,1,0,nn)
c**              write(*,'(" 4th: ylap4: bLapY1: rect=",e12.4," curv=",e12.4)') bLapY1,aa8(4,5,0,nn)
c**            end if
c**  
c**            aa8(4,0,0,nn)= aLapX0
c**            aa8(4,1,0,nn)= bLapY0
c**            aa8(4,4,0,nn)= aLapX1
c**            aa8(4,5,0,nn)= bLapY1
c**  
c**            aa8(4,2,0,nn)=-cLapX0
c**            aa8(4,3,0,nn)=-dLapY0
c**            aa8(4,6,0,nn)=-cLapX1
c**            aa8(4,7,0,nn)=-dLapY1
c**  
c**            ! 5  [ {(Delta v).x - (Delta u).y}/eps ] =0  -> [ {(v.xxx+v.xyy)-(u.xxy+u.yyy)}/eps ] = 0
c**  
c**            aa8(5,0,0,nn)=-bLapY0/eps1
c**            aa8(5,1,0,nn)= aLapX0/eps1
c**            aa8(5,4,0,nn)=-bLapY1/eps1
c**            aa8(5,5,0,nn)= aLapX1/eps1
c**  
c**            aa8(5,2,0,nn)= dLapY0/eps2
c**            aa8(5,3,0,nn)=-cLapX0/eps2
c**            aa8(5,6,0,nn)= dLapY1/eps2
c**            aa8(5,7,0,nn)=-cLapX1/eps2
c**  
c**  
c**             ! 6  [ n.Delta^2 u/eps ] = 0
c**  
c**             ! assign rx,ry,rxx,rxy,... 
c**             setJacobian(rv1,axis1,axis1p1)
c**             dr0=dr1(axis1)
c**             ds0=dr1(axis1p1)
c**             aLapSq0 = lapSqCoeff4a(is,dr0,ds0)
c**             aLapSq1 = lapSqCoeff4b(is,dr0,ds0)
c**  
c**             if( debug.gt.3 )then
c**               aa8(6,0,0,nn) = -(4./(dx1(axis1)**4) +4./(dx1(0)**2*dx1(1)**2) )
c**               aa8(6,4,0,nn) =   1./(dx1(axis1)**4)
c**               write(*,'(" 4th: lapSq: aLapSq0: rect=",e12.4," curv=",e12.4)') aLapSq0,aa8(6,0,0,nn)
c**               write(*,'(" 4th: lapSq: aLapSq1: rect=",e12.4," curv=",e12.4)') aLapSq1,aa8(6,4,0,nn)
c**             end if
c**  
c**             aa8(6,0,0,nn) = an1*aLapSq0/eps1
c**             aa8(6,1,0,nn) = an2*aLapSq0/eps1
c**             aa8(6,4,0,nn) = an1*aLapSq1/eps1
c**             aa8(6,5,0,nn) = an2*aLapSq1/eps1
c**  
c**             setJacobian(rv2,axis2,axis2p1)
c**             dr0=dr2(axis2)
c**             ds0=dr2(axis2p1)
c**             bLapSq0 = lapSqCoeff4a(js,dr0,ds0)
c**             bLapSq1 = lapSqCoeff4b(js,dr0,ds0)
c**  
c**             aa8(6,2,0,nn) = -an1*bLapSq0/eps2
c**             aa8(6,3,0,nn) = -an2*bLapSq0/eps2
c**             aa8(6,6,0,nn) = -an1*bLapSq1/eps2
c**             aa8(6,7,0,nn) = -an2*bLapSq1/eps2
c**  
c**             ! 7  [ tau.Delta^2 v/eps^2 ] = 0 
c**             aa8(7,0,0,nn) = tau1*aLapSq0/eps1**2
c**             aa8(7,1,0,nn) = tau2*aLapSq0/eps1**2
c**             aa8(7,4,0,nn) = tau1*aLapSq1/eps1**2
c**             aa8(7,5,0,nn) = tau2*aLapSq1/eps1**2
c**  
c**             aa8(7,2,0,nn) = -tau1*bLapSq0/eps2**2
c**             aa8(7,3,0,nn) = -tau2*bLapSq0/eps2**2
c**             aa8(7,6,0,nn) = -tau1*bLapSq1/eps2**2
c**             aa8(7,7,0,nn) = -tau2*bLapSq1/eps2**2
c**  
c**             ! save a copy of the matrix
c**             do n2=0,7
c**             do n1=0,7
c**               aa8(n1,n2,1,nn)=aa8(n1,n2,0,nn)
c**             end do
c**             end do
c**  
c**             ! solve A Q = F
c**             ! factor the matrix
c**             numberOfEquations=8
c**             call dgeco( aa8(0,0,0,nn), numberOfEquations, numberOfEquations, ipvt8(0,nn),rcond,work(0))
c**
c**             if( debug.gt.3 ) write(*,'(" --> 4cth: i1,i2=",2i4," rcond=",e10.2)') i1,i2,rcond
c**             ! '
c**           end if
c**
c**
c**           q(0) = u1(i1-is1,i2-is2,i3,ex)
c**           q(1) = u1(i1-is1,i2-is2,i3,ey)
c**           q(2) = u2(j1-js1,j2-js2,j3,ex)
c**           q(3) = u2(j1-js1,j2-js2,j3,ey)
c**
c**           q(4) = u1(i1-2*is1,i2-2*is2,i3,ex)
c**           q(5) = u1(i1-2*is1,i2-2*is2,i3,ey)
c**           q(6) = u2(j1-2*js1,j2-2*js2,j3,ex)
c**           q(7) = u2(j1-2*js1,j2-2*js2,j3,ey)
c**
c**       if( debug.gt.3 ) write(*,'(" --> 4cth: i1,i2=",2i4," q=",8e10.2)') i1,i2,q(0),q(1),q(2),q(3),q(4),q(5),q(6),q(7)
c**
c**           ! subtract off the contributions from the initial (wrong) values at the ghost points:
c**           do n=0,7
c**             f(n) = (aa8(n,0,1,nn)*q(0)+aa8(n,1,1,nn)*q(1)+aa8(n,2,1,nn)*q(2)+aa8(n,3,1,nn)*q(3)+\
c**                     aa8(n,4,1,nn)*q(4)+aa8(n,5,1,nn)*q(5)+aa8(n,6,1,nn)*q(6)+aa8(n,7,1,nn)*q(7)) - f(n)
c**           end do
c**
c**                                ! '
c**
c**           ! solve A Q = F
c**           job=0
c**           numberOfEquations=8
c**           call dgesl( aa8(0,0,0,nn), numberOfEquations, numberOfEquations, ipvt8(0,nn), f(0), job)
c**
c**       if( debug.gt.3 ) write(*,'(" --> 4cth: i1,i2=",2i4," f(solve)=",8e10.2)') i1,i2,f(0),f(1),f(2),f(3),f(4),f(5),f(6),f(7)
c**           ! '
c**
c**           if( .true. )then
c**           u1(i1-is1,i2-is2,i3,ex)=f(0)
c**           u1(i1-is1,i2-is2,i3,ey)=f(1)
c**           u2(j1-js1,j2-js2,j3,ex)=f(2)
c**           u2(j1-js1,j2-js2,j3,ey)=f(3)
c**
c**           u1(i1-2*is1,i2-2*is2,i3,ex)=f(4)
c**           u1(i1-2*is1,i2-2*is2,i3,ey)=f(5)
c**           u2(j1-2*js1,j2-2*js2,j3,ex)=f(6)
c**           u2(j1-2*js1,j2-2*js2,j3,ey)=f(7)
c**           end if
c**
c**          if( debug.gt.3 )then ! re-evaluate
c**
c**           ! compute the maximum change in the solution for this iteration
c**           do n=0,7
c**             err=max(err,abs(q(n)-f(n)))
c**           end do
c**
c**           ulap1=u1Laplacian42(i1,i2,i3,ex)
c**           vlap1=u1Laplacian42(i1,i2,i3,ey)
c**           ulap2=u2Laplacian42(j1,j2,j3,ex)
c**           vlap2=u2Laplacian42(j1,j2,j3,ey)
c**
c**           ulapSq1=u1LapSq22(i1,i2,i3,ex)
c**           vlapSq1=u1LapSq22(i1,i2,i3,ey)
c**           ulapSq2=u2LapSq22(j1,j2,j3,ex)
c**           vlapSq2=u2LapSq22(j1,j2,j3,ey)
c**
c**           f(0)=(u1x42(i1,i2,i3,ex)+u1y42(i1,i2,i3,ey)) - \
c**                (u2x42(j1,j2,j3,ex)+u2y42(j1,j2,j3,ey))
c**
c**           f(1)=(an1*ulap1+an2*vlap1) - \
c**                (an1*ulap2+an2*vlap2)
c**
c**           f(2)=(u1x42(i1,i2,i3,ey)-u1y42(i1,i2,i3,ex)) - \
c**                (u2x42(j1,j2,j3,ey)-u2y42(j1,j2,j3,ex))
c**           
c**           f(3)=(tau1*ulap1+tau2*vlap1)/eps1 - \
c**                (tau1*ulap2+tau2*vlap2)/eps2
c**    
c**           ! These next we can do to 2nd order -- these need a value on the first ghost line --
c**           f(4)=(u1xxx22(i1,i2,i3,ex)+u1xyy22(i1,i2,i3,ex)+u1xxy22(i1,i2,i3,ey)+u1yyy22(i1,i2,i3,ey)) - \
c**                (u2xxx22(j1,j2,j3,ex)+u2xyy22(j1,j2,j3,ex)+u2xxy22(j1,j2,j3,ey)+u2yyy22(j1,j2,j3,ey))
c**
c**           f(5)=((u1xxx22(i1,i2,i3,ey)+u1xyy22(i1,i2,i3,ey))-(u1xxy22(i1,i2,i3,ex)+u1yyy22(i1,i2,i3,ex)))/eps1 - \
c**                ((u2xxx22(j1,j2,j3,ey)+u2xyy22(j1,j2,j3,ey))-(u2xxy22(j1,j2,j3,ex)+u2yyy22(j1,j2,j3,ex)))/eps2
c**
c**           f(6)=(an1*ulapSq1+an2*vlapSq1)/eps1 - \
c**                (an1*ulapSq2+an2*vlapSq2)/eps2
c**
c**           f(7)=(tau1*ulapSq1+tau2*vlapSq1)/eps1**2 - \
c**                (tau1*ulapSq2+tau2*vlapSq2)/eps2**2
c**
c**    
c**           if( debug.gt.3 ) write(*,'(" --> 4cth: i1,i2=",2i4," f(re-eval)=",8e10.2)') i1,i2,f(0),f(1),f(2),f(3),f(4),f(5),f(6),f(7)
c**             ! '
c**          end if
c**
c**           ! ******************************************************
c**           ! solve for Hz
c**           !  [ w.n/eps ] = 0
c**           !  [ lap(w)/eps ] = 0
c**           !  [ lap(w).n/eps**2 ] = 0
c**           !  [ lapSq(w)/eps**2 ] = 0
c**
c**           ! first evaluate the equations we want to solve with the wrong values at the ghost points:
c**           wlap1=u1Laplacian42(i1,i2,i3,hz)
c**           wlap2=u2Laplacian42(j1,j2,j3,hz)
c**
c**           wlapSq1=u1LapSq22(i1,i2,i3,hz)
c**           wlapSq2=u2LapSq22(j1,j2,j3,hz)
c**
c**           f(0)=(an1*u1x42(i1,i2,i3,hz)+an2*u1y42(i1,i2,i3,hz))/eps1 - \
c**                (an1*u2x42(j1,j2,j3,hz)+an2*u2y42(j1,j2,j3,hz))/eps2
c**
c**           f(1)=wlap1/eps1 - \
c**                wlap2/eps2
c**
c**           ! These next we can do to 2nd order -- these need a value on the first ghost line --
c**           f(2)=(an1*(u1xxx22(i1,i2,i3,hz)+u1xyy22(i1,i2,i3,hz))+an2*(u1xxy22(i1,i2,i3,hz)+u1yyy22(i1,i2,i3,hz)))/eps1**2 - \
c**                (an1*(u2xxx22(j1,j2,j3,hz)+u2xyy22(j1,j2,j3,hz))+an2*(u2xxy22(j1,j2,j3,hz)+u2yyy22(j1,j2,j3,hz)))/eps2**2
c**
c**           f(3)=wlapSq1/eps1**2 - \
c**                wlapSq2/eps2**2
c**
c**           if( .false. .or. (initialized.eq.0 .and. it.eq.1) )then
c**             ! form the matrix for computing Hz (and save factor for later use)
c**
c**             ! 1: [ w.n/eps ] = 0
c**             a0 = (an1*rsxy1(i1,i2,i3,axis1,0)+an2*rsxy1(i1,i2,i3,axis1,1))*dr114(axis1)/eps1
c**             b0 = (an1*rsxy2(j1,j2,j3,axis2,0)+an2*rsxy2(j1,j2,j3,axis2,1))*dr214(axis2)/eps2
c**             aa4(0,0,0,nn) = -is*8.*a0
c**             aa4(0,2,0,nn) =  is*   a0
c**             aa4(0,1,0,nn) =  js*8.*b0
c**             aa4(0,3,0,nn) = -js*   b0
c**  
c**             ! 2: [ lap(w)/eps ] = 0 
c**             aa4(1,0,0,nn) = aLap0/eps1
c**             aa4(1,2,0,nn) = aLap1/eps1
c**             aa4(1,1,0,nn) =-bLap0/eps2
c**             aa4(1,3,0,nn) =-bLap1/eps2
c**  
c**             ! 3  [ (an1*(w.xx+w.yy).x + an2.(w.xx+w.yy).y)/eps**2 ] = 0
c**             aa4(2,0,0,nn)= (an1*aLapX0+an2*bLapY0)/eps1**2
c**             aa4(2,2,0,nn)= (an1*aLapX1+an2*bLapY1)/eps1**2
c**             aa4(2,1,0,nn)=-(an1*cLapX0+an2*dLapY0)/eps2**2
c**             aa4(2,3,0,nn)=-(an1*cLapX1+an2*dLapY1)/eps2**2
c**  
c**             ! 4 [ lapSq(w)/eps**2 ] = 0 
c**             aa4(3,0,0,nn) = aLapSq0/eps1**2
c**             aa4(3,2,0,nn) = aLapSq1/eps1**2
c**             aa4(3,1,0,nn) =-bLapSq0/eps2**2
c**             aa4(3,3,0,nn) =-bLapSq1/eps2**2
c**
c**             ! save a copy of the matrix
c**             do n2=0,3
c**             do n1=0,3
c**               aa4(n1,n2,1,nn)=aa4(n1,n2,0,nn)
c**             end do
c**             end do
c**  
c**             ! factor the matrix
c**             numberOfEquations=4
c**             call dgeco( aa4(0,0,0,nn), numberOfEquations, numberOfEquations, ipvt4(0,nn),rcond,work(0))
c**           end if
c**
c**           q(0) = u1(i1-is1,i2-is2,i3,hz)
c**           q(1) = u2(j1-js1,j2-js2,j3,hz)
c**           q(2) = u1(i1-2*is1,i2-2*is2,i3,hz)
c**           q(3) = u2(j1-2*js1,j2-2*js2,j3,hz)
c**
c**           ! subtract off the contributions from the wrong values at the ghost points:
c**           do n=0,3
c**             f(n) = (aa4(n,0,1,nn)*q(0)+aa4(n,1,1,nn)*q(1)+aa4(n,2,1,nn)*q(2)+aa4(n,3,1,nn)*q(3)) - f(n)
c**           end do
c**           ! solve
c**           numberOfEquations=4
c**           job=0
c**           call dgesl( aa4(0,0,0,nn), numberOfEquations, numberOfEquations, ipvt4(0,nn), f(0), job)
c**
c**           u1(i1-is1,i2-is2,i3,hz)=f(0)
c**           u2(j1-js1,j2-js2,j3,hz)=f(1)
c**           u1(i1-2*is1,i2-2*is2,i3,hz)=f(2)
c**           u2(j1-2*js1,j2-2*js2,j3,hz)=f(3)
c**
c**          if( debug.gt.3 )then ! re-evaluate
c**
c**
c**           wlap1=u1Laplacian42(i1,i2,i3,hz)
c**           wlap2=u2Laplacian42(j1,j2,j3,hz)
c**
c**           wlapSq1=u1LapSq22(i1,i2,i3,hz)
c**           wlapSq2=u2LapSq22(j1,j2,j3,hz)
c**
c**           ! first evaluate the equations we want to solve with the wrong values at the ghost points:
c**           f(0)=(an1*u1x42(i1,i2,i3,hz)+an2*u1y42(i1,i2,i3,hz))/eps1 - \
c**                (an1*u2x42(j1,j2,j3,hz)+an2*u2y42(j1,j2,j3,hz))/eps2
c**
c**           f(1)=wlap1/eps1 - \
c**                wlap2/eps2
c**
c**           ! These next we can do to 2nd order -- these need a value on the first ghost line --
c**           f(2)=(an1*(u1xxx22(i1,i2,i3,hz)+u1xyy22(i1,i2,i3,hz))+an2*(u1xxy22(i1,i2,i3,hz)+u1yyy22(i1,i2,i3,hz)))/eps1**2 - \
c**                (an1*(u2xxx22(j1,j2,j3,hz)+u2xyy22(j1,j2,j3,hz))+an2*(u2xxy22(j1,j2,j3,hz)+u2yyy22(j1,j2,j3,hz)))/eps2**2
c**
c**           f(3)=wlapSq1/eps1**2 - \
c**                wlapSq2/eps2**2
c**    
c**           if( debug.gt.3 ) write(*,'(" --> 4cth: i1,i2=",2i4," hz-f(re-eval)=",8e10.2)') i1,i2,f(0),f(1),f(2),f(3)
c**             ! '
c**          end if
c**
c**
c**
c**           ! ***********************
c**
c**           ! u1(i1-is1,i2-is2,i3,hz)=u2(j1+js1,j2+js2,j3,hz) 
c**           ! u2(j1-js1,j2-js2,j3,hz)=u1(i1+is1,i2+is2,i3,hz)
c**           ! u1(i1-2*is1,i2-2*is2,i3,hz)=u2(j1+2*js1,j2+2*js2,j3,hz) 
c**           ! u2(j1-2*js1,j2-2*js2,j3,hz)=u1(i1+2*is1,i2+2*is2,i3,hz)
c**
c**
c**
c**         endLoops2d()
c**         ! =============== end loops =======================
c**      
c**         periodicUpdate2d(u1,boundaryCondition1,gridIndexRange1,side1,axis1)
c**         periodicUpdate2d(u2,boundaryCondition2,gridIndexRange2,side2,axis2)
c**
c**           if( debug.gt.3 )then 
c**             write(*,'(" ***it=",i2," max-diff = ",e11.2)') it,err
c**               ! '
c**           end if
c**         end do ! ************** end iteration **************
c**
c**
c**         ! now make sure that div(u)=0 etc.
c**         if( .false. )then
c**         beginLoops2d() ! =============== start loops =======================
c**
c**           ! 0  [ u.x + v.y ] = 0
c**c           a8(0,0) = -is*8.*rsxy1(i1,i2,i3,axis1,0)*dr114(axis1)     ! coeff of u1(-1) from [u.x+v.y] 
c**c           a8(0,1) = -is*8.*rsxy1(i1,i2,i3,axis1,1)*dr114(axis1)     ! coeff of v1(-1) from [u.x+v.y] 
c**c           a8(0,4) =  is*   rsxy1(i1,i2,i3,axis1,0)*dr114(axis1)     ! u1(-2)
c**c           a8(0,5) =  is*   rsxy1(i1,i2,i3,axis1,1)*dr114(axis1)     ! v1(-2) 
c**
c**c           a8(0,2) =  js*8.*rsxy2(j1,j2,j3,axis2,0)*dr214(axis2)     ! coeff of u2(-1) from [u.x+v.y] 
c**c           a8(0,3) =  js*8.*rsxy2(j1,j2,j3,axis2,1)*dr214(axis2)  
c**c           a8(0,6) = -js*   rsxy2(j1,j2,j3,axis2,0)*dr214(axis2) 
c**c           a8(0,7) = -js*   rsxy2(j1,j2,j3,axis2,1)*dr214(axis2) 
c**
c**           ! first evaluate the equations we want to solve with the wrong values at the ghost points:
c**           divu=u1x42(i1,i2,i3,ex)+u1y42(i1,i2,i3,ey)
c**           a0=is*   rsxy1(i1,i2,i3,axis1,0)*dr114(axis1)
c**           a1=is*   rsxy1(i1,i2,i3,axis1,1)*dr114(axis1)
c**           aNormSq=a0**2+a1**2
c**           ! now project:  a.uNew = a.uOld - div  ->  (div-a.uOld)+a.uNew = div(uNew) = 0
c**           u1(i1-2*is1,i2-2*is2,i3,ex)=u1(i1-2*is1,i2-2*is2,i3,ex)-divu*a0/aNormSq
c**           u1(i1-2*is1,i2-2*is2,i3,ey)=u1(i1-2*is1,i2-2*is2,i3,ey)-divu*a1/aNormSq
c**
c**           divu=u2x42(j1,j2,j3,ex)+u2y42(j1,j2,j3,ey)
c**           a0=js*   rsxy2(j1,j2,j3,axis2,0)*dr214(axis2) 
c**           a1=js*   rsxy2(j1,j2,j3,axis2,1)*dr214(axis2) 
c**           aNormSq=a0**2+a1**2
c**
c**           u2(j1-2*js1,j2-2*js2,j3,ex)=u2(j1-2*js1,j2-2*js2,j3,ex)-divu*a0/aNormSq
c**           u2(j1-2*js1,j2-2*js2,j3,ey)=u2(j1-2*js1,j2-2*js2,j3,ey)-divu*a1/aNormSq
c**
c**           if( debug.gt.3 )then
c**             divu=u1x42(i1,i2,i3,ex)+u1y42(i1,i2,i3,ey)
c**              write(*,'(" --> 4cth: eval div1,div2=",2e10.2)') u1x42(i1,i2,i3,ex)+u1y42(i1,i2,i3,ey),u2x42(j1,j2,j3,ex)+u2y42(j1,j2,j3,ey)
c**           end if
c**         endLoops2d()
c**         end if
c**       else
c**         stop 3214
c**       end if
c**
c**      else  
c**         ! 3D
c**        stop 6676
c**      end if

      return
      end
