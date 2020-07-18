! This file automatically generated from nrbcUtil.bf with bpp.
! *******************************************************************************
!   Non-reflecting Boundary Condition Utility functions
! *******************************************************************************

! Here are macros that define the planeWave solution
! -*- mode: f90; -*-

! **************************************************
! Here are macros that define the:
!      planeWave solution 
! **************************************************

! ======================================================================
!  Slow start function 
!    tba = length of slow start interval (<0 mean no slow start)
! ======================================================================

! cubic ramp
! tba=max(REAL_EPSILON,tb-ta);
! dta=t-ta;
      
! This (cubic) ramp has 1-derivative zero at t=0 and t=tba

! This ramp has 3-derivatives zero at t=0 and t=1
! This is from ramp.maple
! r=-84*t**5+35*t**4-20*t**7+70*t**6
! rt=-420*t**4+140*t**3-140*t**6+420*t**5
! rtt=-1680*t**3+420*t**2-840*t**5+2100*t**4
! rttt=-5040*t**2+840*t-4200*t**4+8400*t**3


! This ramp has 4-derivatives zero at t=0 and t=1
! This is from ramp.maple
! r=126*(t)**5-315*(t)**8+70*(t)**9-420*(t)**6+540*(t)**7
! rt=630*(t)**4-2520*(t)**7+630*(t)**8-2520*(t)**5+3780*(t)**6
! rtt=2520*(t)**3-17640*(t)**6+5040*(t)**7-12600*(t)**4+22680*(t)**5
! rttt=7560*(t)**2-105840*(t)**5+35280*(t)**6-50400*(t)**3+113400*(t)**4


! ============================================================
!  Initialize parameters for the boundary forcing
!   tba: slow start time interval -- no slow start if this is negative
! ===========================================================

! **************** Here is the new generic plane wave solution *******************

! component n=ex,ey,ez, hx,hy,hz (assumes ex=0)
! one time derivative:
! two time derivatives:
! three time derivatives:

! *************** Here is the 2D planeWave solution ******************************


! one time derivative:

! two time derivatives:

! three time derivatives:

! four time derivatives:

! Here are the slow start versions

! one time derivative:

! two time derivatives:

! three time derivatives:

! four time derivatives:


! **************** Here is the 3D planeWave solution ***************************************



! one time derivative:


! two time derivatives:


! three time derivatives:


! four time derivatives:


! Here are the slow start versions


! one time derivative:


! two time derivatives:

! three time derivatives:

! four time derivatives:


! -------------------------------------------------------------------
! Helper function: Return minus the second time derivative
! -------------------------------------------------------------------


! --------------------------------------------------------------------
! Evaluate the plane wave in 2D
! 
!  x,y,t (input) : point to evaluate at 
!  numberOfTimeDerivatives : evaluate this time derivative
!  ubc(.)  (output) : ubc(ex), etc. 
! --------------------------------------------------------------------


! --------------------------------------------------------------------
! Evaluate the plane wave in 3D
! 
!  x,y,z,t (input) : point to evaluate at 
!  numberOfTimeDerivatives : evaluate this time derivative
!  ubc(.)  (output) : ubc(ex), etc. 
! --------------------------------------------------------------------




! ======================================================================
! Loop over faces of the grid and assign points near that boundary
! ======================================================================


! These next formulae must match the one in getInitialConditions.bC


!===============================================================================================
! Macro:
!   Define a Gaussian Plane Wave incident field in 2D   
!===============================================================================================

!===============================================================================================
! Macro:
!   Define a Gaussian Plane Wave incident field in 3D   
!   Changed: May 16, 2020 -- use pwc() coefficients 
!===============================================================================================

! ===================================================================================
! --- Subtract/add the incident wave, before/after applying the non-reflecting BC ---
! OP : "+" or "-" 
! ADJUST : YES or NO to adjust for the bounding box or not
! ===================================================================================


      subroutine adjustForIncident( nd, nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,
     & gridIndexRange, um, u, un, mask,rsxy, xy,icBoundingBox,bc, 
     & boundaryCondition, ipar, rpar, ierr )
! ===================================================================================
!  Non-reflecting BC utility routine: 
!      Subtract/add an incident field near boundaries before/after a non-reflecting BC is applied.
!
!  gridType : 0=rectangular, 1=curvilinear
!  useForcing : 1=use f for RHS to BC
!  side,axis : 0:1 and 0:2
!
!  um : solution at time t-2*dt
!  u : solution at time t-dt
!  un : solution at time t
!  icBoundingBox : i we are given an initial condition bounding box (with positive volume) then only adjust points in this box
!
! ===================================================================================

      implicit none

      integer nd, nd1a,nd1b,nd2a,nd2b,nd3a,nd3b, ndf1a,ndf1b,ndf2a,
     & ndf2b,ndf3a,ndf3b,n1a,n1b,n2a,n2b,n3a,n3b, ndc, bc,ierr

      real um(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
      real u(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
      real un(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
      integer mask(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
      real rsxy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1,0:nd-1)
      real xy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1)
      integer gridIndexRange(0:1,0:2)
      real icBoundingBox(0:1,0:2)

      integer ipar(0:*),boundaryCondition(0:1,0:2)
      real rpar(0:*),pwc(0:5)

!     --- local variables ----

      integer side,axis,gridType,orderOfAccuracy,orderOfExtrapolation,
     & useForcing,ex,ey,ez,hx,hy,hz,useWhereMask,grid,debug,side1,
     & side2,side3,solveForAllFields
      real dx(0:2),dr(0:2),xa(0:2)
      real t,ep,dt,c
      real dxa,dya,dza
      integer axisp1,axisp2,i1,i2,i3,is1,is2,is3,js1,js2,js3,ks1,ks2,
     & ks3,is
      integer extra,extra1a,extra1b,extra2a,extra2b,extra3a,extra3b,
     & numberOfGhostPoints
      integer edgeDirection,sidea,sideb,sidec,bc1,bc2,numberLinesForPML

      real p0,p2,q0,q2,c1abcem2,c2abcem2
      real an1,an2,an3,aNorm,epsX

      real rx0,ry0,rz0 , rxx0,ryy0, rzz0
      real dr0,cxt,cxx,cyy,czz,cm1,g,bxx,byy,bzz
      real rxNorm, rxNormSq, Dn2, Lu, ur0,urr0, unr0, unrr0
      real ux0,uy0,uz0, uxx0,uyy0,uzz0
      real unx0,uny0,unz0, unxx0,unyy0,unzz0
      real t0,t1,t2

      real eps,mu,kx,ky,kz,slowStartInterval,twoPi,cc

      integer range(0:1,0:2), gid(0:1,0:2), dim(0:1,0:2), halfWidth, s,
     & a
      integer adjustForIncidentField,adjustThreeLevels
      logical adjustForBoundingBox
      real x,y,z

      ! parameters for tanh() in smooth transition for IC bounding box:
      real amp, beta, nv(0:2), xv0(0:2)
      integer smoothBoundingBox

      ! There are different possible incident fields
      integer planeWaveIncidentField, gaussianPlaneWaveIncidentField
      parameter( planeWaveIncidentField=0, 
     & gaussianPlaneWaveIncidentField=1 )
      integer incidentFieldType

      real xi,x0GP,y0GP,z0GP,uex,uey,uez,uhz,betaGP,expxi

      ! boundary conditions parameters
! define BC parameters for fortran routines
! boundary conditions
      integer dirichlet,perfectElectricalConductor,
     & perfectMagneticConductor,planeWaveBoundaryCondition,
     & interfaceBC,symmetryBoundaryCondition,abcEM2,abcPML,abc3,abc4,
     & abc5,rbcNonLocal,rbcLocal,characteristic,absorbing,lastBC
      parameter( dirichlet=1,perfectElectricalConductor=2,
     & perfectMagneticConductor=3,planeWaveBoundaryCondition=4,
     & symmetryBoundaryCondition=5,interfaceBC=6,abcEM2=7,abcPML=8,
     & abc3=9,abc4=10,abc5=11,rbcNonLocal=12,rbcLocal=13,
     & characteristic=14,absorbing=15,lastBC=15 )

      integer rectangular,curvilinear
      parameter(rectangular=0,curvilinear=1)


!$$$c     --- start statement function ----
!$$$      integer kd,m,n
!$$$      real rx,ry,rz,sx,sy,sz,tx,ty,tz
!$$$      ! include 'declareDiffOrder2f.h'
!$$$      ! include 'declareDiffOrder4f.h'
!$$$c*      declareDifferenceOrder2(u,RX)
!$$$c*      declareDifferenceOrder4(u,RX)
!$$$      declareDifferenceOrder2(u,RX)
!$$$      declareDifferenceOrder2(un,none)
!$$$
!$$$      declareDifferenceOrder4(u,RX)
!$$$#Include "declareJacobianDerivatives.h"
!$$$
!$$$c.......statement functions for jacobian
!$$$      rx(i1,i2,i3)=rsxy(i1,i2,i3,0,0)
!$$$      ry(i1,i2,i3)=rsxy(i1,i2,i3,0,1)
!$$$      rz(i1,i2,i3)=rsxy(i1,i2,i3,0,2)
!$$$      sx(i1,i2,i3)=rsxy(i1,i2,i3,1,0)
!$$$      sy(i1,i2,i3)=rsxy(i1,i2,i3,1,1)
!$$$      sz(i1,i2,i3)=rsxy(i1,i2,i3,1,2)
!$$$      tx(i1,i2,i3)=rsxy(i1,i2,i3,2,0)
!$$$      ty(i1,i2,i3)=rsxy(i1,i2,i3,2,1)
!$$$      tz(i1,i2,i3)=rsxy(i1,i2,i3,2,2)
!$$$
!$$$
!$$$c     The next macro call will define the difference approximation statement functions
!$$$      defineDifferenceOrder2Components1(u,RX)
!$$$      defineDifferenceOrder2Components1(un,none)
!$$$      defineDifferenceOrder4Components1(u,RX)
!$$$
!$$$#Include "jacobianDerivatives.h"

!............... end statement functions

      ierr=0

      side                 =ipar(0)
      axis                 =ipar(1)
      n1a                  =ipar(2)
      n1b                  =ipar(3)
      n2a                  =ipar(4)
      n2b                  =ipar(5)
      n3a                  =ipar(6)
      n3b                  =ipar(7)
      gridType             =ipar(8)
      orderOfAccuracy      =ipar(9)
      orderOfExtrapolation =ipar(10)
      useForcing           =ipar(11)
      ex                   =ipar(12)
      ey                   =ipar(13)
      ez                   =ipar(14)
      hx                   =ipar(15)
      hy                   =ipar(16)
      hz                   =ipar(17)
      useWhereMask         =ipar(18)
      grid                 =ipar(19)
      debug                =ipar(20)

      adjustForIncidentField=ipar(25)
      numberLinesForPML    =ipar(26)
      adjustThreeLevels    =ipar(27)

      halfWidth            =ipar(31) ! *new* June 20, 2016 *wdh*
      smoothBoundingBox    =ipar(35) ! smooth the IC at the bounding box edge

      incidentFieldType    =ipar(37)  ! defines type of incident field -- plane wave, Gaussian pulse, ...

      solveForAllFields    =ipar(38)

      dx(0)                =rpar(0)
      dx(1)                =rpar(1)
      dx(2)                =rpar(2)
      dr(0)                =rpar(3)
      dr(1)                =rpar(4)
      dr(2)                =rpar(5)
      t                    =rpar(6)
      ep                   =rpar(7) ! pointer for exact solution
      dt                   =rpar(8)
      c                    =rpar(9)

      eps                  =rpar(10)
      mu                   =rpar(11)
      kx                   =rpar(12)  ! for plane wave forcing
      ky                   =rpar(13)
      kz                   =rpar(14)
      slowStartInterval    =rpar(15)

      pwc(0)               =rpar(20) ! coeffs. for plane wave
      pwc(1)               =rpar(21)
      pwc(2)               =rpar(22)
      pwc(3)               =rpar(23)
      pwc(4)               =rpar(24)
      pwc(5)               =rpar(25)

      xa(0)                =rpar(26)  ! for rectangular grids
      xa(1)                =rpar(27)
      xa(2)                =rpar(28)

      ! parameters for tanh() in smooth transition for IC bounding box:
      beta = rpar(29)
      nv(0)  =rpar(30)
      nv(1)  =rpar(31)
      nv(2)  =rpar(32)
      xv0(0) =rpar(33)
      xv0(1) =rpar(34)
      xv0(2) =rpar(35)


      ! If we are given an initial condition bounding box then only adjust points in this box
      if( icBoundingBox(0,0) .lt. icBoundingBox(1,0) )then
        adjustForBoundingBox=.true.
        ! write(*,'(" adjustForIncident: adjustForBoundingBox, grid=",i4," t=",e9.3)') grid,t
      else
        adjustForBoundingBox=.false.
      end if

      ! adjustForBoundingBox=.false.

      ! Gaussian pulse fix me: 
      if( incidentFieldType.eq.gaussianPlaneWaveIncidentField )then
       ! incidentFieldType=1
        x0GP  =rpar(50)
        y0GP  =rpar(51)
        z0GP  =rpar(52)
        betaGP=rpar(53)
        if( t.le.2*dt )then
          write(*,'(" adjustForIncident: 
     & incidentFieldType=gaussianPlaneWaveIncidentField")')
          write(*,'("    GPW: (x0,y0,z0)=(",e10.2,",",e10.2,",",e10.2,
     & ") beta=",e10.2)') x0GP,y0GP,z0GP,betaGP
          write(*,'(" n1a,n1b,n2a,n2b,n3a,n3b=",6i3)') n1a,n1b,n2a,n2b,
     & n3a,n3b
          write(*,'(" solveForAllFields=",i2," adjustForBoundingBox=",
     & l2)') solveForAllFields,adjustForBoundingBox
          write(*,'(" adjustForIncident: icbb=[",e8.2,",",e8.2,"][",
     & e8.2,",",e8.2,"][",e8.2,",",e8.2,"]")') icBoundingBox(0,0),
     & icBoundingBox(1,0),icBoundingBox(0,1),icBoundingBox(1,1),
     & icBoundingBox(0,2),icBoundingBox(1,2)
        end if
      end if

      if( abs(pwc(0))+abs(pwc(1))+abs(pwc(2)) .eq. 0. )then
        ! sanity check
        stop 12345
      end if
      if( debug.gt.1 )then
        write(*,'(" adjustForIncident",i2,": **START** grid=",i4," 
     & side,axis=",2i2)') adjustForIncidentField,grid,side,axis
        ! ' 
      end if


      ! for plane wave forcing 
      twoPi=8.*atan2(1.,1.)
      cc= c*sqrt( kx*kx+ky*ky+kz*kz )

      if( adjustForIncidentField.eq.-1 )then
        ! --- Subtract off the incident wave, before applying the non-reflecting BC ---
        if( adjustForBoundingBox )then
           if( debug.gt.4 )then
             write(*,'(" adjust -: kx,ky,kz,eps,cc",5e10.3)') kx,ky,kz,
     & eps,cc
           end if
           amp=1.
            ! adjust this many points near the boundary (needs to include width of extrapolation too!)
            ! halfWidth = orderOfAccuracy/2 ** now passed in**
            ! restrict all bounds to fit in the array dimensions: 
            dim(0,0)=nd1a
            dim(1,0)=nd1b
            dim(0,1)=nd2a
            dim(1,1)=nd2b
            dim(0,2)=nd3a
            dim(1,2)=nd3b
            ! gid(size,axis) : holds bounds on the points that still may need adjusting 
            !                : as we assign points near faces, gid will be reduced in size so that we 
            !                : do not adjust points multiple times.  
            gid(0,2)=0
            gid(1,2)=0
            range(0,2)=0
            range(1,2)=0
            do axis=0,nd-1
            do side=0,1
             if( boundaryCondition(side,axis).eq.abcPML )then
              gid(side,axis)=dim(side,axis)
             else
              gid(side,axis)=gridIndexRange(side,axis)-halfWidth*(1-2*
     & side)
             end if
            end do
            end do
            do axis=0,nd-1
            do side=0,1
             if( boundaryCondition(side,axis).ge.abcEM2 .and. 
     & boundaryCondition(side,axis).le.rbcLocal )then
              do a=0,nd-1
              do s=0,1
               range(s,a)=gid(s,a)
              end do
              end do
              if( boundaryCondition(side,axis).eq.abcPML )then
                ! we need to adjust more lines for a PML
                range(  side,axis)=dim(side,axis)
                range(1-side,axis)=gridIndexRange(side,axis)+(
     & numberLinesForPML+halfWidth)*(1-2*side) ! check this
              else
                range(0,axis)=gridIndexRange(side,axis)-halfWidth
                range(1,axis)=gridIndexRange(side,axis)+halfWidth
              end if
              ! restrict all bounds to fit in the array dimensions: 
              do a=0,nd-1
              do s=0,1
               range(s,a)=max(dim(0,a),min(dim(1,a), range(s,a) ))
              end do
              end do
              ! reduce gid bounds for future tangential directions so pts are not assigned twice
              ! gid(side,axis)=gridIndexRange(side,axis)+(halfWidth+1)*(1-2*side)
              gid(side,axis)=range(1-side,axis)+(1)*(1-2*side)
              n1a=range(0,0)
              n1b=range(1,0)
              n2a=range(0,1)
              n2b=range(1,1)
              n3a=range(0,2)
              n3b=range(1,2)
              if( debug.gt.0 )then
                write(*,'(" abc:adjustIncindent: side,axis,bc=",3i2," 
     & n1a,n1b,...=",3(i3,1x))') side,axis,boundaryCondition(side,
     & axis),n1a,n1b,n2a,n2b,n3a,n3b
                write(*,'(" abc:adjustIncindent: gid=",4i4)') gid(0,0),
     & gid(1,0),gid(0,1),gid(1,1)
                ! ' 
              end if
             if( nd.eq.2 )then
               do i3=n3a,n3b
               do i2=n2a,n2b
               do i1=n1a,n1b
                 if( gridType.eq.rectangular )then
                   x = xa(0)+i1*dx(0)
                   y = xa(1)+i2*dx(1)
                 else
                   x=xy(i1,i2,i3,0)
                   y=xy(i1,i2,i3,1)
                 end if
                 ! if( debug.gt.1 )then
                 !  t0=planeWave2Dhz0(x,y,t-dt)
                 !  write(*,'("nrbc: adjust -: i=",2i3," Hz,true=",2e10.3)') i1,i2,u(i1,i2,i3,hz),t0
                 ! end if
                   if( x.ge.icBoundingBox(0,0) .and. 
     & x.le.icBoundingBox(1,0) .and. y.ge.icBoundingBox(0,1) .and. 
     & y.le.icBoundingBox(1,1) )then
                     if( smoothBoundingBox.eq.1 )then
                       amp = (.5*(1.-tanh(beta*twoPi*(nv(0)*((x)-xv0(0)
     & )+nv(1)*((y)-xv0(1))-cc*(t)))))
                     else
                       amp=1.
                     end if
                   ! if( amp.gt. 1.e-9 )then
                   !   if( x.gt.1.5 .and. y.gt.-.05 .and. y.lt.0.5 )then
                   !     write(*,'(" x,amp=",2e10.2)') x,amp
                   !   end if
                 if( incidentFieldType .eq. planeWaveIncidentField )
     & then
                   ! --- plane wave incident field ---
                   u(i1,i2,i3,ex) = u(i1,i2,i3,ex) - amp*sin(twoPi*(kx*
     & (x)+ky*(y)-cc*(t-dt)))*pwc(0)
                   u(i1,i2,i3,ey) = u(i1,i2,i3,ey) - amp*sin(twoPi*(kx*
     & (x)+ky*(y)-cc*(t-dt)))*pwc(1)
                   u(i1,i2,i3,hz) = u(i1,i2,i3,hz) - amp*sin(twoPi*(kx*
     & (x)+ky*(y)-cc*(t-dt)))*pwc(5)
                   un(i1,i2,i3,ex)= un(i1,i2,i3,ex) - amp*sin(twoPi*(
     & kx*(x)+ky*(y)-cc*(t)))*pwc(0)
                   un(i1,i2,i3,ey)= un(i1,i2,i3,ey) - amp*sin(twoPi*(
     & kx*(x)+ky*(y)-cc*(t)))*pwc(1)
                   un(i1,i2,i3,hz)= un(i1,i2,i3,hz) - amp*sin(twoPi*(
     & kx*(x)+ky*(y)-cc*(t)))*pwc(5)
                   if( adjustThreeLevels.eq.1 )then
                    um(i1,i2,i3,ex) = um(i1,i2,i3,ex) - amp*sin(twoPi*(
     & kx*(x)+ky*(y)-cc*(t-2.*dt)))*pwc(0)
                    um(i1,i2,i3,ey) = um(i1,i2,i3,ey) - amp*sin(twoPi*(
     & kx*(x)+ky*(y)-cc*(t-2.*dt)))*pwc(1)
                    um(i1,i2,i3,hz) = um(i1,i2,i3,hz) - amp*sin(twoPi*(
     & kx*(x)+ky*(y)-cc*(t-2.*dt)))*pwc(5)
                   end if
                 else if( incidentFieldType .eq. 
     & gaussianPlaneWaveIncidentField )then
                   ! --- plane wave incident field ---
                     xi = kx*(x-x0GP) + ky*(y-y0GP) - cc*(t-dt)
                     expxi = amp*exp(-betaGP*xi*xi )
                     if( solveForAllFields.eq.0 )then
                       ! *new* way June 16, 2020 *wdh*
                       u(i1,i2,i3,ex) = u(i1,i2,i3,ex) - expxi*pwc(0)
                       u(i1,i2,i3,ey) = u(i1,i2,i3,ey) - expxi*pwc(1)
                       u(i1,i2,i3,hz) = u(i1,i2,i3,hz) - expxi*pwc(5)
                       ! *old way*
                       ! uhz = expxi
                       ! uex = uhz*(-ky/(eps*cc))
                       ! uey = uhz*( kx/(eps*cc))
                       ! u(i1,i2,i3,ex) = u(i1,i2,i3,ex) - uex
                       ! u(i1,i2,i3,ey) = u(i1,i2,i3,ey) - uey
                       ! u(i1,i2,i3,hz) = u(i1,i2,i3,hz) - uhz
                     else
                       u(i1,i2,i3,ex) = u(i1,i2,i3,ex) - expxi*pwc(0)
                       u(i1,i2,i3,ey) = u(i1,i2,i3,ey) - expxi*pwc(1)
                       u(i1,i2,i3,ez) = u(i1,i2,i3,ez) - expxi*pwc(2)
                       u(i1,i2,i3,hx) = u(i1,i2,i3,hx) - expxi*pwc(3)
                       u(i1,i2,i3,hy) = u(i1,i2,i3,hy) - expxi*pwc(4)
                       u(i1,i2,i3,hz) = u(i1,i2,i3,hz) - expxi*pwc(5)
                     end if
                     xi = kx*(x-x0GP) + ky*(y-y0GP) - cc*(t)
                     expxi = amp*exp(-betaGP*xi*xi )
                     if( solveForAllFields.eq.0 )then
                       ! *new* way June 16, 2020 *wdh*
                       un(i1,i2,i3,ex) = un(i1,i2,i3,ex) - expxi*pwc(0)
                       un(i1,i2,i3,ey) = un(i1,i2,i3,ey) - expxi*pwc(1)
                       un(i1,i2,i3,hz) = un(i1,i2,i3,hz) - expxi*pwc(5)
                       ! *old way*
                       ! uhz = expxi
                       ! uex = uhz*(-ky/(eps*cc))
                       ! uey = uhz*( kx/(eps*cc))
                       ! un(i1,i2,i3,ex) = un(i1,i2,i3,ex) - uex
                       ! un(i1,i2,i3,ey) = un(i1,i2,i3,ey) - uey
                       ! un(i1,i2,i3,hz) = un(i1,i2,i3,hz) - uhz
                     else
                       un(i1,i2,i3,ex) = un(i1,i2,i3,ex) - expxi*pwc(0)
                       un(i1,i2,i3,ey) = un(i1,i2,i3,ey) - expxi*pwc(1)
                       un(i1,i2,i3,ez) = un(i1,i2,i3,ez) - expxi*pwc(2)
                       un(i1,i2,i3,hx) = un(i1,i2,i3,hx) - expxi*pwc(3)
                       un(i1,i2,i3,hy) = un(i1,i2,i3,hy) - expxi*pwc(4)
                       un(i1,i2,i3,hz) = un(i1,i2,i3,hz) - expxi*pwc(5)
                     end if
                   if( adjustThreeLevels.eq.1 )then
                       xi = kx*(x-x0GP) + ky*(y-y0GP) - cc*(t-2.*dt)
                       expxi = amp*exp(-betaGP*xi*xi )
                       if( solveForAllFields.eq.0 )then
                         ! *new* way June 16, 2020 *wdh*
                         um(i1,i2,i3,ex) = um(i1,i2,i3,ex) - expxi*pwc(
     & 0)
                         um(i1,i2,i3,ey) = um(i1,i2,i3,ey) - expxi*pwc(
     & 1)
                         um(i1,i2,i3,hz) = um(i1,i2,i3,hz) - expxi*pwc(
     & 5)
                         ! *old way*
                         ! uhz = expxi
                         ! uex = uhz*(-ky/(eps*cc))
                         ! uey = uhz*( kx/(eps*cc))
                         ! um(i1,i2,i3,ex) = um(i1,i2,i3,ex) - uex
                         ! um(i1,i2,i3,ey) = um(i1,i2,i3,ey) - uey
                         ! um(i1,i2,i3,hz) = um(i1,i2,i3,hz) - uhz
                       else
                         um(i1,i2,i3,ex) = um(i1,i2,i3,ex) - expxi*pwc(
     & 0)
                         um(i1,i2,i3,ey) = um(i1,i2,i3,ey) - expxi*pwc(
     & 1)
                         um(i1,i2,i3,ez) = um(i1,i2,i3,ez) - expxi*pwc(
     & 2)
                         um(i1,i2,i3,hx) = um(i1,i2,i3,hx) - expxi*pwc(
     & 3)
                         um(i1,i2,i3,hy) = um(i1,i2,i3,hy) - expxi*pwc(
     & 4)
                         um(i1,i2,i3,hz) = um(i1,i2,i3,hz) - expxi*pwc(
     & 5)
                       end if
                   end if
                 else
                   write(*,'("adjust Inc. : ERROR: unknown incident 
     & field type")')
                   stop 88991
                 end if
                  endif
               end do
               end do
               end do
             else
               ! --- THREE DIMENSIONS ---
               !if( incidentFieldType .eq. gaussianPlaneWaveIncidentField )then
               !  write(*,'(" - YES : alter fields for gaussianPlaneWaveIncidentField ")')
               !end if 
               do i3=n3a,n3b
               do i2=n2a,n2b
               do i1=n1a,n1b
                 if( gridType.eq.rectangular )then
                   x = xa(0)+i1*dx(0)
                   y = xa(1)+i2*dx(1)
                   z = xa(2)+i3*dx(2)
                 else
                   x=xy(i1,i2,i3,0)
                   y=xy(i1,i2,i3,1)
                   z=xy(i1,i2,i3,2)
                 end if
                 if( x.ge.icBoundingBox(0,0) .and. x.le.icBoundingBox(
     & 1,0) .and. y.ge.icBoundingBox(0,1) .and. y.le.icBoundingBox(1,
     & 1) .and. z.ge.icBoundingBox(0,2) .and. z.le.icBoundingBox(1,2) 
     & )then
                   if( smoothBoundingBox.eq.1 )then
                     amp = (.5*(1.-tanh(beta*twoPi*(nv(0)*((x)-xv0(0))+
     & nv(1)*((y)-xv0(1))+nv(2)*((z)-xv0(2))-cc*(t)))))
                   else
                     amp=1.
                   end if
                 if( incidentFieldType .eq. planeWaveIncidentField )
     & then
                   u(i1,i2,i3,ex) = u(i1,i2,i3,ex) - amp*sin(twoPi*(kx*
     & (x)+ky*(y)+kz*(z)-cc*(t-dt)))*pwc(0)
                   u(i1,i2,i3,ey) = u(i1,i2,i3,ey) - amp*sin(twoPi*(kx*
     & (x)+ky*(y)+kz*(z)-cc*(t-dt)))*pwc(1)
                   u(i1,i2,i3,ez) = u(i1,i2,i3,ez) - amp*sin(twoPi*(kx*
     & (x)+ky*(y)+kz*(z)-cc*(t-dt)))*pwc(2)
                   un(i1,i2,i3,ex)=un(i1,i2,i3,ex) - amp*sin(twoPi*(kx*
     & (x)+ky*(y)+kz*(z)-cc*(t)))*pwc(0)
                   un(i1,i2,i3,ey)=un(i1,i2,i3,ey) - amp*sin(twoPi*(kx*
     & (x)+ky*(y)+kz*(z)-cc*(t)))*pwc(1)
                   un(i1,i2,i3,ez)=un(i1,i2,i3,ez) - amp*sin(twoPi*(kx*
     & (x)+ky*(y)+kz*(z)-cc*(t)))*pwc(2)
                   if( adjustThreeLevels.eq.1 )then
                    um(i1,i2,i3,ex)=um(i1,i2,i3,ex) - amp*sin(twoPi*(
     & kx*(x)+ky*(y)+kz*(z)-cc*(t-2.*dt)))*pwc(0)
                    um(i1,i2,i3,ey)=um(i1,i2,i3,ey) - amp*sin(twoPi*(
     & kx*(x)+ky*(y)+kz*(z)-cc*(t-2.*dt)))*pwc(1)
                    um(i1,i2,i3,ez)=um(i1,i2,i3,ez) - amp*sin(twoPi*(
     & kx*(x)+ky*(y)+kz*(z)-cc*(t-2.*dt)))*pwc(2)
                   end if
                 else if( incidentFieldType .eq. 
     & gaussianPlaneWaveIncidentField )then
                   ! --- Gaussian plane wave incident field ---  *wdh* Aug 18, 2019
                     xi = kx*(x-x0GP) + ky*(y-y0GP) + kz*(y-z0GP) - cc*
     & (t-dt)
                     expxi = amp*exp(-betaGP*xi*xi )
                     u(i1,i2,i3,ex) = u(i1,i2,i3,ex) - expxi*pwc(0)
                     u(i1,i2,i3,ey) = u(i1,i2,i3,ey) - expxi*pwc(1)
                     u(i1,i2,i3,ez) = u(i1,i2,i3,ez) - expxi*pwc(2)
                     if( solveForAllFields.eq.1 )then
                       u(i1,i2,i3,hx) = u(i1,i2,i3,hx) - expxi*pwc(3)
                       u(i1,i2,i3,hy) = u(i1,i2,i3,hy) - expxi*pwc(4)
                       u(i1,i2,i3,hz) = u(i1,i2,i3,hz) - expxi*pwc(5)
                     end if
                     xi = kx*(x-x0GP) + ky*(y-y0GP) + kz*(y-z0GP) - cc*
     & (t)
                     expxi = amp*exp(-betaGP*xi*xi )
                     un(i1,i2,i3,ex) = un(i1,i2,i3,ex) - expxi*pwc(0)
                     un(i1,i2,i3,ey) = un(i1,i2,i3,ey) - expxi*pwc(1)
                     un(i1,i2,i3,ez) = un(i1,i2,i3,ez) - expxi*pwc(2)
                     if( solveForAllFields.eq.1 )then
                       un(i1,i2,i3,hx) = un(i1,i2,i3,hx) - expxi*pwc(3)
                       un(i1,i2,i3,hy) = un(i1,i2,i3,hy) - expxi*pwc(4)
                       un(i1,i2,i3,hz) = un(i1,i2,i3,hz) - expxi*pwc(5)
                     end if
                   if( adjustThreeLevels.eq.1 )then
                       xi = kx*(x-x0GP) + ky*(y-y0GP) + kz*(y-z0GP) - 
     & cc*(t-2.*dt)
                       expxi = amp*exp(-betaGP*xi*xi )
                       um(i1,i2,i3,ex) = um(i1,i2,i3,ex) - expxi*pwc(0)
                       um(i1,i2,i3,ey) = um(i1,i2,i3,ey) - expxi*pwc(1)
                       um(i1,i2,i3,ez) = um(i1,i2,i3,ez) - expxi*pwc(2)
                       if( solveForAllFields.eq.1 )then
                         um(i1,i2,i3,hx) = um(i1,i2,i3,hx) - expxi*pwc(
     & 3)
                         um(i1,i2,i3,hy) = um(i1,i2,i3,hy) - expxi*pwc(
     & 4)
                         um(i1,i2,i3,hz) = um(i1,i2,i3,hz) - expxi*pwc(
     & 5)
                       end if
                   end if
                 else
                   write(*,'("adjust Inc. : ERROR: unknown incident 
     & field type")')
                   stop 88992
                 end if
                 endif
               end do
               end do
               end do
             end if
             end if
            end do
            end do
        else
           if( debug.gt.4 )then
             write(*,'(" adjust -: kx,ky,kz,eps,cc",5e10.3)') kx,ky,kz,
     & eps,cc
           end if
           amp=1.
            ! adjust this many points near the boundary (needs to include width of extrapolation too!)
            ! halfWidth = orderOfAccuracy/2 ** now passed in**
            ! restrict all bounds to fit in the array dimensions: 
            dim(0,0)=nd1a
            dim(1,0)=nd1b
            dim(0,1)=nd2a
            dim(1,1)=nd2b
            dim(0,2)=nd3a
            dim(1,2)=nd3b
            ! gid(size,axis) : holds bounds on the points that still may need adjusting 
            !                : as we assign points near faces, gid will be reduced in size so that we 
            !                : do not adjust points multiple times.  
            gid(0,2)=0
            gid(1,2)=0
            range(0,2)=0
            range(1,2)=0
            do axis=0,nd-1
            do side=0,1
             if( boundaryCondition(side,axis).eq.abcPML )then
              gid(side,axis)=dim(side,axis)
             else
              gid(side,axis)=gridIndexRange(side,axis)-halfWidth*(1-2*
     & side)
             end if
            end do
            end do
            do axis=0,nd-1
            do side=0,1
             if( boundaryCondition(side,axis).ge.abcEM2 .and. 
     & boundaryCondition(side,axis).le.rbcLocal )then
              do a=0,nd-1
              do s=0,1
               range(s,a)=gid(s,a)
              end do
              end do
              if( boundaryCondition(side,axis).eq.abcPML )then
                ! we need to adjust more lines for a PML
                range(  side,axis)=dim(side,axis)
                range(1-side,axis)=gridIndexRange(side,axis)+(
     & numberLinesForPML+halfWidth)*(1-2*side) ! check this
              else
                range(0,axis)=gridIndexRange(side,axis)-halfWidth
                range(1,axis)=gridIndexRange(side,axis)+halfWidth
              end if
              ! restrict all bounds to fit in the array dimensions: 
              do a=0,nd-1
              do s=0,1
               range(s,a)=max(dim(0,a),min(dim(1,a), range(s,a) ))
              end do
              end do
              ! reduce gid bounds for future tangential directions so pts are not assigned twice
              ! gid(side,axis)=gridIndexRange(side,axis)+(halfWidth+1)*(1-2*side)
              gid(side,axis)=range(1-side,axis)+(1)*(1-2*side)
              n1a=range(0,0)
              n1b=range(1,0)
              n2a=range(0,1)
              n2b=range(1,1)
              n3a=range(0,2)
              n3b=range(1,2)
              if( debug.gt.0 )then
                write(*,'(" abc:adjustIncindent: side,axis,bc=",3i2," 
     & n1a,n1b,...=",3(i3,1x))') side,axis,boundaryCondition(side,
     & axis),n1a,n1b,n2a,n2b,n3a,n3b
                write(*,'(" abc:adjustIncindent: gid=",4i4)') gid(0,0),
     & gid(1,0),gid(0,1),gid(1,1)
                ! ' 
              end if
             if( nd.eq.2 )then
               do i3=n3a,n3b
               do i2=n2a,n2b
               do i1=n1a,n1b
                 if( gridType.eq.rectangular )then
                   x = xa(0)+i1*dx(0)
                   y = xa(1)+i2*dx(1)
                 else
                   x=xy(i1,i2,i3,0)
                   y=xy(i1,i2,i3,1)
                 end if
                 ! if( debug.gt.1 )then
                 !  t0=planeWave2Dhz0(x,y,t-dt)
                 !  write(*,'("nrbc: adjust -: i=",2i3," Hz,true=",2e10.3)') i1,i2,u(i1,i2,i3,hz),t0
                 ! end if
                 if( incidentFieldType .eq. planeWaveIncidentField )
     & then
                   ! --- plane wave incident field ---
                   u(i1,i2,i3,ex) = u(i1,i2,i3,ex) - amp*sin(twoPi*(kx*
     & (x)+ky*(y)-cc*(t-dt)))*pwc(0)
                   u(i1,i2,i3,ey) = u(i1,i2,i3,ey) - amp*sin(twoPi*(kx*
     & (x)+ky*(y)-cc*(t-dt)))*pwc(1)
                   u(i1,i2,i3,hz) = u(i1,i2,i3,hz) - amp*sin(twoPi*(kx*
     & (x)+ky*(y)-cc*(t-dt)))*pwc(5)
                   un(i1,i2,i3,ex)= un(i1,i2,i3,ex) - amp*sin(twoPi*(
     & kx*(x)+ky*(y)-cc*(t)))*pwc(0)
                   un(i1,i2,i3,ey)= un(i1,i2,i3,ey) - amp*sin(twoPi*(
     & kx*(x)+ky*(y)-cc*(t)))*pwc(1)
                   un(i1,i2,i3,hz)= un(i1,i2,i3,hz) - amp*sin(twoPi*(
     & kx*(x)+ky*(y)-cc*(t)))*pwc(5)
                   if( adjustThreeLevels.eq.1 )then
                    um(i1,i2,i3,ex) = um(i1,i2,i3,ex) - amp*sin(twoPi*(
     & kx*(x)+ky*(y)-cc*(t-2.*dt)))*pwc(0)
                    um(i1,i2,i3,ey) = um(i1,i2,i3,ey) - amp*sin(twoPi*(
     & kx*(x)+ky*(y)-cc*(t-2.*dt)))*pwc(1)
                    um(i1,i2,i3,hz) = um(i1,i2,i3,hz) - amp*sin(twoPi*(
     & kx*(x)+ky*(y)-cc*(t-2.*dt)))*pwc(5)
                   end if
                 else if( incidentFieldType .eq. 
     & gaussianPlaneWaveIncidentField )then
                   ! --- plane wave incident field ---
                     xi = kx*(x-x0GP) + ky*(y-y0GP) - cc*(t-dt)
                     expxi = amp*exp(-betaGP*xi*xi )
                     if( solveForAllFields.eq.0 )then
                       ! *new* way June 16, 2020 *wdh*
                       u(i1,i2,i3,ex) = u(i1,i2,i3,ex) - expxi*pwc(0)
                       u(i1,i2,i3,ey) = u(i1,i2,i3,ey) - expxi*pwc(1)
                       u(i1,i2,i3,hz) = u(i1,i2,i3,hz) - expxi*pwc(5)
                       ! *old way*
                       ! uhz = expxi
                       ! uex = uhz*(-ky/(eps*cc))
                       ! uey = uhz*( kx/(eps*cc))
                       ! u(i1,i2,i3,ex) = u(i1,i2,i3,ex) - uex
                       ! u(i1,i2,i3,ey) = u(i1,i2,i3,ey) - uey
                       ! u(i1,i2,i3,hz) = u(i1,i2,i3,hz) - uhz
                     else
                       u(i1,i2,i3,ex) = u(i1,i2,i3,ex) - expxi*pwc(0)
                       u(i1,i2,i3,ey) = u(i1,i2,i3,ey) - expxi*pwc(1)
                       u(i1,i2,i3,ez) = u(i1,i2,i3,ez) - expxi*pwc(2)
                       u(i1,i2,i3,hx) = u(i1,i2,i3,hx) - expxi*pwc(3)
                       u(i1,i2,i3,hy) = u(i1,i2,i3,hy) - expxi*pwc(4)
                       u(i1,i2,i3,hz) = u(i1,i2,i3,hz) - expxi*pwc(5)
                     end if
                     xi = kx*(x-x0GP) + ky*(y-y0GP) - cc*(t)
                     expxi = amp*exp(-betaGP*xi*xi )
                     if( solveForAllFields.eq.0 )then
                       ! *new* way June 16, 2020 *wdh*
                       un(i1,i2,i3,ex) = un(i1,i2,i3,ex) - expxi*pwc(0)
                       un(i1,i2,i3,ey) = un(i1,i2,i3,ey) - expxi*pwc(1)
                       un(i1,i2,i3,hz) = un(i1,i2,i3,hz) - expxi*pwc(5)
                       ! *old way*
                       ! uhz = expxi
                       ! uex = uhz*(-ky/(eps*cc))
                       ! uey = uhz*( kx/(eps*cc))
                       ! un(i1,i2,i3,ex) = un(i1,i2,i3,ex) - uex
                       ! un(i1,i2,i3,ey) = un(i1,i2,i3,ey) - uey
                       ! un(i1,i2,i3,hz) = un(i1,i2,i3,hz) - uhz
                     else
                       un(i1,i2,i3,ex) = un(i1,i2,i3,ex) - expxi*pwc(0)
                       un(i1,i2,i3,ey) = un(i1,i2,i3,ey) - expxi*pwc(1)
                       un(i1,i2,i3,ez) = un(i1,i2,i3,ez) - expxi*pwc(2)
                       un(i1,i2,i3,hx) = un(i1,i2,i3,hx) - expxi*pwc(3)
                       un(i1,i2,i3,hy) = un(i1,i2,i3,hy) - expxi*pwc(4)
                       un(i1,i2,i3,hz) = un(i1,i2,i3,hz) - expxi*pwc(5)
                     end if
                   if( adjustThreeLevels.eq.1 )then
                       xi = kx*(x-x0GP) + ky*(y-y0GP) - cc*(t-2.*dt)
                       expxi = amp*exp(-betaGP*xi*xi )
                       if( solveForAllFields.eq.0 )then
                         ! *new* way June 16, 2020 *wdh*
                         um(i1,i2,i3,ex) = um(i1,i2,i3,ex) - expxi*pwc(
     & 0)
                         um(i1,i2,i3,ey) = um(i1,i2,i3,ey) - expxi*pwc(
     & 1)
                         um(i1,i2,i3,hz) = um(i1,i2,i3,hz) - expxi*pwc(
     & 5)
                         ! *old way*
                         ! uhz = expxi
                         ! uex = uhz*(-ky/(eps*cc))
                         ! uey = uhz*( kx/(eps*cc))
                         ! um(i1,i2,i3,ex) = um(i1,i2,i3,ex) - uex
                         ! um(i1,i2,i3,ey) = um(i1,i2,i3,ey) - uey
                         ! um(i1,i2,i3,hz) = um(i1,i2,i3,hz) - uhz
                       else
                         um(i1,i2,i3,ex) = um(i1,i2,i3,ex) - expxi*pwc(
     & 0)
                         um(i1,i2,i3,ey) = um(i1,i2,i3,ey) - expxi*pwc(
     & 1)
                         um(i1,i2,i3,ez) = um(i1,i2,i3,ez) - expxi*pwc(
     & 2)
                         um(i1,i2,i3,hx) = um(i1,i2,i3,hx) - expxi*pwc(
     & 3)
                         um(i1,i2,i3,hy) = um(i1,i2,i3,hy) - expxi*pwc(
     & 4)
                         um(i1,i2,i3,hz) = um(i1,i2,i3,hz) - expxi*pwc(
     & 5)
                       end if
                   end if
                 else
                   write(*,'("adjust Inc. : ERROR: unknown incident 
     & field type")')
                   stop 88991
                 end if
               end do
               end do
               end do
             else
               ! --- THREE DIMENSIONS ---
               !if( incidentFieldType .eq. gaussianPlaneWaveIncidentField )then
               !  write(*,'(" - NO : alter fields for gaussianPlaneWaveIncidentField ")')
               !end if 
               do i3=n3a,n3b
               do i2=n2a,n2b
               do i1=n1a,n1b
                 if( gridType.eq.rectangular )then
                   x = xa(0)+i1*dx(0)
                   y = xa(1)+i2*dx(1)
                   z = xa(2)+i3*dx(2)
                 else
                   x=xy(i1,i2,i3,0)
                   y=xy(i1,i2,i3,1)
                   z=xy(i1,i2,i3,2)
                 end if
                 if( incidentFieldType .eq. planeWaveIncidentField )
     & then
                   u(i1,i2,i3,ex) = u(i1,i2,i3,ex) - amp*sin(twoPi*(kx*
     & (x)+ky*(y)+kz*(z)-cc*(t-dt)))*pwc(0)
                   u(i1,i2,i3,ey) = u(i1,i2,i3,ey) - amp*sin(twoPi*(kx*
     & (x)+ky*(y)+kz*(z)-cc*(t-dt)))*pwc(1)
                   u(i1,i2,i3,ez) = u(i1,i2,i3,ez) - amp*sin(twoPi*(kx*
     & (x)+ky*(y)+kz*(z)-cc*(t-dt)))*pwc(2)
                   un(i1,i2,i3,ex)=un(i1,i2,i3,ex) - amp*sin(twoPi*(kx*
     & (x)+ky*(y)+kz*(z)-cc*(t)))*pwc(0)
                   un(i1,i2,i3,ey)=un(i1,i2,i3,ey) - amp*sin(twoPi*(kx*
     & (x)+ky*(y)+kz*(z)-cc*(t)))*pwc(1)
                   un(i1,i2,i3,ez)=un(i1,i2,i3,ez) - amp*sin(twoPi*(kx*
     & (x)+ky*(y)+kz*(z)-cc*(t)))*pwc(2)
                   if( adjustThreeLevels.eq.1 )then
                    um(i1,i2,i3,ex)=um(i1,i2,i3,ex) - amp*sin(twoPi*(
     & kx*(x)+ky*(y)+kz*(z)-cc*(t-2.*dt)))*pwc(0)
                    um(i1,i2,i3,ey)=um(i1,i2,i3,ey) - amp*sin(twoPi*(
     & kx*(x)+ky*(y)+kz*(z)-cc*(t-2.*dt)))*pwc(1)
                    um(i1,i2,i3,ez)=um(i1,i2,i3,ez) - amp*sin(twoPi*(
     & kx*(x)+ky*(y)+kz*(z)-cc*(t-2.*dt)))*pwc(2)
                   end if
                 else if( incidentFieldType .eq. 
     & gaussianPlaneWaveIncidentField )then
                   ! --- Gaussian plane wave incident field ---  *wdh* Aug 18, 2019
                     xi = kx*(x-x0GP) + ky*(y-y0GP) + kz*(y-z0GP) - cc*
     & (t-dt)
                     expxi = amp*exp(-betaGP*xi*xi )
                     u(i1,i2,i3,ex) = u(i1,i2,i3,ex) - expxi*pwc(0)
                     u(i1,i2,i3,ey) = u(i1,i2,i3,ey) - expxi*pwc(1)
                     u(i1,i2,i3,ez) = u(i1,i2,i3,ez) - expxi*pwc(2)
                     if( solveForAllFields.eq.1 )then
                       u(i1,i2,i3,hx) = u(i1,i2,i3,hx) - expxi*pwc(3)
                       u(i1,i2,i3,hy) = u(i1,i2,i3,hy) - expxi*pwc(4)
                       u(i1,i2,i3,hz) = u(i1,i2,i3,hz) - expxi*pwc(5)
                     end if
                     xi = kx*(x-x0GP) + ky*(y-y0GP) + kz*(y-z0GP) - cc*
     & (t)
                     expxi = amp*exp(-betaGP*xi*xi )
                     un(i1,i2,i3,ex) = un(i1,i2,i3,ex) - expxi*pwc(0)
                     un(i1,i2,i3,ey) = un(i1,i2,i3,ey) - expxi*pwc(1)
                     un(i1,i2,i3,ez) = un(i1,i2,i3,ez) - expxi*pwc(2)
                     if( solveForAllFields.eq.1 )then
                       un(i1,i2,i3,hx) = un(i1,i2,i3,hx) - expxi*pwc(3)
                       un(i1,i2,i3,hy) = un(i1,i2,i3,hy) - expxi*pwc(4)
                       un(i1,i2,i3,hz) = un(i1,i2,i3,hz) - expxi*pwc(5)
                     end if
                   if( adjustThreeLevels.eq.1 )then
                       xi = kx*(x-x0GP) + ky*(y-y0GP) + kz*(y-z0GP) - 
     & cc*(t-2.*dt)
                       expxi = amp*exp(-betaGP*xi*xi )
                       um(i1,i2,i3,ex) = um(i1,i2,i3,ex) - expxi*pwc(0)
                       um(i1,i2,i3,ey) = um(i1,i2,i3,ey) - expxi*pwc(1)
                       um(i1,i2,i3,ez) = um(i1,i2,i3,ez) - expxi*pwc(2)
                       if( solveForAllFields.eq.1 )then
                         um(i1,i2,i3,hx) = um(i1,i2,i3,hx) - expxi*pwc(
     & 3)
                         um(i1,i2,i3,hy) = um(i1,i2,i3,hy) - expxi*pwc(
     & 4)
                         um(i1,i2,i3,hz) = um(i1,i2,i3,hz) - expxi*pwc(
     & 5)
                       end if
                   end if
                 else
                   write(*,'("adjust Inc. : ERROR: unknown incident 
     & field type")')
                   stop 88992
                 end if
               end do
               end do
               end do
             end if
             end if
            end do
            end do
        end if
      else if( adjustForIncidentField.eq.1 )then
        ! --- Add back the incident wave, before applying the non-reflecting BC ---
        if( adjustForBoundingBox )then
           if( debug.gt.4 )then
             write(*,'(" adjust +: kx,ky,kz,eps,cc",5e10.3)') kx,ky,kz,
     & eps,cc
           end if
           amp=1.
            ! adjust this many points near the boundary (needs to include width of extrapolation too!)
            ! halfWidth = orderOfAccuracy/2 ** now passed in**
            ! restrict all bounds to fit in the array dimensions: 
            dim(0,0)=nd1a
            dim(1,0)=nd1b
            dim(0,1)=nd2a
            dim(1,1)=nd2b
            dim(0,2)=nd3a
            dim(1,2)=nd3b
            ! gid(size,axis) : holds bounds on the points that still may need adjusting 
            !                : as we assign points near faces, gid will be reduced in size so that we 
            !                : do not adjust points multiple times.  
            gid(0,2)=0
            gid(1,2)=0
            range(0,2)=0
            range(1,2)=0
            do axis=0,nd-1
            do side=0,1
             if( boundaryCondition(side,axis).eq.abcPML )then
              gid(side,axis)=dim(side,axis)
             else
              gid(side,axis)=gridIndexRange(side,axis)-halfWidth*(1-2*
     & side)
             end if
            end do
            end do
            do axis=0,nd-1
            do side=0,1
             if( boundaryCondition(side,axis).ge.abcEM2 .and. 
     & boundaryCondition(side,axis).le.rbcLocal )then
              do a=0,nd-1
              do s=0,1
               range(s,a)=gid(s,a)
              end do
              end do
              if( boundaryCondition(side,axis).eq.abcPML )then
                ! we need to adjust more lines for a PML
                range(  side,axis)=dim(side,axis)
                range(1-side,axis)=gridIndexRange(side,axis)+(
     & numberLinesForPML+halfWidth)*(1-2*side) ! check this
              else
                range(0,axis)=gridIndexRange(side,axis)-halfWidth
                range(1,axis)=gridIndexRange(side,axis)+halfWidth
              end if
              ! restrict all bounds to fit in the array dimensions: 
              do a=0,nd-1
              do s=0,1
               range(s,a)=max(dim(0,a),min(dim(1,a), range(s,a) ))
              end do
              end do
              ! reduce gid bounds for future tangential directions so pts are not assigned twice
              ! gid(side,axis)=gridIndexRange(side,axis)+(halfWidth+1)*(1-2*side)
              gid(side,axis)=range(1-side,axis)+(1)*(1-2*side)
              n1a=range(0,0)
              n1b=range(1,0)
              n2a=range(0,1)
              n2b=range(1,1)
              n3a=range(0,2)
              n3b=range(1,2)
              if( debug.gt.0 )then
                write(*,'(" abc:adjustIncindent: side,axis,bc=",3i2," 
     & n1a,n1b,...=",3(i3,1x))') side,axis,boundaryCondition(side,
     & axis),n1a,n1b,n2a,n2b,n3a,n3b
                write(*,'(" abc:adjustIncindent: gid=",4i4)') gid(0,0),
     & gid(1,0),gid(0,1),gid(1,1)
                ! ' 
              end if
             if( nd.eq.2 )then
               do i3=n3a,n3b
               do i2=n2a,n2b
               do i1=n1a,n1b
                 if( gridType.eq.rectangular )then
                   x = xa(0)+i1*dx(0)
                   y = xa(1)+i2*dx(1)
                 else
                   x=xy(i1,i2,i3,0)
                   y=xy(i1,i2,i3,1)
                 end if
                 ! if( debug.gt.1 )then
                 !  t0=planeWave2Dhz0(x,y,t-dt)
                 !  write(*,'("nrbc: adjust +: i=",2i3," Hz,true=",2e10.3)') i1,i2,u(i1,i2,i3,hz),t0
                 ! end if
                   if( x.ge.icBoundingBox(0,0) .and. 
     & x.le.icBoundingBox(1,0) .and. y.ge.icBoundingBox(0,1) .and. 
     & y.le.icBoundingBox(1,1) )then
                     if( smoothBoundingBox.eq.1 )then
                       amp = (.5*(1.-tanh(beta*twoPi*(nv(0)*((x)-xv0(0)
     & )+nv(1)*((y)-xv0(1))-cc*(t)))))
                     else
                       amp=1.
                     end if
                   ! if( amp.gt. 1.e-9 )then
                   !   if( x.gt.1.5 .and. y.gt.-.05 .and. y.lt.0.5 )then
                   !     write(*,'(" x,amp=",2e10.2)') x,amp
                   !   end if
                 if( incidentFieldType .eq. planeWaveIncidentField )
     & then
                   ! --- plane wave incident field ---
                   u(i1,i2,i3,ex) = u(i1,i2,i3,ex) + amp*sin(twoPi*(kx*
     & (x)+ky*(y)-cc*(t-dt)))*pwc(0)
                   u(i1,i2,i3,ey) = u(i1,i2,i3,ey) + amp*sin(twoPi*(kx*
     & (x)+ky*(y)-cc*(t-dt)))*pwc(1)
                   u(i1,i2,i3,hz) = u(i1,i2,i3,hz) + amp*sin(twoPi*(kx*
     & (x)+ky*(y)-cc*(t-dt)))*pwc(5)
                   un(i1,i2,i3,ex)= un(i1,i2,i3,ex) + amp*sin(twoPi*(
     & kx*(x)+ky*(y)-cc*(t)))*pwc(0)
                   un(i1,i2,i3,ey)= un(i1,i2,i3,ey) + amp*sin(twoPi*(
     & kx*(x)+ky*(y)-cc*(t)))*pwc(1)
                   un(i1,i2,i3,hz)= un(i1,i2,i3,hz) + amp*sin(twoPi*(
     & kx*(x)+ky*(y)-cc*(t)))*pwc(5)
                   if( adjustThreeLevels.eq.1 )then
                    um(i1,i2,i3,ex) = um(i1,i2,i3,ex) + amp*sin(twoPi*(
     & kx*(x)+ky*(y)-cc*(t-2.*dt)))*pwc(0)
                    um(i1,i2,i3,ey) = um(i1,i2,i3,ey) + amp*sin(twoPi*(
     & kx*(x)+ky*(y)-cc*(t-2.*dt)))*pwc(1)
                    um(i1,i2,i3,hz) = um(i1,i2,i3,hz) + amp*sin(twoPi*(
     & kx*(x)+ky*(y)-cc*(t-2.*dt)))*pwc(5)
                   end if
                 else if( incidentFieldType .eq. 
     & gaussianPlaneWaveIncidentField )then
                   ! --- plane wave incident field ---
                     xi = kx*(x-x0GP) + ky*(y-y0GP) - cc*(t-dt)
                     expxi = amp*exp(-betaGP*xi*xi )
                     if( solveForAllFields.eq.0 )then
                       ! *new* way June 16, 2020 *wdh*
                       u(i1,i2,i3,ex) = u(i1,i2,i3,ex) + expxi*pwc(0)
                       u(i1,i2,i3,ey) = u(i1,i2,i3,ey) + expxi*pwc(1)
                       u(i1,i2,i3,hz) = u(i1,i2,i3,hz) + expxi*pwc(5)
                       ! *old way*
                       ! uhz = expxi
                       ! uex = uhz*(-ky/(eps*cc))
                       ! uey = uhz*( kx/(eps*cc))
                       ! u(i1,i2,i3,ex) = u(i1,i2,i3,ex) + uex
                       ! u(i1,i2,i3,ey) = u(i1,i2,i3,ey) + uey
                       ! u(i1,i2,i3,hz) = u(i1,i2,i3,hz) + uhz
                     else
                       u(i1,i2,i3,ex) = u(i1,i2,i3,ex) + expxi*pwc(0)
                       u(i1,i2,i3,ey) = u(i1,i2,i3,ey) + expxi*pwc(1)
                       u(i1,i2,i3,ez) = u(i1,i2,i3,ez) + expxi*pwc(2)
                       u(i1,i2,i3,hx) = u(i1,i2,i3,hx) + expxi*pwc(3)
                       u(i1,i2,i3,hy) = u(i1,i2,i3,hy) + expxi*pwc(4)
                       u(i1,i2,i3,hz) = u(i1,i2,i3,hz) + expxi*pwc(5)
                     end if
                     xi = kx*(x-x0GP) + ky*(y-y0GP) - cc*(t)
                     expxi = amp*exp(-betaGP*xi*xi )
                     if( solveForAllFields.eq.0 )then
                       ! *new* way June 16, 2020 *wdh*
                       un(i1,i2,i3,ex) = un(i1,i2,i3,ex) + expxi*pwc(0)
                       un(i1,i2,i3,ey) = un(i1,i2,i3,ey) + expxi*pwc(1)
                       un(i1,i2,i3,hz) = un(i1,i2,i3,hz) + expxi*pwc(5)
                       ! *old way*
                       ! uhz = expxi
                       ! uex = uhz*(-ky/(eps*cc))
                       ! uey = uhz*( kx/(eps*cc))
                       ! un(i1,i2,i3,ex) = un(i1,i2,i3,ex) + uex
                       ! un(i1,i2,i3,ey) = un(i1,i2,i3,ey) + uey
                       ! un(i1,i2,i3,hz) = un(i1,i2,i3,hz) + uhz
                     else
                       un(i1,i2,i3,ex) = un(i1,i2,i3,ex) + expxi*pwc(0)
                       un(i1,i2,i3,ey) = un(i1,i2,i3,ey) + expxi*pwc(1)
                       un(i1,i2,i3,ez) = un(i1,i2,i3,ez) + expxi*pwc(2)
                       un(i1,i2,i3,hx) = un(i1,i2,i3,hx) + expxi*pwc(3)
                       un(i1,i2,i3,hy) = un(i1,i2,i3,hy) + expxi*pwc(4)
                       un(i1,i2,i3,hz) = un(i1,i2,i3,hz) + expxi*pwc(5)
                     end if
                   if( adjustThreeLevels.eq.1 )then
                       xi = kx*(x-x0GP) + ky*(y-y0GP) - cc*(t-2.*dt)
                       expxi = amp*exp(-betaGP*xi*xi )
                       if( solveForAllFields.eq.0 )then
                         ! *new* way June 16, 2020 *wdh*
                         um(i1,i2,i3,ex) = um(i1,i2,i3,ex) + expxi*pwc(
     & 0)
                         um(i1,i2,i3,ey) = um(i1,i2,i3,ey) + expxi*pwc(
     & 1)
                         um(i1,i2,i3,hz) = um(i1,i2,i3,hz) + expxi*pwc(
     & 5)
                         ! *old way*
                         ! uhz = expxi
                         ! uex = uhz*(-ky/(eps*cc))
                         ! uey = uhz*( kx/(eps*cc))
                         ! um(i1,i2,i3,ex) = um(i1,i2,i3,ex) + uex
                         ! um(i1,i2,i3,ey) = um(i1,i2,i3,ey) + uey
                         ! um(i1,i2,i3,hz) = um(i1,i2,i3,hz) + uhz
                       else
                         um(i1,i2,i3,ex) = um(i1,i2,i3,ex) + expxi*pwc(
     & 0)
                         um(i1,i2,i3,ey) = um(i1,i2,i3,ey) + expxi*pwc(
     & 1)
                         um(i1,i2,i3,ez) = um(i1,i2,i3,ez) + expxi*pwc(
     & 2)
                         um(i1,i2,i3,hx) = um(i1,i2,i3,hx) + expxi*pwc(
     & 3)
                         um(i1,i2,i3,hy) = um(i1,i2,i3,hy) + expxi*pwc(
     & 4)
                         um(i1,i2,i3,hz) = um(i1,i2,i3,hz) + expxi*pwc(
     & 5)
                       end if
                   end if
                 else
                   write(*,'("adjust Inc. : ERROR: unknown incident 
     & field type")')
                   stop 88991
                 end if
                  endif
               end do
               end do
               end do
             else
               ! --- THREE DIMENSIONS ---
               !if( incidentFieldType .eq. gaussianPlaneWaveIncidentField )then
               !  write(*,'(" + YES : alter fields for gaussianPlaneWaveIncidentField ")')
               !end if 
               do i3=n3a,n3b
               do i2=n2a,n2b
               do i1=n1a,n1b
                 if( gridType.eq.rectangular )then
                   x = xa(0)+i1*dx(0)
                   y = xa(1)+i2*dx(1)
                   z = xa(2)+i3*dx(2)
                 else
                   x=xy(i1,i2,i3,0)
                   y=xy(i1,i2,i3,1)
                   z=xy(i1,i2,i3,2)
                 end if
                 if( x.ge.icBoundingBox(0,0) .and. x.le.icBoundingBox(
     & 1,0) .and. y.ge.icBoundingBox(0,1) .and. y.le.icBoundingBox(1,
     & 1) .and. z.ge.icBoundingBox(0,2) .and. z.le.icBoundingBox(1,2) 
     & )then
                   if( smoothBoundingBox.eq.1 )then
                     amp = (.5*(1.-tanh(beta*twoPi*(nv(0)*((x)-xv0(0))+
     & nv(1)*((y)-xv0(1))+nv(2)*((z)-xv0(2))-cc*(t)))))
                   else
                     amp=1.
                   end if
                 if( incidentFieldType .eq. planeWaveIncidentField )
     & then
                   u(i1,i2,i3,ex) = u(i1,i2,i3,ex) + amp*sin(twoPi*(kx*
     & (x)+ky*(y)+kz*(z)-cc*(t-dt)))*pwc(0)
                   u(i1,i2,i3,ey) = u(i1,i2,i3,ey) + amp*sin(twoPi*(kx*
     & (x)+ky*(y)+kz*(z)-cc*(t-dt)))*pwc(1)
                   u(i1,i2,i3,ez) = u(i1,i2,i3,ez) + amp*sin(twoPi*(kx*
     & (x)+ky*(y)+kz*(z)-cc*(t-dt)))*pwc(2)
                   un(i1,i2,i3,ex)=un(i1,i2,i3,ex) + amp*sin(twoPi*(kx*
     & (x)+ky*(y)+kz*(z)-cc*(t)))*pwc(0)
                   un(i1,i2,i3,ey)=un(i1,i2,i3,ey) + amp*sin(twoPi*(kx*
     & (x)+ky*(y)+kz*(z)-cc*(t)))*pwc(1)
                   un(i1,i2,i3,ez)=un(i1,i2,i3,ez) + amp*sin(twoPi*(kx*
     & (x)+ky*(y)+kz*(z)-cc*(t)))*pwc(2)
                   if( adjustThreeLevels.eq.1 )then
                    um(i1,i2,i3,ex)=um(i1,i2,i3,ex) + amp*sin(twoPi*(
     & kx*(x)+ky*(y)+kz*(z)-cc*(t-2.*dt)))*pwc(0)
                    um(i1,i2,i3,ey)=um(i1,i2,i3,ey) + amp*sin(twoPi*(
     & kx*(x)+ky*(y)+kz*(z)-cc*(t-2.*dt)))*pwc(1)
                    um(i1,i2,i3,ez)=um(i1,i2,i3,ez) + amp*sin(twoPi*(
     & kx*(x)+ky*(y)+kz*(z)-cc*(t-2.*dt)))*pwc(2)
                   end if
                 else if( incidentFieldType .eq. 
     & gaussianPlaneWaveIncidentField )then
                   ! --- Gaussian plane wave incident field ---  *wdh* Aug 18, 2019
                     xi = kx*(x-x0GP) + ky*(y-y0GP) + kz*(y-z0GP) - cc*
     & (t-dt)
                     expxi = amp*exp(-betaGP*xi*xi )
                     u(i1,i2,i3,ex) = u(i1,i2,i3,ex) + expxi*pwc(0)
                     u(i1,i2,i3,ey) = u(i1,i2,i3,ey) + expxi*pwc(1)
                     u(i1,i2,i3,ez) = u(i1,i2,i3,ez) + expxi*pwc(2)
                     if( solveForAllFields.eq.1 )then
                       u(i1,i2,i3,hx) = u(i1,i2,i3,hx) + expxi*pwc(3)
                       u(i1,i2,i3,hy) = u(i1,i2,i3,hy) + expxi*pwc(4)
                       u(i1,i2,i3,hz) = u(i1,i2,i3,hz) + expxi*pwc(5)
                     end if
                     xi = kx*(x-x0GP) + ky*(y-y0GP) + kz*(y-z0GP) - cc*
     & (t)
                     expxi = amp*exp(-betaGP*xi*xi )
                     un(i1,i2,i3,ex) = un(i1,i2,i3,ex) + expxi*pwc(0)
                     un(i1,i2,i3,ey) = un(i1,i2,i3,ey) + expxi*pwc(1)
                     un(i1,i2,i3,ez) = un(i1,i2,i3,ez) + expxi*pwc(2)
                     if( solveForAllFields.eq.1 )then
                       un(i1,i2,i3,hx) = un(i1,i2,i3,hx) + expxi*pwc(3)
                       un(i1,i2,i3,hy) = un(i1,i2,i3,hy) + expxi*pwc(4)
                       un(i1,i2,i3,hz) = un(i1,i2,i3,hz) + expxi*pwc(5)
                     end if
                   if( adjustThreeLevels.eq.1 )then
                       xi = kx*(x-x0GP) + ky*(y-y0GP) + kz*(y-z0GP) - 
     & cc*(t-2.*dt)
                       expxi = amp*exp(-betaGP*xi*xi )
                       um(i1,i2,i3,ex) = um(i1,i2,i3,ex) + expxi*pwc(0)
                       um(i1,i2,i3,ey) = um(i1,i2,i3,ey) + expxi*pwc(1)
                       um(i1,i2,i3,ez) = um(i1,i2,i3,ez) + expxi*pwc(2)
                       if( solveForAllFields.eq.1 )then
                         um(i1,i2,i3,hx) = um(i1,i2,i3,hx) + expxi*pwc(
     & 3)
                         um(i1,i2,i3,hy) = um(i1,i2,i3,hy) + expxi*pwc(
     & 4)
                         um(i1,i2,i3,hz) = um(i1,i2,i3,hz) + expxi*pwc(
     & 5)
                       end if
                   end if
                 else
                   write(*,'("adjust Inc. : ERROR: unknown incident 
     & field type")')
                   stop 88992
                 end if
                 endif
               end do
               end do
               end do
             end if
             end if
            end do
            end do
        else
           if( debug.gt.4 )then
             write(*,'(" adjust +: kx,ky,kz,eps,cc",5e10.3)') kx,ky,kz,
     & eps,cc
           end if
           amp=1.
            ! adjust this many points near the boundary (needs to include width of extrapolation too!)
            ! halfWidth = orderOfAccuracy/2 ** now passed in**
            ! restrict all bounds to fit in the array dimensions: 
            dim(0,0)=nd1a
            dim(1,0)=nd1b
            dim(0,1)=nd2a
            dim(1,1)=nd2b
            dim(0,2)=nd3a
            dim(1,2)=nd3b
            ! gid(size,axis) : holds bounds on the points that still may need adjusting 
            !                : as we assign points near faces, gid will be reduced in size so that we 
            !                : do not adjust points multiple times.  
            gid(0,2)=0
            gid(1,2)=0
            range(0,2)=0
            range(1,2)=0
            do axis=0,nd-1
            do side=0,1
             if( boundaryCondition(side,axis).eq.abcPML )then
              gid(side,axis)=dim(side,axis)
             else
              gid(side,axis)=gridIndexRange(side,axis)-halfWidth*(1-2*
     & side)
             end if
            end do
            end do
            do axis=0,nd-1
            do side=0,1
             if( boundaryCondition(side,axis).ge.abcEM2 .and. 
     & boundaryCondition(side,axis).le.rbcLocal )then
              do a=0,nd-1
              do s=0,1
               range(s,a)=gid(s,a)
              end do
              end do
              if( boundaryCondition(side,axis).eq.abcPML )then
                ! we need to adjust more lines for a PML
                range(  side,axis)=dim(side,axis)
                range(1-side,axis)=gridIndexRange(side,axis)+(
     & numberLinesForPML+halfWidth)*(1-2*side) ! check this
              else
                range(0,axis)=gridIndexRange(side,axis)-halfWidth
                range(1,axis)=gridIndexRange(side,axis)+halfWidth
              end if
              ! restrict all bounds to fit in the array dimensions: 
              do a=0,nd-1
              do s=0,1
               range(s,a)=max(dim(0,a),min(dim(1,a), range(s,a) ))
              end do
              end do
              ! reduce gid bounds for future tangential directions so pts are not assigned twice
              ! gid(side,axis)=gridIndexRange(side,axis)+(halfWidth+1)*(1-2*side)
              gid(side,axis)=range(1-side,axis)+(1)*(1-2*side)
              n1a=range(0,0)
              n1b=range(1,0)
              n2a=range(0,1)
              n2b=range(1,1)
              n3a=range(0,2)
              n3b=range(1,2)
              if( debug.gt.0 )then
                write(*,'(" abc:adjustIncindent: side,axis,bc=",3i2," 
     & n1a,n1b,...=",3(i3,1x))') side,axis,boundaryCondition(side,
     & axis),n1a,n1b,n2a,n2b,n3a,n3b
                write(*,'(" abc:adjustIncindent: gid=",4i4)') gid(0,0),
     & gid(1,0),gid(0,1),gid(1,1)
                ! ' 
              end if
             if( nd.eq.2 )then
               do i3=n3a,n3b
               do i2=n2a,n2b
               do i1=n1a,n1b
                 if( gridType.eq.rectangular )then
                   x = xa(0)+i1*dx(0)
                   y = xa(1)+i2*dx(1)
                 else
                   x=xy(i1,i2,i3,0)
                   y=xy(i1,i2,i3,1)
                 end if
                 ! if( debug.gt.1 )then
                 !  t0=planeWave2Dhz0(x,y,t-dt)
                 !  write(*,'("nrbc: adjust +: i=",2i3," Hz,true=",2e10.3)') i1,i2,u(i1,i2,i3,hz),t0
                 ! end if
                 if( incidentFieldType .eq. planeWaveIncidentField )
     & then
                   ! --- plane wave incident field ---
                   u(i1,i2,i3,ex) = u(i1,i2,i3,ex) + amp*sin(twoPi*(kx*
     & (x)+ky*(y)-cc*(t-dt)))*pwc(0)
                   u(i1,i2,i3,ey) = u(i1,i2,i3,ey) + amp*sin(twoPi*(kx*
     & (x)+ky*(y)-cc*(t-dt)))*pwc(1)
                   u(i1,i2,i3,hz) = u(i1,i2,i3,hz) + amp*sin(twoPi*(kx*
     & (x)+ky*(y)-cc*(t-dt)))*pwc(5)
                   un(i1,i2,i3,ex)= un(i1,i2,i3,ex) + amp*sin(twoPi*(
     & kx*(x)+ky*(y)-cc*(t)))*pwc(0)
                   un(i1,i2,i3,ey)= un(i1,i2,i3,ey) + amp*sin(twoPi*(
     & kx*(x)+ky*(y)-cc*(t)))*pwc(1)
                   un(i1,i2,i3,hz)= un(i1,i2,i3,hz) + amp*sin(twoPi*(
     & kx*(x)+ky*(y)-cc*(t)))*pwc(5)
                   if( adjustThreeLevels.eq.1 )then
                    um(i1,i2,i3,ex) = um(i1,i2,i3,ex) + amp*sin(twoPi*(
     & kx*(x)+ky*(y)-cc*(t-2.*dt)))*pwc(0)
                    um(i1,i2,i3,ey) = um(i1,i2,i3,ey) + amp*sin(twoPi*(
     & kx*(x)+ky*(y)-cc*(t-2.*dt)))*pwc(1)
                    um(i1,i2,i3,hz) = um(i1,i2,i3,hz) + amp*sin(twoPi*(
     & kx*(x)+ky*(y)-cc*(t-2.*dt)))*pwc(5)
                   end if
                 else if( incidentFieldType .eq. 
     & gaussianPlaneWaveIncidentField )then
                   ! --- plane wave incident field ---
                     xi = kx*(x-x0GP) + ky*(y-y0GP) - cc*(t-dt)
                     expxi = amp*exp(-betaGP*xi*xi )
                     if( solveForAllFields.eq.0 )then
                       ! *new* way June 16, 2020 *wdh*
                       u(i1,i2,i3,ex) = u(i1,i2,i3,ex) + expxi*pwc(0)
                       u(i1,i2,i3,ey) = u(i1,i2,i3,ey) + expxi*pwc(1)
                       u(i1,i2,i3,hz) = u(i1,i2,i3,hz) + expxi*pwc(5)
                       ! *old way*
                       ! uhz = expxi
                       ! uex = uhz*(-ky/(eps*cc))
                       ! uey = uhz*( kx/(eps*cc))
                       ! u(i1,i2,i3,ex) = u(i1,i2,i3,ex) + uex
                       ! u(i1,i2,i3,ey) = u(i1,i2,i3,ey) + uey
                       ! u(i1,i2,i3,hz) = u(i1,i2,i3,hz) + uhz
                     else
                       u(i1,i2,i3,ex) = u(i1,i2,i3,ex) + expxi*pwc(0)
                       u(i1,i2,i3,ey) = u(i1,i2,i3,ey) + expxi*pwc(1)
                       u(i1,i2,i3,ez) = u(i1,i2,i3,ez) + expxi*pwc(2)
                       u(i1,i2,i3,hx) = u(i1,i2,i3,hx) + expxi*pwc(3)
                       u(i1,i2,i3,hy) = u(i1,i2,i3,hy) + expxi*pwc(4)
                       u(i1,i2,i3,hz) = u(i1,i2,i3,hz) + expxi*pwc(5)
                     end if
                     xi = kx*(x-x0GP) + ky*(y-y0GP) - cc*(t)
                     expxi = amp*exp(-betaGP*xi*xi )
                     if( solveForAllFields.eq.0 )then
                       ! *new* way June 16, 2020 *wdh*
                       un(i1,i2,i3,ex) = un(i1,i2,i3,ex) + expxi*pwc(0)
                       un(i1,i2,i3,ey) = un(i1,i2,i3,ey) + expxi*pwc(1)
                       un(i1,i2,i3,hz) = un(i1,i2,i3,hz) + expxi*pwc(5)
                       ! *old way*
                       ! uhz = expxi
                       ! uex = uhz*(-ky/(eps*cc))
                       ! uey = uhz*( kx/(eps*cc))
                       ! un(i1,i2,i3,ex) = un(i1,i2,i3,ex) + uex
                       ! un(i1,i2,i3,ey) = un(i1,i2,i3,ey) + uey
                       ! un(i1,i2,i3,hz) = un(i1,i2,i3,hz) + uhz
                     else
                       un(i1,i2,i3,ex) = un(i1,i2,i3,ex) + expxi*pwc(0)
                       un(i1,i2,i3,ey) = un(i1,i2,i3,ey) + expxi*pwc(1)
                       un(i1,i2,i3,ez) = un(i1,i2,i3,ez) + expxi*pwc(2)
                       un(i1,i2,i3,hx) = un(i1,i2,i3,hx) + expxi*pwc(3)
                       un(i1,i2,i3,hy) = un(i1,i2,i3,hy) + expxi*pwc(4)
                       un(i1,i2,i3,hz) = un(i1,i2,i3,hz) + expxi*pwc(5)
                     end if
                   if( adjustThreeLevels.eq.1 )then
                       xi = kx*(x-x0GP) + ky*(y-y0GP) - cc*(t-2.*dt)
                       expxi = amp*exp(-betaGP*xi*xi )
                       if( solveForAllFields.eq.0 )then
                         ! *new* way June 16, 2020 *wdh*
                         um(i1,i2,i3,ex) = um(i1,i2,i3,ex) + expxi*pwc(
     & 0)
                         um(i1,i2,i3,ey) = um(i1,i2,i3,ey) + expxi*pwc(
     & 1)
                         um(i1,i2,i3,hz) = um(i1,i2,i3,hz) + expxi*pwc(
     & 5)
                         ! *old way*
                         ! uhz = expxi
                         ! uex = uhz*(-ky/(eps*cc))
                         ! uey = uhz*( kx/(eps*cc))
                         ! um(i1,i2,i3,ex) = um(i1,i2,i3,ex) + uex
                         ! um(i1,i2,i3,ey) = um(i1,i2,i3,ey) + uey
                         ! um(i1,i2,i3,hz) = um(i1,i2,i3,hz) + uhz
                       else
                         um(i1,i2,i3,ex) = um(i1,i2,i3,ex) + expxi*pwc(
     & 0)
                         um(i1,i2,i3,ey) = um(i1,i2,i3,ey) + expxi*pwc(
     & 1)
                         um(i1,i2,i3,ez) = um(i1,i2,i3,ez) + expxi*pwc(
     & 2)
                         um(i1,i2,i3,hx) = um(i1,i2,i3,hx) + expxi*pwc(
     & 3)
                         um(i1,i2,i3,hy) = um(i1,i2,i3,hy) + expxi*pwc(
     & 4)
                         um(i1,i2,i3,hz) = um(i1,i2,i3,hz) + expxi*pwc(
     & 5)
                       end if
                   end if
                 else
                   write(*,'("adjust Inc. : ERROR: unknown incident 
     & field type")')
                   stop 88991
                 end if
               end do
               end do
               end do
             else
               ! --- THREE DIMENSIONS ---
               !if( incidentFieldType .eq. gaussianPlaneWaveIncidentField )then
               !  write(*,'(" + NO : alter fields for gaussianPlaneWaveIncidentField ")')
               !end if 
               do i3=n3a,n3b
               do i2=n2a,n2b
               do i1=n1a,n1b
                 if( gridType.eq.rectangular )then
                   x = xa(0)+i1*dx(0)
                   y = xa(1)+i2*dx(1)
                   z = xa(2)+i3*dx(2)
                 else
                   x=xy(i1,i2,i3,0)
                   y=xy(i1,i2,i3,1)
                   z=xy(i1,i2,i3,2)
                 end if
                 if( incidentFieldType .eq. planeWaveIncidentField )
     & then
                   u(i1,i2,i3,ex) = u(i1,i2,i3,ex) + amp*sin(twoPi*(kx*
     & (x)+ky*(y)+kz*(z)-cc*(t-dt)))*pwc(0)
                   u(i1,i2,i3,ey) = u(i1,i2,i3,ey) + amp*sin(twoPi*(kx*
     & (x)+ky*(y)+kz*(z)-cc*(t-dt)))*pwc(1)
                   u(i1,i2,i3,ez) = u(i1,i2,i3,ez) + amp*sin(twoPi*(kx*
     & (x)+ky*(y)+kz*(z)-cc*(t-dt)))*pwc(2)
                   un(i1,i2,i3,ex)=un(i1,i2,i3,ex) + amp*sin(twoPi*(kx*
     & (x)+ky*(y)+kz*(z)-cc*(t)))*pwc(0)
                   un(i1,i2,i3,ey)=un(i1,i2,i3,ey) + amp*sin(twoPi*(kx*
     & (x)+ky*(y)+kz*(z)-cc*(t)))*pwc(1)
                   un(i1,i2,i3,ez)=un(i1,i2,i3,ez) + amp*sin(twoPi*(kx*
     & (x)+ky*(y)+kz*(z)-cc*(t)))*pwc(2)
                   if( adjustThreeLevels.eq.1 )then
                    um(i1,i2,i3,ex)=um(i1,i2,i3,ex) + amp*sin(twoPi*(
     & kx*(x)+ky*(y)+kz*(z)-cc*(t-2.*dt)))*pwc(0)
                    um(i1,i2,i3,ey)=um(i1,i2,i3,ey) + amp*sin(twoPi*(
     & kx*(x)+ky*(y)+kz*(z)-cc*(t-2.*dt)))*pwc(1)
                    um(i1,i2,i3,ez)=um(i1,i2,i3,ez) + amp*sin(twoPi*(
     & kx*(x)+ky*(y)+kz*(z)-cc*(t-2.*dt)))*pwc(2)
                   end if
                 else if( incidentFieldType .eq. 
     & gaussianPlaneWaveIncidentField )then
                   ! --- Gaussian plane wave incident field ---  *wdh* Aug 18, 2019
                     xi = kx*(x-x0GP) + ky*(y-y0GP) + kz*(y-z0GP) - cc*
     & (t-dt)
                     expxi = amp*exp(-betaGP*xi*xi )
                     u(i1,i2,i3,ex) = u(i1,i2,i3,ex) + expxi*pwc(0)
                     u(i1,i2,i3,ey) = u(i1,i2,i3,ey) + expxi*pwc(1)
                     u(i1,i2,i3,ez) = u(i1,i2,i3,ez) + expxi*pwc(2)
                     if( solveForAllFields.eq.1 )then
                       u(i1,i2,i3,hx) = u(i1,i2,i3,hx) + expxi*pwc(3)
                       u(i1,i2,i3,hy) = u(i1,i2,i3,hy) + expxi*pwc(4)
                       u(i1,i2,i3,hz) = u(i1,i2,i3,hz) + expxi*pwc(5)
                     end if
                     xi = kx*(x-x0GP) + ky*(y-y0GP) + kz*(y-z0GP) - cc*
     & (t)
                     expxi = amp*exp(-betaGP*xi*xi )
                     un(i1,i2,i3,ex) = un(i1,i2,i3,ex) + expxi*pwc(0)
                     un(i1,i2,i3,ey) = un(i1,i2,i3,ey) + expxi*pwc(1)
                     un(i1,i2,i3,ez) = un(i1,i2,i3,ez) + expxi*pwc(2)
                     if( solveForAllFields.eq.1 )then
                       un(i1,i2,i3,hx) = un(i1,i2,i3,hx) + expxi*pwc(3)
                       un(i1,i2,i3,hy) = un(i1,i2,i3,hy) + expxi*pwc(4)
                       un(i1,i2,i3,hz) = un(i1,i2,i3,hz) + expxi*pwc(5)
                     end if
                   if( adjustThreeLevels.eq.1 )then
                       xi = kx*(x-x0GP) + ky*(y-y0GP) + kz*(y-z0GP) - 
     & cc*(t-2.*dt)
                       expxi = amp*exp(-betaGP*xi*xi )
                       um(i1,i2,i3,ex) = um(i1,i2,i3,ex) + expxi*pwc(0)
                       um(i1,i2,i3,ey) = um(i1,i2,i3,ey) + expxi*pwc(1)
                       um(i1,i2,i3,ez) = um(i1,i2,i3,ez) + expxi*pwc(2)
                       if( solveForAllFields.eq.1 )then
                         um(i1,i2,i3,hx) = um(i1,i2,i3,hx) + expxi*pwc(
     & 3)
                         um(i1,i2,i3,hy) = um(i1,i2,i3,hy) + expxi*pwc(
     & 4)
                         um(i1,i2,i3,hz) = um(i1,i2,i3,hz) + expxi*pwc(
     & 5)
                       end if
                   end if
                 else
                   write(*,'("adjust Inc. : ERROR: unknown incident 
     & field type")')
                   stop 88992
                 end if
               end do
               end do
               end do
             end if
             end if
            end do
            end do
        end if
      else
        write(*,'(" adjustForIncident: unexpected value for 
     & adjustForIncidentField=",i6)') adjustForIncidentField
        ! ' 
        stop 20031
      end if

      return
      end
