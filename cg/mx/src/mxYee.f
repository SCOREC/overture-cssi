! This file automatically generated from mxYee.bf with bpp.
c *******************************************************************************
c     Define the Yee Approximation for Cartesian Grids
c *******************************************************************************

c These next include files will define the macros that will define the difference approximations
c The actual macro is called below
c#Include "defineDiffOrder2f.h"
c#Include "defineDiffOrder4f.h"


c Here are macros that define the planeWave solution
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





c ------------ macros for the plane material interface -------------------------
c Incident + reflected: 
c Transmitted: 

! =============================================================================
! Evaluate the plane material interface and return the results in ue(..)
!  Input: x0,y0,z0, xp,yp,zp
!  Input: te,th : eval E and H at these times
!  Input: DIM : 2 or 3 
!  Assumes that the initialize routine has been called
! =============================================================================




c ************************************************************************************************
c  This macro is used for looping over the faces of a grid to assign booundary conditions
c
c extra: extra points to assign
c          Case 1: extra=numberOfGhostPoints -- for assigning extended boundaries
c          Case 2: extra=-1 -- for assigning ghost points but not including extended boundaries
c numberOfGhostPoints : number of ghost points (1 for 2nd order, 2 for fourth-order ...)
c ***********************************************************************************************


c ========================================================================
c Begin loop over edges in 3D
c ========================================================================



! ===============================================================================
!  Assign Dirichlet boundary conditions to the E field
!   OPTION : TZ, planeWave, planeScat, planeMat
! ===============================================================================


! ===============================================================================
!  Macro: 
!      Assign PEC boundary conditions to the stair-step boundary of bodies.
! 
!   OPTION : TZ or planeWave or planeScat 
! ===============================================================================







! ==================================================================================
! Compute the material properties for E or H computations
!
! FIELD : E or H
! DIM   : 2 or 3x
! ==================================================================================

! ====================================================================================
!  Return the material speed of light at the cell center
! ====================================================================================





! ========================================================================
! Apply the second order accurate Mur version of the EM2 absorbing BC
!  AXIS : R or S or T 
! ========================================================================

! ==========================================================================================================
!  Set the ghost point values for tangential Magnetic components on a face (side,axis)
!   e2,e3 : tangential components to set 
!   side1,dir1 : a face adjacent to (side,axis)
!   j3,m3a,m3b : a loop over the third direction orthogonal to axis and dir1
!   
! ==========================================================================================================

! =============================================================================================
! Apply the second order accurate Mur version of the EM2 absorbing BC for the Magnetic field
!  AXIS : R or S or T 
! ==========================================================================================


! vertex coordinates: 

! 

      subroutine mxYee( nd, nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,ndf1a,ndf1b,
     & ndf2a,ndf2b,ndf3a,ndf3b,gridIndexRange, um, u, un, f, media, 
     & epsv, muv, sigmaEv, sigmaHv, mask, uKnown, boundaryCondition, 
     & ipar, rpar, ierr )
c ====================================================================================================
c
c   ---- Yee Approximation for Cartesian Grids ---
c
c  gridType : 0=rectangular, 1=curvilinear
c  useForcing : 1=use f for RHS to BC
c
c  um (input)  : old solution     E(t-dt), H(t-dt/2)  (used by ABC)
c  u  (input)  : current solution E(t   ), H(t+dt/2)
c  un (output) : new solution     E(t+dt), H(t+3*dt/2)
c
c  uKnown : holds known solution for some cases
c 
c  option = ipar[0] : 
c       option = 1 : advance E
c       option = 2 : advance H 
c       option = 1+2 : advance E and H 
c
c maskBodies : if maskBodies=1 then the mask(i1,i2,i3) defines cells that are inside a body
c =====================================================================================================

      implicit none

      integer nd, nd1a,nd1b,nd2a,nd2b,nd3a,nd3b, ndf1a,ndf1b,ndf2a,
     & ndf2b,ndf3a,ndf3b,n1a,n1b,n2a,n2b,n3a,n3b, ndc, bc,ierr

      real um(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
      real u(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
      real un(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
      real f(ndf1a:ndf1b,ndf2a:ndf2b,ndf3a:ndf3b,0:*)
      real uKnown(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
      real epsv(0:*),muv(0:*),sigmaEv(0:*),sigmaHv(0:*)
      integer media(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
      integer mask(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
      integer gridIndexRange(0:1,0:2)

      integer ipar(0:*),boundaryCondition(0:1,0:2)
      real rpar(0:*)

c     --- local variables ----

      integer side,axis,useForcing,ex,ey,ez,hx,hy,hz,grid,debug,side1,
     & side2,side3,forcingOption,option
      real dx(0:2),t,ep,dt,c
      integer addSourceTerm,initialConditionOption,maskBodies,
     & knownSolutionOption,useTwilightZoneMaterials
      integer axisp1,axisp2,i1,i2,i3,is1,is2,is3,js1,js2,js3,ks1,ks2,
     & ks3,is
      integer ip1,ip2,ip3,ig1,ig2,ig3,ghost1,ghost2,ghost3
      integer extra,extra1a,extra1b,extra2a,extra2b,extra3a,extra3b,
     & numberOfGhostPoints
      integer edgeDirection,sidea,sideb,sidec,bc1,bc2,bc3

      real eps,mu,kx,ky,kz,slowStartInterval,twoPi,cc

      integer nra(0:2)
      real xv(0:2),xa(0:2),pwc(0:5)
      real ue(0:9),uet(0:9),uex(0:9),uey(0:9),uez(0:9)
      ! real sigmaE,sigmaH,dc,ca,cbx,cby,cbz,cj
      real cf
      real x0,y0,z0,xp,yp,zp,te,tep,tef,th,thp,thf
      integer et1,et2,hn,ht1,ht2,bc0
      real coste,sinte,costh,sinth
      integer epsc,muc,sigmaEc,sigmaHc

      real xh,yh,zh,ue00(0:5),ue10(0:5),ue01(0:5),ue11(0:5)

      integer m00,m10,m01,m11
      integer m000,m100,m010,m110,m001,m101,m011,m111

      real eps00,eps10,eps01,eps11
      real eps000,eps100,eps010,eps110,eps001,eps101,eps011,eps111

      real sigmaE00,sigmaE10,sigmaE01,sigmaE11
      real sigmaE000,sigmaE100,sigmaE010,sigmaE110,sigmaE001,sigmaE101,
     & sigmaE011,sigmaE111

      real mu00,mu10,mu01,mu11
      real mu000,mu100,mu010,mu110,mu001,mu101,mu011,mu111

      real sigmaH00,sigmaH10,sigmaH01,sigmaH11
      real sigmaH000,sigmaH100,sigmaH010,sigmaH110,sigmaH001,sigmaH101,
     & sigmaH011,sigmaH111

      real epsEx,epsEy,epsEz,sigmaEx,sigmaEy,sigmaEz
      real muHx,muHy,muHz,sigmaHx,sigmaHy,sigmaHz

      real dcEx,caEx,cbxEx,cbyEx,cbzEx,cjEx
      real dcEy,caEy,cbxEy,cbyEy,cbzEy,cjEy
      real dcEz,caEz,cbxEz,cbyEz,cbzEz,cjEz

      real dcHx,caHx,cbxHx,cbyHx,cbzHx,cjHx
      real dcHy,caHy,cbxHy,cbyHy,cbzHy,cjHy
      real dcHz,caHz,cbxHz,cbyHz,cbzHz,cjHz

      integer i1a,i1b,i2a,i2b,i3a,i3b,dir1,dir2,dir3,i
      integer i1p,i2p,i3p
      real ctx1,ctx2,ctx3,ctx4

      ! for plane material interfaces:
      integer numberOfPMC
      real pmc(0:40)
      real xPMI(0:2),nPMI(0:2)

      integer myid
      integer debugFile
      character*20 debugFileName

      ! boundary conditions parameters
! define BC parameters for fortran routines
! boundary conditions
      integer dirichlet,perfectElectricalConductor,
     & perfectMagneticConductor,planeWaveBoundaryCondition,
     & interfaceBC,symmetryBoundaryCondition,abcEM2,abcPML,abc3,abc4,
     & abc5,rbcNonLocal,rbcLocal,lastBC
      parameter( dirichlet=1,perfectElectricalConductor=2,
     & perfectMagneticConductor=3,planeWaveBoundaryCondition=4,
     & symmetryBoundaryCondition=5,interfaceBC=6,abcEM2=7,abcPML=8,
     & abc3=9,abc4=10,abc5=11,rbcNonLocal=12,rbcLocal=13,lastBC=13 )

      integer rectangular,curvilinear
      parameter(rectangular=0,curvilinear=1)

      ! initial condition parameters
!
! define initial condition options (see Maxwell.h :  enum InitialConditionEnum)
!
      integer defaultInitialCondition,planeWaveInitialCondition,
     & gaussianPlaneWave,gaussianPulseInitialCondition,
     & squareEigenfunctionInitialCondition, 
     & annulusEigenfunctionInitialCondition,zeroInitialCondition,
     & planeWaveScatteredFieldInitialCondition,
     & planeMaterialInterfaceInitialCondition,
     & gaussianIntegralInitialCondition,twilightZoneInitialConditions,
     & userDefinedInitialConditionOption

      parameter( defaultInitialCondition              =0,      
     & planeWaveInitialCondition                =1,      
     & gaussianPlaneWave                        =2,      
     & gaussianPulseInitialCondition            =3,      
     & squareEigenfunctionInitialCondition      =4,      
     & annulusEigenfunctionInitialCondition     =5,      
     & zeroInitialCondition                     =6,      
     & planeWaveScatteredFieldInitialCondition  =7,      
     & planeMaterialInterfaceInitialCondition   =8,      
     & gaussianIntegralInitialCondition         =9,   
     & twilightZoneInitialConditions        =10,
     & userDefinedInitialConditionOption        =11)

      ! forcing options
      ! forcingOptions -- these should match ForcingEnum in Maxwell.h 
      integer noForcing,magneticSinusoidalPointSource,gaussianSource,
     & twilightZoneForcing, gaussianChargeSource, 
     & userDefinedForcingOption
      integer noBoundaryForcing,planeWaveBoundaryForcing,
     & chirpedPlaneWaveBoundaryForcing
      parameter(noForcing                =0,
     & magneticSinusoidalPointSource =1,gaussianSource                
     & =2,twilightZoneForcing           =3,    gaussianChargeSource   
     &        =4,userDefinedForcingOption      =5 )
      ! boundary forcing options when solved directly for the scattered field:
      parameter( noBoundaryForcing              =0,   
     & planeWaveBoundaryForcing       =1,
     & chirpedPlaneWaveBoundaryForcing=2 )

      ! known solutions

      integer noKnownSolution,twilightZoneKnownSolution,
     & planeWaveKnownSolution,gaussianPlaneWaveKnownSolution,
     & gaussianIntegralKnownSolution,
     & planeMaterialInterfaceKnownSolution,
     & scatteringFromADiskKnownSolution,
     & scatteringFromADielectricDiskKnownSolution,
     & scatteringFromASphereKnownSolution,
     & scatteringFromADielectricSphereKnownSolution,
     & squareEigenfunctionKnownSolution,
     & annulusEigenfunctionKnownSolution,
     & eigenfunctionsOfASphereKnownSolution,userDefinedKnownSolution

      parameter( noKnownSolution=0,twilightZoneKnownSolution=1,
     & planeWaveKnownSolution=2,gaussianPlaneWaveKnownSolution=3,
     & gaussianIntegralKnownSolution=4,
     & planeMaterialInterfaceKnownSolution=5,
     & scatteringFromADiskKnownSolution=6,
     & scatteringFromADielectricDiskKnownSolution=7,
     & scatteringFromASphereKnownSolution=8,
     & scatteringFromADielectricSphereKnownSolution=9,
     & squareEigenfunctionKnownSolution=10,
     & annulusEigenfunctionKnownSolution=11,
     & eigenfunctionsOfASphereKnownSolution=12,
     & userDefinedKnownSolution=13 )

c     --- start statement function ----
      integer kd,m,n

c............... end statement functions
      save debugFile
      data debugFile/-1/

      ierr=0

      option                =ipar( 0)
      ex                    =ipar( 1)
      ey                    =ipar( 2)
      ez                    =ipar( 3)
      hx                    =ipar( 4)
      hy                    =ipar( 5)
      hz                    =ipar( 6)
      epsc                  =ipar( 7)  ! component location for TZ values of eps
      muc                   =ipar( 8)
      sigmaEc               =ipar( 9)
      sigmaHc               =ipar(10)
      grid                  =ipar(11)
      debug                 =ipar(12)
      addSourceTerm         =ipar(13)
      useForcing            =ipar(14)
      forcingOption         =ipar(15)
      nra(0)                =ipar(16)  ! for computing x,y,z coordinates
      nra(1)                =ipar(17)
      nra(2)                =ipar(18)
      initialConditionOption=ipar(19)
      maskBodies            =ipar(20)
      knownSolutionOption   =ipar(21)
      useTwilightZoneMaterials=ipar(22)
      myid                  =ipar(23)
      numberOfPMC           =ipar(24)

      dx(0)                =rpar(0)
      dx(1)                =rpar(1)
      dx(2)                =rpar(2)
      t                    =rpar(3)
      ep                   =rpar(4) ! pointer for exact solution
      dt                   =rpar(5)
      c                    =rpar(6)
      kx                   =rpar(7)  ! for plane wave forcing
      ky                   =rpar(8)
      kz                   =rpar(9)
      slowStartInterval    =rpar(10)
      xa(0)                =rpar(11) ! for computing x,y,z coordinates
      xa(1)                =rpar(12)
      xa(2)                =rpar(13)
      pwc(0)               =rpar(14) ! for the plane wave solution
      pwc(1)               =rpar(15)
      pwc(2)               =rpar(16)
      pwc(3)               =rpar(17)
      pwc(4)               =rpar(18)
      pwc(5)               =rpar(19)

      ! plane material interface:
      do i=0,numberOfPMC-1
        pmc(i)            =rpar(i+20)
      end do
      do i=0,2
        xPMI(i)=pmc(28+i)
        nPMI(i)=pmc(30+i)
      end do

      if( debug.gt.0 .and. debugFile.lt.0 )then
        ! open debug files
        debugFile=10
        if( myid.lt.10 )then
          write(debugFileName,'("mxYee",i1,".fdebug")') myid
        else
          write(debugFileName,'("mxYee",i4.4,".fdebug")') myid
        end if
        ! write(*,*) 'mxYee: myid=',myid,' open debug file:',debugFileName
        open(debugFile,file=debugFileName,status='unknown',
     & form='formatted')
        ! '
        ! INQUIRE(FILE=filen, EXIST=filex)
      end if


      if( debug.gt.1 )then
       write(debugFile,'(/," mxYee: **START** useForcing=",i4," 
     & forcingOption=",i4)') useForcing,forcingOption

       write(debugFile,'(" mxYee: option,ex,ey,ez,hx,hy,hz=",7i4)') 
     & option,ex,ey,ez,hx,hy,hz
       write(debugFile,'(" mxYee: t,dt,dx(0),dx(1),dx(2)=",3e10.2)') t,
     & dt,dx(0),dx(1),dx(2)
       write(debugFile,'(" mxYee: xa(0),xa(1),xa(2)=",3e10.2)') xa(0),
     & xa(1),xa(2)
       write(debugFile,'(" mxYee: nra(0),nra(1),nra(2)=",3i6)') nra(0),
     & nra(1),nra(2)
       write(debugFile,'(" mxYee: epsv(0),muv(0),sigmaEv(0),sigmaHv(0)
     & =",4f8.4)') epsv(0),muv(0),sigmaEv(0),sigmaHv(0)
       write(debugFile,'(" mxYee: gridIndexRange=",3(2i6,2x))') ((
     & gridIndexRange(side,axis),side=0,1),axis=0,2)
       write(debugFile,'(" mxYee: boundaryCondition=",3(2i6,2x))') ((
     & boundaryCondition(side,axis),side=0,1),axis=0,2)
       ! '
      end if

      if( ex.lt.0 .or. ey.lt.0 .or. hz.lt.0 .or. epsc.lt.0 .or. 
     & muc.lt.0 .or. sigmaEc.lt.0 .or. sigmaHc.lt.0 .or.(nd.eq.3 
     & .and. ( ez.lt.0 .or. hx.lt.0 .or. hy.lt.0) ) )then
        write(*,'(" mxYee:ERROR: invalid component number (must be >=0)
     & :")')
        write(*,'(" ex,ey,hz =",3i5)') ex,ey,hz
        write(*,'(" epsc,muc,sigmaEc,sigmaHc =",4i4)') epsc,muc,
     & sigmaEc,sigmaHc
        stop 1110
      end if


      ! for plane wave forcing 
      twoPi=8.*atan2(1.,1.)
      cc= c*sqrt( kx*kx+ky*ky+kz*kz )

!      epsX=1.e-30 ! fix this ***

      te = t         ! E lives at this time (u)
      th = t+.5*dt   ! H lives at this time

      tep = te+dt    ! new E will be at this time (un)
      thp = th+dt    ! new H will be at this time

      tef = te+.5*dt ! forcing for E is centered at this time
      thf = th+.5*dt ! forcing for H is centered at this time

      n1a=gridIndexRange(0,0)
      n1b=gridIndexRange(1,0)
      n2a=gridIndexRange(0,1)
      n2b=gridIndexRange(1,1)
      n3a=gridIndexRange(0,2)
      n3b=gridIndexRange(1,2)


      ! if( initialConditionOption .eq. planeWaveScatteredFieldInitialCondition )then
      coste = cos(-twoPi*cc*tep)
      sinte = sin(-twoPi*cc*tep)
      costh = cos(-twoPi*cc*thp)
      sinth = sin(-twoPi*cc*thp)

      if( useForcing.ne.0 )then
        cf=1.
      else
        cf=0.  ! turn off the forcing
      end if


      if( debug.ge.16 .and. forcingOption.eq.twilightZoneForcing )then
      if( nd.eq.2 )then
       ! --- print errors ---
       write(debugFile,'(/," mxYee: ERRORS at start E, t=",e10.2)') te
       z0=0.
       do i3=n3a,n3b
       do i2=n2a,n2b
       do i1=n1a,n1b
        x0 = (xa(0)+dx(0)*(i1-nra(0)))
        y0 = (xa(1)+dx(1)*(i2-nra(1)))
        xp = x0 + .5*dx(0)
        yp = y0 + .5*dx(1)

        call ogDeriv(ep,0,0,0,0,xp,y0,z0,te,ex,ue (ex))
        call ogDeriv(ep,0,0,0,0,x0,yp,z0,te,ey,ue (ey))
        write(debugFile,'(" mxYee: t=",e10.2," i1,i2=",2i4," Ex,Ey =",
     & 2e10.2," err=",2e10.2)') te,i1,i2,u(i1,i2,i3,ex),u(i1,i2,i3,ey)
     & ,u(i1,i2,i3,ex)-ue(ex),u(i1,i2,i3,ey)-ue(ey)
          ! ' 
       end do
       end do
       end do
      end if
      end if


      ! --- Advance E ---
      if( nd.eq.2 )then
        ! --- Advance E : 2D ---
       do i3=n3a,n3b
       do i2=n2a,n2b
       do i1=n1a,n1b

        ! m = media(i1,i2,i3)  ! pointer to media properties
        ! eps = epsv(m)
        ! mu  = muv(m)
        ! sigmaE = sigmaEv(m)

        ! this next macro computes the properties epsEx, epsEy, sigmaEx, sigmaEy 
         ! Here we assume the material properties are cell centered
          ! ----- 2D ----
          m00 = media(i1  ,i2  ,i3)  ! pointer to media properties in cell (i1+.5,i2+.5)
           m10 = media(i1-1,i2  ,i3)  ! pointer to media properties in cell (i1-.5,i2+.5)
           m01 = media(i1  ,i2-1,i3)  ! pointer to media properties in cell (i1+.5,i2-.5)
           if( m00.eq.m10 .and. m00.eq.m01 )then
             ! locally uniform material properties: 
             epsEx = epsv(m00)
             epsEy = epsEx
             sigmaEx = sigmaEv(m00)
             sigmaEy = sigmaEx
           else
             ! material properties vary, we need to average cell centered values
            eps00    = epsv(m00)
            eps10    = epsv(m10)
            eps01    = epsv(m01)
            sigmaE00 = sigmaEv(m00)
            sigmaE10 = sigmaEv(m10)
            sigmaE01 = sigmaEv(m01)
            epsEx = .5*( eps00 + eps01 )  ! epsEx = eps(i1+.5,i2) = ave of eps(i1+.5,i2+.5) and eps(i1+.5,i2-.5)
            epsEy = .5*( eps00 + eps10 )
            sigmaEx = .5*( sigmaE00 + sigmaE01 )
            sigmaEy = .5*( sigmaE00 + sigmaE10 )
           end if

        ! Ex(i1ph,i2,i3,te) : edge centered
        dcEx =  1.+sigmaEx*dt/(2.*epsEx)
        caEx = (1.-sigmaEx*dt/(2.*epsEx))/dcEx
        cbyEx = (dt/(epsEx*dx(1)))/dcEx
        cjEx  =  cf*dt/(epsEx*dcEx)
        un(i1,i2,i3,ex) = caEx*u(i1,i2,i3,ex) + cbyEx*( u(i1,i2,i3,hz)-
     & u(i1,i2-1,i3,hz) )- cjEx*f(i1,i2,i3,ex)

        ! Ey(i1,i2ph,i3,te)
        dcEy =  1.+sigmaEy*dt/(2.*epsEy)
        caEy = (1.-sigmaEy*dt/(2.*epsEy))/dcEy
        cbxEy = (dt/(epsEy*dx(0)))/dcEy
        cjEy  =  cf*dt/(epsEy*dcEy)
        un(i1,i2,i3,ey) = caEy*u(i1,i2,i3,ey) - cbxEy*( u(i1,i2,i3,hz)-
     & u(i1-1,i2,i3,hz) )- cjEy*f(i1,i2,i3,ey)

        if( .false. )then
          write(debugFile,'(" mxYee:advE i1,i2=",2i4," old Ex,Ey =",
     & 2e10.2," new=",2e10.2)') i1,i2,u(i1,i2,i3,ex),u(i1,i2,i3,ey),
     & un(i1,i2,i3,ex),un(i1,i2,i3,ey)
          ! ' 
        end if

       end do
       end do
       end do

      else
        ! -- Advance E : 3D ---
       do i3=n3a,n3b
       do i2=n2a,n2b
       do i1=n1a,n1b

         ! Here we assume the material properties are cell centered
          ! ----- 3D ----
          m000 = media(i1  ,i2  ,i3  )  ! pointer to media properties in cell (i1+.5,i2+.5,i3+.5)
          m100 = media(i1-1,i2  ,i3  )  ! pointer to media properties in cell (i1-.5,i2+.5,i3+.5)
          m010 = media(i1  ,i2-1,i3  )  ! pointer to media properties in cell (i1+.5,i2-.5,i3+.5)
          m001 = media(i1  ,i2  ,i3-1)  ! pointer to media properties in cell (i1+.5,i2+.5,i3-.5)
           m110 = media(i1-1,i2-1,i3  )  ! pointer to media properties in cell (i1-.5,i2-.5,i3+.5)
           m101 = media(i1-1,i2  ,i3-1)  ! pointer to media properties in cell (i1-.5,i2+.5,i3-.5)
           m011 = media(i1  ,i2-1,i3-1)  ! pointer to media properties in cell (i1+.5,i2-.5,i3-.5)
           if( m000.eq.m100 .and. m000.eq.m010 .and. m000.eq.m001 
     & .and. m000.eq.m110 .and. m000.eq.m101 .and. m000.eq.m011 )then
             ! locally uniform material properties: 
            epsEx = epsv(m000)
            epsEy = epsEx
            epsEz = epsEx
            sigmaEx = sigmaEv(m000)
            sigmaEy = sigmaEx
            sigmaEz = sigmaEx
           else
             ! material properties vary, we need to average cell centered values
            eps000    = epsv(m000)
            eps100    = epsv(m100)
            eps010    = epsv(m010)
            eps110    = epsv(m110)
            eps001    = epsv(m001)
            eps101    = epsv(m101)
            eps011    = epsv(m011)
            sigmaE000    = sigmaEv(m000)
            sigmaE100    = sigmaEv(m100)
            sigmaE010    = sigmaEv(m010)
            sigmaE110    = sigmaEv(m110)
            sigmaE001    = sigmaEv(m001)
            sigmaE101    = sigmaEv(m101)
            sigmaE011    = sigmaEv(m011)
            epsEx = .25*( eps000 + eps010 + eps001 + eps011 )  ! epsEx = eps(i1+.5,i2,i3) = ave of 4 neighbouring cells
            epsEy = .25*( eps000 + eps100 + eps001 + eps101 )
            epsEz = .25*( eps000 + eps100 + eps010 + eps110 )
            sigmaEx = .25*( sigmaE000 + sigmaE010 + sigmaE001 + 
     & sigmaE011 )
            sigmaEy = .25*( sigmaE000 + sigmaE100 + sigmaE001 + 
     & sigmaE101 )
            sigmaEz = .25*( sigmaE000 + sigmaE100 + sigmaE010 + 
     & sigmaE110 )
           end if

        ! m = media(i1,i2,i3)  ! pointer to media properties
        ! eps = epsv(m)
        ! mu  = muv(m)
        ! sigmaE = sigmaEv(m)

        ! Ex(i1ph,i2,i3,nph) : edge centered
        dcEx =  1.+sigmaEx*dt/(2.*epsEx)
        caEx = (1.-sigmaEx*dt/(2.*epsEx))/dcEx
        cbyEx= (dt/(epsEx*dx(1)))/dcEx
        cbzEx= (dt/(epsEx*dx(2)))/dcEx
        cjEx =  cf*dt/(epsEx*dcEx)
        un(i1,i2,i3,ex) = caEx*u(i1,i2,i3,ex) + cbyEx*( u(i1,i2,i3,hz)-
     & u(i1,i2-1,i3,hz) )- cbzEx*( u(i1,i2,i3,hy)-u(i1,i2,i3-1,hy) )- 
     &  cjEx*f(i1,i2,i3,ex)

        ! Ey(i1,i2ph,i3,nph)
        dcEy =  1.+sigmaEy*dt/(2.*epsEy)
        caEy = (1.-sigmaEy*dt/(2.*epsEy))/dcEy
        cbxEy= (dt/(epsEy*dx(0)))/dcEy
        cbzEy= (dt/(epsEy*dx(2)))/dcEy
        cjEy =  cf*dt/(epsEy*dcEy)
        un(i1,i2,i3,ey) = caEy*u(i1,i2,i3,ey) + cbzEy*( u(i1,i2,i3,hx)-
     & u(i1,i2,i3-1,hx) )- cbxEy*( u(i1,i2,i3,hz)-u(i1-1,i2,i3,hz) )- 
     & cjEy*f(i1,i2,i3,ey)

        ! Ez(i1,i2,i3ph,nph)
        dcEz =  1.+sigmaEz*dt/(2.*epsEz)
        caEz = (1.-sigmaEz*dt/(2.*epsEz))/dcEz
        cbxEz= (dt/(epsEz*dx(0)))/dcEz
        cbyEz= (dt/(epsEz*dx(1)))/dcEz
        cjEz =  cf*dt/(epsEz*dcEz)
        un(i1,i2,i3,ez) = caEz*u(i1,i2,i3,ez) + cbxEz*( u(i1,i2,i3,hy)-
     & u(i1-1,i2,i3,hy) )- cbyEz*( u(i1,i2,i3,hx)-u(i1,i2-1,i3,hx) )- 
     &  cjEz*f(i1,i2,i3,ez)

       end do
       end do
       end do
      end if ! nd



      ! ============== TZ Forcing ===============
      if( forcingOption.eq.twilightZoneForcing )then
      if( nd.eq.2 )then
       ! --- Add TZ Forcing to E : 2D ---
       z0=0.
       do i3=n3a,n3b
       do i2=n2a,n2b
       do i1=n1a,n1b

        ! m = media(i1,i2,i3)  ! pointer to media properties
        ! eps = epsv(m)
        ! mu  = muv(m)
        ! sigmaE = sigmaEv(m)

        x0 = (xa(0)+dx(0)*(i1-nra(0)))
        y0 = (xa(1)+dx(1)*(i2-nra(1)))
        xp = x0 + .5*dx(0)
        yp = y0 + .5*dx(1)

        if( useTwilightZoneMaterials.eq.1 )then
         ! Get eps.mu etc. from TZ :
         call ogDeriv(ep,0,0,0,0,xp,y0,z0,tef,epsc,epsEx)
         call ogDeriv(ep,0,0,0,0,xp,y0,z0,tef,sigmaEc,sigmaEx)

         call ogDeriv(ep,0,0,0,0,x0,yp,z0,tef,epsc,epsEy)
         call ogDeriv(ep,0,0,0,0,x0,yp,z0,tef,sigmaEc,sigmaEy)
        else
           ! Here we assume the material properties are cell centered
            ! ----- 2D ----
            m00 = media(i1  ,i2  ,i3)  ! pointer to media properties in cell (i1+.5,i2+.5)
             m10 = media(i1-1,i2  ,i3)  ! pointer to media properties in cell (i1-.5,i2+.5)
             m01 = media(i1  ,i2-1,i3)  ! pointer to media properties in cell (i1+.5,i2-.5)
             if( m00.eq.m10 .and. m00.eq.m01 )then
               ! locally uniform material properties: 
               epsEx = epsv(m00)
               epsEy = epsEx
               sigmaEx = sigmaEv(m00)
               sigmaEy = sigmaEx
             else
               ! material properties vary, we need to average cell centered values
              eps00    = epsv(m00)
              eps10    = epsv(m10)
              eps01    = epsv(m01)
              sigmaE00 = sigmaEv(m00)
              sigmaE10 = sigmaEv(m10)
              sigmaE01 = sigmaEv(m01)
              epsEx = .5*( eps00 + eps01 )  ! epsEx = eps(i1+.5,i2) = ave of eps(i1+.5,i2+.5) and eps(i1+.5,i2-.5)
              epsEy = .5*( eps00 + eps10 )
              sigmaEx = .5*( sigmaE00 + sigmaE01 )
              sigmaEy = .5*( sigmaE00 + sigmaE10 )
             end if
        end if

        call ogDeriv(ep,0,0,0,0,xp,y0,z0,tef,ex,ue (ex))
        call ogDeriv(ep,1,0,0,0,xp,y0,z0,tef,ex,uet(ex))
        call ogDeriv(ep,0,0,1,0,xp,y0,z0,tef,hz,uey(hz))

        un(i1,i2,i3,ex) = un(i1,i2,i3,ex) + dt*( uet(ex)  - ( uey(hz) -
     & sigmaEx*ue(ex) )/epsEx )

        call ogDeriv(ep,0,0,0,0,x0,yp,z0,tef,ey,ue (ey))
        call ogDeriv(ep,1,0,0,0,x0,yp,z0,tef,ey,uet(ey))
        call ogDeriv(ep,0,1,0,0,x0,yp,z0,tef,hz,uex(hz))

        un(i1,i2,i3,ey) = un(i1,i2,i3,ey) + dt*( uet(ey)  - (-uex(hz) -
     & sigmaEy*ue(ey) )/epsEy )

        ! write(debugFile,'(" mxYee:advE:TZ: i1,i2=",2i4," epsEx,epsEy =",2(1pe10.2))') i1,i2,epsEx,epsEy

       end do
       end do
       end do

      else
        ! -- Add TZ Forcing to E : 3D ---
       do i3=n3a,n3b
       do i2=n2a,n2b
       do i1=n1a,n1b

        ! m = media(i1,i2,i3)  ! pointer to media properties
        ! eps = epsv(m)
        ! mu  = muv(m)
        ! sigmaE = sigmaEv(m)


        x0 = (xa(0)+dx(0)*(i1-nra(0)))
        y0 = (xa(1)+dx(1)*(i2-nra(1)))
        z0 = (xa(2)+dx(2)*(i3-nra(2)))

        xp = x0 + .5*dx(0)
        yp = y0 + .5*dx(1)
        zp = z0 + .5*dx(2)

        if( useTwilightZoneMaterials.eq.1 )then
         ! Get eps.mu etc. from TZ :
         call ogDeriv(ep,0,0,0,0,xp,y0,z0,tef,epsc,epsEx)
         call ogDeriv(ep,0,0,0,0,xp,y0,z0,tef,sigmaEc,sigmaEx)

         call ogDeriv(ep,0,0,0,0,x0,yp,z0,tef,epsc,epsEy)
         call ogDeriv(ep,0,0,0,0,x0,yp,z0,tef,sigmaEc,sigmaEy)

         call ogDeriv(ep,0,0,0,0,x0,y0,zp,tef,epsc,epsEz)
         call ogDeriv(ep,0,0,0,0,x0,y0,zp,tef,sigmaEc,sigmaEz)
        else
          ! Here we assume the material properties are cell centered
           ! ----- 3D ----
           m000 = media(i1  ,i2  ,i3  )  ! pointer to media properties in cell (i1+.5,i2+.5,i3+.5)
           m100 = media(i1-1,i2  ,i3  )  ! pointer to media properties in cell (i1-.5,i2+.5,i3+.5)
           m010 = media(i1  ,i2-1,i3  )  ! pointer to media properties in cell (i1+.5,i2-.5,i3+.5)
           m001 = media(i1  ,i2  ,i3-1)  ! pointer to media properties in cell (i1+.5,i2+.5,i3-.5)
            m110 = media(i1-1,i2-1,i3  )  ! pointer to media properties in cell (i1-.5,i2-.5,i3+.5)
            m101 = media(i1-1,i2  ,i3-1)  ! pointer to media properties in cell (i1-.5,i2+.5,i3-.5)
            m011 = media(i1  ,i2-1,i3-1)  ! pointer to media properties in cell (i1+.5,i2-.5,i3-.5)
            if( m000.eq.m100 .and. m000.eq.m010 .and. m000.eq.m001 
     & .and. m000.eq.m110 .and. m000.eq.m101 .and. m000.eq.m011 )then
              ! locally uniform material properties: 
             epsEx = epsv(m000)
             epsEy = epsEx
             epsEz = epsEx
             sigmaEx = sigmaEv(m000)
             sigmaEy = sigmaEx
             sigmaEz = sigmaEx
            else
              ! material properties vary, we need to average cell centered values
             eps000    = epsv(m000)
             eps100    = epsv(m100)
             eps010    = epsv(m010)
             eps110    = epsv(m110)
             eps001    = epsv(m001)
             eps101    = epsv(m101)
             eps011    = epsv(m011)
             sigmaE000    = sigmaEv(m000)
             sigmaE100    = sigmaEv(m100)
             sigmaE010    = sigmaEv(m010)
             sigmaE110    = sigmaEv(m110)
             sigmaE001    = sigmaEv(m001)
             sigmaE101    = sigmaEv(m101)
             sigmaE011    = sigmaEv(m011)
             epsEx = .25*( eps000 + eps010 + eps001 + eps011 )  ! epsEx = eps(i1+.5,i2,i3) = ave of 4 neighbouring cells
             epsEy = .25*( eps000 + eps100 + eps001 + eps101 )
             epsEz = .25*( eps000 + eps100 + eps010 + eps110 )
             sigmaEx = .25*( sigmaE000 + sigmaE010 + sigmaE001 + 
     & sigmaE011 )
             sigmaEy = .25*( sigmaE000 + sigmaE100 + sigmaE001 + 
     & sigmaE101 )
             sigmaEz = .25*( sigmaE000 + sigmaE100 + sigmaE010 + 
     & sigmaE110 )
            end if
        end if

        ! note: TZ values for Ex equation are centered t (xp,y0,z0)
        call ogDeriv(ep,0,0,0,0,xp,y0,z0,tef,ex,ue (ex))
        call ogDeriv(ep,1,0,0,0,xp,y0,z0,tef,ex,uet(ex))
        call ogDeriv(ep,0,0,1,0,xp,y0,z0,tef,hz,uey(hz))
        call ogDeriv(ep,0,0,0,1,xp,y0,z0,tef,hy,uez(hy))

        un(i1,i2,i3,ex) = un(i1,i2,i3,ex) + dt*( uet(ex)  - ( uey(hz) -
     & uez(hy) -sigmaEx*ue(ex) )/epsEx )

        call ogDeriv(ep,0,0,0,0,x0,yp,z0,tef,ey,ue (ey))
        call ogDeriv(ep,1,0,0,0,x0,yp,z0,tef,ey,uet(ey))
        call ogDeriv(ep,0,0,0,1,x0,yp,z0,tef,hx,uez(hx))
        call ogDeriv(ep,0,1,0,0,x0,yp,z0,tef,hz,uex(hz))

        un(i1,i2,i3,ey) = un(i1,i2,i3,ey) + dt*( uet(ey)  - ( uez(hx) -
     & uex(hz) -sigmaEy*ue(ey) )/epsEy )

        call ogDeriv(ep,0,0,0,0,x0,y0,zp,tef,ez,ue (ez))
        call ogDeriv(ep,1,0,0,0,x0,y0,zp,tef,ez,uet(ez))
        call ogDeriv(ep,0,0,1,0,x0,y0,zp,tef,hx,uey(hx))
        call ogDeriv(ep,0,1,0,0,x0,y0,zp,tef,hy,uex(hy))

        un(i1,i2,i3,ez) = un(i1,i2,i3,ez) + dt*( uet(ez)  - ( uex(hy) -
     & uey(hx) -sigmaEz*ue(ez) )/epsEz )

       end do
       end do
       end do
      end if ! nd
      end if ! TZ



      ! ---- Apply BC's to E (un) at t=tep ------
      extra=0
      numberOfGhostPoints=0
       extra1a=extra
       extra1b=extra
       extra2a=extra
       extra2b=extra
       if( nd.eq.3 )then
         extra3a=extra
         extra3b=extra
       else
         extra3a=0
         extra3b=0
       end if
       if( boundaryCondition(0,0).lt.0 )then
         extra1a=max(0,extra1a) ! over-ride extra=-1 : assign ends in periodic directions
         extra1b=extra1a
       else
         if( boundaryCondition(0,0).eq.0 )then
           extra1a=numberOfGhostPoints  ! include interpolation points since we assign ghost points outside these
         end if
         if( boundaryCondition(1,0).eq.0 )then
           extra1b=numberOfGhostPoints
         end if
       end if
       if( boundaryCondition(0,1).lt.0 )then
        extra2a=max(0,extra2a) ! over-ride extra=-1 : assign ends in periodic directions
        extra2b=extra2a
       else
         if( boundaryCondition(0,1).eq.0 )then
           extra2a=numberOfGhostPoints
         end if
         if( boundaryCondition(1,1).eq.0 )then
           extra2b=numberOfGhostPoints
         end if
       end if
       if(  nd.eq.3 .and. boundaryCondition(0,2).lt.0 )then
        extra3a=max(0,extra3a) ! over-ride extra=-1 : assign ends in periodic directions
        extra3b=extra3a
       else
         if( boundaryCondition(0,2).eq.0 )then
           extra3a=numberOfGhostPoints
         end if
         if( boundaryCondition(1,2).eq.0 )then
           extra3b=numberOfGhostPoints
         end if
       end if
       do axis=0,nd-1
       do side=0,1
         n1a=gridIndexRange(0,0)-extra1a
         n1b=gridIndexRange(1,0)+extra1b
         n2a=gridIndexRange(0,1)-extra2a
         n2b=gridIndexRange(1,1)+extra2b
         n3a=gridIndexRange(0,2)-extra3a
         n3b=gridIndexRange(1,2)+extra3b
         if( axis.eq.0 )then
           n1a=gridIndexRange(side,axis)
           n1b=gridIndexRange(side,axis)
         else if( axis.eq.1 )then
           n2a=gridIndexRange(side,axis)
           n2b=gridIndexRange(side,axis)
         else
           n3a=gridIndexRange(side,axis)
           n3b=gridIndexRange(side,axis)
         end if
         is1=0
         is2=0
         is3=0
         if( axis.eq.0 )then
           is1=1-2*side
         else if( axis.eq.1 )then
           is2=1-2*side
         else if( axis.eq.2 )then
           is3=1-2*side
         else
           stop 5
         end if
         axisp1=mod(axis+1,nd)
         axisp2=mod(axis+2,nd)
      !    ! (js1,js2,js3) used to compute tangential derivatives
      !    js1=0
      !    js2=0
      !    js3=0
      !    if( axisp1.eq.0 )then
      !      js1=1-2*side
      !    else if( axisp1.eq.1 )then
      !      js2=1-2*side
      !    else if( axisp1.eq.2 )then
      !      js3=1-2*side
      !    else
      !      stop 5
      !    end if
      ! 
      !    ! (ks1,ks2,ks3) used to compute second tangential derivative
      !    ks1=0
      !    ks2=0
      !    ks3=0
      !    if( axisp2.eq.0 )then
      !      ks1=1-2*side
      !    else if( axisp2.eq.1 )then
      !      ks2=1-2*side
      !    else if( axisp2.eq.2 )then
      !      ks3=1-2*side
      !    else
      !      stop 5
      !    end if
       bc0 =boundaryCondition(side,axis)

       ! write(*,'(" mxYee:E BC side,axis,bc0 =",3i4)') side,axis,bc0

       et1=ex + axisp1 ! tangential component 1
       et2=ex + axisp2 ! tangential component 2
       if( bc0.eq.perfectElectricalConductor )then
         ! PEC : set tangential components to zero
         if( nd.eq.2 )then
          do i3=n3a,n3b
          do i2=n2a,n2b
          do i1=n1a,n1b
           un(i1,i2,i3,et1)=0.
          end do
          end do
          end do
         else
          do i3=n3a,n3b
          do i2=n2a,n2b
          do i1=n1a,n1b
           un(i1,i2,i3,et1)=0.
           un(i1,i2,i3,et2)=0.
          end do
          end do
          end do
         end if
         if( forcingOption.eq.twilightZoneForcing )then
            if( nd.eq.2 )then
             z0=0.
             do i3=n3a,n3b
             do i2=n2a,n2b
             do i1=n1a,n1b
              xv(0) = (xa(0)+dx(0)*(i1-nra(0)))
              xv(1) = (xa(1)+dx(1)*(i2-nra(1)))
              xv(axisp1) = xv(axisp1) + dx(axisp1)*.5
                call ogDeriv(ep,0,0,0,0,xv(0),xv(1),z0,tep,et1,ue(et1))
              un(i1,i2,i3,et1)=ue(et1)
             end do
             end do
             end do
            else ! 3D
             do i3=n3a,n3b
             do i2=n2a,n2b
             do i1=n1a,n1b
              xv(0) = (xa(0)+dx(0)*(i1-nra(0)))
              xv(1) = (xa(1)+dx(1)*(i2-nra(1)))
              xv(2) = (xa(2)+dx(2)*(i3-nra(2)))
              xv(axisp1) = xv(axisp1) + dx(axisp1)*.5
               call ogDeriv(ep,0,0,0,0,xv(0),xv(1),xv(2),tep,et1,ue(
     & et1))
              xv(axisp1) = xv(axisp1) - dx(axisp1)*.5
              xv(axisp2) = xv(axisp2) + dx(axisp2)*.5
               call ogDeriv(ep,0,0,0,0,xv(0),xv(1),xv(2),tep,et2,ue(
     & et2))
              un(i1,i2,i3,et1)=ue(et1)
              un(i1,i2,i3,et2)=ue(et2)
            end do
            end do
            end do
           end if
         end if
       else if( bc0.eq.symmetryBoundaryCondition )then
         ! Symmetry : tangential components of E are odd symmetry 
         if( nd.eq.2 )then
          do i3=n3a,n3b
          do i2=n2a,n2b
          do i1=n1a,n1b
           un(i1,i2,i3,et1)=0.
          end do
          end do
          end do
         else
          do i3=n3a,n3b
          do i2=n2a,n2b
          do i1=n1a,n1b
           ! *** I don't think these are needed since they are not used by H ??
           ! -- Assign symmetry conditions on tau.E : these are used when advancing tau.H
           ! -- note tau.E components on on the boundary
           !un(i1-is1,i2-is2,i3-is3,et1)=2.*un(i1,i2,i3,et1)-un(i1+is1,i2+is2,i3+is3,et1)
          end do
          end do
          end do
         end if

         if( forcingOption.eq.twilightZoneForcing )then
            if( nd.eq.2 )then
             z0=0.
             do i3=n3a,n3b
             do i2=n2a,n2b
             do i1=n1a,n1b
              xv(0) = (xa(0)+dx(0)*(i1-nra(0)))
              xv(1) = (xa(1)+dx(1)*(i2-nra(1)))
              xv(axisp1) = xv(axisp1) + dx(axisp1)*.5
                call ogDeriv(ep,0,0,0,0,xv(0),xv(1),z0,tep,et1,ue(et1))
              un(i1,i2,i3,et1)=ue(et1)
             end do
             end do
             end do
            else ! 3D
             do i3=n3a,n3b
             do i2=n2a,n2b
             do i1=n1a,n1b
              xv(0) = (xa(0)+dx(0)*(i1-nra(0)))
              xv(1) = (xa(1)+dx(1)*(i2-nra(1)))
              xv(2) = (xa(2)+dx(2)*(i3-nra(2)))
              xv(axisp1) = xv(axisp1) + dx(axisp1)*.5
               call ogDeriv(ep,0,0,0,0,xv(0),xv(1),xv(2),tep,et1,ue(
     & et1))
              xv(axisp1) = xv(axisp1) - dx(axisp1)*.5
              xv(axisp2) = xv(axisp2) + dx(axisp2)*.5
               call ogDeriv(ep,0,0,0,0,xv(0),xv(1),xv(2),tep,et2,ue(
     & et2))
              un(i1,i2,i3,et1)=ue(et1)
              un(i1,i2,i3,et2)=ue(et2)
            end do
            end do
            end do
           end if
         end if

       else if( bc0.eq.dirichlet )then

        if( forcingOption.eq.twilightZoneForcing )then
           if( nd.eq.2 )then
            z0=0.
            do i3=n3a,n3b
            do i2=n2a,n2b
            do i1=n1a,n1b
             xv(0) = (xa(0)+dx(0)*(i1-nra(0)))
             xv(1) = (xa(1)+dx(1)*(i2-nra(1)))
             xv(axisp1) = xv(axisp1) + dx(axisp1)*.5
               call ogDeriv(ep,0,0,0,0,xv(0),xv(1),z0,tep,et1,ue(et1))
             un(i1,i2,i3,et1)=ue(et1)
            end do
            end do
            end do
           else ! 3D
            do i3=n3a,n3b
            do i2=n2a,n2b
            do i1=n1a,n1b
             xv(0) = (xa(0)+dx(0)*(i1-nra(0)))
             xv(1) = (xa(1)+dx(1)*(i2-nra(1)))
             xv(2) = (xa(2)+dx(2)*(i3-nra(2)))
             xv(axisp1) = xv(axisp1) + dx(axisp1)*.5
              call ogDeriv(ep,0,0,0,0,xv(0),xv(1),xv(2),tep,et1,ue(et1)
     & )
             xv(axisp1) = xv(axisp1) - dx(axisp1)*.5
             xv(axisp2) = xv(axisp2) + dx(axisp2)*.5
              call ogDeriv(ep,0,0,0,0,xv(0),xv(1),xv(2),tep,et2,ue(et2)
     & )
             un(i1,i2,i3,et1)=ue(et1)
             un(i1,i2,i3,et2)=ue(et2)
           end do
           end do
           end do
          end if

        else if( initialConditionOption.eq.planeWaveInitialCondition )
     & then
           if( nd.eq.2 )then
            z0=0.
            do i3=n3a,n3b
            do i2=n2a,n2b
            do i1=n1a,n1b
             xv(0) = (xa(0)+dx(0)*(i1-nra(0)))
             xv(1) = (xa(1)+dx(1)*(i2-nra(1)))
             xv(axisp1) = xv(axisp1) + dx(axisp1)*.5
               ue(et1) = sin(twoPi*(kx*(xv(0))+ky*(xv(1))+kz*(z0)-cc*(
     & tep)))*pwc(et1)
             un(i1,i2,i3,et1)=ue(et1)
            end do
            end do
            end do
           else ! 3D
            do i3=n3a,n3b
            do i2=n2a,n2b
            do i1=n1a,n1b
             xv(0) = (xa(0)+dx(0)*(i1-nra(0)))
             xv(1) = (xa(1)+dx(1)*(i2-nra(1)))
             xv(2) = (xa(2)+dx(2)*(i3-nra(2)))
             xv(axisp1) = xv(axisp1) + dx(axisp1)*.5
              ue(et1) = sin(twoPi*(kx*(xv(0))+ky*(xv(1))+kz*(xv(2))-cc*
     & (tep)))*pwc(et1)
             xv(axisp1) = xv(axisp1) - dx(axisp1)*.5
             xv(axisp2) = xv(axisp2) + dx(axisp2)*.5
              ue(et2) = sin(twoPi*(kx*(xv(0))+ky*(xv(1))+kz*(xv(2))-cc*
     & (tep)))*pwc(et2)
             un(i1,i2,i3,et1)=ue(et1)
             un(i1,i2,i3,et2)=ue(et2)
           end do
           end do
           end do
          end if

        else if( initialConditionOption .eq. 
     & planeMaterialInterfaceInitialCondition )then
           if( nd.eq.2 )then
            z0=0.
            do i3=n3a,n3b
            do i2=n2a,n2b
            do i1=n1a,n1b
             xv(0) = (xa(0)+dx(0)*(i1-nra(0)))
             xv(1) = (xa(1)+dx(1)*(i2-nra(1)))
             xv(axisp1) = xv(axisp1) + dx(axisp1)*.5
              x0 = (xa(0)+dx(0)*(i1-nra(0)))
              y0 = (xa(1)+dx(1)*(i2-nra(1)))
              xp = x0 + .5*dx(0)
              yp = y0 + .5*dx(1)
               ! **  We check the location of each component ***
                if( nPMI(0)*(xp-xPMI(0)) + nPMI(1)*(y0-xPMI(1)) .le. 
     & 0. )then
                 ! incident + reflected wave: 
                 ue(ex) = (pmc(0)*cos(twoPi*(pmc(19)*(xp-pmc(28))+pmc(
     & 20)*(y0-pmc(29))+pmc(21)*(z0-pmc(30))-pmc(18)*(tep)))+pmc(1)*
     & cos(twoPi*(pmc(22)*(xp-pmc(28))+pmc(23)*(y0-pmc(29))+pmc(24)*(
     & z0-pmc(30))-pmc(18)*(tep))))
                else
                 ! transmitted wave: 
                 ue(ex) = (pmc(12)*cos(twoPi*(pmc(25)*(xp-pmc(28))+pmc(
     & 26)*(y0-pmc(29))+pmc(27)*(z0-pmc(30))-pmc(18)*(tep))))
                end if
                if( nPMI(0)*(x0-xPMI(0)) + nPMI(1)*(yp-xPMI(1)) .le. 
     & 0. )then
                 ! incident + reflected wave: 
                 ue(ey) = (pmc(2)*cos(twoPi*(pmc(19)*(x0-pmc(28))+pmc(
     & 20)*(yp-pmc(29))+pmc(21)*(z0-pmc(30))-pmc(18)*(tep)))+pmc(3)*
     & cos(twoPi*(pmc(22)*(x0-pmc(28))+pmc(23)*(yp-pmc(29))+pmc(24)*(
     & z0-pmc(30))-pmc(18)*(tep))))
                else
                 ! transmitted wave: 
                 ue(ey) = (pmc(13)*cos(twoPi*(pmc(25)*(x0-pmc(28))+pmc(
     & 26)*(yp-pmc(29))+pmc(27)*(z0-pmc(30))-pmc(18)*(tep))))
                end if
                if( nPMI(0)*(xp-xPMI(0)) + nPMI(1)*(yp-xPMI(1)) .le. 
     & 0. )then
                 ! incident + reflected wave: 
                 ue(hz) = (pmc(10)*cos(twoPi*(pmc(19)*(xp-pmc(28))+pmc(
     & 20)*(yp-pmc(29))+pmc(21)*(z0-pmc(30))-pmc(18)*(thp)))+pmc(11)*
     & cos(twoPi*(pmc(22)*(xp-pmc(28))+pmc(23)*(yp-pmc(29))+pmc(24)*(
     & z0-pmc(30))-pmc(18)*(thp))))
                else
                 ! transmitted wave: 
                 ue(hz) = (pmc(17)*cos(twoPi*(pmc(25)*(xp-pmc(28))+pmc(
     & 26)*(yp-pmc(29))+pmc(27)*(z0-pmc(30))-pmc(18)*(thp))))
                end if
             un(i1,i2,i3,et1)=ue(et1)
            end do
            end do
            end do
           else ! 3D
            do i3=n3a,n3b
            do i2=n2a,n2b
            do i1=n1a,n1b
             xv(0) = (xa(0)+dx(0)*(i1-nra(0)))
             xv(1) = (xa(1)+dx(1)*(i2-nra(1)))
             xv(2) = (xa(2)+dx(2)*(i3-nra(2)))
             xv(axisp1) = xv(axisp1) + dx(axisp1)*.5
              x0 = (xa(0)+dx(0)*(i1-nra(0)))
              y0 = (xa(1)+dx(1)*(i2-nra(1)))
              z0 = (xa(2)+dx(2)*(i3-nra(2)))
              xp = x0 + .5*dx(0)
              yp = y0 + .5*dx(1)
              zp = z0 + .5*dx(2)
               ! **  We check the location of each component ***
                if( nPMI(0)*(xp-xPMI(0)) + nPMI(1)*(y0-xPMI(1)) + nPMI(
     & 2)*(z0-xPMI(2)) .le. 0. )then
                 ue(ex) = (pmc(0)*cos(twoPi*(pmc(19)*(xp-pmc(28))+pmc(
     & 20)*(y0-pmc(29))+pmc(21)*(z0-pmc(30))-pmc(18)*(tep)))+pmc(1)*
     & cos(twoPi*(pmc(22)*(xp-pmc(28))+pmc(23)*(y0-pmc(29))+pmc(24)*(
     & z0-pmc(30))-pmc(18)*(tep))))
                else
                 ue(ex) = (pmc(12)*cos(twoPi*(pmc(25)*(xp-pmc(28))+pmc(
     & 26)*(y0-pmc(29))+pmc(27)*(z0-pmc(30))-pmc(18)*(tep))))
                end if
                if( nPMI(0)*(x0-xPMI(0)) + nPMI(1)*(yp-xPMI(1)) + nPMI(
     & 2)*(z0-xPMI(2)) .le. 0. )then
                 ue(ey) = (pmc(2)*cos(twoPi*(pmc(19)*(x0-pmc(28))+pmc(
     & 20)*(yp-pmc(29))+pmc(21)*(z0-pmc(30))-pmc(18)*(tep)))+pmc(3)*
     & cos(twoPi*(pmc(22)*(x0-pmc(28))+pmc(23)*(yp-pmc(29))+pmc(24)*(
     & z0-pmc(30))-pmc(18)*(tep))))
                else
                 ue(ey) = (pmc(13)*cos(twoPi*(pmc(25)*(x0-pmc(28))+pmc(
     & 26)*(yp-pmc(29))+pmc(27)*(z0-pmc(30))-pmc(18)*(tep))))
                end if
                if( nPMI(0)*(x0-xPMI(0)) + nPMI(1)*(y0-xPMI(1)) + nPMI(
     & 2)*(zp-xPMI(2)) .le. 0. )then
                 ue(ez) = (pmc(4)*cos(twoPi*(pmc(19)*(x0-pmc(28))+pmc(
     & 20)*(y0-pmc(29))+pmc(21)*(zp-pmc(30))-pmc(18)*(tep)))+pmc(5)*
     & cos(twoPi*(pmc(22)*(x0-pmc(28))+pmc(23)*(y0-pmc(29))+pmc(24)*(
     & zp-pmc(30))-pmc(18)*(tep))))
                else
                 ue(ez) = (pmc(14)*cos(twoPi*(pmc(25)*(x0-pmc(28))+pmc(
     & 26)*(y0-pmc(29))+pmc(27)*(zp-pmc(30))-pmc(18)*(tep))))
                end if
                if( nPMI(0)*(x0-xPMI(0)) + nPMI(1)*(yp-xPMI(1)) + nPMI(
     & 2)*(zp-xPMI(2)) .le. 0. )then
                 ue(hx) = (pmc(6)*cos(twoPi*(pmc(19)*(x0-pmc(28))+pmc(
     & 20)*(yp-pmc(29))+pmc(21)*(zp-pmc(30))-pmc(18)*(thp)))+pmc(7)*
     & cos(twoPi*(pmc(22)*(x0-pmc(28))+pmc(23)*(yp-pmc(29))+pmc(24)*(
     & zp-pmc(30))-pmc(18)*(thp))))
                else
                 ue(hx) = (pmc(15)*cos(twoPi*(pmc(25)*(x0-pmc(28))+pmc(
     & 26)*(yp-pmc(29))+pmc(27)*(zp-pmc(30))-pmc(18)*(thp))))
                end if
                if( nPMI(0)*(xp-xPMI(0)) + nPMI(1)*(y0-xPMI(1)) + nPMI(
     & 2)*(zp-xPMI(2)) .le. 0. )then
                 ue(hy) = (pmc(8)*cos(twoPi*(pmc(19)*(xp-pmc(28))+pmc(
     & 20)*(y0-pmc(29))+pmc(21)*(zp-pmc(30))-pmc(18)*(thp)))+pmc(9)*
     & cos(twoPi*(pmc(22)*(xp-pmc(28))+pmc(23)*(y0-pmc(29))+pmc(24)*(
     & zp-pmc(30))-pmc(18)*(thp))))
                else
                 ue(hy) = (pmc(16)*cos(twoPi*(pmc(25)*(xp-pmc(28))+pmc(
     & 26)*(y0-pmc(29))+pmc(27)*(zp-pmc(30))-pmc(18)*(thp))))
                end if
                if( nPMI(0)*(xp-xPMI(0)) + nPMI(1)*(yp-xPMI(1)) + nPMI(
     & 2)*(z0-xPMI(2)) .le. 0. )then
                 ue(hz) = (pmc(10)*cos(twoPi*(pmc(19)*(xp-pmc(28))+pmc(
     & 20)*(yp-pmc(29))+pmc(21)*(z0-pmc(30))-pmc(18)*(thp)))+pmc(11)*
     & cos(twoPi*(pmc(22)*(xp-pmc(28))+pmc(23)*(yp-pmc(29))+pmc(24)*(
     & z0-pmc(30))-pmc(18)*(thp))))
                else
                 ue(hz) = (pmc(17)*cos(twoPi*(pmc(25)*(xp-pmc(28))+pmc(
     & 26)*(yp-pmc(29))+pmc(27)*(z0-pmc(30))-pmc(18)*(thp))))
                end if
             xv(axisp1) = xv(axisp1) - dx(axisp1)*.5
             xv(axisp2) = xv(axisp2) + dx(axisp2)*.5
               ! this is done above
             un(i1,i2,i3,et1)=ue(et1)
             un(i1,i2,i3,et2)=ue(et2)
           end do
           end do
           end do
          end if

        else if( knownSolutionOption .eq. 
     & scatteringFromADiskKnownSolution .or.knownSolutionOption .eq. 
     & scatteringFromADielectricDiskKnownSolution 
     & .or.knownSolutionOption .eq. 
     & scatteringFromASphereKnownSolution .or.knownSolutionOption 
     & .eq. scatteringFromADielectricSphereKnownSolution )then

           if( nd.eq.2 )then
            z0=0.
            do i3=n3a,n3b
            do i2=n2a,n2b
            do i1=n1a,n1b
             xv(0) = (xa(0)+dx(0)*(i1-nra(0)))
             xv(1) = (xa(1)+dx(1)*(i2-nra(1)))
             xv(axisp1) = xv(axisp1) + dx(axisp1)*.5
               ue(et1) = uKnown(i1,i2,i3,et1)*sinte+uKnown(i1,i2,i3,
     & et1+3)*coste
             un(i1,i2,i3,et1)=ue(et1)
            end do
            end do
            end do
           else ! 3D
            do i3=n3a,n3b
            do i2=n2a,n2b
            do i1=n1a,n1b
             xv(0) = (xa(0)+dx(0)*(i1-nra(0)))
             xv(1) = (xa(1)+dx(1)*(i2-nra(1)))
             xv(2) = (xa(2)+dx(2)*(i3-nra(2)))
             xv(axisp1) = xv(axisp1) + dx(axisp1)*.5
               ue(et1) = uKnown(i1,i2,i3,et1)*sinte+uKnown(i1,i2,i3,
     & et1+6)*coste
             xv(axisp1) = xv(axisp1) - dx(axisp1)*.5
             xv(axisp2) = xv(axisp2) + dx(axisp2)*.5
               ue(et2) = uKnown(i1,i2,i3,et2)*sinte+uKnown(i1,i2,i3,
     & et2+6)*coste
             un(i1,i2,i3,et1)=ue(et1)
             un(i1,i2,i3,et2)=ue(et2)
           end do
           end do
           end do
          end if

        end if

       else if( bc0.eq.abcEM2 )then

        ! --absorbing BC -- 

         if( axis.eq.0 )then
             if( nd.eq.2 )then
              ! ***************** we need ghost tangential values ***************
              i1=n1a
              i2=n2a
              i3=n3a
                ! -- apply symmetry BC's on all adjacent faces for now: *NOTE* apply to solution u 
                i1p=i1+is1
                if( boundaryCondition(0,axisp1).gt.0 )then
                  i2=n2a
                  u(i1 ,i2-1,i3,ey)=u(i1 ,i2,i3,ey)
                  u(i1p,i2-1,i3,ey)=u(i1p,i2,i3,ey)
                end if
                if( boundaryCondition(1,axisp1).gt.0 )then
                  i2=n2b-1
                  u(i1 ,i2+1,i3,ey)=u(i1 ,i2,i3,ey)
                  u(i1p,i2+1,i3,ey)=u(i1p,i2,i3,ey)
                end if
              ! --- Now apply the ABC ---
              do i3=n3a,n3b
              do i2=n2a,n2b
              do i1=n1a,n1b
                 m000 = media(i1  ,i2  ,i3)  ! pointer to media properties in cell (i1+.5,i2+.5,i3+.5)
                 c= 1./sqrt(epsv(m000)*muv(m000))
               ctx1 = (c*dt-dx(axis))/(c*dt+dx(axis))
               ctx2 = 2.*dx(axis)/(c*dt+dx(axis))
               ctx3 = (c*dt)**2*dx(axis)/( 2.*dx(axisp1)**2*( c*dt+dx(
     & axis)) )
                i1p=i1+is1
                un(i1,i2,i3,ey)=-um(i1p,i2,i3,ey) + ctx1*( un(i1p,i2,
     & i3,ey)+um(i1,i2,i3,ey) ) +ctx2*( u(i1p,i2,i3,ey)+u(i1,i2,i3,ey)
     &  ) +ctx3*( u(i1p,i2-1,i3,ey)-2.*u(i1p,i2,i3,ey)+u(i1p,i2+1,i3,
     & ey) +u(i1 ,i2-1,i3,ey)-2.*u(i1 ,i2,i3,ey)+u(i1 ,i2+1,i3,ey) )
              end do
              end do
              end do
             else
              ! ---- 3D -----
              ! ***************** we need ghost tangential values ***************
              i1=n1a
              i2=n2a
              i3=n3a
               ! -- apply symmetry BC's on all adjacent faces for now: *NOTE* apply to solution u 
               ! we could define a macro for below that permuted 1->2->3->1 and also permuted the args (i1,i2a,i3)->(i2,i3a,i1)->(i1,i2,i3a)
               i1p=i1+is1
               dir2=1
               do side2=0,1
                is=1-2*side2
                if( boundaryCondition(side2,dir2).gt.0 )then
                 i2=gridIndexRange(side2,dir2)
                 i2a=i2-is-side2
                 i2b=i2-is
                 do i3=n3a,n3b
                  ! symmetry condition in y-direction, normal=even, tangential=odd
                  u(i1,i2a,i3,ey)=u(i1,i2a+is,i3,ey)
                  u(i1p,i2a,i3,ey)=u(i1p,i2a+is,i3,ey)
                  u(i1,i2b,i3,ez)=2.*u(i1,i2b+is,i3,ez)-u(i1,i2b+2*(is)
     & ,i3,ez)
                  u(i1p,i2b,i3,ez)=2.*u(i1p,i2b+is,i3,ez)-u(i1p,i2b+2*(
     & is),i3,ez)
                 end do
                end if
               end do
               dir3=2
               do side3=0,1
                is=1-2*side3
                if( boundaryCondition(side3,dir3).gt.0 )then
                 i3=gridIndexRange(side3,dir3)
                 i3a=i3-is-side3
                 i3b=i3-is
                 do i2=n2a,n2b
                  ! symmetry condition in z-direction, normal=even, tangential=odd
                  u(i1,i2,i3a,ez)=u(i1,i2,i3a+is,ez)
                  u(i1p,i2,i3a,ez)=u(i1p,i2,i3a+is,ez)
                  u(i1,i2,i3b,ey)=2.*u(i1,i2,i3b+is,ey)-u(i1,i2,i3b+2*(
     & is),ey)
                  u(i1p,i2,i3b,ey)=2.*u(i1p,i2,i3b+is,ey)-u(i1p,i2,i3b+
     & 2*(is),ey)
                 end do
                end if
               end do
              ! ---- evaluate the ABC ---------------
              do i3=n3a,n3b
              do i2=n2a,n2b
              do i1=n1a,n1b
                 m000 = media(i1  ,i2  ,i3)  ! pointer to media properties in cell (i1+.5,i2+.5,i3+.5)
                 c= 1./sqrt(epsv(m000)*muv(m000))
               ctx1 = (c*dt-dx(axis))/(c*dt+dx(axis))
               ctx2 = 2.*dx(axis)/(c*dt+dx(axis))
                ctx3 = (c*dt)**2*dx(axis)/( 2.*dx(1)**2*( c*dt+dx(axis)
     & ) )
                ctx4 = (c*dt)**2*dx(axis)/( 2.*dx(2)**2*( c*dt+dx(axis)
     & ) )
                i1p=i1+is1
                un(i1,i2,i3,ey)=-um(i1p,i2,i3,ey) + ctx1*( un(i1p,i2,
     & i3,ey)+um(i1,i2,i3,ey) ) +ctx2*( u(i1p,i2,i3,ey)+u(i1,i2,i3,ey)
     &  ) +ctx3*( u(i1p,i2-1,i3,ey)-2.*u(i1p,i2,i3,ey)+u(i1p,i2+1,i3,
     & ey) +u(i1 ,i2-1,i3,ey)-2.*u(i1 ,i2,i3,ey)+u(i1 ,i2+1,i3,ey) )+
     & ctx4*( u(i1p,i2,i3-1,ey)-2.*u(i1p,i2,i3,ey)+u(i1p,i2,i3+1,ey) +
     & u(i1 ,i2,i3-1,ey)-2.*u(i1 ,i2,i3,ey)+u(i1 ,i2,i3+1,ey) )
                ! this is the same formula as above with ey -> ez
                un(i1,i2,i3,ez)=-um(i1p,i2,i3,ez) + ctx1*( un(i1p,i2,
     & i3,ez)+um(i1,i2,i3,ez) ) +ctx2*( u(i1p,i2,i3,ez)+u(i1,i2,i3,ez)
     &  ) +ctx3*( u(i1p,i2-1,i3,ez)-2.*u(i1p,i2,i3,ez)+u(i1p,i2+1,i3,
     & ez) +u(i1 ,i2-1,i3,ez)-2.*u(i1 ,i2,i3,ez)+u(i1 ,i2+1,i3,ez) )+
     & ctx4*( u(i1p,i2,i3-1,ez)-2.*u(i1p,i2,i3,ez)+u(i1p,i2,i3+1,ez) +
     & u(i1 ,i2,i3-1,ez)-2.*u(i1 ,i2,i3,ez)+u(i1 ,i2,i3+1,ez) )
              end do
              end do
              end do
             end if
         else if( axis.eq.1 )then
             if( nd.eq.2 )then
              ! ***************** we need ghost tangential values ***************
              i1=n1a
              i2=n2a
              i3=n3a
                ! -- symmetry BC's for now:
                i2p=i2+is2
                if( boundaryCondition(0,axisp1).gt.0 )then
                  i1=n1a
                  u(i1-1,i2 ,i3,ex)=u(i1,i2 ,i3,ex)
                  u(i1-1,i2p,i3,ex)=u(i1,i2p,i3,ex)
                end if
                if( boundaryCondition(1,axisp1).gt.0 )then
                  i1=n1b-1
                  u(i1+1,i2 ,i3,ex)=u(i1,i2 ,i3,ex)
                  u(i1+1,i2p,i3,ex)=u(i1,i2p,i3,ex)
                end if
              ! --- Now apply the ABC ---
              do i3=n3a,n3b
              do i2=n2a,n2b
              do i1=n1a,n1b
                 m000 = media(i1  ,i2  ,i3)  ! pointer to media properties in cell (i1+.5,i2+.5,i3+.5)
                 c= 1./sqrt(epsv(m000)*muv(m000))
               ctx1 = (c*dt-dx(axis))/(c*dt+dx(axis))
               ctx2 = 2.*dx(axis)/(c*dt+dx(axis))
               ctx3 = (c*dt)**2*dx(axis)/( 2.*dx(axisp1)**2*( c*dt+dx(
     & axis)) )
                i2p=i2+is2
                un(i1,i2,i3,ex)=-um(i1,i2p,i3,ex) + ctx1*( un(i1,i2p,
     & i3,ex)+um(i1,i2,i3,ex) ) +ctx2*( u(i1,i2p,i3,ex)+u(i1,i2,i3,ex)
     &  ) +ctx3*( u(i1-1,i2p,i3,ex)-2.*u(i1,i2p,i3,ex)+u(i1+1,i2p,i3,
     & ex) +u(i1-1,i2p,i3,ex)-2.*u(i1,i2p,i3,ex)+u(i1+1,i2p,i3,ex) )
              end do
              end do
              end do
             else
              ! ---- 3D -----
              ! ***************** we need ghost tangential values ***************
              i1=n1a
              i2=n2a
              i3=n3a
               i2p=i2+is2
               dir3=2
               do side3=0,1
                is=1-2*side3
                if( boundaryCondition(side3,dir3).gt.0 )then
                 i3=gridIndexRange(side3,dir3)
                 i3a=i3-is-side3
                 i3b=i3-is
                 do i1=n1a,n1b
                  ! symmetry condition in z-direction, normal=even, tangential=odd
                  u(i1,i2,i3a,ez)=u(i1,i2,i3a+is,ez)
                  u(i1,i2p,i3a,ez)=u(i1,i2p,i3a+is,ez)
                  u(i1,i2,i3b,ex)=2.*u(i1,i2,i3b+is,ex)-u(i1,i2,i3b+2*(
     & is),ex)
                  u(i1,i2p,i3b,ex)=2.*u(i1,i2p,i3b+is,ex)-u(i1,i2p,i3b+
     & 2*(is),ex)
                 end do
                end if
               end do
               dir1=0
               do side1=0,1
                is=1-2*side1
                if( boundaryCondition(side1,dir1).gt.0 )then
                 i1=gridIndexRange(side1,dir1)
                 i1a=i1-is-side1
                 i1b=i1-is
                 do i3=n3a,n3b
                  ! symmetry condition in x-direction, normal=even, tangential=odd
                  u(i1a,i2,i3,ex)=u(i1a+is,i2,i3,ex)
                  u(i1a,i2p,i3,ex)=u(i1a+is,i2p,i3,ex)
                  u(i1b,i2,i3,ez)=2.*u(i1b+is,i2,i3,ez)-u(i1b+2*(is),
     & i2,i3,ez)
                  u(i1b,i2p,i3,ez)=2.*u(i1b+is,i2p,i3,ez)-u(i1b+2*(is),
     & i2p,i3,ez)
                 end do
                end if
               end do
              ! ---- evaluate the ABC ---------------
              do i3=n3a,n3b
              do i2=n2a,n2b
              do i1=n1a,n1b
                 m000 = media(i1  ,i2  ,i3)  ! pointer to media properties in cell (i1+.5,i2+.5,i3+.5)
                 c= 1./sqrt(epsv(m000)*muv(m000))
               ctx1 = (c*dt-dx(axis))/(c*dt+dx(axis))
               ctx2 = 2.*dx(axis)/(c*dt+dx(axis))
                ctx3 = (c*dt)**2*dx(axis)/( 2.*dx(0)**2*( c*dt+dx(axis)
     & ) )
                ctx4 = (c*dt)**2*dx(axis)/( 2.*dx(2)**2*( c*dt+dx(axis)
     & ) )
                i2p=i2+is2
                un(i1,i2,i3,ex)=-um(i1,i2p,i3,ex) + ctx1*( un(i1,i2p,
     & i3,ex)+um(i1,i2,i3,ex) ) +ctx2*( u(i1,i2p,i3,ex)+u(i1,i2,i3,ex)
     &  ) +ctx3*( u(i1-1,i2p,i3,ex)-2.*u(i1,i2p,i3,ex)+u(i1+1,i2p,i3,
     & ex) +u(i1-1,i2 ,i3,ex)-2.*u(i1,i2 ,i3,ex)+u(i1+1,i2 ,i3,ex) )+
     & ctx4*( u(i1,i2p,i3-1,ex)-2.*u(i1,i2p,i3,ex)+u(i1,i2p,i3+1,ex) +
     & u(i1,i2 ,i3-1,ex)-2.*u(i1,i2 ,i3,ex)+u(i1,i2 ,i3+1,ex) )
                ! this is the same formula as above with ex -> ez
                un(i1,i2,i3,ez)=-um(i1,i2p,i3,ez) + ctx1*( un(i1,i2p,
     & i3,ez)+um(i1,i2,i3,ez) ) +ctx2*( u(i1,i2p,i3,ez)+u(i1,i2,i3,ez)
     &  ) +ctx3*( u(i1-1,i2p,i3,ez)-2.*u(i1,i2p,i3,ez)+u(i1+1,i2p,i3,
     & ez) +u(i1-1,i2 ,i3,ez)-2.*u(i1,i2 ,i3,ez)+u(i1+1,i2 ,i3,ez) )+
     & ctx4*( u(i1,i2p,i3-1,ez)-2.*u(i1,i2p,i3,ez)+u(i1,i2p,i3+1,ez) +
     & u(i1,i2 ,i3-1,ez)-2.*u(i1,i2 ,i3,ez)+u(i1,i2 ,i3+1,ez) )
              end do
              end do
              end do
             end if
         else
             if( nd.eq.2 )then
              ! ***************** we need ghost tangential values ***************
              i1=n1a
              i2=n2a
              i3=n3a
                ! -- symmetry BC's for now:
                i2p=i2+is2
                if( boundaryCondition(0,axisp1).gt.0 )then
                  i1=n1a
                  u(i1-1,i2 ,i3,ex)=u(i1,i2 ,i3,ex)
                  u(i1-1,i2p,i3,ex)=u(i1,i2p,i3,ex)
                end if
                if( boundaryCondition(1,axisp1).gt.0 )then
                  i1=n1b-1
                  u(i1+1,i2 ,i3,ex)=u(i1,i2 ,i3,ex)
                  u(i1+1,i2p,i3,ex)=u(i1,i2p,i3,ex)
                end if
              ! --- Now apply the ABC ---
              do i3=n3a,n3b
              do i2=n2a,n2b
              do i1=n1a,n1b
                 m000 = media(i1  ,i2  ,i3)  ! pointer to media properties in cell (i1+.5,i2+.5,i3+.5)
                 c= 1./sqrt(epsv(m000)*muv(m000))
               ctx1 = (c*dt-dx(axis))/(c*dt+dx(axis))
               ctx2 = 2.*dx(axis)/(c*dt+dx(axis))
               ctx3 = (c*dt)**2*dx(axis)/( 2.*dx(axisp1)**2*( c*dt+dx(
     & axis)) )
                 ! unknown T
                 stop 8008
              end do
              end do
              end do
             else
              ! ---- 3D -----
              ! ***************** we need ghost tangential values ***************
              i1=n1a
              i2=n2a
              i3=n3a
               i3p=i3+is3
               dir1=0
               do side1=0,1
                is=1-2*side1
                if( boundaryCondition(side1,dir1).gt.0 )then
                 i1=gridIndexRange(side1,dir1)
                 i1a=i1-is-side1
                 i1b=i1-is
                 do i2=n2a,n2b
                  ! symmetry condition in x-direction, normal=even, tangential=odd
                  u(i1a,i2,i3,ex)=u(i1a+is,i2,i3,ex)
                  u(i1a,i2,i3p,ex)=u(i1a+is,i2,i3p,ex)
                  u(i1b,i2,i3,ey)=2.*u(i1b+is,i2,i3,ey)-u(i1b+2*(is),
     & i2,i3,ey)
                  u(i1b,i2,i3p,ey)=2.*u(i1b+is,i2,i3p,ey)-u(i1b+2*(is),
     & i2,i3p,ey)
                 end do
                end if
               end do
               dir2=1
               do side2=0,1
                is=1-2*side2
                if( boundaryCondition(side2,dir2).gt.0 )then
                 i2=gridIndexRange(side2,dir2)
                 i2a=i2-is-side2
                 i2b=i2-is
                 do i1=n1a,n1b
                  ! symmetry condition in y-direction, normal=even, tangential=odd
                  u(i1,i2a,i3,ey)=u(i1,i2a+is,i3,ey)
                  u(i1,i2a,i3p,ey)=u(i1,i2a+is,i3p,ey)
                  u(i1,i2b,i3,ex)=2.*u(i1,i2b+is,i3,ex)-u(i1,i2b+2*(is)
     & ,i3,ex)
                  u(i1,i2b,i3p,ex)=2.*u(i1,i2b+is,i3p,ex)-u(i1,i2b+2*(
     & is),i3p,ex)
                 end do
                end if
               end do
              ! ---- evaluate the ABC ---------------
              do i3=n3a,n3b
              do i2=n2a,n2b
              do i1=n1a,n1b
                 m000 = media(i1  ,i2  ,i3)  ! pointer to media properties in cell (i1+.5,i2+.5,i3+.5)
                 c= 1./sqrt(epsv(m000)*muv(m000))
               ctx1 = (c*dt-dx(axis))/(c*dt+dx(axis))
               ctx2 = 2.*dx(axis)/(c*dt+dx(axis))
                ctx3 = (c*dt)**2*dx(axis)/( 2.*dx(0)**2*( c*dt+dx(axis)
     & ) )
                ctx4 = (c*dt)**2*dx(axis)/( 2.*dx(1)**2*( c*dt+dx(axis)
     & ) )
                i3p=i3+is3
                un(i1,i2,i3,ex)=-um(i1,i2,i3p,ex) + ctx1*( un(i1,i2,
     & i3p,ex)+um(i1,i2,i3,ex) ) +ctx2*( u(i1,i2,i3p,ex)+u(i1,i2,i3,
     & ex) ) +ctx3*( u(i1-1,i2,i3p,ex)-2.*u(i1,i2,i3p,ex)+u(i1+1,i2,
     & i3p,ex) +u(i1-1,i2,i3 ,ex)-2.*u(i1,i2,i3 ,ex)+u(i1+1,i2,i3 ,ex)
     &  )+ctx4*( u(i1,i2-1,i3p,ex)-2.*u(i1,i2,i3p,ex)+u(i1,i2+1,i3p,
     & ex) +u(i1,i2-1,i3 ,ex)-2.*u(i1,i2,i3 ,ex)+u(i1,i2+1,i3 ,ex) )
                ! this is the same formula as above with ey -> ez
                un(i1,i2,i3,ey)=-um(i1,i2,i3p,ey) + ctx1*( un(i1,i2,
     & i3p,ey)+um(i1,i2,i3,ey) ) +ctx2*( u(i1,i2,i3p,ey)+u(i1,i2,i3,
     & ey) ) +ctx3*( u(i1-1,i2,i3p,ey)-2.*u(i1,i2,i3p,ey)+u(i1+1,i2,
     & i3p,ey) +u(i1-1,i2,i3 ,ey)-2.*u(i1,i2,i3 ,ey)+u(i1+1,i2,i3 ,ey)
     &  )+ctx4*( u(i1,i2-1,i3p,ey)-2.*u(i1,i2,i3p,ey)+u(i1,i2+1,i3p,
     & ey) +u(i1,i2-1,i3 ,ey)-2.*u(i1,i2,i3 ,ey)+u(i1,i2+1,i3 ,ey) )
              end do
              end do
              end do
             end if
         end if

!         ! For now just extrapolate in space time ******** fix me ********
!          if( nd.eq.2 )then
!           beginLoops()
!            ! un(i1,i2,i3,et1)=extrap3(un,i1+is1,i2+is2,i3+is3,et1,is1,is2,is3)
!            un(i1,i2,i3,et1)=u(i1+is1,i2+is2,i3+is3,et1)
!           endLoops()
!          else
!           beginLoops()
!            ! un(i1,i2,i3,et1)=extrap3(un,i1+is1,i2+is2,i3+is3,et1,is1,is2,is3)
!            ! un(i1,i2,i3,et2)=extrap3(un,i1+is1,i2+is2,i3+is3,et2,is1,is2,is3)
!            un(i1,i2,i3,et1)=u(i1+is1,i2+is2,i3+is3,et1)
!            un(i1,i2,i3,et2)=u(i1+is1,i2+is2,i3+is3,et2)
!           endLoops()
!          end if

       else if( bc0.gt.0 )then
         write(*,'(" mxYee:ERROR: unknown BC for side,axis,bc0 =",3i4)
     & ') side,axis,bc0
         stop 7701
       end if
       end do
       end do
       ! reset these values
       n1a=gridIndexRange(0,0)
       n1b=gridIndexRange(1,0)
       n2a=gridIndexRange(0,1)
       n2b=gridIndexRange(1,1)
       n3a=gridIndexRange(0,2)
       n3b=gridIndexRange(1,2)


      ! ------------------------------
      ! ---------- Bodies ------------
      ! ------------------------------
      if( maskBodies.eq.1 )then
       if( debug.gt.0 )then
         write(debugFile,'(" mxYee: apply embedded BC to masked bodies,
     &  t=",e10.2)') tep
       end if
       ! ' 
       ! mask(i1,i2,i3) != 0 : this cell is masked out -- assume a perfect conductor

       if( forcingOption.eq.twilightZoneForcing )then
         ! Assign the TZ solution on the stair-step boundaries of bodies
          if( nd.eq.2 )then
           z0=0.
           do i3=n3a,n3b
           do i2=n2a,n2b
           do i1=n1a,n1b
            if( mask(i1,i2,i3).ne.0 )then
              ! this cell in inside the body
              x0 = (xa(0)+dx(0)*(i1-nra(0)))
              y0 = (xa(1)+dx(1)*(i2-nra(1)))
              xh = x0 + .5*dx(0)
              yh = y0 + .5*dx(1)
              xp = x0 + dx(0)
              yp = y0 + dx(1)
             ! this is a bit inefficient since we set most values twice, once for each cell
               call ogDeriv(ep,0,0,0,0,xh,y0,z0,tep,ex,ue00(ex))
               call ogDeriv(ep,0,0,0,0,xh,yp,z0,tep,ex,ue10(ex))
               call ogDeriv(ep,0,0,0,0,x0,yh,z0,tep,ey,ue00(ey))
               call ogDeriv(ep,0,0,0,0,xp,yh,z0,tep,ey,ue10(ey))
             un(i1  ,i2  ,i3,ex)=ue00(ex)
             un(i1  ,i2+1,i3,ex)=ue10(ex)
             un(i1  ,i2  ,i3,ey)=ue00(ey)
             un(i1+1,i2  ,i3,ey)=ue10(ey)
            end if
           end do
           end do
           end do
          else
           ! --- 3D ---
           do i3=n3a,n3b
           do i2=n2a,n2b
           do i1=n1a,n1b
             if( mask(i1,i2,i3).ne.0 )then
              x0 = (xa(0)+dx(0)*(i1-nra(0)))
              y0 = (xa(1)+dx(1)*(i2-nra(1)))
              z0 = (xa(2)+dx(2)*(i3-nra(2)))
              xh = x0 + .5*dx(0)
              yh = y0 + .5*dx(1)
              zh = z0 + .5*dx(2)
              xp = x0 + dx(0)
              yp = y0 + dx(1)
              zp = z0 + dx(2)
               call ogDeriv(ep,0,0,0,0,xh,y0,z0,tep,ex,ue00(ex))
               call ogDeriv(ep,0,0,0,0,xh,yp,z0,tep,ex,ue10(ex))
               call ogDeriv(ep,0,0,0,0,xh,y0,zp,tep,ex,ue01(ex))
               call ogDeriv(ep,0,0,0,0,xh,yp,zp,tep,ex,ue11(ex))
               call ogDeriv(ep,0,0,0,0,x0,yh,z0,tep,ey,ue00(ey))
               call ogDeriv(ep,0,0,0,0,xp,yh,z0,tep,ey,ue10(ey))
               call ogDeriv(ep,0,0,0,0,x0,yh,zp,tep,ey,ue01(ey))
               call ogDeriv(ep,0,0,0,0,xp,yh,zp,tep,ey,ue11(ey))
               call ogDeriv(ep,0,0,0,0,x0,y0,zh,tep,ez,ue00(ez))
               call ogDeriv(ep,0,0,0,0,xp,y0,zh,tep,ez,ue10(ez))
               call ogDeriv(ep,0,0,0,0,x0,yp,zh,tep,ez,ue01(ez))
               call ogDeriv(ep,0,0,0,0,xp,yp,zh,tep,ez,ue11(ez))
              un(i1  ,i2  ,i3  ,ex)=ue00(ex)
              un(i1  ,i2+1,i3  ,ex)=ue10(ex)
              un(i1  ,i2  ,i3+1,ex)=ue01(ex)
              un(i1  ,i2+1,i3+1,ex)=ue11(ex)
              un(i1  ,i2  ,i3  ,ey)=ue00(ey)
              un(i1+1,i2  ,i3  ,ey)=ue10(ey)
              un(i1  ,i2  ,i3+1,ey)=ue01(ey)
              un(i1+1,i2  ,i3+1,ey)=ue11(ey)
              un(i1  ,i2  ,i3  ,ez)=ue00(ez)
              un(i1+1,i2  ,i3  ,ez)=ue10(ez)
              un(i1  ,i2+1,i3  ,ez)=ue01(ez)
              un(i1+1,i2+1,i3  ,ez)=ue11(ez)
             end if
           end do
           end do
           end do
          end if

       else if( knownSolutionOption .eq. 
     & scatteringFromADiskKnownSolution .or.knownSolutionOption .eq. 
     & scatteringFromASphereKnownSolution )then

        if( .true. )then
          ! ----- Assign the plane wave  solution on the stair-step boundaries of bodies
          ! write(debugFile,'("mxYee: set plane wave on stair-step")') 
           if( nd.eq.2 )then
            z0=0.
            do i3=n3a,n3b
            do i2=n2a,n2b
            do i1=n1a,n1b
             if( mask(i1,i2,i3).ne.0 )then
               ! this cell in inside the body
               x0 = (xa(0)+dx(0)*(i1-nra(0)))
               y0 = (xa(1)+dx(1)*(i2-nra(1)))
               xh = x0 + .5*dx(0)
               yh = y0 + .5*dx(1)
               xp = x0 + dx(0)
               yp = y0 + dx(1)
              ! this is a bit inefficient since we set most values twice, once for each cell
                ! note that we use minus here since we have subtracted off the plane wave solution
                ue00(ex) = -sin(twoPi*(kx*(xh)+ky*(y0)+kz*(z0)-cc*(tep)
     & ))*pwc(ex)
                ue10(ex) = -sin(twoPi*(kx*(xh)+ky*(yp)+kz*(z0)-cc*(tep)
     & ))*pwc(ex)
                ue00(ey) = -sin(twoPi*(kx*(x0)+ky*(yh)+kz*(z0)-cc*(tep)
     & ))*pwc(ey)
                ue10(ey) = -sin(twoPi*(kx*(xp)+ky*(yh)+kz*(z0)-cc*(tep)
     & ))*pwc(ey)
               ! write(*,'("planeWave stair-step: ",4e10.2)') ue00(ex),ue10(ex),ue00(ey),ue10(ey) 
              un(i1  ,i2  ,i3,ex)=ue00(ex)
              un(i1  ,i2+1,i3,ex)=ue10(ex)
              un(i1  ,i2  ,i3,ey)=ue00(ey)
              un(i1+1,i2  ,i3,ey)=ue10(ey)
             end if
            end do
            end do
            end do
           else
            ! --- 3D ---
            do i3=n3a,n3b
            do i2=n2a,n2b
            do i1=n1a,n1b
              if( mask(i1,i2,i3).ne.0 )then
               x0 = (xa(0)+dx(0)*(i1-nra(0)))
               y0 = (xa(1)+dx(1)*(i2-nra(1)))
               z0 = (xa(2)+dx(2)*(i3-nra(2)))
               xh = x0 + .5*dx(0)
               yh = y0 + .5*dx(1)
               zh = z0 + .5*dx(2)
               xp = x0 + dx(0)
               yp = y0 + dx(1)
               zp = z0 + dx(2)
                ! note that we use minus here since we have subtracted off the plane wave solution
                ue00(ex) = -sin(twoPi*(kx*(xh)+ky*(y0)+kz*(z0)-cc*(tep)
     & ))*pwc(ex)
                ue10(ex) = -sin(twoPi*(kx*(xh)+ky*(yp)+kz*(z0)-cc*(tep)
     & ))*pwc(ex)
                ue01(ex) = -sin(twoPi*(kx*(xh)+ky*(y0)+kz*(zp)-cc*(tep)
     & ))*pwc(ex)
                ue11(ex) = -sin(twoPi*(kx*(xh)+ky*(yp)+kz*(zp)-cc*(tep)
     & ))*pwc(ex)
                ue00(ey) = -sin(twoPi*(kx*(x0)+ky*(yh)+kz*(z0)-cc*(tep)
     & ))*pwc(ey)
                ue10(ey) = -sin(twoPi*(kx*(xp)+ky*(yh)+kz*(z0)-cc*(tep)
     & ))*pwc(ey)
                ue01(ey) = -sin(twoPi*(kx*(x0)+ky*(yh)+kz*(zp)-cc*(tep)
     & ))*pwc(ey)
                ue11(ey) = -sin(twoPi*(kx*(xp)+ky*(yh)+kz*(zp)-cc*(tep)
     & ))*pwc(ey)
                ue00(ez) = -sin(twoPi*(kx*(x0)+ky*(y0)+kz*(zh)-cc*(tep)
     & ))*pwc(ez)
                ue10(ez) = -sin(twoPi*(kx*(xp)+ky*(y0)+kz*(zh)-cc*(tep)
     & ))*pwc(ez)
                ue01(ez) = -sin(twoPi*(kx*(x0)+ky*(yp)+kz*(zh)-cc*(tep)
     & ))*pwc(ez)
                ue11(ez) = -sin(twoPi*(kx*(xp)+ky*(yp)+kz*(zh)-cc*(tep)
     & ))*pwc(ez)
               un(i1  ,i2  ,i3  ,ex)=ue00(ex)
               un(i1  ,i2+1,i3  ,ex)=ue10(ex)
               un(i1  ,i2  ,i3+1,ex)=ue01(ex)
               un(i1  ,i2+1,i3+1,ex)=ue11(ex)
               un(i1  ,i2  ,i3  ,ey)=ue00(ey)
               un(i1+1,i2  ,i3  ,ey)=ue10(ey)
               un(i1  ,i2  ,i3+1,ey)=ue01(ey)
               un(i1+1,i2  ,i3+1,ey)=ue11(ey)
               un(i1  ,i2  ,i3  ,ez)=ue00(ez)
               un(i1+1,i2  ,i3  ,ez)=ue10(ez)
               un(i1  ,i2+1,i3  ,ez)=ue01(ez)
               un(i1+1,i2+1,i3  ,ez)=ue11(ez)
              end if
            end do
            end do
            end do
           end if
        else
        ! ---- For testing set the boundary values to the known scattered solution ----
           if( nd.eq.2 )then
            z0=0.
            do i3=n3a,n3b
            do i2=n2a,n2b
            do i1=n1a,n1b
             if( mask(i1,i2,i3).ne.0 )then
               ! this cell in inside the body
               x0 = (xa(0)+dx(0)*(i1-nra(0)))
               y0 = (xa(1)+dx(1)*(i2-nra(1)))
               xh = x0 + .5*dx(0)
               yh = y0 + .5*dx(1)
               xp = x0 + dx(0)
               yp = y0 + dx(1)
              ! this is a bit inefficient since we set most values twice, once for each cell
                ue00(ex) = uKnown(i1  ,i2  ,i3,ex)*sinte+uKnown(i1  ,
     & i2  ,i3,ex+3)*coste
                ue10(ex) = uKnown(i1  ,i2+1,i3,ex)*sinte+uKnown(i1  ,
     & i2+1,i3,ex+3)*coste
                ue00(ey) = uKnown(i1  ,i2  ,i3,ey)*sinte+uKnown(i1  ,
     & i2  ,i3,ey+3)*coste
                ue10(ey) = uKnown(i1+1,i2  ,i3,ey)*sinte+uKnown(i1+1,
     & i2  ,i3,ey+3)*coste
              un(i1  ,i2  ,i3,ex)=ue00(ex)
              un(i1  ,i2+1,i3,ex)=ue10(ex)
              un(i1  ,i2  ,i3,ey)=ue00(ey)
              un(i1+1,i2  ,i3,ey)=ue10(ey)
             end if
            end do
            end do
            end do
           else
            ! --- 3D ---
            do i3=n3a,n3b
            do i2=n2a,n2b
            do i1=n1a,n1b
              if( mask(i1,i2,i3).ne.0 )then
               x0 = (xa(0)+dx(0)*(i1-nra(0)))
               y0 = (xa(1)+dx(1)*(i2-nra(1)))
               z0 = (xa(2)+dx(2)*(i3-nra(2)))
               xh = x0 + .5*dx(0)
               yh = y0 + .5*dx(1)
               zh = z0 + .5*dx(2)
               xp = x0 + dx(0)
               yp = y0 + dx(1)
               zp = z0 + dx(2)
                stop 1164 ! -- finish me ---
                ue00(ex)=uKnown(i1  ,i2  ,i3  ,ex)*sinte+uKnown(i1  ,
     & i2  ,i3  ,ex+6)*coste
                ue10(ex)=uKnown(i1  ,i2+1,i3  ,ex)*sinte+uKnown(i1  ,
     & i2+1,i3  ,ex+6)*coste
                ue01(ex)=uKnown(i1  ,i2  ,i3+1,ex)*sinte+uKnown(i1  ,
     & i2  ,i3+1,ex+6)*coste
                ue11(ex)=uKnown(i1  ,i2+1,i3+1,ex)*sinte+uKnown(i1  ,
     & i2+1,i3+1,ex+6)*coste
                ue00(ey)=uKnown(i1  ,i2  ,i3  ,ey)*sinte+uKnown(i1  ,
     & i2  ,i3  ,ey+6)*coste
                ue10(ey)=uKnown(i1  ,i2+1,i3  ,ey)*sinte+uKnown(i1  ,
     & i2+1,i3  ,ey+6)*coste
                ue01(ey)=uKnown(i1  ,i2  ,i3+1,ey)*sinte+uKnown(i1  ,
     & i2  ,i3+1,ey+6)*coste
                ue11(ey)=uKnown(i1  ,i2+1,i3+1,ey)*sinte+uKnown(i1  ,
     & i2+1,i3+1,ey+6)*coste
                ue00(ez)=uKnown(i1  ,i2  ,i3  ,ez)*sinte+uKnown(i1  ,
     & i2  ,i3  ,ez+6)*coste
                ue10(ez)=uKnown(i1  ,i2+1,i3  ,ez)*sinte+uKnown(i1  ,
     & i2+1,i3  ,ez+6)*coste
                ue01(ez)=uKnown(i1  ,i2  ,i3+1,ez)*sinte+uKnown(i1  ,
     & i2  ,i3+1,ez+6)*coste
                ue11(ez)=uKnown(i1  ,i2+1,i3+1,ez)*sinte+uKnown(i1  ,
     & i2+1,i3+1,ez+6)*coste
               un(i1  ,i2  ,i3  ,ex)=ue00(ex)
               un(i1  ,i2+1,i3  ,ex)=ue10(ex)
               un(i1  ,i2  ,i3+1,ex)=ue01(ex)
               un(i1  ,i2+1,i3+1,ex)=ue11(ex)
               un(i1  ,i2  ,i3  ,ey)=ue00(ey)
               un(i1+1,i2  ,i3  ,ey)=ue10(ey)
               un(i1  ,i2  ,i3+1,ey)=ue01(ey)
               un(i1+1,i2  ,i3+1,ey)=ue11(ey)
               un(i1  ,i2  ,i3  ,ez)=ue00(ez)
               un(i1+1,i2  ,i3  ,ez)=ue10(ez)
               un(i1  ,i2+1,i3  ,ez)=ue01(ez)
               un(i1+1,i2+1,i3  ,ez)=ue11(ez)
              end if
            end do
            end do
            end do
           end if
        end if

       else

        ! ----- PEC body  -----

        if( nd.eq.2 )then
         do i3=n3a,n3b
         do i2=n2a,n2b
         do i1=n1a,n1b
           if( mask(i1,i2,i3).ne.0 )then
             un(i1  ,i2  ,i3,ex)=0.
             un(i1  ,i2+1,i3,ex)=0.
             un(i1  ,i2  ,i3,ey)=0.
             un(i1+1,i2  ,i3,ey)=0.
           end if
         end do
         end do
         end do
        else
         do i3=n3a,n3b
         do i2=n2a,n2b
         do i1=n1a,n1b
           if( mask(i1,i2,i3).ne.0 )then
             un(i1  ,i2  ,i3  ,ex)=0.
             un(i1  ,i2+1,i3  ,ex)=0.
             un(i1  ,i2  ,i3+1,ex)=0.
             un(i1  ,i2+1,i3+1,ex)=0.

             un(i1  ,i2  ,i3  ,ey)=0.
             un(i1+1,i2  ,i3  ,ey)=0.
             un(i1  ,i2  ,i3+1,ey)=0.
             un(i1+1,i2  ,i3+1,ey)=0.

             un(i1  ,i2  ,i3  ,ez)=0.
             un(i1+1,i2  ,i3  ,ez)=0.
             un(i1  ,i2+1,i3  ,ez)=0.
             un(i1+1,i2+1,i3  ,ez)=0.

           end if
         end do
         end do
         end do
        end if

       end if

      end if


      if( debug.ge.16 .and. forcingOption.eq.twilightZoneForcing )then
      if( nd.eq.2 )then
       ! --- print errors ---
       write(debugFile,'(/," mxYee: After update E, t=",e10.2)') tep
       z0=0.
       do i3=n3a,n3b
       do i2=n2a,n2b
       do i1=n1a,n1b
        x0 = (xa(0)+dx(0)*(i1-nra(0)))
        y0 = (xa(1)+dx(1)*(i2-nra(1)))
        xp = x0 + .5*dx(0)
        yp = y0 + .5*dx(1)

        call ogDeriv(ep,0,0,0,0,xp,y0,z0,tep,ex,ue (ex))
        call ogDeriv(ep,0,0,0,0,x0,yp,z0,tep,ey,ue (ey))
        write(debugFile,'(" mxYee: t=",e10.2," i1,i2=",2i4," Ex,Ey =",
     & 2e10.2," err=",2e10.2)') tep,i1,i2,un(i1,i2,i3,ex),un(i1,i2,i3,
     & ey),un(i1,i2,i3,ex)-ue(ex),un(i1,i2,i3,ey)-ue(ey)
          ! ' 
       end do
       end do
       end do
      end if
      end if

      ! -----------------
      ! --- Advance H ---
      ! -----------------
      if( nd.eq.2 )then
        ! --- Advance H : 2D ---
       do i3=n3a,n3b
       do i2=n2a,n2b
       do i1=n1a,n1b

        !m = media(i1,i2,i3)  ! pointer to media properties
        !eps = epsv(m)
        !mu  = muv(m)
        !sigmaH = sigmaHv(m)

         ! Here we assume the material properties are cell centered
          ! ----- 2D ----
          m00 = media(i1  ,i2  ,i3)  ! pointer to media properties in cell (i1+.5,i2+.5)
            ! Hz is cell centered already: 
            muHz  = muv(m00)
            sigmaHz = sigmaHv(m00)

        dcHz =  1.+sigmaHz*dt/(2.*muHz)
        caHz = (1.-sigmaHz*dt/(2.*muHz))/dcHz

        cbxHz= (dt/(muHz*dx(0)))/dcHz
        cbyHz= (dt/(muHz*dx(1)))/dcHz
        cjHz =  cf*dt/(muHz*dcHz)

        ! Hz(i1ph,i2ph,i3,np1) : face centered
        un(i1,i2,i3,hz) = caHz*u(i1,i2,i3,hz) - cbxHz*( un(i1+1,i2,i3,
     & ey)-un(i1,i2,i3,ey) )+ cbyHz*( un(i1,i2+1,i3,ex)-un(i1,i2,i3,
     & ex) )-  cjHz*f(i1,i2,i3,hz)

       end do
       end do
       end do

      else
        ! --- Advance H : 3D ---
       do i3=n3a,n3b
       do i2=n2a,n2b
       do i1=n1a,n1b

        ! m = media(i1,i2,i3)  ! pointer to media properties
        ! eps = epsv(m)
        ! mu  = muv(m)
        ! sigmaH = sigmaHv(m)

         ! Here we assume the material properties are cell centered
          ! ----- 3D ----
          m000 = media(i1  ,i2  ,i3  )  ! pointer to media properties in cell (i1+.5,i2+.5,i3+.5)
          m100 = media(i1-1,i2  ,i3  )  ! pointer to media properties in cell (i1-.5,i2+.5,i3+.5)
          m010 = media(i1  ,i2-1,i3  )  ! pointer to media properties in cell (i1+.5,i2-.5,i3+.5)
          m001 = media(i1  ,i2  ,i3-1)  ! pointer to media properties in cell (i1+.5,i2+.5,i3-.5)
           ! H lives on a face so we only need to average the two adjacent cell values
           if( m000.eq.m100 .and. m000.eq.m010 .and. m000.eq.m001 )then
             ! locally uniform material properties: 
             muHx = muv(m000)
             muHy = muHx
             muHz = muHx
             sigmaHx = sigmaHv(m000)
             sigmaHy = sigmaHx
             sigmaHz = sigmaHx
           else
            mu000 = muv(m000)
            mu100 = muv(m100)
            mu010 = muv(m010)
            mu001 = muv(m001)
            muHx = .5*( mu000 + mu100 )
            muHy = .5*( mu000 + mu010 )
            muHz = .5*( mu000 + mu001 )
            sigmaH000 = sigmaHv(m000)
            sigmaH100 = sigmaHv(m100)
            sigmaH010 = sigmaHv(m010)
            sigmaH001 = sigmaHv(m001)
            sigmaHx = .5*( sigmaH000 + sigmaH100 )
            sigmaHy = .5*( sigmaH000 + sigmaH010 )
            sigmaHz = .5*( sigmaH000 + sigmaH001 )
           end if

        ! Hx(i1,i2ph,i3ph,np1) : face centered
        dcHx =  1.+sigmaHx*dt/(2.*muHx)
        caHx = (1.-sigmaHx*dt/(2.*muHx))/dcHx
        cbyHx= (dt/(muHx*dx(1)))/dcHx
        cbzHx= (dt/(muHx*dx(2)))/dcHx
        cjHx =  cf*dt/(muHx*dcHx)
        un(i1,i2,i3,hx) = caHx*u(i1,i2,i3,hx) - cbyHx*( un(i1,i2+1,i3,
     & ez)-un(i1,i2,i3,ez) )+ cbzHx*( un(i1,i2,i3+1,ey)-un(i1,i2,i3,
     & ey) )-  cjHx*f(i1,i2,i3,hx)

        ! Hy(i1ph,i2,i3ph,np1) : face centered
        dcHy =  1.+sigmaHy*dt/(2.*muHy)
        caHy = (1.-sigmaHy*dt/(2.*muHy))/dcHy
        cbxHy= (dt/(muHy*dx(0)))/dcHy
        cbzHy= (dt/(muHy*dx(2)))/dcHy
        cjHy =  cf*dt/(muHy*dcHy)
        un(i1,i2,i3,hy) = caHy*u(i1,i2,i3,hy) - cbzHy*( un(i1,i2,i3+1,
     & ex)-un(i1,i2,i3,ex) )+ cbxHy*( un(i1+1,i2,i3,ez)-un(i1,i2,i3,
     & ez) )-  cjHy*f(i1,i2,i3,hy)

        ! Hz(i1ph,i2ph,i3,np1) : face centered
        dcHz =  1.+sigmaHz*dt/(2.*muHz)
        caHz = (1.-sigmaHz*dt/(2.*muHz))/dcHz
        cbxHz= (dt/(muHz*dx(0)))/dcHz
        cbyHz= (dt/(muHz*dx(1)))/dcHz
        cjHz =  cf*dt/(muHz*dcHz)
        un(i1,i2,i3,hz) = caHz*u(i1,i2,i3,hz) - cbxHz*( un(i1+1,i2,i3,
     & ey)-un(i1,i2,i3,ey) )+ cbyHz*( un(i1,i2+1,i3,ex)-un(i1,i2,i3,
     & ex) )-  cjHz*f(i1,i2,i3,hz)

       end do
       end do
       end do
      end if ! end nd

      ! ============== TZ Forcing : H , center forcing at t+dt/2 ===========
      if( forcingOption.eq.twilightZoneForcing )then
      if( nd.eq.2 )then
       ! --- Add TZ Forcing to H: 2D ---
       z0=0.
       do i3=n3a,n3b
       do i2=n2a,n2b
       do i1=n1a,n1b

        ! m = media(i1,i2,i3)  ! pointer to media properties
        ! eps = epsv(m)
        ! mu  = muv(m)
        ! sigmaH = sigmaHv(m)

        x0 = (xa(0)+dx(0)*(i1-nra(0)))
        y0 = (xa(1)+dx(1)*(i2-nra(1)))

        xp = x0 + .5*dx(0)
        yp = y0 + .5*dx(1)

        if( useTwilightZoneMaterials.eq.1 )then
         ! Get eps.mu etc. from TZ :
         call ogDeriv(ep,0,0,0,0,xp,yp,z0,thf,muc,muHz)
         call ogDeriv(ep,0,0,0,0,xp,yp,z0,thf,sigmaHc,sigmaHz)
        else
          ! Here we assume the material properties are cell centered
           ! ----- 2D ----
           m00 = media(i1  ,i2  ,i3)  ! pointer to media properties in cell (i1+.5,i2+.5)
             ! Hz is cell centered already: 
             muHz  = muv(m00)
             sigmaHz = sigmaHv(m00)
        end if

        ! TZ values for Hz are evaluated at (xp,yp,z0)
        call ogDeriv(ep,0,0,0,0,xp,yp,z0,thf,hz,ue (hz))
        call ogDeriv(ep,1,0,0,0,xp,yp,z0,thf,hz,uet(hz))
        call ogDeriv(ep,0,1,0,0,xp,yp,z0,thf,ey,uex(ey))
        call ogDeriv(ep,0,0,1,0,xp,yp,z0,thf,ex,uey(ex))

        un(i1,i2,i3,hz) = un(i1,i2,i3,hz) + dt*( uet(hz)  - ( -uex(ey) 
     & + uey(ex) -sigmaHz*ue(hz) )/muHz )

       end do
       end do
       end do

      else
        ! -- Add TZ Forcing to H: 3D ---
       do i3=n3a,n3b
       do i2=n2a,n2b
       do i1=n1a,n1b

        ! m = media(i1,i2,i3)  ! pointer to media properties
        ! eps = epsv(m)
        ! mu  = muv(m)
        ! sigmaH = sigmaHv(m)

        x0 = (xa(0)+dx(0)*(i1-nra(0)))
        y0 = (xa(1)+dx(1)*(i2-nra(1)))
        z0 = (xa(2)+dx(2)*(i3-nra(2)))

        xp = x0 + .5*dx(0)
        yp = y0 + .5*dx(1)
        zp = z0 + .5*dx(2)

        if( useTwilightZoneMaterials.eq.1 )then
         ! Get eps.mu etc. from TZ :
         call ogDeriv(ep,0,0,0,0,x0,yp,zp,thf,muc,muHx)
         call ogDeriv(ep,0,0,0,0,x0,yp,zp,thf,sigmaHc,sigmaHx)
         call ogDeriv(ep,0,0,0,0,xp,y0,zp,thf,muc,muHy)
         call ogDeriv(ep,0,0,0,0,xp,y0,zp,thf,sigmaHc,sigmaHy)
         call ogDeriv(ep,0,0,0,0,xp,yp,z0,thf,muc,muHz)
         call ogDeriv(ep,0,0,0,0,xp,yp,z0,thf,sigmaHc,sigmaHz)
        else
          ! Here we assume the material properties are cell centered
           ! ----- 3D ----
           m000 = media(i1  ,i2  ,i3  )  ! pointer to media properties in cell (i1+.5,i2+.5,i3+.5)
           m100 = media(i1-1,i2  ,i3  )  ! pointer to media properties in cell (i1-.5,i2+.5,i3+.5)
           m010 = media(i1  ,i2-1,i3  )  ! pointer to media properties in cell (i1+.5,i2-.5,i3+.5)
           m001 = media(i1  ,i2  ,i3-1)  ! pointer to media properties in cell (i1+.5,i2+.5,i3-.5)
            ! H lives on a face so we only need to average the two adjacent cell values
            if( m000.eq.m100 .and. m000.eq.m010 .and. m000.eq.m001 )
     & then
              ! locally uniform material properties: 
              muHx = muv(m000)
              muHy = muHx
              muHz = muHx
              sigmaHx = sigmaHv(m000)
              sigmaHy = sigmaHx
              sigmaHz = sigmaHx
            else
             mu000 = muv(m000)
             mu100 = muv(m100)
             mu010 = muv(m010)
             mu001 = muv(m001)
             muHx = .5*( mu000 + mu100 )
             muHy = .5*( mu000 + mu010 )
             muHz = .5*( mu000 + mu001 )
             sigmaH000 = sigmaHv(m000)
             sigmaH100 = sigmaHv(m100)
             sigmaH010 = sigmaHv(m010)
             sigmaH001 = sigmaHv(m001)
             sigmaHx = .5*( sigmaH000 + sigmaH100 )
             sigmaHy = .5*( sigmaH000 + sigmaH010 )
             sigmaHz = .5*( sigmaH000 + sigmaH001 )
            end if
        end if

        call ogDeriv(ep,0,0,0,0,x0,yp,zp,thf,hx,ue (hx))
        call ogDeriv(ep,1,0,0,0,x0,yp,zp,thf,hx,uet(hx))
        call ogDeriv(ep,0,0,1,0,x0,yp,zp,thf,ez,uey(ez))
        call ogDeriv(ep,0,0,0,1,x0,yp,zp,thf,ey,uez(ey))

        un(i1,i2,i3,hx) = un(i1,i2,i3,hx) + dt*( uet(hx)  - (-uey(ez) +
     & uez(ey) -sigmaHx*ue(hx) )/muHx )

        call ogDeriv(ep,0,0,0,0,xp,y0,zp,thf,hy,ue (hy))
        call ogDeriv(ep,1,0,0,0,xp,y0,zp,thf,hy,uet(hy))
        call ogDeriv(ep,0,0,0,1,xp,y0,zp,thf,ex,uez(ex))
        call ogDeriv(ep,0,1,0,0,xp,y0,zp,thf,ez,uex(ez))

        un(i1,i2,i3,hy) = un(i1,i2,i3,hy) + dt*( uet(hy)  - (-uez(ex) +
     & uex(ez) -sigmaHy*ue(hy) )/muHy )

        call ogDeriv(ep,0,0,0,0,xp,yp,z0,thf,hz,ue (hz))
        call ogDeriv(ep,1,0,0,0,xp,yp,z0,thf,hz,uet(hz))
        call ogDeriv(ep,0,0,1,0,xp,yp,z0,thf,ex,uey(ex))
        call ogDeriv(ep,0,1,0,0,xp,yp,z0,thf,ey,uex(ey))

        un(i1,i2,i3,hz) = un(i1,i2,i3,hz) + dt*( uet(hz)  - (-uex(ey) +
     & uey(ex) -sigmaHz*ue(hz) )/muHz )

       end do
       end do
       end do
      end if ! nd
      end if ! TZ



      ! ---- Apply BC's to H (un) at t=thp---------
      extra=0
      numberOfGhostPoints=0
       extra1a=extra
       extra1b=extra
       extra2a=extra
       extra2b=extra
       if( nd.eq.3 )then
         extra3a=extra
         extra3b=extra
       else
         extra3a=0
         extra3b=0
       end if
       if( boundaryCondition(0,0).lt.0 )then
         extra1a=max(0,extra1a) ! over-ride extra=-1 : assign ends in periodic directions
         extra1b=extra1a
       else
         if( boundaryCondition(0,0).eq.0 )then
           extra1a=numberOfGhostPoints  ! include interpolation points since we assign ghost points outside these
         end if
         if( boundaryCondition(1,0).eq.0 )then
           extra1b=numberOfGhostPoints
         end if
       end if
       if( boundaryCondition(0,1).lt.0 )then
        extra2a=max(0,extra2a) ! over-ride extra=-1 : assign ends in periodic directions
        extra2b=extra2a
       else
         if( boundaryCondition(0,1).eq.0 )then
           extra2a=numberOfGhostPoints
         end if
         if( boundaryCondition(1,1).eq.0 )then
           extra2b=numberOfGhostPoints
         end if
       end if
       if(  nd.eq.3 .and. boundaryCondition(0,2).lt.0 )then
        extra3a=max(0,extra3a) ! over-ride extra=-1 : assign ends in periodic directions
        extra3b=extra3a
       else
         if( boundaryCondition(0,2).eq.0 )then
           extra3a=numberOfGhostPoints
         end if
         if( boundaryCondition(1,2).eq.0 )then
           extra3b=numberOfGhostPoints
         end if
       end if
       do axis=0,nd-1
       do side=0,1
         n1a=gridIndexRange(0,0)-extra1a
         n1b=gridIndexRange(1,0)+extra1b
         n2a=gridIndexRange(0,1)-extra2a
         n2b=gridIndexRange(1,1)+extra2b
         n3a=gridIndexRange(0,2)-extra3a
         n3b=gridIndexRange(1,2)+extra3b
         if( axis.eq.0 )then
           n1a=gridIndexRange(side,axis)
           n1b=gridIndexRange(side,axis)
         else if( axis.eq.1 )then
           n2a=gridIndexRange(side,axis)
           n2b=gridIndexRange(side,axis)
         else
           n3a=gridIndexRange(side,axis)
           n3b=gridIndexRange(side,axis)
         end if
         is1=0
         is2=0
         is3=0
         if( axis.eq.0 )then
           is1=1-2*side
         else if( axis.eq.1 )then
           is2=1-2*side
         else if( axis.eq.2 )then
           is3=1-2*side
         else
           stop 5
         end if
         axisp1=mod(axis+1,nd)
         axisp2=mod(axis+2,nd)
      !    ! (js1,js2,js3) used to compute tangential derivatives
      !    js1=0
      !    js2=0
      !    js3=0
      !    if( axisp1.eq.0 )then
      !      js1=1-2*side
      !    else if( axisp1.eq.1 )then
      !      js2=1-2*side
      !    else if( axisp1.eq.2 )then
      !      js3=1-2*side
      !    else
      !      stop 5
      !    end if
      ! 
      !    ! (ks1,ks2,ks3) used to compute second tangential derivative
      !    ks1=0
      !    ks2=0
      !    ks3=0
      !    if( axisp2.eq.0 )then
      !      ks1=1-2*side
      !    else if( axisp2.eq.1 )then
      !      ks2=1-2*side
      !    else if( axisp2.eq.2 )then
      !      ks3=1-2*side
      !    else
      !      stop 5
      !    end if
       bc0 =boundaryCondition(side,axis)

       hn = hx + axis ! normal component of H
       ht1 = hx + axisp1
       ht2 = hx + axisp2

       if( bc0.eq.perfectElectricalConductor )then
         ! PEC : n.H = 0 
         if( nd.eq.2 )then
          do i3=n3a,n3b
          do i2=n2a,n2b
          do i1=n1a,n1b
           un(i1,i2,i3,hn)=0.
          end do
          end do
          end do
         else
          do i3=n3a,n3b
          do i2=n2a,n2b
          do i1=n1a,n1b
           un(i1,i2,i3,hn)=0.
          end do
          end do
          end do
         end if

       else if( bc0.eq.symmetryBoundaryCondition )then
         ! Symmetry : even symmetry on tau.H
         if( nd.eq.2 )then
          do i3=n3a,n3b
          do i2=n2a,n2b
          do i1=n1a,n1b
           ig1 = i1-(is1+1)/2
           ig2 = i2-(is2+1)/2
           un(ig1,ig2,i3,ht1)=un(ig1+is1,ig2+is2,i3,ht1)
          end do
          end do
          end do
         else
          do i3=n3a,n3b
          do i2=n2a,n2b
          do i1=n1a,n1b
           ! -- Assign symmetry conditions on tangent.H : these are used when advancing tangent.E on the boundary

           ! set ghost cell value equal to first cell inside :
           ! Note: left-side:  (is+1)/2 = 1 
           !       right-side: (is+1)/2 = 0 
           ig1 = i1-(is1+1)/2
           ig2 = i2-(is2+1)/2
           ig3 = i3-(is3+1)/2
           un(ig1,ig2,ig3,ht1)=un(ig1+is1,ig2+is2,ig3+is3,ht1)
           un(ig1,ig2,ig3,ht2)=un(ig1+is1,ig2+is2,ig3+is3,ht2)

          end do
          end do
          end do
         end if

       else if( bc0.eq.dirichlet )then
         if( forcingOption.eq.twilightZoneForcing )then
          if( nd.eq.2 )then
            ! nothing to do here
          else
           do i3=n3a,n3b
           do i2=n2a,n2b
           do i1=n1a,n1b

            xv(0) = (xa(0)+dx(0)*(i1-nra(0)))
            xv(1) = (xa(1)+dx(1)*(i2-nra(1)))
            xv(2) = (xa(2)+dx(2)*(i3-nra(2)))

            xv(axisp1)=xv(axisp1)+.5*dx(axisp1)
            xv(axisp2)=xv(axisp2)+.5*dx(axisp2)

            call ogDeriv(ep,0,0,0,0,xv(0),xv(1),xv(2),thp,hn,ue(hn))

            un(i1,i2,i3,hn)=ue(hn)

           end do
           end do
           end do
          end if

         else

         end if

       else if( bc0.eq.abcEM2 )then

        ! --absorbing BC -- 

         if( axis.eq.0 )then
             if( nd.eq.2 )then
               ! nothing to do 
             else
              ! ---- 3D -----
              ! ***************** we need ghost tangential values ***************
              i1=n1a
              i2=n2a
              i3=n3a
                ! -- apply symmetry BC's on all adjacent faces for now: *NOTE* apply to solution u 
                i1p=i1+is1
                 if( boundaryCondition(0,1).gt.0 )then
                  do i3=n3a,n3b
                   u(i1,n2a-1,i3,hn)=u(i1,n2a,i3,hn)
                   u(i1p,n2a-1,i3,hn)=u(i1p,n2a,i3,hn)
                  end do
                 end if
                 if( boundaryCondition(1,1).gt.0 )then
                  do i3=n3a,n3b
                   u(i1,n2b,i3,hn)=u(i1,n2b-1,i3,hn)
                   u(i1p,n2b,i3,hn)=u(i1p,n2b-1,i3,hn)
                  end do
                 end if
                 if( boundaryCondition(0,2).gt.0 )then
                  do i2=n2a,n2b
                   u(i1,i2,n3a-1,hn)=u(i1,i2,n3a,hn)
                   u(i1p,i2,n3a-1,hn)=u(i1p,i2,n3a,hn)
                  end do
                 end if
                 if( boundaryCondition(1,2).gt.0 )then
                  do i2=n2a,n2b
                   u(i1,i2,n3b,hn)=u(i1,i2,n3b-1,hn)
                   u(i1p,i2,n3b,hn)=u(i1p,i2,n3b-1,hn)
                  end do
                 end if
              ! ---- evaluate the ABC ---------------
              do i3=n3a,n3b
              do i2=n2a,n2b
              do i1=n1a,n1b
                 m000 = media(i1  ,i2  ,i3)  ! pointer to media properties in cell (i1+.5,i2+.5,i3+.5)
                 c= 1./sqrt(epsv(m000)*muv(m000))
               ctx1 = (c*dt-dx(axis))/(c*dt+dx(axis))
               ctx2 = 2.*dx(axis)/(c*dt+dx(axis))
                ctx3 = (c*dt)**2*dx(axis)/( 2.*dx(1)**2*( c*dt+dx(axis)
     & ) )
                ctx4 = (c*dt)**2*dx(axis)/( 2.*dx(2)**2*( c*dt+dx(axis)
     & ) )
                i1p=i1+is1
                un(i1,i2,i3,hn)=-um(i1p,i2,i3,hn) + ctx1*( un(i1p,i2,
     & i3,hn)+um(i1,i2,i3,hn) ) +ctx2*( u(i1p,i2,i3,hn)+u(i1,i2,i3,hn)
     &  ) +ctx3*( u(i1p,i2-1,i3,hn)-2.*u(i1p,i2,i3,hn)+u(i1p,i2+1,i3,
     & hn) +u(i1 ,i2-1,i3,hn)-2.*u(i1 ,i2,i3,hn)+u(i1 ,i2+1,i3,hn) )+
     & ctx4*( u(i1p,i2,i3-1,hn)-2.*u(i1p,i2,i3,hn)+u(i1p,i2,i3+1,hn) +
     & u(i1 ,i2,i3-1,hn)-2.*u(i1 ,i2,i3,hn)+u(i1 ,i2,i3+1,hn) )
              end do
              end do
              end do
             end if
         else if( axis.eq.1 )then
             if( nd.eq.2 )then
               ! nothing to do 
             else
              ! ---- 3D -----
              ! ***************** we need ghost tangential values ***************
              i1=n1a
              i2=n2a
              i3=n3a
                i2p=i2+is2
                 if( boundaryCondition(0,2).gt.0 )then
                  do i1=n1a,n1b
                   u(i1,i2,n3a-1,hn)=u(i1,i2,n3a,hn)
                   u(i1,i2p,n3a-1,hn)=u(i1,i2p,n3a,hn)
                  end do
                 end if
                 if( boundaryCondition(1,2).gt.0 )then
                  do i1=n1a,n1b
                   u(i1,i2,n3b,hn)=u(i1,i2,n3b-1,hn)
                   u(i1,i2p,n3b,hn)=u(i1,i2p,n3b-1,hn)
                  end do
                 end if
                 if( boundaryCondition(0,0).gt.0 )then
                  do i3=n3a,n3b
                   u(n1a-1,i2,i3,hn)=u(n1a,i2,i3,hn)
                   u(n1a-1,i2p,i3,hn)=u(n1a,i2p,i3,hn)
                  end do
                 end if
                 if( boundaryCondition(1,0).gt.0 )then
                  do i3=n3a,n3b
                   u(n1b,i2,i3,hn)=u(n1b-1,i2,i3,hn)
                   u(n1b,i2p,i3,hn)=u(n1b-1,i2p,i3,hn)
                  end do
                 end if
              ! ---- evaluate the ABC ---------------
              do i3=n3a,n3b
              do i2=n2a,n2b
              do i1=n1a,n1b
                 m000 = media(i1  ,i2  ,i3)  ! pointer to media properties in cell (i1+.5,i2+.5,i3+.5)
                 c= 1./sqrt(epsv(m000)*muv(m000))
               ctx1 = (c*dt-dx(axis))/(c*dt+dx(axis))
               ctx2 = 2.*dx(axis)/(c*dt+dx(axis))
                ctx3 = (c*dt)**2*dx(axis)/( 2.*dx(0)**2*( c*dt+dx(axis)
     & ) )
                ctx4 = (c*dt)**2*dx(axis)/( 2.*dx(2)**2*( c*dt+dx(axis)
     & ) )
                i2p=i2+is2
                un(i1,i2,i3,hn)=-um(i1,i2p,i3,hn) + ctx1*( un(i1,i2p,
     & i3,hn)+um(i1,i2,i3,hn) ) +ctx2*( u(i1,i2p,i3,hn)+u(i1,i2,i3,hn)
     &  ) +ctx3*( u(i1-1,i2p,i3,hn)-2.*u(i1,i2p,i3,hn)+u(i1+1,i2p,i3,
     & hn) +u(i1-1,i2 ,i3,hn)-2.*u(i1,i2 ,i3,hn)+u(i1+1,i2 ,i3,hn) )+
     & ctx4*( u(i1,i2p,i3-1,hn)-2.*u(i1,i2p,i3,hn)+u(i1,i2p,i3+1,hn) +
     & u(i1,i2 ,i3-1,hn)-2.*u(i1,i2 ,i3,hn)+u(i1,i2 ,i3+1,hn) )
              end do
              end do
              end do
             end if
         else
             if( nd.eq.2 )then
               ! nothing to do 
             else
              ! ---- 3D -----
              ! ***************** we need ghost tangential values ***************
              i1=n1a
              i2=n2a
              i3=n3a
                i3p=i3+is3
                 if( boundaryCondition(0,0).gt.0 )then
                  do i2=n2a,n2b
                   u(n1a-1,i2,i3,hn)=u(n1a,i2,i3,hn)
                   u(n1a-1,i2,i3p,hn)=u(n1a,i2,i3p,hn)
                  end do
                 end if
                 if( boundaryCondition(1,0).gt.0 )then
                  do i2=n2a,n2b
                   u(n1b,i2,i3,hn)=u(n1b-1,i2,i3,hn)
                   u(n1b,i2,i3p,hn)=u(n1b-1,i2,i3p,hn)
                  end do
                 end if
                 if( boundaryCondition(0,1).gt.0 )then
                  do i1=n1a,n1b
                   u(i1,n2a-1,i3,hn)=u(i1,n2a,i3,hn)
                   u(i1,n2a-1,i3p,hn)=u(i1,n2a,i3p,hn)
                  end do
                 end if
                 if( boundaryCondition(1,1).gt.0 )then
                  do i1=n1a,n1b
                   u(i1,n2b,i3,hn)=u(i1,n2b-1,i3,hn)
                   u(i1,n2b,i3p,hn)=u(i1,n2b-1,i3p,hn)
                  end do
                 end if
              ! ---- evaluate the ABC ---------------
              do i3=n3a,n3b
              do i2=n2a,n2b
              do i1=n1a,n1b
                 m000 = media(i1  ,i2  ,i3)  ! pointer to media properties in cell (i1+.5,i2+.5,i3+.5)
                 c= 1./sqrt(epsv(m000)*muv(m000))
               ctx1 = (c*dt-dx(axis))/(c*dt+dx(axis))
               ctx2 = 2.*dx(axis)/(c*dt+dx(axis))
                ctx3 = (c*dt)**2*dx(axis)/( 2.*dx(0)**2*( c*dt+dx(axis)
     & ) )
                ctx4 = (c*dt)**2*dx(axis)/( 2.*dx(1)**2*( c*dt+dx(axis)
     & ) )
                i3p=i3+is3
                un(i1,i2,i3,hn)=-um(i1,i2,i3p,hn) + ctx1*( un(i1,i2,
     & i3p,hn)+um(i1,i2,i3,hn) ) +ctx2*( u(i1,i2,i3p,hn)+u(i1,i2,i3,
     & hn) ) +ctx3*( u(i1-1,i2,i3p,hn)-2.*u(i1,i2,i3p,hn)+u(i1+1,i2,
     & i3p,hn) +u(i1-1,i2,i3 ,hn)-2.*u(i1,i2,i3 ,hn)+u(i1+1,i2,i3 ,hn)
     &  )+ctx4*( u(i1,i2-1,i3p,hn)-2.*u(i1,i2,i3p,hn)+u(i1,i2+1,i3p,
     & hn) +u(i1,i2-1,i3 ,hn)-2.*u(i1,i2,i3 ,hn)+u(i1,i2+1,i3 ,hn) )
              end do
              end do
              end do
             end if
         end if

       end if
       end do
       end do
       ! reset these values
       n1a=gridIndexRange(0,0)
       n1b=gridIndexRange(1,0)
       n2a=gridIndexRange(0,1)
       n2b=gridIndexRange(1,1)
       n3a=gridIndexRange(0,2)
       n3b=gridIndexRange(1,2)

      if( debug.ge.16 .and. forcingOption.eq.twilightZoneForcing )then
      if( nd.eq.2 )then
       ! --- print errors ---
       write(debugFile,'(/," mxYee: After update H, t=",e10.2)') thp
       z0=0.
       do i3=n3a,n3b
       do i2=n2a,n2b
       do i1=n1a,n1b
        x0 = (xa(0)+dx(0)*(i1-nra(0)))
        y0 = (xa(1)+dx(1)*(i2-nra(1)))
        xp = x0 + .5*dx(0)
        yp = y0 + .5*dx(1)

        call ogDeriv(ep,0,0,0,0,xp,yp,z0,thp,hz,ue (hz))
        write(debugFile,'(" mxYee: t=",e10.2," i1,i2=",2i4," Hz=",
     & e10.2," err=",e10.2)') tep,i1,i2,un(i1,i2,i3,hz),un(i1,i2,i3,
     & hz)-ue(hz)
          ! ' 
       end do
       end do
       end do
      end if
      end if

      return
      end





! =====================================================================================
! Macro to compute the errors in the Yee scheme
!  OPTION : TZ, planeWave, planeScat, planeMat
! =====================================================================================





      subroutine mxYeeIcErr( nd, nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,
     & gridIndexRange, u, v, mask, uKnown, boundaryCondition, ipar, 
     & rpar, ierr )
c ===================================================================================
c
c   ---- Initial Conditions and errors for the Yee Approximation for Cartesian Grids ---
c
c
c  option = ipar[0] : 
c       option = 0 : initial conditions
c              = 1 : compute errors 
c              = 2 : compute div(E) in 2d and 3d and div(H) in 3d
c              = 3 : extrapolate u and fill in v with node centered values of the field
c              = 4 : compute  v(nCurlE) = curl( Eknown )/( mu*omega ) 
c  uKnown : holds known solution for some cases
c 
c Output:
c  rpar[20+i]=err(i)  i=ex,ey,...
c
c option==2:
c     divEMax =rpar[30];
c     divHMax =rpar[31];
c     gradEMax=rpar[32]; 
c     gradHMax=rpar[33]; 
c option==4 
c 
c ===================================================================================

      implicit none

      integer nd, nd1a,nd1b,nd2a,nd2b,nd3a,nd3b, ndf1a,ndf1b,ndf2a,
     & ndf2b,ndf3a,ndf3b,n1a,n1b,n2a,n2b,n3a,n3b, ndc, bc,ierr

      real u(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
      real v(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
      real uKnown(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
      integer mask(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
      integer gridIndexRange(0:1,0:2)

      integer ipar(0:*),boundaryCondition(0:1,0:2)
      real rpar(0:*)

c     --- local variables ----

      integer option, forcingOption,initialConditionOption,
     & numberOfComponents,knownSolutionOption,maskBodies
      integer side,axis,ex,ey,ez,hx,hy,hz,grid,debug,side1,side2,side3
      real dx(0:2),t,ep,dt,c,divEMax,divHMax,gradEMax,gradHMax
      integer axisp1,axisp2,i1,i2,i3,i,is1,is2,is3
      integer extra,extra1a,extra1b,extra2a,extra2b,extra3a,extra3b,
     & numberOfGhostPoints
      real eps,mu,kx,ky,kz,slowStartInterval,twoPi,cc

      integer nra(0:2)
      real xv(0:2),xa(0:2),pwc(0:5)
      real ue(0:9),uet(0:9),uex(0:9),uey(0:9),uez(0:9), err(0:5), 
     & uNorm(0:5)
c      real sigmaE,sigmaH,dc,ca,cbx,cby,cbz,cj
      real x0,y0,z0,xp,yp,zp,te,th
      integer en,et1,et2,hn,bc0,ed,hd,nCurlE,n0,n1,n2,n3,n4,n5,ex0,ey0,
     & ez0
      real coste,sinte,costh,sinth,fact

      ! for plane material interfaces:
      integer numberOfPMC
      real pmc(0:40)  !
      real xPMI(0:2),nPMI(0:2)

      integer myid
      integer debugFile
      character*20 debugFileName

      ! boundary condition parameters
! define BC parameters for fortran routines
! boundary conditions
      integer dirichlet,perfectElectricalConductor,
     & perfectMagneticConductor,planeWaveBoundaryCondition,
     & interfaceBC,symmetryBoundaryCondition,abcEM2,abcPML,abc3,abc4,
     & abc5,rbcNonLocal,rbcLocal,lastBC
      parameter( dirichlet=1,perfectElectricalConductor=2,
     & perfectMagneticConductor=3,planeWaveBoundaryCondition=4,
     & symmetryBoundaryCondition=5,interfaceBC=6,abcEM2=7,abcPML=8,
     & abc3=9,abc4=10,abc5=11,rbcNonLocal=12,rbcLocal=13,lastBC=13 )

      ! initial condition parameters
!
! define initial condition options (see Maxwell.h :  enum InitialConditionEnum)
!
      integer defaultInitialCondition,planeWaveInitialCondition,
     & gaussianPlaneWave,gaussianPulseInitialCondition,
     & squareEigenfunctionInitialCondition, 
     & annulusEigenfunctionInitialCondition,zeroInitialCondition,
     & planeWaveScatteredFieldInitialCondition,
     & planeMaterialInterfaceInitialCondition,
     & gaussianIntegralInitialCondition,twilightZoneInitialConditions,
     & userDefinedInitialConditionOption

      parameter( defaultInitialCondition              =0,      
     & planeWaveInitialCondition                =1,      
     & gaussianPlaneWave                        =2,      
     & gaussianPulseInitialCondition            =3,      
     & squareEigenfunctionInitialCondition      =4,      
     & annulusEigenfunctionInitialCondition     =5,      
     & zeroInitialCondition                     =6,      
     & planeWaveScatteredFieldInitialCondition  =7,      
     & planeMaterialInterfaceInitialCondition   =8,      
     & gaussianIntegralInitialCondition         =9,   
     & twilightZoneInitialConditions        =10,
     & userDefinedInitialConditionOption        =11)

      ! forcing options
      ! forcingOptions -- these should match ForcingEnum in Maxwell.h 
      integer noForcing,magneticSinusoidalPointSource,gaussianSource,
     & twilightZoneForcing, gaussianChargeSource, 
     & userDefinedForcingOption
      integer noBoundaryForcing,planeWaveBoundaryForcing,
     & chirpedPlaneWaveBoundaryForcing
      parameter(noForcing                =0,
     & magneticSinusoidalPointSource =1,gaussianSource                
     & =2,twilightZoneForcing           =3,    gaussianChargeSource   
     &        =4,userDefinedForcingOption      =5 )
      ! boundary forcing options when solved directly for the scattered field:
      parameter( noBoundaryForcing              =0,   
     & planeWaveBoundaryForcing       =1,
     & chirpedPlaneWaveBoundaryForcing=2 )

      ! known solutions

      integer noKnownSolution,twilightZoneKnownSolution,
     & planeWaveKnownSolution,gaussianPlaneWaveKnownSolution,
     & gaussianIntegralKnownSolution,
     & planeMaterialInterfaceKnownSolution,
     & scatteringFromADiskKnownSolution,
     & scatteringFromADielectricDiskKnownSolution,
     & scatteringFromASphereKnownSolution,
     & scatteringFromADielectricSphereKnownSolution,
     & squareEigenfunctionKnownSolution,
     & annulusEigenfunctionKnownSolution,
     & eigenfunctionsOfASphereKnownSolution,userDefinedKnownSolution

      parameter( noKnownSolution=0,twilightZoneKnownSolution=1,
     & planeWaveKnownSolution=2,gaussianPlaneWaveKnownSolution=3,
     & gaussianIntegralKnownSolution=4,
     & planeMaterialInterfaceKnownSolution=5,
     & scatteringFromADiskKnownSolution=6,
     & scatteringFromADielectricDiskKnownSolution=7,
     & scatteringFromASphereKnownSolution=8,
     & scatteringFromADielectricSphereKnownSolution=9,
     & squareEigenfunctionKnownSolution=10,
     & annulusEigenfunctionKnownSolution=11,
     & eigenfunctionsOfASphereKnownSolution=12,
     & userDefinedKnownSolution=13 )

      integer rectangular,curvilinear
      parameter(rectangular=0,curvilinear=1)

c     --- start statement function ----
      integer kd,m,n

c............... end statement functions
      save debugFile
      data debugFile/-1/

      ierr=0

      option                =ipar( 0)
      ex                    =ipar( 1)
      ey                    =ipar( 2)
      ez                    =ipar( 3)
      hx                    =ipar( 4)
      hy                    =ipar( 5)
      hz                    =ipar( 6)
      grid                  =ipar( 7)
      debug                 =ipar( 8)
      initialConditionOption=ipar( 9)
      forcingOption         =ipar(10)
      nra(0)                =ipar(11)  ! for computing x,y,z coordinates
      nra(1)                =ipar(12)
      nra(2)                =ipar(13)
      ed                    =ipar(14) ! save div(E) in this component of v for option=2
      hd                    =ipar(15) ! save div(H) in this component of v for option=2
      knownSolutionOption   =ipar(16)
      maskBodies            =ipar(17)
      myid                  =ipar(18)
      numberOfPMC           =ipar(19)

      dx(0)                =rpar(0)
      dx(1)                =rpar(1)
      dx(2)                =rpar(2)
      t                    =rpar(3)
      ep                   =rpar(4)  ! pointer for exact solution
      dt                   =rpar(5)
      c                    =rpar(6)
      kx                   =rpar(7)  ! for plane wave forcing
      ky                   =rpar(8)
      kz                   =rpar(9)
      slowStartInterval    =rpar(10)
      xa(0)                =rpar(11) ! for computing x,y,z coordinates
      xa(1)                =rpar(12)
      xa(2)                =rpar(13)
      pwc(0)               =rpar(14) ! for the plane wave solution
      pwc(1)               =rpar(15)
      pwc(2)               =rpar(16)
      pwc(3)               =rpar(17)
      pwc(4)               =rpar(18)
      pwc(5)               =rpar(19)

      ! plane material interface:
      if( numberOfPMC.ge.40 )then
       ! fix dimension of pmc
       stop 4410
      end if
      do i=0,numberOfPMC-1
        pmc(i)            =rpar(i+20)
      end do
      do i=0,2
        xPMI(i)=pmc(28+i)
        nPMI(i)=pmc(30+i)
      end do
      ! write(*,'(" mxYee: pmc=",(i2,1x,f10.3,2x))') (i,pmc(i),i=0,numberOfPMC-1)

      if( debug.gt.0 .and. debugFile.lt.0 )then
        ! open debug files
        debugFile=10
        if( myid.lt.10 )then
          write(debugFileName,'("mxYee",i1,".fdebug")') myid
        else
          write(debugFileName,'("mxYee",i4.4,".fdebug")') myid
        end if
        ! write(*,*) 'mxYeeIcErr: myid=',myid,' open debug file:',debugFileName
        open(debugFile,file=debugFileName,status='unknown',
     & form='formatted')
        ! '
        ! INQUIRE(FILE=filen, EXIST=filex)
      end if

      if( debug.gt.1 )then
!        write(debugFile,'(" mxYeeIcErr: **START** initialConditionOption,forcingOption=",i4') initialConditionOption,forcingOption
        write(debugFile,1000) initialConditionOption,forcingOption

      end if
 1000 format(" mxYeeIcErr: **START** initialConditionOption,
     & forcingOption=",i4)
      ! for plane wave forcing 
      twoPi=8.*atan2(1.,1.)
      cc= c*sqrt( kx*kx+ky*ky+kz*kz )  ! this should really be the frequency

      if( debug.gt.2 )then
       write(debugFile,'(" mxYeeIcErr: option,ex,ey,ez,hx,hy,hz=",7i4)
     & ') option,ex,ey,ez,hx,hy,hz
       write(debugFile,'(" mxYeeIcErr: t,dt,dx(0),dx(1),dx(2)=",5e10.2)
     & ') t,dt,dx(0),dx(1),dx(2)
       write(debugFile,'(" mxYeeIcErr: xa(0),xa(1),xa(2)=",3e10.2)') 
     & xa(0),xa(1),xa(2)
       write(debugFile,'(" mxYeeIcErr: pwc(0:5)=",6e10.2," cc=",e10.2)
     & ') (pwc(i),i=0,5),cc
       write(debugFile,'(" mxYeeIcErr: nra(0),nra(1),nra(2)=",3i6)') 
     & nra(0),nra(1),nra(2)
       write(debugFile,'(" mxYeeIcErr: gridIndexRange=",3(2i6,2x))') ((
     & gridIndexRange(side,axis),side=0,1),axis=0,2)
       ! '
      end if

!      epsX=1.e-30 ! fix this ***


      te = t        ! E lives at this time
      th = t+.5*dt  ! H lives at this time

      ! write(debugFile,'(" mxYeeIcErr: t,dt,cc,c,kx,ky,kz=",7e10.2)') t,dt,cc,c,kx,ky,kz

      n1a=gridIndexRange(0,0)
      n1b=gridIndexRange(1,0)
      n2a=gridIndexRange(0,1)
      n2b=gridIndexRange(1,1)
      n3a=gridIndexRange(0,2)
      n3b=gridIndexRange(1,2)

      numberOfComponents = 3*(nd-1)

      coste = cos(-twoPi*cc*te)
      sinte = sin(-twoPi*cc*te)
      costh = cos(-twoPi*cc*th)
      sinth = sin(-twoPi*cc*th)

      if( option.eq.0 )then
        ! ************************************************************
        ! ***************** Initial Conditions ***********************
        ! ************************************************************

        ! ============== TZ Forcing ===============
        if( initialConditionOption.eq.twilightZoneInitialConditions 
     & .or. forcingOption.eq.twilightZoneForcing )then
         if( nd.eq.2 )then
          z0=0.
          do i3=n3a,n3b
          do i2=n2a,n2b
          do i1=n1a,n1b
           x0 = (xa(0)+dx(0)*(i1-nra(0)))
           y0 = (xa(1)+dx(1)*(i2-nra(1)))
           xp = x0 + .5*dx(0)
           yp = y0 + .5*dx(1)
           call ogDeriv(ep,0,0,0,0,xp,y0,z0,te,ex,ue(ex))
           call ogDeriv(ep,0,0,0,0,x0,yp,z0,te,ey,ue(ey))
           call ogDeriv(ep,0,0,0,0,xp,yp,z0,th,hz,ue(hz))
           u(i1,i2,i3,ex) = ue(ex)
           u(i1,i2,i3,ey) = ue(ey)
           u(i1,i2,i3,hz) = ue(hz)
          end do
          end do
          end do

         else
           ! -- 3D TZ Initial condition ---
          do i3=n3a,n3b
          do i2=n2a,n2b
          do i1=n1a,n1b
           x0 = (xa(0)+dx(0)*(i1-nra(0)))
           y0 = (xa(1)+dx(1)*(i2-nra(1)))
           z0 = (xa(2)+dx(2)*(i3-nra(2)))
           xp = x0 + .5*dx(0)
           yp = y0 + .5*dx(1)
           zp = z0 + .5*dx(2)
           call ogDeriv(ep,0,0,0,0,xp,y0,z0,te,ex,ue(ex))
           call ogDeriv(ep,0,0,0,0,x0,yp,z0,te,ey,ue(ey))
           call ogDeriv(ep,0,0,0,0,x0,y0,zp,te,ez,ue(ez))

           call ogDeriv(ep,0,0,0,0,x0,yp,zp,th,hx,ue(hx))
           call ogDeriv(ep,0,0,0,0,xp,y0,zp,th,hy,ue(hy))
           call ogDeriv(ep,0,0,0,0,xp,yp,z0,th,hz,ue(hz))

           u(i1,i2,i3,ex) = ue(ex)
           u(i1,i2,i3,ey) = ue(ey)
           u(i1,i2,i3,ez) = ue(ez)

           u(i1,i2,i3,hx) = ue(hx)
           u(i1,i2,i3,hy) = ue(hy)
           u(i1,i2,i3,hz) = ue(hz)

       ! write(debugFile,'(" mxYeeIcE: i1,i2,i3=",3i4," Ex,Ey,Ez =",3e10.2)') i1,i2,i3,u(i1,i2,i3,ex),u(i1,i2,i3,ey),u(i1,i2,i3,ez)
          end do
          end do
          end do
        end if ! nd

       else if( initialConditionOption.eq.planeWaveInitialCondition )
     & then
         ! --------- plane wave ----------------------
         if( nd.eq.2 )then
          z0=0.
          do i3=n3a,n3b
          do i2=n2a,n2b
          do i1=n1a,n1b
           x0 = (xa(0)+dx(0)*(i1-nra(0)))
           y0 = (xa(1)+dx(1)*(i2-nra(1)))
           xp = x0 + .5*dx(0)
           yp = y0 + .5*dx(1)

           u(i1,i2,i3,ex) = sin(twoPi*(kx*(xp)+ky*(y0)-cc*(te)))*pwc(0)
           u(i1,i2,i3,ey) = sin(twoPi*(kx*(x0)+ky*(yp)-cc*(te)))*pwc(1)
           u(i1,i2,i3,hz) = sin(twoPi*(kx*(xp)+ky*(yp)-cc*(th)))*pwc(5)
          end do
          end do
          end do
         else
          do i3=n3a,n3b
          do i2=n2a,n2b
          do i1=n1a,n1b
           x0 = (xa(0)+dx(0)*(i1-nra(0)))
           y0 = (xa(1)+dx(1)*(i2-nra(1)))
           z0 = (xa(2)+dx(2)*(i3-nra(2)))
           xp = x0 + .5*dx(0)
           yp = y0 + .5*dx(1)
           zp = z0 + .5*dx(2)

           u(i1,i2,i3,ex) = sin(twoPi*(kx*(xp)+ky*(y0)+kz*(z0)-cc*(te))
     & )*pwc(0)
           u(i1,i2,i3,ey) = sin(twoPi*(kx*(x0)+ky*(yp)+kz*(z0)-cc*(te))
     & )*pwc(1)
           u(i1,i2,i3,ez) = sin(twoPi*(kx*(x0)+ky*(y0)+kz*(zp)-cc*(te))
     & )*pwc(2)

           u(i1,i2,i3,hx) = sin(twoPi*(kx*(x0)+ky*(yp)+kz*(zp)-cc*(th))
     & )*pwc(3)
           u(i1,i2,i3,hy) = sin(twoPi*(kx*(xp)+ky*(y0)+kz*(zp)-cc*(th))
     & )*pwc(4)
           u(i1,i2,i3,hz) = sin(twoPi*(kx*(xp)+ky*(yp)+kz*(z0)-cc*(th))
     & )*pwc(5)

          end do
          end do
          end do
        end if ! nd


       else if( initialConditionOption .eq. 
     & planeWaveScatteredFieldInitialCondition )then

         if( debug.gt.0 )then
           write(debugFile,'("mxYeeIcErr:Set known solution as IC, t,
     & dt=",2e10.2)') t,dt
           write(debugFile,'("mxYeeIcErr: coste,sinte,costh,sinth",
     & 4e10.2)') coste,sinte,costh,sinth
         end if

         ! ' 
         if( nd.eq.2 )then
          do i3=n3a,n3b
          do i2=n2a,n2b
          do i1=n1a,n1b
           u(i1,i2,i3,ex) = uKnown(i1,i2,i3,ex)*sinte+uKnown(i1,i2,i3,
     & ex+3)*coste
           u(i1,i2,i3,ey) = uKnown(i1,i2,i3,ey)*sinte+uKnown(i1,i2,i3,
     & ey+3)*coste
           u(i1,i2,i3,hz) = uKnown(i1,i2,i3,hz)*sinth+uKnown(i1,i2,i3,
     & hz+3)*costh

! -----------------
!           x0 = XV(i1,i2,i3) 
!           y0 = YV(i1,i2,i3) 
!           xp = x0 + .5*dx(0)
!           yp = y0 + .5*dx(1)
!
!      write(debugFile,'("mxYeeIcErr: known,plane=",6e12.4)') u(i1,i2,i3,ex),planeWave2Dex0(xp,y0,te),u(i1,i2,i3,ey),planeWave2Dey0(x0,yp,te),u(i1,i2,i3,hz),planeWave2Dhz0(xp,yp,th)
!     -----------                
          end do
          end do
          end do
         else
          do i3=n3a,n3b
          do i2=n2a,n2b
          do i1=n1a,n1b
           u(i1,i2,i3,ex) = uKnown(i1,i2,i3,ex)*sinte+uKnown(i1,i2,i3,
     & ex+6)*coste
           u(i1,i2,i3,ey) = uKnown(i1,i2,i3,ey)*sinte+uKnown(i1,i2,i3,
     & ey+6)*coste
           u(i1,i2,i3,ez) = uKnown(i1,i2,i3,ez)*sinte+uKnown(i1,i2,i3,
     & ez+6)*coste

           u(i1,i2,i3,hx) = uKnown(i1,i2,i3,hx)*sinth+uKnown(i1,i2,i3,
     & hx+6)*costh
           u(i1,i2,i3,hy) = uKnown(i1,i2,i3,hy)*sinth+uKnown(i1,i2,i3,
     & hy+6)*costh
           u(i1,i2,i3,hz) = uKnown(i1,i2,i3,hz)*sinth+uKnown(i1,i2,i3,
     & hz+6)*costh
          end do
          end do
          end do

         end if


       else if( initialConditionOption .eq. 
     & planeMaterialInterfaceInitialCondition )then


         if( nd.eq.2 )then
          z0=0.
          do i3=n3a,n3b
          do i2=n2a,n2b
          do i1=n1a,n1b
           x0 = (xa(0)+dx(0)*(i1-nra(0)))
           y0 = (xa(1)+dx(1)*(i2-nra(1)))
           xp = x0 + .5*dx(0)
           yp = y0 + .5*dx(1)

            ! **  We check the location of each component ***
             if( nPMI(0)*(xp-xPMI(0)) + nPMI(1)*(y0-xPMI(1)) .le. 0. )
     & then
              ! incident + reflected wave: 
              ue(ex) = (pmc(0)*cos(twoPi*(pmc(19)*(xp-pmc(28))+pmc(20)*
     & (y0-pmc(29))+pmc(21)*(z0-pmc(30))-pmc(18)*(te)))+pmc(1)*cos(
     & twoPi*(pmc(22)*(xp-pmc(28))+pmc(23)*(y0-pmc(29))+pmc(24)*(z0-
     & pmc(30))-pmc(18)*(te))))
             else
              ! transmitted wave: 
              ue(ex) = (pmc(12)*cos(twoPi*(pmc(25)*(xp-pmc(28))+pmc(26)
     & *(y0-pmc(29))+pmc(27)*(z0-pmc(30))-pmc(18)*(te))))
             end if
             if( nPMI(0)*(x0-xPMI(0)) + nPMI(1)*(yp-xPMI(1)) .le. 0. )
     & then
              ! incident + reflected wave: 
              ue(ey) = (pmc(2)*cos(twoPi*(pmc(19)*(x0-pmc(28))+pmc(20)*
     & (yp-pmc(29))+pmc(21)*(z0-pmc(30))-pmc(18)*(te)))+pmc(3)*cos(
     & twoPi*(pmc(22)*(x0-pmc(28))+pmc(23)*(yp-pmc(29))+pmc(24)*(z0-
     & pmc(30))-pmc(18)*(te))))
             else
              ! transmitted wave: 
              ue(ey) = (pmc(13)*cos(twoPi*(pmc(25)*(x0-pmc(28))+pmc(26)
     & *(yp-pmc(29))+pmc(27)*(z0-pmc(30))-pmc(18)*(te))))
             end if
             if( nPMI(0)*(xp-xPMI(0)) + nPMI(1)*(yp-xPMI(1)) .le. 0. )
     & then
              ! incident + reflected wave: 
              ue(hz) = (pmc(10)*cos(twoPi*(pmc(19)*(xp-pmc(28))+pmc(20)
     & *(yp-pmc(29))+pmc(21)*(z0-pmc(30))-pmc(18)*(th)))+pmc(11)*cos(
     & twoPi*(pmc(22)*(xp-pmc(28))+pmc(23)*(yp-pmc(29))+pmc(24)*(z0-
     & pmc(30))-pmc(18)*(th))))
             else
              ! transmitted wave: 
              ue(hz) = (pmc(17)*cos(twoPi*(pmc(25)*(xp-pmc(28))+pmc(26)
     & *(yp-pmc(29))+pmc(27)*(z0-pmc(30))-pmc(18)*(th))))
             end if
           u(i1,i2,i3,ex)=ue(ex)
           u(i1,i2,i3,ey)=ue(ey)
           u(i1,i2,i3,hz)=ue(hz)

          end do
          end do
          end do

         else ! 3D

          do i3=n3a,n3b
          do i2=n2a,n2b
          do i1=n1a,n1b
           x0 = (xa(0)+dx(0)*(i1-nra(0)))
           y0 = (xa(1)+dx(1)*(i2-nra(1)))
           z0 = (xa(2)+dx(2)*(i3-nra(2)))
           xp = x0 + .5*dx(0)
           yp = y0 + .5*dx(1)
           zp = z0 + .5*dx(2)

            ! **  We check the location of each component ***
             if( nPMI(0)*(xp-xPMI(0)) + nPMI(1)*(y0-xPMI(1)) + nPMI(2)*
     & (z0-xPMI(2)) .le. 0. )then
              ue(ex) = (pmc(0)*cos(twoPi*(pmc(19)*(xp-pmc(28))+pmc(20)*
     & (y0-pmc(29))+pmc(21)*(z0-pmc(30))-pmc(18)*(te)))+pmc(1)*cos(
     & twoPi*(pmc(22)*(xp-pmc(28))+pmc(23)*(y0-pmc(29))+pmc(24)*(z0-
     & pmc(30))-pmc(18)*(te))))
             else
              ue(ex) = (pmc(12)*cos(twoPi*(pmc(25)*(xp-pmc(28))+pmc(26)
     & *(y0-pmc(29))+pmc(27)*(z0-pmc(30))-pmc(18)*(te))))
             end if
             if( nPMI(0)*(x0-xPMI(0)) + nPMI(1)*(yp-xPMI(1)) + nPMI(2)*
     & (z0-xPMI(2)) .le. 0. )then
              ue(ey) = (pmc(2)*cos(twoPi*(pmc(19)*(x0-pmc(28))+pmc(20)*
     & (yp-pmc(29))+pmc(21)*(z0-pmc(30))-pmc(18)*(te)))+pmc(3)*cos(
     & twoPi*(pmc(22)*(x0-pmc(28))+pmc(23)*(yp-pmc(29))+pmc(24)*(z0-
     & pmc(30))-pmc(18)*(te))))
             else
              ue(ey) = (pmc(13)*cos(twoPi*(pmc(25)*(x0-pmc(28))+pmc(26)
     & *(yp-pmc(29))+pmc(27)*(z0-pmc(30))-pmc(18)*(te))))
             end if
             if( nPMI(0)*(x0-xPMI(0)) + nPMI(1)*(y0-xPMI(1)) + nPMI(2)*
     & (zp-xPMI(2)) .le. 0. )then
              ue(ez) = (pmc(4)*cos(twoPi*(pmc(19)*(x0-pmc(28))+pmc(20)*
     & (y0-pmc(29))+pmc(21)*(zp-pmc(30))-pmc(18)*(te)))+pmc(5)*cos(
     & twoPi*(pmc(22)*(x0-pmc(28))+pmc(23)*(y0-pmc(29))+pmc(24)*(zp-
     & pmc(30))-pmc(18)*(te))))
             else
              ue(ez) = (pmc(14)*cos(twoPi*(pmc(25)*(x0-pmc(28))+pmc(26)
     & *(y0-pmc(29))+pmc(27)*(zp-pmc(30))-pmc(18)*(te))))
             end if
             if( nPMI(0)*(x0-xPMI(0)) + nPMI(1)*(yp-xPMI(1)) + nPMI(2)*
     & (zp-xPMI(2)) .le. 0. )then
              ue(hx) = (pmc(6)*cos(twoPi*(pmc(19)*(x0-pmc(28))+pmc(20)*
     & (yp-pmc(29))+pmc(21)*(zp-pmc(30))-pmc(18)*(th)))+pmc(7)*cos(
     & twoPi*(pmc(22)*(x0-pmc(28))+pmc(23)*(yp-pmc(29))+pmc(24)*(zp-
     & pmc(30))-pmc(18)*(th))))
             else
              ue(hx) = (pmc(15)*cos(twoPi*(pmc(25)*(x0-pmc(28))+pmc(26)
     & *(yp-pmc(29))+pmc(27)*(zp-pmc(30))-pmc(18)*(th))))
             end if
             if( nPMI(0)*(xp-xPMI(0)) + nPMI(1)*(y0-xPMI(1)) + nPMI(2)*
     & (zp-xPMI(2)) .le. 0. )then
              ue(hy) = (pmc(8)*cos(twoPi*(pmc(19)*(xp-pmc(28))+pmc(20)*
     & (y0-pmc(29))+pmc(21)*(zp-pmc(30))-pmc(18)*(th)))+pmc(9)*cos(
     & twoPi*(pmc(22)*(xp-pmc(28))+pmc(23)*(y0-pmc(29))+pmc(24)*(zp-
     & pmc(30))-pmc(18)*(th))))
             else
              ue(hy) = (pmc(16)*cos(twoPi*(pmc(25)*(xp-pmc(28))+pmc(26)
     & *(y0-pmc(29))+pmc(27)*(zp-pmc(30))-pmc(18)*(th))))
             end if
             if( nPMI(0)*(xp-xPMI(0)) + nPMI(1)*(yp-xPMI(1)) + nPMI(2)*
     & (z0-xPMI(2)) .le. 0. )then
              ue(hz) = (pmc(10)*cos(twoPi*(pmc(19)*(xp-pmc(28))+pmc(20)
     & *(yp-pmc(29))+pmc(21)*(z0-pmc(30))-pmc(18)*(th)))+pmc(11)*cos(
     & twoPi*(pmc(22)*(xp-pmc(28))+pmc(23)*(yp-pmc(29))+pmc(24)*(z0-
     & pmc(30))-pmc(18)*(th))))
             else
              ue(hz) = (pmc(17)*cos(twoPi*(pmc(25)*(xp-pmc(28))+pmc(26)
     & *(yp-pmc(29))+pmc(27)*(z0-pmc(30))-pmc(18)*(th))))
             end if

           u(i1,i2,i3,ex)=ue(ex)
           u(i1,i2,i3,ey)=ue(ey)
           u(i1,i2,i3,ez)=ue(ez)

           u(i1,i2,i3,hx)=ue(hx)
           u(i1,i2,i3,hy)=ue(hy)
           u(i1,i2,i3,hz)=ue(hz)

          end do
          end do
          end do

        end if ! nd


       else
         write(*,'("mxYeeIcErr:ERROR: unknown initialConditionOption = 
     & ",i6)') initialConditionOption
          ! ' 
         stop 9026

       end if ! TZ

      end if

      if( option.eq.1 )then
        ! ************************************************************
        ! ***************** Errors ***********************************
        ! ************************************************************

       do i=0,numberOfComponents-1
         err(i)=0.
         uNorm(i)=0.
       end do

       if( forcingOption.eq.twilightZoneForcing )then
          if( nd.eq.2 )then
           z0=0.
           do i3=n3a,n3b
           do i2=n2a,n2b
           do i1=n1a,n1b
            x0 = (xa(0)+dx(0)*(i1-nra(0)))
            y0 = (xa(1)+dx(1)*(i2-nra(1)))
            xp = x0 + .5*dx(0)
            yp = y0 + .5*dx(1)
              call ogDeriv(ep,0,0,0,0,xp,y0,z0,te,ex,ue(ex))
              call ogDeriv(ep,0,0,0,0,x0,yp,z0,te,ey,ue(ey))
              call ogDeriv(ep,0,0,0,0,xp,yp,z0,th,hz,ue(hz))
            v(i1,i2,i3,ex)=u(i1,i2,i3,ex) - ue(ex)
            v(i1,i2,i3,ey)=u(i1,i2,i3,ey) - ue(ey)
            v(i1,i2,i3,hz)=u(i1,i2,i3,hz) - ue(hz)
            ! Ex is not valid at i1=n1b  (mask=0 is a valid cell)
            if( i1.lt.n1b .and. (maskBodies.eq.0 .or. mask(i1,i2,i3)
     & .eq.0 .or. mask(i1,i2-1,i3).eq.0) )then
              err(ex)=max(err(ex),abs(v(i1,i2,i3,ex)))
              uNorm(ex)=max(uNorm(ex),abs(u(i1,i2,i3,ex)))
            else
              v(i1,i2,i3,ex)=0.  ! could set to v(i1-1,i2,i3,ex)
            end if
            ! Ey is not valid at i2=n2b
            if( i2.lt.n2b .and. (maskBodies.eq.0 .or. mask(i1,i2,i3)
     & .eq.0 .or. mask(i1-1,i2,i3).eq.0) )then
              err(ey)=max(err(ey),abs(v(i1,i2,i3,ey)))
              uNorm(ey)=max(uNorm(ey),abs(u(i1,i2,i3,ey)))
            else
              v(i1,i2,i3,ey)=0.
            end if
            if( i1.lt.n1b .and. i2.lt.n2b .and. (maskBodies.eq.0 .or. 
     & mask(i1,i2,i3).eq.0) )then
              err(hz)=max(err(hz),abs(v(i1,i2,i3,hz)))
              uNorm(hz)=max(uNorm(hz),abs(u(i1,i2,i3,hz)))
            else
              v(i1,i2,i3,hz)=0.
            end if
           end do
           end do
           end do
          else
            ! -- Errors in 3D ---
           do i3=n3a,n3b
           do i2=n2a,n2b
           do i1=n1a,n1b
            x0 = (xa(0)+dx(0)*(i1-nra(0)))
            y0 = (xa(1)+dx(1)*(i2-nra(1)))
            z0 = (xa(2)+dx(2)*(i3-nra(2)))
            xp = x0 + .5*dx(0)
            yp = y0 + .5*dx(1)
            zp = z0 + .5*dx(2)
             call ogDeriv(ep,0,0,0,0,xp,y0,z0,te,ex,ue(ex))
             call ogDeriv(ep,0,0,0,0,x0,yp,z0,te,ey,ue(ey))
             call ogDeriv(ep,0,0,0,0,x0,y0,zp,te,ez,ue(ez))
             call ogDeriv(ep,0,0,0,0,x0,yp,zp,th,hx,ue(hx))
             call ogDeriv(ep,0,0,0,0,xp,y0,zp,th,hy,ue(hy))
             call ogDeriv(ep,0,0,0,0,xp,yp,z0,th,hz,ue(hz))
            v(i1,i2,i3,ex) =u(i1,i2,i3,ex) - ue(ex)
            v(i1,i2,i3,ey) =u(i1,i2,i3,ey) - ue(ey)
            v(i1,i2,i3,ez) =u(i1,i2,i3,ez) - ue(ez)
            v(i1,i2,i3,hx) =u(i1,i2,i3,hx) - ue(hx)
            v(i1,i2,i3,hy) =u(i1,i2,i3,hy) - ue(hy)
            v(i1,i2,i3,hz) =u(i1,i2,i3,hz) - ue(hz)
            ! Ex is not valid at i1=n1b, the edge for Ex is valid if any of the 4 cells next to the edge are valid
            if( i1.lt.n1b .and. (maskBodies.eq.0 .or. mask(i1,i2,i3  )
     & .eq.0 .or. mask(i1,i2-1,i3  ).eq.0 .or. mask(i1,i2,i3-1).eq.0 
     & .or. mask(i1,i2-1,i3-1).eq.0 ) )then
              err(ex)=max(err(ex),abs(v(i1,i2,i3,ex)))
              uNorm(ex)=max(uNorm(ex),abs(u(i1,i2,i3,ex)))
            else
              v(i1,i2,i3,ex)=0.  ! could set to v(i1-1,i2,i3,ex)
            end if
            ! Ey is not valid at i2=n2b, the edge for Ey is valid if any of the 4 cells next to the edge are valid
            if( i2.lt.n2b .and. (maskBodies.eq.0 .or. mask(i1,i2,i3  )
     & .eq.0 .or. mask(i1-1,i2,i3  ).eq.0 .or. mask(i1,i2,i3-1).eq.0 
     & .or. mask(i1-1,i2,i3-1).eq.0 ) )then
              err(ey)=max(err(ey),abs(v(i1,i2,i3,ey)))
              uNorm(ey)=max(uNorm(ey),abs(u(i1,i2,i3,ey)))
            else
              v(i1,i2,i3,ey)=0.
            end if
            ! Ez is not valid at i3=n3b, the edge for Ez is valid if any of the 4 cells next to the edge are valid
            if( i3.lt.n3b .and. (maskBodies.eq.0 .or. mask(i1,i2  ,i3)
     & .eq.0 .or. mask(i1-1,i2  ,i3).eq.0 .or. mask(i1,i2-1,i3).eq.0 
     & .or. mask(i1-1,i2-1,i3).eq.0 ) )then
              err(ez)=max(err(ez),abs(v(i1,i2,i3,ez)))
              uNorm(ez)=max(uNorm(ez),abs(u(i1,i2,i3,ez)))
            else
              v(i1,i2,i3,ez)=0.
            end if
            !  The face Hz is valid if any of the 2 cells next to the face are valid
            if( i2.lt.n2b .and. i3.lt.n3b .and. (maskBodies.eq.0 .or. 
     & mask(i1,i2,i3).eq.0 .or. mask(i1-1,i2,i3).eq.0)  )then
              err(hx)=max(err(hx),abs(v(i1,i2,i3,hx)))
              uNorm(hx)=max(uNorm(hx),abs(u(i1,i2,i3,hx)))
            else
              v(i1,i2,i3,hx)=0.
            end if
            if( i1.lt.n1b .and. i3.lt.n3b .and. (maskBodies.eq.0 .or. 
     & mask(i1,i2,i3).eq.0 .or. mask(i1,i2-1,i3).eq.0)  )then
              err(hy)=max(err(hy),abs(v(i1,i2,i3,hy)))
              uNorm(hy)=max(uNorm(hy),abs(u(i1,i2,i3,hy)))
            else
              v(i1,i2,i3,hy)=0.
            end if
            if( i1.lt.n1b .and. i2.lt.n2b .and. (maskBodies.eq.0 .or. 
     & mask(i1,i2,i3).eq.0 .or. mask(i1,i2,i3-1).eq.0)  )then
              err(hz)=max(err(hz),abs(v(i1,i2,i3,hz)))
              uNorm(hz)=max(uNorm(hz),abs(u(i1,i2,i3,hz)))
            else
              v(i1,i2,i3,hz)=0.
            end if
           end do
           end do
           end do
          end if ! nd

       else if( initialConditionOption.eq.planeWaveInitialCondition )
     & then
          if( nd.eq.2 )then
           z0=0.
           do i3=n3a,n3b
           do i2=n2a,n2b
           do i1=n1a,n1b
            x0 = (xa(0)+dx(0)*(i1-nra(0)))
            y0 = (xa(1)+dx(1)*(i2-nra(1)))
            xp = x0 + .5*dx(0)
            yp = y0 + .5*dx(1)
              ue(ex) = sin(twoPi*(kx*(xp)+ky*(y0)-cc*(te)))*pwc(0)
              ue(ey) = sin(twoPi*(kx*(x0)+ky*(yp)-cc*(te)))*pwc(1)
              ue(hz) = sin(twoPi*(kx*(xp)+ky*(yp)-cc*(th)))*pwc(5)
              ! write(debugFile,'(" planeWave: Ex,exact=",2e10.2)') u(i1,i2,i3,ex),ue(ex)
            v(i1,i2,i3,ex)=u(i1,i2,i3,ex) - ue(ex)
            v(i1,i2,i3,ey)=u(i1,i2,i3,ey) - ue(ey)
            v(i1,i2,i3,hz)=u(i1,i2,i3,hz) - ue(hz)
            ! Ex is not valid at i1=n1b  (mask=0 is a valid cell)
            if( i1.lt.n1b .and. (maskBodies.eq.0 .or. mask(i1,i2,i3)
     & .eq.0 .or. mask(i1,i2-1,i3).eq.0) )then
              err(ex)=max(err(ex),abs(v(i1,i2,i3,ex)))
              uNorm(ex)=max(uNorm(ex),abs(u(i1,i2,i3,ex)))
            else
              v(i1,i2,i3,ex)=0.  ! could set to v(i1-1,i2,i3,ex)
            end if
            ! Ey is not valid at i2=n2b
            if( i2.lt.n2b .and. (maskBodies.eq.0 .or. mask(i1,i2,i3)
     & .eq.0 .or. mask(i1-1,i2,i3).eq.0) )then
              err(ey)=max(err(ey),abs(v(i1,i2,i3,ey)))
              uNorm(ey)=max(uNorm(ey),abs(u(i1,i2,i3,ey)))
            else
              v(i1,i2,i3,ey)=0.
            end if
            if( i1.lt.n1b .and. i2.lt.n2b .and. (maskBodies.eq.0 .or. 
     & mask(i1,i2,i3).eq.0) )then
              err(hz)=max(err(hz),abs(v(i1,i2,i3,hz)))
              uNorm(hz)=max(uNorm(hz),abs(u(i1,i2,i3,hz)))
            else
              v(i1,i2,i3,hz)=0.
            end if
           end do
           end do
           end do
          else
            ! -- Errors in 3D ---
           do i3=n3a,n3b
           do i2=n2a,n2b
           do i1=n1a,n1b
            x0 = (xa(0)+dx(0)*(i1-nra(0)))
            y0 = (xa(1)+dx(1)*(i2-nra(1)))
            z0 = (xa(2)+dx(2)*(i3-nra(2)))
            xp = x0 + .5*dx(0)
            yp = y0 + .5*dx(1)
            zp = z0 + .5*dx(2)
             ue(ex) = sin(twoPi*(kx*(xp)+ky*(y0)+kz*(z0)-cc*(te)))*pwc(
     & 0)
             ue(ey) = sin(twoPi*(kx*(x0)+ky*(yp)+kz*(z0)-cc*(te)))*pwc(
     & 1)
             ue(ez) = sin(twoPi*(kx*(x0)+ky*(y0)+kz*(zp)-cc*(te)))*pwc(
     & 2)
             ue(hx) = sin(twoPi*(kx*(x0)+ky*(yp)+kz*(zp)-cc*(th)))*pwc(
     & 3)
             ue(hy) = sin(twoPi*(kx*(xp)+ky*(y0)+kz*(zp)-cc*(th)))*pwc(
     & 4)
             ue(hz) = sin(twoPi*(kx*(xp)+ky*(yp)+kz*(z0)-cc*(th)))*pwc(
     & 5)
            v(i1,i2,i3,ex) =u(i1,i2,i3,ex) - ue(ex)
            v(i1,i2,i3,ey) =u(i1,i2,i3,ey) - ue(ey)
            v(i1,i2,i3,ez) =u(i1,i2,i3,ez) - ue(ez)
            v(i1,i2,i3,hx) =u(i1,i2,i3,hx) - ue(hx)
            v(i1,i2,i3,hy) =u(i1,i2,i3,hy) - ue(hy)
            v(i1,i2,i3,hz) =u(i1,i2,i3,hz) - ue(hz)
            ! Ex is not valid at i1=n1b, the edge for Ex is valid if any of the 4 cells next to the edge are valid
            if( i1.lt.n1b .and. (maskBodies.eq.0 .or. mask(i1,i2,i3  )
     & .eq.0 .or. mask(i1,i2-1,i3  ).eq.0 .or. mask(i1,i2,i3-1).eq.0 
     & .or. mask(i1,i2-1,i3-1).eq.0 ) )then
              err(ex)=max(err(ex),abs(v(i1,i2,i3,ex)))
              uNorm(ex)=max(uNorm(ex),abs(u(i1,i2,i3,ex)))
            else
              v(i1,i2,i3,ex)=0.  ! could set to v(i1-1,i2,i3,ex)
            end if
            ! Ey is not valid at i2=n2b, the edge for Ey is valid if any of the 4 cells next to the edge are valid
            if( i2.lt.n2b .and. (maskBodies.eq.0 .or. mask(i1,i2,i3  )
     & .eq.0 .or. mask(i1-1,i2,i3  ).eq.0 .or. mask(i1,i2,i3-1).eq.0 
     & .or. mask(i1-1,i2,i3-1).eq.0 ) )then
              err(ey)=max(err(ey),abs(v(i1,i2,i3,ey)))
              uNorm(ey)=max(uNorm(ey),abs(u(i1,i2,i3,ey)))
            else
              v(i1,i2,i3,ey)=0.
            end if
            ! Ez is not valid at i3=n3b, the edge for Ez is valid if any of the 4 cells next to the edge are valid
            if( i3.lt.n3b .and. (maskBodies.eq.0 .or. mask(i1,i2  ,i3)
     & .eq.0 .or. mask(i1-1,i2  ,i3).eq.0 .or. mask(i1,i2-1,i3).eq.0 
     & .or. mask(i1-1,i2-1,i3).eq.0 ) )then
              err(ez)=max(err(ez),abs(v(i1,i2,i3,ez)))
              uNorm(ez)=max(uNorm(ez),abs(u(i1,i2,i3,ez)))
            else
              v(i1,i2,i3,ez)=0.
            end if
            !  The face Hz is valid if any of the 2 cells next to the face are valid
            if( i2.lt.n2b .and. i3.lt.n3b .and. (maskBodies.eq.0 .or. 
     & mask(i1,i2,i3).eq.0 .or. mask(i1-1,i2,i3).eq.0)  )then
              err(hx)=max(err(hx),abs(v(i1,i2,i3,hx)))
              uNorm(hx)=max(uNorm(hx),abs(u(i1,i2,i3,hx)))
            else
              v(i1,i2,i3,hx)=0.
            end if
            if( i1.lt.n1b .and. i3.lt.n3b .and. (maskBodies.eq.0 .or. 
     & mask(i1,i2,i3).eq.0 .or. mask(i1,i2-1,i3).eq.0)  )then
              err(hy)=max(err(hy),abs(v(i1,i2,i3,hy)))
              uNorm(hy)=max(uNorm(hy),abs(u(i1,i2,i3,hy)))
            else
              v(i1,i2,i3,hy)=0.
            end if
            if( i1.lt.n1b .and. i2.lt.n2b .and. (maskBodies.eq.0 .or. 
     & mask(i1,i2,i3).eq.0 .or. mask(i1,i2,i3-1).eq.0)  )then
              err(hz)=max(err(hz),abs(v(i1,i2,i3,hz)))
              uNorm(hz)=max(uNorm(hz),abs(u(i1,i2,i3,hz)))
            else
              v(i1,i2,i3,hz)=0.
            end if
           end do
           end do
           end do
          end if ! nd

       else if( knownSolutionOption .eq. 
     & scatteringFromADiskKnownSolution .or.knownSolutionOption .eq. 
     & scatteringFromADielectricDiskKnownSolution 
     & .or.knownSolutionOption .eq. 
     & scatteringFromASphereKnownSolution .or.knownSolutionOption 
     & .eq. scatteringFromADielectricSphereKnownSolution )then

          if( nd.eq.2 )then
           z0=0.
           do i3=n3a,n3b
           do i2=n2a,n2b
           do i1=n1a,n1b
            x0 = (xa(0)+dx(0)*(i1-nra(0)))
            y0 = (xa(1)+dx(1)*(i2-nra(1)))
            xp = x0 + .5*dx(0)
            yp = y0 + .5*dx(1)
              ue(ex) = uKnown(i1,i2,i3,ex)*sinte+uKnown(i1,i2,i3,ex+3)*
     & coste
              ue(ey) = uKnown(i1,i2,i3,ey)*sinte+uKnown(i1,i2,i3,ey+3)*
     & coste
              ue(hz) = uKnown(i1,i2,i3,hz)*sinth+uKnown(i1,i2,i3,hz+3)*
     & costh
            v(i1,i2,i3,ex)=u(i1,i2,i3,ex) - ue(ex)
            v(i1,i2,i3,ey)=u(i1,i2,i3,ey) - ue(ey)
            v(i1,i2,i3,hz)=u(i1,i2,i3,hz) - ue(hz)
            ! Ex is not valid at i1=n1b  (mask=0 is a valid cell)
            if( i1.lt.n1b .and. (maskBodies.eq.0 .or. mask(i1,i2,i3)
     & .eq.0 .or. mask(i1,i2-1,i3).eq.0) )then
              err(ex)=max(err(ex),abs(v(i1,i2,i3,ex)))
              uNorm(ex)=max(uNorm(ex),abs(u(i1,i2,i3,ex)))
            else
              v(i1,i2,i3,ex)=0.  ! could set to v(i1-1,i2,i3,ex)
            end if
            ! Ey is not valid at i2=n2b
            if( i2.lt.n2b .and. (maskBodies.eq.0 .or. mask(i1,i2,i3)
     & .eq.0 .or. mask(i1-1,i2,i3).eq.0) )then
              err(ey)=max(err(ey),abs(v(i1,i2,i3,ey)))
              uNorm(ey)=max(uNorm(ey),abs(u(i1,i2,i3,ey)))
            else
              v(i1,i2,i3,ey)=0.
            end if
            if( i1.lt.n1b .and. i2.lt.n2b .and. (maskBodies.eq.0 .or. 
     & mask(i1,i2,i3).eq.0) )then
              err(hz)=max(err(hz),abs(v(i1,i2,i3,hz)))
              uNorm(hz)=max(uNorm(hz),abs(u(i1,i2,i3,hz)))
            else
              v(i1,i2,i3,hz)=0.
            end if
           end do
           end do
           end do
          else
            ! -- Errors in 3D ---
           do i3=n3a,n3b
           do i2=n2a,n2b
           do i1=n1a,n1b
            x0 = (xa(0)+dx(0)*(i1-nra(0)))
            y0 = (xa(1)+dx(1)*(i2-nra(1)))
            z0 = (xa(2)+dx(2)*(i3-nra(2)))
            xp = x0 + .5*dx(0)
            yp = y0 + .5*dx(1)
            zp = z0 + .5*dx(2)
             ue(ex) = uKnown(i1,i2,i3,ex)*sinte+uKnown(i1,i2,i3,ex+6)*
     & coste
             ue(ey) = uKnown(i1,i2,i3,ey)*sinte+uKnown(i1,i2,i3,ey+6)*
     & coste
             ue(ez) = uKnown(i1,i2,i3,ez)*sinte+uKnown(i1,i2,i3,ez+6)*
     & coste
             ue(hx) = uKnown(i1,i2,i3,hx)*sinth+uKnown(i1,i2,i3,hx+6)*
     & costh
             ue(hy) = uKnown(i1,i2,i3,hy)*sinth+uKnown(i1,i2,i3,hy+6)*
     & costh
             ue(hz) = uKnown(i1,i2,i3,hz)*sinth+uKnown(i1,i2,i3,hz+6)*
     & costh
            v(i1,i2,i3,ex) =u(i1,i2,i3,ex) - ue(ex)
            v(i1,i2,i3,ey) =u(i1,i2,i3,ey) - ue(ey)
            v(i1,i2,i3,ez) =u(i1,i2,i3,ez) - ue(ez)
            v(i1,i2,i3,hx) =u(i1,i2,i3,hx) - ue(hx)
            v(i1,i2,i3,hy) =u(i1,i2,i3,hy) - ue(hy)
            v(i1,i2,i3,hz) =u(i1,i2,i3,hz) - ue(hz)
            ! Ex is not valid at i1=n1b, the edge for Ex is valid if any of the 4 cells next to the edge are valid
            if( i1.lt.n1b .and. (maskBodies.eq.0 .or. mask(i1,i2,i3  )
     & .eq.0 .or. mask(i1,i2-1,i3  ).eq.0 .or. mask(i1,i2,i3-1).eq.0 
     & .or. mask(i1,i2-1,i3-1).eq.0 ) )then
              err(ex)=max(err(ex),abs(v(i1,i2,i3,ex)))
              uNorm(ex)=max(uNorm(ex),abs(u(i1,i2,i3,ex)))
            else
              v(i1,i2,i3,ex)=0.  ! could set to v(i1-1,i2,i3,ex)
            end if
            ! Ey is not valid at i2=n2b, the edge for Ey is valid if any of the 4 cells next to the edge are valid
            if( i2.lt.n2b .and. (maskBodies.eq.0 .or. mask(i1,i2,i3  )
     & .eq.0 .or. mask(i1-1,i2,i3  ).eq.0 .or. mask(i1,i2,i3-1).eq.0 
     & .or. mask(i1-1,i2,i3-1).eq.0 ) )then
              err(ey)=max(err(ey),abs(v(i1,i2,i3,ey)))
              uNorm(ey)=max(uNorm(ey),abs(u(i1,i2,i3,ey)))
            else
              v(i1,i2,i3,ey)=0.
            end if
            ! Ez is not valid at i3=n3b, the edge for Ez is valid if any of the 4 cells next to the edge are valid
            if( i3.lt.n3b .and. (maskBodies.eq.0 .or. mask(i1,i2  ,i3)
     & .eq.0 .or. mask(i1-1,i2  ,i3).eq.0 .or. mask(i1,i2-1,i3).eq.0 
     & .or. mask(i1-1,i2-1,i3).eq.0 ) )then
              err(ez)=max(err(ez),abs(v(i1,i2,i3,ez)))
              uNorm(ez)=max(uNorm(ez),abs(u(i1,i2,i3,ez)))
            else
              v(i1,i2,i3,ez)=0.
            end if
            !  The face Hz is valid if any of the 2 cells next to the face are valid
            if( i2.lt.n2b .and. i3.lt.n3b .and. (maskBodies.eq.0 .or. 
     & mask(i1,i2,i3).eq.0 .or. mask(i1-1,i2,i3).eq.0)  )then
              err(hx)=max(err(hx),abs(v(i1,i2,i3,hx)))
              uNorm(hx)=max(uNorm(hx),abs(u(i1,i2,i3,hx)))
            else
              v(i1,i2,i3,hx)=0.
            end if
            if( i1.lt.n1b .and. i3.lt.n3b .and. (maskBodies.eq.0 .or. 
     & mask(i1,i2,i3).eq.0 .or. mask(i1,i2-1,i3).eq.0)  )then
              err(hy)=max(err(hy),abs(v(i1,i2,i3,hy)))
              uNorm(hy)=max(uNorm(hy),abs(u(i1,i2,i3,hy)))
            else
              v(i1,i2,i3,hy)=0.
            end if
            if( i1.lt.n1b .and. i2.lt.n2b .and. (maskBodies.eq.0 .or. 
     & mask(i1,i2,i3).eq.0 .or. mask(i1,i2,i3-1).eq.0)  )then
              err(hz)=max(err(hz),abs(v(i1,i2,i3,hz)))
              uNorm(hz)=max(uNorm(hz),abs(u(i1,i2,i3,hz)))
            else
              v(i1,i2,i3,hz)=0.
            end if
           end do
           end do
           end do
          end if ! nd

       else if( initialConditionOption .eq. 
     & planeMaterialInterfaceInitialCondition )then

          if( nd.eq.2 )then
           z0=0.
           do i3=n3a,n3b
           do i2=n2a,n2b
           do i1=n1a,n1b
            x0 = (xa(0)+dx(0)*(i1-nra(0)))
            y0 = (xa(1)+dx(1)*(i2-nra(1)))
            xp = x0 + .5*dx(0)
            yp = y0 + .5*dx(1)
               ! **  We check the location of each component ***
                if( nPMI(0)*(xp-xPMI(0)) + nPMI(1)*(y0-xPMI(1)) .le. 
     & 0. )then
                 ! incident + reflected wave: 
                 ue(ex) = (pmc(0)*cos(twoPi*(pmc(19)*(xp-pmc(28))+pmc(
     & 20)*(y0-pmc(29))+pmc(21)*(z0-pmc(30))-pmc(18)*(te)))+pmc(1)*
     & cos(twoPi*(pmc(22)*(xp-pmc(28))+pmc(23)*(y0-pmc(29))+pmc(24)*(
     & z0-pmc(30))-pmc(18)*(te))))
                else
                 ! transmitted wave: 
                 ue(ex) = (pmc(12)*cos(twoPi*(pmc(25)*(xp-pmc(28))+pmc(
     & 26)*(y0-pmc(29))+pmc(27)*(z0-pmc(30))-pmc(18)*(te))))
                end if
                if( nPMI(0)*(x0-xPMI(0)) + nPMI(1)*(yp-xPMI(1)) .le. 
     & 0. )then
                 ! incident + reflected wave: 
                 ue(ey) = (pmc(2)*cos(twoPi*(pmc(19)*(x0-pmc(28))+pmc(
     & 20)*(yp-pmc(29))+pmc(21)*(z0-pmc(30))-pmc(18)*(te)))+pmc(3)*
     & cos(twoPi*(pmc(22)*(x0-pmc(28))+pmc(23)*(yp-pmc(29))+pmc(24)*(
     & z0-pmc(30))-pmc(18)*(te))))
                else
                 ! transmitted wave: 
                 ue(ey) = (pmc(13)*cos(twoPi*(pmc(25)*(x0-pmc(28))+pmc(
     & 26)*(yp-pmc(29))+pmc(27)*(z0-pmc(30))-pmc(18)*(te))))
                end if
                if( nPMI(0)*(xp-xPMI(0)) + nPMI(1)*(yp-xPMI(1)) .le. 
     & 0. )then
                 ! incident + reflected wave: 
                 ue(hz) = (pmc(10)*cos(twoPi*(pmc(19)*(xp-pmc(28))+pmc(
     & 20)*(yp-pmc(29))+pmc(21)*(z0-pmc(30))-pmc(18)*(th)))+pmc(11)*
     & cos(twoPi*(pmc(22)*(xp-pmc(28))+pmc(23)*(yp-pmc(29))+pmc(24)*(
     & z0-pmc(30))-pmc(18)*(th))))
                else
                 ! transmitted wave: 
                 ue(hz) = (pmc(17)*cos(twoPi*(pmc(25)*(xp-pmc(28))+pmc(
     & 26)*(yp-pmc(29))+pmc(27)*(z0-pmc(30))-pmc(18)*(th))))
                end if
            v(i1,i2,i3,ex)=u(i1,i2,i3,ex) - ue(ex)
            v(i1,i2,i3,ey)=u(i1,i2,i3,ey) - ue(ey)
            v(i1,i2,i3,hz)=u(i1,i2,i3,hz) - ue(hz)
            ! Ex is not valid at i1=n1b  (mask=0 is a valid cell)
            if( i1.lt.n1b .and. (maskBodies.eq.0 .or. mask(i1,i2,i3)
     & .eq.0 .or. mask(i1,i2-1,i3).eq.0) )then
              err(ex)=max(err(ex),abs(v(i1,i2,i3,ex)))
              uNorm(ex)=max(uNorm(ex),abs(u(i1,i2,i3,ex)))
            else
              v(i1,i2,i3,ex)=0.  ! could set to v(i1-1,i2,i3,ex)
            end if
            ! Ey is not valid at i2=n2b
            if( i2.lt.n2b .and. (maskBodies.eq.0 .or. mask(i1,i2,i3)
     & .eq.0 .or. mask(i1-1,i2,i3).eq.0) )then
              err(ey)=max(err(ey),abs(v(i1,i2,i3,ey)))
              uNorm(ey)=max(uNorm(ey),abs(u(i1,i2,i3,ey)))
            else
              v(i1,i2,i3,ey)=0.
            end if
            if( i1.lt.n1b .and. i2.lt.n2b .and. (maskBodies.eq.0 .or. 
     & mask(i1,i2,i3).eq.0) )then
              err(hz)=max(err(hz),abs(v(i1,i2,i3,hz)))
              uNorm(hz)=max(uNorm(hz),abs(u(i1,i2,i3,hz)))
            else
              v(i1,i2,i3,hz)=0.
            end if
           end do
           end do
           end do
          else
            ! -- Errors in 3D ---
           do i3=n3a,n3b
           do i2=n2a,n2b
           do i1=n1a,n1b
            x0 = (xa(0)+dx(0)*(i1-nra(0)))
            y0 = (xa(1)+dx(1)*(i2-nra(1)))
            z0 = (xa(2)+dx(2)*(i3-nra(2)))
            xp = x0 + .5*dx(0)
            yp = y0 + .5*dx(1)
            zp = z0 + .5*dx(2)
               ! **  We check the location of each component ***
                if( nPMI(0)*(xp-xPMI(0)) + nPMI(1)*(y0-xPMI(1)) + nPMI(
     & 2)*(z0-xPMI(2)) .le. 0. )then
                 ue(ex) = (pmc(0)*cos(twoPi*(pmc(19)*(xp-pmc(28))+pmc(
     & 20)*(y0-pmc(29))+pmc(21)*(z0-pmc(30))-pmc(18)*(te)))+pmc(1)*
     & cos(twoPi*(pmc(22)*(xp-pmc(28))+pmc(23)*(y0-pmc(29))+pmc(24)*(
     & z0-pmc(30))-pmc(18)*(te))))
                else
                 ue(ex) = (pmc(12)*cos(twoPi*(pmc(25)*(xp-pmc(28))+pmc(
     & 26)*(y0-pmc(29))+pmc(27)*(z0-pmc(30))-pmc(18)*(te))))
                end if
                if( nPMI(0)*(x0-xPMI(0)) + nPMI(1)*(yp-xPMI(1)) + nPMI(
     & 2)*(z0-xPMI(2)) .le. 0. )then
                 ue(ey) = (pmc(2)*cos(twoPi*(pmc(19)*(x0-pmc(28))+pmc(
     & 20)*(yp-pmc(29))+pmc(21)*(z0-pmc(30))-pmc(18)*(te)))+pmc(3)*
     & cos(twoPi*(pmc(22)*(x0-pmc(28))+pmc(23)*(yp-pmc(29))+pmc(24)*(
     & z0-pmc(30))-pmc(18)*(te))))
                else
                 ue(ey) = (pmc(13)*cos(twoPi*(pmc(25)*(x0-pmc(28))+pmc(
     & 26)*(yp-pmc(29))+pmc(27)*(z0-pmc(30))-pmc(18)*(te))))
                end if
                if( nPMI(0)*(x0-xPMI(0)) + nPMI(1)*(y0-xPMI(1)) + nPMI(
     & 2)*(zp-xPMI(2)) .le. 0. )then
                 ue(ez) = (pmc(4)*cos(twoPi*(pmc(19)*(x0-pmc(28))+pmc(
     & 20)*(y0-pmc(29))+pmc(21)*(zp-pmc(30))-pmc(18)*(te)))+pmc(5)*
     & cos(twoPi*(pmc(22)*(x0-pmc(28))+pmc(23)*(y0-pmc(29))+pmc(24)*(
     & zp-pmc(30))-pmc(18)*(te))))
                else
                 ue(ez) = (pmc(14)*cos(twoPi*(pmc(25)*(x0-pmc(28))+pmc(
     & 26)*(y0-pmc(29))+pmc(27)*(zp-pmc(30))-pmc(18)*(te))))
                end if
                if( nPMI(0)*(x0-xPMI(0)) + nPMI(1)*(yp-xPMI(1)) + nPMI(
     & 2)*(zp-xPMI(2)) .le. 0. )then
                 ue(hx) = (pmc(6)*cos(twoPi*(pmc(19)*(x0-pmc(28))+pmc(
     & 20)*(yp-pmc(29))+pmc(21)*(zp-pmc(30))-pmc(18)*(th)))+pmc(7)*
     & cos(twoPi*(pmc(22)*(x0-pmc(28))+pmc(23)*(yp-pmc(29))+pmc(24)*(
     & zp-pmc(30))-pmc(18)*(th))))
                else
                 ue(hx) = (pmc(15)*cos(twoPi*(pmc(25)*(x0-pmc(28))+pmc(
     & 26)*(yp-pmc(29))+pmc(27)*(zp-pmc(30))-pmc(18)*(th))))
                end if
                if( nPMI(0)*(xp-xPMI(0)) + nPMI(1)*(y0-xPMI(1)) + nPMI(
     & 2)*(zp-xPMI(2)) .le. 0. )then
                 ue(hy) = (pmc(8)*cos(twoPi*(pmc(19)*(xp-pmc(28))+pmc(
     & 20)*(y0-pmc(29))+pmc(21)*(zp-pmc(30))-pmc(18)*(th)))+pmc(9)*
     & cos(twoPi*(pmc(22)*(xp-pmc(28))+pmc(23)*(y0-pmc(29))+pmc(24)*(
     & zp-pmc(30))-pmc(18)*(th))))
                else
                 ue(hy) = (pmc(16)*cos(twoPi*(pmc(25)*(xp-pmc(28))+pmc(
     & 26)*(y0-pmc(29))+pmc(27)*(zp-pmc(30))-pmc(18)*(th))))
                end if
                if( nPMI(0)*(xp-xPMI(0)) + nPMI(1)*(yp-xPMI(1)) + nPMI(
     & 2)*(z0-xPMI(2)) .le. 0. )then
                 ue(hz) = (pmc(10)*cos(twoPi*(pmc(19)*(xp-pmc(28))+pmc(
     & 20)*(yp-pmc(29))+pmc(21)*(z0-pmc(30))-pmc(18)*(th)))+pmc(11)*
     & cos(twoPi*(pmc(22)*(xp-pmc(28))+pmc(23)*(yp-pmc(29))+pmc(24)*(
     & z0-pmc(30))-pmc(18)*(th))))
                else
                 ue(hz) = (pmc(17)*cos(twoPi*(pmc(25)*(xp-pmc(28))+pmc(
     & 26)*(yp-pmc(29))+pmc(27)*(z0-pmc(30))-pmc(18)*(th))))
                end if
            v(i1,i2,i3,ex) =u(i1,i2,i3,ex) - ue(ex)
            v(i1,i2,i3,ey) =u(i1,i2,i3,ey) - ue(ey)
            v(i1,i2,i3,ez) =u(i1,i2,i3,ez) - ue(ez)
            v(i1,i2,i3,hx) =u(i1,i2,i3,hx) - ue(hx)
            v(i1,i2,i3,hy) =u(i1,i2,i3,hy) - ue(hy)
            v(i1,i2,i3,hz) =u(i1,i2,i3,hz) - ue(hz)
            ! Ex is not valid at i1=n1b, the edge for Ex is valid if any of the 4 cells next to the edge are valid
            if( i1.lt.n1b .and. (maskBodies.eq.0 .or. mask(i1,i2,i3  )
     & .eq.0 .or. mask(i1,i2-1,i3  ).eq.0 .or. mask(i1,i2,i3-1).eq.0 
     & .or. mask(i1,i2-1,i3-1).eq.0 ) )then
              err(ex)=max(err(ex),abs(v(i1,i2,i3,ex)))
              uNorm(ex)=max(uNorm(ex),abs(u(i1,i2,i3,ex)))
            else
              v(i1,i2,i3,ex)=0.  ! could set to v(i1-1,i2,i3,ex)
            end if
            ! Ey is not valid at i2=n2b, the edge for Ey is valid if any of the 4 cells next to the edge are valid
            if( i2.lt.n2b .and. (maskBodies.eq.0 .or. mask(i1,i2,i3  )
     & .eq.0 .or. mask(i1-1,i2,i3  ).eq.0 .or. mask(i1,i2,i3-1).eq.0 
     & .or. mask(i1-1,i2,i3-1).eq.0 ) )then
              err(ey)=max(err(ey),abs(v(i1,i2,i3,ey)))
              uNorm(ey)=max(uNorm(ey),abs(u(i1,i2,i3,ey)))
            else
              v(i1,i2,i3,ey)=0.
            end if
            ! Ez is not valid at i3=n3b, the edge for Ez is valid if any of the 4 cells next to the edge are valid
            if( i3.lt.n3b .and. (maskBodies.eq.0 .or. mask(i1,i2  ,i3)
     & .eq.0 .or. mask(i1-1,i2  ,i3).eq.0 .or. mask(i1,i2-1,i3).eq.0 
     & .or. mask(i1-1,i2-1,i3).eq.0 ) )then
              err(ez)=max(err(ez),abs(v(i1,i2,i3,ez)))
              uNorm(ez)=max(uNorm(ez),abs(u(i1,i2,i3,ez)))
            else
              v(i1,i2,i3,ez)=0.
            end if
            !  The face Hz is valid if any of the 2 cells next to the face are valid
            if( i2.lt.n2b .and. i3.lt.n3b .and. (maskBodies.eq.0 .or. 
     & mask(i1,i2,i3).eq.0 .or. mask(i1-1,i2,i3).eq.0)  )then
              err(hx)=max(err(hx),abs(v(i1,i2,i3,hx)))
              uNorm(hx)=max(uNorm(hx),abs(u(i1,i2,i3,hx)))
            else
              v(i1,i2,i3,hx)=0.
            end if
            if( i1.lt.n1b .and. i3.lt.n3b .and. (maskBodies.eq.0 .or. 
     & mask(i1,i2,i3).eq.0 .or. mask(i1,i2-1,i3).eq.0)  )then
              err(hy)=max(err(hy),abs(v(i1,i2,i3,hy)))
              uNorm(hy)=max(uNorm(hy),abs(u(i1,i2,i3,hy)))
            else
              v(i1,i2,i3,hy)=0.
            end if
            if( i1.lt.n1b .and. i2.lt.n2b .and. (maskBodies.eq.0 .or. 
     & mask(i1,i2,i3).eq.0 .or. mask(i1,i2,i3-1).eq.0)  )then
              err(hz)=max(err(hz),abs(v(i1,i2,i3,hz)))
              uNorm(hz)=max(uNorm(hz),abs(u(i1,i2,i3,hz)))
            else
              v(i1,i2,i3,hz)=0.
            end if
           end do
           end do
           end do
          end if ! nd

       else
         write(*,'("mxYeeIcErr:ERROR: unknown forcing option = ",i6)') 
     & forcingOption
         write(*,'("mxYeeIcErr:INFO: knownSolutionOption = ",i6)') 
     & knownSolutionOption
         stop 9027
       end if

       ! write(*,'(" mxYeeIcErr: err  =",6e10.2)') (err(i),i=0,numberOfComponents-1)
       ! write(*,'(" mxYeeIcErr: uNorm=",6e10.2)') (uNorm(i),i=0,numberOfComponents-1)

       do i=0,numberOfComponents-1
         rpar(20+i)=err(i)
         rpar(20+i+numberOfComponents)=uNorm(i)
       end do
      end if


      if( option.eq.2 )then
       ! ************************************************************
       ! ***************** Divergence         ***********************
       ! ************************************************************

        ! *** note: this should be div( eps*E ) ********* fix me 

        ! write(debugFile,'(" mxYeeIcErr: compute div, ed,hd=",2i3)') ed,hd

        divEMax=0.
        divHMax=0.
        gradEMax=0.
        gradHMax=0.

        ! div(E) lives at the nodes
        if( ed.ge.0 )then
         if( nd.eq.2 )then
          n1a=gridIndexRange(0,0)+1
          n2a=gridIndexRange(0,1)+1
          n1b=gridIndexRange(1,0)-1
          n2b=gridIndexRange(1,1)-1
          do i3=n3a,n3b
          do i2=n2a,n2b
          do i1=n1a,n1b
           v(i1,i2,i3,ed) = (u(i1,i2,i3,ex)-u(i1-1,i2,i3,ex))/dx(0) + (
     & u(i1,i2,i3,ey)-u(i1,i2-1,i3,ey))/dx(1)
           divEMax=max(divEMax,abs(v(i1,i2,i3,ed)))
           gradEMax=max(gradEMax,abs((u(i1,i2,i3,ex)-u(i1-1,i2,i3,ex))
     & /dx(0)),abs((u(i1,i2,i3,ey)-u(i1,i2-1,i3,ey))/dx(1)))
          end do
          end do
          end do
         else
          n1a=gridIndexRange(0,0)+1
          n2a=gridIndexRange(0,1)+1
          n3a=gridIndexRange(0,2)+1
          n1b=gridIndexRange(1,0)-1
          n2b=gridIndexRange(1,1)-1
          n3b=gridIndexRange(1,2)-1
          do i3=n3a,n3b
          do i2=n2a,n2b
          do i1=n1a,n1b
           v(i1,i2,i3,ed) = (u(i1,i2,i3,ex)-u(i1-1,i2,i3,ex))/dx(0) + (
     & u(i1,i2,i3,ey)-u(i1,i2-1,i3,ey))/dx(1) + (u(i1,i2,i3,ez)-u(i1,
     & i2,i3-1,ez))/dx(2)
           divEMax=max(divEMax,abs(v(i1,i2,i3,ed)))
           gradEMax=max(gradEMax,abs((u(i1,i2,i3,ex)-u(i1-1,i2,i3,ex))
     & /dx(0)),abs((u(i1,i2,i3,ey)-u(i1,i2-1,i3,ey))/dx(1)),abs((u(i1,
     & i2,i3,ez)-u(i1,i2,i3-1,ez))/dx(2)))
          end do
          end do
          end do
         end if ! nd
        end if
        ! div(H) is cell centered
        if( hd.ge.0 .and. nd.eq.3 )then
          n1b=gridIndexRange(1,0)-1
          n2b=gridIndexRange(1,1)-1
          n3b=gridIndexRange(1,2)-1
          do i3=n3a,n3b
          do i2=n2a,n2b
          do i1=n1a,n1b
           v(i1,i2,i3,hd) = (u(i1+1,i2,i3,hx)-u(i1,i2,i3,hx))/dx(0) + (
     & u(i1,i2+1,i3,hy)-u(i1,i2,i3,hy))/dx(1) + (u(i1,i2,i3+1,hz)-u(
     & i1,i2,i3,hz))/dx(2)
           divHMax=max(divHMax,abs(v(i1,i2,i3,ed)))
           gradEMax=max(gradEMax,abs((u(i1+1,i2,i3,hx)-u(i1,i2,i3,hx))
     & /dx(0)),abs((u(i1,i2+1,i3,hy)-u(i1,i2,i3,hy))/dx(1)),abs((u(i1,
     & i2,i3+1,hz)-u(i1,i2,i3,hz))/dx(2)))
          end do
          end do
          end do
        end if ! nd
        ! reset 
        n1a=gridIndexRange(0,0)
        n1b=gridIndexRange(1,0)
        n2a=gridIndexRange(0,1)
        n2b=gridIndexRange(1,1)
        n3a=gridIndexRange(0,2)
        n3b=gridIndexRange(1,2)

        rpar(30)=divEMax
        rpar(31)=divHMax
        rpar(32)=gradEMax
        rpar(33)=gradHMax

      end if


      if( option.eq.3 )then
       ! determine a node centered version of the fields for plotting 


       !      x----- -- x
       !      |         |
       !      Ey  Hz    |
       !      |         |
       !      x-- Ex -- x
       !  (i1,i2)

c$$$       ! *******************************************************************************************************
c$$$       if( .false. )then ! *** we cannot do this -- changes values we now may use ---
c$$$       ! Step I : extrapolate boundary points
c$$$       extra=1 ! assign an extra pt in the tangential direction so that we get corners assigned (eventually)
c$$$       numberOfGhostPoints=0
c$$$       beginLoopOverSides(extra,numberOfGhostPoints)
c$$$        en = ex+axis
c$$$        et1 = ex + axisp1
c$$$        et2 = ex + axisp2
c$$$        if( nd.eq.2 )then
c$$$         ! --- 2d ---
c$$$         if( side.eq.0 )then
c$$$           beginLoops()
c$$$            u(i1-is1,i2-is2,i3,ex)=extrap2(u,i1,i2,i3,ex,is1,is2,is3)
c$$$            u(i1-is1,i2-is2,i3,ey)=extrap2(u,i1,i2,i3,ey,is1,is2,is3)
c$$$            u(i1-is1,i2-is2,i3,hz)=extrap2(u,i1,i2,i3,hz,is1,is2,is3)
c$$$           endLoops()
c$$$         else
c$$$           ! on the right side we extrap differently for Hz and the normal component of the field
c$$$           beginLoops()
c$$$            u(i1-is1,i2-is2,i3,et1)=extrap2(u,i1,i2,i3,et1,is1,is2,is3)
c$$$            u(i1,i2,i3,en)=extrap2(u,i1+is1,i2+is2,i3,en,is1,is2,is3)
c$$$            u(i1,i2,i3,hz)=extrap2(u,i1+is1,i2+is2,i3,hz,is1,is2,is3)
c$$$           endLoops()
c$$$         end if
c$$$        else
c$$$         ! --- 3d ---
c$$$         if( side.eq.0 )then
c$$$           beginLoops()
c$$$            u(i1-is1,i2-is2,i3-is3,ex)=extrap2(u,i1,i2,i3,ex,is1,is2,is3)
c$$$            u(i1-is1,i2-is2,i3-is3,ey)=extrap2(u,i1,i2,i3,ey,is1,is2,is3)
c$$$            u(i1-is1,i2-is2,i3-is3,ez)=extrap2(u,i1,i2,i3,ez,is1,is2,is3)
c$$$
c$$$            u(i1-is1,i2-is2,i3-is3,hx)=extrap2(u,i1,i2,i3,hx,is1,is2,is3)
c$$$            u(i1-is1,i2-is2,i3-is3,hy)=extrap2(u,i1,i2,i3,hy,is1,is2,is3)
c$$$            u(i1-is1,i2-is2,i3-is3,hz)=extrap2(u,i1,i2,i3,hz,is1,is2,is3)
c$$$           endLoops()
c$$$         else
c$$$           ! on the right side we extrap differently for Hz and the normal component of the field
c$$$           beginLoops()
c$$$            u(i1-is1,i2-is2,i3-is3,et1)=extrap2(u,i1,i2,i3,et1,is1,is2,is3)
c$$$            u(i1-is1,i2-is2,i3-is3,et2)=extrap2(u,i1,i2,i3,et2,is1,is2,is3)
c$$$            u(i1,i2,i3,en)=extrap2(u,i1+is1,i2+is2,i3+is3,en,is1,is2,is3)
c$$$
c$$$            u(i1,i2,i3,hx)=extrap2(u,i1+is1,i2+is2,i3+is3,hx,is1,is2,is3)
c$$$            u(i1,i2,i3,hy)=extrap2(u,i1+is1,i2+is2,i3+is3,hy,is1,is2,is3)
c$$$            u(i1,i2,i3,hz)=extrap2(u,i1+is1,i2+is2,i3+is3,hz,is1,is2,is3)
c$$$           endLoops()
c$$$         end if
c$$$        end if
c$$$       endLoopOverSides() 
c$$$       end if
       ! *******************************************************************************************************

       ! Step I : interpolate to nodes
       !   E : located on edges
       !   H : located on cell centres

       if( nd.eq.2 )then
        do i3=n3a,n3b
        do i2=n2a,n2b
        do i1=n1a,n1b
         v(i1,i2,i3,ex) = .5*(u(i1-1,i2,i3,ex) + u(i1,i2,i3,ex))
         v(i1,i2,i3,ey) = .5*(u(i1,i2-1,i3,ey) + u(i1,i2,i3,ey))
         v(i1,i2,i3,hz) = .25*(u(i1-1,i2-1,i3,hz) + u(i1,i2-1,i3,hz) + 
     & u(i1-1,i2,i3,hz) + u(i1,i2,i3,hz))
        end do
        end do
        end do

        ! Now fixup values on the boundaries that were not averaged correctly
        ! For the E-field only the normal component is wrong
        extra=0
        numberOfGhostPoints=0
         extra1a=extra
         extra1b=extra
         extra2a=extra
         extra2b=extra
         if( nd.eq.3 )then
           extra3a=extra
           extra3b=extra
         else
           extra3a=0
           extra3b=0
         end if
         if( boundaryCondition(0,0).lt.0 )then
           extra1a=max(0,extra1a) ! over-ride extra=-1 : assign ends in periodic directions
           extra1b=extra1a
         else
           if( boundaryCondition(0,0).eq.0 )then
             extra1a=numberOfGhostPoints  ! include interpolation points since we assign ghost points outside these
           end if
           if( boundaryCondition(1,0).eq.0 )then
             extra1b=numberOfGhostPoints
           end if
         end if
         if( boundaryCondition(0,1).lt.0 )then
          extra2a=max(0,extra2a) ! over-ride extra=-1 : assign ends in periodic directions
          extra2b=extra2a
         else
           if( boundaryCondition(0,1).eq.0 )then
             extra2a=numberOfGhostPoints
           end if
           if( boundaryCondition(1,1).eq.0 )then
             extra2b=numberOfGhostPoints
           end if
         end if
         if(  nd.eq.3 .and. boundaryCondition(0,2).lt.0 )then
          extra3a=max(0,extra3a) ! over-ride extra=-1 : assign ends in periodic directions
          extra3b=extra3a
         else
           if( boundaryCondition(0,2).eq.0 )then
             extra3a=numberOfGhostPoints
           end if
           if( boundaryCondition(1,2).eq.0 )then
             extra3b=numberOfGhostPoints
           end if
         end if
         do axis=0,nd-1
         do side=0,1
           n1a=gridIndexRange(0,0)-extra1a
           n1b=gridIndexRange(1,0)+extra1b
           n2a=gridIndexRange(0,1)-extra2a
           n2b=gridIndexRange(1,1)+extra2b
           n3a=gridIndexRange(0,2)-extra3a
           n3b=gridIndexRange(1,2)+extra3b
           if( axis.eq.0 )then
             n1a=gridIndexRange(side,axis)
             n1b=gridIndexRange(side,axis)
           else if( axis.eq.1 )then
             n2a=gridIndexRange(side,axis)
             n2b=gridIndexRange(side,axis)
           else
             n3a=gridIndexRange(side,axis)
             n3b=gridIndexRange(side,axis)
           end if
           is1=0
           is2=0
           is3=0
           if( axis.eq.0 )then
             is1=1-2*side
           else if( axis.eq.1 )then
             is2=1-2*side
           else if( axis.eq.2 )then
             is3=1-2*side
           else
             stop 5
           end if
           axisp1=mod(axis+1,nd)
           axisp2=mod(axis+2,nd)
        !    ! (js1,js2,js3) used to compute tangential derivatives
        !    js1=0
        !    js2=0
        !    js3=0
        !    if( axisp1.eq.0 )then
        !      js1=1-2*side
        !    else if( axisp1.eq.1 )then
        !      js2=1-2*side
        !    else if( axisp1.eq.2 )then
        !      js3=1-2*side
        !    else
        !      stop 5
        !    end if
        ! 
        !    ! (ks1,ks2,ks3) used to compute second tangential derivative
        !    ks1=0
        !    ks2=0
        !    ks3=0
        !    if( axisp2.eq.0 )then
        !      ks1=1-2*side
        !    else if( axisp2.eq.1 )then
        !      ks2=1-2*side
        !    else if( axisp2.eq.2 )then
        !      ks3=1-2*side
        !    else
        !      stop 5
        !    end if
         if( boundaryCondition(side,axis).gt.0 )then
          en = ex+axis
          do i3=n3a,n3b
          do i2=n2a,n2b
          do i1=n1a,n1b
            v(i1,i2,i3,en)=(2.*v(i1+is1,i2+is2,i3,en)-v(i1+is1+is1,i2+
     & is2+is2,i3+is3,en))
            v(i1,i2,i3,hz)=(2.*v(i1+is1,i2+is2,i3,hz)-v(i1+is1+is1,i2+
     & is2+is2,i3+is3,hz))
          end do
          end do
          end do
         end if
         end do
         end do
         ! reset these values
         n1a=gridIndexRange(0,0)
         n1b=gridIndexRange(1,0)
         n2a=gridIndexRange(0,1)
         n2b=gridIndexRange(1,1)
         n3a=gridIndexRange(0,2)
         n3b=gridIndexRange(1,2)

       else
         ! -- 3D  ---
        do i3=n3a,n3b
        do i2=n2a,n2b
        do i1=n1a,n1b
         v(i1,i2,i3,ex) = .5*(u(i1-1,i2,i3,ex) + u(i1,i2,i3,ex))
         v(i1,i2,i3,ey) = .5*(u(i1,i2-1,i3,ey) + u(i1,i2,i3,ey))
         v(i1,i2,i3,ez) = .5*(u(i1,i2,i3-1,ez) + u(i1,i2,i3,ez))

         v(i1,i2,i3,hx) =.125*(u(i1-1,i2-1,i3-1,hx) + u(i1,i2-1,i3-1,
     & hx) + u(i1-1,i2,i3-1,hx) + u(i1,i2,i3-1,hx) +u(i1-1,i2-1,i3  ,
     & hx) + u(i1,i2-1,i3  ,hx) + u(i1-1,i2,i3  ,hx) + u(i1,i2,i3  ,
     & hx))
         v(i1,i2,i3,hy) =.125*(u(i1-1,i2-1,i3-1,hy) + u(i1,i2-1,i3-1,
     & hy) + u(i1-1,i2,i3-1,hy) + u(i1,i2,i3-1,hy) +u(i1-1,i2-1,i3  ,
     & hy) + u(i1,i2-1,i3  ,hy) + u(i1-1,i2,i3  ,hy) + u(i1,i2,i3  ,
     & hy))
         v(i1,i2,i3,hz) =.125*(u(i1-1,i2-1,i3-1,hz) + u(i1,i2-1,i3-1,
     & hz) + u(i1-1,i2,i3-1,hz) + u(i1,i2,i3-1,hz) +u(i1-1,i2-1,i3  ,
     & hz) + u(i1,i2-1,i3  ,hz) + u(i1-1,i2,i3  ,hz) + u(i1,i2,i3  ,
     & hz))
        end do
        end do
        end do

        ! Now fixup values on the boundaries that were not averaged correctly
        ! For the E-field only the normal component is wrong
        ! Eventually the corner points should become correct
        extra=0
        numberOfGhostPoints=0
         extra1a=extra
         extra1b=extra
         extra2a=extra
         extra2b=extra
         if( nd.eq.3 )then
           extra3a=extra
           extra3b=extra
         else
           extra3a=0
           extra3b=0
         end if
         if( boundaryCondition(0,0).lt.0 )then
           extra1a=max(0,extra1a) ! over-ride extra=-1 : assign ends in periodic directions
           extra1b=extra1a
         else
           if( boundaryCondition(0,0).eq.0 )then
             extra1a=numberOfGhostPoints  ! include interpolation points since we assign ghost points outside these
           end if
           if( boundaryCondition(1,0).eq.0 )then
             extra1b=numberOfGhostPoints
           end if
         end if
         if( boundaryCondition(0,1).lt.0 )then
          extra2a=max(0,extra2a) ! over-ride extra=-1 : assign ends in periodic directions
          extra2b=extra2a
         else
           if( boundaryCondition(0,1).eq.0 )then
             extra2a=numberOfGhostPoints
           end if
           if( boundaryCondition(1,1).eq.0 )then
             extra2b=numberOfGhostPoints
           end if
         end if
         if(  nd.eq.3 .and. boundaryCondition(0,2).lt.0 )then
          extra3a=max(0,extra3a) ! over-ride extra=-1 : assign ends in periodic directions
          extra3b=extra3a
         else
           if( boundaryCondition(0,2).eq.0 )then
             extra3a=numberOfGhostPoints
           end if
           if( boundaryCondition(1,2).eq.0 )then
             extra3b=numberOfGhostPoints
           end if
         end if
         do axis=0,nd-1
         do side=0,1
           n1a=gridIndexRange(0,0)-extra1a
           n1b=gridIndexRange(1,0)+extra1b
           n2a=gridIndexRange(0,1)-extra2a
           n2b=gridIndexRange(1,1)+extra2b
           n3a=gridIndexRange(0,2)-extra3a
           n3b=gridIndexRange(1,2)+extra3b
           if( axis.eq.0 )then
             n1a=gridIndexRange(side,axis)
             n1b=gridIndexRange(side,axis)
           else if( axis.eq.1 )then
             n2a=gridIndexRange(side,axis)
             n2b=gridIndexRange(side,axis)
           else
             n3a=gridIndexRange(side,axis)
             n3b=gridIndexRange(side,axis)
           end if
           is1=0
           is2=0
           is3=0
           if( axis.eq.0 )then
             is1=1-2*side
           else if( axis.eq.1 )then
             is2=1-2*side
           else if( axis.eq.2 )then
             is3=1-2*side
           else
             stop 5
           end if
           axisp1=mod(axis+1,nd)
           axisp2=mod(axis+2,nd)
        !    ! (js1,js2,js3) used to compute tangential derivatives
        !    js1=0
        !    js2=0
        !    js3=0
        !    if( axisp1.eq.0 )then
        !      js1=1-2*side
        !    else if( axisp1.eq.1 )then
        !      js2=1-2*side
        !    else if( axisp1.eq.2 )then
        !      js3=1-2*side
        !    else
        !      stop 5
        !    end if
        ! 
        !    ! (ks1,ks2,ks3) used to compute second tangential derivative
        !    ks1=0
        !    ks2=0
        !    ks3=0
        !    if( axisp2.eq.0 )then
        !      ks1=1-2*side
        !    else if( axisp2.eq.1 )then
        !      ks2=1-2*side
        !    else if( axisp2.eq.2 )then
        !      ks3=1-2*side
        !    else
        !      stop 5
        !    end if
         if( boundaryCondition(side,axis).gt.0 )then
          en = ex+axis
          do i3=n3a,n3b
          do i2=n2a,n2b
          do i1=n1a,n1b
           v(i1,i2,i3,en)=(2.*v(i1+is1,i2+is2,i3+is3,en)-v(i1+is1+is1,
     & i2+is2+is2,i3+is3+is3,en))
           v(i1,i2,i3,hx)=(2.*v(i1+is1,i2+is2,i3+is3,hx)-v(i1+is1+is1,
     & i2+is2+is2,i3+is3+is3,hx))
           v(i1,i2,i3,hy)=(2.*v(i1+is1,i2+is2,i3+is3,hy)-v(i1+is1+is1,
     & i2+is2+is2,i3+is3+is3,hy))
           v(i1,i2,i3,hz)=(2.*v(i1+is1,i2+is2,i3+is3,hz)-v(i1+is1+is1,
     & i2+is2+is2,i3+is3+is3,hz))
          end do
          end do
          end do
         end if
         end do
         end do
         ! reset these values
         n1a=gridIndexRange(0,0)
         n1b=gridIndexRange(1,0)
         n2a=gridIndexRange(0,1)
         n2b=gridIndexRange(1,1)
         n3a=gridIndexRange(0,2)
         n3b=gridIndexRange(1,2)

       end if ! nd


      end if

      if( option.eq.4 )then

       ! plot curl(E)

       nCurlE= ed ! save results here
       n0=nCurlE
       write(*,'(" mxYeeIC: compute curl( Eknown ) nCurlE=",i2)') 
     & nCurlE

       n1=n0+1
       n2=n1+1
       n3=n2+1
       n4=n3+1
       n5=n4+1
       ! compute the curl of the imaginary part of the known solution
       if( nd.eq.2 )then
        ! --- Compute curl( Eknown )
        ex0=ex+3
        ey0=ey+3
        do i3=n3a,n3b
        do i2=n2a,n2b
        do i1=n1a,n1b

        v(i1,i2,i3,n0) =  ( uKnown(i1+1,i2,i3,ey0)-uKnown(i1,i2,i3,ey0)
     &  )/dx(0)-( uKnown(i1,i2+1,i3,ex0)-uKnown(i1,i2,i3,ex0) )/dx(1)

        end do
        end do
        end do

       else
        ! --- Compute curl( Eknown )
        ! compute the curl of the imaginary part of the known solution
        fact = -1./(twoPi*cc)
        ex0=ex+6
        ey0=ey+6
        ez0=ez+6
        do i3=n3a,n3b
        do i2=n2a,n2b
        do i1=n1a,n1b

        v(i1,i2,i3,n0) =(  ( uKnown(i1,i2+1,i3,ez )-uKnown(i1,i2,i3,ez 
     & ) )/dx(1)-( uKnown(i1,i2,i3+1,ey )-uKnown(i1,i2,i3,ey ) )/dx(2)
     & )*fact

        v(i1,i2,i3,n1) =(  ( uKnown(i1,i2,i3+1,ex )-uKnown(i1,i2,i3,ex 
     & ) )/dx(2)-( uKnown(i1+1,i2,i3,ez )-uKnown(i1,i2,i3,ez ) )/dx(0)
     & )*fact

        v(i1,i2,i3,n2) =(  ( uKnown(i1+1,i2,i3,ey )-uKnown(i1,i2,i3,ey 
     & ) )/dx(0)-( uKnown(i1,i2+1,i3,ex )-uKnown(i1,i2,i3,ex ) )/dx(1)
     & )*fact

        v(i1,i2,i3,n3) =(  ( uKnown(i1,i2+1,i3,ez0)-uKnown(i1,i2,i3,
     & ez0) )/dx(1)-( uKnown(i1,i2,i3+1,ey0)-uKnown(i1,i2,i3,ey0) )
     & /dx(2))*fact

        v(i1,i2,i3,n4) =(  ( uKnown(i1,i2,i3+1,ex0)-uKnown(i1,i2,i3,
     & ex0) )/dx(2)-( uKnown(i1+1,i2,i3,ez0)-uKnown(i1,i2,i3,ez0) )
     & /dx(0))*fact

        v(i1,i2,i3,n5) =(  ( uKnown(i1+1,i2,i3,ey0)-uKnown(i1,i2,i3,
     & ey0) )/dx(0)-( uKnown(i1,i2+1,i3,ex0)-uKnown(i1,i2,i3,ex0) )
     & /dx(1))*fact

        end do
        end do
        end do
       end if ! end nd


      end if

      return
      end
