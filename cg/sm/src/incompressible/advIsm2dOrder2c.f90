! This file automatically generated from advOptIsm.bf90 with bpp.
        subroutine advIsm2dOrder2c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,rsxy,xy,um,u,un,f,vt,vn,v,v1,v2,v3,v4,vt1,vt2,vt3,vt4,gridIndexRange,boundaryCondition,ipar,rpar,ierr )
       !======================================================================
       !   Advance a time step for the equations of Solid Mechanics (linear elasticity for now)
       ! 
       ! nd : number of space dimensions
       !
       !  option:  0  : update solution modified equation
       !           1  : add upwind dissipation
       !           2  : compute v.t for MOL schemes
       !
       !======================================================================
         implicit none
          integer nd, n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b
          real um(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
          real u(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
          real un(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
          real f(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
          real vt(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
          real vn(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
          real  v(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
          real v1(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
          real v2(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
          real v3(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
          real v4(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
          real vt1(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
          real vt2(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
          real vt3(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
          real vt4(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
          real rsxy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1,0:nd-1)
          real xy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1)
          integer mask(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
          integer boundaryCondition(0:1,0:2),gridIndexRange(0:1,0:2),ierr
          integer ipar(0:*)
          real rpar(0:*)
        !     ---- local variables -----
        ! real v(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)     ! ************************************ LOCALLY ALLOCATED
         real uTemp(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1)    ! ************************************ LOCALLY ALLOCATED
         integer c,i1,i2,i3,n,gridType,orderOfAccuracy,orderOfAccuracyInTime,debug,computeUt,dir
         integer addForcing,twilightZone,option,upwindSOS
         integer useWhereMask,useWhereMaskSave,grid,useVariableDissipation
         integer useConservative,combineDissipationWithAdvance
         integer u1c,u2c,u3c,pc
         integer v1c,v2c,v3c
         integer materialFormat,myid
         logical useLowerOrderUpwindingOnBoundaries
         integer useSosupDissipation,preComputeUpwindUt,correction
         real adSosup,uDotFactor,sosupParameter, adxSosup(0:2) 
         real adSosup2,adSosup4,adSosup6,upwindScaleFactor 
         integer m1a,m1b,m2a,m2b,m3a,m3b,ec,numGhost
         integer side,axis,axisp1,axisp2,is1,is2,is3,is
         integer nn1a,nn1b,nn2a,nn2b,nn3a,nn3b
         real cs,dtSqByRho,rhoi
         real cmol1,cmol2,cmol3,cmol4  
         real cld,cldamp  
         real cc,dt,dy,dz,cdt,cdtdx,cdtdy,cdtdz,adc,adcdt,add,adddt,dtOld,cu,cum
         real dt4by12
         real kx,ky,kz
         real t,ep
         real dx(0:2),dr(0:2)
         real ue,ve,pe
         real uett,uexx,ueyy,uezz
         real vett,vexx,veyy,vezz
         real wett,wexx,weyy,wezz
         real pex, pey, pez
         real time0,time1
         real c0,csq,dtsq,cdtsq,cdtsq12,lap(0:20)
         ! boundary conditions parameters and interfaceType values
! define BC parameters for fortran routines
! boundary conditions
!123456789012345678901234567890123456789012345678901234567890123456789
      integer interpolation,displacementBC,tractionBC
      integer slipWall,symmetry,interfaceBC
      integer abcEM2,abcPML,abc3,abc4,abc5,rbcNonLocal,rbcLocal,lastBC
      integer dirichletBoundaryCondition
      parameter( interpolation=0,displacementBC=1,tractionBC=2)
      parameter( slipWall=3,symmetry=4 )
      parameter( interfaceBC=5,abcEM2=6,abcPML=7,abc3=8,abc4=9 )
      parameter( abc5=10,rbcNonLocal=11,rbcLocal=12 )
      parameter( dirichletBoundaryCondition=13 )
      parameter( lastBC=14 )
! define interfaceType values for fortran routines
      integer noInterface                     ! no interface conditions are imposed
      integer heatFluxInterface               ! [ T.n ] = g
      integer tractionInterface               ! [ n.tau ] = g 
      integer tractionAndHeatFluxInterface
      parameter( noInterface=0, heatFluxInterface=1 )
      parameter( tractionInterface=2,tractionAndHeatFluxInterface=3 )
         integer rectangular,curvilinear
         parameter( rectangular=0, curvilinear=1 )
         integer timeSteppingMethod
         integer defaultTimeStepping,adamsSymmetricOrder3,rungeKuttaFourthOrder,stoermerTimeStepping,modifiedEquationTimeStepping
         parameter(defaultTimeStepping=0,adamsSymmetricOrder3=1,rungeKuttaFourthOrder=2,stoermerTimeStepping=3,modifiedEquationTimeStepping=4)
        !...........start statement function
         integer kd,m
         real rx,ry,rz,sx,sy,sz,tx,ty,tz
        ! include 'declareDiffOrder2f.h'
        ! include 'declareDiffOrder4f.h'
          real d12
          real d22
          real h12
          real h22
          real rxr2
          real rxs2
          real rxt2
          real rxrr2
          real rxss2
          real rxrs2
          real ryr2
          real rys2
          real ryt2
          real ryrr2
          real ryss2
          real ryrs2
          real rzr2
          real rzs2
          real rzt2
          real rzrr2
          real rzss2
          real rzrs2
          real sxr2
          real sxs2
          real sxt2
          real sxrr2
          real sxss2
          real sxrs2
          real syr2
          real sys2
          real syt2
          real syrr2
          real syss2
          real syrs2
          real szr2
          real szs2
          real szt2
          real szrr2
          real szss2
          real szrs2
          real txr2
          real txs2
          real txt2
          real txrr2
          real txss2
          real txrs2
          real tyr2
          real tys2
          real tyt2
          real tyrr2
          real tyss2
          real tyrs2
          real tzr2
          real tzs2
          real tzt2
          real tzrr2
          real tzss2
          real tzrs2
          real rxx21
          real rxx22
          real rxy22
          real rxx23
          real rxy23
          real rxz23
          real ryx22
          real ryy22
          real ryx23
          real ryy23
          real ryz23
          real rzx22
          real rzy22
          real rzx23
          real rzy23
          real rzz23
          real sxx22
          real sxy22
          real sxx23
          real sxy23
          real sxz23
          real syx22
          real syy22
          real syx23
          real syy23
          real syz23
          real szx22
          real szy22
          real szx23
          real szy23
          real szz23
          real txx22
          real txy22
          real txx23
          real txy23
          real txz23
          real tyx22
          real tyy22
          real tyx23
          real tyy23
          real tyz23
          real tzx22
          real tzy22
          real tzx23
          real tzy23
          real tzz23
          real ur2
          real us2
          real ut2
          real urr2
          real uss2
          real urs2
          real utt2
          real urt2
          real ust2
          real urrr2
          real usss2
          real uttt2
          real ux21
          real uy21
          real uz21
          real ux22
          real uy22
          real uz22
          real ux23
          real uy23
          real uz23
          real uxx21
          real uyy21
          real uxy21
          real uxz21
          real uyz21
          real uzz21
          real ulaplacian21
          real uxx22
          real uyy22
          real uxy22
          real uxz22
          real uyz22
          real uzz22
          real ulaplacian22
          real uxx23
          real uyy23
          real uzz23
          real uxy23
          real uxz23
          real uyz23
          real ulaplacian23
          real ux23r
          real uy23r
          real uz23r
          real uxx23r
          real uyy23r
          real uxy23r
          real uzz23r
          real uxz23r
          real uyz23r
          real ux21r
          real uy21r
          real uz21r
          real uxx21r
          real uyy21r
          real uzz21r
          real uxy21r
          real uxz21r
          real uyz21r
          real ulaplacian21r
          real ux22r
          real uy22r
          real uz22r
          real uxx22r
          real uyy22r
          real uzz22r
          real uxy22r
          real uxz22r
          real uyz22r
          real ulaplacian22r
          real ulaplacian23r
          real uxxx22r
          real uyyy22r
          real uxxy22r
          real uxyy22r
          real uxxxx22r
          real uyyyy22r
          real uxxyy22r
          real uxxx23r
          real uyyy23r
          real uzzz23r
          real uxxy23r
          real uxxz23r
          real uxyy23r
          real uyyz23r
          real uxzz23r
          real uyzz23r
          real uxxxx23r
          real uyyyy23r
          real uzzzz23r
          real uxxyy23r
          real uxxzz23r
          real uyyzz23r
          real uLapSq22r
          real uLapSq23r
          real unr2
          real uns2
          real unt2
          real unrr2
          real unss2
          real unrs2
          real untt2
          real unrt2
          real unst2
          real unrrr2
          real unsss2
          real unttt2
          real unx21
          real uny21
          real unz21
          real unx22
          real uny22
          real unz22
          real unx23
          real uny23
          real unz23
          real unxx21
          real unyy21
          real unxy21
          real unxz21
          real unyz21
          real unzz21
          real unlaplacian21
          real unxx22
          real unyy22
          real unxy22
          real unxz22
          real unyz22
          real unzz22
          real unlaplacian22
          real unxx23
          real unyy23
          real unzz23
          real unxy23
          real unxz23
          real unyz23
          real unlaplacian23
          real unx23r
          real uny23r
          real unz23r
          real unxx23r
          real unyy23r
          real unxy23r
          real unzz23r
          real unxz23r
          real unyz23r
          real unx21r
          real uny21r
          real unz21r
          real unxx21r
          real unyy21r
          real unzz21r
          real unxy21r
          real unxz21r
          real unyz21r
          real unlaplacian21r
          real unx22r
          real uny22r
          real unz22r
          real unxx22r
          real unyy22r
          real unzz22r
          real unxy22r
          real unxz22r
          real unyz22r
          real unlaplacian22r
          real unlaplacian23r
          real unxxx22r
          real unyyy22r
          real unxxy22r
          real unxyy22r
          real unxxxx22r
          real unyyyy22r
          real unxxyy22r
          real unxxx23r
          real unyyy23r
          real unzzz23r
          real unxxy23r
          real unxxz23r
          real unxyy23r
          real unyyz23r
          real unxzz23r
          real unyzz23r
          real unxxxx23r
          real unyyyy23r
          real unzzzz23r
          real unxxyy23r
          real unxxzz23r
          real unyyzz23r
          real unLapSq22r
          real unLapSq23r
          real umr2
          real ums2
          real umt2
          real umrr2
          real umss2
          real umrs2
          real umtt2
          real umrt2
          real umst2
          real umrrr2
          real umsss2
          real umttt2
          real umx21
          real umy21
          real umz21
          real umx22
          real umy22
          real umz22
          real umx23
          real umy23
          real umz23
          real umxx21
          real umyy21
          real umxy21
          real umxz21
          real umyz21
          real umzz21
          real umlaplacian21
          real umxx22
          real umyy22
          real umxy22
          real umxz22
          real umyz22
          real umzz22
          real umlaplacian22
          real umxx23
          real umyy23
          real umzz23
          real umxy23
          real umxz23
          real umyz23
          real umlaplacian23
          real umx23r
          real umy23r
          real umz23r
          real umxx23r
          real umyy23r
          real umxy23r
          real umzz23r
          real umxz23r
          real umyz23r
          real umx21r
          real umy21r
          real umz21r
          real umxx21r
          real umyy21r
          real umzz21r
          real umxy21r
          real umxz21r
          real umyz21r
          real umlaplacian21r
          real umx22r
          real umy22r
          real umz22r
          real umxx22r
          real umyy22r
          real umzz22r
          real umxy22r
          real umxz22r
          real umyz22r
          real umlaplacian22r
          real umlaplacian23r
          real umxxx22r
          real umyyy22r
          real umxxy22r
          real umxyy22r
          real umxxxx22r
          real umyyyy22r
          real umxxyy22r
          real umxxx23r
          real umyyy23r
          real umzzz23r
          real umxxy23r
          real umxxz23r
          real umxyy23r
          real umyyz23r
          real umxzz23r
          real umyzz23r
          real umxxxx23r
          real umyyyy23r
          real umzzzz23r
          real umxxyy23r
          real umxxzz23r
          real umyyzz23r
          real umLapSq22r
          real umLapSq23r
         ! declareDifferenceOrder2(v,none)
          real d14
          real d24
          real h41
          real h42
          real rxr4
          real rxs4
          real rxt4
          real ryr4
          real rys4
          real ryt4
          real rzr4
          real rzs4
          real rzt4
          real sxr4
          real sxs4
          real sxt4
          real syr4
          real sys4
          real syt4
          real szr4
          real szs4
          real szt4
          real txr4
          real txs4
          real txt4
          real tyr4
          real tys4
          real tyt4
          real tzr4
          real tzs4
          real tzt4
          real rxx41
          real rxx42
          real rxy42
          real rxx43
          real rxy43
          real rxz43
          real ryx42
          real ryy42
          real ryx43
          real ryy43
          real ryz43
          real rzx42
          real rzy42
          real rzx43
          real rzy43
          real rzz43
          real sxx42
          real sxy42
          real sxx43
          real sxy43
          real sxz43
          real syx42
          real syy42
          real syx43
          real syy43
          real syz43
          real szx42
          real szy42
          real szx43
          real szy43
          real szz43
          real txx42
          real txy42
          real txx43
          real txy43
          real txz43
          real tyx42
          real tyy42
          real tyx43
          real tyy43
          real tyz43
          real tzx42
          real tzy42
          real tzx43
          real tzy43
          real tzz43
          real ur4
          real us4
          real ut4
          real urr4
          real uss4
          real utt4
          real urs4
          real urt4
          real ust4
          real ux41
          real uy41
          real uz41
          real ux42
          real uy42
          real uz42
          real ux43
          real uy43
          real uz43
          real uxx41
          real uyy41
          real uxy41
          real uxz41
          real uyz41
          real uzz41
          real ulaplacian41
          real uxx42
          real uyy42
          real uxy42
          real uxz42
          real uyz42
          real uzz42
          real ulaplacian42
          real uxx43
          real uyy43
          real uzz43
          real uxy43
          real uxz43
          real uyz43
          real ulaplacian43
          real ux43r
          real uy43r
          real uz43r
          real uxx43r
          real uyy43r
          real uzz43r
          real uxy43r
          real uxz43r
          real uyz43r
          real ux41r
          real uy41r
          real uz41r
          real uxx41r
          real uyy41r
          real uzz41r
          real uxy41r
          real uxz41r
          real uyz41r
          real ulaplacian41r
          real ux42r
          real uy42r
          real uz42r
          real uxx42r
          real uyy42r
          real uzz42r
          real uxy42r
          real uxz42r
          real uyz42r
          real ulaplacian42r
          real ulaplacian43r
          real unr4
          real uns4
          real unt4
          real unrr4
          real unss4
          real untt4
          real unrs4
          real unrt4
          real unst4
          real unx41
          real uny41
          real unz41
          real unx42
          real uny42
          real unz42
          real unx43
          real uny43
          real unz43
          real unxx41
          real unyy41
          real unxy41
          real unxz41
          real unyz41
          real unzz41
          real unlaplacian41
          real unxx42
          real unyy42
          real unxy42
          real unxz42
          real unyz42
          real unzz42
          real unlaplacian42
          real unxx43
          real unyy43
          real unzz43
          real unxy43
          real unxz43
          real unyz43
          real unlaplacian43
          real unx43r
          real uny43r
          real unz43r
          real unxx43r
          real unyy43r
          real unzz43r
          real unxy43r
          real unxz43r
          real unyz43r
          real unx41r
          real uny41r
          real unz41r
          real unxx41r
          real unyy41r
          real unzz41r
          real unxy41r
          real unxz41r
          real unyz41r
          real unlaplacian41r
          real unx42r
          real uny42r
          real unz42r
          real unxx42r
          real unyy42r
          real unzz42r
          real unxy42r
          real unxz42r
          real unyz42r
          real unlaplacian42r
          real unlaplacian43r
          real umr4
          real ums4
          real umt4
          real umrr4
          real umss4
          real umtt4
          real umrs4
          real umrt4
          real umst4
          real umx41
          real umy41
          real umz41
          real umx42
          real umy42
          real umz42
          real umx43
          real umy43
          real umz43
          real umxx41
          real umyy41
          real umxy41
          real umxz41
          real umyz41
          real umzz41
          real umlaplacian41
          real umxx42
          real umyy42
          real umxy42
          real umxz42
          real umyz42
          real umzz42
          real umlaplacian42
          real umxx43
          real umyy43
          real umzz43
          real umxy43
          real umxz43
          real umyz43
          real umlaplacian43
          real umx43r
          real umy43r
          real umz43r
          real umxx43r
          real umyy43r
          real umzz43r
          real umxy43r
          real umxz43r
          real umyz43r
          real umx41r
          real umy41r
          real umz41r
          real umxx41r
          real umyy41r
          real umzz41r
          real umxy41r
          real umxz41r
          real umyz41r
          real umlaplacian41r
          real umx42r
          real umy42r
          real umz42r
          real umxx42r
          real umyy42r
          real umzz42r
          real umxy42r
          real umxz42r
          real umyz42r
          real umlaplacian42r
          real umlaplacian43r
         real mu,rho
         ! real cdt4by360,cdt6by20160
         ! real lap2d2,lap3d2,lap2d4,lap3d4,lap2d6,lap3d6,lap2d8,lap3d8,lap2d2Pow2,lap3d2Pow2,lap2d2Pow3,lap3d2Pow3,!      lap2d2Pow4,lap3d2Pow4,lap2d4Pow2,lap3d4Pow2,lap2d4Pow3,lap3d4Pow3,lap2d6Pow2,lap3d6Pow2
         ! real du,fd22d,fd23d,fd42d,fd43d,fd62d,fd63d,fd82d,fd83d
         real DztU
        ! real unxx22r,unyy22r,unxy22r,unx22r
        !.......statement functions for jacobian
         rx(i1,i2,i3)=rsxy(i1,i2,i3,0,0)
         ry(i1,i2,i3)=rsxy(i1,i2,i3,0,1)
         rz(i1,i2,i3)=rsxy(i1,i2,i3,0,2)
         sx(i1,i2,i3)=rsxy(i1,i2,i3,1,0)
         sy(i1,i2,i3)=rsxy(i1,i2,i3,1,1)
         sz(i1,i2,i3)=rsxy(i1,i2,i3,1,2)
         tx(i1,i2,i3)=rsxy(i1,i2,i3,2,0)
         ty(i1,i2,i3)=rsxy(i1,i2,i3,2,1)
         tz(i1,i2,i3)=rsxy(i1,i2,i3,2,2)
        ! D-zero in time (really undivided)
          DztU(i1,i2,i3,n) = (un(i1,i2,i3,n)-um(i1,i2,i3,n)) 
        !     The next macro call will define the difference approximation statement functions
         d12(kd) = 1./(2.*dr(kd))
         d22(kd) = 1./(dr(kd)**2)
         ur2(i1,i2,i3,kd)=(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))*d12(0)
         us2(i1,i2,i3,kd)=(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))*d12(1)
         ut2(i1,i2,i3,kd)=(u(i1,i2,i3+1,kd)-u(i1,i2,i3-1,kd))*d12(2)
         urr2(i1,i2,i3,kd)=(-2.*u(i1,i2,i3,kd)+(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd)) )*d22(0)
         uss2(i1,i2,i3,kd)=(-2.*u(i1,i2,i3,kd)+(u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd)) )*d22(1)
         urs2(i1,i2,i3,kd)=(ur2(i1,i2+1,i3,kd)-ur2(i1,i2-1,i3,kd))*d12(1)
         utt2(i1,i2,i3,kd)=(-2.*u(i1,i2,i3,kd)+(u(i1,i2,i3+1,kd)+u(i1,i2,i3-1,kd)) )*d22(2)
         urt2(i1,i2,i3,kd)=(ur2(i1,i2,i3+1,kd)-ur2(i1,i2,i3-1,kd))*d12(2)
         ust2(i1,i2,i3,kd)=(us2(i1,i2,i3+1,kd)-us2(i1,i2,i3-1,kd))*d12(2)
         urrr2(i1,i2,i3,kd)=(-2.*(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))+(u(i1+2,i2,i3,kd)-u(i1-2,i2,i3,kd)) )*d22(0)*d12(0)
         usss2(i1,i2,i3,kd)=(-2.*(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))+(u(i1,i2+2,i3,kd)-u(i1,i2-2,i3,kd)) )*d22(1)*d12(1)
         uttt2(i1,i2,i3,kd)=(-2.*(u(i1,i2,i3+1,kd)-u(i1,i2,i3-1,kd))+(u(i1,i2,i3+2,kd)-u(i1,i2,i3-2,kd)) )*d22(2)*d12(2)
         rxr2(i1,i2,i3)=(rx(i1+1,i2,i3)-rx(i1-1,i2,i3))*d12(0)
         rxs2(i1,i2,i3)=(rx(i1,i2+1,i3)-rx(i1,i2-1,i3))*d12(1)
         rxt2(i1,i2,i3)=(rx(i1,i2,i3+1)-rx(i1,i2,i3-1))*d12(2)
         rxrr2(i1,i2,i3)=(-2.*rx(i1,i2,i3)+(rx(i1+1,i2,i3)+rx(i1-1,i2,i3)) )*d22(0)
         rxss2(i1,i2,i3)=(-2.*rx(i1,i2,i3)+(rx(i1,i2+1,i3)+rx(i1,i2-1,i3)) )*d22(1)
         rxrs2(i1,i2,i3)=(rxr2(i1,i2+1,i3)-rxr2(i1,i2-1,i3))*d12(1)
         ryr2(i1,i2,i3)=(ry(i1+1,i2,i3)-ry(i1-1,i2,i3))*d12(0)
         rys2(i1,i2,i3)=(ry(i1,i2+1,i3)-ry(i1,i2-1,i3))*d12(1)
         ryt2(i1,i2,i3)=(ry(i1,i2,i3+1)-ry(i1,i2,i3-1))*d12(2)
         ryrr2(i1,i2,i3)=(-2.*ry(i1,i2,i3)+(ry(i1+1,i2,i3)+ry(i1-1,i2,i3)) )*d22(0)
         ryss2(i1,i2,i3)=(-2.*ry(i1,i2,i3)+(ry(i1,i2+1,i3)+ry(i1,i2-1,i3)) )*d22(1)
         ryrs2(i1,i2,i3)=(ryr2(i1,i2+1,i3)-ryr2(i1,i2-1,i3))*d12(1)
         rzr2(i1,i2,i3)=(rz(i1+1,i2,i3)-rz(i1-1,i2,i3))*d12(0)
         rzs2(i1,i2,i3)=(rz(i1,i2+1,i3)-rz(i1,i2-1,i3))*d12(1)
         rzt2(i1,i2,i3)=(rz(i1,i2,i3+1)-rz(i1,i2,i3-1))*d12(2)
         rzrr2(i1,i2,i3)=(-2.*rz(i1,i2,i3)+(rz(i1+1,i2,i3)+rz(i1-1,i2,i3)) )*d22(0)
         rzss2(i1,i2,i3)=(-2.*rz(i1,i2,i3)+(rz(i1,i2+1,i3)+rz(i1,i2-1,i3)) )*d22(1)
         rzrs2(i1,i2,i3)=(rzr2(i1,i2+1,i3)-rzr2(i1,i2-1,i3))*d12(1)
         sxr2(i1,i2,i3)=(sx(i1+1,i2,i3)-sx(i1-1,i2,i3))*d12(0)
         sxs2(i1,i2,i3)=(sx(i1,i2+1,i3)-sx(i1,i2-1,i3))*d12(1)
         sxt2(i1,i2,i3)=(sx(i1,i2,i3+1)-sx(i1,i2,i3-1))*d12(2)
         sxrr2(i1,i2,i3)=(-2.*sx(i1,i2,i3)+(sx(i1+1,i2,i3)+sx(i1-1,i2,i3)) )*d22(0)
         sxss2(i1,i2,i3)=(-2.*sx(i1,i2,i3)+(sx(i1,i2+1,i3)+sx(i1,i2-1,i3)) )*d22(1)
         sxrs2(i1,i2,i3)=(sxr2(i1,i2+1,i3)-sxr2(i1,i2-1,i3))*d12(1)
         syr2(i1,i2,i3)=(sy(i1+1,i2,i3)-sy(i1-1,i2,i3))*d12(0)
         sys2(i1,i2,i3)=(sy(i1,i2+1,i3)-sy(i1,i2-1,i3))*d12(1)
         syt2(i1,i2,i3)=(sy(i1,i2,i3+1)-sy(i1,i2,i3-1))*d12(2)
         syrr2(i1,i2,i3)=(-2.*sy(i1,i2,i3)+(sy(i1+1,i2,i3)+sy(i1-1,i2,i3)) )*d22(0)
         syss2(i1,i2,i3)=(-2.*sy(i1,i2,i3)+(sy(i1,i2+1,i3)+sy(i1,i2-1,i3)) )*d22(1)
         syrs2(i1,i2,i3)=(syr2(i1,i2+1,i3)-syr2(i1,i2-1,i3))*d12(1)
         szr2(i1,i2,i3)=(sz(i1+1,i2,i3)-sz(i1-1,i2,i3))*d12(0)
         szs2(i1,i2,i3)=(sz(i1,i2+1,i3)-sz(i1,i2-1,i3))*d12(1)
         szt2(i1,i2,i3)=(sz(i1,i2,i3+1)-sz(i1,i2,i3-1))*d12(2)
         szrr2(i1,i2,i3)=(-2.*sz(i1,i2,i3)+(sz(i1+1,i2,i3)+sz(i1-1,i2,i3)) )*d22(0)
         szss2(i1,i2,i3)=(-2.*sz(i1,i2,i3)+(sz(i1,i2+1,i3)+sz(i1,i2-1,i3)) )*d22(1)
         szrs2(i1,i2,i3)=(szr2(i1,i2+1,i3)-szr2(i1,i2-1,i3))*d12(1)
         txr2(i1,i2,i3)=(tx(i1+1,i2,i3)-tx(i1-1,i2,i3))*d12(0)
         txs2(i1,i2,i3)=(tx(i1,i2+1,i3)-tx(i1,i2-1,i3))*d12(1)
         txt2(i1,i2,i3)=(tx(i1,i2,i3+1)-tx(i1,i2,i3-1))*d12(2)
         txrr2(i1,i2,i3)=(-2.*tx(i1,i2,i3)+(tx(i1+1,i2,i3)+tx(i1-1,i2,i3)) )*d22(0)
         txss2(i1,i2,i3)=(-2.*tx(i1,i2,i3)+(tx(i1,i2+1,i3)+tx(i1,i2-1,i3)) )*d22(1)
         txrs2(i1,i2,i3)=(txr2(i1,i2+1,i3)-txr2(i1,i2-1,i3))*d12(1)
         tyr2(i1,i2,i3)=(ty(i1+1,i2,i3)-ty(i1-1,i2,i3))*d12(0)
         tys2(i1,i2,i3)=(ty(i1,i2+1,i3)-ty(i1,i2-1,i3))*d12(1)
         tyt2(i1,i2,i3)=(ty(i1,i2,i3+1)-ty(i1,i2,i3-1))*d12(2)
         tyrr2(i1,i2,i3)=(-2.*ty(i1,i2,i3)+(ty(i1+1,i2,i3)+ty(i1-1,i2,i3)) )*d22(0)
         tyss2(i1,i2,i3)=(-2.*ty(i1,i2,i3)+(ty(i1,i2+1,i3)+ty(i1,i2-1,i3)) )*d22(1)
         tyrs2(i1,i2,i3)=(tyr2(i1,i2+1,i3)-tyr2(i1,i2-1,i3))*d12(1)
         tzr2(i1,i2,i3)=(tz(i1+1,i2,i3)-tz(i1-1,i2,i3))*d12(0)
         tzs2(i1,i2,i3)=(tz(i1,i2+1,i3)-tz(i1,i2-1,i3))*d12(1)
         tzt2(i1,i2,i3)=(tz(i1,i2,i3+1)-tz(i1,i2,i3-1))*d12(2)
         tzrr2(i1,i2,i3)=(-2.*tz(i1,i2,i3)+(tz(i1+1,i2,i3)+tz(i1-1,i2,i3)) )*d22(0)
         tzss2(i1,i2,i3)=(-2.*tz(i1,i2,i3)+(tz(i1,i2+1,i3)+tz(i1,i2-1,i3)) )*d22(1)
         tzrs2(i1,i2,i3)=(tzr2(i1,i2+1,i3)-tzr2(i1,i2-1,i3))*d12(1)
         ux21(i1,i2,i3,kd)= rx(i1,i2,i3)*ur2(i1,i2,i3,kd)
         uy21(i1,i2,i3,kd)=0
         uz21(i1,i2,i3,kd)=0
         ux22(i1,i2,i3,kd)= rx(i1,i2,i3)*ur2(i1,i2,i3,kd)+sx(i1,i2,i3)*us2(i1,i2,i3,kd)
         uy22(i1,i2,i3,kd)= ry(i1,i2,i3)*ur2(i1,i2,i3,kd)+sy(i1,i2,i3)*us2(i1,i2,i3,kd)
         uz22(i1,i2,i3,kd)=0
         ux23(i1,i2,i3,kd)=rx(i1,i2,i3)*ur2(i1,i2,i3,kd)+sx(i1,i2,i3)*us2(i1,i2,i3,kd)+tx(i1,i2,i3)*ut2(i1,i2,i3,kd)
         uy23(i1,i2,i3,kd)=ry(i1,i2,i3)*ur2(i1,i2,i3,kd)+sy(i1,i2,i3)*us2(i1,i2,i3,kd)+ty(i1,i2,i3)*ut2(i1,i2,i3,kd)
         uz23(i1,i2,i3,kd)=rz(i1,i2,i3)*ur2(i1,i2,i3,kd)+sz(i1,i2,i3)*us2(i1,i2,i3,kd)+tz(i1,i2,i3)*ut2(i1,i2,i3,kd)
         rxx21(i1,i2,i3)= rx(i1,i2,i3)*rxr2(i1,i2,i3)
         rxx22(i1,i2,i3)= rx(i1,i2,i3)*rxr2(i1,i2,i3)+sx(i1,i2,i3)*rxs2(i1,i2,i3)
         rxy22(i1,i2,i3)= ry(i1,i2,i3)*rxr2(i1,i2,i3)+sy(i1,i2,i3)*rxs2(i1,i2,i3)
         rxx23(i1,i2,i3)=rx(i1,i2,i3)*rxr2(i1,i2,i3)+sx(i1,i2,i3)*rxs2(i1,i2,i3)+tx(i1,i2,i3)*rxt2(i1,i2,i3)
         rxy23(i1,i2,i3)=ry(i1,i2,i3)*rxr2(i1,i2,i3)+sy(i1,i2,i3)*rxs2(i1,i2,i3)+ty(i1,i2,i3)*rxt2(i1,i2,i3)
         rxz23(i1,i2,i3)=rz(i1,i2,i3)*rxr2(i1,i2,i3)+sz(i1,i2,i3)*rxs2(i1,i2,i3)+tz(i1,i2,i3)*rxt2(i1,i2,i3)
         ryx22(i1,i2,i3)= rx(i1,i2,i3)*ryr2(i1,i2,i3)+sx(i1,i2,i3)*rys2(i1,i2,i3)
         ryy22(i1,i2,i3)= ry(i1,i2,i3)*ryr2(i1,i2,i3)+sy(i1,i2,i3)*rys2(i1,i2,i3)
         ryx23(i1,i2,i3)=rx(i1,i2,i3)*ryr2(i1,i2,i3)+sx(i1,i2,i3)*rys2(i1,i2,i3)+tx(i1,i2,i3)*ryt2(i1,i2,i3)
         ryy23(i1,i2,i3)=ry(i1,i2,i3)*ryr2(i1,i2,i3)+sy(i1,i2,i3)*rys2(i1,i2,i3)+ty(i1,i2,i3)*ryt2(i1,i2,i3)
         ryz23(i1,i2,i3)=rz(i1,i2,i3)*ryr2(i1,i2,i3)+sz(i1,i2,i3)*rys2(i1,i2,i3)+tz(i1,i2,i3)*ryt2(i1,i2,i3)
         rzx22(i1,i2,i3)= rx(i1,i2,i3)*rzr2(i1,i2,i3)+sx(i1,i2,i3)*rzs2(i1,i2,i3)
         rzy22(i1,i2,i3)= ry(i1,i2,i3)*rzr2(i1,i2,i3)+sy(i1,i2,i3)*rzs2(i1,i2,i3)
         rzx23(i1,i2,i3)=rx(i1,i2,i3)*rzr2(i1,i2,i3)+sx(i1,i2,i3)*rzs2(i1,i2,i3)+tx(i1,i2,i3)*rzt2(i1,i2,i3)
         rzy23(i1,i2,i3)=ry(i1,i2,i3)*rzr2(i1,i2,i3)+sy(i1,i2,i3)*rzs2(i1,i2,i3)+ty(i1,i2,i3)*rzt2(i1,i2,i3)
         rzz23(i1,i2,i3)=rz(i1,i2,i3)*rzr2(i1,i2,i3)+sz(i1,i2,i3)*rzs2(i1,i2,i3)+tz(i1,i2,i3)*rzt2(i1,i2,i3)
         sxx22(i1,i2,i3)= rx(i1,i2,i3)*sxr2(i1,i2,i3)+sx(i1,i2,i3)*sxs2(i1,i2,i3)
         sxy22(i1,i2,i3)= ry(i1,i2,i3)*sxr2(i1,i2,i3)+sy(i1,i2,i3)*sxs2(i1,i2,i3)
         sxx23(i1,i2,i3)=rx(i1,i2,i3)*sxr2(i1,i2,i3)+sx(i1,i2,i3)*sxs2(i1,i2,i3)+tx(i1,i2,i3)*sxt2(i1,i2,i3)
         sxy23(i1,i2,i3)=ry(i1,i2,i3)*sxr2(i1,i2,i3)+sy(i1,i2,i3)*sxs2(i1,i2,i3)+ty(i1,i2,i3)*sxt2(i1,i2,i3)
         sxz23(i1,i2,i3)=rz(i1,i2,i3)*sxr2(i1,i2,i3)+sz(i1,i2,i3)*sxs2(i1,i2,i3)+tz(i1,i2,i3)*sxt2(i1,i2,i3)
         syx22(i1,i2,i3)= rx(i1,i2,i3)*syr2(i1,i2,i3)+sx(i1,i2,i3)*sys2(i1,i2,i3)
         syy22(i1,i2,i3)= ry(i1,i2,i3)*syr2(i1,i2,i3)+sy(i1,i2,i3)*sys2(i1,i2,i3)
         syx23(i1,i2,i3)=rx(i1,i2,i3)*syr2(i1,i2,i3)+sx(i1,i2,i3)*sys2(i1,i2,i3)+tx(i1,i2,i3)*syt2(i1,i2,i3)
         syy23(i1,i2,i3)=ry(i1,i2,i3)*syr2(i1,i2,i3)+sy(i1,i2,i3)*sys2(i1,i2,i3)+ty(i1,i2,i3)*syt2(i1,i2,i3)
         syz23(i1,i2,i3)=rz(i1,i2,i3)*syr2(i1,i2,i3)+sz(i1,i2,i3)*sys2(i1,i2,i3)+tz(i1,i2,i3)*syt2(i1,i2,i3)
         szx22(i1,i2,i3)= rx(i1,i2,i3)*szr2(i1,i2,i3)+sx(i1,i2,i3)*szs2(i1,i2,i3)
         szy22(i1,i2,i3)= ry(i1,i2,i3)*szr2(i1,i2,i3)+sy(i1,i2,i3)*szs2(i1,i2,i3)
         szx23(i1,i2,i3)=rx(i1,i2,i3)*szr2(i1,i2,i3)+sx(i1,i2,i3)*szs2(i1,i2,i3)+tx(i1,i2,i3)*szt2(i1,i2,i3)
         szy23(i1,i2,i3)=ry(i1,i2,i3)*szr2(i1,i2,i3)+sy(i1,i2,i3)*szs2(i1,i2,i3)+ty(i1,i2,i3)*szt2(i1,i2,i3)
         szz23(i1,i2,i3)=rz(i1,i2,i3)*szr2(i1,i2,i3)+sz(i1,i2,i3)*szs2(i1,i2,i3)+tz(i1,i2,i3)*szt2(i1,i2,i3)
         txx22(i1,i2,i3)= rx(i1,i2,i3)*txr2(i1,i2,i3)+sx(i1,i2,i3)*txs2(i1,i2,i3)
         txy22(i1,i2,i3)= ry(i1,i2,i3)*txr2(i1,i2,i3)+sy(i1,i2,i3)*txs2(i1,i2,i3)
         txx23(i1,i2,i3)=rx(i1,i2,i3)*txr2(i1,i2,i3)+sx(i1,i2,i3)*txs2(i1,i2,i3)+tx(i1,i2,i3)*txt2(i1,i2,i3)
         txy23(i1,i2,i3)=ry(i1,i2,i3)*txr2(i1,i2,i3)+sy(i1,i2,i3)*txs2(i1,i2,i3)+ty(i1,i2,i3)*txt2(i1,i2,i3)
         txz23(i1,i2,i3)=rz(i1,i2,i3)*txr2(i1,i2,i3)+sz(i1,i2,i3)*txs2(i1,i2,i3)+tz(i1,i2,i3)*txt2(i1,i2,i3)
         tyx22(i1,i2,i3)= rx(i1,i2,i3)*tyr2(i1,i2,i3)+sx(i1,i2,i3)*tys2(i1,i2,i3)
         tyy22(i1,i2,i3)= ry(i1,i2,i3)*tyr2(i1,i2,i3)+sy(i1,i2,i3)*tys2(i1,i2,i3)
         tyx23(i1,i2,i3)=rx(i1,i2,i3)*tyr2(i1,i2,i3)+sx(i1,i2,i3)*tys2(i1,i2,i3)+tx(i1,i2,i3)*tyt2(i1,i2,i3)
         tyy23(i1,i2,i3)=ry(i1,i2,i3)*tyr2(i1,i2,i3)+sy(i1,i2,i3)*tys2(i1,i2,i3)+ty(i1,i2,i3)*tyt2(i1,i2,i3)
         tyz23(i1,i2,i3)=rz(i1,i2,i3)*tyr2(i1,i2,i3)+sz(i1,i2,i3)*tys2(i1,i2,i3)+tz(i1,i2,i3)*tyt2(i1,i2,i3)
         tzx22(i1,i2,i3)= rx(i1,i2,i3)*tzr2(i1,i2,i3)+sx(i1,i2,i3)*tzs2(i1,i2,i3)
         tzy22(i1,i2,i3)= ry(i1,i2,i3)*tzr2(i1,i2,i3)+sy(i1,i2,i3)*tzs2(i1,i2,i3)
         tzx23(i1,i2,i3)=rx(i1,i2,i3)*tzr2(i1,i2,i3)+sx(i1,i2,i3)*tzs2(i1,i2,i3)+tx(i1,i2,i3)*tzt2(i1,i2,i3)
         tzy23(i1,i2,i3)=ry(i1,i2,i3)*tzr2(i1,i2,i3)+sy(i1,i2,i3)*tzs2(i1,i2,i3)+ty(i1,i2,i3)*tzt2(i1,i2,i3)
         tzz23(i1,i2,i3)=rz(i1,i2,i3)*tzr2(i1,i2,i3)+sz(i1,i2,i3)*tzs2(i1,i2,i3)+tz(i1,i2,i3)*tzt2(i1,i2,i3)
         uxx21(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*urr2(i1,i2,i3,kd)+(rxx22(i1,i2,i3))*ur2(i1,i2,i3,kd)
         uyy21(i1,i2,i3,kd)=0
         uxy21(i1,i2,i3,kd)=0
         uxz21(i1,i2,i3,kd)=0
         uyz21(i1,i2,i3,kd)=0
         uzz21(i1,i2,i3,kd)=0
         ulaplacian21(i1,i2,i3,kd)=uxx21(i1,i2,i3,kd)
         uxx22(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*urr2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3))*urs2(i1,i2,i3,kd)+(sx(i1,i2,i3)**2)*uss2(i1,i2,i3,kd)+(rxx22(i1,i2,i3))*ur2(i1,i2,i3,kd)+(sxx22(i1,i2,i3))*us2(i1,i2,i3,kd)
         uyy22(i1,i2,i3,kd)=(ry(i1,i2,i3)**2)*urr2(i1,i2,i3,kd)+2.*(ry(i1,i2,i3)*sy(i1,i2,i3))*urs2(i1,i2,i3,kd)+(sy(i1,i2,i3)**2)*uss2(i1,i2,i3,kd)+(ryy22(i1,i2,i3))*ur2(i1,i2,i3,kd)+(syy22(i1,i2,i3))*us2(i1,i2,i3,kd)
         uxy22(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*urr2(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*urs2(i1,i2,i3,kd)+sx(i1,i2,i3)*sy(i1,i2,i3)*uss2(i1,i2,i3,kd)+rxy22(i1,i2,i3)*ur2(i1,i2,i3,kd)+sxy22(i1,i2,i3)*us2(i1,i2,i3,kd)
         uxz22(i1,i2,i3,kd)=0
         uyz22(i1,i2,i3,kd)=0
         uzz22(i1,i2,i3,kd)=0
         ulaplacian22(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2)*urr2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*sy(i1,i2,i3))*urs2(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**2)*uss2(i1,i2,i3,kd)+(rxx22(i1,i2,i3)+ryy22(i1,i2,i3))*ur2(i1,i2,i3,kd)+(sxx22(i1,i2,i3)+syy22(i1,i2,i3))*us2(i1,i2,i3,kd)
         uxx23(i1,i2,i3,kd)=rx(i1,i2,i3)**2*urr2(i1,i2,i3,kd)+sx(i1,i2,i3)**2*uss2(i1,i2,i3,kd)+tx(i1,i2,i3)**2*utt2(i1,i2,i3,kd)+2.*rx(i1,i2,i3)*sx(i1,i2,i3)*urs2(i1,i2,i3,kd)+2.*rx(i1,i2,i3)*tx(i1,i2,i3)*urt2(i1,i2,i3,kd)+2.*sx(i1,i2,i3)*tx(i1,i2,i3)*ust2(i1,i2,i3,kd)+rxx23(i1,i2,i3)*ur2(i1,i2,i3,kd)+sxx23(i1,i2,i3)*us2(i1,i2,i3,kd)+txx23(i1,i2,i3)*ut2(i1,i2,i3,kd)
         uyy23(i1,i2,i3,kd)=ry(i1,i2,i3)**2*urr2(i1,i2,i3,kd)+sy(i1,i2,i3)**2*uss2(i1,i2,i3,kd)+ty(i1,i2,i3)**2*utt2(i1,i2,i3,kd)+2.*ry(i1,i2,i3)*sy(i1,i2,i3)*urs2(i1,i2,i3,kd)+2.*ry(i1,i2,i3)*ty(i1,i2,i3)*urt2(i1,i2,i3,kd)+2.*sy(i1,i2,i3)*ty(i1,i2,i3)*ust2(i1,i2,i3,kd)+ryy23(i1,i2,i3)*ur2(i1,i2,i3,kd)+syy23(i1,i2,i3)*us2(i1,i2,i3,kd)+tyy23(i1,i2,i3)*ut2(i1,i2,i3,kd)
         uzz23(i1,i2,i3,kd)=rz(i1,i2,i3)**2*urr2(i1,i2,i3,kd)+sz(i1,i2,i3)**2*uss2(i1,i2,i3,kd)+tz(i1,i2,i3)**2*utt2(i1,i2,i3,kd)+2.*rz(i1,i2,i3)*sz(i1,i2,i3)*urs2(i1,i2,i3,kd)+2.*rz(i1,i2,i3)*tz(i1,i2,i3)*urt2(i1,i2,i3,kd)+2.*sz(i1,i2,i3)*tz(i1,i2,i3)*ust2(i1,i2,i3,kd)+rzz23(i1,i2,i3)*ur2(i1,i2,i3,kd)+szz23(i1,i2,i3)*us2(i1,i2,i3,kd)+tzz23(i1,i2,i3)*ut2(i1,i2,i3,kd)
         uxy23(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*urr2(i1,i2,i3,kd)+sx(i1,i2,i3)*sy(i1,i2,i3)*uss2(i1,i2,i3,kd)+tx(i1,i2,i3)*ty(i1,i2,i3)*utt2(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*urs2(i1,i2,i3,kd)+(rx(i1,i2,i3)*ty(i1,i2,i3)+ry(i1,i2,i3)*tx(i1,i2,i3))*urt2(i1,i2,i3,kd)+(sx(i1,i2,i3)*ty(i1,i2,i3)+sy(i1,i2,i3)*tx(i1,i2,i3))*ust2(i1,i2,i3,kd)+rxy23(i1,i2,i3)*ur2(i1,i2,i3,kd)+sxy23(i1,i2,i3)*us2(i1,i2,i3,kd)+txy23(i1,i2,i3)*ut2(i1,i2,i3,kd)
         uxz23(i1,i2,i3,kd)=rx(i1,i2,i3)*rz(i1,i2,i3)*urr2(i1,i2,i3,kd)+sx(i1,i2,i3)*sz(i1,i2,i3)*uss2(i1,i2,i3,kd)+tx(i1,i2,i3)*tz(i1,i2,i3)*utt2(i1,i2,i3,kd)+(rx(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,i2,i3)*sx(i1,i2,i3))*urs2(i1,i2,i3,kd)+(rx(i1,i2,i3)*tz(i1,i2,i3)+rz(i1,i2,i3)*tx(i1,i2,i3))*urt2(i1,i2,i3,kd)+(sx(i1,i2,i3)*tz(i1,i2,i3)+sz(i1,i2,i3)*tx(i1,i2,i3))*ust2(i1,i2,i3,kd)+rxz23(i1,i2,i3)*ur2(i1,i2,i3,kd)+sxz23(i1,i2,i3)*us2(i1,i2,i3,kd)+txz23(i1,i2,i3)*ut2(i1,i2,i3,kd)
         uyz23(i1,i2,i3,kd)=ry(i1,i2,i3)*rz(i1,i2,i3)*urr2(i1,i2,i3,kd)+sy(i1,i2,i3)*sz(i1,i2,i3)*uss2(i1,i2,i3,kd)+ty(i1,i2,i3)*tz(i1,i2,i3)*utt2(i1,i2,i3,kd)+(ry(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,i2,i3)*sy(i1,i2,i3))*urs2(i1,i2,i3,kd)+(ry(i1,i2,i3)*tz(i1,i2,i3)+rz(i1,i2,i3)*ty(i1,i2,i3))*urt2(i1,i2,i3,kd)+(sy(i1,i2,i3)*tz(i1,i2,i3)+sz(i1,i2,i3)*ty(i1,i2,i3))*ust2(i1,i2,i3,kd)+ryz23(i1,i2,i3)*ur2(i1,i2,i3,kd)+syz23(i1,i2,i3)*us2(i1,i2,i3,kd)+tyz23(i1,i2,i3)*ut2(i1,i2,i3,kd)
         ulaplacian23(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2+rz(i1,i2,i3)**2)*urr2(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**2+sz(i1,i2,i3)**2)*uss2(i1,i2,i3,kd)+(tx(i1,i2,i3)**2+ty(i1,i2,i3)**2+tz(i1,i2,i3)**2)*utt2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*sy(i1,i2,i3)+rz(i1,i2,i3)*sz(i1,i2,i3))*urs2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*tx(i1,i2,i3)+ ry(i1,i2,i3)*ty(i1,i2,i3)+rz(i1,i2,i3)*tz(i1,i2,i3))*urt2(i1,i2,i3,kd)+2.*(sx(i1,i2,i3)*tx(i1,i2,i3)+ sy(i1,i2,i3)*ty(i1,i2,i3)+sz(i1,i2,i3)*tz(i1,i2,i3))*ust2(i1,i2,i3,kd)+(rxx23(i1,i2,i3)+ryy23(i1,i2,i3)+rzz23(i1,i2,i3))*ur2(i1,i2,i3,kd)+(sxx23(i1,i2,i3)+syy23(i1,i2,i3)+szz23(i1,i2,i3))*us2(i1,i2,i3,kd)+(txx23(i1,i2,i3)+tyy23(i1,i2,i3)+tzz23(i1,i2,i3))*ut2(i1,i2,i3,kd)
         !============================================================================================
         ! Define derivatives for a rectangular grid
         !
         !============================================================================================
         h12(kd) = 1./(2.*dx(kd))
         h22(kd) = 1./(dx(kd)**2)
         ux23r(i1,i2,i3,kd)=(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))*h12(0)
         uy23r(i1,i2,i3,kd)=(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))*h12(1)
         uz23r(i1,i2,i3,kd)=(u(i1,i2,i3+1,kd)-u(i1,i2,i3-1,kd))*h12(2)
         uxx23r(i1,i2,i3,kd)=(-2.*u(i1,i2,i3,kd)+(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd)) )*h22(0)
         uyy23r(i1,i2,i3,kd)=(-2.*u(i1,i2,i3,kd)+(u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd)) )*h22(1)
         uxy23r(i1,i2,i3,kd)=(ux23r(i1,i2+1,i3,kd)-ux23r(i1,i2-1,i3,kd))*h12(1)
         uzz23r(i1,i2,i3,kd)=(-2.*u(i1,i2,i3,kd)+(u(i1,i2,i3+1,kd)+u(i1,i2,i3-1,kd)) )*h22(2)
         uxz23r(i1,i2,i3,kd)=(ux23r(i1,i2,i3+1,kd)-ux23r(i1,i2,i3-1,kd))*h12(2)
         uyz23r(i1,i2,i3,kd)=(uy23r(i1,i2,i3+1,kd)-uy23r(i1,i2,i3-1,kd))*h12(2)
         ux21r(i1,i2,i3,kd)= ux23r(i1,i2,i3,kd)
         uy21r(i1,i2,i3,kd)= uy23r(i1,i2,i3,kd)
         uz21r(i1,i2,i3,kd)= uz23r(i1,i2,i3,kd)
         uxx21r(i1,i2,i3,kd)= uxx23r(i1,i2,i3,kd)
         uyy21r(i1,i2,i3,kd)= uyy23r(i1,i2,i3,kd)
         uzz21r(i1,i2,i3,kd)= uzz23r(i1,i2,i3,kd)
         uxy21r(i1,i2,i3,kd)= uxy23r(i1,i2,i3,kd)
         uxz21r(i1,i2,i3,kd)= uxz23r(i1,i2,i3,kd)
         uyz21r(i1,i2,i3,kd)= uyz23r(i1,i2,i3,kd)
         ulaplacian21r(i1,i2,i3,kd)=uxx23r(i1,i2,i3,kd)
         ux22r(i1,i2,i3,kd)= ux23r(i1,i2,i3,kd)
         uy22r(i1,i2,i3,kd)= uy23r(i1,i2,i3,kd)
         uz22r(i1,i2,i3,kd)= uz23r(i1,i2,i3,kd)
         uxx22r(i1,i2,i3,kd)= uxx23r(i1,i2,i3,kd)
         uyy22r(i1,i2,i3,kd)= uyy23r(i1,i2,i3,kd)
         uzz22r(i1,i2,i3,kd)= uzz23r(i1,i2,i3,kd)
         uxy22r(i1,i2,i3,kd)= uxy23r(i1,i2,i3,kd)
         uxz22r(i1,i2,i3,kd)= uxz23r(i1,i2,i3,kd)
         uyz22r(i1,i2,i3,kd)= uyz23r(i1,i2,i3,kd)
         ulaplacian22r(i1,i2,i3,kd)=uxx23r(i1,i2,i3,kd)+uyy23r(i1,i2,i3,kd)
         ulaplacian23r(i1,i2,i3,kd)=uxx23r(i1,i2,i3,kd)+uyy23r(i1,i2,i3,kd)+uzz23r(i1,i2,i3,kd)
         uxxx22r(i1,i2,i3,kd)=(-2.*(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))+(u(i1+2,i2,i3,kd)-u(i1-2,i2,i3,kd)) )*h22(0)*h12(0)
         uyyy22r(i1,i2,i3,kd)=(-2.*(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))+(u(i1,i2+2,i3,kd)-u(i1,i2-2,i3,kd)) )*h22(1)*h12(1)
         uxxy22r(i1,i2,i3,kd)=( uxx22r(i1,i2+1,i3,kd)-uxx22r(i1,i2-1,i3,kd))/(2.*dx(1))
         uxyy22r(i1,i2,i3,kd)=( uyy22r(i1+1,i2,i3,kd)-uyy22r(i1-1,i2,i3,kd))/(2.*dx(0))
         uxxxx22r(i1,i2,i3,kd)=(6.*u(i1,i2,i3,kd)-4.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd))+(u(i1+2,i2,i3,kd)+u(i1-2,i2,i3,kd)) )/(dx(0)**4)
         uyyyy22r(i1,i2,i3,kd)=(6.*u(i1,i2,i3,kd)-4.*(u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd))+(u(i1,i2+2,i3,kd)+u(i1,i2-2,i3,kd)) )/(dx(1)**4)
         uxxyy22r(i1,i2,i3,kd)=( 4.*u(i1,i2,i3,kd)     -2.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd)+u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd))   +   (u(i1+1,i2+1,i3,kd)+u(i1-1,i2+1,i3,kd)+u(i1+1,i2-1,i3,kd)+u(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)
         ! 2D laplacian squared = u.xxxx + 2 u.xxyy + u.yyyy
         uLapSq22r(i1,i2,i3,kd)= ( 6.*u(i1,i2,i3,kd)   - 4.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd))    +(u(i1+2,i2,i3,kd)+u(i1-2,i2,i3,kd)) )/(dx(0)**4) +( 6.*u(i1,i2,i3,kd)    -4.*(u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd))    +(u(i1,i2+2,i3,kd)+u(i1,i2-2,i3,kd)) )/(dx(1)**4)  +( 8.*u(i1,i2,i3,kd)     -4.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd)+u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd))   +2.*(u(i1+1,i2+1,i3,kd)+u(i1-1,i2+1,i3,kd)+u(i1+1,i2-1,i3,kd)+u(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)
         uxxx23r(i1,i2,i3,kd)=(-2.*(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))+(u(i1+2,i2,i3,kd)-u(i1-2,i2,i3,kd)) )*h22(0)*h12(0)
         uyyy23r(i1,i2,i3,kd)=(-2.*(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))+(u(i1,i2+2,i3,kd)-u(i1,i2-2,i3,kd)) )*h22(1)*h12(1)
         uzzz23r(i1,i2,i3,kd)=(-2.*(u(i1,i2,i3+1,kd)-u(i1,i2,i3-1,kd))+(u(i1,i2,i3+2,kd)-u(i1,i2,i3-2,kd)) )*h22(1)*h12(2)
         uxxy23r(i1,i2,i3,kd)=( uxx22r(i1,i2+1,i3,kd)-uxx22r(i1,i2-1,i3,kd))/(2.*dx(1))
         uxyy23r(i1,i2,i3,kd)=( uyy22r(i1+1,i2,i3,kd)-uyy22r(i1-1,i2,i3,kd))/(2.*dx(0))
         uxxz23r(i1,i2,i3,kd)=( uxx22r(i1,i2,i3+1,kd)-uxx22r(i1,i2,i3-1,kd))/(2.*dx(2))
         uyyz23r(i1,i2,i3,kd)=( uyy22r(i1,i2,i3+1,kd)-uyy22r(i1,i2,i3-1,kd))/(2.*dx(2))
         uxzz23r(i1,i2,i3,kd)=( uzz22r(i1+1,i2,i3,kd)-uzz22r(i1-1,i2,i3,kd))/(2.*dx(0))
         uyzz23r(i1,i2,i3,kd)=( uzz22r(i1,i2+1,i3,kd)-uzz22r(i1,i2-1,i3,kd))/(2.*dx(1))
         uxxxx23r(i1,i2,i3,kd)=(6.*u(i1,i2,i3,kd)-4.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd))+(u(i1+2,i2,i3,kd)+u(i1-2,i2,i3,kd)) )/(dx(0)**4)
         uyyyy23r(i1,i2,i3,kd)=(6.*u(i1,i2,i3,kd)-4.*(u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd))+(u(i1,i2+2,i3,kd)+u(i1,i2-2,i3,kd)) )/(dx(1)**4)
         uzzzz23r(i1,i2,i3,kd)=(6.*u(i1,i2,i3,kd)-4.*(u(i1,i2,i3+1,kd)+u(i1,i2,i3-1,kd))+(u(i1,i2,i3+2,kd)+u(i1,i2,i3-2,kd)) )/(dx(2)**4)
         uxxyy23r(i1,i2,i3,kd)=( 4.*u(i1,i2,i3,kd)     -2.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd)+u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd))   +   (u(i1+1,i2+1,i3,kd)+u(i1-1,i2+1,i3,kd)+u(i1+1,i2-1,i3,kd)+u(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)
         uxxzz23r(i1,i2,i3,kd)=( 4.*u(i1,i2,i3,kd)     -2.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd)+u(i1,i2,i3+1,kd)+u(i1,i2,i3-1,kd))   +   (u(i1+1,i2,i3+1,kd)+u(i1-1,i2,i3+1,kd)+u(i1+1,i2,i3-1,kd)+u(i1-1,i2,i3-1,kd)) )/(dx(0)**2*dx(2)**2)
         uyyzz23r(i1,i2,i3,kd)=( 4.*u(i1,i2,i3,kd)     -2.*(u(i1,i2+1,i3,kd)  +u(i1,i2-1,i3,kd)+  u(i1,i2  ,i3+1,kd)+u(i1,i2  ,i3-1,kd))   +   (u(i1,i2+1,i3+1,kd)+u(i1,i2-1,i3+1,kd)+u(i1,i2+1,i3-1,kd)+u(i1,i2-1,i3-1,kd)) )/(dx(1)**2*dx(2)**2)
         ! 3D laplacian squared = u.xxxx + u.yyyy + u.zzzz + 2 (u.xxyy + u.xxzz + u.yyzz )
         uLapSq23r(i1,i2,i3,kd)= ( 6.*u(i1,i2,i3,kd)   - 4.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd))    +(u(i1+2,i2,i3,kd)+u(i1-2,i2,i3,kd)) )/(dx(0)**4) +( 6.*u(i1,i2,i3,kd)    -4.*(u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd))    +(u(i1,i2+2,i3,kd)+u(i1,i2-2,i3,kd)) )/(dx(1)**4)  +( 6.*u(i1,i2,i3,kd)    -4.*(u(i1,i2,i3+1,kd)+u(i1,i2,i3-1,kd))    +(u(i1,i2,i3+2,kd)+u(i1,i2,i3-2,kd)) )/(dx(2)**4)  +( 8.*u(i1,i2,i3,kd)     -4.*(u(i1+1,i2,i3,kd)  +u(i1-1,i2,i3,kd)  +u(i1  ,i2+1,i3,kd)+u(i1  ,i2-1,i3,kd))   +2.*(u(i1+1,i2+1,i3,kd)+u(i1-1,i2+1,i3,kd)+u(i1+1,i2-1,i3,kd)+u(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)+( 8.*u(i1,i2,i3,kd)     -4.*(u(i1+1,i2,i3,kd)  +u(i1-1,i2,i3,kd)  +u(i1  ,i2,i3+1,kd)+u(i1  ,i2,i3-1,kd))   +2.*(u(i1+1,i2,i3+1,kd)+u(i1-1,i2,i3+1,kd)+u(i1+1,i2,i3-1,kd)+u(i1-1,i2,i3-1,kd)) )/(dx(0)**2*dx(2)**2)+( 8.*u(i1,i2,i3,kd)     -4.*(u(i1,i2+1,i3,kd)  +u(i1,i2-1,i3,kd)  +u(i1,i2  ,i3+1,kd)+u(i1,i2  ,i3-1,kd))   +2.*(u(i1,i2+1,i3+1,kd)+u(i1,i2-1,i3+1,kd)+u(i1,i2+1,i3-1,kd)+u(i1,i2-1,i3-1,kd)) )/(dx(1)**2*dx(2)**2)
         d14(kd) = 1./(12.*dr(kd))
         d24(kd) = 1./(12.*dr(kd)**2)
         ur4(i1,i2,i3,kd)=(8.*(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))-(u(i1+2,i2,i3,kd)-u(i1-2,i2,i3,kd)))*d14(0)
         us4(i1,i2,i3,kd)=(8.*(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))-(u(i1,i2+2,i3,kd)-u(i1,i2-2,i3,kd)))*d14(1)
         ut4(i1,i2,i3,kd)=(8.*(u(i1,i2,i3+1,kd)-u(i1,i2,i3-1,kd))-(u(i1,i2,i3+2,kd)-u(i1,i2,i3-2,kd)))*d14(2)
         urr4(i1,i2,i3,kd)=(-30.*u(i1,i2,i3,kd)+16.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd))-(u(i1+2,i2,i3,kd)+u(i1-2,i2,i3,kd)) )*d24(0)
         uss4(i1,i2,i3,kd)=(-30.*u(i1,i2,i3,kd)+16.*(u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd))-(u(i1,i2+2,i3,kd)+u(i1,i2-2,i3,kd)) )*d24(1)
         utt4(i1,i2,i3,kd)=(-30.*u(i1,i2,i3,kd)+16.*(u(i1,i2,i3+1,kd)+u(i1,i2,i3-1,kd))-(u(i1,i2,i3+2,kd)+u(i1,i2,i3-2,kd)) )*d24(2)
         urs4(i1,i2,i3,kd)=(8.*(ur4(i1,i2+1,i3,kd)-ur4(i1,i2-1,i3,kd))-(ur4(i1,i2+2,i3,kd)-ur4(i1,i2-2,i3,kd)))*d14(1)
         urt4(i1,i2,i3,kd)=(8.*(ur4(i1,i2,i3+1,kd)-ur4(i1,i2,i3-1,kd))-(ur4(i1,i2,i3+2,kd)-ur4(i1,i2,i3-2,kd)))*d14(2)
         ust4(i1,i2,i3,kd)=(8.*(us4(i1,i2,i3+1,kd)-us4(i1,i2,i3-1,kd))-(us4(i1,i2,i3+2,kd)-us4(i1,i2,i3-2,kd)))*d14(2)
         rxr4(i1,i2,i3)=(8.*(rx(i1+1,i2,i3)-rx(i1-1,i2,i3))-(rx(i1+2,i2,i3)-rx(i1-2,i2,i3)))*d14(0)
         rxs4(i1,i2,i3)=(8.*(rx(i1,i2+1,i3)-rx(i1,i2-1,i3))-(rx(i1,i2+2,i3)-rx(i1,i2-2,i3)))*d14(1)
         rxt4(i1,i2,i3)=(8.*(rx(i1,i2,i3+1)-rx(i1,i2,i3-1))-(rx(i1,i2,i3+2)-rx(i1,i2,i3-2)))*d14(2)
         ryr4(i1,i2,i3)=(8.*(ry(i1+1,i2,i3)-ry(i1-1,i2,i3))-(ry(i1+2,i2,i3)-ry(i1-2,i2,i3)))*d14(0)
         rys4(i1,i2,i3)=(8.*(ry(i1,i2+1,i3)-ry(i1,i2-1,i3))-(ry(i1,i2+2,i3)-ry(i1,i2-2,i3)))*d14(1)
         ryt4(i1,i2,i3)=(8.*(ry(i1,i2,i3+1)-ry(i1,i2,i3-1))-(ry(i1,i2,i3+2)-ry(i1,i2,i3-2)))*d14(2)
         rzr4(i1,i2,i3)=(8.*(rz(i1+1,i2,i3)-rz(i1-1,i2,i3))-(rz(i1+2,i2,i3)-rz(i1-2,i2,i3)))*d14(0)
         rzs4(i1,i2,i3)=(8.*(rz(i1,i2+1,i3)-rz(i1,i2-1,i3))-(rz(i1,i2+2,i3)-rz(i1,i2-2,i3)))*d14(1)
         rzt4(i1,i2,i3)=(8.*(rz(i1,i2,i3+1)-rz(i1,i2,i3-1))-(rz(i1,i2,i3+2)-rz(i1,i2,i3-2)))*d14(2)
         sxr4(i1,i2,i3)=(8.*(sx(i1+1,i2,i3)-sx(i1-1,i2,i3))-(sx(i1+2,i2,i3)-sx(i1-2,i2,i3)))*d14(0)
         sxs4(i1,i2,i3)=(8.*(sx(i1,i2+1,i3)-sx(i1,i2-1,i3))-(sx(i1,i2+2,i3)-sx(i1,i2-2,i3)))*d14(1)
         sxt4(i1,i2,i3)=(8.*(sx(i1,i2,i3+1)-sx(i1,i2,i3-1))-(sx(i1,i2,i3+2)-sx(i1,i2,i3-2)))*d14(2)
         syr4(i1,i2,i3)=(8.*(sy(i1+1,i2,i3)-sy(i1-1,i2,i3))-(sy(i1+2,i2,i3)-sy(i1-2,i2,i3)))*d14(0)
         sys4(i1,i2,i3)=(8.*(sy(i1,i2+1,i3)-sy(i1,i2-1,i3))-(sy(i1,i2+2,i3)-sy(i1,i2-2,i3)))*d14(1)
         syt4(i1,i2,i3)=(8.*(sy(i1,i2,i3+1)-sy(i1,i2,i3-1))-(sy(i1,i2,i3+2)-sy(i1,i2,i3-2)))*d14(2)
         szr4(i1,i2,i3)=(8.*(sz(i1+1,i2,i3)-sz(i1-1,i2,i3))-(sz(i1+2,i2,i3)-sz(i1-2,i2,i3)))*d14(0)
         szs4(i1,i2,i3)=(8.*(sz(i1,i2+1,i3)-sz(i1,i2-1,i3))-(sz(i1,i2+2,i3)-sz(i1,i2-2,i3)))*d14(1)
         szt4(i1,i2,i3)=(8.*(sz(i1,i2,i3+1)-sz(i1,i2,i3-1))-(sz(i1,i2,i3+2)-sz(i1,i2,i3-2)))*d14(2)
         txr4(i1,i2,i3)=(8.*(tx(i1+1,i2,i3)-tx(i1-1,i2,i3))-(tx(i1+2,i2,i3)-tx(i1-2,i2,i3)))*d14(0)
         txs4(i1,i2,i3)=(8.*(tx(i1,i2+1,i3)-tx(i1,i2-1,i3))-(tx(i1,i2+2,i3)-tx(i1,i2-2,i3)))*d14(1)
         txt4(i1,i2,i3)=(8.*(tx(i1,i2,i3+1)-tx(i1,i2,i3-1))-(tx(i1,i2,i3+2)-tx(i1,i2,i3-2)))*d14(2)
         tyr4(i1,i2,i3)=(8.*(ty(i1+1,i2,i3)-ty(i1-1,i2,i3))-(ty(i1+2,i2,i3)-ty(i1-2,i2,i3)))*d14(0)
         tys4(i1,i2,i3)=(8.*(ty(i1,i2+1,i3)-ty(i1,i2-1,i3))-(ty(i1,i2+2,i3)-ty(i1,i2-2,i3)))*d14(1)
         tyt4(i1,i2,i3)=(8.*(ty(i1,i2,i3+1)-ty(i1,i2,i3-1))-(ty(i1,i2,i3+2)-ty(i1,i2,i3-2)))*d14(2)
         tzr4(i1,i2,i3)=(8.*(tz(i1+1,i2,i3)-tz(i1-1,i2,i3))-(tz(i1+2,i2,i3)-tz(i1-2,i2,i3)))*d14(0)
         tzs4(i1,i2,i3)=(8.*(tz(i1,i2+1,i3)-tz(i1,i2-1,i3))-(tz(i1,i2+2,i3)-tz(i1,i2-2,i3)))*d14(1)
         tzt4(i1,i2,i3)=(8.*(tz(i1,i2,i3+1)-tz(i1,i2,i3-1))-(tz(i1,i2,i3+2)-tz(i1,i2,i3-2)))*d14(2)
         ux41(i1,i2,i3,kd)= rx(i1,i2,i3)*ur4(i1,i2,i3,kd)
         uy41(i1,i2,i3,kd)=0
         uz41(i1,i2,i3,kd)=0
         ux42(i1,i2,i3,kd)= rx(i1,i2,i3)*ur4(i1,i2,i3,kd)+sx(i1,i2,i3)*us4(i1,i2,i3,kd)
         uy42(i1,i2,i3,kd)= ry(i1,i2,i3)*ur4(i1,i2,i3,kd)+sy(i1,i2,i3)*us4(i1,i2,i3,kd)
         uz42(i1,i2,i3,kd)=0
         ux43(i1,i2,i3,kd)=rx(i1,i2,i3)*ur4(i1,i2,i3,kd)+sx(i1,i2,i3)*us4(i1,i2,i3,kd)+tx(i1,i2,i3)*ut4(i1,i2,i3,kd)
         uy43(i1,i2,i3,kd)=ry(i1,i2,i3)*ur4(i1,i2,i3,kd)+sy(i1,i2,i3)*us4(i1,i2,i3,kd)+ty(i1,i2,i3)*ut4(i1,i2,i3,kd)
         uz43(i1,i2,i3,kd)=rz(i1,i2,i3)*ur4(i1,i2,i3,kd)+sz(i1,i2,i3)*us4(i1,i2,i3,kd)+tz(i1,i2,i3)*ut4(i1,i2,i3,kd)
         rxx41(i1,i2,i3)= rx(i1,i2,i3)*rxr4(i1,i2,i3)
         rxx42(i1,i2,i3)= rx(i1,i2,i3)*rxr4(i1,i2,i3)+sx(i1,i2,i3)*rxs4(i1,i2,i3)
         rxy42(i1,i2,i3)= ry(i1,i2,i3)*rxr4(i1,i2,i3)+sy(i1,i2,i3)*rxs4(i1,i2,i3)
         rxx43(i1,i2,i3)=rx(i1,i2,i3)*rxr4(i1,i2,i3)+sx(i1,i2,i3)*rxs4(i1,i2,i3)+tx(i1,i2,i3)*rxt4(i1,i2,i3)
         rxy43(i1,i2,i3)=ry(i1,i2,i3)*rxr4(i1,i2,i3)+sy(i1,i2,i3)*rxs4(i1,i2,i3)+ty(i1,i2,i3)*rxt4(i1,i2,i3)
         rxz43(i1,i2,i3)=rz(i1,i2,i3)*rxr4(i1,i2,i3)+sz(i1,i2,i3)*rxs4(i1,i2,i3)+tz(i1,i2,i3)*rxt4(i1,i2,i3)
         ryx42(i1,i2,i3)= rx(i1,i2,i3)*ryr4(i1,i2,i3)+sx(i1,i2,i3)*rys4(i1,i2,i3)
         ryy42(i1,i2,i3)= ry(i1,i2,i3)*ryr4(i1,i2,i3)+sy(i1,i2,i3)*rys4(i1,i2,i3)
         ryx43(i1,i2,i3)=rx(i1,i2,i3)*ryr4(i1,i2,i3)+sx(i1,i2,i3)*rys4(i1,i2,i3)+tx(i1,i2,i3)*ryt4(i1,i2,i3)
         ryy43(i1,i2,i3)=ry(i1,i2,i3)*ryr4(i1,i2,i3)+sy(i1,i2,i3)*rys4(i1,i2,i3)+ty(i1,i2,i3)*ryt4(i1,i2,i3)
         ryz43(i1,i2,i3)=rz(i1,i2,i3)*ryr4(i1,i2,i3)+sz(i1,i2,i3)*rys4(i1,i2,i3)+tz(i1,i2,i3)*ryt4(i1,i2,i3)
         rzx42(i1,i2,i3)= rx(i1,i2,i3)*rzr4(i1,i2,i3)+sx(i1,i2,i3)*rzs4(i1,i2,i3)
         rzy42(i1,i2,i3)= ry(i1,i2,i3)*rzr4(i1,i2,i3)+sy(i1,i2,i3)*rzs4(i1,i2,i3)
         rzx43(i1,i2,i3)=rx(i1,i2,i3)*rzr4(i1,i2,i3)+sx(i1,i2,i3)*rzs4(i1,i2,i3)+tx(i1,i2,i3)*rzt4(i1,i2,i3)
         rzy43(i1,i2,i3)=ry(i1,i2,i3)*rzr4(i1,i2,i3)+sy(i1,i2,i3)*rzs4(i1,i2,i3)+ty(i1,i2,i3)*rzt4(i1,i2,i3)
         rzz43(i1,i2,i3)=rz(i1,i2,i3)*rzr4(i1,i2,i3)+sz(i1,i2,i3)*rzs4(i1,i2,i3)+tz(i1,i2,i3)*rzt4(i1,i2,i3)
         sxx42(i1,i2,i3)= rx(i1,i2,i3)*sxr4(i1,i2,i3)+sx(i1,i2,i3)*sxs4(i1,i2,i3)
         sxy42(i1,i2,i3)= ry(i1,i2,i3)*sxr4(i1,i2,i3)+sy(i1,i2,i3)*sxs4(i1,i2,i3)
         sxx43(i1,i2,i3)=rx(i1,i2,i3)*sxr4(i1,i2,i3)+sx(i1,i2,i3)*sxs4(i1,i2,i3)+tx(i1,i2,i3)*sxt4(i1,i2,i3)
         sxy43(i1,i2,i3)=ry(i1,i2,i3)*sxr4(i1,i2,i3)+sy(i1,i2,i3)*sxs4(i1,i2,i3)+ty(i1,i2,i3)*sxt4(i1,i2,i3)
         sxz43(i1,i2,i3)=rz(i1,i2,i3)*sxr4(i1,i2,i3)+sz(i1,i2,i3)*sxs4(i1,i2,i3)+tz(i1,i2,i3)*sxt4(i1,i2,i3)
         syx42(i1,i2,i3)= rx(i1,i2,i3)*syr4(i1,i2,i3)+sx(i1,i2,i3)*sys4(i1,i2,i3)
         syy42(i1,i2,i3)= ry(i1,i2,i3)*syr4(i1,i2,i3)+sy(i1,i2,i3)*sys4(i1,i2,i3)
         syx43(i1,i2,i3)=rx(i1,i2,i3)*syr4(i1,i2,i3)+sx(i1,i2,i3)*sys4(i1,i2,i3)+tx(i1,i2,i3)*syt4(i1,i2,i3)
         syy43(i1,i2,i3)=ry(i1,i2,i3)*syr4(i1,i2,i3)+sy(i1,i2,i3)*sys4(i1,i2,i3)+ty(i1,i2,i3)*syt4(i1,i2,i3)
         syz43(i1,i2,i3)=rz(i1,i2,i3)*syr4(i1,i2,i3)+sz(i1,i2,i3)*sys4(i1,i2,i3)+tz(i1,i2,i3)*syt4(i1,i2,i3)
         szx42(i1,i2,i3)= rx(i1,i2,i3)*szr4(i1,i2,i3)+sx(i1,i2,i3)*szs4(i1,i2,i3)
         szy42(i1,i2,i3)= ry(i1,i2,i3)*szr4(i1,i2,i3)+sy(i1,i2,i3)*szs4(i1,i2,i3)
         szx43(i1,i2,i3)=rx(i1,i2,i3)*szr4(i1,i2,i3)+sx(i1,i2,i3)*szs4(i1,i2,i3)+tx(i1,i2,i3)*szt4(i1,i2,i3)
         szy43(i1,i2,i3)=ry(i1,i2,i3)*szr4(i1,i2,i3)+sy(i1,i2,i3)*szs4(i1,i2,i3)+ty(i1,i2,i3)*szt4(i1,i2,i3)
         szz43(i1,i2,i3)=rz(i1,i2,i3)*szr4(i1,i2,i3)+sz(i1,i2,i3)*szs4(i1,i2,i3)+tz(i1,i2,i3)*szt4(i1,i2,i3)
         txx42(i1,i2,i3)= rx(i1,i2,i3)*txr4(i1,i2,i3)+sx(i1,i2,i3)*txs4(i1,i2,i3)
         txy42(i1,i2,i3)= ry(i1,i2,i3)*txr4(i1,i2,i3)+sy(i1,i2,i3)*txs4(i1,i2,i3)
         txx43(i1,i2,i3)=rx(i1,i2,i3)*txr4(i1,i2,i3)+sx(i1,i2,i3)*txs4(i1,i2,i3)+tx(i1,i2,i3)*txt4(i1,i2,i3)
         txy43(i1,i2,i3)=ry(i1,i2,i3)*txr4(i1,i2,i3)+sy(i1,i2,i3)*txs4(i1,i2,i3)+ty(i1,i2,i3)*txt4(i1,i2,i3)
         txz43(i1,i2,i3)=rz(i1,i2,i3)*txr4(i1,i2,i3)+sz(i1,i2,i3)*txs4(i1,i2,i3)+tz(i1,i2,i3)*txt4(i1,i2,i3)
         tyx42(i1,i2,i3)= rx(i1,i2,i3)*tyr4(i1,i2,i3)+sx(i1,i2,i3)*tys4(i1,i2,i3)
         tyy42(i1,i2,i3)= ry(i1,i2,i3)*tyr4(i1,i2,i3)+sy(i1,i2,i3)*tys4(i1,i2,i3)
         tyx43(i1,i2,i3)=rx(i1,i2,i3)*tyr4(i1,i2,i3)+sx(i1,i2,i3)*tys4(i1,i2,i3)+tx(i1,i2,i3)*tyt4(i1,i2,i3)
         tyy43(i1,i2,i3)=ry(i1,i2,i3)*tyr4(i1,i2,i3)+sy(i1,i2,i3)*tys4(i1,i2,i3)+ty(i1,i2,i3)*tyt4(i1,i2,i3)
         tyz43(i1,i2,i3)=rz(i1,i2,i3)*tyr4(i1,i2,i3)+sz(i1,i2,i3)*tys4(i1,i2,i3)+tz(i1,i2,i3)*tyt4(i1,i2,i3)
         tzx42(i1,i2,i3)= rx(i1,i2,i3)*tzr4(i1,i2,i3)+sx(i1,i2,i3)*tzs4(i1,i2,i3)
         tzy42(i1,i2,i3)= ry(i1,i2,i3)*tzr4(i1,i2,i3)+sy(i1,i2,i3)*tzs4(i1,i2,i3)
         tzx43(i1,i2,i3)=rx(i1,i2,i3)*tzr4(i1,i2,i3)+sx(i1,i2,i3)*tzs4(i1,i2,i3)+tx(i1,i2,i3)*tzt4(i1,i2,i3)
         tzy43(i1,i2,i3)=ry(i1,i2,i3)*tzr4(i1,i2,i3)+sy(i1,i2,i3)*tzs4(i1,i2,i3)+ty(i1,i2,i3)*tzt4(i1,i2,i3)
         tzz43(i1,i2,i3)=rz(i1,i2,i3)*tzr4(i1,i2,i3)+sz(i1,i2,i3)*tzs4(i1,i2,i3)+tz(i1,i2,i3)*tzt4(i1,i2,i3)
         uxx41(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*urr4(i1,i2,i3,kd)+(rxx42(i1,i2,i3))*ur4(i1,i2,i3,kd)
         uyy41(i1,i2,i3,kd)=0
         uxy41(i1,i2,i3,kd)=0
         uxz41(i1,i2,i3,kd)=0
         uyz41(i1,i2,i3,kd)=0
         uzz41(i1,i2,i3,kd)=0
         ulaplacian41(i1,i2,i3,kd)=uxx41(i1,i2,i3,kd)
         uxx42(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*urr4(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3))*urs4(i1,i2,i3,kd)+(sx(i1,i2,i3)**2)*uss4(i1,i2,i3,kd)+(rxx42(i1,i2,i3))*ur4(i1,i2,i3,kd)+(sxx42(i1,i2,i3))*us4(i1,i2,i3,kd)
         uyy42(i1,i2,i3,kd)=(ry(i1,i2,i3)**2)*urr4(i1,i2,i3,kd)+2.*(ry(i1,i2,i3)*sy(i1,i2,i3))*urs4(i1,i2,i3,kd)+(sy(i1,i2,i3)**2)*uss4(i1,i2,i3,kd)+(ryy42(i1,i2,i3))*ur4(i1,i2,i3,kd)+(syy42(i1,i2,i3))*us4(i1,i2,i3,kd)
         uxy42(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*urr4(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*urs4(i1,i2,i3,kd)+sx(i1,i2,i3)*sy(i1,i2,i3)*uss4(i1,i2,i3,kd)+rxy42(i1,i2,i3)*ur4(i1,i2,i3,kd)+sxy42(i1,i2,i3)*us4(i1,i2,i3,kd)
         uxz42(i1,i2,i3,kd)=0
         uyz42(i1,i2,i3,kd)=0
         uzz42(i1,i2,i3,kd)=0
         ulaplacian42(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2)*urr4(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*sy(i1,i2,i3))*urs4(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**2)*uss4(i1,i2,i3,kd)+(rxx42(i1,i2,i3)+ryy42(i1,i2,i3))*ur4(i1,i2,i3,kd)+(sxx42(i1,i2,i3)+syy42(i1,i2,i3))*us4(i1,i2,i3,kd)
         uxx43(i1,i2,i3,kd)=rx(i1,i2,i3)**2*urr4(i1,i2,i3,kd)+sx(i1,i2,i3)**2*uss4(i1,i2,i3,kd)+tx(i1,i2,i3)**2*utt4(i1,i2,i3,kd)+2.*rx(i1,i2,i3)*sx(i1,i2,i3)*urs4(i1,i2,i3,kd)+2.*rx(i1,i2,i3)*tx(i1,i2,i3)*urt4(i1,i2,i3,kd)+2.*sx(i1,i2,i3)*tx(i1,i2,i3)*ust4(i1,i2,i3,kd)+rxx43(i1,i2,i3)*ur4(i1,i2,i3,kd)+sxx43(i1,i2,i3)*us4(i1,i2,i3,kd)+txx43(i1,i2,i3)*ut4(i1,i2,i3,kd)
         uyy43(i1,i2,i3,kd)=ry(i1,i2,i3)**2*urr4(i1,i2,i3,kd)+sy(i1,i2,i3)**2*uss4(i1,i2,i3,kd)+ty(i1,i2,i3)**2*utt4(i1,i2,i3,kd)+2.*ry(i1,i2,i3)*sy(i1,i2,i3)*urs4(i1,i2,i3,kd)+2.*ry(i1,i2,i3)*ty(i1,i2,i3)*urt4(i1,i2,i3,kd)+2.*sy(i1,i2,i3)*ty(i1,i2,i3)*ust4(i1,i2,i3,kd)+ryy43(i1,i2,i3)*ur4(i1,i2,i3,kd)+syy43(i1,i2,i3)*us4(i1,i2,i3,kd)+tyy43(i1,i2,i3)*ut4(i1,i2,i3,kd)
         uzz43(i1,i2,i3,kd)=rz(i1,i2,i3)**2*urr4(i1,i2,i3,kd)+sz(i1,i2,i3)**2*uss4(i1,i2,i3,kd)+tz(i1,i2,i3)**2*utt4(i1,i2,i3,kd)+2.*rz(i1,i2,i3)*sz(i1,i2,i3)*urs4(i1,i2,i3,kd)+2.*rz(i1,i2,i3)*tz(i1,i2,i3)*urt4(i1,i2,i3,kd)+2.*sz(i1,i2,i3)*tz(i1,i2,i3)*ust4(i1,i2,i3,kd)+rzz43(i1,i2,i3)*ur4(i1,i2,i3,kd)+szz43(i1,i2,i3)*us4(i1,i2,i3,kd)+tzz43(i1,i2,i3)*ut4(i1,i2,i3,kd)
         uxy43(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*urr4(i1,i2,i3,kd)+sx(i1,i2,i3)*sy(i1,i2,i3)*uss4(i1,i2,i3,kd)+tx(i1,i2,i3)*ty(i1,i2,i3)*utt4(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*urs4(i1,i2,i3,kd)+(rx(i1,i2,i3)*ty(i1,i2,i3)+ry(i1,i2,i3)*tx(i1,i2,i3))*urt4(i1,i2,i3,kd)+(sx(i1,i2,i3)*ty(i1,i2,i3)+sy(i1,i2,i3)*tx(i1,i2,i3))*ust4(i1,i2,i3,kd)+rxy43(i1,i2,i3)*ur4(i1,i2,i3,kd)+sxy43(i1,i2,i3)*us4(i1,i2,i3,kd)+txy43(i1,i2,i3)*ut4(i1,i2,i3,kd)
         uxz43(i1,i2,i3,kd)=rx(i1,i2,i3)*rz(i1,i2,i3)*urr4(i1,i2,i3,kd)+sx(i1,i2,i3)*sz(i1,i2,i3)*uss4(i1,i2,i3,kd)+tx(i1,i2,i3)*tz(i1,i2,i3)*utt4(i1,i2,i3,kd)+(rx(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,i2,i3)*sx(i1,i2,i3))*urs4(i1,i2,i3,kd)+(rx(i1,i2,i3)*tz(i1,i2,i3)+rz(i1,i2,i3)*tx(i1,i2,i3))*urt4(i1,i2,i3,kd)+(sx(i1,i2,i3)*tz(i1,i2,i3)+sz(i1,i2,i3)*tx(i1,i2,i3))*ust4(i1,i2,i3,kd)+rxz43(i1,i2,i3)*ur4(i1,i2,i3,kd)+sxz43(i1,i2,i3)*us4(i1,i2,i3,kd)+txz43(i1,i2,i3)*ut4(i1,i2,i3,kd)
         uyz43(i1,i2,i3,kd)=ry(i1,i2,i3)*rz(i1,i2,i3)*urr4(i1,i2,i3,kd)+sy(i1,i2,i3)*sz(i1,i2,i3)*uss4(i1,i2,i3,kd)+ty(i1,i2,i3)*tz(i1,i2,i3)*utt4(i1,i2,i3,kd)+(ry(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,i2,i3)*sy(i1,i2,i3))*urs4(i1,i2,i3,kd)+(ry(i1,i2,i3)*tz(i1,i2,i3)+rz(i1,i2,i3)*ty(i1,i2,i3))*urt4(i1,i2,i3,kd)+(sy(i1,i2,i3)*tz(i1,i2,i3)+sz(i1,i2,i3)*ty(i1,i2,i3))*ust4(i1,i2,i3,kd)+ryz43(i1,i2,i3)*ur4(i1,i2,i3,kd)+syz43(i1,i2,i3)*us4(i1,i2,i3,kd)+tyz43(i1,i2,i3)*ut4(i1,i2,i3,kd)
         ulaplacian43(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2+rz(i1,i2,i3)**2)*urr4(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**2+sz(i1,i2,i3)**2)*uss4(i1,i2,i3,kd)+(tx(i1,i2,i3)**2+ty(i1,i2,i3)**2+tz(i1,i2,i3)**2)*utt4(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*sy(i1,i2,i3)+rz(i1,i2,i3)*sz(i1,i2,i3))*urs4(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*tx(i1,i2,i3)+ ry(i1,i2,i3)*ty(i1,i2,i3)+rz(i1,i2,i3)*tz(i1,i2,i3))*urt4(i1,i2,i3,kd)+2.*(sx(i1,i2,i3)*tx(i1,i2,i3)+ sy(i1,i2,i3)*ty(i1,i2,i3)+sz(i1,i2,i3)*tz(i1,i2,i3))*ust4(i1,i2,i3,kd)+(rxx43(i1,i2,i3)+ryy43(i1,i2,i3)+rzz43(i1,i2,i3))*ur4(i1,i2,i3,kd)+(sxx43(i1,i2,i3)+syy43(i1,i2,i3)+szz43(i1,i2,i3))*us4(i1,i2,i3,kd)+(txx43(i1,i2,i3)+tyy43(i1,i2,i3)+tzz43(i1,i2,i3))*ut4(i1,i2,i3,kd)
         !============================================================================================
         ! Define derivatives for a rectangular grid
         !
         !============================================================================================
         h41(kd) = 1./(12.*dx(kd))
         h42(kd) = 1./(12.*dx(kd)**2)
         ux43r(i1,i2,i3,kd)=(8.*(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))-(u(i1+2,i2,i3,kd)-u(i1-2,i2,i3,kd)))*h41(0)
         uy43r(i1,i2,i3,kd)=(8.*(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))-(u(i1,i2+2,i3,kd)-u(i1,i2-2,i3,kd)))*h41(1)
         uz43r(i1,i2,i3,kd)=(8.*(u(i1,i2,i3+1,kd)-u(i1,i2,i3-1,kd))-(u(i1,i2,i3+2,kd)-u(i1,i2,i3-2,kd)))*h41(2)
         uxx43r(i1,i2,i3,kd)=( -30.*u(i1,i2,i3,kd)+16.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd))-(u(i1+2,i2,i3,kd)+u(i1-2,i2,i3,kd)) )*h42(0) 
         uyy43r(i1,i2,i3,kd)=( -30.*u(i1,i2,i3,kd)+16.*(u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd))-(u(i1,i2+2,i3,kd)+u(i1,i2-2,i3,kd)) )*h42(1) 
         uzz43r(i1,i2,i3,kd)=( -30.*u(i1,i2,i3,kd)+16.*(u(i1,i2,i3+1,kd)+u(i1,i2,i3-1,kd))-(u(i1,i2,i3+2,kd)+u(i1,i2,i3-2,kd)) )*h42(2)
         uxy43r(i1,i2,i3,kd)=( (u(i1+2,i2+2,i3,kd)-u(i1-2,i2+2,i3,kd)- u(i1+2,i2-2,i3,kd)+u(i1-2,i2-2,i3,kd)) +8.*(u(i1-1,i2+2,i3,kd)-u(i1-1,i2-2,i3,kd)-u(i1+1,i2+2,i3,kd)+u(i1+1,i2-2,i3,kd) +u(i1+2,i2-1,i3,kd)-u(i1-2,i2-1,i3,kd)-u(i1+2,i2+1,i3,kd)+u(i1-2,i2+1,i3,kd))+64.*(u(i1+1,i2+1,i3,kd)-u(i1-1,i2+1,i3,kd)- u(i1+1,i2-1,i3,kd)+u(i1-1,i2-1,i3,kd)))*(h41(0)*h41(1))
         uxz43r(i1,i2,i3,kd)=( (u(i1+2,i2,i3+2,kd)-u(i1-2,i2,i3+2,kd)-u(i1+2,i2,i3-2,kd)+u(i1-2,i2,i3-2,kd)) +8.*(u(i1-1,i2,i3+2,kd)-u(i1-1,i2,i3-2,kd)-u(i1+1,i2,i3+2,kd)+u(i1+1,i2,i3-2,kd) +u(i1+2,i2,i3-1,kd)-u(i1-2,i2,i3-1,kd)- u(i1+2,i2,i3+1,kd)+u(i1-2,i2,i3+1,kd)) +64.*(u(i1+1,i2,i3+1,kd)-u(i1-1,i2,i3+1,kd)-u(i1+1,i2,i3-1,kd)+u(i1-1,i2,i3-1,kd)) )*(h41(0)*h41(2))
         uyz43r(i1,i2,i3,kd)=( (u(i1,i2+2,i3+2,kd)-u(i1,i2-2,i3+2,kd)-u(i1,i2+2,i3-2,kd)+u(i1,i2-2,i3-2,kd)) +8.*(u(i1,i2-1,i3+2,kd)-u(i1,i2-1,i3-2,kd)-u(i1,i2+1,i3+2,kd)+u(i1,i2+1,i3-2,kd) +u(i1,i2+2,i3-1,kd)-u(i1,i2-2,i3-1,kd)-u(i1,i2+2,i3+1,kd)+u(i1,i2-2,i3+1,kd)) +64.*(u(i1,i2+1,i3+1,kd)-u(i1,i2-1,i3+1,kd)-u(i1,i2+1,i3-1,kd)+u(i1,i2-1,i3-1,kd)) )*(h41(1)*h41(2))
         ux41r(i1,i2,i3,kd)= ux43r(i1,i2,i3,kd)
         uy41r(i1,i2,i3,kd)= uy43r(i1,i2,i3,kd)
         uz41r(i1,i2,i3,kd)= uz43r(i1,i2,i3,kd)
         uxx41r(i1,i2,i3,kd)= uxx43r(i1,i2,i3,kd)
         uyy41r(i1,i2,i3,kd)= uyy43r(i1,i2,i3,kd)
         uzz41r(i1,i2,i3,kd)= uzz43r(i1,i2,i3,kd)
         uxy41r(i1,i2,i3,kd)= uxy43r(i1,i2,i3,kd)
         uxz41r(i1,i2,i3,kd)= uxz43r(i1,i2,i3,kd)
         uyz41r(i1,i2,i3,kd)= uyz43r(i1,i2,i3,kd)
         ulaplacian41r(i1,i2,i3,kd)=uxx43r(i1,i2,i3,kd)
         ux42r(i1,i2,i3,kd)= ux43r(i1,i2,i3,kd)
         uy42r(i1,i2,i3,kd)= uy43r(i1,i2,i3,kd)
         uz42r(i1,i2,i3,kd)= uz43r(i1,i2,i3,kd)
         uxx42r(i1,i2,i3,kd)= uxx43r(i1,i2,i3,kd)
         uyy42r(i1,i2,i3,kd)= uyy43r(i1,i2,i3,kd)
         uzz42r(i1,i2,i3,kd)= uzz43r(i1,i2,i3,kd)
         uxy42r(i1,i2,i3,kd)= uxy43r(i1,i2,i3,kd)
         uxz42r(i1,i2,i3,kd)= uxz43r(i1,i2,i3,kd)
         uyz42r(i1,i2,i3,kd)= uyz43r(i1,i2,i3,kd)
         ulaplacian42r(i1,i2,i3,kd)=uxx43r(i1,i2,i3,kd)+uyy43r(i1,i2,i3,kd)
         ulaplacian43r(i1,i2,i3,kd)=uxx43r(i1,i2,i3,kd)+uyy43r(i1,i2,i3,kd)+uzz43r(i1,i2,i3,kd)
         unr2(i1,i2,i3,kd)=(un(i1+1,i2,i3,kd)-un(i1-1,i2,i3,kd))*d12(0)
         uns2(i1,i2,i3,kd)=(un(i1,i2+1,i3,kd)-un(i1,i2-1,i3,kd))*d12(1)
         unt2(i1,i2,i3,kd)=(un(i1,i2,i3+1,kd)-un(i1,i2,i3-1,kd))*d12(2)
         unrr2(i1,i2,i3,kd)=(-2.*un(i1,i2,i3,kd)+(un(i1+1,i2,i3,kd)+un(i1-1,i2,i3,kd)) )*d22(0)
         unss2(i1,i2,i3,kd)=(-2.*un(i1,i2,i3,kd)+(un(i1,i2+1,i3,kd)+un(i1,i2-1,i3,kd)) )*d22(1)
         unrs2(i1,i2,i3,kd)=(unr2(i1,i2+1,i3,kd)-unr2(i1,i2-1,i3,kd))*d12(1)
         untt2(i1,i2,i3,kd)=(-2.*un(i1,i2,i3,kd)+(un(i1,i2,i3+1,kd)+un(i1,i2,i3-1,kd)) )*d22(2)
         unrt2(i1,i2,i3,kd)=(unr2(i1,i2,i3+1,kd)-unr2(i1,i2,i3-1,kd))*d12(2)
         unst2(i1,i2,i3,kd)=(uns2(i1,i2,i3+1,kd)-uns2(i1,i2,i3-1,kd))*d12(2)
         unrrr2(i1,i2,i3,kd)=(-2.*(un(i1+1,i2,i3,kd)-un(i1-1,i2,i3,kd))+(un(i1+2,i2,i3,kd)-un(i1-2,i2,i3,kd)) )*d22(0)*d12(0)
         unsss2(i1,i2,i3,kd)=(-2.*(un(i1,i2+1,i3,kd)-un(i1,i2-1,i3,kd))+(un(i1,i2+2,i3,kd)-un(i1,i2-2,i3,kd)) )*d22(1)*d12(1)
         unttt2(i1,i2,i3,kd)=(-2.*(un(i1,i2,i3+1,kd)-un(i1,i2,i3-1,kd))+(un(i1,i2,i3+2,kd)-un(i1,i2,i3-2,kd)) )*d22(2)*d12(2)
         unx21(i1,i2,i3,kd)= rx(i1,i2,i3)*unr2(i1,i2,i3,kd)
         uny21(i1,i2,i3,kd)=0
         unz21(i1,i2,i3,kd)=0
         unx22(i1,i2,i3,kd)= rx(i1,i2,i3)*unr2(i1,i2,i3,kd)+sx(i1,i2,i3)*uns2(i1,i2,i3,kd)
         uny22(i1,i2,i3,kd)= ry(i1,i2,i3)*unr2(i1,i2,i3,kd)+sy(i1,i2,i3)*uns2(i1,i2,i3,kd)
         unz22(i1,i2,i3,kd)=0
         unx23(i1,i2,i3,kd)=rx(i1,i2,i3)*unr2(i1,i2,i3,kd)+sx(i1,i2,i3)*uns2(i1,i2,i3,kd)+tx(i1,i2,i3)*unt2(i1,i2,i3,kd)
         uny23(i1,i2,i3,kd)=ry(i1,i2,i3)*unr2(i1,i2,i3,kd)+sy(i1,i2,i3)*uns2(i1,i2,i3,kd)+ty(i1,i2,i3)*unt2(i1,i2,i3,kd)
         unz23(i1,i2,i3,kd)=rz(i1,i2,i3)*unr2(i1,i2,i3,kd)+sz(i1,i2,i3)*uns2(i1,i2,i3,kd)+tz(i1,i2,i3)*unt2(i1,i2,i3,kd)
         unxx21(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*unrr2(i1,i2,i3,kd)+(rxx22(i1,i2,i3))*unr2(i1,i2,i3,kd)
         unyy21(i1,i2,i3,kd)=0
         unxy21(i1,i2,i3,kd)=0
         unxz21(i1,i2,i3,kd)=0
         unyz21(i1,i2,i3,kd)=0
         unzz21(i1,i2,i3,kd)=0
         unlaplacian21(i1,i2,i3,kd)=unxx21(i1,i2,i3,kd)
         unxx22(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*unrr2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3))*unrs2(i1,i2,i3,kd)+(sx(i1,i2,i3)**2)*unss2(i1,i2,i3,kd)+(rxx22(i1,i2,i3))*unr2(i1,i2,i3,kd)+(sxx22(i1,i2,i3))*uns2(i1,i2,i3,kd)
         unyy22(i1,i2,i3,kd)=(ry(i1,i2,i3)**2)*unrr2(i1,i2,i3,kd)+2.*(ry(i1,i2,i3)*sy(i1,i2,i3))*unrs2(i1,i2,i3,kd)+(sy(i1,i2,i3)**2)*unss2(i1,i2,i3,kd)+(ryy22(i1,i2,i3))*unr2(i1,i2,i3,kd)+(syy22(i1,i2,i3))*uns2(i1,i2,i3,kd)
         unxy22(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*unrr2(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*unrs2(i1,i2,i3,kd)+sx(i1,i2,i3)*sy(i1,i2,i3)*unss2(i1,i2,i3,kd)+rxy22(i1,i2,i3)*unr2(i1,i2,i3,kd)+sxy22(i1,i2,i3)*uns2(i1,i2,i3,kd)
         unxz22(i1,i2,i3,kd)=0
         unyz22(i1,i2,i3,kd)=0
         unzz22(i1,i2,i3,kd)=0
         unlaplacian22(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2)*unrr2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*sy(i1,i2,i3))*unrs2(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**2)*unss2(i1,i2,i3,kd)+(rxx22(i1,i2,i3)+ryy22(i1,i2,i3))*unr2(i1,i2,i3,kd)+(sxx22(i1,i2,i3)+syy22(i1,i2,i3))*uns2(i1,i2,i3,kd)
         unxx23(i1,i2,i3,kd)=rx(i1,i2,i3)**2*unrr2(i1,i2,i3,kd)+sx(i1,i2,i3)**2*unss2(i1,i2,i3,kd)+tx(i1,i2,i3)**2*untt2(i1,i2,i3,kd)+2.*rx(i1,i2,i3)*sx(i1,i2,i3)*unrs2(i1,i2,i3,kd)+2.*rx(i1,i2,i3)*tx(i1,i2,i3)*unrt2(i1,i2,i3,kd)+2.*sx(i1,i2,i3)*tx(i1,i2,i3)*unst2(i1,i2,i3,kd)+rxx23(i1,i2,i3)*unr2(i1,i2,i3,kd)+sxx23(i1,i2,i3)*uns2(i1,i2,i3,kd)+txx23(i1,i2,i3)*unt2(i1,i2,i3,kd)
         unyy23(i1,i2,i3,kd)=ry(i1,i2,i3)**2*unrr2(i1,i2,i3,kd)+sy(i1,i2,i3)**2*unss2(i1,i2,i3,kd)+ty(i1,i2,i3)**2*untt2(i1,i2,i3,kd)+2.*ry(i1,i2,i3)*sy(i1,i2,i3)*unrs2(i1,i2,i3,kd)+2.*ry(i1,i2,i3)*ty(i1,i2,i3)*unrt2(i1,i2,i3,kd)+2.*sy(i1,i2,i3)*ty(i1,i2,i3)*unst2(i1,i2,i3,kd)+ryy23(i1,i2,i3)*unr2(i1,i2,i3,kd)+syy23(i1,i2,i3)*uns2(i1,i2,i3,kd)+tyy23(i1,i2,i3)*unt2(i1,i2,i3,kd)
         unzz23(i1,i2,i3,kd)=rz(i1,i2,i3)**2*unrr2(i1,i2,i3,kd)+sz(i1,i2,i3)**2*unss2(i1,i2,i3,kd)+tz(i1,i2,i3)**2*untt2(i1,i2,i3,kd)+2.*rz(i1,i2,i3)*sz(i1,i2,i3)*unrs2(i1,i2,i3,kd)+2.*rz(i1,i2,i3)*tz(i1,i2,i3)*unrt2(i1,i2,i3,kd)+2.*sz(i1,i2,i3)*tz(i1,i2,i3)*unst2(i1,i2,i3,kd)+rzz23(i1,i2,i3)*unr2(i1,i2,i3,kd)+szz23(i1,i2,i3)*uns2(i1,i2,i3,kd)+tzz23(i1,i2,i3)*unt2(i1,i2,i3,kd)
         unxy23(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*unrr2(i1,i2,i3,kd)+sx(i1,i2,i3)*sy(i1,i2,i3)*unss2(i1,i2,i3,kd)+tx(i1,i2,i3)*ty(i1,i2,i3)*untt2(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*unrs2(i1,i2,i3,kd)+(rx(i1,i2,i3)*ty(i1,i2,i3)+ry(i1,i2,i3)*tx(i1,i2,i3))*unrt2(i1,i2,i3,kd)+(sx(i1,i2,i3)*ty(i1,i2,i3)+sy(i1,i2,i3)*tx(i1,i2,i3))*unst2(i1,i2,i3,kd)+rxy23(i1,i2,i3)*unr2(i1,i2,i3,kd)+sxy23(i1,i2,i3)*uns2(i1,i2,i3,kd)+txy23(i1,i2,i3)*unt2(i1,i2,i3,kd)
         unxz23(i1,i2,i3,kd)=rx(i1,i2,i3)*rz(i1,i2,i3)*unrr2(i1,i2,i3,kd)+sx(i1,i2,i3)*sz(i1,i2,i3)*unss2(i1,i2,i3,kd)+tx(i1,i2,i3)*tz(i1,i2,i3)*untt2(i1,i2,i3,kd)+(rx(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,i2,i3)*sx(i1,i2,i3))*unrs2(i1,i2,i3,kd)+(rx(i1,i2,i3)*tz(i1,i2,i3)+rz(i1,i2,i3)*tx(i1,i2,i3))*unrt2(i1,i2,i3,kd)+(sx(i1,i2,i3)*tz(i1,i2,i3)+sz(i1,i2,i3)*tx(i1,i2,i3))*unst2(i1,i2,i3,kd)+rxz23(i1,i2,i3)*unr2(i1,i2,i3,kd)+sxz23(i1,i2,i3)*uns2(i1,i2,i3,kd)+txz23(i1,i2,i3)*unt2(i1,i2,i3,kd)
         unyz23(i1,i2,i3,kd)=ry(i1,i2,i3)*rz(i1,i2,i3)*unrr2(i1,i2,i3,kd)+sy(i1,i2,i3)*sz(i1,i2,i3)*unss2(i1,i2,i3,kd)+ty(i1,i2,i3)*tz(i1,i2,i3)*untt2(i1,i2,i3,kd)+(ry(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,i2,i3)*sy(i1,i2,i3))*unrs2(i1,i2,i3,kd)+(ry(i1,i2,i3)*tz(i1,i2,i3)+rz(i1,i2,i3)*ty(i1,i2,i3))*unrt2(i1,i2,i3,kd)+(sy(i1,i2,i3)*tz(i1,i2,i3)+sz(i1,i2,i3)*ty(i1,i2,i3))*unst2(i1,i2,i3,kd)+ryz23(i1,i2,i3)*unr2(i1,i2,i3,kd)+syz23(i1,i2,i3)*uns2(i1,i2,i3,kd)+tyz23(i1,i2,i3)*unt2(i1,i2,i3,kd)
         unlaplacian23(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2+rz(i1,i2,i3)**2)*unrr2(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**2+sz(i1,i2,i3)**2)*unss2(i1,i2,i3,kd)+(tx(i1,i2,i3)**2+ty(i1,i2,i3)**2+tz(i1,i2,i3)**2)*untt2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*sy(i1,i2,i3)+rz(i1,i2,i3)*sz(i1,i2,i3))*unrs2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*tx(i1,i2,i3)+ ry(i1,i2,i3)*ty(i1,i2,i3)+rz(i1,i2,i3)*tz(i1,i2,i3))*unrt2(i1,i2,i3,kd)+2.*(sx(i1,i2,i3)*tx(i1,i2,i3)+ sy(i1,i2,i3)*ty(i1,i2,i3)+sz(i1,i2,i3)*tz(i1,i2,i3))*unst2(i1,i2,i3,kd)+(rxx23(i1,i2,i3)+ryy23(i1,i2,i3)+rzz23(i1,i2,i3))*unr2(i1,i2,i3,kd)+(sxx23(i1,i2,i3)+syy23(i1,i2,i3)+szz23(i1,i2,i3))*uns2(i1,i2,i3,kd)+(txx23(i1,i2,i3)+tyy23(i1,i2,i3)+tzz23(i1,i2,i3))*unt2(i1,i2,i3,kd)
         !============================================================================================
         ! Define derivatives for a rectangular grid
         !
         !============================================================================================
         unx23r(i1,i2,i3,kd)=(un(i1+1,i2,i3,kd)-un(i1-1,i2,i3,kd))*h12(0)
         uny23r(i1,i2,i3,kd)=(un(i1,i2+1,i3,kd)-un(i1,i2-1,i3,kd))*h12(1)
         unz23r(i1,i2,i3,kd)=(un(i1,i2,i3+1,kd)-un(i1,i2,i3-1,kd))*h12(2)
         unxx23r(i1,i2,i3,kd)=(-2.*un(i1,i2,i3,kd)+(un(i1+1,i2,i3,kd)+un(i1-1,i2,i3,kd)) )*h22(0)
         unyy23r(i1,i2,i3,kd)=(-2.*un(i1,i2,i3,kd)+(un(i1,i2+1,i3,kd)+un(i1,i2-1,i3,kd)) )*h22(1)
         unxy23r(i1,i2,i3,kd)=(unx23r(i1,i2+1,i3,kd)-unx23r(i1,i2-1,i3,kd))*h12(1)
         unzz23r(i1,i2,i3,kd)=(-2.*un(i1,i2,i3,kd)+(un(i1,i2,i3+1,kd)+un(i1,i2,i3-1,kd)) )*h22(2)
         unxz23r(i1,i2,i3,kd)=(unx23r(i1,i2,i3+1,kd)-unx23r(i1,i2,i3-1,kd))*h12(2)
         unyz23r(i1,i2,i3,kd)=(uny23r(i1,i2,i3+1,kd)-uny23r(i1,i2,i3-1,kd))*h12(2)
         unx21r(i1,i2,i3,kd)= unx23r(i1,i2,i3,kd)
         uny21r(i1,i2,i3,kd)= uny23r(i1,i2,i3,kd)
         unz21r(i1,i2,i3,kd)= unz23r(i1,i2,i3,kd)
         unxx21r(i1,i2,i3,kd)= unxx23r(i1,i2,i3,kd)
         unyy21r(i1,i2,i3,kd)= unyy23r(i1,i2,i3,kd)
         unzz21r(i1,i2,i3,kd)= unzz23r(i1,i2,i3,kd)
         unxy21r(i1,i2,i3,kd)= unxy23r(i1,i2,i3,kd)
         unxz21r(i1,i2,i3,kd)= unxz23r(i1,i2,i3,kd)
         unyz21r(i1,i2,i3,kd)= unyz23r(i1,i2,i3,kd)
         unlaplacian21r(i1,i2,i3,kd)=unxx23r(i1,i2,i3,kd)
         unx22r(i1,i2,i3,kd)= unx23r(i1,i2,i3,kd)
         uny22r(i1,i2,i3,kd)= uny23r(i1,i2,i3,kd)
         unz22r(i1,i2,i3,kd)= unz23r(i1,i2,i3,kd)
         unxx22r(i1,i2,i3,kd)= unxx23r(i1,i2,i3,kd)
         unyy22r(i1,i2,i3,kd)= unyy23r(i1,i2,i3,kd)
         unzz22r(i1,i2,i3,kd)= unzz23r(i1,i2,i3,kd)
         unxy22r(i1,i2,i3,kd)= unxy23r(i1,i2,i3,kd)
         unxz22r(i1,i2,i3,kd)= unxz23r(i1,i2,i3,kd)
         unyz22r(i1,i2,i3,kd)= unyz23r(i1,i2,i3,kd)
         unlaplacian22r(i1,i2,i3,kd)=unxx23r(i1,i2,i3,kd)+unyy23r(i1,i2,i3,kd)
         unlaplacian23r(i1,i2,i3,kd)=unxx23r(i1,i2,i3,kd)+unyy23r(i1,i2,i3,kd)+unzz23r(i1,i2,i3,kd)
         unxxx22r(i1,i2,i3,kd)=(-2.*(un(i1+1,i2,i3,kd)-un(i1-1,i2,i3,kd))+(un(i1+2,i2,i3,kd)-un(i1-2,i2,i3,kd)) )*h22(0)*h12(0)
         unyyy22r(i1,i2,i3,kd)=(-2.*(un(i1,i2+1,i3,kd)-un(i1,i2-1,i3,kd))+(un(i1,i2+2,i3,kd)-un(i1,i2-2,i3,kd)) )*h22(1)*h12(1)
         unxxy22r(i1,i2,i3,kd)=( unxx22r(i1,i2+1,i3,kd)-unxx22r(i1,i2-1,i3,kd))/(2.*dx(1))
         unxyy22r(i1,i2,i3,kd)=( unyy22r(i1+1,i2,i3,kd)-unyy22r(i1-1,i2,i3,kd))/(2.*dx(0))
         unxxxx22r(i1,i2,i3,kd)=(6.*un(i1,i2,i3,kd)-4.*(un(i1+1,i2,i3,kd)+un(i1-1,i2,i3,kd))+(un(i1+2,i2,i3,kd)+un(i1-2,i2,i3,kd)) )/(dx(0)**4)
         unyyyy22r(i1,i2,i3,kd)=(6.*un(i1,i2,i3,kd)-4.*(un(i1,i2+1,i3,kd)+un(i1,i2-1,i3,kd))+(un(i1,i2+2,i3,kd)+un(i1,i2-2,i3,kd)) )/(dx(1)**4)
         unxxyy22r(i1,i2,i3,kd)=( 4.*un(i1,i2,i3,kd)     -2.*(un(i1+1,i2,i3,kd)+un(i1-1,i2,i3,kd)+un(i1,i2+1,i3,kd)+un(i1,i2-1,i3,kd))   +   (un(i1+1,i2+1,i3,kd)+un(i1-1,i2+1,i3,kd)+un(i1+1,i2-1,i3,kd)+un(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)
         ! 2D laplacian squared = un.xxxx + 2 un.xxyy + un.yyyy
         unLapSq22r(i1,i2,i3,kd)= ( 6.*un(i1,i2,i3,kd)   - 4.*(un(i1+1,i2,i3,kd)+un(i1-1,i2,i3,kd))    +(un(i1+2,i2,i3,kd)+un(i1-2,i2,i3,kd)) )/(dx(0)**4) +( 6.*un(i1,i2,i3,kd)    -4.*(un(i1,i2+1,i3,kd)+un(i1,i2-1,i3,kd))    +(un(i1,i2+2,i3,kd)+un(i1,i2-2,i3,kd)) )/(dx(1)**4)  +( 8.*un(i1,i2,i3,kd)     -4.*(un(i1+1,i2,i3,kd)+un(i1-1,i2,i3,kd)+un(i1,i2+1,i3,kd)+un(i1,i2-1,i3,kd))   +2.*(un(i1+1,i2+1,i3,kd)+un(i1-1,i2+1,i3,kd)+un(i1+1,i2-1,i3,kd)+un(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)
         unxxx23r(i1,i2,i3,kd)=(-2.*(un(i1+1,i2,i3,kd)-un(i1-1,i2,i3,kd))+(un(i1+2,i2,i3,kd)-un(i1-2,i2,i3,kd)) )*h22(0)*h12(0)
         unyyy23r(i1,i2,i3,kd)=(-2.*(un(i1,i2+1,i3,kd)-un(i1,i2-1,i3,kd))+(un(i1,i2+2,i3,kd)-un(i1,i2-2,i3,kd)) )*h22(1)*h12(1)
         unzzz23r(i1,i2,i3,kd)=(-2.*(un(i1,i2,i3+1,kd)-un(i1,i2,i3-1,kd))+(un(i1,i2,i3+2,kd)-un(i1,i2,i3-2,kd)) )*h22(1)*h12(2)
         unxxy23r(i1,i2,i3,kd)=( unxx22r(i1,i2+1,i3,kd)-unxx22r(i1,i2-1,i3,kd))/(2.*dx(1))
         unxyy23r(i1,i2,i3,kd)=( unyy22r(i1+1,i2,i3,kd)-unyy22r(i1-1,i2,i3,kd))/(2.*dx(0))
         unxxz23r(i1,i2,i3,kd)=( unxx22r(i1,i2,i3+1,kd)-unxx22r(i1,i2,i3-1,kd))/(2.*dx(2))
         unyyz23r(i1,i2,i3,kd)=( unyy22r(i1,i2,i3+1,kd)-unyy22r(i1,i2,i3-1,kd))/(2.*dx(2))
         unxzz23r(i1,i2,i3,kd)=( unzz22r(i1+1,i2,i3,kd)-unzz22r(i1-1,i2,i3,kd))/(2.*dx(0))
         unyzz23r(i1,i2,i3,kd)=( unzz22r(i1,i2+1,i3,kd)-unzz22r(i1,i2-1,i3,kd))/(2.*dx(1))
         unxxxx23r(i1,i2,i3,kd)=(6.*un(i1,i2,i3,kd)-4.*(un(i1+1,i2,i3,kd)+un(i1-1,i2,i3,kd))+(un(i1+2,i2,i3,kd)+un(i1-2,i2,i3,kd)) )/(dx(0)**4)
         unyyyy23r(i1,i2,i3,kd)=(6.*un(i1,i2,i3,kd)-4.*(un(i1,i2+1,i3,kd)+un(i1,i2-1,i3,kd))+(un(i1,i2+2,i3,kd)+un(i1,i2-2,i3,kd)) )/(dx(1)**4)
         unzzzz23r(i1,i2,i3,kd)=(6.*un(i1,i2,i3,kd)-4.*(un(i1,i2,i3+1,kd)+un(i1,i2,i3-1,kd))+(un(i1,i2,i3+2,kd)+un(i1,i2,i3-2,kd)) )/(dx(2)**4)
         unxxyy23r(i1,i2,i3,kd)=( 4.*un(i1,i2,i3,kd)     -2.*(un(i1+1,i2,i3,kd)+un(i1-1,i2,i3,kd)+un(i1,i2+1,i3,kd)+un(i1,i2-1,i3,kd))   +   (un(i1+1,i2+1,i3,kd)+un(i1-1,i2+1,i3,kd)+un(i1+1,i2-1,i3,kd)+un(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)
         unxxzz23r(i1,i2,i3,kd)=( 4.*un(i1,i2,i3,kd)     -2.*(un(i1+1,i2,i3,kd)+un(i1-1,i2,i3,kd)+un(i1,i2,i3+1,kd)+un(i1,i2,i3-1,kd))   +   (un(i1+1,i2,i3+1,kd)+un(i1-1,i2,i3+1,kd)+un(i1+1,i2,i3-1,kd)+un(i1-1,i2,i3-1,kd)) )/(dx(0)**2*dx(2)**2)
         unyyzz23r(i1,i2,i3,kd)=( 4.*un(i1,i2,i3,kd)     -2.*(un(i1,i2+1,i3,kd)  +un(i1,i2-1,i3,kd)+  un(i1,i2  ,i3+1,kd)+un(i1,i2  ,i3-1,kd))   +   (un(i1,i2+1,i3+1,kd)+un(i1,i2-1,i3+1,kd)+un(i1,i2+1,i3-1,kd)+un(i1,i2-1,i3-1,kd)) )/(dx(1)**2*dx(2)**2)
         ! 3D laplacian squared = un.xxxx + un.yyyy + un.zzzz + 2 (un.xxyy + un.xxzz + un.yyzz )
         unLapSq23r(i1,i2,i3,kd)= ( 6.*un(i1,i2,i3,kd)   - 4.*(un(i1+1,i2,i3,kd)+un(i1-1,i2,i3,kd))    +(un(i1+2,i2,i3,kd)+un(i1-2,i2,i3,kd)) )/(dx(0)**4) +( 6.*un(i1,i2,i3,kd)    -4.*(un(i1,i2+1,i3,kd)+un(i1,i2-1,i3,kd))    +(un(i1,i2+2,i3,kd)+un(i1,i2-2,i3,kd)) )/(dx(1)**4)  +( 6.*un(i1,i2,i3,kd)    -4.*(un(i1,i2,i3+1,kd)+un(i1,i2,i3-1,kd))    +(un(i1,i2,i3+2,kd)+un(i1,i2,i3-2,kd)) )/(dx(2)**4)  +( 8.*un(i1,i2,i3,kd)     -4.*(un(i1+1,i2,i3,kd)  +un(i1-1,i2,i3,kd)  +un(i1  ,i2+1,i3,kd)+un(i1  ,i2-1,i3,kd))   +2.*(un(i1+1,i2+1,i3,kd)+un(i1-1,i2+1,i3,kd)+un(i1+1,i2-1,i3,kd)+un(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)+( 8.*un(i1,i2,i3,kd)     -4.*(un(i1+1,i2,i3,kd)  +un(i1-1,i2,i3,kd)  +un(i1  ,i2,i3+1,kd)+un(i1  ,i2,i3-1,kd))   +2.*(un(i1+1,i2,i3+1,kd)+un(i1-1,i2,i3+1,kd)+un(i1+1,i2,i3-1,kd)+un(i1-1,i2,i3-1,kd)) )/(dx(0)**2*dx(2)**2)+( 8.*un(i1,i2,i3,kd)     -4.*(un(i1,i2+1,i3,kd)  +un(i1,i2-1,i3,kd)  +un(i1,i2  ,i3+1,kd)+un(i1,i2  ,i3-1,kd))   +2.*(un(i1,i2+1,i3+1,kd)+un(i1,i2-1,i3+1,kd)+un(i1,i2+1,i3-1,kd)+un(i1,i2-1,i3-1,kd)) )/(dx(1)**2*dx(2)**2)
         unr4(i1,i2,i3,kd)=(8.*(un(i1+1,i2,i3,kd)-un(i1-1,i2,i3,kd))-(un(i1+2,i2,i3,kd)-un(i1-2,i2,i3,kd)))*d14(0)
         uns4(i1,i2,i3,kd)=(8.*(un(i1,i2+1,i3,kd)-un(i1,i2-1,i3,kd))-(un(i1,i2+2,i3,kd)-un(i1,i2-2,i3,kd)))*d14(1)
         unt4(i1,i2,i3,kd)=(8.*(un(i1,i2,i3+1,kd)-un(i1,i2,i3-1,kd))-(un(i1,i2,i3+2,kd)-un(i1,i2,i3-2,kd)))*d14(2)
         unrr4(i1,i2,i3,kd)=(-30.*un(i1,i2,i3,kd)+16.*(un(i1+1,i2,i3,kd)+un(i1-1,i2,i3,kd))-(un(i1+2,i2,i3,kd)+un(i1-2,i2,i3,kd)) )*d24(0)
         unss4(i1,i2,i3,kd)=(-30.*un(i1,i2,i3,kd)+16.*(un(i1,i2+1,i3,kd)+un(i1,i2-1,i3,kd))-(un(i1,i2+2,i3,kd)+un(i1,i2-2,i3,kd)) )*d24(1)
         untt4(i1,i2,i3,kd)=(-30.*un(i1,i2,i3,kd)+16.*(un(i1,i2,i3+1,kd)+un(i1,i2,i3-1,kd))-(un(i1,i2,i3+2,kd)+un(i1,i2,i3-2,kd)) )*d24(2)
         unrs4(i1,i2,i3,kd)=(8.*(unr4(i1,i2+1,i3,kd)-unr4(i1,i2-1,i3,kd))-(unr4(i1,i2+2,i3,kd)-unr4(i1,i2-2,i3,kd)))*d14(1)
         unrt4(i1,i2,i3,kd)=(8.*(unr4(i1,i2,i3+1,kd)-unr4(i1,i2,i3-1,kd))-(unr4(i1,i2,i3+2,kd)-unr4(i1,i2,i3-2,kd)))*d14(2)
         unst4(i1,i2,i3,kd)=(8.*(uns4(i1,i2,i3+1,kd)-uns4(i1,i2,i3-1,kd))-(uns4(i1,i2,i3+2,kd)-uns4(i1,i2,i3-2,kd)))*d14(2)
         unx41(i1,i2,i3,kd)= rx(i1,i2,i3)*unr4(i1,i2,i3,kd)
         uny41(i1,i2,i3,kd)=0
         unz41(i1,i2,i3,kd)=0
         unx42(i1,i2,i3,kd)= rx(i1,i2,i3)*unr4(i1,i2,i3,kd)+sx(i1,i2,i3)*uns4(i1,i2,i3,kd)
         uny42(i1,i2,i3,kd)= ry(i1,i2,i3)*unr4(i1,i2,i3,kd)+sy(i1,i2,i3)*uns4(i1,i2,i3,kd)
         unz42(i1,i2,i3,kd)=0
         unx43(i1,i2,i3,kd)=rx(i1,i2,i3)*unr4(i1,i2,i3,kd)+sx(i1,i2,i3)*uns4(i1,i2,i3,kd)+tx(i1,i2,i3)*unt4(i1,i2,i3,kd)
         uny43(i1,i2,i3,kd)=ry(i1,i2,i3)*unr4(i1,i2,i3,kd)+sy(i1,i2,i3)*uns4(i1,i2,i3,kd)+ty(i1,i2,i3)*unt4(i1,i2,i3,kd)
         unz43(i1,i2,i3,kd)=rz(i1,i2,i3)*unr4(i1,i2,i3,kd)+sz(i1,i2,i3)*uns4(i1,i2,i3,kd)+tz(i1,i2,i3)*unt4(i1,i2,i3,kd)
         unxx41(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*unrr4(i1,i2,i3,kd)+(rxx42(i1,i2,i3))*unr4(i1,i2,i3,kd)
         unyy41(i1,i2,i3,kd)=0
         unxy41(i1,i2,i3,kd)=0
         unxz41(i1,i2,i3,kd)=0
         unyz41(i1,i2,i3,kd)=0
         unzz41(i1,i2,i3,kd)=0
         unlaplacian41(i1,i2,i3,kd)=unxx41(i1,i2,i3,kd)
         unxx42(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*unrr4(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3))*unrs4(i1,i2,i3,kd)+(sx(i1,i2,i3)**2)*unss4(i1,i2,i3,kd)+(rxx42(i1,i2,i3))*unr4(i1,i2,i3,kd)+(sxx42(i1,i2,i3))*uns4(i1,i2,i3,kd)
         unyy42(i1,i2,i3,kd)=(ry(i1,i2,i3)**2)*unrr4(i1,i2,i3,kd)+2.*(ry(i1,i2,i3)*sy(i1,i2,i3))*unrs4(i1,i2,i3,kd)+(sy(i1,i2,i3)**2)*unss4(i1,i2,i3,kd)+(ryy42(i1,i2,i3))*unr4(i1,i2,i3,kd)+(syy42(i1,i2,i3))*uns4(i1,i2,i3,kd)
         unxy42(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*unrr4(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*unrs4(i1,i2,i3,kd)+sx(i1,i2,i3)*sy(i1,i2,i3)*unss4(i1,i2,i3,kd)+rxy42(i1,i2,i3)*unr4(i1,i2,i3,kd)+sxy42(i1,i2,i3)*uns4(i1,i2,i3,kd)
         unxz42(i1,i2,i3,kd)=0
         unyz42(i1,i2,i3,kd)=0
         unzz42(i1,i2,i3,kd)=0
         unlaplacian42(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2)*unrr4(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*sy(i1,i2,i3))*unrs4(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**2)*unss4(i1,i2,i3,kd)+(rxx42(i1,i2,i3)+ryy42(i1,i2,i3))*unr4(i1,i2,i3,kd)+(sxx42(i1,i2,i3)+syy42(i1,i2,i3))*uns4(i1,i2,i3,kd)
         unxx43(i1,i2,i3,kd)=rx(i1,i2,i3)**2*unrr4(i1,i2,i3,kd)+sx(i1,i2,i3)**2*unss4(i1,i2,i3,kd)+tx(i1,i2,i3)**2*untt4(i1,i2,i3,kd)+2.*rx(i1,i2,i3)*sx(i1,i2,i3)*unrs4(i1,i2,i3,kd)+2.*rx(i1,i2,i3)*tx(i1,i2,i3)*unrt4(i1,i2,i3,kd)+2.*sx(i1,i2,i3)*tx(i1,i2,i3)*unst4(i1,i2,i3,kd)+rxx43(i1,i2,i3)*unr4(i1,i2,i3,kd)+sxx43(i1,i2,i3)*uns4(i1,i2,i3,kd)+txx43(i1,i2,i3)*unt4(i1,i2,i3,kd)
         unyy43(i1,i2,i3,kd)=ry(i1,i2,i3)**2*unrr4(i1,i2,i3,kd)+sy(i1,i2,i3)**2*unss4(i1,i2,i3,kd)+ty(i1,i2,i3)**2*untt4(i1,i2,i3,kd)+2.*ry(i1,i2,i3)*sy(i1,i2,i3)*unrs4(i1,i2,i3,kd)+2.*ry(i1,i2,i3)*ty(i1,i2,i3)*unrt4(i1,i2,i3,kd)+2.*sy(i1,i2,i3)*ty(i1,i2,i3)*unst4(i1,i2,i3,kd)+ryy43(i1,i2,i3)*unr4(i1,i2,i3,kd)+syy43(i1,i2,i3)*uns4(i1,i2,i3,kd)+tyy43(i1,i2,i3)*unt4(i1,i2,i3,kd)
         unzz43(i1,i2,i3,kd)=rz(i1,i2,i3)**2*unrr4(i1,i2,i3,kd)+sz(i1,i2,i3)**2*unss4(i1,i2,i3,kd)+tz(i1,i2,i3)**2*untt4(i1,i2,i3,kd)+2.*rz(i1,i2,i3)*sz(i1,i2,i3)*unrs4(i1,i2,i3,kd)+2.*rz(i1,i2,i3)*tz(i1,i2,i3)*unrt4(i1,i2,i3,kd)+2.*sz(i1,i2,i3)*tz(i1,i2,i3)*unst4(i1,i2,i3,kd)+rzz43(i1,i2,i3)*unr4(i1,i2,i3,kd)+szz43(i1,i2,i3)*uns4(i1,i2,i3,kd)+tzz43(i1,i2,i3)*unt4(i1,i2,i3,kd)
         unxy43(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*unrr4(i1,i2,i3,kd)+sx(i1,i2,i3)*sy(i1,i2,i3)*unss4(i1,i2,i3,kd)+tx(i1,i2,i3)*ty(i1,i2,i3)*untt4(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*unrs4(i1,i2,i3,kd)+(rx(i1,i2,i3)*ty(i1,i2,i3)+ry(i1,i2,i3)*tx(i1,i2,i3))*unrt4(i1,i2,i3,kd)+(sx(i1,i2,i3)*ty(i1,i2,i3)+sy(i1,i2,i3)*tx(i1,i2,i3))*unst4(i1,i2,i3,kd)+rxy43(i1,i2,i3)*unr4(i1,i2,i3,kd)+sxy43(i1,i2,i3)*uns4(i1,i2,i3,kd)+txy43(i1,i2,i3)*unt4(i1,i2,i3,kd)
         unxz43(i1,i2,i3,kd)=rx(i1,i2,i3)*rz(i1,i2,i3)*unrr4(i1,i2,i3,kd)+sx(i1,i2,i3)*sz(i1,i2,i3)*unss4(i1,i2,i3,kd)+tx(i1,i2,i3)*tz(i1,i2,i3)*untt4(i1,i2,i3,kd)+(rx(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,i2,i3)*sx(i1,i2,i3))*unrs4(i1,i2,i3,kd)+(rx(i1,i2,i3)*tz(i1,i2,i3)+rz(i1,i2,i3)*tx(i1,i2,i3))*unrt4(i1,i2,i3,kd)+(sx(i1,i2,i3)*tz(i1,i2,i3)+sz(i1,i2,i3)*tx(i1,i2,i3))*unst4(i1,i2,i3,kd)+rxz43(i1,i2,i3)*unr4(i1,i2,i3,kd)+sxz43(i1,i2,i3)*uns4(i1,i2,i3,kd)+txz43(i1,i2,i3)*unt4(i1,i2,i3,kd)
         unyz43(i1,i2,i3,kd)=ry(i1,i2,i3)*rz(i1,i2,i3)*unrr4(i1,i2,i3,kd)+sy(i1,i2,i3)*sz(i1,i2,i3)*unss4(i1,i2,i3,kd)+ty(i1,i2,i3)*tz(i1,i2,i3)*untt4(i1,i2,i3,kd)+(ry(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,i2,i3)*sy(i1,i2,i3))*unrs4(i1,i2,i3,kd)+(ry(i1,i2,i3)*tz(i1,i2,i3)+rz(i1,i2,i3)*ty(i1,i2,i3))*unrt4(i1,i2,i3,kd)+(sy(i1,i2,i3)*tz(i1,i2,i3)+sz(i1,i2,i3)*ty(i1,i2,i3))*unst4(i1,i2,i3,kd)+ryz43(i1,i2,i3)*unr4(i1,i2,i3,kd)+syz43(i1,i2,i3)*uns4(i1,i2,i3,kd)+tyz43(i1,i2,i3)*unt4(i1,i2,i3,kd)
         unlaplacian43(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2+rz(i1,i2,i3)**2)*unrr4(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**2+sz(i1,i2,i3)**2)*unss4(i1,i2,i3,kd)+(tx(i1,i2,i3)**2+ty(i1,i2,i3)**2+tz(i1,i2,i3)**2)*untt4(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*sy(i1,i2,i3)+rz(i1,i2,i3)*sz(i1,i2,i3))*unrs4(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*tx(i1,i2,i3)+ ry(i1,i2,i3)*ty(i1,i2,i3)+rz(i1,i2,i3)*tz(i1,i2,i3))*unrt4(i1,i2,i3,kd)+2.*(sx(i1,i2,i3)*tx(i1,i2,i3)+ sy(i1,i2,i3)*ty(i1,i2,i3)+sz(i1,i2,i3)*tz(i1,i2,i3))*unst4(i1,i2,i3,kd)+(rxx43(i1,i2,i3)+ryy43(i1,i2,i3)+rzz43(i1,i2,i3))*unr4(i1,i2,i3,kd)+(sxx43(i1,i2,i3)+syy43(i1,i2,i3)+szz43(i1,i2,i3))*uns4(i1,i2,i3,kd)+(txx43(i1,i2,i3)+tyy43(i1,i2,i3)+tzz43(i1,i2,i3))*unt4(i1,i2,i3,kd)
         !============================================================================================
         ! Define derivatives for a rectangular grid
         !
         !============================================================================================
         unx43r(i1,i2,i3,kd)=(8.*(un(i1+1,i2,i3,kd)-un(i1-1,i2,i3,kd))-(un(i1+2,i2,i3,kd)-un(i1-2,i2,i3,kd)))*h41(0)
         uny43r(i1,i2,i3,kd)=(8.*(un(i1,i2+1,i3,kd)-un(i1,i2-1,i3,kd))-(un(i1,i2+2,i3,kd)-un(i1,i2-2,i3,kd)))*h41(1)
         unz43r(i1,i2,i3,kd)=(8.*(un(i1,i2,i3+1,kd)-un(i1,i2,i3-1,kd))-(un(i1,i2,i3+2,kd)-un(i1,i2,i3-2,kd)))*h41(2)
         unxx43r(i1,i2,i3,kd)=( -30.*un(i1,i2,i3,kd)+16.*(un(i1+1,i2,i3,kd)+un(i1-1,i2,i3,kd))-(un(i1+2,i2,i3,kd)+un(i1-2,i2,i3,kd)) )*h42(0) 
         unyy43r(i1,i2,i3,kd)=( -30.*un(i1,i2,i3,kd)+16.*(un(i1,i2+1,i3,kd)+un(i1,i2-1,i3,kd))-(un(i1,i2+2,i3,kd)+un(i1,i2-2,i3,kd)) )*h42(1) 
         unzz43r(i1,i2,i3,kd)=( -30.*un(i1,i2,i3,kd)+16.*(un(i1,i2,i3+1,kd)+un(i1,i2,i3-1,kd))-(un(i1,i2,i3+2,kd)+un(i1,i2,i3-2,kd)) )*h42(2)
         unxy43r(i1,i2,i3,kd)=( (un(i1+2,i2+2,i3,kd)-un(i1-2,i2+2,i3,kd)- un(i1+2,i2-2,i3,kd)+un(i1-2,i2-2,i3,kd)) +8.*(un(i1-1,i2+2,i3,kd)-un(i1-1,i2-2,i3,kd)-un(i1+1,i2+2,i3,kd)+un(i1+1,i2-2,i3,kd) +un(i1+2,i2-1,i3,kd)-un(i1-2,i2-1,i3,kd)-un(i1+2,i2+1,i3,kd)+un(i1-2,i2+1,i3,kd))+64.*(un(i1+1,i2+1,i3,kd)-un(i1-1,i2+1,i3,kd)- un(i1+1,i2-1,i3,kd)+un(i1-1,i2-1,i3,kd)))*(h41(0)*h41(1))
         unxz43r(i1,i2,i3,kd)=( (un(i1+2,i2,i3+2,kd)-un(i1-2,i2,i3+2,kd)-un(i1+2,i2,i3-2,kd)+un(i1-2,i2,i3-2,kd)) +8.*(un(i1-1,i2,i3+2,kd)-un(i1-1,i2,i3-2,kd)-un(i1+1,i2,i3+2,kd)+un(i1+1,i2,i3-2,kd) +un(i1+2,i2,i3-1,kd)-un(i1-2,i2,i3-1,kd)- un(i1+2,i2,i3+1,kd)+un(i1-2,i2,i3+1,kd)) +64.*(un(i1+1,i2,i3+1,kd)-un(i1-1,i2,i3+1,kd)-un(i1+1,i2,i3-1,kd)+un(i1-1,i2,i3-1,kd)) )*(h41(0)*h41(2))
         unyz43r(i1,i2,i3,kd)=( (un(i1,i2+2,i3+2,kd)-un(i1,i2-2,i3+2,kd)-un(i1,i2+2,i3-2,kd)+un(i1,i2-2,i3-2,kd)) +8.*(un(i1,i2-1,i3+2,kd)-un(i1,i2-1,i3-2,kd)-un(i1,i2+1,i3+2,kd)+un(i1,i2+1,i3-2,kd) +un(i1,i2+2,i3-1,kd)-un(i1,i2-2,i3-1,kd)-un(i1,i2+2,i3+1,kd)+un(i1,i2-2,i3+1,kd)) +64.*(un(i1,i2+1,i3+1,kd)-un(i1,i2-1,i3+1,kd)-un(i1,i2+1,i3-1,kd)+un(i1,i2-1,i3-1,kd)) )*(h41(1)*h41(2))
         unx41r(i1,i2,i3,kd)= unx43r(i1,i2,i3,kd)
         uny41r(i1,i2,i3,kd)= uny43r(i1,i2,i3,kd)
         unz41r(i1,i2,i3,kd)= unz43r(i1,i2,i3,kd)
         unxx41r(i1,i2,i3,kd)= unxx43r(i1,i2,i3,kd)
         unyy41r(i1,i2,i3,kd)= unyy43r(i1,i2,i3,kd)
         unzz41r(i1,i2,i3,kd)= unzz43r(i1,i2,i3,kd)
         unxy41r(i1,i2,i3,kd)= unxy43r(i1,i2,i3,kd)
         unxz41r(i1,i2,i3,kd)= unxz43r(i1,i2,i3,kd)
         unyz41r(i1,i2,i3,kd)= unyz43r(i1,i2,i3,kd)
         unlaplacian41r(i1,i2,i3,kd)=unxx43r(i1,i2,i3,kd)
         unx42r(i1,i2,i3,kd)= unx43r(i1,i2,i3,kd)
         uny42r(i1,i2,i3,kd)= uny43r(i1,i2,i3,kd)
         unz42r(i1,i2,i3,kd)= unz43r(i1,i2,i3,kd)
         unxx42r(i1,i2,i3,kd)= unxx43r(i1,i2,i3,kd)
         unyy42r(i1,i2,i3,kd)= unyy43r(i1,i2,i3,kd)
         unzz42r(i1,i2,i3,kd)= unzz43r(i1,i2,i3,kd)
         unxy42r(i1,i2,i3,kd)= unxy43r(i1,i2,i3,kd)
         unxz42r(i1,i2,i3,kd)= unxz43r(i1,i2,i3,kd)
         unyz42r(i1,i2,i3,kd)= unyz43r(i1,i2,i3,kd)
         unlaplacian42r(i1,i2,i3,kd)=unxx43r(i1,i2,i3,kd)+unyy43r(i1,i2,i3,kd)
         unlaplacian43r(i1,i2,i3,kd)=unxx43r(i1,i2,i3,kd)+unyy43r(i1,i2,i3,kd)+unzz43r(i1,i2,i3,kd)
         umr2(i1,i2,i3,kd)=(um(i1+1,i2,i3,kd)-um(i1-1,i2,i3,kd))*d12(0)
         ums2(i1,i2,i3,kd)=(um(i1,i2+1,i3,kd)-um(i1,i2-1,i3,kd))*d12(1)
         umt2(i1,i2,i3,kd)=(um(i1,i2,i3+1,kd)-um(i1,i2,i3-1,kd))*d12(2)
         umrr2(i1,i2,i3,kd)=(-2.*um(i1,i2,i3,kd)+(um(i1+1,i2,i3,kd)+um(i1-1,i2,i3,kd)) )*d22(0)
         umss2(i1,i2,i3,kd)=(-2.*um(i1,i2,i3,kd)+(um(i1,i2+1,i3,kd)+um(i1,i2-1,i3,kd)) )*d22(1)
         umrs2(i1,i2,i3,kd)=(umr2(i1,i2+1,i3,kd)-umr2(i1,i2-1,i3,kd))*d12(1)
         umtt2(i1,i2,i3,kd)=(-2.*um(i1,i2,i3,kd)+(um(i1,i2,i3+1,kd)+um(i1,i2,i3-1,kd)) )*d22(2)
         umrt2(i1,i2,i3,kd)=(umr2(i1,i2,i3+1,kd)-umr2(i1,i2,i3-1,kd))*d12(2)
         umst2(i1,i2,i3,kd)=(ums2(i1,i2,i3+1,kd)-ums2(i1,i2,i3-1,kd))*d12(2)
         umrrr2(i1,i2,i3,kd)=(-2.*(um(i1+1,i2,i3,kd)-um(i1-1,i2,i3,kd))+(um(i1+2,i2,i3,kd)-um(i1-2,i2,i3,kd)) )*d22(0)*d12(0)
         umsss2(i1,i2,i3,kd)=(-2.*(um(i1,i2+1,i3,kd)-um(i1,i2-1,i3,kd))+(um(i1,i2+2,i3,kd)-um(i1,i2-2,i3,kd)) )*d22(1)*d12(1)
         umttt2(i1,i2,i3,kd)=(-2.*(um(i1,i2,i3+1,kd)-um(i1,i2,i3-1,kd))+(um(i1,i2,i3+2,kd)-um(i1,i2,i3-2,kd)) )*d22(2)*d12(2)
         umx21(i1,i2,i3,kd)= rx(i1,i2,i3)*umr2(i1,i2,i3,kd)
         umy21(i1,i2,i3,kd)=0
         umz21(i1,i2,i3,kd)=0
         umx22(i1,i2,i3,kd)= rx(i1,i2,i3)*umr2(i1,i2,i3,kd)+sx(i1,i2,i3)*ums2(i1,i2,i3,kd)
         umy22(i1,i2,i3,kd)= ry(i1,i2,i3)*umr2(i1,i2,i3,kd)+sy(i1,i2,i3)*ums2(i1,i2,i3,kd)
         umz22(i1,i2,i3,kd)=0
         umx23(i1,i2,i3,kd)=rx(i1,i2,i3)*umr2(i1,i2,i3,kd)+sx(i1,i2,i3)*ums2(i1,i2,i3,kd)+tx(i1,i2,i3)*umt2(i1,i2,i3,kd)
         umy23(i1,i2,i3,kd)=ry(i1,i2,i3)*umr2(i1,i2,i3,kd)+sy(i1,i2,i3)*ums2(i1,i2,i3,kd)+ty(i1,i2,i3)*umt2(i1,i2,i3,kd)
         umz23(i1,i2,i3,kd)=rz(i1,i2,i3)*umr2(i1,i2,i3,kd)+sz(i1,i2,i3)*ums2(i1,i2,i3,kd)+tz(i1,i2,i3)*umt2(i1,i2,i3,kd)
         umxx21(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*umrr2(i1,i2,i3,kd)+(rxx22(i1,i2,i3))*umr2(i1,i2,i3,kd)
         umyy21(i1,i2,i3,kd)=0
         umxy21(i1,i2,i3,kd)=0
         umxz21(i1,i2,i3,kd)=0
         umyz21(i1,i2,i3,kd)=0
         umzz21(i1,i2,i3,kd)=0
         umlaplacian21(i1,i2,i3,kd)=umxx21(i1,i2,i3,kd)
         umxx22(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*umrr2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3))*umrs2(i1,i2,i3,kd)+(sx(i1,i2,i3)**2)*umss2(i1,i2,i3,kd)+(rxx22(i1,i2,i3))*umr2(i1,i2,i3,kd)+(sxx22(i1,i2,i3))*ums2(i1,i2,i3,kd)
         umyy22(i1,i2,i3,kd)=(ry(i1,i2,i3)**2)*umrr2(i1,i2,i3,kd)+2.*(ry(i1,i2,i3)*sy(i1,i2,i3))*umrs2(i1,i2,i3,kd)+(sy(i1,i2,i3)**2)*umss2(i1,i2,i3,kd)+(ryy22(i1,i2,i3))*umr2(i1,i2,i3,kd)+(syy22(i1,i2,i3))*ums2(i1,i2,i3,kd)
         umxy22(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*umrr2(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*umrs2(i1,i2,i3,kd)+sx(i1,i2,i3)*sy(i1,i2,i3)*umss2(i1,i2,i3,kd)+rxy22(i1,i2,i3)*umr2(i1,i2,i3,kd)+sxy22(i1,i2,i3)*ums2(i1,i2,i3,kd)
         umxz22(i1,i2,i3,kd)=0
         umyz22(i1,i2,i3,kd)=0
         umzz22(i1,i2,i3,kd)=0
         umlaplacian22(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2)*umrr2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*sy(i1,i2,i3))*umrs2(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**2)*umss2(i1,i2,i3,kd)+(rxx22(i1,i2,i3)+ryy22(i1,i2,i3))*umr2(i1,i2,i3,kd)+(sxx22(i1,i2,i3)+syy22(i1,i2,i3))*ums2(i1,i2,i3,kd)
         umxx23(i1,i2,i3,kd)=rx(i1,i2,i3)**2*umrr2(i1,i2,i3,kd)+sx(i1,i2,i3)**2*umss2(i1,i2,i3,kd)+tx(i1,i2,i3)**2*umtt2(i1,i2,i3,kd)+2.*rx(i1,i2,i3)*sx(i1,i2,i3)*umrs2(i1,i2,i3,kd)+2.*rx(i1,i2,i3)*tx(i1,i2,i3)*umrt2(i1,i2,i3,kd)+2.*sx(i1,i2,i3)*tx(i1,i2,i3)*umst2(i1,i2,i3,kd)+rxx23(i1,i2,i3)*umr2(i1,i2,i3,kd)+sxx23(i1,i2,i3)*ums2(i1,i2,i3,kd)+txx23(i1,i2,i3)*umt2(i1,i2,i3,kd)
         umyy23(i1,i2,i3,kd)=ry(i1,i2,i3)**2*umrr2(i1,i2,i3,kd)+sy(i1,i2,i3)**2*umss2(i1,i2,i3,kd)+ty(i1,i2,i3)**2*umtt2(i1,i2,i3,kd)+2.*ry(i1,i2,i3)*sy(i1,i2,i3)*umrs2(i1,i2,i3,kd)+2.*ry(i1,i2,i3)*ty(i1,i2,i3)*umrt2(i1,i2,i3,kd)+2.*sy(i1,i2,i3)*ty(i1,i2,i3)*umst2(i1,i2,i3,kd)+ryy23(i1,i2,i3)*umr2(i1,i2,i3,kd)+syy23(i1,i2,i3)*ums2(i1,i2,i3,kd)+tyy23(i1,i2,i3)*umt2(i1,i2,i3,kd)
         umzz23(i1,i2,i3,kd)=rz(i1,i2,i3)**2*umrr2(i1,i2,i3,kd)+sz(i1,i2,i3)**2*umss2(i1,i2,i3,kd)+tz(i1,i2,i3)**2*umtt2(i1,i2,i3,kd)+2.*rz(i1,i2,i3)*sz(i1,i2,i3)*umrs2(i1,i2,i3,kd)+2.*rz(i1,i2,i3)*tz(i1,i2,i3)*umrt2(i1,i2,i3,kd)+2.*sz(i1,i2,i3)*tz(i1,i2,i3)*umst2(i1,i2,i3,kd)+rzz23(i1,i2,i3)*umr2(i1,i2,i3,kd)+szz23(i1,i2,i3)*ums2(i1,i2,i3,kd)+tzz23(i1,i2,i3)*umt2(i1,i2,i3,kd)
         umxy23(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*umrr2(i1,i2,i3,kd)+sx(i1,i2,i3)*sy(i1,i2,i3)*umss2(i1,i2,i3,kd)+tx(i1,i2,i3)*ty(i1,i2,i3)*umtt2(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*umrs2(i1,i2,i3,kd)+(rx(i1,i2,i3)*ty(i1,i2,i3)+ry(i1,i2,i3)*tx(i1,i2,i3))*umrt2(i1,i2,i3,kd)+(sx(i1,i2,i3)*ty(i1,i2,i3)+sy(i1,i2,i3)*tx(i1,i2,i3))*umst2(i1,i2,i3,kd)+rxy23(i1,i2,i3)*umr2(i1,i2,i3,kd)+sxy23(i1,i2,i3)*ums2(i1,i2,i3,kd)+txy23(i1,i2,i3)*umt2(i1,i2,i3,kd)
         umxz23(i1,i2,i3,kd)=rx(i1,i2,i3)*rz(i1,i2,i3)*umrr2(i1,i2,i3,kd)+sx(i1,i2,i3)*sz(i1,i2,i3)*umss2(i1,i2,i3,kd)+tx(i1,i2,i3)*tz(i1,i2,i3)*umtt2(i1,i2,i3,kd)+(rx(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,i2,i3)*sx(i1,i2,i3))*umrs2(i1,i2,i3,kd)+(rx(i1,i2,i3)*tz(i1,i2,i3)+rz(i1,i2,i3)*tx(i1,i2,i3))*umrt2(i1,i2,i3,kd)+(sx(i1,i2,i3)*tz(i1,i2,i3)+sz(i1,i2,i3)*tx(i1,i2,i3))*umst2(i1,i2,i3,kd)+rxz23(i1,i2,i3)*umr2(i1,i2,i3,kd)+sxz23(i1,i2,i3)*ums2(i1,i2,i3,kd)+txz23(i1,i2,i3)*umt2(i1,i2,i3,kd)
         umyz23(i1,i2,i3,kd)=ry(i1,i2,i3)*rz(i1,i2,i3)*umrr2(i1,i2,i3,kd)+sy(i1,i2,i3)*sz(i1,i2,i3)*umss2(i1,i2,i3,kd)+ty(i1,i2,i3)*tz(i1,i2,i3)*umtt2(i1,i2,i3,kd)+(ry(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,i2,i3)*sy(i1,i2,i3))*umrs2(i1,i2,i3,kd)+(ry(i1,i2,i3)*tz(i1,i2,i3)+rz(i1,i2,i3)*ty(i1,i2,i3))*umrt2(i1,i2,i3,kd)+(sy(i1,i2,i3)*tz(i1,i2,i3)+sz(i1,i2,i3)*ty(i1,i2,i3))*umst2(i1,i2,i3,kd)+ryz23(i1,i2,i3)*umr2(i1,i2,i3,kd)+syz23(i1,i2,i3)*ums2(i1,i2,i3,kd)+tyz23(i1,i2,i3)*umt2(i1,i2,i3,kd)
         umlaplacian23(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2+rz(i1,i2,i3)**2)*umrr2(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**2+sz(i1,i2,i3)**2)*umss2(i1,i2,i3,kd)+(tx(i1,i2,i3)**2+ty(i1,i2,i3)**2+tz(i1,i2,i3)**2)*umtt2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*sy(i1,i2,i3)+rz(i1,i2,i3)*sz(i1,i2,i3))*umrs2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*tx(i1,i2,i3)+ ry(i1,i2,i3)*ty(i1,i2,i3)+rz(i1,i2,i3)*tz(i1,i2,i3))*umrt2(i1,i2,i3,kd)+2.*(sx(i1,i2,i3)*tx(i1,i2,i3)+ sy(i1,i2,i3)*ty(i1,i2,i3)+sz(i1,i2,i3)*tz(i1,i2,i3))*umst2(i1,i2,i3,kd)+(rxx23(i1,i2,i3)+ryy23(i1,i2,i3)+rzz23(i1,i2,i3))*umr2(i1,i2,i3,kd)+(sxx23(i1,i2,i3)+syy23(i1,i2,i3)+szz23(i1,i2,i3))*ums2(i1,i2,i3,kd)+(txx23(i1,i2,i3)+tyy23(i1,i2,i3)+tzz23(i1,i2,i3))*umt2(i1,i2,i3,kd)
         !============================================================================================
         ! Define derivatives for a rectangular grid
         !
         !============================================================================================
         umx23r(i1,i2,i3,kd)=(um(i1+1,i2,i3,kd)-um(i1-1,i2,i3,kd))*h12(0)
         umy23r(i1,i2,i3,kd)=(um(i1,i2+1,i3,kd)-um(i1,i2-1,i3,kd))*h12(1)
         umz23r(i1,i2,i3,kd)=(um(i1,i2,i3+1,kd)-um(i1,i2,i3-1,kd))*h12(2)
         umxx23r(i1,i2,i3,kd)=(-2.*um(i1,i2,i3,kd)+(um(i1+1,i2,i3,kd)+um(i1-1,i2,i3,kd)) )*h22(0)
         umyy23r(i1,i2,i3,kd)=(-2.*um(i1,i2,i3,kd)+(um(i1,i2+1,i3,kd)+um(i1,i2-1,i3,kd)) )*h22(1)
         umxy23r(i1,i2,i3,kd)=(umx23r(i1,i2+1,i3,kd)-umx23r(i1,i2-1,i3,kd))*h12(1)
         umzz23r(i1,i2,i3,kd)=(-2.*um(i1,i2,i3,kd)+(um(i1,i2,i3+1,kd)+um(i1,i2,i3-1,kd)) )*h22(2)
         umxz23r(i1,i2,i3,kd)=(umx23r(i1,i2,i3+1,kd)-umx23r(i1,i2,i3-1,kd))*h12(2)
         umyz23r(i1,i2,i3,kd)=(umy23r(i1,i2,i3+1,kd)-umy23r(i1,i2,i3-1,kd))*h12(2)
         umx21r(i1,i2,i3,kd)= umx23r(i1,i2,i3,kd)
         umy21r(i1,i2,i3,kd)= umy23r(i1,i2,i3,kd)
         umz21r(i1,i2,i3,kd)= umz23r(i1,i2,i3,kd)
         umxx21r(i1,i2,i3,kd)= umxx23r(i1,i2,i3,kd)
         umyy21r(i1,i2,i3,kd)= umyy23r(i1,i2,i3,kd)
         umzz21r(i1,i2,i3,kd)= umzz23r(i1,i2,i3,kd)
         umxy21r(i1,i2,i3,kd)= umxy23r(i1,i2,i3,kd)
         umxz21r(i1,i2,i3,kd)= umxz23r(i1,i2,i3,kd)
         umyz21r(i1,i2,i3,kd)= umyz23r(i1,i2,i3,kd)
         umlaplacian21r(i1,i2,i3,kd)=umxx23r(i1,i2,i3,kd)
         umx22r(i1,i2,i3,kd)= umx23r(i1,i2,i3,kd)
         umy22r(i1,i2,i3,kd)= umy23r(i1,i2,i3,kd)
         umz22r(i1,i2,i3,kd)= umz23r(i1,i2,i3,kd)
         umxx22r(i1,i2,i3,kd)= umxx23r(i1,i2,i3,kd)
         umyy22r(i1,i2,i3,kd)= umyy23r(i1,i2,i3,kd)
         umzz22r(i1,i2,i3,kd)= umzz23r(i1,i2,i3,kd)
         umxy22r(i1,i2,i3,kd)= umxy23r(i1,i2,i3,kd)
         umxz22r(i1,i2,i3,kd)= umxz23r(i1,i2,i3,kd)
         umyz22r(i1,i2,i3,kd)= umyz23r(i1,i2,i3,kd)
         umlaplacian22r(i1,i2,i3,kd)=umxx23r(i1,i2,i3,kd)+umyy23r(i1,i2,i3,kd)
         umlaplacian23r(i1,i2,i3,kd)=umxx23r(i1,i2,i3,kd)+umyy23r(i1,i2,i3,kd)+umzz23r(i1,i2,i3,kd)
         umxxx22r(i1,i2,i3,kd)=(-2.*(um(i1+1,i2,i3,kd)-um(i1-1,i2,i3,kd))+(um(i1+2,i2,i3,kd)-um(i1-2,i2,i3,kd)) )*h22(0)*h12(0)
         umyyy22r(i1,i2,i3,kd)=(-2.*(um(i1,i2+1,i3,kd)-um(i1,i2-1,i3,kd))+(um(i1,i2+2,i3,kd)-um(i1,i2-2,i3,kd)) )*h22(1)*h12(1)
         umxxy22r(i1,i2,i3,kd)=( umxx22r(i1,i2+1,i3,kd)-umxx22r(i1,i2-1,i3,kd))/(2.*dx(1))
         umxyy22r(i1,i2,i3,kd)=( umyy22r(i1+1,i2,i3,kd)-umyy22r(i1-1,i2,i3,kd))/(2.*dx(0))
         umxxxx22r(i1,i2,i3,kd)=(6.*um(i1,i2,i3,kd)-4.*(um(i1+1,i2,i3,kd)+um(i1-1,i2,i3,kd))+(um(i1+2,i2,i3,kd)+um(i1-2,i2,i3,kd)) )/(dx(0)**4)
         umyyyy22r(i1,i2,i3,kd)=(6.*um(i1,i2,i3,kd)-4.*(um(i1,i2+1,i3,kd)+um(i1,i2-1,i3,kd))+(um(i1,i2+2,i3,kd)+um(i1,i2-2,i3,kd)) )/(dx(1)**4)
         umxxyy22r(i1,i2,i3,kd)=( 4.*um(i1,i2,i3,kd)     -2.*(um(i1+1,i2,i3,kd)+um(i1-1,i2,i3,kd)+um(i1,i2+1,i3,kd)+um(i1,i2-1,i3,kd))   +   (um(i1+1,i2+1,i3,kd)+um(i1-1,i2+1,i3,kd)+um(i1+1,i2-1,i3,kd)+um(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)
         ! 2D laplacian squared = um.xxxx + 2 um.xxyy + um.yyyy
         umLapSq22r(i1,i2,i3,kd)= ( 6.*um(i1,i2,i3,kd)   - 4.*(um(i1+1,i2,i3,kd)+um(i1-1,i2,i3,kd))    +(um(i1+2,i2,i3,kd)+um(i1-2,i2,i3,kd)) )/(dx(0)**4) +( 6.*um(i1,i2,i3,kd)    -4.*(um(i1,i2+1,i3,kd)+um(i1,i2-1,i3,kd))    +(um(i1,i2+2,i3,kd)+um(i1,i2-2,i3,kd)) )/(dx(1)**4)  +( 8.*um(i1,i2,i3,kd)     -4.*(um(i1+1,i2,i3,kd)+um(i1-1,i2,i3,kd)+um(i1,i2+1,i3,kd)+um(i1,i2-1,i3,kd))   +2.*(um(i1+1,i2+1,i3,kd)+um(i1-1,i2+1,i3,kd)+um(i1+1,i2-1,i3,kd)+um(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)
         umxxx23r(i1,i2,i3,kd)=(-2.*(um(i1+1,i2,i3,kd)-um(i1-1,i2,i3,kd))+(um(i1+2,i2,i3,kd)-um(i1-2,i2,i3,kd)) )*h22(0)*h12(0)
         umyyy23r(i1,i2,i3,kd)=(-2.*(um(i1,i2+1,i3,kd)-um(i1,i2-1,i3,kd))+(um(i1,i2+2,i3,kd)-um(i1,i2-2,i3,kd)) )*h22(1)*h12(1)
         umzzz23r(i1,i2,i3,kd)=(-2.*(um(i1,i2,i3+1,kd)-um(i1,i2,i3-1,kd))+(um(i1,i2,i3+2,kd)-um(i1,i2,i3-2,kd)) )*h22(1)*h12(2)
         umxxy23r(i1,i2,i3,kd)=( umxx22r(i1,i2+1,i3,kd)-umxx22r(i1,i2-1,i3,kd))/(2.*dx(1))
         umxyy23r(i1,i2,i3,kd)=( umyy22r(i1+1,i2,i3,kd)-umyy22r(i1-1,i2,i3,kd))/(2.*dx(0))
         umxxz23r(i1,i2,i3,kd)=( umxx22r(i1,i2,i3+1,kd)-umxx22r(i1,i2,i3-1,kd))/(2.*dx(2))
         umyyz23r(i1,i2,i3,kd)=( umyy22r(i1,i2,i3+1,kd)-umyy22r(i1,i2,i3-1,kd))/(2.*dx(2))
         umxzz23r(i1,i2,i3,kd)=( umzz22r(i1+1,i2,i3,kd)-umzz22r(i1-1,i2,i3,kd))/(2.*dx(0))
         umyzz23r(i1,i2,i3,kd)=( umzz22r(i1,i2+1,i3,kd)-umzz22r(i1,i2-1,i3,kd))/(2.*dx(1))
         umxxxx23r(i1,i2,i3,kd)=(6.*um(i1,i2,i3,kd)-4.*(um(i1+1,i2,i3,kd)+um(i1-1,i2,i3,kd))+(um(i1+2,i2,i3,kd)+um(i1-2,i2,i3,kd)) )/(dx(0)**4)
         umyyyy23r(i1,i2,i3,kd)=(6.*um(i1,i2,i3,kd)-4.*(um(i1,i2+1,i3,kd)+um(i1,i2-1,i3,kd))+(um(i1,i2+2,i3,kd)+um(i1,i2-2,i3,kd)) )/(dx(1)**4)
         umzzzz23r(i1,i2,i3,kd)=(6.*um(i1,i2,i3,kd)-4.*(um(i1,i2,i3+1,kd)+um(i1,i2,i3-1,kd))+(um(i1,i2,i3+2,kd)+um(i1,i2,i3-2,kd)) )/(dx(2)**4)
         umxxyy23r(i1,i2,i3,kd)=( 4.*um(i1,i2,i3,kd)     -2.*(um(i1+1,i2,i3,kd)+um(i1-1,i2,i3,kd)+um(i1,i2+1,i3,kd)+um(i1,i2-1,i3,kd))   +   (um(i1+1,i2+1,i3,kd)+um(i1-1,i2+1,i3,kd)+um(i1+1,i2-1,i3,kd)+um(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)
         umxxzz23r(i1,i2,i3,kd)=( 4.*um(i1,i2,i3,kd)     -2.*(um(i1+1,i2,i3,kd)+um(i1-1,i2,i3,kd)+um(i1,i2,i3+1,kd)+um(i1,i2,i3-1,kd))   +   (um(i1+1,i2,i3+1,kd)+um(i1-1,i2,i3+1,kd)+um(i1+1,i2,i3-1,kd)+um(i1-1,i2,i3-1,kd)) )/(dx(0)**2*dx(2)**2)
         umyyzz23r(i1,i2,i3,kd)=( 4.*um(i1,i2,i3,kd)     -2.*(um(i1,i2+1,i3,kd)  +um(i1,i2-1,i3,kd)+  um(i1,i2  ,i3+1,kd)+um(i1,i2  ,i3-1,kd))   +   (um(i1,i2+1,i3+1,kd)+um(i1,i2-1,i3+1,kd)+um(i1,i2+1,i3-1,kd)+um(i1,i2-1,i3-1,kd)) )/(dx(1)**2*dx(2)**2)
         ! 3D laplacian squared = um.xxxx + um.yyyy + um.zzzz + 2 (um.xxyy + um.xxzz + um.yyzz )
         umLapSq23r(i1,i2,i3,kd)= ( 6.*um(i1,i2,i3,kd)   - 4.*(um(i1+1,i2,i3,kd)+um(i1-1,i2,i3,kd))    +(um(i1+2,i2,i3,kd)+um(i1-2,i2,i3,kd)) )/(dx(0)**4) +( 6.*um(i1,i2,i3,kd)    -4.*(um(i1,i2+1,i3,kd)+um(i1,i2-1,i3,kd))    +(um(i1,i2+2,i3,kd)+um(i1,i2-2,i3,kd)) )/(dx(1)**4)  +( 6.*um(i1,i2,i3,kd)    -4.*(um(i1,i2,i3+1,kd)+um(i1,i2,i3-1,kd))    +(um(i1,i2,i3+2,kd)+um(i1,i2,i3-2,kd)) )/(dx(2)**4)  +( 8.*um(i1,i2,i3,kd)     -4.*(um(i1+1,i2,i3,kd)  +um(i1-1,i2,i3,kd)  +um(i1  ,i2+1,i3,kd)+um(i1  ,i2-1,i3,kd))   +2.*(um(i1+1,i2+1,i3,kd)+um(i1-1,i2+1,i3,kd)+um(i1+1,i2-1,i3,kd)+um(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)+( 8.*um(i1,i2,i3,kd)     -4.*(um(i1+1,i2,i3,kd)  +um(i1-1,i2,i3,kd)  +um(i1  ,i2,i3+1,kd)+um(i1  ,i2,i3-1,kd))   +2.*(um(i1+1,i2,i3+1,kd)+um(i1-1,i2,i3+1,kd)+um(i1+1,i2,i3-1,kd)+um(i1-1,i2,i3-1,kd)) )/(dx(0)**2*dx(2)**2)+( 8.*um(i1,i2,i3,kd)     -4.*(um(i1,i2+1,i3,kd)  +um(i1,i2-1,i3,kd)  +um(i1,i2  ,i3+1,kd)+um(i1,i2  ,i3-1,kd))   +2.*(um(i1,i2+1,i3+1,kd)+um(i1,i2-1,i3+1,kd)+um(i1,i2+1,i3-1,kd)+um(i1,i2-1,i3-1,kd)) )/(dx(1)**2*dx(2)**2)
         umr4(i1,i2,i3,kd)=(8.*(um(i1+1,i2,i3,kd)-um(i1-1,i2,i3,kd))-(um(i1+2,i2,i3,kd)-um(i1-2,i2,i3,kd)))*d14(0)
         ums4(i1,i2,i3,kd)=(8.*(um(i1,i2+1,i3,kd)-um(i1,i2-1,i3,kd))-(um(i1,i2+2,i3,kd)-um(i1,i2-2,i3,kd)))*d14(1)
         umt4(i1,i2,i3,kd)=(8.*(um(i1,i2,i3+1,kd)-um(i1,i2,i3-1,kd))-(um(i1,i2,i3+2,kd)-um(i1,i2,i3-2,kd)))*d14(2)
         umrr4(i1,i2,i3,kd)=(-30.*um(i1,i2,i3,kd)+16.*(um(i1+1,i2,i3,kd)+um(i1-1,i2,i3,kd))-(um(i1+2,i2,i3,kd)+um(i1-2,i2,i3,kd)) )*d24(0)
         umss4(i1,i2,i3,kd)=(-30.*um(i1,i2,i3,kd)+16.*(um(i1,i2+1,i3,kd)+um(i1,i2-1,i3,kd))-(um(i1,i2+2,i3,kd)+um(i1,i2-2,i3,kd)) )*d24(1)
         umtt4(i1,i2,i3,kd)=(-30.*um(i1,i2,i3,kd)+16.*(um(i1,i2,i3+1,kd)+um(i1,i2,i3-1,kd))-(um(i1,i2,i3+2,kd)+um(i1,i2,i3-2,kd)) )*d24(2)
         umrs4(i1,i2,i3,kd)=(8.*(umr4(i1,i2+1,i3,kd)-umr4(i1,i2-1,i3,kd))-(umr4(i1,i2+2,i3,kd)-umr4(i1,i2-2,i3,kd)))*d14(1)
         umrt4(i1,i2,i3,kd)=(8.*(umr4(i1,i2,i3+1,kd)-umr4(i1,i2,i3-1,kd))-(umr4(i1,i2,i3+2,kd)-umr4(i1,i2,i3-2,kd)))*d14(2)
         umst4(i1,i2,i3,kd)=(8.*(ums4(i1,i2,i3+1,kd)-ums4(i1,i2,i3-1,kd))-(ums4(i1,i2,i3+2,kd)-ums4(i1,i2,i3-2,kd)))*d14(2)
         umx41(i1,i2,i3,kd)= rx(i1,i2,i3)*umr4(i1,i2,i3,kd)
         umy41(i1,i2,i3,kd)=0
         umz41(i1,i2,i3,kd)=0
         umx42(i1,i2,i3,kd)= rx(i1,i2,i3)*umr4(i1,i2,i3,kd)+sx(i1,i2,i3)*ums4(i1,i2,i3,kd)
         umy42(i1,i2,i3,kd)= ry(i1,i2,i3)*umr4(i1,i2,i3,kd)+sy(i1,i2,i3)*ums4(i1,i2,i3,kd)
         umz42(i1,i2,i3,kd)=0
         umx43(i1,i2,i3,kd)=rx(i1,i2,i3)*umr4(i1,i2,i3,kd)+sx(i1,i2,i3)*ums4(i1,i2,i3,kd)+tx(i1,i2,i3)*umt4(i1,i2,i3,kd)
         umy43(i1,i2,i3,kd)=ry(i1,i2,i3)*umr4(i1,i2,i3,kd)+sy(i1,i2,i3)*ums4(i1,i2,i3,kd)+ty(i1,i2,i3)*umt4(i1,i2,i3,kd)
         umz43(i1,i2,i3,kd)=rz(i1,i2,i3)*umr4(i1,i2,i3,kd)+sz(i1,i2,i3)*ums4(i1,i2,i3,kd)+tz(i1,i2,i3)*umt4(i1,i2,i3,kd)
         umxx41(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*umrr4(i1,i2,i3,kd)+(rxx42(i1,i2,i3))*umr4(i1,i2,i3,kd)
         umyy41(i1,i2,i3,kd)=0
         umxy41(i1,i2,i3,kd)=0
         umxz41(i1,i2,i3,kd)=0
         umyz41(i1,i2,i3,kd)=0
         umzz41(i1,i2,i3,kd)=0
         umlaplacian41(i1,i2,i3,kd)=umxx41(i1,i2,i3,kd)
         umxx42(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*umrr4(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3))*umrs4(i1,i2,i3,kd)+(sx(i1,i2,i3)**2)*umss4(i1,i2,i3,kd)+(rxx42(i1,i2,i3))*umr4(i1,i2,i3,kd)+(sxx42(i1,i2,i3))*ums4(i1,i2,i3,kd)
         umyy42(i1,i2,i3,kd)=(ry(i1,i2,i3)**2)*umrr4(i1,i2,i3,kd)+2.*(ry(i1,i2,i3)*sy(i1,i2,i3))*umrs4(i1,i2,i3,kd)+(sy(i1,i2,i3)**2)*umss4(i1,i2,i3,kd)+(ryy42(i1,i2,i3))*umr4(i1,i2,i3,kd)+(syy42(i1,i2,i3))*ums4(i1,i2,i3,kd)
         umxy42(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*umrr4(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*umrs4(i1,i2,i3,kd)+sx(i1,i2,i3)*sy(i1,i2,i3)*umss4(i1,i2,i3,kd)+rxy42(i1,i2,i3)*umr4(i1,i2,i3,kd)+sxy42(i1,i2,i3)*ums4(i1,i2,i3,kd)
         umxz42(i1,i2,i3,kd)=0
         umyz42(i1,i2,i3,kd)=0
         umzz42(i1,i2,i3,kd)=0
         umlaplacian42(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2)*umrr4(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*sy(i1,i2,i3))*umrs4(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**2)*umss4(i1,i2,i3,kd)+(rxx42(i1,i2,i3)+ryy42(i1,i2,i3))*umr4(i1,i2,i3,kd)+(sxx42(i1,i2,i3)+syy42(i1,i2,i3))*ums4(i1,i2,i3,kd)
         umxx43(i1,i2,i3,kd)=rx(i1,i2,i3)**2*umrr4(i1,i2,i3,kd)+sx(i1,i2,i3)**2*umss4(i1,i2,i3,kd)+tx(i1,i2,i3)**2*umtt4(i1,i2,i3,kd)+2.*rx(i1,i2,i3)*sx(i1,i2,i3)*umrs4(i1,i2,i3,kd)+2.*rx(i1,i2,i3)*tx(i1,i2,i3)*umrt4(i1,i2,i3,kd)+2.*sx(i1,i2,i3)*tx(i1,i2,i3)*umst4(i1,i2,i3,kd)+rxx43(i1,i2,i3)*umr4(i1,i2,i3,kd)+sxx43(i1,i2,i3)*ums4(i1,i2,i3,kd)+txx43(i1,i2,i3)*umt4(i1,i2,i3,kd)
         umyy43(i1,i2,i3,kd)=ry(i1,i2,i3)**2*umrr4(i1,i2,i3,kd)+sy(i1,i2,i3)**2*umss4(i1,i2,i3,kd)+ty(i1,i2,i3)**2*umtt4(i1,i2,i3,kd)+2.*ry(i1,i2,i3)*sy(i1,i2,i3)*umrs4(i1,i2,i3,kd)+2.*ry(i1,i2,i3)*ty(i1,i2,i3)*umrt4(i1,i2,i3,kd)+2.*sy(i1,i2,i3)*ty(i1,i2,i3)*umst4(i1,i2,i3,kd)+ryy43(i1,i2,i3)*umr4(i1,i2,i3,kd)+syy43(i1,i2,i3)*ums4(i1,i2,i3,kd)+tyy43(i1,i2,i3)*umt4(i1,i2,i3,kd)
         umzz43(i1,i2,i3,kd)=rz(i1,i2,i3)**2*umrr4(i1,i2,i3,kd)+sz(i1,i2,i3)**2*umss4(i1,i2,i3,kd)+tz(i1,i2,i3)**2*umtt4(i1,i2,i3,kd)+2.*rz(i1,i2,i3)*sz(i1,i2,i3)*umrs4(i1,i2,i3,kd)+2.*rz(i1,i2,i3)*tz(i1,i2,i3)*umrt4(i1,i2,i3,kd)+2.*sz(i1,i2,i3)*tz(i1,i2,i3)*umst4(i1,i2,i3,kd)+rzz43(i1,i2,i3)*umr4(i1,i2,i3,kd)+szz43(i1,i2,i3)*ums4(i1,i2,i3,kd)+tzz43(i1,i2,i3)*umt4(i1,i2,i3,kd)
         umxy43(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*umrr4(i1,i2,i3,kd)+sx(i1,i2,i3)*sy(i1,i2,i3)*umss4(i1,i2,i3,kd)+tx(i1,i2,i3)*ty(i1,i2,i3)*umtt4(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*umrs4(i1,i2,i3,kd)+(rx(i1,i2,i3)*ty(i1,i2,i3)+ry(i1,i2,i3)*tx(i1,i2,i3))*umrt4(i1,i2,i3,kd)+(sx(i1,i2,i3)*ty(i1,i2,i3)+sy(i1,i2,i3)*tx(i1,i2,i3))*umst4(i1,i2,i3,kd)+rxy43(i1,i2,i3)*umr4(i1,i2,i3,kd)+sxy43(i1,i2,i3)*ums4(i1,i2,i3,kd)+txy43(i1,i2,i3)*umt4(i1,i2,i3,kd)
         umxz43(i1,i2,i3,kd)=rx(i1,i2,i3)*rz(i1,i2,i3)*umrr4(i1,i2,i3,kd)+sx(i1,i2,i3)*sz(i1,i2,i3)*umss4(i1,i2,i3,kd)+tx(i1,i2,i3)*tz(i1,i2,i3)*umtt4(i1,i2,i3,kd)+(rx(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,i2,i3)*sx(i1,i2,i3))*umrs4(i1,i2,i3,kd)+(rx(i1,i2,i3)*tz(i1,i2,i3)+rz(i1,i2,i3)*tx(i1,i2,i3))*umrt4(i1,i2,i3,kd)+(sx(i1,i2,i3)*tz(i1,i2,i3)+sz(i1,i2,i3)*tx(i1,i2,i3))*umst4(i1,i2,i3,kd)+rxz43(i1,i2,i3)*umr4(i1,i2,i3,kd)+sxz43(i1,i2,i3)*ums4(i1,i2,i3,kd)+txz43(i1,i2,i3)*umt4(i1,i2,i3,kd)
         umyz43(i1,i2,i3,kd)=ry(i1,i2,i3)*rz(i1,i2,i3)*umrr4(i1,i2,i3,kd)+sy(i1,i2,i3)*sz(i1,i2,i3)*umss4(i1,i2,i3,kd)+ty(i1,i2,i3)*tz(i1,i2,i3)*umtt4(i1,i2,i3,kd)+(ry(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,i2,i3)*sy(i1,i2,i3))*umrs4(i1,i2,i3,kd)+(ry(i1,i2,i3)*tz(i1,i2,i3)+rz(i1,i2,i3)*ty(i1,i2,i3))*umrt4(i1,i2,i3,kd)+(sy(i1,i2,i3)*tz(i1,i2,i3)+sz(i1,i2,i3)*ty(i1,i2,i3))*umst4(i1,i2,i3,kd)+ryz43(i1,i2,i3)*umr4(i1,i2,i3,kd)+syz43(i1,i2,i3)*ums4(i1,i2,i3,kd)+tyz43(i1,i2,i3)*umt4(i1,i2,i3,kd)
         umlaplacian43(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2+rz(i1,i2,i3)**2)*umrr4(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**2+sz(i1,i2,i3)**2)*umss4(i1,i2,i3,kd)+(tx(i1,i2,i3)**2+ty(i1,i2,i3)**2+tz(i1,i2,i3)**2)*umtt4(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*sy(i1,i2,i3)+rz(i1,i2,i3)*sz(i1,i2,i3))*umrs4(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*tx(i1,i2,i3)+ ry(i1,i2,i3)*ty(i1,i2,i3)+rz(i1,i2,i3)*tz(i1,i2,i3))*umrt4(i1,i2,i3,kd)+2.*(sx(i1,i2,i3)*tx(i1,i2,i3)+ sy(i1,i2,i3)*ty(i1,i2,i3)+sz(i1,i2,i3)*tz(i1,i2,i3))*umst4(i1,i2,i3,kd)+(rxx43(i1,i2,i3)+ryy43(i1,i2,i3)+rzz43(i1,i2,i3))*umr4(i1,i2,i3,kd)+(sxx43(i1,i2,i3)+syy43(i1,i2,i3)+szz43(i1,i2,i3))*ums4(i1,i2,i3,kd)+(txx43(i1,i2,i3)+tyy43(i1,i2,i3)+tzz43(i1,i2,i3))*umt4(i1,i2,i3,kd)
         !============================================================================================
         ! Define derivatives for a rectangular grid
         !
         !============================================================================================
         umx43r(i1,i2,i3,kd)=(8.*(um(i1+1,i2,i3,kd)-um(i1-1,i2,i3,kd))-(um(i1+2,i2,i3,kd)-um(i1-2,i2,i3,kd)))*h41(0)
         umy43r(i1,i2,i3,kd)=(8.*(um(i1,i2+1,i3,kd)-um(i1,i2-1,i3,kd))-(um(i1,i2+2,i3,kd)-um(i1,i2-2,i3,kd)))*h41(1)
         umz43r(i1,i2,i3,kd)=(8.*(um(i1,i2,i3+1,kd)-um(i1,i2,i3-1,kd))-(um(i1,i2,i3+2,kd)-um(i1,i2,i3-2,kd)))*h41(2)
         umxx43r(i1,i2,i3,kd)=( -30.*um(i1,i2,i3,kd)+16.*(um(i1+1,i2,i3,kd)+um(i1-1,i2,i3,kd))-(um(i1+2,i2,i3,kd)+um(i1-2,i2,i3,kd)) )*h42(0) 
         umyy43r(i1,i2,i3,kd)=( -30.*um(i1,i2,i3,kd)+16.*(um(i1,i2+1,i3,kd)+um(i1,i2-1,i3,kd))-(um(i1,i2+2,i3,kd)+um(i1,i2-2,i3,kd)) )*h42(1) 
         umzz43r(i1,i2,i3,kd)=( -30.*um(i1,i2,i3,kd)+16.*(um(i1,i2,i3+1,kd)+um(i1,i2,i3-1,kd))-(um(i1,i2,i3+2,kd)+um(i1,i2,i3-2,kd)) )*h42(2)
         umxy43r(i1,i2,i3,kd)=( (um(i1+2,i2+2,i3,kd)-um(i1-2,i2+2,i3,kd)- um(i1+2,i2-2,i3,kd)+um(i1-2,i2-2,i3,kd)) +8.*(um(i1-1,i2+2,i3,kd)-um(i1-1,i2-2,i3,kd)-um(i1+1,i2+2,i3,kd)+um(i1+1,i2-2,i3,kd) +um(i1+2,i2-1,i3,kd)-um(i1-2,i2-1,i3,kd)-um(i1+2,i2+1,i3,kd)+um(i1-2,i2+1,i3,kd))+64.*(um(i1+1,i2+1,i3,kd)-um(i1-1,i2+1,i3,kd)- um(i1+1,i2-1,i3,kd)+um(i1-1,i2-1,i3,kd)))*(h41(0)*h41(1))
         umxz43r(i1,i2,i3,kd)=( (um(i1+2,i2,i3+2,kd)-um(i1-2,i2,i3+2,kd)-um(i1+2,i2,i3-2,kd)+um(i1-2,i2,i3-2,kd)) +8.*(um(i1-1,i2,i3+2,kd)-um(i1-1,i2,i3-2,kd)-um(i1+1,i2,i3+2,kd)+um(i1+1,i2,i3-2,kd) +um(i1+2,i2,i3-1,kd)-um(i1-2,i2,i3-1,kd)- um(i1+2,i2,i3+1,kd)+um(i1-2,i2,i3+1,kd)) +64.*(um(i1+1,i2,i3+1,kd)-um(i1-1,i2,i3+1,kd)-um(i1+1,i2,i3-1,kd)+um(i1-1,i2,i3-1,kd)) )*(h41(0)*h41(2))
         umyz43r(i1,i2,i3,kd)=( (um(i1,i2+2,i3+2,kd)-um(i1,i2-2,i3+2,kd)-um(i1,i2+2,i3-2,kd)+um(i1,i2-2,i3-2,kd)) +8.*(um(i1,i2-1,i3+2,kd)-um(i1,i2-1,i3-2,kd)-um(i1,i2+1,i3+2,kd)+um(i1,i2+1,i3-2,kd) +um(i1,i2+2,i3-1,kd)-um(i1,i2-2,i3-1,kd)-um(i1,i2+2,i3+1,kd)+um(i1,i2-2,i3+1,kd)) +64.*(um(i1,i2+1,i3+1,kd)-um(i1,i2-1,i3+1,kd)-um(i1,i2+1,i3-1,kd)+um(i1,i2-1,i3-1,kd)) )*(h41(1)*h41(2))
         umx41r(i1,i2,i3,kd)= umx43r(i1,i2,i3,kd)
         umy41r(i1,i2,i3,kd)= umy43r(i1,i2,i3,kd)
         umz41r(i1,i2,i3,kd)= umz43r(i1,i2,i3,kd)
         umxx41r(i1,i2,i3,kd)= umxx43r(i1,i2,i3,kd)
         umyy41r(i1,i2,i3,kd)= umyy43r(i1,i2,i3,kd)
         umzz41r(i1,i2,i3,kd)= umzz43r(i1,i2,i3,kd)
         umxy41r(i1,i2,i3,kd)= umxy43r(i1,i2,i3,kd)
         umxz41r(i1,i2,i3,kd)= umxz43r(i1,i2,i3,kd)
         umyz41r(i1,i2,i3,kd)= umyz43r(i1,i2,i3,kd)
         umlaplacian41r(i1,i2,i3,kd)=umxx43r(i1,i2,i3,kd)
         umx42r(i1,i2,i3,kd)= umx43r(i1,i2,i3,kd)
         umy42r(i1,i2,i3,kd)= umy43r(i1,i2,i3,kd)
         umz42r(i1,i2,i3,kd)= umz43r(i1,i2,i3,kd)
         umxx42r(i1,i2,i3,kd)= umxx43r(i1,i2,i3,kd)
         umyy42r(i1,i2,i3,kd)= umyy43r(i1,i2,i3,kd)
         umzz42r(i1,i2,i3,kd)= umzz43r(i1,i2,i3,kd)
         umxy42r(i1,i2,i3,kd)= umxy43r(i1,i2,i3,kd)
         umxz42r(i1,i2,i3,kd)= umxz43r(i1,i2,i3,kd)
         umyz42r(i1,i2,i3,kd)= umyz43r(i1,i2,i3,kd)
         umlaplacian42r(i1,i2,i3,kd)=umxx43r(i1,i2,i3,kd)+umyy43r(i1,i2,i3,kd)
         umlaplacian43r(i1,i2,i3,kd)=umxx43r(i1,i2,i3,kd)+umyy43r(i1,i2,i3,kd)+umzz43r(i1,i2,i3,kd)
         !** defineDifferenceOrder2Components1(v,none)
         !** defineDifferenceOrder4Components1(v,none)
         ! *************************************************
         ! *********2nd-order in space and time*************
         ! *************************************************
         ! --- 2D ---
         !...........end   statement functions
         ! write(*,*) 'Inside advIsm...'
         dt    =rpar(0)
         dx(0) =rpar(1)
         dx(1) =rpar(2)
         dx(2) =rpar(3)
         adc   =rpar(4)  ! coefficient of artificial dissipation
         dr(0) =rpar(5)
         dr(1) =rpar(6)
         dr(2) =rpar(7)
         rho   =rpar(8)
         mu    =rpar(9) 
         kx    =rpar(10) 
         ky    =rpar(11) 
         kz    =rpar(12) 
         ep    =rpar(13)
         t     =rpar(14)
         dtOld =rpar(15) ! dt used on the previous time step 
         rpar(20)=0.  ! return the time used for adding dissipation
         dy=dx(1)  ! Are these needed?
         dz=dx(2)
         ! timeForArtificialDissipation=rpar(6) ! return value
         option                        = ipar(0)   ! 0=update solution, 1=add upwind dissipation, 2=compute v.t for MOL schemes
         gridType                      = ipar(1)
         orderOfAccuracy               = ipar(2)
         orderOfAccuracyInTime         = ipar(3)
         addForcing                    = ipar(4)
         twilightZone                  = ipar(5)
         u1c                           = ipar(6)
         u2c                           = ipar(7)
         u3c                           = ipar(8)
         useWhereMask                  = ipar(9)
         timeSteppingMethod            = ipar(10)
         useVariableDissipation        = ipar(11)
         useConservative               = ipar(12)   
         combineDissipationWithAdvance = ipar(13)
         debug                         = ipar(14)
         computeUt                     = ipar(15)  
         materialFormat                = ipar(16)   ! 0=const, 1=piece-wise const, 2=varaiable
         myid                          = ipar(17)
         pc                            = ipar(18)
         upwindSOS                     = ipar(19)
         correction                    = ipar(20)
         ! velocity components : should we pass these in?
         v1c = 0 
         v2c = 1
         v3c = 2
         if( orderOfAccuracy == 2 )then
           useLowerOrderUpwindingOnBoundaries=.false.
         else
           useLowerOrderUpwindingOnBoundaries=.true.
         end if
         cu=  2.     ! coeff. of u(t) in the time-step formula
         cum=-1.     ! coeff. of u(t-dtOld)
         csq=cc**2
         dtsq=dt**2
         cdt=cc*dt
         cs = sqrt(mu/rho)
         cdtSq= (cs*dt)**2
         dtSqByRho = dt**2/rho
         cc = cs 
         rhoi = 1./rho
         ! --- coefficient of linear damping ----   ****TESTING****
         ! This multiplies    dt*( D0t(u) )
         cld = 0.
         cldamp = cld * ( cs/(2.*dt) ) *dt**2   ! for O(1) damping 
         ! cldamp = cldamp * dt   ! make damping O(dt)
         if( option.eq.1 )then 
          useSosupDissipation = 1
         else
          useSosupDissipation = 0
         end if
         preComputeUpwindUt=0 ! What should this be ?
         sosupParameter=1 ! sosupParameter=gamma in sosup scheme  0<= gamma <=1   0=centered scheme
         if( useSosupDissipation.ne.0 )then
           if( .true. )then
             ! -- new way ---
             ! Coefficients in the sosup dissipation from stability analysis 
             if( .true. .or. preComputeUpwindUt.eq.1 )then
                ! We need to reduce the upwind coefficient for stability if we pre-compute uDot: (see CgWave documentation)
               upwindScaleFactor= 1./sqrt(1.*nd)
             else
                upwindScaleFactor= 1.
             end if 
             adSosup2 = -cc*dt*upwindScaleFactor/( 4. )   ! coeff of (h^2 D+D-) 
             adSosup4 = -cc*dt*upwindScaleFactor/( 8. )   ! coeff of (h^2 D+D-)^2 : normal dissipation for order=2
             adSosup = -cc*dt*upwindScaleFactor/( 2**(orderOfAccuracy+1)  )
             ! if( preComputeUpwindUt.eq.1 )then
             !  ! We need to reduce the upwind coefficient for stability if we pre-compute uDot: (see CgWave documentation)
             !  adSosup =adSosup /sqrt(1.*nd)
             ! end if       
           else
             ! Coefficients in the sosup dissipation
             adSosup2 = -cc*dt/( sqrt(1.*nd)*4. )
             adSosup4 =-cc*dt*1./8.
             if( preComputeUpwindUt.eq.1 )then
              ! We need to reduce the upwind coefficient for stability if we pre-compute uDot: (see CgWave documentation)
              adSosup4=adSosup4/sqrt(1.*nd)
             end if     
             if( orderOfAccuracy.eq.2 )then
              adSosup=-cc*dt*1./8.
              if( preComputeUpwindUt.eq.1 )then
                ! We need to reduce the upwind coefficient for stability if we pre-compute uDot: (see CgWave documentation)
                adSosup=adSosup/sqrt(1.*nd)
              end if 
             else if( orderOfAccuracy.eq.4 )then 
               adSosup= -cc*dt*5./288.
               if( .false. )then 
                  adSosup = adSosup*.5 ! ****TEST****
               end if
             else if( orderOfAccuracy.eq.6 )then 
               adSosup=-cc*dt*31./8640.
             else
               stop 1005
             end if
           end if
           if( t<=2*dt )then
             write(*,'("advIsm: adSousp=",1pe12.4)') adSosup
           end if 
           uDotFactor=.5  ! By default uDot is D-zero and so we scale (un-um) by .5 --> .5*(un-um)/(dt)
           ! sosupParameter=gamma in sosup scheme  0<= gamma <=1   0=centered scheme
           adSosup=sosupParameter*adSosup
           ! Note: these next values are only used for rectangular grids. (curvilinear grid values are computed in the loops)
           adxSosup(0)=  uDotFactor*adSosup/dx(0)
           adxSosup(1)=  uDotFactor*adSosup/dx(1)
           adxSosup(2)=  uDotFactor*adSosup/dx(2)    
         end if
        ! c1dtsq=c1*dtsq
        ! c2dtsq=c2*dtsq
        if( t<=2*dt )then
          write(*,'(" advIsm: u1c,u2c,pc=",3i3," t=",1pe12.3)') u1c,u2c,pc,t
          write(*,'("## debug=",i3," addForcing=",i4," twilightZone=",i4," option=",i2," useWhereMask=",i2)') debug,addForcing,twilightZone,option,useWhereMask
          if( option==3 )then
            write(*,'("  MethodOfLines update: correction=",i3)') correction
          end if
        end if
        if( dtOld.le.0 )then
          write(*,'(" advIsm:ERROR : dtOld<=0 ")')
          stop 8167
        end if
        if( dt.ne.dtOld )then
          write(*,'(" advIsm:INFO: dt=",e12.4," <> dtOld=",e12.4," diff=",e9.2)') dt,dtOld,dt-dtOld
          if( orderOfAccuracy.ne.2 )then
            write(*,'(" advIsm:ERROR: variable dt not implemented for orderOfAccuracy=",i4)') orderOfAccuracy
            ! '
            stop 8168
          end if
          ! adjust the coefficients for a variable time step : this is locally second order accurate
          cu= 1.+dt/dtOld     ! coeff. of u(t) in the time-step formula
          cum=-dt/dtOld       ! coeff. of u(t-dtOld)
          ! c1dtsq=c1*dt*(dt+dtOld)*.5
          ! c2dtsq=c2*dt*(dt+dtOld)*.5
        end if
        dt4by12=dtsq*dtsq/12.
        call ovtime( time0 )
          ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          ! ++++++++++++++++++++ DEFINE OPERATORS +++++++++++++++++++++++
          ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            ! ***** CURVILINEAR *******
         ! write(*,*) 'Inside advIsm: order=2 dim=2'
         if( option==0 )then
           ! ------------ UPDATE SOLUTION FOR Modified Equation SOS Scheme ----------------
           if( correction==0 )then
             if( twilightZone==1 )then
                   do i3=n3a,n3b
                   do i2=n2a,n2b
                   do i1=n1a,n1b
                   un(i1,i2,i3,u1c) = 2.*u(i1,i2,i3,u1c) - um(i1,i2,i3,u1c) - dtSqByRho*ux22(i1,i2,i3,pc) + cdtSq*( ulaplacian22(i1,i2,i3,u1c) )
                   un(i1,i2,i3,u2c) = 2.*u(i1,i2,i3,u2c) - um(i1,i2,i3,u2c) - dtSqByRho*uy22(i1,i2,i3,pc) + cdtSq*( ulaplacian22(i1,i2,i3,u2c) )
                       call ogDeriv(ep,2,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,u1c,uett)
                       call ogDeriv(ep,0,2,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,u1c,uexx)
                       call ogDeriv(ep,0,0,2,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,u1c,ueyy)
                       call ogDeriv(ep,2,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,u2c,vett)
                       call ogDeriv(ep,0,2,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,u2c,vexx)
                       call ogDeriv(ep,0,0,2,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,u2c,veyy) 
                       call ogDeriv(ep,0,1,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,pc,pex)                           
                       call ogDeriv(ep,0,0,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,pc,pey)
                       un(i1,i2,i3,u1c) =  un(i1,i2,i3,u1c) + dtSq*uett + dtSqByRho*pex - cdtSq*( uexx + ueyy )
                       un(i1,i2,i3,u2c) =  un(i1,i2,i3,u2c) + dtSq*vett + dtSqByRho*pey - cdtSq*( vexx + veyy )
                  ! un(i1,i2,i3,u1c)=sm22ru(i1,i2,i3) + dtsq*f(i1,i2,i3,u1c)
                   ! un(i1,i2,i3,u2c)=sm22rv(i1,i2,i3) + dtsq*f(i1,i2,i3,u2c) 
                   ! FILL IN p with EXACT Solution for testing
                   ! call ogDeriv(ep,0,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t+dt,pc,pe)
                   ! un(i1,i2,i3,pc)=pe
                   end do
                   end do
                   end do
             else
                   do i3=n3a,n3b
                   do i2=n2a,n2b
                   do i1=n1a,n1b
                   un(i1,i2,i3,u1c) = 2.*u(i1,i2,i3,u1c) - um(i1,i2,i3,u1c) - dtSqByRho*ux22(i1,i2,i3,pc) + cdtSq*( ulaplacian22(i1,i2,i3,u1c) )
                   un(i1,i2,i3,u2c) = 2.*u(i1,i2,i3,u2c) - um(i1,i2,i3,u2c) - dtSqByRho*uy22(i1,i2,i3,pc) + cdtSq*( ulaplacian22(i1,i2,i3,u2c) )
                  ! un(i1,i2,i3,u1c)=sm22ru(i1,i2,i3) + dtsq*f(i1,i2,i3,u1c)
                   ! un(i1,i2,i3,u2c)=sm22rv(i1,i2,i3) + dtsq*f(i1,i2,i3,u2c) 
                   ! FILL IN p with EXACT Solution for testing
                   ! call ogDeriv(ep,0,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t+dt,pc,pe)
                   ! un(i1,i2,i3,pc)=pe
                   end do
                   end do
                   end do
             end if
           else
             ! ----- Modified equation corrector stage ----
             write(*,*) 'advOptISM: Corrector called!'
             if( twilightZone==1 )then
                   do i3=n3a,n3b
                   do i2=n2a,n2b
                   do i1=n1a,n1b
                   uTemp(i1,i2,i3,u1c) = 2.*u(i1,i2,i3,u1c) - um(i1,i2,i3,u1c) - (.5*dtSqByRho)*(unx22(i1,i2,i3,pc) + umx22(i1,i2,i3,pc)) + (.5*cdtSq)*( unlaplacian22(i1,i2,i3,u1c) + umlaplacian22(i1,i2,i3,u1c) )
                   uTemp(i1,i2,i3,u2c) = 2.*u(i1,i2,i3,u2c) - um(i1,i2,i3,u2c) - (.5*dtSqByRho)*(uny22(i1,i2,i3,pc) + umy22(i1,i2,i3,pc)) + (.5*cdtSq)*( unlaplacian22(i1,i2,i3,u2c) + umlaplacian22(i1,i2,i3,u2c) )
                       call ogDeriv(ep,2,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,u1c,uett)
                       call ogDeriv(ep,0,2,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,u1c,uexx)
                       call ogDeriv(ep,0,0,2,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,u1c,ueyy)
                       call ogDeriv(ep,2,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,u2c,vett)
                       call ogDeriv(ep,0,2,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,u2c,vexx)
                       call ogDeriv(ep,0,0,2,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,u2c,veyy) 
                       call ogDeriv(ep,0,1,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,pc,pex)                           
                       call ogDeriv(ep,0,0,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,pc,pey)
                       uTemp(i1,i2,i3,u1c) =  uTemp(i1,i2,i3,u1c) + dtSq*uett + dtSqByRho*pex - cdtSq*( uexx + ueyy )
                       uTemp(i1,i2,i3,u2c) =  uTemp(i1,i2,i3,u2c) + dtSq*vett + dtSqByRho*pey - cdtSq*( vexx + veyy )
                  ! un(i1,i2,i3,u1c)=sm22ru(i1,i2,i3) + dtsq*f(i1,i2,i3,u1c)
                   ! un(i1,i2,i3,u2c)=sm22rv(i1,i2,i3) + dtsq*f(i1,i2,i3,u2c) 
                   ! FILL IN p with EXACT Solution for testing
                   ! call ogDeriv(ep,0,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t+dt,pc,pe)
                   ! un(i1,i2,i3,pc)=pe
                   end do
                   end do
                   end do
                   do i3=n3a,n3b
                   do i2=n2a,n2b
                   do i1=n1a,n1b
                   un(i1,i2,i3,u1c) = uTemp(i1,i2,i3,u1c)
                   un(i1,i2,i3,u2c) = uTemp(i1,i2,i3,u2c)
                   end do
                   end do
                   end do
             else
                   do i3=n3a,n3b
                   do i2=n2a,n2b
                   do i1=n1a,n1b
                   uTemp(i1,i2,i3,u1c) = 2.*u(i1,i2,i3,u1c) - um(i1,i2,i3,u1c) - (.5*dtSqByRho)*(unx22(i1,i2,i3,pc) + umx22(i1,i2,i3,pc)) + (.5*cdtSq)*( unlaplacian22(i1,i2,i3,u1c) + umlaplacian22(i1,i2,i3,u1c) )
                   uTemp(i1,i2,i3,u2c) = 2.*u(i1,i2,i3,u2c) - um(i1,i2,i3,u2c) - (.5*dtSqByRho)*(uny22(i1,i2,i3,pc) + umy22(i1,i2,i3,pc)) + (.5*cdtSq)*( unlaplacian22(i1,i2,i3,u2c) + umlaplacian22(i1,i2,i3,u2c) )
                  ! un(i1,i2,i3,u1c)=sm22ru(i1,i2,i3) + dtsq*f(i1,i2,i3,u1c)
                   ! un(i1,i2,i3,u2c)=sm22rv(i1,i2,i3) + dtsq*f(i1,i2,i3,u2c) 
                   ! FILL IN p with EXACT Solution for testing
                   ! call ogDeriv(ep,0,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t+dt,pc,pe)
                   ! un(i1,i2,i3,pc)=pe
                   end do
                   end do
                   end do
                   do i3=n3a,n3b
                   do i2=n2a,n2b
                   do i1=n1a,n1b
                   un(i1,i2,i3,u1c) = uTemp(i1,i2,i3,u1c)
                   un(i1,i2,i3,u2c) = uTemp(i1,i2,i3,u2c)
                   end do
                   end do
                   end do
             end if      
           end if
         elseif( option==1 )then 
           !   ------------ Add upwind dissipation ---
           if( useSosupDissipation.eq.1 .and. preComputeUpwindUt.eq.1 )then
             if( option.ne.1 )then
               write(*,'("advIsmOpt:ERROR: useSosupDissipation.eq.1 BUT option.ne.1")')
               stop 6663
             end if
                ! precompute "uDot" = dt*du/dt used in the dissipation and store in v 
                ! we need uDot at enough ghost points for the dissipation operator 
                if( debug.gt.3 .and. t.le.3.*dt )then
                  write(*,'(" advWave: add UPWIND DISSIPATION: Evaluate v= uDot")') 
                end if
                numGhost=orderOfAccuracy/2
                numGhost=numGhost+1
                m1a=n1a-numGhost
                m1b=n1b+numGhost
                m2a=n2a-numGhost
                m2b=n2b+numGhost
                if( nd.eq.2 )then
                 m3a=n3a
                 m3b=n3b
                else
                  m3a=n3a-numGhost
                  m3b=n3b+numGhost
                end if
                ! write(*,'(" numGhost=",i2," m1a,m1b,m2a,m2b=",4i4)') numGhost,m1a,m1b,m2a,m2b
                ! We need v at ghost outside interpolation points -- do not use mask here
                  do i3=m3a,m3b
                  do i2=m2a,m2b
                  do i1=m1a,m1b
                  v(i1,i2,i3,0)=un(i1,i2,i3,0)-um(i1,i2,i3,0)
                  ! v(i1,i2,i3,0)=u(i1,i2,i3,0)-un(i1,i2,i3,0)
                  end do
                  end do
                  end do
           end if 
           ! ---- add upwind dissipation -----
           ! preComputeUpwindUt : true=precompute Ut in upwind dissipation,  (uses v=uDot computed above)
           !                      false=compute Ut inline in Gauss-Seidel fashion 
           if( preComputeUpwindUt.eq.1 )then
             ! precompute Ut in upwind dissipation,  (uses v=uDot computed above)
               ! Note: these next values are only used for rectangular grids. (curvilinear grid values are computed in the loops)
               adxSosup(0)=  uDotFactor*adSosup/dx(0)
               adxSosup(1)=  uDotFactor*adSosup/dx(1)
               adxSosup(2)=  uDotFactor*adSosup/dx(2) 
               if( debug.gt.3 .or. t.lt. 2.*dt )then
                 write(*,'("addUpwindDiss: UPWIND DISS using u.t=v dim=2 order=2 grid=curvilinear... t=",e10.2)') t
                 write(*,'(" adxSosup=",3e12.4)') adxSosup(0), adxSosup(1),adxSosup(2)
                 write(*,'(" n1a,n1b,n2a,n2b,n3a,n3b=",6i3)')  n1a,n1b,n2a,n2b,n3a,n3b
               end if
               m=0 ! component number 
               ec = 0 ! component number
               ! We use lower-order upwind on boundaries for traction BCs so skip the boundaries with the high-order upwind
               ! Upwinding is not needed on displacement BCs where u and v are set so skip the boundaries too
               nn1a=n1a; nn1b=n1b; 
               nn2a=n2a; nn2b=n2b; 
               nn3a=n3a; nn3b=n3b;
               if( useLowerOrderUpwindingOnBoundaries )then
                 if( boundaryCondition(0,0)==tractionBC .or. boundaryCondition(0,0)==displacementBC )then
                   nn1a=nn1a+1
                 end if
                 if( boundaryCondition(1,0)==tractionBC .or. boundaryCondition(1,0)==displacementBC )then
                   nn1b=nn1b-1
                 end if 
                 if( boundaryCondition(0,1)==tractionBC .or. boundaryCondition(0,1)==displacementBC )then
                   nn2a=nn2a+1
                 end if
                 if( boundaryCondition(1,1)==tractionBC .or. boundaryCondition(1,1)==displacementBC )then
                   nn2b=nn2b-1
                 end if 
                 if( nd==3 )then
                   if( boundaryCondition(0,2)==tractionBC .or. boundaryCondition(0,2)==displacementBC )then
                     nn3a=nn3a+1
                   end if
                   if( boundaryCondition(1,2)==tractionBC .or. boundaryCondition(1,2)==displacementBC )then
                     nn3b=nn3b-1
                   end if    
                 end if  
               end if
                 do i3=nn3a,nn3b
                 do i2=nn2a,nn2b
                 do i1=nn1a,nn1b
                  ! --- SECOND 2 ---
                    ! --- TWO DIMENSIONS ---
                       if( mask(i1,i2,i3)>0 )then
                          do dir=0,1
                            ! diss-coeff ~= 1/(change in x along direction r(dir) )
                            ! Assuming a nearly orthogonal grid gives ||dx|| = || grad_x(r_i) || / dr_i 
                            adxSosup(dir) = adSosup*uDotFactor*sqrt( rsxy(i1,i2,i3,dir,0)**2 + rsxy(i1,i2,i3,dir,1)**2 )/dr(dir) 
                          end do
                         ! if( mask(i1-1,i2-1,i3)>0 .and. mask(i1  ,i2-1,i3)>0 .and. mask(i1+1,i2-1,i3)>0 .and. !     mask(i1-1,i2  ,i3)>0 .and. mask(i1  ,i2  ,i3)>0 .and. mask(i1+1,i2  ,i3)>0 .and. !     mask(i1-1,i2+1,i3)>0 .and. mask(i1  ,i2+1,i3)>0 .and. mask(i1+1,i2+1,i3)>0  )then
                         ! this point is not next to an interp point
                         uTemp(i1,i2,i3,u1c)=un(i1,i2,i3,u1c)+(+6.*v(i1,i2,i3,u1c)-4.*(v(i1+1,i2,i3,u1c)+v(i1-1,i2,i3,u1c))+(v(i1+2,i2,i3,u1c)+v(i1-2,i2,i3,u1c)))*adxSosup(0)+(+6.*v(i1,i2,i3,u1c)-4.*(v(i1,i2+1,i3,u1c)+v(i1,i2-1,i3,u1c))+(v(i1,i2+2,i3,u1c)+v(i1,i2-2,i3,u1c)))*adxSosup(1) - cldamp*v(i1,i2,i3,u1c)
                         uTemp(i1,i2,i3,u2c)=un(i1,i2,i3,u2c)+(+6.*v(i1,i2,i3,u2c)-4.*(v(i1+1,i2,i3,u2c)+v(i1-1,i2,i3,u2c))+(v(i1+2,i2,i3,u2c)+v(i1-2,i2,i3,u2c)))*adxSosup(0)+(+6.*v(i1,i2,i3,u2c)-4.*(v(i1,i2+1,i3,u2c)+v(i1,i2-1,i3,u2c))+(v(i1,i2+2,i3,u2c)+v(i1,i2-2,i3,u2c)))*adxSosup(1) - cldamp*v(i1,i2,i3,u2c)
                       end if
                     ! else if( mask(i1,i2,i3)>0 )then
                     !   ! this point must be next to an interp point
                     !   un(i1,i2,i3,u1c)=un(i1,i2,i3,u1c)+sosupDiss2d2(v,i1,i2,i3,u1c)
                     !   un(i1,i2,i3,u2c)=un(i1,i2,i3,u2c)+sosupDiss2d2(v,i1,i2,i3,u2c)          
                     ! end if
                 end do
                 end do
                 end do
               ! endLoopsMask()
               if( useLowerOrderUpwindingOnBoundaries )then
                 ! ------ Use lower-order dissipation on the boundary ------
                   adSosup = adSosup2  ! use this coeff
                 ! Note: these next values are only used for rectangular grids. (curvilinear grid values are computed in the loops)
                 adxSosup(0)=  uDotFactor*adSosup/dx(0)
                 adxSosup(1)=  uDotFactor*adSosup/dx(1)
                 adxSosup(2)=  uDotFactor*adSosup/dx(2) 
                  ! extra1a=extra
                  ! extra1b=extra
                  ! extra2a=extra
                  ! extra2b=extra
                  ! if( nd.eq.3 )then
                  !   extra3a=extra
                  !   extra3b=extra
                  ! else
                  !   extra3a=0
                  !   extra3b=0
                  ! end if
                  ! if( boundaryCondition(0,0).lt.0 )then
                  !   extra1a=max(0,extra1a) ! over-ride extra=-1 : assign ends in periodic directions (or internal parallel boundaries)
                  ! else if( boundaryCondition(0,0).eq.0 )then
                  !   extra1a=numberOfGhostPoints  ! include interpolation points since we assign ghost points outside these
                  ! end if
                  ! ! **NOTE** the bc on the right may be negative even it is not on the left (for parallel)
                  ! if( boundaryCondition(1,0).lt.0 )then
                  !   extra1b=max(0,extra1b) ! over-ride extra=-1 : assign ends in periodic directions
                  ! else if( boundaryCondition(1,0).eq.0 )then
                  !   extra1b=numberOfGhostPoints
                  ! end if
                  ! if( boundaryCondition(0,1).lt.0 )then
                  !   extra2a=max(0,extra2a) ! over-ride extra=-1 : assign ends in periodic directions (or internal parallel boundaries)
                  ! else if( boundaryCondition(0,1).eq.0 )then
                  !   extra2a=numberOfGhostPoints  ! include interpolation points since we assign ghost points outside these
                  ! end if
                  ! ! **NOTE** the bc on the right may be negative even it is not on the left (for parallel)
                  ! if( boundaryCondition(1,1).lt.0 )then
                  !   extra2b=max(0,extra2b) ! over-ride extra=-1 : assign ends in periodic directions
                  ! else if( boundaryCondition(1,1).eq.0 )then
                  !   extra2b=numberOfGhostPoints
                  ! end if
                  ! if(  nd.eq.3 )then
                  !  if( boundaryCondition(0,2).lt.0 )then
                  !    extra3a=max(0,extra3a) ! over-ride extra=-1 : assign ends in periodic directions (or internal parallel boundaries)
                  !  else if( boundaryCondition(0,2).eq.0 )then
                  !    extra3a=numberOfGhostPoints  ! include interpolation points since we assign ghost points outside these
                  !  end if
                  !  ! **NOTE** the bc on the right may be negative even it is not on the left (for parallel)
                  !  if( boundaryCondition(1,2).lt.0 )then
                  !    extra3b=max(0,extra3b) ! over-ride extra=-1 : assign ends in periodic directions
                  !  else if( boundaryCondition(1,2).eq.0 )then
                  !    extra3b=numberOfGhostPoints
                  !  end if
                  ! end if
                  do axis=0,nd-1
                  do side=0,1
                    if( boundaryCondition(side,axis).gt.0 )then
                      ! write(*,'(" bcOpt: side,axis,bc=",3i2)') side,axis,boundaryCondition(side,axis)
                      nn1a=gridIndexRange(0,0)
                      nn1b=gridIndexRange(1,0)
                      nn2a=gridIndexRange(0,1)
                      nn2b=gridIndexRange(1,1)
                      nn3a=gridIndexRange(0,2)
                      nn3b=gridIndexRange(1,2)
                      if( axis.eq.0 )then
                        nn1a=gridIndexRange(side,axis)
                        nn1b=gridIndexRange(side,axis)
                      else if( axis.eq.1 )then
                        nn2a=gridIndexRange(side,axis)
                        nn2b=gridIndexRange(side,axis)
                      else
                        nn3a=gridIndexRange(side,axis)
                        nn3b=gridIndexRange(side,axis)
                      end if
                      ! nn1a=gridIndexRange(0,0)-extra1a
                      ! nn1b=gridIndexRange(1,0)+extra1b
                      ! nn2a=gridIndexRange(0,1)-extra2a
                      ! nn2b=gridIndexRange(1,1)+extra2b
                      ! nn3a=gridIndexRange(0,2)-extra3a
                      ! nn3b=gridIndexRange(1,2)+extra3b
                      ! if( axis.eq.0 )then
                      !   nn1a=gridIndexRange(side,axis)
                      !   nn1b=gridIndexRange(side,axis)
                      ! else if( axis.eq.1 )then
                      !   nn2a=gridIndexRange(side,axis)
                      !   nn2b=gridIndexRange(side,axis)
                      ! else
                      !   nn3a=gridIndexRange(side,axis)
                      !   nn3b=gridIndexRange(side,axis)
                      ! end if
                      is=1-2*side
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
                      i3=n3a
                      if( debug.gt.7 )then
                        write(*,'(" advOptIsm: grid,side,axis=",3i3,", loop bounds: nn1a,nn1b,nn2a,nn2b,nn3a,nn3b=",6i3)') grid,side,axis,nn1a,nn1b,nn2a,nn2b,nn3a,nn3b
                      end if
                    end if ! if bc>0 
                    ! ! On interfaces we should use the bcf array values even for TZ since then
                    ! ! we get a coupling at the interface: 
                    ! !   bcf = n.sigma(fluid) + [ n.sigma_e(solid) - n.sigma_e(fluid) ]
                    ! if( interfaceType(side,axis,grid).eq.noInterface )then
                    !   assignTwilightZone=twilightZone
                    ! else
                    !   assignTwilightZone=0  ! this will turn off the use of TZ
                    ! end if
                  ! if( boundaryCondition(side,axis)==tractionBC .or. boundaryCondition(side,axis)==displacementBC )then
                  if( boundaryCondition(side,axis)==tractionBC )then
                    if( t.le.2*dt )then
                      write(*,'("Use lower order upwind dissipation on boundaries: side,axis=",2i2," adSosup=",1pe12.4)') side,axis,adSosup
                      write(*,'("  Boundary: nn1a,nn1b,nn2a,nn2b=",4i5)') nn1a,nn1b,nn2a,nn2b
                    end if      
                      do i3=nn3a,nn3b
                      do i2=nn2a,nn2b
                      do i1=nn1a,nn1b
                        ! --- SECOND 2 ---
                          ! --- TWO DIMENSIONS ---
                             do dir=0,1
                               ! diss-coeff ~= 1/(change in x along direction r(dir) )
                               ! Assuming a nearly orthogonal grid gives ||dx|| = || grad_x(r_i) || / dr_i 
                               adxSosup(dir) = adSosup*uDotFactor*sqrt( rsxy(i1,i2,i3,dir,0)**2 + rsxy(i1,i2,i3,dir,1)**2 )/dr(dir) 
                             end do
                            uTemp(i1,i2,i3,u1c)=un(i1,i2,i3,u1c)+(+2.*v(i1,i2,i3,u1c)-(v(i1+1,i2,i3,u1c)+v(i1-1,i2,i3,u1c)))*adxSosup(0)+(+2.*v(i1,i2,i3,u1c)-(v(i1,i2+1,i3,u1c)+v(i1,i2-1,i3,u1c)))*adxSosup(1)
                            uTemp(i1,i2,i3,u2c)=un(i1,i2,i3,u2c)+(+2.*v(i1,i2,i3,u2c)-(v(i1+1,i2,i3,u2c)+v(i1-1,i2,i3,u2c)))*adxSosup(0)+(+2.*v(i1,i2,i3,u2c)-(v(i1,i2+1,i3,u2c)+v(i1,i2-1,i3,u2c)))*adxSosup(1)
                      end do
                      end do
                      end do
                   end if ! end if bc==
                  end do ! end side
                  end do ! end axis
               end if
               ! -- Copy temp solution back to un ----
                 do i3=n3a,n3b
                 do i2=n2a,n2b
                 do i1=n1a,n1b
                 un(i1,i2,i3,u1c)=uTemp(i1,i2,i3,u1c)
                 un(i1,i2,i3,u2c)=uTemp(i1,i2,i3,u2c)
                 end do
                 end do
                 end do
           else
             ! compute Ut inline in Gauss-Seidel fashion (this is more stable)
               ! Note: these next values are only used for rectangular grids. (curvilinear grid values are computed in the loops)
               adxSosup(0)=  uDotFactor*adSosup/dx(0)
               adxSosup(1)=  uDotFactor*adSosup/dx(1)
               adxSosup(2)=  uDotFactor*adSosup/dx(2) 
               if( debug.gt.3 .or. t.lt. 2.*dt )then
                 write(*,'("addUpwindDiss: UPWIND DISS using u.t=Dztu dim=2 order=2 grid=curvilinear... t=",e10.2)') t
                 write(*,'(" adxSosup=",3e12.4)') adxSosup(0), adxSosup(1),adxSosup(2)
                 write(*,'(" n1a,n1b,n2a,n2b,n3a,n3b=",6i3)')  n1a,n1b,n2a,n2b,n3a,n3b
               end if
               m=0 ! component number 
               ec = 0 ! component number
               ! We use lower-order upwind on boundaries for traction BCs so skip the boundaries with the high-order upwind
               ! Upwinding is not needed on displacement BCs where u and v are set so skip the boundaries too
               nn1a=n1a; nn1b=n1b; 
               nn2a=n2a; nn2b=n2b; 
               nn3a=n3a; nn3b=n3b;
               if( useLowerOrderUpwindingOnBoundaries )then
                 if( boundaryCondition(0,0)==tractionBC .or. boundaryCondition(0,0)==displacementBC )then
                   nn1a=nn1a+1
                 end if
                 if( boundaryCondition(1,0)==tractionBC .or. boundaryCondition(1,0)==displacementBC )then
                   nn1b=nn1b-1
                 end if 
                 if( boundaryCondition(0,1)==tractionBC .or. boundaryCondition(0,1)==displacementBC )then
                   nn2a=nn2a+1
                 end if
                 if( boundaryCondition(1,1)==tractionBC .or. boundaryCondition(1,1)==displacementBC )then
                   nn2b=nn2b-1
                 end if 
                 if( nd==3 )then
                   if( boundaryCondition(0,2)==tractionBC .or. boundaryCondition(0,2)==displacementBC )then
                     nn3a=nn3a+1
                   end if
                   if( boundaryCondition(1,2)==tractionBC .or. boundaryCondition(1,2)==displacementBC )then
                     nn3b=nn3b-1
                   end if    
                 end if  
               end if
                 do i3=nn3a,nn3b
                 do i2=nn2a,nn2b
                 do i1=nn1a,nn1b
                  ! --- SECOND 2 ---
                    ! --- TWO DIMENSIONS ---
                       if( mask(i1,i2,i3)>0 )then
                          do dir=0,1
                            ! diss-coeff ~= 1/(change in x along direction r(dir) )
                            ! Assuming a nearly orthogonal grid gives ||dx|| = || grad_x(r_i) || / dr_i 
                            adxSosup(dir) = adSosup*uDotFactor*sqrt( rsxy(i1,i2,i3,dir,0)**2 + rsxy(i1,i2,i3,dir,1)**2 )/dr(dir) 
                          end do
                         ! if( mask(i1-1,i2-1,i3)>0 .and. mask(i1  ,i2-1,i3)>0 .and. mask(i1+1,i2-1,i3)>0 .and. !     mask(i1-1,i2  ,i3)>0 .and. mask(i1  ,i2  ,i3)>0 .and. mask(i1+1,i2  ,i3)>0 .and. !     mask(i1-1,i2+1,i3)>0 .and. mask(i1  ,i2+1,i3)>0 .and. mask(i1+1,i2+1,i3)>0  )then
                         ! this point is not next to an interp point
                         uTemp(i1,i2,i3,u1c)=un(i1,i2,i3,u1c)+(+6.*Dztu(i1,i2,i3,u1c)-4.*(Dztu(i1+1,i2,i3,u1c)+Dztu(i1-1,i2,i3,u1c))+(Dztu(i1+2,i2,i3,u1c)+Dztu(i1-2,i2,i3,u1c)))*adxSosup(0)+(+6.*Dztu(i1,i2,i3,u1c)-4.*(Dztu(i1,i2+1,i3,u1c)+Dztu(i1,i2-1,i3,u1c))+(Dztu(i1,i2+2,i3,u1c)+Dztu(i1,i2-2,i3,u1c)))*adxSosup(1) - cldamp*Dztu(i1,i2,i3,u1c)
                         uTemp(i1,i2,i3,u2c)=un(i1,i2,i3,u2c)+(+6.*Dztu(i1,i2,i3,u2c)-4.*(Dztu(i1+1,i2,i3,u2c)+Dztu(i1-1,i2,i3,u2c))+(Dztu(i1+2,i2,i3,u2c)+Dztu(i1-2,i2,i3,u2c)))*adxSosup(0)+(+6.*Dztu(i1,i2,i3,u2c)-4.*(Dztu(i1,i2+1,i3,u2c)+Dztu(i1,i2-1,i3,u2c))+(Dztu(i1,i2+2,i3,u2c)+Dztu(i1,i2-2,i3,u2c)))*adxSosup(1) - cldamp*Dztu(i1,i2,i3,u2c)
                       end if
                     ! else if( mask(i1,i2,i3)>0 )then
                     !   ! this point must be next to an interp point
                     !   un(i1,i2,i3,u1c)=un(i1,i2,i3,u1c)+sosupDiss2d2(Dztu,i1,i2,i3,u1c)
                     !   un(i1,i2,i3,u2c)=un(i1,i2,i3,u2c)+sosupDiss2d2(Dztu,i1,i2,i3,u2c)          
                     ! end if
                 end do
                 end do
                 end do
               ! endLoopsMask()
               if( useLowerOrderUpwindingOnBoundaries )then
                 ! ------ Use lower-order dissipation on the boundary ------
                   adSosup = adSosup2  ! use this coeff
                 ! Note: these next values are only used for rectangular grids. (curvilinear grid values are computed in the loops)
                 adxSosup(0)=  uDotFactor*adSosup/dx(0)
                 adxSosup(1)=  uDotFactor*adSosup/dx(1)
                 adxSosup(2)=  uDotFactor*adSosup/dx(2) 
                  ! extra1a=extra
                  ! extra1b=extra
                  ! extra2a=extra
                  ! extra2b=extra
                  ! if( nd.eq.3 )then
                  !   extra3a=extra
                  !   extra3b=extra
                  ! else
                  !   extra3a=0
                  !   extra3b=0
                  ! end if
                  ! if( boundaryCondition(0,0).lt.0 )then
                  !   extra1a=max(0,extra1a) ! over-ride extra=-1 : assign ends in periodic directions (or internal parallel boundaries)
                  ! else if( boundaryCondition(0,0).eq.0 )then
                  !   extra1a=numberOfGhostPoints  ! include interpolation points since we assign ghost points outside these
                  ! end if
                  ! ! **NOTE** the bc on the right may be negative even it is not on the left (for parallel)
                  ! if( boundaryCondition(1,0).lt.0 )then
                  !   extra1b=max(0,extra1b) ! over-ride extra=-1 : assign ends in periodic directions
                  ! else if( boundaryCondition(1,0).eq.0 )then
                  !   extra1b=numberOfGhostPoints
                  ! end if
                  ! if( boundaryCondition(0,1).lt.0 )then
                  !   extra2a=max(0,extra2a) ! over-ride extra=-1 : assign ends in periodic directions (or internal parallel boundaries)
                  ! else if( boundaryCondition(0,1).eq.0 )then
                  !   extra2a=numberOfGhostPoints  ! include interpolation points since we assign ghost points outside these
                  ! end if
                  ! ! **NOTE** the bc on the right may be negative even it is not on the left (for parallel)
                  ! if( boundaryCondition(1,1).lt.0 )then
                  !   extra2b=max(0,extra2b) ! over-ride extra=-1 : assign ends in periodic directions
                  ! else if( boundaryCondition(1,1).eq.0 )then
                  !   extra2b=numberOfGhostPoints
                  ! end if
                  ! if(  nd.eq.3 )then
                  !  if( boundaryCondition(0,2).lt.0 )then
                  !    extra3a=max(0,extra3a) ! over-ride extra=-1 : assign ends in periodic directions (or internal parallel boundaries)
                  !  else if( boundaryCondition(0,2).eq.0 )then
                  !    extra3a=numberOfGhostPoints  ! include interpolation points since we assign ghost points outside these
                  !  end if
                  !  ! **NOTE** the bc on the right may be negative even it is not on the left (for parallel)
                  !  if( boundaryCondition(1,2).lt.0 )then
                  !    extra3b=max(0,extra3b) ! over-ride extra=-1 : assign ends in periodic directions
                  !  else if( boundaryCondition(1,2).eq.0 )then
                  !    extra3b=numberOfGhostPoints
                  !  end if
                  ! end if
                  do axis=0,nd-1
                  do side=0,1
                    if( boundaryCondition(side,axis).gt.0 )then
                      ! write(*,'(" bcOpt: side,axis,bc=",3i2)') side,axis,boundaryCondition(side,axis)
                      nn1a=gridIndexRange(0,0)
                      nn1b=gridIndexRange(1,0)
                      nn2a=gridIndexRange(0,1)
                      nn2b=gridIndexRange(1,1)
                      nn3a=gridIndexRange(0,2)
                      nn3b=gridIndexRange(1,2)
                      if( axis.eq.0 )then
                        nn1a=gridIndexRange(side,axis)
                        nn1b=gridIndexRange(side,axis)
                      else if( axis.eq.1 )then
                        nn2a=gridIndexRange(side,axis)
                        nn2b=gridIndexRange(side,axis)
                      else
                        nn3a=gridIndexRange(side,axis)
                        nn3b=gridIndexRange(side,axis)
                      end if
                      ! nn1a=gridIndexRange(0,0)-extra1a
                      ! nn1b=gridIndexRange(1,0)+extra1b
                      ! nn2a=gridIndexRange(0,1)-extra2a
                      ! nn2b=gridIndexRange(1,1)+extra2b
                      ! nn3a=gridIndexRange(0,2)-extra3a
                      ! nn3b=gridIndexRange(1,2)+extra3b
                      ! if( axis.eq.0 )then
                      !   nn1a=gridIndexRange(side,axis)
                      !   nn1b=gridIndexRange(side,axis)
                      ! else if( axis.eq.1 )then
                      !   nn2a=gridIndexRange(side,axis)
                      !   nn2b=gridIndexRange(side,axis)
                      ! else
                      !   nn3a=gridIndexRange(side,axis)
                      !   nn3b=gridIndexRange(side,axis)
                      ! end if
                      is=1-2*side
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
                      i3=n3a
                      if( debug.gt.7 )then
                        write(*,'(" advOptIsm: grid,side,axis=",3i3,", loop bounds: nn1a,nn1b,nn2a,nn2b,nn3a,nn3b=",6i3)') grid,side,axis,nn1a,nn1b,nn2a,nn2b,nn3a,nn3b
                      end if
                    end if ! if bc>0 
                    ! ! On interfaces we should use the bcf array values even for TZ since then
                    ! ! we get a coupling at the interface: 
                    ! !   bcf = n.sigma(fluid) + [ n.sigma_e(solid) - n.sigma_e(fluid) ]
                    ! if( interfaceType(side,axis,grid).eq.noInterface )then
                    !   assignTwilightZone=twilightZone
                    ! else
                    !   assignTwilightZone=0  ! this will turn off the use of TZ
                    ! end if
                  ! if( boundaryCondition(side,axis)==tractionBC .or. boundaryCondition(side,axis)==displacementBC )then
                  if( boundaryCondition(side,axis)==tractionBC )then
                    if( t.le.2*dt )then
                      write(*,'("Use lower order upwind dissipation on boundaries: side,axis=",2i2," adSosup=",1pe12.4)') side,axis,adSosup
                      write(*,'("  Boundary: nn1a,nn1b,nn2a,nn2b=",4i5)') nn1a,nn1b,nn2a,nn2b
                    end if      
                      do i3=nn3a,nn3b
                      do i2=nn2a,nn2b
                      do i1=nn1a,nn1b
                        ! --- SECOND 2 ---
                          ! --- TWO DIMENSIONS ---
                             do dir=0,1
                               ! diss-coeff ~= 1/(change in x along direction r(dir) )
                               ! Assuming a nearly orthogonal grid gives ||dx|| = || grad_x(r_i) || / dr_i 
                               adxSosup(dir) = adSosup*uDotFactor*sqrt( rsxy(i1,i2,i3,dir,0)**2 + rsxy(i1,i2,i3,dir,1)**2 )/dr(dir) 
                             end do
                            uTemp(i1,i2,i3,u1c)=un(i1,i2,i3,u1c)+(+2.*Dztu(i1,i2,i3,u1c)-(Dztu(i1+1,i2,i3,u1c)+Dztu(i1-1,i2,i3,u1c)))*adxSosup(0)+(+2.*Dztu(i1,i2,i3,u1c)-(Dztu(i1,i2+1,i3,u1c)+Dztu(i1,i2-1,i3,u1c)))*adxSosup(1)
                            uTemp(i1,i2,i3,u2c)=un(i1,i2,i3,u2c)+(+2.*Dztu(i1,i2,i3,u2c)-(Dztu(i1+1,i2,i3,u2c)+Dztu(i1-1,i2,i3,u2c)))*adxSosup(0)+(+2.*Dztu(i1,i2,i3,u2c)-(Dztu(i1,i2+1,i3,u2c)+Dztu(i1,i2-1,i3,u2c)))*adxSosup(1)
                      end do
                      end do
                      end do
                   end if ! end if bc==
                  end do ! end side
                  end do ! end axis
               end if
               ! -- Copy temp solution back to un ----
                 do i3=n3a,n3b
                 do i2=n2a,n2b
                 do i1=n1a,n1b
                 un(i1,i2,i3,u1c)=uTemp(i1,i2,i3,u1c)
                 un(i1,i2,i3,u2c)=uTemp(i1,i2,i3,u2c)
                 end do
                 end do
                 end do
           end if
         elseif( option==2 )then
           ! --- evaluate the time derivative of v for method-of-lines schemes ---
           if( twilightZone==1 )then
                 do i3=n3a,n3b
                 do i2=n2a,n2b
                 do i1=n1a,n1b
                 vt(i1,i2,i3,u1c) =  -rhoi*ux22(i1,i2,i3,pc) + cs*( ulaplacian22(i1,i2,i3,u1c) )
                 vt(i1,i2,i3,u2c) =  -rhoi*uy22(i1,i2,i3,pc) + cs*( ulaplacian22(i1,i2,i3,u2c) )
                     call ogDeriv(ep,2,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,u1c,uett)
                     call ogDeriv(ep,0,2,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,u1c,uexx)
                     call ogDeriv(ep,0,0,2,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,u1c,ueyy)
                     call ogDeriv(ep,2,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,u2c,vett)
                     call ogDeriv(ep,0,2,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,u2c,vexx)
                     call ogDeriv(ep,0,0,2,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,u2c,veyy) 
                     call ogDeriv(ep,0,1,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,pc,pex)                           
                     call ogDeriv(ep,0,0,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,pc,pey)
                     vt(i1,i2,i3,u1c) =  vt(i1,i2,i3,u1c) + uett + rhoi*pex - cs*( uexx + ueyy )
                     vt(i1,i2,i3,u2c) =  vt(i1,i2,i3,u2c) + vett + rhoi*pey - cs*( vexx + veyy )
                 end do
                 end do
                 end do
           else
                 do i3=n3a,n3b
                 do i2=n2a,n2b
                 do i1=n1a,n1b
                 vt(i1,i2,i3,u1c) =  -rhoi*ux22(i1,i2,i3,pc) + cs*( ulaplacian22(i1,i2,i3,u1c) )
                 vt(i1,i2,i3,u2c) =  -rhoi*uy22(i1,i2,i3,pc) + cs*( ulaplacian22(i1,i2,i3,u2c) )
                 end do
                 end do
                 end do
           end if
         elseif( option==3 )then
           ! -----METHOD-OF-LINES UPDATE -----
           if( orderOfAccuracyInTime==2 )then
             if( correction==0 )then
               ! AB2 predictor
               cmol1 =  (3./2.)*dt
               cmol2 = -(1./2.)*dt
             else
               ! AM2 corrector
               cmol1 = (1./2.)*dt
               cmol2 = (1./2.)*dt
             end if
               ! Update u first since v1 may be over-written by vn 
                 do i3=n3a,n3b
                 do i2=n2a,n2b
                 do i1=n1a,n1b
                 un(i1,i2,i3,u1c) = u(i1,i2,i3,u1c) + cmol1*v1(i1,i2,i3,v1c) + cmol2*v2(i1,i2,i3,v1c)
                 un(i1,i2,i3,u2c) = u(i1,i2,i3,u2c) + cmol1*v1(i1,i2,i3,v2c) + cmol2*v2(i1,i2,i3,v2c)
                 end do
                 end do
                 end do
                 do i3=n3a,n3b
                 do i2=n2a,n2b
                 do i1=n1a,n1b
                 vn(i1,i2,i3,v1c) = v(i1,i2,i3,v1c) + cmol1*vt1(i1,i2,i3,u1c) + cmol2*vt2(i1,i2,i3,u1c)
                 vn(i1,i2,i3,v2c) = v(i1,i2,i3,v2c) + cmol1*vt1(i1,i2,i3,u2c) + cmol2*vt2(i1,i2,i3,u2c)
                 end do
                 end do
                 end do
           else if( orderOfAccuracyInTime==4 )then
             if( correction==0 )then
               ! AB3 predictor
               cmol1 =  (23./12.)*dt
               cmol2 = -(16./12.)*dt
               cmol3 =  ( 5./12.)*dt
                   do i3=n3a,n3b
                   do i2=n2a,n2b
                   do i1=n1a,n1b
                   un(i1,i2,i3,u1c) = u(i1,i2,i3,u1c) + cmol1*v1(i1,i2,i3,v1c) + cmol2*v2(i1,i2,i3,v1c) + cmol3*v3(i1,i2,i3,v1c)
                   un(i1,i2,i3,u2c) = u(i1,i2,i3,u2c) + cmol1*v1(i1,i2,i3,v2c) + cmol2*v2(i1,i2,i3,v2c) + cmol3*v3(i1,i2,i3,v2c)
                   end do
                   end do
                   end do
                   do i3=n3a,n3b
                   do i2=n2a,n2b
                   do i1=n1a,n1b
                   vn(i1,i2,i3,v1c) = v(i1,i2,i3,v1c) + cmol1*vt1(i1,i2,i3,u1c) + cmol2*vt2(i1,i2,i3,u1c) + cmol3*vt3(i1,i2,i3,u1c)
                   vn(i1,i2,i3,v2c) = v(i1,i2,i3,v2c) + cmol1*vt1(i1,i2,i3,u2c) + cmol2*vt2(i1,i2,i3,u2c) + cmol3*vt3(i1,i2,i3,u2c)
                   end do
                   end do
                   end do
             else
               ! AM4 corrector
               cmol1 =  ( 9./24.)*dt
               cmol2 =  (19./24.)*dt
               cmol3 = -( 5./24.)*dt
               cmol4 =  ( 1./24.)*dt
                   do i3=n3a,n3b
                   do i2=n2a,n2b
                   do i1=n1a,n1b
                   un(i1,i2,i3,u1c) = u(i1,i2,i3,u1c) + cmol1*v1(i1,i2,i3,v1c) + cmol2*v2(i1,i2,i3,v1c) + cmol3*v3(i1,i2,i3,v1c) + cmol4*v4(i1,i2,i3,v1c)
                   un(i1,i2,i3,u2c) = u(i1,i2,i3,u2c) + cmol1*v1(i1,i2,i3,v2c) + cmol2*v2(i1,i2,i3,v2c) + cmol3*v3(i1,i2,i3,v2c) + cmol4*v4(i1,i2,i3,v2c)
                   end do
                   end do
                   end do
                   do i3=n3a,n3b
                   do i2=n2a,n2b
                   do i1=n1a,n1b
                   vn(i1,i2,i3,v1c) = v(i1,i2,i3,v1c) + cmol1*vt1(i1,i2,i3,u1c) + cmol2*vt2(i1,i2,i3,u1c) + cmol3*vt3(i1,i2,i3,u1c) + cmol4*vt4(i1,i2,i3,u1c)
                   vn(i1,i2,i3,v2c) = v(i1,i2,i3,v2c) + cmol1*vt1(i1,i2,i3,u2c) + cmol2*vt2(i1,i2,i3,u2c) + cmol3*vt3(i1,i2,i3,u2c) + cmol4*vt4(i1,i2,i3,u2c)
                   end do
                   end do
                   end do
             end if
           else
             write(*,*) 'advIsmOpt: not implemented for orderOfAccuracyInTime=',orderOfAccuracyInTime
             stop 6666
           end if
         else
           write(*,*) 'advIsmOpt: Unknown option=',option
           stop 2222
         end if
        return
        end
