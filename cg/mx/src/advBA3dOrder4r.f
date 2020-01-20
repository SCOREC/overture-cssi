! This file automatically generated from advBA.bf with bpp.
       subroutine advBA3dOrder4r(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,
     & nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,rsxy,  um,u,un,f,fa, K0i, 
     & matMask, pm,p,pt,xy, etax,etay,etaz, bc, dis, varDis, ipar, 
     & rpar, ierr )
      !======================================================================
      !   Advance a time step for Maxwells equations
      !     OPTIMIZED version for rectangular grids.
      ! nd : number of space dimensions
      !
      ! ipar(0)  = option : option=0 - Maxwell+Artificial diffusion
      !                           =1 - AD only
      !
      !  dis(i1,i2,i3) : temp space to hold artificial dissipation
      !  varDis(i1,i2,i3) : coefficient of the variable artificial dissipation
      !======================================================================
       implicit none
       integer nd, n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,
     & nd3b,nd4a,nd4b
       real um(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
       real u(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
       real un(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
       real f(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
       real fa(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b,0:*)  ! forcings at different times
       ! Polarization vectors
       real pm(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
       real p(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
       real pt(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
       real xy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
       ! Super-grid layer functions 
       real etax(nd1a:nd1b), etay(nd2a:nd2b), etaz(nd3a:nd3b)
      ! real vvt2(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
      ! real ut3(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
      ! real vvt4(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
      !  real ut5(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
      !  real ut6(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
      !  real ut7(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
       real dis(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
       real varDis(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
       real rsxy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1,0:nd-1)
       integer mask(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
       real K0i(0:5,0:5,0:*)  ! material matrix
       integer matMask(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
       integer bc(0:1,0:2),ierr
       integer ipar(0:*)
       real rpar(0:*)
      !     ---- local variables -----
       integer c,i1,i2,i3,n,gridType,orderOfAccuracy,orderInTime,k1,k2,
     & debug,numComp
       integer addForcing,orderOfDissipation,option
       integer useWhereMask,useWhereMaskSave,solveForE,solveForH,grid,
     & useVariableDissipation
       integer useCurvilinearOpt,useConservative,
     & combineDissipationWithAdvance,useDivergenceCleaning
       integer useNewForcingMethod,numberOfForcingFunctions,fcur,fnext,
     & fprev
       integer ex,ey,ez, hx,hy,hz, solveForAllFields, useSuperGrid, 
     & useAbsorbingLayer(0:2)
       real t,cc,dt,dy,dz,cdt,cdtdx,cdtdy,cdtdz,adc,adcdt,add,adddt
       real dt4by12
       real eps,mu,sigmaE,sigmaH,kx,ky,kz,divergenceCleaningCoefficient
       logical addDissipation, updateInterior, methodOfLines
       real fv(0:10) , curl(0:5)
       real ep
       integer maxRegions,NpMax
       parameter( maxRegions=100,NpMax=10 )  ! FIX ME
       integer numberOfMaterialRegions, mr
       real Ki(0:2,0:2,0:maxRegions) ! 3x3 material matrix for TEz polarization
       integer Np(6,6,0:maxRegions-1)
       real gdmPar(4,NpMax,6,6,0:maxRegions-1), ptSum(0:5)
       real a0,a1,b0,b1
       integer ec,pc,qc, pct,qct
       real e0,e0t,p0,p0t,q0,q0t
       real ev(0:5),evt(0:5)
       integer numPolarizationTerms
       integer maxNumPolarizationTerms
       parameter( maxNumPolarizationTerms=200 )
       real pv(0:maxNumPolarizationTerms), pvt(
     & 0:maxNumPolarizationTerms) !
       real fp(0:maxNumPolarizationTerms)
       real dx(0:2),dr(0:2)
       real adxSosup(0:2), sigma1
       real dx2i,dy2i,dz2i,dxsqi,dysqi,dzsqi,dxi,dyi,dzi
       real dx12i,dy12i,dz12i,dxsq12i,dysq12i,dzsq12i,dxy4i,dxz4i,dyz4,
     & time0,time1
       real dxi4,dyi4,dzi4,dxdyi2,dxdzi2,dydzi2
       real uLap(-1:1,-1:1,0:5),uLapSq(0:5)
       real uLaprr2,uLapss2,uLaprs2,uLapr2,uLaps2
       real c0,c1,csq,dtsq,cdtsq,cdtsq12,cdtSqBy12, csqdt
      !  real c0,c1,csq,dtsq,cdtsq,cdtsq12,lap(0:20),cdtSqBy12, csqdt
      !  real c40,c41,c42,c43
      !  real c60,c61,c62,c63,c64,c65
      !  real c80,c81,c82,c83,c84,c85,c86,c87
      ! real c00lap2d6,c10lap2d6,c01lap2d6,c20lap2d6,c02lap2d6,c30lap2d6,c03lap2d6
      !  real c00lap2d8,c10lap2d8,c01lap2d8,c20lap2d8,c02lap2d8,c30lap2d8,c03lap2d8,c40lap2d8,c04lap2d8
      !  real c000lap3d6,c100lap3d6,c010lap3d6,c001lap3d6,!                  c200lap3d6,c020lap3d6,c002lap3d6,!                  c300lap3d6,c030lap3d6,c003lap3d6
      !  real c000lap3d8,c100lap3d8,c010lap3d8,c001lap3d8,!                  c200lap3d8,c020lap3d8,c002lap3d8,!                  c300lap3d8,c030lap3d8,c003lap3d8,!                  c400lap3d8,c040lap3d8,c004lap3d8
       integer rectangular,curvilinear
       parameter( rectangular=0, curvilinear=1 )
       integer timeSteppingMethod
       integer defaultTimeStepping,adamsSymmetricOrder3,rungeKutta,
     & stoermerTimeStepping,modifiedEquationTimeStepping
       parameter(defaultTimeStepping=0,adamsSymmetricOrder3=1,
     & rungeKutta=2,stoermerTimeStepping=3,
     & modifiedEquationTimeStepping=4)
       integer materialType
       integer isotropic,bianisotropic
       parameter( isotropic=0, bianisotropic=1 )
      !...........start statement function
       integer kd,m
       real rx,ry,rz,sx,sy,sz,tx,ty,tz
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
        real vr2
        real vs2
        real vt2
        real vrr2
        real vss2
        real vrs2
        real vtt2
        real vrt2
        real vst2
        real vrrr2
        real vsss2
        real vttt2
        real vx21
        real vy21
        real vz21
        real vx22
        real vy22
        real vz22
        real vx23
        real vy23
        real vz23
        real vxx21
        real vyy21
        real vxy21
        real vxz21
        real vyz21
        real vzz21
        real vlaplacian21
        real vxx22
        real vyy22
        real vxy22
        real vxz22
        real vyz22
        real vzz22
        real vlaplacian22
        real vxx23
        real vyy23
        real vzz23
        real vxy23
        real vxz23
        real vyz23
        real vlaplacian23
        real vx23r
        real vy23r
        real vz23r
        real vxx23r
        real vyy23r
        real vxy23r
        real vzz23r
        real vxz23r
        real vyz23r
        real vx21r
        real vy21r
        real vz21r
        real vxx21r
        real vyy21r
        real vzz21r
        real vxy21r
        real vxz21r
        real vyz21r
        real vlaplacian21r
        real vx22r
        real vy22r
        real vz22r
        real vxx22r
        real vyy22r
        real vzz22r
        real vxy22r
        real vxz22r
        real vyz22r
        real vlaplacian22r
        real vlaplacian23r
        real vxxx22r
        real vyyy22r
        real vxxy22r
        real vxyy22r
        real vxxxx22r
        real vyyyy22r
        real vxxyy22r
        real vxxx23r
        real vyyy23r
        real vzzz23r
        real vxxy23r
        real vxxz23r
        real vxyy23r
        real vyyz23r
        real vxzz23r
        real vyzz23r
        real vxxxx23r
        real vyyyy23r
        real vzzzz23r
        real vxxyy23r
        real vxxzz23r
        real vyyzz23r
        real vLapSq22r
        real vLapSq23r
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
        real ffr2
        real ffs2
        real fft2
        real ffrr2
        real ffss2
        real ffrs2
        real fftt2
        real ffrt2
        real ffst2
        real ffrrr2
        real ffsss2
        real ffttt2
        real ffx21
        real ffy21
        real ffz21
        real ffx22
        real ffy22
        real ffz22
        real ffx23
        real ffy23
        real ffz23
        real ffxx21
        real ffyy21
        real ffxy21
        real ffxz21
        real ffyz21
        real ffzz21
        real fflaplacian21
        real ffxx22
        real ffyy22
        real ffxy22
        real ffxz22
        real ffyz22
        real ffzz22
        real fflaplacian22
        real ffxx23
        real ffyy23
        real ffzz23
        real ffxy23
        real ffxz23
        real ffyz23
        real fflaplacian23
        real ffx23r
        real ffy23r
        real ffz23r
        real ffxx23r
        real ffyy23r
        real ffxy23r
        real ffzz23r
        real ffxz23r
        real ffyz23r
        real ffx21r
        real ffy21r
        real ffz21r
        real ffxx21r
        real ffyy21r
        real ffzz21r
        real ffxy21r
        real ffxz21r
        real ffyz21r
        real fflaplacian21r
        real ffx22r
        real ffy22r
        real ffz22r
        real ffxx22r
        real ffyy22r
        real ffzz22r
        real ffxy22r
        real ffxz22r
        real ffyz22r
        real fflaplacian22r
        real fflaplacian23r
        real ffxxx22r
        real ffyyy22r
        real ffxxy22r
        real ffxyy22r
        real ffxxxx22r
        real ffyyyy22r
        real ffxxyy22r
        real ffxxx23r
        real ffyyy23r
        real ffzzz23r
        real ffxxy23r
        real ffxxz23r
        real ffxyy23r
        real ffyyz23r
        real ffxzz23r
        real ffyzz23r
        real ffxxxx23r
        real ffyyyy23r
        real ffzzzz23r
        real ffxxyy23r
        real ffxxzz23r
        real ffyyzz23r
        real ffLapSq22r
        real ffLapSq23r
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
        real vr4
        real vs4
        real vt4
        real vrr4
        real vss4
        real vtt4
        real vrs4
        real vrt4
        real vst4
        real vx41
        real vy41
        real vz41
        real vx42
        real vy42
        real vz42
        real vx43
        real vy43
        real vz43
        real vxx41
        real vyy41
        real vxy41
        real vxz41
        real vyz41
        real vzz41
        real vlaplacian41
        real vxx42
        real vyy42
        real vxy42
        real vxz42
        real vyz42
        real vzz42
        real vlaplacian42
        real vxx43
        real vyy43
        real vzz43
        real vxy43
        real vxz43
        real vyz43
        real vlaplacian43
        real vx43r
        real vy43r
        real vz43r
        real vxx43r
        real vyy43r
        real vzz43r
        real vxy43r
        real vxz43r
        real vyz43r
        real vx41r
        real vy41r
        real vz41r
        real vxx41r
        real vyy41r
        real vzz41r
        real vxy41r
        real vxz41r
        real vyz41r
        real vlaplacian41r
        real vx42r
        real vy42r
        real vz42r
        real vxx42r
        real vyy42r
        real vzz42r
        real vxy42r
        real vxz42r
        real vyz42r
        real vlaplacian42r
        real vlaplacian43r
       real maxwell2dr,maxwell3dr,maxwellr44,maxwellr66,maxwellr88
       real maxwellc22,maxwellc44,maxwellc66,maxwellc88
       real maxwell2dr44me,maxwell2dr66me,maxwell2dr88me
       real maxwell3dr44me,maxwell3dr66me,maxwell3dr88me
       real maxwellc44me,maxwellc66me,maxwellc88me
       real max2dc44me,max2dc44me2,max3dc44me
       real mxdc2d2Ex,mxdc2d2Ey,mxdc2d4Ex,mxdc2d4Ey, mxdc2d4cEx,
     & mxdc2d4cEy
       real mxdc2d2cEx,mxdc2d2cEy
       real mxdc3d2Ex,mxdc3d2Ey,mxdc3d2Ez,mxdc3d2Hx,mxdc3d2Hy,mxdc3d2Hz
       real mxdc3d2cEx,mxdc3d2cEy,mxdc3d2cEz,mxdc3d2cHx,mxdc3d2cHy,
     & mxdc3d2cHz
       real mxdc2d4cConsEx,mxdc2d4cConsEy,mxdc2d4cConsEz
       real mxdc3d4Ex,mxdc3d4Ey,mxdc3d4Ez,mxdc3d4Hx,mxdc3d4Hy,mxdc3d4Hz
      ! real vr2,vs2,vrr2,vss2,vrs2,vLaplacian22
       real cdt4by360,cdt6by20160
       real lap2d2,lap3d2,lap2d4,lap3d4,lap2d6,lap3d6,lap2d8,lap3d8,
     & lap2d2Pow2,lap3d2Pow2,lap2d2Pow3,lap3d2Pow3,lap2d2Pow4,
     & lap3d2Pow4,lap2d4Pow2,lap3d4Pow2,lap2d4Pow3,lap3d4Pow3,
     & lap2d6Pow2,lap3d6Pow2
       real lap2d2m,lap3d2m
       real du,fd22d,fd23d,fd42d,fd43d,fd62d,fd63d,fd82d,fd83d
       ! forcing correction functions: 
       real lap2d2f,f2drme44, lap3d2f, f3drme44, f2dcme44, f3dcme44, ff
       ! div cleaning: 
       real dc,dcp,cdc0,cdc1,cdcxx,cdcyy,cdczz,cdcEdx,cdcEdy,cdcEdz,
     & cdcHdx,cdcHdy,cdcHdz,cdcf
       real cdcE,cdcELap,cdcELapsq,cdcELapm,cdcHzxLap,cdcHzyLap
       real cdcH,cdcHLap,cdcHLapsq,cdcHLapm
       ! dispersion
       integer dispersionModel,pxc,pyc,pzc,qxc,qyc,qzc,rxc,ryc,rzc
       integer numTerms1(0:maxRegions),ecIndex1(
     & maxNumPolarizationTerms,0:maxRegions),qcIndex1(
     & maxNumPolarizationTerms,0:maxRegions)
       integer numTerms2(0:maxRegions),ecIndex2(
     & maxNumPolarizationTerms,0:maxRegions),pcIndex2(
     & maxNumPolarizationTerms,0:maxRegions)
       real a0v(maxNumPolarizationTerms,0:maxRegions),a1v(
     & maxNumPolarizationTerms,0:maxRegions)
       real b0v(maxNumPolarizationTerms,0:maxRegions),b1v(
     & maxNumPolarizationTerms,0:maxRegions)
       ! real uv(0:5), unv(0:5)
       logical useOpt
      ! Dispersion models
       integer noDispersion,drude,gdm
       parameter( noDispersion=0, drude=1, gdm=2 )
       integer forcingOption
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
       ! boundary conditions parameters
       ! #Include "bcDefineFortranInclude.h"
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
      !     The next macro will define the difference approximation statement functions
       d12(kd) = 1./(2.*dr(kd))
       d22(kd) = 1./(dr(kd)**2)
       ur2(i1,i2,i3,kd)=(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))*d12(0)
       us2(i1,i2,i3,kd)=(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))*d12(1)
       ut2(i1,i2,i3,kd)=(u(i1,i2,i3+1,kd)-u(i1,i2,i3-1,kd))*d12(2)
       urr2(i1,i2,i3,kd)=(-2.*u(i1,i2,i3,kd)+(u(i1+1,i2,i3,kd)+u(i1-1,
     & i2,i3,kd)) )*d22(0)
       uss2(i1,i2,i3,kd)=(-2.*u(i1,i2,i3,kd)+(u(i1,i2+1,i3,kd)+u(i1,i2-
     & 1,i3,kd)) )*d22(1)
       urs2(i1,i2,i3,kd)=(ur2(i1,i2+1,i3,kd)-ur2(i1,i2-1,i3,kd))*d12(1)
       utt2(i1,i2,i3,kd)=(-2.*u(i1,i2,i3,kd)+(u(i1,i2,i3+1,kd)+u(i1,i2,
     & i3-1,kd)) )*d22(2)
       urt2(i1,i2,i3,kd)=(ur2(i1,i2,i3+1,kd)-ur2(i1,i2,i3-1,kd))*d12(2)
       ust2(i1,i2,i3,kd)=(us2(i1,i2,i3+1,kd)-us2(i1,i2,i3-1,kd))*d12(2)
       urrr2(i1,i2,i3,kd)=(-2.*(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))+(u(
     & i1+2,i2,i3,kd)-u(i1-2,i2,i3,kd)) )*d22(0)*d12(0)
       usss2(i1,i2,i3,kd)=(-2.*(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))+(u(
     & i1,i2+2,i3,kd)-u(i1,i2-2,i3,kd)) )*d22(1)*d12(1)
       uttt2(i1,i2,i3,kd)=(-2.*(u(i1,i2,i3+1,kd)-u(i1,i2,i3-1,kd))+(u(
     & i1,i2,i3+2,kd)-u(i1,i2,i3-2,kd)) )*d22(2)*d12(2)
       rxr2(i1,i2,i3)=(rx(i1+1,i2,i3)-rx(i1-1,i2,i3))*d12(0)
       rxs2(i1,i2,i3)=(rx(i1,i2+1,i3)-rx(i1,i2-1,i3))*d12(1)
       rxt2(i1,i2,i3)=(rx(i1,i2,i3+1)-rx(i1,i2,i3-1))*d12(2)
       rxrr2(i1,i2,i3)=(-2.*rx(i1,i2,i3)+(rx(i1+1,i2,i3)+rx(i1-1,i2,i3)
     & ) )*d22(0)
       rxss2(i1,i2,i3)=(-2.*rx(i1,i2,i3)+(rx(i1,i2+1,i3)+rx(i1,i2-1,i3)
     & ) )*d22(1)
       rxrs2(i1,i2,i3)=(rxr2(i1,i2+1,i3)-rxr2(i1,i2-1,i3))*d12(1)
       ryr2(i1,i2,i3)=(ry(i1+1,i2,i3)-ry(i1-1,i2,i3))*d12(0)
       rys2(i1,i2,i3)=(ry(i1,i2+1,i3)-ry(i1,i2-1,i3))*d12(1)
       ryt2(i1,i2,i3)=(ry(i1,i2,i3+1)-ry(i1,i2,i3-1))*d12(2)
       ryrr2(i1,i2,i3)=(-2.*ry(i1,i2,i3)+(ry(i1+1,i2,i3)+ry(i1-1,i2,i3)
     & ) )*d22(0)
       ryss2(i1,i2,i3)=(-2.*ry(i1,i2,i3)+(ry(i1,i2+1,i3)+ry(i1,i2-1,i3)
     & ) )*d22(1)
       ryrs2(i1,i2,i3)=(ryr2(i1,i2+1,i3)-ryr2(i1,i2-1,i3))*d12(1)
       rzr2(i1,i2,i3)=(rz(i1+1,i2,i3)-rz(i1-1,i2,i3))*d12(0)
       rzs2(i1,i2,i3)=(rz(i1,i2+1,i3)-rz(i1,i2-1,i3))*d12(1)
       rzt2(i1,i2,i3)=(rz(i1,i2,i3+1)-rz(i1,i2,i3-1))*d12(2)
       rzrr2(i1,i2,i3)=(-2.*rz(i1,i2,i3)+(rz(i1+1,i2,i3)+rz(i1-1,i2,i3)
     & ) )*d22(0)
       rzss2(i1,i2,i3)=(-2.*rz(i1,i2,i3)+(rz(i1,i2+1,i3)+rz(i1,i2-1,i3)
     & ) )*d22(1)
       rzrs2(i1,i2,i3)=(rzr2(i1,i2+1,i3)-rzr2(i1,i2-1,i3))*d12(1)
       sxr2(i1,i2,i3)=(sx(i1+1,i2,i3)-sx(i1-1,i2,i3))*d12(0)
       sxs2(i1,i2,i3)=(sx(i1,i2+1,i3)-sx(i1,i2-1,i3))*d12(1)
       sxt2(i1,i2,i3)=(sx(i1,i2,i3+1)-sx(i1,i2,i3-1))*d12(2)
       sxrr2(i1,i2,i3)=(-2.*sx(i1,i2,i3)+(sx(i1+1,i2,i3)+sx(i1-1,i2,i3)
     & ) )*d22(0)
       sxss2(i1,i2,i3)=(-2.*sx(i1,i2,i3)+(sx(i1,i2+1,i3)+sx(i1,i2-1,i3)
     & ) )*d22(1)
       sxrs2(i1,i2,i3)=(sxr2(i1,i2+1,i3)-sxr2(i1,i2-1,i3))*d12(1)
       syr2(i1,i2,i3)=(sy(i1+1,i2,i3)-sy(i1-1,i2,i3))*d12(0)
       sys2(i1,i2,i3)=(sy(i1,i2+1,i3)-sy(i1,i2-1,i3))*d12(1)
       syt2(i1,i2,i3)=(sy(i1,i2,i3+1)-sy(i1,i2,i3-1))*d12(2)
       syrr2(i1,i2,i3)=(-2.*sy(i1,i2,i3)+(sy(i1+1,i2,i3)+sy(i1-1,i2,i3)
     & ) )*d22(0)
       syss2(i1,i2,i3)=(-2.*sy(i1,i2,i3)+(sy(i1,i2+1,i3)+sy(i1,i2-1,i3)
     & ) )*d22(1)
       syrs2(i1,i2,i3)=(syr2(i1,i2+1,i3)-syr2(i1,i2-1,i3))*d12(1)
       szr2(i1,i2,i3)=(sz(i1+1,i2,i3)-sz(i1-1,i2,i3))*d12(0)
       szs2(i1,i2,i3)=(sz(i1,i2+1,i3)-sz(i1,i2-1,i3))*d12(1)
       szt2(i1,i2,i3)=(sz(i1,i2,i3+1)-sz(i1,i2,i3-1))*d12(2)
       szrr2(i1,i2,i3)=(-2.*sz(i1,i2,i3)+(sz(i1+1,i2,i3)+sz(i1-1,i2,i3)
     & ) )*d22(0)
       szss2(i1,i2,i3)=(-2.*sz(i1,i2,i3)+(sz(i1,i2+1,i3)+sz(i1,i2-1,i3)
     & ) )*d22(1)
       szrs2(i1,i2,i3)=(szr2(i1,i2+1,i3)-szr2(i1,i2-1,i3))*d12(1)
       txr2(i1,i2,i3)=(tx(i1+1,i2,i3)-tx(i1-1,i2,i3))*d12(0)
       txs2(i1,i2,i3)=(tx(i1,i2+1,i3)-tx(i1,i2-1,i3))*d12(1)
       txt2(i1,i2,i3)=(tx(i1,i2,i3+1)-tx(i1,i2,i3-1))*d12(2)
       txrr2(i1,i2,i3)=(-2.*tx(i1,i2,i3)+(tx(i1+1,i2,i3)+tx(i1-1,i2,i3)
     & ) )*d22(0)
       txss2(i1,i2,i3)=(-2.*tx(i1,i2,i3)+(tx(i1,i2+1,i3)+tx(i1,i2-1,i3)
     & ) )*d22(1)
       txrs2(i1,i2,i3)=(txr2(i1,i2+1,i3)-txr2(i1,i2-1,i3))*d12(1)
       tyr2(i1,i2,i3)=(ty(i1+1,i2,i3)-ty(i1-1,i2,i3))*d12(0)
       tys2(i1,i2,i3)=(ty(i1,i2+1,i3)-ty(i1,i2-1,i3))*d12(1)
       tyt2(i1,i2,i3)=(ty(i1,i2,i3+1)-ty(i1,i2,i3-1))*d12(2)
       tyrr2(i1,i2,i3)=(-2.*ty(i1,i2,i3)+(ty(i1+1,i2,i3)+ty(i1-1,i2,i3)
     & ) )*d22(0)
       tyss2(i1,i2,i3)=(-2.*ty(i1,i2,i3)+(ty(i1,i2+1,i3)+ty(i1,i2-1,i3)
     & ) )*d22(1)
       tyrs2(i1,i2,i3)=(tyr2(i1,i2+1,i3)-tyr2(i1,i2-1,i3))*d12(1)
       tzr2(i1,i2,i3)=(tz(i1+1,i2,i3)-tz(i1-1,i2,i3))*d12(0)
       tzs2(i1,i2,i3)=(tz(i1,i2+1,i3)-tz(i1,i2-1,i3))*d12(1)
       tzt2(i1,i2,i3)=(tz(i1,i2,i3+1)-tz(i1,i2,i3-1))*d12(2)
       tzrr2(i1,i2,i3)=(-2.*tz(i1,i2,i3)+(tz(i1+1,i2,i3)+tz(i1-1,i2,i3)
     & ) )*d22(0)
       tzss2(i1,i2,i3)=(-2.*tz(i1,i2,i3)+(tz(i1,i2+1,i3)+tz(i1,i2-1,i3)
     & ) )*d22(1)
       tzrs2(i1,i2,i3)=(tzr2(i1,i2+1,i3)-tzr2(i1,i2-1,i3))*d12(1)
       ux21(i1,i2,i3,kd)= rx(i1,i2,i3)*ur2(i1,i2,i3,kd)
       uy21(i1,i2,i3,kd)=0
       uz21(i1,i2,i3,kd)=0
       ux22(i1,i2,i3,kd)= rx(i1,i2,i3)*ur2(i1,i2,i3,kd)+sx(i1,i2,i3)*
     & us2(i1,i2,i3,kd)
       uy22(i1,i2,i3,kd)= ry(i1,i2,i3)*ur2(i1,i2,i3,kd)+sy(i1,i2,i3)*
     & us2(i1,i2,i3,kd)
       uz22(i1,i2,i3,kd)=0
       ux23(i1,i2,i3,kd)=rx(i1,i2,i3)*ur2(i1,i2,i3,kd)+sx(i1,i2,i3)*
     & us2(i1,i2,i3,kd)+tx(i1,i2,i3)*ut2(i1,i2,i3,kd)
       uy23(i1,i2,i3,kd)=ry(i1,i2,i3)*ur2(i1,i2,i3,kd)+sy(i1,i2,i3)*
     & us2(i1,i2,i3,kd)+ty(i1,i2,i3)*ut2(i1,i2,i3,kd)
       uz23(i1,i2,i3,kd)=rz(i1,i2,i3)*ur2(i1,i2,i3,kd)+sz(i1,i2,i3)*
     & us2(i1,i2,i3,kd)+tz(i1,i2,i3)*ut2(i1,i2,i3,kd)
       rxx21(i1,i2,i3)= rx(i1,i2,i3)*rxr2(i1,i2,i3)
       rxx22(i1,i2,i3)= rx(i1,i2,i3)*rxr2(i1,i2,i3)+sx(i1,i2,i3)*rxs2(
     & i1,i2,i3)
       rxy22(i1,i2,i3)= ry(i1,i2,i3)*rxr2(i1,i2,i3)+sy(i1,i2,i3)*rxs2(
     & i1,i2,i3)
       rxx23(i1,i2,i3)=rx(i1,i2,i3)*rxr2(i1,i2,i3)+sx(i1,i2,i3)*rxs2(
     & i1,i2,i3)+tx(i1,i2,i3)*rxt2(i1,i2,i3)
       rxy23(i1,i2,i3)=ry(i1,i2,i3)*rxr2(i1,i2,i3)+sy(i1,i2,i3)*rxs2(
     & i1,i2,i3)+ty(i1,i2,i3)*rxt2(i1,i2,i3)
       rxz23(i1,i2,i3)=rz(i1,i2,i3)*rxr2(i1,i2,i3)+sz(i1,i2,i3)*rxs2(
     & i1,i2,i3)+tz(i1,i2,i3)*rxt2(i1,i2,i3)
       ryx22(i1,i2,i3)= rx(i1,i2,i3)*ryr2(i1,i2,i3)+sx(i1,i2,i3)*rys2(
     & i1,i2,i3)
       ryy22(i1,i2,i3)= ry(i1,i2,i3)*ryr2(i1,i2,i3)+sy(i1,i2,i3)*rys2(
     & i1,i2,i3)
       ryx23(i1,i2,i3)=rx(i1,i2,i3)*ryr2(i1,i2,i3)+sx(i1,i2,i3)*rys2(
     & i1,i2,i3)+tx(i1,i2,i3)*ryt2(i1,i2,i3)
       ryy23(i1,i2,i3)=ry(i1,i2,i3)*ryr2(i1,i2,i3)+sy(i1,i2,i3)*rys2(
     & i1,i2,i3)+ty(i1,i2,i3)*ryt2(i1,i2,i3)
       ryz23(i1,i2,i3)=rz(i1,i2,i3)*ryr2(i1,i2,i3)+sz(i1,i2,i3)*rys2(
     & i1,i2,i3)+tz(i1,i2,i3)*ryt2(i1,i2,i3)
       rzx22(i1,i2,i3)= rx(i1,i2,i3)*rzr2(i1,i2,i3)+sx(i1,i2,i3)*rzs2(
     & i1,i2,i3)
       rzy22(i1,i2,i3)= ry(i1,i2,i3)*rzr2(i1,i2,i3)+sy(i1,i2,i3)*rzs2(
     & i1,i2,i3)
       rzx23(i1,i2,i3)=rx(i1,i2,i3)*rzr2(i1,i2,i3)+sx(i1,i2,i3)*rzs2(
     & i1,i2,i3)+tx(i1,i2,i3)*rzt2(i1,i2,i3)
       rzy23(i1,i2,i3)=ry(i1,i2,i3)*rzr2(i1,i2,i3)+sy(i1,i2,i3)*rzs2(
     & i1,i2,i3)+ty(i1,i2,i3)*rzt2(i1,i2,i3)
       rzz23(i1,i2,i3)=rz(i1,i2,i3)*rzr2(i1,i2,i3)+sz(i1,i2,i3)*rzs2(
     & i1,i2,i3)+tz(i1,i2,i3)*rzt2(i1,i2,i3)
       sxx22(i1,i2,i3)= rx(i1,i2,i3)*sxr2(i1,i2,i3)+sx(i1,i2,i3)*sxs2(
     & i1,i2,i3)
       sxy22(i1,i2,i3)= ry(i1,i2,i3)*sxr2(i1,i2,i3)+sy(i1,i2,i3)*sxs2(
     & i1,i2,i3)
       sxx23(i1,i2,i3)=rx(i1,i2,i3)*sxr2(i1,i2,i3)+sx(i1,i2,i3)*sxs2(
     & i1,i2,i3)+tx(i1,i2,i3)*sxt2(i1,i2,i3)
       sxy23(i1,i2,i3)=ry(i1,i2,i3)*sxr2(i1,i2,i3)+sy(i1,i2,i3)*sxs2(
     & i1,i2,i3)+ty(i1,i2,i3)*sxt2(i1,i2,i3)
       sxz23(i1,i2,i3)=rz(i1,i2,i3)*sxr2(i1,i2,i3)+sz(i1,i2,i3)*sxs2(
     & i1,i2,i3)+tz(i1,i2,i3)*sxt2(i1,i2,i3)
       syx22(i1,i2,i3)= rx(i1,i2,i3)*syr2(i1,i2,i3)+sx(i1,i2,i3)*sys2(
     & i1,i2,i3)
       syy22(i1,i2,i3)= ry(i1,i2,i3)*syr2(i1,i2,i3)+sy(i1,i2,i3)*sys2(
     & i1,i2,i3)
       syx23(i1,i2,i3)=rx(i1,i2,i3)*syr2(i1,i2,i3)+sx(i1,i2,i3)*sys2(
     & i1,i2,i3)+tx(i1,i2,i3)*syt2(i1,i2,i3)
       syy23(i1,i2,i3)=ry(i1,i2,i3)*syr2(i1,i2,i3)+sy(i1,i2,i3)*sys2(
     & i1,i2,i3)+ty(i1,i2,i3)*syt2(i1,i2,i3)
       syz23(i1,i2,i3)=rz(i1,i2,i3)*syr2(i1,i2,i3)+sz(i1,i2,i3)*sys2(
     & i1,i2,i3)+tz(i1,i2,i3)*syt2(i1,i2,i3)
       szx22(i1,i2,i3)= rx(i1,i2,i3)*szr2(i1,i2,i3)+sx(i1,i2,i3)*szs2(
     & i1,i2,i3)
       szy22(i1,i2,i3)= ry(i1,i2,i3)*szr2(i1,i2,i3)+sy(i1,i2,i3)*szs2(
     & i1,i2,i3)
       szx23(i1,i2,i3)=rx(i1,i2,i3)*szr2(i1,i2,i3)+sx(i1,i2,i3)*szs2(
     & i1,i2,i3)+tx(i1,i2,i3)*szt2(i1,i2,i3)
       szy23(i1,i2,i3)=ry(i1,i2,i3)*szr2(i1,i2,i3)+sy(i1,i2,i3)*szs2(
     & i1,i2,i3)+ty(i1,i2,i3)*szt2(i1,i2,i3)
       szz23(i1,i2,i3)=rz(i1,i2,i3)*szr2(i1,i2,i3)+sz(i1,i2,i3)*szs2(
     & i1,i2,i3)+tz(i1,i2,i3)*szt2(i1,i2,i3)
       txx22(i1,i2,i3)= rx(i1,i2,i3)*txr2(i1,i2,i3)+sx(i1,i2,i3)*txs2(
     & i1,i2,i3)
       txy22(i1,i2,i3)= ry(i1,i2,i3)*txr2(i1,i2,i3)+sy(i1,i2,i3)*txs2(
     & i1,i2,i3)
       txx23(i1,i2,i3)=rx(i1,i2,i3)*txr2(i1,i2,i3)+sx(i1,i2,i3)*txs2(
     & i1,i2,i3)+tx(i1,i2,i3)*txt2(i1,i2,i3)
       txy23(i1,i2,i3)=ry(i1,i2,i3)*txr2(i1,i2,i3)+sy(i1,i2,i3)*txs2(
     & i1,i2,i3)+ty(i1,i2,i3)*txt2(i1,i2,i3)
       txz23(i1,i2,i3)=rz(i1,i2,i3)*txr2(i1,i2,i3)+sz(i1,i2,i3)*txs2(
     & i1,i2,i3)+tz(i1,i2,i3)*txt2(i1,i2,i3)
       tyx22(i1,i2,i3)= rx(i1,i2,i3)*tyr2(i1,i2,i3)+sx(i1,i2,i3)*tys2(
     & i1,i2,i3)
       tyy22(i1,i2,i3)= ry(i1,i2,i3)*tyr2(i1,i2,i3)+sy(i1,i2,i3)*tys2(
     & i1,i2,i3)
       tyx23(i1,i2,i3)=rx(i1,i2,i3)*tyr2(i1,i2,i3)+sx(i1,i2,i3)*tys2(
     & i1,i2,i3)+tx(i1,i2,i3)*tyt2(i1,i2,i3)
       tyy23(i1,i2,i3)=ry(i1,i2,i3)*tyr2(i1,i2,i3)+sy(i1,i2,i3)*tys2(
     & i1,i2,i3)+ty(i1,i2,i3)*tyt2(i1,i2,i3)
       tyz23(i1,i2,i3)=rz(i1,i2,i3)*tyr2(i1,i2,i3)+sz(i1,i2,i3)*tys2(
     & i1,i2,i3)+tz(i1,i2,i3)*tyt2(i1,i2,i3)
       tzx22(i1,i2,i3)= rx(i1,i2,i3)*tzr2(i1,i2,i3)+sx(i1,i2,i3)*tzs2(
     & i1,i2,i3)
       tzy22(i1,i2,i3)= ry(i1,i2,i3)*tzr2(i1,i2,i3)+sy(i1,i2,i3)*tzs2(
     & i1,i2,i3)
       tzx23(i1,i2,i3)=rx(i1,i2,i3)*tzr2(i1,i2,i3)+sx(i1,i2,i3)*tzs2(
     & i1,i2,i3)+tx(i1,i2,i3)*tzt2(i1,i2,i3)
       tzy23(i1,i2,i3)=ry(i1,i2,i3)*tzr2(i1,i2,i3)+sy(i1,i2,i3)*tzs2(
     & i1,i2,i3)+ty(i1,i2,i3)*tzt2(i1,i2,i3)
       tzz23(i1,i2,i3)=rz(i1,i2,i3)*tzr2(i1,i2,i3)+sz(i1,i2,i3)*tzs2(
     & i1,i2,i3)+tz(i1,i2,i3)*tzt2(i1,i2,i3)
       uxx21(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*urr2(i1,i2,i3,kd)+(rxx22(
     & i1,i2,i3))*ur2(i1,i2,i3,kd)
       uyy21(i1,i2,i3,kd)=0
       uxy21(i1,i2,i3,kd)=0
       uxz21(i1,i2,i3,kd)=0
       uyz21(i1,i2,i3,kd)=0
       uzz21(i1,i2,i3,kd)=0
       ulaplacian21(i1,i2,i3,kd)=uxx21(i1,i2,i3,kd)
       uxx22(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*urr2(i1,i2,i3,kd)+2.*(rx(
     & i1,i2,i3)*sx(i1,i2,i3))*urs2(i1,i2,i3,kd)+(sx(i1,i2,i3)**2)*
     & uss2(i1,i2,i3,kd)+(rxx22(i1,i2,i3))*ur2(i1,i2,i3,kd)+(sxx22(i1,
     & i2,i3))*us2(i1,i2,i3,kd)
       uyy22(i1,i2,i3,kd)=(ry(i1,i2,i3)**2)*urr2(i1,i2,i3,kd)+2.*(ry(
     & i1,i2,i3)*sy(i1,i2,i3))*urs2(i1,i2,i3,kd)+(sy(i1,i2,i3)**2)*
     & uss2(i1,i2,i3,kd)+(ryy22(i1,i2,i3))*ur2(i1,i2,i3,kd)+(syy22(i1,
     & i2,i3))*us2(i1,i2,i3,kd)
       uxy22(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*urr2(i1,i2,i3,kd)+(
     & rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*urs2(i1,
     & i2,i3,kd)+sx(i1,i2,i3)*sy(i1,i2,i3)*uss2(i1,i2,i3,kd)+rxy22(i1,
     & i2,i3)*ur2(i1,i2,i3,kd)+sxy22(i1,i2,i3)*us2(i1,i2,i3,kd)
       uxz22(i1,i2,i3,kd)=0
       uyz22(i1,i2,i3,kd)=0
       uzz22(i1,i2,i3,kd)=0
       ulaplacian22(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2)*
     & urr2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*
     & sy(i1,i2,i3))*urs2(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**
     & 2)*uss2(i1,i2,i3,kd)+(rxx22(i1,i2,i3)+ryy22(i1,i2,i3))*ur2(i1,
     & i2,i3,kd)+(sxx22(i1,i2,i3)+syy22(i1,i2,i3))*us2(i1,i2,i3,kd)
       uxx23(i1,i2,i3,kd)=rx(i1,i2,i3)**2*urr2(i1,i2,i3,kd)+sx(i1,i2,
     & i3)**2*uss2(i1,i2,i3,kd)+tx(i1,i2,i3)**2*utt2(i1,i2,i3,kd)+2.*
     & rx(i1,i2,i3)*sx(i1,i2,i3)*urs2(i1,i2,i3,kd)+2.*rx(i1,i2,i3)*tx(
     & i1,i2,i3)*urt2(i1,i2,i3,kd)+2.*sx(i1,i2,i3)*tx(i1,i2,i3)*ust2(
     & i1,i2,i3,kd)+rxx23(i1,i2,i3)*ur2(i1,i2,i3,kd)+sxx23(i1,i2,i3)*
     & us2(i1,i2,i3,kd)+txx23(i1,i2,i3)*ut2(i1,i2,i3,kd)
       uyy23(i1,i2,i3,kd)=ry(i1,i2,i3)**2*urr2(i1,i2,i3,kd)+sy(i1,i2,
     & i3)**2*uss2(i1,i2,i3,kd)+ty(i1,i2,i3)**2*utt2(i1,i2,i3,kd)+2.*
     & ry(i1,i2,i3)*sy(i1,i2,i3)*urs2(i1,i2,i3,kd)+2.*ry(i1,i2,i3)*ty(
     & i1,i2,i3)*urt2(i1,i2,i3,kd)+2.*sy(i1,i2,i3)*ty(i1,i2,i3)*ust2(
     & i1,i2,i3,kd)+ryy23(i1,i2,i3)*ur2(i1,i2,i3,kd)+syy23(i1,i2,i3)*
     & us2(i1,i2,i3,kd)+tyy23(i1,i2,i3)*ut2(i1,i2,i3,kd)
       uzz23(i1,i2,i3,kd)=rz(i1,i2,i3)**2*urr2(i1,i2,i3,kd)+sz(i1,i2,
     & i3)**2*uss2(i1,i2,i3,kd)+tz(i1,i2,i3)**2*utt2(i1,i2,i3,kd)+2.*
     & rz(i1,i2,i3)*sz(i1,i2,i3)*urs2(i1,i2,i3,kd)+2.*rz(i1,i2,i3)*tz(
     & i1,i2,i3)*urt2(i1,i2,i3,kd)+2.*sz(i1,i2,i3)*tz(i1,i2,i3)*ust2(
     & i1,i2,i3,kd)+rzz23(i1,i2,i3)*ur2(i1,i2,i3,kd)+szz23(i1,i2,i3)*
     & us2(i1,i2,i3,kd)+tzz23(i1,i2,i3)*ut2(i1,i2,i3,kd)
       uxy23(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*urr2(i1,i2,i3,kd)+
     & sx(i1,i2,i3)*sy(i1,i2,i3)*uss2(i1,i2,i3,kd)+tx(i1,i2,i3)*ty(i1,
     & i2,i3)*utt2(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,
     & i3)*sx(i1,i2,i3))*urs2(i1,i2,i3,kd)+(rx(i1,i2,i3)*ty(i1,i2,i3)+
     & ry(i1,i2,i3)*tx(i1,i2,i3))*urt2(i1,i2,i3,kd)+(sx(i1,i2,i3)*ty(
     & i1,i2,i3)+sy(i1,i2,i3)*tx(i1,i2,i3))*ust2(i1,i2,i3,kd)+rxy23(
     & i1,i2,i3)*ur2(i1,i2,i3,kd)+sxy23(i1,i2,i3)*us2(i1,i2,i3,kd)+
     & txy23(i1,i2,i3)*ut2(i1,i2,i3,kd)
       uxz23(i1,i2,i3,kd)=rx(i1,i2,i3)*rz(i1,i2,i3)*urr2(i1,i2,i3,kd)+
     & sx(i1,i2,i3)*sz(i1,i2,i3)*uss2(i1,i2,i3,kd)+tx(i1,i2,i3)*tz(i1,
     & i2,i3)*utt2(i1,i2,i3,kd)+(rx(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,i2,
     & i3)*sx(i1,i2,i3))*urs2(i1,i2,i3,kd)+(rx(i1,i2,i3)*tz(i1,i2,i3)+
     & rz(i1,i2,i3)*tx(i1,i2,i3))*urt2(i1,i2,i3,kd)+(sx(i1,i2,i3)*tz(
     & i1,i2,i3)+sz(i1,i2,i3)*tx(i1,i2,i3))*ust2(i1,i2,i3,kd)+rxz23(
     & i1,i2,i3)*ur2(i1,i2,i3,kd)+sxz23(i1,i2,i3)*us2(i1,i2,i3,kd)+
     & txz23(i1,i2,i3)*ut2(i1,i2,i3,kd)
       uyz23(i1,i2,i3,kd)=ry(i1,i2,i3)*rz(i1,i2,i3)*urr2(i1,i2,i3,kd)+
     & sy(i1,i2,i3)*sz(i1,i2,i3)*uss2(i1,i2,i3,kd)+ty(i1,i2,i3)*tz(i1,
     & i2,i3)*utt2(i1,i2,i3,kd)+(ry(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,i2,
     & i3)*sy(i1,i2,i3))*urs2(i1,i2,i3,kd)+(ry(i1,i2,i3)*tz(i1,i2,i3)+
     & rz(i1,i2,i3)*ty(i1,i2,i3))*urt2(i1,i2,i3,kd)+(sy(i1,i2,i3)*tz(
     & i1,i2,i3)+sz(i1,i2,i3)*ty(i1,i2,i3))*ust2(i1,i2,i3,kd)+ryz23(
     & i1,i2,i3)*ur2(i1,i2,i3,kd)+syz23(i1,i2,i3)*us2(i1,i2,i3,kd)+
     & tyz23(i1,i2,i3)*ut2(i1,i2,i3,kd)
       ulaplacian23(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2+rz(
     & i1,i2,i3)**2)*urr2(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**
     & 2+sz(i1,i2,i3)**2)*uss2(i1,i2,i3,kd)+(tx(i1,i2,i3)**2+ty(i1,i2,
     & i3)**2+tz(i1,i2,i3)**2)*utt2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(
     & i1,i2,i3)+ ry(i1,i2,i3)*sy(i1,i2,i3)+rz(i1,i2,i3)*sz(i1,i2,i3))
     & *urs2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*tx(i1,i2,i3)+ ry(i1,i2,i3)*
     & ty(i1,i2,i3)+rz(i1,i2,i3)*tz(i1,i2,i3))*urt2(i1,i2,i3,kd)+2.*(
     & sx(i1,i2,i3)*tx(i1,i2,i3)+ sy(i1,i2,i3)*ty(i1,i2,i3)+sz(i1,i2,
     & i3)*tz(i1,i2,i3))*ust2(i1,i2,i3,kd)+(rxx23(i1,i2,i3)+ryy23(i1,
     & i2,i3)+rzz23(i1,i2,i3))*ur2(i1,i2,i3,kd)+(sxx23(i1,i2,i3)+
     & syy23(i1,i2,i3)+szz23(i1,i2,i3))*us2(i1,i2,i3,kd)+(txx23(i1,i2,
     & i3)+tyy23(i1,i2,i3)+tzz23(i1,i2,i3))*ut2(i1,i2,i3,kd)
c============================================================================================
c Define derivatives for a rectangular grid
c
c============================================================================================
       h12(kd) = 1./(2.*dx(kd))
       h22(kd) = 1./(dx(kd)**2)
       ux23r(i1,i2,i3,kd)=(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))*h12(0)
       uy23r(i1,i2,i3,kd)=(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))*h12(1)
       uz23r(i1,i2,i3,kd)=(u(i1,i2,i3+1,kd)-u(i1,i2,i3-1,kd))*h12(2)
       uxx23r(i1,i2,i3,kd)=(-2.*u(i1,i2,i3,kd)+(u(i1+1,i2,i3,kd)+u(i1-
     & 1,i2,i3,kd)) )*h22(0)
       uyy23r(i1,i2,i3,kd)=(-2.*u(i1,i2,i3,kd)+(u(i1,i2+1,i3,kd)+u(i1,
     & i2-1,i3,kd)) )*h22(1)
       uxy23r(i1,i2,i3,kd)=(ux23r(i1,i2+1,i3,kd)-ux23r(i1,i2-1,i3,kd))*
     & h12(1)
       uzz23r(i1,i2,i3,kd)=(-2.*u(i1,i2,i3,kd)+(u(i1,i2,i3+1,kd)+u(i1,
     & i2,i3-1,kd)) )*h22(2)
       uxz23r(i1,i2,i3,kd)=(ux23r(i1,i2,i3+1,kd)-ux23r(i1,i2,i3-1,kd))*
     & h12(2)
       uyz23r(i1,i2,i3,kd)=(uy23r(i1,i2,i3+1,kd)-uy23r(i1,i2,i3-1,kd))*
     & h12(2)
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
       ulaplacian22r(i1,i2,i3,kd)=uxx23r(i1,i2,i3,kd)+uyy23r(i1,i2,i3,
     & kd)
       ulaplacian23r(i1,i2,i3,kd)=uxx23r(i1,i2,i3,kd)+uyy23r(i1,i2,i3,
     & kd)+uzz23r(i1,i2,i3,kd)
       uxxx22r(i1,i2,i3,kd)=(-2.*(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))+(
     & u(i1+2,i2,i3,kd)-u(i1-2,i2,i3,kd)) )*h22(0)*h12(0)
       uyyy22r(i1,i2,i3,kd)=(-2.*(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))+(
     & u(i1,i2+2,i3,kd)-u(i1,i2-2,i3,kd)) )*h22(1)*h12(1)
       uxxy22r(i1,i2,i3,kd)=( uxx22r(i1,i2+1,i3,kd)-uxx22r(i1,i2-1,i3,
     & kd))/(2.*dx(1))
       uxyy22r(i1,i2,i3,kd)=( uyy22r(i1+1,i2,i3,kd)-uyy22r(i1-1,i2,i3,
     & kd))/(2.*dx(0))
       uxxxx22r(i1,i2,i3,kd)=(6.*u(i1,i2,i3,kd)-4.*(u(i1+1,i2,i3,kd)+u(
     & i1-1,i2,i3,kd))+(u(i1+2,i2,i3,kd)+u(i1-2,i2,i3,kd)) )/(dx(0)**
     & 4)
       uyyyy22r(i1,i2,i3,kd)=(6.*u(i1,i2,i3,kd)-4.*(u(i1,i2+1,i3,kd)+u(
     & i1,i2-1,i3,kd))+(u(i1,i2+2,i3,kd)+u(i1,i2-2,i3,kd)) )/(dx(1)**
     & 4)
       uxxyy22r(i1,i2,i3,kd)=( 4.*u(i1,i2,i3,kd)     -2.*(u(i1+1,i2,i3,
     & kd)+u(i1-1,i2,i3,kd)+u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd))   +   (
     & u(i1+1,i2+1,i3,kd)+u(i1-1,i2+1,i3,kd)+u(i1+1,i2-1,i3,kd)+u(i1-
     & 1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)
       ! 2D laplacian squared = u.xxxx + 2 u.xxyy + u.yyyy
       uLapSq22r(i1,i2,i3,kd)= ( 6.*u(i1,i2,i3,kd)   - 4.*(u(i1+1,i2,
     & i3,kd)+u(i1-1,i2,i3,kd))    +(u(i1+2,i2,i3,kd)+u(i1-2,i2,i3,kd)
     & ) )/(dx(0)**4) +( 6.*u(i1,i2,i3,kd)    -4.*(u(i1,i2+1,i3,kd)+u(
     & i1,i2-1,i3,kd))    +(u(i1,i2+2,i3,kd)+u(i1,i2-2,i3,kd)) )/(dx(
     & 1)**4)  +( 8.*u(i1,i2,i3,kd)     -4.*(u(i1+1,i2,i3,kd)+u(i1-1,
     & i2,i3,kd)+u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd))   +2.*(u(i1+1,i2+
     & 1,i3,kd)+u(i1-1,i2+1,i3,kd)+u(i1+1,i2-1,i3,kd)+u(i1-1,i2-1,i3,
     & kd)) )/(dx(0)**2*dx(1)**2)
       uxxx23r(i1,i2,i3,kd)=(-2.*(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))+(
     & u(i1+2,i2,i3,kd)-u(i1-2,i2,i3,kd)) )*h22(0)*h12(0)
       uyyy23r(i1,i2,i3,kd)=(-2.*(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))+(
     & u(i1,i2+2,i3,kd)-u(i1,i2-2,i3,kd)) )*h22(1)*h12(1)
       uzzz23r(i1,i2,i3,kd)=(-2.*(u(i1,i2,i3+1,kd)-u(i1,i2,i3-1,kd))+(
     & u(i1,i2,i3+2,kd)-u(i1,i2,i3-2,kd)) )*h22(1)*h12(2)
       uxxy23r(i1,i2,i3,kd)=( uxx22r(i1,i2+1,i3,kd)-uxx22r(i1,i2-1,i3,
     & kd))/(2.*dx(1))
       uxyy23r(i1,i2,i3,kd)=( uyy22r(i1+1,i2,i3,kd)-uyy22r(i1-1,i2,i3,
     & kd))/(2.*dx(0))
       uxxz23r(i1,i2,i3,kd)=( uxx22r(i1,i2,i3+1,kd)-uxx22r(i1,i2,i3-1,
     & kd))/(2.*dx(2))
       uyyz23r(i1,i2,i3,kd)=( uyy22r(i1,i2,i3+1,kd)-uyy22r(i1,i2,i3-1,
     & kd))/(2.*dx(2))
       uxzz23r(i1,i2,i3,kd)=( uzz22r(i1+1,i2,i3,kd)-uzz22r(i1-1,i2,i3,
     & kd))/(2.*dx(0))
       uyzz23r(i1,i2,i3,kd)=( uzz22r(i1,i2+1,i3,kd)-uzz22r(i1,i2-1,i3,
     & kd))/(2.*dx(1))
       uxxxx23r(i1,i2,i3,kd)=(6.*u(i1,i2,i3,kd)-4.*(u(i1+1,i2,i3,kd)+u(
     & i1-1,i2,i3,kd))+(u(i1+2,i2,i3,kd)+u(i1-2,i2,i3,kd)) )/(dx(0)**
     & 4)
       uyyyy23r(i1,i2,i3,kd)=(6.*u(i1,i2,i3,kd)-4.*(u(i1,i2+1,i3,kd)+u(
     & i1,i2-1,i3,kd))+(u(i1,i2+2,i3,kd)+u(i1,i2-2,i3,kd)) )/(dx(1)**
     & 4)
       uzzzz23r(i1,i2,i3,kd)=(6.*u(i1,i2,i3,kd)-4.*(u(i1,i2,i3+1,kd)+u(
     & i1,i2,i3-1,kd))+(u(i1,i2,i3+2,kd)+u(i1,i2,i3-2,kd)) )/(dx(2)**
     & 4)
       uxxyy23r(i1,i2,i3,kd)=( 4.*u(i1,i2,i3,kd)     -2.*(u(i1+1,i2,i3,
     & kd)+u(i1-1,i2,i3,kd)+u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd))   +   (
     & u(i1+1,i2+1,i3,kd)+u(i1-1,i2+1,i3,kd)+u(i1+1,i2-1,i3,kd)+u(i1-
     & 1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)
       uxxzz23r(i1,i2,i3,kd)=( 4.*u(i1,i2,i3,kd)     -2.*(u(i1+1,i2,i3,
     & kd)+u(i1-1,i2,i3,kd)+u(i1,i2,i3+1,kd)+u(i1,i2,i3-1,kd))   +   (
     & u(i1+1,i2,i3+1,kd)+u(i1-1,i2,i3+1,kd)+u(i1+1,i2,i3-1,kd)+u(i1-
     & 1,i2,i3-1,kd)) )/(dx(0)**2*dx(2)**2)
       uyyzz23r(i1,i2,i3,kd)=( 4.*u(i1,i2,i3,kd)     -2.*(u(i1,i2+1,i3,
     & kd)  +u(i1,i2-1,i3,kd)+  u(i1,i2  ,i3+1,kd)+u(i1,i2  ,i3-1,kd))
     &    +   (u(i1,i2+1,i3+1,kd)+u(i1,i2-1,i3+1,kd)+u(i1,i2+1,i3-1,
     & kd)+u(i1,i2-1,i3-1,kd)) )/(dx(1)**2*dx(2)**2)
       ! 3D laplacian squared = u.xxxx + u.yyyy + u.zzzz + 2 (u.xxyy + u.xxzz + u.yyzz )
       uLapSq23r(i1,i2,i3,kd)= ( 6.*u(i1,i2,i3,kd)   - 4.*(u(i1+1,i2,
     & i3,kd)+u(i1-1,i2,i3,kd))    +(u(i1+2,i2,i3,kd)+u(i1-2,i2,i3,kd)
     & ) )/(dx(0)**4) +( 6.*u(i1,i2,i3,kd)    -4.*(u(i1,i2+1,i3,kd)+u(
     & i1,i2-1,i3,kd))    +(u(i1,i2+2,i3,kd)+u(i1,i2-2,i3,kd)) )/(dx(
     & 1)**4)  +( 6.*u(i1,i2,i3,kd)    -4.*(u(i1,i2,i3+1,kd)+u(i1,i2,
     & i3-1,kd))    +(u(i1,i2,i3+2,kd)+u(i1,i2,i3-2,kd)) )/(dx(2)**4) 
     &  +( 8.*u(i1,i2,i3,kd)     -4.*(u(i1+1,i2,i3,kd)  +u(i1-1,i2,i3,
     & kd)  +u(i1  ,i2+1,i3,kd)+u(i1  ,i2-1,i3,kd))   +2.*(u(i1+1,i2+
     & 1,i3,kd)+u(i1-1,i2+1,i3,kd)+u(i1+1,i2-1,i3,kd)+u(i1-1,i2-1,i3,
     & kd)) )/(dx(0)**2*dx(1)**2)+( 8.*u(i1,i2,i3,kd)     -4.*(u(i1+1,
     & i2,i3,kd)  +u(i1-1,i2,i3,kd)  +u(i1  ,i2,i3+1,kd)+u(i1  ,i2,i3-
     & 1,kd))   +2.*(u(i1+1,i2,i3+1,kd)+u(i1-1,i2,i3+1,kd)+u(i1+1,i2,
     & i3-1,kd)+u(i1-1,i2,i3-1,kd)) )/(dx(0)**2*dx(2)**2)+( 8.*u(i1,
     & i2,i3,kd)     -4.*(u(i1,i2+1,i3,kd)  +u(i1,i2-1,i3,kd)  +u(i1,
     & i2  ,i3+1,kd)+u(i1,i2  ,i3-1,kd))   +2.*(u(i1,i2+1,i3+1,kd)+u(
     & i1,i2-1,i3+1,kd)+u(i1,i2+1,i3-1,kd)+u(i1,i2-1,i3-1,kd)) )/(dx(
     & 1)**2*dx(2)**2)
       d14(kd) = 1./(12.*dr(kd))
       d24(kd) = 1./(12.*dr(kd)**2)
       ur4(i1,i2,i3,kd)=(8.*(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))-(u(i1+
     & 2,i2,i3,kd)-u(i1-2,i2,i3,kd)))*d14(0)
       us4(i1,i2,i3,kd)=(8.*(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))-(u(i1,
     & i2+2,i3,kd)-u(i1,i2-2,i3,kd)))*d14(1)
       ut4(i1,i2,i3,kd)=(8.*(u(i1,i2,i3+1,kd)-u(i1,i2,i3-1,kd))-(u(i1,
     & i2,i3+2,kd)-u(i1,i2,i3-2,kd)))*d14(2)
       urr4(i1,i2,i3,kd)=(-30.*u(i1,i2,i3,kd)+16.*(u(i1+1,i2,i3,kd)+u(
     & i1-1,i2,i3,kd))-(u(i1+2,i2,i3,kd)+u(i1-2,i2,i3,kd)) )*d24(0)
       uss4(i1,i2,i3,kd)=(-30.*u(i1,i2,i3,kd)+16.*(u(i1,i2+1,i3,kd)+u(
     & i1,i2-1,i3,kd))-(u(i1,i2+2,i3,kd)+u(i1,i2-2,i3,kd)) )*d24(1)
       utt4(i1,i2,i3,kd)=(-30.*u(i1,i2,i3,kd)+16.*(u(i1,i2,i3+1,kd)+u(
     & i1,i2,i3-1,kd))-(u(i1,i2,i3+2,kd)+u(i1,i2,i3-2,kd)) )*d24(2)
       urs4(i1,i2,i3,kd)=(8.*(ur4(i1,i2+1,i3,kd)-ur4(i1,i2-1,i3,kd))-(
     & ur4(i1,i2+2,i3,kd)-ur4(i1,i2-2,i3,kd)))*d14(1)
       urt4(i1,i2,i3,kd)=(8.*(ur4(i1,i2,i3+1,kd)-ur4(i1,i2,i3-1,kd))-(
     & ur4(i1,i2,i3+2,kd)-ur4(i1,i2,i3-2,kd)))*d14(2)
       ust4(i1,i2,i3,kd)=(8.*(us4(i1,i2,i3+1,kd)-us4(i1,i2,i3-1,kd))-(
     & us4(i1,i2,i3+2,kd)-us4(i1,i2,i3-2,kd)))*d14(2)
       rxr4(i1,i2,i3)=(8.*(rx(i1+1,i2,i3)-rx(i1-1,i2,i3))-(rx(i1+2,i2,
     & i3)-rx(i1-2,i2,i3)))*d14(0)
       rxs4(i1,i2,i3)=(8.*(rx(i1,i2+1,i3)-rx(i1,i2-1,i3))-(rx(i1,i2+2,
     & i3)-rx(i1,i2-2,i3)))*d14(1)
       rxt4(i1,i2,i3)=(8.*(rx(i1,i2,i3+1)-rx(i1,i2,i3-1))-(rx(i1,i2,i3+
     & 2)-rx(i1,i2,i3-2)))*d14(2)
       ryr4(i1,i2,i3)=(8.*(ry(i1+1,i2,i3)-ry(i1-1,i2,i3))-(ry(i1+2,i2,
     & i3)-ry(i1-2,i2,i3)))*d14(0)
       rys4(i1,i2,i3)=(8.*(ry(i1,i2+1,i3)-ry(i1,i2-1,i3))-(ry(i1,i2+2,
     & i3)-ry(i1,i2-2,i3)))*d14(1)
       ryt4(i1,i2,i3)=(8.*(ry(i1,i2,i3+1)-ry(i1,i2,i3-1))-(ry(i1,i2,i3+
     & 2)-ry(i1,i2,i3-2)))*d14(2)
       rzr4(i1,i2,i3)=(8.*(rz(i1+1,i2,i3)-rz(i1-1,i2,i3))-(rz(i1+2,i2,
     & i3)-rz(i1-2,i2,i3)))*d14(0)
       rzs4(i1,i2,i3)=(8.*(rz(i1,i2+1,i3)-rz(i1,i2-1,i3))-(rz(i1,i2+2,
     & i3)-rz(i1,i2-2,i3)))*d14(1)
       rzt4(i1,i2,i3)=(8.*(rz(i1,i2,i3+1)-rz(i1,i2,i3-1))-(rz(i1,i2,i3+
     & 2)-rz(i1,i2,i3-2)))*d14(2)
       sxr4(i1,i2,i3)=(8.*(sx(i1+1,i2,i3)-sx(i1-1,i2,i3))-(sx(i1+2,i2,
     & i3)-sx(i1-2,i2,i3)))*d14(0)
       sxs4(i1,i2,i3)=(8.*(sx(i1,i2+1,i3)-sx(i1,i2-1,i3))-(sx(i1,i2+2,
     & i3)-sx(i1,i2-2,i3)))*d14(1)
       sxt4(i1,i2,i3)=(8.*(sx(i1,i2,i3+1)-sx(i1,i2,i3-1))-(sx(i1,i2,i3+
     & 2)-sx(i1,i2,i3-2)))*d14(2)
       syr4(i1,i2,i3)=(8.*(sy(i1+1,i2,i3)-sy(i1-1,i2,i3))-(sy(i1+2,i2,
     & i3)-sy(i1-2,i2,i3)))*d14(0)
       sys4(i1,i2,i3)=(8.*(sy(i1,i2+1,i3)-sy(i1,i2-1,i3))-(sy(i1,i2+2,
     & i3)-sy(i1,i2-2,i3)))*d14(1)
       syt4(i1,i2,i3)=(8.*(sy(i1,i2,i3+1)-sy(i1,i2,i3-1))-(sy(i1,i2,i3+
     & 2)-sy(i1,i2,i3-2)))*d14(2)
       szr4(i1,i2,i3)=(8.*(sz(i1+1,i2,i3)-sz(i1-1,i2,i3))-(sz(i1+2,i2,
     & i3)-sz(i1-2,i2,i3)))*d14(0)
       szs4(i1,i2,i3)=(8.*(sz(i1,i2+1,i3)-sz(i1,i2-1,i3))-(sz(i1,i2+2,
     & i3)-sz(i1,i2-2,i3)))*d14(1)
       szt4(i1,i2,i3)=(8.*(sz(i1,i2,i3+1)-sz(i1,i2,i3-1))-(sz(i1,i2,i3+
     & 2)-sz(i1,i2,i3-2)))*d14(2)
       txr4(i1,i2,i3)=(8.*(tx(i1+1,i2,i3)-tx(i1-1,i2,i3))-(tx(i1+2,i2,
     & i3)-tx(i1-2,i2,i3)))*d14(0)
       txs4(i1,i2,i3)=(8.*(tx(i1,i2+1,i3)-tx(i1,i2-1,i3))-(tx(i1,i2+2,
     & i3)-tx(i1,i2-2,i3)))*d14(1)
       txt4(i1,i2,i3)=(8.*(tx(i1,i2,i3+1)-tx(i1,i2,i3-1))-(tx(i1,i2,i3+
     & 2)-tx(i1,i2,i3-2)))*d14(2)
       tyr4(i1,i2,i3)=(8.*(ty(i1+1,i2,i3)-ty(i1-1,i2,i3))-(ty(i1+2,i2,
     & i3)-ty(i1-2,i2,i3)))*d14(0)
       tys4(i1,i2,i3)=(8.*(ty(i1,i2+1,i3)-ty(i1,i2-1,i3))-(ty(i1,i2+2,
     & i3)-ty(i1,i2-2,i3)))*d14(1)
       tyt4(i1,i2,i3)=(8.*(ty(i1,i2,i3+1)-ty(i1,i2,i3-1))-(ty(i1,i2,i3+
     & 2)-ty(i1,i2,i3-2)))*d14(2)
       tzr4(i1,i2,i3)=(8.*(tz(i1+1,i2,i3)-tz(i1-1,i2,i3))-(tz(i1+2,i2,
     & i3)-tz(i1-2,i2,i3)))*d14(0)
       tzs4(i1,i2,i3)=(8.*(tz(i1,i2+1,i3)-tz(i1,i2-1,i3))-(tz(i1,i2+2,
     & i3)-tz(i1,i2-2,i3)))*d14(1)
       tzt4(i1,i2,i3)=(8.*(tz(i1,i2,i3+1)-tz(i1,i2,i3-1))-(tz(i1,i2,i3+
     & 2)-tz(i1,i2,i3-2)))*d14(2)
       ux41(i1,i2,i3,kd)= rx(i1,i2,i3)*ur4(i1,i2,i3,kd)
       uy41(i1,i2,i3,kd)=0
       uz41(i1,i2,i3,kd)=0
       ux42(i1,i2,i3,kd)= rx(i1,i2,i3)*ur4(i1,i2,i3,kd)+sx(i1,i2,i3)*
     & us4(i1,i2,i3,kd)
       uy42(i1,i2,i3,kd)= ry(i1,i2,i3)*ur4(i1,i2,i3,kd)+sy(i1,i2,i3)*
     & us4(i1,i2,i3,kd)
       uz42(i1,i2,i3,kd)=0
       ux43(i1,i2,i3,kd)=rx(i1,i2,i3)*ur4(i1,i2,i3,kd)+sx(i1,i2,i3)*
     & us4(i1,i2,i3,kd)+tx(i1,i2,i3)*ut4(i1,i2,i3,kd)
       uy43(i1,i2,i3,kd)=ry(i1,i2,i3)*ur4(i1,i2,i3,kd)+sy(i1,i2,i3)*
     & us4(i1,i2,i3,kd)+ty(i1,i2,i3)*ut4(i1,i2,i3,kd)
       uz43(i1,i2,i3,kd)=rz(i1,i2,i3)*ur4(i1,i2,i3,kd)+sz(i1,i2,i3)*
     & us4(i1,i2,i3,kd)+tz(i1,i2,i3)*ut4(i1,i2,i3,kd)
       rxx41(i1,i2,i3)= rx(i1,i2,i3)*rxr4(i1,i2,i3)
       rxx42(i1,i2,i3)= rx(i1,i2,i3)*rxr4(i1,i2,i3)+sx(i1,i2,i3)*rxs4(
     & i1,i2,i3)
       rxy42(i1,i2,i3)= ry(i1,i2,i3)*rxr4(i1,i2,i3)+sy(i1,i2,i3)*rxs4(
     & i1,i2,i3)
       rxx43(i1,i2,i3)=rx(i1,i2,i3)*rxr4(i1,i2,i3)+sx(i1,i2,i3)*rxs4(
     & i1,i2,i3)+tx(i1,i2,i3)*rxt4(i1,i2,i3)
       rxy43(i1,i2,i3)=ry(i1,i2,i3)*rxr4(i1,i2,i3)+sy(i1,i2,i3)*rxs4(
     & i1,i2,i3)+ty(i1,i2,i3)*rxt4(i1,i2,i3)
       rxz43(i1,i2,i3)=rz(i1,i2,i3)*rxr4(i1,i2,i3)+sz(i1,i2,i3)*rxs4(
     & i1,i2,i3)+tz(i1,i2,i3)*rxt4(i1,i2,i3)
       ryx42(i1,i2,i3)= rx(i1,i2,i3)*ryr4(i1,i2,i3)+sx(i1,i2,i3)*rys4(
     & i1,i2,i3)
       ryy42(i1,i2,i3)= ry(i1,i2,i3)*ryr4(i1,i2,i3)+sy(i1,i2,i3)*rys4(
     & i1,i2,i3)
       ryx43(i1,i2,i3)=rx(i1,i2,i3)*ryr4(i1,i2,i3)+sx(i1,i2,i3)*rys4(
     & i1,i2,i3)+tx(i1,i2,i3)*ryt4(i1,i2,i3)
       ryy43(i1,i2,i3)=ry(i1,i2,i3)*ryr4(i1,i2,i3)+sy(i1,i2,i3)*rys4(
     & i1,i2,i3)+ty(i1,i2,i3)*ryt4(i1,i2,i3)
       ryz43(i1,i2,i3)=rz(i1,i2,i3)*ryr4(i1,i2,i3)+sz(i1,i2,i3)*rys4(
     & i1,i2,i3)+tz(i1,i2,i3)*ryt4(i1,i2,i3)
       rzx42(i1,i2,i3)= rx(i1,i2,i3)*rzr4(i1,i2,i3)+sx(i1,i2,i3)*rzs4(
     & i1,i2,i3)
       rzy42(i1,i2,i3)= ry(i1,i2,i3)*rzr4(i1,i2,i3)+sy(i1,i2,i3)*rzs4(
     & i1,i2,i3)
       rzx43(i1,i2,i3)=rx(i1,i2,i3)*rzr4(i1,i2,i3)+sx(i1,i2,i3)*rzs4(
     & i1,i2,i3)+tx(i1,i2,i3)*rzt4(i1,i2,i3)
       rzy43(i1,i2,i3)=ry(i1,i2,i3)*rzr4(i1,i2,i3)+sy(i1,i2,i3)*rzs4(
     & i1,i2,i3)+ty(i1,i2,i3)*rzt4(i1,i2,i3)
       rzz43(i1,i2,i3)=rz(i1,i2,i3)*rzr4(i1,i2,i3)+sz(i1,i2,i3)*rzs4(
     & i1,i2,i3)+tz(i1,i2,i3)*rzt4(i1,i2,i3)
       sxx42(i1,i2,i3)= rx(i1,i2,i3)*sxr4(i1,i2,i3)+sx(i1,i2,i3)*sxs4(
     & i1,i2,i3)
       sxy42(i1,i2,i3)= ry(i1,i2,i3)*sxr4(i1,i2,i3)+sy(i1,i2,i3)*sxs4(
     & i1,i2,i3)
       sxx43(i1,i2,i3)=rx(i1,i2,i3)*sxr4(i1,i2,i3)+sx(i1,i2,i3)*sxs4(
     & i1,i2,i3)+tx(i1,i2,i3)*sxt4(i1,i2,i3)
       sxy43(i1,i2,i3)=ry(i1,i2,i3)*sxr4(i1,i2,i3)+sy(i1,i2,i3)*sxs4(
     & i1,i2,i3)+ty(i1,i2,i3)*sxt4(i1,i2,i3)
       sxz43(i1,i2,i3)=rz(i1,i2,i3)*sxr4(i1,i2,i3)+sz(i1,i2,i3)*sxs4(
     & i1,i2,i3)+tz(i1,i2,i3)*sxt4(i1,i2,i3)
       syx42(i1,i2,i3)= rx(i1,i2,i3)*syr4(i1,i2,i3)+sx(i1,i2,i3)*sys4(
     & i1,i2,i3)
       syy42(i1,i2,i3)= ry(i1,i2,i3)*syr4(i1,i2,i3)+sy(i1,i2,i3)*sys4(
     & i1,i2,i3)
       syx43(i1,i2,i3)=rx(i1,i2,i3)*syr4(i1,i2,i3)+sx(i1,i2,i3)*sys4(
     & i1,i2,i3)+tx(i1,i2,i3)*syt4(i1,i2,i3)
       syy43(i1,i2,i3)=ry(i1,i2,i3)*syr4(i1,i2,i3)+sy(i1,i2,i3)*sys4(
     & i1,i2,i3)+ty(i1,i2,i3)*syt4(i1,i2,i3)
       syz43(i1,i2,i3)=rz(i1,i2,i3)*syr4(i1,i2,i3)+sz(i1,i2,i3)*sys4(
     & i1,i2,i3)+tz(i1,i2,i3)*syt4(i1,i2,i3)
       szx42(i1,i2,i3)= rx(i1,i2,i3)*szr4(i1,i2,i3)+sx(i1,i2,i3)*szs4(
     & i1,i2,i3)
       szy42(i1,i2,i3)= ry(i1,i2,i3)*szr4(i1,i2,i3)+sy(i1,i2,i3)*szs4(
     & i1,i2,i3)
       szx43(i1,i2,i3)=rx(i1,i2,i3)*szr4(i1,i2,i3)+sx(i1,i2,i3)*szs4(
     & i1,i2,i3)+tx(i1,i2,i3)*szt4(i1,i2,i3)
       szy43(i1,i2,i3)=ry(i1,i2,i3)*szr4(i1,i2,i3)+sy(i1,i2,i3)*szs4(
     & i1,i2,i3)+ty(i1,i2,i3)*szt4(i1,i2,i3)
       szz43(i1,i2,i3)=rz(i1,i2,i3)*szr4(i1,i2,i3)+sz(i1,i2,i3)*szs4(
     & i1,i2,i3)+tz(i1,i2,i3)*szt4(i1,i2,i3)
       txx42(i1,i2,i3)= rx(i1,i2,i3)*txr4(i1,i2,i3)+sx(i1,i2,i3)*txs4(
     & i1,i2,i3)
       txy42(i1,i2,i3)= ry(i1,i2,i3)*txr4(i1,i2,i3)+sy(i1,i2,i3)*txs4(
     & i1,i2,i3)
       txx43(i1,i2,i3)=rx(i1,i2,i3)*txr4(i1,i2,i3)+sx(i1,i2,i3)*txs4(
     & i1,i2,i3)+tx(i1,i2,i3)*txt4(i1,i2,i3)
       txy43(i1,i2,i3)=ry(i1,i2,i3)*txr4(i1,i2,i3)+sy(i1,i2,i3)*txs4(
     & i1,i2,i3)+ty(i1,i2,i3)*txt4(i1,i2,i3)
       txz43(i1,i2,i3)=rz(i1,i2,i3)*txr4(i1,i2,i3)+sz(i1,i2,i3)*txs4(
     & i1,i2,i3)+tz(i1,i2,i3)*txt4(i1,i2,i3)
       tyx42(i1,i2,i3)= rx(i1,i2,i3)*tyr4(i1,i2,i3)+sx(i1,i2,i3)*tys4(
     & i1,i2,i3)
       tyy42(i1,i2,i3)= ry(i1,i2,i3)*tyr4(i1,i2,i3)+sy(i1,i2,i3)*tys4(
     & i1,i2,i3)
       tyx43(i1,i2,i3)=rx(i1,i2,i3)*tyr4(i1,i2,i3)+sx(i1,i2,i3)*tys4(
     & i1,i2,i3)+tx(i1,i2,i3)*tyt4(i1,i2,i3)
       tyy43(i1,i2,i3)=ry(i1,i2,i3)*tyr4(i1,i2,i3)+sy(i1,i2,i3)*tys4(
     & i1,i2,i3)+ty(i1,i2,i3)*tyt4(i1,i2,i3)
       tyz43(i1,i2,i3)=rz(i1,i2,i3)*tyr4(i1,i2,i3)+sz(i1,i2,i3)*tys4(
     & i1,i2,i3)+tz(i1,i2,i3)*tyt4(i1,i2,i3)
       tzx42(i1,i2,i3)= rx(i1,i2,i3)*tzr4(i1,i2,i3)+sx(i1,i2,i3)*tzs4(
     & i1,i2,i3)
       tzy42(i1,i2,i3)= ry(i1,i2,i3)*tzr4(i1,i2,i3)+sy(i1,i2,i3)*tzs4(
     & i1,i2,i3)
       tzx43(i1,i2,i3)=rx(i1,i2,i3)*tzr4(i1,i2,i3)+sx(i1,i2,i3)*tzs4(
     & i1,i2,i3)+tx(i1,i2,i3)*tzt4(i1,i2,i3)
       tzy43(i1,i2,i3)=ry(i1,i2,i3)*tzr4(i1,i2,i3)+sy(i1,i2,i3)*tzs4(
     & i1,i2,i3)+ty(i1,i2,i3)*tzt4(i1,i2,i3)
       tzz43(i1,i2,i3)=rz(i1,i2,i3)*tzr4(i1,i2,i3)+sz(i1,i2,i3)*tzs4(
     & i1,i2,i3)+tz(i1,i2,i3)*tzt4(i1,i2,i3)
       uxx41(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*urr4(i1,i2,i3,kd)+(rxx42(
     & i1,i2,i3))*ur4(i1,i2,i3,kd)
       uyy41(i1,i2,i3,kd)=0
       uxy41(i1,i2,i3,kd)=0
       uxz41(i1,i2,i3,kd)=0
       uyz41(i1,i2,i3,kd)=0
       uzz41(i1,i2,i3,kd)=0
       ulaplacian41(i1,i2,i3,kd)=uxx41(i1,i2,i3,kd)
       uxx42(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*urr4(i1,i2,i3,kd)+2.*(rx(
     & i1,i2,i3)*sx(i1,i2,i3))*urs4(i1,i2,i3,kd)+(sx(i1,i2,i3)**2)*
     & uss4(i1,i2,i3,kd)+(rxx42(i1,i2,i3))*ur4(i1,i2,i3,kd)+(sxx42(i1,
     & i2,i3))*us4(i1,i2,i3,kd)
       uyy42(i1,i2,i3,kd)=(ry(i1,i2,i3)**2)*urr4(i1,i2,i3,kd)+2.*(ry(
     & i1,i2,i3)*sy(i1,i2,i3))*urs4(i1,i2,i3,kd)+(sy(i1,i2,i3)**2)*
     & uss4(i1,i2,i3,kd)+(ryy42(i1,i2,i3))*ur4(i1,i2,i3,kd)+(syy42(i1,
     & i2,i3))*us4(i1,i2,i3,kd)
       uxy42(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*urr4(i1,i2,i3,kd)+(
     & rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*urs4(i1,
     & i2,i3,kd)+sx(i1,i2,i3)*sy(i1,i2,i3)*uss4(i1,i2,i3,kd)+rxy42(i1,
     & i2,i3)*ur4(i1,i2,i3,kd)+sxy42(i1,i2,i3)*us4(i1,i2,i3,kd)
       uxz42(i1,i2,i3,kd)=0
       uyz42(i1,i2,i3,kd)=0
       uzz42(i1,i2,i3,kd)=0
       ulaplacian42(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2)*
     & urr4(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*
     & sy(i1,i2,i3))*urs4(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**
     & 2)*uss4(i1,i2,i3,kd)+(rxx42(i1,i2,i3)+ryy42(i1,i2,i3))*ur4(i1,
     & i2,i3,kd)+(sxx42(i1,i2,i3)+syy42(i1,i2,i3))*us4(i1,i2,i3,kd)
       uxx43(i1,i2,i3,kd)=rx(i1,i2,i3)**2*urr4(i1,i2,i3,kd)+sx(i1,i2,
     & i3)**2*uss4(i1,i2,i3,kd)+tx(i1,i2,i3)**2*utt4(i1,i2,i3,kd)+2.*
     & rx(i1,i2,i3)*sx(i1,i2,i3)*urs4(i1,i2,i3,kd)+2.*rx(i1,i2,i3)*tx(
     & i1,i2,i3)*urt4(i1,i2,i3,kd)+2.*sx(i1,i2,i3)*tx(i1,i2,i3)*ust4(
     & i1,i2,i3,kd)+rxx43(i1,i2,i3)*ur4(i1,i2,i3,kd)+sxx43(i1,i2,i3)*
     & us4(i1,i2,i3,kd)+txx43(i1,i2,i3)*ut4(i1,i2,i3,kd)
       uyy43(i1,i2,i3,kd)=ry(i1,i2,i3)**2*urr4(i1,i2,i3,kd)+sy(i1,i2,
     & i3)**2*uss4(i1,i2,i3,kd)+ty(i1,i2,i3)**2*utt4(i1,i2,i3,kd)+2.*
     & ry(i1,i2,i3)*sy(i1,i2,i3)*urs4(i1,i2,i3,kd)+2.*ry(i1,i2,i3)*ty(
     & i1,i2,i3)*urt4(i1,i2,i3,kd)+2.*sy(i1,i2,i3)*ty(i1,i2,i3)*ust4(
     & i1,i2,i3,kd)+ryy43(i1,i2,i3)*ur4(i1,i2,i3,kd)+syy43(i1,i2,i3)*
     & us4(i1,i2,i3,kd)+tyy43(i1,i2,i3)*ut4(i1,i2,i3,kd)
       uzz43(i1,i2,i3,kd)=rz(i1,i2,i3)**2*urr4(i1,i2,i3,kd)+sz(i1,i2,
     & i3)**2*uss4(i1,i2,i3,kd)+tz(i1,i2,i3)**2*utt4(i1,i2,i3,kd)+2.*
     & rz(i1,i2,i3)*sz(i1,i2,i3)*urs4(i1,i2,i3,kd)+2.*rz(i1,i2,i3)*tz(
     & i1,i2,i3)*urt4(i1,i2,i3,kd)+2.*sz(i1,i2,i3)*tz(i1,i2,i3)*ust4(
     & i1,i2,i3,kd)+rzz43(i1,i2,i3)*ur4(i1,i2,i3,kd)+szz43(i1,i2,i3)*
     & us4(i1,i2,i3,kd)+tzz43(i1,i2,i3)*ut4(i1,i2,i3,kd)
       uxy43(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*urr4(i1,i2,i3,kd)+
     & sx(i1,i2,i3)*sy(i1,i2,i3)*uss4(i1,i2,i3,kd)+tx(i1,i2,i3)*ty(i1,
     & i2,i3)*utt4(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,
     & i3)*sx(i1,i2,i3))*urs4(i1,i2,i3,kd)+(rx(i1,i2,i3)*ty(i1,i2,i3)+
     & ry(i1,i2,i3)*tx(i1,i2,i3))*urt4(i1,i2,i3,kd)+(sx(i1,i2,i3)*ty(
     & i1,i2,i3)+sy(i1,i2,i3)*tx(i1,i2,i3))*ust4(i1,i2,i3,kd)+rxy43(
     & i1,i2,i3)*ur4(i1,i2,i3,kd)+sxy43(i1,i2,i3)*us4(i1,i2,i3,kd)+
     & txy43(i1,i2,i3)*ut4(i1,i2,i3,kd)
       uxz43(i1,i2,i3,kd)=rx(i1,i2,i3)*rz(i1,i2,i3)*urr4(i1,i2,i3,kd)+
     & sx(i1,i2,i3)*sz(i1,i2,i3)*uss4(i1,i2,i3,kd)+tx(i1,i2,i3)*tz(i1,
     & i2,i3)*utt4(i1,i2,i3,kd)+(rx(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,i2,
     & i3)*sx(i1,i2,i3))*urs4(i1,i2,i3,kd)+(rx(i1,i2,i3)*tz(i1,i2,i3)+
     & rz(i1,i2,i3)*tx(i1,i2,i3))*urt4(i1,i2,i3,kd)+(sx(i1,i2,i3)*tz(
     & i1,i2,i3)+sz(i1,i2,i3)*tx(i1,i2,i3))*ust4(i1,i2,i3,kd)+rxz43(
     & i1,i2,i3)*ur4(i1,i2,i3,kd)+sxz43(i1,i2,i3)*us4(i1,i2,i3,kd)+
     & txz43(i1,i2,i3)*ut4(i1,i2,i3,kd)
       uyz43(i1,i2,i3,kd)=ry(i1,i2,i3)*rz(i1,i2,i3)*urr4(i1,i2,i3,kd)+
     & sy(i1,i2,i3)*sz(i1,i2,i3)*uss4(i1,i2,i3,kd)+ty(i1,i2,i3)*tz(i1,
     & i2,i3)*utt4(i1,i2,i3,kd)+(ry(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,i2,
     & i3)*sy(i1,i2,i3))*urs4(i1,i2,i3,kd)+(ry(i1,i2,i3)*tz(i1,i2,i3)+
     & rz(i1,i2,i3)*ty(i1,i2,i3))*urt4(i1,i2,i3,kd)+(sy(i1,i2,i3)*tz(
     & i1,i2,i3)+sz(i1,i2,i3)*ty(i1,i2,i3))*ust4(i1,i2,i3,kd)+ryz43(
     & i1,i2,i3)*ur4(i1,i2,i3,kd)+syz43(i1,i2,i3)*us4(i1,i2,i3,kd)+
     & tyz43(i1,i2,i3)*ut4(i1,i2,i3,kd)
       ulaplacian43(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2+rz(
     & i1,i2,i3)**2)*urr4(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**
     & 2+sz(i1,i2,i3)**2)*uss4(i1,i2,i3,kd)+(tx(i1,i2,i3)**2+ty(i1,i2,
     & i3)**2+tz(i1,i2,i3)**2)*utt4(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(
     & i1,i2,i3)+ ry(i1,i2,i3)*sy(i1,i2,i3)+rz(i1,i2,i3)*sz(i1,i2,i3))
     & *urs4(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*tx(i1,i2,i3)+ ry(i1,i2,i3)*
     & ty(i1,i2,i3)+rz(i1,i2,i3)*tz(i1,i2,i3))*urt4(i1,i2,i3,kd)+2.*(
     & sx(i1,i2,i3)*tx(i1,i2,i3)+ sy(i1,i2,i3)*ty(i1,i2,i3)+sz(i1,i2,
     & i3)*tz(i1,i2,i3))*ust4(i1,i2,i3,kd)+(rxx43(i1,i2,i3)+ryy43(i1,
     & i2,i3)+rzz43(i1,i2,i3))*ur4(i1,i2,i3,kd)+(sxx43(i1,i2,i3)+
     & syy43(i1,i2,i3)+szz43(i1,i2,i3))*us4(i1,i2,i3,kd)+(txx43(i1,i2,
     & i3)+tyy43(i1,i2,i3)+tzz43(i1,i2,i3))*ut4(i1,i2,i3,kd)
c============================================================================================
c Define derivatives for a rectangular grid
c
c============================================================================================
       h41(kd) = 1./(12.*dx(kd))
       h42(kd) = 1./(12.*dx(kd)**2)
       ux43r(i1,i2,i3,kd)=(8.*(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))-(u(
     & i1+2,i2,i3,kd)-u(i1-2,i2,i3,kd)))*h41(0)
       uy43r(i1,i2,i3,kd)=(8.*(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))-(u(
     & i1,i2+2,i3,kd)-u(i1,i2-2,i3,kd)))*h41(1)
       uz43r(i1,i2,i3,kd)=(8.*(u(i1,i2,i3+1,kd)-u(i1,i2,i3-1,kd))-(u(
     & i1,i2,i3+2,kd)-u(i1,i2,i3-2,kd)))*h41(2)
       uxx43r(i1,i2,i3,kd)=( -30.*u(i1,i2,i3,kd)+16.*(u(i1+1,i2,i3,kd)+
     & u(i1-1,i2,i3,kd))-(u(i1+2,i2,i3,kd)+u(i1-2,i2,i3,kd)) )*h42(0)
       uyy43r(i1,i2,i3,kd)=( -30.*u(i1,i2,i3,kd)+16.*(u(i1,i2+1,i3,kd)+
     & u(i1,i2-1,i3,kd))-(u(i1,i2+2,i3,kd)+u(i1,i2-2,i3,kd)) )*h42(1)
       uzz43r(i1,i2,i3,kd)=( -30.*u(i1,i2,i3,kd)+16.*(u(i1,i2,i3+1,kd)+
     & u(i1,i2,i3-1,kd))-(u(i1,i2,i3+2,kd)+u(i1,i2,i3-2,kd)) )*h42(2)
       uxy43r(i1,i2,i3,kd)=( (u(i1+2,i2+2,i3,kd)-u(i1-2,i2+2,i3,kd)- u(
     & i1+2,i2-2,i3,kd)+u(i1-2,i2-2,i3,kd)) +8.*(u(i1-1,i2+2,i3,kd)-u(
     & i1-1,i2-2,i3,kd)-u(i1+1,i2+2,i3,kd)+u(i1+1,i2-2,i3,kd) +u(i1+2,
     & i2-1,i3,kd)-u(i1-2,i2-1,i3,kd)-u(i1+2,i2+1,i3,kd)+u(i1-2,i2+1,
     & i3,kd))+64.*(u(i1+1,i2+1,i3,kd)-u(i1-1,i2+1,i3,kd)- u(i1+1,i2-
     & 1,i3,kd)+u(i1-1,i2-1,i3,kd)))*(h41(0)*h41(1))
       uxz43r(i1,i2,i3,kd)=( (u(i1+2,i2,i3+2,kd)-u(i1-2,i2,i3+2,kd)-u(
     & i1+2,i2,i3-2,kd)+u(i1-2,i2,i3-2,kd)) +8.*(u(i1-1,i2,i3+2,kd)-u(
     & i1-1,i2,i3-2,kd)-u(i1+1,i2,i3+2,kd)+u(i1+1,i2,i3-2,kd) +u(i1+2,
     & i2,i3-1,kd)-u(i1-2,i2,i3-1,kd)- u(i1+2,i2,i3+1,kd)+u(i1-2,i2,
     & i3+1,kd)) +64.*(u(i1+1,i2,i3+1,kd)-u(i1-1,i2,i3+1,kd)-u(i1+1,
     & i2,i3-1,kd)+u(i1-1,i2,i3-1,kd)) )*(h41(0)*h41(2))
       uyz43r(i1,i2,i3,kd)=( (u(i1,i2+2,i3+2,kd)-u(i1,i2-2,i3+2,kd)-u(
     & i1,i2+2,i3-2,kd)+u(i1,i2-2,i3-2,kd)) +8.*(u(i1,i2-1,i3+2,kd)-u(
     & i1,i2-1,i3-2,kd)-u(i1,i2+1,i3+2,kd)+u(i1,i2+1,i3-2,kd) +u(i1,
     & i2+2,i3-1,kd)-u(i1,i2-2,i3-1,kd)-u(i1,i2+2,i3+1,kd)+u(i1,i2-2,
     & i3+1,kd)) +64.*(u(i1,i2+1,i3+1,kd)-u(i1,i2-1,i3+1,kd)-u(i1,i2+
     & 1,i3-1,kd)+u(i1,i2-1,i3-1,kd)) )*(h41(1)*h41(2))
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
       ulaplacian42r(i1,i2,i3,kd)=uxx43r(i1,i2,i3,kd)+uyy43r(i1,i2,i3,
     & kd)
       ulaplacian43r(i1,i2,i3,kd)=uxx43r(i1,i2,i3,kd)+uyy43r(i1,i2,i3,
     & kd)+uzz43r(i1,i2,i3,kd)
       unr2(i1,i2,i3,kd)=(un(i1+1,i2,i3,kd)-un(i1-1,i2,i3,kd))*d12(0)
       uns2(i1,i2,i3,kd)=(un(i1,i2+1,i3,kd)-un(i1,i2-1,i3,kd))*d12(1)
       unt2(i1,i2,i3,kd)=(un(i1,i2,i3+1,kd)-un(i1,i2,i3-1,kd))*d12(2)
       unrr2(i1,i2,i3,kd)=(-2.*un(i1,i2,i3,kd)+(un(i1+1,i2,i3,kd)+un(
     & i1-1,i2,i3,kd)) )*d22(0)
       unss2(i1,i2,i3,kd)=(-2.*un(i1,i2,i3,kd)+(un(i1,i2+1,i3,kd)+un(
     & i1,i2-1,i3,kd)) )*d22(1)
       unrs2(i1,i2,i3,kd)=(unr2(i1,i2+1,i3,kd)-unr2(i1,i2-1,i3,kd))*
     & d12(1)
       untt2(i1,i2,i3,kd)=(-2.*un(i1,i2,i3,kd)+(un(i1,i2,i3+1,kd)+un(
     & i1,i2,i3-1,kd)) )*d22(2)
       unrt2(i1,i2,i3,kd)=(unr2(i1,i2,i3+1,kd)-unr2(i1,i2,i3-1,kd))*
     & d12(2)
       unst2(i1,i2,i3,kd)=(uns2(i1,i2,i3+1,kd)-uns2(i1,i2,i3-1,kd))*
     & d12(2)
       unrrr2(i1,i2,i3,kd)=(-2.*(un(i1+1,i2,i3,kd)-un(i1-1,i2,i3,kd))+(
     & un(i1+2,i2,i3,kd)-un(i1-2,i2,i3,kd)) )*d22(0)*d12(0)
       unsss2(i1,i2,i3,kd)=(-2.*(un(i1,i2+1,i3,kd)-un(i1,i2-1,i3,kd))+(
     & un(i1,i2+2,i3,kd)-un(i1,i2-2,i3,kd)) )*d22(1)*d12(1)
       unttt2(i1,i2,i3,kd)=(-2.*(un(i1,i2,i3+1,kd)-un(i1,i2,i3-1,kd))+(
     & un(i1,i2,i3+2,kd)-un(i1,i2,i3-2,kd)) )*d22(2)*d12(2)
       unx21(i1,i2,i3,kd)= rx(i1,i2,i3)*unr2(i1,i2,i3,kd)
       uny21(i1,i2,i3,kd)=0
       unz21(i1,i2,i3,kd)=0
       unx22(i1,i2,i3,kd)= rx(i1,i2,i3)*unr2(i1,i2,i3,kd)+sx(i1,i2,i3)*
     & uns2(i1,i2,i3,kd)
       uny22(i1,i2,i3,kd)= ry(i1,i2,i3)*unr2(i1,i2,i3,kd)+sy(i1,i2,i3)*
     & uns2(i1,i2,i3,kd)
       unz22(i1,i2,i3,kd)=0
       unx23(i1,i2,i3,kd)=rx(i1,i2,i3)*unr2(i1,i2,i3,kd)+sx(i1,i2,i3)*
     & uns2(i1,i2,i3,kd)+tx(i1,i2,i3)*unt2(i1,i2,i3,kd)
       uny23(i1,i2,i3,kd)=ry(i1,i2,i3)*unr2(i1,i2,i3,kd)+sy(i1,i2,i3)*
     & uns2(i1,i2,i3,kd)+ty(i1,i2,i3)*unt2(i1,i2,i3,kd)
       unz23(i1,i2,i3,kd)=rz(i1,i2,i3)*unr2(i1,i2,i3,kd)+sz(i1,i2,i3)*
     & uns2(i1,i2,i3,kd)+tz(i1,i2,i3)*unt2(i1,i2,i3,kd)
       unxx21(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*unrr2(i1,i2,i3,kd)+(rxx22(
     & i1,i2,i3))*unr2(i1,i2,i3,kd)
       unyy21(i1,i2,i3,kd)=0
       unxy21(i1,i2,i3,kd)=0
       unxz21(i1,i2,i3,kd)=0
       unyz21(i1,i2,i3,kd)=0
       unzz21(i1,i2,i3,kd)=0
       unlaplacian21(i1,i2,i3,kd)=unxx21(i1,i2,i3,kd)
       unxx22(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*unrr2(i1,i2,i3,kd)+2.*(rx(
     & i1,i2,i3)*sx(i1,i2,i3))*unrs2(i1,i2,i3,kd)+(sx(i1,i2,i3)**2)*
     & unss2(i1,i2,i3,kd)+(rxx22(i1,i2,i3))*unr2(i1,i2,i3,kd)+(sxx22(
     & i1,i2,i3))*uns2(i1,i2,i3,kd)
       unyy22(i1,i2,i3,kd)=(ry(i1,i2,i3)**2)*unrr2(i1,i2,i3,kd)+2.*(ry(
     & i1,i2,i3)*sy(i1,i2,i3))*unrs2(i1,i2,i3,kd)+(sy(i1,i2,i3)**2)*
     & unss2(i1,i2,i3,kd)+(ryy22(i1,i2,i3))*unr2(i1,i2,i3,kd)+(syy22(
     & i1,i2,i3))*uns2(i1,i2,i3,kd)
       unxy22(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*unrr2(i1,i2,i3,kd)
     & +(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*unrs2(
     & i1,i2,i3,kd)+sx(i1,i2,i3)*sy(i1,i2,i3)*unss2(i1,i2,i3,kd)+
     & rxy22(i1,i2,i3)*unr2(i1,i2,i3,kd)+sxy22(i1,i2,i3)*uns2(i1,i2,
     & i3,kd)
       unxz22(i1,i2,i3,kd)=0
       unyz22(i1,i2,i3,kd)=0
       unzz22(i1,i2,i3,kd)=0
       unlaplacian22(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2)*
     & unrr2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*
     & sy(i1,i2,i3))*unrs2(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)*
     & *2)*unss2(i1,i2,i3,kd)+(rxx22(i1,i2,i3)+ryy22(i1,i2,i3))*unr2(
     & i1,i2,i3,kd)+(sxx22(i1,i2,i3)+syy22(i1,i2,i3))*uns2(i1,i2,i3,
     & kd)
       unxx23(i1,i2,i3,kd)=rx(i1,i2,i3)**2*unrr2(i1,i2,i3,kd)+sx(i1,i2,
     & i3)**2*unss2(i1,i2,i3,kd)+tx(i1,i2,i3)**2*untt2(i1,i2,i3,kd)+
     & 2.*rx(i1,i2,i3)*sx(i1,i2,i3)*unrs2(i1,i2,i3,kd)+2.*rx(i1,i2,i3)
     & *tx(i1,i2,i3)*unrt2(i1,i2,i3,kd)+2.*sx(i1,i2,i3)*tx(i1,i2,i3)*
     & unst2(i1,i2,i3,kd)+rxx23(i1,i2,i3)*unr2(i1,i2,i3,kd)+sxx23(i1,
     & i2,i3)*uns2(i1,i2,i3,kd)+txx23(i1,i2,i3)*unt2(i1,i2,i3,kd)
       unyy23(i1,i2,i3,kd)=ry(i1,i2,i3)**2*unrr2(i1,i2,i3,kd)+sy(i1,i2,
     & i3)**2*unss2(i1,i2,i3,kd)+ty(i1,i2,i3)**2*untt2(i1,i2,i3,kd)+
     & 2.*ry(i1,i2,i3)*sy(i1,i2,i3)*unrs2(i1,i2,i3,kd)+2.*ry(i1,i2,i3)
     & *ty(i1,i2,i3)*unrt2(i1,i2,i3,kd)+2.*sy(i1,i2,i3)*ty(i1,i2,i3)*
     & unst2(i1,i2,i3,kd)+ryy23(i1,i2,i3)*unr2(i1,i2,i3,kd)+syy23(i1,
     & i2,i3)*uns2(i1,i2,i3,kd)+tyy23(i1,i2,i3)*unt2(i1,i2,i3,kd)
       unzz23(i1,i2,i3,kd)=rz(i1,i2,i3)**2*unrr2(i1,i2,i3,kd)+sz(i1,i2,
     & i3)**2*unss2(i1,i2,i3,kd)+tz(i1,i2,i3)**2*untt2(i1,i2,i3,kd)+
     & 2.*rz(i1,i2,i3)*sz(i1,i2,i3)*unrs2(i1,i2,i3,kd)+2.*rz(i1,i2,i3)
     & *tz(i1,i2,i3)*unrt2(i1,i2,i3,kd)+2.*sz(i1,i2,i3)*tz(i1,i2,i3)*
     & unst2(i1,i2,i3,kd)+rzz23(i1,i2,i3)*unr2(i1,i2,i3,kd)+szz23(i1,
     & i2,i3)*uns2(i1,i2,i3,kd)+tzz23(i1,i2,i3)*unt2(i1,i2,i3,kd)
       unxy23(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*unrr2(i1,i2,i3,kd)
     & +sx(i1,i2,i3)*sy(i1,i2,i3)*unss2(i1,i2,i3,kd)+tx(i1,i2,i3)*ty(
     & i1,i2,i3)*untt2(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,
     & i2,i3)*sx(i1,i2,i3))*unrs2(i1,i2,i3,kd)+(rx(i1,i2,i3)*ty(i1,i2,
     & i3)+ry(i1,i2,i3)*tx(i1,i2,i3))*unrt2(i1,i2,i3,kd)+(sx(i1,i2,i3)
     & *ty(i1,i2,i3)+sy(i1,i2,i3)*tx(i1,i2,i3))*unst2(i1,i2,i3,kd)+
     & rxy23(i1,i2,i3)*unr2(i1,i2,i3,kd)+sxy23(i1,i2,i3)*uns2(i1,i2,
     & i3,kd)+txy23(i1,i2,i3)*unt2(i1,i2,i3,kd)
       unxz23(i1,i2,i3,kd)=rx(i1,i2,i3)*rz(i1,i2,i3)*unrr2(i1,i2,i3,kd)
     & +sx(i1,i2,i3)*sz(i1,i2,i3)*unss2(i1,i2,i3,kd)+tx(i1,i2,i3)*tz(
     & i1,i2,i3)*untt2(i1,i2,i3,kd)+(rx(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,
     & i2,i3)*sx(i1,i2,i3))*unrs2(i1,i2,i3,kd)+(rx(i1,i2,i3)*tz(i1,i2,
     & i3)+rz(i1,i2,i3)*tx(i1,i2,i3))*unrt2(i1,i2,i3,kd)+(sx(i1,i2,i3)
     & *tz(i1,i2,i3)+sz(i1,i2,i3)*tx(i1,i2,i3))*unst2(i1,i2,i3,kd)+
     & rxz23(i1,i2,i3)*unr2(i1,i2,i3,kd)+sxz23(i1,i2,i3)*uns2(i1,i2,
     & i3,kd)+txz23(i1,i2,i3)*unt2(i1,i2,i3,kd)
       unyz23(i1,i2,i3,kd)=ry(i1,i2,i3)*rz(i1,i2,i3)*unrr2(i1,i2,i3,kd)
     & +sy(i1,i2,i3)*sz(i1,i2,i3)*unss2(i1,i2,i3,kd)+ty(i1,i2,i3)*tz(
     & i1,i2,i3)*untt2(i1,i2,i3,kd)+(ry(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,
     & i2,i3)*sy(i1,i2,i3))*unrs2(i1,i2,i3,kd)+(ry(i1,i2,i3)*tz(i1,i2,
     & i3)+rz(i1,i2,i3)*ty(i1,i2,i3))*unrt2(i1,i2,i3,kd)+(sy(i1,i2,i3)
     & *tz(i1,i2,i3)+sz(i1,i2,i3)*ty(i1,i2,i3))*unst2(i1,i2,i3,kd)+
     & ryz23(i1,i2,i3)*unr2(i1,i2,i3,kd)+syz23(i1,i2,i3)*uns2(i1,i2,
     & i3,kd)+tyz23(i1,i2,i3)*unt2(i1,i2,i3,kd)
       unlaplacian23(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2+rz(
     & i1,i2,i3)**2)*unrr2(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)*
     & *2+sz(i1,i2,i3)**2)*unss2(i1,i2,i3,kd)+(tx(i1,i2,i3)**2+ty(i1,
     & i2,i3)**2+tz(i1,i2,i3)**2)*untt2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*
     & sx(i1,i2,i3)+ ry(i1,i2,i3)*sy(i1,i2,i3)+rz(i1,i2,i3)*sz(i1,i2,
     & i3))*unrs2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*tx(i1,i2,i3)+ ry(i1,
     & i2,i3)*ty(i1,i2,i3)+rz(i1,i2,i3)*tz(i1,i2,i3))*unrt2(i1,i2,i3,
     & kd)+2.*(sx(i1,i2,i3)*tx(i1,i2,i3)+ sy(i1,i2,i3)*ty(i1,i2,i3)+
     & sz(i1,i2,i3)*tz(i1,i2,i3))*unst2(i1,i2,i3,kd)+(rxx23(i1,i2,i3)+
     & ryy23(i1,i2,i3)+rzz23(i1,i2,i3))*unr2(i1,i2,i3,kd)+(sxx23(i1,
     & i2,i3)+syy23(i1,i2,i3)+szz23(i1,i2,i3))*uns2(i1,i2,i3,kd)+(
     & txx23(i1,i2,i3)+tyy23(i1,i2,i3)+tzz23(i1,i2,i3))*unt2(i1,i2,i3,
     & kd)
c============================================================================================
c Define derivatives for a rectangular grid
c
c============================================================================================
       unx23r(i1,i2,i3,kd)=(un(i1+1,i2,i3,kd)-un(i1-1,i2,i3,kd))*h12(0)
       uny23r(i1,i2,i3,kd)=(un(i1,i2+1,i3,kd)-un(i1,i2-1,i3,kd))*h12(1)
       unz23r(i1,i2,i3,kd)=(un(i1,i2,i3+1,kd)-un(i1,i2,i3-1,kd))*h12(2)
       unxx23r(i1,i2,i3,kd)=(-2.*un(i1,i2,i3,kd)+(un(i1+1,i2,i3,kd)+un(
     & i1-1,i2,i3,kd)) )*h22(0)
       unyy23r(i1,i2,i3,kd)=(-2.*un(i1,i2,i3,kd)+(un(i1,i2+1,i3,kd)+un(
     & i1,i2-1,i3,kd)) )*h22(1)
       unxy23r(i1,i2,i3,kd)=(unx23r(i1,i2+1,i3,kd)-unx23r(i1,i2-1,i3,
     & kd))*h12(1)
       unzz23r(i1,i2,i3,kd)=(-2.*un(i1,i2,i3,kd)+(un(i1,i2,i3+1,kd)+un(
     & i1,i2,i3-1,kd)) )*h22(2)
       unxz23r(i1,i2,i3,kd)=(unx23r(i1,i2,i3+1,kd)-unx23r(i1,i2,i3-1,
     & kd))*h12(2)
       unyz23r(i1,i2,i3,kd)=(uny23r(i1,i2,i3+1,kd)-uny23r(i1,i2,i3-1,
     & kd))*h12(2)
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
       unlaplacian22r(i1,i2,i3,kd)=unxx23r(i1,i2,i3,kd)+unyy23r(i1,i2,
     & i3,kd)
       unlaplacian23r(i1,i2,i3,kd)=unxx23r(i1,i2,i3,kd)+unyy23r(i1,i2,
     & i3,kd)+unzz23r(i1,i2,i3,kd)
       unxxx22r(i1,i2,i3,kd)=(-2.*(un(i1+1,i2,i3,kd)-un(i1-1,i2,i3,kd))
     & +(un(i1+2,i2,i3,kd)-un(i1-2,i2,i3,kd)) )*h22(0)*h12(0)
       unyyy22r(i1,i2,i3,kd)=(-2.*(un(i1,i2+1,i3,kd)-un(i1,i2-1,i3,kd))
     & +(un(i1,i2+2,i3,kd)-un(i1,i2-2,i3,kd)) )*h22(1)*h12(1)
       unxxy22r(i1,i2,i3,kd)=( unxx22r(i1,i2+1,i3,kd)-unxx22r(i1,i2-1,
     & i3,kd))/(2.*dx(1))
       unxyy22r(i1,i2,i3,kd)=( unyy22r(i1+1,i2,i3,kd)-unyy22r(i1-1,i2,
     & i3,kd))/(2.*dx(0))
       unxxxx22r(i1,i2,i3,kd)=(6.*un(i1,i2,i3,kd)-4.*(un(i1+1,i2,i3,kd)
     & +un(i1-1,i2,i3,kd))+(un(i1+2,i2,i3,kd)+un(i1-2,i2,i3,kd)) )/(
     & dx(0)**4)
       unyyyy22r(i1,i2,i3,kd)=(6.*un(i1,i2,i3,kd)-4.*(un(i1,i2+1,i3,kd)
     & +un(i1,i2-1,i3,kd))+(un(i1,i2+2,i3,kd)+un(i1,i2-2,i3,kd)) )/(
     & dx(1)**4)
       unxxyy22r(i1,i2,i3,kd)=( 4.*un(i1,i2,i3,kd)     -2.*(un(i1+1,i2,
     & i3,kd)+un(i1-1,i2,i3,kd)+un(i1,i2+1,i3,kd)+un(i1,i2-1,i3,kd))  
     &  +   (un(i1+1,i2+1,i3,kd)+un(i1-1,i2+1,i3,kd)+un(i1+1,i2-1,i3,
     & kd)+un(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)
       ! 2D laplacian squared = un.xxxx + 2 un.xxyy + un.yyyy
       unLapSq22r(i1,i2,i3,kd)= ( 6.*un(i1,i2,i3,kd)   - 4.*(un(i1+1,
     & i2,i3,kd)+un(i1-1,i2,i3,kd))    +(un(i1+2,i2,i3,kd)+un(i1-2,i2,
     & i3,kd)) )/(dx(0)**4) +( 6.*un(i1,i2,i3,kd)    -4.*(un(i1,i2+1,
     & i3,kd)+un(i1,i2-1,i3,kd))    +(un(i1,i2+2,i3,kd)+un(i1,i2-2,i3,
     & kd)) )/(dx(1)**4)  +( 8.*un(i1,i2,i3,kd)     -4.*(un(i1+1,i2,
     & i3,kd)+un(i1-1,i2,i3,kd)+un(i1,i2+1,i3,kd)+un(i1,i2-1,i3,kd))  
     &  +2.*(un(i1+1,i2+1,i3,kd)+un(i1-1,i2+1,i3,kd)+un(i1+1,i2-1,i3,
     & kd)+un(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)
       unxxx23r(i1,i2,i3,kd)=(-2.*(un(i1+1,i2,i3,kd)-un(i1-1,i2,i3,kd))
     & +(un(i1+2,i2,i3,kd)-un(i1-2,i2,i3,kd)) )*h22(0)*h12(0)
       unyyy23r(i1,i2,i3,kd)=(-2.*(un(i1,i2+1,i3,kd)-un(i1,i2-1,i3,kd))
     & +(un(i1,i2+2,i3,kd)-un(i1,i2-2,i3,kd)) )*h22(1)*h12(1)
       unzzz23r(i1,i2,i3,kd)=(-2.*(un(i1,i2,i3+1,kd)-un(i1,i2,i3-1,kd))
     & +(un(i1,i2,i3+2,kd)-un(i1,i2,i3-2,kd)) )*h22(1)*h12(2)
       unxxy23r(i1,i2,i3,kd)=( unxx22r(i1,i2+1,i3,kd)-unxx22r(i1,i2-1,
     & i3,kd))/(2.*dx(1))
       unxyy23r(i1,i2,i3,kd)=( unyy22r(i1+1,i2,i3,kd)-unyy22r(i1-1,i2,
     & i3,kd))/(2.*dx(0))
       unxxz23r(i1,i2,i3,kd)=( unxx22r(i1,i2,i3+1,kd)-unxx22r(i1,i2,i3-
     & 1,kd))/(2.*dx(2))
       unyyz23r(i1,i2,i3,kd)=( unyy22r(i1,i2,i3+1,kd)-unyy22r(i1,i2,i3-
     & 1,kd))/(2.*dx(2))
       unxzz23r(i1,i2,i3,kd)=( unzz22r(i1+1,i2,i3,kd)-unzz22r(i1-1,i2,
     & i3,kd))/(2.*dx(0))
       unyzz23r(i1,i2,i3,kd)=( unzz22r(i1,i2+1,i3,kd)-unzz22r(i1,i2-1,
     & i3,kd))/(2.*dx(1))
       unxxxx23r(i1,i2,i3,kd)=(6.*un(i1,i2,i3,kd)-4.*(un(i1+1,i2,i3,kd)
     & +un(i1-1,i2,i3,kd))+(un(i1+2,i2,i3,kd)+un(i1-2,i2,i3,kd)) )/(
     & dx(0)**4)
       unyyyy23r(i1,i2,i3,kd)=(6.*un(i1,i2,i3,kd)-4.*(un(i1,i2+1,i3,kd)
     & +un(i1,i2-1,i3,kd))+(un(i1,i2+2,i3,kd)+un(i1,i2-2,i3,kd)) )/(
     & dx(1)**4)
       unzzzz23r(i1,i2,i3,kd)=(6.*un(i1,i2,i3,kd)-4.*(un(i1,i2,i3+1,kd)
     & +un(i1,i2,i3-1,kd))+(un(i1,i2,i3+2,kd)+un(i1,i2,i3-2,kd)) )/(
     & dx(2)**4)
       unxxyy23r(i1,i2,i3,kd)=( 4.*un(i1,i2,i3,kd)     -2.*(un(i1+1,i2,
     & i3,kd)+un(i1-1,i2,i3,kd)+un(i1,i2+1,i3,kd)+un(i1,i2-1,i3,kd))  
     &  +   (un(i1+1,i2+1,i3,kd)+un(i1-1,i2+1,i3,kd)+un(i1+1,i2-1,i3,
     & kd)+un(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)
       unxxzz23r(i1,i2,i3,kd)=( 4.*un(i1,i2,i3,kd)     -2.*(un(i1+1,i2,
     & i3,kd)+un(i1-1,i2,i3,kd)+un(i1,i2,i3+1,kd)+un(i1,i2,i3-1,kd))  
     &  +   (un(i1+1,i2,i3+1,kd)+un(i1-1,i2,i3+1,kd)+un(i1+1,i2,i3-1,
     & kd)+un(i1-1,i2,i3-1,kd)) )/(dx(0)**2*dx(2)**2)
       unyyzz23r(i1,i2,i3,kd)=( 4.*un(i1,i2,i3,kd)     -2.*(un(i1,i2+1,
     & i3,kd)  +un(i1,i2-1,i3,kd)+  un(i1,i2  ,i3+1,kd)+un(i1,i2  ,i3-
     & 1,kd))   +   (un(i1,i2+1,i3+1,kd)+un(i1,i2-1,i3+1,kd)+un(i1,i2+
     & 1,i3-1,kd)+un(i1,i2-1,i3-1,kd)) )/(dx(1)**2*dx(2)**2)
       ! 3D laplacian squared = un.xxxx + un.yyyy + un.zzzz + 2 (un.xxyy + un.xxzz + un.yyzz )
       unLapSq23r(i1,i2,i3,kd)= ( 6.*un(i1,i2,i3,kd)   - 4.*(un(i1+1,
     & i2,i3,kd)+un(i1-1,i2,i3,kd))    +(un(i1+2,i2,i3,kd)+un(i1-2,i2,
     & i3,kd)) )/(dx(0)**4) +( 6.*un(i1,i2,i3,kd)    -4.*(un(i1,i2+1,
     & i3,kd)+un(i1,i2-1,i3,kd))    +(un(i1,i2+2,i3,kd)+un(i1,i2-2,i3,
     & kd)) )/(dx(1)**4)  +( 6.*un(i1,i2,i3,kd)    -4.*(un(i1,i2,i3+1,
     & kd)+un(i1,i2,i3-1,kd))    +(un(i1,i2,i3+2,kd)+un(i1,i2,i3-2,kd)
     & ) )/(dx(2)**4)  +( 8.*un(i1,i2,i3,kd)     -4.*(un(i1+1,i2,i3,
     & kd)  +un(i1-1,i2,i3,kd)  +un(i1  ,i2+1,i3,kd)+un(i1  ,i2-1,i3,
     & kd))   +2.*(un(i1+1,i2+1,i3,kd)+un(i1-1,i2+1,i3,kd)+un(i1+1,i2-
     & 1,i3,kd)+un(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)+( 8.*un(i1,
     & i2,i3,kd)     -4.*(un(i1+1,i2,i3,kd)  +un(i1-1,i2,i3,kd)  +un(
     & i1  ,i2,i3+1,kd)+un(i1  ,i2,i3-1,kd))   +2.*(un(i1+1,i2,i3+1,
     & kd)+un(i1-1,i2,i3+1,kd)+un(i1+1,i2,i3-1,kd)+un(i1-1,i2,i3-1,kd)
     & ) )/(dx(0)**2*dx(2)**2)+( 8.*un(i1,i2,i3,kd)     -4.*(un(i1,i2+
     & 1,i3,kd)  +un(i1,i2-1,i3,kd)  +un(i1,i2  ,i3+1,kd)+un(i1,i2  ,
     & i3-1,kd))   +2.*(un(i1,i2+1,i3+1,kd)+un(i1,i2-1,i3+1,kd)+un(i1,
     & i2+1,i3-1,kd)+un(i1,i2-1,i3-1,kd)) )/(dx(1)**2*dx(2)**2)
       unr4(i1,i2,i3,kd)=(8.*(un(i1+1,i2,i3,kd)-un(i1-1,i2,i3,kd))-(un(
     & i1+2,i2,i3,kd)-un(i1-2,i2,i3,kd)))*d14(0)
       uns4(i1,i2,i3,kd)=(8.*(un(i1,i2+1,i3,kd)-un(i1,i2-1,i3,kd))-(un(
     & i1,i2+2,i3,kd)-un(i1,i2-2,i3,kd)))*d14(1)
       unt4(i1,i2,i3,kd)=(8.*(un(i1,i2,i3+1,kd)-un(i1,i2,i3-1,kd))-(un(
     & i1,i2,i3+2,kd)-un(i1,i2,i3-2,kd)))*d14(2)
       unrr4(i1,i2,i3,kd)=(-30.*un(i1,i2,i3,kd)+16.*(un(i1+1,i2,i3,kd)+
     & un(i1-1,i2,i3,kd))-(un(i1+2,i2,i3,kd)+un(i1-2,i2,i3,kd)) )*d24(
     & 0)
       unss4(i1,i2,i3,kd)=(-30.*un(i1,i2,i3,kd)+16.*(un(i1,i2+1,i3,kd)+
     & un(i1,i2-1,i3,kd))-(un(i1,i2+2,i3,kd)+un(i1,i2-2,i3,kd)) )*d24(
     & 1)
       untt4(i1,i2,i3,kd)=(-30.*un(i1,i2,i3,kd)+16.*(un(i1,i2,i3+1,kd)+
     & un(i1,i2,i3-1,kd))-(un(i1,i2,i3+2,kd)+un(i1,i2,i3-2,kd)) )*d24(
     & 2)
       unrs4(i1,i2,i3,kd)=(8.*(unr4(i1,i2+1,i3,kd)-unr4(i1,i2-1,i3,kd))
     & -(unr4(i1,i2+2,i3,kd)-unr4(i1,i2-2,i3,kd)))*d14(1)
       unrt4(i1,i2,i3,kd)=(8.*(unr4(i1,i2,i3+1,kd)-unr4(i1,i2,i3-1,kd))
     & -(unr4(i1,i2,i3+2,kd)-unr4(i1,i2,i3-2,kd)))*d14(2)
       unst4(i1,i2,i3,kd)=(8.*(uns4(i1,i2,i3+1,kd)-uns4(i1,i2,i3-1,kd))
     & -(uns4(i1,i2,i3+2,kd)-uns4(i1,i2,i3-2,kd)))*d14(2)
       unx41(i1,i2,i3,kd)= rx(i1,i2,i3)*unr4(i1,i2,i3,kd)
       uny41(i1,i2,i3,kd)=0
       unz41(i1,i2,i3,kd)=0
       unx42(i1,i2,i3,kd)= rx(i1,i2,i3)*unr4(i1,i2,i3,kd)+sx(i1,i2,i3)*
     & uns4(i1,i2,i3,kd)
       uny42(i1,i2,i3,kd)= ry(i1,i2,i3)*unr4(i1,i2,i3,kd)+sy(i1,i2,i3)*
     & uns4(i1,i2,i3,kd)
       unz42(i1,i2,i3,kd)=0
       unx43(i1,i2,i3,kd)=rx(i1,i2,i3)*unr4(i1,i2,i3,kd)+sx(i1,i2,i3)*
     & uns4(i1,i2,i3,kd)+tx(i1,i2,i3)*unt4(i1,i2,i3,kd)
       uny43(i1,i2,i3,kd)=ry(i1,i2,i3)*unr4(i1,i2,i3,kd)+sy(i1,i2,i3)*
     & uns4(i1,i2,i3,kd)+ty(i1,i2,i3)*unt4(i1,i2,i3,kd)
       unz43(i1,i2,i3,kd)=rz(i1,i2,i3)*unr4(i1,i2,i3,kd)+sz(i1,i2,i3)*
     & uns4(i1,i2,i3,kd)+tz(i1,i2,i3)*unt4(i1,i2,i3,kd)
       unxx41(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*unrr4(i1,i2,i3,kd)+(rxx42(
     & i1,i2,i3))*unr4(i1,i2,i3,kd)
       unyy41(i1,i2,i3,kd)=0
       unxy41(i1,i2,i3,kd)=0
       unxz41(i1,i2,i3,kd)=0
       unyz41(i1,i2,i3,kd)=0
       unzz41(i1,i2,i3,kd)=0
       unlaplacian41(i1,i2,i3,kd)=unxx41(i1,i2,i3,kd)
       unxx42(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*unrr4(i1,i2,i3,kd)+2.*(rx(
     & i1,i2,i3)*sx(i1,i2,i3))*unrs4(i1,i2,i3,kd)+(sx(i1,i2,i3)**2)*
     & unss4(i1,i2,i3,kd)+(rxx42(i1,i2,i3))*unr4(i1,i2,i3,kd)+(sxx42(
     & i1,i2,i3))*uns4(i1,i2,i3,kd)
       unyy42(i1,i2,i3,kd)=(ry(i1,i2,i3)**2)*unrr4(i1,i2,i3,kd)+2.*(ry(
     & i1,i2,i3)*sy(i1,i2,i3))*unrs4(i1,i2,i3,kd)+(sy(i1,i2,i3)**2)*
     & unss4(i1,i2,i3,kd)+(ryy42(i1,i2,i3))*unr4(i1,i2,i3,kd)+(syy42(
     & i1,i2,i3))*uns4(i1,i2,i3,kd)
       unxy42(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*unrr4(i1,i2,i3,kd)
     & +(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*unrs4(
     & i1,i2,i3,kd)+sx(i1,i2,i3)*sy(i1,i2,i3)*unss4(i1,i2,i3,kd)+
     & rxy42(i1,i2,i3)*unr4(i1,i2,i3,kd)+sxy42(i1,i2,i3)*uns4(i1,i2,
     & i3,kd)
       unxz42(i1,i2,i3,kd)=0
       unyz42(i1,i2,i3,kd)=0
       unzz42(i1,i2,i3,kd)=0
       unlaplacian42(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2)*
     & unrr4(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*
     & sy(i1,i2,i3))*unrs4(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)*
     & *2)*unss4(i1,i2,i3,kd)+(rxx42(i1,i2,i3)+ryy42(i1,i2,i3))*unr4(
     & i1,i2,i3,kd)+(sxx42(i1,i2,i3)+syy42(i1,i2,i3))*uns4(i1,i2,i3,
     & kd)
       unxx43(i1,i2,i3,kd)=rx(i1,i2,i3)**2*unrr4(i1,i2,i3,kd)+sx(i1,i2,
     & i3)**2*unss4(i1,i2,i3,kd)+tx(i1,i2,i3)**2*untt4(i1,i2,i3,kd)+
     & 2.*rx(i1,i2,i3)*sx(i1,i2,i3)*unrs4(i1,i2,i3,kd)+2.*rx(i1,i2,i3)
     & *tx(i1,i2,i3)*unrt4(i1,i2,i3,kd)+2.*sx(i1,i2,i3)*tx(i1,i2,i3)*
     & unst4(i1,i2,i3,kd)+rxx43(i1,i2,i3)*unr4(i1,i2,i3,kd)+sxx43(i1,
     & i2,i3)*uns4(i1,i2,i3,kd)+txx43(i1,i2,i3)*unt4(i1,i2,i3,kd)
       unyy43(i1,i2,i3,kd)=ry(i1,i2,i3)**2*unrr4(i1,i2,i3,kd)+sy(i1,i2,
     & i3)**2*unss4(i1,i2,i3,kd)+ty(i1,i2,i3)**2*untt4(i1,i2,i3,kd)+
     & 2.*ry(i1,i2,i3)*sy(i1,i2,i3)*unrs4(i1,i2,i3,kd)+2.*ry(i1,i2,i3)
     & *ty(i1,i2,i3)*unrt4(i1,i2,i3,kd)+2.*sy(i1,i2,i3)*ty(i1,i2,i3)*
     & unst4(i1,i2,i3,kd)+ryy43(i1,i2,i3)*unr4(i1,i2,i3,kd)+syy43(i1,
     & i2,i3)*uns4(i1,i2,i3,kd)+tyy43(i1,i2,i3)*unt4(i1,i2,i3,kd)
       unzz43(i1,i2,i3,kd)=rz(i1,i2,i3)**2*unrr4(i1,i2,i3,kd)+sz(i1,i2,
     & i3)**2*unss4(i1,i2,i3,kd)+tz(i1,i2,i3)**2*untt4(i1,i2,i3,kd)+
     & 2.*rz(i1,i2,i3)*sz(i1,i2,i3)*unrs4(i1,i2,i3,kd)+2.*rz(i1,i2,i3)
     & *tz(i1,i2,i3)*unrt4(i1,i2,i3,kd)+2.*sz(i1,i2,i3)*tz(i1,i2,i3)*
     & unst4(i1,i2,i3,kd)+rzz43(i1,i2,i3)*unr4(i1,i2,i3,kd)+szz43(i1,
     & i2,i3)*uns4(i1,i2,i3,kd)+tzz43(i1,i2,i3)*unt4(i1,i2,i3,kd)
       unxy43(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*unrr4(i1,i2,i3,kd)
     & +sx(i1,i2,i3)*sy(i1,i2,i3)*unss4(i1,i2,i3,kd)+tx(i1,i2,i3)*ty(
     & i1,i2,i3)*untt4(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,
     & i2,i3)*sx(i1,i2,i3))*unrs4(i1,i2,i3,kd)+(rx(i1,i2,i3)*ty(i1,i2,
     & i3)+ry(i1,i2,i3)*tx(i1,i2,i3))*unrt4(i1,i2,i3,kd)+(sx(i1,i2,i3)
     & *ty(i1,i2,i3)+sy(i1,i2,i3)*tx(i1,i2,i3))*unst4(i1,i2,i3,kd)+
     & rxy43(i1,i2,i3)*unr4(i1,i2,i3,kd)+sxy43(i1,i2,i3)*uns4(i1,i2,
     & i3,kd)+txy43(i1,i2,i3)*unt4(i1,i2,i3,kd)
       unxz43(i1,i2,i3,kd)=rx(i1,i2,i3)*rz(i1,i2,i3)*unrr4(i1,i2,i3,kd)
     & +sx(i1,i2,i3)*sz(i1,i2,i3)*unss4(i1,i2,i3,kd)+tx(i1,i2,i3)*tz(
     & i1,i2,i3)*untt4(i1,i2,i3,kd)+(rx(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,
     & i2,i3)*sx(i1,i2,i3))*unrs4(i1,i2,i3,kd)+(rx(i1,i2,i3)*tz(i1,i2,
     & i3)+rz(i1,i2,i3)*tx(i1,i2,i3))*unrt4(i1,i2,i3,kd)+(sx(i1,i2,i3)
     & *tz(i1,i2,i3)+sz(i1,i2,i3)*tx(i1,i2,i3))*unst4(i1,i2,i3,kd)+
     & rxz43(i1,i2,i3)*unr4(i1,i2,i3,kd)+sxz43(i1,i2,i3)*uns4(i1,i2,
     & i3,kd)+txz43(i1,i2,i3)*unt4(i1,i2,i3,kd)
       unyz43(i1,i2,i3,kd)=ry(i1,i2,i3)*rz(i1,i2,i3)*unrr4(i1,i2,i3,kd)
     & +sy(i1,i2,i3)*sz(i1,i2,i3)*unss4(i1,i2,i3,kd)+ty(i1,i2,i3)*tz(
     & i1,i2,i3)*untt4(i1,i2,i3,kd)+(ry(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,
     & i2,i3)*sy(i1,i2,i3))*unrs4(i1,i2,i3,kd)+(ry(i1,i2,i3)*tz(i1,i2,
     & i3)+rz(i1,i2,i3)*ty(i1,i2,i3))*unrt4(i1,i2,i3,kd)+(sy(i1,i2,i3)
     & *tz(i1,i2,i3)+sz(i1,i2,i3)*ty(i1,i2,i3))*unst4(i1,i2,i3,kd)+
     & ryz43(i1,i2,i3)*unr4(i1,i2,i3,kd)+syz43(i1,i2,i3)*uns4(i1,i2,
     & i3,kd)+tyz43(i1,i2,i3)*unt4(i1,i2,i3,kd)
       unlaplacian43(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2+rz(
     & i1,i2,i3)**2)*unrr4(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)*
     & *2+sz(i1,i2,i3)**2)*unss4(i1,i2,i3,kd)+(tx(i1,i2,i3)**2+ty(i1,
     & i2,i3)**2+tz(i1,i2,i3)**2)*untt4(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*
     & sx(i1,i2,i3)+ ry(i1,i2,i3)*sy(i1,i2,i3)+rz(i1,i2,i3)*sz(i1,i2,
     & i3))*unrs4(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*tx(i1,i2,i3)+ ry(i1,
     & i2,i3)*ty(i1,i2,i3)+rz(i1,i2,i3)*tz(i1,i2,i3))*unrt4(i1,i2,i3,
     & kd)+2.*(sx(i1,i2,i3)*tx(i1,i2,i3)+ sy(i1,i2,i3)*ty(i1,i2,i3)+
     & sz(i1,i2,i3)*tz(i1,i2,i3))*unst4(i1,i2,i3,kd)+(rxx43(i1,i2,i3)+
     & ryy43(i1,i2,i3)+rzz43(i1,i2,i3))*unr4(i1,i2,i3,kd)+(sxx43(i1,
     & i2,i3)+syy43(i1,i2,i3)+szz43(i1,i2,i3))*uns4(i1,i2,i3,kd)+(
     & txx43(i1,i2,i3)+tyy43(i1,i2,i3)+tzz43(i1,i2,i3))*unt4(i1,i2,i3,
     & kd)
c============================================================================================
c Define derivatives for a rectangular grid
c
c============================================================================================
       unx43r(i1,i2,i3,kd)=(8.*(un(i1+1,i2,i3,kd)-un(i1-1,i2,i3,kd))-(
     & un(i1+2,i2,i3,kd)-un(i1-2,i2,i3,kd)))*h41(0)
       uny43r(i1,i2,i3,kd)=(8.*(un(i1,i2+1,i3,kd)-un(i1,i2-1,i3,kd))-(
     & un(i1,i2+2,i3,kd)-un(i1,i2-2,i3,kd)))*h41(1)
       unz43r(i1,i2,i3,kd)=(8.*(un(i1,i2,i3+1,kd)-un(i1,i2,i3-1,kd))-(
     & un(i1,i2,i3+2,kd)-un(i1,i2,i3-2,kd)))*h41(2)
       unxx43r(i1,i2,i3,kd)=( -30.*un(i1,i2,i3,kd)+16.*(un(i1+1,i2,i3,
     & kd)+un(i1-1,i2,i3,kd))-(un(i1+2,i2,i3,kd)+un(i1-2,i2,i3,kd)) )*
     & h42(0)
       unyy43r(i1,i2,i3,kd)=( -30.*un(i1,i2,i3,kd)+16.*(un(i1,i2+1,i3,
     & kd)+un(i1,i2-1,i3,kd))-(un(i1,i2+2,i3,kd)+un(i1,i2-2,i3,kd)) )*
     & h42(1)
       unzz43r(i1,i2,i3,kd)=( -30.*un(i1,i2,i3,kd)+16.*(un(i1,i2,i3+1,
     & kd)+un(i1,i2,i3-1,kd))-(un(i1,i2,i3+2,kd)+un(i1,i2,i3-2,kd)) )*
     & h42(2)
       unxy43r(i1,i2,i3,kd)=( (un(i1+2,i2+2,i3,kd)-un(i1-2,i2+2,i3,kd)-
     &  un(i1+2,i2-2,i3,kd)+un(i1-2,i2-2,i3,kd)) +8.*(un(i1-1,i2+2,i3,
     & kd)-un(i1-1,i2-2,i3,kd)-un(i1+1,i2+2,i3,kd)+un(i1+1,i2-2,i3,kd)
     &  +un(i1+2,i2-1,i3,kd)-un(i1-2,i2-1,i3,kd)-un(i1+2,i2+1,i3,kd)+
     & un(i1-2,i2+1,i3,kd))+64.*(un(i1+1,i2+1,i3,kd)-un(i1-1,i2+1,i3,
     & kd)- un(i1+1,i2-1,i3,kd)+un(i1-1,i2-1,i3,kd)))*(h41(0)*h41(1))
       unxz43r(i1,i2,i3,kd)=( (un(i1+2,i2,i3+2,kd)-un(i1-2,i2,i3+2,kd)-
     & un(i1+2,i2,i3-2,kd)+un(i1-2,i2,i3-2,kd)) +8.*(un(i1-1,i2,i3+2,
     & kd)-un(i1-1,i2,i3-2,kd)-un(i1+1,i2,i3+2,kd)+un(i1+1,i2,i3-2,kd)
     &  +un(i1+2,i2,i3-1,kd)-un(i1-2,i2,i3-1,kd)- un(i1+2,i2,i3+1,kd)+
     & un(i1-2,i2,i3+1,kd)) +64.*(un(i1+1,i2,i3+1,kd)-un(i1-1,i2,i3+1,
     & kd)-un(i1+1,i2,i3-1,kd)+un(i1-1,i2,i3-1,kd)) )*(h41(0)*h41(2))
       unyz43r(i1,i2,i3,kd)=( (un(i1,i2+2,i3+2,kd)-un(i1,i2-2,i3+2,kd)-
     & un(i1,i2+2,i3-2,kd)+un(i1,i2-2,i3-2,kd)) +8.*(un(i1,i2-1,i3+2,
     & kd)-un(i1,i2-1,i3-2,kd)-un(i1,i2+1,i3+2,kd)+un(i1,i2+1,i3-2,kd)
     &  +un(i1,i2+2,i3-1,kd)-un(i1,i2-2,i3-1,kd)-un(i1,i2+2,i3+1,kd)+
     & un(i1,i2-2,i3+1,kd)) +64.*(un(i1,i2+1,i3+1,kd)-un(i1,i2-1,i3+1,
     & kd)-un(i1,i2+1,i3-1,kd)+un(i1,i2-1,i3-1,kd)) )*(h41(1)*h41(2))
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
       unlaplacian42r(i1,i2,i3,kd)=unxx43r(i1,i2,i3,kd)+unyy43r(i1,i2,
     & i3,kd)
       unlaplacian43r(i1,i2,i3,kd)=unxx43r(i1,i2,i3,kd)+unyy43r(i1,i2,
     & i3,kd)+unzz43r(i1,i2,i3,kd)
       ! efineDifferenceOrder2Components1(v,none)
       ! defineDifferenceOrder4Components1(v,none)
       ! defineDifferenceOrder2Components1(um,none)
       ! defineDifferenceOrder2Components1(ff,none)
       !    *** 2nd order ***
       lap2d2(i1,i2,i3,c)=(u(i1+1,i2,i3,c)-2.*u(i1,i2,i3,c)+u(i1-1,i2,
     & i3,c))*dxsqi+(u(i1,i2+1,i3,c)-2.*u(i1,i2,i3,c)+u(i1,i2-1,i3,c))
     & *dysqi
      !-  lap2d2m(i1,i2,i3,c)=(um(i1+1,i2,i3,c)-2.*um(i1,i2,i3,c)+um(i1-1,i2,i3,c))*dxsqi!-                     +(um(i1,i2+1,i3,c)-2.*um(i1,i2,i3,c)+um(i1,i2-1,i3,c))*dysqi
      !- 
      !-  lap3d2(i1,i2,i3,c)=(u(i1+1,i2,i3,c)-2.*u(i1,i2,i3,c)+u(i1-1,i2,i3,c))*dxsqi!-                    +(u(i1,i2+1,i3,c)-2.*u(i1,i2,i3,c)+u(i1,i2-1,i3,c))*dysqi!-                    +(u(i1,i2,i3+1,c)-2.*u(i1,i2,i3,c)+u(i1,i2,i3-1,c))*dzsqi
      !- 
      !-  lap3d2m(i1,i2,i3,c)=(um(i1+1,i2,i3,c)-2.*um(i1,i2,i3,c)+um(i1-1,i2,i3,c))*dxsqi!-                     +(um(i1,i2+1,i3,c)-2.*um(i1,i2,i3,c)+um(i1,i2-1,i3,c))*dysqi!-                     +(um(i1,i2,i3+1,c)-2.*um(i1,i2,i3,c)+um(i1,i2,i3-1,c))*dzsqi
      !- 
      !-  lap2d2f(i1,i2,i3,c,m)=(fa(i1+1,i2,i3,c,m)-2.*fa(i1,i2,i3,c,m)+fa(i1-1,i2,i3,c,m))*dxsqi!-                       +(fa(i1,i2+1,i3,c,m)-2.*fa(i1,i2,i3,c,m)+fa(i1,i2-1,i3,c,m))*dysqi
      !- 
      !-  lap3d2f(i1,i2,i3,c,m)=(fa(i1+1,i2,i3,c,m)-2.*fa(i1,i2,i3,c,m)+fa(i1-1,i2,i3,c,m))*dxsqi!-                       +(fa(i1,i2+1,i3,c,m)-2.*fa(i1,i2,i3,c,m)+fa(i1,i2-1,i3,c,m))*dysqi!-                       +(fa(i1,i2,i3+1,c,m)-2.*fa(i1,i2,i3,c,m)+fa(i1,i2,i3-1,c,m))*dzsqi
      !- 
      !- 
      ! 2D laplacian squared = u.xxxx + 2 u.xxyy + u.yyyy
      lap2d2Pow2(i1,i2,i3,c)= ( 6.*u(i1,i2,i3,c)   - 4.*(u(i1+1,i2,i3,
     & c)+u(i1-1,i2,i3,c))    +(u(i1+2,i2,i3,c)+u(i1-2,i2,i3,c)) )*
     & dxi4 +( 6.*u(i1,i2,i3,c)    -4.*(u(i1,i2+1,i3,c)+u(i1,i2-1,i3,
     & c))    +(u(i1,i2+2,i3,c)+u(i1,i2-2,i3,c)) )*dyi4  +( 8.*u(i1,
     & i2,i3,c)     -4.*(u(i1+1,i2,i3,c)+u(i1-1,i2,i3,c)+u(i1,i2+1,i3,
     & c)+u(i1,i2-1,i3,c))   +2.*(u(i1+1,i2+1,i3,c)+u(i1-1,i2+1,i3,c)+
     & u(i1+1,i2-1,i3,c)+u(i1-1,i2-1,i3,c)) )*dxdyi2
      !- 
      !-  ! 3D laplacian squared = u.xxxx + u.yyyy + u.zzzz + 2 (u.xxyy + u.xxzz + u.yyzz )
      !-  lap3d2Pow2(i1,i2,i3,c)= ( 6.*u(i1,i2,i3,c)   !-    - 4.*(u(i1+1,i2,i3,c)+u(i1-1,i2,i3,c))    !-        +(u(i1+2,i2,i3,c)+u(i1-2,i2,i3,c)) )*dxi4 !-   +(  +6.*u(i1,i2,i3,c)    !-     -4.*(u(i1,i2+1,i3,c)+u(i1,i2-1,i3,c))    !-        +(u(i1,i2+2,i3,c)+u(i1,i2-2,i3,c)) )*dyi4!-   +(  +6.*u(i1,i2,i3,c)    !-     -4.*(u(i1,i2,i3+1,c)+u(i1,i2,i3-1,c))    !-        +(u(i1,i2,i3+2,c)+u(i1,i2,i3-2,c)) )*dzi4!-    +(8.*u(i1,i2,i3,c)     !-     -4.*(u(i1+1,i2,i3,c)+u(i1-1,i2,i3,c)+u(i1,i2+1,i3,c)+u(i1,i2-1,i3,c))   !-     +2.*(u(i1+1,i2+1,i3,c)+u(i1-1,i2+1,i3,c)+u(i1+1,i2-1,i3,c)+u(i1-1,i2-1,i3,c)) )*dxdyi2 !-    +(8.*u(i1,i2,i3,c)     !-     -4.*(u(i1+1,i2,i3,c)+u(i1-1,i2,i3,c)+u(i1,i2,i3+1,c)+u(i1,i2,i3-1,c))   !-     +2.*(u(i1+1,i2,i3+1,c)+u(i1-1,i2,i3+1,c)+u(i1+1,i2,i3-1,c)+u(i1-1,i2,i3-1,c)) )*dxdzi2 !-    +(8.*u(i1,i2,i3,c)     !-     -4.*(u(i1,i2+1,i3,c)+u(i1,i2-1,i3,c)+u(i1,i2,i3+1,c)+u(i1,i2,i3-1,c))   !-     +2.*(u(i1,i2+1,i3+1,c)+u(i1,i2-1,i3+1,c)+u(i1,i2+1,i3-1,c)+u(i1,i2-1,i3-1,c)) )*dydzi2 
      !- 
      !-  lap2d2Pow3(i1,i2,i3,c)=LAP2D2(lap2d2Pow2,i1,i2,i3,c)
      !- 
      !-  lap3d2Pow3(i1,i2,i3,c)=LAP3D2(lap3d2Pow2,i1,i2,i3,c)
      !- 
      !-  lap2d2Pow4(i1,i2,i3,c)=LAP2D2POW2(lap2d2Pow2,i1,i2,i3,c)
      !-  lap3d2Pow4(i1,i2,i3,c)=LAP3D2POW2(lap3d2Pow2,i1,i2,i3,c)
      !-  
      !    ** 4th order ****
      lap2d4(i1,i2,i3,c)=( -30.*u(i1,i2,i3,c)     +16.*(u(i1+1,i2,i3,c)
     & +u(i1-1,i2,i3,c))     -(u(i1+2,i2,i3,c)+u(i1-2,i2,i3,c)) )*
     & dxsq12i + ( -30.*u(i1,i2,i3,c)     +16.*(u(i1,i2+1,i3,c)+u(i1,
     & i2-1,i3,c))     -(u(i1,i2+2,i3,c)+u(i1,i2-2,i3,c)) )*dysq12i
      !- 
      !-  lap3d4(i1,i2,i3,c)=lap2d4(i1,i2,i3,c)+ !-   ( -30.*u(i1,i2,i3,c)      !-    +16.*(u(i1,i2,i3+1,c)+u(i1,i2,i3-1,c))      !-        -(u(i1,i2,i3+2,c)+u(i1,i2,i3-2,c)) )*dzsq12i 
      !- 
      !-  lap2d4Pow2(i1,i2,i3,c)=LAP2D4(lap2d4,i1,i2,i3,c)
      !-  lap3d4Pow2(i1,i2,i3,c)=LAP3D4(lap3d4,i1,i2,i3,c)
      !- 
      !-  lap2d4Pow3(i1,i2,i3,c)=LAP2D4(lap2d4Pow2,i1,i2,i3,c)
      !-  lap3d4Pow3(i1,i2,i3,c)=LAP3D4(lap3d4Pow2,i1,i2,i3,c)
      !- 
      !- !     *** 6th order ***
      !- 
      !-  lap2d6(i1,i2,i3,c)= !-           c00lap2d6*u(i1,i2,i3,c)     !-          +c10lap2d6*(u(i1+1,i2,i3,c)+u(i1-1,i2,i3,c)) !-          +c01lap2d6*(u(i1,i2+1,i3,c)+u(i1,i2-1,i3,c)) !-          +c20lap2d6*(u(i1+2,i2,i3,c)+u(i1-2,i2,i3,c)) !-          +c02lap2d6*(u(i1,i2+2,i3,c)+u(i1,i2-2,i3,c)) !-          +c30lap2d6*(u(i1+3,i2,i3,c)+u(i1-3,i2,i3,c)) !-          +c03lap2d6*(u(i1,i2+3,i3,c)+u(i1,i2-3,i3,c)) 
      !- 
      !-  lap3d6(i1,i2,i3,c)=!-           c000lap3d6*u(i1,i2,i3,c) !-          +c100lap3d6*(u(i1+1,i2,i3,c)+u(i1-1,i2,i3,c)) !-          +c010lap3d6*(u(i1,i2+1,i3,c)+u(i1,i2-1,i3,c)) !-          +c001lap3d6*(u(i1,i2,i3+1,c)+u(i1,i2,i3-1,c)) !-          +c200lap3d6*(u(i1+2,i2,i3,c)+u(i1-2,i2,i3,c)) !-          +c020lap3d6*(u(i1,i2+2,i3,c)+u(i1,i2-2,i3,c)) !-          +c002lap3d6*(u(i1,i2,i3+2,c)+u(i1,i2,i3-2,c)) !-          +c300lap3d6*(u(i1+3,i2,i3,c)+u(i1-3,i2,i3,c)) !-          +c030lap3d6*(u(i1,i2+3,i3,c)+u(i1,i2-3,i3,c)) !-          +c003lap3d6*(u(i1,i2,i3+3,c)+u(i1,i2,i3-3,c))
      !- 
      !-  lap2d6Pow2(i1,i2,i3,c)=LAP2D6(lap2d6,i1,i2,i3,c)
      !-  lap3d6Pow2(i1,i2,i3,c)=LAP3D6(lap3d6,i1,i2,i3,c)
      !- 
      !- 
      !- !     *** 8th order ***
      !- 
      !-  lap2d8(i1,i2,i3,c)=c00lap2d8*u(i1,i2,i3,c)      !-           +c10lap2d8*(u(i1+1,i2,i3,c)+u(i1-1,i2,i3,c))     !-           +c01lap2d8*(u(i1,i2+1,i3,c)+u(i1,i2-1,i3,c)) !-           +c20lap2d8*(u(i1+2,i2,i3,c)+u(i1-2,i2,i3,c))  !-           +c02lap2d8*(u(i1,i2+2,i3,c)+u(i1,i2-2,i3,c)) !-           +c30lap2d8*(u(i1+3,i2,i3,c)+u(i1-3,i2,i3,c))  !-           +c03lap2d8*(u(i1,i2+3,i3,c)+u(i1,i2-3,i3,c)) !-           +c40lap2d8*(u(i1+4,i2,i3,c)+u(i1-4,i2,i3,c))  !-           +c04lap2d8*(u(i1,i2+4,i3,c)+u(i1,i2-4,i3,c))
      !- 
      !-  lap3d8(i1,i2,i3,c)=c000lap3d8*u(i1,i2,i3,c)      !-           +c100lap3d8*(u(i1+1,i2,i3,c)+u(i1-1,i2,i3,c))     !-           +c010lap3d8*(u(i1,i2+1,i3,c)+u(i1,i2-1,i3,c)) !-           +c001lap3d8*(u(i1,i2,i3+1,c)+u(i1,i2,i3-1,c)) !-           +c200lap3d8*(u(i1+2,i2,i3,c)+u(i1-2,i2,i3,c))  !-           +c020lap3d8*(u(i1,i2+2,i3,c)+u(i1,i2-2,i3,c)) !-           +c002lap3d8*(u(i1,i2,i3+2,c)+u(i1,i2,i3-2,c)) !-           +c300lap3d8*(u(i1+3,i2,i3,c)+u(i1-3,i2,i3,c))  !-           +c030lap3d8*(u(i1,i2+3,i3,c)+u(i1,i2-3,i3,c)) !-           +c003lap3d8*(u(i1,i2,i3+3,c)+u(i1,i2,i3-3,c)) !-           +c400lap3d8*(u(i1+4,i2,i3,c)+u(i1-4,i2,i3,c))  !-           +c040lap3d8*(u(i1,i2+4,i3,c)+u(i1,i2-4,i3,c)) !-           +c004lap3d8*(u(i1,i2,i3+4,c)+u(i1,i2,i3-4,c))
      !- 
      !- ! ******* artificial dissipation ******
      !-  du(i1,i2,i3,c)=u(i1,i2,i3,c)-um(i1,i2,i3,c)
      !- 
      !- !      (2nd difference)
      !-  fd22d(i1,i2,i3,c)= !-  (     ( du(i1-1,i2,i3,c)+du(i1+1,i2,i3,c)+du(i1,i2-1,i3,c)+du(i1,i2+1,i3,c) ) !-   -4.*du(i1,i2,i3,c) )
      !- !
      !-  fd23d(i1,i2,i3,c)=!-  (     ( du(i1-1,i2,i3,c)+du(i1+1,i2,i3,c)+du(i1,i2-1,i3,c)+du(i1,i2+1,i3,c)+du(i1,i2,i3-1,c)+du(i1,i2,i3+1,c) ) !-    -6.*du(i1,i2,i3,c) )
      !- 
      !- !     -(fourth difference)
      !-  fd42d(i1,i2,i3,c)= !-  (    -( du(i1-2,i2,i3,c)+du(i1+2,i2,i3,c)+du(i1,i2-2,i3,c)+du(i1,i2+2,i3,c) ) !-    +4.*( du(i1-1,i2,i3,c)+du(i1+1,i2,i3,c)+du(i1,i2-1,i3,c)+du(i1,i2+1,i3,c) ) !-   -12.*du(i1,i2,i3,c) )
      !- !
      !-  fd43d(i1,i2,i3,c)=!-  (    -( du(i1-2,i2,i3,c)+du(i1+2,i2,i3,c)+du(i1,i2-2,i3,c)+du(i1,i2+2,i3,c)+du(i1,i2,i3-2,c)+du(i1,i2,i3+2,c) ) !-    +4.*( du(i1-1,i2,i3,c)+du(i1+1,i2,i3,c)+du(i1,i2-1,i3,c)+du(i1,i2+1,i3,c)+du(i1,i2,i3-1,c)+du(i1,i2,i3+1,c) ) !-   -18.*du(i1,i2,i3,c) )
      !- 
      !-  ! (sixth  difference)
      !-  fd62d(i1,i2,i3,c)= !-  (     ( du(i1-3,i2,i3,c)+du(i1+3,i2,i3,c)+du(i1,i2-3,i3,c)+du(i1,i2+3,i3,c) ) !-    -6.*( du(i1-2,i2,i3,c)+du(i1+2,i2,i3,c)+du(i1,i2-2,i3,c)+du(i1,i2+2,i3,c) ) !-   +15.*( du(i1-1,i2,i3,c)+du(i1+1,i2,i3,c)+du(i1,i2-1,i3,c)+du(i1,i2+1,i3,c) ) !-   -40.*du(i1,i2,i3,c) )
      !- 
      !-  fd63d(i1,i2,i3,c)=!-  (     ( du(i1-3,i2,i3,c)+du(i1+3,i2,i3,c)+du(i1,i2-3,i3,c)+du(i1,i2+3,i3,c)+du(i1,i2,i3-3,c)+du(i1,i2,i3+3,c) ) !-    -6.*( du(i1-2,i2,i3,c)+du(i1+2,i2,i3,c)+du(i1,i2-2,i3,c)+du(i1,i2+2,i3,c)+du(i1,i2,i3-2,c)+du(i1,i2,i3+2,c) ) !-   +15.*( du(i1-1,i2,i3,c)+du(i1+1,i2,i3,c)+du(i1,i2-1,i3,c)+du(i1,i2+1,i3,c)+du(i1,i2,i3-1,c)+du(i1,i2,i3+1,c) ) !-   -60.*du(i1,i2,i3,c) )
      !- 
      !-  ! -(eighth  difference)
      !-  fd82d(i1,i2,i3,c)= !-  (    -( du(i1-4,i2,i3,c)+du(i1+4,i2,i3,c)+du(i1,i2-4,i3,c)+du(i1,i2+4,i3,c) ) !-    +8.*( du(i1-3,i2,i3,c)+du(i1+3,i2,i3,c)+du(i1,i2-3,i3,c)+du(i1,i2+3,i3,c) ) !-   -28.*( du(i1-2,i2,i3,c)+du(i1+2,i2,i3,c)+du(i1,i2-2,i3,c)+du(i1,i2+2,i3,c) ) !-   +56.*( du(i1-1,i2,i3,c)+du(i1+1,i2,i3,c)+du(i1,i2-1,i3,c)+du(i1,i2+1,i3,c) ) !-  -140.*du(i1,i2,i3,c) )
      !- 
      !-  fd83d(i1,i2,i3,c)=!-  (    -( du(i1-4,i2,i3,c)+du(i1+4,i2,i3,c)+du(i1,i2-4,i3,c)+du(i1,i2+4,i3,c)+du(i1,i2,i3-4,c)+du(i1,i2,i3+4,c) ) !-    +8.*( du(i1-3,i2,i3,c)+du(i1+3,i2,i3,c)+du(i1,i2-3,i3,c)+du(i1,i2+3,i3,c)+du(i1,i2,i3-3,c)+du(i1,i2,i3+3,c) ) !-   -28.*( du(i1-2,i2,i3,c)+du(i1+2,i2,i3,c)+du(i1,i2-2,i3,c)+du(i1,i2+2,i3,c)+du(i1,i2,i3-2,c)+du(i1,i2,i3+2,c) ) !-   +56.*( du(i1-1,i2,i3,c)+du(i1+1,i2,i3,c)+du(i1,i2-1,i3,c)+du(i1,i2+1,i3,c)+du(i1,i2,i3-1,c)+du(i1,i2,i3+1,c) ) !-  -210.*du(i1,i2,i3,c) )
      !...........end   statement functions
       useOpt=.true. ! is true, use new optimized version
       debug=0
       cc    =rpar(0)  ! this is c
       dt    =rpar(1)
       dx(0) =rpar(2)
       dx(1) =rpar(3)
       dx(2) =rpar(4)
       adc   =rpar(5)  ! coefficient of artificial dissipation
       add   =rpar(6)  ! coefficient of divergence damping
       dr(0) =rpar(7)
       dr(1) =rpar(8)
       dr(2) =rpar(9)
       eps   =rpar(10)
       mu    =rpar(11)
       kx    =rpar(12)
       ky    =rpar(13)
       kz    =rpar(14)
       sigmaE=rpar(15)  ! electric conductivity (for lossy materials, complex index of refraction)
       sigmaH=rpar(16)  ! magnetic conductivity
       divergenceCleaningCoefficient=rpar(17)
       t     =rpar(18)
       ep    =rpar(19)  ! for TZ forcing
       rpar(20)=0.  ! return the time used for adding dissipation
       dy=dx(1)  ! Are these needed?
       dz=dx(2)
       ! timeForArtificialDissipation=rpar(6) ! return value
       option             =ipar(0)
       gridType           =ipar(1)
       orderOfAccuracy    =ipar(2)
       orderInTime        =ipar(3)
       addForcing         =ipar(4)
       orderOfDissipation =ipar(5)
       ex                 =ipar(6)
       ey                 =ipar(7)
       ez                 =ipar(8)
       hx                 =ipar(9)
       hy                 =ipar(10)
       hz                 =ipar(11)
       solveForE          =ipar(12)
       solveForH          =ipar(13)
       useWhereMask       =ipar(14)
       timeSteppingMethod =ipar(15)
       useVariableDissipation        =ipar(16)
       useCurvilinearOpt             =ipar(17)
       useConservative               =ipar(18)
       combineDissipationWithAdvance =ipar(19)
       useDivergenceCleaning         =ipar(20)
       useNewForcingMethod           =ipar(21)
       numberOfForcingFunctions      =ipar(22)
       fcur                          =ipar(23)
       dispersionModel     =ipar(24)
       pxc                 =ipar(25)
       pyc                 =ipar(26)
       pzc                 =ipar(27)
       ! qxc                 =ipar(28)
       grid                =ipar(29)
       ! qzc                 =ipar(30)
       ! rxc                 =ipar(31)
       ! ryc                 =ipar(32)
       ! rzc                 =ipar(33)
       ! useSosupDissipation,
       ! sosupDissipationOption,
       updateInterior        =ipar(36).eq.1
       addDissipation        =ipar(37).eq.1
       forcingOption         =ipar(39)
       ! addDissipation=.true. if we add the dissipation in the dis(i1,i2,i3,c) array
       !  if combineDissipationWithAdvance.ne.0 we compute the dissipation on the fly in the time step
       !  rather than pre-computing it in diss(i1,i2,i3,c)
       ! OLD addDissipation = adc.gt.0. .and. combineDissipationWithAdvance.eq.0
       ! adcdt=adc*dt
       solveForAllFields       = ipar(40)
       materialType            = ipar(41)
       numberOfMaterialRegions = ipar(42)
       useSuperGrid            = ipar(43)
       useAbsorbingLayer(0)    = ipar(44)
       useAbsorbingLayer(1)    = ipar(45)
       useAbsorbingLayer(2)    = ipar(46)
       if( nd.eq.2 .and. solveForAllFields.eq.0 )then
         numComp=3  ! TEZ has 3 components
       else
         numComp=6
       end if
       if( t.le.2*dt )then
         write(*,*) 'Inside advBA3dOrder4r...'
         write(*,'("addForcing=",i2, " solveForAllFields=",i2," 
     & useSuperGrid=",i2)') addForcing,solveForAllFields
         write(*,'(" useSuperGrid=",i2," useAbsorbingLayer(0:2)=",3i2)
     & ') useSuperGrid,(useAbsorbingLayer(n),n=0,2)
         write(*,'("addDissipation=",l2, " updateInterior=",l2)') 
     & addDissipation,updateInterior
         write(*,'("dispersionModel=",i2)') dispersionModel
         if( timeSteppingMethod.eq.modifiedEquationTimeStepping )then
            write(*,*) 'timeSteppingMethod =  
     & modifiedEquationTimeStepping'
         else if( timeSteppingMethod.eq.rungeKutta ) then
            write(*,*) 'timeSteppingMethod = rungeKutta'
         end if
         if( .false. .and. useSuperGrid.ne.0 )then
           write(*,'(" etax=",100(f6.3,1x))') (etax(i1),i1=nd1a,nd1b)
        end if
       end if
       if( numberOfMaterialRegions.gt.maxRegions )then
          stop 1002
       end if
       ! 3x3 Material matrix for TEz polarization
       ! We use ex=0,ey=1 and hz=5 entries in K0i(0:5,0:5) 
       do mr=0,numberOfMaterialRegions-1
         Ki(0,0,mr) = K0i(0,0,mr)
         Ki(0,1,mr) = K0i(0,1,mr)
         Ki(0,2,mr) = K0i(0,5,mr)
         Ki(1,0,mr) = K0i(1,0,mr)
         Ki(1,1,mr) = K0i(1,1,mr)
         Ki(1,2,mr) = K0i(1,5,mr)
         Ki(2,0,mr) = K0i(5,0,mr)
         Ki(2,1,mr) = K0i(5,1,mr)
         Ki(2,2,mr) = K0i(5,5,mr)
       end do
       if( t.lt.dt )then
         write(*,*) 'materialType=',materialType
         write(*,*) 'numberOfMaterialRegions=',numberOfMaterialRegions
         do mr=0,numberOfMaterialRegions-1
           write(*,*) 'Material region=',mr
           write(*,'("K0i=",6("[",6(f6.3,1x),"]",/,4x))') ((K0i(i1,i2,
     & mr),i1=0,5),i2=0,5)
           if( solveForAllFields .eq. 0 )then
             write(*,'("Ki =",3("[",3(f6.3,1x),"]",/,4x))') ((Ki(i1,i2,
     & mr),i1=0,2),i2=0,2)
           end if
         end do
         if( numberOfMaterialRegions>1 .and. n1b-n1a.lt.20 )then
           do i3=n3a,n3b
           do i2=n2a,n2b
             write(*,*) 'matMask=',(matMask(i1,i2,i3),i1=n1a,n1b)
           end do
           end do
         end if
       end if
       if( dispersionModel.ne.noDispersion )then
         ! get the BA gdm parameters
         !  gdmPar(4,NpMax,6,6,0:maxRegions-1)
         call getBAGDMParameters( grid, gdmPar,Np, NpMax,maxRegions )
         ! count the total number of polarization terms: 
         numPolarizationTerms=0
         do mr=0,numberOfMaterialRegions-1
           do k1=1,6
             do k2=1,6
               numPolarizationTerms = numPolarizationTerms + Np(k1,k2,
     & mr)
             end do
           end do
         end do
         numPolarizationTerms = numPolarizationTerms*2  ! we store p and p.t
         if( numPolarizationTerms > maxNumPolarizationTerms )then
            write(*,'("advBA: ERROR: numPolarizationTerms > 
     & maxNumPolarizationTerms")');
            stop 1234
         end if
         if( t.eq.0. .and. dispersionModel.ne.noDispersion )then
           ! ---- Dispersive Maxwell ----
           write(*,'("--advOpt-- dispersionModel=",i4," 
     & numPolarizationTerms=",i6)') dispersionModel,
     & numPolarizationTerms
           !write(*,'("--advOpt-- GDM: numberOfPolarizationVectors=",i4," alphaP=",e8.2)') numberOfPolarizationVectors,alphaP
           !write(*,'("--advOpt-- GDM: alphaP,a0,a1,b0,b1=",5(1p,e10.2))') alphaP,a0,a1,b0,b1
           !do iv=0,numberOfPolarizationVectors-1
           !  write(*,'("--advOpt-- GDM: eqn=",i3," a0,a1,b0,b1=",4(1p,e10.2))') iv,a0v(iv),a1v(iv),b0v(iv),b1v(iv)
           !end do
           if( .false. )then
            do mr=0,numberOfMaterialRegions-1
              write(*,'("BA-GDM: material region mr=",i2)') mr
              ! write(*,'((5x,6i3))') ((Np(k1,k2,mr),k1=1,6),k2=1,6)
              do k1=1,6
              do k2=1,6
                 if( Np(k1,k2,mr) > 0 )then
                   write(*,'("  K(",i1,",",i1,") : Np=",i3," :")') k1,
     & k2,Np(k1,k2,mr)
                   do n=1,Np(k1,k2,mr)
                     write(*,'("    n=",i3," [a0,a1,b0,b1]=",4(e9.3,1x)
     & )') n,(gdmPar(m,n,k1,k2,mr),m=1,4)
                   end do
                 end if
              end do
              end do
            end do
           end if
          ! stop 3333
        end if
       end if
       methodOfLines=.false.
       if( timeSteppingMethod.ne.modifiedEquationTimeStepping )then
          methodOfLines=.true.
       end if
       fprev = mod(fcur-1+numberOfForcingFunctions,max(1,
     & numberOfForcingFunctions))
       fnext = mod(fcur+1                         ,max(1,
     & numberOfForcingFunctions))
       csq=cc**2
       dtsq=dt**2
       cdt=cc*dt
       cdtsq=(cc**2)*(dt**2)
       cdtsq12=cdtsq*cdtsq/12.  ! c^4 dt^4 /14
       cdt4by360=(cdt)**4/360.
       cdt6by20160=cdt**6/(8.*7.*6.*5.*4.*3.)
       csqdt = cc**2*dt
       cdtSqBy12= cdtsq/12.   ! c^2*dt*2/12
       dt4by12=dtsq*dtsq/12.
       cdtdx = (cc*dt/dx(0))**2
       cdtdy = (cc*dt/dy)**2
       cdtdz = (cc*dt/dz)**2
       dxsqi=1./(dx(0)**2)
       dysqi=1./(dy**2)
       dzsqi=1./(dz**2)
       dxsq12i=1./(12.*dx(0)**2)
       dysq12i=1./(12.*dy**2)
       dzsq12i=1./(12.*dz**2)
       dxi4=1./(dx(0)**4)
       dyi4=1./(dy**4)
       dxdyi2=1./(dx(0)*dx(0)*dy*dy)
       dzi4=1./(dz**4)
       dxdzi2=1./(dx(0)*dx(0)*dz*dz)
       dydzi2=1./(dy*dy*dz*dz)
       do m=0,10
          fv(m)=0.  ! temp for forcing
       end do
       if( t.eq.0. .and. dispersionModel.ne.noDispersion )then
          write(*,'("--advOpt-- dispersionModel=",i4," px,py,pz=",3i2)
     & ') dispersionModel,pxc,pyc,pzc
       end if
      ! write(*,'(" advMaxwell: timeSteppingMethod=",i2)') timeSteppingMethod
       if( timeSteppingMethod.eq.defaultTimeStepping )then
        write(*,'(" advMaxwell:ERROR: 
     & timeSteppingMethod=defaultTimeStepping -- this should be set")
     & ')
          ! '
        stop 83322
       end if
       if( .not. addDissipation )then
         ! **new way ** Dec 17, 2019 
         ! All methods are done here 
             if( dispersionModel.eq.noDispersion )then
               if( useSuperGrid.eq.0  ) then
                   if( t.lt.2*dt )then
                     write(*,'("advBA: advance BA dim=3, order=4 
     & grid=rectangular polar=NONE... t=",e10.2)') t
                   end if
                   mr=0
                   if( .not.methodOfLines )then
                     ! --- TAYLOR TIME-STEPPING --- 
                     stop 111
                   else
                     ! --- METHOD OF LINES (RK) ---
                       do i3=n3a,n3b
                       do i2=n2a,n2b
                       do i1=n1a,n1b
                         if( mask(i1,i2,i3).gt.0 )then
                       if( addForcing.ne.0 )then  ! do this for now *fix me*
                           do m=0,5
                            fv(m)=f(i1,i2,i3,m)
                           end do
                       end if
                       if( numberOfMaterialRegions.gt.1 )then
                         mr = matMask(i1,i2,i3)
                         if( mr.lt.0 .or. 
     & mr.ge.numberOfMaterialRegions )then  ! do this for now
                            stop 9999
                         end if
                       end if
                        ! --- 3D -----
                         curl(0) =   uy43r(i1,i2,i3,hz)- uz43r(i1,i2,
     & i3,hy)
                         curl(1) =   uz43r(i1,i2,i3,hx)- ux43r(i1,i2,
     & i3,hz)
                         curl(2) =   ux43r(i1,i2,i3,hy)- uy43r(i1,i2,
     & i3,hx)
                         curl(3) =-( uy43r(i1,i2,i3,ez)- uz43r(i1,i2,
     & i3,ey))
                         curl(4) =-( uz43r(i1,i2,i3,ex)- ux43r(i1,i2,
     & i3,ez))
                         curl(5) =-( ux43r(i1,i2,i3,ey)- uy43r(i1,i2,
     & i3,ex))
                         do m=0,5
                           un(i1,i2,i3,m)= K0i(m,0,mr)*(curl(0)+fv(0)) 
     & + K0i(m,1,mr)*(curl(1)+fv(1)) + K0i(m,2,mr)*(curl(2)+fv(2)) + 
     & K0i(m,3,mr)*(curl(3)+fv(3)) + K0i(m,4,mr)*(curl(4)+fv(4)) + 
     & K0i(m,5,mr)*(curl(5)+fv(5))
                         end do
                         end if
                       end do
                       end do
                       end do
                     if( .false. .or. debug > 15 )then
                       stop  4747
                    end if
                   end if ! end MOL
               else
                 ! --- SUPERGRID ---
                 if( t.le.3*dt )then
                   write(*,'(" USE SUPERGRID...")' )
                 end if
                     !  --- THREE-DIMENSIONS ---
                     if( useAbsorbingLayer(0).eq.1 .and. 
     & useAbsorbingLayer(1).eq.1 .and. useAbsorbingLayer(2).eq.1 )then
                         if( t.lt.2*dt )then
                           write(*,'("advBA: advance BA dim=3, order=4 
     & grid=rectangular polar=NONE... t=",e10.2)') t
                         end if
                         mr=0
                         if( .not.methodOfLines )then
                           ! --- TAYLOR TIME-STEPPING --- 
                           stop 111
                         else
                           ! --- METHOD OF LINES (RK) ---
                             do i3=n3a,n3b
                             do i2=n2a,n2b
                             do i1=n1a,n1b
                               if( mask(i1,i2,i3).gt.0 )then
                             if( addForcing.ne.0 )then  ! do this for now *fix me*
                                 do m=0,5
                                  fv(m)=f(i1,i2,i3,m)
                                 end do
                             end if
                             if( numberOfMaterialRegions.gt.1 )then
                               mr = matMask(i1,i2,i3)
                               if( mr.lt.0 .or. 
     & mr.ge.numberOfMaterialRegions )then  ! do this for now
                                  stop 9999
                               end if
                             end if
                              ! --- 3D -----
                               curl(0) =  etay(i2)* uy43r(i1,i2,i3,hz)-
     & etaz(i3)* uz43r(i1,i2,i3,hy)
                               curl(1) =  etaz(i3)* uz43r(i1,i2,i3,hx)-
     & etax(i1)* ux43r(i1,i2,i3,hz)
                               curl(2) =  etax(i1)* ux43r(i1,i2,i3,hy)-
     & etay(i2)* uy43r(i1,i2,i3,hx)
                               curl(3) =-(etay(i2)* uy43r(i1,i2,i3,ez)-
     & etaz(i3)* uz43r(i1,i2,i3,ey))
                               curl(4) =-(etaz(i3)* uz43r(i1,i2,i3,ex)-
     & etax(i1)* ux43r(i1,i2,i3,ez))
                               curl(5) =-(etax(i1)* ux43r(i1,i2,i3,ey)-
     & etay(i2)* uy43r(i1,i2,i3,ex))
                               do m=0,5
                                 un(i1,i2,i3,m)= K0i(m,0,mr)*(curl(0)+
     & fv(0)) + K0i(m,1,mr)*(curl(1)+fv(1)) + K0i(m,2,mr)*(curl(2)+fv(
     & 2)) + K0i(m,3,mr)*(curl(3)+fv(3)) + K0i(m,4,mr)*(curl(4)+fv(4))
     &  + K0i(m,5,mr)*(curl(5)+fv(5))
                               end do
                               end if
                             end do
                             end do
                             end do
                           if( .false. .or. debug > 15 )then
                             stop  4747
                          end if
                         end if ! end MOL
                     else if( useAbsorbingLayer(0).eq.1 .and. 
     & useAbsorbingLayer(1).eq.1 .and. useAbsorbingLayer(2).eq.0 )then
                         if( t.lt.2*dt )then
                           write(*,'("advBA: advance BA dim=3, order=4 
     & grid=rectangular polar=NONE... t=",e10.2)') t
                         end if
                         mr=0
                         if( .not.methodOfLines )then
                           ! --- TAYLOR TIME-STEPPING --- 
                           stop 111
                         else
                           ! --- METHOD OF LINES (RK) ---
                             do i3=n3a,n3b
                             do i2=n2a,n2b
                             do i1=n1a,n1b
                               if( mask(i1,i2,i3).gt.0 )then
                             if( addForcing.ne.0 )then  ! do this for now *fix me*
                                 do m=0,5
                                  fv(m)=f(i1,i2,i3,m)
                                 end do
                             end if
                             if( numberOfMaterialRegions.gt.1 )then
                               mr = matMask(i1,i2,i3)
                               if( mr.lt.0 .or. 
     & mr.ge.numberOfMaterialRegions )then  ! do this for now
                                  stop 9999
                               end if
                             end if
                              ! --- 3D -----
                               curl(0) =  etay(i2)* uy43r(i1,i2,i3,hz)-
     &  uz43r(i1,i2,i3,hy)
                               curl(1) =   uz43r(i1,i2,i3,hx)-etax(i1)*
     &  ux43r(i1,i2,i3,hz)
                               curl(2) =  etax(i1)* ux43r(i1,i2,i3,hy)-
     & etay(i2)* uy43r(i1,i2,i3,hx)
                               curl(3) =-(etay(i2)* uy43r(i1,i2,i3,ez)-
     &  uz43r(i1,i2,i3,ey))
                               curl(4) =-( uz43r(i1,i2,i3,ex)-etax(i1)*
     &  ux43r(i1,i2,i3,ez))
                               curl(5) =-(etax(i1)* ux43r(i1,i2,i3,ey)-
     & etay(i2)* uy43r(i1,i2,i3,ex))
                               do m=0,5
                                 un(i1,i2,i3,m)= K0i(m,0,mr)*(curl(0)+
     & fv(0)) + K0i(m,1,mr)*(curl(1)+fv(1)) + K0i(m,2,mr)*(curl(2)+fv(
     & 2)) + K0i(m,3,mr)*(curl(3)+fv(3)) + K0i(m,4,mr)*(curl(4)+fv(4))
     &  + K0i(m,5,mr)*(curl(5)+fv(5))
                               end do
                               end if
                             end do
                             end do
                             end do
                           if( .false. .or. debug > 15 )then
                             stop  4747
                          end if
                         end if ! end MOL
                     else if( useAbsorbingLayer(0).eq.1 .and. 
     & useAbsorbingLayer(1).eq.0 .and. useAbsorbingLayer(2).eq.1 )then
                         if( t.lt.2*dt )then
                           write(*,'("advBA: advance BA dim=3, order=4 
     & grid=rectangular polar=NONE... t=",e10.2)') t
                         end if
                         mr=0
                         if( .not.methodOfLines )then
                           ! --- TAYLOR TIME-STEPPING --- 
                           stop 111
                         else
                           ! --- METHOD OF LINES (RK) ---
                             do i3=n3a,n3b
                             do i2=n2a,n2b
                             do i1=n1a,n1b
                               if( mask(i1,i2,i3).gt.0 )then
                             if( addForcing.ne.0 )then  ! do this for now *fix me*
                                 do m=0,5
                                  fv(m)=f(i1,i2,i3,m)
                                 end do
                             end if
                             if( numberOfMaterialRegions.gt.1 )then
                               mr = matMask(i1,i2,i3)
                               if( mr.lt.0 .or. 
     & mr.ge.numberOfMaterialRegions )then  ! do this for now
                                  stop 9999
                               end if
                             end if
                              ! --- 3D -----
                               curl(0) =   uy43r(i1,i2,i3,hz)-etaz(i3)*
     &  uz43r(i1,i2,i3,hy)
                               curl(1) =  etaz(i3)* uz43r(i1,i2,i3,hx)-
     & etax(i1)* ux43r(i1,i2,i3,hz)
                               curl(2) =  etax(i1)* ux43r(i1,i2,i3,hy)-
     &  uy43r(i1,i2,i3,hx)
                               curl(3) =-( uy43r(i1,i2,i3,ez)-etaz(i3)*
     &  uz43r(i1,i2,i3,ey))
                               curl(4) =-(etaz(i3)* uz43r(i1,i2,i3,ex)-
     & etax(i1)* ux43r(i1,i2,i3,ez))
                               curl(5) =-(etax(i1)* ux43r(i1,i2,i3,ey)-
     &  uy43r(i1,i2,i3,ex))
                               do m=0,5
                                 un(i1,i2,i3,m)= K0i(m,0,mr)*(curl(0)+
     & fv(0)) + K0i(m,1,mr)*(curl(1)+fv(1)) + K0i(m,2,mr)*(curl(2)+fv(
     & 2)) + K0i(m,3,mr)*(curl(3)+fv(3)) + K0i(m,4,mr)*(curl(4)+fv(4))
     &  + K0i(m,5,mr)*(curl(5)+fv(5))
                               end do
                               end if
                             end do
                             end do
                             end do
                           if( .false. .or. debug > 15 )then
                             stop  4747
                          end if
                         end if ! end MOL
                     else if( useAbsorbingLayer(0).eq.0 .and. 
     & useAbsorbingLayer(1).eq.1 .and. useAbsorbingLayer(2).eq.1 )then
                         if( t.lt.2*dt )then
                           write(*,'("advBA: advance BA dim=3, order=4 
     & grid=rectangular polar=NONE... t=",e10.2)') t
                         end if
                         mr=0
                         if( .not.methodOfLines )then
                           ! --- TAYLOR TIME-STEPPING --- 
                           stop 111
                         else
                           ! --- METHOD OF LINES (RK) ---
                             do i3=n3a,n3b
                             do i2=n2a,n2b
                             do i1=n1a,n1b
                               if( mask(i1,i2,i3).gt.0 )then
                             if( addForcing.ne.0 )then  ! do this for now *fix me*
                                 do m=0,5
                                  fv(m)=f(i1,i2,i3,m)
                                 end do
                             end if
                             if( numberOfMaterialRegions.gt.1 )then
                               mr = matMask(i1,i2,i3)
                               if( mr.lt.0 .or. 
     & mr.ge.numberOfMaterialRegions )then  ! do this for now
                                  stop 9999
                               end if
                             end if
                              ! --- 3D -----
                               curl(0) =  etay(i2)* uy43r(i1,i2,i3,hz)-
     & etaz(i3)* uz43r(i1,i2,i3,hy)
                               curl(1) =  etaz(i3)* uz43r(i1,i2,i3,hx)-
     &  ux43r(i1,i2,i3,hz)
                               curl(2) =   ux43r(i1,i2,i3,hy)-etay(i2)*
     &  uy43r(i1,i2,i3,hx)
                               curl(3) =-(etay(i2)* uy43r(i1,i2,i3,ez)-
     & etaz(i3)* uz43r(i1,i2,i3,ey))
                               curl(4) =-(etaz(i3)* uz43r(i1,i2,i3,ex)-
     &  ux43r(i1,i2,i3,ez))
                               curl(5) =-( ux43r(i1,i2,i3,ey)-etay(i2)*
     &  uy43r(i1,i2,i3,ex))
                               do m=0,5
                                 un(i1,i2,i3,m)= K0i(m,0,mr)*(curl(0)+
     & fv(0)) + K0i(m,1,mr)*(curl(1)+fv(1)) + K0i(m,2,mr)*(curl(2)+fv(
     & 2)) + K0i(m,3,mr)*(curl(3)+fv(3)) + K0i(m,4,mr)*(curl(4)+fv(4))
     &  + K0i(m,5,mr)*(curl(5)+fv(5))
                               end do
                               end if
                             end do
                             end do
                             end do
                           if( .false. .or. debug > 15 )then
                             stop  4747
                          end if
                         end if ! end MOL
                     else if( useAbsorbingLayer(0).eq.1 .and. 
     & useAbsorbingLayer(1).eq.0 .and. useAbsorbingLayer(2).eq.0 )then
                         if( t.lt.2*dt )then
                           write(*,'("advBA: advance BA dim=3, order=4 
     & grid=rectangular polar=NONE... t=",e10.2)') t
                         end if
                         mr=0
                         if( .not.methodOfLines )then
                           ! --- TAYLOR TIME-STEPPING --- 
                           stop 111
                         else
                           ! --- METHOD OF LINES (RK) ---
                             do i3=n3a,n3b
                             do i2=n2a,n2b
                             do i1=n1a,n1b
                               if( mask(i1,i2,i3).gt.0 )then
                             if( addForcing.ne.0 )then  ! do this for now *fix me*
                                 do m=0,5
                                  fv(m)=f(i1,i2,i3,m)
                                 end do
                             end if
                             if( numberOfMaterialRegions.gt.1 )then
                               mr = matMask(i1,i2,i3)
                               if( mr.lt.0 .or. 
     & mr.ge.numberOfMaterialRegions )then  ! do this for now
                                  stop 9999
                               end if
                             end if
                              ! --- 3D -----
                               curl(0) =   uy43r(i1,i2,i3,hz)- uz43r(
     & i1,i2,i3,hy)
                               curl(1) =   uz43r(i1,i2,i3,hx)-etax(i1)*
     &  ux43r(i1,i2,i3,hz)
                               curl(2) =  etax(i1)* ux43r(i1,i2,i3,hy)-
     &  uy43r(i1,i2,i3,hx)
                               curl(3) =-( uy43r(i1,i2,i3,ez)- uz43r(
     & i1,i2,i3,ey))
                               curl(4) =-( uz43r(i1,i2,i3,ex)-etax(i1)*
     &  ux43r(i1,i2,i3,ez))
                               curl(5) =-(etax(i1)* ux43r(i1,i2,i3,ey)-
     &  uy43r(i1,i2,i3,ex))
                               do m=0,5
                                 un(i1,i2,i3,m)= K0i(m,0,mr)*(curl(0)+
     & fv(0)) + K0i(m,1,mr)*(curl(1)+fv(1)) + K0i(m,2,mr)*(curl(2)+fv(
     & 2)) + K0i(m,3,mr)*(curl(3)+fv(3)) + K0i(m,4,mr)*(curl(4)+fv(4))
     &  + K0i(m,5,mr)*(curl(5)+fv(5))
                               end do
                               end if
                             end do
                             end do
                             end do
                           if( .false. .or. debug > 15 )then
                             stop  4747
                          end if
                         end if ! end MOL
                     else if( useAbsorbingLayer(0).eq.0 .and. 
     & useAbsorbingLayer(1).eq.1 .and. useAbsorbingLayer(2).eq.0 )then
                         if( t.lt.2*dt )then
                           write(*,'("advBA: advance BA dim=3, order=4 
     & grid=rectangular polar=NONE... t=",e10.2)') t
                         end if
                         mr=0
                         if( .not.methodOfLines )then
                           ! --- TAYLOR TIME-STEPPING --- 
                           stop 111
                         else
                           ! --- METHOD OF LINES (RK) ---
                             do i3=n3a,n3b
                             do i2=n2a,n2b
                             do i1=n1a,n1b
                               if( mask(i1,i2,i3).gt.0 )then
                             if( addForcing.ne.0 )then  ! do this for now *fix me*
                                 do m=0,5
                                  fv(m)=f(i1,i2,i3,m)
                                 end do
                             end if
                             if( numberOfMaterialRegions.gt.1 )then
                               mr = matMask(i1,i2,i3)
                               if( mr.lt.0 .or. 
     & mr.ge.numberOfMaterialRegions )then  ! do this for now
                                  stop 9999
                               end if
                             end if
                              ! --- 3D -----
                               curl(0) =  etay(i2)* uy43r(i1,i2,i3,hz)-
     &  uz43r(i1,i2,i3,hy)
                               curl(1) =   uz43r(i1,i2,i3,hx)- ux43r(
     & i1,i2,i3,hz)
                               curl(2) =   ux43r(i1,i2,i3,hy)-etay(i2)*
     &  uy43r(i1,i2,i3,hx)
                               curl(3) =-(etay(i2)* uy43r(i1,i2,i3,ez)-
     &  uz43r(i1,i2,i3,ey))
                               curl(4) =-( uz43r(i1,i2,i3,ex)- ux43r(
     & i1,i2,i3,ez))
                               curl(5) =-( ux43r(i1,i2,i3,ey)-etay(i2)*
     &  uy43r(i1,i2,i3,ex))
                               do m=0,5
                                 un(i1,i2,i3,m)= K0i(m,0,mr)*(curl(0)+
     & fv(0)) + K0i(m,1,mr)*(curl(1)+fv(1)) + K0i(m,2,mr)*(curl(2)+fv(
     & 2)) + K0i(m,3,mr)*(curl(3)+fv(3)) + K0i(m,4,mr)*(curl(4)+fv(4))
     &  + K0i(m,5,mr)*(curl(5)+fv(5))
                               end do
                               end if
                             end do
                             end do
                             end do
                           if( .false. .or. debug > 15 )then
                             stop  4747
                          end if
                         end if ! end MOL
                     else if( useAbsorbingLayer(0).eq.0 .and. 
     & useAbsorbingLayer(1).eq.0 .and. useAbsorbingLayer(2).eq.1 )then
                         if( t.lt.2*dt )then
                           write(*,'("advBA: advance BA dim=3, order=4 
     & grid=rectangular polar=NONE... t=",e10.2)') t
                         end if
                         mr=0
                         if( .not.methodOfLines )then
                           ! --- TAYLOR TIME-STEPPING --- 
                           stop 111
                         else
                           ! --- METHOD OF LINES (RK) ---
                             do i3=n3a,n3b
                             do i2=n2a,n2b
                             do i1=n1a,n1b
                               if( mask(i1,i2,i3).gt.0 )then
                             if( addForcing.ne.0 )then  ! do this for now *fix me*
                                 do m=0,5
                                  fv(m)=f(i1,i2,i3,m)
                                 end do
                             end if
                             if( numberOfMaterialRegions.gt.1 )then
                               mr = matMask(i1,i2,i3)
                               if( mr.lt.0 .or. 
     & mr.ge.numberOfMaterialRegions )then  ! do this for now
                                  stop 9999
                               end if
                             end if
                              ! --- 3D -----
                               curl(0) =   uy43r(i1,i2,i3,hz)-etaz(i3)*
     &  uz43r(i1,i2,i3,hy)
                               curl(1) =  etaz(i3)* uz43r(i1,i2,i3,hx)-
     &  ux43r(i1,i2,i3,hz)
                               curl(2) =   ux43r(i1,i2,i3,hy)- uy43r(
     & i1,i2,i3,hx)
                               curl(3) =-( uy43r(i1,i2,i3,ez)-etaz(i3)*
     &  uz43r(i1,i2,i3,ey))
                               curl(4) =-(etaz(i3)* uz43r(i1,i2,i3,ex)-
     &  ux43r(i1,i2,i3,ez))
                               curl(5) =-( ux43r(i1,i2,i3,ey)- uy43r(
     & i1,i2,i3,ex))
                               do m=0,5
                                 un(i1,i2,i3,m)= K0i(m,0,mr)*(curl(0)+
     & fv(0)) + K0i(m,1,mr)*(curl(1)+fv(1)) + K0i(m,2,mr)*(curl(2)+fv(
     & 2)) + K0i(m,3,mr)*(curl(3)+fv(3)) + K0i(m,4,mr)*(curl(4)+fv(4))
     &  + K0i(m,5,mr)*(curl(5)+fv(5))
                               end do
                               end if
                             end do
                             end do
                             end do
                           if( .false. .or. debug > 15 )then
                             stop  4747
                          end if
                         end if ! end MOL
                     else
                         if( t.lt.2*dt )then
                           write(*,'("advBA: advance BA dim=3, order=4 
     & grid=rectangular polar=NONE... t=",e10.2)') t
                         end if
                         mr=0
                         if( .not.methodOfLines )then
                           ! --- TAYLOR TIME-STEPPING --- 
                           stop 111
                         else
                           ! --- METHOD OF LINES (RK) ---
                             do i3=n3a,n3b
                             do i2=n2a,n2b
                             do i1=n1a,n1b
                               if( mask(i1,i2,i3).gt.0 )then
                             if( addForcing.ne.0 )then  ! do this for now *fix me*
                                 do m=0,5
                                  fv(m)=f(i1,i2,i3,m)
                                 end do
                             end if
                             if( numberOfMaterialRegions.gt.1 )then
                               mr = matMask(i1,i2,i3)
                               if( mr.lt.0 .or. 
     & mr.ge.numberOfMaterialRegions )then  ! do this for now
                                  stop 9999
                               end if
                             end if
                              ! --- 3D -----
                               curl(0) =   uy43r(i1,i2,i3,hz)- uz43r(
     & i1,i2,i3,hy)
                               curl(1) =   uz43r(i1,i2,i3,hx)- ux43r(
     & i1,i2,i3,hz)
                               curl(2) =   ux43r(i1,i2,i3,hy)- uy43r(
     & i1,i2,i3,hx)
                               curl(3) =-( uy43r(i1,i2,i3,ez)- uz43r(
     & i1,i2,i3,ey))
                               curl(4) =-( uz43r(i1,i2,i3,ex)- ux43r(
     & i1,i2,i3,ez))
                               curl(5) =-( ux43r(i1,i2,i3,ey)- uy43r(
     & i1,i2,i3,ex))
                               do m=0,5
                                 un(i1,i2,i3,m)= K0i(m,0,mr)*(curl(0)+
     & fv(0)) + K0i(m,1,mr)*(curl(1)+fv(1)) + K0i(m,2,mr)*(curl(2)+fv(
     & 2)) + K0i(m,3,mr)*(curl(3)+fv(3)) + K0i(m,4,mr)*(curl(4)+fv(4))
     &  + K0i(m,5,mr)*(curl(5)+fv(5))
                               end do
                               end if
                             end do
                             end do
                             end do
                           if( .false. .or. debug > 15 )then
                             stop  4747
                          end if
                         end if ! end MOL
                     end if
               end if
             else
               if( useSuperGrid.eq.0  ) then
                   if( t.lt.2*dt )then
                     write(*,'("advBA: advance BA GDM dim=3 order=4 
     & grid=rectangular polar=NONE... t=",e10.2)') t
                   end if
                   ! ---- Precompute some indirection arrays to make the dispersion loops go faster ----
                   do mr=0,numberOfMaterialRegions-1
                     m=0
                     pc=0
                     qc=1
                     do k1=1,6
                       ec=k1-1 ! E or H component
                       do k2=1,6
                         do n=1,Np(k1,k2,mr)
                           m=m+1
                           ecIndex1(m,mr)=ec
                           qcIndex1(m,mr)=qc
                           ! ptSum(ec) = ptSum(ec) + p(i1,i2,i3,qc) - pv(qc)
                           qc=qc+2
                         end do
                       end do
                       ! subtract off P.t = sum Q_m
                       ! curl(ec) = curl(ec) - ptSum(ec)
                     end do
                     numTerms1(mr)=m
                     if( numTerms1(mr).gt.maxNumPolarizationTerms )then
                       write(*,'(" ERROR numTerms1=",i6," too big")') 
     & numTerms1(mr)
                       stop 1616
                     end if
                     pc=0  ! P
                     qc=1  ! Q = P.t
                     m=0
                     do k1=1,6
                       do k2=1,6
                         ec=k2-1 ! This GDM term involves this E or H component
                         do n=1,Np(k1,k2,mr)
                           a0 = gdmPar(1,n,k1,k2,mr)
                           a1 = gdmPar(2,n,k1,k2,mr)
                           b0 = gdmPar(3,n,k1,k2,mr)
                           b1 = gdmPar(4,n,k1,k2,mr)
                           pct=pc+numComp  ! p.t is stored in un here
                           qct=qc+numComp  ! q.t is stored in un here
                           m=m+1
                           ecIndex2(m,mr)=ec
                           pcIndex2(m,mr)=pc
                           a0v(m,mr)=a0
                           a1v(m,mr)=a1
                           b0v(m,mr)=b0
                           b1v(m,mr)=b1
                           ! un(i1,i2,i3,pct) = p(i1,i2,i3,qc) + fp(pc) 
                           ! un(i1,i2,i3,qct) = a0*u(i1,i2,i3,ec) + a1*un(i1,i2,i3,ec) - b0*p(i1,i2,i3,pc)- b1*p(i1,i2,i3,qc) + fp(qc) 
                           pc=pc+2
                           qc=qc+2
                         end do
                       end do
                     end do
                     numTerms2(mr)=m
                     if( numTerms2(mr).gt.maxNumPolarizationTerms )then
                       write(*,'(" ERROR numTerms2=",i6," too big")') 
     & numTerms2(mr)
                       stop 1616
                     end if
                   end do ! end do mr
                   mr=0
                   if( .not.methodOfLines )then
                     ! --- TAYLOR TIME-STEPPING --- 
                     stop 111
                   else
                     ! --- METHOD OF LINES (RK) ---
                     ! zero out some forcing terms 
                     do m=0,numPolarizationTerms-1
                        pv(m)=0.
                       fp(m)=0.
                     end do
                       do i3=n3a,n3b
                       do i2=n2a,n2b
                       do i1=n1a,n1b
                         if( mask(i1,i2,i3).gt.0 )then
                       if( addForcing.ne.0 )then  ! do this for now *fix me*
                           do m=0,5
                            fv(m)=f(i1,i2,i3,m)
                           end do
                       end if
                       if( numberOfMaterialRegions.gt.1 )then
                         mr = matMask(i1,i2,i3)
                         if( mr.lt.0 .or. 
     & mr.ge.numberOfMaterialRegions )then  ! do this for now
                            stop 9999
                         end if
                       end if
                       if( forcingOption.eq.twilightZoneForcing )then
                         if( nd.eq.2 )then
                           do m=0,numComp-1
                             ec = m
                               call ogDeriv(ep, 0,0,0,0, xy(i1,i2,i3,0)
     & ,xy(i1,i2,i3,1),0.,t, ec,ev(m) )
                               call ogDeriv(ep, 1,0,0,0, xy(i1,i2,i3,0)
     & ,xy(i1,i2,i3,1),0.,t, ec,evt(m) )
                           end do
                           ! eval the polarization terms and time derivatives 
                           do m=0,numPolarizationTerms-1
                             pc = m+numComp  ! TZ index
                               call ogDeriv(ep, 0,0,0,0, xy(i1,i2,i3,0)
     & ,xy(i1,i2,i3,1),0.,t, pc,pv(m) )
                               call ogDeriv(ep, 1,0,0,0, xy(i1,i2,i3,0)
     & ,xy(i1,i2,i3,1),0.,t, pc,pvt(m) )
                           end do
                         else
                           do m=0,numComp-1
                             ec = m
                               call ogDeriv(ep, 0,0,0,0, xy(i1,i2,i3,0)
     & ,xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,ev(m) )
                               call ogDeriv(ep, 1,0,0,0, xy(i1,i2,i3,0)
     & ,xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,evt(m) )
                           end do
                           ! eval the polarization terms and time derivatives 
                           do m=0,numPolarizationTerms-1
                             pc = m+numComp  ! TZ index
                               call ogDeriv(ep, 0,0,0,0, xy(i1,i2,i3,0)
     & ,xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, pc,pv(m) )
                               call ogDeriv(ep, 1,0,0,0, xy(i1,i2,i3,0)
     & ,xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, pc,pvt(m) )
                           end do
                        end if
                        !write(*,'(" i1,i2=",2i3," xy=",2(f5.2,1x))') i1,i2,xy(i1,i2,i3,0),xy(i1,i2,i3,1)
                        !write(*,'(" i1,i2=",2i3," ev=",6(f6.3,1x))') i1,i2,(ev(m),m=0,5)
                        !write(*,'(" i1,i2=",2i3," pv=",10(f6.3,1x))') i1,i2,(pv(m),m=0,numPolarizationTerms-1)
                         ! FIX ME -- use opt version here too
                         pc=0  ! P
                         qc=1  ! Q = P.t
                         do k1=1,6
                           do k2=1,6
                             ec=k2-1 ! E or H component
                             do n=1,Np(k1,k2,mr)
                               a0 = gdmPar(1,n,k1,k2,mr)
                               a1 = gdmPar(2,n,k1,k2,mr)
                               b0 = gdmPar(3,n,k1,k2,mr)
                               b1 = gdmPar(4,n,k1,k2,mr)
                               ! TZ forcing for polarization equations: 
                               fp(pc) = pvt(pc) - pv(qc)
                               fp(qc) = pvt(qc) - (  a0*ev(ec) + a1*
     & evt(ec) - b0*pv(pc)- b1*pv(qc) )
                               pc=pc+2
                               qc=qc+2
                             end do
                           end do
                         end do
                      end if
                       ! compute components of the curl(H) and -curl(E)
                        ! --- 3D -----
                         curl(0) =   uy43r(i1,i2,i3,hz)- uz43r(i1,i2,
     & i3,hy)
                         curl(1) =   uz43r(i1,i2,i3,hx)- ux43r(i1,i2,
     & i3,hz)
                         curl(2) =   ux43r(i1,i2,i3,hy)- uy43r(i1,i2,
     & i3,hx)
                         curl(3) =-( uy43r(i1,i2,i3,ez)- uz43r(i1,i2,
     & i3,ey))
                         curl(4) =-( uz43r(i1,i2,i3,ex)- ux43r(i1,i2,
     & i3,ez))
                         curl(5) =-( ux43r(i1,i2,i3,ey)- uy43r(i1,i2,
     & i3,ex))
                       if( debug.gt.3 )then
                         write(*,'("----- i1,i2=",2i3)') i1,i2
                       end if
                       ! ---- Compute q = p.t = SUM_k2 SUM_n  p(i1,i2,i3, n,k1,k2, qc )
                       ! opt version 
                         do m=0,5
                           ptSum(m)=0
                         end do
                       do m=1,numTerms1(mr)
                          ec = ecIndex1(m,mr)
                          qc = qcIndex1(m,mr)
                          ! upiv(m) = p(i1,i2,i3,qc-1)  ! didn't seem to help
                          ! uqiv(m) = p(i1,i2,i3,qc )
                          ! uptSum(ec) = ptSum(ec) + qiv(m) - pv(qc)
                          ptSum(ec) = ptSum(ec) + p(i1,i2,i3,qc) - pv(
     & qc)
                       end do
                         do m=0,5
                           curl(m) = curl(m) - ptSum(m)
                         end do
                       ! if( debug.gt.3 )then
                       !   write(*,'("   ptSum=",6(f6.3,1x))') (ptSum(m),m=0,5)
                       ! end if
                         do m=0,5
                           un(i1,i2,i3,m)= K0i(m,0,mr)*(curl(0)+fv(0)) 
     & + K0i(m,1,mr)*(curl(1)+fv(1)) + K0i(m,2,mr)*(curl(2)+fv(2)) + 
     & K0i(m,3,mr)*(curl(3)+fv(3)) + K0i(m,4,mr)*(curl(4)+fv(4)) + 
     & K0i(m,5,mr)*(curl(5)+fv(5))
                         end do
                       ! if( .false. )then
                       !   write(*,'(" i1,i2=",2i3," ut=",6(f6.3,1x))') i1,i2,(un(i1,i2,i3,m),m=0,5)
                       ! end if
                       ! --- compute time derivatives of P and Q
                       ! p.t = q
                       ! q.t = a0*E + a1*Et - b0*p - b1*q   
                       ! or 
                       ! q.t = a0*H + a1*Ht - b0*p - b1*q   
                       ! optimized version
                       do m=1,numTerms2(mr)
                         ec = ecIndex2(m,mr)
                         pc = pcIndex2(m,mr)
                         qc = pc+1
                         pct=pc+numComp  ! p.t is stored in un here
                         qct=qc+numComp  ! q.t is stored in un here
                         a0 = a0v(m,mr)
                         a1 = a1v(m,mr)
                         b0 = b0v(m,mr)
                         b1 = b1v(m,mr)
                         un(i1,i2,i3,pct) = p(i1,i2,i3,qc) + fp(pc)
                         un(i1,i2,i3,qct) = a0*u(i1,i2,i3,ec) + a1*un(
     & i1,i2,i3,ec) - b0*p(i1,i2,i3,pc)- b1*p(i1,i2,i3,qc) + fp(qc)
                         ! uun(i1,i2,i3,pct) = qiv(m) + fp(pc) 
                         ! uun(i1,i2,i3,qct) = a0*uv(ec) + a1*unv(ec) - b0*piv(m)- b1*qiv(m) + fp(qc) 
                       end do
                         end if
                       end do
                       end do
                       end do
                     if( .false. .or. debug > 15 )then
                       stop  4444
                    end if
                   end if
               else
                 ! --- SUPERGRID ---
                 if( t.le.3*dt )then
                   write(*,'(" USE SUPERGRID...")' )
                 end if
                     !  --- THREE-DIMENSIONS ---
                     if( useAbsorbingLayer(0).eq.1 .and. 
     & useAbsorbingLayer(1).eq.1 .and. useAbsorbingLayer(2).eq.1 )then
                         if( t.lt.2*dt )then
                           write(*,'("advBA: advance BA GDM dim=3 
     & order=4 grid=rectangular polar=NONE... t=",e10.2)') t
                         end if
                         ! ---- Precompute some indirection arrays to make the dispersion loops go faster ----
                         do mr=0,numberOfMaterialRegions-1
                           m=0
                           pc=0
                           qc=1
                           do k1=1,6
                             ec=k1-1 ! E or H component
                             do k2=1,6
                               do n=1,Np(k1,k2,mr)
                                 m=m+1
                                 ecIndex1(m,mr)=ec
                                 qcIndex1(m,mr)=qc
                                 ! ptSum(ec) = ptSum(ec) + p(i1,i2,i3,qc) - pv(qc)
                                 qc=qc+2
                               end do
                             end do
                             ! subtract off P.t = sum Q_m
                             ! curl(ec) = curl(ec) - ptSum(ec)
                           end do
                           numTerms1(mr)=m
                           if( numTerms1(mr)
     & .gt.maxNumPolarizationTerms )then
                             write(*,'(" ERROR numTerms1=",i6," too 
     & big")') numTerms1(mr)
                             stop 1616
                           end if
                           pc=0  ! P
                           qc=1  ! Q = P.t
                           m=0
                           do k1=1,6
                             do k2=1,6
                               ec=k2-1 ! This GDM term involves this E or H component
                               do n=1,Np(k1,k2,mr)
                                 a0 = gdmPar(1,n,k1,k2,mr)
                                 a1 = gdmPar(2,n,k1,k2,mr)
                                 b0 = gdmPar(3,n,k1,k2,mr)
                                 b1 = gdmPar(4,n,k1,k2,mr)
                                 pct=pc+numComp  ! p.t is stored in un here
                                 qct=qc+numComp  ! q.t is stored in un here
                                 m=m+1
                                 ecIndex2(m,mr)=ec
                                 pcIndex2(m,mr)=pc
                                 a0v(m,mr)=a0
                                 a1v(m,mr)=a1
                                 b0v(m,mr)=b0
                                 b1v(m,mr)=b1
                                 ! un(i1,i2,i3,pct) = p(i1,i2,i3,qc) + fp(pc) 
                                 ! un(i1,i2,i3,qct) = a0*u(i1,i2,i3,ec) + a1*un(i1,i2,i3,ec) - b0*p(i1,i2,i3,pc)- b1*p(i1,i2,i3,qc) + fp(qc) 
                                 pc=pc+2
                                 qc=qc+2
                               end do
                             end do
                           end do
                           numTerms2(mr)=m
                           if( numTerms2(mr)
     & .gt.maxNumPolarizationTerms )then
                             write(*,'(" ERROR numTerms2=",i6," too 
     & big")') numTerms2(mr)
                             stop 1616
                           end if
                         end do ! end do mr
                         mr=0
                         if( .not.methodOfLines )then
                           ! --- TAYLOR TIME-STEPPING --- 
                           stop 111
                         else
                           ! --- METHOD OF LINES (RK) ---
                           ! zero out some forcing terms 
                           do m=0,numPolarizationTerms-1
                              pv(m)=0.
                             fp(m)=0.
                           end do
                             do i3=n3a,n3b
                             do i2=n2a,n2b
                             do i1=n1a,n1b
                               if( mask(i1,i2,i3).gt.0 )then
                             if( addForcing.ne.0 )then  ! do this for now *fix me*
                                 do m=0,5
                                  fv(m)=f(i1,i2,i3,m)
                                 end do
                             end if
                             if( numberOfMaterialRegions.gt.1 )then
                               mr = matMask(i1,i2,i3)
                               if( mr.lt.0 .or. 
     & mr.ge.numberOfMaterialRegions )then  ! do this for now
                                  stop 9999
                               end if
                             end if
                             if( forcingOption.eq.twilightZoneForcing )
     & then
                               if( nd.eq.2 )then
                                 do m=0,numComp-1
                                   ec = m
                                     call ogDeriv(ep, 0,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,ev(m) )
                                     call ogDeriv(ep, 1,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,evt(m) )
                                 end do
                                 ! eval the polarization terms and time derivatives 
                                 do m=0,numPolarizationTerms-1
                                   pc = m+numComp  ! TZ index
                                     call ogDeriv(ep, 0,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),0.,t, pc,pv(m) )
                                     call ogDeriv(ep, 1,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),0.,t, pc,pvt(m) )
                                 end do
                               else
                                 do m=0,numComp-1
                                   ec = m
                                     call ogDeriv(ep, 0,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,ev(m) )
                                     call ogDeriv(ep, 1,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,evt(m) )
                                 end do
                                 ! eval the polarization terms and time derivatives 
                                 do m=0,numPolarizationTerms-1
                                   pc = m+numComp  ! TZ index
                                     call ogDeriv(ep, 0,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, pc,pv(m) )
                                     call ogDeriv(ep, 1,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, pc,pvt(m) )
                                 end do
                              end if
                              !write(*,'(" i1,i2=",2i3," xy=",2(f5.2,1x))') i1,i2,xy(i1,i2,i3,0),xy(i1,i2,i3,1)
                              !write(*,'(" i1,i2=",2i3," ev=",6(f6.3,1x))') i1,i2,(ev(m),m=0,5)
                              !write(*,'(" i1,i2=",2i3," pv=",10(f6.3,1x))') i1,i2,(pv(m),m=0,numPolarizationTerms-1)
                               ! FIX ME -- use opt version here too
                               pc=0  ! P
                               qc=1  ! Q = P.t
                               do k1=1,6
                                 do k2=1,6
                                   ec=k2-1 ! E or H component
                                   do n=1,Np(k1,k2,mr)
                                     a0 = gdmPar(1,n,k1,k2,mr)
                                     a1 = gdmPar(2,n,k1,k2,mr)
                                     b0 = gdmPar(3,n,k1,k2,mr)
                                     b1 = gdmPar(4,n,k1,k2,mr)
                                     ! TZ forcing for polarization equations: 
                                     fp(pc) = pvt(pc) - pv(qc)
                                     fp(qc) = pvt(qc) - (  a0*ev(ec) + 
     & a1*evt(ec) - b0*pv(pc)- b1*pv(qc) )
                                     pc=pc+2
                                     qc=qc+2
                                   end do
                                 end do
                               end do
                            end if
                             ! compute components of the curl(H) and -curl(E)
                              ! --- 3D -----
                               curl(0) =  etay(i2)* uy43r(i1,i2,i3,hz)-
     & etaz(i3)* uz43r(i1,i2,i3,hy)
                               curl(1) =  etaz(i3)* uz43r(i1,i2,i3,hx)-
     & etax(i1)* ux43r(i1,i2,i3,hz)
                               curl(2) =  etax(i1)* ux43r(i1,i2,i3,hy)-
     & etay(i2)* uy43r(i1,i2,i3,hx)
                               curl(3) =-(etay(i2)* uy43r(i1,i2,i3,ez)-
     & etaz(i3)* uz43r(i1,i2,i3,ey))
                               curl(4) =-(etaz(i3)* uz43r(i1,i2,i3,ex)-
     & etax(i1)* ux43r(i1,i2,i3,ez))
                               curl(5) =-(etax(i1)* ux43r(i1,i2,i3,ey)-
     & etay(i2)* uy43r(i1,i2,i3,ex))
                             if( debug.gt.3 )then
                               write(*,'("----- i1,i2=",2i3)') i1,i2
                             end if
                             ! ---- Compute q = p.t = SUM_k2 SUM_n  p(i1,i2,i3, n,k1,k2, qc )
                             ! opt version 
                               do m=0,5
                                 ptSum(m)=0
                               end do
                             do m=1,numTerms1(mr)
                                ec = ecIndex1(m,mr)
                                qc = qcIndex1(m,mr)
                                ! upiv(m) = p(i1,i2,i3,qc-1)  ! didn't seem to help
                                ! uqiv(m) = p(i1,i2,i3,qc )
                                ! uptSum(ec) = ptSum(ec) + qiv(m) - pv(qc)
                                ptSum(ec) = ptSum(ec) + p(i1,i2,i3,qc) 
     & - pv(qc)
                             end do
                               do m=0,5
                                 curl(m) = curl(m) - ptSum(m)
                               end do
                             ! if( debug.gt.3 )then
                             !   write(*,'("   ptSum=",6(f6.3,1x))') (ptSum(m),m=0,5)
                             ! end if
                               do m=0,5
                                 un(i1,i2,i3,m)= K0i(m,0,mr)*(curl(0)+
     & fv(0)) + K0i(m,1,mr)*(curl(1)+fv(1)) + K0i(m,2,mr)*(curl(2)+fv(
     & 2)) + K0i(m,3,mr)*(curl(3)+fv(3)) + K0i(m,4,mr)*(curl(4)+fv(4))
     &  + K0i(m,5,mr)*(curl(5)+fv(5))
                               end do
                             ! if( .false. )then
                             !   write(*,'(" i1,i2=",2i3," ut=",6(f6.3,1x))') i1,i2,(un(i1,i2,i3,m),m=0,5)
                             ! end if
                             ! --- compute time derivatives of P and Q
                             ! p.t = q
                             ! q.t = a0*E + a1*Et - b0*p - b1*q   
                             ! or 
                             ! q.t = a0*H + a1*Ht - b0*p - b1*q   
                             ! optimized version
                             do m=1,numTerms2(mr)
                               ec = ecIndex2(m,mr)
                               pc = pcIndex2(m,mr)
                               qc = pc+1
                               pct=pc+numComp  ! p.t is stored in un here
                               qct=qc+numComp  ! q.t is stored in un here
                               a0 = a0v(m,mr)
                               a1 = a1v(m,mr)
                               b0 = b0v(m,mr)
                               b1 = b1v(m,mr)
                               un(i1,i2,i3,pct) = p(i1,i2,i3,qc) + fp(
     & pc)
                               un(i1,i2,i3,qct) = a0*u(i1,i2,i3,ec) + 
     & a1*un(i1,i2,i3,ec) - b0*p(i1,i2,i3,pc)- b1*p(i1,i2,i3,qc) + fp(
     & qc)
                               ! uun(i1,i2,i3,pct) = qiv(m) + fp(pc) 
                               ! uun(i1,i2,i3,qct) = a0*uv(ec) + a1*unv(ec) - b0*piv(m)- b1*qiv(m) + fp(qc) 
                             end do
                               end if
                             end do
                             end do
                             end do
                           if( .false. .or. debug > 15 )then
                             stop  4444
                          end if
                         end if
                     else if( useAbsorbingLayer(0).eq.1 .and. 
     & useAbsorbingLayer(1).eq.1 .and. useAbsorbingLayer(2).eq.0 )then
                         if( t.lt.2*dt )then
                           write(*,'("advBA: advance BA GDM dim=3 
     & order=4 grid=rectangular polar=NONE... t=",e10.2)') t
                         end if
                         ! ---- Precompute some indirection arrays to make the dispersion loops go faster ----
                         do mr=0,numberOfMaterialRegions-1
                           m=0
                           pc=0
                           qc=1
                           do k1=1,6
                             ec=k1-1 ! E or H component
                             do k2=1,6
                               do n=1,Np(k1,k2,mr)
                                 m=m+1
                                 ecIndex1(m,mr)=ec
                                 qcIndex1(m,mr)=qc
                                 ! ptSum(ec) = ptSum(ec) + p(i1,i2,i3,qc) - pv(qc)
                                 qc=qc+2
                               end do
                             end do
                             ! subtract off P.t = sum Q_m
                             ! curl(ec) = curl(ec) - ptSum(ec)
                           end do
                           numTerms1(mr)=m
                           if( numTerms1(mr)
     & .gt.maxNumPolarizationTerms )then
                             write(*,'(" ERROR numTerms1=",i6," too 
     & big")') numTerms1(mr)
                             stop 1616
                           end if
                           pc=0  ! P
                           qc=1  ! Q = P.t
                           m=0
                           do k1=1,6
                             do k2=1,6
                               ec=k2-1 ! This GDM term involves this E or H component
                               do n=1,Np(k1,k2,mr)
                                 a0 = gdmPar(1,n,k1,k2,mr)
                                 a1 = gdmPar(2,n,k1,k2,mr)
                                 b0 = gdmPar(3,n,k1,k2,mr)
                                 b1 = gdmPar(4,n,k1,k2,mr)
                                 pct=pc+numComp  ! p.t is stored in un here
                                 qct=qc+numComp  ! q.t is stored in un here
                                 m=m+1
                                 ecIndex2(m,mr)=ec
                                 pcIndex2(m,mr)=pc
                                 a0v(m,mr)=a0
                                 a1v(m,mr)=a1
                                 b0v(m,mr)=b0
                                 b1v(m,mr)=b1
                                 ! un(i1,i2,i3,pct) = p(i1,i2,i3,qc) + fp(pc) 
                                 ! un(i1,i2,i3,qct) = a0*u(i1,i2,i3,ec) + a1*un(i1,i2,i3,ec) - b0*p(i1,i2,i3,pc)- b1*p(i1,i2,i3,qc) + fp(qc) 
                                 pc=pc+2
                                 qc=qc+2
                               end do
                             end do
                           end do
                           numTerms2(mr)=m
                           if( numTerms2(mr)
     & .gt.maxNumPolarizationTerms )then
                             write(*,'(" ERROR numTerms2=",i6," too 
     & big")') numTerms2(mr)
                             stop 1616
                           end if
                         end do ! end do mr
                         mr=0
                         if( .not.methodOfLines )then
                           ! --- TAYLOR TIME-STEPPING --- 
                           stop 111
                         else
                           ! --- METHOD OF LINES (RK) ---
                           ! zero out some forcing terms 
                           do m=0,numPolarizationTerms-1
                              pv(m)=0.
                             fp(m)=0.
                           end do
                             do i3=n3a,n3b
                             do i2=n2a,n2b
                             do i1=n1a,n1b
                               if( mask(i1,i2,i3).gt.0 )then
                             if( addForcing.ne.0 )then  ! do this for now *fix me*
                                 do m=0,5
                                  fv(m)=f(i1,i2,i3,m)
                                 end do
                             end if
                             if( numberOfMaterialRegions.gt.1 )then
                               mr = matMask(i1,i2,i3)
                               if( mr.lt.0 .or. 
     & mr.ge.numberOfMaterialRegions )then  ! do this for now
                                  stop 9999
                               end if
                             end if
                             if( forcingOption.eq.twilightZoneForcing )
     & then
                               if( nd.eq.2 )then
                                 do m=0,numComp-1
                                   ec = m
                                     call ogDeriv(ep, 0,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,ev(m) )
                                     call ogDeriv(ep, 1,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,evt(m) )
                                 end do
                                 ! eval the polarization terms and time derivatives 
                                 do m=0,numPolarizationTerms-1
                                   pc = m+numComp  ! TZ index
                                     call ogDeriv(ep, 0,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),0.,t, pc,pv(m) )
                                     call ogDeriv(ep, 1,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),0.,t, pc,pvt(m) )
                                 end do
                               else
                                 do m=0,numComp-1
                                   ec = m
                                     call ogDeriv(ep, 0,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,ev(m) )
                                     call ogDeriv(ep, 1,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,evt(m) )
                                 end do
                                 ! eval the polarization terms and time derivatives 
                                 do m=0,numPolarizationTerms-1
                                   pc = m+numComp  ! TZ index
                                     call ogDeriv(ep, 0,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, pc,pv(m) )
                                     call ogDeriv(ep, 1,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, pc,pvt(m) )
                                 end do
                              end if
                              !write(*,'(" i1,i2=",2i3," xy=",2(f5.2,1x))') i1,i2,xy(i1,i2,i3,0),xy(i1,i2,i3,1)
                              !write(*,'(" i1,i2=",2i3," ev=",6(f6.3,1x))') i1,i2,(ev(m),m=0,5)
                              !write(*,'(" i1,i2=",2i3," pv=",10(f6.3,1x))') i1,i2,(pv(m),m=0,numPolarizationTerms-1)
                               ! FIX ME -- use opt version here too
                               pc=0  ! P
                               qc=1  ! Q = P.t
                               do k1=1,6
                                 do k2=1,6
                                   ec=k2-1 ! E or H component
                                   do n=1,Np(k1,k2,mr)
                                     a0 = gdmPar(1,n,k1,k2,mr)
                                     a1 = gdmPar(2,n,k1,k2,mr)
                                     b0 = gdmPar(3,n,k1,k2,mr)
                                     b1 = gdmPar(4,n,k1,k2,mr)
                                     ! TZ forcing for polarization equations: 
                                     fp(pc) = pvt(pc) - pv(qc)
                                     fp(qc) = pvt(qc) - (  a0*ev(ec) + 
     & a1*evt(ec) - b0*pv(pc)- b1*pv(qc) )
                                     pc=pc+2
                                     qc=qc+2
                                   end do
                                 end do
                               end do
                            end if
                             ! compute components of the curl(H) and -curl(E)
                              ! --- 3D -----
                               curl(0) =  etay(i2)* uy43r(i1,i2,i3,hz)-
     &  uz43r(i1,i2,i3,hy)
                               curl(1) =   uz43r(i1,i2,i3,hx)-etax(i1)*
     &  ux43r(i1,i2,i3,hz)
                               curl(2) =  etax(i1)* ux43r(i1,i2,i3,hy)-
     & etay(i2)* uy43r(i1,i2,i3,hx)
                               curl(3) =-(etay(i2)* uy43r(i1,i2,i3,ez)-
     &  uz43r(i1,i2,i3,ey))
                               curl(4) =-( uz43r(i1,i2,i3,ex)-etax(i1)*
     &  ux43r(i1,i2,i3,ez))
                               curl(5) =-(etax(i1)* ux43r(i1,i2,i3,ey)-
     & etay(i2)* uy43r(i1,i2,i3,ex))
                             if( debug.gt.3 )then
                               write(*,'("----- i1,i2=",2i3)') i1,i2
                             end if
                             ! ---- Compute q = p.t = SUM_k2 SUM_n  p(i1,i2,i3, n,k1,k2, qc )
                             ! opt version 
                               do m=0,5
                                 ptSum(m)=0
                               end do
                             do m=1,numTerms1(mr)
                                ec = ecIndex1(m,mr)
                                qc = qcIndex1(m,mr)
                                ! upiv(m) = p(i1,i2,i3,qc-1)  ! didn't seem to help
                                ! uqiv(m) = p(i1,i2,i3,qc )
                                ! uptSum(ec) = ptSum(ec) + qiv(m) - pv(qc)
                                ptSum(ec) = ptSum(ec) + p(i1,i2,i3,qc) 
     & - pv(qc)
                             end do
                               do m=0,5
                                 curl(m) = curl(m) - ptSum(m)
                               end do
                             ! if( debug.gt.3 )then
                             !   write(*,'("   ptSum=",6(f6.3,1x))') (ptSum(m),m=0,5)
                             ! end if
                               do m=0,5
                                 un(i1,i2,i3,m)= K0i(m,0,mr)*(curl(0)+
     & fv(0)) + K0i(m,1,mr)*(curl(1)+fv(1)) + K0i(m,2,mr)*(curl(2)+fv(
     & 2)) + K0i(m,3,mr)*(curl(3)+fv(3)) + K0i(m,4,mr)*(curl(4)+fv(4))
     &  + K0i(m,5,mr)*(curl(5)+fv(5))
                               end do
                             ! if( .false. )then
                             !   write(*,'(" i1,i2=",2i3," ut=",6(f6.3,1x))') i1,i2,(un(i1,i2,i3,m),m=0,5)
                             ! end if
                             ! --- compute time derivatives of P and Q
                             ! p.t = q
                             ! q.t = a0*E + a1*Et - b0*p - b1*q   
                             ! or 
                             ! q.t = a0*H + a1*Ht - b0*p - b1*q   
                             ! optimized version
                             do m=1,numTerms2(mr)
                               ec = ecIndex2(m,mr)
                               pc = pcIndex2(m,mr)
                               qc = pc+1
                               pct=pc+numComp  ! p.t is stored in un here
                               qct=qc+numComp  ! q.t is stored in un here
                               a0 = a0v(m,mr)
                               a1 = a1v(m,mr)
                               b0 = b0v(m,mr)
                               b1 = b1v(m,mr)
                               un(i1,i2,i3,pct) = p(i1,i2,i3,qc) + fp(
     & pc)
                               un(i1,i2,i3,qct) = a0*u(i1,i2,i3,ec) + 
     & a1*un(i1,i2,i3,ec) - b0*p(i1,i2,i3,pc)- b1*p(i1,i2,i3,qc) + fp(
     & qc)
                               ! uun(i1,i2,i3,pct) = qiv(m) + fp(pc) 
                               ! uun(i1,i2,i3,qct) = a0*uv(ec) + a1*unv(ec) - b0*piv(m)- b1*qiv(m) + fp(qc) 
                             end do
                               end if
                             end do
                             end do
                             end do
                           if( .false. .or. debug > 15 )then
                             stop  4444
                          end if
                         end if
                     else if( useAbsorbingLayer(0).eq.1 .and. 
     & useAbsorbingLayer(1).eq.0 .and. useAbsorbingLayer(2).eq.1 )then
                         if( t.lt.2*dt )then
                           write(*,'("advBA: advance BA GDM dim=3 
     & order=4 grid=rectangular polar=NONE... t=",e10.2)') t
                         end if
                         ! ---- Precompute some indirection arrays to make the dispersion loops go faster ----
                         do mr=0,numberOfMaterialRegions-1
                           m=0
                           pc=0
                           qc=1
                           do k1=1,6
                             ec=k1-1 ! E or H component
                             do k2=1,6
                               do n=1,Np(k1,k2,mr)
                                 m=m+1
                                 ecIndex1(m,mr)=ec
                                 qcIndex1(m,mr)=qc
                                 ! ptSum(ec) = ptSum(ec) + p(i1,i2,i3,qc) - pv(qc)
                                 qc=qc+2
                               end do
                             end do
                             ! subtract off P.t = sum Q_m
                             ! curl(ec) = curl(ec) - ptSum(ec)
                           end do
                           numTerms1(mr)=m
                           if( numTerms1(mr)
     & .gt.maxNumPolarizationTerms )then
                             write(*,'(" ERROR numTerms1=",i6," too 
     & big")') numTerms1(mr)
                             stop 1616
                           end if
                           pc=0  ! P
                           qc=1  ! Q = P.t
                           m=0
                           do k1=1,6
                             do k2=1,6
                               ec=k2-1 ! This GDM term involves this E or H component
                               do n=1,Np(k1,k2,mr)
                                 a0 = gdmPar(1,n,k1,k2,mr)
                                 a1 = gdmPar(2,n,k1,k2,mr)
                                 b0 = gdmPar(3,n,k1,k2,mr)
                                 b1 = gdmPar(4,n,k1,k2,mr)
                                 pct=pc+numComp  ! p.t is stored in un here
                                 qct=qc+numComp  ! q.t is stored in un here
                                 m=m+1
                                 ecIndex2(m,mr)=ec
                                 pcIndex2(m,mr)=pc
                                 a0v(m,mr)=a0
                                 a1v(m,mr)=a1
                                 b0v(m,mr)=b0
                                 b1v(m,mr)=b1
                                 ! un(i1,i2,i3,pct) = p(i1,i2,i3,qc) + fp(pc) 
                                 ! un(i1,i2,i3,qct) = a0*u(i1,i2,i3,ec) + a1*un(i1,i2,i3,ec) - b0*p(i1,i2,i3,pc)- b1*p(i1,i2,i3,qc) + fp(qc) 
                                 pc=pc+2
                                 qc=qc+2
                               end do
                             end do
                           end do
                           numTerms2(mr)=m
                           if( numTerms2(mr)
     & .gt.maxNumPolarizationTerms )then
                             write(*,'(" ERROR numTerms2=",i6," too 
     & big")') numTerms2(mr)
                             stop 1616
                           end if
                         end do ! end do mr
                         mr=0
                         if( .not.methodOfLines )then
                           ! --- TAYLOR TIME-STEPPING --- 
                           stop 111
                         else
                           ! --- METHOD OF LINES (RK) ---
                           ! zero out some forcing terms 
                           do m=0,numPolarizationTerms-1
                              pv(m)=0.
                             fp(m)=0.
                           end do
                             do i3=n3a,n3b
                             do i2=n2a,n2b
                             do i1=n1a,n1b
                               if( mask(i1,i2,i3).gt.0 )then
                             if( addForcing.ne.0 )then  ! do this for now *fix me*
                                 do m=0,5
                                  fv(m)=f(i1,i2,i3,m)
                                 end do
                             end if
                             if( numberOfMaterialRegions.gt.1 )then
                               mr = matMask(i1,i2,i3)
                               if( mr.lt.0 .or. 
     & mr.ge.numberOfMaterialRegions )then  ! do this for now
                                  stop 9999
                               end if
                             end if
                             if( forcingOption.eq.twilightZoneForcing )
     & then
                               if( nd.eq.2 )then
                                 do m=0,numComp-1
                                   ec = m
                                     call ogDeriv(ep, 0,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,ev(m) )
                                     call ogDeriv(ep, 1,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,evt(m) )
                                 end do
                                 ! eval the polarization terms and time derivatives 
                                 do m=0,numPolarizationTerms-1
                                   pc = m+numComp  ! TZ index
                                     call ogDeriv(ep, 0,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),0.,t, pc,pv(m) )
                                     call ogDeriv(ep, 1,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),0.,t, pc,pvt(m) )
                                 end do
                               else
                                 do m=0,numComp-1
                                   ec = m
                                     call ogDeriv(ep, 0,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,ev(m) )
                                     call ogDeriv(ep, 1,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,evt(m) )
                                 end do
                                 ! eval the polarization terms and time derivatives 
                                 do m=0,numPolarizationTerms-1
                                   pc = m+numComp  ! TZ index
                                     call ogDeriv(ep, 0,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, pc,pv(m) )
                                     call ogDeriv(ep, 1,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, pc,pvt(m) )
                                 end do
                              end if
                              !write(*,'(" i1,i2=",2i3," xy=",2(f5.2,1x))') i1,i2,xy(i1,i2,i3,0),xy(i1,i2,i3,1)
                              !write(*,'(" i1,i2=",2i3," ev=",6(f6.3,1x))') i1,i2,(ev(m),m=0,5)
                              !write(*,'(" i1,i2=",2i3," pv=",10(f6.3,1x))') i1,i2,(pv(m),m=0,numPolarizationTerms-1)
                               ! FIX ME -- use opt version here too
                               pc=0  ! P
                               qc=1  ! Q = P.t
                               do k1=1,6
                                 do k2=1,6
                                   ec=k2-1 ! E or H component
                                   do n=1,Np(k1,k2,mr)
                                     a0 = gdmPar(1,n,k1,k2,mr)
                                     a1 = gdmPar(2,n,k1,k2,mr)
                                     b0 = gdmPar(3,n,k1,k2,mr)
                                     b1 = gdmPar(4,n,k1,k2,mr)
                                     ! TZ forcing for polarization equations: 
                                     fp(pc) = pvt(pc) - pv(qc)
                                     fp(qc) = pvt(qc) - (  a0*ev(ec) + 
     & a1*evt(ec) - b0*pv(pc)- b1*pv(qc) )
                                     pc=pc+2
                                     qc=qc+2
                                   end do
                                 end do
                               end do
                            end if
                             ! compute components of the curl(H) and -curl(E)
                              ! --- 3D -----
                               curl(0) =   uy43r(i1,i2,i3,hz)-etaz(i3)*
     &  uz43r(i1,i2,i3,hy)
                               curl(1) =  etaz(i3)* uz43r(i1,i2,i3,hx)-
     & etax(i1)* ux43r(i1,i2,i3,hz)
                               curl(2) =  etax(i1)* ux43r(i1,i2,i3,hy)-
     &  uy43r(i1,i2,i3,hx)
                               curl(3) =-( uy43r(i1,i2,i3,ez)-etaz(i3)*
     &  uz43r(i1,i2,i3,ey))
                               curl(4) =-(etaz(i3)* uz43r(i1,i2,i3,ex)-
     & etax(i1)* ux43r(i1,i2,i3,ez))
                               curl(5) =-(etax(i1)* ux43r(i1,i2,i3,ey)-
     &  uy43r(i1,i2,i3,ex))
                             if( debug.gt.3 )then
                               write(*,'("----- i1,i2=",2i3)') i1,i2
                             end if
                             ! ---- Compute q = p.t = SUM_k2 SUM_n  p(i1,i2,i3, n,k1,k2, qc )
                             ! opt version 
                               do m=0,5
                                 ptSum(m)=0
                               end do
                             do m=1,numTerms1(mr)
                                ec = ecIndex1(m,mr)
                                qc = qcIndex1(m,mr)
                                ! upiv(m) = p(i1,i2,i3,qc-1)  ! didn't seem to help
                                ! uqiv(m) = p(i1,i2,i3,qc )
                                ! uptSum(ec) = ptSum(ec) + qiv(m) - pv(qc)
                                ptSum(ec) = ptSum(ec) + p(i1,i2,i3,qc) 
     & - pv(qc)
                             end do
                               do m=0,5
                                 curl(m) = curl(m) - ptSum(m)
                               end do
                             ! if( debug.gt.3 )then
                             !   write(*,'("   ptSum=",6(f6.3,1x))') (ptSum(m),m=0,5)
                             ! end if
                               do m=0,5
                                 un(i1,i2,i3,m)= K0i(m,0,mr)*(curl(0)+
     & fv(0)) + K0i(m,1,mr)*(curl(1)+fv(1)) + K0i(m,2,mr)*(curl(2)+fv(
     & 2)) + K0i(m,3,mr)*(curl(3)+fv(3)) + K0i(m,4,mr)*(curl(4)+fv(4))
     &  + K0i(m,5,mr)*(curl(5)+fv(5))
                               end do
                             ! if( .false. )then
                             !   write(*,'(" i1,i2=",2i3," ut=",6(f6.3,1x))') i1,i2,(un(i1,i2,i3,m),m=0,5)
                             ! end if
                             ! --- compute time derivatives of P and Q
                             ! p.t = q
                             ! q.t = a0*E + a1*Et - b0*p - b1*q   
                             ! or 
                             ! q.t = a0*H + a1*Ht - b0*p - b1*q   
                             ! optimized version
                             do m=1,numTerms2(mr)
                               ec = ecIndex2(m,mr)
                               pc = pcIndex2(m,mr)
                               qc = pc+1
                               pct=pc+numComp  ! p.t is stored in un here
                               qct=qc+numComp  ! q.t is stored in un here
                               a0 = a0v(m,mr)
                               a1 = a1v(m,mr)
                               b0 = b0v(m,mr)
                               b1 = b1v(m,mr)
                               un(i1,i2,i3,pct) = p(i1,i2,i3,qc) + fp(
     & pc)
                               un(i1,i2,i3,qct) = a0*u(i1,i2,i3,ec) + 
     & a1*un(i1,i2,i3,ec) - b0*p(i1,i2,i3,pc)- b1*p(i1,i2,i3,qc) + fp(
     & qc)
                               ! uun(i1,i2,i3,pct) = qiv(m) + fp(pc) 
                               ! uun(i1,i2,i3,qct) = a0*uv(ec) + a1*unv(ec) - b0*piv(m)- b1*qiv(m) + fp(qc) 
                             end do
                               end if
                             end do
                             end do
                             end do
                           if( .false. .or. debug > 15 )then
                             stop  4444
                          end if
                         end if
                     else if( useAbsorbingLayer(0).eq.0 .and. 
     & useAbsorbingLayer(1).eq.1 .and. useAbsorbingLayer(2).eq.1 )then
                         if( t.lt.2*dt )then
                           write(*,'("advBA: advance BA GDM dim=3 
     & order=4 grid=rectangular polar=NONE... t=",e10.2)') t
                         end if
                         ! ---- Precompute some indirection arrays to make the dispersion loops go faster ----
                         do mr=0,numberOfMaterialRegions-1
                           m=0
                           pc=0
                           qc=1
                           do k1=1,6
                             ec=k1-1 ! E or H component
                             do k2=1,6
                               do n=1,Np(k1,k2,mr)
                                 m=m+1
                                 ecIndex1(m,mr)=ec
                                 qcIndex1(m,mr)=qc
                                 ! ptSum(ec) = ptSum(ec) + p(i1,i2,i3,qc) - pv(qc)
                                 qc=qc+2
                               end do
                             end do
                             ! subtract off P.t = sum Q_m
                             ! curl(ec) = curl(ec) - ptSum(ec)
                           end do
                           numTerms1(mr)=m
                           if( numTerms1(mr)
     & .gt.maxNumPolarizationTerms )then
                             write(*,'(" ERROR numTerms1=",i6," too 
     & big")') numTerms1(mr)
                             stop 1616
                           end if
                           pc=0  ! P
                           qc=1  ! Q = P.t
                           m=0
                           do k1=1,6
                             do k2=1,6
                               ec=k2-1 ! This GDM term involves this E or H component
                               do n=1,Np(k1,k2,mr)
                                 a0 = gdmPar(1,n,k1,k2,mr)
                                 a1 = gdmPar(2,n,k1,k2,mr)
                                 b0 = gdmPar(3,n,k1,k2,mr)
                                 b1 = gdmPar(4,n,k1,k2,mr)
                                 pct=pc+numComp  ! p.t is stored in un here
                                 qct=qc+numComp  ! q.t is stored in un here
                                 m=m+1
                                 ecIndex2(m,mr)=ec
                                 pcIndex2(m,mr)=pc
                                 a0v(m,mr)=a0
                                 a1v(m,mr)=a1
                                 b0v(m,mr)=b0
                                 b1v(m,mr)=b1
                                 ! un(i1,i2,i3,pct) = p(i1,i2,i3,qc) + fp(pc) 
                                 ! un(i1,i2,i3,qct) = a0*u(i1,i2,i3,ec) + a1*un(i1,i2,i3,ec) - b0*p(i1,i2,i3,pc)- b1*p(i1,i2,i3,qc) + fp(qc) 
                                 pc=pc+2
                                 qc=qc+2
                               end do
                             end do
                           end do
                           numTerms2(mr)=m
                           if( numTerms2(mr)
     & .gt.maxNumPolarizationTerms )then
                             write(*,'(" ERROR numTerms2=",i6," too 
     & big")') numTerms2(mr)
                             stop 1616
                           end if
                         end do ! end do mr
                         mr=0
                         if( .not.methodOfLines )then
                           ! --- TAYLOR TIME-STEPPING --- 
                           stop 111
                         else
                           ! --- METHOD OF LINES (RK) ---
                           ! zero out some forcing terms 
                           do m=0,numPolarizationTerms-1
                              pv(m)=0.
                             fp(m)=0.
                           end do
                             do i3=n3a,n3b
                             do i2=n2a,n2b
                             do i1=n1a,n1b
                               if( mask(i1,i2,i3).gt.0 )then
                             if( addForcing.ne.0 )then  ! do this for now *fix me*
                                 do m=0,5
                                  fv(m)=f(i1,i2,i3,m)
                                 end do
                             end if
                             if( numberOfMaterialRegions.gt.1 )then
                               mr = matMask(i1,i2,i3)
                               if( mr.lt.0 .or. 
     & mr.ge.numberOfMaterialRegions )then  ! do this for now
                                  stop 9999
                               end if
                             end if
                             if( forcingOption.eq.twilightZoneForcing )
     & then
                               if( nd.eq.2 )then
                                 do m=0,numComp-1
                                   ec = m
                                     call ogDeriv(ep, 0,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,ev(m) )
                                     call ogDeriv(ep, 1,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,evt(m) )
                                 end do
                                 ! eval the polarization terms and time derivatives 
                                 do m=0,numPolarizationTerms-1
                                   pc = m+numComp  ! TZ index
                                     call ogDeriv(ep, 0,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),0.,t, pc,pv(m) )
                                     call ogDeriv(ep, 1,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),0.,t, pc,pvt(m) )
                                 end do
                               else
                                 do m=0,numComp-1
                                   ec = m
                                     call ogDeriv(ep, 0,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,ev(m) )
                                     call ogDeriv(ep, 1,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,evt(m) )
                                 end do
                                 ! eval the polarization terms and time derivatives 
                                 do m=0,numPolarizationTerms-1
                                   pc = m+numComp  ! TZ index
                                     call ogDeriv(ep, 0,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, pc,pv(m) )
                                     call ogDeriv(ep, 1,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, pc,pvt(m) )
                                 end do
                              end if
                              !write(*,'(" i1,i2=",2i3," xy=",2(f5.2,1x))') i1,i2,xy(i1,i2,i3,0),xy(i1,i2,i3,1)
                              !write(*,'(" i1,i2=",2i3," ev=",6(f6.3,1x))') i1,i2,(ev(m),m=0,5)
                              !write(*,'(" i1,i2=",2i3," pv=",10(f6.3,1x))') i1,i2,(pv(m),m=0,numPolarizationTerms-1)
                               ! FIX ME -- use opt version here too
                               pc=0  ! P
                               qc=1  ! Q = P.t
                               do k1=1,6
                                 do k2=1,6
                                   ec=k2-1 ! E or H component
                                   do n=1,Np(k1,k2,mr)
                                     a0 = gdmPar(1,n,k1,k2,mr)
                                     a1 = gdmPar(2,n,k1,k2,mr)
                                     b0 = gdmPar(3,n,k1,k2,mr)
                                     b1 = gdmPar(4,n,k1,k2,mr)
                                     ! TZ forcing for polarization equations: 
                                     fp(pc) = pvt(pc) - pv(qc)
                                     fp(qc) = pvt(qc) - (  a0*ev(ec) + 
     & a1*evt(ec) - b0*pv(pc)- b1*pv(qc) )
                                     pc=pc+2
                                     qc=qc+2
                                   end do
                                 end do
                               end do
                            end if
                             ! compute components of the curl(H) and -curl(E)
                              ! --- 3D -----
                               curl(0) =  etay(i2)* uy43r(i1,i2,i3,hz)-
     & etaz(i3)* uz43r(i1,i2,i3,hy)
                               curl(1) =  etaz(i3)* uz43r(i1,i2,i3,hx)-
     &  ux43r(i1,i2,i3,hz)
                               curl(2) =   ux43r(i1,i2,i3,hy)-etay(i2)*
     &  uy43r(i1,i2,i3,hx)
                               curl(3) =-(etay(i2)* uy43r(i1,i2,i3,ez)-
     & etaz(i3)* uz43r(i1,i2,i3,ey))
                               curl(4) =-(etaz(i3)* uz43r(i1,i2,i3,ex)-
     &  ux43r(i1,i2,i3,ez))
                               curl(5) =-( ux43r(i1,i2,i3,ey)-etay(i2)*
     &  uy43r(i1,i2,i3,ex))
                             if( debug.gt.3 )then
                               write(*,'("----- i1,i2=",2i3)') i1,i2
                             end if
                             ! ---- Compute q = p.t = SUM_k2 SUM_n  p(i1,i2,i3, n,k1,k2, qc )
                             ! opt version 
                               do m=0,5
                                 ptSum(m)=0
                               end do
                             do m=1,numTerms1(mr)
                                ec = ecIndex1(m,mr)
                                qc = qcIndex1(m,mr)
                                ! upiv(m) = p(i1,i2,i3,qc-1)  ! didn't seem to help
                                ! uqiv(m) = p(i1,i2,i3,qc )
                                ! uptSum(ec) = ptSum(ec) + qiv(m) - pv(qc)
                                ptSum(ec) = ptSum(ec) + p(i1,i2,i3,qc) 
     & - pv(qc)
                             end do
                               do m=0,5
                                 curl(m) = curl(m) - ptSum(m)
                               end do
                             ! if( debug.gt.3 )then
                             !   write(*,'("   ptSum=",6(f6.3,1x))') (ptSum(m),m=0,5)
                             ! end if
                               do m=0,5
                                 un(i1,i2,i3,m)= K0i(m,0,mr)*(curl(0)+
     & fv(0)) + K0i(m,1,mr)*(curl(1)+fv(1)) + K0i(m,2,mr)*(curl(2)+fv(
     & 2)) + K0i(m,3,mr)*(curl(3)+fv(3)) + K0i(m,4,mr)*(curl(4)+fv(4))
     &  + K0i(m,5,mr)*(curl(5)+fv(5))
                               end do
                             ! if( .false. )then
                             !   write(*,'(" i1,i2=",2i3," ut=",6(f6.3,1x))') i1,i2,(un(i1,i2,i3,m),m=0,5)
                             ! end if
                             ! --- compute time derivatives of P and Q
                             ! p.t = q
                             ! q.t = a0*E + a1*Et - b0*p - b1*q   
                             ! or 
                             ! q.t = a0*H + a1*Ht - b0*p - b1*q   
                             ! optimized version
                             do m=1,numTerms2(mr)
                               ec = ecIndex2(m,mr)
                               pc = pcIndex2(m,mr)
                               qc = pc+1
                               pct=pc+numComp  ! p.t is stored in un here
                               qct=qc+numComp  ! q.t is stored in un here
                               a0 = a0v(m,mr)
                               a1 = a1v(m,mr)
                               b0 = b0v(m,mr)
                               b1 = b1v(m,mr)
                               un(i1,i2,i3,pct) = p(i1,i2,i3,qc) + fp(
     & pc)
                               un(i1,i2,i3,qct) = a0*u(i1,i2,i3,ec) + 
     & a1*un(i1,i2,i3,ec) - b0*p(i1,i2,i3,pc)- b1*p(i1,i2,i3,qc) + fp(
     & qc)
                               ! uun(i1,i2,i3,pct) = qiv(m) + fp(pc) 
                               ! uun(i1,i2,i3,qct) = a0*uv(ec) + a1*unv(ec) - b0*piv(m)- b1*qiv(m) + fp(qc) 
                             end do
                               end if
                             end do
                             end do
                             end do
                           if( .false. .or. debug > 15 )then
                             stop  4444
                          end if
                         end if
                     else if( useAbsorbingLayer(0).eq.1 .and. 
     & useAbsorbingLayer(1).eq.0 .and. useAbsorbingLayer(2).eq.0 )then
                         if( t.lt.2*dt )then
                           write(*,'("advBA: advance BA GDM dim=3 
     & order=4 grid=rectangular polar=NONE... t=",e10.2)') t
                         end if
                         ! ---- Precompute some indirection arrays to make the dispersion loops go faster ----
                         do mr=0,numberOfMaterialRegions-1
                           m=0
                           pc=0
                           qc=1
                           do k1=1,6
                             ec=k1-1 ! E or H component
                             do k2=1,6
                               do n=1,Np(k1,k2,mr)
                                 m=m+1
                                 ecIndex1(m,mr)=ec
                                 qcIndex1(m,mr)=qc
                                 ! ptSum(ec) = ptSum(ec) + p(i1,i2,i3,qc) - pv(qc)
                                 qc=qc+2
                               end do
                             end do
                             ! subtract off P.t = sum Q_m
                             ! curl(ec) = curl(ec) - ptSum(ec)
                           end do
                           numTerms1(mr)=m
                           if( numTerms1(mr)
     & .gt.maxNumPolarizationTerms )then
                             write(*,'(" ERROR numTerms1=",i6," too 
     & big")') numTerms1(mr)
                             stop 1616
                           end if
                           pc=0  ! P
                           qc=1  ! Q = P.t
                           m=0
                           do k1=1,6
                             do k2=1,6
                               ec=k2-1 ! This GDM term involves this E or H component
                               do n=1,Np(k1,k2,mr)
                                 a0 = gdmPar(1,n,k1,k2,mr)
                                 a1 = gdmPar(2,n,k1,k2,mr)
                                 b0 = gdmPar(3,n,k1,k2,mr)
                                 b1 = gdmPar(4,n,k1,k2,mr)
                                 pct=pc+numComp  ! p.t is stored in un here
                                 qct=qc+numComp  ! q.t is stored in un here
                                 m=m+1
                                 ecIndex2(m,mr)=ec
                                 pcIndex2(m,mr)=pc
                                 a0v(m,mr)=a0
                                 a1v(m,mr)=a1
                                 b0v(m,mr)=b0
                                 b1v(m,mr)=b1
                                 ! un(i1,i2,i3,pct) = p(i1,i2,i3,qc) + fp(pc) 
                                 ! un(i1,i2,i3,qct) = a0*u(i1,i2,i3,ec) + a1*un(i1,i2,i3,ec) - b0*p(i1,i2,i3,pc)- b1*p(i1,i2,i3,qc) + fp(qc) 
                                 pc=pc+2
                                 qc=qc+2
                               end do
                             end do
                           end do
                           numTerms2(mr)=m
                           if( numTerms2(mr)
     & .gt.maxNumPolarizationTerms )then
                             write(*,'(" ERROR numTerms2=",i6," too 
     & big")') numTerms2(mr)
                             stop 1616
                           end if
                         end do ! end do mr
                         mr=0
                         if( .not.methodOfLines )then
                           ! --- TAYLOR TIME-STEPPING --- 
                           stop 111
                         else
                           ! --- METHOD OF LINES (RK) ---
                           ! zero out some forcing terms 
                           do m=0,numPolarizationTerms-1
                              pv(m)=0.
                             fp(m)=0.
                           end do
                             do i3=n3a,n3b
                             do i2=n2a,n2b
                             do i1=n1a,n1b
                               if( mask(i1,i2,i3).gt.0 )then
                             if( addForcing.ne.0 )then  ! do this for now *fix me*
                                 do m=0,5
                                  fv(m)=f(i1,i2,i3,m)
                                 end do
                             end if
                             if( numberOfMaterialRegions.gt.1 )then
                               mr = matMask(i1,i2,i3)
                               if( mr.lt.0 .or. 
     & mr.ge.numberOfMaterialRegions )then  ! do this for now
                                  stop 9999
                               end if
                             end if
                             if( forcingOption.eq.twilightZoneForcing )
     & then
                               if( nd.eq.2 )then
                                 do m=0,numComp-1
                                   ec = m
                                     call ogDeriv(ep, 0,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,ev(m) )
                                     call ogDeriv(ep, 1,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,evt(m) )
                                 end do
                                 ! eval the polarization terms and time derivatives 
                                 do m=0,numPolarizationTerms-1
                                   pc = m+numComp  ! TZ index
                                     call ogDeriv(ep, 0,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),0.,t, pc,pv(m) )
                                     call ogDeriv(ep, 1,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),0.,t, pc,pvt(m) )
                                 end do
                               else
                                 do m=0,numComp-1
                                   ec = m
                                     call ogDeriv(ep, 0,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,ev(m) )
                                     call ogDeriv(ep, 1,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,evt(m) )
                                 end do
                                 ! eval the polarization terms and time derivatives 
                                 do m=0,numPolarizationTerms-1
                                   pc = m+numComp  ! TZ index
                                     call ogDeriv(ep, 0,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, pc,pv(m) )
                                     call ogDeriv(ep, 1,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, pc,pvt(m) )
                                 end do
                              end if
                              !write(*,'(" i1,i2=",2i3," xy=",2(f5.2,1x))') i1,i2,xy(i1,i2,i3,0),xy(i1,i2,i3,1)
                              !write(*,'(" i1,i2=",2i3," ev=",6(f6.3,1x))') i1,i2,(ev(m),m=0,5)
                              !write(*,'(" i1,i2=",2i3," pv=",10(f6.3,1x))') i1,i2,(pv(m),m=0,numPolarizationTerms-1)
                               ! FIX ME -- use opt version here too
                               pc=0  ! P
                               qc=1  ! Q = P.t
                               do k1=1,6
                                 do k2=1,6
                                   ec=k2-1 ! E or H component
                                   do n=1,Np(k1,k2,mr)
                                     a0 = gdmPar(1,n,k1,k2,mr)
                                     a1 = gdmPar(2,n,k1,k2,mr)
                                     b0 = gdmPar(3,n,k1,k2,mr)
                                     b1 = gdmPar(4,n,k1,k2,mr)
                                     ! TZ forcing for polarization equations: 
                                     fp(pc) = pvt(pc) - pv(qc)
                                     fp(qc) = pvt(qc) - (  a0*ev(ec) + 
     & a1*evt(ec) - b0*pv(pc)- b1*pv(qc) )
                                     pc=pc+2
                                     qc=qc+2
                                   end do
                                 end do
                               end do
                            end if
                             ! compute components of the curl(H) and -curl(E)
                              ! --- 3D -----
                               curl(0) =   uy43r(i1,i2,i3,hz)- uz43r(
     & i1,i2,i3,hy)
                               curl(1) =   uz43r(i1,i2,i3,hx)-etax(i1)*
     &  ux43r(i1,i2,i3,hz)
                               curl(2) =  etax(i1)* ux43r(i1,i2,i3,hy)-
     &  uy43r(i1,i2,i3,hx)
                               curl(3) =-( uy43r(i1,i2,i3,ez)- uz43r(
     & i1,i2,i3,ey))
                               curl(4) =-( uz43r(i1,i2,i3,ex)-etax(i1)*
     &  ux43r(i1,i2,i3,ez))
                               curl(5) =-(etax(i1)* ux43r(i1,i2,i3,ey)-
     &  uy43r(i1,i2,i3,ex))
                             if( debug.gt.3 )then
                               write(*,'("----- i1,i2=",2i3)') i1,i2
                             end if
                             ! ---- Compute q = p.t = SUM_k2 SUM_n  p(i1,i2,i3, n,k1,k2, qc )
                             ! opt version 
                               do m=0,5
                                 ptSum(m)=0
                               end do
                             do m=1,numTerms1(mr)
                                ec = ecIndex1(m,mr)
                                qc = qcIndex1(m,mr)
                                ! upiv(m) = p(i1,i2,i3,qc-1)  ! didn't seem to help
                                ! uqiv(m) = p(i1,i2,i3,qc )
                                ! uptSum(ec) = ptSum(ec) + qiv(m) - pv(qc)
                                ptSum(ec) = ptSum(ec) + p(i1,i2,i3,qc) 
     & - pv(qc)
                             end do
                               do m=0,5
                                 curl(m) = curl(m) - ptSum(m)
                               end do
                             ! if( debug.gt.3 )then
                             !   write(*,'("   ptSum=",6(f6.3,1x))') (ptSum(m),m=0,5)
                             ! end if
                               do m=0,5
                                 un(i1,i2,i3,m)= K0i(m,0,mr)*(curl(0)+
     & fv(0)) + K0i(m,1,mr)*(curl(1)+fv(1)) + K0i(m,2,mr)*(curl(2)+fv(
     & 2)) + K0i(m,3,mr)*(curl(3)+fv(3)) + K0i(m,4,mr)*(curl(4)+fv(4))
     &  + K0i(m,5,mr)*(curl(5)+fv(5))
                               end do
                             ! if( .false. )then
                             !   write(*,'(" i1,i2=",2i3," ut=",6(f6.3,1x))') i1,i2,(un(i1,i2,i3,m),m=0,5)
                             ! end if
                             ! --- compute time derivatives of P and Q
                             ! p.t = q
                             ! q.t = a0*E + a1*Et - b0*p - b1*q   
                             ! or 
                             ! q.t = a0*H + a1*Ht - b0*p - b1*q   
                             ! optimized version
                             do m=1,numTerms2(mr)
                               ec = ecIndex2(m,mr)
                               pc = pcIndex2(m,mr)
                               qc = pc+1
                               pct=pc+numComp  ! p.t is stored in un here
                               qct=qc+numComp  ! q.t is stored in un here
                               a0 = a0v(m,mr)
                               a1 = a1v(m,mr)
                               b0 = b0v(m,mr)
                               b1 = b1v(m,mr)
                               un(i1,i2,i3,pct) = p(i1,i2,i3,qc) + fp(
     & pc)
                               un(i1,i2,i3,qct) = a0*u(i1,i2,i3,ec) + 
     & a1*un(i1,i2,i3,ec) - b0*p(i1,i2,i3,pc)- b1*p(i1,i2,i3,qc) + fp(
     & qc)
                               ! uun(i1,i2,i3,pct) = qiv(m) + fp(pc) 
                               ! uun(i1,i2,i3,qct) = a0*uv(ec) + a1*unv(ec) - b0*piv(m)- b1*qiv(m) + fp(qc) 
                             end do
                               end if
                             end do
                             end do
                             end do
                           if( .false. .or. debug > 15 )then
                             stop  4444
                          end if
                         end if
                     else if( useAbsorbingLayer(0).eq.0 .and. 
     & useAbsorbingLayer(1).eq.1 .and. useAbsorbingLayer(2).eq.0 )then
                         if( t.lt.2*dt )then
                           write(*,'("advBA: advance BA GDM dim=3 
     & order=4 grid=rectangular polar=NONE... t=",e10.2)') t
                         end if
                         ! ---- Precompute some indirection arrays to make the dispersion loops go faster ----
                         do mr=0,numberOfMaterialRegions-1
                           m=0
                           pc=0
                           qc=1
                           do k1=1,6
                             ec=k1-1 ! E or H component
                             do k2=1,6
                               do n=1,Np(k1,k2,mr)
                                 m=m+1
                                 ecIndex1(m,mr)=ec
                                 qcIndex1(m,mr)=qc
                                 ! ptSum(ec) = ptSum(ec) + p(i1,i2,i3,qc) - pv(qc)
                                 qc=qc+2
                               end do
                             end do
                             ! subtract off P.t = sum Q_m
                             ! curl(ec) = curl(ec) - ptSum(ec)
                           end do
                           numTerms1(mr)=m
                           if( numTerms1(mr)
     & .gt.maxNumPolarizationTerms )then
                             write(*,'(" ERROR numTerms1=",i6," too 
     & big")') numTerms1(mr)
                             stop 1616
                           end if
                           pc=0  ! P
                           qc=1  ! Q = P.t
                           m=0
                           do k1=1,6
                             do k2=1,6
                               ec=k2-1 ! This GDM term involves this E or H component
                               do n=1,Np(k1,k2,mr)
                                 a0 = gdmPar(1,n,k1,k2,mr)
                                 a1 = gdmPar(2,n,k1,k2,mr)
                                 b0 = gdmPar(3,n,k1,k2,mr)
                                 b1 = gdmPar(4,n,k1,k2,mr)
                                 pct=pc+numComp  ! p.t is stored in un here
                                 qct=qc+numComp  ! q.t is stored in un here
                                 m=m+1
                                 ecIndex2(m,mr)=ec
                                 pcIndex2(m,mr)=pc
                                 a0v(m,mr)=a0
                                 a1v(m,mr)=a1
                                 b0v(m,mr)=b0
                                 b1v(m,mr)=b1
                                 ! un(i1,i2,i3,pct) = p(i1,i2,i3,qc) + fp(pc) 
                                 ! un(i1,i2,i3,qct) = a0*u(i1,i2,i3,ec) + a1*un(i1,i2,i3,ec) - b0*p(i1,i2,i3,pc)- b1*p(i1,i2,i3,qc) + fp(qc) 
                                 pc=pc+2
                                 qc=qc+2
                               end do
                             end do
                           end do
                           numTerms2(mr)=m
                           if( numTerms2(mr)
     & .gt.maxNumPolarizationTerms )then
                             write(*,'(" ERROR numTerms2=",i6," too 
     & big")') numTerms2(mr)
                             stop 1616
                           end if
                         end do ! end do mr
                         mr=0
                         if( .not.methodOfLines )then
                           ! --- TAYLOR TIME-STEPPING --- 
                           stop 111
                         else
                           ! --- METHOD OF LINES (RK) ---
                           ! zero out some forcing terms 
                           do m=0,numPolarizationTerms-1
                              pv(m)=0.
                             fp(m)=0.
                           end do
                             do i3=n3a,n3b
                             do i2=n2a,n2b
                             do i1=n1a,n1b
                               if( mask(i1,i2,i3).gt.0 )then
                             if( addForcing.ne.0 )then  ! do this for now *fix me*
                                 do m=0,5
                                  fv(m)=f(i1,i2,i3,m)
                                 end do
                             end if
                             if( numberOfMaterialRegions.gt.1 )then
                               mr = matMask(i1,i2,i3)
                               if( mr.lt.0 .or. 
     & mr.ge.numberOfMaterialRegions )then  ! do this for now
                                  stop 9999
                               end if
                             end if
                             if( forcingOption.eq.twilightZoneForcing )
     & then
                               if( nd.eq.2 )then
                                 do m=0,numComp-1
                                   ec = m
                                     call ogDeriv(ep, 0,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,ev(m) )
                                     call ogDeriv(ep, 1,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,evt(m) )
                                 end do
                                 ! eval the polarization terms and time derivatives 
                                 do m=0,numPolarizationTerms-1
                                   pc = m+numComp  ! TZ index
                                     call ogDeriv(ep, 0,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),0.,t, pc,pv(m) )
                                     call ogDeriv(ep, 1,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),0.,t, pc,pvt(m) )
                                 end do
                               else
                                 do m=0,numComp-1
                                   ec = m
                                     call ogDeriv(ep, 0,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,ev(m) )
                                     call ogDeriv(ep, 1,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,evt(m) )
                                 end do
                                 ! eval the polarization terms and time derivatives 
                                 do m=0,numPolarizationTerms-1
                                   pc = m+numComp  ! TZ index
                                     call ogDeriv(ep, 0,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, pc,pv(m) )
                                     call ogDeriv(ep, 1,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, pc,pvt(m) )
                                 end do
                              end if
                              !write(*,'(" i1,i2=",2i3," xy=",2(f5.2,1x))') i1,i2,xy(i1,i2,i3,0),xy(i1,i2,i3,1)
                              !write(*,'(" i1,i2=",2i3," ev=",6(f6.3,1x))') i1,i2,(ev(m),m=0,5)
                              !write(*,'(" i1,i2=",2i3," pv=",10(f6.3,1x))') i1,i2,(pv(m),m=0,numPolarizationTerms-1)
                               ! FIX ME -- use opt version here too
                               pc=0  ! P
                               qc=1  ! Q = P.t
                               do k1=1,6
                                 do k2=1,6
                                   ec=k2-1 ! E or H component
                                   do n=1,Np(k1,k2,mr)
                                     a0 = gdmPar(1,n,k1,k2,mr)
                                     a1 = gdmPar(2,n,k1,k2,mr)
                                     b0 = gdmPar(3,n,k1,k2,mr)
                                     b1 = gdmPar(4,n,k1,k2,mr)
                                     ! TZ forcing for polarization equations: 
                                     fp(pc) = pvt(pc) - pv(qc)
                                     fp(qc) = pvt(qc) - (  a0*ev(ec) + 
     & a1*evt(ec) - b0*pv(pc)- b1*pv(qc) )
                                     pc=pc+2
                                     qc=qc+2
                                   end do
                                 end do
                               end do
                            end if
                             ! compute components of the curl(H) and -curl(E)
                              ! --- 3D -----
                               curl(0) =  etay(i2)* uy43r(i1,i2,i3,hz)-
     &  uz43r(i1,i2,i3,hy)
                               curl(1) =   uz43r(i1,i2,i3,hx)- ux43r(
     & i1,i2,i3,hz)
                               curl(2) =   ux43r(i1,i2,i3,hy)-etay(i2)*
     &  uy43r(i1,i2,i3,hx)
                               curl(3) =-(etay(i2)* uy43r(i1,i2,i3,ez)-
     &  uz43r(i1,i2,i3,ey))
                               curl(4) =-( uz43r(i1,i2,i3,ex)- ux43r(
     & i1,i2,i3,ez))
                               curl(5) =-( ux43r(i1,i2,i3,ey)-etay(i2)*
     &  uy43r(i1,i2,i3,ex))
                             if( debug.gt.3 )then
                               write(*,'("----- i1,i2=",2i3)') i1,i2
                             end if
                             ! ---- Compute q = p.t = SUM_k2 SUM_n  p(i1,i2,i3, n,k1,k2, qc )
                             ! opt version 
                               do m=0,5
                                 ptSum(m)=0
                               end do
                             do m=1,numTerms1(mr)
                                ec = ecIndex1(m,mr)
                                qc = qcIndex1(m,mr)
                                ! upiv(m) = p(i1,i2,i3,qc-1)  ! didn't seem to help
                                ! uqiv(m) = p(i1,i2,i3,qc )
                                ! uptSum(ec) = ptSum(ec) + qiv(m) - pv(qc)
                                ptSum(ec) = ptSum(ec) + p(i1,i2,i3,qc) 
     & - pv(qc)
                             end do
                               do m=0,5
                                 curl(m) = curl(m) - ptSum(m)
                               end do
                             ! if( debug.gt.3 )then
                             !   write(*,'("   ptSum=",6(f6.3,1x))') (ptSum(m),m=0,5)
                             ! end if
                               do m=0,5
                                 un(i1,i2,i3,m)= K0i(m,0,mr)*(curl(0)+
     & fv(0)) + K0i(m,1,mr)*(curl(1)+fv(1)) + K0i(m,2,mr)*(curl(2)+fv(
     & 2)) + K0i(m,3,mr)*(curl(3)+fv(3)) + K0i(m,4,mr)*(curl(4)+fv(4))
     &  + K0i(m,5,mr)*(curl(5)+fv(5))
                               end do
                             ! if( .false. )then
                             !   write(*,'(" i1,i2=",2i3," ut=",6(f6.3,1x))') i1,i2,(un(i1,i2,i3,m),m=0,5)
                             ! end if
                             ! --- compute time derivatives of P and Q
                             ! p.t = q
                             ! q.t = a0*E + a1*Et - b0*p - b1*q   
                             ! or 
                             ! q.t = a0*H + a1*Ht - b0*p - b1*q   
                             ! optimized version
                             do m=1,numTerms2(mr)
                               ec = ecIndex2(m,mr)
                               pc = pcIndex2(m,mr)
                               qc = pc+1
                               pct=pc+numComp  ! p.t is stored in un here
                               qct=qc+numComp  ! q.t is stored in un here
                               a0 = a0v(m,mr)
                               a1 = a1v(m,mr)
                               b0 = b0v(m,mr)
                               b1 = b1v(m,mr)
                               un(i1,i2,i3,pct) = p(i1,i2,i3,qc) + fp(
     & pc)
                               un(i1,i2,i3,qct) = a0*u(i1,i2,i3,ec) + 
     & a1*un(i1,i2,i3,ec) - b0*p(i1,i2,i3,pc)- b1*p(i1,i2,i3,qc) + fp(
     & qc)
                               ! uun(i1,i2,i3,pct) = qiv(m) + fp(pc) 
                               ! uun(i1,i2,i3,qct) = a0*uv(ec) + a1*unv(ec) - b0*piv(m)- b1*qiv(m) + fp(qc) 
                             end do
                               end if
                             end do
                             end do
                             end do
                           if( .false. .or. debug > 15 )then
                             stop  4444
                          end if
                         end if
                     else if( useAbsorbingLayer(0).eq.0 .and. 
     & useAbsorbingLayer(1).eq.0 .and. useAbsorbingLayer(2).eq.1 )then
                         if( t.lt.2*dt )then
                           write(*,'("advBA: advance BA GDM dim=3 
     & order=4 grid=rectangular polar=NONE... t=",e10.2)') t
                         end if
                         ! ---- Precompute some indirection arrays to make the dispersion loops go faster ----
                         do mr=0,numberOfMaterialRegions-1
                           m=0
                           pc=0
                           qc=1
                           do k1=1,6
                             ec=k1-1 ! E or H component
                             do k2=1,6
                               do n=1,Np(k1,k2,mr)
                                 m=m+1
                                 ecIndex1(m,mr)=ec
                                 qcIndex1(m,mr)=qc
                                 ! ptSum(ec) = ptSum(ec) + p(i1,i2,i3,qc) - pv(qc)
                                 qc=qc+2
                               end do
                             end do
                             ! subtract off P.t = sum Q_m
                             ! curl(ec) = curl(ec) - ptSum(ec)
                           end do
                           numTerms1(mr)=m
                           if( numTerms1(mr)
     & .gt.maxNumPolarizationTerms )then
                             write(*,'(" ERROR numTerms1=",i6," too 
     & big")') numTerms1(mr)
                             stop 1616
                           end if
                           pc=0  ! P
                           qc=1  ! Q = P.t
                           m=0
                           do k1=1,6
                             do k2=1,6
                               ec=k2-1 ! This GDM term involves this E or H component
                               do n=1,Np(k1,k2,mr)
                                 a0 = gdmPar(1,n,k1,k2,mr)
                                 a1 = gdmPar(2,n,k1,k2,mr)
                                 b0 = gdmPar(3,n,k1,k2,mr)
                                 b1 = gdmPar(4,n,k1,k2,mr)
                                 pct=pc+numComp  ! p.t is stored in un here
                                 qct=qc+numComp  ! q.t is stored in un here
                                 m=m+1
                                 ecIndex2(m,mr)=ec
                                 pcIndex2(m,mr)=pc
                                 a0v(m,mr)=a0
                                 a1v(m,mr)=a1
                                 b0v(m,mr)=b0
                                 b1v(m,mr)=b1
                                 ! un(i1,i2,i3,pct) = p(i1,i2,i3,qc) + fp(pc) 
                                 ! un(i1,i2,i3,qct) = a0*u(i1,i2,i3,ec) + a1*un(i1,i2,i3,ec) - b0*p(i1,i2,i3,pc)- b1*p(i1,i2,i3,qc) + fp(qc) 
                                 pc=pc+2
                                 qc=qc+2
                               end do
                             end do
                           end do
                           numTerms2(mr)=m
                           if( numTerms2(mr)
     & .gt.maxNumPolarizationTerms )then
                             write(*,'(" ERROR numTerms2=",i6," too 
     & big")') numTerms2(mr)
                             stop 1616
                           end if
                         end do ! end do mr
                         mr=0
                         if( .not.methodOfLines )then
                           ! --- TAYLOR TIME-STEPPING --- 
                           stop 111
                         else
                           ! --- METHOD OF LINES (RK) ---
                           ! zero out some forcing terms 
                           do m=0,numPolarizationTerms-1
                              pv(m)=0.
                             fp(m)=0.
                           end do
                             do i3=n3a,n3b
                             do i2=n2a,n2b
                             do i1=n1a,n1b
                               if( mask(i1,i2,i3).gt.0 )then
                             if( addForcing.ne.0 )then  ! do this for now *fix me*
                                 do m=0,5
                                  fv(m)=f(i1,i2,i3,m)
                                 end do
                             end if
                             if( numberOfMaterialRegions.gt.1 )then
                               mr = matMask(i1,i2,i3)
                               if( mr.lt.0 .or. 
     & mr.ge.numberOfMaterialRegions )then  ! do this for now
                                  stop 9999
                               end if
                             end if
                             if( forcingOption.eq.twilightZoneForcing )
     & then
                               if( nd.eq.2 )then
                                 do m=0,numComp-1
                                   ec = m
                                     call ogDeriv(ep, 0,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,ev(m) )
                                     call ogDeriv(ep, 1,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,evt(m) )
                                 end do
                                 ! eval the polarization terms and time derivatives 
                                 do m=0,numPolarizationTerms-1
                                   pc = m+numComp  ! TZ index
                                     call ogDeriv(ep, 0,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),0.,t, pc,pv(m) )
                                     call ogDeriv(ep, 1,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),0.,t, pc,pvt(m) )
                                 end do
                               else
                                 do m=0,numComp-1
                                   ec = m
                                     call ogDeriv(ep, 0,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,ev(m) )
                                     call ogDeriv(ep, 1,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,evt(m) )
                                 end do
                                 ! eval the polarization terms and time derivatives 
                                 do m=0,numPolarizationTerms-1
                                   pc = m+numComp  ! TZ index
                                     call ogDeriv(ep, 0,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, pc,pv(m) )
                                     call ogDeriv(ep, 1,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, pc,pvt(m) )
                                 end do
                              end if
                              !write(*,'(" i1,i2=",2i3," xy=",2(f5.2,1x))') i1,i2,xy(i1,i2,i3,0),xy(i1,i2,i3,1)
                              !write(*,'(" i1,i2=",2i3," ev=",6(f6.3,1x))') i1,i2,(ev(m),m=0,5)
                              !write(*,'(" i1,i2=",2i3," pv=",10(f6.3,1x))') i1,i2,(pv(m),m=0,numPolarizationTerms-1)
                               ! FIX ME -- use opt version here too
                               pc=0  ! P
                               qc=1  ! Q = P.t
                               do k1=1,6
                                 do k2=1,6
                                   ec=k2-1 ! E or H component
                                   do n=1,Np(k1,k2,mr)
                                     a0 = gdmPar(1,n,k1,k2,mr)
                                     a1 = gdmPar(2,n,k1,k2,mr)
                                     b0 = gdmPar(3,n,k1,k2,mr)
                                     b1 = gdmPar(4,n,k1,k2,mr)
                                     ! TZ forcing for polarization equations: 
                                     fp(pc) = pvt(pc) - pv(qc)
                                     fp(qc) = pvt(qc) - (  a0*ev(ec) + 
     & a1*evt(ec) - b0*pv(pc)- b1*pv(qc) )
                                     pc=pc+2
                                     qc=qc+2
                                   end do
                                 end do
                               end do
                            end if
                             ! compute components of the curl(H) and -curl(E)
                              ! --- 3D -----
                               curl(0) =   uy43r(i1,i2,i3,hz)-etaz(i3)*
     &  uz43r(i1,i2,i3,hy)
                               curl(1) =  etaz(i3)* uz43r(i1,i2,i3,hx)-
     &  ux43r(i1,i2,i3,hz)
                               curl(2) =   ux43r(i1,i2,i3,hy)- uy43r(
     & i1,i2,i3,hx)
                               curl(3) =-( uy43r(i1,i2,i3,ez)-etaz(i3)*
     &  uz43r(i1,i2,i3,ey))
                               curl(4) =-(etaz(i3)* uz43r(i1,i2,i3,ex)-
     &  ux43r(i1,i2,i3,ez))
                               curl(5) =-( ux43r(i1,i2,i3,ey)- uy43r(
     & i1,i2,i3,ex))
                             if( debug.gt.3 )then
                               write(*,'("----- i1,i2=",2i3)') i1,i2
                             end if
                             ! ---- Compute q = p.t = SUM_k2 SUM_n  p(i1,i2,i3, n,k1,k2, qc )
                             ! opt version 
                               do m=0,5
                                 ptSum(m)=0
                               end do
                             do m=1,numTerms1(mr)
                                ec = ecIndex1(m,mr)
                                qc = qcIndex1(m,mr)
                                ! upiv(m) = p(i1,i2,i3,qc-1)  ! didn't seem to help
                                ! uqiv(m) = p(i1,i2,i3,qc )
                                ! uptSum(ec) = ptSum(ec) + qiv(m) - pv(qc)
                                ptSum(ec) = ptSum(ec) + p(i1,i2,i3,qc) 
     & - pv(qc)
                             end do
                               do m=0,5
                                 curl(m) = curl(m) - ptSum(m)
                               end do
                             ! if( debug.gt.3 )then
                             !   write(*,'("   ptSum=",6(f6.3,1x))') (ptSum(m),m=0,5)
                             ! end if
                               do m=0,5
                                 un(i1,i2,i3,m)= K0i(m,0,mr)*(curl(0)+
     & fv(0)) + K0i(m,1,mr)*(curl(1)+fv(1)) + K0i(m,2,mr)*(curl(2)+fv(
     & 2)) + K0i(m,3,mr)*(curl(3)+fv(3)) + K0i(m,4,mr)*(curl(4)+fv(4))
     &  + K0i(m,5,mr)*(curl(5)+fv(5))
                               end do
                             ! if( .false. )then
                             !   write(*,'(" i1,i2=",2i3," ut=",6(f6.3,1x))') i1,i2,(un(i1,i2,i3,m),m=0,5)
                             ! end if
                             ! --- compute time derivatives of P and Q
                             ! p.t = q
                             ! q.t = a0*E + a1*Et - b0*p - b1*q   
                             ! or 
                             ! q.t = a0*H + a1*Ht - b0*p - b1*q   
                             ! optimized version
                             do m=1,numTerms2(mr)
                               ec = ecIndex2(m,mr)
                               pc = pcIndex2(m,mr)
                               qc = pc+1
                               pct=pc+numComp  ! p.t is stored in un here
                               qct=qc+numComp  ! q.t is stored in un here
                               a0 = a0v(m,mr)
                               a1 = a1v(m,mr)
                               b0 = b0v(m,mr)
                               b1 = b1v(m,mr)
                               un(i1,i2,i3,pct) = p(i1,i2,i3,qc) + fp(
     & pc)
                               un(i1,i2,i3,qct) = a0*u(i1,i2,i3,ec) + 
     & a1*un(i1,i2,i3,ec) - b0*p(i1,i2,i3,pc)- b1*p(i1,i2,i3,qc) + fp(
     & qc)
                               ! uun(i1,i2,i3,pct) = qiv(m) + fp(pc) 
                               ! uun(i1,i2,i3,qct) = a0*uv(ec) + a1*unv(ec) - b0*piv(m)- b1*qiv(m) + fp(qc) 
                             end do
                               end if
                             end do
                             end do
                             end do
                           if( .false. .or. debug > 15 )then
                             stop  4444
                          end if
                         end if
                     else
                         if( t.lt.2*dt )then
                           write(*,'("advBA: advance BA GDM dim=3 
     & order=4 grid=rectangular polar=NONE... t=",e10.2)') t
                         end if
                         ! ---- Precompute some indirection arrays to make the dispersion loops go faster ----
                         do mr=0,numberOfMaterialRegions-1
                           m=0
                           pc=0
                           qc=1
                           do k1=1,6
                             ec=k1-1 ! E or H component
                             do k2=1,6
                               do n=1,Np(k1,k2,mr)
                                 m=m+1
                                 ecIndex1(m,mr)=ec
                                 qcIndex1(m,mr)=qc
                                 ! ptSum(ec) = ptSum(ec) + p(i1,i2,i3,qc) - pv(qc)
                                 qc=qc+2
                               end do
                             end do
                             ! subtract off P.t = sum Q_m
                             ! curl(ec) = curl(ec) - ptSum(ec)
                           end do
                           numTerms1(mr)=m
                           if( numTerms1(mr)
     & .gt.maxNumPolarizationTerms )then
                             write(*,'(" ERROR numTerms1=",i6," too 
     & big")') numTerms1(mr)
                             stop 1616
                           end if
                           pc=0  ! P
                           qc=1  ! Q = P.t
                           m=0
                           do k1=1,6
                             do k2=1,6
                               ec=k2-1 ! This GDM term involves this E or H component
                               do n=1,Np(k1,k2,mr)
                                 a0 = gdmPar(1,n,k1,k2,mr)
                                 a1 = gdmPar(2,n,k1,k2,mr)
                                 b0 = gdmPar(3,n,k1,k2,mr)
                                 b1 = gdmPar(4,n,k1,k2,mr)
                                 pct=pc+numComp  ! p.t is stored in un here
                                 qct=qc+numComp  ! q.t is stored in un here
                                 m=m+1
                                 ecIndex2(m,mr)=ec
                                 pcIndex2(m,mr)=pc
                                 a0v(m,mr)=a0
                                 a1v(m,mr)=a1
                                 b0v(m,mr)=b0
                                 b1v(m,mr)=b1
                                 ! un(i1,i2,i3,pct) = p(i1,i2,i3,qc) + fp(pc) 
                                 ! un(i1,i2,i3,qct) = a0*u(i1,i2,i3,ec) + a1*un(i1,i2,i3,ec) - b0*p(i1,i2,i3,pc)- b1*p(i1,i2,i3,qc) + fp(qc) 
                                 pc=pc+2
                                 qc=qc+2
                               end do
                             end do
                           end do
                           numTerms2(mr)=m
                           if( numTerms2(mr)
     & .gt.maxNumPolarizationTerms )then
                             write(*,'(" ERROR numTerms2=",i6," too 
     & big")') numTerms2(mr)
                             stop 1616
                           end if
                         end do ! end do mr
                         mr=0
                         if( .not.methodOfLines )then
                           ! --- TAYLOR TIME-STEPPING --- 
                           stop 111
                         else
                           ! --- METHOD OF LINES (RK) ---
                           ! zero out some forcing terms 
                           do m=0,numPolarizationTerms-1
                              pv(m)=0.
                             fp(m)=0.
                           end do
                             do i3=n3a,n3b
                             do i2=n2a,n2b
                             do i1=n1a,n1b
                               if( mask(i1,i2,i3).gt.0 )then
                             if( addForcing.ne.0 )then  ! do this for now *fix me*
                                 do m=0,5
                                  fv(m)=f(i1,i2,i3,m)
                                 end do
                             end if
                             if( numberOfMaterialRegions.gt.1 )then
                               mr = matMask(i1,i2,i3)
                               if( mr.lt.0 .or. 
     & mr.ge.numberOfMaterialRegions )then  ! do this for now
                                  stop 9999
                               end if
                             end if
                             if( forcingOption.eq.twilightZoneForcing )
     & then
                               if( nd.eq.2 )then
                                 do m=0,numComp-1
                                   ec = m
                                     call ogDeriv(ep, 0,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,ev(m) )
                                     call ogDeriv(ep, 1,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,evt(m) )
                                 end do
                                 ! eval the polarization terms and time derivatives 
                                 do m=0,numPolarizationTerms-1
                                   pc = m+numComp  ! TZ index
                                     call ogDeriv(ep, 0,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),0.,t, pc,pv(m) )
                                     call ogDeriv(ep, 1,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),0.,t, pc,pvt(m) )
                                 end do
                               else
                                 do m=0,numComp-1
                                   ec = m
                                     call ogDeriv(ep, 0,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,ev(m) )
                                     call ogDeriv(ep, 1,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,evt(m) )
                                 end do
                                 ! eval the polarization terms and time derivatives 
                                 do m=0,numPolarizationTerms-1
                                   pc = m+numComp  ! TZ index
                                     call ogDeriv(ep, 0,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, pc,pv(m) )
                                     call ogDeriv(ep, 1,0,0,0, xy(i1,
     & i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, pc,pvt(m) )
                                 end do
                              end if
                              !write(*,'(" i1,i2=",2i3," xy=",2(f5.2,1x))') i1,i2,xy(i1,i2,i3,0),xy(i1,i2,i3,1)
                              !write(*,'(" i1,i2=",2i3," ev=",6(f6.3,1x))') i1,i2,(ev(m),m=0,5)
                              !write(*,'(" i1,i2=",2i3," pv=",10(f6.3,1x))') i1,i2,(pv(m),m=0,numPolarizationTerms-1)
                               ! FIX ME -- use opt version here too
                               pc=0  ! P
                               qc=1  ! Q = P.t
                               do k1=1,6
                                 do k2=1,6
                                   ec=k2-1 ! E or H component
                                   do n=1,Np(k1,k2,mr)
                                     a0 = gdmPar(1,n,k1,k2,mr)
                                     a1 = gdmPar(2,n,k1,k2,mr)
                                     b0 = gdmPar(3,n,k1,k2,mr)
                                     b1 = gdmPar(4,n,k1,k2,mr)
                                     ! TZ forcing for polarization equations: 
                                     fp(pc) = pvt(pc) - pv(qc)
                                     fp(qc) = pvt(qc) - (  a0*ev(ec) + 
     & a1*evt(ec) - b0*pv(pc)- b1*pv(qc) )
                                     pc=pc+2
                                     qc=qc+2
                                   end do
                                 end do
                               end do
                            end if
                             ! compute components of the curl(H) and -curl(E)
                              ! --- 3D -----
                               curl(0) =   uy43r(i1,i2,i3,hz)- uz43r(
     & i1,i2,i3,hy)
                               curl(1) =   uz43r(i1,i2,i3,hx)- ux43r(
     & i1,i2,i3,hz)
                               curl(2) =   ux43r(i1,i2,i3,hy)- uy43r(
     & i1,i2,i3,hx)
                               curl(3) =-( uy43r(i1,i2,i3,ez)- uz43r(
     & i1,i2,i3,ey))
                               curl(4) =-( uz43r(i1,i2,i3,ex)- ux43r(
     & i1,i2,i3,ez))
                               curl(5) =-( ux43r(i1,i2,i3,ey)- uy43r(
     & i1,i2,i3,ex))
                             if( debug.gt.3 )then
                               write(*,'("----- i1,i2=",2i3)') i1,i2
                             end if
                             ! ---- Compute q = p.t = SUM_k2 SUM_n  p(i1,i2,i3, n,k1,k2, qc )
                             ! opt version 
                               do m=0,5
                                 ptSum(m)=0
                               end do
                             do m=1,numTerms1(mr)
                                ec = ecIndex1(m,mr)
                                qc = qcIndex1(m,mr)
                                ! upiv(m) = p(i1,i2,i3,qc-1)  ! didn't seem to help
                                ! uqiv(m) = p(i1,i2,i3,qc )
                                ! uptSum(ec) = ptSum(ec) + qiv(m) - pv(qc)
                                ptSum(ec) = ptSum(ec) + p(i1,i2,i3,qc) 
     & - pv(qc)
                             end do
                               do m=0,5
                                 curl(m) = curl(m) - ptSum(m)
                               end do
                             ! if( debug.gt.3 )then
                             !   write(*,'("   ptSum=",6(f6.3,1x))') (ptSum(m),m=0,5)
                             ! end if
                               do m=0,5
                                 un(i1,i2,i3,m)= K0i(m,0,mr)*(curl(0)+
     & fv(0)) + K0i(m,1,mr)*(curl(1)+fv(1)) + K0i(m,2,mr)*(curl(2)+fv(
     & 2)) + K0i(m,3,mr)*(curl(3)+fv(3)) + K0i(m,4,mr)*(curl(4)+fv(4))
     &  + K0i(m,5,mr)*(curl(5)+fv(5))
                               end do
                             ! if( .false. )then
                             !   write(*,'(" i1,i2=",2i3," ut=",6(f6.3,1x))') i1,i2,(un(i1,i2,i3,m),m=0,5)
                             ! end if
                             ! --- compute time derivatives of P and Q
                             ! p.t = q
                             ! q.t = a0*E + a1*Et - b0*p - b1*q   
                             ! or 
                             ! q.t = a0*H + a1*Ht - b0*p - b1*q   
                             ! optimized version
                             do m=1,numTerms2(mr)
                               ec = ecIndex2(m,mr)
                               pc = pcIndex2(m,mr)
                               qc = pc+1
                               pct=pc+numComp  ! p.t is stored in un here
                               qct=qc+numComp  ! q.t is stored in un here
                               a0 = a0v(m,mr)
                               a1 = a1v(m,mr)
                               b0 = b0v(m,mr)
                               b1 = b1v(m,mr)
                               un(i1,i2,i3,pct) = p(i1,i2,i3,qc) + fp(
     & pc)
                               un(i1,i2,i3,qct) = a0*u(i1,i2,i3,ec) + 
     & a1*un(i1,i2,i3,ec) - b0*p(i1,i2,i3,pc)- b1*p(i1,i2,i3,qc) + fp(
     & qc)
                               ! uun(i1,i2,i3,pct) = qiv(m) + fp(pc) 
                               ! uun(i1,i2,i3,qct) = a0*uv(ec) + a1*unv(ec) - b0*piv(m)- b1*qiv(m) + fp(qc) 
                             end do
                               end if
                             end do
                             end do
                             end do
                           if( .false. .or. debug > 15 )then
                             stop  4444
                          end if
                         end if
                     end if
               end if
             end if
       end if
       if( gridType.eq.rectangular )then
      !       **********************************************
      !       *************** rectangular ******************
      !       **********************************************
          ! ------------------------------------------------------------------------------
          !    3D : 4th order  (rectangular)
          ! ------------------------------------------------------------------------------
          if( addDissipation )then
            if( t.le.3*dt )then
              write(*,'("advBA: order=4: addDissipation=",l2," adc=",
     & e10.2)') addDissipation,adc
            end if
            sigma1=adc
            adxSosup(0)=sigma1/(2*64.)
            adxSosup(1)=sigma1/(2*64.)
            adxSosup(2)=sigma1/(2*64.)
              do i3=n3a,n3b
              do i2=n2a,n2b
              do i1=n1a,n1b
                if( mask(i1,i2,i3).gt.0 )then
              do m=ex,hz
                un(i1,i2,i3,m)=u(i1,i2,i3,m) + (-20.*u(i1,i2,i3,m)+15.*
     & (u(i1+1,i2,i3,m)+u(i1-1,i2,i3,m))-6.*(u(i1+2,i2,i3,m)+u(i1-2,
     & i2,i3,m))+(u(i1+3,i2,i3,m)+u(i1-3,i2,i3,m)))*adxSosup(0)+(-20.*
     & u(i1,i2,i3,m)+15.*(u(i1,i2+1,i3,m)+u(i1,i2-1,i3,m))-6.*(u(i1,
     & i2+2,i3,m)+u(i1,i2-2,i3,m))+(u(i1,i2+3,i3,m)+u(i1,i2-3,i3,m)))*
     & adxSosup(1)+(-20.*u(i1,i2,i3,m)+15.*(u(i1,i2,i3+1,m)+u(i1,i2,
     & i3-1,m))-6.*(u(i1,i2,i3+2,m)+u(i1,i2,i3-2,m))+(u(i1,i2,i3+3,m)+
     & u(i1,i2,i3-3,m)))*adxSosup(2)
              end do
                end if
              end do
              end do
              end do
          end if
      !    else if( materialType.eq.isotropic )then
      !      ! updateMx3dOrder4()
      !      stop 7676
      !    else if( dispersionModel.eq.noDispersion )then
      !      updateBA3dOrder4()
      !    else
      !      if( useOpt )then
      !        if( useSuperGrid.eq.0  ) then
      !          updateBAGDMOpt(3,4,rectangular,,,)
      !        else
      !          updateBAGDMOpt(3,4,rectangular,etax(i1)*,etay(i2)*,etaz(i3)*)
      !        end if 
      !      else 
      !        updateBAGDM(3,4,rectangular)
      !      end if 
      !    end if 
        ! End if rectnagular 
       else
       end if
       return
       end
