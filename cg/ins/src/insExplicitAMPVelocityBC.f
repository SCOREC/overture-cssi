! This file automatically generated from insExplicitAMPVelocityBC.bf with bpp.
!=============================================================================================================
!
!   Routines to assign AMP velocity Boundary conditions for explicit time-stepping
!
! Notes:
!   May 3, 2018 -- initial version
!============================================================================================================
!

c Use this next macro to declare the statement functions that are defined below
c To include derivatives of rx use OPTION=RX


c Define statement functions for difference approximations of order 2 
c To include derivatives of rx use OPTION=RX
c To include derivatives of rx use OPTION=RX






! ************************************************************************************************
!  This macro is used for looping over the faces of a grid to assign booundary conditions
!
! extra: extra points to assign
!          Case 1: extra=numberOfGhostPoints -- for assigning extended boundaries
!          Case 2: extra=-1 -- for assigning ghost points but not including extended boundaries
! numberOfGhostPoints : number of ghost points (1 for 2nd order, 2 for fourth-order ...)
!
!
! Output:
!  n1a,n1b,n2a,n2b,n3a,n3b : from gridIndexRange
!  nn1a,nn1b,nn2a,nn2b,nn3a,nn3b : includes "extra" points
! 
! ***********************************************************************************************










! ==========================================================================================================
! Apply the TRACTION FREE boundary condition to determine the velocity on the ghost points
!   RECTANGULAR GRID CASE
! 
! ORDER: 2 or 4
! DIR = x, y, z
! DIM = 2 or 3 dimensions
!
! ==========================================================================================================


! ==========================================================================================================
! Apply the TRACTION boundary condition to determine the velocity on the ghost points
!  Curvilinear grid case
! 
!   tv.sigma.nv = tv.solidTraction
!    div( v ) = 0
!  sigma_{ij} = mu*( (v_i)_j + (v_j)_i )
!
! ORDER: 2 or 4
! DIR = r,s,t
! GRIDTYPE: rectangular, curvilinear
!
! ==========================================================================================================
! ====================== END TRACTION 2D CURVILINEAR =========================


! ==========================================================================================================
! Apply the TRACTION boundary condition to determine the velocity on the ghost points
!  Curvilinear grid case **THREE-DIMENSIONS***
! 
! ORDER: 2 or 4
! DIR = r,s,t
! GRIDTYPE: rectangular, curvilinear
!
! ==========================================================================================================
! ====================== END TRACTION 3D Order=2 CURVILINEAR =========================




! ==========================================================================================================
!   ** NEW ** June 26, 2018
!
! Apply the AMP VELOCITY boundary condition to determine the velocity on the ghost points
!  Curvilinear grid case
! 
! ORDER: 2 or 4
! DIR = r,s,t
! GRIDTYPE: rectangular, curvilinear
!
! ==========================================================================================================
! ====================== END AMP 2D CURVILINEAR =========================






! ==========================================================================================================
!   ** OLD**
! Apply the AMP VELOCITY boundary condition to determine the velocity on the ghost points
!  Curvilinear grid case
! 
! ORDER: 2 or 4
! DIR = r,s,t
! GRIDTYPE: rectangular, curvilinear
!
! ==========================================================================================================
! ====================== END *OLD* AMP 2D CURVILINEAR =========================




      subroutine insExplicitAMPVelocityBC(bcOption, nd,nd1a,nd1b,nd2a,
     & nd2b,nd3a,nd3b,nd4a,nd4b,
     & ipar,rpar, u,un, mask, x,rsxy, nds1a,nds1b,nds2a,nds2b,nds3a,
     & nds3b, solidVelocity,
     & solidTraction, fluidVelocity, bc, gridIndexRange, ierr )
!=============================================================================================================
!     Assign traction (free surface) Boundary conditions
!
! u : fill in interface conditions to u, on input contains solution on boundary and interior at time t
!     and contains a predicted value for p at all values
! un : solution at time t-dt
!
!============================================================================================================
      implicit none
      integer nd, nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b
      real ep ! holds pointer to OGFunction
      real u(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
      real un(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
      real x(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1)
      real rsxy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1,0:nd-1)

      integer nds1a,nds1b,nds2a,nds2b,nds3a,nds3b,nds4a,nds4b
      real solidVelocity(nds1a:nds1b,nds2a:nds2b,nds3a:nds3b,0:nd-1)
      real solidTraction(nds1a:nds1b,nds2a:nds2b,nds3a:nds3b,0:nd-1)
      real fluidVelocity(nds1a:nds1b,nds2a:nds2b,nds3a:nds3b,0:nd-1)
      real rpar(0:*),dt

      integer mask(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
      integer gridIndexRange(0:1,0:2), bc(0:1,0:2)
      integer ipar(0:*),ierr

!.......local
      ! temp work space -- *fix me*
      real uNew(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)

      integer numberOfProcessors,debug,myid
      integer kd,kd3,i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b,c,nr0,nr1
      integer is,j1,j2,j3,side,axis,twilightZone,bcOption
      integer i1m,i2m,i3m,i1p,i2p,i3p
      integer pc,uc,vc,wc,sc,grid,orderOfAccuracy,gridIsMoving,
     & useWhereMask,tc,assignTemperature
      integer gridType,gridIsImplicit,implicitMethod,implicitOption,
     & isAxisymmetric
      integer use2ndOrderAD,use4thOrderAD,advectPassiveScalar
      integer nr(0:1,0:2)
      ! integer bcOptionWallNormal
      integer bc1,bc2,extrapOrder,ks1,kd1,ks2,kd2,is1,is2,is3
      integer axisp1,axisp2,extra,extra1a,extra1b,extra2a,extra2b,
     & extra3a,extra3b
      integer nn1a,nn1b,nn2a,nn2b,nn3a,nn3b
      integer numberOfGhostPoints

      real t,nu,mu,zp,zs,thermalExpansivity,te,zsScaled,rho,realMin
      ! real cd42,adCoeff4, cd22, adCoeff2
      real dr(0:2),dx(0:2),d14v(0:2),d24v(0:2), gravity(0:2)
      ! real vy,vxy,ux,uxy
      real f1,f2,f3,f4,det,alpha,beta,rxsqr,ajs,gamma,alphaMu
      real an1,an2,an3,aNormi, t1,t2,t3, b1,b2,b3, crxd,cryd,crzd, ux,
     & uy, vx,vy,wx,wy,uz,vz,wz, um1,vm1,wm1
      real csf1,csf2,csf3,epsX
      real rxd,ryd,rzd,sxd,syd,szd,txd,tyd,tzd,rxsqd
      real rxxd,ryyd,sxxd,syyd
      real rxi,ryi,rzi,sxi,syi,szi,txi,tyi,tzi,rxxi,ryyi,rzzi
      real a11,a12,a13,a21,a22,a23,a31,a32,a33
      ! real div,div11,div12,div13
      real dive,f1e,f2e,f3e,fMax,res1,res2,res3,resMax
      real c11,c12,c13,c21,c22,c23,c31,c32,c33
      real tvx,tvy,tvz, tvxe,tvye,tvze

      ! variables to hold the exact solution:
      real ue,uxe,uye,uze,uxxe,uyye,uzze,ute
      real ve,vxe,vye,vze,vxxe,vyye,vzze,vte
      real we,wxe,wye,wze,wxxe,wyye,wzze,wte
      real pe,pxe,pye,pze,pxxe,pyye,pzze,pte
      real uxxhe,uyyhe,vxxhe,vyyhe,uthe,vthe
      real phe,pxhe,pyhe,pzhe,pxxhe,pyyhe,pzzhe
      real une,vne,wne

      real uxx,uyy,uzz, vxx,vyy,vzz, wxx,wyy,wzz
      real unxx,vnxx,wnxx, unyy,vnyy,wnyy, unzz,vnzz,wnzz
      real px,py,pz, pnx,pny,pnz
      real rxx,ryy,sxx,syy
      real ume,vme,wme
      integer useJacobiUpdate,checkResults,correctionStage
      real maxDiv,char,maxChar,tDotTaun,maxRes1,maxRes2,res,solidT1,
     & solidT2
      real nDotU,nDotUs,nTauTilden,pStar,pcSwitch,signForTraction
      real nDotUe,nTauTildene,f1tz,f2tz,f3tz,f4tz,zf

      real uExtrap2,uExtrap3,epsu, u1a,u2a, uLim, clim

      real div,div11,div12,div13,dn,tHalf,nuDtBy2,scale
      real rv(0:3),cLap0,cLapg,rxSq,rySq,rzSq,sxSq,sySq,szSq,txSq,tySq,
     & tzSq

      integer body,stateChoice,n,knownSolution,addedMassVelocityBC
      real val(0:10)

      ! enum DeformingBodyStateOptionEnum -- this should match the enum Parameters.h
      integer boundaryPosition,boundaryVelocity,boundaryAcceleration,
     & boundaryTraction,boundaryTractionRate
      parameter(
     &   boundaryPosition    =0,
     &   boundaryVelocity    =1,
     &   boundaryAcceleration=2,
     &   boundaryTraction    =3,
     &   boundaryTractionRate=4 )


      integer numberOfEquations,job
      real a4(1:4,1:4),bv(1:4),ipvt(1:4),rcond,work(1:4), us(1:4)


!..................
      integer rectangular,curvilinear
      parameter( rectangular=0, curvilinear=1 )

      integer doubleDiv,divAndDivN
      parameter( doubleDiv=0, divAndDivN=1 )
      integer
     &     noSlipWall,
     &     inflowWithVelocityGiven,
     &     slipWall,
     &     outflow,
     &     convectiveOutflow,
     &     tractionFree,
     &     inflowWithPandTV,
     &     dirichletBoundaryCondition,
     &     symmetry,
     &     axisymmetric,
     &     freeSurfaceBoundaryCondition,
     &     penaltyBoundaryCondition
      parameter( noSlipWall=1,inflowWithVelocityGiven=2,
     & slipWall=4,outflow=5,convectiveOutflow=14,tractionFree=15,
     & inflowWithPandTV=3,
     &  dirichletBoundaryCondition=12,
     &  symmetry=11,axisymmetric=13, penaltyBoundaryCondition=100,
     &  freeSurfaceBoundaryCondition=31 )

      ! outflowOption values:
      integer extrapolateOutflow,neumannAtOuflow
      parameter( extrapolateOutflow=0,neumannAtOuflow=1 )

      ! declare variables for difference approximations
      ! include 'declareDiffOrder4f.h'
      ! declareDifferenceOrder4(u,RX)

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

! .............. begin statement functions
      real divBCr2d,divBCs2d, divBCr3d,divBCs3d,divBCt3d
      real rx,ry,rz,sx,sy,sz,tx,ty,tz
      real insbfu2d,insbfv2d,insbfu3d,insbfv3d,insbfw3d,ogf
      real delta42,delta43, delta22, delta23

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

!     The next macro call will define the difference approximation statement functions
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
      rxrr2(i1,i2,i3)=(-2.*rx(i1,i2,i3)+(rx(i1+1,i2,i3)+rx(i1-1,i2,i3))
     &  )*d22(0)
      rxss2(i1,i2,i3)=(-2.*rx(i1,i2,i3)+(rx(i1,i2+1,i3)+rx(i1,i2-1,i3))
     &  )*d22(1)
      rxrs2(i1,i2,i3)=(rxr2(i1,i2+1,i3)-rxr2(i1,i2-1,i3))*d12(1)
      ryr2(i1,i2,i3)=(ry(i1+1,i2,i3)-ry(i1-1,i2,i3))*d12(0)
      rys2(i1,i2,i3)=(ry(i1,i2+1,i3)-ry(i1,i2-1,i3))*d12(1)
      ryt2(i1,i2,i3)=(ry(i1,i2,i3+1)-ry(i1,i2,i3-1))*d12(2)
      ryrr2(i1,i2,i3)=(-2.*ry(i1,i2,i3)+(ry(i1+1,i2,i3)+ry(i1-1,i2,i3))
     &  )*d22(0)
      ryss2(i1,i2,i3)=(-2.*ry(i1,i2,i3)+(ry(i1,i2+1,i3)+ry(i1,i2-1,i3))
     &  )*d22(1)
      ryrs2(i1,i2,i3)=(ryr2(i1,i2+1,i3)-ryr2(i1,i2-1,i3))*d12(1)
      rzr2(i1,i2,i3)=(rz(i1+1,i2,i3)-rz(i1-1,i2,i3))*d12(0)
      rzs2(i1,i2,i3)=(rz(i1,i2+1,i3)-rz(i1,i2-1,i3))*d12(1)
      rzt2(i1,i2,i3)=(rz(i1,i2,i3+1)-rz(i1,i2,i3-1))*d12(2)
      rzrr2(i1,i2,i3)=(-2.*rz(i1,i2,i3)+(rz(i1+1,i2,i3)+rz(i1-1,i2,i3))
     &  )*d22(0)
      rzss2(i1,i2,i3)=(-2.*rz(i1,i2,i3)+(rz(i1,i2+1,i3)+rz(i1,i2-1,i3))
     &  )*d22(1)
      rzrs2(i1,i2,i3)=(rzr2(i1,i2+1,i3)-rzr2(i1,i2-1,i3))*d12(1)
      sxr2(i1,i2,i3)=(sx(i1+1,i2,i3)-sx(i1-1,i2,i3))*d12(0)
      sxs2(i1,i2,i3)=(sx(i1,i2+1,i3)-sx(i1,i2-1,i3))*d12(1)
      sxt2(i1,i2,i3)=(sx(i1,i2,i3+1)-sx(i1,i2,i3-1))*d12(2)
      sxrr2(i1,i2,i3)=(-2.*sx(i1,i2,i3)+(sx(i1+1,i2,i3)+sx(i1-1,i2,i3))
     &  )*d22(0)
      sxss2(i1,i2,i3)=(-2.*sx(i1,i2,i3)+(sx(i1,i2+1,i3)+sx(i1,i2-1,i3))
     &  )*d22(1)
      sxrs2(i1,i2,i3)=(sxr2(i1,i2+1,i3)-sxr2(i1,i2-1,i3))*d12(1)
      syr2(i1,i2,i3)=(sy(i1+1,i2,i3)-sy(i1-1,i2,i3))*d12(0)
      sys2(i1,i2,i3)=(sy(i1,i2+1,i3)-sy(i1,i2-1,i3))*d12(1)
      syt2(i1,i2,i3)=(sy(i1,i2,i3+1)-sy(i1,i2,i3-1))*d12(2)
      syrr2(i1,i2,i3)=(-2.*sy(i1,i2,i3)+(sy(i1+1,i2,i3)+sy(i1-1,i2,i3))
     &  )*d22(0)
      syss2(i1,i2,i3)=(-2.*sy(i1,i2,i3)+(sy(i1,i2+1,i3)+sy(i1,i2-1,i3))
     &  )*d22(1)
      syrs2(i1,i2,i3)=(syr2(i1,i2+1,i3)-syr2(i1,i2-1,i3))*d12(1)
      szr2(i1,i2,i3)=(sz(i1+1,i2,i3)-sz(i1-1,i2,i3))*d12(0)
      szs2(i1,i2,i3)=(sz(i1,i2+1,i3)-sz(i1,i2-1,i3))*d12(1)
      szt2(i1,i2,i3)=(sz(i1,i2,i3+1)-sz(i1,i2,i3-1))*d12(2)
      szrr2(i1,i2,i3)=(-2.*sz(i1,i2,i3)+(sz(i1+1,i2,i3)+sz(i1-1,i2,i3))
     &  )*d22(0)
      szss2(i1,i2,i3)=(-2.*sz(i1,i2,i3)+(sz(i1,i2+1,i3)+sz(i1,i2-1,i3))
     &  )*d22(1)
      szrs2(i1,i2,i3)=(szr2(i1,i2+1,i3)-szr2(i1,i2-1,i3))*d12(1)
      txr2(i1,i2,i3)=(tx(i1+1,i2,i3)-tx(i1-1,i2,i3))*d12(0)
      txs2(i1,i2,i3)=(tx(i1,i2+1,i3)-tx(i1,i2-1,i3))*d12(1)
      txt2(i1,i2,i3)=(tx(i1,i2,i3+1)-tx(i1,i2,i3-1))*d12(2)
      txrr2(i1,i2,i3)=(-2.*tx(i1,i2,i3)+(tx(i1+1,i2,i3)+tx(i1-1,i2,i3))
     &  )*d22(0)
      txss2(i1,i2,i3)=(-2.*tx(i1,i2,i3)+(tx(i1,i2+1,i3)+tx(i1,i2-1,i3))
     &  )*d22(1)
      txrs2(i1,i2,i3)=(txr2(i1,i2+1,i3)-txr2(i1,i2-1,i3))*d12(1)
      tyr2(i1,i2,i3)=(ty(i1+1,i2,i3)-ty(i1-1,i2,i3))*d12(0)
      tys2(i1,i2,i3)=(ty(i1,i2+1,i3)-ty(i1,i2-1,i3))*d12(1)
      tyt2(i1,i2,i3)=(ty(i1,i2,i3+1)-ty(i1,i2,i3-1))*d12(2)
      tyrr2(i1,i2,i3)=(-2.*ty(i1,i2,i3)+(ty(i1+1,i2,i3)+ty(i1-1,i2,i3))
     &  )*d22(0)
      tyss2(i1,i2,i3)=(-2.*ty(i1,i2,i3)+(ty(i1,i2+1,i3)+ty(i1,i2-1,i3))
     &  )*d22(1)
      tyrs2(i1,i2,i3)=(tyr2(i1,i2+1,i3)-tyr2(i1,i2-1,i3))*d12(1)
      tzr2(i1,i2,i3)=(tz(i1+1,i2,i3)-tz(i1-1,i2,i3))*d12(0)
      tzs2(i1,i2,i3)=(tz(i1,i2+1,i3)-tz(i1,i2-1,i3))*d12(1)
      tzt2(i1,i2,i3)=(tz(i1,i2,i3+1)-tz(i1,i2,i3-1))*d12(2)
      tzrr2(i1,i2,i3)=(-2.*tz(i1,i2,i3)+(tz(i1+1,i2,i3)+tz(i1-1,i2,i3))
     &  )*d22(0)
      tzss2(i1,i2,i3)=(-2.*tz(i1,i2,i3)+(tz(i1,i2+1,i3)+tz(i1,i2-1,i3))
     &  )*d22(1)
      tzrs2(i1,i2,i3)=(tzr2(i1,i2+1,i3)-tzr2(i1,i2-1,i3))*d12(1)
      ux21(i1,i2,i3,kd)= rx(i1,i2,i3)*ur2(i1,i2,i3,kd)
      uy21(i1,i2,i3,kd)=0
      uz21(i1,i2,i3,kd)=0
      ux22(i1,i2,i3,kd)= rx(i1,i2,i3)*ur2(i1,i2,i3,kd)+sx(i1,i2,i3)*
     & us2(i1,i2,i3,kd)
      uy22(i1,i2,i3,kd)= ry(i1,i2,i3)*ur2(i1,i2,i3,kd)+sy(i1,i2,i3)*
     & us2(i1,i2,i3,kd)
      uz22(i1,i2,i3,kd)=0
      ux23(i1,i2,i3,kd)=rx(i1,i2,i3)*ur2(i1,i2,i3,kd)+sx(i1,i2,i3)*us2(
     & i1,i2,i3,kd)+tx(i1,i2,i3)*ut2(i1,i2,i3,kd)
      uy23(i1,i2,i3,kd)=ry(i1,i2,i3)*ur2(i1,i2,i3,kd)+sy(i1,i2,i3)*us2(
     & i1,i2,i3,kd)+ty(i1,i2,i3)*ut2(i1,i2,i3,kd)
      uz23(i1,i2,i3,kd)=rz(i1,i2,i3)*ur2(i1,i2,i3,kd)+sz(i1,i2,i3)*us2(
     & i1,i2,i3,kd)+tz(i1,i2,i3)*ut2(i1,i2,i3,kd)
      rxx21(i1,i2,i3)= rx(i1,i2,i3)*rxr2(i1,i2,i3)
      rxx22(i1,i2,i3)= rx(i1,i2,i3)*rxr2(i1,i2,i3)+sx(i1,i2,i3)*rxs2(
     & i1,i2,i3)
      rxy22(i1,i2,i3)= ry(i1,i2,i3)*rxr2(i1,i2,i3)+sy(i1,i2,i3)*rxs2(
     & i1,i2,i3)
      rxx23(i1,i2,i3)=rx(i1,i2,i3)*rxr2(i1,i2,i3)+sx(i1,i2,i3)*rxs2(i1,
     & i2,i3)+tx(i1,i2,i3)*rxt2(i1,i2,i3)
      rxy23(i1,i2,i3)=ry(i1,i2,i3)*rxr2(i1,i2,i3)+sy(i1,i2,i3)*rxs2(i1,
     & i2,i3)+ty(i1,i2,i3)*rxt2(i1,i2,i3)
      rxz23(i1,i2,i3)=rz(i1,i2,i3)*rxr2(i1,i2,i3)+sz(i1,i2,i3)*rxs2(i1,
     & i2,i3)+tz(i1,i2,i3)*rxt2(i1,i2,i3)
      ryx22(i1,i2,i3)= rx(i1,i2,i3)*ryr2(i1,i2,i3)+sx(i1,i2,i3)*rys2(
     & i1,i2,i3)
      ryy22(i1,i2,i3)= ry(i1,i2,i3)*ryr2(i1,i2,i3)+sy(i1,i2,i3)*rys2(
     & i1,i2,i3)
      ryx23(i1,i2,i3)=rx(i1,i2,i3)*ryr2(i1,i2,i3)+sx(i1,i2,i3)*rys2(i1,
     & i2,i3)+tx(i1,i2,i3)*ryt2(i1,i2,i3)
      ryy23(i1,i2,i3)=ry(i1,i2,i3)*ryr2(i1,i2,i3)+sy(i1,i2,i3)*rys2(i1,
     & i2,i3)+ty(i1,i2,i3)*ryt2(i1,i2,i3)
      ryz23(i1,i2,i3)=rz(i1,i2,i3)*ryr2(i1,i2,i3)+sz(i1,i2,i3)*rys2(i1,
     & i2,i3)+tz(i1,i2,i3)*ryt2(i1,i2,i3)
      rzx22(i1,i2,i3)= rx(i1,i2,i3)*rzr2(i1,i2,i3)+sx(i1,i2,i3)*rzs2(
     & i1,i2,i3)
      rzy22(i1,i2,i3)= ry(i1,i2,i3)*rzr2(i1,i2,i3)+sy(i1,i2,i3)*rzs2(
     & i1,i2,i3)
      rzx23(i1,i2,i3)=rx(i1,i2,i3)*rzr2(i1,i2,i3)+sx(i1,i2,i3)*rzs2(i1,
     & i2,i3)+tx(i1,i2,i3)*rzt2(i1,i2,i3)
      rzy23(i1,i2,i3)=ry(i1,i2,i3)*rzr2(i1,i2,i3)+sy(i1,i2,i3)*rzs2(i1,
     & i2,i3)+ty(i1,i2,i3)*rzt2(i1,i2,i3)
      rzz23(i1,i2,i3)=rz(i1,i2,i3)*rzr2(i1,i2,i3)+sz(i1,i2,i3)*rzs2(i1,
     & i2,i3)+tz(i1,i2,i3)*rzt2(i1,i2,i3)
      sxx22(i1,i2,i3)= rx(i1,i2,i3)*sxr2(i1,i2,i3)+sx(i1,i2,i3)*sxs2(
     & i1,i2,i3)
      sxy22(i1,i2,i3)= ry(i1,i2,i3)*sxr2(i1,i2,i3)+sy(i1,i2,i3)*sxs2(
     & i1,i2,i3)
      sxx23(i1,i2,i3)=rx(i1,i2,i3)*sxr2(i1,i2,i3)+sx(i1,i2,i3)*sxs2(i1,
     & i2,i3)+tx(i1,i2,i3)*sxt2(i1,i2,i3)
      sxy23(i1,i2,i3)=ry(i1,i2,i3)*sxr2(i1,i2,i3)+sy(i1,i2,i3)*sxs2(i1,
     & i2,i3)+ty(i1,i2,i3)*sxt2(i1,i2,i3)
      sxz23(i1,i2,i3)=rz(i1,i2,i3)*sxr2(i1,i2,i3)+sz(i1,i2,i3)*sxs2(i1,
     & i2,i3)+tz(i1,i2,i3)*sxt2(i1,i2,i3)
      syx22(i1,i2,i3)= rx(i1,i2,i3)*syr2(i1,i2,i3)+sx(i1,i2,i3)*sys2(
     & i1,i2,i3)
      syy22(i1,i2,i3)= ry(i1,i2,i3)*syr2(i1,i2,i3)+sy(i1,i2,i3)*sys2(
     & i1,i2,i3)
      syx23(i1,i2,i3)=rx(i1,i2,i3)*syr2(i1,i2,i3)+sx(i1,i2,i3)*sys2(i1,
     & i2,i3)+tx(i1,i2,i3)*syt2(i1,i2,i3)
      syy23(i1,i2,i3)=ry(i1,i2,i3)*syr2(i1,i2,i3)+sy(i1,i2,i3)*sys2(i1,
     & i2,i3)+ty(i1,i2,i3)*syt2(i1,i2,i3)
      syz23(i1,i2,i3)=rz(i1,i2,i3)*syr2(i1,i2,i3)+sz(i1,i2,i3)*sys2(i1,
     & i2,i3)+tz(i1,i2,i3)*syt2(i1,i2,i3)
      szx22(i1,i2,i3)= rx(i1,i2,i3)*szr2(i1,i2,i3)+sx(i1,i2,i3)*szs2(
     & i1,i2,i3)
      szy22(i1,i2,i3)= ry(i1,i2,i3)*szr2(i1,i2,i3)+sy(i1,i2,i3)*szs2(
     & i1,i2,i3)
      szx23(i1,i2,i3)=rx(i1,i2,i3)*szr2(i1,i2,i3)+sx(i1,i2,i3)*szs2(i1,
     & i2,i3)+tx(i1,i2,i3)*szt2(i1,i2,i3)
      szy23(i1,i2,i3)=ry(i1,i2,i3)*szr2(i1,i2,i3)+sy(i1,i2,i3)*szs2(i1,
     & i2,i3)+ty(i1,i2,i3)*szt2(i1,i2,i3)
      szz23(i1,i2,i3)=rz(i1,i2,i3)*szr2(i1,i2,i3)+sz(i1,i2,i3)*szs2(i1,
     & i2,i3)+tz(i1,i2,i3)*szt2(i1,i2,i3)
      txx22(i1,i2,i3)= rx(i1,i2,i3)*txr2(i1,i2,i3)+sx(i1,i2,i3)*txs2(
     & i1,i2,i3)
      txy22(i1,i2,i3)= ry(i1,i2,i3)*txr2(i1,i2,i3)+sy(i1,i2,i3)*txs2(
     & i1,i2,i3)
      txx23(i1,i2,i3)=rx(i1,i2,i3)*txr2(i1,i2,i3)+sx(i1,i2,i3)*txs2(i1,
     & i2,i3)+tx(i1,i2,i3)*txt2(i1,i2,i3)
      txy23(i1,i2,i3)=ry(i1,i2,i3)*txr2(i1,i2,i3)+sy(i1,i2,i3)*txs2(i1,
     & i2,i3)+ty(i1,i2,i3)*txt2(i1,i2,i3)
      txz23(i1,i2,i3)=rz(i1,i2,i3)*txr2(i1,i2,i3)+sz(i1,i2,i3)*txs2(i1,
     & i2,i3)+tz(i1,i2,i3)*txt2(i1,i2,i3)
      tyx22(i1,i2,i3)= rx(i1,i2,i3)*tyr2(i1,i2,i3)+sx(i1,i2,i3)*tys2(
     & i1,i2,i3)
      tyy22(i1,i2,i3)= ry(i1,i2,i3)*tyr2(i1,i2,i3)+sy(i1,i2,i3)*tys2(
     & i1,i2,i3)
      tyx23(i1,i2,i3)=rx(i1,i2,i3)*tyr2(i1,i2,i3)+sx(i1,i2,i3)*tys2(i1,
     & i2,i3)+tx(i1,i2,i3)*tyt2(i1,i2,i3)
      tyy23(i1,i2,i3)=ry(i1,i2,i3)*tyr2(i1,i2,i3)+sy(i1,i2,i3)*tys2(i1,
     & i2,i3)+ty(i1,i2,i3)*tyt2(i1,i2,i3)
      tyz23(i1,i2,i3)=rz(i1,i2,i3)*tyr2(i1,i2,i3)+sz(i1,i2,i3)*tys2(i1,
     & i2,i3)+tz(i1,i2,i3)*tyt2(i1,i2,i3)
      tzx22(i1,i2,i3)= rx(i1,i2,i3)*tzr2(i1,i2,i3)+sx(i1,i2,i3)*tzs2(
     & i1,i2,i3)
      tzy22(i1,i2,i3)= ry(i1,i2,i3)*tzr2(i1,i2,i3)+sy(i1,i2,i3)*tzs2(
     & i1,i2,i3)
      tzx23(i1,i2,i3)=rx(i1,i2,i3)*tzr2(i1,i2,i3)+sx(i1,i2,i3)*tzs2(i1,
     & i2,i3)+tx(i1,i2,i3)*tzt2(i1,i2,i3)
      tzy23(i1,i2,i3)=ry(i1,i2,i3)*tzr2(i1,i2,i3)+sy(i1,i2,i3)*tzs2(i1,
     & i2,i3)+ty(i1,i2,i3)*tzt2(i1,i2,i3)
      tzz23(i1,i2,i3)=rz(i1,i2,i3)*tzr2(i1,i2,i3)+sz(i1,i2,i3)*tzs2(i1,
     & i2,i3)+tz(i1,i2,i3)*tzt2(i1,i2,i3)
      uxx21(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*urr2(i1,i2,i3,kd)+(rxx22(i1,
     & i2,i3))*ur2(i1,i2,i3,kd)
      uyy21(i1,i2,i3,kd)=0
      uxy21(i1,i2,i3,kd)=0
      uxz21(i1,i2,i3,kd)=0
      uyz21(i1,i2,i3,kd)=0
      uzz21(i1,i2,i3,kd)=0
      ulaplacian21(i1,i2,i3,kd)=uxx21(i1,i2,i3,kd)
      uxx22(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*urr2(i1,i2,i3,kd)+2.*(rx(i1,
     & i2,i3)*sx(i1,i2,i3))*urs2(i1,i2,i3,kd)+(sx(i1,i2,i3)**2)*uss2(
     & i1,i2,i3,kd)+(rxx22(i1,i2,i3))*ur2(i1,i2,i3,kd)+(sxx22(i1,i2,
     & i3))*us2(i1,i2,i3,kd)
      uyy22(i1,i2,i3,kd)=(ry(i1,i2,i3)**2)*urr2(i1,i2,i3,kd)+2.*(ry(i1,
     & i2,i3)*sy(i1,i2,i3))*urs2(i1,i2,i3,kd)+(sy(i1,i2,i3)**2)*uss2(
     & i1,i2,i3,kd)+(ryy22(i1,i2,i3))*ur2(i1,i2,i3,kd)+(syy22(i1,i2,
     & i3))*us2(i1,i2,i3,kd)
      uxy22(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*urr2(i1,i2,i3,kd)+(
     & rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*urs2(i1,
     & i2,i3,kd)+sx(i1,i2,i3)*sy(i1,i2,i3)*uss2(i1,i2,i3,kd)+rxy22(i1,
     & i2,i3)*ur2(i1,i2,i3,kd)+sxy22(i1,i2,i3)*us2(i1,i2,i3,kd)
      uxz22(i1,i2,i3,kd)=0
      uyz22(i1,i2,i3,kd)=0
      uzz22(i1,i2,i3,kd)=0
      ulaplacian22(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2)*urr2(
     & i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*sy(i1,
     & i2,i3))*urs2(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**2)*
     & uss2(i1,i2,i3,kd)+(rxx22(i1,i2,i3)+ryy22(i1,i2,i3))*ur2(i1,i2,
     & i3,kd)+(sxx22(i1,i2,i3)+syy22(i1,i2,i3))*us2(i1,i2,i3,kd)
      uxx23(i1,i2,i3,kd)=rx(i1,i2,i3)**2*urr2(i1,i2,i3,kd)+sx(i1,i2,i3)
     & **2*uss2(i1,i2,i3,kd)+tx(i1,i2,i3)**2*utt2(i1,i2,i3,kd)+2.*rx(
     & i1,i2,i3)*sx(i1,i2,i3)*urs2(i1,i2,i3,kd)+2.*rx(i1,i2,i3)*tx(i1,
     & i2,i3)*urt2(i1,i2,i3,kd)+2.*sx(i1,i2,i3)*tx(i1,i2,i3)*ust2(i1,
     & i2,i3,kd)+rxx23(i1,i2,i3)*ur2(i1,i2,i3,kd)+sxx23(i1,i2,i3)*us2(
     & i1,i2,i3,kd)+txx23(i1,i2,i3)*ut2(i1,i2,i3,kd)
      uyy23(i1,i2,i3,kd)=ry(i1,i2,i3)**2*urr2(i1,i2,i3,kd)+sy(i1,i2,i3)
     & **2*uss2(i1,i2,i3,kd)+ty(i1,i2,i3)**2*utt2(i1,i2,i3,kd)+2.*ry(
     & i1,i2,i3)*sy(i1,i2,i3)*urs2(i1,i2,i3,kd)+2.*ry(i1,i2,i3)*ty(i1,
     & i2,i3)*urt2(i1,i2,i3,kd)+2.*sy(i1,i2,i3)*ty(i1,i2,i3)*ust2(i1,
     & i2,i3,kd)+ryy23(i1,i2,i3)*ur2(i1,i2,i3,kd)+syy23(i1,i2,i3)*us2(
     & i1,i2,i3,kd)+tyy23(i1,i2,i3)*ut2(i1,i2,i3,kd)
      uzz23(i1,i2,i3,kd)=rz(i1,i2,i3)**2*urr2(i1,i2,i3,kd)+sz(i1,i2,i3)
     & **2*uss2(i1,i2,i3,kd)+tz(i1,i2,i3)**2*utt2(i1,i2,i3,kd)+2.*rz(
     & i1,i2,i3)*sz(i1,i2,i3)*urs2(i1,i2,i3,kd)+2.*rz(i1,i2,i3)*tz(i1,
     & i2,i3)*urt2(i1,i2,i3,kd)+2.*sz(i1,i2,i3)*tz(i1,i2,i3)*ust2(i1,
     & i2,i3,kd)+rzz23(i1,i2,i3)*ur2(i1,i2,i3,kd)+szz23(i1,i2,i3)*us2(
     & i1,i2,i3,kd)+tzz23(i1,i2,i3)*ut2(i1,i2,i3,kd)
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
      ulaplacian23(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2+rz(i1,
     & i2,i3)**2)*urr2(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**2+
     & sz(i1,i2,i3)**2)*uss2(i1,i2,i3,kd)+(tx(i1,i2,i3)**2+ty(i1,i2,
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
      uxx23r(i1,i2,i3,kd)=(-2.*u(i1,i2,i3,kd)+(u(i1+1,i2,i3,kd)+u(i1-1,
     & i2,i3,kd)) )*h22(0)
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
      uxxx22r(i1,i2,i3,kd)=(-2.*(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))+(u(
     & i1+2,i2,i3,kd)-u(i1-2,i2,i3,kd)) )*h22(0)*h12(0)
      uyyy22r(i1,i2,i3,kd)=(-2.*(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))+(u(
     & i1,i2+2,i3,kd)-u(i1,i2-2,i3,kd)) )*h22(1)*h12(1)
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
      uLapSq22r(i1,i2,i3,kd)= ( 6.*u(i1,i2,i3,kd)   - 4.*(u(i1+1,i2,i3,
     & kd)+u(i1-1,i2,i3,kd))    +(u(i1+2,i2,i3,kd)+u(i1-2,i2,i3,kd)) )
     & /(dx(0)**4) +( 6.*u(i1,i2,i3,kd)    -4.*(u(i1,i2+1,i3,kd)+u(i1,
     & i2-1,i3,kd))    +(u(i1,i2+2,i3,kd)+u(i1,i2-2,i3,kd)) )/(dx(1)**
     & 4)  +( 8.*u(i1,i2,i3,kd)     -4.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,
     & i3,kd)+u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd))   +2.*(u(i1+1,i2+1,
     & i3,kd)+u(i1-1,i2+1,i3,kd)+u(i1+1,i2-1,i3,kd)+u(i1-1,i2-1,i3,kd)
     & ) )/(dx(0)**2*dx(1)**2)
      uxxx23r(i1,i2,i3,kd)=(-2.*(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))+(u(
     & i1+2,i2,i3,kd)-u(i1-2,i2,i3,kd)) )*h22(0)*h12(0)
      uyyy23r(i1,i2,i3,kd)=(-2.*(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))+(u(
     & i1,i2+2,i3,kd)-u(i1,i2-2,i3,kd)) )*h22(1)*h12(1)
      uzzz23r(i1,i2,i3,kd)=(-2.*(u(i1,i2,i3+1,kd)-u(i1,i2,i3-1,kd))+(u(
     & i1,i2,i3+2,kd)-u(i1,i2,i3-2,kd)) )*h22(1)*h12(2)
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
      uLapSq23r(i1,i2,i3,kd)= ( 6.*u(i1,i2,i3,kd)   - 4.*(u(i1+1,i2,i3,
     & kd)+u(i1-1,i2,i3,kd))    +(u(i1+2,i2,i3,kd)+u(i1-2,i2,i3,kd)) )
     & /(dx(0)**4) +( 6.*u(i1,i2,i3,kd)    -4.*(u(i1,i2+1,i3,kd)+u(i1,
     & i2-1,i3,kd))    +(u(i1,i2+2,i3,kd)+u(i1,i2-2,i3,kd)) )/(dx(1)**
     & 4)  +( 6.*u(i1,i2,i3,kd)    -4.*(u(i1,i2,i3+1,kd)+u(i1,i2,i3-1,
     & kd))    +(u(i1,i2,i3+2,kd)+u(i1,i2,i3-2,kd)) )/(dx(2)**4)  +( 
     & 8.*u(i1,i2,i3,kd)     -4.*(u(i1+1,i2,i3,kd)  +u(i1-1,i2,i3,kd) 
     &  +u(i1  ,i2+1,i3,kd)+u(i1  ,i2-1,i3,kd))   +2.*(u(i1+1,i2+1,i3,
     & kd)+u(i1-1,i2+1,i3,kd)+u(i1+1,i2-1,i3,kd)+u(i1-1,i2-1,i3,kd)) )
     & /(dx(0)**2*dx(1)**2)+( 8.*u(i1,i2,i3,kd)     -4.*(u(i1+1,i2,i3,
     & kd)  +u(i1-1,i2,i3,kd)  +u(i1  ,i2,i3+1,kd)+u(i1  ,i2,i3-1,kd))
     &    +2.*(u(i1+1,i2,i3+1,kd)+u(i1-1,i2,i3+1,kd)+u(i1+1,i2,i3-1,
     & kd)+u(i1-1,i2,i3-1,kd)) )/(dx(0)**2*dx(2)**2)+( 8.*u(i1,i2,i3,
     & kd)     -4.*(u(i1,i2+1,i3,kd)  +u(i1,i2-1,i3,kd)  +u(i1,i2  ,
     & i3+1,kd)+u(i1,i2  ,i3-1,kd))   +2.*(u(i1,i2+1,i3+1,kd)+u(i1,i2-
     & 1,i3+1,kd)+u(i1,i2+1,i3-1,kd)+u(i1,i2-1,i3-1,kd)) )/(dx(1)**2*
     & dx(2)**2)
      unr2(i1,i2,i3,kd)=(un(i1+1,i2,i3,kd)-un(i1-1,i2,i3,kd))*d12(0)
      uns2(i1,i2,i3,kd)=(un(i1,i2+1,i3,kd)-un(i1,i2-1,i3,kd))*d12(1)
      unt2(i1,i2,i3,kd)=(un(i1,i2,i3+1,kd)-un(i1,i2,i3-1,kd))*d12(2)
      unrr2(i1,i2,i3,kd)=(-2.*un(i1,i2,i3,kd)+(un(i1+1,i2,i3,kd)+un(i1-
     & 1,i2,i3,kd)) )*d22(0)
      unss2(i1,i2,i3,kd)=(-2.*un(i1,i2,i3,kd)+(un(i1,i2+1,i3,kd)+un(i1,
     & i2-1,i3,kd)) )*d22(1)
      unrs2(i1,i2,i3,kd)=(unr2(i1,i2+1,i3,kd)-unr2(i1,i2-1,i3,kd))*d12(
     & 1)
      untt2(i1,i2,i3,kd)=(-2.*un(i1,i2,i3,kd)+(un(i1,i2,i3+1,kd)+un(i1,
     & i2,i3-1,kd)) )*d22(2)
      unrt2(i1,i2,i3,kd)=(unr2(i1,i2,i3+1,kd)-unr2(i1,i2,i3-1,kd))*d12(
     & 2)
      unst2(i1,i2,i3,kd)=(uns2(i1,i2,i3+1,kd)-uns2(i1,i2,i3-1,kd))*d12(
     & 2)
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
      unxy22(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*unrr2(i1,i2,i3,kd)+
     & (rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*unrs2(i1,
     & i2,i3,kd)+sx(i1,i2,i3)*sy(i1,i2,i3)*unss2(i1,i2,i3,kd)+rxy22(
     & i1,i2,i3)*unr2(i1,i2,i3,kd)+sxy22(i1,i2,i3)*uns2(i1,i2,i3,kd)
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
      unxy23(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*unrr2(i1,i2,i3,kd)+
     & sx(i1,i2,i3)*sy(i1,i2,i3)*unss2(i1,i2,i3,kd)+tx(i1,i2,i3)*ty(
     & i1,i2,i3)*untt2(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,
     & i2,i3)*sx(i1,i2,i3))*unrs2(i1,i2,i3,kd)+(rx(i1,i2,i3)*ty(i1,i2,
     & i3)+ry(i1,i2,i3)*tx(i1,i2,i3))*unrt2(i1,i2,i3,kd)+(sx(i1,i2,i3)
     & *ty(i1,i2,i3)+sy(i1,i2,i3)*tx(i1,i2,i3))*unst2(i1,i2,i3,kd)+
     & rxy23(i1,i2,i3)*unr2(i1,i2,i3,kd)+sxy23(i1,i2,i3)*uns2(i1,i2,
     & i3,kd)+txy23(i1,i2,i3)*unt2(i1,i2,i3,kd)
      unxz23(i1,i2,i3,kd)=rx(i1,i2,i3)*rz(i1,i2,i3)*unrr2(i1,i2,i3,kd)+
     & sx(i1,i2,i3)*sz(i1,i2,i3)*unss2(i1,i2,i3,kd)+tx(i1,i2,i3)*tz(
     & i1,i2,i3)*untt2(i1,i2,i3,kd)+(rx(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,
     & i2,i3)*sx(i1,i2,i3))*unrs2(i1,i2,i3,kd)+(rx(i1,i2,i3)*tz(i1,i2,
     & i3)+rz(i1,i2,i3)*tx(i1,i2,i3))*unrt2(i1,i2,i3,kd)+(sx(i1,i2,i3)
     & *tz(i1,i2,i3)+sz(i1,i2,i3)*tx(i1,i2,i3))*unst2(i1,i2,i3,kd)+
     & rxz23(i1,i2,i3)*unr2(i1,i2,i3,kd)+sxz23(i1,i2,i3)*uns2(i1,i2,
     & i3,kd)+txz23(i1,i2,i3)*unt2(i1,i2,i3,kd)
      unyz23(i1,i2,i3,kd)=ry(i1,i2,i3)*rz(i1,i2,i3)*unrr2(i1,i2,i3,kd)+
     & sy(i1,i2,i3)*sz(i1,i2,i3)*unss2(i1,i2,i3,kd)+ty(i1,i2,i3)*tz(
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
      unxy23r(i1,i2,i3,kd)=(unx23r(i1,i2+1,i3,kd)-unx23r(i1,i2-1,i3,kd)
     & )*h12(1)
      unzz23r(i1,i2,i3,kd)=(-2.*un(i1,i2,i3,kd)+(un(i1,i2,i3+1,kd)+un(
     & i1,i2,i3-1,kd)) )*h22(2)
      unxz23r(i1,i2,i3,kd)=(unx23r(i1,i2,i3+1,kd)-unx23r(i1,i2,i3-1,kd)
     & )*h12(2)
      unyz23r(i1,i2,i3,kd)=(uny23r(i1,i2,i3+1,kd)-uny23r(i1,i2,i3-1,kd)
     & )*h12(2)
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
      unxxx22r(i1,i2,i3,kd)=(-2.*(un(i1+1,i2,i3,kd)-un(i1-1,i2,i3,kd))+
     & (un(i1+2,i2,i3,kd)-un(i1-2,i2,i3,kd)) )*h22(0)*h12(0)
      unyyy22r(i1,i2,i3,kd)=(-2.*(un(i1,i2+1,i3,kd)-un(i1,i2-1,i3,kd))+
     & (un(i1,i2+2,i3,kd)-un(i1,i2-2,i3,kd)) )*h22(1)*h12(1)
      unxxy22r(i1,i2,i3,kd)=( unxx22r(i1,i2+1,i3,kd)-unxx22r(i1,i2-1,
     & i3,kd))/(2.*dx(1))
      unxyy22r(i1,i2,i3,kd)=( unyy22r(i1+1,i2,i3,kd)-unyy22r(i1-1,i2,
     & i3,kd))/(2.*dx(0))
      unxxxx22r(i1,i2,i3,kd)=(6.*un(i1,i2,i3,kd)-4.*(un(i1+1,i2,i3,kd)+
     & un(i1-1,i2,i3,kd))+(un(i1+2,i2,i3,kd)+un(i1-2,i2,i3,kd)) )/(dx(
     & 0)**4)
      unyyyy22r(i1,i2,i3,kd)=(6.*un(i1,i2,i3,kd)-4.*(un(i1,i2+1,i3,kd)+
     & un(i1,i2-1,i3,kd))+(un(i1,i2+2,i3,kd)+un(i1,i2-2,i3,kd)) )/(dx(
     & 1)**4)
      unxxyy22r(i1,i2,i3,kd)=( 4.*un(i1,i2,i3,kd)     -2.*(un(i1+1,i2,
     & i3,kd)+un(i1-1,i2,i3,kd)+un(i1,i2+1,i3,kd)+un(i1,i2-1,i3,kd))  
     &  +   (un(i1+1,i2+1,i3,kd)+un(i1-1,i2+1,i3,kd)+un(i1+1,i2-1,i3,
     & kd)+un(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)
      ! 2D laplacian squared = un.xxxx + 2 un.xxyy + un.yyyy
      unLapSq22r(i1,i2,i3,kd)= ( 6.*un(i1,i2,i3,kd)   - 4.*(un(i1+1,i2,
     & i3,kd)+un(i1-1,i2,i3,kd))    +(un(i1+2,i2,i3,kd)+un(i1-2,i2,i3,
     & kd)) )/(dx(0)**4) +( 6.*un(i1,i2,i3,kd)    -4.*(un(i1,i2+1,i3,
     & kd)+un(i1,i2-1,i3,kd))    +(un(i1,i2+2,i3,kd)+un(i1,i2-2,i3,kd)
     & ) )/(dx(1)**4)  +( 8.*un(i1,i2,i3,kd)     -4.*(un(i1+1,i2,i3,
     & kd)+un(i1-1,i2,i3,kd)+un(i1,i2+1,i3,kd)+un(i1,i2-1,i3,kd))   +
     & 2.*(un(i1+1,i2+1,i3,kd)+un(i1-1,i2+1,i3,kd)+un(i1+1,i2-1,i3,kd)
     & +un(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)
      unxxx23r(i1,i2,i3,kd)=(-2.*(un(i1+1,i2,i3,kd)-un(i1-1,i2,i3,kd))+
     & (un(i1+2,i2,i3,kd)-un(i1-2,i2,i3,kd)) )*h22(0)*h12(0)
      unyyy23r(i1,i2,i3,kd)=(-2.*(un(i1,i2+1,i3,kd)-un(i1,i2-1,i3,kd))+
     & (un(i1,i2+2,i3,kd)-un(i1,i2-2,i3,kd)) )*h22(1)*h12(1)
      unzzz23r(i1,i2,i3,kd)=(-2.*(un(i1,i2,i3+1,kd)-un(i1,i2,i3-1,kd))+
     & (un(i1,i2,i3+2,kd)-un(i1,i2,i3-2,kd)) )*h22(1)*h12(2)
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
      unxxxx23r(i1,i2,i3,kd)=(6.*un(i1,i2,i3,kd)-4.*(un(i1+1,i2,i3,kd)+
     & un(i1-1,i2,i3,kd))+(un(i1+2,i2,i3,kd)+un(i1-2,i2,i3,kd)) )/(dx(
     & 0)**4)
      unyyyy23r(i1,i2,i3,kd)=(6.*un(i1,i2,i3,kd)-4.*(un(i1,i2+1,i3,kd)+
     & un(i1,i2-1,i3,kd))+(un(i1,i2+2,i3,kd)+un(i1,i2-2,i3,kd)) )/(dx(
     & 1)**4)
      unzzzz23r(i1,i2,i3,kd)=(6.*un(i1,i2,i3,kd)-4.*(un(i1,i2,i3+1,kd)+
     & un(i1,i2,i3-1,kd))+(un(i1,i2,i3+2,kd)+un(i1,i2,i3-2,kd)) )/(dx(
     & 2)**4)
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
      unLapSq23r(i1,i2,i3,kd)= ( 6.*un(i1,i2,i3,kd)   - 4.*(un(i1+1,i2,
     & i3,kd)+un(i1-1,i2,i3,kd))    +(un(i1+2,i2,i3,kd)+un(i1-2,i2,i3,
     & kd)) )/(dx(0)**4) +( 6.*un(i1,i2,i3,kd)    -4.*(un(i1,i2+1,i3,
     & kd)+un(i1,i2-1,i3,kd))    +(un(i1,i2+2,i3,kd)+un(i1,i2-2,i3,kd)
     & ) )/(dx(1)**4)  +( 6.*un(i1,i2,i3,kd)    -4.*(un(i1,i2,i3+1,kd)
     & +un(i1,i2,i3-1,kd))    +(un(i1,i2,i3+2,kd)+un(i1,i2,i3-2,kd)) )
     & /(dx(2)**4)  +( 8.*un(i1,i2,i3,kd)     -4.*(un(i1+1,i2,i3,kd)  
     & +un(i1-1,i2,i3,kd)  +un(i1  ,i2+1,i3,kd)+un(i1  ,i2-1,i3,kd))  
     &  +2.*(un(i1+1,i2+1,i3,kd)+un(i1-1,i2+1,i3,kd)+un(i1+1,i2-1,i3,
     & kd)+un(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)+( 8.*un(i1,i2,
     & i3,kd)     -4.*(un(i1+1,i2,i3,kd)  +un(i1-1,i2,i3,kd)  +un(i1  
     & ,i2,i3+1,kd)+un(i1  ,i2,i3-1,kd))   +2.*(un(i1+1,i2,i3+1,kd)+
     & un(i1-1,i2,i3+1,kd)+un(i1+1,i2,i3-1,kd)+un(i1-1,i2,i3-1,kd)) )
     & /(dx(0)**2*dx(2)**2)+( 8.*un(i1,i2,i3,kd)     -4.*(un(i1,i2+1,
     & i3,kd)  +un(i1,i2-1,i3,kd)  +un(i1,i2  ,i3+1,kd)+un(i1,i2  ,i3-
     & 1,kd))   +2.*(un(i1,i2+1,i3+1,kd)+un(i1,i2-1,i3+1,kd)+un(i1,i2+
     & 1,i3-1,kd)+un(i1,i2-1,i3-1,kd)) )/(dx(1)**2*dx(2)**2)

! .............. end statement functions


      ierr=0

      pc                 =ipar(0)
      uc                 =ipar(1)
      vc                 =ipar(2)
      wc                 =ipar(3)
      sc                 =ipar(4)
      grid               =ipar(5)
      gridType           =ipar(6)
      orderOfAccuracy    =ipar(7)
      gridIsMoving       =ipar(8)
      useWhereMask       =ipar(9)
      gridIsImplicit     =ipar(10)
      implicitMethod     =ipar(11)
      implicitOption     =ipar(12)
      isAxisymmetric     =ipar(13)
      twilightZone       =ipar(14)
      numberOfProcessors =ipar(15)
      debug              =ipar(16)
      myid               =ipar(17)
      assignTemperature  =ipar(18)
      tc                 =ipar(19)
      side               =ipar(20)
      axis               =ipar(21)
      knownSolution      =ipar(22)
      addedMassVelocityBC=ipar(23)
      correctionStage    =ipar(24)

      dx(0)             =rpar(0)
      dx(1)             =rpar(1)
      dx(2)             =rpar(2)
      dr(0)             =rpar(3)
      dr(1)             =rpar(4)
      dr(2)             =rpar(5)
      nu                =rpar(6)
      t                 =rpar(7)
      zp                =rpar(8)
      zs                =rpar(9)
      alpha             =rpar(10) ! added-mass ratio for fluid projection
      mu                =rpar(11)
      ajs               =rpar(12)
      gravity(0)        =rpar(13) ! not used
      gravity(1)        =rpar(14) ! not used
      gravity(2)        =rpar(15) ! not used
      thermalExpansivity=rpar(16) ! not used
      ep                =rpar(17) ! pointer for exact solution
      dt                =rpar(18)
      realMin           =rpar(19)
      rho               =rpar(20) ! fluid density

      epsX = 1.e-30 ! fix me -- pass in

      useJacobiUpdate=1 ! =1 : first compute all new values, then replace

      if( correctionStage.gt.0 )then
        pcSwitch=1.
      else
        pcSwitch=0.
      end if

      checkResults= 0 ! =1 : for debugging, check results

      if( .true. )then
        write(*,'("Inside insExplicitAMPVelocityBC nd=",i2," tz=",i2,",
     & t=",e9.2," mu=",e8.2)') nd,twilightZone,t,mu
        write(*,'("  +++  gridType=",i2," orderOfAccuracy=",i2," (side,
     & axis)=",2i3)') gridType,orderOfAccuracy,side,axis
        write(*,'("  +++  zp=",e10.2," zs=",e10.2," alpha=",e10.2," 
     & dt=",e10.2)') zp,zs,alpha,dt
        write(*,'("  +++  addedMassVelocityBC=",i2," rho=",e12.4)') 
     & addedMassVelocityBC,rho
        write(*,'("  +++  useJacobiUpdate=",i2," pcSwitch=",f3.1)') 
     & useJacobiUpdate,pcSwitch
      end if

      if( mu.le.0. )then
        write(*,'("insExplicitAMPVelocityBC: mu=",e9.2, " must be 
     & POSITIVE!")') mu
        stop 6590
      end if

      extra=0
      numberOfGhostPoints=1
      ! ================= START LOOP OVER SIDES ===============================
       ! *NOTE: extra is not used yet -- keep for future 
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
       if( bc(0,0).lt.0 )then
         extra1a=max(0,extra1a) ! over-ride extra=-1 : assign ends in periodic directions (or internal parallel boundaries)
       else if( bc(0,0).eq.0 )then
         extra1a=numberOfGhostPoints  ! include interpolation points since we assign ghost points outside these
       end if
       ! **NOTE** the bc on the right may be negative even it is not on the left (for parallel)
       if( bc(1,0).lt.0 )then
         extra1b=max(0,extra1b) ! over-ride extra=-1 : assign ends in periodic directions
       else if( bc(1,0).eq.0 )then
         extra1b=numberOfGhostPoints
       end if
       if( bc(0,1).lt.0 )then
         extra2a=max(0,extra2a) ! over-ride extra=-1 : assign ends in periodic directions (or internal parallel boundaries)
       else if( bc(0,1).eq.0 )then
         extra2a=numberOfGhostPoints  ! include interpolation points since we assign ghost points outside these
       end if
       ! **NOTE** the bc on the right may be negative even it is not on the left (for parallel)
       if( bc(1,1).lt.0 )then
         extra2b=max(0,extra2b) ! over-ride extra=-1 : assign ends in periodic directions
       else if( bc(1,1).eq.0 )then
         extra2b=numberOfGhostPoints
       end if
       if(  nd.eq.3 )then
        if( bc(0,2).lt.0 )then
          extra3a=max(0,extra3a) ! over-ride extra=-1 : assign ends in periodic directions (or internal parallel boundaries)
        else if( bc(0,2).eq.0 )then
          extra3a=numberOfGhostPoints  ! include interpolation points since we assign ghost points outside these
        end if
        ! **NOTE** the bc on the right may be negative even it is not on the left (for parallel)
        if( bc(1,2).lt.0 )then
          extra3b=max(0,extra3b) ! over-ride extra=-1 : assign ends in periodic directions
        else if( bc(1,2).eq.0 )then
          extra3b=numberOfGhostPoints
        end if
       end if
      !! do axis=0,nd-1
      !! do side=0,1
      !!   if( bc(side,axis).eq.tractionFree .or. bc(side,axis).eq.freeSurfaceBoundaryCondition )then 
         if( bc(side,axis).eq.noSlipWall )then
           ! write(*,'(" insTractionBC: nd,side,axis,bc=",4i4)') nd,side,axis,bc(side,axis)
           n1a=gridIndexRange(0,0)
           n1b=gridIndexRange(1,0)
           n2a=gridIndexRange(0,1)
           n2b=gridIndexRange(1,1)
           n3a=gridIndexRange(0,2)
           n3b=gridIndexRange(1,2)
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
           nn1a=gridIndexRange(0,0)-extra1a
           nn1b=gridIndexRange(1,0)+extra1b
           nn2a=gridIndexRange(0,1)-extra2a
           nn2b=gridIndexRange(1,1)+extra2b
           nn3a=gridIndexRange(0,2)-extra3a
           nn3b=gridIndexRange(1,2)+extra3b
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
             write(*,'(" insExplicitAMPVelocityBC: grid,side,axis=",
     & 3i3,", loop bounds: n1a,n1b,n2a,n2b,n3a,n3b=",6i3)') grid,side,
     & axis,n1a,n1b,n2a,n2b,n3a,n3b
           end if

       if( nd.eq.2 )then
         ! --- 2D TRACTION  ----
        if( gridType.eq.0 )then
          ! --- RECTANGULAR ----
          write(*,'("insExplicitAMPVelocityBC: finish me")')
          stop 20202

          if( orderOfAccuracy.eq.2 )then
            if( axis.eq.0 )then
              ! tractionRectangular(2,x,2)
            else if( axis.eq.1 )then
              ! tractionRectangular(2,y,2)
            else
              stop 4444
            end if
          else if( orderOfAccuracy.eq.4 )then
            stop 3333
          else
           stop 1234
          end if

        else if( gridType.eq.1 )then
          ! --- CURVILINEAR ----

         if( orderOfAccuracy.eq.2 )then

          if( addedMassVelocityBC.eq.1 )then
            ! *new* version June 25, 2018
            write(*,'("START *NEW* AMP INTERFACE VELOCITY BC LOOPS 
     & CURVILINEAR ")')
            ! First project interface values?
            if( .false. )then
              write(*,'(" INS AMP VELOCITY BC -- pre-project velocity")
     & ')
              do i3=nn3a,nn3b
              do i2=nn2a,nn2b
              do i1=nn1a,nn1b
                ! **FIX ME** There is one alpha(zp) and one alpha(zs)
                u(i1,i2,i3,uc)=alpha*u(i1,i2,i3,uc)+ (1.-alpha)*
     & solidVelocity(i1,i2,i3,0)
                u(i1,i2,i3,vc)=alpha*u(i1,i2,i3,vc)+ (1.-alpha)*
     & solidVelocity(i1,i2,i3,1)
              end do
              end do
              end do
            end if
            maxDiv=0.
            maxChar=0.
            do i3=n3a,n3b
            do i2=n2a,n2b
            do i1=n1a,n1b
             i1m=i1-is1
             i2m=i2-is2
             i3m=i3-is3
             i1p=i1+is1
             i2p=i2+is2
             i3p=i3+is3
             ! *************** TRACTION BC CURVILINEAR GRIDS ****************
             ! (rxd,ryd) : direction of the normal to r(axis)=const
             rxd = rsxy(i1,i2,i3,  axis,0)
             ryd = rsxy(i1,i2,i3,  axis,1)
             sxd = rsxy(i1,i2,i3,axisp1,0)
             syd = rsxy(i1,i2,i3,axisp1,1)
             if( axis.eq.0 )then
               rxxd  = rxx22(i1,i2,i3)
               ryyd  = ryy22(i1,i2,i3)
               sxxd  = sxx22(i1,i2,i3)
               syyd  = syy22(i1,i2,i3)
             else
               rxxd  = sxx22(i1,i2,i3)
               ryyd  = syy22(i1,i2,i3)
               sxxd  = rxx22(i1,i2,i3)
               syyd  = ryy22(i1,i2,i3)
             end if
              an1 = rsxy(i1,i2,i3,axis,0)
              an2 = rsxy(i1,i2,i3,axis,1)
              aNormi = -is/max(epsX,sqrt(an1**2 + an2**2))
              an1=an1*aNormi
              an2=an2*aNormi
             ! tangent
             t1=-an2
             t2= an1
             ux = ux22(i1,i2,i3,uc)
             uy = uy22(i1,i2,i3,uc)
             vx = ux22(i1,i2,i3,vc)
             vy = uy22(i1,i2,i3,vc)
             uxx = uxx22(i1,i2,i3,uc)
             uyy = uyy22(i1,i2,i3,uc)
             vxx = uxx22(i1,i2,i3,vc)
             vyy = uyy22(i1,i2,i3,vc)
             px  = ux22(i1,i2,i3,pc)
             py  = uy22(i1,i2,i3,pc)
             unxx = unxx22(i1,i2,i3,uc)
             unyy = unyy22(i1,i2,i3,uc)
             vnxx = unxx22(i1,i2,i3,vc)
             vnyy = unyy22(i1,i2,i3,vc)
             pnx  = unx22(i1,i2,i3,pc)
             pny  = uny22(i1,i2,i3,pc)
             if( .FALSE. .and. twilightZone.eq.1 )then
               ! Test -- set exact for px,py 
               call ogDeriv(ep,0,1,0,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & t-dt,pc,pxe )
               call ogDeriv(ep,0,0,1,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & t-dt,pc,pye )
               pnx = pxe
               pny = pye
               call ogDeriv(ep,0,1,0,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & t,pc,pxe )
               call ogDeriv(ep,0,0,1,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & t,pc,pye )
               px = pxe
               py = pye
             end if
             ! Us minus of solid traction since the solid normal is reversed
             solidT1 = -solidTraction(i1,i2,i3,0)
             solidT2 = -solidTraction(i1,i2,i3,1)
             ! ========== Equation 1:   ==========
             !        2*mu*div(v) + zp*n.v = pcSwitch*( n.sigmaS.n - ( -p + n.tauTilde.n ) ) + zp*n.vS 
             !
             !   tauTilde = tau - 2*mu*div(v)*I 
             !   pcSwitch = 0 for predictor stage, 1 for corrector stage 
             !
             ! First evaluate the equation using current ghost values 
             !   f1 = ux+vy = a11*u(-1) + a12*v(-1) + a13*u(0) + a14*v(0) = rest
             !   rest = f1(uCurrent) - a11*uCurrent(-1) + a12*vCurrent(-1)
             nDotU = an1*u(i1,i2,i3,uc)+an2*u(i1,i2,i3,vc)
             nDotUs= an1*solidVelocity(i1,i2,i3,0) + an2*solidVelocity(
     & i1,i2,i3,1)
             nTauTilden = mu*( an1*( -2.*vy*an1 + (uy+vx)*an2) + an2*( 
     & (uy+vx)*an1 -2.*ux*an2 ) )
             pStar = u(i1,i2,i3,pc) ! current value for p
             signForTraction=1.
             ! "equation" 1: (all terms on the LHS)
             ! f1 = 2.*mu*(ux+vy) 
             f1 = 2.*mu*(ux+vy) + zp*nDotU - zp*nDotUs + pcSwitch*
     & signForTraction*( -pStar+nTauTilden  -( an1*solidT1+an2*
     & solidT2 ) )
             if( twilightZone.eq.1 )then
               ! adjust for TZ 
               call ogDeriv(ep,0,0,0,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & t,uc,ue )
               call ogDeriv(ep,0,0,0,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & t,vc,ve )
               call ogDeriv(ep,0,0,0,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & t,pc,pe )
               call ogDeriv(ep,0,1,0,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & t,uc,uxe )
               call ogDeriv(ep,0,1,0,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & t,vc,vxe )
               call ogDeriv(ep,0,0,1,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & t,uc,uye )
               call ogDeriv(ep,0,0,1,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & t,vc,vye )
               nDotUe = an1*ue+an2*ve
               nTauTildene = mu*( an1*( -2*vye*an1 + (uye+vxe)*an2) + 
     & an2*( (uye+vxe)*an1 -2.*uxe*an2 ) )
               f1tz = 2.*mu*(uxe+vye) + zp*nDotUe + pcSwitch*(-pe + 
     & nTauTildene )
               ! f1tz = 2.*mu*(uxe+vye) 
               f1 = f1 - f1tz
             end if
             ! crxd = coeff of u(-1) in u.x  
             ! cryd = coeff of u(-1) in u.y  
             crxd=-is*rxd/(2.*dr(axis))
             cryd=-is*ryd/(2.*dr(axis))
             a4(1,1) = 2.*mu*crxd             ! coeff of u(-1) in eqn 1
             a4(1,2) = 2.*mu*cryd             ! coeff of v(-1) in eqn 1
             a4(1,3) = zp*an1                 ! coeff of u( 0) in eqn 1
             a4(1,4) = zp*an2                 ! coeff of v( 0) in eqn 1
             ! RHS: 
             bv(1) = a4(1,1)*u(i1m,i2m,i3m,uc) + a4(1,2)*u(i1m,i2m,i3m,
     & vc) + a4(1,3)*u(i1 ,i2 ,i3 ,uc) + a4(1,4)*u(i1 ,i2 ,i3 ,vc) - 
     & f1
             ! assume the TZ solution is divergence free so we do not need to add a TZ correction to f1
             ! =========== Equation 2:  tangential traction equation =================
             !    (1/mu)* tv.tau.n + (zs/mu)*tv.v = (1/mu)*tv.(-solidTraction) + (zs/mu)*tv.solidVelocity
             ! 
             ! First evaluate the tangential traction equation using current ghost values 
             !  f2 = (1/mu) * tv.tauv.nv 
             !     =  2*ux t1*n1 + (uy+vx)*(t1*n2+t2*n1) + 2* t2*n2* vy 
             !     = csf1*ux + csf2*(uy+vx) + csf3*vy
             !     = a21*u(-1) + a22*v(-1) + .... = f2
            ! zsScaled =zs/(1.+10./zs) ! set to zero to turn off "zs" terms
             zsScaled =zs
             csf1= 2.*t1*an1
             csf2=(t1*an2+t2*an1)
             csf3= 2.*t2*an2
             a4(2,1) = csf1*crxd + csf2*cryd  ! coeff of u(-1)
             a4(2,2) = csf2*crxd + csf3*cryd  ! coeff of v(-1)
             a4(2,3) = (zsScaled/mu)*t1       ! coeff of u( 0)
             a4(2,4) = (zsScaled/mu)*t2       ! coeff of v( 0)
             ! Eqn 2: (leave off (zs/mu)*tv.v since it will be otherwise just subtracted below) 
             f2 = csf1*ux + csf2*(uy+vx) + csf3*vy - (1./mu)*( t1*
     & solidT1                   + t2*solidT2                   ) - (
     & signForTraction*zsScaled/mu)*( t1*solidVelocity(i1,i2,i3,0) + 
     & t2*solidVelocity(i1,i2,i3,1) )
             if( twilightZone.eq.1 )then
               call ogDeriv(ep,0,0,0,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & t,uc,ue )
               call ogDeriv(ep,0,0,0,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & t,vc,ve )
               call ogDeriv(ep,0,0,0,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & t-dt,uc,une )
               call ogDeriv(ep,0,0,0,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & t-dt,vc,vne )
               call ogDeriv(ep,0,1,0,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & t,uc,uxe )
               call ogDeriv(ep,0,1,0,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & t,vc,vxe )
               call ogDeriv(ep,0,0,1,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & t,uc,uye )
               call ogDeriv(ep,0,0,1,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & t,vc,vye )
               ! Adjust for TZ:  
               f2 = f2 - ( (csf1*uxe + csf2*(uye+vxe) + csf3*vye) + (
     & zsScaled/mu)*( t1*ue+t2*ve ) )
             end if
             bv(2) = a4(2,1)*u(i1m,i2m,i3m,uc) + a4(2,2)*u(i1m,i2m,i3m,
     & vc) - f2
             ! =========== Equation 3: tangential velocity equation ================
             !          tv.v - (nu*dt)/2*tv.Delta v = tv.rv 
             !          rv = vn - dt*( .5*grad(p) + .5*grad(pn) ) + (nu*dt/2)*(Delta un )
             !  pn : solution at the old time 
             !  (on the predictor stange we assume "u" holds a predicted pressure) 
             nuDtby2 = .5*nu*dt
             ! Delta = (rx^2+ry^2)*urr + (sx^2+sy^2)*uss + (rxx+ryy)*ur + (sxx+syy)*us 
             rxSq = rxd**2 + ryd**2
             sxSq = sxd**2 + syd**2
             cLapg = rxSq/dr(axis)**2 -is*(rxxd+ryyd)/(2.*dr(axis))   ! coeff of ghost in Delta()
             cLap0 = -2.*( rxSq/dr(axis)**2 + sxSq/dr(axisp1)**2 )    ! coeff of bndry in Delta()
             scale=1.  ! Set to zero to turn on interior equation **TEMP**
             a4(3,1) =    -scale*nuDtby2*cLapg*t1   ! coeff of u_{-1}
             a4(3,2) =    -scale*nuDtby2*cLapg*t2   ! coeff of v_{-1}
             a4(3,3) = t1 -scale*nuDtby2*cLap0*t1   ! coeff of u_{0}
             a4(3,4) = t2 -scale*nuDtby2*cLap0*t2   ! coeff of v_{0}
             rv(0) = un(i1,i2,i3,uc) -(.5*dt/rho)*( px + pnx ) + 
     & nuDtby2*( unxx + unyy + uxx+uyy )
             rv(1) = un(i1,i2,i3,vc) -(.5*dt/rho)*( py + pny ) + 
     & nuDtby2*( vnxx + vnyy + vxx+vyy )
             ! adjust RHS for evaluating derivatives with wrong values on the ghost and boundary
             !rv(0) = rv(0) - nuDtby2*( cLapg*u(i1m,i2m,i3m,uc) + cLap0*u(i1,i2,i3,uc) )
             !rv(1) = rv(1) - nuDtby2*( cLapg*u(i1m,i2m,i3m,vc) + cLap0*u(i1,i2,i3,vc) )
             if( twilightZone.eq.1 )then
               ! Eval TZ forcing at t-dt/2 (centered in [t-dt,t]
               thalf=t-.5*dt
               call ogDeriv(ep,1,0,0,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & tHalf,uc,uthe )
               call ogDeriv(ep,1,0,0,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & tHalf,vc,vthe )
               call ogDeriv(ep,0,2,0,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & tHalf,uc,uxxhe )
               call ogDeriv(ep,0,0,2,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & tHalf,uc,uyyhe )
               call ogDeriv(ep,0,2,0,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & tHalf,vc,vxxhe )
               call ogDeriv(ep,0,0,2,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & tHalf,vc,vyyhe )
               call ogDeriv(ep,0,1,0,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & tHalf,pc,pxhe )
               call ogDeriv(ep,0,0,1,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & tHalf,pc,pyhe )
               rv(0) = rv(0) + dt*uthe + (dt/rho)*( pxhe ) - (nu*dt)*( 
     & uxxhe + uyyhe )
               rv(1) = rv(1) + dt*vthe + (dt/rho)*( pyhe ) - (nu*dt)*( 
     & vxxhe + vyyhe )
             end if
             bv(3) = t1*(-u(i1,i2,i3,uc)+rv(0)) + t2*(-u(i1,i2,i3,vc)+
     & rv(1))
             ! adjust RHS for guessed values at boundary and ghost points:
             bv(3) =   a4(3,1)*u(i1m,i2m,i3m,uc) + a4(3,2)*u(i1m,i2m,
     & i3m,vc) + a4(3,3)*u(i1 ,i2 ,i3 ,uc) + a4(3,4)*u(i1 ,i2 ,i3 ,vc)
     &   + bv(3)
             ! =========== Equation 4:  normal velocity equation + divergence ================
             !          n.( v- rv) - gamma*div(v) = 0 
             !          rv = vn - dt*( (1/2)*grad(p) + .5*grad(pn) ) + (nu*dt/2)*(Delta un )
             !  pn : solution at the old time 
             !  (on the predictor stage we assume "u" holds a predicted pressure) 
             gamma = zp*dt/rho
             beta =  .5*nu*dt
             a4(4,1) =     -beta*cLapg*an1 -gamma*crxd   ! coeff of u_{-1}
             a4(4,2) =     -beta*cLapg*an2 -gamma*cryd   ! coeff of v_{-1}
             a4(4,3) = an1 -beta*cLap0*an1               ! coeff of u_{0}
             a4(4,4) = an2 -beta*cLap0*an2               ! coeff of v_{0}
             ! RHS: 
             bv(4)= an1*( -u(i1,i2,i3,uc) + rv(0) ) + an2*( -u(i1,i2,
     & i3,vc) +rv(1) ) + gamma*( ux + vy )
             ! adjust RHS for guessed values at boundary and ghost points:
             bv(4) =   a4(4,1)*u(i1m,i2m,i3m,uc) + a4(4,2)*u(i1m,i2m,
     & i3m,vc) + a4(4,3)*u(i1 ,i2 ,i3 ,uc) + a4(4,4)*u(i1 ,i2 ,i3 ,vc)
     &   + bv(4)
             if( twilightZone.eq.1 )then
               ! write(*,'(" an1*rv(0) + an2*rv(1) - an1*ue + an2*ve=",1pe12.3)') an1*rv(0) + an2*rv(1)- (an1*ue + an2*ve)
               ! Adjust for TZ:  
               ! With TZ, the RHS above should equal (-1+alpha)*( n.ve ) so we need to add (1-alpha)*n.ve
               !   n.v = alpha*n.ve + (1-alpha)n.ve = n.ve 
               bv(4) = bv(4) -gamma*( uxe + vye )
             end if
             ! write(*,'(" i1,i2=",2i3," rxd,ryd=",2f8.4," sxd,syd=",2f10.4)') i1,i2,rxd,ryd,sxd,syd
             ! write(*,'(" i1,i2=",2i3," a11,a12=",2f8.4," a21,a22=",2f10.4)') i1,i2,a11,a12,a21,a22
             ! write(*,'(" i1,i2=",2i3," solidTraction=",2e10.2)') i1,i2,solidTraction(i1,i2,i3,0),solidTraction(i1,i2,i3,1)
             ! Solve
             !   [a11 a12 a13 a14 ][ u(-1)] = [ b1 ]
             !   [a21 a22 a23 a24 ][ v(-1)] = [ b2 ]
             !   [a31 a32 a33 a34 ][ u( 0)] = [ b3 ]
             !   [a41 a42 a43 a44 ][ v( 0)] = [ b4 ]
             job=0
             numberOfEquations=4
             call dgeco( a4(1,1), numberOfEquations, numberOfEquations,
     &  ipvt(1),rcond,work(1))
             call dgesl( a4(1,1), numberOfEquations, numberOfEquations,
     &  ipvt(1), bv(1), job)
             ! *********************** CHECK RESULTS ***************************
             if( checkResults.eq.1 )then
               ! check that the equations are all satisfied
               write(*,'("AMP EXPLICIT VELOCITY BC: i1,i2=",2i3," zs=",
     & e16.8," alpha=",e16.8," beta=",e10.2," rcond=",e10.2)') i1,i2,
     & zs,alpha,beta,rcond
               write(*,'(" normal=",2e14.4," tangent=",2e14.4)') an1,
     & an2,t1,t2
               write(*,'(" px,py=",2e14.4," old: px,py=",2e14.4)') px,
     & py,pnx,pny
               write(*,'(" vSolid=",2e16.6," solidTraction=",2e16.6)') 
     & solidVelocity(i1,i2,i3,0),solidVelocity(i1,i2,i3,1),
     & solidTraction(i1,i2,i3,0),solidTraction(i1,i2,i3,1)
               write(*,'(" old: uxx,uyy,vxx,vyy==",4e16.6)') unxx,unyy,
     & vnxx,vnyy
               ! save current values
               us(1)=u(i1m,i2m,i3m,uc)
               us(2)=u(i1m,i2m,i3m,vc)
               us(3)=u(i1 ,i2 ,i3 ,uc)
               us(4)=u(i1 ,i2 ,i3 ,vc)
               u(i1m,i2m,i3m,uc)=bv(1)
               u(i1m,i2m,i3m,vc)=bv(2)
               u(i1 ,i2 ,i3 ,uc)=bv(3)
               u(i1 ,i2 ,i3 ,vc)=bv(4)
               write(*,'(" Old: u,v=",2e16.6," ghost=",2e16.6)') un(i1 
     & ,i2 ,i3 ,uc),un(i1 ,i2 ,i3 ,vc),un(i1m,i2m,i3m,uc),un(i1m,i2m,
     & i3m,vc)
               write(*,'(" New: u,v=",2e16.6," ghost=",2e16.6)') u(i1 ,
     & i2 ,i3 ,uc),u(i1 ,i2 ,i3 ,vc),u(i1m,i2m,i3m,uc),u(i1m,i2m,i3m,
     & vc)
               ux = ux22(i1,i2,i3,uc)
               uy = uy22(i1,i2,i3,uc)
               vx = ux22(i1,i2,i3,vc)
               vy = uy22(i1,i2,i3,vc)
               uxx = uxx22(i1,i2,i3,uc)
               uyy = uyy22(i1,i2,i3,uc)
               vxx = uxx22(i1,i2,i3,vc)
               vyy = uyy22(i1,i2,i3,vc)
               ! =======================  EQN 1 ==============================================
               ! 2*mu*div(v) + zp*n.v = n.sigmaS.n + zp*n.vS - pcSwitch*( -p + n.tauTilde.n ) 
               nDotU = an1*u(i1,i2,i3,uc)+an2*u(i1,i2,i3,vc)
               nDotUs= an1*solidVelocity(i1,i2,i3,0) + an2*
     & solidVelocity(i1,i2,i3,1)
               nTauTilden = mu*( an1*( -2.*vy*an1 + (uy+vx)*an2) + an2*
     & ( (uy+vx)*an1 -2.*ux*an2 ) )
               pStar = u(i1,i2,i3,pc) ! current value for p
               ! "equation" 1: (all terms on the LHS)
               f1 = 2.*mu*(ux+vy) + zp*an1*nDotU  -( an1*solidT1+an2*
     & solidT2 + zp*nDotUs ) + pcSwitch*(-pStar + nTauTilden )
               ! f1 = 2.*mu*(ux+vy) 
               if( twilightZone.eq.1 )then
                 ! adjust for TZ 
                 f1 = f1 - f1tz
               end if
               res=abs(f1)
               maxRes1=max(maxRes1,res)
               ! div=ux+vy
               ! maxDiv=max(maxDiv,abs(div))
               write(*,'("EQN1: After: 2*mu*div(v) + zp*n.v ...=",
     & e10.2)') res
               ! =======================   EQN 2 ==============================================
               ! (1/mu)* tv.sigma.n + (zs/mu)*tv.v = (1/mu)*tv.solidTraction + (zs/mu)*tv.solidVelocity
               tDotTaun =  mu*( t1*(2.*ux*an1+(uy+vx)*an2) + t2*((uy+
     & vx)*an1+2.*vy*an2) )
               char = (1./mu)*(  tDotTaun + zs*(t1*u(i1 ,i2 ,i3 ,uc)+
     & t2*u(i1 ,i2 ,i3 ,vc)) )
               f2 =  (1./mu)*( t1*solidT1 + t2*solidT2 ) + (zs/mu)*( 
     & t1*solidVelocity(i1,i2,i3,0) + t2*solidVelocity(i1,i2,i3,1) )
               ! Adjust for TZ ? June 25, 2018
               if( twilightZone.eq.1 )then
                call ogDeriv(ep,0,0,0,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & t,uc,ue )
                call ogDeriv(ep,0,0,0,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & t,vc,ve )
                call ogDeriv(ep,0,1,0,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & t,uc,uxe )
                call ogDeriv(ep,0,1,0,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & t,vc,vxe )
                call ogDeriv(ep,0,0,1,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & t,uc,uye )
                call ogDeriv(ep,0,0,1,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & t,vc,vye )
                ! Adjust for TZ:  
                f2 = f2 + (csf1*uxe + csf2*(uye+vxe) + csf3*vye) + (
     & zsScaled/mu)*( t1*ue+t2*ve )
               end if
               res = abs(char-f2)
               maxChar=max(maxChar,res)
               write(*,'("EQN2: After: ( t.tau.n + zs*t.v ) ...=",
     & e10.2,"  (char,f2=",2e10.2,")")') res,char,f2
               ! =======================   EQN 3 ==============================================
               !
               ! Here is the solution on the boundary from the interior equations
               rv(0) = un(i1,i2,i3,uc) -(.5*dt/rho)*( px + pnx ) + 
     & nuDtby2*( unxx + unyy + uxx+uyy )
               rv(1) = un(i1,i2,i3,vc) -(.5*dt/rho)*( py + pny ) + 
     & nuDtby2*( vnxx + vnyy + vxx+vyy )
               if( twilightZone.eq.1 )then
                 rv(0) = rv(0) + dt*uthe + (dt/rho)*( pxhe ) - (nu*dt)*
     & ( uxxhe + uyyhe )
                 rv(1) = rv(1) + dt*vthe + (dt/rho)*( pyhe ) - (nu*dt)*
     & ( vxxhe + vyyhe )
               end if
               res = abs( t1*rv(0) + t2*rv(1) - (t1*u(i1 ,i2 ,i3 ,uc) +
     &  t2*u(i1 ,i2 ,i3 ,vc) ) )
               maxRes1=max(maxRes1,res)
               write(*,'("EQN3: After: res in tv.( EQN) ) ...=",e10.2)
     & ') res
               ! =========== Equation 4:  normal velocity equation + divergence  ================
               !     n.( v-rv ) - gamma*div(v) = 0 
               res = an1*( -u(i1,i2,i3,uc) + rv(0) ) + an2*( -u(i1,i2,
     & i3,vc) +rv(1) ) + gamma*( ux + vy )
               if( twilightZone.eq.1 )then
                 ! Adjust for TZ:  
                 res = res  -gamma*( uxe + vye )
               end if
               ! res = an1*un(i1,i2,i3,uc) + an2*un(i1,i2,i3,vc) -res 
               maxRes2=max(maxRes2,res)
               write(*,'("EQN4: After: res in nv.( EQN) ) ...=",e10.2)
     & ') res
               ! reset
               u(i1m,i2m,i3m,uc)=us(1)
               u(i1m,i2m,i3m,vc)=us(2)
               u(i1 ,i2 ,i3 ,uc)=us(3)
               u(i1 ,i2 ,i3 ,vc)=us(4)
             end if
             if( .false. .and. knownSolution.ne.0 )then
              body=0
              stateChoice=boundaryVelocity
              n=0
              call evalUserDefinedDeformingBodyKnownSolution(  body,
     & stateChoice,t,grid,i1m,i2m,i3m,n,val(0) )
              call evalUserDefinedDeformingBodyKnownSolution(  body,
     & stateChoice,t,grid,i1 ,i2 ,i3 ,n,val(2) )
              write(*,'(" i1,i2=",2i3," alpha=",e10.2," beta=",e10.2," 
     & rcond=",e10.2)') i1,i2,alpha,beta,rcond
              write(*,'("  insExpAmpBC: (knownSolution): v-err=",4(
     & 1pe10.2)," |u|,|v|=",2e10.2)') bv(1)-val(0),bv(2)-val(1),bv(3)-
     & val(2),bv(4)-val(3),abs(val(0)),abs(val(1))
              ! set exact 
              if( .true. )then
                write(*,'("  insExpAmpBC: **TEMP** set EXACT velocity 
     & on the interface")')
                bv(1)=val(0)
                bv(2)=val(1)
                bv(3)=val(2)
                bv(4)=val(3)
              end if
             end if
             if( useJacobiUpdate.eq.1 )then
               uNew(i1m,i2m,i3m,uc)=bv(1)
               uNew(i1m,i2m,i3m,vc)=bv(2)
               uNew(i1 ,i2 ,i3 ,uc)=bv(3)
               uNew(i1 ,i2 ,i3 ,vc)=bv(4)
             else
               u(i1m,i2m,i3m,uc)=bv(1)
               u(i1m,i2m,i3m,vc)=bv(2)
               u(i1 ,i2 ,i3 ,uc)=bv(3)
               u(i1 ,i2 ,i3 ,vc)=bv(4)
             end if
            end do
            end do
            end do
            ! Now copy the copy into u 
            if( useJacobiUpdate.eq.1 )then
             do i3=n3a,n3b
             do i2=n2a,n2b
             do i1=n1a,n1b
              i1m=i1-is1
              i2m=i2-is2
              i3m=i3-is3
              u(i1m,i2m,i3m,uc)=uNew(i1m,i2m,i3m,uc)
              u(i1m,i2m,i3m,vc)=uNew(i1m,i2m,i3m,vc)
              u(i1 ,i2 ,i3 ,uc)=uNew(i1 ,i2 ,i3 ,uc)
              u(i1 ,i2 ,i3 ,vc)=uNew(i1 ,i2 ,i3 ,vc)
             end do
             end do
             end do
            end if
            if( checkResults.eq.1 )then
              write(*,'(" maxDiv=",e10.2," maxChar=",e10.2," maxRes1=",
     & e10.2," maxRes2=",e10.2)') maxDiv,maxChar,maxRes1,maxRes2
              stop 1111
            end if
          else if( addedMassVelocityBC.eq.0 )then
            ! *old* version
            write(*,'("START *OLD* AMP INTERFACE VELOCITY BC LOOPS 
     & CURVILINEAR ")')
            ! First project interface values?
            if( .false. )then
              write(*,'(" INS AMP VELOCITY BC -- pre-project velocity")
     & ')
              do i3=nn3a,nn3b
              do i2=nn2a,nn2b
              do i1=nn1a,nn1b
                ! **FIX ME** There is one alpha(zp) and one alpha(zs)
                u(i1,i2,i3,uc)=alpha*u(i1,i2,i3,uc)+ (1.-alpha)*
     & solidVelocity(i1,i2,i3,0)
                u(i1,i2,i3,vc)=alpha*u(i1,i2,i3,vc)+ (1.-alpha)*
     & solidVelocity(i1,i2,i3,1)
              end do
              end do
              end do
            end if
            maxDiv=0.
            maxChar=0.
            do i3=n3a,n3b
            do i2=n2a,n2b
            do i1=n1a,n1b
             i1m=i1-is1
             i2m=i2-is2
             i3m=i3-is3
             i1p=i1+is1
             i2p=i2+is2
             i3p=i3+is3
             ! *************** TRACTION BC CURVILINEAR GRIDS ****************
             ! (rxd,ryd) : direction of the normal to r(axis)=const
             rxd = rsxy(i1,i2,i3,  axis,0)
             ryd = rsxy(i1,i2,i3,  axis,1)
             sxd = rsxy(i1,i2,i3,axisp1,0)
             syd = rsxy(i1,i2,i3,axisp1,1)
             if( axis.eq.0 )then
               rxxd  = rxx22(i1,i2,i3)
               ryyd  = ryy22(i1,i2,i3)
               sxxd  = sxx22(i1,i2,i3)
               syyd  = syy22(i1,i2,i3)
             else
               rxxd  = sxx22(i1,i2,i3)
               ryyd  = syy22(i1,i2,i3)
               sxxd  = rxx22(i1,i2,i3)
               syyd  = ryy22(i1,i2,i3)
             end if
              an1 = rsxy(i1,i2,i3,axis,0)
              an2 = rsxy(i1,i2,i3,axis,1)
              aNormi = -is/max(epsX,sqrt(an1**2 + an2**2))
              an1=an1*aNormi
              an2=an2*aNormi
             ! tangent
             t1=-an2
             t2= an1
             ux = ux22(i1,i2,i3,uc)
             uy = uy22(i1,i2,i3,uc)
             vx = ux22(i1,i2,i3,vc)
             vy = uy22(i1,i2,i3,vc)
             uxx = uxx22(i1,i2,i3,uc)
             uyy = uyy22(i1,i2,i3,uc)
             vxx = uxx22(i1,i2,i3,vc)
             vyy = uyy22(i1,i2,i3,vc)
             px  = ux22(i1,i2,i3,pc)
             py  = uy22(i1,i2,i3,pc)
             unxx = unxx22(i1,i2,i3,uc)
             unyy = unyy22(i1,i2,i3,uc)
             vnxx = unxx22(i1,i2,i3,vc)
             vnyy = unyy22(i1,i2,i3,vc)
             pnx  = unx22(i1,i2,i3,pc)
             pny  = uny22(i1,i2,i3,pc)
             if( .FALSE. .and. twilightZone.eq.1 )then
               ! Test -- set exact for px,py ****TEMP****
               call ogDeriv(ep,0,1,0,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & t-dt,pc,pxe )
               call ogDeriv(ep,0,0,1,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & t-dt,pc,pye )
               pnx = pxe
               pny = pye
               call ogDeriv(ep,0,1,0,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & t,pc,pxe )
               call ogDeriv(ep,0,0,1,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & t,pc,pye )
               px = pxe
               py = pye
             end if
             ! Us minus of solid traction since the solid normal is reversed
             solidT1 = -solidTraction(i1,i2,i3,0)
             solidT2 = -solidTraction(i1,i2,i3,1)
             ! ========== Equation 1: divergence=0  ==========
             ! crxd = coeff of u(-1) in u.x  
             ! cryd = coeff of u(-1) in u.y  
             crxd=-is*rxd/(2.*dr(axis))
             cryd=-is*ryd/(2.*dr(axis))
             ! First evaluate div(u) using current ghost values 
             !   f1 = ux+vy = a11*u(-1) + a12*v(-1) + a13*u(0) + a14*v(0) = rest
             !   rest = f1(uCurrent) - a11*uCurrent(-1) + a12*vCurrent(-1)
             f1 = ux+vy  ! "equation" 1
             a4(1,1) = -is*rxd/(2.*dr(axis))  ! coeff of u(-1) in div(v)=0
             a4(1,2) = -is*ryd/(2.*dr(axis))  ! coeff of v(-1) in div(v)=0
             a4(1,3) = 0.                     ! coeff of u( 0) in div(v)=0
             a4(1,4) = 0.                     ! coeff of v( 0) in div(v)=0
             bv(1) = a4(1,1)*u(i1m,i2m,i3m,uc) + a4(1,2)*u(i1m,i2m,i3m,
     & vc) - f1  ! RHS
             ! assume the TZ solution is divergence free so we do not need to add a TZ correction to f1
             ! =========== Equation 2:  tangential traction equation =================
             !    (1/mu)* tv.sigma.n + (zs/mu)*tv.v = (1/mu)*tv.(-solidTraction) + (zs/mu)*tv.solidVelocity
             ! 
             ! First evaluate the tangential traction equation using current ghost values 
             !  f2 = (1/mu) * tv.tauv.nv 
             !     =  2*ux t1*n1 + (uy+vx)*(t1*n2+t2*n1) + 2* t2*n2* vy 
             !     = csf1*ux + csf2*(uy+vx) + csf3*vy
             !     = a21*u(-1) + a22*v(-1) + .... = f2
            ! zsScaled =zs/(1.+10./zs) ! set to zero to turn off "zs" terms
             zsScaled =zs
             csf1= 2.*t1*an1
             csf2=(t1*an2+t2*an1)
             csf3= 2.*t2*an2
             a4(2,1) = csf1*crxd + csf2*cryd  ! coeff of u(-1)
             a4(2,2) = csf2*crxd + csf3*cryd  ! coeff of v(-1)
             a4(2,3) = (zsScaled/mu)*t1       ! coeff of u( 0)
             a4(2,4) = (zsScaled/mu)*t2       ! coeff of v( 0)
             ! Eqn 2: (leave off (zs/mu)*tv.v since it will be otherwise just subtracted below) 
             f2 = csf1*ux + csf2*(uy+vx) + csf3*vy - (1./mu)*( t1*
     & solidT1                   + t2*solidT2                   ) - (
     & zsScaled/mu)*( t1*solidVelocity(i1,i2,i3,0) + t2*solidVelocity(
     & i1,i2,i3,1) )
             if( twilightZone.eq.1 )then
               call ogDeriv(ep,0,0,0,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & t,uc,ue )
               call ogDeriv(ep,0,0,0,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & t,vc,ve )
               call ogDeriv(ep,0,0,0,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & t-dt,uc,une )
               call ogDeriv(ep,0,0,0,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & t-dt,vc,vne )
               call ogDeriv(ep,0,1,0,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & t,uc,uxe )
               call ogDeriv(ep,0,1,0,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & t,vc,vxe )
               call ogDeriv(ep,0,0,1,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & t,uc,uye )
               call ogDeriv(ep,0,0,1,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & t,vc,vye )
               ! Adjust for TZ:  
               f2tz = (csf1*uxe + csf2*(uye+vxe) + csf3*vye) + (
     & zsScaled/mu)*( t1*ue+t2*ve )
               f2 = f2 - f2tz
             end if
             bv(2) = a4(2,1)*u(i1m,i2m,i3m,uc) + a4(2,2)*u(i1m,i2m,i3m,
     & vc) - f2
             ! =========== Equation 3: tangential velocity equation ================
             !          tv.v - (nu*dt)/2*tv.Delta v = tv.rv 
             !          rv = vn - dt*( .5*grad(p) + .5*grad(pn) ) + (nu*dt/2)*(Delta un )
             !  pn : solution at the old time 
             !  (on the predictor stange we assume "u" holds a predicted pressure) 
             nuDtby2 = .5*nu*dt
             ! Delta = (rx^2+ry^2)*urr + (sx^2+sy^2)*uss + (rxx+ryy)*ur + (sxx+syy)*us 
             rxSq = rxd**2 + ryd**2
             sxSq = sxd**2 + syd**2
             cLapg = rxSq/dr(axis)**2 -is*(rxxd+ryyd)/(2.*dr(axis))   ! coeff of ghost in Delta()
             cLap0 = -2.*( rxSq/dr(axis)**2 + sxSq/dr(axisp1)**2 )    ! coeff of bndry in Delta()
             scale=1.  ! Set to zero to turn on interior equation **TEMP**
             a4(3,1) =    -scale*nuDtby2*cLapg*t1   ! coeff of u_{-1}
             a4(3,2) =    -scale*nuDtby2*cLapg*t2   ! coeff of v_{-1}
             a4(3,3) = t1 -scale*nuDtby2*cLap0*t1   ! coeff of u_{0}
             a4(3,4) = t2 -scale*nuDtby2*cLap0*t2   ! coeff of v_{0}
             rv(0) = un(i1,i2,i3,uc) -(.5*dt/rho)*( px + pnx ) + 
     & nuDtby2*( unxx + unyy + uxx+uyy )
             rv(1) = un(i1,i2,i3,vc) -(.5*dt/rho)*( py + pny ) + 
     & nuDtby2*( vnxx + vnyy + vxx+vyy )
             ! adjust RHS for evaluating derivatives with wrong values on the ghost and boundary
             !rv(0) = rv(0) - nuDtby2*( cLapg*u(i1m,i2m,i3m,uc) + cLap0*u(i1,i2,i3,uc) )
             !rv(1) = rv(1) - nuDtby2*( cLapg*u(i1m,i2m,i3m,vc) + cLap0*u(i1,i2,i3,vc) )
             if( twilightZone.eq.1 )then
               ! Eval TZ forcing at t-dt/2 (centered in [t-dt,t]
               thalf=t-.5*dt
               call ogDeriv(ep,1,0,0,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & tHalf,uc,uthe )
               call ogDeriv(ep,1,0,0,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & tHalf,vc,vthe )
               call ogDeriv(ep,0,2,0,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & tHalf,uc,uxxhe )
               call ogDeriv(ep,0,0,2,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & tHalf,uc,uyyhe )
               call ogDeriv(ep,0,2,0,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & tHalf,vc,vxxhe )
               call ogDeriv(ep,0,0,2,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & tHalf,vc,vyyhe )
               call ogDeriv(ep,0,1,0,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & tHalf,pc,pxhe )
               call ogDeriv(ep,0,0,1,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & tHalf,pc,pyhe )
               rv(0) = rv(0) + dt*uthe + (dt/rho)*( pxhe ) - (nu*dt)*( 
     & uxxhe + uyyhe )
               rv(1) = rv(1) + dt*vthe + (dt/rho)*( pyhe ) - (nu*dt)*( 
     & vxxhe + vyyhe )
             end if
             bv(3) = t1*(-u(i1,i2,i3,uc)+rv(0)) + t2*(-u(i1,i2,i3,vc)+
     & rv(1))
             ! adjust RHS for guessed values at boundary and ghost points:
             bv(3) =   a4(3,1)*u(i1m,i2m,i3m,uc) + a4(3,2)*u(i1m,i2m,
     & i3m,vc) + a4(3,3)*u(i1 ,i2 ,i3 ,uc) + a4(3,4)*u(i1 ,i2 ,i3 ,vc)
     &   + bv(3)
             ! =========== Equation 4:  Impedance average of normal velocity equation ================
             !          n.v = alphaMu*n.( rv ) + (1-alphaMu)*n.vSolid + (pcSwitch/(zf+zp))*[ n.SigmaS.n - (-p+n.tauTilde.n) ]
             !          alphaMu = zf/(zf+zp), zf=2*mu*rho/(zp*dt) 
             !          pcSwitch = 0:predictor, 1:corrector
             !          rv = vn - dt*( (1/2)*grad(p) + .5*grad(pn) ) + (nu*dt/2)*(Delta un )
             !   pn : solution at the old time 
             !  (on the predictor stange we assume "u" holds a predicted pressure) 
             ! 
             nDotU = an1*u(i1,i2,i3,uc)+an2*u(i1,i2,i3,vc)
             nDotUs= an1*solidVelocity(i1,i2,i3,0) + an2*solidVelocity(
     & i1,i2,i3,1)
             nTauTilden = mu*( an1*( -2.*vy*an1 + (uy+vx)*an2 ) + an2*(
     &  (uy+vx)*an1 -2.*ux*an2 ) )
             pStar = u(i1,i2,i3,pc) ! current value for p
             scale=1. ! Set to zero to turn off interior equation
             zf = 2.*mu*rho/(zp*dt)
             alphaMu = zf/(zf+zp)
             ! alphaMu=alpha
             beta =  scale*alphaMu*.5*nu*dt
             a4(4,1) =     -beta*cLapg*an1     ! coeff of u_{-1}
             a4(4,2) =     -beta*cLapg*an2     ! coeff of v_{-1}
             a4(4,3) = an1 -beta*cLap0*an1     ! coeff of u_{0}
             a4(4,4) = an2 -beta*cLap0*an2     ! coeff of v_{0}
             ! RHS: 
             bv(4)= -nDotU + scale*alphaMu*( an1*rv(0) + an2*rv(1) ) + 
     & (1.-alphaMu)*nDotUs + (pcSwitch/(zf+zp))*( ( an1*solidT1+an2*
     & solidT2 ) -(-pStar+nTauTilden) )
             ! adjust RHS for guessed values at boundary and ghost points:
             bv(4) =   a4(4,1)*u(i1m,i2m,i3m,uc) + a4(4,2)*u(i1m,i2m,
     & i3m,vc) + a4(4,3)*u(i1 ,i2 ,i3 ,uc) + a4(4,4)*u(i1 ,i2 ,i3 ,vc)
     &   + bv(4)
             if( twilightZone.eq.1 )then
               ! write(*,'(" an1*rv(0) + an2*rv(1) - an1*ue + an2*ve=",1pe12.3)') an1*rv(0) + an2*rv(1)- (an1*ue + an2*ve)
               ! Adjust for TZ:  
               ! With TZ, the RHS above should equal (-1+alpha)*( n.ve ) so we need to add (1-alpha)*n.ve
               !   n.v = alpha*n.ve + (1-alpha)n.ve = n.ve 
               nDotUe = an1*ue+an2*ve
               nTauTildene = mu*( an1*( -2*vye*an1 + (uye+vxe)*an2) + 
     & an2*( (uye+vxe)*an1 -2.*uxe*an2 ) )
               f4tz = (1.-scale*alphaMu)*nDotUe + zp*nDotUe + (
     & pcSwitch/(zf+zp))*(-pe + nTauTildene )
               bv(4) = bv(4) + f4tz
             end if
             ! write(*,'(" i1,i2=",2i3," rxd,ryd=",2f8.4," sxd,syd=",2f10.4)') i1,i2,rxd,ryd,sxd,syd
             ! write(*,'(" i1,i2=",2i3," a11,a12=",2f8.4," a21,a22=",2f10.4)') i1,i2,a11,a12,a21,a22
             ! write(*,'(" i1,i2=",2i3," solidTraction=",2e10.2)') i1,i2,solidTraction(i1,i2,i3,0),solidTraction(i1,i2,i3,1)
             ! Solve
             !   [a11 a12 a13 a14 ][ u(-1)] = [ b1 ]
             !   [a21 a22 a23 a24 ][ v(-1)] = [ b2 ]
             !   [a31 a32 a33 a34 ][ u( 0)] = [ b3 ]
             !   [a41 a42 a43 a44 ][ v( 0)] = [ b4 ]
             job=0
             numberOfEquations=4
             call dgeco( a4(1,1), numberOfEquations, numberOfEquations,
     &  ipvt(1),rcond,work(1))
             call dgesl( a4(1,1), numberOfEquations, numberOfEquations,
     &  ipvt(1), bv(1), job)
             ! --------------------- CHECK RESULTS -------------------
             if( checkResults.eq.1 )then
               ! check that the equations are all satisfied
               write(*,'("AMP EXPLICIT VELOCITY BC: i1,i2=",2i3," zs=",
     & e16.8," alpha=",e16.8," beta=",e10.2," rcond=",e10.2)') i1,i2,
     & zs,alpha,beta,rcond
               write(*,'(" normal=",2e14.4," tangent=",2e14.4)') an1,
     & an2,t1,t2
               write(*,'(" px,py=",2e14.4," old: px,py=",2e14.4)') px,
     & py,pnx,pny
               write(*,'(" vSolid=",2e16.6," solidTraction=",2e16.6)') 
     & solidVelocity(i1,i2,i3,0),solidVelocity(i1,i2,i3,1),
     & solidTraction(i1,i2,i3,0),solidTraction(i1,i2,i3,1)
               write(*,'(" old: uxx,uyy,vxx,vyy==",4e16.6)') unxx,unyy,
     & vnxx,vnyy
               ! save current values
               us(1)=u(i1m,i2m,i3m,uc)
               us(2)=u(i1m,i2m,i3m,vc)
               us(3)=u(i1 ,i2 ,i3 ,uc)
               us(4)=u(i1 ,i2 ,i3 ,vc)
               u(i1m,i2m,i3m,uc)=bv(1)
               u(i1m,i2m,i3m,vc)=bv(2)
               u(i1 ,i2 ,i3 ,uc)=bv(3)
               u(i1 ,i2 ,i3 ,vc)=bv(4)
               write(*,'(" Old: u,v=",2e16.6," ghost=",2e16.6)') un(i1 
     & ,i2 ,i3 ,uc),un(i1 ,i2 ,i3 ,vc),un(i1m,i2m,i3m,uc),un(i1m,i2m,
     & i3m,vc)
               write(*,'(" New: u,v=",2e16.6," ghost=",2e16.6)') u(i1 ,
     & i2 ,i3 ,uc),u(i1 ,i2 ,i3 ,vc),u(i1m,i2m,i3m,uc),u(i1m,i2m,i3m,
     & vc)
               ux = ux22(i1,i2,i3,uc)
               uy = uy22(i1,i2,i3,uc)
               vx = ux22(i1,i2,i3,vc)
               vy = uy22(i1,i2,i3,vc)
               uxx = uxx22(i1,i2,i3,uc)
               uyy = uyy22(i1,i2,i3,uc)
               vxx = uxx22(i1,i2,i3,vc)
               vyy = uyy22(i1,i2,i3,vc)
               ! =======================   EQN 1 ==============================================
               div=ux+vy
               maxDiv=max(maxDiv,abs(div))
               write(*,'("EQN1: After: divergence=",e10.2)') div
               ! =======================   EQN 2 ==============================================
               ! (1/mu)* tv.sigma.n + (zs/mu)*tv.v = (1/mu)*tv.solidTraction + (zs/mu)*tv.solidVelocity
               tDotTaun =  mu*( t1*(2.*ux*an1+(uy+vx)*an2) + t2*((uy+
     & vx)*an1+2.*vy*an2) )
               char = (1./mu)*(  tDotTaun + zs*(t1*u(i1 ,i2 ,i3 ,uc)+
     & t2*u(i1 ,i2 ,i3 ,vc)) )
               f2 =  (1./mu)*( t1*solidT1 + t2*solidT2 ) + (zs/mu)*( 
     & t1*solidVelocity(i1,i2,i3,0) + t2*solidVelocity(i1,i2,i3,1) )
               ! Adjust for TZ ? June 25, 2018
               if( twilightZone.eq.1 )then
                ! Adjust for TZ:  
                f2 = f2 + f2tz
               end if
               res = abs(char-f2)
               maxChar=max(maxChar,res)
               write(*,'("EQN2: After: ( t.tau.n + zs*t.v ) ...=",
     & e10.2,"  (char,f2=",2e10.2,")")') res,char,f2
               ! =======================   EQN 3 ==============================================
               ! Here is the solution on the boundary from the interior equations
               rv(0) = un(i1,i2,i3,uc) -(.5*dt)*( px + pnx ) + nuDtby2*
     & ( unxx + unyy + uxx+uyy )
               rv(1) = un(i1,i2,i3,vc) -(.5*dt)*( py + pny ) + nuDtby2*
     & ( vnxx + vnyy + vxx+vyy )
               if( twilightZone.eq.1 )then
                 rv(0) = rv(0) + dt*uthe + (dt/rho)*( pxhe ) - (nu*dt)*
     & ( uxxhe + uyyhe )
                 rv(1) = rv(1) + dt*vthe + (dt/rho)*( pyhe ) - (nu*dt)*
     & ( vxxhe + vyyhe )
               end if
               res = abs( t1*rv(0) + t2*rv(1) - (t1*u(i1 ,i2 ,i3 ,uc) +
     &  t2*u(i1 ,i2 ,i3 ,vc) ) )
               maxRes1=max(maxRes1,res)
               write(*,'("EQN3: After: res in tv.( EQN) ) ...=",e10.2)
     & ') res
               ! =========== Equation 4:  Impedance average of normal velocity equation ================
               !          n.v = alphaMu*n.( rv ) + (1-alphaMu)*n.vSolid + (pcSwitch/(zf+zp))*[ n.SigmaS.n - (-p+n.tauTilde.n) ]
               nDotU = an1*u(i1,i2,i3,uc)+an2*u(i1,i2,i3,vc)
               nDotUs= an1*solidVelocity(i1,i2,i3,0) + an2*
     & solidVelocity(i1,i2,i3,1)
               nTauTilden = mu*( an1*( -2.*vy*an1 + (uy+vx)*an2 ) + 
     & an2*( (uy+vx)*an1 -2.*ux*an2 ) )
               res =  -nDotU + scale*alphaMu*( an1*rv(0) + an2*rv(1) ) 
     & + (1.-alphaMu)*nDotUs + (pcSwitch/(zf+zp))*( ( an1*solidT1+an2*
     & solidT2 ) -(-pStar+nTauTilden) )
               if( twilightZone.eq.1 )then
                 ! Adjust for TZ:  
                 res = res  + f4tz
               end if
               maxRes2=max(maxRes2,res)
               write(*,'("EQN4: After: res in nv.( EQN) ) ...=",e10.2)
     & ') res
               ! reset
               u(i1m,i2m,i3m,uc)=us(1)
               u(i1m,i2m,i3m,vc)=us(2)
               u(i1 ,i2 ,i3 ,uc)=us(3)
               u(i1 ,i2 ,i3 ,vc)=us(4)
             end if
             if( .false. .and. knownSolution.ne.0 )then
              body=0
              stateChoice=boundaryVelocity
              n=0
              call evalUserDefinedDeformingBodyKnownSolution(  body,
     & stateChoice,t,grid,i1m,i2m,i3m,n,val(0) )
              call evalUserDefinedDeformingBodyKnownSolution(  body,
     & stateChoice,t,grid,i1 ,i2 ,i3 ,n,val(2) )
              write(*,'(" i1,i2=",2i3," alpha=",e10.2," beta=",e10.2," 
     & rcond=",e10.2)') i1,i2,alpha,beta,rcond
              write(*,'("  insExpAmpBC: (knownSolution): v-err=",4(
     & 1pe10.2)," |u|,|v|=",2e10.2)') bv(1)-val(0),bv(2)-val(1),bv(3)-
     & val(2),bv(4)-val(3),abs(val(0)),abs(val(1))
              ! set exact 
              if( .true. )then
                write(*,'("  insExpAmpBC: **TEMP** set EXACT velocity 
     & on the interface")')
                bv(1)=val(0)
                bv(2)=val(1)
                bv(3)=val(2)
                bv(4)=val(3)
              end if
             end if
             if( .FALSE. )then
             write(*,'(" i1,i2=",2i3," alpha=",e10.2," beta=",e10.2," 
     & rcond=",e10.2)') i1,i2,alpha,beta,rcond
             write(*,'(" bv=",4e10.2)') bv(1),bv(2),bv(3),bv(4)
             if( twilightZone.eq.1 )then
               call ogDeriv(ep,0,2,0,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & t,uc,uxxe )
               call ogDeriv(ep,0,0,2,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & t,uc,uyye )
               call ogDeriv(ep,0,2,0,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & t,vc,vxxe )
               call ogDeriv(ep,0,0,2,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & t,vc,vyye )
               write(*,'(" uxx-uxxe,uyy-uyye,vxx-vxxe,vyy-vyye=",
     & 4e10.2)')  uxx-uxxe,uyy-uyye,vxx-vxxe,vyy-vyye
               call ogDeriv(ep,0,2,0,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & t-dt,uc,uxxe )
               call ogDeriv(ep,0,0,2,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & t-dt,uc,uyye )
               call ogDeriv(ep,0,2,0,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & t-dt,vc,vxxe )
               call ogDeriv(ep,0,0,2,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & t-dt,vc,vyye )
               write(*,'(" unxx-unxxe,unyy-unyye,vnxx-vnxxe,vnyy-
     & vnyye=",4e10.2)')  unxx-uxxe,unyy-uyye,vnxx-vxxe,vnyy-vyye
               write(*,'(" rxxd,ryyd=",2e10.2)') rxxd,ryyd
               call ogDeriv(ep,0,1,0,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & t,pc,pxe )
               call ogDeriv(ep,0,0,1,0,x(i1,i2,i3,0),x(i1,i2,i3,1),0.,
     & t,pc,pye )
               write(*,'(" px-pxe,py-pye=",2e10.2)') px-pxe, py-pye
               write(*,'(" Old time: un-une=",e10.2," vn-vne=",e10.2)')
     &  un(i1,i2,i3,uc)-une,un(i1,i2,i3,vc)-vne
             end if
             write(*,'(" solidVelocity=",2e10.2)') solidVelocity(i1,i2,
     & i3,0),solidVelocity(i1,i2,i3,0)
             write(*,'(" -> predicted u(-1),v(-1),u(0),v(0)=",4e12.4)')
     &  u(i1m,i2m,i3m,uc),u(i1m,i2m,i3m,vc),u(i1 ,i2 ,i3 ,uc),u(i1 ,
     & i2 ,i3 ,vc)
             write(*,'(" -> new       u(-1),v(-1),u(0),v(0)=",4e12.4)')
     &  bv(1),bv(2),bv(3),bv(4)
             write(*,'(" -> DIFF      u(-1),v(-1),u(0),v(0)=",4e12.4)')
     &   u(i1m,i2m,i3m,uc)-bv(1),u(i1m,i2m,i3m,vc)-bv(2),u(i1 ,i2 ,i3 
     & ,uc)-bv(3),u(i1 ,i2 ,i3 ,vc)-bv(4)
             end if
             if( useJacobiUpdate.eq.1 )then
               uNew(i1m,i2m,i3m,uc)=bv(1)
               uNew(i1m,i2m,i3m,vc)=bv(2)
               uNew(i1 ,i2 ,i3 ,uc)=bv(3)
               uNew(i1 ,i2 ,i3 ,vc)=bv(4)
             else
               u(i1m,i2m,i3m,uc)=bv(1)
               u(i1m,i2m,i3m,vc)=bv(2)
               u(i1 ,i2 ,i3 ,uc)=bv(3)
               u(i1 ,i2 ,i3 ,vc)=bv(4)
             end if
            end do
            end do
            end do
            ! Now copy the copy into u 
            if( useJacobiUpdate.eq.1 )then
             do i3=n3a,n3b
             do i2=n2a,n2b
             do i1=n1a,n1b
              i1m=i1-is1
              i2m=i2-is2
              i3m=i3-is3
              u(i1m,i2m,i3m,uc)=uNew(i1m,i2m,i3m,uc)
              u(i1m,i2m,i3m,vc)=uNew(i1m,i2m,i3m,vc)
              u(i1 ,i2 ,i3 ,uc)=uNew(i1 ,i2 ,i3 ,uc)
              u(i1 ,i2 ,i3 ,vc)=uNew(i1 ,i2 ,i3 ,vc)
             end do
             end do
             end do
            end if
            if( checkResults.eq.1 )then
              write(*,'(" maxDiv=",e10.2," maxChar=",e10.2," maxRes1=",
     & e10.2," maxRes2=",e10.2)') maxDiv,maxChar,maxRes1,maxRes2
              stop 1111
            end if
          else
            write(*,'("insExplicitAMPVelocityBC: finish me")')
            stop 2222
            ! tractionCurvilinear2dOrder2()
          end if
         else
           stop 4455
         end if

        else

          stop 1111
        end if

       else
         ! ---- 3D TRACTION ----
         write(*,'("insExplicitAMPVelocityBC: finish me")')
         stop 3333
        if( gridType.eq.0 )then
          ! --- RECTANGULAR ----
          if( orderOfAccuracy.eq.2 )then
            if( axis.eq.0 )then
              ! tractionRectangular(2,x,3)

            else if( axis.eq.1 )then
              write(*,'("CALL TRACTION FREE Y")')
              ! tractionRectangular(2,y,3)

            else if( axis.eq.2 )then
               ! finish me for solid traction
               stop 4321
              ! write(*,'("START TRACTION LOOPS RECTANGULAR ")') 
              f1=0.
              f2=0.
              f3=0.
              do i3=n3a,n3b
              do i2=n2a,n2b
              do i1=n1a,n1b
               i1m=i1-is1
               i2m=i2-is2
               i3m=i3-is3
               i1p=i1+is1
               i2p=i2+is2
               i3p=i3+is3
                ! --------------------------------------------------------------
                ! ----------------- 3D Traction Rectangular --------------------
                ! --------------------------------------------------------------
                 ! uz =0 
                 ! vz = 0
                 ! wz = -(ux+vy)
                 if( twilightZone.eq.1 )then
                   ! assume the TZ solution is divergence free
                   call ogDeriv(ep,0,0,0,1,x(i1,i2,i3,0),x(i1,i2,i3,1),
     & x(i1,i2,i3,2),t,uc,uze )
                   call ogDeriv(ep,0,0,0,1,x(i1,i2,i3,0),x(i1,i2,i3,1),
     & x(i1,i2,i3,2),t,vc,vze )
                   f1 = -is*2.*dx(axis)*uze
                   f2 = -is*2.*dx(axis)*vze
                 end if
                 u(i1m,i2m,i3m,uc)= u(i1p,i2p,i3p,uc) + f1
                 u(i1m,i2m,i3m,vc)= u(i1p,i2p,i3p,vc) + f2
                 u(i1m,i2m,i3m,wc)= u(i1p,i2p,i3p,wc) + is*2.*dx(axis)*
     & (ux23r(i1,i2,i3,uc)+uy23r(i1,i2,i3,vc))
              end do
              end do
              end do
            else
              stop 4444
            end if
          else if( orderOfAccuracy.eq.4 )then
            stop 3333
          else
           stop 1234
          end if

        else if( gridType.eq.1 )then
          ! --- CURVILINEAR ----

          ! tractionCurvilinear3dOrder2()

        else
          ! unknown gridType 
          stop 2222
        end if
       end if


        end if ! end if noSlipWall
      !! end do ! end side
      !! end do ! end axis
      ! ================= END LOOP OVER SIDES ===============================


      return
      end


