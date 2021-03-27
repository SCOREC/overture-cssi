! This file automatically generated from getDiv.bf with bpp.
       subroutine getDiv3dOrder2c(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,
     & nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,rsxy,  u,p,v, divD, divB, 
     & K0, matMask, ipar, rpar, ierr )
      !======================================================================
      !   Compute the divergence 
      !
      !   v : tempoary storage to hold D and B for the BA Maxwell
      !   divD(i1,i2,i3,0:1) : return div(E) ( or div(D) for BA) here 
      !   divB(i1,i2,i3,0:1) : return div(H) ( or Div(B) for BA) here 
      !======================================================================
       implicit none
       integer nd, n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,
     & nd3b,nd4a,nd4b
       real u(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
       real v(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
       real divD(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
       real divB(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
       ! Polarization vectors
       real p(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
       real rsxy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1,0:nd-1)
       integer mask(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
       real K0(0:5,0:5,0:*)  ! material matrix
       integer matMask(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
       integer ierr
       integer ipar(0:*)
       real rpar(0:*)
       !     ---- local variables -----
       integer c,i1,i2,i3,n,gridType,orderOfAccuracy,debug
       integer m1a,m1b,m2a,m2b,m3a,m3b,numGhost,k1,k2
       integer saveDivergence
       real t,dt
       real d1x,d2x,d3x, d1y,d2y,d3y, d1z,d2z,d3z
       real b1x,b2x,b3x, b1y,b2y,b3y, b1z,b2z,b3z
       real dDiv, bDiv
       real divDmax,divBmax, gradDmax, gradBmax
       integer ex,ey,ez, hx,hy,hz, solveForAllFields, grid
       integer maxRegions,NpMax
       parameter( maxRegions=100,NpMax=10 )  ! FIX ME
       integer numberOfMaterialRegions, mr
       real K03(0:2,0:2,0:maxRegions) ! 3x3 material matrix for TEz polarization
       integer Np(6,6,0:maxRegions-1)
       real gdmPar(4,NpMax,6,6,0:maxRegions-1), ptSum(0:5)
       integer ec,pc,qc, pct,qct
       integer numPolarizationTerms
       integer maxNumPolarizationTerms
       parameter( maxNumPolarizationTerms=200 )
       real dx(0:2),dr(0:2)
       integer rectangular,curvilinear
       parameter( rectangular=0, curvilinear=1 )
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
       ! dispersion
       integer dispersionModel
       integer numTerms1(0:maxRegions),ecIndex1(
     & maxNumPolarizationTerms,0:maxRegions),pcIndex1(
     & maxNumPolarizationTerms,0:maxRegions)
      ! Dispersion models
       integer noDispersion,drude,gdm
       parameter( noDispersion=0, drude=1, gdm=2 )
      !  integer forcingOption   
      !  ! forcing options
      !  #Include "forcingDefineFortranInclude.h"
       integer method,nfdtd,bamx
       parameter( nfdtd=5, bamx=7 )
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
       vr2(i1,i2,i3,kd)=(v(i1+1,i2,i3,kd)-v(i1-1,i2,i3,kd))*d12(0)
       vs2(i1,i2,i3,kd)=(v(i1,i2+1,i3,kd)-v(i1,i2-1,i3,kd))*d12(1)
       vt2(i1,i2,i3,kd)=(v(i1,i2,i3+1,kd)-v(i1,i2,i3-1,kd))*d12(2)
       vrr2(i1,i2,i3,kd)=(-2.*v(i1,i2,i3,kd)+(v(i1+1,i2,i3,kd)+v(i1-1,
     & i2,i3,kd)) )*d22(0)
       vss2(i1,i2,i3,kd)=(-2.*v(i1,i2,i3,kd)+(v(i1,i2+1,i3,kd)+v(i1,i2-
     & 1,i3,kd)) )*d22(1)
       vrs2(i1,i2,i3,kd)=(vr2(i1,i2+1,i3,kd)-vr2(i1,i2-1,i3,kd))*d12(1)
       vtt2(i1,i2,i3,kd)=(-2.*v(i1,i2,i3,kd)+(v(i1,i2,i3+1,kd)+v(i1,i2,
     & i3-1,kd)) )*d22(2)
       vrt2(i1,i2,i3,kd)=(vr2(i1,i2,i3+1,kd)-vr2(i1,i2,i3-1,kd))*d12(2)
       vst2(i1,i2,i3,kd)=(vs2(i1,i2,i3+1,kd)-vs2(i1,i2,i3-1,kd))*d12(2)
       vrrr2(i1,i2,i3,kd)=(-2.*(v(i1+1,i2,i3,kd)-v(i1-1,i2,i3,kd))+(v(
     & i1+2,i2,i3,kd)-v(i1-2,i2,i3,kd)) )*d22(0)*d12(0)
       vsss2(i1,i2,i3,kd)=(-2.*(v(i1,i2+1,i3,kd)-v(i1,i2-1,i3,kd))+(v(
     & i1,i2+2,i3,kd)-v(i1,i2-2,i3,kd)) )*d22(1)*d12(1)
       vttt2(i1,i2,i3,kd)=(-2.*(v(i1,i2,i3+1,kd)-v(i1,i2,i3-1,kd))+(v(
     & i1,i2,i3+2,kd)-v(i1,i2,i3-2,kd)) )*d22(2)*d12(2)
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
       vx21(i1,i2,i3,kd)= rx(i1,i2,i3)*vr2(i1,i2,i3,kd)
       vy21(i1,i2,i3,kd)=0
       vz21(i1,i2,i3,kd)=0
       vx22(i1,i2,i3,kd)= rx(i1,i2,i3)*vr2(i1,i2,i3,kd)+sx(i1,i2,i3)*
     & vs2(i1,i2,i3,kd)
       vy22(i1,i2,i3,kd)= ry(i1,i2,i3)*vr2(i1,i2,i3,kd)+sy(i1,i2,i3)*
     & vs2(i1,i2,i3,kd)
       vz22(i1,i2,i3,kd)=0
       vx23(i1,i2,i3,kd)=rx(i1,i2,i3)*vr2(i1,i2,i3,kd)+sx(i1,i2,i3)*
     & vs2(i1,i2,i3,kd)+tx(i1,i2,i3)*vt2(i1,i2,i3,kd)
       vy23(i1,i2,i3,kd)=ry(i1,i2,i3)*vr2(i1,i2,i3,kd)+sy(i1,i2,i3)*
     & vs2(i1,i2,i3,kd)+ty(i1,i2,i3)*vt2(i1,i2,i3,kd)
       vz23(i1,i2,i3,kd)=rz(i1,i2,i3)*vr2(i1,i2,i3,kd)+sz(i1,i2,i3)*
     & vs2(i1,i2,i3,kd)+tz(i1,i2,i3)*vt2(i1,i2,i3,kd)
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
       vxx21(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*vrr2(i1,i2,i3,kd)+(rxx22(
     & i1,i2,i3))*vr2(i1,i2,i3,kd)
       vyy21(i1,i2,i3,kd)=0
       vxy21(i1,i2,i3,kd)=0
       vxz21(i1,i2,i3,kd)=0
       vyz21(i1,i2,i3,kd)=0
       vzz21(i1,i2,i3,kd)=0
       vlaplacian21(i1,i2,i3,kd)=vxx21(i1,i2,i3,kd)
       vxx22(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*vrr2(i1,i2,i3,kd)+2.*(rx(
     & i1,i2,i3)*sx(i1,i2,i3))*vrs2(i1,i2,i3,kd)+(sx(i1,i2,i3)**2)*
     & vss2(i1,i2,i3,kd)+(rxx22(i1,i2,i3))*vr2(i1,i2,i3,kd)+(sxx22(i1,
     & i2,i3))*vs2(i1,i2,i3,kd)
       vyy22(i1,i2,i3,kd)=(ry(i1,i2,i3)**2)*vrr2(i1,i2,i3,kd)+2.*(ry(
     & i1,i2,i3)*sy(i1,i2,i3))*vrs2(i1,i2,i3,kd)+(sy(i1,i2,i3)**2)*
     & vss2(i1,i2,i3,kd)+(ryy22(i1,i2,i3))*vr2(i1,i2,i3,kd)+(syy22(i1,
     & i2,i3))*vs2(i1,i2,i3,kd)
       vxy22(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*vrr2(i1,i2,i3,kd)+(
     & rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*vrs2(i1,
     & i2,i3,kd)+sx(i1,i2,i3)*sy(i1,i2,i3)*vss2(i1,i2,i3,kd)+rxy22(i1,
     & i2,i3)*vr2(i1,i2,i3,kd)+sxy22(i1,i2,i3)*vs2(i1,i2,i3,kd)
       vxz22(i1,i2,i3,kd)=0
       vyz22(i1,i2,i3,kd)=0
       vzz22(i1,i2,i3,kd)=0
       vlaplacian22(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2)*
     & vrr2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*
     & sy(i1,i2,i3))*vrs2(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**
     & 2)*vss2(i1,i2,i3,kd)+(rxx22(i1,i2,i3)+ryy22(i1,i2,i3))*vr2(i1,
     & i2,i3,kd)+(sxx22(i1,i2,i3)+syy22(i1,i2,i3))*vs2(i1,i2,i3,kd)
       vxx23(i1,i2,i3,kd)=rx(i1,i2,i3)**2*vrr2(i1,i2,i3,kd)+sx(i1,i2,
     & i3)**2*vss2(i1,i2,i3,kd)+tx(i1,i2,i3)**2*vtt2(i1,i2,i3,kd)+2.*
     & rx(i1,i2,i3)*sx(i1,i2,i3)*vrs2(i1,i2,i3,kd)+2.*rx(i1,i2,i3)*tx(
     & i1,i2,i3)*vrt2(i1,i2,i3,kd)+2.*sx(i1,i2,i3)*tx(i1,i2,i3)*vst2(
     & i1,i2,i3,kd)+rxx23(i1,i2,i3)*vr2(i1,i2,i3,kd)+sxx23(i1,i2,i3)*
     & vs2(i1,i2,i3,kd)+txx23(i1,i2,i3)*vt2(i1,i2,i3,kd)
       vyy23(i1,i2,i3,kd)=ry(i1,i2,i3)**2*vrr2(i1,i2,i3,kd)+sy(i1,i2,
     & i3)**2*vss2(i1,i2,i3,kd)+ty(i1,i2,i3)**2*vtt2(i1,i2,i3,kd)+2.*
     & ry(i1,i2,i3)*sy(i1,i2,i3)*vrs2(i1,i2,i3,kd)+2.*ry(i1,i2,i3)*ty(
     & i1,i2,i3)*vrt2(i1,i2,i3,kd)+2.*sy(i1,i2,i3)*ty(i1,i2,i3)*vst2(
     & i1,i2,i3,kd)+ryy23(i1,i2,i3)*vr2(i1,i2,i3,kd)+syy23(i1,i2,i3)*
     & vs2(i1,i2,i3,kd)+tyy23(i1,i2,i3)*vt2(i1,i2,i3,kd)
       vzz23(i1,i2,i3,kd)=rz(i1,i2,i3)**2*vrr2(i1,i2,i3,kd)+sz(i1,i2,
     & i3)**2*vss2(i1,i2,i3,kd)+tz(i1,i2,i3)**2*vtt2(i1,i2,i3,kd)+2.*
     & rz(i1,i2,i3)*sz(i1,i2,i3)*vrs2(i1,i2,i3,kd)+2.*rz(i1,i2,i3)*tz(
     & i1,i2,i3)*vrt2(i1,i2,i3,kd)+2.*sz(i1,i2,i3)*tz(i1,i2,i3)*vst2(
     & i1,i2,i3,kd)+rzz23(i1,i2,i3)*vr2(i1,i2,i3,kd)+szz23(i1,i2,i3)*
     & vs2(i1,i2,i3,kd)+tzz23(i1,i2,i3)*vt2(i1,i2,i3,kd)
       vxy23(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*vrr2(i1,i2,i3,kd)+
     & sx(i1,i2,i3)*sy(i1,i2,i3)*vss2(i1,i2,i3,kd)+tx(i1,i2,i3)*ty(i1,
     & i2,i3)*vtt2(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,
     & i3)*sx(i1,i2,i3))*vrs2(i1,i2,i3,kd)+(rx(i1,i2,i3)*ty(i1,i2,i3)+
     & ry(i1,i2,i3)*tx(i1,i2,i3))*vrt2(i1,i2,i3,kd)+(sx(i1,i2,i3)*ty(
     & i1,i2,i3)+sy(i1,i2,i3)*tx(i1,i2,i3))*vst2(i1,i2,i3,kd)+rxy23(
     & i1,i2,i3)*vr2(i1,i2,i3,kd)+sxy23(i1,i2,i3)*vs2(i1,i2,i3,kd)+
     & txy23(i1,i2,i3)*vt2(i1,i2,i3,kd)
       vxz23(i1,i2,i3,kd)=rx(i1,i2,i3)*rz(i1,i2,i3)*vrr2(i1,i2,i3,kd)+
     & sx(i1,i2,i3)*sz(i1,i2,i3)*vss2(i1,i2,i3,kd)+tx(i1,i2,i3)*tz(i1,
     & i2,i3)*vtt2(i1,i2,i3,kd)+(rx(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,i2,
     & i3)*sx(i1,i2,i3))*vrs2(i1,i2,i3,kd)+(rx(i1,i2,i3)*tz(i1,i2,i3)+
     & rz(i1,i2,i3)*tx(i1,i2,i3))*vrt2(i1,i2,i3,kd)+(sx(i1,i2,i3)*tz(
     & i1,i2,i3)+sz(i1,i2,i3)*tx(i1,i2,i3))*vst2(i1,i2,i3,kd)+rxz23(
     & i1,i2,i3)*vr2(i1,i2,i3,kd)+sxz23(i1,i2,i3)*vs2(i1,i2,i3,kd)+
     & txz23(i1,i2,i3)*vt2(i1,i2,i3,kd)
       vyz23(i1,i2,i3,kd)=ry(i1,i2,i3)*rz(i1,i2,i3)*vrr2(i1,i2,i3,kd)+
     & sy(i1,i2,i3)*sz(i1,i2,i3)*vss2(i1,i2,i3,kd)+ty(i1,i2,i3)*tz(i1,
     & i2,i3)*vtt2(i1,i2,i3,kd)+(ry(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,i2,
     & i3)*sy(i1,i2,i3))*vrs2(i1,i2,i3,kd)+(ry(i1,i2,i3)*tz(i1,i2,i3)+
     & rz(i1,i2,i3)*ty(i1,i2,i3))*vrt2(i1,i2,i3,kd)+(sy(i1,i2,i3)*tz(
     & i1,i2,i3)+sz(i1,i2,i3)*ty(i1,i2,i3))*vst2(i1,i2,i3,kd)+ryz23(
     & i1,i2,i3)*vr2(i1,i2,i3,kd)+syz23(i1,i2,i3)*vs2(i1,i2,i3,kd)+
     & tyz23(i1,i2,i3)*vt2(i1,i2,i3,kd)
       vlaplacian23(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2+rz(
     & i1,i2,i3)**2)*vrr2(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**
     & 2+sz(i1,i2,i3)**2)*vss2(i1,i2,i3,kd)+(tx(i1,i2,i3)**2+ty(i1,i2,
     & i3)**2+tz(i1,i2,i3)**2)*vtt2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(
     & i1,i2,i3)+ ry(i1,i2,i3)*sy(i1,i2,i3)+rz(i1,i2,i3)*sz(i1,i2,i3))
     & *vrs2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*tx(i1,i2,i3)+ ry(i1,i2,i3)*
     & ty(i1,i2,i3)+rz(i1,i2,i3)*tz(i1,i2,i3))*vrt2(i1,i2,i3,kd)+2.*(
     & sx(i1,i2,i3)*tx(i1,i2,i3)+ sy(i1,i2,i3)*ty(i1,i2,i3)+sz(i1,i2,
     & i3)*tz(i1,i2,i3))*vst2(i1,i2,i3,kd)+(rxx23(i1,i2,i3)+ryy23(i1,
     & i2,i3)+rzz23(i1,i2,i3))*vr2(i1,i2,i3,kd)+(sxx23(i1,i2,i3)+
     & syy23(i1,i2,i3)+szz23(i1,i2,i3))*vs2(i1,i2,i3,kd)+(txx23(i1,i2,
     & i3)+tyy23(i1,i2,i3)+tzz23(i1,i2,i3))*vt2(i1,i2,i3,kd)
c============================================================================================
c Define derivatives for a rectangular grid
c
c============================================================================================
       h12(kd) = 1./(2.*dx(kd))
       h22(kd) = 1./(dx(kd)**2)
       vx23r(i1,i2,i3,kd)=(v(i1+1,i2,i3,kd)-v(i1-1,i2,i3,kd))*h12(0)
       vy23r(i1,i2,i3,kd)=(v(i1,i2+1,i3,kd)-v(i1,i2-1,i3,kd))*h12(1)
       vz23r(i1,i2,i3,kd)=(v(i1,i2,i3+1,kd)-v(i1,i2,i3-1,kd))*h12(2)
       vxx23r(i1,i2,i3,kd)=(-2.*v(i1,i2,i3,kd)+(v(i1+1,i2,i3,kd)+v(i1-
     & 1,i2,i3,kd)) )*h22(0)
       vyy23r(i1,i2,i3,kd)=(-2.*v(i1,i2,i3,kd)+(v(i1,i2+1,i3,kd)+v(i1,
     & i2-1,i3,kd)) )*h22(1)
       vxy23r(i1,i2,i3,kd)=(vx23r(i1,i2+1,i3,kd)-vx23r(i1,i2-1,i3,kd))*
     & h12(1)
       vzz23r(i1,i2,i3,kd)=(-2.*v(i1,i2,i3,kd)+(v(i1,i2,i3+1,kd)+v(i1,
     & i2,i3-1,kd)) )*h22(2)
       vxz23r(i1,i2,i3,kd)=(vx23r(i1,i2,i3+1,kd)-vx23r(i1,i2,i3-1,kd))*
     & h12(2)
       vyz23r(i1,i2,i3,kd)=(vy23r(i1,i2,i3+1,kd)-vy23r(i1,i2,i3-1,kd))*
     & h12(2)
       vx21r(i1,i2,i3,kd)= vx23r(i1,i2,i3,kd)
       vy21r(i1,i2,i3,kd)= vy23r(i1,i2,i3,kd)
       vz21r(i1,i2,i3,kd)= vz23r(i1,i2,i3,kd)
       vxx21r(i1,i2,i3,kd)= vxx23r(i1,i2,i3,kd)
       vyy21r(i1,i2,i3,kd)= vyy23r(i1,i2,i3,kd)
       vzz21r(i1,i2,i3,kd)= vzz23r(i1,i2,i3,kd)
       vxy21r(i1,i2,i3,kd)= vxy23r(i1,i2,i3,kd)
       vxz21r(i1,i2,i3,kd)= vxz23r(i1,i2,i3,kd)
       vyz21r(i1,i2,i3,kd)= vyz23r(i1,i2,i3,kd)
       vlaplacian21r(i1,i2,i3,kd)=vxx23r(i1,i2,i3,kd)
       vx22r(i1,i2,i3,kd)= vx23r(i1,i2,i3,kd)
       vy22r(i1,i2,i3,kd)= vy23r(i1,i2,i3,kd)
       vz22r(i1,i2,i3,kd)= vz23r(i1,i2,i3,kd)
       vxx22r(i1,i2,i3,kd)= vxx23r(i1,i2,i3,kd)
       vyy22r(i1,i2,i3,kd)= vyy23r(i1,i2,i3,kd)
       vzz22r(i1,i2,i3,kd)= vzz23r(i1,i2,i3,kd)
       vxy22r(i1,i2,i3,kd)= vxy23r(i1,i2,i3,kd)
       vxz22r(i1,i2,i3,kd)= vxz23r(i1,i2,i3,kd)
       vyz22r(i1,i2,i3,kd)= vyz23r(i1,i2,i3,kd)
       vlaplacian22r(i1,i2,i3,kd)=vxx23r(i1,i2,i3,kd)+vyy23r(i1,i2,i3,
     & kd)
       vlaplacian23r(i1,i2,i3,kd)=vxx23r(i1,i2,i3,kd)+vyy23r(i1,i2,i3,
     & kd)+vzz23r(i1,i2,i3,kd)
       vxxx22r(i1,i2,i3,kd)=(-2.*(v(i1+1,i2,i3,kd)-v(i1-1,i2,i3,kd))+(
     & v(i1+2,i2,i3,kd)-v(i1-2,i2,i3,kd)) )*h22(0)*h12(0)
       vyyy22r(i1,i2,i3,kd)=(-2.*(v(i1,i2+1,i3,kd)-v(i1,i2-1,i3,kd))+(
     & v(i1,i2+2,i3,kd)-v(i1,i2-2,i3,kd)) )*h22(1)*h12(1)
       vxxy22r(i1,i2,i3,kd)=( vxx22r(i1,i2+1,i3,kd)-vxx22r(i1,i2-1,i3,
     & kd))/(2.*dx(1))
       vxyy22r(i1,i2,i3,kd)=( vyy22r(i1+1,i2,i3,kd)-vyy22r(i1-1,i2,i3,
     & kd))/(2.*dx(0))
       vxxxx22r(i1,i2,i3,kd)=(6.*v(i1,i2,i3,kd)-4.*(v(i1+1,i2,i3,kd)+v(
     & i1-1,i2,i3,kd))+(v(i1+2,i2,i3,kd)+v(i1-2,i2,i3,kd)) )/(dx(0)**
     & 4)
       vyyyy22r(i1,i2,i3,kd)=(6.*v(i1,i2,i3,kd)-4.*(v(i1,i2+1,i3,kd)+v(
     & i1,i2-1,i3,kd))+(v(i1,i2+2,i3,kd)+v(i1,i2-2,i3,kd)) )/(dx(1)**
     & 4)
       vxxyy22r(i1,i2,i3,kd)=( 4.*v(i1,i2,i3,kd)     -2.*(v(i1+1,i2,i3,
     & kd)+v(i1-1,i2,i3,kd)+v(i1,i2+1,i3,kd)+v(i1,i2-1,i3,kd))   +   (
     & v(i1+1,i2+1,i3,kd)+v(i1-1,i2+1,i3,kd)+v(i1+1,i2-1,i3,kd)+v(i1-
     & 1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)
       ! 2D laplacian squared = v.xxxx + 2 v.xxyy + v.yyyy
       vLapSq22r(i1,i2,i3,kd)= ( 6.*v(i1,i2,i3,kd)   - 4.*(v(i1+1,i2,
     & i3,kd)+v(i1-1,i2,i3,kd))    +(v(i1+2,i2,i3,kd)+v(i1-2,i2,i3,kd)
     & ) )/(dx(0)**4) +( 6.*v(i1,i2,i3,kd)    -4.*(v(i1,i2+1,i3,kd)+v(
     & i1,i2-1,i3,kd))    +(v(i1,i2+2,i3,kd)+v(i1,i2-2,i3,kd)) )/(dx(
     & 1)**4)  +( 8.*v(i1,i2,i3,kd)     -4.*(v(i1+1,i2,i3,kd)+v(i1-1,
     & i2,i3,kd)+v(i1,i2+1,i3,kd)+v(i1,i2-1,i3,kd))   +2.*(v(i1+1,i2+
     & 1,i3,kd)+v(i1-1,i2+1,i3,kd)+v(i1+1,i2-1,i3,kd)+v(i1-1,i2-1,i3,
     & kd)) )/(dx(0)**2*dx(1)**2)
       vxxx23r(i1,i2,i3,kd)=(-2.*(v(i1+1,i2,i3,kd)-v(i1-1,i2,i3,kd))+(
     & v(i1+2,i2,i3,kd)-v(i1-2,i2,i3,kd)) )*h22(0)*h12(0)
       vyyy23r(i1,i2,i3,kd)=(-2.*(v(i1,i2+1,i3,kd)-v(i1,i2-1,i3,kd))+(
     & v(i1,i2+2,i3,kd)-v(i1,i2-2,i3,kd)) )*h22(1)*h12(1)
       vzzz23r(i1,i2,i3,kd)=(-2.*(v(i1,i2,i3+1,kd)-v(i1,i2,i3-1,kd))+(
     & v(i1,i2,i3+2,kd)-v(i1,i2,i3-2,kd)) )*h22(1)*h12(2)
       vxxy23r(i1,i2,i3,kd)=( vxx22r(i1,i2+1,i3,kd)-vxx22r(i1,i2-1,i3,
     & kd))/(2.*dx(1))
       vxyy23r(i1,i2,i3,kd)=( vyy22r(i1+1,i2,i3,kd)-vyy22r(i1-1,i2,i3,
     & kd))/(2.*dx(0))
       vxxz23r(i1,i2,i3,kd)=( vxx22r(i1,i2,i3+1,kd)-vxx22r(i1,i2,i3-1,
     & kd))/(2.*dx(2))
       vyyz23r(i1,i2,i3,kd)=( vyy22r(i1,i2,i3+1,kd)-vyy22r(i1,i2,i3-1,
     & kd))/(2.*dx(2))
       vxzz23r(i1,i2,i3,kd)=( vzz22r(i1+1,i2,i3,kd)-vzz22r(i1-1,i2,i3,
     & kd))/(2.*dx(0))
       vyzz23r(i1,i2,i3,kd)=( vzz22r(i1,i2+1,i3,kd)-vzz22r(i1,i2-1,i3,
     & kd))/(2.*dx(1))
       vxxxx23r(i1,i2,i3,kd)=(6.*v(i1,i2,i3,kd)-4.*(v(i1+1,i2,i3,kd)+v(
     & i1-1,i2,i3,kd))+(v(i1+2,i2,i3,kd)+v(i1-2,i2,i3,kd)) )/(dx(0)**
     & 4)
       vyyyy23r(i1,i2,i3,kd)=(6.*v(i1,i2,i3,kd)-4.*(v(i1,i2+1,i3,kd)+v(
     & i1,i2-1,i3,kd))+(v(i1,i2+2,i3,kd)+v(i1,i2-2,i3,kd)) )/(dx(1)**
     & 4)
       vzzzz23r(i1,i2,i3,kd)=(6.*v(i1,i2,i3,kd)-4.*(v(i1,i2,i3+1,kd)+v(
     & i1,i2,i3-1,kd))+(v(i1,i2,i3+2,kd)+v(i1,i2,i3-2,kd)) )/(dx(2)**
     & 4)
       vxxyy23r(i1,i2,i3,kd)=( 4.*v(i1,i2,i3,kd)     -2.*(v(i1+1,i2,i3,
     & kd)+v(i1-1,i2,i3,kd)+v(i1,i2+1,i3,kd)+v(i1,i2-1,i3,kd))   +   (
     & v(i1+1,i2+1,i3,kd)+v(i1-1,i2+1,i3,kd)+v(i1+1,i2-1,i3,kd)+v(i1-
     & 1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)
       vxxzz23r(i1,i2,i3,kd)=( 4.*v(i1,i2,i3,kd)     -2.*(v(i1+1,i2,i3,
     & kd)+v(i1-1,i2,i3,kd)+v(i1,i2,i3+1,kd)+v(i1,i2,i3-1,kd))   +   (
     & v(i1+1,i2,i3+1,kd)+v(i1-1,i2,i3+1,kd)+v(i1+1,i2,i3-1,kd)+v(i1-
     & 1,i2,i3-1,kd)) )/(dx(0)**2*dx(2)**2)
       vyyzz23r(i1,i2,i3,kd)=( 4.*v(i1,i2,i3,kd)     -2.*(v(i1,i2+1,i3,
     & kd)  +v(i1,i2-1,i3,kd)+  v(i1,i2  ,i3+1,kd)+v(i1,i2  ,i3-1,kd))
     &    +   (v(i1,i2+1,i3+1,kd)+v(i1,i2-1,i3+1,kd)+v(i1,i2+1,i3-1,
     & kd)+v(i1,i2-1,i3-1,kd)) )/(dx(1)**2*dx(2)**2)
       ! 3D laplacian squared = v.xxxx + v.yyyy + v.zzzz + 2 (v.xxyy + v.xxzz + v.yyzz )
       vLapSq23r(i1,i2,i3,kd)= ( 6.*v(i1,i2,i3,kd)   - 4.*(v(i1+1,i2,
     & i3,kd)+v(i1-1,i2,i3,kd))    +(v(i1+2,i2,i3,kd)+v(i1-2,i2,i3,kd)
     & ) )/(dx(0)**4) +( 6.*v(i1,i2,i3,kd)    -4.*(v(i1,i2+1,i3,kd)+v(
     & i1,i2-1,i3,kd))    +(v(i1,i2+2,i3,kd)+v(i1,i2-2,i3,kd)) )/(dx(
     & 1)**4)  +( 6.*v(i1,i2,i3,kd)    -4.*(v(i1,i2,i3+1,kd)+v(i1,i2,
     & i3-1,kd))    +(v(i1,i2,i3+2,kd)+v(i1,i2,i3-2,kd)) )/(dx(2)**4) 
     &  +( 8.*v(i1,i2,i3,kd)     -4.*(v(i1+1,i2,i3,kd)  +v(i1-1,i2,i3,
     & kd)  +v(i1  ,i2+1,i3,kd)+v(i1  ,i2-1,i3,kd))   +2.*(v(i1+1,i2+
     & 1,i3,kd)+v(i1-1,i2+1,i3,kd)+v(i1+1,i2-1,i3,kd)+v(i1-1,i2-1,i3,
     & kd)) )/(dx(0)**2*dx(1)**2)+( 8.*v(i1,i2,i3,kd)     -4.*(v(i1+1,
     & i2,i3,kd)  +v(i1-1,i2,i3,kd)  +v(i1  ,i2,i3+1,kd)+v(i1  ,i2,i3-
     & 1,kd))   +2.*(v(i1+1,i2,i3+1,kd)+v(i1-1,i2,i3+1,kd)+v(i1+1,i2,
     & i3-1,kd)+v(i1-1,i2,i3-1,kd)) )/(dx(0)**2*dx(2)**2)+( 8.*v(i1,
     & i2,i3,kd)     -4.*(v(i1,i2+1,i3,kd)  +v(i1,i2-1,i3,kd)  +v(i1,
     & i2  ,i3+1,kd)+v(i1,i2  ,i3-1,kd))   +2.*(v(i1,i2+1,i3+1,kd)+v(
     & i1,i2-1,i3+1,kd)+v(i1,i2+1,i3-1,kd)+v(i1,i2-1,i3-1,kd)) )/(dx(
     & 1)**2*dx(2)**2)
       d14(kd) = 1./(12.*dr(kd))
       d24(kd) = 1./(12.*dr(kd)**2)
       vr4(i1,i2,i3,kd)=(8.*(v(i1+1,i2,i3,kd)-v(i1-1,i2,i3,kd))-(v(i1+
     & 2,i2,i3,kd)-v(i1-2,i2,i3,kd)))*d14(0)
       vs4(i1,i2,i3,kd)=(8.*(v(i1,i2+1,i3,kd)-v(i1,i2-1,i3,kd))-(v(i1,
     & i2+2,i3,kd)-v(i1,i2-2,i3,kd)))*d14(1)
       vt4(i1,i2,i3,kd)=(8.*(v(i1,i2,i3+1,kd)-v(i1,i2,i3-1,kd))-(v(i1,
     & i2,i3+2,kd)-v(i1,i2,i3-2,kd)))*d14(2)
       vrr4(i1,i2,i3,kd)=(-30.*v(i1,i2,i3,kd)+16.*(v(i1+1,i2,i3,kd)+v(
     & i1-1,i2,i3,kd))-(v(i1+2,i2,i3,kd)+v(i1-2,i2,i3,kd)) )*d24(0)
       vss4(i1,i2,i3,kd)=(-30.*v(i1,i2,i3,kd)+16.*(v(i1,i2+1,i3,kd)+v(
     & i1,i2-1,i3,kd))-(v(i1,i2+2,i3,kd)+v(i1,i2-2,i3,kd)) )*d24(1)
       vtt4(i1,i2,i3,kd)=(-30.*v(i1,i2,i3,kd)+16.*(v(i1,i2,i3+1,kd)+v(
     & i1,i2,i3-1,kd))-(v(i1,i2,i3+2,kd)+v(i1,i2,i3-2,kd)) )*d24(2)
       vrs4(i1,i2,i3,kd)=(8.*(vr4(i1,i2+1,i3,kd)-vr4(i1,i2-1,i3,kd))-(
     & vr4(i1,i2+2,i3,kd)-vr4(i1,i2-2,i3,kd)))*d14(1)
       vrt4(i1,i2,i3,kd)=(8.*(vr4(i1,i2,i3+1,kd)-vr4(i1,i2,i3-1,kd))-(
     & vr4(i1,i2,i3+2,kd)-vr4(i1,i2,i3-2,kd)))*d14(2)
       vst4(i1,i2,i3,kd)=(8.*(vs4(i1,i2,i3+1,kd)-vs4(i1,i2,i3-1,kd))-(
     & vs4(i1,i2,i3+2,kd)-vs4(i1,i2,i3-2,kd)))*d14(2)
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
       vx41(i1,i2,i3,kd)= rx(i1,i2,i3)*vr4(i1,i2,i3,kd)
       vy41(i1,i2,i3,kd)=0
       vz41(i1,i2,i3,kd)=0
       vx42(i1,i2,i3,kd)= rx(i1,i2,i3)*vr4(i1,i2,i3,kd)+sx(i1,i2,i3)*
     & vs4(i1,i2,i3,kd)
       vy42(i1,i2,i3,kd)= ry(i1,i2,i3)*vr4(i1,i2,i3,kd)+sy(i1,i2,i3)*
     & vs4(i1,i2,i3,kd)
       vz42(i1,i2,i3,kd)=0
       vx43(i1,i2,i3,kd)=rx(i1,i2,i3)*vr4(i1,i2,i3,kd)+sx(i1,i2,i3)*
     & vs4(i1,i2,i3,kd)+tx(i1,i2,i3)*vt4(i1,i2,i3,kd)
       vy43(i1,i2,i3,kd)=ry(i1,i2,i3)*vr4(i1,i2,i3,kd)+sy(i1,i2,i3)*
     & vs4(i1,i2,i3,kd)+ty(i1,i2,i3)*vt4(i1,i2,i3,kd)
       vz43(i1,i2,i3,kd)=rz(i1,i2,i3)*vr4(i1,i2,i3,kd)+sz(i1,i2,i3)*
     & vs4(i1,i2,i3,kd)+tz(i1,i2,i3)*vt4(i1,i2,i3,kd)
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
       vxx41(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*vrr4(i1,i2,i3,kd)+(rxx42(
     & i1,i2,i3))*vr4(i1,i2,i3,kd)
       vyy41(i1,i2,i3,kd)=0
       vxy41(i1,i2,i3,kd)=0
       vxz41(i1,i2,i3,kd)=0
       vyz41(i1,i2,i3,kd)=0
       vzz41(i1,i2,i3,kd)=0
       vlaplacian41(i1,i2,i3,kd)=vxx41(i1,i2,i3,kd)
       vxx42(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*vrr4(i1,i2,i3,kd)+2.*(rx(
     & i1,i2,i3)*sx(i1,i2,i3))*vrs4(i1,i2,i3,kd)+(sx(i1,i2,i3)**2)*
     & vss4(i1,i2,i3,kd)+(rxx42(i1,i2,i3))*vr4(i1,i2,i3,kd)+(sxx42(i1,
     & i2,i3))*vs4(i1,i2,i3,kd)
       vyy42(i1,i2,i3,kd)=(ry(i1,i2,i3)**2)*vrr4(i1,i2,i3,kd)+2.*(ry(
     & i1,i2,i3)*sy(i1,i2,i3))*vrs4(i1,i2,i3,kd)+(sy(i1,i2,i3)**2)*
     & vss4(i1,i2,i3,kd)+(ryy42(i1,i2,i3))*vr4(i1,i2,i3,kd)+(syy42(i1,
     & i2,i3))*vs4(i1,i2,i3,kd)
       vxy42(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*vrr4(i1,i2,i3,kd)+(
     & rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*vrs4(i1,
     & i2,i3,kd)+sx(i1,i2,i3)*sy(i1,i2,i3)*vss4(i1,i2,i3,kd)+rxy42(i1,
     & i2,i3)*vr4(i1,i2,i3,kd)+sxy42(i1,i2,i3)*vs4(i1,i2,i3,kd)
       vxz42(i1,i2,i3,kd)=0
       vyz42(i1,i2,i3,kd)=0
       vzz42(i1,i2,i3,kd)=0
       vlaplacian42(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2)*
     & vrr4(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*
     & sy(i1,i2,i3))*vrs4(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**
     & 2)*vss4(i1,i2,i3,kd)+(rxx42(i1,i2,i3)+ryy42(i1,i2,i3))*vr4(i1,
     & i2,i3,kd)+(sxx42(i1,i2,i3)+syy42(i1,i2,i3))*vs4(i1,i2,i3,kd)
       vxx43(i1,i2,i3,kd)=rx(i1,i2,i3)**2*vrr4(i1,i2,i3,kd)+sx(i1,i2,
     & i3)**2*vss4(i1,i2,i3,kd)+tx(i1,i2,i3)**2*vtt4(i1,i2,i3,kd)+2.*
     & rx(i1,i2,i3)*sx(i1,i2,i3)*vrs4(i1,i2,i3,kd)+2.*rx(i1,i2,i3)*tx(
     & i1,i2,i3)*vrt4(i1,i2,i3,kd)+2.*sx(i1,i2,i3)*tx(i1,i2,i3)*vst4(
     & i1,i2,i3,kd)+rxx43(i1,i2,i3)*vr4(i1,i2,i3,kd)+sxx43(i1,i2,i3)*
     & vs4(i1,i2,i3,kd)+txx43(i1,i2,i3)*vt4(i1,i2,i3,kd)
       vyy43(i1,i2,i3,kd)=ry(i1,i2,i3)**2*vrr4(i1,i2,i3,kd)+sy(i1,i2,
     & i3)**2*vss4(i1,i2,i3,kd)+ty(i1,i2,i3)**2*vtt4(i1,i2,i3,kd)+2.*
     & ry(i1,i2,i3)*sy(i1,i2,i3)*vrs4(i1,i2,i3,kd)+2.*ry(i1,i2,i3)*ty(
     & i1,i2,i3)*vrt4(i1,i2,i3,kd)+2.*sy(i1,i2,i3)*ty(i1,i2,i3)*vst4(
     & i1,i2,i3,kd)+ryy43(i1,i2,i3)*vr4(i1,i2,i3,kd)+syy43(i1,i2,i3)*
     & vs4(i1,i2,i3,kd)+tyy43(i1,i2,i3)*vt4(i1,i2,i3,kd)
       vzz43(i1,i2,i3,kd)=rz(i1,i2,i3)**2*vrr4(i1,i2,i3,kd)+sz(i1,i2,
     & i3)**2*vss4(i1,i2,i3,kd)+tz(i1,i2,i3)**2*vtt4(i1,i2,i3,kd)+2.*
     & rz(i1,i2,i3)*sz(i1,i2,i3)*vrs4(i1,i2,i3,kd)+2.*rz(i1,i2,i3)*tz(
     & i1,i2,i3)*vrt4(i1,i2,i3,kd)+2.*sz(i1,i2,i3)*tz(i1,i2,i3)*vst4(
     & i1,i2,i3,kd)+rzz43(i1,i2,i3)*vr4(i1,i2,i3,kd)+szz43(i1,i2,i3)*
     & vs4(i1,i2,i3,kd)+tzz43(i1,i2,i3)*vt4(i1,i2,i3,kd)
       vxy43(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*vrr4(i1,i2,i3,kd)+
     & sx(i1,i2,i3)*sy(i1,i2,i3)*vss4(i1,i2,i3,kd)+tx(i1,i2,i3)*ty(i1,
     & i2,i3)*vtt4(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,
     & i3)*sx(i1,i2,i3))*vrs4(i1,i2,i3,kd)+(rx(i1,i2,i3)*ty(i1,i2,i3)+
     & ry(i1,i2,i3)*tx(i1,i2,i3))*vrt4(i1,i2,i3,kd)+(sx(i1,i2,i3)*ty(
     & i1,i2,i3)+sy(i1,i2,i3)*tx(i1,i2,i3))*vst4(i1,i2,i3,kd)+rxy43(
     & i1,i2,i3)*vr4(i1,i2,i3,kd)+sxy43(i1,i2,i3)*vs4(i1,i2,i3,kd)+
     & txy43(i1,i2,i3)*vt4(i1,i2,i3,kd)
       vxz43(i1,i2,i3,kd)=rx(i1,i2,i3)*rz(i1,i2,i3)*vrr4(i1,i2,i3,kd)+
     & sx(i1,i2,i3)*sz(i1,i2,i3)*vss4(i1,i2,i3,kd)+tx(i1,i2,i3)*tz(i1,
     & i2,i3)*vtt4(i1,i2,i3,kd)+(rx(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,i2,
     & i3)*sx(i1,i2,i3))*vrs4(i1,i2,i3,kd)+(rx(i1,i2,i3)*tz(i1,i2,i3)+
     & rz(i1,i2,i3)*tx(i1,i2,i3))*vrt4(i1,i2,i3,kd)+(sx(i1,i2,i3)*tz(
     & i1,i2,i3)+sz(i1,i2,i3)*tx(i1,i2,i3))*vst4(i1,i2,i3,kd)+rxz43(
     & i1,i2,i3)*vr4(i1,i2,i3,kd)+sxz43(i1,i2,i3)*vs4(i1,i2,i3,kd)+
     & txz43(i1,i2,i3)*vt4(i1,i2,i3,kd)
       vyz43(i1,i2,i3,kd)=ry(i1,i2,i3)*rz(i1,i2,i3)*vrr4(i1,i2,i3,kd)+
     & sy(i1,i2,i3)*sz(i1,i2,i3)*vss4(i1,i2,i3,kd)+ty(i1,i2,i3)*tz(i1,
     & i2,i3)*vtt4(i1,i2,i3,kd)+(ry(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,i2,
     & i3)*sy(i1,i2,i3))*vrs4(i1,i2,i3,kd)+(ry(i1,i2,i3)*tz(i1,i2,i3)+
     & rz(i1,i2,i3)*ty(i1,i2,i3))*vrt4(i1,i2,i3,kd)+(sy(i1,i2,i3)*tz(
     & i1,i2,i3)+sz(i1,i2,i3)*ty(i1,i2,i3))*vst4(i1,i2,i3,kd)+ryz43(
     & i1,i2,i3)*vr4(i1,i2,i3,kd)+syz43(i1,i2,i3)*vs4(i1,i2,i3,kd)+
     & tyz43(i1,i2,i3)*vt4(i1,i2,i3,kd)
       vlaplacian43(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2+rz(
     & i1,i2,i3)**2)*vrr4(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**
     & 2+sz(i1,i2,i3)**2)*vss4(i1,i2,i3,kd)+(tx(i1,i2,i3)**2+ty(i1,i2,
     & i3)**2+tz(i1,i2,i3)**2)*vtt4(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(
     & i1,i2,i3)+ ry(i1,i2,i3)*sy(i1,i2,i3)+rz(i1,i2,i3)*sz(i1,i2,i3))
     & *vrs4(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*tx(i1,i2,i3)+ ry(i1,i2,i3)*
     & ty(i1,i2,i3)+rz(i1,i2,i3)*tz(i1,i2,i3))*vrt4(i1,i2,i3,kd)+2.*(
     & sx(i1,i2,i3)*tx(i1,i2,i3)+ sy(i1,i2,i3)*ty(i1,i2,i3)+sz(i1,i2,
     & i3)*tz(i1,i2,i3))*vst4(i1,i2,i3,kd)+(rxx43(i1,i2,i3)+ryy43(i1,
     & i2,i3)+rzz43(i1,i2,i3))*vr4(i1,i2,i3,kd)+(sxx43(i1,i2,i3)+
     & syy43(i1,i2,i3)+szz43(i1,i2,i3))*vs4(i1,i2,i3,kd)+(txx43(i1,i2,
     & i3)+tyy43(i1,i2,i3)+tzz43(i1,i2,i3))*vt4(i1,i2,i3,kd)
c============================================================================================
c Define derivatives for a rectangular grid
c
c============================================================================================
       h41(kd) = 1./(12.*dx(kd))
       h42(kd) = 1./(12.*dx(kd)**2)
       vx43r(i1,i2,i3,kd)=(8.*(v(i1+1,i2,i3,kd)-v(i1-1,i2,i3,kd))-(v(
     & i1+2,i2,i3,kd)-v(i1-2,i2,i3,kd)))*h41(0)
       vy43r(i1,i2,i3,kd)=(8.*(v(i1,i2+1,i3,kd)-v(i1,i2-1,i3,kd))-(v(
     & i1,i2+2,i3,kd)-v(i1,i2-2,i3,kd)))*h41(1)
       vz43r(i1,i2,i3,kd)=(8.*(v(i1,i2,i3+1,kd)-v(i1,i2,i3-1,kd))-(v(
     & i1,i2,i3+2,kd)-v(i1,i2,i3-2,kd)))*h41(2)
       vxx43r(i1,i2,i3,kd)=( -30.*v(i1,i2,i3,kd)+16.*(v(i1+1,i2,i3,kd)+
     & v(i1-1,i2,i3,kd))-(v(i1+2,i2,i3,kd)+v(i1-2,i2,i3,kd)) )*h42(0)
       vyy43r(i1,i2,i3,kd)=( -30.*v(i1,i2,i3,kd)+16.*(v(i1,i2+1,i3,kd)+
     & v(i1,i2-1,i3,kd))-(v(i1,i2+2,i3,kd)+v(i1,i2-2,i3,kd)) )*h42(1)
       vzz43r(i1,i2,i3,kd)=( -30.*v(i1,i2,i3,kd)+16.*(v(i1,i2,i3+1,kd)+
     & v(i1,i2,i3-1,kd))-(v(i1,i2,i3+2,kd)+v(i1,i2,i3-2,kd)) )*h42(2)
       vxy43r(i1,i2,i3,kd)=( (v(i1+2,i2+2,i3,kd)-v(i1-2,i2+2,i3,kd)- v(
     & i1+2,i2-2,i3,kd)+v(i1-2,i2-2,i3,kd)) +8.*(v(i1-1,i2+2,i3,kd)-v(
     & i1-1,i2-2,i3,kd)-v(i1+1,i2+2,i3,kd)+v(i1+1,i2-2,i3,kd) +v(i1+2,
     & i2-1,i3,kd)-v(i1-2,i2-1,i3,kd)-v(i1+2,i2+1,i3,kd)+v(i1-2,i2+1,
     & i3,kd))+64.*(v(i1+1,i2+1,i3,kd)-v(i1-1,i2+1,i3,kd)- v(i1+1,i2-
     & 1,i3,kd)+v(i1-1,i2-1,i3,kd)))*(h41(0)*h41(1))
       vxz43r(i1,i2,i3,kd)=( (v(i1+2,i2,i3+2,kd)-v(i1-2,i2,i3+2,kd)-v(
     & i1+2,i2,i3-2,kd)+v(i1-2,i2,i3-2,kd)) +8.*(v(i1-1,i2,i3+2,kd)-v(
     & i1-1,i2,i3-2,kd)-v(i1+1,i2,i3+2,kd)+v(i1+1,i2,i3-2,kd) +v(i1+2,
     & i2,i3-1,kd)-v(i1-2,i2,i3-1,kd)- v(i1+2,i2,i3+1,kd)+v(i1-2,i2,
     & i3+1,kd)) +64.*(v(i1+1,i2,i3+1,kd)-v(i1-1,i2,i3+1,kd)-v(i1+1,
     & i2,i3-1,kd)+v(i1-1,i2,i3-1,kd)) )*(h41(0)*h41(2))
       vyz43r(i1,i2,i3,kd)=( (v(i1,i2+2,i3+2,kd)-v(i1,i2-2,i3+2,kd)-v(
     & i1,i2+2,i3-2,kd)+v(i1,i2-2,i3-2,kd)) +8.*(v(i1,i2-1,i3+2,kd)-v(
     & i1,i2-1,i3-2,kd)-v(i1,i2+1,i3+2,kd)+v(i1,i2+1,i3-2,kd) +v(i1,
     & i2+2,i3-1,kd)-v(i1,i2-2,i3-1,kd)-v(i1,i2+2,i3+1,kd)+v(i1,i2-2,
     & i3+1,kd)) +64.*(v(i1,i2+1,i3+1,kd)-v(i1,i2-1,i3+1,kd)-v(i1,i2+
     & 1,i3-1,kd)+v(i1,i2-1,i3-1,kd)) )*(h41(1)*h41(2))
       vx41r(i1,i2,i3,kd)= vx43r(i1,i2,i3,kd)
       vy41r(i1,i2,i3,kd)= vy43r(i1,i2,i3,kd)
       vz41r(i1,i2,i3,kd)= vz43r(i1,i2,i3,kd)
       vxx41r(i1,i2,i3,kd)= vxx43r(i1,i2,i3,kd)
       vyy41r(i1,i2,i3,kd)= vyy43r(i1,i2,i3,kd)
       vzz41r(i1,i2,i3,kd)= vzz43r(i1,i2,i3,kd)
       vxy41r(i1,i2,i3,kd)= vxy43r(i1,i2,i3,kd)
       vxz41r(i1,i2,i3,kd)= vxz43r(i1,i2,i3,kd)
       vyz41r(i1,i2,i3,kd)= vyz43r(i1,i2,i3,kd)
       vlaplacian41r(i1,i2,i3,kd)=vxx43r(i1,i2,i3,kd)
       vx42r(i1,i2,i3,kd)= vx43r(i1,i2,i3,kd)
       vy42r(i1,i2,i3,kd)= vy43r(i1,i2,i3,kd)
       vz42r(i1,i2,i3,kd)= vz43r(i1,i2,i3,kd)
       vxx42r(i1,i2,i3,kd)= vxx43r(i1,i2,i3,kd)
       vyy42r(i1,i2,i3,kd)= vyy43r(i1,i2,i3,kd)
       vzz42r(i1,i2,i3,kd)= vzz43r(i1,i2,i3,kd)
       vxy42r(i1,i2,i3,kd)= vxy43r(i1,i2,i3,kd)
       vxz42r(i1,i2,i3,kd)= vxz43r(i1,i2,i3,kd)
       vyz42r(i1,i2,i3,kd)= vyz43r(i1,i2,i3,kd)
       vlaplacian42r(i1,i2,i3,kd)=vxx43r(i1,i2,i3,kd)+vyy43r(i1,i2,i3,
     & kd)
       vlaplacian43r(i1,i2,i3,kd)=vxx43r(i1,i2,i3,kd)+vyy43r(i1,i2,i3,
     & kd)+vzz43r(i1,i2,i3,kd)
      !...........end   statement functions
       dx(0) =rpar(0)
       dx(1) =rpar(1)
       dx(2) =rpar(2)
       dr(0) =rpar(3)
       dr(1) =rpar(4)
       dr(2) =rpar(5)
       t     =rpar(6)
       ! This are provided on input and output (max is taken over all grids)
       divDmax = rpar(7)  ! max(div(D)) returned here
       gradDmax= rpar(8)  ! max(grad(D)) returned here
       divBmax = rpar(9)  ! max(div(B)) returned here
       gradBmax= rpar(10) ! max(grad(B)) returned here
       dt      = rpar(11) ! added March 26, 2021 -- was missing
       saveDivergence     =ipar(0)
       method             =ipar(1)
       gridType           =ipar(2)
       orderOfAccuracy    =ipar(3)
       ex                 =ipar(4)
       ey                 =ipar(5)
       ez                 =ipar(6)
       hx                 =ipar(7)
       hy                 =ipar(8)
       hz                 =ipar(9)
       dispersionModel    =ipar(10)
       grid               =ipar(11)
       solveForAllFields  =ipar(12)
       numberOfMaterialRegions = ipar(13)
       debug              =ipar(14)
       if( t.le.2*dt .and. debug.gt.0 )then
         write(*,*) 'Inside getDiv3dOrder2c...'
         write(*,'("solveForAllFields=",i2)') solveForAllFields
         write(*,'("dispersionModel=",i2)') dispersionModel
       end if
       if( method.eq.bamx )then
         if( numberOfMaterialRegions.gt.maxRegions )then
           write(*,*) 'getDiv: Error: numberOfMaterialRegions=',
     & numberOfMaterialRegions,' is bigger than maxRegions=',
     & maxRegions
           write(*,*) 'FIX ME BILL!'
           stop 1002
         end if
         ! 3x3 Material matrix for TEz polarization
         ! We use ex=0,ey=1 and hz=5 entries in K0i(0:5,0:5) 
         do mr=0,numberOfMaterialRegions-1
           K03(0,0,mr) = K0(0,0,mr)
           K03(0,1,mr) = K0(0,1,mr)
           K03(0,2,mr) = K0(0,5,mr)
           K03(1,0,mr) = K0(1,0,mr)
           K03(1,1,mr) = K0(1,1,mr)
           K03(1,2,mr) = K0(1,5,mr)
           K03(2,0,mr) = K0(5,0,mr)
           K03(2,1,mr) = K0(5,1,mr)
           K03(2,2,mr) = K0(5,5,mr)
         end do
       end if
       if( method.eq.bamx .and. t.lt.dt )then
         write(*,*) 'numberOfMaterialRegions=',numberOfMaterialRegions
         do mr=0,numberOfMaterialRegions-1
           write(*,*) 'Material region=',mr
           write(*,'("K0=",6("[",6(f6.3,1x),"]",/,4x))') ((K0(i1,i2,mr)
     & ,i1=0,5),i2=0,5)
           if( solveForAllFields .eq. 0 )then
             write(*,'("K03 =",3("[",3(f6.3,1x),"]",/,4x))') ((K03(i1,
     & i2,mr),i1=0,2),i2=0,2)
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
       if( method.eq.bamx .and. dispersionModel.ne.noDispersion )then
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
           if( debug.gt.1 )then
             write(*,'("--getDiv-- dispersionModel=",i4," 
     & numPolarizationTerms=",i6)') dispersionModel,
     & numPolarizationTerms
           end if
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
       if( t.eq.0. .and. dispersionModel.ne.noDispersion .and. debug>1 
     & )then
          write(*,'("--getDiv-- dispersionModel=",i4)') dispersionModel
       end if
        ! We need to evaluate D and B at extra points so that we can take the divergence
        numGhost = orderOfAccuracy/2
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
        ! divDmax=0.
        ! divBmax=0.
        ! gradDmax=0.
        ! gradBmax=0.
        if( method.ne.bamx )then
         ! --- Isotropic Maxwell ----
            if( t.lt.2*dt )then
              write(*,'("computeDiv dim=3 order=2 grid=curvilinear ... 
     & t=",e10.2)') t
            end if
            ! FINISH ME -- divMax, gradmax ...
            stop 1234
            ! Loop bounds -- exclude absorbing layers ? 
              do i3=n3a,n3b
              do i2=n2a,n2b
              do i1=n1a,n1b
                if( mask(i1,i2,i3).gt.0 )then
               ! --- 3D -----
                  divD(i1,i2,i3) = vx23(i1,i2,i3,ex) + vy23(i1,i2,i3,
     & ey) + vz23(i1,i2,i3,ez)
             divDmax=max(divDmax,divD(i1,i2,i3))
                end if
              end do
              end do
              end do
        else
          ! --- BA EQUATIONS ---
          if( dispersionModel.eq.noDispersion )then
            if( nd.eq.2 .and. solveForAllFields.eq.0 )then
                if( t.lt.2*dt )then
                  write(*,'("computeInducedFields polar=TEZ... t=",
     & e10.2)') t
                end if
                ! loop bounds -- include ghost points 
                mr=0  ! default value for a single material
                  do i3=m3a,m3b
                  do i2=m2a,m2b
                  do i1=m1a,m1b
                    if( mask(i1,i2,i3).gt.0 )then
                  if( numberOfMaterialRegions.gt.1 )then
                    mr = matMask(i1,i2,i3)
                    if( mr.lt.0 .or. mr.ge.numberOfMaterialRegions )
     & then  ! do this for now
                       stop 9999
                    end if
                  end if
                  ! -- Compute [D,B] = K0*U + P 
                    do m=0,2
                      v(i1,i2,i3,m)= K03(m,0,mr)*u(i1,i2,i3,0) + K03(m,
     & 1,mr)*u(i1,i2,i3,1) + K03(m,2,mr)*u(i1,i2,i3,2)
                    end do
                    end if
                  end do
                  end do
                  end do
                if( t.lt.2*dt )then
                  write(*,'("computeDiv (BA) dim=3 order=2 
     & grid=curvilinear polar=TEZ... t=",e10.2)') t
                end if
                ! Loop bounds -- exclude super-grid layer ? 
                  do i3=n3a,n3b
                  do i2=n2a,n2b
                  do i1=n1a,n1b
                    if( mask(i1,i2,i3).gt.0 )then
                   ! --- 3D -----
                     ! -- curvilinear grid --
                       d1x = vx23(i1,i2,i3,ex)
                       d2x = vx23(i1,i2,i3,ey)
                       d3x = vx23(i1,i2,i3,ez)
                       d1y = vy23(i1,i2,i3,ex)
                       d2y = vy23(i1,i2,i3,ey)
                       d3y = vy23(i1,i2,i3,ez)
                       d1z = vz23(i1,i2,i3,ex)
                       d2z = vz23(i1,i2,i3,ey)
                       d3z = vz23(i1,i2,i3,ez)
                       b1x = vx23(i1,i2,i3,hx)
                       b2x = vx23(i1,i2,i3,hy)
                       b3x = vx23(i1,i2,i3,hz)
                       b1y = vy23(i1,i2,i3,hx)
                       b2y = vy23(i1,i2,i3,hy)
                       b3y = vy23(i1,i2,i3,hz)
                       b1z = vz23(i1,i2,i3,hx)
                       b2z = vz23(i1,i2,i3,hy)
                       b3z = vz23(i1,i2,i3,hz)
                   dDiv = d1x + d2y + d3z
                   if( saveDivergence.eq.1 )then
                     divD(i1,i2,i3) = dDiv
                   end if
                   divDmax=max(divDmax,abs(dDiv))
                   bDiv = b1x + b2y + b3z
                   if( saveDivergence.eq.1 )then
                     divB(i1,i2,i3)=bDiv
                   end if
                   divBmax = max(divBmax,abs(bDiv))
                   gradDmax = max( gradDmax,abs(d1x),abs(d2x),abs(d3x),
     & abs(d1y),abs(d2y),abs(d3y),abs(d1z),abs(d2z),abs(d3z))
                   gradBmax = max( gradBmax,abs(b1x),abs(b2x),abs(b3x),
     & abs(b1y),abs(b2y),abs(b3y),abs(b1z),abs(b2z),abs(b3z))
                    end if
                  end do
                  end do
                  end do
            else
                if( t.lt.2*dt )then
                  write(*,'("computeInducedFields polar=NONE... t=",
     & e10.2)') t
                end if
                ! loop bounds -- include ghost points 
                mr=0  ! default value for a single material
                  do i3=m3a,m3b
                  do i2=m2a,m2b
                  do i1=m1a,m1b
                    if( mask(i1,i2,i3).gt.0 )then
                  if( numberOfMaterialRegions.gt.1 )then
                    mr = matMask(i1,i2,i3)
                    if( mr.lt.0 .or. mr.ge.numberOfMaterialRegions )
     & then  ! do this for now
                       stop 9999
                    end if
                  end if
                  ! -- Compute [D,B] = K0*U + P 
                    do m=0,5
                      v(i1,i2,i3,m)= K0(m,0,mr)*u(i1,i2,i3,0) + K0(m,1,
     & mr)*u(i1,i2,i3,1) + K0(m,2,mr)*u(i1,i2,i3,2) + K0(m,3,mr)*u(i1,
     & i2,i3,3) + K0(m,4,mr)*u(i1,i2,i3,4) + K0(m,5,mr)*u(i1,i2,i3,5)
                    end do
                    end if
                  end do
                  end do
                  end do
                if( t.lt.2*dt )then
                  write(*,'("computeDiv (BA) dim=3 order=2 
     & grid=curvilinear polar=NONE... t=",e10.2)') t
                end if
                ! Loop bounds -- exclude super-grid layer ? 
                  do i3=n3a,n3b
                  do i2=n2a,n2b
                  do i1=n1a,n1b
                    if( mask(i1,i2,i3).gt.0 )then
                   ! --- 3D -----
                     ! -- curvilinear grid --
                       d1x = vx23(i1,i2,i3,ex)
                       d2x = vx23(i1,i2,i3,ey)
                       d3x = vx23(i1,i2,i3,ez)
                       d1y = vy23(i1,i2,i3,ex)
                       d2y = vy23(i1,i2,i3,ey)
                       d3y = vy23(i1,i2,i3,ez)
                       d1z = vz23(i1,i2,i3,ex)
                       d2z = vz23(i1,i2,i3,ey)
                       d3z = vz23(i1,i2,i3,ez)
                       b1x = vx23(i1,i2,i3,hx)
                       b2x = vx23(i1,i2,i3,hy)
                       b3x = vx23(i1,i2,i3,hz)
                       b1y = vy23(i1,i2,i3,hx)
                       b2y = vy23(i1,i2,i3,hy)
                       b3y = vy23(i1,i2,i3,hz)
                       b1z = vz23(i1,i2,i3,hx)
                       b2z = vz23(i1,i2,i3,hy)
                       b3z = vz23(i1,i2,i3,hz)
                   dDiv = d1x + d2y + d3z
                   if( saveDivergence.eq.1 )then
                     divD(i1,i2,i3) = dDiv
                   end if
                   divDmax=max(divDmax,abs(dDiv))
                   bDiv = b1x + b2y + b3z
                   if( saveDivergence.eq.1 )then
                     divB(i1,i2,i3)=bDiv
                   end if
                   divBmax = max(divBmax,abs(bDiv))
                   gradDmax = max( gradDmax,abs(d1x),abs(d2x),abs(d3x),
     & abs(d1y),abs(d2y),abs(d3y),abs(d1z),abs(d2z),abs(d3z))
                   gradBmax = max( gradBmax,abs(b1x),abs(b2x),abs(b3x),
     & abs(b1y),abs(b2y),abs(b3y),abs(b1z),abs(b2z),abs(b3z))
                    end if
                  end do
                  end do
                  end do
            end if
          else
            if( nd.eq.2 .and. solveForAllFields.eq.0 )then
                if( t.lt.2*dt )then
                  write(*,'("computeInducedFieldsGDM polar=TEZ... t=",
     & e10.2)') t
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
                        pcIndex1(m,mr)=pc ! ** check me
                        pc=pc+2
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
                end do ! end do mr
                mr=0
                ! loop bounds -- include ghost points **finish me**
                ! include ghost 
                mr=0  ! default value for a single material
                  do i3=m3a,m3b
                  do i2=m2a,m2b
                  do i1=m1a,m1b
                    if( mask(i1,i2,i3).gt.0 )then
                  if( numberOfMaterialRegions.gt.1 )then
                    mr = matMask(i1,i2,i3)
                    if( mr.lt.0 .or. mr.ge.numberOfMaterialRegions )
     & then  ! do this for now
                       stop 9999
                    end if
                  end if
                  ! ---- Compute p = SUM_k2 SUM_n  p(i1,i2,i3, n,k1,k2, pc )
                    do m=0,2
                      ptSum(m)=0
                    end do
                  do m=1,numTerms1(mr)
                     ec = ecIndex1(m,mr)
                     pc = pcIndex1(m,mr)
                     ptSum(ec) = ptSum(ec) + p(i1,i2,i3,pc)
                  end do
                  ! -- Compute [D,B] = K0*U + P 
                    do m=0,2
                      v(i1,i2,i3,m)= K03(m,0,mr)*u(i1,i2,i3,0) + K03(m,
     & 1,mr)*u(i1,i2,i3,1) + K03(m,2,mr)*u(i1,i2,i3,2) + ptSum(m)
                    end do
                    end if
                  end do
                  end do
                  end do
                if( t.lt.2*dt )then
                  write(*,'("computeDiv (BA) dim=3 order=2 
     & grid=curvilinear polar=TEZ... t=",e10.2)') t
                end if
                ! Loop bounds -- exclude super-grid layer ? 
                  do i3=n3a,n3b
                  do i2=n2a,n2b
                  do i1=n1a,n1b
                    if( mask(i1,i2,i3).gt.0 )then
                   ! --- 3D -----
                     ! -- curvilinear grid --
                       d1x = vx23(i1,i2,i3,ex)
                       d2x = vx23(i1,i2,i3,ey)
                       d3x = vx23(i1,i2,i3,ez)
                       d1y = vy23(i1,i2,i3,ex)
                       d2y = vy23(i1,i2,i3,ey)
                       d3y = vy23(i1,i2,i3,ez)
                       d1z = vz23(i1,i2,i3,ex)
                       d2z = vz23(i1,i2,i3,ey)
                       d3z = vz23(i1,i2,i3,ez)
                       b1x = vx23(i1,i2,i3,hx)
                       b2x = vx23(i1,i2,i3,hy)
                       b3x = vx23(i1,i2,i3,hz)
                       b1y = vy23(i1,i2,i3,hx)
                       b2y = vy23(i1,i2,i3,hy)
                       b3y = vy23(i1,i2,i3,hz)
                       b1z = vz23(i1,i2,i3,hx)
                       b2z = vz23(i1,i2,i3,hy)
                       b3z = vz23(i1,i2,i3,hz)
                   dDiv = d1x + d2y + d3z
                   if( saveDivergence.eq.1 )then
                     divD(i1,i2,i3) = dDiv
                   end if
                   divDmax=max(divDmax,abs(dDiv))
                   bDiv = b1x + b2y + b3z
                   if( saveDivergence.eq.1 )then
                     divB(i1,i2,i3)=bDiv
                   end if
                   divBmax = max(divBmax,abs(bDiv))
                   gradDmax = max( gradDmax,abs(d1x),abs(d2x),abs(d3x),
     & abs(d1y),abs(d2y),abs(d3y),abs(d1z),abs(d2z),abs(d3z))
                   gradBmax = max( gradBmax,abs(b1x),abs(b2x),abs(b3x),
     & abs(b1y),abs(b2y),abs(b3y),abs(b1z),abs(b2z),abs(b3z))
                    end if
                  end do
                  end do
                  end do
            else
                if( t.lt.2*dt )then
                  write(*,'("computeInducedFieldsGDM polar=NONE... t=",
     & e10.2)') t
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
                        pcIndex1(m,mr)=pc ! ** check me
                        pc=pc+2
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
                end do ! end do mr
                mr=0
                ! loop bounds -- include ghost points **finish me**
                ! include ghost 
                mr=0  ! default value for a single material
                  do i3=m3a,m3b
                  do i2=m2a,m2b
                  do i1=m1a,m1b
                    if( mask(i1,i2,i3).gt.0 )then
                  if( numberOfMaterialRegions.gt.1 )then
                    mr = matMask(i1,i2,i3)
                    if( mr.lt.0 .or. mr.ge.numberOfMaterialRegions )
     & then  ! do this for now
                       stop 9999
                    end if
                  end if
                  ! ---- Compute p = SUM_k2 SUM_n  p(i1,i2,i3, n,k1,k2, pc )
                    do m=0,5
                      ptSum(m)=0
                    end do
                  do m=1,numTerms1(mr)
                     ec = ecIndex1(m,mr)
                     pc = pcIndex1(m,mr)
                     ptSum(ec) = ptSum(ec) + p(i1,i2,i3,pc)
                  end do
                  ! -- Compute [D,B] = K0*U + P 
                    do m=0,5
                      v(i1,i2,i3,m)= K0(m,0,mr)*u(i1,i2,i3,0) + K0(m,1,
     & mr)*u(i1,i2,i3,1) + K0(m,2,mr)*u(i1,i2,i3,2) + K0(m,3,mr)*u(i1,
     & i2,i3,3) + K0(m,4,mr)*u(i1,i2,i3,4) + K0(m,5,mr)*u(i1,i2,i3,5) 
     & + ptSum(m)
                    end do
                    end if
                  end do
                  end do
                  end do
                if( t.lt.2*dt )then
                  write(*,'("computeDiv (BA) dim=3 order=2 
     & grid=curvilinear polar=NONE... t=",e10.2)') t
                end if
                ! Loop bounds -- exclude super-grid layer ? 
                  do i3=n3a,n3b
                  do i2=n2a,n2b
                  do i1=n1a,n1b
                    if( mask(i1,i2,i3).gt.0 )then
                   ! --- 3D -----
                     ! -- curvilinear grid --
                       d1x = vx23(i1,i2,i3,ex)
                       d2x = vx23(i1,i2,i3,ey)
                       d3x = vx23(i1,i2,i3,ez)
                       d1y = vy23(i1,i2,i3,ex)
                       d2y = vy23(i1,i2,i3,ey)
                       d3y = vy23(i1,i2,i3,ez)
                       d1z = vz23(i1,i2,i3,ex)
                       d2z = vz23(i1,i2,i3,ey)
                       d3z = vz23(i1,i2,i3,ez)
                       b1x = vx23(i1,i2,i3,hx)
                       b2x = vx23(i1,i2,i3,hy)
                       b3x = vx23(i1,i2,i3,hz)
                       b1y = vy23(i1,i2,i3,hx)
                       b2y = vy23(i1,i2,i3,hy)
                       b3y = vy23(i1,i2,i3,hz)
                       b1z = vz23(i1,i2,i3,hx)
                       b2z = vz23(i1,i2,i3,hy)
                       b3z = vz23(i1,i2,i3,hz)
                   dDiv = d1x + d2y + d3z
                   if( saveDivergence.eq.1 )then
                     divD(i1,i2,i3) = dDiv
                   end if
                   divDmax=max(divDmax,abs(dDiv))
                   bDiv = b1x + b2y + b3z
                   if( saveDivergence.eq.1 )then
                     divB(i1,i2,i3)=bDiv
                   end if
                   divBmax = max(divBmax,abs(bDiv))
                   gradDmax = max( gradDmax,abs(d1x),abs(d2x),abs(d3x),
     & abs(d1y),abs(d2y),abs(d3y),abs(d1z),abs(d2z),abs(d3z))
                   gradBmax = max( gradBmax,abs(b1x),abs(b2x),abs(b3x),
     & abs(b1y),abs(b2y),abs(b3y),abs(b1z),abs(b2z),abs(b3z))
                    end if
                  end do
                  end do
                  end do
            end if
          end if
        end if
        ! return max values in rpar
        rpar(7) =divDmax
        rpar(8) =gradDMax
        rpar(9) =divBmax
        rpar(10)=gradBMax
        return
        end
