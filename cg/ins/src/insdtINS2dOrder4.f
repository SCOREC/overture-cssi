! This file automatically generated from insdtINS.bf with bpp.
         subroutine insdtINS2dOrder4(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,
     & nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,radiusInverse, 
     &  u,uu, ut,uti,gv,dw,  bc, ipar, rpar, ierr )
        !======================================================================
        !   Compute du/dt for the incompressible NS on rectangular grids
        !     OPTIMIZED version for rectangular grids.
        ! nd : number of space dimensions
        !
        ! gv : gridVelocity for moving grids
        ! uu : for moving grids uu is a workspace to hold u-gv, otherwise uu==u
        ! dw : distance to the wall for some turbulence models
        !======================================================================
         implicit none
         integer nd, n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,
     & nd3b,nd4a,nd4b
         real u(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
         real uu(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
         real ut(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
         real uti(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
         real gv(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1)
         real dw(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
         real xy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1)
         real rsxy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1,0:nd-1)
         real radiusInverse(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
         integer mask(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
         integer bc(0:1,0:2),ierr
         integer ipar(0:*)
         real rpar(0:*)
         !     ---- local variables -----
         integer c,i1,i2,i3,kd,kd3,orderOfAccuracy,gridIsMoving,
     & useWhereMask
         integer gridIsImplicit,implicitOption,implicitMethod,
     & isAxisymmetric,use2ndOrderAD,use4thOrderAD
         integer rc,pc,uc,vc,wc,sc,nc,kc,ec,tc,grid,m,
     & advectPassiveScalar,vsc
         real nu,dt,nuPassiveScalar,adcPassiveScalar
         real dxi,dyi,dzi,dri,dsi,dti,dr2i,ds2i,dt2i
         real ad21,ad22,ad41,ad42,cd22,cd42,adc
         real ad21n,ad22n,ad41n,ad42n,cd22n,cd42n
         real yy,ri
         integer materialFormat
         real t
         integer gridType
         integer rectangular,curvilinear
         parameter( rectangular=0, curvilinear=1 )
         integer turbulenceModel,noTurbulenceModel
         integer baldwinLomax,spalartAllmaras,kEpsilon,kOmega,
     & largeEddySimulation
         parameter (noTurbulenceModel=0,baldwinLomax=1,kEpsilon=2,
     & kOmega=3,spalartAllmaras=4,largeEddySimulation=5 )
         integer pdeModel,standardModel,BoussinesqModel,
     & viscoPlasticModel,twoPhaseFlowModel
         parameter( standardModel=0,BoussinesqModel=1,
     & viscoPlasticModel=2,twoPhaseFlowModel=3 )
         integer upwindOrder,debug,orderOneFix, augmentPlot
         integer advectionOption, centeredAdvection,upwindAdvection,
     & bwenoAdvection
         parameter( centeredAdvection=0, upwindAdvection=1, 
     & bwenoAdvection=2 )
         real au,agu(0:25,0:25) ! for holdings upwind approximations to (a.grad)u
         real drl,wb(0:10) ! for holding the bweno weights
         integer s1,s2,s3,var,vard
         real aur,aus,aguP1,aguP2,aum,aup,gvU,gvV,gvW
         real wpl,wpr,wml,wmr,Fpl,Fpr,Fml,Fmr,Fp,Fm ! for holdings bweno flux and weights
         real Apl,Apr,Bpl,Bpr,Aml,Amr,Bml,Bmr
         real betapl,betapr,betaml,betamr
         real maxLocal
         real ep,ap1,ap2,am1,am2,wp1,wp2,wm1,wm2
         real bp1,bp2,bm1,bm2,wplw,wprw,wmlw,wmrw,up,um
         real aumax,expAd
         !***************! for holdings bweno flux and weights
         integer computeAllTerms,doNotComputeImplicitTerms,
     & computeImplicitTermsSeparately,computeAllWithWeightedImplicit
         parameter( computeAllTerms=0,doNotComputeImplicitTerms=1,
     & computeImplicitTermsSeparately=2,
     & computeAllWithWeightedImplicit=3 )
         real rx,ry,rz,sx,sy,sz,tx,ty,tz
         real dr(0:2), dx(0:2)
         ! for SPAL TM
         real n0,n0x,n0y,n0z
         real cb1, cb2, cv1, sigma, sigmai, kappa, cw1, cw2, cw3, 
     & cw3e6, cv1e3, cd0, cr0
         real chi,chi3,fnu1,fnu2,s,r,g,fw,dKappaSq,nSqBydSq,dd
         real nuT,nuTx,nuTy,nuTz,nuTd
         real u0,u0x,u0y,u0z
         real v0,v0x,v0y,v0z
         real w0,w0x,w0y,w0z
         ! for k-epsilon
         real k0,k0x,k0y,k0z, e0,e0x,e0y,e0z
         real nuP,prod
         real cMu,cEps1,cEps2,sigmaEpsI,sigmaKI
         ! for visco-plastic
         ! real nuVP,etaVP,yieldStressVP,exponentVP,epsVP
         ! real eDotNorm,exp0
         ! real u0xx,u0xy,u0xz,u0yy,u0yz,u0zz
         ! real v0xx,v0xy,v0xz,v0yy,v0yz,v0zz
         ! real w0xx,w0xy,w0xz,w0yy,w0yz,w0zz
         real delta22,delta23,delta42,delta43
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
         !  --- begin statement functions
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
         d14(kd) = 1./(12.*dr(kd))
         d24(kd) = 1./(12.*dr(kd)**2)
         ur4(i1,i2,i3,kd)=(8.*(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))-(u(
     & i1+2,i2,i3,kd)-u(i1-2,i2,i3,kd)))*d14(0)
         us4(i1,i2,i3,kd)=(8.*(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))-(u(
     & i1,i2+2,i3,kd)-u(i1,i2-2,i3,kd)))*d14(1)
         ut4(i1,i2,i3,kd)=(8.*(u(i1,i2,i3+1,kd)-u(i1,i2,i3-1,kd))-(u(
     & i1,i2,i3+2,kd)-u(i1,i2,i3-2,kd)))*d14(2)
         urr4(i1,i2,i3,kd)=(-30.*u(i1,i2,i3,kd)+16.*(u(i1+1,i2,i3,kd)+
     & u(i1-1,i2,i3,kd))-(u(i1+2,i2,i3,kd)+u(i1-2,i2,i3,kd)) )*d24(0)
         uss4(i1,i2,i3,kd)=(-30.*u(i1,i2,i3,kd)+16.*(u(i1,i2+1,i3,kd)+
     & u(i1,i2-1,i3,kd))-(u(i1,i2+2,i3,kd)+u(i1,i2-2,i3,kd)) )*d24(1)
         utt4(i1,i2,i3,kd)=(-30.*u(i1,i2,i3,kd)+16.*(u(i1,i2,i3+1,kd)+
     & u(i1,i2,i3-1,kd))-(u(i1,i2,i3+2,kd)+u(i1,i2,i3-2,kd)) )*d24(2)
         urs4(i1,i2,i3,kd)=(8.*(ur4(i1,i2+1,i3,kd)-ur4(i1,i2-1,i3,kd))-
     & (ur4(i1,i2+2,i3,kd)-ur4(i1,i2-2,i3,kd)))*d14(1)
         urt4(i1,i2,i3,kd)=(8.*(ur4(i1,i2,i3+1,kd)-ur4(i1,i2,i3-1,kd))-
     & (ur4(i1,i2,i3+2,kd)-ur4(i1,i2,i3-2,kd)))*d14(2)
         ust4(i1,i2,i3,kd)=(8.*(us4(i1,i2,i3+1,kd)-us4(i1,i2,i3-1,kd))-
     & (us4(i1,i2,i3+2,kd)-us4(i1,i2,i3-2,kd)))*d14(2)
         rxr4(i1,i2,i3)=(8.*(rx(i1+1,i2,i3)-rx(i1-1,i2,i3))-(rx(i1+2,
     & i2,i3)-rx(i1-2,i2,i3)))*d14(0)
         rxs4(i1,i2,i3)=(8.*(rx(i1,i2+1,i3)-rx(i1,i2-1,i3))-(rx(i1,i2+
     & 2,i3)-rx(i1,i2-2,i3)))*d14(1)
         rxt4(i1,i2,i3)=(8.*(rx(i1,i2,i3+1)-rx(i1,i2,i3-1))-(rx(i1,i2,
     & i3+2)-rx(i1,i2,i3-2)))*d14(2)
         ryr4(i1,i2,i3)=(8.*(ry(i1+1,i2,i3)-ry(i1-1,i2,i3))-(ry(i1+2,
     & i2,i3)-ry(i1-2,i2,i3)))*d14(0)
         rys4(i1,i2,i3)=(8.*(ry(i1,i2+1,i3)-ry(i1,i2-1,i3))-(ry(i1,i2+
     & 2,i3)-ry(i1,i2-2,i3)))*d14(1)
         ryt4(i1,i2,i3)=(8.*(ry(i1,i2,i3+1)-ry(i1,i2,i3-1))-(ry(i1,i2,
     & i3+2)-ry(i1,i2,i3-2)))*d14(2)
         rzr4(i1,i2,i3)=(8.*(rz(i1+1,i2,i3)-rz(i1-1,i2,i3))-(rz(i1+2,
     & i2,i3)-rz(i1-2,i2,i3)))*d14(0)
         rzs4(i1,i2,i3)=(8.*(rz(i1,i2+1,i3)-rz(i1,i2-1,i3))-(rz(i1,i2+
     & 2,i3)-rz(i1,i2-2,i3)))*d14(1)
         rzt4(i1,i2,i3)=(8.*(rz(i1,i2,i3+1)-rz(i1,i2,i3-1))-(rz(i1,i2,
     & i3+2)-rz(i1,i2,i3-2)))*d14(2)
         sxr4(i1,i2,i3)=(8.*(sx(i1+1,i2,i3)-sx(i1-1,i2,i3))-(sx(i1+2,
     & i2,i3)-sx(i1-2,i2,i3)))*d14(0)
         sxs4(i1,i2,i3)=(8.*(sx(i1,i2+1,i3)-sx(i1,i2-1,i3))-(sx(i1,i2+
     & 2,i3)-sx(i1,i2-2,i3)))*d14(1)
         sxt4(i1,i2,i3)=(8.*(sx(i1,i2,i3+1)-sx(i1,i2,i3-1))-(sx(i1,i2,
     & i3+2)-sx(i1,i2,i3-2)))*d14(2)
         syr4(i1,i2,i3)=(8.*(sy(i1+1,i2,i3)-sy(i1-1,i2,i3))-(sy(i1+2,
     & i2,i3)-sy(i1-2,i2,i3)))*d14(0)
         sys4(i1,i2,i3)=(8.*(sy(i1,i2+1,i3)-sy(i1,i2-1,i3))-(sy(i1,i2+
     & 2,i3)-sy(i1,i2-2,i3)))*d14(1)
         syt4(i1,i2,i3)=(8.*(sy(i1,i2,i3+1)-sy(i1,i2,i3-1))-(sy(i1,i2,
     & i3+2)-sy(i1,i2,i3-2)))*d14(2)
         szr4(i1,i2,i3)=(8.*(sz(i1+1,i2,i3)-sz(i1-1,i2,i3))-(sz(i1+2,
     & i2,i3)-sz(i1-2,i2,i3)))*d14(0)
         szs4(i1,i2,i3)=(8.*(sz(i1,i2+1,i3)-sz(i1,i2-1,i3))-(sz(i1,i2+
     & 2,i3)-sz(i1,i2-2,i3)))*d14(1)
         szt4(i1,i2,i3)=(8.*(sz(i1,i2,i3+1)-sz(i1,i2,i3-1))-(sz(i1,i2,
     & i3+2)-sz(i1,i2,i3-2)))*d14(2)
         txr4(i1,i2,i3)=(8.*(tx(i1+1,i2,i3)-tx(i1-1,i2,i3))-(tx(i1+2,
     & i2,i3)-tx(i1-2,i2,i3)))*d14(0)
         txs4(i1,i2,i3)=(8.*(tx(i1,i2+1,i3)-tx(i1,i2-1,i3))-(tx(i1,i2+
     & 2,i3)-tx(i1,i2-2,i3)))*d14(1)
         txt4(i1,i2,i3)=(8.*(tx(i1,i2,i3+1)-tx(i1,i2,i3-1))-(tx(i1,i2,
     & i3+2)-tx(i1,i2,i3-2)))*d14(2)
         tyr4(i1,i2,i3)=(8.*(ty(i1+1,i2,i3)-ty(i1-1,i2,i3))-(ty(i1+2,
     & i2,i3)-ty(i1-2,i2,i3)))*d14(0)
         tys4(i1,i2,i3)=(8.*(ty(i1,i2+1,i3)-ty(i1,i2-1,i3))-(ty(i1,i2+
     & 2,i3)-ty(i1,i2-2,i3)))*d14(1)
         tyt4(i1,i2,i3)=(8.*(ty(i1,i2,i3+1)-ty(i1,i2,i3-1))-(ty(i1,i2,
     & i3+2)-ty(i1,i2,i3-2)))*d14(2)
         tzr4(i1,i2,i3)=(8.*(tz(i1+1,i2,i3)-tz(i1-1,i2,i3))-(tz(i1+2,
     & i2,i3)-tz(i1-2,i2,i3)))*d14(0)
         tzs4(i1,i2,i3)=(8.*(tz(i1,i2+1,i3)-tz(i1,i2-1,i3))-(tz(i1,i2+
     & 2,i3)-tz(i1,i2-2,i3)))*d14(1)
         tzt4(i1,i2,i3)=(8.*(tz(i1,i2,i3+1)-tz(i1,i2,i3-1))-(tz(i1,i2,
     & i3+2)-tz(i1,i2,i3-2)))*d14(2)
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
         rxx42(i1,i2,i3)= rx(i1,i2,i3)*rxr4(i1,i2,i3)+sx(i1,i2,i3)*
     & rxs4(i1,i2,i3)
         rxy42(i1,i2,i3)= ry(i1,i2,i3)*rxr4(i1,i2,i3)+sy(i1,i2,i3)*
     & rxs4(i1,i2,i3)
         rxx43(i1,i2,i3)=rx(i1,i2,i3)*rxr4(i1,i2,i3)+sx(i1,i2,i3)*rxs4(
     & i1,i2,i3)+tx(i1,i2,i3)*rxt4(i1,i2,i3)
         rxy43(i1,i2,i3)=ry(i1,i2,i3)*rxr4(i1,i2,i3)+sy(i1,i2,i3)*rxs4(
     & i1,i2,i3)+ty(i1,i2,i3)*rxt4(i1,i2,i3)
         rxz43(i1,i2,i3)=rz(i1,i2,i3)*rxr4(i1,i2,i3)+sz(i1,i2,i3)*rxs4(
     & i1,i2,i3)+tz(i1,i2,i3)*rxt4(i1,i2,i3)
         ryx42(i1,i2,i3)= rx(i1,i2,i3)*ryr4(i1,i2,i3)+sx(i1,i2,i3)*
     & rys4(i1,i2,i3)
         ryy42(i1,i2,i3)= ry(i1,i2,i3)*ryr4(i1,i2,i3)+sy(i1,i2,i3)*
     & rys4(i1,i2,i3)
         ryx43(i1,i2,i3)=rx(i1,i2,i3)*ryr4(i1,i2,i3)+sx(i1,i2,i3)*rys4(
     & i1,i2,i3)+tx(i1,i2,i3)*ryt4(i1,i2,i3)
         ryy43(i1,i2,i3)=ry(i1,i2,i3)*ryr4(i1,i2,i3)+sy(i1,i2,i3)*rys4(
     & i1,i2,i3)+ty(i1,i2,i3)*ryt4(i1,i2,i3)
         ryz43(i1,i2,i3)=rz(i1,i2,i3)*ryr4(i1,i2,i3)+sz(i1,i2,i3)*rys4(
     & i1,i2,i3)+tz(i1,i2,i3)*ryt4(i1,i2,i3)
         rzx42(i1,i2,i3)= rx(i1,i2,i3)*rzr4(i1,i2,i3)+sx(i1,i2,i3)*
     & rzs4(i1,i2,i3)
         rzy42(i1,i2,i3)= ry(i1,i2,i3)*rzr4(i1,i2,i3)+sy(i1,i2,i3)*
     & rzs4(i1,i2,i3)
         rzx43(i1,i2,i3)=rx(i1,i2,i3)*rzr4(i1,i2,i3)+sx(i1,i2,i3)*rzs4(
     & i1,i2,i3)+tx(i1,i2,i3)*rzt4(i1,i2,i3)
         rzy43(i1,i2,i3)=ry(i1,i2,i3)*rzr4(i1,i2,i3)+sy(i1,i2,i3)*rzs4(
     & i1,i2,i3)+ty(i1,i2,i3)*rzt4(i1,i2,i3)
         rzz43(i1,i2,i3)=rz(i1,i2,i3)*rzr4(i1,i2,i3)+sz(i1,i2,i3)*rzs4(
     & i1,i2,i3)+tz(i1,i2,i3)*rzt4(i1,i2,i3)
         sxx42(i1,i2,i3)= rx(i1,i2,i3)*sxr4(i1,i2,i3)+sx(i1,i2,i3)*
     & sxs4(i1,i2,i3)
         sxy42(i1,i2,i3)= ry(i1,i2,i3)*sxr4(i1,i2,i3)+sy(i1,i2,i3)*
     & sxs4(i1,i2,i3)
         sxx43(i1,i2,i3)=rx(i1,i2,i3)*sxr4(i1,i2,i3)+sx(i1,i2,i3)*sxs4(
     & i1,i2,i3)+tx(i1,i2,i3)*sxt4(i1,i2,i3)
         sxy43(i1,i2,i3)=ry(i1,i2,i3)*sxr4(i1,i2,i3)+sy(i1,i2,i3)*sxs4(
     & i1,i2,i3)+ty(i1,i2,i3)*sxt4(i1,i2,i3)
         sxz43(i1,i2,i3)=rz(i1,i2,i3)*sxr4(i1,i2,i3)+sz(i1,i2,i3)*sxs4(
     & i1,i2,i3)+tz(i1,i2,i3)*sxt4(i1,i2,i3)
         syx42(i1,i2,i3)= rx(i1,i2,i3)*syr4(i1,i2,i3)+sx(i1,i2,i3)*
     & sys4(i1,i2,i3)
         syy42(i1,i2,i3)= ry(i1,i2,i3)*syr4(i1,i2,i3)+sy(i1,i2,i3)*
     & sys4(i1,i2,i3)
         syx43(i1,i2,i3)=rx(i1,i2,i3)*syr4(i1,i2,i3)+sx(i1,i2,i3)*sys4(
     & i1,i2,i3)+tx(i1,i2,i3)*syt4(i1,i2,i3)
         syy43(i1,i2,i3)=ry(i1,i2,i3)*syr4(i1,i2,i3)+sy(i1,i2,i3)*sys4(
     & i1,i2,i3)+ty(i1,i2,i3)*syt4(i1,i2,i3)
         syz43(i1,i2,i3)=rz(i1,i2,i3)*syr4(i1,i2,i3)+sz(i1,i2,i3)*sys4(
     & i1,i2,i3)+tz(i1,i2,i3)*syt4(i1,i2,i3)
         szx42(i1,i2,i3)= rx(i1,i2,i3)*szr4(i1,i2,i3)+sx(i1,i2,i3)*
     & szs4(i1,i2,i3)
         szy42(i1,i2,i3)= ry(i1,i2,i3)*szr4(i1,i2,i3)+sy(i1,i2,i3)*
     & szs4(i1,i2,i3)
         szx43(i1,i2,i3)=rx(i1,i2,i3)*szr4(i1,i2,i3)+sx(i1,i2,i3)*szs4(
     & i1,i2,i3)+tx(i1,i2,i3)*szt4(i1,i2,i3)
         szy43(i1,i2,i3)=ry(i1,i2,i3)*szr4(i1,i2,i3)+sy(i1,i2,i3)*szs4(
     & i1,i2,i3)+ty(i1,i2,i3)*szt4(i1,i2,i3)
         szz43(i1,i2,i3)=rz(i1,i2,i3)*szr4(i1,i2,i3)+sz(i1,i2,i3)*szs4(
     & i1,i2,i3)+tz(i1,i2,i3)*szt4(i1,i2,i3)
         txx42(i1,i2,i3)= rx(i1,i2,i3)*txr4(i1,i2,i3)+sx(i1,i2,i3)*
     & txs4(i1,i2,i3)
         txy42(i1,i2,i3)= ry(i1,i2,i3)*txr4(i1,i2,i3)+sy(i1,i2,i3)*
     & txs4(i1,i2,i3)
         txx43(i1,i2,i3)=rx(i1,i2,i3)*txr4(i1,i2,i3)+sx(i1,i2,i3)*txs4(
     & i1,i2,i3)+tx(i1,i2,i3)*txt4(i1,i2,i3)
         txy43(i1,i2,i3)=ry(i1,i2,i3)*txr4(i1,i2,i3)+sy(i1,i2,i3)*txs4(
     & i1,i2,i3)+ty(i1,i2,i3)*txt4(i1,i2,i3)
         txz43(i1,i2,i3)=rz(i1,i2,i3)*txr4(i1,i2,i3)+sz(i1,i2,i3)*txs4(
     & i1,i2,i3)+tz(i1,i2,i3)*txt4(i1,i2,i3)
         tyx42(i1,i2,i3)= rx(i1,i2,i3)*tyr4(i1,i2,i3)+sx(i1,i2,i3)*
     & tys4(i1,i2,i3)
         tyy42(i1,i2,i3)= ry(i1,i2,i3)*tyr4(i1,i2,i3)+sy(i1,i2,i3)*
     & tys4(i1,i2,i3)
         tyx43(i1,i2,i3)=rx(i1,i2,i3)*tyr4(i1,i2,i3)+sx(i1,i2,i3)*tys4(
     & i1,i2,i3)+tx(i1,i2,i3)*tyt4(i1,i2,i3)
         tyy43(i1,i2,i3)=ry(i1,i2,i3)*tyr4(i1,i2,i3)+sy(i1,i2,i3)*tys4(
     & i1,i2,i3)+ty(i1,i2,i3)*tyt4(i1,i2,i3)
         tyz43(i1,i2,i3)=rz(i1,i2,i3)*tyr4(i1,i2,i3)+sz(i1,i2,i3)*tys4(
     & i1,i2,i3)+tz(i1,i2,i3)*tyt4(i1,i2,i3)
         tzx42(i1,i2,i3)= rx(i1,i2,i3)*tzr4(i1,i2,i3)+sx(i1,i2,i3)*
     & tzs4(i1,i2,i3)
         tzy42(i1,i2,i3)= ry(i1,i2,i3)*tzr4(i1,i2,i3)+sy(i1,i2,i3)*
     & tzs4(i1,i2,i3)
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
         uxy42(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*urr4(i1,i2,i3,kd)
     & +(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*urs4(i1,
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
         uxy43(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*urr4(i1,i2,i3,kd)
     & +sx(i1,i2,i3)*sy(i1,i2,i3)*uss4(i1,i2,i3,kd)+tx(i1,i2,i3)*ty(
     & i1,i2,i3)*utt4(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,
     & i2,i3)*sx(i1,i2,i3))*urs4(i1,i2,i3,kd)+(rx(i1,i2,i3)*ty(i1,i2,
     & i3)+ry(i1,i2,i3)*tx(i1,i2,i3))*urt4(i1,i2,i3,kd)+(sx(i1,i2,i3)*
     & ty(i1,i2,i3)+sy(i1,i2,i3)*tx(i1,i2,i3))*ust4(i1,i2,i3,kd)+
     & rxy43(i1,i2,i3)*ur4(i1,i2,i3,kd)+sxy43(i1,i2,i3)*us4(i1,i2,i3,
     & kd)+txy43(i1,i2,i3)*ut4(i1,i2,i3,kd)
         uxz43(i1,i2,i3,kd)=rx(i1,i2,i3)*rz(i1,i2,i3)*urr4(i1,i2,i3,kd)
     & +sx(i1,i2,i3)*sz(i1,i2,i3)*uss4(i1,i2,i3,kd)+tx(i1,i2,i3)*tz(
     & i1,i2,i3)*utt4(i1,i2,i3,kd)+(rx(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,
     & i2,i3)*sx(i1,i2,i3))*urs4(i1,i2,i3,kd)+(rx(i1,i2,i3)*tz(i1,i2,
     & i3)+rz(i1,i2,i3)*tx(i1,i2,i3))*urt4(i1,i2,i3,kd)+(sx(i1,i2,i3)*
     & tz(i1,i2,i3)+sz(i1,i2,i3)*tx(i1,i2,i3))*ust4(i1,i2,i3,kd)+
     & rxz43(i1,i2,i3)*ur4(i1,i2,i3,kd)+sxz43(i1,i2,i3)*us4(i1,i2,i3,
     & kd)+txz43(i1,i2,i3)*ut4(i1,i2,i3,kd)
         uyz43(i1,i2,i3,kd)=ry(i1,i2,i3)*rz(i1,i2,i3)*urr4(i1,i2,i3,kd)
     & +sy(i1,i2,i3)*sz(i1,i2,i3)*uss4(i1,i2,i3,kd)+ty(i1,i2,i3)*tz(
     & i1,i2,i3)*utt4(i1,i2,i3,kd)+(ry(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,
     & i2,i3)*sy(i1,i2,i3))*urs4(i1,i2,i3,kd)+(ry(i1,i2,i3)*tz(i1,i2,
     & i3)+rz(i1,i2,i3)*ty(i1,i2,i3))*urt4(i1,i2,i3,kd)+(sy(i1,i2,i3)*
     & tz(i1,i2,i3)+sz(i1,i2,i3)*ty(i1,i2,i3))*ust4(i1,i2,i3,kd)+
     & ryz43(i1,i2,i3)*ur4(i1,i2,i3,kd)+syz43(i1,i2,i3)*us4(i1,i2,i3,
     & kd)+tyz43(i1,i2,i3)*ut4(i1,i2,i3,kd)
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
         uxx43r(i1,i2,i3,kd)=( -30.*u(i1,i2,i3,kd)+16.*(u(i1+1,i2,i3,
     & kd)+u(i1-1,i2,i3,kd))-(u(i1+2,i2,i3,kd)+u(i1-2,i2,i3,kd)) )*
     & h42(0)
         uyy43r(i1,i2,i3,kd)=( -30.*u(i1,i2,i3,kd)+16.*(u(i1,i2+1,i3,
     & kd)+u(i1,i2-1,i3,kd))-(u(i1,i2+2,i3,kd)+u(i1,i2-2,i3,kd)) )*
     & h42(1)
         uzz43r(i1,i2,i3,kd)=( -30.*u(i1,i2,i3,kd)+16.*(u(i1,i2,i3+1,
     & kd)+u(i1,i2,i3-1,kd))-(u(i1,i2,i3+2,kd)+u(i1,i2,i3-2,kd)) )*
     & h42(2)
         uxy43r(i1,i2,i3,kd)=( (u(i1+2,i2+2,i3,kd)-u(i1-2,i2+2,i3,kd)- 
     & u(i1+2,i2-2,i3,kd)+u(i1-2,i2-2,i3,kd)) +8.*(u(i1-1,i2+2,i3,kd)-
     & u(i1-1,i2-2,i3,kd)-u(i1+1,i2+2,i3,kd)+u(i1+1,i2-2,i3,kd) +u(i1+
     & 2,i2-1,i3,kd)-u(i1-2,i2-1,i3,kd)-u(i1+2,i2+1,i3,kd)+u(i1-2,i2+
     & 1,i3,kd))+64.*(u(i1+1,i2+1,i3,kd)-u(i1-1,i2+1,i3,kd)- u(i1+1,
     & i2-1,i3,kd)+u(i1-1,i2-1,i3,kd)))*(h41(0)*h41(1))
         uxz43r(i1,i2,i3,kd)=( (u(i1+2,i2,i3+2,kd)-u(i1-2,i2,i3+2,kd)-
     & u(i1+2,i2,i3-2,kd)+u(i1-2,i2,i3-2,kd)) +8.*(u(i1-1,i2,i3+2,kd)-
     & u(i1-1,i2,i3-2,kd)-u(i1+1,i2,i3+2,kd)+u(i1+1,i2,i3-2,kd) +u(i1+
     & 2,i2,i3-1,kd)-u(i1-2,i2,i3-1,kd)- u(i1+2,i2,i3+1,kd)+u(i1-2,i2,
     & i3+1,kd)) +64.*(u(i1+1,i2,i3+1,kd)-u(i1-1,i2,i3+1,kd)-u(i1+1,
     & i2,i3-1,kd)+u(i1-1,i2,i3-1,kd)) )*(h41(0)*h41(2))
         uyz43r(i1,i2,i3,kd)=( (u(i1,i2+2,i3+2,kd)-u(i1,i2-2,i3+2,kd)-
     & u(i1,i2+2,i3-2,kd)+u(i1,i2-2,i3-2,kd)) +8.*(u(i1,i2-1,i3+2,kd)-
     & u(i1,i2-1,i3-2,kd)-u(i1,i2+1,i3+2,kd)+u(i1,i2+1,i3-2,kd) +u(i1,
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
         ulaplacian42r(i1,i2,i3,kd)=uxx43r(i1,i2,i3,kd)+uyy43r(i1,i2,
     & i3,kd)
         ulaplacian43r(i1,i2,i3,kd)=uxx43r(i1,i2,i3,kd)+uyy43r(i1,i2,
     & i3,kd)+uzz43r(i1,i2,i3,kd)
        !    --- For 2nd order 2D artificial diffusion ---
         delta22(c)=    (u(i1+1,i2,i3,c)-4.*u(i1,i2,i3,c)+u(i1-1,i2,i3,
     & c)  +u(i1,i2+1,i3,c)                 +u(i1,i2-1,i3,c))
        !    --- For 2nd order 3D artificial diffusion ---
         delta23(c)= (u(i1+1,i2,i3,c)-6.*u(i1,i2,i3,c)+u(i1-1,i2,i3,c) 
     &   +u(i1,i2+1,i3,c)                   +u(i1,i2-1,i3,c)  +u(i1,
     & i2,i3+1,c)                   +u(i1,i2,i3-1,c))
        !     ---For fourth-order artificial diffusion in 2D
         delta42(c)= (   -u(i1+2,i2,i3,c)-u(i1-2,i2,i3,c)   -u(i1,i2+2,
     & i3,c)-u(i1,i2-2,i3,c)   +4.*(u(i1+1,i2,i3,c)+u(i1-1,i2,i3,c)   
     & +u(i1,i2+1,i3,c)+u(i1,i2-1,i3,c))  -12.*u(i1,i2,i3,c) )
        !     ---For fourth-order artificial diffusion in 3D
         delta43(c)= (   -u(i1+2,i2,i3,c)-u(i1-2,i2,i3,c)  -u(i1,i2+2,
     & i3,c)-u(i1,i2-2,i3,c)  -u(i1,i2,i3+2,c)-u(i1,i2,i3-2,c)  +4.*(
     & u(i1+1,i2,i3,c)+u(i1-1,i2,i3,c)  +u(i1,i2+1,i3,c)+u(i1,i2-1,i3,
     & c)  +u(i1,i2,i3+1,c)+u(i1,i2,i3-1,c)) -18.*u(i1,i2,i3,c) )
        !     --- end statement functions
         ierr=0
         ! write(*,'("Inside insdt: gridType=",i2)') gridType
         pc                 =ipar(0)
         uc                 =ipar(1)
         vc                 =ipar(2)
         wc                 =ipar(3)
         nc                 =ipar(4)
         sc                 =ipar(5)
         tc                 =ipar(6)  ! **new**
         grid               =ipar(7)
         orderOfAccuracy    =ipar(8)
         gridIsMoving       =ipar(9)
         useWhereMask       =ipar(10)
         gridIsImplicit     =ipar(11)
         implicitMethod     =ipar(12)
         implicitOption     =ipar(13)
         isAxisymmetric     =ipar(14)
         use2ndOrderAD      =ipar(15)
         use4thOrderAD      =ipar(16)
         advectPassiveScalar=ipar(17)
         gridType           =ipar(18)
         turbulenceModel    =ipar(19)
         pdeModel           =ipar(20)
         vsc                =ipar(21)
         rc                 =ipar(22)
         debug              =ipar(23)
         materialFormat     =ipar(24)
         advectionOption    =ipar(25)  ! *new* 2017/01/27
         upwindOrder        =ipar(26)
         dr(0)             =rpar(0)
         dr(1)             =rpar(1)
         dr(2)             =rpar(2)
         dx(0)             =rpar(3)
         dx(1)             =rpar(4)
         dx(2)             =rpar(5)
         nu                =rpar(6)
         ad21              =rpar(7)
         ad22              =rpar(8)
         ad41              =rpar(9)
         ad42              =rpar(10)
         nuPassiveScalar   =rpar(11)
         adcPassiveScalar  =rpar(12)
         ad21n             =rpar(13)
         ad22n             =rpar(14)
         ad41n             =rpar(15)
         ad42n             =rpar(16)
        !       gravity(0)        =rpar(18)
        !      gravity(1)        =rpar(19)
        !      gravity(2)        =rpar(20)
        !      thermalExpansivity=rpar(21)
        !      adcBoussinesq     =rpar(22) ! coefficient of artificial diffusion for Boussinesq T equation 
        !      kThermal          =rpar(23)
         t                 =rpar(24)
        ! nuVP              =rpar(24)  ! for visco-plastic
         ! etaVP             =rpar(25)
         ! yieldStressVP     =rpar(26)
         ! exponentVP        =rpar(27)
         ! epsVP             =rpar(28)
        ! write(*,'(" insdt: eta,yield,exp,eps=",4e10.2)') etaVP,yieldStressVP,exponentVP,epsVP
         kc=nc
         ec=kc+1
         if( orderOfAccuracy.ne.2 .and. orderOfAccuracy.ne.4 )then
           write(*,'("insdt:ERROR orderOfAccuracy=",i6)') 
     & orderOfAccuracy
           stop 1
         end if
         if( gridType.ne.rectangular .and. gridType.ne.curvilinear )
     & then
           write(*,'("insdt:ERROR gridType=",i6)') gridType
           stop 2
         end if
         if( uc.lt.0 .or. vc.lt.0 .or. (nd.eq.3 .and. wc.lt.0) )then
           write(*,'("insdt:ERROR uc,vc,ws=",3i6)') uc,vc,wc
           stop 4
         end if
        !      write(*,'("insdt: turbulenceModel=",2i6)') turbulenceModel
        !      write(*,'("insdt: nd,uc,vc,wc,kc=",2i6)') nd,uc,vc,wc,kc
         if( turbulenceModel.eq.kEpsilon .and. (kc.lt.uc+nd .or. 
     & kc.gt.1000) )then
           write(*,'("insdt:ERROR in kc: nd,uc,vc,wc,kc=",2i6)') nd,uc,
     & vc,wc,kc
           stop 5
         end if
        ! ** these are needed by self-adjoint terms **fix**
         dxi=1./dx(0)
         dyi=1./dx(1)
         dzi=1./dx(2)
         dri=1./dr(0)
         dsi=1./dr(1)
         dti=1./dr(2)
         dr2i=1./(2.*dr(0))
         ds2i=1./(2.*dr(1))
         dt2i=1./(2.*dr(2))
         if( turbulenceModel.eq.spalartAllmaras )then
           call getSpalartAllmarasParameters(cb1, cb2, cv1, sigma, 
     & sigmai, kappa, cw1, cw2, cw3, cw3e6, cv1e3, cd0, cr0)
         else if( turbulenceModel.eq.kEpsilon )then
          ! write(*,'(" insdt: k-epsilon: nc,kc,ec=",3i3)') nc,kc,ec
           call getKEpsilonParameters( cMu,cEps1,cEps2,sigmaEpsI,
     & sigmaKI )
           !  write(*,'(" insdt: cMu,cEps1,cEps2,sigmaEpsI,sigmaKI=",5f8.3)') cMu,cEps1,cEps2,sigmaEpsI,sigmaKI
         else if( turbulenceModel.eq.largeEddySimulation )then
           ! do nothing
         else if( turbulenceModel.ne.noTurbulenceModel )then
           stop 88
         end if
         adc=adcPassiveScalar ! coefficient of linear artificial diffusion
         cd22=ad22/(nd**2)
         cd42=ad42/(nd**2)
        !     *********************************      
        !     ********MAIN LOOPS***************      
        !     *********************************      
         if( gridType.eq.rectangular )then
          if( isAxisymmetric.eq.0 )then
           if( gridIsImplicit.eq.0 )then
            ! --- explicit time-stepping ---
            if( advectionOption.eq.centeredAdvection )then
              if( useWhereMask.ne.0 )then
               do i3=n3a,n3b
               do i2=n2a,n2b
               do i1=n1a,n1b
                if( mask(i1,i2,i3).gt.0 )then
                   ! INS, no AD
                    ut(i1,i2,i3,uc)= -uu(i1,i2,i3,uc)*ux42r(i1,i2,i3,
     & uc)-uu(i1,i2,i3,vc)*uy42r(i1,i2,i3,uc)-ux42r(i1,i2,i3,pc)+nu*
     & ulaplacian42r(i1,i2,i3,uc)
                    ut(i1,i2,i3,vc)= -uu(i1,i2,i3,uc)*ux42r(i1,i2,i3,
     & vc)-uu(i1,i2,i3,vc)*uy42r(i1,i2,i3,vc)-uy42r(i1,i2,i3,pc)+nu*
     & ulaplacian42r(i1,i2,i3,vc)
                end if
               end do
               end do
               end do
              else
               do i3=n3a,n3b
               do i2=n2a,n2b
               do i1=n1a,n1b
                  ! INS, no AD
                   ut(i1,i2,i3,uc)= -uu(i1,i2,i3,uc)*ux42r(i1,i2,i3,uc)
     & -uu(i1,i2,i3,vc)*uy42r(i1,i2,i3,uc)-ux42r(i1,i2,i3,pc)+nu*
     & ulaplacian42r(i1,i2,i3,uc)
                   ut(i1,i2,i3,vc)= -uu(i1,i2,i3,uc)*ux42r(i1,i2,i3,vc)
     & -uu(i1,i2,i3,vc)*uy42r(i1,i2,i3,vc)-uy42r(i1,i2,i3,pc)+nu*
     & ulaplacian42r(i1,i2,i3,vc)
               end do
               end do
               end do
              end if
            else if( advectionOption.eq.upwindAdvection )then
              ! --- upwind ---
              if( useWhereMask.ne.0 )then
               do i3=n3a,n3b
               do i2=n2a,n2b
               do i1=n1a,n1b
                if( mask(i1,i2,i3).gt.0 )then
                  ! INS, no AD
                    ! --- upwind approximations ---
                      ! --- upwind scheme ---
                      ! for testing output this next message:
                      ! if( t.le. 0. )then
                      !    write(*,'(" getAdvection upwind scheme (7)")') 
                      ! end if
                        ! --- CARTESIAN GRID ---
                        !- agu(uc,uc)=UU(uc)*UX(uc)
                        !- agu(vc,uc)=UU(vc)*UY(uc)
                        !- agu(uc,vc)=UU(uc)*UX(vc)
                        !- agu(vc,vc)=UU(vc)*UY(vc)
                        ! -- first order upwind --
                        if( upwindOrder.eq.1 )then
                         au = uu(i1,i2,i3,uc)
                         if( au.gt.0. )then
                           ! u*ux = u*D-x(u)
                           agu(uc,uc)= au*(u(i1,i2,i3,uc)-u(i1-1,i2,i3,
     & uc))/(dx(0))
                           ! u*vx = u*D-x(v)
                           agu(uc,vc)= au*(u(i1,i2,i3,vc)-u(i1-1,i2,i3,
     & vc))/(dx(0))
                         else
                           ! u*ux = u*D+x(u)
                           agu(uc,uc)= au*(u(i1+1,i2,i3,uc)-u(i1,i2,i3,
     & uc))/(dx(0))
                           ! u*vx = u*D+x(v)
                           agu(uc,vc)= au*(u(i1+1,i2,i3,vc)-u(i1,i2,i3,
     & vc))/(dx(0))
                         end if
                         au = uu(i1,i2,i3,vc)
                         if( au.gt.0. )then
                           ! v*uy = v*D-y(u)
                           agu(vc,uc)= au*(u(i1,i2,i3,uc)-u(i1,i2-1,i3,
     & uc))/(dx(1))
                           ! v*vy = v*D-y(v)
                           agu(vc,vc)= au*(u(i1,i2,i3,vc)-u(i1,i2-1,i3,
     & vc))/(dx(1))
                         else
                           ! v*uy = v*D+y(u)
                           agu(vc,uc)= au*(u(i1,i2+1,i3,uc)-u(i1,i2,i3,
     & uc))/(dx(1))
                           ! v*vy = v*D+y(v) 
                           agu(vc,vc)= au*(u(i1,i2+1,i3,vc)-u(i1,i2,i3,
     & vc))/(dx(1))
                         end if
                        elseif( upwindOrder.eq.2 )then
                          ! write(*,'(" finish me, upwindOrder=",i2)') upwindOrder
                          ! stop 222
                           au = uu(i1,i2,i3,uc)
                           if( au.gt.0. )then
                              ! u*ux = u*D-x(u)
                              agu(uc,uc)= au*(3.*u(i1,i2,i3,uc)-4.*u(
     & i1-1,i2,i3,uc)+u(i1-2,i2,i3,uc))/(2.*dx(0))
                              ! u*vx = u*D-x(v)
                              agu(uc,vc)= au*(3.*u(i1,i2,i3,vc)-4.*u(
     & i1-1,i2,i3,vc)+u(i1-2,i2,i3,vc))/(2.*dx(0))
                           else
                              ! u*ux = u*D+x(u)
                              agu(uc,uc)= au*(-u(i1+2,i2,i3,uc)+4.*u(
     & i1+1,i2,i3,uc)-3.*u(i1,i2,i3,uc))/(2.*dx(0))
                              ! u*vx = u*D+x(v)
                              agu(uc,vc)= au*(-u(i1+2,i2,i3,vc)+4.*u(
     & i1+1,i2,i3,vc)-3.*u(i1,i2,i3,vc))/(2.*dx(0))
                           end if
                           au = uu(i1,i2,i3,vc)
                           if( au.gt.0. )then
                              ! v*uy = v*D-y(u)
                              agu(vc,uc)= au*(3.*u(i1,i2,i3,uc)-4.*u(
     & i1,i2-1,i3,uc)+u(i1,i2-2,i3,uc))/(2.*dx(1))
                              ! v*vy = v*D-y(v)
                              agu(vc,vc)= au*(3.*u(i1,i2,i3,vc)-4.*u(
     & i1,i2-1,i3,vc)+u(i1,i2-2,i3,vc))/(2.*dx(1))
                           else
                              ! v*uy = v*D+y(u)
                              agu(vc,uc)= au*(-u(i1,i2+2,i3,uc)+4.*u(
     & i1,i2+1,i3,uc)-3.*u(i1,i2,i3,uc))/(2.*dx(1))
                              ! v*vy = v*D+y(v) 
                              agu(vc,vc)= au*(-u(i1,i2+2,i3,vc)+4.*u(
     & i1,i2+1,i3,vc)-3.*u(i1,i2,i3,vc))/(2.*dx(1))
                           end if
                        elseif( upwindOrder.eq.3 )then
                           au = uu(i1,i2,i3,uc)
                           if( au.gt.0. )then
                              ! u*ux = u*D-x(u)
                              agu(uc,uc)= au*(2.*u(i1+1,i2,i3,uc)+3.*u(
     & i1,i2,i3,uc)-6.*u(i1-1,i2,i3,uc)+u(i1-2,i2,i3,uc))/(6.*dx(0))
                              ! u*vx = u*D-x(v)
                              agu(uc,vc)= au*(2.*u(i1+1,i2,i3,vc)+3.*u(
     & i1,i2,i3,vc)-6.*u(i1-1,i2,i3,vc)+u(i1-2,i2,i3,vc))/(6.*dx(0))
                           else
                              ! u*ux = u*D+x(u)
                              agu(uc,uc)= au*(-u(i1+2,i2,i3,uc)+6.*u(
     & i1+1,i2,i3,uc)-3.*u(i1,i2,i3,uc)-2.*u(i1-1,i2,i3,uc))/(6.*dx(0)
     & )
                              ! u*vx = u*D+x(v)
                              agu(uc,vc)= au*(-u(i1+2,i2,i3,vc)+6.*u(
     & i1+1,i2,i3,vc)-3.*u(i1,i2,i3,vc)-2.*u(i1-1,i2,i3,vc))/(6.*dx(0)
     & )
                           end if
                           au = uu(i1,i2,i3,vc)
                           if( au.gt.0. )then
                              ! v*uy = v*D-y(u)
                              agu(vc,uc)= au*(2.*u(i1,i2+1,i3,uc)+3.*u(
     & i1,i2,i3,uc)-6.*u(i1,i2-1,i3,uc)+u(i1,i2-2,i3,uc))/(6.*dx(1))
                              ! v*vy = v*D-y(v)
                              agu(vc,vc)= au*(2.*u(i1,i2+1,i3,vc)+3.*u(
     & i1,i2,i3,vc)-6.*u(i1,i2-1,i3,vc)+u(i1,i2-2,i3,vc))/(6.*dx(1))
                           else
                              ! v*uy = v*D+y(u)
                              agu(vc,uc)= au*(-u(i1,i2+2,i3,uc)+6.*u(
     & i1,i2+1,i3,uc)-3.*u(i1,i2,i3,uc)-2.*u(i1,i2-1,i3,uc))/(6.*dx(1)
     & )
                              ! v*vy = v*D+y(v) 
                              agu(vc,vc)= au*(-u(i1,i2+2,i3,vc)+6.*u(
     & i1,i2+1,i3,vc)-3.*u(i1,i2,i3,vc)-2.*u(i1,i2-1,i3,vc))/(6.*dx(1)
     & )
                           end if
                        end if
                    ! #If ("UPWIND" == "CENTERED") 
                    !  ! -- centered advection ---
                    !  ! write(*,'(" getAdvection -- centered")')
                    !  #If "2" == "2"
                    !   ! -- 2D --
                    !   agu(uc,uc)=UU(uc)*UX(uc)
                    !   agu(vc,uc)=UU(vc)*UY(uc)
                    !   agu(uc,vc)=UU(uc)*UX(vc)
                    !   agu(vc,vc)=UU(vc)*UY(vc)
                    !  #Else
                    !   ! -- 3D -- *check me*
                    !   agu(uc,uc)=UU(uc)*UX(uc)
                    !   agu(vc,uc)=UU(vc)*UY(uc)
                    !   agu(wc,uc)=UU(wc)*UY(uc)
                    !   agu(uc,vc)=UU(uc)*UX(vc)
                    !   agu(vc,vc)=UU(vc)*UY(vc)
                    !   agu(wc,vc)=UU(wc)*UY(vc)
                    !   agu(uc,wc)=UU(uc)*UX(wc)
                    !   agu(vc,wc)=UU(vc)*UY(wc)
                    !   agu(wc,wc)=UU(wc)*UY(wc)
                    !  #End
                    ! #Elif "UPWIND" == "UPWIND" 
                    !   ! --- upwind scheme ---
                    !   ! for testing output this next message:
                    !   if( t.le. 0. )then
                    !     write(*,'(" getAdvection upwind scheme (7)")') 
                    !   end if
                    !   getUpwindAdvection(u,i1,i2,i3,NONE,2,4,rectangular, agu)
                    ! #Elif "UPWIND" == "BWENO" 
                    !   ! --- Bweno scheme ---
                    !   ! for testing output this message:
                    !   if( t.le. 0. )then
                    !      write(*,'(" getAdvection BWENO scheme (7)")') 
                    !   end if
                    !   getBwenoAdvection(u,i1,i2,i3,NONE,2,4,rectangular, agu)
                    ! #Else
                    !   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
                    !   stop 999
                    ! #End
                       ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-ux42r(
     & i1,i2,i3,pc)+nu*ulaplacian42r(i1,i2,i3,uc)
                       ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-uy42r(
     & i1,i2,i3,pc)+nu*ulaplacian42r(i1,i2,i3,vc)
                 !       u(i1,i2,i3,3)= agu(uc,3) ! au in r direction 
                 !       u(i1,i2,i3,4)= agu(uc,8) ! au in s direction
                 !       u(i1,i2,i3,5)= agu(uc,4)  ! wpl uur
                 !       u(i1,i2,i3,6)= agu(uc,5)  ! wpr uur
                 !       u(i1,i2,i3,7)= agu(uc,6)  ! wml uur
                 !       u(i1,i2,i3,8)= agu(uc,7)  ! wmr uur
                 !       u(i1,i2,i3,9)= agu(uc,9)  ! wpl uus
                 !       u(i1,i2,i3,10)= agu(uc,10)  ! wpr uus
                 !       u(i1,i2,i3,11)= agu(uc,11)  ! wml uus
                 !       u(i1,i2,i3,12)= agu(uc,12)  ! wmr uus
                 ! !
                 !       u(i1,i2,i3,13)= agu(vc,4)  ! wpl uvr
                 !       u(i1,i2,i3,14)= agu(vc,5)  ! wpr uvr
                 !       u(i1,i2,i3,15)= agu(vc,6)  ! wml uvr
                 !       u(i1,i2,i3,16)= agu(vc,7)  ! wmr uvr
                 !       u(i1,i2,i3,17)= agu(vc,9)  ! wpl uvs
                 !       u(i1,i2,i3,18)= agu(vc,10)  ! wpr uvs
                 !       u(i1,i2,i3,19)= agu(vc,11)  ! wml uvs
                 !       u(i1,i2,i3,20)= agu(vc,12)  ! wmr uvs
                 !       u(i1,i2,i3,21)= agu(uc,13)  ! au at where fix is turned on 
                 !       u(i1,i2,i3,22)= agu(vc,13)  ! au at where fix is turned on 
                end if
               end do
               end do
               end do
              else
               do i3=n3a,n3b
               do i2=n2a,n2b
               do i1=n1a,n1b
                 ! INS, no AD
                   ! --- upwind approximations ---
                     ! --- upwind scheme ---
                     ! for testing output this next message:
                     ! if( t.le. 0. )then
                     !    write(*,'(" getAdvection upwind scheme (7)")') 
                     ! end if
                       ! --- CARTESIAN GRID ---
                       !- agu(uc,uc)=UU(uc)*UX(uc)
                       !- agu(vc,uc)=UU(vc)*UY(uc)
                       !- agu(uc,vc)=UU(uc)*UX(vc)
                       !- agu(vc,vc)=UU(vc)*UY(vc)
                       ! -- first order upwind --
                       if( upwindOrder.eq.1 )then
                        au = uu(i1,i2,i3,uc)
                        if( au.gt.0. )then
                          ! u*ux = u*D-x(u)
                          agu(uc,uc)= au*(u(i1,i2,i3,uc)-u(i1-1,i2,i3,
     & uc))/(dx(0))
                          ! u*vx = u*D-x(v)
                          agu(uc,vc)= au*(u(i1,i2,i3,vc)-u(i1-1,i2,i3,
     & vc))/(dx(0))
                        else
                          ! u*ux = u*D+x(u)
                          agu(uc,uc)= au*(u(i1+1,i2,i3,uc)-u(i1,i2,i3,
     & uc))/(dx(0))
                          ! u*vx = u*D+x(v)
                          agu(uc,vc)= au*(u(i1+1,i2,i3,vc)-u(i1,i2,i3,
     & vc))/(dx(0))
                        end if
                        au = uu(i1,i2,i3,vc)
                        if( au.gt.0. )then
                          ! v*uy = v*D-y(u)
                          agu(vc,uc)= au*(u(i1,i2,i3,uc)-u(i1,i2-1,i3,
     & uc))/(dx(1))
                          ! v*vy = v*D-y(v)
                          agu(vc,vc)= au*(u(i1,i2,i3,vc)-u(i1,i2-1,i3,
     & vc))/(dx(1))
                        else
                          ! v*uy = v*D+y(u)
                          agu(vc,uc)= au*(u(i1,i2+1,i3,uc)-u(i1,i2,i3,
     & uc))/(dx(1))
                          ! v*vy = v*D+y(v) 
                          agu(vc,vc)= au*(u(i1,i2+1,i3,vc)-u(i1,i2,i3,
     & vc))/(dx(1))
                        end if
                       elseif( upwindOrder.eq.2 )then
                         ! write(*,'(" finish me, upwindOrder=",i2)') upwindOrder
                         ! stop 222
                          au = uu(i1,i2,i3,uc)
                          if( au.gt.0. )then
                             ! u*ux = u*D-x(u)
                             agu(uc,uc)= au*(3.*u(i1,i2,i3,uc)-4.*u(i1-
     & 1,i2,i3,uc)+u(i1-2,i2,i3,uc))/(2.*dx(0))
                             ! u*vx = u*D-x(v)
                             agu(uc,vc)= au*(3.*u(i1,i2,i3,vc)-4.*u(i1-
     & 1,i2,i3,vc)+u(i1-2,i2,i3,vc))/(2.*dx(0))
                          else
                             ! u*ux = u*D+x(u)
                             agu(uc,uc)= au*(-u(i1+2,i2,i3,uc)+4.*u(i1+
     & 1,i2,i3,uc)-3.*u(i1,i2,i3,uc))/(2.*dx(0))
                             ! u*vx = u*D+x(v)
                             agu(uc,vc)= au*(-u(i1+2,i2,i3,vc)+4.*u(i1+
     & 1,i2,i3,vc)-3.*u(i1,i2,i3,vc))/(2.*dx(0))
                          end if
                          au = uu(i1,i2,i3,vc)
                          if( au.gt.0. )then
                             ! v*uy = v*D-y(u)
                             agu(vc,uc)= au*(3.*u(i1,i2,i3,uc)-4.*u(i1,
     & i2-1,i3,uc)+u(i1,i2-2,i3,uc))/(2.*dx(1))
                             ! v*vy = v*D-y(v)
                             agu(vc,vc)= au*(3.*u(i1,i2,i3,vc)-4.*u(i1,
     & i2-1,i3,vc)+u(i1,i2-2,i3,vc))/(2.*dx(1))
                          else
                             ! v*uy = v*D+y(u)
                             agu(vc,uc)= au*(-u(i1,i2+2,i3,uc)+4.*u(i1,
     & i2+1,i3,uc)-3.*u(i1,i2,i3,uc))/(2.*dx(1))
                             ! v*vy = v*D+y(v) 
                             agu(vc,vc)= au*(-u(i1,i2+2,i3,vc)+4.*u(i1,
     & i2+1,i3,vc)-3.*u(i1,i2,i3,vc))/(2.*dx(1))
                          end if
                       elseif( upwindOrder.eq.3 )then
                          au = uu(i1,i2,i3,uc)
                          if( au.gt.0. )then
                             ! u*ux = u*D-x(u)
                             agu(uc,uc)= au*(2.*u(i1+1,i2,i3,uc)+3.*u(
     & i1,i2,i3,uc)-6.*u(i1-1,i2,i3,uc)+u(i1-2,i2,i3,uc))/(6.*dx(0))
                             ! u*vx = u*D-x(v)
                             agu(uc,vc)= au*(2.*u(i1+1,i2,i3,vc)+3.*u(
     & i1,i2,i3,vc)-6.*u(i1-1,i2,i3,vc)+u(i1-2,i2,i3,vc))/(6.*dx(0))
                          else
                             ! u*ux = u*D+x(u)
                             agu(uc,uc)= au*(-u(i1+2,i2,i3,uc)+6.*u(i1+
     & 1,i2,i3,uc)-3.*u(i1,i2,i3,uc)-2.*u(i1-1,i2,i3,uc))/(6.*dx(0))
                             ! u*vx = u*D+x(v)
                             agu(uc,vc)= au*(-u(i1+2,i2,i3,vc)+6.*u(i1+
     & 1,i2,i3,vc)-3.*u(i1,i2,i3,vc)-2.*u(i1-1,i2,i3,vc))/(6.*dx(0))
                          end if
                          au = uu(i1,i2,i3,vc)
                          if( au.gt.0. )then
                             ! v*uy = v*D-y(u)
                             agu(vc,uc)= au*(2.*u(i1,i2+1,i3,uc)+3.*u(
     & i1,i2,i3,uc)-6.*u(i1,i2-1,i3,uc)+u(i1,i2-2,i3,uc))/(6.*dx(1))
                             ! v*vy = v*D-y(v)
                             agu(vc,vc)= au*(2.*u(i1,i2+1,i3,vc)+3.*u(
     & i1,i2,i3,vc)-6.*u(i1,i2-1,i3,vc)+u(i1,i2-2,i3,vc))/(6.*dx(1))
                          else
                             ! v*uy = v*D+y(u)
                             agu(vc,uc)= au*(-u(i1,i2+2,i3,uc)+6.*u(i1,
     & i2+1,i3,uc)-3.*u(i1,i2,i3,uc)-2.*u(i1,i2-1,i3,uc))/(6.*dx(1))
                             ! v*vy = v*D+y(v) 
                             agu(vc,vc)= au*(-u(i1,i2+2,i3,vc)+6.*u(i1,
     & i2+1,i3,vc)-3.*u(i1,i2,i3,vc)-2.*u(i1,i2-1,i3,vc))/(6.*dx(1))
                          end if
                       end if
                   ! #If ("UPWIND" == "CENTERED") 
                   !  ! -- centered advection ---
                   !  ! write(*,'(" getAdvection -- centered")')
                   !  #If "2" == "2"
                   !   ! -- 2D --
                   !   agu(uc,uc)=UU(uc)*UX(uc)
                   !   agu(vc,uc)=UU(vc)*UY(uc)
                   !   agu(uc,vc)=UU(uc)*UX(vc)
                   !   agu(vc,vc)=UU(vc)*UY(vc)
                   !  #Else
                   !   ! -- 3D -- *check me*
                   !   agu(uc,uc)=UU(uc)*UX(uc)
                   !   agu(vc,uc)=UU(vc)*UY(uc)
                   !   agu(wc,uc)=UU(wc)*UY(uc)
                   !   agu(uc,vc)=UU(uc)*UX(vc)
                   !   agu(vc,vc)=UU(vc)*UY(vc)
                   !   agu(wc,vc)=UU(wc)*UY(vc)
                   !   agu(uc,wc)=UU(uc)*UX(wc)
                   !   agu(vc,wc)=UU(vc)*UY(wc)
                   !   agu(wc,wc)=UU(wc)*UY(wc)
                   !  #End
                   ! #Elif "UPWIND" == "UPWIND" 
                   !   ! --- upwind scheme ---
                   !   ! for testing output this next message:
                   !   if( t.le. 0. )then
                   !     write(*,'(" getAdvection upwind scheme (7)")') 
                   !   end if
                   !   getUpwindAdvection(u,i1,i2,i3,NONE,2,4,rectangular, agu)
                   ! #Elif "UPWIND" == "BWENO" 
                   !   ! --- Bweno scheme ---
                   !   ! for testing output this message:
                   !   if( t.le. 0. )then
                   !      write(*,'(" getAdvection BWENO scheme (7)")') 
                   !   end if
                   !   getBwenoAdvection(u,i1,i2,i3,NONE,2,4,rectangular, agu)
                   ! #Else
                   !   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
                   !   stop 999
                   ! #End
                      ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-ux42r(
     & i1,i2,i3,pc)+nu*ulaplacian42r(i1,i2,i3,uc)
                      ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-uy42r(
     & i1,i2,i3,pc)+nu*ulaplacian42r(i1,i2,i3,vc)
                !       u(i1,i2,i3,3)= agu(uc,3) ! au in r direction 
                !       u(i1,i2,i3,4)= agu(uc,8) ! au in s direction
                !       u(i1,i2,i3,5)= agu(uc,4)  ! wpl uur
                !       u(i1,i2,i3,6)= agu(uc,5)  ! wpr uur
                !       u(i1,i2,i3,7)= agu(uc,6)  ! wml uur
                !       u(i1,i2,i3,8)= agu(uc,7)  ! wmr uur
                !       u(i1,i2,i3,9)= agu(uc,9)  ! wpl uus
                !       u(i1,i2,i3,10)= agu(uc,10)  ! wpr uus
                !       u(i1,i2,i3,11)= agu(uc,11)  ! wml uus
                !       u(i1,i2,i3,12)= agu(uc,12)  ! wmr uus
                ! !
                !       u(i1,i2,i3,13)= agu(vc,4)  ! wpl uvr
                !       u(i1,i2,i3,14)= agu(vc,5)  ! wpr uvr
                !       u(i1,i2,i3,15)= agu(vc,6)  ! wml uvr
                !       u(i1,i2,i3,16)= agu(vc,7)  ! wmr uvr
                !       u(i1,i2,i3,17)= agu(vc,9)  ! wpl uvs
                !       u(i1,i2,i3,18)= agu(vc,10)  ! wpr uvs
                !       u(i1,i2,i3,19)= agu(vc,11)  ! wml uvs
                !       u(i1,i2,i3,20)= agu(vc,12)  ! wmr uvs
                !       u(i1,i2,i3,21)= agu(uc,13)  ! au at where fix is turned on 
                !       u(i1,i2,i3,22)= agu(vc,13)  ! au at where fix is turned on 
               end do
               end do
               end do
              end if
            else if( advectionOption.eq.bwenoAdvection )then
              ! --- bweno ---
              if( useWhereMask.ne.0 )then
               do i3=n3a,n3b
               do i2=n2a,n2b
               do i1=n1a,n1b
                if( mask(i1,i2,i3).gt.0 )then
                  ! INS, no AD
                    ! --- upwind approximations ---
                      ! --- Bweno scheme ---
                      ! for testing output this message:
                      ! if( t.le. 0. )then
                      !    write(*,'(" getAdvection BWENO scheme (7)")') 
                      ! end if
                         gvU = u(i1,i2,i3,uc) - uu(i1,i2,i3,uc)
                         gvV = u(i1,i2,i3,vc) - uu(i1,i2,i3,vc)
                         gvW = u(i1,i2,i3,wc) - uu(i1,i2,i3,wc)
                         orderOneFix = 1
                         augmentPlot = 1
                         if( upwindOrder.eq.4 )then
                        ! -- fourth order bweno --
                            au = uu(i1,i2,i3,uc)
                            aum = u(i1-1,i2,i3,uc) - gvU
                            aup = u(i1+1,i2,i3,uc) - gvU
                            !usage
                            !s1   is the shifted grid in r1 direction
                            !var  is the variable taking derivative (uc,vc,wc)
                            !vard is the variable of direction
                            !agu is only passed in for augmented plot
                            !expAd 
                            ! u*ux
                            s1=1
                            s2=0
                            s3=0
                            var = uc
                            vard= uc
                            drl = dx(0)
                              Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*u(
     & i1   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
                              Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    u(
     & i1-s1,i2-s2,i3-s3,var)
                              Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) - 
     & 2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
                              Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) - 
     &    u(i1   ,i2   ,i3   ,var)
                              Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,
     & i2-s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
                              Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var)
                              Amr   = Apl
                              Bmr   = Bpl
                              betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 
     & 1./4.*Bpl**2
                              betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 
     & 1./4.*Bpr**2
                              betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 
     & 1./4.*Bml**2
                              betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 
     & 1./4.*Bmr**2
                              ep     = 1e-15
                              ap1    = (1./2.)/(ep + betapl)**(4)
                              ap2    = (1./2.)/(ep + betapr)**(4)
                              am1    = (1./2.)/(ep + betaml)**(4)
                              am2    = (1./2.)/(ep + betamr)**(4)
                              !calculate weights
                              wp1 = ap1/(ap1+ap2)
                              wp2 = ap2/(ap1+ap2)
                              wm1 = am1/(am1+am2)
                              wm2 = am2/(am1+am2)
                              !mapping of the weights
                              bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                              bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                              bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                              bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                              wplw = bp1/(bp1+bp2)
                              wprw = bp2/(bp1+bp2)
                              wmlw = bm1/(bm1+bm2)
                              wmrw = bm2/(bm1+bm2)
                              !get face value of u
                              up = au
                              um = au
                              !picking side for the weights
                              if (up.ge.0) then
                                 wpl = max(wplw,wprw)
                                 wpr = min(wplw,wprw)
                              else
                                 wpl = min(wplw,wprw)
                                 wpr = max(wplw,wprw)
                              end if
                              if (um.ge.0) then
                                 wml = max(wmlw,wmrw)
                                 wmr = min(wmlw,wmrw)
                              else
                                 wml = min(wmlw,wmrw)
                                 wmr = max(wmlw,wmrw)
                              end if
                              expAd = 0.
                              if (orderOneFix.eq.1) then
                                 ! if (augmentPlot.eq.1) then
                                 !    agu(uc,13)=0
                                 ! end if
                                 if (aum*aup .le. 0) then
                                    aumax = max(abs(au),abs(aup),abs(
     & aum))
                                    expAd = (1./12.)*aumax*( u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(i1,
     & i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*s2,
     & i3+2*s3,var))
                                    wpl = 0.5
                                    wpr = 0.5
                                    wml = 0.5
                                    wmr = 0.5
                                    ! plot au at where fix is turned on
                                    ! if (augmentPlot.eq.1) then
                                    !    agu(uc,13)=au
                                    ! end if
                                 end if
                              end if
                              ! wpl = 0.5
                              ! wpr = 0.5
                              ! wml = 0.5
                              ! wmr = 0.5
                              Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,var) 
     & + 5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
                              Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 5.*
     & u(i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
                              Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*
     & s3,var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
                              Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,var)
     &  + 5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
                              Fp  = wpl*Fpl + wpr*Fpr
                              Fm  = wml*Fml + wmr*Fmr
                                 agu(vard,var)= au*(Fp - Fm)/drl + 
     & expAd/drl
                            !     if ( i1 .eq. 14 .and. i2.eq.16 ) then
                                    !     if ( i1 .eq. 14) then
                                    ! write(*,'(" betapl,betapr,betaml,betamr = ",4F10.8)')  betapl,betapr,betaml,betamr
                                    ! write(*,'(" ap1,ap2,am1,am2 = ",4e10.8)')  ap1,ap2,am1,am2
                                    ! write(*,'(" wp1,wpr,wm1,wm2 = ",4F10.8)')  wp1,wp2,wm1,wm2
                                    ! write(*,'(" bp1,bp2,bm1,bm2 = ",4F10.8)')  bp1,bp2,bm1,bm2
                                    ! write(*,'(" wplw,wprw,wmlw,wmrw = ",4F10.8)')  wplw,wprw,wmlw,wmrw
                                    ! write(*,'(" wpl,wpr,wml,wmr = ",4F10.8)')  wpl,wpr,wml,wmr
                                    ! write(*,'(" agu = ",F10.8)') agu(vard,var)
                            !     end if
                            !u*vx
                            s1=1
                            s2=0
                            s3=0
                            var = vc
                            vard= uc
                            drl = dx(0)
                              Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*u(
     & i1   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
                              Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    u(
     & i1-s1,i2-s2,i3-s3,var)
                              Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) - 
     & 2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
                              Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) - 
     &    u(i1   ,i2   ,i3   ,var)
                              Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,
     & i2-s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
                              Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var)
                              Amr   = Apl
                              Bmr   = Bpl
                              betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 
     & 1./4.*Bpl**2
                              betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 
     & 1./4.*Bpr**2
                              betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 
     & 1./4.*Bml**2
                              betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 
     & 1./4.*Bmr**2
                              ep     = 1e-15
                              ap1    = (1./2.)/(ep + betapl)**(4)
                              ap2    = (1./2.)/(ep + betapr)**(4)
                              am1    = (1./2.)/(ep + betaml)**(4)
                              am2    = (1./2.)/(ep + betamr)**(4)
                              !calculate weights
                              wp1 = ap1/(ap1+ap2)
                              wp2 = ap2/(ap1+ap2)
                              wm1 = am1/(am1+am2)
                              wm2 = am2/(am1+am2)
                              !mapping of the weights
                              bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                              bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                              bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                              bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                              wplw = bp1/(bp1+bp2)
                              wprw = bp2/(bp1+bp2)
                              wmlw = bm1/(bm1+bm2)
                              wmrw = bm2/(bm1+bm2)
                              !get face value of u
                              up = au
                              um = au
                              !picking side for the weights
                              if (up.ge.0) then
                                 wpl = max(wplw,wprw)
                                 wpr = min(wplw,wprw)
                              else
                                 wpl = min(wplw,wprw)
                                 wpr = max(wplw,wprw)
                              end if
                              if (um.ge.0) then
                                 wml = max(wmlw,wmrw)
                                 wmr = min(wmlw,wmrw)
                              else
                                 wml = min(wmlw,wmrw)
                                 wmr = max(wmlw,wmrw)
                              end if
                              expAd = 0.
                              if (orderOneFix.eq.1) then
                                 ! if (augmentPlot.eq.1) then
                                 !    agu(uc,13)=0
                                 ! end if
                                 if (aum*aup .le. 0) then
                                    aumax = max(abs(au),abs(aup),abs(
     & aum))
                                    expAd = (1./12.)*aumax*( u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(i1,
     & i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*s2,
     & i3+2*s3,var))
                                    wpl = 0.5
                                    wpr = 0.5
                                    wml = 0.5
                                    wmr = 0.5
                                    ! plot au at where fix is turned on
                                    ! if (augmentPlot.eq.1) then
                                    !    agu(uc,13)=au
                                    ! end if
                                 end if
                              end if
                              ! wpl = 0.5
                              ! wpr = 0.5
                              ! wml = 0.5
                              ! wmr = 0.5
                              Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,var) 
     & + 5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
                              Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 5.*
     & u(i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
                              Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*
     & s3,var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
                              Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,var)
     &  + 5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
                              Fp  = wpl*Fpl + wpr*Fpr
                              Fm  = wml*Fml + wmr*Fmr
                                 agu(vard,var)= au*(Fp - Fm)/drl + 
     & expAd/drl
                            !     if ( i1 .eq. 14 .and. i2.eq.16 ) then
                                    !     if ( i1 .eq. 14) then
                                    ! write(*,'(" betapl,betapr,betaml,betamr = ",4F10.8)')  betapl,betapr,betaml,betamr
                                    ! write(*,'(" ap1,ap2,am1,am2 = ",4e10.8)')  ap1,ap2,am1,am2
                                    ! write(*,'(" wp1,wpr,wm1,wm2 = ",4F10.8)')  wp1,wp2,wm1,wm2
                                    ! write(*,'(" bp1,bp2,bm1,bm2 = ",4F10.8)')  bp1,bp2,bm1,bm2
                                    ! write(*,'(" wplw,wprw,wmlw,wmrw = ",4F10.8)')  wplw,wprw,wmlw,wmrw
                                    ! write(*,'(" wpl,wpr,wml,wmr = ",4F10.8)')  wpl,wpr,wml,wmr
                                    ! write(*,'(" agu = ",F10.8)') agu(vard,var)
                            !     end if
                            !v*uy
                            au = uu(i1,i2,i3,vc)
                            aum = u(i1,i2-1,i3,vc) - gvV
                            aup = u(i1,i2+1,i3,vc) - gvV
                            s1=0
                            s2=1
                            s3=0
                            var = uc
                            vard= vc
                            drl = dx(1)
                              Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*u(
     & i1   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
                              Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    u(
     & i1-s1,i2-s2,i3-s3,var)
                              Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) - 
     & 2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
                              Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) - 
     &    u(i1   ,i2   ,i3   ,var)
                              Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,
     & i2-s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
                              Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var)
                              Amr   = Apl
                              Bmr   = Bpl
                              betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 
     & 1./4.*Bpl**2
                              betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 
     & 1./4.*Bpr**2
                              betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 
     & 1./4.*Bml**2
                              betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 
     & 1./4.*Bmr**2
                              ep     = 1e-15
                              ap1    = (1./2.)/(ep + betapl)**(4)
                              ap2    = (1./2.)/(ep + betapr)**(4)
                              am1    = (1./2.)/(ep + betaml)**(4)
                              am2    = (1./2.)/(ep + betamr)**(4)
                              !calculate weights
                              wp1 = ap1/(ap1+ap2)
                              wp2 = ap2/(ap1+ap2)
                              wm1 = am1/(am1+am2)
                              wm2 = am2/(am1+am2)
                              !mapping of the weights
                              bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                              bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                              bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                              bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                              wplw = bp1/(bp1+bp2)
                              wprw = bp2/(bp1+bp2)
                              wmlw = bm1/(bm1+bm2)
                              wmrw = bm2/(bm1+bm2)
                              !get face value of u
                              up = au
                              um = au
                              !picking side for the weights
                              if (up.ge.0) then
                                 wpl = max(wplw,wprw)
                                 wpr = min(wplw,wprw)
                              else
                                 wpl = min(wplw,wprw)
                                 wpr = max(wplw,wprw)
                              end if
                              if (um.ge.0) then
                                 wml = max(wmlw,wmrw)
                                 wmr = min(wmlw,wmrw)
                              else
                                 wml = min(wmlw,wmrw)
                                 wmr = max(wmlw,wmrw)
                              end if
                              expAd = 0.
                              if (orderOneFix.eq.1) then
                                 ! if (augmentPlot.eq.1) then
                                 !    agu(uc,13)=0
                                 ! end if
                                 if (aum*aup .le. 0) then
                                    aumax = max(abs(au),abs(aup),abs(
     & aum))
                                    expAd = (1./12.)*aumax*( u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(i1,
     & i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*s2,
     & i3+2*s3,var))
                                    wpl = 0.5
                                    wpr = 0.5
                                    wml = 0.5
                                    wmr = 0.5
                                    ! plot au at where fix is turned on
                                    ! if (augmentPlot.eq.1) then
                                    !    agu(uc,13)=au
                                    ! end if
                                 end if
                              end if
                              ! wpl = 0.5
                              ! wpr = 0.5
                              ! wml = 0.5
                              ! wmr = 0.5
                              Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,var) 
     & + 5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
                              Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 5.*
     & u(i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
                              Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*
     & s3,var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
                              Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,var)
     &  + 5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
                              Fp  = wpl*Fpl + wpr*Fpr
                              Fm  = wml*Fml + wmr*Fmr
                                 agu(vard,var)= au*(Fp - Fm)/drl + 
     & expAd/drl
                            !     if ( i1 .eq. 14 .and. i2.eq.16 ) then
                                    !     if ( i1 .eq. 14) then
                                    ! write(*,'(" betapl,betapr,betaml,betamr = ",4F10.8)')  betapl,betapr,betaml,betamr
                                    ! write(*,'(" ap1,ap2,am1,am2 = ",4e10.8)')  ap1,ap2,am1,am2
                                    ! write(*,'(" wp1,wpr,wm1,wm2 = ",4F10.8)')  wp1,wp2,wm1,wm2
                                    ! write(*,'(" bp1,bp2,bm1,bm2 = ",4F10.8)')  bp1,bp2,bm1,bm2
                                    ! write(*,'(" wplw,wprw,wmlw,wmrw = ",4F10.8)')  wplw,wprw,wmlw,wmrw
                                    ! write(*,'(" wpl,wpr,wml,wmr = ",4F10.8)')  wpl,wpr,wml,wmr
                                    ! write(*,'(" agu = ",F10.8)') agu(vard,var)
                            !     end if
                            !v*vy
                            s1=0
                            s2=1
                            s3=0
                            var = vc
                            vard= vc
                            drl = dx(1)
                              Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*u(
     & i1   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
                              Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    u(
     & i1-s1,i2-s2,i3-s3,var)
                              Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) - 
     & 2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
                              Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) - 
     &    u(i1   ,i2   ,i3   ,var)
                              Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,
     & i2-s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
                              Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var)
                              Amr   = Apl
                              Bmr   = Bpl
                              betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 
     & 1./4.*Bpl**2
                              betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 
     & 1./4.*Bpr**2
                              betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 
     & 1./4.*Bml**2
                              betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 
     & 1./4.*Bmr**2
                              ep     = 1e-15
                              ap1    = (1./2.)/(ep + betapl)**(4)
                              ap2    = (1./2.)/(ep + betapr)**(4)
                              am1    = (1./2.)/(ep + betaml)**(4)
                              am2    = (1./2.)/(ep + betamr)**(4)
                              !calculate weights
                              wp1 = ap1/(ap1+ap2)
                              wp2 = ap2/(ap1+ap2)
                              wm1 = am1/(am1+am2)
                              wm2 = am2/(am1+am2)
                              !mapping of the weights
                              bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                              bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                              bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                              bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                              wplw = bp1/(bp1+bp2)
                              wprw = bp2/(bp1+bp2)
                              wmlw = bm1/(bm1+bm2)
                              wmrw = bm2/(bm1+bm2)
                              !get face value of u
                              up = au
                              um = au
                              !picking side for the weights
                              if (up.ge.0) then
                                 wpl = max(wplw,wprw)
                                 wpr = min(wplw,wprw)
                              else
                                 wpl = min(wplw,wprw)
                                 wpr = max(wplw,wprw)
                              end if
                              if (um.ge.0) then
                                 wml = max(wmlw,wmrw)
                                 wmr = min(wmlw,wmrw)
                              else
                                 wml = min(wmlw,wmrw)
                                 wmr = max(wmlw,wmrw)
                              end if
                              expAd = 0.
                              if (orderOneFix.eq.1) then
                                 ! if (augmentPlot.eq.1) then
                                 !    agu(uc,13)=0
                                 ! end if
                                 if (aum*aup .le. 0) then
                                    aumax = max(abs(au),abs(aup),abs(
     & aum))
                                    expAd = (1./12.)*aumax*( u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(i1,
     & i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*s2,
     & i3+2*s3,var))
                                    wpl = 0.5
                                    wpr = 0.5
                                    wml = 0.5
                                    wmr = 0.5
                                    ! plot au at where fix is turned on
                                    ! if (augmentPlot.eq.1) then
                                    !    agu(uc,13)=au
                                    ! end if
                                 end if
                              end if
                              ! wpl = 0.5
                              ! wpr = 0.5
                              ! wml = 0.5
                              ! wmr = 0.5
                              Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,var) 
     & + 5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
                              Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 5.*
     & u(i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
                              Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*
     & s3,var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
                              Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,var)
     &  + 5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
                              Fp  = wpl*Fpl + wpr*Fpr
                              Fm  = wml*Fml + wmr*Fmr
                                 agu(vard,var)= au*(Fp - Fm)/drl + 
     & expAd/drl
                            !     if ( i1 .eq. 14 .and. i2.eq.16 ) then
                                    !     if ( i1 .eq. 14) then
                                    ! write(*,'(" betapl,betapr,betaml,betamr = ",4F10.8)')  betapl,betapr,betaml,betamr
                                    ! write(*,'(" ap1,ap2,am1,am2 = ",4e10.8)')  ap1,ap2,am1,am2
                                    ! write(*,'(" wp1,wpr,wm1,wm2 = ",4F10.8)')  wp1,wp2,wm1,wm2
                                    ! write(*,'(" bp1,bp2,bm1,bm2 = ",4F10.8)')  bp1,bp2,bm1,bm2
                                    ! write(*,'(" wplw,wprw,wmlw,wmrw = ",4F10.8)')  wplw,wprw,wmlw,wmrw
                                    ! write(*,'(" wpl,wpr,wml,wmr = ",4F10.8)')  wpl,wpr,wml,wmr
                                    ! write(*,'(" agu = ",F10.8)') agu(vard,var)
                            !     end if
                           ! u*ux-------------------------------------------------------------------------
                           !calculate smooth factor
                           ! Apl   = u(i1+1,i2,i3,uc) - 2.*u(i1  ,i2,i3,uc)   + u(i1-1,i2,i3,uc)
                           ! Bpl   = u(i1+1,i2,i3,uc) -    u(i1-1,i2,i3,uc)
                           ! Apr   = u(i1+2,i2,i3,uc) - 2.*u(i1+1,i2,i3,uc) + u(i1,i2,i3,uc)
                           ! Bpr   = u(i1+2,i2,i3,uc) -    u(i1  ,i2,i3,uc)
                           ! Aml   = u(i1,i2,i3,uc) - 2.*u(i1-1,i2,i3,uc) + u(i1-2,i2,i3,uc)
                           ! Bml   = u(i1,i2,i3,uc) -    u(i1-2,i2,i3,uc)
                           ! Amr   = Apl
                           ! Bmr   = Bpl
                           ! betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 1./4.*Bpl**2
                           ! betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 1./4.*Bpr**2
                           ! betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 1./4.*Bml**2
                           ! betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 1./4.*Bmr**2
                           ! ep     = 1e-15
                           ! ap1    = (1./2.)/(ep + betapl)**(4)
                           ! ap2    = (1./2.)/(ep + betapr)**(4)
                           ! am1    = (1./2.)/(ep + betaml)**(4)
                           ! am2    = (1./2.)/(ep + betamr)**(4)
                           ! !calculate weights
                           ! wp1 = ap1/(ap1+ap2)
                           ! wp2 = ap2/(ap1+ap2)
                           ! wm1 = am1/(am1+am2)
                           ! wm2 = am2/(am1+am2)
                           ! !mapping of the weights
                           ! bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                           ! bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                           ! bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                           ! bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                           ! wplw = bp1/(bp1+bp2)
                           ! wprw = bp2/(bp1+bp2)
                           ! wmlw = bm1/(bm1+bm2)
                           ! wmrw = bm2/(bm1+bm2)
                           ! !get face value of u
                           ! up = au
                           ! um = au
                           ! !picking side for the weights
                           !  if (up.ge.0) then
                           !    wpl = max(wplw,wprw)
                           !    wpr = min(wplw,wprw)
                           ! else
                           !    wpl = min(wplw,wprw)
                           !    wpr = max(wplw,wprw)
                           ! end if
                           ! if (um.ge.0) then
                           !    wml = max(wmlw,wmrw)
                           !    wmr = min(wmlw,wmrw)
                           ! else
                           !    wml = min(wmlw,wmrw)
                           !    wmr = max(wmlw,wmrw)
                           ! end if
                           ! expAd = 0.
                           ! if (orderOneFix.eq.1) then 
                           !    if (augmentPlot.eq.1) then
                           !       agu(uc,13)=0
                           !    end if
                           !    if (aum*aup .le. 0) then
                           !       aumax = max(abs(au),abs(aup),abs(aum))
                           !       expAd = (1./12.)*aumax*( u(i1-2,i2,i3,uc) - 4.*u(i1-1,i2,i3,uc) + 6.*u(i1,i2,i3,uc) - 4.*u(i1+1,i2,i3,uc) +  u(i1+2,i2,i3,uc))
                           !       wpl = 0.5
                           !       wpr = 0.5
                           !       wml = 0.5
                           !       wmr = 0.5
                           !       ! plot au at where fix is turned on
                           !       if (augmentPlot.eq.1) then
                           !          agu(uc,13)=au
                           !       end if
                           !    end if
                           ! end if
                           ! wpl = 0.5
                           ! wpr = 0.5
                           ! wml = 0.5
                           ! wmr = 0.5
                           ! Fpl = 1./6.*(  -u(i1-1,i2,i3,uc) + 5.*u(i1,i2,i3,uc)   + 2.*u(i1+1,i2,i3,uc))
                           ! Fpr = 1./6.*( 2.*u(i1,i2,i3,uc)   + 5.*u(i1+1,i2,i3,uc) -   u(i1+2,i2,i3,uc))
                           ! Fml = 1./6.*(  -u(i1-2,i2,i3,uc) + 5.*u(i1-1,i2,i3,uc)   + 2.*u(i1,i2,i3,uc))
                           ! Fmr = 1./6.*( 2.*u(i1-1,i2,i3,uc) + 5.*u(i1,i2,i3,uc) -       u(i1+1,i2,i3,uc))
                           ! Fp  = wpl*Fpl + wpr*Fpr
                           ! Fm  = wml*Fml + wmr*Fmr
                           ! agu(uc,uc)= au*(Fp - Fm)/(dx(0)) + expAd*dx(0)**(-1)
                           ! if (augmentPlot.eq.1) then
                           !    agu(uc,3)=au
                           !    agu(uc,4)=wpl
                           !    agu(uc,5)=wpr
                           !    agu(uc,6)=wml
                           !    agu(uc,7)=wmr
                           ! end if
                           ! u*vx-------------------------------------------------------------------------
                           !calculate smooth factor
                           ! Apl   = u(i1+1,i2,i3,vc) - 2.*u(i1,i2,i3,vc)   + u(i1-1,i2,i3,vc)
                           ! Bpl   = u(i1+1,i2,i3,vc) -    u(i1-1,i2,i3,vc)
                           ! Apr   = u(i1+2,i2,i3,vc) - 2.*u(i1+1,i2,i3,vc) + u(i1,i2,i3,vc)
                           ! Bpr   = u(i1+2,i2,i3,vc) -    u(i1,i2,i3,vc)
                           ! Aml   = u(i1,i2,i3,vc) - 2.*u(i1-1,i2,i3,vc) + u(i1-2,i2,i3,vc)
                           ! Bml   = u(i1,i2,i3,vc) -    u(i1-2,i2,i3,vc)
                           ! Amr   = Apl
                           ! Bmr   = Bpl
                           ! betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 1./4.*Bpl**2
                           ! betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 1./4.*Bpr**2
                           ! betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 1./4.*Bml**2
                           ! betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 1./4.*Bmr**2
                           ! ep     = 1e-15
                           ! ap1    = (1./2.)/(ep + betapl)**(4)
                           ! ap2    = (1./2.)/(ep + betapr)**(4)
                           ! am1    = (1./2.)/(ep + betaml)**(4)
                           ! am2    = (1./2.)/(ep + betamr)**(4)
                           ! !calculate weights
                           ! wp1 = ap1/(ap1+ap2)
                           ! wp2 = ap2/(ap1+ap2)
                           ! wm1 = am1/(am1+am2)
                           ! wm2 = am2/(am1+am2)
                           ! !mapping of the weights
                           ! bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                           ! bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                           ! bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                           ! bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                           ! wplw = bp1/(bp1+bp2)
                           ! wprw = bp2/(bp1+bp2)
                           ! wmlw = bm1/(bm1+bm2)
                           ! wmrw = bm2/(bm1+bm2)
                           ! !get face value of u
                           ! ! up = 0.5*(u(i1  ,i2,i3,uc)+u(i1+1,i2,i3,uc))
                           ! ! um = 0.5*(u(i1-1,i2,i3,uc)+u(i1  ,i2,i3,uc))
                           ! up = au
                           ! um = au
                           ! !picking side for the weights
                           ! if (up.ge.0) then
                           !    wpl = max(wplw,wprw)
                           !    wpr = min(wplw,wprw)
                           ! else
                           !    wpl = min(wplw,wprw)
                           !    wpr = max(wplw,wprw)
                           ! end if
                           ! if (um.ge.0) then
                           !    wml = max(wmlw,wmrw)
                           !    wmr = min(wmlw,wmrw)
                           ! else
                           !    wml = min(wmlw,wmrw)
                           !    wmr = max(wmlw,wmrw)
                           ! end if
                           ! expAd = 0.
                           ! if (orderOneFix.eq.1) then
                           !    if (aum*aup .le. 0) then
                           !       aumax = max(abs(au),abs(aup),abs(aum))
                           !       expAd = (1./12.)*aumax*( u(i1-2,i2,i3,vc) - 4.*u(i1-1,i2,i3,vc) + 6.*u(i1,i2,i3,vc) - 4.*u(i1+1,i2,i3,vc) +  u(i1+2,i2,i3,vc))
                           !       wpl = 0.5
                           !       wpr = 0.5
                           !       wml = 0.5
                           !       wmr = 0.5
                           !    end if
                           ! end if
                           ! Fpl = 1./6.*(  -u(i1-1,i2,i3,vc) + 5.*u(i1,i2,i3,vc)   + 2.*u(i1+1,i2,i3,vc))
                           ! Fpr = 1./6.*( 2.*u(i1,i2,i3,vc)   + 5.*u(i1+1,i2,i3,vc) -   u(i1+2,i2,i3,vc))
                           ! Fml = 1./6.*(  -u(i1-2,i2,i3,vc) + 5.*u(i1-1,i2,i3,vc)   + 2.*u(i1,i2,i3,vc))
                           ! Fmr = 1./6.*( 2.*u(i1-1,i2,i3,vc) + 5.*u(i1,i2,i3,vc) -       u(i1+1,i2,i3,vc))
                           ! Fp  = wpl*Fpl + wpr*Fpr
                           ! Fm  = wml*Fml + wmr*Fmr
                           ! agu(uc,vc)= au*(Fp - Fm)/(dx(0)) + expAd*dx(0)**(-1) 
                           ! if (augmentPlot .eq. 1) then 
                           !    agu(vc,4)=wpl
                           !    agu(vc,5)=wpr
                           !    agu(vc,6)=wml
                           !    agu(vc,7)=wmr
                           ! end if
                           !v*uy-------------------------------------------------------------------------
                           ! au = u(i1,i2,i3,vc)
                           ! aum = u(i1,i2-1,i3,vc)
                           ! aup = u(i1,i2+1,i3,vc)
                           !calculate smooth factor
                           ! Apl   = u(i1,i2+1,i3,uc) - 2.*u(i1,i2,i3,uc)   + u(i1,i2-1,i3,uc)
                           ! Bpl   = u(i1,i2+1,i3,uc) -    u(i1,i2-1,i3,uc)
                           ! Apr   = u(i1,i2+2,i3,uc) - 2.*u(i1,i2+1,i3,uc) + u(i1,i2,i3,uc)
                           ! Bpr   = u(i1,i2+2,i3,uc) -    u(i1,i2,i3,uc)
                           ! Aml   = u(i1,i2,i3,uc) - 2.*u(i1,i2-1,i3,uc) + u(i1,i2-2,i3,uc)
                           ! Bml   = u(i1,i2,i3,uc) -    u(i1,i2-2,i3,uc)
                           ! Amr   = Apl
                           ! Bmr   = Bpl
                           ! betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 1./4.*Bpl**2
                           ! betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 1./4.*Bpr**2
                           ! betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 1./4.*Bml**2
                           ! betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 1./4.*Bmr**2
                           ! ep     = 1e-15
                           ! ap1    = (1./2.)/(ep + betapl)**(4)
                           ! ap2    = (1./2.)/(ep + betapr)**(4)
                           ! am1    = (1./2.)/(ep + betaml)**(4)
                           ! am2    = (1./2.)/(ep + betamr)**(4)
                           ! !calculate weights
                           ! wp1 = ap1/(ap1+ap2)
                           ! wp2 = ap2/(ap1+ap2)
                           ! wm1 = am1/(am1+am2)
                           ! wm2 = am2/(am1+am2)
                           ! !mapping of the weights
                           ! bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                           ! bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                           ! bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                           ! bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                           ! wplw = bp1/(bp1+bp2)
                           ! wprw = bp2/(bp1+bp2)
                           ! wmlw = bm1/(bm1+bm2)
                           ! wmrw = bm2/(bm1+bm2)
                           ! !get face value of v
                           ! ! up = 0.5*(u(i1  ,i2+1,i3,vc)+u(i1,i2 ,i3,vc))
                           ! ! um = 0.5*(u(i1  ,i2  ,i3,vc)+u(i1,i2-1,i3,vc))
                           ! up = au
                           ! um = au
                           ! !picking side for the weights
                           !  if (up.ge.0) then
                           !    wpl = max(wplw,wprw)
                           !    wpr = min(wplw,wprw)
                           ! else
                           !    wpl = min(wplw,wprw)
                           !    wpr = max(wplw,wprw)
                           ! end if
                           ! if (um.ge.0) then
                           !    wml = max(wmlw,wmrw)
                           !    wmr = min(wmlw,wmrw)
                           ! else
                           !    wml = min(wmlw,wmrw)
                           !    wmr = max(wmlw,wmrw)
                           ! end if
                           ! expAd = 0.
                           ! if (orderOneFix.eq.1) then 
                           !    if (augmentPlot .eq. 1) then 
                           !       agu(vc,13)=0
                           !    end if
                           !    if (aum*aup .le. 0) then
                           !       aumax = max(abs(au),abs(aup),abs(aum))
                           !       expAd = (1./12.)*aumax*( u(i1,i2-2,i3,uc) - 4.*u(i1,i2-1,i3,uc) + 6.*u(i1,i2,i3,uc) - 4.*u(i1,i2+1,i3,uc) +  u(i1,i2+2,i3,uc))
                           !       wpl = 0.5
                           !       wpr = 0.5
                           !       wml = 0.5
                           !       wmr = 0.5
                           !       if (augmentPlot .eq. 1) then 
                           !          agu(vc,13)=au
                           !       end if
                           !    end if
                           ! end if
                           ! Fpl = 1./6.*(  -u(i1,i2-1,i3,uc) + 5.*u(i1,i2,i3,uc)   + 2.*u(i1,i2+1,i3,uc))
                           ! Fpr = 1./6.*( 2.*u(i1,i2,i3,uc)   + 5.*u(i1,i2+1,i3,uc) -   u(i1,i2+2,i3,uc))
                           ! Fml = 1./6.*(  -u(i1,i2-2,i3,uc) + 5.*u(i1,i2-1,i3,uc)   + 2.*u(i1,i2,i3,uc))
                           ! Fmr = 1./6.*( 2.*u(i1,i2-1,i3,uc) + 5.*u(i1,i2,i3,uc) -       u(i1,i2+1,i3,uc))
                           ! Fp  = wpl*Fpl + wpr*Fpr
                           ! Fm  = wml*Fml + wmr*Fmr
                           ! agu(vc,uc)= au*(Fp - Fm)/(dx(1)) + expAd*dx(1)**(-1)
                           ! if (augmentPlot .eq. 1) then 
                           !    agu(uc,8)=au
                           !    agu(uc,9) =wpl
                           !    agu(uc,10)=wpr
                           !    agu(uc,11)=wml
                           !    agu(uc,12)=wmr
                           ! end if
                           !v*vy---------------------------------------------------------------
                           ! Apl   = u(i1,i2+1,i3,vc) - 2.*u(i1,i2,i3,vc)   + u(i1,i2-1,i3,vc)
                           ! Bpl   = u(i1,i2+1,i3,vc) -    u(i1,i2-1,i3,vc)
                           ! Apr   = u(i1,i2+2,i3,vc) - 2.*u(i1,i2+1,i3,vc) + u(i1,i2,i3,vc)
                           ! Bpr   = u(i1,i2+2,i3,vc) -    u(i1,i2,i3,vc)
                           ! Aml   = u(i1,i2,i3,vc) - 2.*u(i1,i2-1,i3,vc) + u(i1,i2-2,i3,vc)
                           ! Bml   = u(i1,i2,i3,vc) -    u(i1,i2-2,i3,vc)
                           ! Amr   = Apl
                           ! Bmr   = Bpl
                           ! betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 1./4.*Bpl**2
                           ! betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 1./4.*Bpr**2
                           ! betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 1./4.*Bml**2
                           ! betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 1./4.*Bmr**2
                           ! ep     = 1e-15
                           ! ap1    = (1./2.)/(ep + betapl)**(4)
                           ! ap2    = (1./2.)/(ep + betapr)**(4)
                           ! am1    = (1./2.)/(ep + betaml)**(4)
                           ! am2    = (1./2.)/(ep + betamr)**(4)
                           ! !calculate weights
                           ! wp1 = ap1/(ap1+ap2)
                           ! wp2 = ap2/(ap1+ap2)
                           ! wm1 = am1/(am1+am2)
                           ! wm2 = am2/(am1+am2)
                           ! !mapping of the weights
                           ! bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                           ! bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                           ! bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                           ! bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                           ! wplw = bp1/(bp1+bp2)
                           ! wprw = bp2/(bp1+bp2)
                           ! wmlw = bm1/(bm1+bm2)
                           ! wmrw = bm2/(bm1+bm2)
                           ! !get face value of v
                           !  ! up = 0.5*(u(i1  ,i2,i3,vc)  + u(i1,i2+1,i3,vc))
                           !  ! um = 0.5*(u(i1-1,i2,i3,vc)  + u(i1,i2  ,i3,vc))
                           ! up = au
                           ! um = au
                           ! !picking side for the weights
                           !  if (up.ge.0) then
                           !    wpl = max(wplw,wprw)
                           !    wpr = min(wplw,wprw)
                           ! else
                           !    wpl = min(wplw,wprw)
                           !    wpr = max(wplw,wprw)
                           ! end if
                           ! if (um.ge.0) then
                           !    wml = max(wmlw,wmrw)
                           !    wmr = min(wmlw,wmrw)
                           ! else
                           !    wml = min(wmlw,wmrw)
                           !    wmr = max(wmlw,wmrw)
                           ! end if
                           ! expAd=0.
                           ! if (orderOneFix.eq.1) then
                           !    if (aum*aup .le. 0) then
                           !       aumax = max(abs(au),abs(aup),abs(aum))
                           !       expAd = (1./12.)*aumax*( u(i1,i2-2,i3,vc) - 4.*u(i1,i2-1,i3,vc) + 6.*u(i1,i2,i3,vc) - 4.*u(i1,i2+1,i3,vc) +  u(i1,i2+2,i3,vc))
                           !       wpl = 0.5
                           !       wpr = 0.5
                           !       wml = 0.5
                           !       wmr = 0.5
                           !    end if
                           ! end if
                           ! Fpl = 1./6.*(  -u(i1,i2-1,i3,vc) + 5.*u(i1,i2,i3,vc)   + 2.*u(i1,i2+1,i3,vc))
                           ! Fpr = 1./6.*( 2.*u(i1,i2,i3,vc)   + 5.*u(i1,i2+1,i3,vc) -   u(i1,i2+2,i3,vc))
                           ! Fml = 1./6.*(  -u(i1,i2-2,i3,vc) + 5.*u(i1,i2-1,i3,vc)   + 2.*u(i1,i2,i3,vc))
                           ! Fmr = 1./6.*( 2.*u(i1,i2-1,i3,vc) + 5.*u(i1,i2,i3,vc) -       u(i1,i2+1,i3,vc))
                           ! Fp  = wpl*Fpl + wpr*Fpr
                           ! Fm  = wml*Fml + wmr*Fmr
                           ! agu(vc,vc)= au*(Fp - Fm)/(dx(1)) + expAd*dx(1)**(-1)
                           ! if (augmentPlot .eq. 1) then  
                           !    agu(vc,9) =wpl
                           !    agu(vc,10)=wpr
                           !    agu(vc,11)=wml
                           !    agu(vc,12)=wmr
                           ! end if
                        else ! if (upwindOrder.eq.4)
                           write(*,'(" getBewnoAdvection: only 4th 
     & order is avaliable now")' )
                           stop 777
                        end if
                    ! #If ("BWENO" == "CENTERED") 
                    !  ! -- centered advection ---
                    !  ! write(*,'(" getAdvection -- centered")')
                    !  #If "2" == "2"
                    !   ! -- 2D --
                    !   agu(uc,uc)=UU(uc)*UX(uc)
                    !   agu(vc,uc)=UU(vc)*UY(uc)
                    !   agu(uc,vc)=UU(uc)*UX(vc)
                    !   agu(vc,vc)=UU(vc)*UY(vc)
                    !  #Else
                    !   ! -- 3D -- *check me*
                    !   agu(uc,uc)=UU(uc)*UX(uc)
                    !   agu(vc,uc)=UU(vc)*UY(uc)
                    !   agu(wc,uc)=UU(wc)*UY(uc)
                    !   agu(uc,vc)=UU(uc)*UX(vc)
                    !   agu(vc,vc)=UU(vc)*UY(vc)
                    !   agu(wc,vc)=UU(wc)*UY(vc)
                    !   agu(uc,wc)=UU(uc)*UX(wc)
                    !   agu(vc,wc)=UU(vc)*UY(wc)
                    !   agu(wc,wc)=UU(wc)*UY(wc)
                    !  #End
                    ! #Elif "BWENO" == "BWENO" 
                    !   ! --- upwind scheme ---
                    !   ! for testing output this next message:
                    !   if( t.le. 0. )then
                    !     write(*,'(" getAdvection upwind scheme (7)")') 
                    !   end if
                    !   getUpwindAdvection(u,i1,i2,i3,NONE,2,4,rectangular, agu)
                    ! #Elif "BWENO" == "BWENO" 
                    !   ! --- Bweno scheme ---
                    !   ! for testing output this message:
                    !   if( t.le. 0. )then
                    !      write(*,'(" getAdvection BWENO scheme (7)")') 
                    !   end if
                    !   getBwenoAdvection(u,i1,i2,i3,NONE,2,4,rectangular, agu)
                    ! #Else
                    !   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
                    !   stop 999
                    ! #End
                       ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-ux42r(
     & i1,i2,i3,pc)+nu*ulaplacian42r(i1,i2,i3,uc)
                       ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-uy42r(
     & i1,i2,i3,pc)+nu*ulaplacian42r(i1,i2,i3,vc)
                 !       u(i1,i2,i3,3)= agu(uc,3) ! au in r direction 
                 !       u(i1,i2,i3,4)= agu(uc,8) ! au in s direction
                 !       u(i1,i2,i3,5)= agu(uc,4)  ! wpl uur
                 !       u(i1,i2,i3,6)= agu(uc,5)  ! wpr uur
                 !       u(i1,i2,i3,7)= agu(uc,6)  ! wml uur
                 !       u(i1,i2,i3,8)= agu(uc,7)  ! wmr uur
                 !       u(i1,i2,i3,9)= agu(uc,9)  ! wpl uus
                 !       u(i1,i2,i3,10)= agu(uc,10)  ! wpr uus
                 !       u(i1,i2,i3,11)= agu(uc,11)  ! wml uus
                 !       u(i1,i2,i3,12)= agu(uc,12)  ! wmr uus
                 ! !
                 !       u(i1,i2,i3,13)= agu(vc,4)  ! wpl uvr
                 !       u(i1,i2,i3,14)= agu(vc,5)  ! wpr uvr
                 !       u(i1,i2,i3,15)= agu(vc,6)  ! wml uvr
                 !       u(i1,i2,i3,16)= agu(vc,7)  ! wmr uvr
                 !       u(i1,i2,i3,17)= agu(vc,9)  ! wpl uvs
                 !       u(i1,i2,i3,18)= agu(vc,10)  ! wpr uvs
                 !       u(i1,i2,i3,19)= agu(vc,11)  ! wml uvs
                 !       u(i1,i2,i3,20)= agu(vc,12)  ! wmr uvs
                 !       u(i1,i2,i3,21)= agu(uc,13)  ! au at where fix is turned on 
                 !       u(i1,i2,i3,22)= agu(vc,13)  ! au at where fix is turned on 
                end if
               end do
               end do
               end do
              else
               do i3=n3a,n3b
               do i2=n2a,n2b
               do i1=n1a,n1b
                 ! INS, no AD
                   ! --- upwind approximations ---
                     ! --- Bweno scheme ---
                     ! for testing output this message:
                     ! if( t.le. 0. )then
                     !    write(*,'(" getAdvection BWENO scheme (7)")') 
                     ! end if
                        gvU = u(i1,i2,i3,uc) - uu(i1,i2,i3,uc)
                        gvV = u(i1,i2,i3,vc) - uu(i1,i2,i3,vc)
                        gvW = u(i1,i2,i3,wc) - uu(i1,i2,i3,wc)
                        orderOneFix = 1
                        augmentPlot = 1
                        if( upwindOrder.eq.4 )then
                       ! -- fourth order bweno --
                           au = uu(i1,i2,i3,uc)
                           aum = u(i1-1,i2,i3,uc) - gvU
                           aup = u(i1+1,i2,i3,uc) - gvU
                           !usage
                           !s1   is the shifted grid in r1 direction
                           !var  is the variable taking derivative (uc,vc,wc)
                           !vard is the variable of direction
                           !agu is only passed in for augmented plot
                           !expAd 
                           ! u*ux
                           s1=1
                           s2=0
                           s3=0
                           var = uc
                           vard= uc
                           drl = dx(0)
                             Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*u(
     & i1   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
                             Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    u(
     & i1-s1,i2-s2,i3-s3,var)
                             Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) - 
     & 2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
                             Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) -  
     &   u(i1   ,i2   ,i3   ,var)
                             Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,i2-
     & s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
                             Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,i2-
     & 2*s2,i3-2*s3,var)
                             Amr   = Apl
                             Bmr   = Bpl
                             betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 
     & 1./4.*Bpl**2
                             betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 
     & 1./4.*Bpr**2
                             betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 
     & 1./4.*Bml**2
                             betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 
     & 1./4.*Bmr**2
                             ep     = 1e-15
                             ap1    = (1./2.)/(ep + betapl)**(4)
                             ap2    = (1./2.)/(ep + betapr)**(4)
                             am1    = (1./2.)/(ep + betaml)**(4)
                             am2    = (1./2.)/(ep + betamr)**(4)
                             !calculate weights
                             wp1 = ap1/(ap1+ap2)
                             wp2 = ap2/(ap1+ap2)
                             wm1 = am1/(am1+am2)
                             wm2 = am2/(am1+am2)
                             !mapping of the weights
                             bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                             bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                             bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                             bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                             wplw = bp1/(bp1+bp2)
                             wprw = bp2/(bp1+bp2)
                             wmlw = bm1/(bm1+bm2)
                             wmrw = bm2/(bm1+bm2)
                             !get face value of u
                             up = au
                             um = au
                             !picking side for the weights
                             if (up.ge.0) then
                                wpl = max(wplw,wprw)
                                wpr = min(wplw,wprw)
                             else
                                wpl = min(wplw,wprw)
                                wpr = max(wplw,wprw)
                             end if
                             if (um.ge.0) then
                                wml = max(wmlw,wmrw)
                                wmr = min(wmlw,wmrw)
                             else
                                wml = min(wmlw,wmrw)
                                wmr = max(wmlw,wmrw)
                             end if
                             expAd = 0.
                             if (orderOneFix.eq.1) then
                                ! if (augmentPlot.eq.1) then
                                !    agu(uc,13)=0
                                ! end if
                                if (aum*aup .le. 0) then
                                   aumax = max(abs(au),abs(aup),abs(
     & aum))
                                   expAd = (1./12.)*aumax*( u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(i1,
     & i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*s2,
     & i3+2*s3,var))
                                   wpl = 0.5
                                   wpr = 0.5
                                   wml = 0.5
                                   wmr = 0.5
                                   ! plot au at where fix is turned on
                                   ! if (augmentPlot.eq.1) then
                                   !    agu(uc,13)=au
                                   ! end if
                                end if
                             end if
                             ! wpl = 0.5
                             ! wpr = 0.5
                             ! wml = 0.5
                             ! wmr = 0.5
                             Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,var) +
     &  5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
                             Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 5.*u(
     & i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
                             Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*s3,
     & var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
                             Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,var) 
     & + 5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
                             Fp  = wpl*Fpl + wpr*Fpr
                             Fm  = wml*Fml + wmr*Fmr
                                agu(vard,var)= au*(Fp - Fm)/drl + 
     & expAd/drl
                           !     if ( i1 .eq. 14 .and. i2.eq.16 ) then
                                   !     if ( i1 .eq. 14) then
                                   ! write(*,'(" betapl,betapr,betaml,betamr = ",4F10.8)')  betapl,betapr,betaml,betamr
                                   ! write(*,'(" ap1,ap2,am1,am2 = ",4e10.8)')  ap1,ap2,am1,am2
                                   ! write(*,'(" wp1,wpr,wm1,wm2 = ",4F10.8)')  wp1,wp2,wm1,wm2
                                   ! write(*,'(" bp1,bp2,bm1,bm2 = ",4F10.8)')  bp1,bp2,bm1,bm2
                                   ! write(*,'(" wplw,wprw,wmlw,wmrw = ",4F10.8)')  wplw,wprw,wmlw,wmrw
                                   ! write(*,'(" wpl,wpr,wml,wmr = ",4F10.8)')  wpl,wpr,wml,wmr
                                   ! write(*,'(" agu = ",F10.8)') agu(vard,var)
                           !     end if
                           !u*vx
                           s1=1
                           s2=0
                           s3=0
                           var = vc
                           vard= uc
                           drl = dx(0)
                             Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*u(
     & i1   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
                             Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    u(
     & i1-s1,i2-s2,i3-s3,var)
                             Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) - 
     & 2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
                             Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) -  
     &   u(i1   ,i2   ,i3   ,var)
                             Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,i2-
     & s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
                             Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,i2-
     & 2*s2,i3-2*s3,var)
                             Amr   = Apl
                             Bmr   = Bpl
                             betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 
     & 1./4.*Bpl**2
                             betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 
     & 1./4.*Bpr**2
                             betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 
     & 1./4.*Bml**2
                             betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 
     & 1./4.*Bmr**2
                             ep     = 1e-15
                             ap1    = (1./2.)/(ep + betapl)**(4)
                             ap2    = (1./2.)/(ep + betapr)**(4)
                             am1    = (1./2.)/(ep + betaml)**(4)
                             am2    = (1./2.)/(ep + betamr)**(4)
                             !calculate weights
                             wp1 = ap1/(ap1+ap2)
                             wp2 = ap2/(ap1+ap2)
                             wm1 = am1/(am1+am2)
                             wm2 = am2/(am1+am2)
                             !mapping of the weights
                             bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                             bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                             bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                             bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                             wplw = bp1/(bp1+bp2)
                             wprw = bp2/(bp1+bp2)
                             wmlw = bm1/(bm1+bm2)
                             wmrw = bm2/(bm1+bm2)
                             !get face value of u
                             up = au
                             um = au
                             !picking side for the weights
                             if (up.ge.0) then
                                wpl = max(wplw,wprw)
                                wpr = min(wplw,wprw)
                             else
                                wpl = min(wplw,wprw)
                                wpr = max(wplw,wprw)
                             end if
                             if (um.ge.0) then
                                wml = max(wmlw,wmrw)
                                wmr = min(wmlw,wmrw)
                             else
                                wml = min(wmlw,wmrw)
                                wmr = max(wmlw,wmrw)
                             end if
                             expAd = 0.
                             if (orderOneFix.eq.1) then
                                ! if (augmentPlot.eq.1) then
                                !    agu(uc,13)=0
                                ! end if
                                if (aum*aup .le. 0) then
                                   aumax = max(abs(au),abs(aup),abs(
     & aum))
                                   expAd = (1./12.)*aumax*( u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(i1,
     & i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*s2,
     & i3+2*s3,var))
                                   wpl = 0.5
                                   wpr = 0.5
                                   wml = 0.5
                                   wmr = 0.5
                                   ! plot au at where fix is turned on
                                   ! if (augmentPlot.eq.1) then
                                   !    agu(uc,13)=au
                                   ! end if
                                end if
                             end if
                             ! wpl = 0.5
                             ! wpr = 0.5
                             ! wml = 0.5
                             ! wmr = 0.5
                             Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,var) +
     &  5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
                             Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 5.*u(
     & i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
                             Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*s3,
     & var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
                             Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,var) 
     & + 5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
                             Fp  = wpl*Fpl + wpr*Fpr
                             Fm  = wml*Fml + wmr*Fmr
                                agu(vard,var)= au*(Fp - Fm)/drl + 
     & expAd/drl
                           !     if ( i1 .eq. 14 .and. i2.eq.16 ) then
                                   !     if ( i1 .eq. 14) then
                                   ! write(*,'(" betapl,betapr,betaml,betamr = ",4F10.8)')  betapl,betapr,betaml,betamr
                                   ! write(*,'(" ap1,ap2,am1,am2 = ",4e10.8)')  ap1,ap2,am1,am2
                                   ! write(*,'(" wp1,wpr,wm1,wm2 = ",4F10.8)')  wp1,wp2,wm1,wm2
                                   ! write(*,'(" bp1,bp2,bm1,bm2 = ",4F10.8)')  bp1,bp2,bm1,bm2
                                   ! write(*,'(" wplw,wprw,wmlw,wmrw = ",4F10.8)')  wplw,wprw,wmlw,wmrw
                                   ! write(*,'(" wpl,wpr,wml,wmr = ",4F10.8)')  wpl,wpr,wml,wmr
                                   ! write(*,'(" agu = ",F10.8)') agu(vard,var)
                           !     end if
                           !v*uy
                           au = uu(i1,i2,i3,vc)
                           aum = u(i1,i2-1,i3,vc) - gvV
                           aup = u(i1,i2+1,i3,vc) - gvV
                           s1=0
                           s2=1
                           s3=0
                           var = uc
                           vard= vc
                           drl = dx(1)
                             Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*u(
     & i1   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
                             Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    u(
     & i1-s1,i2-s2,i3-s3,var)
                             Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) - 
     & 2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
                             Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) -  
     &   u(i1   ,i2   ,i3   ,var)
                             Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,i2-
     & s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
                             Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,i2-
     & 2*s2,i3-2*s3,var)
                             Amr   = Apl
                             Bmr   = Bpl
                             betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 
     & 1./4.*Bpl**2
                             betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 
     & 1./4.*Bpr**2
                             betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 
     & 1./4.*Bml**2
                             betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 
     & 1./4.*Bmr**2
                             ep     = 1e-15
                             ap1    = (1./2.)/(ep + betapl)**(4)
                             ap2    = (1./2.)/(ep + betapr)**(4)
                             am1    = (1./2.)/(ep + betaml)**(4)
                             am2    = (1./2.)/(ep + betamr)**(4)
                             !calculate weights
                             wp1 = ap1/(ap1+ap2)
                             wp2 = ap2/(ap1+ap2)
                             wm1 = am1/(am1+am2)
                             wm2 = am2/(am1+am2)
                             !mapping of the weights
                             bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                             bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                             bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                             bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                             wplw = bp1/(bp1+bp2)
                             wprw = bp2/(bp1+bp2)
                             wmlw = bm1/(bm1+bm2)
                             wmrw = bm2/(bm1+bm2)
                             !get face value of u
                             up = au
                             um = au
                             !picking side for the weights
                             if (up.ge.0) then
                                wpl = max(wplw,wprw)
                                wpr = min(wplw,wprw)
                             else
                                wpl = min(wplw,wprw)
                                wpr = max(wplw,wprw)
                             end if
                             if (um.ge.0) then
                                wml = max(wmlw,wmrw)
                                wmr = min(wmlw,wmrw)
                             else
                                wml = min(wmlw,wmrw)
                                wmr = max(wmlw,wmrw)
                             end if
                             expAd = 0.
                             if (orderOneFix.eq.1) then
                                ! if (augmentPlot.eq.1) then
                                !    agu(uc,13)=0
                                ! end if
                                if (aum*aup .le. 0) then
                                   aumax = max(abs(au),abs(aup),abs(
     & aum))
                                   expAd = (1./12.)*aumax*( u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(i1,
     & i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*s2,
     & i3+2*s3,var))
                                   wpl = 0.5
                                   wpr = 0.5
                                   wml = 0.5
                                   wmr = 0.5
                                   ! plot au at where fix is turned on
                                   ! if (augmentPlot.eq.1) then
                                   !    agu(uc,13)=au
                                   ! end if
                                end if
                             end if
                             ! wpl = 0.5
                             ! wpr = 0.5
                             ! wml = 0.5
                             ! wmr = 0.5
                             Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,var) +
     &  5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
                             Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 5.*u(
     & i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
                             Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*s3,
     & var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
                             Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,var) 
     & + 5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
                             Fp  = wpl*Fpl + wpr*Fpr
                             Fm  = wml*Fml + wmr*Fmr
                                agu(vard,var)= au*(Fp - Fm)/drl + 
     & expAd/drl
                           !     if ( i1 .eq. 14 .and. i2.eq.16 ) then
                                   !     if ( i1 .eq. 14) then
                                   ! write(*,'(" betapl,betapr,betaml,betamr = ",4F10.8)')  betapl,betapr,betaml,betamr
                                   ! write(*,'(" ap1,ap2,am1,am2 = ",4e10.8)')  ap1,ap2,am1,am2
                                   ! write(*,'(" wp1,wpr,wm1,wm2 = ",4F10.8)')  wp1,wp2,wm1,wm2
                                   ! write(*,'(" bp1,bp2,bm1,bm2 = ",4F10.8)')  bp1,bp2,bm1,bm2
                                   ! write(*,'(" wplw,wprw,wmlw,wmrw = ",4F10.8)')  wplw,wprw,wmlw,wmrw
                                   ! write(*,'(" wpl,wpr,wml,wmr = ",4F10.8)')  wpl,wpr,wml,wmr
                                   ! write(*,'(" agu = ",F10.8)') agu(vard,var)
                           !     end if
                           !v*vy
                           s1=0
                           s2=1
                           s3=0
                           var = vc
                           vard= vc
                           drl = dx(1)
                             Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*u(
     & i1   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
                             Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    u(
     & i1-s1,i2-s2,i3-s3,var)
                             Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) - 
     & 2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
                             Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) -  
     &   u(i1   ,i2   ,i3   ,var)
                             Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,i2-
     & s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
                             Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,i2-
     & 2*s2,i3-2*s3,var)
                             Amr   = Apl
                             Bmr   = Bpl
                             betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 
     & 1./4.*Bpl**2
                             betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 
     & 1./4.*Bpr**2
                             betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 
     & 1./4.*Bml**2
                             betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 
     & 1./4.*Bmr**2
                             ep     = 1e-15
                             ap1    = (1./2.)/(ep + betapl)**(4)
                             ap2    = (1./2.)/(ep + betapr)**(4)
                             am1    = (1./2.)/(ep + betaml)**(4)
                             am2    = (1./2.)/(ep + betamr)**(4)
                             !calculate weights
                             wp1 = ap1/(ap1+ap2)
                             wp2 = ap2/(ap1+ap2)
                             wm1 = am1/(am1+am2)
                             wm2 = am2/(am1+am2)
                             !mapping of the weights
                             bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                             bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                             bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                             bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                             wplw = bp1/(bp1+bp2)
                             wprw = bp2/(bp1+bp2)
                             wmlw = bm1/(bm1+bm2)
                             wmrw = bm2/(bm1+bm2)
                             !get face value of u
                             up = au
                             um = au
                             !picking side for the weights
                             if (up.ge.0) then
                                wpl = max(wplw,wprw)
                                wpr = min(wplw,wprw)
                             else
                                wpl = min(wplw,wprw)
                                wpr = max(wplw,wprw)
                             end if
                             if (um.ge.0) then
                                wml = max(wmlw,wmrw)
                                wmr = min(wmlw,wmrw)
                             else
                                wml = min(wmlw,wmrw)
                                wmr = max(wmlw,wmrw)
                             end if
                             expAd = 0.
                             if (orderOneFix.eq.1) then
                                ! if (augmentPlot.eq.1) then
                                !    agu(uc,13)=0
                                ! end if
                                if (aum*aup .le. 0) then
                                   aumax = max(abs(au),abs(aup),abs(
     & aum))
                                   expAd = (1./12.)*aumax*( u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(i1,
     & i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*s2,
     & i3+2*s3,var))
                                   wpl = 0.5
                                   wpr = 0.5
                                   wml = 0.5
                                   wmr = 0.5
                                   ! plot au at where fix is turned on
                                   ! if (augmentPlot.eq.1) then
                                   !    agu(uc,13)=au
                                   ! end if
                                end if
                             end if
                             ! wpl = 0.5
                             ! wpr = 0.5
                             ! wml = 0.5
                             ! wmr = 0.5
                             Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,var) +
     &  5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
                             Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 5.*u(
     & i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
                             Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*s3,
     & var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
                             Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,var) 
     & + 5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
                             Fp  = wpl*Fpl + wpr*Fpr
                             Fm  = wml*Fml + wmr*Fmr
                                agu(vard,var)= au*(Fp - Fm)/drl + 
     & expAd/drl
                           !     if ( i1 .eq. 14 .and. i2.eq.16 ) then
                                   !     if ( i1 .eq. 14) then
                                   ! write(*,'(" betapl,betapr,betaml,betamr = ",4F10.8)')  betapl,betapr,betaml,betamr
                                   ! write(*,'(" ap1,ap2,am1,am2 = ",4e10.8)')  ap1,ap2,am1,am2
                                   ! write(*,'(" wp1,wpr,wm1,wm2 = ",4F10.8)')  wp1,wp2,wm1,wm2
                                   ! write(*,'(" bp1,bp2,bm1,bm2 = ",4F10.8)')  bp1,bp2,bm1,bm2
                                   ! write(*,'(" wplw,wprw,wmlw,wmrw = ",4F10.8)')  wplw,wprw,wmlw,wmrw
                                   ! write(*,'(" wpl,wpr,wml,wmr = ",4F10.8)')  wpl,wpr,wml,wmr
                                   ! write(*,'(" agu = ",F10.8)') agu(vard,var)
                           !     end if
                          ! u*ux-------------------------------------------------------------------------
                          !calculate smooth factor
                          ! Apl   = u(i1+1,i2,i3,uc) - 2.*u(i1  ,i2,i3,uc)   + u(i1-1,i2,i3,uc)
                          ! Bpl   = u(i1+1,i2,i3,uc) -    u(i1-1,i2,i3,uc)
                          ! Apr   = u(i1+2,i2,i3,uc) - 2.*u(i1+1,i2,i3,uc) + u(i1,i2,i3,uc)
                          ! Bpr   = u(i1+2,i2,i3,uc) -    u(i1  ,i2,i3,uc)
                          ! Aml   = u(i1,i2,i3,uc) - 2.*u(i1-1,i2,i3,uc) + u(i1-2,i2,i3,uc)
                          ! Bml   = u(i1,i2,i3,uc) -    u(i1-2,i2,i3,uc)
                          ! Amr   = Apl
                          ! Bmr   = Bpl
                          ! betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 1./4.*Bpl**2
                          ! betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 1./4.*Bpr**2
                          ! betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 1./4.*Bml**2
                          ! betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 1./4.*Bmr**2
                          ! ep     = 1e-15
                          ! ap1    = (1./2.)/(ep + betapl)**(4)
                          ! ap2    = (1./2.)/(ep + betapr)**(4)
                          ! am1    = (1./2.)/(ep + betaml)**(4)
                          ! am2    = (1./2.)/(ep + betamr)**(4)
                          ! !calculate weights
                          ! wp1 = ap1/(ap1+ap2)
                          ! wp2 = ap2/(ap1+ap2)
                          ! wm1 = am1/(am1+am2)
                          ! wm2 = am2/(am1+am2)
                          ! !mapping of the weights
                          ! bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                          ! bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                          ! bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                          ! bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                          ! wplw = bp1/(bp1+bp2)
                          ! wprw = bp2/(bp1+bp2)
                          ! wmlw = bm1/(bm1+bm2)
                          ! wmrw = bm2/(bm1+bm2)
                          ! !get face value of u
                          ! up = au
                          ! um = au
                          ! !picking side for the weights
                          !  if (up.ge.0) then
                          !    wpl = max(wplw,wprw)
                          !    wpr = min(wplw,wprw)
                          ! else
                          !    wpl = min(wplw,wprw)
                          !    wpr = max(wplw,wprw)
                          ! end if
                          ! if (um.ge.0) then
                          !    wml = max(wmlw,wmrw)
                          !    wmr = min(wmlw,wmrw)
                          ! else
                          !    wml = min(wmlw,wmrw)
                          !    wmr = max(wmlw,wmrw)
                          ! end if
                          ! expAd = 0.
                          ! if (orderOneFix.eq.1) then 
                          !    if (augmentPlot.eq.1) then
                          !       agu(uc,13)=0
                          !    end if
                          !    if (aum*aup .le. 0) then
                          !       aumax = max(abs(au),abs(aup),abs(aum))
                          !       expAd = (1./12.)*aumax*( u(i1-2,i2,i3,uc) - 4.*u(i1-1,i2,i3,uc) + 6.*u(i1,i2,i3,uc) - 4.*u(i1+1,i2,i3,uc) +  u(i1+2,i2,i3,uc))
                          !       wpl = 0.5
                          !       wpr = 0.5
                          !       wml = 0.5
                          !       wmr = 0.5
                          !       ! plot au at where fix is turned on
                          !       if (augmentPlot.eq.1) then
                          !          agu(uc,13)=au
                          !       end if
                          !    end if
                          ! end if
                          ! wpl = 0.5
                          ! wpr = 0.5
                          ! wml = 0.5
                          ! wmr = 0.5
                          ! Fpl = 1./6.*(  -u(i1-1,i2,i3,uc) + 5.*u(i1,i2,i3,uc)   + 2.*u(i1+1,i2,i3,uc))
                          ! Fpr = 1./6.*( 2.*u(i1,i2,i3,uc)   + 5.*u(i1+1,i2,i3,uc) -   u(i1+2,i2,i3,uc))
                          ! Fml = 1./6.*(  -u(i1-2,i2,i3,uc) + 5.*u(i1-1,i2,i3,uc)   + 2.*u(i1,i2,i3,uc))
                          ! Fmr = 1./6.*( 2.*u(i1-1,i2,i3,uc) + 5.*u(i1,i2,i3,uc) -       u(i1+1,i2,i3,uc))
                          ! Fp  = wpl*Fpl + wpr*Fpr
                          ! Fm  = wml*Fml + wmr*Fmr
                          ! agu(uc,uc)= au*(Fp - Fm)/(dx(0)) + expAd*dx(0)**(-1)
                          ! if (augmentPlot.eq.1) then
                          !    agu(uc,3)=au
                          !    agu(uc,4)=wpl
                          !    agu(uc,5)=wpr
                          !    agu(uc,6)=wml
                          !    agu(uc,7)=wmr
                          ! end if
                          ! u*vx-------------------------------------------------------------------------
                          !calculate smooth factor
                          ! Apl   = u(i1+1,i2,i3,vc) - 2.*u(i1,i2,i3,vc)   + u(i1-1,i2,i3,vc)
                          ! Bpl   = u(i1+1,i2,i3,vc) -    u(i1-1,i2,i3,vc)
                          ! Apr   = u(i1+2,i2,i3,vc) - 2.*u(i1+1,i2,i3,vc) + u(i1,i2,i3,vc)
                          ! Bpr   = u(i1+2,i2,i3,vc) -    u(i1,i2,i3,vc)
                          ! Aml   = u(i1,i2,i3,vc) - 2.*u(i1-1,i2,i3,vc) + u(i1-2,i2,i3,vc)
                          ! Bml   = u(i1,i2,i3,vc) -    u(i1-2,i2,i3,vc)
                          ! Amr   = Apl
                          ! Bmr   = Bpl
                          ! betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 1./4.*Bpl**2
                          ! betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 1./4.*Bpr**2
                          ! betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 1./4.*Bml**2
                          ! betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 1./4.*Bmr**2
                          ! ep     = 1e-15
                          ! ap1    = (1./2.)/(ep + betapl)**(4)
                          ! ap2    = (1./2.)/(ep + betapr)**(4)
                          ! am1    = (1./2.)/(ep + betaml)**(4)
                          ! am2    = (1./2.)/(ep + betamr)**(4)
                          ! !calculate weights
                          ! wp1 = ap1/(ap1+ap2)
                          ! wp2 = ap2/(ap1+ap2)
                          ! wm1 = am1/(am1+am2)
                          ! wm2 = am2/(am1+am2)
                          ! !mapping of the weights
                          ! bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                          ! bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                          ! bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                          ! bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                          ! wplw = bp1/(bp1+bp2)
                          ! wprw = bp2/(bp1+bp2)
                          ! wmlw = bm1/(bm1+bm2)
                          ! wmrw = bm2/(bm1+bm2)
                          ! !get face value of u
                          ! ! up = 0.5*(u(i1  ,i2,i3,uc)+u(i1+1,i2,i3,uc))
                          ! ! um = 0.5*(u(i1-1,i2,i3,uc)+u(i1  ,i2,i3,uc))
                          ! up = au
                          ! um = au
                          ! !picking side for the weights
                          ! if (up.ge.0) then
                          !    wpl = max(wplw,wprw)
                          !    wpr = min(wplw,wprw)
                          ! else
                          !    wpl = min(wplw,wprw)
                          !    wpr = max(wplw,wprw)
                          ! end if
                          ! if (um.ge.0) then
                          !    wml = max(wmlw,wmrw)
                          !    wmr = min(wmlw,wmrw)
                          ! else
                          !    wml = min(wmlw,wmrw)
                          !    wmr = max(wmlw,wmrw)
                          ! end if
                          ! expAd = 0.
                          ! if (orderOneFix.eq.1) then
                          !    if (aum*aup .le. 0) then
                          !       aumax = max(abs(au),abs(aup),abs(aum))
                          !       expAd = (1./12.)*aumax*( u(i1-2,i2,i3,vc) - 4.*u(i1-1,i2,i3,vc) + 6.*u(i1,i2,i3,vc) - 4.*u(i1+1,i2,i3,vc) +  u(i1+2,i2,i3,vc))
                          !       wpl = 0.5
                          !       wpr = 0.5
                          !       wml = 0.5
                          !       wmr = 0.5
                          !    end if
                          ! end if
                          ! Fpl = 1./6.*(  -u(i1-1,i2,i3,vc) + 5.*u(i1,i2,i3,vc)   + 2.*u(i1+1,i2,i3,vc))
                          ! Fpr = 1./6.*( 2.*u(i1,i2,i3,vc)   + 5.*u(i1+1,i2,i3,vc) -   u(i1+2,i2,i3,vc))
                          ! Fml = 1./6.*(  -u(i1-2,i2,i3,vc) + 5.*u(i1-1,i2,i3,vc)   + 2.*u(i1,i2,i3,vc))
                          ! Fmr = 1./6.*( 2.*u(i1-1,i2,i3,vc) + 5.*u(i1,i2,i3,vc) -       u(i1+1,i2,i3,vc))
                          ! Fp  = wpl*Fpl + wpr*Fpr
                          ! Fm  = wml*Fml + wmr*Fmr
                          ! agu(uc,vc)= au*(Fp - Fm)/(dx(0)) + expAd*dx(0)**(-1) 
                          ! if (augmentPlot .eq. 1) then 
                          !    agu(vc,4)=wpl
                          !    agu(vc,5)=wpr
                          !    agu(vc,6)=wml
                          !    agu(vc,7)=wmr
                          ! end if
                          !v*uy-------------------------------------------------------------------------
                          ! au = u(i1,i2,i3,vc)
                          ! aum = u(i1,i2-1,i3,vc)
                          ! aup = u(i1,i2+1,i3,vc)
                          !calculate smooth factor
                          ! Apl   = u(i1,i2+1,i3,uc) - 2.*u(i1,i2,i3,uc)   + u(i1,i2-1,i3,uc)
                          ! Bpl   = u(i1,i2+1,i3,uc) -    u(i1,i2-1,i3,uc)
                          ! Apr   = u(i1,i2+2,i3,uc) - 2.*u(i1,i2+1,i3,uc) + u(i1,i2,i3,uc)
                          ! Bpr   = u(i1,i2+2,i3,uc) -    u(i1,i2,i3,uc)
                          ! Aml   = u(i1,i2,i3,uc) - 2.*u(i1,i2-1,i3,uc) + u(i1,i2-2,i3,uc)
                          ! Bml   = u(i1,i2,i3,uc) -    u(i1,i2-2,i3,uc)
                          ! Amr   = Apl
                          ! Bmr   = Bpl
                          ! betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 1./4.*Bpl**2
                          ! betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 1./4.*Bpr**2
                          ! betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 1./4.*Bml**2
                          ! betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 1./4.*Bmr**2
                          ! ep     = 1e-15
                          ! ap1    = (1./2.)/(ep + betapl)**(4)
                          ! ap2    = (1./2.)/(ep + betapr)**(4)
                          ! am1    = (1./2.)/(ep + betaml)**(4)
                          ! am2    = (1./2.)/(ep + betamr)**(4)
                          ! !calculate weights
                          ! wp1 = ap1/(ap1+ap2)
                          ! wp2 = ap2/(ap1+ap2)
                          ! wm1 = am1/(am1+am2)
                          ! wm2 = am2/(am1+am2)
                          ! !mapping of the weights
                          ! bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                          ! bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                          ! bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                          ! bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                          ! wplw = bp1/(bp1+bp2)
                          ! wprw = bp2/(bp1+bp2)
                          ! wmlw = bm1/(bm1+bm2)
                          ! wmrw = bm2/(bm1+bm2)
                          ! !get face value of v
                          ! ! up = 0.5*(u(i1  ,i2+1,i3,vc)+u(i1,i2 ,i3,vc))
                          ! ! um = 0.5*(u(i1  ,i2  ,i3,vc)+u(i1,i2-1,i3,vc))
                          ! up = au
                          ! um = au
                          ! !picking side for the weights
                          !  if (up.ge.0) then
                          !    wpl = max(wplw,wprw)
                          !    wpr = min(wplw,wprw)
                          ! else
                          !    wpl = min(wplw,wprw)
                          !    wpr = max(wplw,wprw)
                          ! end if
                          ! if (um.ge.0) then
                          !    wml = max(wmlw,wmrw)
                          !    wmr = min(wmlw,wmrw)
                          ! else
                          !    wml = min(wmlw,wmrw)
                          !    wmr = max(wmlw,wmrw)
                          ! end if
                          ! expAd = 0.
                          ! if (orderOneFix.eq.1) then 
                          !    if (augmentPlot .eq. 1) then 
                          !       agu(vc,13)=0
                          !    end if
                          !    if (aum*aup .le. 0) then
                          !       aumax = max(abs(au),abs(aup),abs(aum))
                          !       expAd = (1./12.)*aumax*( u(i1,i2-2,i3,uc) - 4.*u(i1,i2-1,i3,uc) + 6.*u(i1,i2,i3,uc) - 4.*u(i1,i2+1,i3,uc) +  u(i1,i2+2,i3,uc))
                          !       wpl = 0.5
                          !       wpr = 0.5
                          !       wml = 0.5
                          !       wmr = 0.5
                          !       if (augmentPlot .eq. 1) then 
                          !          agu(vc,13)=au
                          !       end if
                          !    end if
                          ! end if
                          ! Fpl = 1./6.*(  -u(i1,i2-1,i3,uc) + 5.*u(i1,i2,i3,uc)   + 2.*u(i1,i2+1,i3,uc))
                          ! Fpr = 1./6.*( 2.*u(i1,i2,i3,uc)   + 5.*u(i1,i2+1,i3,uc) -   u(i1,i2+2,i3,uc))
                          ! Fml = 1./6.*(  -u(i1,i2-2,i3,uc) + 5.*u(i1,i2-1,i3,uc)   + 2.*u(i1,i2,i3,uc))
                          ! Fmr = 1./6.*( 2.*u(i1,i2-1,i3,uc) + 5.*u(i1,i2,i3,uc) -       u(i1,i2+1,i3,uc))
                          ! Fp  = wpl*Fpl + wpr*Fpr
                          ! Fm  = wml*Fml + wmr*Fmr
                          ! agu(vc,uc)= au*(Fp - Fm)/(dx(1)) + expAd*dx(1)**(-1)
                          ! if (augmentPlot .eq. 1) then 
                          !    agu(uc,8)=au
                          !    agu(uc,9) =wpl
                          !    agu(uc,10)=wpr
                          !    agu(uc,11)=wml
                          !    agu(uc,12)=wmr
                          ! end if
                          !v*vy---------------------------------------------------------------
                          ! Apl   = u(i1,i2+1,i3,vc) - 2.*u(i1,i2,i3,vc)   + u(i1,i2-1,i3,vc)
                          ! Bpl   = u(i1,i2+1,i3,vc) -    u(i1,i2-1,i3,vc)
                          ! Apr   = u(i1,i2+2,i3,vc) - 2.*u(i1,i2+1,i3,vc) + u(i1,i2,i3,vc)
                          ! Bpr   = u(i1,i2+2,i3,vc) -    u(i1,i2,i3,vc)
                          ! Aml   = u(i1,i2,i3,vc) - 2.*u(i1,i2-1,i3,vc) + u(i1,i2-2,i3,vc)
                          ! Bml   = u(i1,i2,i3,vc) -    u(i1,i2-2,i3,vc)
                          ! Amr   = Apl
                          ! Bmr   = Bpl
                          ! betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 1./4.*Bpl**2
                          ! betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 1./4.*Bpr**2
                          ! betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 1./4.*Bml**2
                          ! betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 1./4.*Bmr**2
                          ! ep     = 1e-15
                          ! ap1    = (1./2.)/(ep + betapl)**(4)
                          ! ap2    = (1./2.)/(ep + betapr)**(4)
                          ! am1    = (1./2.)/(ep + betaml)**(4)
                          ! am2    = (1./2.)/(ep + betamr)**(4)
                          ! !calculate weights
                          ! wp1 = ap1/(ap1+ap2)
                          ! wp2 = ap2/(ap1+ap2)
                          ! wm1 = am1/(am1+am2)
                          ! wm2 = am2/(am1+am2)
                          ! !mapping of the weights
                          ! bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                          ! bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                          ! bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                          ! bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                          ! wplw = bp1/(bp1+bp2)
                          ! wprw = bp2/(bp1+bp2)
                          ! wmlw = bm1/(bm1+bm2)
                          ! wmrw = bm2/(bm1+bm2)
                          ! !get face value of v
                          !  ! up = 0.5*(u(i1  ,i2,i3,vc)  + u(i1,i2+1,i3,vc))
                          !  ! um = 0.5*(u(i1-1,i2,i3,vc)  + u(i1,i2  ,i3,vc))
                          ! up = au
                          ! um = au
                          ! !picking side for the weights
                          !  if (up.ge.0) then
                          !    wpl = max(wplw,wprw)
                          !    wpr = min(wplw,wprw)
                          ! else
                          !    wpl = min(wplw,wprw)
                          !    wpr = max(wplw,wprw)
                          ! end if
                          ! if (um.ge.0) then
                          !    wml = max(wmlw,wmrw)
                          !    wmr = min(wmlw,wmrw)
                          ! else
                          !    wml = min(wmlw,wmrw)
                          !    wmr = max(wmlw,wmrw)
                          ! end if
                          ! expAd=0.
                          ! if (orderOneFix.eq.1) then
                          !    if (aum*aup .le. 0) then
                          !       aumax = max(abs(au),abs(aup),abs(aum))
                          !       expAd = (1./12.)*aumax*( u(i1,i2-2,i3,vc) - 4.*u(i1,i2-1,i3,vc) + 6.*u(i1,i2,i3,vc) - 4.*u(i1,i2+1,i3,vc) +  u(i1,i2+2,i3,vc))
                          !       wpl = 0.5
                          !       wpr = 0.5
                          !       wml = 0.5
                          !       wmr = 0.5
                          !    end if
                          ! end if
                          ! Fpl = 1./6.*(  -u(i1,i2-1,i3,vc) + 5.*u(i1,i2,i3,vc)   + 2.*u(i1,i2+1,i3,vc))
                          ! Fpr = 1./6.*( 2.*u(i1,i2,i3,vc)   + 5.*u(i1,i2+1,i3,vc) -   u(i1,i2+2,i3,vc))
                          ! Fml = 1./6.*(  -u(i1,i2-2,i3,vc) + 5.*u(i1,i2-1,i3,vc)   + 2.*u(i1,i2,i3,vc))
                          ! Fmr = 1./6.*( 2.*u(i1,i2-1,i3,vc) + 5.*u(i1,i2,i3,vc) -       u(i1,i2+1,i3,vc))
                          ! Fp  = wpl*Fpl + wpr*Fpr
                          ! Fm  = wml*Fml + wmr*Fmr
                          ! agu(vc,vc)= au*(Fp - Fm)/(dx(1)) + expAd*dx(1)**(-1)
                          ! if (augmentPlot .eq. 1) then  
                          !    agu(vc,9) =wpl
                          !    agu(vc,10)=wpr
                          !    agu(vc,11)=wml
                          !    agu(vc,12)=wmr
                          ! end if
                       else ! if (upwindOrder.eq.4)
                          write(*,'(" getBewnoAdvection: only 4th 
     & order is avaliable now")' )
                          stop 777
                       end if
                   ! #If ("BWENO" == "CENTERED") 
                   !  ! -- centered advection ---
                   !  ! write(*,'(" getAdvection -- centered")')
                   !  #If "2" == "2"
                   !   ! -- 2D --
                   !   agu(uc,uc)=UU(uc)*UX(uc)
                   !   agu(vc,uc)=UU(vc)*UY(uc)
                   !   agu(uc,vc)=UU(uc)*UX(vc)
                   !   agu(vc,vc)=UU(vc)*UY(vc)
                   !  #Else
                   !   ! -- 3D -- *check me*
                   !   agu(uc,uc)=UU(uc)*UX(uc)
                   !   agu(vc,uc)=UU(vc)*UY(uc)
                   !   agu(wc,uc)=UU(wc)*UY(uc)
                   !   agu(uc,vc)=UU(uc)*UX(vc)
                   !   agu(vc,vc)=UU(vc)*UY(vc)
                   !   agu(wc,vc)=UU(wc)*UY(vc)
                   !   agu(uc,wc)=UU(uc)*UX(wc)
                   !   agu(vc,wc)=UU(vc)*UY(wc)
                   !   agu(wc,wc)=UU(wc)*UY(wc)
                   !  #End
                   ! #Elif "BWENO" == "BWENO" 
                   !   ! --- upwind scheme ---
                   !   ! for testing output this next message:
                   !   if( t.le. 0. )then
                   !     write(*,'(" getAdvection upwind scheme (7)")') 
                   !   end if
                   !   getUpwindAdvection(u,i1,i2,i3,NONE,2,4,rectangular, agu)
                   ! #Elif "BWENO" == "BWENO" 
                   !   ! --- Bweno scheme ---
                   !   ! for testing output this message:
                   !   if( t.le. 0. )then
                   !      write(*,'(" getAdvection BWENO scheme (7)")') 
                   !   end if
                   !   getBwenoAdvection(u,i1,i2,i3,NONE,2,4,rectangular, agu)
                   ! #Else
                   !   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
                   !   stop 999
                   ! #End
                      ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-ux42r(
     & i1,i2,i3,pc)+nu*ulaplacian42r(i1,i2,i3,uc)
                      ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-uy42r(
     & i1,i2,i3,pc)+nu*ulaplacian42r(i1,i2,i3,vc)
                !       u(i1,i2,i3,3)= agu(uc,3) ! au in r direction 
                !       u(i1,i2,i3,4)= agu(uc,8) ! au in s direction
                !       u(i1,i2,i3,5)= agu(uc,4)  ! wpl uur
                !       u(i1,i2,i3,6)= agu(uc,5)  ! wpr uur
                !       u(i1,i2,i3,7)= agu(uc,6)  ! wml uur
                !       u(i1,i2,i3,8)= agu(uc,7)  ! wmr uur
                !       u(i1,i2,i3,9)= agu(uc,9)  ! wpl uus
                !       u(i1,i2,i3,10)= agu(uc,10)  ! wpr uus
                !       u(i1,i2,i3,11)= agu(uc,11)  ! wml uus
                !       u(i1,i2,i3,12)= agu(uc,12)  ! wmr uus
                ! !
                !       u(i1,i2,i3,13)= agu(vc,4)  ! wpl uvr
                !       u(i1,i2,i3,14)= agu(vc,5)  ! wpr uvr
                !       u(i1,i2,i3,15)= agu(vc,6)  ! wml uvr
                !       u(i1,i2,i3,16)= agu(vc,7)  ! wmr uvr
                !       u(i1,i2,i3,17)= agu(vc,9)  ! wpl uvs
                !       u(i1,i2,i3,18)= agu(vc,10)  ! wpr uvs
                !       u(i1,i2,i3,19)= agu(vc,11)  ! wml uvs
                !       u(i1,i2,i3,20)= agu(vc,12)  ! wmr uvs
                !       u(i1,i2,i3,21)= agu(uc,13)  ! au at where fix is turned on 
                !       u(i1,i2,i3,22)= agu(vc,13)  ! au at where fix is turned on 
               end do
               end do
               end do
              end if
            else
              write(*,'(" unknown advectionOption")')
              stop 1010
            end if
           else ! gridIsImplicit
            ! ---- implicit time-stepping ---
            if( advectionOption.eq.centeredAdvection )then
             if( implicitOption .eq.computeImplicitTermsSeparately )
     & then
               if( useWhereMask.ne.0 )then
                do i3=n3a,n3b
                do i2=n2a,n2b
                do i1=n1a,n1b
                 if( mask(i1,i2,i3).gt.0 )then
                    ! explicit terms only, no diffusion
                     ut(i1,i2,i3,uc)= -uu(i1,i2,i3,uc)*ux42r(i1,i2,i3,
     & uc)-uu(i1,i2,i3,vc)*uy42r(i1,i2,i3,uc)-ux42r(i1,i2,i3,pc)
                     ut(i1,i2,i3,vc)= -uu(i1,i2,i3,uc)*ux42r(i1,i2,i3,
     & vc)-uu(i1,i2,i3,vc)*uy42r(i1,i2,i3,vc)-uy42r(i1,i2,i3,pc)
                    ! include implicit terms - diffusion
                     uti(i1,i2,i3,uc)= nu*ulaplacian42r(i1,i2,i3,uc)
                     uti(i1,i2,i3,vc)= nu*ulaplacian42r(i1,i2,i3,vc)
                 end if
                end do
                end do
                end do
               else
                do i3=n3a,n3b
                do i2=n2a,n2b
                do i1=n1a,n1b
                   ! explicit terms only, no diffusion
                    ut(i1,i2,i3,uc)= -uu(i1,i2,i3,uc)*ux42r(i1,i2,i3,
     & uc)-uu(i1,i2,i3,vc)*uy42r(i1,i2,i3,uc)-ux42r(i1,i2,i3,pc)
                    ut(i1,i2,i3,vc)= -uu(i1,i2,i3,uc)*ux42r(i1,i2,i3,
     & vc)-uu(i1,i2,i3,vc)*uy42r(i1,i2,i3,vc)-uy42r(i1,i2,i3,pc)
                   ! include implicit terms - diffusion
                    uti(i1,i2,i3,uc)= nu*ulaplacian42r(i1,i2,i3,uc)
                    uti(i1,i2,i3,vc)= nu*ulaplacian42r(i1,i2,i3,vc)
                end do
                end do
                end do
               end if
             else if( implicitOption.eq.doNotComputeImplicitTerms )then
               if( useWhereMask.ne.0 )then
                do i3=n3a,n3b
                do i2=n2a,n2b
                do i1=n1a,n1b
                 if( mask(i1,i2,i3).gt.0 )then
                    ! explicit terms only, no diffusion
                     ut(i1,i2,i3,uc)= -uu(i1,i2,i3,uc)*ux42r(i1,i2,i3,
     & uc)-uu(i1,i2,i3,vc)*uy42r(i1,i2,i3,uc)-ux42r(i1,i2,i3,pc)
                     ut(i1,i2,i3,vc)= -uu(i1,i2,i3,uc)*ux42r(i1,i2,i3,
     & vc)-uu(i1,i2,i3,vc)*uy42r(i1,i2,i3,vc)-uy42r(i1,i2,i3,pc)
                 end if
                end do
                end do
                end do
               else
                do i3=n3a,n3b
                do i2=n2a,n2b
                do i1=n1a,n1b
                   ! explicit terms only, no diffusion
                    ut(i1,i2,i3,uc)= -uu(i1,i2,i3,uc)*ux42r(i1,i2,i3,
     & uc)-uu(i1,i2,i3,vc)*uy42r(i1,i2,i3,uc)-ux42r(i1,i2,i3,pc)
                    ut(i1,i2,i3,vc)= -uu(i1,i2,i3,uc)*ux42r(i1,i2,i3,
     & vc)-uu(i1,i2,i3,vc)*uy42r(i1,i2,i3,vc)-uy42r(i1,i2,i3,pc)
                end do
                end do
                end do
               end if
             else
              write(*,*)'insdt: Unknown implicitOption=',implicitOption
              stop 5
             end if  ! end implicitOption
            else if( advectionOption.eq.upwindAdvection )then
              ! --- upwind ---
             if( implicitOption .eq.computeImplicitTermsSeparately )
     & then
               if( useWhereMask.ne.0 )then
                do i3=n3a,n3b
                do i2=n2a,n2b
                do i1=n1a,n1b
                 if( mask(i1,i2,i3).gt.0 )then
                    ! explicit terms only, no diffusion
                      ! upwind approximation/bweno
                        ! --- upwind scheme ---
                        ! for testing output this next message:
                        ! if( t.le. 0. )then
                        !    write(*,'(" getAdvection upwind scheme (7)")') 
                        ! end if
                          ! --- CARTESIAN GRID ---
                          !- agu(uc,uc)=UU(uc)*UX(uc)
                          !- agu(vc,uc)=UU(vc)*UY(uc)
                          !- agu(uc,vc)=UU(uc)*UX(vc)
                          !- agu(vc,vc)=UU(vc)*UY(vc)
                          ! -- first order upwind --
                          if( upwindOrder.eq.1 )then
                           au = uu(i1,i2,i3,uc)
                           if( au.gt.0. )then
                             ! u*ux = u*D-x(u)
                             agu(uc,uc)= au*(u(i1,i2,i3,uc)-u(i1-1,i2,
     & i3,uc))/(dx(0))
                             ! u*vx = u*D-x(v)
                             agu(uc,vc)= au*(u(i1,i2,i3,vc)-u(i1-1,i2,
     & i3,vc))/(dx(0))
                           else
                             ! u*ux = u*D+x(u)
                             agu(uc,uc)= au*(u(i1+1,i2,i3,uc)-u(i1,i2,
     & i3,uc))/(dx(0))
                             ! u*vx = u*D+x(v)
                             agu(uc,vc)= au*(u(i1+1,i2,i3,vc)-u(i1,i2,
     & i3,vc))/(dx(0))
                           end if
                           au = uu(i1,i2,i3,vc)
                           if( au.gt.0. )then
                             ! v*uy = v*D-y(u)
                             agu(vc,uc)= au*(u(i1,i2,i3,uc)-u(i1,i2-1,
     & i3,uc))/(dx(1))
                             ! v*vy = v*D-y(v)
                             agu(vc,vc)= au*(u(i1,i2,i3,vc)-u(i1,i2-1,
     & i3,vc))/(dx(1))
                           else
                             ! v*uy = v*D+y(u)
                             agu(vc,uc)= au*(u(i1,i2+1,i3,uc)-u(i1,i2,
     & i3,uc))/(dx(1))
                             ! v*vy = v*D+y(v) 
                             agu(vc,vc)= au*(u(i1,i2+1,i3,vc)-u(i1,i2,
     & i3,vc))/(dx(1))
                           end if
                          elseif( upwindOrder.eq.2 )then
                            ! write(*,'(" finish me, upwindOrder=",i2)') upwindOrder
                            ! stop 222
                             au = uu(i1,i2,i3,uc)
                             if( au.gt.0. )then
                                ! u*ux = u*D-x(u)
                                agu(uc,uc)= au*(3.*u(i1,i2,i3,uc)-4.*u(
     & i1-1,i2,i3,uc)+u(i1-2,i2,i3,uc))/(2.*dx(0))
                                ! u*vx = u*D-x(v)
                                agu(uc,vc)= au*(3.*u(i1,i2,i3,vc)-4.*u(
     & i1-1,i2,i3,vc)+u(i1-2,i2,i3,vc))/(2.*dx(0))
                             else
                                ! u*ux = u*D+x(u)
                                agu(uc,uc)= au*(-u(i1+2,i2,i3,uc)+4.*u(
     & i1+1,i2,i3,uc)-3.*u(i1,i2,i3,uc))/(2.*dx(0))
                                ! u*vx = u*D+x(v)
                                agu(uc,vc)= au*(-u(i1+2,i2,i3,vc)+4.*u(
     & i1+1,i2,i3,vc)-3.*u(i1,i2,i3,vc))/(2.*dx(0))
                             end if
                             au = uu(i1,i2,i3,vc)
                             if( au.gt.0. )then
                                ! v*uy = v*D-y(u)
                                agu(vc,uc)= au*(3.*u(i1,i2,i3,uc)-4.*u(
     & i1,i2-1,i3,uc)+u(i1,i2-2,i3,uc))/(2.*dx(1))
                                ! v*vy = v*D-y(v)
                                agu(vc,vc)= au*(3.*u(i1,i2,i3,vc)-4.*u(
     & i1,i2-1,i3,vc)+u(i1,i2-2,i3,vc))/(2.*dx(1))
                             else
                                ! v*uy = v*D+y(u)
                                agu(vc,uc)= au*(-u(i1,i2+2,i3,uc)+4.*u(
     & i1,i2+1,i3,uc)-3.*u(i1,i2,i3,uc))/(2.*dx(1))
                                ! v*vy = v*D+y(v) 
                                agu(vc,vc)= au*(-u(i1,i2+2,i3,vc)+4.*u(
     & i1,i2+1,i3,vc)-3.*u(i1,i2,i3,vc))/(2.*dx(1))
                             end if
                          elseif( upwindOrder.eq.3 )then
                             au = uu(i1,i2,i3,uc)
                             if( au.gt.0. )then
                                ! u*ux = u*D-x(u)
                                agu(uc,uc)= au*(2.*u(i1+1,i2,i3,uc)+3.*
     & u(i1,i2,i3,uc)-6.*u(i1-1,i2,i3,uc)+u(i1-2,i2,i3,uc))/(6.*dx(0))
                                ! u*vx = u*D-x(v)
                                agu(uc,vc)= au*(2.*u(i1+1,i2,i3,vc)+3.*
     & u(i1,i2,i3,vc)-6.*u(i1-1,i2,i3,vc)+u(i1-2,i2,i3,vc))/(6.*dx(0))
                             else
                                ! u*ux = u*D+x(u)
                                agu(uc,uc)= au*(-u(i1+2,i2,i3,uc)+6.*u(
     & i1+1,i2,i3,uc)-3.*u(i1,i2,i3,uc)-2.*u(i1-1,i2,i3,uc))/(6.*dx(0)
     & )
                                ! u*vx = u*D+x(v)
                                agu(uc,vc)= au*(-u(i1+2,i2,i3,vc)+6.*u(
     & i1+1,i2,i3,vc)-3.*u(i1,i2,i3,vc)-2.*u(i1-1,i2,i3,vc))/(6.*dx(0)
     & )
                             end if
                             au = uu(i1,i2,i3,vc)
                             if( au.gt.0. )then
                                ! v*uy = v*D-y(u)
                                agu(vc,uc)= au*(2.*u(i1,i2+1,i3,uc)+3.*
     & u(i1,i2,i3,uc)-6.*u(i1,i2-1,i3,uc)+u(i1,i2-2,i3,uc))/(6.*dx(1))
                                ! v*vy = v*D-y(v)
                                agu(vc,vc)= au*(2.*u(i1,i2+1,i3,vc)+3.*
     & u(i1,i2,i3,vc)-6.*u(i1,i2-1,i3,vc)+u(i1,i2-2,i3,vc))/(6.*dx(1))
                             else
                                ! v*uy = v*D+y(u)
                                agu(vc,uc)= au*(-u(i1,i2+2,i3,uc)+6.*u(
     & i1,i2+1,i3,uc)-3.*u(i1,i2,i3,uc)-2.*u(i1,i2-1,i3,uc))/(6.*dx(1)
     & )
                                ! v*vy = v*D+y(v) 
                                agu(vc,vc)= au*(-u(i1,i2+2,i3,vc)+6.*u(
     & i1,i2+1,i3,vc)-3.*u(i1,i2,i3,vc)-2.*u(i1,i2-1,i3,vc))/(6.*dx(1)
     & )
                             end if
                          end if
                      ! #If ("UPWIND" == "CENTERED") 
                      !  ! -- centered advection ---
                      !  ! write(*,'(" getAdvection -- centered")')
                      !  #If "2" == "2"
                      !   ! -- 2D --
                      !   agu(uc,uc)=UU(uc)*UX(uc)
                      !   agu(vc,uc)=UU(vc)*UY(uc)
                      !   agu(uc,vc)=UU(uc)*UX(vc)
                      !   agu(vc,vc)=UU(vc)*UY(vc)
                      !  #Else
                      !   ! -- 3D -- *check me*
                      !   agu(uc,uc)=UU(uc)*UX(uc)
                      !   agu(vc,uc)=UU(vc)*UY(uc)
                      !   agu(wc,uc)=UU(wc)*UY(uc)
                      !   agu(uc,vc)=UU(uc)*UX(vc)
                      !   agu(vc,vc)=UU(vc)*UY(vc)
                      !   agu(wc,vc)=UU(wc)*UY(vc)
                      !   agu(uc,wc)=UU(uc)*UX(wc)
                      !   agu(vc,wc)=UU(vc)*UY(wc)
                      !   agu(wc,wc)=UU(wc)*UY(wc)
                      !  #End
                      ! #Elif "UPWIND" == "UPWIND" 
                      !   ! --- upwind scheme ---
                      !   ! for testing output this next message:
                      !   if( t.le. 0. )then
                      !     write(*,'(" getAdvection upwind scheme (7)")') 
                      !   end if
                      !   getUpwindAdvection(u,i1,i2,i3,NONE,2,4,rectangular, agu)
                      ! #Elif "UPWIND" == "BWENO" 
                      !   ! --- Bweno scheme ---
                      !   ! for testing output this message:
                      !   if( t.le. 0. )then
                      !      write(*,'(" getAdvection BWENO scheme (7)")') 
                      !   end if
                      !   getBwenoAdvection(u,i1,i2,i3,NONE,2,4,rectangular, agu)
                      ! #Else
                      !   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
                      !   stop 999
                      ! #End
                        ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-
     & ux42r(i1,i2,i3,pc)
                        ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-
     & uy42r(i1,i2,i3,pc)
                    ! include implicit terms - diffusion
                     uti(i1,i2,i3,uc)= nu*ulaplacian42r(i1,i2,i3,uc)
                     uti(i1,i2,i3,vc)= nu*ulaplacian42r(i1,i2,i3,vc)
                 end if
                end do
                end do
                end do
               else
                do i3=n3a,n3b
                do i2=n2a,n2b
                do i1=n1a,n1b
                   ! explicit terms only, no diffusion
                     ! upwind approximation/bweno
                       ! --- upwind scheme ---
                       ! for testing output this next message:
                       ! if( t.le. 0. )then
                       !    write(*,'(" getAdvection upwind scheme (7)")') 
                       ! end if
                         ! --- CARTESIAN GRID ---
                         !- agu(uc,uc)=UU(uc)*UX(uc)
                         !- agu(vc,uc)=UU(vc)*UY(uc)
                         !- agu(uc,vc)=UU(uc)*UX(vc)
                         !- agu(vc,vc)=UU(vc)*UY(vc)
                         ! -- first order upwind --
                         if( upwindOrder.eq.1 )then
                          au = uu(i1,i2,i3,uc)
                          if( au.gt.0. )then
                            ! u*ux = u*D-x(u)
                            agu(uc,uc)= au*(u(i1,i2,i3,uc)-u(i1-1,i2,
     & i3,uc))/(dx(0))
                            ! u*vx = u*D-x(v)
                            agu(uc,vc)= au*(u(i1,i2,i3,vc)-u(i1-1,i2,
     & i3,vc))/(dx(0))
                          else
                            ! u*ux = u*D+x(u)
                            agu(uc,uc)= au*(u(i1+1,i2,i3,uc)-u(i1,i2,
     & i3,uc))/(dx(0))
                            ! u*vx = u*D+x(v)
                            agu(uc,vc)= au*(u(i1+1,i2,i3,vc)-u(i1,i2,
     & i3,vc))/(dx(0))
                          end if
                          au = uu(i1,i2,i3,vc)
                          if( au.gt.0. )then
                            ! v*uy = v*D-y(u)
                            agu(vc,uc)= au*(u(i1,i2,i3,uc)-u(i1,i2-1,
     & i3,uc))/(dx(1))
                            ! v*vy = v*D-y(v)
                            agu(vc,vc)= au*(u(i1,i2,i3,vc)-u(i1,i2-1,
     & i3,vc))/(dx(1))
                          else
                            ! v*uy = v*D+y(u)
                            agu(vc,uc)= au*(u(i1,i2+1,i3,uc)-u(i1,i2,
     & i3,uc))/(dx(1))
                            ! v*vy = v*D+y(v) 
                            agu(vc,vc)= au*(u(i1,i2+1,i3,vc)-u(i1,i2,
     & i3,vc))/(dx(1))
                          end if
                         elseif( upwindOrder.eq.2 )then
                           ! write(*,'(" finish me, upwindOrder=",i2)') upwindOrder
                           ! stop 222
                            au = uu(i1,i2,i3,uc)
                            if( au.gt.0. )then
                               ! u*ux = u*D-x(u)
                               agu(uc,uc)= au*(3.*u(i1,i2,i3,uc)-4.*u(
     & i1-1,i2,i3,uc)+u(i1-2,i2,i3,uc))/(2.*dx(0))
                               ! u*vx = u*D-x(v)
                               agu(uc,vc)= au*(3.*u(i1,i2,i3,vc)-4.*u(
     & i1-1,i2,i3,vc)+u(i1-2,i2,i3,vc))/(2.*dx(0))
                            else
                               ! u*ux = u*D+x(u)
                               agu(uc,uc)= au*(-u(i1+2,i2,i3,uc)+4.*u(
     & i1+1,i2,i3,uc)-3.*u(i1,i2,i3,uc))/(2.*dx(0))
                               ! u*vx = u*D+x(v)
                               agu(uc,vc)= au*(-u(i1+2,i2,i3,vc)+4.*u(
     & i1+1,i2,i3,vc)-3.*u(i1,i2,i3,vc))/(2.*dx(0))
                            end if
                            au = uu(i1,i2,i3,vc)
                            if( au.gt.0. )then
                               ! v*uy = v*D-y(u)
                               agu(vc,uc)= au*(3.*u(i1,i2,i3,uc)-4.*u(
     & i1,i2-1,i3,uc)+u(i1,i2-2,i3,uc))/(2.*dx(1))
                               ! v*vy = v*D-y(v)
                               agu(vc,vc)= au*(3.*u(i1,i2,i3,vc)-4.*u(
     & i1,i2-1,i3,vc)+u(i1,i2-2,i3,vc))/(2.*dx(1))
                            else
                               ! v*uy = v*D+y(u)
                               agu(vc,uc)= au*(-u(i1,i2+2,i3,uc)+4.*u(
     & i1,i2+1,i3,uc)-3.*u(i1,i2,i3,uc))/(2.*dx(1))
                               ! v*vy = v*D+y(v) 
                               agu(vc,vc)= au*(-u(i1,i2+2,i3,vc)+4.*u(
     & i1,i2+1,i3,vc)-3.*u(i1,i2,i3,vc))/(2.*dx(1))
                            end if
                         elseif( upwindOrder.eq.3 )then
                            au = uu(i1,i2,i3,uc)
                            if( au.gt.0. )then
                               ! u*ux = u*D-x(u)
                               agu(uc,uc)= au*(2.*u(i1+1,i2,i3,uc)+3.*
     & u(i1,i2,i3,uc)-6.*u(i1-1,i2,i3,uc)+u(i1-2,i2,i3,uc))/(6.*dx(0))
                               ! u*vx = u*D-x(v)
                               agu(uc,vc)= au*(2.*u(i1+1,i2,i3,vc)+3.*
     & u(i1,i2,i3,vc)-6.*u(i1-1,i2,i3,vc)+u(i1-2,i2,i3,vc))/(6.*dx(0))
                            else
                               ! u*ux = u*D+x(u)
                               agu(uc,uc)= au*(-u(i1+2,i2,i3,uc)+6.*u(
     & i1+1,i2,i3,uc)-3.*u(i1,i2,i3,uc)-2.*u(i1-1,i2,i3,uc))/(6.*dx(0)
     & )
                               ! u*vx = u*D+x(v)
                               agu(uc,vc)= au*(-u(i1+2,i2,i3,vc)+6.*u(
     & i1+1,i2,i3,vc)-3.*u(i1,i2,i3,vc)-2.*u(i1-1,i2,i3,vc))/(6.*dx(0)
     & )
                            end if
                            au = uu(i1,i2,i3,vc)
                            if( au.gt.0. )then
                               ! v*uy = v*D-y(u)
                               agu(vc,uc)= au*(2.*u(i1,i2+1,i3,uc)+3.*
     & u(i1,i2,i3,uc)-6.*u(i1,i2-1,i3,uc)+u(i1,i2-2,i3,uc))/(6.*dx(1))
                               ! v*vy = v*D-y(v)
                               agu(vc,vc)= au*(2.*u(i1,i2+1,i3,vc)+3.*
     & u(i1,i2,i3,vc)-6.*u(i1,i2-1,i3,vc)+u(i1,i2-2,i3,vc))/(6.*dx(1))
                            else
                               ! v*uy = v*D+y(u)
                               agu(vc,uc)= au*(-u(i1,i2+2,i3,uc)+6.*u(
     & i1,i2+1,i3,uc)-3.*u(i1,i2,i3,uc)-2.*u(i1,i2-1,i3,uc))/(6.*dx(1)
     & )
                               ! v*vy = v*D+y(v) 
                               agu(vc,vc)= au*(-u(i1,i2+2,i3,vc)+6.*u(
     & i1,i2+1,i3,vc)-3.*u(i1,i2,i3,vc)-2.*u(i1,i2-1,i3,vc))/(6.*dx(1)
     & )
                            end if
                         end if
                     ! #If ("UPWIND" == "CENTERED") 
                     !  ! -- centered advection ---
                     !  ! write(*,'(" getAdvection -- centered")')
                     !  #If "2" == "2"
                     !   ! -- 2D --
                     !   agu(uc,uc)=UU(uc)*UX(uc)
                     !   agu(vc,uc)=UU(vc)*UY(uc)
                     !   agu(uc,vc)=UU(uc)*UX(vc)
                     !   agu(vc,vc)=UU(vc)*UY(vc)
                     !  #Else
                     !   ! -- 3D -- *check me*
                     !   agu(uc,uc)=UU(uc)*UX(uc)
                     !   agu(vc,uc)=UU(vc)*UY(uc)
                     !   agu(wc,uc)=UU(wc)*UY(uc)
                     !   agu(uc,vc)=UU(uc)*UX(vc)
                     !   agu(vc,vc)=UU(vc)*UY(vc)
                     !   agu(wc,vc)=UU(wc)*UY(vc)
                     !   agu(uc,wc)=UU(uc)*UX(wc)
                     !   agu(vc,wc)=UU(vc)*UY(wc)
                     !   agu(wc,wc)=UU(wc)*UY(wc)
                     !  #End
                     ! #Elif "UPWIND" == "UPWIND" 
                     !   ! --- upwind scheme ---
                     !   ! for testing output this next message:
                     !   if( t.le. 0. )then
                     !     write(*,'(" getAdvection upwind scheme (7)")') 
                     !   end if
                     !   getUpwindAdvection(u,i1,i2,i3,NONE,2,4,rectangular, agu)
                     ! #Elif "UPWIND" == "BWENO" 
                     !   ! --- Bweno scheme ---
                     !   ! for testing output this message:
                     !   if( t.le. 0. )then
                     !      write(*,'(" getAdvection BWENO scheme (7)")') 
                     !   end if
                     !   getBwenoAdvection(u,i1,i2,i3,NONE,2,4,rectangular, agu)
                     ! #Else
                     !   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
                     !   stop 999
                     ! #End
                       ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-ux42r(
     & i1,i2,i3,pc)
                       ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-uy42r(
     & i1,i2,i3,pc)
                   ! include implicit terms - diffusion
                    uti(i1,i2,i3,uc)= nu*ulaplacian42r(i1,i2,i3,uc)
                    uti(i1,i2,i3,vc)= nu*ulaplacian42r(i1,i2,i3,vc)
                end do
                end do
                end do
               end if
             else if( implicitOption.eq.doNotComputeImplicitTerms )then
               if( useWhereMask.ne.0 )then
                do i3=n3a,n3b
                do i2=n2a,n2b
                do i1=n1a,n1b
                 if( mask(i1,i2,i3).gt.0 )then
                    ! explicit terms only, no diffusion
                      ! upwind approximation/bweno
                        ! --- upwind scheme ---
                        ! for testing output this next message:
                        ! if( t.le. 0. )then
                        !    write(*,'(" getAdvection upwind scheme (7)")') 
                        ! end if
                          ! --- CARTESIAN GRID ---
                          !- agu(uc,uc)=UU(uc)*UX(uc)
                          !- agu(vc,uc)=UU(vc)*UY(uc)
                          !- agu(uc,vc)=UU(uc)*UX(vc)
                          !- agu(vc,vc)=UU(vc)*UY(vc)
                          ! -- first order upwind --
                          if( upwindOrder.eq.1 )then
                           au = uu(i1,i2,i3,uc)
                           if( au.gt.0. )then
                             ! u*ux = u*D-x(u)
                             agu(uc,uc)= au*(u(i1,i2,i3,uc)-u(i1-1,i2,
     & i3,uc))/(dx(0))
                             ! u*vx = u*D-x(v)
                             agu(uc,vc)= au*(u(i1,i2,i3,vc)-u(i1-1,i2,
     & i3,vc))/(dx(0))
                           else
                             ! u*ux = u*D+x(u)
                             agu(uc,uc)= au*(u(i1+1,i2,i3,uc)-u(i1,i2,
     & i3,uc))/(dx(0))
                             ! u*vx = u*D+x(v)
                             agu(uc,vc)= au*(u(i1+1,i2,i3,vc)-u(i1,i2,
     & i3,vc))/(dx(0))
                           end if
                           au = uu(i1,i2,i3,vc)
                           if( au.gt.0. )then
                             ! v*uy = v*D-y(u)
                             agu(vc,uc)= au*(u(i1,i2,i3,uc)-u(i1,i2-1,
     & i3,uc))/(dx(1))
                             ! v*vy = v*D-y(v)
                             agu(vc,vc)= au*(u(i1,i2,i3,vc)-u(i1,i2-1,
     & i3,vc))/(dx(1))
                           else
                             ! v*uy = v*D+y(u)
                             agu(vc,uc)= au*(u(i1,i2+1,i3,uc)-u(i1,i2,
     & i3,uc))/(dx(1))
                             ! v*vy = v*D+y(v) 
                             agu(vc,vc)= au*(u(i1,i2+1,i3,vc)-u(i1,i2,
     & i3,vc))/(dx(1))
                           end if
                          elseif( upwindOrder.eq.2 )then
                            ! write(*,'(" finish me, upwindOrder=",i2)') upwindOrder
                            ! stop 222
                             au = uu(i1,i2,i3,uc)
                             if( au.gt.0. )then
                                ! u*ux = u*D-x(u)
                                agu(uc,uc)= au*(3.*u(i1,i2,i3,uc)-4.*u(
     & i1-1,i2,i3,uc)+u(i1-2,i2,i3,uc))/(2.*dx(0))
                                ! u*vx = u*D-x(v)
                                agu(uc,vc)= au*(3.*u(i1,i2,i3,vc)-4.*u(
     & i1-1,i2,i3,vc)+u(i1-2,i2,i3,vc))/(2.*dx(0))
                             else
                                ! u*ux = u*D+x(u)
                                agu(uc,uc)= au*(-u(i1+2,i2,i3,uc)+4.*u(
     & i1+1,i2,i3,uc)-3.*u(i1,i2,i3,uc))/(2.*dx(0))
                                ! u*vx = u*D+x(v)
                                agu(uc,vc)= au*(-u(i1+2,i2,i3,vc)+4.*u(
     & i1+1,i2,i3,vc)-3.*u(i1,i2,i3,vc))/(2.*dx(0))
                             end if
                             au = uu(i1,i2,i3,vc)
                             if( au.gt.0. )then
                                ! v*uy = v*D-y(u)
                                agu(vc,uc)= au*(3.*u(i1,i2,i3,uc)-4.*u(
     & i1,i2-1,i3,uc)+u(i1,i2-2,i3,uc))/(2.*dx(1))
                                ! v*vy = v*D-y(v)
                                agu(vc,vc)= au*(3.*u(i1,i2,i3,vc)-4.*u(
     & i1,i2-1,i3,vc)+u(i1,i2-2,i3,vc))/(2.*dx(1))
                             else
                                ! v*uy = v*D+y(u)
                                agu(vc,uc)= au*(-u(i1,i2+2,i3,uc)+4.*u(
     & i1,i2+1,i3,uc)-3.*u(i1,i2,i3,uc))/(2.*dx(1))
                                ! v*vy = v*D+y(v) 
                                agu(vc,vc)= au*(-u(i1,i2+2,i3,vc)+4.*u(
     & i1,i2+1,i3,vc)-3.*u(i1,i2,i3,vc))/(2.*dx(1))
                             end if
                          elseif( upwindOrder.eq.3 )then
                             au = uu(i1,i2,i3,uc)
                             if( au.gt.0. )then
                                ! u*ux = u*D-x(u)
                                agu(uc,uc)= au*(2.*u(i1+1,i2,i3,uc)+3.*
     & u(i1,i2,i3,uc)-6.*u(i1-1,i2,i3,uc)+u(i1-2,i2,i3,uc))/(6.*dx(0))
                                ! u*vx = u*D-x(v)
                                agu(uc,vc)= au*(2.*u(i1+1,i2,i3,vc)+3.*
     & u(i1,i2,i3,vc)-6.*u(i1-1,i2,i3,vc)+u(i1-2,i2,i3,vc))/(6.*dx(0))
                             else
                                ! u*ux = u*D+x(u)
                                agu(uc,uc)= au*(-u(i1+2,i2,i3,uc)+6.*u(
     & i1+1,i2,i3,uc)-3.*u(i1,i2,i3,uc)-2.*u(i1-1,i2,i3,uc))/(6.*dx(0)
     & )
                                ! u*vx = u*D+x(v)
                                agu(uc,vc)= au*(-u(i1+2,i2,i3,vc)+6.*u(
     & i1+1,i2,i3,vc)-3.*u(i1,i2,i3,vc)-2.*u(i1-1,i2,i3,vc))/(6.*dx(0)
     & )
                             end if
                             au = uu(i1,i2,i3,vc)
                             if( au.gt.0. )then
                                ! v*uy = v*D-y(u)
                                agu(vc,uc)= au*(2.*u(i1,i2+1,i3,uc)+3.*
     & u(i1,i2,i3,uc)-6.*u(i1,i2-1,i3,uc)+u(i1,i2-2,i3,uc))/(6.*dx(1))
                                ! v*vy = v*D-y(v)
                                agu(vc,vc)= au*(2.*u(i1,i2+1,i3,vc)+3.*
     & u(i1,i2,i3,vc)-6.*u(i1,i2-1,i3,vc)+u(i1,i2-2,i3,vc))/(6.*dx(1))
                             else
                                ! v*uy = v*D+y(u)
                                agu(vc,uc)= au*(-u(i1,i2+2,i3,uc)+6.*u(
     & i1,i2+1,i3,uc)-3.*u(i1,i2,i3,uc)-2.*u(i1,i2-1,i3,uc))/(6.*dx(1)
     & )
                                ! v*vy = v*D+y(v) 
                                agu(vc,vc)= au*(-u(i1,i2+2,i3,vc)+6.*u(
     & i1,i2+1,i3,vc)-3.*u(i1,i2,i3,vc)-2.*u(i1,i2-1,i3,vc))/(6.*dx(1)
     & )
                             end if
                          end if
                      ! #If ("UPWIND" == "CENTERED") 
                      !  ! -- centered advection ---
                      !  ! write(*,'(" getAdvection -- centered")')
                      !  #If "2" == "2"
                      !   ! -- 2D --
                      !   agu(uc,uc)=UU(uc)*UX(uc)
                      !   agu(vc,uc)=UU(vc)*UY(uc)
                      !   agu(uc,vc)=UU(uc)*UX(vc)
                      !   agu(vc,vc)=UU(vc)*UY(vc)
                      !  #Else
                      !   ! -- 3D -- *check me*
                      !   agu(uc,uc)=UU(uc)*UX(uc)
                      !   agu(vc,uc)=UU(vc)*UY(uc)
                      !   agu(wc,uc)=UU(wc)*UY(uc)
                      !   agu(uc,vc)=UU(uc)*UX(vc)
                      !   agu(vc,vc)=UU(vc)*UY(vc)
                      !   agu(wc,vc)=UU(wc)*UY(vc)
                      !   agu(uc,wc)=UU(uc)*UX(wc)
                      !   agu(vc,wc)=UU(vc)*UY(wc)
                      !   agu(wc,wc)=UU(wc)*UY(wc)
                      !  #End
                      ! #Elif "UPWIND" == "UPWIND" 
                      !   ! --- upwind scheme ---
                      !   ! for testing output this next message:
                      !   if( t.le. 0. )then
                      !     write(*,'(" getAdvection upwind scheme (7)")') 
                      !   end if
                      !   getUpwindAdvection(u,i1,i2,i3,NONE,2,4,rectangular, agu)
                      ! #Elif "UPWIND" == "BWENO" 
                      !   ! --- Bweno scheme ---
                      !   ! for testing output this message:
                      !   if( t.le. 0. )then
                      !      write(*,'(" getAdvection BWENO scheme (7)")') 
                      !   end if
                      !   getBwenoAdvection(u,i1,i2,i3,NONE,2,4,rectangular, agu)
                      ! #Else
                      !   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
                      !   stop 999
                      ! #End
                        ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-
     & ux42r(i1,i2,i3,pc)
                        ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-
     & uy42r(i1,i2,i3,pc)
                 end if
                end do
                end do
                end do
               else
                do i3=n3a,n3b
                do i2=n2a,n2b
                do i1=n1a,n1b
                   ! explicit terms only, no diffusion
                     ! upwind approximation/bweno
                       ! --- upwind scheme ---
                       ! for testing output this next message:
                       ! if( t.le. 0. )then
                       !    write(*,'(" getAdvection upwind scheme (7)")') 
                       ! end if
                         ! --- CARTESIAN GRID ---
                         !- agu(uc,uc)=UU(uc)*UX(uc)
                         !- agu(vc,uc)=UU(vc)*UY(uc)
                         !- agu(uc,vc)=UU(uc)*UX(vc)
                         !- agu(vc,vc)=UU(vc)*UY(vc)
                         ! -- first order upwind --
                         if( upwindOrder.eq.1 )then
                          au = uu(i1,i2,i3,uc)
                          if( au.gt.0. )then
                            ! u*ux = u*D-x(u)
                            agu(uc,uc)= au*(u(i1,i2,i3,uc)-u(i1-1,i2,
     & i3,uc))/(dx(0))
                            ! u*vx = u*D-x(v)
                            agu(uc,vc)= au*(u(i1,i2,i3,vc)-u(i1-1,i2,
     & i3,vc))/(dx(0))
                          else
                            ! u*ux = u*D+x(u)
                            agu(uc,uc)= au*(u(i1+1,i2,i3,uc)-u(i1,i2,
     & i3,uc))/(dx(0))
                            ! u*vx = u*D+x(v)
                            agu(uc,vc)= au*(u(i1+1,i2,i3,vc)-u(i1,i2,
     & i3,vc))/(dx(0))
                          end if
                          au = uu(i1,i2,i3,vc)
                          if( au.gt.0. )then
                            ! v*uy = v*D-y(u)
                            agu(vc,uc)= au*(u(i1,i2,i3,uc)-u(i1,i2-1,
     & i3,uc))/(dx(1))
                            ! v*vy = v*D-y(v)
                            agu(vc,vc)= au*(u(i1,i2,i3,vc)-u(i1,i2-1,
     & i3,vc))/(dx(1))
                          else
                            ! v*uy = v*D+y(u)
                            agu(vc,uc)= au*(u(i1,i2+1,i3,uc)-u(i1,i2,
     & i3,uc))/(dx(1))
                            ! v*vy = v*D+y(v) 
                            agu(vc,vc)= au*(u(i1,i2+1,i3,vc)-u(i1,i2,
     & i3,vc))/(dx(1))
                          end if
                         elseif( upwindOrder.eq.2 )then
                           ! write(*,'(" finish me, upwindOrder=",i2)') upwindOrder
                           ! stop 222
                            au = uu(i1,i2,i3,uc)
                            if( au.gt.0. )then
                               ! u*ux = u*D-x(u)
                               agu(uc,uc)= au*(3.*u(i1,i2,i3,uc)-4.*u(
     & i1-1,i2,i3,uc)+u(i1-2,i2,i3,uc))/(2.*dx(0))
                               ! u*vx = u*D-x(v)
                               agu(uc,vc)= au*(3.*u(i1,i2,i3,vc)-4.*u(
     & i1-1,i2,i3,vc)+u(i1-2,i2,i3,vc))/(2.*dx(0))
                            else
                               ! u*ux = u*D+x(u)
                               agu(uc,uc)= au*(-u(i1+2,i2,i3,uc)+4.*u(
     & i1+1,i2,i3,uc)-3.*u(i1,i2,i3,uc))/(2.*dx(0))
                               ! u*vx = u*D+x(v)
                               agu(uc,vc)= au*(-u(i1+2,i2,i3,vc)+4.*u(
     & i1+1,i2,i3,vc)-3.*u(i1,i2,i3,vc))/(2.*dx(0))
                            end if
                            au = uu(i1,i2,i3,vc)
                            if( au.gt.0. )then
                               ! v*uy = v*D-y(u)
                               agu(vc,uc)= au*(3.*u(i1,i2,i3,uc)-4.*u(
     & i1,i2-1,i3,uc)+u(i1,i2-2,i3,uc))/(2.*dx(1))
                               ! v*vy = v*D-y(v)
                               agu(vc,vc)= au*(3.*u(i1,i2,i3,vc)-4.*u(
     & i1,i2-1,i3,vc)+u(i1,i2-2,i3,vc))/(2.*dx(1))
                            else
                               ! v*uy = v*D+y(u)
                               agu(vc,uc)= au*(-u(i1,i2+2,i3,uc)+4.*u(
     & i1,i2+1,i3,uc)-3.*u(i1,i2,i3,uc))/(2.*dx(1))
                               ! v*vy = v*D+y(v) 
                               agu(vc,vc)= au*(-u(i1,i2+2,i3,vc)+4.*u(
     & i1,i2+1,i3,vc)-3.*u(i1,i2,i3,vc))/(2.*dx(1))
                            end if
                         elseif( upwindOrder.eq.3 )then
                            au = uu(i1,i2,i3,uc)
                            if( au.gt.0. )then
                               ! u*ux = u*D-x(u)
                               agu(uc,uc)= au*(2.*u(i1+1,i2,i3,uc)+3.*
     & u(i1,i2,i3,uc)-6.*u(i1-1,i2,i3,uc)+u(i1-2,i2,i3,uc))/(6.*dx(0))
                               ! u*vx = u*D-x(v)
                               agu(uc,vc)= au*(2.*u(i1+1,i2,i3,vc)+3.*
     & u(i1,i2,i3,vc)-6.*u(i1-1,i2,i3,vc)+u(i1-2,i2,i3,vc))/(6.*dx(0))
                            else
                               ! u*ux = u*D+x(u)
                               agu(uc,uc)= au*(-u(i1+2,i2,i3,uc)+6.*u(
     & i1+1,i2,i3,uc)-3.*u(i1,i2,i3,uc)-2.*u(i1-1,i2,i3,uc))/(6.*dx(0)
     & )
                               ! u*vx = u*D+x(v)
                               agu(uc,vc)= au*(-u(i1+2,i2,i3,vc)+6.*u(
     & i1+1,i2,i3,vc)-3.*u(i1,i2,i3,vc)-2.*u(i1-1,i2,i3,vc))/(6.*dx(0)
     & )
                            end if
                            au = uu(i1,i2,i3,vc)
                            if( au.gt.0. )then
                               ! v*uy = v*D-y(u)
                               agu(vc,uc)= au*(2.*u(i1,i2+1,i3,uc)+3.*
     & u(i1,i2,i3,uc)-6.*u(i1,i2-1,i3,uc)+u(i1,i2-2,i3,uc))/(6.*dx(1))
                               ! v*vy = v*D-y(v)
                               agu(vc,vc)= au*(2.*u(i1,i2+1,i3,vc)+3.*
     & u(i1,i2,i3,vc)-6.*u(i1,i2-1,i3,vc)+u(i1,i2-2,i3,vc))/(6.*dx(1))
                            else
                               ! v*uy = v*D+y(u)
                               agu(vc,uc)= au*(-u(i1,i2+2,i3,uc)+6.*u(
     & i1,i2+1,i3,uc)-3.*u(i1,i2,i3,uc)-2.*u(i1,i2-1,i3,uc))/(6.*dx(1)
     & )
                               ! v*vy = v*D+y(v) 
                               agu(vc,vc)= au*(-u(i1,i2+2,i3,vc)+6.*u(
     & i1,i2+1,i3,vc)-3.*u(i1,i2,i3,vc)-2.*u(i1,i2-1,i3,vc))/(6.*dx(1)
     & )
                            end if
                         end if
                     ! #If ("UPWIND" == "CENTERED") 
                     !  ! -- centered advection ---
                     !  ! write(*,'(" getAdvection -- centered")')
                     !  #If "2" == "2"
                     !   ! -- 2D --
                     !   agu(uc,uc)=UU(uc)*UX(uc)
                     !   agu(vc,uc)=UU(vc)*UY(uc)
                     !   agu(uc,vc)=UU(uc)*UX(vc)
                     !   agu(vc,vc)=UU(vc)*UY(vc)
                     !  #Else
                     !   ! -- 3D -- *check me*
                     !   agu(uc,uc)=UU(uc)*UX(uc)
                     !   agu(vc,uc)=UU(vc)*UY(uc)
                     !   agu(wc,uc)=UU(wc)*UY(uc)
                     !   agu(uc,vc)=UU(uc)*UX(vc)
                     !   agu(vc,vc)=UU(vc)*UY(vc)
                     !   agu(wc,vc)=UU(wc)*UY(vc)
                     !   agu(uc,wc)=UU(uc)*UX(wc)
                     !   agu(vc,wc)=UU(vc)*UY(wc)
                     !   agu(wc,wc)=UU(wc)*UY(wc)
                     !  #End
                     ! #Elif "UPWIND" == "UPWIND" 
                     !   ! --- upwind scheme ---
                     !   ! for testing output this next message:
                     !   if( t.le. 0. )then
                     !     write(*,'(" getAdvection upwind scheme (7)")') 
                     !   end if
                     !   getUpwindAdvection(u,i1,i2,i3,NONE,2,4,rectangular, agu)
                     ! #Elif "UPWIND" == "BWENO" 
                     !   ! --- Bweno scheme ---
                     !   ! for testing output this message:
                     !   if( t.le. 0. )then
                     !      write(*,'(" getAdvection BWENO scheme (7)")') 
                     !   end if
                     !   getBwenoAdvection(u,i1,i2,i3,NONE,2,4,rectangular, agu)
                     ! #Else
                     !   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
                     !   stop 999
                     ! #End
                       ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-ux42r(
     & i1,i2,i3,pc)
                       ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-uy42r(
     & i1,i2,i3,pc)
                end do
                end do
                end do
               end if
             else
              write(*,*)'insdt: Unknown implicitOption=',implicitOption
              stop 6
             end if  ! end implicitOption
            else if( advectionOption.eq.bwenoAdvection )then
              ! --- bweno ---
             if( implicitOption .eq.computeImplicitTermsSeparately )
     & then
               if( useWhereMask.ne.0 )then
                do i3=n3a,n3b
                do i2=n2a,n2b
                do i1=n1a,n1b
                 if( mask(i1,i2,i3).gt.0 )then
                    ! explicit terms only, no diffusion
                      ! upwind approximation/bweno
                        ! --- Bweno scheme ---
                        ! for testing output this message:
                        ! if( t.le. 0. )then
                        !    write(*,'(" getAdvection BWENO scheme (7)")') 
                        ! end if
                           gvU = u(i1,i2,i3,uc) - uu(i1,i2,i3,uc)
                           gvV = u(i1,i2,i3,vc) - uu(i1,i2,i3,vc)
                           gvW = u(i1,i2,i3,wc) - uu(i1,i2,i3,wc)
                           orderOneFix = 1
                           augmentPlot = 1
                           if( upwindOrder.eq.4 )then
                          ! -- fourth order bweno --
                              au = uu(i1,i2,i3,uc)
                              aum = u(i1-1,i2,i3,uc) - gvU
                              aup = u(i1+1,i2,i3,uc) - gvU
                              !usage
                              !s1   is the shifted grid in r1 direction
                              !var  is the variable taking derivative (uc,vc,wc)
                              !vard is the variable of direction
                              !agu is only passed in for augmented plot
                              !expAd 
                              ! u*ux
                              s1=1
                              s2=0
                              s3=0
                              var = uc
                              vard= uc
                              drl = dx(0)
                                Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*
     & u(i1   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
                                Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    
     & u(i1-s1,i2-s2,i3-s3,var)
                                Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) 
     & - 2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
                                Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) 
     & -    u(i1   ,i2   ,i3   ,var)
                                Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,
     & i2-s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
                                Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var)
                                Amr   = Apl
                                Bmr   = Bpl
                                betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 
     & 1./4.*Bpl**2
                                betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 
     & 1./4.*Bpr**2
                                betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 
     & 1./4.*Bml**2
                                betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 
     & 1./4.*Bmr**2
                                ep     = 1e-15
                                ap1    = (1./2.)/(ep + betapl)**(4)
                                ap2    = (1./2.)/(ep + betapr)**(4)
                                am1    = (1./2.)/(ep + betaml)**(4)
                                am2    = (1./2.)/(ep + betamr)**(4)
                                !calculate weights
                                wp1 = ap1/(ap1+ap2)
                                wp2 = ap2/(ap1+ap2)
                                wm1 = am1/(am1+am2)
                                wm2 = am2/(am1+am2)
                                !mapping of the weights
                                bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                                bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                                bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                                bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                                wplw = bp1/(bp1+bp2)
                                wprw = bp2/(bp1+bp2)
                                wmlw = bm1/(bm1+bm2)
                                wmrw = bm2/(bm1+bm2)
                                !get face value of u
                                up = au
                                um = au
                                !picking side for the weights
                                if (up.ge.0) then
                                   wpl = max(wplw,wprw)
                                   wpr = min(wplw,wprw)
                                else
                                   wpl = min(wplw,wprw)
                                   wpr = max(wplw,wprw)
                                end if
                                if (um.ge.0) then
                                   wml = max(wmlw,wmrw)
                                   wmr = min(wmlw,wmrw)
                                else
                                   wml = min(wmlw,wmrw)
                                   wmr = max(wmlw,wmrw)
                                end if
                                expAd = 0.
                                if (orderOneFix.eq.1) then
                                   ! if (augmentPlot.eq.1) then
                                   !    agu(uc,13)=0
                                   ! end if
                                   if (aum*aup .le. 0) then
                                      aumax = max(abs(au),abs(aup),abs(
     & aum))
                                      expAd = (1./12.)*aumax*( u(i1-2*
     & s1,i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(
     & i1,i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*
     & s2,i3+2*s3,var))
                                      wpl = 0.5
                                      wpr = 0.5
                                      wml = 0.5
                                      wmr = 0.5
                                      ! plot au at where fix is turned on
                                      ! if (augmentPlot.eq.1) then
                                      !    agu(uc,13)=au
                                      ! end if
                                   end if
                                end if
                                ! wpl = 0.5
                                ! wpr = 0.5
                                ! wml = 0.5
                                ! wmr = 0.5
                                Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,
     & var) + 5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
                                Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 
     & 5.*u(i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
                                Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*
     & s3,var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
                                Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,
     & var) + 5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
                                Fp  = wpl*Fpl + wpr*Fpr
                                Fm  = wml*Fml + wmr*Fmr
                                   agu(vard,var)= au*(Fp - Fm)/drl + 
     & expAd/drl
                              !     if ( i1 .eq. 14 .and. i2.eq.16 ) then
                                      !     if ( i1 .eq. 14) then
                                      ! write(*,'(" betapl,betapr,betaml,betamr = ",4F10.8)')  betapl,betapr,betaml,betamr
                                      ! write(*,'(" ap1,ap2,am1,am2 = ",4e10.8)')  ap1,ap2,am1,am2
                                      ! write(*,'(" wp1,wpr,wm1,wm2 = ",4F10.8)')  wp1,wp2,wm1,wm2
                                      ! write(*,'(" bp1,bp2,bm1,bm2 = ",4F10.8)')  bp1,bp2,bm1,bm2
                                      ! write(*,'(" wplw,wprw,wmlw,wmrw = ",4F10.8)')  wplw,wprw,wmlw,wmrw
                                      ! write(*,'(" wpl,wpr,wml,wmr = ",4F10.8)')  wpl,wpr,wml,wmr
                                      ! write(*,'(" agu = ",F10.8)') agu(vard,var)
                              !     end if
                              !u*vx
                              s1=1
                              s2=0
                              s3=0
                              var = vc
                              vard= uc
                              drl = dx(0)
                                Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*
     & u(i1   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
                                Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    
     & u(i1-s1,i2-s2,i3-s3,var)
                                Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) 
     & - 2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
                                Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) 
     & -    u(i1   ,i2   ,i3   ,var)
                                Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,
     & i2-s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
                                Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var)
                                Amr   = Apl
                                Bmr   = Bpl
                                betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 
     & 1./4.*Bpl**2
                                betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 
     & 1./4.*Bpr**2
                                betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 
     & 1./4.*Bml**2
                                betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 
     & 1./4.*Bmr**2
                                ep     = 1e-15
                                ap1    = (1./2.)/(ep + betapl)**(4)
                                ap2    = (1./2.)/(ep + betapr)**(4)
                                am1    = (1./2.)/(ep + betaml)**(4)
                                am2    = (1./2.)/(ep + betamr)**(4)
                                !calculate weights
                                wp1 = ap1/(ap1+ap2)
                                wp2 = ap2/(ap1+ap2)
                                wm1 = am1/(am1+am2)
                                wm2 = am2/(am1+am2)
                                !mapping of the weights
                                bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                                bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                                bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                                bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                                wplw = bp1/(bp1+bp2)
                                wprw = bp2/(bp1+bp2)
                                wmlw = bm1/(bm1+bm2)
                                wmrw = bm2/(bm1+bm2)
                                !get face value of u
                                up = au
                                um = au
                                !picking side for the weights
                                if (up.ge.0) then
                                   wpl = max(wplw,wprw)
                                   wpr = min(wplw,wprw)
                                else
                                   wpl = min(wplw,wprw)
                                   wpr = max(wplw,wprw)
                                end if
                                if (um.ge.0) then
                                   wml = max(wmlw,wmrw)
                                   wmr = min(wmlw,wmrw)
                                else
                                   wml = min(wmlw,wmrw)
                                   wmr = max(wmlw,wmrw)
                                end if
                                expAd = 0.
                                if (orderOneFix.eq.1) then
                                   ! if (augmentPlot.eq.1) then
                                   !    agu(uc,13)=0
                                   ! end if
                                   if (aum*aup .le. 0) then
                                      aumax = max(abs(au),abs(aup),abs(
     & aum))
                                      expAd = (1./12.)*aumax*( u(i1-2*
     & s1,i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(
     & i1,i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*
     & s2,i3+2*s3,var))
                                      wpl = 0.5
                                      wpr = 0.5
                                      wml = 0.5
                                      wmr = 0.5
                                      ! plot au at where fix is turned on
                                      ! if (augmentPlot.eq.1) then
                                      !    agu(uc,13)=au
                                      ! end if
                                   end if
                                end if
                                ! wpl = 0.5
                                ! wpr = 0.5
                                ! wml = 0.5
                                ! wmr = 0.5
                                Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,
     & var) + 5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
                                Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 
     & 5.*u(i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
                                Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*
     & s3,var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
                                Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,
     & var) + 5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
                                Fp  = wpl*Fpl + wpr*Fpr
                                Fm  = wml*Fml + wmr*Fmr
                                   agu(vard,var)= au*(Fp - Fm)/drl + 
     & expAd/drl
                              !     if ( i1 .eq. 14 .and. i2.eq.16 ) then
                                      !     if ( i1 .eq. 14) then
                                      ! write(*,'(" betapl,betapr,betaml,betamr = ",4F10.8)')  betapl,betapr,betaml,betamr
                                      ! write(*,'(" ap1,ap2,am1,am2 = ",4e10.8)')  ap1,ap2,am1,am2
                                      ! write(*,'(" wp1,wpr,wm1,wm2 = ",4F10.8)')  wp1,wp2,wm1,wm2
                                      ! write(*,'(" bp1,bp2,bm1,bm2 = ",4F10.8)')  bp1,bp2,bm1,bm2
                                      ! write(*,'(" wplw,wprw,wmlw,wmrw = ",4F10.8)')  wplw,wprw,wmlw,wmrw
                                      ! write(*,'(" wpl,wpr,wml,wmr = ",4F10.8)')  wpl,wpr,wml,wmr
                                      ! write(*,'(" agu = ",F10.8)') agu(vard,var)
                              !     end if
                              !v*uy
                              au = uu(i1,i2,i3,vc)
                              aum = u(i1,i2-1,i3,vc) - gvV
                              aup = u(i1,i2+1,i3,vc) - gvV
                              s1=0
                              s2=1
                              s3=0
                              var = uc
                              vard= vc
                              drl = dx(1)
                                Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*
     & u(i1   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
                                Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    
     & u(i1-s1,i2-s2,i3-s3,var)
                                Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) 
     & - 2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
                                Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) 
     & -    u(i1   ,i2   ,i3   ,var)
                                Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,
     & i2-s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
                                Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var)
                                Amr   = Apl
                                Bmr   = Bpl
                                betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 
     & 1./4.*Bpl**2
                                betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 
     & 1./4.*Bpr**2
                                betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 
     & 1./4.*Bml**2
                                betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 
     & 1./4.*Bmr**2
                                ep     = 1e-15
                                ap1    = (1./2.)/(ep + betapl)**(4)
                                ap2    = (1./2.)/(ep + betapr)**(4)
                                am1    = (1./2.)/(ep + betaml)**(4)
                                am2    = (1./2.)/(ep + betamr)**(4)
                                !calculate weights
                                wp1 = ap1/(ap1+ap2)
                                wp2 = ap2/(ap1+ap2)
                                wm1 = am1/(am1+am2)
                                wm2 = am2/(am1+am2)
                                !mapping of the weights
                                bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                                bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                                bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                                bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                                wplw = bp1/(bp1+bp2)
                                wprw = bp2/(bp1+bp2)
                                wmlw = bm1/(bm1+bm2)
                                wmrw = bm2/(bm1+bm2)
                                !get face value of u
                                up = au
                                um = au
                                !picking side for the weights
                                if (up.ge.0) then
                                   wpl = max(wplw,wprw)
                                   wpr = min(wplw,wprw)
                                else
                                   wpl = min(wplw,wprw)
                                   wpr = max(wplw,wprw)
                                end if
                                if (um.ge.0) then
                                   wml = max(wmlw,wmrw)
                                   wmr = min(wmlw,wmrw)
                                else
                                   wml = min(wmlw,wmrw)
                                   wmr = max(wmlw,wmrw)
                                end if
                                expAd = 0.
                                if (orderOneFix.eq.1) then
                                   ! if (augmentPlot.eq.1) then
                                   !    agu(uc,13)=0
                                   ! end if
                                   if (aum*aup .le. 0) then
                                      aumax = max(abs(au),abs(aup),abs(
     & aum))
                                      expAd = (1./12.)*aumax*( u(i1-2*
     & s1,i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(
     & i1,i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*
     & s2,i3+2*s3,var))
                                      wpl = 0.5
                                      wpr = 0.5
                                      wml = 0.5
                                      wmr = 0.5
                                      ! plot au at where fix is turned on
                                      ! if (augmentPlot.eq.1) then
                                      !    agu(uc,13)=au
                                      ! end if
                                   end if
                                end if
                                ! wpl = 0.5
                                ! wpr = 0.5
                                ! wml = 0.5
                                ! wmr = 0.5
                                Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,
     & var) + 5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
                                Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 
     & 5.*u(i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
                                Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*
     & s3,var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
                                Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,
     & var) + 5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
                                Fp  = wpl*Fpl + wpr*Fpr
                                Fm  = wml*Fml + wmr*Fmr
                                   agu(vard,var)= au*(Fp - Fm)/drl + 
     & expAd/drl
                              !     if ( i1 .eq. 14 .and. i2.eq.16 ) then
                                      !     if ( i1 .eq. 14) then
                                      ! write(*,'(" betapl,betapr,betaml,betamr = ",4F10.8)')  betapl,betapr,betaml,betamr
                                      ! write(*,'(" ap1,ap2,am1,am2 = ",4e10.8)')  ap1,ap2,am1,am2
                                      ! write(*,'(" wp1,wpr,wm1,wm2 = ",4F10.8)')  wp1,wp2,wm1,wm2
                                      ! write(*,'(" bp1,bp2,bm1,bm2 = ",4F10.8)')  bp1,bp2,bm1,bm2
                                      ! write(*,'(" wplw,wprw,wmlw,wmrw = ",4F10.8)')  wplw,wprw,wmlw,wmrw
                                      ! write(*,'(" wpl,wpr,wml,wmr = ",4F10.8)')  wpl,wpr,wml,wmr
                                      ! write(*,'(" agu = ",F10.8)') agu(vard,var)
                              !     end if
                              !v*vy
                              s1=0
                              s2=1
                              s3=0
                              var = vc
                              vard= vc
                              drl = dx(1)
                                Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*
     & u(i1   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
                                Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    
     & u(i1-s1,i2-s2,i3-s3,var)
                                Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) 
     & - 2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
                                Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) 
     & -    u(i1   ,i2   ,i3   ,var)
                                Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,
     & i2-s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
                                Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var)
                                Amr   = Apl
                                Bmr   = Bpl
                                betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 
     & 1./4.*Bpl**2
                                betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 
     & 1./4.*Bpr**2
                                betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 
     & 1./4.*Bml**2
                                betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 
     & 1./4.*Bmr**2
                                ep     = 1e-15
                                ap1    = (1./2.)/(ep + betapl)**(4)
                                ap2    = (1./2.)/(ep + betapr)**(4)
                                am1    = (1./2.)/(ep + betaml)**(4)
                                am2    = (1./2.)/(ep + betamr)**(4)
                                !calculate weights
                                wp1 = ap1/(ap1+ap2)
                                wp2 = ap2/(ap1+ap2)
                                wm1 = am1/(am1+am2)
                                wm2 = am2/(am1+am2)
                                !mapping of the weights
                                bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                                bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                                bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                                bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                                wplw = bp1/(bp1+bp2)
                                wprw = bp2/(bp1+bp2)
                                wmlw = bm1/(bm1+bm2)
                                wmrw = bm2/(bm1+bm2)
                                !get face value of u
                                up = au
                                um = au
                                !picking side for the weights
                                if (up.ge.0) then
                                   wpl = max(wplw,wprw)
                                   wpr = min(wplw,wprw)
                                else
                                   wpl = min(wplw,wprw)
                                   wpr = max(wplw,wprw)
                                end if
                                if (um.ge.0) then
                                   wml = max(wmlw,wmrw)
                                   wmr = min(wmlw,wmrw)
                                else
                                   wml = min(wmlw,wmrw)
                                   wmr = max(wmlw,wmrw)
                                end if
                                expAd = 0.
                                if (orderOneFix.eq.1) then
                                   ! if (augmentPlot.eq.1) then
                                   !    agu(uc,13)=0
                                   ! end if
                                   if (aum*aup .le. 0) then
                                      aumax = max(abs(au),abs(aup),abs(
     & aum))
                                      expAd = (1./12.)*aumax*( u(i1-2*
     & s1,i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(
     & i1,i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*
     & s2,i3+2*s3,var))
                                      wpl = 0.5
                                      wpr = 0.5
                                      wml = 0.5
                                      wmr = 0.5
                                      ! plot au at where fix is turned on
                                      ! if (augmentPlot.eq.1) then
                                      !    agu(uc,13)=au
                                      ! end if
                                   end if
                                end if
                                ! wpl = 0.5
                                ! wpr = 0.5
                                ! wml = 0.5
                                ! wmr = 0.5
                                Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,
     & var) + 5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
                                Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 
     & 5.*u(i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
                                Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*
     & s3,var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
                                Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,
     & var) + 5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
                                Fp  = wpl*Fpl + wpr*Fpr
                                Fm  = wml*Fml + wmr*Fmr
                                   agu(vard,var)= au*(Fp - Fm)/drl + 
     & expAd/drl
                              !     if ( i1 .eq. 14 .and. i2.eq.16 ) then
                                      !     if ( i1 .eq. 14) then
                                      ! write(*,'(" betapl,betapr,betaml,betamr = ",4F10.8)')  betapl,betapr,betaml,betamr
                                      ! write(*,'(" ap1,ap2,am1,am2 = ",4e10.8)')  ap1,ap2,am1,am2
                                      ! write(*,'(" wp1,wpr,wm1,wm2 = ",4F10.8)')  wp1,wp2,wm1,wm2
                                      ! write(*,'(" bp1,bp2,bm1,bm2 = ",4F10.8)')  bp1,bp2,bm1,bm2
                                      ! write(*,'(" wplw,wprw,wmlw,wmrw = ",4F10.8)')  wplw,wprw,wmlw,wmrw
                                      ! write(*,'(" wpl,wpr,wml,wmr = ",4F10.8)')  wpl,wpr,wml,wmr
                                      ! write(*,'(" agu = ",F10.8)') agu(vard,var)
                              !     end if
                             ! u*ux-------------------------------------------------------------------------
                             !calculate smooth factor
                             ! Apl   = u(i1+1,i2,i3,uc) - 2.*u(i1  ,i2,i3,uc)   + u(i1-1,i2,i3,uc)
                             ! Bpl   = u(i1+1,i2,i3,uc) -    u(i1-1,i2,i3,uc)
                             ! Apr   = u(i1+2,i2,i3,uc) - 2.*u(i1+1,i2,i3,uc) + u(i1,i2,i3,uc)
                             ! Bpr   = u(i1+2,i2,i3,uc) -    u(i1  ,i2,i3,uc)
                             ! Aml   = u(i1,i2,i3,uc) - 2.*u(i1-1,i2,i3,uc) + u(i1-2,i2,i3,uc)
                             ! Bml   = u(i1,i2,i3,uc) -    u(i1-2,i2,i3,uc)
                             ! Amr   = Apl
                             ! Bmr   = Bpl
                             ! betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 1./4.*Bpl**2
                             ! betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 1./4.*Bpr**2
                             ! betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 1./4.*Bml**2
                             ! betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 1./4.*Bmr**2
                             ! ep     = 1e-15
                             ! ap1    = (1./2.)/(ep + betapl)**(4)
                             ! ap2    = (1./2.)/(ep + betapr)**(4)
                             ! am1    = (1./2.)/(ep + betaml)**(4)
                             ! am2    = (1./2.)/(ep + betamr)**(4)
                             ! !calculate weights
                             ! wp1 = ap1/(ap1+ap2)
                             ! wp2 = ap2/(ap1+ap2)
                             ! wm1 = am1/(am1+am2)
                             ! wm2 = am2/(am1+am2)
                             ! !mapping of the weights
                             ! bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                             ! bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                             ! bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                             ! bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                             ! wplw = bp1/(bp1+bp2)
                             ! wprw = bp2/(bp1+bp2)
                             ! wmlw = bm1/(bm1+bm2)
                             ! wmrw = bm2/(bm1+bm2)
                             ! !get face value of u
                             ! up = au
                             ! um = au
                             ! !picking side for the weights
                             !  if (up.ge.0) then
                             !    wpl = max(wplw,wprw)
                             !    wpr = min(wplw,wprw)
                             ! else
                             !    wpl = min(wplw,wprw)
                             !    wpr = max(wplw,wprw)
                             ! end if
                             ! if (um.ge.0) then
                             !    wml = max(wmlw,wmrw)
                             !    wmr = min(wmlw,wmrw)
                             ! else
                             !    wml = min(wmlw,wmrw)
                             !    wmr = max(wmlw,wmrw)
                             ! end if
                             ! expAd = 0.
                             ! if (orderOneFix.eq.1) then 
                             !    if (augmentPlot.eq.1) then
                             !       agu(uc,13)=0
                             !    end if
                             !    if (aum*aup .le. 0) then
                             !       aumax = max(abs(au),abs(aup),abs(aum))
                             !       expAd = (1./12.)*aumax*( u(i1-2,i2,i3,uc) - 4.*u(i1-1,i2,i3,uc) + 6.*u(i1,i2,i3,uc) - 4.*u(i1+1,i2,i3,uc) +  u(i1+2,i2,i3,uc))
                             !       wpl = 0.5
                             !       wpr = 0.5
                             !       wml = 0.5
                             !       wmr = 0.5
                             !       ! plot au at where fix is turned on
                             !       if (augmentPlot.eq.1) then
                             !          agu(uc,13)=au
                             !       end if
                             !    end if
                             ! end if
                             ! wpl = 0.5
                             ! wpr = 0.5
                             ! wml = 0.5
                             ! wmr = 0.5
                             ! Fpl = 1./6.*(  -u(i1-1,i2,i3,uc) + 5.*u(i1,i2,i3,uc)   + 2.*u(i1+1,i2,i3,uc))
                             ! Fpr = 1./6.*( 2.*u(i1,i2,i3,uc)   + 5.*u(i1+1,i2,i3,uc) -   u(i1+2,i2,i3,uc))
                             ! Fml = 1./6.*(  -u(i1-2,i2,i3,uc) + 5.*u(i1-1,i2,i3,uc)   + 2.*u(i1,i2,i3,uc))
                             ! Fmr = 1./6.*( 2.*u(i1-1,i2,i3,uc) + 5.*u(i1,i2,i3,uc) -       u(i1+1,i2,i3,uc))
                             ! Fp  = wpl*Fpl + wpr*Fpr
                             ! Fm  = wml*Fml + wmr*Fmr
                             ! agu(uc,uc)= au*(Fp - Fm)/(dx(0)) + expAd*dx(0)**(-1)
                             ! if (augmentPlot.eq.1) then
                             !    agu(uc,3)=au
                             !    agu(uc,4)=wpl
                             !    agu(uc,5)=wpr
                             !    agu(uc,6)=wml
                             !    agu(uc,7)=wmr
                             ! end if
                             ! u*vx-------------------------------------------------------------------------
                             !calculate smooth factor
                             ! Apl   = u(i1+1,i2,i3,vc) - 2.*u(i1,i2,i3,vc)   + u(i1-1,i2,i3,vc)
                             ! Bpl   = u(i1+1,i2,i3,vc) -    u(i1-1,i2,i3,vc)
                             ! Apr   = u(i1+2,i2,i3,vc) - 2.*u(i1+1,i2,i3,vc) + u(i1,i2,i3,vc)
                             ! Bpr   = u(i1+2,i2,i3,vc) -    u(i1,i2,i3,vc)
                             ! Aml   = u(i1,i2,i3,vc) - 2.*u(i1-1,i2,i3,vc) + u(i1-2,i2,i3,vc)
                             ! Bml   = u(i1,i2,i3,vc) -    u(i1-2,i2,i3,vc)
                             ! Amr   = Apl
                             ! Bmr   = Bpl
                             ! betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 1./4.*Bpl**2
                             ! betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 1./4.*Bpr**2
                             ! betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 1./4.*Bml**2
                             ! betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 1./4.*Bmr**2
                             ! ep     = 1e-15
                             ! ap1    = (1./2.)/(ep + betapl)**(4)
                             ! ap2    = (1./2.)/(ep + betapr)**(4)
                             ! am1    = (1./2.)/(ep + betaml)**(4)
                             ! am2    = (1./2.)/(ep + betamr)**(4)
                             ! !calculate weights
                             ! wp1 = ap1/(ap1+ap2)
                             ! wp2 = ap2/(ap1+ap2)
                             ! wm1 = am1/(am1+am2)
                             ! wm2 = am2/(am1+am2)
                             ! !mapping of the weights
                             ! bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                             ! bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                             ! bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                             ! bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                             ! wplw = bp1/(bp1+bp2)
                             ! wprw = bp2/(bp1+bp2)
                             ! wmlw = bm1/(bm1+bm2)
                             ! wmrw = bm2/(bm1+bm2)
                             ! !get face value of u
                             ! ! up = 0.5*(u(i1  ,i2,i3,uc)+u(i1+1,i2,i3,uc))
                             ! ! um = 0.5*(u(i1-1,i2,i3,uc)+u(i1  ,i2,i3,uc))
                             ! up = au
                             ! um = au
                             ! !picking side for the weights
                             ! if (up.ge.0) then
                             !    wpl = max(wplw,wprw)
                             !    wpr = min(wplw,wprw)
                             ! else
                             !    wpl = min(wplw,wprw)
                             !    wpr = max(wplw,wprw)
                             ! end if
                             ! if (um.ge.0) then
                             !    wml = max(wmlw,wmrw)
                             !    wmr = min(wmlw,wmrw)
                             ! else
                             !    wml = min(wmlw,wmrw)
                             !    wmr = max(wmlw,wmrw)
                             ! end if
                             ! expAd = 0.
                             ! if (orderOneFix.eq.1) then
                             !    if (aum*aup .le. 0) then
                             !       aumax = max(abs(au),abs(aup),abs(aum))
                             !       expAd = (1./12.)*aumax*( u(i1-2,i2,i3,vc) - 4.*u(i1-1,i2,i3,vc) + 6.*u(i1,i2,i3,vc) - 4.*u(i1+1,i2,i3,vc) +  u(i1+2,i2,i3,vc))
                             !       wpl = 0.5
                             !       wpr = 0.5
                             !       wml = 0.5
                             !       wmr = 0.5
                             !    end if
                             ! end if
                             ! Fpl = 1./6.*(  -u(i1-1,i2,i3,vc) + 5.*u(i1,i2,i3,vc)   + 2.*u(i1+1,i2,i3,vc))
                             ! Fpr = 1./6.*( 2.*u(i1,i2,i3,vc)   + 5.*u(i1+1,i2,i3,vc) -   u(i1+2,i2,i3,vc))
                             ! Fml = 1./6.*(  -u(i1-2,i2,i3,vc) + 5.*u(i1-1,i2,i3,vc)   + 2.*u(i1,i2,i3,vc))
                             ! Fmr = 1./6.*( 2.*u(i1-1,i2,i3,vc) + 5.*u(i1,i2,i3,vc) -       u(i1+1,i2,i3,vc))
                             ! Fp  = wpl*Fpl + wpr*Fpr
                             ! Fm  = wml*Fml + wmr*Fmr
                             ! agu(uc,vc)= au*(Fp - Fm)/(dx(0)) + expAd*dx(0)**(-1) 
                             ! if (augmentPlot .eq. 1) then 
                             !    agu(vc,4)=wpl
                             !    agu(vc,5)=wpr
                             !    agu(vc,6)=wml
                             !    agu(vc,7)=wmr
                             ! end if
                             !v*uy-------------------------------------------------------------------------
                             ! au = u(i1,i2,i3,vc)
                             ! aum = u(i1,i2-1,i3,vc)
                             ! aup = u(i1,i2+1,i3,vc)
                             !calculate smooth factor
                             ! Apl   = u(i1,i2+1,i3,uc) - 2.*u(i1,i2,i3,uc)   + u(i1,i2-1,i3,uc)
                             ! Bpl   = u(i1,i2+1,i3,uc) -    u(i1,i2-1,i3,uc)
                             ! Apr   = u(i1,i2+2,i3,uc) - 2.*u(i1,i2+1,i3,uc) + u(i1,i2,i3,uc)
                             ! Bpr   = u(i1,i2+2,i3,uc) -    u(i1,i2,i3,uc)
                             ! Aml   = u(i1,i2,i3,uc) - 2.*u(i1,i2-1,i3,uc) + u(i1,i2-2,i3,uc)
                             ! Bml   = u(i1,i2,i3,uc) -    u(i1,i2-2,i3,uc)
                             ! Amr   = Apl
                             ! Bmr   = Bpl
                             ! betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 1./4.*Bpl**2
                             ! betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 1./4.*Bpr**2
                             ! betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 1./4.*Bml**2
                             ! betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 1./4.*Bmr**2
                             ! ep     = 1e-15
                             ! ap1    = (1./2.)/(ep + betapl)**(4)
                             ! ap2    = (1./2.)/(ep + betapr)**(4)
                             ! am1    = (1./2.)/(ep + betaml)**(4)
                             ! am2    = (1./2.)/(ep + betamr)**(4)
                             ! !calculate weights
                             ! wp1 = ap1/(ap1+ap2)
                             ! wp2 = ap2/(ap1+ap2)
                             ! wm1 = am1/(am1+am2)
                             ! wm2 = am2/(am1+am2)
                             ! !mapping of the weights
                             ! bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                             ! bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                             ! bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                             ! bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                             ! wplw = bp1/(bp1+bp2)
                             ! wprw = bp2/(bp1+bp2)
                             ! wmlw = bm1/(bm1+bm2)
                             ! wmrw = bm2/(bm1+bm2)
                             ! !get face value of v
                             ! ! up = 0.5*(u(i1  ,i2+1,i3,vc)+u(i1,i2 ,i3,vc))
                             ! ! um = 0.5*(u(i1  ,i2  ,i3,vc)+u(i1,i2-1,i3,vc))
                             ! up = au
                             ! um = au
                             ! !picking side for the weights
                             !  if (up.ge.0) then
                             !    wpl = max(wplw,wprw)
                             !    wpr = min(wplw,wprw)
                             ! else
                             !    wpl = min(wplw,wprw)
                             !    wpr = max(wplw,wprw)
                             ! end if
                             ! if (um.ge.0) then
                             !    wml = max(wmlw,wmrw)
                             !    wmr = min(wmlw,wmrw)
                             ! else
                             !    wml = min(wmlw,wmrw)
                             !    wmr = max(wmlw,wmrw)
                             ! end if
                             ! expAd = 0.
                             ! if (orderOneFix.eq.1) then 
                             !    if (augmentPlot .eq. 1) then 
                             !       agu(vc,13)=0
                             !    end if
                             !    if (aum*aup .le. 0) then
                             !       aumax = max(abs(au),abs(aup),abs(aum))
                             !       expAd = (1./12.)*aumax*( u(i1,i2-2,i3,uc) - 4.*u(i1,i2-1,i3,uc) + 6.*u(i1,i2,i3,uc) - 4.*u(i1,i2+1,i3,uc) +  u(i1,i2+2,i3,uc))
                             !       wpl = 0.5
                             !       wpr = 0.5
                             !       wml = 0.5
                             !       wmr = 0.5
                             !       if (augmentPlot .eq. 1) then 
                             !          agu(vc,13)=au
                             !       end if
                             !    end if
                             ! end if
                             ! Fpl = 1./6.*(  -u(i1,i2-1,i3,uc) + 5.*u(i1,i2,i3,uc)   + 2.*u(i1,i2+1,i3,uc))
                             ! Fpr = 1./6.*( 2.*u(i1,i2,i3,uc)   + 5.*u(i1,i2+1,i3,uc) -   u(i1,i2+2,i3,uc))
                             ! Fml = 1./6.*(  -u(i1,i2-2,i3,uc) + 5.*u(i1,i2-1,i3,uc)   + 2.*u(i1,i2,i3,uc))
                             ! Fmr = 1./6.*( 2.*u(i1,i2-1,i3,uc) + 5.*u(i1,i2,i3,uc) -       u(i1,i2+1,i3,uc))
                             ! Fp  = wpl*Fpl + wpr*Fpr
                             ! Fm  = wml*Fml + wmr*Fmr
                             ! agu(vc,uc)= au*(Fp - Fm)/(dx(1)) + expAd*dx(1)**(-1)
                             ! if (augmentPlot .eq. 1) then 
                             !    agu(uc,8)=au
                             !    agu(uc,9) =wpl
                             !    agu(uc,10)=wpr
                             !    agu(uc,11)=wml
                             !    agu(uc,12)=wmr
                             ! end if
                             !v*vy---------------------------------------------------------------
                             ! Apl   = u(i1,i2+1,i3,vc) - 2.*u(i1,i2,i3,vc)   + u(i1,i2-1,i3,vc)
                             ! Bpl   = u(i1,i2+1,i3,vc) -    u(i1,i2-1,i3,vc)
                             ! Apr   = u(i1,i2+2,i3,vc) - 2.*u(i1,i2+1,i3,vc) + u(i1,i2,i3,vc)
                             ! Bpr   = u(i1,i2+2,i3,vc) -    u(i1,i2,i3,vc)
                             ! Aml   = u(i1,i2,i3,vc) - 2.*u(i1,i2-1,i3,vc) + u(i1,i2-2,i3,vc)
                             ! Bml   = u(i1,i2,i3,vc) -    u(i1,i2-2,i3,vc)
                             ! Amr   = Apl
                             ! Bmr   = Bpl
                             ! betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 1./4.*Bpl**2
                             ! betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 1./4.*Bpr**2
                             ! betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 1./4.*Bml**2
                             ! betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 1./4.*Bmr**2
                             ! ep     = 1e-15
                             ! ap1    = (1./2.)/(ep + betapl)**(4)
                             ! ap2    = (1./2.)/(ep + betapr)**(4)
                             ! am1    = (1./2.)/(ep + betaml)**(4)
                             ! am2    = (1./2.)/(ep + betamr)**(4)
                             ! !calculate weights
                             ! wp1 = ap1/(ap1+ap2)
                             ! wp2 = ap2/(ap1+ap2)
                             ! wm1 = am1/(am1+am2)
                             ! wm2 = am2/(am1+am2)
                             ! !mapping of the weights
                             ! bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                             ! bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                             ! bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                             ! bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                             ! wplw = bp1/(bp1+bp2)
                             ! wprw = bp2/(bp1+bp2)
                             ! wmlw = bm1/(bm1+bm2)
                             ! wmrw = bm2/(bm1+bm2)
                             ! !get face value of v
                             !  ! up = 0.5*(u(i1  ,i2,i3,vc)  + u(i1,i2+1,i3,vc))
                             !  ! um = 0.5*(u(i1-1,i2,i3,vc)  + u(i1,i2  ,i3,vc))
                             ! up = au
                             ! um = au
                             ! !picking side for the weights
                             !  if (up.ge.0) then
                             !    wpl = max(wplw,wprw)
                             !    wpr = min(wplw,wprw)
                             ! else
                             !    wpl = min(wplw,wprw)
                             !    wpr = max(wplw,wprw)
                             ! end if
                             ! if (um.ge.0) then
                             !    wml = max(wmlw,wmrw)
                             !    wmr = min(wmlw,wmrw)
                             ! else
                             !    wml = min(wmlw,wmrw)
                             !    wmr = max(wmlw,wmrw)
                             ! end if
                             ! expAd=0.
                             ! if (orderOneFix.eq.1) then
                             !    if (aum*aup .le. 0) then
                             !       aumax = max(abs(au),abs(aup),abs(aum))
                             !       expAd = (1./12.)*aumax*( u(i1,i2-2,i3,vc) - 4.*u(i1,i2-1,i3,vc) + 6.*u(i1,i2,i3,vc) - 4.*u(i1,i2+1,i3,vc) +  u(i1,i2+2,i3,vc))
                             !       wpl = 0.5
                             !       wpr = 0.5
                             !       wml = 0.5
                             !       wmr = 0.5
                             !    end if
                             ! end if
                             ! Fpl = 1./6.*(  -u(i1,i2-1,i3,vc) + 5.*u(i1,i2,i3,vc)   + 2.*u(i1,i2+1,i3,vc))
                             ! Fpr = 1./6.*( 2.*u(i1,i2,i3,vc)   + 5.*u(i1,i2+1,i3,vc) -   u(i1,i2+2,i3,vc))
                             ! Fml = 1./6.*(  -u(i1,i2-2,i3,vc) + 5.*u(i1,i2-1,i3,vc)   + 2.*u(i1,i2,i3,vc))
                             ! Fmr = 1./6.*( 2.*u(i1,i2-1,i3,vc) + 5.*u(i1,i2,i3,vc) -       u(i1,i2+1,i3,vc))
                             ! Fp  = wpl*Fpl + wpr*Fpr
                             ! Fm  = wml*Fml + wmr*Fmr
                             ! agu(vc,vc)= au*(Fp - Fm)/(dx(1)) + expAd*dx(1)**(-1)
                             ! if (augmentPlot .eq. 1) then  
                             !    agu(vc,9) =wpl
                             !    agu(vc,10)=wpr
                             !    agu(vc,11)=wml
                             !    agu(vc,12)=wmr
                             ! end if
                          else ! if (upwindOrder.eq.4)
                             write(*,'(" getBewnoAdvection: only 4th 
     & order is avaliable now")' )
                             stop 777
                          end if
                      ! #If ("BWENO" == "CENTERED") 
                      !  ! -- centered advection ---
                      !  ! write(*,'(" getAdvection -- centered")')
                      !  #If "2" == "2"
                      !   ! -- 2D --
                      !   agu(uc,uc)=UU(uc)*UX(uc)
                      !   agu(vc,uc)=UU(vc)*UY(uc)
                      !   agu(uc,vc)=UU(uc)*UX(vc)
                      !   agu(vc,vc)=UU(vc)*UY(vc)
                      !  #Else
                      !   ! -- 3D -- *check me*
                      !   agu(uc,uc)=UU(uc)*UX(uc)
                      !   agu(vc,uc)=UU(vc)*UY(uc)
                      !   agu(wc,uc)=UU(wc)*UY(uc)
                      !   agu(uc,vc)=UU(uc)*UX(vc)
                      !   agu(vc,vc)=UU(vc)*UY(vc)
                      !   agu(wc,vc)=UU(wc)*UY(vc)
                      !   agu(uc,wc)=UU(uc)*UX(wc)
                      !   agu(vc,wc)=UU(vc)*UY(wc)
                      !   agu(wc,wc)=UU(wc)*UY(wc)
                      !  #End
                      ! #Elif "BWENO" == "BWENO" 
                      !   ! --- upwind scheme ---
                      !   ! for testing output this next message:
                      !   if( t.le. 0. )then
                      !     write(*,'(" getAdvection upwind scheme (7)")') 
                      !   end if
                      !   getUpwindAdvection(u,i1,i2,i3,NONE,2,4,rectangular, agu)
                      ! #Elif "BWENO" == "BWENO" 
                      !   ! --- Bweno scheme ---
                      !   ! for testing output this message:
                      !   if( t.le. 0. )then
                      !      write(*,'(" getAdvection BWENO scheme (7)")') 
                      !   end if
                      !   getBwenoAdvection(u,i1,i2,i3,NONE,2,4,rectangular, agu)
                      ! #Else
                      !   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
                      !   stop 999
                      ! #End
                        ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-
     & ux42r(i1,i2,i3,pc)
                        ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-
     & uy42r(i1,i2,i3,pc)
                    ! include implicit terms - diffusion
                     uti(i1,i2,i3,uc)= nu*ulaplacian42r(i1,i2,i3,uc)
                     uti(i1,i2,i3,vc)= nu*ulaplacian42r(i1,i2,i3,vc)
                 end if
                end do
                end do
                end do
               else
                do i3=n3a,n3b
                do i2=n2a,n2b
                do i1=n1a,n1b
                   ! explicit terms only, no diffusion
                     ! upwind approximation/bweno
                       ! --- Bweno scheme ---
                       ! for testing output this message:
                       ! if( t.le. 0. )then
                       !    write(*,'(" getAdvection BWENO scheme (7)")') 
                       ! end if
                          gvU = u(i1,i2,i3,uc) - uu(i1,i2,i3,uc)
                          gvV = u(i1,i2,i3,vc) - uu(i1,i2,i3,vc)
                          gvW = u(i1,i2,i3,wc) - uu(i1,i2,i3,wc)
                          orderOneFix = 1
                          augmentPlot = 1
                          if( upwindOrder.eq.4 )then
                         ! -- fourth order bweno --
                             au = uu(i1,i2,i3,uc)
                             aum = u(i1-1,i2,i3,uc) - gvU
                             aup = u(i1+1,i2,i3,uc) - gvU
                             !usage
                             !s1   is the shifted grid in r1 direction
                             !var  is the variable taking derivative (uc,vc,wc)
                             !vard is the variable of direction
                             !agu is only passed in for augmented plot
                             !expAd 
                             ! u*ux
                             s1=1
                             s2=0
                             s3=0
                             var = uc
                             vard= uc
                             drl = dx(0)
                               Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*u(
     & i1   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
                               Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    u(
     & i1-s1,i2-s2,i3-s3,var)
                               Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) -
     &  2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
                               Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) -
     &     u(i1   ,i2   ,i3   ,var)
                               Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,
     & i2-s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
                               Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var)
                               Amr   = Apl
                               Bmr   = Bpl
                               betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 
     & 1./4.*Bpl**2
                               betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 
     & 1./4.*Bpr**2
                               betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 
     & 1./4.*Bml**2
                               betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 
     & 1./4.*Bmr**2
                               ep     = 1e-15
                               ap1    = (1./2.)/(ep + betapl)**(4)
                               ap2    = (1./2.)/(ep + betapr)**(4)
                               am1    = (1./2.)/(ep + betaml)**(4)
                               am2    = (1./2.)/(ep + betamr)**(4)
                               !calculate weights
                               wp1 = ap1/(ap1+ap2)
                               wp2 = ap2/(ap1+ap2)
                               wm1 = am1/(am1+am2)
                               wm2 = am2/(am1+am2)
                               !mapping of the weights
                               bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                               bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                               bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                               bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                               wplw = bp1/(bp1+bp2)
                               wprw = bp2/(bp1+bp2)
                               wmlw = bm1/(bm1+bm2)
                               wmrw = bm2/(bm1+bm2)
                               !get face value of u
                               up = au
                               um = au
                               !picking side for the weights
                               if (up.ge.0) then
                                  wpl = max(wplw,wprw)
                                  wpr = min(wplw,wprw)
                               else
                                  wpl = min(wplw,wprw)
                                  wpr = max(wplw,wprw)
                               end if
                               if (um.ge.0) then
                                  wml = max(wmlw,wmrw)
                                  wmr = min(wmlw,wmrw)
                               else
                                  wml = min(wmlw,wmrw)
                                  wmr = max(wmlw,wmrw)
                               end if
                               expAd = 0.
                               if (orderOneFix.eq.1) then
                                  ! if (augmentPlot.eq.1) then
                                  !    agu(uc,13)=0
                                  ! end if
                                  if (aum*aup .le. 0) then
                                     aumax = max(abs(au),abs(aup),abs(
     & aum))
                                     expAd = (1./12.)*aumax*( u(i1-2*
     & s1,i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(
     & i1,i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*
     & s2,i3+2*s3,var))
                                     wpl = 0.5
                                     wpr = 0.5
                                     wml = 0.5
                                     wmr = 0.5
                                     ! plot au at where fix is turned on
                                     ! if (augmentPlot.eq.1) then
                                     !    agu(uc,13)=au
                                     ! end if
                                  end if
                               end if
                               ! wpl = 0.5
                               ! wpr = 0.5
                               ! wml = 0.5
                               ! wmr = 0.5
                               Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,var)
     &  + 5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
                               Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 5.*
     & u(i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
                               Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*
     & s3,var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
                               Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,
     & var) + 5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
                               Fp  = wpl*Fpl + wpr*Fpr
                               Fm  = wml*Fml + wmr*Fmr
                                  agu(vard,var)= au*(Fp - Fm)/drl + 
     & expAd/drl
                             !     if ( i1 .eq. 14 .and. i2.eq.16 ) then
                                     !     if ( i1 .eq. 14) then
                                     ! write(*,'(" betapl,betapr,betaml,betamr = ",4F10.8)')  betapl,betapr,betaml,betamr
                                     ! write(*,'(" ap1,ap2,am1,am2 = ",4e10.8)')  ap1,ap2,am1,am2
                                     ! write(*,'(" wp1,wpr,wm1,wm2 = ",4F10.8)')  wp1,wp2,wm1,wm2
                                     ! write(*,'(" bp1,bp2,bm1,bm2 = ",4F10.8)')  bp1,bp2,bm1,bm2
                                     ! write(*,'(" wplw,wprw,wmlw,wmrw = ",4F10.8)')  wplw,wprw,wmlw,wmrw
                                     ! write(*,'(" wpl,wpr,wml,wmr = ",4F10.8)')  wpl,wpr,wml,wmr
                                     ! write(*,'(" agu = ",F10.8)') agu(vard,var)
                             !     end if
                             !u*vx
                             s1=1
                             s2=0
                             s3=0
                             var = vc
                             vard= uc
                             drl = dx(0)
                               Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*u(
     & i1   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
                               Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    u(
     & i1-s1,i2-s2,i3-s3,var)
                               Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) -
     &  2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
                               Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) -
     &     u(i1   ,i2   ,i3   ,var)
                               Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,
     & i2-s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
                               Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var)
                               Amr   = Apl
                               Bmr   = Bpl
                               betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 
     & 1./4.*Bpl**2
                               betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 
     & 1./4.*Bpr**2
                               betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 
     & 1./4.*Bml**2
                               betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 
     & 1./4.*Bmr**2
                               ep     = 1e-15
                               ap1    = (1./2.)/(ep + betapl)**(4)
                               ap2    = (1./2.)/(ep + betapr)**(4)
                               am1    = (1./2.)/(ep + betaml)**(4)
                               am2    = (1./2.)/(ep + betamr)**(4)
                               !calculate weights
                               wp1 = ap1/(ap1+ap2)
                               wp2 = ap2/(ap1+ap2)
                               wm1 = am1/(am1+am2)
                               wm2 = am2/(am1+am2)
                               !mapping of the weights
                               bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                               bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                               bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                               bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                               wplw = bp1/(bp1+bp2)
                               wprw = bp2/(bp1+bp2)
                               wmlw = bm1/(bm1+bm2)
                               wmrw = bm2/(bm1+bm2)
                               !get face value of u
                               up = au
                               um = au
                               !picking side for the weights
                               if (up.ge.0) then
                                  wpl = max(wplw,wprw)
                                  wpr = min(wplw,wprw)
                               else
                                  wpl = min(wplw,wprw)
                                  wpr = max(wplw,wprw)
                               end if
                               if (um.ge.0) then
                                  wml = max(wmlw,wmrw)
                                  wmr = min(wmlw,wmrw)
                               else
                                  wml = min(wmlw,wmrw)
                                  wmr = max(wmlw,wmrw)
                               end if
                               expAd = 0.
                               if (orderOneFix.eq.1) then
                                  ! if (augmentPlot.eq.1) then
                                  !    agu(uc,13)=0
                                  ! end if
                                  if (aum*aup .le. 0) then
                                     aumax = max(abs(au),abs(aup),abs(
     & aum))
                                     expAd = (1./12.)*aumax*( u(i1-2*
     & s1,i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(
     & i1,i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*
     & s2,i3+2*s3,var))
                                     wpl = 0.5
                                     wpr = 0.5
                                     wml = 0.5
                                     wmr = 0.5
                                     ! plot au at where fix is turned on
                                     ! if (augmentPlot.eq.1) then
                                     !    agu(uc,13)=au
                                     ! end if
                                  end if
                               end if
                               ! wpl = 0.5
                               ! wpr = 0.5
                               ! wml = 0.5
                               ! wmr = 0.5
                               Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,var)
     &  + 5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
                               Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 5.*
     & u(i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
                               Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*
     & s3,var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
                               Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,
     & var) + 5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
                               Fp  = wpl*Fpl + wpr*Fpr
                               Fm  = wml*Fml + wmr*Fmr
                                  agu(vard,var)= au*(Fp - Fm)/drl + 
     & expAd/drl
                             !     if ( i1 .eq. 14 .and. i2.eq.16 ) then
                                     !     if ( i1 .eq. 14) then
                                     ! write(*,'(" betapl,betapr,betaml,betamr = ",4F10.8)')  betapl,betapr,betaml,betamr
                                     ! write(*,'(" ap1,ap2,am1,am2 = ",4e10.8)')  ap1,ap2,am1,am2
                                     ! write(*,'(" wp1,wpr,wm1,wm2 = ",4F10.8)')  wp1,wp2,wm1,wm2
                                     ! write(*,'(" bp1,bp2,bm1,bm2 = ",4F10.8)')  bp1,bp2,bm1,bm2
                                     ! write(*,'(" wplw,wprw,wmlw,wmrw = ",4F10.8)')  wplw,wprw,wmlw,wmrw
                                     ! write(*,'(" wpl,wpr,wml,wmr = ",4F10.8)')  wpl,wpr,wml,wmr
                                     ! write(*,'(" agu = ",F10.8)') agu(vard,var)
                             !     end if
                             !v*uy
                             au = uu(i1,i2,i3,vc)
                             aum = u(i1,i2-1,i3,vc) - gvV
                             aup = u(i1,i2+1,i3,vc) - gvV
                             s1=0
                             s2=1
                             s3=0
                             var = uc
                             vard= vc
                             drl = dx(1)
                               Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*u(
     & i1   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
                               Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    u(
     & i1-s1,i2-s2,i3-s3,var)
                               Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) -
     &  2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
                               Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) -
     &     u(i1   ,i2   ,i3   ,var)
                               Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,
     & i2-s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
                               Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var)
                               Amr   = Apl
                               Bmr   = Bpl
                               betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 
     & 1./4.*Bpl**2
                               betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 
     & 1./4.*Bpr**2
                               betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 
     & 1./4.*Bml**2
                               betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 
     & 1./4.*Bmr**2
                               ep     = 1e-15
                               ap1    = (1./2.)/(ep + betapl)**(4)
                               ap2    = (1./2.)/(ep + betapr)**(4)
                               am1    = (1./2.)/(ep + betaml)**(4)
                               am2    = (1./2.)/(ep + betamr)**(4)
                               !calculate weights
                               wp1 = ap1/(ap1+ap2)
                               wp2 = ap2/(ap1+ap2)
                               wm1 = am1/(am1+am2)
                               wm2 = am2/(am1+am2)
                               !mapping of the weights
                               bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                               bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                               bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                               bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                               wplw = bp1/(bp1+bp2)
                               wprw = bp2/(bp1+bp2)
                               wmlw = bm1/(bm1+bm2)
                               wmrw = bm2/(bm1+bm2)
                               !get face value of u
                               up = au
                               um = au
                               !picking side for the weights
                               if (up.ge.0) then
                                  wpl = max(wplw,wprw)
                                  wpr = min(wplw,wprw)
                               else
                                  wpl = min(wplw,wprw)
                                  wpr = max(wplw,wprw)
                               end if
                               if (um.ge.0) then
                                  wml = max(wmlw,wmrw)
                                  wmr = min(wmlw,wmrw)
                               else
                                  wml = min(wmlw,wmrw)
                                  wmr = max(wmlw,wmrw)
                               end if
                               expAd = 0.
                               if (orderOneFix.eq.1) then
                                  ! if (augmentPlot.eq.1) then
                                  !    agu(uc,13)=0
                                  ! end if
                                  if (aum*aup .le. 0) then
                                     aumax = max(abs(au),abs(aup),abs(
     & aum))
                                     expAd = (1./12.)*aumax*( u(i1-2*
     & s1,i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(
     & i1,i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*
     & s2,i3+2*s3,var))
                                     wpl = 0.5
                                     wpr = 0.5
                                     wml = 0.5
                                     wmr = 0.5
                                     ! plot au at where fix is turned on
                                     ! if (augmentPlot.eq.1) then
                                     !    agu(uc,13)=au
                                     ! end if
                                  end if
                               end if
                               ! wpl = 0.5
                               ! wpr = 0.5
                               ! wml = 0.5
                               ! wmr = 0.5
                               Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,var)
     &  + 5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
                               Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 5.*
     & u(i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
                               Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*
     & s3,var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
                               Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,
     & var) + 5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
                               Fp  = wpl*Fpl + wpr*Fpr
                               Fm  = wml*Fml + wmr*Fmr
                                  agu(vard,var)= au*(Fp - Fm)/drl + 
     & expAd/drl
                             !     if ( i1 .eq. 14 .and. i2.eq.16 ) then
                                     !     if ( i1 .eq. 14) then
                                     ! write(*,'(" betapl,betapr,betaml,betamr = ",4F10.8)')  betapl,betapr,betaml,betamr
                                     ! write(*,'(" ap1,ap2,am1,am2 = ",4e10.8)')  ap1,ap2,am1,am2
                                     ! write(*,'(" wp1,wpr,wm1,wm2 = ",4F10.8)')  wp1,wp2,wm1,wm2
                                     ! write(*,'(" bp1,bp2,bm1,bm2 = ",4F10.8)')  bp1,bp2,bm1,bm2
                                     ! write(*,'(" wplw,wprw,wmlw,wmrw = ",4F10.8)')  wplw,wprw,wmlw,wmrw
                                     ! write(*,'(" wpl,wpr,wml,wmr = ",4F10.8)')  wpl,wpr,wml,wmr
                                     ! write(*,'(" agu = ",F10.8)') agu(vard,var)
                             !     end if
                             !v*vy
                             s1=0
                             s2=1
                             s3=0
                             var = vc
                             vard= vc
                             drl = dx(1)
                               Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*u(
     & i1   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
                               Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    u(
     & i1-s1,i2-s2,i3-s3,var)
                               Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) -
     &  2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
                               Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) -
     &     u(i1   ,i2   ,i3   ,var)
                               Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,
     & i2-s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
                               Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var)
                               Amr   = Apl
                               Bmr   = Bpl
                               betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 
     & 1./4.*Bpl**2
                               betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 
     & 1./4.*Bpr**2
                               betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 
     & 1./4.*Bml**2
                               betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 
     & 1./4.*Bmr**2
                               ep     = 1e-15
                               ap1    = (1./2.)/(ep + betapl)**(4)
                               ap2    = (1./2.)/(ep + betapr)**(4)
                               am1    = (1./2.)/(ep + betaml)**(4)
                               am2    = (1./2.)/(ep + betamr)**(4)
                               !calculate weights
                               wp1 = ap1/(ap1+ap2)
                               wp2 = ap2/(ap1+ap2)
                               wm1 = am1/(am1+am2)
                               wm2 = am2/(am1+am2)
                               !mapping of the weights
                               bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                               bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                               bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                               bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                               wplw = bp1/(bp1+bp2)
                               wprw = bp2/(bp1+bp2)
                               wmlw = bm1/(bm1+bm2)
                               wmrw = bm2/(bm1+bm2)
                               !get face value of u
                               up = au
                               um = au
                               !picking side for the weights
                               if (up.ge.0) then
                                  wpl = max(wplw,wprw)
                                  wpr = min(wplw,wprw)
                               else
                                  wpl = min(wplw,wprw)
                                  wpr = max(wplw,wprw)
                               end if
                               if (um.ge.0) then
                                  wml = max(wmlw,wmrw)
                                  wmr = min(wmlw,wmrw)
                               else
                                  wml = min(wmlw,wmrw)
                                  wmr = max(wmlw,wmrw)
                               end if
                               expAd = 0.
                               if (orderOneFix.eq.1) then
                                  ! if (augmentPlot.eq.1) then
                                  !    agu(uc,13)=0
                                  ! end if
                                  if (aum*aup .le. 0) then
                                     aumax = max(abs(au),abs(aup),abs(
     & aum))
                                     expAd = (1./12.)*aumax*( u(i1-2*
     & s1,i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(
     & i1,i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*
     & s2,i3+2*s3,var))
                                     wpl = 0.5
                                     wpr = 0.5
                                     wml = 0.5
                                     wmr = 0.5
                                     ! plot au at where fix is turned on
                                     ! if (augmentPlot.eq.1) then
                                     !    agu(uc,13)=au
                                     ! end if
                                  end if
                               end if
                               ! wpl = 0.5
                               ! wpr = 0.5
                               ! wml = 0.5
                               ! wmr = 0.5
                               Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,var)
     &  + 5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
                               Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 5.*
     & u(i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
                               Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*
     & s3,var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
                               Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,
     & var) + 5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
                               Fp  = wpl*Fpl + wpr*Fpr
                               Fm  = wml*Fml + wmr*Fmr
                                  agu(vard,var)= au*(Fp - Fm)/drl + 
     & expAd/drl
                             !     if ( i1 .eq. 14 .and. i2.eq.16 ) then
                                     !     if ( i1 .eq. 14) then
                                     ! write(*,'(" betapl,betapr,betaml,betamr = ",4F10.8)')  betapl,betapr,betaml,betamr
                                     ! write(*,'(" ap1,ap2,am1,am2 = ",4e10.8)')  ap1,ap2,am1,am2
                                     ! write(*,'(" wp1,wpr,wm1,wm2 = ",4F10.8)')  wp1,wp2,wm1,wm2
                                     ! write(*,'(" bp1,bp2,bm1,bm2 = ",4F10.8)')  bp1,bp2,bm1,bm2
                                     ! write(*,'(" wplw,wprw,wmlw,wmrw = ",4F10.8)')  wplw,wprw,wmlw,wmrw
                                     ! write(*,'(" wpl,wpr,wml,wmr = ",4F10.8)')  wpl,wpr,wml,wmr
                                     ! write(*,'(" agu = ",F10.8)') agu(vard,var)
                             !     end if
                            ! u*ux-------------------------------------------------------------------------
                            !calculate smooth factor
                            ! Apl   = u(i1+1,i2,i3,uc) - 2.*u(i1  ,i2,i3,uc)   + u(i1-1,i2,i3,uc)
                            ! Bpl   = u(i1+1,i2,i3,uc) -    u(i1-1,i2,i3,uc)
                            ! Apr   = u(i1+2,i2,i3,uc) - 2.*u(i1+1,i2,i3,uc) + u(i1,i2,i3,uc)
                            ! Bpr   = u(i1+2,i2,i3,uc) -    u(i1  ,i2,i3,uc)
                            ! Aml   = u(i1,i2,i3,uc) - 2.*u(i1-1,i2,i3,uc) + u(i1-2,i2,i3,uc)
                            ! Bml   = u(i1,i2,i3,uc) -    u(i1-2,i2,i3,uc)
                            ! Amr   = Apl
                            ! Bmr   = Bpl
                            ! betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 1./4.*Bpl**2
                            ! betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 1./4.*Bpr**2
                            ! betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 1./4.*Bml**2
                            ! betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 1./4.*Bmr**2
                            ! ep     = 1e-15
                            ! ap1    = (1./2.)/(ep + betapl)**(4)
                            ! ap2    = (1./2.)/(ep + betapr)**(4)
                            ! am1    = (1./2.)/(ep + betaml)**(4)
                            ! am2    = (1./2.)/(ep + betamr)**(4)
                            ! !calculate weights
                            ! wp1 = ap1/(ap1+ap2)
                            ! wp2 = ap2/(ap1+ap2)
                            ! wm1 = am1/(am1+am2)
                            ! wm2 = am2/(am1+am2)
                            ! !mapping of the weights
                            ! bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                            ! bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                            ! bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                            ! bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                            ! wplw = bp1/(bp1+bp2)
                            ! wprw = bp2/(bp1+bp2)
                            ! wmlw = bm1/(bm1+bm2)
                            ! wmrw = bm2/(bm1+bm2)
                            ! !get face value of u
                            ! up = au
                            ! um = au
                            ! !picking side for the weights
                            !  if (up.ge.0) then
                            !    wpl = max(wplw,wprw)
                            !    wpr = min(wplw,wprw)
                            ! else
                            !    wpl = min(wplw,wprw)
                            !    wpr = max(wplw,wprw)
                            ! end if
                            ! if (um.ge.0) then
                            !    wml = max(wmlw,wmrw)
                            !    wmr = min(wmlw,wmrw)
                            ! else
                            !    wml = min(wmlw,wmrw)
                            !    wmr = max(wmlw,wmrw)
                            ! end if
                            ! expAd = 0.
                            ! if (orderOneFix.eq.1) then 
                            !    if (augmentPlot.eq.1) then
                            !       agu(uc,13)=0
                            !    end if
                            !    if (aum*aup .le. 0) then
                            !       aumax = max(abs(au),abs(aup),abs(aum))
                            !       expAd = (1./12.)*aumax*( u(i1-2,i2,i3,uc) - 4.*u(i1-1,i2,i3,uc) + 6.*u(i1,i2,i3,uc) - 4.*u(i1+1,i2,i3,uc) +  u(i1+2,i2,i3,uc))
                            !       wpl = 0.5
                            !       wpr = 0.5
                            !       wml = 0.5
                            !       wmr = 0.5
                            !       ! plot au at where fix is turned on
                            !       if (augmentPlot.eq.1) then
                            !          agu(uc,13)=au
                            !       end if
                            !    end if
                            ! end if
                            ! wpl = 0.5
                            ! wpr = 0.5
                            ! wml = 0.5
                            ! wmr = 0.5
                            ! Fpl = 1./6.*(  -u(i1-1,i2,i3,uc) + 5.*u(i1,i2,i3,uc)   + 2.*u(i1+1,i2,i3,uc))
                            ! Fpr = 1./6.*( 2.*u(i1,i2,i3,uc)   + 5.*u(i1+1,i2,i3,uc) -   u(i1+2,i2,i3,uc))
                            ! Fml = 1./6.*(  -u(i1-2,i2,i3,uc) + 5.*u(i1-1,i2,i3,uc)   + 2.*u(i1,i2,i3,uc))
                            ! Fmr = 1./6.*( 2.*u(i1-1,i2,i3,uc) + 5.*u(i1,i2,i3,uc) -       u(i1+1,i2,i3,uc))
                            ! Fp  = wpl*Fpl + wpr*Fpr
                            ! Fm  = wml*Fml + wmr*Fmr
                            ! agu(uc,uc)= au*(Fp - Fm)/(dx(0)) + expAd*dx(0)**(-1)
                            ! if (augmentPlot.eq.1) then
                            !    agu(uc,3)=au
                            !    agu(uc,4)=wpl
                            !    agu(uc,5)=wpr
                            !    agu(uc,6)=wml
                            !    agu(uc,7)=wmr
                            ! end if
                            ! u*vx-------------------------------------------------------------------------
                            !calculate smooth factor
                            ! Apl   = u(i1+1,i2,i3,vc) - 2.*u(i1,i2,i3,vc)   + u(i1-1,i2,i3,vc)
                            ! Bpl   = u(i1+1,i2,i3,vc) -    u(i1-1,i2,i3,vc)
                            ! Apr   = u(i1+2,i2,i3,vc) - 2.*u(i1+1,i2,i3,vc) + u(i1,i2,i3,vc)
                            ! Bpr   = u(i1+2,i2,i3,vc) -    u(i1,i2,i3,vc)
                            ! Aml   = u(i1,i2,i3,vc) - 2.*u(i1-1,i2,i3,vc) + u(i1-2,i2,i3,vc)
                            ! Bml   = u(i1,i2,i3,vc) -    u(i1-2,i2,i3,vc)
                            ! Amr   = Apl
                            ! Bmr   = Bpl
                            ! betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 1./4.*Bpl**2
                            ! betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 1./4.*Bpr**2
                            ! betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 1./4.*Bml**2
                            ! betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 1./4.*Bmr**2
                            ! ep     = 1e-15
                            ! ap1    = (1./2.)/(ep + betapl)**(4)
                            ! ap2    = (1./2.)/(ep + betapr)**(4)
                            ! am1    = (1./2.)/(ep + betaml)**(4)
                            ! am2    = (1./2.)/(ep + betamr)**(4)
                            ! !calculate weights
                            ! wp1 = ap1/(ap1+ap2)
                            ! wp2 = ap2/(ap1+ap2)
                            ! wm1 = am1/(am1+am2)
                            ! wm2 = am2/(am1+am2)
                            ! !mapping of the weights
                            ! bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                            ! bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                            ! bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                            ! bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                            ! wplw = bp1/(bp1+bp2)
                            ! wprw = bp2/(bp1+bp2)
                            ! wmlw = bm1/(bm1+bm2)
                            ! wmrw = bm2/(bm1+bm2)
                            ! !get face value of u
                            ! ! up = 0.5*(u(i1  ,i2,i3,uc)+u(i1+1,i2,i3,uc))
                            ! ! um = 0.5*(u(i1-1,i2,i3,uc)+u(i1  ,i2,i3,uc))
                            ! up = au
                            ! um = au
                            ! !picking side for the weights
                            ! if (up.ge.0) then
                            !    wpl = max(wplw,wprw)
                            !    wpr = min(wplw,wprw)
                            ! else
                            !    wpl = min(wplw,wprw)
                            !    wpr = max(wplw,wprw)
                            ! end if
                            ! if (um.ge.0) then
                            !    wml = max(wmlw,wmrw)
                            !    wmr = min(wmlw,wmrw)
                            ! else
                            !    wml = min(wmlw,wmrw)
                            !    wmr = max(wmlw,wmrw)
                            ! end if
                            ! expAd = 0.
                            ! if (orderOneFix.eq.1) then
                            !    if (aum*aup .le. 0) then
                            !       aumax = max(abs(au),abs(aup),abs(aum))
                            !       expAd = (1./12.)*aumax*( u(i1-2,i2,i3,vc) - 4.*u(i1-1,i2,i3,vc) + 6.*u(i1,i2,i3,vc) - 4.*u(i1+1,i2,i3,vc) +  u(i1+2,i2,i3,vc))
                            !       wpl = 0.5
                            !       wpr = 0.5
                            !       wml = 0.5
                            !       wmr = 0.5
                            !    end if
                            ! end if
                            ! Fpl = 1./6.*(  -u(i1-1,i2,i3,vc) + 5.*u(i1,i2,i3,vc)   + 2.*u(i1+1,i2,i3,vc))
                            ! Fpr = 1./6.*( 2.*u(i1,i2,i3,vc)   + 5.*u(i1+1,i2,i3,vc) -   u(i1+2,i2,i3,vc))
                            ! Fml = 1./6.*(  -u(i1-2,i2,i3,vc) + 5.*u(i1-1,i2,i3,vc)   + 2.*u(i1,i2,i3,vc))
                            ! Fmr = 1./6.*( 2.*u(i1-1,i2,i3,vc) + 5.*u(i1,i2,i3,vc) -       u(i1+1,i2,i3,vc))
                            ! Fp  = wpl*Fpl + wpr*Fpr
                            ! Fm  = wml*Fml + wmr*Fmr
                            ! agu(uc,vc)= au*(Fp - Fm)/(dx(0)) + expAd*dx(0)**(-1) 
                            ! if (augmentPlot .eq. 1) then 
                            !    agu(vc,4)=wpl
                            !    agu(vc,5)=wpr
                            !    agu(vc,6)=wml
                            !    agu(vc,7)=wmr
                            ! end if
                            !v*uy-------------------------------------------------------------------------
                            ! au = u(i1,i2,i3,vc)
                            ! aum = u(i1,i2-1,i3,vc)
                            ! aup = u(i1,i2+1,i3,vc)
                            !calculate smooth factor
                            ! Apl   = u(i1,i2+1,i3,uc) - 2.*u(i1,i2,i3,uc)   + u(i1,i2-1,i3,uc)
                            ! Bpl   = u(i1,i2+1,i3,uc) -    u(i1,i2-1,i3,uc)
                            ! Apr   = u(i1,i2+2,i3,uc) - 2.*u(i1,i2+1,i3,uc) + u(i1,i2,i3,uc)
                            ! Bpr   = u(i1,i2+2,i3,uc) -    u(i1,i2,i3,uc)
                            ! Aml   = u(i1,i2,i3,uc) - 2.*u(i1,i2-1,i3,uc) + u(i1,i2-2,i3,uc)
                            ! Bml   = u(i1,i2,i3,uc) -    u(i1,i2-2,i3,uc)
                            ! Amr   = Apl
                            ! Bmr   = Bpl
                            ! betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 1./4.*Bpl**2
                            ! betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 1./4.*Bpr**2
                            ! betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 1./4.*Bml**2
                            ! betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 1./4.*Bmr**2
                            ! ep     = 1e-15
                            ! ap1    = (1./2.)/(ep + betapl)**(4)
                            ! ap2    = (1./2.)/(ep + betapr)**(4)
                            ! am1    = (1./2.)/(ep + betaml)**(4)
                            ! am2    = (1./2.)/(ep + betamr)**(4)
                            ! !calculate weights
                            ! wp1 = ap1/(ap1+ap2)
                            ! wp2 = ap2/(ap1+ap2)
                            ! wm1 = am1/(am1+am2)
                            ! wm2 = am2/(am1+am2)
                            ! !mapping of the weights
                            ! bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                            ! bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                            ! bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                            ! bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                            ! wplw = bp1/(bp1+bp2)
                            ! wprw = bp2/(bp1+bp2)
                            ! wmlw = bm1/(bm1+bm2)
                            ! wmrw = bm2/(bm1+bm2)
                            ! !get face value of v
                            ! ! up = 0.5*(u(i1  ,i2+1,i3,vc)+u(i1,i2 ,i3,vc))
                            ! ! um = 0.5*(u(i1  ,i2  ,i3,vc)+u(i1,i2-1,i3,vc))
                            ! up = au
                            ! um = au
                            ! !picking side for the weights
                            !  if (up.ge.0) then
                            !    wpl = max(wplw,wprw)
                            !    wpr = min(wplw,wprw)
                            ! else
                            !    wpl = min(wplw,wprw)
                            !    wpr = max(wplw,wprw)
                            ! end if
                            ! if (um.ge.0) then
                            !    wml = max(wmlw,wmrw)
                            !    wmr = min(wmlw,wmrw)
                            ! else
                            !    wml = min(wmlw,wmrw)
                            !    wmr = max(wmlw,wmrw)
                            ! end if
                            ! expAd = 0.
                            ! if (orderOneFix.eq.1) then 
                            !    if (augmentPlot .eq. 1) then 
                            !       agu(vc,13)=0
                            !    end if
                            !    if (aum*aup .le. 0) then
                            !       aumax = max(abs(au),abs(aup),abs(aum))
                            !       expAd = (1./12.)*aumax*( u(i1,i2-2,i3,uc) - 4.*u(i1,i2-1,i3,uc) + 6.*u(i1,i2,i3,uc) - 4.*u(i1,i2+1,i3,uc) +  u(i1,i2+2,i3,uc))
                            !       wpl = 0.5
                            !       wpr = 0.5
                            !       wml = 0.5
                            !       wmr = 0.5
                            !       if (augmentPlot .eq. 1) then 
                            !          agu(vc,13)=au
                            !       end if
                            !    end if
                            ! end if
                            ! Fpl = 1./6.*(  -u(i1,i2-1,i3,uc) + 5.*u(i1,i2,i3,uc)   + 2.*u(i1,i2+1,i3,uc))
                            ! Fpr = 1./6.*( 2.*u(i1,i2,i3,uc)   + 5.*u(i1,i2+1,i3,uc) -   u(i1,i2+2,i3,uc))
                            ! Fml = 1./6.*(  -u(i1,i2-2,i3,uc) + 5.*u(i1,i2-1,i3,uc)   + 2.*u(i1,i2,i3,uc))
                            ! Fmr = 1./6.*( 2.*u(i1,i2-1,i3,uc) + 5.*u(i1,i2,i3,uc) -       u(i1,i2+1,i3,uc))
                            ! Fp  = wpl*Fpl + wpr*Fpr
                            ! Fm  = wml*Fml + wmr*Fmr
                            ! agu(vc,uc)= au*(Fp - Fm)/(dx(1)) + expAd*dx(1)**(-1)
                            ! if (augmentPlot .eq. 1) then 
                            !    agu(uc,8)=au
                            !    agu(uc,9) =wpl
                            !    agu(uc,10)=wpr
                            !    agu(uc,11)=wml
                            !    agu(uc,12)=wmr
                            ! end if
                            !v*vy---------------------------------------------------------------
                            ! Apl   = u(i1,i2+1,i3,vc) - 2.*u(i1,i2,i3,vc)   + u(i1,i2-1,i3,vc)
                            ! Bpl   = u(i1,i2+1,i3,vc) -    u(i1,i2-1,i3,vc)
                            ! Apr   = u(i1,i2+2,i3,vc) - 2.*u(i1,i2+1,i3,vc) + u(i1,i2,i3,vc)
                            ! Bpr   = u(i1,i2+2,i3,vc) -    u(i1,i2,i3,vc)
                            ! Aml   = u(i1,i2,i3,vc) - 2.*u(i1,i2-1,i3,vc) + u(i1,i2-2,i3,vc)
                            ! Bml   = u(i1,i2,i3,vc) -    u(i1,i2-2,i3,vc)
                            ! Amr   = Apl
                            ! Bmr   = Bpl
                            ! betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 1./4.*Bpl**2
                            ! betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 1./4.*Bpr**2
                            ! betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 1./4.*Bml**2
                            ! betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 1./4.*Bmr**2
                            ! ep     = 1e-15
                            ! ap1    = (1./2.)/(ep + betapl)**(4)
                            ! ap2    = (1./2.)/(ep + betapr)**(4)
                            ! am1    = (1./2.)/(ep + betaml)**(4)
                            ! am2    = (1./2.)/(ep + betamr)**(4)
                            ! !calculate weights
                            ! wp1 = ap1/(ap1+ap2)
                            ! wp2 = ap2/(ap1+ap2)
                            ! wm1 = am1/(am1+am2)
                            ! wm2 = am2/(am1+am2)
                            ! !mapping of the weights
                            ! bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                            ! bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                            ! bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                            ! bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                            ! wplw = bp1/(bp1+bp2)
                            ! wprw = bp2/(bp1+bp2)
                            ! wmlw = bm1/(bm1+bm2)
                            ! wmrw = bm2/(bm1+bm2)
                            ! !get face value of v
                            !  ! up = 0.5*(u(i1  ,i2,i3,vc)  + u(i1,i2+1,i3,vc))
                            !  ! um = 0.5*(u(i1-1,i2,i3,vc)  + u(i1,i2  ,i3,vc))
                            ! up = au
                            ! um = au
                            ! !picking side for the weights
                            !  if (up.ge.0) then
                            !    wpl = max(wplw,wprw)
                            !    wpr = min(wplw,wprw)
                            ! else
                            !    wpl = min(wplw,wprw)
                            !    wpr = max(wplw,wprw)
                            ! end if
                            ! if (um.ge.0) then
                            !    wml = max(wmlw,wmrw)
                            !    wmr = min(wmlw,wmrw)
                            ! else
                            !    wml = min(wmlw,wmrw)
                            !    wmr = max(wmlw,wmrw)
                            ! end if
                            ! expAd=0.
                            ! if (orderOneFix.eq.1) then
                            !    if (aum*aup .le. 0) then
                            !       aumax = max(abs(au),abs(aup),abs(aum))
                            !       expAd = (1./12.)*aumax*( u(i1,i2-2,i3,vc) - 4.*u(i1,i2-1,i3,vc) + 6.*u(i1,i2,i3,vc) - 4.*u(i1,i2+1,i3,vc) +  u(i1,i2+2,i3,vc))
                            !       wpl = 0.5
                            !       wpr = 0.5
                            !       wml = 0.5
                            !       wmr = 0.5
                            !    end if
                            ! end if
                            ! Fpl = 1./6.*(  -u(i1,i2-1,i3,vc) + 5.*u(i1,i2,i3,vc)   + 2.*u(i1,i2+1,i3,vc))
                            ! Fpr = 1./6.*( 2.*u(i1,i2,i3,vc)   + 5.*u(i1,i2+1,i3,vc) -   u(i1,i2+2,i3,vc))
                            ! Fml = 1./6.*(  -u(i1,i2-2,i3,vc) + 5.*u(i1,i2-1,i3,vc)   + 2.*u(i1,i2,i3,vc))
                            ! Fmr = 1./6.*( 2.*u(i1,i2-1,i3,vc) + 5.*u(i1,i2,i3,vc) -       u(i1,i2+1,i3,vc))
                            ! Fp  = wpl*Fpl + wpr*Fpr
                            ! Fm  = wml*Fml + wmr*Fmr
                            ! agu(vc,vc)= au*(Fp - Fm)/(dx(1)) + expAd*dx(1)**(-1)
                            ! if (augmentPlot .eq. 1) then  
                            !    agu(vc,9) =wpl
                            !    agu(vc,10)=wpr
                            !    agu(vc,11)=wml
                            !    agu(vc,12)=wmr
                            ! end if
                         else ! if (upwindOrder.eq.4)
                            write(*,'(" getBewnoAdvection: only 4th 
     & order is avaliable now")' )
                            stop 777
                         end if
                     ! #If ("BWENO" == "CENTERED") 
                     !  ! -- centered advection ---
                     !  ! write(*,'(" getAdvection -- centered")')
                     !  #If "2" == "2"
                     !   ! -- 2D --
                     !   agu(uc,uc)=UU(uc)*UX(uc)
                     !   agu(vc,uc)=UU(vc)*UY(uc)
                     !   agu(uc,vc)=UU(uc)*UX(vc)
                     !   agu(vc,vc)=UU(vc)*UY(vc)
                     !  #Else
                     !   ! -- 3D -- *check me*
                     !   agu(uc,uc)=UU(uc)*UX(uc)
                     !   agu(vc,uc)=UU(vc)*UY(uc)
                     !   agu(wc,uc)=UU(wc)*UY(uc)
                     !   agu(uc,vc)=UU(uc)*UX(vc)
                     !   agu(vc,vc)=UU(vc)*UY(vc)
                     !   agu(wc,vc)=UU(wc)*UY(vc)
                     !   agu(uc,wc)=UU(uc)*UX(wc)
                     !   agu(vc,wc)=UU(vc)*UY(wc)
                     !   agu(wc,wc)=UU(wc)*UY(wc)
                     !  #End
                     ! #Elif "BWENO" == "BWENO" 
                     !   ! --- upwind scheme ---
                     !   ! for testing output this next message:
                     !   if( t.le. 0. )then
                     !     write(*,'(" getAdvection upwind scheme (7)")') 
                     !   end if
                     !   getUpwindAdvection(u,i1,i2,i3,NONE,2,4,rectangular, agu)
                     ! #Elif "BWENO" == "BWENO" 
                     !   ! --- Bweno scheme ---
                     !   ! for testing output this message:
                     !   if( t.le. 0. )then
                     !      write(*,'(" getAdvection BWENO scheme (7)")') 
                     !   end if
                     !   getBwenoAdvection(u,i1,i2,i3,NONE,2,4,rectangular, agu)
                     ! #Else
                     !   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
                     !   stop 999
                     ! #End
                       ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-ux42r(
     & i1,i2,i3,pc)
                       ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-uy42r(
     & i1,i2,i3,pc)
                   ! include implicit terms - diffusion
                    uti(i1,i2,i3,uc)= nu*ulaplacian42r(i1,i2,i3,uc)
                    uti(i1,i2,i3,vc)= nu*ulaplacian42r(i1,i2,i3,vc)
                end do
                end do
                end do
               end if
             else if( implicitOption.eq.doNotComputeImplicitTerms )then
               if( useWhereMask.ne.0 )then
                do i3=n3a,n3b
                do i2=n2a,n2b
                do i1=n1a,n1b
                 if( mask(i1,i2,i3).gt.0 )then
                    ! explicit terms only, no diffusion
                      ! upwind approximation/bweno
                        ! --- Bweno scheme ---
                        ! for testing output this message:
                        ! if( t.le. 0. )then
                        !    write(*,'(" getAdvection BWENO scheme (7)")') 
                        ! end if
                           gvU = u(i1,i2,i3,uc) - uu(i1,i2,i3,uc)
                           gvV = u(i1,i2,i3,vc) - uu(i1,i2,i3,vc)
                           gvW = u(i1,i2,i3,wc) - uu(i1,i2,i3,wc)
                           orderOneFix = 1
                           augmentPlot = 1
                           if( upwindOrder.eq.4 )then
                          ! -- fourth order bweno --
                              au = uu(i1,i2,i3,uc)
                              aum = u(i1-1,i2,i3,uc) - gvU
                              aup = u(i1+1,i2,i3,uc) - gvU
                              !usage
                              !s1   is the shifted grid in r1 direction
                              !var  is the variable taking derivative (uc,vc,wc)
                              !vard is the variable of direction
                              !agu is only passed in for augmented plot
                              !expAd 
                              ! u*ux
                              s1=1
                              s2=0
                              s3=0
                              var = uc
                              vard= uc
                              drl = dx(0)
                                Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*
     & u(i1   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
                                Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    
     & u(i1-s1,i2-s2,i3-s3,var)
                                Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) 
     & - 2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
                                Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) 
     & -    u(i1   ,i2   ,i3   ,var)
                                Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,
     & i2-s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
                                Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var)
                                Amr   = Apl
                                Bmr   = Bpl
                                betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 
     & 1./4.*Bpl**2
                                betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 
     & 1./4.*Bpr**2
                                betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 
     & 1./4.*Bml**2
                                betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 
     & 1./4.*Bmr**2
                                ep     = 1e-15
                                ap1    = (1./2.)/(ep + betapl)**(4)
                                ap2    = (1./2.)/(ep + betapr)**(4)
                                am1    = (1./2.)/(ep + betaml)**(4)
                                am2    = (1./2.)/(ep + betamr)**(4)
                                !calculate weights
                                wp1 = ap1/(ap1+ap2)
                                wp2 = ap2/(ap1+ap2)
                                wm1 = am1/(am1+am2)
                                wm2 = am2/(am1+am2)
                                !mapping of the weights
                                bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                                bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                                bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                                bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                                wplw = bp1/(bp1+bp2)
                                wprw = bp2/(bp1+bp2)
                                wmlw = bm1/(bm1+bm2)
                                wmrw = bm2/(bm1+bm2)
                                !get face value of u
                                up = au
                                um = au
                                !picking side for the weights
                                if (up.ge.0) then
                                   wpl = max(wplw,wprw)
                                   wpr = min(wplw,wprw)
                                else
                                   wpl = min(wplw,wprw)
                                   wpr = max(wplw,wprw)
                                end if
                                if (um.ge.0) then
                                   wml = max(wmlw,wmrw)
                                   wmr = min(wmlw,wmrw)
                                else
                                   wml = min(wmlw,wmrw)
                                   wmr = max(wmlw,wmrw)
                                end if
                                expAd = 0.
                                if (orderOneFix.eq.1) then
                                   ! if (augmentPlot.eq.1) then
                                   !    agu(uc,13)=0
                                   ! end if
                                   if (aum*aup .le. 0) then
                                      aumax = max(abs(au),abs(aup),abs(
     & aum))
                                      expAd = (1./12.)*aumax*( u(i1-2*
     & s1,i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(
     & i1,i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*
     & s2,i3+2*s3,var))
                                      wpl = 0.5
                                      wpr = 0.5
                                      wml = 0.5
                                      wmr = 0.5
                                      ! plot au at where fix is turned on
                                      ! if (augmentPlot.eq.1) then
                                      !    agu(uc,13)=au
                                      ! end if
                                   end if
                                end if
                                ! wpl = 0.5
                                ! wpr = 0.5
                                ! wml = 0.5
                                ! wmr = 0.5
                                Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,
     & var) + 5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
                                Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 
     & 5.*u(i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
                                Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*
     & s3,var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
                                Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,
     & var) + 5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
                                Fp  = wpl*Fpl + wpr*Fpr
                                Fm  = wml*Fml + wmr*Fmr
                                   agu(vard,var)= au*(Fp - Fm)/drl + 
     & expAd/drl
                              !     if ( i1 .eq. 14 .and. i2.eq.16 ) then
                                      !     if ( i1 .eq. 14) then
                                      ! write(*,'(" betapl,betapr,betaml,betamr = ",4F10.8)')  betapl,betapr,betaml,betamr
                                      ! write(*,'(" ap1,ap2,am1,am2 = ",4e10.8)')  ap1,ap2,am1,am2
                                      ! write(*,'(" wp1,wpr,wm1,wm2 = ",4F10.8)')  wp1,wp2,wm1,wm2
                                      ! write(*,'(" bp1,bp2,bm1,bm2 = ",4F10.8)')  bp1,bp2,bm1,bm2
                                      ! write(*,'(" wplw,wprw,wmlw,wmrw = ",4F10.8)')  wplw,wprw,wmlw,wmrw
                                      ! write(*,'(" wpl,wpr,wml,wmr = ",4F10.8)')  wpl,wpr,wml,wmr
                                      ! write(*,'(" agu = ",F10.8)') agu(vard,var)
                              !     end if
                              !u*vx
                              s1=1
                              s2=0
                              s3=0
                              var = vc
                              vard= uc
                              drl = dx(0)
                                Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*
     & u(i1   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
                                Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    
     & u(i1-s1,i2-s2,i3-s3,var)
                                Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) 
     & - 2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
                                Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) 
     & -    u(i1   ,i2   ,i3   ,var)
                                Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,
     & i2-s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
                                Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var)
                                Amr   = Apl
                                Bmr   = Bpl
                                betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 
     & 1./4.*Bpl**2
                                betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 
     & 1./4.*Bpr**2
                                betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 
     & 1./4.*Bml**2
                                betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 
     & 1./4.*Bmr**2
                                ep     = 1e-15
                                ap1    = (1./2.)/(ep + betapl)**(4)
                                ap2    = (1./2.)/(ep + betapr)**(4)
                                am1    = (1./2.)/(ep + betaml)**(4)
                                am2    = (1./2.)/(ep + betamr)**(4)
                                !calculate weights
                                wp1 = ap1/(ap1+ap2)
                                wp2 = ap2/(ap1+ap2)
                                wm1 = am1/(am1+am2)
                                wm2 = am2/(am1+am2)
                                !mapping of the weights
                                bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                                bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                                bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                                bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                                wplw = bp1/(bp1+bp2)
                                wprw = bp2/(bp1+bp2)
                                wmlw = bm1/(bm1+bm2)
                                wmrw = bm2/(bm1+bm2)
                                !get face value of u
                                up = au
                                um = au
                                !picking side for the weights
                                if (up.ge.0) then
                                   wpl = max(wplw,wprw)
                                   wpr = min(wplw,wprw)
                                else
                                   wpl = min(wplw,wprw)
                                   wpr = max(wplw,wprw)
                                end if
                                if (um.ge.0) then
                                   wml = max(wmlw,wmrw)
                                   wmr = min(wmlw,wmrw)
                                else
                                   wml = min(wmlw,wmrw)
                                   wmr = max(wmlw,wmrw)
                                end if
                                expAd = 0.
                                if (orderOneFix.eq.1) then
                                   ! if (augmentPlot.eq.1) then
                                   !    agu(uc,13)=0
                                   ! end if
                                   if (aum*aup .le. 0) then
                                      aumax = max(abs(au),abs(aup),abs(
     & aum))
                                      expAd = (1./12.)*aumax*( u(i1-2*
     & s1,i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(
     & i1,i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*
     & s2,i3+2*s3,var))
                                      wpl = 0.5
                                      wpr = 0.5
                                      wml = 0.5
                                      wmr = 0.5
                                      ! plot au at where fix is turned on
                                      ! if (augmentPlot.eq.1) then
                                      !    agu(uc,13)=au
                                      ! end if
                                   end if
                                end if
                                ! wpl = 0.5
                                ! wpr = 0.5
                                ! wml = 0.5
                                ! wmr = 0.5
                                Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,
     & var) + 5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
                                Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 
     & 5.*u(i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
                                Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*
     & s3,var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
                                Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,
     & var) + 5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
                                Fp  = wpl*Fpl + wpr*Fpr
                                Fm  = wml*Fml + wmr*Fmr
                                   agu(vard,var)= au*(Fp - Fm)/drl + 
     & expAd/drl
                              !     if ( i1 .eq. 14 .and. i2.eq.16 ) then
                                      !     if ( i1 .eq. 14) then
                                      ! write(*,'(" betapl,betapr,betaml,betamr = ",4F10.8)')  betapl,betapr,betaml,betamr
                                      ! write(*,'(" ap1,ap2,am1,am2 = ",4e10.8)')  ap1,ap2,am1,am2
                                      ! write(*,'(" wp1,wpr,wm1,wm2 = ",4F10.8)')  wp1,wp2,wm1,wm2
                                      ! write(*,'(" bp1,bp2,bm1,bm2 = ",4F10.8)')  bp1,bp2,bm1,bm2
                                      ! write(*,'(" wplw,wprw,wmlw,wmrw = ",4F10.8)')  wplw,wprw,wmlw,wmrw
                                      ! write(*,'(" wpl,wpr,wml,wmr = ",4F10.8)')  wpl,wpr,wml,wmr
                                      ! write(*,'(" agu = ",F10.8)') agu(vard,var)
                              !     end if
                              !v*uy
                              au = uu(i1,i2,i3,vc)
                              aum = u(i1,i2-1,i3,vc) - gvV
                              aup = u(i1,i2+1,i3,vc) - gvV
                              s1=0
                              s2=1
                              s3=0
                              var = uc
                              vard= vc
                              drl = dx(1)
                                Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*
     & u(i1   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
                                Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    
     & u(i1-s1,i2-s2,i3-s3,var)
                                Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) 
     & - 2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
                                Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) 
     & -    u(i1   ,i2   ,i3   ,var)
                                Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,
     & i2-s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
                                Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var)
                                Amr   = Apl
                                Bmr   = Bpl
                                betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 
     & 1./4.*Bpl**2
                                betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 
     & 1./4.*Bpr**2
                                betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 
     & 1./4.*Bml**2
                                betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 
     & 1./4.*Bmr**2
                                ep     = 1e-15
                                ap1    = (1./2.)/(ep + betapl)**(4)
                                ap2    = (1./2.)/(ep + betapr)**(4)
                                am1    = (1./2.)/(ep + betaml)**(4)
                                am2    = (1./2.)/(ep + betamr)**(4)
                                !calculate weights
                                wp1 = ap1/(ap1+ap2)
                                wp2 = ap2/(ap1+ap2)
                                wm1 = am1/(am1+am2)
                                wm2 = am2/(am1+am2)
                                !mapping of the weights
                                bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                                bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                                bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                                bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                                wplw = bp1/(bp1+bp2)
                                wprw = bp2/(bp1+bp2)
                                wmlw = bm1/(bm1+bm2)
                                wmrw = bm2/(bm1+bm2)
                                !get face value of u
                                up = au
                                um = au
                                !picking side for the weights
                                if (up.ge.0) then
                                   wpl = max(wplw,wprw)
                                   wpr = min(wplw,wprw)
                                else
                                   wpl = min(wplw,wprw)
                                   wpr = max(wplw,wprw)
                                end if
                                if (um.ge.0) then
                                   wml = max(wmlw,wmrw)
                                   wmr = min(wmlw,wmrw)
                                else
                                   wml = min(wmlw,wmrw)
                                   wmr = max(wmlw,wmrw)
                                end if
                                expAd = 0.
                                if (orderOneFix.eq.1) then
                                   ! if (augmentPlot.eq.1) then
                                   !    agu(uc,13)=0
                                   ! end if
                                   if (aum*aup .le. 0) then
                                      aumax = max(abs(au),abs(aup),abs(
     & aum))
                                      expAd = (1./12.)*aumax*( u(i1-2*
     & s1,i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(
     & i1,i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*
     & s2,i3+2*s3,var))
                                      wpl = 0.5
                                      wpr = 0.5
                                      wml = 0.5
                                      wmr = 0.5
                                      ! plot au at where fix is turned on
                                      ! if (augmentPlot.eq.1) then
                                      !    agu(uc,13)=au
                                      ! end if
                                   end if
                                end if
                                ! wpl = 0.5
                                ! wpr = 0.5
                                ! wml = 0.5
                                ! wmr = 0.5
                                Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,
     & var) + 5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
                                Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 
     & 5.*u(i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
                                Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*
     & s3,var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
                                Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,
     & var) + 5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
                                Fp  = wpl*Fpl + wpr*Fpr
                                Fm  = wml*Fml + wmr*Fmr
                                   agu(vard,var)= au*(Fp - Fm)/drl + 
     & expAd/drl
                              !     if ( i1 .eq. 14 .and. i2.eq.16 ) then
                                      !     if ( i1 .eq. 14) then
                                      ! write(*,'(" betapl,betapr,betaml,betamr = ",4F10.8)')  betapl,betapr,betaml,betamr
                                      ! write(*,'(" ap1,ap2,am1,am2 = ",4e10.8)')  ap1,ap2,am1,am2
                                      ! write(*,'(" wp1,wpr,wm1,wm2 = ",4F10.8)')  wp1,wp2,wm1,wm2
                                      ! write(*,'(" bp1,bp2,bm1,bm2 = ",4F10.8)')  bp1,bp2,bm1,bm2
                                      ! write(*,'(" wplw,wprw,wmlw,wmrw = ",4F10.8)')  wplw,wprw,wmlw,wmrw
                                      ! write(*,'(" wpl,wpr,wml,wmr = ",4F10.8)')  wpl,wpr,wml,wmr
                                      ! write(*,'(" agu = ",F10.8)') agu(vard,var)
                              !     end if
                              !v*vy
                              s1=0
                              s2=1
                              s3=0
                              var = vc
                              vard= vc
                              drl = dx(1)
                                Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*
     & u(i1   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
                                Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    
     & u(i1-s1,i2-s2,i3-s3,var)
                                Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) 
     & - 2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
                                Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) 
     & -    u(i1   ,i2   ,i3   ,var)
                                Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,
     & i2-s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
                                Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var)
                                Amr   = Apl
                                Bmr   = Bpl
                                betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 
     & 1./4.*Bpl**2
                                betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 
     & 1./4.*Bpr**2
                                betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 
     & 1./4.*Bml**2
                                betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 
     & 1./4.*Bmr**2
                                ep     = 1e-15
                                ap1    = (1./2.)/(ep + betapl)**(4)
                                ap2    = (1./2.)/(ep + betapr)**(4)
                                am1    = (1./2.)/(ep + betaml)**(4)
                                am2    = (1./2.)/(ep + betamr)**(4)
                                !calculate weights
                                wp1 = ap1/(ap1+ap2)
                                wp2 = ap2/(ap1+ap2)
                                wm1 = am1/(am1+am2)
                                wm2 = am2/(am1+am2)
                                !mapping of the weights
                                bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                                bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                                bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                                bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                                wplw = bp1/(bp1+bp2)
                                wprw = bp2/(bp1+bp2)
                                wmlw = bm1/(bm1+bm2)
                                wmrw = bm2/(bm1+bm2)
                                !get face value of u
                                up = au
                                um = au
                                !picking side for the weights
                                if (up.ge.0) then
                                   wpl = max(wplw,wprw)
                                   wpr = min(wplw,wprw)
                                else
                                   wpl = min(wplw,wprw)
                                   wpr = max(wplw,wprw)
                                end if
                                if (um.ge.0) then
                                   wml = max(wmlw,wmrw)
                                   wmr = min(wmlw,wmrw)
                                else
                                   wml = min(wmlw,wmrw)
                                   wmr = max(wmlw,wmrw)
                                end if
                                expAd = 0.
                                if (orderOneFix.eq.1) then
                                   ! if (augmentPlot.eq.1) then
                                   !    agu(uc,13)=0
                                   ! end if
                                   if (aum*aup .le. 0) then
                                      aumax = max(abs(au),abs(aup),abs(
     & aum))
                                      expAd = (1./12.)*aumax*( u(i1-2*
     & s1,i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(
     & i1,i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*
     & s2,i3+2*s3,var))
                                      wpl = 0.5
                                      wpr = 0.5
                                      wml = 0.5
                                      wmr = 0.5
                                      ! plot au at where fix is turned on
                                      ! if (augmentPlot.eq.1) then
                                      !    agu(uc,13)=au
                                      ! end if
                                   end if
                                end if
                                ! wpl = 0.5
                                ! wpr = 0.5
                                ! wml = 0.5
                                ! wmr = 0.5
                                Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,
     & var) + 5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
                                Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 
     & 5.*u(i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
                                Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*
     & s3,var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
                                Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,
     & var) + 5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
                                Fp  = wpl*Fpl + wpr*Fpr
                                Fm  = wml*Fml + wmr*Fmr
                                   agu(vard,var)= au*(Fp - Fm)/drl + 
     & expAd/drl
                              !     if ( i1 .eq. 14 .and. i2.eq.16 ) then
                                      !     if ( i1 .eq. 14) then
                                      ! write(*,'(" betapl,betapr,betaml,betamr = ",4F10.8)')  betapl,betapr,betaml,betamr
                                      ! write(*,'(" ap1,ap2,am1,am2 = ",4e10.8)')  ap1,ap2,am1,am2
                                      ! write(*,'(" wp1,wpr,wm1,wm2 = ",4F10.8)')  wp1,wp2,wm1,wm2
                                      ! write(*,'(" bp1,bp2,bm1,bm2 = ",4F10.8)')  bp1,bp2,bm1,bm2
                                      ! write(*,'(" wplw,wprw,wmlw,wmrw = ",4F10.8)')  wplw,wprw,wmlw,wmrw
                                      ! write(*,'(" wpl,wpr,wml,wmr = ",4F10.8)')  wpl,wpr,wml,wmr
                                      ! write(*,'(" agu = ",F10.8)') agu(vard,var)
                              !     end if
                             ! u*ux-------------------------------------------------------------------------
                             !calculate smooth factor
                             ! Apl   = u(i1+1,i2,i3,uc) - 2.*u(i1  ,i2,i3,uc)   + u(i1-1,i2,i3,uc)
                             ! Bpl   = u(i1+1,i2,i3,uc) -    u(i1-1,i2,i3,uc)
                             ! Apr   = u(i1+2,i2,i3,uc) - 2.*u(i1+1,i2,i3,uc) + u(i1,i2,i3,uc)
                             ! Bpr   = u(i1+2,i2,i3,uc) -    u(i1  ,i2,i3,uc)
                             ! Aml   = u(i1,i2,i3,uc) - 2.*u(i1-1,i2,i3,uc) + u(i1-2,i2,i3,uc)
                             ! Bml   = u(i1,i2,i3,uc) -    u(i1-2,i2,i3,uc)
                             ! Amr   = Apl
                             ! Bmr   = Bpl
                             ! betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 1./4.*Bpl**2
                             ! betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 1./4.*Bpr**2
                             ! betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 1./4.*Bml**2
                             ! betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 1./4.*Bmr**2
                             ! ep     = 1e-15
                             ! ap1    = (1./2.)/(ep + betapl)**(4)
                             ! ap2    = (1./2.)/(ep + betapr)**(4)
                             ! am1    = (1./2.)/(ep + betaml)**(4)
                             ! am2    = (1./2.)/(ep + betamr)**(4)
                             ! !calculate weights
                             ! wp1 = ap1/(ap1+ap2)
                             ! wp2 = ap2/(ap1+ap2)
                             ! wm1 = am1/(am1+am2)
                             ! wm2 = am2/(am1+am2)
                             ! !mapping of the weights
                             ! bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                             ! bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                             ! bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                             ! bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                             ! wplw = bp1/(bp1+bp2)
                             ! wprw = bp2/(bp1+bp2)
                             ! wmlw = bm1/(bm1+bm2)
                             ! wmrw = bm2/(bm1+bm2)
                             ! !get face value of u
                             ! up = au
                             ! um = au
                             ! !picking side for the weights
                             !  if (up.ge.0) then
                             !    wpl = max(wplw,wprw)
                             !    wpr = min(wplw,wprw)
                             ! else
                             !    wpl = min(wplw,wprw)
                             !    wpr = max(wplw,wprw)
                             ! end if
                             ! if (um.ge.0) then
                             !    wml = max(wmlw,wmrw)
                             !    wmr = min(wmlw,wmrw)
                             ! else
                             !    wml = min(wmlw,wmrw)
                             !    wmr = max(wmlw,wmrw)
                             ! end if
                             ! expAd = 0.
                             ! if (orderOneFix.eq.1) then 
                             !    if (augmentPlot.eq.1) then
                             !       agu(uc,13)=0
                             !    end if
                             !    if (aum*aup .le. 0) then
                             !       aumax = max(abs(au),abs(aup),abs(aum))
                             !       expAd = (1./12.)*aumax*( u(i1-2,i2,i3,uc) - 4.*u(i1-1,i2,i3,uc) + 6.*u(i1,i2,i3,uc) - 4.*u(i1+1,i2,i3,uc) +  u(i1+2,i2,i3,uc))
                             !       wpl = 0.5
                             !       wpr = 0.5
                             !       wml = 0.5
                             !       wmr = 0.5
                             !       ! plot au at where fix is turned on
                             !       if (augmentPlot.eq.1) then
                             !          agu(uc,13)=au
                             !       end if
                             !    end if
                             ! end if
                             ! wpl = 0.5
                             ! wpr = 0.5
                             ! wml = 0.5
                             ! wmr = 0.5
                             ! Fpl = 1./6.*(  -u(i1-1,i2,i3,uc) + 5.*u(i1,i2,i3,uc)   + 2.*u(i1+1,i2,i3,uc))
                             ! Fpr = 1./6.*( 2.*u(i1,i2,i3,uc)   + 5.*u(i1+1,i2,i3,uc) -   u(i1+2,i2,i3,uc))
                             ! Fml = 1./6.*(  -u(i1-2,i2,i3,uc) + 5.*u(i1-1,i2,i3,uc)   + 2.*u(i1,i2,i3,uc))
                             ! Fmr = 1./6.*( 2.*u(i1-1,i2,i3,uc) + 5.*u(i1,i2,i3,uc) -       u(i1+1,i2,i3,uc))
                             ! Fp  = wpl*Fpl + wpr*Fpr
                             ! Fm  = wml*Fml + wmr*Fmr
                             ! agu(uc,uc)= au*(Fp - Fm)/(dx(0)) + expAd*dx(0)**(-1)
                             ! if (augmentPlot.eq.1) then
                             !    agu(uc,3)=au
                             !    agu(uc,4)=wpl
                             !    agu(uc,5)=wpr
                             !    agu(uc,6)=wml
                             !    agu(uc,7)=wmr
                             ! end if
                             ! u*vx-------------------------------------------------------------------------
                             !calculate smooth factor
                             ! Apl   = u(i1+1,i2,i3,vc) - 2.*u(i1,i2,i3,vc)   + u(i1-1,i2,i3,vc)
                             ! Bpl   = u(i1+1,i2,i3,vc) -    u(i1-1,i2,i3,vc)
                             ! Apr   = u(i1+2,i2,i3,vc) - 2.*u(i1+1,i2,i3,vc) + u(i1,i2,i3,vc)
                             ! Bpr   = u(i1+2,i2,i3,vc) -    u(i1,i2,i3,vc)
                             ! Aml   = u(i1,i2,i3,vc) - 2.*u(i1-1,i2,i3,vc) + u(i1-2,i2,i3,vc)
                             ! Bml   = u(i1,i2,i3,vc) -    u(i1-2,i2,i3,vc)
                             ! Amr   = Apl
                             ! Bmr   = Bpl
                             ! betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 1./4.*Bpl**2
                             ! betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 1./4.*Bpr**2
                             ! betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 1./4.*Bml**2
                             ! betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 1./4.*Bmr**2
                             ! ep     = 1e-15
                             ! ap1    = (1./2.)/(ep + betapl)**(4)
                             ! ap2    = (1./2.)/(ep + betapr)**(4)
                             ! am1    = (1./2.)/(ep + betaml)**(4)
                             ! am2    = (1./2.)/(ep + betamr)**(4)
                             ! !calculate weights
                             ! wp1 = ap1/(ap1+ap2)
                             ! wp2 = ap2/(ap1+ap2)
                             ! wm1 = am1/(am1+am2)
                             ! wm2 = am2/(am1+am2)
                             ! !mapping of the weights
                             ! bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                             ! bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                             ! bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                             ! bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                             ! wplw = bp1/(bp1+bp2)
                             ! wprw = bp2/(bp1+bp2)
                             ! wmlw = bm1/(bm1+bm2)
                             ! wmrw = bm2/(bm1+bm2)
                             ! !get face value of u
                             ! ! up = 0.5*(u(i1  ,i2,i3,uc)+u(i1+1,i2,i3,uc))
                             ! ! um = 0.5*(u(i1-1,i2,i3,uc)+u(i1  ,i2,i3,uc))
                             ! up = au
                             ! um = au
                             ! !picking side for the weights
                             ! if (up.ge.0) then
                             !    wpl = max(wplw,wprw)
                             !    wpr = min(wplw,wprw)
                             ! else
                             !    wpl = min(wplw,wprw)
                             !    wpr = max(wplw,wprw)
                             ! end if
                             ! if (um.ge.0) then
                             !    wml = max(wmlw,wmrw)
                             !    wmr = min(wmlw,wmrw)
                             ! else
                             !    wml = min(wmlw,wmrw)
                             !    wmr = max(wmlw,wmrw)
                             ! end if
                             ! expAd = 0.
                             ! if (orderOneFix.eq.1) then
                             !    if (aum*aup .le. 0) then
                             !       aumax = max(abs(au),abs(aup),abs(aum))
                             !       expAd = (1./12.)*aumax*( u(i1-2,i2,i3,vc) - 4.*u(i1-1,i2,i3,vc) + 6.*u(i1,i2,i3,vc) - 4.*u(i1+1,i2,i3,vc) +  u(i1+2,i2,i3,vc))
                             !       wpl = 0.5
                             !       wpr = 0.5
                             !       wml = 0.5
                             !       wmr = 0.5
                             !    end if
                             ! end if
                             ! Fpl = 1./6.*(  -u(i1-1,i2,i3,vc) + 5.*u(i1,i2,i3,vc)   + 2.*u(i1+1,i2,i3,vc))
                             ! Fpr = 1./6.*( 2.*u(i1,i2,i3,vc)   + 5.*u(i1+1,i2,i3,vc) -   u(i1+2,i2,i3,vc))
                             ! Fml = 1./6.*(  -u(i1-2,i2,i3,vc) + 5.*u(i1-1,i2,i3,vc)   + 2.*u(i1,i2,i3,vc))
                             ! Fmr = 1./6.*( 2.*u(i1-1,i2,i3,vc) + 5.*u(i1,i2,i3,vc) -       u(i1+1,i2,i3,vc))
                             ! Fp  = wpl*Fpl + wpr*Fpr
                             ! Fm  = wml*Fml + wmr*Fmr
                             ! agu(uc,vc)= au*(Fp - Fm)/(dx(0)) + expAd*dx(0)**(-1) 
                             ! if (augmentPlot .eq. 1) then 
                             !    agu(vc,4)=wpl
                             !    agu(vc,5)=wpr
                             !    agu(vc,6)=wml
                             !    agu(vc,7)=wmr
                             ! end if
                             !v*uy-------------------------------------------------------------------------
                             ! au = u(i1,i2,i3,vc)
                             ! aum = u(i1,i2-1,i3,vc)
                             ! aup = u(i1,i2+1,i3,vc)
                             !calculate smooth factor
                             ! Apl   = u(i1,i2+1,i3,uc) - 2.*u(i1,i2,i3,uc)   + u(i1,i2-1,i3,uc)
                             ! Bpl   = u(i1,i2+1,i3,uc) -    u(i1,i2-1,i3,uc)
                             ! Apr   = u(i1,i2+2,i3,uc) - 2.*u(i1,i2+1,i3,uc) + u(i1,i2,i3,uc)
                             ! Bpr   = u(i1,i2+2,i3,uc) -    u(i1,i2,i3,uc)
                             ! Aml   = u(i1,i2,i3,uc) - 2.*u(i1,i2-1,i3,uc) + u(i1,i2-2,i3,uc)
                             ! Bml   = u(i1,i2,i3,uc) -    u(i1,i2-2,i3,uc)
                             ! Amr   = Apl
                             ! Bmr   = Bpl
                             ! betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 1./4.*Bpl**2
                             ! betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 1./4.*Bpr**2
                             ! betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 1./4.*Bml**2
                             ! betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 1./4.*Bmr**2
                             ! ep     = 1e-15
                             ! ap1    = (1./2.)/(ep + betapl)**(4)
                             ! ap2    = (1./2.)/(ep + betapr)**(4)
                             ! am1    = (1./2.)/(ep + betaml)**(4)
                             ! am2    = (1./2.)/(ep + betamr)**(4)
                             ! !calculate weights
                             ! wp1 = ap1/(ap1+ap2)
                             ! wp2 = ap2/(ap1+ap2)
                             ! wm1 = am1/(am1+am2)
                             ! wm2 = am2/(am1+am2)
                             ! !mapping of the weights
                             ! bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                             ! bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                             ! bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                             ! bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                             ! wplw = bp1/(bp1+bp2)
                             ! wprw = bp2/(bp1+bp2)
                             ! wmlw = bm1/(bm1+bm2)
                             ! wmrw = bm2/(bm1+bm2)
                             ! !get face value of v
                             ! ! up = 0.5*(u(i1  ,i2+1,i3,vc)+u(i1,i2 ,i3,vc))
                             ! ! um = 0.5*(u(i1  ,i2  ,i3,vc)+u(i1,i2-1,i3,vc))
                             ! up = au
                             ! um = au
                             ! !picking side for the weights
                             !  if (up.ge.0) then
                             !    wpl = max(wplw,wprw)
                             !    wpr = min(wplw,wprw)
                             ! else
                             !    wpl = min(wplw,wprw)
                             !    wpr = max(wplw,wprw)
                             ! end if
                             ! if (um.ge.0) then
                             !    wml = max(wmlw,wmrw)
                             !    wmr = min(wmlw,wmrw)
                             ! else
                             !    wml = min(wmlw,wmrw)
                             !    wmr = max(wmlw,wmrw)
                             ! end if
                             ! expAd = 0.
                             ! if (orderOneFix.eq.1) then 
                             !    if (augmentPlot .eq. 1) then 
                             !       agu(vc,13)=0
                             !    end if
                             !    if (aum*aup .le. 0) then
                             !       aumax = max(abs(au),abs(aup),abs(aum))
                             !       expAd = (1./12.)*aumax*( u(i1,i2-2,i3,uc) - 4.*u(i1,i2-1,i3,uc) + 6.*u(i1,i2,i3,uc) - 4.*u(i1,i2+1,i3,uc) +  u(i1,i2+2,i3,uc))
                             !       wpl = 0.5
                             !       wpr = 0.5
                             !       wml = 0.5
                             !       wmr = 0.5
                             !       if (augmentPlot .eq. 1) then 
                             !          agu(vc,13)=au
                             !       end if
                             !    end if
                             ! end if
                             ! Fpl = 1./6.*(  -u(i1,i2-1,i3,uc) + 5.*u(i1,i2,i3,uc)   + 2.*u(i1,i2+1,i3,uc))
                             ! Fpr = 1./6.*( 2.*u(i1,i2,i3,uc)   + 5.*u(i1,i2+1,i3,uc) -   u(i1,i2+2,i3,uc))
                             ! Fml = 1./6.*(  -u(i1,i2-2,i3,uc) + 5.*u(i1,i2-1,i3,uc)   + 2.*u(i1,i2,i3,uc))
                             ! Fmr = 1./6.*( 2.*u(i1,i2-1,i3,uc) + 5.*u(i1,i2,i3,uc) -       u(i1,i2+1,i3,uc))
                             ! Fp  = wpl*Fpl + wpr*Fpr
                             ! Fm  = wml*Fml + wmr*Fmr
                             ! agu(vc,uc)= au*(Fp - Fm)/(dx(1)) + expAd*dx(1)**(-1)
                             ! if (augmentPlot .eq. 1) then 
                             !    agu(uc,8)=au
                             !    agu(uc,9) =wpl
                             !    agu(uc,10)=wpr
                             !    agu(uc,11)=wml
                             !    agu(uc,12)=wmr
                             ! end if
                             !v*vy---------------------------------------------------------------
                             ! Apl   = u(i1,i2+1,i3,vc) - 2.*u(i1,i2,i3,vc)   + u(i1,i2-1,i3,vc)
                             ! Bpl   = u(i1,i2+1,i3,vc) -    u(i1,i2-1,i3,vc)
                             ! Apr   = u(i1,i2+2,i3,vc) - 2.*u(i1,i2+1,i3,vc) + u(i1,i2,i3,vc)
                             ! Bpr   = u(i1,i2+2,i3,vc) -    u(i1,i2,i3,vc)
                             ! Aml   = u(i1,i2,i3,vc) - 2.*u(i1,i2-1,i3,vc) + u(i1,i2-2,i3,vc)
                             ! Bml   = u(i1,i2,i3,vc) -    u(i1,i2-2,i3,vc)
                             ! Amr   = Apl
                             ! Bmr   = Bpl
                             ! betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 1./4.*Bpl**2
                             ! betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 1./4.*Bpr**2
                             ! betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 1./4.*Bml**2
                             ! betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 1./4.*Bmr**2
                             ! ep     = 1e-15
                             ! ap1    = (1./2.)/(ep + betapl)**(4)
                             ! ap2    = (1./2.)/(ep + betapr)**(4)
                             ! am1    = (1./2.)/(ep + betaml)**(4)
                             ! am2    = (1./2.)/(ep + betamr)**(4)
                             ! !calculate weights
                             ! wp1 = ap1/(ap1+ap2)
                             ! wp2 = ap2/(ap1+ap2)
                             ! wm1 = am1/(am1+am2)
                             ! wm2 = am2/(am1+am2)
                             ! !mapping of the weights
                             ! bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                             ! bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                             ! bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                             ! bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                             ! wplw = bp1/(bp1+bp2)
                             ! wprw = bp2/(bp1+bp2)
                             ! wmlw = bm1/(bm1+bm2)
                             ! wmrw = bm2/(bm1+bm2)
                             ! !get face value of v
                             !  ! up = 0.5*(u(i1  ,i2,i3,vc)  + u(i1,i2+1,i3,vc))
                             !  ! um = 0.5*(u(i1-1,i2,i3,vc)  + u(i1,i2  ,i3,vc))
                             ! up = au
                             ! um = au
                             ! !picking side for the weights
                             !  if (up.ge.0) then
                             !    wpl = max(wplw,wprw)
                             !    wpr = min(wplw,wprw)
                             ! else
                             !    wpl = min(wplw,wprw)
                             !    wpr = max(wplw,wprw)
                             ! end if
                             ! if (um.ge.0) then
                             !    wml = max(wmlw,wmrw)
                             !    wmr = min(wmlw,wmrw)
                             ! else
                             !    wml = min(wmlw,wmrw)
                             !    wmr = max(wmlw,wmrw)
                             ! end if
                             ! expAd=0.
                             ! if (orderOneFix.eq.1) then
                             !    if (aum*aup .le. 0) then
                             !       aumax = max(abs(au),abs(aup),abs(aum))
                             !       expAd = (1./12.)*aumax*( u(i1,i2-2,i3,vc) - 4.*u(i1,i2-1,i3,vc) + 6.*u(i1,i2,i3,vc) - 4.*u(i1,i2+1,i3,vc) +  u(i1,i2+2,i3,vc))
                             !       wpl = 0.5
                             !       wpr = 0.5
                             !       wml = 0.5
                             !       wmr = 0.5
                             !    end if
                             ! end if
                             ! Fpl = 1./6.*(  -u(i1,i2-1,i3,vc) + 5.*u(i1,i2,i3,vc)   + 2.*u(i1,i2+1,i3,vc))
                             ! Fpr = 1./6.*( 2.*u(i1,i2,i3,vc)   + 5.*u(i1,i2+1,i3,vc) -   u(i1,i2+2,i3,vc))
                             ! Fml = 1./6.*(  -u(i1,i2-2,i3,vc) + 5.*u(i1,i2-1,i3,vc)   + 2.*u(i1,i2,i3,vc))
                             ! Fmr = 1./6.*( 2.*u(i1,i2-1,i3,vc) + 5.*u(i1,i2,i3,vc) -       u(i1,i2+1,i3,vc))
                             ! Fp  = wpl*Fpl + wpr*Fpr
                             ! Fm  = wml*Fml + wmr*Fmr
                             ! agu(vc,vc)= au*(Fp - Fm)/(dx(1)) + expAd*dx(1)**(-1)
                             ! if (augmentPlot .eq. 1) then  
                             !    agu(vc,9) =wpl
                             !    agu(vc,10)=wpr
                             !    agu(vc,11)=wml
                             !    agu(vc,12)=wmr
                             ! end if
                          else ! if (upwindOrder.eq.4)
                             write(*,'(" getBewnoAdvection: only 4th 
     & order is avaliable now")' )
                             stop 777
                          end if
                      ! #If ("BWENO" == "CENTERED") 
                      !  ! -- centered advection ---
                      !  ! write(*,'(" getAdvection -- centered")')
                      !  #If "2" == "2"
                      !   ! -- 2D --
                      !   agu(uc,uc)=UU(uc)*UX(uc)
                      !   agu(vc,uc)=UU(vc)*UY(uc)
                      !   agu(uc,vc)=UU(uc)*UX(vc)
                      !   agu(vc,vc)=UU(vc)*UY(vc)
                      !  #Else
                      !   ! -- 3D -- *check me*
                      !   agu(uc,uc)=UU(uc)*UX(uc)
                      !   agu(vc,uc)=UU(vc)*UY(uc)
                      !   agu(wc,uc)=UU(wc)*UY(uc)
                      !   agu(uc,vc)=UU(uc)*UX(vc)
                      !   agu(vc,vc)=UU(vc)*UY(vc)
                      !   agu(wc,vc)=UU(wc)*UY(vc)
                      !   agu(uc,wc)=UU(uc)*UX(wc)
                      !   agu(vc,wc)=UU(vc)*UY(wc)
                      !   agu(wc,wc)=UU(wc)*UY(wc)
                      !  #End
                      ! #Elif "BWENO" == "BWENO" 
                      !   ! --- upwind scheme ---
                      !   ! for testing output this next message:
                      !   if( t.le. 0. )then
                      !     write(*,'(" getAdvection upwind scheme (7)")') 
                      !   end if
                      !   getUpwindAdvection(u,i1,i2,i3,NONE,2,4,rectangular, agu)
                      ! #Elif "BWENO" == "BWENO" 
                      !   ! --- Bweno scheme ---
                      !   ! for testing output this message:
                      !   if( t.le. 0. )then
                      !      write(*,'(" getAdvection BWENO scheme (7)")') 
                      !   end if
                      !   getBwenoAdvection(u,i1,i2,i3,NONE,2,4,rectangular, agu)
                      ! #Else
                      !   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
                      !   stop 999
                      ! #End
                        ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-
     & ux42r(i1,i2,i3,pc)
                        ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-
     & uy42r(i1,i2,i3,pc)
                 end if
                end do
                end do
                end do
               else
                do i3=n3a,n3b
                do i2=n2a,n2b
                do i1=n1a,n1b
                   ! explicit terms only, no diffusion
                     ! upwind approximation/bweno
                       ! --- Bweno scheme ---
                       ! for testing output this message:
                       ! if( t.le. 0. )then
                       !    write(*,'(" getAdvection BWENO scheme (7)")') 
                       ! end if
                          gvU = u(i1,i2,i3,uc) - uu(i1,i2,i3,uc)
                          gvV = u(i1,i2,i3,vc) - uu(i1,i2,i3,vc)
                          gvW = u(i1,i2,i3,wc) - uu(i1,i2,i3,wc)
                          orderOneFix = 1
                          augmentPlot = 1
                          if( upwindOrder.eq.4 )then
                         ! -- fourth order bweno --
                             au = uu(i1,i2,i3,uc)
                             aum = u(i1-1,i2,i3,uc) - gvU
                             aup = u(i1+1,i2,i3,uc) - gvU
                             !usage
                             !s1   is the shifted grid in r1 direction
                             !var  is the variable taking derivative (uc,vc,wc)
                             !vard is the variable of direction
                             !agu is only passed in for augmented plot
                             !expAd 
                             ! u*ux
                             s1=1
                             s2=0
                             s3=0
                             var = uc
                             vard= uc
                             drl = dx(0)
                               Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*u(
     & i1   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
                               Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    u(
     & i1-s1,i2-s2,i3-s3,var)
                               Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) -
     &  2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
                               Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) -
     &     u(i1   ,i2   ,i3   ,var)
                               Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,
     & i2-s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
                               Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var)
                               Amr   = Apl
                               Bmr   = Bpl
                               betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 
     & 1./4.*Bpl**2
                               betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 
     & 1./4.*Bpr**2
                               betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 
     & 1./4.*Bml**2
                               betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 
     & 1./4.*Bmr**2
                               ep     = 1e-15
                               ap1    = (1./2.)/(ep + betapl)**(4)
                               ap2    = (1./2.)/(ep + betapr)**(4)
                               am1    = (1./2.)/(ep + betaml)**(4)
                               am2    = (1./2.)/(ep + betamr)**(4)
                               !calculate weights
                               wp1 = ap1/(ap1+ap2)
                               wp2 = ap2/(ap1+ap2)
                               wm1 = am1/(am1+am2)
                               wm2 = am2/(am1+am2)
                               !mapping of the weights
                               bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                               bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                               bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                               bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                               wplw = bp1/(bp1+bp2)
                               wprw = bp2/(bp1+bp2)
                               wmlw = bm1/(bm1+bm2)
                               wmrw = bm2/(bm1+bm2)
                               !get face value of u
                               up = au
                               um = au
                               !picking side for the weights
                               if (up.ge.0) then
                                  wpl = max(wplw,wprw)
                                  wpr = min(wplw,wprw)
                               else
                                  wpl = min(wplw,wprw)
                                  wpr = max(wplw,wprw)
                               end if
                               if (um.ge.0) then
                                  wml = max(wmlw,wmrw)
                                  wmr = min(wmlw,wmrw)
                               else
                                  wml = min(wmlw,wmrw)
                                  wmr = max(wmlw,wmrw)
                               end if
                               expAd = 0.
                               if (orderOneFix.eq.1) then
                                  ! if (augmentPlot.eq.1) then
                                  !    agu(uc,13)=0
                                  ! end if
                                  if (aum*aup .le. 0) then
                                     aumax = max(abs(au),abs(aup),abs(
     & aum))
                                     expAd = (1./12.)*aumax*( u(i1-2*
     & s1,i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(
     & i1,i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*
     & s2,i3+2*s3,var))
                                     wpl = 0.5
                                     wpr = 0.5
                                     wml = 0.5
                                     wmr = 0.5
                                     ! plot au at where fix is turned on
                                     ! if (augmentPlot.eq.1) then
                                     !    agu(uc,13)=au
                                     ! end if
                                  end if
                               end if
                               ! wpl = 0.5
                               ! wpr = 0.5
                               ! wml = 0.5
                               ! wmr = 0.5
                               Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,var)
     &  + 5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
                               Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 5.*
     & u(i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
                               Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*
     & s3,var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
                               Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,
     & var) + 5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
                               Fp  = wpl*Fpl + wpr*Fpr
                               Fm  = wml*Fml + wmr*Fmr
                                  agu(vard,var)= au*(Fp - Fm)/drl + 
     & expAd/drl
                             !     if ( i1 .eq. 14 .and. i2.eq.16 ) then
                                     !     if ( i1 .eq. 14) then
                                     ! write(*,'(" betapl,betapr,betaml,betamr = ",4F10.8)')  betapl,betapr,betaml,betamr
                                     ! write(*,'(" ap1,ap2,am1,am2 = ",4e10.8)')  ap1,ap2,am1,am2
                                     ! write(*,'(" wp1,wpr,wm1,wm2 = ",4F10.8)')  wp1,wp2,wm1,wm2
                                     ! write(*,'(" bp1,bp2,bm1,bm2 = ",4F10.8)')  bp1,bp2,bm1,bm2
                                     ! write(*,'(" wplw,wprw,wmlw,wmrw = ",4F10.8)')  wplw,wprw,wmlw,wmrw
                                     ! write(*,'(" wpl,wpr,wml,wmr = ",4F10.8)')  wpl,wpr,wml,wmr
                                     ! write(*,'(" agu = ",F10.8)') agu(vard,var)
                             !     end if
                             !u*vx
                             s1=1
                             s2=0
                             s3=0
                             var = vc
                             vard= uc
                             drl = dx(0)
                               Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*u(
     & i1   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
                               Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    u(
     & i1-s1,i2-s2,i3-s3,var)
                               Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) -
     &  2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
                               Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) -
     &     u(i1   ,i2   ,i3   ,var)
                               Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,
     & i2-s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
                               Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var)
                               Amr   = Apl
                               Bmr   = Bpl
                               betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 
     & 1./4.*Bpl**2
                               betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 
     & 1./4.*Bpr**2
                               betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 
     & 1./4.*Bml**2
                               betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 
     & 1./4.*Bmr**2
                               ep     = 1e-15
                               ap1    = (1./2.)/(ep + betapl)**(4)
                               ap2    = (1./2.)/(ep + betapr)**(4)
                               am1    = (1./2.)/(ep + betaml)**(4)
                               am2    = (1./2.)/(ep + betamr)**(4)
                               !calculate weights
                               wp1 = ap1/(ap1+ap2)
                               wp2 = ap2/(ap1+ap2)
                               wm1 = am1/(am1+am2)
                               wm2 = am2/(am1+am2)
                               !mapping of the weights
                               bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                               bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                               bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                               bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                               wplw = bp1/(bp1+bp2)
                               wprw = bp2/(bp1+bp2)
                               wmlw = bm1/(bm1+bm2)
                               wmrw = bm2/(bm1+bm2)
                               !get face value of u
                               up = au
                               um = au
                               !picking side for the weights
                               if (up.ge.0) then
                                  wpl = max(wplw,wprw)
                                  wpr = min(wplw,wprw)
                               else
                                  wpl = min(wplw,wprw)
                                  wpr = max(wplw,wprw)
                               end if
                               if (um.ge.0) then
                                  wml = max(wmlw,wmrw)
                                  wmr = min(wmlw,wmrw)
                               else
                                  wml = min(wmlw,wmrw)
                                  wmr = max(wmlw,wmrw)
                               end if
                               expAd = 0.
                               if (orderOneFix.eq.1) then
                                  ! if (augmentPlot.eq.1) then
                                  !    agu(uc,13)=0
                                  ! end if
                                  if (aum*aup .le. 0) then
                                     aumax = max(abs(au),abs(aup),abs(
     & aum))
                                     expAd = (1./12.)*aumax*( u(i1-2*
     & s1,i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(
     & i1,i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*
     & s2,i3+2*s3,var))
                                     wpl = 0.5
                                     wpr = 0.5
                                     wml = 0.5
                                     wmr = 0.5
                                     ! plot au at where fix is turned on
                                     ! if (augmentPlot.eq.1) then
                                     !    agu(uc,13)=au
                                     ! end if
                                  end if
                               end if
                               ! wpl = 0.5
                               ! wpr = 0.5
                               ! wml = 0.5
                               ! wmr = 0.5
                               Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,var)
     &  + 5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
                               Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 5.*
     & u(i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
                               Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*
     & s3,var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
                               Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,
     & var) + 5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
                               Fp  = wpl*Fpl + wpr*Fpr
                               Fm  = wml*Fml + wmr*Fmr
                                  agu(vard,var)= au*(Fp - Fm)/drl + 
     & expAd/drl
                             !     if ( i1 .eq. 14 .and. i2.eq.16 ) then
                                     !     if ( i1 .eq. 14) then
                                     ! write(*,'(" betapl,betapr,betaml,betamr = ",4F10.8)')  betapl,betapr,betaml,betamr
                                     ! write(*,'(" ap1,ap2,am1,am2 = ",4e10.8)')  ap1,ap2,am1,am2
                                     ! write(*,'(" wp1,wpr,wm1,wm2 = ",4F10.8)')  wp1,wp2,wm1,wm2
                                     ! write(*,'(" bp1,bp2,bm1,bm2 = ",4F10.8)')  bp1,bp2,bm1,bm2
                                     ! write(*,'(" wplw,wprw,wmlw,wmrw = ",4F10.8)')  wplw,wprw,wmlw,wmrw
                                     ! write(*,'(" wpl,wpr,wml,wmr = ",4F10.8)')  wpl,wpr,wml,wmr
                                     ! write(*,'(" agu = ",F10.8)') agu(vard,var)
                             !     end if
                             !v*uy
                             au = uu(i1,i2,i3,vc)
                             aum = u(i1,i2-1,i3,vc) - gvV
                             aup = u(i1,i2+1,i3,vc) - gvV
                             s1=0
                             s2=1
                             s3=0
                             var = uc
                             vard= vc
                             drl = dx(1)
                               Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*u(
     & i1   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
                               Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    u(
     & i1-s1,i2-s2,i3-s3,var)
                               Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) -
     &  2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
                               Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) -
     &     u(i1   ,i2   ,i3   ,var)
                               Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,
     & i2-s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
                               Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var)
                               Amr   = Apl
                               Bmr   = Bpl
                               betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 
     & 1./4.*Bpl**2
                               betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 
     & 1./4.*Bpr**2
                               betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 
     & 1./4.*Bml**2
                               betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 
     & 1./4.*Bmr**2
                               ep     = 1e-15
                               ap1    = (1./2.)/(ep + betapl)**(4)
                               ap2    = (1./2.)/(ep + betapr)**(4)
                               am1    = (1./2.)/(ep + betaml)**(4)
                               am2    = (1./2.)/(ep + betamr)**(4)
                               !calculate weights
                               wp1 = ap1/(ap1+ap2)
                               wp2 = ap2/(ap1+ap2)
                               wm1 = am1/(am1+am2)
                               wm2 = am2/(am1+am2)
                               !mapping of the weights
                               bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                               bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                               bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                               bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                               wplw = bp1/(bp1+bp2)
                               wprw = bp2/(bp1+bp2)
                               wmlw = bm1/(bm1+bm2)
                               wmrw = bm2/(bm1+bm2)
                               !get face value of u
                               up = au
                               um = au
                               !picking side for the weights
                               if (up.ge.0) then
                                  wpl = max(wplw,wprw)
                                  wpr = min(wplw,wprw)
                               else
                                  wpl = min(wplw,wprw)
                                  wpr = max(wplw,wprw)
                               end if
                               if (um.ge.0) then
                                  wml = max(wmlw,wmrw)
                                  wmr = min(wmlw,wmrw)
                               else
                                  wml = min(wmlw,wmrw)
                                  wmr = max(wmlw,wmrw)
                               end if
                               expAd = 0.
                               if (orderOneFix.eq.1) then
                                  ! if (augmentPlot.eq.1) then
                                  !    agu(uc,13)=0
                                  ! end if
                                  if (aum*aup .le. 0) then
                                     aumax = max(abs(au),abs(aup),abs(
     & aum))
                                     expAd = (1./12.)*aumax*( u(i1-2*
     & s1,i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(
     & i1,i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*
     & s2,i3+2*s3,var))
                                     wpl = 0.5
                                     wpr = 0.5
                                     wml = 0.5
                                     wmr = 0.5
                                     ! plot au at where fix is turned on
                                     ! if (augmentPlot.eq.1) then
                                     !    agu(uc,13)=au
                                     ! end if
                                  end if
                               end if
                               ! wpl = 0.5
                               ! wpr = 0.5
                               ! wml = 0.5
                               ! wmr = 0.5
                               Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,var)
     &  + 5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
                               Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 5.*
     & u(i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
                               Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*
     & s3,var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
                               Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,
     & var) + 5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
                               Fp  = wpl*Fpl + wpr*Fpr
                               Fm  = wml*Fml + wmr*Fmr
                                  agu(vard,var)= au*(Fp - Fm)/drl + 
     & expAd/drl
                             !     if ( i1 .eq. 14 .and. i2.eq.16 ) then
                                     !     if ( i1 .eq. 14) then
                                     ! write(*,'(" betapl,betapr,betaml,betamr = ",4F10.8)')  betapl,betapr,betaml,betamr
                                     ! write(*,'(" ap1,ap2,am1,am2 = ",4e10.8)')  ap1,ap2,am1,am2
                                     ! write(*,'(" wp1,wpr,wm1,wm2 = ",4F10.8)')  wp1,wp2,wm1,wm2
                                     ! write(*,'(" bp1,bp2,bm1,bm2 = ",4F10.8)')  bp1,bp2,bm1,bm2
                                     ! write(*,'(" wplw,wprw,wmlw,wmrw = ",4F10.8)')  wplw,wprw,wmlw,wmrw
                                     ! write(*,'(" wpl,wpr,wml,wmr = ",4F10.8)')  wpl,wpr,wml,wmr
                                     ! write(*,'(" agu = ",F10.8)') agu(vard,var)
                             !     end if
                             !v*vy
                             s1=0
                             s2=1
                             s3=0
                             var = vc
                             vard= vc
                             drl = dx(1)
                               Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*u(
     & i1   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
                               Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    u(
     & i1-s1,i2-s2,i3-s3,var)
                               Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) -
     &  2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
                               Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) -
     &     u(i1   ,i2   ,i3   ,var)
                               Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,
     & i2-s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
                               Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var)
                               Amr   = Apl
                               Bmr   = Bpl
                               betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 
     & 1./4.*Bpl**2
                               betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 
     & 1./4.*Bpr**2
                               betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 
     & 1./4.*Bml**2
                               betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 
     & 1./4.*Bmr**2
                               ep     = 1e-15
                               ap1    = (1./2.)/(ep + betapl)**(4)
                               ap2    = (1./2.)/(ep + betapr)**(4)
                               am1    = (1./2.)/(ep + betaml)**(4)
                               am2    = (1./2.)/(ep + betamr)**(4)
                               !calculate weights
                               wp1 = ap1/(ap1+ap2)
                               wp2 = ap2/(ap1+ap2)
                               wm1 = am1/(am1+am2)
                               wm2 = am2/(am1+am2)
                               !mapping of the weights
                               bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                               bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                               bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                               bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                               wplw = bp1/(bp1+bp2)
                               wprw = bp2/(bp1+bp2)
                               wmlw = bm1/(bm1+bm2)
                               wmrw = bm2/(bm1+bm2)
                               !get face value of u
                               up = au
                               um = au
                               !picking side for the weights
                               if (up.ge.0) then
                                  wpl = max(wplw,wprw)
                                  wpr = min(wplw,wprw)
                               else
                                  wpl = min(wplw,wprw)
                                  wpr = max(wplw,wprw)
                               end if
                               if (um.ge.0) then
                                  wml = max(wmlw,wmrw)
                                  wmr = min(wmlw,wmrw)
                               else
                                  wml = min(wmlw,wmrw)
                                  wmr = max(wmlw,wmrw)
                               end if
                               expAd = 0.
                               if (orderOneFix.eq.1) then
                                  ! if (augmentPlot.eq.1) then
                                  !    agu(uc,13)=0
                                  ! end if
                                  if (aum*aup .le. 0) then
                                     aumax = max(abs(au),abs(aup),abs(
     & aum))
                                     expAd = (1./12.)*aumax*( u(i1-2*
     & s1,i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(
     & i1,i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*
     & s2,i3+2*s3,var))
                                     wpl = 0.5
                                     wpr = 0.5
                                     wml = 0.5
                                     wmr = 0.5
                                     ! plot au at where fix is turned on
                                     ! if (augmentPlot.eq.1) then
                                     !    agu(uc,13)=au
                                     ! end if
                                  end if
                               end if
                               ! wpl = 0.5
                               ! wpr = 0.5
                               ! wml = 0.5
                               ! wmr = 0.5
                               Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,var)
     &  + 5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
                               Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 5.*
     & u(i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
                               Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*
     & s3,var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
                               Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,
     & var) + 5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
                               Fp  = wpl*Fpl + wpr*Fpr
                               Fm  = wml*Fml + wmr*Fmr
                                  agu(vard,var)= au*(Fp - Fm)/drl + 
     & expAd/drl
                             !     if ( i1 .eq. 14 .and. i2.eq.16 ) then
                                     !     if ( i1 .eq. 14) then
                                     ! write(*,'(" betapl,betapr,betaml,betamr = ",4F10.8)')  betapl,betapr,betaml,betamr
                                     ! write(*,'(" ap1,ap2,am1,am2 = ",4e10.8)')  ap1,ap2,am1,am2
                                     ! write(*,'(" wp1,wpr,wm1,wm2 = ",4F10.8)')  wp1,wp2,wm1,wm2
                                     ! write(*,'(" bp1,bp2,bm1,bm2 = ",4F10.8)')  bp1,bp2,bm1,bm2
                                     ! write(*,'(" wplw,wprw,wmlw,wmrw = ",4F10.8)')  wplw,wprw,wmlw,wmrw
                                     ! write(*,'(" wpl,wpr,wml,wmr = ",4F10.8)')  wpl,wpr,wml,wmr
                                     ! write(*,'(" agu = ",F10.8)') agu(vard,var)
                             !     end if
                            ! u*ux-------------------------------------------------------------------------
                            !calculate smooth factor
                            ! Apl   = u(i1+1,i2,i3,uc) - 2.*u(i1  ,i2,i3,uc)   + u(i1-1,i2,i3,uc)
                            ! Bpl   = u(i1+1,i2,i3,uc) -    u(i1-1,i2,i3,uc)
                            ! Apr   = u(i1+2,i2,i3,uc) - 2.*u(i1+1,i2,i3,uc) + u(i1,i2,i3,uc)
                            ! Bpr   = u(i1+2,i2,i3,uc) -    u(i1  ,i2,i3,uc)
                            ! Aml   = u(i1,i2,i3,uc) - 2.*u(i1-1,i2,i3,uc) + u(i1-2,i2,i3,uc)
                            ! Bml   = u(i1,i2,i3,uc) -    u(i1-2,i2,i3,uc)
                            ! Amr   = Apl
                            ! Bmr   = Bpl
                            ! betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 1./4.*Bpl**2
                            ! betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 1./4.*Bpr**2
                            ! betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 1./4.*Bml**2
                            ! betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 1./4.*Bmr**2
                            ! ep     = 1e-15
                            ! ap1    = (1./2.)/(ep + betapl)**(4)
                            ! ap2    = (1./2.)/(ep + betapr)**(4)
                            ! am1    = (1./2.)/(ep + betaml)**(4)
                            ! am2    = (1./2.)/(ep + betamr)**(4)
                            ! !calculate weights
                            ! wp1 = ap1/(ap1+ap2)
                            ! wp2 = ap2/(ap1+ap2)
                            ! wm1 = am1/(am1+am2)
                            ! wm2 = am2/(am1+am2)
                            ! !mapping of the weights
                            ! bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                            ! bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                            ! bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                            ! bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                            ! wplw = bp1/(bp1+bp2)
                            ! wprw = bp2/(bp1+bp2)
                            ! wmlw = bm1/(bm1+bm2)
                            ! wmrw = bm2/(bm1+bm2)
                            ! !get face value of u
                            ! up = au
                            ! um = au
                            ! !picking side for the weights
                            !  if (up.ge.0) then
                            !    wpl = max(wplw,wprw)
                            !    wpr = min(wplw,wprw)
                            ! else
                            !    wpl = min(wplw,wprw)
                            !    wpr = max(wplw,wprw)
                            ! end if
                            ! if (um.ge.0) then
                            !    wml = max(wmlw,wmrw)
                            !    wmr = min(wmlw,wmrw)
                            ! else
                            !    wml = min(wmlw,wmrw)
                            !    wmr = max(wmlw,wmrw)
                            ! end if
                            ! expAd = 0.
                            ! if (orderOneFix.eq.1) then 
                            !    if (augmentPlot.eq.1) then
                            !       agu(uc,13)=0
                            !    end if
                            !    if (aum*aup .le. 0) then
                            !       aumax = max(abs(au),abs(aup),abs(aum))
                            !       expAd = (1./12.)*aumax*( u(i1-2,i2,i3,uc) - 4.*u(i1-1,i2,i3,uc) + 6.*u(i1,i2,i3,uc) - 4.*u(i1+1,i2,i3,uc) +  u(i1+2,i2,i3,uc))
                            !       wpl = 0.5
                            !       wpr = 0.5
                            !       wml = 0.5
                            !       wmr = 0.5
                            !       ! plot au at where fix is turned on
                            !       if (augmentPlot.eq.1) then
                            !          agu(uc,13)=au
                            !       end if
                            !    end if
                            ! end if
                            ! wpl = 0.5
                            ! wpr = 0.5
                            ! wml = 0.5
                            ! wmr = 0.5
                            ! Fpl = 1./6.*(  -u(i1-1,i2,i3,uc) + 5.*u(i1,i2,i3,uc)   + 2.*u(i1+1,i2,i3,uc))
                            ! Fpr = 1./6.*( 2.*u(i1,i2,i3,uc)   + 5.*u(i1+1,i2,i3,uc) -   u(i1+2,i2,i3,uc))
                            ! Fml = 1./6.*(  -u(i1-2,i2,i3,uc) + 5.*u(i1-1,i2,i3,uc)   + 2.*u(i1,i2,i3,uc))
                            ! Fmr = 1./6.*( 2.*u(i1-1,i2,i3,uc) + 5.*u(i1,i2,i3,uc) -       u(i1+1,i2,i3,uc))
                            ! Fp  = wpl*Fpl + wpr*Fpr
                            ! Fm  = wml*Fml + wmr*Fmr
                            ! agu(uc,uc)= au*(Fp - Fm)/(dx(0)) + expAd*dx(0)**(-1)
                            ! if (augmentPlot.eq.1) then
                            !    agu(uc,3)=au
                            !    agu(uc,4)=wpl
                            !    agu(uc,5)=wpr
                            !    agu(uc,6)=wml
                            !    agu(uc,7)=wmr
                            ! end if
                            ! u*vx-------------------------------------------------------------------------
                            !calculate smooth factor
                            ! Apl   = u(i1+1,i2,i3,vc) - 2.*u(i1,i2,i3,vc)   + u(i1-1,i2,i3,vc)
                            ! Bpl   = u(i1+1,i2,i3,vc) -    u(i1-1,i2,i3,vc)
                            ! Apr   = u(i1+2,i2,i3,vc) - 2.*u(i1+1,i2,i3,vc) + u(i1,i2,i3,vc)
                            ! Bpr   = u(i1+2,i2,i3,vc) -    u(i1,i2,i3,vc)
                            ! Aml   = u(i1,i2,i3,vc) - 2.*u(i1-1,i2,i3,vc) + u(i1-2,i2,i3,vc)
                            ! Bml   = u(i1,i2,i3,vc) -    u(i1-2,i2,i3,vc)
                            ! Amr   = Apl
                            ! Bmr   = Bpl
                            ! betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 1./4.*Bpl**2
                            ! betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 1./4.*Bpr**2
                            ! betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 1./4.*Bml**2
                            ! betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 1./4.*Bmr**2
                            ! ep     = 1e-15
                            ! ap1    = (1./2.)/(ep + betapl)**(4)
                            ! ap2    = (1./2.)/(ep + betapr)**(4)
                            ! am1    = (1./2.)/(ep + betaml)**(4)
                            ! am2    = (1./2.)/(ep + betamr)**(4)
                            ! !calculate weights
                            ! wp1 = ap1/(ap1+ap2)
                            ! wp2 = ap2/(ap1+ap2)
                            ! wm1 = am1/(am1+am2)
                            ! wm2 = am2/(am1+am2)
                            ! !mapping of the weights
                            ! bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                            ! bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                            ! bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                            ! bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                            ! wplw = bp1/(bp1+bp2)
                            ! wprw = bp2/(bp1+bp2)
                            ! wmlw = bm1/(bm1+bm2)
                            ! wmrw = bm2/(bm1+bm2)
                            ! !get face value of u
                            ! ! up = 0.5*(u(i1  ,i2,i3,uc)+u(i1+1,i2,i3,uc))
                            ! ! um = 0.5*(u(i1-1,i2,i3,uc)+u(i1  ,i2,i3,uc))
                            ! up = au
                            ! um = au
                            ! !picking side for the weights
                            ! if (up.ge.0) then
                            !    wpl = max(wplw,wprw)
                            !    wpr = min(wplw,wprw)
                            ! else
                            !    wpl = min(wplw,wprw)
                            !    wpr = max(wplw,wprw)
                            ! end if
                            ! if (um.ge.0) then
                            !    wml = max(wmlw,wmrw)
                            !    wmr = min(wmlw,wmrw)
                            ! else
                            !    wml = min(wmlw,wmrw)
                            !    wmr = max(wmlw,wmrw)
                            ! end if
                            ! expAd = 0.
                            ! if (orderOneFix.eq.1) then
                            !    if (aum*aup .le. 0) then
                            !       aumax = max(abs(au),abs(aup),abs(aum))
                            !       expAd = (1./12.)*aumax*( u(i1-2,i2,i3,vc) - 4.*u(i1-1,i2,i3,vc) + 6.*u(i1,i2,i3,vc) - 4.*u(i1+1,i2,i3,vc) +  u(i1+2,i2,i3,vc))
                            !       wpl = 0.5
                            !       wpr = 0.5
                            !       wml = 0.5
                            !       wmr = 0.5
                            !    end if
                            ! end if
                            ! Fpl = 1./6.*(  -u(i1-1,i2,i3,vc) + 5.*u(i1,i2,i3,vc)   + 2.*u(i1+1,i2,i3,vc))
                            ! Fpr = 1./6.*( 2.*u(i1,i2,i3,vc)   + 5.*u(i1+1,i2,i3,vc) -   u(i1+2,i2,i3,vc))
                            ! Fml = 1./6.*(  -u(i1-2,i2,i3,vc) + 5.*u(i1-1,i2,i3,vc)   + 2.*u(i1,i2,i3,vc))
                            ! Fmr = 1./6.*( 2.*u(i1-1,i2,i3,vc) + 5.*u(i1,i2,i3,vc) -       u(i1+1,i2,i3,vc))
                            ! Fp  = wpl*Fpl + wpr*Fpr
                            ! Fm  = wml*Fml + wmr*Fmr
                            ! agu(uc,vc)= au*(Fp - Fm)/(dx(0)) + expAd*dx(0)**(-1) 
                            ! if (augmentPlot .eq. 1) then 
                            !    agu(vc,4)=wpl
                            !    agu(vc,5)=wpr
                            !    agu(vc,6)=wml
                            !    agu(vc,7)=wmr
                            ! end if
                            !v*uy-------------------------------------------------------------------------
                            ! au = u(i1,i2,i3,vc)
                            ! aum = u(i1,i2-1,i3,vc)
                            ! aup = u(i1,i2+1,i3,vc)
                            !calculate smooth factor
                            ! Apl   = u(i1,i2+1,i3,uc) - 2.*u(i1,i2,i3,uc)   + u(i1,i2-1,i3,uc)
                            ! Bpl   = u(i1,i2+1,i3,uc) -    u(i1,i2-1,i3,uc)
                            ! Apr   = u(i1,i2+2,i3,uc) - 2.*u(i1,i2+1,i3,uc) + u(i1,i2,i3,uc)
                            ! Bpr   = u(i1,i2+2,i3,uc) -    u(i1,i2,i3,uc)
                            ! Aml   = u(i1,i2,i3,uc) - 2.*u(i1,i2-1,i3,uc) + u(i1,i2-2,i3,uc)
                            ! Bml   = u(i1,i2,i3,uc) -    u(i1,i2-2,i3,uc)
                            ! Amr   = Apl
                            ! Bmr   = Bpl
                            ! betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 1./4.*Bpl**2
                            ! betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 1./4.*Bpr**2
                            ! betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 1./4.*Bml**2
                            ! betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 1./4.*Bmr**2
                            ! ep     = 1e-15
                            ! ap1    = (1./2.)/(ep + betapl)**(4)
                            ! ap2    = (1./2.)/(ep + betapr)**(4)
                            ! am1    = (1./2.)/(ep + betaml)**(4)
                            ! am2    = (1./2.)/(ep + betamr)**(4)
                            ! !calculate weights
                            ! wp1 = ap1/(ap1+ap2)
                            ! wp2 = ap2/(ap1+ap2)
                            ! wm1 = am1/(am1+am2)
                            ! wm2 = am2/(am1+am2)
                            ! !mapping of the weights
                            ! bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                            ! bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                            ! bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                            ! bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                            ! wplw = bp1/(bp1+bp2)
                            ! wprw = bp2/(bp1+bp2)
                            ! wmlw = bm1/(bm1+bm2)
                            ! wmrw = bm2/(bm1+bm2)
                            ! !get face value of v
                            ! ! up = 0.5*(u(i1  ,i2+1,i3,vc)+u(i1,i2 ,i3,vc))
                            ! ! um = 0.5*(u(i1  ,i2  ,i3,vc)+u(i1,i2-1,i3,vc))
                            ! up = au
                            ! um = au
                            ! !picking side for the weights
                            !  if (up.ge.0) then
                            !    wpl = max(wplw,wprw)
                            !    wpr = min(wplw,wprw)
                            ! else
                            !    wpl = min(wplw,wprw)
                            !    wpr = max(wplw,wprw)
                            ! end if
                            ! if (um.ge.0) then
                            !    wml = max(wmlw,wmrw)
                            !    wmr = min(wmlw,wmrw)
                            ! else
                            !    wml = min(wmlw,wmrw)
                            !    wmr = max(wmlw,wmrw)
                            ! end if
                            ! expAd = 0.
                            ! if (orderOneFix.eq.1) then 
                            !    if (augmentPlot .eq. 1) then 
                            !       agu(vc,13)=0
                            !    end if
                            !    if (aum*aup .le. 0) then
                            !       aumax = max(abs(au),abs(aup),abs(aum))
                            !       expAd = (1./12.)*aumax*( u(i1,i2-2,i3,uc) - 4.*u(i1,i2-1,i3,uc) + 6.*u(i1,i2,i3,uc) - 4.*u(i1,i2+1,i3,uc) +  u(i1,i2+2,i3,uc))
                            !       wpl = 0.5
                            !       wpr = 0.5
                            !       wml = 0.5
                            !       wmr = 0.5
                            !       if (augmentPlot .eq. 1) then 
                            !          agu(vc,13)=au
                            !       end if
                            !    end if
                            ! end if
                            ! Fpl = 1./6.*(  -u(i1,i2-1,i3,uc) + 5.*u(i1,i2,i3,uc)   + 2.*u(i1,i2+1,i3,uc))
                            ! Fpr = 1./6.*( 2.*u(i1,i2,i3,uc)   + 5.*u(i1,i2+1,i3,uc) -   u(i1,i2+2,i3,uc))
                            ! Fml = 1./6.*(  -u(i1,i2-2,i3,uc) + 5.*u(i1,i2-1,i3,uc)   + 2.*u(i1,i2,i3,uc))
                            ! Fmr = 1./6.*( 2.*u(i1,i2-1,i3,uc) + 5.*u(i1,i2,i3,uc) -       u(i1,i2+1,i3,uc))
                            ! Fp  = wpl*Fpl + wpr*Fpr
                            ! Fm  = wml*Fml + wmr*Fmr
                            ! agu(vc,uc)= au*(Fp - Fm)/(dx(1)) + expAd*dx(1)**(-1)
                            ! if (augmentPlot .eq. 1) then 
                            !    agu(uc,8)=au
                            !    agu(uc,9) =wpl
                            !    agu(uc,10)=wpr
                            !    agu(uc,11)=wml
                            !    agu(uc,12)=wmr
                            ! end if
                            !v*vy---------------------------------------------------------------
                            ! Apl   = u(i1,i2+1,i3,vc) - 2.*u(i1,i2,i3,vc)   + u(i1,i2-1,i3,vc)
                            ! Bpl   = u(i1,i2+1,i3,vc) -    u(i1,i2-1,i3,vc)
                            ! Apr   = u(i1,i2+2,i3,vc) - 2.*u(i1,i2+1,i3,vc) + u(i1,i2,i3,vc)
                            ! Bpr   = u(i1,i2+2,i3,vc) -    u(i1,i2,i3,vc)
                            ! Aml   = u(i1,i2,i3,vc) - 2.*u(i1,i2-1,i3,vc) + u(i1,i2-2,i3,vc)
                            ! Bml   = u(i1,i2,i3,vc) -    u(i1,i2-2,i3,vc)
                            ! Amr   = Apl
                            ! Bmr   = Bpl
                            ! betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 1./4.*Bpl**2
                            ! betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 1./4.*Bpr**2
                            ! betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 1./4.*Bml**2
                            ! betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 1./4.*Bmr**2
                            ! ep     = 1e-15
                            ! ap1    = (1./2.)/(ep + betapl)**(4)
                            ! ap2    = (1./2.)/(ep + betapr)**(4)
                            ! am1    = (1./2.)/(ep + betaml)**(4)
                            ! am2    = (1./2.)/(ep + betamr)**(4)
                            ! !calculate weights
                            ! wp1 = ap1/(ap1+ap2)
                            ! wp2 = ap2/(ap1+ap2)
                            ! wm1 = am1/(am1+am2)
                            ! wm2 = am2/(am1+am2)
                            ! !mapping of the weights
                            ! bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                            ! bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                            ! bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                            ! bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                            ! wplw = bp1/(bp1+bp2)
                            ! wprw = bp2/(bp1+bp2)
                            ! wmlw = bm1/(bm1+bm2)
                            ! wmrw = bm2/(bm1+bm2)
                            ! !get face value of v
                            !  ! up = 0.5*(u(i1  ,i2,i3,vc)  + u(i1,i2+1,i3,vc))
                            !  ! um = 0.5*(u(i1-1,i2,i3,vc)  + u(i1,i2  ,i3,vc))
                            ! up = au
                            ! um = au
                            ! !picking side for the weights
                            !  if (up.ge.0) then
                            !    wpl = max(wplw,wprw)
                            !    wpr = min(wplw,wprw)
                            ! else
                            !    wpl = min(wplw,wprw)
                            !    wpr = max(wplw,wprw)
                            ! end if
                            ! if (um.ge.0) then
                            !    wml = max(wmlw,wmrw)
                            !    wmr = min(wmlw,wmrw)
                            ! else
                            !    wml = min(wmlw,wmrw)
                            !    wmr = max(wmlw,wmrw)
                            ! end if
                            ! expAd=0.
                            ! if (orderOneFix.eq.1) then
                            !    if (aum*aup .le. 0) then
                            !       aumax = max(abs(au),abs(aup),abs(aum))
                            !       expAd = (1./12.)*aumax*( u(i1,i2-2,i3,vc) - 4.*u(i1,i2-1,i3,vc) + 6.*u(i1,i2,i3,vc) - 4.*u(i1,i2+1,i3,vc) +  u(i1,i2+2,i3,vc))
                            !       wpl = 0.5
                            !       wpr = 0.5
                            !       wml = 0.5
                            !       wmr = 0.5
                            !    end if
                            ! end if
                            ! Fpl = 1./6.*(  -u(i1,i2-1,i3,vc) + 5.*u(i1,i2,i3,vc)   + 2.*u(i1,i2+1,i3,vc))
                            ! Fpr = 1./6.*( 2.*u(i1,i2,i3,vc)   + 5.*u(i1,i2+1,i3,vc) -   u(i1,i2+2,i3,vc))
                            ! Fml = 1./6.*(  -u(i1,i2-2,i3,vc) + 5.*u(i1,i2-1,i3,vc)   + 2.*u(i1,i2,i3,vc))
                            ! Fmr = 1./6.*( 2.*u(i1,i2-1,i3,vc) + 5.*u(i1,i2,i3,vc) -       u(i1,i2+1,i3,vc))
                            ! Fp  = wpl*Fpl + wpr*Fpr
                            ! Fm  = wml*Fml + wmr*Fmr
                            ! agu(vc,vc)= au*(Fp - Fm)/(dx(1)) + expAd*dx(1)**(-1)
                            ! if (augmentPlot .eq. 1) then  
                            !    agu(vc,9) =wpl
                            !    agu(vc,10)=wpr
                            !    agu(vc,11)=wml
                            !    agu(vc,12)=wmr
                            ! end if
                         else ! if (upwindOrder.eq.4)
                            write(*,'(" getBewnoAdvection: only 4th 
     & order is avaliable now")' )
                            stop 777
                         end if
                     ! #If ("BWENO" == "CENTERED") 
                     !  ! -- centered advection ---
                     !  ! write(*,'(" getAdvection -- centered")')
                     !  #If "2" == "2"
                     !   ! -- 2D --
                     !   agu(uc,uc)=UU(uc)*UX(uc)
                     !   agu(vc,uc)=UU(vc)*UY(uc)
                     !   agu(uc,vc)=UU(uc)*UX(vc)
                     !   agu(vc,vc)=UU(vc)*UY(vc)
                     !  #Else
                     !   ! -- 3D -- *check me*
                     !   agu(uc,uc)=UU(uc)*UX(uc)
                     !   agu(vc,uc)=UU(vc)*UY(uc)
                     !   agu(wc,uc)=UU(wc)*UY(uc)
                     !   agu(uc,vc)=UU(uc)*UX(vc)
                     !   agu(vc,vc)=UU(vc)*UY(vc)
                     !   agu(wc,vc)=UU(wc)*UY(vc)
                     !   agu(uc,wc)=UU(uc)*UX(wc)
                     !   agu(vc,wc)=UU(vc)*UY(wc)
                     !   agu(wc,wc)=UU(wc)*UY(wc)
                     !  #End
                     ! #Elif "BWENO" == "BWENO" 
                     !   ! --- upwind scheme ---
                     !   ! for testing output this next message:
                     !   if( t.le. 0. )then
                     !     write(*,'(" getAdvection upwind scheme (7)")') 
                     !   end if
                     !   getUpwindAdvection(u,i1,i2,i3,NONE,2,4,rectangular, agu)
                     ! #Elif "BWENO" == "BWENO" 
                     !   ! --- Bweno scheme ---
                     !   ! for testing output this message:
                     !   if( t.le. 0. )then
                     !      write(*,'(" getAdvection BWENO scheme (7)")') 
                     !   end if
                     !   getBwenoAdvection(u,i1,i2,i3,NONE,2,4,rectangular, agu)
                     ! #Else
                     !   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
                     !   stop 999
                     ! #End
                       ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-ux42r(
     & i1,i2,i3,pc)
                       ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-uy42r(
     & i1,i2,i3,pc)
                end do
                end do
                end do
               end if
             else
              write(*,*)'insdt: Unknown implicitOption=',implicitOption
              stop 7
             end if  ! end implicitOption
            else
              write(*,'(" unknown advectionOption")')
              stop 1010
            end if
           end if
          else if( isAxisymmetric.eq.1 )then
           if( advectionOption.ne.centeredAdvection )then
             write(*,*) 'insdt.h : finish me for axisymmetric'
             stop 2020
           end if
           ! **** axisymmetric case ****
           if( gridIsImplicit.eq.0 )then
            ! explicit
            if( useWhereMask.ne.0 )then
             do i3=n3a,n3b
             do i2=n2a,n2b
             do i1=n1a,n1b
              if( mask(i1,i2,i3).gt.0 )then
                 ! INS, no AD
                  ut(i1,i2,i3,uc)= -uu(i1,i2,i3,uc)*ux42r(i1,i2,i3,uc)-
     & uu(i1,i2,i3,vc)*uy42r(i1,i2,i3,uc)-ux42r(i1,i2,i3,pc)+nu*
     & ulaplacian42r(i1,i2,i3,uc)
                  ut(i1,i2,i3,vc)= -uu(i1,i2,i3,uc)*ux42r(i1,i2,i3,vc)-
     & uu(i1,i2,i3,vc)*uy42r(i1,i2,i3,vc)-uy42r(i1,i2,i3,pc)+nu*
     & ulaplacian42r(i1,i2,i3,vc)
                  ! -- add on axisymmetric corrections ---
                  ri=radiusInverse(i1,i2,i3)
                  if( ri.ne.0. )then
                    ut(i1,i2,i3,uc)=ut(i1,i2,i3,uc)+nu*( uy42r(i1,i2,
     & i3,uc)*ri )
                    ut(i1,i2,i3,vc)=ut(i1,i2,i3,vc)+nu*( (uy42r(i1,i2,
     & i3,vc)-uu(i1,i2,i3,vc)*ri)*ri )
                  else
                    ut(i1,i2,i3,uc)=ut(i1,i2,i3,uc)+nu*( uyy42r(i1,i2,
     & i3,uc) )
                    ut(i1,i2,i3,vc)=ut(i1,i2,i3,vc)+nu*( .5*uyy42r(i1,
     & i2,i3,vc) )
                  end if
              end if
             end do
             end do
             end do
            else
             do i3=n3a,n3b
             do i2=n2a,n2b
             do i1=n1a,n1b
                ! INS, no AD
                 ut(i1,i2,i3,uc)= -uu(i1,i2,i3,uc)*ux42r(i1,i2,i3,uc)-
     & uu(i1,i2,i3,vc)*uy42r(i1,i2,i3,uc)-ux42r(i1,i2,i3,pc)+nu*
     & ulaplacian42r(i1,i2,i3,uc)
                 ut(i1,i2,i3,vc)= -uu(i1,i2,i3,uc)*ux42r(i1,i2,i3,vc)-
     & uu(i1,i2,i3,vc)*uy42r(i1,i2,i3,vc)-uy42r(i1,i2,i3,pc)+nu*
     & ulaplacian42r(i1,i2,i3,vc)
                 ! -- add on axisymmetric corrections ---
                 ri=radiusInverse(i1,i2,i3)
                 if( ri.ne.0. )then
                   ut(i1,i2,i3,uc)=ut(i1,i2,i3,uc)+nu*( uy42r(i1,i2,i3,
     & uc)*ri )
                   ut(i1,i2,i3,vc)=ut(i1,i2,i3,vc)+nu*( (uy42r(i1,i2,
     & i3,vc)-uu(i1,i2,i3,vc)*ri)*ri )
                 else
                   ut(i1,i2,i3,uc)=ut(i1,i2,i3,uc)+nu*( uyy42r(i1,i2,
     & i3,uc) )
                   ut(i1,i2,i3,vc)=ut(i1,i2,i3,vc)+nu*( .5*uyy42r(i1,
     & i2,i3,vc) )
                 end if
             end do
             end do
             end do
            end if
           else ! gridIsImplicit
            ! ***** implicit *******
            if( implicitOption .eq.computeImplicitTermsSeparately )then
              if( useWhereMask.ne.0 )then
               do i3=n3a,n3b
               do i2=n2a,n2b
               do i1=n1a,n1b
                if( mask(i1,i2,i3).gt.0 )then
                   ! explicit terms only, no diffusion
                    ut(i1,i2,i3,uc)= -uu(i1,i2,i3,uc)*ux42r(i1,i2,i3,
     & uc)-uu(i1,i2,i3,vc)*uy42r(i1,i2,i3,uc)-ux42r(i1,i2,i3,pc)
                    ut(i1,i2,i3,vc)= -uu(i1,i2,i3,uc)*ux42r(i1,i2,i3,
     & vc)-uu(i1,i2,i3,vc)*uy42r(i1,i2,i3,vc)-uy42r(i1,i2,i3,pc)
                   ! include implicit terms - diffusion
                    uti(i1,i2,i3,uc)= nu*ulaplacian42r(i1,i2,i3,uc)
                    uti(i1,i2,i3,vc)= nu*ulaplacian42r(i1,i2,i3,vc)
                    ri=radiusInverse(i1,i2,i3)
                    if( ri.ne.0. )then
                      uti(i1,i2,i3,uc)=uti(i1,i2,i3,uc)+nu*( uy42r(i1,
     & i2,i3,uc)*ri )
                      uti(i1,i2,i3,vc)=uti(i1,i2,i3,vc)+nu*( (uy42r(i1,
     & i2,i3,vc)-uu(i1,i2,i3,vc)*ri)*ri )
                    else
                      uti(i1,i2,i3,uc)=uti(i1,i2,i3,uc)+nu*( uyy42r(i1,
     & i2,i3,uc) )
                      uti(i1,i2,i3,vc)=uti(i1,i2,i3,vc)+nu*( .5*uyy42r(
     & i1,i2,i3,vc) )
                    end if
                end if
               end do
               end do
               end do
              else
               do i3=n3a,n3b
               do i2=n2a,n2b
               do i1=n1a,n1b
                  ! explicit terms only, no diffusion
                   ut(i1,i2,i3,uc)= -uu(i1,i2,i3,uc)*ux42r(i1,i2,i3,uc)
     & -uu(i1,i2,i3,vc)*uy42r(i1,i2,i3,uc)-ux42r(i1,i2,i3,pc)
                   ut(i1,i2,i3,vc)= -uu(i1,i2,i3,uc)*ux42r(i1,i2,i3,vc)
     & -uu(i1,i2,i3,vc)*uy42r(i1,i2,i3,vc)-uy42r(i1,i2,i3,pc)
                  ! include implicit terms - diffusion
                   uti(i1,i2,i3,uc)= nu*ulaplacian42r(i1,i2,i3,uc)
                   uti(i1,i2,i3,vc)= nu*ulaplacian42r(i1,i2,i3,vc)
                   ri=radiusInverse(i1,i2,i3)
                   if( ri.ne.0. )then
                     uti(i1,i2,i3,uc)=uti(i1,i2,i3,uc)+nu*( uy42r(i1,
     & i2,i3,uc)*ri )
                     uti(i1,i2,i3,vc)=uti(i1,i2,i3,vc)+nu*( (uy42r(i1,
     & i2,i3,vc)-uu(i1,i2,i3,vc)*ri)*ri )
                   else
                     uti(i1,i2,i3,uc)=uti(i1,i2,i3,uc)+nu*( uyy42r(i1,
     & i2,i3,uc) )
                     uti(i1,i2,i3,vc)=uti(i1,i2,i3,vc)+nu*( .5*uyy42r(
     & i1,i2,i3,vc) )
                   end if
               end do
               end do
               end do
              end if
            else if( implicitOption.eq.doNotComputeImplicitTerms )then
              if( useWhereMask.ne.0 )then
               do i3=n3a,n3b
               do i2=n2a,n2b
               do i1=n1a,n1b
                if( mask(i1,i2,i3).gt.0 )then
                   ! explicit terms only, no diffusion
                    ut(i1,i2,i3,uc)= -uu(i1,i2,i3,uc)*ux42r(i1,i2,i3,
     & uc)-uu(i1,i2,i3,vc)*uy42r(i1,i2,i3,uc)-ux42r(i1,i2,i3,pc)
                    ut(i1,i2,i3,vc)= -uu(i1,i2,i3,uc)*ux42r(i1,i2,i3,
     & vc)-uu(i1,i2,i3,vc)*uy42r(i1,i2,i3,vc)-uy42r(i1,i2,i3,pc)
                end if
               end do
               end do
               end do
              else
               do i3=n3a,n3b
               do i2=n2a,n2b
               do i1=n1a,n1b
                  ! explicit terms only, no diffusion
                   ut(i1,i2,i3,uc)= -uu(i1,i2,i3,uc)*ux42r(i1,i2,i3,uc)
     & -uu(i1,i2,i3,vc)*uy42r(i1,i2,i3,uc)-ux42r(i1,i2,i3,pc)
                   ut(i1,i2,i3,vc)= -uu(i1,i2,i3,uc)*ux42r(i1,i2,i3,vc)
     & -uu(i1,i2,i3,vc)*uy42r(i1,i2,i3,vc)-uy42r(i1,i2,i3,pc)
               end do
               end do
               end do
              end if
            else
             write(*,*)'insdt: Unknown implicitOption=',implicitOption
             stop 5
            end if  ! end implicitOption
           end if
          else
            stop 88733
          end if
         else if( gridType.eq.curvilinear )then
          if( isAxisymmetric.eq.0 )then
           if( gridIsImplicit.eq.0 )then
            ! --- explicit time-stepping ---
            if( advectionOption.eq.centeredAdvection )then
              if( useWhereMask.ne.0 )then
               do i3=n3a,n3b
               do i2=n2a,n2b
               do i1=n1a,n1b
                if( mask(i1,i2,i3).gt.0 )then
                   ! INS, no AD
                    ut(i1,i2,i3,uc)= -uu(i1,i2,i3,uc)*ux42(i1,i2,i3,uc)
     & -uu(i1,i2,i3,vc)*uy42(i1,i2,i3,uc)-ux42(i1,i2,i3,pc)+nu*
     & ulaplacian42(i1,i2,i3,uc)
                    ut(i1,i2,i3,vc)= -uu(i1,i2,i3,uc)*ux42(i1,i2,i3,vc)
     & -uu(i1,i2,i3,vc)*uy42(i1,i2,i3,vc)-uy42(i1,i2,i3,pc)+nu*
     & ulaplacian42(i1,i2,i3,vc)
                end if
               end do
               end do
               end do
              else
               do i3=n3a,n3b
               do i2=n2a,n2b
               do i1=n1a,n1b
                  ! INS, no AD
                   ut(i1,i2,i3,uc)= -uu(i1,i2,i3,uc)*ux42(i1,i2,i3,uc)-
     & uu(i1,i2,i3,vc)*uy42(i1,i2,i3,uc)-ux42(i1,i2,i3,pc)+nu*
     & ulaplacian42(i1,i2,i3,uc)
                   ut(i1,i2,i3,vc)= -uu(i1,i2,i3,uc)*ux42(i1,i2,i3,vc)-
     & uu(i1,i2,i3,vc)*uy42(i1,i2,i3,vc)-uy42(i1,i2,i3,pc)+nu*
     & ulaplacian42(i1,i2,i3,vc)
               end do
               end do
               end do
              end if
            else if( advectionOption.eq.upwindAdvection )then
              ! --- upwind ---
              if( useWhereMask.ne.0 )then
               do i3=n3a,n3b
               do i2=n2a,n2b
               do i1=n1a,n1b
                if( mask(i1,i2,i3).gt.0 )then
                  ! INS, no AD
                    ! --- upwind approximations ---
                      ! --- upwind scheme ---
                      ! for testing output this next message:
                      ! if( t.le. 0. )then
                      !    write(*,'(" getAdvection upwind scheme (7)")') 
                      ! end if
                        ! write(*,'(" we in HERE")' )
                        ! stop 7171
                        ! --- CURVILINEAR GRID ---
                            if( upwindOrder.eq.1 )then
                               au   = rsxy(i1,i2,i3,0,0)*uu(i1,i2,i3,
     & uc)+rsxy(i1,i2,i3,0,1)*uu(i1,i2,i3,vc)
                               if( au.gt.0. )then
                                  ! u*ux
                                  agu(uc,uc)= au*(u(i1,i2,i3,uc)-u(i1-
     & 1,i2,i3,uc))/(dr(0))
                                  ! u*vx
                                  agu(uc,vc)= au*(u(i1,i2,i3,vc)-u(i1-
     & 1,i2,i3,vc))/(dr(0))
                               else
                                  ! u*ux
                                  agu(uc,uc)= au*(u(i1+1,i2,i3,uc)-u(
     & i1,i2,i3,uc))/(dr(0))
                                  ! u*vx
                                  agu(uc,vc)= au*(u(i1+1,i2,i3,vc)-u(
     & i1,i2,i3,vc))/(dr(0))
                               end if
                               au = rsxy(i1,i2,i3,1,0)*uu(i1,i2,i3,uc)+
     & rsxy(i1,i2,i3,1,1)*uu(i1,i2,i3,vc)
                               if( au.gt.0. )then
                                  ! v*uy
                                  agu(vc,uc)= au*(u(i1,i2,i3,uc)-u(i1,
     & i2-1,i3,uc))/(dr(1))
                                  ! v*vy
                                  agu(vc,vc)= au*(u(i1,i2,i3,vc)-u(i1,
     & i2-1,i3,vc))/(dr(1))
                               else
                                  ! v*uy
                                  agu(vc,uc)= au*(u(i1,i2+1,i3,uc)-u(
     & i1,i2,i3,uc))/(dr(1))
                                  ! v*vy
                                  agu(vc,vc)= au*(u(i1,i2+1,i3,vc)-u(
     & i1,i2,i3,vc))/(dr(1))
                               end if
                            elseif( upwindOrder.eq.2 )then
                               au = rsxy(i1,i2,i3,0,0)*uu(i1,i2,i3,uc)+
     & rsxy(i1,i2,i3,0,1)*uu(i1,i2,i3,vc)
                               if( au.gt.0. )then
                                ! u*ux = u*D-x(u)
                                  agu(uc,uc)= au*(3.*u(i1,i2,i3,uc)-4.*
     & u(i1-1,i2,i3,uc)+u(i1-2,i2,i3,uc))/(2.*dr(0))
                                ! u*vx = u*D-x(v)
                                  agu(uc,vc)= au*(3.*u(i1,i2,i3,vc)-4.*
     & u(i1-1,i2,i3,vc)+u(i1-2,i2,i3,vc))/(2.*dr(0))
                               else
                                ! u*ux = u*D+x(u)
                                  agu(uc,uc)= au*(-u(i1+2,i2,i3,uc)+4.*
     & u(i1+1,i2,i3,uc)-3.*u(i1,i2,i3,uc))/(2.*dr(0))
                                ! u*vx = u*D+x(v)
                                  agu(uc,vc)= au*(-u(i1+2,i2,i3,vc)+4.*
     & u(i1+1,i2,i3,vc)-3.*u(i1,i2,i3,vc))/(2.*dr(0))
                               end if
                               au = rsxy(i1,i2,i3,1,0)*uu(i1,i2,i3,uc)+
     & rsxy(i1,i2,i3,1,1)*uu(i1,i2,i3,vc)
                               if( au.gt.0. )then
                                ! v*uy = v*D-y(u)
                                  agu(vc,uc)= au*(3.*u(i1,i2,i3,uc)-4.*
     & u(i1,i2-1,i3,uc)+u(i1,i2-2,i3,uc))/(2.*dr(1))
                                ! v*vy = v*D-y(v)
                                  agu(vc,vc)= au*(3.*u(i1,i2,i3,vc)-4.*
     & u(i1,i2-1,i3,vc)+u(i1,i2-2,i3,vc))/(2.*dr(1))
                               else
                                ! v*uy = v*D+y(u)
                                  agu(vc,uc)= au*(-u(i1,i2+2,i3,uc)+4.*
     & u(i1,i2+1,i3,uc)-3.*u(i1,i2,i3,uc))/(2.*dr(1))
                                ! v*vy = v*D+y(v) 
                                  agu(vc,vc)= au*(-u(i1,i2+2,i3,vc)+4.*
     & u(i1,i2+1,i3,vc)-3.*u(i1,i2,i3,vc))/(2.*dr(1))
                               end if
                            elseif( upwindOrder.eq.3 )then
                               au = rsxy(i1,i2,i3,0,0)*uu(i1,i2,i3,uc)+
     & rsxy(i1,i2,i3,0,1)*uu(i1,i2,i3,vc)
                               if( au.gt.0. )then
                              ! u*ux = u*D-x(u)
                                  agu(uc,uc)= au*(2.*u(i1+1,i2,i3,uc)+
     & 3.*u(i1,i2,i3,uc)-6.*u(i1-1,i2,i3,uc)+u(i1-2,i2,i3,uc))/(6.*dr(
     & 0))
                              ! u*vx = u*D-x(v)
                                  agu(uc,vc)= au*(2.*u(i1+1,i2,i3,vc)+
     & 3.*u(i1,i2,i3,vc)-6.*u(i1-1,i2,i3,vc)+u(i1-2,i2,i3,vc))/(6.*dr(
     & 0))
                               else
                              ! u*ux = u*D+x(u)
                                  agu(uc,uc)= au*(-u(i1+2,i2,i3,uc)+6.*
     & u(i1+1,i2,i3,uc)-3.*u(i1,i2,i3,uc)-2.*u(i1-1,i2,i3,uc))/(6.*dr(
     & 0))
                              ! u*vx = u*D+x(v)
                                  agu(uc,vc)= au*(-u(i1+2,i2,i3,vc)+6.*
     & u(i1+1,i2,i3,vc)-3.*u(i1,i2,i3,vc)-2.*u(i1-1,i2,i3,vc))/(6.*dr(
     & 0))
                               end if
                               au = rsxy(i1,i2,i3,1,0)*uu(i1,i2,i3,uc)+
     & rsxy(i1,i2,i3,1,1)*uu(i1,i2,i3,vc)
                               if( au.gt.0. )then
                              ! v*uy = v*D-y(u)
                                  agu(vc,uc)= au*(2.*u(i1,i2+1,i3,uc)+
     & 3.*u(i1,i2,i3,uc)-6.*u(i1,i2-1,i3,uc)+u(i1,i2-2,i3,uc))/(6.*dr(
     & 1))
                              ! v*vy = v*D-y(v)
                                  agu(vc,vc)= au*(2.*u(i1,i2+1,i3,vc)+
     & 3.*u(i1,i2,i3,vc)-6.*u(i1,i2-1,i3,vc)+u(i1,i2-2,i3,vc))/(6.*dr(
     & 1))
                               else
                              ! v*uy = v*D+y(u)
                                  agu(vc,uc)= au*(-u(i1,i2+2,i3,uc)+6.*
     & u(i1,i2+1,i3,uc)-3.*u(i1,i2,i3,uc)-2.*u(i1,i2-1,i3,uc))/(6.*dr(
     & 1))
                              ! v*vy = v*D+y(v) 
                                  agu(vc,vc)= au*(-u(i1,i2+2,i3,vc)+6.*
     & u(i1,i2+1,i3,vc)-3.*u(i1,i2,i3,vc)-2.*u(i1,i2-1,i3,vc))/(6.*dr(
     & 1))
                               end if
                            end if
                      !
                      !      
                    ! #If ("UPWIND" == "CENTERED") 
                    !  ! -- centered advection ---
                    !  ! write(*,'(" getAdvection -- centered")')
                    !  #If "2" == "2"
                    !   ! -- 2D --
                    !   agu(uc,uc)=UU(uc)*UX(uc)
                    !   agu(vc,uc)=UU(vc)*UY(uc)
                    !   agu(uc,vc)=UU(uc)*UX(vc)
                    !   agu(vc,vc)=UU(vc)*UY(vc)
                    !  #Else
                    !   ! -- 3D -- *check me*
                    !   agu(uc,uc)=UU(uc)*UX(uc)
                    !   agu(vc,uc)=UU(vc)*UY(uc)
                    !   agu(wc,uc)=UU(wc)*UY(uc)
                    !   agu(uc,vc)=UU(uc)*UX(vc)
                    !   agu(vc,vc)=UU(vc)*UY(vc)
                    !   agu(wc,vc)=UU(wc)*UY(vc)
                    !   agu(uc,wc)=UU(uc)*UX(wc)
                    !   agu(vc,wc)=UU(vc)*UY(wc)
                    !   agu(wc,wc)=UU(wc)*UY(wc)
                    !  #End
                    ! #Elif "UPWIND" == "UPWIND" 
                    !   ! --- upwind scheme ---
                    !   ! for testing output this next message:
                    !   if( t.le. 0. )then
                    !     write(*,'(" getAdvection upwind scheme (7)")') 
                    !   end if
                    !   getUpwindAdvection(u,i1,i2,i3,NONE,2,4,curvilinear, agu)
                    ! #Elif "UPWIND" == "BWENO" 
                    !   ! --- Bweno scheme ---
                    !   ! for testing output this message:
                    !   if( t.le. 0. )then
                    !      write(*,'(" getAdvection BWENO scheme (7)")') 
                    !   end if
                    !   getBwenoAdvection(u,i1,i2,i3,NONE,2,4,curvilinear, agu)
                    ! #Else
                    !   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
                    !   stop 999
                    ! #End
                       ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-ux42(
     & i1,i2,i3,pc)+nu*ulaplacian42(i1,i2,i3,uc)
                       ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-uy42(
     & i1,i2,i3,pc)+nu*ulaplacian42(i1,i2,i3,vc)
                 !       u(i1,i2,i3,3)= agu(uc,3) ! au in r direction 
                 !       u(i1,i2,i3,4)= agu(uc,8) ! au in s direction
                 !       u(i1,i2,i3,5)= agu(uc,4)  ! wpl uur
                 !       u(i1,i2,i3,6)= agu(uc,5)  ! wpr uur
                 !       u(i1,i2,i3,7)= agu(uc,6)  ! wml uur
                 !       u(i1,i2,i3,8)= agu(uc,7)  ! wmr uur
                 !       u(i1,i2,i3,9)= agu(uc,9)  ! wpl uus
                 !       u(i1,i2,i3,10)= agu(uc,10)  ! wpr uus
                 !       u(i1,i2,i3,11)= agu(uc,11)  ! wml uus
                 !       u(i1,i2,i3,12)= agu(uc,12)  ! wmr uus
                 ! !
                 !       u(i1,i2,i3,13)= agu(vc,4)  ! wpl uvr
                 !       u(i1,i2,i3,14)= agu(vc,5)  ! wpr uvr
                 !       u(i1,i2,i3,15)= agu(vc,6)  ! wml uvr
                 !       u(i1,i2,i3,16)= agu(vc,7)  ! wmr uvr
                 !       u(i1,i2,i3,17)= agu(vc,9)  ! wpl uvs
                 !       u(i1,i2,i3,18)= agu(vc,10)  ! wpr uvs
                 !       u(i1,i2,i3,19)= agu(vc,11)  ! wml uvs
                 !       u(i1,i2,i3,20)= agu(vc,12)  ! wmr uvs
                 !       u(i1,i2,i3,21)= agu(uc,13)  ! au at where fix is turned on 
                 !       u(i1,i2,i3,22)= agu(vc,13)  ! au at where fix is turned on 
                end if
               end do
               end do
               end do
              else
               do i3=n3a,n3b
               do i2=n2a,n2b
               do i1=n1a,n1b
                 ! INS, no AD
                   ! --- upwind approximations ---
                     ! --- upwind scheme ---
                     ! for testing output this next message:
                     ! if( t.le. 0. )then
                     !    write(*,'(" getAdvection upwind scheme (7)")') 
                     ! end if
                       ! write(*,'(" we in HERE")' )
                       ! stop 7171
                       ! --- CURVILINEAR GRID ---
                           if( upwindOrder.eq.1 )then
                              au   = rsxy(i1,i2,i3,0,0)*uu(i1,i2,i3,uc)
     & +rsxy(i1,i2,i3,0,1)*uu(i1,i2,i3,vc)
                              if( au.gt.0. )then
                                 ! u*ux
                                 agu(uc,uc)= au*(u(i1,i2,i3,uc)-u(i1-1,
     & i2,i3,uc))/(dr(0))
                                 ! u*vx
                                 agu(uc,vc)= au*(u(i1,i2,i3,vc)-u(i1-1,
     & i2,i3,vc))/(dr(0))
                              else
                                 ! u*ux
                                 agu(uc,uc)= au*(u(i1+1,i2,i3,uc)-u(i1,
     & i2,i3,uc))/(dr(0))
                                 ! u*vx
                                 agu(uc,vc)= au*(u(i1+1,i2,i3,vc)-u(i1,
     & i2,i3,vc))/(dr(0))
                              end if
                              au = rsxy(i1,i2,i3,1,0)*uu(i1,i2,i3,uc)+
     & rsxy(i1,i2,i3,1,1)*uu(i1,i2,i3,vc)
                              if( au.gt.0. )then
                                 ! v*uy
                                 agu(vc,uc)= au*(u(i1,i2,i3,uc)-u(i1,
     & i2-1,i3,uc))/(dr(1))
                                 ! v*vy
                                 agu(vc,vc)= au*(u(i1,i2,i3,vc)-u(i1,
     & i2-1,i3,vc))/(dr(1))
                              else
                                 ! v*uy
                                 agu(vc,uc)= au*(u(i1,i2+1,i3,uc)-u(i1,
     & i2,i3,uc))/(dr(1))
                                 ! v*vy
                                 agu(vc,vc)= au*(u(i1,i2+1,i3,vc)-u(i1,
     & i2,i3,vc))/(dr(1))
                              end if
                           elseif( upwindOrder.eq.2 )then
                              au = rsxy(i1,i2,i3,0,0)*uu(i1,i2,i3,uc)+
     & rsxy(i1,i2,i3,0,1)*uu(i1,i2,i3,vc)
                              if( au.gt.0. )then
                               ! u*ux = u*D-x(u)
                                 agu(uc,uc)= au*(3.*u(i1,i2,i3,uc)-4.*
     & u(i1-1,i2,i3,uc)+u(i1-2,i2,i3,uc))/(2.*dr(0))
                               ! u*vx = u*D-x(v)
                                 agu(uc,vc)= au*(3.*u(i1,i2,i3,vc)-4.*
     & u(i1-1,i2,i3,vc)+u(i1-2,i2,i3,vc))/(2.*dr(0))
                              else
                               ! u*ux = u*D+x(u)
                                 agu(uc,uc)= au*(-u(i1+2,i2,i3,uc)+4.*
     & u(i1+1,i2,i3,uc)-3.*u(i1,i2,i3,uc))/(2.*dr(0))
                               ! u*vx = u*D+x(v)
                                 agu(uc,vc)= au*(-u(i1+2,i2,i3,vc)+4.*
     & u(i1+1,i2,i3,vc)-3.*u(i1,i2,i3,vc))/(2.*dr(0))
                              end if
                              au = rsxy(i1,i2,i3,1,0)*uu(i1,i2,i3,uc)+
     & rsxy(i1,i2,i3,1,1)*uu(i1,i2,i3,vc)
                              if( au.gt.0. )then
                               ! v*uy = v*D-y(u)
                                 agu(vc,uc)= au*(3.*u(i1,i2,i3,uc)-4.*
     & u(i1,i2-1,i3,uc)+u(i1,i2-2,i3,uc))/(2.*dr(1))
                               ! v*vy = v*D-y(v)
                                 agu(vc,vc)= au*(3.*u(i1,i2,i3,vc)-4.*
     & u(i1,i2-1,i3,vc)+u(i1,i2-2,i3,vc))/(2.*dr(1))
                              else
                               ! v*uy = v*D+y(u)
                                 agu(vc,uc)= au*(-u(i1,i2+2,i3,uc)+4.*
     & u(i1,i2+1,i3,uc)-3.*u(i1,i2,i3,uc))/(2.*dr(1))
                               ! v*vy = v*D+y(v) 
                                 agu(vc,vc)= au*(-u(i1,i2+2,i3,vc)+4.*
     & u(i1,i2+1,i3,vc)-3.*u(i1,i2,i3,vc))/(2.*dr(1))
                              end if
                           elseif( upwindOrder.eq.3 )then
                              au = rsxy(i1,i2,i3,0,0)*uu(i1,i2,i3,uc)+
     & rsxy(i1,i2,i3,0,1)*uu(i1,i2,i3,vc)
                              if( au.gt.0. )then
                             ! u*ux = u*D-x(u)
                                 agu(uc,uc)= au*(2.*u(i1+1,i2,i3,uc)+
     & 3.*u(i1,i2,i3,uc)-6.*u(i1-1,i2,i3,uc)+u(i1-2,i2,i3,uc))/(6.*dr(
     & 0))
                             ! u*vx = u*D-x(v)
                                 agu(uc,vc)= au*(2.*u(i1+1,i2,i3,vc)+
     & 3.*u(i1,i2,i3,vc)-6.*u(i1-1,i2,i3,vc)+u(i1-2,i2,i3,vc))/(6.*dr(
     & 0))
                              else
                             ! u*ux = u*D+x(u)
                                 agu(uc,uc)= au*(-u(i1+2,i2,i3,uc)+6.*
     & u(i1+1,i2,i3,uc)-3.*u(i1,i2,i3,uc)-2.*u(i1-1,i2,i3,uc))/(6.*dr(
     & 0))
                             ! u*vx = u*D+x(v)
                                 agu(uc,vc)= au*(-u(i1+2,i2,i3,vc)+6.*
     & u(i1+1,i2,i3,vc)-3.*u(i1,i2,i3,vc)-2.*u(i1-1,i2,i3,vc))/(6.*dr(
     & 0))
                              end if
                              au = rsxy(i1,i2,i3,1,0)*uu(i1,i2,i3,uc)+
     & rsxy(i1,i2,i3,1,1)*uu(i1,i2,i3,vc)
                              if( au.gt.0. )then
                             ! v*uy = v*D-y(u)
                                 agu(vc,uc)= au*(2.*u(i1,i2+1,i3,uc)+
     & 3.*u(i1,i2,i3,uc)-6.*u(i1,i2-1,i3,uc)+u(i1,i2-2,i3,uc))/(6.*dr(
     & 1))
                             ! v*vy = v*D-y(v)
                                 agu(vc,vc)= au*(2.*u(i1,i2+1,i3,vc)+
     & 3.*u(i1,i2,i3,vc)-6.*u(i1,i2-1,i3,vc)+u(i1,i2-2,i3,vc))/(6.*dr(
     & 1))
                              else
                             ! v*uy = v*D+y(u)
                                 agu(vc,uc)= au*(-u(i1,i2+2,i3,uc)+6.*
     & u(i1,i2+1,i3,uc)-3.*u(i1,i2,i3,uc)-2.*u(i1,i2-1,i3,uc))/(6.*dr(
     & 1))
                             ! v*vy = v*D+y(v) 
                                 agu(vc,vc)= au*(-u(i1,i2+2,i3,vc)+6.*
     & u(i1,i2+1,i3,vc)-3.*u(i1,i2,i3,vc)-2.*u(i1,i2-1,i3,vc))/(6.*dr(
     & 1))
                              end if
                           end if
                     !
                     !      
                   ! #If ("UPWIND" == "CENTERED") 
                   !  ! -- centered advection ---
                   !  ! write(*,'(" getAdvection -- centered")')
                   !  #If "2" == "2"
                   !   ! -- 2D --
                   !   agu(uc,uc)=UU(uc)*UX(uc)
                   !   agu(vc,uc)=UU(vc)*UY(uc)
                   !   agu(uc,vc)=UU(uc)*UX(vc)
                   !   agu(vc,vc)=UU(vc)*UY(vc)
                   !  #Else
                   !   ! -- 3D -- *check me*
                   !   agu(uc,uc)=UU(uc)*UX(uc)
                   !   agu(vc,uc)=UU(vc)*UY(uc)
                   !   agu(wc,uc)=UU(wc)*UY(uc)
                   !   agu(uc,vc)=UU(uc)*UX(vc)
                   !   agu(vc,vc)=UU(vc)*UY(vc)
                   !   agu(wc,vc)=UU(wc)*UY(vc)
                   !   agu(uc,wc)=UU(uc)*UX(wc)
                   !   agu(vc,wc)=UU(vc)*UY(wc)
                   !   agu(wc,wc)=UU(wc)*UY(wc)
                   !  #End
                   ! #Elif "UPWIND" == "UPWIND" 
                   !   ! --- upwind scheme ---
                   !   ! for testing output this next message:
                   !   if( t.le. 0. )then
                   !     write(*,'(" getAdvection upwind scheme (7)")') 
                   !   end if
                   !   getUpwindAdvection(u,i1,i2,i3,NONE,2,4,curvilinear, agu)
                   ! #Elif "UPWIND" == "BWENO" 
                   !   ! --- Bweno scheme ---
                   !   ! for testing output this message:
                   !   if( t.le. 0. )then
                   !      write(*,'(" getAdvection BWENO scheme (7)")') 
                   !   end if
                   !   getBwenoAdvection(u,i1,i2,i3,NONE,2,4,curvilinear, agu)
                   ! #Else
                   !   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
                   !   stop 999
                   ! #End
                      ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-ux42(
     & i1,i2,i3,pc)+nu*ulaplacian42(i1,i2,i3,uc)
                      ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-uy42(
     & i1,i2,i3,pc)+nu*ulaplacian42(i1,i2,i3,vc)
                !       u(i1,i2,i3,3)= agu(uc,3) ! au in r direction 
                !       u(i1,i2,i3,4)= agu(uc,8) ! au in s direction
                !       u(i1,i2,i3,5)= agu(uc,4)  ! wpl uur
                !       u(i1,i2,i3,6)= agu(uc,5)  ! wpr uur
                !       u(i1,i2,i3,7)= agu(uc,6)  ! wml uur
                !       u(i1,i2,i3,8)= agu(uc,7)  ! wmr uur
                !       u(i1,i2,i3,9)= agu(uc,9)  ! wpl uus
                !       u(i1,i2,i3,10)= agu(uc,10)  ! wpr uus
                !       u(i1,i2,i3,11)= agu(uc,11)  ! wml uus
                !       u(i1,i2,i3,12)= agu(uc,12)  ! wmr uus
                ! !
                !       u(i1,i2,i3,13)= agu(vc,4)  ! wpl uvr
                !       u(i1,i2,i3,14)= agu(vc,5)  ! wpr uvr
                !       u(i1,i2,i3,15)= agu(vc,6)  ! wml uvr
                !       u(i1,i2,i3,16)= agu(vc,7)  ! wmr uvr
                !       u(i1,i2,i3,17)= agu(vc,9)  ! wpl uvs
                !       u(i1,i2,i3,18)= agu(vc,10)  ! wpr uvs
                !       u(i1,i2,i3,19)= agu(vc,11)  ! wml uvs
                !       u(i1,i2,i3,20)= agu(vc,12)  ! wmr uvs
                !       u(i1,i2,i3,21)= agu(uc,13)  ! au at where fix is turned on 
                !       u(i1,i2,i3,22)= agu(vc,13)  ! au at where fix is turned on 
               end do
               end do
               end do
              end if
            else if( advectionOption.eq.bwenoAdvection )then
              ! --- bweno ---
              if( useWhereMask.ne.0 )then
               do i3=n3a,n3b
               do i2=n2a,n2b
               do i1=n1a,n1b
                if( mask(i1,i2,i3).gt.0 )then
                  ! INS, no AD
                    ! --- upwind approximations ---
                      ! --- Bweno scheme ---
                      ! for testing output this message:
                      ! if( t.le. 0. )then
                      !    write(*,'(" getAdvection BWENO scheme (7)")') 
                      ! end if
                         gvU = u(i1,i2,i3,uc) - uu(i1,i2,i3,uc)
                         gvV = u(i1,i2,i3,vc) - uu(i1,i2,i3,vc)
                         gvW = u(i1,i2,i3,wc) - uu(i1,i2,i3,wc)
                        !bweno curvulinear
                        orderOneFix = 1
                        augmentPlot = 0
                        if( upwindOrder.eq.4 )then
                           ! u*ux-------------------------------------------------------------------------------
                           ! au = rsxy(i1,i2,i3,0,0)*u(i1,i2,i3,uc)+rsxy(i1,i2,i3,0,1)*u(i1,i2,i3,vc) + rsxy(i1,i2,i3,0,2)*u(i1,i2,i3,wc)
                           ! aum   = rsxy(i1-1,i2,i3,0,0)*u(i1-1,i2,i3,uc)+rsxy(i1-1,i2,i3,0,1)*u(i1-1,i2,i3,vc)+ rsxy(i1-1,i2,i3,0,2)*u(i1-1,i2,i3,wc)
                           ! aup   = rsxy(i1+1,i2,i3,0,0)*u(i1+1,i2,i3,uc)+rsxy(i1+1,i2,i3,0,1)*u(i1+1,i2,i3,vc)+ rsxy(i1+1,i2,i3,0,2)*u(i1+1,i2,i3,wc)
                           au    = rsxy(i1  ,i2,i3,0,0)*uu(i1  ,i2,i3,
     & uc) + rsxy(i1  ,i2,i3,0,1) * uu(i1  ,i2,i3,vc)
                           aum   = rsxy(i1-1,i2,i3,0,0)*(u(i1-1,i2,i3,
     & uc) - gvU ) + rsxy(i1-1,i2,i3,0,1) * (u(i1-1,i2,i3,vc)- gvV )
                           aup   = rsxy(i1+1,i2,i3,0,0)*(u(i1+1,i2,i3,
     & uc) - gvU ) + rsxy(i1+1,i2,i3,0,1) * (u(i1+1,i2,i3,vc)- gvV )
                           s1=1
                           s2=0
                           s3=0
                           var = uc
                           vard= uc
                           drl = dr(0)
                             Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*u(
     & i1   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
                             Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    u(
     & i1-s1,i2-s2,i3-s3,var)
                             Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) - 
     & 2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
                             Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) -  
     &   u(i1   ,i2   ,i3   ,var)
                             Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,i2-
     & s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
                             Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,i2-
     & 2*s2,i3-2*s3,var)
                             Amr   = Apl
                             Bmr   = Bpl
                             betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 
     & 1./4.*Bpl**2
                             betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 
     & 1./4.*Bpr**2
                             betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 
     & 1./4.*Bml**2
                             betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 
     & 1./4.*Bmr**2
                             ep     = 1e-15
                             ap1    = (1./2.)/(ep + betapl)**(4)
                             ap2    = (1./2.)/(ep + betapr)**(4)
                             am1    = (1./2.)/(ep + betaml)**(4)
                             am2    = (1./2.)/(ep + betamr)**(4)
                             !calculate weights
                             wp1 = ap1/(ap1+ap2)
                             wp2 = ap2/(ap1+ap2)
                             wm1 = am1/(am1+am2)
                             wm2 = am2/(am1+am2)
                             !mapping of the weights
                             bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                             bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                             bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                             bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                             wplw = bp1/(bp1+bp2)
                             wprw = bp2/(bp1+bp2)
                             wmlw = bm1/(bm1+bm2)
                             wmrw = bm2/(bm1+bm2)
                             !get face value of u
                             up = au
                             um = au
                             !picking side for the weights
                             if (up.ge.0) then
                                wpl = max(wplw,wprw)
                                wpr = min(wplw,wprw)
                             else
                                wpl = min(wplw,wprw)
                                wpr = max(wplw,wprw)
                             end if
                             if (um.ge.0) then
                                wml = max(wmlw,wmrw)
                                wmr = min(wmlw,wmrw)
                             else
                                wml = min(wmlw,wmrw)
                                wmr = max(wmlw,wmrw)
                             end if
                             expAd = 0.
                             if (orderOneFix.eq.1) then
                                ! if (augmentPlot.eq.1) then
                                !    agu(uc,13)=0
                                ! end if
                                if (aum*aup .le. 0) then
                                   aumax = max(abs(au),abs(aup),abs(
     & aum))
                                   expAd = (1./12.)*aumax*( u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(i1,
     & i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*s2,
     & i3+2*s3,var))
                                   wpl = 0.5
                                   wpr = 0.5
                                   wml = 0.5
                                   wmr = 0.5
                                   ! plot au at where fix is turned on
                                   ! if (augmentPlot.eq.1) then
                                   !    agu(uc,13)=au
                                   ! end if
                                end if
                             end if
                             ! wpl = 0.5
                             ! wpr = 0.5
                             ! wml = 0.5
                             ! wmr = 0.5
                             Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,var) +
     &  5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
                             Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 5.*u(
     & i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
                             Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*s3,
     & var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
                             Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,var) 
     & + 5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
                             Fp  = wpl*Fpl + wpr*Fpr
                             Fm  = wml*Fml + wmr*Fmr
                                agu(vard,var)= au*(Fp - Fm)/drl + 
     & expAd/drl
                           !     if ( i1 .eq. 14 .and. i2.eq.16 ) then
                                   !     if ( i1 .eq. 14) then
                                   ! write(*,'(" betapl,betapr,betaml,betamr = ",4F10.8)')  betapl,betapr,betaml,betamr
                                   ! write(*,'(" ap1,ap2,am1,am2 = ",4e10.8)')  ap1,ap2,am1,am2
                                   ! write(*,'(" wp1,wpr,wm1,wm2 = ",4F10.8)')  wp1,wp2,wm1,wm2
                                   ! write(*,'(" bp1,bp2,bm1,bm2 = ",4F10.8)')  bp1,bp2,bm1,bm2
                                   ! write(*,'(" wplw,wprw,wmlw,wmrw = ",4F10.8)')  wplw,wprw,wmlw,wmrw
                                   ! write(*,'(" wpl,wpr,wml,wmr = ",4F10.8)')  wpl,wpr,wml,wmr
                                   ! write(*,'(" agu = ",F10.8)') agu(vard,var)
                           !     end if
                           !u*vx
                           s1=1
                           s2=0
                           s3=0
                           var = vc
                           vard= uc
                           drl = dr(0)
                             Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*u(
     & i1   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
                             Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    u(
     & i1-s1,i2-s2,i3-s3,var)
                             Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) - 
     & 2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
                             Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) -  
     &   u(i1   ,i2   ,i3   ,var)
                             Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,i2-
     & s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
                             Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,i2-
     & 2*s2,i3-2*s3,var)
                             Amr   = Apl
                             Bmr   = Bpl
                             betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 
     & 1./4.*Bpl**2
                             betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 
     & 1./4.*Bpr**2
                             betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 
     & 1./4.*Bml**2
                             betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 
     & 1./4.*Bmr**2
                             ep     = 1e-15
                             ap1    = (1./2.)/(ep + betapl)**(4)
                             ap2    = (1./2.)/(ep + betapr)**(4)
                             am1    = (1./2.)/(ep + betaml)**(4)
                             am2    = (1./2.)/(ep + betamr)**(4)
                             !calculate weights
                             wp1 = ap1/(ap1+ap2)
                             wp2 = ap2/(ap1+ap2)
                             wm1 = am1/(am1+am2)
                             wm2 = am2/(am1+am2)
                             !mapping of the weights
                             bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                             bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                             bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                             bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                             wplw = bp1/(bp1+bp2)
                             wprw = bp2/(bp1+bp2)
                             wmlw = bm1/(bm1+bm2)
                             wmrw = bm2/(bm1+bm2)
                             !get face value of u
                             up = au
                             um = au
                             !picking side for the weights
                             if (up.ge.0) then
                                wpl = max(wplw,wprw)
                                wpr = min(wplw,wprw)
                             else
                                wpl = min(wplw,wprw)
                                wpr = max(wplw,wprw)
                             end if
                             if (um.ge.0) then
                                wml = max(wmlw,wmrw)
                                wmr = min(wmlw,wmrw)
                             else
                                wml = min(wmlw,wmrw)
                                wmr = max(wmlw,wmrw)
                             end if
                             expAd = 0.
                             if (orderOneFix.eq.1) then
                                ! if (augmentPlot.eq.1) then
                                !    agu(uc,13)=0
                                ! end if
                                if (aum*aup .le. 0) then
                                   aumax = max(abs(au),abs(aup),abs(
     & aum))
                                   expAd = (1./12.)*aumax*( u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(i1,
     & i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*s2,
     & i3+2*s3,var))
                                   wpl = 0.5
                                   wpr = 0.5
                                   wml = 0.5
                                   wmr = 0.5
                                   ! plot au at where fix is turned on
                                   ! if (augmentPlot.eq.1) then
                                   !    agu(uc,13)=au
                                   ! end if
                                end if
                             end if
                             ! wpl = 0.5
                             ! wpr = 0.5
                             ! wml = 0.5
                             ! wmr = 0.5
                             Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,var) +
     &  5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
                             Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 5.*u(
     & i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
                             Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*s3,
     & var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
                             Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,var) 
     & + 5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
                             Fp  = wpl*Fpl + wpr*Fpr
                             Fm  = wml*Fml + wmr*Fmr
                                agu(vard,var)= au*(Fp - Fm)/drl + 
     & expAd/drl
                           !     if ( i1 .eq. 14 .and. i2.eq.16 ) then
                                   !     if ( i1 .eq. 14) then
                                   ! write(*,'(" betapl,betapr,betaml,betamr = ",4F10.8)')  betapl,betapr,betaml,betamr
                                   ! write(*,'(" ap1,ap2,am1,am2 = ",4e10.8)')  ap1,ap2,am1,am2
                                   ! write(*,'(" wp1,wpr,wm1,wm2 = ",4F10.8)')  wp1,wp2,wm1,wm2
                                   ! write(*,'(" bp1,bp2,bm1,bm2 = ",4F10.8)')  bp1,bp2,bm1,bm2
                                   ! write(*,'(" wplw,wprw,wmlw,wmrw = ",4F10.8)')  wplw,wprw,wmlw,wmrw
                                   ! write(*,'(" wpl,wpr,wml,wmr = ",4F10.8)')  wpl,wpr,wml,wmr
                                   ! write(*,'(" agu = ",F10.8)') agu(vard,var)
                           !     end if
                           !v*uy
                           au    = rsxy(i1,i2  ,i3,1,0)*uu(i1,i2  ,i3,
     & uc)+rsxy(i1,i2,i3,1,1)*uu(i1,i2,i3,vc)
                           aum   = rsxy(i1,i2-1,i3,1,0)*(u(i1,i2-1,i3,
     & uc) - gvU) + rsxy(i1,i2-1,i3,1,1)*(u(i1,i2-1,i3,vc) - gvV)
                           aup   = rsxy(i1,i2+1,i3,1,0)*(u(i1,i2+1,i3,
     & uc) - gvU) + rsxy(i1,i2+1,i3,1,1)*(u(i1,i2+1,i3,vc) - gvV)
                           s1=0
                           s2=1
                           s3=0
                           var = uc
                           vard= vc
                           drl = dr(1)
                             Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*u(
     & i1   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
                             Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    u(
     & i1-s1,i2-s2,i3-s3,var)
                             Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) - 
     & 2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
                             Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) -  
     &   u(i1   ,i2   ,i3   ,var)
                             Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,i2-
     & s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
                             Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,i2-
     & 2*s2,i3-2*s3,var)
                             Amr   = Apl
                             Bmr   = Bpl
                             betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 
     & 1./4.*Bpl**2
                             betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 
     & 1./4.*Bpr**2
                             betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 
     & 1./4.*Bml**2
                             betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 
     & 1./4.*Bmr**2
                             ep     = 1e-15
                             ap1    = (1./2.)/(ep + betapl)**(4)
                             ap2    = (1./2.)/(ep + betapr)**(4)
                             am1    = (1./2.)/(ep + betaml)**(4)
                             am2    = (1./2.)/(ep + betamr)**(4)
                             !calculate weights
                             wp1 = ap1/(ap1+ap2)
                             wp2 = ap2/(ap1+ap2)
                             wm1 = am1/(am1+am2)
                             wm2 = am2/(am1+am2)
                             !mapping of the weights
                             bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                             bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                             bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                             bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                             wplw = bp1/(bp1+bp2)
                             wprw = bp2/(bp1+bp2)
                             wmlw = bm1/(bm1+bm2)
                             wmrw = bm2/(bm1+bm2)
                             !get face value of u
                             up = au
                             um = au
                             !picking side for the weights
                             if (up.ge.0) then
                                wpl = max(wplw,wprw)
                                wpr = min(wplw,wprw)
                             else
                                wpl = min(wplw,wprw)
                                wpr = max(wplw,wprw)
                             end if
                             if (um.ge.0) then
                                wml = max(wmlw,wmrw)
                                wmr = min(wmlw,wmrw)
                             else
                                wml = min(wmlw,wmrw)
                                wmr = max(wmlw,wmrw)
                             end if
                             expAd = 0.
                             if (orderOneFix.eq.1) then
                                ! if (augmentPlot.eq.1) then
                                !    agu(uc,13)=0
                                ! end if
                                if (aum*aup .le. 0) then
                                   aumax = max(abs(au),abs(aup),abs(
     & aum))
                                   expAd = (1./12.)*aumax*( u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(i1,
     & i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*s2,
     & i3+2*s3,var))
                                   wpl = 0.5
                                   wpr = 0.5
                                   wml = 0.5
                                   wmr = 0.5
                                   ! plot au at where fix is turned on
                                   ! if (augmentPlot.eq.1) then
                                   !    agu(uc,13)=au
                                   ! end if
                                end if
                             end if
                             ! wpl = 0.5
                             ! wpr = 0.5
                             ! wml = 0.5
                             ! wmr = 0.5
                             Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,var) +
     &  5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
                             Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 5.*u(
     & i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
                             Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*s3,
     & var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
                             Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,var) 
     & + 5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
                             Fp  = wpl*Fpl + wpr*Fpr
                             Fm  = wml*Fml + wmr*Fmr
                                agu(vard,var)= au*(Fp - Fm)/drl + 
     & expAd/drl
                           !     if ( i1 .eq. 14 .and. i2.eq.16 ) then
                                   !     if ( i1 .eq. 14) then
                                   ! write(*,'(" betapl,betapr,betaml,betamr = ",4F10.8)')  betapl,betapr,betaml,betamr
                                   ! write(*,'(" ap1,ap2,am1,am2 = ",4e10.8)')  ap1,ap2,am1,am2
                                   ! write(*,'(" wp1,wpr,wm1,wm2 = ",4F10.8)')  wp1,wp2,wm1,wm2
                                   ! write(*,'(" bp1,bp2,bm1,bm2 = ",4F10.8)')  bp1,bp2,bm1,bm2
                                   ! write(*,'(" wplw,wprw,wmlw,wmrw = ",4F10.8)')  wplw,wprw,wmlw,wmrw
                                   ! write(*,'(" wpl,wpr,wml,wmr = ",4F10.8)')  wpl,wpr,wml,wmr
                                   ! write(*,'(" agu = ",F10.8)') agu(vard,var)
                           !     end if
                           !v*vy
                           s1=0
                           s2=1
                           s3=0
                           var = vc
                           vard= vc
                           drl = dr(1)
                             Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*u(
     & i1   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
                             Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    u(
     & i1-s1,i2-s2,i3-s3,var)
                             Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) - 
     & 2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
                             Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) -  
     &   u(i1   ,i2   ,i3   ,var)
                             Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,i2-
     & s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
                             Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,i2-
     & 2*s2,i3-2*s3,var)
                             Amr   = Apl
                             Bmr   = Bpl
                             betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 
     & 1./4.*Bpl**2
                             betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 
     & 1./4.*Bpr**2
                             betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 
     & 1./4.*Bml**2
                             betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 
     & 1./4.*Bmr**2
                             ep     = 1e-15
                             ap1    = (1./2.)/(ep + betapl)**(4)
                             ap2    = (1./2.)/(ep + betapr)**(4)
                             am1    = (1./2.)/(ep + betaml)**(4)
                             am2    = (1./2.)/(ep + betamr)**(4)
                             !calculate weights
                             wp1 = ap1/(ap1+ap2)
                             wp2 = ap2/(ap1+ap2)
                             wm1 = am1/(am1+am2)
                             wm2 = am2/(am1+am2)
                             !mapping of the weights
                             bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                             bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                             bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                             bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                             wplw = bp1/(bp1+bp2)
                             wprw = bp2/(bp1+bp2)
                             wmlw = bm1/(bm1+bm2)
                             wmrw = bm2/(bm1+bm2)
                             !get face value of u
                             up = au
                             um = au
                             !picking side for the weights
                             if (up.ge.0) then
                                wpl = max(wplw,wprw)
                                wpr = min(wplw,wprw)
                             else
                                wpl = min(wplw,wprw)
                                wpr = max(wplw,wprw)
                             end if
                             if (um.ge.0) then
                                wml = max(wmlw,wmrw)
                                wmr = min(wmlw,wmrw)
                             else
                                wml = min(wmlw,wmrw)
                                wmr = max(wmlw,wmrw)
                             end if
                             expAd = 0.
                             if (orderOneFix.eq.1) then
                                ! if (augmentPlot.eq.1) then
                                !    agu(uc,13)=0
                                ! end if
                                if (aum*aup .le. 0) then
                                   aumax = max(abs(au),abs(aup),abs(
     & aum))
                                   expAd = (1./12.)*aumax*( u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(i1,
     & i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*s2,
     & i3+2*s3,var))
                                   wpl = 0.5
                                   wpr = 0.5
                                   wml = 0.5
                                   wmr = 0.5
                                   ! plot au at where fix is turned on
                                   ! if (augmentPlot.eq.1) then
                                   !    agu(uc,13)=au
                                   ! end if
                                end if
                             end if
                             ! wpl = 0.5
                             ! wpr = 0.5
                             ! wml = 0.5
                             ! wmr = 0.5
                             Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,var) +
     &  5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
                             Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 5.*u(
     & i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
                             Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*s3,
     & var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
                             Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,var) 
     & + 5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
                             Fp  = wpl*Fpl + wpr*Fpr
                             Fm  = wml*Fml + wmr*Fmr
                                agu(vard,var)= au*(Fp - Fm)/drl + 
     & expAd/drl
                           !     if ( i1 .eq. 14 .and. i2.eq.16 ) then
                                   !     if ( i1 .eq. 14) then
                                   ! write(*,'(" betapl,betapr,betaml,betamr = ",4F10.8)')  betapl,betapr,betaml,betamr
                                   ! write(*,'(" ap1,ap2,am1,am2 = ",4e10.8)')  ap1,ap2,am1,am2
                                   ! write(*,'(" wp1,wpr,wm1,wm2 = ",4F10.8)')  wp1,wp2,wm1,wm2
                                   ! write(*,'(" bp1,bp2,bm1,bm2 = ",4F10.8)')  bp1,bp2,bm1,bm2
                                   ! write(*,'(" wplw,wprw,wmlw,wmrw = ",4F10.8)')  wplw,wprw,wmlw,wmrw
                                   ! write(*,'(" wpl,wpr,wml,wmr = ",4F10.8)')  wpl,wpr,wml,wmr
                                   ! write(*,'(" agu = ",F10.8)') agu(vard,var)
                           !     end if
                        else ! if (upwindOrder.eq.4)
                           write(*,'(" getBewnoAdvection: only 4th 
     & order is avaliable now")' )
                           stop 777
                        end if
                    ! #If ("BWENO" == "CENTERED") 
                    !  ! -- centered advection ---
                    !  ! write(*,'(" getAdvection -- centered")')
                    !  #If "2" == "2"
                    !   ! -- 2D --
                    !   agu(uc,uc)=UU(uc)*UX(uc)
                    !   agu(vc,uc)=UU(vc)*UY(uc)
                    !   agu(uc,vc)=UU(uc)*UX(vc)
                    !   agu(vc,vc)=UU(vc)*UY(vc)
                    !  #Else
                    !   ! -- 3D -- *check me*
                    !   agu(uc,uc)=UU(uc)*UX(uc)
                    !   agu(vc,uc)=UU(vc)*UY(uc)
                    !   agu(wc,uc)=UU(wc)*UY(uc)
                    !   agu(uc,vc)=UU(uc)*UX(vc)
                    !   agu(vc,vc)=UU(vc)*UY(vc)
                    !   agu(wc,vc)=UU(wc)*UY(vc)
                    !   agu(uc,wc)=UU(uc)*UX(wc)
                    !   agu(vc,wc)=UU(vc)*UY(wc)
                    !   agu(wc,wc)=UU(wc)*UY(wc)
                    !  #End
                    ! #Elif "BWENO" == "BWENO" 
                    !   ! --- upwind scheme ---
                    !   ! for testing output this next message:
                    !   if( t.le. 0. )then
                    !     write(*,'(" getAdvection upwind scheme (7)")') 
                    !   end if
                    !   getUpwindAdvection(u,i1,i2,i3,NONE,2,4,curvilinear, agu)
                    ! #Elif "BWENO" == "BWENO" 
                    !   ! --- Bweno scheme ---
                    !   ! for testing output this message:
                    !   if( t.le. 0. )then
                    !      write(*,'(" getAdvection BWENO scheme (7)")') 
                    !   end if
                    !   getBwenoAdvection(u,i1,i2,i3,NONE,2,4,curvilinear, agu)
                    ! #Else
                    !   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
                    !   stop 999
                    ! #End
                       ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-ux42(
     & i1,i2,i3,pc)+nu*ulaplacian42(i1,i2,i3,uc)
                       ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-uy42(
     & i1,i2,i3,pc)+nu*ulaplacian42(i1,i2,i3,vc)
                 !       u(i1,i2,i3,3)= agu(uc,3) ! au in r direction 
                 !       u(i1,i2,i3,4)= agu(uc,8) ! au in s direction
                 !       u(i1,i2,i3,5)= agu(uc,4)  ! wpl uur
                 !       u(i1,i2,i3,6)= agu(uc,5)  ! wpr uur
                 !       u(i1,i2,i3,7)= agu(uc,6)  ! wml uur
                 !       u(i1,i2,i3,8)= agu(uc,7)  ! wmr uur
                 !       u(i1,i2,i3,9)= agu(uc,9)  ! wpl uus
                 !       u(i1,i2,i3,10)= agu(uc,10)  ! wpr uus
                 !       u(i1,i2,i3,11)= agu(uc,11)  ! wml uus
                 !       u(i1,i2,i3,12)= agu(uc,12)  ! wmr uus
                 ! !
                 !       u(i1,i2,i3,13)= agu(vc,4)  ! wpl uvr
                 !       u(i1,i2,i3,14)= agu(vc,5)  ! wpr uvr
                 !       u(i1,i2,i3,15)= agu(vc,6)  ! wml uvr
                 !       u(i1,i2,i3,16)= agu(vc,7)  ! wmr uvr
                 !       u(i1,i2,i3,17)= agu(vc,9)  ! wpl uvs
                 !       u(i1,i2,i3,18)= agu(vc,10)  ! wpr uvs
                 !       u(i1,i2,i3,19)= agu(vc,11)  ! wml uvs
                 !       u(i1,i2,i3,20)= agu(vc,12)  ! wmr uvs
                 !       u(i1,i2,i3,21)= agu(uc,13)  ! au at where fix is turned on 
                 !       u(i1,i2,i3,22)= agu(vc,13)  ! au at where fix is turned on 
                end if
               end do
               end do
               end do
              else
               do i3=n3a,n3b
               do i2=n2a,n2b
               do i1=n1a,n1b
                 ! INS, no AD
                   ! --- upwind approximations ---
                     ! --- Bweno scheme ---
                     ! for testing output this message:
                     ! if( t.le. 0. )then
                     !    write(*,'(" getAdvection BWENO scheme (7)")') 
                     ! end if
                        gvU = u(i1,i2,i3,uc) - uu(i1,i2,i3,uc)
                        gvV = u(i1,i2,i3,vc) - uu(i1,i2,i3,vc)
                        gvW = u(i1,i2,i3,wc) - uu(i1,i2,i3,wc)
                       !bweno curvulinear
                       orderOneFix = 1
                       augmentPlot = 0
                       if( upwindOrder.eq.4 )then
                          ! u*ux-------------------------------------------------------------------------------
                          ! au = rsxy(i1,i2,i3,0,0)*u(i1,i2,i3,uc)+rsxy(i1,i2,i3,0,1)*u(i1,i2,i3,vc) + rsxy(i1,i2,i3,0,2)*u(i1,i2,i3,wc)
                          ! aum   = rsxy(i1-1,i2,i3,0,0)*u(i1-1,i2,i3,uc)+rsxy(i1-1,i2,i3,0,1)*u(i1-1,i2,i3,vc)+ rsxy(i1-1,i2,i3,0,2)*u(i1-1,i2,i3,wc)
                          ! aup   = rsxy(i1+1,i2,i3,0,0)*u(i1+1,i2,i3,uc)+rsxy(i1+1,i2,i3,0,1)*u(i1+1,i2,i3,vc)+ rsxy(i1+1,i2,i3,0,2)*u(i1+1,i2,i3,wc)
                          au    = rsxy(i1  ,i2,i3,0,0)*uu(i1  ,i2,i3,
     & uc) + rsxy(i1  ,i2,i3,0,1) * uu(i1  ,i2,i3,vc)
                          aum   = rsxy(i1-1,i2,i3,0,0)*(u(i1-1,i2,i3,
     & uc) - gvU ) + rsxy(i1-1,i2,i3,0,1) * (u(i1-1,i2,i3,vc)- gvV )
                          aup   = rsxy(i1+1,i2,i3,0,0)*(u(i1+1,i2,i3,
     & uc) - gvU ) + rsxy(i1+1,i2,i3,0,1) * (u(i1+1,i2,i3,vc)- gvV )
                          s1=1
                          s2=0
                          s3=0
                          var = uc
                          vard= uc
                          drl = dr(0)
                            Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*u(i1 
     &   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
                            Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    u(i1-
     & s1,i2-s2,i3-s3,var)
                            Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) - 
     & 2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
                            Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) -   
     &  u(i1   ,i2   ,i3   ,var)
                            Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,i2-
     & s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
                            Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,i2-
     & 2*s2,i3-2*s3,var)
                            Amr   = Apl
                            Bmr   = Bpl
                            betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 
     & 1./4.*Bpl**2
                            betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 
     & 1./4.*Bpr**2
                            betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 
     & 1./4.*Bml**2
                            betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 
     & 1./4.*Bmr**2
                            ep     = 1e-15
                            ap1    = (1./2.)/(ep + betapl)**(4)
                            ap2    = (1./2.)/(ep + betapr)**(4)
                            am1    = (1./2.)/(ep + betaml)**(4)
                            am2    = (1./2.)/(ep + betamr)**(4)
                            !calculate weights
                            wp1 = ap1/(ap1+ap2)
                            wp2 = ap2/(ap1+ap2)
                            wm1 = am1/(am1+am2)
                            wm2 = am2/(am1+am2)
                            !mapping of the weights
                            bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                            bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                            bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                            bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                            wplw = bp1/(bp1+bp2)
                            wprw = bp2/(bp1+bp2)
                            wmlw = bm1/(bm1+bm2)
                            wmrw = bm2/(bm1+bm2)
                            !get face value of u
                            up = au
                            um = au
                            !picking side for the weights
                            if (up.ge.0) then
                               wpl = max(wplw,wprw)
                               wpr = min(wplw,wprw)
                            else
                               wpl = min(wplw,wprw)
                               wpr = max(wplw,wprw)
                            end if
                            if (um.ge.0) then
                               wml = max(wmlw,wmrw)
                               wmr = min(wmlw,wmrw)
                            else
                               wml = min(wmlw,wmrw)
                               wmr = max(wmlw,wmrw)
                            end if
                            expAd = 0.
                            if (orderOneFix.eq.1) then
                               ! if (augmentPlot.eq.1) then
                               !    agu(uc,13)=0
                               ! end if
                               if (aum*aup .le. 0) then
                                  aumax = max(abs(au),abs(aup),abs(aum)
     & )
                                  expAd = (1./12.)*aumax*( u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(i1,
     & i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*s2,
     & i3+2*s3,var))
                                  wpl = 0.5
                                  wpr = 0.5
                                  wml = 0.5
                                  wmr = 0.5
                                  ! plot au at where fix is turned on
                                  ! if (augmentPlot.eq.1) then
                                  !    agu(uc,13)=au
                                  ! end if
                               end if
                            end if
                            ! wpl = 0.5
                            ! wpr = 0.5
                            ! wml = 0.5
                            ! wmr = 0.5
                            Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,var) + 
     & 5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
                            Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 5.*u(
     & i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
                            Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*s3,
     & var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
                            Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,var) +
     &  5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
                            Fp  = wpl*Fpl + wpr*Fpr
                            Fm  = wml*Fml + wmr*Fmr
                               agu(vard,var)= au*(Fp - Fm)/drl + 
     & expAd/drl
                          !     if ( i1 .eq. 14 .and. i2.eq.16 ) then
                                  !     if ( i1 .eq. 14) then
                                  ! write(*,'(" betapl,betapr,betaml,betamr = ",4F10.8)')  betapl,betapr,betaml,betamr
                                  ! write(*,'(" ap1,ap2,am1,am2 = ",4e10.8)')  ap1,ap2,am1,am2
                                  ! write(*,'(" wp1,wpr,wm1,wm2 = ",4F10.8)')  wp1,wp2,wm1,wm2
                                  ! write(*,'(" bp1,bp2,bm1,bm2 = ",4F10.8)')  bp1,bp2,bm1,bm2
                                  ! write(*,'(" wplw,wprw,wmlw,wmrw = ",4F10.8)')  wplw,wprw,wmlw,wmrw
                                  ! write(*,'(" wpl,wpr,wml,wmr = ",4F10.8)')  wpl,wpr,wml,wmr
                                  ! write(*,'(" agu = ",F10.8)') agu(vard,var)
                          !     end if
                          !u*vx
                          s1=1
                          s2=0
                          s3=0
                          var = vc
                          vard= uc
                          drl = dr(0)
                            Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*u(i1 
     &   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
                            Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    u(i1-
     & s1,i2-s2,i3-s3,var)
                            Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) - 
     & 2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
                            Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) -   
     &  u(i1   ,i2   ,i3   ,var)
                            Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,i2-
     & s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
                            Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,i2-
     & 2*s2,i3-2*s3,var)
                            Amr   = Apl
                            Bmr   = Bpl
                            betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 
     & 1./4.*Bpl**2
                            betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 
     & 1./4.*Bpr**2
                            betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 
     & 1./4.*Bml**2
                            betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 
     & 1./4.*Bmr**2
                            ep     = 1e-15
                            ap1    = (1./2.)/(ep + betapl)**(4)
                            ap2    = (1./2.)/(ep + betapr)**(4)
                            am1    = (1./2.)/(ep + betaml)**(4)
                            am2    = (1./2.)/(ep + betamr)**(4)
                            !calculate weights
                            wp1 = ap1/(ap1+ap2)
                            wp2 = ap2/(ap1+ap2)
                            wm1 = am1/(am1+am2)
                            wm2 = am2/(am1+am2)
                            !mapping of the weights
                            bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                            bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                            bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                            bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                            wplw = bp1/(bp1+bp2)
                            wprw = bp2/(bp1+bp2)
                            wmlw = bm1/(bm1+bm2)
                            wmrw = bm2/(bm1+bm2)
                            !get face value of u
                            up = au
                            um = au
                            !picking side for the weights
                            if (up.ge.0) then
                               wpl = max(wplw,wprw)
                               wpr = min(wplw,wprw)
                            else
                               wpl = min(wplw,wprw)
                               wpr = max(wplw,wprw)
                            end if
                            if (um.ge.0) then
                               wml = max(wmlw,wmrw)
                               wmr = min(wmlw,wmrw)
                            else
                               wml = min(wmlw,wmrw)
                               wmr = max(wmlw,wmrw)
                            end if
                            expAd = 0.
                            if (orderOneFix.eq.1) then
                               ! if (augmentPlot.eq.1) then
                               !    agu(uc,13)=0
                               ! end if
                               if (aum*aup .le. 0) then
                                  aumax = max(abs(au),abs(aup),abs(aum)
     & )
                                  expAd = (1./12.)*aumax*( u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(i1,
     & i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*s2,
     & i3+2*s3,var))
                                  wpl = 0.5
                                  wpr = 0.5
                                  wml = 0.5
                                  wmr = 0.5
                                  ! plot au at where fix is turned on
                                  ! if (augmentPlot.eq.1) then
                                  !    agu(uc,13)=au
                                  ! end if
                               end if
                            end if
                            ! wpl = 0.5
                            ! wpr = 0.5
                            ! wml = 0.5
                            ! wmr = 0.5
                            Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,var) + 
     & 5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
                            Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 5.*u(
     & i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
                            Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*s3,
     & var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
                            Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,var) +
     &  5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
                            Fp  = wpl*Fpl + wpr*Fpr
                            Fm  = wml*Fml + wmr*Fmr
                               agu(vard,var)= au*(Fp - Fm)/drl + 
     & expAd/drl
                          !     if ( i1 .eq. 14 .and. i2.eq.16 ) then
                                  !     if ( i1 .eq. 14) then
                                  ! write(*,'(" betapl,betapr,betaml,betamr = ",4F10.8)')  betapl,betapr,betaml,betamr
                                  ! write(*,'(" ap1,ap2,am1,am2 = ",4e10.8)')  ap1,ap2,am1,am2
                                  ! write(*,'(" wp1,wpr,wm1,wm2 = ",4F10.8)')  wp1,wp2,wm1,wm2
                                  ! write(*,'(" bp1,bp2,bm1,bm2 = ",4F10.8)')  bp1,bp2,bm1,bm2
                                  ! write(*,'(" wplw,wprw,wmlw,wmrw = ",4F10.8)')  wplw,wprw,wmlw,wmrw
                                  ! write(*,'(" wpl,wpr,wml,wmr = ",4F10.8)')  wpl,wpr,wml,wmr
                                  ! write(*,'(" agu = ",F10.8)') agu(vard,var)
                          !     end if
                          !v*uy
                          au    = rsxy(i1,i2  ,i3,1,0)*uu(i1,i2  ,i3,
     & uc)+rsxy(i1,i2,i3,1,1)*uu(i1,i2,i3,vc)
                          aum   = rsxy(i1,i2-1,i3,1,0)*(u(i1,i2-1,i3,
     & uc) - gvU) + rsxy(i1,i2-1,i3,1,1)*(u(i1,i2-1,i3,vc) - gvV)
                          aup   = rsxy(i1,i2+1,i3,1,0)*(u(i1,i2+1,i3,
     & uc) - gvU) + rsxy(i1,i2+1,i3,1,1)*(u(i1,i2+1,i3,vc) - gvV)
                          s1=0
                          s2=1
                          s3=0
                          var = uc
                          vard= vc
                          drl = dr(1)
                            Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*u(i1 
     &   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
                            Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    u(i1-
     & s1,i2-s2,i3-s3,var)
                            Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) - 
     & 2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
                            Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) -   
     &  u(i1   ,i2   ,i3   ,var)
                            Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,i2-
     & s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
                            Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,i2-
     & 2*s2,i3-2*s3,var)
                            Amr   = Apl
                            Bmr   = Bpl
                            betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 
     & 1./4.*Bpl**2
                            betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 
     & 1./4.*Bpr**2
                            betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 
     & 1./4.*Bml**2
                            betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 
     & 1./4.*Bmr**2
                            ep     = 1e-15
                            ap1    = (1./2.)/(ep + betapl)**(4)
                            ap2    = (1./2.)/(ep + betapr)**(4)
                            am1    = (1./2.)/(ep + betaml)**(4)
                            am2    = (1./2.)/(ep + betamr)**(4)
                            !calculate weights
                            wp1 = ap1/(ap1+ap2)
                            wp2 = ap2/(ap1+ap2)
                            wm1 = am1/(am1+am2)
                            wm2 = am2/(am1+am2)
                            !mapping of the weights
                            bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                            bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                            bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                            bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                            wplw = bp1/(bp1+bp2)
                            wprw = bp2/(bp1+bp2)
                            wmlw = bm1/(bm1+bm2)
                            wmrw = bm2/(bm1+bm2)
                            !get face value of u
                            up = au
                            um = au
                            !picking side for the weights
                            if (up.ge.0) then
                               wpl = max(wplw,wprw)
                               wpr = min(wplw,wprw)
                            else
                               wpl = min(wplw,wprw)
                               wpr = max(wplw,wprw)
                            end if
                            if (um.ge.0) then
                               wml = max(wmlw,wmrw)
                               wmr = min(wmlw,wmrw)
                            else
                               wml = min(wmlw,wmrw)
                               wmr = max(wmlw,wmrw)
                            end if
                            expAd = 0.
                            if (orderOneFix.eq.1) then
                               ! if (augmentPlot.eq.1) then
                               !    agu(uc,13)=0
                               ! end if
                               if (aum*aup .le. 0) then
                                  aumax = max(abs(au),abs(aup),abs(aum)
     & )
                                  expAd = (1./12.)*aumax*( u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(i1,
     & i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*s2,
     & i3+2*s3,var))
                                  wpl = 0.5
                                  wpr = 0.5
                                  wml = 0.5
                                  wmr = 0.5
                                  ! plot au at where fix is turned on
                                  ! if (augmentPlot.eq.1) then
                                  !    agu(uc,13)=au
                                  ! end if
                               end if
                            end if
                            ! wpl = 0.5
                            ! wpr = 0.5
                            ! wml = 0.5
                            ! wmr = 0.5
                            Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,var) + 
     & 5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
                            Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 5.*u(
     & i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
                            Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*s3,
     & var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
                            Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,var) +
     &  5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
                            Fp  = wpl*Fpl + wpr*Fpr
                            Fm  = wml*Fml + wmr*Fmr
                               agu(vard,var)= au*(Fp - Fm)/drl + 
     & expAd/drl
                          !     if ( i1 .eq. 14 .and. i2.eq.16 ) then
                                  !     if ( i1 .eq. 14) then
                                  ! write(*,'(" betapl,betapr,betaml,betamr = ",4F10.8)')  betapl,betapr,betaml,betamr
                                  ! write(*,'(" ap1,ap2,am1,am2 = ",4e10.8)')  ap1,ap2,am1,am2
                                  ! write(*,'(" wp1,wpr,wm1,wm2 = ",4F10.8)')  wp1,wp2,wm1,wm2
                                  ! write(*,'(" bp1,bp2,bm1,bm2 = ",4F10.8)')  bp1,bp2,bm1,bm2
                                  ! write(*,'(" wplw,wprw,wmlw,wmrw = ",4F10.8)')  wplw,wprw,wmlw,wmrw
                                  ! write(*,'(" wpl,wpr,wml,wmr = ",4F10.8)')  wpl,wpr,wml,wmr
                                  ! write(*,'(" agu = ",F10.8)') agu(vard,var)
                          !     end if
                          !v*vy
                          s1=0
                          s2=1
                          s3=0
                          var = vc
                          vard= vc
                          drl = dr(1)
                            Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*u(i1 
     &   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
                            Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    u(i1-
     & s1,i2-s2,i3-s3,var)
                            Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) - 
     & 2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
                            Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) -   
     &  u(i1   ,i2   ,i3   ,var)
                            Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,i2-
     & s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
                            Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,i2-
     & 2*s2,i3-2*s3,var)
                            Amr   = Apl
                            Bmr   = Bpl
                            betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 
     & 1./4.*Bpl**2
                            betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 
     & 1./4.*Bpr**2
                            betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 
     & 1./4.*Bml**2
                            betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 
     & 1./4.*Bmr**2
                            ep     = 1e-15
                            ap1    = (1./2.)/(ep + betapl)**(4)
                            ap2    = (1./2.)/(ep + betapr)**(4)
                            am1    = (1./2.)/(ep + betaml)**(4)
                            am2    = (1./2.)/(ep + betamr)**(4)
                            !calculate weights
                            wp1 = ap1/(ap1+ap2)
                            wp2 = ap2/(ap1+ap2)
                            wm1 = am1/(am1+am2)
                            wm2 = am2/(am1+am2)
                            !mapping of the weights
                            bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                            bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                            bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                            bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                            wplw = bp1/(bp1+bp2)
                            wprw = bp2/(bp1+bp2)
                            wmlw = bm1/(bm1+bm2)
                            wmrw = bm2/(bm1+bm2)
                            !get face value of u
                            up = au
                            um = au
                            !picking side for the weights
                            if (up.ge.0) then
                               wpl = max(wplw,wprw)
                               wpr = min(wplw,wprw)
                            else
                               wpl = min(wplw,wprw)
                               wpr = max(wplw,wprw)
                            end if
                            if (um.ge.0) then
                               wml = max(wmlw,wmrw)
                               wmr = min(wmlw,wmrw)
                            else
                               wml = min(wmlw,wmrw)
                               wmr = max(wmlw,wmrw)
                            end if
                            expAd = 0.
                            if (orderOneFix.eq.1) then
                               ! if (augmentPlot.eq.1) then
                               !    agu(uc,13)=0
                               ! end if
                               if (aum*aup .le. 0) then
                                  aumax = max(abs(au),abs(aup),abs(aum)
     & )
                                  expAd = (1./12.)*aumax*( u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(i1,
     & i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*s2,
     & i3+2*s3,var))
                                  wpl = 0.5
                                  wpr = 0.5
                                  wml = 0.5
                                  wmr = 0.5
                                  ! plot au at where fix is turned on
                                  ! if (augmentPlot.eq.1) then
                                  !    agu(uc,13)=au
                                  ! end if
                               end if
                            end if
                            ! wpl = 0.5
                            ! wpr = 0.5
                            ! wml = 0.5
                            ! wmr = 0.5
                            Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,var) + 
     & 5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
                            Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 5.*u(
     & i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
                            Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*s3,
     & var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
                            Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,var) +
     &  5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
                            Fp  = wpl*Fpl + wpr*Fpr
                            Fm  = wml*Fml + wmr*Fmr
                               agu(vard,var)= au*(Fp - Fm)/drl + 
     & expAd/drl
                          !     if ( i1 .eq. 14 .and. i2.eq.16 ) then
                                  !     if ( i1 .eq. 14) then
                                  ! write(*,'(" betapl,betapr,betaml,betamr = ",4F10.8)')  betapl,betapr,betaml,betamr
                                  ! write(*,'(" ap1,ap2,am1,am2 = ",4e10.8)')  ap1,ap2,am1,am2
                                  ! write(*,'(" wp1,wpr,wm1,wm2 = ",4F10.8)')  wp1,wp2,wm1,wm2
                                  ! write(*,'(" bp1,bp2,bm1,bm2 = ",4F10.8)')  bp1,bp2,bm1,bm2
                                  ! write(*,'(" wplw,wprw,wmlw,wmrw = ",4F10.8)')  wplw,wprw,wmlw,wmrw
                                  ! write(*,'(" wpl,wpr,wml,wmr = ",4F10.8)')  wpl,wpr,wml,wmr
                                  ! write(*,'(" agu = ",F10.8)') agu(vard,var)
                          !     end if
                       else ! if (upwindOrder.eq.4)
                          write(*,'(" getBewnoAdvection: only 4th 
     & order is avaliable now")' )
                          stop 777
                       end if
                   ! #If ("BWENO" == "CENTERED") 
                   !  ! -- centered advection ---
                   !  ! write(*,'(" getAdvection -- centered")')
                   !  #If "2" == "2"
                   !   ! -- 2D --
                   !   agu(uc,uc)=UU(uc)*UX(uc)
                   !   agu(vc,uc)=UU(vc)*UY(uc)
                   !   agu(uc,vc)=UU(uc)*UX(vc)
                   !   agu(vc,vc)=UU(vc)*UY(vc)
                   !  #Else
                   !   ! -- 3D -- *check me*
                   !   agu(uc,uc)=UU(uc)*UX(uc)
                   !   agu(vc,uc)=UU(vc)*UY(uc)
                   !   agu(wc,uc)=UU(wc)*UY(uc)
                   !   agu(uc,vc)=UU(uc)*UX(vc)
                   !   agu(vc,vc)=UU(vc)*UY(vc)
                   !   agu(wc,vc)=UU(wc)*UY(vc)
                   !   agu(uc,wc)=UU(uc)*UX(wc)
                   !   agu(vc,wc)=UU(vc)*UY(wc)
                   !   agu(wc,wc)=UU(wc)*UY(wc)
                   !  #End
                   ! #Elif "BWENO" == "BWENO" 
                   !   ! --- upwind scheme ---
                   !   ! for testing output this next message:
                   !   if( t.le. 0. )then
                   !     write(*,'(" getAdvection upwind scheme (7)")') 
                   !   end if
                   !   getUpwindAdvection(u,i1,i2,i3,NONE,2,4,curvilinear, agu)
                   ! #Elif "BWENO" == "BWENO" 
                   !   ! --- Bweno scheme ---
                   !   ! for testing output this message:
                   !   if( t.le. 0. )then
                   !      write(*,'(" getAdvection BWENO scheme (7)")') 
                   !   end if
                   !   getBwenoAdvection(u,i1,i2,i3,NONE,2,4,curvilinear, agu)
                   ! #Else
                   !   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
                   !   stop 999
                   ! #End
                      ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-ux42(
     & i1,i2,i3,pc)+nu*ulaplacian42(i1,i2,i3,uc)
                      ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-uy42(
     & i1,i2,i3,pc)+nu*ulaplacian42(i1,i2,i3,vc)
                !       u(i1,i2,i3,3)= agu(uc,3) ! au in r direction 
                !       u(i1,i2,i3,4)= agu(uc,8) ! au in s direction
                !       u(i1,i2,i3,5)= agu(uc,4)  ! wpl uur
                !       u(i1,i2,i3,6)= agu(uc,5)  ! wpr uur
                !       u(i1,i2,i3,7)= agu(uc,6)  ! wml uur
                !       u(i1,i2,i3,8)= agu(uc,7)  ! wmr uur
                !       u(i1,i2,i3,9)= agu(uc,9)  ! wpl uus
                !       u(i1,i2,i3,10)= agu(uc,10)  ! wpr uus
                !       u(i1,i2,i3,11)= agu(uc,11)  ! wml uus
                !       u(i1,i2,i3,12)= agu(uc,12)  ! wmr uus
                ! !
                !       u(i1,i2,i3,13)= agu(vc,4)  ! wpl uvr
                !       u(i1,i2,i3,14)= agu(vc,5)  ! wpr uvr
                !       u(i1,i2,i3,15)= agu(vc,6)  ! wml uvr
                !       u(i1,i2,i3,16)= agu(vc,7)  ! wmr uvr
                !       u(i1,i2,i3,17)= agu(vc,9)  ! wpl uvs
                !       u(i1,i2,i3,18)= agu(vc,10)  ! wpr uvs
                !       u(i1,i2,i3,19)= agu(vc,11)  ! wml uvs
                !       u(i1,i2,i3,20)= agu(vc,12)  ! wmr uvs
                !       u(i1,i2,i3,21)= agu(uc,13)  ! au at where fix is turned on 
                !       u(i1,i2,i3,22)= agu(vc,13)  ! au at where fix is turned on 
               end do
               end do
               end do
              end if
            else
              write(*,'(" unknown advectionOption")')
              stop 1010
            end if
           else ! gridIsImplicit
            ! ---- implicit time-stepping ---
            if( advectionOption.eq.centeredAdvection )then
             if( implicitOption .eq.computeImplicitTermsSeparately )
     & then
               if( useWhereMask.ne.0 )then
                do i3=n3a,n3b
                do i2=n2a,n2b
                do i1=n1a,n1b
                 if( mask(i1,i2,i3).gt.0 )then
                    ! explicit terms only, no diffusion
                     ut(i1,i2,i3,uc)= -uu(i1,i2,i3,uc)*ux42(i1,i2,i3,
     & uc)-uu(i1,i2,i3,vc)*uy42(i1,i2,i3,uc)-ux42(i1,i2,i3,pc)
                     ut(i1,i2,i3,vc)= -uu(i1,i2,i3,uc)*ux42(i1,i2,i3,
     & vc)-uu(i1,i2,i3,vc)*uy42(i1,i2,i3,vc)-uy42(i1,i2,i3,pc)
                    ! include implicit terms - diffusion
                     uti(i1,i2,i3,uc)= nu*ulaplacian42(i1,i2,i3,uc)
                     uti(i1,i2,i3,vc)= nu*ulaplacian42(i1,i2,i3,vc)
                 end if
                end do
                end do
                end do
               else
                do i3=n3a,n3b
                do i2=n2a,n2b
                do i1=n1a,n1b
                   ! explicit terms only, no diffusion
                    ut(i1,i2,i3,uc)= -uu(i1,i2,i3,uc)*ux42(i1,i2,i3,uc)
     & -uu(i1,i2,i3,vc)*uy42(i1,i2,i3,uc)-ux42(i1,i2,i3,pc)
                    ut(i1,i2,i3,vc)= -uu(i1,i2,i3,uc)*ux42(i1,i2,i3,vc)
     & -uu(i1,i2,i3,vc)*uy42(i1,i2,i3,vc)-uy42(i1,i2,i3,pc)
                   ! include implicit terms - diffusion
                    uti(i1,i2,i3,uc)= nu*ulaplacian42(i1,i2,i3,uc)
                    uti(i1,i2,i3,vc)= nu*ulaplacian42(i1,i2,i3,vc)
                end do
                end do
                end do
               end if
             else if( implicitOption.eq.doNotComputeImplicitTerms )then
               if( useWhereMask.ne.0 )then
                do i3=n3a,n3b
                do i2=n2a,n2b
                do i1=n1a,n1b
                 if( mask(i1,i2,i3).gt.0 )then
                    ! explicit terms only, no diffusion
                     ut(i1,i2,i3,uc)= -uu(i1,i2,i3,uc)*ux42(i1,i2,i3,
     & uc)-uu(i1,i2,i3,vc)*uy42(i1,i2,i3,uc)-ux42(i1,i2,i3,pc)
                     ut(i1,i2,i3,vc)= -uu(i1,i2,i3,uc)*ux42(i1,i2,i3,
     & vc)-uu(i1,i2,i3,vc)*uy42(i1,i2,i3,vc)-uy42(i1,i2,i3,pc)
                 end if
                end do
                end do
                end do
               else
                do i3=n3a,n3b
                do i2=n2a,n2b
                do i1=n1a,n1b
                   ! explicit terms only, no diffusion
                    ut(i1,i2,i3,uc)= -uu(i1,i2,i3,uc)*ux42(i1,i2,i3,uc)
     & -uu(i1,i2,i3,vc)*uy42(i1,i2,i3,uc)-ux42(i1,i2,i3,pc)
                    ut(i1,i2,i3,vc)= -uu(i1,i2,i3,uc)*ux42(i1,i2,i3,vc)
     & -uu(i1,i2,i3,vc)*uy42(i1,i2,i3,vc)-uy42(i1,i2,i3,pc)
                end do
                end do
                end do
               end if
             else
              write(*,*)'insdt: Unknown implicitOption=',implicitOption
              stop 5
             end if  ! end implicitOption
            else if( advectionOption.eq.upwindAdvection )then
              ! --- upwind ---
             if( implicitOption .eq.computeImplicitTermsSeparately )
     & then
               if( useWhereMask.ne.0 )then
                do i3=n3a,n3b
                do i2=n2a,n2b
                do i1=n1a,n1b
                 if( mask(i1,i2,i3).gt.0 )then
                    ! explicit terms only, no diffusion
                      ! upwind approximation/bweno
                        ! --- upwind scheme ---
                        ! for testing output this next message:
                        ! if( t.le. 0. )then
                        !    write(*,'(" getAdvection upwind scheme (7)")') 
                        ! end if
                          ! write(*,'(" we in HERE")' )
                          ! stop 7171
                          ! --- CURVILINEAR GRID ---
                              if( upwindOrder.eq.1 )then
                                 au   = rsxy(i1,i2,i3,0,0)*uu(i1,i2,i3,
     & uc)+rsxy(i1,i2,i3,0,1)*uu(i1,i2,i3,vc)
                                 if( au.gt.0. )then
                                    ! u*ux
                                    agu(uc,uc)= au*(u(i1,i2,i3,uc)-u(
     & i1-1,i2,i3,uc))/(dr(0))
                                    ! u*vx
                                    agu(uc,vc)= au*(u(i1,i2,i3,vc)-u(
     & i1-1,i2,i3,vc))/(dr(0))
                                 else
                                    ! u*ux
                                    agu(uc,uc)= au*(u(i1+1,i2,i3,uc)-u(
     & i1,i2,i3,uc))/(dr(0))
                                    ! u*vx
                                    agu(uc,vc)= au*(u(i1+1,i2,i3,vc)-u(
     & i1,i2,i3,vc))/(dr(0))
                                 end if
                                 au = rsxy(i1,i2,i3,1,0)*uu(i1,i2,i3,
     & uc)+rsxy(i1,i2,i3,1,1)*uu(i1,i2,i3,vc)
                                 if( au.gt.0. )then
                                    ! v*uy
                                    agu(vc,uc)= au*(u(i1,i2,i3,uc)-u(
     & i1,i2-1,i3,uc))/(dr(1))
                                    ! v*vy
                                    agu(vc,vc)= au*(u(i1,i2,i3,vc)-u(
     & i1,i2-1,i3,vc))/(dr(1))
                                 else
                                    ! v*uy
                                    agu(vc,uc)= au*(u(i1,i2+1,i3,uc)-u(
     & i1,i2,i3,uc))/(dr(1))
                                    ! v*vy
                                    agu(vc,vc)= au*(u(i1,i2+1,i3,vc)-u(
     & i1,i2,i3,vc))/(dr(1))
                                 end if
                              elseif( upwindOrder.eq.2 )then
                                 au = rsxy(i1,i2,i3,0,0)*uu(i1,i2,i3,
     & uc)+rsxy(i1,i2,i3,0,1)*uu(i1,i2,i3,vc)
                                 if( au.gt.0. )then
                                  ! u*ux = u*D-x(u)
                                    agu(uc,uc)= au*(3.*u(i1,i2,i3,uc)-
     & 4.*u(i1-1,i2,i3,uc)+u(i1-2,i2,i3,uc))/(2.*dr(0))
                                  ! u*vx = u*D-x(v)
                                    agu(uc,vc)= au*(3.*u(i1,i2,i3,vc)-
     & 4.*u(i1-1,i2,i3,vc)+u(i1-2,i2,i3,vc))/(2.*dr(0))
                                 else
                                  ! u*ux = u*D+x(u)
                                    agu(uc,uc)= au*(-u(i1+2,i2,i3,uc)+
     & 4.*u(i1+1,i2,i3,uc)-3.*u(i1,i2,i3,uc))/(2.*dr(0))
                                  ! u*vx = u*D+x(v)
                                    agu(uc,vc)= au*(-u(i1+2,i2,i3,vc)+
     & 4.*u(i1+1,i2,i3,vc)-3.*u(i1,i2,i3,vc))/(2.*dr(0))
                                 end if
                                 au = rsxy(i1,i2,i3,1,0)*uu(i1,i2,i3,
     & uc)+rsxy(i1,i2,i3,1,1)*uu(i1,i2,i3,vc)
                                 if( au.gt.0. )then
                                  ! v*uy = v*D-y(u)
                                    agu(vc,uc)= au*(3.*u(i1,i2,i3,uc)-
     & 4.*u(i1,i2-1,i3,uc)+u(i1,i2-2,i3,uc))/(2.*dr(1))
                                  ! v*vy = v*D-y(v)
                                    agu(vc,vc)= au*(3.*u(i1,i2,i3,vc)-
     & 4.*u(i1,i2-1,i3,vc)+u(i1,i2-2,i3,vc))/(2.*dr(1))
                                 else
                                  ! v*uy = v*D+y(u)
                                    agu(vc,uc)= au*(-u(i1,i2+2,i3,uc)+
     & 4.*u(i1,i2+1,i3,uc)-3.*u(i1,i2,i3,uc))/(2.*dr(1))
                                  ! v*vy = v*D+y(v) 
                                    agu(vc,vc)= au*(-u(i1,i2+2,i3,vc)+
     & 4.*u(i1,i2+1,i3,vc)-3.*u(i1,i2,i3,vc))/(2.*dr(1))
                                 end if
                              elseif( upwindOrder.eq.3 )then
                                 au = rsxy(i1,i2,i3,0,0)*uu(i1,i2,i3,
     & uc)+rsxy(i1,i2,i3,0,1)*uu(i1,i2,i3,vc)
                                 if( au.gt.0. )then
                                ! u*ux = u*D-x(u)
                                    agu(uc,uc)= au*(2.*u(i1+1,i2,i3,uc)
     & +3.*u(i1,i2,i3,uc)-6.*u(i1-1,i2,i3,uc)+u(i1-2,i2,i3,uc))/(6.*
     & dr(0))
                                ! u*vx = u*D-x(v)
                                    agu(uc,vc)= au*(2.*u(i1+1,i2,i3,vc)
     & +3.*u(i1,i2,i3,vc)-6.*u(i1-1,i2,i3,vc)+u(i1-2,i2,i3,vc))/(6.*
     & dr(0))
                                 else
                                ! u*ux = u*D+x(u)
                                    agu(uc,uc)= au*(-u(i1+2,i2,i3,uc)+
     & 6.*u(i1+1,i2,i3,uc)-3.*u(i1,i2,i3,uc)-2.*u(i1-1,i2,i3,uc))/(6.*
     & dr(0))
                                ! u*vx = u*D+x(v)
                                    agu(uc,vc)= au*(-u(i1+2,i2,i3,vc)+
     & 6.*u(i1+1,i2,i3,vc)-3.*u(i1,i2,i3,vc)-2.*u(i1-1,i2,i3,vc))/(6.*
     & dr(0))
                                 end if
                                 au = rsxy(i1,i2,i3,1,0)*uu(i1,i2,i3,
     & uc)+rsxy(i1,i2,i3,1,1)*uu(i1,i2,i3,vc)
                                 if( au.gt.0. )then
                                ! v*uy = v*D-y(u)
                                    agu(vc,uc)= au*(2.*u(i1,i2+1,i3,uc)
     & +3.*u(i1,i2,i3,uc)-6.*u(i1,i2-1,i3,uc)+u(i1,i2-2,i3,uc))/(6.*
     & dr(1))
                                ! v*vy = v*D-y(v)
                                    agu(vc,vc)= au*(2.*u(i1,i2+1,i3,vc)
     & +3.*u(i1,i2,i3,vc)-6.*u(i1,i2-1,i3,vc)+u(i1,i2-2,i3,vc))/(6.*
     & dr(1))
                                 else
                                ! v*uy = v*D+y(u)
                                    agu(vc,uc)= au*(-u(i1,i2+2,i3,uc)+
     & 6.*u(i1,i2+1,i3,uc)-3.*u(i1,i2,i3,uc)-2.*u(i1,i2-1,i3,uc))/(6.*
     & dr(1))
                                ! v*vy = v*D+y(v) 
                                    agu(vc,vc)= au*(-u(i1,i2+2,i3,vc)+
     & 6.*u(i1,i2+1,i3,vc)-3.*u(i1,i2,i3,vc)-2.*u(i1,i2-1,i3,vc))/(6.*
     & dr(1))
                                 end if
                              end if
                        !
                        !      
                      ! #If ("UPWIND" == "CENTERED") 
                      !  ! -- centered advection ---
                      !  ! write(*,'(" getAdvection -- centered")')
                      !  #If "2" == "2"
                      !   ! -- 2D --
                      !   agu(uc,uc)=UU(uc)*UX(uc)
                      !   agu(vc,uc)=UU(vc)*UY(uc)
                      !   agu(uc,vc)=UU(uc)*UX(vc)
                      !   agu(vc,vc)=UU(vc)*UY(vc)
                      !  #Else
                      !   ! -- 3D -- *check me*
                      !   agu(uc,uc)=UU(uc)*UX(uc)
                      !   agu(vc,uc)=UU(vc)*UY(uc)
                      !   agu(wc,uc)=UU(wc)*UY(uc)
                      !   agu(uc,vc)=UU(uc)*UX(vc)
                      !   agu(vc,vc)=UU(vc)*UY(vc)
                      !   agu(wc,vc)=UU(wc)*UY(vc)
                      !   agu(uc,wc)=UU(uc)*UX(wc)
                      !   agu(vc,wc)=UU(vc)*UY(wc)
                      !   agu(wc,wc)=UU(wc)*UY(wc)
                      !  #End
                      ! #Elif "UPWIND" == "UPWIND" 
                      !   ! --- upwind scheme ---
                      !   ! for testing output this next message:
                      !   if( t.le. 0. )then
                      !     write(*,'(" getAdvection upwind scheme (7)")') 
                      !   end if
                      !   getUpwindAdvection(u,i1,i2,i3,NONE,2,4,curvilinear, agu)
                      ! #Elif "UPWIND" == "BWENO" 
                      !   ! --- Bweno scheme ---
                      !   ! for testing output this message:
                      !   if( t.le. 0. )then
                      !      write(*,'(" getAdvection BWENO scheme (7)")') 
                      !   end if
                      !   getBwenoAdvection(u,i1,i2,i3,NONE,2,4,curvilinear, agu)
                      ! #Else
                      !   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
                      !   stop 999
                      ! #End
                        ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-ux42(
     & i1,i2,i3,pc)
                        ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-uy42(
     & i1,i2,i3,pc)
                    ! include implicit terms - diffusion
                     uti(i1,i2,i3,uc)= nu*ulaplacian42(i1,i2,i3,uc)
                     uti(i1,i2,i3,vc)= nu*ulaplacian42(i1,i2,i3,vc)
                 end if
                end do
                end do
                end do
               else
                do i3=n3a,n3b
                do i2=n2a,n2b
                do i1=n1a,n1b
                   ! explicit terms only, no diffusion
                     ! upwind approximation/bweno
                       ! --- upwind scheme ---
                       ! for testing output this next message:
                       ! if( t.le. 0. )then
                       !    write(*,'(" getAdvection upwind scheme (7)")') 
                       ! end if
                         ! write(*,'(" we in HERE")' )
                         ! stop 7171
                         ! --- CURVILINEAR GRID ---
                             if( upwindOrder.eq.1 )then
                                au   = rsxy(i1,i2,i3,0,0)*uu(i1,i2,i3,
     & uc)+rsxy(i1,i2,i3,0,1)*uu(i1,i2,i3,vc)
                                if( au.gt.0. )then
                                   ! u*ux
                                   agu(uc,uc)= au*(u(i1,i2,i3,uc)-u(i1-
     & 1,i2,i3,uc))/(dr(0))
                                   ! u*vx
                                   agu(uc,vc)= au*(u(i1,i2,i3,vc)-u(i1-
     & 1,i2,i3,vc))/(dr(0))
                                else
                                   ! u*ux
                                   agu(uc,uc)= au*(u(i1+1,i2,i3,uc)-u(
     & i1,i2,i3,uc))/(dr(0))
                                   ! u*vx
                                   agu(uc,vc)= au*(u(i1+1,i2,i3,vc)-u(
     & i1,i2,i3,vc))/(dr(0))
                                end if
                                au = rsxy(i1,i2,i3,1,0)*uu(i1,i2,i3,uc)
     & +rsxy(i1,i2,i3,1,1)*uu(i1,i2,i3,vc)
                                if( au.gt.0. )then
                                   ! v*uy
                                   agu(vc,uc)= au*(u(i1,i2,i3,uc)-u(i1,
     & i2-1,i3,uc))/(dr(1))
                                   ! v*vy
                                   agu(vc,vc)= au*(u(i1,i2,i3,vc)-u(i1,
     & i2-1,i3,vc))/(dr(1))
                                else
                                   ! v*uy
                                   agu(vc,uc)= au*(u(i1,i2+1,i3,uc)-u(
     & i1,i2,i3,uc))/(dr(1))
                                   ! v*vy
                                   agu(vc,vc)= au*(u(i1,i2+1,i3,vc)-u(
     & i1,i2,i3,vc))/(dr(1))
                                end if
                             elseif( upwindOrder.eq.2 )then
                                au = rsxy(i1,i2,i3,0,0)*uu(i1,i2,i3,uc)
     & +rsxy(i1,i2,i3,0,1)*uu(i1,i2,i3,vc)
                                if( au.gt.0. )then
                                 ! u*ux = u*D-x(u)
                                   agu(uc,uc)= au*(3.*u(i1,i2,i3,uc)-
     & 4.*u(i1-1,i2,i3,uc)+u(i1-2,i2,i3,uc))/(2.*dr(0))
                                 ! u*vx = u*D-x(v)
                                   agu(uc,vc)= au*(3.*u(i1,i2,i3,vc)-
     & 4.*u(i1-1,i2,i3,vc)+u(i1-2,i2,i3,vc))/(2.*dr(0))
                                else
                                 ! u*ux = u*D+x(u)
                                   agu(uc,uc)= au*(-u(i1+2,i2,i3,uc)+
     & 4.*u(i1+1,i2,i3,uc)-3.*u(i1,i2,i3,uc))/(2.*dr(0))
                                 ! u*vx = u*D+x(v)
                                   agu(uc,vc)= au*(-u(i1+2,i2,i3,vc)+
     & 4.*u(i1+1,i2,i3,vc)-3.*u(i1,i2,i3,vc))/(2.*dr(0))
                                end if
                                au = rsxy(i1,i2,i3,1,0)*uu(i1,i2,i3,uc)
     & +rsxy(i1,i2,i3,1,1)*uu(i1,i2,i3,vc)
                                if( au.gt.0. )then
                                 ! v*uy = v*D-y(u)
                                   agu(vc,uc)= au*(3.*u(i1,i2,i3,uc)-
     & 4.*u(i1,i2-1,i3,uc)+u(i1,i2-2,i3,uc))/(2.*dr(1))
                                 ! v*vy = v*D-y(v)
                                   agu(vc,vc)= au*(3.*u(i1,i2,i3,vc)-
     & 4.*u(i1,i2-1,i3,vc)+u(i1,i2-2,i3,vc))/(2.*dr(1))
                                else
                                 ! v*uy = v*D+y(u)
                                   agu(vc,uc)= au*(-u(i1,i2+2,i3,uc)+
     & 4.*u(i1,i2+1,i3,uc)-3.*u(i1,i2,i3,uc))/(2.*dr(1))
                                 ! v*vy = v*D+y(v) 
                                   agu(vc,vc)= au*(-u(i1,i2+2,i3,vc)+
     & 4.*u(i1,i2+1,i3,vc)-3.*u(i1,i2,i3,vc))/(2.*dr(1))
                                end if
                             elseif( upwindOrder.eq.3 )then
                                au = rsxy(i1,i2,i3,0,0)*uu(i1,i2,i3,uc)
     & +rsxy(i1,i2,i3,0,1)*uu(i1,i2,i3,vc)
                                if( au.gt.0. )then
                               ! u*ux = u*D-x(u)
                                   agu(uc,uc)= au*(2.*u(i1+1,i2,i3,uc)+
     & 3.*u(i1,i2,i3,uc)-6.*u(i1-1,i2,i3,uc)+u(i1-2,i2,i3,uc))/(6.*dr(
     & 0))
                               ! u*vx = u*D-x(v)
                                   agu(uc,vc)= au*(2.*u(i1+1,i2,i3,vc)+
     & 3.*u(i1,i2,i3,vc)-6.*u(i1-1,i2,i3,vc)+u(i1-2,i2,i3,vc))/(6.*dr(
     & 0))
                                else
                               ! u*ux = u*D+x(u)
                                   agu(uc,uc)= au*(-u(i1+2,i2,i3,uc)+
     & 6.*u(i1+1,i2,i3,uc)-3.*u(i1,i2,i3,uc)-2.*u(i1-1,i2,i3,uc))/(6.*
     & dr(0))
                               ! u*vx = u*D+x(v)
                                   agu(uc,vc)= au*(-u(i1+2,i2,i3,vc)+
     & 6.*u(i1+1,i2,i3,vc)-3.*u(i1,i2,i3,vc)-2.*u(i1-1,i2,i3,vc))/(6.*
     & dr(0))
                                end if
                                au = rsxy(i1,i2,i3,1,0)*uu(i1,i2,i3,uc)
     & +rsxy(i1,i2,i3,1,1)*uu(i1,i2,i3,vc)
                                if( au.gt.0. )then
                               ! v*uy = v*D-y(u)
                                   agu(vc,uc)= au*(2.*u(i1,i2+1,i3,uc)+
     & 3.*u(i1,i2,i3,uc)-6.*u(i1,i2-1,i3,uc)+u(i1,i2-2,i3,uc))/(6.*dr(
     & 1))
                               ! v*vy = v*D-y(v)
                                   agu(vc,vc)= au*(2.*u(i1,i2+1,i3,vc)+
     & 3.*u(i1,i2,i3,vc)-6.*u(i1,i2-1,i3,vc)+u(i1,i2-2,i3,vc))/(6.*dr(
     & 1))
                                else
                               ! v*uy = v*D+y(u)
                                   agu(vc,uc)= au*(-u(i1,i2+2,i3,uc)+
     & 6.*u(i1,i2+1,i3,uc)-3.*u(i1,i2,i3,uc)-2.*u(i1,i2-1,i3,uc))/(6.*
     & dr(1))
                               ! v*vy = v*D+y(v) 
                                   agu(vc,vc)= au*(-u(i1,i2+2,i3,vc)+
     & 6.*u(i1,i2+1,i3,vc)-3.*u(i1,i2,i3,vc)-2.*u(i1,i2-1,i3,vc))/(6.*
     & dr(1))
                                end if
                             end if
                       !
                       !      
                     ! #If ("UPWIND" == "CENTERED") 
                     !  ! -- centered advection ---
                     !  ! write(*,'(" getAdvection -- centered")')
                     !  #If "2" == "2"
                     !   ! -- 2D --
                     !   agu(uc,uc)=UU(uc)*UX(uc)
                     !   agu(vc,uc)=UU(vc)*UY(uc)
                     !   agu(uc,vc)=UU(uc)*UX(vc)
                     !   agu(vc,vc)=UU(vc)*UY(vc)
                     !  #Else
                     !   ! -- 3D -- *check me*
                     !   agu(uc,uc)=UU(uc)*UX(uc)
                     !   agu(vc,uc)=UU(vc)*UY(uc)
                     !   agu(wc,uc)=UU(wc)*UY(uc)
                     !   agu(uc,vc)=UU(uc)*UX(vc)
                     !   agu(vc,vc)=UU(vc)*UY(vc)
                     !   agu(wc,vc)=UU(wc)*UY(vc)
                     !   agu(uc,wc)=UU(uc)*UX(wc)
                     !   agu(vc,wc)=UU(vc)*UY(wc)
                     !   agu(wc,wc)=UU(wc)*UY(wc)
                     !  #End
                     ! #Elif "UPWIND" == "UPWIND" 
                     !   ! --- upwind scheme ---
                     !   ! for testing output this next message:
                     !   if( t.le. 0. )then
                     !     write(*,'(" getAdvection upwind scheme (7)")') 
                     !   end if
                     !   getUpwindAdvection(u,i1,i2,i3,NONE,2,4,curvilinear, agu)
                     ! #Elif "UPWIND" == "BWENO" 
                     !   ! --- Bweno scheme ---
                     !   ! for testing output this message:
                     !   if( t.le. 0. )then
                     !      write(*,'(" getAdvection BWENO scheme (7)")') 
                     !   end if
                     !   getBwenoAdvection(u,i1,i2,i3,NONE,2,4,curvilinear, agu)
                     ! #Else
                     !   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
                     !   stop 999
                     ! #End
                       ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-ux42(
     & i1,i2,i3,pc)
                       ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-uy42(
     & i1,i2,i3,pc)
                   ! include implicit terms - diffusion
                    uti(i1,i2,i3,uc)= nu*ulaplacian42(i1,i2,i3,uc)
                    uti(i1,i2,i3,vc)= nu*ulaplacian42(i1,i2,i3,vc)
                end do
                end do
                end do
               end if
             else if( implicitOption.eq.doNotComputeImplicitTerms )then
               if( useWhereMask.ne.0 )then
                do i3=n3a,n3b
                do i2=n2a,n2b
                do i1=n1a,n1b
                 if( mask(i1,i2,i3).gt.0 )then
                    ! explicit terms only, no diffusion
                      ! upwind approximation/bweno
                        ! --- upwind scheme ---
                        ! for testing output this next message:
                        ! if( t.le. 0. )then
                        !    write(*,'(" getAdvection upwind scheme (7)")') 
                        ! end if
                          ! write(*,'(" we in HERE")' )
                          ! stop 7171
                          ! --- CURVILINEAR GRID ---
                              if( upwindOrder.eq.1 )then
                                 au   = rsxy(i1,i2,i3,0,0)*uu(i1,i2,i3,
     & uc)+rsxy(i1,i2,i3,0,1)*uu(i1,i2,i3,vc)
                                 if( au.gt.0. )then
                                    ! u*ux
                                    agu(uc,uc)= au*(u(i1,i2,i3,uc)-u(
     & i1-1,i2,i3,uc))/(dr(0))
                                    ! u*vx
                                    agu(uc,vc)= au*(u(i1,i2,i3,vc)-u(
     & i1-1,i2,i3,vc))/(dr(0))
                                 else
                                    ! u*ux
                                    agu(uc,uc)= au*(u(i1+1,i2,i3,uc)-u(
     & i1,i2,i3,uc))/(dr(0))
                                    ! u*vx
                                    agu(uc,vc)= au*(u(i1+1,i2,i3,vc)-u(
     & i1,i2,i3,vc))/(dr(0))
                                 end if
                                 au = rsxy(i1,i2,i3,1,0)*uu(i1,i2,i3,
     & uc)+rsxy(i1,i2,i3,1,1)*uu(i1,i2,i3,vc)
                                 if( au.gt.0. )then
                                    ! v*uy
                                    agu(vc,uc)= au*(u(i1,i2,i3,uc)-u(
     & i1,i2-1,i3,uc))/(dr(1))
                                    ! v*vy
                                    agu(vc,vc)= au*(u(i1,i2,i3,vc)-u(
     & i1,i2-1,i3,vc))/(dr(1))
                                 else
                                    ! v*uy
                                    agu(vc,uc)= au*(u(i1,i2+1,i3,uc)-u(
     & i1,i2,i3,uc))/(dr(1))
                                    ! v*vy
                                    agu(vc,vc)= au*(u(i1,i2+1,i3,vc)-u(
     & i1,i2,i3,vc))/(dr(1))
                                 end if
                              elseif( upwindOrder.eq.2 )then
                                 au = rsxy(i1,i2,i3,0,0)*uu(i1,i2,i3,
     & uc)+rsxy(i1,i2,i3,0,1)*uu(i1,i2,i3,vc)
                                 if( au.gt.0. )then
                                  ! u*ux = u*D-x(u)
                                    agu(uc,uc)= au*(3.*u(i1,i2,i3,uc)-
     & 4.*u(i1-1,i2,i3,uc)+u(i1-2,i2,i3,uc))/(2.*dr(0))
                                  ! u*vx = u*D-x(v)
                                    agu(uc,vc)= au*(3.*u(i1,i2,i3,vc)-
     & 4.*u(i1-1,i2,i3,vc)+u(i1-2,i2,i3,vc))/(2.*dr(0))
                                 else
                                  ! u*ux = u*D+x(u)
                                    agu(uc,uc)= au*(-u(i1+2,i2,i3,uc)+
     & 4.*u(i1+1,i2,i3,uc)-3.*u(i1,i2,i3,uc))/(2.*dr(0))
                                  ! u*vx = u*D+x(v)
                                    agu(uc,vc)= au*(-u(i1+2,i2,i3,vc)+
     & 4.*u(i1+1,i2,i3,vc)-3.*u(i1,i2,i3,vc))/(2.*dr(0))
                                 end if
                                 au = rsxy(i1,i2,i3,1,0)*uu(i1,i2,i3,
     & uc)+rsxy(i1,i2,i3,1,1)*uu(i1,i2,i3,vc)
                                 if( au.gt.0. )then
                                  ! v*uy = v*D-y(u)
                                    agu(vc,uc)= au*(3.*u(i1,i2,i3,uc)-
     & 4.*u(i1,i2-1,i3,uc)+u(i1,i2-2,i3,uc))/(2.*dr(1))
                                  ! v*vy = v*D-y(v)
                                    agu(vc,vc)= au*(3.*u(i1,i2,i3,vc)-
     & 4.*u(i1,i2-1,i3,vc)+u(i1,i2-2,i3,vc))/(2.*dr(1))
                                 else
                                  ! v*uy = v*D+y(u)
                                    agu(vc,uc)= au*(-u(i1,i2+2,i3,uc)+
     & 4.*u(i1,i2+1,i3,uc)-3.*u(i1,i2,i3,uc))/(2.*dr(1))
                                  ! v*vy = v*D+y(v) 
                                    agu(vc,vc)= au*(-u(i1,i2+2,i3,vc)+
     & 4.*u(i1,i2+1,i3,vc)-3.*u(i1,i2,i3,vc))/(2.*dr(1))
                                 end if
                              elseif( upwindOrder.eq.3 )then
                                 au = rsxy(i1,i2,i3,0,0)*uu(i1,i2,i3,
     & uc)+rsxy(i1,i2,i3,0,1)*uu(i1,i2,i3,vc)
                                 if( au.gt.0. )then
                                ! u*ux = u*D-x(u)
                                    agu(uc,uc)= au*(2.*u(i1+1,i2,i3,uc)
     & +3.*u(i1,i2,i3,uc)-6.*u(i1-1,i2,i3,uc)+u(i1-2,i2,i3,uc))/(6.*
     & dr(0))
                                ! u*vx = u*D-x(v)
                                    agu(uc,vc)= au*(2.*u(i1+1,i2,i3,vc)
     & +3.*u(i1,i2,i3,vc)-6.*u(i1-1,i2,i3,vc)+u(i1-2,i2,i3,vc))/(6.*
     & dr(0))
                                 else
                                ! u*ux = u*D+x(u)
                                    agu(uc,uc)= au*(-u(i1+2,i2,i3,uc)+
     & 6.*u(i1+1,i2,i3,uc)-3.*u(i1,i2,i3,uc)-2.*u(i1-1,i2,i3,uc))/(6.*
     & dr(0))
                                ! u*vx = u*D+x(v)
                                    agu(uc,vc)= au*(-u(i1+2,i2,i3,vc)+
     & 6.*u(i1+1,i2,i3,vc)-3.*u(i1,i2,i3,vc)-2.*u(i1-1,i2,i3,vc))/(6.*
     & dr(0))
                                 end if
                                 au = rsxy(i1,i2,i3,1,0)*uu(i1,i2,i3,
     & uc)+rsxy(i1,i2,i3,1,1)*uu(i1,i2,i3,vc)
                                 if( au.gt.0. )then
                                ! v*uy = v*D-y(u)
                                    agu(vc,uc)= au*(2.*u(i1,i2+1,i3,uc)
     & +3.*u(i1,i2,i3,uc)-6.*u(i1,i2-1,i3,uc)+u(i1,i2-2,i3,uc))/(6.*
     & dr(1))
                                ! v*vy = v*D-y(v)
                                    agu(vc,vc)= au*(2.*u(i1,i2+1,i3,vc)
     & +3.*u(i1,i2,i3,vc)-6.*u(i1,i2-1,i3,vc)+u(i1,i2-2,i3,vc))/(6.*
     & dr(1))
                                 else
                                ! v*uy = v*D+y(u)
                                    agu(vc,uc)= au*(-u(i1,i2+2,i3,uc)+
     & 6.*u(i1,i2+1,i3,uc)-3.*u(i1,i2,i3,uc)-2.*u(i1,i2-1,i3,uc))/(6.*
     & dr(1))
                                ! v*vy = v*D+y(v) 
                                    agu(vc,vc)= au*(-u(i1,i2+2,i3,vc)+
     & 6.*u(i1,i2+1,i3,vc)-3.*u(i1,i2,i3,vc)-2.*u(i1,i2-1,i3,vc))/(6.*
     & dr(1))
                                 end if
                              end if
                        !
                        !      
                      ! #If ("UPWIND" == "CENTERED") 
                      !  ! -- centered advection ---
                      !  ! write(*,'(" getAdvection -- centered")')
                      !  #If "2" == "2"
                      !   ! -- 2D --
                      !   agu(uc,uc)=UU(uc)*UX(uc)
                      !   agu(vc,uc)=UU(vc)*UY(uc)
                      !   agu(uc,vc)=UU(uc)*UX(vc)
                      !   agu(vc,vc)=UU(vc)*UY(vc)
                      !  #Else
                      !   ! -- 3D -- *check me*
                      !   agu(uc,uc)=UU(uc)*UX(uc)
                      !   agu(vc,uc)=UU(vc)*UY(uc)
                      !   agu(wc,uc)=UU(wc)*UY(uc)
                      !   agu(uc,vc)=UU(uc)*UX(vc)
                      !   agu(vc,vc)=UU(vc)*UY(vc)
                      !   agu(wc,vc)=UU(wc)*UY(vc)
                      !   agu(uc,wc)=UU(uc)*UX(wc)
                      !   agu(vc,wc)=UU(vc)*UY(wc)
                      !   agu(wc,wc)=UU(wc)*UY(wc)
                      !  #End
                      ! #Elif "UPWIND" == "UPWIND" 
                      !   ! --- upwind scheme ---
                      !   ! for testing output this next message:
                      !   if( t.le. 0. )then
                      !     write(*,'(" getAdvection upwind scheme (7)")') 
                      !   end if
                      !   getUpwindAdvection(u,i1,i2,i3,NONE,2,4,curvilinear, agu)
                      ! #Elif "UPWIND" == "BWENO" 
                      !   ! --- Bweno scheme ---
                      !   ! for testing output this message:
                      !   if( t.le. 0. )then
                      !      write(*,'(" getAdvection BWENO scheme (7)")') 
                      !   end if
                      !   getBwenoAdvection(u,i1,i2,i3,NONE,2,4,curvilinear, agu)
                      ! #Else
                      !   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
                      !   stop 999
                      ! #End
                        ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-ux42(
     & i1,i2,i3,pc)
                        ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-uy42(
     & i1,i2,i3,pc)
                 end if
                end do
                end do
                end do
               else
                do i3=n3a,n3b
                do i2=n2a,n2b
                do i1=n1a,n1b
                   ! explicit terms only, no diffusion
                     ! upwind approximation/bweno
                       ! --- upwind scheme ---
                       ! for testing output this next message:
                       ! if( t.le. 0. )then
                       !    write(*,'(" getAdvection upwind scheme (7)")') 
                       ! end if
                         ! write(*,'(" we in HERE")' )
                         ! stop 7171
                         ! --- CURVILINEAR GRID ---
                             if( upwindOrder.eq.1 )then
                                au   = rsxy(i1,i2,i3,0,0)*uu(i1,i2,i3,
     & uc)+rsxy(i1,i2,i3,0,1)*uu(i1,i2,i3,vc)
                                if( au.gt.0. )then
                                   ! u*ux
                                   agu(uc,uc)= au*(u(i1,i2,i3,uc)-u(i1-
     & 1,i2,i3,uc))/(dr(0))
                                   ! u*vx
                                   agu(uc,vc)= au*(u(i1,i2,i3,vc)-u(i1-
     & 1,i2,i3,vc))/(dr(0))
                                else
                                   ! u*ux
                                   agu(uc,uc)= au*(u(i1+1,i2,i3,uc)-u(
     & i1,i2,i3,uc))/(dr(0))
                                   ! u*vx
                                   agu(uc,vc)= au*(u(i1+1,i2,i3,vc)-u(
     & i1,i2,i3,vc))/(dr(0))
                                end if
                                au = rsxy(i1,i2,i3,1,0)*uu(i1,i2,i3,uc)
     & +rsxy(i1,i2,i3,1,1)*uu(i1,i2,i3,vc)
                                if( au.gt.0. )then
                                   ! v*uy
                                   agu(vc,uc)= au*(u(i1,i2,i3,uc)-u(i1,
     & i2-1,i3,uc))/(dr(1))
                                   ! v*vy
                                   agu(vc,vc)= au*(u(i1,i2,i3,vc)-u(i1,
     & i2-1,i3,vc))/(dr(1))
                                else
                                   ! v*uy
                                   agu(vc,uc)= au*(u(i1,i2+1,i3,uc)-u(
     & i1,i2,i3,uc))/(dr(1))
                                   ! v*vy
                                   agu(vc,vc)= au*(u(i1,i2+1,i3,vc)-u(
     & i1,i2,i3,vc))/(dr(1))
                                end if
                             elseif( upwindOrder.eq.2 )then
                                au = rsxy(i1,i2,i3,0,0)*uu(i1,i2,i3,uc)
     & +rsxy(i1,i2,i3,0,1)*uu(i1,i2,i3,vc)
                                if( au.gt.0. )then
                                 ! u*ux = u*D-x(u)
                                   agu(uc,uc)= au*(3.*u(i1,i2,i3,uc)-
     & 4.*u(i1-1,i2,i3,uc)+u(i1-2,i2,i3,uc))/(2.*dr(0))
                                 ! u*vx = u*D-x(v)
                                   agu(uc,vc)= au*(3.*u(i1,i2,i3,vc)-
     & 4.*u(i1-1,i2,i3,vc)+u(i1-2,i2,i3,vc))/(2.*dr(0))
                                else
                                 ! u*ux = u*D+x(u)
                                   agu(uc,uc)= au*(-u(i1+2,i2,i3,uc)+
     & 4.*u(i1+1,i2,i3,uc)-3.*u(i1,i2,i3,uc))/(2.*dr(0))
                                 ! u*vx = u*D+x(v)
                                   agu(uc,vc)= au*(-u(i1+2,i2,i3,vc)+
     & 4.*u(i1+1,i2,i3,vc)-3.*u(i1,i2,i3,vc))/(2.*dr(0))
                                end if
                                au = rsxy(i1,i2,i3,1,0)*uu(i1,i2,i3,uc)
     & +rsxy(i1,i2,i3,1,1)*uu(i1,i2,i3,vc)
                                if( au.gt.0. )then
                                 ! v*uy = v*D-y(u)
                                   agu(vc,uc)= au*(3.*u(i1,i2,i3,uc)-
     & 4.*u(i1,i2-1,i3,uc)+u(i1,i2-2,i3,uc))/(2.*dr(1))
                                 ! v*vy = v*D-y(v)
                                   agu(vc,vc)= au*(3.*u(i1,i2,i3,vc)-
     & 4.*u(i1,i2-1,i3,vc)+u(i1,i2-2,i3,vc))/(2.*dr(1))
                                else
                                 ! v*uy = v*D+y(u)
                                   agu(vc,uc)= au*(-u(i1,i2+2,i3,uc)+
     & 4.*u(i1,i2+1,i3,uc)-3.*u(i1,i2,i3,uc))/(2.*dr(1))
                                 ! v*vy = v*D+y(v) 
                                   agu(vc,vc)= au*(-u(i1,i2+2,i3,vc)+
     & 4.*u(i1,i2+1,i3,vc)-3.*u(i1,i2,i3,vc))/(2.*dr(1))
                                end if
                             elseif( upwindOrder.eq.3 )then
                                au = rsxy(i1,i2,i3,0,0)*uu(i1,i2,i3,uc)
     & +rsxy(i1,i2,i3,0,1)*uu(i1,i2,i3,vc)
                                if( au.gt.0. )then
                               ! u*ux = u*D-x(u)
                                   agu(uc,uc)= au*(2.*u(i1+1,i2,i3,uc)+
     & 3.*u(i1,i2,i3,uc)-6.*u(i1-1,i2,i3,uc)+u(i1-2,i2,i3,uc))/(6.*dr(
     & 0))
                               ! u*vx = u*D-x(v)
                                   agu(uc,vc)= au*(2.*u(i1+1,i2,i3,vc)+
     & 3.*u(i1,i2,i3,vc)-6.*u(i1-1,i2,i3,vc)+u(i1-2,i2,i3,vc))/(6.*dr(
     & 0))
                                else
                               ! u*ux = u*D+x(u)
                                   agu(uc,uc)= au*(-u(i1+2,i2,i3,uc)+
     & 6.*u(i1+1,i2,i3,uc)-3.*u(i1,i2,i3,uc)-2.*u(i1-1,i2,i3,uc))/(6.*
     & dr(0))
                               ! u*vx = u*D+x(v)
                                   agu(uc,vc)= au*(-u(i1+2,i2,i3,vc)+
     & 6.*u(i1+1,i2,i3,vc)-3.*u(i1,i2,i3,vc)-2.*u(i1-1,i2,i3,vc))/(6.*
     & dr(0))
                                end if
                                au = rsxy(i1,i2,i3,1,0)*uu(i1,i2,i3,uc)
     & +rsxy(i1,i2,i3,1,1)*uu(i1,i2,i3,vc)
                                if( au.gt.0. )then
                               ! v*uy = v*D-y(u)
                                   agu(vc,uc)= au*(2.*u(i1,i2+1,i3,uc)+
     & 3.*u(i1,i2,i3,uc)-6.*u(i1,i2-1,i3,uc)+u(i1,i2-2,i3,uc))/(6.*dr(
     & 1))
                               ! v*vy = v*D-y(v)
                                   agu(vc,vc)= au*(2.*u(i1,i2+1,i3,vc)+
     & 3.*u(i1,i2,i3,vc)-6.*u(i1,i2-1,i3,vc)+u(i1,i2-2,i3,vc))/(6.*dr(
     & 1))
                                else
                               ! v*uy = v*D+y(u)
                                   agu(vc,uc)= au*(-u(i1,i2+2,i3,uc)+
     & 6.*u(i1,i2+1,i3,uc)-3.*u(i1,i2,i3,uc)-2.*u(i1,i2-1,i3,uc))/(6.*
     & dr(1))
                               ! v*vy = v*D+y(v) 
                                   agu(vc,vc)= au*(-u(i1,i2+2,i3,vc)+
     & 6.*u(i1,i2+1,i3,vc)-3.*u(i1,i2,i3,vc)-2.*u(i1,i2-1,i3,vc))/(6.*
     & dr(1))
                                end if
                             end if
                       !
                       !      
                     ! #If ("UPWIND" == "CENTERED") 
                     !  ! -- centered advection ---
                     !  ! write(*,'(" getAdvection -- centered")')
                     !  #If "2" == "2"
                     !   ! -- 2D --
                     !   agu(uc,uc)=UU(uc)*UX(uc)
                     !   agu(vc,uc)=UU(vc)*UY(uc)
                     !   agu(uc,vc)=UU(uc)*UX(vc)
                     !   agu(vc,vc)=UU(vc)*UY(vc)
                     !  #Else
                     !   ! -- 3D -- *check me*
                     !   agu(uc,uc)=UU(uc)*UX(uc)
                     !   agu(vc,uc)=UU(vc)*UY(uc)
                     !   agu(wc,uc)=UU(wc)*UY(uc)
                     !   agu(uc,vc)=UU(uc)*UX(vc)
                     !   agu(vc,vc)=UU(vc)*UY(vc)
                     !   agu(wc,vc)=UU(wc)*UY(vc)
                     !   agu(uc,wc)=UU(uc)*UX(wc)
                     !   agu(vc,wc)=UU(vc)*UY(wc)
                     !   agu(wc,wc)=UU(wc)*UY(wc)
                     !  #End
                     ! #Elif "UPWIND" == "UPWIND" 
                     !   ! --- upwind scheme ---
                     !   ! for testing output this next message:
                     !   if( t.le. 0. )then
                     !     write(*,'(" getAdvection upwind scheme (7)")') 
                     !   end if
                     !   getUpwindAdvection(u,i1,i2,i3,NONE,2,4,curvilinear, agu)
                     ! #Elif "UPWIND" == "BWENO" 
                     !   ! --- Bweno scheme ---
                     !   ! for testing output this message:
                     !   if( t.le. 0. )then
                     !      write(*,'(" getAdvection BWENO scheme (7)")') 
                     !   end if
                     !   getBwenoAdvection(u,i1,i2,i3,NONE,2,4,curvilinear, agu)
                     ! #Else
                     !   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
                     !   stop 999
                     ! #End
                       ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-ux42(
     & i1,i2,i3,pc)
                       ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-uy42(
     & i1,i2,i3,pc)
                end do
                end do
                end do
               end if
             else
              write(*,*)'insdt: Unknown implicitOption=',implicitOption
              stop 6
             end if  ! end implicitOption
            else if( advectionOption.eq.bwenoAdvection )then
              ! --- bweno ---
             if( implicitOption .eq.computeImplicitTermsSeparately )
     & then
               if( useWhereMask.ne.0 )then
                do i3=n3a,n3b
                do i2=n2a,n2b
                do i1=n1a,n1b
                 if( mask(i1,i2,i3).gt.0 )then
                    ! explicit terms only, no diffusion
                      ! upwind approximation/bweno
                        ! --- Bweno scheme ---
                        ! for testing output this message:
                        ! if( t.le. 0. )then
                        !    write(*,'(" getAdvection BWENO scheme (7)")') 
                        ! end if
                           gvU = u(i1,i2,i3,uc) - uu(i1,i2,i3,uc)
                           gvV = u(i1,i2,i3,vc) - uu(i1,i2,i3,vc)
                           gvW = u(i1,i2,i3,wc) - uu(i1,i2,i3,wc)
                          !bweno curvulinear
                          orderOneFix = 1
                          augmentPlot = 0
                          if( upwindOrder.eq.4 )then
                             ! u*ux-------------------------------------------------------------------------------
                             ! au = rsxy(i1,i2,i3,0,0)*u(i1,i2,i3,uc)+rsxy(i1,i2,i3,0,1)*u(i1,i2,i3,vc) + rsxy(i1,i2,i3,0,2)*u(i1,i2,i3,wc)
                             ! aum   = rsxy(i1-1,i2,i3,0,0)*u(i1-1,i2,i3,uc)+rsxy(i1-1,i2,i3,0,1)*u(i1-1,i2,i3,vc)+ rsxy(i1-1,i2,i3,0,2)*u(i1-1,i2,i3,wc)
                             ! aup   = rsxy(i1+1,i2,i3,0,0)*u(i1+1,i2,i3,uc)+rsxy(i1+1,i2,i3,0,1)*u(i1+1,i2,i3,vc)+ rsxy(i1+1,i2,i3,0,2)*u(i1+1,i2,i3,wc)
                             au    = rsxy(i1  ,i2,i3,0,0)*uu(i1  ,i2,
     & i3,uc) + rsxy(i1  ,i2,i3,0,1) * uu(i1  ,i2,i3,vc)
                             aum   = rsxy(i1-1,i2,i3,0,0)*(u(i1-1,i2,
     & i3,uc) - gvU ) + rsxy(i1-1,i2,i3,0,1) * (u(i1-1,i2,i3,vc)- gvV 
     & )
                             aup   = rsxy(i1+1,i2,i3,0,0)*(u(i1+1,i2,
     & i3,uc) - gvU ) + rsxy(i1+1,i2,i3,0,1) * (u(i1+1,i2,i3,vc)- gvV 
     & )
                             s1=1
                             s2=0
                             s3=0
                             var = uc
                             vard= uc
                             drl = dr(0)
                               Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*u(
     & i1   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
                               Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    u(
     & i1-s1,i2-s2,i3-s3,var)
                               Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) -
     &  2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
                               Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) -
     &     u(i1   ,i2   ,i3   ,var)
                               Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,
     & i2-s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
                               Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var)
                               Amr   = Apl
                               Bmr   = Bpl
                               betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 
     & 1./4.*Bpl**2
                               betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 
     & 1./4.*Bpr**2
                               betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 
     & 1./4.*Bml**2
                               betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 
     & 1./4.*Bmr**2
                               ep     = 1e-15
                               ap1    = (1./2.)/(ep + betapl)**(4)
                               ap2    = (1./2.)/(ep + betapr)**(4)
                               am1    = (1./2.)/(ep + betaml)**(4)
                               am2    = (1./2.)/(ep + betamr)**(4)
                               !calculate weights
                               wp1 = ap1/(ap1+ap2)
                               wp2 = ap2/(ap1+ap2)
                               wm1 = am1/(am1+am2)
                               wm2 = am2/(am1+am2)
                               !mapping of the weights
                               bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                               bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                               bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                               bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                               wplw = bp1/(bp1+bp2)
                               wprw = bp2/(bp1+bp2)
                               wmlw = bm1/(bm1+bm2)
                               wmrw = bm2/(bm1+bm2)
                               !get face value of u
                               up = au
                               um = au
                               !picking side for the weights
                               if (up.ge.0) then
                                  wpl = max(wplw,wprw)
                                  wpr = min(wplw,wprw)
                               else
                                  wpl = min(wplw,wprw)
                                  wpr = max(wplw,wprw)
                               end if
                               if (um.ge.0) then
                                  wml = max(wmlw,wmrw)
                                  wmr = min(wmlw,wmrw)
                               else
                                  wml = min(wmlw,wmrw)
                                  wmr = max(wmlw,wmrw)
                               end if
                               expAd = 0.
                               if (orderOneFix.eq.1) then
                                  ! if (augmentPlot.eq.1) then
                                  !    agu(uc,13)=0
                                  ! end if
                                  if (aum*aup .le. 0) then
                                     aumax = max(abs(au),abs(aup),abs(
     & aum))
                                     expAd = (1./12.)*aumax*( u(i1-2*
     & s1,i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(
     & i1,i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*
     & s2,i3+2*s3,var))
                                     wpl = 0.5
                                     wpr = 0.5
                                     wml = 0.5
                                     wmr = 0.5
                                     ! plot au at where fix is turned on
                                     ! if (augmentPlot.eq.1) then
                                     !    agu(uc,13)=au
                                     ! end if
                                  end if
                               end if
                               ! wpl = 0.5
                               ! wpr = 0.5
                               ! wml = 0.5
                               ! wmr = 0.5
                               Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,var)
     &  + 5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
                               Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 5.*
     & u(i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
                               Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*
     & s3,var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
                               Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,
     & var) + 5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
                               Fp  = wpl*Fpl + wpr*Fpr
                               Fm  = wml*Fml + wmr*Fmr
                                  agu(vard,var)= au*(Fp - Fm)/drl + 
     & expAd/drl
                             !     if ( i1 .eq. 14 .and. i2.eq.16 ) then
                                     !     if ( i1 .eq. 14) then
                                     ! write(*,'(" betapl,betapr,betaml,betamr = ",4F10.8)')  betapl,betapr,betaml,betamr
                                     ! write(*,'(" ap1,ap2,am1,am2 = ",4e10.8)')  ap1,ap2,am1,am2
                                     ! write(*,'(" wp1,wpr,wm1,wm2 = ",4F10.8)')  wp1,wp2,wm1,wm2
                                     ! write(*,'(" bp1,bp2,bm1,bm2 = ",4F10.8)')  bp1,bp2,bm1,bm2
                                     ! write(*,'(" wplw,wprw,wmlw,wmrw = ",4F10.8)')  wplw,wprw,wmlw,wmrw
                                     ! write(*,'(" wpl,wpr,wml,wmr = ",4F10.8)')  wpl,wpr,wml,wmr
                                     ! write(*,'(" agu = ",F10.8)') agu(vard,var)
                             !     end if
                             !u*vx
                             s1=1
                             s2=0
                             s3=0
                             var = vc
                             vard= uc
                             drl = dr(0)
                               Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*u(
     & i1   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
                               Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    u(
     & i1-s1,i2-s2,i3-s3,var)
                               Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) -
     &  2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
                               Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) -
     &     u(i1   ,i2   ,i3   ,var)
                               Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,
     & i2-s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
                               Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var)
                               Amr   = Apl
                               Bmr   = Bpl
                               betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 
     & 1./4.*Bpl**2
                               betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 
     & 1./4.*Bpr**2
                               betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 
     & 1./4.*Bml**2
                               betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 
     & 1./4.*Bmr**2
                               ep     = 1e-15
                               ap1    = (1./2.)/(ep + betapl)**(4)
                               ap2    = (1./2.)/(ep + betapr)**(4)
                               am1    = (1./2.)/(ep + betaml)**(4)
                               am2    = (1./2.)/(ep + betamr)**(4)
                               !calculate weights
                               wp1 = ap1/(ap1+ap2)
                               wp2 = ap2/(ap1+ap2)
                               wm1 = am1/(am1+am2)
                               wm2 = am2/(am1+am2)
                               !mapping of the weights
                               bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                               bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                               bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                               bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                               wplw = bp1/(bp1+bp2)
                               wprw = bp2/(bp1+bp2)
                               wmlw = bm1/(bm1+bm2)
                               wmrw = bm2/(bm1+bm2)
                               !get face value of u
                               up = au
                               um = au
                               !picking side for the weights
                               if (up.ge.0) then
                                  wpl = max(wplw,wprw)
                                  wpr = min(wplw,wprw)
                               else
                                  wpl = min(wplw,wprw)
                                  wpr = max(wplw,wprw)
                               end if
                               if (um.ge.0) then
                                  wml = max(wmlw,wmrw)
                                  wmr = min(wmlw,wmrw)
                               else
                                  wml = min(wmlw,wmrw)
                                  wmr = max(wmlw,wmrw)
                               end if
                               expAd = 0.
                               if (orderOneFix.eq.1) then
                                  ! if (augmentPlot.eq.1) then
                                  !    agu(uc,13)=0
                                  ! end if
                                  if (aum*aup .le. 0) then
                                     aumax = max(abs(au),abs(aup),abs(
     & aum))
                                     expAd = (1./12.)*aumax*( u(i1-2*
     & s1,i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(
     & i1,i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*
     & s2,i3+2*s3,var))
                                     wpl = 0.5
                                     wpr = 0.5
                                     wml = 0.5
                                     wmr = 0.5
                                     ! plot au at where fix is turned on
                                     ! if (augmentPlot.eq.1) then
                                     !    agu(uc,13)=au
                                     ! end if
                                  end if
                               end if
                               ! wpl = 0.5
                               ! wpr = 0.5
                               ! wml = 0.5
                               ! wmr = 0.5
                               Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,var)
     &  + 5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
                               Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 5.*
     & u(i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
                               Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*
     & s3,var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
                               Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,
     & var) + 5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
                               Fp  = wpl*Fpl + wpr*Fpr
                               Fm  = wml*Fml + wmr*Fmr
                                  agu(vard,var)= au*(Fp - Fm)/drl + 
     & expAd/drl
                             !     if ( i1 .eq. 14 .and. i2.eq.16 ) then
                                     !     if ( i1 .eq. 14) then
                                     ! write(*,'(" betapl,betapr,betaml,betamr = ",4F10.8)')  betapl,betapr,betaml,betamr
                                     ! write(*,'(" ap1,ap2,am1,am2 = ",4e10.8)')  ap1,ap2,am1,am2
                                     ! write(*,'(" wp1,wpr,wm1,wm2 = ",4F10.8)')  wp1,wp2,wm1,wm2
                                     ! write(*,'(" bp1,bp2,bm1,bm2 = ",4F10.8)')  bp1,bp2,bm1,bm2
                                     ! write(*,'(" wplw,wprw,wmlw,wmrw = ",4F10.8)')  wplw,wprw,wmlw,wmrw
                                     ! write(*,'(" wpl,wpr,wml,wmr = ",4F10.8)')  wpl,wpr,wml,wmr
                                     ! write(*,'(" agu = ",F10.8)') agu(vard,var)
                             !     end if
                             !v*uy
                             au    = rsxy(i1,i2  ,i3,1,0)*uu(i1,i2  ,
     & i3,uc)+rsxy(i1,i2,i3,1,1)*uu(i1,i2,i3,vc)
                             aum   = rsxy(i1,i2-1,i3,1,0)*(u(i1,i2-1,
     & i3,uc) - gvU) + rsxy(i1,i2-1,i3,1,1)*(u(i1,i2-1,i3,vc) - gvV)
                             aup   = rsxy(i1,i2+1,i3,1,0)*(u(i1,i2+1,
     & i3,uc) - gvU) + rsxy(i1,i2+1,i3,1,1)*(u(i1,i2+1,i3,vc) - gvV)
                             s1=0
                             s2=1
                             s3=0
                             var = uc
                             vard= vc
                             drl = dr(1)
                               Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*u(
     & i1   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
                               Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    u(
     & i1-s1,i2-s2,i3-s3,var)
                               Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) -
     &  2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
                               Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) -
     &     u(i1   ,i2   ,i3   ,var)
                               Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,
     & i2-s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
                               Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var)
                               Amr   = Apl
                               Bmr   = Bpl
                               betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 
     & 1./4.*Bpl**2
                               betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 
     & 1./4.*Bpr**2
                               betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 
     & 1./4.*Bml**2
                               betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 
     & 1./4.*Bmr**2
                               ep     = 1e-15
                               ap1    = (1./2.)/(ep + betapl)**(4)
                               ap2    = (1./2.)/(ep + betapr)**(4)
                               am1    = (1./2.)/(ep + betaml)**(4)
                               am2    = (1./2.)/(ep + betamr)**(4)
                               !calculate weights
                               wp1 = ap1/(ap1+ap2)
                               wp2 = ap2/(ap1+ap2)
                               wm1 = am1/(am1+am2)
                               wm2 = am2/(am1+am2)
                               !mapping of the weights
                               bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                               bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                               bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                               bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                               wplw = bp1/(bp1+bp2)
                               wprw = bp2/(bp1+bp2)
                               wmlw = bm1/(bm1+bm2)
                               wmrw = bm2/(bm1+bm2)
                               !get face value of u
                               up = au
                               um = au
                               !picking side for the weights
                               if (up.ge.0) then
                                  wpl = max(wplw,wprw)
                                  wpr = min(wplw,wprw)
                               else
                                  wpl = min(wplw,wprw)
                                  wpr = max(wplw,wprw)
                               end if
                               if (um.ge.0) then
                                  wml = max(wmlw,wmrw)
                                  wmr = min(wmlw,wmrw)
                               else
                                  wml = min(wmlw,wmrw)
                                  wmr = max(wmlw,wmrw)
                               end if
                               expAd = 0.
                               if (orderOneFix.eq.1) then
                                  ! if (augmentPlot.eq.1) then
                                  !    agu(uc,13)=0
                                  ! end if
                                  if (aum*aup .le. 0) then
                                     aumax = max(abs(au),abs(aup),abs(
     & aum))
                                     expAd = (1./12.)*aumax*( u(i1-2*
     & s1,i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(
     & i1,i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*
     & s2,i3+2*s3,var))
                                     wpl = 0.5
                                     wpr = 0.5
                                     wml = 0.5
                                     wmr = 0.5
                                     ! plot au at where fix is turned on
                                     ! if (augmentPlot.eq.1) then
                                     !    agu(uc,13)=au
                                     ! end if
                                  end if
                               end if
                               ! wpl = 0.5
                               ! wpr = 0.5
                               ! wml = 0.5
                               ! wmr = 0.5
                               Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,var)
     &  + 5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
                               Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 5.*
     & u(i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
                               Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*
     & s3,var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
                               Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,
     & var) + 5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
                               Fp  = wpl*Fpl + wpr*Fpr
                               Fm  = wml*Fml + wmr*Fmr
                                  agu(vard,var)= au*(Fp - Fm)/drl + 
     & expAd/drl
                             !     if ( i1 .eq. 14 .and. i2.eq.16 ) then
                                     !     if ( i1 .eq. 14) then
                                     ! write(*,'(" betapl,betapr,betaml,betamr = ",4F10.8)')  betapl,betapr,betaml,betamr
                                     ! write(*,'(" ap1,ap2,am1,am2 = ",4e10.8)')  ap1,ap2,am1,am2
                                     ! write(*,'(" wp1,wpr,wm1,wm2 = ",4F10.8)')  wp1,wp2,wm1,wm2
                                     ! write(*,'(" bp1,bp2,bm1,bm2 = ",4F10.8)')  bp1,bp2,bm1,bm2
                                     ! write(*,'(" wplw,wprw,wmlw,wmrw = ",4F10.8)')  wplw,wprw,wmlw,wmrw
                                     ! write(*,'(" wpl,wpr,wml,wmr = ",4F10.8)')  wpl,wpr,wml,wmr
                                     ! write(*,'(" agu = ",F10.8)') agu(vard,var)
                             !     end if
                             !v*vy
                             s1=0
                             s2=1
                             s3=0
                             var = vc
                             vard= vc
                             drl = dr(1)
                               Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*u(
     & i1   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
                               Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    u(
     & i1-s1,i2-s2,i3-s3,var)
                               Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) -
     &  2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
                               Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) -
     &     u(i1   ,i2   ,i3   ,var)
                               Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,
     & i2-s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
                               Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var)
                               Amr   = Apl
                               Bmr   = Bpl
                               betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 
     & 1./4.*Bpl**2
                               betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 
     & 1./4.*Bpr**2
                               betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 
     & 1./4.*Bml**2
                               betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 
     & 1./4.*Bmr**2
                               ep     = 1e-15
                               ap1    = (1./2.)/(ep + betapl)**(4)
                               ap2    = (1./2.)/(ep + betapr)**(4)
                               am1    = (1./2.)/(ep + betaml)**(4)
                               am2    = (1./2.)/(ep + betamr)**(4)
                               !calculate weights
                               wp1 = ap1/(ap1+ap2)
                               wp2 = ap2/(ap1+ap2)
                               wm1 = am1/(am1+am2)
                               wm2 = am2/(am1+am2)
                               !mapping of the weights
                               bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                               bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                               bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                               bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                               wplw = bp1/(bp1+bp2)
                               wprw = bp2/(bp1+bp2)
                               wmlw = bm1/(bm1+bm2)
                               wmrw = bm2/(bm1+bm2)
                               !get face value of u
                               up = au
                               um = au
                               !picking side for the weights
                               if (up.ge.0) then
                                  wpl = max(wplw,wprw)
                                  wpr = min(wplw,wprw)
                               else
                                  wpl = min(wplw,wprw)
                                  wpr = max(wplw,wprw)
                               end if
                               if (um.ge.0) then
                                  wml = max(wmlw,wmrw)
                                  wmr = min(wmlw,wmrw)
                               else
                                  wml = min(wmlw,wmrw)
                                  wmr = max(wmlw,wmrw)
                               end if
                               expAd = 0.
                               if (orderOneFix.eq.1) then
                                  ! if (augmentPlot.eq.1) then
                                  !    agu(uc,13)=0
                                  ! end if
                                  if (aum*aup .le. 0) then
                                     aumax = max(abs(au),abs(aup),abs(
     & aum))
                                     expAd = (1./12.)*aumax*( u(i1-2*
     & s1,i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(
     & i1,i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*
     & s2,i3+2*s3,var))
                                     wpl = 0.5
                                     wpr = 0.5
                                     wml = 0.5
                                     wmr = 0.5
                                     ! plot au at where fix is turned on
                                     ! if (augmentPlot.eq.1) then
                                     !    agu(uc,13)=au
                                     ! end if
                                  end if
                               end if
                               ! wpl = 0.5
                               ! wpr = 0.5
                               ! wml = 0.5
                               ! wmr = 0.5
                               Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,var)
     &  + 5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
                               Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 5.*
     & u(i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
                               Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*
     & s3,var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
                               Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,
     & var) + 5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
                               Fp  = wpl*Fpl + wpr*Fpr
                               Fm  = wml*Fml + wmr*Fmr
                                  agu(vard,var)= au*(Fp - Fm)/drl + 
     & expAd/drl
                             !     if ( i1 .eq. 14 .and. i2.eq.16 ) then
                                     !     if ( i1 .eq. 14) then
                                     ! write(*,'(" betapl,betapr,betaml,betamr = ",4F10.8)')  betapl,betapr,betaml,betamr
                                     ! write(*,'(" ap1,ap2,am1,am2 = ",4e10.8)')  ap1,ap2,am1,am2
                                     ! write(*,'(" wp1,wpr,wm1,wm2 = ",4F10.8)')  wp1,wp2,wm1,wm2
                                     ! write(*,'(" bp1,bp2,bm1,bm2 = ",4F10.8)')  bp1,bp2,bm1,bm2
                                     ! write(*,'(" wplw,wprw,wmlw,wmrw = ",4F10.8)')  wplw,wprw,wmlw,wmrw
                                     ! write(*,'(" wpl,wpr,wml,wmr = ",4F10.8)')  wpl,wpr,wml,wmr
                                     ! write(*,'(" agu = ",F10.8)') agu(vard,var)
                             !     end if
                          else ! if (upwindOrder.eq.4)
                             write(*,'(" getBewnoAdvection: only 4th 
     & order is avaliable now")' )
                             stop 777
                          end if
                      ! #If ("BWENO" == "CENTERED") 
                      !  ! -- centered advection ---
                      !  ! write(*,'(" getAdvection -- centered")')
                      !  #If "2" == "2"
                      !   ! -- 2D --
                      !   agu(uc,uc)=UU(uc)*UX(uc)
                      !   agu(vc,uc)=UU(vc)*UY(uc)
                      !   agu(uc,vc)=UU(uc)*UX(vc)
                      !   agu(vc,vc)=UU(vc)*UY(vc)
                      !  #Else
                      !   ! -- 3D -- *check me*
                      !   agu(uc,uc)=UU(uc)*UX(uc)
                      !   agu(vc,uc)=UU(vc)*UY(uc)
                      !   agu(wc,uc)=UU(wc)*UY(uc)
                      !   agu(uc,vc)=UU(uc)*UX(vc)
                      !   agu(vc,vc)=UU(vc)*UY(vc)
                      !   agu(wc,vc)=UU(wc)*UY(vc)
                      !   agu(uc,wc)=UU(uc)*UX(wc)
                      !   agu(vc,wc)=UU(vc)*UY(wc)
                      !   agu(wc,wc)=UU(wc)*UY(wc)
                      !  #End
                      ! #Elif "BWENO" == "BWENO" 
                      !   ! --- upwind scheme ---
                      !   ! for testing output this next message:
                      !   if( t.le. 0. )then
                      !     write(*,'(" getAdvection upwind scheme (7)")') 
                      !   end if
                      !   getUpwindAdvection(u,i1,i2,i3,NONE,2,4,curvilinear, agu)
                      ! #Elif "BWENO" == "BWENO" 
                      !   ! --- Bweno scheme ---
                      !   ! for testing output this message:
                      !   if( t.le. 0. )then
                      !      write(*,'(" getAdvection BWENO scheme (7)")') 
                      !   end if
                      !   getBwenoAdvection(u,i1,i2,i3,NONE,2,4,curvilinear, agu)
                      ! #Else
                      !   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
                      !   stop 999
                      ! #End
                        ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-ux42(
     & i1,i2,i3,pc)
                        ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-uy42(
     & i1,i2,i3,pc)
                    ! include implicit terms - diffusion
                     uti(i1,i2,i3,uc)= nu*ulaplacian42(i1,i2,i3,uc)
                     uti(i1,i2,i3,vc)= nu*ulaplacian42(i1,i2,i3,vc)
                 end if
                end do
                end do
                end do
               else
                do i3=n3a,n3b
                do i2=n2a,n2b
                do i1=n1a,n1b
                   ! explicit terms only, no diffusion
                     ! upwind approximation/bweno
                       ! --- Bweno scheme ---
                       ! for testing output this message:
                       ! if( t.le. 0. )then
                       !    write(*,'(" getAdvection BWENO scheme (7)")') 
                       ! end if
                          gvU = u(i1,i2,i3,uc) - uu(i1,i2,i3,uc)
                          gvV = u(i1,i2,i3,vc) - uu(i1,i2,i3,vc)
                          gvW = u(i1,i2,i3,wc) - uu(i1,i2,i3,wc)
                         !bweno curvulinear
                         orderOneFix = 1
                         augmentPlot = 0
                         if( upwindOrder.eq.4 )then
                            ! u*ux-------------------------------------------------------------------------------
                            ! au = rsxy(i1,i2,i3,0,0)*u(i1,i2,i3,uc)+rsxy(i1,i2,i3,0,1)*u(i1,i2,i3,vc) + rsxy(i1,i2,i3,0,2)*u(i1,i2,i3,wc)
                            ! aum   = rsxy(i1-1,i2,i3,0,0)*u(i1-1,i2,i3,uc)+rsxy(i1-1,i2,i3,0,1)*u(i1-1,i2,i3,vc)+ rsxy(i1-1,i2,i3,0,2)*u(i1-1,i2,i3,wc)
                            ! aup   = rsxy(i1+1,i2,i3,0,0)*u(i1+1,i2,i3,uc)+rsxy(i1+1,i2,i3,0,1)*u(i1+1,i2,i3,vc)+ rsxy(i1+1,i2,i3,0,2)*u(i1+1,i2,i3,wc)
                            au    = rsxy(i1  ,i2,i3,0,0)*uu(i1  ,i2,i3,
     & uc) + rsxy(i1  ,i2,i3,0,1) * uu(i1  ,i2,i3,vc)
                            aum   = rsxy(i1-1,i2,i3,0,0)*(u(i1-1,i2,i3,
     & uc) - gvU ) + rsxy(i1-1,i2,i3,0,1) * (u(i1-1,i2,i3,vc)- gvV )
                            aup   = rsxy(i1+1,i2,i3,0,0)*(u(i1+1,i2,i3,
     & uc) - gvU ) + rsxy(i1+1,i2,i3,0,1) * (u(i1+1,i2,i3,vc)- gvV )
                            s1=1
                            s2=0
                            s3=0
                            var = uc
                            vard= uc
                            drl = dr(0)
                              Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*u(
     & i1   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
                              Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    u(
     & i1-s1,i2-s2,i3-s3,var)
                              Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) - 
     & 2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
                              Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) - 
     &    u(i1   ,i2   ,i3   ,var)
                              Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,
     & i2-s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
                              Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var)
                              Amr   = Apl
                              Bmr   = Bpl
                              betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 
     & 1./4.*Bpl**2
                              betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 
     & 1./4.*Bpr**2
                              betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 
     & 1./4.*Bml**2
                              betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 
     & 1./4.*Bmr**2
                              ep     = 1e-15
                              ap1    = (1./2.)/(ep + betapl)**(4)
                              ap2    = (1./2.)/(ep + betapr)**(4)
                              am1    = (1./2.)/(ep + betaml)**(4)
                              am2    = (1./2.)/(ep + betamr)**(4)
                              !calculate weights
                              wp1 = ap1/(ap1+ap2)
                              wp2 = ap2/(ap1+ap2)
                              wm1 = am1/(am1+am2)
                              wm2 = am2/(am1+am2)
                              !mapping of the weights
                              bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                              bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                              bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                              bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                              wplw = bp1/(bp1+bp2)
                              wprw = bp2/(bp1+bp2)
                              wmlw = bm1/(bm1+bm2)
                              wmrw = bm2/(bm1+bm2)
                              !get face value of u
                              up = au
                              um = au
                              !picking side for the weights
                              if (up.ge.0) then
                                 wpl = max(wplw,wprw)
                                 wpr = min(wplw,wprw)
                              else
                                 wpl = min(wplw,wprw)
                                 wpr = max(wplw,wprw)
                              end if
                              if (um.ge.0) then
                                 wml = max(wmlw,wmrw)
                                 wmr = min(wmlw,wmrw)
                              else
                                 wml = min(wmlw,wmrw)
                                 wmr = max(wmlw,wmrw)
                              end if
                              expAd = 0.
                              if (orderOneFix.eq.1) then
                                 ! if (augmentPlot.eq.1) then
                                 !    agu(uc,13)=0
                                 ! end if
                                 if (aum*aup .le. 0) then
                                    aumax = max(abs(au),abs(aup),abs(
     & aum))
                                    expAd = (1./12.)*aumax*( u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(i1,
     & i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*s2,
     & i3+2*s3,var))
                                    wpl = 0.5
                                    wpr = 0.5
                                    wml = 0.5
                                    wmr = 0.5
                                    ! plot au at where fix is turned on
                                    ! if (augmentPlot.eq.1) then
                                    !    agu(uc,13)=au
                                    ! end if
                                 end if
                              end if
                              ! wpl = 0.5
                              ! wpr = 0.5
                              ! wml = 0.5
                              ! wmr = 0.5
                              Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,var) 
     & + 5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
                              Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 5.*
     & u(i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
                              Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*
     & s3,var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
                              Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,var)
     &  + 5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
                              Fp  = wpl*Fpl + wpr*Fpr
                              Fm  = wml*Fml + wmr*Fmr
                                 agu(vard,var)= au*(Fp - Fm)/drl + 
     & expAd/drl
                            !     if ( i1 .eq. 14 .and. i2.eq.16 ) then
                                    !     if ( i1 .eq. 14) then
                                    ! write(*,'(" betapl,betapr,betaml,betamr = ",4F10.8)')  betapl,betapr,betaml,betamr
                                    ! write(*,'(" ap1,ap2,am1,am2 = ",4e10.8)')  ap1,ap2,am1,am2
                                    ! write(*,'(" wp1,wpr,wm1,wm2 = ",4F10.8)')  wp1,wp2,wm1,wm2
                                    ! write(*,'(" bp1,bp2,bm1,bm2 = ",4F10.8)')  bp1,bp2,bm1,bm2
                                    ! write(*,'(" wplw,wprw,wmlw,wmrw = ",4F10.8)')  wplw,wprw,wmlw,wmrw
                                    ! write(*,'(" wpl,wpr,wml,wmr = ",4F10.8)')  wpl,wpr,wml,wmr
                                    ! write(*,'(" agu = ",F10.8)') agu(vard,var)
                            !     end if
                            !u*vx
                            s1=1
                            s2=0
                            s3=0
                            var = vc
                            vard= uc
                            drl = dr(0)
                              Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*u(
     & i1   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
                              Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    u(
     & i1-s1,i2-s2,i3-s3,var)
                              Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) - 
     & 2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
                              Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) - 
     &    u(i1   ,i2   ,i3   ,var)
                              Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,
     & i2-s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
                              Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var)
                              Amr   = Apl
                              Bmr   = Bpl
                              betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 
     & 1./4.*Bpl**2
                              betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 
     & 1./4.*Bpr**2
                              betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 
     & 1./4.*Bml**2
                              betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 
     & 1./4.*Bmr**2
                              ep     = 1e-15
                              ap1    = (1./2.)/(ep + betapl)**(4)
                              ap2    = (1./2.)/(ep + betapr)**(4)
                              am1    = (1./2.)/(ep + betaml)**(4)
                              am2    = (1./2.)/(ep + betamr)**(4)
                              !calculate weights
                              wp1 = ap1/(ap1+ap2)
                              wp2 = ap2/(ap1+ap2)
                              wm1 = am1/(am1+am2)
                              wm2 = am2/(am1+am2)
                              !mapping of the weights
                              bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                              bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                              bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                              bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                              wplw = bp1/(bp1+bp2)
                              wprw = bp2/(bp1+bp2)
                              wmlw = bm1/(bm1+bm2)
                              wmrw = bm2/(bm1+bm2)
                              !get face value of u
                              up = au
                              um = au
                              !picking side for the weights
                              if (up.ge.0) then
                                 wpl = max(wplw,wprw)
                                 wpr = min(wplw,wprw)
                              else
                                 wpl = min(wplw,wprw)
                                 wpr = max(wplw,wprw)
                              end if
                              if (um.ge.0) then
                                 wml = max(wmlw,wmrw)
                                 wmr = min(wmlw,wmrw)
                              else
                                 wml = min(wmlw,wmrw)
                                 wmr = max(wmlw,wmrw)
                              end if
                              expAd = 0.
                              if (orderOneFix.eq.1) then
                                 ! if (augmentPlot.eq.1) then
                                 !    agu(uc,13)=0
                                 ! end if
                                 if (aum*aup .le. 0) then
                                    aumax = max(abs(au),abs(aup),abs(
     & aum))
                                    expAd = (1./12.)*aumax*( u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(i1,
     & i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*s2,
     & i3+2*s3,var))
                                    wpl = 0.5
                                    wpr = 0.5
                                    wml = 0.5
                                    wmr = 0.5
                                    ! plot au at where fix is turned on
                                    ! if (augmentPlot.eq.1) then
                                    !    agu(uc,13)=au
                                    ! end if
                                 end if
                              end if
                              ! wpl = 0.5
                              ! wpr = 0.5
                              ! wml = 0.5
                              ! wmr = 0.5
                              Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,var) 
     & + 5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
                              Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 5.*
     & u(i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
                              Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*
     & s3,var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
                              Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,var)
     &  + 5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
                              Fp  = wpl*Fpl + wpr*Fpr
                              Fm  = wml*Fml + wmr*Fmr
                                 agu(vard,var)= au*(Fp - Fm)/drl + 
     & expAd/drl
                            !     if ( i1 .eq. 14 .and. i2.eq.16 ) then
                                    !     if ( i1 .eq. 14) then
                                    ! write(*,'(" betapl,betapr,betaml,betamr = ",4F10.8)')  betapl,betapr,betaml,betamr
                                    ! write(*,'(" ap1,ap2,am1,am2 = ",4e10.8)')  ap1,ap2,am1,am2
                                    ! write(*,'(" wp1,wpr,wm1,wm2 = ",4F10.8)')  wp1,wp2,wm1,wm2
                                    ! write(*,'(" bp1,bp2,bm1,bm2 = ",4F10.8)')  bp1,bp2,bm1,bm2
                                    ! write(*,'(" wplw,wprw,wmlw,wmrw = ",4F10.8)')  wplw,wprw,wmlw,wmrw
                                    ! write(*,'(" wpl,wpr,wml,wmr = ",4F10.8)')  wpl,wpr,wml,wmr
                                    ! write(*,'(" agu = ",F10.8)') agu(vard,var)
                            !     end if
                            !v*uy
                            au    = rsxy(i1,i2  ,i3,1,0)*uu(i1,i2  ,i3,
     & uc)+rsxy(i1,i2,i3,1,1)*uu(i1,i2,i3,vc)
                            aum   = rsxy(i1,i2-1,i3,1,0)*(u(i1,i2-1,i3,
     & uc) - gvU) + rsxy(i1,i2-1,i3,1,1)*(u(i1,i2-1,i3,vc) - gvV)
                            aup   = rsxy(i1,i2+1,i3,1,0)*(u(i1,i2+1,i3,
     & uc) - gvU) + rsxy(i1,i2+1,i3,1,1)*(u(i1,i2+1,i3,vc) - gvV)
                            s1=0
                            s2=1
                            s3=0
                            var = uc
                            vard= vc
                            drl = dr(1)
                              Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*u(
     & i1   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
                              Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    u(
     & i1-s1,i2-s2,i3-s3,var)
                              Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) - 
     & 2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
                              Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) - 
     &    u(i1   ,i2   ,i3   ,var)
                              Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,
     & i2-s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
                              Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var)
                              Amr   = Apl
                              Bmr   = Bpl
                              betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 
     & 1./4.*Bpl**2
                              betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 
     & 1./4.*Bpr**2
                              betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 
     & 1./4.*Bml**2
                              betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 
     & 1./4.*Bmr**2
                              ep     = 1e-15
                              ap1    = (1./2.)/(ep + betapl)**(4)
                              ap2    = (1./2.)/(ep + betapr)**(4)
                              am1    = (1./2.)/(ep + betaml)**(4)
                              am2    = (1./2.)/(ep + betamr)**(4)
                              !calculate weights
                              wp1 = ap1/(ap1+ap2)
                              wp2 = ap2/(ap1+ap2)
                              wm1 = am1/(am1+am2)
                              wm2 = am2/(am1+am2)
                              !mapping of the weights
                              bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                              bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                              bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                              bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                              wplw = bp1/(bp1+bp2)
                              wprw = bp2/(bp1+bp2)
                              wmlw = bm1/(bm1+bm2)
                              wmrw = bm2/(bm1+bm2)
                              !get face value of u
                              up = au
                              um = au
                              !picking side for the weights
                              if (up.ge.0) then
                                 wpl = max(wplw,wprw)
                                 wpr = min(wplw,wprw)
                              else
                                 wpl = min(wplw,wprw)
                                 wpr = max(wplw,wprw)
                              end if
                              if (um.ge.0) then
                                 wml = max(wmlw,wmrw)
                                 wmr = min(wmlw,wmrw)
                              else
                                 wml = min(wmlw,wmrw)
                                 wmr = max(wmlw,wmrw)
                              end if
                              expAd = 0.
                              if (orderOneFix.eq.1) then
                                 ! if (augmentPlot.eq.1) then
                                 !    agu(uc,13)=0
                                 ! end if
                                 if (aum*aup .le. 0) then
                                    aumax = max(abs(au),abs(aup),abs(
     & aum))
                                    expAd = (1./12.)*aumax*( u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(i1,
     & i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*s2,
     & i3+2*s3,var))
                                    wpl = 0.5
                                    wpr = 0.5
                                    wml = 0.5
                                    wmr = 0.5
                                    ! plot au at where fix is turned on
                                    ! if (augmentPlot.eq.1) then
                                    !    agu(uc,13)=au
                                    ! end if
                                 end if
                              end if
                              ! wpl = 0.5
                              ! wpr = 0.5
                              ! wml = 0.5
                              ! wmr = 0.5
                              Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,var) 
     & + 5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
                              Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 5.*
     & u(i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
                              Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*
     & s3,var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
                              Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,var)
     &  + 5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
                              Fp  = wpl*Fpl + wpr*Fpr
                              Fm  = wml*Fml + wmr*Fmr
                                 agu(vard,var)= au*(Fp - Fm)/drl + 
     & expAd/drl
                            !     if ( i1 .eq. 14 .and. i2.eq.16 ) then
                                    !     if ( i1 .eq. 14) then
                                    ! write(*,'(" betapl,betapr,betaml,betamr = ",4F10.8)')  betapl,betapr,betaml,betamr
                                    ! write(*,'(" ap1,ap2,am1,am2 = ",4e10.8)')  ap1,ap2,am1,am2
                                    ! write(*,'(" wp1,wpr,wm1,wm2 = ",4F10.8)')  wp1,wp2,wm1,wm2
                                    ! write(*,'(" bp1,bp2,bm1,bm2 = ",4F10.8)')  bp1,bp2,bm1,bm2
                                    ! write(*,'(" wplw,wprw,wmlw,wmrw = ",4F10.8)')  wplw,wprw,wmlw,wmrw
                                    ! write(*,'(" wpl,wpr,wml,wmr = ",4F10.8)')  wpl,wpr,wml,wmr
                                    ! write(*,'(" agu = ",F10.8)') agu(vard,var)
                            !     end if
                            !v*vy
                            s1=0
                            s2=1
                            s3=0
                            var = vc
                            vard= vc
                            drl = dr(1)
                              Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*u(
     & i1   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
                              Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    u(
     & i1-s1,i2-s2,i3-s3,var)
                              Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) - 
     & 2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
                              Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) - 
     &    u(i1   ,i2   ,i3   ,var)
                              Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,
     & i2-s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
                              Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var)
                              Amr   = Apl
                              Bmr   = Bpl
                              betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 
     & 1./4.*Bpl**2
                              betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 
     & 1./4.*Bpr**2
                              betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 
     & 1./4.*Bml**2
                              betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 
     & 1./4.*Bmr**2
                              ep     = 1e-15
                              ap1    = (1./2.)/(ep + betapl)**(4)
                              ap2    = (1./2.)/(ep + betapr)**(4)
                              am1    = (1./2.)/(ep + betaml)**(4)
                              am2    = (1./2.)/(ep + betamr)**(4)
                              !calculate weights
                              wp1 = ap1/(ap1+ap2)
                              wp2 = ap2/(ap1+ap2)
                              wm1 = am1/(am1+am2)
                              wm2 = am2/(am1+am2)
                              !mapping of the weights
                              bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                              bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                              bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                              bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                              wplw = bp1/(bp1+bp2)
                              wprw = bp2/(bp1+bp2)
                              wmlw = bm1/(bm1+bm2)
                              wmrw = bm2/(bm1+bm2)
                              !get face value of u
                              up = au
                              um = au
                              !picking side for the weights
                              if (up.ge.0) then
                                 wpl = max(wplw,wprw)
                                 wpr = min(wplw,wprw)
                              else
                                 wpl = min(wplw,wprw)
                                 wpr = max(wplw,wprw)
                              end if
                              if (um.ge.0) then
                                 wml = max(wmlw,wmrw)
                                 wmr = min(wmlw,wmrw)
                              else
                                 wml = min(wmlw,wmrw)
                                 wmr = max(wmlw,wmrw)
                              end if
                              expAd = 0.
                              if (orderOneFix.eq.1) then
                                 ! if (augmentPlot.eq.1) then
                                 !    agu(uc,13)=0
                                 ! end if
                                 if (aum*aup .le. 0) then
                                    aumax = max(abs(au),abs(aup),abs(
     & aum))
                                    expAd = (1./12.)*aumax*( u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(i1,
     & i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*s2,
     & i3+2*s3,var))
                                    wpl = 0.5
                                    wpr = 0.5
                                    wml = 0.5
                                    wmr = 0.5
                                    ! plot au at where fix is turned on
                                    ! if (augmentPlot.eq.1) then
                                    !    agu(uc,13)=au
                                    ! end if
                                 end if
                              end if
                              ! wpl = 0.5
                              ! wpr = 0.5
                              ! wml = 0.5
                              ! wmr = 0.5
                              Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,var) 
     & + 5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
                              Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 5.*
     & u(i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
                              Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*
     & s3,var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
                              Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,var)
     &  + 5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
                              Fp  = wpl*Fpl + wpr*Fpr
                              Fm  = wml*Fml + wmr*Fmr
                                 agu(vard,var)= au*(Fp - Fm)/drl + 
     & expAd/drl
                            !     if ( i1 .eq. 14 .and. i2.eq.16 ) then
                                    !     if ( i1 .eq. 14) then
                                    ! write(*,'(" betapl,betapr,betaml,betamr = ",4F10.8)')  betapl,betapr,betaml,betamr
                                    ! write(*,'(" ap1,ap2,am1,am2 = ",4e10.8)')  ap1,ap2,am1,am2
                                    ! write(*,'(" wp1,wpr,wm1,wm2 = ",4F10.8)')  wp1,wp2,wm1,wm2
                                    ! write(*,'(" bp1,bp2,bm1,bm2 = ",4F10.8)')  bp1,bp2,bm1,bm2
                                    ! write(*,'(" wplw,wprw,wmlw,wmrw = ",4F10.8)')  wplw,wprw,wmlw,wmrw
                                    ! write(*,'(" wpl,wpr,wml,wmr = ",4F10.8)')  wpl,wpr,wml,wmr
                                    ! write(*,'(" agu = ",F10.8)') agu(vard,var)
                            !     end if
                         else ! if (upwindOrder.eq.4)
                            write(*,'(" getBewnoAdvection: only 4th 
     & order is avaliable now")' )
                            stop 777
                         end if
                     ! #If ("BWENO" == "CENTERED") 
                     !  ! -- centered advection ---
                     !  ! write(*,'(" getAdvection -- centered")')
                     !  #If "2" == "2"
                     !   ! -- 2D --
                     !   agu(uc,uc)=UU(uc)*UX(uc)
                     !   agu(vc,uc)=UU(vc)*UY(uc)
                     !   agu(uc,vc)=UU(uc)*UX(vc)
                     !   agu(vc,vc)=UU(vc)*UY(vc)
                     !  #Else
                     !   ! -- 3D -- *check me*
                     !   agu(uc,uc)=UU(uc)*UX(uc)
                     !   agu(vc,uc)=UU(vc)*UY(uc)
                     !   agu(wc,uc)=UU(wc)*UY(uc)
                     !   agu(uc,vc)=UU(uc)*UX(vc)
                     !   agu(vc,vc)=UU(vc)*UY(vc)
                     !   agu(wc,vc)=UU(wc)*UY(vc)
                     !   agu(uc,wc)=UU(uc)*UX(wc)
                     !   agu(vc,wc)=UU(vc)*UY(wc)
                     !   agu(wc,wc)=UU(wc)*UY(wc)
                     !  #End
                     ! #Elif "BWENO" == "BWENO" 
                     !   ! --- upwind scheme ---
                     !   ! for testing output this next message:
                     !   if( t.le. 0. )then
                     !     write(*,'(" getAdvection upwind scheme (7)")') 
                     !   end if
                     !   getUpwindAdvection(u,i1,i2,i3,NONE,2,4,curvilinear, agu)
                     ! #Elif "BWENO" == "BWENO" 
                     !   ! --- Bweno scheme ---
                     !   ! for testing output this message:
                     !   if( t.le. 0. )then
                     !      write(*,'(" getAdvection BWENO scheme (7)")') 
                     !   end if
                     !   getBwenoAdvection(u,i1,i2,i3,NONE,2,4,curvilinear, agu)
                     ! #Else
                     !   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
                     !   stop 999
                     ! #End
                       ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-ux42(
     & i1,i2,i3,pc)
                       ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-uy42(
     & i1,i2,i3,pc)
                   ! include implicit terms - diffusion
                    uti(i1,i2,i3,uc)= nu*ulaplacian42(i1,i2,i3,uc)
                    uti(i1,i2,i3,vc)= nu*ulaplacian42(i1,i2,i3,vc)
                end do
                end do
                end do
               end if
             else if( implicitOption.eq.doNotComputeImplicitTerms )then
               if( useWhereMask.ne.0 )then
                do i3=n3a,n3b
                do i2=n2a,n2b
                do i1=n1a,n1b
                 if( mask(i1,i2,i3).gt.0 )then
                    ! explicit terms only, no diffusion
                      ! upwind approximation/bweno
                        ! --- Bweno scheme ---
                        ! for testing output this message:
                        ! if( t.le. 0. )then
                        !    write(*,'(" getAdvection BWENO scheme (7)")') 
                        ! end if
                           gvU = u(i1,i2,i3,uc) - uu(i1,i2,i3,uc)
                           gvV = u(i1,i2,i3,vc) - uu(i1,i2,i3,vc)
                           gvW = u(i1,i2,i3,wc) - uu(i1,i2,i3,wc)
                          !bweno curvulinear
                          orderOneFix = 1
                          augmentPlot = 0
                          if( upwindOrder.eq.4 )then
                             ! u*ux-------------------------------------------------------------------------------
                             ! au = rsxy(i1,i2,i3,0,0)*u(i1,i2,i3,uc)+rsxy(i1,i2,i3,0,1)*u(i1,i2,i3,vc) + rsxy(i1,i2,i3,0,2)*u(i1,i2,i3,wc)
                             ! aum   = rsxy(i1-1,i2,i3,0,0)*u(i1-1,i2,i3,uc)+rsxy(i1-1,i2,i3,0,1)*u(i1-1,i2,i3,vc)+ rsxy(i1-1,i2,i3,0,2)*u(i1-1,i2,i3,wc)
                             ! aup   = rsxy(i1+1,i2,i3,0,0)*u(i1+1,i2,i3,uc)+rsxy(i1+1,i2,i3,0,1)*u(i1+1,i2,i3,vc)+ rsxy(i1+1,i2,i3,0,2)*u(i1+1,i2,i3,wc)
                             au    = rsxy(i1  ,i2,i3,0,0)*uu(i1  ,i2,
     & i3,uc) + rsxy(i1  ,i2,i3,0,1) * uu(i1  ,i2,i3,vc)
                             aum   = rsxy(i1-1,i2,i3,0,0)*(u(i1-1,i2,
     & i3,uc) - gvU ) + rsxy(i1-1,i2,i3,0,1) * (u(i1-1,i2,i3,vc)- gvV 
     & )
                             aup   = rsxy(i1+1,i2,i3,0,0)*(u(i1+1,i2,
     & i3,uc) - gvU ) + rsxy(i1+1,i2,i3,0,1) * (u(i1+1,i2,i3,vc)- gvV 
     & )
                             s1=1
                             s2=0
                             s3=0
                             var = uc
                             vard= uc
                             drl = dr(0)
                               Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*u(
     & i1   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
                               Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    u(
     & i1-s1,i2-s2,i3-s3,var)
                               Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) -
     &  2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
                               Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) -
     &     u(i1   ,i2   ,i3   ,var)
                               Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,
     & i2-s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
                               Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var)
                               Amr   = Apl
                               Bmr   = Bpl
                               betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 
     & 1./4.*Bpl**2
                               betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 
     & 1./4.*Bpr**2
                               betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 
     & 1./4.*Bml**2
                               betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 
     & 1./4.*Bmr**2
                               ep     = 1e-15
                               ap1    = (1./2.)/(ep + betapl)**(4)
                               ap2    = (1./2.)/(ep + betapr)**(4)
                               am1    = (1./2.)/(ep + betaml)**(4)
                               am2    = (1./2.)/(ep + betamr)**(4)
                               !calculate weights
                               wp1 = ap1/(ap1+ap2)
                               wp2 = ap2/(ap1+ap2)
                               wm1 = am1/(am1+am2)
                               wm2 = am2/(am1+am2)
                               !mapping of the weights
                               bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                               bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                               bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                               bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                               wplw = bp1/(bp1+bp2)
                               wprw = bp2/(bp1+bp2)
                               wmlw = bm1/(bm1+bm2)
                               wmrw = bm2/(bm1+bm2)
                               !get face value of u
                               up = au
                               um = au
                               !picking side for the weights
                               if (up.ge.0) then
                                  wpl = max(wplw,wprw)
                                  wpr = min(wplw,wprw)
                               else
                                  wpl = min(wplw,wprw)
                                  wpr = max(wplw,wprw)
                               end if
                               if (um.ge.0) then
                                  wml = max(wmlw,wmrw)
                                  wmr = min(wmlw,wmrw)
                               else
                                  wml = min(wmlw,wmrw)
                                  wmr = max(wmlw,wmrw)
                               end if
                               expAd = 0.
                               if (orderOneFix.eq.1) then
                                  ! if (augmentPlot.eq.1) then
                                  !    agu(uc,13)=0
                                  ! end if
                                  if (aum*aup .le. 0) then
                                     aumax = max(abs(au),abs(aup),abs(
     & aum))
                                     expAd = (1./12.)*aumax*( u(i1-2*
     & s1,i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(
     & i1,i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*
     & s2,i3+2*s3,var))
                                     wpl = 0.5
                                     wpr = 0.5
                                     wml = 0.5
                                     wmr = 0.5
                                     ! plot au at where fix is turned on
                                     ! if (augmentPlot.eq.1) then
                                     !    agu(uc,13)=au
                                     ! end if
                                  end if
                               end if
                               ! wpl = 0.5
                               ! wpr = 0.5
                               ! wml = 0.5
                               ! wmr = 0.5
                               Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,var)
     &  + 5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
                               Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 5.*
     & u(i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
                               Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*
     & s3,var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
                               Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,
     & var) + 5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
                               Fp  = wpl*Fpl + wpr*Fpr
                               Fm  = wml*Fml + wmr*Fmr
                                  agu(vard,var)= au*(Fp - Fm)/drl + 
     & expAd/drl
                             !     if ( i1 .eq. 14 .and. i2.eq.16 ) then
                                     !     if ( i1 .eq. 14) then
                                     ! write(*,'(" betapl,betapr,betaml,betamr = ",4F10.8)')  betapl,betapr,betaml,betamr
                                     ! write(*,'(" ap1,ap2,am1,am2 = ",4e10.8)')  ap1,ap2,am1,am2
                                     ! write(*,'(" wp1,wpr,wm1,wm2 = ",4F10.8)')  wp1,wp2,wm1,wm2
                                     ! write(*,'(" bp1,bp2,bm1,bm2 = ",4F10.8)')  bp1,bp2,bm1,bm2
                                     ! write(*,'(" wplw,wprw,wmlw,wmrw = ",4F10.8)')  wplw,wprw,wmlw,wmrw
                                     ! write(*,'(" wpl,wpr,wml,wmr = ",4F10.8)')  wpl,wpr,wml,wmr
                                     ! write(*,'(" agu = ",F10.8)') agu(vard,var)
                             !     end if
                             !u*vx
                             s1=1
                             s2=0
                             s3=0
                             var = vc
                             vard= uc
                             drl = dr(0)
                               Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*u(
     & i1   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
                               Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    u(
     & i1-s1,i2-s2,i3-s3,var)
                               Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) -
     &  2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
                               Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) -
     &     u(i1   ,i2   ,i3   ,var)
                               Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,
     & i2-s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
                               Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var)
                               Amr   = Apl
                               Bmr   = Bpl
                               betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 
     & 1./4.*Bpl**2
                               betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 
     & 1./4.*Bpr**2
                               betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 
     & 1./4.*Bml**2
                               betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 
     & 1./4.*Bmr**2
                               ep     = 1e-15
                               ap1    = (1./2.)/(ep + betapl)**(4)
                               ap2    = (1./2.)/(ep + betapr)**(4)
                               am1    = (1./2.)/(ep + betaml)**(4)
                               am2    = (1./2.)/(ep + betamr)**(4)
                               !calculate weights
                               wp1 = ap1/(ap1+ap2)
                               wp2 = ap2/(ap1+ap2)
                               wm1 = am1/(am1+am2)
                               wm2 = am2/(am1+am2)
                               !mapping of the weights
                               bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                               bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                               bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                               bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                               wplw = bp1/(bp1+bp2)
                               wprw = bp2/(bp1+bp2)
                               wmlw = bm1/(bm1+bm2)
                               wmrw = bm2/(bm1+bm2)
                               !get face value of u
                               up = au
                               um = au
                               !picking side for the weights
                               if (up.ge.0) then
                                  wpl = max(wplw,wprw)
                                  wpr = min(wplw,wprw)
                               else
                                  wpl = min(wplw,wprw)
                                  wpr = max(wplw,wprw)
                               end if
                               if (um.ge.0) then
                                  wml = max(wmlw,wmrw)
                                  wmr = min(wmlw,wmrw)
                               else
                                  wml = min(wmlw,wmrw)
                                  wmr = max(wmlw,wmrw)
                               end if
                               expAd = 0.
                               if (orderOneFix.eq.1) then
                                  ! if (augmentPlot.eq.1) then
                                  !    agu(uc,13)=0
                                  ! end if
                                  if (aum*aup .le. 0) then
                                     aumax = max(abs(au),abs(aup),abs(
     & aum))
                                     expAd = (1./12.)*aumax*( u(i1-2*
     & s1,i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(
     & i1,i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*
     & s2,i3+2*s3,var))
                                     wpl = 0.5
                                     wpr = 0.5
                                     wml = 0.5
                                     wmr = 0.5
                                     ! plot au at where fix is turned on
                                     ! if (augmentPlot.eq.1) then
                                     !    agu(uc,13)=au
                                     ! end if
                                  end if
                               end if
                               ! wpl = 0.5
                               ! wpr = 0.5
                               ! wml = 0.5
                               ! wmr = 0.5
                               Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,var)
     &  + 5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
                               Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 5.*
     & u(i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
                               Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*
     & s3,var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
                               Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,
     & var) + 5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
                               Fp  = wpl*Fpl + wpr*Fpr
                               Fm  = wml*Fml + wmr*Fmr
                                  agu(vard,var)= au*(Fp - Fm)/drl + 
     & expAd/drl
                             !     if ( i1 .eq. 14 .and. i2.eq.16 ) then
                                     !     if ( i1 .eq. 14) then
                                     ! write(*,'(" betapl,betapr,betaml,betamr = ",4F10.8)')  betapl,betapr,betaml,betamr
                                     ! write(*,'(" ap1,ap2,am1,am2 = ",4e10.8)')  ap1,ap2,am1,am2
                                     ! write(*,'(" wp1,wpr,wm1,wm2 = ",4F10.8)')  wp1,wp2,wm1,wm2
                                     ! write(*,'(" bp1,bp2,bm1,bm2 = ",4F10.8)')  bp1,bp2,bm1,bm2
                                     ! write(*,'(" wplw,wprw,wmlw,wmrw = ",4F10.8)')  wplw,wprw,wmlw,wmrw
                                     ! write(*,'(" wpl,wpr,wml,wmr = ",4F10.8)')  wpl,wpr,wml,wmr
                                     ! write(*,'(" agu = ",F10.8)') agu(vard,var)
                             !     end if
                             !v*uy
                             au    = rsxy(i1,i2  ,i3,1,0)*uu(i1,i2  ,
     & i3,uc)+rsxy(i1,i2,i3,1,1)*uu(i1,i2,i3,vc)
                             aum   = rsxy(i1,i2-1,i3,1,0)*(u(i1,i2-1,
     & i3,uc) - gvU) + rsxy(i1,i2-1,i3,1,1)*(u(i1,i2-1,i3,vc) - gvV)
                             aup   = rsxy(i1,i2+1,i3,1,0)*(u(i1,i2+1,
     & i3,uc) - gvU) + rsxy(i1,i2+1,i3,1,1)*(u(i1,i2+1,i3,vc) - gvV)
                             s1=0
                             s2=1
                             s3=0
                             var = uc
                             vard= vc
                             drl = dr(1)
                               Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*u(
     & i1   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
                               Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    u(
     & i1-s1,i2-s2,i3-s3,var)
                               Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) -
     &  2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
                               Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) -
     &     u(i1   ,i2   ,i3   ,var)
                               Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,
     & i2-s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
                               Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var)
                               Amr   = Apl
                               Bmr   = Bpl
                               betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 
     & 1./4.*Bpl**2
                               betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 
     & 1./4.*Bpr**2
                               betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 
     & 1./4.*Bml**2
                               betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 
     & 1./4.*Bmr**2
                               ep     = 1e-15
                               ap1    = (1./2.)/(ep + betapl)**(4)
                               ap2    = (1./2.)/(ep + betapr)**(4)
                               am1    = (1./2.)/(ep + betaml)**(4)
                               am2    = (1./2.)/(ep + betamr)**(4)
                               !calculate weights
                               wp1 = ap1/(ap1+ap2)
                               wp2 = ap2/(ap1+ap2)
                               wm1 = am1/(am1+am2)
                               wm2 = am2/(am1+am2)
                               !mapping of the weights
                               bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                               bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                               bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                               bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                               wplw = bp1/(bp1+bp2)
                               wprw = bp2/(bp1+bp2)
                               wmlw = bm1/(bm1+bm2)
                               wmrw = bm2/(bm1+bm2)
                               !get face value of u
                               up = au
                               um = au
                               !picking side for the weights
                               if (up.ge.0) then
                                  wpl = max(wplw,wprw)
                                  wpr = min(wplw,wprw)
                               else
                                  wpl = min(wplw,wprw)
                                  wpr = max(wplw,wprw)
                               end if
                               if (um.ge.0) then
                                  wml = max(wmlw,wmrw)
                                  wmr = min(wmlw,wmrw)
                               else
                                  wml = min(wmlw,wmrw)
                                  wmr = max(wmlw,wmrw)
                               end if
                               expAd = 0.
                               if (orderOneFix.eq.1) then
                                  ! if (augmentPlot.eq.1) then
                                  !    agu(uc,13)=0
                                  ! end if
                                  if (aum*aup .le. 0) then
                                     aumax = max(abs(au),abs(aup),abs(
     & aum))
                                     expAd = (1./12.)*aumax*( u(i1-2*
     & s1,i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(
     & i1,i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*
     & s2,i3+2*s3,var))
                                     wpl = 0.5
                                     wpr = 0.5
                                     wml = 0.5
                                     wmr = 0.5
                                     ! plot au at where fix is turned on
                                     ! if (augmentPlot.eq.1) then
                                     !    agu(uc,13)=au
                                     ! end if
                                  end if
                               end if
                               ! wpl = 0.5
                               ! wpr = 0.5
                               ! wml = 0.5
                               ! wmr = 0.5
                               Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,var)
     &  + 5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
                               Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 5.*
     & u(i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
                               Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*
     & s3,var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
                               Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,
     & var) + 5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
                               Fp  = wpl*Fpl + wpr*Fpr
                               Fm  = wml*Fml + wmr*Fmr
                                  agu(vard,var)= au*(Fp - Fm)/drl + 
     & expAd/drl
                             !     if ( i1 .eq. 14 .and. i2.eq.16 ) then
                                     !     if ( i1 .eq. 14) then
                                     ! write(*,'(" betapl,betapr,betaml,betamr = ",4F10.8)')  betapl,betapr,betaml,betamr
                                     ! write(*,'(" ap1,ap2,am1,am2 = ",4e10.8)')  ap1,ap2,am1,am2
                                     ! write(*,'(" wp1,wpr,wm1,wm2 = ",4F10.8)')  wp1,wp2,wm1,wm2
                                     ! write(*,'(" bp1,bp2,bm1,bm2 = ",4F10.8)')  bp1,bp2,bm1,bm2
                                     ! write(*,'(" wplw,wprw,wmlw,wmrw = ",4F10.8)')  wplw,wprw,wmlw,wmrw
                                     ! write(*,'(" wpl,wpr,wml,wmr = ",4F10.8)')  wpl,wpr,wml,wmr
                                     ! write(*,'(" agu = ",F10.8)') agu(vard,var)
                             !     end if
                             !v*vy
                             s1=0
                             s2=1
                             s3=0
                             var = vc
                             vard= vc
                             drl = dr(1)
                               Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*u(
     & i1   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
                               Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    u(
     & i1-s1,i2-s2,i3-s3,var)
                               Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) -
     &  2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
                               Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) -
     &     u(i1   ,i2   ,i3   ,var)
                               Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,
     & i2-s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
                               Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var)
                               Amr   = Apl
                               Bmr   = Bpl
                               betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 
     & 1./4.*Bpl**2
                               betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 
     & 1./4.*Bpr**2
                               betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 
     & 1./4.*Bml**2
                               betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 
     & 1./4.*Bmr**2
                               ep     = 1e-15
                               ap1    = (1./2.)/(ep + betapl)**(4)
                               ap2    = (1./2.)/(ep + betapr)**(4)
                               am1    = (1./2.)/(ep + betaml)**(4)
                               am2    = (1./2.)/(ep + betamr)**(4)
                               !calculate weights
                               wp1 = ap1/(ap1+ap2)
                               wp2 = ap2/(ap1+ap2)
                               wm1 = am1/(am1+am2)
                               wm2 = am2/(am1+am2)
                               !mapping of the weights
                               bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                               bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                               bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                               bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                               wplw = bp1/(bp1+bp2)
                               wprw = bp2/(bp1+bp2)
                               wmlw = bm1/(bm1+bm2)
                               wmrw = bm2/(bm1+bm2)
                               !get face value of u
                               up = au
                               um = au
                               !picking side for the weights
                               if (up.ge.0) then
                                  wpl = max(wplw,wprw)
                                  wpr = min(wplw,wprw)
                               else
                                  wpl = min(wplw,wprw)
                                  wpr = max(wplw,wprw)
                               end if
                               if (um.ge.0) then
                                  wml = max(wmlw,wmrw)
                                  wmr = min(wmlw,wmrw)
                               else
                                  wml = min(wmlw,wmrw)
                                  wmr = max(wmlw,wmrw)
                               end if
                               expAd = 0.
                               if (orderOneFix.eq.1) then
                                  ! if (augmentPlot.eq.1) then
                                  !    agu(uc,13)=0
                                  ! end if
                                  if (aum*aup .le. 0) then
                                     aumax = max(abs(au),abs(aup),abs(
     & aum))
                                     expAd = (1./12.)*aumax*( u(i1-2*
     & s1,i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(
     & i1,i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*
     & s2,i3+2*s3,var))
                                     wpl = 0.5
                                     wpr = 0.5
                                     wml = 0.5
                                     wmr = 0.5
                                     ! plot au at where fix is turned on
                                     ! if (augmentPlot.eq.1) then
                                     !    agu(uc,13)=au
                                     ! end if
                                  end if
                               end if
                               ! wpl = 0.5
                               ! wpr = 0.5
                               ! wml = 0.5
                               ! wmr = 0.5
                               Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,var)
     &  + 5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
                               Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 5.*
     & u(i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
                               Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*
     & s3,var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
                               Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,
     & var) + 5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
                               Fp  = wpl*Fpl + wpr*Fpr
                               Fm  = wml*Fml + wmr*Fmr
                                  agu(vard,var)= au*(Fp - Fm)/drl + 
     & expAd/drl
                             !     if ( i1 .eq. 14 .and. i2.eq.16 ) then
                                     !     if ( i1 .eq. 14) then
                                     ! write(*,'(" betapl,betapr,betaml,betamr = ",4F10.8)')  betapl,betapr,betaml,betamr
                                     ! write(*,'(" ap1,ap2,am1,am2 = ",4e10.8)')  ap1,ap2,am1,am2
                                     ! write(*,'(" wp1,wpr,wm1,wm2 = ",4F10.8)')  wp1,wp2,wm1,wm2
                                     ! write(*,'(" bp1,bp2,bm1,bm2 = ",4F10.8)')  bp1,bp2,bm1,bm2
                                     ! write(*,'(" wplw,wprw,wmlw,wmrw = ",4F10.8)')  wplw,wprw,wmlw,wmrw
                                     ! write(*,'(" wpl,wpr,wml,wmr = ",4F10.8)')  wpl,wpr,wml,wmr
                                     ! write(*,'(" agu = ",F10.8)') agu(vard,var)
                             !     end if
                          else ! if (upwindOrder.eq.4)
                             write(*,'(" getBewnoAdvection: only 4th 
     & order is avaliable now")' )
                             stop 777
                          end if
                      ! #If ("BWENO" == "CENTERED") 
                      !  ! -- centered advection ---
                      !  ! write(*,'(" getAdvection -- centered")')
                      !  #If "2" == "2"
                      !   ! -- 2D --
                      !   agu(uc,uc)=UU(uc)*UX(uc)
                      !   agu(vc,uc)=UU(vc)*UY(uc)
                      !   agu(uc,vc)=UU(uc)*UX(vc)
                      !   agu(vc,vc)=UU(vc)*UY(vc)
                      !  #Else
                      !   ! -- 3D -- *check me*
                      !   agu(uc,uc)=UU(uc)*UX(uc)
                      !   agu(vc,uc)=UU(vc)*UY(uc)
                      !   agu(wc,uc)=UU(wc)*UY(uc)
                      !   agu(uc,vc)=UU(uc)*UX(vc)
                      !   agu(vc,vc)=UU(vc)*UY(vc)
                      !   agu(wc,vc)=UU(wc)*UY(vc)
                      !   agu(uc,wc)=UU(uc)*UX(wc)
                      !   agu(vc,wc)=UU(vc)*UY(wc)
                      !   agu(wc,wc)=UU(wc)*UY(wc)
                      !  #End
                      ! #Elif "BWENO" == "BWENO" 
                      !   ! --- upwind scheme ---
                      !   ! for testing output this next message:
                      !   if( t.le. 0. )then
                      !     write(*,'(" getAdvection upwind scheme (7)")') 
                      !   end if
                      !   getUpwindAdvection(u,i1,i2,i3,NONE,2,4,curvilinear, agu)
                      ! #Elif "BWENO" == "BWENO" 
                      !   ! --- Bweno scheme ---
                      !   ! for testing output this message:
                      !   if( t.le. 0. )then
                      !      write(*,'(" getAdvection BWENO scheme (7)")') 
                      !   end if
                      !   getBwenoAdvection(u,i1,i2,i3,NONE,2,4,curvilinear, agu)
                      ! #Else
                      !   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
                      !   stop 999
                      ! #End
                        ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-ux42(
     & i1,i2,i3,pc)
                        ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-uy42(
     & i1,i2,i3,pc)
                 end if
                end do
                end do
                end do
               else
                do i3=n3a,n3b
                do i2=n2a,n2b
                do i1=n1a,n1b
                   ! explicit terms only, no diffusion
                     ! upwind approximation/bweno
                       ! --- Bweno scheme ---
                       ! for testing output this message:
                       ! if( t.le. 0. )then
                       !    write(*,'(" getAdvection BWENO scheme (7)")') 
                       ! end if
                          gvU = u(i1,i2,i3,uc) - uu(i1,i2,i3,uc)
                          gvV = u(i1,i2,i3,vc) - uu(i1,i2,i3,vc)
                          gvW = u(i1,i2,i3,wc) - uu(i1,i2,i3,wc)
                         !bweno curvulinear
                         orderOneFix = 1
                         augmentPlot = 0
                         if( upwindOrder.eq.4 )then
                            ! u*ux-------------------------------------------------------------------------------
                            ! au = rsxy(i1,i2,i3,0,0)*u(i1,i2,i3,uc)+rsxy(i1,i2,i3,0,1)*u(i1,i2,i3,vc) + rsxy(i1,i2,i3,0,2)*u(i1,i2,i3,wc)
                            ! aum   = rsxy(i1-1,i2,i3,0,0)*u(i1-1,i2,i3,uc)+rsxy(i1-1,i2,i3,0,1)*u(i1-1,i2,i3,vc)+ rsxy(i1-1,i2,i3,0,2)*u(i1-1,i2,i3,wc)
                            ! aup   = rsxy(i1+1,i2,i3,0,0)*u(i1+1,i2,i3,uc)+rsxy(i1+1,i2,i3,0,1)*u(i1+1,i2,i3,vc)+ rsxy(i1+1,i2,i3,0,2)*u(i1+1,i2,i3,wc)
                            au    = rsxy(i1  ,i2,i3,0,0)*uu(i1  ,i2,i3,
     & uc) + rsxy(i1  ,i2,i3,0,1) * uu(i1  ,i2,i3,vc)
                            aum   = rsxy(i1-1,i2,i3,0,0)*(u(i1-1,i2,i3,
     & uc) - gvU ) + rsxy(i1-1,i2,i3,0,1) * (u(i1-1,i2,i3,vc)- gvV )
                            aup   = rsxy(i1+1,i2,i3,0,0)*(u(i1+1,i2,i3,
     & uc) - gvU ) + rsxy(i1+1,i2,i3,0,1) * (u(i1+1,i2,i3,vc)- gvV )
                            s1=1
                            s2=0
                            s3=0
                            var = uc
                            vard= uc
                            drl = dr(0)
                              Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*u(
     & i1   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
                              Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    u(
     & i1-s1,i2-s2,i3-s3,var)
                              Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) - 
     & 2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
                              Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) - 
     &    u(i1   ,i2   ,i3   ,var)
                              Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,
     & i2-s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
                              Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var)
                              Amr   = Apl
                              Bmr   = Bpl
                              betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 
     & 1./4.*Bpl**2
                              betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 
     & 1./4.*Bpr**2
                              betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 
     & 1./4.*Bml**2
                              betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 
     & 1./4.*Bmr**2
                              ep     = 1e-15
                              ap1    = (1./2.)/(ep + betapl)**(4)
                              ap2    = (1./2.)/(ep + betapr)**(4)
                              am1    = (1./2.)/(ep + betaml)**(4)
                              am2    = (1./2.)/(ep + betamr)**(4)
                              !calculate weights
                              wp1 = ap1/(ap1+ap2)
                              wp2 = ap2/(ap1+ap2)
                              wm1 = am1/(am1+am2)
                              wm2 = am2/(am1+am2)
                              !mapping of the weights
                              bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                              bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                              bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                              bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                              wplw = bp1/(bp1+bp2)
                              wprw = bp2/(bp1+bp2)
                              wmlw = bm1/(bm1+bm2)
                              wmrw = bm2/(bm1+bm2)
                              !get face value of u
                              up = au
                              um = au
                              !picking side for the weights
                              if (up.ge.0) then
                                 wpl = max(wplw,wprw)
                                 wpr = min(wplw,wprw)
                              else
                                 wpl = min(wplw,wprw)
                                 wpr = max(wplw,wprw)
                              end if
                              if (um.ge.0) then
                                 wml = max(wmlw,wmrw)
                                 wmr = min(wmlw,wmrw)
                              else
                                 wml = min(wmlw,wmrw)
                                 wmr = max(wmlw,wmrw)
                              end if
                              expAd = 0.
                              if (orderOneFix.eq.1) then
                                 ! if (augmentPlot.eq.1) then
                                 !    agu(uc,13)=0
                                 ! end if
                                 if (aum*aup .le. 0) then
                                    aumax = max(abs(au),abs(aup),abs(
     & aum))
                                    expAd = (1./12.)*aumax*( u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(i1,
     & i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*s2,
     & i3+2*s3,var))
                                    wpl = 0.5
                                    wpr = 0.5
                                    wml = 0.5
                                    wmr = 0.5
                                    ! plot au at where fix is turned on
                                    ! if (augmentPlot.eq.1) then
                                    !    agu(uc,13)=au
                                    ! end if
                                 end if
                              end if
                              ! wpl = 0.5
                              ! wpr = 0.5
                              ! wml = 0.5
                              ! wmr = 0.5
                              Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,var) 
     & + 5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
                              Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 5.*
     & u(i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
                              Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*
     & s3,var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
                              Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,var)
     &  + 5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
                              Fp  = wpl*Fpl + wpr*Fpr
                              Fm  = wml*Fml + wmr*Fmr
                                 agu(vard,var)= au*(Fp - Fm)/drl + 
     & expAd/drl
                            !     if ( i1 .eq. 14 .and. i2.eq.16 ) then
                                    !     if ( i1 .eq. 14) then
                                    ! write(*,'(" betapl,betapr,betaml,betamr = ",4F10.8)')  betapl,betapr,betaml,betamr
                                    ! write(*,'(" ap1,ap2,am1,am2 = ",4e10.8)')  ap1,ap2,am1,am2
                                    ! write(*,'(" wp1,wpr,wm1,wm2 = ",4F10.8)')  wp1,wp2,wm1,wm2
                                    ! write(*,'(" bp1,bp2,bm1,bm2 = ",4F10.8)')  bp1,bp2,bm1,bm2
                                    ! write(*,'(" wplw,wprw,wmlw,wmrw = ",4F10.8)')  wplw,wprw,wmlw,wmrw
                                    ! write(*,'(" wpl,wpr,wml,wmr = ",4F10.8)')  wpl,wpr,wml,wmr
                                    ! write(*,'(" agu = ",F10.8)') agu(vard,var)
                            !     end if
                            !u*vx
                            s1=1
                            s2=0
                            s3=0
                            var = vc
                            vard= uc
                            drl = dr(0)
                              Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*u(
     & i1   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
                              Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    u(
     & i1-s1,i2-s2,i3-s3,var)
                              Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) - 
     & 2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
                              Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) - 
     &    u(i1   ,i2   ,i3   ,var)
                              Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,
     & i2-s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
                              Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var)
                              Amr   = Apl
                              Bmr   = Bpl
                              betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 
     & 1./4.*Bpl**2
                              betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 
     & 1./4.*Bpr**2
                              betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 
     & 1./4.*Bml**2
                              betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 
     & 1./4.*Bmr**2
                              ep     = 1e-15
                              ap1    = (1./2.)/(ep + betapl)**(4)
                              ap2    = (1./2.)/(ep + betapr)**(4)
                              am1    = (1./2.)/(ep + betaml)**(4)
                              am2    = (1./2.)/(ep + betamr)**(4)
                              !calculate weights
                              wp1 = ap1/(ap1+ap2)
                              wp2 = ap2/(ap1+ap2)
                              wm1 = am1/(am1+am2)
                              wm2 = am2/(am1+am2)
                              !mapping of the weights
                              bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                              bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                              bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                              bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                              wplw = bp1/(bp1+bp2)
                              wprw = bp2/(bp1+bp2)
                              wmlw = bm1/(bm1+bm2)
                              wmrw = bm2/(bm1+bm2)
                              !get face value of u
                              up = au
                              um = au
                              !picking side for the weights
                              if (up.ge.0) then
                                 wpl = max(wplw,wprw)
                                 wpr = min(wplw,wprw)
                              else
                                 wpl = min(wplw,wprw)
                                 wpr = max(wplw,wprw)
                              end if
                              if (um.ge.0) then
                                 wml = max(wmlw,wmrw)
                                 wmr = min(wmlw,wmrw)
                              else
                                 wml = min(wmlw,wmrw)
                                 wmr = max(wmlw,wmrw)
                              end if
                              expAd = 0.
                              if (orderOneFix.eq.1) then
                                 ! if (augmentPlot.eq.1) then
                                 !    agu(uc,13)=0
                                 ! end if
                                 if (aum*aup .le. 0) then
                                    aumax = max(abs(au),abs(aup),abs(
     & aum))
                                    expAd = (1./12.)*aumax*( u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(i1,
     & i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*s2,
     & i3+2*s3,var))
                                    wpl = 0.5
                                    wpr = 0.5
                                    wml = 0.5
                                    wmr = 0.5
                                    ! plot au at where fix is turned on
                                    ! if (augmentPlot.eq.1) then
                                    !    agu(uc,13)=au
                                    ! end if
                                 end if
                              end if
                              ! wpl = 0.5
                              ! wpr = 0.5
                              ! wml = 0.5
                              ! wmr = 0.5
                              Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,var) 
     & + 5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
                              Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 5.*
     & u(i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
                              Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*
     & s3,var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
                              Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,var)
     &  + 5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
                              Fp  = wpl*Fpl + wpr*Fpr
                              Fm  = wml*Fml + wmr*Fmr
                                 agu(vard,var)= au*(Fp - Fm)/drl + 
     & expAd/drl
                            !     if ( i1 .eq. 14 .and. i2.eq.16 ) then
                                    !     if ( i1 .eq. 14) then
                                    ! write(*,'(" betapl,betapr,betaml,betamr = ",4F10.8)')  betapl,betapr,betaml,betamr
                                    ! write(*,'(" ap1,ap2,am1,am2 = ",4e10.8)')  ap1,ap2,am1,am2
                                    ! write(*,'(" wp1,wpr,wm1,wm2 = ",4F10.8)')  wp1,wp2,wm1,wm2
                                    ! write(*,'(" bp1,bp2,bm1,bm2 = ",4F10.8)')  bp1,bp2,bm1,bm2
                                    ! write(*,'(" wplw,wprw,wmlw,wmrw = ",4F10.8)')  wplw,wprw,wmlw,wmrw
                                    ! write(*,'(" wpl,wpr,wml,wmr = ",4F10.8)')  wpl,wpr,wml,wmr
                                    ! write(*,'(" agu = ",F10.8)') agu(vard,var)
                            !     end if
                            !v*uy
                            au    = rsxy(i1,i2  ,i3,1,0)*uu(i1,i2  ,i3,
     & uc)+rsxy(i1,i2,i3,1,1)*uu(i1,i2,i3,vc)
                            aum   = rsxy(i1,i2-1,i3,1,0)*(u(i1,i2-1,i3,
     & uc) - gvU) + rsxy(i1,i2-1,i3,1,1)*(u(i1,i2-1,i3,vc) - gvV)
                            aup   = rsxy(i1,i2+1,i3,1,0)*(u(i1,i2+1,i3,
     & uc) - gvU) + rsxy(i1,i2+1,i3,1,1)*(u(i1,i2+1,i3,vc) - gvV)
                            s1=0
                            s2=1
                            s3=0
                            var = uc
                            vard= vc
                            drl = dr(1)
                              Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*u(
     & i1   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
                              Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    u(
     & i1-s1,i2-s2,i3-s3,var)
                              Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) - 
     & 2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
                              Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) - 
     &    u(i1   ,i2   ,i3   ,var)
                              Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,
     & i2-s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
                              Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var)
                              Amr   = Apl
                              Bmr   = Bpl
                              betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 
     & 1./4.*Bpl**2
                              betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 
     & 1./4.*Bpr**2
                              betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 
     & 1./4.*Bml**2
                              betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 
     & 1./4.*Bmr**2
                              ep     = 1e-15
                              ap1    = (1./2.)/(ep + betapl)**(4)
                              ap2    = (1./2.)/(ep + betapr)**(4)
                              am1    = (1./2.)/(ep + betaml)**(4)
                              am2    = (1./2.)/(ep + betamr)**(4)
                              !calculate weights
                              wp1 = ap1/(ap1+ap2)
                              wp2 = ap2/(ap1+ap2)
                              wm1 = am1/(am1+am2)
                              wm2 = am2/(am1+am2)
                              !mapping of the weights
                              bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                              bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                              bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                              bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                              wplw = bp1/(bp1+bp2)
                              wprw = bp2/(bp1+bp2)
                              wmlw = bm1/(bm1+bm2)
                              wmrw = bm2/(bm1+bm2)
                              !get face value of u
                              up = au
                              um = au
                              !picking side for the weights
                              if (up.ge.0) then
                                 wpl = max(wplw,wprw)
                                 wpr = min(wplw,wprw)
                              else
                                 wpl = min(wplw,wprw)
                                 wpr = max(wplw,wprw)
                              end if
                              if (um.ge.0) then
                                 wml = max(wmlw,wmrw)
                                 wmr = min(wmlw,wmrw)
                              else
                                 wml = min(wmlw,wmrw)
                                 wmr = max(wmlw,wmrw)
                              end if
                              expAd = 0.
                              if (orderOneFix.eq.1) then
                                 ! if (augmentPlot.eq.1) then
                                 !    agu(uc,13)=0
                                 ! end if
                                 if (aum*aup .le. 0) then
                                    aumax = max(abs(au),abs(aup),abs(
     & aum))
                                    expAd = (1./12.)*aumax*( u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(i1,
     & i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*s2,
     & i3+2*s3,var))
                                    wpl = 0.5
                                    wpr = 0.5
                                    wml = 0.5
                                    wmr = 0.5
                                    ! plot au at where fix is turned on
                                    ! if (augmentPlot.eq.1) then
                                    !    agu(uc,13)=au
                                    ! end if
                                 end if
                              end if
                              ! wpl = 0.5
                              ! wpr = 0.5
                              ! wml = 0.5
                              ! wmr = 0.5
                              Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,var) 
     & + 5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
                              Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 5.*
     & u(i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
                              Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*
     & s3,var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
                              Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,var)
     &  + 5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
                              Fp  = wpl*Fpl + wpr*Fpr
                              Fm  = wml*Fml + wmr*Fmr
                                 agu(vard,var)= au*(Fp - Fm)/drl + 
     & expAd/drl
                            !     if ( i1 .eq. 14 .and. i2.eq.16 ) then
                                    !     if ( i1 .eq. 14) then
                                    ! write(*,'(" betapl,betapr,betaml,betamr = ",4F10.8)')  betapl,betapr,betaml,betamr
                                    ! write(*,'(" ap1,ap2,am1,am2 = ",4e10.8)')  ap1,ap2,am1,am2
                                    ! write(*,'(" wp1,wpr,wm1,wm2 = ",4F10.8)')  wp1,wp2,wm1,wm2
                                    ! write(*,'(" bp1,bp2,bm1,bm2 = ",4F10.8)')  bp1,bp2,bm1,bm2
                                    ! write(*,'(" wplw,wprw,wmlw,wmrw = ",4F10.8)')  wplw,wprw,wmlw,wmrw
                                    ! write(*,'(" wpl,wpr,wml,wmr = ",4F10.8)')  wpl,wpr,wml,wmr
                                    ! write(*,'(" agu = ",F10.8)') agu(vard,var)
                            !     end if
                            !v*vy
                            s1=0
                            s2=1
                            s3=0
                            var = vc
                            vard= vc
                            drl = dr(1)
                              Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*u(
     & i1   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
                              Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    u(
     & i1-s1,i2-s2,i3-s3,var)
                              Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) - 
     & 2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
                              Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) - 
     &    u(i1   ,i2   ,i3   ,var)
                              Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,
     & i2-s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
                              Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var)
                              Amr   = Apl
                              Bmr   = Bpl
                              betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 
     & 1./4.*Bpl**2
                              betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 
     & 1./4.*Bpr**2
                              betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 
     & 1./4.*Bml**2
                              betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 
     & 1./4.*Bmr**2
                              ep     = 1e-15
                              ap1    = (1./2.)/(ep + betapl)**(4)
                              ap2    = (1./2.)/(ep + betapr)**(4)
                              am1    = (1./2.)/(ep + betaml)**(4)
                              am2    = (1./2.)/(ep + betamr)**(4)
                              !calculate weights
                              wp1 = ap1/(ap1+ap2)
                              wp2 = ap2/(ap1+ap2)
                              wm1 = am1/(am1+am2)
                              wm2 = am2/(am1+am2)
                              !mapping of the weights
                              bp1 = wp1*(3./4.+wp1*(wp1-1./2.))
                              bp2 = wp2*(3./4.+wp2*(wp2-1./2.))
                              bm1 = wm1*(3./4.+wm1*(wm1-1./2.))
                              bm2 = wm2*(3./4.+wm2*(wm2-1./2.))
                              wplw = bp1/(bp1+bp2)
                              wprw = bp2/(bp1+bp2)
                              wmlw = bm1/(bm1+bm2)
                              wmrw = bm2/(bm1+bm2)
                              !get face value of u
                              up = au
                              um = au
                              !picking side for the weights
                              if (up.ge.0) then
                                 wpl = max(wplw,wprw)
                                 wpr = min(wplw,wprw)
                              else
                                 wpl = min(wplw,wprw)
                                 wpr = max(wplw,wprw)
                              end if
                              if (um.ge.0) then
                                 wml = max(wmlw,wmrw)
                                 wmr = min(wmlw,wmrw)
                              else
                                 wml = min(wmlw,wmrw)
                                 wmr = max(wmlw,wmrw)
                              end if
                              expAd = 0.
                              if (orderOneFix.eq.1) then
                                 ! if (augmentPlot.eq.1) then
                                 !    agu(uc,13)=0
                                 ! end if
                                 if (aum*aup .le. 0) then
                                    aumax = max(abs(au),abs(aup),abs(
     & aum))
                                    expAd = (1./12.)*aumax*( u(i1-2*s1,
     & i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(i1,
     & i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*s2,
     & i3+2*s3,var))
                                    wpl = 0.5
                                    wpr = 0.5
                                    wml = 0.5
                                    wmr = 0.5
                                    ! plot au at where fix is turned on
                                    ! if (augmentPlot.eq.1) then
                                    !    agu(uc,13)=au
                                    ! end if
                                 end if
                              end if
                              ! wpl = 0.5
                              ! wpr = 0.5
                              ! wml = 0.5
                              ! wmr = 0.5
                              Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,var) 
     & + 5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
                              Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 5.*
     & u(i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
                              Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*
     & s3,var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
                              Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,var)
     &  + 5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
                              Fp  = wpl*Fpl + wpr*Fpr
                              Fm  = wml*Fml + wmr*Fmr
                                 agu(vard,var)= au*(Fp - Fm)/drl + 
     & expAd/drl
                            !     if ( i1 .eq. 14 .and. i2.eq.16 ) then
                                    !     if ( i1 .eq. 14) then
                                    ! write(*,'(" betapl,betapr,betaml,betamr = ",4F10.8)')  betapl,betapr,betaml,betamr
                                    ! write(*,'(" ap1,ap2,am1,am2 = ",4e10.8)')  ap1,ap2,am1,am2
                                    ! write(*,'(" wp1,wpr,wm1,wm2 = ",4F10.8)')  wp1,wp2,wm1,wm2
                                    ! write(*,'(" bp1,bp2,bm1,bm2 = ",4F10.8)')  bp1,bp2,bm1,bm2
                                    ! write(*,'(" wplw,wprw,wmlw,wmrw = ",4F10.8)')  wplw,wprw,wmlw,wmrw
                                    ! write(*,'(" wpl,wpr,wml,wmr = ",4F10.8)')  wpl,wpr,wml,wmr
                                    ! write(*,'(" agu = ",F10.8)') agu(vard,var)
                            !     end if
                         else ! if (upwindOrder.eq.4)
                            write(*,'(" getBewnoAdvection: only 4th 
     & order is avaliable now")' )
                            stop 777
                         end if
                     ! #If ("BWENO" == "CENTERED") 
                     !  ! -- centered advection ---
                     !  ! write(*,'(" getAdvection -- centered")')
                     !  #If "2" == "2"
                     !   ! -- 2D --
                     !   agu(uc,uc)=UU(uc)*UX(uc)
                     !   agu(vc,uc)=UU(vc)*UY(uc)
                     !   agu(uc,vc)=UU(uc)*UX(vc)
                     !   agu(vc,vc)=UU(vc)*UY(vc)
                     !  #Else
                     !   ! -- 3D -- *check me*
                     !   agu(uc,uc)=UU(uc)*UX(uc)
                     !   agu(vc,uc)=UU(vc)*UY(uc)
                     !   agu(wc,uc)=UU(wc)*UY(uc)
                     !   agu(uc,vc)=UU(uc)*UX(vc)
                     !   agu(vc,vc)=UU(vc)*UY(vc)
                     !   agu(wc,vc)=UU(wc)*UY(vc)
                     !   agu(uc,wc)=UU(uc)*UX(wc)
                     !   agu(vc,wc)=UU(vc)*UY(wc)
                     !   agu(wc,wc)=UU(wc)*UY(wc)
                     !  #End
                     ! #Elif "BWENO" == "BWENO" 
                     !   ! --- upwind scheme ---
                     !   ! for testing output this next message:
                     !   if( t.le. 0. )then
                     !     write(*,'(" getAdvection upwind scheme (7)")') 
                     !   end if
                     !   getUpwindAdvection(u,i1,i2,i3,NONE,2,4,curvilinear, agu)
                     ! #Elif "BWENO" == "BWENO" 
                     !   ! --- Bweno scheme ---
                     !   ! for testing output this message:
                     !   if( t.le. 0. )then
                     !      write(*,'(" getAdvection BWENO scheme (7)")') 
                     !   end if
                     !   getBwenoAdvection(u,i1,i2,i3,NONE,2,4,curvilinear, agu)
                     ! #Else
                     !   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
                     !   stop 999
                     ! #End
                       ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-ux42(
     & i1,i2,i3,pc)
                       ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-uy42(
     & i1,i2,i3,pc)
                end do
                end do
                end do
               end if
             else
              write(*,*)'insdt: Unknown implicitOption=',implicitOption
              stop 7
             end if  ! end implicitOption
            else
              write(*,'(" unknown advectionOption")')
              stop 1010
            end if
           end if
          else if( isAxisymmetric.eq.1 )then
           if( advectionOption.ne.centeredAdvection )then
             write(*,*) 'insdt.h : finish me for axisymmetric'
             stop 2020
           end if
           ! **** axisymmetric case ****
           if( gridIsImplicit.eq.0 )then
            ! explicit
            if( useWhereMask.ne.0 )then
             do i3=n3a,n3b
             do i2=n2a,n2b
             do i1=n1a,n1b
              if( mask(i1,i2,i3).gt.0 )then
                 ! INS, no AD
                  ut(i1,i2,i3,uc)= -uu(i1,i2,i3,uc)*ux42(i1,i2,i3,uc)-
     & uu(i1,i2,i3,vc)*uy42(i1,i2,i3,uc)-ux42(i1,i2,i3,pc)+nu*
     & ulaplacian42(i1,i2,i3,uc)
                  ut(i1,i2,i3,vc)= -uu(i1,i2,i3,uc)*ux42(i1,i2,i3,vc)-
     & uu(i1,i2,i3,vc)*uy42(i1,i2,i3,vc)-uy42(i1,i2,i3,pc)+nu*
     & ulaplacian42(i1,i2,i3,vc)
                  ! -- add on axisymmetric corrections ---
                  ri=radiusInverse(i1,i2,i3)
                  if( ri.ne.0. )then
                    ut(i1,i2,i3,uc)=ut(i1,i2,i3,uc)+nu*( uy42(i1,i2,i3,
     & uc)*ri )
                    ut(i1,i2,i3,vc)=ut(i1,i2,i3,vc)+nu*( (uy42(i1,i2,
     & i3,vc)-uu(i1,i2,i3,vc)*ri)*ri )
                  else
                    ut(i1,i2,i3,uc)=ut(i1,i2,i3,uc)+nu*( uyy42(i1,i2,
     & i3,uc) )
                    ut(i1,i2,i3,vc)=ut(i1,i2,i3,vc)+nu*( .5*uyy42(i1,
     & i2,i3,vc) )
                  end if
              end if
             end do
             end do
             end do
            else
             do i3=n3a,n3b
             do i2=n2a,n2b
             do i1=n1a,n1b
                ! INS, no AD
                 ut(i1,i2,i3,uc)= -uu(i1,i2,i3,uc)*ux42(i1,i2,i3,uc)-
     & uu(i1,i2,i3,vc)*uy42(i1,i2,i3,uc)-ux42(i1,i2,i3,pc)+nu*
     & ulaplacian42(i1,i2,i3,uc)
                 ut(i1,i2,i3,vc)= -uu(i1,i2,i3,uc)*ux42(i1,i2,i3,vc)-
     & uu(i1,i2,i3,vc)*uy42(i1,i2,i3,vc)-uy42(i1,i2,i3,pc)+nu*
     & ulaplacian42(i1,i2,i3,vc)
                 ! -- add on axisymmetric corrections ---
                 ri=radiusInverse(i1,i2,i3)
                 if( ri.ne.0. )then
                   ut(i1,i2,i3,uc)=ut(i1,i2,i3,uc)+nu*( uy42(i1,i2,i3,
     & uc)*ri )
                   ut(i1,i2,i3,vc)=ut(i1,i2,i3,vc)+nu*( (uy42(i1,i2,i3,
     & vc)-uu(i1,i2,i3,vc)*ri)*ri )
                 else
                   ut(i1,i2,i3,uc)=ut(i1,i2,i3,uc)+nu*( uyy42(i1,i2,i3,
     & uc) )
                   ut(i1,i2,i3,vc)=ut(i1,i2,i3,vc)+nu*( .5*uyy42(i1,i2,
     & i3,vc) )
                 end if
             end do
             end do
             end do
            end if
           else ! gridIsImplicit
            ! ***** implicit *******
            if( implicitOption .eq.computeImplicitTermsSeparately )then
              if( useWhereMask.ne.0 )then
               do i3=n3a,n3b
               do i2=n2a,n2b
               do i1=n1a,n1b
                if( mask(i1,i2,i3).gt.0 )then
                   ! explicit terms only, no diffusion
                    ut(i1,i2,i3,uc)= -uu(i1,i2,i3,uc)*ux42(i1,i2,i3,uc)
     & -uu(i1,i2,i3,vc)*uy42(i1,i2,i3,uc)-ux42(i1,i2,i3,pc)
                    ut(i1,i2,i3,vc)= -uu(i1,i2,i3,uc)*ux42(i1,i2,i3,vc)
     & -uu(i1,i2,i3,vc)*uy42(i1,i2,i3,vc)-uy42(i1,i2,i3,pc)
                   ! include implicit terms - diffusion
                    uti(i1,i2,i3,uc)= nu*ulaplacian42(i1,i2,i3,uc)
                    uti(i1,i2,i3,vc)= nu*ulaplacian42(i1,i2,i3,vc)
                    ri=radiusInverse(i1,i2,i3)
                    if( ri.ne.0. )then
                      uti(i1,i2,i3,uc)=uti(i1,i2,i3,uc)+nu*( uy42(i1,
     & i2,i3,uc)*ri )
                      uti(i1,i2,i3,vc)=uti(i1,i2,i3,vc)+nu*( (uy42(i1,
     & i2,i3,vc)-uu(i1,i2,i3,vc)*ri)*ri )
                    else
                      uti(i1,i2,i3,uc)=uti(i1,i2,i3,uc)+nu*( uyy42(i1,
     & i2,i3,uc) )
                      uti(i1,i2,i3,vc)=uti(i1,i2,i3,vc)+nu*( .5*uyy42(
     & i1,i2,i3,vc) )
                    end if
                end if
               end do
               end do
               end do
              else
               do i3=n3a,n3b
               do i2=n2a,n2b
               do i1=n1a,n1b
                  ! explicit terms only, no diffusion
                   ut(i1,i2,i3,uc)= -uu(i1,i2,i3,uc)*ux42(i1,i2,i3,uc)-
     & uu(i1,i2,i3,vc)*uy42(i1,i2,i3,uc)-ux42(i1,i2,i3,pc)
                   ut(i1,i2,i3,vc)= -uu(i1,i2,i3,uc)*ux42(i1,i2,i3,vc)-
     & uu(i1,i2,i3,vc)*uy42(i1,i2,i3,vc)-uy42(i1,i2,i3,pc)
                  ! include implicit terms - diffusion
                   uti(i1,i2,i3,uc)= nu*ulaplacian42(i1,i2,i3,uc)
                   uti(i1,i2,i3,vc)= nu*ulaplacian42(i1,i2,i3,vc)
                   ri=radiusInverse(i1,i2,i3)
                   if( ri.ne.0. )then
                     uti(i1,i2,i3,uc)=uti(i1,i2,i3,uc)+nu*( uy42(i1,i2,
     & i3,uc)*ri )
                     uti(i1,i2,i3,vc)=uti(i1,i2,i3,vc)+nu*( (uy42(i1,
     & i2,i3,vc)-uu(i1,i2,i3,vc)*ri)*ri )
                   else
                     uti(i1,i2,i3,uc)=uti(i1,i2,i3,uc)+nu*( uyy42(i1,
     & i2,i3,uc) )
                     uti(i1,i2,i3,vc)=uti(i1,i2,i3,vc)+nu*( .5*uyy42(
     & i1,i2,i3,vc) )
                   end if
               end do
               end do
               end do
              end if
            else if( implicitOption.eq.doNotComputeImplicitTerms )then
              if( useWhereMask.ne.0 )then
               do i3=n3a,n3b
               do i2=n2a,n2b
               do i1=n1a,n1b
                if( mask(i1,i2,i3).gt.0 )then
                   ! explicit terms only, no diffusion
                    ut(i1,i2,i3,uc)= -uu(i1,i2,i3,uc)*ux42(i1,i2,i3,uc)
     & -uu(i1,i2,i3,vc)*uy42(i1,i2,i3,uc)-ux42(i1,i2,i3,pc)
                    ut(i1,i2,i3,vc)= -uu(i1,i2,i3,uc)*ux42(i1,i2,i3,vc)
     & -uu(i1,i2,i3,vc)*uy42(i1,i2,i3,vc)-uy42(i1,i2,i3,pc)
                end if
               end do
               end do
               end do
              else
               do i3=n3a,n3b
               do i2=n2a,n2b
               do i1=n1a,n1b
                  ! explicit terms only, no diffusion
                   ut(i1,i2,i3,uc)= -uu(i1,i2,i3,uc)*ux42(i1,i2,i3,uc)-
     & uu(i1,i2,i3,vc)*uy42(i1,i2,i3,uc)-ux42(i1,i2,i3,pc)
                   ut(i1,i2,i3,vc)= -uu(i1,i2,i3,uc)*ux42(i1,i2,i3,vc)-
     & uu(i1,i2,i3,vc)*uy42(i1,i2,i3,vc)-uy42(i1,i2,i3,pc)
               end do
               end do
               end do
              end if
            else
             write(*,*)'insdt: Unknown implicitOption=',implicitOption
             stop 5
            end if  ! end implicitOption
           end if
          else
            stop 88733
          end if
         else
           stop 77
         end if
         return
         end
