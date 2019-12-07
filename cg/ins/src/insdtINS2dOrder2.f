! This file automatically generated from insdtINS.bf with bpp.
         subroutine insdtINS2dOrder2(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,
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
         d12(kd) = 1./(2.*dr(kd))
         d22(kd) = 1./(dr(kd)**2)
         ur2(i1,i2,i3,kd)=(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))*d12(0)
         us2(i1,i2,i3,kd)=(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))*d12(1)
         ut2(i1,i2,i3,kd)=(u(i1,i2,i3+1,kd)-u(i1,i2,i3-1,kd))*d12(2)
         urr2(i1,i2,i3,kd)=(-2.*u(i1,i2,i3,kd)+(u(i1+1,i2,i3,kd)+u(i1-
     & 1,i2,i3,kd)) )*d22(0)
         uss2(i1,i2,i3,kd)=(-2.*u(i1,i2,i3,kd)+(u(i1,i2+1,i3,kd)+u(i1,
     & i2-1,i3,kd)) )*d22(1)
         urs2(i1,i2,i3,kd)=(ur2(i1,i2+1,i3,kd)-ur2(i1,i2-1,i3,kd))*d12(
     & 1)
         utt2(i1,i2,i3,kd)=(-2.*u(i1,i2,i3,kd)+(u(i1,i2,i3+1,kd)+u(i1,
     & i2,i3-1,kd)) )*d22(2)
         urt2(i1,i2,i3,kd)=(ur2(i1,i2,i3+1,kd)-ur2(i1,i2,i3-1,kd))*d12(
     & 2)
         ust2(i1,i2,i3,kd)=(us2(i1,i2,i3+1,kd)-us2(i1,i2,i3-1,kd))*d12(
     & 2)
         urrr2(i1,i2,i3,kd)=(-2.*(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))+(
     & u(i1+2,i2,i3,kd)-u(i1-2,i2,i3,kd)) )*d22(0)*d12(0)
         usss2(i1,i2,i3,kd)=(-2.*(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))+(
     & u(i1,i2+2,i3,kd)-u(i1,i2-2,i3,kd)) )*d22(1)*d12(1)
         uttt2(i1,i2,i3,kd)=(-2.*(u(i1,i2,i3+1,kd)-u(i1,i2,i3-1,kd))+(
     & u(i1,i2,i3+2,kd)-u(i1,i2,i3-2,kd)) )*d22(2)*d12(2)
         rxr2(i1,i2,i3)=(rx(i1+1,i2,i3)-rx(i1-1,i2,i3))*d12(0)
         rxs2(i1,i2,i3)=(rx(i1,i2+1,i3)-rx(i1,i2-1,i3))*d12(1)
         rxt2(i1,i2,i3)=(rx(i1,i2,i3+1)-rx(i1,i2,i3-1))*d12(2)
         rxrr2(i1,i2,i3)=(-2.*rx(i1,i2,i3)+(rx(i1+1,i2,i3)+rx(i1-1,i2,
     & i3)) )*d22(0)
         rxss2(i1,i2,i3)=(-2.*rx(i1,i2,i3)+(rx(i1,i2+1,i3)+rx(i1,i2-1,
     & i3)) )*d22(1)
         rxrs2(i1,i2,i3)=(rxr2(i1,i2+1,i3)-rxr2(i1,i2-1,i3))*d12(1)
         ryr2(i1,i2,i3)=(ry(i1+1,i2,i3)-ry(i1-1,i2,i3))*d12(0)
         rys2(i1,i2,i3)=(ry(i1,i2+1,i3)-ry(i1,i2-1,i3))*d12(1)
         ryt2(i1,i2,i3)=(ry(i1,i2,i3+1)-ry(i1,i2,i3-1))*d12(2)
         ryrr2(i1,i2,i3)=(-2.*ry(i1,i2,i3)+(ry(i1+1,i2,i3)+ry(i1-1,i2,
     & i3)) )*d22(0)
         ryss2(i1,i2,i3)=(-2.*ry(i1,i2,i3)+(ry(i1,i2+1,i3)+ry(i1,i2-1,
     & i3)) )*d22(1)
         ryrs2(i1,i2,i3)=(ryr2(i1,i2+1,i3)-ryr2(i1,i2-1,i3))*d12(1)
         rzr2(i1,i2,i3)=(rz(i1+1,i2,i3)-rz(i1-1,i2,i3))*d12(0)
         rzs2(i1,i2,i3)=(rz(i1,i2+1,i3)-rz(i1,i2-1,i3))*d12(1)
         rzt2(i1,i2,i3)=(rz(i1,i2,i3+1)-rz(i1,i2,i3-1))*d12(2)
         rzrr2(i1,i2,i3)=(-2.*rz(i1,i2,i3)+(rz(i1+1,i2,i3)+rz(i1-1,i2,
     & i3)) )*d22(0)
         rzss2(i1,i2,i3)=(-2.*rz(i1,i2,i3)+(rz(i1,i2+1,i3)+rz(i1,i2-1,
     & i3)) )*d22(1)
         rzrs2(i1,i2,i3)=(rzr2(i1,i2+1,i3)-rzr2(i1,i2-1,i3))*d12(1)
         sxr2(i1,i2,i3)=(sx(i1+1,i2,i3)-sx(i1-1,i2,i3))*d12(0)
         sxs2(i1,i2,i3)=(sx(i1,i2+1,i3)-sx(i1,i2-1,i3))*d12(1)
         sxt2(i1,i2,i3)=(sx(i1,i2,i3+1)-sx(i1,i2,i3-1))*d12(2)
         sxrr2(i1,i2,i3)=(-2.*sx(i1,i2,i3)+(sx(i1+1,i2,i3)+sx(i1-1,i2,
     & i3)) )*d22(0)
         sxss2(i1,i2,i3)=(-2.*sx(i1,i2,i3)+(sx(i1,i2+1,i3)+sx(i1,i2-1,
     & i3)) )*d22(1)
         sxrs2(i1,i2,i3)=(sxr2(i1,i2+1,i3)-sxr2(i1,i2-1,i3))*d12(1)
         syr2(i1,i2,i3)=(sy(i1+1,i2,i3)-sy(i1-1,i2,i3))*d12(0)
         sys2(i1,i2,i3)=(sy(i1,i2+1,i3)-sy(i1,i2-1,i3))*d12(1)
         syt2(i1,i2,i3)=(sy(i1,i2,i3+1)-sy(i1,i2,i3-1))*d12(2)
         syrr2(i1,i2,i3)=(-2.*sy(i1,i2,i3)+(sy(i1+1,i2,i3)+sy(i1-1,i2,
     & i3)) )*d22(0)
         syss2(i1,i2,i3)=(-2.*sy(i1,i2,i3)+(sy(i1,i2+1,i3)+sy(i1,i2-1,
     & i3)) )*d22(1)
         syrs2(i1,i2,i3)=(syr2(i1,i2+1,i3)-syr2(i1,i2-1,i3))*d12(1)
         szr2(i1,i2,i3)=(sz(i1+1,i2,i3)-sz(i1-1,i2,i3))*d12(0)
         szs2(i1,i2,i3)=(sz(i1,i2+1,i3)-sz(i1,i2-1,i3))*d12(1)
         szt2(i1,i2,i3)=(sz(i1,i2,i3+1)-sz(i1,i2,i3-1))*d12(2)
         szrr2(i1,i2,i3)=(-2.*sz(i1,i2,i3)+(sz(i1+1,i2,i3)+sz(i1-1,i2,
     & i3)) )*d22(0)
         szss2(i1,i2,i3)=(-2.*sz(i1,i2,i3)+(sz(i1,i2+1,i3)+sz(i1,i2-1,
     & i3)) )*d22(1)
         szrs2(i1,i2,i3)=(szr2(i1,i2+1,i3)-szr2(i1,i2-1,i3))*d12(1)
         txr2(i1,i2,i3)=(tx(i1+1,i2,i3)-tx(i1-1,i2,i3))*d12(0)
         txs2(i1,i2,i3)=(tx(i1,i2+1,i3)-tx(i1,i2-1,i3))*d12(1)
         txt2(i1,i2,i3)=(tx(i1,i2,i3+1)-tx(i1,i2,i3-1))*d12(2)
         txrr2(i1,i2,i3)=(-2.*tx(i1,i2,i3)+(tx(i1+1,i2,i3)+tx(i1-1,i2,
     & i3)) )*d22(0)
         txss2(i1,i2,i3)=(-2.*tx(i1,i2,i3)+(tx(i1,i2+1,i3)+tx(i1,i2-1,
     & i3)) )*d22(1)
         txrs2(i1,i2,i3)=(txr2(i1,i2+1,i3)-txr2(i1,i2-1,i3))*d12(1)
         tyr2(i1,i2,i3)=(ty(i1+1,i2,i3)-ty(i1-1,i2,i3))*d12(0)
         tys2(i1,i2,i3)=(ty(i1,i2+1,i3)-ty(i1,i2-1,i3))*d12(1)
         tyt2(i1,i2,i3)=(ty(i1,i2,i3+1)-ty(i1,i2,i3-1))*d12(2)
         tyrr2(i1,i2,i3)=(-2.*ty(i1,i2,i3)+(ty(i1+1,i2,i3)+ty(i1-1,i2,
     & i3)) )*d22(0)
         tyss2(i1,i2,i3)=(-2.*ty(i1,i2,i3)+(ty(i1,i2+1,i3)+ty(i1,i2-1,
     & i3)) )*d22(1)
         tyrs2(i1,i2,i3)=(tyr2(i1,i2+1,i3)-tyr2(i1,i2-1,i3))*d12(1)
         tzr2(i1,i2,i3)=(tz(i1+1,i2,i3)-tz(i1-1,i2,i3))*d12(0)
         tzs2(i1,i2,i3)=(tz(i1,i2+1,i3)-tz(i1,i2-1,i3))*d12(1)
         tzt2(i1,i2,i3)=(tz(i1,i2,i3+1)-tz(i1,i2,i3-1))*d12(2)
         tzrr2(i1,i2,i3)=(-2.*tz(i1,i2,i3)+(tz(i1+1,i2,i3)+tz(i1-1,i2,
     & i3)) )*d22(0)
         tzss2(i1,i2,i3)=(-2.*tz(i1,i2,i3)+(tz(i1,i2+1,i3)+tz(i1,i2-1,
     & i3)) )*d22(1)
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
         rxx22(i1,i2,i3)= rx(i1,i2,i3)*rxr2(i1,i2,i3)+sx(i1,i2,i3)*
     & rxs2(i1,i2,i3)
         rxy22(i1,i2,i3)= ry(i1,i2,i3)*rxr2(i1,i2,i3)+sy(i1,i2,i3)*
     & rxs2(i1,i2,i3)
         rxx23(i1,i2,i3)=rx(i1,i2,i3)*rxr2(i1,i2,i3)+sx(i1,i2,i3)*rxs2(
     & i1,i2,i3)+tx(i1,i2,i3)*rxt2(i1,i2,i3)
         rxy23(i1,i2,i3)=ry(i1,i2,i3)*rxr2(i1,i2,i3)+sy(i1,i2,i3)*rxs2(
     & i1,i2,i3)+ty(i1,i2,i3)*rxt2(i1,i2,i3)
         rxz23(i1,i2,i3)=rz(i1,i2,i3)*rxr2(i1,i2,i3)+sz(i1,i2,i3)*rxs2(
     & i1,i2,i3)+tz(i1,i2,i3)*rxt2(i1,i2,i3)
         ryx22(i1,i2,i3)= rx(i1,i2,i3)*ryr2(i1,i2,i3)+sx(i1,i2,i3)*
     & rys2(i1,i2,i3)
         ryy22(i1,i2,i3)= ry(i1,i2,i3)*ryr2(i1,i2,i3)+sy(i1,i2,i3)*
     & rys2(i1,i2,i3)
         ryx23(i1,i2,i3)=rx(i1,i2,i3)*ryr2(i1,i2,i3)+sx(i1,i2,i3)*rys2(
     & i1,i2,i3)+tx(i1,i2,i3)*ryt2(i1,i2,i3)
         ryy23(i1,i2,i3)=ry(i1,i2,i3)*ryr2(i1,i2,i3)+sy(i1,i2,i3)*rys2(
     & i1,i2,i3)+ty(i1,i2,i3)*ryt2(i1,i2,i3)
         ryz23(i1,i2,i3)=rz(i1,i2,i3)*ryr2(i1,i2,i3)+sz(i1,i2,i3)*rys2(
     & i1,i2,i3)+tz(i1,i2,i3)*ryt2(i1,i2,i3)
         rzx22(i1,i2,i3)= rx(i1,i2,i3)*rzr2(i1,i2,i3)+sx(i1,i2,i3)*
     & rzs2(i1,i2,i3)
         rzy22(i1,i2,i3)= ry(i1,i2,i3)*rzr2(i1,i2,i3)+sy(i1,i2,i3)*
     & rzs2(i1,i2,i3)
         rzx23(i1,i2,i3)=rx(i1,i2,i3)*rzr2(i1,i2,i3)+sx(i1,i2,i3)*rzs2(
     & i1,i2,i3)+tx(i1,i2,i3)*rzt2(i1,i2,i3)
         rzy23(i1,i2,i3)=ry(i1,i2,i3)*rzr2(i1,i2,i3)+sy(i1,i2,i3)*rzs2(
     & i1,i2,i3)+ty(i1,i2,i3)*rzt2(i1,i2,i3)
         rzz23(i1,i2,i3)=rz(i1,i2,i3)*rzr2(i1,i2,i3)+sz(i1,i2,i3)*rzs2(
     & i1,i2,i3)+tz(i1,i2,i3)*rzt2(i1,i2,i3)
         sxx22(i1,i2,i3)= rx(i1,i2,i3)*sxr2(i1,i2,i3)+sx(i1,i2,i3)*
     & sxs2(i1,i2,i3)
         sxy22(i1,i2,i3)= ry(i1,i2,i3)*sxr2(i1,i2,i3)+sy(i1,i2,i3)*
     & sxs2(i1,i2,i3)
         sxx23(i1,i2,i3)=rx(i1,i2,i3)*sxr2(i1,i2,i3)+sx(i1,i2,i3)*sxs2(
     & i1,i2,i3)+tx(i1,i2,i3)*sxt2(i1,i2,i3)
         sxy23(i1,i2,i3)=ry(i1,i2,i3)*sxr2(i1,i2,i3)+sy(i1,i2,i3)*sxs2(
     & i1,i2,i3)+ty(i1,i2,i3)*sxt2(i1,i2,i3)
         sxz23(i1,i2,i3)=rz(i1,i2,i3)*sxr2(i1,i2,i3)+sz(i1,i2,i3)*sxs2(
     & i1,i2,i3)+tz(i1,i2,i3)*sxt2(i1,i2,i3)
         syx22(i1,i2,i3)= rx(i1,i2,i3)*syr2(i1,i2,i3)+sx(i1,i2,i3)*
     & sys2(i1,i2,i3)
         syy22(i1,i2,i3)= ry(i1,i2,i3)*syr2(i1,i2,i3)+sy(i1,i2,i3)*
     & sys2(i1,i2,i3)
         syx23(i1,i2,i3)=rx(i1,i2,i3)*syr2(i1,i2,i3)+sx(i1,i2,i3)*sys2(
     & i1,i2,i3)+tx(i1,i2,i3)*syt2(i1,i2,i3)
         syy23(i1,i2,i3)=ry(i1,i2,i3)*syr2(i1,i2,i3)+sy(i1,i2,i3)*sys2(
     & i1,i2,i3)+ty(i1,i2,i3)*syt2(i1,i2,i3)
         syz23(i1,i2,i3)=rz(i1,i2,i3)*syr2(i1,i2,i3)+sz(i1,i2,i3)*sys2(
     & i1,i2,i3)+tz(i1,i2,i3)*syt2(i1,i2,i3)
         szx22(i1,i2,i3)= rx(i1,i2,i3)*szr2(i1,i2,i3)+sx(i1,i2,i3)*
     & szs2(i1,i2,i3)
         szy22(i1,i2,i3)= ry(i1,i2,i3)*szr2(i1,i2,i3)+sy(i1,i2,i3)*
     & szs2(i1,i2,i3)
         szx23(i1,i2,i3)=rx(i1,i2,i3)*szr2(i1,i2,i3)+sx(i1,i2,i3)*szs2(
     & i1,i2,i3)+tx(i1,i2,i3)*szt2(i1,i2,i3)
         szy23(i1,i2,i3)=ry(i1,i2,i3)*szr2(i1,i2,i3)+sy(i1,i2,i3)*szs2(
     & i1,i2,i3)+ty(i1,i2,i3)*szt2(i1,i2,i3)
         szz23(i1,i2,i3)=rz(i1,i2,i3)*szr2(i1,i2,i3)+sz(i1,i2,i3)*szs2(
     & i1,i2,i3)+tz(i1,i2,i3)*szt2(i1,i2,i3)
         txx22(i1,i2,i3)= rx(i1,i2,i3)*txr2(i1,i2,i3)+sx(i1,i2,i3)*
     & txs2(i1,i2,i3)
         txy22(i1,i2,i3)= ry(i1,i2,i3)*txr2(i1,i2,i3)+sy(i1,i2,i3)*
     & txs2(i1,i2,i3)
         txx23(i1,i2,i3)=rx(i1,i2,i3)*txr2(i1,i2,i3)+sx(i1,i2,i3)*txs2(
     & i1,i2,i3)+tx(i1,i2,i3)*txt2(i1,i2,i3)
         txy23(i1,i2,i3)=ry(i1,i2,i3)*txr2(i1,i2,i3)+sy(i1,i2,i3)*txs2(
     & i1,i2,i3)+ty(i1,i2,i3)*txt2(i1,i2,i3)
         txz23(i1,i2,i3)=rz(i1,i2,i3)*txr2(i1,i2,i3)+sz(i1,i2,i3)*txs2(
     & i1,i2,i3)+tz(i1,i2,i3)*txt2(i1,i2,i3)
         tyx22(i1,i2,i3)= rx(i1,i2,i3)*tyr2(i1,i2,i3)+sx(i1,i2,i3)*
     & tys2(i1,i2,i3)
         tyy22(i1,i2,i3)= ry(i1,i2,i3)*tyr2(i1,i2,i3)+sy(i1,i2,i3)*
     & tys2(i1,i2,i3)
         tyx23(i1,i2,i3)=rx(i1,i2,i3)*tyr2(i1,i2,i3)+sx(i1,i2,i3)*tys2(
     & i1,i2,i3)+tx(i1,i2,i3)*tyt2(i1,i2,i3)
         tyy23(i1,i2,i3)=ry(i1,i2,i3)*tyr2(i1,i2,i3)+sy(i1,i2,i3)*tys2(
     & i1,i2,i3)+ty(i1,i2,i3)*tyt2(i1,i2,i3)
         tyz23(i1,i2,i3)=rz(i1,i2,i3)*tyr2(i1,i2,i3)+sz(i1,i2,i3)*tys2(
     & i1,i2,i3)+tz(i1,i2,i3)*tyt2(i1,i2,i3)
         tzx22(i1,i2,i3)= rx(i1,i2,i3)*tzr2(i1,i2,i3)+sx(i1,i2,i3)*
     & tzs2(i1,i2,i3)
         tzy22(i1,i2,i3)= ry(i1,i2,i3)*tzr2(i1,i2,i3)+sy(i1,i2,i3)*
     & tzs2(i1,i2,i3)
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
         uxy22(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*urr2(i1,i2,i3,kd)
     & +(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*urs2(i1,
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
         uxy23(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*urr2(i1,i2,i3,kd)
     & +sx(i1,i2,i3)*sy(i1,i2,i3)*uss2(i1,i2,i3,kd)+tx(i1,i2,i3)*ty(
     & i1,i2,i3)*utt2(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,
     & i2,i3)*sx(i1,i2,i3))*urs2(i1,i2,i3,kd)+(rx(i1,i2,i3)*ty(i1,i2,
     & i3)+ry(i1,i2,i3)*tx(i1,i2,i3))*urt2(i1,i2,i3,kd)+(sx(i1,i2,i3)*
     & ty(i1,i2,i3)+sy(i1,i2,i3)*tx(i1,i2,i3))*ust2(i1,i2,i3,kd)+
     & rxy23(i1,i2,i3)*ur2(i1,i2,i3,kd)+sxy23(i1,i2,i3)*us2(i1,i2,i3,
     & kd)+txy23(i1,i2,i3)*ut2(i1,i2,i3,kd)
         uxz23(i1,i2,i3,kd)=rx(i1,i2,i3)*rz(i1,i2,i3)*urr2(i1,i2,i3,kd)
     & +sx(i1,i2,i3)*sz(i1,i2,i3)*uss2(i1,i2,i3,kd)+tx(i1,i2,i3)*tz(
     & i1,i2,i3)*utt2(i1,i2,i3,kd)+(rx(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,
     & i2,i3)*sx(i1,i2,i3))*urs2(i1,i2,i3,kd)+(rx(i1,i2,i3)*tz(i1,i2,
     & i3)+rz(i1,i2,i3)*tx(i1,i2,i3))*urt2(i1,i2,i3,kd)+(sx(i1,i2,i3)*
     & tz(i1,i2,i3)+sz(i1,i2,i3)*tx(i1,i2,i3))*ust2(i1,i2,i3,kd)+
     & rxz23(i1,i2,i3)*ur2(i1,i2,i3,kd)+sxz23(i1,i2,i3)*us2(i1,i2,i3,
     & kd)+txz23(i1,i2,i3)*ut2(i1,i2,i3,kd)
         uyz23(i1,i2,i3,kd)=ry(i1,i2,i3)*rz(i1,i2,i3)*urr2(i1,i2,i3,kd)
     & +sy(i1,i2,i3)*sz(i1,i2,i3)*uss2(i1,i2,i3,kd)+ty(i1,i2,i3)*tz(
     & i1,i2,i3)*utt2(i1,i2,i3,kd)+(ry(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,
     & i2,i3)*sy(i1,i2,i3))*urs2(i1,i2,i3,kd)+(ry(i1,i2,i3)*tz(i1,i2,
     & i3)+rz(i1,i2,i3)*ty(i1,i2,i3))*urt2(i1,i2,i3,kd)+(sy(i1,i2,i3)*
     & tz(i1,i2,i3)+sz(i1,i2,i3)*ty(i1,i2,i3))*ust2(i1,i2,i3,kd)+
     & ryz23(i1,i2,i3)*ur2(i1,i2,i3,kd)+syz23(i1,i2,i3)*us2(i1,i2,i3,
     & kd)+tyz23(i1,i2,i3)*ut2(i1,i2,i3,kd)
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
         uxx23r(i1,i2,i3,kd)=(-2.*u(i1,i2,i3,kd)+(u(i1+1,i2,i3,kd)+u(
     & i1-1,i2,i3,kd)) )*h22(0)
         uyy23r(i1,i2,i3,kd)=(-2.*u(i1,i2,i3,kd)+(u(i1,i2+1,i3,kd)+u(
     & i1,i2-1,i3,kd)) )*h22(1)
         uxy23r(i1,i2,i3,kd)=(ux23r(i1,i2+1,i3,kd)-ux23r(i1,i2-1,i3,kd)
     & )*h12(1)
         uzz23r(i1,i2,i3,kd)=(-2.*u(i1,i2,i3,kd)+(u(i1,i2,i3+1,kd)+u(
     & i1,i2,i3-1,kd)) )*h22(2)
         uxz23r(i1,i2,i3,kd)=(ux23r(i1,i2,i3+1,kd)-ux23r(i1,i2,i3-1,kd)
     & )*h12(2)
         uyz23r(i1,i2,i3,kd)=(uy23r(i1,i2,i3+1,kd)-uy23r(i1,i2,i3-1,kd)
     & )*h12(2)
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
         ulaplacian22r(i1,i2,i3,kd)=uxx23r(i1,i2,i3,kd)+uyy23r(i1,i2,
     & i3,kd)
         ulaplacian23r(i1,i2,i3,kd)=uxx23r(i1,i2,i3,kd)+uyy23r(i1,i2,
     & i3,kd)+uzz23r(i1,i2,i3,kd)
         uxxx22r(i1,i2,i3,kd)=(-2.*(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))+
     & (u(i1+2,i2,i3,kd)-u(i1-2,i2,i3,kd)) )*h22(0)*h12(0)
         uyyy22r(i1,i2,i3,kd)=(-2.*(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))+
     & (u(i1,i2+2,i3,kd)-u(i1,i2-2,i3,kd)) )*h22(1)*h12(1)
         uxxy22r(i1,i2,i3,kd)=( uxx22r(i1,i2+1,i3,kd)-uxx22r(i1,i2-1,
     & i3,kd))/(2.*dx(1))
         uxyy22r(i1,i2,i3,kd)=( uyy22r(i1+1,i2,i3,kd)-uyy22r(i1-1,i2,
     & i3,kd))/(2.*dx(0))
         uxxxx22r(i1,i2,i3,kd)=(6.*u(i1,i2,i3,kd)-4.*(u(i1+1,i2,i3,kd)+
     & u(i1-1,i2,i3,kd))+(u(i1+2,i2,i3,kd)+u(i1-2,i2,i3,kd)) )/(dx(0)*
     & *4)
         uyyyy22r(i1,i2,i3,kd)=(6.*u(i1,i2,i3,kd)-4.*(u(i1,i2+1,i3,kd)+
     & u(i1,i2-1,i3,kd))+(u(i1,i2+2,i3,kd)+u(i1,i2-2,i3,kd)) )/(dx(1)*
     & *4)
         uxxyy22r(i1,i2,i3,kd)=( 4.*u(i1,i2,i3,kd)     -2.*(u(i1+1,i2,
     & i3,kd)+u(i1-1,i2,i3,kd)+u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd))   + 
     &   (u(i1+1,i2+1,i3,kd)+u(i1-1,i2+1,i3,kd)+u(i1+1,i2-1,i3,kd)+u(
     & i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)
         ! 2D laplacian squared = u.xxxx + 2 u.xxyy + u.yyyy
         uLapSq22r(i1,i2,i3,kd)= ( 6.*u(i1,i2,i3,kd)   - 4.*(u(i1+1,i2,
     & i3,kd)+u(i1-1,i2,i3,kd))    +(u(i1+2,i2,i3,kd)+u(i1-2,i2,i3,kd)
     & ) )/(dx(0)**4) +( 6.*u(i1,i2,i3,kd)    -4.*(u(i1,i2+1,i3,kd)+u(
     & i1,i2-1,i3,kd))    +(u(i1,i2+2,i3,kd)+u(i1,i2-2,i3,kd)) )/(dx(
     & 1)**4)  +( 8.*u(i1,i2,i3,kd)     -4.*(u(i1+1,i2,i3,kd)+u(i1-1,
     & i2,i3,kd)+u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd))   +2.*(u(i1+1,i2+
     & 1,i3,kd)+u(i1-1,i2+1,i3,kd)+u(i1+1,i2-1,i3,kd)+u(i1-1,i2-1,i3,
     & kd)) )/(dx(0)**2*dx(1)**2)
         uxxx23r(i1,i2,i3,kd)=(-2.*(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))+
     & (u(i1+2,i2,i3,kd)-u(i1-2,i2,i3,kd)) )*h22(0)*h12(0)
         uyyy23r(i1,i2,i3,kd)=(-2.*(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))+
     & (u(i1,i2+2,i3,kd)-u(i1,i2-2,i3,kd)) )*h22(1)*h12(1)
         uzzz23r(i1,i2,i3,kd)=(-2.*(u(i1,i2,i3+1,kd)-u(i1,i2,i3-1,kd))+
     & (u(i1,i2,i3+2,kd)-u(i1,i2,i3-2,kd)) )*h22(1)*h12(2)
         uxxy23r(i1,i2,i3,kd)=( uxx22r(i1,i2+1,i3,kd)-uxx22r(i1,i2-1,
     & i3,kd))/(2.*dx(1))
         uxyy23r(i1,i2,i3,kd)=( uyy22r(i1+1,i2,i3,kd)-uyy22r(i1-1,i2,
     & i3,kd))/(2.*dx(0))
         uxxz23r(i1,i2,i3,kd)=( uxx22r(i1,i2,i3+1,kd)-uxx22r(i1,i2,i3-
     & 1,kd))/(2.*dx(2))
         uyyz23r(i1,i2,i3,kd)=( uyy22r(i1,i2,i3+1,kd)-uyy22r(i1,i2,i3-
     & 1,kd))/(2.*dx(2))
         uxzz23r(i1,i2,i3,kd)=( uzz22r(i1+1,i2,i3,kd)-uzz22r(i1-1,i2,
     & i3,kd))/(2.*dx(0))
         uyzz23r(i1,i2,i3,kd)=( uzz22r(i1,i2+1,i3,kd)-uzz22r(i1,i2-1,
     & i3,kd))/(2.*dx(1))
         uxxxx23r(i1,i2,i3,kd)=(6.*u(i1,i2,i3,kd)-4.*(u(i1+1,i2,i3,kd)+
     & u(i1-1,i2,i3,kd))+(u(i1+2,i2,i3,kd)+u(i1-2,i2,i3,kd)) )/(dx(0)*
     & *4)
         uyyyy23r(i1,i2,i3,kd)=(6.*u(i1,i2,i3,kd)-4.*(u(i1,i2+1,i3,kd)+
     & u(i1,i2-1,i3,kd))+(u(i1,i2+2,i3,kd)+u(i1,i2-2,i3,kd)) )/(dx(1)*
     & *4)
         uzzzz23r(i1,i2,i3,kd)=(6.*u(i1,i2,i3,kd)-4.*(u(i1,i2,i3+1,kd)+
     & u(i1,i2,i3-1,kd))+(u(i1,i2,i3+2,kd)+u(i1,i2,i3-2,kd)) )/(dx(2)*
     & *4)
         uxxyy23r(i1,i2,i3,kd)=( 4.*u(i1,i2,i3,kd)     -2.*(u(i1+1,i2,
     & i3,kd)+u(i1-1,i2,i3,kd)+u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd))   + 
     &   (u(i1+1,i2+1,i3,kd)+u(i1-1,i2+1,i3,kd)+u(i1+1,i2-1,i3,kd)+u(
     & i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)
         uxxzz23r(i1,i2,i3,kd)=( 4.*u(i1,i2,i3,kd)     -2.*(u(i1+1,i2,
     & i3,kd)+u(i1-1,i2,i3,kd)+u(i1,i2,i3+1,kd)+u(i1,i2,i3-1,kd))   + 
     &   (u(i1+1,i2,i3+1,kd)+u(i1-1,i2,i3+1,kd)+u(i1+1,i2,i3-1,kd)+u(
     & i1-1,i2,i3-1,kd)) )/(dx(0)**2*dx(2)**2)
         uyyzz23r(i1,i2,i3,kd)=( 4.*u(i1,i2,i3,kd)     -2.*(u(i1,i2+1,
     & i3,kd)  +u(i1,i2-1,i3,kd)+  u(i1,i2  ,i3+1,kd)+u(i1,i2  ,i3-1,
     & kd))   +   (u(i1,i2+1,i3+1,kd)+u(i1,i2-1,i3+1,kd)+u(i1,i2+1,i3-
     & 1,kd)+u(i1,i2-1,i3-1,kd)) )/(dx(1)**2*dx(2)**2)
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
                    ut(i1,i2,i3,uc)= -uu(i1,i2,i3,uc)*ux22r(i1,i2,i3,
     & uc)-uu(i1,i2,i3,vc)*uy22r(i1,i2,i3,uc)-ux22r(i1,i2,i3,pc)+nu*
     & ulaplacian22r(i1,i2,i3,uc)
                    ut(i1,i2,i3,vc)= -uu(i1,i2,i3,uc)*ux22r(i1,i2,i3,
     & vc)-uu(i1,i2,i3,vc)*uy22r(i1,i2,i3,vc)-uy22r(i1,i2,i3,pc)+nu*
     & ulaplacian22r(i1,i2,i3,vc)
                end if
               end do
               end do
               end do
              else
               do i3=n3a,n3b
               do i2=n2a,n2b
               do i1=n1a,n1b
                  ! INS, no AD
                   ut(i1,i2,i3,uc)= -uu(i1,i2,i3,uc)*ux22r(i1,i2,i3,uc)
     & -uu(i1,i2,i3,vc)*uy22r(i1,i2,i3,uc)-ux22r(i1,i2,i3,pc)+nu*
     & ulaplacian22r(i1,i2,i3,uc)
                   ut(i1,i2,i3,vc)= -uu(i1,i2,i3,uc)*ux22r(i1,i2,i3,vc)
     & -uu(i1,i2,i3,vc)*uy22r(i1,i2,i3,vc)-uy22r(i1,i2,i3,pc)+nu*
     & ulaplacian22r(i1,i2,i3,vc)
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
                    !   getUpwindAdvection(u,i1,i2,i3,NONE,2,2,rectangular, agu)
                    ! #Elif "UPWIND" == "BWENO" 
                    !   ! --- Bweno scheme ---
                    !   ! for testing output this message:
                    !   if( t.le. 0. )then
                    !      write(*,'(" getAdvection BWENO scheme (7)")') 
                    !   end if
                    !   getBwenoAdvection(u,i1,i2,i3,NONE,2,2,rectangular, agu)
                    ! #Else
                    !   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
                    !   stop 999
                    ! #End
                       ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-ux22r(
     & i1,i2,i3,pc)+nu*ulaplacian22r(i1,i2,i3,uc)
                       ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-uy22r(
     & i1,i2,i3,pc)+nu*ulaplacian22r(i1,i2,i3,vc)
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
                   !   getUpwindAdvection(u,i1,i2,i3,NONE,2,2,rectangular, agu)
                   ! #Elif "UPWIND" == "BWENO" 
                   !   ! --- Bweno scheme ---
                   !   ! for testing output this message:
                   !   if( t.le. 0. )then
                   !      write(*,'(" getAdvection BWENO scheme (7)")') 
                   !   end if
                   !   getBwenoAdvection(u,i1,i2,i3,NONE,2,2,rectangular, agu)
                   ! #Else
                   !   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
                   !   stop 999
                   ! #End
                      ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-ux22r(
     & i1,i2,i3,pc)+nu*ulaplacian22r(i1,i2,i3,uc)
                      ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-uy22r(
     & i1,i2,i3,pc)+nu*ulaplacian22r(i1,i2,i3,vc)
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
                    !   getUpwindAdvection(u,i1,i2,i3,NONE,2,2,rectangular, agu)
                    ! #Elif "BWENO" == "BWENO" 
                    !   ! --- Bweno scheme ---
                    !   ! for testing output this message:
                    !   if( t.le. 0. )then
                    !      write(*,'(" getAdvection BWENO scheme (7)")') 
                    !   end if
                    !   getBwenoAdvection(u,i1,i2,i3,NONE,2,2,rectangular, agu)
                    ! #Else
                    !   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
                    !   stop 999
                    ! #End
                       ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-ux22r(
     & i1,i2,i3,pc)+nu*ulaplacian22r(i1,i2,i3,uc)
                       ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-uy22r(
     & i1,i2,i3,pc)+nu*ulaplacian22r(i1,i2,i3,vc)
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
                   !   getUpwindAdvection(u,i1,i2,i3,NONE,2,2,rectangular, agu)
                   ! #Elif "BWENO" == "BWENO" 
                   !   ! --- Bweno scheme ---
                   !   ! for testing output this message:
                   !   if( t.le. 0. )then
                   !      write(*,'(" getAdvection BWENO scheme (7)")') 
                   !   end if
                   !   getBwenoAdvection(u,i1,i2,i3,NONE,2,2,rectangular, agu)
                   ! #Else
                   !   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
                   !   stop 999
                   ! #End
                      ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-ux22r(
     & i1,i2,i3,pc)+nu*ulaplacian22r(i1,i2,i3,uc)
                      ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-uy22r(
     & i1,i2,i3,pc)+nu*ulaplacian22r(i1,i2,i3,vc)
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
                     ut(i1,i2,i3,uc)= -uu(i1,i2,i3,uc)*ux22r(i1,i2,i3,
     & uc)-uu(i1,i2,i3,vc)*uy22r(i1,i2,i3,uc)-ux22r(i1,i2,i3,pc)
                     ut(i1,i2,i3,vc)= -uu(i1,i2,i3,uc)*ux22r(i1,i2,i3,
     & vc)-uu(i1,i2,i3,vc)*uy22r(i1,i2,i3,vc)-uy22r(i1,i2,i3,pc)
                    ! include implicit terms - diffusion
                     uti(i1,i2,i3,uc)= nu*ulaplacian22r(i1,i2,i3,uc)
                     uti(i1,i2,i3,vc)= nu*ulaplacian22r(i1,i2,i3,vc)
                 end if
                end do
                end do
                end do
               else
                do i3=n3a,n3b
                do i2=n2a,n2b
                do i1=n1a,n1b
                   ! explicit terms only, no diffusion
                    ut(i1,i2,i3,uc)= -uu(i1,i2,i3,uc)*ux22r(i1,i2,i3,
     & uc)-uu(i1,i2,i3,vc)*uy22r(i1,i2,i3,uc)-ux22r(i1,i2,i3,pc)
                    ut(i1,i2,i3,vc)= -uu(i1,i2,i3,uc)*ux22r(i1,i2,i3,
     & vc)-uu(i1,i2,i3,vc)*uy22r(i1,i2,i3,vc)-uy22r(i1,i2,i3,pc)
                   ! include implicit terms - diffusion
                    uti(i1,i2,i3,uc)= nu*ulaplacian22r(i1,i2,i3,uc)
                    uti(i1,i2,i3,vc)= nu*ulaplacian22r(i1,i2,i3,vc)
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
                     ut(i1,i2,i3,uc)= -uu(i1,i2,i3,uc)*ux22r(i1,i2,i3,
     & uc)-uu(i1,i2,i3,vc)*uy22r(i1,i2,i3,uc)-ux22r(i1,i2,i3,pc)
                     ut(i1,i2,i3,vc)= -uu(i1,i2,i3,uc)*ux22r(i1,i2,i3,
     & vc)-uu(i1,i2,i3,vc)*uy22r(i1,i2,i3,vc)-uy22r(i1,i2,i3,pc)
                 end if
                end do
                end do
                end do
               else
                do i3=n3a,n3b
                do i2=n2a,n2b
                do i1=n1a,n1b
                   ! explicit terms only, no diffusion
                    ut(i1,i2,i3,uc)= -uu(i1,i2,i3,uc)*ux22r(i1,i2,i3,
     & uc)-uu(i1,i2,i3,vc)*uy22r(i1,i2,i3,uc)-ux22r(i1,i2,i3,pc)
                    ut(i1,i2,i3,vc)= -uu(i1,i2,i3,uc)*ux22r(i1,i2,i3,
     & vc)-uu(i1,i2,i3,vc)*uy22r(i1,i2,i3,vc)-uy22r(i1,i2,i3,pc)
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
                      !   getUpwindAdvection(u,i1,i2,i3,NONE,2,2,rectangular, agu)
                      ! #Elif "UPWIND" == "BWENO" 
                      !   ! --- Bweno scheme ---
                      !   ! for testing output this message:
                      !   if( t.le. 0. )then
                      !      write(*,'(" getAdvection BWENO scheme (7)")') 
                      !   end if
                      !   getBwenoAdvection(u,i1,i2,i3,NONE,2,2,rectangular, agu)
                      ! #Else
                      !   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
                      !   stop 999
                      ! #End
                        ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-
     & ux22r(i1,i2,i3,pc)
                        ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-
     & uy22r(i1,i2,i3,pc)
                    ! include implicit terms - diffusion
                     uti(i1,i2,i3,uc)= nu*ulaplacian22r(i1,i2,i3,uc)
                     uti(i1,i2,i3,vc)= nu*ulaplacian22r(i1,i2,i3,vc)
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
                     !   getUpwindAdvection(u,i1,i2,i3,NONE,2,2,rectangular, agu)
                     ! #Elif "UPWIND" == "BWENO" 
                     !   ! --- Bweno scheme ---
                     !   ! for testing output this message:
                     !   if( t.le. 0. )then
                     !      write(*,'(" getAdvection BWENO scheme (7)")') 
                     !   end if
                     !   getBwenoAdvection(u,i1,i2,i3,NONE,2,2,rectangular, agu)
                     ! #Else
                     !   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
                     !   stop 999
                     ! #End
                       ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-ux22r(
     & i1,i2,i3,pc)
                       ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-uy22r(
     & i1,i2,i3,pc)
                   ! include implicit terms - diffusion
                    uti(i1,i2,i3,uc)= nu*ulaplacian22r(i1,i2,i3,uc)
                    uti(i1,i2,i3,vc)= nu*ulaplacian22r(i1,i2,i3,vc)
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
                      !   getUpwindAdvection(u,i1,i2,i3,NONE,2,2,rectangular, agu)
                      ! #Elif "UPWIND" == "BWENO" 
                      !   ! --- Bweno scheme ---
                      !   ! for testing output this message:
                      !   if( t.le. 0. )then
                      !      write(*,'(" getAdvection BWENO scheme (7)")') 
                      !   end if
                      !   getBwenoAdvection(u,i1,i2,i3,NONE,2,2,rectangular, agu)
                      ! #Else
                      !   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
                      !   stop 999
                      ! #End
                        ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-
     & ux22r(i1,i2,i3,pc)
                        ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-
     & uy22r(i1,i2,i3,pc)
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
                     !   getUpwindAdvection(u,i1,i2,i3,NONE,2,2,rectangular, agu)
                     ! #Elif "UPWIND" == "BWENO" 
                     !   ! --- Bweno scheme ---
                     !   ! for testing output this message:
                     !   if( t.le. 0. )then
                     !      write(*,'(" getAdvection BWENO scheme (7)")') 
                     !   end if
                     !   getBwenoAdvection(u,i1,i2,i3,NONE,2,2,rectangular, agu)
                     ! #Else
                     !   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
                     !   stop 999
                     ! #End
                       ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-ux22r(
     & i1,i2,i3,pc)
                       ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-uy22r(
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
                      !   getUpwindAdvection(u,i1,i2,i3,NONE,2,2,rectangular, agu)
                      ! #Elif "BWENO" == "BWENO" 
                      !   ! --- Bweno scheme ---
                      !   ! for testing output this message:
                      !   if( t.le. 0. )then
                      !      write(*,'(" getAdvection BWENO scheme (7)")') 
                      !   end if
                      !   getBwenoAdvection(u,i1,i2,i3,NONE,2,2,rectangular, agu)
                      ! #Else
                      !   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
                      !   stop 999
                      ! #End
                        ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-
     & ux22r(i1,i2,i3,pc)
                        ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-
     & uy22r(i1,i2,i3,pc)
                    ! include implicit terms - diffusion
                     uti(i1,i2,i3,uc)= nu*ulaplacian22r(i1,i2,i3,uc)
                     uti(i1,i2,i3,vc)= nu*ulaplacian22r(i1,i2,i3,vc)
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
                     !   getUpwindAdvection(u,i1,i2,i3,NONE,2,2,rectangular, agu)
                     ! #Elif "BWENO" == "BWENO" 
                     !   ! --- Bweno scheme ---
                     !   ! for testing output this message:
                     !   if( t.le. 0. )then
                     !      write(*,'(" getAdvection BWENO scheme (7)")') 
                     !   end if
                     !   getBwenoAdvection(u,i1,i2,i3,NONE,2,2,rectangular, agu)
                     ! #Else
                     !   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
                     !   stop 999
                     ! #End
                       ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-ux22r(
     & i1,i2,i3,pc)
                       ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-uy22r(
     & i1,i2,i3,pc)
                   ! include implicit terms - diffusion
                    uti(i1,i2,i3,uc)= nu*ulaplacian22r(i1,i2,i3,uc)
                    uti(i1,i2,i3,vc)= nu*ulaplacian22r(i1,i2,i3,vc)
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
                      !   getUpwindAdvection(u,i1,i2,i3,NONE,2,2,rectangular, agu)
                      ! #Elif "BWENO" == "BWENO" 
                      !   ! --- Bweno scheme ---
                      !   ! for testing output this message:
                      !   if( t.le. 0. )then
                      !      write(*,'(" getAdvection BWENO scheme (7)")') 
                      !   end if
                      !   getBwenoAdvection(u,i1,i2,i3,NONE,2,2,rectangular, agu)
                      ! #Else
                      !   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
                      !   stop 999
                      ! #End
                        ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-
     & ux22r(i1,i2,i3,pc)
                        ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-
     & uy22r(i1,i2,i3,pc)
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
                     !   getUpwindAdvection(u,i1,i2,i3,NONE,2,2,rectangular, agu)
                     ! #Elif "BWENO" == "BWENO" 
                     !   ! --- Bweno scheme ---
                     !   ! for testing output this message:
                     !   if( t.le. 0. )then
                     !      write(*,'(" getAdvection BWENO scheme (7)")') 
                     !   end if
                     !   getBwenoAdvection(u,i1,i2,i3,NONE,2,2,rectangular, agu)
                     ! #Else
                     !   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
                     !   stop 999
                     ! #End
                       ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-ux22r(
     & i1,i2,i3,pc)
                       ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-uy22r(
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
                  ut(i1,i2,i3,uc)= -uu(i1,i2,i3,uc)*ux22r(i1,i2,i3,uc)-
     & uu(i1,i2,i3,vc)*uy22r(i1,i2,i3,uc)-ux22r(i1,i2,i3,pc)+nu*
     & ulaplacian22r(i1,i2,i3,uc)
                  ut(i1,i2,i3,vc)= -uu(i1,i2,i3,uc)*ux22r(i1,i2,i3,vc)-
     & uu(i1,i2,i3,vc)*uy22r(i1,i2,i3,vc)-uy22r(i1,i2,i3,pc)+nu*
     & ulaplacian22r(i1,i2,i3,vc)
                  ! -- add on axisymmetric corrections ---
                  ri=radiusInverse(i1,i2,i3)
                  if( ri.ne.0. )then
                    ut(i1,i2,i3,uc)=ut(i1,i2,i3,uc)+nu*( uy22r(i1,i2,
     & i3,uc)*ri )
                    ut(i1,i2,i3,vc)=ut(i1,i2,i3,vc)+nu*( (uy22r(i1,i2,
     & i3,vc)-uu(i1,i2,i3,vc)*ri)*ri )
                  else
                    ut(i1,i2,i3,uc)=ut(i1,i2,i3,uc)+nu*( uyy22r(i1,i2,
     & i3,uc) )
                    ut(i1,i2,i3,vc)=ut(i1,i2,i3,vc)+nu*( .5*uyy22r(i1,
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
                 ut(i1,i2,i3,uc)= -uu(i1,i2,i3,uc)*ux22r(i1,i2,i3,uc)-
     & uu(i1,i2,i3,vc)*uy22r(i1,i2,i3,uc)-ux22r(i1,i2,i3,pc)+nu*
     & ulaplacian22r(i1,i2,i3,uc)
                 ut(i1,i2,i3,vc)= -uu(i1,i2,i3,uc)*ux22r(i1,i2,i3,vc)-
     & uu(i1,i2,i3,vc)*uy22r(i1,i2,i3,vc)-uy22r(i1,i2,i3,pc)+nu*
     & ulaplacian22r(i1,i2,i3,vc)
                 ! -- add on axisymmetric corrections ---
                 ri=radiusInverse(i1,i2,i3)
                 if( ri.ne.0. )then
                   ut(i1,i2,i3,uc)=ut(i1,i2,i3,uc)+nu*( uy22r(i1,i2,i3,
     & uc)*ri )
                   ut(i1,i2,i3,vc)=ut(i1,i2,i3,vc)+nu*( (uy22r(i1,i2,
     & i3,vc)-uu(i1,i2,i3,vc)*ri)*ri )
                 else
                   ut(i1,i2,i3,uc)=ut(i1,i2,i3,uc)+nu*( uyy22r(i1,i2,
     & i3,uc) )
                   ut(i1,i2,i3,vc)=ut(i1,i2,i3,vc)+nu*( .5*uyy22r(i1,
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
                    ut(i1,i2,i3,uc)= -uu(i1,i2,i3,uc)*ux22r(i1,i2,i3,
     & uc)-uu(i1,i2,i3,vc)*uy22r(i1,i2,i3,uc)-ux22r(i1,i2,i3,pc)
                    ut(i1,i2,i3,vc)= -uu(i1,i2,i3,uc)*ux22r(i1,i2,i3,
     & vc)-uu(i1,i2,i3,vc)*uy22r(i1,i2,i3,vc)-uy22r(i1,i2,i3,pc)
                   ! include implicit terms - diffusion
                    uti(i1,i2,i3,uc)= nu*ulaplacian22r(i1,i2,i3,uc)
                    uti(i1,i2,i3,vc)= nu*ulaplacian22r(i1,i2,i3,vc)
                    ri=radiusInverse(i1,i2,i3)
                    if( ri.ne.0. )then
                      uti(i1,i2,i3,uc)=uti(i1,i2,i3,uc)+nu*( uy22r(i1,
     & i2,i3,uc)*ri )
                      uti(i1,i2,i3,vc)=uti(i1,i2,i3,vc)+nu*( (uy22r(i1,
     & i2,i3,vc)-uu(i1,i2,i3,vc)*ri)*ri )
                    else
                      uti(i1,i2,i3,uc)=uti(i1,i2,i3,uc)+nu*( uyy22r(i1,
     & i2,i3,uc) )
                      uti(i1,i2,i3,vc)=uti(i1,i2,i3,vc)+nu*( .5*uyy22r(
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
                   ut(i1,i2,i3,uc)= -uu(i1,i2,i3,uc)*ux22r(i1,i2,i3,uc)
     & -uu(i1,i2,i3,vc)*uy22r(i1,i2,i3,uc)-ux22r(i1,i2,i3,pc)
                   ut(i1,i2,i3,vc)= -uu(i1,i2,i3,uc)*ux22r(i1,i2,i3,vc)
     & -uu(i1,i2,i3,vc)*uy22r(i1,i2,i3,vc)-uy22r(i1,i2,i3,pc)
                  ! include implicit terms - diffusion
                   uti(i1,i2,i3,uc)= nu*ulaplacian22r(i1,i2,i3,uc)
                   uti(i1,i2,i3,vc)= nu*ulaplacian22r(i1,i2,i3,vc)
                   ri=radiusInverse(i1,i2,i3)
                   if( ri.ne.0. )then
                     uti(i1,i2,i3,uc)=uti(i1,i2,i3,uc)+nu*( uy22r(i1,
     & i2,i3,uc)*ri )
                     uti(i1,i2,i3,vc)=uti(i1,i2,i3,vc)+nu*( (uy22r(i1,
     & i2,i3,vc)-uu(i1,i2,i3,vc)*ri)*ri )
                   else
                     uti(i1,i2,i3,uc)=uti(i1,i2,i3,uc)+nu*( uyy22r(i1,
     & i2,i3,uc) )
                     uti(i1,i2,i3,vc)=uti(i1,i2,i3,vc)+nu*( .5*uyy22r(
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
                    ut(i1,i2,i3,uc)= -uu(i1,i2,i3,uc)*ux22r(i1,i2,i3,
     & uc)-uu(i1,i2,i3,vc)*uy22r(i1,i2,i3,uc)-ux22r(i1,i2,i3,pc)
                    ut(i1,i2,i3,vc)= -uu(i1,i2,i3,uc)*ux22r(i1,i2,i3,
     & vc)-uu(i1,i2,i3,vc)*uy22r(i1,i2,i3,vc)-uy22r(i1,i2,i3,pc)
                end if
               end do
               end do
               end do
              else
               do i3=n3a,n3b
               do i2=n2a,n2b
               do i1=n1a,n1b
                  ! explicit terms only, no diffusion
                   ut(i1,i2,i3,uc)= -uu(i1,i2,i3,uc)*ux22r(i1,i2,i3,uc)
     & -uu(i1,i2,i3,vc)*uy22r(i1,i2,i3,uc)-ux22r(i1,i2,i3,pc)
                   ut(i1,i2,i3,vc)= -uu(i1,i2,i3,uc)*ux22r(i1,i2,i3,vc)
     & -uu(i1,i2,i3,vc)*uy22r(i1,i2,i3,vc)-uy22r(i1,i2,i3,pc)
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
                    ut(i1,i2,i3,uc)= -uu(i1,i2,i3,uc)*ux22(i1,i2,i3,uc)
     & -uu(i1,i2,i3,vc)*uy22(i1,i2,i3,uc)-ux22(i1,i2,i3,pc)+nu*
     & ulaplacian22(i1,i2,i3,uc)
                    ut(i1,i2,i3,vc)= -uu(i1,i2,i3,uc)*ux22(i1,i2,i3,vc)
     & -uu(i1,i2,i3,vc)*uy22(i1,i2,i3,vc)-uy22(i1,i2,i3,pc)+nu*
     & ulaplacian22(i1,i2,i3,vc)
                end if
               end do
               end do
               end do
              else
               do i3=n3a,n3b
               do i2=n2a,n2b
               do i1=n1a,n1b
                  ! INS, no AD
                   ut(i1,i2,i3,uc)= -uu(i1,i2,i3,uc)*ux22(i1,i2,i3,uc)-
     & uu(i1,i2,i3,vc)*uy22(i1,i2,i3,uc)-ux22(i1,i2,i3,pc)+nu*
     & ulaplacian22(i1,i2,i3,uc)
                   ut(i1,i2,i3,vc)= -uu(i1,i2,i3,uc)*ux22(i1,i2,i3,vc)-
     & uu(i1,i2,i3,vc)*uy22(i1,i2,i3,vc)-uy22(i1,i2,i3,pc)+nu*
     & ulaplacian22(i1,i2,i3,vc)
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
                    !   getUpwindAdvection(u,i1,i2,i3,NONE,2,2,curvilinear, agu)
                    ! #Elif "UPWIND" == "BWENO" 
                    !   ! --- Bweno scheme ---
                    !   ! for testing output this message:
                    !   if( t.le. 0. )then
                    !      write(*,'(" getAdvection BWENO scheme (7)")') 
                    !   end if
                    !   getBwenoAdvection(u,i1,i2,i3,NONE,2,2,curvilinear, agu)
                    ! #Else
                    !   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
                    !   stop 999
                    ! #End
                       ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-ux22(
     & i1,i2,i3,pc)+nu*ulaplacian22(i1,i2,i3,uc)
                       ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-uy22(
     & i1,i2,i3,pc)+nu*ulaplacian22(i1,i2,i3,vc)
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
                   !   getUpwindAdvection(u,i1,i2,i3,NONE,2,2,curvilinear, agu)
                   ! #Elif "UPWIND" == "BWENO" 
                   !   ! --- Bweno scheme ---
                   !   ! for testing output this message:
                   !   if( t.le. 0. )then
                   !      write(*,'(" getAdvection BWENO scheme (7)")') 
                   !   end if
                   !   getBwenoAdvection(u,i1,i2,i3,NONE,2,2,curvilinear, agu)
                   ! #Else
                   !   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
                   !   stop 999
                   ! #End
                      ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-ux22(
     & i1,i2,i3,pc)+nu*ulaplacian22(i1,i2,i3,uc)
                      ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-uy22(
     & i1,i2,i3,pc)+nu*ulaplacian22(i1,i2,i3,vc)
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
                    !   getUpwindAdvection(u,i1,i2,i3,NONE,2,2,curvilinear, agu)
                    ! #Elif "BWENO" == "BWENO" 
                    !   ! --- Bweno scheme ---
                    !   ! for testing output this message:
                    !   if( t.le. 0. )then
                    !      write(*,'(" getAdvection BWENO scheme (7)")') 
                    !   end if
                    !   getBwenoAdvection(u,i1,i2,i3,NONE,2,2,curvilinear, agu)
                    ! #Else
                    !   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
                    !   stop 999
                    ! #End
                       ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-ux22(
     & i1,i2,i3,pc)+nu*ulaplacian22(i1,i2,i3,uc)
                       ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-uy22(
     & i1,i2,i3,pc)+nu*ulaplacian22(i1,i2,i3,vc)
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
                   !   getUpwindAdvection(u,i1,i2,i3,NONE,2,2,curvilinear, agu)
                   ! #Elif "BWENO" == "BWENO" 
                   !   ! --- Bweno scheme ---
                   !   ! for testing output this message:
                   !   if( t.le. 0. )then
                   !      write(*,'(" getAdvection BWENO scheme (7)")') 
                   !   end if
                   !   getBwenoAdvection(u,i1,i2,i3,NONE,2,2,curvilinear, agu)
                   ! #Else
                   !   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
                   !   stop 999
                   ! #End
                      ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-ux22(
     & i1,i2,i3,pc)+nu*ulaplacian22(i1,i2,i3,uc)
                      ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-uy22(
     & i1,i2,i3,pc)+nu*ulaplacian22(i1,i2,i3,vc)
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
                     ut(i1,i2,i3,uc)= -uu(i1,i2,i3,uc)*ux22(i1,i2,i3,
     & uc)-uu(i1,i2,i3,vc)*uy22(i1,i2,i3,uc)-ux22(i1,i2,i3,pc)
                     ut(i1,i2,i3,vc)= -uu(i1,i2,i3,uc)*ux22(i1,i2,i3,
     & vc)-uu(i1,i2,i3,vc)*uy22(i1,i2,i3,vc)-uy22(i1,i2,i3,pc)
                    ! include implicit terms - diffusion
                     uti(i1,i2,i3,uc)= nu*ulaplacian22(i1,i2,i3,uc)
                     uti(i1,i2,i3,vc)= nu*ulaplacian22(i1,i2,i3,vc)
                 end if
                end do
                end do
                end do
               else
                do i3=n3a,n3b
                do i2=n2a,n2b
                do i1=n1a,n1b
                   ! explicit terms only, no diffusion
                    ut(i1,i2,i3,uc)= -uu(i1,i2,i3,uc)*ux22(i1,i2,i3,uc)
     & -uu(i1,i2,i3,vc)*uy22(i1,i2,i3,uc)-ux22(i1,i2,i3,pc)
                    ut(i1,i2,i3,vc)= -uu(i1,i2,i3,uc)*ux22(i1,i2,i3,vc)
     & -uu(i1,i2,i3,vc)*uy22(i1,i2,i3,vc)-uy22(i1,i2,i3,pc)
                   ! include implicit terms - diffusion
                    uti(i1,i2,i3,uc)= nu*ulaplacian22(i1,i2,i3,uc)
                    uti(i1,i2,i3,vc)= nu*ulaplacian22(i1,i2,i3,vc)
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
                     ut(i1,i2,i3,uc)= -uu(i1,i2,i3,uc)*ux22(i1,i2,i3,
     & uc)-uu(i1,i2,i3,vc)*uy22(i1,i2,i3,uc)-ux22(i1,i2,i3,pc)
                     ut(i1,i2,i3,vc)= -uu(i1,i2,i3,uc)*ux22(i1,i2,i3,
     & vc)-uu(i1,i2,i3,vc)*uy22(i1,i2,i3,vc)-uy22(i1,i2,i3,pc)
                 end if
                end do
                end do
                end do
               else
                do i3=n3a,n3b
                do i2=n2a,n2b
                do i1=n1a,n1b
                   ! explicit terms only, no diffusion
                    ut(i1,i2,i3,uc)= -uu(i1,i2,i3,uc)*ux22(i1,i2,i3,uc)
     & -uu(i1,i2,i3,vc)*uy22(i1,i2,i3,uc)-ux22(i1,i2,i3,pc)
                    ut(i1,i2,i3,vc)= -uu(i1,i2,i3,uc)*ux22(i1,i2,i3,vc)
     & -uu(i1,i2,i3,vc)*uy22(i1,i2,i3,vc)-uy22(i1,i2,i3,pc)
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
                      !   getUpwindAdvection(u,i1,i2,i3,NONE,2,2,curvilinear, agu)
                      ! #Elif "UPWIND" == "BWENO" 
                      !   ! --- Bweno scheme ---
                      !   ! for testing output this message:
                      !   if( t.le. 0. )then
                      !      write(*,'(" getAdvection BWENO scheme (7)")') 
                      !   end if
                      !   getBwenoAdvection(u,i1,i2,i3,NONE,2,2,curvilinear, agu)
                      ! #Else
                      !   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
                      !   stop 999
                      ! #End
                        ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-ux22(
     & i1,i2,i3,pc)
                        ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-uy22(
     & i1,i2,i3,pc)
                    ! include implicit terms - diffusion
                     uti(i1,i2,i3,uc)= nu*ulaplacian22(i1,i2,i3,uc)
                     uti(i1,i2,i3,vc)= nu*ulaplacian22(i1,i2,i3,vc)
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
                     !   getUpwindAdvection(u,i1,i2,i3,NONE,2,2,curvilinear, agu)
                     ! #Elif "UPWIND" == "BWENO" 
                     !   ! --- Bweno scheme ---
                     !   ! for testing output this message:
                     !   if( t.le. 0. )then
                     !      write(*,'(" getAdvection BWENO scheme (7)")') 
                     !   end if
                     !   getBwenoAdvection(u,i1,i2,i3,NONE,2,2,curvilinear, agu)
                     ! #Else
                     !   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
                     !   stop 999
                     ! #End
                       ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-ux22(
     & i1,i2,i3,pc)
                       ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-uy22(
     & i1,i2,i3,pc)
                   ! include implicit terms - diffusion
                    uti(i1,i2,i3,uc)= nu*ulaplacian22(i1,i2,i3,uc)
                    uti(i1,i2,i3,vc)= nu*ulaplacian22(i1,i2,i3,vc)
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
                      !   getUpwindAdvection(u,i1,i2,i3,NONE,2,2,curvilinear, agu)
                      ! #Elif "UPWIND" == "BWENO" 
                      !   ! --- Bweno scheme ---
                      !   ! for testing output this message:
                      !   if( t.le. 0. )then
                      !      write(*,'(" getAdvection BWENO scheme (7)")') 
                      !   end if
                      !   getBwenoAdvection(u,i1,i2,i3,NONE,2,2,curvilinear, agu)
                      ! #Else
                      !   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
                      !   stop 999
                      ! #End
                        ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-ux22(
     & i1,i2,i3,pc)
                        ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-uy22(
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
                     !   getUpwindAdvection(u,i1,i2,i3,NONE,2,2,curvilinear, agu)
                     ! #Elif "UPWIND" == "BWENO" 
                     !   ! --- Bweno scheme ---
                     !   ! for testing output this message:
                     !   if( t.le. 0. )then
                     !      write(*,'(" getAdvection BWENO scheme (7)")') 
                     !   end if
                     !   getBwenoAdvection(u,i1,i2,i3,NONE,2,2,curvilinear, agu)
                     ! #Else
                     !   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
                     !   stop 999
                     ! #End
                       ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-ux22(
     & i1,i2,i3,pc)
                       ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-uy22(
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
                      !   getUpwindAdvection(u,i1,i2,i3,NONE,2,2,curvilinear, agu)
                      ! #Elif "BWENO" == "BWENO" 
                      !   ! --- Bweno scheme ---
                      !   ! for testing output this message:
                      !   if( t.le. 0. )then
                      !      write(*,'(" getAdvection BWENO scheme (7)")') 
                      !   end if
                      !   getBwenoAdvection(u,i1,i2,i3,NONE,2,2,curvilinear, agu)
                      ! #Else
                      !   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
                      !   stop 999
                      ! #End
                        ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-ux22(
     & i1,i2,i3,pc)
                        ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-uy22(
     & i1,i2,i3,pc)
                    ! include implicit terms - diffusion
                     uti(i1,i2,i3,uc)= nu*ulaplacian22(i1,i2,i3,uc)
                     uti(i1,i2,i3,vc)= nu*ulaplacian22(i1,i2,i3,vc)
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
                     !   getUpwindAdvection(u,i1,i2,i3,NONE,2,2,curvilinear, agu)
                     ! #Elif "BWENO" == "BWENO" 
                     !   ! --- Bweno scheme ---
                     !   ! for testing output this message:
                     !   if( t.le. 0. )then
                     !      write(*,'(" getAdvection BWENO scheme (7)")') 
                     !   end if
                     !   getBwenoAdvection(u,i1,i2,i3,NONE,2,2,curvilinear, agu)
                     ! #Else
                     !   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
                     !   stop 999
                     ! #End
                       ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-ux22(
     & i1,i2,i3,pc)
                       ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-uy22(
     & i1,i2,i3,pc)
                   ! include implicit terms - diffusion
                    uti(i1,i2,i3,uc)= nu*ulaplacian22(i1,i2,i3,uc)
                    uti(i1,i2,i3,vc)= nu*ulaplacian22(i1,i2,i3,vc)
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
                      !   getUpwindAdvection(u,i1,i2,i3,NONE,2,2,curvilinear, agu)
                      ! #Elif "BWENO" == "BWENO" 
                      !   ! --- Bweno scheme ---
                      !   ! for testing output this message:
                      !   if( t.le. 0. )then
                      !      write(*,'(" getAdvection BWENO scheme (7)")') 
                      !   end if
                      !   getBwenoAdvection(u,i1,i2,i3,NONE,2,2,curvilinear, agu)
                      ! #Else
                      !   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
                      !   stop 999
                      ! #End
                        ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-ux22(
     & i1,i2,i3,pc)
                        ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-uy22(
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
                     !   getUpwindAdvection(u,i1,i2,i3,NONE,2,2,curvilinear, agu)
                     ! #Elif "BWENO" == "BWENO" 
                     !   ! --- Bweno scheme ---
                     !   ! for testing output this message:
                     !   if( t.le. 0. )then
                     !      write(*,'(" getAdvection BWENO scheme (7)")') 
                     !   end if
                     !   getBwenoAdvection(u,i1,i2,i3,NONE,2,2,curvilinear, agu)
                     ! #Else
                     !   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
                     !   stop 999
                     ! #End
                       ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-ux22(
     & i1,i2,i3,pc)
                       ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-uy22(
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
                  ut(i1,i2,i3,uc)= -uu(i1,i2,i3,uc)*ux22(i1,i2,i3,uc)-
     & uu(i1,i2,i3,vc)*uy22(i1,i2,i3,uc)-ux22(i1,i2,i3,pc)+nu*
     & ulaplacian22(i1,i2,i3,uc)
                  ut(i1,i2,i3,vc)= -uu(i1,i2,i3,uc)*ux22(i1,i2,i3,vc)-
     & uu(i1,i2,i3,vc)*uy22(i1,i2,i3,vc)-uy22(i1,i2,i3,pc)+nu*
     & ulaplacian22(i1,i2,i3,vc)
                  ! -- add on axisymmetric corrections ---
                  ri=radiusInverse(i1,i2,i3)
                  if( ri.ne.0. )then
                    ut(i1,i2,i3,uc)=ut(i1,i2,i3,uc)+nu*( uy22(i1,i2,i3,
     & uc)*ri )
                    ut(i1,i2,i3,vc)=ut(i1,i2,i3,vc)+nu*( (uy22(i1,i2,
     & i3,vc)-uu(i1,i2,i3,vc)*ri)*ri )
                  else
                    ut(i1,i2,i3,uc)=ut(i1,i2,i3,uc)+nu*( uyy22(i1,i2,
     & i3,uc) )
                    ut(i1,i2,i3,vc)=ut(i1,i2,i3,vc)+nu*( .5*uyy22(i1,
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
                 ut(i1,i2,i3,uc)= -uu(i1,i2,i3,uc)*ux22(i1,i2,i3,uc)-
     & uu(i1,i2,i3,vc)*uy22(i1,i2,i3,uc)-ux22(i1,i2,i3,pc)+nu*
     & ulaplacian22(i1,i2,i3,uc)
                 ut(i1,i2,i3,vc)= -uu(i1,i2,i3,uc)*ux22(i1,i2,i3,vc)-
     & uu(i1,i2,i3,vc)*uy22(i1,i2,i3,vc)-uy22(i1,i2,i3,pc)+nu*
     & ulaplacian22(i1,i2,i3,vc)
                 ! -- add on axisymmetric corrections ---
                 ri=radiusInverse(i1,i2,i3)
                 if( ri.ne.0. )then
                   ut(i1,i2,i3,uc)=ut(i1,i2,i3,uc)+nu*( uy22(i1,i2,i3,
     & uc)*ri )
                   ut(i1,i2,i3,vc)=ut(i1,i2,i3,vc)+nu*( (uy22(i1,i2,i3,
     & vc)-uu(i1,i2,i3,vc)*ri)*ri )
                 else
                   ut(i1,i2,i3,uc)=ut(i1,i2,i3,uc)+nu*( uyy22(i1,i2,i3,
     & uc) )
                   ut(i1,i2,i3,vc)=ut(i1,i2,i3,vc)+nu*( .5*uyy22(i1,i2,
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
                    ut(i1,i2,i3,uc)= -uu(i1,i2,i3,uc)*ux22(i1,i2,i3,uc)
     & -uu(i1,i2,i3,vc)*uy22(i1,i2,i3,uc)-ux22(i1,i2,i3,pc)
                    ut(i1,i2,i3,vc)= -uu(i1,i2,i3,uc)*ux22(i1,i2,i3,vc)
     & -uu(i1,i2,i3,vc)*uy22(i1,i2,i3,vc)-uy22(i1,i2,i3,pc)
                   ! include implicit terms - diffusion
                    uti(i1,i2,i3,uc)= nu*ulaplacian22(i1,i2,i3,uc)
                    uti(i1,i2,i3,vc)= nu*ulaplacian22(i1,i2,i3,vc)
                    ri=radiusInverse(i1,i2,i3)
                    if( ri.ne.0. )then
                      uti(i1,i2,i3,uc)=uti(i1,i2,i3,uc)+nu*( uy22(i1,
     & i2,i3,uc)*ri )
                      uti(i1,i2,i3,vc)=uti(i1,i2,i3,vc)+nu*( (uy22(i1,
     & i2,i3,vc)-uu(i1,i2,i3,vc)*ri)*ri )
                    else
                      uti(i1,i2,i3,uc)=uti(i1,i2,i3,uc)+nu*( uyy22(i1,
     & i2,i3,uc) )
                      uti(i1,i2,i3,vc)=uti(i1,i2,i3,vc)+nu*( .5*uyy22(
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
                   ut(i1,i2,i3,uc)= -uu(i1,i2,i3,uc)*ux22(i1,i2,i3,uc)-
     & uu(i1,i2,i3,vc)*uy22(i1,i2,i3,uc)-ux22(i1,i2,i3,pc)
                   ut(i1,i2,i3,vc)= -uu(i1,i2,i3,uc)*ux22(i1,i2,i3,vc)-
     & uu(i1,i2,i3,vc)*uy22(i1,i2,i3,vc)-uy22(i1,i2,i3,pc)
                  ! include implicit terms - diffusion
                   uti(i1,i2,i3,uc)= nu*ulaplacian22(i1,i2,i3,uc)
                   uti(i1,i2,i3,vc)= nu*ulaplacian22(i1,i2,i3,vc)
                   ri=radiusInverse(i1,i2,i3)
                   if( ri.ne.0. )then
                     uti(i1,i2,i3,uc)=uti(i1,i2,i3,uc)+nu*( uy22(i1,i2,
     & i3,uc)*ri )
                     uti(i1,i2,i3,vc)=uti(i1,i2,i3,vc)+nu*( (uy22(i1,
     & i2,i3,vc)-uu(i1,i2,i3,vc)*ri)*ri )
                   else
                     uti(i1,i2,i3,uc)=uti(i1,i2,i3,uc)+nu*( uyy22(i1,
     & i2,i3,uc) )
                     uti(i1,i2,i3,vc)=uti(i1,i2,i3,vc)+nu*( .5*uyy22(
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
                    ut(i1,i2,i3,uc)= -uu(i1,i2,i3,uc)*ux22(i1,i2,i3,uc)
     & -uu(i1,i2,i3,vc)*uy22(i1,i2,i3,uc)-ux22(i1,i2,i3,pc)
                    ut(i1,i2,i3,vc)= -uu(i1,i2,i3,uc)*ux22(i1,i2,i3,vc)
     & -uu(i1,i2,i3,vc)*uy22(i1,i2,i3,vc)-uy22(i1,i2,i3,pc)
                end if
               end do
               end do
               end do
              else
               do i3=n3a,n3b
               do i2=n2a,n2b
               do i1=n1a,n1b
                  ! explicit terms only, no diffusion
                   ut(i1,i2,i3,uc)= -uu(i1,i2,i3,uc)*ux22(i1,i2,i3,uc)-
     & uu(i1,i2,i3,vc)*uy22(i1,i2,i3,uc)-ux22(i1,i2,i3,pc)
                   ut(i1,i2,i3,vc)= -uu(i1,i2,i3,uc)*ux22(i1,i2,i3,vc)-
     & uu(i1,i2,i3,vc)*uy22(i1,i2,i3,vc)-uy22(i1,i2,i3,pc)
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
