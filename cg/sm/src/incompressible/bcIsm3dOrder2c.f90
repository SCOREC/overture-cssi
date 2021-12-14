! This file automatically generated from bcOptIsm.bf90 with bpp.
        subroutine bcIsm3dOrder2c( nd,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,gridIndexRange,u,mask,rsxy,xy,ndMatProp,matIndex,matValpc,matVal,boundaryCondition,addBoundaryForcing,interfaceType,dim,bcf00,bcf10,bcf01,bcf11,bcf02,bcf12,bcf0,bcOffset,ipar,rpar,ierr )
       ! ===================================================================================
       !  Boundary conditions for solid mechanics
       !
       !  gridType : 0=rectangular, 1=curvilinear
       !
       !  c2= mu/rho, c1=(mu+lambda)/rho;
       ! 
       ! The forcing for the boundary conditions can be accessed in two ways. One can either 
       ! use the arrays: 
       !       bcf00(i1,i2,i3,m), bcf10(i1,i2,i3,m), bcf01(i1,i2,i3,m), bcf11(i1,i2,i3,m), 
       !       bcf02(i1,i2,i3,m), bcf12(i1,i2,i3,m)
       ! which provide values for the 6 different faces in 6 different arrays. One can also
       ! access the same values using the single statement function
       !         bcf(side,axis,i1,i2,i3,m)
       ! which is defined below. 
       ! ===================================================================================
         implicit none
           integer nd, nd1a,nd1b,nd2a,nd2b,nd3a,nd3b, ierr
           real u(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
           integer mask(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
           real rsxy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1,0:nd-1)
           real xy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1)
           integer gridIndexRange(0:1,0:2),boundaryCondition(0:1,0:2)
           integer addBoundaryForcing(0:1,0:2)
           integer interfaceType(0:1,0:2,0:*)
           integer dim(0:1,0:2,0:1,0:2)
           real bcf00(dim(0,0,0,0):dim(1,0,0,0), dim(0,1,0,0):dim(1,1,0,0), dim(0,2,0,0):dim(1,2,0,0),0:*)
           real bcf10(dim(0,0,1,0):dim(1,0,1,0), dim(0,1,1,0):dim(1,1,1,0), dim(0,2,1,0):dim(1,2,1,0),0:*)
           real bcf01(dim(0,0,0,1):dim(1,0,0,1), dim(0,1,0,1):dim(1,1,0,1), dim(0,2,0,1):dim(1,2,0,1),0:*)
           real bcf11(dim(0,0,1,1):dim(1,0,1,1), dim(0,1,1,1):dim(1,1,1,1), dim(0,2,1,1):dim(1,2,1,1),0:*)
           real bcf02(dim(0,0,0,2):dim(1,0,0,2), dim(0,1,0,2):dim(1,1,0,2), dim(0,2,0,2):dim(1,2,0,2),0:*)
           real bcf12(dim(0,0,1,2):dim(1,0,1,2), dim(0,1,1,2):dim(1,1,1,2), dim(0,2,1,2):dim(1,2,1,2),0:*)
           real bcf0(0:*)
           integer*8 bcOffset(0:1,0:2)
           integer ipar(0:*)
           real rpar(0:*)
           ! -- Declare arrays for variable material properties --
           include '../declareVarMatProp.h'
           integer rectangular,curvilinear
           parameter(rectangular=0,curvilinear=1)
             ! work space arrays that must be saved from call to call:
       !**      real aa2(0:1,0:1,0:1,0:*),aa4(0:3,0:3,0:1,0:*),aa8(0:7,0:7,0:1,0:*)
       !**      integer ipvt2(0:1,0:*), ipvt4(0:3,0:*), ipvt8(0:7,0:*)
       !     --- local variables ----
         integer side,axis,grid,gridType,orderOfAccuracy,orderOfExtrapolation,twilightZone,assignTwilightZone,u1c,u2c,u3c,useWhereMask,debug,nn,n1,n2,pc,upwindSOS,useCurlCurlBoundaryCondition
         integer displacementTractionCorner
         real dx(0:2),dr(0:2)
         real t,ep,dt,c1,c2,rho,mu,lambda,alpha,beta
         integer axisp1,axisp2,i1,i2,i3,is1,is2,is3,j1,j2,j3,js1,js2,js3,ks1,ks2,ks3,is,js,it,nit
         integer option,initialized
         integer ghost,numGhost,numberOfGhostPoints,ghost1,ghost2
         integer side1,side2,side3, sidea,sideb,bc1,bc2, edgeDirection
         integer n1a,n1b,n2a,n2b,n3a,n3b
         integer nn1a,nn1b,nn2a,nn2b,nn3a,nn3b
         integer extra1a,extra1b,extra2a,extra2b,extra3a,extra3b
         logical checkResiduals
         integer bcSide,bcAdjacent
         logical adjustEnds
         integer leftRight,dir,axisp,gid(0:1,0:2)
         real urv(0:5),usv(0:5),utv(0:5)
         integer iw1,iw2,iw3
         real delta 
         integer hw1,hw2,hw3   
         real a11,a12,a21,a22,det,b0,b1,b2
         real a0,a1,cc0,cc1,d0,d1,dr0,ds0
         real aNormSq,divu,uAve
         real epsRatio,an1,an2,an3,aNorm,aNormi,ua,ub,nDotU,t1,t2,t3,tn1,tn2,tn3
         real epsx,tmp
         real tau11,tau12,tau13,tau21,tau22,tau23,tau31,tau32,tau33
         real ux,uy,uz,vx,vy,vz,wx,wy,wz
         real ux0,uy0,uz0,vx0,vy0,vz0,wx0,wy0,wz0
         real uxx0,uxy0,uxz0,uyy0,uyz0,uzz0
         real vxx0,vxy0,vxz0,vyy0,vyz0,vzz0
         real wxx0,wxy0,wxz0,wyy0,wyz0,wzz0
         real u0,v0,w0, u1,v1,w1
         real um,up,vm,vp,wm,wp
         real tau1,tau2,tau3,clap1,clap2,ulap1,vlap1,wlap1,ulap2,vlap2,wlap2,an1Cartesian,an2Cartesian
         real u1x,u1y,u1z, u2x,u2y,u2z, u3x,u3y,u3z
         real trac1,trac2,trac3
         integer numberOfEquations,job
         real a2(0:1,0:1),a3(0:2,0:2),a4(0:3,0:3),a6(0:5,0:5),a8(0:7,0:7),q(0:11),f(0:11),f0(0:11),rcond,work(0:11)
         integer ipvt(0:11)
         real resMax,fe(0:11),g(0:11)
         real err
         real Au(1:2,1:2), Av(1:2,1:2),cu11,cu12,cu21,cu22,cv11,cv12,cv21,cv22,f1e,f2e,f1,f2
         real res1,res2,resTol
         real u1e, u1ex, u1ey, u1exx, u1exy, u1eyy, u1exxx, u1exxy, u1exyy, u1eyyy
         real u2e, u2ex, u2ey, u2exx, u2exy, u2eyy, u2exxx, u2exxy, u2exyy, u2eyyy
         real pex,pey,pexx,pexy,peyy
         ! shift and components used in discrete delta 
         integer shift1(0:11), shift2(0:11), shift3(0:11), comp(0:11)
         integer startGhost
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
         !     --- start statement function ----
         real bcf
         integer kd,m,n
         real uxOneSided
         real dr12
         real dr22
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
         real urrs2
         real urss2
         real urrt2
         real usst2
         real urtt2
         real ustt2
         real urrrr2
         real ussss2
         real utttt2
         real urrss2
         real urrtt2
         real usstt2
         real urrrs2
         real ursss2
         real urrrt2
         real ussst2
         real urttt2
         real usttt2
         real rsxyr2
         real rsxys2
         real rsxyt2
         real rsxyrr2
         real rsxyss2
         real rsxyrs2
         real rsxytt2
         real rsxyrt2
         real rsxyst2
         real rsxyrrr2
         real rsxysss2
         real rsxyttt2
         real rsxyrrs2
         real rsxyrss2
         real rsxyrrt2
         real rsxysst2
         real rsxyrtt2
         real rsxystt2
         real rsxyrrrr2
         real rsxyssss2
         real rsxytttt2
         real rsxyrrss2
         real rsxyrrtt2
         real rsxysstt2
         real ux21
         real uy21
         real uz21
         real ux22
         real uy22
         real uz22
         real ux23
         real uy23
         real uz23
         real rsxyx21
         real rsxyx22
         real rsxyy22
         real rsxyx23
         real rsxyy23
         real rsxyz23
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
         real rsxyxx22
         real rsxyyy22
         real rsxyxy22
         real rsxyxxx22
         real rsxyxxy22
         real rsxyxyy22
         real rsxyyyy22
         real uxxx22
         real uxxy22
         real uxyy22
         real uyyy22
         real uxxxx22
         real uxxxy22
         real uxxyy22
         real uxyyy22
         real uyyyy22
         real uLapSq22
         real uxx23
         real uyy23
         real uzz23
         real uxy23
         real uxz23
         real uyz23
         real ulaplacian23
         real dx12
         real dx22
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
         real uLapSq22r
         real uxxx23r
         real uyyy23r
         real uzzz23r
         real uxxy23r
         real uxyy23r
         real uxxz23r
         real uyyz23r
         real uxzz23r
         real uyzz23r
         real uxxxx23r
         real uyyyy23r
         real uzzzz23r
         real uxxyy23r
         real uxxzz23r
         real uyyzz23r
         real dr14
         real dr24
         real ur4
         real us4
         real ut4
         real urr4
         real uss4
         real utt4
         real urs4
         real urt4
         real ust4
         real rsxyr4
         real rsxys4
         real rsxyt4
         real ux41
         real uy41
         real uz41
         real ux42
         real uy42
         real uz42
         real ux43
         real uy43
         real uz43
         real rsxyx41
         real rsxyx42
         real rsxyy42
         real rsxyx43
         real rsxyy43
         real rsxyz43
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
         real dx41
         real dx42
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
         !     The next macro call will define the difference approximation statement functions
         dr12(kd) = 1./(2.*dr(kd))
         dr22(kd) = 1./(dr(kd)**2)
         ur2(i1,i2,i3,kd)=(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))*dr12(0)
         us2(i1,i2,i3,kd)=(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))*dr12(1)
         ut2(i1,i2,i3,kd)=(u(i1,i2,i3+1,kd)-u(i1,i2,i3-1,kd))*dr12(2)
         urr2(i1,i2,i3,kd)=(-2.*u(i1,i2,i3,kd)+(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd)) )*dr22(0)
         uss2(i1,i2,i3,kd)=(-2.*u(i1,i2,i3,kd)+(u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd)) )*dr22(1)
         urs2(i1,i2,i3,kd)=(ur2(i1,i2+1,i3,kd)-ur2(i1,i2-1,i3,kd))*dr12(1)
         utt2(i1,i2,i3,kd)=(-2.*u(i1,i2,i3,kd)+(u(i1,i2,i3+1,kd)+u(i1,i2,i3-1,kd)) )*dr22(2)
         urt2(i1,i2,i3,kd)=(ur2(i1,i2,i3+1,kd)-ur2(i1,i2,i3-1,kd))*dr12(2)
         ust2(i1,i2,i3,kd)=(us2(i1,i2,i3+1,kd)-us2(i1,i2,i3-1,kd))*dr12(2)
         urrr2(i1,i2,i3,kd)=(-2.*(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))+(u(i1+2,i2,i3,kd)-u(i1-2,i2,i3,kd)) )*dr22(0)*dr12(0)
         usss2(i1,i2,i3,kd)=(-2.*(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))+(u(i1,i2+2,i3,kd)-u(i1,i2-2,i3,kd)) )*dr22(1)*dr12(1)
         uttt2(i1,i2,i3,kd)=(-2.*(u(i1,i2,i3+1,kd)-u(i1,i2,i3-1,kd))+(u(i1,i2,i3+2,kd)-u(i1,i2,i3-2,kd)) )*dr22(1)*dr12(2)
         urrs2(i1,i2,i3,kd)=( urr2(i1,i2+1,i3,kd)-urr2(i1,i2-1,i3,kd))/(2.*dr(1))
         urss2(i1,i2,i3,kd)=( uss2(i1+1,i2,i3,kd)-uss2(i1-1,i2,i3,kd))/(2.*dr(0))
         urrt2(i1,i2,i3,kd)=( urr2(i1,i2,i3+1,kd)-urr2(i1,i2,i3-1,kd))/(2.*dr(2))
         usst2(i1,i2,i3,kd)=( uss2(i1,i2,i3+1,kd)-uss2(i1,i2,i3-1,kd))/(2.*dr(2))
         urtt2(i1,i2,i3,kd)=( utt2(i1+1,i2,i3,kd)-utt2(i1-1,i2,i3,kd))/(2.*dr(0))
         ustt2(i1,i2,i3,kd)=( utt2(i1,i2+1,i3,kd)-utt2(i1,i2-1,i3,kd))/(2.*dr(1))
         urrrr2(i1,i2,i3,kd)=(6.*u(i1,i2,i3,kd)-4.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd))+(u(i1+2,i2,i3,kd)+u(i1-2,i2,i3,kd)) )/(dr(0)**4)
         ussss2(i1,i2,i3,kd)=(6.*u(i1,i2,i3,kd)-4.*(u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd))+(u(i1,i2+2,i3,kd)+u(i1,i2-2,i3,kd)) )/(dr(1)**4)
         utttt2(i1,i2,i3,kd)=(6.*u(i1,i2,i3,kd)-4.*(u(i1,i2,i3+1,kd)+u(i1,i2,i3-1,kd))+(u(i1,i2,i3+2,kd)+u(i1,i2,i3-2,kd)) )/(dr(2)**4)
         urrss2(i1,i2,i3,kd)=( 4.*u(i1,i2,i3,kd)-2.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd)+u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd))+   (u(i1+1,i2+1,i3,kd)+u(i1-1,i2+1,i3,kd)+u(i1+1,i2-1,i3,kd)+u(i1-1,i2-1,i3,kd)) )/(dr(0)**2*dr(1)**2)
         urrtt2(i1,i2,i3,kd)=( 4.*u(i1,i2,i3,kd)-2.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd)+u(i1,i2,i3+1,kd)+u(i1,i2,i3-1,kd))+   (u(i1+1,i2,i3+1,kd)+u(i1-1,i2,i3+1,kd)+u(i1+1,i2,i3-1,kd)+u(i1-1,i2,i3-1,kd)) )/(dr(0)**2*dr(2)**2)
         usstt2(i1,i2,i3,kd)=( 4.*u(i1,i2,i3,kd)-2.*(u(i1,i2+1,i3,kd)  +u(i1,i2-1,i3,kd)+  u(i1,i2  ,i3+1,kd)+u(i1,i2  ,i3-1,kd))+   (u(i1,i2+1,i3+1,kd)+u(i1,i2-1,i3+1,kd)+u(i1,i2+1,i3-1,kd)+u(i1,i2-1,i3-1,kd)) )/(dr(1)**2*dr(2)**2)
         urrrs2(i1,i2,i3,kd)=( urrr2(i1,i2+1,i3,kd)-urrr2(i1,i2-1,i3,kd))/(2.*dr(1))
         ursss2(i1,i2,i3,kd)=( usss2(i1+1,i2,i3,kd)-usss2(i1-1,i2,i3,kd))/(2.*dr(0))
         urrrt2(i1,i2,i3,kd)=( urrr2(i1,i2,i3+1,kd)-urrr2(i1,i2,i3-1,kd))/(2.*dr(2))
         ussst2(i1,i2,i3,kd)=( usss2(i1,i2,i3+1,kd)-usss2(i1,i2,i3-1,kd))/(2.*dr(2))
         urttt2(i1,i2,i3,kd)=( uttt2(i1+1,i2,i3,kd)-uttt2(i1-1,i2,i3,kd))/(2.*dr(0))
         usttt2(i1,i2,i3,kd)=( uttt2(i1,i2+1,i3,kd)-uttt2(i1,i2-1,i3,kd))/(2.*dr(1))
         rsxyr2(i1,i2,i3,m,n)=(rsxy(i1+1,i2,i3,m,n)-rsxy(i1-1,i2,i3,m,n))*dr12(0)
         rsxys2(i1,i2,i3,m,n)=(rsxy(i1,i2+1,i3,m,n)-rsxy(i1,i2-1,i3,m,n))*dr12(1)
         rsxyt2(i1,i2,i3,m,n)=(rsxy(i1,i2,i3+1,m,n)-rsxy(i1,i2,i3-1,m,n))*dr12(2)
         rsxyrr2(i1,i2,i3,m,n)=(-2.*rsxy(i1,i2,i3,m,n)+(rsxy(i1+1,i2,i3,m,n)+rsxy(i1-1,i2,i3,m,n)) )*dr22(0)
         rsxyss2(i1,i2,i3,m,n)=(-2.*rsxy(i1,i2,i3,m,n)+(rsxy(i1,i2+1,i3,m,n)+rsxy(i1,i2-1,i3,m,n)) )*dr22(1)
         rsxyrs2(i1,i2,i3,m,n)=(rsxyr2(i1,i2+1,i3,m,n)-rsxyr2(i1,i2-1,i3,m,n))*dr12(1)
         rsxytt2(i1,i2,i3,m,n)=(-2.*rsxy(i1,i2,i3,m,n)+(rsxy(i1,i2,i3+1,m,n)+rsxy(i1,i2,i3-1,m,n)) )*dr22(2)
         rsxyrt2(i1,i2,i3,m,n)=(rsxyr2(i1,i2,i3+1,m,n)-rsxyr2(i1,i2,i3-1,m,n))*dr12(2)
         rsxyst2(i1,i2,i3,m,n)=(rsxys2(i1,i2,i3+1,m,n)-rsxys2(i1,i2,i3-1,m,n))*dr12(2)
         rsxyrrr2(i1,i2,i3,m,n)=(-2.*(rsxy(i1+1,i2,i3,m,n)-rsxy(i1-1,i2,i3,m,n))+(rsxy(i1+2,i2,i3,m,n)-rsxy(i1-2,i2,i3,m,n)) )*dr22(0)*dr12(0)
         rsxysss2(i1,i2,i3,m,n)=(-2.*(rsxy(i1,i2+1,i3,m,n)-rsxy(i1,i2-1,i3,m,n))+(rsxy(i1,i2+2,i3,m,n)-rsxy(i1,i2-2,i3,m,n)) )*dr22(1)*dr12(1)
         rsxyttt2(i1,i2,i3,m,n)=(-2.*(rsxy(i1,i2,i3+1,m,n)-rsxy(i1,i2,i3-1,m,n))+(rsxy(i1,i2,i3+2,m,n)-rsxy(i1,i2,i3-2,m,n)) )*dr22(1)*dr12(2)
         rsxyrrs2(i1,i2,i3,m,n)=( rsxyrr2(i1,i2+1,i3,m,n)-rsxyrr2(i1,i2-1,i3,m,n))/(2.*dr(1))
         rsxyrss2(i1,i2,i3,m,n)=( rsxyss2(i1+1,i2,i3,m,n)-rsxyss2(i1-1,i2,i3,m,n))/(2.*dr(0))
         rsxyrrt2(i1,i2,i3,m,n)=( rsxyrr2(i1,i2,i3+1,m,n)-rsxyrr2(i1,i2,i3-1,m,n))/(2.*dr(2))
         rsxysst2(i1,i2,i3,m,n)=( rsxyss2(i1,i2,i3+1,m,n)-rsxyss2(i1,i2,i3-1,m,n))/(2.*dr(2))
         rsxyrtt2(i1,i2,i3,m,n)=( rsxytt2(i1+1,i2,i3,m,n)-rsxytt2(i1-1,i2,i3,m,n))/(2.*dr(0))
         rsxystt2(i1,i2,i3,m,n)=( rsxytt2(i1,i2+1,i3,m,n)-rsxytt2(i1,i2-1,i3,m,n))/(2.*dr(1))
         rsxyrrrr2(i1,i2,i3,m,n)=(6.*rsxy(i1,i2,i3,m,n)-4.*(rsxy(i1+1,i2,i3,m,n)+rsxy(i1-1,i2,i3,m,n))+(rsxy(i1+2,i2,i3,m,n)+rsxy(i1-2,i2,i3,m,n)) )/(dr(0)**4)
         rsxyssss2(i1,i2,i3,m,n)=(6.*rsxy(i1,i2,i3,m,n)-4.*(rsxy(i1,i2+1,i3,m,n)+rsxy(i1,i2-1,i3,m,n))+(rsxy(i1,i2+2,i3,m,n)+rsxy(i1,i2-2,i3,m,n)) )/(dr(1)**4)
         rsxytttt2(i1,i2,i3,m,n)=(6.*rsxy(i1,i2,i3,m,n)-4.*(rsxy(i1,i2,i3+1,m,n)+rsxy(i1,i2,i3-1,m,n))+(rsxy(i1,i2,i3+2,m,n)+rsxy(i1,i2,i3-2,m,n)) )/(dr(2)**4)
         rsxyrrss2(i1,i2,i3,m,n)=( 4.*rsxy(i1,i2,i3,m,n)-2.*(rsxy(i1+1,i2,i3,m,n)+rsxy(i1-1,i2,i3,m,n)+rsxy(i1,i2+1,i3,m,n)+rsxy(i1,i2-1,i3,m,n))+   (rsxy(i1+1,i2+1,i3,m,n)+rsxy(i1-1,i2+1,i3,m,n)+rsxy(i1+1,i2-1,i3,m,n)+rsxy(i1-1,i2-1,i3,m,n)) )/(dr(0)**2*dr(1)**2)
         rsxyrrtt2(i1,i2,i3,m,n)=( 4.*rsxy(i1,i2,i3,m,n)-2.*(rsxy(i1+1,i2,i3,m,n)+rsxy(i1-1,i2,i3,m,n)+rsxy(i1,i2,i3+1,m,n)+rsxy(i1,i2,i3-1,m,n))+   (rsxy(i1+1,i2,i3+1,m,n)+rsxy(i1-1,i2,i3+1,m,n)+rsxy(i1+1,i2,i3-1,m,n)+rsxy(i1-1,i2,i3-1,m,n)) )/(dr(0)**2*dr(2)**2)
         rsxysstt2(i1,i2,i3,m,n)=( 4.*rsxy(i1,i2,i3,m,n)-2.*(rsxy(i1,i2+1,i3,m,n)  +rsxy(i1,i2-1,i3,m,n)+  rsxy(i1,i2  ,i3+1,m,n)+rsxy(i1,i2  ,i3-1,m,n))+   (rsxy(i1,i2+1,i3+1,m,n)+rsxy(i1,i2-1,i3+1,m,n)+rsxy(i1,i2+1,i3-1,m,n)+rsxy(i1,i2-1,i3-1,m,n)) )/(dr(1)**2*dr(2)**2)
         ux21(i1,i2,i3,kd)= rsxy(i1,i2,i3,0,0)*ur2(i1,i2,i3,kd)
         uy21(i1,i2,i3,kd)=0
         uz21(i1,i2,i3,kd)=0
         ux22(i1,i2,i3,kd)= rsxy(i1,i2,i3,0,0)*ur2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*us2(i1,i2,i3,kd)
         uy22(i1,i2,i3,kd)= rsxy(i1,i2,i3,0,1)*ur2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,1)*us2(i1,i2,i3,kd)
         uz22(i1,i2,i3,kd)=0
         ux23(i1,i2,i3,kd)=rsxy(i1,i2,i3,0,0)*ur2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*us2(i1,i2,i3,kd)+rsxy(i1,i2,i3,2,0)*ut2(i1,i2,i3,kd)
         uy23(i1,i2,i3,kd)=rsxy(i1,i2,i3,0,1)*ur2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,1)*us2(i1,i2,i3,kd)+rsxy(i1,i2,i3,2,1)*ut2(i1,i2,i3,kd)
         uz23(i1,i2,i3,kd)=rsxy(i1,i2,i3,0,2)*ur2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,2)*us2(i1,i2,i3,kd)+rsxy(i1,i2,i3,2,2)*ut2(i1,i2,i3,kd)
         rsxyx21(i1,i2,i3,m,n)= rsxy(i1,i2,i3,0,0)*rsxyr2(i1,i2,i3,m,n)
         rsxyx22(i1,i2,i3,m,n)= rsxy(i1,i2,i3,0,0)*rsxyr2(i1,i2,i3,m,n)+rsxy(i1,i2,i3,1,0)*rsxys2(i1,i2,i3,m,n)
         rsxyy22(i1,i2,i3,m,n)= rsxy(i1,i2,i3,0,1)*rsxyr2(i1,i2,i3,m,n)+rsxy(i1,i2,i3,1,1)*rsxys2(i1,i2,i3,m,n)
         rsxyx23(i1,i2,i3,m,n)=rsxy(i1,i2,i3,0,0)*rsxyr2(i1,i2,i3,m,n)+rsxy(i1,i2,i3,1,0)*rsxys2(i1,i2,i3,m,n)+rsxy(i1,i2,i3,2,0)*rsxyt2(i1,i2,i3,m,n)
         rsxyy23(i1,i2,i3,m,n)=rsxy(i1,i2,i3,0,1)*rsxyr2(i1,i2,i3,m,n)+rsxy(i1,i2,i3,1,1)*rsxys2(i1,i2,i3,m,n)+rsxy(i1,i2,i3,2,1)*rsxyt2(i1,i2,i3,m,n)
         rsxyz23(i1,i2,i3,m,n)=rsxy(i1,i2,i3,0,2)*rsxyr2(i1,i2,i3,m,n)+rsxy(i1,i2,i3,1,2)*rsxys2(i1,i2,i3,m,n)+rsxy(i1,i2,i3,2,2)*rsxyt2(i1,i2,i3,m,n)
         uxx21(i1,i2,i3,kd)=(rsxy(i1,i2,i3,0,0)**2)*urr2(i1,i2,i3,kd)+(rsxyx22(i1,i2,i3,0,0))*ur2(i1,i2,i3,kd)
         uyy21(i1,i2,i3,kd)=0
         uxy21(i1,i2,i3,kd)=0
         uxz21(i1,i2,i3,kd)=0
         uyz21(i1,i2,i3,kd)=0
         uzz21(i1,i2,i3,kd)=0
         ulaplacian21(i1,i2,i3,kd)=uxx21(i1,i2,i3,kd)
         uxx22(i1,i2,i3,kd)=(rsxy(i1,i2,i3,0,0)**2)*urr2(i1,i2,i3,kd)+2.*(rsxy(i1,i2,i3,0,0)*rsxy(i1,i2,i3,1,0))*urs2(i1,i2,i3,kd)+(rsxy(i1,i2,i3,1,0)**2)*uss2(i1,i2,i3,kd)+(rsxyx22(i1,i2,i3,0,0))*ur2(i1,i2,i3,kd)+(rsxyx22(i1,i2,i3,1,0))*us2(i1,i2,i3,kd)
         uyy22(i1,i2,i3,kd)=(rsxy(i1,i2,i3,0,1)**2)*urr2(i1,i2,i3,kd)+2.*(rsxy(i1,i2,i3,0,1)*rsxy(i1,i2,i3,1,1))*urs2(i1,i2,i3,kd)+(rsxy(i1,i2,i3,1,1)**2)*uss2(i1,i2,i3,kd)+(rsxyy22(i1,i2,i3,0,1))*ur2(i1,i2,i3,kd)+(rsxyy22(i1,i2,i3,1,1))*us2(i1,i2,i3,kd)
         uxy22(i1,i2,i3,kd)=rsxy(i1,i2,i3,0,0)*rsxy(i1,i2,i3,0,1)*urr2(i1,i2,i3,kd)+(rsxy(i1,i2,i3,0,0)*rsxy(i1,i2,i3,1,1)+rsxy(i1,i2,i3,0,1)*rsxy(i1,i2,i3,1,0))*urs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*rsxy(i1,i2,i3,1,1)*uss2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,0,1)*ur2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,1)*us2(i1,i2,i3,kd)
         uxz22(i1,i2,i3,kd)=0
         uyz22(i1,i2,i3,kd)=0
         uzz22(i1,i2,i3,kd)=0
         ulaplacian22(i1,i2,i3,kd)=(rsxy(i1,i2,i3,0,0)**2+rsxy(i1,i2,i3,0,1)**2)*urr2(i1,i2,i3,kd)+2.*(rsxy(i1,i2,i3,0,0)*rsxy(i1,i2,i3,1,0)+ rsxy(i1,i2,i3,0,1)*rsxy(i1,i2,i3,1,1))*urs2(i1,i2,i3,kd)+(rsxy(i1,i2,i3,1,0)**2+rsxy(i1,i2,i3,1,1)**2)*uss2(i1,i2,i3,kd)+(rsxyx22(i1,i2,i3,0,0)+rsxyy22(i1,i2,i3,0,1))*ur2(i1,i2,i3,kd)+(rsxyx22(i1,i2,i3,1,0)+rsxyy22(i1,i2,i3,1,1))*us2(i1,i2,i3,kd)
         ! ..... start: 3rd and 4th derivatives, 2D ....
         rsxyxx22(i1,i2,i3,m,n)=(rsxy(i1,i2,i3,0,0)**2)*rsxyrr2(i1,i2,i3,m,n)+2.*(rsxy(i1,i2,i3,0,0)*rsxy(i1,i2,i3,1,0))*rsxyrs2(i1,i2,i3,m,n)+(rsxy(i1,i2,i3,1,0)**2)*rsxyss2(i1,i2,i3,m,n)+(rsxyx22(i1,i2,i3,0,0))*rsxyr2(i1,i2,i3,m,n)+(rsxyx22(i1,i2,i3,1,0))*rsxys2(i1,i2,i3,m,n)
         rsxyyy22(i1,i2,i3,m,n)=(rsxy(i1,i2,i3,0,1)**2)*rsxyrr2(i1,i2,i3,m,n)+2.*(rsxy(i1,i2,i3,0,1)*rsxy(i1,i2,i3,1,1))*rsxyrs2(i1,i2,i3,m,n)+(rsxy(i1,i2,i3,1,1)**2)*rsxyss2(i1,i2,i3,m,n)+(rsxyy22(i1,i2,i3,0,1))*rsxyr2(i1,i2,i3,m,n)+(rsxyy22(i1,i2,i3,1,1))*rsxys2(i1,i2,i3,m,n)
         rsxyxy22(i1,i2,i3,m,n)=rsxy(i1,i2,i3,0,0)*rsxy(i1,i2,i3,0,1)*rsxyrr2(i1,i2,i3,m,n)+(rsxy(i1,i2,i3,0,0)*rsxy(i1,i2,i3,1,1)+rsxy(i1,i2,i3,0,1)*rsxy(i1,i2,i3,1,0))*rsxyrs2(i1,i2,i3,m,n)+rsxy(i1,i2,i3,1,0)*rsxy(i1,i2,i3,1,1)*rsxyss2(i1,i2,i3,m,n)+rsxyx22(i1,i2,i3,0,1)*rsxyr2(i1,i2,i3,m,n)+rsxyx22(i1,i2,i3,1,1)*rsxys2(i1,i2,i3,m,n)
         rsxyxxx22(i1,i2,i3,m,n)=rsxyxx22(i1,i2,i3,0,0)*rsxyr2(i1,i2,i3,m,n)+rsxyxx22(i1,i2,i3,1,0)*rsxys2(i1,i2,i3,m,n)+rsxyx22(i1,i2,i3,0,0)*(rsxy(i1,i2,i3,0,0)*rsxyrr2(i1,i2,i3,m,n)+rsxy(i1,i2,i3,1,0)*rsxyrs2(i1,i2,i3,m,n))+rsxy(i1,i2,i3,0,0)*(rsxyx22(i1,i2,i3,0,0)*rsxyrr2(i1,i2,i3,m,n)+rsxyx22(i1,i2,i3,1,0)*rsxyrs2(i1,i2,i3,m,n))+rsxyx22(i1,i2,i3,1,0)*(rsxy(i1,i2,i3,0,0)*rsxyrs2(i1,i2,i3,m,n)+rsxy(i1,i2,i3,1,0)*rsxyss2(i1,i2,i3,m,n))+rsxy(i1,i2,i3,1,0)*(rsxyx22(i1,i2,i3,0,0)*rsxyrs2(i1,i2,i3,m,n)+rsxyx22(i1,i2,i3,1,0)*rsxyss2(i1,i2,i3,m,n))+rsxy(i1,i2,i3,0,0)*(rsxyx22(i1,i2,i3,0,0)*rsxyrr2(i1,i2,i3,m,n)+rsxyx22(i1,i2,i3,1,0)*rsxyrs2(i1,i2,i3,m,n)+rsxy(i1,i2,i3,0,0)*(rsxy(i1,i2,i3,0,0)*rsxyrrr2(i1,i2,i3,m,n)+rsxy(i1,i2,i3,1,0)*rsxyrrs2(i1,i2,i3,m,n))+rsxy(i1,i2,i3,1,0)*(rsxy(i1,i2,i3,0,0)*rsxyrrs2(i1,i2,i3,m,n)+rsxy(i1,i2,i3,1,0)*rsxyrss2(i1,i2,i3,m,n)))+rsxy(i1,i2,i3,1,0)*(rsxyx22(i1,i2,i3,0,0)*rsxyrs2(i1,i2,i3,m,n)+rsxyx22(i1,i2,i3,1,0)*rsxyss2(i1,i2,i3,m,n)+rsxy(i1,i2,i3,0,0)*(rsxy(i1,i2,i3,0,0)*rsxyrrs2(i1,i2,i3,m,n)+rsxy(i1,i2,i3,1,0)*rsxyrss2(i1,i2,i3,m,n))+rsxy(i1,i2,i3,1,0)*(rsxy(i1,i2,i3,0,0)*rsxyrss2(i1,i2,i3,m,n)+rsxy(i1,i2,i3,1,0)*rsxysss2(i1,i2,i3,m,n)))
         rsxyxxy22(i1,i2,i3,m,n)=rsxyxy22(i1,i2,i3,0,0)*rsxyr2(i1,i2,i3,m,n)+rsxyxy22(i1,i2,i3,1,0)*rsxys2(i1,i2,i3,m,n)+rsxyx22(i1,i2,i3,0,1)*(rsxy(i1,i2,i3,0,0)*rsxyrr2(i1,i2,i3,m,n)+rsxy(i1,i2,i3,1,0)*rsxyrs2(i1,i2,i3,m,n))+rsxy(i1,i2,i3,0,0)*(rsxyx22(i1,i2,i3,0,1)*rsxyrr2(i1,i2,i3,m,n)+rsxyx22(i1,i2,i3,1,1)*rsxyrs2(i1,i2,i3,m,n))+rsxyx22(i1,i2,i3,1,1)*(rsxy(i1,i2,i3,0,0)*rsxyrs2(i1,i2,i3,m,n)+rsxy(i1,i2,i3,1,0)*rsxyss2(i1,i2,i3,m,n))+rsxy(i1,i2,i3,1,0)*(rsxyx22(i1,i2,i3,0,1)*rsxyrs2(i1,i2,i3,m,n)+rsxyx22(i1,i2,i3,1,1)*rsxyss2(i1,i2,i3,m,n))+rsxy(i1,i2,i3,0,1)*(rsxyx22(i1,i2,i3,0,0)*rsxyrr2(i1,i2,i3,m,n)+rsxyx22(i1,i2,i3,1,0)*rsxyrs2(i1,i2,i3,m,n)+rsxy(i1,i2,i3,0,0)*(rsxy(i1,i2,i3,0,0)*rsxyrrr2(i1,i2,i3,m,n)+rsxy(i1,i2,i3,1,0)*rsxyrrs2(i1,i2,i3,m,n))+rsxy(i1,i2,i3,1,0)*(rsxy(i1,i2,i3,0,0)*rsxyrrs2(i1,i2,i3,m,n)+rsxy(i1,i2,i3,1,0)*rsxyrss2(i1,i2,i3,m,n)))+rsxy(i1,i2,i3,1,1)*(rsxyx22(i1,i2,i3,0,0)*rsxyrs2(i1,i2,i3,m,n)+rsxyx22(i1,i2,i3,1,0)*rsxyss2(i1,i2,i3,m,n)+rsxy(i1,i2,i3,0,0)*(rsxy(i1,i2,i3,0,0)*rsxyrrs2(i1,i2,i3,m,n)+rsxy(i1,i2,i3,1,0)*rsxyrss2(i1,i2,i3,m,n))+rsxy(i1,i2,i3,1,0)*(rsxy(i1,i2,i3,0,0)*rsxyrss2(i1,i2,i3,m,n)+rsxy(i1,i2,i3,1,0)*rsxysss2(i1,i2,i3,m,n)))
         rsxyxyy22(i1,i2,i3,m,n)=rsxyyy22(i1,i2,i3,0,0)*rsxyr2(i1,i2,i3,m,n)+rsxyyy22(i1,i2,i3,1,0)*rsxys2(i1,i2,i3,m,n)+rsxyy22(i1,i2,i3,0,1)*(rsxy(i1,i2,i3,0,0)*rsxyrr2(i1,i2,i3,m,n)+rsxy(i1,i2,i3,1,0)*rsxyrs2(i1,i2,i3,m,n))+rsxy(i1,i2,i3,0,1)*(rsxyx22(i1,i2,i3,0,1)*rsxyrr2(i1,i2,i3,m,n)+rsxyx22(i1,i2,i3,1,1)*rsxyrs2(i1,i2,i3,m,n))+rsxyy22(i1,i2,i3,1,1)*(rsxy(i1,i2,i3,0,0)*rsxyrs2(i1,i2,i3,m,n)+rsxy(i1,i2,i3,1,0)*rsxyss2(i1,i2,i3,m,n))+rsxy(i1,i2,i3,1,1)*(rsxyx22(i1,i2,i3,0,1)*rsxyrs2(i1,i2,i3,m,n)+rsxyx22(i1,i2,i3,1,1)*rsxyss2(i1,i2,i3,m,n))+rsxy(i1,i2,i3,0,1)*(rsxyx22(i1,i2,i3,0,1)*rsxyrr2(i1,i2,i3,m,n)+rsxyx22(i1,i2,i3,1,1)*rsxyrs2(i1,i2,i3,m,n)+rsxy(i1,i2,i3,0,1)*(rsxy(i1,i2,i3,0,0)*rsxyrrr2(i1,i2,i3,m,n)+rsxy(i1,i2,i3,1,0)*rsxyrrs2(i1,i2,i3,m,n))+rsxy(i1,i2,i3,1,1)*(rsxy(i1,i2,i3,0,0)*rsxyrrs2(i1,i2,i3,m,n)+rsxy(i1,i2,i3,1,0)*rsxyrss2(i1,i2,i3,m,n)))+rsxy(i1,i2,i3,1,1)*(rsxyx22(i1,i2,i3,0,1)*rsxyrs2(i1,i2,i3,m,n)+rsxyx22(i1,i2,i3,1,1)*rsxyss2(i1,i2,i3,m,n)+rsxy(i1,i2,i3,0,1)*(rsxy(i1,i2,i3,0,0)*rsxyrrs2(i1,i2,i3,m,n)+rsxy(i1,i2,i3,1,0)*rsxyrss2(i1,i2,i3,m,n))+rsxy(i1,i2,i3,1,1)*(rsxy(i1,i2,i3,0,0)*rsxyrss2(i1,i2,i3,m,n)+rsxy(i1,i2,i3,1,0)*rsxysss2(i1,i2,i3,m,n)))
         rsxyyyy22(i1,i2,i3,m,n)=rsxyyy22(i1,i2,i3,0,1)*rsxyr2(i1,i2,i3,m,n)+rsxyyy22(i1,i2,i3,1,1)*rsxys2(i1,i2,i3,m,n)+rsxyy22(i1,i2,i3,0,1)*(rsxy(i1,i2,i3,0,1)*rsxyrr2(i1,i2,i3,m,n)+rsxy(i1,i2,i3,1,1)*rsxyrs2(i1,i2,i3,m,n))+rsxy(i1,i2,i3,0,1)*(rsxyy22(i1,i2,i3,0,1)*rsxyrr2(i1,i2,i3,m,n)+rsxyy22(i1,i2,i3,1,1)*rsxyrs2(i1,i2,i3,m,n))+rsxyy22(i1,i2,i3,1,1)*(rsxy(i1,i2,i3,0,1)*rsxyrs2(i1,i2,i3,m,n)+rsxy(i1,i2,i3,1,1)*rsxyss2(i1,i2,i3,m,n))+rsxy(i1,i2,i3,1,1)*(rsxyy22(i1,i2,i3,0,1)*rsxyrs2(i1,i2,i3,m,n)+rsxyy22(i1,i2,i3,1,1)*rsxyss2(i1,i2,i3,m,n))+rsxy(i1,i2,i3,0,1)*(rsxyy22(i1,i2,i3,0,1)*rsxyrr2(i1,i2,i3,m,n)+rsxyy22(i1,i2,i3,1,1)*rsxyrs2(i1,i2,i3,m,n)+rsxy(i1,i2,i3,0,1)*(rsxy(i1,i2,i3,0,1)*rsxyrrr2(i1,i2,i3,m,n)+rsxy(i1,i2,i3,1,1)*rsxyrrs2(i1,i2,i3,m,n))+rsxy(i1,i2,i3,1,1)*(rsxy(i1,i2,i3,0,1)*rsxyrrs2(i1,i2,i3,m,n)+rsxy(i1,i2,i3,1,1)*rsxyrss2(i1,i2,i3,m,n)))+rsxy(i1,i2,i3,1,1)*(rsxyy22(i1,i2,i3,0,1)*rsxyrs2(i1,i2,i3,m,n)+rsxyy22(i1,i2,i3,1,1)*rsxyss2(i1,i2,i3,m,n)+rsxy(i1,i2,i3,0,1)*(rsxy(i1,i2,i3,0,1)*rsxyrrs2(i1,i2,i3,m,n)+rsxy(i1,i2,i3,1,1)*rsxyrss2(i1,i2,i3,m,n))+rsxy(i1,i2,i3,1,1)*(rsxy(i1,i2,i3,0,1)*rsxyrss2(i1,i2,i3,m,n)+rsxy(i1,i2,i3,1,1)*rsxysss2(i1,i2,i3,m,n)))
         uxxx22(i1,i2,i3,kd)=rsxyxx22(i1,i2,i3,0,0)*ur2(i1,i2,i3,kd)+rsxyxx22(i1,i2,i3,1,0)*us2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,0,0)*(rsxy(i1,i2,i3,0,0)*urr2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urs2(i1,i2,i3,kd))+rsxy(i1,i2,i3,0,0)*(rsxyx22(i1,i2,i3,0,0)*urr2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,0)*urs2(i1,i2,i3,kd))+rsxyx22(i1,i2,i3,1,0)*(rsxy(i1,i2,i3,0,0)*urs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*uss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,0)*(rsxyx22(i1,i2,i3,0,0)*urs2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,0)*uss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,0,0)*(rsxyx22(i1,i2,i3,0,0)*urr2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,0)*urs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,0,0)*(rsxy(i1,i2,i3,0,0)*urrr2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urrs2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,0)*(rsxy(i1,i2,i3,0,0)*urrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urss2(i1,i2,i3,kd)))+rsxy(i1,i2,i3,1,0)*(rsxyx22(i1,i2,i3,0,0)*urs2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,0)*uss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,0,0)*(rsxy(i1,i2,i3,0,0)*urrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,0)*(rsxy(i1,i2,i3,0,0)*urss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*usss2(i1,i2,i3,kd)))
         uxxy22(i1,i2,i3,kd)=rsxyxy22(i1,i2,i3,0,0)*ur2(i1,i2,i3,kd)+rsxyxy22(i1,i2,i3,1,0)*us2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,0,1)*(rsxy(i1,i2,i3,0,0)*urr2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urs2(i1,i2,i3,kd))+rsxy(i1,i2,i3,0,0)*(rsxyx22(i1,i2,i3,0,1)*urr2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,1)*urs2(i1,i2,i3,kd))+rsxyx22(i1,i2,i3,1,1)*(rsxy(i1,i2,i3,0,0)*urs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*uss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,0)*(rsxyx22(i1,i2,i3,0,1)*urs2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,1)*uss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,0,1)*(rsxyx22(i1,i2,i3,0,0)*urr2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,0)*urs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,0,0)*(rsxy(i1,i2,i3,0,0)*urrr2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urrs2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,0)*(rsxy(i1,i2,i3,0,0)*urrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urss2(i1,i2,i3,kd)))+rsxy(i1,i2,i3,1,1)*(rsxyx22(i1,i2,i3,0,0)*urs2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,0)*uss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,0,0)*(rsxy(i1,i2,i3,0,0)*urrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,0)*(rsxy(i1,i2,i3,0,0)*urss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*usss2(i1,i2,i3,kd)))
         uxyy22(i1,i2,i3,kd)=rsxyyy22(i1,i2,i3,0,0)*ur2(i1,i2,i3,kd)+rsxyyy22(i1,i2,i3,1,0)*us2(i1,i2,i3,kd)+rsxyy22(i1,i2,i3,0,1)*(rsxy(i1,i2,i3,0,0)*urr2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urs2(i1,i2,i3,kd))+rsxy(i1,i2,i3,0,1)*(rsxyx22(i1,i2,i3,0,1)*urr2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,1)*urs2(i1,i2,i3,kd))+rsxyy22(i1,i2,i3,1,1)*(rsxy(i1,i2,i3,0,0)*urs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*uss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,1)*(rsxyx22(i1,i2,i3,0,1)*urs2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,1)*uss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,0,1)*(rsxyx22(i1,i2,i3,0,1)*urr2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,1)*urs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,0,1)*(rsxy(i1,i2,i3,0,0)*urrr2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urrs2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,1)*(rsxy(i1,i2,i3,0,0)*urrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urss2(i1,i2,i3,kd)))+rsxy(i1,i2,i3,1,1)*(rsxyx22(i1,i2,i3,0,1)*urs2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,1)*uss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,0,1)*(rsxy(i1,i2,i3,0,0)*urrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,1)*(rsxy(i1,i2,i3,0,0)*urss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*usss2(i1,i2,i3,kd)))
         uyyy22(i1,i2,i3,kd)=rsxyyy22(i1,i2,i3,0,1)*ur2(i1,i2,i3,kd)+rsxyyy22(i1,i2,i3,1,1)*us2(i1,i2,i3,kd)+rsxyy22(i1,i2,i3,0,1)*(rsxy(i1,i2,i3,0,1)*urr2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,1)*urs2(i1,i2,i3,kd))+rsxy(i1,i2,i3,0,1)*(rsxyy22(i1,i2,i3,0,1)*urr2(i1,i2,i3,kd)+rsxyy22(i1,i2,i3,1,1)*urs2(i1,i2,i3,kd))+rsxyy22(i1,i2,i3,1,1)*(rsxy(i1,i2,i3,0,1)*urs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,1)*uss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,1)*(rsxyy22(i1,i2,i3,0,1)*urs2(i1,i2,i3,kd)+rsxyy22(i1,i2,i3,1,1)*uss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,0,1)*(rsxyy22(i1,i2,i3,0,1)*urr2(i1,i2,i3,kd)+rsxyy22(i1,i2,i3,1,1)*urs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,0,1)*(rsxy(i1,i2,i3,0,1)*urrr2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,1)*urrs2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,1)*(rsxy(i1,i2,i3,0,1)*urrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,1)*urss2(i1,i2,i3,kd)))+rsxy(i1,i2,i3,1,1)*(rsxyy22(i1,i2,i3,0,1)*urs2(i1,i2,i3,kd)+rsxyy22(i1,i2,i3,1,1)*uss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,0,1)*(rsxy(i1,i2,i3,0,1)*urrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,1)*urss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,1)*(rsxy(i1,i2,i3,0,1)*urss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,1)*usss2(i1,i2,i3,kd)))
         uxxxx22(i1,i2,i3,kd)=rsxyxxx22(i1,i2,i3,0,0)*ur2(i1,i2,i3,kd)+rsxyxxx22(i1,i2,i3,1,0)*us2(i1,i2,i3,kd)+rsxyxx22(i1,i2,i3,0,0)*(rsxy(i1,i2,i3,0,0)*urr2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urs2(i1,i2,i3,kd))+2*rsxyx22(i1,i2,i3,0,0)*(rsxyx22(i1,i2,i3,0,0)*urr2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,0)*urs2(i1,i2,i3,kd))+rsxy(i1,i2,i3,0,0)*(rsxyxx22(i1,i2,i3,0,0)*urr2(i1,i2,i3,kd)+rsxyxx22(i1,i2,i3,1,0)*urs2(i1,i2,i3,kd))+rsxyxx22(i1,i2,i3,1,0)*(rsxy(i1,i2,i3,0,0)*urs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*uss2(i1,i2,i3,kd))+2*rsxyx22(i1,i2,i3,1,0)*(rsxyx22(i1,i2,i3,0,0)*urs2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,0)*uss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,0)*(rsxyxx22(i1,i2,i3,0,0)*urs2(i1,i2,i3,kd)+rsxyxx22(i1,i2,i3,1,0)*uss2(i1,i2,i3,kd))+rsxyx22(i1,i2,i3,0,0)*(rsxyx22(i1,i2,i3,0,0)*urr2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,0)*urs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,0,0)*(rsxy(i1,i2,i3,0,0)*urrr2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urrs2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,0)*(rsxy(i1,i2,i3,0,0)*urrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urss2(i1,i2,i3,kd)))+rsxy(i1,i2,i3,0,0)*(rsxyxx22(i1,i2,i3,0,0)*urr2(i1,i2,i3,kd)+rsxyxx22(i1,i2,i3,1,0)*urs2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,0,0)*(rsxy(i1,i2,i3,0,0)*urrr2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urrs2(i1,i2,i3,kd))+rsxy(i1,i2,i3,0,0)*(rsxyx22(i1,i2,i3,0,0)*urrr2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,0)*urrs2(i1,i2,i3,kd))+rsxyx22(i1,i2,i3,1,0)*(rsxy(i1,i2,i3,0,0)*urrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,0)*(rsxyx22(i1,i2,i3,0,0)*urrs2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,0)*urss2(i1,i2,i3,kd)))+rsxyx22(i1,i2,i3,1,0)*(rsxyx22(i1,i2,i3,0,0)*urs2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,0)*uss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,0,0)*(rsxy(i1,i2,i3,0,0)*urrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,0)*(rsxy(i1,i2,i3,0,0)*urss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*usss2(i1,i2,i3,kd)))+rsxy(i1,i2,i3,1,0)*(rsxyxx22(i1,i2,i3,0,0)*urs2(i1,i2,i3,kd)+rsxyxx22(i1,i2,i3,1,0)*uss2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,0,0)*(rsxy(i1,i2,i3,0,0)*urrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,0,0)*(rsxyx22(i1,i2,i3,0,0)*urrs2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,0)*urss2(i1,i2,i3,kd))+rsxyx22(i1,i2,i3,1,0)*(rsxy(i1,i2,i3,0,0)*urss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*usss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,0)*(rsxyx22(i1,i2,i3,0,0)*urss2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,0)*usss2(i1,i2,i3,kd)))+rsxy(i1,i2,i3,0,0)*(rsxyxx22(i1,i2,i3,0,0)*urr2(i1,i2,i3,kd)+rsxyxx22(i1,i2,i3,1,0)*urs2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,0,0)*(rsxy(i1,i2,i3,0,0)*urrr2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urrs2(i1,i2,i3,kd))+rsxy(i1,i2,i3,0,0)*(rsxyx22(i1,i2,i3,0,0)*urrr2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,0)*urrs2(i1,i2,i3,kd))+rsxyx22(i1,i2,i3,1,0)*(rsxy(i1,i2,i3,0,0)*urrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,0)*(rsxyx22(i1,i2,i3,0,0)*urrs2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,0)*urss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,0,0)*(rsxyx22(i1,i2,i3,0,0)*urrr2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,0)*urrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,0,0)*(rsxy(i1,i2,i3,0,0)*urrrr2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urrrs2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,0)*(rsxy(i1,i2,i3,0,0)*urrrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urrss2(i1,i2,i3,kd)))+rsxy(i1,i2,i3,1,0)*(rsxyx22(i1,i2,i3,0,0)*urrs2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,0)*urss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,0,0)*(rsxy(i1,i2,i3,0,0)*urrrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urrss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,0)*(rsxy(i1,i2,i3,0,0)*urrss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*ursss2(i1,i2,i3,kd))))+rsxy(i1,i2,i3,1,0)*(rsxyxx22(i1,i2,i3,0,0)*urs2(i1,i2,i3,kd)+rsxyxx22(i1,i2,i3,1,0)*uss2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,0,0)*(rsxy(i1,i2,i3,0,0)*urrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,0,0)*(rsxyx22(i1,i2,i3,0,0)*urrs2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,0)*urss2(i1,i2,i3,kd))+rsxyx22(i1,i2,i3,1,0)*(rsxy(i1,i2,i3,0,0)*urss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*usss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,0)*(rsxyx22(i1,i2,i3,0,0)*urss2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,0)*usss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,0,0)*(rsxyx22(i1,i2,i3,0,0)*urrs2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,0)*urss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,0,0)*(rsxy(i1,i2,i3,0,0)*urrrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urrss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,0)*(rsxy(i1,i2,i3,0,0)*urrss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*ursss2(i1,i2,i3,kd)))+rsxy(i1,i2,i3,1,0)*(rsxyx22(i1,i2,i3,0,0)*urss2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,0)*usss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,0,0)*(rsxy(i1,i2,i3,0,0)*urrss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*ursss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,0)*(rsxy(i1,i2,i3,0,0)*ursss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*ussss2(i1,i2,i3,kd))))
         uxxxy22(i1,i2,i3,kd)=rsxyxxy22(i1,i2,i3,0,0)*ur2(i1,i2,i3,kd)+rsxyxxy22(i1,i2,i3,1,0)*us2(i1,i2,i3,kd)+rsxyxy22(i1,i2,i3,0,0)*(rsxy(i1,i2,i3,0,0)*urr2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urs2(i1,i2,i3,kd))+rsxyx22(i1,i2,i3,0,0)*(rsxyx22(i1,i2,i3,0,1)*urr2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,1)*urs2(i1,i2,i3,kd))+rsxyx22(i1,i2,i3,0,1)*(rsxyx22(i1,i2,i3,0,0)*urr2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,0)*urs2(i1,i2,i3,kd))+rsxy(i1,i2,i3,0,0)*(rsxyxy22(i1,i2,i3,0,0)*urr2(i1,i2,i3,kd)+rsxyxy22(i1,i2,i3,1,0)*urs2(i1,i2,i3,kd))+rsxyxy22(i1,i2,i3,1,0)*(rsxy(i1,i2,i3,0,0)*urs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*uss2(i1,i2,i3,kd))+rsxyx22(i1,i2,i3,1,0)*(rsxyx22(i1,i2,i3,0,1)*urs2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,1)*uss2(i1,i2,i3,kd))+rsxyx22(i1,i2,i3,1,1)*(rsxyx22(i1,i2,i3,0,0)*urs2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,0)*uss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,0)*(rsxyxy22(i1,i2,i3,0,0)*urs2(i1,i2,i3,kd)+rsxyxy22(i1,i2,i3,1,0)*uss2(i1,i2,i3,kd))+rsxyx22(i1,i2,i3,0,1)*(rsxyx22(i1,i2,i3,0,0)*urr2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,0)*urs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,0,0)*(rsxy(i1,i2,i3,0,0)*urrr2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urrs2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,0)*(rsxy(i1,i2,i3,0,0)*urrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urss2(i1,i2,i3,kd)))+rsxy(i1,i2,i3,0,0)*(rsxyxy22(i1,i2,i3,0,0)*urr2(i1,i2,i3,kd)+rsxyxy22(i1,i2,i3,1,0)*urs2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,0,1)*(rsxy(i1,i2,i3,0,0)*urrr2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urrs2(i1,i2,i3,kd))+rsxy(i1,i2,i3,0,0)*(rsxyx22(i1,i2,i3,0,1)*urrr2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,1)*urrs2(i1,i2,i3,kd))+rsxyx22(i1,i2,i3,1,1)*(rsxy(i1,i2,i3,0,0)*urrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,0)*(rsxyx22(i1,i2,i3,0,1)*urrs2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,1)*urss2(i1,i2,i3,kd)))+rsxyx22(i1,i2,i3,1,1)*(rsxyx22(i1,i2,i3,0,0)*urs2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,0)*uss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,0,0)*(rsxy(i1,i2,i3,0,0)*urrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,0)*(rsxy(i1,i2,i3,0,0)*urss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*usss2(i1,i2,i3,kd)))+rsxy(i1,i2,i3,1,0)*(rsxyxy22(i1,i2,i3,0,0)*urs2(i1,i2,i3,kd)+rsxyxy22(i1,i2,i3,1,0)*uss2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,0,1)*(rsxy(i1,i2,i3,0,0)*urrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,0,0)*(rsxyx22(i1,i2,i3,0,1)*urrs2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,1)*urss2(i1,i2,i3,kd))+rsxyx22(i1,i2,i3,1,1)*(rsxy(i1,i2,i3,0,0)*urss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*usss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,0)*(rsxyx22(i1,i2,i3,0,1)*urss2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,1)*usss2(i1,i2,i3,kd)))+rsxy(i1,i2,i3,0,1)*(rsxyxx22(i1,i2,i3,0,0)*urr2(i1,i2,i3,kd)+rsxyxx22(i1,i2,i3,1,0)*urs2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,0,0)*(rsxy(i1,i2,i3,0,0)*urrr2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urrs2(i1,i2,i3,kd))+rsxy(i1,i2,i3,0,0)*(rsxyx22(i1,i2,i3,0,0)*urrr2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,0)*urrs2(i1,i2,i3,kd))+rsxyx22(i1,i2,i3,1,0)*(rsxy(i1,i2,i3,0,0)*urrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,0)*(rsxyx22(i1,i2,i3,0,0)*urrs2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,0)*urss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,0,0)*(rsxyx22(i1,i2,i3,0,0)*urrr2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,0)*urrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,0,0)*(rsxy(i1,i2,i3,0,0)*urrrr2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urrrs2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,0)*(rsxy(i1,i2,i3,0,0)*urrrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urrss2(i1,i2,i3,kd)))+rsxy(i1,i2,i3,1,0)*(rsxyx22(i1,i2,i3,0,0)*urrs2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,0)*urss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,0,0)*(rsxy(i1,i2,i3,0,0)*urrrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urrss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,0)*(rsxy(i1,i2,i3,0,0)*urrss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*ursss2(i1,i2,i3,kd))))+rsxy(i1,i2,i3,1,1)*(rsxyxx22(i1,i2,i3,0,0)*urs2(i1,i2,i3,kd)+rsxyxx22(i1,i2,i3,1,0)*uss2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,0,0)*(rsxy(i1,i2,i3,0,0)*urrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,0,0)*(rsxyx22(i1,i2,i3,0,0)*urrs2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,0)*urss2(i1,i2,i3,kd))+rsxyx22(i1,i2,i3,1,0)*(rsxy(i1,i2,i3,0,0)*urss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*usss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,0)*(rsxyx22(i1,i2,i3,0,0)*urss2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,0)*usss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,0,0)*(rsxyx22(i1,i2,i3,0,0)*urrs2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,0)*urss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,0,0)*(rsxy(i1,i2,i3,0,0)*urrrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urrss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,0)*(rsxy(i1,i2,i3,0,0)*urrss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*ursss2(i1,i2,i3,kd)))+rsxy(i1,i2,i3,1,0)*(rsxyx22(i1,i2,i3,0,0)*urss2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,0)*usss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,0,0)*(rsxy(i1,i2,i3,0,0)*urrss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*ursss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,0)*(rsxy(i1,i2,i3,0,0)*ursss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*ussss2(i1,i2,i3,kd))))
         uxxyy22(i1,i2,i3,kd)=rsxyxyy22(i1,i2,i3,0,0)*ur2(i1,i2,i3,kd)+rsxyxyy22(i1,i2,i3,1,0)*us2(i1,i2,i3,kd)+rsxyyy22(i1,i2,i3,0,0)*(rsxy(i1,i2,i3,0,0)*urr2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urs2(i1,i2,i3,kd))+2*rsxyx22(i1,i2,i3,0,1)*(rsxyx22(i1,i2,i3,0,1)*urr2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,1)*urs2(i1,i2,i3,kd))+rsxy(i1,i2,i3,0,0)*(rsxyyy22(i1,i2,i3,0,0)*urr2(i1,i2,i3,kd)+rsxyyy22(i1,i2,i3,1,0)*urs2(i1,i2,i3,kd))+rsxyyy22(i1,i2,i3,1,0)*(rsxy(i1,i2,i3,0,0)*urs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*uss2(i1,i2,i3,kd))+2*rsxyx22(i1,i2,i3,1,1)*(rsxyx22(i1,i2,i3,0,1)*urs2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,1)*uss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,0)*(rsxyyy22(i1,i2,i3,0,0)*urs2(i1,i2,i3,kd)+rsxyyy22(i1,i2,i3,1,0)*uss2(i1,i2,i3,kd))+rsxyy22(i1,i2,i3,0,1)*(rsxyx22(i1,i2,i3,0,0)*urr2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,0)*urs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,0,0)*(rsxy(i1,i2,i3,0,0)*urrr2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urrs2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,0)*(rsxy(i1,i2,i3,0,0)*urrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urss2(i1,i2,i3,kd)))+rsxy(i1,i2,i3,0,1)*(rsxyxy22(i1,i2,i3,0,0)*urr2(i1,i2,i3,kd)+rsxyxy22(i1,i2,i3,1,0)*urs2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,0,1)*(rsxy(i1,i2,i3,0,0)*urrr2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urrs2(i1,i2,i3,kd))+rsxy(i1,i2,i3,0,0)*(rsxyx22(i1,i2,i3,0,1)*urrr2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,1)*urrs2(i1,i2,i3,kd))+rsxyx22(i1,i2,i3,1,1)*(rsxy(i1,i2,i3,0,0)*urrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,0)*(rsxyx22(i1,i2,i3,0,1)*urrs2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,1)*urss2(i1,i2,i3,kd)))+rsxyy22(i1,i2,i3,1,1)*(rsxyx22(i1,i2,i3,0,0)*urs2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,0)*uss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,0,0)*(rsxy(i1,i2,i3,0,0)*urrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,0)*(rsxy(i1,i2,i3,0,0)*urss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*usss2(i1,i2,i3,kd)))+rsxy(i1,i2,i3,1,1)*(rsxyxy22(i1,i2,i3,0,0)*urs2(i1,i2,i3,kd)+rsxyxy22(i1,i2,i3,1,0)*uss2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,0,1)*(rsxy(i1,i2,i3,0,0)*urrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,0,0)*(rsxyx22(i1,i2,i3,0,1)*urrs2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,1)*urss2(i1,i2,i3,kd))+rsxyx22(i1,i2,i3,1,1)*(rsxy(i1,i2,i3,0,0)*urss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*usss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,0)*(rsxyx22(i1,i2,i3,0,1)*urss2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,1)*usss2(i1,i2,i3,kd)))+rsxy(i1,i2,i3,0,1)*(rsxyxy22(i1,i2,i3,0,0)*urr2(i1,i2,i3,kd)+rsxyxy22(i1,i2,i3,1,0)*urs2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,0,1)*(rsxy(i1,i2,i3,0,0)*urrr2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urrs2(i1,i2,i3,kd))+rsxy(i1,i2,i3,0,0)*(rsxyx22(i1,i2,i3,0,1)*urrr2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,1)*urrs2(i1,i2,i3,kd))+rsxyx22(i1,i2,i3,1,1)*(rsxy(i1,i2,i3,0,0)*urrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,0)*(rsxyx22(i1,i2,i3,0,1)*urrs2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,1)*urss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,0,1)*(rsxyx22(i1,i2,i3,0,0)*urrr2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,0)*urrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,0,0)*(rsxy(i1,i2,i3,0,0)*urrrr2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urrrs2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,0)*(rsxy(i1,i2,i3,0,0)*urrrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urrss2(i1,i2,i3,kd)))+rsxy(i1,i2,i3,1,1)*(rsxyx22(i1,i2,i3,0,0)*urrs2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,0)*urss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,0,0)*(rsxy(i1,i2,i3,0,0)*urrrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urrss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,0)*(rsxy(i1,i2,i3,0,0)*urrss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*ursss2(i1,i2,i3,kd))))+rsxy(i1,i2,i3,1,1)*(rsxyxy22(i1,i2,i3,0,0)*urs2(i1,i2,i3,kd)+rsxyxy22(i1,i2,i3,1,0)*uss2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,0,1)*(rsxy(i1,i2,i3,0,0)*urrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,0,0)*(rsxyx22(i1,i2,i3,0,1)*urrs2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,1)*urss2(i1,i2,i3,kd))+rsxyx22(i1,i2,i3,1,1)*(rsxy(i1,i2,i3,0,0)*urss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*usss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,0)*(rsxyx22(i1,i2,i3,0,1)*urss2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,1)*usss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,0,1)*(rsxyx22(i1,i2,i3,0,0)*urrs2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,0)*urss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,0,0)*(rsxy(i1,i2,i3,0,0)*urrrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urrss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,0)*(rsxy(i1,i2,i3,0,0)*urrss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*ursss2(i1,i2,i3,kd)))+rsxy(i1,i2,i3,1,1)*(rsxyx22(i1,i2,i3,0,0)*urss2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,0)*usss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,0,0)*(rsxy(i1,i2,i3,0,0)*urrss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*ursss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,0)*(rsxy(i1,i2,i3,0,0)*ursss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*ussss2(i1,i2,i3,kd))))
         uxyyy22(i1,i2,i3,kd)=rsxyyyy22(i1,i2,i3,0,0)*ur2(i1,i2,i3,kd)+rsxyyyy22(i1,i2,i3,1,0)*us2(i1,i2,i3,kd)+rsxyyy22(i1,i2,i3,0,1)*(rsxy(i1,i2,i3,0,0)*urr2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urs2(i1,i2,i3,kd))+2*rsxyy22(i1,i2,i3,0,1)*(rsxyx22(i1,i2,i3,0,1)*urr2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,1)*urs2(i1,i2,i3,kd))+rsxy(i1,i2,i3,0,1)*(rsxyyy22(i1,i2,i3,0,0)*urr2(i1,i2,i3,kd)+rsxyyy22(i1,i2,i3,1,0)*urs2(i1,i2,i3,kd))+rsxyyy22(i1,i2,i3,1,1)*(rsxy(i1,i2,i3,0,0)*urs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*uss2(i1,i2,i3,kd))+2*rsxyy22(i1,i2,i3,1,1)*(rsxyx22(i1,i2,i3,0,1)*urs2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,1)*uss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,1)*(rsxyyy22(i1,i2,i3,0,0)*urs2(i1,i2,i3,kd)+rsxyyy22(i1,i2,i3,1,0)*uss2(i1,i2,i3,kd))+rsxyy22(i1,i2,i3,0,1)*(rsxyx22(i1,i2,i3,0,1)*urr2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,1)*urs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,0,1)*(rsxy(i1,i2,i3,0,0)*urrr2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urrs2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,1)*(rsxy(i1,i2,i3,0,0)*urrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urss2(i1,i2,i3,kd)))+rsxy(i1,i2,i3,0,1)*(rsxyyy22(i1,i2,i3,0,0)*urr2(i1,i2,i3,kd)+rsxyyy22(i1,i2,i3,1,0)*urs2(i1,i2,i3,kd)+rsxyy22(i1,i2,i3,0,1)*(rsxy(i1,i2,i3,0,0)*urrr2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urrs2(i1,i2,i3,kd))+rsxy(i1,i2,i3,0,1)*(rsxyx22(i1,i2,i3,0,1)*urrr2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,1)*urrs2(i1,i2,i3,kd))+rsxyy22(i1,i2,i3,1,1)*(rsxy(i1,i2,i3,0,0)*urrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,1)*(rsxyx22(i1,i2,i3,0,1)*urrs2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,1)*urss2(i1,i2,i3,kd)))+rsxyy22(i1,i2,i3,1,1)*(rsxyx22(i1,i2,i3,0,1)*urs2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,1)*uss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,0,1)*(rsxy(i1,i2,i3,0,0)*urrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,1)*(rsxy(i1,i2,i3,0,0)*urss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*usss2(i1,i2,i3,kd)))+rsxy(i1,i2,i3,1,1)*(rsxyyy22(i1,i2,i3,0,0)*urs2(i1,i2,i3,kd)+rsxyyy22(i1,i2,i3,1,0)*uss2(i1,i2,i3,kd)+rsxyy22(i1,i2,i3,0,1)*(rsxy(i1,i2,i3,0,0)*urrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,0,1)*(rsxyx22(i1,i2,i3,0,1)*urrs2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,1)*urss2(i1,i2,i3,kd))+rsxyy22(i1,i2,i3,1,1)*(rsxy(i1,i2,i3,0,0)*urss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*usss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,1)*(rsxyx22(i1,i2,i3,0,1)*urss2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,1)*usss2(i1,i2,i3,kd)))+rsxy(i1,i2,i3,0,1)*(rsxyyy22(i1,i2,i3,0,0)*urr2(i1,i2,i3,kd)+rsxyyy22(i1,i2,i3,1,0)*urs2(i1,i2,i3,kd)+rsxyy22(i1,i2,i3,0,1)*(rsxy(i1,i2,i3,0,0)*urrr2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urrs2(i1,i2,i3,kd))+rsxy(i1,i2,i3,0,1)*(rsxyx22(i1,i2,i3,0,1)*urrr2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,1)*urrs2(i1,i2,i3,kd))+rsxyy22(i1,i2,i3,1,1)*(rsxy(i1,i2,i3,0,0)*urrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,1)*(rsxyx22(i1,i2,i3,0,1)*urrs2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,1)*urss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,0,1)*(rsxyx22(i1,i2,i3,0,1)*urrr2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,1)*urrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,0,1)*(rsxy(i1,i2,i3,0,0)*urrrr2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urrrs2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,1)*(rsxy(i1,i2,i3,0,0)*urrrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urrss2(i1,i2,i3,kd)))+rsxy(i1,i2,i3,1,1)*(rsxyx22(i1,i2,i3,0,1)*urrs2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,1)*urss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,0,1)*(rsxy(i1,i2,i3,0,0)*urrrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urrss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,1)*(rsxy(i1,i2,i3,0,0)*urrss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*ursss2(i1,i2,i3,kd))))+rsxy(i1,i2,i3,1,1)*(rsxyyy22(i1,i2,i3,0,0)*urs2(i1,i2,i3,kd)+rsxyyy22(i1,i2,i3,1,0)*uss2(i1,i2,i3,kd)+rsxyy22(i1,i2,i3,0,1)*(rsxy(i1,i2,i3,0,0)*urrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,0,1)*(rsxyx22(i1,i2,i3,0,1)*urrs2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,1)*urss2(i1,i2,i3,kd))+rsxyy22(i1,i2,i3,1,1)*(rsxy(i1,i2,i3,0,0)*urss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*usss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,1)*(rsxyx22(i1,i2,i3,0,1)*urss2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,1)*usss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,0,1)*(rsxyx22(i1,i2,i3,0,1)*urrs2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,1)*urss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,0,1)*(rsxy(i1,i2,i3,0,0)*urrrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*urrss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,1)*(rsxy(i1,i2,i3,0,0)*urrss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*ursss2(i1,i2,i3,kd)))+rsxy(i1,i2,i3,1,1)*(rsxyx22(i1,i2,i3,0,1)*urss2(i1,i2,i3,kd)+rsxyx22(i1,i2,i3,1,1)*usss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,0,1)*(rsxy(i1,i2,i3,0,0)*urrss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*ursss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,1)*(rsxy(i1,i2,i3,0,0)*ursss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*ussss2(i1,i2,i3,kd))))
         uyyyy22(i1,i2,i3,kd)=rsxyyyy22(i1,i2,i3,0,1)*ur2(i1,i2,i3,kd)+rsxyyyy22(i1,i2,i3,1,1)*us2(i1,i2,i3,kd)+rsxyyy22(i1,i2,i3,0,1)*(rsxy(i1,i2,i3,0,1)*urr2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,1)*urs2(i1,i2,i3,kd))+2*rsxyy22(i1,i2,i3,0,1)*(rsxyy22(i1,i2,i3,0,1)*urr2(i1,i2,i3,kd)+rsxyy22(i1,i2,i3,1,1)*urs2(i1,i2,i3,kd))+rsxy(i1,i2,i3,0,1)*(rsxyyy22(i1,i2,i3,0,1)*urr2(i1,i2,i3,kd)+rsxyyy22(i1,i2,i3,1,1)*urs2(i1,i2,i3,kd))+rsxyyy22(i1,i2,i3,1,1)*(rsxy(i1,i2,i3,0,1)*urs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,1)*uss2(i1,i2,i3,kd))+2*rsxyy22(i1,i2,i3,1,1)*(rsxyy22(i1,i2,i3,0,1)*urs2(i1,i2,i3,kd)+rsxyy22(i1,i2,i3,1,1)*uss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,1)*(rsxyyy22(i1,i2,i3,0,1)*urs2(i1,i2,i3,kd)+rsxyyy22(i1,i2,i3,1,1)*uss2(i1,i2,i3,kd))+rsxyy22(i1,i2,i3,0,1)*(rsxyy22(i1,i2,i3,0,1)*urr2(i1,i2,i3,kd)+rsxyy22(i1,i2,i3,1,1)*urs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,0,1)*(rsxy(i1,i2,i3,0,1)*urrr2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,1)*urrs2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,1)*(rsxy(i1,i2,i3,0,1)*urrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,1)*urss2(i1,i2,i3,kd)))+rsxy(i1,i2,i3,0,1)*(rsxyyy22(i1,i2,i3,0,1)*urr2(i1,i2,i3,kd)+rsxyyy22(i1,i2,i3,1,1)*urs2(i1,i2,i3,kd)+rsxyy22(i1,i2,i3,0,1)*(rsxy(i1,i2,i3,0,1)*urrr2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,1)*urrs2(i1,i2,i3,kd))+rsxy(i1,i2,i3,0,1)*(rsxyy22(i1,i2,i3,0,1)*urrr2(i1,i2,i3,kd)+rsxyy22(i1,i2,i3,1,1)*urrs2(i1,i2,i3,kd))+rsxyy22(i1,i2,i3,1,1)*(rsxy(i1,i2,i3,0,1)*urrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,1)*urss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,1)*(rsxyy22(i1,i2,i3,0,1)*urrs2(i1,i2,i3,kd)+rsxyy22(i1,i2,i3,1,1)*urss2(i1,i2,i3,kd)))+rsxyy22(i1,i2,i3,1,1)*(rsxyy22(i1,i2,i3,0,1)*urs2(i1,i2,i3,kd)+rsxyy22(i1,i2,i3,1,1)*uss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,0,1)*(rsxy(i1,i2,i3,0,1)*urrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,1)*urss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,1)*(rsxy(i1,i2,i3,0,1)*urss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,1)*usss2(i1,i2,i3,kd)))+rsxy(i1,i2,i3,1,1)*(rsxyyy22(i1,i2,i3,0,1)*urs2(i1,i2,i3,kd)+rsxyyy22(i1,i2,i3,1,1)*uss2(i1,i2,i3,kd)+rsxyy22(i1,i2,i3,0,1)*(rsxy(i1,i2,i3,0,1)*urrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,1)*urss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,0,1)*(rsxyy22(i1,i2,i3,0,1)*urrs2(i1,i2,i3,kd)+rsxyy22(i1,i2,i3,1,1)*urss2(i1,i2,i3,kd))+rsxyy22(i1,i2,i3,1,1)*(rsxy(i1,i2,i3,0,1)*urss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,1)*usss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,1)*(rsxyy22(i1,i2,i3,0,1)*urss2(i1,i2,i3,kd)+rsxyy22(i1,i2,i3,1,1)*usss2(i1,i2,i3,kd)))+rsxy(i1,i2,i3,0,1)*(rsxyyy22(i1,i2,i3,0,1)*urr2(i1,i2,i3,kd)+rsxyyy22(i1,i2,i3,1,1)*urs2(i1,i2,i3,kd)+rsxyy22(i1,i2,i3,0,1)*(rsxy(i1,i2,i3,0,1)*urrr2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,1)*urrs2(i1,i2,i3,kd))+rsxy(i1,i2,i3,0,1)*(rsxyy22(i1,i2,i3,0,1)*urrr2(i1,i2,i3,kd)+rsxyy22(i1,i2,i3,1,1)*urrs2(i1,i2,i3,kd))+rsxyy22(i1,i2,i3,1,1)*(rsxy(i1,i2,i3,0,1)*urrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,1)*urss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,1)*(rsxyy22(i1,i2,i3,0,1)*urrs2(i1,i2,i3,kd)+rsxyy22(i1,i2,i3,1,1)*urss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,0,1)*(rsxyy22(i1,i2,i3,0,1)*urrr2(i1,i2,i3,kd)+rsxyy22(i1,i2,i3,1,1)*urrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,0,1)*(rsxy(i1,i2,i3,0,1)*urrrr2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,1)*urrrs2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,1)*(rsxy(i1,i2,i3,0,1)*urrrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,1)*urrss2(i1,i2,i3,kd)))+rsxy(i1,i2,i3,1,1)*(rsxyy22(i1,i2,i3,0,1)*urrs2(i1,i2,i3,kd)+rsxyy22(i1,i2,i3,1,1)*urss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,0,1)*(rsxy(i1,i2,i3,0,1)*urrrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,1)*urrss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,1)*(rsxy(i1,i2,i3,0,1)*urrss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,1)*ursss2(i1,i2,i3,kd))))+rsxy(i1,i2,i3,1,1)*(rsxyyy22(i1,i2,i3,0,1)*urs2(i1,i2,i3,kd)+rsxyyy22(i1,i2,i3,1,1)*uss2(i1,i2,i3,kd)+rsxyy22(i1,i2,i3,0,1)*(rsxy(i1,i2,i3,0,1)*urrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,1)*urss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,0,1)*(rsxyy22(i1,i2,i3,0,1)*urrs2(i1,i2,i3,kd)+rsxyy22(i1,i2,i3,1,1)*urss2(i1,i2,i3,kd))+rsxyy22(i1,i2,i3,1,1)*(rsxy(i1,i2,i3,0,1)*urss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,1)*usss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,1)*(rsxyy22(i1,i2,i3,0,1)*urss2(i1,i2,i3,kd)+rsxyy22(i1,i2,i3,1,1)*usss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,0,1)*(rsxyy22(i1,i2,i3,0,1)*urrs2(i1,i2,i3,kd)+rsxyy22(i1,i2,i3,1,1)*urss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,0,1)*(rsxy(i1,i2,i3,0,1)*urrrs2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,1)*urrss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,1)*(rsxy(i1,i2,i3,0,1)*urrss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,1)*ursss2(i1,i2,i3,kd)))+rsxy(i1,i2,i3,1,1)*(rsxyy22(i1,i2,i3,0,1)*urss2(i1,i2,i3,kd)+rsxyy22(i1,i2,i3,1,1)*usss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,0,1)*(rsxy(i1,i2,i3,0,1)*urrss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,1)*ursss2(i1,i2,i3,kd))+rsxy(i1,i2,i3,1,1)*(rsxy(i1,i2,i3,0,1)*ursss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,1)*ussss2(i1,i2,i3,kd))))
         uLapSq22(i1,i2,i3,kd)=uxxxx22(i1,i2,i3,kd)+uyyyy22(i1,i2,i3,kd)+2.*uxxyy22(i1,i2,i3,kd)
         ! ..... end: 3rd and 4th derivatives, 2D ....
         uxx23(i1,i2,i3,kd)=rsxy(i1,i2,i3,0,0)**2*urr2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)**2*uss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,2,0)**2*utt2(i1,i2,i3,kd)+2.*rsxy(i1,i2,i3,0,0)*rsxy(i1,i2,i3,1,0)*urs2(i1,i2,i3,kd)+2.*rsxy(i1,i2,i3,0,0)*rsxy(i1,i2,i3,2,0)*urt2(i1,i2,i3,kd)+2.*rsxy(i1,i2,i3,1,0)*rsxy(i1,i2,i3,2,0)*ust2(i1,i2,i3,kd)+rsxyx23(i1,i2,i3,0,0)*ur2(i1,i2,i3,kd)+rsxyx23(i1,i2,i3,1,0)*us2(i1,i2,i3,kd)+rsxyx23(i1,i2,i3,2,0)*ut2(i1,i2,i3,kd)
         uyy23(i1,i2,i3,kd)=rsxy(i1,i2,i3,0,1)**2*urr2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,1)**2*uss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,2,1)**2*utt2(i1,i2,i3,kd)+2.*rsxy(i1,i2,i3,0,1)*rsxy(i1,i2,i3,1,1)*urs2(i1,i2,i3,kd)+2.*rsxy(i1,i2,i3,0,1)*rsxy(i1,i2,i3,2,1)*urt2(i1,i2,i3,kd)+2.*rsxy(i1,i2,i3,1,1)*rsxy(i1,i2,i3,2,1)*ust2(i1,i2,i3,kd)+rsxyy23(i1,i2,i3,0,1)*ur2(i1,i2,i3,kd)+rsxyy23(i1,i2,i3,1,1)*us2(i1,i2,i3,kd)+rsxyy23(i1,i2,i3,2,1)*ut2(i1,i2,i3,kd)
         uzz23(i1,i2,i3,kd)=rsxy(i1,i2,i3,0,2)**2*urr2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,2)**2*uss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,2,2)**2*utt2(i1,i2,i3,kd)+2.*rsxy(i1,i2,i3,0,2)*rsxy(i1,i2,i3,1,2)*urs2(i1,i2,i3,kd)+2.*rsxy(i1,i2,i3,0,2)*rsxy(i1,i2,i3,2,2)*urt2(i1,i2,i3,kd)+2.*rsxy(i1,i2,i3,1,2)*rsxy(i1,i2,i3,2,2)*ust2(i1,i2,i3,kd)+rsxyz23(i1,i2,i3,0,2)*ur2(i1,i2,i3,kd)+rsxyz23(i1,i2,i3,1,2)*us2(i1,i2,i3,kd)+rsxyz23(i1,i2,i3,2,2)*ut2(i1,i2,i3,kd)
         uxy23(i1,i2,i3,kd)=rsxy(i1,i2,i3,0,0)*rsxy(i1,i2,i3,0,1)*urr2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*rsxy(i1,i2,i3,1,1)*uss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,2,0)*rsxy(i1,i2,i3,2,1)*utt2(i1,i2,i3,kd)+(rsxy(i1,i2,i3,0,0)*rsxy(i1,i2,i3,1,1)+rsxy(i1,i2,i3,0,1)*rsxy(i1,i2,i3,1,0))*urs2(i1,i2,i3,kd)+(rsxy(i1,i2,i3,0,0)*rsxy(i1,i2,i3,2,1)+rsxy(i1,i2,i3,0,1)*rsxy(i1,i2,i3,2,0))*urt2(i1,i2,i3,kd)+(rsxy(i1,i2,i3,1,0)*rsxy(i1,i2,i3,2,1)+rsxy(i1,i2,i3,1,1)*rsxy(i1,i2,i3,2,0))*ust2(i1,i2,i3,kd)+rsxyx23(i1,i2,i3,0,1)*ur2(i1,i2,i3,kd)+rsxyx23(i1,i2,i3,1,1)*us2(i1,i2,i3,kd)+rsxyx23(i1,i2,i3,2,1)*ut2(i1,i2,i3,kd)
         uxz23(i1,i2,i3,kd)=rsxy(i1,i2,i3,0,0)*rsxy(i1,i2,i3,0,2)*urr2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*rsxy(i1,i2,i3,1,2)*uss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,2,0)*rsxy(i1,i2,i3,2,2)*utt2(i1,i2,i3,kd)+(rsxy(i1,i2,i3,0,0)*rsxy(i1,i2,i3,1,2)+rsxy(i1,i2,i3,0,2)*rsxy(i1,i2,i3,1,0))*urs2(i1,i2,i3,kd)+(rsxy(i1,i2,i3,0,0)*rsxy(i1,i2,i3,2,2)+rsxy(i1,i2,i3,0,2)*rsxy(i1,i2,i3,2,0))*urt2(i1,i2,i3,kd)+(rsxy(i1,i2,i3,1,0)*rsxy(i1,i2,i3,2,2)+rsxy(i1,i2,i3,1,2)*rsxy(i1,i2,i3,2,0))*ust2(i1,i2,i3,kd)+rsxyx23(i1,i2,i3,0,2)*ur2(i1,i2,i3,kd)+rsxyx23(i1,i2,i3,1,2)*us2(i1,i2,i3,kd)+rsxyx23(i1,i2,i3,2,2)*ut2(i1,i2,i3,kd)
         uyz23(i1,i2,i3,kd)=rsxy(i1,i2,i3,0,1)*rsxy(i1,i2,i3,0,2)*urr2(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,1)*rsxy(i1,i2,i3,1,2)*uss2(i1,i2,i3,kd)+rsxy(i1,i2,i3,2,1)*rsxy(i1,i2,i3,2,2)*utt2(i1,i2,i3,kd)+(rsxy(i1,i2,i3,0,1)*rsxy(i1,i2,i3,1,2)+rsxy(i1,i2,i3,0,2)*rsxy(i1,i2,i3,1,1))*urs2(i1,i2,i3,kd)+(rsxy(i1,i2,i3,0,1)*rsxy(i1,i2,i3,2,2)+rsxy(i1,i2,i3,0,2)*rsxy(i1,i2,i3,2,1))*urt2(i1,i2,i3,kd)+(rsxy(i1,i2,i3,1,1)*rsxy(i1,i2,i3,2,2)+rsxy(i1,i2,i3,1,2)*rsxy(i1,i2,i3,2,1))*ust2(i1,i2,i3,kd)+rsxyy23(i1,i2,i3,0,2)*ur2(i1,i2,i3,kd)+rsxyy23(i1,i2,i3,1,2)*us2(i1,i2,i3,kd)+rsxyy23(i1,i2,i3,2,2)*ut2(i1,i2,i3,kd)
         ulaplacian23(i1,i2,i3,kd)=(rsxy(i1,i2,i3,0,0)**2+rsxy(i1,i2,i3,0,1)**2+rsxy(i1,i2,i3,0,2)**2)*urr2(i1,i2,i3,kd)+(rsxy(i1,i2,i3,1,0)**2+rsxy(i1,i2,i3,1,1)**2+rsxy(i1,i2,i3,1,2)**2)*uss2(i1,i2,i3,kd)+(rsxy(i1,i2,i3,2,0)**2+rsxy(i1,i2,i3,2,1)**2+rsxy(i1,i2,i3,2,2)**2)*utt2(i1,i2,i3,kd)+2.*(rsxy(i1,i2,i3,0,0)*rsxy(i1,i2,i3,1,0)+ rsxy(i1,i2,i3,0,1)*rsxy(i1,i2,i3,1,1)+rsxy(i1,i2,i3,0,2)*rsxy(i1,i2,i3,1,2))*urs2(i1,i2,i3,kd)+2.*(rsxy(i1,i2,i3,0,0)*rsxy(i1,i2,i3,2,0)+ rsxy(i1,i2,i3,0,1)*rsxy(i1,i2,i3,2,1)+rsxy(i1,i2,i3,0,2)*rsxy(i1,i2,i3,2,2))*urt2(i1,i2,i3,kd)+2.*(rsxy(i1,i2,i3,1,0)*rsxy(i1,i2,i3,2,0)+ rsxy(i1,i2,i3,1,1)*rsxy(i1,i2,i3,2,1)+rsxy(i1,i2,i3,1,2)*rsxy(i1,i2,i3,2,2))*ust2(i1,i2,i3,kd)+(rsxyx23(i1,i2,i3,0,0)+rsxyy23(i1,i2,i3,0,1)+rsxyz23(i1,i2,i3,0,2))*ur2(i1,i2,i3,kd)+(rsxyx23(i1,i2,i3,1,0)+rsxyy23(i1,i2,i3,1,1)+rsxyz23(i1,i2,i3,1,2))*us2(i1,i2,i3,kd)+(rsxyx23(i1,i2,i3,2,0)+rsxyy23(i1,i2,i3,2,1)+rsxyz23(i1,i2,i3,2,2))*ut2(i1,i2,i3,kd)
         !============================================================================================
         ! Define derivatives for a rectangular grid
         !
         !============================================================================================
         dx12(kd) = 1./(2.*dx(kd))
         dx22(kd) = 1./(dx(kd)**2)
         ux23r(i1,i2,i3,kd)=(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))*dx12(0)
         uy23r(i1,i2,i3,kd)=(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))*dx12(1)
         uz23r(i1,i2,i3,kd)=(u(i1,i2,i3+1,kd)-u(i1,i2,i3-1,kd))*dx12(2)
         uxx23r(i1,i2,i3,kd)=(-2.*u(i1,i2,i3,kd)+(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd)) )*dx22(0)
         uyy23r(i1,i2,i3,kd)=(-2.*u(i1,i2,i3,kd)+(u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd)) )*dx22(1)
         uxy23r(i1,i2,i3,kd)=(ux23r(i1,i2+1,i3,kd)-ux23r(i1,i2-1,i3,kd))*dx12(1)
         uzz23r(i1,i2,i3,kd)=(-2.*u(i1,i2,i3,kd)+(u(i1,i2,i3+1,kd)+u(i1,i2,i3-1,kd)) )*dx22(2)
         uxz23r(i1,i2,i3,kd)=(ux23r(i1,i2,i3+1,kd)-ux23r(i1,i2,i3-1,kd))*dx12(2)
         uyz23r(i1,i2,i3,kd)=(uy23r(i1,i2,i3+1,kd)-uy23r(i1,i2,i3-1,kd))*dx12(2)
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
         uxxx22r(i1,i2,i3,kd)=(-2.*(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))+(u(i1+2,i2,i3,kd)-u(i1-2,i2,i3,kd)) )*dx22(0)*dx12(0)
         uyyy22r(i1,i2,i3,kd)=(-2.*(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))+(u(i1,i2+2,i3,kd)-u(i1,i2-2,i3,kd)) )*dx22(1)*dx12(1)
         uxxy22r(i1,i2,i3,kd)=( uxx22r(i1,i2+1,i3,kd)-uxx22r(i1,i2-1,i3,kd))/(2.*dx(1))
         uxyy22r(i1,i2,i3,kd)=( uyy22r(i1+1,i2,i3,kd)-uyy22r(i1-1,i2,i3,kd))/(2.*dx(0))
         uxxxx22r(i1,i2,i3,kd)=(6.*u(i1,i2,i3,kd)-4.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd)) +(u(i1+2,i2,i3,kd)+u(i1-2,i2,i3,kd)) )/(dx(0)**4)
         uyyyy22r(i1,i2,i3,kd)=(6.*u(i1,i2,i3,kd)-4.*(u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd)) +(u(i1,i2+2,i3,kd)+u(i1,i2-2,i3,kd)) )/(dx(1)**4)
         uxxyy22r(i1,i2,i3,kd)=( 4.*u(i1,i2,i3,kd)-2.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd)+u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd))+   (u(i1+1,i2+1,i3,kd)+u(i1-1,i2+1,i3,kd)+u(i1+1,i2-1,i3,kd)+u(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)
         uLapSq22r(i1,i2,i3,kd)= ( 6.*u(i1,i2,i3,kd)- 4.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd))+(u(i1+2,i2,i3,kd)+u(i1-2,i2,i3,kd)) )/(dx(0)**4)+( 6.*u(i1,i2,i3,kd)-4.*(u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd)) +(u(i1,i2+2,i3,kd)+u(i1,i2-2,i3,kd)) )/(dx(1)**4)+( 8.*u(i1,i2,i3,kd)-4.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd)+u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd))+2.*(u(i1+1,i2+1,i3,kd)+u(i1-1,i2+1,i3,kd)+u(i1+1,i2-1,i3,kd)+u(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)
         uxxx23r(i1,i2,i3,kd)=(-2.*(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))+(u(i1+2,i2,i3,kd)-u(i1-2,i2,i3,kd)) )*dx22(0)*dx12(0)
         uyyy23r(i1,i2,i3,kd)=(-2.*(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))+(u(i1,i2+2,i3,kd)-u(i1,i2-2,i3,kd)) )*dx22(1)*dx12(1)
         uzzz23r(i1,i2,i3,kd)=(-2.*(u(i1,i2,i3+1,kd)-u(i1,i2,i3-1,kd))+(u(i1,i2,i3+2,kd)-u(i1,i2,i3-2,kd)) )*dx22(1)*dx12(2)
         uxxy23r(i1,i2,i3,kd)=( uxx22r(i1,i2+1,i3,kd)-uxx22r(i1,i2-1,i3,kd))/(2.*dx(1))
         uxyy23r(i1,i2,i3,kd)=( uyy22r(i1+1,i2,i3,kd)-uyy22r(i1-1,i2,i3,kd))/(2.*dx(0))
         uxxz23r(i1,i2,i3,kd)=( uxx22r(i1,i2,i3+1,kd)-uxx22r(i1,i2,i3-1,kd))/(2.*dx(2))
         uyyz23r(i1,i2,i3,kd)=( uyy22r(i1,i2,i3+1,kd)-uyy22r(i1,i2,i3-1,kd))/(2.*dx(2))
         uxzz23r(i1,i2,i3,kd)=( uzz22r(i1+1,i2,i3,kd)-uzz22r(i1-1,i2,i3,kd))/(2.*dx(0))
         uyzz23r(i1,i2,i3,kd)=( uzz22r(i1,i2+1,i3,kd)-uzz22r(i1,i2-1,i3,kd))/(2.*dx(1))
         uxxxx23r(i1,i2,i3,kd)=(6.*u(i1,i2,i3,kd)-4.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd))+(u(i1+2,i2,i3,kd)+u(i1-2,i2,i3,kd)) )/(dx(0)**4)
         uyyyy23r(i1,i2,i3,kd)=(6.*u(i1,i2,i3,kd)-4.*(u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd))+(u(i1,i2+2,i3,kd)+u(i1,i2-2,i3,kd)) )/(dx(1)**4)
         uzzzz23r(i1,i2,i3,kd)=(6.*u(i1,i2,i3,kd)-4.*(u(i1,i2,i3+1,kd)+u(i1,i2,i3-1,kd))+(u(i1,i2,i3+2,kd)+u(i1,i2,i3-2,kd)) )/(dx(2)**4)
         uxxyy23r(i1,i2,i3,kd)=( 4.*u(i1,i2,i3,kd)-2.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd)+u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd))+   (u(i1+1,i2+1,i3,kd)+u(i1-1,i2+1,i3,kd)+u(i1+1,i2-1,i3,kd)+u(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)
         uxxzz23r(i1,i2,i3,kd)=( 4.*u(i1,i2,i3,kd)-2.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd)+u(i1,i2,i3+1,kd)+u(i1,i2,i3-1,kd))+   (u(i1+1,i2,i3+1,kd)+u(i1-1,i2,i3+1,kd)+u(i1+1,i2,i3-1,kd)+u(i1-1,i2,i3-1,kd)) )/(dx(0)**2*dx(2)**2)
         uyyzz23r(i1,i2,i3,kd)=( 4.*u(i1,i2,i3,kd)-2.*(u(i1,i2+1,i3,kd)  +u(i1,i2-1,i3,kd)+  u(i1,i2  ,i3+1,kd)+u(i1,i2  ,i3-1,kd))+   (u(i1,i2+1,i3+1,kd)+u(i1,i2-1,i3+1,kd)+u(i1,i2+1,i3-1,kd)+u(i1,i2-1,i3-1,kd)) )/(dx(1)**2*dx(2)**2)
         dr14(kd) = 1./(12.*dr(kd))
         dr24(kd) = 1./(12.*dr(kd)**2)
         ur4(i1,i2,i3,kd)=(8.*(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))-(u(i1+2,i2,i3,kd)-u(i1-2,i2,i3,kd)))*dr14(0)
         us4(i1,i2,i3,kd)=(8.*(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))-(u(i1,i2+2,i3,kd)-u(i1,i2-2,i3,kd)))*dr14(1)
         ut4(i1,i2,i3,kd)=(8.*(u(i1,i2,i3+1,kd)-u(i1,i2,i3-1,kd))-(u(i1,i2,i3+2,kd)-u(i1,i2,i3-2,kd)))*dr14(2)
         urr4(i1,i2,i3,kd)=(-30.*u(i1,i2,i3,kd)+16.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd))-(u(i1+2,i2,i3,kd)+u(i1-2,i2,i3,kd)) )*dr24(0)
         uss4(i1,i2,i3,kd)=(-30.*u(i1,i2,i3,kd)+16.*(u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd))-(u(i1,i2+2,i3,kd)+u(i1,i2-2,i3,kd)) )*dr24(1)
         utt4(i1,i2,i3,kd)=(-30.*u(i1,i2,i3,kd)+16.*(u(i1,i2,i3+1,kd)+u(i1,i2,i3-1,kd))-(u(i1,i2,i3+2,kd)+u(i1,i2,i3-2,kd)) )*dr24(2)
         urs4(i1,i2,i3,kd)=(8.*(ur4(i1,i2+1,i3,kd)-ur4(i1,i2-1,i3,kd))-(ur4(i1,i2+2,i3,kd)-ur4(i1,i2-2,i3,kd)))*dr14(1)
         urt4(i1,i2,i3,kd)=(8.*(ur4(i1,i2,i3+1,kd)-ur4(i1,i2,i3-1,kd))-(ur4(i1,i2,i3+2,kd)-ur4(i1,i2,i3-2,kd)))*dr14(2)
         ust4(i1,i2,i3,kd)=(8.*(us4(i1,i2,i3+1,kd)-us4(i1,i2,i3-1,kd))-(us4(i1,i2,i3+2,kd)-us4(i1,i2,i3-2,kd)))*dr14(2)
         rsxyr4(i1,i2,i3,m,n)=(8.*(rsxy(i1+1,i2,i3,m,n)-rsxy(i1-1,i2,i3,m,n))-(rsxy(i1+2,i2,i3,m,n)-rsxy(i1-2,i2,i3,m,n)))*dr14(0)
         rsxys4(i1,i2,i3,m,n)=(8.*(rsxy(i1,i2+1,i3,m,n)-rsxy(i1,i2-1,i3,m,n))-(rsxy(i1,i2+2,i3,m,n)-rsxy(i1,i2-2,i3,m,n)))*dr14(1)
         rsxyt4(i1,i2,i3,m,n)=(8.*(rsxy(i1,i2,i3+1,m,n)-rsxy(i1,i2,i3-1,m,n))-(rsxy(i1,i2,i3+2,m,n)-rsxy(i1,i2,i3-2,m,n)))*dr14(2)
         ux41(i1,i2,i3,kd)= rsxy(i1,i2,i3,0,0)*ur4(i1,i2,i3,kd)
         uy41(i1,i2,i3,kd)=0
         uz41(i1,i2,i3,kd)=0
         ux42(i1,i2,i3,kd)= rsxy(i1,i2,i3,0,0)*ur4(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*us4(i1,i2,i3,kd)
         uy42(i1,i2,i3,kd)= rsxy(i1,i2,i3,0,1)*ur4(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,1)*us4(i1,i2,i3,kd)
         uz42(i1,i2,i3,kd)=0
         ux43(i1,i2,i3,kd)=rsxy(i1,i2,i3,0,0)*ur4(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*us4(i1,i2,i3,kd)+rsxy(i1,i2,i3,2,0)*ut4(i1,i2,i3,kd)
         uy43(i1,i2,i3,kd)=rsxy(i1,i2,i3,0,1)*ur4(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,1)*us4(i1,i2,i3,kd)+rsxy(i1,i2,i3,2,1)*ut4(i1,i2,i3,kd)
         uz43(i1,i2,i3,kd)=rsxy(i1,i2,i3,0,2)*ur4(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,2)*us4(i1,i2,i3,kd)+rsxy(i1,i2,i3,2,2)*ut4(i1,i2,i3,kd)
         rsxyx41(i1,i2,i3,m,n)= rsxy(i1,i2,i3,0,0)*rsxyr4(i1,i2,i3,m,n)
         rsxyx42(i1,i2,i3,m,n)= rsxy(i1,i2,i3,0,0)*rsxyr4(i1,i2,i3,m,n)+rsxy(i1,i2,i3,1,0)*rsxys4(i1,i2,i3,m,n)
         rsxyy42(i1,i2,i3,m,n)= rsxy(i1,i2,i3,0,1)*rsxyr4(i1,i2,i3,m,n)+rsxy(i1,i2,i3,1,1)*rsxys4(i1,i2,i3,m,n)
         rsxyx43(i1,i2,i3,m,n)=rsxy(i1,i2,i3,0,0)*rsxyr4(i1,i2,i3,m,n)+rsxy(i1,i2,i3,1,0)*rsxys4(i1,i2,i3,m,n)+rsxy(i1,i2,i3,2,0)*rsxyt4(i1,i2,i3,m,n)
         rsxyy43(i1,i2,i3,m,n)=rsxy(i1,i2,i3,0,1)*rsxyr4(i1,i2,i3,m,n)+rsxy(i1,i2,i3,1,1)*rsxys4(i1,i2,i3,m,n)+rsxy(i1,i2,i3,2,1)*rsxyt4(i1,i2,i3,m,n)
         rsxyz43(i1,i2,i3,m,n)=rsxy(i1,i2,i3,0,2)*rsxyr4(i1,i2,i3,m,n)+rsxy(i1,i2,i3,1,2)*rsxys4(i1,i2,i3,m,n)+rsxy(i1,i2,i3,2,2)*rsxyt4(i1,i2,i3,m,n)
         uxx41(i1,i2,i3,kd)=(rsxy(i1,i2,i3,0,0)**2)*urr4(i1,i2,i3,kd)+(rsxyx42(i1,i2,i3,0,0))*ur4(i1,i2,i3,kd)
         uyy41(i1,i2,i3,kd)=0
         uxy41(i1,i2,i3,kd)=0
         uxz41(i1,i2,i3,kd)=0
         uyz41(i1,i2,i3,kd)=0
         uzz41(i1,i2,i3,kd)=0
         ulaplacian41(i1,i2,i3,kd)=uxx41(i1,i2,i3,kd)
         uxx42(i1,i2,i3,kd)=(rsxy(i1,i2,i3,0,0)**2)*urr4(i1,i2,i3,kd)+2.*(rsxy(i1,i2,i3,0,0)*rsxy(i1,i2,i3,1,0))*urs4(i1,i2,i3,kd)+(rsxy(i1,i2,i3,1,0)**2)*uss4(i1,i2,i3,kd)+(rsxyx42(i1,i2,i3,0,0))*ur4(i1,i2,i3,kd)+(rsxyx42(i1,i2,i3,1,0))*us4(i1,i2,i3,kd)
         uyy42(i1,i2,i3,kd)=(rsxy(i1,i2,i3,0,1)**2)*urr4(i1,i2,i3,kd)+2.*(rsxy(i1,i2,i3,0,1)*rsxy(i1,i2,i3,1,1))*urs4(i1,i2,i3,kd)+(rsxy(i1,i2,i3,1,1)**2)*uss4(i1,i2,i3,kd)+(rsxyy42(i1,i2,i3,0,1))*ur4(i1,i2,i3,kd)+(rsxyy42(i1,i2,i3,1,1))*us4(i1,i2,i3,kd)
         uxy42(i1,i2,i3,kd)=rsxy(i1,i2,i3,0,0)*rsxy(i1,i2,i3,0,1)*urr4(i1,i2,i3,kd)+(rsxy(i1,i2,i3,0,0)*rsxy(i1,i2,i3,1,1)+rsxy(i1,i2,i3,0,1)*rsxy(i1,i2,i3,1,0))*urs4(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*rsxy(i1,i2,i3,1,1)*uss4(i1,i2,i3,kd)+rsxyx42(i1,i2,i3,0,1)*ur4(i1,i2,i3,kd)+rsxyx42(i1,i2,i3,1,1)*us4(i1,i2,i3,kd)
         uxz42(i1,i2,i3,kd)=0
         uyz42(i1,i2,i3,kd)=0
         uzz42(i1,i2,i3,kd)=0
         ulaplacian42(i1,i2,i3,kd)=(rsxy(i1,i2,i3,0,0)**2+rsxy(i1,i2,i3,0,1)**2)*urr4(i1,i2,i3,kd)+2.*(rsxy(i1,i2,i3,0,0)*rsxy(i1,i2,i3,1,0)+ rsxy(i1,i2,i3,0,1)*rsxy(i1,i2,i3,1,1))*urs4(i1,i2,i3,kd)+(rsxy(i1,i2,i3,1,0)**2+rsxy(i1,i2,i3,1,1)**2)*uss4(i1,i2,i3,kd)+(rsxyx42(i1,i2,i3,0,0)+rsxyy42(i1,i2,i3,0,1))*ur4(i1,i2,i3,kd)+(rsxyx42(i1,i2,i3,1,0)+rsxyy42(i1,i2,i3,1,1))*us4(i1,i2,i3,kd)
         uxx43(i1,i2,i3,kd)=rsxy(i1,i2,i3,0,0)**2*urr4(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)**2*uss4(i1,i2,i3,kd)+rsxy(i1,i2,i3,2,0)**2*utt4(i1,i2,i3,kd)+2.*rsxy(i1,i2,i3,0,0)*rsxy(i1,i2,i3,1,0)*urs4(i1,i2,i3,kd)+2.*rsxy(i1,i2,i3,0,0)*rsxy(i1,i2,i3,2,0)*urt4(i1,i2,i3,kd)+2.*rsxy(i1,i2,i3,1,0)*rsxy(i1,i2,i3,2,0)*ust4(i1,i2,i3,kd)+rsxyx43(i1,i2,i3,0,0)*ur4(i1,i2,i3,kd)+rsxyx43(i1,i2,i3,1,0)*us4(i1,i2,i3,kd)+rsxyx43(i1,i2,i3,2,0)*ut4(i1,i2,i3,kd)
         uyy43(i1,i2,i3,kd)=rsxy(i1,i2,i3,0,1)**2*urr4(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,1)**2*uss4(i1,i2,i3,kd)+rsxy(i1,i2,i3,2,1)**2*utt4(i1,i2,i3,kd)+2.*rsxy(i1,i2,i3,0,1)*rsxy(i1,i2,i3,1,1)*urs4(i1,i2,i3,kd)+2.*rsxy(i1,i2,i3,0,1)*rsxy(i1,i2,i3,2,1)*urt4(i1,i2,i3,kd)+2.*rsxy(i1,i2,i3,1,1)*rsxy(i1,i2,i3,2,1)*ust4(i1,i2,i3,kd)+rsxyy43(i1,i2,i3,0,1)*ur4(i1,i2,i3,kd)+rsxyy43(i1,i2,i3,1,1)*us4(i1,i2,i3,kd)+rsxyy43(i1,i2,i3,2,1)*ut4(i1,i2,i3,kd)
         uzz43(i1,i2,i3,kd)=rsxy(i1,i2,i3,0,2)**2*urr4(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,2)**2*uss4(i1,i2,i3,kd)+rsxy(i1,i2,i3,2,2)**2*utt4(i1,i2,i3,kd)+2.*rsxy(i1,i2,i3,0,2)*rsxy(i1,i2,i3,1,2)*urs4(i1,i2,i3,kd)+2.*rsxy(i1,i2,i3,0,2)*rsxy(i1,i2,i3,2,2)*urt4(i1,i2,i3,kd)+2.*rsxy(i1,i2,i3,1,2)*rsxy(i1,i2,i3,2,2)*ust4(i1,i2,i3,kd)+rsxyz43(i1,i2,i3,0,2)*ur4(i1,i2,i3,kd)+rsxyz43(i1,i2,i3,1,2)*us4(i1,i2,i3,kd)+rsxyz43(i1,i2,i3,2,2)*ut4(i1,i2,i3,kd)
         uxy43(i1,i2,i3,kd)=rsxy(i1,i2,i3,0,0)*rsxy(i1,i2,i3,0,1)*urr4(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*rsxy(i1,i2,i3,1,1)*uss4(i1,i2,i3,kd)+rsxy(i1,i2,i3,2,0)*rsxy(i1,i2,i3,2,1)*utt4(i1,i2,i3,kd)+(rsxy(i1,i2,i3,0,0)*rsxy(i1,i2,i3,1,1)+rsxy(i1,i2,i3,0,1)*rsxy(i1,i2,i3,1,0))*urs4(i1,i2,i3,kd)+(rsxy(i1,i2,i3,0,0)*rsxy(i1,i2,i3,2,1)+rsxy(i1,i2,i3,0,1)*rsxy(i1,i2,i3,2,0))*urt4(i1,i2,i3,kd)+(rsxy(i1,i2,i3,1,0)*rsxy(i1,i2,i3,2,1)+rsxy(i1,i2,i3,1,1)*rsxy(i1,i2,i3,2,0))*ust4(i1,i2,i3,kd)+rsxyx43(i1,i2,i3,0,1)*ur4(i1,i2,i3,kd)+rsxyx43(i1,i2,i3,1,1)*us4(i1,i2,i3,kd)+rsxyx43(i1,i2,i3,2,1)*ut4(i1,i2,i3,kd)
         uxz43(i1,i2,i3,kd)=rsxy(i1,i2,i3,0,0)*rsxy(i1,i2,i3,0,2)*urr4(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,0)*rsxy(i1,i2,i3,1,2)*uss4(i1,i2,i3,kd)+rsxy(i1,i2,i3,2,0)*rsxy(i1,i2,i3,2,2)*utt4(i1,i2,i3,kd)+(rsxy(i1,i2,i3,0,0)*rsxy(i1,i2,i3,1,2)+rsxy(i1,i2,i3,0,2)*rsxy(i1,i2,i3,1,0))*urs4(i1,i2,i3,kd)+(rsxy(i1,i2,i3,0,0)*rsxy(i1,i2,i3,2,2)+rsxy(i1,i2,i3,0,2)*rsxy(i1,i2,i3,2,0))*urt4(i1,i2,i3,kd)+(rsxy(i1,i2,i3,1,0)*rsxy(i1,i2,i3,2,2)+rsxy(i1,i2,i3,1,2)*rsxy(i1,i2,i3,2,0))*ust4(i1,i2,i3,kd)+rsxyx43(i1,i2,i3,0,2)*ur4(i1,i2,i3,kd)+rsxyx43(i1,i2,i3,1,2)*us4(i1,i2,i3,kd)+rsxyx43(i1,i2,i3,2,2)*ut4(i1,i2,i3,kd)
         uyz43(i1,i2,i3,kd)=rsxy(i1,i2,i3,0,1)*rsxy(i1,i2,i3,0,2)*urr4(i1,i2,i3,kd)+rsxy(i1,i2,i3,1,1)*rsxy(i1,i2,i3,1,2)*uss4(i1,i2,i3,kd)+rsxy(i1,i2,i3,2,1)*rsxy(i1,i2,i3,2,2)*utt4(i1,i2,i3,kd)+(rsxy(i1,i2,i3,0,1)*rsxy(i1,i2,i3,1,2)+rsxy(i1,i2,i3,0,2)*rsxy(i1,i2,i3,1,1))*urs4(i1,i2,i3,kd)+(rsxy(i1,i2,i3,0,1)*rsxy(i1,i2,i3,2,2)+rsxy(i1,i2,i3,0,2)*rsxy(i1,i2,i3,2,1))*urt4(i1,i2,i3,kd)+(rsxy(i1,i2,i3,1,1)*rsxy(i1,i2,i3,2,2)+rsxy(i1,i2,i3,1,2)*rsxy(i1,i2,i3,2,1))*ust4(i1,i2,i3,kd)+rsxyy43(i1,i2,i3,0,2)*ur4(i1,i2,i3,kd)+rsxyy43(i1,i2,i3,1,2)*us4(i1,i2,i3,kd)+rsxyy43(i1,i2,i3,2,2)*ut4(i1,i2,i3,kd)
         ulaplacian43(i1,i2,i3,kd)=(rsxy(i1,i2,i3,0,0)**2+rsxy(i1,i2,i3,0,1)**2+rsxy(i1,i2,i3,0,2)**2)*urr4(i1,i2,i3,kd)+(rsxy(i1,i2,i3,1,0)**2+rsxy(i1,i2,i3,1,1)**2+rsxy(i1,i2,i3,1,2)**2)*uss4(i1,i2,i3,kd)+(rsxy(i1,i2,i3,2,0)**2+rsxy(i1,i2,i3,2,1)**2+rsxy(i1,i2,i3,2,2)**2)*utt4(i1,i2,i3,kd)+2.*(rsxy(i1,i2,i3,0,0)*rsxy(i1,i2,i3,1,0)+ rsxy(i1,i2,i3,0,1)*rsxy(i1,i2,i3,1,1)+rsxy(i1,i2,i3,0,2)*rsxy(i1,i2,i3,1,2))*urs4(i1,i2,i3,kd)+2.*(rsxy(i1,i2,i3,0,0)*rsxy(i1,i2,i3,2,0)+ rsxy(i1,i2,i3,0,1)*rsxy(i1,i2,i3,2,1)+rsxy(i1,i2,i3,0,2)*rsxy(i1,i2,i3,2,2))*urt4(i1,i2,i3,kd)+2.*(rsxy(i1,i2,i3,1,0)*rsxy(i1,i2,i3,2,0)+ rsxy(i1,i2,i3,1,1)*rsxy(i1,i2,i3,2,1)+rsxy(i1,i2,i3,1,2)*rsxy(i1,i2,i3,2,2))*ust4(i1,i2,i3,kd)+(rsxyx43(i1,i2,i3,0,0)+rsxyy43(i1,i2,i3,0,1)+rsxyz43(i1,i2,i3,0,2))*ur4(i1,i2,i3,kd)+(rsxyx43(i1,i2,i3,1,0)+rsxyy43(i1,i2,i3,1,1)+rsxyz43(i1,i2,i3,1,2))*us4(i1,i2,i3,kd)+(rsxyx43(i1,i2,i3,2,0)+rsxyy43(i1,i2,i3,2,1)+rsxyz43(i1,i2,i3,2,2))*ut4(i1,i2,i3,kd)
         !============================================================================================
         ! Define derivatives for a rectangular grid
         !
         !============================================================================================
         dx41(kd) = 1./(12.*dx(kd))
         dx42(kd) = 1./(12.*dx(kd)**2)
         ux43r(i1,i2,i3,kd)=(8.*(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))-(u(i1+2,i2,i3,kd)-u(i1-2,i2,i3,kd)))*dx41(0)
         uy43r(i1,i2,i3,kd)=(8.*(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))-(u(i1,i2+2,i3,kd)-u(i1,i2-2,i3,kd)))*dx41(1)
         uz43r(i1,i2,i3,kd)=(8.*(u(i1,i2,i3+1,kd)-u(i1,i2,i3-1,kd))-(u(i1,i2,i3+2,kd)-u(i1,i2,i3-2,kd)))*dx41(2)
         uxx43r(i1,i2,i3,kd)=( -30.*u(i1,i2,i3,kd)+16.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd))-(u(i1+2,i2,i3,kd)+u(i1-2,i2,i3,kd)) )*dx42(0) 
         uyy43r(i1,i2,i3,kd)=( -30.*u(i1,i2,i3,kd)+16.*(u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd))-(u(i1,i2+2,i3,kd)+u(i1,i2-2,i3,kd)) )*dx42(1) 
         uzz43r(i1,i2,i3,kd)=( -30.*u(i1,i2,i3,kd)+16.*(u(i1,i2,i3+1,kd)+u(i1,i2,i3-1,kd))-(u(i1,i2,i3+2,kd)+u(i1,i2,i3-2,kd)) )*dx42(2)
         uxy43r(i1,i2,i3,kd)=( (u(i1+2,i2+2,i3,kd)-u(i1-2,i2+2,i3,kd)- u(i1+2,i2-2,i3,kd)+u(i1-2,i2-2,i3,kd)) +8.*(u(i1-1,i2+2,i3,kd)-u(i1-1,i2-2,i3,kd)-u(i1+1,i2+2,i3,kd)+u(i1+1,i2-2,i3,kd) +u(i1+2,i2-1,i3,kd)-u(i1-2,i2-1,i3,kd)-u(i1+2,i2+1,i3,kd)+u(i1-2,i2+1,i3,kd))+64.*(u(i1+1,i2+1,i3,kd)-u(i1-1,i2+1,i3,kd)- u(i1+1,i2-1,i3,kd)+u(i1-1,i2-1,i3,kd)))*(dx41(0)*dx41(1))
         uxz43r(i1,i2,i3,kd)=( (u(i1+2,i2,i3+2,kd)-u(i1-2,i2,i3+2,kd)-u(i1+2,i2,i3-2,kd)+u(i1-2,i2,i3-2,kd)) +8.*(u(i1-1,i2,i3+2,kd)-u(i1-1,i2,i3-2,kd)-u(i1+1,i2,i3+2,kd)+u(i1+1,i2,i3-2,kd) +u(i1+2,i2,i3-1,kd)-u(i1-2,i2,i3-1,kd)- u(i1+2,i2,i3+1,kd)+u(i1-2,i2,i3+1,kd)) +64.*(u(i1+1,i2,i3+1,kd)-u(i1-1,i2,i3+1,kd)-u(i1+1,i2,i3-1,kd)+u(i1-1,i2,i3-1,kd)) )*(dx41(0)*dx41(2))
         uyz43r(i1,i2,i3,kd)=( (u(i1,i2+2,i3+2,kd)-u(i1,i2-2,i3+2,kd)-u(i1,i2+2,i3-2,kd)+u(i1,i2-2,i3-2,kd)) +8.*(u(i1,i2-1,i3+2,kd)-u(i1,i2-1,i3-2,kd)-u(i1,i2+1,i3+2,kd)+u(i1,i2+1,i3-2,kd) +u(i1,i2+2,i3-1,kd)-u(i1,i2-2,i3-1,kd)-u(i1,i2+2,i3+1,kd)+u(i1,i2-2,i3+1,kd)) +64.*(u(i1,i2+1,i3+1,kd)-u(i1,i2-1,i3+1,kd)-u(i1,i2+1,i3-1,kd)+u(i1,i2-1,i3-1,kd)) )*(dx41(1)*dx41(2))
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
         ! 4th-order 1 sided derivative  extrap=(1 5 10 10 5 1)
         uxOneSided(i1,i2,i3,m)=-(10./3.)*u(i1,i2,i3,m)+6.*u(i1+is1,i2+is2,i3+is3,m)-2.*u(i1+2*is1,i2+2*is2,i3+2*is3,m)+(1./3.)*u(i1+3*is1,i2+3*is2,i3+3*is3,m)
         ! Here is the the generic boundary condition forcing array. It uses the bcOffset(side,axis) values as an
         ! an offset from the bcf0 array to access the bcf10, bcf01, bcf11, ... arrays
         bcf(side,axis,i1,i2,i3,m) = bcf0(bcOffset(side,axis) + (i1-dim(0,0,side,axis)+(dim(1,0,side,axis)-dim(0,0,side,axis)+1)* (i2-dim(0,1,side,axis)+(dim(1,1,side,axis)-dim(0,1,side,axis)+1)* (i3-dim(0,2,side,axis)+(dim(1,2,side,axis)-dim(0,2,side,axis)+1)*(m)))))
         !............... end statement functions
         ierr=0
         nd                           = ipar(0)
         grid                         = ipar(1)
         u1c                          = ipar(2)
         u2c                          = ipar(3)
         u3c                          = ipar(4)
         gridType                     = ipar(5)
         orderOfAccuracy              = ipar(6)
         orderOfExtrapolation         = ipar(7)
         twilightZone                 = ipar(8)
         useWhereMask                 = ipar(9)
         debug                        = ipar(10)
         materialFormat               = ipar(15)
         pc                           = ipar(16)
         upwindSOS                    = ipar(17)
         useCurlCurlBoundaryCondition = ipar(18);
         dx(0)                = rpar(0)
         dx(1)                = rpar(1)
         dx(2)                = rpar(2)
         dr(0)                = rpar(3)
         dr(1)                = rpar(4)
         dr(2)                = rpar(5)
         t                    = rpar(6)
         ep                   = rpar(7) ! pointer for exact solution
         dt                   = rpar(8)
         rho                  = rpar(9)
         mu                   = rpar(10)
         lambda               = rpar(11)
         c1                   = rpar(12)
         c2                   = rpar(13)
         job=0  ! *wdh* 090101
         numGhost=orderOfAccuracy/2
         if( upwindSOS==1 )then
           ! upwinding uses an extra ghost line
           numGhost=numGhost+1
         end if 
         ! check residuals for the first few steps
         checkResiduals=.false.
         if( t <=2.*dt )then
           checkResiduals=.true.
         end if
         if( t<=2*dt .or. debug.gt.3 )then
           write(*,'("# bcOptIsm: t=",1pe12.3)') t
           write(*,'(" bcOptIsm: rho,mu=",2f10.5," gridType=",i2," upwindSOS=",i2," numGhost=",i2)') rho,mu,gridType,upwindSOS,numGhost
           write(*,'(" bcOptIsm: u1c,u2c,u3c,pc=",4i3," twilightZone=",i2," orderOfAccuracy=",i2)') u1c,u2c,u3c,pc,twilightZone,orderOfAccuracy
           write(*,'(" boundaryCondition=",6i4)') ((boundaryCondition(side,axis),side=0,1),axis=0,nd-1)
           write(*,'(" addBoundaryForcing=",6i4)') ((addBoundaryForcing(side,axis),side=0,1),axis=0,nd-1)
         end if
         if( debug.gt.7 )then
           write(*,'(" bcOptIsm: **START** grid=",i4," u1c,u2c,u3c=",3i2)') grid,u1c,u2c,u3c
              ! '
         end if
         if( debug.gt.7 )then
          n1a=gridIndexRange(0,0)
          n1b=gridIndexRange(1,0)
          n2a=gridIndexRange(0,1)
          n2b=gridIndexRange(1,1)
          n3a=gridIndexRange(0,2)
          n3b=gridIndexRange(1,2)
          write(*,'(" bcOptIsm: grid=",i3,",n1a,n1b,n2a,n2b,n3a,n3b=",6i3)') grid,n1a,n1b,n2a,n2b,n3a,n3b
          ! write(*,*) 'bcOptIsm: u=',((((u(i1,i2,i3,m),m=0,nd-1),i1=n1a,n1b),i2=n2a,n2b),i3=n3a,n3b)
         end if
         if( materialFormat.ne.constantMaterialProperties )then
           write(*,'(" ***bcOptIsm:ERROR: Finish me for variable material")')
           stop 7736
         end if
         epsx=1.e-20  ! fix this 
         resTol = 1.e-9
         ! write(*,*) 'bcOffset:',bcOffset
         ! ! assign corners and edges (3d)
         ! if( .false. )then ! *wdh* 071027 -- turn this off for now ---
         !   if( orderOfAccuracy.eq.2 .and. nd.eq.2 )then
         !     ! *** fix this for traction BCs with forcing 
         !    ! For interfaces and TZ it is ok to just set the corners using TZ. This is maybe cheating a bit.
         !    if( gridType.eq.rectangular )then
         !     if( twilightZone.eq.0 )then
         !       assignCorners2d(2,rectangular,none)
         !     else
         !       assignCorners2d(2,rectangular,twilightZone)
         !     end if
         !    else
         !     if( twilightZone.eq.0 )then
         !       assignCorners2d(2,curvilinear,none)
         !     else
         !       assignCorners2d(2,curvilinear,twilightZone)
         !     end if
         !    end if      
         !   else if( orderOfAccuracy.eq.2 .and. nd.eq.3 )then
         !     !$$$       if( gridType.eq.rectangular )then
         !     !$$$        if( twilightZone.eq.0 )then
         !     !$$$          assignCorners3d(2,rectangular,none)
         !     !$$$        else
         !     !$$$          assignCorners2d(2,rectangular,twilightZone)
         !     !$$$        end if
         !     !$$$       else
         !     !$$$        if( twilightZone.eq.0 )then
         !     !$$$          assignCorners3d(2,curvilinear,none)
         !     !$$$        else
         !     !$$$          assignCorners3d(2,curvilinear,twilightZone)
         !     !$$$        end if
         !     !$$$       end if      
         !   else
         !      stop 5533
         !   end if
         ! end if
         ! if( .false. )then
         !   ! check the boundary forcing arrays: check that bcf(side,axis,i1,i2,i3,m) agrees with bcf00, bcf10, ...
         !   write(*,*) dim
         !   beginLoopOverSides(numGhost,numGhost,none)
         !    write(*,'(" BCF: side,axis=",2i3," bcOffset(side,axis)=",i8)') side,axis,bcOffset(side,axis)
         !    if( addBoundaryForcing(side,axis).ne.0 )then
         !      beginLoops3d()
         !        if( side.eq.0 .and. axis.eq.0 )then
         !          do m=0,nd-1
         !            tmp = bcf00(i1,i2,i3,m) - bcf(side,axis,i1,i2,i3,m)
         !            write(*,'(" BCF(0,0): i=",3i3," f=",e8.2,2x,e8.2," diff=",e8.2)') i1,i2,i3,bcf00(i1,i2,i3,m),bcf(side,axis,i1,i2,i3,m),tmp
         !            ! '
         !          end do
         !        else if( side.eq.1 .and. axis.eq.0 )then
         !          do m=0,nd-1
         !            tmp = bcf10(i1,i2,i3,m) - bcf(side,axis,i1,i2,i3,m)
         !            write(*,'(" BCF(1,0): i=",3i3," f=",e8.2,2x,e8.2," diff=",e8.2)') i1,i2,i3,bcf10(i1,i2,i3,m),bcf(side,axis,i1,i2,i3,m),tmp
         !            ! '
         !          end do
         !        else if( side.eq.0 .and. axis.eq.1 )then
         !          do m=0,nd-1
         !            tmp = bcf01(i1,i2,i3,m) - bcf(side,axis,i1,i2,i3,m)
         !            write(*,'(" BCF(0,1): i=",3i3," f=",e8.2,2x,e8.2," diff=",e8.2)') i1,i2,i3,bcf01(i1,i2,i3,m),bcf(side,axis,i1,i2,i3,m),tmp
         !            ! '
         !          end do
         !        else if( side.eq.1 .and. axis.eq.1 )then
         !          do m=0,nd-1
         !            tmp = bcf11(i1,i2,i3,m) - bcf(side,axis,i1,i2,i3,m)
         !            write(*,'(" BCF(1,1): i=",3i3," f=",e8.2,2x,e8.2," diff=",e8.2)') i1,i2,i3,bcf11(i1,i2,i3,m),bcf(side,axis,i1,i2,i3,m),tmp
         !            ! '
         !          end do
         !        else if( side.eq.0 .and. axis.eq.2 )then
         !          do m=0,nd-1
         !            tmp = bcf02(i1,i2,i3,m) - bcf(side,axis,i1,i2,i3,m)
         !            write(*,'(" BCF(0,2): i=",3i3," f=",e8.2,2x,e8.2," diff=",e8.2)') i1,i2,i3,bcf02(i1,i2,i3,m),bcf(side,axis,i1,i2,i3,m),tmp
         !            ! '
         !          end do
         !        else if( side.eq.1 .and. axis.eq.2 )then
         !          do m=0,nd-1
         !            tmp = bcf12(i1,i2,i3,m) - bcf(side,axis,i1,i2,i3,m)
         !            write(*,'(" BCF(1,2): i=",3i3," f=",e8.2,2x,e8.2," diff=",e8.2)') i1,i2,i3,bcf12(i1,i2,i3,m),bcf(side,axis,i1,i2,i3,m),tmp
         !            ! '
         !          end do
         !        end if
         !      endLoops3d()
         !    end if
         !   endLoopOverSides()
         ! end if
         ! --- NEW WAY -----
         ! ---------------------------------------------------------------------
         ! ----------- STAGE 1, ASSIGN DIRICHLET TYPE CONDITIONS ---------------
         ! --------------------------------------------------------------------- 
            extra1a=numGhost
            extra1b=numGhost
            extra2a=numGhost
            extra2b=numGhost
            if( nd.eq.3 )then
              extra3a=numGhost
              extra3b=numGhost
            else
              extra3a=0
              extra3b=0
            end if
            if( boundaryCondition(0,0).lt.0 )then
              extra1a=max(0,extra1a) ! over-ride numGhost=-1 : assign ends in periodic directions (or internal parallel boundaries)
            else if( boundaryCondition(0,0).eq.0 )then
              extra1a=numGhost  ! include interpolation points since we assign ghost points outside these
            end if
            ! **NOTE** the bc on the right may be negative even it is not on the left (for parallel)
            if( boundaryCondition(1,0).lt.0 )then
              extra1b=max(0,extra1b) ! over-ride numGhost=-1 : assign ends in periodic directions
            else if( boundaryCondition(1,0).eq.0 )then
              extra1b=numGhost
            end if
            if( boundaryCondition(0,1).lt.0 )then
              extra2a=max(0,extra2a) ! over-ride numGhost=-1 : assign ends in periodic directions (or internal parallel boundaries)
            else if( boundaryCondition(0,1).eq.0 )then
              extra2a=numGhost  ! include interpolation points since we assign ghost points outside these
            end if
            ! **NOTE** the bc on the right may be negative even it is not on the left (for parallel)
            if( boundaryCondition(1,1).lt.0 )then
              extra2b=max(0,extra2b) ! over-ride numGhost=-1 : assign ends in periodic directions
            else if( boundaryCondition(1,1).eq.0 )then
              extra2b=numGhost
            end if
            if(  nd.eq.3 )then
             if( boundaryCondition(0,2).lt.0 )then
               extra3a=max(0,extra3a) ! over-ride numGhost=-1 : assign ends in periodic directions (or internal parallel boundaries)
             else if( boundaryCondition(0,2).eq.0 )then
               extra3a=numGhost  ! include interpolation points since we assign ghost points outside these
             end if
             ! **NOTE** the bc on the right may be negative even it is not on the left (for parallel)
             if( boundaryCondition(1,2).lt.0 )then
               extra3b=max(0,extra3b) ! over-ride numGhost=-1 : assign ends in periodic directions
             else if( boundaryCondition(1,2).eq.0 )then
               extra3b=numGhost
             end if
            end if
            do axis=0,nd-1
            do side=0,1
              if( boundaryCondition(side,axis).gt.0 )then
                ! write(*,'(" bcOpt: side,axis,bc=",3i2)') side,axis,boundaryCondition(side,axis)
                ! get gid(0;1,0:2) -- adjusted for priority at corners
                 do leftRight=0,1
                   do dir=0,2
                     gid(leftRight,dir) = gridIndexRange(leftRight,dir)
                   end do
                 end do       
                axisp1=mod(axis+1,nd)
                axisp2=mod(axis+2,nd)
                i3=n3a
                n1a=gid(0,0)
                n1b=gid(1,0)
                n2a=gid(0,1)
                n2b=gid(1,1)
                n3a=gid(0,2)
                n3b=gid(1,2)
                if( axis.eq.0 )then
                  n1a=gid(side,axis)
                  n1b=gid(side,axis)
                else if( axis.eq.1 )then
                  n2a=gid(side,axis)
                  n2b=gid(side,axis)
                else
                  n3a=gid(side,axis)
                  n3b=gid(side,axis)
                end if
                ! if( boundaryCondition(side,axis)==tractionBC )then
                !   if( boundaryCondition(axisp1,0)==dirichletBoundaryCondition )then
                !     n2a = n2a+1; ! skip boundary
                !   end if
                !   if( boundaryCondition(axisp1,1)==dirichletBoundaryCondition )then
                !     n2b = n2b-1; ! skip boundary
                !   end if 
                ! end if
                nn1a=gid(0,0)-extra1a
                nn1b=gid(1,0)+extra1b
                nn2a=gid(0,1)-extra2a
                nn2b=gid(1,1)+extra2b
                nn3a=gid(0,2)-extra3a
                nn3b=gid(1,2)+extra3b
                if( axis.eq.0 )then
                  nn1a=gid(side,axis)
                  nn1b=gid(side,axis)
                else if( axis.eq.1 )then
                  nn2a=gid(side,axis)
                  nn2b=gid(side,axis)
                else
                  nn3a=gid(side,axis)
                  nn3b=gid(side,axis)
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
                if( debug.gt.7 )then
                  write(*,'(" bcOpt: grid,side,axis=",3i3,",loop bounds: n1a,n1b,n2a,n2b,n3a,n3b=",6i3)') grid,side,axis,n1a,n1b,n2a,n2b,n3a,n3b
                  write(*,'("                                                       gridIndexRange=",6i3)') ((gridIndexRange(m,dir),m=0,1),dir=0,nd-1)
                  write(*,'("                                                                  gid=",6i3)') ((gid(m,dir),m=0,1),dir=0,nd-1)
                end if
              end if ! if bc>0 
              ! On interfaces we should use the bcf array values even for TZ since then
              ! we get a coupling at the interface: 
              !   bcf = n.sigma(fluid) + [ n.sigma_e(solid) - n.sigma_e(fluid) ]
              if( interfaceType(side,axis,grid).eq.noInterface )then
                assignTwilightZone=twilightZone
              else
                assignTwilightZone=0  ! this will turn off the use of TZ
              end if
             if( boundaryCondition(side,axis).eq.displacementBC )then
               ! Note: addBoundaryForcing =1 for TZ as well
               if( twilightZone==1 )then
                  if( .false. )then
                    write(*,*) 'Displacement BC -- assign extended boundaries:'
                    write(*,'(" bcOpt: grid,side,axis=",3i3,",loop bounds: nn1a,nn1b,nn2a,nn2b,nn3a,nn3b=",6i3)') grid,side,axis,nn1a,nn1b,nn2a,nn2b,nn3a,nn3b
                    write(*,'("                                                           gridIndexRange=",6i3)') ((gridIndexRange(m,dir),m=0,1),dir=0,nd-1)
                    write(*,'("                                                                      gid=",6i3)') ((gid(m,dir),m=0,1),dir=0,nd-1)
                  end if
                   ! ------ Assigned extended boundaries with dirichlet values ------
                    do i3=nn3a,nn3b
                    do i2=nn2a,nn2b
                    do i1=nn1a,nn1b
                          call ogf3d(ep,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,u0,v0,w0)
                         u(i1,i2,i3,u1c)=u0
                         u(i1,i2,i3,u2c)=v0
                         u(i1,i2,i3,u3c)=w0       
                    end do
                    end do
                    end do
               else if( .true. .or. addBoundaryForcing(side,axis)==0 )then  ! FOR NOW SKIP addBoundaryFORCING
                  if( .false. )then
                    write(*,*) 'Displacement BC -- assign extended boundaries:'
                    write(*,'(" bcOpt: grid,side,axis=",3i3,",loop bounds: nn1a,nn1b,nn2a,nn2b,nn3a,nn3b=",6i3)') grid,side,axis,nn1a,nn1b,nn2a,nn2b,nn3a,nn3b
                    write(*,'("                                                           gridIndexRange=",6i3)') ((gridIndexRange(m,dir),m=0,1),dir=0,nd-1)
                    write(*,'("                                                                      gid=",6i3)') ((gid(m,dir),m=0,1),dir=0,nd-1)
                  end if
                   ! ------ Assigned extended boundaries with dirichlet values ------
                    do i3=nn3a,nn3b
                    do i2=nn2a,nn2b
                    do i1=nn1a,nn1b
                       u(i1,i2,i3,u1c)=0.
                       u(i1,i2,i3,u2c)=0.
                         u(i1,i2,i3,u3c)=0.
                    end do
                    end do
                    end do
               else 
                  if( .false. )then
                    write(*,*) 'Displacement BC -- assign extended boundaries:'
                    write(*,'(" bcOpt: grid,side,axis=",3i3,",loop bounds: nn1a,nn1b,nn2a,nn2b,nn3a,nn3b=",6i3)') grid,side,axis,nn1a,nn1b,nn2a,nn2b,nn3a,nn3b
                    write(*,'("                                                           gridIndexRange=",6i3)') ((gridIndexRange(m,dir),m=0,1),dir=0,nd-1)
                    write(*,'("                                                                      gid=",6i3)') ((gid(m,dir),m=0,1),dir=0,nd-1)
                  end if
                   ! ------ Assigned extended boundaries with dirichlet values ------
                    do i3=nn3a,nn3b
                    do i2=nn2a,nn2b
                    do i1=nn1a,nn1b
                       ! Use forcing in array bcf(..)  
                       stop 1234 
                    end do
                    end do
                    end do
               end if
            else if( boundaryCondition(side,axis).eq.tractionBC )then 
              ! Extrapolation now done above
            else if( boundaryCondition(side,axis).eq.slipWall )then
              ! set n.u = given on the boundary 
              an1=0.  ! (an1,an2) = outward normal 
              an2=0.
              if( axis.eq.0 )then
                an1=2*side-1
              else
                an2=2*side-1
              end if
              if( addBoundaryForcing(side,axis).eq.0 )then
                ! no forcing 
                 i3=n3a
                 do i2=n2a,n2b
                 do i1=n1a,n1b
                 u1 = u(i1,i2,i3,u1c)
                 v1 = u(i1,i2,i3,u2c)
                 nDotU = an1*u1 + an2*v1  
                 u(i1,i2,i3,u1c)=u1 - nDotU*an1
                 u(i1,i2,i3,u2c)=v1 - nDotU*an2
                 end do
                 end do
               else if( assignTwilightZone.eq.0 )then
                 ! include forcing terms 
                 ! n.u = n.g 
                 i3=n3a
                 do i2=n2a,n2b
                 do i1=n1a,n1b
                 u1 = u(i1,i2,i3,u1c)
                 v1 = u(i1,i2,i3,u2c)
                 nDotU = an1*(u1-bcf(side,axis,i1,i2,i3,u1c)) + an2*(v1-bcf(side,axis,i1,i2,i3,u2c))
                 u(i1,i2,i3,u1c)=u1 - nDotU*an1
                 u(i1,i2,i3,u2c)=v1 - nDotU*an2
                 end do
                 end do
               else
                ! Twilight-zone: 
                !   n.u = n.ue
                 i3=n3a
                 do i2=n2a,n2b
                 do i1=n1a,n1b
                  call ogf2d(ep,xy(i1,i2,i3,0),xy(i1,i2,i3,1),t,u0,v0)
                 u1 = u(i1,i2,i3,u1c)
                 v1 = u(i1,i2,i3,u2c)
                 nDotU = an1*(u1-u0) + an2*(v1-v0)
                 u(i1,i2,i3,u1c)=u1 - nDotU*an1
                 u(i1,i2,i3,u2c)=v1 - nDotU*an2
                 end do
                 end do
               end if
               ! extrap values to the ghost line 
                i3=n3a
                do i2=nn2a,nn2b
                do i1=nn1a,nn1b
                u(i1-is1,i2-is2,i3,u1c)=(3.*u(i1,i2,i3,u1c)-3.*u(i1+is1,i2+is2,i3+is3,u1c)+u(i1+2*is1,i2+2*is2,i3+2*is3,u1c))
                u(i1-is1,i2-is2,i3,u2c)=(3.*u(i1,i2,i3,u2c)-3.*u(i1+is1,i2+is2,i3+is3,u2c)+u(i1+2*is1,i2+2*is2,i3+2*is3,u2c))
                end do
                end do
            else if( boundaryCondition(side,axis).eq.dirichletBoundaryCondition )then
               ! Set exact values on boundary and ghost
               if( twilightZone.eq.1 )then
                  do i3=nn3a,nn3b
                  do i2=nn2a,nn2b
                  do i1=nn1a,nn1b
                      call ogf3d(ep,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,u0,v0,w0)
                     u(i1,i2,i3,u1c)=u0
                     u(i1,i2,i3,u2c)=v0
                     u(i1,i2,i3,u3c)=w0
                   do ghost=1,numGhost
                     j1 = i1 - ghost*is1
                     j2 = i2 - ghost*is2
                       j3 = i3 - ghost*is3
                        call ogf3d(ep,xy(j1,j2,j3,0),xy(j1,j2,j3,1),xy(j1,j2,j3,2),t,u0,v0,w0)
                       u(j1,j2,j3,u1c)=u0
                       u(j1,j2,j3,u2c)=v0              
                       u(j1,j2,j3,u3c)=w0              
                   end do
                  end do
                  end do
                  end do
               end if
            else if( boundaryCondition(side,axis).eq.symmetry )then
               write(*,*) 'finish me for symmetry BC'
              stop 2856
            else if( boundaryCondition(side,axis).gt.0 )then
              stop 1193
            end if
            end do ! end side
            end do ! end axis
         ! ---------------------------------------------------------------
         ! Stage 1(b) extrap values to ghost points (may be needed at corners)
         ! ---------------------------------------------------------------
          extra1a=numGhost
          extra1b=numGhost
          extra2a=numGhost
          extra2b=numGhost
          if( nd.eq.3 )then
            extra3a=numGhost
            extra3b=numGhost
          else
            extra3a=0
            extra3b=0
          end if
          if( boundaryCondition(0,0).lt.0 )then
            extra1a=max(0,extra1a) ! over-ride numGhost=-1 : assign ends in periodic directions (or internal parallel boundaries)
          else if( boundaryCondition(0,0).eq.0 )then
            extra1a=numGhost  ! include interpolation points since we assign ghost points outside these
          end if
          ! **NOTE** the bc on the right may be negative even it is not on the left (for parallel)
          if( boundaryCondition(1,0).lt.0 )then
            extra1b=max(0,extra1b) ! over-ride numGhost=-1 : assign ends in periodic directions
          else if( boundaryCondition(1,0).eq.0 )then
            extra1b=numGhost
          end if
          if( boundaryCondition(0,1).lt.0 )then
            extra2a=max(0,extra2a) ! over-ride numGhost=-1 : assign ends in periodic directions (or internal parallel boundaries)
          else if( boundaryCondition(0,1).eq.0 )then
            extra2a=numGhost  ! include interpolation points since we assign ghost points outside these
          end if
          ! **NOTE** the bc on the right may be negative even it is not on the left (for parallel)
          if( boundaryCondition(1,1).lt.0 )then
            extra2b=max(0,extra2b) ! over-ride numGhost=-1 : assign ends in periodic directions
          else if( boundaryCondition(1,1).eq.0 )then
            extra2b=numGhost
          end if
          if(  nd.eq.3 )then
           if( boundaryCondition(0,2).lt.0 )then
             extra3a=max(0,extra3a) ! over-ride numGhost=-1 : assign ends in periodic directions (or internal parallel boundaries)
           else if( boundaryCondition(0,2).eq.0 )then
             extra3a=numGhost  ! include interpolation points since we assign ghost points outside these
           end if
           ! **NOTE** the bc on the right may be negative even it is not on the left (for parallel)
           if( boundaryCondition(1,2).lt.0 )then
             extra3b=max(0,extra3b) ! over-ride numGhost=-1 : assign ends in periodic directions
           else if( boundaryCondition(1,2).eq.0 )then
             extra3b=numGhost
           end if
          end if
          do axis=0,nd-1
          do side=0,1
            if( boundaryCondition(side,axis).gt.0 )then
              ! write(*,'(" bcOpt: side,axis,bc=",3i2)') side,axis,boundaryCondition(side,axis)
              ! get gid(0;1,0:2) -- adjusted for priority at corners
                  do leftRight=0,1
                    do dir=0,2
                      gid(leftRight,dir) = gridIndexRange(leftRight,dir)
                    end do
                  end do
                  bcSide = boundaryCondition(side,axis)
                  if( bcSide.eq.tractionBC .or. bcSide.eq.displacementBC )then
                    ! adjust gridIndexRange at ends of the boundary 
                    !  Note: traction-traction corner assigns the extended ghost on the boundary
                    do dir=1,nd-1
                      axisp = mod(axis+dir,nd)  ! tangential direction
                      do leftRight=0,1
                        bcAdjacent = boundaryCondition(leftRight,axisp)
                        adjustEnds=.false.
                        if( bcSide.eq.displacementBC .and. ( bcAdjacent.eq.displacementBC .or. bcAdjacent.eq.dirichletBoundaryCondition ) )then
                          ! extended boundary has been set, no need to fill ghost on ends
                          adjustEnds=.true.
                        end if 
                        if( bcSide.eq.tractionBC .and. ( bcAdjacent.eq.displacementBC .or. bcAdjacent.eq.dirichletBoundaryCondition ) )then  
                          ! extended boundary has been set, no need to fill ghost on ends
                          adjustEnds=.true.
                        end if
                        ! if( boundaryCondition(leftRight,axisp)==dirichletBoundaryCondition .or. !     boundaryCondition(leftRight,axisp)==displacementBC             .or. !     boundaryCondition(leftRight,axisp)==tractionBC      )then
                        if( adjustEnds )then
                          ! shift to avoid the corner
                          gid(leftRight,axisp) = gid(leftRight,axisp) + (1-2*leftRight)
                        end if
                      end do
                    end do
                  end if
              axisp1=mod(axis+1,nd)
              axisp2=mod(axis+2,nd)
              i3=n3a
              n1a=gid(0,0)
              n1b=gid(1,0)
              n2a=gid(0,1)
              n2b=gid(1,1)
              n3a=gid(0,2)
              n3b=gid(1,2)
              if( axis.eq.0 )then
                n1a=gid(side,axis)
                n1b=gid(side,axis)
              else if( axis.eq.1 )then
                n2a=gid(side,axis)
                n2b=gid(side,axis)
              else
                n3a=gid(side,axis)
                n3b=gid(side,axis)
              end if
              ! if( boundaryCondition(side,axis)==tractionBC )then
              !   if( boundaryCondition(axisp1,0)==dirichletBoundaryCondition )then
              !     n2a = n2a+1; ! skip boundary
              !   end if
              !   if( boundaryCondition(axisp1,1)==dirichletBoundaryCondition )then
              !     n2b = n2b-1; ! skip boundary
              !   end if 
              ! end if
              nn1a=gid(0,0)-extra1a
              nn1b=gid(1,0)+extra1b
              nn2a=gid(0,1)-extra2a
              nn2b=gid(1,1)+extra2b
              nn3a=gid(0,2)-extra3a
              nn3b=gid(1,2)+extra3b
              if( axis.eq.0 )then
                nn1a=gid(side,axis)
                nn1b=gid(side,axis)
              else if( axis.eq.1 )then
                nn2a=gid(side,axis)
                nn2b=gid(side,axis)
              else
                nn3a=gid(side,axis)
                nn3b=gid(side,axis)
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
              if( debug.gt.7 )then
                write(*,'(" bcOpt: grid,side,axis=",3i3,",loop bounds: n1a,n1b,n2a,n2b,n3a,n3b=",6i3)') grid,side,axis,n1a,n1b,n2a,n2b,n3a,n3b
                write(*,'("                                                       gridIndexRange=",6i3)') ((gridIndexRange(m,dir),m=0,1),dir=0,nd-1)
                write(*,'("                                                                  gid=",6i3)') ((gid(m,dir),m=0,1),dir=0,nd-1)
              end if
            end if ! if bc>0 
            ! On interfaces we should use the bcf array values even for TZ since then
            ! we get a coupling at the interface: 
            !   bcf = n.sigma(fluid) + [ n.sigma_e(solid) - n.sigma_e(fluid) ]
            if( interfaceType(side,axis,grid).eq.noInterface )then
              assignTwilightZone=twilightZone
            else
              assignTwilightZone=0  ! this will turn off the use of TZ
            end if
           if( boundaryCondition(side,axis)==displacementBC .or. boundaryCondition(side,axis)==tractionBC )then
             ! write(*,'("bcOptIsm: extrap ghost on displacementBC or tractionBC ")')
               !  if( boundaryCondition(side,axis).eq.tractionBC )then
               !     ! Set exact values on boundary and ghost
               !     if( twilightZone.eq.1 )then
               !       beginLoops2d()
               !         OGF2D(i1,i2,i3,t,u0,v0)
               !         write(*,'(" traction: boundary i1,i2=",2i3," error=",2(1pe10.2))') i1,i2,u(i1,i2,i3,u1c)-u0,u(i1,i2,i3,u2c)-v0
               !       endLoops2d()
               !     end if        
               ! end if
                do i3=nn3a,nn3b
                do i2=nn2a,nn2b
                do i1=nn1a,nn1b
               ! beginLoops3d()
                 do ghost=1,numGhost
                   j1 = i1 - ghost*is1
                   j2 = i2 - ghost*is2
                   j3 = i3 - ghost*is3        
                   u(j1,j2,j3,u1c)=(3.*u(j1+is1,j2+is2,j3+is3,u1c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u1c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u1c))
                   u(j1,j2,j3,u2c)=(3.*u(j1+is1,j2+is2,j3+is3,u2c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u2c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u2c))
                   u(j1,j2,j3,u3c)=(3.*u(j1+is1,j2+is2,j3+is3,u3c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u3c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u3c))
                 end do
                end do
                end do
                end do
           end if
          end do ! end side
          end do ! end axis
         ! assign values on boundary again to get extended boundaries -- *fix me*
            extra1a=numGhost
            extra1b=numGhost
            extra2a=numGhost
            extra2b=numGhost
            if( nd.eq.3 )then
              extra3a=numGhost
              extra3b=numGhost
            else
              extra3a=0
              extra3b=0
            end if
            if( boundaryCondition(0,0).lt.0 )then
              extra1a=max(0,extra1a) ! over-ride numGhost=-1 : assign ends in periodic directions (or internal parallel boundaries)
            else if( boundaryCondition(0,0).eq.0 )then
              extra1a=numGhost  ! include interpolation points since we assign ghost points outside these
            end if
            ! **NOTE** the bc on the right may be negative even it is not on the left (for parallel)
            if( boundaryCondition(1,0).lt.0 )then
              extra1b=max(0,extra1b) ! over-ride numGhost=-1 : assign ends in periodic directions
            else if( boundaryCondition(1,0).eq.0 )then
              extra1b=numGhost
            end if
            if( boundaryCondition(0,1).lt.0 )then
              extra2a=max(0,extra2a) ! over-ride numGhost=-1 : assign ends in periodic directions (or internal parallel boundaries)
            else if( boundaryCondition(0,1).eq.0 )then
              extra2a=numGhost  ! include interpolation points since we assign ghost points outside these
            end if
            ! **NOTE** the bc on the right may be negative even it is not on the left (for parallel)
            if( boundaryCondition(1,1).lt.0 )then
              extra2b=max(0,extra2b) ! over-ride numGhost=-1 : assign ends in periodic directions
            else if( boundaryCondition(1,1).eq.0 )then
              extra2b=numGhost
            end if
            if(  nd.eq.3 )then
             if( boundaryCondition(0,2).lt.0 )then
               extra3a=max(0,extra3a) ! over-ride numGhost=-1 : assign ends in periodic directions (or internal parallel boundaries)
             else if( boundaryCondition(0,2).eq.0 )then
               extra3a=numGhost  ! include interpolation points since we assign ghost points outside these
             end if
             ! **NOTE** the bc on the right may be negative even it is not on the left (for parallel)
             if( boundaryCondition(1,2).lt.0 )then
               extra3b=max(0,extra3b) ! over-ride numGhost=-1 : assign ends in periodic directions
             else if( boundaryCondition(1,2).eq.0 )then
               extra3b=numGhost
             end if
            end if
            do axis=0,nd-1
            do side=0,1
              if( boundaryCondition(side,axis).gt.0 )then
                ! write(*,'(" bcOpt: side,axis,bc=",3i2)') side,axis,boundaryCondition(side,axis)
                ! get gid(0;1,0:2) -- adjusted for priority at corners
                 do leftRight=0,1
                   do dir=0,2
                     gid(leftRight,dir) = gridIndexRange(leftRight,dir)
                   end do
                 end do       
                axisp1=mod(axis+1,nd)
                axisp2=mod(axis+2,nd)
                i3=n3a
                n1a=gid(0,0)
                n1b=gid(1,0)
                n2a=gid(0,1)
                n2b=gid(1,1)
                n3a=gid(0,2)
                n3b=gid(1,2)
                if( axis.eq.0 )then
                  n1a=gid(side,axis)
                  n1b=gid(side,axis)
                else if( axis.eq.1 )then
                  n2a=gid(side,axis)
                  n2b=gid(side,axis)
                else
                  n3a=gid(side,axis)
                  n3b=gid(side,axis)
                end if
                ! if( boundaryCondition(side,axis)==tractionBC )then
                !   if( boundaryCondition(axisp1,0)==dirichletBoundaryCondition )then
                !     n2a = n2a+1; ! skip boundary
                !   end if
                !   if( boundaryCondition(axisp1,1)==dirichletBoundaryCondition )then
                !     n2b = n2b-1; ! skip boundary
                !   end if 
                ! end if
                nn1a=gid(0,0)-extra1a
                nn1b=gid(1,0)+extra1b
                nn2a=gid(0,1)-extra2a
                nn2b=gid(1,1)+extra2b
                nn3a=gid(0,2)-extra3a
                nn3b=gid(1,2)+extra3b
                if( axis.eq.0 )then
                  nn1a=gid(side,axis)
                  nn1b=gid(side,axis)
                else if( axis.eq.1 )then
                  nn2a=gid(side,axis)
                  nn2b=gid(side,axis)
                else
                  nn3a=gid(side,axis)
                  nn3b=gid(side,axis)
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
                if( debug.gt.7 )then
                  write(*,'(" bcOpt: grid,side,axis=",3i3,",loop bounds: n1a,n1b,n2a,n2b,n3a,n3b=",6i3)') grid,side,axis,n1a,n1b,n2a,n2b,n3a,n3b
                  write(*,'("                                                       gridIndexRange=",6i3)') ((gridIndexRange(m,dir),m=0,1),dir=0,nd-1)
                  write(*,'("                                                                  gid=",6i3)') ((gid(m,dir),m=0,1),dir=0,nd-1)
                end if
              end if ! if bc>0 
              ! On interfaces we should use the bcf array values even for TZ since then
              ! we get a coupling at the interface: 
              !   bcf = n.sigma(fluid) + [ n.sigma_e(solid) - n.sigma_e(fluid) ]
              if( interfaceType(side,axis,grid).eq.noInterface )then
                assignTwilightZone=twilightZone
              else
                assignTwilightZone=0  ! this will turn off the use of TZ
              end if
             if( boundaryCondition(side,axis).eq.displacementBC )then
               ! Note: addBoundaryForcing =1 for TZ as well
               if( twilightZone==1 )then
                  if( .false. )then
                    write(*,*) 'Displacement BC -- assign extended boundaries:'
                    write(*,'(" bcOpt: grid,side,axis=",3i3,",loop bounds: nn1a,nn1b,nn2a,nn2b,nn3a,nn3b=",6i3)') grid,side,axis,nn1a,nn1b,nn2a,nn2b,nn3a,nn3b
                    write(*,'("                                                           gridIndexRange=",6i3)') ((gridIndexRange(m,dir),m=0,1),dir=0,nd-1)
                    write(*,'("                                                                      gid=",6i3)') ((gid(m,dir),m=0,1),dir=0,nd-1)
                  end if
                   ! ------ Assigned extended boundaries with dirichlet values ------
                    do i3=nn3a,nn3b
                    do i2=nn2a,nn2b
                    do i1=nn1a,nn1b
                          call ogf3d(ep,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,u0,v0,w0)
                         u(i1,i2,i3,u1c)=u0
                         u(i1,i2,i3,u2c)=v0
                         u(i1,i2,i3,u3c)=w0       
                    end do
                    end do
                    end do
               else if( .true. .or. addBoundaryForcing(side,axis)==0 )then  ! FOR NOW SKIP addBoundaryFORCING
                  if( .false. )then
                    write(*,*) 'Displacement BC -- assign extended boundaries:'
                    write(*,'(" bcOpt: grid,side,axis=",3i3,",loop bounds: nn1a,nn1b,nn2a,nn2b,nn3a,nn3b=",6i3)') grid,side,axis,nn1a,nn1b,nn2a,nn2b,nn3a,nn3b
                    write(*,'("                                                           gridIndexRange=",6i3)') ((gridIndexRange(m,dir),m=0,1),dir=0,nd-1)
                    write(*,'("                                                                      gid=",6i3)') ((gid(m,dir),m=0,1),dir=0,nd-1)
                  end if
                   ! ------ Assigned extended boundaries with dirichlet values ------
                    do i3=nn3a,nn3b
                    do i2=nn2a,nn2b
                    do i1=nn1a,nn1b
                       u(i1,i2,i3,u1c)=0.
                       u(i1,i2,i3,u2c)=0.
                         u(i1,i2,i3,u3c)=0.
                    end do
                    end do
                    end do
               else 
                  if( .false. )then
                    write(*,*) 'Displacement BC -- assign extended boundaries:'
                    write(*,'(" bcOpt: grid,side,axis=",3i3,",loop bounds: nn1a,nn1b,nn2a,nn2b,nn3a,nn3b=",6i3)') grid,side,axis,nn1a,nn1b,nn2a,nn2b,nn3a,nn3b
                    write(*,'("                                                           gridIndexRange=",6i3)') ((gridIndexRange(m,dir),m=0,1),dir=0,nd-1)
                    write(*,'("                                                                      gid=",6i3)') ((gid(m,dir),m=0,1),dir=0,nd-1)
                  end if
                   ! ------ Assigned extended boundaries with dirichlet values ------
                    do i3=nn3a,nn3b
                    do i2=nn2a,nn2b
                    do i1=nn1a,nn1b
                       ! Use forcing in array bcf(..)  
                       stop 1234 
                    end do
                    end do
                    end do
               end if
            else if( boundaryCondition(side,axis).eq.tractionBC )then 
              ! Extrapolation now done above
            else if( boundaryCondition(side,axis).eq.slipWall )then
              ! set n.u = given on the boundary 
              an1=0.  ! (an1,an2) = outward normal 
              an2=0.
              if( axis.eq.0 )then
                an1=2*side-1
              else
                an2=2*side-1
              end if
              if( addBoundaryForcing(side,axis).eq.0 )then
                ! no forcing 
                 i3=n3a
                 do i2=n2a,n2b
                 do i1=n1a,n1b
                 u1 = u(i1,i2,i3,u1c)
                 v1 = u(i1,i2,i3,u2c)
                 nDotU = an1*u1 + an2*v1  
                 u(i1,i2,i3,u1c)=u1 - nDotU*an1
                 u(i1,i2,i3,u2c)=v1 - nDotU*an2
                 end do
                 end do
               else if( assignTwilightZone.eq.0 )then
                 ! include forcing terms 
                 ! n.u = n.g 
                 i3=n3a
                 do i2=n2a,n2b
                 do i1=n1a,n1b
                 u1 = u(i1,i2,i3,u1c)
                 v1 = u(i1,i2,i3,u2c)
                 nDotU = an1*(u1-bcf(side,axis,i1,i2,i3,u1c)) + an2*(v1-bcf(side,axis,i1,i2,i3,u2c))
                 u(i1,i2,i3,u1c)=u1 - nDotU*an1
                 u(i1,i2,i3,u2c)=v1 - nDotU*an2
                 end do
                 end do
               else
                ! Twilight-zone: 
                !   n.u = n.ue
                 i3=n3a
                 do i2=n2a,n2b
                 do i1=n1a,n1b
                  call ogf2d(ep,xy(i1,i2,i3,0),xy(i1,i2,i3,1),t,u0,v0)
                 u1 = u(i1,i2,i3,u1c)
                 v1 = u(i1,i2,i3,u2c)
                 nDotU = an1*(u1-u0) + an2*(v1-v0)
                 u(i1,i2,i3,u1c)=u1 - nDotU*an1
                 u(i1,i2,i3,u2c)=v1 - nDotU*an2
                 end do
                 end do
               end if
               ! extrap values to the ghost line 
                i3=n3a
                do i2=nn2a,nn2b
                do i1=nn1a,nn1b
                u(i1-is1,i2-is2,i3,u1c)=(3.*u(i1,i2,i3,u1c)-3.*u(i1+is1,i2+is2,i3+is3,u1c)+u(i1+2*is1,i2+2*is2,i3+2*is3,u1c))
                u(i1-is1,i2-is2,i3,u2c)=(3.*u(i1,i2,i3,u2c)-3.*u(i1+is1,i2+is2,i3+is3,u2c)+u(i1+2*is1,i2+2*is2,i3+2*is3,u2c))
                end do
                end do
            else if( boundaryCondition(side,axis).eq.dirichletBoundaryCondition )then
               ! Set exact values on boundary and ghost
               if( twilightZone.eq.1 )then
                  do i3=nn3a,nn3b
                  do i2=nn2a,nn2b
                  do i1=nn1a,nn1b
                      call ogf3d(ep,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,u0,v0,w0)
                     u(i1,i2,i3,u1c)=u0
                     u(i1,i2,i3,u2c)=v0
                     u(i1,i2,i3,u3c)=w0
                   do ghost=1,numGhost
                     j1 = i1 - ghost*is1
                     j2 = i2 - ghost*is2
                       j3 = i3 - ghost*is3
                        call ogf3d(ep,xy(j1,j2,j3,0),xy(j1,j2,j3,1),xy(j1,j2,j3,2),t,u0,v0,w0)
                       u(j1,j2,j3,u1c)=u0
                       u(j1,j2,j3,u2c)=v0              
                       u(j1,j2,j3,u3c)=w0              
                   end do
                  end do
                  end do
                  end do
               end if
            else if( boundaryCondition(side,axis).eq.symmetry )then
               write(*,*) 'finish me for symmetry BC'
              stop 2856
            else if( boundaryCondition(side,axis).gt.0 )then
              stop 1193
            end if
            end do ! end side
            end do ! end axis
         ! ---------------------------------------------------------------------
         ! ----------- STAGE 2, ASSIGN NEUMANN TYPE CONDITIONS ---------------
         ! --------------------------------------------------------------------- 
          extra1a=numGhost
          extra1b=numGhost
          extra2a=numGhost
          extra2b=numGhost
          if( nd.eq.3 )then
            extra3a=numGhost
            extra3b=numGhost
          else
            extra3a=0
            extra3b=0
          end if
          if( boundaryCondition(0,0).lt.0 )then
            extra1a=max(0,extra1a) ! over-ride numGhost=-1 : assign ends in periodic directions (or internal parallel boundaries)
          else if( boundaryCondition(0,0).eq.0 )then
            extra1a=numGhost  ! include interpolation points since we assign ghost points outside these
          end if
          ! **NOTE** the bc on the right may be negative even it is not on the left (for parallel)
          if( boundaryCondition(1,0).lt.0 )then
            extra1b=max(0,extra1b) ! over-ride numGhost=-1 : assign ends in periodic directions
          else if( boundaryCondition(1,0).eq.0 )then
            extra1b=numGhost
          end if
          if( boundaryCondition(0,1).lt.0 )then
            extra2a=max(0,extra2a) ! over-ride numGhost=-1 : assign ends in periodic directions (or internal parallel boundaries)
          else if( boundaryCondition(0,1).eq.0 )then
            extra2a=numGhost  ! include interpolation points since we assign ghost points outside these
          end if
          ! **NOTE** the bc on the right may be negative even it is not on the left (for parallel)
          if( boundaryCondition(1,1).lt.0 )then
            extra2b=max(0,extra2b) ! over-ride numGhost=-1 : assign ends in periodic directions
          else if( boundaryCondition(1,1).eq.0 )then
            extra2b=numGhost
          end if
          if(  nd.eq.3 )then
           if( boundaryCondition(0,2).lt.0 )then
             extra3a=max(0,extra3a) ! over-ride numGhost=-1 : assign ends in periodic directions (or internal parallel boundaries)
           else if( boundaryCondition(0,2).eq.0 )then
             extra3a=numGhost  ! include interpolation points since we assign ghost points outside these
           end if
           ! **NOTE** the bc on the right may be negative even it is not on the left (for parallel)
           if( boundaryCondition(1,2).lt.0 )then
             extra3b=max(0,extra3b) ! over-ride numGhost=-1 : assign ends in periodic directions
           else if( boundaryCondition(1,2).eq.0 )then
             extra3b=numGhost
           end if
          end if
          do axis=0,nd-1
          do side=0,1
            if( boundaryCondition(side,axis).gt.0 )then
              ! write(*,'(" bcOpt: side,axis,bc=",3i2)') side,axis,boundaryCondition(side,axis)
              ! get gid(0;1,0:2) -- adjusted for priority at corners
                  do leftRight=0,1
                    do dir=0,2
                      gid(leftRight,dir) = gridIndexRange(leftRight,dir)
                    end do
                  end do
                  bcSide = boundaryCondition(side,axis)
                  if( bcSide.eq.tractionBC .or. bcSide.eq.displacementBC )then
                    ! adjust gridIndexRange at ends of the boundary 
                    !  Note: traction-traction corner assigns the extended ghost on the boundary
                    do dir=1,nd-1
                      axisp = mod(axis+dir,nd)  ! tangential direction
                      do leftRight=0,1
                        bcAdjacent = boundaryCondition(leftRight,axisp)
                        adjustEnds=.false.
                        if( bcSide.eq.displacementBC .and. ( bcAdjacent.eq.displacementBC .or. bcAdjacent.eq.dirichletBoundaryCondition ) )then
                          ! extended boundary has been set, no need to fill ghost on ends
                          adjustEnds=.true.
                        end if 
                        if( bcSide.eq.tractionBC .and. ( bcAdjacent.eq.displacementBC .or. bcAdjacent.eq.dirichletBoundaryCondition ) )then  
                          ! extended boundary has been set, no need to fill ghost on ends
                          adjustEnds=.true.
                        end if
                        ! if( boundaryCondition(leftRight,axisp)==dirichletBoundaryCondition .or. !     boundaryCondition(leftRight,axisp)==displacementBC             .or. !     boundaryCondition(leftRight,axisp)==tractionBC      )then
                        if( adjustEnds )then
                          ! shift to avoid the corner
                          gid(leftRight,axisp) = gid(leftRight,axisp) + (1-2*leftRight)
                        end if
                      end do
                    end do
                  end if
              axisp1=mod(axis+1,nd)
              axisp2=mod(axis+2,nd)
              i3=n3a
              n1a=gid(0,0)
              n1b=gid(1,0)
              n2a=gid(0,1)
              n2b=gid(1,1)
              n3a=gid(0,2)
              n3b=gid(1,2)
              if( axis.eq.0 )then
                n1a=gid(side,axis)
                n1b=gid(side,axis)
              else if( axis.eq.1 )then
                n2a=gid(side,axis)
                n2b=gid(side,axis)
              else
                n3a=gid(side,axis)
                n3b=gid(side,axis)
              end if
              ! if( boundaryCondition(side,axis)==tractionBC )then
              !   if( boundaryCondition(axisp1,0)==dirichletBoundaryCondition )then
              !     n2a = n2a+1; ! skip boundary
              !   end if
              !   if( boundaryCondition(axisp1,1)==dirichletBoundaryCondition )then
              !     n2b = n2b-1; ! skip boundary
              !   end if 
              ! end if
              nn1a=gid(0,0)-extra1a
              nn1b=gid(1,0)+extra1b
              nn2a=gid(0,1)-extra2a
              nn2b=gid(1,1)+extra2b
              nn3a=gid(0,2)-extra3a
              nn3b=gid(1,2)+extra3b
              if( axis.eq.0 )then
                nn1a=gid(side,axis)
                nn1b=gid(side,axis)
              else if( axis.eq.1 )then
                nn2a=gid(side,axis)
                nn2b=gid(side,axis)
              else
                nn3a=gid(side,axis)
                nn3b=gid(side,axis)
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
              if( debug.gt.7 )then
                write(*,'(" bcOpt: grid,side,axis=",3i3,",loop bounds: n1a,n1b,n2a,n2b,n3a,n3b=",6i3)') grid,side,axis,n1a,n1b,n2a,n2b,n3a,n3b
                write(*,'("                                                       gridIndexRange=",6i3)') ((gridIndexRange(m,dir),m=0,1),dir=0,nd-1)
                write(*,'("                                                                  gid=",6i3)') ((gid(m,dir),m=0,1),dir=0,nd-1)
              end if
            end if ! if bc>0 
            ! On interfaces we should use the bcf array values even for TZ since then
            ! we get a coupling at the interface: 
            !   bcf = n.sigma(fluid) + [ n.sigma_e(solid) - n.sigma_e(fluid) ]
            if( interfaceType(side,axis,grid).eq.noInterface )then
              assignTwilightZone=twilightZone
            else
              assignTwilightZone=0  ! this will turn off the use of TZ
            end if
           if( boundaryCondition(side,axis).eq.displacementBC )then 
             ! --------- DISPLACEMENT GHOST -----
             ! Note: addBoundaryForcing =1 for TZ as well
             if( twilightZone==1 )then
                 ! ----------- Assign ghost points ----------
                 fe(0)=0.; fe(1)=0.; fe(2)=0.; fe(3)=0.;   ! holds forcing 
                 numberOfEquations=4;                      ! number of ghost points we solve for   
                if( .false. )then
                  write(*,*) 'Displacement BC -- assign GHOST:'
                  write(*,'(" bcOpt: grid,side,axis=",3i3,",loop bounds: n1a,n1b,n2a,n2b,n3a,n3b=",6i3)') grid,side,axis,n1a,n1b,n2a,n2b,n3a,n3b
                  write(*,'("                                                     gridIndexRange=",6i3)') ((gridIndexRange(m,dir),m=0,1),dir=0,nd-1)
                  write(*,'("                                                                gid=",6i3)') ((gid(m,dir),m=0,1),dir=0,nd-1)
                end if
                  do i3=n3a,n3b
                  do i2=n2a,n2b
                  do i1=n1a,n1b
                   !   Set div(u) = 0 
                   !      u.x + v.y + w.z = 0
                       ! --- displacement ghost : curvilinear  -----
                          an1 = rsxy(i1,i2,i3,axis,0)
                          an2 = rsxy(i1,i2,i3,axis,1)
                          an3 = rsxy(i1,i2,i3,axis,2)
                          aNormi = -is/max(epsx,sqrt(an1**2 + an2**2+ an3**2))
                          an1=an1*aNormi
                          an2=an2*aNormi
                          an3=an3*aNormi
                         ! --- 3D Displacement ghost curvlinear ----
                         ! here are  the equations we mean to satisfy: 
                         !   u.x + v.y + w.z = 0
                         !   Extrapolate t.uv
                         ! This can be written as: 
                         !   div(u)* n  + (I - n n^T)*( u - extrap( u ) ) 
                         ! 
                         ! first ghost:
                         j1 = i1-is1; j2 = i2-is2; j3 = i3-is3;
                         ux = ux23(i1,i2,i3,u1c)
                         vy = uy23(i1,i2,i3,u2c)          
                         wz = uz23(i1,i2,i3,u3c)   
                         divu = ux + vy + wz - fe(0) 
                         g(0) = u(j1,j2,j3,u1c)-(3.*u(j1+is1,j2+is2,j3+is3,u1c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u1c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u1c))
                         g(1) = u(j1,j2,j3,u2c)-(3.*u(j1+is1,j2+is2,j3+is3,u2c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u2c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u2c))                 
                         g(2) = u(j1,j2,j3,u3c)-(3.*u(j1+is1,j2+is2,j3+is3,u3c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u3c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u3c))    
                         f(0) = divu*an1 + (1.-an1*an1)*g(0) + (  -an1*an2)*g(1) +(  -an1*an3)*g(2)               
                         f(1) = divu*an2 + (  -an2*an1)*g(0) + (1.-an2*an2)*g(1) +(  -an2*an3)*g(2)               
                         f(2) = divu*an3 + (  -an3*an1)*g(0) + (  -an3*an2)*g(1) +(1.-an3*an3)*g(2)    
                         a3(0,0)=-is*( rsxy(i1,i2,i3,axis,0) )/(2.*dr(axis)) *an1 + (1.-an1*an1)  ! coeff of u1(-1) in f(0)
                         a3(0,1)=-is*( rsxy(i1,i2,i3,axis,1) )/(2.*dr(axis)) *an1 + (  -an1*an2)  ! coeff of u2(-1) in f(0)
                         a3(0,2)=-is*( rsxy(i1,i2,i3,axis,2) )/(2.*dr(axis)) *an1 + (  -an1*an3)  ! coeff of u3(-1) in f(0)
                         a3(1,0)=-is*( rsxy(i1,i2,i3,axis,0) )/(2.*dr(axis)) *an2 + (  -an2*an1)  ! coeff of u1(-1) in f(1)
                         a3(1,1)=-is*( rsxy(i1,i2,i3,axis,1) )/(2.*dr(axis)) *an2 + (1.-an2*an2)
                         a3(1,2)=-is*( rsxy(i1,i2,i3,axis,2) )/(2.*dr(axis)) *an2 + (  -an2*an3)
                         a3(2,0)=-is*( rsxy(i1,i2,i3,axis,0) )/(2.*dr(axis)) *an3 + (  -an3*an1)
                         a3(2,1)=-is*( rsxy(i1,i2,i3,axis,1) )/(2.*dr(axis)) *an3 + (  -an3*an2)
                         a3(2,2)=-is*( rsxy(i1,i2,i3,axis,2) )/(2.*dr(axis)) *an3 + (1.-an3*an3)                    
                         ! here are the wrong ghostpoint values
                         q(0) = u(i1-is1,i2-is2,i3-is3,u1c)
                         q(1) = u(i1-is1,i2-is2,i3-is3,u2c)
                         q(2) = u(i1-is1,i2-is2,i3-is3,u3c)
                         ! subtract off the contributions from the wrong values at the ghost points:
                         do n=0,2
                           f(n) = ( a3(n,0)*q(0) +a3(n,1)*q(1)+ a3(n,2)*q(2) ) - f(n)
                         end do
                         ! could optimize and do this by hand 
                         call dgeco( a3(0,0), 3, 3, ipvt(0),rcond,work(0))
                         call dgesl( a3(0,0), 3, 3, ipvt(0), f(0), job)
                         u(i1-is1,i2-is2,i3-is3,u1c)=f(0)
                         u(i1-is1,i2-is2,i3-is3,u2c)=f(1)
                         u(i1-is1,i2-is2,i3-is3,u3c)=f(2)
                         if( .true. .or. checkResiduals )then 
                           ! re-evaluate and check residuals
                           ux = ux23(i1,i2,i3,u1c)
                           vy = uy23(i1,i2,i3,u2c)          
                           wz = uz23(i1,i2,i3,u3c)   
                           divu = ux + vy + wz - fe(0) 
                           g(0) = u(j1,j2,j3,u1c)-(3.*u(j1+is1,j2+is2,j3+is3,u1c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u1c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u1c))
                           g(1) = u(j1,j2,j3,u2c)-(3.*u(j1+is1,j2+is2,j3+is3,u2c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u2c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u2c))                 
                           g(2) = u(j1,j2,j3,u3c)-(3.*u(j1+is1,j2+is2,j3+is3,u3c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u3c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u3c))    
                           f(0) = divu*an1 + (1.-an1*an1)*g(0) + (  -an1*an2)*g(1) +(  -an1*an3)*g(2)               
                           f(1) = divu*an2 + (  -an2*an1)*g(0) + (1.-an2*an2)*g(1) +(  -an2*an3)*g(2)               
                           f(2) = divu*an3 + (  -an3*an1)*g(0) + ( - an3*an2)*g(1) +(1.-an3*an3)*g(2)               
                           resMax = max(abs(f(0)),abs(f(1)),abs(f(2)))      
                           if( resMax > resTol )then
                             write(*,'(" --> displacement BC:ERROR curvilinear fill ghost 3D: residuals =",3(1pe12.4))') f(0),f(1),f(2)
                             stop 3232
                           end if
                         end if
                     ! --- Extrapolate any extra ghost --- 
                     do ghost=2,numGhost
                       ! (j1,j2,j3) is the ghost point index
                       j1 = i1 - is1*ghost 
                       j2 = i2 - is2*ghost 
                       j3 = i3 - is3*ghost 
                       u(j1,j2,j3,u1c)=(3.*u(j1+is1,j2+is2,j3+is3,u1c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u1c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u1c))
                       u(j1,j2,j3,u2c)=(3.*u(j1+is1,j2+is2,j3+is3,u2c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u2c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u2c))
                         u(j1,j2,j3,u3c)=(3.*u(j1+is1,j2+is2,j3+is3,u3c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u3c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u3c))
                     end do 
                  end do
                  end do
                  end do
                 ! write(*,'("Traction BC -- stop here for now")')
                 ! stop 9999
             else if( .true. .or. addBoundaryForcing(side,axis)==0 )then  ! FOR NOW SKIP addBoundaryFORCING
                 ! ----------- Assign ghost points ----------
                 fe(0)=0.; fe(1)=0.; fe(2)=0.; fe(3)=0.;   ! holds forcing 
                 numberOfEquations=4;                      ! number of ghost points we solve for   
                if( .false. )then
                  write(*,*) 'Displacement BC -- assign GHOST:'
                  write(*,'(" bcOpt: grid,side,axis=",3i3,",loop bounds: n1a,n1b,n2a,n2b,n3a,n3b=",6i3)') grid,side,axis,n1a,n1b,n2a,n2b,n3a,n3b
                  write(*,'("                                                     gridIndexRange=",6i3)') ((gridIndexRange(m,dir),m=0,1),dir=0,nd-1)
                  write(*,'("                                                                gid=",6i3)') ((gid(m,dir),m=0,1),dir=0,nd-1)
                end if
                  do i3=n3a,n3b
                  do i2=n2a,n2b
                  do i1=n1a,n1b
                   !   Set div(u) = 0 
                   !      u.x + v.y + w.z = 0
                       ! --- displacement ghost : curvilinear  -----
                          an1 = rsxy(i1,i2,i3,axis,0)
                          an2 = rsxy(i1,i2,i3,axis,1)
                          an3 = rsxy(i1,i2,i3,axis,2)
                          aNormi = -is/max(epsx,sqrt(an1**2 + an2**2+ an3**2))
                          an1=an1*aNormi
                          an2=an2*aNormi
                          an3=an3*aNormi
                         ! --- 3D Displacement ghost curvlinear ----
                         ! here are  the equations we mean to satisfy: 
                         !   u.x + v.y + w.z = 0
                         !   Extrapolate t.uv
                         ! This can be written as: 
                         !   div(u)* n  + (I - n n^T)*( u - extrap( u ) ) 
                         ! 
                         ! first ghost:
                         j1 = i1-is1; j2 = i2-is2; j3 = i3-is3;
                         ux = ux23(i1,i2,i3,u1c)
                         vy = uy23(i1,i2,i3,u2c)          
                         wz = uz23(i1,i2,i3,u3c)   
                         divu = ux + vy + wz - fe(0) 
                         g(0) = u(j1,j2,j3,u1c)-(3.*u(j1+is1,j2+is2,j3+is3,u1c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u1c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u1c))
                         g(1) = u(j1,j2,j3,u2c)-(3.*u(j1+is1,j2+is2,j3+is3,u2c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u2c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u2c))                 
                         g(2) = u(j1,j2,j3,u3c)-(3.*u(j1+is1,j2+is2,j3+is3,u3c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u3c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u3c))    
                         f(0) = divu*an1 + (1.-an1*an1)*g(0) + (  -an1*an2)*g(1) +(  -an1*an3)*g(2)               
                         f(1) = divu*an2 + (  -an2*an1)*g(0) + (1.-an2*an2)*g(1) +(  -an2*an3)*g(2)               
                         f(2) = divu*an3 + (  -an3*an1)*g(0) + (  -an3*an2)*g(1) +(1.-an3*an3)*g(2)    
                         a3(0,0)=-is*( rsxy(i1,i2,i3,axis,0) )/(2.*dr(axis)) *an1 + (1.-an1*an1)  ! coeff of u1(-1) in f(0)
                         a3(0,1)=-is*( rsxy(i1,i2,i3,axis,1) )/(2.*dr(axis)) *an1 + (  -an1*an2)  ! coeff of u2(-1) in f(0)
                         a3(0,2)=-is*( rsxy(i1,i2,i3,axis,2) )/(2.*dr(axis)) *an1 + (  -an1*an3)  ! coeff of u3(-1) in f(0)
                         a3(1,0)=-is*( rsxy(i1,i2,i3,axis,0) )/(2.*dr(axis)) *an2 + (  -an2*an1)  ! coeff of u1(-1) in f(1)
                         a3(1,1)=-is*( rsxy(i1,i2,i3,axis,1) )/(2.*dr(axis)) *an2 + (1.-an2*an2)
                         a3(1,2)=-is*( rsxy(i1,i2,i3,axis,2) )/(2.*dr(axis)) *an2 + (  -an2*an3)
                         a3(2,0)=-is*( rsxy(i1,i2,i3,axis,0) )/(2.*dr(axis)) *an3 + (  -an3*an1)
                         a3(2,1)=-is*( rsxy(i1,i2,i3,axis,1) )/(2.*dr(axis)) *an3 + (  -an3*an2)
                         a3(2,2)=-is*( rsxy(i1,i2,i3,axis,2) )/(2.*dr(axis)) *an3 + (1.-an3*an3)                    
                         ! here are the wrong ghostpoint values
                         q(0) = u(i1-is1,i2-is2,i3-is3,u1c)
                         q(1) = u(i1-is1,i2-is2,i3-is3,u2c)
                         q(2) = u(i1-is1,i2-is2,i3-is3,u3c)
                         ! subtract off the contributions from the wrong values at the ghost points:
                         do n=0,2
                           f(n) = ( a3(n,0)*q(0) +a3(n,1)*q(1)+ a3(n,2)*q(2) ) - f(n)
                         end do
                         ! could optimize and do this by hand 
                         call dgeco( a3(0,0), 3, 3, ipvt(0),rcond,work(0))
                         call dgesl( a3(0,0), 3, 3, ipvt(0), f(0), job)
                         u(i1-is1,i2-is2,i3-is3,u1c)=f(0)
                         u(i1-is1,i2-is2,i3-is3,u2c)=f(1)
                         u(i1-is1,i2-is2,i3-is3,u3c)=f(2)
                         if( .true. .or. checkResiduals )then 
                           ! re-evaluate and check residuals
                           ux = ux23(i1,i2,i3,u1c)
                           vy = uy23(i1,i2,i3,u2c)          
                           wz = uz23(i1,i2,i3,u3c)   
                           divu = ux + vy + wz - fe(0) 
                           g(0) = u(j1,j2,j3,u1c)-(3.*u(j1+is1,j2+is2,j3+is3,u1c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u1c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u1c))
                           g(1) = u(j1,j2,j3,u2c)-(3.*u(j1+is1,j2+is2,j3+is3,u2c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u2c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u2c))                 
                           g(2) = u(j1,j2,j3,u3c)-(3.*u(j1+is1,j2+is2,j3+is3,u3c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u3c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u3c))    
                           f(0) = divu*an1 + (1.-an1*an1)*g(0) + (  -an1*an2)*g(1) +(  -an1*an3)*g(2)               
                           f(1) = divu*an2 + (  -an2*an1)*g(0) + (1.-an2*an2)*g(1) +(  -an2*an3)*g(2)               
                           f(2) = divu*an3 + (  -an3*an1)*g(0) + ( - an3*an2)*g(1) +(1.-an3*an3)*g(2)               
                           resMax = max(abs(f(0)),abs(f(1)),abs(f(2)))      
                           if( resMax > resTol )then
                             write(*,'(" --> displacement BC:ERROR curvilinear fill ghost 3D: residuals =",3(1pe12.4))') f(0),f(1),f(2)
                             stop 3232
                           end if
                         end if
                     ! --- Extrapolate any extra ghost --- 
                     do ghost=2,numGhost
                       ! (j1,j2,j3) is the ghost point index
                       j1 = i1 - is1*ghost 
                       j2 = i2 - is2*ghost 
                       j3 = i3 - is3*ghost 
                       u(j1,j2,j3,u1c)=(3.*u(j1+is1,j2+is2,j3+is3,u1c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u1c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u1c))
                       u(j1,j2,j3,u2c)=(3.*u(j1+is1,j2+is2,j3+is3,u2c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u2c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u2c))
                         u(j1,j2,j3,u3c)=(3.*u(j1+is1,j2+is2,j3+is3,u3c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u3c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u3c))
                     end do 
                  end do
                  end do
                  end do
                 ! write(*,'("Traction BC -- stop here for now")')
                 ! stop 9999
             else 
                 ! ----------- Assign ghost points ----------
                 fe(0)=0.; fe(1)=0.; fe(2)=0.; fe(3)=0.;   ! holds forcing 
                 numberOfEquations=4;                      ! number of ghost points we solve for   
                if( .false. )then
                  write(*,*) 'Displacement BC -- assign GHOST:'
                  write(*,'(" bcOpt: grid,side,axis=",3i3,",loop bounds: n1a,n1b,n2a,n2b,n3a,n3b=",6i3)') grid,side,axis,n1a,n1b,n2a,n2b,n3a,n3b
                  write(*,'("                                                     gridIndexRange=",6i3)') ((gridIndexRange(m,dir),m=0,1),dir=0,nd-1)
                  write(*,'("                                                                gid=",6i3)') ((gid(m,dir),m=0,1),dir=0,nd-1)
                end if
                  do i3=n3a,n3b
                  do i2=n2a,n2b
                  do i1=n1a,n1b
                   !   Set div(u) = 0 
                   !      u.x + v.y + w.z = 0
                       ! --- displacement ghost : curvilinear  -----
                          an1 = rsxy(i1,i2,i3,axis,0)
                          an2 = rsxy(i1,i2,i3,axis,1)
                          an3 = rsxy(i1,i2,i3,axis,2)
                          aNormi = -is/max(epsx,sqrt(an1**2 + an2**2+ an3**2))
                          an1=an1*aNormi
                          an2=an2*aNormi
                          an3=an3*aNormi
                        if( assignTwilightZone.eq.1 )then
                           call ogDeriv(ep,0,1,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,u1c,ux0)
                           call ogDeriv(ep,0,0,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,u2c,vy0)
                           call ogDeriv(ep,0,0,0,1,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,u3c,wz0)            
                           fe(0) = ux0 + vy0 + wz0
                         fe(1) = 0.
                        end if
                         ! --- 3D Displacement ghost curvlinear ----
                         ! here are  the equations we mean to satisfy: 
                         !   u.x + v.y + w.z = 0
                         !   Extrapolate t.uv
                         ! This can be written as: 
                         !   div(u)* n  + (I - n n^T)*( u - extrap( u ) ) 
                         ! 
                         ! first ghost:
                         j1 = i1-is1; j2 = i2-is2; j3 = i3-is3;
                         ux = ux23(i1,i2,i3,u1c)
                         vy = uy23(i1,i2,i3,u2c)          
                         wz = uz23(i1,i2,i3,u3c)   
                         divu = ux + vy + wz - fe(0) 
                         g(0) = u(j1,j2,j3,u1c)-(3.*u(j1+is1,j2+is2,j3+is3,u1c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u1c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u1c))
                         g(1) = u(j1,j2,j3,u2c)-(3.*u(j1+is1,j2+is2,j3+is3,u2c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u2c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u2c))                 
                         g(2) = u(j1,j2,j3,u3c)-(3.*u(j1+is1,j2+is2,j3+is3,u3c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u3c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u3c))    
                         f(0) = divu*an1 + (1.-an1*an1)*g(0) + (  -an1*an2)*g(1) +(  -an1*an3)*g(2)               
                         f(1) = divu*an2 + (  -an2*an1)*g(0) + (1.-an2*an2)*g(1) +(  -an2*an3)*g(2)               
                         f(2) = divu*an3 + (  -an3*an1)*g(0) + (  -an3*an2)*g(1) +(1.-an3*an3)*g(2)    
                         a3(0,0)=-is*( rsxy(i1,i2,i3,axis,0) )/(2.*dr(axis)) *an1 + (1.-an1*an1)  ! coeff of u1(-1) in f(0)
                         a3(0,1)=-is*( rsxy(i1,i2,i3,axis,1) )/(2.*dr(axis)) *an1 + (  -an1*an2)  ! coeff of u2(-1) in f(0)
                         a3(0,2)=-is*( rsxy(i1,i2,i3,axis,2) )/(2.*dr(axis)) *an1 + (  -an1*an3)  ! coeff of u3(-1) in f(0)
                         a3(1,0)=-is*( rsxy(i1,i2,i3,axis,0) )/(2.*dr(axis)) *an2 + (  -an2*an1)  ! coeff of u1(-1) in f(1)
                         a3(1,1)=-is*( rsxy(i1,i2,i3,axis,1) )/(2.*dr(axis)) *an2 + (1.-an2*an2)
                         a3(1,2)=-is*( rsxy(i1,i2,i3,axis,2) )/(2.*dr(axis)) *an2 + (  -an2*an3)
                         a3(2,0)=-is*( rsxy(i1,i2,i3,axis,0) )/(2.*dr(axis)) *an3 + (  -an3*an1)
                         a3(2,1)=-is*( rsxy(i1,i2,i3,axis,1) )/(2.*dr(axis)) *an3 + (  -an3*an2)
                         a3(2,2)=-is*( rsxy(i1,i2,i3,axis,2) )/(2.*dr(axis)) *an3 + (1.-an3*an3)                    
                         ! here are the wrong ghostpoint values
                         q(0) = u(i1-is1,i2-is2,i3-is3,u1c)
                         q(1) = u(i1-is1,i2-is2,i3-is3,u2c)
                         q(2) = u(i1-is1,i2-is2,i3-is3,u3c)
                         ! subtract off the contributions from the wrong values at the ghost points:
                         do n=0,2
                           f(n) = ( a3(n,0)*q(0) +a3(n,1)*q(1)+ a3(n,2)*q(2) ) - f(n)
                         end do
                         ! could optimize and do this by hand 
                         call dgeco( a3(0,0), 3, 3, ipvt(0),rcond,work(0))
                         call dgesl( a3(0,0), 3, 3, ipvt(0), f(0), job)
                         u(i1-is1,i2-is2,i3-is3,u1c)=f(0)
                         u(i1-is1,i2-is2,i3-is3,u2c)=f(1)
                         u(i1-is1,i2-is2,i3-is3,u3c)=f(2)
                         if( .true. .or. checkResiduals )then 
                           ! re-evaluate and check residuals
                           ux = ux23(i1,i2,i3,u1c)
                           vy = uy23(i1,i2,i3,u2c)          
                           wz = uz23(i1,i2,i3,u3c)   
                           divu = ux + vy + wz - fe(0) 
                           g(0) = u(j1,j2,j3,u1c)-(3.*u(j1+is1,j2+is2,j3+is3,u1c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u1c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u1c))
                           g(1) = u(j1,j2,j3,u2c)-(3.*u(j1+is1,j2+is2,j3+is3,u2c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u2c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u2c))                 
                           g(2) = u(j1,j2,j3,u3c)-(3.*u(j1+is1,j2+is2,j3+is3,u3c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u3c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u3c))    
                           f(0) = divu*an1 + (1.-an1*an1)*g(0) + (  -an1*an2)*g(1) +(  -an1*an3)*g(2)               
                           f(1) = divu*an2 + (  -an2*an1)*g(0) + (1.-an2*an2)*g(1) +(  -an2*an3)*g(2)               
                           f(2) = divu*an3 + (  -an3*an1)*g(0) + ( - an3*an2)*g(1) +(1.-an3*an3)*g(2)               
                           resMax = max(abs(f(0)),abs(f(1)),abs(f(2)))      
                           if( resMax > resTol )then
                             write(*,'(" --> displacement BC:ERROR curvilinear fill ghost 3D: residuals =",3(1pe12.4))') f(0),f(1),f(2)
                             stop 3232
                           end if
                         end if
                     ! --- Extrapolate any extra ghost --- 
                     do ghost=2,numGhost
                       ! (j1,j2,j3) is the ghost point index
                       j1 = i1 - is1*ghost 
                       j2 = i2 - is2*ghost 
                       j3 = i3 - is3*ghost 
                       u(j1,j2,j3,u1c)=(3.*u(j1+is1,j2+is2,j3+is3,u1c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u1c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u1c))
                       u(j1,j2,j3,u2c)=(3.*u(j1+is1,j2+is2,j3+is3,u2c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u2c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u2c))
                         u(j1,j2,j3,u3c)=(3.*u(j1+is1,j2+is2,j3+is3,u3c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u3c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u3c))
                     end do 
                  end do
                  end do
                  end do
                 ! write(*,'("Traction BC -- stop here for now")')
                 ! stop 9999
             end if
          else if( boundaryCondition(side,axis).eq.tractionBC )then 
            ! ------ TRACTION BC ----
            !  if( addBoundaryForcing(side,axis).eq.0 .and. assignTwilightZone==0 )then   *** FIX ME ***
            ! if( addBoundaryForcing(side,axis)==0 .and. twilightZone==0 )then   ! FOR NOW SKIP addBoundaryFORCING  *** FIX ME *****
            if( twilightZone==0 )then   ! FOR NOW SKIP addBoundaryFORCING  *** FIX ME ***** addBoundaryForcing =1 for an exact solution 
                fe(0)=0.; fe(1)=0.;  fe(2)=0.; fe(3)=0.;   ! holds forcing 
                numberOfEquations=4;      ! number of ghost points we solve for   
                 do i3=n3a,n3b
                 do i2=n2a,n2b
                 do i1=n1a,n1b
                  ! here is the normal
                       an1 = rsxy(i1,i2,i3,axis,0)
                       an2 = rsxy(i1,i2,i3,axis,1)
                       an3 = rsxy(i1,i2,i3,axis,2)
                       aNormi = -is/max(epsx,sqrt(an1**2 + an2**2+ an3**2))
                       an1=an1*aNormi
                       an2=an2*aNormi
                       an3=an3*aNormi
                  ! *** FIX ME : USE THIS BELOW FOR RECTANGULAR CASE####################################
                  ! #If "2" eq "2"
                  !   #If "noForcing" eq "twilightZone" 
                  !     #If "3" eq "2"
                  !       OGDERIV2D(0,1,0,0,i1,i2,i3,t,ux0,vx0)
                  !       OGDERIV2D(0,0,1,0,i1,i2,i3,t,uy0,vy0)
                  !       fe(0) = ux0 + vy0
                  !       fe(1) = 2*an1*an2*( ux0 - vy0 ) + ( an2**2 -an1**2 )*( uy0 + vx0 )
                  !     #Else
                  !       stop 333
                  !     #End
                  !   #End
                  ! #End
                    ! ----- TRACTION CURVILINEAR ----
                    ! --- here are  the equations we mean to satisfy ----
                    !  --> evaluate with the wrong values in the ghost  
                      ! --------- 2 2 -----------
                        ! --------------- TRACTION, 2=2, 3=3, CURVILINEAR -----------------
                        ! write(*,*) '--------------- TRACTION, order=2, dim=3, CURVILINEAR -----------------'          
                        numberOfEquations = 3 
                        ! Here are the realtive index positions of the 3 ghost we solve for 
                        shift1(0) = -is1; shift2(0) = -is2; shift3(0) = -is3; comp(0)=u1c; 
                        shift1(1) = -is1; shift2(1) = -is2; shift3(1) = -is3; comp(1)=u2c; 
                        shift1(2) = -is1; shift2(2) = -is2; shift3(2) = -is3; comp(2)=u3c; 
                        ! --- eval any forcing ---
                          !    (div(u)) n + mu* ( I - n n^T) epsm * n 
                        ! write(*,'("traction:3d: i1,i2,i3=",3i3," fe=",3(1pe12.4,1x))') i1,i2,i3,fe(0),fe(1),fe(2)
                        ! Fill in the matrix using the discrete delta approach
                          ! ! hw1 = half stencil width
                          ! hw1=orderOfAccuracy/2
                          ! hw2=hw1
                          ! if( nd.eq.2 )then
                          !   hw3=0
                          ! else
                          !   hw3=hw1
                          ! end if
                          !  write(*,'("EVAL-COEFF: i1,i2,i3=",3i3," hw1,hw2,hw3=",3i2)') i1,i2,i3,hw1,hw2,hw3
                          ! First eval equartions with no pertutbation --> save in f0 
                            !    (div(u)) n + mu* ( I - n n^T) epsm * n 
                            u1x = ux23(i1,i2,i3,u1c); u1y = uy23(i1,i2,i3,u1c); u1z = uz23(i1,i2,i3,u1c);
                            u2x = ux23(i1,i2,i3,u2c); u2y = uy23(i1,i2,i3,u2c); u2z = uz23(i1,i2,i3,u2c);
                            u3x = ux23(i1,i2,i3,u3c); u3y = uy23(i1,i2,i3,u3c); u3z = uz23(i1,i2,i3,u3c);
                            divu = u1x + u2y + u3z 
                            ! traction = epsm*n : 
                            trac1 = mu*( ( u1x + u1x )*an1 + ( u1y + u2x )*an2 + ( u1z + u3x )*an3 )
                            trac2 = mu*( ( u2x + u1y )*an1 + ( u2y + u2y )*an2 + ( u2z + u3y )*an3 )
                            trac3 = mu*( ( u3x + u1z )*an1 + ( u3y + u2z )*an2 + ( u3z + u3z )*an3 )
                            f(0) = divu*an1 +  (1.-an1*an1)*trac1 -    an1*an2 *trac2     - an1*an3 *trac3  - fe(0)
                            f(1) = divu*an2 +     -an2*an1 *trac1 +(1.-an2*an2)*trac2     - an2*an3 *trac3  - fe(1)
                            f(2) = divu*an3 +     -an3*an1 *trac1 +   -an3*an2 *trac2 + (1.-an3*an3)*trac3  - fe(2)         
                          do n1=0,numberOfEquations-1
                           f0(n1)=f(n1)
                          end do
                          delta=1.  ! perturb E by this amount 
                          do n2=0,numberOfEquations-1
                            ! pertub one component: 
                              u(i1+shift1(n2),i2+shift2(n2),i3+shift3(n2),comp(n2))=u(i1+shift1(n2),i2+shift2(n2),i3+shift3(n2),comp(n2))+(delta)
                              !    (div(u)) n + mu* ( I - n n^T) epsm * n 
                              u1x = ux23(i1,i2,i3,u1c); u1y = uy23(i1,i2,i3,u1c); u1z = uz23(i1,i2,i3,u1c);
                              u2x = ux23(i1,i2,i3,u2c); u2y = uy23(i1,i2,i3,u2c); u2z = uz23(i1,i2,i3,u2c);
                              u3x = ux23(i1,i2,i3,u3c); u3y = uy23(i1,i2,i3,u3c); u3z = uz23(i1,i2,i3,u3c);
                              divu = u1x + u2y + u3z 
                              ! traction = epsm*n : 
                              trac1 = mu*( ( u1x + u1x )*an1 + ( u1y + u2x )*an2 + ( u1z + u3x )*an3 )
                              trac2 = mu*( ( u2x + u1y )*an1 + ( u2y + u2y )*an2 + ( u2z + u3y )*an3 )
                              trac3 = mu*( ( u3x + u1z )*an1 + ( u3y + u2z )*an2 + ( u3z + u3z )*an3 )
                              f(0) = divu*an1 +  (1.-an1*an1)*trac1 -    an1*an2 *trac2     - an1*an3 *trac3  - fe(0)
                              f(1) = divu*an2 +     -an2*an1 *trac1 +(1.-an2*an2)*trac2     - an2*an3 *trac3  - fe(1)
                              f(2) = divu*an3 +     -an3*an1 *trac1 +   -an3*an2 *trac2 + (1.-an3*an3)*trac3  - fe(2)         
                            ! compute the difference
                            do n1=0,numberOfEquations-1
                              f(n1)=f(n1)-f0(n1)
                              a3(n1,n2) = f(n1)
                            end do
                            ! reset pertubation
                              u(i1+shift1(n2),i2+shift2(n2),i3+shift3(n2),comp(n2))=u(i1+shift1(n2),i2+shift2(n2),i3+shift3(n2),comp(n2))+(-delta)
                          end do 
                          ! restore f -- not needed here
                          ! evalTractionCurvilinear3d()
                        ! print the matrix: 
                        ! write(*,'("traction 3d curv: i1,i2=",2i3,/," a3=[",3(1pe12.4,1x),"; ...")') i1,i2,((a3(n,n2),n2=0,2),n=0,2) 
                        ! evaluate equations with wrong values in the ghost        
                          !    (div(u)) n + mu* ( I - n n^T) epsm * n 
                          u1x = ux23(i1,i2,i3,u1c); u1y = uy23(i1,i2,i3,u1c); u1z = uz23(i1,i2,i3,u1c);
                          u2x = ux23(i1,i2,i3,u2c); u2y = uy23(i1,i2,i3,u2c); u2z = uz23(i1,i2,i3,u2c);
                          u3x = ux23(i1,i2,i3,u3c); u3y = uy23(i1,i2,i3,u3c); u3z = uz23(i1,i2,i3,u3c);
                          divu = u1x + u2y + u3z 
                          ! traction = epsm*n : 
                          trac1 = mu*( ( u1x + u1x )*an1 + ( u1y + u2x )*an2 + ( u1z + u3x )*an3 )
                          trac2 = mu*( ( u2x + u1y )*an1 + ( u2y + u2y )*an2 + ( u2z + u3y )*an3 )
                          trac3 = mu*( ( u3x + u1z )*an1 + ( u3y + u2z )*an2 + ( u3z + u3z )*an3 )
                          f(0) = divu*an1 +  (1.-an1*an1)*trac1 -    an1*an2 *trac2     - an1*an3 *trac3  - fe(0)
                          f(1) = divu*an2 +     -an2*an1 *trac1 +(1.-an2*an2)*trac2     - an2*an3 *trac3  - fe(1)
                          f(2) = divu*an3 +     -an3*an1 *trac1 +   -an3*an2 *trac2 + (1.-an3*an3)*trac3  - fe(2)         
                        ! adjust rhs for wrong values on ghost 
                        do n=0,numberOfEquations-1
                          f(n) = -f(n) 
                          do n2=0,numberOfEquations-1
                            f(n) = f(n) + a3(n,n2)*u(i1+shift1(n2),i2+shift2(n2),i3+shift3(n2),comp(n2))
                          end do
                        end do
                        ! solve for ghost
                        call dgeco( a3(0,0), numberOfEquations, numberOfEquations, ipvt(0),rcond,work(0) )
                        call dgesl( a3(0,0), numberOfEquations, numberOfEquations, ipvt(0), f(0), job )
                        ! write(*,'(" traction curvilinear 3d: rcond=",1pe10.2)') rcond
                        ! assign values 
                        do n2=0,numberOfEquations-1
                          u(i1+shift1(n2),i2+shift2(n2),i3+shift3(n2),comp(n2)) = f(n2)
                        end do    
                        ! call ogDeriv(ep,0,0,0,0,xy(i1-is1,i2-is2,i3-is3,0),xy(i1-is1,i2-is2,i3-is3,1),xy(i1-is1,i2-is2,i3-is3,2),t,u1c,u1e)
                        ! call ogDeriv(ep,0,0,0,0,xy(i1-is1,i2-is2,i3-is3,0),xy(i1-is1,i2-is2,i3-is3,1),xy(i1-is1,i2-is2,i3-is3,2),t,u2c,u2e)
                        ! write(*,'("traction:3d: i1,i2,i3=",3i3," u1(-1)=",1pe12.4," err=",1pe10.2," u2(-1)=",1pe12.4," err=",1pe10.2)') !        i1,i2,i3,f(0),f(0)-u1e,f(1),f(1)-u2e           
                        ! ---- check residuals ---   
                        if( checkResiduals )then
                          ! re-evaluate the equations
                            !    (div(u)) n + mu* ( I - n n^T) epsm * n 
                            u1x = ux23(i1,i2,i3,u1c); u1y = uy23(i1,i2,i3,u1c); u1z = uz23(i1,i2,i3,u1c);
                            u2x = ux23(i1,i2,i3,u2c); u2y = uy23(i1,i2,i3,u2c); u2z = uz23(i1,i2,i3,u2c);
                            u3x = ux23(i1,i2,i3,u3c); u3y = uy23(i1,i2,i3,u3c); u3z = uz23(i1,i2,i3,u3c);
                            divu = u1x + u2y + u3z 
                            ! traction = epsm*n : 
                            trac1 = mu*( ( u1x + u1x )*an1 + ( u1y + u2x )*an2 + ( u1z + u3x )*an3 )
                            trac2 = mu*( ( u2x + u1y )*an1 + ( u2y + u2y )*an2 + ( u2z + u3y )*an3 )
                            trac3 = mu*( ( u3x + u1z )*an1 + ( u3y + u2z )*an2 + ( u3z + u3z )*an3 )
                            f(0) = divu*an1 +  (1.-an1*an1)*trac1 -    an1*an2 *trac2     - an1*an3 *trac3  - fe(0)
                            f(1) = divu*an2 +     -an2*an1 *trac1 +(1.-an2*an2)*trac2     - an2*an3 *trac3  - fe(1)
                            f(2) = divu*an3 +     -an3*an1 *trac1 +   -an3*an2 *trac2 + (1.-an3*an3)*trac3  - fe(2)         
                          resMax=0.
                          do n=0,numberOfEquations-1
                            resMax = max(resMax,abs(f(n)))
                          end do
                          if( resMax>resTol )then
                            write(*,'("Traction 3d curvilinear: ERROR residuals are large =",4(1pe12.4,1x))') (f(n),n=0,numberOfEquations-1)
                          else
                            ! write(*,'("Traction 3d curvilinear resMax=",1pe8.2)') resMax
                          end if
                        end if          
                      ! --- Extrapolate any extra ghost --- 
                      do ghost=3,numGhost
                        ! (j1,j2,j3) is the ghost point index
                        j1 = i1 - is1*ghost 
                        j2 = i2 - is2*ghost 
                        j3 = i3 - is3*ghost 
                        u(j1,j2,j3,u1c)=(3.*u(j1+is1,j2+is2,j3+is3,u1c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u1c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u1c))
                        u(j1,j2,j3,u2c)=(3.*u(j1+is1,j2+is2,j3+is3,u2c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u2c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u2c))
                          u(j1,j2,j3,u3c)=(3.*u(j1+is1,j2+is2,j3+is3,u3c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u3c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u3c))
                      end do
                 end do
                 end do
                 end do
            else if( addBoundaryForcing(side,axis)==1 .and. twilightZone==0 )then
                fe(0)=0.; fe(1)=0.;  fe(2)=0.; fe(3)=0.;   ! holds forcing 
                numberOfEquations=4;      ! number of ghost points we solve for   
                 do i3=n3a,n3b
                 do i2=n2a,n2b
                 do i1=n1a,n1b
                  ! here is the normal
                       an1 = rsxy(i1,i2,i3,axis,0)
                       an2 = rsxy(i1,i2,i3,axis,1)
                       an3 = rsxy(i1,i2,i3,axis,2)
                       aNormi = -is/max(epsx,sqrt(an1**2 + an2**2+ an3**2))
                       an1=an1*aNormi
                       an2=an2*aNormi
                       an3=an3*aNormi
                  ! *** FIX ME : USE THIS BELOW FOR RECTANGULAR CASE####################################
                  ! #If "2" eq "2"
                  !   #If "forcing" eq "twilightZone" 
                  !     #If "3" eq "2"
                  !       OGDERIV2D(0,1,0,0,i1,i2,i3,t,ux0,vx0)
                  !       OGDERIV2D(0,0,1,0,i1,i2,i3,t,uy0,vy0)
                  !       fe(0) = ux0 + vy0
                  !       fe(1) = 2*an1*an2*( ux0 - vy0 ) + ( an2**2 -an1**2 )*( uy0 + vx0 )
                  !     #Else
                  !       stop 333
                  !     #End
                  !   #End
                  ! #End
                    ! ----- TRACTION CURVILINEAR ----
                    ! --- here are  the equations we mean to satisfy ----
                    !  --> evaluate with the wrong values in the ghost  
                      ! --------- 2 2 -----------
                        ! --------------- TRACTION, 2=2, 3=3, CURVILINEAR -----------------
                        ! write(*,*) '--------------- TRACTION, order=2, dim=3, CURVILINEAR -----------------'          
                        numberOfEquations = 3 
                        ! Here are the realtive index positions of the 3 ghost we solve for 
                        shift1(0) = -is1; shift2(0) = -is2; shift3(0) = -is3; comp(0)=u1c; 
                        shift1(1) = -is1; shift2(1) = -is2; shift3(1) = -is3; comp(1)=u2c; 
                        shift1(2) = -is1; shift2(2) = -is2; shift3(2) = -is3; comp(2)=u3c; 
                        ! --- eval any forcing ---
                          !    (div(u)) n + mu* ( I - n n^T) epsm * n 
                        ! write(*,'("traction:3d: i1,i2,i3=",3i3," fe=",3(1pe12.4,1x))') i1,i2,i3,fe(0),fe(1),fe(2)
                        ! Fill in the matrix using the discrete delta approach
                          ! ! hw1 = half stencil width
                          ! hw1=orderOfAccuracy/2
                          ! hw2=hw1
                          ! if( nd.eq.2 )then
                          !   hw3=0
                          ! else
                          !   hw3=hw1
                          ! end if
                          !  write(*,'("EVAL-COEFF: i1,i2,i3=",3i3," hw1,hw2,hw3=",3i2)') i1,i2,i3,hw1,hw2,hw3
                          ! First eval equartions with no pertutbation --> save in f0 
                            !    (div(u)) n + mu* ( I - n n^T) epsm * n 
                            u1x = ux23(i1,i2,i3,u1c); u1y = uy23(i1,i2,i3,u1c); u1z = uz23(i1,i2,i3,u1c);
                            u2x = ux23(i1,i2,i3,u2c); u2y = uy23(i1,i2,i3,u2c); u2z = uz23(i1,i2,i3,u2c);
                            u3x = ux23(i1,i2,i3,u3c); u3y = uy23(i1,i2,i3,u3c); u3z = uz23(i1,i2,i3,u3c);
                            divu = u1x + u2y + u3z 
                            ! traction = epsm*n : 
                            trac1 = mu*( ( u1x + u1x )*an1 + ( u1y + u2x )*an2 + ( u1z + u3x )*an3 )
                            trac2 = mu*( ( u2x + u1y )*an1 + ( u2y + u2y )*an2 + ( u2z + u3y )*an3 )
                            trac3 = mu*( ( u3x + u1z )*an1 + ( u3y + u2z )*an2 + ( u3z + u3z )*an3 )
                            f(0) = divu*an1 +  (1.-an1*an1)*trac1 -    an1*an2 *trac2     - an1*an3 *trac3  - fe(0)
                            f(1) = divu*an2 +     -an2*an1 *trac1 +(1.-an2*an2)*trac2     - an2*an3 *trac3  - fe(1)
                            f(2) = divu*an3 +     -an3*an1 *trac1 +   -an3*an2 *trac2 + (1.-an3*an3)*trac3  - fe(2)         
                          do n1=0,numberOfEquations-1
                           f0(n1)=f(n1)
                          end do
                          delta=1.  ! perturb E by this amount 
                          do n2=0,numberOfEquations-1
                            ! pertub one component: 
                              u(i1+shift1(n2),i2+shift2(n2),i3+shift3(n2),comp(n2))=u(i1+shift1(n2),i2+shift2(n2),i3+shift3(n2),comp(n2))+(delta)
                              !    (div(u)) n + mu* ( I - n n^T) epsm * n 
                              u1x = ux23(i1,i2,i3,u1c); u1y = uy23(i1,i2,i3,u1c); u1z = uz23(i1,i2,i3,u1c);
                              u2x = ux23(i1,i2,i3,u2c); u2y = uy23(i1,i2,i3,u2c); u2z = uz23(i1,i2,i3,u2c);
                              u3x = ux23(i1,i2,i3,u3c); u3y = uy23(i1,i2,i3,u3c); u3z = uz23(i1,i2,i3,u3c);
                              divu = u1x + u2y + u3z 
                              ! traction = epsm*n : 
                              trac1 = mu*( ( u1x + u1x )*an1 + ( u1y + u2x )*an2 + ( u1z + u3x )*an3 )
                              trac2 = mu*( ( u2x + u1y )*an1 + ( u2y + u2y )*an2 + ( u2z + u3y )*an3 )
                              trac3 = mu*( ( u3x + u1z )*an1 + ( u3y + u2z )*an2 + ( u3z + u3z )*an3 )
                              f(0) = divu*an1 +  (1.-an1*an1)*trac1 -    an1*an2 *trac2     - an1*an3 *trac3  - fe(0)
                              f(1) = divu*an2 +     -an2*an1 *trac1 +(1.-an2*an2)*trac2     - an2*an3 *trac3  - fe(1)
                              f(2) = divu*an3 +     -an3*an1 *trac1 +   -an3*an2 *trac2 + (1.-an3*an3)*trac3  - fe(2)         
                            ! compute the difference
                            do n1=0,numberOfEquations-1
                              f(n1)=f(n1)-f0(n1)
                              a3(n1,n2) = f(n1)
                            end do
                            ! reset pertubation
                              u(i1+shift1(n2),i2+shift2(n2),i3+shift3(n2),comp(n2))=u(i1+shift1(n2),i2+shift2(n2),i3+shift3(n2),comp(n2))+(-delta)
                          end do 
                          ! restore f -- not needed here
                          ! evalTractionCurvilinear3d()
                        ! print the matrix: 
                        ! write(*,'("traction 3d curv: i1,i2=",2i3,/," a3=[",3(1pe12.4,1x),"; ...")') i1,i2,((a3(n,n2),n2=0,2),n=0,2) 
                        ! evaluate equations with wrong values in the ghost        
                          !    (div(u)) n + mu* ( I - n n^T) epsm * n 
                          u1x = ux23(i1,i2,i3,u1c); u1y = uy23(i1,i2,i3,u1c); u1z = uz23(i1,i2,i3,u1c);
                          u2x = ux23(i1,i2,i3,u2c); u2y = uy23(i1,i2,i3,u2c); u2z = uz23(i1,i2,i3,u2c);
                          u3x = ux23(i1,i2,i3,u3c); u3y = uy23(i1,i2,i3,u3c); u3z = uz23(i1,i2,i3,u3c);
                          divu = u1x + u2y + u3z 
                          ! traction = epsm*n : 
                          trac1 = mu*( ( u1x + u1x )*an1 + ( u1y + u2x )*an2 + ( u1z + u3x )*an3 )
                          trac2 = mu*( ( u2x + u1y )*an1 + ( u2y + u2y )*an2 + ( u2z + u3y )*an3 )
                          trac3 = mu*( ( u3x + u1z )*an1 + ( u3y + u2z )*an2 + ( u3z + u3z )*an3 )
                          f(0) = divu*an1 +  (1.-an1*an1)*trac1 -    an1*an2 *trac2     - an1*an3 *trac3  - fe(0)
                          f(1) = divu*an2 +     -an2*an1 *trac1 +(1.-an2*an2)*trac2     - an2*an3 *trac3  - fe(1)
                          f(2) = divu*an3 +     -an3*an1 *trac1 +   -an3*an2 *trac2 + (1.-an3*an3)*trac3  - fe(2)         
                        ! adjust rhs for wrong values on ghost 
                        do n=0,numberOfEquations-1
                          f(n) = -f(n) 
                          do n2=0,numberOfEquations-1
                            f(n) = f(n) + a3(n,n2)*u(i1+shift1(n2),i2+shift2(n2),i3+shift3(n2),comp(n2))
                          end do
                        end do
                        ! solve for ghost
                        call dgeco( a3(0,0), numberOfEquations, numberOfEquations, ipvt(0),rcond,work(0) )
                        call dgesl( a3(0,0), numberOfEquations, numberOfEquations, ipvt(0), f(0), job )
                        ! write(*,'(" traction curvilinear 3d: rcond=",1pe10.2)') rcond
                        ! assign values 
                        do n2=0,numberOfEquations-1
                          u(i1+shift1(n2),i2+shift2(n2),i3+shift3(n2),comp(n2)) = f(n2)
                        end do    
                        ! call ogDeriv(ep,0,0,0,0,xy(i1-is1,i2-is2,i3-is3,0),xy(i1-is1,i2-is2,i3-is3,1),xy(i1-is1,i2-is2,i3-is3,2),t,u1c,u1e)
                        ! call ogDeriv(ep,0,0,0,0,xy(i1-is1,i2-is2,i3-is3,0),xy(i1-is1,i2-is2,i3-is3,1),xy(i1-is1,i2-is2,i3-is3,2),t,u2c,u2e)
                        ! write(*,'("traction:3d: i1,i2,i3=",3i3," u1(-1)=",1pe12.4," err=",1pe10.2," u2(-1)=",1pe12.4," err=",1pe10.2)') !        i1,i2,i3,f(0),f(0)-u1e,f(1),f(1)-u2e           
                        ! ---- check residuals ---   
                        if( checkResiduals )then
                          ! re-evaluate the equations
                            !    (div(u)) n + mu* ( I - n n^T) epsm * n 
                            u1x = ux23(i1,i2,i3,u1c); u1y = uy23(i1,i2,i3,u1c); u1z = uz23(i1,i2,i3,u1c);
                            u2x = ux23(i1,i2,i3,u2c); u2y = uy23(i1,i2,i3,u2c); u2z = uz23(i1,i2,i3,u2c);
                            u3x = ux23(i1,i2,i3,u3c); u3y = uy23(i1,i2,i3,u3c); u3z = uz23(i1,i2,i3,u3c);
                            divu = u1x + u2y + u3z 
                            ! traction = epsm*n : 
                            trac1 = mu*( ( u1x + u1x )*an1 + ( u1y + u2x )*an2 + ( u1z + u3x )*an3 )
                            trac2 = mu*( ( u2x + u1y )*an1 + ( u2y + u2y )*an2 + ( u2z + u3y )*an3 )
                            trac3 = mu*( ( u3x + u1z )*an1 + ( u3y + u2z )*an2 + ( u3z + u3z )*an3 )
                            f(0) = divu*an1 +  (1.-an1*an1)*trac1 -    an1*an2 *trac2     - an1*an3 *trac3  - fe(0)
                            f(1) = divu*an2 +     -an2*an1 *trac1 +(1.-an2*an2)*trac2     - an2*an3 *trac3  - fe(1)
                            f(2) = divu*an3 +     -an3*an1 *trac1 +   -an3*an2 *trac2 + (1.-an3*an3)*trac3  - fe(2)         
                          resMax=0.
                          do n=0,numberOfEquations-1
                            resMax = max(resMax,abs(f(n)))
                          end do
                          if( resMax>resTol )then
                            write(*,'("Traction 3d curvilinear: ERROR residuals are large =",4(1pe12.4,1x))') (f(n),n=0,numberOfEquations-1)
                          else
                            ! write(*,'("Traction 3d curvilinear resMax=",1pe8.2)') resMax
                          end if
                        end if          
                      ! --- Extrapolate any extra ghost --- 
                      do ghost=3,numGhost
                        ! (j1,j2,j3) is the ghost point index
                        j1 = i1 - is1*ghost 
                        j2 = i2 - is2*ghost 
                        j3 = i3 - is3*ghost 
                        u(j1,j2,j3,u1c)=(3.*u(j1+is1,j2+is2,j3+is3,u1c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u1c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u1c))
                        u(j1,j2,j3,u2c)=(3.*u(j1+is1,j2+is2,j3+is3,u2c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u2c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u2c))
                          u(j1,j2,j3,u3c)=(3.*u(j1+is1,j2+is2,j3+is3,u3c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u3c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u3c))
                      end do
                 end do
                 end do
                 end do
            else
                fe(0)=0.; fe(1)=0.;  fe(2)=0.; fe(3)=0.;   ! holds forcing 
                numberOfEquations=4;      ! number of ghost points we solve for   
                 do i3=n3a,n3b
                 do i2=n2a,n2b
                 do i1=n1a,n1b
                  ! here is the normal
                       an1 = rsxy(i1,i2,i3,axis,0)
                       an2 = rsxy(i1,i2,i3,axis,1)
                       an3 = rsxy(i1,i2,i3,axis,2)
                       aNormi = -is/max(epsx,sqrt(an1**2 + an2**2+ an3**2))
                       an1=an1*aNormi
                       an2=an2*aNormi
                       an3=an3*aNormi
                  ! *** FIX ME : USE THIS BELOW FOR RECTANGULAR CASE####################################
                  ! #If "2" eq "2"
                  !   #If "twilightZone" eq "twilightZone" 
                  !     #If "3" eq "2"
                  !       OGDERIV2D(0,1,0,0,i1,i2,i3,t,ux0,vx0)
                  !       OGDERIV2D(0,0,1,0,i1,i2,i3,t,uy0,vy0)
                  !       fe(0) = ux0 + vy0
                  !       fe(1) = 2*an1*an2*( ux0 - vy0 ) + ( an2**2 -an1**2 )*( uy0 + vx0 )
                  !     #Else
                  !       stop 333
                  !     #End
                  !   #End
                  ! #End
                      call ogDeriv(ep,0,1,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,u1c,ux0)
                      call ogDeriv(ep,0,0,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,u1c,uy0)
                      call ogDeriv(ep,0,0,0,1,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,u1c,uz0)
                      call ogDeriv(ep,0,1,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,u2c,vx0)
                      call ogDeriv(ep,0,0,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,u2c,vy0)
                      call ogDeriv(ep,0,0,0,1,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,u2c,vz0) 
                      call ogDeriv(ep,0,1,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,u3c,wx0)
                      call ogDeriv(ep,0,0,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,u3c,wy0)
                      call ogDeriv(ep,0,0,0,1,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,u3c,wz0)               
                    ! ----- TRACTION CURVILINEAR ----
                    ! --- here are  the equations we mean to satisfy ----
                    !  --> evaluate with the wrong values in the ghost  
                      ! --------- 2 2 -----------
                        ! --------------- TRACTION, 2=2, 3=3, CURVILINEAR -----------------
                        ! write(*,*) '--------------- TRACTION, order=2, dim=3, CURVILINEAR -----------------'          
                        numberOfEquations = 3 
                        ! Here are the realtive index positions of the 3 ghost we solve for 
                        shift1(0) = -is1; shift2(0) = -is2; shift3(0) = -is3; comp(0)=u1c; 
                        shift1(1) = -is1; shift2(1) = -is2; shift3(1) = -is3; comp(1)=u2c; 
                        shift1(2) = -is1; shift2(2) = -is2; shift3(2) = -is3; comp(2)=u3c; 
                        ! --- eval any forcing ---
                          !    (div(u)) n + mu* ( I - n n^T) epsm * n 
                            call ogDeriv(ep,0,1,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,u1c,u1x)
                            call ogDeriv(ep,0,0,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,u1c,u1y)
                            call ogDeriv(ep,0,0,0,1,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,u1c,u1z)
                            call ogDeriv(ep,0,1,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,u2c,u2x)
                            call ogDeriv(ep,0,0,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,u2c,u2y)
                            call ogDeriv(ep,0,0,0,1,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,u2c,u2z)
                            call ogDeriv(ep,0,1,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,u3c,u3x)
                            call ogDeriv(ep,0,0,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,u3c,u3y)
                            call ogDeriv(ep,0,0,0,1,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,u3c,u3z)        
                            divu = u1x + u2y + u3z 
                            trac1 = mu*( ( u1x + u1x )*an1 + ( u1y + u2x )*an2 + ( u1z + u3x )*an3 )
                            trac2 = mu*( ( u2x + u1y )*an1 + ( u2y + u2y )*an2 + ( u2z + u3y )*an3 )
                            trac3 = mu*( ( u3x + u1z )*an1 + ( u3y + u2z )*an2 + ( u3z + u3z )*an3 )
                            fe(0) = divu*an1 +  (1.-an1*an1)*trac1 -    an1*an2 *trac2     - an1*an3 *trac3  
                            fe(1) = divu*an2 +     -an2*an1 *trac1 +(1.-an2*an2)*trac2     - an2*an3 *trac3  
                            fe(2) = divu*an3 +     -an3*an1 *trac1 +   -an3*an2 *trac2 + (1.-an3*an3)*trac3      
                        ! write(*,'("traction:3d: i1,i2,i3=",3i3," fe=",3(1pe12.4,1x))') i1,i2,i3,fe(0),fe(1),fe(2)
                        ! Fill in the matrix using the discrete delta approach
                          ! ! hw1 = half stencil width
                          ! hw1=orderOfAccuracy/2
                          ! hw2=hw1
                          ! if( nd.eq.2 )then
                          !   hw3=0
                          ! else
                          !   hw3=hw1
                          ! end if
                          !  write(*,'("EVAL-COEFF: i1,i2,i3=",3i3," hw1,hw2,hw3=",3i2)') i1,i2,i3,hw1,hw2,hw3
                          ! First eval equartions with no pertutbation --> save in f0 
                            !    (div(u)) n + mu* ( I - n n^T) epsm * n 
                            u1x = ux23(i1,i2,i3,u1c); u1y = uy23(i1,i2,i3,u1c); u1z = uz23(i1,i2,i3,u1c);
                            u2x = ux23(i1,i2,i3,u2c); u2y = uy23(i1,i2,i3,u2c); u2z = uz23(i1,i2,i3,u2c);
                            u3x = ux23(i1,i2,i3,u3c); u3y = uy23(i1,i2,i3,u3c); u3z = uz23(i1,i2,i3,u3c);
                            divu = u1x + u2y + u3z 
                            ! traction = epsm*n : 
                            trac1 = mu*( ( u1x + u1x )*an1 + ( u1y + u2x )*an2 + ( u1z + u3x )*an3 )
                            trac2 = mu*( ( u2x + u1y )*an1 + ( u2y + u2y )*an2 + ( u2z + u3y )*an3 )
                            trac3 = mu*( ( u3x + u1z )*an1 + ( u3y + u2z )*an2 + ( u3z + u3z )*an3 )
                            f(0) = divu*an1 +  (1.-an1*an1)*trac1 -    an1*an2 *trac2     - an1*an3 *trac3  - fe(0)
                            f(1) = divu*an2 +     -an2*an1 *trac1 +(1.-an2*an2)*trac2     - an2*an3 *trac3  - fe(1)
                            f(2) = divu*an3 +     -an3*an1 *trac1 +   -an3*an2 *trac2 + (1.-an3*an3)*trac3  - fe(2)         
                          do n1=0,numberOfEquations-1
                           f0(n1)=f(n1)
                          end do
                          delta=1.  ! perturb E by this amount 
                          do n2=0,numberOfEquations-1
                            ! pertub one component: 
                              u(i1+shift1(n2),i2+shift2(n2),i3+shift3(n2),comp(n2))=u(i1+shift1(n2),i2+shift2(n2),i3+shift3(n2),comp(n2))+(delta)
                              !    (div(u)) n + mu* ( I - n n^T) epsm * n 
                              u1x = ux23(i1,i2,i3,u1c); u1y = uy23(i1,i2,i3,u1c); u1z = uz23(i1,i2,i3,u1c);
                              u2x = ux23(i1,i2,i3,u2c); u2y = uy23(i1,i2,i3,u2c); u2z = uz23(i1,i2,i3,u2c);
                              u3x = ux23(i1,i2,i3,u3c); u3y = uy23(i1,i2,i3,u3c); u3z = uz23(i1,i2,i3,u3c);
                              divu = u1x + u2y + u3z 
                              ! traction = epsm*n : 
                              trac1 = mu*( ( u1x + u1x )*an1 + ( u1y + u2x )*an2 + ( u1z + u3x )*an3 )
                              trac2 = mu*( ( u2x + u1y )*an1 + ( u2y + u2y )*an2 + ( u2z + u3y )*an3 )
                              trac3 = mu*( ( u3x + u1z )*an1 + ( u3y + u2z )*an2 + ( u3z + u3z )*an3 )
                              f(0) = divu*an1 +  (1.-an1*an1)*trac1 -    an1*an2 *trac2     - an1*an3 *trac3  - fe(0)
                              f(1) = divu*an2 +     -an2*an1 *trac1 +(1.-an2*an2)*trac2     - an2*an3 *trac3  - fe(1)
                              f(2) = divu*an3 +     -an3*an1 *trac1 +   -an3*an2 *trac2 + (1.-an3*an3)*trac3  - fe(2)         
                            ! compute the difference
                            do n1=0,numberOfEquations-1
                              f(n1)=f(n1)-f0(n1)
                              a3(n1,n2) = f(n1)
                            end do
                            ! reset pertubation
                              u(i1+shift1(n2),i2+shift2(n2),i3+shift3(n2),comp(n2))=u(i1+shift1(n2),i2+shift2(n2),i3+shift3(n2),comp(n2))+(-delta)
                          end do 
                          ! restore f -- not needed here
                          ! evalTractionCurvilinear3d()
                        ! print the matrix: 
                        ! write(*,'("traction 3d curv: i1,i2=",2i3,/," a3=[",3(1pe12.4,1x),"; ...")') i1,i2,((a3(n,n2),n2=0,2),n=0,2) 
                        ! evaluate equations with wrong values in the ghost        
                          !    (div(u)) n + mu* ( I - n n^T) epsm * n 
                          u1x = ux23(i1,i2,i3,u1c); u1y = uy23(i1,i2,i3,u1c); u1z = uz23(i1,i2,i3,u1c);
                          u2x = ux23(i1,i2,i3,u2c); u2y = uy23(i1,i2,i3,u2c); u2z = uz23(i1,i2,i3,u2c);
                          u3x = ux23(i1,i2,i3,u3c); u3y = uy23(i1,i2,i3,u3c); u3z = uz23(i1,i2,i3,u3c);
                          divu = u1x + u2y + u3z 
                          ! traction = epsm*n : 
                          trac1 = mu*( ( u1x + u1x )*an1 + ( u1y + u2x )*an2 + ( u1z + u3x )*an3 )
                          trac2 = mu*( ( u2x + u1y )*an1 + ( u2y + u2y )*an2 + ( u2z + u3y )*an3 )
                          trac3 = mu*( ( u3x + u1z )*an1 + ( u3y + u2z )*an2 + ( u3z + u3z )*an3 )
                          f(0) = divu*an1 +  (1.-an1*an1)*trac1 -    an1*an2 *trac2     - an1*an3 *trac3  - fe(0)
                          f(1) = divu*an2 +     -an2*an1 *trac1 +(1.-an2*an2)*trac2     - an2*an3 *trac3  - fe(1)
                          f(2) = divu*an3 +     -an3*an1 *trac1 +   -an3*an2 *trac2 + (1.-an3*an3)*trac3  - fe(2)         
                        ! adjust rhs for wrong values on ghost 
                        do n=0,numberOfEquations-1
                          f(n) = -f(n) 
                          do n2=0,numberOfEquations-1
                            f(n) = f(n) + a3(n,n2)*u(i1+shift1(n2),i2+shift2(n2),i3+shift3(n2),comp(n2))
                          end do
                        end do
                        ! solve for ghost
                        call dgeco( a3(0,0), numberOfEquations, numberOfEquations, ipvt(0),rcond,work(0) )
                        call dgesl( a3(0,0), numberOfEquations, numberOfEquations, ipvt(0), f(0), job )
                        ! write(*,'(" traction curvilinear 3d: rcond=",1pe10.2)') rcond
                        ! assign values 
                        do n2=0,numberOfEquations-1
                          u(i1+shift1(n2),i2+shift2(n2),i3+shift3(n2),comp(n2)) = f(n2)
                        end do    
                        ! call ogDeriv(ep,0,0,0,0,xy(i1-is1,i2-is2,i3-is3,0),xy(i1-is1,i2-is2,i3-is3,1),xy(i1-is1,i2-is2,i3-is3,2),t,u1c,u1e)
                        ! call ogDeriv(ep,0,0,0,0,xy(i1-is1,i2-is2,i3-is3,0),xy(i1-is1,i2-is2,i3-is3,1),xy(i1-is1,i2-is2,i3-is3,2),t,u2c,u2e)
                        ! write(*,'("traction:3d: i1,i2,i3=",3i3," u1(-1)=",1pe12.4," err=",1pe10.2," u2(-1)=",1pe12.4," err=",1pe10.2)') !        i1,i2,i3,f(0),f(0)-u1e,f(1),f(1)-u2e           
                        ! ---- check residuals ---   
                        if( checkResiduals )then
                          ! re-evaluate the equations
                            !    (div(u)) n + mu* ( I - n n^T) epsm * n 
                            u1x = ux23(i1,i2,i3,u1c); u1y = uy23(i1,i2,i3,u1c); u1z = uz23(i1,i2,i3,u1c);
                            u2x = ux23(i1,i2,i3,u2c); u2y = uy23(i1,i2,i3,u2c); u2z = uz23(i1,i2,i3,u2c);
                            u3x = ux23(i1,i2,i3,u3c); u3y = uy23(i1,i2,i3,u3c); u3z = uz23(i1,i2,i3,u3c);
                            divu = u1x + u2y + u3z 
                            ! traction = epsm*n : 
                            trac1 = mu*( ( u1x + u1x )*an1 + ( u1y + u2x )*an2 + ( u1z + u3x )*an3 )
                            trac2 = mu*( ( u2x + u1y )*an1 + ( u2y + u2y )*an2 + ( u2z + u3y )*an3 )
                            trac3 = mu*( ( u3x + u1z )*an1 + ( u3y + u2z )*an2 + ( u3z + u3z )*an3 )
                            f(0) = divu*an1 +  (1.-an1*an1)*trac1 -    an1*an2 *trac2     - an1*an3 *trac3  - fe(0)
                            f(1) = divu*an2 +     -an2*an1 *trac1 +(1.-an2*an2)*trac2     - an2*an3 *trac3  - fe(1)
                            f(2) = divu*an3 +     -an3*an1 *trac1 +   -an3*an2 *trac2 + (1.-an3*an3)*trac3  - fe(2)         
                          resMax=0.
                          do n=0,numberOfEquations-1
                            resMax = max(resMax,abs(f(n)))
                          end do
                          if( resMax>resTol )then
                            write(*,'("Traction 3d curvilinear: ERROR residuals are large =",4(1pe12.4,1x))') (f(n),n=0,numberOfEquations-1)
                          else
                            ! write(*,'("Traction 3d curvilinear resMax=",1pe8.2)') resMax
                          end if
                        end if          
                      ! --- Extrapolate any extra ghost --- 
                      do ghost=3,numGhost
                        ! (j1,j2,j3) is the ghost point index
                        j1 = i1 - is1*ghost 
                        j2 = i2 - is2*ghost 
                        j3 = i3 - is3*ghost 
                        u(j1,j2,j3,u1c)=(3.*u(j1+is1,j2+is2,j3+is3,u1c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u1c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u1c))
                        u(j1,j2,j3,u2c)=(3.*u(j1+is1,j2+is2,j3+is3,u2c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u2c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u2c))
                          u(j1,j2,j3,u3c)=(3.*u(j1+is1,j2+is2,j3+is3,u3c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u3c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u3c))
                      end do
                 end do
                 end do
                 end do
            end if
          else if( boundaryCondition(side,axis).eq.slipWall )then 
            write(*,*) 'bcOptIsm: slipWall - finish me'
            stop 1234
            alpha=lambda/(lambda+2.*mu)
            beta =1./(lambda+2.*mu)
            if( axis.eq.0 )then
              ! v.x = -u.y       
             if( addBoundaryForcing(side,axis).eq.0 )then
              ! no forcing 
               i3=n3a
               do i2=n2a,n2b
               do i1=n1a,n1b
               if( mask(i1,i2,i3).gt.0 )then
               u(i1-is1,i2-is2,i3,u2c)=u(i1+is1,i2+is2,i3,u2c)+is1*dx(0)*2.*      uy22r(i1,i2,i3,u1c)
               end if
               end do
               end do
             else if( assignTwilightZone.eq.0 )then
               ! include forcing terms 
               i3=n3a
               do i2=n2a,n2b
               do i1=n1a,n1b
               if( mask(i1,i2,i3).gt.0 )then
               u(i1-is1,i2-is2,i3,u2c)=u(i1+is1,i2+is2,i3,u2c)+dx(0)*2.*(is1*uy22r(i1,i2,i3,u1c)    + (1./mu)*bcf(side,axis,i1,i2,i3,u2c) )
               end if
               end do
               end do
             else 
              ! Twilight-zone: 
              ! u.x = -alpha*v.y + ue.x -alpha*ve.y 
               i3=n3a
               do i2=n2a,n2b
               do i1=n1a,n1b
               if( mask(i1,i2,i3).gt.0 )then
                 call ogDeriv2(ep, 0,1,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, u1c,ux0, u2c,vx0)
                 call ogDeriv2(ep, 0,0,1,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, u1c,uy0, u2c,vy0)
               u(i1-is1,i2-is2,i3,u2c)=u(i1+is1,i2+is2,i3,u2c)-is1*dx(0)*2.*(-uy22r(i1,i2,i3,u1c) + vx0 + uy0 )
               !     write(*,'("i1,i2=",2i3," ux0,vx0,uy0,vy0=",4e10.2)') i1,i2, ux0,vx0,uy0,vy0                
               end if
               end do
               end do
             end if
            else if( axis.eq.1 )then
            ! u.y = - v.x
             if( addBoundaryForcing(side,axis).eq.0 )then
               i3=n3a
               do i2=n2a,n2b
               do i1=n1a,n1b
               if( mask(i1,i2,i3).gt.0 )then
               u(i1-is1,i2-is2,i3,u1c)=u(i1+is1,i2+is2,i3,u1c)+is2*dx(1)*2.*      ux22r(i1,i2,i3,u2c)
               end if
               end do
               end do
             else if( assignTwilightZone.eq.0 )then
               ! include forcing terms
               i3=n3a
               do i2=n2a,n2b
               do i1=n1a,n1b
               if( mask(i1,i2,i3).gt.0 )then
               u(i1-is1,i2-is2,i3,u1c)=u(i1+is1,i2+is2,i3,u1c)+dx(1)*2.*(is2*ux22r(i1,i2,i3,u2c)    + (1./mu)*bcf(side,axis,i1,i2,i3,u1c) )
               end if
               end do
               end do
             else
              ! Twilight-zone: 
               i3=n3a
               do i2=n2a,n2b
               do i1=n1a,n1b
               if( mask(i1,i2,i3).gt.0 )then
                 call ogDeriv2(ep, 0,1,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, u1c,ux0, u2c,vx0)
                 call ogDeriv2(ep, 0,0,1,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, u1c,uy0, u2c,vy0)
               u(i1-is1,i2-is2,i3,u1c)=u(i1+is1,i2+is2,i3,u1c)-is2*dx(1)*2.*(      -ux22r(i1,i2,i3,u2c)+uy0+vx0)
               end if
               end do
               end do
             end if
            end if   
          end if ! end bc 
          end do ! end side
          end do ! end axis
         ! ---------------------------------------------------------------------
         ! ----------- STAGE 3, ASSIGN CORNERS AND EDGES -----------------------
         ! ---------------------------------------------------------------------
           if( twilightZone==1 )then
               ! numberOfGhostPoints=orderOfAccuracy/2
               ! ----- Assign the edges -------
               ! write(*,*) ' finish me: assign edges in 3D'
                 do edgeDirection=0,2 ! direction parallel to the edge
                   do sidea=0,1
                   do sideb=0,1
                     if( edgeDirection.eq.0 )then
                       side1=0
                       side2=sidea
                       side3=sideb
                     else if( edgeDirection.eq.1 )then
                       side1=sideb 
                       side2=0
                       side3=sidea
                     else
                       side1=sidea
                       side2=sideb
                       side3=0
                     end if
                     is1=1-2*(side1)
                     is2=1-2*(side2)
                     is3=1-2*(side3)
                     if( edgeDirection.eq.2 )then
                      is3=0
                      n1a=gridIndexRange(side1,0)
                      n1b=gridIndexRange(side1,0)
                      n2a=gridIndexRange(side2,1)
                      n2b=gridIndexRange(side2,1)
                      n3a=gridIndexRange(0,2)
                      n3b=gridIndexRange(1,2)
                      bc1=boundaryCondition(side1,0)
                      bc2=boundaryCondition(side2,1)
                     else if( edgeDirection.eq.1 )then
                      is2=0
                      n1a=gridIndexRange(side1,0)
                      n1b=gridIndexRange(side1,0)
                      n2a=gridIndexRange(    0,1)
                      n2b=gridIndexRange(    1,1)
                      n3a=gridIndexRange(side3,2)
                      n3b=gridIndexRange(side3,2)
                      bc1=boundaryCondition(side1,0)
                      bc2=boundaryCondition(side3,2)
                     else 
                      is1=0  
                      n1a=gridIndexRange(    0,0)
                      n1b=gridIndexRange(    1,0)
                      n2a=gridIndexRange(side2,1)
                      n2b=gridIndexRange(side2,1)
                      n3a=gridIndexRange(side3,2)
                      n3b=gridIndexRange(side3,2)
                      bc1=boundaryCondition(side2,1)
                      bc2=boundaryCondition(side3,2)
                     end if
                     ! ************ assign corner points outside edges ***********************
                     if( bc1.gt.0 .and. bc2.gt.0 )then
                       ! write(*,'(" assign edges: edgeDirection=",i2," sidea,sideb=",2i3)') edgeDirection,sidea,sideb
                       do ghost1=1,numGhost
                       do ghost2=1,numGhost
                         ! shift to ghost point "(m1,m2)"
                         if( edgeDirection.eq.2 )then 
                           js1=is1*ghost1  
                           js2=is2*ghost2
                           js3=0
                         else if( edgeDirection.eq.1 )then 
                           js1=is1*ghost1  
                           js2=0
                           js3=is3*ghost2
                         else 
                           js1=0
                           js2=is2*ghost1
                           js3=is3*ghost2
                         end if      
                          do i3=n3a,n3b
                          do i2=n2a,n2b
                          do i1=n1a,n1b
                           j1 = i1-js1; j2=i2-js2; j3=i3-js3;  ! ghost point 
                           u(j1,j2,j3,u1c) = (3.*u(j1+is1,j2+is2,j3+is3,u1c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u1c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u1c))
                           u(j1,j2,j3,u2c) = (3.*u(j1+is1,j2+is2,j3+is3,u2c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u2c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u2c))
                           u(j1,j2,j3,u3c) = (3.*u(j1+is1,j2+is2,j3+is3,u3c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u3c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u3c))
                          end do
                          end do
                          end do
                        end do ! ghost1
                        end do ! ghost2
                     end if ! end if bc1>0 and bc2>0
                   end do ! end do sideb
                   end do ! end do sidea
                 end do ! end do edgeDirection
               do side3=0,1
               do side2=0,1
               do side1=0,1
                 ! ----- assign ghost values outside the corner (vertex) -----
                 i1=gridIndexRange(side1,0)
                 i2=gridIndexRange(side2,1)
                 i3=gridIndexRange(side3,2)
                 is1=1-2*side1
                 is2=1-2*side2
                 is3=1-2*side3
                 if( boundaryCondition(side1,0).gt.0 .and.boundaryCondition(side2,1).gt.0 .and.boundaryCondition(side3,2).gt.0 )then
                  do ghost=1,numGhost
                    j1 = i1 - ghost*is1
                    j2 = i2 - ghost*is2
                    j3 = i3 - ghost*is3
                    u(j1,j2,j3,u1c) = (3.*u(j1+is1,j2+is2,j3+is3,u1c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u1c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u1c)) 
                    u(j1,j2,j3,u2c) = (3.*u(j1+is1,j2+is2,j3+is3,u2c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u2c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u2c)) 
                    u(j1,j2,j3,u3c) = (3.*u(j1+is1,j2+is2,j3+is3,u3c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u3c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u3c)) 
                  end do
                 end if
               end do
               end do
               end do
           else
               ! numberOfGhostPoints=orderOfAccuracy/2
               ! ----- Assign the edges -------
               ! write(*,*) ' finish me: assign edges in 3D'
                 do edgeDirection=0,2 ! direction parallel to the edge
                   do sidea=0,1
                   do sideb=0,1
                     if( edgeDirection.eq.0 )then
                       side1=0
                       side2=sidea
                       side3=sideb
                     else if( edgeDirection.eq.1 )then
                       side1=sideb 
                       side2=0
                       side3=sidea
                     else
                       side1=sidea
                       side2=sideb
                       side3=0
                     end if
                     is1=1-2*(side1)
                     is2=1-2*(side2)
                     is3=1-2*(side3)
                     if( edgeDirection.eq.2 )then
                      is3=0
                      n1a=gridIndexRange(side1,0)
                      n1b=gridIndexRange(side1,0)
                      n2a=gridIndexRange(side2,1)
                      n2b=gridIndexRange(side2,1)
                      n3a=gridIndexRange(0,2)
                      n3b=gridIndexRange(1,2)
                      bc1=boundaryCondition(side1,0)
                      bc2=boundaryCondition(side2,1)
                     else if( edgeDirection.eq.1 )then
                      is2=0
                      n1a=gridIndexRange(side1,0)
                      n1b=gridIndexRange(side1,0)
                      n2a=gridIndexRange(    0,1)
                      n2b=gridIndexRange(    1,1)
                      n3a=gridIndexRange(side3,2)
                      n3b=gridIndexRange(side3,2)
                      bc1=boundaryCondition(side1,0)
                      bc2=boundaryCondition(side3,2)
                     else 
                      is1=0  
                      n1a=gridIndexRange(    0,0)
                      n1b=gridIndexRange(    1,0)
                      n2a=gridIndexRange(side2,1)
                      n2b=gridIndexRange(side2,1)
                      n3a=gridIndexRange(side3,2)
                      n3b=gridIndexRange(side3,2)
                      bc1=boundaryCondition(side2,1)
                      bc2=boundaryCondition(side3,2)
                     end if
                     ! ************ assign corner points outside edges ***********************
                     if( bc1.gt.0 .and. bc2.gt.0 )then
                       ! write(*,'(" assign edges: edgeDirection=",i2," sidea,sideb=",2i3)') edgeDirection,sidea,sideb
                       do ghost1=1,numGhost
                       do ghost2=1,numGhost
                         ! shift to ghost point "(m1,m2)"
                         if( edgeDirection.eq.2 )then 
                           js1=is1*ghost1  
                           js2=is2*ghost2
                           js3=0
                         else if( edgeDirection.eq.1 )then 
                           js1=is1*ghost1  
                           js2=0
                           js3=is3*ghost2
                         else 
                           js1=0
                           js2=is2*ghost1
                           js3=is3*ghost2
                         end if      
                          do i3=n3a,n3b
                          do i2=n2a,n2b
                          do i1=n1a,n1b
                           j1 = i1-js1; j2=i2-js2; j3=i3-js3;  ! ghost point 
                           u(j1,j2,j3,u1c) = (3.*u(j1+is1,j2+is2,j3+is3,u1c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u1c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u1c))
                           u(j1,j2,j3,u2c) = (3.*u(j1+is1,j2+is2,j3+is3,u2c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u2c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u2c))
                           u(j1,j2,j3,u3c) = (3.*u(j1+is1,j2+is2,j3+is3,u3c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u3c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u3c))
                          end do
                          end do
                          end do
                        end do ! ghost1
                        end do ! ghost2
                     end if ! end if bc1>0 and bc2>0
                   end do ! end do sideb
                   end do ! end do sidea
                 end do ! end do edgeDirection
               do side3=0,1
               do side2=0,1
               do side1=0,1
                 ! ----- assign ghost values outside the corner (vertex) -----
                 i1=gridIndexRange(side1,0)
                 i2=gridIndexRange(side2,1)
                 i3=gridIndexRange(side3,2)
                 is1=1-2*side1
                 is2=1-2*side2
                 is3=1-2*side3
                 if( boundaryCondition(side1,0).gt.0 .and.boundaryCondition(side2,1).gt.0 .and.boundaryCondition(side3,2).gt.0 )then
                  do ghost=1,numGhost
                    j1 = i1 - ghost*is1
                    j2 = i2 - ghost*is2
                    j3 = i3 - ghost*is3
                    u(j1,j2,j3,u1c) = (3.*u(j1+is1,j2+is2,j3+is3,u1c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u1c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u1c)) 
                    u(j1,j2,j3,u2c) = (3.*u(j1+is1,j2+is2,j3+is3,u2c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u2c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u2c)) 
                    u(j1,j2,j3,u3c) = (3.*u(j1+is1,j2+is2,j3+is3,u3c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u3c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u3c)) 
                  end do
                 end if
               end do
               end do
               end do
           end if    
         return
         end
