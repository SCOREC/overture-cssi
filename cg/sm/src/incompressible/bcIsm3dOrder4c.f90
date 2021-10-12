! This file automatically generated from bcOptIsm.bf90 with bpp.
        subroutine bcIsm3dOrder4c( nd,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,gridIndexRange,u,mask,rsxy,xy,ndMatProp,matIndex,matValpc,matVal,boundaryCondition,addBoundaryForcing,interfaceType,dim,bcf00,bcf10,bcf01,bcf11,bcf02,bcf12,bcf0,bcOffset,ipar,rpar,ierr )
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
             real dx(0:2),dr(0:2)
             real t,ep,dt,c1,c2,rho,mu,lambda,alpha,beta
             integer axisp1,axisp2,i1,i2,i3,is1,is2,is3,j1,j2,j3,js1,js2,js3,ks1,ks2,ks3,is,js,it,nit
             integer option,initialized
             integer ghost,numGhost,numberOfGhostPoints
             integer side1,side2
             integer n1a,n1b,n2a,n2b,n3a,n3b
             integer nn1a,nn1b,nn2a,nn2b,nn3a,nn3b
             integer extra1a,extra1b,extra2a,extra2b,extra3a,extra3b
             real urv(0:5),usv(0:5),utv(0:5)
             integer iw1,iw2,iw3      
             real a11,a12,a21,a22,det,b0,b1,b2
             real a0,a1,cc0,cc1,d0,d1,dr0,ds0
             real aNormSq,divu,uAve
             real epsRatio,an1,an2,an3,aNorm,aNormi,ua,ub,nDotU,t1,t2,t3
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
       !**   real ulapSq1,vlapSq1,ulapSq2,vlapSq2,wlapSq1,wlapSq2
       !**      integer np1a,np1b,np2a,np2b,np3a,np3b,diff(0:2)
       !**      real rx,ry,rxx,rxy,ryy,rxxx,rxxy,rxyy,ryyy,rxxxx,rxxyy,ryyyy
       !**      real sx,sy,sxx,sxy,syy,sxxx,sxxy,sxyy,syyy,sxxxx,sxxyy,syyyy
       !**      real rv1x(0:2),rv1y(0:2),rv1xx(0:2),rv1xy(0:2),rv1yy(0:2),rv1xxx(0:2),rv1xxy(0:2),rv1xyy(0:2),rv1yyy(0:2),!**           rv1xxxx(0:2),rv1xxyy(0:2),rv1yyyy(0:2)
       !**      real sv1x(0:2),sv1y(0:2),sv1xx(0:2),sv1xy(0:2),sv1yy(0:2),sv1xxx(0:2),sv1xxy(0:2),sv1xyy(0:2),sv1yyy(0:2),!**           sv1xxxx(0:2),sv1xxyy(0:2),sv1yyyy(0:2)
       !**      real rv2x(0:2),rv2y(0:2),rv2xx(0:2),rv2xy(0:2),rv2yy(0:2),rv2xxx(0:2),rv2xxy(0:2),rv2xyy(0:2),rv2yyy(0:2),!**           rv2xxxx(0:2),rv2xxyy(0:2),rv2yyyy(0:2)
       !**      real sv2x(0:2),sv2y(0:2),sv2xx(0:2),sv2xy(0:2),sv2yy(0:2),sv2xxx(0:2),sv2xxy(0:2),sv2xyy(0:2),sv2yyy(0:2),!**           sv2xxxx(0:2),sv2xxyy(0:2),sv2yyyy(0:2)
             integer numberOfEquations,job
             real a2(0:1,0:1),a3(0:2,0:2),a4(0:3,0:3),a8(0:7,0:7),q(0:11),f(0:11),rcond,work(0:11)
             integer ipvt(0:11)
             real resMax,fe(0:3)
             real err
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
             ! integer rectangular,curvilinear
             ! parameter(!   rectangular=0,!   curvilinear=1)
       !     --- start statement function ----
             real bcf
             integer kd,m,n
             real uxOneSided
       !     real rx,ry,rz,sx,sy,sz,tx,ty,tz
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
             if( t<=2*dt .or. debug.gt.3 )then
               write(*,'("# bcOptIsm: t=",1pe12.3)') t
               write(*,'(" bcOptIsm: rho,mu=",2f10.5," gridType=",i2," upwindSOS=",i2," numGhost=",i2)') rho,mu,gridType,upwindSOS,numGhost
               write(*,'(" bcOptIsm: u1c,u2c,u3c,pc=",4i3," twilightZone=",i2," orderOfAccuracy=",i2)') u1c,u2c,u3c,pc,twilightZone,orderOfAccuracy
               write(*,'(" boundaryCondition=",6i4)') ((boundaryCondition(side,axis),side=0,1),axis=0,nd)
               write(*,'(" addBoundaryForcing=",6i4)') ((addBoundaryForcing(side,axis),side=0,1),axis=0,nd)
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
             ! write(*,*) 'bcOffset:',bcOffset
             ! assign corners and edges (3d)
             if( .false. )then ! *wdh* 071027 -- turn this off for now ---
               if( orderOfAccuracy.eq.2 .and. nd.eq.2 )then
                 ! *** fix this for traction BCs with forcing 
                ! For interfaces and TZ it is ok to just set the corners using TZ. This is maybe cheating a bit.
                if( gridType.eq.rectangular )then
                 if( twilightZone.eq.0 )then
                     axis=0
                     axisp1=1
                     i3=gridIndexRange(0,2)
                     numberOfGhostPoints=orderOfAccuracy/2
                     do side1=0,1
                     do side2=0,1
                     if( boundaryCondition(side1,0).eq.tractionBC .and.boundaryCondition(side2,1).eq.tractionBC )then
                       i1=gridIndexRange(side1,0) ! (i1,i2,i3)=corner point
                       i2=gridIndexRange(side2,1)
                       ! write(*,'("bcOpt: assign corner side1,side2,i1,i2,i3=",2i2,3i5)') side1,side2,i1,i2,i3
                       is1=1-2*side1
                       is2=1-2*side2
                   !   dra=dr(0)*is1
                   !   dsa=dr(1)*is2
                       ! First assign normal component of the displacement:
                       ! u.x=u.xxx=0 --> u is even in x
                       ! v.y=v.yyy=0 --> v is even in y
                       do m=1,numberOfGhostPoints
                         js1=is1*m  ! shift to ghost point "m"
                         js2=is2*m
                           u(i1-js1,i2,i3,u1c)=u(i1+js1,i2,i3,u1c)
                           u(i1,i2-js2,i3,u2c)=u(i1,i2+js2,i3,u2c)
                       end do 
                       ! Now assign the tangential components of the displacement 
                       alpha=lambda/(lambda+2.*mu)  
                         js1=is1
                         js2=is2
                           ! u.yy = alpha*u.xx
                           ! v.xx = alpha*v.yy
                           u(i1,i2-js2,i3,u1c)=2.*u(i1,i2,i3,u1c)-u(i1,i2+js2,i3,u1c) +dx(1)**2*alpha*uxx22r(i1,i2,i3,u1c)
                           u(i1-js1,i2,i3,u2c)=2.*u(i1,i2,i3,u2c)-u(i1+js1,i2,i3,u2c) +dx(0)**2*alpha*uyy22r(i1,i2,i3,u2c)
                       ! Now do corner (C) points
                       ! Taylor series: 
                       !   u(-x,-y) = u(x,y) - 2*x*u.x(0,0) - 2*y*u.y(0,0) + O( h^3 )
                     ! ** u(i1-is1,i2-is2,i3,u1c)=u(i1+is1,i2+is2,i3,u1c) -2.*is1*dx(0)*ux22r(i1,i2,i3,u1c) -2.*is2*dx(1)*uy22r(i1,i2,i3,u1c)
                     ! ** u(i1-is1,i2-is2,i3,u2c)=u(i1+is1,i2+is2,i3,u2c) -2.*is1*dx(0)*ux22r(i1,i2,i3,u2c) -2.*is2*dx(1)*uy22r(i1,i2,i3,u2c)
                       ! This version uses u.xy = - v.xx, v.xy = - u.yy
                       u(i1-is1,i2-is2,i3,u1c)=2.*u(i1,i2,i3,u1c) - u(i1+is1,i2+is2,i3,u1c) +dx(0)**2*uxx22r(i1,i2,i3,u1c) - 2.*dx(0)*dx(1)*uxx22r(i1,i2,i3,u2c) +dx(1)**2*uyy22r(i1,i2,i3,u1c)
                       u(i1-is1,i2-is2,i3,u2c)=2.*u(i1,i2,i3,u2c) - u(i1+is1,i2+is2,i3,u2c) +dx(0)**2*uxx22r(i1,i2,i3,u2c) - 2.*dx(0)*dx(1)*uyy22r(i1,i2,i3,u1c) +dx(1)**2*uyy22r(i1,i2,i3,u2c)
                     else if( (boundaryCondition(side1,0).eq.tractionBC .and. boundaryCondition(side2,1).eq.displacementBC) .or.(boundaryCondition(side1,0).eq.displacementBC  .and. boundaryCondition(side2,1).eq.tractionBC) )then 
                       ! displacementBC next to stress free
                       stop 2311
                     else if( boundaryCondition(side1,0).eq.displacementBC .and. boundaryCondition(side2,1).eq.displacementBC )then
                       ! displacementBC next to displacementBC
                       ! do we need to do anything in this case ? *wdh* 071012
                       ! stop 2312
                     else if( boundaryCondition(side1,0).gt.0 .and. boundaryCondition(side2,1).gt.0 )then
                       ! unknown 
                       stop 2313
                     end if
                     end do
                     end do
                 else
                     axis=0
                     axisp1=1
                     i3=gridIndexRange(0,2)
                     numberOfGhostPoints=orderOfAccuracy/2
                     do side1=0,1
                     do side2=0,1
                     if( boundaryCondition(side1,0).eq.tractionBC .and.boundaryCondition(side2,1).eq.tractionBC )then
                       i1=gridIndexRange(side1,0) ! (i1,i2,i3)=corner point
                       i2=gridIndexRange(side2,1)
                       ! write(*,'("bcOpt: assign corner side1,side2,i1,i2,i3=",2i2,3i5)') side1,side2,i1,i2,i3
                       is1=1-2*side1
                       is2=1-2*side2
                   !   dra=dr(0)*is1
                   !   dsa=dr(1)*is2
                       ! First assign normal component of the displacement:
                       ! u.x=u.xxx=0 --> u is even in x
                       ! v.y=v.yyy=0 --> v is even in y
                       do m=1,numberOfGhostPoints
                         js1=is1*m  ! shift to ghost point "m"
                         js2=is2*m
                           u(i1-js1,i2,i3,u1c)=u(i1+js1,i2,i3,u1c)
                           u(i1,i2-js2,i3,u2c)=u(i1,i2+js2,i3,u2c)
                              call ogf2d(ep,xy(i1-js1,i2,i3,0),xy(i1-js1,i2,i3,1),t,um,vm)
                              call ogf2d(ep,xy(i1+js1,i2,i3,0),xy(i1+js1,i2,i3,1),t,up,vp)
                             u(i1-js1,i2,i3,u1c)=u(i1-js1,i2,i3,u1c) + um-up
                              call ogf2d(ep,xy(i1,i2-js2,i3,0),xy(i1,i2-js2,i3,1),t,um,vm)
                              call ogf2d(ep,xy(i1,i2+js2,i3,0),xy(i1,i2+js2,i3,1),t,up,vp)
                             u(i1,i2-js2,i3,u2c)=u(i1,i2-js2,i3,u2c) + vm-vp
                       end do 
                       ! Now assign the tangential components of the displacement 
                       alpha=lambda/(lambda+2.*mu)  
                         js1=is1
                         js2=is2
                           ! u.yy = alpha*u.xx
                           ! v.xx = alpha*v.yy
                           u(i1,i2-js2,i3,u1c)=2.*u(i1,i2,i3,u1c)-u(i1,i2+js2,i3,u1c) +dx(1)**2*alpha*uxx22r(i1,i2,i3,u1c)
                           u(i1-js1,i2,i3,u2c)=2.*u(i1,i2,i3,u2c)-u(i1+js1,i2,i3,u2c) +dx(0)**2*alpha*uyy22r(i1,i2,i3,u2c)
                               call ogDeriv2(ep, 0,2,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, u1c,uxx0, u2c,vxx0)
                               call ogDeriv2(ep, 0,0,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, u1c,uyy0, u2c,vyy0)
                             u(i1,i2-js2,i3,u1c)=u(i1,i2-js2,i3,u1c) +dx(1)**2*( -alpha*uxx0+ uyy0)
                             u(i1-js1,i2,i3,u2c)=u(i1-js1,i2,i3,u2c) +dx(0)**2*( -alpha*vyy0+ vxx0)
                             if( debug.gt.0 )then
                                call ogf2d(ep,xy(i1-js1,i2,i3,0),xy(i1-js1,i2,i3,1),t,um,vm)
                                call ogf2d(ep,xy(i1,i2-js2,i3,0),xy(i1,i2-js2,i3,1),t,up,vp)
                               write(*,'(" bcOpt:corner: i1,i2=",2i4," uerr,verr=",4e10.2)') i1,i2,u(i1-js1,i2,i3,u1c)-um,u(i1-js1,i2,i3,u2c)-vm,u(i1,i2-js2,i3,u1c)-up,u(i1,i2-js2,i3,u2c)-vp
                               ! '
                             end if
                       ! Now do corner (C) points
                       ! Taylor series: 
                       !   u(-x,-y) = u(x,y) - 2*x*u.x(0,0) - 2*y*u.y(0,0) + O( h^3 )
                     ! ** u(i1-is1,i2-is2,i3,u1c)=u(i1+is1,i2+is2,i3,u1c) -2.*is1*dx(0)*ux22r(i1,i2,i3,u1c) -2.*is2*dx(1)*uy22r(i1,i2,i3,u1c)
                     ! ** u(i1-is1,i2-is2,i3,u2c)=u(i1+is1,i2+is2,i3,u2c) -2.*is1*dx(0)*ux22r(i1,i2,i3,u2c) -2.*is2*dx(1)*uy22r(i1,i2,i3,u2c)
                       ! This version uses u.xy = - v.xx, v.xy = - u.yy
                       u(i1-is1,i2-is2,i3,u1c)=2.*u(i1,i2,i3,u1c) - u(i1+is1,i2+is2,i3,u1c) +dx(0)**2*uxx22r(i1,i2,i3,u1c) - 2.*dx(0)*dx(1)*uxx22r(i1,i2,i3,u2c) +dx(1)**2*uyy22r(i1,i2,i3,u1c)
                       u(i1-is1,i2-is2,i3,u2c)=2.*u(i1,i2,i3,u2c) - u(i1+is1,i2+is2,i3,u2c) +dx(0)**2*uxx22r(i1,i2,i3,u2c) - 2.*dx(0)*dx(1)*uyy22r(i1,i2,i3,u1c) +dx(1)**2*uyy22r(i1,i2,i3,u2c)
                     else if( (boundaryCondition(side1,0).eq.tractionBC .and. boundaryCondition(side2,1).eq.displacementBC) .or.(boundaryCondition(side1,0).eq.displacementBC  .and. boundaryCondition(side2,1).eq.tractionBC) )then 
                       ! displacementBC next to stress free
                       stop 2311
                     else if( boundaryCondition(side1,0).eq.displacementBC .and. boundaryCondition(side2,1).eq.displacementBC )then
                       ! displacementBC next to displacementBC
                       ! do we need to do anything in this case ? *wdh* 071012
                       ! stop 2312
                     else if( boundaryCondition(side1,0).gt.0 .and. boundaryCondition(side2,1).gt.0 )then
                       ! unknown 
                       stop 2313
                     end if
                     end do
                     end do
                 end if
                else
                 if( twilightZone.eq.0 )then
                     axis=0
                     axisp1=1
                     i3=gridIndexRange(0,2)
                     numberOfGhostPoints=orderOfAccuracy/2
                     do side1=0,1
                     do side2=0,1
                     if( boundaryCondition(side1,0).eq.tractionBC .and.boundaryCondition(side2,1).eq.tractionBC )then
                       i1=gridIndexRange(side1,0) ! (i1,i2,i3)=corner point
                       i2=gridIndexRange(side2,1)
                       ! write(*,'("bcOpt: assign corner side1,side2,i1,i2,i3=",2i2,3i5)') side1,side2,i1,i2,i3
                       is1=1-2*side1
                       is2=1-2*side2
                   !   dra=dr(0)*is1
                   !   dsa=dr(1)*is2
                       ! First assign normal component of the displacement:
                       ! u.x=u.xxx=0 --> u is even in x
                       ! v.y=v.yyy=0 --> v is even in y
                       do m=1,numberOfGhostPoints
                         js1=is1*m  ! shift to ghost point "m"
                         js2=is2*m
                           stop 1117
                       end do 
                       ! Now assign the tangential components of the displacement 
                       alpha=lambda/(lambda+2.*mu)  
                         js1=is1
                         js2=is2
                           stop 1117
                       ! Now do corner (C) points
                       ! Taylor series: 
                       !   u(-x,-y) = u(x,y) - 2*x*u.x(0,0) - 2*y*u.y(0,0) + O( h^3 )
                     ! ** u(i1-is1,i2-is2,i3,u1c)=u(i1+is1,i2+is2,i3,u1c) -2.*is1*dx(0)*ux22r(i1,i2,i3,u1c) -2.*is2*dx(1)*uy22r(i1,i2,i3,u1c)
                     ! ** u(i1-is1,i2-is2,i3,u2c)=u(i1+is1,i2+is2,i3,u2c) -2.*is1*dx(0)*ux22r(i1,i2,i3,u2c) -2.*is2*dx(1)*uy22r(i1,i2,i3,u2c)
                       ! This version uses u.xy = - v.xx, v.xy = - u.yy
                       u(i1-is1,i2-is2,i3,u1c)=2.*u(i1,i2,i3,u1c) - u(i1+is1,i2+is2,i3,u1c) +dx(0)**2*uxx22r(i1,i2,i3,u1c) - 2.*dx(0)*dx(1)*uxx22r(i1,i2,i3,u2c) +dx(1)**2*uyy22r(i1,i2,i3,u1c)
                       u(i1-is1,i2-is2,i3,u2c)=2.*u(i1,i2,i3,u2c) - u(i1+is1,i2+is2,i3,u2c) +dx(0)**2*uxx22r(i1,i2,i3,u2c) - 2.*dx(0)*dx(1)*uyy22r(i1,i2,i3,u1c) +dx(1)**2*uyy22r(i1,i2,i3,u2c)
                     else if( (boundaryCondition(side1,0).eq.tractionBC .and. boundaryCondition(side2,1).eq.displacementBC) .or.(boundaryCondition(side1,0).eq.displacementBC  .and. boundaryCondition(side2,1).eq.tractionBC) )then 
                       ! displacementBC next to stress free
                       stop 2311
                     else if( boundaryCondition(side1,0).eq.displacementBC .and. boundaryCondition(side2,1).eq.displacementBC )then
                       ! displacementBC next to displacementBC
                       ! do we need to do anything in this case ? *wdh* 071012
                       ! stop 2312
                     else if( boundaryCondition(side1,0).gt.0 .and. boundaryCondition(side2,1).gt.0 )then
                       ! unknown 
                       stop 2313
                     end if
                     end do
                     end do
                 else
                     axis=0
                     axisp1=1
                     i3=gridIndexRange(0,2)
                     numberOfGhostPoints=orderOfAccuracy/2
                     do side1=0,1
                     do side2=0,1
                     if( boundaryCondition(side1,0).eq.tractionBC .and.boundaryCondition(side2,1).eq.tractionBC )then
                       i1=gridIndexRange(side1,0) ! (i1,i2,i3)=corner point
                       i2=gridIndexRange(side2,1)
                       ! write(*,'("bcOpt: assign corner side1,side2,i1,i2,i3=",2i2,3i5)') side1,side2,i1,i2,i3
                       is1=1-2*side1
                       is2=1-2*side2
                   !   dra=dr(0)*is1
                   !   dsa=dr(1)*is2
                       ! First assign normal component of the displacement:
                       ! u.x=u.xxx=0 --> u is even in x
                       ! v.y=v.yyy=0 --> v is even in y
                       do m=1,numberOfGhostPoints
                         js1=is1*m  ! shift to ghost point "m"
                         js2=is2*m
                           stop 1117
                       end do 
                       ! Now assign the tangential components of the displacement 
                       alpha=lambda/(lambda+2.*mu)  
                         js1=is1
                         js2=is2
                           stop 1117
                       ! Now do corner (C) points
                       ! Taylor series: 
                       !   u(-x,-y) = u(x,y) - 2*x*u.x(0,0) - 2*y*u.y(0,0) + O( h^3 )
                     ! ** u(i1-is1,i2-is2,i3,u1c)=u(i1+is1,i2+is2,i3,u1c) -2.*is1*dx(0)*ux22r(i1,i2,i3,u1c) -2.*is2*dx(1)*uy22r(i1,i2,i3,u1c)
                     ! ** u(i1-is1,i2-is2,i3,u2c)=u(i1+is1,i2+is2,i3,u2c) -2.*is1*dx(0)*ux22r(i1,i2,i3,u2c) -2.*is2*dx(1)*uy22r(i1,i2,i3,u2c)
                       ! This version uses u.xy = - v.xx, v.xy = - u.yy
                       u(i1-is1,i2-is2,i3,u1c)=2.*u(i1,i2,i3,u1c) - u(i1+is1,i2+is2,i3,u1c) +dx(0)**2*uxx22r(i1,i2,i3,u1c) - 2.*dx(0)*dx(1)*uxx22r(i1,i2,i3,u2c) +dx(1)**2*uyy22r(i1,i2,i3,u1c)
                       u(i1-is1,i2-is2,i3,u2c)=2.*u(i1,i2,i3,u2c) - u(i1+is1,i2+is2,i3,u2c) +dx(0)**2*uxx22r(i1,i2,i3,u2c) - 2.*dx(0)*dx(1)*uyy22r(i1,i2,i3,u1c) +dx(1)**2*uyy22r(i1,i2,i3,u2c)
                     else if( (boundaryCondition(side1,0).eq.tractionBC .and. boundaryCondition(side2,1).eq.displacementBC) .or.(boundaryCondition(side1,0).eq.displacementBC  .and. boundaryCondition(side2,1).eq.tractionBC) )then 
                       ! displacementBC next to stress free
                       stop 2311
                     else if( boundaryCondition(side1,0).eq.displacementBC .and. boundaryCondition(side2,1).eq.displacementBC )then
                       ! displacementBC next to displacementBC
                       ! do we need to do anything in this case ? *wdh* 071012
                       ! stop 2312
                     else if( boundaryCondition(side1,0).gt.0 .and. boundaryCondition(side2,1).gt.0 )then
                       ! unknown 
                       stop 2313
                     end if
                     end do
                     end do
                 end if
                end if      
               else if( orderOfAccuracy.eq.2 .and. nd.eq.3 )then
                 !$$$       if( gridType.eq.rectangular )then
                 !$$$        if( twilightZone.eq.0 )then
                 !$$$          assignCorners3d(2,rectangular,none)
                 !$$$        else
                 !$$$          assignCorners2d(2,rectangular,twilightZone)
                 !$$$        end if
                 !$$$       else
                 !$$$        if( twilightZone.eq.0 )then
                 !$$$          assignCorners3d(2,curvilinear,none)
                 !$$$        else
                 !$$$          assignCorners3d(2,curvilinear,twilightZone)
                 !$$$        end if
                 !$$$       end if      
               else
                  stop 5533
               end if
             end if
             if( .false. )then
               ! check the boundary forcing arrays: check that bcf(side,axis,i1,i2,i3,m) agrees with bcf00, bcf10, ...
               write(*,*) dim
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
               !*      ! (js1,js2,js3) used to compute tangential derivatives
               !*      js1=0
               !*      js2=0
               !*      js3=0
               !*      if( axisp1.eq.0 )then
               !*        js1=1-2*side
               !*      else if( axisp1.eq.1 )then
               !*        js2=1-2*side
               !*      else if( axisp1.eq.2 )then
               !*        js3=1-2*side
               !*      else
               !*        stop 5
               !*      end if
               !* 
               !*      ! (ks1,ks2,ks3) used to compute second tangential derivative
               !*      ks1=0
               !*      ks2=0
               !*      ks3=0
               !*      if( axisp2.eq.0 )then
               !*        ks1=1-2*side
               !*      else if( axisp2.eq.1 )then
               !*        ks2=1-2*side
               !*      else if( axisp2.eq.2 )then
               !*        ks3=1-2*side
               !*      else
               !*        stop 5
               !*      end if
                    if( debug.gt.7 )then
                      write(*,'(" bcOpt: grid,side,axis=",3i3,", loop bounds: n1a,n1b,n2a,n2b,n3a,n3b=",6i3)') grid,side,axis,n1a,n1b,n2a,n2b,n3a,n3b
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
                write(*,'(" BCF: side,axis=",2i3," bcOffset(side,axis)=",i8)') side,axis,bcOffset(side,axis)
                if( addBoundaryForcing(side,axis).ne.0 )then
                   do i3=n3a,n3b
                   do i2=n2a,n2b
                   do i1=n1a,n1b
                    if( side.eq.0 .and. axis.eq.0 )then
                      do m=0,nd-1
                        tmp = bcf00(i1,i2,i3,m) - bcf(side,axis,i1,i2,i3,m)
                        write(*,'(" BCF(0,0): i=",3i3," f=",e8.2,2x,e8.2," diff=",e8.2)') i1,i2,i3,bcf00(i1,i2,i3,m),bcf(side,axis,i1,i2,i3,m),tmp
                        ! '
                      end do
                    else if( side.eq.1 .and. axis.eq.0 )then
                      do m=0,nd-1
                        tmp = bcf10(i1,i2,i3,m) - bcf(side,axis,i1,i2,i3,m)
                        write(*,'(" BCF(1,0): i=",3i3," f=",e8.2,2x,e8.2," diff=",e8.2)') i1,i2,i3,bcf10(i1,i2,i3,m),bcf(side,axis,i1,i2,i3,m),tmp
                        ! '
                      end do
                    else if( side.eq.0 .and. axis.eq.1 )then
                      do m=0,nd-1
                        tmp = bcf01(i1,i2,i3,m) - bcf(side,axis,i1,i2,i3,m)
                        write(*,'(" BCF(0,1): i=",3i3," f=",e8.2,2x,e8.2," diff=",e8.2)') i1,i2,i3,bcf01(i1,i2,i3,m),bcf(side,axis,i1,i2,i3,m),tmp
                        ! '
                      end do
                    else if( side.eq.1 .and. axis.eq.1 )then
                      do m=0,nd-1
                        tmp = bcf11(i1,i2,i3,m) - bcf(side,axis,i1,i2,i3,m)
                        write(*,'(" BCF(1,1): i=",3i3," f=",e8.2,2x,e8.2," diff=",e8.2)') i1,i2,i3,bcf11(i1,i2,i3,m),bcf(side,axis,i1,i2,i3,m),tmp
                        ! '
                      end do
                    else if( side.eq.0 .and. axis.eq.2 )then
                      do m=0,nd-1
                        tmp = bcf02(i1,i2,i3,m) - bcf(side,axis,i1,i2,i3,m)
                        write(*,'(" BCF(0,2): i=",3i3," f=",e8.2,2x,e8.2," diff=",e8.2)') i1,i2,i3,bcf02(i1,i2,i3,m),bcf(side,axis,i1,i2,i3,m),tmp
                        ! '
                      end do
                    else if( side.eq.1 .and. axis.eq.2 )then
                      do m=0,nd-1
                        tmp = bcf12(i1,i2,i3,m) - bcf(side,axis,i1,i2,i3,m)
                        write(*,'(" BCF(1,2): i=",3i3," f=",e8.2,2x,e8.2," diff=",e8.2)') i1,i2,i3,bcf12(i1,i2,i3,m),bcf(side,axis,i1,i2,i3,m),tmp
                        ! '
                      end do
                    end if
                   end do
                   end do
                   end do
                end if
                end do ! end side
                end do ! end axis
             end if
             ! ---------------------------------------------------------------
             ! first extrap values to ghost points (may be needed at corners)
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
             !*      ! (js1,js2,js3) used to compute tangential derivatives
             !*      js1=0
             !*      js2=0
             !*      js3=0
             !*      if( axisp1.eq.0 )then
             !*        js1=1-2*side
             !*      else if( axisp1.eq.1 )then
             !*        js2=1-2*side
             !*      else if( axisp1.eq.2 )then
             !*        js3=1-2*side
             !*      else
             !*        stop 5
             !*      end if
             !* 
             !*      ! (ks1,ks2,ks3) used to compute second tangential derivative
             !*      ks1=0
             !*      ks2=0
             !*      ks3=0
             !*      if( axisp2.eq.0 )then
             !*        ks1=1-2*side
             !*      else if( axisp2.eq.1 )then
             !*        ks2=1-2*side
             !*      else if( axisp2.eq.2 )then
             !*        ks3=1-2*side
             !*      else
             !*        stop 5
             !*      end if
                  if( debug.gt.7 )then
                    write(*,'(" bcOpt: grid,side,axis=",3i3,", loop bounds: n1a,n1b,n2a,n2b,n3a,n3b=",6i3)') grid,side,axis,n1a,n1b,n2a,n2b,n3a,n3b
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
                 if( orderOfAccuracy==2 )then
                   if( nd==2 )then
                      i3=n3a
                      do i2=nn2a,nn2b
                      do i1=nn1a,nn1b
                       u(i1-is1,i2-is2,i3,u1c)=(3.*u(i1,i2,i3,u1c)-3.*u(i1+is1,i2+is2,i3+is3,u1c)+u(i1+2*is1,i2+2*is2,i3+2*is3,u1c))
                       u(i1-is1,i2-is2,i3,u2c)=(3.*u(i1,i2,i3,u2c)-3.*u(i1+is1,i2+is2,i3+is3,u2c)+u(i1+2*is1,i2+2*is2,i3+2*is3,u2c))
                      end do
                      end do
                   else
                      i3=n3a
                      do i2=nn2a,nn2b
                      do i1=nn1a,nn1b
                       u(i1-is1,i2-is2,i3-is3,u1c)=(3.*u(i1,i2,i3,u1c)-3.*u(i1+is1,i2+is2,i3+is3,u1c)+u(i1+2*is1,i2+2*is2,i3+2*is3,u1c))
                       u(i1-is1,i2-is2,i3-is3,u2c)=(3.*u(i1,i2,i3,u2c)-3.*u(i1+is1,i2+is2,i3+is3,u2c)+u(i1+2*is1,i2+2*is2,i3+2*is3,u2c))
                       u(i1-is1,i2-is2,i3-is3,u3c)=(3.*u(i1,i2,i3,u3c)-3.*u(i1+is1,i2+is2,i3+is3,u3c)+u(i1+2*is1,i2+2*is2,i3+2*is3,u3c))
                      end do
                      end do
                   end if
                 elseif( orderOfAccuracy==4 )then
                   if( nd==2 )then
                      i3=n3a
                      do i2=nn2a,nn2b
                      do i1=nn1a,nn1b
                       u(i1-is1,i2-is2,i3,u1c)=(5.*u(i1,i2,i3,u1c)-10.*u(i1+is1,i2+is2,i3+is3,u1c)+10.*u(i1+2*is1,i2+2*is2,i3+2*is3,u1c)-5.*u(i1+3*is1,i2+3*is2,i3+3*is3,u1c)+u(i1+4*is1,i2+4*is2,i3+4*is3,u1c))
                       u(i1-is1,i2-is2,i3,u2c)=(5.*u(i1,i2,i3,u2c)-10.*u(i1+is1,i2+is2,i3+is3,u2c)+10.*u(i1+2*is1,i2+2*is2,i3+2*is3,u2c)-5.*u(i1+3*is1,i2+3*is2,i3+3*is3,u2c)+u(i1+4*is1,i2+4*is2,i3+4*is3,u2c))
                      end do
                      end do
                   else
                      i3=n3a
                      do i2=nn2a,nn2b
                      do i1=nn1a,nn1b
                       u(i1-is1,i2-is2,i3-is3,u1c)=(5.*u(i1,i2,i3,u1c)-10.*u(i1+is1,i2+is2,i3+is3,u1c)+10.*u(i1+2*is1,i2+2*is2,i3+2*is3,u1c)-5.*u(i1+3*is1,i2+3*is2,i3+3*is3,u1c)+u(i1+4*is1,i2+4*is2,i3+4*is3,u1c))
                       u(i1-is1,i2-is2,i3-is3,u2c)=(5.*u(i1,i2,i3,u2c)-10.*u(i1+is1,i2+is2,i3+is3,u2c)+10.*u(i1+2*is1,i2+2*is2,i3+2*is3,u2c)-5.*u(i1+3*is1,i2+3*is2,i3+3*is3,u2c)+u(i1+4*is1,i2+4*is2,i3+4*is3,u2c))
                       u(i1-is1,i2-is2,i3-is3,u3c)=(5.*u(i1,i2,i3,u3c)-10.*u(i1+is1,i2+is2,i3+is3,u3c)+10.*u(i1+2*is1,i2+2*is2,i3+2*is3,u3c)-5.*u(i1+3*is1,i2+3*is2,i3+3*is3,u3c)+u(i1+4*is1,i2+4*is2,i3+4*is3,u3c))
                      end do
                      end do
                   end if
                 else
                   write(*,*) 'finish me for orderOfAccuracy=',orderOfAccuracy
                   stop 6666
                 end if
               end if
              end do ! end side
              end do ! end axis
             if( nd.eq.2 )then
               ! *********************************** 
               ! **************** 2D ***************
               ! *********************************** 
              if( orderOfAccuracy.eq.2 .and. gridType.eq.rectangular )then
               ! ---------------------------------------------------------
               ! ----------- STAGE 1, 4=2, RECTANGULAR ---------------
               ! ---------------------------------------------------------
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
               !*      ! (js1,js2,js3) used to compute tangential derivatives
               !*      js1=0
               !*      js2=0
               !*      js3=0
               !*      if( axisp1.eq.0 )then
               !*        js1=1-2*side
               !*      else if( axisp1.eq.1 )then
               !*        js2=1-2*side
               !*      else if( axisp1.eq.2 )then
               !*        js3=1-2*side
               !*      else
               !*        stop 5
               !*      end if
               !* 
               !*      ! (ks1,ks2,ks3) used to compute second tangential derivative
               !*      ks1=0
               !*      ks2=0
               !*      ks3=0
               !*      if( axisp2.eq.0 )then
               !*        ks1=1-2*side
               !*      else if( axisp2.eq.1 )then
               !*        ks2=1-2*side
               !*      else if( axisp2.eq.2 )then
               !*        ks3=1-2*side
               !*      else
               !*        stop 5
               !*      end if
                    if( debug.gt.7 )then
                      write(*,'(" bcOpt: grid,side,axis=",3i3,", loop bounds: n1a,n1b,n2a,n2b,n3a,n3b=",6i3)') grid,side,axis,n1a,n1b,n2a,n2b,n3a,n3b
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
                    ! old: displacementBC2dMacro(forcing,2)
                      ! ------ Assigned extended boundaries with dirichlet values ------
                       do i3=nn3a,nn3b
                       do i2=nn2a,nn2b
                       do i1=nn1a,nn1b
                             call ogf2d(ep,xy(i1,i2,i3,0),xy(i1,i2,i3,1),t,u0,v0)
                            u(i1,i2,i3,u1c)=u0
                            u(i1,i2,i3,u2c)=v0
                       end do
                       end do
                       end do
                      ! ----------- Assign ghost points ----------
                      fe(0)=0.; fe(1)=0; ! forcings
                       do i3=n3a,n3b
                       do i2=n2a,n2b
                       do i1=n1a,n1b
                        !   Set div(u) = 0 
                        !      u.x + v.y + w.z = 0
                            if( axis==0 )then
                                ! Twilight-zone: 
                                ! u.x = -v.y + ue.x -ve.y 
                                    call ogDeriv2(ep, 0,1,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, u1c,ux0, u2c,vx0)
                                    call ogDeriv2(ep, 0,0,1,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, u1c,uy0, u2c,vy0)
                                  u(i1-is1,i2-is2,i3,u1c)=u(i1+is1,i2+is2,i3,u1c) - is1*dx(0)*2.*( -uy22r(i1,i2,i3,u2c) + ux0 +vy0 )
                                  u(i1-is1,i2-is2,i3,u2c)=(3.*u(i1,i2,i3,u2c)-3.*u(i1+is1,i2+is2,i3+is3,u2c)+u(i1+2*is1,i2+2*is2,i3+2*is3,u2c))
                            else if( axis==1 )then
                                ! Twilight-zone: 
                                    call ogDeriv2(ep, 0,1,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, u1c,ux0, u2c,vx0)
                                    call ogDeriv2(ep, 0,0,1,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, u1c,uy0, u2c,vy0)
                                  u(i1-is1,i2-is2,i3,u1c)=(3.*u(i1,i2,i3,u1c)-3.*u(i1+is1,i2+is2,i3+is3,u1c)+u(i1+2*is1,i2+2*is2,i3+2*is3,u1c))
                                  u(i1-is1,i2-is2,i3,u2c)=u(i1+is1,i2+is2,i3,u2c) -is2*dx(1)*2.*( -ux22r(i1,i2,i3,u1c) + vy0 + ux0 )
                            else
                              ! axis==2 
                              ! w.z = - ( u.x + v.y )
                              stop 1234
                            end if
                          ! --- Extrapolate any extra ghost --- 
                          do ghost=2,numGhost
                            ! (j1,j2,j3) is the ghost point index
                            j1 = i1 - is1*ghost 
                            j2 = i2 - is2*ghost 
                            j3 = i3 - is3*ghost 
                            u(j1,j2,j3,u1c)=(3.*u(j1+is1,j2+is2,j3+is3,u1c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u1c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u1c))
                            u(j1,j2,j3,u2c)=(3.*u(j1+is1,j2+is2,j3+is3,u2c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u2c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u2c))
                          end do 
                       end do
                       end do
                       end do
                      ! write(*,'("Traction BC -- stop here for now")')
                      ! stop 9999
                  else if( .true. .or. addBoundaryForcing(side,axis)==0 )then  ! FOR NOW SKIP addBoundaryFORCING
                    ! displacementBC2dMacro(noForcing,2)
                      ! ------ Assigned extended boundaries with dirichlet values ------
                       do i3=nn3a,nn3b
                       do i2=nn2a,nn2b
                       do i1=nn1a,nn1b
                          u(i1,i2,i3,u1c)=0.
                          u(i1,i2,i3,u2c)=0.
                       end do
                       end do
                       end do
                      ! ----------- Assign ghost points ----------
                      fe(0)=0.; fe(1)=0; ! forcings
                       do i3=n3a,n3b
                       do i2=n2a,n2b
                       do i1=n1a,n1b
                        !   Set div(u) = 0 
                        !      u.x + v.y + w.z = 0
                            if( axis==0 )then
                                  ! (1) u.x = - v.y
                                  ! (2) Extrap v 
                                  u(i1-is1,i2-is2,i3,u1c)=u(i1+is1,i2+is2,i3,u1c) + is1*dx(0)*2.*uy22r(i1,i2,i3,u2c)
                                  u(i1-is1,i2-is2,i3,u2c)=(3.*u(i1,i2,i3,u2c)-3.*u(i1+is1,i2+is2,i3+is3,u2c)+u(i1+2*is1,i2+2*is2,i3+2*is3,u2c))
                                  ! write(*,'(" DBC: i1,i2=",2i3," div(u)=",1pe12.4)') i1,i2, ux22r(i1,i2,i3,u1c) + uy22r(i1,i2,i3,u2c)
                            else if( axis==1 )then
                                ! (1) extrap u
                                ! (2) v.y = - u.x
                                  u(i1-is1,i2-is2,i3,u1c)=(3.*u(i1,i2,i3,u1c)-3.*u(i1+is1,i2+is2,i3+is3,u1c)+u(i1+2*is1,i2+2*is2,i3+2*is3,u1c))
                                  u(i1-is1,i2-is2,i3,u2c)=u(i1+is1,i2+is2,i3,u2c) + is2*dx(1)*2.*ux22r(i1,i2,i3,u1c)
                            else
                              ! axis==2 
                              ! w.z = - ( u.x + v.y )
                              stop 1234
                            end if
                          ! --- Extrapolate any extra ghost --- 
                          do ghost=2,numGhost
                            ! (j1,j2,j3) is the ghost point index
                            j1 = i1 - is1*ghost 
                            j2 = i2 - is2*ghost 
                            j3 = i3 - is3*ghost 
                            u(j1,j2,j3,u1c)=(3.*u(j1+is1,j2+is2,j3+is3,u1c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u1c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u1c))
                            u(j1,j2,j3,u2c)=(3.*u(j1+is1,j2+is2,j3+is3,u2c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u2c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u2c))
                          end do 
                       end do
                       end do
                       end do
                      ! write(*,'("Traction BC -- stop here for now")')
                      ! stop 9999
                  else 
                      ! ------ Assigned extended boundaries with dirichlet values ------
                       do i3=nn3a,nn3b
                       do i2=nn2a,nn2b
                       do i1=nn1a,nn1b
                          ! Use forcing in array bcf(..)  
                          stop 1234 
                       end do
                       end do
                       end do
                      ! ----------- Assign ghost points ----------
                      fe(0)=0.; fe(1)=0; ! forcings
                       do i3=n3a,n3b
                       do i2=n2a,n2b
                       do i1=n1a,n1b
                        !   Set div(u) = 0 
                        !      u.x + v.y + w.z = 0
                            if( axis==0 )then
                                ! Use forcing in array bcf(..)  
                                stop 1234     
                            else if( axis==1 )then
                                ! Use forcing in array bcf(..)  
                                stop 1234     
                            else
                              ! axis==2 
                              ! w.z = - ( u.x + v.y )
                              stop 1234
                            end if
                          ! --- Extrapolate any extra ghost --- 
                          do ghost=2,numGhost
                            ! (j1,j2,j3) is the ghost point index
                            j1 = i1 - is1*ghost 
                            j2 = i2 - is2*ghost 
                            j3 = i3 - is3*ghost 
                            u(j1,j2,j3,u1c)=(3.*u(j1+is1,j2+is2,j3+is3,u1c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u1c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u1c))
                            u(j1,j2,j3,u2c)=(3.*u(j1+is1,j2+is2,j3+is3,u2c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u2c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u2c))
                          end do 
                       end do
                       end do
                       end do
                      ! write(*,'("Traction BC -- stop here for now")')
                      ! stop 9999
                  end if
                else if( boundaryCondition(side,axis).eq.tractionBC )then 
                  ! Extrapolation now done above
                  ! ! first extrap values to ghost points (may be needed at corners)
                  ! ! *wdh* 081117 -- this is still needed
                  ! write(*,'("bcOptIsm: extrap ghost on traction BC")')
                  ! beginGhostLoops2d()
                  !   u(i1-is1,i2-is2,i3,u1c)=extrap3(u,i1,i2,i3,u1c,is1,is2,is3)
                  !   u(i1-is1,i2-is2,i3,u2c)=extrap3(u,i1,i2,i3,u2c,is1,is2,is3)
                  ! endLoops2d()
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
                      i3=n3a
                      do i2=n2a,n2b
                      do i1=n1a,n1b
                        call ogf2d(ep,xy(i1,i2,i3,0),xy(i1,i2,i3,1),t,u0,v0)
                       u(i1,i2,i3,u1c)=u0
                       u(i1,i2,i3,u2c)=v0
                       do ghost=1,numGhost
                         j1 = i1 - ghost*is1
                         j2 = i2 - ghost*is2
                          call ogf2d(ep,xy(j1,j2,i3,0),xy(j1,j2,i3,1),t,u0,v0)
                         u(j1,j2,i3,u1c)=u0
                         u(j1,j2,i3,u2c)=v0
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
               ! *********** now apply BC's that assign the ghost values *********
               ! ---------------------------------------------------------
               ! ----------- STAGE 2, 4=2, RECTANGULAR ---------------
               ! ---------------------------------------------------------        
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
               !*      ! (js1,js2,js3) used to compute tangential derivatives
               !*      js1=0
               !*      js2=0
               !*      js3=0
               !*      if( axisp1.eq.0 )then
               !*        js1=1-2*side
               !*      else if( axisp1.eq.1 )then
               !*        js2=1-2*side
               !*      else if( axisp1.eq.2 )then
               !*        js3=1-2*side
               !*      else
               !*        stop 5
               !*      end if
               !* 
               !*      ! (ks1,ks2,ks3) used to compute second tangential derivative
               !*      ks1=0
               !*      ks2=0
               !*      ks3=0
               !*      if( axisp2.eq.0 )then
               !*        ks1=1-2*side
               !*      else if( axisp2.eq.1 )then
               !*        ks2=1-2*side
               !*      else if( axisp2.eq.2 )then
               !*        ks3=1-2*side
               !*      else
               !*        stop 5
               !*      end if
                    if( debug.gt.7 )then
                      write(*,'(" bcOpt: grid,side,axis=",3i3,", loop bounds: n1a,n1b,n2a,n2b,n3a,n3b=",6i3)') grid,side,axis,n1a,n1b,n2a,n2b,n3a,n3b
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
                if( boundaryCondition(side,axis).eq.tractionBC )then 
                  ! ------ TRACTION BC ----
                  ! Compressible case was
                  !  (lambda+2*mu) u_x + lambda*v_y = f   left/right
                  !  mu*( u_y + v_x ) = 0                 left/right
                  ! 
                  ! Incompressible:
                  !    -p + 2*mu*u_x = f1  (for pressure eqn)
                  !    u_x + v_y =0 
                  !    mu*( u_y + v_x ) =0  (same)
                  ! alpha=lambda/(lambda+2.*mu)
                  ! beta =1./(lambda+2.*mu)
                  ! write(*,'(" tractionBC: side,axis=",2i3," assignTwilightZone=",i2," numGhost=",i2)') side,axis,assignTwilightZone,numGhost
                  alpha=1. ! old 
                  beta=1.  ! old
                  if( axis.eq.0 )then
                    ! u.x = -v.y  
                    ! v.x = -u.y       
                    if( addBoundaryForcing(side,axis).eq.0 .and. assignTwilightZone==0 )then
                     ! no forcing 
                      i3=n3a
                      do i2=n2a,n2b
                      do i1=n1a,n1b
                      if( mask(i1,i2,i3).gt.0 )then
                       do ghost=1,numGhost
                         js1 = ghost*is1
                         js2 = ghost*is2  
                         ! Note: Use wider stencil on farther ghost  
                         !  u(-2)  = u(2) + (2 * 2*dx) ( ... )         
                         u(i1-js1,i2-js2,i3,u1c)=u(i1+js1,i2+js2,i3,u1c) +js1*dx(0)*2.*uy22r(i1,i2,i3,u2c)
                         u(i1-js1,i2-js2,i3,u2c)=u(i1+js1,i2+js2,i3,u2c) +js1*dx(0)*2.*uy22r(i1,i2,i3,u1c)
                       end do
                      end if
                      end do
                      end do
                    else if( assignTwilightZone.eq.0 )then
                      ! include forcing terms 
                       i3=n3a
                       do i2=n2a,n2b
                       do i1=n1a,n1b
                       if( mask(i1,i2,i3).gt.0 )then
                        do ghost=1,numGhost
                          js1 = ghost*is1
                          js2 = ghost*is2                
                          u(i1-js1,i2-js2,i3,u1c)=u(i1+js1,i2+js2,i3,u1c) +dx(0)*2.*( js1*uy22r(i1,i2,i3,u2c) +bcf(side,axis,i1,i2,i3,u1c) )
                          u(i1-js1,i2-js2,i3,u2c)=u(i1+js1,i2+js2,i3,u2c) +dx(0)*2.*( js1*uy22r(i1,i2,i3,u1c) + (1./mu)*bcf(side,axis,i1,i2,i3,u2c) )
                        end do
                       end if
                       end do
                       end do
                    else
                      ! Twilight-zone: 
                      ! u.x = -v.y + ue.x -ve.y 
                      ! write(*,'(" bcOptIsm: assign traction values for axis==0 and TZ")')
                       i3=n3a
                       do i2=n2a,n2b
                       do i1=n1a,n1b
                       if( mask(i1,i2,i3).gt.0 )then
                          call ogDeriv2(ep, 0,1,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, u1c,ux0, u2c,vx0)
                          call ogDeriv2(ep, 0,0,1,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, u1c,uy0, u2c,vy0)
                        do ghost=1,numGhost
                          js1 = ghost*is1
                          js2 = ghost*is2                
                          u(i1-js1,i2-js2,i3,u1c)=u(i1+js1,i2+js2,i3,u1c) -js1*dx(0)*2.*(-uy22r(i1,i2,i3,u2c) + ux0 + vy0 )
                          u(i1-js1,i2-js2,i3,u2c)=u(i1+js1,i2+js2,i3,u2c) -js1*dx(0)*2.*(-uy22r(i1,i2,i3,u1c) + vx0 + uy0 )
                           call ogf2d(ep,xy(i1-js1,i2-js2,i3,0),xy(i1-js1,i2-js2,i3,1),t,u0,v0)
                          ! write(*,'("tractionBC: i1,i2=",2i3," j1,j2=",2i3," v(ghost)=",1pe12.4," ve=",1pe12.4)') i1,i2,i1-js1,i2-js2,u(i1-js1,i2-js2,i3,u2c),v0 
                          !     write(*,'("i1,i2=",2i3," ux0,vx0,uy0,vy0=",4e10.2)') i1,i2, ux0,vx0,uy0,vy0 
                        end do                
                       end if
                       end do
                       end do
                    end if
                  else ! axis==1 
                    ! u.y = - v.x
                    ! v.y = - u.x 
                    if( addBoundaryForcing(side,axis).eq.0 )then
                      i3=n3a
                      do i2=n2a,n2b
                      do i1=n1a,n1b
                      if( mask(i1,i2,i3).gt.0 )then
                       do ghost=1,numGhost
                         js1 = ghost*is1
                         js2 = ghost*is2                
                         u(i1-js1,i2-js2,i3,u1c)=u(i1+js1,i2+js2,i3,u1c) +js2*dx(1)*2.*ux22r(i1,i2,i3,u2c)
                         u(i1-js1,i2-js2,i3,u2c)=u(i1+js1,i2+js2,i3,u2c) +js2*dx(1)*2.*ux22r(i1,i2,i3,u1c)
                       end do
                      end if
                      end do
                      end do
                    else if( assignTwilightZone.eq.0 )then
                      ! include forcing terms
                       i3=n3a
                       do i2=n2a,n2b
                       do i1=n1a,n1b
                       if( mask(i1,i2,i3).gt.0 )then
                       do ghost=1,numGhost
                         js1 = ghost*is1
                         js2 = ghost*is2                   
                         u(i1-js1,i2-js2,i3,u1c)=u(i1+js1,i2+js2,i3,u1c) +dx(1)*2.*(js2*ux22r(i1,i2,i3,u2c) + (1./mu)*bcf(side,axis,i1,i2,i3,u1c) )
                         u(i1-js1,i2-js2,i3,u2c)=u(i1+js1,i2+js2,i3,u2c) +dx(1)*2.*(js2*ux22r(i1,i2,i3,u1c) + bcf(side,axis,i1,i2,i3,u2c) )
                       end do
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
                        do ghost=1,numGhost
                          js1 = ghost*is1
                          js2 = ghost*is2                     
                          u(i1-js1,i2-js2,i3,u1c)=u(i1+js1,i2+js2,i3,u1c) -js2*dx(1)*2.*(-ux22r(i1,i2,i3,u2c) +uy0+vx0)
                          u(i1-js1,i2-js2,i3,u2c)=u(i1+js1,i2+js2,i3,u2c) -js2*dx(1)*2.*(-ux22r(i1,i2,i3,u1c) +vy0+ux0)
                        end do 
                       end if
                       end do
                       end do
                    end if
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
              else if( orderOfAccuracy.eq.2 .and. gridType.eq.curvilinear )then
               ! ---------------------------------------------------------
               ! ----------- STAGE 1, 4=2, CURVILINEAR ---------------
               ! ---------------------------------------------------------        
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
               !*      ! (js1,js2,js3) used to compute tangential derivatives
               !*      js1=0
               !*      js2=0
               !*      js3=0
               !*      if( axisp1.eq.0 )then
               !*        js1=1-2*side
               !*      else if( axisp1.eq.1 )then
               !*        js2=1-2*side
               !*      else if( axisp1.eq.2 )then
               !*        js3=1-2*side
               !*      else
               !*        stop 5
               !*      end if
               !* 
               !*      ! (ks1,ks2,ks3) used to compute second tangential derivative
               !*      ks1=0
               !*      ks2=0
               !*      ks3=0
               !*      if( axisp2.eq.0 )then
               !*        ks1=1-2*side
               !*      else if( axisp2.eq.1 )then
               !*        ks2=1-2*side
               !*      else if( axisp2.eq.2 )then
               !*        ks3=1-2*side
               !*      else
               !*        stop 5
               !*      end if
                    if( debug.gt.7 )then
                      write(*,'(" bcOpt: grid,side,axis=",3i3,", loop bounds: n1a,n1b,n2a,n2b,n3a,n3b=",6i3)') grid,side,axis,n1a,n1b,n2a,n2b,n3a,n3b
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
                      ! ------ Assigned extended boundaries with dirichlet values ------
                       do i3=nn3a,nn3b
                       do i2=nn2a,nn2b
                       do i1=nn1a,nn1b
                             call ogf2d(ep,xy(i1,i2,i3,0),xy(i1,i2,i3,1),t,u0,v0)
                            u(i1,i2,i3,u1c)=u0
                            u(i1,i2,i3,u2c)=v0
                       end do
                       end do
                       end do
                      ! ----------- Assign ghost points ----------
                      fe(0)=0.; fe(1)=0; ! forcings
                       do i3=n3a,n3b
                       do i2=n2a,n2b
                       do i1=n1a,n1b
                        !   Set div(u) = 0 
                        !      u.x + v.y + w.z = 0
                            ! --- traction ghost : curvilinear  -----
                             an1 = rsxy(i1,i2,i3,axis,0)
                             an2 = rsxy(i1,i2,i3,axis,1)
                             aNormi = -is/max(epsx,sqrt(an1**2 + an2**2))
                             an1=an1*aNormi
                             an2=an2*aNormi
                            ! tangent
                            t1=-an2
                            t2= an1
                            ! here are  the equations we mean to satisfy: 
                            !   u.x + v.y = 0
                            !   Extrapolate t.uv
                            ux = ux22(i1,i2,i3,u1c)
                            vy = uy22(i1,i2,i3,u2c)
                            j1 = i1-is1; j2=i2-is2; j3=i3;
                            f(0) = ux + vy - fe(0) 
                            f(1) =  t1*( u(j1,j2,j3,u1c)-(3.*u(j1+is1,j2+is2,j3+is3,u1c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u1c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u1c)) ) +t2*( u(j1,j2,j3,u2c)-(3.*u(j1+is1,j2+is2,j3+is3,u2c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u2c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u2c)) )
                            !  u.x + v.y =  r.x u.r + s.x u.s  + r.y v.r + s.y v.s
                            !            =  
                            !  [ a2(0,0) a2(0,1) ][ u(-1) ] =  RHS
                            !  [ a2(1,0) a2(1,1) ][ v(-1) ]   
                            a2(0,0)=-is*( rsxy(i1,i2,i3,axis,0) )/(2.*dr(axis))
                            a2(0,1)=-is*( rsxy(i1,i2,i3,axis,1) )/(2.*dr(axis))
                            a2(1,0)= t1 
                            a2(1,1)= t2
                            ! here are the wrong ghostpoint values
                            q(0) = u(i1-is1,i2-is2,i3,u1c)
                            q(1) = u(i1-is1,i2-is2,i3,u2c)
                            ! subtract off the contributions from the wrong values at the ghost points:
                            do n=0,1
                              f(n) = (a2(n,0)*q(0)+a2(n,1)*q(1)) - f(n)
                            end do
                            ! could optimize and do this by hand 
                            call dgeco( a2(0,0), 2, 2, ipvt(0),rcond,work(0))
                            call dgesl( a2(0,0), 2, 2, ipvt(0), f(0), job)
                            u(i1-is1,i2-is2,i3,u1c)=f(0)
                            u(i1-is1,i2-is2,i3,u2c)=f(1)
                            if( debug.gt.3 )then 
                              ! re-evaluate and check residuals
                              ux=ux22(i1,i2,i3,u1c)
                              vy=uy22(i1,i2,i3,u2c)
                              f(0) = ux + vy - fe(0) 
                              f(1) =  t1*( u(j1,j2,j3,u1c)-(3.*u(j1+is1,j2+is2,j3+is3,u1c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u1c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u1c)) ) +t2*( u(j1,j2,j3,u2c)-(3.*u(j1+is1,j2+is2,j3+is3,u2c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u2c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u2c)) )      
                              resMax = max(abs(f(0)),abs(f(1)))      
                              write(*,'(" --> displacement BC: curvilinear fill ghost: residuals =",2(1pe12.4))') f(0),f(1)
                            end if
                          ! --- Extrapolate any extra ghost --- 
                          do ghost=2,numGhost
                            ! (j1,j2,j3) is the ghost point index
                            j1 = i1 - is1*ghost 
                            j2 = i2 - is2*ghost 
                            j3 = i3 - is3*ghost 
                            u(j1,j2,j3,u1c)=(3.*u(j1+is1,j2+is2,j3+is3,u1c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u1c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u1c))
                            u(j1,j2,j3,u2c)=(3.*u(j1+is1,j2+is2,j3+is3,u2c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u2c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u2c))
                          end do 
                       end do
                       end do
                       end do
                      ! write(*,'("Traction BC -- stop here for now")')
                      ! stop 9999
                  else if( .true. .or. addBoundaryForcing(side,axis)==0 )then  ! FOR NOW SKIP addBoundaryFORCING
                      ! ------ Assigned extended boundaries with dirichlet values ------
                       do i3=nn3a,nn3b
                       do i2=nn2a,nn2b
                       do i1=nn1a,nn1b
                          u(i1,i2,i3,u1c)=0.
                          u(i1,i2,i3,u2c)=0.
                       end do
                       end do
                       end do
                      ! ----------- Assign ghost points ----------
                      fe(0)=0.; fe(1)=0; ! forcings
                       do i3=n3a,n3b
                       do i2=n2a,n2b
                       do i1=n1a,n1b
                        !   Set div(u) = 0 
                        !      u.x + v.y + w.z = 0
                            ! --- traction ghost : curvilinear  -----
                             an1 = rsxy(i1,i2,i3,axis,0)
                             an2 = rsxy(i1,i2,i3,axis,1)
                             aNormi = -is/max(epsx,sqrt(an1**2 + an2**2))
                             an1=an1*aNormi
                             an2=an2*aNormi
                            ! tangent
                            t1=-an2
                            t2= an1
                            ! here are  the equations we mean to satisfy: 
                            !   u.x + v.y = 0
                            !   Extrapolate t.uv
                            ux = ux22(i1,i2,i3,u1c)
                            vy = uy22(i1,i2,i3,u2c)
                            j1 = i1-is1; j2=i2-is2; j3=i3;
                            f(0) = ux + vy - fe(0) 
                            f(1) =  t1*( u(j1,j2,j3,u1c)-(3.*u(j1+is1,j2+is2,j3+is3,u1c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u1c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u1c)) ) +t2*( u(j1,j2,j3,u2c)-(3.*u(j1+is1,j2+is2,j3+is3,u2c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u2c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u2c)) )
                            !  u.x + v.y =  r.x u.r + s.x u.s  + r.y v.r + s.y v.s
                            !            =  
                            !  [ a2(0,0) a2(0,1) ][ u(-1) ] =  RHS
                            !  [ a2(1,0) a2(1,1) ][ v(-1) ]   
                            a2(0,0)=-is*( rsxy(i1,i2,i3,axis,0) )/(2.*dr(axis))
                            a2(0,1)=-is*( rsxy(i1,i2,i3,axis,1) )/(2.*dr(axis))
                            a2(1,0)= t1 
                            a2(1,1)= t2
                            ! here are the wrong ghostpoint values
                            q(0) = u(i1-is1,i2-is2,i3,u1c)
                            q(1) = u(i1-is1,i2-is2,i3,u2c)
                            ! subtract off the contributions from the wrong values at the ghost points:
                            do n=0,1
                              f(n) = (a2(n,0)*q(0)+a2(n,1)*q(1)) - f(n)
                            end do
                            ! could optimize and do this by hand 
                            call dgeco( a2(0,0), 2, 2, ipvt(0),rcond,work(0))
                            call dgesl( a2(0,0), 2, 2, ipvt(0), f(0), job)
                            u(i1-is1,i2-is2,i3,u1c)=f(0)
                            u(i1-is1,i2-is2,i3,u2c)=f(1)
                            if( debug.gt.3 )then 
                              ! re-evaluate and check residuals
                              ux=ux22(i1,i2,i3,u1c)
                              vy=uy22(i1,i2,i3,u2c)
                              f(0) = ux + vy - fe(0) 
                              f(1) =  t1*( u(j1,j2,j3,u1c)-(3.*u(j1+is1,j2+is2,j3+is3,u1c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u1c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u1c)) ) +t2*( u(j1,j2,j3,u2c)-(3.*u(j1+is1,j2+is2,j3+is3,u2c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u2c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u2c)) )      
                              resMax = max(abs(f(0)),abs(f(1)))      
                              write(*,'(" --> displacement BC: curvilinear fill ghost: residuals =",2(1pe12.4))') f(0),f(1)
                            end if
                          ! --- Extrapolate any extra ghost --- 
                          do ghost=2,numGhost
                            ! (j1,j2,j3) is the ghost point index
                            j1 = i1 - is1*ghost 
                            j2 = i2 - is2*ghost 
                            j3 = i3 - is3*ghost 
                            u(j1,j2,j3,u1c)=(3.*u(j1+is1,j2+is2,j3+is3,u1c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u1c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u1c))
                            u(j1,j2,j3,u2c)=(3.*u(j1+is1,j2+is2,j3+is3,u2c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u2c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u2c))
                          end do 
                       end do
                       end do
                       end do
                      ! write(*,'("Traction BC -- stop here for now")')
                      ! stop 9999
                  else 
                      ! ------ Assigned extended boundaries with dirichlet values ------
                       do i3=nn3a,nn3b
                       do i2=nn2a,nn2b
                       do i1=nn1a,nn1b
                          ! Use forcing in array bcf(..)  
                          stop 1234 
                       end do
                       end do
                       end do
                      ! ----------- Assign ghost points ----------
                      fe(0)=0.; fe(1)=0; ! forcings
                       do i3=n3a,n3b
                       do i2=n2a,n2b
                       do i1=n1a,n1b
                        !   Set div(u) = 0 
                        !      u.x + v.y + w.z = 0
                            ! --- traction ghost : curvilinear  -----
                             an1 = rsxy(i1,i2,i3,axis,0)
                             an2 = rsxy(i1,i2,i3,axis,1)
                             aNormi = -is/max(epsx,sqrt(an1**2 + an2**2))
                             an1=an1*aNormi
                             an2=an2*aNormi
                            ! tangent
                            t1=-an2
                            t2= an1
                             if( assignTwilightZone.eq.1 )then
                              call ogDeriv(ep,0,1,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,u1c,ux0)
                              call ogDeriv(ep,0,0,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,u2c,vy0)
                              fe(0) = ux0 + vy0
                              fe(1) = 0.
                             end if
                            ! here are  the equations we mean to satisfy: 
                            !   u.x + v.y = 0
                            !   Extrapolate t.uv
                            ux = ux22(i1,i2,i3,u1c)
                            vy = uy22(i1,i2,i3,u2c)
                            j1 = i1-is1; j2=i2-is2; j3=i3;
                            f(0) = ux + vy - fe(0) 
                            f(1) =  t1*( u(j1,j2,j3,u1c)-(3.*u(j1+is1,j2+is2,j3+is3,u1c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u1c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u1c)) ) +t2*( u(j1,j2,j3,u2c)-(3.*u(j1+is1,j2+is2,j3+is3,u2c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u2c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u2c)) )
                            !  u.x + v.y =  r.x u.r + s.x u.s  + r.y v.r + s.y v.s
                            !            =  
                            !  [ a2(0,0) a2(0,1) ][ u(-1) ] =  RHS
                            !  [ a2(1,0) a2(1,1) ][ v(-1) ]   
                            a2(0,0)=-is*( rsxy(i1,i2,i3,axis,0) )/(2.*dr(axis))
                            a2(0,1)=-is*( rsxy(i1,i2,i3,axis,1) )/(2.*dr(axis))
                            a2(1,0)= t1 
                            a2(1,1)= t2
                            ! here are the wrong ghostpoint values
                            q(0) = u(i1-is1,i2-is2,i3,u1c)
                            q(1) = u(i1-is1,i2-is2,i3,u2c)
                            ! subtract off the contributions from the wrong values at the ghost points:
                            do n=0,1
                              f(n) = (a2(n,0)*q(0)+a2(n,1)*q(1)) - f(n)
                            end do
                            ! could optimize and do this by hand 
                            call dgeco( a2(0,0), 2, 2, ipvt(0),rcond,work(0))
                            call dgesl( a2(0,0), 2, 2, ipvt(0), f(0), job)
                            u(i1-is1,i2-is2,i3,u1c)=f(0)
                            u(i1-is1,i2-is2,i3,u2c)=f(1)
                            if( debug.gt.3 )then 
                              ! re-evaluate and check residuals
                              ux=ux22(i1,i2,i3,u1c)
                              vy=uy22(i1,i2,i3,u2c)
                              f(0) = ux + vy - fe(0) 
                              f(1) =  t1*( u(j1,j2,j3,u1c)-(3.*u(j1+is1,j2+is2,j3+is3,u1c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u1c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u1c)) ) +t2*( u(j1,j2,j3,u2c)-(3.*u(j1+is1,j2+is2,j3+is3,u2c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u2c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u2c)) )      
                              resMax = max(abs(f(0)),abs(f(1)))      
                              write(*,'(" --> displacement BC: curvilinear fill ghost: residuals =",2(1pe12.4))') f(0),f(1)
                            end if
                          ! --- Extrapolate any extra ghost --- 
                          do ghost=2,numGhost
                            ! (j1,j2,j3) is the ghost point index
                            j1 = i1 - is1*ghost 
                            j2 = i2 - is2*ghost 
                            j3 = i3 - is3*ghost 
                            u(j1,j2,j3,u1c)=(3.*u(j1+is1,j2+is2,j3+is3,u1c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u1c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u1c))
                            u(j1,j2,j3,u2c)=(3.*u(j1+is1,j2+is2,j3+is3,u2c)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,u2c)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,u2c))
                          end do 
                       end do
                       end do
                       end do
                      ! write(*,'("Traction BC -- stop here for now")')
                      ! stop 9999
                  end if
                  ! ! For now displacement BC's are done by the calling program
                  ! if( addBoundaryForcing(side,axis).eq.0 )then
                  !  displacementBC2dMacro(noForcing,2)
                  ! else
                  !  displacementBC2dMacro(forcing,2)
                  ! end if
                 else if( boundaryCondition(side,axis).eq.tractionBC )then 
                  ! first extrap values to ghost points (may be needed at corners)
                   i3=n3a
                   do i2=nn2a,nn2b
                   do i1=nn1a,nn1b
                    u(i1-is1,i2-is2,i3,u1c)=(3.*u(i1,i2,i3,u1c)-3.*u(i1+is1,i2+is2,i3+is3,u1c)+u(i1+2*is1,i2+2*is2,i3+2*is3,u1c))
                    u(i1-is1,i2-is2,i3,u2c)=(3.*u(i1,i2,i3,u2c)-3.*u(i1+is1,i2+is2,i3+is3,u2c)+u(i1+2*is1,i2+2*is2,i3+2*is3,u2c))
                   end do
                   end do
                 else if( boundaryCondition(side,axis).eq.slipWall )then
                  if( addBoundaryForcing(side,axis).eq.0 )then
                    ! no forcing 
                     i3=n3a
                     do i2=n2a,n2b
                     do i1=n1a,n1b
                     if( mask(i1,i2,i3).gt.0 )then
                     ! (an1,an2) : outward normal is=1-2*side
                      an1 = rsxy(i1,i2,i3,axis,0)
                      an2 = rsxy(i1,i2,i3,axis,1)
                      aNormi = -is/max(epsx,sqrt(an1**2 + an2**2))
                      an1=an1*aNormi
                      an2=an2*aNormi
                     u1 = u(i1,i2,i3,u1c)
                     v1 = u(i1,i2,i3,u2c)
                     nDotU = an1*u1 + an2*v1  
                     u(i1,i2,i3,u1c)=u1 - nDotU*an1
                     u(i1,i2,i3,u2c)=v1 - nDotU*an2
                     end if
                     end do
                     end do
                   else if( assignTwilightZone.eq.0 )then
                     ! include forcing terms 
                     ! n.u = n.g 
                     i3=n3a
                     do i2=n2a,n2b
                     do i1=n1a,n1b
                     if( mask(i1,i2,i3).gt.0 )then
                      an1 = rsxy(i1,i2,i3,axis,0)
                      an2 = rsxy(i1,i2,i3,axis,1)
                      aNormi = -is/max(epsx,sqrt(an1**2 + an2**2))
                      an1=an1*aNormi
                      an2=an2*aNormi
                     u1 = u(i1,i2,i3,u1c)
                     v1 = u(i1,i2,i3,u2c)
                     nDotU = an1*(u1-bcf(side,axis,i1,i2,i3,u1c)) + an2*(v1-bcf(side,axis,i1,i2,i3,u2c))
                     u(i1,i2,i3,u1c)=u1 - nDotU*an1
                     u(i1,i2,i3,u2c)=v1 - nDotU*an2
                     end if
                     end do
                     end do
                   else
                    ! Twilight-zone: 
                    !   n.u = n.ue
                     i3=n3a
                     do i2=n2a,n2b
                     do i1=n1a,n1b
                     if( mask(i1,i2,i3).gt.0 )then
                      an1 = rsxy(i1,i2,i3,axis,0)
                      an2 = rsxy(i1,i2,i3,axis,1)
                      aNormi = -is/max(epsx,sqrt(an1**2 + an2**2))
                      an1=an1*aNormi
                      an2=an2*aNormi
                      call ogf2d(ep,xy(i1,i2,i3,0),xy(i1,i2,i3,1),t,u0,v0)
                     u1 = u(i1,i2,i3,u1c)
                     v1 = u(i1,i2,i3,u2c)
                     nDotU = an1*(u1-u0) + an2*(v1-v0)
                     u(i1,i2,i3,u1c)=u1 - nDotU*an1
                     u(i1,i2,i3,u2c)=v1 - nDotU*an2
                     end if
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
                      i3=n3a
                      do i2=n2a,n2b
                      do i1=n1a,n1b
                        call ogf2d(ep,xy(i1,i2,i3,0),xy(i1,i2,i3,1),t,u0,v0)
                       u(i1,i2,i3,u1c)=u0
                       u(i1,i2,i3,u2c)=v0
                       do ghost=1,numGhost
                         j1 = i1 - ghost*is1
                         j2 = i2 - ghost*is2
                          call ogf2d(ep,xy(j1,j2,i3,0),xy(j1,j2,i3,1),t,u0,v0)
                         u(j1,j2,i3,u1c)=u0
                         u(j1,j2,i3,u2c)=v0
                       end do
                      end do
                      end do
                   end if
                 ! else if( boundaryCondition(side,axis).eq.dirichletBoundaryCondition )then
                 !   ! Set exact values on boundary and ghost
                 !   if( twilightZone.eq.1 )then
                 !     beginLoops2d()
                 !       OGF2D(i1,i2,i3,t,u0,v0)
                 !       u(i1,i2,i3,u1c)=u0
                 !       u(i1,i2,i3,u2c)=v0
                 !       OGF2D(i1-is1,i2-is2,i3,t,u0,v0)
                 !       u(i1-is1,i2-is2,i3,u1c)=u0
                 !       u(i1-is1,i2-is2,i3,u2c)=v0
                 !     endLoops2d()
                 !   end if
                 else if( boundaryCondition(side,axis).eq.symmetry )then
                   write(*,*) 'finish me for symmetry BC'
                  stop 2856
                 else if( boundaryCondition(side,axis).gt.0 )then
                  stop 1193
                 end if
                end do ! end side
                end do ! end axis
               ! ** now apply BC's that assign the ghost values *********
               ! ---------------------------------------------------------
               ! ----------- STAGE 2, 4=2, CURVILINEAR ---------------
               ! ---------------------------------------------------------                
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
               !*      ! (js1,js2,js3) used to compute tangential derivatives
               !*      js1=0
               !*      js2=0
               !*      js3=0
               !*      if( axisp1.eq.0 )then
               !*        js1=1-2*side
               !*      else if( axisp1.eq.1 )then
               !*        js2=1-2*side
               !*      else if( axisp1.eq.2 )then
               !*        js3=1-2*side
               !*      else
               !*        stop 5
               !*      end if
               !* 
               !*      ! (ks1,ks2,ks3) used to compute second tangential derivative
               !*      ks1=0
               !*      ks2=0
               !*      ks3=0
               !*      if( axisp2.eq.0 )then
               !*        ks1=1-2*side
               !*      else if( axisp2.eq.1 )then
               !*        ks2=1-2*side
               !*      else if( axisp2.eq.2 )then
               !*        ks3=1-2*side
               !*      else
               !*        stop 5
               !*      end if
                    if( debug.gt.7 )then
                      write(*,'(" bcOpt: grid,side,axis=",3i3,", loop bounds: n1a,n1b,n2a,n2b,n3a,n3b=",6i3)') grid,side,axis,n1a,n1b,n2a,n2b,n3a,n3b
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
                if( boundaryCondition(side,axis).eq.tractionBC )then 
                  if( addBoundaryForcing(side,axis).eq.0 )then
                     fe(0)=0.; fe(1)=0.;  ! holds forcing 
                      i3=n3a
                      do i2=n2a,n2b
                      do i1=n1a,n1b
                       ! here is the normal (assumed to be the same on both sides)
                        an1 = rsxy(i1,i2,i3,axis,0)
                        an2 = rsxy(i1,i2,i3,axis,1)
                        aNormi = -is/max(epsx,sqrt(an1**2 + an2**2))
                        an1=an1*aNormi
                        an2=an2*aNormi
                       ! --- here are  the equations we mean to satisfy ----
                       !  --> evaluate with the wrong values in the ghost  
                       ux=ux22(i1,i2,i3,u1c)
                       uy=uy22(i1,i2,i3,u1c)
                       vx=ux22(i1,i2,i3,u2c)
                       vy=uy22(i1,i2,i3,u2c) 
                       f(0) = ux + vy                                                - fe(0)
                       f(1) = 2.*an1*an2*( ux - vy ) + ( an2**2 -an1**2 )*( uy + vx ) - fe(1)
                       !  [ a2(0,0) a2(0,1) ][ u(-1) ] =  RHS
                       !  [ a2(1,0) a2(1,1) ][ v(-1) ]   
                       ! u_x + v_y = rx*ur + ry*vy  + sx*us + sy*vs 
                       a2(0,0)=-is*( rsxy(i1,i2,i3,axis,0) )/(2.*dr(axis) )  ! coeff of u(-1) in u_x + v_y = 0 
                       a2(0,1)=-is*( rsxy(i1,i2,i3,axis,1) )/(2.*dr(axis) )  ! coeff of v(-1) in u_x + v_y = 0
                       ! 2*an1*an2*( ux - vy ) + (an2**2 -an1**2 )*( uy + vx )
                       !   ux = rx*ur + sx*us   vx = rx*vr + sx*vs   
                       !   uy = ry*ur + sy*us   vy = ry*vr + sy*vs
                       a2(1,0)=-is*( 2.*an1*an2*(  rsxy(i1,i2,i3,axis,0) ) + ( an2**2 -an1**2 )*( rsxy(i1,i2,i3,axis,1) ) )/(2.*dr(axis))   ! coeff of u(-1)
                       a2(1,1)=-is*( 2.*an1*an2*( -rsxy(i1,i2,i3,axis,1) ) + ( an2**2 -an1**2 )*( rsxy(i1,i2,i3,axis,0) ) )/(2.*dr(axis))   ! coeff of v(-1)
                       ! here are the wrong ghost point values
                       q(0) = u(i1-is1,i2-is2,i3,u1c)
                       q(1) = u(i1-is1,i2-is2,i3,u2c)
                       ! subtract off the contributions from the wrong values at the ghost points:
                       do n=0,1
                         f(n) = (a2(n,0)*q(0)+a2(n,1)*q(1)) - f(n)
                       end do
                       ! *optimize me* solve this by hand 
                       call dgeco( a2(0,0), 2, 2, ipvt(0),rcond,work(0))
                       call dgesl( a2(0,0), 2, 2, ipvt(0), f(0), job)
                       u(i1-is1,i2-is2,i3,u1c)=f(0)
                       u(i1-is1,i2-is2,i3,u2c)=f(1)
                       if( debug.gt.3 )then ! re-evaluate
                         ux=ux22(i1,i2,i3,u1c)
                         uy=uy22(i1,i2,i3,u1c)
                         vx=ux22(i1,i2,i3,u2c)
                         vy=uy22(i1,i2,i3,u2c) 
                         f(0) = ux + vy                                                - fe(0)
                         f(1) = 2*an1*an2*( ux - vy ) + ( an2**2 -an1**2 )*( uy + vx ) - fe(1)  
                         resMax = max(abs(f(0)),abs(f(1)))      
                         write(*,'(" --> traction BC: curvilinear fill ghost: residuals =",2(1pe12.4))') f(0),f(1)
                       end if
                       ! ------ get second ghost by using a wide formula ----
                       !  axis=0: use
                       !     D0r with dr -> 2*dr
                       !  axis=1: use
                       !     D0s with ds -> 2*ds
                       if( numGhost.gt.1 )then
                         ghost=2;                          ! fill this ghost line 
                         js1=ghost*is1; js2=ghost*is2;     ! shift from boundary point (i1,i2,i3) to ghost point
                         if( axis.eq.0 )then
                           iw1=ghost; iw2=1;               ! use wide stencil for i1
                         else
                           iw1=1;     iw2=ghost;           ! use wide stencil for i2
                         end if
                           urv(u1c) = (u(i1+iw1,i2,i3,u1c) - u(i1-iw1,i2,i3,u1c))/(2.*iw1*dr(0))
                           urv(u2c) = (u(i1+iw1,i2,i3,u2c) - u(i1-iw1,i2,i3,u2c))/(2.*iw1*dr(0))
                           usv(u1c) = (u(i1,i2+iw2,i3,u1c) - u(i1,i2-iw2,i3,u1c))/(2.*iw2*dr(1))
                           usv(u2c) = (u(i1,i2+iw2,i3,u2c) - u(i1,i2-iw2,i3,u2c))/(2.*iw2*dr(1))
                           ux=rsxy(i1,i2,i3,0,0)*urv(u1c) + rsxy(i1,i2,i3,1,0)*usv(u1c)
                           vx=rsxy(i1,i2,i3,0,0)*urv(u2c) + rsxy(i1,i2,i3,1,0)*usv(u2c)
                           uy=rsxy(i1,i2,i3,0,1)*urv(u1c) + rsxy(i1,i2,i3,1,1)*usv(u1c)
                           vy=rsxy(i1,i2,i3,0,1)*urv(u2c) + rsxy(i1,i2,i3,1,1)*usv(u2c)
                         ! --- here are  the equations we mean to satisfy ----
                         !  --> evaluate with the wrong values in the ghost 
                         f(0) = ux + vy                                                - fe(0)
                         f(1) = 2*an1*an2*( ux - vy ) + ( an2**2 -an1**2 )*( uy + vx ) - fe(1)
                         !  [ a2(0,0) a2(0,1) ][ u(-1) ] =  RHS
                         !  [ a2(1,0) a2(1,1) ][ v(-1) ]   
                         ! u_x + v_y = rx*ur + ry*vy  + sx*us + sy*vs 
                         a2(0,0)=-is*( rsxy(i1,i2,i3,axis,0) )/(2.*ghost*dr(axis) )  ! coeff of u(-ghost) in u_x + v_y = 0 
                         a2(0,1)=-is*( rsxy(i1,i2,i3,axis,1) )/(2.*ghost*dr(axis) )  ! coeff of v(-ghost) in u_x + v_y = 0
                         ! 2*an1*an2*( ux - vy ) + (an2**2 -an1**2 )*( uy + vx )
                         !   ux = rx*ur + sx*us   vx = rx*vr + sx*vs   
                         !   uy = ry*ur + sy*us   vy = ry*vr + sy*vs
                         a2(1,0)=-is*( 2*an1*an2*(  rsxy(i1,i2,i3,axis,0) ) + ( an2**2 -an1**2 )*( rsxy(i1,i2,i3,axis,1) ) )/(2.*ghost*dr(axis))  ! coeff of u(-ghost)
                         a2(1,1)=-is*( 2*an1*an2*( -rsxy(i1,i2,i3,axis,1) ) + ( an2**2 -an1**2 )*( rsxy(i1,i2,i3,axis,0) ) )/(2.*ghost*dr(axis))  ! coeff of v(-ghost)
                         ! here are the wrong ghost point values
                         q(0) = u(i1-js1,i2-js2,i3,u1c)
                         q(1) = u(i1-js1,i2-js2,i3,u2c)
                         ! subtract off the contributions from the wrong values at the ghost points:
                         do n=0,1
                           f(n) = (a2(n,0)*q(0)+a2(n,1)*q(1)) - f(n)
                         end do
                         ! *optimize me* solve this by hand 
                         call dgeco( a2(0,0), 2, 2, ipvt(0),rcond,work(0))
                         call dgesl( a2(0,0), 2, 2, ipvt(0), f(0), job)
                         u(i1-js1,i2-js2,i3,u1c)=f(0)
                         u(i1-js1,i2-js2,i3,u2c)=f(1)
                         if( debug.gt.3 )then ! re-evaluate
                             urv(u1c) = (u(i1+iw1,i2,i3,u1c) - u(i1-iw1,i2,i3,u1c))/(2.*iw1*dr(0))
                             urv(u2c) = (u(i1+iw1,i2,i3,u2c) - u(i1-iw1,i2,i3,u2c))/(2.*iw1*dr(0))
                             usv(u1c) = (u(i1,i2+iw2,i3,u1c) - u(i1,i2-iw2,i3,u1c))/(2.*iw2*dr(1))
                             usv(u2c) = (u(i1,i2+iw2,i3,u2c) - u(i1,i2-iw2,i3,u2c))/(2.*iw2*dr(1))
                             ux=rsxy(i1,i2,i3,0,0)*urv(u1c) + rsxy(i1,i2,i3,1,0)*usv(u1c)
                             vx=rsxy(i1,i2,i3,0,0)*urv(u2c) + rsxy(i1,i2,i3,1,0)*usv(u2c)
                             uy=rsxy(i1,i2,i3,0,1)*urv(u1c) + rsxy(i1,i2,i3,1,1)*usv(u1c)
                             vy=rsxy(i1,i2,i3,0,1)*urv(u2c) + rsxy(i1,i2,i3,1,1)*usv(u2c)
                           f(0) = ux + vy                                                - fe(0)
                           f(1) = 2*an1*an2*( ux - vy ) + ( an2**2 -an1**2 )*( uy + vx ) - fe(1)  
                           resMax = max(abs(f(0)),abs(f(1)))      
                           write(*,'(" --> traction BC: curvilinear fill ghost=2 using wide approx: residuals =",2(1pe12.4))') f(0),f(1)
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
                       end do  
                      end do
                      end do
                  else
                     fe(0)=0.; fe(1)=0.;  ! holds forcing 
                      i3=n3a
                      do i2=n2a,n2b
                      do i1=n1a,n1b
                       ! here is the normal (assumed to be the same on both sides)
                        an1 = rsxy(i1,i2,i3,axis,0)
                        an2 = rsxy(i1,i2,i3,axis,1)
                        aNormi = -is/max(epsx,sqrt(an1**2 + an2**2))
                        an1=an1*aNormi
                        an2=an2*aNormi
                         if( assignTwilightZone.eq.1 )then
                            call ogDeriv2(ep, 0,1,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, u1c,ux0, u2c,vx0)
                            call ogDeriv2(ep, 0,0,1,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, u1c,uy0, u2c,vy0)
                          fe(0) = ux0 + vy0
                          fe(1) = 2*an1*an2*( ux0 - vy0 ) + ( an2**2 -an1**2 )*( uy0 + vx0 )
                         end if
                       ! --- here are  the equations we mean to satisfy ----
                       !  --> evaluate with the wrong values in the ghost  
                       ux=ux22(i1,i2,i3,u1c)
                       uy=uy22(i1,i2,i3,u1c)
                       vx=ux22(i1,i2,i3,u2c)
                       vy=uy22(i1,i2,i3,u2c) 
                       f(0) = ux + vy                                                - fe(0)
                       f(1) = 2.*an1*an2*( ux - vy ) + ( an2**2 -an1**2 )*( uy + vx ) - fe(1)
                       !  [ a2(0,0) a2(0,1) ][ u(-1) ] =  RHS
                       !  [ a2(1,0) a2(1,1) ][ v(-1) ]   
                       ! u_x + v_y = rx*ur + ry*vy  + sx*us + sy*vs 
                       a2(0,0)=-is*( rsxy(i1,i2,i3,axis,0) )/(2.*dr(axis) )  ! coeff of u(-1) in u_x + v_y = 0 
                       a2(0,1)=-is*( rsxy(i1,i2,i3,axis,1) )/(2.*dr(axis) )  ! coeff of v(-1) in u_x + v_y = 0
                       ! 2*an1*an2*( ux - vy ) + (an2**2 -an1**2 )*( uy + vx )
                       !   ux = rx*ur + sx*us   vx = rx*vr + sx*vs   
                       !   uy = ry*ur + sy*us   vy = ry*vr + sy*vs
                       a2(1,0)=-is*( 2.*an1*an2*(  rsxy(i1,i2,i3,axis,0) ) + ( an2**2 -an1**2 )*( rsxy(i1,i2,i3,axis,1) ) )/(2.*dr(axis))   ! coeff of u(-1)
                       a2(1,1)=-is*( 2.*an1*an2*( -rsxy(i1,i2,i3,axis,1) ) + ( an2**2 -an1**2 )*( rsxy(i1,i2,i3,axis,0) ) )/(2.*dr(axis))   ! coeff of v(-1)
                       ! here are the wrong ghost point values
                       q(0) = u(i1-is1,i2-is2,i3,u1c)
                       q(1) = u(i1-is1,i2-is2,i3,u2c)
                       ! subtract off the contributions from the wrong values at the ghost points:
                       do n=0,1
                         f(n) = (a2(n,0)*q(0)+a2(n,1)*q(1)) - f(n)
                       end do
                       ! *optimize me* solve this by hand 
                       call dgeco( a2(0,0), 2, 2, ipvt(0),rcond,work(0))
                       call dgesl( a2(0,0), 2, 2, ipvt(0), f(0), job)
                       u(i1-is1,i2-is2,i3,u1c)=f(0)
                       u(i1-is1,i2-is2,i3,u2c)=f(1)
                       if( debug.gt.3 )then ! re-evaluate
                         ux=ux22(i1,i2,i3,u1c)
                         uy=uy22(i1,i2,i3,u1c)
                         vx=ux22(i1,i2,i3,u2c)
                         vy=uy22(i1,i2,i3,u2c) 
                         f(0) = ux + vy                                                - fe(0)
                         f(1) = 2*an1*an2*( ux - vy ) + ( an2**2 -an1**2 )*( uy + vx ) - fe(1)  
                         resMax = max(abs(f(0)),abs(f(1)))      
                         write(*,'(" --> traction BC: curvilinear fill ghost: residuals =",2(1pe12.4))') f(0),f(1)
                       end if
                       ! ------ get second ghost by using a wide formula ----
                       !  axis=0: use
                       !     D0r with dr -> 2*dr
                       !  axis=1: use
                       !     D0s with ds -> 2*ds
                       if( numGhost.gt.1 )then
                         ghost=2;                          ! fill this ghost line 
                         js1=ghost*is1; js2=ghost*is2;     ! shift from boundary point (i1,i2,i3) to ghost point
                         if( axis.eq.0 )then
                           iw1=ghost; iw2=1;               ! use wide stencil for i1
                         else
                           iw1=1;     iw2=ghost;           ! use wide stencil for i2
                         end if
                           urv(u1c) = (u(i1+iw1,i2,i3,u1c) - u(i1-iw1,i2,i3,u1c))/(2.*iw1*dr(0))
                           urv(u2c) = (u(i1+iw1,i2,i3,u2c) - u(i1-iw1,i2,i3,u2c))/(2.*iw1*dr(0))
                           usv(u1c) = (u(i1,i2+iw2,i3,u1c) - u(i1,i2-iw2,i3,u1c))/(2.*iw2*dr(1))
                           usv(u2c) = (u(i1,i2+iw2,i3,u2c) - u(i1,i2-iw2,i3,u2c))/(2.*iw2*dr(1))
                           ux=rsxy(i1,i2,i3,0,0)*urv(u1c) + rsxy(i1,i2,i3,1,0)*usv(u1c)
                           vx=rsxy(i1,i2,i3,0,0)*urv(u2c) + rsxy(i1,i2,i3,1,0)*usv(u2c)
                           uy=rsxy(i1,i2,i3,0,1)*urv(u1c) + rsxy(i1,i2,i3,1,1)*usv(u1c)
                           vy=rsxy(i1,i2,i3,0,1)*urv(u2c) + rsxy(i1,i2,i3,1,1)*usv(u2c)
                         ! --- here are  the equations we mean to satisfy ----
                         !  --> evaluate with the wrong values in the ghost 
                         f(0) = ux + vy                                                - fe(0)
                         f(1) = 2*an1*an2*( ux - vy ) + ( an2**2 -an1**2 )*( uy + vx ) - fe(1)
                         !  [ a2(0,0) a2(0,1) ][ u(-1) ] =  RHS
                         !  [ a2(1,0) a2(1,1) ][ v(-1) ]   
                         ! u_x + v_y = rx*ur + ry*vy  + sx*us + sy*vs 
                         a2(0,0)=-is*( rsxy(i1,i2,i3,axis,0) )/(2.*ghost*dr(axis) )  ! coeff of u(-ghost) in u_x + v_y = 0 
                         a2(0,1)=-is*( rsxy(i1,i2,i3,axis,1) )/(2.*ghost*dr(axis) )  ! coeff of v(-ghost) in u_x + v_y = 0
                         ! 2*an1*an2*( ux - vy ) + (an2**2 -an1**2 )*( uy + vx )
                         !   ux = rx*ur + sx*us   vx = rx*vr + sx*vs   
                         !   uy = ry*ur + sy*us   vy = ry*vr + sy*vs
                         a2(1,0)=-is*( 2*an1*an2*(  rsxy(i1,i2,i3,axis,0) ) + ( an2**2 -an1**2 )*( rsxy(i1,i2,i3,axis,1) ) )/(2.*ghost*dr(axis))  ! coeff of u(-ghost)
                         a2(1,1)=-is*( 2*an1*an2*( -rsxy(i1,i2,i3,axis,1) ) + ( an2**2 -an1**2 )*( rsxy(i1,i2,i3,axis,0) ) )/(2.*ghost*dr(axis))  ! coeff of v(-ghost)
                         ! here are the wrong ghost point values
                         q(0) = u(i1-js1,i2-js2,i3,u1c)
                         q(1) = u(i1-js1,i2-js2,i3,u2c)
                         ! subtract off the contributions from the wrong values at the ghost points:
                         do n=0,1
                           f(n) = (a2(n,0)*q(0)+a2(n,1)*q(1)) - f(n)
                         end do
                         ! *optimize me* solve this by hand 
                         call dgeco( a2(0,0), 2, 2, ipvt(0),rcond,work(0))
                         call dgesl( a2(0,0), 2, 2, ipvt(0), f(0), job)
                         u(i1-js1,i2-js2,i3,u1c)=f(0)
                         u(i1-js1,i2-js2,i3,u2c)=f(1)
                         if( debug.gt.3 )then ! re-evaluate
                             urv(u1c) = (u(i1+iw1,i2,i3,u1c) - u(i1-iw1,i2,i3,u1c))/(2.*iw1*dr(0))
                             urv(u2c) = (u(i1+iw1,i2,i3,u2c) - u(i1-iw1,i2,i3,u2c))/(2.*iw1*dr(0))
                             usv(u1c) = (u(i1,i2+iw2,i3,u1c) - u(i1,i2-iw2,i3,u1c))/(2.*iw2*dr(1))
                             usv(u2c) = (u(i1,i2+iw2,i3,u2c) - u(i1,i2-iw2,i3,u2c))/(2.*iw2*dr(1))
                             ux=rsxy(i1,i2,i3,0,0)*urv(u1c) + rsxy(i1,i2,i3,1,0)*usv(u1c)
                             vx=rsxy(i1,i2,i3,0,0)*urv(u2c) + rsxy(i1,i2,i3,1,0)*usv(u2c)
                             uy=rsxy(i1,i2,i3,0,1)*urv(u1c) + rsxy(i1,i2,i3,1,1)*usv(u1c)
                             vy=rsxy(i1,i2,i3,0,1)*urv(u2c) + rsxy(i1,i2,i3,1,1)*usv(u2c)
                           f(0) = ux + vy                                                - fe(0)
                           f(1) = 2*an1*an2*( ux - vy ) + ( an2**2 -an1**2 )*( uy + vx ) - fe(1)  
                           resMax = max(abs(f(0)),abs(f(1)))      
                           write(*,'(" --> traction BC: curvilinear fill ghost=2 using wide approx: residuals =",2(1pe12.4))') f(0),f(1)
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
                       end do  
                      end do
                      end do
                  end if
                else if( boundaryCondition(side,axis).eq.slipWall )then
                  if( addBoundaryForcing(side,axis).eq.0 )then
                   alpha=lambda+2*mu
                    i3=n3a
                    do i2=n2a,n2b
                    do i1=n1a,n1b
                   !  if( debug.gt.2 )then ! re-evaluate
                   ! 
                   !    OGF2D(i1-is1,i2-is2,i3,t,u0,v0)
                   !    write(*,'(" --> slipWall:curv:START i1,i2=",2i4," ghost: (u,v)=",2e10.2," (ue,ve)=",2e10.2)') i1,i2,u(i1-is1,i2-is2,i3,u1c),u(i1-is1,i2-is2,i3,u2c),u0,v0
                   !    ! ' 
                   !  end if
                    ! Solve 
                    !       n.u(-1) = given (by extrapolation)
                    !       n.tauv.t = 0 
                    !    -->   A uv.r + B uv.s = 0
                    ! here is the normal 
                     an1 = rsxy(i1,i2,i3,axis,0)
                     an2 = rsxy(i1,i2,i3,axis,1)
                     aNormi = -is/max(epsx,sqrt(an1**2 + an2**2))
                     an1=an1*aNormi
                     an2=an2*aNormi
                    ! tangent: 
                    t1 =-an2
                    t2 = an1 
                     ux=ux22(i1,i2,i3,u1c)
                     uy=uy22(i1,i2,i3,u1c)
                     vx=ux22(i1,i2,i3,u2c)
                     vy=uy22(i1,i2,i3,u2c)
                     ! components of the stress tensor: 
                     tau11 = alpha*ux + lambda*vy
                     tau12 = mu*( uy + vx )
                     tau21 = tau12
                     tau22 = lambda*ux + alpha*vy
                     ! f(m) holds the current residuals in the equations we mean to satisfy:  
                     ! f(0) = an1*u(i1-is1,i2-is2,i3,u1c)+an2*u(i1-is1,i2-is2,i3,u2c) - ( an1*um + an2*vm )
                     f(0) = 0.
                     f(1) = t1*(an1*tau11+an2*tau21) + t2*(an1*tau12+an2*tau22)
                    !  [ a2(0,0) a2(0,1) ][ u(-1) ] =  RHS
                    !  [ a2(1,0) a2(1,1) ][ v(-1) ]   
                    !a2(0,0)=-is*( an1*alpha *rsxy(i1,i2,i3,axis,0)+an2*mu*rsxy(i1,i2,i3,axis,1) )/(2.*dr(axis))
                    !a2(0,1)=-is*( an1*lambda*rsxy(i1,i2,i3,axis,1)+an2*mu*rsxy(i1,i2,i3,axis,0) )/(2.*dr(axis))
                    !
                    !a2(1,0)=-is*( an1*mu*rsxy(i1,i2,i3,axis,1)+an2*lambda*rsxy(i1,i2,i3,axis,0) )/(2.*dr(axis))
                    !a2(1,1)=-is*( an1*mu*rsxy(i1,i2,i3,axis,0)+an2*alpha *rsxy(i1,i2,i3,axis,1) )/(2.*dr(axis))
                    a2(0,0)=an1
                    a2(0,1)=an2
                    a2(1,0)=-is*( t1*( an1*alpha *rsxy(i1,i2,i3,axis,0)+an2*mu*rsxy(i1,i2,i3,axis,1) )+t2*( an1*mu*rsxy(i1,i2,i3,axis,1)+an2*lambda*rsxy(i1,i2,i3,axis,0) ) )/(2.*dr(axis))
                    a2(1,1)=-is*( t1*( an1*lambda*rsxy(i1,i2,i3,axis,1)+an2*mu*rsxy(i1,i2,i3,axis,0) )+t2*( an1*mu*rsxy(i1,i2,i3,axis,0)+an2*alpha *rsxy(i1,i2,i3,axis,1) ) )/(2.*dr(axis))
                    ! here are the wrong ghostpoint values
                    q(0) = u(i1-is1,i2-is2,i3,u1c)
                    q(1) = u(i1-is1,i2-is2,i3,u2c)
                    ! subtract off the contributions from the wrong values at the ghost points:
                    do n=0,1
                      f(n) = (a2(n,0)*q(0)+a2(n,1)*q(1)) - f(n)
                    end do
                    ! write(*,'(" --> slipWall:curv:Before solve q=",2e10.2," f=",2e10.2)') q(0),q(1),f(0),f(1)
                    ! write(*,'(" --> slipWall:curv:Before solve a2=",4e10.2)') a2(0,0),a2(0,1),a2(1,0),a2(1,1)
                    call dgeco( a2(0,0), 2, 2, ipvt(0),rcond,work(0))
                    call dgesl( a2(0,0), 2, 2, ipvt(0), f(0), job)
                    u(i1-is1,i2-is2,i3,u1c)=f(0)
                    u(i1-is1,i2-is2,i3,u2c)=f(1)
                    if( debug.gt.3 )then ! re-evaluate
                   !    OGF2D(i1-is1,i2-is2,i3,t,u0,v0)
                   !    write(*,'(" --> slipWall:curv: i1,i2=",2i4," (u,v)=",2e10.2," (ue,ve)=",2e10.2)') i1,i2,u(i1-is1,i2-is2,i3,u1c),u(i1-is1,i2-is2,i3,u2c),u0,v0
                   !    u(i1-is1,i2-is2,i3,u1c)=u0
                   !    u(i1-is1,i2-is2,i3,u2c)=v0
                       ux=ux22(i1,i2,i3,u1c)
                       uy=uy22(i1,i2,i3,u1c)
                       vx=ux22(i1,i2,i3,u2c)
                       vy=uy22(i1,i2,i3,u2c)
                       ! components of the stress tensor: 
                       tau11 = alpha*ux + lambda*vy
                       tau12 = mu*( uy + vx )
                       tau21 = tau12
                       tau22 = lambda*ux + alpha*vy
                       ! f(m) holds the current residuals in the equations we mean to satisfy:  
                       ! f(0) = an1*u(i1-is1,i2-is2,i3,u1c)+an2*u(i1-is1,i2-is2,i3,u2c) - ( an1*um + an2*vm )
                       f(0) = 0.
                       f(1) = t1*(an1*tau11+an2*tau21) + t2*(an1*tau12+an2*tau22)
                      write(*,'(" --> slipWall:curv: i1,i2=",2i4," f=",2e10.2," an1,an2=",2f6.2," is,is1,is2=",3i2)') i1,i2,f(0),f(1),an1,an2,is,is1,is2
                      !'
                      ! write(*,'(" --> slipWall:curv: i1,i2,i3=",3i4," dr=",2e10.2)') i1,i2,i3,dr(0),dr(1)
                      ! write(*,'(" --> slipWall:curv: i1,i2,i3=",3i4," assignTwilightZone,addBoundaryForcing=",2i3,", rcond=",e10.2)') i1,i2,i3,assignTwilightZone,addBoundaryForcing(side,axis),rcond
                      ! write(*,'(" --> slipWall:curv: i1,i2=",2i4," rsxy=",4e10.2)') i1,i2,rsxy(i1,i2,i3,0,0),rsxy(i1,i2,i3,1,0),rsxy(i1,i2,i3,0,1),rsxy(i1,i2,i3,1,1)
                        ! '
                    end if
                    end do
                    end do
                  else
                   alpha=lambda+2*mu
                    i3=n3a
                    do i2=n2a,n2b
                    do i1=n1a,n1b
                   !  if( debug.gt.2 )then ! re-evaluate
                   ! 
                   !    OGF2D(i1-is1,i2-is2,i3,t,u0,v0)
                   !    write(*,'(" --> slipWall:curv:START i1,i2=",2i4," ghost: (u,v)=",2e10.2," (ue,ve)=",2e10.2)') i1,i2,u(i1-is1,i2-is2,i3,u1c),u(i1-is1,i2-is2,i3,u2c),u0,v0
                   !    ! ' 
                   !  end if
                    ! Solve 
                    !       n.u(-1) = given (by extrapolation)
                    !       n.tauv.t = 0 
                    !    -->   A uv.r + B uv.s = 0
                    ! here is the normal 
                     an1 = rsxy(i1,i2,i3,axis,0)
                     an2 = rsxy(i1,i2,i3,axis,1)
                     aNormi = -is/max(epsx,sqrt(an1**2 + an2**2))
                     an1=an1*aNormi
                     an2=an2*aNormi
                    ! tangent: 
                    t1 =-an2
                    t2 = an1 
                     ux=ux22(i1,i2,i3,u1c)
                     uy=uy22(i1,i2,i3,u1c)
                     vx=ux22(i1,i2,i3,u2c)
                     vy=uy22(i1,i2,i3,u2c)
                     ! components of the stress tensor: 
                     tau11 = alpha*ux + lambda*vy
                     tau12 = mu*( uy + vx )
                     tau21 = tau12
                     tau22 = lambda*ux + alpha*vy
                     ! f(m) holds the current residuals in the equations we mean to satisfy:  
                     ! f(0) = an1*u(i1-is1,i2-is2,i3,u1c)+an2*u(i1-is1,i2-is2,i3,u2c) - ( an1*um + an2*vm )
                     f(0) = 0.
                     f(1) = t1*(an1*tau11+an2*tau21) + t2*(an1*tau12+an2*tau22)
                      if( assignTwilightZone.eq.0 )then
                       ! forced case: what should this be ??  do nothing for now 
                       !f(0) = f(0) + is*bcf(side,axis,i1,i2,i3,u1c)
                       !f(1) = f(1) + is*bcf(side,axis,i1,i2,i3,u2c)
                      else
                       ! We could do the following: 
                       !OGF2D(i1-is1,i2-is2,i3,t,um,vm)
                       !f(0) = an1*u(i1-is1,i2-is2,i3,u1c)+an2*u(i1-is1,i2-is2,i3,u2c) - ( an1*um + an2*vm )
                         call ogDeriv2(ep, 0,1,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, u1c,ux0, u2c,vx0)
                         call ogDeriv2(ep, 0,0,1,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, u1c,uy0, u2c,vy0)
                       f(1) = f(1) - t1*( an1*(alpha*ux0 + lambda*vy0) + an2*(mu*(uy0+vx0) ) )- t2*( an1*(mu*(uy0+vx0) )          + an2*(lambda*ux0 + alpha*vy0) )
                      end if
                    !  [ a2(0,0) a2(0,1) ][ u(-1) ] =  RHS
                    !  [ a2(1,0) a2(1,1) ][ v(-1) ]   
                    !a2(0,0)=-is*( an1*alpha *rsxy(i1,i2,i3,axis,0)+an2*mu*rsxy(i1,i2,i3,axis,1) )/(2.*dr(axis))
                    !a2(0,1)=-is*( an1*lambda*rsxy(i1,i2,i3,axis,1)+an2*mu*rsxy(i1,i2,i3,axis,0) )/(2.*dr(axis))
                    !
                    !a2(1,0)=-is*( an1*mu*rsxy(i1,i2,i3,axis,1)+an2*lambda*rsxy(i1,i2,i3,axis,0) )/(2.*dr(axis))
                    !a2(1,1)=-is*( an1*mu*rsxy(i1,i2,i3,axis,0)+an2*alpha *rsxy(i1,i2,i3,axis,1) )/(2.*dr(axis))
                    a2(0,0)=an1
                    a2(0,1)=an2
                    a2(1,0)=-is*( t1*( an1*alpha *rsxy(i1,i2,i3,axis,0)+an2*mu*rsxy(i1,i2,i3,axis,1) )+t2*( an1*mu*rsxy(i1,i2,i3,axis,1)+an2*lambda*rsxy(i1,i2,i3,axis,0) ) )/(2.*dr(axis))
                    a2(1,1)=-is*( t1*( an1*lambda*rsxy(i1,i2,i3,axis,1)+an2*mu*rsxy(i1,i2,i3,axis,0) )+t2*( an1*mu*rsxy(i1,i2,i3,axis,0)+an2*alpha *rsxy(i1,i2,i3,axis,1) ) )/(2.*dr(axis))
                    ! here are the wrong ghostpoint values
                    q(0) = u(i1-is1,i2-is2,i3,u1c)
                    q(1) = u(i1-is1,i2-is2,i3,u2c)
                    ! subtract off the contributions from the wrong values at the ghost points:
                    do n=0,1
                      f(n) = (a2(n,0)*q(0)+a2(n,1)*q(1)) - f(n)
                    end do
                    ! write(*,'(" --> slipWall:curv:Before solve q=",2e10.2," f=",2e10.2)') q(0),q(1),f(0),f(1)
                    ! write(*,'(" --> slipWall:curv:Before solve a2=",4e10.2)') a2(0,0),a2(0,1),a2(1,0),a2(1,1)
                    call dgeco( a2(0,0), 2, 2, ipvt(0),rcond,work(0))
                    call dgesl( a2(0,0), 2, 2, ipvt(0), f(0), job)
                    u(i1-is1,i2-is2,i3,u1c)=f(0)
                    u(i1-is1,i2-is2,i3,u2c)=f(1)
                    if( debug.gt.3 )then ! re-evaluate
                   !    OGF2D(i1-is1,i2-is2,i3,t,u0,v0)
                   !    write(*,'(" --> slipWall:curv: i1,i2=",2i4," (u,v)=",2e10.2," (ue,ve)=",2e10.2)') i1,i2,u(i1-is1,i2-is2,i3,u1c),u(i1-is1,i2-is2,i3,u2c),u0,v0
                   !    u(i1-is1,i2-is2,i3,u1c)=u0
                   !    u(i1-is1,i2-is2,i3,u2c)=v0
                       ux=ux22(i1,i2,i3,u1c)
                       uy=uy22(i1,i2,i3,u1c)
                       vx=ux22(i1,i2,i3,u2c)
                       vy=uy22(i1,i2,i3,u2c)
                       ! components of the stress tensor: 
                       tau11 = alpha*ux + lambda*vy
                       tau12 = mu*( uy + vx )
                       tau21 = tau12
                       tau22 = lambda*ux + alpha*vy
                       ! f(m) holds the current residuals in the equations we mean to satisfy:  
                       ! f(0) = an1*u(i1-is1,i2-is2,i3,u1c)+an2*u(i1-is1,i2-is2,i3,u2c) - ( an1*um + an2*vm )
                       f(0) = 0.
                       f(1) = t1*(an1*tau11+an2*tau21) + t2*(an1*tau12+an2*tau22)
                        if( assignTwilightZone.eq.0 )then
                         ! forced case: what should this be ??  do nothing for now 
                         !f(0) = f(0) + is*bcf(side,axis,i1,i2,i3,u1c)
                         !f(1) = f(1) + is*bcf(side,axis,i1,i2,i3,u2c)
                        else
                         ! We could do the following: 
                         !OGF2D(i1-is1,i2-is2,i3,t,um,vm)
                         !f(0) = an1*u(i1-is1,i2-is2,i3,u1c)+an2*u(i1-is1,i2-is2,i3,u2c) - ( an1*um + an2*vm )
                           call ogDeriv2(ep, 0,1,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, u1c,ux0, u2c,vx0)
                           call ogDeriv2(ep, 0,0,1,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, u1c,uy0, u2c,vy0)
                         f(1) = f(1) - t1*( an1*(alpha*ux0 + lambda*vy0) + an2*(mu*(uy0+vx0) ) )- t2*( an1*(mu*(uy0+vx0) )          + an2*(lambda*ux0 + alpha*vy0) )
                        end if
                      write(*,'(" --> slipWall:curv: i1,i2=",2i4," f=",2e10.2," an1,an2=",2f6.2," is,is1,is2=",3i2)') i1,i2,f(0),f(1),an1,an2,is,is1,is2
                      !'
                      ! write(*,'(" --> slipWall:curv: i1,i2,i3=",3i4," dr=",2e10.2)') i1,i2,i3,dr(0),dr(1)
                      ! write(*,'(" --> slipWall:curv: i1,i2,i3=",3i4," assignTwilightZone,addBoundaryForcing=",2i3,", rcond=",e10.2)') i1,i2,i3,assignTwilightZone,addBoundaryForcing(side,axis),rcond
                      ! write(*,'(" --> slipWall:curv: i1,i2=",2i4," rsxy=",4e10.2)') i1,i2,rsxy(i1,i2,i3,0,0),rsxy(i1,i2,i3,1,0),rsxy(i1,i2,i3,0,1),rsxy(i1,i2,i3,1,1)
                        ! '
                    end if
                    end do
                    end do
                  end if
                end if
                end do ! end side
                end do ! end axis
               ! -----------------------------------
               ! -----------4th Order---------------
               ! -----------------------------------
              else if( orderOfAccuracy.eq.4 .and. gridType.eq.rectangular )then
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
               !*      ! (js1,js2,js3) used to compute tangential derivatives
               !*      js1=0
               !*      js2=0
               !*      js3=0
               !*      if( axisp1.eq.0 )then
               !*        js1=1-2*side
               !*      else if( axisp1.eq.1 )then
               !*        js2=1-2*side
               !*      else if( axisp1.eq.2 )then
               !*        js3=1-2*side
               !*      else
               !*        stop 5
               !*      end if
               !* 
               !*      ! (ks1,ks2,ks3) used to compute second tangential derivative
               !*      ks1=0
               !*      ks2=0
               !*      ks3=0
               !*      if( axisp2.eq.0 )then
               !*        ks1=1-2*side
               !*      else if( axisp2.eq.1 )then
               !*        ks2=1-2*side
               !*      else if( axisp2.eq.2 )then
               !*        ks3=1-2*side
               !*      else
               !*        stop 5
               !*      end if
                    if( debug.gt.7 )then
                      write(*,'(" bcOpt: grid,side,axis=",3i3,", loop bounds: n1a,n1b,n2a,n2b,n3a,n3b=",6i3)') grid,side,axis,n1a,n1b,n2a,n2b,n3a,n3b
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
                  write(*,*) 'bcOptIsm: order=4 finish me'
                  stop 4444
                             ! For now displacement BC's are done by the calling program
                  if( addBoundaryForcing(side,axis).eq.0 )then
                    stop 2222
                     do i3=nn3a,nn3b
                     do i2=nn2a,nn2b
                     do i1=nn1a,nn1b
                        u(i1,i2,i3,u1c)=0.
                        u(i1,i2,i3,u2c)=0.
                      ! Ghost points: 
                      !   u.x + v.y = 0
                        u(i1-is1,i2-is2,i3,u1c)=(5.*u(i1,i2,i3,u1c)-10.*u(i1+is1,i2+is2,i3+is3,u1c)+10.*u(i1+2*is1,i2+2*is2,i3+2*is3,u1c)-5.*u(i1+3*is1,i2+3*is2,i3+3*is3,u1c)+u(i1+4*is1,i2+4*is2,i3+4*is3,u1c))
                        u(i1-is1,i2-is2,i3,u2c)=(5.*u(i1,i2,i3,u2c)-10.*u(i1+is1,i2+is2,i3+is3,u2c)+10.*u(i1+2*is1,i2+2*is2,i3+2*is3,u2c)-5.*u(i1+3*is1,i2+3*is2,i3+3*is3,u2c)+u(i1+4*is1,i2+4*is2,i3+4*is3,u2c))
                     end do
                     end do
                     end do
                  else
                    stop 2222
                     do i3=nn3a,nn3b
                     do i2=nn2a,nn2b
                     do i1=nn1a,nn1b
                         call ogf2d(ep,xy(i1,i2,i3,0),xy(i1,i2,i3,1),t,u0,v0)
                        u(i1,i2,i3,u1c)=u0
                        u(i1,i2,i3,u2c)=v0
                      ! Ghost points: 
                      !   u.x + v.y = 0
                        u(i1-is1,i2-is2,i3,u1c)=(5.*u(i1,i2,i3,u1c)-10.*u(i1+is1,i2+is2,i3+is3,u1c)+10.*u(i1+2*is1,i2+2*is2,i3+2*is3,u1c)-5.*u(i1+3*is1,i2+3*is2,i3+3*is3,u1c)+u(i1+4*is1,i2+4*is2,i3+4*is3,u1c))
                        u(i1-is1,i2-is2,i3,u2c)=(5.*u(i1,i2,i3,u2c)-10.*u(i1+is1,i2+is2,i3+is3,u2c)+10.*u(i1+2*is1,i2+2*is2,i3+2*is3,u2c)-5.*u(i1+3*is1,i2+3*is2,i3+3*is3,u2c)+u(i1+4*is1,i2+4*is2,i3+4*is3,u2c))
                     end do
                     end do
                     end do
                  end if
                else if( boundaryCondition(side,axis).eq.tractionBC )then 
                  write(*,*) 'bcOptIsm: order=4 finish me'
                  stop 4444
                  ! first extrap values to ghost points (may be needed at corners)
                   i3=n3a
                   do i2=nn2a,nn2b
                   do i1=nn1a,nn1b
                    u(i1-is1,i2-is2,i3,u1c)=(5.*u(i1,i2,i3,u1c)-10.*u(i1+is1,i2+is2,i3+is3,u1c)+10.*u(i1+2*is1,i2+2*is2,i3+2*is3,u1c)-5.*u(i1+3*is1,i2+3*is2,i3+3*is3,u1c)+u(i1+4*is1,i2+4*is2,i3+4*is3,u1c))
                    u(i1-is1,i2-is2,i3,u2c)=(5.*u(i1,i2,i3,u2c)-10.*u(i1+is1,i2+is2,i3+is3,u2c)+10.*u(i1+2*is1,i2+2*is2,i3+2*is3,u2c)-5.*u(i1+3*is1,i2+3*is2,i3+3*is3,u2c)+u(i1+4*is1,i2+4*is2,i3+4*is3,u2c))
                    u(i1-2*is1,i2-2*is2,i3,u1c)=(5.*u(i1-is1,i2-is2,i3,u1c)-10.*u(i1-is1+is1,i2-is2+is2,i3+is3,u1c)+10.*u(i1-is1+2*is1,i2-is2+2*is2,i3+2*is3,u1c)-5.*u(i1-is1+3*is1,i2-is2+3*is2,i3+3*is3,u1c)+u(i1-is1+4*is1,i2-is2+4*is2,i3+4*is3,u1c))
                    u(i1-2*is1,i2-2*is2,i3,u2c)=(5.*u(i1-is1,i2-is2,i3,u2c)-10.*u(i1-is1+is1,i2-is2+is2,i3+is3,u2c)+10.*u(i1-is1+2*is1,i2-is2+2*is2,i3+2*is3,u2c)-5.*u(i1-is1+3*is1,i2-is2+3*is2,i3+3*is3,u2c)+u(i1-is1+4*is1,i2-is2+4*is2,i3+4*is3,u2c))
                   end do
                   end do
                else if( boundaryCondition(side,axis).eq.dirichletBoundaryCondition )then
                   ! Set exact values on boundary and ghost
                   if( twilightZone.eq.1 )then
                      i3=n3a
                      do i2=n2a,n2b
                      do i1=n1a,n1b
                        call ogf2d(ep,xy(i1,i2,i3,0),xy(i1,i2,i3,1),t,u0,v0)
                       u(i1,i2,i3,u1c)=u0
                       u(i1,i2,i3,u2c)=v0
                       do ghost=1,numGhost
                         j1 = i1 - ghost*is1
                         j2 = i2 - ghost*is2
                          call ogf2d(ep,xy(j1,j2,i3,0),xy(j1,j2,i3,1),t,u0,v0)
                         u(j1,j2,i3,u1c)=u0
                         u(j1,j2,i3,u2c)=v0
                       end do
                      end do
                      end do
                   else
                     stop 777
                   end if
                else if( boundaryCondition(side,axis).eq.symmetry )then
                   ! finish me 
                   stop 7676
                else if( boundaryCondition(side,axis).gt.0 )then
                  stop 1193
                end if
                end do ! end side
                end do ! end axis
               ! ** now apply BC's that assign the ghost values *********
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
               !*      ! (js1,js2,js3) used to compute tangential derivatives
               !*      js1=0
               !*      js2=0
               !*      js3=0
               !*      if( axisp1.eq.0 )then
               !*        js1=1-2*side
               !*      else if( axisp1.eq.1 )then
               !*        js2=1-2*side
               !*      else if( axisp1.eq.2 )then
               !*        js3=1-2*side
               !*      else
               !*        stop 5
               !*      end if
               !* 
               !*      ! (ks1,ks2,ks3) used to compute second tangential derivative
               !*      ks1=0
               !*      ks2=0
               !*      ks3=0
               !*      if( axisp2.eq.0 )then
               !*        ks1=1-2*side
               !*      else if( axisp2.eq.1 )then
               !*        ks2=1-2*side
               !*      else if( axisp2.eq.2 )then
               !*        ks3=1-2*side
               !*      else
               !*        stop 5
               !*      end if
                    if( debug.gt.7 )then
                      write(*,'(" bcOpt: grid,side,axis=",3i3,", loop bounds: n1a,n1b,n2a,n2b,n3a,n3b=",6i3)') grid,side,axis,n1a,n1b,n2a,n2b,n3a,n3b
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
                if( boundaryCondition(side,axis).eq.tractionBC )then 
                  alpha=lambda/(lambda+2.*mu)
                  if( axis.eq.0 )then
                    ! u.x = -alpha*v.y  
                    ! v.x = -u.y       
                   if( addBoundaryForcing(side,axis).eq.0 )then
                     i3=n3a
                     do i2=n2a,n2b
                     do i1=n1a,n1b
                     if( mask(i1,i2,i3).gt.0 )then
                     u(i1-is1,i2-is2,i3,u1c)=uxOneSided(i1,i2,i3,u1c)+is1*dx(0)*4.*alpha*uy42r(i1,i2,i3,u2c)
                     u(i1-is1,i2-is2,i3,u2c)=uxOneSided(i1,i2,i3,u2c)+is1*dx(0)*4.*      uy42r(i1,i2,i3,u1c)
                     end if
                     end do
                     end do
                   else if( assignTwilightZone.eq.0 )then
                     ! finish me
                     stop 1609
                   else
                    ! u.x = -alpha*v.y + ue.x -alpha*ve.y 
                     i3=n3a
                     do i2=n2a,n2b
                     do i1=n1a,n1b
                     if( mask(i1,i2,i3).gt.0 )then
                       call ogDeriv2(ep, 0,1,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, u1c,ux0, u2c,vx0)
                       call ogDeriv2(ep, 0,0,1,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, u1c,uy0, u2c,vy0)
                     u(i1-is1,i2-is2,i3,u1c)=uxOneSided(i1,i2,i3,u1c)-is1*dx(0)*4.*(-alpha*uy22r(i1,i2,i3,u2c) + ux0 +alpha*vy0 )
                     u(i1-is1,i2-is2,i3,u2c)=uxOneSided(i1,i2,i3,u2c)-is1*dx(0)*4.*(-uy22r(i1,i2,i3,u1c) + vx0 + uy0 )
                     !     write(*,'("i1,i2=",2i3," ux0,vx0,uy0,vy0=",4e10.2)') i1,i2, ux0,vx0,uy0,vy0                
                     end if
                     end do
                     end do
                   end if
                  else
                    ! u.y = - v.x
                    ! v.y = -alpha*u.x 
                   if( addBoundaryForcing(side,axis).eq.0 )then
                     i3=n3a
                     do i2=n2a,n2b
                     do i1=n1a,n1b
                     if( mask(i1,i2,i3).gt.0 )then
                     u(i1-is1,i2-is2,i3,u1c)=uxOneSided(i1,i2,i3,u1c)+is2*dx(1)*4.*      ux42r(i1,i2,i3,u2c)
                     u(i1-is1,i2-is2,i3,u2c)=uxOneSided(i1,i2,i3,u2c)+is2*dx(1)*4.*alpha*ux42r(i1,i2,i3,u1c)
                     end if
                     end do
                     end do
                   else if( assignTwilightZone.eq.0 )then
                     ! finish me
                     stop 1633
                   else
                     i3=n3a
                     do i2=n2a,n2b
                     do i1=n1a,n1b
                     if( mask(i1,i2,i3).gt.0 )then
                       call ogDeriv2(ep, 0,1,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, u1c,ux0, u2c,vx0)
                       call ogDeriv2(ep, 0,0,1,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, u1c,uy0, u2c,vy0)
                     u(i1-is1,i2-is2,i3,u1c)=uxOneSided(i1,i2,i3,u1c)-is2*dx(1)*4.*(      -ux42r(i1,i2,i3,u2c)+uy0+vx0)
                     u(i1-is1,i2-is2,i3,u2c)=uxOneSided(i1,i2,i3,u2c)-is2*dx(1)*4.*(-alpha*ux42r(i1,i2,i3,u1c)+vy0+alpha*ux0)
                     end if
                     end do
                     end do
                   end if
                  end if      
                  ! Now extrap second ghost line
                   i3=n3a
                   do i2=nn2a,nn2b
                   do i1=nn1a,nn1b
                    u(i1-2*is1,i2-2*is2,i3,u1c)=(5.*u(i1-is1,i2-is2,i3,u1c)-10.*u(i1-is1+is1,i2-is2+is2,i3+is3,u1c)+10.*u(i1-is1+2*is1,i2-is2+2*is2,i3+2*is3,u1c)-5.*u(i1-is1+3*is1,i2-is2+3*is2,i3+3*is3,u1c)+u(i1-is1+4*is1,i2-is2+4*is2,i3+4*is3,u1c))
                    u(i1-2*is1,i2-2*is2,i3,u2c)=(5.*u(i1-is1,i2-is2,i3,u2c)-10.*u(i1-is1+is1,i2-is2+is2,i3+is3,u2c)+10.*u(i1-is1+2*is1,i2-is2+2*is2,i3+2*is3,u2c)-5.*u(i1-is1+3*is1,i2-is2+3*is2,i3+3*is3,u2c)+u(i1-is1+4*is1,i2-is2+4*is2,i3+4*is3,u2c))
                   end do
                   end do
                end if
                end do ! end side
                end do ! end axis
              else if( orderOfAccuracy.eq.4 .and. gridType.eq.curvilinear )then
               ! *********************************************
               ! ************* 2d Curvilinear ****************
               ! *********************************************
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
               !*      ! (js1,js2,js3) used to compute tangential derivatives
               !*      js1=0
               !*      js2=0
               !*      js3=0
               !*      if( axisp1.eq.0 )then
               !*        js1=1-2*side
               !*      else if( axisp1.eq.1 )then
               !*        js2=1-2*side
               !*      else if( axisp1.eq.2 )then
               !*        js3=1-2*side
               !*      else
               !*        stop 5
               !*      end if
               !* 
               !*      ! (ks1,ks2,ks3) used to compute second tangential derivative
               !*      ks1=0
               !*      ks2=0
               !*      ks3=0
               !*      if( axisp2.eq.0 )then
               !*        ks1=1-2*side
               !*      else if( axisp2.eq.1 )then
               !*        ks2=1-2*side
               !*      else if( axisp2.eq.2 )then
               !*        ks3=1-2*side
               !*      else
               !*        stop 5
               !*      end if
                    if( debug.gt.7 )then
                      write(*,'(" bcOpt: grid,side,axis=",3i3,", loop bounds: n1a,n1b,n2a,n2b,n3a,n3b=",6i3)') grid,side,axis,n1a,n1b,n2a,n2b,n3a,n3b
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
                  ! For now displacement BC's are done by the calling program
                  if( addBoundaryForcing(side,axis).eq.0 )then
                    stop 2222
                     do i3=nn3a,nn3b
                     do i2=nn2a,nn2b
                     do i1=nn1a,nn1b
                        u(i1,i2,i3,u1c)=0.
                        u(i1,i2,i3,u2c)=0.
                      ! Ghost points: 
                      !   u.x + v.y = 0
                        u(i1-is1,i2-is2,i3,u1c)=(5.*u(i1,i2,i3,u1c)-10.*u(i1+is1,i2+is2,i3+is3,u1c)+10.*u(i1+2*is1,i2+2*is2,i3+2*is3,u1c)-5.*u(i1+3*is1,i2+3*is2,i3+3*is3,u1c)+u(i1+4*is1,i2+4*is2,i3+4*is3,u1c))
                        u(i1-is1,i2-is2,i3,u2c)=(5.*u(i1,i2,i3,u2c)-10.*u(i1+is1,i2+is2,i3+is3,u2c)+10.*u(i1+2*is1,i2+2*is2,i3+2*is3,u2c)-5.*u(i1+3*is1,i2+3*is2,i3+3*is3,u2c)+u(i1+4*is1,i2+4*is2,i3+4*is3,u2c))
                     end do
                     end do
                     end do
                  else
                    stop 2222
                     do i3=nn3a,nn3b
                     do i2=nn2a,nn2b
                     do i1=nn1a,nn1b
                         call ogf2d(ep,xy(i1,i2,i3,0),xy(i1,i2,i3,1),t,u0,v0)
                        u(i1,i2,i3,u1c)=u0
                        u(i1,i2,i3,u2c)=v0
                      ! Ghost points: 
                      !   u.x + v.y = 0
                        u(i1-is1,i2-is2,i3,u1c)=(5.*u(i1,i2,i3,u1c)-10.*u(i1+is1,i2+is2,i3+is3,u1c)+10.*u(i1+2*is1,i2+2*is2,i3+2*is3,u1c)-5.*u(i1+3*is1,i2+3*is2,i3+3*is3,u1c)+u(i1+4*is1,i2+4*is2,i3+4*is3,u1c))
                        u(i1-is1,i2-is2,i3,u2c)=(5.*u(i1,i2,i3,u2c)-10.*u(i1+is1,i2+is2,i3+is3,u2c)+10.*u(i1+2*is1,i2+2*is2,i3+2*is3,u2c)-5.*u(i1+3*is1,i2+3*is2,i3+3*is3,u2c)+u(i1+4*is1,i2+4*is2,i3+4*is3,u2c))
                     end do
                     end do
                     end do
                  end if
                 else if( boundaryCondition(side,axis).eq.tractionBC )then 
                  ! first extrap values to ghost points (may be needed at corners)
                   i3=n3a
                   do i2=nn2a,nn2b
                   do i1=nn1a,nn1b
                    u(i1-is1,i2-is2,i3,u1c)=(5.*u(i1,i2,i3,u1c)-10.*u(i1+is1,i2+is2,i3+is3,u1c)+10.*u(i1+2*is1,i2+2*is2,i3+2*is3,u1c)-5.*u(i1+3*is1,i2+3*is2,i3+3*is3,u1c)+u(i1+4*is1,i2+4*is2,i3+4*is3,u1c))
                    u(i1-is1,i2-is2,i3,u2c)=(5.*u(i1,i2,i3,u2c)-10.*u(i1+is1,i2+is2,i3+is3,u2c)+10.*u(i1+2*is1,i2+2*is2,i3+2*is3,u2c)-5.*u(i1+3*is1,i2+3*is2,i3+3*is3,u2c)+u(i1+4*is1,i2+4*is2,i3+4*is3,u2c))
                   end do
                   end do
                 else if( boundaryCondition(side,axis).eq.dirichletBoundaryCondition .or.boundaryCondition(side,axis).eq.symmetry )then
                  write(*,*) 'finish me BC'
                  stop 2222
                 else if( boundaryCondition(side,axis).gt.0 )then
                  stop 1193
                 end if
                end do ! end side
                end do ! end axis
               ! ** now apply BC's that assign the ghost values *********
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
               !*      ! (js1,js2,js3) used to compute tangential derivatives
               !*      js1=0
               !*      js2=0
               !*      js3=0
               !*      if( axisp1.eq.0 )then
               !*        js1=1-2*side
               !*      else if( axisp1.eq.1 )then
               !*        js2=1-2*side
               !*      else if( axisp1.eq.2 )then
               !*        js3=1-2*side
               !*      else
               !*        stop 5
               !*      end if
               !* 
               !*      ! (ks1,ks2,ks3) used to compute second tangential derivative
               !*      ks1=0
               !*      ks2=0
               !*      ks3=0
               !*      if( axisp2.eq.0 )then
               !*        ks1=1-2*side
               !*      else if( axisp2.eq.1 )then
               !*        ks2=1-2*side
               !*      else if( axisp2.eq.2 )then
               !*        ks3=1-2*side
               !*      else
               !*        stop 5
               !*      end if
                    if( debug.gt.7 )then
                      write(*,'(" bcOpt: grid,side,axis=",3i3,", loop bounds: n1a,n1b,n2a,n2b,n3a,n3b=",6i3)') grid,side,axis,n1a,n1b,n2a,n2b,n3a,n3b
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
                if( boundaryCondition(side,axis).eq.tractionBC )then 
                  if( addBoundaryForcing(side,axis).eq.0 )then
                     fe(0)=0.; fe(1)=0.;  ! holds forcing 
                      i3=n3a
                      do i2=n2a,n2b
                      do i1=n1a,n1b
                       ! here is the normal (assumed to be the same on both sides)
                        an1 = rsxy(i1,i2,i3,axis,0)
                        an2 = rsxy(i1,i2,i3,axis,1)
                        aNormi = -is/max(epsx,sqrt(an1**2 + an2**2))
                        an1=an1*aNormi
                        an2=an2*aNormi
                       ! --- here are  the equations we mean to satisfy ----
                       !  --> evaluate with the wrong values in the ghost  
                       ux=ux22(i1,i2,i3,u1c)
                       uy=uy22(i1,i2,i3,u1c)
                       vx=ux22(i1,i2,i3,u2c)
                       vy=uy22(i1,i2,i3,u2c) 
                       f(0) = ux + vy                                                - fe(0)
                       f(1) = 2.*an1*an2*( ux - vy ) + ( an2**2 -an1**2 )*( uy + vx ) - fe(1)
                       !  [ a2(0,0) a2(0,1) ][ u(-1) ] =  RHS
                       !  [ a2(1,0) a2(1,1) ][ v(-1) ]   
                       ! u_x + v_y = rx*ur + ry*vy  + sx*us + sy*vs 
                       a2(0,0)=-is*( rsxy(i1,i2,i3,axis,0) )/(2.*dr(axis) )  ! coeff of u(-1) in u_x + v_y = 0 
                       a2(0,1)=-is*( rsxy(i1,i2,i3,axis,1) )/(2.*dr(axis) )  ! coeff of v(-1) in u_x + v_y = 0
                       ! 2*an1*an2*( ux - vy ) + (an2**2 -an1**2 )*( uy + vx )
                       !   ux = rx*ur + sx*us   vx = rx*vr + sx*vs   
                       !   uy = ry*ur + sy*us   vy = ry*vr + sy*vs
                       a2(1,0)=-is*( 2.*an1*an2*(  rsxy(i1,i2,i3,axis,0) ) + ( an2**2 -an1**2 )*( rsxy(i1,i2,i3,axis,1) ) )/(2.*dr(axis))   ! coeff of u(-1)
                       a2(1,1)=-is*( 2.*an1*an2*( -rsxy(i1,i2,i3,axis,1) ) + ( an2**2 -an1**2 )*( rsxy(i1,i2,i3,axis,0) ) )/(2.*dr(axis))   ! coeff of v(-1)
                       ! here are the wrong ghost point values
                       q(0) = u(i1-is1,i2-is2,i3,u1c)
                       q(1) = u(i1-is1,i2-is2,i3,u2c)
                       ! subtract off the contributions from the wrong values at the ghost points:
                       do n=0,1
                         f(n) = (a2(n,0)*q(0)+a2(n,1)*q(1)) - f(n)
                       end do
                       ! *optimize me* solve this by hand 
                       call dgeco( a2(0,0), 2, 2, ipvt(0),rcond,work(0))
                       call dgesl( a2(0,0), 2, 2, ipvt(0), f(0), job)
                       u(i1-is1,i2-is2,i3,u1c)=f(0)
                       u(i1-is1,i2-is2,i3,u2c)=f(1)
                       if( debug.gt.3 )then ! re-evaluate
                         ux=ux22(i1,i2,i3,u1c)
                         uy=uy22(i1,i2,i3,u1c)
                         vx=ux22(i1,i2,i3,u2c)
                         vy=uy22(i1,i2,i3,u2c) 
                         f(0) = ux + vy                                                - fe(0)
                         f(1) = 2*an1*an2*( ux - vy ) + ( an2**2 -an1**2 )*( uy + vx ) - fe(1)  
                         resMax = max(abs(f(0)),abs(f(1)))      
                         write(*,'(" --> traction BC: curvilinear fill ghost: residuals =",2(1pe12.4))') f(0),f(1)
                       end if
                       ! ------ get second ghost by using a wide formula ----
                       !  axis=0: use
                       !     D0r with dr -> 2*dr
                       !  axis=1: use
                       !     D0s with ds -> 2*ds
                       if( numGhost.gt.1 )then
                         ghost=2;                          ! fill this ghost line 
                         js1=ghost*is1; js2=ghost*is2;     ! shift from boundary point (i1,i2,i3) to ghost point
                         if( axis.eq.0 )then
                           iw1=ghost; iw2=1;               ! use wide stencil for i1
                         else
                           iw1=1;     iw2=ghost;           ! use wide stencil for i2
                         end if
                           urv(u1c) = (u(i1+iw1,i2,i3,u1c) - u(i1-iw1,i2,i3,u1c))/(2.*iw1*dr(0))
                           urv(u2c) = (u(i1+iw1,i2,i3,u2c) - u(i1-iw1,i2,i3,u2c))/(2.*iw1*dr(0))
                           usv(u1c) = (u(i1,i2+iw2,i3,u1c) - u(i1,i2-iw2,i3,u1c))/(2.*iw2*dr(1))
                           usv(u2c) = (u(i1,i2+iw2,i3,u2c) - u(i1,i2-iw2,i3,u2c))/(2.*iw2*dr(1))
                           ux=rsxy(i1,i2,i3,0,0)*urv(u1c) + rsxy(i1,i2,i3,1,0)*usv(u1c)
                           vx=rsxy(i1,i2,i3,0,0)*urv(u2c) + rsxy(i1,i2,i3,1,0)*usv(u2c)
                           uy=rsxy(i1,i2,i3,0,1)*urv(u1c) + rsxy(i1,i2,i3,1,1)*usv(u1c)
                           vy=rsxy(i1,i2,i3,0,1)*urv(u2c) + rsxy(i1,i2,i3,1,1)*usv(u2c)
                         ! --- here are  the equations we mean to satisfy ----
                         !  --> evaluate with the wrong values in the ghost 
                         f(0) = ux + vy                                                - fe(0)
                         f(1) = 2*an1*an2*( ux - vy ) + ( an2**2 -an1**2 )*( uy + vx ) - fe(1)
                         !  [ a2(0,0) a2(0,1) ][ u(-1) ] =  RHS
                         !  [ a2(1,0) a2(1,1) ][ v(-1) ]   
                         ! u_x + v_y = rx*ur + ry*vy  + sx*us + sy*vs 
                         a2(0,0)=-is*( rsxy(i1,i2,i3,axis,0) )/(2.*ghost*dr(axis) )  ! coeff of u(-ghost) in u_x + v_y = 0 
                         a2(0,1)=-is*( rsxy(i1,i2,i3,axis,1) )/(2.*ghost*dr(axis) )  ! coeff of v(-ghost) in u_x + v_y = 0
                         ! 2*an1*an2*( ux - vy ) + (an2**2 -an1**2 )*( uy + vx )
                         !   ux = rx*ur + sx*us   vx = rx*vr + sx*vs   
                         !   uy = ry*ur + sy*us   vy = ry*vr + sy*vs
                         a2(1,0)=-is*( 2*an1*an2*(  rsxy(i1,i2,i3,axis,0) ) + ( an2**2 -an1**2 )*( rsxy(i1,i2,i3,axis,1) ) )/(2.*ghost*dr(axis))  ! coeff of u(-ghost)
                         a2(1,1)=-is*( 2*an1*an2*( -rsxy(i1,i2,i3,axis,1) ) + ( an2**2 -an1**2 )*( rsxy(i1,i2,i3,axis,0) ) )/(2.*ghost*dr(axis))  ! coeff of v(-ghost)
                         ! here are the wrong ghost point values
                         q(0) = u(i1-js1,i2-js2,i3,u1c)
                         q(1) = u(i1-js1,i2-js2,i3,u2c)
                         ! subtract off the contributions from the wrong values at the ghost points:
                         do n=0,1
                           f(n) = (a2(n,0)*q(0)+a2(n,1)*q(1)) - f(n)
                         end do
                         ! *optimize me* solve this by hand 
                         call dgeco( a2(0,0), 2, 2, ipvt(0),rcond,work(0))
                         call dgesl( a2(0,0), 2, 2, ipvt(0), f(0), job)
                         u(i1-js1,i2-js2,i3,u1c)=f(0)
                         u(i1-js1,i2-js2,i3,u2c)=f(1)
                         if( debug.gt.3 )then ! re-evaluate
                             urv(u1c) = (u(i1+iw1,i2,i3,u1c) - u(i1-iw1,i2,i3,u1c))/(2.*iw1*dr(0))
                             urv(u2c) = (u(i1+iw1,i2,i3,u2c) - u(i1-iw1,i2,i3,u2c))/(2.*iw1*dr(0))
                             usv(u1c) = (u(i1,i2+iw2,i3,u1c) - u(i1,i2-iw2,i3,u1c))/(2.*iw2*dr(1))
                             usv(u2c) = (u(i1,i2+iw2,i3,u2c) - u(i1,i2-iw2,i3,u2c))/(2.*iw2*dr(1))
                             ux=rsxy(i1,i2,i3,0,0)*urv(u1c) + rsxy(i1,i2,i3,1,0)*usv(u1c)
                             vx=rsxy(i1,i2,i3,0,0)*urv(u2c) + rsxy(i1,i2,i3,1,0)*usv(u2c)
                             uy=rsxy(i1,i2,i3,0,1)*urv(u1c) + rsxy(i1,i2,i3,1,1)*usv(u1c)
                             vy=rsxy(i1,i2,i3,0,1)*urv(u2c) + rsxy(i1,i2,i3,1,1)*usv(u2c)
                           f(0) = ux + vy                                                - fe(0)
                           f(1) = 2*an1*an2*( ux - vy ) + ( an2**2 -an1**2 )*( uy + vx ) - fe(1)  
                           resMax = max(abs(f(0)),abs(f(1)))      
                           write(*,'(" --> traction BC: curvilinear fill ghost=2 using wide approx: residuals =",2(1pe12.4))') f(0),f(1)
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
                       end do  
                      end do
                      end do
                  else
                     fe(0)=0.; fe(1)=0.;  ! holds forcing 
                      i3=n3a
                      do i2=n2a,n2b
                      do i1=n1a,n1b
                       ! here is the normal (assumed to be the same on both sides)
                        an1 = rsxy(i1,i2,i3,axis,0)
                        an2 = rsxy(i1,i2,i3,axis,1)
                        aNormi = -is/max(epsx,sqrt(an1**2 + an2**2))
                        an1=an1*aNormi
                        an2=an2*aNormi
                         if( assignTwilightZone.eq.1 )then
                            call ogDeriv2(ep, 0,1,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, u1c,ux0, u2c,vx0)
                            call ogDeriv2(ep, 0,0,1,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, u1c,uy0, u2c,vy0)
                          fe(0) = ux0 + vy0
                          fe(1) = 2*an1*an2*( ux0 - vy0 ) + ( an2**2 -an1**2 )*( uy0 + vx0 )
                         end if
                       ! --- here are  the equations we mean to satisfy ----
                       !  --> evaluate with the wrong values in the ghost  
                       ux=ux22(i1,i2,i3,u1c)
                       uy=uy22(i1,i2,i3,u1c)
                       vx=ux22(i1,i2,i3,u2c)
                       vy=uy22(i1,i2,i3,u2c) 
                       f(0) = ux + vy                                                - fe(0)
                       f(1) = 2.*an1*an2*( ux - vy ) + ( an2**2 -an1**2 )*( uy + vx ) - fe(1)
                       !  [ a2(0,0) a2(0,1) ][ u(-1) ] =  RHS
                       !  [ a2(1,0) a2(1,1) ][ v(-1) ]   
                       ! u_x + v_y = rx*ur + ry*vy  + sx*us + sy*vs 
                       a2(0,0)=-is*( rsxy(i1,i2,i3,axis,0) )/(2.*dr(axis) )  ! coeff of u(-1) in u_x + v_y = 0 
                       a2(0,1)=-is*( rsxy(i1,i2,i3,axis,1) )/(2.*dr(axis) )  ! coeff of v(-1) in u_x + v_y = 0
                       ! 2*an1*an2*( ux - vy ) + (an2**2 -an1**2 )*( uy + vx )
                       !   ux = rx*ur + sx*us   vx = rx*vr + sx*vs   
                       !   uy = ry*ur + sy*us   vy = ry*vr + sy*vs
                       a2(1,0)=-is*( 2.*an1*an2*(  rsxy(i1,i2,i3,axis,0) ) + ( an2**2 -an1**2 )*( rsxy(i1,i2,i3,axis,1) ) )/(2.*dr(axis))   ! coeff of u(-1)
                       a2(1,1)=-is*( 2.*an1*an2*( -rsxy(i1,i2,i3,axis,1) ) + ( an2**2 -an1**2 )*( rsxy(i1,i2,i3,axis,0) ) )/(2.*dr(axis))   ! coeff of v(-1)
                       ! here are the wrong ghost point values
                       q(0) = u(i1-is1,i2-is2,i3,u1c)
                       q(1) = u(i1-is1,i2-is2,i3,u2c)
                       ! subtract off the contributions from the wrong values at the ghost points:
                       do n=0,1
                         f(n) = (a2(n,0)*q(0)+a2(n,1)*q(1)) - f(n)
                       end do
                       ! *optimize me* solve this by hand 
                       call dgeco( a2(0,0), 2, 2, ipvt(0),rcond,work(0))
                       call dgesl( a2(0,0), 2, 2, ipvt(0), f(0), job)
                       u(i1-is1,i2-is2,i3,u1c)=f(0)
                       u(i1-is1,i2-is2,i3,u2c)=f(1)
                       if( debug.gt.3 )then ! re-evaluate
                         ux=ux22(i1,i2,i3,u1c)
                         uy=uy22(i1,i2,i3,u1c)
                         vx=ux22(i1,i2,i3,u2c)
                         vy=uy22(i1,i2,i3,u2c) 
                         f(0) = ux + vy                                                - fe(0)
                         f(1) = 2*an1*an2*( ux - vy ) + ( an2**2 -an1**2 )*( uy + vx ) - fe(1)  
                         resMax = max(abs(f(0)),abs(f(1)))      
                         write(*,'(" --> traction BC: curvilinear fill ghost: residuals =",2(1pe12.4))') f(0),f(1)
                       end if
                       ! ------ get second ghost by using a wide formula ----
                       !  axis=0: use
                       !     D0r with dr -> 2*dr
                       !  axis=1: use
                       !     D0s with ds -> 2*ds
                       if( numGhost.gt.1 )then
                         ghost=2;                          ! fill this ghost line 
                         js1=ghost*is1; js2=ghost*is2;     ! shift from boundary point (i1,i2,i3) to ghost point
                         if( axis.eq.0 )then
                           iw1=ghost; iw2=1;               ! use wide stencil for i1
                         else
                           iw1=1;     iw2=ghost;           ! use wide stencil for i2
                         end if
                           urv(u1c) = (u(i1+iw1,i2,i3,u1c) - u(i1-iw1,i2,i3,u1c))/(2.*iw1*dr(0))
                           urv(u2c) = (u(i1+iw1,i2,i3,u2c) - u(i1-iw1,i2,i3,u2c))/(2.*iw1*dr(0))
                           usv(u1c) = (u(i1,i2+iw2,i3,u1c) - u(i1,i2-iw2,i3,u1c))/(2.*iw2*dr(1))
                           usv(u2c) = (u(i1,i2+iw2,i3,u2c) - u(i1,i2-iw2,i3,u2c))/(2.*iw2*dr(1))
                           ux=rsxy(i1,i2,i3,0,0)*urv(u1c) + rsxy(i1,i2,i3,1,0)*usv(u1c)
                           vx=rsxy(i1,i2,i3,0,0)*urv(u2c) + rsxy(i1,i2,i3,1,0)*usv(u2c)
                           uy=rsxy(i1,i2,i3,0,1)*urv(u1c) + rsxy(i1,i2,i3,1,1)*usv(u1c)
                           vy=rsxy(i1,i2,i3,0,1)*urv(u2c) + rsxy(i1,i2,i3,1,1)*usv(u2c)
                         ! --- here are  the equations we mean to satisfy ----
                         !  --> evaluate with the wrong values in the ghost 
                         f(0) = ux + vy                                                - fe(0)
                         f(1) = 2*an1*an2*( ux - vy ) + ( an2**2 -an1**2 )*( uy + vx ) - fe(1)
                         !  [ a2(0,0) a2(0,1) ][ u(-1) ] =  RHS
                         !  [ a2(1,0) a2(1,1) ][ v(-1) ]   
                         ! u_x + v_y = rx*ur + ry*vy  + sx*us + sy*vs 
                         a2(0,0)=-is*( rsxy(i1,i2,i3,axis,0) )/(2.*ghost*dr(axis) )  ! coeff of u(-ghost) in u_x + v_y = 0 
                         a2(0,1)=-is*( rsxy(i1,i2,i3,axis,1) )/(2.*ghost*dr(axis) )  ! coeff of v(-ghost) in u_x + v_y = 0
                         ! 2*an1*an2*( ux - vy ) + (an2**2 -an1**2 )*( uy + vx )
                         !   ux = rx*ur + sx*us   vx = rx*vr + sx*vs   
                         !   uy = ry*ur + sy*us   vy = ry*vr + sy*vs
                         a2(1,0)=-is*( 2*an1*an2*(  rsxy(i1,i2,i3,axis,0) ) + ( an2**2 -an1**2 )*( rsxy(i1,i2,i3,axis,1) ) )/(2.*ghost*dr(axis))  ! coeff of u(-ghost)
                         a2(1,1)=-is*( 2*an1*an2*( -rsxy(i1,i2,i3,axis,1) ) + ( an2**2 -an1**2 )*( rsxy(i1,i2,i3,axis,0) ) )/(2.*ghost*dr(axis))  ! coeff of v(-ghost)
                         ! here are the wrong ghost point values
                         q(0) = u(i1-js1,i2-js2,i3,u1c)
                         q(1) = u(i1-js1,i2-js2,i3,u2c)
                         ! subtract off the contributions from the wrong values at the ghost points:
                         do n=0,1
                           f(n) = (a2(n,0)*q(0)+a2(n,1)*q(1)) - f(n)
                         end do
                         ! *optimize me* solve this by hand 
                         call dgeco( a2(0,0), 2, 2, ipvt(0),rcond,work(0))
                         call dgesl( a2(0,0), 2, 2, ipvt(0), f(0), job)
                         u(i1-js1,i2-js2,i3,u1c)=f(0)
                         u(i1-js1,i2-js2,i3,u2c)=f(1)
                         if( debug.gt.3 )then ! re-evaluate
                             urv(u1c) = (u(i1+iw1,i2,i3,u1c) - u(i1-iw1,i2,i3,u1c))/(2.*iw1*dr(0))
                             urv(u2c) = (u(i1+iw1,i2,i3,u2c) - u(i1-iw1,i2,i3,u2c))/(2.*iw1*dr(0))
                             usv(u1c) = (u(i1,i2+iw2,i3,u1c) - u(i1,i2-iw2,i3,u1c))/(2.*iw2*dr(1))
                             usv(u2c) = (u(i1,i2+iw2,i3,u2c) - u(i1,i2-iw2,i3,u2c))/(2.*iw2*dr(1))
                             ux=rsxy(i1,i2,i3,0,0)*urv(u1c) + rsxy(i1,i2,i3,1,0)*usv(u1c)
                             vx=rsxy(i1,i2,i3,0,0)*urv(u2c) + rsxy(i1,i2,i3,1,0)*usv(u2c)
                             uy=rsxy(i1,i2,i3,0,1)*urv(u1c) + rsxy(i1,i2,i3,1,1)*usv(u1c)
                             vy=rsxy(i1,i2,i3,0,1)*urv(u2c) + rsxy(i1,i2,i3,1,1)*usv(u2c)
                           f(0) = ux + vy                                                - fe(0)
                           f(1) = 2*an1*an2*( ux - vy ) + ( an2**2 -an1**2 )*( uy + vx ) - fe(1)  
                           resMax = max(abs(f(0)),abs(f(1)))      
                           write(*,'(" --> traction BC: curvilinear fill ghost=2 using wide approx: residuals =",2(1pe12.4))') f(0),f(1)
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
                       end do  
                      end do
                      end do
                  end if
                end if
                end do ! end side
                end do ! end axis
              else 
                ! un-known nd and orderOfAccuracy
                stop 6663
              end if
             else if( nd.eq.3 )then
              !    *************************
              !    ********** 3D ***********
              !    *************************
              if( orderOfAccuracy.eq.2 .and. gridType.eq.rectangular )then
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
               !*      ! (js1,js2,js3) used to compute tangential derivatives
               !*      js1=0
               !*      js2=0
               !*      js3=0
               !*      if( axisp1.eq.0 )then
               !*        js1=1-2*side
               !*      else if( axisp1.eq.1 )then
               !*        js2=1-2*side
               !*      else if( axisp1.eq.2 )then
               !*        js3=1-2*side
               !*      else
               !*        stop 5
               !*      end if
               !* 
               !*      ! (ks1,ks2,ks3) used to compute second tangential derivative
               !*      ks1=0
               !*      ks2=0
               !*      ks3=0
               !*      if( axisp2.eq.0 )then
               !*        ks1=1-2*side
               !*      else if( axisp2.eq.1 )then
               !*        ks2=1-2*side
               !*      else if( axisp2.eq.2 )then
               !*        ks3=1-2*side
               !*      else
               !*        stop 5
               !*      end if
                    if( debug.gt.7 )then
                      write(*,'(" bcOpt: grid,side,axis=",3i3,", loop bounds: n1a,n1b,n2a,n2b,n3a,n3b=",6i3)') grid,side,axis,n1a,n1b,n2a,n2b,n3a,n3b
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
                  if( addBoundaryForcing(side,axis).eq.0 )then
                    stop 333
                     do i3=nn3a,nn3b
                     do i2=nn2a,nn2b
                     do i1=nn1a,nn1b
                        u(i1,i2,i3,u1c)=0.
                        u(i1,i2,i3,u2c)=0.
                        u(i1,i2,i3,u3c)=0.
                        u(i1-is1,i2-is2,i3-is3,u1c)=(3.*u(i1,i2,i3,u1c)-3.*u(i1+is1,i2+is2,i3+is3,u1c)+u(i1+2*is1,i2+2*is2,i3+2*is3,u1c))
                        u(i1-is1,i2-is2,i3-is3,u2c)=(3.*u(i1,i2,i3,u2c)-3.*u(i1+is1,i2+is2,i3+is3,u2c)+u(i1+2*is1,i2+2*is2,i3+2*is3,u2c))
                        u(i1-is1,i2-is2,i3-is3,u3c)=(3.*u(i1,i2,i3,u3c)-3.*u(i1+is1,i2+is2,i3+is3,u3c)+u(i1+2*is1,i2+2*is2,i3+2*is3,u3c))
                     end do
                     end do
                     end do
                  else
                    stop 333
                     do i3=nn3a,nn3b
                     do i2=nn2a,nn2b
                     do i1=nn1a,nn1b
                         call ogf3d(ep,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,u0,v0,w0)
                        u(i1,i2,i3,u1c)=u0
                        u(i1,i2,i3,u2c)=v0
                        u(i1,i2,i3,u3c)=w0
                        u(i1-is1,i2-is2,i3-is3,u1c)=(3.*u(i1,i2,i3,u1c)-3.*u(i1+is1,i2+is2,i3+is3,u1c)+u(i1+2*is1,i2+2*is2,i3+2*is3,u1c))
                        u(i1-is1,i2-is2,i3-is3,u2c)=(3.*u(i1,i2,i3,u2c)-3.*u(i1+is1,i2+is2,i3+is3,u2c)+u(i1+2*is1,i2+2*is2,i3+2*is3,u2c))
                        u(i1-is1,i2-is2,i3-is3,u3c)=(3.*u(i1,i2,i3,u3c)-3.*u(i1+is1,i2+is2,i3+is3,u3c)+u(i1+2*is1,i2+2*is2,i3+2*is3,u3c))
                     end do
                     end do
                     end do
                  end if
                else if( boundaryCondition(side,axis).eq.tractionBC )then 
                  ! first extrap values to ghost points (may be needed at corners)
                   do i3=nn3a,nn3b
                   do i2=nn2a,nn2b
                   do i1=nn1a,nn1b
                    u(i1-is1,i2-is2,i3-is3,u1c)=(3.*u(i1,i2,i3,u1c)-3.*u(i1+is1,i2+is2,i3+is3,u1c)+u(i1+2*is1,i2+2*is2,i3+2*is3,u1c))
                    u(i1-is1,i2-is2,i3-is3,u2c)=(3.*u(i1,i2,i3,u2c)-3.*u(i1+is1,i2+is2,i3+is3,u2c)+u(i1+2*is1,i2+2*is2,i3+2*is3,u2c))
                    u(i1-is1,i2-is2,i3-is3,u3c)=(3.*u(i1,i2,i3,u3c)-3.*u(i1+is1,i2+is2,i3+is3,u3c)+u(i1+2*is1,i2+2*is2,i3+2*is3,u3c))
                   end do
                   end do
                   end do
                else if( boundaryCondition(side,axis).eq.dirichletBoundaryCondition .or.boundaryCondition(side,axis).eq.symmetry )then
                  ! do nothing here
                else if( boundaryCondition(side,axis).gt.0 )then
                  stop 1193
                end if
                end do ! end side
                end do ! end axis
               ! ** now apply BC's that assign the ghost values *********
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
               !*      ! (js1,js2,js3) used to compute tangential derivatives
               !*      js1=0
               !*      js2=0
               !*      js3=0
               !*      if( axisp1.eq.0 )then
               !*        js1=1-2*side
               !*      else if( axisp1.eq.1 )then
               !*        js2=1-2*side
               !*      else if( axisp1.eq.2 )then
               !*        js3=1-2*side
               !*      else
               !*        stop 5
               !*      end if
               !* 
               !*      ! (ks1,ks2,ks3) used to compute second tangential derivative
               !*      ks1=0
               !*      ks2=0
               !*      ks3=0
               !*      if( axisp2.eq.0 )then
               !*        ks1=1-2*side
               !*      else if( axisp2.eq.1 )then
               !*        ks2=1-2*side
               !*      else if( axisp2.eq.2 )then
               !*        ks3=1-2*side
               !*      else
               !*        stop 5
               !*      end if
                    if( debug.gt.7 )then
                      write(*,'(" bcOpt: grid,side,axis=",3i3,", loop bounds: n1a,n1b,n2a,n2b,n3a,n3b=",6i3)') grid,side,axis,n1a,n1b,n2a,n2b,n3a,n3b
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
                if( boundaryCondition(side,axis).eq.tractionBC )then 
                  if( addBoundaryForcing(side,axis).eq.0 )then
                   alpha=lambda/(lambda+2.*mu)
                   beta=1./(lambda+2.*mu)
                   if( axis.eq.0 )then
                     ! u.x = -alpha*(v.y+w.z)
                     ! v.x = -u.y  
                     ! w.x = -u.z
                     do i3=n3a,n3b
                     do i2=n2a,n2b
                     do i1=n1a,n1b
                     if( mask(i1,i2,i3).gt.0 )then
                     vy=uy23r(i1,i2,i3,u2c)
                     wz=uz23r(i1,i2,i3,u3c)
                     uy=uy23r(i1,i2,i3,u1c)
                     uz=uz23r(i1,i2,i3,u1c)
                      u(i1-is1,i2-is2,i3-is3,u1c)=u(i1+is1,i2+is2,i3+is3,u1c)-is1*dx(0)*2.*(-alpha*(vy+wz))
                      u(i1-is1,i2-is2,i3-is3,u2c)=u(i1+is1,i2+is2,i3+is3,u2c)-is1*dx(0)*2.*(-uy)
                      u(i1-is1,i2-is2,i3-is3,u3c)=u(i1+is1,i2+is2,i3+is3,u3c)-is1*dx(0)*2.*(-uz)
                     end if
                     end do
                     end do
                     end do
                   else if( axis.eq.1 )then
                   ! u.y = - v.x
                   ! v.y = -alpha*(u.x+w.z)
                   ! w.y = - v.z
                     do i3=n3a,n3b
                     do i2=n2a,n2b
                     do i1=n1a,n1b
                     if( mask(i1,i2,i3).gt.0 )then
                     vx=ux23r(i1,i2,i3,u2c)
                     ux=ux23r(i1,i2,i3,u1c)
                     wz=uz23r(i1,i2,i3,u3c)
                     vz=uz23r(i1,i2,i3,u2c)
                      u(i1-is1,i2-is2,i3-is3,u1c)=u(i1+is1,i2+is2,i3+is3,u1c)-is2*dx(1)*2.*(-vx)
                      u(i1-is1,i2-is2,i3-is3,u2c)=u(i1+is1,i2+is2,i3+is3,u2c)-is2*dx(1)*2.*(-alpha*(ux+wz))
                      u(i1-is1,i2-is2,i3-is3,u3c)=u(i1+is1,i2+is2,i3+is3,u3c)-is2*dx(1)*2.*(-vz)
                     end if
                     end do
                     end do
                     end do
                   else 
                   ! u.z = - w.x
                   ! v.z = - w.y
                   ! w.z = -alpha*(u.x+v.y)
                     do i3=n3a,n3b
                     do i2=n2a,n2b
                     do i1=n1a,n1b
                     if( mask(i1,i2,i3).gt.0 )then
                     wx=ux23r(i1,i2,i3,u3c)
                     wy=uy23r(i1,i2,i3,u3c)
                     ux=ux23r(i1,i2,i3,u1c)
                     vy=uy23r(i1,i2,i3,u2c)
                      u(i1-is1,i2-is2,i3-is3,u1c)=u(i1+is1,i2+is2,i3+is3,u1c)-is3*dx(2)*2.*(-wx)
                      u(i1-is1,i2-is2,i3-is3,u2c)=u(i1+is1,i2+is2,i3+is3,u2c)-is3*dx(2)*2.*(-wy)
                      u(i1-is1,i2-is2,i3-is3,u3c)=u(i1+is1,i2+is2,i3+is3,u3c)-is3*dx(2)*2.*(-alpha*(ux+vy))
                     end if
                     end do
                     end do
                     end do
                   end if      
                  else
                   alpha=lambda/(lambda+2.*mu)
                   beta=1./(lambda+2.*mu)
                   if( axis.eq.0 )then
                     ! u.x = -alpha*(v.y+w.z)
                     ! v.x = -u.y  
                     ! w.x = -u.z
                     do i3=n3a,n3b
                     do i2=n2a,n2b
                     do i1=n1a,n1b
                     if( mask(i1,i2,i3).gt.0 )then
                     vy=uy23r(i1,i2,i3,u2c)
                     wz=uz23r(i1,i2,i3,u3c)
                     uy=uy23r(i1,i2,i3,u1c)
                     uz=uz23r(i1,i2,i3,u1c)
                      if( assignTwilightZone.eq.0 )then
                       u(i1-is1,i2-is2,i3-is3,u1c)=u(i1+is1,i2+is2,i3+is3,u1c)+dx(0)*2.*(is1*alpha*(vy+wz)+ beta*bcf(side,axis,i1,i2,i3,u1c) )
                       u(i1-is1,i2-is2,i3-is3,u2c)=u(i1+is1,i2+is2,i3+is3,u2c)+dx(0)*2.*(is1*uy+         (1./mu)*bcf(side,axis,i1,i2,i3,u2c) )
                       u(i1-is1,i2-is2,i3-is3,u3c)=u(i1+is1,i2+is2,i3+is3,u3c)+dx(0)*2.*(is1*uz+         (1./mu)*bcf(side,axis,i1,i2,i3,u3c) )
                      else
                         call ogDeriv3(ep, 0,1,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, u1c,ux0, u2c,vx0, u3c,wx0)
                         call ogDeriv3(ep, 0,0,1,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, u1c,uy0, u2c,vy0, u3c,wy0)
                         call ogDeriv3(ep, 0,0,0,1, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, u1c,uz0, u2c,vz0, u3c,wz0)
                       u(i1-is1,i2-is2,i3-is3,u1c)=u(i1+is1,i2+is2,i3+is3,u1c)-is1*dx(0)*2.*(-alpha*(vy+wz)+ux0+alpha*(vy0+wz0))
                       u(i1-is1,i2-is2,i3-is3,u2c)=u(i1+is1,i2+is2,i3+is3,u2c)-is1*dx(0)*2.*(-uy           +vx0+uy0)
                       u(i1-is1,i2-is2,i3-is3,u3c)=u(i1+is1,i2+is2,i3+is3,u3c)-is1*dx(0)*2.*(-uz           +wx0+uz0)
                      end if
                     end if
                     end do
                     end do
                     end do
                   else if( axis.eq.1 )then
                   ! u.y = - v.x
                   ! v.y = -alpha*(u.x+w.z)
                   ! w.y = - v.z
                     do i3=n3a,n3b
                     do i2=n2a,n2b
                     do i1=n1a,n1b
                     if( mask(i1,i2,i3).gt.0 )then
                     vx=ux23r(i1,i2,i3,u2c)
                     ux=ux23r(i1,i2,i3,u1c)
                     wz=uz23r(i1,i2,i3,u3c)
                     vz=uz23r(i1,i2,i3,u2c)
                      if( assignTwilightZone.eq.0 )then
                       u(i1-is1,i2-is2,i3-is3,u1c)=u(i1+is1,i2+is2,i3+is3,u1c)+dx(1)*2.*(is2*vx +         (1./mu)*bcf(side,axis,i1,i2,i3,u1c))
                       u(i1-is1,i2-is2,i3-is3,u2c)=u(i1+is1,i2+is2,i3+is3,u2c)+dx(1)*2.*(is2*alpha*(ux+wz) + beta*bcf(side,axis,i1,i2,i3,u2c))
                       u(i1-is1,i2-is2,i3-is3,u3c)=u(i1+is1,i2+is2,i3+is3,u3c)+dx(1)*2.*(is2*vz +         (1./mu)*bcf(side,axis,i1,i2,i3,u3c) )
                      else
                         call ogDeriv3(ep, 0,1,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, u1c,ux0, u2c,vx0, u3c,wx0)
                         call ogDeriv3(ep, 0,0,1,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, u1c,uy0, u2c,vy0, u3c,wy0)
                         call ogDeriv3(ep, 0,0,0,1, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, u1c,uz0, u2c,vz0, u3c,wz0)
                       u(i1-is1,i2-is2,i3-is3,u1c)=u(i1+is1,i2+is2,i3+is3,u1c)-is2*dx(1)*2.*(-vx            +uy0+vx0)
                       u(i1-is1,i2-is2,i3-is3,u2c)=u(i1+is1,i2+is2,i3+is3,u2c)-is2*dx(1)*2.*(-alpha*(ux+wz) +vy0+alpha*(ux0+wz0))
                       u(i1-is1,i2-is2,i3-is3,u3c)=u(i1+is1,i2+is2,i3+is3,u3c)-is2*dx(1)*2.*(-vz            +wy0+vz0)
                      end if
                     end if
                     end do
                     end do
                     end do
                   else 
                   ! u.z = - w.x
                   ! v.z = - w.y
                   ! w.z = -alpha*(u.x+v.y)
                     do i3=n3a,n3b
                     do i2=n2a,n2b
                     do i1=n1a,n1b
                     if( mask(i1,i2,i3).gt.0 )then
                     wx=ux23r(i1,i2,i3,u3c)
                     wy=uy23r(i1,i2,i3,u3c)
                     ux=ux23r(i1,i2,i3,u1c)
                     vy=uy23r(i1,i2,i3,u2c)
                      if( assignTwilightZone.eq.0 )then
                       u(i1-is1,i2-is2,i3-is3,u1c)=u(i1+is1,i2+is2,i3+is3,u1c)+dx(2)*2.*( is3*wx +         (1./mu)*bcf(side,axis,i1,i2,i3,u1c))
                       u(i1-is1,i2-is2,i3-is3,u2c)=u(i1+is1,i2+is2,i3+is3,u2c)+dx(2)*2.*( is3*wy +         (1./mu)*bcf(side,axis,i1,i2,i3,u2c))
                       u(i1-is1,i2-is2,i3-is3,u3c)=u(i1+is1,i2+is2,i3+is3,u3c)+dx(2)*2.*(is3*alpha*(ux+vy) + beta*bcf(side,axis,i1,i2,i3,u3c))
                      else
                         call ogDeriv3(ep, 0,1,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, u1c,ux0, u2c,vx0, u3c,wx0)
                         call ogDeriv3(ep, 0,0,1,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, u1c,uy0, u2c,vy0, u3c,wy0)
                         call ogDeriv3(ep, 0,0,0,1, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, u1c,uz0, u2c,vz0, u3c,wz0)
                       u(i1-is1,i2-is2,i3-is3,u1c)=u(i1+is1,i2+is2,i3+is3,u1c)-is3*dx(2)*2.*(-wx            +uz0+wx0)
                       u(i1-is1,i2-is2,i3-is3,u2c)=u(i1+is1,i2+is2,i3+is3,u2c)-is3*dx(2)*2.*(-wy            +vz0+wy0)
                       u(i1-is1,i2-is2,i3-is3,u3c)=u(i1+is1,i2+is2,i3+is3,u3c)-is3*dx(2)*2.*(-alpha*(ux+vy) +wz0+alpha*(ux0+vy0))
                      end if
                     end if
                     end do
                     end do
                     end do
                   end if      
                  end if
                end if
                end do ! end side
                end do ! end axis
              else if( orderOfAccuracy.eq.2 .and. gridType.eq.curvilinear )then
               ! *********************************************
               ! ************* 3d Curvilinear ****************
               ! *********************************************
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
               !*      ! (js1,js2,js3) used to compute tangential derivatives
               !*      js1=0
               !*      js2=0
               !*      js3=0
               !*      if( axisp1.eq.0 )then
               !*        js1=1-2*side
               !*      else if( axisp1.eq.1 )then
               !*        js2=1-2*side
               !*      else if( axisp1.eq.2 )then
               !*        js3=1-2*side
               !*      else
               !*        stop 5
               !*      end if
               !* 
               !*      ! (ks1,ks2,ks3) used to compute second tangential derivative
               !*      ks1=0
               !*      ks2=0
               !*      ks3=0
               !*      if( axisp2.eq.0 )then
               !*        ks1=1-2*side
               !*      else if( axisp2.eq.1 )then
               !*        ks2=1-2*side
               !*      else if( axisp2.eq.2 )then
               !*        ks3=1-2*side
               !*      else
               !*        stop 5
               !*      end if
                    if( debug.gt.7 )then
                      write(*,'(" bcOpt: grid,side,axis=",3i3,", loop bounds: n1a,n1b,n2a,n2b,n3a,n3b=",6i3)') grid,side,axis,n1a,n1b,n2a,n2b,n3a,n3b
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
                  ! note: we can assign ghost pts in tangential dir too:
                  if( addBoundaryForcing(side,axis).eq.0 )then
                    stop 333
                     do i3=nn3a,nn3b
                     do i2=nn2a,nn2b
                     do i1=nn1a,nn1b
                        u(i1,i2,i3,u1c)=0.
                        u(i1,i2,i3,u2c)=0.
                        u(i1,i2,i3,u3c)=0.
                        u(i1-is1,i2-is2,i3-is3,u1c)=(3.*u(i1,i2,i3,u1c)-3.*u(i1+is1,i2+is2,i3+is3,u1c)+u(i1+2*is1,i2+2*is2,i3+2*is3,u1c))
                        u(i1-is1,i2-is2,i3-is3,u2c)=(3.*u(i1,i2,i3,u2c)-3.*u(i1+is1,i2+is2,i3+is3,u2c)+u(i1+2*is1,i2+2*is2,i3+2*is3,u2c))
                        u(i1-is1,i2-is2,i3-is3,u3c)=(3.*u(i1,i2,i3,u3c)-3.*u(i1+is1,i2+is2,i3+is3,u3c)+u(i1+2*is1,i2+2*is2,i3+2*is3,u3c))
                     end do
                     end do
                     end do
                  else
                    stop 333
                     do i3=nn3a,nn3b
                     do i2=nn2a,nn2b
                     do i1=nn1a,nn1b
                         call ogf3d(ep,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,u0,v0,w0)
                        u(i1,i2,i3,u1c)=u0
                        u(i1,i2,i3,u2c)=v0
                        u(i1,i2,i3,u3c)=w0
                        u(i1-is1,i2-is2,i3-is3,u1c)=(3.*u(i1,i2,i3,u1c)-3.*u(i1+is1,i2+is2,i3+is3,u1c)+u(i1+2*is1,i2+2*is2,i3+2*is3,u1c))
                        u(i1-is1,i2-is2,i3-is3,u2c)=(3.*u(i1,i2,i3,u2c)-3.*u(i1+is1,i2+is2,i3+is3,u2c)+u(i1+2*is1,i2+2*is2,i3+2*is3,u2c))
                        u(i1-is1,i2-is2,i3-is3,u3c)=(3.*u(i1,i2,i3,u3c)-3.*u(i1+is1,i2+is2,i3+is3,u3c)+u(i1+2*is1,i2+2*is2,i3+2*is3,u3c))
                     end do
                     end do
                     end do
                  end if
                 else if( boundaryCondition(side,axis).eq.tractionBC )then 
                  ! first extrap values to ghost points (may be needed at corners)
                   do i3=nn3a,nn3b
                   do i2=nn2a,nn2b
                   do i1=nn1a,nn1b
                    u(i1-is1,i2-is2,i3-is3,u1c)=(3.*u(i1,i2,i3,u1c)-3.*u(i1+is1,i2+is2,i3+is3,u1c)+u(i1+2*is1,i2+2*is2,i3+2*is3,u1c))
                    u(i1-is1,i2-is2,i3-is3,u2c)=(3.*u(i1,i2,i3,u2c)-3.*u(i1+is1,i2+is2,i3+is3,u2c)+u(i1+2*is1,i2+2*is2,i3+2*is3,u2c))
                    u(i1-is1,i2-is2,i3-is3,u3c)=(3.*u(i1,i2,i3,u3c)-3.*u(i1+is1,i2+is2,i3+is3,u3c)+u(i1+2*is1,i2+2*is2,i3+2*is3,u3c))
                   end do
                   end do
                   end do
                 else if( boundaryCondition(side,axis).eq.dirichletBoundaryCondition .or.boundaryCondition(side,axis).eq.symmetry )then
                   write(*,*) 'finish me BC'
                   stop 3333
                 else if( boundaryCondition(side,axis).gt.0 )then
                  stop 1193
                 end if
                end do ! end side
                end do ! end axis
               ! ** now apply BC's that assign the ghost values *********
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
               !*      ! (js1,js2,js3) used to compute tangential derivatives
               !*      js1=0
               !*      js2=0
               !*      js3=0
               !*      if( axisp1.eq.0 )then
               !*        js1=1-2*side
               !*      else if( axisp1.eq.1 )then
               !*        js2=1-2*side
               !*      else if( axisp1.eq.2 )then
               !*        js3=1-2*side
               !*      else
               !*        stop 5
               !*      end if
               !* 
               !*      ! (ks1,ks2,ks3) used to compute second tangential derivative
               !*      ks1=0
               !*      ks2=0
               !*      ks3=0
               !*      if( axisp2.eq.0 )then
               !*        ks1=1-2*side
               !*      else if( axisp2.eq.1 )then
               !*        ks2=1-2*side
               !*      else if( axisp2.eq.2 )then
               !*        ks3=1-2*side
               !*      else
               !*        stop 5
               !*      end if
                    if( debug.gt.7 )then
                      write(*,'(" bcOpt: grid,side,axis=",3i3,", loop bounds: n1a,n1b,n2a,n2b,n3a,n3b=",6i3)') grid,side,axis,n1a,n1b,n2a,n2b,n3a,n3b
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
                if( boundaryCondition(side,axis).eq.tractionBC )then 
                  if( addBoundaryForcing(side,axis).eq.0 )then
                     write(*,*) 'bcOptIsm -- finish me for 3d curvlinear'
                     stop 6678
                   alpha=lambda+2*mu
                    do i3=n3a,n3b
                    do i2=n2a,n2b
                    do i1=n1a,n1b
                    if( mask(i1,i2,i3).gt.0 )then
                    ! Solve n.tauv = 0 
                    !    -->   A uv.r + B uv.s = 0
                    ! here is the normal (assumed to be the same on both sides)
                    an1=rsxy(i1,i2,i3,axis,0)   ! normal (an1,an2,an3)
                    an2=rsxy(i1,i2,i3,axis,1)
                    an3=rsxy(i1,i2,i3,axis,2)
                    aNorm=max(epsx,sqrt(an1**2+an2**2+an3**2))
                    an1=an1/aNorm
                    an2=an2/aNorm
                    an3=an3/aNorm
                    ux=ux23(i1,i2,i3,u1c)
                    uy=uy23(i1,i2,i3,u1c)
                    uz=uz23(i1,i2,i3,u1c)
                    vx=ux23(i1,i2,i3,u2c)
                    vy=uy23(i1,i2,i3,u2c)
                    vz=uz23(i1,i2,i3,u2c)
                    wx=ux23(i1,i2,i3,u3c)
                    wy=uy23(i1,i2,i3,u3c)
                    wz=uz23(i1,i2,i3,u3c)
                    tau11 = alpha*ux + lambda*(vy+wz)
                    tau12 = mu*( uy + vx )
                    tau13 = mu*( uz + wx )
                    tau21 = tau12
                    tau22 = lambda*(ux+wz) + alpha*vy
                    tau23 = mu*( vz + wy )
                    tau31 = tau13
                    tau32 = tau23
                    tau33 = lambda*(ux+vy) + alpha*wz
                    ! here are  the equations we mean to satisfy:  
                    f(0) = an1*tau11+an2*tau21+an3*tau31
                    f(1) = an1*tau12+an2*tau22+an3*tau32
                    f(2) = an1*tau13+an2*tau23+an3*tau33
                    !  [ a3(0,0) a3(0,1) a3(0,2) ][ u(-1) ] =  RHS
                    !  [ a3(1,0) a3(1,1) a3(1,2) ][ v(-1) ]   
                    !  [ a3(2,0) a3(2,1) a3(2,2) ][ w(-1) ]   
                    a3(0,0)=-is*( an1*alpha* rsxy(i1,i2,i3,axis,0)+an2*mu*    rsxy(i1,i2,i3,axis,1)+an3*mu*    rsxy(i1,i2,i3,axis,2) )/(2.*dr(axis))
                    a3(0,1)=-is*( an1*lambda*rsxy(i1,i2,i3,axis,1)+an2*mu*    rsxy(i1,i2,i3,axis,0) )/(2.*dr(axis))
                    a3(0,2)=-is*( an1*lambda*rsxy(i1,i2,i3,axis,2)+an3*mu*    rsxy(i1,i2,i3,axis,0) )/(2.*dr(axis))
                    a3(1,0)=-is*( an1*mu*    rsxy(i1,i2,i3,axis,1)+an2*lambda*rsxy(i1,i2,i3,axis,0) )/(2.*dr(axis))
                    a3(1,1)=-is*( an1*mu*    rsxy(i1,i2,i3,axis,0)+an2*alpha* rsxy(i1,i2,i3,axis,1)+an3*mu*    rsxy(i1,i2,i3,axis,2) )/(2.*dr(axis))
                    a3(1,2)=-is*( an2*lambda*rsxy(i1,i2,i3,axis,2)+an3*mu*    rsxy(i1,i2,i3,axis,1) )/(2.*dr(axis))
                    a3(2,0)=-is*( an1*mu*    rsxy(i1,i2,i3,axis,2)+an3*lambda*rsxy(i1,i2,i3,axis,0) )/(2.*dr(axis))
                    a3(2,1)=-is*( an2*mu*    rsxy(i1,i2,i3,axis,2)+an3*lambda*rsxy(i1,i2,i3,axis,1) )/(2.*dr(axis))
                    a3(2,2)=-is*( an1*mu*    rsxy(i1,i2,i3,axis,0)+an2*mu*    rsxy(i1,i2,i3,axis,1)+an3*alpha* rsxy(i1,i2,i3,axis,2) )/(2.*dr(axis))
                    ! here are the wrong ghostpoint values
                    q(0) = u(i1-is1,i2-is2,i3-is3,u1c)
                    q(1) = u(i1-is1,i2-is2,i3-is3,u2c)
                    q(2) = u(i1-is1,i2-is2,i3-is3,u3c)
                    ! subtract off the contributions from the wrong values at the ghost points:
                    do n=0,2
                      f(n) = (a3(n,0)*q(0)+a3(n,1)*q(1)+a3(n,2)*q(2)) - f(n)
                    end do
                    call dgeco( a3(0,0), 3, 3, ipvt(0),rcond,work(0))
                    call dgesl( a3(0,0), 3, 3, ipvt(0), f(0), job)
                    u(i1-is1,i2-is2,i3-is3,u1c)=f(0)
                    u(i1-is1,i2-is2,i3-is3,u2c)=f(1)
                    u(i1-is1,i2-is2,i3-is3,u3c)=f(2)
                    if( debug.gt.0 )then ! re-evaluate
                      ux=ux23(i1,i2,i3,u1c)
                      uy=uy23(i1,i2,i3,u1c)
                      uz=uz23(i1,i2,i3,u1c)
                      vx=ux23(i1,i2,i3,u2c)
                      vy=uy23(i1,i2,i3,u2c)
                      vz=uz23(i1,i2,i3,u2c)
                      wx=ux23(i1,i2,i3,u3c)
                      wy=uy23(i1,i2,i3,u3c)
                      wz=uz23(i1,i2,i3,u3c)
                      tau11 = alpha*ux + lambda*(vy+wz)
                      tau12 = mu*( uy + vx )
                      tau13 = mu*( uz + wx )
                      tau21 = tau12
                      tau22 = lambda*(ux+wz) + alpha*vy
                      tau23 = mu*( vz + wy )
                      tau31 = tau13
                      tau32 = tau23
                      tau33 = lambda*(ux+vy) + alpha*wz
                      ! here are  the equations we mean to satisfy:  
                      f(0) = an1*tau11+an2*tau21+an3*tau31
                      f(1) = an1*tau12+an2*tau22+an3*tau32
                      f(2) = an1*tau13+an2*tau23+an3*tau33
                      write(*,'(" --> bc: (",i1,",",i1,") i1,i2,i3=",3i4," n.tau=",4e10.2)') side,axis,i1,i2,i3,f(0),f(1),f(2)
                        ! '
                    end if
                    end if
                    end do
                    end do
                    end do
                  else
                     write(*,*) 'bcOptIsm -- finish me for 3d curvlinear'
                     stop 6678
                   alpha=lambda+2*mu
                    do i3=n3a,n3b
                    do i2=n2a,n2b
                    do i1=n1a,n1b
                    if( mask(i1,i2,i3).gt.0 )then
                    ! Solve n.tauv = 0 
                    !    -->   A uv.r + B uv.s = 0
                    ! here is the normal (assumed to be the same on both sides)
                    an1=rsxy(i1,i2,i3,axis,0)   ! normal (an1,an2,an3)
                    an2=rsxy(i1,i2,i3,axis,1)
                    an3=rsxy(i1,i2,i3,axis,2)
                    aNorm=max(epsx,sqrt(an1**2+an2**2+an3**2))
                    an1=an1/aNorm
                    an2=an2/aNorm
                    an3=an3/aNorm
                    ux=ux23(i1,i2,i3,u1c)
                    uy=uy23(i1,i2,i3,u1c)
                    uz=uz23(i1,i2,i3,u1c)
                    vx=ux23(i1,i2,i3,u2c)
                    vy=uy23(i1,i2,i3,u2c)
                    vz=uz23(i1,i2,i3,u2c)
                    wx=ux23(i1,i2,i3,u3c)
                    wy=uy23(i1,i2,i3,u3c)
                    wz=uz23(i1,i2,i3,u3c)
                    tau11 = alpha*ux + lambda*(vy+wz)
                    tau12 = mu*( uy + vx )
                    tau13 = mu*( uz + wx )
                    tau21 = tau12
                    tau22 = lambda*(ux+wz) + alpha*vy
                    tau23 = mu*( vz + wy )
                    tau31 = tau13
                    tau32 = tau23
                    tau33 = lambda*(ux+vy) + alpha*wz
                    ! here are  the equations we mean to satisfy:  
                    f(0) = an1*tau11+an2*tau21+an3*tau31
                    f(1) = an1*tau12+an2*tau22+an3*tau32
                    f(2) = an1*tau13+an2*tau23+an3*tau33
                     if( assignTwilightZone.eq.0 )then
                      ! forced case: solve n.tau - f = 0 
                      ! *wdh* 080523 multiply bcf by is since (an1,an2) is not the outward normal
                      f(0) = f(0) + is*bcf(side,axis,i1,i2,i3,u1c)
                      f(1) = f(1) + is*bcf(side,axis,i1,i2,i3,u2c)
                      f(2) = f(2) + is*bcf(side,axis,i1,i2,i3,u3c)
                     else
                        call ogDeriv3(ep, 0,1,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, u1c,ux0, u2c,vx0, u3c,wx0)
                        call ogDeriv3(ep, 0,0,1,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, u1c,uy0, u2c,vy0, u3c,wy0)
                        call ogDeriv3(ep, 0,0,0,1, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, u1c,uz0, u2c,vz0, u3c,wz0)
                      f(0) = f(0) - ( an1*(alpha*ux0+lambda*(vy0+wz0))+an2*(mu*(uy0+vx0))+an3*(mu*(uz0+wx0)) )
                      f(1) = f(1) - ( an1*(mu*(uy0+vx0))+an2*(lambda*(ux0+wz0)+alpha*vy0)+an3*(mu*(vz0+wy0)) )
                      f(2) = f(2) - ( an1*(mu*(uz0+wx0))+an2*(mu*(vz0+wy0))+an3*(lambda*(ux0+vy0)+alpha*wz0) )
                     end if
                    !  [ a3(0,0) a3(0,1) a3(0,2) ][ u(-1) ] =  RHS
                    !  [ a3(1,0) a3(1,1) a3(1,2) ][ v(-1) ]   
                    !  [ a3(2,0) a3(2,1) a3(2,2) ][ w(-1) ]   
                    a3(0,0)=-is*( an1*alpha* rsxy(i1,i2,i3,axis,0)+an2*mu*    rsxy(i1,i2,i3,axis,1)+an3*mu*    rsxy(i1,i2,i3,axis,2) )/(2.*dr(axis))
                    a3(0,1)=-is*( an1*lambda*rsxy(i1,i2,i3,axis,1)+an2*mu*    rsxy(i1,i2,i3,axis,0) )/(2.*dr(axis))
                    a3(0,2)=-is*( an1*lambda*rsxy(i1,i2,i3,axis,2)+an3*mu*    rsxy(i1,i2,i3,axis,0) )/(2.*dr(axis))
                    a3(1,0)=-is*( an1*mu*    rsxy(i1,i2,i3,axis,1)+an2*lambda*rsxy(i1,i2,i3,axis,0) )/(2.*dr(axis))
                    a3(1,1)=-is*( an1*mu*    rsxy(i1,i2,i3,axis,0)+an2*alpha* rsxy(i1,i2,i3,axis,1)+an3*mu*    rsxy(i1,i2,i3,axis,2) )/(2.*dr(axis))
                    a3(1,2)=-is*( an2*lambda*rsxy(i1,i2,i3,axis,2)+an3*mu*    rsxy(i1,i2,i3,axis,1) )/(2.*dr(axis))
                    a3(2,0)=-is*( an1*mu*    rsxy(i1,i2,i3,axis,2)+an3*lambda*rsxy(i1,i2,i3,axis,0) )/(2.*dr(axis))
                    a3(2,1)=-is*( an2*mu*    rsxy(i1,i2,i3,axis,2)+an3*lambda*rsxy(i1,i2,i3,axis,1) )/(2.*dr(axis))
                    a3(2,2)=-is*( an1*mu*    rsxy(i1,i2,i3,axis,0)+an2*mu*    rsxy(i1,i2,i3,axis,1)+an3*alpha* rsxy(i1,i2,i3,axis,2) )/(2.*dr(axis))
                    ! here are the wrong ghostpoint values
                    q(0) = u(i1-is1,i2-is2,i3-is3,u1c)
                    q(1) = u(i1-is1,i2-is2,i3-is3,u2c)
                    q(2) = u(i1-is1,i2-is2,i3-is3,u3c)
                    ! subtract off the contributions from the wrong values at the ghost points:
                    do n=0,2
                      f(n) = (a3(n,0)*q(0)+a3(n,1)*q(1)+a3(n,2)*q(2)) - f(n)
                    end do
                    call dgeco( a3(0,0), 3, 3, ipvt(0),rcond,work(0))
                    call dgesl( a3(0,0), 3, 3, ipvt(0), f(0), job)
                    u(i1-is1,i2-is2,i3-is3,u1c)=f(0)
                    u(i1-is1,i2-is2,i3-is3,u2c)=f(1)
                    u(i1-is1,i2-is2,i3-is3,u3c)=f(2)
                    if( debug.gt.0 )then ! re-evaluate
                      ux=ux23(i1,i2,i3,u1c)
                      uy=uy23(i1,i2,i3,u1c)
                      uz=uz23(i1,i2,i3,u1c)
                      vx=ux23(i1,i2,i3,u2c)
                      vy=uy23(i1,i2,i3,u2c)
                      vz=uz23(i1,i2,i3,u2c)
                      wx=ux23(i1,i2,i3,u3c)
                      wy=uy23(i1,i2,i3,u3c)
                      wz=uz23(i1,i2,i3,u3c)
                      tau11 = alpha*ux + lambda*(vy+wz)
                      tau12 = mu*( uy + vx )
                      tau13 = mu*( uz + wx )
                      tau21 = tau12
                      tau22 = lambda*(ux+wz) + alpha*vy
                      tau23 = mu*( vz + wy )
                      tau31 = tau13
                      tau32 = tau23
                      tau33 = lambda*(ux+vy) + alpha*wz
                      ! here are  the equations we mean to satisfy:  
                      f(0) = an1*tau11+an2*tau21+an3*tau31
                      f(1) = an1*tau12+an2*tau22+an3*tau32
                      f(2) = an1*tau13+an2*tau23+an3*tau33
                       if( assignTwilightZone.eq.0 )then
                        ! forced case: solve n.tau - f = 0 
                        ! *wdh* 080523 multiply bcf by is since (an1,an2) is not the outward normal
                        f(0) = f(0) + is*bcf(side,axis,i1,i2,i3,u1c)
                        f(1) = f(1) + is*bcf(side,axis,i1,i2,i3,u2c)
                        f(2) = f(2) + is*bcf(side,axis,i1,i2,i3,u3c)
                       else
                          call ogDeriv3(ep, 0,1,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, u1c,ux0, u2c,vx0, u3c,wx0)
                          call ogDeriv3(ep, 0,0,1,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, u1c,uy0, u2c,vy0, u3c,wy0)
                          call ogDeriv3(ep, 0,0,0,1, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, u1c,uz0, u2c,vz0, u3c,wz0)
                        f(0) = f(0) - ( an1*(alpha*ux0+lambda*(vy0+wz0))+an2*(mu*(uy0+vx0))+an3*(mu*(uz0+wx0)) )
                        f(1) = f(1) - ( an1*(mu*(uy0+vx0))+an2*(lambda*(ux0+wz0)+alpha*vy0)+an3*(mu*(vz0+wy0)) )
                        f(2) = f(2) - ( an1*(mu*(uz0+wx0))+an2*(mu*(vz0+wy0))+an3*(lambda*(ux0+vy0)+alpha*wz0) )
                       end if 
                      write(*,'(" --> bc: (",i1,",",i1,") i1,i2,i3=",3i4," n.tau=",4e10.2)') side,axis,i1,i2,i3,f(0),f(1),f(2)
                        ! '
                    end if
                    end if
                    end do
                    end do
                    end do
                  end if
                end if
                end do ! end side
                end do ! end axis
              else 
                ! un-known nd and orderOfAccuracy
                stop 6663
              end if
             else
               ! unknown nd 
               stop 8826 
             end if 
             return
             end
