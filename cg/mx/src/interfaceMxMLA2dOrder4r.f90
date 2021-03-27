! This file automatically generated from interfaceOpt.bf90 with bpp.
 subroutine interfaceMxMLA2dOrder4r( nd, nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,gridIndexRange1, u1,u1n,u1m, wk1, mask1,rsxy1, xy1, p1,p1n,p1m, q1,q1n,q1m, boundaryCondition1, md1a,md1b,md2a,md2b,md3a,md3b,gridIndexRange2, u2,u2n,u2m, wk2, mask2,rsxy2, xy2, p2,p2n,p2m, q2,q2n,q2m, boundaryCondition2, ipar, rpar, aa2,aa4,aa8, ipvt2,ipvt4,ipvt8, ierr )
 ! ===================================================================================
 !  Interface boundary conditions for Maxwells Equations in 3D.
 !
 !  gridType : 0=rectangular, 1=curvilinear
 !
 !  u1: solution on the domain 1 ("left") of the interface
 !  u2: solution on the domain 2 ("right") of the interface
 !
 !  u1n,u1m: past time solutions on domain 1 
 !  u2n,u2m: past time solutions on domain 2
 !
 !  p1: polarization vectors on domain 1 (for dispersive models)
 !  p2: polarization vectors on domain 2 (for dispersive models)
 ! 
 !  wk1,wk2 : work-space for jacobi update of ghost values 
 !
 !  aa2,aa4,aa8 : real work space arrays that must be saved from call to call
 !  ipvt2,ipvt4,ipvt8: integer work space arrays that must be saved from call to call
 ! ===================================================================================
 implicit none
       integer nd, nd1a,nd1b,nd2a,nd2b,nd3a,nd3b, md1a,md1b,md2a,md2b,md3a,md3b, n1a,n1b,n2a,n2b,n3a,n3b,  m1a,m1b,m2a,m2b,m3a,m3b,  ierr
       ! ------- arrays for domain 1 -------------
       real u1(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
       real u1n(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
       real u1m(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
       real wk1(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
       ! polarization vectors
       real p1(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
       real p1n(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
       real p1m(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
       ! nonlinear variables 
       real q1(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
       real q1n(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
       real q1m(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
       integer mask1(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
       real rsxy1(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1,0:nd-1)
       real xy1(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1)
       integer gridIndexRange1(0:1,0:2),boundaryCondition1(0:1,0:2)
       ! ------- arrays for domain 2 -------------
       real u2(md1a:md1b,md2a:md2b,md3a:md3b,0:*)
       real u2n(md1a:md1b,md2a:md2b,md3a:md3b,0:*)
       real u2m(md1a:md1b,md2a:md2b,md3a:md3b,0:*)
       real wk2(md1a:md1b,md2a:md2b,md3a:md3b,0:*)
       ! polarization vectors
       real p2(md1a:md1b,md2a:md2b,md3a:md3b,0:*)
       real p2n(md1a:md1b,md2a:md2b,md3a:md3b,0:*)
       real p2m(md1a:md1b,md2a:md2b,md3a:md3b,0:*)
       ! nonlinear variables 
       real q2(md1a:md1b,md2a:md2b,md3a:md3b,0:*)
       real q2n(md1a:md1b,md2a:md2b,md3a:md3b,0:*)
       real q2m(md1a:md1b,md2a:md2b,md3a:md3b,0:*)
       integer mask2(md1a:md1b,md2a:md2b,md3a:md3b)
       real rsxy2(md1a:md1b,md2a:md2b,md3a:md3b,0:nd-1,0:nd-1)
       real xy2(md1a:md1b,md2a:md2b,md3a:md3b,0:nd-1)
       integer gridIndexRange2(0:1,0:2),boundaryCondition2(0:1,0:2)
       integer ipar(0:*),numberOfInterfaceIterationsUsed
       real rpar(0:*)
       ! work space arrays that must be saved from call to call:
       real aa2(0:1,0:1,0:1,0:*),aa4(0:3,0:3,0:1,0:*),aa8(0:7,0:7,0:1,0:*)
       integer ipvt2(0:1,0:*), ipvt4(0:3,0:*), ipvt8(0:7,0:*)
 !     --- local variables ----
       integer side1,axis1,grid1,side2,axis2,grid2,gridType,orderOfAccuracy,orderOfExtrapolation,useForcing,ex,ey,ez,hx,hy,hz,useWhereMask,debug,solveForE,solveForH,axis1p1,axis1p2,axis2p1,axis2p2,nn,n1,n2,twilightZone
       real dx1(0:2),dr1(0:2),dx2(0:2),dr2(0:2)
       real dr1a(0:2), dr2a(0:2)
       real dxxx1by2i, dxxx2by2i, dxx1by12i, dxx2by12i
       real t,ep,dt,eps1,mu1,c1,eps2,mu2,c2,epsmu1,epsmu2,eta1,eta2,eta1i,eta2i
       real cSq1,cSq2,dtSq
       real absoluteErrorTolerance,relativeErrorTolerance,omega
       integer axisp1,axisp2,i1,i2,i3,is1,is2,is3,j1,j2,j3,js1,js2,js3,ks1,ks2,ks3,is,js,it,nit,k1,k2,k3,i,j
       integer i1m,i2m,i3m,j1m,j2m,j3m,i1p,i2p,i3p,j1p,j2p,j3p
       integer ii1,ii2,ii3,numberOfIterations
       integer interfaceOption,interfaceEquationsOption,initialized,forcingOption
       integer assignInterfaceValues,assignInterfaceGhostValues,setDivergenceAtInterfaces
       integer numGhost,numGhost2,giveDiv
       ! grid1
       integer ns1a,ns1b,ns2a,ns2b,ns3a,ns3b  ! save n1a,n1b,...
       integer ne1a,ne1b,ne2a,ne2b,ne3a,ne3b  ! one extra ghost for 4th-order Stage I order2 
       integer nn1a,nn1b,nn2a,nn2b,nn3a,nn3b
       integer ng1a,ng1b,ng2a,ng2b,ng3a,ng3b  ! for ghost too    
       ! grid2
       integer ms1a,ms1b,ms2a,ms2b,ms3a,ms3b  ! save m1a,m1b,...
       integer me1a,me1b,me2a,me2b,me3a,me3b  ! one extra ghost for 4th-order Stage I order2 
       integer mm1a,mm1b,mm2a,mm2b,mm3a,mm3b
       integer mg1a,mg1b,mg2a,mg2b,mg3a,mg3b  ! for ghost too 
       integer m1,m2,mm
       integer internalGhostBC
       real rx1,ry1,rz1,rx2,ry2,rz2
       real aLap0,aLap1,bLap0,bLap1,aLapX0,aLapX1,bLapY0,bLapY1,cLapX0,cLapX1,dLapY0,dLapY1,aLapSq0,aLapSq1,bLapSq0,bLapSq1
       real a0,a1,b0,b1,cc0,cc1,d0,d1,dr0,ds0
       real aNormSq,divu,dive
       real epsRatio,an1,an2,an3,aNorm,ua,ub,uc,nDotU
       real epsx
       real tau1,tau2,clap1,clap2,u1Lap,v1Lap,w1Lap,u2Lap,v2Lap,w2Lap,an1Cartesian,an2Cartesian,an3Cartesian
       real u1LapSq,v1LapSq,u2LapSq,v2LapSq,w1LapSq,w2LapSq
       integer np1a,np1b,np2a,np2b,np3a,np3b,diff(0:2)
       real rx,ry,rxx,rxy,ryy,rxxx,rxxy,rxyy,ryyy,rxxxx,rxxyy,ryyyy
       real sx,sy,sxx,sxy,syy,sxxx,sxxy,sxyy,syyy,sxxxx,sxxyy,syyyy
 !     real rv1x(0:2),rv1y(0:2),rv1xx(0:2),rv1xy(0:2),rv1yy(0:2),rv1xxx(0:2),rv1xxy(0:2),rv1xyy(0:2),rv1yyy(0:2),!          rv1xxxx(0:2),rv1xxyy(0:2),rv1yyyy(0:2)
 !     real sv1x(0:2),sv1y(0:2),sv1xx(0:2),sv1xy(0:2),sv1yy(0:2),sv1xxx(0:2),sv1xxy(0:2),sv1xyy(0:2),sv1yyy(0:2),!          sv1xxxx(0:2),sv1xxyy(0:2),sv1yyyy(0:2)
 !     real rv2x(0:2),rv2y(0:2),rv2xx(0:2),rv2xy(0:2),rv2yy(0:2),rv2xxx(0:2),rv2xxy(0:2),rv2xyy(0:2),rv2yyy(0:2),!          rv2xxxx(0:2),rv2xxyy(0:2),rv2yyyy(0:2)
 !     real sv2x(0:2),sv2y(0:2),sv2xx(0:2),sv2xy(0:2),sv2yy(0:2),sv2xxx(0:2),sv2xxy(0:2),sv2xyy(0:2),sv2yyy(0:2),!          sv2xxxx(0:2),sv2xxyy(0:2),sv2yyyy(0:2)
       integer numberOfEquations,job
       real a2(0:1,0:1),a4(0:3,0:3),a6(0:5,0:5),a8(0:7,0:7),a12(0:11,0:11),q(0:11),qOld(0:11),f(0:11),rcond,work(0:11)
       real a2f(0:1,0:1), a4f(0:3,0:3), a8f(0:7,0:7)
       integer ipvt(0:11), ipivot2(0:1), ipivot4(0:3), ipivot8(0:7)
       real a2h(0:1,0:1), a2hf(0:1,0:1), a4h(0:3,0:3), a4hf(0:3,0:3)
       integer ipivot2h(0:1), ipivot4h(0:3)
       logical first
       real err,errOld,errRatio,ratioAve
       integer debugFile,myid,parallel
       character*20 debugFileName
       ! for new evaluation method:
       real u1x,u1y,u1z,u1xx,u1xy,u1yy,u1xz,u1yz,u1zz
       real u2x,u2y,u2z,u2xx,u2xy,u2yy,u2xz,u2yz,u2zz
       real e1x,e1y,e1z,e1xx,e1xy,e1yy,e1xz,e1yz,e1zz,e1Lap
       real e2x,e2y,e2z,e2xx,e2xy,e2yy,e2xz,e2yz,e2zz,e2Lap
       real v1x,v1y,v1z,v1xx,v1xy,v1yy,v1xz,v1yz,v1zz
       real v2x,v2y,v2z,v2xx,v2xy,v2yy,v2xz,v2yz,v2zz
       real w1x,w1y,w1z,w1xx,w1xy,w1yy,w1xz,w1yz,w1zz
       real w2x,w2y,w2z,w2xx,w2xy,w2yy,w2xz,w2yz,w2zz
       real u1nx,u1ny,u1nz,u1nxx,u1nxy,u1nyy,u1nxz,u1nyz,u1nzz,u1nLap
       real v1nx,v1ny,v1nz,v1nxx,v1nxy,v1nyy,v1nxz,v1nyz,v1nzz,v1nLap
       real u2nx,u2ny,u2nz,u2nxx,u2nxy,u2nyy,u2nxz,u2nyz,u2nzz,u2nLap
       real v2nx,v2ny,v2nz,v2nxx,v2nxy,v2nyy,v2nxz,v2nyz,v2nzz,v2nLap
       real p1x,p1y,p1z,p1xx,p1xy,p1yy,p1xz,p1yz,p1zz,p1Lap
       real p1nx,p1ny,p1nz,p1nxx,p1nxy,p1nyy,p1nxz,p1nyz,p1nzz,p1nLap
       real p2x,p2y,p2z,p2xx,p2xy,p2yy,p2xz,p2yz,p2zz,p2Lap
       real p2nx,p2ny,p2nz,p2nxx,p2nxy,p2nyy,p2nxz,p2nyz,p2nzz,p2nLap
       real u1xxx,u1xxy,u1xyy,u1yyy, u1xxz,u1xzz,u1zzz, u1yyz, u1yzz
       real u2xxx,u2xxy,u2xyy,u2yyy, u2xxz,u2xzz,u2zzz, u2yyz, u2yzz
       real v1xxx,v1xxy,v1xyy,v1yyy, v1xxz,v1xzz,v1zzz, v1yyz, v1yzz
       real v2xxx,v2xxy,v2xyy,v2yyy, v2xxz,v2xzz,v2zzz, v2yyz, v2yzz
       real w1xxx,w1xxy,w1xyy,w1yyy, w1xxz,w1xzz,w1zzz, w1yyz, w1yzz
       real w2xxx,w2xxy,w2xyy,w2yyy, w2xxz,w2xzz,w2zzz, w2yyz, w2yzz
       real u1xxxx,u1xxyy,u1yyyy, u1xxzz,u1zzzz, u1yyzz
       real u2xxxx,u2xxyy,u2yyyy, u2xxzz,u2zzzz, u2yyzz
       real v1xxxx,v1xxyy,v1yyyy, v1xxzz,v1zzzz, v1yyzz
       real v2xxxx,v2xxyy,v2yyyy, v2xxzz,v2zzzz, v2yyzz
       real w1xxxx,w1xxyy,w1yyyy, w1xxzz,w1zzzz, w1yyzz
       real w2xxxx,w2xxyy,w2yyyy, w2xxzz,w2zzzz, w2yyzz
       real rxx1(0:2,0:2,0:2), rxx2(0:2,0:2,0:2)
       real dx112(0:2),dx122(0:2),dx212(0:2),dx222(0:2),dx141(0:2),dx142(0:2),dx241(0:2),dx242(0:2)
       real dr114(0:2),dr214(0:2)
       real cem1,divE1,curlE1x,curlE1y,curlE1z,nDotCurlE1,nDotLapE1
       real cem2,divE2,curlE2x,curlE2y,curlE2z,nDotCurlE2,nDotLapE2
       real c1x,c1y,c1z
       real c2x,c2y,c2z
       ! these are for the exact solution from TZ flow: 
       real ue,ve,we, we0,we1,we2,we3
       real uex,uey,uez, vex,vey,vez, wex,wey,wez, hex,hey,hez, evv(0:11), errv(0:11), maxErr
       real uexx,ueyy,uezz, vexx,veyy,vezz, wexx,weyy,wezz
       real ueLap, veLap, weLap
       real curlEex,curlEey,curlEez,nDotCurlEe,nDotLapEe
       real uexxx,uexxy,uexyy,ueyyy
       real vexxx,vexxy,vexyy,veyyy
       real wexxx,wexxy,wexyy,weyyy
       real uexxxx,uexxyy,ueyyyy,ueLapSq
       real vexxxx,vexxyy,veyyyy,veLapSq
       real wexxxx,wexxyy,weyyyy,weLapSq
       integer iv(0:2),ksv(0:2),sidea
       ! for impedance projection 
       integer useImpedanceInterfaceProjection
       real  ex1,ey1,ez1, hz1, ex2,ey2,ez2, hz2, nDotE1, nDotE2, epsNDotEI,  nDotEI, nDotE1I, nDotE2I
       real  exI, eyI, ezI, hzI, g1,g2, nDotEe
       ! boundary conditions parameters
! define BC parameters for fortran routines
! boundary conditions
      integer dirichlet,perfectElectricalConductor,perfectMagneticConductor,planeWaveBoundaryCondition,interfaceBC,symmetryBoundaryCondition,abcEM2,abcPML,abc3,abc4,abc5,rbcNonLocal,rbcLocal,characteristic,absorbing,lastBC
      parameter( dirichlet=1,perfectElectricalConductor=2,perfectMagneticConductor=3,planeWaveBoundaryCondition=4,symmetryBoundaryCondition=5,interfaceBC=6,abcEM2=7,abcPML=8,abc3=9,abc4=10,abc5=11,rbcNonLocal=12,rbcLocal=13,characteristic=14,absorbing=15,lastBC=15 )
       integer rectangular,curvilinear
       parameter(rectangular=0,curvilinear=1)
       ! Dispersion models
       integer noDispersion,drude,gdm
       parameter( noDispersion=0, drude=1, gdm=2 )
       ! Nonlinear models
       integer noNonlinearModel,multilevelAtomic
       parameter( noNonlinearModel=0, multilevelAtomic=1 )
       ! integer nonlinearModel
       ! forcing options
      ! forcingOptions -- these should match ForcingEnum in Maxwell.h 
      integer noForcing,magneticSinusoidalPointSource,gaussianSource,twilightZoneForcing,	gaussianChargeSource, userDefinedForcingOption
	integer noBoundaryForcing,planeWaveBoundaryForcing,chirpedPlaneWaveBoundaryForcing
      parameter(noForcing                =0,magneticSinusoidalPointSource =1,gaussianSource                =2,twilightZoneForcing           =3,	   gaussianChargeSource          =4,userDefinedForcingOption      =5 )
      ! boundary forcing options when solved directly for the scattered field:
      parameter( noBoundaryForcing              =0,		 planeWaveBoundaryForcing       =1,chirpedPlaneWaveBoundaryForcing=2 )
       ! Known solution options

      integer noKnownSolution,twilightZoneKnownSolution,planeWaveKnownSolution,gaussianPlaneWaveKnownSolution,gaussianIntegralKnownSolution,planeMaterialInterfaceKnownSolution,scatteringFromADiskKnownSolution,scatteringFromADielectricDiskKnownSolution,scatteringFromASphereKnownSolution,scatteringFromADielectricSphereKnownSolution,squareEigenfunctionKnownSolution,annulusEigenfunctionKnownSolution,eigenfunctionsOfASphereKnownSolution,userDefinedKnownSolution

      parameter( noKnownSolution=0,twilightZoneKnownSolution=1,planeWaveKnownSolution=2,gaussianPlaneWaveKnownSolution=3,gaussianIntegralKnownSolution=4,planeMaterialInterfaceKnownSolution=5,scatteringFromADiskKnownSolution=6,scatteringFromADielectricDiskKnownSolution=7,scatteringFromASphereKnownSolution=8,scatteringFromADielectricSphereKnownSolution=9,squareEigenfunctionKnownSolution=10,annulusEigenfunctionKnownSolution=11,eigenfunctionsOfASphereKnownSolution=12,userDefinedKnownSolution=13 )   
       integer knownSolutionOption
       integer useJacobiUpdate, useUnifiedInterfaceMacros
       integer dispersionModel1, dispersionModel2, dispersive, jv, pxc
       integer gdmParOption
       integer maxNumberOfParameters,maxNumberOfPolarizationVectors,maxPolarizationComponents
       parameter( maxNumberOfParameters=4, maxNumberOfPolarizationVectors=20, maxPolarizationComponents=20*3 )
       integer numberOfPolarizationVectors1
       real gdmPar1(0:maxNumberOfParameters-1,0:maxNumberOfPolarizationVectors-1)
       real a0v1,a1v1,b0v1,b1v1
       integer numberOfPolarizationVectors2
       real gdmPar2(0:maxNumberOfParameters-1,0:maxNumberOfPolarizationVectors-1)
       real a0v2,a1v2,b0v2,b1v2
       integer numberOfTimeDerivatives
       real evals(0:2),evalse(0:2),dpdm
       real pvals(0:maxPolarizationComponents-1),pvalse(0:maxPolarizationComponents-1)
       real pvalsm(0:maxPolarizationComponents-1),pvalsp(0:maxPolarizationComponents-1)
       real pvc(0:maxNumberOfPolarizationVectors-1),pvm(0:maxNumberOfPolarizationVectors-1),fpv(0:maxNumberOfPolarizationVectors-1)
       real rhspv(0:maxNumberOfPolarizationVectors-1)
       real betav(0:maxNumberOfPolarizationVectors-1)
       real beta,fe,evm,pSum,rhsE,rhsP,tm
       real uet,uett
       integer addForcing
       ! variables for dispersive models
       real alphaP1,alphaP2
       real fp1(0:2), fp2(0:2)
       integer ec,pc,pce
       real pv,pvt,pvtt
       real ev,evt,evtt
       real Bk, Ck
       real Csum, beta1,beta2, c2Sum
       real fev1(0:2), fpSum1(0:2), fpv1(0:2,0:maxNumberOfPolarizationVectors-1)
       real fev2(0:2), fpSum2(0:2), fpv2(0:2,0:maxNumberOfPolarizationVectors-1)
       real es(0:2), est(0:2), estt(0:2), esxx(0:2), esyy(0:2), eszz(0:2)
       real pe(0:2), pet(0:2), pett(0:2), pettSum1(0:2), pettSum2(0:2)
       real p0,p0t,p0tt
       real p1v(0:2), p2v(0:2), D1v(0:2), D2v(0:2)
       real p1ev(0:2), p2ev(0:2)
       real nDotD1,nDotD2, nDotP1, nDotP2, nDotP1e,nDotP2e, nDotDe, nDotDI
       real nDotFp1, nDotFp2
       real betac1,betac2
       ! variables for 4th-order GDM interface
       real c2PttLE,c2PttE,c2PttEm,c2PttP,c2PttPm,c2PttfE,c2PttfP
       real c4PttLE,c4PttE,c4PttEm,c4PttP,c4PttPm,c4PttfE,c4PttfP
       real c4PttLLE,c4PttLP,c4PttLEm,c4PttLPm,c4PttLfE,c4PttLfP
       real c4PttfEt,c4PttfEtt,c4PttfPt,c4PttfPtt
       real c2EtLE,c2EtE,c2EtEm,c2EtP,c2EtPm,c2EtfE,c2EtfP
       real c2PtLE,c2PtE,c2PtEm,c2PtP,c2PtPm,c2PtfE,c2PtfP
       real c2PttttLE,c2PttttE,c2PttttEm,c2PttttP,c2PttttPm,c2PttttLLE,c2PttttLEm,c2PttttLP,c2PttttLPm
       real c2PttttfE,c2PttttfP,c2PttttfEt,c2PttttfPt,c2PttttfEtt,c2PttttfPtt,c2PttttLfe,c2PttttLfP
       real c2PttEsum1, c2PttEsum2, c2PttLEsum1, c2PttLEsum2
       real c4PttLEsum1, c4PttLLEsum1, c4PttLEsum2, c4PttLLEsum2
       real c2PttttLEsum1, c2PttttLLEsum1, c2PttttLEsum2, c2PttttLLEsum2
       real coeffLap1,coeffLapSq1,coeffLap2,coeffLapSq2
       integer ismooth,nsmooth
       real eqnCoeff,eqnCoeffb 
       real curl1um1,curl1vm1,curl1um2,curl1vm2
       real curl2um1,curl2vm1,curl2um2,curl2vm2
       real evn,pvn,evec(0:2),pvec(0:2,0:maxNumberOfPolarizationVectors-1)
       real ptv(0:2,0:maxNumberOfPolarizationVectors-1),pttv(0:2,0:maxNumberOfPolarizationVectors-1),ptttv(0:2,0:maxNumberOfPolarizationVectors-1),pttttv(0:2,0:maxNumberOfPolarizationVectors-1)
       real alpha,d4,LPn,LP,LPm
       real LE(0:2)
       real LE1(0:2),LE1m(0:2),LLE1(0:2), LEx1(0:2),LEy1(0:2)
       real LE2(0:2),LE2m(0:2),LLE2(0:2), LEx2(0:2),LEy2(0:2) 
       real LfE1(0:2),fEt1(0:2),fEtt1(0:2)
       real LfP1(0:2,0:maxNumberOfPolarizationVectors-1),fPt1(0:2,0:maxNumberOfPolarizationVectors-1)
       real fPtt1(0:2,0:maxNumberOfPolarizationVectors-1)
       real LfE2(0:2),fEt2(0:2),fEtt2(0:2)
       real LfP2(0:2,0:maxNumberOfPolarizationVectors-1),fPt2(0:2,0:maxNumberOfPolarizationVectors-1)
       real fPtt2(0:2,0:maxNumberOfPolarizationVectors-1)
       real pevtt1(0:2,0:maxNumberOfPolarizationVectors-1)
       real pevtt2(0:2,0:maxNumberOfPolarizationVectors-1)
       real pevttL1(0:2,0:maxNumberOfPolarizationVectors-1)
       real pevttL2(0:2,0:maxNumberOfPolarizationVectors-1)
       real pevttx1(0:2,0:maxNumberOfPolarizationVectors-1),pevtty1(0:2,0:maxNumberOfPolarizationVectors-1)
       real pevttx2(0:2,0:maxNumberOfPolarizationVectors-1),pevtty2(0:2,0:maxNumberOfPolarizationVectors-1)
       real pevttSum1(0:2),pevttxSum1(0:2),pevttySum1(0:2)
       real pevttSum2(0:2),pevttxSum2(0:2),pevttySum2(0:2)
       real pevttt1(0:2,0:maxNumberOfPolarizationVectors-1),pevttt2(0:2,0:maxNumberOfPolarizationVectors-1)
       real pevtttt1(0:2,0:maxNumberOfPolarizationVectors-1),pevtttt2(0:2,0:maxNumberOfPolarizationVectors-1)
       real pevttLSum1(0:2),pevttLSum2(0:2),pevttttSum1(0:2),pevttttSum2(0:2)
       real petttSum
       real nDotPevttSum1, nDotPevttSum2
       real x1,y1,z1
       real esttt(0:2),estxx(0:2),estyy(0:2),esttxx(0:2),esttyy(0:2),esxxxx(0:2),esxxyy(0:2),esyyyy(0:2)
       real esttx(0:2),estty(0:2),esxxx(0:2),esxyy(0:2),esxxy(0:2),esyyy(0:2),esx(0:2),esy(0:2),estx(0:2),esty(0:2)
       real pettt(0:2),pexx(0:2),peyy(0:2),petxx(0:2),petyy(0:2),pettxx(0:2),pettyy(0:2)
       real pettx(0:2),petty(0:2),pex(0:2),pey(0:2),petx(0:2),pety(0:2)
       real estttt(0:2),petttt(0:2)
       real esL,estL,esttL,esLL,esLx,esLy
       real peL,petL,pettL,ptta,ptttta
       real pvx,pvy,pvz, pvnx,pvny,pvnz, pttxa,pttya,pttza, Lptta
       real fPttx1(0:2), fPtty1(0:2), fPttx2(0:2), fPtty2(0:2), fLPtt1(0:2), fLPtt2(0:2), fPtttt1(0:2), fPtttt2(0:2)
       real evx1(0:2),evy1(0:2),evnx1(0:2),evny1(0:2)
       real evx2(0:2),evy2(0:2),evnx2(0:2),evny2(0:2)
       real fevx1(0:2),fevy1(0:2),fevz1(0:2)
       real fevx2(0:2),fevy2(0:2),fevz2(0:2)
       real fpvx1(0:2,0:maxNumberOfPolarizationVectors-1),fpvy1(0:2,0:maxNumberOfPolarizationVectors-1)
       real fpvx2(0:2,0:maxNumberOfPolarizationVectors-1),fpvy2(0:2,0:maxNumberOfPolarizationVectors-1)
       ! For checking coefficients by delta approach: 
       real u1s(-5:5,-5:5,-5:5,0:2), u2s(-5:5,-5:5,-5:5,0:2)
       real coeffDiff, dsBig
       integer checkCoeff, saveCoeff, coeffFile, k, avoidInterfaceIterations
       integer hw1,hw2,hw3
       real f0(0:11),delta
       ! for saving coefficients
       real am(12,500)
       integer ieqn,ja
       integer i1a,i1b,i1c,i2a,i2b,i2c,i3a,i3b,i3c
       integer j1a,j1b,j1c,j2a,j2b,j2c,j3a,j3b,j3c
       ! ----- multilevel atomic model -----
       integer useNonlinearModel
       integer nonlinearModel1,numberOfAtomicLevels1,maxPar1,numPolar1
       integer nonlinearModel2,numberOfAtomicLevels2,maxPar2,numPolar2
       integer na,nce
       real q0,q0t,q0tt,q0ttt,q0tttt,q0xx,q0yy,q0x,q0y
       real pnec1,prc1,peptc1
       real pnec2,prc2,peptc2
       integer maxPar
       parameter( maxPar=20 )
       real nlPar1(0:maxPar-1,0:maxPar-1,0:2)
       real nlPar2(0:maxPar-1,0:maxPar-1,0:2)
       real fnv(0:maxPar-1),fnv1(0:maxPar-1),fnv2(0:maxPar-1)
       real fntv(0:maxPar-1),fntv1(0:maxPar-1),fntv2(0:maxPar-1)
       real fnttv(0:maxPar-1),fnttv1(0:maxPar-1),fnttv2(0:maxPar-1)
       real fntttv(0:maxPar-1),fntttv1(0:maxPar-1),fntttv2(0:maxPar-1)
       real qv(0:maxPar-1),qvn(0:maxPar-1),qvm(0:maxPar-1)
       real qt(0:maxPar-1),qtt(0:maxPar-1),qttt(0:maxPar-1),qtttt(0:maxPar-1)
       real qex(0:maxPar-1),qey(0:maxPar-1),qeLap(0:maxPar-1)
       real q1x,q1y,q1Lap,q1xx,q1yy,q1xy
       real q2x,q2y,q2Lap,q2xx,q2yy,q2xy
       real qvx,qvy,qvLap
       real evx0,evy0,pv0,ev0,evLap
       real checkParallelFortranArrayReal, checkParallelFortranArrayInt, pdiff
       integer numParallelGhost
       integer nd4a,nd4b,n4a,n4b
       integer md4a,md4b,m4a,m4b
       character*180 label
 !     --- start statement function ----
 ! .......statement functions for GDM parameters
       a0v1(jv) = gdmPar1(0,jv)
       a1v1(jv) = gdmPar1(1,jv)
       b0v1(jv) = gdmPar1(2,jv)
       b1v1(jv) = gdmPar1(3,jv)
       a0v2(jv) = gdmPar2(0,jv)
       a1v2(jv) = gdmPar2(1,jv)
       b0v2(jv) = gdmPar2(2,jv)
       b1v2(jv) = gdmPar2(3,jv)
       ! ..... statement functions for multilevel atomic model
       ! pnec  = polarizationNECoefficients
       ! prc   = populationRelaxationCoefficients
       ! peptc = populationEPtCoefficients
       pnec1(m1,m2)  = nlPar1(m1,m2,0)
       prc1(m1,m2)   = nlPar1(m1,m2,1)
       peptc1(m1,m2) = nlPar1(m1,m2,2)
       pnec2(m1,m2)  = nlPar2(m1,m2,0)
       prc2(m1,m2)   = nlPar2(m1,m2,1)
       peptc2(m1,m2) = nlPar2(m1,m2,2)
       integer kd,m,n
 !     real rx,ry,rz,sx,sy,sz,tx,ty,tz
 !*      declareDifferenceNewOrder2(u1,rsxy1,dr1,dx1,RX)
 !*      declareDifferenceNewOrder2(u2,rsxy2,dr2,dx2,RX)
 !*      declareDifferenceNewOrder4(u1,rsxy1,dr1,dx1,RX)
 !*      declareDifferenceNewOrder4(u2,rsxy2,dr2,dx2,RX)
 !.......statement functions for jacobian
 !     rx(i1,i2,i3)=rsxy1(i1,i2,i3,0,0)
 !     ry(i1,i2,i3)=rsxy1(i1,i2,i3,0,1)
 !     rz(i1,i2,i3)=rsxy1(i1,i2,i3,0,2)
 !     sx(i1,i2,i3)=rsxy1(i1,i2,i3,1,0)
 !     sy(i1,i2,i3)=rsxy1(i1,i2,i3,1,1)
 !     sz(i1,i2,i3)=rsxy1(i1,i2,i3,1,2)
 !     tx(i1,i2,i3)=rsxy1(i1,i2,i3,2,0)
 !     ty(i1,i2,i3)=rsxy1(i1,i2,i3,2,1)
 !     tz(i1,i2,i3)=rsxy1(i1,i2,i3,2,2) 
 !     The next macro call will define the difference approximation statement functions
 !*      defineDifferenceNewOrder2Components1(u1,rsxy1,dr1,dx1,RX)
 !*      defineDifferenceNewOrder2Components1(u2,rsxy2,dr2,dx2,RX)
 !*      defineDifferenceNewOrder4Components1(u1,rsxy1,dr1,dx1,RX)
 !*      defineDifferenceNewOrder4Components1(u2,rsxy2,dr2,dx2,RX)
       real t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,t86,t87,t88,t89,t90,t91,t92,t93,t94,t95,t96,t97,t98,t99,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t130,t131,t132,t133,t134,t135,t136,t137,t138,t139,t140,t141,t142,t143,t144,t145,t146,t147,t148,t149,t150,t151,t152,t153,t154,t155,t156,t157,t158,t159,t160,t161,t162,t163,t164,t165,t166,t167,t168,t169,t170,t171,t172,t173,t174,t175,t176,t177,t178,t179,t180,t181,t182,t183,t184,t185,t186,t187,t188,t189,t190,t191,t192,t193,t194,t195,t196,t197,t198,t199,t200,t201,t202,t203,t204,t205,t206,t207,t208,t209,t210,t211,t212,t213,t214,t215,t216,t217,t218,t219,t220,t221,t222,t223,t224,t225,t226,t227,t228,t229,t230,t231,t232,t233,t234,t235,t236,t237,t238,t239,t240,t241,t242,t243,t244,t245,t246,t247,t248,t249,t250,t251,t252,t253,t254,t255,t256,t257,t258,t259,t260,t261,t262,t263,t264,t265,t266,t267,t268,t269,t270,t271,t272,t273,t274,t275,t276,t277,t278,t279,t280,t281,t282,t283,t284,t285,t286,t287,t288,t289,t290,t291,t292,t293,t294,t295,t296,t297,t298,t299,t300,t301,t302,t303,t304,t305,t306,t307,t308,t309,t310,t311,t312,t313,t314,t315,t316,t317,t318,t319,t320,t321,t322,t323,t324,t325,t326,t327,t328,t329,t330,t331,t332,t333,t334,t335,t336,t337,t338,t339,t340,t341,t342,t343,t344,t345,t346,t347,t348,t349,t350,t351,t352,t353,t354,t355,t356,t357,t358,t359,t360,t361,t362,t363,t364,t365,t366,t367,t368,t369,t370,t371,t372,t373,t374,t375,t376,t377,t378,t379,t380,t381,t382,t383,t384,t385,t386,t387,t388,t389,t390,t391,t392,t393,t394,t395,t396,t397,t398,t399,t400,t401,t402,t403,t404,t405,t406,t407,t408,t409,t410,t411,t412,t413,t414,t415,t416,t417,t418,t419,t420,t421,t422,t423,t424,t425,t426,t427,t428,t429,t430,t431,t432,t433,t434,t435,t436,t437,t438,t439,t440,t441,t442,t443,t444,t445,t446,t447,t448,t449,t450,t451,t452,t453,t454,t455,t456,t457,t458,t459,t460,t461,t462,t463,t464,t465,t466,t467,t468,t469,t470,t471,t472,t473,t474,t475,t476,t477,t478,t479,t480,t481,t482,t483,t484,t485,t486,t487,t488,t489,t490,t491,t492,t493,t494,t495,t496,t497,t498,t499,t500,t501,t502,t503,t504,t505,t506,t507,t508,t509,t510,t511,t512,t513,t514,t515,t516,t517,t518,t519,t520,t521,t522,t523,t524,t525,t526,t527,t528,t529,t530,t531,t532,t533,t534,t535,t536,t537,t538,t539,t540,t541,t542,t543,t544,t545,t546,t547,t548,t549,t550,t551,t552,t553,t554,t555,t556,t557,t558,t559,t560,t561,t562,t563,t564,t565,t566,t567,t568,t569,t570,t571,t572,t573,t574,t575,t576,t577,t578,t579,t580,t581,t582,t583,t584,t585,t586,t587,t588,t589,t590,t591,t592,t593,t594,t595,t596,t597,t598,t599,t600,t601,t602,t603,t604,t605,t606,t607,t608,t609,t610,t611,t612,t613,t614,t615,t616,t617,t618,t619,t620,t621,t622,t623,t624,t625,t626,t627,t628,t629,t630,t631,t632,t633,t634,t635,t636,t637,t638,t639,t640,t641,t642,t643,t644,t645,t646,t647,t648,t649,t650,t651,t652,t653,t654,t655,t656,t657,t658,t659,t660,t661,t662,t663,t664,t665,t666,t667,t668,t669,t670,t671,t672,t673,t674,t675,t676,t677,t678,t679,t680,t681,t682,t683,t684,t685,t686,t687,t688,t689,t690,t691,t692,t693,t694,t695,t696,t697,t698,t699,t700,t701,t702,t703,t704,t705,t706,t707,t708,t709,t710,t711,t712,t713,t714,t715,t716,t717,t718,t719,t720,t721,t722,t723,t724,t725,t726,t727,t728,t729,t730,t731,t732,t733,t734,t735,t736,t737,t738,t739,t740,t741,t742,t743,t744,t745,t746,t747,t748,t749,t750,t751,t752,t753,t754,t755,t756,t757,t758,t759,t760,t761,t762,t763,t764,t765,t766,t767,t768,t769,t770,t771,t772,t773,t774,t775,t776,t777,t778,t779,t780,t781,t782,t783,t784,t785,t786,t787,t788,t789,t790,t791,t792,t793,t794,t795,t796,t797,t798,t799,t800,t801,t802,t803,t804,t805,t806,t807,t808,t809,t810,t811,t812,t813,t814,t815,t816,t817,t818,t819,t820,t821,t822,t823,t824,t825,t826,t827,t828,t829,t830,t831,t832,t833,t834,t835,t836,t837,t838,t839,t840,t841,t842,t843,t844,t845,t846,t847,t848,t849,t850,t851,t852,t853,t854,t855,t856,t857,t858,t859,t860,t861,t862,t863,t864,t865,t866,t867,t868,t869,t870,t871,t872,t873,t874,t875,t876,t877,t878,t879,t880,t881,t882,t883,t884,t885,t886,t887,t888,t889,t890,t891,t892,t893,t894,t895,t896,t897,t898,t899,t900,t901,t902,t903,t904,t905,t906,t907,t908,t909,t910,t911,t912,t913,t914,t915,t916,t917,t918,t919,t920,t921,t922,t923,t924,t925,t926,t927,t928,t929,t930,t931,t932,t933,t934,t935,t936,t937,t938,t939,t940,t941,t942,t943,t944,t945,t946,t947,t948,t949,t950,t951,t952,t953,t954,t955,t956,t957,t958,t959,t960,t961,t962,t963,t964,t965,t966,t967,t968,t969,t970,t971,t972,t973,t974,t975,t976,t977,t978,t979,t980,t981,t982,t983,t984,t985,t986,t987,t988,t989,t990,t991,t992,t993,t994,t995,t996,t997,t998,t999,t1000,t1001,t1002,t1003,t1004,t1005,t1006,t1007,t1008,t1009,t1010,t1011,t1012,t1013,t1014,t1015,t1016,t1017,t1018,t1019,t1020,t1021,t1022,t1023,t1024,t1025,t1026,t1027,t1028,t1029,t1030,t1031,t1032,t1033,t1034,t1035,t1036,t1037,t1038,t1039,t1040,t1041,t1042,t1043,t1044,t1045,t1046,t1047,t1048,t1049,t1050,t1051,t1052,t1053,t1054,t1055,t1056,t1057,t1058,t1059,t1060,t1061,t1062,t1063,t1064,t1065,t1066,t1067,t1068,t1069,t1070,t1071,t1072,t1073,t1074,t1075,t1076,t1077,t1078,t1079,t1080,t1081,t1082,t1083,t1084,t1085,t1086,t1087,t1088,t1089,t1090,t1091,t1092,t1093,t1094,t1095,t1096,t1097,t1098,t1099,t1100,t1101,t1102,t1103,t1104,t1105,t1106,t1107,t1108,t1109,t1110,t1111,t1112,t1113,t1114,t1115,t1116,t1117,t1118,t1119,t1120,t1121,t1122,t1123,t1124,t1125,t1126,t1127,t1128,t1129,t1130,t1131,t1132,t1133,t1134,t1135,t1136,t1137,t1138,t1139,t1140,t1141,t1142,t1143,t1144,t1145,t1146,t1147,t1148,t1149,t1150,t1151,t1152,t1153,t1154,t1155,t1156,t1157,t1158,t1159,t1160,t1161,t1162,t1163,t1164,t1165,t1166,t1167,t1168,t1169,t1170,t1171,t1172,t1173,t1174,t1175,t1176,t1177,t1178,t1179,t1180,t1181,t1182,t1183,t1184,t1185,t1186,t1187,t1188,t1189,t1190,t1191,t1192,t1193,t1194,t1195,t1196,t1197,t1198,t1199,t1200,t1201,t1202,t1203,t1204,t1205,t1206,t1207,t1208,t1209,t1210,t1211,t1212,t1213,t1214,t1215,t1216,t1217,t1218,t1219,t1220,t1221,t1222,t1223,t1224,t1225,t1226,t1227,t1228,t1229,t1230,t1231,t1232,t1233,t1234,t1235,t1236,t1237,t1238,t1239,t1240,t1241,t1242,t1243,t1244,t1245,t1246,t1247,t1248,t1249,t1250,t1251,t1252,t1253,t1254,t1255,t1256,t1257,t1258,t1259,t1260,t1261,t1262,t1263,t1264,t1265,t1266,t1267,t1268,t1269,t1270,t1271,t1272,t1273,t1274,t1275,t1276,t1277,t1278,t1279,t1280,t1281,t1282,t1283,t1284,t1285,t1286,t1287,t1288,t1289,t1290,t1291,t1292,t1293,t1294,t1295,t1296,t1297,t1298,t1299,t1300,t1301,t1302,t1303,t1304,t1305,t1306,t1307,t1308,t1309,t1310,t1311,t1312,t1313,t1314,t1315,t1316,t1317,t1318,t1319,t1320,t1321,t1322,t1323,t1324,t1325,t1326,t1327,t1328,t1329,t1330,t1331,t1332,t1333,t1334,t1335,t1336,t1337,t1338,t1339,t1340,t1341,t1342,t1343,t1344,t1345,t1346,t1347,t1348,t1349,t1350,t1351,t1352,t1353,t1354,t1355,t1356,t1357,t1358,t1359,t1360,t1361,t1362,t1363,t1364,t1365,t1366,t1367,t1368,t1369,t1370,t1371,t1372,t1373,t1374,t1375,t1376,t1377,t1378,t1379,t1380,t1381,t1382,t1383,t1384,t1385,t1386,t1387,t1388,t1389,t1390,t1391,t1392,t1393,t1394,t1395,t1396,t1397,t1398,t1399,t1400,t1401,t1402,t1403,t1404,t1405,t1406,t1407,t1408,t1409,t1410,t1411,t1412,t1413,t1414,t1415,t1416,t1417,t1418,t1419,t1420,t1421,t1422,t1423,t1424,t1425,t1426,t1427,t1428,t1429,t1430,t1431,t1432,t1433,t1434,t1435,t1436,t1437,t1438,t1439,t1440,t1441,t1442,t1443,t1444,t1445,t1446,t1447,t1448,t1449,t1450,t1451,t1452,t1453,t1454,t1455,t1456,t1457,t1458,t1459,t1460,t1461,t1462,t1463,t1464,t1465,t1466,t1467,t1468,t1469,t1470,t1471,t1472,t1473,t1474,t1475,t1476,t1477,t1478,t1479,t1480,t1481,t1482,t1483,t1484,t1485,t1486,t1487,t1488,t1489,t1490,t1491,t1492,t1493,t1494,t1495,t1496,t1497,t1498,t1499,t1500,t1501,t1502,t1503,t1504,t1505,t1506,t1507,t1508,t1509,t1510,t1511,t1512,t1513,t1514,t1515,t1516,t1517,t1518,t1519,t1520,t1521,t1522,t1523,t1524,t1525,t1526,t1527,t1528,t1529,t1530,t1531,t1532,t1533,t1534,t1535,t1536,t1537,t1538,t1539,t1540,t1541,t1542,t1543,t1544,t1545,t1546,t1547,t1548,t1549,t1550,t1551,t1552,t1553,t1554,t1555,t1556,t1557,t1558,t1559,t1560,t1561,t1562,t1563,t1564,t1565,t1566,t1567,t1568,t1569,t1570,t1571,t1572,t1573,t1574,t1575,t1576,t1577,t1578,t1579,t1580,t1581,t1582,t1583,t1584,t1585,t1586,t1587,t1588,t1589,t1590,t1591,t1592,t1593,t1594,t1595,t1596,t1597,t1598,t1599,t1600,t1601,t1602,t1603,t1604,t1605,t1606,t1607,t1608,t1609,t1610,t1611,t1612,t1613,t1614,t1615,t1616,t1617,t1618,t1619,t1620,t1621,t1622,t1623,t1624,t1625,t1626,t1627,t1628,t1629,t1630,t1631,t1632,t1633,t1634,t1635,t1636,t1637,t1638,t1639,t1640,t1641,t1642,t1643,t1644,t1645,t1646,t1647,t1648,t1649,t1650,t1651,t1652,t1653,t1654,t1655,t1656,t1657,t1658,t1659,t1660,t1661,t1662,t1663,t1664,t1665,t1666,t1667,t1668,t1669,t1670,t1671,t1672,t1673,t1674,t1675,t1676,t1677,t1678,t1679,t1680,t1681,t1682,t1683,t1684,t1685,t1686,t1687,t1688,t1689,t1690,t1691,t1692,t1693,t1694,t1695,t1696,t1697,t1698,t1699,t1700,t1701,t1702,t1703,t1704,t1705,t1706,t1707,t1708,t1709,t1710,t1711,t1712,t1713,t1714,t1715,t1716,t1717,t1718,t1719,t1720,t1721,t1722,t1723,t1724,t1725,t1726,t1727,t1728,t1729,t1730,t1731,t1732,t1733,t1734,t1735,t1736,t1737,t1738,t1739,t1740,t1741,t1742,t1743,t1744,t1745,t1746,t1747,t1748,t1749,t1750,t1751,t1752,t1753,t1754,t1755,t1756,t1757,t1758,t1759,t1760,t1761,t1762,t1763,t1764,t1765,t1766,t1767,t1768,t1769,t1770,t1771,t1772,t1773,t1774,t1775,t1776,t1777,t1778,t1779,t1780,t1781,t1782,t1783,t1784,t1785,t1786,t1787,t1788,t1789,t1790,t1791,t1792,t1793,t1794,t1795,t1796,t1797,t1798,t1799,t1800,t1801,t1802,t1803,t1804,t1805,t1806,t1807,t1808,t1809,t1810,t1811,t1812,t1813,t1814,t1815,t1816,t1817,t1818,t1819,t1820,t1821,t1822,t1823,t1824,t1825,t1826,t1827,t1828,t1829,t1830,t1831,t1832,t1833,t1834,t1835,t1836,t1837,t1838,t1839,t1840,t1841,t1842,t1843,t1844,t1845,t1846,t1847,t1848,t1849,t1850,t1851,t1852,t1853,t1854,t1855,t1856,t1857,t1858,t1859,t1860,t1861,t1862,t1863,t1864,t1865,t1866,t1867,t1868,t1869,t1870,t1871,t1872,t1873,t1874,t1875,t1876,t1877,t1878,t1879,t1880,t1881,t1882,t1883,t1884,t1885,t1886,t1887,t1888,t1889,t1890,t1891,t1892,t1893,t1894,t1895,t1896,t1897,t1898,t1899,t1900,t1901,t1902,t1903,t1904,t1905,t1906,t1907,t1908,t1909,t1910,t1911,t1912,t1913,t1914,t1915,t1916,t1917,t1918,t1919,t1920,t1921,t1922,t1923,t1924,t1925,t1926,t1927,t1928,t1929,t1930,t1931,t1932,t1933,t1934,t1935,t1936,t1937,t1938,t1939,t1940,t1941,t1942,t1943,t1944,t1945,t1946,t1947,t1948,t1949,t1950,t1951,t1952,t1953,t1954,t1955,t1956,t1957,t1958,t1959,t1960,t1961,t1962,t1963,t1964,t1965,t1966,t1967,t1968,t1969,t1970,t1971,t1972,t1973,t1974,t1975,t1976,t1977,t1978,t1979,t1980,t1981,t1982,t1983,t1984,t1985,t1986,t1987,t1988,t1989,t1990,t1991,t1992,t1993,t1994,t1995,t1996,t1997,t1998,t1999,t2000,t2001,t2002,t2003,t2004,t2005,t2006,t2007,t2008,t2009,t2010,t2011,t2012,t2013,t2014,t2015,t2016,t2017,t2018,t2019,t2020,t2021,t2022,t2023,t2024,t2025,t2026,t2027,t2028,t2029,t2030,t2031,t2032,t2033,t2034,t2035,t2036,t2037,t2038,t2039,t2040,t2041,t2042,t2043,t2044,t2045,t2046,t2047,t2048,t2049,t2050,t2051,t2052,t2053,t2054,t2055,t2056,t2057,t2058,t2059,t2060,t2061,t2062,t2063,t2064,t2065,t2066,t2067,t2068,t2069,t2070,t2071,t2072,t2073,t2074,t2075,t2076,t2077,t2078,t2079,t2080,t2081,t2082,t2083,t2084,t2085,t2086,t2087,t2088,t2089,t2090,t2091,t2092,t2093,t2094,t2095,t2096,t2097,t2098,t2099,t2100,t2101,t2102,t2103,t2104,t2105,t2106,t2107,t2108,t2109,t2110,t2111,t2112,t2113,t2114,t2115,t2116,t2117,t2118,t2119,t2120,t2121,t2122,t2123,t2124,t2125,t2126,t2127,t2128,t2129,t2130,t2131,t2132,t2133,t2134,t2135,t2136,t2137,t2138,t2139,t2140,t2141,t2142,t2143,t2144,t2145,t2146,t2147,t2148,t2149,t2150,t2151,t2152,t2153,t2154,t2155,t2156,t2157,t2158,t2159,t2160,t2161,t2162,t2163,t2164,t2165,t2166,t2167,t2168,t2169,t2170,t2171,t2172,t2173,t2174,t2175,t2176,t2177,t2178,t2179,t2180,t2181,t2182,t2183,t2184,t2185,t2186,t2187,t2188,t2189,t2190,t2191,t2192,t2193,t2194,t2195,t2196,t2197,t2198,t2199,t2200,t2201,t2202,t2203,t2204,t2205,t2206,t2207,t2208,t2209,t2210,t2211,t2212,t2213,t2214,t2215,t2216,t2217,t2218,t2219,t2220,t2221,t2222,t2223,t2224,t2225,t2226,t2227,t2228,t2229,t2230,t2231,t2232,t2233,t2234,t2235,t2236,t2237,t2238,t2239,t2240,t2241,t2242,t2243,t2244,t2245,t2246,t2247,t2248,t2249,t2250,t2251,t2252,t2253,t2254,t2255,t2256,t2257,t2258,t2259,t2260,t2261,t2262,t2263,t2264,t2265,t2266,t2267,t2268,t2269,t2270,t2271,t2272,t2273,t2274,t2275,t2276,t2277,t2278,t2279,t2280,t2281,t2282,t2283,t2284,t2285,t2286,t2287,t2288,t2289,t2290,t2291,t2292,t2293,t2294,t2295,t2296,t2297,t2298,t2299,t2300,t2301,t2302,t2303,t2304,t2305,t2306,t2307,t2308,t2309,t2310,t2311,t2312,t2313,t2314,t2315,t2316,t2317,t2318,t2319,t2320,t2321,t2322,t2323,t2324,t2325,t2326,t2327,t2328,t2329,t2330,t2331,t2332,t2333,t2334,t2335,t2336,t2337,t2338,t2339,t2340,t2341,t2342,t2343,t2344,t2345,t2346,t2347,t2348,t2349,t2350,t2351,t2352,t2353,t2354,t2355,t2356,t2357,t2358,t2359,t2360,t2361,t2362,t2363,t2364,t2365,t2366,t2367,t2368,t2369,t2370,t2371,t2372,t2373,t2374,t2375,t2376,t2377,t2378,t2379,t2380,t2381,t2382,t2383,t2384,t2385,t2386,t2387,t2388,t2389,t2390,t2391,t2392,t2393,t2394,t2395,t2396,t2397,t2398,t2399,t2400,t2401,t2402,t2403,t2404,t2405,t2406,t2407,t2408,t2409,t2410,t2411,t2412,t2413,t2414,t2415,t2416,t2417,t2418,t2419,t2420,t2421,t2422,t2423,t2424,t2425,t2426,t2427,t2428,t2429,t2430,t2431,t2432,t2433,t2434,t2435,t2436,t2437,t2438,t2439,t2440,t2441,t2442,t2443,t2444,t2445,t2446,t2447,t2448,t2449,t2450,t2451,t2452,t2453,t2454,t2455,t2456,t2457,t2458,t2459,t2460,t2461,t2462,t2463,t2464,t2465,t2466,t2467,t2468,t2469,t2470,t2471,t2472,t2473,t2474,t2475,t2476,t2477,t2478,t2479,t2480,t2481,t2482,t2483,t2484,t2485,t2486,t2487,t2488,t2489,t2490,t2491,t2492,t2493,t2494,t2495,t2496,t2497,t2498,t2499,t2500,t2501,t2502,t2503,t2504,t2505,t2506,t2507,t2508,t2509,t2510,t2511,t2512,t2513,t2514,t2515,t2516,t2517,t2518,t2519,t2520,t2521,t2522,t2523,t2524,t2525,t2526,t2527,t2528,t2529,t2530,t2531,t2532,t2533,t2534,t2535,t2536,t2537,t2538,t2539,t2540,t2541,t2542,t2543,t2544,t2545,t2546,t2547,t2548,t2549,t2550,t2551,t2552,t2553,t2554,t2555,t2556,t2557,t2558,t2559,t2560,t2561,t2562,t2563,t2564,t2565,t2566,t2567,t2568,t2569,t2570,t2571,t2572,t2573,t2574,t2575,t2576,t2577,t2578,t2579,t2580,t2581,t2582,t2583,t2584,t2585,t2586,t2587,t2588,t2589,t2590,t2591,t2592,t2593,t2594,t2595,t2596,t2597,t2598,t2599,t2600,t2601,t2602,t2603,t2604,t2605,t2606,t2607,t2608,t2609,t2610,t2611,t2612,t2613,t2614,t2615,t2616,t2617,t2618,t2619,t2620,t2621,t2622,t2623,t2624,t2625,t2626,t2627,t2628,t2629,t2630,t2631,t2632,t2633,t2634,t2635,t2636,t2637,t2638,t2639,t2640,t2641,t2642,t2643,t2644,t2645,t2646,t2647,t2648,t2649,t2650,t2651,t2652,t2653,t2654,t2655,t2656,t2657,t2658,t2659,t2660,t2661,t2662,t2663,t2664,t2665,t2666,t2667,t2668,t2669,t2670,t2671,t2672,t2673,t2674,t2675,t2676,t2677,t2678,t2679,t2680,t2681,t2682,t2683,t2684,t2685,t2686,t2687,t2688,t2689,t2690,t2691,t2692,t2693,t2694,t2695,t2696,t2697,t2698,t2699,t2700,t2701,t2702,t2703,t2704,t2705,t2706,t2707,t2708,t2709,t2710,t2711,t2712,t2713,t2714,t2715,t2716,t2717,t2718,t2719,t2720,t2721,t2722,t2723,t2724,t2725,t2726,t2727,t2728,t2729,t2730,t2731,t2732,t2733,t2734,t2735,t2736,t2737,t2738,t2739,t2740,t2741,t2742,t2743,t2744,t2745,t2746,t2747,t2748,t2749,t2750,t2751,t2752,t2753,t2754,t2755,t2756,t2757,t2758,t2759,t2760,t2761,t2762,t2763,t2764,t2765,t2766,t2767,t2768,t2769,t2770,t2771,t2772,t2773,t2774,t2775,t2776,t2777,t2778,t2779,t2780,t2781,t2782,t2783,t2784,t2785,t2786,t2787,t2788,t2789,t2790,t2791,t2792,t2793,t2794,t2795,t2796,t2797,t2798,t2799,t2800,t2801,t2802,t2803,t2804,t2805,t2806,t2807,t2808,t2809,t2810,t2811,t2812,t2813,t2814,t2815,t2816,t2817,t2818,t2819,t2820,t2821,t2822,t2823,t2824,t2825,t2826,t2827,t2828,t2829,t2830,t2831,t2832,t2833,t2834,t2835,t2836,t2837,t2838,t2839,t2840,t2841,t2842,t2843,t2844,t2845,t2846,t2847,t2848,t2849,t2850,t2851,t2852,t2853,t2854,t2855,t2856,t2857,t2858,t2859,t2860,t2861,t2862,t2863,t2864,t2865,t2866,t2867,t2868,t2869,t2870,t2871,t2872,t2873,t2874,t2875,t2876,t2877,t2878,t2879,t2880,t2881,t2882,t2883,t2884,t2885,t2886,t2887,t2888,t2889,t2890,t2891,t2892,t2893,t2894,t2895,t2896,t2897,t2898,t2899,t2900,t2901,t2902,t2903,t2904,t2905,t2906,t2907,t2908,t2909,t2910,t2911,t2912,t2913,t2914,t2915,t2916,t2917,t2918,t2919,t2920,t2921,t2922,t2923,t2924,t2925,t2926,t2927,t2928,t2929,t2930,t2931,t2932,t2933,t2934,t2935,t2936,t2937,t2938,t2939,t2940,t2941,t2942,t2943,t2944,t2945,t2946,t2947,t2948,t2949,t2950,t2951,t2952,t2953,t2954,t2955,t2956,t2957,t2958,t2959,t2960,t2961,t2962,t2963,t2964,t2965,t2966,t2967,t2968,t2969,t2970,t2971,t2972,t2973,t2974,t2975,t2976,t2977,t2978,t2979,t2980,t2981,t2982,t2983,t2984,t2985,t2986,t2987,t2988,t2989,t2990,t2991,t2992,t2993,t2994,t2995,t2996,t2997,t2998,t2999,t3000
       real uu1,uu1r,uu1s,uu1t,uu1rr,uu1rs,uu1ss,uu1rt,uu1st,uu1tt,uu1rrr,uu1rrs,uu1rss,uu1sss,uu1rrt,uu1rst,uu1sst,uu1rtt,uu1stt,uu1ttt,uu1rrrr,uu1rrrs,uu1rrss,uu1rsss,uu1ssss,uu1rrrt,uu1rrst,uu1rsst,uu1ssst,uu1rrtt,uu1rstt,uu1sstt,uu1rttt,uu1sttt,uu1tttt,uu1rrrrr,uu1rrrrs,uu1rrrss,uu1rrsss,uu1rssss,uu1sssss,uu1rrrrt,uu1rrrst,uu1rrsst,uu1rssst,uu1sssst,uu1rrrtt,uu1rrstt,uu1rsstt,uu1ssstt,uu1rrttt,uu1rsttt,uu1ssttt,uu1rtttt,uu1stttt,uu1ttttt,uu1rrrrrr,uu1rrrrrs,uu1rrrrss,uu1rrrsss,uu1rrssss,uu1rsssss,uu1ssssss,uu1rrrrrt,uu1rrrrst,uu1rrrsst,uu1rrssst,uu1rsssst,uu1ssssst,uu1rrrrtt,uu1rrrstt,uu1rrsstt,uu1rssstt,uu1sssstt,uu1rrrttt,uu1rrsttt,uu1rssttt,uu1sssttt,uu1rrtttt,uu1rstttt,uu1sstttt,uu1rttttt,uu1sttttt,uu1tttttt
       real uu2,uu2r,uu2s,uu2t,uu2rr,uu2rs,uu2ss,uu2rt,uu2st,uu2tt,uu2rrr,uu2rrs,uu2rss,uu2sss,uu2rrt,uu2rst,uu2sst,uu2rtt,uu2stt,uu2ttt,uu2rrrr,uu2rrrs,uu2rrss,uu2rsss,uu2ssss,uu2rrrt,uu2rrst,uu2rsst,uu2ssst,uu2rrtt,uu2rstt,uu2sstt,uu2rttt,uu2sttt,uu2tttt,uu2rrrrr,uu2rrrrs,uu2rrrss,uu2rrsss,uu2rssss,uu2sssss,uu2rrrrt,uu2rrrst,uu2rrsst,uu2rssst,uu2sssst,uu2rrrtt,uu2rrstt,uu2rsstt,uu2ssstt,uu2rrttt,uu2rsttt,uu2ssttt,uu2rtttt,uu2stttt,uu2ttttt,uu2rrrrrr,uu2rrrrrs,uu2rrrrss,uu2rrrsss,uu2rrssss,uu2rsssss,uu2ssssss,uu2rrrrrt,uu2rrrrst,uu2rrrsst,uu2rrssst,uu2rsssst,uu2ssssst,uu2rrrrtt,uu2rrrstt,uu2rrsstt,uu2rssstt,uu2sssstt,uu2rrrttt,uu2rrsttt,uu2rssttt,uu2sssttt,uu2rrtttt,uu2rstttt,uu2sstttt,uu2rttttt,uu2sttttt,uu2tttttt
       real pp1,pp1r,pp1s,pp1t,pp1rr,pp1rs,pp1ss,pp1rt,pp1st,pp1tt,pp1rrr,pp1rrs,pp1rss,pp1sss,pp1rrt,pp1rst,pp1sst,pp1rtt,pp1stt,pp1ttt,pp1rrrr,pp1rrrs,pp1rrss,pp1rsss,pp1ssss,pp1rrrt,pp1rrst,pp1rsst,pp1ssst,pp1rrtt,pp1rstt,pp1sstt,pp1rttt,pp1sttt,pp1tttt,pp1rrrrr,pp1rrrrs,pp1rrrss,pp1rrsss,pp1rssss,pp1sssss,pp1rrrrt,pp1rrrst,pp1rrsst,pp1rssst,pp1sssst,pp1rrrtt,pp1rrstt,pp1rsstt,pp1ssstt,pp1rrttt,pp1rsttt,pp1ssttt,pp1rtttt,pp1stttt,pp1ttttt,pp1rrrrrr,pp1rrrrrs,pp1rrrrss,pp1rrrsss,pp1rrssss,pp1rsssss,pp1ssssss,pp1rrrrrt,pp1rrrrst,pp1rrrsst,pp1rrssst,pp1rsssst,pp1ssssst,pp1rrrrtt,pp1rrrstt,pp1rrsstt,pp1rssstt,pp1sssstt,pp1rrrttt,pp1rrsttt,pp1rssttt,pp1sssttt,pp1rrtttt,pp1rstttt,pp1sstttt,pp1rttttt,pp1sttttt,pp1tttttt
       real pp2,pp2r,pp2s,pp2t,pp2rr,pp2rs,pp2ss,pp2rt,pp2st,pp2tt,pp2rrr,pp2rrs,pp2rss,pp2sss,pp2rrt,pp2rst,pp2sst,pp2rtt,pp2stt,pp2ttt,pp2rrrr,pp2rrrs,pp2rrss,pp2rsss,pp2ssss,pp2rrrt,pp2rrst,pp2rsst,pp2ssst,pp2rrtt,pp2rstt,pp2sstt,pp2rttt,pp2sttt,pp2tttt,pp2rrrrr,pp2rrrrs,pp2rrrss,pp2rrsss,pp2rssss,pp2sssss,pp2rrrrt,pp2rrrst,pp2rrsst,pp2rssst,pp2sssst,pp2rrrtt,pp2rrstt,pp2rsstt,pp2ssstt,pp2rrttt,pp2rsttt,pp2ssttt,pp2rtttt,pp2stttt,pp2ttttt,pp2rrrrrr,pp2rrrrrs,pp2rrrrss,pp2rrrsss,pp2rrssss,pp2rsssss,pp2ssssss,pp2rrrrrt,pp2rrrrst,pp2rrrsst,pp2rrssst,pp2rsssst,pp2ssssst,pp2rrrrtt,pp2rrrstt,pp2rrsstt,pp2rssstt,pp2sssstt,pp2rrrttt,pp2rrsttt,pp2rssttt,pp2sssttt,pp2rrtttt,pp2rstttt,pp2sstttt,pp2rttttt,pp2sttttt,pp2tttttt
       real qq1,qq1r,qq1s,qq1t,qq1rr,qq1rs,qq1ss,qq1rt,qq1st,qq1tt,qq1rrr,qq1rrs,qq1rss,qq1sss,qq1rrt,qq1rst,qq1sst,qq1rtt,qq1stt,qq1ttt,qq1rrrr,qq1rrrs,qq1rrss,qq1rsss,qq1ssss,qq1rrrt,qq1rrst,qq1rsst,qq1ssst,qq1rrtt,qq1rstt,qq1sstt,qq1rttt,qq1sttt,qq1tttt,qq1rrrrr,qq1rrrrs,qq1rrrss,qq1rrsss,qq1rssss,qq1sssss,qq1rrrrt,qq1rrrst,qq1rrsst,qq1rssst,qq1sssst,qq1rrrtt,qq1rrstt,qq1rsstt,qq1ssstt,qq1rrttt,qq1rsttt,qq1ssttt,qq1rtttt,qq1stttt,qq1ttttt,qq1rrrrrr,qq1rrrrrs,qq1rrrrss,qq1rrrsss,qq1rrssss,qq1rsssss,qq1ssssss,qq1rrrrrt,qq1rrrrst,qq1rrrsst,qq1rrssst,qq1rsssst,qq1ssssst,qq1rrrrtt,qq1rrrstt,qq1rrsstt,qq1rssstt,qq1sssstt,qq1rrrttt,qq1rrsttt,qq1rssttt,qq1sssttt,qq1rrtttt,qq1rstttt,qq1sstttt,qq1rttttt,qq1sttttt,qq1tttttt
       real qq2,qq2r,qq2s,qq2t,qq2rr,qq2rs,qq2ss,qq2rt,qq2st,qq2tt,qq2rrr,qq2rrs,qq2rss,qq2sss,qq2rrt,qq2rst,qq2sst,qq2rtt,qq2stt,qq2ttt,qq2rrrr,qq2rrrs,qq2rrss,qq2rsss,qq2ssss,qq2rrrt,qq2rrst,qq2rsst,qq2ssst,qq2rrtt,qq2rstt,qq2sstt,qq2rttt,qq2sttt,qq2tttt,qq2rrrrr,qq2rrrrs,qq2rrrss,qq2rrsss,qq2rssss,qq2sssss,qq2rrrrt,qq2rrrst,qq2rrsst,qq2rssst,qq2sssst,qq2rrrtt,qq2rrstt,qq2rsstt,qq2ssstt,qq2rrttt,qq2rsttt,qq2ssttt,qq2rtttt,qq2stttt,qq2ttttt,qq2rrrrrr,qq2rrrrrs,qq2rrrrss,qq2rrrsss,qq2rrssss,qq2rsssss,qq2ssssss,qq2rrrrrt,qq2rrrrst,qq2rrrsst,qq2rrssst,qq2rsssst,qq2ssssst,qq2rrrrtt,qq2rrrstt,qq2rrsstt,qq2rssstt,qq2sssstt,qq2rrrttt,qq2rrsttt,qq2rssttt,qq2sssttt,qq2rrtttt,qq2rstttt,qq2sstttt,qq2rttttt,qq2sttttt,qq2tttttt
       real vv1,vv1r,vv1s,vv1t,vv1rr,vv1rs,vv1ss,vv1rt,vv1st,vv1tt,vv1rrr,vv1rrs,vv1rss,vv1sss,vv1rrt,vv1rst,vv1sst,vv1rtt,vv1stt,vv1ttt,vv1rrrr,vv1rrrs,vv1rrss,vv1rsss,vv1ssss,vv1rrrt,vv1rrst,vv1rsst,vv1ssst,vv1rrtt,vv1rstt,vv1sstt,vv1rttt,vv1sttt,vv1tttt,vv1rrrrr,vv1rrrrs,vv1rrrss,vv1rrsss,vv1rssss,vv1sssss,vv1rrrrt,vv1rrrst,vv1rrsst,vv1rssst,vv1sssst,vv1rrrtt,vv1rrstt,vv1rsstt,vv1ssstt,vv1rrttt,vv1rsttt,vv1ssttt,vv1rtttt,vv1stttt,vv1ttttt,vv1rrrrrr,vv1rrrrrs,vv1rrrrss,vv1rrrsss,vv1rrssss,vv1rsssss,vv1ssssss,vv1rrrrrt,vv1rrrrst,vv1rrrsst,vv1rrssst,vv1rsssst,vv1ssssst,vv1rrrrtt,vv1rrrstt,vv1rrsstt,vv1rssstt,vv1sssstt,vv1rrrttt,vv1rrsttt,vv1rssttt,vv1sssttt,vv1rrtttt,vv1rstttt,vv1sstttt,vv1rttttt,vv1sttttt,vv1tttttt
       real vv2,vv2r,vv2s,vv2t,vv2rr,vv2rs,vv2ss,vv2rt,vv2st,vv2tt,vv2rrr,vv2rrs,vv2rss,vv2sss,vv2rrt,vv2rst,vv2sst,vv2rtt,vv2stt,vv2ttt,vv2rrrr,vv2rrrs,vv2rrss,vv2rsss,vv2ssss,vv2rrrt,vv2rrst,vv2rsst,vv2ssst,vv2rrtt,vv2rstt,vv2sstt,vv2rttt,vv2sttt,vv2tttt,vv2rrrrr,vv2rrrrs,vv2rrrss,vv2rrsss,vv2rssss,vv2sssss,vv2rrrrt,vv2rrrst,vv2rrsst,vv2rssst,vv2sssst,vv2rrrtt,vv2rrstt,vv2rsstt,vv2ssstt,vv2rrttt,vv2rsttt,vv2ssttt,vv2rtttt,vv2stttt,vv2ttttt,vv2rrrrrr,vv2rrrrrs,vv2rrrrss,vv2rrrsss,vv2rrssss,vv2rsssss,vv2ssssss,vv2rrrrrt,vv2rrrrst,vv2rrrsst,vv2rrssst,vv2rsssst,vv2ssssst,vv2rrrrtt,vv2rrrstt,vv2rrsstt,vv2rssstt,vv2sssstt,vv2rrrttt,vv2rrsttt,vv2rssttt,vv2sssttt,vv2rrtttt,vv2rstttt,vv2sstttt,vv2rttttt,vv2sttttt,vv2tttttt
       real ww1,ww1r,ww1s,ww1t,ww1rr,ww1rs,ww1ss,ww1rt,ww1st,ww1tt,ww1rrr,ww1rrs,ww1rss,ww1sss,ww1rrt,ww1rst,ww1sst,ww1rtt,ww1stt,ww1ttt,ww1rrrr,ww1rrrs,ww1rrss,ww1rsss,ww1ssss,ww1rrrt,ww1rrst,ww1rsst,ww1ssst,ww1rrtt,ww1rstt,ww1sstt,ww1rttt,ww1sttt,ww1tttt,ww1rrrrr,ww1rrrrs,ww1rrrss,ww1rrsss,ww1rssss,ww1sssss,ww1rrrrt,ww1rrrst,ww1rrsst,ww1rssst,ww1sssst,ww1rrrtt,ww1rrstt,ww1rsstt,ww1ssstt,ww1rrttt,ww1rsttt,ww1ssttt,ww1rtttt,ww1stttt,ww1ttttt,ww1rrrrrr,ww1rrrrrs,ww1rrrrss,ww1rrrsss,ww1rrssss,ww1rsssss,ww1ssssss,ww1rrrrrt,ww1rrrrst,ww1rrrsst,ww1rrssst,ww1rsssst,ww1ssssst,ww1rrrrtt,ww1rrrstt,ww1rrsstt,ww1rssstt,ww1sssstt,ww1rrrttt,ww1rrsttt,ww1rssttt,ww1sssttt,ww1rrtttt,ww1rstttt,ww1sstttt,ww1rttttt,ww1sttttt,ww1tttttt
       real ww2,ww2r,ww2s,ww2t,ww2rr,ww2rs,ww2ss,ww2rt,ww2st,ww2tt,ww2rrr,ww2rrs,ww2rss,ww2sss,ww2rrt,ww2rst,ww2sst,ww2rtt,ww2stt,ww2ttt,ww2rrrr,ww2rrrs,ww2rrss,ww2rsss,ww2ssss,ww2rrrt,ww2rrst,ww2rsst,ww2ssst,ww2rrtt,ww2rstt,ww2sstt,ww2rttt,ww2sttt,ww2tttt,ww2rrrrr,ww2rrrrs,ww2rrrss,ww2rrsss,ww2rssss,ww2sssss,ww2rrrrt,ww2rrrst,ww2rrsst,ww2rssst,ww2sssst,ww2rrrtt,ww2rrstt,ww2rsstt,ww2ssstt,ww2rrttt,ww2rsttt,ww2ssttt,ww2rtttt,ww2stttt,ww2ttttt,ww2rrrrrr,ww2rrrrrs,ww2rrrrss,ww2rrrsss,ww2rrssss,ww2rsssss,ww2ssssss,ww2rrrrrt,ww2rrrrst,ww2rrrsst,ww2rrssst,ww2rsssst,ww2ssssst,ww2rrrrtt,ww2rrrstt,ww2rrsstt,ww2rssstt,ww2sssstt,ww2rrrttt,ww2rrsttt,ww2rssttt,ww2sssttt,ww2rrtttt,ww2rstttt,ww2sstttt,ww2rttttt,ww2sttttt,ww2tttttt
        real aj1rx,aj1rxr,aj1rxs,aj1rxt,aj1rxrr,aj1rxrs,aj1rxss,aj1rxrt,aj1rxst,aj1rxtt,aj1rxrrr,aj1rxrrs,aj1rxrss,aj1rxsss,aj1rxrrt,aj1rxrst,aj1rxsst,aj1rxrtt,aj1rxstt,aj1rxttt,aj1rxrrrr,aj1rxrrrs,aj1rxrrss,aj1rxrsss,aj1rxssss,aj1rxrrrt,aj1rxrrst,aj1rxrsst,aj1rxssst,aj1rxrrtt,aj1rxrstt,aj1rxsstt,aj1rxrttt,aj1rxsttt,aj1rxtttt,aj1rxrrrrr,aj1rxrrrrs,aj1rxrrrss,aj1rxrrsss,aj1rxrssss,aj1rxsssss,aj1rxrrrrt,aj1rxrrrst,aj1rxrrsst,aj1rxrssst,aj1rxsssst,aj1rxrrrtt,aj1rxrrstt,aj1rxrsstt,aj1rxssstt,aj1rxrrttt,aj1rxrsttt,aj1rxssttt,aj1rxrtttt,aj1rxstttt,aj1rxttttt,aj1rxrrrrrr,aj1rxrrrrrs,aj1rxrrrrss,aj1rxrrrsss,aj1rxrrssss,aj1rxrsssss,aj1rxssssss,aj1rxrrrrrt,aj1rxrrrrst,aj1rxrrrsst,aj1rxrrssst,aj1rxrsssst,aj1rxssssst,aj1rxrrrrtt,aj1rxrrrstt,aj1rxrrsstt,aj1rxrssstt,aj1rxsssstt,aj1rxrrrttt,aj1rxrrsttt,aj1rxrssttt,aj1rxsssttt,aj1rxrrtttt,aj1rxrstttt,aj1rxsstttt,aj1rxrttttt,aj1rxsttttt,aj1rxtttttt
        real aj1sx,aj1sxr,aj1sxs,aj1sxt,aj1sxrr,aj1sxrs,aj1sxss,aj1sxrt,aj1sxst,aj1sxtt,aj1sxrrr,aj1sxrrs,aj1sxrss,aj1sxsss,aj1sxrrt,aj1sxrst,aj1sxsst,aj1sxrtt,aj1sxstt,aj1sxttt,aj1sxrrrr,aj1sxrrrs,aj1sxrrss,aj1sxrsss,aj1sxssss,aj1sxrrrt,aj1sxrrst,aj1sxrsst,aj1sxssst,aj1sxrrtt,aj1sxrstt,aj1sxsstt,aj1sxrttt,aj1sxsttt,aj1sxtttt,aj1sxrrrrr,aj1sxrrrrs,aj1sxrrrss,aj1sxrrsss,aj1sxrssss,aj1sxsssss,aj1sxrrrrt,aj1sxrrrst,aj1sxrrsst,aj1sxrssst,aj1sxsssst,aj1sxrrrtt,aj1sxrrstt,aj1sxrsstt,aj1sxssstt,aj1sxrrttt,aj1sxrsttt,aj1sxssttt,aj1sxrtttt,aj1sxstttt,aj1sxttttt,aj1sxrrrrrr,aj1sxrrrrrs,aj1sxrrrrss,aj1sxrrrsss,aj1sxrrssss,aj1sxrsssss,aj1sxssssss,aj1sxrrrrrt,aj1sxrrrrst,aj1sxrrrsst,aj1sxrrssst,aj1sxrsssst,aj1sxssssst,aj1sxrrrrtt,aj1sxrrrstt,aj1sxrrsstt,aj1sxrssstt,aj1sxsssstt,aj1sxrrrttt,aj1sxrrsttt,aj1sxrssttt,aj1sxsssttt,aj1sxrrtttt,aj1sxrstttt,aj1sxsstttt,aj1sxrttttt,aj1sxsttttt,aj1sxtttttt
        real aj1ry,aj1ryr,aj1rys,aj1ryt,aj1ryrr,aj1ryrs,aj1ryss,aj1ryrt,aj1ryst,aj1rytt,aj1ryrrr,aj1ryrrs,aj1ryrss,aj1rysss,aj1ryrrt,aj1ryrst,aj1rysst,aj1ryrtt,aj1rystt,aj1ryttt,aj1ryrrrr,aj1ryrrrs,aj1ryrrss,aj1ryrsss,aj1ryssss,aj1ryrrrt,aj1ryrrst,aj1ryrsst,aj1ryssst,aj1ryrrtt,aj1ryrstt,aj1rysstt,aj1ryrttt,aj1rysttt,aj1rytttt,aj1ryrrrrr,aj1ryrrrrs,aj1ryrrrss,aj1ryrrsss,aj1ryrssss,aj1rysssss,aj1ryrrrrt,aj1ryrrrst,aj1ryrrsst,aj1ryrssst,aj1rysssst,aj1ryrrrtt,aj1ryrrstt,aj1ryrsstt,aj1ryssstt,aj1ryrrttt,aj1ryrsttt,aj1ryssttt,aj1ryrtttt,aj1rystttt,aj1ryttttt,aj1ryrrrrrr,aj1ryrrrrrs,aj1ryrrrrss,aj1ryrrrsss,aj1ryrrssss,aj1ryrsssss,aj1ryssssss,aj1ryrrrrrt,aj1ryrrrrst,aj1ryrrrsst,aj1ryrrssst,aj1ryrsssst,aj1ryssssst,aj1ryrrrrtt,aj1ryrrrstt,aj1ryrrsstt,aj1ryrssstt,aj1rysssstt,aj1ryrrrttt,aj1ryrrsttt,aj1ryrssttt,aj1rysssttt,aj1ryrrtttt,aj1ryrstttt,aj1rysstttt,aj1ryrttttt,aj1rysttttt,aj1rytttttt
        real aj1sy,aj1syr,aj1sys,aj1syt,aj1syrr,aj1syrs,aj1syss,aj1syrt,aj1syst,aj1sytt,aj1syrrr,aj1syrrs,aj1syrss,aj1sysss,aj1syrrt,aj1syrst,aj1sysst,aj1syrtt,aj1systt,aj1syttt,aj1syrrrr,aj1syrrrs,aj1syrrss,aj1syrsss,aj1syssss,aj1syrrrt,aj1syrrst,aj1syrsst,aj1syssst,aj1syrrtt,aj1syrstt,aj1sysstt,aj1syrttt,aj1systtt,aj1sytttt,aj1syrrrrr,aj1syrrrrs,aj1syrrrss,aj1syrrsss,aj1syrssss,aj1sysssss,aj1syrrrrt,aj1syrrrst,aj1syrrsst,aj1syrssst,aj1sysssst,aj1syrrrtt,aj1syrrstt,aj1syrsstt,aj1syssstt,aj1syrrttt,aj1syrsttt,aj1syssttt,aj1syrtttt,aj1systttt,aj1syttttt,aj1syrrrrrr,aj1syrrrrrs,aj1syrrrrss,aj1syrrrsss,aj1syrrssss,aj1syrsssss,aj1syssssss,aj1syrrrrrt,aj1syrrrrst,aj1syrrrsst,aj1syrrssst,aj1syrsssst,aj1syssssst,aj1syrrrrtt,aj1syrrrstt,aj1syrrsstt,aj1syrssstt,aj1sysssstt,aj1syrrrttt,aj1syrrsttt,aj1syrssttt,aj1sysssttt,aj1syrrtttt,aj1syrstttt,aj1sysstttt,aj1syrttttt,aj1systtttt,aj1sytttttt
        real aj1rxx,aj1rxy,aj1rxz,aj1rxxx,aj1rxxy,aj1rxyy,aj1rxxz,aj1rxyz,aj1rxzz,aj1rxxxx,aj1rxxxy,aj1rxxyy,aj1rxyyy,aj1rxxxz,aj1rxxyz,aj1rxyyz,aj1rxxzz,aj1rxyzz,aj1rxzzz,aj1rxxxxx,aj1rxxxxy,aj1rxxxyy,aj1rxxyyy,aj1rxyyyy,aj1rxxxxz,aj1rxxxyz,aj1rxxyyz,aj1rxyyyz,aj1rxxxzz,aj1rxxyzz,aj1rxyyzz,aj1rxxzzz,aj1rxyzzz,aj1rxzzzz,aj1rxxxxxx,aj1rxxxxxy,aj1rxxxxyy,aj1rxxxyyy,aj1rxxyyyy,aj1rxyyyyy,aj1rxxxxxz,aj1rxxxxyz,aj1rxxxyyz,aj1rxxyyyz,aj1rxyyyyz,aj1rxxxxzz,aj1rxxxyzz,aj1rxxyyzz,aj1rxyyyzz,aj1rxxxzzz,aj1rxxyzzz,aj1rxyyzzz,aj1rxxzzzz,aj1rxyzzzz,aj1rxzzzzz,aj1rxxxxxxx,aj1rxxxxxxy,aj1rxxxxxyy,aj1rxxxxyyy,aj1rxxxyyyy,aj1rxxyyyyy,aj1rxyyyyyy,aj1rxxxxxxz,aj1rxxxxxyz,aj1rxxxxyyz,aj1rxxxyyyz,aj1rxxyyyyz,aj1rxyyyyyz,aj1rxxxxxzz,aj1rxxxxyzz,aj1rxxxyyzz,aj1rxxyyyzz,aj1rxyyyyzz,aj1rxxxxzzz,aj1rxxxyzzz,aj1rxxyyzzz,aj1rxyyyzzz,aj1rxxxzzzz,aj1rxxyzzzz,aj1rxyyzzzz,aj1rxxzzzzz,aj1rxyzzzzz,aj1rxzzzzzz
        real aj1sxx,aj1sxy,aj1sxz,aj1sxxx,aj1sxxy,aj1sxyy,aj1sxxz,aj1sxyz,aj1sxzz,aj1sxxxx,aj1sxxxy,aj1sxxyy,aj1sxyyy,aj1sxxxz,aj1sxxyz,aj1sxyyz,aj1sxxzz,aj1sxyzz,aj1sxzzz,aj1sxxxxx,aj1sxxxxy,aj1sxxxyy,aj1sxxyyy,aj1sxyyyy,aj1sxxxxz,aj1sxxxyz,aj1sxxyyz,aj1sxyyyz,aj1sxxxzz,aj1sxxyzz,aj1sxyyzz,aj1sxxzzz,aj1sxyzzz,aj1sxzzzz,aj1sxxxxxx,aj1sxxxxxy,aj1sxxxxyy,aj1sxxxyyy,aj1sxxyyyy,aj1sxyyyyy,aj1sxxxxxz,aj1sxxxxyz,aj1sxxxyyz,aj1sxxyyyz,aj1sxyyyyz,aj1sxxxxzz,aj1sxxxyzz,aj1sxxyyzz,aj1sxyyyzz,aj1sxxxzzz,aj1sxxyzzz,aj1sxyyzzz,aj1sxxzzzz,aj1sxyzzzz,aj1sxzzzzz,aj1sxxxxxxx,aj1sxxxxxxy,aj1sxxxxxyy,aj1sxxxxyyy,aj1sxxxyyyy,aj1sxxyyyyy,aj1sxyyyyyy,aj1sxxxxxxz,aj1sxxxxxyz,aj1sxxxxyyz,aj1sxxxyyyz,aj1sxxyyyyz,aj1sxyyyyyz,aj1sxxxxxzz,aj1sxxxxyzz,aj1sxxxyyzz,aj1sxxyyyzz,aj1sxyyyyzz,aj1sxxxxzzz,aj1sxxxyzzz,aj1sxxyyzzz,aj1sxyyyzzz,aj1sxxxzzzz,aj1sxxyzzzz,aj1sxyyzzzz,aj1sxxzzzzz,aj1sxyzzzzz,aj1sxzzzzzz
        real aj1ryx,aj1ryy,aj1ryz,aj1ryxx,aj1ryxy,aj1ryyy,aj1ryxz,aj1ryyz,aj1ryzz,aj1ryxxx,aj1ryxxy,aj1ryxyy,aj1ryyyy,aj1ryxxz,aj1ryxyz,aj1ryyyz,aj1ryxzz,aj1ryyzz,aj1ryzzz,aj1ryxxxx,aj1ryxxxy,aj1ryxxyy,aj1ryxyyy,aj1ryyyyy,aj1ryxxxz,aj1ryxxyz,aj1ryxyyz,aj1ryyyyz,aj1ryxxzz,aj1ryxyzz,aj1ryyyzz,aj1ryxzzz,aj1ryyzzz,aj1ryzzzz,aj1ryxxxxx,aj1ryxxxxy,aj1ryxxxyy,aj1ryxxyyy,aj1ryxyyyy,aj1ryyyyyy,aj1ryxxxxz,aj1ryxxxyz,aj1ryxxyyz,aj1ryxyyyz,aj1ryyyyyz,aj1ryxxxzz,aj1ryxxyzz,aj1ryxyyzz,aj1ryyyyzz,aj1ryxxzzz,aj1ryxyzzz,aj1ryyyzzz,aj1ryxzzzz,aj1ryyzzzz,aj1ryzzzzz,aj1ryxxxxxx,aj1ryxxxxxy,aj1ryxxxxyy,aj1ryxxxyyy,aj1ryxxyyyy,aj1ryxyyyyy,aj1ryyyyyyy,aj1ryxxxxxz,aj1ryxxxxyz,aj1ryxxxyyz,aj1ryxxyyyz,aj1ryxyyyyz,aj1ryyyyyyz,aj1ryxxxxzz,aj1ryxxxyzz,aj1ryxxyyzz,aj1ryxyyyzz,aj1ryyyyyzz,aj1ryxxxzzz,aj1ryxxyzzz,aj1ryxyyzzz,aj1ryyyyzzz,aj1ryxxzzzz,aj1ryxyzzzz,aj1ryyyzzzz,aj1ryxzzzzz,aj1ryyzzzzz,aj1ryzzzzzz
        real aj1syx,aj1syy,aj1syz,aj1syxx,aj1syxy,aj1syyy,aj1syxz,aj1syyz,aj1syzz,aj1syxxx,aj1syxxy,aj1syxyy,aj1syyyy,aj1syxxz,aj1syxyz,aj1syyyz,aj1syxzz,aj1syyzz,aj1syzzz,aj1syxxxx,aj1syxxxy,aj1syxxyy,aj1syxyyy,aj1syyyyy,aj1syxxxz,aj1syxxyz,aj1syxyyz,aj1syyyyz,aj1syxxzz,aj1syxyzz,aj1syyyzz,aj1syxzzz,aj1syyzzz,aj1syzzzz,aj1syxxxxx,aj1syxxxxy,aj1syxxxyy,aj1syxxyyy,aj1syxyyyy,aj1syyyyyy,aj1syxxxxz,aj1syxxxyz,aj1syxxyyz,aj1syxyyyz,aj1syyyyyz,aj1syxxxzz,aj1syxxyzz,aj1syxyyzz,aj1syyyyzz,aj1syxxzzz,aj1syxyzzz,aj1syyyzzz,aj1syxzzzz,aj1syyzzzz,aj1syzzzzz,aj1syxxxxxx,aj1syxxxxxy,aj1syxxxxyy,aj1syxxxyyy,aj1syxxyyyy,aj1syxyyyyy,aj1syyyyyyy,aj1syxxxxxz,aj1syxxxxyz,aj1syxxxyyz,aj1syxxyyyz,aj1syxyyyyz,aj1syyyyyyz,aj1syxxxxzz,aj1syxxxyzz,aj1syxxyyzz,aj1syxyyyzz,aj1syyyyyzz,aj1syxxxzzz,aj1syxxyzzz,aj1syxyyzzz,aj1syyyyzzz,aj1syxxzzzz,aj1syxyzzzz,aj1syyyzzzz,aj1syxzzzzz,aj1syyzzzzz,aj1syzzzzzz
        real aj1rz,aj1rzr,aj1rzs,aj1rzt,aj1rzrr,aj1rzrs,aj1rzss,aj1rzrt,aj1rzst,aj1rztt,aj1rzrrr,aj1rzrrs,aj1rzrss,aj1rzsss,aj1rzrrt,aj1rzrst,aj1rzsst,aj1rzrtt,aj1rzstt,aj1rzttt,aj1rzrrrr,aj1rzrrrs,aj1rzrrss,aj1rzrsss,aj1rzssss,aj1rzrrrt,aj1rzrrst,aj1rzrsst,aj1rzssst,aj1rzrrtt,aj1rzrstt,aj1rzsstt,aj1rzrttt,aj1rzsttt,aj1rztttt,aj1rzrrrrr,aj1rzrrrrs,aj1rzrrrss,aj1rzrrsss,aj1rzrssss,aj1rzsssss,aj1rzrrrrt,aj1rzrrrst,aj1rzrrsst,aj1rzrssst,aj1rzsssst,aj1rzrrrtt,aj1rzrrstt,aj1rzrsstt,aj1rzssstt,aj1rzrrttt,aj1rzrsttt,aj1rzssttt,aj1rzrtttt,aj1rzstttt,aj1rzttttt,aj1rzrrrrrr,aj1rzrrrrrs,aj1rzrrrrss,aj1rzrrrsss,aj1rzrrssss,aj1rzrsssss,aj1rzssssss,aj1rzrrrrrt,aj1rzrrrrst,aj1rzrrrsst,aj1rzrrssst,aj1rzrsssst,aj1rzssssst,aj1rzrrrrtt,aj1rzrrrstt,aj1rzrrsstt,aj1rzrssstt,aj1rzsssstt,aj1rzrrrttt,aj1rzrrsttt,aj1rzrssttt,aj1rzsssttt,aj1rzrrtttt,aj1rzrstttt,aj1rzsstttt,aj1rzrttttt,aj1rzsttttt,aj1rztttttt
        real aj1sz,aj1szr,aj1szs,aj1szt,aj1szrr,aj1szrs,aj1szss,aj1szrt,aj1szst,aj1sztt,aj1szrrr,aj1szrrs,aj1szrss,aj1szsss,aj1szrrt,aj1szrst,aj1szsst,aj1szrtt,aj1szstt,aj1szttt,aj1szrrrr,aj1szrrrs,aj1szrrss,aj1szrsss,aj1szssss,aj1szrrrt,aj1szrrst,aj1szrsst,aj1szssst,aj1szrrtt,aj1szrstt,aj1szsstt,aj1szrttt,aj1szsttt,aj1sztttt,aj1szrrrrr,aj1szrrrrs,aj1szrrrss,aj1szrrsss,aj1szrssss,aj1szsssss,aj1szrrrrt,aj1szrrrst,aj1szrrsst,aj1szrssst,aj1szsssst,aj1szrrrtt,aj1szrrstt,aj1szrsstt,aj1szssstt,aj1szrrttt,aj1szrsttt,aj1szssttt,aj1szrtttt,aj1szstttt,aj1szttttt,aj1szrrrrrr,aj1szrrrrrs,aj1szrrrrss,aj1szrrrsss,aj1szrrssss,aj1szrsssss,aj1szssssss,aj1szrrrrrt,aj1szrrrrst,aj1szrrrsst,aj1szrrssst,aj1szrsssst,aj1szssssst,aj1szrrrrtt,aj1szrrrstt,aj1szrrsstt,aj1szrssstt,aj1szsssstt,aj1szrrrttt,aj1szrrsttt,aj1szrssttt,aj1szsssttt,aj1szrrtttt,aj1szrstttt,aj1szsstttt,aj1szrttttt,aj1szsttttt,aj1sztttttt
        real aj1tx,aj1txr,aj1txs,aj1txt,aj1txrr,aj1txrs,aj1txss,aj1txrt,aj1txst,aj1txtt,aj1txrrr,aj1txrrs,aj1txrss,aj1txsss,aj1txrrt,aj1txrst,aj1txsst,aj1txrtt,aj1txstt,aj1txttt,aj1txrrrr,aj1txrrrs,aj1txrrss,aj1txrsss,aj1txssss,aj1txrrrt,aj1txrrst,aj1txrsst,aj1txssst,aj1txrrtt,aj1txrstt,aj1txsstt,aj1txrttt,aj1txsttt,aj1txtttt,aj1txrrrrr,aj1txrrrrs,aj1txrrrss,aj1txrrsss,aj1txrssss,aj1txsssss,aj1txrrrrt,aj1txrrrst,aj1txrrsst,aj1txrssst,aj1txsssst,aj1txrrrtt,aj1txrrstt,aj1txrsstt,aj1txssstt,aj1txrrttt,aj1txrsttt,aj1txssttt,aj1txrtttt,aj1txstttt,aj1txttttt,aj1txrrrrrr,aj1txrrrrrs,aj1txrrrrss,aj1txrrrsss,aj1txrrssss,aj1txrsssss,aj1txssssss,aj1txrrrrrt,aj1txrrrrst,aj1txrrrsst,aj1txrrssst,aj1txrsssst,aj1txssssst,aj1txrrrrtt,aj1txrrrstt,aj1txrrsstt,aj1txrssstt,aj1txsssstt,aj1txrrrttt,aj1txrrsttt,aj1txrssttt,aj1txsssttt,aj1txrrtttt,aj1txrstttt,aj1txsstttt,aj1txrttttt,aj1txsttttt,aj1txtttttt
        real aj1ty,aj1tyr,aj1tys,aj1tyt,aj1tyrr,aj1tyrs,aj1tyss,aj1tyrt,aj1tyst,aj1tytt,aj1tyrrr,aj1tyrrs,aj1tyrss,aj1tysss,aj1tyrrt,aj1tyrst,aj1tysst,aj1tyrtt,aj1tystt,aj1tyttt,aj1tyrrrr,aj1tyrrrs,aj1tyrrss,aj1tyrsss,aj1tyssss,aj1tyrrrt,aj1tyrrst,aj1tyrsst,aj1tyssst,aj1tyrrtt,aj1tyrstt,aj1tysstt,aj1tyrttt,aj1tysttt,aj1tytttt,aj1tyrrrrr,aj1tyrrrrs,aj1tyrrrss,aj1tyrrsss,aj1tyrssss,aj1tysssss,aj1tyrrrrt,aj1tyrrrst,aj1tyrrsst,aj1tyrssst,aj1tysssst,aj1tyrrrtt,aj1tyrrstt,aj1tyrsstt,aj1tyssstt,aj1tyrrttt,aj1tyrsttt,aj1tyssttt,aj1tyrtttt,aj1tystttt,aj1tyttttt,aj1tyrrrrrr,aj1tyrrrrrs,aj1tyrrrrss,aj1tyrrrsss,aj1tyrrssss,aj1tyrsssss,aj1tyssssss,aj1tyrrrrrt,aj1tyrrrrst,aj1tyrrrsst,aj1tyrrssst,aj1tyrsssst,aj1tyssssst,aj1tyrrrrtt,aj1tyrrrstt,aj1tyrrsstt,aj1tyrssstt,aj1tysssstt,aj1tyrrrttt,aj1tyrrsttt,aj1tyrssttt,aj1tysssttt,aj1tyrrtttt,aj1tyrstttt,aj1tysstttt,aj1tyrttttt,aj1tysttttt,aj1tytttttt
        real aj1tz,aj1tzr,aj1tzs,aj1tzt,aj1tzrr,aj1tzrs,aj1tzss,aj1tzrt,aj1tzst,aj1tztt,aj1tzrrr,aj1tzrrs,aj1tzrss,aj1tzsss,aj1tzrrt,aj1tzrst,aj1tzsst,aj1tzrtt,aj1tzstt,aj1tzttt,aj1tzrrrr,aj1tzrrrs,aj1tzrrss,aj1tzrsss,aj1tzssss,aj1tzrrrt,aj1tzrrst,aj1tzrsst,aj1tzssst,aj1tzrrtt,aj1tzrstt,aj1tzsstt,aj1tzrttt,aj1tzsttt,aj1tztttt,aj1tzrrrrr,aj1tzrrrrs,aj1tzrrrss,aj1tzrrsss,aj1tzrssss,aj1tzsssss,aj1tzrrrrt,aj1tzrrrst,aj1tzrrsst,aj1tzrssst,aj1tzsssst,aj1tzrrrtt,aj1tzrrstt,aj1tzrsstt,aj1tzssstt,aj1tzrrttt,aj1tzrsttt,aj1tzssttt,aj1tzrtttt,aj1tzstttt,aj1tzttttt,aj1tzrrrrrr,aj1tzrrrrrs,aj1tzrrrrss,aj1tzrrrsss,aj1tzrrssss,aj1tzrsssss,aj1tzssssss,aj1tzrrrrrt,aj1tzrrrrst,aj1tzrrrsst,aj1tzrrssst,aj1tzrsssst,aj1tzssssst,aj1tzrrrrtt,aj1tzrrrstt,aj1tzrrsstt,aj1tzrssstt,aj1tzsssstt,aj1tzrrrttt,aj1tzrrsttt,aj1tzrssttt,aj1tzsssttt,aj1tzrrtttt,aj1tzrstttt,aj1tzsstttt,aj1tzrttttt,aj1tzsttttt,aj1tztttttt
        real aj1rzx,aj1rzy,aj1rzz,aj1rzxx,aj1rzxy,aj1rzyy,aj1rzxz,aj1rzyz,aj1rzzz,aj1rzxxx,aj1rzxxy,aj1rzxyy,aj1rzyyy,aj1rzxxz,aj1rzxyz,aj1rzyyz,aj1rzxzz,aj1rzyzz,aj1rzzzz,aj1rzxxxx,aj1rzxxxy,aj1rzxxyy,aj1rzxyyy,aj1rzyyyy,aj1rzxxxz,aj1rzxxyz,aj1rzxyyz,aj1rzyyyz,aj1rzxxzz,aj1rzxyzz,aj1rzyyzz,aj1rzxzzz,aj1rzyzzz,aj1rzzzzz,aj1rzxxxxx,aj1rzxxxxy,aj1rzxxxyy,aj1rzxxyyy,aj1rzxyyyy,aj1rzyyyyy,aj1rzxxxxz,aj1rzxxxyz,aj1rzxxyyz,aj1rzxyyyz,aj1rzyyyyz,aj1rzxxxzz,aj1rzxxyzz,aj1rzxyyzz,aj1rzyyyzz,aj1rzxxzzz,aj1rzxyzzz,aj1rzyyzzz,aj1rzxzzzz,aj1rzyzzzz,aj1rzzzzzz,aj1rzxxxxxx,aj1rzxxxxxy,aj1rzxxxxyy,aj1rzxxxyyy,aj1rzxxyyyy,aj1rzxyyyyy,aj1rzyyyyyy,aj1rzxxxxxz,aj1rzxxxxyz,aj1rzxxxyyz,aj1rzxxyyyz,aj1rzxyyyyz,aj1rzyyyyyz,aj1rzxxxxzz,aj1rzxxxyzz,aj1rzxxyyzz,aj1rzxyyyzz,aj1rzyyyyzz,aj1rzxxxzzz,aj1rzxxyzzz,aj1rzxyyzzz,aj1rzyyyzzz,aj1rzxxzzzz,aj1rzxyzzzz,aj1rzyyzzzz,aj1rzxzzzzz,aj1rzyzzzzz,aj1rzzzzzzz
        real aj1szx,aj1szy,aj1szz,aj1szxx,aj1szxy,aj1szyy,aj1szxz,aj1szyz,aj1szzz,aj1szxxx,aj1szxxy,aj1szxyy,aj1szyyy,aj1szxxz,aj1szxyz,aj1szyyz,aj1szxzz,aj1szyzz,aj1szzzz,aj1szxxxx,aj1szxxxy,aj1szxxyy,aj1szxyyy,aj1szyyyy,aj1szxxxz,aj1szxxyz,aj1szxyyz,aj1szyyyz,aj1szxxzz,aj1szxyzz,aj1szyyzz,aj1szxzzz,aj1szyzzz,aj1szzzzz,aj1szxxxxx,aj1szxxxxy,aj1szxxxyy,aj1szxxyyy,aj1szxyyyy,aj1szyyyyy,aj1szxxxxz,aj1szxxxyz,aj1szxxyyz,aj1szxyyyz,aj1szyyyyz,aj1szxxxzz,aj1szxxyzz,aj1szxyyzz,aj1szyyyzz,aj1szxxzzz,aj1szxyzzz,aj1szyyzzz,aj1szxzzzz,aj1szyzzzz,aj1szzzzzz,aj1szxxxxxx,aj1szxxxxxy,aj1szxxxxyy,aj1szxxxyyy,aj1szxxyyyy,aj1szxyyyyy,aj1szyyyyyy,aj1szxxxxxz,aj1szxxxxyz,aj1szxxxyyz,aj1szxxyyyz,aj1szxyyyyz,aj1szyyyyyz,aj1szxxxxzz,aj1szxxxyzz,aj1szxxyyzz,aj1szxyyyzz,aj1szyyyyzz,aj1szxxxzzz,aj1szxxyzzz,aj1szxyyzzz,aj1szyyyzzz,aj1szxxzzzz,aj1szxyzzzz,aj1szyyzzzz,aj1szxzzzzz,aj1szyzzzzz,aj1szzzzzzz
        real aj1txx,aj1txy,aj1txz,aj1txxx,aj1txxy,aj1txyy,aj1txxz,aj1txyz,aj1txzz,aj1txxxx,aj1txxxy,aj1txxyy,aj1txyyy,aj1txxxz,aj1txxyz,aj1txyyz,aj1txxzz,aj1txyzz,aj1txzzz,aj1txxxxx,aj1txxxxy,aj1txxxyy,aj1txxyyy,aj1txyyyy,aj1txxxxz,aj1txxxyz,aj1txxyyz,aj1txyyyz,aj1txxxzz,aj1txxyzz,aj1txyyzz,aj1txxzzz,aj1txyzzz,aj1txzzzz,aj1txxxxxx,aj1txxxxxy,aj1txxxxyy,aj1txxxyyy,aj1txxyyyy,aj1txyyyyy,aj1txxxxxz,aj1txxxxyz,aj1txxxyyz,aj1txxyyyz,aj1txyyyyz,aj1txxxxzz,aj1txxxyzz,aj1txxyyzz,aj1txyyyzz,aj1txxxzzz,aj1txxyzzz,aj1txyyzzz,aj1txxzzzz,aj1txyzzzz,aj1txzzzzz,aj1txxxxxxx,aj1txxxxxxy,aj1txxxxxyy,aj1txxxxyyy,aj1txxxyyyy,aj1txxyyyyy,aj1txyyyyyy,aj1txxxxxxz,aj1txxxxxyz,aj1txxxxyyz,aj1txxxyyyz,aj1txxyyyyz,aj1txyyyyyz,aj1txxxxxzz,aj1txxxxyzz,aj1txxxyyzz,aj1txxyyyzz,aj1txyyyyzz,aj1txxxxzzz,aj1txxxyzzz,aj1txxyyzzz,aj1txyyyzzz,aj1txxxzzzz,aj1txxyzzzz,aj1txyyzzzz,aj1txxzzzzz,aj1txyzzzzz,aj1txzzzzzz
        real aj1tyx,aj1tyy,aj1tyz,aj1tyxx,aj1tyxy,aj1tyyy,aj1tyxz,aj1tyyz,aj1tyzz,aj1tyxxx,aj1tyxxy,aj1tyxyy,aj1tyyyy,aj1tyxxz,aj1tyxyz,aj1tyyyz,aj1tyxzz,aj1tyyzz,aj1tyzzz,aj1tyxxxx,aj1tyxxxy,aj1tyxxyy,aj1tyxyyy,aj1tyyyyy,aj1tyxxxz,aj1tyxxyz,aj1tyxyyz,aj1tyyyyz,aj1tyxxzz,aj1tyxyzz,aj1tyyyzz,aj1tyxzzz,aj1tyyzzz,aj1tyzzzz,aj1tyxxxxx,aj1tyxxxxy,aj1tyxxxyy,aj1tyxxyyy,aj1tyxyyyy,aj1tyyyyyy,aj1tyxxxxz,aj1tyxxxyz,aj1tyxxyyz,aj1tyxyyyz,aj1tyyyyyz,aj1tyxxxzz,aj1tyxxyzz,aj1tyxyyzz,aj1tyyyyzz,aj1tyxxzzz,aj1tyxyzzz,aj1tyyyzzz,aj1tyxzzzz,aj1tyyzzzz,aj1tyzzzzz,aj1tyxxxxxx,aj1tyxxxxxy,aj1tyxxxxyy,aj1tyxxxyyy,aj1tyxxyyyy,aj1tyxyyyyy,aj1tyyyyyyy,aj1tyxxxxxz,aj1tyxxxxyz,aj1tyxxxyyz,aj1tyxxyyyz,aj1tyxyyyyz,aj1tyyyyyyz,aj1tyxxxxzz,aj1tyxxxyzz,aj1tyxxyyzz,aj1tyxyyyzz,aj1tyyyyyzz,aj1tyxxxzzz,aj1tyxxyzzz,aj1tyxyyzzz,aj1tyyyyzzz,aj1tyxxzzzz,aj1tyxyzzzz,aj1tyyyzzzz,aj1tyxzzzzz,aj1tyyzzzzz,aj1tyzzzzzz
        real aj1tzx,aj1tzy,aj1tzz,aj1tzxx,aj1tzxy,aj1tzyy,aj1tzxz,aj1tzyz,aj1tzzz,aj1tzxxx,aj1tzxxy,aj1tzxyy,aj1tzyyy,aj1tzxxz,aj1tzxyz,aj1tzyyz,aj1tzxzz,aj1tzyzz,aj1tzzzz,aj1tzxxxx,aj1tzxxxy,aj1tzxxyy,aj1tzxyyy,aj1tzyyyy,aj1tzxxxz,aj1tzxxyz,aj1tzxyyz,aj1tzyyyz,aj1tzxxzz,aj1tzxyzz,aj1tzyyzz,aj1tzxzzz,aj1tzyzzz,aj1tzzzzz,aj1tzxxxxx,aj1tzxxxxy,aj1tzxxxyy,aj1tzxxyyy,aj1tzxyyyy,aj1tzyyyyy,aj1tzxxxxz,aj1tzxxxyz,aj1tzxxyyz,aj1tzxyyyz,aj1tzyyyyz,aj1tzxxxzz,aj1tzxxyzz,aj1tzxyyzz,aj1tzyyyzz,aj1tzxxzzz,aj1tzxyzzz,aj1tzyyzzz,aj1tzxzzzz,aj1tzyzzzz,aj1tzzzzzz,aj1tzxxxxxx,aj1tzxxxxxy,aj1tzxxxxyy,aj1tzxxxyyy,aj1tzxxyyyy,aj1tzxyyyyy,aj1tzyyyyyy,aj1tzxxxxxz,aj1tzxxxxyz,aj1tzxxxyyz,aj1tzxxyyyz,aj1tzxyyyyz,aj1tzyyyyyz,aj1tzxxxxzz,aj1tzxxxyzz,aj1tzxxyyzz,aj1tzxyyyzz,aj1tzyyyyzz,aj1tzxxxzzz,aj1tzxxyzzz,aj1tzxyyzzz,aj1tzyyyzzz,aj1tzxxzzzz,aj1tzxyzzzz,aj1tzyyzzzz,aj1tzxzzzzz,aj1tzyzzzzz,aj1tzzzzzzz
        real aj2rx,aj2rxr,aj2rxs,aj2rxt,aj2rxrr,aj2rxrs,aj2rxss,aj2rxrt,aj2rxst,aj2rxtt,aj2rxrrr,aj2rxrrs,aj2rxrss,aj2rxsss,aj2rxrrt,aj2rxrst,aj2rxsst,aj2rxrtt,aj2rxstt,aj2rxttt,aj2rxrrrr,aj2rxrrrs,aj2rxrrss,aj2rxrsss,aj2rxssss,aj2rxrrrt,aj2rxrrst,aj2rxrsst,aj2rxssst,aj2rxrrtt,aj2rxrstt,aj2rxsstt,aj2rxrttt,aj2rxsttt,aj2rxtttt,aj2rxrrrrr,aj2rxrrrrs,aj2rxrrrss,aj2rxrrsss,aj2rxrssss,aj2rxsssss,aj2rxrrrrt,aj2rxrrrst,aj2rxrrsst,aj2rxrssst,aj2rxsssst,aj2rxrrrtt,aj2rxrrstt,aj2rxrsstt,aj2rxssstt,aj2rxrrttt,aj2rxrsttt,aj2rxssttt,aj2rxrtttt,aj2rxstttt,aj2rxttttt,aj2rxrrrrrr,aj2rxrrrrrs,aj2rxrrrrss,aj2rxrrrsss,aj2rxrrssss,aj2rxrsssss,aj2rxssssss,aj2rxrrrrrt,aj2rxrrrrst,aj2rxrrrsst,aj2rxrrssst,aj2rxrsssst,aj2rxssssst,aj2rxrrrrtt,aj2rxrrrstt,aj2rxrrsstt,aj2rxrssstt,aj2rxsssstt,aj2rxrrrttt,aj2rxrrsttt,aj2rxrssttt,aj2rxsssttt,aj2rxrrtttt,aj2rxrstttt,aj2rxsstttt,aj2rxrttttt,aj2rxsttttt,aj2rxtttttt
        real aj2sx,aj2sxr,aj2sxs,aj2sxt,aj2sxrr,aj2sxrs,aj2sxss,aj2sxrt,aj2sxst,aj2sxtt,aj2sxrrr,aj2sxrrs,aj2sxrss,aj2sxsss,aj2sxrrt,aj2sxrst,aj2sxsst,aj2sxrtt,aj2sxstt,aj2sxttt,aj2sxrrrr,aj2sxrrrs,aj2sxrrss,aj2sxrsss,aj2sxssss,aj2sxrrrt,aj2sxrrst,aj2sxrsst,aj2sxssst,aj2sxrrtt,aj2sxrstt,aj2sxsstt,aj2sxrttt,aj2sxsttt,aj2sxtttt,aj2sxrrrrr,aj2sxrrrrs,aj2sxrrrss,aj2sxrrsss,aj2sxrssss,aj2sxsssss,aj2sxrrrrt,aj2sxrrrst,aj2sxrrsst,aj2sxrssst,aj2sxsssst,aj2sxrrrtt,aj2sxrrstt,aj2sxrsstt,aj2sxssstt,aj2sxrrttt,aj2sxrsttt,aj2sxssttt,aj2sxrtttt,aj2sxstttt,aj2sxttttt,aj2sxrrrrrr,aj2sxrrrrrs,aj2sxrrrrss,aj2sxrrrsss,aj2sxrrssss,aj2sxrsssss,aj2sxssssss,aj2sxrrrrrt,aj2sxrrrrst,aj2sxrrrsst,aj2sxrrssst,aj2sxrsssst,aj2sxssssst,aj2sxrrrrtt,aj2sxrrrstt,aj2sxrrsstt,aj2sxrssstt,aj2sxsssstt,aj2sxrrrttt,aj2sxrrsttt,aj2sxrssttt,aj2sxsssttt,aj2sxrrtttt,aj2sxrstttt,aj2sxsstttt,aj2sxrttttt,aj2sxsttttt,aj2sxtttttt
        real aj2ry,aj2ryr,aj2rys,aj2ryt,aj2ryrr,aj2ryrs,aj2ryss,aj2ryrt,aj2ryst,aj2rytt,aj2ryrrr,aj2ryrrs,aj2ryrss,aj2rysss,aj2ryrrt,aj2ryrst,aj2rysst,aj2ryrtt,aj2rystt,aj2ryttt,aj2ryrrrr,aj2ryrrrs,aj2ryrrss,aj2ryrsss,aj2ryssss,aj2ryrrrt,aj2ryrrst,aj2ryrsst,aj2ryssst,aj2ryrrtt,aj2ryrstt,aj2rysstt,aj2ryrttt,aj2rysttt,aj2rytttt,aj2ryrrrrr,aj2ryrrrrs,aj2ryrrrss,aj2ryrrsss,aj2ryrssss,aj2rysssss,aj2ryrrrrt,aj2ryrrrst,aj2ryrrsst,aj2ryrssst,aj2rysssst,aj2ryrrrtt,aj2ryrrstt,aj2ryrsstt,aj2ryssstt,aj2ryrrttt,aj2ryrsttt,aj2ryssttt,aj2ryrtttt,aj2rystttt,aj2ryttttt,aj2ryrrrrrr,aj2ryrrrrrs,aj2ryrrrrss,aj2ryrrrsss,aj2ryrrssss,aj2ryrsssss,aj2ryssssss,aj2ryrrrrrt,aj2ryrrrrst,aj2ryrrrsst,aj2ryrrssst,aj2ryrsssst,aj2ryssssst,aj2ryrrrrtt,aj2ryrrrstt,aj2ryrrsstt,aj2ryrssstt,aj2rysssstt,aj2ryrrrttt,aj2ryrrsttt,aj2ryrssttt,aj2rysssttt,aj2ryrrtttt,aj2ryrstttt,aj2rysstttt,aj2ryrttttt,aj2rysttttt,aj2rytttttt
        real aj2sy,aj2syr,aj2sys,aj2syt,aj2syrr,aj2syrs,aj2syss,aj2syrt,aj2syst,aj2sytt,aj2syrrr,aj2syrrs,aj2syrss,aj2sysss,aj2syrrt,aj2syrst,aj2sysst,aj2syrtt,aj2systt,aj2syttt,aj2syrrrr,aj2syrrrs,aj2syrrss,aj2syrsss,aj2syssss,aj2syrrrt,aj2syrrst,aj2syrsst,aj2syssst,aj2syrrtt,aj2syrstt,aj2sysstt,aj2syrttt,aj2systtt,aj2sytttt,aj2syrrrrr,aj2syrrrrs,aj2syrrrss,aj2syrrsss,aj2syrssss,aj2sysssss,aj2syrrrrt,aj2syrrrst,aj2syrrsst,aj2syrssst,aj2sysssst,aj2syrrrtt,aj2syrrstt,aj2syrsstt,aj2syssstt,aj2syrrttt,aj2syrsttt,aj2syssttt,aj2syrtttt,aj2systttt,aj2syttttt,aj2syrrrrrr,aj2syrrrrrs,aj2syrrrrss,aj2syrrrsss,aj2syrrssss,aj2syrsssss,aj2syssssss,aj2syrrrrrt,aj2syrrrrst,aj2syrrrsst,aj2syrrssst,aj2syrsssst,aj2syssssst,aj2syrrrrtt,aj2syrrrstt,aj2syrrsstt,aj2syrssstt,aj2sysssstt,aj2syrrrttt,aj2syrrsttt,aj2syrssttt,aj2sysssttt,aj2syrrtttt,aj2syrstttt,aj2sysstttt,aj2syrttttt,aj2systtttt,aj2sytttttt
        real aj2rxx,aj2rxy,aj2rxz,aj2rxxx,aj2rxxy,aj2rxyy,aj2rxxz,aj2rxyz,aj2rxzz,aj2rxxxx,aj2rxxxy,aj2rxxyy,aj2rxyyy,aj2rxxxz,aj2rxxyz,aj2rxyyz,aj2rxxzz,aj2rxyzz,aj2rxzzz,aj2rxxxxx,aj2rxxxxy,aj2rxxxyy,aj2rxxyyy,aj2rxyyyy,aj2rxxxxz,aj2rxxxyz,aj2rxxyyz,aj2rxyyyz,aj2rxxxzz,aj2rxxyzz,aj2rxyyzz,aj2rxxzzz,aj2rxyzzz,aj2rxzzzz,aj2rxxxxxx,aj2rxxxxxy,aj2rxxxxyy,aj2rxxxyyy,aj2rxxyyyy,aj2rxyyyyy,aj2rxxxxxz,aj2rxxxxyz,aj2rxxxyyz,aj2rxxyyyz,aj2rxyyyyz,aj2rxxxxzz,aj2rxxxyzz,aj2rxxyyzz,aj2rxyyyzz,aj2rxxxzzz,aj2rxxyzzz,aj2rxyyzzz,aj2rxxzzzz,aj2rxyzzzz,aj2rxzzzzz,aj2rxxxxxxx,aj2rxxxxxxy,aj2rxxxxxyy,aj2rxxxxyyy,aj2rxxxyyyy,aj2rxxyyyyy,aj2rxyyyyyy,aj2rxxxxxxz,aj2rxxxxxyz,aj2rxxxxyyz,aj2rxxxyyyz,aj2rxxyyyyz,aj2rxyyyyyz,aj2rxxxxxzz,aj2rxxxxyzz,aj2rxxxyyzz,aj2rxxyyyzz,aj2rxyyyyzz,aj2rxxxxzzz,aj2rxxxyzzz,aj2rxxyyzzz,aj2rxyyyzzz,aj2rxxxzzzz,aj2rxxyzzzz,aj2rxyyzzzz,aj2rxxzzzzz,aj2rxyzzzzz,aj2rxzzzzzz
        real aj2sxx,aj2sxy,aj2sxz,aj2sxxx,aj2sxxy,aj2sxyy,aj2sxxz,aj2sxyz,aj2sxzz,aj2sxxxx,aj2sxxxy,aj2sxxyy,aj2sxyyy,aj2sxxxz,aj2sxxyz,aj2sxyyz,aj2sxxzz,aj2sxyzz,aj2sxzzz,aj2sxxxxx,aj2sxxxxy,aj2sxxxyy,aj2sxxyyy,aj2sxyyyy,aj2sxxxxz,aj2sxxxyz,aj2sxxyyz,aj2sxyyyz,aj2sxxxzz,aj2sxxyzz,aj2sxyyzz,aj2sxxzzz,aj2sxyzzz,aj2sxzzzz,aj2sxxxxxx,aj2sxxxxxy,aj2sxxxxyy,aj2sxxxyyy,aj2sxxyyyy,aj2sxyyyyy,aj2sxxxxxz,aj2sxxxxyz,aj2sxxxyyz,aj2sxxyyyz,aj2sxyyyyz,aj2sxxxxzz,aj2sxxxyzz,aj2sxxyyzz,aj2sxyyyzz,aj2sxxxzzz,aj2sxxyzzz,aj2sxyyzzz,aj2sxxzzzz,aj2sxyzzzz,aj2sxzzzzz,aj2sxxxxxxx,aj2sxxxxxxy,aj2sxxxxxyy,aj2sxxxxyyy,aj2sxxxyyyy,aj2sxxyyyyy,aj2sxyyyyyy,aj2sxxxxxxz,aj2sxxxxxyz,aj2sxxxxyyz,aj2sxxxyyyz,aj2sxxyyyyz,aj2sxyyyyyz,aj2sxxxxxzz,aj2sxxxxyzz,aj2sxxxyyzz,aj2sxxyyyzz,aj2sxyyyyzz,aj2sxxxxzzz,aj2sxxxyzzz,aj2sxxyyzzz,aj2sxyyyzzz,aj2sxxxzzzz,aj2sxxyzzzz,aj2sxyyzzzz,aj2sxxzzzzz,aj2sxyzzzzz,aj2sxzzzzzz
        real aj2ryx,aj2ryy,aj2ryz,aj2ryxx,aj2ryxy,aj2ryyy,aj2ryxz,aj2ryyz,aj2ryzz,aj2ryxxx,aj2ryxxy,aj2ryxyy,aj2ryyyy,aj2ryxxz,aj2ryxyz,aj2ryyyz,aj2ryxzz,aj2ryyzz,aj2ryzzz,aj2ryxxxx,aj2ryxxxy,aj2ryxxyy,aj2ryxyyy,aj2ryyyyy,aj2ryxxxz,aj2ryxxyz,aj2ryxyyz,aj2ryyyyz,aj2ryxxzz,aj2ryxyzz,aj2ryyyzz,aj2ryxzzz,aj2ryyzzz,aj2ryzzzz,aj2ryxxxxx,aj2ryxxxxy,aj2ryxxxyy,aj2ryxxyyy,aj2ryxyyyy,aj2ryyyyyy,aj2ryxxxxz,aj2ryxxxyz,aj2ryxxyyz,aj2ryxyyyz,aj2ryyyyyz,aj2ryxxxzz,aj2ryxxyzz,aj2ryxyyzz,aj2ryyyyzz,aj2ryxxzzz,aj2ryxyzzz,aj2ryyyzzz,aj2ryxzzzz,aj2ryyzzzz,aj2ryzzzzz,aj2ryxxxxxx,aj2ryxxxxxy,aj2ryxxxxyy,aj2ryxxxyyy,aj2ryxxyyyy,aj2ryxyyyyy,aj2ryyyyyyy,aj2ryxxxxxz,aj2ryxxxxyz,aj2ryxxxyyz,aj2ryxxyyyz,aj2ryxyyyyz,aj2ryyyyyyz,aj2ryxxxxzz,aj2ryxxxyzz,aj2ryxxyyzz,aj2ryxyyyzz,aj2ryyyyyzz,aj2ryxxxzzz,aj2ryxxyzzz,aj2ryxyyzzz,aj2ryyyyzzz,aj2ryxxzzzz,aj2ryxyzzzz,aj2ryyyzzzz,aj2ryxzzzzz,aj2ryyzzzzz,aj2ryzzzzzz
        real aj2syx,aj2syy,aj2syz,aj2syxx,aj2syxy,aj2syyy,aj2syxz,aj2syyz,aj2syzz,aj2syxxx,aj2syxxy,aj2syxyy,aj2syyyy,aj2syxxz,aj2syxyz,aj2syyyz,aj2syxzz,aj2syyzz,aj2syzzz,aj2syxxxx,aj2syxxxy,aj2syxxyy,aj2syxyyy,aj2syyyyy,aj2syxxxz,aj2syxxyz,aj2syxyyz,aj2syyyyz,aj2syxxzz,aj2syxyzz,aj2syyyzz,aj2syxzzz,aj2syyzzz,aj2syzzzz,aj2syxxxxx,aj2syxxxxy,aj2syxxxyy,aj2syxxyyy,aj2syxyyyy,aj2syyyyyy,aj2syxxxxz,aj2syxxxyz,aj2syxxyyz,aj2syxyyyz,aj2syyyyyz,aj2syxxxzz,aj2syxxyzz,aj2syxyyzz,aj2syyyyzz,aj2syxxzzz,aj2syxyzzz,aj2syyyzzz,aj2syxzzzz,aj2syyzzzz,aj2syzzzzz,aj2syxxxxxx,aj2syxxxxxy,aj2syxxxxyy,aj2syxxxyyy,aj2syxxyyyy,aj2syxyyyyy,aj2syyyyyyy,aj2syxxxxxz,aj2syxxxxyz,aj2syxxxyyz,aj2syxxyyyz,aj2syxyyyyz,aj2syyyyyyz,aj2syxxxxzz,aj2syxxxyzz,aj2syxxyyzz,aj2syxyyyzz,aj2syyyyyzz,aj2syxxxzzz,aj2syxxyzzz,aj2syxyyzzz,aj2syyyyzzz,aj2syxxzzzz,aj2syxyzzzz,aj2syyyzzzz,aj2syxzzzzz,aj2syyzzzzz,aj2syzzzzzz
        real aj2rz,aj2rzr,aj2rzs,aj2rzt,aj2rzrr,aj2rzrs,aj2rzss,aj2rzrt,aj2rzst,aj2rztt,aj2rzrrr,aj2rzrrs,aj2rzrss,aj2rzsss,aj2rzrrt,aj2rzrst,aj2rzsst,aj2rzrtt,aj2rzstt,aj2rzttt,aj2rzrrrr,aj2rzrrrs,aj2rzrrss,aj2rzrsss,aj2rzssss,aj2rzrrrt,aj2rzrrst,aj2rzrsst,aj2rzssst,aj2rzrrtt,aj2rzrstt,aj2rzsstt,aj2rzrttt,aj2rzsttt,aj2rztttt,aj2rzrrrrr,aj2rzrrrrs,aj2rzrrrss,aj2rzrrsss,aj2rzrssss,aj2rzsssss,aj2rzrrrrt,aj2rzrrrst,aj2rzrrsst,aj2rzrssst,aj2rzsssst,aj2rzrrrtt,aj2rzrrstt,aj2rzrsstt,aj2rzssstt,aj2rzrrttt,aj2rzrsttt,aj2rzssttt,aj2rzrtttt,aj2rzstttt,aj2rzttttt,aj2rzrrrrrr,aj2rzrrrrrs,aj2rzrrrrss,aj2rzrrrsss,aj2rzrrssss,aj2rzrsssss,aj2rzssssss,aj2rzrrrrrt,aj2rzrrrrst,aj2rzrrrsst,aj2rzrrssst,aj2rzrsssst,aj2rzssssst,aj2rzrrrrtt,aj2rzrrrstt,aj2rzrrsstt,aj2rzrssstt,aj2rzsssstt,aj2rzrrrttt,aj2rzrrsttt,aj2rzrssttt,aj2rzsssttt,aj2rzrrtttt,aj2rzrstttt,aj2rzsstttt,aj2rzrttttt,aj2rzsttttt,aj2rztttttt
        real aj2sz,aj2szr,aj2szs,aj2szt,aj2szrr,aj2szrs,aj2szss,aj2szrt,aj2szst,aj2sztt,aj2szrrr,aj2szrrs,aj2szrss,aj2szsss,aj2szrrt,aj2szrst,aj2szsst,aj2szrtt,aj2szstt,aj2szttt,aj2szrrrr,aj2szrrrs,aj2szrrss,aj2szrsss,aj2szssss,aj2szrrrt,aj2szrrst,aj2szrsst,aj2szssst,aj2szrrtt,aj2szrstt,aj2szsstt,aj2szrttt,aj2szsttt,aj2sztttt,aj2szrrrrr,aj2szrrrrs,aj2szrrrss,aj2szrrsss,aj2szrssss,aj2szsssss,aj2szrrrrt,aj2szrrrst,aj2szrrsst,aj2szrssst,aj2szsssst,aj2szrrrtt,aj2szrrstt,aj2szrsstt,aj2szssstt,aj2szrrttt,aj2szrsttt,aj2szssttt,aj2szrtttt,aj2szstttt,aj2szttttt,aj2szrrrrrr,aj2szrrrrrs,aj2szrrrrss,aj2szrrrsss,aj2szrrssss,aj2szrsssss,aj2szssssss,aj2szrrrrrt,aj2szrrrrst,aj2szrrrsst,aj2szrrssst,aj2szrsssst,aj2szssssst,aj2szrrrrtt,aj2szrrrstt,aj2szrrsstt,aj2szrssstt,aj2szsssstt,aj2szrrrttt,aj2szrrsttt,aj2szrssttt,aj2szsssttt,aj2szrrtttt,aj2szrstttt,aj2szsstttt,aj2szrttttt,aj2szsttttt,aj2sztttttt
        real aj2tx,aj2txr,aj2txs,aj2txt,aj2txrr,aj2txrs,aj2txss,aj2txrt,aj2txst,aj2txtt,aj2txrrr,aj2txrrs,aj2txrss,aj2txsss,aj2txrrt,aj2txrst,aj2txsst,aj2txrtt,aj2txstt,aj2txttt,aj2txrrrr,aj2txrrrs,aj2txrrss,aj2txrsss,aj2txssss,aj2txrrrt,aj2txrrst,aj2txrsst,aj2txssst,aj2txrrtt,aj2txrstt,aj2txsstt,aj2txrttt,aj2txsttt,aj2txtttt,aj2txrrrrr,aj2txrrrrs,aj2txrrrss,aj2txrrsss,aj2txrssss,aj2txsssss,aj2txrrrrt,aj2txrrrst,aj2txrrsst,aj2txrssst,aj2txsssst,aj2txrrrtt,aj2txrrstt,aj2txrsstt,aj2txssstt,aj2txrrttt,aj2txrsttt,aj2txssttt,aj2txrtttt,aj2txstttt,aj2txttttt,aj2txrrrrrr,aj2txrrrrrs,aj2txrrrrss,aj2txrrrsss,aj2txrrssss,aj2txrsssss,aj2txssssss,aj2txrrrrrt,aj2txrrrrst,aj2txrrrsst,aj2txrrssst,aj2txrsssst,aj2txssssst,aj2txrrrrtt,aj2txrrrstt,aj2txrrsstt,aj2txrssstt,aj2txsssstt,aj2txrrrttt,aj2txrrsttt,aj2txrssttt,aj2txsssttt,aj2txrrtttt,aj2txrstttt,aj2txsstttt,aj2txrttttt,aj2txsttttt,aj2txtttttt
        real aj2ty,aj2tyr,aj2tys,aj2tyt,aj2tyrr,aj2tyrs,aj2tyss,aj2tyrt,aj2tyst,aj2tytt,aj2tyrrr,aj2tyrrs,aj2tyrss,aj2tysss,aj2tyrrt,aj2tyrst,aj2tysst,aj2tyrtt,aj2tystt,aj2tyttt,aj2tyrrrr,aj2tyrrrs,aj2tyrrss,aj2tyrsss,aj2tyssss,aj2tyrrrt,aj2tyrrst,aj2tyrsst,aj2tyssst,aj2tyrrtt,aj2tyrstt,aj2tysstt,aj2tyrttt,aj2tysttt,aj2tytttt,aj2tyrrrrr,aj2tyrrrrs,aj2tyrrrss,aj2tyrrsss,aj2tyrssss,aj2tysssss,aj2tyrrrrt,aj2tyrrrst,aj2tyrrsst,aj2tyrssst,aj2tysssst,aj2tyrrrtt,aj2tyrrstt,aj2tyrsstt,aj2tyssstt,aj2tyrrttt,aj2tyrsttt,aj2tyssttt,aj2tyrtttt,aj2tystttt,aj2tyttttt,aj2tyrrrrrr,aj2tyrrrrrs,aj2tyrrrrss,aj2tyrrrsss,aj2tyrrssss,aj2tyrsssss,aj2tyssssss,aj2tyrrrrrt,aj2tyrrrrst,aj2tyrrrsst,aj2tyrrssst,aj2tyrsssst,aj2tyssssst,aj2tyrrrrtt,aj2tyrrrstt,aj2tyrrsstt,aj2tyrssstt,aj2tysssstt,aj2tyrrrttt,aj2tyrrsttt,aj2tyrssttt,aj2tysssttt,aj2tyrrtttt,aj2tyrstttt,aj2tysstttt,aj2tyrttttt,aj2tysttttt,aj2tytttttt
        real aj2tz,aj2tzr,aj2tzs,aj2tzt,aj2tzrr,aj2tzrs,aj2tzss,aj2tzrt,aj2tzst,aj2tztt,aj2tzrrr,aj2tzrrs,aj2tzrss,aj2tzsss,aj2tzrrt,aj2tzrst,aj2tzsst,aj2tzrtt,aj2tzstt,aj2tzttt,aj2tzrrrr,aj2tzrrrs,aj2tzrrss,aj2tzrsss,aj2tzssss,aj2tzrrrt,aj2tzrrst,aj2tzrsst,aj2tzssst,aj2tzrrtt,aj2tzrstt,aj2tzsstt,aj2tzrttt,aj2tzsttt,aj2tztttt,aj2tzrrrrr,aj2tzrrrrs,aj2tzrrrss,aj2tzrrsss,aj2tzrssss,aj2tzsssss,aj2tzrrrrt,aj2tzrrrst,aj2tzrrsst,aj2tzrssst,aj2tzsssst,aj2tzrrrtt,aj2tzrrstt,aj2tzrsstt,aj2tzssstt,aj2tzrrttt,aj2tzrsttt,aj2tzssttt,aj2tzrtttt,aj2tzstttt,aj2tzttttt,aj2tzrrrrrr,aj2tzrrrrrs,aj2tzrrrrss,aj2tzrrrsss,aj2tzrrssss,aj2tzrsssss,aj2tzssssss,aj2tzrrrrrt,aj2tzrrrrst,aj2tzrrrsst,aj2tzrrssst,aj2tzrsssst,aj2tzssssst,aj2tzrrrrtt,aj2tzrrrstt,aj2tzrrsstt,aj2tzrssstt,aj2tzsssstt,aj2tzrrrttt,aj2tzrrsttt,aj2tzrssttt,aj2tzsssttt,aj2tzrrtttt,aj2tzrstttt,aj2tzsstttt,aj2tzrttttt,aj2tzsttttt,aj2tztttttt
        real aj2rzx,aj2rzy,aj2rzz,aj2rzxx,aj2rzxy,aj2rzyy,aj2rzxz,aj2rzyz,aj2rzzz,aj2rzxxx,aj2rzxxy,aj2rzxyy,aj2rzyyy,aj2rzxxz,aj2rzxyz,aj2rzyyz,aj2rzxzz,aj2rzyzz,aj2rzzzz,aj2rzxxxx,aj2rzxxxy,aj2rzxxyy,aj2rzxyyy,aj2rzyyyy,aj2rzxxxz,aj2rzxxyz,aj2rzxyyz,aj2rzyyyz,aj2rzxxzz,aj2rzxyzz,aj2rzyyzz,aj2rzxzzz,aj2rzyzzz,aj2rzzzzz,aj2rzxxxxx,aj2rzxxxxy,aj2rzxxxyy,aj2rzxxyyy,aj2rzxyyyy,aj2rzyyyyy,aj2rzxxxxz,aj2rzxxxyz,aj2rzxxyyz,aj2rzxyyyz,aj2rzyyyyz,aj2rzxxxzz,aj2rzxxyzz,aj2rzxyyzz,aj2rzyyyzz,aj2rzxxzzz,aj2rzxyzzz,aj2rzyyzzz,aj2rzxzzzz,aj2rzyzzzz,aj2rzzzzzz,aj2rzxxxxxx,aj2rzxxxxxy,aj2rzxxxxyy,aj2rzxxxyyy,aj2rzxxyyyy,aj2rzxyyyyy,aj2rzyyyyyy,aj2rzxxxxxz,aj2rzxxxxyz,aj2rzxxxyyz,aj2rzxxyyyz,aj2rzxyyyyz,aj2rzyyyyyz,aj2rzxxxxzz,aj2rzxxxyzz,aj2rzxxyyzz,aj2rzxyyyzz,aj2rzyyyyzz,aj2rzxxxzzz,aj2rzxxyzzz,aj2rzxyyzzz,aj2rzyyyzzz,aj2rzxxzzzz,aj2rzxyzzzz,aj2rzyyzzzz,aj2rzxzzzzz,aj2rzyzzzzz,aj2rzzzzzzz
        real aj2szx,aj2szy,aj2szz,aj2szxx,aj2szxy,aj2szyy,aj2szxz,aj2szyz,aj2szzz,aj2szxxx,aj2szxxy,aj2szxyy,aj2szyyy,aj2szxxz,aj2szxyz,aj2szyyz,aj2szxzz,aj2szyzz,aj2szzzz,aj2szxxxx,aj2szxxxy,aj2szxxyy,aj2szxyyy,aj2szyyyy,aj2szxxxz,aj2szxxyz,aj2szxyyz,aj2szyyyz,aj2szxxzz,aj2szxyzz,aj2szyyzz,aj2szxzzz,aj2szyzzz,aj2szzzzz,aj2szxxxxx,aj2szxxxxy,aj2szxxxyy,aj2szxxyyy,aj2szxyyyy,aj2szyyyyy,aj2szxxxxz,aj2szxxxyz,aj2szxxyyz,aj2szxyyyz,aj2szyyyyz,aj2szxxxzz,aj2szxxyzz,aj2szxyyzz,aj2szyyyzz,aj2szxxzzz,aj2szxyzzz,aj2szyyzzz,aj2szxzzzz,aj2szyzzzz,aj2szzzzzz,aj2szxxxxxx,aj2szxxxxxy,aj2szxxxxyy,aj2szxxxyyy,aj2szxxyyyy,aj2szxyyyyy,aj2szyyyyyy,aj2szxxxxxz,aj2szxxxxyz,aj2szxxxyyz,aj2szxxyyyz,aj2szxyyyyz,aj2szyyyyyz,aj2szxxxxzz,aj2szxxxyzz,aj2szxxyyzz,aj2szxyyyzz,aj2szyyyyzz,aj2szxxxzzz,aj2szxxyzzz,aj2szxyyzzz,aj2szyyyzzz,aj2szxxzzzz,aj2szxyzzzz,aj2szyyzzzz,aj2szxzzzzz,aj2szyzzzzz,aj2szzzzzzz
        real aj2txx,aj2txy,aj2txz,aj2txxx,aj2txxy,aj2txyy,aj2txxz,aj2txyz,aj2txzz,aj2txxxx,aj2txxxy,aj2txxyy,aj2txyyy,aj2txxxz,aj2txxyz,aj2txyyz,aj2txxzz,aj2txyzz,aj2txzzz,aj2txxxxx,aj2txxxxy,aj2txxxyy,aj2txxyyy,aj2txyyyy,aj2txxxxz,aj2txxxyz,aj2txxyyz,aj2txyyyz,aj2txxxzz,aj2txxyzz,aj2txyyzz,aj2txxzzz,aj2txyzzz,aj2txzzzz,aj2txxxxxx,aj2txxxxxy,aj2txxxxyy,aj2txxxyyy,aj2txxyyyy,aj2txyyyyy,aj2txxxxxz,aj2txxxxyz,aj2txxxyyz,aj2txxyyyz,aj2txyyyyz,aj2txxxxzz,aj2txxxyzz,aj2txxyyzz,aj2txyyyzz,aj2txxxzzz,aj2txxyzzz,aj2txyyzzz,aj2txxzzzz,aj2txyzzzz,aj2txzzzzz,aj2txxxxxxx,aj2txxxxxxy,aj2txxxxxyy,aj2txxxxyyy,aj2txxxyyyy,aj2txxyyyyy,aj2txyyyyyy,aj2txxxxxxz,aj2txxxxxyz,aj2txxxxyyz,aj2txxxyyyz,aj2txxyyyyz,aj2txyyyyyz,aj2txxxxxzz,aj2txxxxyzz,aj2txxxyyzz,aj2txxyyyzz,aj2txyyyyzz,aj2txxxxzzz,aj2txxxyzzz,aj2txxyyzzz,aj2txyyyzzz,aj2txxxzzzz,aj2txxyzzzz,aj2txyyzzzz,aj2txxzzzzz,aj2txyzzzzz,aj2txzzzzzz
        real aj2tyx,aj2tyy,aj2tyz,aj2tyxx,aj2tyxy,aj2tyyy,aj2tyxz,aj2tyyz,aj2tyzz,aj2tyxxx,aj2tyxxy,aj2tyxyy,aj2tyyyy,aj2tyxxz,aj2tyxyz,aj2tyyyz,aj2tyxzz,aj2tyyzz,aj2tyzzz,aj2tyxxxx,aj2tyxxxy,aj2tyxxyy,aj2tyxyyy,aj2tyyyyy,aj2tyxxxz,aj2tyxxyz,aj2tyxyyz,aj2tyyyyz,aj2tyxxzz,aj2tyxyzz,aj2tyyyzz,aj2tyxzzz,aj2tyyzzz,aj2tyzzzz,aj2tyxxxxx,aj2tyxxxxy,aj2tyxxxyy,aj2tyxxyyy,aj2tyxyyyy,aj2tyyyyyy,aj2tyxxxxz,aj2tyxxxyz,aj2tyxxyyz,aj2tyxyyyz,aj2tyyyyyz,aj2tyxxxzz,aj2tyxxyzz,aj2tyxyyzz,aj2tyyyyzz,aj2tyxxzzz,aj2tyxyzzz,aj2tyyyzzz,aj2tyxzzzz,aj2tyyzzzz,aj2tyzzzzz,aj2tyxxxxxx,aj2tyxxxxxy,aj2tyxxxxyy,aj2tyxxxyyy,aj2tyxxyyyy,aj2tyxyyyyy,aj2tyyyyyyy,aj2tyxxxxxz,aj2tyxxxxyz,aj2tyxxxyyz,aj2tyxxyyyz,aj2tyxyyyyz,aj2tyyyyyyz,aj2tyxxxxzz,aj2tyxxxyzz,aj2tyxxyyzz,aj2tyxyyyzz,aj2tyyyyyzz,aj2tyxxxzzz,aj2tyxxyzzz,aj2tyxyyzzz,aj2tyyyyzzz,aj2tyxxzzzz,aj2tyxyzzzz,aj2tyyyzzzz,aj2tyxzzzzz,aj2tyyzzzzz,aj2tyzzzzzz
        real aj2tzx,aj2tzy,aj2tzz,aj2tzxx,aj2tzxy,aj2tzyy,aj2tzxz,aj2tzyz,aj2tzzz,aj2tzxxx,aj2tzxxy,aj2tzxyy,aj2tzyyy,aj2tzxxz,aj2tzxyz,aj2tzyyz,aj2tzxzz,aj2tzyzz,aj2tzzzz,aj2tzxxxx,aj2tzxxxy,aj2tzxxyy,aj2tzxyyy,aj2tzyyyy,aj2tzxxxz,aj2tzxxyz,aj2tzxyyz,aj2tzyyyz,aj2tzxxzz,aj2tzxyzz,aj2tzyyzz,aj2tzxzzz,aj2tzyzzz,aj2tzzzzz,aj2tzxxxxx,aj2tzxxxxy,aj2tzxxxyy,aj2tzxxyyy,aj2tzxyyyy,aj2tzyyyyy,aj2tzxxxxz,aj2tzxxxyz,aj2tzxxyyz,aj2tzxyyyz,aj2tzyyyyz,aj2tzxxxzz,aj2tzxxyzz,aj2tzxyyzz,aj2tzyyyzz,aj2tzxxzzz,aj2tzxyzzz,aj2tzyyzzz,aj2tzxzzzz,aj2tzyzzzz,aj2tzzzzzz,aj2tzxxxxxx,aj2tzxxxxxy,aj2tzxxxxyy,aj2tzxxxyyy,aj2tzxxyyyy,aj2tzxyyyyy,aj2tzyyyyyy,aj2tzxxxxxz,aj2tzxxxxyz,aj2tzxxxyyz,aj2tzxxyyyz,aj2tzxyyyyz,aj2tzyyyyyz,aj2tzxxxxzz,aj2tzxxxyzz,aj2tzxxyyzz,aj2tzxyyyzz,aj2tzyyyyzz,aj2tzxxxzzz,aj2tzxxyzzz,aj2tzxyyzzz,aj2tzyyyzzz,aj2tzxxzzzz,aj2tzxyzzzz,aj2tzyyzzzz,aj2tzxzzzzz,aj2tzyzzzzz,aj2tzzzzzzz
 !............... end statement functions
       ierr=0
       side1                =ipar(0)
       axis1                =ipar(1)
       grid1                =ipar(2)
       n1a                  =ipar(3)
       n1b                  =ipar(4)
       n2a                  =ipar(5)
       n2b                  =ipar(6)
       n3a                  =ipar(7)
       n3b                  =ipar(8)
       side2                =ipar(9)
       axis2                =ipar(10)
       grid2                =ipar(11)
       m1a                  =ipar(12)
       m1b                  =ipar(13)
       m2a                  =ipar(14)
       m2b                  =ipar(15)
       m3a                  =ipar(16)
       m3b                  =ipar(17)
       gridType             =ipar(18)
       orderOfAccuracy      =ipar(19)
       orderOfExtrapolation =ipar(20)  ! maximum allowable order of extrapolation
       useForcing           =ipar(21)
       ex                   =ipar(22)
       ey                   =ipar(23)
       ez                   =ipar(24)
       hx                   =ipar(25)
       hy                   =ipar(26)
       hz                   =ipar(27)
       solveForE            =ipar(28)
       solveForH            =ipar(29)
       useWhereMask         =ipar(30)
       debug                =ipar(31)
       nit                  =ipar(32)
       interfaceOption      =ipar(33)
       initialized          =ipar(34)
       myid                 =ipar(35)
       parallel             =ipar(36)
       forcingOption        =ipar(37)
       interfaceEquationsOption  =ipar(38)
       assignInterfaceValues     =ipar(39)
       assignInterfaceGhostValues=ipar(40)
       setDivergenceAtInterfaces =ipar(41)
       ! *new* *wdh* June 28, 2016
       !  useImpedanceInterfaceProjection=0: OLD way 
       !  useImpedanceInterfaceProjection=1: new way using impedance weighting (see maxwell.pdf)
       useImpedanceInterfaceProjection=ipar(42)
       ! numberOfInterfaceIterationsUsed = ipar(43)  ! returned value 
       ipar(43)=0
       dispersionModel1    = ipar(44)
       dispersionModel2    = ipar(45)
       pxc                 = ipar(46)
       knownSolutionOption = ipar(47)
       useJacobiUpdate     = ipar(48)
       nonlinearModel1     = ipar(49)
       nonlinearModel2     = ipar(50)
       numParallelGhost    = ipar(51)
       internalGhostBC     = ipar(52)  ! bc value for internal parallel boundaries 
       useUnifiedInterfaceMacros = ipar(53)
       dx1(0)                =rpar(0)
       dx1(1)                =rpar(1)
       dx1(2)                =rpar(2)
       dr1(0)                =rpar(3)
       dr1(1)                =rpar(4)
       dr1(2)                =rpar(5)
       dx2(0)                =rpar(6)
       dx2(1)                =rpar(7)
       dx2(2)                =rpar(8)
       dr2(0)                =rpar(9)
       dr2(1)                =rpar(10)
       dr2(2)                =rpar(11)
       t                    =rpar(12)
       ep                   =rpar(13) ! pointer for exact solution
       dt                   =rpar(14)
       eps1                 =rpar(15)
       mu1                  =rpar(16)
       c1                   =rpar(17)
       eps2                 =rpar(18)
       mu2                  =rpar(19)
       c2                   =rpar(20)
       omega                =rpar(21)
       ! rpar(22) : averageInterfaceConvergenceRate : return value 
       ! rpar(23) : maxFinalResidual : return value 
       relativeErrorTolerance = rpar(24)
       absoluteErrorTolerance = rpar(25)
       ! absoluteErrorTolerance=1.e-10  ! fix me -- need a relative tol
       ! absoluteErrorTolerance=1.e-15  ! fix me -- need a relative tol
       cSq1=c1**2
       cSq2=c2**2
       dtSq=dt**2
       eta1=sqrt(mu1/eps1) ! electrical impedance
       eta2=sqrt(mu2/eps2) ! electrical impedance
       eta1i=1./eta1
       eta2i=1./eta2
       epsmu1=eps1*mu1
       epsmu2=eps2*mu2
       twilightZone=useForcing
       ! For updating P on ghost using PDE: 
       addForcing=0 ! fix me *******************************
       if( twilightZone.ne.0 )then
         addForcing=twilightZoneForcing
       end if
       avoidInterfaceIterations=1  ! option to avoid interface iterations. *new* July 6 2019 
       dsBig = 1.e20               ! Large value for tangential grid spacing to zero out mixed derivatives
       if( avoidInterfaceIterations.eq.1 )then     
         nit=1
         omega=1.    ! do not under-relax iterations
         axis1p1=mod(axis1+1,3)
         axis1p2=mod(axis1+2,3)
         axis2p1=mod(axis2+1,3)
         axis2p2=mod(axis2+2,3)
         ! dr1a(0:2) = sets grid spacing in tangential directions to dsBig 
         dr1a(axis1  )=dr1(axis1)
         dr1a(axis1p1)=dsBig
         dr1a(axis1p2)=dsBig
         ! dr2a(0:2) = sets grid spacing in tangential directions to dsBig 
         dr2a(axis2  )=dr2(axis2)
         dr2a(axis2p1)=dsBig
         dr2a(axis2p2)=dsBig
       end if
       checkCoeff=0 ! 0 ! 1  ! if non-zero then check coefficients in interface equations using delta approach
       if( t.gt.1.5*dt )then
          checkCoeff=0    ! turn off for t>0 
       end if 
       coeffDiff=0.
       saveCoeff=0  ! 1 : save coefficients in the interface equations to a file
       coeffFile=20  ! file for saving coefficients
       ieqn=-1      ! counts equations
       if( saveCoeff.eq.1 .and. initialized.eq.0 .and. assignInterfaceGhostValues.eq.1 )then
         open (coeffFile,file='interfaceCoeff.dat',status='unknown',form='formatted')
         write(coeffFile,'("! Coefficients in the interface equations")') 
         write(coeffFile,'("! nd,orderOfAccuracy")') 
         write(coeffFile,'(2(i6,1x))') nd,orderOfAccuracy
         write(coeffFile,'(9(i6,1x))') side1,axis1,grid1, n1a,n1b,n2a,n2b,n3a,n3b
         write(coeffFile,'(9(i6,1x))') side2,axis2,grid2, m1a,m1b,m2a,m2b,m3a,m3b
       end if
       debugFile=10
       if( initialized.eq.0 .and. debug.gt.0 )then
         ! open debug files
         ! open (debugFile,file=filen,status='unknown',form='formatted')
         if( myid.lt.10 )then
           write(debugFileName,'("mxi",i1,".fdebug")') myid
         else
           write(debugFileName,'("mxi",i2,".fdebug")') myid
         end if
         write(*,*) 'interface3d: myid=',myid,' open debug file:',debugFileName
         open (debugFile,file=debugFileName,status='unknown',form='formatted')
         ! '
         ! INQUIRE(FILE=filen, EXIST=filex)
       end if
       if( t.le. 1.5*dt .and. debug.gt.0 )then
         write(*,'(" +++++++++cgmx interface3d t=",e9.2," dt=",e9.2," nit=",i3," ++++++++")') t,dt,nit
            ! '
         write(*,'("  ... nd=",i2," gridType=",i2," order=",i2," debug=",i3,", ex=",i2)') nd,gridType,orderOfAccuracy,debug,ex
         write(*,'("  ... assignInterface=",i2," assignGhost=",i2)') assignInterfaceValues,assignInterfaceGhostValues
         write(*,'("  ... useJacobiUpdate=",i2," numParallelGhost=",i2)') useJacobiUpdate,numParallelGhost
         write(*,'("  ... setDivergenceAtInterfaces=",i2)') setDivergenceAtInterfaces
         write(*,'("  ... useImpedanceInterfaceProjection=",i2)') useImpedanceInterfaceProjection
         write(*,'("  ... avoidInterfaceIterations=",i2)') avoidInterfaceIterations
         write(*,'("  ... useUnifiedInterfaceMacros=",i2)') useUnifiedInterfaceMacros
         write(*,'("  ... interface its (4th-order) relativeTol=",e12.3," absoluteTol=",e12.3)') relativeErrorTolerance,absoluteErrorTolerance
         write(*,'("  ... bc1=",6i6)') ((boundaryCondition1(i1,i2),i1=0,1),i2=0,nd-1)
         write(*,'("  ... bc2=",6i6)') ((boundaryCondition2(i1,i2),i1=0,1),i2=0,nd-1)
         write(*,'("  ... internalGhostBC=",i6)') internalGhostBC
       end if
       if( t.lt.1.5*dt .and. debug.gt.0 )then
         write(debugFile,'(" +++++++++cgmx interface3d t=",e9.2," ++++++++")') t
            ! '
         write(debugFile,'(" interface3d: nd=",i2," gridType=",i2)') nd,gridType
         write(debugFile,'("  ... nd=",i2," gridType=",i2," order=",i2," debug=",i3)') nd,gridType,orderOfAccuracy,debug
         write(debugFile,'("  ... assignInterface=",i2," assignGhost=",i2)') assignInterfaceValues,assignInterfaceGhostValues
       end if
       if( abs(c1*c1-1./(mu1*eps1)).gt. 1.e-10 )then
         write(debugFile,'(" interface3d:ERROR: c1,eps1,mu1=",3e10.2," not consistent")') c1,eps1,mu1
            ! '
         stop 11
       end if
       if( abs(c2*c2-1./(mu2*eps2)).gt. 1.e-10 )then
         write(debugFile,'(" interface3d:ERROR: c2,eps2,mu2=",3e10.2," not consistent")') c2,eps2,mu2
            ! '
         stop 11
       end if
       if( .false. )then
         write(debugFile,'(" interface3d: eps1,eps2=",2f10.5," c1,c2=",2f10.5)') eps1,eps2,c1,c2
            ! '
       end if
       if( nit.lt.0 .or. nit.gt.100 )then
         write(debugFile,'(" interfaceBC: ERROR: nit=",i9)') nit
         nit=max(1,min(100,nit))
       end if
       if( debug.gt.1 )then
         write(debugFile,'("********************************************************************** ")')
         write(debugFile,'(" interface3d: **START** t=",e10.2)') t
         write(debugFile,'(" interface3d: **START** grid1=",i4," side1,axis1=",2i2," bc=",6i3)') grid1,side1,axis1,boundaryCondition1(0,0),boundaryCondition1(1,0),boundaryCondition1(0,1),boundaryCondition1(1,1),boundaryCondition1(0,2),boundaryCondition1(1,2)
            ! '
         write(debugFile,'(" interface3d: **START** grid2=",i4," side2,axis2=",2i2," bc=",6i3)') grid2,side2,axis2,boundaryCondition2(0,0),boundaryCondition2(1,0),boundaryCondition2(0,1),boundaryCondition2(1,1),boundaryCondition2(0,2),boundaryCondition2(1,2)
            ! '
         write(debugFile,'("n1a,n1b,...=",6i5)') n1a,n1b,n2a,n2b,n3a,n3b
         write(debugFile,'("m1a,m1b,...=",6i5)') m1a,m1b,m2a,m2b,m3a,m3b
       end if
       if( debug.gt.8 )then
        write(debugFile,'("start u1=",(3i4,1x,3e11.2))') (((i1,i2,i3,(u1(i1,i2,i3,m),m=0,2),i1=nd1a,nd1b),i2=nd2a,nd2b),i3=nd3a,nd3b)
        write(debugFile,'("start u2=",(3i4,1x,3e11.2))') (((i1,i2,i3,(u2(i1,i2,i3,m),m=0,2),i1=md1a,md1b),i2=md2a,md2b),i3=md3a,md3b)
       end if
       dispersive=0 ! >0 implies at least one side is dispersive
       gdmParOption=1 ! scale a0 and a1 parameters by eps
       if( dispersionModel1.ne.noDispersion )then
         dispersive=dispersive+1
         ! get the gdm parameters
         !   gdmPar(0:3,jv) = (a0,a1,b0,b1) 
         call getGDMParameters( grid1,alphaP1,gdmPar1,numberOfPolarizationVectors1, maxNumberOfParameters,maxNumberOfPolarizationVectors,gdmParOption )
         if( t.le. 1.5*dt .and. debug.gt.0 )then
           ! ---- Dispersive Maxwell ----
           write(*,'("--interface3d-- dispersionModel1=",i4," grid1=",i4," pxc=",i4)') dispersionModel1,grid1,pxc
           write(*,'("--interface3d-- GDM: numberOfPolarizationVectors1=",i4," alphaP1=",e8.2)') numberOfPolarizationVectors1,alphaP1
           do jv=0,numberOfPolarizationVectors1-1
             write(*,'("--interface3d-- GDM: eqn=",i3," a0,a1,b0,b1=",4(1p,e10.2))') jv,a0v1(jv),a1v1(jv),b0v1(jv),b1v1(jv)
           end do 
         end if
       else
         numberOfPolarizationVectors1=0
         alphaP1 = 0.
       end if
       if( dispersionModel2.ne.noDispersion )then
         dispersive=dispersive+1
         ! get the gdm parameters
         !   gdmPar(0:3,jv) = (a0,a1,b0,b1) 
         call getGDMParameters( grid2,alphaP2,gdmPar2,numberOfPolarizationVectors2, maxNumberOfParameters,maxNumberOfPolarizationVectors,gdmParOption )
         if( t.le. 1.5*dt .and. debug.gt.0 )then
           ! ---- Dispersive Maxwell ----
           write(*,'("--interface3d-- dispersionModel2=",i4," grid2=",i4)') dispersionModel2,grid2
           write(*,'("--interface3d-- GDM: numberOfPolarizationVectors2=",i4," alphaP2=",e8.2)') numberOfPolarizationVectors2,alphaP2
           do jv=0,numberOfPolarizationVectors2-1
             write(*,'("--interface3d-- GDM: eqn=",i3," a0,a1,b0,b1=",4(1p,e10.2))') jv,a0v2(jv),a1v2(jv),b0v2(jv),b1v2(jv)
           end do 
         end if
         ! write(*,'(" interface: FINISH ME")') 
         ! stop 1111
       else
         numberOfPolarizationVectors2=0
         alphaP2 = 0.
       end if
       useNonlinearModel=0
       if( nonlinearModel1 .ne. noNonlinearModel )then
         useNonlinearModel=1
         call getMultilevelAtomicParameters( grid1, nlPar1, maxPar, maxPar, numPolar1, numberOfAtomicLevels1 )
         if( numPolar1.ne.numberOfPolarizationVectors1 )then
           write(*,'(" interfaceOpt:ERROR: numberOfPolarizationVectors1 does not match numPolar1 from nonlinear model!!")')
           stop 8888
         end if
         if( t.le. 1.5*dt .and. debug.gt.0 )then
           write(*,'("--interfaceOpt-- nonlinearModel1=",i4," (1=multilevelAtomic)")') nonlinearModel1
           write(*,'("multilevelAtomic: numberOfPolarizationVectors1=",i4,"  numberOfAtomicLevels1=",i4)') numberOfPolarizationVectors1, numberOfAtomicLevels1
           write(*,'("polarizationNECoefficients1:")')
           do m1=0,numberOfPolarizationVectors1-1
             write(*,'( 10(e12.3,1x) )') (pnec1(m1,m2),m2=0,numberOfAtomicLevels1-1)
           end do 
           write(*,'("populationRelaxationCoefficients1:")')
           do m1=0,numberOfAtomicLevels1-1
             write(*,'( 10(e12.3,1x) )') (prc1(m1,m2),m2=0,numberOfAtomicLevels1-1)
           end do 
           write(*,'("populationEPtCoefficients1:")')
           do m1=0,numberOfAtomicLevels1-1
             write(*,'( 10(e12.3,1x) )') (peptc1(m1,m2),m2=0,numberOfPolarizationVectors1-1)
           end do 
         end if
       else
         numberOfAtomicLevels1=0
       end if 
       if( nonlinearModel2 .ne. noNonlinearModel )then
         useNonlinearModel=1
         call getMultilevelAtomicParameters( grid2, nlPar2, maxPar, maxPar, numPolar2, numberOfAtomicLevels2 )
         if( numPolar2.ne.numberOfPolarizationVectors2 )then
           write(*,'(" interfaceOpt:ERROR: numberOfPolarizationVectors2 does not match numPolar2 from nonlinear model!!")')
           stop 9999
         end if
         if( t.le. 1.5*dt .and. debug.gt.0 )then
           write(*,'("--interfaceOpt-- nonlinearModel2=",i4," (1=multilevelAtomic)")') nonlinearModel2
           write(*,'("multilevelAtomic: numberOfPolarizationVectors2=",i4,"  numberOfAtomicLevels2=",i4)') numberOfPolarizationVectors2, numberOfAtomicLevels2
           write(*,'("polarizationNECoefficients2:")')
           do m1=0,numberOfPolarizationVectors2-1
             write(*,'( 10(e12.3,1x) )') (pnec2(m1,m2),m2=0,numberOfAtomicLevels2-1)
           end do 
           write(*,'("populationRelaxationCoefficients2:")')
           do m1=0,numberOfAtomicLevels2-1
             write(*,'( 10(e12.3,1x) )') (prc2(m1,m2),m2=0,numberOfAtomicLevels2-1)
           end do 
           write(*,'("populationEPtCoefficients2:")')
           do m1=0,numberOfAtomicLevels2-1
             write(*,'( 10(e12.3,1x) )') (peptc2(m1,m2),m2=0,numberOfPolarizationVectors2-1)
           end do 
         end if
       else
         numberOfAtomicLevels2 = 0
       end if 
       epsx=1.e-20  ! fix this 
       ! --- init various variables and loop bounds ----
             do kd=0,nd-1
              dx112(kd) = 1./(2.*dx1(kd))
              dx122(kd) = 1./(dx1(kd)**2)
              dx212(kd) = 1./(2.*dx2(kd))
              dx222(kd) = 1./(dx2(kd)**2)
              dx141(kd) = 1./(12.*dx1(kd))
              dx142(kd) = 1./(12.*dx1(kd)**2)
              dx241(kd) = 1./(12.*dx2(kd))
              dx242(kd) = 1./(12.*dx2(kd)**2)
              dr114(kd) = 1./(12.*dr1(kd))
              dr214(kd) = 1./(12.*dr2(kd))
             end do
             numGhost=orderOfAccuracy/2
             giveDiv=0   ! set to 1 to give div(u) on both sides, rather than setting the jump in div(u)
             ! save n1a,n1b,... for 2-stage fourth-order scheme
             ns1a=n1a
             ns1b=n1b
             ns2a=n2a
             ns2b=n2b
             ns3a=n3a
             ns3b=n3b
             ms1a=m1a
             ms1b=m1b
             ms2a=m2a
             ms2b=m2b
             ms3a=m3a
             ms3b=m3b
             ! For 2nd-order Stage 1 of fourth-order scheme 
             if( nd.eq.2 )then
               numGhost2=1 ! we always do at least this many 
             else
               numGhost2=0 ! FIX 3D Order 4 -- do this for backward compatibility
             end if 
             if( numParallelGhost.eq.3 )then
               numGhost2=1  ! num ghost for 2nd-order update (Stage I of 4th order update)
             else if( numParallelGhost.gt.3 )then
               numGhost2=2  ! include an extra ghost -- is this needed ? 
             else
               if( orderOfAccuracy.eq.4 .and. numParallelGhost.gt.0 )then
                 if( t.le. 1.5*dt .and. debug.gt.0 )then
                   write(*,'(/,"---------------------------------------------------------------")')
                   write(*,'("interface3d:WARNING: orderOfAccuracy=",i2," but numParallelGhost=",i3)') orderOfAccuracy,numParallelGhost
                   write(*,'("interface3d: Choose numParallelGhost==3 to make answers match to np=1 results")') 
                   write(*,'("---------------------------------------------------------------",/)')
                 end if 
               end if
             end if 
             ne1a=n1a
             ne1b=n1b
             ne2a=n2a
             ne2b=n2b
             ne3a=n3a
             ne3b=n3b
             me1a=m1a
             me1b=m1b
             me2a=m2a
             me2b=m2b
             me3a=m3a
             me3b=m3b
             ! bounds for loops that include ghost points in the tangential directions:
             nn1a=n1a
             nn1b=n1b
             nn2a=n2a
             nn2b=n2b
             nn3a=n3a
             nn3b=n3b
             mm1a=m1a
             mm1b=m1b
             mm2a=m2a
             mm2b=m2b
             mm3a=m3a
             mm3b=m3b
             i3=n3a
             j3=m3a
             axis1p1=mod(axis1+1,nd)
             axis1p2=mod(axis1+2,nd)
             axis2p1=mod(axis2+1,nd)
             axis2p2=mod(axis2+2,nd)
             is1=0
             is2=0
             is3=0
             if( axis1.ne.0 )then
               ! include ghost lines in tangential periodic (and parallel) directions (for extrapolating)
               ! *wdh* Also include ghost on interpolation boundaries 2015/06/29 
               if( boundaryCondition1(0,0).le.0 )then ! parallel ghost may only have bc<0 on one side
                 nn1a=nn1a-numGhost
                 if( boundaryCondition2(0,0).gt.0 )then
                   write(*,'("interface3d: bc is inconsistent")')
                   stop 178
                 end if
               end if
               if( boundaryCondition1(1,0).le.0 )then ! parallel ghost may only have bc<0 on one side
                 nn1b=nn1b+numGhost
                 if( boundaryCondition2(1,0).gt.0 )then
                   write(*,'("interface3d: bc is inconsistent")')
                   stop 179
                 end if
               end if
               ! -- bounds for order 2 stage I of fourth-order scheme 
               if( .true. .or. boundaryCondition1(0,0).eq.internalGhostBC )then
                 ! parallel ghost: 
                 ne1a=ne1a-numGhost2
               end if 
               if( .true. .or. boundaryCondition1(1,0).eq.internalGhostBC )then
                 ! parallel ghost:  
                 ne1b=ne1b+numGhost2
               end if 
             end if
             if( axis1.ne.1 )then
               ! include ghost lines in tangential periodic (and parallel) directions (for extrapolating)
               if( boundaryCondition1(0,1).le.0 )then
                 nn2a=nn2a-numGhost
                 if( boundaryCondition2(0,1).gt.0 )then
                   write(*,'("interface3d: bc is inconsistent")')
                   stop 180
                 end if
               end if
               if( boundaryCondition1(1,1).le.0 )then
                 nn2b=nn2b+numGhost
                 if( boundaryCondition2(1,1).gt.0 )then
                   write(*,'("interface3d: bc is inconsistent")')
                   stop 181
                 end if
               end if
               ! -- bounds for order 2 stage I of fourth-order scheme 
               if( .true. .or. boundaryCondition1(0,1).eq.internalGhostBC )then
                 ! adjust for parallel ghost 
                 ne2a=ne2a-numGhost2
               end if 
               if( .true. .or. boundaryCondition1(1,1).eq.internalGhostBC )then
                 ! adjust for parallel ghost 
                 ne2b=ne2b+numGhost2
               end if
             end if
             if( nd.eq.3 .and. axis1.ne.2 )then
               ! include ghost lines in tangential periodic (and parallel) directions (for extrapolating)
               if( boundaryCondition1(0,2).le.0 )then
                 nn3a=nn3a-numGhost
                 if( boundaryCondition2(0,2).gt.0 )then
                   write(*,'("interface3d: bc is inconsistent")')
                   stop 182
                 end if
               end if
               if( boundaryCondition1(1,2).le.0 )then
                 nn3b=nn3b+numGhost
                 if( boundaryCondition2(1,2).gt.0 )then
                   write(*,'("interface3d: bc is inconsistent")')
                   stop 183
                 end if
               end if
               ! -- bounds for order 2 stage I of fourth-order scheme 
               if( .true. .or. boundaryCondition1(0,2).eq.internalGhostBC )then
                 ! adjust for parallel ghost 
                 ne3a=ne3a-numGhost2
               end if
               if( .true. .or. boundaryCondition1(1,2).eq.internalGhostBC )then
                 ! adjust for parallel ghost 
                 ne3b=ne3b+numGhost2
               end if 
             end if
             if( axis1.eq.0 ) then
               is1=1-2*side1
               an1Cartesian=1. ! normal for a cartesian grid
               an2Cartesian=0.
               an3Cartesian=0.
             else if( axis1.eq.1 )then
               is2=1-2*side1
               an1Cartesian=0.
               an2Cartesian=1.
               an3Cartesian=0.
             else if( axis1.eq.2 )then
               is3=1-2*side1
               an1Cartesian=0.
               an2Cartesian=0.
               an3Cartesian=1.
             else
               stop 5528
             end if
             js1=0
             js2=0
             js3=0
             if( axis2.ne.0 )then
               if( boundaryCondition2(0,0).le.0 )then
                 mm1a=mm1a-numGhost
               end if
               if( boundaryCondition2(1,0).le.0 )then
                 mm1b=mm1b+numGhost
               end if
               if( .true. .or. boundaryCondition2(0,0).eq.internalGhostBC )then
                 me1a=me1a-numGhost2
               end if
               if( .true. .or. boundaryCondition2(1,0).eq.internalGhostBC )then
                 me1b=me1b+numGhost2
               end if
             end if
             if( axis2.ne.1 )then
               if( boundaryCondition2(0,1).le.0 )then
                 mm2a=mm2a-numGhost
               end if
               if( boundaryCondition2(1,1).le.0 )then
                 mm2b=mm2b+numGhost
               end if
               if( .true. .or. boundaryCondition2(0,1).eq.internalGhostBC )then
                 me2a=me2a-numGhost2
               end if
               if( .true. .or. boundaryCondition2(1,1).eq.internalGhostBC )then
                 me2b=me2b+numGhost2
               end if
             end if
             if( nd.eq.3 .and. axis2.ne.2 )then
               if( boundaryCondition2(0,2).le.0 )then
                 mm3a=mm3a-numGhost
               end if
               if( boundaryCondition2(1,2).le.0 )then
                 mm3b=mm3b+numGhost
               end if
               if( .true. .or. boundaryCondition2(0,2).eq.internalGhostBC )then
                 me3a=me3a-numGhost2
               end if
               if( .true. .or. boundaryCondition2(1,2).eq.internalGhostBC )then
                 me3b=me3b+numGhost2
               end if
             end if
             if( axis2.eq.0 ) then
               js1=1-2*side2
             else if( axis2.eq.1 ) then
               js2=1-2*side2
             else  if( axis2.eq.2 ) then
               js3=1-2*side2
             else
               stop 3384
             end if
             is=1-2*side1
             js=1-2*side2
             if( t.le. 1.5*dt .and. debug.gt.0 )then
               write(debugFile,'("myid=",i3," nd1a,nd1b,...  =",6i5)') myid,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b
               write(debugFile,'("myid=",i3," n1a,n1b,...    =",6i5)') myid,n1a,n1b,n2a,n2b,n3a,n3b
               write(debugFile,'("myid=",i3," ne1a,ne1b,...  =",6i5)') myid,ne1a,ne1b,ne2a,ne2b,ne3a,ne3b
               write(debugFile,'("myid=",i3," md1a,md1b,...  =",6i5)') myid,md1a,md1b,md2a,md2b,md3a,md3b
               write(debugFile,'("myid=",i3," m1a,m1b,...    =",6i5)') myid,m1a,m1b,m2a,m2b,m3a,m3b
               write(debugFile,'("myid=",i3," me1a,me1b,...  =",6i5)') myid,me1a,me1b,me2a,me2b,me3a,me3b
             end if
             if( debug.gt.1 )then
               write(debugFile,'("nn1a,nn1b,...=",6i5)') nn1a,nn1b,nn2a,nn2b,nn3a,nn3b
               write(debugFile,'("mm1a,mm1b,...=",6i5)') mm1a,mm1b,mm2a,mm2b,mm3a,mm3b
             end if
             if( orderOfAccuracy.eq.2 .and. orderOfExtrapolation.lt.3 )then
               write(debugFile,'(" ERROR: interface3d: orderOfExtrapolation<3 ")')
               stop 7716
             end if
             if( orderOfAccuracy.eq.4 .and. orderOfExtrapolation.lt.4 )then
               write(debugFile,'(" ERROR: interface3d: orderOfExtrapolation<4 ")')
               stop 7716
             end if
       rx1=0.
       ry1=0.
       rz1=0.
       if( axis1.eq.0 )then
         rx1=1.
       else if( axis1.eq.1 )then
         ry1=1.
       else 
         rz1=1.
       endif
       rx2=0.
       ry2=0.
       rz2=0.
       if( axis2.eq.0 )then
         rx2=1.
       else if( axis2.eq.1 )then
         ry2=1.
       else 
         rz2=1.
       endif
       ! first time through check that the mask's are consistent
       ! For now we require the masks to both be positive at the same points on the interface
       ! We assign pts where both mask1 and mask2 are discretization pts.
       ! If mask1>0 and mask2<0 then we just leave the extrapolated values in u1 and u2 .  
       if( initialized.eq.0 )then
        if( nd.eq.2 )then
         ! check the consistency of the mask arrays
          i3=n3a
          j3=m3a
          j2=m2a
          do i2=n2a,n2b
           j1=m1a
           do i1=n1a,n1b
           m1 = mask1(i1,i2,i3)
           m2 = mask2(j1,j2,j3)
           if( (m1.gt.0 .and. m2.eq.0) .or. (m1.eq.0 .and. m2.gt.0) )then
             write(debugFile,'(" interface3d:ERROR: mask1 and mask2 do not agree. One is >0 and one =0 ")')
              ! '
             stop 1111
           end if 
            j1=j1+1
           end do
           j2=j2+1
          end do
        else if( nd.eq.3 )then
         ! check the consistency of the mask arrays
          j3=m3a
          do i3=n3a,n3b
          j2=m2a
          do i2=n2a,n2b
          j1=m1a
          do i1=n1a,n1b
           m1 = mask1(i1,i2,i3)
           m2 = mask2(j1,j2,j3)
           if( (m1.gt.0 .and. m2.eq.0) .or. (m1.eq.0 .and. m2.gt.0) )then
             write(debugFile,'(" interface3d:ERROR: mask1 and mask2 do not agree. One is >0 and one =0")')
              ! '
             stop 1111
           end if 
            j1=j1+1
           end do
           j2=j2+1
          end do
           j3=j3+1
          end do
        end if
        if( debug.gt.0 )then
          write(debugFile,'("cgmx:interface3d: The mask arrays for grid1=",i3," and grid2=",i3," were found to be consistent")') grid1,grid2
          ! ' 
        end if
       end if
       if( nd.eq.2 .and. orderOfAccuracy.eq.2 .and. gridType.eq.rectangular )then
         ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         ! End rectangular, 2D, order 2
         ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        else if( nd.eq.2 .and. orderOfAccuracy.eq.2 .and. gridType.eq.curvilinear )then
         ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         ! End curvilinear, 2D, order 2
         ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        else if( nd.eq.2 .and. orderOfAccuracy.eq.4 .and. gridType.eq.rectangular )then
         ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          ! ************************************
          ! ***** 2d rectangular 4th-order *****
          ! ************************************
            ! --------------- 4th Order Rectangular ---------------
            !- ! **TEMP**
            !- alphaP1=0.
            !- alphaP2=0.
            !-write(*,'(" ***TEMP Setting alphaP1=alphaP2=0")') 
            ! if( useForcing.ne.0 )then 
            !   ! finish me 
            !   stop 7716
            ! end if
            ! ***** fix these for [mu] != 0 ****
            if( mu1.ne.mu2 )then
              stop 9924
            end if
            ! ---- first satisfy the jump conditions on the boundary --------
            !    [ eps n.u ] = 0
            !    [ tau.u ] = 0
            if( assignInterfaceValues.eq.1 )then
              if( dispersive.gt.0 )then
               ! ****** DISPERSIVE CASE *******
                if( .true. .and. useImpedanceInterfaceProjection.eq.1 )then
                  ! --------------------------------------------------------------
                  ! ------ Use impedance weighting to project the interface ------
                  ! --------------------------------------------------------------
                  !  (see maxwell.pdf)
                  !
                  if( t.le.3*dt .and. debug.gt.0  )then
                    write(*,'("cgmx:interface3d: PROJECT INTERFACE DISPERSIVE in 2D")')
                  end if
                  ! --- init polarization vectors in case one side is non-dispersive --- 
                  do n=0,nd-1   
                    p1v(n)=0.
                    p2v(n)=0.
                  end do
                  g1=0.
                  g2=0.
                  ! --- LOOP over the interface ---
                   i3=n3a
                   j3=m3a
                   j2=mm2a
                   do i2=nn2a,nn2b
                    j1=mm1a
                    do i1=nn1a,nn1b
                    ! if( mask1(i1,i2,i3).gt.0 .and. mask2(j1,j2,j3).gt.0 )then
                    ! *wdh* 2015/08/14 -- project interpolation points too
                    if( mask1(i1,i2,i3).ne.0 .and. mask2(j1,j2,j3).ne.0 )then
                       an1=an1Cartesian
                       an2=an2Cartesian
                    ! left state:
                    ex1=u1(i1,i2,i3,ex)
                    ey1=u1(i1,i2,i3,ey)
                    hz1=u1(i1,i2,i3,hz)
                    if( dispersionModel1.ne.noDispersion )then
                      ! eval sum of Pv vectors 
                      do n=0,nd-1
                        p1v(n)=0.
                        do jv=0,numberOfPolarizationVectors1-1
                          pc = jv*nd
                          p1v(n)=p1v(n) + p1(i1,i2,i3,pc+n)
                        end do
                      end do
                    end if
                    ! displacement vector 
                    D1v(0) = eps1*(ex1 + alphaP1*p1v(0))
                    D1v(1) = eps1*(ey1 + alphaP1*p1v(1))
                    ! right state:
                    ex2=u2(j1,j2,j3,ex)
                    ey2=u2(j1,j2,j3,ey)
                    hz2=u2(j1,j2,j3,hz)
                    if( dispersionModel2.ne.noDispersion )then
                      ! eval sum of Pv vectors 
                      do n=0,nd-1
                        p2v(n)=0.
                        do jv=0,numberOfPolarizationVectors2-1
                          pc = jv*nd
                          p2v(n)=p2v(n) + p2(j1,j2,j3,pc+n)
                        end do
                      end do
                    end if
                    ! displacement vector 
                    D2v(0) = eps2*(ex2 + alphaP2*p2v(0))
                    D2v(1) = eps2*(ey2 + alphaP2*p2v(1))
                    ! normal components 
                    ! nDotE1 = an1*ex1+an2*ey1
                    ! nDotE2 = an1*ex2+an2*ey2
                    nDotD1 = an1*D1v(0)+an2*D1v(1)
                    nDotD2 = an1*D2v(0)+an2*D2v(1)
                    nDotP1 = an1*p1v(0)+an2*p1v(1)
                    nDotP2 = an1*p2v(0)+an2*p2v(1)
                    ! Interface displacement-vector is an impedance weighted average
                    nDotDI = ( eta1*nDotD1 + eta2*nDotD2 )/(eta1+eta2)
                    ! Normal components of n.E on the left and right:
                    nDotE1I = nDotDI/eps1 -alphaP1*nDotP1
                    nDotE2I = nDotDI/eps2 -alphaP2*nDotP2
                    ! nDotE1I = nDotE1 + (1./(1.*eps1))*(nDotDI-nDotD1)
                    ! nDotE2I = nDotE2 + (1./(1.*eps2))*(nDotDI-nDotD2)
                    if( twilightZone.eq.1)then
                     ! --- adjust for TZ forcing ---
                     ! *NOTE* Here we assume the exact solution for E is the same on both sides
                     ! *NOTE* exact soltuion for P may be different 
                     ! See notes...
                     !    gm = Ee - D_e^I /epsm + alphaP*Pe
                     !    D_e^I = (eta1*D1e + eta2*D2e )/( eta1+eta2 )
                     call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, ue )
                     call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, ve )
                     nDotEe = an1*ue+an2*ve
                     ! Compute sum of exact polarization vectors: 
                     do n=0,nd-1
                      p1ev(n)=0. 
                      do jv=0,numberOfPolarizationVectors1-1
                        ! NOTE: the TZ component number is offset by pxc
                        pc = pxc + jv*nd
                        call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pe(n)  )
                        p1ev(n) = p1ev(n) + pe(n)
                      end do
                      p2ev(n)=0. 
                      do jv=0,numberOfPolarizationVectors2-1
                        ! NOTE: the TZ component number is offset by pxc
                        pc = pxc + jv*nd
                        call ogderiv(ep, 0,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pe(n)  )
                        p2ev(n) = p2ev(n) + pe(n)
                      end do
                     end do
                     ! normal components 
                     nDotP1e = an1*p1ev(0) + an2*p1ev(1)
                     nDotP2e = an1*p2ev(0) + an2*p2ev(1)
                     nDotDe = ( eta1*eps1*(nDotEe+alphaP1*nDotP1e) + eta2*eps2*(nDotEe+alphaP2*nDotP2e) )/(eta1+eta2)
                     ! epsNDotEI + g1 = eps1*nDotE 
                     ! epsNDotEI + g2 = eps2*nDotE 
                     g1= nDotEe - nDotDe/eps1 + alphaP1*nDotP1e 
                     g2= nDotEe - nDotDe/eps2 + alphaP2*nDotP2e 
                     ! Adjust n.E 
                     nDotE1I = nDotE1I+g1  ! nDotE for interface on left 
                     nDotE2I = nDotE2I+g2  ! nDotE for interface on right 
                    end if
                    ! inverse impedance average of tangential components  (set values for full vector and then correct below)
                    exI = ( eta1i*ex1 + eta2i*ex2 )/( eta1i + eta2i )
                    eyI = ( eta1i*ey1 + eta2i*ey2 )/( eta1i + eta2i )
                    ! hz : impedance weighted average 
                    hzI = ( eta1*hz1 + eta2*hz2 )/(eta1+eta2) 
                    nDotEI= an1*exI+an2*eyI  ! we need to subtract off normal component of (exI,eyI) 
                    u1(i1,i2,i3,ex) = exI + (nDotE1I - nDotEI)*an1
                    u1(i1,i2,i3,ey) = eyI + (nDotE1I - nDotEI)*an2
                    u1(i1,i2,i3,hz) = hzI 
                    u2(j1,j2,j3,ex) = exI + (nDotE2I - nDotEI)*an1
                    u2(j1,j2,j3,ey) = eyI + (nDotE2I - nDotEI)*an2
                    u2(j1,j2,j3,hz) = hzI 
                    if( .false. .and. twilightZone.eq.1 )then
                        ! check errors on the boundary
                        do n=0,nd-1
                          call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, es(n)   ) 
                        end do
                        f(0) =  u1(i1,i2,i3,ex) -es(0)
                        f(1) =  u1(i1,i2,i3,ey) -es(1)
                        do n=0,nd-1
                          call ogderiv(ep, 0,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, est(n)   ) 
                        end do
                        f(2) =  u2(j1,j2,j3,ex) -est(0)
                        f(3) =  u2(j1,j2,j3,ey) -est(1)
                        write(*,'(" PROJECT-INTERFACE ERRORS =",4e10.2)') f(0),f(1),f(2),f(3) 
                    end if
                     end if
                     j1=j1+1
                    end do
                    j2=j2+1
                   end do
                else 
                  write(*,'("cgmx:interface3d: FINISH ME -- TWO-D dispersive bndry jump conditions")')
                end if
              else
               ! ****** NON-DISPERSIVE CASE *******
                if( useImpedanceInterfaceProjection.eq.1 )then
                  ! --------------------------------------------------------------
                  ! ------ Use impedance weighting to project the interface ------
                  ! --------------------------------------------------------------
                  !  (see maxwell.pdf)
                  !
                   i3=n3a
                   j3=m3a
                   j2=mm2a
                   do i2=nn2a,nn2b
                    j1=mm1a
                    do i1=nn1a,nn1b
                    ! if( mask1(i1,i2,i3).gt.0 .and. mask2(j1,j2,j3).gt.0 )then
                    ! *wdh* 2015/08/14 -- project interpolation points too
                    if( mask1(i1,i2,i3).ne.0 .and. mask2(j1,j2,j3).ne.0 )then
                    ! eps2 n.u2 = eps1 n.u1
                    !     tau.u2 = tau.u1
                       an1=an1Cartesian
                       an2=an2Cartesian
                    ! left state:
                    ex1=u1(i1,i2,i3,ex)
                    ey1=u1(i1,i2,i3,ey)
                    hz1=u1(i1,i2,i3,hz)
                    ! right state:
                    ex2=u2(j1,j2,j3,ex)
                    ey2=u2(j1,j2,j3,ey)
                    hz2=u2(j1,j2,j3,hz)
                    ! normal components 
                    nDotE1 = an1*ex1+an2*ey1
                    nDotE2 = an1*ex2+an2*ey2
                    ! The interface value of (eps n.E) is an impedance average of (eps* n.E)
                    !
                    !      (eps n.E)_I = [ eta1*( eps1*n.E1 ) +  eta2*( eps2*n.E2 ) ]/[ eta1 + eta2 ]
                    !
                    ! We then set 
                    !   eps1*nDotE1_I = (eps n.E)_I 
                    !   eps2*nDotE2_I = (eps n.E)_I 
                    epsNDotEI = ( eta1*(eps1*nDotE1) + eta2*(eps2*nDotE2) )/(eta1+eta2) ! (eps * n. E)_I 
                    if( twilightZone.eq.1)then
                     ! adjust for TZ forcing (here we assume the exact solution is the same on both sides)
                     call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, ue )
                     call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, ve )
                     nDotEe = an1*ue+an2*ve
                     ! epsNDotEI + g1 = eps1*nDotE 
                     ! epsNDotEI + g2 = eps2*nDotE 
                     g1= (eps1-eps2) * nDotEe *eta2/(eta1+eta2)
                     g2= (eps2-eps1) * nDotEe *eta1/(eta1+eta2)
                    else
                      g1=0.
                      g2=0. 
                    end if
                    nDotE1I = (epsNDotEI+g1)/eps1  ! nDotE for interface on left 
                    nDotE2I = (epsNDotEI+g2)/eps2  ! nDotE for interface on right 
                    ! inverse impedance average of tangential components  (do full vector and correct below)
                    exI = ( eta1i*ex1 + eta2i*ex2 )/( eta1i + eta2i )
                    eyI = ( eta1i*ey1 + eta2i*ey2 )/( eta1i + eta2i )
                    ! hz : impedance weighted average 
                    hzI = ( eta1*hz1 + eta2*hz2 )/(eta1+eta2) 
                    nDotEI= an1*exI+an2*eyI  ! we need to subtract off normal component of (exI,eyI) 
                    u1(i1,i2,i3,ex) = exI + (nDotE1I - nDotEI)*an1
                    u1(i1,i2,i3,ey) = eyI + (nDotE1I - nDotEI)*an2
                    u1(i1,i2,i3,hz) = hzI 
                    u2(j1,j2,j3,ex) = exI + (nDotE2I - nDotEI)*an1
                    u2(j1,j2,j3,ey) = eyI + (nDotE2I - nDotEI)*an2
                    u2(j1,j2,j3,hz) = hzI 
                     end if
                     j1=j1+1
                    end do
                    j2=j2+1
                   end do
                else if( eps1.lt.eps2 )then
                  epsRatio=eps1/eps2
                   i3=n3a
                   j3=m3a
                   j2=mm2a
                   do i2=nn2a,nn2b
                    j1=mm1a
                    do i1=nn1a,nn1b
                    ! if( mask1(i1,i2,i3).gt.0 .and. mask2(j1,j2,j3).gt.0 )then
                    ! *wdh* 2015/08/14 -- project interpolation points too
                    if( mask1(i1,i2,i3).ne.0 .and. mask2(j1,j2,j3).ne.0 )then
                    ! eps2 n.u2 = eps1 n.u1
                    !     tau.u2 = tau.u1
                     an1=an1Cartesian
                     an2=an2Cartesian
                    ua=u1(i1,i2,i3,ex)
                    ub=u1(i1,i2,i3,ey)
                    nDotU = an1*ua+an2*ub
                    if( twilightZone.eq.1 )then
                     ! adjust for TZ forcing (here we assume the exact solution is the same on both sides)
                     call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, ue )
                     call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, ve )
                     nDotU = nDotU - (an1*ue+an2*ve)
                    end if
                    ! u2 equals u1 but with normal component = eps1/eps2*(n.u1)
                    u2(j1,j2,j3,ex) = ua + (nDotU*epsRatio - nDotU)*an1
                    u2(j1,j2,j3,ey) = ub + (nDotU*epsRatio - nDotU)*an2
                    u2(j1,j2,j3,hz) = u1(i1,i2,i3,hz)
                     end if
                     j1=j1+1
                    end do
                    j2=j2+1
                   end do
                else
                  epsRatio=eps2/eps1
                   i3=n3a
                   j3=m3a
                   j2=mm2a
                   do i2=nn2a,nn2b
                    j1=mm1a
                    do i1=nn1a,nn1b
                    ! if( mask1(i1,i2,i3).gt.0 .and. mask2(j1,j2,j3).gt.0 )then
                    ! *wdh* 2015/08/14 -- project interpolation points too
                    if( mask1(i1,i2,i3).ne.0 .and. mask2(j1,j2,j3).ne.0 )then
                    ! eps2 n.u2 = eps1 n.u1
                    !     tau.u2 = tau.u1
                     an1=an1Cartesian
                     an2=an2Cartesian
                    ua=u2(j1,j2,j3,ex)
                    ub=u2(j1,j2,j3,ey)
                    nDotU = an1*ua+an2*ub
                    if( twilightZone.eq.1 )then
                     ! adjust for TZ forcing (here we assume the exact solution is the same on both sides)
                     call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, ue )
                     call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, ve )
              ! write(*,'(" jump: x,y=",2e10.2," ua,ue=",2e10.2," ub,ve=",2e10.2)') xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),ua,ue,ub,ve
                     nDotU = nDotU - (an1*ue+an2*ve)
                    end if
                    u1(i1,i2,i3,ex) = ua + (nDotU*epsRatio - nDotU)*an1
                    u1(i1,i2,i3,ey) = ub + (nDotU*epsRatio - nDotU)*an2
                    u1(i1,i2,i3,hz) = u2(j1,j2,j3,hz)
                     end if
                     j1=j1+1
                    end do
                    j2=j2+1
                   end do
                end if
              end if
            end if
            ! here are the real jump conditions for the ghost points
            ! 0  [ u.x + v.y ] = 0
            ! 1  [ u.xx + u.yy ] = 0
            ! 2  [ v.x - u.y ] =0 
            ! 3  [ (v.xx+v.yy)/eps ] = 0
            ! 4  [ (u.xx+u.yy).x + (v.xx+v.yy).y ] = 0  OR [ (u.xx).x + (v.xx).y ] = 0 OR  [ (u.yy).x + (v.yy).y ] = 0 
            ! 5  [ {(Delta v).x - (Delta u).y}/eps ] =0  -> [ {(v.xxx+v.xyy)-(u.xxy+u.yyy)}/eps ] = 0
            ! 6  [ Delta^2 u/eps ] = 0
            ! 7  [ Delta^2 v/eps^2 ] = 0 
            ! ----- assign ghost using jump conditions -----
            if( assignInterfaceGhostValues.eq.1 )then
              ! initialization step: assign two ghost lines by extrapolation
              if( .true. )then
                i3=n3a
                j3=m3a
                j2=mm2a
                do i2=nn2a,nn2b
                 j1=mm1a
                 do i1=nn1a,nn1b
                ! ---- extrapolate to order 5 -----
                ! extrap to order 5 so exact for degree 4 *wdh* 2015/06/29
                u1(i1-is1,i2-is2,i3,ex)=(5.*u1(i1,i2,i3,ex)-10.*u1(i1+is1,i2+is2,i3+is3,ex)+10.*u1(i1+2*is1,i2+2*is2,i3+2*is3,ex)-5.*u1(i1+3*is1,i2+3*is2,i3+3*is3,ex)+u1(i1+4*is1,i2+4*is2,i3+4*is3,ex))
                u1(i1-is1,i2-is2,i3,ey)=(5.*u1(i1,i2,i3,ey)-10.*u1(i1+is1,i2+is2,i3+is3,ey)+10.*u1(i1+2*is1,i2+2*is2,i3+2*is3,ey)-5.*u1(i1+3*is1,i2+3*is2,i3+3*is3,ey)+u1(i1+4*is1,i2+4*is2,i3+4*is3,ey))
                u1(i1-is1,i2-is2,i3,hz)=(5.*u1(i1,i2,i3,hz)-10.*u1(i1+is1,i2+is2,i3+is3,hz)+10.*u1(i1+2*is1,i2+2*is2,i3+2*is3,hz)-5.*u1(i1+3*is1,i2+3*is2,i3+3*is3,hz)+u1(i1+4*is1,i2+4*is2,i3+4*is3,hz))
                u2(j1-js1,j2-js2,j3,ex)=(5.*u2(j1,j2,j3,ex)-10.*u2(j1+js1,j2+js2,j3+js3,ex)+10.*u2(j1+2*js1,j2+2*js2,j3+2*js3,ex)-5.*u2(j1+3*js1,j2+3*js2,j3+3*js3,ex)+u2(j1+4*js1,j2+4*js2,j3+4*js3,ex))
                u2(j1-js1,j2-js2,j3,ey)=(5.*u2(j1,j2,j3,ey)-10.*u2(j1+js1,j2+js2,j3+js3,ey)+10.*u2(j1+2*js1,j2+2*js2,j3+2*js3,ey)-5.*u2(j1+3*js1,j2+3*js2,j3+3*js3,ey)+u2(j1+4*js1,j2+4*js2,j3+4*js3,ey))
                u2(j1-js1,j2-js2,j3,hz)=(5.*u2(j1,j2,j3,hz)-10.*u2(j1+js1,j2+js2,j3+js3,hz)+10.*u2(j1+2*js1,j2+2*js2,j3+2*js3,hz)-5.*u2(j1+3*js1,j2+3*js2,j3+3*js3,hz)+u2(j1+4*js1,j2+4*js2,j3+4*js3,hz))
                ! --- also extrap 2nd line for now
                u1(i1-2*is1,i2-2*is2,i3,ex)=(5.*u1(i1-is1,i2-is2,i3,ex)-10.*u1(i1-is1+is1,i2-is2+is2,i3+is3,ex)+10.*u1(i1-is1+2*is1,i2-is2+2*is2,i3+2*is3,ex)-5.*u1(i1-is1+3*is1,i2-is2+3*is2,i3+3*is3,ex)+u1(i1-is1+4*is1,i2-is2+4*is2,i3+4*is3,ex))
                u1(i1-2*is1,i2-2*is2,i3,ey)=(5.*u1(i1-is1,i2-is2,i3,ey)-10.*u1(i1-is1+is1,i2-is2+is2,i3+is3,ey)+10.*u1(i1-is1+2*is1,i2-is2+2*is2,i3+2*is3,ey)-5.*u1(i1-is1+3*is1,i2-is2+3*is2,i3+3*is3,ey)+u1(i1-is1+4*is1,i2-is2+4*is2,i3+4*is3,ey))
                u1(i1-2*is1,i2-2*is2,i3,hz)=(5.*u1(i1-is1,i2-is2,i3,hz)-10.*u1(i1-is1+is1,i2-is2+is2,i3+is3,hz)+10.*u1(i1-is1+2*is1,i2-is2+2*is2,i3+2*is3,hz)-5.*u1(i1-is1+3*is1,i2-is2+3*is2,i3+3*is3,hz)+u1(i1-is1+4*is1,i2-is2+4*is2,i3+4*is3,hz))
                u2(j1-2*js1,j2-2*js2,j3,ex)=(5.*u2(j1-js1,j2-js2,j3,ex)-10.*u2(j1-js1+js1,j2-js2+js2,j3+js3,ex)+10.*u2(j1-js1+2*js1,j2-js2+2*js2,j3+2*js3,ex)-5.*u2(j1-js1+3*js1,j2-js2+3*js2,j3+3*js3,ex)+u2(j1-js1+4*js1,j2-js2+4*js2,j3+4*js3,ex))
                u2(j1-2*js1,j2-2*js2,j3,ey)=(5.*u2(j1-js1,j2-js2,j3,ey)-10.*u2(j1-js1+js1,j2-js2+js2,j3+js3,ey)+10.*u2(j1-js1+2*js1,j2-js2+2*js2,j3+2*js3,ey)-5.*u2(j1-js1+3*js1,j2-js2+3*js2,j3+3*js3,ey)+u2(j1-js1+4*js1,j2-js2+4*js2,j3+4*js3,ey))
                u2(j1-2*js1,j2-2*js2,j3,hz)=(5.*u2(j1-js1,j2-js2,j3,hz)-10.*u2(j1-js1+js1,j2-js2+js2,j3+js3,hz)+10.*u2(j1-js1+2*js1,j2-js2+2*js2,j3+2*js3,hz)-5.*u2(j1-js1+3*js1,j2-js2+3*js2,j3+3*js3,hz)+u2(j1-js1+4*js1,j2-js2+4*js2,j3+4*js3,hz))
                if( dispersive.ne.noDispersion )then
                 do jv=0,numberOfPolarizationVectors1-1
                   do n=0,nd-1
                     pc = n + jv*nd 
                     p1(i1  -is1,i2  -is2,i3,pc)=(5.*p1(i1,i2,i3,pc)-10.*p1(i1+is1,i2+is2,i3+is3,pc)+10.*p1(i1+2*is1,i2+2*is2,i3+2*is3,pc)-5.*p1(i1+3*is1,i2+3*is2,i3+3*is3,pc)+p1(i1+4*is1,i2+4*is2,i3+4*is3,pc))
                     p1(i1-2*is1,i2-2*is2,i3,pc)=(5.*p1(i1-is1,i2-is2,i3,pc)-10.*p1(i1-is1+is1,i2-is2+is2,i3+is3,pc)+10.*p1(i1-is1+2*is1,i2-is2+2*is2,i3+2*is3,pc)-5.*p1(i1-is1+3*is1,i2-is2+3*is2,i3+3*is3,pc)+p1(i1-is1+4*is1,i2-is2+4*is2,i3+4*is3,pc))
                     if( .false. .and. twilightZone.eq.1 )then
                       ! *** TEST ****
                       ! call ogderiv(ep, 0,0,0,0, xy1(i1      ,i2      ,i3,0),xy1(i1      ,i2      ,i3,1),0.,t,pc+pxc, p1(i1      ,i2      ,i3,pc)   )
                       call ogderiv(ep, 0,0,0,0, xy1(i1-  is1,i2-  is2,i3,0),xy1(i1-  is1,i2-  is2,i3,1),0.,t,pc+pxc, p1(i1-  is1,i2-  is2,i3,pc)   )
                       call ogderiv(ep, 0,0,0,0, xy1(i1-2*is1,i2-2*is2,i3,0),xy1(i1-2*is1,i2-2*is2,i3,1),0.,t,pc+pxc, p1(i1-2*is1,i2-2*is2,i3,pc)   )
                     end if
                   end do
                 end do
                 do jv=0,numberOfPolarizationVectors2-1
                   do n=0,nd-1
                     pc = n + jv*nd 
                     p2(j1  -js1,j2  -js2,j3,pc)=(5.*p2(j1,j2,j3,pc)-10.*p2(j1+js1,j2+js2,j3+js3,pc)+10.*p2(j1+2*js1,j2+2*js2,j3+2*js3,pc)-5.*p2(j1+3*js1,j2+3*js2,j3+3*js3,pc)+p2(j1+4*js1,j2+4*js2,j3+4*js3,pc))
                     p2(j1-2*js1,j2-2*js2,j3,pc)=(5.*p2(j1-js1,j2-js2,j3,pc)-10.*p2(j1-js1+js1,j2-js2+js2,j3+js3,pc)+10.*p2(j1-js1+2*js1,j2-js2+2*js2,j3+2*js3,pc)-5.*p2(j1-js1+3*js1,j2-js2+3*js2,j3+3*js3,pc)+p2(j1-js1+4*js1,j2-js2+4*js2,j3+4*js3,pc))
                     if( .false. .and. twilightZone.eq.1 )then
                       ! *** TEST ****
                       ! call ogderiv(ep, 0,0,0,0, xy2(j1      ,j2      ,j3,0),xy2(j1      ,j2      ,j3,1),0.,t,pc+pxc, p2(j1      ,j2      ,j3,pc)   )
                       call ogderiv(ep, 0,0,0,0, xy2(j1-  js1,j2-  js2,j3,0),xy2(j1-  js1,j2-  js2,j3,1),0.,t,pc+pxc, p2(j1-  js1,j2-  js2,j3,pc)   )
                       call ogderiv(ep, 0,0,0,0, xy2(j1-2*js1,j2-2*js2,j3,0),xy2(j1-2*js1,j2-2*js2,j3,1),0.,t,pc+pxc, p2(j1-2*js1,j2-2*js2,j3,pc)   )
                     end if
                   end do
                 end do
                end if
                if (nonlinearModel1 .ne. noNonlinearModel) then
                  do na = 0,numberOfAtomicLevels1-1
                    q1(i1  -is1,i2  -is2,i3,na)=(5.*q1(i1,i2,i3,na)-10.*q1(i1+is1,i2+is2,i3+is3,na)+10.*q1(i1+2*is1,i2+2*is2,i3+2*is3,na)-5.*q1(i1+3*is1,i2+3*is2,i3+3*is3,na)+q1(i1+4*is1,i2+4*is2,i3+4*is3,na))
                    q1(i1-2*is1,i2-2*is2,i3,na)=(5.*q1(i1-is1,i2-is2,i3,na)-10.*q1(i1-is1+is1,i2-is2+is2,i3+is3,na)+10.*q1(i1-is1+2*is1,i2-is2+2*is2,i3+2*is3,na)-5.*q1(i1-is1+3*is1,i2-is2+3*is2,i3+3*is3,na)+q1(i1-is1+4*is1,i2-is2+4*is2,i3+4*is3,na))
                  enddo
                end if
                if (nonlinearModel2 .ne. noNonlinearModel) then
                  do na = 0,numberOfAtomicLevels2-1
                    q2(j1  -js1,j2  -js2,j3,na)=(5.*q2(j1,j2,j3,na)-10.*q2(j1+js1,j2+js2,j3+js3,na)+10.*q2(j1+2*js1,j2+2*js2,j3+2*js3,na)-5.*q2(j1+3*js1,j2+3*js2,j3+3*js3,na)+q2(j1+4*js1,j2+4*js2,j3+4*js3,na))
                    q2(j1-2*js1,j2-2*js2,j3,na)=(5.*q2(j1-js1,j2-js2,j3,na)-10.*q2(j1-js1+js1,j2-js2+js2,j3+js3,na)+10.*q2(j1-js1+2*js1,j2-js2+2*js2,j3+2*js3,na)-5.*q2(j1-js1+3*js1,j2-js2+3*js2,j3+3*js3,na)+q2(j1-js1+4*js1,j2-js2+4*js2,j3+4*js3,na))
                  enddo
                endif
                  j1=j1+1
                 end do
                 j2=j2+1
                end do
              else
               ! ---- extrapolate to order 4 -----
                i3=n3a
                j3=m3a
                j2=mm2a
                do i2=nn2a,nn2b
                 j1=mm1a
                 do i1=nn1a,nn1b
                u1(i1-is1,i2-is2,i3,ex)=(4.*u1(i1,i2,i3,ex)-6.*u1(i1+is1,i2+is2,i3+is3,ex)+4.*u1(i1+2*is1,i2+2*is2,i3+2*is3,ex)-u1(i1+3*is1,i2+3*is2,i3+3*is3,ex))
                u1(i1-is1,i2-is2,i3,ey)=(4.*u1(i1,i2,i3,ey)-6.*u1(i1+is1,i2+is2,i3+is3,ey)+4.*u1(i1+2*is1,i2+2*is2,i3+2*is3,ey)-u1(i1+3*is1,i2+3*is2,i3+3*is3,ey))
                u1(i1-is1,i2-is2,i3,hz)=(4.*u1(i1,i2,i3,hz)-6.*u1(i1+is1,i2+is2,i3+is3,hz)+4.*u1(i1+2*is1,i2+2*is2,i3+2*is3,hz)-u1(i1+3*is1,i2+3*is2,i3+3*is3,hz))
                u2(j1-js1,j2-js2,j3,ex)=(4.*u2(j1,j2,j3,ex)-6.*u2(j1+js1,j2+js2,j3+js3,ex)+4.*u2(j1+2*js1,j2+2*js2,j3+2*js3,ex)-u2(j1+3*js1,j2+3*js2,j3+3*js3,ex))
                u2(j1-js1,j2-js2,j3,ey)=(4.*u2(j1,j2,j3,ey)-6.*u2(j1+js1,j2+js2,j3+js3,ey)+4.*u2(j1+2*js1,j2+2*js2,j3+2*js3,ey)-u2(j1+3*js1,j2+3*js2,j3+3*js3,ey))
                u2(j1-js1,j2-js2,j3,hz)=(4.*u2(j1,j2,j3,hz)-6.*u2(j1+js1,j2+js2,j3+js3,hz)+4.*u2(j1+2*js1,j2+2*js2,j3+2*js3,hz)-u2(j1+3*js1,j2+3*js2,j3+3*js3,hz))
                ! --- also extrap 2nd line for now
                u1(i1-2*is1,i2-2*is2,i3,ex)=(4.*u1(i1-is1,i2-is2,i3,ex)-6.*u1(i1-is1+is1,i2-is2+is2,i3+is3,ex)+4.*u1(i1-is1+2*is1,i2-is2+2*is2,i3+2*is3,ex)-u1(i1-is1+3*is1,i2-is2+3*is2,i3+3*is3,ex))
                u1(i1-2*is1,i2-2*is2,i3,ey)=(4.*u1(i1-is1,i2-is2,i3,ey)-6.*u1(i1-is1+is1,i2-is2+is2,i3+is3,ey)+4.*u1(i1-is1+2*is1,i2-is2+2*is2,i3+2*is3,ey)-u1(i1-is1+3*is1,i2-is2+3*is2,i3+3*is3,ey))
                u1(i1-2*is1,i2-2*is2,i3,hz)=(4.*u1(i1-is1,i2-is2,i3,hz)-6.*u1(i1-is1+is1,i2-is2+is2,i3+is3,hz)+4.*u1(i1-is1+2*is1,i2-is2+2*is2,i3+2*is3,hz)-u1(i1-is1+3*is1,i2-is2+3*is2,i3+3*is3,hz))
                u2(j1-2*js1,j2-2*js2,j3,ex)=(4.*u2(j1-js1,j2-js2,j3,ex)-6.*u2(j1-js1+js1,j2-js2+js2,j3+js3,ex)+4.*u2(j1-js1+2*js1,j2-js2+2*js2,j3+2*js3,ex)-u2(j1-js1+3*js1,j2-js2+3*js2,j3+3*js3,ex))
                u2(j1-2*js1,j2-2*js2,j3,ey)=(4.*u2(j1-js1,j2-js2,j3,ey)-6.*u2(j1-js1+js1,j2-js2+js2,j3+js3,ey)+4.*u2(j1-js1+2*js1,j2-js2+2*js2,j3+2*js3,ey)-u2(j1-js1+3*js1,j2-js2+3*js2,j3+3*js3,ey))
                u2(j1-2*js1,j2-2*js2,j3,hz)=(4.*u2(j1-js1,j2-js2,j3,hz)-6.*u2(j1-js1+js1,j2-js2+js2,j3+js3,hz)+4.*u2(j1-js1+2*js1,j2-js2+2*js2,j3+2*js3,hz)-u2(j1-js1+3*js1,j2-js2+3*js2,j3+3*js3,hz))
                if( dispersive.ne.noDispersion )then
                 do jv=0,numberOfPolarizationVectors1-1
                   do n=0,nd-1
                     pc = n + jv*nd 
                     p1(i1  -is1,i2  -is2,i3,pc)=(4.*p1(i1,i2,i3,pc)-6.*p1(i1+is1,i2+is2,i3+is3,pc)+4.*p1(i1+2*is1,i2+2*is2,i3+2*is3,pc)-p1(i1+3*is1,i2+3*is2,i3+3*is3,pc))
                     p1(i1-2*is1,i2-2*is2,i3,pc)=(4.*p1(i1-is1,i2-is2,i3,pc)-6.*p1(i1-is1+is1,i2-is2+is2,i3+is3,pc)+4.*p1(i1-is1+2*is1,i2-is2+2*is2,i3+2*is3,pc)-p1(i1-is1+3*is1,i2-is2+3*is2,i3+3*is3,pc))
                   end do
                 end do
                 do jv=0,numberOfPolarizationVectors2-1
                   do n=0,nd-1
                     pc = n + jv*nd 
                     p2(j1  -js1,j2  -js2,j3,pc)=(4.*p2(j1,j2,j3,pc)-6.*p2(j1+js1,j2+js2,j3+js3,pc)+4.*p2(j1+2*js1,j2+2*js2,j3+2*js3,pc)-p2(j1+3*js1,j2+3*js2,j3+3*js3,pc))
                     p2(j1-2*js1,j2-2*js2,j3,pc)=(4.*p2(j1-js1,j2-js2,j3,pc)-6.*p2(j1-js1+js1,j2-js2+js2,j3+js3,pc)+4.*p2(j1-js1+2*js1,j2-js2+2*js2,j3+2*js3,pc)-p2(j1-js1+3*js1,j2-js2+3*js2,j3+3*js3,pc))
                   end do
                 end do
                end if
                if (nonlinearModel1 .ne. noNonlinearModel) then
                  do na = 0,numberOfAtomicLevels1-1
                    q1(i1  -is1,i2  -is2,i3,na)=(4.*q1(i1,i2,i3,na)-6.*q1(i1+is1,i2+is2,i3+is3,na)+4.*q1(i1+2*is1,i2+2*is2,i3+2*is3,na)-q1(i1+3*is1,i2+3*is2,i3+3*is3,na))
                    q1(i1-2*is1,i2-2*is2,i3,na)=(4.*q1(i1-is1,i2-is2,i3,na)-6.*q1(i1-is1+is1,i2-is2+is2,i3+is3,na)+4.*q1(i1-is1+2*is1,i2-is2+2*is2,i3+2*is3,na)-q1(i1-is1+3*is1,i2-is2+3*is2,i3+3*is3,na))
                  enddo
                end if
                if (nonlinearModel2 .ne. noNonlinearModel) then
                  do na = 0,numberOfAtomicLevels2-1
                    q2(j1  -js1,j2  -js2,j3,na)=(4.*q2(j1,j2,j3,na)-6.*q2(j1+js1,j2+js2,j3+js3,na)+4.*q2(j1+2*js1,j2+2*js2,j3+2*js3,na)-q2(j1+3*js1,j2+3*js2,j3+3*js3,na))
                    q2(j1-2*js1,j2-2*js2,j3,na)=(4.*q2(j1-js1,j2-js2,j3,na)-6.*q2(j1-js1+js1,j2-js2+js2,j3+js3,na)+4.*q2(j1-js1+2*js1,j2-js2+2*js2,j3+2*js3,na)-q2(j1-js1+3*js1,j2-js2+3*js2,j3+3*js3,na))
                  enddo
                endif
                  j1=j1+1
                 end do
                 j2=j2+1
                end do
              end if
              orderOfAccuracy=2 ! temporarily set (for checkCoeff only?)
              if( t.lt. 3.*dt .and. debug.gt.0 )then
                ! This next perl command will cause macro derivatives to be computed to order=2
                 write(*,'(" **** STAGE I: ASSIGN GHOST TO 2ND ORDER FOR 4TH ORDER RECTANGUAR INTERFACE (residuals Order=2)****")')
              end if 
              ! in parallel we add extra points in the tangential direction on parallel boundaries
              ! (otherwise we would use extrapolated values which is probably ok) 
                ! grid1: nn1a,nn1b, etc includes extra ghost in tangential directions
                n1a=ne1a
                n1b=ne1b
                n2a=ne2a
                n2b=ne2b
                n3a=ne3a
                n3b=ne3b
                ! grid2
                m1a=me1a
                m1b=me1b
                m2a=me2a
                m2b=me2b
                m3a=me3a
                m3b=me3b
             ! Macro to assign ghost values:
             if( dispersive.eq.0 )then
             else if( useUnifiedInterfaceMacros.eq.1 ) then
              ! Use Qing's unified interface macros
               ! ****************************************************
               ! ***********  2D, ORDER=2, RECTANGULAR **************
               ! ****************************************************
              if( t.le.3.*dt .and. debug.gt.0 )then
                write(*,'("Interface>>>","22rectangle-unified-dispersive")')
              end if
              ! For rectangular, both sides must axis axis1==axis2: 
              if( axis1.ne.axis2 )then
                stop 8826
              end if
              ! normal and tangent (for TZ forcing)
              an1=an1Cartesian
              an2=an2Cartesian
              tau1=-an2
              tau2= an1
              ! make sure the normal and tangent are set
              if( abs( an1**2 + an2**2 -1. )>1.e-10 .or. abs( tau1**2 + tau2**2 -1. )>1.e-10 )then
                write(*,'("gdm22r - ERROR: incorrect an1,an2, tau1,tau2=",4(1pe9.2))') an1,an2,tau1,tau2
                stop 6666
              end if
              if( abs(abs(an1)-1.) > 1.e-10 )then
                write(*,'("gdm22r - ERROR: only implemented for an1=1, an2=0:  an1,an2, tau1,tau2=",4(1pe9.2))') an1,an2,tau1,tau2
                stop 6667 
              end if 
              ! 
              ! Solve for the unknowns u1(-1),v1(-1),u2(-1),v2(-1)
              !     
              !       A [ U ] = A [ U(old) ] - [ f ]
              !
              !               [ u1(-1) ]
              !       [ U ] = [ v1(-1) ]
              !               [ u2(-1) ]
              !               [ v2(-1) ]
              !             
              ! --- initialize some forcing functions ---
              ! forcing functions for E and P
              do n=0,nd-1
                fev1(n)=0.
                fev2(n)=0.
                if (dispersionModel1 .ne. noDispersion) then
                  do jv=0,numberOfPolarizationVectors1-1
                    fpv1(n,jv)=0.
                  end do
                endif
                if (dispersionModel2 .ne. noDispersion) then
                  do jv=0,numberOfPolarizationVectors2-1
                    fpv2(n,jv)=0.
                  end do
                endif
              end do
              ! forcing functions for N
              if (nonlinearModel1 .ne. noNonlinearModel) then
                do jv = 0,numberOfAtomicLevels1-1
                    fnv1(jv) = 0.
                    fntv1(jv) = 0.
                enddo
              endif
              if (nonlinearModel2 .ne. noNonlinearModel) then
                do jv = 0,numberOfAtomicLevels2-1
                    fnv2(jv) = 0.
                    fntv2(jv) = 0.
                enddo
              endif
              ! print *, "-----------Now using unified (RECTANGULAR)---------------"
              ! ----------------- START LOOP OVER INTERFACE -------------------------
               i3=n3a
               j3=m3a
               j2=m2a
               do i2=n2a,n2b
                j1=m1a
                do i1=n1a,n1b
                 if( mask1(i1,i2,i3).gt.0 .and. mask2(j1,j2,j3).gt.0 )then
                ! first evaluate the equations we want to solve with the wrong values at the ghost points:
                 ! NOTE: the jacobian derivatives can be computed once for all components
                   uu1=u1(i1,i2,i3,ex) ! in the rectangular case just eval the solution
                    u1x = (-u1(i1-1,i2,i3,ex)+u1(i1+1,i2,i3,ex))/(2.*dx1(0))
                    u1y = (-u1(i1,i2-1,i3,ex)+u1(i1,i2+1,i3,ex))/(2.*dx1(1))
                    u1xx = (u1(i1-1,i2,i3,ex)-2.*u1(i1,i2,i3,ex)+u1(i1+1,i2,i3,ex))/(dx1(0)**2)
                    u1yy = (u1(i1,i2-1,i3,ex)-2.*u1(i1,i2,i3,ex)+u1(i1,i2+1,i3,ex))/(dx1(1)**2)
                  u1Lap = u1xx+ u1yy
                   vv1=u1(i1,i2,i3,ey) ! in the rectangular case just eval the solution
                    v1x = (-u1(i1-1,i2,i3,ey)+u1(i1+1,i2,i3,ey))/(2.*dx1(0))
                    v1y = (-u1(i1,i2-1,i3,ey)+u1(i1,i2+1,i3,ey))/(2.*dx1(1))
                    v1xx = (u1(i1-1,i2,i3,ey)-2.*u1(i1,i2,i3,ey)+u1(i1+1,i2,i3,ey))/(dx1(0)**2)
                    v1yy = (u1(i1,i2-1,i3,ey)-2.*u1(i1,i2,i3,ey)+u1(i1,i2+1,i3,ey))/(dx1(1)**2)
                  v1Lap = v1xx+ v1yy
                 ! NOTE: the jacobian derivatives can be computed once for all components
                   uu2=u2(j1,j2,j3,ex) ! in the rectangular case just eval the solution
                    u2x = (-u2(j1-1,j2,j3,ex)+u2(j1+1,j2,j3,ex))/(2.*dx2(0))
                    u2y = (-u2(j1,j2-1,j3,ex)+u2(j1,j2+1,j3,ex))/(2.*dx2(1))
                    u2xx = (u2(j1-1,j2,j3,ex)-2.*u2(j1,j2,j3,ex)+u2(j1+1,j2,j3,ex))/(dx2(0)**2)
                    u2yy = (u2(j1,j2-1,j3,ex)-2.*u2(j1,j2,j3,ex)+u2(j1,j2+1,j3,ex))/(dx2(1)**2)
                  u2Lap = u2xx+ u2yy
                   vv2=u2(j1,j2,j3,ey) ! in the rectangular case just eval the solution
                    v2x = (-u2(j1-1,j2,j3,ey)+u2(j1+1,j2,j3,ey))/(2.*dx2(0))
                    v2y = (-u2(j1,j2-1,j3,ey)+u2(j1,j2+1,j3,ey))/(2.*dx2(1))
                    v2xx = (u2(j1-1,j2,j3,ey)-2.*u2(j1,j2,j3,ey)+u2(j1+1,j2,j3,ey))/(dx2(0)**2)
                    v2yy = (u2(j1,j2-1,j3,ey)-2.*u2(j1,j2,j3,ey)+u2(j1,j2+1,j3,ey))/(dx2(1)**2)
                  v2Lap = v2xx+ v2yy
                ! Evaluate TZ forcing for dispersive equations in 2D 
                  if( twilightZone.eq.1 )then
                      !------------------------
                      ! nonlinear MLA
                      !------------------------
                      if( dispersionModel1.ne.noDispersion .and. nonlinearModel1 .ne. noNonlinearModel) then
                        !-----------------------
                        ! dimension loops for E and P
                        !-----------------------
                        ! t is at new time now
                        nce = pxc+nd*numberOfPolarizationVectors1
                        do n=0,nd-1
                          fpSum1(n)=0.
                          pettSum1(n)=0.
                          call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, es(n)   ) 
                          call ogderiv(ep, 1,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, est(n)  )
                          call ogderiv(ep, 2,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, estt(n) )
                          call ogderiv(ep, 0,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esxx(n) )
                          call ogderiv(ep, 0,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esyy(n) )
                          do jv=0,numberOfPolarizationVectors1-1
                            ! The TZ component is offset by pxc
                            pc = pxc + jv*nd
                            call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pe(n)   )
                            call ogderiv(ep, 1,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pet(n)  )
                            call ogderiv(ep, 2,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pett(n) )
                            ! Normal TZ forcing for P_{n,jv} equation: 
                            ! fpv1(n,jv) = pett(n) + b1v1(jv)*pet(n) + b0v1(jv)*pe(n) - a0v1(jv)*es(n) - a1v1(jv)*est(n)
                            ! left hand side of gdm equations
                            fpv1(n,jv) = pett(n) + b1v1(jv)*pet(n) + b0v1(jv)*pe(n)
                            do na = 0,numberOfAtomicLevels1-1
                              call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, nce+na, q0  )
                              fpv1(n,jv) = fpv1(n,jv) - pnec1(jv,na)*q0*es(n) ! adding \Delta N*E
                              ! print *, 'pnec1',pnec1(jv,na)
                            enddo
                            ! Keep sum: 
                            fpSum1(n)  = fpSum1(n)  + fpv1(n,jv)
                            pettSum1(n) = pettSum1(n) + pett(n) 
                          end do 
                          ! TZ forcing for E_{n} equation:
                          ! E_tt - c1^2 Delta E + alphaP1*Ptt  = 
                          fev1(n) = estt(n) - c1**2*( esxx(n) + esyy(n) ) + alphaP1*pettSum1(n)
                        end do
                        !--------------------------------
                        ! outside of dimension loop for N
                        !--------------------------------
                        do na=0,numberOfAtomicLevels1-1
                          ! na-th level
                          call ogderiv(ep, 1,0,0,0,xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, nce+na, q0t )
                          call ogderiv(ep, 2,0,0,0,xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, nce+na, q0tt)
                          ! initialize
                          fnv1(na)  = q0t ! forcing for \partial_tN_\ell = alpha_{\ell,k}N_k+\beta_{\ell,m}E\cdot\partial_tP_k
                          fntv1(na) = q0tt ! next derivative
                          ! relaxation (alpha_{\ell,m})
                          do jv=0,numberOfAtomicLevels1-1
                            call ogderiv(ep, 0,0,0,0,xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, nce+jv, q0  )
                            call ogderiv(ep, 1,0,0,0,xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, nce+jv, q0t )
                            fnv1(na)  = fnv1(na)  - prc1(na,jv)*q0
                            fntv1(na) = fntv1(na) - prc1(na,jv)*q0t
                          enddo
                          ! dot product (\beta_{\ell,k})
                          do n=0,nd-1 ! loop over dim
                            ! corresponding polarization vector
                            do jv=0,numberOfPolarizationVectors1-1  
                              pc = pxc + jv*nd
                              call ogderiv(ep, 1,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pet(n)  )
                              call ogderiv(ep, 2,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pett(n) )
                              fnv1(na)  = fnv1(na) - peptc1(na,jv)*es(n)*pet(n)
                              fntv1(na) = fntv1(na) - peptc1(na,jv)*est(n)*pet(n) - peptc1(na,jv)*es(n)*pett(n)
                            enddo
                          enddo
                        enddo
                      !------------------------
                      ! GDM
                      !------------------------
                      elseif (dispersionModel1.ne.noDispersion .and. nonlinearModel1.eq.noNonlinearModel) then
                        do n=0,nd-1
                          fpSum1(n)=0.
                          pettSum1(n)=0.
                          call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, es(n)   ) 
                          call ogderiv(ep, 1,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, est(n)  )
                          call ogderiv(ep, 2,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, estt(n) )
                          call ogderiv(ep, 0,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esxx(n) )
                          call ogderiv(ep, 0,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esyy(n) )
                          do jv=0,numberOfPolarizationVectors1-1
                            ! The TZ component is offset by pxc
                            pc = pxc + jv*nd
                            call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pe(n)   )
                            call ogderiv(ep, 1,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pet(n)  )
                            call ogderiv(ep, 2,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pett(n) )
                            ! Normal TZ forcing for P_{n,jv} equation: 
                            fpv1(n,jv) = pett(n) + b1v1(jv)*pet(n) + b0v1(jv)*pe(n) - a0v1(jv)*es(n) - a1v1(jv)*est(n)
                            ! Keep sum: 
                            fpSum1(n)  = fpSum1(n)  + fpv1(n,jv)
                            pettSum1(n) = pettSum1(n) + pett(n) 
                          end do 
                          ! TZ forcing for E_{n} equation:
                          ! E_tt - c1^2 Delta E + alphaP1*Ptt  = 
                          fev1(n) = estt(n) - c1**2*( esxx(n) + esyy(n) ) + alphaP1*pettSum1(n)
                        end do
                      !------------------------
                      ! no dispersion
                      !------------------------
                      elseif (dispersionModel1.eq.noDispersion .and. nonlinearModel1.eq.noNonlinearModel) then
                        do n=0,nd-1
                          fpSum1(n)=0.
                          pettSum1(n)=0.
                          call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, es(n)   ) 
                          call ogderiv(ep, 1,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, est(n)  )
                          call ogderiv(ep, 2,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, estt(n) )
                          call ogderiv(ep, 0,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esxx(n) )
                          call ogderiv(ep, 0,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esyy(n) )
                          fev1(n) = estt(n) - c1**2*( esxx(n) + esyy(n) )
                        end do
                      end if
                      !------------------------
                      ! nonlinear MLA
                      !------------------------
                      if( dispersionModel2.ne.noDispersion .and. nonlinearModel2 .ne. noNonlinearModel) then
                        !-----------------------
                        ! dimension loops for E and P
                        !-----------------------
                        ! t is at new time now
                        nce = pxc+nd*numberOfPolarizationVectors2
                        do n=0,nd-1
                          fpSum2(n)=0.
                          pettSum2(n)=0.
                          call ogderiv(ep, 0,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, es(n)   ) 
                          call ogderiv(ep, 1,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, est(n)  )
                          call ogderiv(ep, 2,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, estt(n) )
                          call ogderiv(ep, 0,2,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esxx(n) )
                          call ogderiv(ep, 0,0,2,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esyy(n) )
                          do jv=0,numberOfPolarizationVectors2-1
                            ! The TZ component is offset by pxc
                            pc = pxc + jv*nd
                            call ogderiv(ep, 0,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pe(n)   )
                            call ogderiv(ep, 1,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pet(n)  )
                            call ogderiv(ep, 2,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pett(n) )
                            ! Normal TZ forcing for P_{n,jv} equation: 
                            ! fpv2(n,jv) = pett(n) + b1v2(jv)*pet(n) + b0v2(jv)*pe(n) - a0v2(jv)*es(n) - a1v2(jv)*est(n)
                            ! left hand side of gdm equations
                            fpv2(n,jv) = pett(n) + b1v2(jv)*pet(n) + b0v2(jv)*pe(n)
                            do na = 0,numberOfAtomicLevels2-1
                              call ogderiv(ep, 0,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t, nce+na, q0  )
                              fpv2(n,jv) = fpv2(n,jv) - pnec2(jv,na)*q0*es(n) ! adding \Delta N*E
                              ! print *, 'pnec2',pnec2(jv,na)
                            enddo
                            ! Keep sum: 
                            fpSum2(n)  = fpSum2(n)  + fpv2(n,jv)
                            pettSum2(n) = pettSum2(n) + pett(n) 
                          end do 
                          ! TZ forcing for E_{n} equation:
                          ! E_tt - c2^2 Delta E + alphaP2*Ptt  = 
                          fev2(n) = estt(n) - c2**2*( esxx(n) + esyy(n) ) + alphaP2*pettSum2(n)
                        end do
                        !--------------------------------
                        ! outside of dimension loop for N
                        !--------------------------------
                        do na=0,numberOfAtomicLevels2-1
                          ! na-th level
                          call ogderiv(ep, 1,0,0,0,xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t, nce+na, q0t )
                          call ogderiv(ep, 2,0,0,0,xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t, nce+na, q0tt)
                          ! initialize
                          fnv2(na)  = q0t ! forcing for \partial_tN_\ell = alpha_{\ell,k}N_k+\beta_{\ell,m}E\cdot\partial_tP_k
                          fntv2(na) = q0tt ! next derivative
                          ! relaxation (alpha_{\ell,m})
                          do jv=0,numberOfAtomicLevels2-1
                            call ogderiv(ep, 0,0,0,0,xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t, nce+jv, q0  )
                            call ogderiv(ep, 1,0,0,0,xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t, nce+jv, q0t )
                            fnv2(na)  = fnv2(na)  - prc2(na,jv)*q0
                            fntv2(na) = fntv2(na) - prc2(na,jv)*q0t
                          enddo
                          ! dot product (\beta_{\ell,k})
                          do n=0,nd-1 ! loop over dim
                            ! corresponding polarization vector
                            do jv=0,numberOfPolarizationVectors2-1  
                              pc = pxc + jv*nd
                              call ogderiv(ep, 1,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pet(n)  )
                              call ogderiv(ep, 2,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pett(n) )
                              fnv2(na)  = fnv2(na) - peptc2(na,jv)*es(n)*pet(n)
                              fntv2(na) = fntv2(na) - peptc2(na,jv)*est(n)*pet(n) - peptc2(na,jv)*es(n)*pett(n)
                            enddo
                          enddo
                        enddo
                      !------------------------
                      ! GDM
                      !------------------------
                      elseif (dispersionModel2.ne.noDispersion .and. nonlinearModel2.eq.noNonlinearModel) then
                        do n=0,nd-1
                          fpSum2(n)=0.
                          pettSum2(n)=0.
                          call ogderiv(ep, 0,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, es(n)   ) 
                          call ogderiv(ep, 1,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, est(n)  )
                          call ogderiv(ep, 2,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, estt(n) )
                          call ogderiv(ep, 0,2,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esxx(n) )
                          call ogderiv(ep, 0,0,2,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esyy(n) )
                          do jv=0,numberOfPolarizationVectors2-1
                            ! The TZ component is offset by pxc
                            pc = pxc + jv*nd
                            call ogderiv(ep, 0,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pe(n)   )
                            call ogderiv(ep, 1,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pet(n)  )
                            call ogderiv(ep, 2,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pett(n) )
                            ! Normal TZ forcing for P_{n,jv} equation: 
                            fpv2(n,jv) = pett(n) + b1v2(jv)*pet(n) + b0v2(jv)*pe(n) - a0v2(jv)*es(n) - a1v2(jv)*est(n)
                            ! Keep sum: 
                            fpSum2(n)  = fpSum2(n)  + fpv2(n,jv)
                            pettSum2(n) = pettSum2(n) + pett(n) 
                          end do 
                          ! TZ forcing for E_{n} equation:
                          ! E_tt - c2^2 Delta E + alphaP2*Ptt  = 
                          fev2(n) = estt(n) - c2**2*( esxx(n) + esyy(n) ) + alphaP2*pettSum2(n)
                        end do
                      !------------------------
                      ! no dispersion
                      !------------------------
                      elseif (dispersionModel2.eq.noDispersion .and. nonlinearModel2.eq.noNonlinearModel) then
                        do n=0,nd-1
                          fpSum2(n)=0.
                          pettSum2(n)=0.
                          call ogderiv(ep, 0,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, es(n)   ) 
                          call ogderiv(ep, 1,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, est(n)  )
                          call ogderiv(ep, 2,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, estt(n) )
                          call ogderiv(ep, 0,2,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esxx(n) )
                          call ogderiv(ep, 0,0,2,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esyy(n) )
                          fev2(n) = estt(n) - c2**2*( esxx(n) + esyy(n) )
                        end do
                      end if
                  end if
                ! eval dispersive forcings for domain 1
                  ! no dispersion
                  do n=0,nd-1
                    fp1(n)=0.
                  end do
                  !------------------------
                  ! nonlinear MLA
                  !------------------------
                  if( dispersionModel1.ne.noDispersion .and. nonlinearModel1.ne.noNonlinearModel) then
                    nce = pxc+nd*numberOfPolarizationVectors1
                    do n=0,nd-1
                      do jv=0,numberOfPolarizationVectors1-1
                        pc = n + jv*nd 
                        ec = ex +n
                        ! in time order: p1m,p1n,p1
                        ! pvm   =  p1(i1,i2,i3,pc) ! this one needs to be replaced due to rank difference
                        ! pv  =  p1n(i1,i2,i3,pc)
                        ! evm    =  u1(i1,i2,i3,ec)
                        ! ev   =  u1n(i1,i2,i3,ec)
                        pvn = 2.*p1(i1,i2,i3,pc)-p1n(i1,i2,i3,pc) + 0.5*dt*b1v1(jv)*p1n(i1,i2,i3,pc) - dtsq*b0v1(jv)*p1(i1,i2,i3,pc) + dtsq*fpv1(n,jv)
                        do na = 0,numberOfAtomicLevels1-1 ! \Delta N^n*E^n
                          pvn = pvn + dtsq*pnec1(jv,na)*q1(i1,i2,i3,na)*u1(i1,i2,i3,ec)
                        enddo
                        pvn= pvn/( 1.+.5*dt*b1v1(jv) )
                        ! #If "p1" eq "p1"
                        ! call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t+dt,pxc+jv*nd+n, pe(n)   )
                        ! #Else
                        ! call ogderiv(ep, 0,0,0,0, xy2(i1,i2,i3,0),xy2(i1,i2,i3,1),0.,t+dt,pxc+jv*nd+n, pe(n)   )
                        ! #End
                        ! print *, '+++++++pvn',pvn,'exact',pe(n),'diff',pvn-pe(n)
                        ! marching from t^{n+1} to t^{n+2}
                        ! #If "p1" eq "p1"
                        ! call ogderiv(ep, 2,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pett(n)   )
                        ! #Else
                        ! call ogderiv(ep, 2,0,0,0, xy2(i1,i2,i3,0),xy2(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pett(n)   )
                        ! #End
                        fp1(n) = fp1(n) + (pvn-2.*p1(i1,i2,i3,pc)+p1n(i1,i2,i3,pc))/dtsq
                        ! fp1(n) = fp1(n) + pett(n)
                        ! print *, '------pvtt',(pvn-2.*p1(i1,i2,i3,pc)+p1n(i1,i2,i3,pc))/dtsq,'exact',pett(n),'diff',(pvn-2.*p1(i1,i2,i3,pc)+p1n(i1,i2,i3,pc))/dtsq-pett(n)
                        ! print*,'dtsq',dtsq,'diff1',pvn-2.*p1(i1,i2,i3,pc)+p1n(i1,i2,i3,pc),'diff2',p1(i1,i2,i3,pc)-2.*p1n(i1,i2,i3,pc)+p1m(i1,i2,i3,pc)
                        ! write(*,'(" numberOfAtomicLevels1=",i3)') numberOfAtomicLevels1
                        ! write(*,'(" i1,i2,i3=",3i3)') i1,i2,i3
                        ! na=numberOfAtomicLevels1-1
                        ! write(*,'(" q1,q1n,q1m=",3e12.2)') q1(i1,i2,i3,na),q1n(i1,i2,i3,na),q1m(i1,i2,i3,na)
                        ! write(*,'(" pc=",i3," pvn,p1,p1n,p1m=",4e12.2)') pc, pvn,p1(i1,i2,i3,pc),p1n(i1,i2,i3,pc),p1m(i1,i2,i3,pc)
                        ! ! write(*,'(" dt=",e12.2," pv,pvt,pvtt, ev,evt,evtt=",6e12.2)') dt,pv,pvt,pvtt, ev,evt,evtt
                        ! ! write(*,'(" jv=",i2," a0,a1,b0,b1=",4e12.2," Bk,Ck=",2e12.2)') jv,a0v1(jv),a1v1(jv),b0v1(jv),b1v1(jv),Bk,Ck
                        ! write(*,'(" n=",i2," fev1(n)=",e12.2," fp1(n)=",e12.2," fpv1(n,jv)=",e12.2)') n,fev1(n),fp1(n),fpv1(n,jv)
                      end do
                    end do
                    beta1 = 1.
                  !-------------------
                  ! GDM
                  !-------------------
                  elseif (dispersionModel1.ne.noDispersion .and. nonlinearModel1.eq.noNonlinearModel) then
                    Csum=0.
                    do jv=0,numberOfPolarizationVectors1-1
                      a0=a0v1(jv)
                      a1=a1v1(jv)
                      b0=b0v1(jv)
                      b1=b1v1(jv)
                      alpha=alphaP1
                      !  Second-order coefficients: 
                      !  Ptt = c2PttLE*LE + c2PttE*E + c2PttEm*Em + c2PttP*P + c2PttPm*Pm + c2PttfE*fE + c2PttfP*fP
      ! File created by Dropbox/GDM/maple/interface.maple
      d1=2.+(a1*alpha+b1)*dt

      ! --------------- Here is Et to second-order ------------

      ! Et = -1/2/dt*(2*a0*alpha*dt^2-2*b1*dt-4)/d1*E-1/2/dt*(2*b1*dt+4)/d1*Em-1/2/dt*(-b1*dt^3-2*dt^2)/d1*LE-1/2/dt*(-2*alpha*b0*dt^2-2*alpha*b1*dt)/d1*P-alpha*b1/d1*Pm-1/2/dt*(-b1*dt^3-2*dt^2)/d1*fE-dt*alpha*fP/d1
      ! Et = c2EtLE*LE + c2EtE*E + c2EtEm*Em + c2EtP*P + c2EtPm*Pm + c2EtfE*fE + c2EtfP*fP
      !   LE = c^2Delta(E) 
      c2EtLE=(.5*b1*dt**2+1.0*dt)/d1
      c2EtE=(-a0*alpha*dt**2+b1*dt+2.)/dt/d1
      c2EtEm=(-b1*dt-2.)/dt/d1
      c2EtP=alpha*(b0*dt+b1)/d1
      c2EtPm=-1.*alpha*b1/d1
      c2EtfE=(.5*b1*dt**2+1.0*dt)/d1
      c2EtfP=-1.*dt*alpha/d1
      ! --------------- Here is Pt to second-order ------------

      ! Pt = 1/2*(2*a0*dt^2+2*a1*dt)/dt/d1*E-a1/d1*Em+1/2*a1*dt^2/d1*LE+1/2*(2*a1*alpha*dt-2*b0*dt^2+4)/dt/d1*P+1/2*(-2*a1*alpha*dt-4)/dt/d1*Pm+1/2*a1*dt^2/d1*fE+fP*dt/d1
      ! Pt = c2PtLE*LE + c2PtE*E + c2PtEm*Em + c2PtP*P + c2PtPm*Pm + c2PtfE*fE + c2PtfP*fP
      !   LE = c^2Delta(E) 
      c2PtLE=.5*a1*dt**2/d1
      c2PtE=(1.*a0*dt+1.*a1)/d1
      c2PtEm=-1.*a1/d1
      c2PtP=(alpha*a1*dt-b0*dt**2+2.)/dt/d1
      c2PtPm=(-alpha*a1*dt-2.)/dt/d1
      c2PtfE=.5*a1*dt**2/d1
      c2PtfP=1.*dt/d1
      ! --------------- Here is Ptt to second-order ------------

      ! Ptt = a1*dt/d1*LE+((-a0*a1*alpha-a0*b1)*dt^2+d1*a0*dt+2*a1)/dt/d1*E-2*a1/dt/d1*Em+((a1*alpha*b0+b0*b1)*dt^2-d1*b0*dt-2*b1)/dt/d1*P+2*b1/dt/d1*Pm+a1*dt/d1*fE+((-a1*alpha-b1)*dt^2+dt*d1)/dt/d1*fP
      ! Ptt = c2PttLE*LE + c2PttE*E + c2PttEm*Em + c2PttP*P + c2PttPm*Pm + c2PttfE*fE + c2PttfP*fP
      c2PttLE=1.*a1*dt/d1
      c2PttE=((-1.*a1*alpha-1.*b1)*a0*dt**2+d1*a0*dt+2.*a1)/dt/d1
      c2PttEm=-2.*a1/dt/d1
      c2PttP=(b0*(a1*alpha+b1)*dt**2-1.*d1*b0*dt-2.*b1)/dt/d1
      c2PttPm=2.*b1/dt/d1
      c2PttfE=1.*a1*dt/d1
      c2PttfP=(1.*d1+(-1.*a1*alpha-1.*b1)*dt)/d1
                      Csum = Csum + c2PttLE
                      ! Bk = 1 + .5*dt*( b1v1(jv) + alphaP1*a1v1(jv) )
                      ! Ck = (1./Bk)*a1v1(jv)*dt*.5
                      ! Csum = Csum + Ck 
                      do n=0,nd-1
                        pc = n + jv*nd 
                        ec = ex +n
                        ! P at new time t+dt
                        ! Pt, Ptt at time t
                        pv   =  p1(i1,i2,i3,pc)
                        pvn  =  p1n(i1,i2,i3,pc)
                        ev    =  u1(i1,i2,i3,ec)
                        evn   =  u1n(i1,i2,i3,ec)
                        ! Ptt = c2PttLE*LE + c2PttE*E + c2PttEm*Em + c2PttP*P + c2PttPm*Pm + c2PttfE*fE + c2PttfP*fP
                        ! Levae off term: c2PttLE*LE(n)
                        fp1(n) = fp1(n) + c2PttE*ev + c2PttEm*evn + c2PttP*pv + c2PttPm*pvn + c2PttfE*fev1(n) + c2PttfP*fpv1(n,jv)
                        ! write(*,'(" i1,i2,i3=",3i3)') i1,i2,i3
                        ! write(*,'(" pc=",i3," p1,p1n,p1m=",3e12.2)') pc, p1(i1,i2,i3,pc),p1n(i1,i2,i3,pc),p1m(i1,i2,i3,pc)
                        ! write(*,'(" dt=",e12.2," pv,pvt,pvtt, ev,evt,evtt=",6e12.2)') dt,pv,pvt,pvtt, ev,evt,evtt
                        ! write(*,'(" jv=",i2," a0,a1,b0,b1=",4e12.2," Bk,Ck=",2e12.2)') jv,a0v1(jv),a1v1(jv),b0v1(jv),b1v1(jv),Bk,Ck
                        ! write(*,'(" n=",i2," fev1(n)=",e12.2," fp1(n)=",e12.2," fpv1(n,jv)=",e12.2)') n,fev1(n),fp1(n),fpv1(n,jv)
                      end do
                    end do
                    ! we could precompute D
                    beta1 = 1. -alphaP1*Csum
                  else
                    beta1 = 1.
                  end if
                ! eval dispersive forcings for domain 2
                  ! no dispersion
                  do n=0,nd-1
                    fp2(n)=0.
                  end do
                  !------------------------
                  ! nonlinear MLA
                  !------------------------
                  if( dispersionModel2.ne.noDispersion .and. nonlinearModel2.ne.noNonlinearModel) then
                    nce = pxc+nd*numberOfPolarizationVectors2
                    do n=0,nd-1
                      do jv=0,numberOfPolarizationVectors2-1
                        pc = n + jv*nd 
                        ec = ex +n
                        ! in time order: p2m,p2n,p2
                        ! pvm   =  p2(j1,j2,j3,pc) ! this one needs to be replaced due to rank difference
                        ! pv  =  p2n(j1,j2,j3,pc)
                        ! evm    =  u2(j1,j2,j3,ec)
                        ! ev   =  u2n(j1,j2,j3,ec)
                        pvn = 2.*p2(j1,j2,j3,pc)-p2n(j1,j2,j3,pc) + 0.5*dt*b1v2(jv)*p2n(j1,j2,j3,pc) - dtsq*b0v2(jv)*p2(j1,j2,j3,pc) + dtsq*fpv2(n,jv)
                        do na = 0,numberOfAtomicLevels2-1 ! \Delta N^n*E^n
                          pvn = pvn + dtsq*pnec2(jv,na)*q2(j1,j2,j3,na)*u2(j1,j2,j3,ec)
                        enddo
                        pvn= pvn/( 1.+.5*dt*b1v2(jv) )
                        ! #If "p2" eq "p1"
                        ! call ogderiv(ep, 0,0,0,0, xy1(j1,j2,j3,0),xy1(j1,j2,j3,1),0.,t+dt,pxc+jv*nd+n, pe(n)   )
                        ! #Else
                        ! call ogderiv(ep, 0,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t+dt,pxc+jv*nd+n, pe(n)   )
                        ! #End
                        ! print *, '+++++++pvn',pvn,'exact',pe(n),'diff',pvn-pe(n)
                        ! marching from t^{n+1} to t^{n+2}
                        ! #If "p2" eq "p1"
                        ! call ogderiv(ep, 2,0,0,0, xy1(j1,j2,j3,0),xy1(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pett(n)   )
                        ! #Else
                        ! call ogderiv(ep, 2,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pett(n)   )
                        ! #End
                        fp2(n) = fp2(n) + (pvn-2.*p2(j1,j2,j3,pc)+p2n(j1,j2,j3,pc))/dtsq
                        ! fp2(n) = fp2(n) + pett(n)
                        ! print *, '------pvtt',(pvn-2.*p2(j1,j2,j3,pc)+p2n(j1,j2,j3,pc))/dtsq,'exact',pett(n),'diff',(pvn-2.*p2(j1,j2,j3,pc)+p2n(j1,j2,j3,pc))/dtsq-pett(n)
                        ! print*,'dtsq',dtsq,'diff1',pvn-2.*p2(j1,j2,j3,pc)+p2n(j1,j2,j3,pc),'diff2',p2(j1,j2,j3,pc)-2.*p2n(j1,j2,j3,pc)+p2m(j1,j2,j3,pc)
                        ! write(*,'(" numberOfAtomicLevels2=",i3)') numberOfAtomicLevels2
                        ! write(*,'(" j1,j2,j3=",3i3)') j1,j2,j3
                        ! na=numberOfAtomicLevels2-1
                        ! write(*,'(" q2,q2n,q2m=",3e12.2)') q2(j1,j2,j3,na),q2n(j1,j2,j3,na),q2m(j1,j2,j3,na)
                        ! write(*,'(" pc=",i3," pvn,p2,p2n,p2m=",4e12.2)') pc, pvn,p2(j1,j2,j3,pc),p2n(j1,j2,j3,pc),p2m(j1,j2,j3,pc)
                        ! ! write(*,'(" dt=",e12.2," pv,pvt,pvtt, ev,evt,evtt=",6e12.2)') dt,pv,pvt,pvtt, ev,evt,evtt
                        ! ! write(*,'(" jv=",i2," a0,a1,b0,b1=",4e12.2," Bk,Ck=",2e12.2)') jv,a0v2(jv),a1v2(jv),b0v2(jv),b1v2(jv),Bk,Ck
                        ! write(*,'(" n=",i2," fev2(n)=",e12.2," fp2(n)=",e12.2," fpv2(n,jv)=",e12.2)') n,fev2(n),fp2(n),fpv2(n,jv)
                      end do
                    end do
                    beta2 = 1.
                  !-------------------
                  ! GDM
                  !-------------------
                  elseif (dispersionModel2.ne.noDispersion .and. nonlinearModel2.eq.noNonlinearModel) then
                    Csum=0.
                    do jv=0,numberOfPolarizationVectors2-1
                      a0=a0v2(jv)
                      a1=a1v2(jv)
                      b0=b0v2(jv)
                      b1=b1v2(jv)
                      alpha=alphaP2
                      !  Second-order coefficients: 
                      !  Ptt = c2PttLE*LE + c2PttE*E + c2PttEm*Em + c2PttP*P + c2PttPm*Pm + c2PttfE*fE + c2PttfP*fP
      ! File created by Dropbox/GDM/maple/interface.maple
      d1=2.+(a1*alpha+b1)*dt

      ! --------------- Here is Et to second-order ------------

      ! Et = -1/2/dt*(2*a0*alpha*dt^2-2*b1*dt-4)/d1*E-1/2/dt*(2*b1*dt+4)/d1*Em-1/2/dt*(-b1*dt^3-2*dt^2)/d1*LE-1/2/dt*(-2*alpha*b0*dt^2-2*alpha*b1*dt)/d1*P-alpha*b1/d1*Pm-1/2/dt*(-b1*dt^3-2*dt^2)/d1*fE-dt*alpha*fP/d1
      ! Et = c2EtLE*LE + c2EtE*E + c2EtEm*Em + c2EtP*P + c2EtPm*Pm + c2EtfE*fE + c2EtfP*fP
      !   LE = c^2Delta(E) 
      c2EtLE=(.5*b1*dt**2+1.0*dt)/d1
      c2EtE=(-a0*alpha*dt**2+b1*dt+2.)/dt/d1
      c2EtEm=(-b1*dt-2.)/dt/d1
      c2EtP=alpha*(b0*dt+b1)/d1
      c2EtPm=-1.*alpha*b1/d1
      c2EtfE=(.5*b1*dt**2+1.0*dt)/d1
      c2EtfP=-1.*dt*alpha/d1
      ! --------------- Here is Pt to second-order ------------

      ! Pt = 1/2*(2*a0*dt^2+2*a1*dt)/dt/d1*E-a1/d1*Em+1/2*a1*dt^2/d1*LE+1/2*(2*a1*alpha*dt-2*b0*dt^2+4)/dt/d1*P+1/2*(-2*a1*alpha*dt-4)/dt/d1*Pm+1/2*a1*dt^2/d1*fE+fP*dt/d1
      ! Pt = c2PtLE*LE + c2PtE*E + c2PtEm*Em + c2PtP*P + c2PtPm*Pm + c2PtfE*fE + c2PtfP*fP
      !   LE = c^2Delta(E) 
      c2PtLE=.5*a1*dt**2/d1
      c2PtE=(1.*a0*dt+1.*a1)/d1
      c2PtEm=-1.*a1/d1
      c2PtP=(alpha*a1*dt-b0*dt**2+2.)/dt/d1
      c2PtPm=(-alpha*a1*dt-2.)/dt/d1
      c2PtfE=.5*a1*dt**2/d1
      c2PtfP=1.*dt/d1
      ! --------------- Here is Ptt to second-order ------------

      ! Ptt = a1*dt/d1*LE+((-a0*a1*alpha-a0*b1)*dt^2+d1*a0*dt+2*a1)/dt/d1*E-2*a1/dt/d1*Em+((a1*alpha*b0+b0*b1)*dt^2-d1*b0*dt-2*b1)/dt/d1*P+2*b1/dt/d1*Pm+a1*dt/d1*fE+((-a1*alpha-b1)*dt^2+dt*d1)/dt/d1*fP
      ! Ptt = c2PttLE*LE + c2PttE*E + c2PttEm*Em + c2PttP*P + c2PttPm*Pm + c2PttfE*fE + c2PttfP*fP
      c2PttLE=1.*a1*dt/d1
      c2PttE=((-1.*a1*alpha-1.*b1)*a0*dt**2+d1*a0*dt+2.*a1)/dt/d1
      c2PttEm=-2.*a1/dt/d1
      c2PttP=(b0*(a1*alpha+b1)*dt**2-1.*d1*b0*dt-2.*b1)/dt/d1
      c2PttPm=2.*b1/dt/d1
      c2PttfE=1.*a1*dt/d1
      c2PttfP=(1.*d1+(-1.*a1*alpha-1.*b1)*dt)/d1
                      Csum = Csum + c2PttLE
                      ! Bk = 1 + .5*dt*( b1v2(jv) + alphaP2*a1v2(jv) )
                      ! Ck = (1./Bk)*a1v2(jv)*dt*.5
                      ! Csum = Csum + Ck 
                      do n=0,nd-1
                        pc = n + jv*nd 
                        ec = ex +n
                        ! P at new time t+dt
                        ! Pt, Ptt at time t
                        pv   =  p2(j1,j2,j3,pc)
                        pvn  =  p2n(j1,j2,j3,pc)
                        ev    =  u2(j1,j2,j3,ec)
                        evn   =  u2n(j1,j2,j3,ec)
                        ! Ptt = c2PttLE*LE + c2PttE*E + c2PttEm*Em + c2PttP*P + c2PttPm*Pm + c2PttfE*fE + c2PttfP*fP
                        ! Levae off term: c2PttLE*LE(n)
                        fp2(n) = fp2(n) + c2PttE*ev + c2PttEm*evn + c2PttP*pv + c2PttPm*pvn + c2PttfE*fev2(n) + c2PttfP*fpv2(n,jv)
                        ! write(*,'(" j1,j2,j3=",3i3)') j1,j2,j3
                        ! write(*,'(" pc=",i3," p2,p2n,p2m=",3e12.2)') pc, p2(j1,j2,j3,pc),p2n(j1,j2,j3,pc),p2m(j1,j2,j3,pc)
                        ! write(*,'(" dt=",e12.2," pv,pvt,pvtt, ev,evt,evtt=",6e12.2)') dt,pv,pvt,pvtt, ev,evt,evtt
                        ! write(*,'(" jv=",i2," a0,a1,b0,b1=",4e12.2," Bk,Ck=",2e12.2)') jv,a0v2(jv),a1v2(jv),b0v2(jv),b1v2(jv),Bk,Ck
                        ! write(*,'(" n=",i2," fev2(n)=",e12.2," fp2(n)=",e12.2," fpv2(n,jv)=",e12.2)') n,fev2(n),fp2(n),fpv2(n,jv)
                      end do
                    end do
                    ! we could precompute D
                    beta2 = 1. -alphaP2*Csum
                  else
                    beta2 = 1.
                  end if
                 f(0)=(u1x+v1y) - (u2x+v2y)
                 f(1)=( an1*u1Lap +an2*v1Lap )/mu1 - ( an1*u2Lap +an2*v2Lap )/mu2 
                 f(2)=(v1x-u1y)/mu1 - (v2x-u2y)/mu2
                 f(3)=( ( tau1*u1Lap +tau2*v1Lap )*beta1/epsmu1 - alphaP1*(tau1*fp1(0)+tau2*fp1(1)) ) - ( ( tau1*u2Lap +tau2*v2Lap )*beta2/epsmu2 - alphaP2*(tau1*fp2(0)+tau2*fp2(1)) )
                 if( twilightZone.eq.1 )then
                    ! For now we assume mu1=mu2 and TZ solutions are the same on both sides.
                    ! call ogderiv(ep, 0,1,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, uex  )
                    call ogderiv(ep, 0,0,1,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, uey  )
                    call ogderiv(ep, 0,1,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, vex  )
                    ! call ogderiv(ep, 0,0,1,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, vey  )
                    call ogderiv(ep, 0,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, uexx )
                    call ogderiv(ep, 0,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, ueyy )
                    call ogderiv(ep, 0,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, vexx )
                    call ogderiv(ep, 0,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, veyy )
                    ueLap = uexx + ueyy
                    veLap = vexx + veyy
                    f(1) = f(1) - ( ueLap      )*(1./mu1 - 1./mu2)
                    f(2) = f(2) - ( vex - uey  )*(1./mu1 - 1./mu2)
                    ! f(3) = [ tv.E.tt] = [ tv.( c^2*Delta(E) - alphaP*P.tt + fev) ] 
                    !      = [ tv.( c^2*Delta(E) - alphaP*P.tt] + [ tv.fev ]
                    ! -- add on the jump in the forcing ---
                    f(3) = f(3) + ( tau1*(fev1(0)-fev2(0)) + tau2*(fev1(1)-fev2(1)) )
                   ! For now we assume mu1=mu2 and TZ solutions are the same on both sides.
                   ! f(3) = [ tv.E.tt] = [ tv.( c^2*Delta(E) - alphaP*P.tt + fev) ] 
                   !      = [ tv.( c^2*Delta(E) - alphaP*P.tt] + [ tv.fev ]
                   ! -- add on the jump in the forcing ---
                   ! f(3) = f(3) + ( tau1*(fev1(0)-fev2(0)) + tau2*(fev1(1)-fev2(1)) )
                  !-    ! f(3) = [ tv.( c^2*Delta(E) - alphaP*P.tt ] - [ tv.( c^2*Delta(E^e) - alphaP*(P^e).tt ] =0 
                  !     call ogderiv(ep, 0,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, uexx )
                  !     call ogderiv(ep, 0,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, ueyy )
                  !     call ogderiv(ep, 0,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, vexx )
                  !     call ogderiv(ep, 0,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, veyy )
                  !     ueLap = uexx + ueyy
                  !     veLap = vexx + veyy
                  !     f(3) = f(3) - ( tau1*ueLap +tau2*veLap )*(beta1/epsmu1-beta2/epsmu2)
                  !     f(3) = f(3) +  alphaP1*( tau1*pettSum1(0) + tau2*pettSum1(1) ) - !                    alphaP2*( tau1*pettSum2(0) + tau2*pettSum2(1) )
                  ! !-    end if 
                   ! write(debugFile,'(" u1Lap,ueLap=",2e10.2," v1Lap,veLap=",2e10.2)') u1Lap,ueLap,v1Lap,veLap
                   ! write(debugFile,'(" u1Lap,ueLap=",2e10.2," v1Lap,veLap=",2e10.2)') u1Lap,ueLap,v1Lap,veLap
                 end if
                if( axis1.eq.0 .and. axis2.eq.0 )then
                  ! Interface equations for a boundary at x = 0 or x=1
                  ! ---- EQUATION 0 -----
                  a4(0,0) = -is1/(2.*dx1(axis1))    ! coeff of u1(-1) from [u.x+v.y] 
                  a4(0,1) = 0.                      ! coeff of v1(-1) from [u.x+v.y] 
                  a4(0,2) =  js1/(2.*dx2(axis2))    ! coeff of u2(-1) from [u.x+v.y] 
                  a4(0,3) = 0.                      ! coeff of v2(-1) from [u.x+v.y]
                  ! ---- EQUATION 1 -----    
                  a4(1,0) = 1./(dx1(axis1)**2)/mu1   ! coeff of u1(-1) from [(u.xx + u.yy)/mu]
                  a4(1,1) = 0. 
                  a4(1,2) =-1./(dx2(axis2)**2)/mu2   ! coeff of u2(-1) from [(u.xx + u.yy)/mu]
                  a4(1,3) = 0. 
                  ! ---- EQUATION 2 -----
                  a4(2,0) = 0.
                  a4(2,1) = -is1/(2.*dx1(axis1))/mu1    ! coeff of v1(-1) from [(v.x - u.y)/mu] 
                  a4(2,2) = 0.
                  a4(2,3) =  js1/(2.*dx2(axis2))/mu2    ! coeff of v2(-1) from [(v.x - u.y)/mu]
                  ! ---- EQUATION 3 -----    
                  ! The coefficient of Delta(E) in this equation is altered due to Ptt term 
                  ! write(*,'(" beta1,beta2=",2e10.2," fp1,fp2=",2e10.2)') beta1,beta2,fp1(1),fp2(1)
                  a4(3,0) = 0.                      
                  a4(3,1) = (beta1/epsmu1)/(dx1(axis1)**2) ! coeff of v1(-1) from [beta*c^2*(v.xx+v.yy)]
                  a4(3,2) = 0. 
                  a4(3,3) =-(beta2/epsmu2)/(dx2(axis2)**2) ! coeff of v2(-1) from [beta*c^2*(v.xx+v.yy)]
                else
                  write(*,*) 'gdm22r finish me for axis1.ne.0 or axis2.ne.0 ...'
                  stop 767
                end if
                 if( debug>7 .and. i2.le.0 )then
                   write(*,*) "gdm22r: Matrix a4"
                   do n=0,3
                     write(*,'(4(1pe10.2))') (a4(n,nn),nn=0,3)
                   end do 
                 end if 
                 q(0) = u1(i1-is1,i2-is2,i3,ex)
                 q(1) = u1(i1-is1,i2-is2,i3,ey)
                 q(2) = u2(j1-js1,j2-js2,j3,ex)
                 q(3) = u2(j1-js1,j2-js2,j3,ey)
                 ! subtract off the contributions from the wrong values at the ghost points:
                 do n=0,3
                   f(n) = (a4(n,0)*q(0)+a4(n,1)*q(1)+a4(n,2)*q(2)+a4(n,3)*q(3)) - f(n)
                 end do
                 if( debug>0 .and. twilightZone.eq.1 )then
                    write(*,'(" unified22r --> i1,i2=",2i4," RHS f=",4f12.7)') i1,i2,f(0),f(1),f(2),f(3)
                    write(*,'("                              q=",4f12.7)') q(0),q(1),q(2),q(3)
                  end if
                 ! solve A Q = F
                 ! factor the matrix
                 ! numberOfEquations=4
                 call dgeco( a4(0,0), 4, 4, ipvt(0),rcond,work(0))
                 ! solve
                 ! write(debugFile,'(" --> i1,i2=",2i4," rcond=",e10.2)') i1,i2,rcond
                 job=0
                 call dgesl( a4(0,0), 4, 4, ipvt(0), f(0), job)
                 ! write(debugFile,'(" --> i1,i2=",2i4," f(solve)=",4f8.3)') i1,i2,f(0),f(1),f(2),f(3)
                 u1(i1-is1,i2-is2,i3,ex)=f(0)
                 u1(i1-is1,i2-is2,i3,ey)=f(1)
                 u2(j1-js1,j2-js2,j3,ex)=f(2)
                 u2(j1-js1,j2-js2,j3,ey)=f(3)
                if( debug>7 .and. twilightZone.eq.1 )then
                  ! check errors
                  call ogderiv(ep, 0,0,0,0, xy1(i1-is1,i2-is2,i3,0),xy1(i1-is1,i2-is2,i3,1),0.,t, ex, evv(0) )
                  call ogderiv(ep, 0,0,0,0, xy1(i1-is1,i2-is2,i3,0),xy1(i1-is1,i2-is2,i3,1),0.,t, ey, evv(1) )
                  call ogderiv(ep, 0,0,0,0, xy2(j1-js1,j2-js2,j3,0),xy2(j1-js1,j2-js2,j3,1),0.,t, ex, evv(2) )
                  call ogderiv(ep, 0,0,0,0, xy2(j1-js1,j2-js2,j3,0),xy2(j1-js1,j2-js2,j3,1),0.,t, ey, evv(3) )
                  maxErr=0.
                  do n=0,3
                    maxErr =max(maxErr,abs(evv(n)-f(n)))
                  end do
                  write(*,'("gdm22r: i1,i2,i3=",3i4," t=",e9.3)') i1,i2,i3,t 
                  write(*,'("     true= ",4(1pe10.2))') (evv(n),n=0,3)
                  write(*,'("      err= ",4(1pe8.1)," -> maxErr=",e8.1)') (abs(evv(n)-f(n)),n=0,3),maxErr
                end if
                 ! -------------------------------------------------------
                 ! No need to solve for Hz as it is just an ODE
                 ! -------------------------------------------------------
                  end if
                  j1=j1+1
                 end do
                 j2=j2+1
                end do
               ! stop 7777
             else if( useNonlinearModel.eq.0 )then
               ! dispersive case
             else
               ! nonlinear dispersive case
                  ! ****************************************************
                  ! ***********  2D, ORDER=2, RECTANGULAR **************
                  ! ****************************************************
                 if( t.le.3.*dt .and. debug.gt.0 )then
                   write(*,'("Interface>>>","22rectangle-nonlinear-MLA")')
                 end if
                 ! For rectangular, both sides must axis axis1==axis2: 
                 if( axis1.ne.axis2 )then
                   stop 8826
                 end if
                 ! 
                 ! Solve for the unknowns u1(-1),v1(-1),u2(-1),v2(-1)
                 !     
                 !       A [ U ] = A [ U(old) ] - [ f ]
                 !
                 !               [ u1(-1) ]
                 !       [ U ] = [ v1(-1) ]
                 !               [ u2(-1) ]
                 !               [ v2(-1) ]
                 !             
                 ! --- initialize some forcing functions ---
                 ! forcing functions for E and P
                 do n=0,nd-1
                   fev1(n)=0.
                   fev2(n)=0.
                   if (dispersionModel1 .ne. noDispersion) then
                     do jv=0,numberOfPolarizationVectors1-1
                       fpv1(n,jv)=0.
                     end do
                   endif
                   if (dispersionModel2 .ne. noDispersion) then
                     do jv=0,numberOfPolarizationVectors2-1
                       fpv2(n,jv)=0.
                     end do
                   endif
                 end do
                 ! forcing functions for N
                 if (nonlinearModel1 .ne. noNonlinearModel) then
                   do jv = 0,numberOfAtomicLevels1-1
                       fnv1(jv) = 0.
                       fntv1(jv) = 0.
                   enddo
                 endif
                 if (nonlinearModel2 .ne. noNonlinearModel) then
                   do jv = 0,numberOfAtomicLevels2-1
                       fnv2(jv) = 0.
                       fntv2(jv) = 0.
                   enddo
                 endif
                 ! print *, "-----------Now using MLA (RECTANGULAR)---------------"
                 ! ----------------- START LOOP OVER INTERFACE -------------------------
                  i3=n3a
                  j3=m3a
                  j2=m2a
                  do i2=n2a,n2b
                   j1=m1a
                   do i1=n1a,n1b
                    if( mask1(i1,i2,i3).gt.0 .and. mask2(j1,j2,j3).gt.0 )then
                   ! first evaluate the equations we want to solve with the wrong values at the ghost points:
                    ! NOTE: the jacobian derivatives can be computed once for all components
                      uu1=u1(i1,i2,i3,ex) ! in the rectangular case just eval the solution
                       u1x = (-u1(i1-1,i2,i3,ex)+u1(i1+1,i2,i3,ex))/(2.*dx1(0))
                       u1y = (-u1(i1,i2-1,i3,ex)+u1(i1,i2+1,i3,ex))/(2.*dx1(1))
                       u1xx = (u1(i1-1,i2,i3,ex)-2.*u1(i1,i2,i3,ex)+u1(i1+1,i2,i3,ex))/(dx1(0)**2)
                       u1yy = (u1(i1,i2-1,i3,ex)-2.*u1(i1,i2,i3,ex)+u1(i1,i2+1,i3,ex))/(dx1(1)**2)
                     u1Lap = u1xx+ u1yy
                      vv1=u1(i1,i2,i3,ey) ! in the rectangular case just eval the solution
                       v1x = (-u1(i1-1,i2,i3,ey)+u1(i1+1,i2,i3,ey))/(2.*dx1(0))
                       v1y = (-u1(i1,i2-1,i3,ey)+u1(i1,i2+1,i3,ey))/(2.*dx1(1))
                       v1xx = (u1(i1-1,i2,i3,ey)-2.*u1(i1,i2,i3,ey)+u1(i1+1,i2,i3,ey))/(dx1(0)**2)
                       v1yy = (u1(i1,i2-1,i3,ey)-2.*u1(i1,i2,i3,ey)+u1(i1,i2+1,i3,ey))/(dx1(1)**2)
                     v1Lap = v1xx+ v1yy
                    ! NOTE: the jacobian derivatives can be computed once for all components
                      uu2=u2(j1,j2,j3,ex) ! in the rectangular case just eval the solution
                       u2x = (-u2(j1-1,j2,j3,ex)+u2(j1+1,j2,j3,ex))/(2.*dx2(0))
                       u2y = (-u2(j1,j2-1,j3,ex)+u2(j1,j2+1,j3,ex))/(2.*dx2(1))
                       u2xx = (u2(j1-1,j2,j3,ex)-2.*u2(j1,j2,j3,ex)+u2(j1+1,j2,j3,ex))/(dx2(0)**2)
                       u2yy = (u2(j1,j2-1,j3,ex)-2.*u2(j1,j2,j3,ex)+u2(j1,j2+1,j3,ex))/(dx2(1)**2)
                     u2Lap = u2xx+ u2yy
                      vv2=u2(j1,j2,j3,ey) ! in the rectangular case just eval the solution
                       v2x = (-u2(j1-1,j2,j3,ey)+u2(j1+1,j2,j3,ey))/(2.*dx2(0))
                       v2y = (-u2(j1,j2-1,j3,ey)+u2(j1,j2+1,j3,ey))/(2.*dx2(1))
                       v2xx = (u2(j1-1,j2,j3,ey)-2.*u2(j1,j2,j3,ey)+u2(j1+1,j2,j3,ey))/(dx2(0)**2)
                       v2yy = (u2(j1,j2-1,j3,ey)-2.*u2(j1,j2,j3,ey)+u2(j1,j2+1,j3,ey))/(dx2(1)**2)
                     v2Lap = v2xx+ v2yy
                   ! Evaluate TZ forcing for dispersive equations in 2D 
                     if( twilightZone.eq.1 )then
                         if( dispersionModel1.ne.noDispersion .and. nonlinearModel1 .ne. noNonlinearModel) then
                           !-----------------------
                           ! dimension loops for E and P
                           !-----------------------
                           ! t is at new time now
                           nce = pxc+nd*numberOfPolarizationVectors1
                           do n=0,nd-1
                             fpSum1(n)=0.
                             pettSum1(n)=0.
                             call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, es(n)   ) 
                             call ogderiv(ep, 1,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, est(n)  )
                             call ogderiv(ep, 2,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, estt(n) )
                             call ogderiv(ep, 0,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esxx(n) )
                             call ogderiv(ep, 0,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esyy(n) )
                             do jv=0,numberOfPolarizationVectors1-1
                               ! The TZ component is offset by pxc
                               pc = pxc + jv*nd
                               call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pe(n)   )
                               call ogderiv(ep, 1,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pet(n)  )
                               call ogderiv(ep, 2,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pett(n) )
                               ! Normal TZ forcing for P_{n,jv} equation: 
                               ! fpv1(n,jv) = pett(n) + b1v1(jv)*pet(n) + b0v1(jv)*pe(n) - a0v(jv)*es(n) - a1v(jv)*est(n)
                               ! left hand side of gdm equations
                               fpv1(n,jv) = pett(n) + b1v1(jv)*pet(n) + b0v1(jv)*pe(n)
                               do na = 0,numberOfAtomicLevels1-1
                                 call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, nce+na, q0  )
                                 fpv1(n,jv) = fpv1(n,jv) - pnec1(jv,na)*q0*es(n) ! adding \Delta N*E
                                 ! print *, 'pnec1',pnec1(jv,na)
                               enddo
                               ! Keep sum: 
                               fpSum1(n)  = fpSum1(n)  + fpv1(n,jv)
                               pettSum1(n) = pettSum1(n) + pett(n) 
                             end do 
                             ! TZ forcing for E_{n} equation:
                             ! E_tt - c1^2 Delta E + alphaP1*Ptt  = 
                             fev1(n) = estt(n) - c1**2*( esxx(n) + esyy(n) ) + alphaP1*pettSum1(n)
                           end do
                           !--------------------------------
                           ! outside of dimension loop for N
                           !--------------------------------
                           do na=0,numberOfAtomicLevels1-1
                             ! na-th level
                             call ogderiv(ep, 1,0,0,0,xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, nce+na, q0t )
                             call ogderiv(ep, 2,0,0,0,xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, nce+na, q0tt)
                             ! initialize
                             fnv1(na)  = q0t ! forcing for \partial_tN_\ell = alpha_{\ell,k}N_k+\beta_{\ell,m}E\cdot\partial_tP_k
                             fntv1(na) = q0tt ! next derivative
                             ! relaxation (alpha_{\ell,m})
                             do jv=0,numberOfAtomicLevels1-1
                               call ogderiv(ep, 0,0,0,0,xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, nce+jv, q0  )
                               call ogderiv(ep, 1,0,0,0,xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, nce+jv, q0t )
                               fnv1(na)  = fnv1(na)  - prc1(na,jv)*q0
                               fntv1(na) = fntv1(na) - prc1(na,jv)*q0t
                             enddo
                             ! dot product (\beta_{\ell,k})
                             do n=0,nd-1 ! loop over dim
                               ! corresponding polarization vector
                               do jv=0,numberOfPolarizationVectors1-1  
                                 pc = pxc + jv*nd
                                 call ogderiv(ep, 1,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pet(n)  )
                                 call ogderiv(ep, 2,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pett(n) )
                                 fnv1(na)  = fnv1(na) - peptc1(na,jv)*es(n)*pet(n)
                                 fntv1(na) = fntv1(na) - peptc1(na,jv)*est(n)*pet(n) - peptc1(na,jv)*es(n)*pett(n)
                               enddo
                             enddo
                           enddo
                         end if
                         if( dispersionModel2.ne.noDispersion .and. nonlinearModel2 .ne. noNonlinearModel) then
                           !-----------------------
                           ! dimension loops for E and P
                           !-----------------------
                           ! t is at new time now
                           nce = pxc+nd*numberOfPolarizationVectors2
                           do n=0,nd-1
                             fpSum2(n)=0.
                             pettSum2(n)=0.
                             call ogderiv(ep, 0,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, es(n)   ) 
                             call ogderiv(ep, 1,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, est(n)  )
                             call ogderiv(ep, 2,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, estt(n) )
                             call ogderiv(ep, 0,2,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esxx(n) )
                             call ogderiv(ep, 0,0,2,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esyy(n) )
                             do jv=0,numberOfPolarizationVectors2-1
                               ! The TZ component is offset by pxc
                               pc = pxc + jv*nd
                               call ogderiv(ep, 0,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pe(n)   )
                               call ogderiv(ep, 1,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pet(n)  )
                               call ogderiv(ep, 2,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pett(n) )
                               ! Normal TZ forcing for P_{n,jv} equation: 
                               ! fpv2(n,jv) = pett(n) + b1v2(jv)*pet(n) + b0v2(jv)*pe(n) - a0v(jv)*es(n) - a1v(jv)*est(n)
                               ! left hand side of gdm equations
                               fpv2(n,jv) = pett(n) + b1v2(jv)*pet(n) + b0v2(jv)*pe(n)
                               do na = 0,numberOfAtomicLevels2-1
                                 call ogderiv(ep, 0,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t, nce+na, q0  )
                                 fpv2(n,jv) = fpv2(n,jv) - pnec2(jv,na)*q0*es(n) ! adding \Delta N*E
                                 ! print *, 'pnec2',pnec2(jv,na)
                               enddo
                               ! Keep sum: 
                               fpSum2(n)  = fpSum2(n)  + fpv2(n,jv)
                               pettSum2(n) = pettSum2(n) + pett(n) 
                             end do 
                             ! TZ forcing for E_{n} equation:
                             ! E_tt - c2^2 Delta E + alphaP2*Ptt  = 
                             fev2(n) = estt(n) - c2**2*( esxx(n) + esyy(n) ) + alphaP2*pettSum2(n)
                           end do
                           !--------------------------------
                           ! outside of dimension loop for N
                           !--------------------------------
                           do na=0,numberOfAtomicLevels2-1
                             ! na-th level
                             call ogderiv(ep, 1,0,0,0,xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t, nce+na, q0t )
                             call ogderiv(ep, 2,0,0,0,xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t, nce+na, q0tt)
                             ! initialize
                             fnv2(na)  = q0t ! forcing for \partial_tN_\ell = alpha_{\ell,k}N_k+\beta_{\ell,m}E\cdot\partial_tP_k
                             fntv2(na) = q0tt ! next derivative
                             ! relaxation (alpha_{\ell,m})
                             do jv=0,numberOfAtomicLevels2-1
                               call ogderiv(ep, 0,0,0,0,xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t, nce+jv, q0  )
                               call ogderiv(ep, 1,0,0,0,xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t, nce+jv, q0t )
                               fnv2(na)  = fnv2(na)  - prc2(na,jv)*q0
                               fntv2(na) = fntv2(na) - prc2(na,jv)*q0t
                             enddo
                             ! dot product (\beta_{\ell,k})
                             do n=0,nd-1 ! loop over dim
                               ! corresponding polarization vector
                               do jv=0,numberOfPolarizationVectors2-1  
                                 pc = pxc + jv*nd
                                 call ogderiv(ep, 1,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pet(n)  )
                                 call ogderiv(ep, 2,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pett(n) )
                                 fnv2(na)  = fnv2(na) - peptc2(na,jv)*es(n)*pet(n)
                                 fntv2(na) = fntv2(na) - peptc2(na,jv)*est(n)*pet(n) - peptc2(na,jv)*es(n)*pett(n)
                               enddo
                             enddo
                           enddo
                         end if
                     end if
                   ! eval dispersive forcings for domain 1
                     ! no dispersion
                     do n=0,nd-1
                       fp1(n)=0.
                     end do
                     ! nonlinear MLA
                     if( dispersionModel1.ne.noDispersion .and. nonlinearModel1.ne.noNonlinearModel) then
                       nce = pxc+nd*numberOfPolarizationVectors1
                       ! print *, '----------------',pxc,n
                       do n=0,nd-1
                         ! b0=b0v1(jv)
                         ! b1=b1v1(jv)
                         ! alpha=alphaP1
                         do jv=0,numberOfPolarizationVectors1-1
                           pc = n + jv*nd 
                           ec = ex +n
                           ! in time order: p1m,p1n,p1
                           ! pvm   =  p1(i1,i2,i3,pc) ! this one needs to be replaced due to rank difference
                           ! pv  =  p1n(i1,i2,i3,pc)
                           ! evm    =  u1(i1,i2,i3,ec)
                           ! ev   =  u1n(i1,i2,i3,ec)
                           pvn = 2.*p1(i1,i2,i3,pc)-p1n(i1,i2,i3,pc) + 0.5*dt*b1v1(jv)*p1n(i1,i2,i3,pc) - dtsq*b0v1(jv)*p1(i1,i2,i3,pc) + dtsq*fpv1(n,jv)
                           do na = 0,numberOfAtomicLevels1-1 ! \Delta N^n*E^n
                             pvn = pvn + dtsq*pnec1(jv,na)*q1(i1,i2,i3,na)*u1(i1,i2,i3,ec)
                           enddo
                           pvn= pvn/( 1.+.5*dt*b1v1(jv) )
                           ! #If "p1" eq "p1"
                           ! call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t+dt,pxc+jv*nd+n, pe(n)   )
                           ! #Else
                           ! call ogderiv(ep, 0,0,0,0, xy2(i1,i2,i3,0),xy2(i1,i2,i3,1),0.,t+dt,pxc+jv*nd+n, pe(n)   )
                           ! #End
                           ! print *, '+++++++pvn',pvn,'exact',pe(n),'diff',pvn-pe(n)
                           ! marching from t^{n+1} to t^{n+2}
                           ! #If "p1" eq "p1"
                           ! call ogderiv(ep, 2,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pett(n)   )
                           ! #Else
                           ! call ogderiv(ep, 2,0,0,0, xy2(i1,i2,i3,0),xy2(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pett(n)   )
                           ! #End
                           fp1(n) = fp1(n) + (pvn-2.*p1(i1,i2,i3,pc)+p1n(i1,i2,i3,pc))/dtsq
                           ! fp1(n) = fp1(n) + pett(n)
                           ! print *, '------pvtt',(pvn-2.*p1(i1,i2,i3,pc)+p1n(i1,i2,i3,pc))/dtsq,'exact',pett(n),'diff',(pvn-2.*p1(i1,i2,i3,pc)+p1n(i1,i2,i3,pc))/dtsq-pett(n)
                           ! print*,'dtsq',dtsq,'diff1',pvn-2.*p1(i1,i2,i3,pc)+p1n(i1,i2,i3,pc),'diff2',p1(i1,i2,i3,pc)-2.*p1n(i1,i2,i3,pc)+p1m(i1,i2,i3,pc)
                           ! write(*,'(" numberOfAtomicLevels1=",i3)') numberOfAtomicLevels1
                           ! write(*,'(" i1,i2,i3=",3i3)') i1,i2,i3
                           ! na=numberOfAtomicLevels1-1
                           ! write(*,'(" q1,q1n,q1m=",3e12.2)') q1(i1,i2,i3,na),q1n(i1,i2,i3,na),q1m(i1,i2,i3,na)
                           ! write(*,'(" pc=",i3," pvn,p1,p1n,p1m=",4e12.2)') pc, pvn,p1(i1,i2,i3,pc),p1n(i1,i2,i3,pc),p1m(i1,i2,i3,pc)
                           ! ! write(*,'(" dt=",e12.2," pv,pvt,pvtt, ev,evt,evtt=",6e12.2)') dt,pv,pvt,pvtt, ev,evt,evtt
                           ! ! write(*,'(" jv=",i2," a0,a1,b0,b1=",4e12.2," Bk,Ck=",2e12.2)') jv,a0v(jv),a1v(jv),b0v1(jv),b1v1(jv),Bk,Ck
                           ! write(*,'(" n=",i2," fev1(n)=",e12.2," fp1(n)=",e12.2," fpv1(n,jv)=",e12.2)') n,fev1(n),fp1(n),fpv1(n,jv)
                         end do
                       end do
                       ! we could precompute D
                       beta1 = 1.
                     else
                       beta1 = 1.
                     end if
                   ! eval dispersive forcings for domain 2
                     ! no dispersion
                     do n=0,nd-1
                       fp2(n)=0.
                     end do
                     ! nonlinear MLA
                     if( dispersionModel2.ne.noDispersion .and. nonlinearModel2.ne.noNonlinearModel) then
                       nce = pxc+nd*numberOfPolarizationVectors2
                       ! print *, '----------------',pxc,n
                       do n=0,nd-1
                         ! b0=b0v2(jv)
                         ! b1=b1v2(jv)
                         ! alpha=alphaP2
                         do jv=0,numberOfPolarizationVectors2-1
                           pc = n + jv*nd 
                           ec = ex +n
                           ! in time order: p2m,p2n,p2
                           ! pvm   =  p2(j1,j2,j3,pc) ! this one needs to be replaced due to rank difference
                           ! pv  =  p2n(j1,j2,j3,pc)
                           ! evm    =  u2(j1,j2,j3,ec)
                           ! ev   =  u2n(j1,j2,j3,ec)
                           pvn = 2.*p2(j1,j2,j3,pc)-p2n(j1,j2,j3,pc) + 0.5*dt*b1v2(jv)*p2n(j1,j2,j3,pc) - dtsq*b0v2(jv)*p2(j1,j2,j3,pc) + dtsq*fpv2(n,jv)
                           do na = 0,numberOfAtomicLevels2-1 ! \Delta N^n*E^n
                             pvn = pvn + dtsq*pnec2(jv,na)*q2(j1,j2,j3,na)*u2(j1,j2,j3,ec)
                           enddo
                           pvn= pvn/( 1.+.5*dt*b1v2(jv) )
                           ! #If "p2" eq "p1"
                           ! call ogderiv(ep, 0,0,0,0, xy1(j1,j2,j3,0),xy1(j1,j2,j3,1),0.,t+dt,pxc+jv*nd+n, pe(n)   )
                           ! #Else
                           ! call ogderiv(ep, 0,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t+dt,pxc+jv*nd+n, pe(n)   )
                           ! #End
                           ! print *, '+++++++pvn',pvn,'exact',pe(n),'diff',pvn-pe(n)
                           ! marching from t^{n+1} to t^{n+2}
                           ! #If "p2" eq "p1"
                           ! call ogderiv(ep, 2,0,0,0, xy1(j1,j2,j3,0),xy1(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pett(n)   )
                           ! #Else
                           ! call ogderiv(ep, 2,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pett(n)   )
                           ! #End
                           fp2(n) = fp2(n) + (pvn-2.*p2(j1,j2,j3,pc)+p2n(j1,j2,j3,pc))/dtsq
                           ! fp2(n) = fp2(n) + pett(n)
                           ! print *, '------pvtt',(pvn-2.*p2(j1,j2,j3,pc)+p2n(j1,j2,j3,pc))/dtsq,'exact',pett(n),'diff',(pvn-2.*p2(j1,j2,j3,pc)+p2n(j1,j2,j3,pc))/dtsq-pett(n)
                           ! print*,'dtsq',dtsq,'diff1',pvn-2.*p2(j1,j2,j3,pc)+p2n(j1,j2,j3,pc),'diff2',p2(j1,j2,j3,pc)-2.*p2n(j1,j2,j3,pc)+p2m(j1,j2,j3,pc)
                           ! write(*,'(" numberOfAtomicLevels2=",i3)') numberOfAtomicLevels2
                           ! write(*,'(" j1,j2,j3=",3i3)') j1,j2,j3
                           ! na=numberOfAtomicLevels2-1
                           ! write(*,'(" q2,q2n,q2m=",3e12.2)') q2(j1,j2,j3,na),q2n(j1,j2,j3,na),q2m(j1,j2,j3,na)
                           ! write(*,'(" pc=",i3," pvn,p2,p2n,p2m=",4e12.2)') pc, pvn,p2(j1,j2,j3,pc),p2n(j1,j2,j3,pc),p2m(j1,j2,j3,pc)
                           ! ! write(*,'(" dt=",e12.2," pv,pvt,pvtt, ev,evt,evtt=",6e12.2)') dt,pv,pvt,pvtt, ev,evt,evtt
                           ! ! write(*,'(" jv=",i2," a0,a1,b0,b1=",4e12.2," Bk,Ck=",2e12.2)') jv,a0v(jv),a1v(jv),b0v2(jv),b1v2(jv),Bk,Ck
                           ! write(*,'(" n=",i2," fev2(n)=",e12.2," fp2(n)=",e12.2," fpv2(n,jv)=",e12.2)') n,fev2(n),fp2(n),fpv2(n,jv)
                         end do
                       end do
                       ! we could precompute D
                       beta2 = 1.
                     else
                       beta2 = 1.
                     end if
                   if( axis1.eq.0 )then ! vertical interfaces
                     ! Interface equations for a boundary at x = 0 or x=1
                     ! ---- EQUATION 0 -----
                     ! [ u.x + v.y ] = 0
                     ! NOTE: if mu==mu2 then we do not need TZ forcing for this eqn:
                     f(0)=(u1x+v1y) - (u2x+v2y)
                     a4(0,0) = -is1/(2.*dx1(axis1))    ! coeff of u1(-1) from [u.x+v.y] 
                     a4(0,1) = 0.                      ! coeff of v1(-1) from [u.x+v.y] 
                     a4(0,2) =  js1/(2.*dx2(axis2))    ! coeff of u2(-1) from [u.x+v.y] 
                     a4(0,3) = 0.                      ! coeff of v2(-1) from [u.x+v.y]
                     ! ---- EQUATION 1 -----
                     ! [ (1/mu)* tv,.( curl(E) ) ] = 0  
                     ! NOTE: if mu==mu2 then we do not need TZ forcing for this eqn:
                     f(1)=(v1x-u1y)/mu1 - (v2x-u2y)/mu2
                     a4(1,0) = 0.
                     a4(1,1) = -is1/(2.*dx1(axis1))    ! coeff of v1(-1) from [v.x - u.y] 
                     a4(1,2) = 0.
                     a4(1,3) =  js1/(2.*dx2(axis2))    ! coeff of v2(-1) from [v.x - u.y]
                     ! ---- EQUATION 2 -----
                     ! [ (1/mu)* nv.( Delta(E) ) ]=0 (normal component)
                     ! NOTE: if mu1==mu2 then we do not need TZ forcing for this eqn (TZ forcing canceled due to nonzero jump conditions)
                     f(2)=( (u1xx+u1yy)/mu1  ) - ( (u2xx+u2yy)/mu2  )
                     a4(2,0) = 1./(dx1(axis1)**2)/mu1   ! coeff of u1(-1) from [(u.xx + u.yy)/mu]
                     a4(2,1) = 0. 
                     a4(2,2) =-1./(dx2(axis2)**2)/mu2   ! coeff of u2(-1) from [(u.xx + u.yy)/mu]
                     a4(2,3) = 0. 
                     ! ---- EQUATION 3 -----    
                     ! [ tv.( c^2*Delta(E) -alphaP*P_tt) ] = 0 (tangential component)
                     ! The coefficient of Delta(E) in this equation is altered due to Ptt term (not true for MLA)
                     f(3)=( (v1xx+v1yy)*beta1/epsmu1 -alphaP1*fp1(1) + fev1(1)) - ( (v2xx+v2yy)*beta2/epsmu2 -alphaP2*fp2(1) + fev2(1))
                     ! TEST 
                     if( .false. )then
                       f(3)=( (v1xx+v1yy)*beta1/epsmu1 ) - ( (v2xx+v2yy)*beta2/epsmu2 )
                     end if
                     a4(3,0) = 0.                      
                     a4(3,1) = (beta1/epsmu1)/(dx1(axis1)**2) ! coeff of v1(-1) from [beta*c^2*(v.xx+v.yy)]
                     a4(3,2) = 0. 
                     a4(3,3) =-(beta2/epsmu2)/(dx2(axis2)**2) ! coeff of v2(-1) from [beta*c^2*(v.xx+v.yy)]
                     ! print *, 'E TZ forcing (x)',fev1(0),fev2(0),'E TZ forcing (y)',fev1(1),fev2(1)
                     ! print *, '============eps:', eps1,eps2, 'mu',mu1,mu2, 'epsmu',epsmu1,epsmu2,'beta',beta1,beta2,'alphaP',alphaP1,alphaP2
                   else ! ---------- horizontal interfaces ---------------
                     ! Interface equations for a boundary at y = 0 or y=1
                     ! Switch u <-> v,  x<-> y in above equations 
                     ! ---- EQUATION 0 -----
                     f(0)=(v1y+u1x) - (v2y+u2x)
                     a4(0,0) = 0.                      ! coeff of u1(-1) from [u.x+v.y] 
                     a4(0,1) = -is1/(2.*dx1(axis1))    ! coeff of v1(-1) from [u.x+v.y] 
                     a4(0,2) = 0.                      ! coeff of u2(-1) from [u.x+v.y] 
                     a4(0,3) = js1/(2.*dx2(axis2))     ! coeff of v2(-1) from [u.x+v.y]
                     ! ---- EQUATION 1 -----
                     f(1)=(u1y-v1x)/mu1 - (u2y-v2x)/mu2
                     a4(1,0) = -is1/(2.*dx1(axis1))
                     a4(1,1) = 0.
                     a4(1,2) =  js1/(2.*dx2(axis2))  
                     a4(1,3) = 0.
                     ! ---- EQUATION 2 -----    
                     f(2)=( (v1xx+v1yy)/mu1 ) - ( (v2xx+v2yy)/mu2 )
                     a4(2,0) = 0.
                     a4(2,1) = 1./(dx1(axis1)**2)/mu1  
                     a4(2,2) = 0.
                     a4(2,3) =-1./(dx2(axis2)**2)/mu2 
                     ! ---- EQUATION 3 -----    
                     ! The coefficient of Delta(E) in this equation is altered due to Ptt term 
                     f(3)=( (u1xx+u1yy)*beta1/epsmu1 -alphaP1*fp1(0) +fev1(0) ) - ( (u2xx+u2yy)*beta2/epsmu2 -alphaP2*fp2(0) +fev2(0) )
                     a4(3,0) = (beta1/epsmu1)/(dx1(axis1)**2)
                     a4(3,1) = 0.
                     a4(3,2) =-(beta2/epsmu2)/(dx2(axis2)**2) 
                     a4(3,3) = 0.
                   end if
                    q(0) = u1(i1-is1,i2-is2,i3,ex)
                    q(1) = u1(i1-is1,i2-is2,i3,ey)
                    q(2) = u2(j1-js1,j2-js2,j3,ex)
                    q(3) = u2(j1-js1,j2-js2,j3,ey)
                    if( .false. .or. debug.gt.4 )then 
                      write(*,'("BEFORE: --> i1,i2=",2i4," j1,j2=",2i4," f()=",4e10.2)') i1,i2,j1,j2,f(0),f(1),f(2),f(3)
                      write(*,'("     beta1,beta2=",2e10.2," fp1=",2e10.2," fp2=",2e10.2)') beta1,beta2,fp1(0),fp1(1),fp2(0),fp2(1)
                      write(*,'("     mu1,mu2=",2e10.2," v1y,u1x,v2y,u2x=",4e10.2)') mu1,mu2,v1y,u1x,v2y,u2x
                    end if
                    ! subtract off the contributions from the wrong values at the ghost points:
                    do n=0,3
                      f(n) = (a4(n,0)*q(0)+a4(n,1)*q(1)+a4(n,2)*q(2)+a4(n,3)*q(3)) - f(n)
                    end do
                    ! write(debugFile,'(" --> i1,i2=",2i4," f(subtract)=",4f8.3)') i1,i2,f(0),f(1),f(2),f(3)
                    ! solve A Q = F
                    ! factor the matrix
                    numberOfEquations=4
                    call dgeco( a4(0,0), numberOfEquations, numberOfEquations, ipvt(0),rcond,work(0))
                    ! solve
                    ! write(debugFile,'(" --> i1,i2=",2i4," rcond=",e10.2)') i1,i2,rcond
                    job=0
                    call dgesl( a4(0,0), numberOfEquations, numberOfEquations, ipvt(0), f(0), job)
                    ! write(debugFile,'(" --> i1,i2=",2i4," f(solve)=",4f8.3)') i1,i2,f(0),f(1),f(2),f(3)
                    u1(i1-is1,i2-is2,i3,ex)=f(0)
                    u1(i1-is1,i2-is2,i3,ey)=f(1)
                    u2(j1-js1,j2-js2,j3,ex)=f(2)
                    u2(j1-js1,j2-js2,j3,ey)=f(3)
                    if( .false. .or. debug.gt.4 )then 
                      ! CHECK: re-evaluate the jump conditions
                       ! NOTE: the jacobian derivatives can be computed once for all components
                         uu1=u1(i1,i2,i3,ex) ! in the rectangular case just eval the solution
                          u1x = (-u1(i1-1,i2,i3,ex)+u1(i1+1,i2,i3,ex))/(2.*dx1(0))
                          u1y = (-u1(i1,i2-1,i3,ex)+u1(i1,i2+1,i3,ex))/(2.*dx1(1))
                          u1xx = (u1(i1-1,i2,i3,ex)-2.*u1(i1,i2,i3,ex)+u1(i1+1,i2,i3,ex))/(dx1(0)**2)
                          u1yy = (u1(i1,i2-1,i3,ex)-2.*u1(i1,i2,i3,ex)+u1(i1,i2+1,i3,ex))/(dx1(1)**2)
                        u1Lap = u1xx+ u1yy
                         vv1=u1(i1,i2,i3,ey) ! in the rectangular case just eval the solution
                          v1x = (-u1(i1-1,i2,i3,ey)+u1(i1+1,i2,i3,ey))/(2.*dx1(0))
                          v1y = (-u1(i1,i2-1,i3,ey)+u1(i1,i2+1,i3,ey))/(2.*dx1(1))
                          v1xx = (u1(i1-1,i2,i3,ey)-2.*u1(i1,i2,i3,ey)+u1(i1+1,i2,i3,ey))/(dx1(0)**2)
                          v1yy = (u1(i1,i2-1,i3,ey)-2.*u1(i1,i2,i3,ey)+u1(i1,i2+1,i3,ey))/(dx1(1)**2)
                        v1Lap = v1xx+ v1yy
                       ! NOTE: the jacobian derivatives can be computed once for all components
                         uu2=u2(j1,j2,j3,ex) ! in the rectangular case just eval the solution
                          u2x = (-u2(j1-1,j2,j3,ex)+u2(j1+1,j2,j3,ex))/(2.*dx2(0))
                          u2y = (-u2(j1,j2-1,j3,ex)+u2(j1,j2+1,j3,ex))/(2.*dx2(1))
                          u2xx = (u2(j1-1,j2,j3,ex)-2.*u2(j1,j2,j3,ex)+u2(j1+1,j2,j3,ex))/(dx2(0)**2)
                          u2yy = (u2(j1,j2-1,j3,ex)-2.*u2(j1,j2,j3,ex)+u2(j1,j2+1,j3,ex))/(dx2(1)**2)
                        u2Lap = u2xx+ u2yy
                         vv2=u2(j1,j2,j3,ey) ! in the rectangular case just eval the solution
                          v2x = (-u2(j1-1,j2,j3,ey)+u2(j1+1,j2,j3,ey))/(2.*dx2(0))
                          v2y = (-u2(j1,j2-1,j3,ey)+u2(j1,j2+1,j3,ey))/(2.*dx2(1))
                          v2xx = (u2(j1-1,j2,j3,ey)-2.*u2(j1,j2,j3,ey)+u2(j1+1,j2,j3,ey))/(dx2(0)**2)
                          v2yy = (u2(j1,j2-1,j3,ey)-2.*u2(j1,j2,j3,ey)+u2(j1,j2+1,j3,ey))/(dx2(1)**2)
                        v2Lap = v2xx+ v2yy
                      if( axis1.eq.0 )then
                         f(0)=(u1x+v1y) - (u2x+v2y)
                         f(1)=(v1x-u1y)/mu1 - (v2x-u2y)/mu2
                         f(2)=( (u1xx+u1yy)/mu1 ) - ( (u2xx+u2yy)/mu2 )
                         f(3)=( (v1xx+v1yy)*beta1/epsmu1 -alphaP1*fp1(1) ) - ( (v2xx+v2yy)*beta2/epsmu2 -alphaP2*fp2(1))
                        ! TEST 
                         if( .false. )then
                           f(3)=( (v1xx+v1yy)*beta1/epsmu1 ) - ( (v2xx+v2yy)*beta2/epsmu2 )
                         end if
                       else
                         f(0)=(v1y+u1x) - (v2y+u2x)
                         f(1)=(u1y-v1x)/mu1 - (u2y-v2x)/mu2
                         f(2)=( (v1xx+v1yy)/mu1 ) - ( (v2xx+v2yy)/mu2 )    
                         f(3)=( (u1xx+u1yy)*beta1/epsmu1 -alphaP1*fp1(0) ) - ( (u2xx+u2yy)*beta2/epsmu2 -alphaP2*fp2(0))
                       end if 
                       write(*,'("AFTER: --> i1,i2=",2i4," j1,j2=",2i4," f(re-eval)=",4e10.2)') i1,i2,j1,j2,f(0),f(1),f(2),f(3)
                       if( twilightZone.eq.1 )then
                         ! check errors in the ghost 
                           k1=i1-is1
                           k2=i2-is2
                           k3=i3
                           do n=0,nd-1
                             call ogderiv(ep, 0,0,0,0, xy1(k1,k2,k3,0),xy1(k1,k2,k3,1),0.,t,ex+n, es(n)   ) 
                           end do
                           f(0) =  u1(i1-is1,i2-is2,i3,ex) -es(0)
                           f(1) =  u1(i1-is1,i2-is2,i3,ey) -es(1)
                           k1=j1-js1
                           k2=j2-js2
                           k3=j3
                           do n=0,nd-1
                             call ogderiv(ep, 0,0,0,0, xy2(k1,k2,k3,0),xy2(k1,k2,k3,1),0.,t,ex+n, est(n)   ) 
                           end do
                           f(2) =  u2(j1-js1,j2-js2,j3,ex) -est(0)
                           f(3) =  u2(j1-js1,j2-js2,j3,ey) -est(1)
                           write(*,'(" ghost err =",4e10.2)') f(0),f(1),f(2),f(3) 
                       end if
                    end if
                    ! -------------------------------------------------------
                    ! No need to solve for Hz as it is just an ODE
                    ! -------------------------------------------------------
                     end if
                     j1=j1+1
                    end do
                    j2=j2+1
                   end do
                  ! stop 7777
             end if    
               ! grid1: ns1a,ns1b, ...  are saved values of n1a,n1b,...
               n1a=ns1a
               n1b=ns1b
               n2a=ns2a
               n2b=ns2b
               n3a=ns3a
               n3b=ns3b
               ! grid2
               m1a=ms1a
               m1b=ms1b
               m2a=ms2a
               m2b=ms2b
               m3a=ms3a
               m3b=ms3b
             orderOfAccuracy=4 ! reset 
              ! This next perl command will reset the macro derivatives to be computed to order=4
             ! Macro to assign ghost values:
             if( dispersive.eq.0 )then
             else if( useUnifiedInterfaceMacros.eq.1 ) then
              ! Use Qing's unified interface macros
               ! ****************************************************************
               ! ***********  Unified, 2D, ORDER=4, RECTANGULAR **************
               ! ****************************************************************
               if( t.le.3.*dt .and. debug.gt.0 )then
                 write(*,'("Interface>>>","24r-unified-dispersive")')
               end if
               if( t.le.5*dt .or. debug.gt.3 )then
                 if( it.le.2 )then
                   write(*,'("macro: assignUnifiedInterfaceGhost24r : it=",i6," t,dt=",2e10.2)') it,t,dt
                 end if
               end if
               ! normal and tangent (for TZ forcing)
               an1=an1Cartesian
               an2=an2Cartesian
               tau1=-an2
               tau2= an1
              ! make sure the normal and tangent are set
              if( abs( an1**2 + an2**2 -1. )>1.e-10 .or. abs( tau1**2 + tau2**2 -1. )>1.e-10 )then
                write(*,'("unified-24r - ERROR: incorrect an1,an2, tau1,tau2=",4(1pe9.2))') an1,an2,tau1,tau2
                stop 6666
              end if
              ! --- initialize some forcing functions ---
              do n=0,nd-1
                fev1(n)=0.
                LfE1(n)=0.
                fEt1(n)=0.
                fEtt1(n)=0.
                fev2(n)=0.
                LfE2(n)=0.
                fEt2(n)=0.
                fEtt2(n)=0.
                fevx1(n)=0.
                fevy1(n)=0.
                fevx2(n)=0.
                fevy2(n)=0.
                if (dispersionModel1 .ne. 0) then
                  do jv=0,numberOfPolarizationVectors1-1
                    fpv1(n,jv)=0.
                    LfP1(n,jv)=0.
                    fPt1(n,jv)=0.
                    fPtt1(n,jv)=0.
                    fpvx1(n,jv)=0.
                    fpvy1(n,jv)=0.
                  end do
                endif
                if (dispersionModel2 .ne. 0) then
                  do jv=0,numberOfPolarizationVectors2-1
                    fpv2(n,jv)=0.
                    LfP2(n,jv)=0.
                    fPt2(n,jv)=0.
                    fPtt2(n,jv)=0.
                    fpvx2(n,jv)=0.
                    fpvy2(n,jv)=0.
                  end do
                endif
              end do
              ! forcing functions for N
              if (nonlinearModel1 .ne. 0) then
                do jv = 0,numberOfAtomicLevels1-1
                    fnv1(jv) = 0.
                    fntv1(jv) = 0.
                enddo
              endif
              if (nonlinearModel2 .ne. 0) then
                do jv = 0,numberOfAtomicLevels2-1
                    fnv2(jv) = 0.
                    fntv2(jv) = 0.
                enddo
              endif
               ! write(*,'("p1=",(15(e10.2,1x)))') (((p1(i1,i2,i3,0),i1=nd1a,nd1b),i2=nd2a,nd2b),i3=nd3a,nd3b)
               ! =============== start loops ======================
                i3=n3a
                j3=m3a
                j2=m2a
                do i2=n2a,n2b
                 j1=m1a
                 do i1=n1a,n1b
                  if( mask1(i1,i2,i3).gt.0 .and. mask2(j1,j2,j3).gt.0 )then
                 ! Evaluate the jump conditions using the wrong values at the ghost points 
                  ! Evaluate TZ forcing for dispersive equations in 2D 
                    if( twilightZone.eq.1 )then
                      !-----------------------
                      ! MLA
                      !-----------------------
                      if( dispersionModel1.ne.noDispersion .and. nonlinearModel1.ne.noNonlinearModel) then
                        nce = pxc+nd*numberOfPolarizationVectors1
                        !-----------------------------
                        ! dimension loop for P and E
                        !-----------------------------
                        nce = pxc+nd*numberOfPolarizationVectors1
                        do n=0,nd-1
                          fpSum1(n)=0.
                          pevttSum1(n)=0.
                          pevttxSum1(n)=0.
                          pevttySum1(n)=0.
                          pevttLSum1(n)=0.
                          pevttttSum1(n)=0.
                          petttSum=0.
                          call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, es(n)   ) 
                          call ogderiv(ep, 1,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, est(n)  )
                          call ogderiv(ep, 2,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, estt(n) )
                          call ogderiv(ep, 0,1,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esx(n) )
                          call ogderiv(ep, 0,0,1,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esy(n) )
                          call ogderiv(ep, 0,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esxx(n) )
                          call ogderiv(ep, 0,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esyy(n) )
                          call ogderiv(ep, 1,1,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, estx(n) )
                          call ogderiv(ep, 1,0,1,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esty(n) )
                          call ogderiv(ep, 2,1,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esttx(n) )
                          call ogderiv(ep, 2,0,1,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, estty(n) )
                          call ogderiv(ep, 0,3,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esxxx(n) )
                          call ogderiv(ep, 0,2,1,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esxxy(n) )
                          call ogderiv(ep, 0,1,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esxyy(n) )
                          call ogderiv(ep, 0,0,3,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esyyy(n) )
                          call ogderiv(ep, 3,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esttt(n) )
                          call ogderiv(ep, 4,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, estttt(n) )
                          call ogderiv(ep, 1,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, estxx(n) )
                          call ogderiv(ep, 1,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, estyy(n) )
                          call ogderiv(ep, 2,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esttxx(n) )
                          call ogderiv(ep, 2,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esttyy(n) )
                          call ogderiv(ep, 0,4,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esxxxx(n) )
                          call ogderiv(ep, 0,2,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esxxyy(n) )
                          call ogderiv(ep, 0,0,4,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esyyyy(n) )
                          ! L = c1^2*Delta
                          esL  = ( esxx(n)   + esyy(n) ) ! deleted c1^2
                          estL = ( estxx(n)  + estyy(n) )
                          esttL= ( esttxx(n) + esttyy(n) )
                          esLx  = ( esxxx(n)   + esxyy(n) )
                          esLy  = ( esxxy(n)   + esyyy(n) )
                          ! L^2 : 
                          esLL = ( esxxxx(n)  + 2.*esxxyy(n) + esyyyy(n) ) ! deleted c1^4
                          do jv=0,numberOfPolarizationVectors1-1
                            ! The TZ component is offset by pxc
                            pc = pxc + jv*nd
                            call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pe(n)   )
                            call ogderiv(ep, 1,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pet(n)  )
                            call ogderiv(ep, 2,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pett(n) )
                            call ogderiv(ep, 3,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pettt(n) )
                            call ogderiv(ep, 4,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, petttt(n) )
                            call ogderiv(ep, 0,1,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pex(n) )
                            call ogderiv(ep, 0,0,1,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pey(n) )
                            call ogderiv(ep, 0,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pexx(n) )
                            call ogderiv(ep, 0,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, peyy(n) )
                            call ogderiv(ep, 1,1,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, petx(n) )
                            call ogderiv(ep, 1,0,1,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pety(n) )
                            call ogderiv(ep, 1,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, petxx(n) )
                            call ogderiv(ep, 1,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, petyy(n) )
                            call ogderiv(ep, 2,1,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pettx(n) )
                            call ogderiv(ep, 2,0,1,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, petty(n) )
                            call ogderiv(ep, 2,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pettxx(n) )
                            call ogderiv(ep, 2,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pettyy(n) )
                            peL  = ( pexx(n)   + peyy(n) ) ! deleted c1^2
                            petL = ( petxx(n)  + petyy(n) )
                            pettL= ( pettxx(n) + pettyy(n) )
                            fpv1(n,jv) = pett(n)   + b1v1(jv)*pet(n)   + b0v1(jv)*pe(n)
                            fPt1(n,jv) = pettt(n)  + b1v1(jv)*pett(n)  + b0v1(jv)*pet(n)
                            fPtt1(n,jv)= petttt(n) + b1v1(jv)*pettt(n) + b0v1(jv)*pett(n)
                            LfP1(n,jv) = pettL     + b1v1(jv)*petL     + b0v1(jv)*peL
                            fpvx1(n,jv)= pettx(n)  + b1v1(jv)*petx(n)  + b0v1(jv)*pex(n)
                            fpvy1(n,jv)= petty(n)  + b1v1(jv)*pety(n)  + b0v1(jv)*pey(n)
                            do na = 0,numberOfAtomicLevels1-1
                              call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, nce+na, q0  )
                              call ogderiv(ep, 1,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, nce+na, q0t  )
                              call ogderiv(ep, 2,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, nce+na, q0tt  )
                              call ogderiv(ep, 0,1,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, nce+na, q0x  )
                              call ogderiv(ep, 0,0,1,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, nce+na, q0y  )
                              call ogderiv(ep, 0,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, nce+na, q0xx  )
                              call ogderiv(ep, 0,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, nce+na, q0yy  )
                              fpv1(n,jv) = fpv1(n,jv) - pnec1(jv,na)*q0*es(n) ! adding \Delta N*E
                              fPt1(n,jv) = fPt1(n,jv) - pnec1(jv,na)*q0t*es(n) - pnec1(jv,na)*q0*est(n)
                              fPtt1(n,jv) = fPtt1(n,jv) - pnec1(jv,na)*q0tt*es(n)- 2.0*pnec1(jv,na)*q0t*est(n) - pnec1(jv,na)*q0*estt(n)
                              LfP1(n,jv) = LfP1(n,jv) - pnec1(jv,na)*(q0xx*es(n)+2.*q0x*esx(n)+q0*esxx(n) + q0yy*es(n)+2.*q0y*esy(n)+q0*esyy(n))
                              fpvx1(n,jv) = fpvx1(n,jv) - pnec1(jv,na)*q0x*es(n) - pnec1(jv,na)*q0*esx(n)
                              fpvy1(n,jv) = fpvy1(n,jv) - pnec1(jv,na)*q0y*es(n) - pnec1(jv,na)*q0*esy(n)
                            enddo
                            ! print *,'xxxxxxxxxxxxxxxxxxxxxxx'
                            ! print *, 'FOR P TZ FORCING'
                            ! print *, n,jv, fpv1(n,jv),fPt1(n,jv),fPtt1(n,jv),LfP1(n,jv),fpvx1(n,jv),fpvy1(n,jv)
                            ! print *,'xxxxxxxxxxxxxxxxxxxxxxx'
                            ! Normal TZ forcing for P_{n,jv} equation: 
                            ! fpv1(n,jv) = pett(n)   + b1v1(jv)*pet(n)   + b0v1(jv)*pe(n)   - a0v1(jv)*es(n)   - a1v1(jv)*est(n)
                            ! fPt1(n,jv) = pettt(n)  + b1v1(jv)*pett(n)  + b0v1(jv)*pet(n)  - a0v1(jv)*est(n)  - a1v1(jv)*estt(n)
                            ! fPtt1(n,jv)= petttt(n) + b1v1(jv)*pettt(n) + b0v1(jv)*pett(n) - a0v1(jv)*estt(n) - a1v1(jv)*esttt(n)
                            ! LfP1(n,jv) = pettL     + b1v1(jv)*petL     + b0v1(jv)*peL     - a0v1(jv)*esL     - a1v1(jv)*estL
                            ! fpvx1(n,jv)= pettx(n)  + b1v1(jv)*petx(n)  + b0v1(jv)*pex(n)  - a0v1(jv)*esx(n)  - a1v1(jv)*estx(n)
                            ! fpvy1(n,jv)= petty(n)  + b1v1(jv)*pety(n)  + b0v1(jv)*pey(n)  - a0v1(jv)*esy(n)  - a1v1(jv)*esty(n)
                            ! write(*,'(" n=",i4," LfP1=",e10.4," pettL,petL,peL,esL,estL=",5e12.4)') n,LfP1(n,jv),pettL,petL,peL,esL,estL
                            ! write(*,'(" pe,pet,pett,pettt,petttt=",5e12.4)') pe(n),pet(n),pett(n),pettt(n),petttt(n)
                            ! write(*,'("TZ: n,jv=",2i4," pex,pey,pexx,peyy=",4(1pe12.4))') n,jv,pex(n),pey(n),pexx(n),peyy(n)
                            ! Save ptt for checking later
                            pevtt1(n,jv)=pett(n)
                            pevttx1(n,jv)=pettx(n)
                            pevtty1(n,jv)=petty(n)
                            pevttt1(n,jv)=pettt(n)
                            pevtttt1(n,jv)=petttt(n)
                            pevttL1(n,jv) = pettL
                            pevttLSum1(n) = pevttLSum1(n)  + c1**2*pettL ! added c1**2 to be consistence with GDM for eval jumps
                            pevttttSum1(n)= pevttttSum1(n) + petttt(n) 
                            ! Keep some sums: 
                            fpSum1(n)   = fpSum1(n)  + fpv1(n,jv)
                            pevttSum1(n)  = pevttSum1(n)  + pett(n) 
                            pevttxSum1(n) = pevttxSum1(n) + pettx(n)
                            pevttySum1(n) = pevttySum1(n) + petty(n)
                            petttSum  = petttSum  + pettt(n) 
                          end do 
                          ! TZ forcing for E_{n} equation:
                          ! E_tt - c1^2 Delta E + alphaP1*Ptt  = 
                          fev1(n) = estt(n)   - c1**2*esL   + alphaP1*pevttSum1(n)
                          fEt1(n) = esttt(n)  - c1**2*estL  + alphaP1*petttSum
                          fEtt1(n)= estttt(n) - c1**2*esttL + alphaP1*pevttttSum1(n)
                          fevx1(n) = esttx(n) - c1**2*esLx   + alphaP1*pevttxSum1(n)
                          fevy1(n) = estty(n) - c1**2*esLy   + alphaP1*pevttySum1(n)
                          ! write(*,'("--> fEtt1=",e10.2," estttt,esttL,pettttSum=",3e10.2)')  fEtt1(n),estttt(n),esttL,pettttSum
                          LfE1(n) = esttL     - c1**2*esLL  + alphaP1*pevttLSum1(n)/c1**2 ! new: divided c1**2
                        end do
                        !--------------------------------
                        ! outside of dimension loop for N
                        !--------------------------------
                        do na=0,numberOfAtomicLevels1-1
                          ! na-th level
                          call ogderiv(ep, 1,0,0,0,xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, nce+na, q0t )
                          call ogderiv(ep, 2,0,0,0,xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, nce+na, q0tt)
                          call ogderiv(ep, 3,0,0,0,xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, nce+na, q0ttt)
                          call ogderiv(ep, 4,0,0,0,xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, nce+na, q0tttt)
                          ! initialize
                          fnv1(na)  = q0t ! forcing for \partial_tN_\ell = alpha_{\ell,k}N_k+\beta_{\ell,m}E\cdot\partial_tP_k
                          fntv1(na) = q0tt ! next derivative
                          fnttv1(na) = q0ttt
                          fntttv1(na) = q0tttt
                          ! relaxation (alpha_{\ell,m})
                          do jv=0,numberOfAtomicLevels1-1
                            call ogderiv(ep, 0,0,0,0,xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, nce+jv, q0 )
                            call ogderiv(ep, 1,0,0,0,xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, nce+jv, q0t)
                            call ogderiv(ep, 2,0,0,0,xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, nce+jv, q0tt)
                            call ogderiv(ep, 3,0,0,0,xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, nce+jv, q0ttt)
                            fnv1(na)  = fnv1(na)  - prc1(na,jv)*q0
                            fntv1(na) = fntv1(na) - prc1(na,jv)*q0t
                            fnttv1(na) = fnttv1(na) - prc1(na,jv)*q0tt
                            fntttv1(na) = fntttv1(na) - prc1(na,jv)*q0ttt
                          enddo
                          ! dot product (\beta_{\ell,k})
                          do n=0,nd-1 ! loop over dim
                            call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, es(n)   ) 
                            call ogderiv(ep, 1,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, est(n)  )
                            call ogderiv(ep, 2,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, estt(n) )
                            call ogderiv(ep, 3,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esttt(n) )
                            ! corresponding polarization vector
                            do jv=0,numberOfPolarizationVectors1-1 
                              pc = pxc + jv*nd
                              call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pe(n)   )
                              call ogderiv(ep, 1,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pet(n)  )
                              call ogderiv(ep, 2,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pett(n) )
                              call ogderiv(ep, 3,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pettt(n) )
                              call ogderiv(ep, 4,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, petttt(n) ) 
                              fnv1(na)  = fnv1(na) - peptc1(na,jv)*es(n)*pet(n)
                              fntv1(na) = fntv1(na) - peptc1(na,jv)*est(n)*pet(n) - peptc1(na,jv)*es(n)*pett(n)
                              fnttv1(na) = fnttv1(na) - peptc1(na,jv)*estt(n)*pet(n) - 2.d0*peptc1(na,jv)*est(n)*pett(n) - peptc1(na,jv)*es(n)*pettt(n)
                              fntttv1(na) = fntttv1(na) - peptc1(na,jv)*esttt(n)*pet(n) - 3.d0*peptc1(na,jv)*estt(n)*pett(n) - 3.d0*peptc1(na,jv)*est(n)*pettt(n) - peptc1(na,jv)*es(n)*petttt(n)
                            enddo
                          enddo
                        enddo
                      !-----------------------
                      ! GDM
                      !-----------------------
                      elseif( dispersionModel1.ne.noDispersion .and. nonlinearModel1.eq.noNonlinearModel) then
                        do n=0,nd-1
                          fpSum1(n)=0.
                          pevttSum1(n)=0.
                          pevttxSum1(n)=0.
                          pevttySum1(n)=0.
                          pevttLSum1(n)=0.
                          pevttttSum1(n)=0.
                          petttSum=0.
                          call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, es(n)   ) 
                          call ogderiv(ep, 1,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, est(n)  )
                          call ogderiv(ep, 2,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, estt(n) )
                          call ogderiv(ep, 0,1,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esx(n) )
                          call ogderiv(ep, 0,0,1,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esy(n) )
                          call ogderiv(ep, 0,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esxx(n) )
                          call ogderiv(ep, 0,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esyy(n) )
                          call ogderiv(ep, 1,1,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, estx(n) )
                          call ogderiv(ep, 1,0,1,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esty(n) )
                          call ogderiv(ep, 2,1,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esttx(n) )
                          call ogderiv(ep, 2,0,1,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, estty(n) )
                          call ogderiv(ep, 0,3,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esxxx(n) )
                          call ogderiv(ep, 0,2,1,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esxxy(n) )
                          call ogderiv(ep, 0,1,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esxyy(n) )
                          call ogderiv(ep, 0,0,3,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esyyy(n) )
                          call ogderiv(ep, 3,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esttt(n) )
                          call ogderiv(ep, 4,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, estttt(n) )
                          call ogderiv(ep, 1,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, estxx(n) )
                          call ogderiv(ep, 1,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, estyy(n) )
                          call ogderiv(ep, 2,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esttxx(n) )
                          call ogderiv(ep, 2,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esttyy(n) )
                          call ogderiv(ep, 0,4,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esxxxx(n) )
                          call ogderiv(ep, 0,2,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esxxyy(n) )
                          call ogderiv(ep, 0,0,4,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esyyyy(n) )
                          ! L = c1^2*Delta
                          esL  = c1**2*( esxx(n)   + esyy(n) )
                          estL = c1**2*( estxx(n)  + estyy(n) )
                          esttL= c1**2*( esttxx(n) + esttyy(n) )
                          esLx  = c1**2*( esxxx(n)   + esxyy(n) )
                          esLy  = c1**2*( esxxy(n)   + esyyy(n) )
                          ! L^2 : 
                          esLL = c1**4*( esxxxx(n)  + 2.*esxxyy(n) + esyyyy(n) )
                          do jv=0,numberOfPolarizationVectors1-1
                            ! The TZ component is offset by pxc
                            pc = pxc + jv*nd
                            call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pe(n)   )
                            call ogderiv(ep, 1,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pet(n)  )
                            call ogderiv(ep, 2,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pett(n) )
                            call ogderiv(ep, 3,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pettt(n) )
                            call ogderiv(ep, 4,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, petttt(n) )
                            call ogderiv(ep, 0,1,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pex(n) )
                            call ogderiv(ep, 0,0,1,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pey(n) )
                            call ogderiv(ep, 0,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pexx(n) )
                            call ogderiv(ep, 0,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, peyy(n) )
                            call ogderiv(ep, 1,1,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, petx(n) )
                            call ogderiv(ep, 1,0,1,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pety(n) )
                            call ogderiv(ep, 1,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, petxx(n) )
                            call ogderiv(ep, 1,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, petyy(n) )
                            call ogderiv(ep, 2,1,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pettx(n) )
                            call ogderiv(ep, 2,0,1,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, petty(n) )
                            call ogderiv(ep, 2,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pettxx(n) )
                            call ogderiv(ep, 2,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pettyy(n) )
                            peL  = c1**2*( pexx(n)   + peyy(n) )
                            petL = c1**2*( petxx(n)  + petyy(n) )
                            pettL= c1**2*( pettxx(n) + pettyy(n) )
                            ! Normal TZ forcing for P_{n,jv} equation: 
                            fpv1(n,jv) = pett(n)   + b1v1(jv)*pet(n)   + b0v1(jv)*pe(n)   - a0v1(jv)*es(n)   - a1v1(jv)*est(n)
                            fPt1(n,jv) = pettt(n)  + b1v1(jv)*pett(n)  + b0v1(jv)*pet(n)  - a0v1(jv)*est(n)  - a1v1(jv)*estt(n)
                            fPtt1(n,jv)= petttt(n) + b1v1(jv)*pettt(n) + b0v1(jv)*pett(n) - a0v1(jv)*estt(n) - a1v1(jv)*esttt(n)
                            LfP1(n,jv) = pettL     + b1v1(jv)*petL     + b0v1(jv)*peL     - a0v1(jv)*esL     - a1v1(jv)*estL
                            fpvx1(n,jv)= pettx(n)  + b1v1(jv)*petx(n)  + b0v1(jv)*pex(n)  - a0v1(jv)*esx(n)  - a1v1(jv)*estx(n)
                            fpvy1(n,jv)= petty(n)  + b1v1(jv)*pety(n)  + b0v1(jv)*pey(n)  - a0v1(jv)*esy(n)  - a1v1(jv)*esty(n)
                            ! write(*,'(" n=",i4," LfP1=",e10.4," pettL,petL,peL,esL,estL=",5e12.4)') n,LfP1(n,jv),pettL,petL,peL,esL,estL
                            ! write(*,'(" pe,pet,pett,pettt,petttt=",5e12.4)') pe(n),pet(n),pett(n),pettt(n),petttt(n)
                            ! write(*,'("TZ: n,jv=",2i4," pex,pey,pexx,peyy=",4(1pe12.4))') n,jv,pex(n),pey(n),pexx(n),peyy(n)
                            ! Save ptt for checking later
                            pevtt1(n,jv)=pett(n)
                            pevttx1(n,jv)=pettx(n)
                            pevtty1(n,jv)=petty(n)
                            pevtttt1(n,jv)=petttt(n)
                            pevttLSum1(n) = pevttLSum1(n)  + pettL
                            pevttttSum1(n)= pevttttSum1(n) + petttt(n) 
                            ! Keep some sums: 
                            fpSum1(n)   = fpSum1(n)  + fpv1(n,jv)
                            pevttSum1(n)  = pevttSum1(n)  + pett(n) 
                            pevttxSum1(n) = pevttxSum1(n) + pettx(n)
                            pevttySum1(n) = pevttySum1(n) + petty(n)
                            petttSum  = petttSum  + pettt(n) 
                          end do 
                          ! TZ forcing for E_{n} equation:
                          ! E_tt - c1^2 Delta E + alphaP1*Ptt  = 
                          fev1(n) = estt(n)   - esL   + alphaP1*pevttSum1(n)
                          fEt1(n) = esttt(n)  - estL  + alphaP1*petttSum
                          fEtt1(n)= estttt(n) - esttL + alphaP1*pevttttSum1(n)
                          fevx1(n) = esttx(n) - esLx   + alphaP1*pevttxSum1(n)
                          fevy1(n) = estty(n) - esLy   + alphaP1*pevttySum1(n)
                          ! write(*,'("--> fEtt1=",e10.2," estttt,esttL,pettttSum=",3e10.2)')  fEtt1(n),estttt(n),esttL,pettttSum
                          LfE1(n) = esttL     - esLL  + alphaP1*pevttLSum1(n)
                       end do
                      !-----------------------
                      ! no dispersion
                      !-----------------------
                      elseif( dispersionModel1.eq.noDispersion .and. nonlinearModel1.eq.noNonlinearModel) then
                        do n=0,nd-1
                          fpSum1(n)=0.
                          pevttSum1(n)=0.
                          pevttxSum1(n)=0.
                          pevttySum1(n)=0.
                          pevttLSum1(n)=0.
                          pevttttSum1(n)=0.
                        enddo
                      end if
                      !-----------------------
                      ! MLA
                      !-----------------------
                      if( dispersionModel2.ne.noDispersion .and. nonlinearModel2.ne.noNonlinearModel) then
                        nce = pxc+nd*numberOfPolarizationVectors2
                        !-----------------------------
                        ! dimension loop for P and E
                        !-----------------------------
                        nce = pxc+nd*numberOfPolarizationVectors2
                        do n=0,nd-1
                          fpSum2(n)=0.
                          pevttSum2(n)=0.
                          pevttxSum2(n)=0.
                          pevttySum2(n)=0.
                          pevttLSum2(n)=0.
                          pevttttSum2(n)=0.
                          petttSum=0.
                          call ogderiv(ep, 0,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, es(n)   ) 
                          call ogderiv(ep, 1,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, est(n)  )
                          call ogderiv(ep, 2,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, estt(n) )
                          call ogderiv(ep, 0,1,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esx(n) )
                          call ogderiv(ep, 0,0,1,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esy(n) )
                          call ogderiv(ep, 0,2,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esxx(n) )
                          call ogderiv(ep, 0,0,2,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esyy(n) )
                          call ogderiv(ep, 1,1,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, estx(n) )
                          call ogderiv(ep, 1,0,1,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esty(n) )
                          call ogderiv(ep, 2,1,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esttx(n) )
                          call ogderiv(ep, 2,0,1,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, estty(n) )
                          call ogderiv(ep, 0,3,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esxxx(n) )
                          call ogderiv(ep, 0,2,1,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esxxy(n) )
                          call ogderiv(ep, 0,1,2,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esxyy(n) )
                          call ogderiv(ep, 0,0,3,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esyyy(n) )
                          call ogderiv(ep, 3,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esttt(n) )
                          call ogderiv(ep, 4,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, estttt(n) )
                          call ogderiv(ep, 1,2,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, estxx(n) )
                          call ogderiv(ep, 1,0,2,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, estyy(n) )
                          call ogderiv(ep, 2,2,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esttxx(n) )
                          call ogderiv(ep, 2,0,2,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esttyy(n) )
                          call ogderiv(ep, 0,4,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esxxxx(n) )
                          call ogderiv(ep, 0,2,2,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esxxyy(n) )
                          call ogderiv(ep, 0,0,4,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esyyyy(n) )
                          ! L = c2^2*Delta
                          esL  = ( esxx(n)   + esyy(n) ) ! deleted c2^2
                          estL = ( estxx(n)  + estyy(n) )
                          esttL= ( esttxx(n) + esttyy(n) )
                          esLx  = ( esxxx(n)   + esxyy(n) )
                          esLy  = ( esxxy(n)   + esyyy(n) )
                          ! L^2 : 
                          esLL = ( esxxxx(n)  + 2.*esxxyy(n) + esyyyy(n) ) ! deleted c2^4
                          do jv=0,numberOfPolarizationVectors2-1
                            ! The TZ component is offset by pxc
                            pc = pxc + jv*nd
                            call ogderiv(ep, 0,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pe(n)   )
                            call ogderiv(ep, 1,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pet(n)  )
                            call ogderiv(ep, 2,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pett(n) )
                            call ogderiv(ep, 3,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pettt(n) )
                            call ogderiv(ep, 4,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, petttt(n) )
                            call ogderiv(ep, 0,1,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pex(n) )
                            call ogderiv(ep, 0,0,1,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pey(n) )
                            call ogderiv(ep, 0,2,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pexx(n) )
                            call ogderiv(ep, 0,0,2,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, peyy(n) )
                            call ogderiv(ep, 1,1,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, petx(n) )
                            call ogderiv(ep, 1,0,1,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pety(n) )
                            call ogderiv(ep, 1,2,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, petxx(n) )
                            call ogderiv(ep, 1,0,2,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, petyy(n) )
                            call ogderiv(ep, 2,1,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pettx(n) )
                            call ogderiv(ep, 2,0,1,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, petty(n) )
                            call ogderiv(ep, 2,2,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pettxx(n) )
                            call ogderiv(ep, 2,0,2,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pettyy(n) )
                            peL  = ( pexx(n)   + peyy(n) ) ! deleted c2^2
                            petL = ( petxx(n)  + petyy(n) )
                            pettL= ( pettxx(n) + pettyy(n) )
                            fpv2(n,jv) = pett(n)   + b1v2(jv)*pet(n)   + b0v2(jv)*pe(n)
                            fPt2(n,jv) = pettt(n)  + b1v2(jv)*pett(n)  + b0v2(jv)*pet(n)
                            fPtt2(n,jv)= petttt(n) + b1v2(jv)*pettt(n) + b0v2(jv)*pett(n)
                            LfP2(n,jv) = pettL     + b1v2(jv)*petL     + b0v2(jv)*peL
                            fpvx2(n,jv)= pettx(n)  + b1v2(jv)*petx(n)  + b0v2(jv)*pex(n)
                            fpvy2(n,jv)= petty(n)  + b1v2(jv)*pety(n)  + b0v2(jv)*pey(n)
                            do na = 0,numberOfAtomicLevels2-1
                              call ogderiv(ep, 0,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t, nce+na, q0  )
                              call ogderiv(ep, 1,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t, nce+na, q0t  )
                              call ogderiv(ep, 2,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t, nce+na, q0tt  )
                              call ogderiv(ep, 0,1,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t, nce+na, q0x  )
                              call ogderiv(ep, 0,0,1,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t, nce+na, q0y  )
                              call ogderiv(ep, 0,2,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t, nce+na, q0xx  )
                              call ogderiv(ep, 0,0,2,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t, nce+na, q0yy  )
                              fpv2(n,jv) = fpv2(n,jv) - pnec2(jv,na)*q0*es(n) ! adding \Delta N*E
                              fPt2(n,jv) = fPt2(n,jv) - pnec2(jv,na)*q0t*es(n) - pnec2(jv,na)*q0*est(n)
                              fPtt2(n,jv) = fPtt2(n,jv) - pnec2(jv,na)*q0tt*es(n)- 2.0*pnec2(jv,na)*q0t*est(n) - pnec2(jv,na)*q0*estt(n)
                              LfP2(n,jv) = LfP2(n,jv) - pnec2(jv,na)*(q0xx*es(n)+2.*q0x*esx(n)+q0*esxx(n) + q0yy*es(n)+2.*q0y*esy(n)+q0*esyy(n))
                              fpvx2(n,jv) = fpvx2(n,jv) - pnec2(jv,na)*q0x*es(n) - pnec2(jv,na)*q0*esx(n)
                              fpvy2(n,jv) = fpvy2(n,jv) - pnec2(jv,na)*q0y*es(n) - pnec2(jv,na)*q0*esy(n)
                            enddo
                            ! print *,'xxxxxxxxxxxxxxxxxxxxxxx'
                            ! print *, 'FOR P TZ FORCING'
                            ! print *, n,jv, fpv2(n,jv),fPt2(n,jv),fPtt2(n,jv),LfP2(n,jv),fpvx2(n,jv),fpvy2(n,jv)
                            ! print *,'xxxxxxxxxxxxxxxxxxxxxxx'
                            ! Normal TZ forcing for P_{n,jv} equation: 
                            ! fpv2(n,jv) = pett(n)   + b1v2(jv)*pet(n)   + b0v2(jv)*pe(n)   - a0v2(jv)*es(n)   - a1v2(jv)*est(n)
                            ! fPt2(n,jv) = pettt(n)  + b1v2(jv)*pett(n)  + b0v2(jv)*pet(n)  - a0v2(jv)*est(n)  - a1v2(jv)*estt(n)
                            ! fPtt2(n,jv)= petttt(n) + b1v2(jv)*pettt(n) + b0v2(jv)*pett(n) - a0v2(jv)*estt(n) - a1v2(jv)*esttt(n)
                            ! LfP2(n,jv) = pettL     + b1v2(jv)*petL     + b0v2(jv)*peL     - a0v2(jv)*esL     - a1v2(jv)*estL
                            ! fpvx2(n,jv)= pettx(n)  + b1v2(jv)*petx(n)  + b0v2(jv)*pex(n)  - a0v2(jv)*esx(n)  - a1v2(jv)*estx(n)
                            ! fpvy2(n,jv)= petty(n)  + b1v2(jv)*pety(n)  + b0v2(jv)*pey(n)  - a0v2(jv)*esy(n)  - a1v2(jv)*esty(n)
                            ! write(*,'(" n=",i4," LfP2=",e10.4," pettL,petL,peL,esL,estL=",5e12.4)') n,LfP2(n,jv),pettL,petL,peL,esL,estL
                            ! write(*,'(" pe,pet,pett,pettt,petttt=",5e12.4)') pe(n),pet(n),pett(n),pettt(n),petttt(n)
                            ! write(*,'("TZ: n,jv=",2i4," pex,pey,pexx,peyy=",4(1pe12.4))') n,jv,pex(n),pey(n),pexx(n),peyy(n)
                            ! Save ptt for checking later
                            pevtt2(n,jv)=pett(n)
                            pevttx2(n,jv)=pettx(n)
                            pevtty2(n,jv)=petty(n)
                            pevttt2(n,jv)=pettt(n)
                            pevtttt2(n,jv)=petttt(n)
                            pevttL2(n,jv) = pettL
                            pevttLSum2(n) = pevttLSum2(n)  + c2**2*pettL ! added c2**2 to be consistence with GDM for eval jumps
                            pevttttSum2(n)= pevttttSum2(n) + petttt(n) 
                            ! Keep some sums: 
                            fpSum2(n)   = fpSum2(n)  + fpv2(n,jv)
                            pevttSum2(n)  = pevttSum2(n)  + pett(n) 
                            pevttxSum2(n) = pevttxSum2(n) + pettx(n)
                            pevttySum2(n) = pevttySum2(n) + petty(n)
                            petttSum  = petttSum  + pettt(n) 
                          end do 
                          ! TZ forcing for E_{n} equation:
                          ! E_tt - c2^2 Delta E + alphaP2*Ptt  = 
                          fev2(n) = estt(n)   - c2**2*esL   + alphaP2*pevttSum2(n)
                          fEt2(n) = esttt(n)  - c2**2*estL  + alphaP2*petttSum
                          fEtt2(n)= estttt(n) - c2**2*esttL + alphaP2*pevttttSum2(n)
                          fevx2(n) = esttx(n) - c2**2*esLx   + alphaP2*pevttxSum2(n)
                          fevy2(n) = estty(n) - c2**2*esLy   + alphaP2*pevttySum2(n)
                          ! write(*,'("--> fEtt2=",e10.2," estttt,esttL,pettttSum=",3e10.2)')  fEtt2(n),estttt(n),esttL,pettttSum
                          LfE2(n) = esttL     - c2**2*esLL  + alphaP2*pevttLSum2(n)/c2**2 ! new: divided c2**2
                        end do
                        !--------------------------------
                        ! outside of dimension loop for N
                        !--------------------------------
                        do na=0,numberOfAtomicLevels2-1
                          ! na-th level
                          call ogderiv(ep, 1,0,0,0,xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t, nce+na, q0t )
                          call ogderiv(ep, 2,0,0,0,xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t, nce+na, q0tt)
                          call ogderiv(ep, 3,0,0,0,xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t, nce+na, q0ttt)
                          call ogderiv(ep, 4,0,0,0,xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t, nce+na, q0tttt)
                          ! initialize
                          fnv2(na)  = q0t ! forcing for \partial_tN_\ell = alpha_{\ell,k}N_k+\beta_{\ell,m}E\cdot\partial_tP_k
                          fntv2(na) = q0tt ! next derivative
                          fnttv2(na) = q0ttt
                          fntttv2(na) = q0tttt
                          ! relaxation (alpha_{\ell,m})
                          do jv=0,numberOfAtomicLevels2-1
                            call ogderiv(ep, 0,0,0,0,xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t, nce+jv, q0 )
                            call ogderiv(ep, 1,0,0,0,xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t, nce+jv, q0t)
                            call ogderiv(ep, 2,0,0,0,xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t, nce+jv, q0tt)
                            call ogderiv(ep, 3,0,0,0,xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t, nce+jv, q0ttt)
                            fnv2(na)  = fnv2(na)  - prc2(na,jv)*q0
                            fntv2(na) = fntv2(na) - prc2(na,jv)*q0t
                            fnttv2(na) = fnttv2(na) - prc2(na,jv)*q0tt
                            fntttv2(na) = fntttv2(na) - prc2(na,jv)*q0ttt
                          enddo
                          ! dot product (\beta_{\ell,k})
                          do n=0,nd-1 ! loop over dim
                            call ogderiv(ep, 0,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, es(n)   ) 
                            call ogderiv(ep, 1,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, est(n)  )
                            call ogderiv(ep, 2,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, estt(n) )
                            call ogderiv(ep, 3,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esttt(n) )
                            ! corresponding polarization vector
                            do jv=0,numberOfPolarizationVectors2-1 
                              pc = pxc + jv*nd
                              call ogderiv(ep, 0,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pe(n)   )
                              call ogderiv(ep, 1,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pet(n)  )
                              call ogderiv(ep, 2,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pett(n) )
                              call ogderiv(ep, 3,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pettt(n) )
                              call ogderiv(ep, 4,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, petttt(n) ) 
                              fnv2(na)  = fnv2(na) - peptc2(na,jv)*es(n)*pet(n)
                              fntv2(na) = fntv2(na) - peptc2(na,jv)*est(n)*pet(n) - peptc2(na,jv)*es(n)*pett(n)
                              fnttv2(na) = fnttv2(na) - peptc2(na,jv)*estt(n)*pet(n) - 2.d0*peptc2(na,jv)*est(n)*pett(n) - peptc2(na,jv)*es(n)*pettt(n)
                              fntttv2(na) = fntttv2(na) - peptc2(na,jv)*esttt(n)*pet(n) - 3.d0*peptc2(na,jv)*estt(n)*pett(n) - 3.d0*peptc2(na,jv)*est(n)*pettt(n) - peptc2(na,jv)*es(n)*petttt(n)
                            enddo
                          enddo
                        enddo
                      !-----------------------
                      ! GDM
                      !-----------------------
                      elseif( dispersionModel2.ne.noDispersion .and. nonlinearModel2.eq.noNonlinearModel) then
                        do n=0,nd-1
                          fpSum2(n)=0.
                          pevttSum2(n)=0.
                          pevttxSum2(n)=0.
                          pevttySum2(n)=0.
                          pevttLSum2(n)=0.
                          pevttttSum2(n)=0.
                          petttSum=0.
                          call ogderiv(ep, 0,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, es(n)   ) 
                          call ogderiv(ep, 1,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, est(n)  )
                          call ogderiv(ep, 2,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, estt(n) )
                          call ogderiv(ep, 0,1,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esx(n) )
                          call ogderiv(ep, 0,0,1,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esy(n) )
                          call ogderiv(ep, 0,2,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esxx(n) )
                          call ogderiv(ep, 0,0,2,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esyy(n) )
                          call ogderiv(ep, 1,1,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, estx(n) )
                          call ogderiv(ep, 1,0,1,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esty(n) )
                          call ogderiv(ep, 2,1,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esttx(n) )
                          call ogderiv(ep, 2,0,1,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, estty(n) )
                          call ogderiv(ep, 0,3,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esxxx(n) )
                          call ogderiv(ep, 0,2,1,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esxxy(n) )
                          call ogderiv(ep, 0,1,2,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esxyy(n) )
                          call ogderiv(ep, 0,0,3,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esyyy(n) )
                          call ogderiv(ep, 3,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esttt(n) )
                          call ogderiv(ep, 4,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, estttt(n) )
                          call ogderiv(ep, 1,2,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, estxx(n) )
                          call ogderiv(ep, 1,0,2,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, estyy(n) )
                          call ogderiv(ep, 2,2,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esttxx(n) )
                          call ogderiv(ep, 2,0,2,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esttyy(n) )
                          call ogderiv(ep, 0,4,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esxxxx(n) )
                          call ogderiv(ep, 0,2,2,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esxxyy(n) )
                          call ogderiv(ep, 0,0,4,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esyyyy(n) )
                          ! L = c2^2*Delta
                          esL  = c2**2*( esxx(n)   + esyy(n) )
                          estL = c2**2*( estxx(n)  + estyy(n) )
                          esttL= c2**2*( esttxx(n) + esttyy(n) )
                          esLx  = c2**2*( esxxx(n)   + esxyy(n) )
                          esLy  = c2**2*( esxxy(n)   + esyyy(n) )
                          ! L^2 : 
                          esLL = c2**4*( esxxxx(n)  + 2.*esxxyy(n) + esyyyy(n) )
                          do jv=0,numberOfPolarizationVectors2-1
                            ! The TZ component is offset by pxc
                            pc = pxc + jv*nd
                            call ogderiv(ep, 0,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pe(n)   )
                            call ogderiv(ep, 1,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pet(n)  )
                            call ogderiv(ep, 2,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pett(n) )
                            call ogderiv(ep, 3,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pettt(n) )
                            call ogderiv(ep, 4,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, petttt(n) )
                            call ogderiv(ep, 0,1,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pex(n) )
                            call ogderiv(ep, 0,0,1,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pey(n) )
                            call ogderiv(ep, 0,2,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pexx(n) )
                            call ogderiv(ep, 0,0,2,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, peyy(n) )
                            call ogderiv(ep, 1,1,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, petx(n) )
                            call ogderiv(ep, 1,0,1,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pety(n) )
                            call ogderiv(ep, 1,2,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, petxx(n) )
                            call ogderiv(ep, 1,0,2,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, petyy(n) )
                            call ogderiv(ep, 2,1,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pettx(n) )
                            call ogderiv(ep, 2,0,1,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, petty(n) )
                            call ogderiv(ep, 2,2,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pettxx(n) )
                            call ogderiv(ep, 2,0,2,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pettyy(n) )
                            peL  = c2**2*( pexx(n)   + peyy(n) )
                            petL = c2**2*( petxx(n)  + petyy(n) )
                            pettL= c2**2*( pettxx(n) + pettyy(n) )
                            ! Normal TZ forcing for P_{n,jv} equation: 
                            fpv2(n,jv) = pett(n)   + b1v2(jv)*pet(n)   + b0v2(jv)*pe(n)   - a0v2(jv)*es(n)   - a1v2(jv)*est(n)
                            fPt2(n,jv) = pettt(n)  + b1v2(jv)*pett(n)  + b0v2(jv)*pet(n)  - a0v2(jv)*est(n)  - a1v2(jv)*estt(n)
                            fPtt2(n,jv)= petttt(n) + b1v2(jv)*pettt(n) + b0v2(jv)*pett(n) - a0v2(jv)*estt(n) - a1v2(jv)*esttt(n)
                            LfP2(n,jv) = pettL     + b1v2(jv)*petL     + b0v2(jv)*peL     - a0v2(jv)*esL     - a1v2(jv)*estL
                            fpvx2(n,jv)= pettx(n)  + b1v2(jv)*petx(n)  + b0v2(jv)*pex(n)  - a0v2(jv)*esx(n)  - a1v2(jv)*estx(n)
                            fpvy2(n,jv)= petty(n)  + b1v2(jv)*pety(n)  + b0v2(jv)*pey(n)  - a0v2(jv)*esy(n)  - a1v2(jv)*esty(n)
                            ! write(*,'(" n=",i4," LfP2=",e10.4," pettL,petL,peL,esL,estL=",5e12.4)') n,LfP2(n,jv),pettL,petL,peL,esL,estL
                            ! write(*,'(" pe,pet,pett,pettt,petttt=",5e12.4)') pe(n),pet(n),pett(n),pettt(n),petttt(n)
                            ! write(*,'("TZ: n,jv=",2i4," pex,pey,pexx,peyy=",4(1pe12.4))') n,jv,pex(n),pey(n),pexx(n),peyy(n)
                            ! Save ptt for checking later
                            pevtt2(n,jv)=pett(n)
                            pevttx2(n,jv)=pettx(n)
                            pevtty2(n,jv)=petty(n)
                            pevtttt2(n,jv)=petttt(n)
                            pevttLSum2(n) = pevttLSum2(n)  + pettL
                            pevttttSum2(n)= pevttttSum2(n) + petttt(n) 
                            ! Keep some sums: 
                            fpSum2(n)   = fpSum2(n)  + fpv2(n,jv)
                            pevttSum2(n)  = pevttSum2(n)  + pett(n) 
                            pevttxSum2(n) = pevttxSum2(n) + pettx(n)
                            pevttySum2(n) = pevttySum2(n) + petty(n)
                            petttSum  = petttSum  + pettt(n) 
                          end do 
                          ! TZ forcing for E_{n} equation:
                          ! E_tt - c2^2 Delta E + alphaP2*Ptt  = 
                          fev2(n) = estt(n)   - esL   + alphaP2*pevttSum2(n)
                          fEt2(n) = esttt(n)  - estL  + alphaP2*petttSum
                          fEtt2(n)= estttt(n) - esttL + alphaP2*pevttttSum2(n)
                          fevx2(n) = esttx(n) - esLx   + alphaP2*pevttxSum2(n)
                          fevy2(n) = estty(n) - esLy   + alphaP2*pevttySum2(n)
                          ! write(*,'("--> fEtt2=",e10.2," estttt,esttL,pettttSum=",3e10.2)')  fEtt2(n),estttt(n),esttL,pettttSum
                          LfE2(n) = esttL     - esLL  + alphaP2*pevttLSum2(n)
                       end do
                      !-----------------------
                      ! no dispersion
                      !-----------------------
                      elseif( dispersionModel2.eq.noDispersion .and. nonlinearModel2.eq.noNonlinearModel) then
                        do n=0,nd-1
                          fpSum2(n)=0.
                          pevttSum2(n)=0.
                          pevttxSum2(n)=0.
                          pevttySum2(n)=0.
                          pevttLSum2(n)=0.
                          pevttttSum2(n)=0.
                        enddo
                      end if
                    end if
                   ! These derivatives are computed to 2nd-order accuracy
                   ! NOTE: the jacobian derivatives can be computed once for all components
                     uu1=u1(i1,i2,i3,ex) ! in the rectangular case just eval the solution
                      u1xxx = (-u1(i1-2,i2,i3,ex)+2.*u1(i1-1,i2,i3,ex)-2.*u1(i1+1,i2,i3,ex)+u1(i1+2,i2,i3,ex))/(2.*dx1(0)**3)
                      u1xxy = ((-u1(i1-1,i2-1,i3,ex)+u1(i1-1,i2+1,i3,ex))/(2.*dx1(1))-2.*(-u1(i1,i2-1,i3,ex)+u1(i1,i2+1,i3,ex))/(2.*dx1(1))+(-u1(i1+1,i2-1,i3,ex)+u1(i1+1,i2+1,i3,ex))/(2.*dx1(1)))/(dx1(0)**2)
                      u1xyy = (-(u1(i1-1,i2-1,i3,ex)-2.*u1(i1-1,i2,i3,ex)+u1(i1-1,i2+1,i3,ex))/(dx1(1)**2)+(u1(i1+1,i2-1,i3,ex)-2.*u1(i1+1,i2,i3,ex)+u1(i1+1,i2+1,i3,ex))/(dx1(1)**2))/(2.*dx1(0))
                      u1yyy = (-u1(i1,i2-2,i3,ex)+2.*u1(i1,i2-1,i3,ex)-2.*u1(i1,i2+1,i3,ex)+u1(i1,i2+2,i3,ex))/(2.*dx1(1)**3)
                      u1xxxx = (u1(i1-2,i2,i3,ex)-4.*u1(i1-1,i2,i3,ex)+6.*u1(i1,i2,i3,ex)-4.*u1(i1+1,i2,i3,ex)+u1(i1+2,i2,i3,ex))/(dx1(0)**4)
                      u1xxyy = ((u1(i1-1,i2-1,i3,ex)-2.*u1(i1-1,i2,i3,ex)+u1(i1-1,i2+1,i3,ex))/(dx1(1)**2)-2.*(u1(i1,i2-1,i3,ex)-2.*u1(i1,i2,i3,ex)+u1(i1,i2+1,i3,ex))/(dx1(1)**2)+(u1(i1+1,i2-1,i3,ex)-2.*u1(i1+1,i2,i3,ex)+u1(i1+1,i2+1,i3,ex))/(dx1(1)**2))/(dx1(0)**2)
                      u1yyyy = (u1(i1,i2-2,i3,ex)-4.*u1(i1,i2-1,i3,ex)+6.*u1(i1,i2,i3,ex)-4.*u1(i1,i2+1,i3,ex)+u1(i1,i2+2,i3,ex))/(dx1(1)**4)
                    u1LapSq = u1xxxx +2.* u1xxyy + u1yyyy
                     vv1=u1(i1,i2,i3,ey) ! in the rectangular case just eval the solution
                      v1xxx = (-u1(i1-2,i2,i3,ey)+2.*u1(i1-1,i2,i3,ey)-2.*u1(i1+1,i2,i3,ey)+u1(i1+2,i2,i3,ey))/(2.*dx1(0)**3)
                      v1xxy = ((-u1(i1-1,i2-1,i3,ey)+u1(i1-1,i2+1,i3,ey))/(2.*dx1(1))-2.*(-u1(i1,i2-1,i3,ey)+u1(i1,i2+1,i3,ey))/(2.*dx1(1))+(-u1(i1+1,i2-1,i3,ey)+u1(i1+1,i2+1,i3,ey))/(2.*dx1(1)))/(dx1(0)**2)
                      v1xyy = (-(u1(i1-1,i2-1,i3,ey)-2.*u1(i1-1,i2,i3,ey)+u1(i1-1,i2+1,i3,ey))/(dx1(1)**2)+(u1(i1+1,i2-1,i3,ey)-2.*u1(i1+1,i2,i3,ey)+u1(i1+1,i2+1,i3,ey))/(dx1(1)**2))/(2.*dx1(0))
                      v1yyy = (-u1(i1,i2-2,i3,ey)+2.*u1(i1,i2-1,i3,ey)-2.*u1(i1,i2+1,i3,ey)+u1(i1,i2+2,i3,ey))/(2.*dx1(1)**3)
                      v1xxxx = (u1(i1-2,i2,i3,ey)-4.*u1(i1-1,i2,i3,ey)+6.*u1(i1,i2,i3,ey)-4.*u1(i1+1,i2,i3,ey)+u1(i1+2,i2,i3,ey))/(dx1(0)**4)
                      v1xxyy = ((u1(i1-1,i2-1,i3,ey)-2.*u1(i1-1,i2,i3,ey)+u1(i1-1,i2+1,i3,ey))/(dx1(1)**2)-2.*(u1(i1,i2-1,i3,ey)-2.*u1(i1,i2,i3,ey)+u1(i1,i2+1,i3,ey))/(dx1(1)**2)+(u1(i1+1,i2-1,i3,ey)-2.*u1(i1+1,i2,i3,ey)+u1(i1+1,i2+1,i3,ey))/(dx1(1)**2))/(dx1(0)**2)
                      v1yyyy = (u1(i1,i2-2,i3,ey)-4.*u1(i1,i2-1,i3,ey)+6.*u1(i1,i2,i3,ey)-4.*u1(i1,i2+1,i3,ey)+u1(i1,i2+2,i3,ey))/(dx1(1)**4)
                    v1LapSq = v1xxxx +2.* v1xxyy + v1yyyy
                   ! NOTE: the jacobian derivatives can be computed once for all components
                     uu2=u2(j1,j2,j3,ex) ! in the rectangular case just eval the solution
                      u2xxx = (-u2(j1-2,j2,j3,ex)+2.*u2(j1-1,j2,j3,ex)-2.*u2(j1+1,j2,j3,ex)+u2(j1+2,j2,j3,ex))/(2.*dx2(0)**3)
                      u2xxy = ((-u2(j1-1,j2-1,j3,ex)+u2(j1-1,j2+1,j3,ex))/(2.*dx2(1))-2.*(-u2(j1,j2-1,j3,ex)+u2(j1,j2+1,j3,ex))/(2.*dx2(1))+(-u2(j1+1,j2-1,j3,ex)+u2(j1+1,j2+1,j3,ex))/(2.*dx2(1)))/(dx2(0)**2)
                      u2xyy = (-(u2(j1-1,j2-1,j3,ex)-2.*u2(j1-1,j2,j3,ex)+u2(j1-1,j2+1,j3,ex))/(dx2(1)**2)+(u2(j1+1,j2-1,j3,ex)-2.*u2(j1+1,j2,j3,ex)+u2(j1+1,j2+1,j3,ex))/(dx2(1)**2))/(2.*dx2(0))
                      u2yyy = (-u2(j1,j2-2,j3,ex)+2.*u2(j1,j2-1,j3,ex)-2.*u2(j1,j2+1,j3,ex)+u2(j1,j2+2,j3,ex))/(2.*dx2(1)**3)
                      u2xxxx = (u2(j1-2,j2,j3,ex)-4.*u2(j1-1,j2,j3,ex)+6.*u2(j1,j2,j3,ex)-4.*u2(j1+1,j2,j3,ex)+u2(j1+2,j2,j3,ex))/(dx2(0)**4)
                      u2xxyy = ((u2(j1-1,j2-1,j3,ex)-2.*u2(j1-1,j2,j3,ex)+u2(j1-1,j2+1,j3,ex))/(dx2(1)**2)-2.*(u2(j1,j2-1,j3,ex)-2.*u2(j1,j2,j3,ex)+u2(j1,j2+1,j3,ex))/(dx2(1)**2)+(u2(j1+1,j2-1,j3,ex)-2.*u2(j1+1,j2,j3,ex)+u2(j1+1,j2+1,j3,ex))/(dx2(1)**2))/(dx2(0)**2)
                      u2yyyy = (u2(j1,j2-2,j3,ex)-4.*u2(j1,j2-1,j3,ex)+6.*u2(j1,j2,j3,ex)-4.*u2(j1,j2+1,j3,ex)+u2(j1,j2+2,j3,ex))/(dx2(1)**4)
                    u2LapSq = u2xxxx +2.* u2xxyy + u2yyyy
                     vv2=u2(j1,j2,j3,ey) ! in the rectangular case just eval the solution
                      v2xxx = (-u2(j1-2,j2,j3,ey)+2.*u2(j1-1,j2,j3,ey)-2.*u2(j1+1,j2,j3,ey)+u2(j1+2,j2,j3,ey))/(2.*dx2(0)**3)
                      v2xxy = ((-u2(j1-1,j2-1,j3,ey)+u2(j1-1,j2+1,j3,ey))/(2.*dx2(1))-2.*(-u2(j1,j2-1,j3,ey)+u2(j1,j2+1,j3,ey))/(2.*dx2(1))+(-u2(j1+1,j2-1,j3,ey)+u2(j1+1,j2+1,j3,ey))/(2.*dx2(1)))/(dx2(0)**2)
                      v2xyy = (-(u2(j1-1,j2-1,j3,ey)-2.*u2(j1-1,j2,j3,ey)+u2(j1-1,j2+1,j3,ey))/(dx2(1)**2)+(u2(j1+1,j2-1,j3,ey)-2.*u2(j1+1,j2,j3,ey)+u2(j1+1,j2+1,j3,ey))/(dx2(1)**2))/(2.*dx2(0))
                      v2yyy = (-u2(j1,j2-2,j3,ey)+2.*u2(j1,j2-1,j3,ey)-2.*u2(j1,j2+1,j3,ey)+u2(j1,j2+2,j3,ey))/(2.*dx2(1)**3)
                      v2xxxx = (u2(j1-2,j2,j3,ey)-4.*u2(j1-1,j2,j3,ey)+6.*u2(j1,j2,j3,ey)-4.*u2(j1+1,j2,j3,ey)+u2(j1+2,j2,j3,ey))/(dx2(0)**4)
                      v2xxyy = ((u2(j1-1,j2-1,j3,ey)-2.*u2(j1-1,j2,j3,ey)+u2(j1-1,j2+1,j3,ey))/(dx2(1)**2)-2.*(u2(j1,j2-1,j3,ey)-2.*u2(j1,j2,j3,ey)+u2(j1,j2+1,j3,ey))/(dx2(1)**2)+(u2(j1+1,j2-1,j3,ey)-2.*u2(j1+1,j2,j3,ey)+u2(j1+1,j2+1,j3,ey))/(dx2(1)**2))/(dx2(0)**2)
                      v2yyyy = (u2(j1,j2-2,j3,ey)-4.*u2(j1,j2-1,j3,ey)+6.*u2(j1,j2,j3,ey)-4.*u2(j1,j2+1,j3,ey)+u2(j1,j2+2,j3,ey))/(dx2(1)**4)
                    v2LapSq = v2xxxx +2.* v2xxyy + v2yyyy
                   ! These derivatives are computed to 4th-order accuracy
                   ! NOTE: the jacobian derivatives can be computed once for all components
                     uu1=u1(i1,i2,i3,ex) ! in the rectangular case just eval the solution
                      u1x = (u1(i1-2,i2,i3,ex)-8.*u1(i1-1,i2,i3,ex)+8.*u1(i1+1,i2,i3,ex)-u1(i1+2,i2,i3,ex))/(12.*dx1(0))
                      u1y = (u1(i1,i2-2,i3,ex)-8.*u1(i1,i2-1,i3,ex)+8.*u1(i1,i2+1,i3,ex)-u1(i1,i2+2,i3,ex))/(12.*dx1(1))
                      u1xx = (-u1(i1-2,i2,i3,ex)+16.*u1(i1-1,i2,i3,ex)-30.*u1(i1,i2,i3,ex)+16.*u1(i1+1,i2,i3,ex)-u1(i1+2,i2,i3,ex))/(12.*dx1(0)**2)
                      u1yy = (-u1(i1,i2-2,i3,ex)+16.*u1(i1,i2-1,i3,ex)-30.*u1(i1,i2,i3,ex)+16.*u1(i1,i2+1,i3,ex)-u1(i1,i2+2,i3,ex))/(12.*dx1(1)**2)
                    u1Lap = u1xx+ u1yy
                     vv1=u1(i1,i2,i3,ey) ! in the rectangular case just eval the solution
                      v1x = (u1(i1-2,i2,i3,ey)-8.*u1(i1-1,i2,i3,ey)+8.*u1(i1+1,i2,i3,ey)-u1(i1+2,i2,i3,ey))/(12.*dx1(0))
                      v1y = (u1(i1,i2-2,i3,ey)-8.*u1(i1,i2-1,i3,ey)+8.*u1(i1,i2+1,i3,ey)-u1(i1,i2+2,i3,ey))/(12.*dx1(1))
                      v1xx = (-u1(i1-2,i2,i3,ey)+16.*u1(i1-1,i2,i3,ey)-30.*u1(i1,i2,i3,ey)+16.*u1(i1+1,i2,i3,ey)-u1(i1+2,i2,i3,ey))/(12.*dx1(0)**2)
                      v1yy = (-u1(i1,i2-2,i3,ey)+16.*u1(i1,i2-1,i3,ey)-30.*u1(i1,i2,i3,ey)+16.*u1(i1,i2+1,i3,ey)-u1(i1,i2+2,i3,ey))/(12.*dx1(1)**2)
                    v1Lap = v1xx+ v1yy
                   ! NOTE: the jacobian derivatives can be computed once for all components
                     uu2=u2(j1,j2,j3,ex) ! in the rectangular case just eval the solution
                      u2x = (u2(j1-2,j2,j3,ex)-8.*u2(j1-1,j2,j3,ex)+8.*u2(j1+1,j2,j3,ex)-u2(j1+2,j2,j3,ex))/(12.*dx2(0))
                      u2y = (u2(j1,j2-2,j3,ex)-8.*u2(j1,j2-1,j3,ex)+8.*u2(j1,j2+1,j3,ex)-u2(j1,j2+2,j3,ex))/(12.*dx2(1))
                      u2xx = (-u2(j1-2,j2,j3,ex)+16.*u2(j1-1,j2,j3,ex)-30.*u2(j1,j2,j3,ex)+16.*u2(j1+1,j2,j3,ex)-u2(j1+2,j2,j3,ex))/(12.*dx2(0)**2)
                      u2yy = (-u2(j1,j2-2,j3,ex)+16.*u2(j1,j2-1,j3,ex)-30.*u2(j1,j2,j3,ex)+16.*u2(j1,j2+1,j3,ex)-u2(j1,j2+2,j3,ex))/(12.*dx2(1)**2)
                    u2Lap = u2xx+ u2yy
                     vv2=u2(j1,j2,j3,ey) ! in the rectangular case just eval the solution
                      v2x = (u2(j1-2,j2,j3,ey)-8.*u2(j1-1,j2,j3,ey)+8.*u2(j1+1,j2,j3,ey)-u2(j1+2,j2,j3,ey))/(12.*dx2(0))
                      v2y = (u2(j1,j2-2,j3,ey)-8.*u2(j1,j2-1,j3,ey)+8.*u2(j1,j2+1,j3,ey)-u2(j1,j2+2,j3,ey))/(12.*dx2(1))
                      v2xx = (-u2(j1-2,j2,j3,ey)+16.*u2(j1-1,j2,j3,ey)-30.*u2(j1,j2,j3,ey)+16.*u2(j1+1,j2,j3,ey)-u2(j1+2,j2,j3,ey))/(12.*dx2(0)**2)
                      v2yy = (-u2(j1,j2-2,j3,ey)+16.*u2(j1,j2-1,j3,ey)-30.*u2(j1,j2,j3,ey)+16.*u2(j1,j2+1,j3,ey)-u2(j1,j2+2,j3,ey))/(12.*dx2(1)**2)
                    v2Lap = v2xx+ v2yy
                   ! Store c^2*Delta(E) in a vector 
                   LE1(0)=(c1**2)*u1Lap
                   LE1(1)=(c1**2)*v1Lap
                   LE2(0)=(c2**2)*u2Lap
                   LE2(1)=(c2**2)*v2Lap
                   ! Store L^2(E) 
                   LLE1(0)=(c1**4)*u1LapSq
                   LLE1(1)=(c1**4)*v1LapSq
                   LLE2(0)=(c2**4)*u2LapSq
                   LLE2(1)=(c2**4)*v2LapSq
                   ! Store (LE).x an (LE).y 
                   LEx1(0) = (c1**2)*( u1xxx + u1xyy )
                   LEx1(1) = (c1**2)*( v1xxx + v1xyy )
                   LEy1(0) = (c1**2)*( u1xxy + u1yyy )
                   LEy1(1) = (c1**2)*( v1xxy + v1yyy )
                   LEx2(0) = (c2**2)*( u2xxx + u2xyy )
                   LEx2(1) = (c2**2)*( v2xxx + v2xyy )
                   LEy2(0) = (c2**2)*( u2xxy + u2yyy )
                   LEy2(1) = (c2**2)*( v2xxy + v2yyy )
                   ! We also need derivatives at the old time:
                   ! These next derivatives may only be needed to order2, but use order 4 for now so exact for degree 4
                   ! #perl $ORDER=2;
                     uu1=u1n(i1,i2,i3,ex) ! in the rectangular case just eval the solution
                      u1nx = (u1n(i1-2,i2,i3,ex)-8.*u1n(i1-1,i2,i3,ex)+8.*u1n(i1+1,i2,i3,ex)-u1n(i1+2,i2,i3,ex))/(12.*dx1(0))
                      u1ny = (u1n(i1,i2-2,i3,ex)-8.*u1n(i1,i2-1,i3,ex)+8.*u1n(i1,i2+1,i3,ex)-u1n(i1,i2+2,i3,ex))/(12.*dx1(1))
                      u1nxx = (-u1n(i1-2,i2,i3,ex)+16.*u1n(i1-1,i2,i3,ex)-30.*u1n(i1,i2,i3,ex)+16.*u1n(i1+1,i2,i3,ex)-u1n(i1+2,i2,i3,ex))/(12.*dx1(0)**2)
                      u1nyy = (-u1n(i1,i2-2,i3,ex)+16.*u1n(i1,i2-1,i3,ex)-30.*u1n(i1,i2,i3,ex)+16.*u1n(i1,i2+1,i3,ex)-u1n(i1,i2+2,i3,ex))/(12.*dx1(1)**2)
                    u1nLap = u1nxx+ u1nyy
                     vv1=u1n(i1,i2,i3,ey) ! in the rectangular case just eval the solution
                      v1nx = (u1n(i1-2,i2,i3,ey)-8.*u1n(i1-1,i2,i3,ey)+8.*u1n(i1+1,i2,i3,ey)-u1n(i1+2,i2,i3,ey))/(12.*dx1(0))
                      v1ny = (u1n(i1,i2-2,i3,ey)-8.*u1n(i1,i2-1,i3,ey)+8.*u1n(i1,i2+1,i3,ey)-u1n(i1,i2+2,i3,ey))/(12.*dx1(1))
                      v1nxx = (-u1n(i1-2,i2,i3,ey)+16.*u1n(i1-1,i2,i3,ey)-30.*u1n(i1,i2,i3,ey)+16.*u1n(i1+1,i2,i3,ey)-u1n(i1+2,i2,i3,ey))/(12.*dx1(0)**2)
                      v1nyy = (-u1n(i1,i2-2,i3,ey)+16.*u1n(i1,i2-1,i3,ey)-30.*u1n(i1,i2,i3,ey)+16.*u1n(i1,i2+1,i3,ey)-u1n(i1,i2+2,i3,ey))/(12.*dx1(1)**2)
                    v1nLap = v1nxx+ v1nyy
                   ! Here are c^2*Delta(E) at the old time: 
                   LE1m(0) = (c1**2)*u1nLap
                   LE1m(1) = (c1**2)*v1nLap
                     uu2=u2n(j1,j2,j3,ex) ! in the rectangular case just eval the solution
                      u2nx = (u2n(j1-2,j2,j3,ex)-8.*u2n(j1-1,j2,j3,ex)+8.*u2n(j1+1,j2,j3,ex)-u2n(j1+2,j2,j3,ex))/(12.*dx2(0))
                      u2ny = (u2n(j1,j2-2,j3,ex)-8.*u2n(j1,j2-1,j3,ex)+8.*u2n(j1,j2+1,j3,ex)-u2n(j1,j2+2,j3,ex))/(12.*dx2(1))
                      u2nxx = (-u2n(j1-2,j2,j3,ex)+16.*u2n(j1-1,j2,j3,ex)-30.*u2n(j1,j2,j3,ex)+16.*u2n(j1+1,j2,j3,ex)-u2n(j1+2,j2,j3,ex))/(12.*dx2(0)**2)
                      u2nyy = (-u2n(j1,j2-2,j3,ex)+16.*u2n(j1,j2-1,j3,ex)-30.*u2n(j1,j2,j3,ex)+16.*u2n(j1,j2+1,j3,ex)-u2n(j1,j2+2,j3,ex))/(12.*dx2(1)**2)
                    u2nLap = u2nxx+ u2nyy
                     vv2=u2n(j1,j2,j3,ey) ! in the rectangular case just eval the solution
                      v2nx = (u2n(j1-2,j2,j3,ey)-8.*u2n(j1-1,j2,j3,ey)+8.*u2n(j1+1,j2,j3,ey)-u2n(j1+2,j2,j3,ey))/(12.*dx2(0))
                      v2ny = (u2n(j1,j2-2,j3,ey)-8.*u2n(j1,j2-1,j3,ey)+8.*u2n(j1,j2+1,j3,ey)-u2n(j1,j2+2,j3,ey))/(12.*dx2(1))
                      v2nxx = (-u2n(j1-2,j2,j3,ey)+16.*u2n(j1-1,j2,j3,ey)-30.*u2n(j1,j2,j3,ey)+16.*u2n(j1+1,j2,j3,ey)-u2n(j1+2,j2,j3,ey))/(12.*dx2(0)**2)
                      v2nyy = (-u2n(j1,j2-2,j3,ey)+16.*u2n(j1,j2-1,j3,ey)-30.*u2n(j1,j2,j3,ey)+16.*u2n(j1,j2+1,j3,ey)-u2n(j1,j2+2,j3,ey))/(12.*dx2(1)**2)
                    v2nLap = v2nxx+ v2nyy
                   LE2m(0) = (c2**2)*u2nLap
                   LE2m(1) = (c2**2)*v2nLap
                   evx1(0) = u1x
                   evx1(1) = v1x
                   evy1(0) = u1y
                   evy1(1) = v1y 
                   evnx1(0) = u1nx
                   evnx1(1) = v1nx
                   evny1(0) = u1ny
                   evny1(1) = v1ny 
                   evx2(0) = u2x
                   evx2(1) = v2x
                   evy2(0) = u2y
                   evy2(1) = v2y 
                   evnx2(0) = u2nx
                   evnx2(1) = v2nx
                   evny2(0) = u2ny
                   evny2(1) = v2ny 
                  ! eval nonlinear dispersive forcings for domain 1
                  ! getDispersiveForcingOrder4(LEFT,i1,i2,i3, fp1, fpv1,fev1, p1,p1n,p1m, u1,u1n,u1m, dispersionModel1,numberOfPolarizationVectors1,alphaP1,!             c2PttEsum1,c2PttLEsum1,c4PttLEsum1,c4PttLLEsum1,c2PttttLEsum1,c2PttttLLEsum1,a0v1,a1v1,!             b0v1,b1v1,LE1,LLE1,LE1m,LfE1,LfP1,fEt1,fEtt1,fPt1,fPtt1,pevtt1,pevttx1,pevtty1,pevtttt1,evx1,evy1,evnx1,evny1,!             fevx1,fevy1,fpvx1,fpvy1,LEx1,LEy1,fPttx1,fPtty1,fLPtt1,fPtttt1)
                    ! pre-assign 0 values
                    do n=0,nd-1 ! dispersive forcing in jump conditions
                      fp1(n)=0.
                      fPttx1(n) =0.
                      fPtty1(n) =0.
                      fLPtt1(n) =0.
                      fPtttt1(n)=0.
                    end do
                    !-----------------------------
                    ! MLA
                    !-----------------------------
                    ! only do this for MLA (dispersive and nonlinear multi-level)
                    if( dispersionModel1.ne.noDispersion .and. nonlinearModel1.ne.noNonlinearModel) then
                      nce = pxc+nd*numberOfPolarizationVectors1
                      ! -----------------------------------------
                      ! order 2 (E, P, N) at the interface (fictitious step)
                      !------------------------------------------
                      ! dimension loop for E and P
                      do n=0,nd-1
                        ec = ex +n 
                        ev0  =  u1n(i1,i2,i3,ec)
                        ev   =  u1(i1,i2,i3,ec) ! time where we need to fill in ghost points
                        pSum=0.
                        do jv=0,numberOfPolarizationVectors1-1
                          pc = n + jv*nd  
                          pv0 =  p1n(i1,i2,i3,pc)
                          pv  =  p1(i1,i2,i3,pc)
                          pvn = 2.*pv-p1n(i1,i2,i3,pc) + 0.5*dt*b1v1(jv)*p1n(i1,i2,i3,pc) - dtsq*b0v1(jv)*pv + dtsq*fpv1(n,jv)
                          do na = 0,numberOfAtomicLevels1-1 ! \Delta N^n*E^n
                            pvn = pvn + dtsq*pnec1(jv,na)*q1(i1,i2,i3,na)*ev
                          enddo ! na
                          pvec(n,jv)= pvn/( 1.+.5*dt*b1v1(jv) ) ! time + dt
                          ! #If "p1" eq "p1"
                          ! call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t+dt,pxc+jv*nd+n, pe(n)   )
                          ! ! call ogderiv(ep, 1,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pet(n)   )
                          ! ! call ogderiv(ep, 2,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pett(n)   )
                          ! #Else
                          ! call ogderiv(ep, 0,0,0,0, xy2(i1,i2,i3,0),xy2(i1,i2,i3,1),0.,t+dt,pxc+jv*nd+n, pe(n)   )
                          ! ! call ogderiv(ep, 1,0,0,0, xy2(i1,i2,i3,0),xy2(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pet(n)   )
                          ! ! call ogderiv(ep, 2,0,0,0, xy2(i1,i2,i3,0),xy2(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pett(n)   )
                          ! #End
                          ! print *, '---------Dispersive forcing 2----------'
                          ! print *, pvec(n,jv),pe(n),pvec(n,jv)-pe(n)
                          pSum = pSum + pvec(n,jv) -2.*pv + pv0 ! keep sum
                          ! print *, '++++++Dispersive forcing+++++++++++++++'
                          ! print *, 'check 2nd output P: ', pc,jv,pvec(n,jv)
                        enddo ! jv
                        ! second order update of E
                        evec(n) = (2.*ev-ev0) + cSq1*dtsq*LE1(n)/cSq1 - alphaP1*pSum + dtsq*fev1(n) ! cSq1 is already in LE1
                        ! print *, '++++++Dispersive forcing+++++++++++++++'
                        ! print *, ev,ev0, cSq1, dtsq,LE1(n),alphaP1,pSum,fev1(n)
                        ! print *, 'check 2nd output E: ', ec,evec(n)
                      enddo ! n
                      ! N outside of space loop
                      ! 1st derivative
                      do na=0,numberOfAtomicLevels1-1
                        qt(na) = fnv1(na)
                        do jv = 0,numberOfAtomicLevels1-1
                          qt(na) = qt(na)+prc1(na,jv)*q1(i1,i2,i3,jv)
                        enddo
                        do n=0,nd-1
                          do jv=0,numberOfPolarizationVectors1-1
                            qt(na) = qt(na) + peptc1(na,jv)*u1(i1,i2,i3,ex+n)*(pvec(n,jv)-p1n(i1,i2,i3,n+jv*nd))/(2.*dt)
                          enddo
                        enddo
                      enddo
                      ! 2nd derivative
                      do na=0,numberOfAtomicLevels1-1
                        qtt(na) = fntv1(na)
                        do jv = 0,numberOfAtomicLevels1-1
                          qtt(na) = qtt(na)+prc1(na,jv)*qt(jv)
                        enddo
                        do n=0,nd-1
                          do jv=0,numberOfPolarizationVectors1-1
                            qtt(na) = qtt(na) + peptc1(na,jv)*(evec(n)-u1n(i1,i2,i3,ex+n))/(2.*dt)*(pvec(n,jv)-p1n(i1,i2,i3,n+jv*nd))/(2.*dt)+ peptc1(na,jv)*u1(i1,i2,i3,ex+n)*(pvec(n,jv)-2.*p1(i1,i2,i3,n+jv*nd)+p1n(i1,i2,i3,n+jv*nd))/(dtsq)
                          enddo
                        enddo
                      enddo
                      ! taylor expansion
                      do na=0,numberOfAtomicLevels1-1
                        qv(na) = q1(i1,i2,i3,na) + dt*qt(na) + dtsq/2.*qtt(na)
                      enddo
                      !----------------------------------------
                      ! order 4 update of P at interface (fictitious step)
                      !----------------------------------------
                      ! second order accurate terms
                      do n=0,nd-1
                        do jv = 0,numberOfPolarizationVectors1-1
                          ! #If "p1" eq "p1"
                          ! call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pe(n)   )
                          ! call ogderiv(ep, 1,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pet(n)   )
                          ! call ogderiv(ep, 2,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pett(n)   )
                          ! #Else
                          ! call ogderiv(ep, 0,0,0,0, xy2(i1,i2,i3,0),xy2(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pe(n)   )
                          ! call ogderiv(ep, 1,0,0,0, xy2(i1,i2,i3,0),xy2(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pet(n)   )
                          ! call ogderiv(ep, 2,0,0,0, xy2(i1,i2,i3,0),xy2(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pett(n)   )
                          ! #End
                          ! pvec(n,jv) = pe(n)
                          ! ptv(n,jv) = pet(n)
                          ! pttv(n,jv) = pett(n)
                          ptv(n,jv) = (pvec(n,jv)-p1n(i1,i2,i3,n+jv*nd))/(2.*dt)
                          pttv(n,jv) = (pvec(n,jv)-2.*p1(i1,i2,i3,n+jv*nd)+p1n(i1,i2,i3,n+jv*nd))/dtsq
                          ptttv(n,jv) = -b1v1(jv)*pttv(n,jv)-b0v1(jv)*ptv(n,jv)+fPt1(n,jv)
                          do na = 0,numberOfAtomicLevels1-1 ! update using ODE
                            ptttv(n,jv) = ptttv(n,jv) + pnec1(jv,na)*qt(na)*u1(i1,i2,i3,ex+n) + pnec1(jv,na)*q1(i1,i2,i3,na)*(evec(n)-u1n(i1,i2,i3,ex+n))/(2.*dt)
                          enddo
                          pttttv(n,jv) = -b1v1(jv)*ptttv(n,jv)-b0v1(jv)*pttv(n,jv)+fPtt1(n,jv) ! update using ODE
                          do na = 0,numberOfAtomicLevels1-1
                            pttttv(n,jv) = pttttv(n,jv) + pnec1(jv,na)*qtt(na)*u1(i1,i2,i3,ex+n) + 2.*pnec1(jv,na)*qt(na)*(evec(n)-u1n(i1,i2,i3,ex+n))/(2.*dt) + pnec1(jv,na)*q1(i1,i2,i3,na)*(evec(n)-2.*u1(i1,i2,i3,ex+n)+u1n(i1,i2,i3,ex+n))/dtsq
                            ! print *, '++++++Dispersive forcing+++++++++++++++'
                            ! print *, 'check P derivatives: ', n,jv,ptv(n,jv),pttv(n,jv),ptttv(n,jv),pttttv(n,jv)
                            ! print *, na,pnec1(jv,na),q1(i1,i2,i3,na),qt(na),qtt(na)
                          enddo
                          ! print *, '++++++Dispersive forcing+++++++++++++++'
                          ! print *, 'check P time derivatives: ', n,jv,ptv(n,jv),pttv(n,jv),ptttv(n,jv),pttttv(n,jv)
                          ! print *, fPt1(n,jv),fPtt1(n,jv)
                        enddo
                        ! print *, evec(n),u1(i1,i2,i3,ex+n),u1n(i1,i2,i3,ex+n)
                      enddo
                      ! dimension loop for E and P
                      do n=0,nd-1
                        ec = ex +n 
                        ev   =  u1(i1,i2,i3,ec) ! time where we need to fill in ghost points
                        ! These next derivatives may only be needed to order2, (use order 4 for testing TZ polynomials)
                        ! #perl $ORDER=2; ! should use order 2 since E is only filled at first ghost lines at t
                          ! uu1 in the next statement defines names of intermediate values
                            uu1=u1(i1,i2,i3,ec) ! in the rectangular case just eval the solution
                             e1x = (u1(i1-2,i2,i3,ec)-8.*u1(i1-1,i2,i3,ec)+8.*u1(i1+1,i2,i3,ec)-u1(i1+2,i2,i3,ec))/(12.*dx1(0))
                             e1y = (u1(i1,i2-2,i3,ec)-8.*u1(i1,i2-1,i3,ec)+8.*u1(i1,i2+1,i3,ec)-u1(i1,i2+2,i3,ec))/(12.*dx1(1))
                             e1xx = (-u1(i1-2,i2,i3,ec)+16.*u1(i1-1,i2,i3,ec)-30.*u1(i1,i2,i3,ec)+16.*u1(i1+1,i2,i3,ec)-u1(i1+2,i2,i3,ec))/(12.*dx1(0)**2)
                             e1yy = (-u1(i1,i2-2,i3,ec)+16.*u1(i1,i2-1,i3,ec)-30.*u1(i1,i2,i3,ec)+16.*u1(i1,i2+1,i3,ec)-u1(i1,i2+2,i3,ec))/(12.*dx1(1)**2)
                           e1Lap = e1xx+ e1yy
                          ! write(*,'("LEFT: u1x,u1y,u1xx,u1yy,u1Lap=",5(1pe12.4))') u1x,u1y,u1xx,u1yy,u1Lap
                          evx0  = e1x
                          evy0  = e1y
                          evLap = e1xx+e1yy ! these values use the second order predicted values in the first ghost lines
                            ! print *, '++++++Dispersive forcing+++++++++++++++'
                            ! print *, 'check E spatial derivatives: ', n,evx0,evy0,evLap
                        do jv=0,numberOfPolarizationVectors1-1
                          pc = n + jv*nd 
                          ! These next derivatives may only be needed to order2, (use order 4 for testing TZ polynomials)
                          ! #perl $ORDER=2;
                            ! uu1 in the next statement defines names of intermediate values
                              pp1=p1(i1,i2,i3,pc) ! in the rectangular case just eval the solution
                               p1x = (p1(i1-2,i2,i3,pc)-8.*p1(i1-1,i2,i3,pc)+8.*p1(i1+1,i2,i3,pc)-p1(i1+2,i2,i3,pc))/(12.*dx1(0))
                               p1y = (p1(i1,i2-2,i3,pc)-8.*p1(i1,i2-1,i3,pc)+8.*p1(i1,i2+1,i3,pc)-p1(i1,i2+2,i3,pc))/(12.*dx1(1))
                               p1xx = (-p1(i1-2,i2,i3,pc)+16.*p1(i1-1,i2,i3,pc)-30.*p1(i1,i2,i3,pc)+16.*p1(i1+1,i2,i3,pc)-p1(i1+2,i2,i3,pc))/(12.*dx1(0)**2)
                               p1yy = (-p1(i1,i2-2,i3,pc)+16.*p1(i1,i2-1,i3,pc)-30.*p1(i1,i2,i3,pc)+16.*p1(i1,i2+1,i3,pc)-p1(i1,i2+2,i3,pc))/(12.*dx1(1)**2)
                             p1Lap = p1xx+ p1yy
                            ! write(*,'("LEFT: p1x,p1y,p1xx,p1yy,p1Lap=",5(1pe12.4))') p1x,p1y,p1xx,p1yy,p1Lap
                            LP  = p1Lap ! removed c^2
                              pp1=p1n(i1,i2,i3,pc) ! in the rectangular case just eval the solution
                               p1nx = (p1n(i1-2,i2,i3,pc)-8.*p1n(i1-1,i2,i3,pc)+8.*p1n(i1+1,i2,i3,pc)-p1n(i1+2,i2,i3,pc))/(12.*dx1(0))
                               p1ny = (p1n(i1,i2-2,i3,pc)-8.*p1n(i1,i2-1,i3,pc)+8.*p1n(i1,i2+1,i3,pc)-p1n(i1,i2+2,i3,pc))/(12.*dx1(1))
                               p1nxx = (-p1n(i1-2,i2,i3,pc)+16.*p1n(i1-1,i2,i3,pc)-30.*p1n(i1,i2,i3,pc)+16.*p1n(i1+1,i2,i3,pc)-p1n(i1+2,i2,i3,pc))/(12.*dx1(0)**2)
                               p1nyy = (-p1n(i1,i2-2,i3,pc)+16.*p1n(i1,i2-1,i3,pc)-30.*p1n(i1,i2,i3,pc)+16.*p1n(i1,i2+1,i3,pc)-p1n(i1,i2+2,i3,pc))/(12.*dx1(1)**2)
                             p1nLap = p1nxx+ p1nyy
                            ! write(*,'("LEFT: p1nxx,p1nyy,p1nLap=",3e12.4)') p1nxx,p1nyy,p1nLap
                            LPm = p1nLap
                            pvx  = p1x
                            pvy  = p1y
                            pvnx = p1nx
                            pvny = p1ny
                          ! print *, '++++++Dispersive forcing+++++++++++++++'
                          ! print *, 'check P spatial derivatives: ', n,jv,LP,LPm,pvx,pvy,pvnx,pvny
                          pv0 =  p1n(i1,i2,i3,pc)
                          pv  =  p1(i1,i2,i3,pc)
                          pvn = 2.*pv-p1n(i1,i2,i3,pc) + 0.5*dt*b1v1(jv)*p1n(i1,i2,i3,pc) + dt**4/12.*pttttv(n,jv) + dt**4/6.*b1v1(jv)*ptttv(n,jv) - dtsq*b0v1(jv)*pv + dtsq*fpv1(n,jv)
                          ! pvn = 2.*pv-p1n(i1,i2,i3,pc) + 0.5*dt*b1v1(jv)*p1n(i1,i2,i3,pc) - dtsq*b0v1(jv)*pv + dtsq*fpv1(n,jv)
                          do na = 0,numberOfAtomicLevels1-1 ! \Delta N^n*E^n
                            pvn = pvn + dtsq*pnec1(jv,na)*q1(i1,i2,i3,na)*ev
                          enddo ! na
                          pvec(n,jv)= pvn/( 1.+.5*dt*b1v1(jv) ) ! time + dt
                          ! 4th order accurate term
                          fp1(n) = fp1(n) + (pvec(n,jv)-2.*p1(i1,i2,i3,pc)+p1n(i1,i2,i3,pc))/dtsq  - dt**2/12.*pttttv(n,jv)
                          ! #If "p1" eq "p1"
                          !   call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t+dt,pxc+jv*nd+n, pe(n)   )
                          !   call ogderiv(ep, 1,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pet(n)   )
                          !   call ogderiv(ep, 2,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pett(n)   )
                          !   call ogderiv(ep, 4,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pxc+jv*nd+n, petttt(n)   )
                          ! #Else
                          !   call ogderiv(ep, 0,0,0,0, xy2(i1,i2,i3,0),xy2(i1,i2,i3,1),0.,t+dt,pxc+jv*nd+n, pe(n)   )
                          !   call ogderiv(ep, 1,0,0,0, xy2(i1,i2,i3,0),xy2(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pet(n)   )
                          !   call ogderiv(ep, 2,0,0,0, xy2(i1,i2,i3,0),xy2(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pett(n)   )
                          !   call ogderiv(ep, 4,0,0,0, xy2(i1,i2,i3,0),xy2(i1,i2,i3,1),0.,t,pxc+jv*nd+n, petttt(n)   )
                          ! #End
                          ! print *, '---------Dispersive forcing 4----------'
                          ! print *, pvec(n,jv),pe(n),pvec(n,jv)-pe(n)
                          ! print *, (pvec(n,jv)-2.*p1(i1,i2,i3,pc)+p1n(i1,i2,i3,pc))/dtsq  - dt**2/12.*pttttv(n,jv), pett(n),(pvec(n,jv)-2.*p1(i1,i2,i3,pc)+p1n(i1,i2,i3,pc))/dtsq  - dt**2/12.*pttttv(n,jv)-pett(n)
                          ! fp1(n) = fp1(n) + pett(n)
                          ! 2nd order accurate terms
                          fPtttt1(n) = fPtttt1(n) + pttttv(n,jv)
                          ! fPtttt1(n) = fPtttt1(n) + petttt(n)
                          ! print *, '++++++Dispersive forcing+++++++++++++++'
                          ! print *, 'check P time derivatives 2 and 4: ', n,jv,(pvec(n,jv)-2.*p1(i1,i2,i3,pc)+p1n(i1,i2,i3,pc))/dtsq  - dt**4/12.*pttttv(n,jv),pttttv(n,jv)
                          !--------------------------------
                          ! spatial derivatives
                          !--------------------------------
                          ! nce = pxc+nd*numberOfPolarizationVectors1
                          ! N*E
                          ! #perl $ORDER=2;
                          do na=0,numberOfAtomicLevels1-1
                            ! uu1 in the next statement defines names of intermediate values
                              qq1=q1(i1,i2,i3,na) ! in the rectangular case just eval the solution
                               q1x = (q1(i1-2,i2,i3,na)-8.*q1(i1-1,i2,i3,na)+8.*q1(i1+1,i2,i3,na)-q1(i1+2,i2,i3,na))/(12.*dx1(0))
                               q1y = (q1(i1,i2-2,i3,na)-8.*q1(i1,i2-1,i3,na)+8.*q1(i1,i2+1,i3,na)-q1(i1,i2+2,i3,na))/(12.*dx1(1))
                               q1xx = (-q1(i1-2,i2,i3,na)+16.*q1(i1-1,i2,i3,na)-30.*q1(i1,i2,i3,na)+16.*q1(i1+1,i2,i3,na)-q1(i1+2,i2,i3,na))/(12.*dx1(0)**2)
                               q1yy = (-q1(i1,i2-2,i3,na)+16.*q1(i1,i2-1,i3,na)-30.*q1(i1,i2,i3,na)+16.*q1(i1,i2+1,i3,na)-q1(i1,i2+2,i3,na))/(12.*dx1(1)**2)
                             q1Lap = q1xx+ q1yy
                            ! write(*,'("LEFT: p1x,p1y,p1xx,p1yy,p1Lap=",5(1pe12.4))') p1x,p1y,p1xx,p1yy,p1Lap
                            qvx  = q1x
                            qvy  = q1y
                            qvLap  = q1Lap
                            qex(na) = evx0*q1(i1,i2,i3,na)+qvx*ev
                            qey(na) = evy0*q1(i1,i2,i3,na)+qvy*ev
                            qeLap(na) = ev*qvLap+q1(i1,i2,i3,na)*evLap+2.*evx0*qvx+2.*evy0*qvy
                            ! print *, '++++++Dispersive forcing+++++++++++++++'
                            ! print *, 'check N spatial derivatives: ', na,qex(na),qey(na),qeLap(na)
                          enddo
                          ! #If "p1" eq "p1"
                          !   call ogderiv(ep, 2,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pevttxx(n,jv)   )
                          !   call ogderiv(ep, 2,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pevttyy(n,jv)   )
                          !   call ogderiv(ep, 2,1,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pevttx1(n,jv)   )
                          !   call ogderiv(ep, 2,0,1,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pevtty1(n,jv)   )
                          ! #Else
                          !   call ogderiv(ep, 2,2,0,0, xy2(i1,i2,i3,0),xy2(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pevttxx(n,jv)   )
                          !   call ogderiv(ep, 2,0,2,0, xy2(i1,i2,i3,0),xy2(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pevttyy(n,jv)   )
                          !   call ogderiv(ep, 2,1,0,0, xy2(i1,i2,i3,0),xy2(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pevttx1(n,jv)   )
                          !   call ogderiv(ep, 2,0,1,0, xy2(i1,i2,i3,0),xy2(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pevtty1(n,jv)   )
                          ! #End
                          ! laplacian
                          LPn = 2.*LP-LPm + 0.5*dt*b1v1(jv)*LPm - dtsq*b0v1(jv)*LP + dtsq*LfP1(n,jv)
                          do na = 0,numberOfAtomicLevels1-1 ! \Delta N^n*E^n
                            LPn = LPn + dtsq*pnec1(jv,na)*qeLap(na)
                          enddo
                          ! time derivatives
                          fLPtt1(n) = fLPtt1(n) + cSq1*(LPn/(1.+.5*dt*b1v1(jv)) - 2.*LP + LPm)/dtsq ! added c^2 to be consistence with GDM version
                          ! fLPtt1(n) = fLPtt1(n) + pevttL1(n,jv)
                          ! print *,'error of laplacian with # of levels',numberOfAtomicLevels1,numberOfPolarizationVectors1,pevttxx(n,jv)+pevttyy(n,jv)-pevttL1(n,jv)
                          ! print *, '++++++Dispersive forcing+++++++++++++++'
                          ! print *, 'check P lap: ', n,jv,LP,LPm,b1v1(jv),dtsq,b0v1(jv),LfP1(n,jv),fLPtt1(n)
                          ! x
                          pttxa = 2.*pvx-pvnx + 0.5*dt*b1v1(jv)*pvnx - dtsq*b0v1(jv)*pvx + dtsq*fpvx1(n,jv)
                          do na = 0,numberOfAtomicLevels1-1 ! \Delta N^n*E^n
                            pttxa = pttxa + dtsq*pnec1(jv,na)*qex(na)
                          enddo
                          ! time derivatives
                          fPttx1(n) = fPttx1(n) + (pttxa/(1.+.5*dt*b1v1(jv)) - 2.*pvx + pvnx)/dtsq
                          ! fPttx1(n) = fPttx1(n) + pevttx1(n,jv)
                          ! print *, '++++++Dispersive forcing+++++++++++++++'
                          ! print *, 'check Pttx: ', n,jv,pvx,pvnx,fpvx1(n,jv),fPttx1(n)
                          ! y
                          pttya = 2.*pvy-pvny + 0.5*dt*b1v1(jv)*pvny - dtsq*b0v1(jv)*pvy + dtsq*fpvy1(n,jv)
                          do na = 0,numberOfAtomicLevels1-1 ! \Delta N^n*E^n
                            pttya = pttya + dtsq*pnec1(jv,na)*qey(na)
                          enddo
                          ! time derivatives
                          fPtty1(n) = fPtty1(n) + (pttya/(1.+.5*dt*b1v1(jv)) - 2.*pvy + pvny)/dtsq
                          ! fPtty1(n) = fPtty1(n) + pevtty1(n,jv)
                          ! print *, '++++++Dispersive forcing+++++++++++++++'
                          ! print *, 'check Ptty: ', n,jv,n,jv,pvy,pvny,fpvy1(n,jv),fPtty1(n)
                          ! if( twilightZone.eq.1 )then
                          !   write(*,'("")')
                          !   write(*,'("DI4:LEFT: i1,i2=",2i3," jv=",i2," n=",i2)') i1,i2,jv,n
                          !   print *, 'ptt diff',(pvec(n,jv)-2.*p1(i1,i2,i3,pc)+p1n(i1,i2,i3,pc))/dtsq - dt**2/12.*pttttv(n,jv)-pevtt1(n,jv)
                          !   print *, 'pttx diff', (pttxa/(1.+.5*dt*b1v1(jv)) - 2.*pvx + pvnx)/dtsq-pevttx1(n,jv)
                          !   print *, 'ptty diff', (pttya/(1.+.5*dt*b1v1(jv)) - 2.*pvy + pvny)/dtsq-pevtty1(n,jv)
                          !   print *, 'ptttt diff', pttttv(n,jv)-pevtttt1(n,jv)
                          !   print *, 'pttL diff',(LPn/(1.+.5*dt*b1v1(jv)) - 2.*LP + LPm)/dtsq-pevttL1(n,jv)
                          ! end if
                        enddo ! jv
                      enddo ! n
                    !-----------------------------
                    ! GDM
                    !-----------------------------
                    elseif( dispersionModel1.ne.noDispersion .and. nonlinearModel1.eq.noNonlinearModel) then
                      c2PttEsum1=0.
                      c2PttLEsum1=0.
                      c4PttLEsum1=0.
                      c4PttLLEsum1=0.
                      c2PttttLEsum1=0.
                      c2PttttLLEsum1=0.
                      do jv=0,numberOfPolarizationVectors1-1
                        a0=a0v1(jv)
                        a1=a1v1(jv)
                        b0=b0v1(jv)
                        b1=b1v1(jv)
                        alpha=alphaP1
                        ! Second-order coefficients: 
                        ! Ptt = c2PttLE*LE1 + c2PttE*E + c2PttEm*Em + c2PttP*P + c2PttPm*Pm + c2PttfE*fE + c2PttfP*fP
      ! File created by Dropbox/GDM/maple/interface.maple
      d1=2.+(a1*alpha+b1)*dt

      ! --------------- Here is Et to second-order ------------

      ! Et = -1/2/dt*(2*a0*alpha*dt^2-2*b1*dt-4)/d1*E-1/2/dt*(2*b1*dt+4)/d1*Em-1/2/dt*(-b1*dt^3-2*dt^2)/d1*LE-1/2/dt*(-2*alpha*b0*dt^2-2*alpha*b1*dt)/d1*P-alpha*b1/d1*Pm-1/2/dt*(-b1*dt^3-2*dt^2)/d1*fE-dt*alpha*fP/d1
      ! Et = c2EtLE*LE + c2EtE*E + c2EtEm*Em + c2EtP*P + c2EtPm*Pm + c2EtfE*fE + c2EtfP*fP
      !   LE = c^2Delta(E) 
      c2EtLE=(.5*b1*dt**2+1.0*dt)/d1
      c2EtE=(-a0*alpha*dt**2+b1*dt+2.)/dt/d1
      c2EtEm=(-b1*dt-2.)/dt/d1
      c2EtP=alpha*(b0*dt+b1)/d1
      c2EtPm=-1.*alpha*b1/d1
      c2EtfE=(.5*b1*dt**2+1.0*dt)/d1
      c2EtfP=-1.*dt*alpha/d1
      ! --------------- Here is Pt to second-order ------------

      ! Pt = 1/2*(2*a0*dt^2+2*a1*dt)/dt/d1*E-a1/d1*Em+1/2*a1*dt^2/d1*LE+1/2*(2*a1*alpha*dt-2*b0*dt^2+4)/dt/d1*P+1/2*(-2*a1*alpha*dt-4)/dt/d1*Pm+1/2*a1*dt^2/d1*fE+fP*dt/d1
      ! Pt = c2PtLE*LE + c2PtE*E + c2PtEm*Em + c2PtP*P + c2PtPm*Pm + c2PtfE*fE + c2PtfP*fP
      !   LE = c^2Delta(E) 
      c2PtLE=.5*a1*dt**2/d1
      c2PtE=(1.*a0*dt+1.*a1)/d1
      c2PtEm=-1.*a1/d1
      c2PtP=(alpha*a1*dt-b0*dt**2+2.)/dt/d1
      c2PtPm=(-alpha*a1*dt-2.)/dt/d1
      c2PtfE=.5*a1*dt**2/d1
      c2PtfP=1.*dt/d1
      ! --------------- Here is Ptt to second-order ------------

      ! Ptt = a1*dt/d1*LE+((-a0*a1*alpha-a0*b1)*dt^2+d1*a0*dt+2*a1)/dt/d1*E-2*a1/dt/d1*Em+((a1*alpha*b0+b0*b1)*dt^2-d1*b0*dt-2*b1)/dt/d1*P+2*b1/dt/d1*Pm+a1*dt/d1*fE+((-a1*alpha-b1)*dt^2+dt*d1)/dt/d1*fP
      ! Ptt = c2PttLE*LE + c2PttE*E + c2PttEm*Em + c2PttP*P + c2PttPm*Pm + c2PttfE*fE + c2PttfP*fP
      c2PttLE=1.*a1*dt/d1
      c2PttE=((-1.*a1*alpha-1.*b1)*a0*dt**2+d1*a0*dt+2.*a1)/dt/d1
      c2PttEm=-2.*a1/dt/d1
      c2PttP=(b0*(a1*alpha+b1)*dt**2-1.*d1*b0*dt-2.*b1)/dt/d1
      c2PttPm=2.*b1/dt/d1
      c2PttfE=1.*a1*dt/d1
      c2PttfP=(1.*d1+(-1.*a1*alpha-1.*b1)*dt)/d1
                        ! Fourth-order coefficients
                        ! Ptt = c4PttLE*LE1 + c4PttE*E + c4PttEm*Em + c4PttP*P + c4PttPm*Pm + c4PttfE*fE + c4PttfP*fP
                        !    + c4PttLLE*LLE1 + c4PttLP*LP + c4PttLEm*LE1m + c4PttLPm*LPm+ c4PttLfE*LfE1 + c4PttLfP*LfP1
                        !    + c4PttfEt*fEt1 + c4PttfEtt*fEtt1 + c4PttfPt*fPt1+ c4PttfPtt*fPtt1
                        ! #Include interfaceAdeGdmOrder4.h 
                        ! Nov 4, 2018 *new way* 
      ! File created by Dropbox/GDM/maple/interfaceNew.maple
      ! 
      ! --------------- Here is Ptt to fourth-order ------------

      ! Ptt = ((1/2*(-b0*c2PttE*alpha*dt^6-6*c2PttE*alpha*b1*dt^5+(12*a1^2*c2PttLE*alpha^3+24*a1*(c2PttLE*b1-1/2*a1+1/2*c2PtLE*b0-1/2*a0*c2EtLE)*alpha^2+(12*c2PttLE*b1^2+(-12*a0*c2EtLE+12*b0*c2PtLE-12*a1)*b1+12*c2EtE*a1-12*c2PttE)*alpha+12*b0)*dt^4+72*b1*dt^3+144*dt^2)/d4/dt-1/6*dt^2*(-(a0*c2EtLE+(-alpha*c2PttLE+1)*a1-c2PtLE*b0-c2PttLE*b1)*alpha+c2EtE))*a1-(1/2*(-a0*c2PttE*alpha*dt^6-6*c2PttE*alpha*a1*dt^5+(-12*(alpha*c2PttLE-1)*alpha*a1^2+((12*a0*c2EtLE-12*b0*c2PtLE-24*b1*c2PttLE)*alpha+12*b1-12*c2EtE)*a1-12*c2PttLE*b1^2+(12*a0*c2EtLE-12*b0*c2PtLE)*b1+12*a0)*dt^4+72*a1*dt^3)/d4/dt-1/6*dt^2*(a0*c2EtLE+(-alpha*c2PttLE+1)*a1-c2PtLE*b0-c2PttLE*b1))*b1)*LE+((1/2*(-b0*(alpha*c2PttLE-1)*dt^6-6*(alpha*c2PttLE-1)*b1*dt^5+((12*a1*c2EtLE-12*c2PttLE)*alpha+12)*dt^4)/d4/dt-1/6*dt^2*c2EtLE)*a1-1/2*(-a0*(alpha*c2PttLE-1)*dt^6-6*(alpha*c2PttLE-1)*a1*dt^5-12*c2EtLE*a1*dt^4)/d4/dt*b1)*LLE+((1/2*(-b0*c2PttEm*alpha*dt^6-6*c2PttEm*alpha*b1*dt^5+(12*a1*c2EtEm-12*c2PttEm)*alpha*dt^4)/d4/dt-1/6*dt^2*c2EtEm)*a1-1/2*(-a0*alpha*c2PttEm*dt^6-6*a1*alpha*c2PttEm*dt^5-12*a1*c2EtEm*dt^4)/d4/dt*b1)*LEm+((1/2*(-b0*c2PttP*alpha*dt^6-6*c2PttP*alpha*b1*dt^5+(12*a1*c2EtP-12*c2PttP)*alpha*dt^4)/d4/dt-1/6*dt^2*c2EtP)*a1-1/2*(-a0*alpha*c2PttP*dt^6-6*a1*alpha*c2PttP*dt^5-12*a1*c2EtP*dt^4)/d4/dt*b1)*LP+((1/2*(-b0*c2PttPm*alpha*dt^6-6*c2PttPm*alpha*b1*dt^5+(12*a1*c2EtPm-12*c2PttPm)*alpha*dt^4)/d4/dt-1/6*dt^2*c2EtPm)*a1-1/2*(-a0*alpha*c2PttPm*dt^6-6*a1*alpha*c2PttPm*dt^5-12*a1*c2EtPm*dt^4)/d4/dt*b1)*LPm+(a0+(1/2*((12*a1^2*c2PttE*alpha^3+24*a1*(c2PttE*b1+1/2*c2PtE*b0-1/2*c2EtE*a0)*alpha^2+(12*c2PttE*b1^2+(-12*a0*c2EtE+12*b0*c2PtE)*b1)*alpha)*dt^4+(-120*a0*alpha+24*b0)*dt^2+144*b1*dt+288)/d4/dt+1/6*dt^2*(-a1*alpha*c2PttE+a0*c2EtE-b0*c2PtE-b1*c2PttE)*alpha)*a1-(1/2*((-12*c2PttE*alpha^2*a1^2+(12*a0*c2EtE-12*b0*c2PtE-24*b1*c2PttE)*alpha*a1-12*c2PttE*b1^2+(12*a0*c2EtE-12*b0*c2PtE)*b1)*dt^4+144*a0*dt^2+144*a1*dt)/d4/dt-1/6*dt^2*(-a1*alpha*c2PttE+a0*c2EtE-b0*c2PtE-b1*c2PttE))*b1)*E+((1/2*(((12*a1^2*c2PttEm*alpha^3+24*a1*(c2PttEm*b1+1/2*c2PtEm*b0-1/2*c2EtEm*a0)*alpha^2+(12*c2PttEm*b1^2+(-12*a0*c2EtEm+12*b0*c2PtEm)*b1)*alpha)*dt^4+(-12*a0*alpha-12*b0)*dt^2+(72*a1*alpha-72*b1)*dt-144)/d4-1)/dt+1/6*dt^2*(-a1*alpha*c2PttEm+a0*c2EtEm-b0*c2PtEm-b1*c2PttEm)*alpha)*a1-(1/2*((-12*c2PttEm*alpha^2*a1^2+(12*a0*c2EtEm-12*b0*c2PtEm-24*b1*c2PttEm)*alpha*a1-12*c2PttEm*b1^2+(12*a0*c2EtEm-12*b0*c2PtEm)*b1)*dt^4-144*a1*dt)/d4/dt-1/6*dt^2*(-a1*alpha*c2PttEm+a0*c2EtEm-b0*c2PtEm-b1*c2PttEm))*b1)*Em+((1/2*((12*a1^2*c2PttP*alpha^3+24*a1*(c2PttP*b1+1/2*c2PtP*b0-1/2*c2EtP*a0)*alpha^2+(12*c2PttP*b1^2+(-12*a0*c2EtP+12*b0*c2PtP)*b1)*alpha)*dt^4+144*b0*alpha*dt^2+144*b1*alpha*dt)/d4/dt+1/6*dt^2*(-a1*alpha*c2PttP+a0*c2EtP-b0*c2PtP-b1*c2PttP)*alpha)*a1-b0-(1/2*((-12*c2PttP*alpha^2*a1^2+(12*a0*c2EtP-12*b0*c2PtP-24*b1*c2PttP)*alpha*a1-12*c2PttP*b1^2+(12*a0*c2EtP-12*b0*c2PtP)*b1)*dt^4+(24*a0*alpha-120*b0)*dt^2+144*alpha*a1*dt+288)/d4/dt-1/6*dt^2*(-a1*alpha*c2PttP+a0*c2EtP-b0*c2PtP-b1*c2PttP))*b1)*P+((1/2*((12*a1^2*c2PttPm*alpha^3+24*a1*(c2PttPm*b1+1/2*c2PtPm*b0-1/2*c2EtPm*a0)*alpha^2+(12*c2PttPm*b1^2+(-12*a0*c2EtPm+12*b0*c2PtPm)*b1)*alpha)*dt^4-144*b1*alpha*dt)/d4/dt+1/6*dt^2*(-a1*alpha*c2PttPm+a0*c2EtPm-b0*c2PtPm-b1*c2PttPm)*alpha)*a1-(1/2*(((-12*c2PttPm*alpha^2*a1^2+(12*a0*c2EtPm-12*b0*c2PtPm-24*b1*c2PttPm)*alpha*a1-12*c2PttPm*b1^2+(12*a0*c2EtPm-12*b0*c2PtPm)*b1)*dt^4+(-12*a0*alpha-12*b0)*dt^2+(-72*a1*alpha+72*b1)*dt-144)/d4-1)/dt-1/6*dt^2*(-a1*alpha*c2PttPm+a0*c2EtPm-b0*c2PtPm-b1*c2PttPm))*b1)*Pm+((1/2*((12*a1^2*c2PttfE*alpha^3+24*a1*(c2PttfE*b1-1/2*a1+1/2*c2PtfE*b0-1/2*c2EtfE*a0)*alpha^2+(12*c2PttfE*b1^2+(-12*a0*c2EtfE+12*b0*c2PtfE-12*a1)*b1)*alpha+12*b0)*dt^4+72*b1*dt^3+144*dt^2)/d4/dt+1/6*dt^2*(c2EtfE*a0+(-alpha*c2PttfE+1)*a1-c2PtfE*b0-c2PttfE*b1)*alpha)*a1-(1/2*((-12*(alpha*c2PttfE-1)*alpha*a1^2+((12*a0*c2EtfE-12*b0*c2PtfE-24*b1*c2PttfE)*alpha+12*b1)*a1-12*c2PttfE*b1^2+(12*a0*c2EtfE-12*b0*c2PtfE)*b1+12*a0)*dt^4+72*a1*dt^3)/d4/dt-1/6*dt^2*(c2EtfE*a0+(-alpha*c2PttfE+1)*a1-c2PtfE*b0-c2PttfE*b1))*b1)*fE+((1/2*(-b0*(alpha*c2PttfE-1)*dt^6-6*(alpha*c2PttfE-1)*b1*dt^5+((12*a1*c2EtfE-12*c2PttfE)*alpha+12)*dt^4)/d4/dt-1/6*dt^2*c2EtfE)*a1-1/2*(-a0*(alpha*c2PttfE-1)*dt^6-6*(alpha*c2PttfE-1)*a1*dt^5-12*c2EtfE*a1*dt^4)/d4/dt*b1)*LfE+((6*alpha*a1*dt^3/d4-1/6*dt^2)*a1+6*a1*dt^3/d4*b1)*fEt+(1/2*(b0*dt^6+6*b1*dt^5+12*dt^4)/d4/dt*a1-1/2*(a0*dt^6+6*a1*dt^5)/d4/dt*b1)*fEtt+((1/2*((12*a1^2*alpha^3*c2PttfP+24*a1*(b1*c2PttfP+1/2*b0*c2PtfP-1/2*a0*c2EtfP)*alpha^2+(12*b1^2*c2PttfP+(-12*a0*c2EtfP+12*b0*c2PtfP)*b1)*alpha)*dt^4-144*alpha*dt^2)/d4/dt+1/6*dt^2*(-a1*alpha*c2PttfP+a0*c2EtfP-b0*c2PtfP-b1*c2PttfP)*alpha)*a1-(1/2*((-12*a1^2*alpha^2*c2PttfP+(12*a0*c2EtfP-12*b0*c2PtfP-24*b1*c2PttfP)*alpha*a1-12*b1^2*c2PttfP+(12*a0*c2EtfP-12*b0*c2PtfP)*b1)*dt^4+144*dt^2)/d4/dt-1/6*dt^2*(-a1*alpha*c2PttfP+a0*c2EtfP-b0*c2PtfP-b1*c2PttfP))*b1+1)*fP+((1/2*(-alpha*b0*c2PttfP*dt^6-6*alpha*b1*c2PttfP*dt^5+(12*a1*c2EtfP-12*c2PttfP)*alpha*dt^4)/d4/dt-1/6*dt^2*c2EtfP)*a1-1/2*(-a0*alpha*c2PttfP*dt^6-6*a1*alpha*c2PttfP*dt^5-12*a1*c2EtfP*dt^4)/d4/dt*b1)*LfP+((1/2*(-12*a1*alpha^2-12*alpha*b1)*dt^3/d4+1/6*alpha*dt^2)*a1-(1/2*(12*a1*alpha+12*b1)*dt^3/d4-1/6*dt^2)*b1)*fPt+(-6*alpha*a1*dt^3/d4-6*dt^3/d4*b1)*fPtt
      ! Ptt = c4PttLE*LE + c4PttE*E + c4PttEm*Em + c4PttP*P + c4PttPm*Pm + c4PttfE*fE + c4PttfP*fP
      !    + c4PttLLE*LLE + c4PttLP*LP + c4PttLEm*LEm + c4PttLPm*LPm+ c4PttLfE*LfE + c4PttLfP*LfP
      !    + c4PttfEt*fEt + c4PttfEtt*fEtt + c4PttfPt*fPt+ c4PttfPtt*fPtt
      d4=144.+(12.*a0*alpha+12.*b0)*dt**2+(72.*alpha*a1+72.*b1)*dt

      c4PttLE=-.166666666666666667*(-432.*a1-3.*a0*alpha*b1*c2PttE*dt**4+3.*a1*alpha*b0*c2PttE*dt**4+36.*a0*a1**2*alpha**2*c2EtLE*dt**2-36.*a1**2*alpha**2*b0*c2PtLE*dt**2-108.*a1**2*alpha**2*b1*c2PttLE*dt**2-108.*a1*alpha*b1**2*c2PttLE*dt**2+a1**2*alpha**2*c2PttLE*d4*dt-a0*b1*c2EtLE*d4*dt+b0*b1*c2PtLE*d4*dt-36.*b1**3*c2PttLE*dt**2-36.*a1*b0*dt**2+36.*a1**3*alpha**2*dt**2+36.*a1*b1**2*dt**2-36.*a1**3*alpha**3*c2PttLE*dt**2+36.*a0*b1**2*c2EtLE*dt**2+72.*a1**2*alpha*b1*dt**2-36.*a1**2*alpha*c2EtE*dt**2-36.*b0*b1**2*c2PtLE*dt**2+36.*a1*alpha*c2PttE*dt**2-36.*a1*b1*c2EtE*dt**2-a1**2*alpha*d4*dt+b1**2*c2PttLE*d4*dt-a1*b1*d4*dt+a1*c2EtE*d4*dt+72.*a0*a1*alpha*b1*c2EtLE*dt**2-72.*a1*alpha*b0*b1*c2PtLE*dt**2-a0*a1*alpha*c2EtLE*d4*dt+a1*alpha*b0*c2PtLE*d4*dt+2.*a1*alpha*b1*c2PttLE*d4*dt+36.*a0*b1*dt**2)*dt/d4
      c4PttLLE=.166666666666666667*dt**2*(3.*a0*alpha*b1*c2PttLE*dt**3-3.*a1*alpha*b0*c2PttLE*dt**3-3.*a0*b1*dt**3+36.*a1**2*alpha*c2EtLE*dt+3.*a1*b0*dt**3-36.*a1*alpha*c2PttLE*dt+36.*a1*b1*c2EtLE*dt-a1*c2EtLE*d4+36.*a1*dt)/d4
      c4PttE=((6.00000000000000001*c2PttE*a1**3*alpha**3+(-6.00000000000000001*c2EtE*a0+6.00000000000000001*c2PtE*b0+18.*c2PttE*b1)*alpha**2*a1**2+(18.*c2PttE*b1**2+(-12.*c2EtE*a0+12.*c2PtE*b0)*b1)*alpha*a1+6.00000000000000001*c2PttE*b1**3+(-6.00000000000000001*c2EtE*a0+6.00000000000000001*c2PtE*b0)*b1**2)*dt**4+(-.166666666666666667*c2PttE*a1**2*d4*alpha**2+(-.333333333333333334*c2PttE*b1*d4+(.166666666666666667*c2EtE*a0-.166666666666666667*c2PtE*b0)*d4)*alpha*a1-.166666666666666667*c2PttE*b1**2*d4+(.166666666666666667*c2EtE*a0-.166666666666666667*c2PtE*b0)*d4*b1)*dt**3+((-60.0000000000000001*a0*alpha+12.*b0)*a1-72.0000000000000001*a0*b1)*dt**2+a0*d4*dt+144.*a1)/d4/dt
      c4PttEm=((6.00000000000000001*c2PttEm*a1**3*alpha**3+(-6.00000000000000001*c2EtEm*a0+6.00000000000000001*c2PtEm*b0+18.*c2PttEm*b1)*alpha**2*a1**2+(18.*c2PttEm*b1**2+(-12.*c2EtEm*a0+12.*c2PtEm*b0)*b1)*alpha*a1+6.00000000000000001*c2PttEm*b1**3+(-6.00000000000000001*c2EtEm*a0+6.00000000000000001*c2PtEm*b0)*b1**2)*dt**4+(-.166666666666666667*c2PttEm*a1**2*d4*alpha**2+(-.333333333333333334*c2PttEm*b1*d4+(.166666666666666667*c2EtEm*a0-.166666666666666667*c2PtEm*b0)*d4)*alpha*a1-.166666666666666667*c2PttEm*b1**2*d4+(.166666666666666667*c2EtEm*a0-.166666666666666667*c2PtEm*b0)*d4*b1)*dt**3+(-6.00000000000000001*a0*alpha-6.00000000000000001*b0)*a1*dt**2+(36.0000000000000001*a1**2*alpha+36.0000000000000001*a1*b1)*dt+(-.500000000000000001*d4-72.0000000000000001)*a1)/d4/dt
      c4PttP=((6.00000000000000001*c2PttP*b1**3+(18.*a1*alpha*c2PttP-6.00000000000000001*c2EtP*a0+6.00000000000000001*c2PtP*b0)*b1**2+(18.*c2PttP*alpha**2*a1**2+(-12.*c2EtP*a0+12.*c2PtP*b0)*a1*alpha)*b1+6.00000000000000001*c2PttP*a1**3*alpha**3+(-6.00000000000000001*c2EtP*a0+6.00000000000000001*c2PtP*b0)*a1**2*alpha**2)*dt**4+(-.166666666666666667*c2PttP*b1**2*d4+(-.333333333333333334*c2PttP*a1*d4*alpha+(.166666666666666667*c2EtP*a0-.166666666666666667*c2PtP*b0)*d4)*b1-.166666666666666667*c2PttP*a1**2*d4*alpha**2+(.166666666666666667*c2EtP*a0-.166666666666666667*c2PtP*b0)*d4*a1*alpha)*dt**3+((-12.*a0*alpha+60.0000000000000001*b0)*b1+72.0000000000000001*a1*b0*alpha)*dt**2-1.*b0*d4*dt-144.*b1)/d4/dt
      c4PttPm=((6.00000000000000001*c2PttPm*b1**3+(18.*a1*alpha*c2PttPm-6.00000000000000001*c2EtPm*a0+6.00000000000000001*c2PtPm*b0)*b1**2+(18.*c2PttPm*alpha**2*a1**2+(-12.*c2EtPm*a0+12.*c2PtPm*b0)*a1*alpha)*b1+6.00000000000000001*c2PttPm*a1**3*alpha**3+(-6.00000000000000001*c2EtPm*a0+6.00000000000000001*c2PtPm*b0)*a1**2*alpha**2)*dt**4+(-.166666666666666667*c2PttPm*b1**2*d4+(-.333333333333333334*c2PttPm*a1*d4*alpha+(.166666666666666667*c2EtPm*a0-.166666666666666667*c2PtPm*b0)*d4)*b1-.166666666666666667*c2PttPm*a1**2*d4*alpha**2+(.166666666666666667*c2EtPm*a0-.166666666666666667*c2PtPm*b0)*d4*a1*alpha)*dt**3+(6.00000000000000001*a0*alpha+6.00000000000000001*b0)*b1*dt**2+(-36.0000000000000001*a1*b1*alpha-36.0000000000000001*b1**2)*dt+(.500000000000000001*d4+72.0000000000000001)*b1)/d4/dt
      c4PttLP=.166666666666666667*dt**2*(3.*a0*alpha*b1*c2PttP*dt**3-3.*a1*alpha*b0*c2PttP*dt**3+36.*a1**2*alpha*c2EtP*dt-36.*a1*alpha*c2PttP*dt+36.*a1*b1*c2EtP*dt-a1*c2EtP*d4)/d4
      c4PttLEm=.166666666666666667*dt**2*(3.*a0*alpha*b1*c2PttEm*dt**3-3.*a1*alpha*b0*c2PttEm*dt**3+36.*a1**2*alpha*c2EtEm*dt-36.*a1*alpha*c2PttEm*dt+36.*a1*b1*c2EtEm*dt-a1*c2EtEm*d4)/d4
      c4PttLPm=.166666666666666667*dt**2*(3.*a0*alpha*b1*c2PttPm*dt**3-3.*a1*alpha*b0*c2PttPm*dt**3+36.*a1**2*alpha*c2EtPm*dt-36.*a1*alpha*c2PttPm*dt+36.*a1*b1*c2EtPm*dt-a1*c2EtPm*d4)/d4
      c4PttfE=-.166666666666666667*dt*(-432.*a1-36.*a1*b0*dt**2+36.*a1**3*alpha**2*dt**2+36.*a1*b1**2*dt**2+72.*a1**2*alpha*b1*dt**2-a1**2*alpha*d4*dt-a1*b1*d4*dt-36.*b1**3*c2PttfE*dt**2+72.*a0*a1*alpha*b1*c2EtfE*dt**2-72.*a1*alpha*b0*b1*c2PtfE*dt**2-a0*a1*alpha*c2EtfE*d4*dt+a1*alpha*b0*c2PtfE*d4*dt+2.*a1*alpha*b1*c2PttfE*d4*dt+36.*a0*a1**2*alpha**2*c2EtfE*dt**2-36.*a1**2*alpha**2*b0*c2PtfE*dt**2-108.*a1**2*alpha**2*b1*c2PttfE*dt**2+a1**2*alpha**2*c2PttfE*d4*dt-108.*a1*alpha*b1**2*c2PttfE*dt**2-a0*b1*c2EtfE*d4*dt+b0*b1*c2PtfE*d4*dt-36.*a1**3*alpha**3*c2PttfE*dt**2+36.*a0*b1**2*c2EtfE*dt**2-36.*b0*b1**2*c2PtfE*dt**2+b1**2*c2PttfE*d4*dt+36.*a0*b1*dt**2)/d4
      c4PttfP=((6.00000000000000001*c2PttfP*a1**3*alpha**3+(-6.00000000000000001*a0*c2EtfP+6.00000000000000001*b0*c2PtfP+18.*b1*c2PttfP)*alpha**2*a1**2+(18.*b1**2*c2PttfP+(-12.*a0*c2EtfP+12.*b0*c2PtfP)*b1)*alpha*a1+6.00000000000000001*c2PttfP*b1**3+(-6.00000000000000001*a0*c2EtfP+6.00000000000000001*b0*c2PtfP)*b1**2)*dt**3+(-.166666666666666667*c2PttfP*a1**2*d4*alpha**2+(-.333333333333333334*c2PttfP*b1*d4+(.166666666666666667*a0*c2EtfP-.166666666666666667*b0*c2PtfP)*d4)*alpha*a1-.166666666666666667*c2PttfP*b1**2*d4+(.166666666666666667*a0*c2EtfP-.166666666666666667*b0*c2PtfP)*d4*b1)*dt**2+(-72.0000000000000001*alpha*a1-72.0000000000000001*b1)*dt+d4)/d4
      c4PttfEt=.166666666666666667*a1*dt**2*(36.*alpha*a1*dt+36.*b1*dt-d4)/d4
      c4PttfPt=-.166666666666666667*dt**2*(36.*alpha*a1*dt+36.*b1*dt-d4)*(a1*alpha+b1)/d4
      c4PttLfE=.166666666666666667*dt**2*(3.*a0*alpha*b1*c2PttfE*dt**3-3.*a1*alpha*b0*c2PttfE*dt**3-3.*a0*b1*dt**3+36.*a1**2*alpha*c2EtfE*dt+3.*a1*b0*dt**3-36.*a1*alpha*c2PttfE*dt+36.*a1*b1*c2EtfE*dt-a1*c2EtfE*d4+36.*a1*dt)/d4
      c4PttLfP=.166666666666666667*dt**2*(3.*a0*alpha*b1*c2PttfP*dt**3-3.*a1*alpha*b0*c2PttfP*dt**3+36.*a1**2*alpha*c2EtfP*dt-36.*a1*alpha*c2PttfP*dt+36.*a1*b1*c2EtfP*dt-a1*c2EtfP*d4)/d4
      c4PttfEtt=(-.5*a0*b1*dt**2+.5*a1*b0*dt**2+6.0*a1)*dt**3/d4
      c4PttfPtt=(-6.*alpha*a1-6.*b1)*dt**3/d4

      ! --------------- Here is Ptttt to second-order ------------

      ! Ptttt = ((-alpha*c2PttLE+1)*a0+(-(a0*c2EtLE+(-alpha*c2PttLE+1)*a1-c2PtLE*b0-c2PttLE*b1)*alpha+c2EtE)*a1-c2PttLE*b0-(a0*c2EtLE+(-alpha*c2PttLE+1)*a1-c2PtLE*b0-c2PttLE*b1)*b1)*LE+c2EtLE*a1*LLE+c2EtEm*a1*LEm+c2EtP*a1*LP+c2EtPm*a1*LPm+(-c2PttE*alpha*a0-(-a1*alpha*c2PttE+a0*c2EtE-b0*c2PtE-b1*c2PttE)*alpha*a1-c2PttE*b0-(-a1*alpha*c2PttE+a0*c2EtE-b0*c2PtE-b1*c2PttE)*b1)*E+(-alpha*c2PttEm*a0-(-a1*alpha*c2PttEm+a0*c2EtEm-b0*c2PtEm-b1*c2PttEm)*alpha*a1-c2PttEm*b0-(-a1*alpha*c2PttEm+a0*c2EtEm-b0*c2PtEm-b1*c2PttEm)*b1)*Em+(-alpha*c2PttP*a0-(-a1*alpha*c2PttP+a0*c2EtP-b0*c2PtP-b1*c2PttP)*alpha*a1-c2PttP*b0-(-a1*alpha*c2PttP+a0*c2EtP-b0*c2PtP-b1*c2PttP)*b1)*P+(-alpha*c2PttPm*a0-(-a1*alpha*c2PttPm+a0*c2EtPm-b0*c2PtPm-b1*c2PttPm)*alpha*a1-c2PttPm*b0-(-a1*alpha*c2PttPm+a0*c2EtPm-b0*c2PtPm-b1*c2PttPm)*b1)*Pm+((-alpha*c2PttfE+1)*a0-(c2EtfE*a0+(-alpha*c2PttfE+1)*a1-c2PtfE*b0-c2PttfE*b1)*alpha*a1-c2PttfE*b0-(c2EtfE*a0+(-alpha*c2PttfE+1)*a1-c2PtfE*b0-c2PttfE*b1)*b1)*fE+c2EtfE*a1*LfE+a1*fEt+(-c2PttfP*alpha*a0-(-a1*alpha*c2PttfP+a0*c2EtfP-b0*c2PtfP-b1*c2PttfP)*alpha*a1-c2PttfP*b0-(-a1*alpha*c2PttfP+a0*c2EtfP-b0*c2PtfP-b1*c2PttfP)*b1)*fP+c2EtfP*a1*LfP+(-a1*alpha-b1)*fPt+fPtt
      ! Ptttt = c2PttLE*LE + c2PttE*E + c2PttEm*Em + c2PttP*P + c2PttPm*Pm + c2PttfE*fE + c2PttfP*fP
      !      + c2PttLLE*LLE + c2PttLP*LP + c2PttLEm*LEm + c2PttLPm*LPm+ c2PttLfE*LfE + c2PttLfP*LfP
      !      + c2PttfEt*fEt + c2PttfEtt*fEtt + c2PttfPt*fPt+ c2PttfPtt*fPtt
      c2PttttLE=(alpha**2*c2PttLE-1.*alpha)*a1**2+((-1.*a0*c2EtLE+c2PtLE*b0+2.*c2PttLE*b1)*alpha-1.*b1+c2EtE)*a1-1.*a0*alpha*c2PttLE+(b1**2-1.*b0)*c2PttLE+(-1.*a0*c2EtLE+c2PtLE*b0)*b1+a0
      c2PttttLLE=1.*c2EtLE*a1
      c2PttttE=(alpha**2*a1**2+(2.*a1*b1-1.*a0)*alpha+b1**2-1.*b0)*c2PttE+(-1.*c2EtE*a0+c2PtE*b0)*a1*alpha+(-1.*c2EtE*a0+c2PtE*b0)*b1
      c2PttttEm=(alpha**2*a1**2+(2.*a1*b1-1.*a0)*alpha+b1**2-1.*b0)*c2PttEm+(-1.*c2EtEm*a0+c2PtEm*b0)*a1*alpha+(-1.*c2EtEm*a0+c2PtEm*b0)*b1
      c2PttttP=(alpha**2*a1**2+(2.*a1*b1-1.*a0)*alpha+b1**2-1.*b0)*c2PttP+(-1.*c2EtP*a0+c2PtP*b0)*a1*alpha+(-1.*c2EtP*a0+c2PtP*b0)*b1
      c2PttttPm=(alpha**2*a1**2+(2.*a1*b1-1.*a0)*alpha+b1**2-1.*b0)*c2PttPm+(-1.*c2EtPm*a0+c2PtPm*b0)*a1*alpha+(-1.*c2EtPm*a0+c2PtPm*b0)*b1
      c2PttttLP=1.*c2EtP*a1
      c2PttttLEm=1.*c2EtEm*a1
      c2PttttLPm=1.*c2EtPm*a1
      c2PttttfE=a1**2*alpha**2*c2PttfE+(-1.*a1**2+(-1.*c2EtfE*a0+c2PtfE*b0+2.*c2PttfE*b1)*a1-1.*a0*c2PttfE)*alpha-1.*a1*b1+(b1**2-1.*b0)*c2PttfE+(-1.*c2EtfE*a0+c2PtfE*b0)*b1+a0
      c2PttttfP=(alpha**2*a1**2+(2.*a1*b1-1.*a0)*alpha+b1**2-1.*b0)*c2PttfP+(-1.*a0*c2EtfP+b0*c2PtfP)*a1*alpha+(-1.*a0*c2EtfP+b0*c2PtfP)*b1
      c2PttttfEt=1.*a1
      c2PttttfPt=-1.*alpha*a1-1.*b1
      c2PttttLfE=1.*c2EtfE*a1
      c2PttttLfP=1.*c2EtfP*a1
      c2PttttfEtt=0.
      c2PttttfPtt=1.
                        if( .false. .and. twilightZone.eq.1 )then
                          write(*,'(" LEFT: alpha,dt,d4,d1=",4e12.4)') alpha,dt,d4,d1
                          write(*,'(" a0,a1,b0,b1=",4e12.4)') a0,a1,b0,b1
                          write(*,'(" c2EtE,c2PtE,c2PttfP=",3e12.4)') c2EtE,c2PtE,c2PttfP
                          write(*,'(" c4PttLE,c4PttE,c4PttEm,c4PttP,c4PttPm,c4PttfE=",6e12.4)') c4PttLE,c4PttE,c4PttEm,c4PttP,c4PttPm,c4PttfE
                          write(*,'(" c4PttfP,c4PttLLE,c4PttLP,c4PttLEm,c4PttLPm,c4PttLfE,c4PttLfP=",7e12.4)') c4PttfP,c4PttLLE,c4PttLP,c4PttLEm,c4PttLPm,c4PttLfE,c4PttLfP
                          write(*,'(" c4PttfEt,c4PttfEtt,c4PttfPt,c4PttfPtt=",4e12.4)') c4PttfEt,c4PttfEtt,c4PttfPt,c4PttfPtt
                        end if
                        ! Coeff of E in P.tt (4th order)
                        c2PttEsum1 = c2PttEsum1 + c2PttE
                        ! Coeff of LE1 in P.tt (4th order)
                        c2PttLEsum1  = c2PttLEsum1 + c2PttLE
                        ! Coeff of LE1 in P.tt (4th order)
                        c4PttLEsum1  = c4PttLEsum1 + c4PttLE
                        ! Coeff of LLE1 in P.tt
                        c4PttLLEsum1 = c4PttLLEsum1 + c4PttLLE
                        ! Coeff of LE1 and LLE1 in P.tttt
                        c2PttttLEsum1 =c2PttttLEsum1 +c2PttttLE
                        c2PttttLLEsum1=c2PttttLLEsum1+c2PttttLLE
                        do n=0,nd-1
                          pc = n + jv*nd 
                          ec = ex +n
                          pv   =  p1(i1,i2,i3,pc)
                          pvn  =  p1n(i1,i2,i3,pc)
                          ev    =  u1(i1,i2,i3,ec)
                          evn   =  u1n(i1,i2,i3,ec)
                          ! Left: u1x,u1y, u1xx, u1yy, u1Lap (ex)
                          !       v1x,v1y, v1xx, v1yy, v1Lap (ey) 
                          ! These next derivatives may only be needed to order2, (use order 4 for testing TZ polynomials)
                          ! perl $ORDER=2;
                            ! uu1 in the next statement defines names of intermediate values
                              uu1=p1(i1,i2,i3,pc) ! in the rectangular case just eval the solution
                               p1x = (p1(i1-2,i2,i3,pc)-8.*p1(i1-1,i2,i3,pc)+8.*p1(i1+1,i2,i3,pc)-p1(i1+2,i2,i3,pc))/(12.*dx1(0))
                               p1y = (p1(i1,i2-2,i3,pc)-8.*p1(i1,i2-1,i3,pc)+8.*p1(i1,i2+1,i3,pc)-p1(i1,i2+2,i3,pc))/(12.*dx1(1))
                               p1xx = (-p1(i1-2,i2,i3,pc)+16.*p1(i1-1,i2,i3,pc)-30.*p1(i1,i2,i3,pc)+16.*p1(i1+1,i2,i3,pc)-p1(i1+2,i2,i3,pc))/(12.*dx1(0)**2)
                               p1yy = (-p1(i1,i2-2,i3,pc)+16.*p1(i1,i2-1,i3,pc)-30.*p1(i1,i2,i3,pc)+16.*p1(i1,i2+1,i3,pc)-p1(i1,i2+2,i3,pc))/(12.*dx1(1)**2)
                             p1Lap = p1xx+ p1yy
                            ! write(*,'("LEFT: p1x,p1y,p1xx,p1yy,p1Lap=",5(1pe12.4))') p1x,p1y,p1xx,p1yy,p1Lap
                            LP  = (c1**2)*p1Lap
                              uu1=p1n(i1,i2,i3,pc) ! in the rectangular case just eval the solution
                               p1nx = (p1n(i1-2,i2,i3,pc)-8.*p1n(i1-1,i2,i3,pc)+8.*p1n(i1+1,i2,i3,pc)-p1n(i1+2,i2,i3,pc))/(12.*dx1(0))
                               p1ny = (p1n(i1,i2-2,i3,pc)-8.*p1n(i1,i2-1,i3,pc)+8.*p1n(i1,i2+1,i3,pc)-p1n(i1,i2+2,i3,pc))/(12.*dx1(1))
                               p1nxx = (-p1n(i1-2,i2,i3,pc)+16.*p1n(i1-1,i2,i3,pc)-30.*p1n(i1,i2,i3,pc)+16.*p1n(i1+1,i2,i3,pc)-p1n(i1+2,i2,i3,pc))/(12.*dx1(0)**2)
                               p1nyy = (-p1n(i1,i2-2,i3,pc)+16.*p1n(i1,i2-1,i3,pc)-30.*p1n(i1,i2,i3,pc)+16.*p1n(i1,i2+1,i3,pc)-p1n(i1,i2+2,i3,pc))/(12.*dx1(1)**2)
                             p1nLap = p1nxx+ p1nyy
                            ! write(*,'("LEFT: p1nxx,p1nyy,p1nLap=",3e12.4)') p1nxx,p1nyy,p1nLap
                            LPm = (c1**2)*p1nLap
                            pvx  = p1x
                            pvy  = p1y
                            pvnx = p1nx
                            pvny = p1ny
                          ! Accumulate: SUM_m Pm,tt
                          ptta = c4PttLE*LE1(n) + c4PttE*ev + c4PttEm*evn + c4PttP*pv + c4PttPm*pvn + c4PttfE*fev1(n) + c4PttfP*fpv1(n,jv) + c4PttLLE*LLE1(n) + c4PttLP*LP + c4PttLEm*LE1m(n) + c4PttLPm*LPm+ c4PttLfE*LfE1(n) + c4PttLfP*LfP1(n,jv) + c4PttfEt*fEt1(n) + c4PttfEtt*fEtt1(n) + c4PttfPt*fPt1(n,jv)+ c4PttfPtt*fPtt1(n,jv)
                          ! ---- Compute fp1 = P.tt 
                          fp1(n) = fp1(n) + ptta
                          ! ----- Compute fPttx1 = (P.tt).x , fPtty1 = (P.tt).y  (second order)
                          pttxa = c2PttLE*LEx1(n) + c2PttE*evx1(n) + c2PttEm*evnx1(n) + c2PttP*pvx + c2PttPm*pvnx + c2PttfE*fevx1(n) + c2PttfP*fpvx1(n,jv)
                          pttya = c2PttLE*LEy1(n) + c2PttE*evy1(n) + c2PttEm*evny1(n) + c2PttP*pvy + c2PttPm*pvny + c2PttfE*fevy1(n) + c2PttfP*fpvy1(n,jv)
                          fPttx1(n) = fPttx1(n) + pttxa
                          fPtty1(n) = fPtty1(n) + pttya
                          ! ----- Compute fLPtt1 = L(P.tt) (second order)
                          Lptta = c2PttLE*LLE1(n) + c2PttE*LE1(n) + c2PttEm*LE1m(n) + c2PttP*LP + c2PttPm*LPm + c2PttfE*LfE1(n) + c2PttfP*LfP1(n,jv)
                          fLPtt1(n) = fLPtt1(n) + Lptta
                          ! ----- Compute fPtttt1 = P.tttt
                          ptttta= c2PttttLE*LE1(n) + c2PttttE*ev + c2PttttEm*evn + c2PttttP*pv + c2PttttPm*pvn + c2PttttfE*fev1(n) + c2PttttfP*fpv1(n,jv) + c2PttttLLE*LLE1(n) + c2PttttLP*LP + c2PttttLEm*LE1m(n) + c2PttttLPm*LPm+ c2PttttLfE*LfE1(n) + c2PttttLfP*LfP1(n,jv) + c2PttttfEt*fEt1(n) + c2PttttfEtt*fEtt1(n) + c2PttttfPt*fPt1(n,jv)+ c2PttttfPtt*fPtt1(n,jv)
                          fPtttt1(n) = fPtttt1(n) + ptttta
                          if( .false. .and. twilightZone.eq.1 )then
                            write(*,'("")')
                            write(*,'("DI4:LEFT: i1,i2=",2i3," jv=",i2," n=",i2," ptta,ptte=",2e12.4)') i1,i2,jv,n,ptta,pevtt1(n,jv)
                            write(*,'("        : pttxa,pttxe=",2(1pe12.4)," pttya,pttye=",2(1pe12.4))') pttxa,pevttx1(n,jv),pttya,pevtty1(n,jv)
                            write(*,'("        : ptttta,ptttte=",2(1pe12.4),/)') ptttta,pevtttt1(n,jv)
                            write(*,'(" c2PttLE,c2PttE,c2PttEm,c2PttP,c2PttPm,c2PttfE,c2PttfP=",7(1pe12.4))') c2PtLE,c2PtE,c2PtEm,c2PtP,c2PtPm,c2PtfE,c2PtfP
                            write(*,'(" LEx1,evx1,evnx1,pvx,pvnx,fevx1,fpvx1=",7(1pe12.4))') LEx1(n),evx1(n),evnx1(n),pvx,pvnx,fevx1(n),fpvx1(n,jv)
                            write(*,'(" LEy1,evy1,evny1,pvy,pvny,fevy1,fpvy1=",7(1pe12.4))') LEy1(n),evy1(n),evny1(n),pvy,pvny,fevy1(n),fpvy1(n,jv)
                            ! write(*,'(" LE1,LLE1,LE1m,LP,LPm=",5e12.4)') LE1(n),LLE1(n),LE1m(n),LP,LPm
                            ! write(*,'(" LfE1,LfP1,fEt1,fEtt1,fPt1,fPtt1=",6e12.4)') LfE1(n),LfP1(n,jv),fEt1(n),fEtt1(n),fPt1(n,jv),fPtt1(n,jv)
                            ! write(*,'(" ev,evn,pv,pvn,fev1,fpv1=",6e12.4)')ev,evn,pv,pvn,fev1(n),fpv1(n,jv)
                          end if
                        end do ! end do n 
                      end do ! enddo jv
                    !-----------------------------
                    ! no dispersion
                    !-----------------------------
                    elseif( dispersionModel1.eq.noDispersion .and. nonlinearModel1.eq.noNonlinearModel) then
                      ! do nothing, dispersive forcing is 0
                    end if
                  ! eval nonlinear dispersive forcings for domain 2
                  ! getDispersiveForcingOrder4(RIGHT,j1,j2,j3, fp2, fpv2,fev2, p2,p2n,p2m, u2,u2n,u2m, dispersionModel2,numberOfPolarizationVectors2,alphaP2,!             c2PttEsum2,c2PttLEsum2,c4PttLEsum2,c4PttLLEsum2,c2PttttLEsum2,c2PttttLLEsum2,a0v2,a1v2,!             b0v2,b1v2,LE2,LLE2,LE2m,LfE2,LfP2,fEt2,fEtt2,fPt2,fPtt2,pevtt2,pevttx2,pevtty2,pevtttt2,evx2,evy2,evnx2,evny2,!             fevx2,fevy2,fpvx2,fpvy2,LEx2,LEy2,fPttx2,fPtty2,fLPtt2,fPtttt2)
                    ! pre-assign 0 values
                    do n=0,nd-1 ! dispersive forcing in jump conditions
                      fp2(n)=0.
                      fPttx2(n) =0.
                      fPtty2(n) =0.
                      fLPtt2(n) =0.
                      fPtttt2(n)=0.
                    end do
                    !-----------------------------
                    ! MLA
                    !-----------------------------
                    ! only do this for MLA (dispersive and nonlinear multi-level)
                    if( dispersionModel2.ne.noDispersion .and. nonlinearModel2.ne.noNonlinearModel) then
                      nce = pxc+nd*numberOfPolarizationVectors2
                      ! -----------------------------------------
                      ! order 2 (E, P, N) at the interface (fictitious step)
                      !------------------------------------------
                      ! dimension loop for E and P
                      do n=0,nd-1
                        ec = ex +n 
                        ev0  =  u2n(j1,j2,j3,ec)
                        ev   =  u2(j1,j2,j3,ec) ! time where we need to fill in ghost points
                        pSum=0.
                        do jv=0,numberOfPolarizationVectors2-1
                          pc = n + jv*nd  
                          pv0 =  p2n(j1,j2,j3,pc)
                          pv  =  p2(j1,j2,j3,pc)
                          pvn = 2.*pv-p2n(j1,j2,j3,pc) + 0.5*dt*b1v2(jv)*p2n(j1,j2,j3,pc) - dtsq*b0v2(jv)*pv + dtsq*fpv2(n,jv)
                          do na = 0,numberOfAtomicLevels2-1 ! \Delta N^n*E^n
                            pvn = pvn + dtsq*pnec2(jv,na)*q2(j1,j2,j3,na)*ev
                          enddo ! na
                          pvec(n,jv)= pvn/( 1.+.5*dt*b1v2(jv) ) ! time + dt
                          ! #If "p2" eq "p1"
                          ! call ogderiv(ep, 0,0,0,0, xy1(j1,j2,j3,0),xy1(j1,j2,j3,1),0.,t+dt,pxc+jv*nd+n, pe(n)   )
                          ! ! call ogderiv(ep, 1,0,0,0, xy1(j1,j2,j3,0),xy1(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pet(n)   )
                          ! ! call ogderiv(ep, 2,0,0,0, xy1(j1,j2,j3,0),xy1(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pett(n)   )
                          ! #Else
                          ! call ogderiv(ep, 0,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t+dt,pxc+jv*nd+n, pe(n)   )
                          ! ! call ogderiv(ep, 1,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pet(n)   )
                          ! ! call ogderiv(ep, 2,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pett(n)   )
                          ! #End
                          ! print *, '---------Dispersive forcing 2----------'
                          ! print *, pvec(n,jv),pe(n),pvec(n,jv)-pe(n)
                          pSum = pSum + pvec(n,jv) -2.*pv + pv0 ! keep sum
                          ! print *, '++++++Dispersive forcing+++++++++++++++'
                          ! print *, 'check 2nd output P: ', pc,jv,pvec(n,jv)
                        enddo ! jv
                        ! second order update of E
                        evec(n) = (2.*ev-ev0) + cSq2*dtsq*LE2(n)/cSq2 - alphaP2*pSum + dtsq*fev2(n) ! cSq2 is already in LE2
                        ! print *, '++++++Dispersive forcing+++++++++++++++'
                        ! print *, ev,ev0, cSq2, dtsq,LE2(n),alphaP2,pSum,fev2(n)
                        ! print *, 'check 2nd output E: ', ec,evec(n)
                      enddo ! n
                      ! N outside of space loop
                      ! 1st derivative
                      do na=0,numberOfAtomicLevels2-1
                        qt(na) = fnv2(na)
                        do jv = 0,numberOfAtomicLevels2-1
                          qt(na) = qt(na)+prc2(na,jv)*q2(j1,j2,j3,jv)
                        enddo
                        do n=0,nd-1
                          do jv=0,numberOfPolarizationVectors2-1
                            qt(na) = qt(na) + peptc2(na,jv)*u2(j1,j2,j3,ex+n)*(pvec(n,jv)-p2n(j1,j2,j3,n+jv*nd))/(2.*dt)
                          enddo
                        enddo
                      enddo
                      ! 2nd derivative
                      do na=0,numberOfAtomicLevels2-1
                        qtt(na) = fntv2(na)
                        do jv = 0,numberOfAtomicLevels2-1
                          qtt(na) = qtt(na)+prc2(na,jv)*qt(jv)
                        enddo
                        do n=0,nd-1
                          do jv=0,numberOfPolarizationVectors2-1
                            qtt(na) = qtt(na) + peptc2(na,jv)*(evec(n)-u2n(j1,j2,j3,ex+n))/(2.*dt)*(pvec(n,jv)-p2n(j1,j2,j3,n+jv*nd))/(2.*dt)+ peptc2(na,jv)*u2(j1,j2,j3,ex+n)*(pvec(n,jv)-2.*p2(j1,j2,j3,n+jv*nd)+p2n(j1,j2,j3,n+jv*nd))/(dtsq)
                          enddo
                        enddo
                      enddo
                      ! taylor expansion
                      do na=0,numberOfAtomicLevels2-1
                        qv(na) = q2(j1,j2,j3,na) + dt*qt(na) + dtsq/2.*qtt(na)
                      enddo
                      !----------------------------------------
                      ! order 4 update of P at interface (fictitious step)
                      !----------------------------------------
                      ! second order accurate terms
                      do n=0,nd-1
                        do jv = 0,numberOfPolarizationVectors2-1
                          ! #If "p2" eq "p1"
                          ! call ogderiv(ep, 0,0,0,0, xy1(j1,j2,j3,0),xy1(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pe(n)   )
                          ! call ogderiv(ep, 1,0,0,0, xy1(j1,j2,j3,0),xy1(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pet(n)   )
                          ! call ogderiv(ep, 2,0,0,0, xy1(j1,j2,j3,0),xy1(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pett(n)   )
                          ! #Else
                          ! call ogderiv(ep, 0,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pe(n)   )
                          ! call ogderiv(ep, 1,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pet(n)   )
                          ! call ogderiv(ep, 2,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pett(n)   )
                          ! #End
                          ! pvec(n,jv) = pe(n)
                          ! ptv(n,jv) = pet(n)
                          ! pttv(n,jv) = pett(n)
                          ptv(n,jv) = (pvec(n,jv)-p2n(j1,j2,j3,n+jv*nd))/(2.*dt)
                          pttv(n,jv) = (pvec(n,jv)-2.*p2(j1,j2,j3,n+jv*nd)+p2n(j1,j2,j3,n+jv*nd))/dtsq
                          ptttv(n,jv) = -b1v2(jv)*pttv(n,jv)-b0v2(jv)*ptv(n,jv)+fPt2(n,jv)
                          do na = 0,numberOfAtomicLevels2-1 ! update using ODE
                            ptttv(n,jv) = ptttv(n,jv) + pnec2(jv,na)*qt(na)*u2(j1,j2,j3,ex+n) + pnec2(jv,na)*q2(j1,j2,j3,na)*(evec(n)-u2n(j1,j2,j3,ex+n))/(2.*dt)
                          enddo
                          pttttv(n,jv) = -b1v2(jv)*ptttv(n,jv)-b0v2(jv)*pttv(n,jv)+fPtt2(n,jv) ! update using ODE
                          do na = 0,numberOfAtomicLevels2-1
                            pttttv(n,jv) = pttttv(n,jv) + pnec2(jv,na)*qtt(na)*u2(j1,j2,j3,ex+n) + 2.*pnec2(jv,na)*qt(na)*(evec(n)-u2n(j1,j2,j3,ex+n))/(2.*dt) + pnec2(jv,na)*q2(j1,j2,j3,na)*(evec(n)-2.*u2(j1,j2,j3,ex+n)+u2n(j1,j2,j3,ex+n))/dtsq
                            ! print *, '++++++Dispersive forcing+++++++++++++++'
                            ! print *, 'check P derivatives: ', n,jv,ptv(n,jv),pttv(n,jv),ptttv(n,jv),pttttv(n,jv)
                            ! print *, na,pnec2(jv,na),q2(j1,j2,j3,na),qt(na),qtt(na)
                          enddo
                          ! print *, '++++++Dispersive forcing+++++++++++++++'
                          ! print *, 'check P time derivatives: ', n,jv,ptv(n,jv),pttv(n,jv),ptttv(n,jv),pttttv(n,jv)
                          ! print *, fPt2(n,jv),fPtt2(n,jv)
                        enddo
                        ! print *, evec(n),u2(j1,j2,j3,ex+n),u2n(j1,j2,j3,ex+n)
                      enddo
                      ! dimension loop for E and P
                      do n=0,nd-1
                        ec = ex +n 
                        ev   =  u2(j1,j2,j3,ec) ! time where we need to fill in ghost points
                        ! These next derivatives may only be needed to order2, (use order 4 for testing TZ polynomials)
                        ! #perl $ORDER=2; ! should use order 2 since E is only filled at first ghost lines at t
                          ! uu2 in the next statement defines names of intermediate values
                            uu2=u2(j1,j2,j3,ec) ! in the rectangular case just eval the solution
                             e2x = (u2(j1-2,j2,j3,ec)-8.*u2(j1-1,j2,j3,ec)+8.*u2(j1+1,j2,j3,ec)-u2(j1+2,j2,j3,ec))/(12.*dx2(0))
                             e2y = (u2(j1,j2-2,j3,ec)-8.*u2(j1,j2-1,j3,ec)+8.*u2(j1,j2+1,j3,ec)-u2(j1,j2+2,j3,ec))/(12.*dx2(1))
                             e2xx = (-u2(j1-2,j2,j3,ec)+16.*u2(j1-1,j2,j3,ec)-30.*u2(j1,j2,j3,ec)+16.*u2(j1+1,j2,j3,ec)-u2(j1+2,j2,j3,ec))/(12.*dx2(0)**2)
                             e2yy = (-u2(j1,j2-2,j3,ec)+16.*u2(j1,j2-1,j3,ec)-30.*u2(j1,j2,j3,ec)+16.*u2(j1,j2+1,j3,ec)-u2(j1,j2+2,j3,ec))/(12.*dx2(1)**2)
                           e2Lap = e2xx+ e2yy
                          ! write(*,'("RIGHT: u2x,u2y,u2xx,u2yy,u2Lap=",5(1pe12.4))') u2x,u2y,u2xx,u2yy,u2Lap
                          evx0  = e2x
                          evy0  = e2y
                          evLap = e2xx+e2yy
                            ! print *, '++++++Dispersive forcing+++++++++++++++'
                            ! print *, 'check E spatial derivatives: ', n,evx0,evy0,evLap
                        do jv=0,numberOfPolarizationVectors2-1
                          pc = n + jv*nd 
                          ! These next derivatives may only be needed to order2, (use order 4 for testing TZ polynomials)
                          ! #perl $ORDER=2;
                            ! uu1 in the next statement defines names of intermediate values
                              pp2=p2(j1,j2,j3,pc) ! in the rectangular case just eval the solution
                               p2x = (p2(j1-2,j2,j3,pc)-8.*p2(j1-1,j2,j3,pc)+8.*p2(j1+1,j2,j3,pc)-p2(j1+2,j2,j3,pc))/(12.*dx1(0))
                               p2y = (p2(j1,j2-2,j3,pc)-8.*p2(j1,j2-1,j3,pc)+8.*p2(j1,j2+1,j3,pc)-p2(j1,j2+2,j3,pc))/(12.*dx1(1))
                               p2xx = (-p2(j1-2,j2,j3,pc)+16.*p2(j1-1,j2,j3,pc)-30.*p2(j1,j2,j3,pc)+16.*p2(j1+1,j2,j3,pc)-p2(j1+2,j2,j3,pc))/(12.*dx1(0)**2)
                               p2yy = (-p2(j1,j2-2,j3,pc)+16.*p2(j1,j2-1,j3,pc)-30.*p2(j1,j2,j3,pc)+16.*p2(j1,j2+1,j3,pc)-p2(j1,j2+2,j3,pc))/(12.*dx1(1)**2)
                             p2Lap = p2xx+ p2yy
                            LP  = p2Lap
                              pp2=p2n(j1,j2,j3,pc) ! in the rectangular case just eval the solution
                               p2nx = (p2n(j1-2,j2,j3,pc)-8.*p2n(j1-1,j2,j3,pc)+8.*p2n(j1+1,j2,j3,pc)-p2n(j1+2,j2,j3,pc))/(12.*dx1(0))
                               p2ny = (p2n(j1,j2-2,j3,pc)-8.*p2n(j1,j2-1,j3,pc)+8.*p2n(j1,j2+1,j3,pc)-p2n(j1,j2+2,j3,pc))/(12.*dx1(1))
                               p2nxx = (-p2n(j1-2,j2,j3,pc)+16.*p2n(j1-1,j2,j3,pc)-30.*p2n(j1,j2,j3,pc)+16.*p2n(j1+1,j2,j3,pc)-p2n(j1+2,j2,j3,pc))/(12.*dx1(0)**2)
                               p2nyy = (-p2n(j1,j2-2,j3,pc)+16.*p2n(j1,j2-1,j3,pc)-30.*p2n(j1,j2,j3,pc)+16.*p2n(j1,j2+1,j3,pc)-p2n(j1,j2+2,j3,pc))/(12.*dx1(1)**2)
                             p2nLap = p2nxx+ p2nyy
                            LPm = p2nLap
                            pvx  = p2x
                            pvy  = p2y
                            pvnx = p2nx
                            pvny = p2ny
                          ! print *, '++++++Dispersive forcing+++++++++++++++'
                          ! print *, 'check P spatial derivatives: ', n,jv,LP,LPm,pvx,pvy,pvnx,pvny
                          pv0 =  p2n(j1,j2,j3,pc)
                          pv  =  p2(j1,j2,j3,pc)
                          pvn = 2.*pv-p2n(j1,j2,j3,pc) + 0.5*dt*b1v2(jv)*p2n(j1,j2,j3,pc) + dt**4/12.*pttttv(n,jv) + dt**4/6.*b1v2(jv)*ptttv(n,jv) - dtsq*b0v2(jv)*pv + dtsq*fpv2(n,jv)
                          ! pvn = 2.*pv-p2n(j1,j2,j3,pc) + 0.5*dt*b1v2(jv)*p2n(j1,j2,j3,pc) - dtsq*b0v2(jv)*pv + dtsq*fpv2(n,jv)
                          do na = 0,numberOfAtomicLevels2-1 ! \Delta N^n*E^n
                            pvn = pvn + dtsq*pnec2(jv,na)*q2(j1,j2,j3,na)*ev
                          enddo ! na
                          pvec(n,jv)= pvn/( 1.+.5*dt*b1v2(jv) ) ! time + dt
                          ! 4th order accurate term
                          fp2(n) = fp2(n) + (pvec(n,jv)-2.*p2(j1,j2,j3,pc)+p2n(j1,j2,j3,pc))/dtsq  - dt**2/12.*pttttv(n,jv)
                          ! #If "p2" eq "p1"
                          !   call ogderiv(ep, 0,0,0,0, xy1(j1,j2,j3,0),xy1(j1,j2,j3,1),0.,t+dt,pxc+jv*nd+n, pe(n)   )
                          !   call ogderiv(ep, 1,0,0,0, xy1(j1,j2,j3,0),xy1(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pet(n)   )
                          !   call ogderiv(ep, 2,0,0,0, xy1(j1,j2,j3,0),xy1(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pett(n)   )
                          !   call ogderiv(ep, 4,0,0,0, xy1(j1,j2,j3,0),xy1(j1,j2,j3,1),0.,t,pxc+jv*nd+n, petttt(n)   )
                          ! #Else
                          !   call ogderiv(ep, 0,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t+dt,pxc+jv*nd+n, pe(n)   )
                          !   call ogderiv(ep, 1,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pet(n)   )
                          !   call ogderiv(ep, 2,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pett(n)   )
                          !   call ogderiv(ep, 4,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pxc+jv*nd+n, petttt(n)   )
                          ! #End
                          ! print *, '---------Dispersive forcing 4----------'
                          ! print *, pvec(n,jv),pe(n),pvec(n,jv)-pe(n)
                          ! print *, (pvec(n,jv)-2.*p2(j1,j2,j3,pc)+p2n(j1,j2,j3,pc))/dtsq  - dt**2/12.*pttttv(n,jv), pett(n),(pvec(n,jv)-2.*p2(j1,j2,j3,pc)+p2n(j1,j2,j3,pc))/dtsq  - dt**2/12.*pttttv(n,jv)-pett(n)
                          ! fp2(n) = fp2(n) + pett(n)
                          ! 2nd order accurate terms
                          fPtttt2(n) = fPtttt2(n) + pttttv(n,jv)
                          ! fPtttt2(n) = fPtttt2(n) + petttt(n)
                          ! print *, '++++++Dispersive forcing+++++++++++++++'
                          ! print *, 'check P time derivatives 2 and 4: ', n,jv,(pvec(n,jv)-2.*p2(j1,j2,j3,pc)+p2n(j1,j2,j3,pc))/dtsq  - dt**4/12.*pttttv(n,jv),pttttv(n,jv)
                          !--------------------------------
                          ! spatial derivatives
                          !--------------------------------
                          ! nce = pxc+nd*numberOfPolarizationVectors2
                          ! N*E
                          ! #perl $ORDER=2;
                          do na=0,numberOfAtomicLevels2-1
                            ! uu2 in the next statement defines names of intermediate values
                              qq2=q2(j1,j2,j3,na) ! in the rectangular case just eval the solution
                               q2x = (q2(j1-2,j2,j3,na)-8.*q2(j1-1,j2,j3,na)+8.*q2(j1+1,j2,j3,na)-q2(j1+2,j2,j3,na))/(12.*dx2(0))
                               q2y = (q2(j1,j2-2,j3,na)-8.*q2(j1,j2-1,j3,na)+8.*q2(j1,j2+1,j3,na)-q2(j1,j2+2,j3,na))/(12.*dx2(1))
                               q2xx = (-q2(j1-2,j2,j3,na)+16.*q2(j1-1,j2,j3,na)-30.*q2(j1,j2,j3,na)+16.*q2(j1+1,j2,j3,na)-q2(j1+2,j2,j3,na))/(12.*dx2(0)**2)
                               q2yy = (-q2(j1,j2-2,j3,na)+16.*q2(j1,j2-1,j3,na)-30.*q2(j1,j2,j3,na)+16.*q2(j1,j2+1,j3,na)-q2(j1,j2+2,j3,na))/(12.*dx2(1)**2)
                             q2Lap = q2xx+ q2yy
                            qvx  = q2x
                            qvy  = q2y
                            qvLap  = q2Lap
                            qex(na) = evx0*q2(j1,j2,j3,na)+qvx*ev
                            qey(na) = evy0*q2(j1,j2,j3,na)+qvy*ev
                            qeLap(na) = ev*qvLap+q2(j1,j2,j3,na)*evLap+2.*evx0*qvx+2.*evy0*qvy
                            ! print *, '++++++Dispersive forcing+++++++++++++++'
                            ! print *, 'check N spatial derivatives: ', na,qex(na),qey(na),qeLap(na)
                          enddo
                          ! #If "p2" eq "p1"
                          !   call ogderiv(ep, 2,2,0,0, xy1(j1,j2,j3,0),xy1(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pevttxx(n,jv)   )
                          !   call ogderiv(ep, 2,0,2,0, xy1(j1,j2,j3,0),xy1(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pevttyy(n,jv)   )
                          !   call ogderiv(ep, 2,1,0,0, xy1(j1,j2,j3,0),xy1(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pevttx2(n,jv)   )
                          !   call ogderiv(ep, 2,0,1,0, xy1(j1,j2,j3,0),xy1(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pevtty2(n,jv)   )
                          ! #Else
                          !   call ogderiv(ep, 2,2,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pevttxx(n,jv)   )
                          !   call ogderiv(ep, 2,0,2,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pevttyy(n,jv)   )
                          !   call ogderiv(ep, 2,1,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pevttx2(n,jv)   )
                          !   call ogderiv(ep, 2,0,1,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pevtty2(n,jv)   )
                          ! #End
                          ! laplacian
                          LPn = 2.*LP-LPm + 0.5*dt*b1v2(jv)*LPm - dtsq*b0v2(jv)*LP + dtsq*LfP2(n,jv)
                          do na = 0,numberOfAtomicLevels2-1 ! \Delta N^n*E^n
                            LPn = LPn + dtsq*pnec2(jv,na)*qeLap(na)
                          enddo
                          ! time derivatives
                          fLPtt2(n) = fLPtt2(n) + cSq2*(LPn/(1.+.5*dt*b1v2(jv)) - 2.*LP + LPm)/dtsq ! added c^2 to be consistence with GDM version
                          ! fLPtt2(n) = fLPtt2(n) + pevttL2(n,jv)
                          ! print *,'error of laplacian with # of levels',numberOfAtomicLevels2,numberOfPolarizationVectors2,pevttxx(n,jv)+pevttyy(n,jv)-pevttL2(n,jv)
                          ! print *, '++++++Dispersive forcing+++++++++++++++'
                          ! print *, 'check P lap: ', n,jv,LP,LPm,b1v2(jv),dtsq,b0v2(jv),LfP2(n,jv),fLPtt2(n)
                          ! x
                          pttxa = 2.*pvx-pvnx + 0.5*dt*b1v2(jv)*pvnx - dtsq*b0v2(jv)*pvx + dtsq*fpvx2(n,jv)
                          do na = 0,numberOfAtomicLevels2-1 ! \Delta N^n*E^n
                            pttxa = pttxa + dtsq*pnec2(jv,na)*qex(na)
                          enddo
                          ! time derivatives
                          fPttx2(n) = fPttx2(n) + (pttxa/(1.+.5*dt*b1v2(jv)) - 2.*pvx + pvnx)/dtsq
                          ! fPttx2(n) = fPttx2(n) + pevttx2(n,jv)
                          ! print *, '++++++Dispersive forcing+++++++++++++++'
                          ! print *, 'check Pttx: ', n,jv,pvx,pvnx,fpvx2(n,jv),fPttx2(n)
                          ! y
                          pttya = 2.*pvy-pvny + 0.5*dt*b1v2(jv)*pvny - dtsq*b0v2(jv)*pvy + dtsq*fpvy2(n,jv)
                          do na = 0,numberOfAtomicLevels2-1 ! \Delta N^n*E^n
                            pttya = pttya + dtsq*pnec2(jv,na)*qey(na)
                          enddo
                          ! time derivatives
                          fPtty2(n) = fPtty2(n) + (pttya/(1.+.5*dt*b1v2(jv)) - 2.*pvy + pvny)/dtsq
                          ! fPtty2(n) = fPtty2(n) + pevtty2(n,jv)
                          ! print *, '++++++Dispersive forcing+++++++++++++++'
                          ! print *, 'check Ptty: ', n,jv,n,jv,pvy,pvny,fpvy2(n,jv),fPtty2(n)
                          ! if( twilightZone.eq.1 )then
                          !   write(*,'("")')
                          !   write(*,'("DI4:RIGHT: j1,j2=",2i3," jv=",i2," n=",i2)') j1,j2,jv,n
                          !   print *, 'ptt diff',(pvec(n,jv)-2.*p2(j1,j2,j3,pc)+p2n(j1,j2,j3,pc))/dtsq - dt**2/12.*pttttv(n,jv)-pevtt2(n,jv)
                          !   print *, 'pttx diff', (pttxa/(1.+.5*dt*b1v2(jv)) - 2.*pvx + pvnx)/dtsq-pevttx2(n,jv)
                          !   print *, 'ptty diff', (pttya/(1.+.5*dt*b1v2(jv)) - 2.*pvy + pvny)/dtsq-pevtty2(n,jv)
                          !   print *, 'ptttt diff', pttttv(n,jv)-pevtttt2(n,jv)
                          !   print *, 'pttL diff',(LPn/(1.+.5*dt*b1v2(jv)) - 2.*LP + LPm)/dtsq-pevttL2(n,jv)
                          ! end if
                        enddo ! jv
                      enddo ! n
                    !-----------------------------
                    ! GDM
                    !-----------------------------
                    elseif( dispersionModel2.ne.noDispersion .and. nonlinearModel2.eq.noNonlinearModel) then
                      c2PttEsum2=0.
                      c2PttLEsum2=0.
                      c4PttLEsum2=0.
                      c4PttLLEsum2=0.
                      c2PttttLEsum2=0.
                      c2PttttLLEsum2=0.
                      do jv=0,numberOfPolarizationVectors2-1
                        a0=a0v2(jv)
                        a1=a1v2(jv)
                        b0=b0v2(jv)
                        b1=b1v2(jv)
                        alpha=alphaP2
                        ! Second-order coefficients: 
                        ! Ptt = c2PttLE*LE2 + c2PttE*E + c2PttEm*Em + c2PttP*P + c2PttPm*Pm + c2PttfE*fE + c2PttfP*fP
      ! File created by Dropbox/GDM/maple/interface.maple
      d1=2.+(a1*alpha+b1)*dt

      ! --------------- Here is Et to second-order ------------

      ! Et = -1/2/dt*(2*a0*alpha*dt^2-2*b1*dt-4)/d1*E-1/2/dt*(2*b1*dt+4)/d1*Em-1/2/dt*(-b1*dt^3-2*dt^2)/d1*LE-1/2/dt*(-2*alpha*b0*dt^2-2*alpha*b1*dt)/d1*P-alpha*b1/d1*Pm-1/2/dt*(-b1*dt^3-2*dt^2)/d1*fE-dt*alpha*fP/d1
      ! Et = c2EtLE*LE + c2EtE*E + c2EtEm*Em + c2EtP*P + c2EtPm*Pm + c2EtfE*fE + c2EtfP*fP
      !   LE = c^2Delta(E) 
      c2EtLE=(.5*b1*dt**2+1.0*dt)/d1
      c2EtE=(-a0*alpha*dt**2+b1*dt+2.)/dt/d1
      c2EtEm=(-b1*dt-2.)/dt/d1
      c2EtP=alpha*(b0*dt+b1)/d1
      c2EtPm=-1.*alpha*b1/d1
      c2EtfE=(.5*b1*dt**2+1.0*dt)/d1
      c2EtfP=-1.*dt*alpha/d1
      ! --------------- Here is Pt to second-order ------------

      ! Pt = 1/2*(2*a0*dt^2+2*a1*dt)/dt/d1*E-a1/d1*Em+1/2*a1*dt^2/d1*LE+1/2*(2*a1*alpha*dt-2*b0*dt^2+4)/dt/d1*P+1/2*(-2*a1*alpha*dt-4)/dt/d1*Pm+1/2*a1*dt^2/d1*fE+fP*dt/d1
      ! Pt = c2PtLE*LE + c2PtE*E + c2PtEm*Em + c2PtP*P + c2PtPm*Pm + c2PtfE*fE + c2PtfP*fP
      !   LE = c^2Delta(E) 
      c2PtLE=.5*a1*dt**2/d1
      c2PtE=(1.*a0*dt+1.*a1)/d1
      c2PtEm=-1.*a1/d1
      c2PtP=(alpha*a1*dt-b0*dt**2+2.)/dt/d1
      c2PtPm=(-alpha*a1*dt-2.)/dt/d1
      c2PtfE=.5*a1*dt**2/d1
      c2PtfP=1.*dt/d1
      ! --------------- Here is Ptt to second-order ------------

      ! Ptt = a1*dt/d1*LE+((-a0*a1*alpha-a0*b1)*dt^2+d1*a0*dt+2*a1)/dt/d1*E-2*a1/dt/d1*Em+((a1*alpha*b0+b0*b1)*dt^2-d1*b0*dt-2*b1)/dt/d1*P+2*b1/dt/d1*Pm+a1*dt/d1*fE+((-a1*alpha-b1)*dt^2+dt*d1)/dt/d1*fP
      ! Ptt = c2PttLE*LE + c2PttE*E + c2PttEm*Em + c2PttP*P + c2PttPm*Pm + c2PttfE*fE + c2PttfP*fP
      c2PttLE=1.*a1*dt/d1
      c2PttE=((-1.*a1*alpha-1.*b1)*a0*dt**2+d1*a0*dt+2.*a1)/dt/d1
      c2PttEm=-2.*a1/dt/d1
      c2PttP=(b0*(a1*alpha+b1)*dt**2-1.*d1*b0*dt-2.*b1)/dt/d1
      c2PttPm=2.*b1/dt/d1
      c2PttfE=1.*a1*dt/d1
      c2PttfP=(1.*d1+(-1.*a1*alpha-1.*b1)*dt)/d1
                        ! Fourth-order coefficients
                        ! Ptt = c4PttLE*LE2 + c4PttE*E + c4PttEm*Em + c4PttP*P + c4PttPm*Pm + c4PttfE*fE + c4PttfP*fP
                        !    + c4PttLLE*LLE2 + c4PttLP*LP + c4PttLEm*LE2m + c4PttLPm*LPm+ c4PttLfE*LfE2 + c4PttLfP*LfP2
                        !    + c4PttfEt*fEt2 + c4PttfEtt*fEtt2 + c4PttfPt*fPt2+ c4PttfPtt*fPtt2
                        ! #Include interfaceAdeGdmOrder4.h 
                        ! Nov 4, 2018 *new way* 
      ! File created by Dropbox/GDM/maple/interfaceNew.maple
      ! 
      ! --------------- Here is Ptt to fourth-order ------------

      ! Ptt = ((1/2*(-b0*c2PttE*alpha*dt^6-6*c2PttE*alpha*b1*dt^5+(12*a1^2*c2PttLE*alpha^3+24*a1*(c2PttLE*b1-1/2*a1+1/2*c2PtLE*b0-1/2*a0*c2EtLE)*alpha^2+(12*c2PttLE*b1^2+(-12*a0*c2EtLE+12*b0*c2PtLE-12*a1)*b1+12*c2EtE*a1-12*c2PttE)*alpha+12*b0)*dt^4+72*b1*dt^3+144*dt^2)/d4/dt-1/6*dt^2*(-(a0*c2EtLE+(-alpha*c2PttLE+1)*a1-c2PtLE*b0-c2PttLE*b1)*alpha+c2EtE))*a1-(1/2*(-a0*c2PttE*alpha*dt^6-6*c2PttE*alpha*a1*dt^5+(-12*(alpha*c2PttLE-1)*alpha*a1^2+((12*a0*c2EtLE-12*b0*c2PtLE-24*b1*c2PttLE)*alpha+12*b1-12*c2EtE)*a1-12*c2PttLE*b1^2+(12*a0*c2EtLE-12*b0*c2PtLE)*b1+12*a0)*dt^4+72*a1*dt^3)/d4/dt-1/6*dt^2*(a0*c2EtLE+(-alpha*c2PttLE+1)*a1-c2PtLE*b0-c2PttLE*b1))*b1)*LE+((1/2*(-b0*(alpha*c2PttLE-1)*dt^6-6*(alpha*c2PttLE-1)*b1*dt^5+((12*a1*c2EtLE-12*c2PttLE)*alpha+12)*dt^4)/d4/dt-1/6*dt^2*c2EtLE)*a1-1/2*(-a0*(alpha*c2PttLE-1)*dt^6-6*(alpha*c2PttLE-1)*a1*dt^5-12*c2EtLE*a1*dt^4)/d4/dt*b1)*LLE+((1/2*(-b0*c2PttEm*alpha*dt^6-6*c2PttEm*alpha*b1*dt^5+(12*a1*c2EtEm-12*c2PttEm)*alpha*dt^4)/d4/dt-1/6*dt^2*c2EtEm)*a1-1/2*(-a0*alpha*c2PttEm*dt^6-6*a1*alpha*c2PttEm*dt^5-12*a1*c2EtEm*dt^4)/d4/dt*b1)*LEm+((1/2*(-b0*c2PttP*alpha*dt^6-6*c2PttP*alpha*b1*dt^5+(12*a1*c2EtP-12*c2PttP)*alpha*dt^4)/d4/dt-1/6*dt^2*c2EtP)*a1-1/2*(-a0*alpha*c2PttP*dt^6-6*a1*alpha*c2PttP*dt^5-12*a1*c2EtP*dt^4)/d4/dt*b1)*LP+((1/2*(-b0*c2PttPm*alpha*dt^6-6*c2PttPm*alpha*b1*dt^5+(12*a1*c2EtPm-12*c2PttPm)*alpha*dt^4)/d4/dt-1/6*dt^2*c2EtPm)*a1-1/2*(-a0*alpha*c2PttPm*dt^6-6*a1*alpha*c2PttPm*dt^5-12*a1*c2EtPm*dt^4)/d4/dt*b1)*LPm+(a0+(1/2*((12*a1^2*c2PttE*alpha^3+24*a1*(c2PttE*b1+1/2*c2PtE*b0-1/2*c2EtE*a0)*alpha^2+(12*c2PttE*b1^2+(-12*a0*c2EtE+12*b0*c2PtE)*b1)*alpha)*dt^4+(-120*a0*alpha+24*b0)*dt^2+144*b1*dt+288)/d4/dt+1/6*dt^2*(-a1*alpha*c2PttE+a0*c2EtE-b0*c2PtE-b1*c2PttE)*alpha)*a1-(1/2*((-12*c2PttE*alpha^2*a1^2+(12*a0*c2EtE-12*b0*c2PtE-24*b1*c2PttE)*alpha*a1-12*c2PttE*b1^2+(12*a0*c2EtE-12*b0*c2PtE)*b1)*dt^4+144*a0*dt^2+144*a1*dt)/d4/dt-1/6*dt^2*(-a1*alpha*c2PttE+a0*c2EtE-b0*c2PtE-b1*c2PttE))*b1)*E+((1/2*(((12*a1^2*c2PttEm*alpha^3+24*a1*(c2PttEm*b1+1/2*c2PtEm*b0-1/2*c2EtEm*a0)*alpha^2+(12*c2PttEm*b1^2+(-12*a0*c2EtEm+12*b0*c2PtEm)*b1)*alpha)*dt^4+(-12*a0*alpha-12*b0)*dt^2+(72*a1*alpha-72*b1)*dt-144)/d4-1)/dt+1/6*dt^2*(-a1*alpha*c2PttEm+a0*c2EtEm-b0*c2PtEm-b1*c2PttEm)*alpha)*a1-(1/2*((-12*c2PttEm*alpha^2*a1^2+(12*a0*c2EtEm-12*b0*c2PtEm-24*b1*c2PttEm)*alpha*a1-12*c2PttEm*b1^2+(12*a0*c2EtEm-12*b0*c2PtEm)*b1)*dt^4-144*a1*dt)/d4/dt-1/6*dt^2*(-a1*alpha*c2PttEm+a0*c2EtEm-b0*c2PtEm-b1*c2PttEm))*b1)*Em+((1/2*((12*a1^2*c2PttP*alpha^3+24*a1*(c2PttP*b1+1/2*c2PtP*b0-1/2*c2EtP*a0)*alpha^2+(12*c2PttP*b1^2+(-12*a0*c2EtP+12*b0*c2PtP)*b1)*alpha)*dt^4+144*b0*alpha*dt^2+144*b1*alpha*dt)/d4/dt+1/6*dt^2*(-a1*alpha*c2PttP+a0*c2EtP-b0*c2PtP-b1*c2PttP)*alpha)*a1-b0-(1/2*((-12*c2PttP*alpha^2*a1^2+(12*a0*c2EtP-12*b0*c2PtP-24*b1*c2PttP)*alpha*a1-12*c2PttP*b1^2+(12*a0*c2EtP-12*b0*c2PtP)*b1)*dt^4+(24*a0*alpha-120*b0)*dt^2+144*alpha*a1*dt+288)/d4/dt-1/6*dt^2*(-a1*alpha*c2PttP+a0*c2EtP-b0*c2PtP-b1*c2PttP))*b1)*P+((1/2*((12*a1^2*c2PttPm*alpha^3+24*a1*(c2PttPm*b1+1/2*c2PtPm*b0-1/2*c2EtPm*a0)*alpha^2+(12*c2PttPm*b1^2+(-12*a0*c2EtPm+12*b0*c2PtPm)*b1)*alpha)*dt^4-144*b1*alpha*dt)/d4/dt+1/6*dt^2*(-a1*alpha*c2PttPm+a0*c2EtPm-b0*c2PtPm-b1*c2PttPm)*alpha)*a1-(1/2*(((-12*c2PttPm*alpha^2*a1^2+(12*a0*c2EtPm-12*b0*c2PtPm-24*b1*c2PttPm)*alpha*a1-12*c2PttPm*b1^2+(12*a0*c2EtPm-12*b0*c2PtPm)*b1)*dt^4+(-12*a0*alpha-12*b0)*dt^2+(-72*a1*alpha+72*b1)*dt-144)/d4-1)/dt-1/6*dt^2*(-a1*alpha*c2PttPm+a0*c2EtPm-b0*c2PtPm-b1*c2PttPm))*b1)*Pm+((1/2*((12*a1^2*c2PttfE*alpha^3+24*a1*(c2PttfE*b1-1/2*a1+1/2*c2PtfE*b0-1/2*c2EtfE*a0)*alpha^2+(12*c2PttfE*b1^2+(-12*a0*c2EtfE+12*b0*c2PtfE-12*a1)*b1)*alpha+12*b0)*dt^4+72*b1*dt^3+144*dt^2)/d4/dt+1/6*dt^2*(c2EtfE*a0+(-alpha*c2PttfE+1)*a1-c2PtfE*b0-c2PttfE*b1)*alpha)*a1-(1/2*((-12*(alpha*c2PttfE-1)*alpha*a1^2+((12*a0*c2EtfE-12*b0*c2PtfE-24*b1*c2PttfE)*alpha+12*b1)*a1-12*c2PttfE*b1^2+(12*a0*c2EtfE-12*b0*c2PtfE)*b1+12*a0)*dt^4+72*a1*dt^3)/d4/dt-1/6*dt^2*(c2EtfE*a0+(-alpha*c2PttfE+1)*a1-c2PtfE*b0-c2PttfE*b1))*b1)*fE+((1/2*(-b0*(alpha*c2PttfE-1)*dt^6-6*(alpha*c2PttfE-1)*b1*dt^5+((12*a1*c2EtfE-12*c2PttfE)*alpha+12)*dt^4)/d4/dt-1/6*dt^2*c2EtfE)*a1-1/2*(-a0*(alpha*c2PttfE-1)*dt^6-6*(alpha*c2PttfE-1)*a1*dt^5-12*c2EtfE*a1*dt^4)/d4/dt*b1)*LfE+((6*alpha*a1*dt^3/d4-1/6*dt^2)*a1+6*a1*dt^3/d4*b1)*fEt+(1/2*(b0*dt^6+6*b1*dt^5+12*dt^4)/d4/dt*a1-1/2*(a0*dt^6+6*a1*dt^5)/d4/dt*b1)*fEtt+((1/2*((12*a1^2*alpha^3*c2PttfP+24*a1*(b1*c2PttfP+1/2*b0*c2PtfP-1/2*a0*c2EtfP)*alpha^2+(12*b1^2*c2PttfP+(-12*a0*c2EtfP+12*b0*c2PtfP)*b1)*alpha)*dt^4-144*alpha*dt^2)/d4/dt+1/6*dt^2*(-a1*alpha*c2PttfP+a0*c2EtfP-b0*c2PtfP-b1*c2PttfP)*alpha)*a1-(1/2*((-12*a1^2*alpha^2*c2PttfP+(12*a0*c2EtfP-12*b0*c2PtfP-24*b1*c2PttfP)*alpha*a1-12*b1^2*c2PttfP+(12*a0*c2EtfP-12*b0*c2PtfP)*b1)*dt^4+144*dt^2)/d4/dt-1/6*dt^2*(-a1*alpha*c2PttfP+a0*c2EtfP-b0*c2PtfP-b1*c2PttfP))*b1+1)*fP+((1/2*(-alpha*b0*c2PttfP*dt^6-6*alpha*b1*c2PttfP*dt^5+(12*a1*c2EtfP-12*c2PttfP)*alpha*dt^4)/d4/dt-1/6*dt^2*c2EtfP)*a1-1/2*(-a0*alpha*c2PttfP*dt^6-6*a1*alpha*c2PttfP*dt^5-12*a1*c2EtfP*dt^4)/d4/dt*b1)*LfP+((1/2*(-12*a1*alpha^2-12*alpha*b1)*dt^3/d4+1/6*alpha*dt^2)*a1-(1/2*(12*a1*alpha+12*b1)*dt^3/d4-1/6*dt^2)*b1)*fPt+(-6*alpha*a1*dt^3/d4-6*dt^3/d4*b1)*fPtt
      ! Ptt = c4PttLE*LE + c4PttE*E + c4PttEm*Em + c4PttP*P + c4PttPm*Pm + c4PttfE*fE + c4PttfP*fP
      !    + c4PttLLE*LLE + c4PttLP*LP + c4PttLEm*LEm + c4PttLPm*LPm+ c4PttLfE*LfE + c4PttLfP*LfP
      !    + c4PttfEt*fEt + c4PttfEtt*fEtt + c4PttfPt*fPt+ c4PttfPtt*fPtt
      d4=144.+(12.*a0*alpha+12.*b0)*dt**2+(72.*alpha*a1+72.*b1)*dt

      c4PttLE=-.166666666666666667*(-432.*a1-3.*a0*alpha*b1*c2PttE*dt**4+3.*a1*alpha*b0*c2PttE*dt**4+36.*a0*a1**2*alpha**2*c2EtLE*dt**2-36.*a1**2*alpha**2*b0*c2PtLE*dt**2-108.*a1**2*alpha**2*b1*c2PttLE*dt**2-108.*a1*alpha*b1**2*c2PttLE*dt**2+a1**2*alpha**2*c2PttLE*d4*dt-a0*b1*c2EtLE*d4*dt+b0*b1*c2PtLE*d4*dt-36.*b1**3*c2PttLE*dt**2-36.*a1*b0*dt**2+36.*a1**3*alpha**2*dt**2+36.*a1*b1**2*dt**2-36.*a1**3*alpha**3*c2PttLE*dt**2+36.*a0*b1**2*c2EtLE*dt**2+72.*a1**2*alpha*b1*dt**2-36.*a1**2*alpha*c2EtE*dt**2-36.*b0*b1**2*c2PtLE*dt**2+36.*a1*alpha*c2PttE*dt**2-36.*a1*b1*c2EtE*dt**2-a1**2*alpha*d4*dt+b1**2*c2PttLE*d4*dt-a1*b1*d4*dt+a1*c2EtE*d4*dt+72.*a0*a1*alpha*b1*c2EtLE*dt**2-72.*a1*alpha*b0*b1*c2PtLE*dt**2-a0*a1*alpha*c2EtLE*d4*dt+a1*alpha*b0*c2PtLE*d4*dt+2.*a1*alpha*b1*c2PttLE*d4*dt+36.*a0*b1*dt**2)*dt/d4
      c4PttLLE=.166666666666666667*dt**2*(3.*a0*alpha*b1*c2PttLE*dt**3-3.*a1*alpha*b0*c2PttLE*dt**3-3.*a0*b1*dt**3+36.*a1**2*alpha*c2EtLE*dt+3.*a1*b0*dt**3-36.*a1*alpha*c2PttLE*dt+36.*a1*b1*c2EtLE*dt-a1*c2EtLE*d4+36.*a1*dt)/d4
      c4PttE=((6.00000000000000001*c2PttE*a1**3*alpha**3+(-6.00000000000000001*c2EtE*a0+6.00000000000000001*c2PtE*b0+18.*c2PttE*b1)*alpha**2*a1**2+(18.*c2PttE*b1**2+(-12.*c2EtE*a0+12.*c2PtE*b0)*b1)*alpha*a1+6.00000000000000001*c2PttE*b1**3+(-6.00000000000000001*c2EtE*a0+6.00000000000000001*c2PtE*b0)*b1**2)*dt**4+(-.166666666666666667*c2PttE*a1**2*d4*alpha**2+(-.333333333333333334*c2PttE*b1*d4+(.166666666666666667*c2EtE*a0-.166666666666666667*c2PtE*b0)*d4)*alpha*a1-.166666666666666667*c2PttE*b1**2*d4+(.166666666666666667*c2EtE*a0-.166666666666666667*c2PtE*b0)*d4*b1)*dt**3+((-60.0000000000000001*a0*alpha+12.*b0)*a1-72.0000000000000001*a0*b1)*dt**2+a0*d4*dt+144.*a1)/d4/dt
      c4PttEm=((6.00000000000000001*c2PttEm*a1**3*alpha**3+(-6.00000000000000001*c2EtEm*a0+6.00000000000000001*c2PtEm*b0+18.*c2PttEm*b1)*alpha**2*a1**2+(18.*c2PttEm*b1**2+(-12.*c2EtEm*a0+12.*c2PtEm*b0)*b1)*alpha*a1+6.00000000000000001*c2PttEm*b1**3+(-6.00000000000000001*c2EtEm*a0+6.00000000000000001*c2PtEm*b0)*b1**2)*dt**4+(-.166666666666666667*c2PttEm*a1**2*d4*alpha**2+(-.333333333333333334*c2PttEm*b1*d4+(.166666666666666667*c2EtEm*a0-.166666666666666667*c2PtEm*b0)*d4)*alpha*a1-.166666666666666667*c2PttEm*b1**2*d4+(.166666666666666667*c2EtEm*a0-.166666666666666667*c2PtEm*b0)*d4*b1)*dt**3+(-6.00000000000000001*a0*alpha-6.00000000000000001*b0)*a1*dt**2+(36.0000000000000001*a1**2*alpha+36.0000000000000001*a1*b1)*dt+(-.500000000000000001*d4-72.0000000000000001)*a1)/d4/dt
      c4PttP=((6.00000000000000001*c2PttP*b1**3+(18.*a1*alpha*c2PttP-6.00000000000000001*c2EtP*a0+6.00000000000000001*c2PtP*b0)*b1**2+(18.*c2PttP*alpha**2*a1**2+(-12.*c2EtP*a0+12.*c2PtP*b0)*a1*alpha)*b1+6.00000000000000001*c2PttP*a1**3*alpha**3+(-6.00000000000000001*c2EtP*a0+6.00000000000000001*c2PtP*b0)*a1**2*alpha**2)*dt**4+(-.166666666666666667*c2PttP*b1**2*d4+(-.333333333333333334*c2PttP*a1*d4*alpha+(.166666666666666667*c2EtP*a0-.166666666666666667*c2PtP*b0)*d4)*b1-.166666666666666667*c2PttP*a1**2*d4*alpha**2+(.166666666666666667*c2EtP*a0-.166666666666666667*c2PtP*b0)*d4*a1*alpha)*dt**3+((-12.*a0*alpha+60.0000000000000001*b0)*b1+72.0000000000000001*a1*b0*alpha)*dt**2-1.*b0*d4*dt-144.*b1)/d4/dt
      c4PttPm=((6.00000000000000001*c2PttPm*b1**3+(18.*a1*alpha*c2PttPm-6.00000000000000001*c2EtPm*a0+6.00000000000000001*c2PtPm*b0)*b1**2+(18.*c2PttPm*alpha**2*a1**2+(-12.*c2EtPm*a0+12.*c2PtPm*b0)*a1*alpha)*b1+6.00000000000000001*c2PttPm*a1**3*alpha**3+(-6.00000000000000001*c2EtPm*a0+6.00000000000000001*c2PtPm*b0)*a1**2*alpha**2)*dt**4+(-.166666666666666667*c2PttPm*b1**2*d4+(-.333333333333333334*c2PttPm*a1*d4*alpha+(.166666666666666667*c2EtPm*a0-.166666666666666667*c2PtPm*b0)*d4)*b1-.166666666666666667*c2PttPm*a1**2*d4*alpha**2+(.166666666666666667*c2EtPm*a0-.166666666666666667*c2PtPm*b0)*d4*a1*alpha)*dt**3+(6.00000000000000001*a0*alpha+6.00000000000000001*b0)*b1*dt**2+(-36.0000000000000001*a1*b1*alpha-36.0000000000000001*b1**2)*dt+(.500000000000000001*d4+72.0000000000000001)*b1)/d4/dt
      c4PttLP=.166666666666666667*dt**2*(3.*a0*alpha*b1*c2PttP*dt**3-3.*a1*alpha*b0*c2PttP*dt**3+36.*a1**2*alpha*c2EtP*dt-36.*a1*alpha*c2PttP*dt+36.*a1*b1*c2EtP*dt-a1*c2EtP*d4)/d4
      c4PttLEm=.166666666666666667*dt**2*(3.*a0*alpha*b1*c2PttEm*dt**3-3.*a1*alpha*b0*c2PttEm*dt**3+36.*a1**2*alpha*c2EtEm*dt-36.*a1*alpha*c2PttEm*dt+36.*a1*b1*c2EtEm*dt-a1*c2EtEm*d4)/d4
      c4PttLPm=.166666666666666667*dt**2*(3.*a0*alpha*b1*c2PttPm*dt**3-3.*a1*alpha*b0*c2PttPm*dt**3+36.*a1**2*alpha*c2EtPm*dt-36.*a1*alpha*c2PttPm*dt+36.*a1*b1*c2EtPm*dt-a1*c2EtPm*d4)/d4
      c4PttfE=-.166666666666666667*dt*(-432.*a1-36.*a1*b0*dt**2+36.*a1**3*alpha**2*dt**2+36.*a1*b1**2*dt**2+72.*a1**2*alpha*b1*dt**2-a1**2*alpha*d4*dt-a1*b1*d4*dt-36.*b1**3*c2PttfE*dt**2+72.*a0*a1*alpha*b1*c2EtfE*dt**2-72.*a1*alpha*b0*b1*c2PtfE*dt**2-a0*a1*alpha*c2EtfE*d4*dt+a1*alpha*b0*c2PtfE*d4*dt+2.*a1*alpha*b1*c2PttfE*d4*dt+36.*a0*a1**2*alpha**2*c2EtfE*dt**2-36.*a1**2*alpha**2*b0*c2PtfE*dt**2-108.*a1**2*alpha**2*b1*c2PttfE*dt**2+a1**2*alpha**2*c2PttfE*d4*dt-108.*a1*alpha*b1**2*c2PttfE*dt**2-a0*b1*c2EtfE*d4*dt+b0*b1*c2PtfE*d4*dt-36.*a1**3*alpha**3*c2PttfE*dt**2+36.*a0*b1**2*c2EtfE*dt**2-36.*b0*b1**2*c2PtfE*dt**2+b1**2*c2PttfE*d4*dt+36.*a0*b1*dt**2)/d4
      c4PttfP=((6.00000000000000001*c2PttfP*a1**3*alpha**3+(-6.00000000000000001*a0*c2EtfP+6.00000000000000001*b0*c2PtfP+18.*b1*c2PttfP)*alpha**2*a1**2+(18.*b1**2*c2PttfP+(-12.*a0*c2EtfP+12.*b0*c2PtfP)*b1)*alpha*a1+6.00000000000000001*c2PttfP*b1**3+(-6.00000000000000001*a0*c2EtfP+6.00000000000000001*b0*c2PtfP)*b1**2)*dt**3+(-.166666666666666667*c2PttfP*a1**2*d4*alpha**2+(-.333333333333333334*c2PttfP*b1*d4+(.166666666666666667*a0*c2EtfP-.166666666666666667*b0*c2PtfP)*d4)*alpha*a1-.166666666666666667*c2PttfP*b1**2*d4+(.166666666666666667*a0*c2EtfP-.166666666666666667*b0*c2PtfP)*d4*b1)*dt**2+(-72.0000000000000001*alpha*a1-72.0000000000000001*b1)*dt+d4)/d4
      c4PttfEt=.166666666666666667*a1*dt**2*(36.*alpha*a1*dt+36.*b1*dt-d4)/d4
      c4PttfPt=-.166666666666666667*dt**2*(36.*alpha*a1*dt+36.*b1*dt-d4)*(a1*alpha+b1)/d4
      c4PttLfE=.166666666666666667*dt**2*(3.*a0*alpha*b1*c2PttfE*dt**3-3.*a1*alpha*b0*c2PttfE*dt**3-3.*a0*b1*dt**3+36.*a1**2*alpha*c2EtfE*dt+3.*a1*b0*dt**3-36.*a1*alpha*c2PttfE*dt+36.*a1*b1*c2EtfE*dt-a1*c2EtfE*d4+36.*a1*dt)/d4
      c4PttLfP=.166666666666666667*dt**2*(3.*a0*alpha*b1*c2PttfP*dt**3-3.*a1*alpha*b0*c2PttfP*dt**3+36.*a1**2*alpha*c2EtfP*dt-36.*a1*alpha*c2PttfP*dt+36.*a1*b1*c2EtfP*dt-a1*c2EtfP*d4)/d4
      c4PttfEtt=(-.5*a0*b1*dt**2+.5*a1*b0*dt**2+6.0*a1)*dt**3/d4
      c4PttfPtt=(-6.*alpha*a1-6.*b1)*dt**3/d4

      ! --------------- Here is Ptttt to second-order ------------

      ! Ptttt = ((-alpha*c2PttLE+1)*a0+(-(a0*c2EtLE+(-alpha*c2PttLE+1)*a1-c2PtLE*b0-c2PttLE*b1)*alpha+c2EtE)*a1-c2PttLE*b0-(a0*c2EtLE+(-alpha*c2PttLE+1)*a1-c2PtLE*b0-c2PttLE*b1)*b1)*LE+c2EtLE*a1*LLE+c2EtEm*a1*LEm+c2EtP*a1*LP+c2EtPm*a1*LPm+(-c2PttE*alpha*a0-(-a1*alpha*c2PttE+a0*c2EtE-b0*c2PtE-b1*c2PttE)*alpha*a1-c2PttE*b0-(-a1*alpha*c2PttE+a0*c2EtE-b0*c2PtE-b1*c2PttE)*b1)*E+(-alpha*c2PttEm*a0-(-a1*alpha*c2PttEm+a0*c2EtEm-b0*c2PtEm-b1*c2PttEm)*alpha*a1-c2PttEm*b0-(-a1*alpha*c2PttEm+a0*c2EtEm-b0*c2PtEm-b1*c2PttEm)*b1)*Em+(-alpha*c2PttP*a0-(-a1*alpha*c2PttP+a0*c2EtP-b0*c2PtP-b1*c2PttP)*alpha*a1-c2PttP*b0-(-a1*alpha*c2PttP+a0*c2EtP-b0*c2PtP-b1*c2PttP)*b1)*P+(-alpha*c2PttPm*a0-(-a1*alpha*c2PttPm+a0*c2EtPm-b0*c2PtPm-b1*c2PttPm)*alpha*a1-c2PttPm*b0-(-a1*alpha*c2PttPm+a0*c2EtPm-b0*c2PtPm-b1*c2PttPm)*b1)*Pm+((-alpha*c2PttfE+1)*a0-(c2EtfE*a0+(-alpha*c2PttfE+1)*a1-c2PtfE*b0-c2PttfE*b1)*alpha*a1-c2PttfE*b0-(c2EtfE*a0+(-alpha*c2PttfE+1)*a1-c2PtfE*b0-c2PttfE*b1)*b1)*fE+c2EtfE*a1*LfE+a1*fEt+(-c2PttfP*alpha*a0-(-a1*alpha*c2PttfP+a0*c2EtfP-b0*c2PtfP-b1*c2PttfP)*alpha*a1-c2PttfP*b0-(-a1*alpha*c2PttfP+a0*c2EtfP-b0*c2PtfP-b1*c2PttfP)*b1)*fP+c2EtfP*a1*LfP+(-a1*alpha-b1)*fPt+fPtt
      ! Ptttt = c2PttLE*LE + c2PttE*E + c2PttEm*Em + c2PttP*P + c2PttPm*Pm + c2PttfE*fE + c2PttfP*fP
      !      + c2PttLLE*LLE + c2PttLP*LP + c2PttLEm*LEm + c2PttLPm*LPm+ c2PttLfE*LfE + c2PttLfP*LfP
      !      + c2PttfEt*fEt + c2PttfEtt*fEtt + c2PttfPt*fPt+ c2PttfPtt*fPtt
      c2PttttLE=(alpha**2*c2PttLE-1.*alpha)*a1**2+((-1.*a0*c2EtLE+c2PtLE*b0+2.*c2PttLE*b1)*alpha-1.*b1+c2EtE)*a1-1.*a0*alpha*c2PttLE+(b1**2-1.*b0)*c2PttLE+(-1.*a0*c2EtLE+c2PtLE*b0)*b1+a0
      c2PttttLLE=1.*c2EtLE*a1
      c2PttttE=(alpha**2*a1**2+(2.*a1*b1-1.*a0)*alpha+b1**2-1.*b0)*c2PttE+(-1.*c2EtE*a0+c2PtE*b0)*a1*alpha+(-1.*c2EtE*a0+c2PtE*b0)*b1
      c2PttttEm=(alpha**2*a1**2+(2.*a1*b1-1.*a0)*alpha+b1**2-1.*b0)*c2PttEm+(-1.*c2EtEm*a0+c2PtEm*b0)*a1*alpha+(-1.*c2EtEm*a0+c2PtEm*b0)*b1
      c2PttttP=(alpha**2*a1**2+(2.*a1*b1-1.*a0)*alpha+b1**2-1.*b0)*c2PttP+(-1.*c2EtP*a0+c2PtP*b0)*a1*alpha+(-1.*c2EtP*a0+c2PtP*b0)*b1
      c2PttttPm=(alpha**2*a1**2+(2.*a1*b1-1.*a0)*alpha+b1**2-1.*b0)*c2PttPm+(-1.*c2EtPm*a0+c2PtPm*b0)*a1*alpha+(-1.*c2EtPm*a0+c2PtPm*b0)*b1
      c2PttttLP=1.*c2EtP*a1
      c2PttttLEm=1.*c2EtEm*a1
      c2PttttLPm=1.*c2EtPm*a1
      c2PttttfE=a1**2*alpha**2*c2PttfE+(-1.*a1**2+(-1.*c2EtfE*a0+c2PtfE*b0+2.*c2PttfE*b1)*a1-1.*a0*c2PttfE)*alpha-1.*a1*b1+(b1**2-1.*b0)*c2PttfE+(-1.*c2EtfE*a0+c2PtfE*b0)*b1+a0
      c2PttttfP=(alpha**2*a1**2+(2.*a1*b1-1.*a0)*alpha+b1**2-1.*b0)*c2PttfP+(-1.*a0*c2EtfP+b0*c2PtfP)*a1*alpha+(-1.*a0*c2EtfP+b0*c2PtfP)*b1
      c2PttttfEt=1.*a1
      c2PttttfPt=-1.*alpha*a1-1.*b1
      c2PttttLfE=1.*c2EtfE*a1
      c2PttttLfP=1.*c2EtfP*a1
      c2PttttfEtt=0.
      c2PttttfPtt=1.
                        if( .false. .and. twilightZone.eq.1 )then
                          write(*,'(" RIGHT: alpha,dt,d4,d1=",4e12.4)') alpha,dt,d4,d1
                          write(*,'(" a0,a1,b0,b1=",4e12.4)') a0,a1,b0,b1
                          write(*,'(" c2EtE,c2PtE,c2PttfP=",3e12.4)') c2EtE,c2PtE,c2PttfP
                          write(*,'(" c4PttLE,c4PttE,c4PttEm,c4PttP,c4PttPm,c4PttfE=",6e12.4)') c4PttLE,c4PttE,c4PttEm,c4PttP,c4PttPm,c4PttfE
                          write(*,'(" c4PttfP,c4PttLLE,c4PttLP,c4PttLEm,c4PttLPm,c4PttLfE,c4PttLfP=",7e12.4)') c4PttfP,c4PttLLE,c4PttLP,c4PttLEm,c4PttLPm,c4PttLfE,c4PttLfP
                          write(*,'(" c4PttfEt,c4PttfEtt,c4PttfPt,c4PttfPtt=",4e12.4)') c4PttfEt,c4PttfEtt,c4PttfPt,c4PttfPtt
                        end if
                        ! Coeff of E in P.tt (4th order)
                        c2PttEsum2 = c2PttEsum2 + c2PttE
                        ! Coeff of LE2 in P.tt (4th order)
                        c2PttLEsum2  = c2PttLEsum2 + c2PttLE
                        ! Coeff of LE2 in P.tt (4th order)
                        c4PttLEsum2  = c4PttLEsum2 + c4PttLE
                        ! Coeff of LLE2 in P.tt
                        c4PttLLEsum2 = c4PttLLEsum2 + c4PttLLE
                        ! Coeff of LE2 and LLE2 in P.tttt
                        c2PttttLEsum2 =c2PttttLEsum2 +c2PttttLE
                        c2PttttLLEsum2=c2PttttLLEsum2+c2PttttLLE
                        do n=0,nd-1
                          pc = n + jv*nd 
                          ec = ex +n
                          pv   =  p2(j1,j2,j3,pc)
                          pvn  =  p2n(j1,j2,j3,pc)
                          ev    =  u2(j1,j2,j3,ec)
                          evn   =  u2n(j1,j2,j3,ec)
                          ! Left: u1x,u1y, u1xx, u1yy, u1Lap (ex)
                          !       v1x,v1y, v1xx, v1yy, v1Lap (ey) 
                          ! These next derivatives may only be needed to order2, (use order 4 for testing TZ polynomials)
                          ! perl $ORDER=2;
                            ! uu1 in the next statement defines names of intermediate values
                              uu2=p2(j1,j2,j3,pc) ! in the rectangular case just eval the solution
                               p2x = (p2(j1-2,j2,j3,pc)-8.*p2(j1-1,j2,j3,pc)+8.*p2(j1+1,j2,j3,pc)-p2(j1+2,j2,j3,pc))/(12.*dx1(0))
                               p2y = (p2(j1,j2-2,j3,pc)-8.*p2(j1,j2-1,j3,pc)+8.*p2(j1,j2+1,j3,pc)-p2(j1,j2+2,j3,pc))/(12.*dx1(1))
                               p2xx = (-p2(j1-2,j2,j3,pc)+16.*p2(j1-1,j2,j3,pc)-30.*p2(j1,j2,j3,pc)+16.*p2(j1+1,j2,j3,pc)-p2(j1+2,j2,j3,pc))/(12.*dx1(0)**2)
                               p2yy = (-p2(j1,j2-2,j3,pc)+16.*p2(j1,j2-1,j3,pc)-30.*p2(j1,j2,j3,pc)+16.*p2(j1,j2+1,j3,pc)-p2(j1,j2+2,j3,pc))/(12.*dx1(1)**2)
                             p2Lap = p2xx+ p2yy
                            LP  = (c2**2)*p2Lap
                              uu2=p2n(j1,j2,j3,pc) ! in the rectangular case just eval the solution
                               p2nx = (p2n(j1-2,j2,j3,pc)-8.*p2n(j1-1,j2,j3,pc)+8.*p2n(j1+1,j2,j3,pc)-p2n(j1+2,j2,j3,pc))/(12.*dx1(0))
                               p2ny = (p2n(j1,j2-2,j3,pc)-8.*p2n(j1,j2-1,j3,pc)+8.*p2n(j1,j2+1,j3,pc)-p2n(j1,j2+2,j3,pc))/(12.*dx1(1))
                               p2nxx = (-p2n(j1-2,j2,j3,pc)+16.*p2n(j1-1,j2,j3,pc)-30.*p2n(j1,j2,j3,pc)+16.*p2n(j1+1,j2,j3,pc)-p2n(j1+2,j2,j3,pc))/(12.*dx1(0)**2)
                               p2nyy = (-p2n(j1,j2-2,j3,pc)+16.*p2n(j1,j2-1,j3,pc)-30.*p2n(j1,j2,j3,pc)+16.*p2n(j1,j2+1,j3,pc)-p2n(j1,j2+2,j3,pc))/(12.*dx1(1)**2)
                             p2nLap = p2nxx+ p2nyy
                            LPm = (c2**2)*p2nLap
                            pvx  = p2x
                            pvy  = p2y
                            pvnx = p2nx
                            pvny = p2ny
                          ! Accumulate: SUM_m Pm,tt
                          ptta = c4PttLE*LE2(n) + c4PttE*ev + c4PttEm*evn + c4PttP*pv + c4PttPm*pvn + c4PttfE*fev2(n) + c4PttfP*fpv2(n,jv) + c4PttLLE*LLE2(n) + c4PttLP*LP + c4PttLEm*LE2m(n) + c4PttLPm*LPm+ c4PttLfE*LfE2(n) + c4PttLfP*LfP2(n,jv) + c4PttfEt*fEt2(n) + c4PttfEtt*fEtt2(n) + c4PttfPt*fPt2(n,jv)+ c4PttfPtt*fPtt2(n,jv)
                          ! ---- Compute fp2 = P.tt 
                          fp2(n) = fp2(n) + ptta
                          ! ----- Compute fPttx2 = (P.tt).x , fPtty2 = (P.tt).y  (second order)
                          pttxa = c2PttLE*LEx2(n) + c2PttE*evx2(n) + c2PttEm*evnx2(n) + c2PttP*pvx + c2PttPm*pvnx + c2PttfE*fevx2(n) + c2PttfP*fpvx2(n,jv)
                          pttya = c2PttLE*LEy2(n) + c2PttE*evy2(n) + c2PttEm*evny2(n) + c2PttP*pvy + c2PttPm*pvny + c2PttfE*fevy2(n) + c2PttfP*fpvy2(n,jv)
                          fPttx2(n) = fPttx2(n) + pttxa
                          fPtty2(n) = fPtty2(n) + pttya
                          ! ----- Compute fLPtt2 = L(P.tt) (second order)
                          Lptta = c2PttLE*LLE2(n) + c2PttE*LE2(n) + c2PttEm*LE2m(n) + c2PttP*LP + c2PttPm*LPm + c2PttfE*LfE2(n) + c2PttfP*LfP2(n,jv)
                          fLPtt2(n) = fLPtt2(n) + Lptta
                          ! ----- Compute fPtttt2 = P.tttt
                          ptttta= c2PttttLE*LE2(n) + c2PttttE*ev + c2PttttEm*evn + c2PttttP*pv + c2PttttPm*pvn + c2PttttfE*fev2(n) + c2PttttfP*fpv2(n,jv) + c2PttttLLE*LLE2(n) + c2PttttLP*LP + c2PttttLEm*LE2m(n) + c2PttttLPm*LPm+ c2PttttLfE*LfE2(n) + c2PttttLfP*LfP2(n,jv) + c2PttttfEt*fEt2(n) + c2PttttfEtt*fEtt2(n) + c2PttttfPt*fPt2(n,jv)+ c2PttttfPtt*fPtt2(n,jv)
                          fPtttt2(n) = fPtttt2(n) + ptttta
                          if( .false. .and. twilightZone.eq.1 )then
                            write(*,'("")')
                            write(*,'("DI4:RIGHT: j1,j2=",2i3," jv=",i2," n=",i2," ptta,ptte=",2e12.4)') j1,j2,jv,n,ptta,pevtt2(n,jv)
                            write(*,'("        : pttxa,pttxe=",2(1pe12.4)," pttya,pttye=",2(1pe12.4))') pttxa,pevttx2(n,jv),pttya,pevtty2(n,jv)
                            write(*,'("        : ptttta,ptttte=",2(1pe12.4),/)') ptttta,pevtttt2(n,jv)
                            write(*,'(" c2PttLE,c2PttE,c2PttEm,c2PttP,c2PttPm,c2PttfE,c2PttfP=",7(1pe12.4))') c2PtLE,c2PtE,c2PtEm,c2PtP,c2PtPm,c2PtfE,c2PtfP
                            write(*,'(" LEx2,evx2,evnx2,pvx,pvnx,fevx2,fpvx2=",7(1pe12.4))') LEx2(n),evx2(n),evnx2(n),pvx,pvnx,fevx2(n),fpvx2(n,jv)
                            write(*,'(" LEy2,evy2,evny2,pvy,pvny,fevy2,fpvy2=",7(1pe12.4))') LEy2(n),evy2(n),evny2(n),pvy,pvny,fevy2(n),fpvy2(n,jv)
                            ! write(*,'(" LE2,LLE2,LE2m,LP,LPm=",5e12.4)') LE2(n),LLE2(n),LE2m(n),LP,LPm
                            ! write(*,'(" LfE2,LfP2,fEt2,fEtt2,fPt2,fPtt2=",6e12.4)') LfE2(n),LfP2(n,jv),fEt2(n),fEtt2(n),fPt2(n,jv),fPtt2(n,jv)
                            ! write(*,'(" ev,evn,pv,pvn,fev2,fpv2=",6e12.4)')ev,evn,pv,pvn,fev2(n),fpv2(n,jv)
                          end if
                        end do ! end do n 
                      end do ! enddo jv
                    !-----------------------------
                    ! no dispersion
                    !-----------------------------
                    elseif( dispersionModel2.eq.noDispersion .and. nonlinearModel2.eq.noNonlinearModel) then
                      ! do nothing, dispersive forcing is 0
                    end if
                  ! first evaluate the equations we want to solve with the wrong values at the ghost points: (assigns f(0:7))
                    f(0)=(u1x+v1y) - (u2x+v2y)
                    ! [ n.Delta( E )/mu ]=0
                    f(1)=(an1*u1Lap+an2*v1Lap)/mu1 - (an1*u2Lap+an2*v2Lap)/mu2
                    ! [ nv X curl( E) /mu ] = 0 
                    f(2)=(v1x-u1y)/mu1 - (v2x-u2y)/mu2
                   ! f(3)=(tau1*u1Lap+tau2*v1Lap)/eps1 - !      (tau1*u2Lap+tau2*v2Lap)/eps2
                    f(3)=( ( tau1*u1Lap +tau2*v1Lap )/epsmu1 - alphaP1*(tau1*fp1(0)+tau2*fp1(1)) ) - ( ( tau1*u2Lap +tau2*v2Lap )/epsmu2 - alphaP2*(tau1*fp2(0)+tau2*fp2(1)) )
                    ! [ Delta( div(E) )*c^2 ] = 0 
                    ! f(4)=(u1xxx+u1xyy+v1xxy+v1yyy)*c1**2/mu1- !      (u2xxx+u2xyy+v2xxy+v2yyy)*c2**2/mu2
                    f(4)=((u1xxx+u1xyy+v1xxy+v1yyy)/epsmu1-alphaP1*(fPttx1(0) + fPtty1(1))) - ((u2xxx+u2xyy+v2xxy+v2yyy)/epsmu2-alphaP2*(fPttx2(0) + fPtty2(1)))
                    ! print *, 'Divergence of Ptt 0 ??', (fPttx1(0) + fPtty1(1)), (fPttx2(0) + fPtty2(1))
                    ! Dispersive:
                    !   [ (c^2/mu)*{(Delta v).x - (Delta u).y} - (alphaP/mu)*( Py.ttx - Px.tty) ] =0 
                    !    fPttx = P.ttx 
                    !    fPtty = P.tty 
                    f(5)= ( ((v1xxx+v1xyy)-(u1xxy+u1yyy))/epsmu1 - alphaP1*(fPttx1(1) - fPtty1(0)) )/mu1 -( ((v2xxx+v2xyy)-(u2xxy+u2yyy))/epsmu2 - alphaP2*(fPttx2(1) - fPtty2(0)) )/mu2
                    ! f(5)= ( ((v1xxx+2.*v1xyy)-(u1yyy))/epsmu1 - alphaP1*(fPttx1(1) - fPtty1(0)) )/mu1 !      -( ((v2xxx+2.*v2xyy)-(u2yyy))/epsmu2 - alphaP2*(fPttx2(1) - fPtty2(0)) )/mu2
                    ! if( nonlinearModel.ne.noNonlinearModel )then ! MLA
                    !   !  [ nv.( c^2*Delta^2(E) - alphaP*Delta(Ptt) )/mu ] = 0 
                    !   ! Note: fLptt = Delta( Ptt ) 
                    !   f(6)= ( (an1*u1LapSq+an2*v1LapSq)/epsmu1  -alphaP1*( an1*fLPtt1(0)+an2*fLPtt1(1) )  )/mu1 !        -( (an1*u2LapSq+an2*v2LapSq)/epsmu2  -alphaP2*( an1*fLPtt2(0)+an2*fLPtt2(1) )  )/mu2 
                    ! else ! GDM
                      !  [ nv.( c^2*Delta^2(E) - alphaP*Delta(Ptt) )/mu ] = 0 
                      ! Note: fLptt = c^2*Delta( Ptt ) 
                      f(6)= ( (an1*u1LapSq+an2*v1LapSq)*c1**2  -alphaP1*( an1*fLPtt1(0)+an2*fLPtt1(1) )*epsmu1  )/mu1 -( (an1*u2LapSq+an2*v2LapSq)*c2**2  -alphaP2*( an1*fLPtt2(0)+an2*fLPtt2(1) )*epsmu2  )/mu2
                    ! end if
                    ! [ tv.( c^4*Delta^2(E) - alphaP*c^2*Delta(P.tt) - alphaP*P.tttt) ]=0 
                    ! if( nonlinearModel.ne.noNonlinearModel )then ! MLA
                    !   f(7)=(tau1*u1LapSq+tau2*v1LapSq)/(epsmu1**2) - !        (tau1*u2LapSq+tau2*v2LapSq)/(epsmu2**2) !        -alphaP1*( ( tau1*fLPtt1(0) + tau2*fLPtt1(1) )/epsmu1 + tau1*fPtttt1(0) + tau2*fPtttt1(1) ) !        +alphaP2*( ( tau1*fLPtt2(0) + tau2*fLPtt2(1) )/epsmu2 + tau1*fPtttt2(0) + tau2*fPtttt2(1) )
                    ! else ! GDM
                      f(7)=(tau1*u1LapSq+tau2*v1LapSq)/epsmu1**2 - (tau1*u2LapSq+tau2*v2LapSq)/epsmu2**2 -alphaP1*( ( tau1*fLPtt1(0) + tau2*fLPtt1(1) ) + tau1*fPtttt1(0) + tau2*fPtttt1(1) ) +alphaP2*( ( tau1*fLPtt2(0) + tau2*fLPtt2(1) ) + tau1*fPtttt2(0) + tau2*fPtttt2(1) )
                    ! endif
                    ! print *, 'Inside jumps order 4:', alphaP1, alphaP2,fp1(0),fp1(1),fp2(0),fp2(1)
                    if( twilightZone.eq.1 )then
                      call ogderiv(ep, 0,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, uexx )
                      call ogderiv(ep, 0,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, ueyy )
                      call ogderiv(ep, 0,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, vexx )
                      call ogderiv(ep, 0,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, veyy )
                      call ogderiv(ep, 0,3,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, uexxx ) 
                      call ogderiv(ep, 0,2,1,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, uexxy )
                      call ogderiv(ep, 0,1,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, uexyy )
                      call ogderiv(ep, 0,0,3,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, ueyyy )
                      call ogderiv(ep, 0,3,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, vexxx )
                      call ogderiv(ep, 0,1,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, vexyy )
                      call ogderiv(ep, 0,2,1,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, vexxy )
                      call ogderiv(ep, 0,0,3,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, veyyy )
                      call ogderiv(ep, 0,4,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, uexxxx )
                      call ogderiv(ep, 0,2,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, uexxyy )
                      call ogderiv(ep, 0,0,4,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, ueyyyy )
                      call ogderiv(ep, 0,4,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, vexxxx )
                      call ogderiv(ep, 0,2,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, vexxyy )
                      call ogderiv(ep, 0,0,4,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, veyyyy )
                      ueLap = uexx + ueyy
                      veLap = vexx + veyy
                      ueLapSq = uexxxx + 2.*uexxyy + ueyyyy
                      veLapSq = vexxxx + 2.*vexxyy + veyyyy
                      ! f(1) = f(1) - ( an1*ueLap + an2*veLap )*(1./mu1 - 1./mu2)
                      f(1) = f(1) - ( ueLap )*(1./mu1 - 1./mu2)
                      f(2) = f(2) - ( vex - uey  )*(1./mu1 - 1./mu2)
                      ! print *, 'in eval tz: ', vex,uey, mu1,mu2
                      f(3) = f(3) - ( tau1*ueLap +tau2*veLap )*(1./epsmu1-1./epsmu2) + alphaP1*(tau1*pevttSum1(0)+tau2*pevttSum1(1)) - alphaP2*(tau1*pevttSum2(0)+tau2*pevttSum2(1))
                                  ! print *,'check dispersive forcing'
                                  ! print *,pevttSum1(0)-fp1(0),pevttSum1(1)-fp1(1),pevttSum2(0)-fp2(0),pevttSum2(1)-fp2(1)
                      f(4) = f(4) - ((uexxx+uexyy+vexxy+veyyy)/epsmu1-alphaP1*(pevttxSum1(0) + pevttySum1(1))) + ((uexxx+uexyy+vexxy+veyyy)/epsmu2-alphaP2*(pevttxSum2(0) + pevttySum2(1)))
                      ! f(4) = f(4) - ( uexxx+uexyy+vexxy+veyyy )*(c1**2/mu1 - c2**2/mu2)
                      f(5) = f(5) - ((vexxx+vexyy)-(uexxy+ueyyy))*(1./(epsmu1*mu1)-1./(epsmu2*mu2)) + (alphaP1/mu1)*( pevttxSum1(1) - pevttySum1(0) ) - (alphaP2/mu2)*( pevttxSum2(1) - pevttySum2(0) )
                     ! f(5) = f(5) - ((vexxx+2.*vexyy)-(ueyyy))*(1./(epsmu1*mu1)-1./(epsmu2*mu2)) !             + (alphaP1/mu1)*( pevttxSum1(1) - pevttySum1(0) ) !             - (alphaP2/mu2)*( pevttxSum2(1) - pevttySum2(0) )
                                 ! print *, fPttx1(1)-pevttxSum1(1),fPtty1(0)-pevttySum1(0),fPttx2(1)-pevttxSum2(1),fPtty2(0)-pevttySum2(0)
                      ! if( nonlinearModel.ne.noNonlinearModel )then ! MLA
                      !   f(6) = f(6) - (an1*ueLapSq+an2*veLapSq)*(1./(epsmu1*mu1) - 1./(epsmu2*mu2) ) !               + (alphaP1/mu1)*( an1*pevttLSum1(0) + an2*pevttLSum1(1) ) !               - (alphaP2/mu2)*( an1*pevttLSum2(0) + an2*pevttLSum2(1) ) 
                      !               ! print *, fLPtt1(0)-pevttLSum1(0),fLPtt1(1)-pevttLSum1(1),fLPtt2(0)-pevttLSum2(0),fLPtt2(1)-pevttLSum2(1)
                      ! else ! GDM
                      f(6) = f(6) - (an1*ueLapSq+an2*veLapSq)*( c1**2/mu1 - c2**2/mu2 ) + (alphaP1/mu1)*( an1*pevttLSum1(0) + an2*pevttLSum1(1) )*epsmu1 - (alphaP2/mu2)*( an1*pevttLSum2(0) + an2*pevttLSum2(1) )*epsmu2
                      ! end if
                      ! if( nonlinearModel.ne.noNonlinearModel )then ! MLA
                      !   f(7) = f(7) - (tau1*ueLapSq+tau2*veLapSq)*(1./(epsmu1**2) - 1./(epsmu2**2)) !               + alphaP1*( ( tau1*pevttLSum1(0) + tau2*pevttLSum1(1) )/epsmu1 + tau1*pevttttSum1(0) + tau2*pevttttSum1(1) ) !               - alphaP2*( ( tau1*pevttLSum2(0) + tau2*pevttLSum2(1) )/epsmu2 + tau1*pevttttSum2(0) + tau2*pevttttSum2(1) ) 
                      !             ! print *, fPtttt1(0)-pevttttSum1(0),fPtttt1(1)-pevttttSum1(1),fPtttt2(0)-pevttttSum2(0),fPtttt2(1)-pevttttSum2(1)
                      ! else ! GDM
                      f(7) = f(7) - (tau1*ueLapSq+tau2*veLapSq)*( c1**4 - c2**4 ) + alphaP1*( ( tau1*pevttLSum1(0) + tau2*pevttLSum1(1) ) + tau1*pevttttSum1(0) + tau2*pevttttSum1(1) ) - alphaP2*( ( tau1*pevttLSum2(0) + tau2*pevttLSum2(1) ) + tau1*pevttttSum2(0) + tau2*pevttttSum2(1) )
                      ! endif
                   end if
                  if( debug.gt.7 ) write(debugFile,'(" --> 4cth: j1,j2=",2i4," u1xx,u1yy,u2xx,u2yy=",4e10.2)') j1,j2,u1xx,u1yy,u2xx,u2yy
                    ! '
                   if( debug.gt.3 ) write(debugFile,'(" --> 4cth: i1,i2=",2i4," f(start)=",8e10.2)') i1,i2,f(0),f(1),f(2),f(3),f(4),f(5),f(6),f(7)
                 ! here is the matrix of coefficients for the unknowns u1(-1),v1(-1),u2(-1),v2(-1)
                 ! Solve:
                 !     
                 !       A [ U ] = A [ U(old) ] - [ f ]
                 ! Equation 0: 
                 ! 0  [ u.x + v.y ] = 0
                 a8(0,0) = -is*8.*rx1*dx141(axis1)     ! coeff of u1(-1) from [u.x+v.y] 
                 a8(0,1) = -is*8.*ry1*dx141(axis1)     ! coeff of v1(-1) from [u.x+v.y] 
                 a8(0,4) =  is*rx1*dx141(axis1)        ! u1(-2)
                 a8(0,5) =  is*ry1*dx141(axis1)        ! v1(-2) 
                 a8(0,2) =  js*8.*rx2*dx241(axis2)     ! coeff of u2(-1) from [u.x+v.y] 
                 a8(0,3) =  js*8.*ry2*dx241(axis2) 
                 a8(0,6) = -js*   rx2*dx241(axis2) 
                 a8(0,7) = -js*   ry2*dx241(axis2) 
                 ! 1  [ (u.xx + u.yy)/mu ] = 0
                 a8(1,0) = 16.*dx142(axis1)/mu1         ! coeff of u1(-1) from [u.xx + u.yy]
                 a8(1,1) = 0. 
                 a8(1,4) =    -dx142(axis1)/mu1         ! coeff of u1(-2) from [u.xx + u.yy]
                 a8(1,5) = 0. 
                 a8(1,2) =-16.*dx242(axis2)/mu2         ! coeff of u2(-1) from [u.xx + u.yy]
                 a8(1,3) = 0. 
                 a8(1,6) =     dx242(axis2)/mu2         ! coeff of u2(-2) from [u.xx + u.yy]
                 a8(1,7) = 0. 
                 ! Equation 2: 
                 curl1um1 =  is*8.*ry1*dx141(axis1)   ! coeff of u(-1) from v.x - u.y 
                 curl1vm1 = -is*8.*rx1*dx141(axis1)   ! coeff of v(-1) from v.x - u.y 
                 curl1um2 = -is*   ry1*dx141(axis1)   ! coeff of u(-2) from v.x - u.y 
                 curl1vm2 =  is*   rx1*dx141(axis1)   ! coeff of v(-2) from v.x - u.y
                 curl2um1 =  js*8.*ry2*dx241(axis2)  ! coeff of u(-1) from v.x - u.y 
                 curl2vm1 = -js*8.*rx2*dx241(axis2)  ! coeff of v(-1) from v.x - u.y 
                 curl2um2 = -js*   ry2*dx241(axis2)  ! coeff of u(-2) from v.x - u.y 
                 curl2vm2 =  js*   rx2*dx241(axis2)  ! coeff of v(-2) from v.x - u.y
                 ! 2  [ (v.x - u.y)/mu ] =0 
                 a8(2,0) =  curl1um1/mu1 
                 a8(2,1) =  curl1vm1/mu1 
                 a8(2,4) =  curl1um2/mu1 
                 a8(2,5) =  curl1vm2/mu1 
                 a8(2,2) = -curl2um1/mu2
                 a8(2,3) = -curl2vm1/mu2
                 a8(2,6) = -curl2um2/mu2
                 a8(2,7) = -curl2vm2/mu2
                 !- a8(2,0) =  is*8.*ry1*dx141(axis1)/mu1 
                 !- a8(2,1) = -is*8.*rx1*dx141(axis1)/mu1     ! coeff of v1(-1) from [v.x - u.y] 
                 !- a8(2,4) = -is*   ry1*dx141(axis1)/mu1 
                 !- a8(2,5) =  is*   rx1*dx141(axis1)/mu1 
                 !- a8(2,2) = -js*8.*ry2*dx241(axis2)/mu2
                 !- a8(2,3) =  js*8.*rx2*dx241(axis2)/mu2
                 !- a8(2,6) =  js*   ry2*dx241(axis2)/mu2
                 !- a8(2,7) = -js*   rx2*dx241(axis2)/mu2
                 aLap0=16.*dx142(axis1)  ! coeff of w1(-1) 
                 aLap1=-1.*dx142(axis1)  ! coeff of w1(-2)
                 bLap0=16.*dx242(axis2)  ! coeff of w2(-1) 
                 bLap1=-1.*dx242(axis2)  ! coeff of w2(-1) 
                 aLapSq0= ( -4./(dx1(axis1)**4) )
                 aLapSq1= (  1./(dx1(axis1)**4) )
                 bLapSq0= ( -4./(dx2(axis2)**4) )
                 bLapSq1= (  1./(dx2(axis2)**4) )
                 ! -------------- Equation 3 -----------------------
                 !   [ tau.{ (uv.xx+uv.yy)/eps -alphaP*P.tt } ] = 0
                 !    P.tt = c4PttLEsum * L(E) + c4PttLLEsum* L^2(E) + ...
                 if( nonlinearModel1.ne.noNonlinearModel )then
                   ! coeff of P is not used, thus set to 0
                   c4PttLEsum1 = 0.
                   c4PttLLEsum1 = 0.
                 endif
                 if( nonlinearModel2.ne.noNonlinearModel )then
                   c4PttLEsum2 = 0.
                   c4PttLLEsum2 = 0.
                 endif
                 a8(3,0) = tau1*( aLap0*( 1./epsmu1 -alphaP1*c4PttLEsum1/epsmu1 ) - aLapSq0*alphaP1*c4PttLLEsum1/epsmu1**2 )
                 a8(3,1) = tau2*( aLap0*( 1./epsmu1 -alphaP1*c4PttLEsum1/epsmu1 ) - aLapSq0*alphaP1*c4PttLLEsum1/epsmu1**2 )
                 a8(3,4) = tau1*( aLap1*( 1./epsmu1 -alphaP1*c4PttLEsum1/epsmu1 ) - aLapSq1*alphaP1*c4PttLLEsum1/epsmu1**2 )
                 a8(3,5) = tau2*( aLap1*( 1./epsmu1 -alphaP1*c4PttLEsum1/epsmu1 ) - aLapSq1*alphaP1*c4PttLLEsum1/epsmu1**2 )
                 a8(3,2) =-tau1*( bLap0*( 1./epsmu2 -alphaP2*c4PttLEsum2/epsmu2 ) - bLapSq0*alphaP2*c4PttLLEsum2/epsmu2**2 )
                 a8(3,3) =-tau2*( bLap0*( 1./epsmu2 -alphaP2*c4PttLEsum2/epsmu2 ) - bLapSq0*alphaP2*c4PttLLEsum2/epsmu2**2 )
                 a8(3,6) =-tau1*( bLap1*( 1./epsmu2 -alphaP2*c4PttLEsum2/epsmu2 ) - bLapSq1*alphaP2*c4PttLLEsum2/epsmu2**2 )
                 a8(3,7) =-tau2*( bLap1*( 1./epsmu2 -alphaP2*c4PttLEsum2/epsmu2 ) - bLapSq1*alphaP2*c4PttLLEsum2/epsmu2**2 )
                 ! -------------- Equation 4 -----------------------
                 !    [ (u.xx+u.yy).x + (v.xx+v.yy).y ] = 0
                 dxxx1by2i = 1./(dx1(axis1)**3*2); 
                 dxxx2by2i = 1./(dx2(axis2)**3*2); 
                 ! Order u1(-1), v1(-1), u2(-1), v2(-1), u1(-2), v1(-2), u2(-2), v2(-2), 
                 a8(4,0)= ( is*rx1 *2.*dxxx1by2i )*c1**2
                 a8(4,1)= 0.
                 a8(4,2)=-( js*rx2* 2.*dxxx2by2i )*c2**2
                 a8(4,3)= 0.
                 a8(4,4)= (-is*rx1    *dxxx1by2i )*c1**2
                 a8(4,5)= 0.
                 a8(4,6)=-(-js*rx2    *dxxx2by2i )*c2**2  
                 a8(4,7)= 0.
                 ! ---------------- Equation 5 (2nd-order) -----------------
                 !   [ ( {(Delta v).x - (Delta u).y}/(epsmu) - alphaP*( Py.ttx - Px.tty) )/mu ] =0 
                 !
                 ! check me: 
                 aLapX0 =  is*rx1*2.*dx122(axis1)*dx112(axis1) 
                 bLapY0 = -is*ry1*2.*dx122(axis1)*dx112(axis1)
                 aLapX1 = -is*rx1   *dx122(axis1)*dx112(axis1)
                 bLapY1 =  is*ry1   *dx122(axis1)*dx112(axis1)
                 cLapX0 =  js*rx2*2.*dx222(axis2)*dx212(axis2) 
                 dLapY0 = -js*ry2*2.*dx222(axis2)*dx212(axis2)
                 cLapX1 = -js*rx2   *dx222(axis2)*dx212(axis2)
                 dLapY1 =  js*ry2   *dx222(axis2)*dx212(axis2)
                 !     P.tt = c2PttLEsum * L(E)
                  if( nonlinearModel1.ne.noNonlinearModel )then
                    ! coeff of P is set to 0
                    c2PttLEsum1 = 0.
                    c2PttEsum1 = 0.
                  endif
                  if( nonlinearModel2.ne.noNonlinearModel )then
                    c2PttLEsum2 = 0.
                    c2PttEsum2 = 0.
                  endif
                 eqnCoeff = ( 1./epsmu1 - alphaP1*c2PttLEsum1/epsmu1 )/mu1 
                 eqnCoeffb = -alphaP1*c2PttEsum1/mu1 ! added sept 16, 2018 
                 a8(5,0)=-bLapY0*eqnCoeff + curl1um1*eqnCoeffb  
                 a8(5,1)= aLapX0*eqnCoeff + curl1vm1*eqnCoeffb
                 a8(5,4)=-bLapY1*eqnCoeff + curl1um2*eqnCoeffb 
                 a8(5,5)= aLapX1*eqnCoeff + curl1vm2*eqnCoeffb
                 eqnCoeff = ( 1./epsmu2 - alphaP2*c2PttLEsum2/epsmu2 )/mu2 
                 eqnCoeffb = -alphaP2*c2PttEsum2/mu2 ! added sept 16, 2018 
                 a8(5,2)=-(-dLapY0*eqnCoeff + curl2um1*eqnCoeffb)
                 a8(5,3)=-( cLapX0*eqnCoeff + curl2vm1*eqnCoeffb)
                 a8(5,6)=-(-dLapY1*eqnCoeff + curl2um2*eqnCoeffb)
                 a8(5,7)=-( cLapX1*eqnCoeff + curl2vm2*eqnCoeffb)
                 ! ------- Equation 6 -----
                 !  [ nv.( c^2*Delta^2(E) - alphaP*Delta(Ptt) )/mu ] = 0 
                 ! 6  [ n.Delta^2 u/eps ] = 0
                 ! use Eqn 6 
                 ! NOTE: LE = c^2*Delta(E) and LLE = (c^4*Delta^2) E 
                 ! Note: the coeff of L(E) in Delta(Ptt) is the coeff of E in Ptt
                 ! Note: the coeff of LL(E) in Delta(Ptt) is the coeff of LE in Ptt
                  if( nonlinearModel1.ne.noNonlinearModel )then
                    c2PttEsum1 = 0.
                    c2PttLEsum1 = 0.
                  endif
                  if( nonlinearModel2.ne.noNonlinearModel )then
                    c2PttEsum2 = 0.
                    c2PttLEsum2 = 0.
                  endif
                 a8(6,0) = an1*( aLapSq0/epsmu1 -alphaP1*( c2PttEsum1*aLap0 + c2PttLEsum1*aLapSq0/epsmu1 ) )/mu1
                 a8(6,1) = an2*( aLapSq0/epsmu1 -alphaP1*( c2PttEsum1*aLap0 + c2PttLEsum1*aLapSq0/epsmu1 ) )/mu1
                 a8(6,4) = an1*( aLapSq1/epsmu1 -alphaP1*( c2PttEsum1*aLap1 + c2PttLEsum1*aLapSq1/epsmu1 ) )/mu1
                 a8(6,5) = an2*( aLapSq1/epsmu1 -alphaP1*( c2PttEsum1*aLap1 + c2PttLEsum1*aLapSq1/epsmu1 ) )/mu1
                 a8(6,2) =-an1*( bLapSq0/epsmu2 -alphaP2*( c2PttEsum2*bLap0 + c2PttLEsum2*bLapSq0/epsmu2 ) )/mu2
                 a8(6,3) =-an2*( bLapSq0/epsmu2 -alphaP2*( c2PttEsum2*bLap0 + c2PttLEsum2*bLapSq0/epsmu2 ) )/mu2
                 a8(6,6) =-an1*( bLapSq1/epsmu2 -alphaP2*( c2PttEsum2*bLap1 + c2PttLEsum2*bLapSq1/epsmu2 ) )/mu2
                 a8(6,7) =-an2*( bLapSq1/epsmu2 -alphaP2*( c2PttEsum2*bLap1 + c2PttLEsum2*bLapSq1/epsmu2 ) )/mu2
                 ! ------- Equation 7 ------
                 ! [ tv.( c^4*Delta^2(E) - alphaP*c^2*Delta(P.tt) - alphaP*P.tttt) ]=0 
                 ! 7  [ tau.Delta^2 v/eps^2 ] = 0 
                 ! Note: the coeff of L(E) in Delta(Ptt) is the coeff of E in Ptt
                 ! Note: the coeff of LL(E) in Delta(Ptt) is the coeff of LE in Ptt
                 if( nonlinearModel1.ne.noNonlinearModel )then
                   c2PttEsum1 = 0.
                   c2PttttLEsum1 = 0.
                   c2PttLEsum1 = 0.
                   c2PttttLLEsum1 = 0.
                 endif
                 if( nonlinearModel2.ne.noNonlinearModel )then
                   c2PttEsum2 = 0.
                   c2PttttLEsum2 = 0.
                   c2PttLEsum2 = 0.
                   c2PttttLLEsum2 = 0.
                 endif
                 coeffLap1   =              -alphaP1*(  c2PttEsum1 + c2PttttLEsum1  )/epsmu1
                 coeffLapSq1 = 1./epsmu1**2 -alphaP1*( c2PttLEsum1 + c2PttttLLEsum1 )/epsmu1**2
                 coeffLap2   =              -alphaP2*(  c2PttEsum2 + c2PttttLEsum2  )/epsmu2
                 coeffLapSq2 = 1./epsmu2**2 -alphaP2*( c2PttLEsum2 + c2PttttLLEsum2 )/epsmu2**2
                 a8(7,0) = tau1*( coeffLapSq1*aLapSq0 + coeffLap1*aLap0 )
                 a8(7,1) = tau2*( coeffLapSq1*aLapSq0 + coeffLap1*aLap0 )
                 a8(7,4) = tau1*( coeffLapSq1*aLapSq1 + coeffLap1*aLap1 )
                 a8(7,5) = tau2*( coeffLapSq1*aLapSq1 + coeffLap1*aLap1 )
                 a8(7,2) =-tau1*( coeffLapSq2*bLapSq0 + coeffLap2*bLap0 )
                 a8(7,3) =-tau2*( coeffLapSq2*bLapSq0 + coeffLap2*bLap0 )
                 a8(7,6) =-tau1*( coeffLapSq2*bLapSq1 + coeffLap2*bLap1 )
                 a8(7,7) =-tau2*( coeffLapSq2*bLapSq1 + coeffLap2*bLap1 )
                 ! if( debug>7 .and. i2.le.0 )then
                 !   write(*,*) "unified24r: Matrix a8"
                 !   do n=0,7
                 !     write(*,'(8(1pe10.2))') (a8(n,nn),nn=0,7)
                 !   end do 
                 ! end if 
                 ! --- check matrix coefficients by delta function approach ----
                 if( checkCoeff.eq.1 )then
                  !! numberOfEquations=8
                  !! checkCoefficients(i1,i2,i3, j1,j2,j3,numberOfEquations,a8,evalDispersiveInterfaceEquations24r )
                 end if
                 ! -- save current (wrong) ghost values in q()
                 q(0) = u1(i1-is1,i2-is2,i3,ex)
                 q(1) = u1(i1-is1,i2-is2,i3,ey)
                 q(2) = u2(j1-js1,j2-js2,j3,ex)
                 q(3) = u2(j1-js1,j2-js2,j3,ey)
                 q(4) = u1(i1-2*is1,i2-2*is2,i3,ex)
                 q(5) = u1(i1-2*is1,i2-2*is2,i3,ey)
                 q(6) = u2(j1-2*js1,j2-2*js2,j3,ex)
                 q(7) = u2(j1-2*js1,j2-2*js2,j3,ey)
                 !- if( debug.gt.4 )then
                 !-   write(debugFile,'(" --> unified24r: i1,i2=",2i4," q=",8e10.2)') i1,i2,(q(n),n=0,7)
                 !- end if
                 if( debug.gt.7 )then
                   write(*,'(" --> before unified24r: i1,i2=",2i4," RHS(A)=",8(1pe13.5))') i1,i2,(f(n),n=0,7)
                 end if
                 ! subtract off the contributions from the initial (wrong) values at the ghost points:
                 do n=0,7
                   f(n) = (a8(n,0)*q(0)+a8(n,1)*q(1)+a8(n,2)*q(2)+a8(n,3)*q(3)+a8(n,4)*q(4)+a8(n,5)*q(5)+a8(n,6)*q(6)+a8(n,7)*q(7)) - f(n)
                 end do
                 if( debug.gt.7 )then
                   write(*,'(" --> after unified24r: i1,i2=",2i4," RHS(A)=",8(1pe13.5))') i1,i2,(f(n),n=0,7)
                 end if
                 ! solve A Q = F
                 ! factor the matrix
                 numberOfEquations=8
                 call dgeco( a8(0,0), numberOfEquations, numberOfEquations, ipivot8(0),rcond,work(0))
                 if( debug.gt.3 ) then
                   write(debugFile,'(" --> unified24r: i1,i2=",2i4," rcond=",e10.2)') i1,i2,rcond
                 end if
                 ! solve A Q = F
                 job=0
                 numberOfEquations=8
                 call dgesl( a8(0,0), numberOfEquations, numberOfEquations, ipivot8(0), f(0), job)
                 if( useJacobiUpdate.eq.0 )then
                   u1(i1-is1,i2-is2,i3,ex)=f(0)
                   u1(i1-is1,i2-is2,i3,ey)=f(1)
                   u2(j1-js1,j2-js2,j3,ex)=f(2)
                   u2(j1-js1,j2-js2,j3,ey)=f(3)
                   u1(i1-2*is1,i2-2*is2,i3,ex)=f(4)
                   u1(i1-2*is1,i2-2*is2,i3,ey)=f(5)
                   u2(j1-2*js1,j2-2*js2,j3,ex)=f(6)
                   u2(j1-2*js1,j2-2*js2,j3,ey)=f(7)
                 else
                   ! Jacobi-update
                   wk1(i1-is1,i2-is2,i3,ex)    =f(0)
                   wk1(i1-is1,i2-is2,i3,ey)    =f(1)
                   wk2(j1-js1,j2-js2,j3,ex)    =f(2)
                   wk2(j1-js1,j2-js2,j3,ey)    =f(3)
                   wk1(i1-2*is1,i2-2*is2,i3,ex)=f(4)
                   wk1(i1-2*is1,i2-2*is2,i3,ey)=f(5)
                   wk2(j1-2*js1,j2-2*js2,j3,ex)=f(6)
                   wk2(j1-2*js1,j2-2*js2,j3,ey)=f(7)
                 end if
                 ! if( debug>3 .and. twilightZone.eq.1 )then
                 !   ! check errors
                 !   call ogderiv(ep, 0,0,0,0, xy1(i1-is1,i2-is2,i3,0),xy1(i1-is1,i2-is2,i3,1),0.,t, ex, evv(0) )
                 !   call ogderiv(ep, 0,0,0,0, xy1(i1-is1,i2-is2,i3,0),xy1(i1-is1,i2-is2,i3,1),0.,t, ey, evv(1) )
                 !   call ogderiv(ep, 0,0,0,0, xy2(j1-js1,j2-js2,j3,0),xy2(j1-js1,j2-js2,j3,1),0.,t, ex, evv(2) )
                 !   call ogderiv(ep, 0,0,0,0, xy2(j1-js1,j2-js2,j3,0),xy2(j1-js1,j2-js2,j3,1),0.,t, ey, evv(3) )
                 !   call ogderiv(ep, 0,0,0,0, xy1(i1-2*is1,i2-2*is2,i3,0),xy1(i1-2*is1,i2-2*is2,i3,1),0.,t, ex, evv(4) )
                 !   call ogderiv(ep, 0,0,0,0, xy1(i1-2*is1,i2-2*is2,i3,0),xy1(i1-2*is1,i2-2*is2,i3,1),0.,t, ey, evv(5) )
                 !   call ogderiv(ep, 0,0,0,0, xy2(j1-2*js1,j2-2*js2,j3,0),xy2(j1-2*js1,j2-2*js2,j3,1),0.,t, ex, evv(6) )
                 !   call ogderiv(ep, 0,0,0,0, xy2(j1-2*js1,j2-2*js2,j3,0),xy2(j1-2*js1,j2-2*js2,j3,1),0.,t, ey, evv(7) )
                 !   write(*,'("gdm24r: Errors in u1(-1), v1(-1), u2(-1), v2(-1), u1(-2), v1(-2), u2(-2), v2(-2)")') 
                 !   maxErr=0.
                 !   do n=0,7
                 !     maxErr =max(maxErr,abs(evv(n)-f(n)))
                 !   end do
                 !   write(*,'("gdm24r: i1,i2=",2i4," err= ",8e8.1," -> maxErr=",e8.1)') i1,i2, (abs(evv(n)-f(n)),n=0,7),maxErr
                 ! end if
                 !-  if( .false. .and. debug.gt.2 )then 
                 !-    ! --- check residuals in the jump conditions ----
                 !-
                 !-    ! Evaluate the jump conditions using the new values at the ghost points 
                 !-    !! evaluateDispersiveInterfaceEquations2dOrder4()
                 !-    write(debugFile,'(" JUMP-residuals: i1,i2=",2i4," f(re-eval)=",8e10.2)') i1,i2,f(0),f(1),f(2),f(3),f(4),f(5),f(6),f(7)
                 !-  end if
                ! ******************************************************
                ! NOTE: This is just the non-dispersive case ***FIX ME FOR GDM ***
                ! 
                ! solve for Hz
                !  [ w.n/eps ] = 0
                !  [ lap(w)/eps ] = 0
                !  [ lap(w).n/eps**2 ] = 0
                !  [ lapSq(w)/eps**2 ] = 0
                 ! first evaluate the equations we want to solve with the wrong values at the ghost points:
                  ! These derivatives are computed to 4th-order accuracy
                    ww1=u1(i1,i2,i3,hz) ! in the rectangular case just eval the solution
                     w1x = (u1(i1-2,i2,i3,hz)-8.*u1(i1-1,i2,i3,hz)+8.*u1(i1+1,i2,i3,hz)-u1(i1+2,i2,i3,hz))/(12.*dx1(0))
                     w1y = (u1(i1,i2-2,i3,hz)-8.*u1(i1,i2-1,i3,hz)+8.*u1(i1,i2+1,i3,hz)-u1(i1,i2+2,i3,hz))/(12.*dx1(1))
                     w1xx = (-u1(i1-2,i2,i3,hz)+16.*u1(i1-1,i2,i3,hz)-30.*u1(i1,i2,i3,hz)+16.*u1(i1+1,i2,i3,hz)-u1(i1+2,i2,i3,hz))/(12.*dx1(0)**2)
                     w1yy = (-u1(i1,i2-2,i3,hz)+16.*u1(i1,i2-1,i3,hz)-30.*u1(i1,i2,i3,hz)+16.*u1(i1,i2+1,i3,hz)-u1(i1,i2+2,i3,hz))/(12.*dx1(1)**2)
                   w1Lap = w1xx+ w1yy
                    ww2=u2(j1,j2,j3,hz) ! in the rectangular case just eval the solution
                     w2x = (u2(j1-2,j2,j3,hz)-8.*u2(j1-1,j2,j3,hz)+8.*u2(j1+1,j2,j3,hz)-u2(j1+2,j2,j3,hz))/(12.*dx2(0))
                     w2y = (u2(j1,j2-2,j3,hz)-8.*u2(j1,j2-1,j3,hz)+8.*u2(j1,j2+1,j3,hz)-u2(j1,j2+2,j3,hz))/(12.*dx2(1))
                     w2xx = (-u2(j1-2,j2,j3,hz)+16.*u2(j1-1,j2,j3,hz)-30.*u2(j1,j2,j3,hz)+16.*u2(j1+1,j2,j3,hz)-u2(j1+2,j2,j3,hz))/(12.*dx2(0)**2)
                     w2yy = (-u2(j1,j2-2,j3,hz)+16.*u2(j1,j2-1,j3,hz)-30.*u2(j1,j2,j3,hz)+16.*u2(j1,j2+1,j3,hz)-u2(j1,j2+2,j3,hz))/(12.*dx2(1)**2)
                   w2Lap = w2xx+ w2yy
                  ! These derivatives are computed to 2nd-order accuracy
                    ww1=u1(i1,i2,i3,hz) ! in the rectangular case just eval the solution
                     w1xxx = (-u1(i1-2,i2,i3,hz)+2.*u1(i1-1,i2,i3,hz)-2.*u1(i1+1,i2,i3,hz)+u1(i1+2,i2,i3,hz))/(2.*dx1(0)**3)
                     w1xxy = ((-u1(i1-1,i2-1,i3,hz)+u1(i1-1,i2+1,i3,hz))/(2.*dx1(1))-2.*(-u1(i1,i2-1,i3,hz)+u1(i1,i2+1,i3,hz))/(2.*dx1(1))+(-u1(i1+1,i2-1,i3,hz)+u1(i1+1,i2+1,i3,hz))/(2.*dx1(1)))/(dx1(0)**2)
                     w1xyy = (-(u1(i1-1,i2-1,i3,hz)-2.*u1(i1-1,i2,i3,hz)+u1(i1-1,i2+1,i3,hz))/(dx1(1)**2)+(u1(i1+1,i2-1,i3,hz)-2.*u1(i1+1,i2,i3,hz)+u1(i1+1,i2+1,i3,hz))/(dx1(1)**2))/(2.*dx1(0))
                     w1yyy = (-u1(i1,i2-2,i3,hz)+2.*u1(i1,i2-1,i3,hz)-2.*u1(i1,i2+1,i3,hz)+u1(i1,i2+2,i3,hz))/(2.*dx1(1)**3)
                     w1xxxx = (u1(i1-2,i2,i3,hz)-4.*u1(i1-1,i2,i3,hz)+6.*u1(i1,i2,i3,hz)-4.*u1(i1+1,i2,i3,hz)+u1(i1+2,i2,i3,hz))/(dx1(0)**4)
                     w1xxyy = ((u1(i1-1,i2-1,i3,hz)-2.*u1(i1-1,i2,i3,hz)+u1(i1-1,i2+1,i3,hz))/(dx1(1)**2)-2.*(u1(i1,i2-1,i3,hz)-2.*u1(i1,i2,i3,hz)+u1(i1,i2+1,i3,hz))/(dx1(1)**2)+(u1(i1+1,i2-1,i3,hz)-2.*u1(i1+1,i2,i3,hz)+u1(i1+1,i2+1,i3,hz))/(dx1(1)**2))/(dx1(0)**2)
                     w1yyyy = (u1(i1,i2-2,i3,hz)-4.*u1(i1,i2-1,i3,hz)+6.*u1(i1,i2,i3,hz)-4.*u1(i1,i2+1,i3,hz)+u1(i1,i2+2,i3,hz))/(dx1(1)**4)
                   w1LapSq = w1xxxx +2.* w1xxyy + w1yyyy
                    ww2=u2(j1,j2,j3,hz) ! in the rectangular case just eval the solution
                     w2xxx = (-u2(j1-2,j2,j3,hz)+2.*u2(j1-1,j2,j3,hz)-2.*u2(j1+1,j2,j3,hz)+u2(j1+2,j2,j3,hz))/(2.*dx2(0)**3)
                     w2xxy = ((-u2(j1-1,j2-1,j3,hz)+u2(j1-1,j2+1,j3,hz))/(2.*dx2(1))-2.*(-u2(j1,j2-1,j3,hz)+u2(j1,j2+1,j3,hz))/(2.*dx2(1))+(-u2(j1+1,j2-1,j3,hz)+u2(j1+1,j2+1,j3,hz))/(2.*dx2(1)))/(dx2(0)**2)
                     w2xyy = (-(u2(j1-1,j2-1,j3,hz)-2.*u2(j1-1,j2,j3,hz)+u2(j1-1,j2+1,j3,hz))/(dx2(1)**2)+(u2(j1+1,j2-1,j3,hz)-2.*u2(j1+1,j2,j3,hz)+u2(j1+1,j2+1,j3,hz))/(dx2(1)**2))/(2.*dx2(0))
                     w2yyy = (-u2(j1,j2-2,j3,hz)+2.*u2(j1,j2-1,j3,hz)-2.*u2(j1,j2+1,j3,hz)+u2(j1,j2+2,j3,hz))/(2.*dx2(1)**3)
                     w2xxxx = (u2(j1-2,j2,j3,hz)-4.*u2(j1-1,j2,j3,hz)+6.*u2(j1,j2,j3,hz)-4.*u2(j1+1,j2,j3,hz)+u2(j1+2,j2,j3,hz))/(dx2(0)**4)
                     w2xxyy = ((u2(j1-1,j2-1,j3,hz)-2.*u2(j1-1,j2,j3,hz)+u2(j1-1,j2+1,j3,hz))/(dx2(1)**2)-2.*(u2(j1,j2-1,j3,hz)-2.*u2(j1,j2,j3,hz)+u2(j1,j2+1,j3,hz))/(dx2(1)**2)+(u2(j1+1,j2-1,j3,hz)-2.*u2(j1+1,j2,j3,hz)+u2(j1+1,j2+1,j3,hz))/(dx2(1)**2))/(dx2(0)**2)
                     w2yyyy = (u2(j1,j2-2,j3,hz)-4.*u2(j1,j2-1,j3,hz)+6.*u2(j1,j2,j3,hz)-4.*u2(j1,j2+1,j3,hz)+u2(j1,j2+2,j3,hz))/(dx2(1)**4)
                   w2LapSq = w2xxxx +2.* w2xxyy + w2yyyy
                  f(0)=(an1*w1x+an2*w1y)/eps1 - (an1*w2x+an2*w2y)/eps2
                  f(1)=w1Lap/eps1 - w2Lap/eps2
                  f(2)=(an1*(w1xxx+w1xyy)+an2*(w1xxy+w1yyy))/eps1**2 - (an1*(w2xxx+w2xyy)+an2*(w2xxy+w2yyy))/eps2**2
                  f(3)=w1LapSq/eps1**2 - w2LapSq/eps2**2
                  if( twilightZone.eq.1 )then
                    call ogderiv(ep, 0,1,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, hz, wex  )
                    call ogderiv(ep, 0,0,1,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, hz, wey  )
                    call ogderiv(ep, 0,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, hz, wexx )
                    call ogderiv(ep, 0,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, hz, weyy )
                    call ogderiv(ep, 0,3,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, hz, wexxx )
                    call ogderiv(ep, 0,2,1,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, hz, wexxy )
                    call ogderiv(ep, 0,1,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, hz, wexyy )
                    call ogderiv(ep, 0,0,3,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, hz, weyyy )
                    call ogderiv(ep, 0,4,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, hz, wexxxx )
                    call ogderiv(ep, 0,2,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, hz, wexxyy )
                    call ogderiv(ep, 0,0,4,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, hz, weyyyy )
                    weLap = wexx + weyy
                    weLapSq = wexxxx + 2.*wexxyy + weyyyy
                    f(0) = f(0) - (an1*wex+an2*wey)*(1./eps1 - 1./eps2)
                    f(1) = f(1) - ( weLap )*(1./eps1 - 1./eps2)
                    f(2) = f(2) - (an1*(wexxx+wexyy)+an2*(wexxy+weyyy))*(1./eps1**2 - 1./eps2**2)
                    f(3) = f(3) - weLapSq*(1./eps1**2 - 1./eps2**2)
                  end if
                 ! form the matrix for computing Hz
                 ! 1: [ w.n/eps ] = 0
                 a0 = dx141(axis1)/eps1
                 b0 = dx241(axis2)/eps2
                 a4h(0,0) = -is*8.*a0
                 a4h(0,2) =  is*   a0
                 a4h(0,1) =  js*8.*b0
                 a4h(0,3) = -js*   b0
                 ! 2: [ lap(w)/eps ] = 0 
                 aLap0=16.*dx142(axis1)  ! coeff of w1(-1) 
                 aLap1=-1.*dx142(axis1)  ! coeff of w1(-2)
                 bLap0=16.*dx242(axis2)  ! coeff of w2(-1) 
                 bLap1=-1.*dx242(axis2)  ! coeff of w2(-1) 
                 a4h(1,0) = aLap0/eps1    ! coeff of w1(-1) 
                 a4h(1,2) = aLap1/eps1    ! coeff of w1(-2)
                 a4h(1,1) =-bLap0/eps2    ! coeff of w2(-1) 
                 a4h(1,3) =-bLap1/eps2    ! coeff of w2(-1) 
                 ! 3:  [ (an1*(w.xx+w.yy).x + an2.(w.xx+w.yy).y)/eps**2 ] = 0
                 !  a4h(2,0)= (an1*aLapX0+an2*bLapY0)/eps1**2  ! coeff of w1(-1) 
                 !  a4h(2,1)=-(an1*cLapX0+an2*dLapY0)/eps2**2  ! coeff of w2(-1)
                 !  a4h(2,2)= (an1*aLapX1+an2*bLapY1)/eps1**2  ! coeff of w1(-2)
                 !  a4h(2,3)=-(an1*cLapX1+an2*dLapY1)/eps2**2  ! coeff of w2(-2)
                 a4h(2,0)=  is*(  1./(dx1(axis1)**3) +1./(dx1(axis1)*dx1(axis1p1)**2) )/eps1**2  ! coeff of w1(-1) aLapX0
                 a4h(2,2)=  is*( -.5/(dx1(axis1)**3)                                  )/eps1**2  ! coeff of w1(-2) aLapX1
                 a4h(2,1)= -js*(  1./(dx2(axis2)**3) +1./(dx2(axis2)*dx2(axis2p1)**2) )/eps2**2  ! coeff of w2(-1) cLapX0 
                 a4h(2,3)= -js*( -.5/(dx2(axis2)**3)                                  )/eps2**2  ! coeff of w2(-2) cLapX1
                 ! 4 [ lapSq(w)/eps**2 ] = 0   [ w_xxxx + 2 * w_xxyy + w_yyyy ]
                 aLapSq0= ( -4./(dx1(axis1)**4) -4./(dx1(axis1)**2 * dx1(axis1p1)**2 ) )
                 aLapSq1= (  1./(dx1(axis1)**4) )
                 bLapSq0= ( -4./(dx2(axis2)**4) -4./(dx2(axis2)**2 * dx2(axis2p1)**2 ) )
                 bLapSq1= (  1./(dx2(axis2)**4) )
                 a4h(3,0) = aLapSq0/eps1**2
                 a4h(3,2) = aLapSq1/eps1**2
                 a4h(3,1) =-bLapSq0/eps2**2
                 a4h(3,3) =-bLapSq1/eps2**2
                 q(0) = u1(i1-is1,i2-is2,i3,hz)
                 q(1) = u2(j1-js1,j2-js2,j3,hz)
                 q(2) = u1(i1-2*is1,i2-2*is2,i3,hz)
                 q(3) = u2(j1-2*js1,j2-2*js2,j3,hz)
                 ! subtract off the contributions from the wrong values at the ghost points:
                 do n=0,3
                   f(n) = (a4h(n,0)*q(0)+a4h(n,1)*q(1)+a4h(n,2)*q(2)+a4h(n,3)*q(3)) - f(n)
                 end do
                 ! write(*,'(" a4h=",4(e9.2,1x))') ((a4h(i,j),j=0,3),i=0,3)
                 ! factor the matrix
                 ! numberOfEquations=4
                 call dgeco( a4h(0,0), 4, 4, ipivot4h(0),rcond,work(0))
                 ! write(*,'("rcond=",e12.4)') rcond
                 ! solve
                 job=0
                 call dgesl( a4h(0,0), 4, 4, ipivot4h(0), f(0), job)
                 if( useJacobiUpdate.eq.0 )then
                   u1(i1-  is1,i2-  is2,i3,hz)=f(0)
                   u2(j1-  js1,j2-  js2,j3,hz)=f(1)
                   u1(i1-2*is1,i2-2*is2,i3,hz)=f(2)
                   u2(j1-2*js1,j2-2*js2,j3,hz)=f(3)
                 else
                   ! Jacobi update -- save answer in work space
                   wk1(i1-  is1,i2-  is2,i3,hz)=f(0)
                   wk2(j1-  js1,j2-  js2,j3,hz)=f(1)
                   wk1(i1-2*is1,i2-2*is2,i3,hz)=f(2)
                   wk2(j1-2*js1,j2-2*js2,j3,hz)=f(3)
                 end if
                  end if
                  j1=j1+1
                 end do
                 j2=j2+1
                end do
               ! =============== end loops =======================
               if( useJacobiUpdate.ne.0 )then
                 ! Jacobi-update: now fill in values 
                  i3=n3a
                  j3=m3a
                  j2=m2a
                  do i2=n2a,n2b
                   j1=m1a
                   do i1=n1a,n1b
                    if( mask1(i1,i2,i3).gt.0 .and. mask2(j1,j2,j3).gt.0 )then
                   u1(i1-is1,i2-is2,i3,ex)=wk1(i1-is1,i2-is2,i3,ex)
                   u1(i1-is1,i2-is2,i3,ey)=wk1(i1-is1,i2-is2,i3,ey)
                   u2(j1-js1,j2-js2,j3,ex)=wk2(j1-js1,j2-js2,j3,ex)
                   u2(j1-js1,j2-js2,j3,ey)=wk2(j1-js1,j2-js2,j3,ey)
                   u1(i1-2*is1,i2-2*is2,i3,ex)=wk1(i1-2*is1,i2-2*is2,i3,ex)
                   u1(i1-2*is1,i2-2*is2,i3,ey)=wk1(i1-2*is1,i2-2*is2,i3,ey)
                   u2(j1-2*js1,j2-2*js2,j3,ex)=wk2(j1-2*js1,j2-2*js2,j3,ex)
                   u2(j1-2*js1,j2-2*js2,j3,ey)=wk2(j1-2*js1,j2-2*js2,j3,ey)
                   u1(i1-  is1,i2-  is2,i3,hz)=wk1(i1-  is1,i2-  is2,i3,hz)
                   u2(j1-  js1,j2-  js2,j3,hz)=wk2(j1-  js1,j2-  js2,j3,hz)
                   u1(i1-2*is1,i2-2*is2,i3,hz)=wk1(i1-2*is1,i2-2*is2,i3,hz)
                   u2(j1-2*js1,j2-2*js2,j3,hz)=wk2(j1-2*js1,j2-2*js2,j3,hz)
                    end if
                    j1=j1+1
                   end do
                   j2=j2+1
                  end do
               end if
             else if( useNonlinearModel.eq.0 )then
               ! dispersive case
             else
               ! nonlinear dispersive case
                  ! ****************************************************************
                  ! ***********  NONLINEAR, 2D, ORDER=4, RECTANGULAR **************
                  ! ****************************************************************
                  if( t.le.3.*dt .and. debug.gt.0 )then
                    write(*,'("Interface>>>","24r-MLA")')
                  end if
                  print *, '--------------Now using nonlinear routines for RECTANGULAR --------------'
                  if( t.le.5*dt .or. debug.gt.3 )then
                    if( it.le.2 )then
                      write(*,'("macro: assignNonlinearInterfaceGhost24r : it=",i6," t,dt=",2e10.2)') it,t,dt
                    end if
                  end if
                  ! normal and tangent (for TZ forcing)
                  an1=an1Cartesian
                  an2=an2Cartesian
                  tau1=-an2
                  tau2= an1
                 ! make sure the normal and tangent are set
                 if( abs( an1**2 + an2**2 -1. )>1.e-10 .or. abs( tau1**2 + tau2**2 -1. )>1.e-10 )then
                   write(*,'("mla24r - ERROR: incorrect an1,an2, tau1,tau2=",4(1pe9.2))') an1,an2,tau1,tau2
                   stop 6666
                 end if
                 ! --- initialize some forcing functions ---
                 do n=0,nd-1
                   fev1(n)=0.
                   LfE1(n)=0.
                   fEt1(n)=0.
                   fEtt1(n)=0.
                   fev2(n)=0.
                   LfE2(n)=0.
                   fEt2(n)=0.
                   fEtt2(n)=0.
                   fevx1(n)=0.
                   fevy1(n)=0.
                   fevx2(n)=0.
                   fevy2(n)=0.
                   if (dispersionModel1 .ne. 0) then
                     do jv=0,numberOfPolarizationVectors1-1
                       fpv1(n,jv)=0.
                       LfP1(n,jv)=0.
                       fPt1(n,jv)=0.
                       fPtt1(n,jv)=0.
                       fpvx1(n,jv)=0.
                       fpvy1(n,jv)=0.
                     end do
                   endif
                   if (dispersionModel2 .ne. 0) then
                     do jv=0,numberOfPolarizationVectors2-1
                       fpv2(n,jv)=0.
                       LfP2(n,jv)=0.
                       fPt2(n,jv)=0.
                       fPtt2(n,jv)=0.
                       fpvx2(n,jv)=0.
                       fpvy2(n,jv)=0.
                     end do
                   endif
                 end do
                 ! forcing functions for N
                 if (nonlinearModel1 .ne. 0) then
                   do jv = 0,numberOfAtomicLevels1-1
                       fnv1(jv) = 0.
                       fntv1(jv) = 0.
                   enddo
                 endif
                 if (nonlinearModel2 .ne. 0) then
                   do jv = 0,numberOfAtomicLevels2-1
                       fnv2(jv) = 0.
                       fntv2(jv) = 0.
                   enddo
                 endif
                  ! write(*,'("p1=",(15(e10.2,1x)))') (((p1(i1,i2,i3,0),i1=nd1a,nd1b),i2=nd2a,nd2b),i3=nd3a,nd3b)
                  ! =============== start loops ======================
                   i3=n3a
                   j3=m3a
                   j2=m2a
                   do i2=n2a,n2b
                    j1=m1a
                    do i1=n1a,n1b
                     if( mask1(i1,i2,i3).gt.0 .and. mask2(j1,j2,j3).gt.0 )then
                    ! Evaluate the jump conditions using the wrong values at the ghost points 
                     ! Evaluate TZ forcing for dispersive equations in 2D 
                       if( twilightZone.eq.1 )then
                         if( dispersionModel1.ne.noDispersion .and. nonlinearModel1.ne.noNonlinearModel) then
                           nce = pxc+nd*numberOfPolarizationVectors1
                           !-----------------------------
                           ! dimension loop for P and E
                           !-----------------------------
                           nce = pxc+nd*numberOfPolarizationVectors1
                           do n=0,nd-1
                             fpSum1(n)=0.
                             pevttSum1(n)=0.
                             pevttxSum1(n)=0.
                             pevttySum1(n)=0.
                             pevttLSum1(n)=0.
                             pevttttSum1(n)=0.
                             petttSum=0.
                             call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, es(n)   ) 
                             call ogderiv(ep, 1,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, est(n)  )
                             call ogderiv(ep, 2,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, estt(n) )
                             call ogderiv(ep, 0,1,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esx(n) )
                             call ogderiv(ep, 0,0,1,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esy(n) )
                             call ogderiv(ep, 0,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esxx(n) )
                             call ogderiv(ep, 0,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esyy(n) )
                             call ogderiv(ep, 1,1,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, estx(n) )
                             call ogderiv(ep, 1,0,1,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esty(n) )
                             call ogderiv(ep, 2,1,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esttx(n) )
                             call ogderiv(ep, 2,0,1,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, estty(n) )
                             call ogderiv(ep, 0,3,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esxxx(n) )
                             call ogderiv(ep, 0,2,1,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esxxy(n) )
                             call ogderiv(ep, 0,1,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esxyy(n) )
                             call ogderiv(ep, 0,0,3,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esyyy(n) )
                             call ogderiv(ep, 3,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esttt(n) )
                             call ogderiv(ep, 4,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, estttt(n) )
                             call ogderiv(ep, 1,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, estxx(n) )
                             call ogderiv(ep, 1,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, estyy(n) )
                             call ogderiv(ep, 2,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esttxx(n) )
                             call ogderiv(ep, 2,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esttyy(n) )
                             call ogderiv(ep, 0,4,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esxxxx(n) )
                             call ogderiv(ep, 0,2,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esxxyy(n) )
                             call ogderiv(ep, 0,0,4,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esyyyy(n) )
                             ! L = c1^2*Delta
                             esL  = ( esxx(n)   + esyy(n) ) ! deleted c1^2
                             estL = ( estxx(n)  + estyy(n) )
                             esttL= ( esttxx(n) + esttyy(n) )
                             esLx  = ( esxxx(n)   + esxyy(n) )
                             esLy  = ( esxxy(n)   + esyyy(n) )
                             ! L^2 : 
                             esLL = ( esxxxx(n)  + 2.*esxxyy(n) + esyyyy(n) ) ! deleted c1^4
                             do jv=0,numberOfPolarizationVectors1-1
                               ! The TZ component is offset by pxc
                               pc = pxc + jv*nd
                               call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pe(n)   )
                               call ogderiv(ep, 1,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pet(n)  )
                               call ogderiv(ep, 2,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pett(n) )
                               call ogderiv(ep, 3,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pettt(n) )
                               call ogderiv(ep, 4,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, petttt(n) )
                               call ogderiv(ep, 0,1,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pex(n) )
                               call ogderiv(ep, 0,0,1,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pey(n) )
                               call ogderiv(ep, 0,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pexx(n) )
                               call ogderiv(ep, 0,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, peyy(n) )
                               call ogderiv(ep, 1,1,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, petx(n) )
                               call ogderiv(ep, 1,0,1,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pety(n) )
                               call ogderiv(ep, 1,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, petxx(n) )
                               call ogderiv(ep, 1,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, petyy(n) )
                               call ogderiv(ep, 2,1,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pettx(n) )
                               call ogderiv(ep, 2,0,1,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, petty(n) )
                               call ogderiv(ep, 2,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pettxx(n) )
                               call ogderiv(ep, 2,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pettyy(n) )
                               peL  = ( pexx(n)   + peyy(n) ) ! deleted c1^2
                               petL = ( petxx(n)  + petyy(n) )
                               pettL= ( pettxx(n) + pettyy(n) )
                               fpv1(n,jv) = pett(n)   + b1v1(jv)*pet(n)   + b0v1(jv)*pe(n)
                               fPt1(n,jv) = pettt(n)  + b1v1(jv)*pett(n)  + b0v1(jv)*pet(n)
                               fPtt1(n,jv)= petttt(n) + b1v1(jv)*pettt(n) + b0v1(jv)*pett(n)
                               LfP1(n,jv) = pettL     + b1v1(jv)*petL     + b0v1(jv)*peL
                               fpvx1(n,jv)= pettx(n)  + b1v1(jv)*petx(n)  + b0v1(jv)*pex(n)
                               fpvy1(n,jv)= petty(n)  + b1v1(jv)*pety(n)  + b0v1(jv)*pey(n)
                               do na = 0,numberOfAtomicLevels1-1
                                 call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, nce+na, q0  )
                                 call ogderiv(ep, 1,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, nce+na, q0t  )
                                 call ogderiv(ep, 2,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, nce+na, q0tt  )
                                 call ogderiv(ep, 0,1,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, nce+na, q0x  )
                                 call ogderiv(ep, 0,0,1,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, nce+na, q0y  )
                                 call ogderiv(ep, 0,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, nce+na, q0xx  )
                                 call ogderiv(ep, 0,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, nce+na, q0yy  )
                                 fpv1(n,jv) = fpv1(n,jv) - pnec1(jv,na)*q0*es(n) ! adding \Delta N*E
                                 fPt1(n,jv) = fPt1(n,jv) - pnec1(jv,na)*q0t*es(n) - pnec1(jv,na)*q0*est(n)
                                 fPtt1(n,jv) = fPtt1(n,jv) - pnec1(jv,na)*q0tt*es(n)- 2.0*pnec1(jv,na)*q0t*est(n) - pnec1(jv,na)*q0*estt(n)
                                 LfP1(n,jv) = LfP1(n,jv) - pnec1(jv,na)*(q0xx*es(n)+2.*q0x*esx(n)+q0*esxx(n) + q0yy*es(n)+2.*q0y*esy(n)+q0*esyy(n))
                                 fpvx1(n,jv) = fpvx1(n,jv) - pnec1(jv,na)*q0x*es(n) - pnec1(jv,na)*q0*esx(n)
                                 fpvy1(n,jv) = fpvy1(n,jv) - pnec1(jv,na)*q0y*es(n) - pnec1(jv,na)*q0*esy(n)
                               enddo
                               ! print *,'xxxxxxxxxxxxxxxxxxxxxxx'
                               ! print *, 'FOR P TZ FORCING'
                               ! print *, n,jv, fpv1(n,jv),fPt1(n,jv),fPtt1(n,jv),LfP1(n,jv),fpvx1(n,jv),fpvy1(n,jv)
                               ! print *,'xxxxxxxxxxxxxxxxxxxxxxx'
                               ! Normal TZ forcing for P_{n,jv} equation: 
                               ! fpv1(n,jv) = pett(n)   + b1v1(jv)*pet(n)   + b0v1(jv)*pe(n)   - a0v(jv)*es(n)   - a1v(jv)*est(n)
                               ! fPt1(n,jv) = pettt(n)  + b1v1(jv)*pett(n)  + b0v1(jv)*pet(n)  - a0v(jv)*est(n)  - a1v(jv)*estt(n)
                               ! fPtt1(n,jv)= petttt(n) + b1v1(jv)*pettt(n) + b0v1(jv)*pett(n) - a0v(jv)*estt(n) - a1v(jv)*esttt(n)
                               ! LfP1(n,jv) = pettL     + b1v1(jv)*petL     + b0v1(jv)*peL     - a0v(jv)*esL     - a1v(jv)*estL
                               ! fpvx1(n,jv)= pettx(n)  + b1v1(jv)*petx(n)  + b0v1(jv)*pex(n)  - a0v(jv)*esx(n)  - a1v(jv)*estx(n)
                               ! fpvy1(n,jv)= petty(n)  + b1v1(jv)*pety(n)  + b0v1(jv)*pey(n)  - a0v(jv)*esy(n)  - a1v(jv)*esty(n)
                               ! write(*,'(" n=",i4," LfP1=",e10.4," pettL,petL,peL,esL,estL=",5e12.4)') n,LfP1(n,jv),pettL,petL,peL,esL,estL
                               ! write(*,'(" pe,pet,pett,pettt,petttt=",5e12.4)') pe(n),pet(n),pett(n),pettt(n),petttt(n)
                               ! write(*,'("TZ: n,jv=",2i4," pex,pey,pexx,peyy=",4(1pe12.4))') n,jv,pex(n),pey(n),pexx(n),peyy(n)
                               ! Save ptt for checking later
                               pevtt1(n,jv)=pett(n)
                               pevttx1(n,jv)=pettx(n)
                               pevtty1(n,jv)=petty(n)
                               pevttt1(n,jv)=pettt(n)
                               pevtttt1(n,jv)=petttt(n)
                               pevttL1(n,jv) = pettL
                               pevttLSum1(n) = pevttLSum1(n)  + pettL
                               pevttttSum1(n)= pevttttSum1(n) + petttt(n) 
                               ! Keep some sums: 
                               fpSum1(n)   = fpSum1(n)  + fpv1(n,jv)
                               pevttSum1(n)  = pevttSum1(n)  + pett(n) 
                               pevttxSum1(n) = pevttxSum1(n) + pettx(n)
                               pevttySum1(n) = pevttySum1(n) + petty(n)
                               petttSum  = petttSum  + pettt(n) 
                             end do 
                             ! TZ forcing for E_{n} equation:
                             ! E_tt - c1^2 Delta E + alphaP1*Ptt  = 
                             fev1(n) = estt(n)   - c1**2*esL   + alphaP1*pevttSum1(n)
                             fEt1(n) = esttt(n)  - c1**2*estL  + alphaP1*petttSum
                             fEtt1(n)= estttt(n) - c1**2*esttL + alphaP1*pevttttSum1(n)
                             fevx1(n) = esttx(n) - c1**2*esLx   + alphaP1*pevttxSum1(n)
                             fevy1(n) = estty(n) - c1**2*esLy   + alphaP1*pevttySum1(n)
                             ! write(*,'("--> fEtt1=",e10.2," estttt,esttL,pettttSum=",3e10.2)')  fEtt1(n),estttt(n),esttL,pettttSum
                             LfE1(n) = esttL     - c1**2*esLL  + alphaP1*pevttLSum1(n)   
                           end do
                           !--------------------------------
                           ! outside of dimension loop for N
                           !--------------------------------
                           do na=0,numberOfAtomicLevels1-1
                             ! na-th level
                             call ogderiv(ep, 1,0,0,0,xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, nce+na, q0t )
                             call ogderiv(ep, 2,0,0,0,xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, nce+na, q0tt)
                             call ogderiv(ep, 3,0,0,0,xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, nce+na, q0ttt)
                             call ogderiv(ep, 4,0,0,0,xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, nce+na, q0tttt)
                             ! initialize
                             fnv1(na)  = q0t ! forcing for \partial_tN_\ell = alpha_{\ell,k}N_k+\beta_{\ell,m}E\cdot\partial_tP_k
                             fntv1(na) = q0tt ! next derivative
                             fnttv1(na) = q0ttt
                             fntttv1(na) = q0tttt
                             ! relaxation (alpha_{\ell,m})
                             do jv=0,numberOfAtomicLevels1-1
                               call ogderiv(ep, 0,0,0,0,xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, nce+jv, q0 )
                               call ogderiv(ep, 1,0,0,0,xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, nce+jv, q0t)
                               call ogderiv(ep, 2,0,0,0,xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, nce+jv, q0tt)
                               call ogderiv(ep, 3,0,0,0,xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, nce+jv, q0ttt)
                               fnv1(na)  = fnv1(na)  - prc1(na,jv)*q0
                               fntv1(na) = fntv1(na) - prc1(na,jv)*q0t
                               fnttv1(na) = fnttv1(na) - prc1(na,jv)*q0tt
                               fntttv1(na) = fntttv1(na) - prc1(na,jv)*q0ttt
                             enddo
                             ! dot product (\beta_{\ell,k})
                             do n=0,nd-1 ! loop over dim
                               call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, es(n)   ) 
                               call ogderiv(ep, 1,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, est(n)  )
                               call ogderiv(ep, 2,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, estt(n) )
                               call ogderiv(ep, 3,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, esttt(n) )
                               ! corresponding polarization vector
                               do jv=0,numberOfPolarizationVectors1-1 
                                 pc = pxc + jv*nd
                                 call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pe(n)   )
                                 call ogderiv(ep, 1,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pet(n)  )
                                 call ogderiv(ep, 2,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pett(n) )
                                 call ogderiv(ep, 3,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pettt(n) )
                                 call ogderiv(ep, 4,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, petttt(n) ) 
                                 fnv1(na)  = fnv1(na) - peptc1(na,jv)*es(n)*pet(n)
                                 fntv1(na) = fntv1(na) - peptc1(na,jv)*est(n)*pet(n) - peptc1(na,jv)*es(n)*pett(n)
                                 fnttv1(na) = fnttv1(na) - peptc1(na,jv)*estt(n)*pet(n) - 2.d0*peptc1(na,jv)*est(n)*pett(n) - peptc1(na,jv)*es(n)*pettt(n)
                                 fntttv1(na) = fntttv1(na) - peptc1(na,jv)*esttt(n)*pet(n) - 3.d0*peptc1(na,jv)*estt(n)*pett(n) - 3.d0*peptc1(na,jv)*est(n)*pettt(n) - peptc1(na,jv)*es(n)*petttt(n)
                               enddo
                             enddo
                           enddo
                         end if
                         if( dispersionModel2.ne.noDispersion .and. nonlinearModel2.ne.noNonlinearModel) then
                           nce = pxc+nd*numberOfPolarizationVectors2
                           !-----------------------------
                           ! dimension loop for P and E
                           !-----------------------------
                           nce = pxc+nd*numberOfPolarizationVectors2
                           do n=0,nd-1
                             fpSum2(n)=0.
                             pevttSum2(n)=0.
                             pevttxSum2(n)=0.
                             pevttySum2(n)=0.
                             pevttLSum2(n)=0.
                             pevttttSum2(n)=0.
                             petttSum=0.
                             call ogderiv(ep, 0,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, es(n)   ) 
                             call ogderiv(ep, 1,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, est(n)  )
                             call ogderiv(ep, 2,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, estt(n) )
                             call ogderiv(ep, 0,1,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esx(n) )
                             call ogderiv(ep, 0,0,1,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esy(n) )
                             call ogderiv(ep, 0,2,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esxx(n) )
                             call ogderiv(ep, 0,0,2,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esyy(n) )
                             call ogderiv(ep, 1,1,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, estx(n) )
                             call ogderiv(ep, 1,0,1,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esty(n) )
                             call ogderiv(ep, 2,1,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esttx(n) )
                             call ogderiv(ep, 2,0,1,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, estty(n) )
                             call ogderiv(ep, 0,3,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esxxx(n) )
                             call ogderiv(ep, 0,2,1,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esxxy(n) )
                             call ogderiv(ep, 0,1,2,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esxyy(n) )
                             call ogderiv(ep, 0,0,3,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esyyy(n) )
                             call ogderiv(ep, 3,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esttt(n) )
                             call ogderiv(ep, 4,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, estttt(n) )
                             call ogderiv(ep, 1,2,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, estxx(n) )
                             call ogderiv(ep, 1,0,2,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, estyy(n) )
                             call ogderiv(ep, 2,2,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esttxx(n) )
                             call ogderiv(ep, 2,0,2,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esttyy(n) )
                             call ogderiv(ep, 0,4,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esxxxx(n) )
                             call ogderiv(ep, 0,2,2,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esxxyy(n) )
                             call ogderiv(ep, 0,0,4,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esyyyy(n) )
                             ! L = c2^2*Delta
                             esL  = ( esxx(n)   + esyy(n) ) ! deleted c2^2
                             estL = ( estxx(n)  + estyy(n) )
                             esttL= ( esttxx(n) + esttyy(n) )
                             esLx  = ( esxxx(n)   + esxyy(n) )
                             esLy  = ( esxxy(n)   + esyyy(n) )
                             ! L^2 : 
                             esLL = ( esxxxx(n)  + 2.*esxxyy(n) + esyyyy(n) ) ! deleted c2^4
                             do jv=0,numberOfPolarizationVectors2-1
                               ! The TZ component is offset by pxc
                               pc = pxc + jv*nd
                               call ogderiv(ep, 0,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pe(n)   )
                               call ogderiv(ep, 1,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pet(n)  )
                               call ogderiv(ep, 2,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pett(n) )
                               call ogderiv(ep, 3,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pettt(n) )
                               call ogderiv(ep, 4,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, petttt(n) )
                               call ogderiv(ep, 0,1,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pex(n) )
                               call ogderiv(ep, 0,0,1,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pey(n) )
                               call ogderiv(ep, 0,2,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pexx(n) )
                               call ogderiv(ep, 0,0,2,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, peyy(n) )
                               call ogderiv(ep, 1,1,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, petx(n) )
                               call ogderiv(ep, 1,0,1,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pety(n) )
                               call ogderiv(ep, 1,2,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, petxx(n) )
                               call ogderiv(ep, 1,0,2,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, petyy(n) )
                               call ogderiv(ep, 2,1,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pettx(n) )
                               call ogderiv(ep, 2,0,1,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, petty(n) )
                               call ogderiv(ep, 2,2,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pettxx(n) )
                               call ogderiv(ep, 2,0,2,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pettyy(n) )
                               peL  = ( pexx(n)   + peyy(n) ) ! deleted c2^2
                               petL = ( petxx(n)  + petyy(n) )
                               pettL= ( pettxx(n) + pettyy(n) )
                               fpv2(n,jv) = pett(n)   + b1v2(jv)*pet(n)   + b0v2(jv)*pe(n)
                               fPt2(n,jv) = pettt(n)  + b1v2(jv)*pett(n)  + b0v2(jv)*pet(n)
                               fPtt2(n,jv)= petttt(n) + b1v2(jv)*pettt(n) + b0v2(jv)*pett(n)
                               LfP2(n,jv) = pettL     + b1v2(jv)*petL     + b0v2(jv)*peL
                               fpvx2(n,jv)= pettx(n)  + b1v2(jv)*petx(n)  + b0v2(jv)*pex(n)
                               fpvy2(n,jv)= petty(n)  + b1v2(jv)*pety(n)  + b0v2(jv)*pey(n)
                               do na = 0,numberOfAtomicLevels2-1
                                 call ogderiv(ep, 0,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t, nce+na, q0  )
                                 call ogderiv(ep, 1,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t, nce+na, q0t  )
                                 call ogderiv(ep, 2,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t, nce+na, q0tt  )
                                 call ogderiv(ep, 0,1,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t, nce+na, q0x  )
                                 call ogderiv(ep, 0,0,1,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t, nce+na, q0y  )
                                 call ogderiv(ep, 0,2,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t, nce+na, q0xx  )
                                 call ogderiv(ep, 0,0,2,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t, nce+na, q0yy  )
                                 fpv2(n,jv) = fpv2(n,jv) - pnec2(jv,na)*q0*es(n) ! adding \Delta N*E
                                 fPt2(n,jv) = fPt2(n,jv) - pnec2(jv,na)*q0t*es(n) - pnec2(jv,na)*q0*est(n)
                                 fPtt2(n,jv) = fPtt2(n,jv) - pnec2(jv,na)*q0tt*es(n)- 2.0*pnec2(jv,na)*q0t*est(n) - pnec2(jv,na)*q0*estt(n)
                                 LfP2(n,jv) = LfP2(n,jv) - pnec2(jv,na)*(q0xx*es(n)+2.*q0x*esx(n)+q0*esxx(n) + q0yy*es(n)+2.*q0y*esy(n)+q0*esyy(n))
                                 fpvx2(n,jv) = fpvx2(n,jv) - pnec2(jv,na)*q0x*es(n) - pnec2(jv,na)*q0*esx(n)
                                 fpvy2(n,jv) = fpvy2(n,jv) - pnec2(jv,na)*q0y*es(n) - pnec2(jv,na)*q0*esy(n)
                               enddo
                               ! print *,'xxxxxxxxxxxxxxxxxxxxxxx'
                               ! print *, 'FOR P TZ FORCING'
                               ! print *, n,jv, fpv2(n,jv),fPt2(n,jv),fPtt2(n,jv),LfP2(n,jv),fpvx2(n,jv),fpvy2(n,jv)
                               ! print *,'xxxxxxxxxxxxxxxxxxxxxxx'
                               ! Normal TZ forcing for P_{n,jv} equation: 
                               ! fpv2(n,jv) = pett(n)   + b1v2(jv)*pet(n)   + b0v2(jv)*pe(n)   - a0v(jv)*es(n)   - a1v(jv)*est(n)
                               ! fPt2(n,jv) = pettt(n)  + b1v2(jv)*pett(n)  + b0v2(jv)*pet(n)  - a0v(jv)*est(n)  - a1v(jv)*estt(n)
                               ! fPtt2(n,jv)= petttt(n) + b1v2(jv)*pettt(n) + b0v2(jv)*pett(n) - a0v(jv)*estt(n) - a1v(jv)*esttt(n)
                               ! LfP2(n,jv) = pettL     + b1v2(jv)*petL     + b0v2(jv)*peL     - a0v(jv)*esL     - a1v(jv)*estL
                               ! fpvx2(n,jv)= pettx(n)  + b1v2(jv)*petx(n)  + b0v2(jv)*pex(n)  - a0v(jv)*esx(n)  - a1v(jv)*estx(n)
                               ! fpvy2(n,jv)= petty(n)  + b1v2(jv)*pety(n)  + b0v2(jv)*pey(n)  - a0v(jv)*esy(n)  - a1v(jv)*esty(n)
                               ! write(*,'(" n=",i4," LfP2=",e10.4," pettL,petL,peL,esL,estL=",5e12.4)') n,LfP2(n,jv),pettL,petL,peL,esL,estL
                               ! write(*,'(" pe,pet,pett,pettt,petttt=",5e12.4)') pe(n),pet(n),pett(n),pettt(n),petttt(n)
                               ! write(*,'("TZ: n,jv=",2i4," pex,pey,pexx,peyy=",4(1pe12.4))') n,jv,pex(n),pey(n),pexx(n),peyy(n)
                               ! Save ptt for checking later
                               pevtt2(n,jv)=pett(n)
                               pevttx2(n,jv)=pettx(n)
                               pevtty2(n,jv)=petty(n)
                               pevttt2(n,jv)=pettt(n)
                               pevtttt2(n,jv)=petttt(n)
                               pevttL2(n,jv) = pettL
                               pevttLSum2(n) = pevttLSum2(n)  + pettL
                               pevttttSum2(n)= pevttttSum2(n) + petttt(n) 
                               ! Keep some sums: 
                               fpSum2(n)   = fpSum2(n)  + fpv2(n,jv)
                               pevttSum2(n)  = pevttSum2(n)  + pett(n) 
                               pevttxSum2(n) = pevttxSum2(n) + pettx(n)
                               pevttySum2(n) = pevttySum2(n) + petty(n)
                               petttSum  = petttSum  + pettt(n) 
                             end do 
                             ! TZ forcing for E_{n} equation:
                             ! E_tt - c2^2 Delta E + alphaP2*Ptt  = 
                             fev2(n) = estt(n)   - c2**2*esL   + alphaP2*pevttSum2(n)
                             fEt2(n) = esttt(n)  - c2**2*estL  + alphaP2*petttSum
                             fEtt2(n)= estttt(n) - c2**2*esttL + alphaP2*pevttttSum2(n)
                             fevx2(n) = esttx(n) - c2**2*esLx   + alphaP2*pevttxSum2(n)
                             fevy2(n) = estty(n) - c2**2*esLy   + alphaP2*pevttySum2(n)
                             ! write(*,'("--> fEtt2=",e10.2," estttt,esttL,pettttSum=",3e10.2)')  fEtt2(n),estttt(n),esttL,pettttSum
                             LfE2(n) = esttL     - c2**2*esLL  + alphaP2*pevttLSum2(n)   
                           end do
                           !--------------------------------
                           ! outside of dimension loop for N
                           !--------------------------------
                           do na=0,numberOfAtomicLevels2-1
                             ! na-th level
                             call ogderiv(ep, 1,0,0,0,xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t, nce+na, q0t )
                             call ogderiv(ep, 2,0,0,0,xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t, nce+na, q0tt)
                             call ogderiv(ep, 3,0,0,0,xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t, nce+na, q0ttt)
                             call ogderiv(ep, 4,0,0,0,xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t, nce+na, q0tttt)
                             ! initialize
                             fnv2(na)  = q0t ! forcing for \partial_tN_\ell = alpha_{\ell,k}N_k+\beta_{\ell,m}E\cdot\partial_tP_k
                             fntv2(na) = q0tt ! next derivative
                             fnttv2(na) = q0ttt
                             fntttv2(na) = q0tttt
                             ! relaxation (alpha_{\ell,m})
                             do jv=0,numberOfAtomicLevels2-1
                               call ogderiv(ep, 0,0,0,0,xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t, nce+jv, q0 )
                               call ogderiv(ep, 1,0,0,0,xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t, nce+jv, q0t)
                               call ogderiv(ep, 2,0,0,0,xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t, nce+jv, q0tt)
                               call ogderiv(ep, 3,0,0,0,xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t, nce+jv, q0ttt)
                               fnv2(na)  = fnv2(na)  - prc2(na,jv)*q0
                               fntv2(na) = fntv2(na) - prc2(na,jv)*q0t
                               fnttv2(na) = fnttv2(na) - prc2(na,jv)*q0tt
                               fntttv2(na) = fntttv2(na) - prc2(na,jv)*q0ttt
                             enddo
                             ! dot product (\beta_{\ell,k})
                             do n=0,nd-1 ! loop over dim
                               call ogderiv(ep, 0,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, es(n)   ) 
                               call ogderiv(ep, 1,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, est(n)  )
                               call ogderiv(ep, 2,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, estt(n) )
                               call ogderiv(ep, 3,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, esttt(n) )
                               ! corresponding polarization vector
                               do jv=0,numberOfPolarizationVectors2-1 
                                 pc = pxc + jv*nd
                                 call ogderiv(ep, 0,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pe(n)   )
                                 call ogderiv(ep, 1,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pet(n)  )
                                 call ogderiv(ep, 2,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pett(n) )
                                 call ogderiv(ep, 3,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pettt(n) )
                                 call ogderiv(ep, 4,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, petttt(n) ) 
                                 fnv2(na)  = fnv2(na) - peptc2(na,jv)*es(n)*pet(n)
                                 fntv2(na) = fntv2(na) - peptc2(na,jv)*est(n)*pet(n) - peptc2(na,jv)*es(n)*pett(n)
                                 fnttv2(na) = fnttv2(na) - peptc2(na,jv)*estt(n)*pet(n) - 2.d0*peptc2(na,jv)*est(n)*pett(n) - peptc2(na,jv)*es(n)*pettt(n)
                                 fntttv2(na) = fntttv2(na) - peptc2(na,jv)*esttt(n)*pet(n) - 3.d0*peptc2(na,jv)*estt(n)*pett(n) - 3.d0*peptc2(na,jv)*est(n)*pettt(n) - peptc2(na,jv)*es(n)*petttt(n)
                               enddo
                             enddo
                           enddo
                         end if
                       end if
                      ! These derivatives are computed to 2nd-order accuracy
                      ! NOTE: the jacobian derivatives can be computed once for all components
                        uu1=u1(i1,i2,i3,ex) ! in the rectangular case just eval the solution
                         u1xxx = (-u1(i1-2,i2,i3,ex)+2.*u1(i1-1,i2,i3,ex)-2.*u1(i1+1,i2,i3,ex)+u1(i1+2,i2,i3,ex))/(2.*dx1(0)**3)
                         u1xxy = ((-u1(i1-1,i2-1,i3,ex)+u1(i1-1,i2+1,i3,ex))/(2.*dx1(1))-2.*(-u1(i1,i2-1,i3,ex)+u1(i1,i2+1,i3,ex))/(2.*dx1(1))+(-u1(i1+1,i2-1,i3,ex)+u1(i1+1,i2+1,i3,ex))/(2.*dx1(1)))/(dx1(0)**2)
                         u1xyy = (-(u1(i1-1,i2-1,i3,ex)-2.*u1(i1-1,i2,i3,ex)+u1(i1-1,i2+1,i3,ex))/(dx1(1)**2)+(u1(i1+1,i2-1,i3,ex)-2.*u1(i1+1,i2,i3,ex)+u1(i1+1,i2+1,i3,ex))/(dx1(1)**2))/(2.*dx1(0))
                         u1yyy = (-u1(i1,i2-2,i3,ex)+2.*u1(i1,i2-1,i3,ex)-2.*u1(i1,i2+1,i3,ex)+u1(i1,i2+2,i3,ex))/(2.*dx1(1)**3)
                         u1xxxx = (u1(i1-2,i2,i3,ex)-4.*u1(i1-1,i2,i3,ex)+6.*u1(i1,i2,i3,ex)-4.*u1(i1+1,i2,i3,ex)+u1(i1+2,i2,i3,ex))/(dx1(0)**4)
                         u1xxyy = ((u1(i1-1,i2-1,i3,ex)-2.*u1(i1-1,i2,i3,ex)+u1(i1-1,i2+1,i3,ex))/(dx1(1)**2)-2.*(u1(i1,i2-1,i3,ex)-2.*u1(i1,i2,i3,ex)+u1(i1,i2+1,i3,ex))/(dx1(1)**2)+(u1(i1+1,i2-1,i3,ex)-2.*u1(i1+1,i2,i3,ex)+u1(i1+1,i2+1,i3,ex))/(dx1(1)**2))/(dx1(0)**2)
                         u1yyyy = (u1(i1,i2-2,i3,ex)-4.*u1(i1,i2-1,i3,ex)+6.*u1(i1,i2,i3,ex)-4.*u1(i1,i2+1,i3,ex)+u1(i1,i2+2,i3,ex))/(dx1(1)**4)
                       u1LapSq = u1xxxx +2.* u1xxyy + u1yyyy
                        vv1=u1(i1,i2,i3,ey) ! in the rectangular case just eval the solution
                         v1xxx = (-u1(i1-2,i2,i3,ey)+2.*u1(i1-1,i2,i3,ey)-2.*u1(i1+1,i2,i3,ey)+u1(i1+2,i2,i3,ey))/(2.*dx1(0)**3)
                         v1xxy = ((-u1(i1-1,i2-1,i3,ey)+u1(i1-1,i2+1,i3,ey))/(2.*dx1(1))-2.*(-u1(i1,i2-1,i3,ey)+u1(i1,i2+1,i3,ey))/(2.*dx1(1))+(-u1(i1+1,i2-1,i3,ey)+u1(i1+1,i2+1,i3,ey))/(2.*dx1(1)))/(dx1(0)**2)
                         v1xyy = (-(u1(i1-1,i2-1,i3,ey)-2.*u1(i1-1,i2,i3,ey)+u1(i1-1,i2+1,i3,ey))/(dx1(1)**2)+(u1(i1+1,i2-1,i3,ey)-2.*u1(i1+1,i2,i3,ey)+u1(i1+1,i2+1,i3,ey))/(dx1(1)**2))/(2.*dx1(0))
                         v1yyy = (-u1(i1,i2-2,i3,ey)+2.*u1(i1,i2-1,i3,ey)-2.*u1(i1,i2+1,i3,ey)+u1(i1,i2+2,i3,ey))/(2.*dx1(1)**3)
                         v1xxxx = (u1(i1-2,i2,i3,ey)-4.*u1(i1-1,i2,i3,ey)+6.*u1(i1,i2,i3,ey)-4.*u1(i1+1,i2,i3,ey)+u1(i1+2,i2,i3,ey))/(dx1(0)**4)
                         v1xxyy = ((u1(i1-1,i2-1,i3,ey)-2.*u1(i1-1,i2,i3,ey)+u1(i1-1,i2+1,i3,ey))/(dx1(1)**2)-2.*(u1(i1,i2-1,i3,ey)-2.*u1(i1,i2,i3,ey)+u1(i1,i2+1,i3,ey))/(dx1(1)**2)+(u1(i1+1,i2-1,i3,ey)-2.*u1(i1+1,i2,i3,ey)+u1(i1+1,i2+1,i3,ey))/(dx1(1)**2))/(dx1(0)**2)
                         v1yyyy = (u1(i1,i2-2,i3,ey)-4.*u1(i1,i2-1,i3,ey)+6.*u1(i1,i2,i3,ey)-4.*u1(i1,i2+1,i3,ey)+u1(i1,i2+2,i3,ey))/(dx1(1)**4)
                       v1LapSq = v1xxxx +2.* v1xxyy + v1yyyy
                      ! NOTE: the jacobian derivatives can be computed once for all components
                        uu2=u2(j1,j2,j3,ex) ! in the rectangular case just eval the solution
                         u2xxx = (-u2(j1-2,j2,j3,ex)+2.*u2(j1-1,j2,j3,ex)-2.*u2(j1+1,j2,j3,ex)+u2(j1+2,j2,j3,ex))/(2.*dx2(0)**3)
                         u2xxy = ((-u2(j1-1,j2-1,j3,ex)+u2(j1-1,j2+1,j3,ex))/(2.*dx2(1))-2.*(-u2(j1,j2-1,j3,ex)+u2(j1,j2+1,j3,ex))/(2.*dx2(1))+(-u2(j1+1,j2-1,j3,ex)+u2(j1+1,j2+1,j3,ex))/(2.*dx2(1)))/(dx2(0)**2)
                         u2xyy = (-(u2(j1-1,j2-1,j3,ex)-2.*u2(j1-1,j2,j3,ex)+u2(j1-1,j2+1,j3,ex))/(dx2(1)**2)+(u2(j1+1,j2-1,j3,ex)-2.*u2(j1+1,j2,j3,ex)+u2(j1+1,j2+1,j3,ex))/(dx2(1)**2))/(2.*dx2(0))
                         u2yyy = (-u2(j1,j2-2,j3,ex)+2.*u2(j1,j2-1,j3,ex)-2.*u2(j1,j2+1,j3,ex)+u2(j1,j2+2,j3,ex))/(2.*dx2(1)**3)
                         u2xxxx = (u2(j1-2,j2,j3,ex)-4.*u2(j1-1,j2,j3,ex)+6.*u2(j1,j2,j3,ex)-4.*u2(j1+1,j2,j3,ex)+u2(j1+2,j2,j3,ex))/(dx2(0)**4)
                         u2xxyy = ((u2(j1-1,j2-1,j3,ex)-2.*u2(j1-1,j2,j3,ex)+u2(j1-1,j2+1,j3,ex))/(dx2(1)**2)-2.*(u2(j1,j2-1,j3,ex)-2.*u2(j1,j2,j3,ex)+u2(j1,j2+1,j3,ex))/(dx2(1)**2)+(u2(j1+1,j2-1,j3,ex)-2.*u2(j1+1,j2,j3,ex)+u2(j1+1,j2+1,j3,ex))/(dx2(1)**2))/(dx2(0)**2)
                         u2yyyy = (u2(j1,j2-2,j3,ex)-4.*u2(j1,j2-1,j3,ex)+6.*u2(j1,j2,j3,ex)-4.*u2(j1,j2+1,j3,ex)+u2(j1,j2+2,j3,ex))/(dx2(1)**4)
                       u2LapSq = u2xxxx +2.* u2xxyy + u2yyyy
                        vv2=u2(j1,j2,j3,ey) ! in the rectangular case just eval the solution
                         v2xxx = (-u2(j1-2,j2,j3,ey)+2.*u2(j1-1,j2,j3,ey)-2.*u2(j1+1,j2,j3,ey)+u2(j1+2,j2,j3,ey))/(2.*dx2(0)**3)
                         v2xxy = ((-u2(j1-1,j2-1,j3,ey)+u2(j1-1,j2+1,j3,ey))/(2.*dx2(1))-2.*(-u2(j1,j2-1,j3,ey)+u2(j1,j2+1,j3,ey))/(2.*dx2(1))+(-u2(j1+1,j2-1,j3,ey)+u2(j1+1,j2+1,j3,ey))/(2.*dx2(1)))/(dx2(0)**2)
                         v2xyy = (-(u2(j1-1,j2-1,j3,ey)-2.*u2(j1-1,j2,j3,ey)+u2(j1-1,j2+1,j3,ey))/(dx2(1)**2)+(u2(j1+1,j2-1,j3,ey)-2.*u2(j1+1,j2,j3,ey)+u2(j1+1,j2+1,j3,ey))/(dx2(1)**2))/(2.*dx2(0))
                         v2yyy = (-u2(j1,j2-2,j3,ey)+2.*u2(j1,j2-1,j3,ey)-2.*u2(j1,j2+1,j3,ey)+u2(j1,j2+2,j3,ey))/(2.*dx2(1)**3)
                         v2xxxx = (u2(j1-2,j2,j3,ey)-4.*u2(j1-1,j2,j3,ey)+6.*u2(j1,j2,j3,ey)-4.*u2(j1+1,j2,j3,ey)+u2(j1+2,j2,j3,ey))/(dx2(0)**4)
                         v2xxyy = ((u2(j1-1,j2-1,j3,ey)-2.*u2(j1-1,j2,j3,ey)+u2(j1-1,j2+1,j3,ey))/(dx2(1)**2)-2.*(u2(j1,j2-1,j3,ey)-2.*u2(j1,j2,j3,ey)+u2(j1,j2+1,j3,ey))/(dx2(1)**2)+(u2(j1+1,j2-1,j3,ey)-2.*u2(j1+1,j2,j3,ey)+u2(j1+1,j2+1,j3,ey))/(dx2(1)**2))/(dx2(0)**2)
                         v2yyyy = (u2(j1,j2-2,j3,ey)-4.*u2(j1,j2-1,j3,ey)+6.*u2(j1,j2,j3,ey)-4.*u2(j1,j2+1,j3,ey)+u2(j1,j2+2,j3,ey))/(dx2(1)**4)
                       v2LapSq = v2xxxx +2.* v2xxyy + v2yyyy
                      ! These derivatives are computed to 4th-order accuracy
                      ! NOTE: the jacobian derivatives can be computed once for all components
                        uu1=u1(i1,i2,i3,ex) ! in the rectangular case just eval the solution
                         u1x = (u1(i1-2,i2,i3,ex)-8.*u1(i1-1,i2,i3,ex)+8.*u1(i1+1,i2,i3,ex)-u1(i1+2,i2,i3,ex))/(12.*dx1(0))
                         u1y = (u1(i1,i2-2,i3,ex)-8.*u1(i1,i2-1,i3,ex)+8.*u1(i1,i2+1,i3,ex)-u1(i1,i2+2,i3,ex))/(12.*dx1(1))
                         u1xx = (-u1(i1-2,i2,i3,ex)+16.*u1(i1-1,i2,i3,ex)-30.*u1(i1,i2,i3,ex)+16.*u1(i1+1,i2,i3,ex)-u1(i1+2,i2,i3,ex))/(12.*dx1(0)**2)
                         u1yy = (-u1(i1,i2-2,i3,ex)+16.*u1(i1,i2-1,i3,ex)-30.*u1(i1,i2,i3,ex)+16.*u1(i1,i2+1,i3,ex)-u1(i1,i2+2,i3,ex))/(12.*dx1(1)**2)
                       u1Lap = u1xx+ u1yy
                        vv1=u1(i1,i2,i3,ey) ! in the rectangular case just eval the solution
                         v1x = (u1(i1-2,i2,i3,ey)-8.*u1(i1-1,i2,i3,ey)+8.*u1(i1+1,i2,i3,ey)-u1(i1+2,i2,i3,ey))/(12.*dx1(0))
                         v1y = (u1(i1,i2-2,i3,ey)-8.*u1(i1,i2-1,i3,ey)+8.*u1(i1,i2+1,i3,ey)-u1(i1,i2+2,i3,ey))/(12.*dx1(1))
                         v1xx = (-u1(i1-2,i2,i3,ey)+16.*u1(i1-1,i2,i3,ey)-30.*u1(i1,i2,i3,ey)+16.*u1(i1+1,i2,i3,ey)-u1(i1+2,i2,i3,ey))/(12.*dx1(0)**2)
                         v1yy = (-u1(i1,i2-2,i3,ey)+16.*u1(i1,i2-1,i3,ey)-30.*u1(i1,i2,i3,ey)+16.*u1(i1,i2+1,i3,ey)-u1(i1,i2+2,i3,ey))/(12.*dx1(1)**2)
                       v1Lap = v1xx+ v1yy
                      ! NOTE: the jacobian derivatives can be computed once for all components
                        uu2=u2(j1,j2,j3,ex) ! in the rectangular case just eval the solution
                         u2x = (u2(j1-2,j2,j3,ex)-8.*u2(j1-1,j2,j3,ex)+8.*u2(j1+1,j2,j3,ex)-u2(j1+2,j2,j3,ex))/(12.*dx2(0))
                         u2y = (u2(j1,j2-2,j3,ex)-8.*u2(j1,j2-1,j3,ex)+8.*u2(j1,j2+1,j3,ex)-u2(j1,j2+2,j3,ex))/(12.*dx2(1))
                         u2xx = (-u2(j1-2,j2,j3,ex)+16.*u2(j1-1,j2,j3,ex)-30.*u2(j1,j2,j3,ex)+16.*u2(j1+1,j2,j3,ex)-u2(j1+2,j2,j3,ex))/(12.*dx2(0)**2)
                         u2yy = (-u2(j1,j2-2,j3,ex)+16.*u2(j1,j2-1,j3,ex)-30.*u2(j1,j2,j3,ex)+16.*u2(j1,j2+1,j3,ex)-u2(j1,j2+2,j3,ex))/(12.*dx2(1)**2)
                       u2Lap = u2xx+ u2yy
                        vv2=u2(j1,j2,j3,ey) ! in the rectangular case just eval the solution
                         v2x = (u2(j1-2,j2,j3,ey)-8.*u2(j1-1,j2,j3,ey)+8.*u2(j1+1,j2,j3,ey)-u2(j1+2,j2,j3,ey))/(12.*dx2(0))
                         v2y = (u2(j1,j2-2,j3,ey)-8.*u2(j1,j2-1,j3,ey)+8.*u2(j1,j2+1,j3,ey)-u2(j1,j2+2,j3,ey))/(12.*dx2(1))
                         v2xx = (-u2(j1-2,j2,j3,ey)+16.*u2(j1-1,j2,j3,ey)-30.*u2(j1,j2,j3,ey)+16.*u2(j1+1,j2,j3,ey)-u2(j1+2,j2,j3,ey))/(12.*dx2(0)**2)
                         v2yy = (-u2(j1,j2-2,j3,ey)+16.*u2(j1,j2-1,j3,ey)-30.*u2(j1,j2,j3,ey)+16.*u2(j1,j2+1,j3,ey)-u2(j1,j2+2,j3,ey))/(12.*dx2(1)**2)
                       v2Lap = v2xx+ v2yy
                      ! Store c^2*Delta(E) in a vector 
                      LE1(0)=(c1**2)*u1Lap
                      LE1(1)=(c1**2)*v1Lap
                      LE2(0)=(c2**2)*u2Lap
                      LE2(1)=(c2**2)*v2Lap
                      ! Store L^2(E) 
                      LLE1(0)=(c1**4)*u1LapSq
                      LLE1(1)=(c1**4)*v1LapSq
                      LLE2(0)=(c2**4)*u2LapSq
                      LLE2(1)=(c2**4)*v2LapSq
                      ! Store (LE).x an (LE).y 
                      LEx1(0) = (c1**2)*( u1xxx + u1xyy )
                      LEx1(1) = (c1**2)*( v1xxx + v1xyy )
                      LEy1(0) = (c1**2)*( u1xxy + u1yyy )
                      LEy1(1) = (c1**2)*( v1xxy + v1yyy )
                      LEx2(0) = (c2**2)*( u2xxx + u2xyy )
                      LEx2(1) = (c2**2)*( v2xxx + v2xyy )
                      LEy2(0) = (c2**2)*( u2xxy + u2yyy )
                      LEy2(1) = (c2**2)*( v2xxy + v2yyy )
                      ! We also need derivatives at the old time:
                      ! These next derivatives may only be needed to order2, but use order 4 for now so exact for degree 4
                      ! #perl $ORDER=2;
                        uu1=u1n(i1,i2,i3,ex) ! in the rectangular case just eval the solution
                         u1nx = (u1n(i1-2,i2,i3,ex)-8.*u1n(i1-1,i2,i3,ex)+8.*u1n(i1+1,i2,i3,ex)-u1n(i1+2,i2,i3,ex))/(12.*dx1(0))
                         u1ny = (u1n(i1,i2-2,i3,ex)-8.*u1n(i1,i2-1,i3,ex)+8.*u1n(i1,i2+1,i3,ex)-u1n(i1,i2+2,i3,ex))/(12.*dx1(1))
                         u1nxx = (-u1n(i1-2,i2,i3,ex)+16.*u1n(i1-1,i2,i3,ex)-30.*u1n(i1,i2,i3,ex)+16.*u1n(i1+1,i2,i3,ex)-u1n(i1+2,i2,i3,ex))/(12.*dx1(0)**2)
                         u1nyy = (-u1n(i1,i2-2,i3,ex)+16.*u1n(i1,i2-1,i3,ex)-30.*u1n(i1,i2,i3,ex)+16.*u1n(i1,i2+1,i3,ex)-u1n(i1,i2+2,i3,ex))/(12.*dx1(1)**2)
                       u1nLap = u1nxx+ u1nyy
                        vv1=u1n(i1,i2,i3,ey) ! in the rectangular case just eval the solution
                         v1nx = (u1n(i1-2,i2,i3,ey)-8.*u1n(i1-1,i2,i3,ey)+8.*u1n(i1+1,i2,i3,ey)-u1n(i1+2,i2,i3,ey))/(12.*dx1(0))
                         v1ny = (u1n(i1,i2-2,i3,ey)-8.*u1n(i1,i2-1,i3,ey)+8.*u1n(i1,i2+1,i3,ey)-u1n(i1,i2+2,i3,ey))/(12.*dx1(1))
                         v1nxx = (-u1n(i1-2,i2,i3,ey)+16.*u1n(i1-1,i2,i3,ey)-30.*u1n(i1,i2,i3,ey)+16.*u1n(i1+1,i2,i3,ey)-u1n(i1+2,i2,i3,ey))/(12.*dx1(0)**2)
                         v1nyy = (-u1n(i1,i2-2,i3,ey)+16.*u1n(i1,i2-1,i3,ey)-30.*u1n(i1,i2,i3,ey)+16.*u1n(i1,i2+1,i3,ey)-u1n(i1,i2+2,i3,ey))/(12.*dx1(1)**2)
                       v1nLap = v1nxx+ v1nyy
                      ! Here are c^2*Delta(E) at the old time: 
                      LE1m(0) = (c1**2)*u1nLap
                      LE1m(1) = (c1**2)*v1nLap
                        uu2=u2n(j1,j2,j3,ex) ! in the rectangular case just eval the solution
                         u2nx = (u2n(j1-2,j2,j3,ex)-8.*u2n(j1-1,j2,j3,ex)+8.*u2n(j1+1,j2,j3,ex)-u2n(j1+2,j2,j3,ex))/(12.*dx2(0))
                         u2ny = (u2n(j1,j2-2,j3,ex)-8.*u2n(j1,j2-1,j3,ex)+8.*u2n(j1,j2+1,j3,ex)-u2n(j1,j2+2,j3,ex))/(12.*dx2(1))
                         u2nxx = (-u2n(j1-2,j2,j3,ex)+16.*u2n(j1-1,j2,j3,ex)-30.*u2n(j1,j2,j3,ex)+16.*u2n(j1+1,j2,j3,ex)-u2n(j1+2,j2,j3,ex))/(12.*dx2(0)**2)
                         u2nyy = (-u2n(j1,j2-2,j3,ex)+16.*u2n(j1,j2-1,j3,ex)-30.*u2n(j1,j2,j3,ex)+16.*u2n(j1,j2+1,j3,ex)-u2n(j1,j2+2,j3,ex))/(12.*dx2(1)**2)
                       u2nLap = u2nxx+ u2nyy
                        vv2=u2n(j1,j2,j3,ey) ! in the rectangular case just eval the solution
                         v2nx = (u2n(j1-2,j2,j3,ey)-8.*u2n(j1-1,j2,j3,ey)+8.*u2n(j1+1,j2,j3,ey)-u2n(j1+2,j2,j3,ey))/(12.*dx2(0))
                         v2ny = (u2n(j1,j2-2,j3,ey)-8.*u2n(j1,j2-1,j3,ey)+8.*u2n(j1,j2+1,j3,ey)-u2n(j1,j2+2,j3,ey))/(12.*dx2(1))
                         v2nxx = (-u2n(j1-2,j2,j3,ey)+16.*u2n(j1-1,j2,j3,ey)-30.*u2n(j1,j2,j3,ey)+16.*u2n(j1+1,j2,j3,ey)-u2n(j1+2,j2,j3,ey))/(12.*dx2(0)**2)
                         v2nyy = (-u2n(j1,j2-2,j3,ey)+16.*u2n(j1,j2-1,j3,ey)-30.*u2n(j1,j2,j3,ey)+16.*u2n(j1,j2+1,j3,ey)-u2n(j1,j2+2,j3,ey))/(12.*dx2(1)**2)
                       v2nLap = v2nxx+ v2nyy
                      LE2m(0) = (c2**2)*u2nLap
                      LE2m(1) = (c2**2)*v2nLap
                      evx1(0) = u1x
                      evx1(1) = v1x
                      evy1(0) = u1y
                      evy1(1) = v1y 
                      evnx1(0) = u1nx
                      evnx1(1) = v1nx
                      evny1(0) = u1ny
                      evny1(1) = v1ny 
                      evx2(0) = u2x
                      evx2(1) = v2x
                      evy2(0) = u2y
                      evy2(1) = v2y 
                      evnx2(0) = u2nx
                      evnx2(1) = v2nx
                      evny2(0) = u2ny
                      evny2(1) = v2ny 
                     ! eval nonlinear dispersive forcings for domain 1
                       ! pre-assign 0 values
                       do n=0,nd-1 ! dispersive forcing in jump conditions
                         fp1(n)=0.
                         fPttx1(n) =0.
                         fPtty1(n) =0.
                         fLPtt1(n) =0.
                         fPtttt1(n)=0.
                       end do
                       ! only do this for MLA (dispersive and nonlinear multi-level)
                       if( dispersionModel1.ne.noDispersion .and. nonlinearModel1.ne.noNonlinearModel) then
                         nce = pxc+nd*numberOfPolarizationVectors1
                         ! -----------------------------------------
                         ! order 2 (E, P, N) at the interface (fictitious step)
                         !------------------------------------------
                         ! dimension loop for E and P
                         do n=0,nd-1
                           ec = ex +n 
                           ev0  =  u1n(i1,i2,i3,ec)
                           ev   =  u1(i1,i2,i3,ec) ! time where we need to fill in ghost points
                           pSum=0.
                           do jv=0,numberOfPolarizationVectors1-1
                             pc = n + jv*nd  
                             pv0 =  p1n(i1,i2,i3,pc)
                             pv  =  p1(i1,i2,i3,pc)
                             pvn = 2.*pv-p1n(i1,i2,i3,pc) + 0.5*dt*b1v1(jv)*p1n(i1,i2,i3,pc) - dtsq*b0v1(jv)*pv + dtsq*fpv1(n,jv)
                             do na = 0,numberOfAtomicLevels1-1 ! \Delta N^n*E^n
                               pvn = pvn + dtsq*pnec1(jv,na)*q1(i1,i2,i3,na)*ev
                             enddo ! na
                             pvec(n,jv)= pvn/( 1.+.5*dt*b1v1(jv) ) ! time + dt
                             ! #If "p1" eq "p1"
                             ! call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t+dt,pxc+jv*nd+n, pe(n)   )
                             ! ! call ogderiv(ep, 1,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pet(n)   )
                             ! ! call ogderiv(ep, 2,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pett(n)   )
                             ! #Else
                             ! call ogderiv(ep, 0,0,0,0, xy2(i1,i2,i3,0),xy2(i1,i2,i3,1),0.,t+dt,pxc+jv*nd+n, pe(n)   )
                             ! ! call ogderiv(ep, 1,0,0,0, xy2(i1,i2,i3,0),xy2(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pet(n)   )
                             ! ! call ogderiv(ep, 2,0,0,0, xy2(i1,i2,i3,0),xy2(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pett(n)   )
                             ! #End
                             ! print *, '---------Dispersive forcing 2----------'
                             ! print *, pvec(n,jv),pe(n),pvec(n,jv)-pe(n)
                             pSum = pSum + pvec(n,jv) -2.*pv + pv0 ! keep sum
                             ! print *, '++++++Dispersive forcing+++++++++++++++'
                             ! print *, 'check 2nd output P: ', pc,jv,pvec(n,jv)
                           enddo ! jv
                           ! second order update of E
                           evec(n) = (2.*ev-ev0) + cSq1*dtsq*LE1(n)/cSq1 - alphaP1*pSum + dtsq*fev1(n) ! cSq1 is already in LE1
                           ! print *, '++++++Dispersive forcing+++++++++++++++'
                           ! print *, ev,ev0, cSq1, dtsq,LE1(n),alphaP1,pSum,fev1(n)
                           ! print *, 'check 2nd output E: ', ec,evec(n)
                         enddo ! n
                         ! N outside of space loop
                         ! 1st derivative
                         do na=0,numberOfAtomicLevels1-1
                           qt(na) = fnv1(na)
                           do jv = 0,numberOfAtomicLevels1-1
                             qt(na) = qt(na)+prc1(na,jv)*q1(i1,i2,i3,jv)
                           enddo
                           do n=0,nd-1
                             do jv=0,numberOfPolarizationVectors1-1
                               qt(na) = qt(na) + peptc1(na,jv)*u1(i1,i2,i3,ex+n)*(pvec(n,jv)-p1n(i1,i2,i3,n+jv*nd))/(2.*dt)
                             enddo
                           enddo
                         enddo
                         ! 2nd derivative
                         do na=0,numberOfAtomicLevels1-1
                           qtt(na) = fntv1(na)
                           do jv = 0,numberOfAtomicLevels1-1
                             qtt(na) = qtt(na)+prc1(na,jv)*qt(jv)
                           enddo
                           do n=0,nd-1
                             do jv=0,numberOfPolarizationVectors1-1
                               qtt(na) = qtt(na) + peptc1(na,jv)*(evec(n)-u1n(i1,i2,i3,ex+n))/(2.*dt)*(pvec(n,jv)-p1n(i1,i2,i3,n+jv*nd))/(2.*dt)+ peptc1(na,jv)*u1(i1,i2,i3,ex+n)*(pvec(n,jv)-2.*p1(i1,i2,i3,n+jv*nd)+p1n(i1,i2,i3,n+jv*nd))/(dtsq)
                             enddo
                           enddo
                         enddo
                         ! taylor expansion
                         do na=0,numberOfAtomicLevels1-1
                           qv(na) = q1(i1,i2,i3,na) + dt*qt(na) + dtsq/2.*qtt(na)
                         enddo
                         !----------------------------------------
                         ! order 4 update of P at interface (fictitious step)
                         !----------------------------------------
                         ! second order accurate terms
                         do n=0,nd-1
                           do jv = 0,numberOfPolarizationVectors1-1
                             ! #If "p1" eq "p1"
                             ! call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pe(n)   )
                             ! call ogderiv(ep, 1,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pet(n)   )
                             ! call ogderiv(ep, 2,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pett(n)   )
                             ! #Else
                             ! call ogderiv(ep, 0,0,0,0, xy2(i1,i2,i3,0),xy2(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pe(n)   )
                             ! call ogderiv(ep, 1,0,0,0, xy2(i1,i2,i3,0),xy2(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pet(n)   )
                             ! call ogderiv(ep, 2,0,0,0, xy2(i1,i2,i3,0),xy2(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pett(n)   )
                             ! #End
                             ! pvec(n,jv) = pe(n)
                             ! ptv(n,jv) = pet(n)
                             ! pttv(n,jv) = pett(n)
                             ptv(n,jv) = (pvec(n,jv)-p1n(i1,i2,i3,n+jv*nd))/(2.*dt)
                             pttv(n,jv) = (pvec(n,jv)-2.*p1(i1,i2,i3,n+jv*nd)+p1n(i1,i2,i3,n+jv*nd))/dtsq
                             ptttv(n,jv) = -b1v1(jv)*pttv(n,jv)-b0v1(jv)*ptv(n,jv)+fPt1(n,jv)
                             do na = 0,numberOfAtomicLevels1-1 ! update using ODE
                               ptttv(n,jv) = ptttv(n,jv) + pnec1(jv,na)*qt(na)*u1(i1,i2,i3,ex+n) + pnec1(jv,na)*q1(i1,i2,i3,na)*(evec(n)-u1n(i1,i2,i3,ex+n))/(2.*dt)
                             enddo
                             pttttv(n,jv) = -b1v1(jv)*ptttv(n,jv)-b0v1(jv)*pttv(n,jv)+fPtt1(n,jv) ! update using ODE
                             do na = 0,numberOfAtomicLevels1-1
                               pttttv(n,jv) = pttttv(n,jv) + pnec1(jv,na)*qtt(na)*u1(i1,i2,i3,ex+n) + 2.*pnec1(jv,na)*qt(na)*(evec(n)-u1n(i1,i2,i3,ex+n))/(2.*dt) + pnec1(jv,na)*q1(i1,i2,i3,na)*(evec(n)-2.*u1(i1,i2,i3,ex+n)+u1n(i1,i2,i3,ex+n))/dtsq
                               ! print *, '++++++Dispersive forcing+++++++++++++++'
                               ! print *, 'check P derivatives: ', n,jv,ptv(n,jv),pttv(n,jv),ptttv(n,jv),pttttv(n,jv)
                               ! print *, na,pnec1(jv,na),q1(i1,i2,i3,na),qt(na),qtt(na)
                             enddo
                             ! print *, '++++++Dispersive forcing+++++++++++++++'
                             ! print *, 'check P time derivatives: ', n,jv,ptv(n,jv),pttv(n,jv),ptttv(n,jv),pttttv(n,jv)
                             ! print *, fPt1(n,jv),fPtt1(n,jv)
                           enddo
                           ! print *, evec(n),u1(i1,i2,i3,ex+n),u1n(i1,i2,i3,ex+n)
                         enddo
                         ! dimension loop for E and P
                         do n=0,nd-1
                           ec = ex +n 
                           ev   =  u1(i1,i2,i3,ec) ! time where we need to fill in ghost points
                           ! These next derivatives may only be needed to order2, (use order 4 for testing TZ polynomials)
                           ! #perl $ORDER=2; ! should use order 2 since E is only filled at first ghost lines at t
                             ! uu1 in the next statement defines names of intermediate values
                               uu1=u1(i1,i2,i3,ec) ! in the rectangular case just eval the solution
                                e1x = (u1(i1-2,i2,i3,ec)-8.*u1(i1-1,i2,i3,ec)+8.*u1(i1+1,i2,i3,ec)-u1(i1+2,i2,i3,ec))/(12.*dx1(0))
                                e1y = (u1(i1,i2-2,i3,ec)-8.*u1(i1,i2-1,i3,ec)+8.*u1(i1,i2+1,i3,ec)-u1(i1,i2+2,i3,ec))/(12.*dx1(1))
                                e1xx = (-u1(i1-2,i2,i3,ec)+16.*u1(i1-1,i2,i3,ec)-30.*u1(i1,i2,i3,ec)+16.*u1(i1+1,i2,i3,ec)-u1(i1+2,i2,i3,ec))/(12.*dx1(0)**2)
                                e1yy = (-u1(i1,i2-2,i3,ec)+16.*u1(i1,i2-1,i3,ec)-30.*u1(i1,i2,i3,ec)+16.*u1(i1,i2+1,i3,ec)-u1(i1,i2+2,i3,ec))/(12.*dx1(1)**2)
                              e1Lap = e1xx+ e1yy
                             ! write(*,'("LEFT: u1x,u1y,u1xx,u1yy,u1Lap=",5(1pe12.4))') u1x,u1y,u1xx,u1yy,u1Lap
                             evx0  = e1x
                             evy0  = e1y
                             evLap = e1xx+e1yy ! these values use the second order predicted values in the first ghost lines
                               ! print *, '++++++Dispersive forcing+++++++++++++++'
                               ! print *, 'check E spatial derivatives: ', n,evx0,evy0,evLap
                           do jv=0,numberOfPolarizationVectors1-1
                             pc = n + jv*nd 
                             ! These next derivatives may only be needed to order2, (use order 4 for testing TZ polynomials)
                             ! #perl $ORDER=2;
                               ! uu1 in the next statement defines names of intermediate values
                                 pp1=p1(i1,i2,i3,pc) ! in the rectangular case just eval the solution
                                  p1x = (p1(i1-2,i2,i3,pc)-8.*p1(i1-1,i2,i3,pc)+8.*p1(i1+1,i2,i3,pc)-p1(i1+2,i2,i3,pc))/(12.*dx1(0))
                                  p1y = (p1(i1,i2-2,i3,pc)-8.*p1(i1,i2-1,i3,pc)+8.*p1(i1,i2+1,i3,pc)-p1(i1,i2+2,i3,pc))/(12.*dx1(1))
                                  p1xx = (-p1(i1-2,i2,i3,pc)+16.*p1(i1-1,i2,i3,pc)-30.*p1(i1,i2,i3,pc)+16.*p1(i1+1,i2,i3,pc)-p1(i1+2,i2,i3,pc))/(12.*dx1(0)**2)
                                  p1yy = (-p1(i1,i2-2,i3,pc)+16.*p1(i1,i2-1,i3,pc)-30.*p1(i1,i2,i3,pc)+16.*p1(i1,i2+1,i3,pc)-p1(i1,i2+2,i3,pc))/(12.*dx1(1)**2)
                                p1Lap = p1xx+ p1yy
                               ! write(*,'("LEFT: p1x,p1y,p1xx,p1yy,p1Lap=",5(1pe12.4))') p1x,p1y,p1xx,p1yy,p1Lap
                               LP  = p1Lap ! removed c^2
                                 pp1=p1n(i1,i2,i3,pc) ! in the rectangular case just eval the solution
                                  p1nx = (p1n(i1-2,i2,i3,pc)-8.*p1n(i1-1,i2,i3,pc)+8.*p1n(i1+1,i2,i3,pc)-p1n(i1+2,i2,i3,pc))/(12.*dx1(0))
                                  p1ny = (p1n(i1,i2-2,i3,pc)-8.*p1n(i1,i2-1,i3,pc)+8.*p1n(i1,i2+1,i3,pc)-p1n(i1,i2+2,i3,pc))/(12.*dx1(1))
                                  p1nxx = (-p1n(i1-2,i2,i3,pc)+16.*p1n(i1-1,i2,i3,pc)-30.*p1n(i1,i2,i3,pc)+16.*p1n(i1+1,i2,i3,pc)-p1n(i1+2,i2,i3,pc))/(12.*dx1(0)**2)
                                  p1nyy = (-p1n(i1,i2-2,i3,pc)+16.*p1n(i1,i2-1,i3,pc)-30.*p1n(i1,i2,i3,pc)+16.*p1n(i1,i2+1,i3,pc)-p1n(i1,i2+2,i3,pc))/(12.*dx1(1)**2)
                                p1nLap = p1nxx+ p1nyy
                               ! write(*,'("LEFT: p1nxx,p1nyy,p1nLap=",3e12.4)') p1nxx,p1nyy,p1nLap
                               LPm = p1nLap
                               pvx  = p1x
                               pvy  = p1y
                               pvnx = p1nx
                               pvny = p1ny
                             ! print *, '++++++Dispersive forcing+++++++++++++++'
                             ! print *, 'check P spatial derivatives: ', n,jv,LP,LPm,pvx,pvy,pvnx,pvny
                             pv0 =  p1n(i1,i2,i3,pc)
                             pv  =  p1(i1,i2,i3,pc)
                             pvn = 2.*pv-p1n(i1,i2,i3,pc) + 0.5*dt*b1v1(jv)*p1n(i1,i2,i3,pc) + dt**4/12.*pttttv(n,jv) + dt**4/6.*b1v1(jv)*ptttv(n,jv) - dtsq*b0v1(jv)*pv + dtsq*fpv1(n,jv)
                             ! pvn = 2.*pv-p1n(i1,i2,i3,pc) + 0.5*dt*b1v1(jv)*p1n(i1,i2,i3,pc) - dtsq*b0v1(jv)*pv + dtsq*fpv1(n,jv)
                             do na = 0,numberOfAtomicLevels1-1 ! \Delta N^n*E^n
                               pvn = pvn + dtsq*pnec1(jv,na)*q1(i1,i2,i3,na)*ev
                             enddo ! na
                             pvec(n,jv)= pvn/( 1.+.5*dt*b1v1(jv) ) ! time + dt
                             ! 4th order accurate term
                             fp1(n) = fp1(n) + (pvec(n,jv)-2.*p1(i1,i2,i3,pc)+p1n(i1,i2,i3,pc))/dtsq  - dt**2/12.*pttttv(n,jv)
                             ! #If "p1" eq "p1"
                             !   call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t+dt,pxc+jv*nd+n, pe(n)   )
                             !   call ogderiv(ep, 1,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pet(n)   )
                             !   call ogderiv(ep, 2,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pett(n)   )
                             !   call ogderiv(ep, 4,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pxc+jv*nd+n, petttt(n)   )
                             ! #Else
                             !   call ogderiv(ep, 0,0,0,0, xy2(i1,i2,i3,0),xy2(i1,i2,i3,1),0.,t+dt,pxc+jv*nd+n, pe(n)   )
                             !   call ogderiv(ep, 1,0,0,0, xy2(i1,i2,i3,0),xy2(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pet(n)   )
                             !   call ogderiv(ep, 2,0,0,0, xy2(i1,i2,i3,0),xy2(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pett(n)   )
                             !   call ogderiv(ep, 4,0,0,0, xy2(i1,i2,i3,0),xy2(i1,i2,i3,1),0.,t,pxc+jv*nd+n, petttt(n)   )
                             ! #End
                             ! print *, '---------Dispersive forcing 4----------'
                             ! print *, pvec(n,jv),pe(n),pvec(n,jv)-pe(n)
                             ! print *, (pvec(n,jv)-2.*p1(i1,i2,i3,pc)+p1n(i1,i2,i3,pc))/dtsq  - dt**2/12.*pttttv(n,jv), pett(n),(pvec(n,jv)-2.*p1(i1,i2,i3,pc)+p1n(i1,i2,i3,pc))/dtsq  - dt**2/12.*pttttv(n,jv)-pett(n)
                             ! fp1(n) = fp1(n) + pett(n)
                             ! 2nd order accurate terms
                             fPtttt1(n) = fPtttt1(n) + pttttv(n,jv)
                             ! fPtttt1(n) = fPtttt1(n) + petttt(n)
                             ! print *, '++++++Dispersive forcing+++++++++++++++'
                             ! print *, 'check P time derivatives 2 and 4: ', n,jv,(pvec(n,jv)-2.*p1(i1,i2,i3,pc)+p1n(i1,i2,i3,pc))/dtsq  - dt**4/12.*pttttv(n,jv),pttttv(n,jv)
                             !--------------------------------
                             ! spatial derivatives
                             !--------------------------------
                             ! nce = pxc+nd*numberOfPolarizationVectors1
                             ! N*E
                             ! #perl $ORDER=2;
                             do na=0,numberOfAtomicLevels1-1
                               ! uu1 in the next statement defines names of intermediate values
                                 qq1=q1(i1,i2,i3,na) ! in the rectangular case just eval the solution
                                  q1x = (q1(i1-2,i2,i3,na)-8.*q1(i1-1,i2,i3,na)+8.*q1(i1+1,i2,i3,na)-q1(i1+2,i2,i3,na))/(12.*dx1(0))
                                  q1y = (q1(i1,i2-2,i3,na)-8.*q1(i1,i2-1,i3,na)+8.*q1(i1,i2+1,i3,na)-q1(i1,i2+2,i3,na))/(12.*dx1(1))
                                  q1xx = (-q1(i1-2,i2,i3,na)+16.*q1(i1-1,i2,i3,na)-30.*q1(i1,i2,i3,na)+16.*q1(i1+1,i2,i3,na)-q1(i1+2,i2,i3,na))/(12.*dx1(0)**2)
                                  q1yy = (-q1(i1,i2-2,i3,na)+16.*q1(i1,i2-1,i3,na)-30.*q1(i1,i2,i3,na)+16.*q1(i1,i2+1,i3,na)-q1(i1,i2+2,i3,na))/(12.*dx1(1)**2)
                                q1Lap = q1xx+ q1yy
                               ! write(*,'("LEFT: p1x,p1y,p1xx,p1yy,p1Lap=",5(1pe12.4))') p1x,p1y,p1xx,p1yy,p1Lap
                               qvx  = q1x
                               qvy  = q1y
                               qvLap  = q1Lap
                               qex(na) = evx0*q1(i1,i2,i3,na)+qvx*ev
                               qey(na) = evy0*q1(i1,i2,i3,na)+qvy*ev
                               qeLap(na) = ev*qvLap+q1(i1,i2,i3,na)*evLap+2.*evx0*qvx+2.*evy0*qvy
                               ! print *, '++++++Dispersive forcing+++++++++++++++'
                               ! print *, 'check N spatial derivatives: ', na,qex(na),qey(na),qeLap(na)
                             enddo
                             ! #If "p1" eq "p1"
                             !   call ogderiv(ep, 2,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pevttxx(n,jv)   )
                             !   call ogderiv(ep, 2,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pevttyy(n,jv)   )
                             !   call ogderiv(ep, 2,1,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pevttx1(n,jv)   )
                             !   call ogderiv(ep, 2,0,1,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pevtty1(n,jv)   )
                             ! #Else
                             !   call ogderiv(ep, 2,2,0,0, xy2(i1,i2,i3,0),xy2(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pevttxx(n,jv)   )
                             !   call ogderiv(ep, 2,0,2,0, xy2(i1,i2,i3,0),xy2(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pevttyy(n,jv)   )
                             !   call ogderiv(ep, 2,1,0,0, xy2(i1,i2,i3,0),xy2(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pevttx1(n,jv)   )
                             !   call ogderiv(ep, 2,0,1,0, xy2(i1,i2,i3,0),xy2(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pevtty1(n,jv)   )
                             ! #End
                             ! laplacian
                             LPn = 2.*LP-LPm + 0.5*dt*b1v1(jv)*LPm - dtsq*b0v1(jv)*LP + dtsq*LfP1(n,jv)
                             do na = 0,numberOfAtomicLevels1-1 ! \Delta N^n*E^n
                               LPn = LPn + dtsq*pnec1(jv,na)*qeLap(na)
                             enddo
                             ! time derivatives
                             fLPtt1(n) = fLPtt1(n) + (LPn/(1.+.5*dt*b1v1(jv)) - 2.*LP + LPm)/dtsq
                             ! fLPtt1(n) = fLPtt1(n) + pevttL1(n,jv)
                             ! print *,'error of laplacian with # of levels',numberOfAtomicLevels1,numberOfPolarizationVectors1,pevttxx(n,jv)+pevttyy(n,jv)-pevttL1(n,jv)
                             ! print *, '++++++Dispersive forcing+++++++++++++++'
                             ! print *, 'check P lap: ', n,jv,LP,LPm,b1v1(jv),dtsq,b0v1(jv),LfP1(n,jv),fLPtt1(n)
                             ! x
                             pttxa = 2.*pvx-pvnx + 0.5*dt*b1v1(jv)*pvnx - dtsq*b0v1(jv)*pvx + dtsq*fpvx1(n,jv)
                             do na = 0,numberOfAtomicLevels1-1 ! \Delta N^n*E^n
                               pttxa = pttxa + dtsq*pnec1(jv,na)*qex(na)
                             enddo
                             ! time derivatives
                             fPttx1(n) = fPttx1(n) + (pttxa/(1.+.5*dt*b1v1(jv)) - 2.*pvx + pvnx)/dtsq
                             ! fPttx1(n) = fPttx1(n) + pevttx1(n,jv)
                             ! print *, '++++++Dispersive forcing+++++++++++++++'
                             ! print *, 'check Pttx: ', n,jv,pvx,pvnx,fpvx1(n,jv),fPttx1(n)
                             ! y
                             pttya = 2.*pvy-pvny + 0.5*dt*b1v1(jv)*pvny - dtsq*b0v1(jv)*pvy + dtsq*fpvy1(n,jv)
                             do na = 0,numberOfAtomicLevels1-1 ! \Delta N^n*E^n
                               pttya = pttya + dtsq*pnec1(jv,na)*qey(na)
                             enddo
                             ! time derivatives
                             fPtty1(n) = fPtty1(n) + (pttya/(1.+.5*dt*b1v1(jv)) - 2.*pvy + pvny)/dtsq
                             ! fPtty1(n) = fPtty1(n) + pevtty1(n,jv)
                             ! print *, '++++++Dispersive forcing+++++++++++++++'
                             ! print *, 'check Ptty: ', n,jv,n,jv,pvy,pvny,fpvy1(n,jv),fPtty1(n)
                             ! if( twilightZone.eq.1 )then
                             !   write(*,'("")')
                             !   write(*,'("DI4:LEFT: i1,i2=",2i3," jv=",i2," n=",i2)') i1,i2,jv,n
                             !   print *, 'ptt diff',(pvec(n,jv)-2.*p1(i1,i2,i3,pc)+p1n(i1,i2,i3,pc))/dtsq - dt**2/12.*pttttv(n,jv)-pevtt1(n,jv)
                             !   print *, 'pttx diff', (pttxa/(1.+.5*dt*b1v1(jv)) - 2.*pvx + pvnx)/dtsq-pevttx1(n,jv)
                             !   print *, 'ptty diff', (pttya/(1.+.5*dt*b1v1(jv)) - 2.*pvy + pvny)/dtsq-pevtty1(n,jv)
                             !   print *, 'ptttt diff', pttttv(n,jv)-pevtttt1(n,jv)
                             !   print *, 'pttL diff',(LPn/(1.+.5*dt*b1v1(jv)) - 2.*LP + LPm)/dtsq-pevttL1(n,jv)
                             ! end if
                           enddo ! jv
                         enddo ! n
                       end if
                     ! eval nonlinear dispersive forcings for domain 2
                       ! pre-assign 0 values
                       do n=0,nd-1 ! dispersive forcing in jump conditions
                         fp2(n)=0.
                         fPttx2(n) =0.
                         fPtty2(n) =0.
                         fLPtt2(n) =0.
                         fPtttt2(n)=0.
                       end do
                       ! only do this for MLA (dispersive and nonlinear multi-level)
                       if( dispersionModel2.ne.noDispersion .and. nonlinearModel2.ne.noNonlinearModel) then
                         nce = pxc+nd*numberOfPolarizationVectors2
                         ! -----------------------------------------
                         ! order 2 (E, P, N) at the interface (fictitious step)
                         !------------------------------------------
                         ! dimension loop for E and P
                         do n=0,nd-1
                           ec = ex +n 
                           ev0  =  u2n(j1,j2,j3,ec)
                           ev   =  u2(j1,j2,j3,ec) ! time where we need to fill in ghost points
                           pSum=0.
                           do jv=0,numberOfPolarizationVectors2-1
                             pc = n + jv*nd  
                             pv0 =  p2n(j1,j2,j3,pc)
                             pv  =  p2(j1,j2,j3,pc)
                             pvn = 2.*pv-p2n(j1,j2,j3,pc) + 0.5*dt*b1v2(jv)*p2n(j1,j2,j3,pc) - dtsq*b0v2(jv)*pv + dtsq*fpv2(n,jv)
                             do na = 0,numberOfAtomicLevels2-1 ! \Delta N^n*E^n
                               pvn = pvn + dtsq*pnec2(jv,na)*q2(j1,j2,j3,na)*ev
                             enddo ! na
                             pvec(n,jv)= pvn/( 1.+.5*dt*b1v2(jv) ) ! time + dt
                             ! #If "p2" eq "p1"
                             ! call ogderiv(ep, 0,0,0,0, xy1(j1,j2,j3,0),xy1(j1,j2,j3,1),0.,t+dt,pxc+jv*nd+n, pe(n)   )
                             ! ! call ogderiv(ep, 1,0,0,0, xy1(j1,j2,j3,0),xy1(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pet(n)   )
                             ! ! call ogderiv(ep, 2,0,0,0, xy1(j1,j2,j3,0),xy1(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pett(n)   )
                             ! #Else
                             ! call ogderiv(ep, 0,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t+dt,pxc+jv*nd+n, pe(n)   )
                             ! ! call ogderiv(ep, 1,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pet(n)   )
                             ! ! call ogderiv(ep, 2,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pett(n)   )
                             ! #End
                             ! print *, '---------Dispersive forcing 2----------'
                             ! print *, pvec(n,jv),pe(n),pvec(n,jv)-pe(n)
                             pSum = pSum + pvec(n,jv) -2.*pv + pv0 ! keep sum
                             ! print *, '++++++Dispersive forcing+++++++++++++++'
                             ! print *, 'check 2nd output P: ', pc,jv,pvec(n,jv)
                           enddo ! jv
                           ! second order update of E
                           evec(n) = (2.*ev-ev0) + cSq2*dtsq*LE2(n)/cSq2 - alphaP2*pSum + dtsq*fev2(n) ! cSq2 is already in LE2
                           ! print *, '++++++Dispersive forcing+++++++++++++++'
                           ! print *, ev,ev0, cSq2, dtsq,LE2(n),alphaP2,pSum,fev2(n)
                           ! print *, 'check 2nd output E: ', ec,evec(n)
                         enddo ! n
                         ! N outside of space loop
                         ! 1st derivative
                         do na=0,numberOfAtomicLevels2-1
                           qt(na) = fnv2(na)
                           do jv = 0,numberOfAtomicLevels2-1
                             qt(na) = qt(na)+prc2(na,jv)*q2(j1,j2,j3,jv)
                           enddo
                           do n=0,nd-1
                             do jv=0,numberOfPolarizationVectors2-1
                               qt(na) = qt(na) + peptc2(na,jv)*u2(j1,j2,j3,ex+n)*(pvec(n,jv)-p2n(j1,j2,j3,n+jv*nd))/(2.*dt)
                             enddo
                           enddo
                         enddo
                         ! 2nd derivative
                         do na=0,numberOfAtomicLevels2-1
                           qtt(na) = fntv2(na)
                           do jv = 0,numberOfAtomicLevels2-1
                             qtt(na) = qtt(na)+prc2(na,jv)*qt(jv)
                           enddo
                           do n=0,nd-1
                             do jv=0,numberOfPolarizationVectors2-1
                               qtt(na) = qtt(na) + peptc2(na,jv)*(evec(n)-u2n(j1,j2,j3,ex+n))/(2.*dt)*(pvec(n,jv)-p2n(j1,j2,j3,n+jv*nd))/(2.*dt)+ peptc2(na,jv)*u2(j1,j2,j3,ex+n)*(pvec(n,jv)-2.*p2(j1,j2,j3,n+jv*nd)+p2n(j1,j2,j3,n+jv*nd))/(dtsq)
                             enddo
                           enddo
                         enddo
                         ! taylor expansion
                         do na=0,numberOfAtomicLevels2-1
                           qv(na) = q2(j1,j2,j3,na) + dt*qt(na) + dtsq/2.*qtt(na)
                         enddo
                         !----------------------------------------
                         ! order 4 update of P at interface (fictitious step)
                         !----------------------------------------
                         ! second order accurate terms
                         do n=0,nd-1
                           do jv = 0,numberOfPolarizationVectors2-1
                             ! #If "p2" eq "p1"
                             ! call ogderiv(ep, 0,0,0,0, xy1(j1,j2,j3,0),xy1(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pe(n)   )
                             ! call ogderiv(ep, 1,0,0,0, xy1(j1,j2,j3,0),xy1(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pet(n)   )
                             ! call ogderiv(ep, 2,0,0,0, xy1(j1,j2,j3,0),xy1(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pett(n)   )
                             ! #Else
                             ! call ogderiv(ep, 0,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pe(n)   )
                             ! call ogderiv(ep, 1,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pet(n)   )
                             ! call ogderiv(ep, 2,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pett(n)   )
                             ! #End
                             ! pvec(n,jv) = pe(n)
                             ! ptv(n,jv) = pet(n)
                             ! pttv(n,jv) = pett(n)
                             ptv(n,jv) = (pvec(n,jv)-p2n(j1,j2,j3,n+jv*nd))/(2.*dt)
                             pttv(n,jv) = (pvec(n,jv)-2.*p2(j1,j2,j3,n+jv*nd)+p2n(j1,j2,j3,n+jv*nd))/dtsq
                             ptttv(n,jv) = -b1v2(jv)*pttv(n,jv)-b0v2(jv)*ptv(n,jv)+fPt2(n,jv)
                             do na = 0,numberOfAtomicLevels2-1 ! update using ODE
                               ptttv(n,jv) = ptttv(n,jv) + pnec2(jv,na)*qt(na)*u2(j1,j2,j3,ex+n) + pnec2(jv,na)*q2(j1,j2,j3,na)*(evec(n)-u2n(j1,j2,j3,ex+n))/(2.*dt)
                             enddo
                             pttttv(n,jv) = -b1v2(jv)*ptttv(n,jv)-b0v2(jv)*pttv(n,jv)+fPtt2(n,jv) ! update using ODE
                             do na = 0,numberOfAtomicLevels2-1
                               pttttv(n,jv) = pttttv(n,jv) + pnec2(jv,na)*qtt(na)*u2(j1,j2,j3,ex+n) + 2.*pnec2(jv,na)*qt(na)*(evec(n)-u2n(j1,j2,j3,ex+n))/(2.*dt) + pnec2(jv,na)*q2(j1,j2,j3,na)*(evec(n)-2.*u2(j1,j2,j3,ex+n)+u2n(j1,j2,j3,ex+n))/dtsq
                               ! print *, '++++++Dispersive forcing+++++++++++++++'
                               ! print *, 'check P derivatives: ', n,jv,ptv(n,jv),pttv(n,jv),ptttv(n,jv),pttttv(n,jv)
                               ! print *, na,pnec2(jv,na),q2(j1,j2,j3,na),qt(na),qtt(na)
                             enddo
                             ! print *, '++++++Dispersive forcing+++++++++++++++'
                             ! print *, 'check P time derivatives: ', n,jv,ptv(n,jv),pttv(n,jv),ptttv(n,jv),pttttv(n,jv)
                             ! print *, fPt2(n,jv),fPtt2(n,jv)
                           enddo
                           ! print *, evec(n),u2(j1,j2,j3,ex+n),u2n(j1,j2,j3,ex+n)
                         enddo
                         ! dimension loop for E and P
                         do n=0,nd-1
                           ec = ex +n 
                           ev   =  u2(j1,j2,j3,ec) ! time where we need to fill in ghost points
                           ! These next derivatives may only be needed to order2, (use order 4 for testing TZ polynomials)
                           ! #perl $ORDER=2; ! should use order 2 since E is only filled at first ghost lines at t
                             ! uu2 in the next statement defines names of intermediate values
                               uu2=u2(j1,j2,j3,ec) ! in the rectangular case just eval the solution
                                e2x = (u2(j1-2,j2,j3,ec)-8.*u2(j1-1,j2,j3,ec)+8.*u2(j1+1,j2,j3,ec)-u2(j1+2,j2,j3,ec))/(12.*dx2(0))
                                e2y = (u2(j1,j2-2,j3,ec)-8.*u2(j1,j2-1,j3,ec)+8.*u2(j1,j2+1,j3,ec)-u2(j1,j2+2,j3,ec))/(12.*dx2(1))
                                e2xx = (-u2(j1-2,j2,j3,ec)+16.*u2(j1-1,j2,j3,ec)-30.*u2(j1,j2,j3,ec)+16.*u2(j1+1,j2,j3,ec)-u2(j1+2,j2,j3,ec))/(12.*dx2(0)**2)
                                e2yy = (-u2(j1,j2-2,j3,ec)+16.*u2(j1,j2-1,j3,ec)-30.*u2(j1,j2,j3,ec)+16.*u2(j1,j2+1,j3,ec)-u2(j1,j2+2,j3,ec))/(12.*dx2(1)**2)
                              e2Lap = e2xx+ e2yy
                             ! write(*,'("RIGHT: u2x,u2y,u2xx,u2yy,u2Lap=",5(1pe12.4))') u2x,u2y,u2xx,u2yy,u2Lap
                             evx0  = e2x
                             evy0  = e2y
                             evLap = e2xx+e2yy
                               ! print *, '++++++Dispersive forcing+++++++++++++++'
                               ! print *, 'check E spatial derivatives: ', n,evx0,evy0,evLap
                           do jv=0,numberOfPolarizationVectors2-1
                             pc = n + jv*nd 
                             ! These next derivatives may only be needed to order2, (use order 4 for testing TZ polynomials)
                             ! #perl $ORDER=2;
                               ! uu1 in the next statement defines names of intermediate values
                                 pp2=p2(j1,j2,j3,pc) ! in the rectangular case just eval the solution
                                  p2x = (p2(j1-2,j2,j3,pc)-8.*p2(j1-1,j2,j3,pc)+8.*p2(j1+1,j2,j3,pc)-p2(j1+2,j2,j3,pc))/(12.*dx1(0))
                                  p2y = (p2(j1,j2-2,j3,pc)-8.*p2(j1,j2-1,j3,pc)+8.*p2(j1,j2+1,j3,pc)-p2(j1,j2+2,j3,pc))/(12.*dx1(1))
                                  p2xx = (-p2(j1-2,j2,j3,pc)+16.*p2(j1-1,j2,j3,pc)-30.*p2(j1,j2,j3,pc)+16.*p2(j1+1,j2,j3,pc)-p2(j1+2,j2,j3,pc))/(12.*dx1(0)**2)
                                  p2yy = (-p2(j1,j2-2,j3,pc)+16.*p2(j1,j2-1,j3,pc)-30.*p2(j1,j2,j3,pc)+16.*p2(j1,j2+1,j3,pc)-p2(j1,j2+2,j3,pc))/(12.*dx1(1)**2)
                                p2Lap = p2xx+ p2yy
                               LP  = p2Lap
                                 pp2=p2n(j1,j2,j3,pc) ! in the rectangular case just eval the solution
                                  p2nx = (p2n(j1-2,j2,j3,pc)-8.*p2n(j1-1,j2,j3,pc)+8.*p2n(j1+1,j2,j3,pc)-p2n(j1+2,j2,j3,pc))/(12.*dx1(0))
                                  p2ny = (p2n(j1,j2-2,j3,pc)-8.*p2n(j1,j2-1,j3,pc)+8.*p2n(j1,j2+1,j3,pc)-p2n(j1,j2+2,j3,pc))/(12.*dx1(1))
                                  p2nxx = (-p2n(j1-2,j2,j3,pc)+16.*p2n(j1-1,j2,j3,pc)-30.*p2n(j1,j2,j3,pc)+16.*p2n(j1+1,j2,j3,pc)-p2n(j1+2,j2,j3,pc))/(12.*dx1(0)**2)
                                  p2nyy = (-p2n(j1,j2-2,j3,pc)+16.*p2n(j1,j2-1,j3,pc)-30.*p2n(j1,j2,j3,pc)+16.*p2n(j1,j2+1,j3,pc)-p2n(j1,j2+2,j3,pc))/(12.*dx1(1)**2)
                                p2nLap = p2nxx+ p2nyy
                               LPm = p2nLap
                               pvx  = p2x
                               pvy  = p2y
                               pvnx = p2nx
                               pvny = p2ny
                             ! print *, '++++++Dispersive forcing+++++++++++++++'
                             ! print *, 'check P spatial derivatives: ', n,jv,LP,LPm,pvx,pvy,pvnx,pvny
                             pv0 =  p2n(j1,j2,j3,pc)
                             pv  =  p2(j1,j2,j3,pc)
                             pvn = 2.*pv-p2n(j1,j2,j3,pc) + 0.5*dt*b1v2(jv)*p2n(j1,j2,j3,pc) + dt**4/12.*pttttv(n,jv) + dt**4/6.*b1v2(jv)*ptttv(n,jv) - dtsq*b0v2(jv)*pv + dtsq*fpv2(n,jv)
                             ! pvn = 2.*pv-p2n(j1,j2,j3,pc) + 0.5*dt*b1v2(jv)*p2n(j1,j2,j3,pc) - dtsq*b0v2(jv)*pv + dtsq*fpv2(n,jv)
                             do na = 0,numberOfAtomicLevels2-1 ! \Delta N^n*E^n
                               pvn = pvn + dtsq*pnec2(jv,na)*q2(j1,j2,j3,na)*ev
                             enddo ! na
                             pvec(n,jv)= pvn/( 1.+.5*dt*b1v2(jv) ) ! time + dt
                             ! 4th order accurate term
                             fp2(n) = fp2(n) + (pvec(n,jv)-2.*p2(j1,j2,j3,pc)+p2n(j1,j2,j3,pc))/dtsq  - dt**2/12.*pttttv(n,jv)
                             ! #If "p2" eq "p1"
                             !   call ogderiv(ep, 0,0,0,0, xy1(j1,j2,j3,0),xy1(j1,j2,j3,1),0.,t+dt,pxc+jv*nd+n, pe(n)   )
                             !   call ogderiv(ep, 1,0,0,0, xy1(j1,j2,j3,0),xy1(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pet(n)   )
                             !   call ogderiv(ep, 2,0,0,0, xy1(j1,j2,j3,0),xy1(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pett(n)   )
                             !   call ogderiv(ep, 4,0,0,0, xy1(j1,j2,j3,0),xy1(j1,j2,j3,1),0.,t,pxc+jv*nd+n, petttt(n)   )
                             ! #Else
                             !   call ogderiv(ep, 0,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t+dt,pxc+jv*nd+n, pe(n)   )
                             !   call ogderiv(ep, 1,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pet(n)   )
                             !   call ogderiv(ep, 2,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pett(n)   )
                             !   call ogderiv(ep, 4,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pxc+jv*nd+n, petttt(n)   )
                             ! #End
                             ! print *, '---------Dispersive forcing 4----------'
                             ! print *, pvec(n,jv),pe(n),pvec(n,jv)-pe(n)
                             ! print *, (pvec(n,jv)-2.*p2(j1,j2,j3,pc)+p2n(j1,j2,j3,pc))/dtsq  - dt**2/12.*pttttv(n,jv), pett(n),(pvec(n,jv)-2.*p2(j1,j2,j3,pc)+p2n(j1,j2,j3,pc))/dtsq  - dt**2/12.*pttttv(n,jv)-pett(n)
                             ! fp2(n) = fp2(n) + pett(n)
                             ! 2nd order accurate terms
                             fPtttt2(n) = fPtttt2(n) + pttttv(n,jv)
                             ! fPtttt2(n) = fPtttt2(n) + petttt(n)
                             ! print *, '++++++Dispersive forcing+++++++++++++++'
                             ! print *, 'check P time derivatives 2 and 4: ', n,jv,(pvec(n,jv)-2.*p2(j1,j2,j3,pc)+p2n(j1,j2,j3,pc))/dtsq  - dt**4/12.*pttttv(n,jv),pttttv(n,jv)
                             !--------------------------------
                             ! spatial derivatives
                             !--------------------------------
                             ! nce = pxc+nd*numberOfPolarizationVectors2
                             ! N*E
                             ! #perl $ORDER=2;
                             do na=0,numberOfAtomicLevels2-1
                               ! uu2 in the next statement defines names of intermediate values
                                 qq2=q2(j1,j2,j3,na) ! in the rectangular case just eval the solution
                                  q2x = (q2(j1-2,j2,j3,na)-8.*q2(j1-1,j2,j3,na)+8.*q2(j1+1,j2,j3,na)-q2(j1+2,j2,j3,na))/(12.*dx2(0))
                                  q2y = (q2(j1,j2-2,j3,na)-8.*q2(j1,j2-1,j3,na)+8.*q2(j1,j2+1,j3,na)-q2(j1,j2+2,j3,na))/(12.*dx2(1))
                                  q2xx = (-q2(j1-2,j2,j3,na)+16.*q2(j1-1,j2,j3,na)-30.*q2(j1,j2,j3,na)+16.*q2(j1+1,j2,j3,na)-q2(j1+2,j2,j3,na))/(12.*dx2(0)**2)
                                  q2yy = (-q2(j1,j2-2,j3,na)+16.*q2(j1,j2-1,j3,na)-30.*q2(j1,j2,j3,na)+16.*q2(j1,j2+1,j3,na)-q2(j1,j2+2,j3,na))/(12.*dx2(1)**2)
                                q2Lap = q2xx+ q2yy
                               qvx  = q2x
                               qvy  = q2y
                               qvLap  = q2Lap
                               qex(na) = evx0*q2(j1,j2,j3,na)+qvx*ev
                               qey(na) = evy0*q2(j1,j2,j3,na)+qvy*ev
                               qeLap(na) = ev*qvLap+q2(j1,j2,j3,na)*evLap+2.*evx0*qvx+2.*evy0*qvy
                               ! print *, '++++++Dispersive forcing+++++++++++++++'
                               ! print *, 'check N spatial derivatives: ', na,qex(na),qey(na),qeLap(na)
                             enddo
                             ! #If "p2" eq "p1"
                             !   call ogderiv(ep, 2,2,0,0, xy1(j1,j2,j3,0),xy1(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pevttxx(n,jv)   )
                             !   call ogderiv(ep, 2,0,2,0, xy1(j1,j2,j3,0),xy1(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pevttyy(n,jv)   )
                             !   call ogderiv(ep, 2,1,0,0, xy1(j1,j2,j3,0),xy1(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pevttx2(n,jv)   )
                             !   call ogderiv(ep, 2,0,1,0, xy1(j1,j2,j3,0),xy1(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pevtty2(n,jv)   )
                             ! #Else
                             !   call ogderiv(ep, 2,2,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pevttxx(n,jv)   )
                             !   call ogderiv(ep, 2,0,2,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pevttyy(n,jv)   )
                             !   call ogderiv(ep, 2,1,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pevttx2(n,jv)   )
                             !   call ogderiv(ep, 2,0,1,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pevtty2(n,jv)   )
                             ! #End
                             ! laplacian
                             LPn = 2.*LP-LPm + 0.5*dt*b1v2(jv)*LPm - dtsq*b0v2(jv)*LP + dtsq*LfP2(n,jv)
                             do na = 0,numberOfAtomicLevels2-1 ! \Delta N^n*E^n
                               LPn = LPn + dtsq*pnec2(jv,na)*qeLap(na)
                             enddo
                             ! time derivatives
                             fLPtt2(n) = fLPtt2(n) + (LPn/(1.+.5*dt*b1v2(jv)) - 2.*LP + LPm)/dtsq
                             ! fLPtt2(n) = fLPtt2(n) + pevttL2(n,jv)
                             ! print *,'error of laplacian with # of levels',numberOfAtomicLevels2,numberOfPolarizationVectors2,pevttxx(n,jv)+pevttyy(n,jv)-pevttL2(n,jv)
                             ! print *, '++++++Dispersive forcing+++++++++++++++'
                             ! print *, 'check P lap: ', n,jv,LP,LPm,b1v2(jv),dtsq,b0v2(jv),LfP2(n,jv),fLPtt2(n)
                             ! x
                             pttxa = 2.*pvx-pvnx + 0.5*dt*b1v2(jv)*pvnx - dtsq*b0v2(jv)*pvx + dtsq*fpvx2(n,jv)
                             do na = 0,numberOfAtomicLevels2-1 ! \Delta N^n*E^n
                               pttxa = pttxa + dtsq*pnec2(jv,na)*qex(na)
                             enddo
                             ! time derivatives
                             fPttx2(n) = fPttx2(n) + (pttxa/(1.+.5*dt*b1v2(jv)) - 2.*pvx + pvnx)/dtsq
                             ! fPttx2(n) = fPttx2(n) + pevttx2(n,jv)
                             ! print *, '++++++Dispersive forcing+++++++++++++++'
                             ! print *, 'check Pttx: ', n,jv,pvx,pvnx,fpvx2(n,jv),fPttx2(n)
                             ! y
                             pttya = 2.*pvy-pvny + 0.5*dt*b1v2(jv)*pvny - dtsq*b0v2(jv)*pvy + dtsq*fpvy2(n,jv)
                             do na = 0,numberOfAtomicLevels2-1 ! \Delta N^n*E^n
                               pttya = pttya + dtsq*pnec2(jv,na)*qey(na)
                             enddo
                             ! time derivatives
                             fPtty2(n) = fPtty2(n) + (pttya/(1.+.5*dt*b1v2(jv)) - 2.*pvy + pvny)/dtsq
                             ! fPtty2(n) = fPtty2(n) + pevtty2(n,jv)
                             ! print *, '++++++Dispersive forcing+++++++++++++++'
                             ! print *, 'check Ptty: ', n,jv,n,jv,pvy,pvny,fpvy2(n,jv),fPtty2(n)
                             ! if( twilightZone.eq.1 )then
                             !   write(*,'("")')
                             !   write(*,'("DI4:RIGHT: j1,j2=",2i3," jv=",i2," n=",i2)') j1,j2,jv,n
                             !   print *, 'ptt diff',(pvec(n,jv)-2.*p2(j1,j2,j3,pc)+p2n(j1,j2,j3,pc))/dtsq - dt**2/12.*pttttv(n,jv)-pevtt2(n,jv)
                             !   print *, 'pttx diff', (pttxa/(1.+.5*dt*b1v2(jv)) - 2.*pvx + pvnx)/dtsq-pevttx2(n,jv)
                             !   print *, 'ptty diff', (pttya/(1.+.5*dt*b1v2(jv)) - 2.*pvy + pvny)/dtsq-pevtty2(n,jv)
                             !   print *, 'ptttt diff', pttttv(n,jv)-pevtttt2(n,jv)
                             !   print *, 'pttL diff',(LPn/(1.+.5*dt*b1v2(jv)) - 2.*LP + LPm)/dtsq-pevttL2(n,jv)
                             ! end if
                           enddo ! jv
                         enddo ! n
                       end if
                     ! first evaluate the equations we want to solve with the wrong values at the ghost points: (assigns f(0:7))
                      f(0)=(u1x+v1y) - (u2x+v2y)
                      ! [ n.Delta( E )/mu ]=0
                      f(1)=(an1*u1Lap+an2*v1Lap)/mu1 - (an1*u2Lap+an2*v2Lap)/mu2
                      ! [ nv X curl( E) /mu ] = 0 
                      f(2)=(v1x-u1y)/mu1 - (v2x-u2y)/mu2
                     ! f(3)=(tau1*u1Lap+tau2*v1Lap)/eps1 - !      (tau1*u2Lap+tau2*v2Lap)/eps2
                      f(3)=( ( tau1*u1Lap +tau2*v1Lap )/epsmu1 - alphaP1*(tau1*fp1(0)+tau2*fp1(1)) ) - ( ( tau1*u2Lap +tau2*v2Lap )/epsmu2 - alphaP2*(tau1*fp2(0)+tau2*fp2(1)) )
                      ! [ Delta( div(E) )*c^2 ] = 0 
                      ! f(4)=(u1xxx+u1xyy+v1xxy+v1yyy)*c1**2- !      (u2xxx+u2xyy+v2xxy+v2yyy)*c2**2
                      f(4)=((u1xxx+u1xyy+v1xxy+v1yyy)/epsmu1-alphaP1*(fPttx1(0) + fPtty1(1))) - ((u2xxx+u2xyy+v2xxy+v2yyy)/epsmu2-alphaP2*(fPttx2(0) + fPtty2(1)))
                      ! print *, 'Divergence of Ptt 0 ??', (fPttx1(0) + fPtty1(1)), (fPttx2(0) + fPtty2(1))
                      ! Dispersive:
                      !   [ (c^2/mu)*{(Delta v).x - (Delta u).y} - (alphaP/mu)*( Py.ttx - Px.tty) ] =0 
                      !    fPttx = P.ttx 
                      !    fPtty = P.tty 
                      f(5)= ( ((v1xxx+v1xyy)-(u1xxy+u1yyy))/epsmu1 - alphaP1*(fPttx1(1) - fPtty1(0)) )/mu1 -( ((v2xxx+v2xyy)-(u2xxy+u2yyy))/epsmu2 - alphaP2*(fPttx2(1) - fPtty2(0)) )/mu2
                      ! f(5)= ( ((v1xxx+2.*v1xyy)-(u1yyy))/epsmu1 - alphaP1*(fPttx1(1) - fPtty1(0)) )/mu1 !      -( ((v2xxx+2.*v2xyy)-(u2yyy))/epsmu2 - alphaP2*(fPttx2(1) - fPtty2(0)) )/mu2
                      if( setDivergenceAtInterfaces.eq.0 )then
                       !  [ nv.( c^2*Delta^2(E) - alphaP*Delta(Ptt) )/mu ] = 0 
                       ! Note: fLptt = Delta( Ptt ) 
                       f(6)= ( (an1*u1LapSq+an2*v1LapSq)/epsmu1  -alphaP1*( an1*fLPtt1(0)+an2*fLPtt1(1) )  )/mu1 -( (an1*u2LapSq+an2*v2LapSq)/epsmu2  -alphaP2*( an1*fLPtt2(0)+an2*fLPtt2(1) )  )/mu2 
                      else
                       f(6)=(u1x+v1y) ! ??
                      end if
                      ! [ tv.( c^4*Delta^2(E) - alphaP*c^2*Delta(P.tt) - alphaP*P.tttt) ]=0 
                      f(7)=(tau1*u1LapSq+tau2*v1LapSq)/(epsmu1**2) - (tau1*u2LapSq+tau2*v2LapSq)/(epsmu2**2) -alphaP1*( ( tau1*fLPtt1(0) + tau2*fLPtt1(1) )/epsmu1 + tau1*fPtttt1(0) + tau2*fPtttt1(1) ) +alphaP2*( ( tau1*fLPtt2(0) + tau2*fLPtt2(1) )/epsmu2 + tau1*fPtttt2(0) + tau2*fPtttt2(1) )
                      ! print *, 'Inside jumps order 4:', alphaP1, alphaP2,fp1(0),fp1(1),fp2(0),fp2(1)
                      if( twilightZone.eq.1 )then
                        call ogderiv(ep, 0,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, uexx )
                        call ogderiv(ep, 0,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, ueyy )
                        call ogderiv(ep, 0,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, vexx )
                        call ogderiv(ep, 0,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, veyy )
                        call ogderiv(ep, 0,3,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, uexxx ) 
                        call ogderiv(ep, 0,2,1,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, uexxy )
                        call ogderiv(ep, 0,1,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, uexyy )
                        call ogderiv(ep, 0,0,3,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, ueyyy )
                        call ogderiv(ep, 0,3,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, vexxx )
                        call ogderiv(ep, 0,1,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, vexyy )
                        call ogderiv(ep, 0,2,1,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, vexxy )
                        call ogderiv(ep, 0,0,3,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, veyyy )
                        call ogderiv(ep, 0,4,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, uexxxx )
                        call ogderiv(ep, 0,2,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, uexxyy )
                        call ogderiv(ep, 0,0,4,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, ueyyyy )
                        call ogderiv(ep, 0,4,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, vexxxx )
                        call ogderiv(ep, 0,2,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, vexxyy )
                        call ogderiv(ep, 0,0,4,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, veyyyy )
                        ueLap = uexx + ueyy
                        veLap = vexx + veyy
                        ueLapSq = uexxxx + 2.*uexxyy + ueyyyy
                        veLapSq = vexxxx + 2.*vexxyy + veyyyy
                        f(1) = f(1) - ( an1*ueLap + an2*veLap )*(1./mu1 - 1./mu2)
                        f(2) = f(2) - ( vex - uey  )*(1./mu1 - 1./mu2)
                        ! print *, 'in eval tz: ', vex,uey, mu1,mu2
                        f(3) = f(3) - ( tau1*ueLap +tau2*veLap )*(1./epsmu1-1./epsmu2) + alphaP1*(tau1*pevttSum1(0)+tau2*pevttSum1(1)) - alphaP2*(tau1*pevttSum2(0)+tau2*pevttSum2(1))
                                    ! print *,'check dispersive forcing'
                                    ! print *,pevttSum1(0)-fp1(0),pevttSum1(1)-fp1(1),pevttSum2(0)-fp2(0),pevttSum2(1)-fp2(1)
                        f(4) = f(4) - ((uexxx+uexyy+vexxy+veyyy)/epsmu1-alphaP1*(pevttxSum1(0) + pevttySum1(1))) + ((uexxx+uexyy+vexxy+veyyy)/epsmu2-alphaP2*(pevttxSum2(0) + pevttySum2(1)))
                        ! f(4) = f(4) - ( uexxx+uexyy+vexxy+veyyy )*(c1**2 - c2**2)
                        f(5) = f(5) - ((vexxx+vexyy)-(uexxy+ueyyy))*(1./(epsmu1*mu1)-1./(epsmu2*mu2)) + (alphaP1/mu1)*( pevttxSum1(1) - pevttySum1(0) ) - (alphaP2/mu2)*( pevttxSum2(1) - pevttySum2(0) )
                        ! f(5) = f(5) - ((vexxx+2.*vexyy)-(ueyyy))*(1./(epsmu1*mu1)-1./(epsmu2*mu2)) !             + (alphaP1/mu1)*( pevttxSum1(1) - pevttySum1(0) ) !             - (alphaP2/mu2)*( pevttxSum2(1) - pevttySum2(0) )
                                    ! print *, fPttx1(1)-pevttxSum1(1),fPtty1(0)-pevttySum1(0),fPttx2(1)-pevttxSum2(1),fPtty2(0)-pevttySum2(0)
                        if( setDivergenceAtInterfaces.eq.0 )then
                          f(6) = f(6) - (an1*ueLapSq+an2*veLapSq)*(1./(epsmu1*mu1) - 1./(epsmu2*mu2) ) + (alphaP1/mu1)*( an1*pevttLSum1(0) + an2*pevttLSum1(1) ) - (alphaP2/mu2)*( an1*pevttLSum2(0) + an2*pevttLSum2(1) ) 
                                      ! print *, fLPtt1(0)-pevttLSum1(0),fLPtt1(1)-pevttLSum1(1),fLPtt2(0)-pevttLSum2(0),fLPtt2(1)-pevttLSum2(1)
                        end if
                        f(7) = f(7) - (tau1*ueLapSq+tau2*veLapSq)*(1./(epsmu1**2) - 1./(epsmu2**2)) + alphaP1*( ( tau1*pevttLSum1(0) + tau2*pevttLSum1(1) )/epsmu1 + tau1*pevttttSum1(0) + tau2*pevttttSum1(1) ) - alphaP2*( ( tau1*pevttLSum2(0) + tau2*pevttLSum2(1) )/epsmu2 + tau1*pevttttSum2(0) + tau2*pevttttSum2(1) ) 
                                    ! print *, fPtttt1(0)-pevttttSum1(0),fPtttt1(1)-pevttttSum1(1),fPtttt2(0)-pevttttSum2(0),fPtttt2(1)-pevttttSum2(1)
                      end if
                     if( debug.gt.7 ) write(debugFile,'(" --> 4cth: j1,j2=",2i4," u1xx,u1yy,u2xx,u2yy=",4e10.2)') j1,j2,u1xx,u1yy,u2xx,u2yy
                       ! '
                      if( debug.gt.3 ) write(debugFile,'(" --> 4cth: i1,i2=",2i4," f(start)=",8e10.2)') i1,i2,f(0),f(1),f(2),f(3),f(4),f(5),f(6),f(7)
                    ! here is the matrix of coefficients for the unknowns u1(-1),v1(-1),u2(-1),v2(-1)
                    ! Solve:
                    !     
                    !       A [ U ] = A [ U(old) ] - [ f ]
                    ! Equation 0: 
                    ! 0  [ u.x + v.y ] = 0
                    a8(0,0) = -is*8.*rx1*dx141(axis1)     ! coeff of u1(-1) from [u.x+v.y] 
                    a8(0,1) = -is*8.*ry1*dx141(axis1)     ! coeff of v1(-1) from [u.x+v.y] 
                    a8(0,4) =  is*rx1*dx141(axis1)        ! u1(-2)
                    a8(0,5) =  is*ry1*dx141(axis1)        ! v1(-2) 
                    a8(0,2) =  js*8.*rx2*dx241(axis2)     ! coeff of u2(-1) from [u.x+v.y] 
                    a8(0,3) =  js*8.*ry2*dx241(axis2) 
                    a8(0,6) = -js*   rx2*dx241(axis2) 
                    a8(0,7) = -js*   ry2*dx241(axis2) 
                    ! 1  [ (u.xx + u.yy)/mu ] = 0
                    a8(1,0) = 16.*dx142(axis1)/mu1         ! coeff of u1(-1) from [u.xx + u.yy]
                    a8(1,1) = 0. 
                    a8(1,4) =    -dx142(axis1)/mu1         ! coeff of u1(-2) from [u.xx + u.yy]
                    a8(1,5) = 0. 
                    a8(1,2) =-16.*dx242(axis2)/mu2         ! coeff of u2(-1) from [u.xx + u.yy]
                    a8(1,3) = 0. 
                    a8(1,6) =     dx242(axis2)/mu2         ! coeff of u2(-2) from [u.xx + u.yy]
                    a8(1,7) = 0. 
                    ! Equation 2: 
                    curl1um1 =  is*8.*ry1*dx141(axis1)   ! coeff of u(-1) from v.x - u.y 
                    curl1vm1 = -is*8.*rx1*dx141(axis1)   ! coeff of v(-1) from v.x - u.y 
                    curl1um2 = -is*   ry1*dx141(axis1)   ! coeff of u(-2) from v.x - u.y 
                    curl1vm2 =  is*   rx1*dx141(axis1)   ! coeff of v(-2) from v.x - u.y
                    curl2um1 =  js*8.*ry2*dx241(axis2)  ! coeff of u(-1) from v.x - u.y 
                    curl2vm1 = -js*8.*rx2*dx241(axis2)  ! coeff of v(-1) from v.x - u.y 
                    curl2um2 = -js*   ry2*dx241(axis2)  ! coeff of u(-2) from v.x - u.y 
                    curl2vm2 =  js*   rx2*dx241(axis2)  ! coeff of v(-2) from v.x - u.y
                    ! 2  [ (v.x - u.y)/mu ] =0 
                    a8(2,0) =  curl1um1/mu1 
                    a8(2,1) =  curl1vm1/mu1 
                    a8(2,4) =  curl1um2/mu1 
                    a8(2,5) =  curl1vm2/mu1 
                    a8(2,2) = -curl2um1/mu2
                    a8(2,3) = -curl2vm1/mu2
                    a8(2,6) = -curl2um2/mu2
                    a8(2,7) = -curl2vm2/mu2
                    !- a8(2,0) =  is*8.*ry1*dx141(axis1)/mu1 
                    !- a8(2,1) = -is*8.*rx1*dx141(axis1)/mu1     ! coeff of v1(-1) from [v.x - u.y] 
                    !- a8(2,4) = -is*   ry1*dx141(axis1)/mu1 
                    !- a8(2,5) =  is*   rx1*dx141(axis1)/mu1 
                    !- a8(2,2) = -js*8.*ry2*dx241(axis2)/mu2
                    !- a8(2,3) =  js*8.*rx2*dx241(axis2)/mu2
                    !- a8(2,6) =  js*   ry2*dx241(axis2)/mu2
                    !- a8(2,7) = -js*   rx2*dx241(axis2)/mu2
                    aLap0=16.*dx142(axis1)  ! coeff of w1(-1) 
                    aLap1=-1.*dx142(axis1)  ! coeff of w1(-2)
                    bLap0=16.*dx242(axis2)  ! coeff of w2(-1) 
                    bLap1=-1.*dx242(axis2)  ! coeff of w2(-1) 
                    aLapSq0= ( -4./(dx1(axis1)**4) )
                    aLapSq1= (  1./(dx1(axis1)**4) )
                    bLapSq0= ( -4./(dx2(axis2)**4) )
                    bLapSq1= (  1./(dx2(axis2)**4) )
                    ! -------------- Equation 3 -----------------------
                    !   [ tau.{ (uv.xx+uv.yy)/eps -alphaP*P.tt } ] = 0
                    !    P.tt = c4PttLEsum * L(E) + c4PttLLEsum* L^2(E) + ...
                    c4PttLEsum1 = 0.
                    c4PttLLEsum1 = 0.
                    c4PttLEsum2 = 0.
                    c4PttLLEsum2 = 0.
                    a8(3,0) = tau1*( aLap0*( 1./epsmu1 -alphaP1*c4PttLEsum1/epsmu1 ) - aLapSq0*alphaP1*c4PttLLEsum1/epsmu1**2 )
                    a8(3,1) = tau2*( aLap0*( 1./epsmu1 -alphaP1*c4PttLEsum1/epsmu1 ) - aLapSq0*alphaP1*c4PttLLEsum1/epsmu1**2 )
                    a8(3,4) = tau1*( aLap1*( 1./epsmu1 -alphaP1*c4PttLEsum1/epsmu1 ) - aLapSq1*alphaP1*c4PttLLEsum1/epsmu1**2 )
                    a8(3,5) = tau2*( aLap1*( 1./epsmu1 -alphaP1*c4PttLEsum1/epsmu1 ) - aLapSq1*alphaP1*c4PttLLEsum1/epsmu1**2 )
                    a8(3,2) =-tau1*( bLap0*( 1./epsmu2 -alphaP2*c4PttLEsum2/epsmu2 ) - bLapSq0*alphaP2*c4PttLLEsum2/epsmu2**2 )
                    a8(3,3) =-tau2*( bLap0*( 1./epsmu2 -alphaP2*c4PttLEsum2/epsmu2 ) - bLapSq0*alphaP2*c4PttLLEsum2/epsmu2**2 )
                    a8(3,6) =-tau1*( bLap1*( 1./epsmu2 -alphaP2*c4PttLEsum2/epsmu2 ) - bLapSq1*alphaP2*c4PttLLEsum2/epsmu2**2 )
                    a8(3,7) =-tau2*( bLap1*( 1./epsmu2 -alphaP2*c4PttLEsum2/epsmu2 ) - bLapSq1*alphaP2*c4PttLLEsum2/epsmu2**2 )
                    ! -------------- Equation 4 -----------------------
                    !    [ (u.xx+u.yy).x + (v.xx+v.yy).y ] = 0
                    dxxx1by2i = 1./(dx1(axis1)**3*2); 
                    dxxx2by2i = 1./(dx2(axis2)**3*2); 
                    ! Order u1(-1), v1(-1), u2(-1), v2(-1), u1(-2), v1(-2), u2(-2), v2(-2), 
                    a8(4,0)= ( is*rx1 *2.*dxxx1by2i )*c1**2
                    a8(4,1)= 0.
                    a8(4,2)=-( js*rx2* 2.*dxxx2by2i )*c2**2
                    a8(4,3)= 0.
                    a8(4,4)= (-is*rx1    *dxxx1by2i )*c1**2
                    a8(4,5)= 0.
                    a8(4,6)=-(-js*rx2    *dxxx2by2i )*c2**2  
                    a8(4,7)= 0.
                    ! ---------------- Equation 5 (2nd-order) -----------------
                    !   [ ( {(Delta v).x - (Delta u).y}/(epsmu) - alphaP*( Py.ttx - Px.tty) )/mu ] =0 
                    !
                    ! check me: 
                    aLapX0 =  is*rx1*2.*dx122(axis1)*dx112(axis1) 
                    bLapY0 = -is*ry1*2.*dx122(axis1)*dx112(axis1)
                    aLapX1 = -is*rx1   *dx122(axis1)*dx112(axis1)
                    bLapY1 =  is*ry1   *dx122(axis1)*dx112(axis1)
                    cLapX0 =  js*rx2*2.*dx222(axis2)*dx212(axis2) 
                    dLapY0 = -js*ry2*2.*dx222(axis2)*dx212(axis2)
                    cLapX1 = -js*rx2   *dx222(axis2)*dx212(axis2)
                    dLapY1 =  js*ry2   *dx222(axis2)*dx212(axis2)
                    !     P.tt = c2PttLEsum * L(E)
                    c2PttLEsum1 = 0.
                    c2PttEsum1 = 0.
                    c2PttLEsum2 = 0.
                    c2PttEsum2 = 0.
                    eqnCoeff = ( 1./epsmu1 - alphaP1*c2PttLEsum1/epsmu1 )/mu1 
                    eqnCoeffb = -alphaP1*c2PttEsum1/mu1 ! added sept 16, 2018 
                    a8(5,0)=-bLapY0*eqnCoeff + curl1um1*eqnCoeffb  
                    a8(5,1)= aLapX0*eqnCoeff + curl1vm1*eqnCoeffb
                    a8(5,4)=-bLapY1*eqnCoeff + curl1um2*eqnCoeffb 
                    a8(5,5)= aLapX1*eqnCoeff + curl1vm2*eqnCoeffb
                    eqnCoeff = ( 1./epsmu2 - alphaP2*c2PttLEsum2/epsmu2 )/mu2 
                    eqnCoeffb = -alphaP2*c2PttEsum2/mu2 ! added sept 16, 2018 
                    a8(5,2)=-(-dLapY0*eqnCoeff + curl2um1*eqnCoeffb)
                    a8(5,3)=-( cLapX0*eqnCoeff + curl2vm1*eqnCoeffb)
                    a8(5,6)=-(-dLapY1*eqnCoeff + curl2um2*eqnCoeffb)
                    a8(5,7)=-( cLapX1*eqnCoeff + curl2vm2*eqnCoeffb)
                    ! ------- Equation 6 -----
                    !  [ nv.( c^2*Delta^2(E) - alphaP*Delta(Ptt) )/mu ] = 0 
                    ! 6  [ n.Delta^2 u/eps ] = 0
                    ! use Eqn 6 
                    ! NOTE: LE = c^2*Delta(E) and LLE = (c^4*Delta^2) E 
                    ! Note: the coeff of L(E) in Delta(Ptt) is the coeff of E in Ptt
                    ! Note: the coeff of LL(E) in Delta(Ptt) is the coeff of LE in Ptt
                    c2PttEsum1 = 0.
                    c2PttLEsum1 = 0.
                    c2PttEsum2 = 0.
                    c2PttLEsum2 = 0.
                    a8(6,0) = an1*( aLapSq0/epsmu1 -alphaP1*( c2PttEsum1*aLap0 + c2PttLEsum1*aLapSq0/epsmu1 ) )/mu1
                    a8(6,1) = an2*( aLapSq0/epsmu1 -alphaP1*( c2PttEsum1*aLap0 + c2PttLEsum1*aLapSq0/epsmu1 ) )/mu1
                    a8(6,4) = an1*( aLapSq1/epsmu1 -alphaP1*( c2PttEsum1*aLap1 + c2PttLEsum1*aLapSq1/epsmu1 ) )/mu1
                    a8(6,5) = an2*( aLapSq1/epsmu1 -alphaP1*( c2PttEsum1*aLap1 + c2PttLEsum1*aLapSq1/epsmu1 ) )/mu1
                    a8(6,2) =-an1*( bLapSq0/epsmu2 -alphaP2*( c2PttEsum2*bLap0 + c2PttLEsum2*bLapSq0/epsmu2 ) )/mu2
                    a8(6,3) =-an2*( bLapSq0/epsmu2 -alphaP2*( c2PttEsum2*bLap0 + c2PttLEsum2*bLapSq0/epsmu2 ) )/mu2
                    a8(6,6) =-an1*( bLapSq1/epsmu2 -alphaP2*( c2PttEsum2*bLap1 + c2PttLEsum2*bLapSq1/epsmu2 ) )/mu2
                    a8(6,7) =-an2*( bLapSq1/epsmu2 -alphaP2*( c2PttEsum2*bLap1 + c2PttLEsum2*bLapSq1/epsmu2 ) )/mu2
                    ! ------- Equation 7 ------
                    ! [ tv.( c^4*Delta^2(E) - alphaP*c^2*Delta(P.tt) - alphaP*P.tttt) ]=0 
                    ! 7  [ tau.Delta^2 v/eps^2 ] = 0 
                    ! Note: the coeff of L(E) in Delta(Ptt) is the coeff of E in Ptt
                    ! Note: the coeff of LL(E) in Delta(Ptt) is the coeff of LE in Ptt
                    c2PttEsum1 = 0.
                    c2PttttLEsum1 = 0.
                    c2PttLEsum1 = 0.
                    c2PttttLLEsum1 = 0.
                    c2PttEsum2 = 0.
                    c2PttttLEsum2 = 0.
                    c2PttLEsum2 = 0.
                    c2PttttLLEsum2 = 0.
                    coeffLap1   =              -alphaP1*(  c2PttEsum1 + c2PttttLEsum1  )/epsmu1
                    coeffLapSq1 = 1./epsmu1**2 -alphaP1*( c2PttLEsum1 + c2PttttLLEsum1 )/epsmu1**2
                    coeffLap2   =              -alphaP2*(  c2PttEsum2 + c2PttttLEsum2  )/epsmu2
                    coeffLapSq2 = 1./epsmu2**2 -alphaP2*( c2PttLEsum2 + c2PttttLLEsum2 )/epsmu2**2
                    a8(7,0) = tau1*( coeffLapSq1*aLapSq0 + coeffLap1*aLap0 )
                    a8(7,1) = tau2*( coeffLapSq1*aLapSq0 + coeffLap1*aLap0 )
                    a8(7,4) = tau1*( coeffLapSq1*aLapSq1 + coeffLap1*aLap1 )
                    a8(7,5) = tau2*( coeffLapSq1*aLapSq1 + coeffLap1*aLap1 )
                    a8(7,2) =-tau1*( coeffLapSq2*bLapSq0 + coeffLap2*bLap0 )
                    a8(7,3) =-tau2*( coeffLapSq2*bLapSq0 + coeffLap2*bLap0 )
                    a8(7,6) =-tau1*( coeffLapSq2*bLapSq1 + coeffLap2*bLap1 )
                    a8(7,7) =-tau2*( coeffLapSq2*bLapSq1 + coeffLap2*bLap1 )
                    if( debug>7 .and. i2.le.0 )then
                      write(*,*) "mla24r: Matrix a8"
                      do n=0,7
                        write(*,'(8(1pe10.2))') (a8(n,nn),nn=0,7)
                      end do 
                    end if 
                    ! --- check matrix coefficients by delta function approach ----
                    if( checkCoeff.eq.1 )then
                     !! numberOfEquations=8
                     !! checkCoefficients(i1,i2,i3, j1,j2,j3,numberOfEquations,a8,evalDispersiveInterfaceEquations24r )
                    end if
                    ! -- save current (wrong) ghost values in q()
                    q(0) = u1(i1-is1,i2-is2,i3,ex)
                    q(1) = u1(i1-is1,i2-is2,i3,ey)
                    q(2) = u2(j1-js1,j2-js2,j3,ex)
                    q(3) = u2(j1-js1,j2-js2,j3,ey)
                    q(4) = u1(i1-2*is1,i2-2*is2,i3,ex)
                    q(5) = u1(i1-2*is1,i2-2*is2,i3,ey)
                    q(6) = u2(j1-2*js1,j2-2*js2,j3,ex)
                    q(7) = u2(j1-2*js1,j2-2*js2,j3,ey)
                    !- if( debug.gt.4 )then
                    !-   write(debugFile,'(" --> mla24r: i1,i2=",2i4," q=",8e10.2)') i1,i2,(q(n),n=0,7)
                    !- end if
                    if( debug.gt.7 )then
                      write(*,'(" --> before mla24r: i1,i2=",2i4," RHS(A)=",8(1pe13.5))') i1,i2,(f(n),n=0,7)
                    end if
                    ! subtract off the contributions from the initial (wrong) values at the ghost points:
                    do n=0,7
                      f(n) = (a8(n,0)*q(0)+a8(n,1)*q(1)+a8(n,2)*q(2)+a8(n,3)*q(3)+a8(n,4)*q(4)+a8(n,5)*q(5)+a8(n,6)*q(6)+a8(n,7)*q(7)) - f(n)
                    end do
                    if( debug.gt.7 )then
                      write(*,'(" --> after mla24r: i1,i2=",2i4," RHS(A)=",8(1pe13.5))') i1,i2,(f(n),n=0,7)
                    end if
                    ! solve A Q = F
                    ! factor the matrix
                    numberOfEquations=8
                    call dgeco( a8(0,0), numberOfEquations, numberOfEquations, ipivot8(0),rcond,work(0))
                    if( debug.gt.3 ) then
                      write(debugFile,'(" --> mla24r: i1,i2=",2i4," rcond=",e10.2)') i1,i2,rcond
                    end if
                    ! solve A Q = F
                    job=0
                    numberOfEquations=8
                    call dgesl( a8(0,0), numberOfEquations, numberOfEquations, ipivot8(0), f(0), job)
                    if( useJacobiUpdate.eq.0 )then
                      u1(i1-is1,i2-is2,i3,ex)=f(0)
                      u1(i1-is1,i2-is2,i3,ey)=f(1)
                      u2(j1-js1,j2-js2,j3,ex)=f(2)
                      u2(j1-js1,j2-js2,j3,ey)=f(3)
                      u1(i1-2*is1,i2-2*is2,i3,ex)=f(4)
                      u1(i1-2*is1,i2-2*is2,i3,ey)=f(5)
                      u2(j1-2*js1,j2-2*js2,j3,ex)=f(6)
                      u2(j1-2*js1,j2-2*js2,j3,ey)=f(7)
                    else
                      ! Jacobi-update
                      wk1(i1-is1,i2-is2,i3,ex)    =f(0)
                      wk1(i1-is1,i2-is2,i3,ey)    =f(1)
                      wk2(j1-js1,j2-js2,j3,ex)    =f(2)
                      wk2(j1-js1,j2-js2,j3,ey)    =f(3)
                      wk1(i1-2*is1,i2-2*is2,i3,ex)=f(4)
                      wk1(i1-2*is1,i2-2*is2,i3,ey)=f(5)
                      wk2(j1-2*js1,j2-2*js2,j3,ex)=f(6)
                      wk2(j1-2*js1,j2-2*js2,j3,ey)=f(7)
                    end if
                    ! if( debug>3 .and. twilightZone.eq.1 )then
                    !   ! check errors
                    !   call ogderiv(ep, 0,0,0,0, xy1(i1-is1,i2-is2,i3,0),xy1(i1-is1,i2-is2,i3,1),0.,t, ex, evv(0) )
                    !   call ogderiv(ep, 0,0,0,0, xy1(i1-is1,i2-is2,i3,0),xy1(i1-is1,i2-is2,i3,1),0.,t, ey, evv(1) )
                    !   call ogderiv(ep, 0,0,0,0, xy2(j1-js1,j2-js2,j3,0),xy2(j1-js1,j2-js2,j3,1),0.,t, ex, evv(2) )
                    !   call ogderiv(ep, 0,0,0,0, xy2(j1-js1,j2-js2,j3,0),xy2(j1-js1,j2-js2,j3,1),0.,t, ey, evv(3) )
                    !   call ogderiv(ep, 0,0,0,0, xy1(i1-2*is1,i2-2*is2,i3,0),xy1(i1-2*is1,i2-2*is2,i3,1),0.,t, ex, evv(4) )
                    !   call ogderiv(ep, 0,0,0,0, xy1(i1-2*is1,i2-2*is2,i3,0),xy1(i1-2*is1,i2-2*is2,i3,1),0.,t, ey, evv(5) )
                    !   call ogderiv(ep, 0,0,0,0, xy2(j1-2*js1,j2-2*js2,j3,0),xy2(j1-2*js1,j2-2*js2,j3,1),0.,t, ex, evv(6) )
                    !   call ogderiv(ep, 0,0,0,0, xy2(j1-2*js1,j2-2*js2,j3,0),xy2(j1-2*js1,j2-2*js2,j3,1),0.,t, ey, evv(7) )
                    !   write(*,'("gdm24r: Errors in u1(-1), v1(-1), u2(-1), v2(-1), u1(-2), v1(-2), u2(-2), v2(-2)")') 
                    !   maxErr=0.
                    !   do n=0,7
                    !     maxErr =max(maxErr,abs(evv(n)-f(n)))
                    !   end do
                    !   write(*,'("gdm24r: i1,i2=",2i4," err= ",8e8.1," -> maxErr=",e8.1)') i1,i2, (abs(evv(n)-f(n)),n=0,7),maxErr
                    ! end if
                    !-  if( .false. .and. debug.gt.2 )then 
                    !-    ! --- check residuals in the jump conditions ----
                    !-
                    !-    ! Evaluate the jump conditions using the new values at the ghost points 
                    !-    !! evaluateDispersiveInterfaceEquations2dOrder4()
                    !-    write(debugFile,'(" JUMP-residuals: i1,i2=",2i4," f(re-eval)=",8e10.2)') i1,i2,f(0),f(1),f(2),f(3),f(4),f(5),f(6),f(7)
                    !-  end if
                   ! ******************************************************
                   ! NOTE: This is just the non-dispersive case ***FIX ME FOR GDM ***
                   ! 
                   ! solve for Hz
                   !  [ w.n/eps ] = 0
                   !  [ lap(w)/eps ] = 0
                   !  [ lap(w).n/eps**2 ] = 0
                   !  [ lapSq(w)/eps**2 ] = 0
                    ! first evaluate the equations we want to solve with the wrong values at the ghost points:
                     ! These derivatives are computed to 4th-order accuracy
                       ww1=u1(i1,i2,i3,hz) ! in the rectangular case just eval the solution
                        w1x = (u1(i1-2,i2,i3,hz)-8.*u1(i1-1,i2,i3,hz)+8.*u1(i1+1,i2,i3,hz)-u1(i1+2,i2,i3,hz))/(12.*dx1(0))
                        w1y = (u1(i1,i2-2,i3,hz)-8.*u1(i1,i2-1,i3,hz)+8.*u1(i1,i2+1,i3,hz)-u1(i1,i2+2,i3,hz))/(12.*dx1(1))
                        w1xx = (-u1(i1-2,i2,i3,hz)+16.*u1(i1-1,i2,i3,hz)-30.*u1(i1,i2,i3,hz)+16.*u1(i1+1,i2,i3,hz)-u1(i1+2,i2,i3,hz))/(12.*dx1(0)**2)
                        w1yy = (-u1(i1,i2-2,i3,hz)+16.*u1(i1,i2-1,i3,hz)-30.*u1(i1,i2,i3,hz)+16.*u1(i1,i2+1,i3,hz)-u1(i1,i2+2,i3,hz))/(12.*dx1(1)**2)
                      w1Lap = w1xx+ w1yy
                       ww2=u2(j1,j2,j3,hz) ! in the rectangular case just eval the solution
                        w2x = (u2(j1-2,j2,j3,hz)-8.*u2(j1-1,j2,j3,hz)+8.*u2(j1+1,j2,j3,hz)-u2(j1+2,j2,j3,hz))/(12.*dx2(0))
                        w2y = (u2(j1,j2-2,j3,hz)-8.*u2(j1,j2-1,j3,hz)+8.*u2(j1,j2+1,j3,hz)-u2(j1,j2+2,j3,hz))/(12.*dx2(1))
                        w2xx = (-u2(j1-2,j2,j3,hz)+16.*u2(j1-1,j2,j3,hz)-30.*u2(j1,j2,j3,hz)+16.*u2(j1+1,j2,j3,hz)-u2(j1+2,j2,j3,hz))/(12.*dx2(0)**2)
                        w2yy = (-u2(j1,j2-2,j3,hz)+16.*u2(j1,j2-1,j3,hz)-30.*u2(j1,j2,j3,hz)+16.*u2(j1,j2+1,j3,hz)-u2(j1,j2+2,j3,hz))/(12.*dx2(1)**2)
                      w2Lap = w2xx+ w2yy
                     ! These derivatives are computed to 2nd-order accuracy
                       ww1=u1(i1,i2,i3,hz) ! in the rectangular case just eval the solution
                        w1xxx = (-u1(i1-2,i2,i3,hz)+2.*u1(i1-1,i2,i3,hz)-2.*u1(i1+1,i2,i3,hz)+u1(i1+2,i2,i3,hz))/(2.*dx1(0)**3)
                        w1xxy = ((-u1(i1-1,i2-1,i3,hz)+u1(i1-1,i2+1,i3,hz))/(2.*dx1(1))-2.*(-u1(i1,i2-1,i3,hz)+u1(i1,i2+1,i3,hz))/(2.*dx1(1))+(-u1(i1+1,i2-1,i3,hz)+u1(i1+1,i2+1,i3,hz))/(2.*dx1(1)))/(dx1(0)**2)
                        w1xyy = (-(u1(i1-1,i2-1,i3,hz)-2.*u1(i1-1,i2,i3,hz)+u1(i1-1,i2+1,i3,hz))/(dx1(1)**2)+(u1(i1+1,i2-1,i3,hz)-2.*u1(i1+1,i2,i3,hz)+u1(i1+1,i2+1,i3,hz))/(dx1(1)**2))/(2.*dx1(0))
                        w1yyy = (-u1(i1,i2-2,i3,hz)+2.*u1(i1,i2-1,i3,hz)-2.*u1(i1,i2+1,i3,hz)+u1(i1,i2+2,i3,hz))/(2.*dx1(1)**3)
                        w1xxxx = (u1(i1-2,i2,i3,hz)-4.*u1(i1-1,i2,i3,hz)+6.*u1(i1,i2,i3,hz)-4.*u1(i1+1,i2,i3,hz)+u1(i1+2,i2,i3,hz))/(dx1(0)**4)
                        w1xxyy = ((u1(i1-1,i2-1,i3,hz)-2.*u1(i1-1,i2,i3,hz)+u1(i1-1,i2+1,i3,hz))/(dx1(1)**2)-2.*(u1(i1,i2-1,i3,hz)-2.*u1(i1,i2,i3,hz)+u1(i1,i2+1,i3,hz))/(dx1(1)**2)+(u1(i1+1,i2-1,i3,hz)-2.*u1(i1+1,i2,i3,hz)+u1(i1+1,i2+1,i3,hz))/(dx1(1)**2))/(dx1(0)**2)
                        w1yyyy = (u1(i1,i2-2,i3,hz)-4.*u1(i1,i2-1,i3,hz)+6.*u1(i1,i2,i3,hz)-4.*u1(i1,i2+1,i3,hz)+u1(i1,i2+2,i3,hz))/(dx1(1)**4)
                      w1LapSq = w1xxxx +2.* w1xxyy + w1yyyy
                       ww2=u2(j1,j2,j3,hz) ! in the rectangular case just eval the solution
                        w2xxx = (-u2(j1-2,j2,j3,hz)+2.*u2(j1-1,j2,j3,hz)-2.*u2(j1+1,j2,j3,hz)+u2(j1+2,j2,j3,hz))/(2.*dx2(0)**3)
                        w2xxy = ((-u2(j1-1,j2-1,j3,hz)+u2(j1-1,j2+1,j3,hz))/(2.*dx2(1))-2.*(-u2(j1,j2-1,j3,hz)+u2(j1,j2+1,j3,hz))/(2.*dx2(1))+(-u2(j1+1,j2-1,j3,hz)+u2(j1+1,j2+1,j3,hz))/(2.*dx2(1)))/(dx2(0)**2)
                        w2xyy = (-(u2(j1-1,j2-1,j3,hz)-2.*u2(j1-1,j2,j3,hz)+u2(j1-1,j2+1,j3,hz))/(dx2(1)**2)+(u2(j1+1,j2-1,j3,hz)-2.*u2(j1+1,j2,j3,hz)+u2(j1+1,j2+1,j3,hz))/(dx2(1)**2))/(2.*dx2(0))
                        w2yyy = (-u2(j1,j2-2,j3,hz)+2.*u2(j1,j2-1,j3,hz)-2.*u2(j1,j2+1,j3,hz)+u2(j1,j2+2,j3,hz))/(2.*dx2(1)**3)
                        w2xxxx = (u2(j1-2,j2,j3,hz)-4.*u2(j1-1,j2,j3,hz)+6.*u2(j1,j2,j3,hz)-4.*u2(j1+1,j2,j3,hz)+u2(j1+2,j2,j3,hz))/(dx2(0)**4)
                        w2xxyy = ((u2(j1-1,j2-1,j3,hz)-2.*u2(j1-1,j2,j3,hz)+u2(j1-1,j2+1,j3,hz))/(dx2(1)**2)-2.*(u2(j1,j2-1,j3,hz)-2.*u2(j1,j2,j3,hz)+u2(j1,j2+1,j3,hz))/(dx2(1)**2)+(u2(j1+1,j2-1,j3,hz)-2.*u2(j1+1,j2,j3,hz)+u2(j1+1,j2+1,j3,hz))/(dx2(1)**2))/(dx2(0)**2)
                        w2yyyy = (u2(j1,j2-2,j3,hz)-4.*u2(j1,j2-1,j3,hz)+6.*u2(j1,j2,j3,hz)-4.*u2(j1,j2+1,j3,hz)+u2(j1,j2+2,j3,hz))/(dx2(1)**4)
                      w2LapSq = w2xxxx +2.* w2xxyy + w2yyyy
                     f(0)=(an1*w1x+an2*w1y)/eps1 - (an1*w2x+an2*w2y)/eps2
                     f(1)=w1Lap/eps1 - w2Lap/eps2
                     f(2)=(an1*(w1xxx+w1xyy)+an2*(w1xxy+w1yyy))/eps1**2 - (an1*(w2xxx+w2xyy)+an2*(w2xxy+w2yyy))/eps2**2
                     f(3)=w1LapSq/eps1**2 - w2LapSq/eps2**2
                     if( twilightZone.eq.1 )then
                       call ogderiv(ep, 0,1,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, hz, wex  )
                       call ogderiv(ep, 0,0,1,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, hz, wey  )
                       call ogderiv(ep, 0,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, hz, wexx )
                       call ogderiv(ep, 0,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, hz, weyy )
                       call ogderiv(ep, 0,3,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, hz, wexxx )
                       call ogderiv(ep, 0,2,1,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, hz, wexxy )
                       call ogderiv(ep, 0,1,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, hz, wexyy )
                       call ogderiv(ep, 0,0,3,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, hz, weyyy )
                       call ogderiv(ep, 0,4,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, hz, wexxxx )
                       call ogderiv(ep, 0,2,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, hz, wexxyy )
                       call ogderiv(ep, 0,0,4,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, hz, weyyyy )
                       weLap = wexx + weyy
                       weLapSq = wexxxx + 2.*wexxyy + weyyyy
                       f(0) = f(0) - (an1*wex+an2*wey)*(1./eps1 - 1./eps2)
                       f(1) = f(1) - ( weLap )*(1./eps1 - 1./eps2)
                       f(2) = f(2) - (an1*(wexxx+wexyy)+an2*(wexxy+weyyy))*(1./eps1**2 - 1./eps2**2)
                       f(3) = f(3) - weLapSq*(1./eps1**2 - 1./eps2**2)
                     end if
                    ! form the matrix for computing Hz
                    ! 1: [ w.n/eps ] = 0
                    a0 = dx141(axis1)/eps1
                    b0 = dx241(axis2)/eps2
                    a4h(0,0) = -is*8.*a0
                    a4h(0,2) =  is*   a0
                    a4h(0,1) =  js*8.*b0
                    a4h(0,3) = -js*   b0
                    ! 2: [ lap(w)/eps ] = 0 
                    aLap0=16.*dx142(axis1)  ! coeff of w1(-1) 
                    aLap1=-1.*dx142(axis1)  ! coeff of w1(-2)
                    bLap0=16.*dx242(axis2)  ! coeff of w2(-1) 
                    bLap1=-1.*dx242(axis2)  ! coeff of w2(-1) 
                    a4h(1,0) = aLap0/eps1    ! coeff of w1(-1) 
                    a4h(1,2) = aLap1/eps1    ! coeff of w1(-2)
                    a4h(1,1) =-bLap0/eps2    ! coeff of w2(-1) 
                    a4h(1,3) =-bLap1/eps2    ! coeff of w2(-1) 
                    ! 3:  [ (an1*(w.xx+w.yy).x + an2.(w.xx+w.yy).y)/eps**2 ] = 0
                    !  a4h(2,0)= (an1*aLapX0+an2*bLapY0)/eps1**2  ! coeff of w1(-1) 
                    !  a4h(2,1)=-(an1*cLapX0+an2*dLapY0)/eps2**2  ! coeff of w2(-1)
                    !  a4h(2,2)= (an1*aLapX1+an2*bLapY1)/eps1**2  ! coeff of w1(-2)
                    !  a4h(2,3)=-(an1*cLapX1+an2*dLapY1)/eps2**2  ! coeff of w2(-2)
                    a4h(2,0)=  is*(  1./(dx1(axis1)**3) +1./(dx1(axis1)*dx1(axis1p1)**2) )/eps1**2  ! coeff of w1(-1) aLapX0
                    a4h(2,2)=  is*( -.5/(dx1(axis1)**3)                                  )/eps1**2  ! coeff of w1(-2) aLapX1
                    a4h(2,1)= -js*(  1./(dx2(axis2)**3) +1./(dx2(axis2)*dx2(axis2p1)**2) )/eps2**2  ! coeff of w2(-1) cLapX0 
                    a4h(2,3)= -js*( -.5/(dx2(axis2)**3)                                  )/eps2**2  ! coeff of w2(-2) cLapX1
                    ! 4 [ lapSq(w)/eps**2 ] = 0   [ w_xxxx + 2 * w_xxyy + w_yyyy ]
                    aLapSq0= ( -4./(dx1(axis1)**4) -4./(dx1(axis1)**2 * dx1(axis1p1)**2 ) )
                    aLapSq1= (  1./(dx1(axis1)**4) )
                    bLapSq0= ( -4./(dx2(axis2)**4) -4./(dx2(axis2)**2 * dx2(axis2p1)**2 ) )
                    bLapSq1= (  1./(dx2(axis2)**4) )
                    a4h(3,0) = aLapSq0/eps1**2
                    a4h(3,2) = aLapSq1/eps1**2
                    a4h(3,1) =-bLapSq0/eps2**2
                    a4h(3,3) =-bLapSq1/eps2**2
                    q(0) = u1(i1-is1,i2-is2,i3,hz)
                    q(1) = u2(j1-js1,j2-js2,j3,hz)
                    q(2) = u1(i1-2*is1,i2-2*is2,i3,hz)
                    q(3) = u2(j1-2*js1,j2-2*js2,j3,hz)
                    ! subtract off the contributions from the wrong values at the ghost points:
                    do n=0,3
                      f(n) = (a4h(n,0)*q(0)+a4h(n,1)*q(1)+a4h(n,2)*q(2)+a4h(n,3)*q(3)) - f(n)
                    end do
                    ! write(*,'(" a4h=",4(e9.2,1x))') ((a4h(i,j),j=0,3),i=0,3)
                    ! factor the matrix
                    ! numberOfEquations=4
                    call dgeco( a4h(0,0), 4, 4, ipivot4h(0),rcond,work(0))
                    ! write(*,'("rcond=",e12.4)') rcond
                    ! solve
                    job=0
                    call dgesl( a4h(0,0), 4, 4, ipivot4h(0), f(0), job)
                    if( useJacobiUpdate.eq.0 )then
                      u1(i1-  is1,i2-  is2,i3,hz)=f(0)
                      u2(j1-  js1,j2-  js2,j3,hz)=f(1)
                      u1(i1-2*is1,i2-2*is2,i3,hz)=f(2)
                      u2(j1-2*js1,j2-2*js2,j3,hz)=f(3)
                    else
                      ! Jacobi update -- save answer in work space
                      wk1(i1-  is1,i2-  is2,i3,hz)=f(0)
                      wk2(j1-  js1,j2-  js2,j3,hz)=f(1)
                      wk1(i1-2*is1,i2-2*is2,i3,hz)=f(2)
                      wk2(j1-2*js1,j2-2*js2,j3,hz)=f(3)
                    end if
                     end if
                     j1=j1+1
                    end do
                    j2=j2+1
                   end do
                  ! =============== end loops =======================
                  if( useJacobiUpdate.ne.0 )then
                    ! Jacobi-update: now fill in values 
                     i3=n3a
                     j3=m3a
                     j2=m2a
                     do i2=n2a,n2b
                      j1=m1a
                      do i1=n1a,n1b
                       if( mask1(i1,i2,i3).gt.0 .and. mask2(j1,j2,j3).gt.0 )then
                      u1(i1-is1,i2-is2,i3,ex)=wk1(i1-is1,i2-is2,i3,ex)
                      u1(i1-is1,i2-is2,i3,ey)=wk1(i1-is1,i2-is2,i3,ey)
                      u2(j1-js1,j2-js2,j3,ex)=wk2(j1-js1,j2-js2,j3,ex)
                      u2(j1-js1,j2-js2,j3,ey)=wk2(j1-js1,j2-js2,j3,ey)
                      u1(i1-2*is1,i2-2*is2,i3,ex)=wk1(i1-2*is1,i2-2*is2,i3,ex)
                      u1(i1-2*is1,i2-2*is2,i3,ey)=wk1(i1-2*is1,i2-2*is2,i3,ey)
                      u2(j1-2*js1,j2-2*js2,j3,ex)=wk2(j1-2*js1,j2-2*js2,j3,ex)
                      u2(j1-2*js1,j2-2*js2,j3,ey)=wk2(j1-2*js1,j2-2*js2,j3,ey)
                      u1(i1-  is1,i2-  is2,i3,hz)=wk1(i1-  is1,i2-  is2,i3,hz)
                      u2(j1-  js1,j2-  js2,j3,hz)=wk2(j1-  js1,j2-  js2,j3,hz)
                      u1(i1-2*is1,i2-2*is2,i3,hz)=wk1(i1-2*is1,i2-2*is2,i3,hz)
                      u2(j1-2*js1,j2-2*js2,j3,hz)=wk2(j1-2*js1,j2-2*js2,j3,hz)
                       end if
                       j1=j1+1
                      end do
                      j2=j2+1
                     end do
                  end if
             end if  
             ! fixup corner points 
             ! if( .false. )then
             !   fixupInterfaceEndValues(2,rectangular,u1,side1,axis1,axis1p1,axis1p2,boundaryCondition1,gridIndexRange1,dx1,dr1)
             !   fixupInterfaceEndValues(2,rectangular,u2,side2,axis2,axis2p1,axis2p2,boundaryCondition2,gridIndexRange2,dx2,dr2)
             ! end if
             ! periodic update
             if( parallel.eq.0 )then
              axisp1=mod(axis1+1,nd)
              if( boundaryCondition1(0,axisp1).lt.0 )then
               ! direction axisp1 is periodic
               diff(axis1)=0
               diff(axisp1)=gridIndexRange1(1,axisp1)-gridIndexRange1(0,axisp1)
               if( side1.eq.0 )then
                 ! assign 4 ghost points outside lower corner
                 np1a=gridIndexRange1(0,0)-2
                 np1b=gridIndexRange1(0,0)-1
                 np2a=gridIndexRange1(0,1)-2
                 np2b=gridIndexRange1(0,1)-1
                 do i3=n3a,n3b
                 do i2=np2a,np2b
                 do i1=np1a,np1b
                 do n=ex,hz
                   ! write(*,'(" periodic i1,i2,i3,n=",4i4)') i1,i2,i3,n
                  u1(i1,i2,i3,n) = u1(i1+diff(0),i2+diff(1),i3,n)
                 end do
                 end do
                 end do
                 end do
                 ! assign 4 ghost points outside upper corner
                 if( axis1.eq.0 )then
                   np2a=gridIndexRange1(1,axisp1)+1
                   np2b=gridIndexRange1(1,axisp1)+2
                 else
                   np1a=gridIndexRange1(1,axisp1)+1
                   np1b=gridIndexRange1(1,axisp1)+2
                 end if
                 do i3=n3a,n3b
                 do i2=np2a,np2b
                 do i1=np1a,np1b
                 do n=ex,hz
                   ! write(*,'(" periodic i1,i2,i3,n=",4i4)') i1,i2,i3,n
                  u1(i1,i2,i3,n) = u1(i1-diff(0),i2-diff(1),i3,n)
                 end do
                 end do
                 end do
                 end do
               else
                 ! assign 4 ghost points outside upper corner
                 np1a=gridIndexRange1(1,0)+1
                 np1b=gridIndexRange1(1,0)+2
                 np2a=gridIndexRange1(1,1)+1
                 np2b=gridIndexRange1(1,1)+2
                 do i3=n3a,n3b
                 do i2=np2a,np2b
                 do i1=np1a,np1b
                 do n=ex,hz
                   ! write(*,'(" periodic i1,i2,i3,n=",4i4)') i1,i2,i3,n
                  u1(i1,i2,i3,n) = u1(i1-diff(0),i2-diff(1),i3,n)
                 end do
                 end do
                 end do
                 end do
                 if( axis1.eq.0 )then
                   np2a=gridIndexRange1(0,axisp1)-2
                   np2b=gridIndexRange1(0,axisp1)-1
                 else
                   np1a=gridIndexRange1(0,axisp1)-2
                   np1b=gridIndexRange1(0,axisp1)-1
                 end if
                 do i3=n3a,n3b
                 do i2=np2a,np2b
                 do i1=np1a,np1b
                 do n=ex,hz
                   ! write(*,'(" periodic i1,i2,i3,n=",4i4)') i1,i2,i3,n
                  u1(i1,i2,i3,n) = u1(i1+diff(0),i2+diff(1),i3,n)
                 end do
                 end do
                 end do
                 end do
               end if
              endif
             end if
             if( parallel.eq.0 )then
              axisp1=mod(axis2+1,nd)
              if( boundaryCondition2(0,axisp1).lt.0 )then
               ! direction axisp1 is periodic
               diff(axis2)=0
               diff(axisp1)=gridIndexRange2(1,axisp1)-gridIndexRange2(0,axisp1)
               if( side2.eq.0 )then
                 ! assign 4 ghost points outside lower corner
                 np1a=gridIndexRange2(0,0)-2
                 np1b=gridIndexRange2(0,0)-1
                 np2a=gridIndexRange2(0,1)-2
                 np2b=gridIndexRange2(0,1)-1
                 do i3=n3a,n3b
                 do i2=np2a,np2b
                 do i1=np1a,np1b
                 do n=ex,hz
                   ! write(*,'(" periodic i1,i2,i3,n=",4i4)') i1,i2,i3,n
                  u2(i1,i2,i3,n) = u2(i1+diff(0),i2+diff(1),i3,n)
                 end do
                 end do
                 end do
                 end do
                 ! assign 4 ghost points outside upper corner
                 if( axis2.eq.0 )then
                   np2a=gridIndexRange2(1,axisp1)+1
                   np2b=gridIndexRange2(1,axisp1)+2
                 else
                   np1a=gridIndexRange2(1,axisp1)+1
                   np1b=gridIndexRange2(1,axisp1)+2
                 end if
                 do i3=n3a,n3b
                 do i2=np2a,np2b
                 do i1=np1a,np1b
                 do n=ex,hz
                   ! write(*,'(" periodic i1,i2,i3,n=",4i4)') i1,i2,i3,n
                  u2(i1,i2,i3,n) = u2(i1-diff(0),i2-diff(1),i3,n)
                 end do
                 end do
                 end do
                 end do
               else
                 ! assign 4 ghost points outside upper corner
                 np1a=gridIndexRange2(1,0)+1
                 np1b=gridIndexRange2(1,0)+2
                 np2a=gridIndexRange2(1,1)+1
                 np2b=gridIndexRange2(1,1)+2
                 do i3=n3a,n3b
                 do i2=np2a,np2b
                 do i1=np1a,np1b
                 do n=ex,hz
                   ! write(*,'(" periodic i1,i2,i3,n=",4i4)') i1,i2,i3,n
                  u2(i1,i2,i3,n) = u2(i1-diff(0),i2-diff(1),i3,n)
                 end do
                 end do
                 end do
                 end do
                 if( axis2.eq.0 )then
                   np2a=gridIndexRange2(0,axisp1)-2
                   np2b=gridIndexRange2(0,axisp1)-1
                 else
                   np1a=gridIndexRange2(0,axisp1)-2
                   np1b=gridIndexRange2(0,axisp1)-1
                 end if
                 do i3=n3a,n3b
                 do i2=np2a,np2b
                 do i1=np1a,np1b
                 do n=ex,hz
                   ! write(*,'(" periodic i1,i2,i3,n=",4i4)') i1,i2,i3,n
                  u2(i1,i2,i3,n) = u2(i1+diff(0),i2+diff(1),i3,n)
                 end do
                 end do
                 end do
                 end do
               end if
              endif
             end if
           end if ! assign ghost
         ! End rectangular, 2D, order 4
         ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        else if( nd.eq.2 .and. orderOfAccuracy.eq.4 .and. gridType.eq.curvilinear )then
         ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         ! End curvilinear, 2D, order 4
         ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        else if( nd.eq.3 .and. orderOfAccuracy.eq.2 .and. gridType.eq.curvilinear )then
         ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         ! End rectangular, 3D, order 2
         ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        else if( nd.eq.3 .and. orderOfAccuracy.eq.4 .and. gridType.eq.curvilinear )then
          ! called elsewhere 
          stop 4562
      !-         call interfaceOpt3dOrder4( nd, nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,!-                               gridIndexRange1, u1,u1n,u1m,wk1, mask1,rsxy1, xy1, p1,p1n,p1m, boundaryCondition1, !-                               md1a,md1b,md2a,md2b,md3a,md3b,!-                               gridIndexRange2, u2,u2n,u2m,wk2, mask2,rsxy2, xy2, p2,p2n,p2m, boundaryCondition2, !-                               ipar, rpar, !-                               aa2,aa4,aa8, ipvt2,ipvt4,ipvt8, !-                               ierr )
        else
          if( nd.eq.3 .and. gridType.eq.rectangular )then
            write(*,'("interface3d: ERROR: 3d rectangular not implemented for interfaces")') 
          else 
            write(debugFile,'("interface3d: ERROR: unknown options nd,order=",2i3)') nd,orderOfAccuracy
            write(*,'("interface3d: ERROR: unknown options nd,order,gridType=",3i3)') nd,orderOfAccuracy,gridType
          end if
          stop 3214
        end if
       return
       end
