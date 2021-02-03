! This file automatically generated from interfaceOpt.bf90 with bpp.
 subroutine interfaceMxMLA2dOrder2c( nd, nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,gridIndexRange1, u1,u1n,u1m, wk1, mask1,rsxy1, xy1, p1,p1n,p1m, q1,q1n,q1m, boundaryCondition1, md1a,md1b,md2a,md2b,md3a,md3b,gridIndexRange2, u2,u2n,u2m, wk2, mask2,rsxy2, xy2, p2,p2n,p2m, q2,q2n,q2m, boundaryCondition2, ipar, rpar, aa2,aa4,aa8, ipvt2,ipvt4,ipvt8, ierr )
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
       integer useJacobiUpdate
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
         write(*,'("  ... useImpedanceInterfaceProjection=",i2)') useImpedanceInterfaceProjection
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
       end if
       useNonlinearModel=0
       if( nonlinearModel1 .ne. noNonlinearModel )then
         write(*,'("--interface3d-- nonlinearModel1=",i4," (1=multilevelAtomic)")') nonlinearModel1
         useNonlinearModel=1
         call getMultilevelAtomicParameters( grid1, nlPar1, maxPar, maxPar, numPolar1, numberOfAtomicLevels1 )
         if( numPolar1.ne.numberOfPolarizationVectors1 )then
           write(*,'(" interface3d:ERROR: numberOfPolarizationVectors1 does not match numPolar1 from nonlinear model!!")')
           stop 8888
         end if
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
       if( nonlinearModel2 .ne. noNonlinearModel )then
         useNonlinearModel=1
         write(*,'("--interface3d-- nonlinearModel2=",i4," (1=multilevelAtomic)")') nonlinearModel2
         call getMultilevelAtomicParameters( grid2, nlPar2, maxPar, maxPar, numPolar2, numberOfAtomicLevels2 )
         if( numPolar2.ne.numberOfPolarizationVectors2 )then
           write(*,'(" interface3d:ERROR: numberOfPolarizationVectors2 does not match numPolar2 from nonlinear model!!")')
           stop 9999
         end if
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
          ! *******************************
          ! ***** 2d curvilinear case *****
          ! *******************************
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
                       an1=rsxy1(i1,i2,i3,axis1,0)   ! normal (an1,an2)
                       an2=rsxy1(i1,i2,i3,axis1,1)
                       aNorm=max(epsx,sqrt(an1**2+an2**2))
                       an1=an1/aNorm
                       an2=an2/aNorm
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
                       an1=rsxy1(i1,i2,i3,axis1,0)   ! normal (an1,an2)
                       an2=rsxy1(i1,i2,i3,axis1,1)
                       aNorm=max(epsx,sqrt(an1**2+an2**2))
                       an1=an1/aNorm
                       an2=an2/aNorm
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
                     an1=rsxy1(i1,i2,i3,axis1,0)   ! normal (an1,an2)
                     an2=rsxy1(i1,i2,i3,axis1,1)
                     aNorm=max(epsx,sqrt(an1**2+an2**2))
                     an1=an1/aNorm
                     an2=an2/aNorm
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
                     an1=rsxy1(i1,i2,i3,axis1,0)
                     an2=rsxy1(i1,i2,i3,axis1,1)
                     aNorm=max(epsx,sqrt(an1**2+an2**2))
                     an1=an1/aNorm
                     an2=an2/aNorm
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
           ! ----- assign ghost using jump conditions -----
           if( assignInterfaceGhostValues.eq.1 )then
            ! initialization step: assign first ghost line by extrapolation
              if( solveForE.ne.0 )then
                i3=n3a
                j3=m3a
                j2=mm2a
                do i2=nn2a,nn2b
                 j1=mm1a
                 do i1=nn1a,nn1b
                 u1(i1-is1,i2-is2,i3,ex)=(3.*u1(i1,i2,i3,ex)-3.*u1(i1+is1,i2+is2,i3+is3,ex)+u1(i1+2*is1,i2+2*is2,i3+2*is3,ex))
                 u1(i1-is1,i2-is2,i3,ey)=(3.*u1(i1,i2,i3,ey)-3.*u1(i1+is1,i2+is2,i3+is3,ey)+u1(i1+2*is1,i2+2*is2,i3+2*is3,ey))
                 u1(i1-is1,i2-is2,i3,hz)=(3.*u1(i1,i2,i3,hz)-3.*u1(i1+is1,i2+is2,i3+is3,hz)+u1(i1+2*is1,i2+2*is2,i3+2*is3,hz))
                 u2(j1-js1,j2-js2,j3,ex)=(3.*u2(j1,j2,j3,ex)-3.*u2(j1+js1,j2+js2,j3+js3,ex)+u2(j1+2*js1,j2+2*js2,j3+2*js3,ex))
                 u2(j1-js1,j2-js2,j3,ey)=(3.*u2(j1,j2,j3,ey)-3.*u2(j1+js1,j2+js2,j3+js3,ey)+u2(j1+2*js1,j2+2*js2,j3+2*js3,ey))
                 u2(j1-js1,j2-js2,j3,hz)=(3.*u2(j1,j2,j3,hz)-3.*u2(j1+js1,j2+js2,j3+js3,hz)+u2(j1+2*js1,j2+2*js2,j3+2*js3,hz))
                 if( dispersive.ne.noDispersion )then
                  do jv=0,numberOfPolarizationVectors1-1
                    do n=0,nd-1
                      pc = n + jv*nd 
                      p1(i1-is1,i2-is2,i3,pc)=(3.*p1(i1,i2,i3,pc)-3.*p1(i1+is1,i2+is2,i3+is3,pc)+p1(i1+2*is1,i2+2*is2,i3+2*is3,pc))
                    end do
                  end do
                  do jv=0,numberOfPolarizationVectors2-1
                    do n=0,nd-1
                      pc = n + jv*nd 
                      p2(j1-js1,j2-js2,j3,pc)=(3.*p2(j1,j2,j3,pc)-3.*p2(j1+js1,j2+js2,j3+js3,pc)+p2(j1+2*js1,j2+2*js2,j3+2*js3,pc))
                    end do
                  end do
                 end if
                if (nonlinearModel1 .ne. noNonlinearModel) then
                  do na = 0,numberOfAtomicLevels1-1
                    q1(i1  -is1,i2  -is2,i3,na)=(3.*q1(i1,i2,i3,na)-3.*q1(i1+is1,i2+is2,i3+is3,na)+q1(i1+2*is1,i2+2*is2,i3+2*is3,na))
                    q1(i1-2*is1,i2-2*is2,i3,na)=(3.*q1(i1-is1,i2-is2,i3,na)-3.*q1(i1-is1+is1,i2-is2+is2,i3+is3,na)+q1(i1-is1+2*is1,i2-is2+2*is2,i3+2*is3,na))
                  enddo
                end if
                if (nonlinearModel2 .ne. noNonlinearModel) then
                  do na = 0,numberOfAtomicLevels2-1
                    q2(j1  -js1,j2  -js2,j3,na)=(3.*q2(j1,j2,j3,na)-3.*q2(j1+js1,j2+js2,j3+js3,na)+q2(j1+2*js1,j2+2*js2,j3+2*js3,na))
                    q2(j1-2*js1,j2-2*js2,j3,na)=(3.*q2(j1-js1,j2-js2,j3,na)-3.*q2(j1-js1+js1,j2-js2+js2,j3+js3,na)+q2(j1-js1+2*js1,j2-js2+2*js2,j3+2*js3,na))
                  enddo
                endif
                  j1=j1+1
                 end do
                 j2=j2+1
                end do
              end if
            ! Macro to assign ghost values:
            if( dispersive.eq.0 )then
            else if( useNonlinearModel.eq.0 )then
              ! dispersive case
            else
              ! nonlinear model
                  ! ****************************************************
                  ! ***********  2D, ORDER=2, CURVILINEAR **************
                  ! ****************************************************
                if( t.le.3.*dt .and. debug.gt.0 )then
                  write(*,'("Interface>>>","22curvilinear-nonlinear-MLA")')
                end if
                ! --- initialize some forcing functions ---
                do n=0,nd-1
                  fev1(n)=0.
                  fev2(n)=0.
                  do jv=0,numberOfPolarizationVectors1-1
                    fpv1(n,jv)=0.
                  end do
                  do jv=0,numberOfPolarizationVectors2-1
                    fpv2(n,jv)=0.
                  end do
                end do
                ! forcing functions for N
                do jv = 0,numberOfAtomicLevels1-1
                    fnv1(jv) = 0.
                    fntv1(jv) = 0.
                enddo
                do jv = 0,numberOfAtomicLevels2-1
                    fnv2(jv) = 0.
                    fntv2(jv) = 0.
                enddo
                ! ----------------- START LOOP OVER INTERFACE -------------------------
                 i3=n3a
                 j3=m3a
                 j2=m2a
                 do i2=n2a,n2b
                  j1=m1a
                  do i1=n1a,n1b
                   if( mask1(i1,i2,i3).gt.0 .and. mask2(j1,j2,j3).gt.0 )then
                  ! here is the normal (assumed to be the same on both sides)
                  an1=rsxy1(i1,i2,i3,axis1,0)   ! normal (an1,an2)
                  an2=rsxy1(i1,i2,i3,axis1,1)
                  aNorm=max(epsx,sqrt(an1**2+an2**2))
                  an1=an1/aNorm
                  an2=an2/aNorm
                  tau1=-an2
                  tau2= an1
                  ! first evaluate the equations we want to solve with the wrong values at the ghost points:
                   ! NOTE: the jacobian derivatives can be computed once for all components
                    ! this next call will define the jacobian and its derivatives (parameteric and spatial)
                    aj1rx = rsxy1(i1,i2,i3,0,0)
                    aj1rxr = (-rsxy1(i1-1,i2,i3,0,0)+rsxy1(i1+1,i2,i3,0,0))/(2.*dr1(0))
                    aj1rxs = (-rsxy1(i1,i2-1,i3,0,0)+rsxy1(i1,i2+1,i3,0,0))/(2.*dr1(1))
                    aj1sx = rsxy1(i1,i2,i3,1,0)
                    aj1sxr = (-rsxy1(i1-1,i2,i3,1,0)+rsxy1(i1+1,i2,i3,1,0))/(2.*dr1(0))
                    aj1sxs = (-rsxy1(i1,i2-1,i3,1,0)+rsxy1(i1,i2+1,i3,1,0))/(2.*dr1(1))
                    aj1ry = rsxy1(i1,i2,i3,0,1)
                    aj1ryr = (-rsxy1(i1-1,i2,i3,0,1)+rsxy1(i1+1,i2,i3,0,1))/(2.*dr1(0))
                    aj1rys = (-rsxy1(i1,i2-1,i3,0,1)+rsxy1(i1,i2+1,i3,0,1))/(2.*dr1(1))
                    aj1sy = rsxy1(i1,i2,i3,1,1)
                    aj1syr = (-rsxy1(i1-1,i2,i3,1,1)+rsxy1(i1+1,i2,i3,1,1))/(2.*dr1(0))
                    aj1sys = (-rsxy1(i1,i2-1,i3,1,1)+rsxy1(i1,i2+1,i3,1,1))/(2.*dr1(1))
                    aj1rxx = aj1rx*aj1rxr+aj1sx*aj1rxs
                    aj1rxy = aj1ry*aj1rxr+aj1sy*aj1rxs
                    aj1sxx = aj1rx*aj1sxr+aj1sx*aj1sxs
                    aj1sxy = aj1ry*aj1sxr+aj1sy*aj1sxs
                    aj1ryx = aj1rx*aj1ryr+aj1sx*aj1rys
                    aj1ryy = aj1ry*aj1ryr+aj1sy*aj1rys
                    aj1syx = aj1rx*aj1syr+aj1sx*aj1sys
                    aj1syy = aj1ry*aj1syr+aj1sy*aj1sys
                     uu1 = u1(i1,i2,i3,ex)
                     uu1r = (-u1(i1-1,i2,i3,ex)+u1(i1+1,i2,i3,ex))/(2.*dr1(0))
                     uu1s = (-u1(i1,i2-1,i3,ex)+u1(i1,i2+1,i3,ex))/(2.*dr1(1))
                     uu1rr = (u1(i1-1,i2,i3,ex)-2.*u1(i1,i2,i3,ex)+u1(i1+1,i2,i3,ex))/(dr1(0)**2)
                     uu1rs = (-(-u1(i1-1,i2-1,i3,ex)+u1(i1-1,i2+1,i3,ex))/(2.*dr1(1))+(-u1(i1+1,i2-1,i3,ex)+u1(i1+1,i2+1,i3,ex))/(2.*dr1(1)))/(2.*dr1(0))
                     uu1ss = (u1(i1,i2-1,i3,ex)-2.*u1(i1,i2,i3,ex)+u1(i1,i2+1,i3,ex))/(dr1(1)**2)
                      u1x = aj1rx*uu1r+aj1sx*uu1s
                      u1y = aj1ry*uu1r+aj1sy*uu1s
                      t1 = aj1rx**2
                      t6 = aj1sx**2
                      u1xx = t1*uu1rr+2*aj1rx*aj1sx*uu1rs+t6*uu1ss+aj1rxx*uu1r+aj1sxx*uu1s
                      t1 = aj1ry**2
                      t6 = aj1sy**2
                      u1yy = t1*uu1rr+2*aj1ry*aj1sy*uu1rs+t6*uu1ss+aj1ryy*uu1r+aj1syy*uu1s
                    u1Lap = u1xx+ u1yy
                     vv1 = u1(i1,i2,i3,ey)
                     vv1r = (-u1(i1-1,i2,i3,ey)+u1(i1+1,i2,i3,ey))/(2.*dr1(0))
                     vv1s = (-u1(i1,i2-1,i3,ey)+u1(i1,i2+1,i3,ey))/(2.*dr1(1))
                     vv1rr = (u1(i1-1,i2,i3,ey)-2.*u1(i1,i2,i3,ey)+u1(i1+1,i2,i3,ey))/(dr1(0)**2)
                     vv1rs = (-(-u1(i1-1,i2-1,i3,ey)+u1(i1-1,i2+1,i3,ey))/(2.*dr1(1))+(-u1(i1+1,i2-1,i3,ey)+u1(i1+1,i2+1,i3,ey))/(2.*dr1(1)))/(2.*dr1(0))
                     vv1ss = (u1(i1,i2-1,i3,ey)-2.*u1(i1,i2,i3,ey)+u1(i1,i2+1,i3,ey))/(dr1(1)**2)
                      v1x = aj1rx*vv1r+aj1sx*vv1s
                      v1y = aj1ry*vv1r+aj1sy*vv1s
                      t1 = aj1rx**2
                      t6 = aj1sx**2
                      v1xx = t1*vv1rr+2*aj1rx*aj1sx*vv1rs+t6*vv1ss+aj1rxx*vv1r+aj1sxx*vv1s
                      t1 = aj1ry**2
                      t6 = aj1sy**2
                      v1yy = t1*vv1rr+2*aj1ry*aj1sy*vv1rs+t6*vv1ss+aj1ryy*vv1r+aj1syy*vv1s
                    v1Lap = v1xx+ v1yy
                   ! NOTE: the jacobian derivatives can be computed once for all components
                    ! this next call will define the jacobian and its derivatives (parameteric and spatial)
                    aj2rx = rsxy2(j1,j2,j3,0,0)
                    aj2rxr = (-rsxy2(j1-1,j2,j3,0,0)+rsxy2(j1+1,j2,j3,0,0))/(2.*dr2(0))
                    aj2rxs = (-rsxy2(j1,j2-1,j3,0,0)+rsxy2(j1,j2+1,j3,0,0))/(2.*dr2(1))
                    aj2sx = rsxy2(j1,j2,j3,1,0)
                    aj2sxr = (-rsxy2(j1-1,j2,j3,1,0)+rsxy2(j1+1,j2,j3,1,0))/(2.*dr2(0))
                    aj2sxs = (-rsxy2(j1,j2-1,j3,1,0)+rsxy2(j1,j2+1,j3,1,0))/(2.*dr2(1))
                    aj2ry = rsxy2(j1,j2,j3,0,1)
                    aj2ryr = (-rsxy2(j1-1,j2,j3,0,1)+rsxy2(j1+1,j2,j3,0,1))/(2.*dr2(0))
                    aj2rys = (-rsxy2(j1,j2-1,j3,0,1)+rsxy2(j1,j2+1,j3,0,1))/(2.*dr2(1))
                    aj2sy = rsxy2(j1,j2,j3,1,1)
                    aj2syr = (-rsxy2(j1-1,j2,j3,1,1)+rsxy2(j1+1,j2,j3,1,1))/(2.*dr2(0))
                    aj2sys = (-rsxy2(j1,j2-1,j3,1,1)+rsxy2(j1,j2+1,j3,1,1))/(2.*dr2(1))
                    aj2rxx = aj2rx*aj2rxr+aj2sx*aj2rxs
                    aj2rxy = aj2ry*aj2rxr+aj2sy*aj2rxs
                    aj2sxx = aj2rx*aj2sxr+aj2sx*aj2sxs
                    aj2sxy = aj2ry*aj2sxr+aj2sy*aj2sxs
                    aj2ryx = aj2rx*aj2ryr+aj2sx*aj2rys
                    aj2ryy = aj2ry*aj2ryr+aj2sy*aj2rys
                    aj2syx = aj2rx*aj2syr+aj2sx*aj2sys
                    aj2syy = aj2ry*aj2syr+aj2sy*aj2sys
                     uu2 = u2(j1,j2,j3,ex)
                     uu2r = (-u2(j1-1,j2,j3,ex)+u2(j1+1,j2,j3,ex))/(2.*dr2(0))
                     uu2s = (-u2(j1,j2-1,j3,ex)+u2(j1,j2+1,j3,ex))/(2.*dr2(1))
                     uu2rr = (u2(j1-1,j2,j3,ex)-2.*u2(j1,j2,j3,ex)+u2(j1+1,j2,j3,ex))/(dr2(0)**2)
                     uu2rs = (-(-u2(j1-1,j2-1,j3,ex)+u2(j1-1,j2+1,j3,ex))/(2.*dr2(1))+(-u2(j1+1,j2-1,j3,ex)+u2(j1+1,j2+1,j3,ex))/(2.*dr2(1)))/(2.*dr2(0))
                     uu2ss = (u2(j1,j2-1,j3,ex)-2.*u2(j1,j2,j3,ex)+u2(j1,j2+1,j3,ex))/(dr2(1)**2)
                      u2x = aj2rx*uu2r+aj2sx*uu2s
                      u2y = aj2ry*uu2r+aj2sy*uu2s
                      t1 = aj2rx**2
                      t6 = aj2sx**2
                      u2xx = t1*uu2rr+2*aj2rx*aj2sx*uu2rs+t6*uu2ss+aj2rxx*uu2r+aj2sxx*uu2s
                      t1 = aj2ry**2
                      t6 = aj2sy**2
                      u2yy = t1*uu2rr+2*aj2ry*aj2sy*uu2rs+t6*uu2ss+aj2ryy*uu2r+aj2syy*uu2s
                    u2Lap = u2xx+ u2yy
                     vv2 = u2(j1,j2,j3,ey)
                     vv2r = (-u2(j1-1,j2,j3,ey)+u2(j1+1,j2,j3,ey))/(2.*dr2(0))
                     vv2s = (-u2(j1,j2-1,j3,ey)+u2(j1,j2+1,j3,ey))/(2.*dr2(1))
                     vv2rr = (u2(j1-1,j2,j3,ey)-2.*u2(j1,j2,j3,ey)+u2(j1+1,j2,j3,ey))/(dr2(0)**2)
                     vv2rs = (-(-u2(j1-1,j2-1,j3,ey)+u2(j1-1,j2+1,j3,ey))/(2.*dr2(1))+(-u2(j1+1,j2-1,j3,ey)+u2(j1+1,j2+1,j3,ey))/(2.*dr2(1)))/(2.*dr2(0))
                     vv2ss = (u2(j1,j2-1,j3,ey)-2.*u2(j1,j2,j3,ey)+u2(j1,j2+1,j3,ey))/(dr2(1)**2)
                      v2x = aj2rx*vv2r+aj2sx*vv2s
                      v2y = aj2ry*vv2r+aj2sy*vv2s
                      t1 = aj2rx**2
                      t6 = aj2sx**2
                      v2xx = t1*vv2rr+2*aj2rx*aj2sx*vv2rs+t6*vv2ss+aj2rxx*vv2r+aj2sxx*vv2s
                      t1 = aj2ry**2
                      t6 = aj2sy**2
                      v2yy = t1*vv2rr+2*aj2ry*aj2sy*vv2rs+t6*vv2ss+aj2ryy*vv2r+aj2syy*vv2s
                    v2Lap = v2xx+ v2yy
                  ! if( .true. .or. debug.gt.4 )then 
                  !    write(*,'(" START  (i1,i2)=",2i3," v=Ey(-1:1,-1:1)",9e14.6)') i1,i2, ((u1(i1+k1,i2+k2,i3,ey),k1=-1,1),k2=-1,1)
                  !    write(*,'(" START    mu1,mu2=",2e10.2," v1y,u1x,v2y,u2x=",4e10.2)') mu1,mu2,v1y,u1x,v2y,u2x
                  !  end if
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
                          call ogderiv(ep, 2,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pxc+jv*nd+n, pett(n)   )
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
                          call ogderiv(ep, 2,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pxc+jv*nd+n, pett(n)   )
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
                  ! Evaulate RHS, f(n),n=0,1,2,3 using current ghost values: 
                   f(0)=(u1x+v1y) - (u2x+v2y)
                   f(1)=( an1*u1Lap +an2*v1Lap )/mu1 - ( an1*u2Lap +an2*v2Lap )/mu2 
                   f(2)=(v1x-u1y)/mu1 - (v2x-u2y)/mu2
                   f(3)=( ( tau1*u1Lap +tau2*v1Lap )*beta1/epsmu1 - alphaP1*(tau1*fp1(0)+tau2*fp1(1)) ) - ( ( tau1*u2Lap +tau2*v2Lap )*beta2/epsmu2 - alphaP2*(tau1*fp2(0)+tau2*fp2(1)) )
                   if( twilightZone.eq.1 )then
                     ! For now we assume mu1=mu2 and TZ solutions are the same on both sides.
                     ! f(3) = [ tv.E.tt] = [ tv.( c^2*Delta(E) - alphaP*P.tt + fev) ] 
                     !      = [ tv.( c^2*Delta(E) - alphaP*P.tt] + [ tv.fev ]
                     ! -- add on the jump in the forcing ---
                     f(3) = f(3) + ( tau1*(fev1(0)-fev2(0)) + tau2*(fev1(1)-fev2(1)) )
                    !-    ! f(3) = [ tv.( c^2*Delta(E) - alphaP*P.tt ] - [ tv.( c^2*Delta(E^e) - alphaP*(P^e).tt ] =0 
                    !-    call ogderiv(ep, 0,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, uexx )
                    !-    call ogderiv(ep, 0,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, ueyy )
                    !-    call ogderiv(ep, 0,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, vexx )
                    !-    call ogderiv(ep, 0,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, veyy )
                    !- 
                    !-    ueLap = uexx + ueyy
                    !-    veLap = vexx + veyy
                    !-    f(3) = f(3) - ( tau1*ueLap +tau2*veLap )*(1./epsmu1-1./epsmu2)
                    !- 
                    !-    f(3) = f(3) +  alphaP1*( tau1*pettSum1(0) + tau2*pettSum1(1) ) - !-                   alphaP2*( tau1*pettSum2(0) + tau2*pettSum2(1) )
                    !-    end if 
                     ! write(debugFile,'(" u1Lap,ueLap=",2e10.2," v1Lap,veLap=",2e10.2)') u1Lap,ueLap,v1Lap,veLap
                     ! write(debugFile,'(" u1Lap,ueLap=",2e10.2," v1Lap,veLap=",2e10.2)') u1Lap,ueLap,v1Lap,veLap
                   end if
                  ! write(debugFile,'(" --> order2-curv: i1,i2=",2i4," f(start)=",4f8.3)') i1,i2,f(0),f(1),f(2),f(3)
                  ! write(debugFile,'(" --> u1(ghost),u1=",4f8.3)') u1(i1-is1,i2-is2,i3,ex),u1(i1,i2,i3,ex)
                  ! write(debugFile,'(" --> u2(ghost),u2=",4f8.3)') u2(j1-js1,j2-js2,j3,ex),u2(j1,j2,j3,ex)
                  ! '
                  ! here is the matrix of coefficients for the unknowns u1(-1),v1(-1),u2(-1),v2(-1)
                  ! Solve:
                  !     
                  !       A [ U ] = A [ U(old) ] - [ f ]
                  ! ---- EQUATION 0 ----- 
                  a4(0,0) = -is*rsxy1(i1,i2,i3,axis1,0)/(2.*dr1(axis1))    ! coeff of u1(-1) from [u.x+v.y] 
                  a4(0,1) = -is*rsxy1(i1,i2,i3,axis1,1)/(2.*dr1(axis1))    ! coeff of v1(-1) from [u.x+v.y] 
                  a4(0,2) =  js*rsxy2(j1,j2,j3,axis2,0)/(2.*dr2(axis2))    ! coeff of u2(-1) from [u.x+v.y] 
                  a4(0,3) =  js*rsxy2(j1,j2,j3,axis2,1)/(2.*dr2(axis2))    ! coeff of v2(-1) from [u.x+v.y] 
                  ! ---- EQUATION 2 ----- 
                  a4(2,0) =  is*rsxy1(i1,i2,i3,axis1,1)/(2.*dr1(axis1))/mu1   ! coeff of u1(-1) from [(v.x - u.y)/mu] 
                  a4(2,1) = -is*rsxy1(i1,i2,i3,axis1,0)/(2.*dr1(axis1))/mu1   ! coeff of v1(-1) from [(v.x - u.y)/mu] 
                  a4(2,2) = -js*rsxy2(j1,j2,j3,axis2,1)/(2.*dr2(axis2))/mu2   ! coeff of u2(-1) from [(v.x - u.y)/mu] 
                  a4(2,3) =  js*rsxy2(j1,j2,j3,axis2,0)/(2.*dr2(axis2))/mu2   ! coeff of v2(-1) from [(v.x - u.y)/mu] 
                  ! coeff of u(-1) from lap = u.xx + u.yy
                  rxx1(0,0,0)=aj1rxx
                  rxx1(1,0,0)=aj1sxx
                  rxx1(0,1,1)=aj1ryy
                  rxx1(1,1,1)=aj1syy
                  rxx2(0,0,0)=aj2rxx
                  rxx2(1,0,0)=aj2sxx
                  rxx2(0,1,1)=aj2ryy
                  rxx2(1,1,1)=aj2syy
                  ! clap1=(rsxy1(i1,i2,i3,axis1,0)**2+rsxy1(i1,i2,i3,axis1,1)**2)/(dr1(axis1)**2) !           -is*(rsxy1x22(i1,i2,i3,axis1,0)+rsxy1y22(i1,i2,i3,axis1,1))/(2.*dr1(axis1))
                  ! clap2=(rsxy2(j1,j2,j3,axis2,0)**2+rsxy2(j1,j2,j3,axis2,1)**2)/(dr2(axis2)**2) !             -js*(rsxy2x22(j1,j2,j3,axis2,0)+rsxy2y22(j1,j2,j3,axis2,1))/(2.*dr2(axis2)) 
                  clap1=(rsxy1(i1,i2,i3,axis1,0)**2+rsxy1(i1,i2,i3,axis1,1)**2)/(dr1(axis1)**2) -is*(rxx1(axis1,0,0)+rxx1(axis1,1,1))/(2.*dr1(axis1))
                  clap2=(rsxy2(j1,j2,j3,axis2,0)**2+rsxy2(j1,j2,j3,axis2,1)**2)/(dr2(axis2)**2) -js*(rxx2(axis2,0,0)+rxx2(axis2,1,1))/(2.*dr2(axis2)) 
                  ! ---- EQUATION 1 ----- 
                  !   [ n.(uv.xx + u.yy)/mu ] = 0
                  a4(1,0) = an1*clap1/mu1
                  a4(1,1) = an2*clap1/mu1
                  a4(1,2) =-an1*clap2/mu2
                  a4(1,3) =-an2*clap2/mu2
                  ! ---- EQUATION 3 ----- 
                  !   [ tau.(uv.xx+uv.yy)*beta/(eps*mu) + ... ] = 0
                  a4(3,0) = tau1*clap1*beta1/epsmu1
                  a4(3,1) = tau2*clap1*beta1/epsmu1
                  a4(3,2) =-tau1*clap2*beta2/epsmu2
                  a4(3,3) =-tau2*clap2*beta2/epsmu2
                   if( .false. .or. debug.gt.4 )then 
                     write(*,'("BEFORE: --> i1,i2=",2i4," j1,j2=",2i4," f()=",4e10.2)') i1,i2,j1,j2,f(0),f(1),f(2),f(3)
                     write(*,'("     beta1,beta2=",2e10.2," fp1=",2e10.2," fp2=",2e10.2)') beta1,beta2,fp1(0),fp1(1),fp2(0),fp2(1)
                     write(*,'("     mu1,mu2=",2e10.2," v1y,u1x,v2y,u2x=",4e10.2)') mu1,mu2,v1y,u1x,v2y,u2x
                   end if
                  q(0) = u1(i1-is1,i2-is2,i3,ex)
                  q(1) = u1(i1-is1,i2-is2,i3,ey)
                  q(2) = u2(j1-js1,j2-js2,j3,ex)
                  q(3) = u2(j1-js1,j2-js2,j3,ey)
                  ! --- check matrix coefficients by delta function approach ----
                  if( checkCoeff.eq.1 )then
                    numberOfEquations=4
                      ! hw1 = half stencil width
                      hw1=orderOfAccuracy/2
                      hw2=hw1
                      if( nd.eq.2 )then
                        hw3=0
                      else
                        hw3=hw1
                      end if
                      write(*,'("CHECK-COEFF: i1,i2,i3=",3i3," hw1,hw2,hw3=",3i2)') i1,i2,i3,hw1,hw2,hw3
                      ! First eval equations with no pertutbation --> save in f0 
                         ! NOTE: the jacobian derivatives can be computed once for all components
                          ! this next call will define the jacobian and its derivatives (parameteric and spatial)
                          aj1rx = rsxy1(i1,i2,i3,0,0)
                          aj1rxr = (-rsxy1(i1-1,i2,i3,0,0)+rsxy1(i1+1,i2,i3,0,0))/(2.*dr1(0))
                          aj1rxs = (-rsxy1(i1,i2-1,i3,0,0)+rsxy1(i1,i2+1,i3,0,0))/(2.*dr1(1))
                          aj1sx = rsxy1(i1,i2,i3,1,0)
                          aj1sxr = (-rsxy1(i1-1,i2,i3,1,0)+rsxy1(i1+1,i2,i3,1,0))/(2.*dr1(0))
                          aj1sxs = (-rsxy1(i1,i2-1,i3,1,0)+rsxy1(i1,i2+1,i3,1,0))/(2.*dr1(1))
                          aj1ry = rsxy1(i1,i2,i3,0,1)
                          aj1ryr = (-rsxy1(i1-1,i2,i3,0,1)+rsxy1(i1+1,i2,i3,0,1))/(2.*dr1(0))
                          aj1rys = (-rsxy1(i1,i2-1,i3,0,1)+rsxy1(i1,i2+1,i3,0,1))/(2.*dr1(1))
                          aj1sy = rsxy1(i1,i2,i3,1,1)
                          aj1syr = (-rsxy1(i1-1,i2,i3,1,1)+rsxy1(i1+1,i2,i3,1,1))/(2.*dr1(0))
                          aj1sys = (-rsxy1(i1,i2-1,i3,1,1)+rsxy1(i1,i2+1,i3,1,1))/(2.*dr1(1))
                          aj1rxx = aj1rx*aj1rxr+aj1sx*aj1rxs
                          aj1rxy = aj1ry*aj1rxr+aj1sy*aj1rxs
                          aj1sxx = aj1rx*aj1sxr+aj1sx*aj1sxs
                          aj1sxy = aj1ry*aj1sxr+aj1sy*aj1sxs
                          aj1ryx = aj1rx*aj1ryr+aj1sx*aj1rys
                          aj1ryy = aj1ry*aj1ryr+aj1sy*aj1rys
                          aj1syx = aj1rx*aj1syr+aj1sx*aj1sys
                          aj1syy = aj1ry*aj1syr+aj1sy*aj1sys
                           uu1 = u1(i1,i2,i3,ex)
                           uu1r = (-u1(i1-1,i2,i3,ex)+u1(i1+1,i2,i3,ex))/(2.*dr1(0))
                           uu1s = (-u1(i1,i2-1,i3,ex)+u1(i1,i2+1,i3,ex))/(2.*dr1(1))
                           uu1rr = (u1(i1-1,i2,i3,ex)-2.*u1(i1,i2,i3,ex)+u1(i1+1,i2,i3,ex))/(dr1(0)**2)
                           uu1rs = (-(-u1(i1-1,i2-1,i3,ex)+u1(i1-1,i2+1,i3,ex))/(2.*dr1(1))+(-u1(i1+1,i2-1,i3,ex)+u1(i1+1,i2+1,i3,ex))/(2.*dr1(1)))/(2.*dr1(0))
                           uu1ss = (u1(i1,i2-1,i3,ex)-2.*u1(i1,i2,i3,ex)+u1(i1,i2+1,i3,ex))/(dr1(1)**2)
                            u1x = aj1rx*uu1r+aj1sx*uu1s
                            u1y = aj1ry*uu1r+aj1sy*uu1s
                            t1 = aj1rx**2
                            t6 = aj1sx**2
                            u1xx = t1*uu1rr+2*aj1rx*aj1sx*uu1rs+t6*uu1ss+aj1rxx*uu1r+aj1sxx*uu1s
                            t1 = aj1ry**2
                            t6 = aj1sy**2
                            u1yy = t1*uu1rr+2*aj1ry*aj1sy*uu1rs+t6*uu1ss+aj1ryy*uu1r+aj1syy*uu1s
                          u1Lap = u1xx+ u1yy
                           vv1 = u1(i1,i2,i3,ey)
                           vv1r = (-u1(i1-1,i2,i3,ey)+u1(i1+1,i2,i3,ey))/(2.*dr1(0))
                           vv1s = (-u1(i1,i2-1,i3,ey)+u1(i1,i2+1,i3,ey))/(2.*dr1(1))
                           vv1rr = (u1(i1-1,i2,i3,ey)-2.*u1(i1,i2,i3,ey)+u1(i1+1,i2,i3,ey))/(dr1(0)**2)
                           vv1rs = (-(-u1(i1-1,i2-1,i3,ey)+u1(i1-1,i2+1,i3,ey))/(2.*dr1(1))+(-u1(i1+1,i2-1,i3,ey)+u1(i1+1,i2+1,i3,ey))/(2.*dr1(1)))/(2.*dr1(0))
                           vv1ss = (u1(i1,i2-1,i3,ey)-2.*u1(i1,i2,i3,ey)+u1(i1,i2+1,i3,ey))/(dr1(1)**2)
                            v1x = aj1rx*vv1r+aj1sx*vv1s
                            v1y = aj1ry*vv1r+aj1sy*vv1s
                            t1 = aj1rx**2
                            t6 = aj1sx**2
                            v1xx = t1*vv1rr+2*aj1rx*aj1sx*vv1rs+t6*vv1ss+aj1rxx*vv1r+aj1sxx*vv1s
                            t1 = aj1ry**2
                            t6 = aj1sy**2
                            v1yy = t1*vv1rr+2*aj1ry*aj1sy*vv1rs+t6*vv1ss+aj1ryy*vv1r+aj1syy*vv1s
                          v1Lap = v1xx+ v1yy
                         ! NOTE: the jacobian derivatives can be computed once for all components
                          ! this next call will define the jacobian and its derivatives (parameteric and spatial)
                          aj2rx = rsxy2(j1,j2,j3,0,0)
                          aj2rxr = (-rsxy2(j1-1,j2,j3,0,0)+rsxy2(j1+1,j2,j3,0,0))/(2.*dr2(0))
                          aj2rxs = (-rsxy2(j1,j2-1,j3,0,0)+rsxy2(j1,j2+1,j3,0,0))/(2.*dr2(1))
                          aj2sx = rsxy2(j1,j2,j3,1,0)
                          aj2sxr = (-rsxy2(j1-1,j2,j3,1,0)+rsxy2(j1+1,j2,j3,1,0))/(2.*dr2(0))
                          aj2sxs = (-rsxy2(j1,j2-1,j3,1,0)+rsxy2(j1,j2+1,j3,1,0))/(2.*dr2(1))
                          aj2ry = rsxy2(j1,j2,j3,0,1)
                          aj2ryr = (-rsxy2(j1-1,j2,j3,0,1)+rsxy2(j1+1,j2,j3,0,1))/(2.*dr2(0))
                          aj2rys = (-rsxy2(j1,j2-1,j3,0,1)+rsxy2(j1,j2+1,j3,0,1))/(2.*dr2(1))
                          aj2sy = rsxy2(j1,j2,j3,1,1)
                          aj2syr = (-rsxy2(j1-1,j2,j3,1,1)+rsxy2(j1+1,j2,j3,1,1))/(2.*dr2(0))
                          aj2sys = (-rsxy2(j1,j2-1,j3,1,1)+rsxy2(j1,j2+1,j3,1,1))/(2.*dr2(1))
                          aj2rxx = aj2rx*aj2rxr+aj2sx*aj2rxs
                          aj2rxy = aj2ry*aj2rxr+aj2sy*aj2rxs
                          aj2sxx = aj2rx*aj2sxr+aj2sx*aj2sxs
                          aj2sxy = aj2ry*aj2sxr+aj2sy*aj2sxs
                          aj2ryx = aj2rx*aj2ryr+aj2sx*aj2rys
                          aj2ryy = aj2ry*aj2ryr+aj2sy*aj2rys
                          aj2syx = aj2rx*aj2syr+aj2sx*aj2sys
                          aj2syy = aj2ry*aj2syr+aj2sy*aj2sys
                           uu2 = u2(j1,j2,j3,ex)
                           uu2r = (-u2(j1-1,j2,j3,ex)+u2(j1+1,j2,j3,ex))/(2.*dr2(0))
                           uu2s = (-u2(j1,j2-1,j3,ex)+u2(j1,j2+1,j3,ex))/(2.*dr2(1))
                           uu2rr = (u2(j1-1,j2,j3,ex)-2.*u2(j1,j2,j3,ex)+u2(j1+1,j2,j3,ex))/(dr2(0)**2)
                           uu2rs = (-(-u2(j1-1,j2-1,j3,ex)+u2(j1-1,j2+1,j3,ex))/(2.*dr2(1))+(-u2(j1+1,j2-1,j3,ex)+u2(j1+1,j2+1,j3,ex))/(2.*dr2(1)))/(2.*dr2(0))
                           uu2ss = (u2(j1,j2-1,j3,ex)-2.*u2(j1,j2,j3,ex)+u2(j1,j2+1,j3,ex))/(dr2(1)**2)
                            u2x = aj2rx*uu2r+aj2sx*uu2s
                            u2y = aj2ry*uu2r+aj2sy*uu2s
                            t1 = aj2rx**2
                            t6 = aj2sx**2
                            u2xx = t1*uu2rr+2*aj2rx*aj2sx*uu2rs+t6*uu2ss+aj2rxx*uu2r+aj2sxx*uu2s
                            t1 = aj2ry**2
                            t6 = aj2sy**2
                            u2yy = t1*uu2rr+2*aj2ry*aj2sy*uu2rs+t6*uu2ss+aj2ryy*uu2r+aj2syy*uu2s
                          u2Lap = u2xx+ u2yy
                           vv2 = u2(j1,j2,j3,ey)
                           vv2r = (-u2(j1-1,j2,j3,ey)+u2(j1+1,j2,j3,ey))/(2.*dr2(0))
                           vv2s = (-u2(j1,j2-1,j3,ey)+u2(j1,j2+1,j3,ey))/(2.*dr2(1))
                           vv2rr = (u2(j1-1,j2,j3,ey)-2.*u2(j1,j2,j3,ey)+u2(j1+1,j2,j3,ey))/(dr2(0)**2)
                           vv2rs = (-(-u2(j1-1,j2-1,j3,ey)+u2(j1-1,j2+1,j3,ey))/(2.*dr2(1))+(-u2(j1+1,j2-1,j3,ey)+u2(j1+1,j2+1,j3,ey))/(2.*dr2(1)))/(2.*dr2(0))
                           vv2ss = (u2(j1,j2-1,j3,ey)-2.*u2(j1,j2,j3,ey)+u2(j1,j2+1,j3,ey))/(dr2(1)**2)
                            v2x = aj2rx*vv2r+aj2sx*vv2s
                            v2y = aj2ry*vv2r+aj2sy*vv2s
                            t1 = aj2rx**2
                            t6 = aj2sx**2
                            v2xx = t1*vv2rr+2*aj2rx*aj2sx*vv2rs+t6*vv2ss+aj2rxx*vv2r+aj2sxx*vv2s
                            t1 = aj2ry**2
                            t6 = aj2sy**2
                            v2yy = t1*vv2rr+2*aj2ry*aj2sy*vv2rs+t6*vv2ss+aj2ryy*vv2r+aj2syy*vv2s
                          v2Lap = v2xx+ v2yy
                         f(0)=(u1x+v1y) - (u2x+v2y)
                         f(1)=( an1*u1Lap +an2*v1Lap )/mu1 - ( an1*u2Lap +an2*v2Lap )/mu2 
                         f(2)=(v1x-u1y)/mu1 - (v2x-u2y)/mu2
                         f(3)=( ( tau1*u1Lap +tau2*v1Lap )*beta1/epsmu1 - alphaP1*(tau1*fp1(0)+tau2*fp1(1)) ) - ( ( tau1*u2Lap +tau2*v2Lap )*beta2/epsmu2 - alphaP2*(tau1*fp2(0)+tau2*fp2(1)) )
                         if( twilightZone.eq.1 )then
                           ! For now we assume mu1=mu2 and TZ solutions are the same on both sides.
                           ! f(3) = [ tv.E.tt] = [ tv.( c^2*Delta(E) - alphaP*P.tt + fev) ] 
                           !      = [ tv.( c^2*Delta(E) - alphaP*P.tt] + [ tv.fev ]
                           ! -- add on the jump in the forcing ---
                           f(3) = f(3) + ( tau1*(fev1(0)-fev2(0)) + tau2*(fev1(1)-fev2(1)) )
                          !-    ! f(3) = [ tv.( c^2*Delta(E) - alphaP*P.tt ] - [ tv.( c^2*Delta(E^e) - alphaP*(P^e).tt ] =0 
                          !-    call ogderiv(ep, 0,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, uexx )
                          !-    call ogderiv(ep, 0,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, ueyy )
                          !-    call ogderiv(ep, 0,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, vexx )
                          !-    call ogderiv(ep, 0,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, veyy )
                          !- 
                          !-    ueLap = uexx + ueyy
                          !-    veLap = vexx + veyy
                          !-    f(3) = f(3) - ( tau1*ueLap +tau2*veLap )*(1./epsmu1-1./epsmu2)
                          !- 
                          !-    f(3) = f(3) +  alphaP1*( tau1*pettSum1(0) + tau2*pettSum1(1) ) - !-                   alphaP2*( tau1*pettSum2(0) + tau2*pettSum2(1) )
                          !-    end if 
                           ! write(debugFile,'(" u1Lap,ueLap=",2e10.2," v1Lap,veLap=",2e10.2)') u1Lap,ueLap,v1Lap,veLap
                           ! write(debugFile,'(" u1Lap,ueLap=",2e10.2," v1Lap,veLap=",2e10.2)') u1Lap,ueLap,v1Lap,veLap
                         end if
                      do n1=0,numberOfEquations-1
                       f0(n1)=f(n1)
                      end do
                      delta=1.  ! perturb E by this amount 
                      do n2=0,numberOfEquations-1
                        ! pertub one component: 
                         if( nd.eq.2 .or. orderOfAccuracy.eq.2 )then
                          if( n2.lt.nd )then
                            u1(i1-is1,i2-is2,i3-is3,ex+n2   )=u1(i1-is1,i2-is2,i3-is3,ex+n2   )+(delta)
                          else if( n2.lt.2*nd )then
                            u2(j1-js1,j2-js2,j3-js3,ex+n2-nd)= u2(j1-js1,j2-js2,j3-js3,ex+n2-nd)+(delta)
                          else if( n2.lt.3*nd )then
                            u1(i1-2*is1,i2-2*is2,i3-2*is3,ex+n2-2*nd)=u1(i1-2*is1,i2-2*is2,i3-2*is3,ex+n2-2*nd)+(delta)
                          else if( n2.lt.4*nd )then
                            u2(j1-2*js1,j2-2*js2,j3-2*js3,ex+n2-3*nd)=u2(j1-2*js1,j2-2*js2,j3-2*js3,ex+n2-3*nd)+(delta)
                          else
                            ! this should not happen
                            stop 6363
                          end if
                         else
                          ! 3D order 4 has a different ordering from 3D order 2
                          !  (Ex1,Ey1,Ez1)(-1) , (Ex1,Ey1,Ez1)(-2), (Ex2,Ey2,Ez2)(-1) , (Ex2,Ey2,Ez2)(-2), 
                          if( n2.lt.nd )then
                            u1(i1-is1,i2-is2,i3-is3,ex+n2   )=u1(i1-is1,i2-is2,i3-is3,ex+n2   )+(delta)
                          else if( n2.lt.2*nd )then
                            u1(i1-2*is1,i2-2*is2,i3-2*is3,ex+n2-nd)=u1(i1-2*is1,i2-2*is2,i3-2*is3,ex+n2-nd)+(delta)
                          else if( n2.lt.3*nd )then
                            u2(j1-js1,j2-js2,j3-js3,ex+n2-2*nd)= u2(j1-js1,j2-js2,j3-js3,ex+n2-2*nd)+(delta)
                          else if( n2.lt.4*nd )then
                            u2(j1-2*js1,j2-2*js2,j3-2*js3,ex+n2-3*nd)=u2(j1-2*js1,j2-2*js2,j3-2*js3,ex+n2-3*nd)+(delta)
                          else
                            ! this should not happen
                            stop 6363
                          end if
                         end if
                           ! NOTE: the jacobian derivatives can be computed once for all components
                            ! this next call will define the jacobian and its derivatives (parameteric and spatial)
                            aj1rx = rsxy1(i1,i2,i3,0,0)
                            aj1rxr = (-rsxy1(i1-1,i2,i3,0,0)+rsxy1(i1+1,i2,i3,0,0))/(2.*dr1(0))
                            aj1rxs = (-rsxy1(i1,i2-1,i3,0,0)+rsxy1(i1,i2+1,i3,0,0))/(2.*dr1(1))
                            aj1sx = rsxy1(i1,i2,i3,1,0)
                            aj1sxr = (-rsxy1(i1-1,i2,i3,1,0)+rsxy1(i1+1,i2,i3,1,0))/(2.*dr1(0))
                            aj1sxs = (-rsxy1(i1,i2-1,i3,1,0)+rsxy1(i1,i2+1,i3,1,0))/(2.*dr1(1))
                            aj1ry = rsxy1(i1,i2,i3,0,1)
                            aj1ryr = (-rsxy1(i1-1,i2,i3,0,1)+rsxy1(i1+1,i2,i3,0,1))/(2.*dr1(0))
                            aj1rys = (-rsxy1(i1,i2-1,i3,0,1)+rsxy1(i1,i2+1,i3,0,1))/(2.*dr1(1))
                            aj1sy = rsxy1(i1,i2,i3,1,1)
                            aj1syr = (-rsxy1(i1-1,i2,i3,1,1)+rsxy1(i1+1,i2,i3,1,1))/(2.*dr1(0))
                            aj1sys = (-rsxy1(i1,i2-1,i3,1,1)+rsxy1(i1,i2+1,i3,1,1))/(2.*dr1(1))
                            aj1rxx = aj1rx*aj1rxr+aj1sx*aj1rxs
                            aj1rxy = aj1ry*aj1rxr+aj1sy*aj1rxs
                            aj1sxx = aj1rx*aj1sxr+aj1sx*aj1sxs
                            aj1sxy = aj1ry*aj1sxr+aj1sy*aj1sxs
                            aj1ryx = aj1rx*aj1ryr+aj1sx*aj1rys
                            aj1ryy = aj1ry*aj1ryr+aj1sy*aj1rys
                            aj1syx = aj1rx*aj1syr+aj1sx*aj1sys
                            aj1syy = aj1ry*aj1syr+aj1sy*aj1sys
                             uu1 = u1(i1,i2,i3,ex)
                             uu1r = (-u1(i1-1,i2,i3,ex)+u1(i1+1,i2,i3,ex))/(2.*dr1(0))
                             uu1s = (-u1(i1,i2-1,i3,ex)+u1(i1,i2+1,i3,ex))/(2.*dr1(1))
                             uu1rr = (u1(i1-1,i2,i3,ex)-2.*u1(i1,i2,i3,ex)+u1(i1+1,i2,i3,ex))/(dr1(0)**2)
                             uu1rs = (-(-u1(i1-1,i2-1,i3,ex)+u1(i1-1,i2+1,i3,ex))/(2.*dr1(1))+(-u1(i1+1,i2-1,i3,ex)+u1(i1+1,i2+1,i3,ex))/(2.*dr1(1)))/(2.*dr1(0))
                             uu1ss = (u1(i1,i2-1,i3,ex)-2.*u1(i1,i2,i3,ex)+u1(i1,i2+1,i3,ex))/(dr1(1)**2)
                              u1x = aj1rx*uu1r+aj1sx*uu1s
                              u1y = aj1ry*uu1r+aj1sy*uu1s
                              t1 = aj1rx**2
                              t6 = aj1sx**2
                              u1xx = t1*uu1rr+2*aj1rx*aj1sx*uu1rs+t6*uu1ss+aj1rxx*uu1r+aj1sxx*uu1s
                              t1 = aj1ry**2
                              t6 = aj1sy**2
                              u1yy = t1*uu1rr+2*aj1ry*aj1sy*uu1rs+t6*uu1ss+aj1ryy*uu1r+aj1syy*uu1s
                            u1Lap = u1xx+ u1yy
                             vv1 = u1(i1,i2,i3,ey)
                             vv1r = (-u1(i1-1,i2,i3,ey)+u1(i1+1,i2,i3,ey))/(2.*dr1(0))
                             vv1s = (-u1(i1,i2-1,i3,ey)+u1(i1,i2+1,i3,ey))/(2.*dr1(1))
                             vv1rr = (u1(i1-1,i2,i3,ey)-2.*u1(i1,i2,i3,ey)+u1(i1+1,i2,i3,ey))/(dr1(0)**2)
                             vv1rs = (-(-u1(i1-1,i2-1,i3,ey)+u1(i1-1,i2+1,i3,ey))/(2.*dr1(1))+(-u1(i1+1,i2-1,i3,ey)+u1(i1+1,i2+1,i3,ey))/(2.*dr1(1)))/(2.*dr1(0))
                             vv1ss = (u1(i1,i2-1,i3,ey)-2.*u1(i1,i2,i3,ey)+u1(i1,i2+1,i3,ey))/(dr1(1)**2)
                              v1x = aj1rx*vv1r+aj1sx*vv1s
                              v1y = aj1ry*vv1r+aj1sy*vv1s
                              t1 = aj1rx**2
                              t6 = aj1sx**2
                              v1xx = t1*vv1rr+2*aj1rx*aj1sx*vv1rs+t6*vv1ss+aj1rxx*vv1r+aj1sxx*vv1s
                              t1 = aj1ry**2
                              t6 = aj1sy**2
                              v1yy = t1*vv1rr+2*aj1ry*aj1sy*vv1rs+t6*vv1ss+aj1ryy*vv1r+aj1syy*vv1s
                            v1Lap = v1xx+ v1yy
                           ! NOTE: the jacobian derivatives can be computed once for all components
                            ! this next call will define the jacobian and its derivatives (parameteric and spatial)
                            aj2rx = rsxy2(j1,j2,j3,0,0)
                            aj2rxr = (-rsxy2(j1-1,j2,j3,0,0)+rsxy2(j1+1,j2,j3,0,0))/(2.*dr2(0))
                            aj2rxs = (-rsxy2(j1,j2-1,j3,0,0)+rsxy2(j1,j2+1,j3,0,0))/(2.*dr2(1))
                            aj2sx = rsxy2(j1,j2,j3,1,0)
                            aj2sxr = (-rsxy2(j1-1,j2,j3,1,0)+rsxy2(j1+1,j2,j3,1,0))/(2.*dr2(0))
                            aj2sxs = (-rsxy2(j1,j2-1,j3,1,0)+rsxy2(j1,j2+1,j3,1,0))/(2.*dr2(1))
                            aj2ry = rsxy2(j1,j2,j3,0,1)
                            aj2ryr = (-rsxy2(j1-1,j2,j3,0,1)+rsxy2(j1+1,j2,j3,0,1))/(2.*dr2(0))
                            aj2rys = (-rsxy2(j1,j2-1,j3,0,1)+rsxy2(j1,j2+1,j3,0,1))/(2.*dr2(1))
                            aj2sy = rsxy2(j1,j2,j3,1,1)
                            aj2syr = (-rsxy2(j1-1,j2,j3,1,1)+rsxy2(j1+1,j2,j3,1,1))/(2.*dr2(0))
                            aj2sys = (-rsxy2(j1,j2-1,j3,1,1)+rsxy2(j1,j2+1,j3,1,1))/(2.*dr2(1))
                            aj2rxx = aj2rx*aj2rxr+aj2sx*aj2rxs
                            aj2rxy = aj2ry*aj2rxr+aj2sy*aj2rxs
                            aj2sxx = aj2rx*aj2sxr+aj2sx*aj2sxs
                            aj2sxy = aj2ry*aj2sxr+aj2sy*aj2sxs
                            aj2ryx = aj2rx*aj2ryr+aj2sx*aj2rys
                            aj2ryy = aj2ry*aj2ryr+aj2sy*aj2rys
                            aj2syx = aj2rx*aj2syr+aj2sx*aj2sys
                            aj2syy = aj2ry*aj2syr+aj2sy*aj2sys
                             uu2 = u2(j1,j2,j3,ex)
                             uu2r = (-u2(j1-1,j2,j3,ex)+u2(j1+1,j2,j3,ex))/(2.*dr2(0))
                             uu2s = (-u2(j1,j2-1,j3,ex)+u2(j1,j2+1,j3,ex))/(2.*dr2(1))
                             uu2rr = (u2(j1-1,j2,j3,ex)-2.*u2(j1,j2,j3,ex)+u2(j1+1,j2,j3,ex))/(dr2(0)**2)
                             uu2rs = (-(-u2(j1-1,j2-1,j3,ex)+u2(j1-1,j2+1,j3,ex))/(2.*dr2(1))+(-u2(j1+1,j2-1,j3,ex)+u2(j1+1,j2+1,j3,ex))/(2.*dr2(1)))/(2.*dr2(0))
                             uu2ss = (u2(j1,j2-1,j3,ex)-2.*u2(j1,j2,j3,ex)+u2(j1,j2+1,j3,ex))/(dr2(1)**2)
                              u2x = aj2rx*uu2r+aj2sx*uu2s
                              u2y = aj2ry*uu2r+aj2sy*uu2s
                              t1 = aj2rx**2
                              t6 = aj2sx**2
                              u2xx = t1*uu2rr+2*aj2rx*aj2sx*uu2rs+t6*uu2ss+aj2rxx*uu2r+aj2sxx*uu2s
                              t1 = aj2ry**2
                              t6 = aj2sy**2
                              u2yy = t1*uu2rr+2*aj2ry*aj2sy*uu2rs+t6*uu2ss+aj2ryy*uu2r+aj2syy*uu2s
                            u2Lap = u2xx+ u2yy
                             vv2 = u2(j1,j2,j3,ey)
                             vv2r = (-u2(j1-1,j2,j3,ey)+u2(j1+1,j2,j3,ey))/(2.*dr2(0))
                             vv2s = (-u2(j1,j2-1,j3,ey)+u2(j1,j2+1,j3,ey))/(2.*dr2(1))
                             vv2rr = (u2(j1-1,j2,j3,ey)-2.*u2(j1,j2,j3,ey)+u2(j1+1,j2,j3,ey))/(dr2(0)**2)
                             vv2rs = (-(-u2(j1-1,j2-1,j3,ey)+u2(j1-1,j2+1,j3,ey))/(2.*dr2(1))+(-u2(j1+1,j2-1,j3,ey)+u2(j1+1,j2+1,j3,ey))/(2.*dr2(1)))/(2.*dr2(0))
                             vv2ss = (u2(j1,j2-1,j3,ey)-2.*u2(j1,j2,j3,ey)+u2(j1,j2+1,j3,ey))/(dr2(1)**2)
                              v2x = aj2rx*vv2r+aj2sx*vv2s
                              v2y = aj2ry*vv2r+aj2sy*vv2s
                              t1 = aj2rx**2
                              t6 = aj2sx**2
                              v2xx = t1*vv2rr+2*aj2rx*aj2sx*vv2rs+t6*vv2ss+aj2rxx*vv2r+aj2sxx*vv2s
                              t1 = aj2ry**2
                              t6 = aj2sy**2
                              v2yy = t1*vv2rr+2*aj2ry*aj2sy*vv2rs+t6*vv2ss+aj2ryy*vv2r+aj2syy*vv2s
                            v2Lap = v2xx+ v2yy
                           f(0)=(u1x+v1y) - (u2x+v2y)
                           f(1)=( an1*u1Lap +an2*v1Lap )/mu1 - ( an1*u2Lap +an2*v2Lap )/mu2 
                           f(2)=(v1x-u1y)/mu1 - (v2x-u2y)/mu2
                           f(3)=( ( tau1*u1Lap +tau2*v1Lap )*beta1/epsmu1 - alphaP1*(tau1*fp1(0)+tau2*fp1(1)) ) - ( ( tau1*u2Lap +tau2*v2Lap )*beta2/epsmu2 - alphaP2*(tau1*fp2(0)+tau2*fp2(1)) )
                           if( twilightZone.eq.1 )then
                             ! For now we assume mu1=mu2 and TZ solutions are the same on both sides.
                             ! f(3) = [ tv.E.tt] = [ tv.( c^2*Delta(E) - alphaP*P.tt + fev) ] 
                             !      = [ tv.( c^2*Delta(E) - alphaP*P.tt] + [ tv.fev ]
                             ! -- add on the jump in the forcing ---
                             f(3) = f(3) + ( tau1*(fev1(0)-fev2(0)) + tau2*(fev1(1)-fev2(1)) )
                            !-    ! f(3) = [ tv.( c^2*Delta(E) - alphaP*P.tt ] - [ tv.( c^2*Delta(E^e) - alphaP*(P^e).tt ] =0 
                            !-    call ogderiv(ep, 0,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, uexx )
                            !-    call ogderiv(ep, 0,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, ueyy )
                            !-    call ogderiv(ep, 0,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, vexx )
                            !-    call ogderiv(ep, 0,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, veyy )
                            !- 
                            !-    ueLap = uexx + ueyy
                            !-    veLap = vexx + veyy
                            !-    f(3) = f(3) - ( tau1*ueLap +tau2*veLap )*(1./epsmu1-1./epsmu2)
                            !- 
                            !-    f(3) = f(3) +  alphaP1*( tau1*pettSum1(0) + tau2*pettSum1(1) ) - !-                   alphaP2*( tau1*pettSum2(0) + tau2*pettSum2(1) )
                            !-    end if 
                             ! write(debugFile,'(" u1Lap,ueLap=",2e10.2," v1Lap,veLap=",2e10.2)') u1Lap,ueLap,v1Lap,veLap
                             ! write(debugFile,'(" u1Lap,ueLap=",2e10.2," v1Lap,veLap=",2e10.2)') u1Lap,ueLap,v1Lap,veLap
                           end if
                        ! compute the difference
                        do n1=0,numberOfEquations-1
                         f(n1)=f(n1)-f0(n1)
                        end do
                        ! write(*,'(" u1x,v1y,u2x,v2y=",4(1pe10.2))') u1x,v1y, u2x,v2y
                        if( .true. )then
                         if( numberOfEquations.eq.4 )then
                          write(*,'(" a4(*,",i1,")=",4(1pe10.2)," diff(*)=",4(1pe8.1) )') n2,(a4(n1,n2),n1=0,numberOfEquations-1),(a4(n1,n2)-f(n1),n1=0,numberOfEquations-1)
                         else if( numberOfEquations.eq.6 )then
                          write(*,'(" a4(*,",i1,")=",6(1pe10.2)," diff(*)=",6(1pe8.1) )') n2,(a4(n1,n2),n1=0,numberOfEquations-1),(a4(n1,n2)-f(n1),n1=0,numberOfEquations-1)
                         else if( numberOfEquations.eq.8 )then
                          write(*,'(" a4(*,",i1,")=",8(1pe10.2)," diff(*)=",8(1pe8.1) )') n2,(a4(n1,n2),n1=0,numberOfEquations-1),(a4(n1,n2)-f(n1),n1=0,numberOfEquations-1)
                         else if( numberOfEquations.eq.12 )then
                          write(*,'(" a4(*,",i1,")=",12(1pe10.2)," diff(*)=",12(1pe8.1) )') n2,(a4(n1,n2),n1=0,numberOfEquations-1),(a4(n1,n2)-f(n1),n1=0,numberOfEquations-1)
                         else 
                           stop 8181
                         end if
                        else
                         if( numberOfEquations.eq.4 )then
                          write(*,'(" a4(*,",i1,")=",4(1pe10.2)," f(*)=",4(1pe10.2)," diff(*)=",4(1pe8.1) )') n2,(a4(n1,n2),n1=0,numberOfEquations-1),(f(n1),n1=0,numberOfEquations-1),(a4(n1,n2)-f(n1),n1=0,numberOfEquations-1)
                         else if( numberOfEquations.eq.6 )then
                          write(*,'(" a4(*,",i1,")=",6(1pe10.2)," f(*)=",6(1pe10.2)," diff(*)=",6(1pe8.1) )') n2,(a4(n1,n2),n1=0,numberOfEquations-1),(f(n1),n1=0,numberOfEquations-1),(a4(n1,n2)-f(n1),n1=0,numberOfEquations-1)
                         else if( numberOfEquations.eq.8 )then
                          write(*,'(" a4(*,",i1,")=",8(1pe10.2)," f(*)=",8(1pe10.2)," diff(*)=",8(1pe8.1) )') n2,(a4(n1,n2),n1=0,numberOfEquations-1),(f(n1),n1=0,numberOfEquations-1),(a4(n1,n2)-f(n1),n1=0,numberOfEquations-1)
                         else if( numberOfEquations.eq.12 )then
                          write(*,'(" a4(*,",i1,")=",12(1pe10.2)," f(*)=",12(1pe10.2)," diff(*)=",12(1pe8.1) )') n2,(a4(n1,n2),n1=0,numberOfEquations-1),(f(n1),n1=0,numberOfEquations-1),(a4(n1,n2)-f(n1),n1=0,numberOfEquations-1)
                         else 
                           stop 8181
                         end if
                        end if
                        do n1=0,numberOfEquations-1
                          coeffDiff = max(coeffDiff,abs(a4(n1,n2)-f(n1)))
                        end do
                        ! reset pertubation
                         if( nd.eq.2 .or. orderOfAccuracy.eq.2 )then
                          if( n2.lt.nd )then
                            u1(i1-is1,i2-is2,i3-is3,ex+n2   )=u1(i1-is1,i2-is2,i3-is3,ex+n2   )+(-delta)
                          else if( n2.lt.2*nd )then
                            u2(j1-js1,j2-js2,j3-js3,ex+n2-nd)= u2(j1-js1,j2-js2,j3-js3,ex+n2-nd)+(-delta)
                          else if( n2.lt.3*nd )then
                            u1(i1-2*is1,i2-2*is2,i3-2*is3,ex+n2-2*nd)=u1(i1-2*is1,i2-2*is2,i3-2*is3,ex+n2-2*nd)+(-delta)
                          else if( n2.lt.4*nd )then
                            u2(j1-2*js1,j2-2*js2,j3-2*js3,ex+n2-3*nd)=u2(j1-2*js1,j2-2*js2,j3-2*js3,ex+n2-3*nd)+(-delta)
                          else
                            ! this should not happen
                            stop 6363
                          end if
                         else
                          ! 3D order 4 has a different ordering from 3D order 2
                          !  (Ex1,Ey1,Ez1)(-1) , (Ex1,Ey1,Ez1)(-2), (Ex2,Ey2,Ez2)(-1) , (Ex2,Ey2,Ez2)(-2), 
                          if( n2.lt.nd )then
                            u1(i1-is1,i2-is2,i3-is3,ex+n2   )=u1(i1-is1,i2-is2,i3-is3,ex+n2   )+(-delta)
                          else if( n2.lt.2*nd )then
                            u1(i1-2*is1,i2-2*is2,i3-2*is3,ex+n2-nd)=u1(i1-2*is1,i2-2*is2,i3-2*is3,ex+n2-nd)+(-delta)
                          else if( n2.lt.3*nd )then
                            u2(j1-js1,j2-js2,j3-js3,ex+n2-2*nd)= u2(j1-js1,j2-js2,j3-js3,ex+n2-2*nd)+(-delta)
                          else if( n2.lt.4*nd )then
                            u2(j1-2*js1,j2-2*js2,j3-2*js3,ex+n2-3*nd)=u2(j1-2*js1,j2-2*js2,j3-2*js3,ex+n2-3*nd)+(-delta)
                          else
                            ! this should not happen
                            stop 6363
                          end if
                         end if
                      end do 
                      ! restore 
                         ! NOTE: the jacobian derivatives can be computed once for all components
                          ! this next call will define the jacobian and its derivatives (parameteric and spatial)
                          aj1rx = rsxy1(i1,i2,i3,0,0)
                          aj1rxr = (-rsxy1(i1-1,i2,i3,0,0)+rsxy1(i1+1,i2,i3,0,0))/(2.*dr1(0))
                          aj1rxs = (-rsxy1(i1,i2-1,i3,0,0)+rsxy1(i1,i2+1,i3,0,0))/(2.*dr1(1))
                          aj1sx = rsxy1(i1,i2,i3,1,0)
                          aj1sxr = (-rsxy1(i1-1,i2,i3,1,0)+rsxy1(i1+1,i2,i3,1,0))/(2.*dr1(0))
                          aj1sxs = (-rsxy1(i1,i2-1,i3,1,0)+rsxy1(i1,i2+1,i3,1,0))/(2.*dr1(1))
                          aj1ry = rsxy1(i1,i2,i3,0,1)
                          aj1ryr = (-rsxy1(i1-1,i2,i3,0,1)+rsxy1(i1+1,i2,i3,0,1))/(2.*dr1(0))
                          aj1rys = (-rsxy1(i1,i2-1,i3,0,1)+rsxy1(i1,i2+1,i3,0,1))/(2.*dr1(1))
                          aj1sy = rsxy1(i1,i2,i3,1,1)
                          aj1syr = (-rsxy1(i1-1,i2,i3,1,1)+rsxy1(i1+1,i2,i3,1,1))/(2.*dr1(0))
                          aj1sys = (-rsxy1(i1,i2-1,i3,1,1)+rsxy1(i1,i2+1,i3,1,1))/(2.*dr1(1))
                          aj1rxx = aj1rx*aj1rxr+aj1sx*aj1rxs
                          aj1rxy = aj1ry*aj1rxr+aj1sy*aj1rxs
                          aj1sxx = aj1rx*aj1sxr+aj1sx*aj1sxs
                          aj1sxy = aj1ry*aj1sxr+aj1sy*aj1sxs
                          aj1ryx = aj1rx*aj1ryr+aj1sx*aj1rys
                          aj1ryy = aj1ry*aj1ryr+aj1sy*aj1rys
                          aj1syx = aj1rx*aj1syr+aj1sx*aj1sys
                          aj1syy = aj1ry*aj1syr+aj1sy*aj1sys
                           uu1 = u1(i1,i2,i3,ex)
                           uu1r = (-u1(i1-1,i2,i3,ex)+u1(i1+1,i2,i3,ex))/(2.*dr1(0))
                           uu1s = (-u1(i1,i2-1,i3,ex)+u1(i1,i2+1,i3,ex))/(2.*dr1(1))
                           uu1rr = (u1(i1-1,i2,i3,ex)-2.*u1(i1,i2,i3,ex)+u1(i1+1,i2,i3,ex))/(dr1(0)**2)
                           uu1rs = (-(-u1(i1-1,i2-1,i3,ex)+u1(i1-1,i2+1,i3,ex))/(2.*dr1(1))+(-u1(i1+1,i2-1,i3,ex)+u1(i1+1,i2+1,i3,ex))/(2.*dr1(1)))/(2.*dr1(0))
                           uu1ss = (u1(i1,i2-1,i3,ex)-2.*u1(i1,i2,i3,ex)+u1(i1,i2+1,i3,ex))/(dr1(1)**2)
                            u1x = aj1rx*uu1r+aj1sx*uu1s
                            u1y = aj1ry*uu1r+aj1sy*uu1s
                            t1 = aj1rx**2
                            t6 = aj1sx**2
                            u1xx = t1*uu1rr+2*aj1rx*aj1sx*uu1rs+t6*uu1ss+aj1rxx*uu1r+aj1sxx*uu1s
                            t1 = aj1ry**2
                            t6 = aj1sy**2
                            u1yy = t1*uu1rr+2*aj1ry*aj1sy*uu1rs+t6*uu1ss+aj1ryy*uu1r+aj1syy*uu1s
                          u1Lap = u1xx+ u1yy
                           vv1 = u1(i1,i2,i3,ey)
                           vv1r = (-u1(i1-1,i2,i3,ey)+u1(i1+1,i2,i3,ey))/(2.*dr1(0))
                           vv1s = (-u1(i1,i2-1,i3,ey)+u1(i1,i2+1,i3,ey))/(2.*dr1(1))
                           vv1rr = (u1(i1-1,i2,i3,ey)-2.*u1(i1,i2,i3,ey)+u1(i1+1,i2,i3,ey))/(dr1(0)**2)
                           vv1rs = (-(-u1(i1-1,i2-1,i3,ey)+u1(i1-1,i2+1,i3,ey))/(2.*dr1(1))+(-u1(i1+1,i2-1,i3,ey)+u1(i1+1,i2+1,i3,ey))/(2.*dr1(1)))/(2.*dr1(0))
                           vv1ss = (u1(i1,i2-1,i3,ey)-2.*u1(i1,i2,i3,ey)+u1(i1,i2+1,i3,ey))/(dr1(1)**2)
                            v1x = aj1rx*vv1r+aj1sx*vv1s
                            v1y = aj1ry*vv1r+aj1sy*vv1s
                            t1 = aj1rx**2
                            t6 = aj1sx**2
                            v1xx = t1*vv1rr+2*aj1rx*aj1sx*vv1rs+t6*vv1ss+aj1rxx*vv1r+aj1sxx*vv1s
                            t1 = aj1ry**2
                            t6 = aj1sy**2
                            v1yy = t1*vv1rr+2*aj1ry*aj1sy*vv1rs+t6*vv1ss+aj1ryy*vv1r+aj1syy*vv1s
                          v1Lap = v1xx+ v1yy
                         ! NOTE: the jacobian derivatives can be computed once for all components
                          ! this next call will define the jacobian and its derivatives (parameteric and spatial)
                          aj2rx = rsxy2(j1,j2,j3,0,0)
                          aj2rxr = (-rsxy2(j1-1,j2,j3,0,0)+rsxy2(j1+1,j2,j3,0,0))/(2.*dr2(0))
                          aj2rxs = (-rsxy2(j1,j2-1,j3,0,0)+rsxy2(j1,j2+1,j3,0,0))/(2.*dr2(1))
                          aj2sx = rsxy2(j1,j2,j3,1,0)
                          aj2sxr = (-rsxy2(j1-1,j2,j3,1,0)+rsxy2(j1+1,j2,j3,1,0))/(2.*dr2(0))
                          aj2sxs = (-rsxy2(j1,j2-1,j3,1,0)+rsxy2(j1,j2+1,j3,1,0))/(2.*dr2(1))
                          aj2ry = rsxy2(j1,j2,j3,0,1)
                          aj2ryr = (-rsxy2(j1-1,j2,j3,0,1)+rsxy2(j1+1,j2,j3,0,1))/(2.*dr2(0))
                          aj2rys = (-rsxy2(j1,j2-1,j3,0,1)+rsxy2(j1,j2+1,j3,0,1))/(2.*dr2(1))
                          aj2sy = rsxy2(j1,j2,j3,1,1)
                          aj2syr = (-rsxy2(j1-1,j2,j3,1,1)+rsxy2(j1+1,j2,j3,1,1))/(2.*dr2(0))
                          aj2sys = (-rsxy2(j1,j2-1,j3,1,1)+rsxy2(j1,j2+1,j3,1,1))/(2.*dr2(1))
                          aj2rxx = aj2rx*aj2rxr+aj2sx*aj2rxs
                          aj2rxy = aj2ry*aj2rxr+aj2sy*aj2rxs
                          aj2sxx = aj2rx*aj2sxr+aj2sx*aj2sxs
                          aj2sxy = aj2ry*aj2sxr+aj2sy*aj2sxs
                          aj2ryx = aj2rx*aj2ryr+aj2sx*aj2rys
                          aj2ryy = aj2ry*aj2ryr+aj2sy*aj2rys
                          aj2syx = aj2rx*aj2syr+aj2sx*aj2sys
                          aj2syy = aj2ry*aj2syr+aj2sy*aj2sys
                           uu2 = u2(j1,j2,j3,ex)
                           uu2r = (-u2(j1-1,j2,j3,ex)+u2(j1+1,j2,j3,ex))/(2.*dr2(0))
                           uu2s = (-u2(j1,j2-1,j3,ex)+u2(j1,j2+1,j3,ex))/(2.*dr2(1))
                           uu2rr = (u2(j1-1,j2,j3,ex)-2.*u2(j1,j2,j3,ex)+u2(j1+1,j2,j3,ex))/(dr2(0)**2)
                           uu2rs = (-(-u2(j1-1,j2-1,j3,ex)+u2(j1-1,j2+1,j3,ex))/(2.*dr2(1))+(-u2(j1+1,j2-1,j3,ex)+u2(j1+1,j2+1,j3,ex))/(2.*dr2(1)))/(2.*dr2(0))
                           uu2ss = (u2(j1,j2-1,j3,ex)-2.*u2(j1,j2,j3,ex)+u2(j1,j2+1,j3,ex))/(dr2(1)**2)
                            u2x = aj2rx*uu2r+aj2sx*uu2s
                            u2y = aj2ry*uu2r+aj2sy*uu2s
                            t1 = aj2rx**2
                            t6 = aj2sx**2
                            u2xx = t1*uu2rr+2*aj2rx*aj2sx*uu2rs+t6*uu2ss+aj2rxx*uu2r+aj2sxx*uu2s
                            t1 = aj2ry**2
                            t6 = aj2sy**2
                            u2yy = t1*uu2rr+2*aj2ry*aj2sy*uu2rs+t6*uu2ss+aj2ryy*uu2r+aj2syy*uu2s
                          u2Lap = u2xx+ u2yy
                           vv2 = u2(j1,j2,j3,ey)
                           vv2r = (-u2(j1-1,j2,j3,ey)+u2(j1+1,j2,j3,ey))/(2.*dr2(0))
                           vv2s = (-u2(j1,j2-1,j3,ey)+u2(j1,j2+1,j3,ey))/(2.*dr2(1))
                           vv2rr = (u2(j1-1,j2,j3,ey)-2.*u2(j1,j2,j3,ey)+u2(j1+1,j2,j3,ey))/(dr2(0)**2)
                           vv2rs = (-(-u2(j1-1,j2-1,j3,ey)+u2(j1-1,j2+1,j3,ey))/(2.*dr2(1))+(-u2(j1+1,j2-1,j3,ey)+u2(j1+1,j2+1,j3,ey))/(2.*dr2(1)))/(2.*dr2(0))
                           vv2ss = (u2(j1,j2-1,j3,ey)-2.*u2(j1,j2,j3,ey)+u2(j1,j2+1,j3,ey))/(dr2(1)**2)
                            v2x = aj2rx*vv2r+aj2sx*vv2s
                            v2y = aj2ry*vv2r+aj2sy*vv2s
                            t1 = aj2rx**2
                            t6 = aj2sx**2
                            v2xx = t1*vv2rr+2*aj2rx*aj2sx*vv2rs+t6*vv2ss+aj2rxx*vv2r+aj2sxx*vv2s
                            t1 = aj2ry**2
                            t6 = aj2sy**2
                            v2yy = t1*vv2rr+2*aj2ry*aj2sy*vv2rs+t6*vv2ss+aj2ryy*vv2r+aj2syy*vv2s
                          v2Lap = v2xx+ v2yy
                         f(0)=(u1x+v1y) - (u2x+v2y)
                         f(1)=( an1*u1Lap +an2*v1Lap )/mu1 - ( an1*u2Lap +an2*v2Lap )/mu2 
                         f(2)=(v1x-u1y)/mu1 - (v2x-u2y)/mu2
                         f(3)=( ( tau1*u1Lap +tau2*v1Lap )*beta1/epsmu1 - alphaP1*(tau1*fp1(0)+tau2*fp1(1)) ) - ( ( tau1*u2Lap +tau2*v2Lap )*beta2/epsmu2 - alphaP2*(tau1*fp2(0)+tau2*fp2(1)) )
                         if( twilightZone.eq.1 )then
                           ! For now we assume mu1=mu2 and TZ solutions are the same on both sides.
                           ! f(3) = [ tv.E.tt] = [ tv.( c^2*Delta(E) - alphaP*P.tt + fev) ] 
                           !      = [ tv.( c^2*Delta(E) - alphaP*P.tt] + [ tv.fev ]
                           ! -- add on the jump in the forcing ---
                           f(3) = f(3) + ( tau1*(fev1(0)-fev2(0)) + tau2*(fev1(1)-fev2(1)) )
                          !-    ! f(3) = [ tv.( c^2*Delta(E) - alphaP*P.tt ] - [ tv.( c^2*Delta(E^e) - alphaP*(P^e).tt ] =0 
                          !-    call ogderiv(ep, 0,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, uexx )
                          !-    call ogderiv(ep, 0,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, ueyy )
                          !-    call ogderiv(ep, 0,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, vexx )
                          !-    call ogderiv(ep, 0,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, veyy )
                          !- 
                          !-    ueLap = uexx + ueyy
                          !-    veLap = vexx + veyy
                          !-    f(3) = f(3) - ( tau1*ueLap +tau2*veLap )*(1./epsmu1-1./epsmu2)
                          !- 
                          !-    f(3) = f(3) +  alphaP1*( tau1*pettSum1(0) + tau2*pettSum1(1) ) - !-                   alphaP2*( tau1*pettSum2(0) + tau2*pettSum2(1) )
                          !-    end if 
                           ! write(debugFile,'(" u1Lap,ueLap=",2e10.2," v1Lap,veLap=",2e10.2)') u1Lap,ueLap,v1Lap,veLap
                           ! write(debugFile,'(" u1Lap,ueLap=",2e10.2," v1Lap,veLap=",2e10.2)') u1Lap,ueLap,v1Lap,veLap
                         end if
                    print*,'Checked coefficients using delta function approach'
                  end if
                  ! write(debugFile,'(" --> xy1=",4f8.3)') xy1(i1,i2,i3,0),xy1(i1,i2,i3,1)
                  ! write(debugFile,'(" --> rsxy1=",4f8.3)') rsxy1(i1,i2,i3,0,0),rsxy1(i1,i2,i3,1,0),rsxy1(i1,i2,i3,0,1),rsxy1(i1,i2,i3,1,1)
                  ! write(debugFile,'(" --> rsxy2=",4f8.3)') rsxy2(j1,j2,j3,0,0),rsxy2(j1,j2,j3,1,0),rsxy2(j1,j2,j3,0,1),rsxy2(j1,j2,j3,1,1)
                  ! write(debugFile,'(" --> rxx1=",2f8.3)') rxx1(axis1,0,0),rxx1(axis1,1,1)
                  ! write(debugFile,'(" --> rxx2=",2f8.3)') rxx2(axis2,0,0),rxx2(axis1,1,1)
                  ! write(debugFile,'(" --> a4(0,.)=",4f8.3)') a4(0,0),a4(0,1),a4(0,2),a4(0,3)
                  ! write(debugFile,'(" --> a4(1,.)=",4f8.3)') a4(1,0),a4(1,1),a4(1,2),a4(1,3)
                  ! write(debugFile,'(" --> a4(2,.)=",4f8.3)') a4(2,0),a4(2,1),a4(2,2),a4(2,3)
                  ! write(debugFile,'(" --> a4(3,.)=",4f8.3)') a4(3,0),a4(3,1),a4(3,2),a4(3,3)
                  ! write(debugFile,'(" --> an1,an2=",2f8.3)') an1,an2
                  ! write(debugFile,'(" --> clap1,clap2=",2f8.3)') clap1,clap2
                  ! subtract off the contributions from the wrong values at the ghost points:
                  do n=0,3
                    f(n) = (a4(n,0)*q(0)+a4(n,1)*q(1)+a4(n,2)*q(2)+a4(n,3)*q(3)) - f(n)
                  end do
                  ! write(debugFile,'(" --> order2-curv: i1,i2=",2i4," f(subtract)=",4f8.3)') i1,i2,f(0),f(1),f(2),f(3)
                  ! solve A Q = F
                  ! factor the matrix
                  numberOfEquations=4
                  call dgeco( a4(0,0), numberOfEquations, numberOfEquations, ipvt(0),rcond,work(0))
                  ! solve
                  ! write(debugFile,'(" --> order2-curv: i1,i2=",2i4," rcond=",e10.2)') i1,i2,rcond
                  job=0
                  call dgesl( a4(0,0), numberOfEquations, numberOfEquations, ipvt(0), f(0), job)
                  ! write(debugFile,'(" --> order2-curv: i1,i2=",2i4," f(solve)=",4f8.3)') i1,i2,f(0),f(1),f(2),f(3)
                  u1(i1-is1,i2-is2,i3,ex)=f(0)
                  u1(i1-is1,i2-is2,i3,ey)=f(1)
                  u2(j1-js1,j2-js2,j3,ex)=f(2)
                  u2(j1-js1,j2-js2,j3,ey)=f(3)
                  if( .false. .or. debug.gt.3 )then ! re-evaluate
                     ! NOTE: the jacobian derivatives can be computed once for all components
                      ! this next call will define the jacobian and its derivatives (parameteric and spatial)
                      aj1rx = rsxy1(i1,i2,i3,0,0)
                      aj1rxr = (-rsxy1(i1-1,i2,i3,0,0)+rsxy1(i1+1,i2,i3,0,0))/(2.*dr1(0))
                      aj1rxs = (-rsxy1(i1,i2-1,i3,0,0)+rsxy1(i1,i2+1,i3,0,0))/(2.*dr1(1))
                      aj1sx = rsxy1(i1,i2,i3,1,0)
                      aj1sxr = (-rsxy1(i1-1,i2,i3,1,0)+rsxy1(i1+1,i2,i3,1,0))/(2.*dr1(0))
                      aj1sxs = (-rsxy1(i1,i2-1,i3,1,0)+rsxy1(i1,i2+1,i3,1,0))/(2.*dr1(1))
                      aj1ry = rsxy1(i1,i2,i3,0,1)
                      aj1ryr = (-rsxy1(i1-1,i2,i3,0,1)+rsxy1(i1+1,i2,i3,0,1))/(2.*dr1(0))
                      aj1rys = (-rsxy1(i1,i2-1,i3,0,1)+rsxy1(i1,i2+1,i3,0,1))/(2.*dr1(1))
                      aj1sy = rsxy1(i1,i2,i3,1,1)
                      aj1syr = (-rsxy1(i1-1,i2,i3,1,1)+rsxy1(i1+1,i2,i3,1,1))/(2.*dr1(0))
                      aj1sys = (-rsxy1(i1,i2-1,i3,1,1)+rsxy1(i1,i2+1,i3,1,1))/(2.*dr1(1))
                      aj1rxx = aj1rx*aj1rxr+aj1sx*aj1rxs
                      aj1rxy = aj1ry*aj1rxr+aj1sy*aj1rxs
                      aj1sxx = aj1rx*aj1sxr+aj1sx*aj1sxs
                      aj1sxy = aj1ry*aj1sxr+aj1sy*aj1sxs
                      aj1ryx = aj1rx*aj1ryr+aj1sx*aj1rys
                      aj1ryy = aj1ry*aj1ryr+aj1sy*aj1rys
                      aj1syx = aj1rx*aj1syr+aj1sx*aj1sys
                      aj1syy = aj1ry*aj1syr+aj1sy*aj1sys
                       uu1 = u1(i1,i2,i3,ex)
                       uu1r = (-u1(i1-1,i2,i3,ex)+u1(i1+1,i2,i3,ex))/(2.*dr1(0))
                       uu1s = (-u1(i1,i2-1,i3,ex)+u1(i1,i2+1,i3,ex))/(2.*dr1(1))
                       uu1rr = (u1(i1-1,i2,i3,ex)-2.*u1(i1,i2,i3,ex)+u1(i1+1,i2,i3,ex))/(dr1(0)**2)
                       uu1rs = (-(-u1(i1-1,i2-1,i3,ex)+u1(i1-1,i2+1,i3,ex))/(2.*dr1(1))+(-u1(i1+1,i2-1,i3,ex)+u1(i1+1,i2+1,i3,ex))/(2.*dr1(1)))/(2.*dr1(0))
                       uu1ss = (u1(i1,i2-1,i3,ex)-2.*u1(i1,i2,i3,ex)+u1(i1,i2+1,i3,ex))/(dr1(1)**2)
                        u1x = aj1rx*uu1r+aj1sx*uu1s
                        u1y = aj1ry*uu1r+aj1sy*uu1s
                        t1 = aj1rx**2
                        t6 = aj1sx**2
                        u1xx = t1*uu1rr+2*aj1rx*aj1sx*uu1rs+t6*uu1ss+aj1rxx*uu1r+aj1sxx*uu1s
                        t1 = aj1ry**2
                        t6 = aj1sy**2
                        u1yy = t1*uu1rr+2*aj1ry*aj1sy*uu1rs+t6*uu1ss+aj1ryy*uu1r+aj1syy*uu1s
                      u1Lap = u1xx+ u1yy
                       vv1 = u1(i1,i2,i3,ey)
                       vv1r = (-u1(i1-1,i2,i3,ey)+u1(i1+1,i2,i3,ey))/(2.*dr1(0))
                       vv1s = (-u1(i1,i2-1,i3,ey)+u1(i1,i2+1,i3,ey))/(2.*dr1(1))
                       vv1rr = (u1(i1-1,i2,i3,ey)-2.*u1(i1,i2,i3,ey)+u1(i1+1,i2,i3,ey))/(dr1(0)**2)
                       vv1rs = (-(-u1(i1-1,i2-1,i3,ey)+u1(i1-1,i2+1,i3,ey))/(2.*dr1(1))+(-u1(i1+1,i2-1,i3,ey)+u1(i1+1,i2+1,i3,ey))/(2.*dr1(1)))/(2.*dr1(0))
                       vv1ss = (u1(i1,i2-1,i3,ey)-2.*u1(i1,i2,i3,ey)+u1(i1,i2+1,i3,ey))/(dr1(1)**2)
                        v1x = aj1rx*vv1r+aj1sx*vv1s
                        v1y = aj1ry*vv1r+aj1sy*vv1s
                        t1 = aj1rx**2
                        t6 = aj1sx**2
                        v1xx = t1*vv1rr+2*aj1rx*aj1sx*vv1rs+t6*vv1ss+aj1rxx*vv1r+aj1sxx*vv1s
                        t1 = aj1ry**2
                        t6 = aj1sy**2
                        v1yy = t1*vv1rr+2*aj1ry*aj1sy*vv1rs+t6*vv1ss+aj1ryy*vv1r+aj1syy*vv1s
                      v1Lap = v1xx+ v1yy
                     ! NOTE: the jacobian derivatives can be computed once for all components
                      ! this next call will define the jacobian and its derivatives (parameteric and spatial)
                      aj2rx = rsxy2(j1,j2,j3,0,0)
                      aj2rxr = (-rsxy2(j1-1,j2,j3,0,0)+rsxy2(j1+1,j2,j3,0,0))/(2.*dr2(0))
                      aj2rxs = (-rsxy2(j1,j2-1,j3,0,0)+rsxy2(j1,j2+1,j3,0,0))/(2.*dr2(1))
                      aj2sx = rsxy2(j1,j2,j3,1,0)
                      aj2sxr = (-rsxy2(j1-1,j2,j3,1,0)+rsxy2(j1+1,j2,j3,1,0))/(2.*dr2(0))
                      aj2sxs = (-rsxy2(j1,j2-1,j3,1,0)+rsxy2(j1,j2+1,j3,1,0))/(2.*dr2(1))
                      aj2ry = rsxy2(j1,j2,j3,0,1)
                      aj2ryr = (-rsxy2(j1-1,j2,j3,0,1)+rsxy2(j1+1,j2,j3,0,1))/(2.*dr2(0))
                      aj2rys = (-rsxy2(j1,j2-1,j3,0,1)+rsxy2(j1,j2+1,j3,0,1))/(2.*dr2(1))
                      aj2sy = rsxy2(j1,j2,j3,1,1)
                      aj2syr = (-rsxy2(j1-1,j2,j3,1,1)+rsxy2(j1+1,j2,j3,1,1))/(2.*dr2(0))
                      aj2sys = (-rsxy2(j1,j2-1,j3,1,1)+rsxy2(j1,j2+1,j3,1,1))/(2.*dr2(1))
                      aj2rxx = aj2rx*aj2rxr+aj2sx*aj2rxs
                      aj2rxy = aj2ry*aj2rxr+aj2sy*aj2rxs
                      aj2sxx = aj2rx*aj2sxr+aj2sx*aj2sxs
                      aj2sxy = aj2ry*aj2sxr+aj2sy*aj2sxs
                      aj2ryx = aj2rx*aj2ryr+aj2sx*aj2rys
                      aj2ryy = aj2ry*aj2ryr+aj2sy*aj2rys
                      aj2syx = aj2rx*aj2syr+aj2sx*aj2sys
                      aj2syy = aj2ry*aj2syr+aj2sy*aj2sys
                       uu2 = u2(j1,j2,j3,ex)
                       uu2r = (-u2(j1-1,j2,j3,ex)+u2(j1+1,j2,j3,ex))/(2.*dr2(0))
                       uu2s = (-u2(j1,j2-1,j3,ex)+u2(j1,j2+1,j3,ex))/(2.*dr2(1))
                       uu2rr = (u2(j1-1,j2,j3,ex)-2.*u2(j1,j2,j3,ex)+u2(j1+1,j2,j3,ex))/(dr2(0)**2)
                       uu2rs = (-(-u2(j1-1,j2-1,j3,ex)+u2(j1-1,j2+1,j3,ex))/(2.*dr2(1))+(-u2(j1+1,j2-1,j3,ex)+u2(j1+1,j2+1,j3,ex))/(2.*dr2(1)))/(2.*dr2(0))
                       uu2ss = (u2(j1,j2-1,j3,ex)-2.*u2(j1,j2,j3,ex)+u2(j1,j2+1,j3,ex))/(dr2(1)**2)
                        u2x = aj2rx*uu2r+aj2sx*uu2s
                        u2y = aj2ry*uu2r+aj2sy*uu2s
                        t1 = aj2rx**2
                        t6 = aj2sx**2
                        u2xx = t1*uu2rr+2*aj2rx*aj2sx*uu2rs+t6*uu2ss+aj2rxx*uu2r+aj2sxx*uu2s
                        t1 = aj2ry**2
                        t6 = aj2sy**2
                        u2yy = t1*uu2rr+2*aj2ry*aj2sy*uu2rs+t6*uu2ss+aj2ryy*uu2r+aj2syy*uu2s
                      u2Lap = u2xx+ u2yy
                       vv2 = u2(j1,j2,j3,ey)
                       vv2r = (-u2(j1-1,j2,j3,ey)+u2(j1+1,j2,j3,ey))/(2.*dr2(0))
                       vv2s = (-u2(j1,j2-1,j3,ey)+u2(j1,j2+1,j3,ey))/(2.*dr2(1))
                       vv2rr = (u2(j1-1,j2,j3,ey)-2.*u2(j1,j2,j3,ey)+u2(j1+1,j2,j3,ey))/(dr2(0)**2)
                       vv2rs = (-(-u2(j1-1,j2-1,j3,ey)+u2(j1-1,j2+1,j3,ey))/(2.*dr2(1))+(-u2(j1+1,j2-1,j3,ey)+u2(j1+1,j2+1,j3,ey))/(2.*dr2(1)))/(2.*dr2(0))
                       vv2ss = (u2(j1,j2-1,j3,ey)-2.*u2(j1,j2,j3,ey)+u2(j1,j2+1,j3,ey))/(dr2(1)**2)
                        v2x = aj2rx*vv2r+aj2sx*vv2s
                        v2y = aj2ry*vv2r+aj2sy*vv2s
                        t1 = aj2rx**2
                        t6 = aj2sx**2
                        v2xx = t1*vv2rr+2*aj2rx*aj2sx*vv2rs+t6*vv2ss+aj2rxx*vv2r+aj2sxx*vv2s
                        t1 = aj2ry**2
                        t6 = aj2sy**2
                        v2yy = t1*vv2rr+2*aj2ry*aj2sy*vv2rs+t6*vv2ss+aj2ryy*vv2r+aj2syy*vv2s
                      v2Lap = v2xx+ v2yy
                     f(0)=(u1x+v1y) - (u2x+v2y)
                     f(1)=( an1*u1Lap +an2*v1Lap )/mu1 - ( an1*u2Lap +an2*v2Lap )/mu2 
                     f(2)=(v1x-u1y)/mu1 - (v2x-u2y)/mu2
                     f(3)=( ( tau1*u1Lap +tau2*v1Lap )*beta1/epsmu1 - alphaP1*(tau1*fp1(0)+tau2*fp1(1)) ) - ( ( tau1*u2Lap +tau2*v2Lap )*beta2/epsmu2 - alphaP2*(tau1*fp2(0)+tau2*fp2(1)) )
                     if( twilightZone.eq.1 )then
                       ! For now we assume mu1=mu2 and TZ solutions are the same on both sides.
                       ! f(3) = [ tv.E.tt] = [ tv.( c^2*Delta(E) - alphaP*P.tt + fev) ] 
                       !      = [ tv.( c^2*Delta(E) - alphaP*P.tt] + [ tv.fev ]
                       ! -- add on the jump in the forcing ---
                       f(3) = f(3) + ( tau1*(fev1(0)-fev2(0)) + tau2*(fev1(1)-fev2(1)) )
                      !-    ! f(3) = [ tv.( c^2*Delta(E) - alphaP*P.tt ] - [ tv.( c^2*Delta(E^e) - alphaP*(P^e).tt ] =0 
                      !-    call ogderiv(ep, 0,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, uexx )
                      !-    call ogderiv(ep, 0,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, ueyy )
                      !-    call ogderiv(ep, 0,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, vexx )
                      !-    call ogderiv(ep, 0,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, veyy )
                      !- 
                      !-    ueLap = uexx + ueyy
                      !-    veLap = vexx + veyy
                      !-    f(3) = f(3) - ( tau1*ueLap +tau2*veLap )*(1./epsmu1-1./epsmu2)
                      !- 
                      !-    f(3) = f(3) +  alphaP1*( tau1*pettSum1(0) + tau2*pettSum1(1) ) - !-                   alphaP2*( tau1*pettSum2(0) + tau2*pettSum2(1) )
                      !-    end if 
                       ! write(debugFile,'(" u1Lap,ueLap=",2e10.2," v1Lap,veLap=",2e10.2)') u1Lap,ueLap,v1Lap,veLap
                       ! write(debugFile,'(" u1Lap,ueLap=",2e10.2," v1Lap,veLap=",2e10.2)') u1Lap,ueLap,v1Lap,veLap
                     end if
                    !write(debugFile,'(" --> order2-curv: xy1(ghost)=",2e11.3)') xy1(i1-is1,i2-is2,i3,0),xy1(i1-is1,i2-is2,i3,1)
                    !write(debugFile,'(" --> order2-curv: xy2(ghost)=",2e11.3)') xy2(j1-js1,j2-js2,j3,0),xy2(j1-js1,j2-js2,j3,1)
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
                  ! -- Hz has already been filled in by extrapolation ----
                   end if
                   j1=j1+1
                  end do
                  j2=j2+1
                 end do
                 if( checkCoeff.eq.1 )then
                   write(*,'("+++++ nonlinearMLA22c: check coeff in interface: max(diff) = ",1pe8.2)') coeffDiff
                 end if
            end if
            ! now make sure that div(u)=0 etc.
            !2          if( .false. )then
            !2         beginLoops2d() ! =============== start loops =======================
            !2
            !2           ! 0  [ u.x + v.y ] = 0
            !2           ! first evaluate the equations we want to solve with the wrong values at the ghost points:
            !2           divu=u1x22(i1,i2,i3,ex)+u1y22(i1,i2,i3,ey)
            !2           a0=-is*rsxy1(i1,i2,i3,axis1,0)*dr112(axis1)
            !2           a1=-is*rsxy1(i1,i2,i3,axis1,1)*dr112(axis1)
            !2           aNormSq=a0**2+a1**2
            !2           ! now project:  a.uNew = a.uOld - div  ->  (div-a.uOld)+a.uNew = div(uNew) = 0
            !2           u1(i1-is1,i2-is2,i3,ex)=u1(i1-is1,i2-is2,i3,ex)-divu*a0/aNormSq
            !2           u1(i1-is1,i2-is2,i3,ey)=u1(i1-is1,i2-is2,i3,ey)-divu*a1/aNormSq
            !2
            !2           divu=u2x22(j1,j2,j3,ex)+u2y22(j1,j2,j3,ey)
            !2           a0=-js*rsxy2(j1,j2,j3,axis2,0)*dr212(axis2) 
            !2           a1=-js*rsxy2(j1,j2,j3,axis2,1)*dr212(axis2) 
            !2           aNormSq=a0**2+a1**2
            !2
            !2           u2(j1-js1,j2-js2,j3,ex)=u2(j1-js1,j2-js2,j3,ex)-divu*a0/aNormSq
            !2           u2(j1-js1,j2-js2,j3,ey)=u2(j1-js1,j2-js2,j3,ey)-divu*a1/aNormSq
            !2
            !2           if( debug.gt.0 )then
            !2             write(debugFile,'(" --> 2cth: eval div1,div2=",2e10.2)') u1x22(i1,i2,i3,ex)+u1y22(i1,i2,i3,ey),u2x22(j1,j2,j3,ex)+u2y22(j1,j2,j3,ey)
            !2           end if
            !2         endLoops2d()
            !2          end if
            ! periodic update **** THIS WON T WORK IN PARALLEL
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
           end if ! end assignInterfaceGhostValues
          else if( .false. .and. orderOfAccuracy.eq.4 )then
            ! for testing -- just assign from the other ghost points
             i3=n3a
             j3=m3a
             j2=m2a
             do i2=n2a,n2b
              j1=m1a
              do i1=n1a,n1b
              u1(i1-is1,i2-is2,i3,ex)=u2(j1+js1,j2+js2,j3,ex)
              u1(i1-is1,i2-is2,i3,ey)=u2(j1+js1,j2+js2,j3,ey)
              u1(i1-is1,i2-is2,i3,hz)=u2(j1+js1,j2+js2,j3,hz) 
              u2(j1-js1,j2-js2,j3,ex)=u1(i1+is1,i2+is2,i3,ex)
              u2(j1-js1,j2-js2,j3,ey)=u1(i1+is1,i2+is2,i3,ey)
              u2(j1-js1,j2-js2,j3,hz)=u1(i1+is1,i2+is2,i3,hz)
              u1(i1-2*is1,i2-2*is2,i3,ex)=u2(j1+2*js1,j2+2*js2,j3,ex)
              u1(i1-2*is1,i2-2*is2,i3,ey)=u2(j1+2*js1,j2+2*js2,j3,ey)
              u1(i1-2*is1,i2-2*is2,i3,hz)=u2(j1+2*js1,j2+2*js2,j3,hz) 
              u2(j1-2*js1,j2-2*js2,j3,ex)=u1(i1+2*is1,i2+2*is2,i3,ex)
              u2(j1-2*js1,j2-2*js2,j3,ey)=u1(i1+2*is1,i2+2*is2,i3,ey)
              u2(j1-2*js1,j2-2*js2,j3,hz)=u1(i1+2*is1,i2+2*is2,i3,hz)
               j1=j1+1
              end do
              j2=j2+1
             end do
         ! End curvilinear, 2D, order 2
         ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        else if( nd.eq.2 .and. orderOfAccuracy.eq.4 .and. gridType.eq.rectangular )then
         ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
