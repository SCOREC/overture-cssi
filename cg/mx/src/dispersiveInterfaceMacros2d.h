!         -*- mode: F90 -*-
! *********************************************************************
! ********* MACROS FOR 2D DISPERSIVE INTERFACE CONDITIONS *************
!    This file is included into interfaceOpt.bf90
! *********************************************************************

! --------------------------------------------------------------------------
! Macro: Assign interface ghost values, DIM=2, ORDER=4, GRID=Rectangular 
!         DISPERSIVE CASE -- GDM 
!
! *new* Feb 01, 2021
! --------------------------------------------------------------------------
#beginMacro assignDispersiveInterfaceGhost24r()

 ! ****************************************************************
 ! ***********  DISPERSIVE, 2D, ORDER=4, RECTANGULAR **************
 ! ****************************************************************

 INFO("24r-GDM")

 if( t.le.5*dt .or. debug.gt.3 )then
   if( it.le.2 )then
     write(*,'("macro: assignDispersiveInterfaceGhost24r : it=",i6," t,dt=",2e10.2)') it,t,dt
   end if
 end if

 ! normal and tangent (for TZ forcing)
 an1=an1Cartesian
 an2=an2Cartesian
 tau1=-an2
 tau2= an1

! make sure the normal and tangent are set
if( abs( an1**2 + an2**2 -1. )>1.e-10 .or. abs( tau1**2 + tau2**2 -1. )>1.e-10 )then
  write(*,'("gdm24r - ERROR: incorrect an1,an2, tau1,tau2=",4(1pe9.2))') an1,an2,tau1,tau2
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

  do jv=0,numberOfPolarizationVectors1-1
    fpv1(n,jv)=0.
    LfP1(n,jv)=0.
    fPt1(n,jv)=0.
    fPtt1(n,jv)=0.

    fpvx1(n,jv)=0.
    fpvy1(n,jv)=0.
  end do
  do jv=0,numberOfPolarizationVectors2-1
    fpv2(n,jv)=0.
    LfP2(n,jv)=0.
    fPt2(n,jv)=0.
    fPtt2(n,jv)=0.

    fpvx2(n,jv)=0.
    fpvy2(n,jv)=0.
  end do
end do

 ! write(*,'("p1=",(15(e10.2,1x)))') (((p1(i1,i2,i3,0),i1=nd1a,nd1b),i2=nd2a,nd2b),i3=nd3a,nd3b)


 ! =============== start loops ======================
 beginLoopsMask2d() 

   ! Evaluate the jump conditions using the wrong values at the ghost points 
   evaluateDispersiveInterfaceEquations2dOrder4()

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
   a8(4,0)= ( is*rx1 *2.*dxxx1by2i )*c1**2/mu1 
   a8(4,1)= 0.
   a8(4,2)=-( js*rx2* 2.*dxxx2by2i )*c2**2/mu2 
   a8(4,3)= 0.
   a8(4,4)= (-is*rx1    *dxxx1by2i )*c1**2/mu1 
   a8(4,5)= 0.
   a8(4,6)=-(-js*rx2    *dxxx2by2i )*c2**2/mu2    
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
     write(*,*) "gdm42r: Matrix a8"
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
   !-   write(debugFile,'(" --> gdm24r: i1,i2=",2i4," q=",8e10.2)') i1,i2,(q(n),n=0,7)
   !- end if
   
   if( debug.gt.7 )then
     write(*,'(" --> gdm24r: i1,i2=",2i4," RHS(A)=",8(1pe13.5))') i1,i2,(f(n),n=0,7)
   end if
   
   ! subtract off the contributions from the initial (wrong) values at the ghost points:
   do n=0,7
     f(n) = (a8(n,0)*q(0)+a8(n,1)*q(1)+a8(n,2)*q(2)+a8(n,3)*q(3)+\
             a8(n,4)*q(4)+a8(n,5)*q(5)+a8(n,6)*q(6)+a8(n,7)*q(7)) - f(n)
   end do


   ! solve A Q = F
   ! factor the matrix
   numberOfEquations=8
   call dgeco( a8(0,0), numberOfEquations, numberOfEquations, ipivot8(0),rcond,work(0))

   if( debug.gt.3 ) then
     write(debugFile,'(" --> gdm24r: i1,i2=",2i4," rcond=",e10.2)') i1,i2,rcond
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

   if( debug>3 .and. twilightZone.eq.1 )then
     ! check errors
     call ogderiv(ep, 0,0,0,0, xy1(i1-is1,i2-is2,i3,0),xy1(i1-is1,i2-is2,i3,1),0.,t, ex, evv(0) )
     call ogderiv(ep, 0,0,0,0, xy1(i1-is1,i2-is2,i3,0),xy1(i1-is1,i2-is2,i3,1),0.,t, ey, evv(1) )
  
     call ogderiv(ep, 0,0,0,0, xy2(j1-js1,j2-js2,j3,0),xy2(j1-js1,j2-js2,j3,1),0.,t, ex, evv(2) )
     call ogderiv(ep, 0,0,0,0, xy2(j1-js1,j2-js2,j3,0),xy2(j1-js1,j2-js2,j3,1),0.,t, ey, evv(3) )
  
     call ogderiv(ep, 0,0,0,0, xy1(i1-2*is1,i2-2*is2,i3,0),xy1(i1-2*is1,i2-2*is2,i3,1),0.,t, ex, evv(4) )
     call ogderiv(ep, 0,0,0,0, xy1(i1-2*is1,i2-2*is2,i3,0),xy1(i1-2*is1,i2-2*is2,i3,1),0.,t, ey, evv(5) )
  
     call ogderiv(ep, 0,0,0,0, xy2(j1-2*js1,j2-2*js2,j3,0),xy2(j1-2*js1,j2-2*js2,j3,1),0.,t, ex, evv(6) )
     call ogderiv(ep, 0,0,0,0, xy2(j1-2*js1,j2-2*js2,j3,0),xy2(j1-2*js1,j2-2*js2,j3,1),0.,t, ey, evv(7) )
  
     write(*,'("gdm24r: Errors in u1(-1), v1(-1), u2(-1), v2(-1), u1(-2), v1(-2), u2(-2), v2(-2)")') 
  
     maxErr=0.
     do n=0,7
       maxErr =max(maxErr,abs(evv(n)-f(n)))
     end do
     write(*,'("gdm24r: i1,i2=",2i4," err= ",8e8.1," -> maxErr=",e8.1)') i1,i2, (abs(evv(n)-f(n)),n=0,7),maxErr
  
  
   end if

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
   evalMagneticDerivs2dOrder4()
   evalMagneticField2dJumpOrder4()

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

 endLoopsMask2d()
 ! =============== end loops =======================
      
 if( useJacobiUpdate.ne.0 )then
   ! Jacobi-update: now fill in values 
   beginLoopsMask2d() 
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
   endLoopsMask2d()
 end if


#endMacro
