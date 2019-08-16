!         -*- mode: F90 -*-
! --------------------------------------------------------------------------
! ******** define some interface macros for 3D *********
! --------------------------------------------------------------------------


! --------------------------------------------------------------------------
!   Apply interface jump conditions, 3D order 2
! --------------------------------------------------------------------------
#beginMacro eval3dJumpOrder2()
 divE1 = u1x+v1y+w1z
 curlE1x = w1y-v1z
 curlE1y = u1z-w1x
 curlE1z = v1x-u1y
 nDotCurlE1=an1*curlE1x+an2*curlE1y+an3*curlE1z
 nDotLapE1 = an1*u1Lap + an2*v1Lap + an3*w1Lap

 divE2 = u2x+v2y+w2z
 curlE2x = w2y-v2z
 curlE2y = u2z-w2x
 curlE2z = v2x-u2y
 nDotCurlE2=an1*curlE2x+an2*curlE2y+an3*curlE2z
 nDotLapE2 = an1*u2Lap + an2*v2Lap + an3*w2Lap


 f(0)=( divE1*an1 + (curlE1x- nDotCurlE1*an1)/mu1 ) - ( divE2*an1 + (curlE2x- nDotCurlE2*an1)/mu2 )
 f(1)=( divE1*an2 + (curlE1y- nDotCurlE1*an2)/mu1 ) - ( divE2*an2 + (curlE2y- nDotCurlE2*an2)/mu2 )
 f(2)=( divE1*an3 + (curlE1z- nDotCurlE1*an3)/mu1 ) - ( divE2*an3 + (curlE2z- nDotCurlE2*an3)/mu2 )


 f(3)=( u1Lap/(epsmu1) + cem1*nDotLapE1*an1 ) - ( u2Lap/(epsmu2) + cem2*nDotLapE2*an1 )
 f(4)=( v1Lap/(epsmu1) + cem1*nDotLapE1*an2 ) - ( v2Lap/(epsmu2) + cem2*nDotLapE2*an2 )
 f(5)=( w1Lap/(epsmu1) + cem1*nDotLapE1*an3 ) - ( w2Lap/(epsmu2) + cem2*nDotLapE2*an3 )

 if( twilightZone.eq.1 )then

   call ogderiv(ep, 0,1,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ex, uex  )
   call ogderiv(ep, 0,0,1,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ex, uey  )
   call ogderiv(ep, 0,0,0,1, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ex, uez  )
   call ogderiv(ep, 0,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ex, uexx )
   call ogderiv(ep, 0,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ex, ueyy )
   call ogderiv(ep, 0,0,0,2, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ex, uezz )

   call ogderiv(ep, 0,1,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ey, vex  )
   call ogderiv(ep, 0,0,1,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ey, vey  )
   call ogderiv(ep, 0,0,0,1, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ey, vez  )
   call ogderiv(ep, 0,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ey, vexx )
   call ogderiv(ep, 0,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ey, veyy )
   call ogderiv(ep, 0,0,0,2, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ey, vezz )

   call ogderiv(ep, 0,1,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ez, wex  )
   call ogderiv(ep, 0,0,1,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ez, wey  )
   call ogderiv(ep, 0,0,0,1, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ez, wez  )
   call ogderiv(ep, 0,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ez, wexx )
   call ogderiv(ep, 0,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ez, weyy )
   call ogderiv(ep, 0,0,0,2, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ez, wezz )

   ueLap = uexx+ueyy+uezz
   veLap = vexx+veyy+vezz
   weLap = wexx+weyy+wezz

   curlEex = wey-vez
   curlEey = uez-wex
   curlEez = vex-uey
   nDotCurlEe=an1*curlEex+an2*curlEey+an3*curlEez
   nDotLapEe=an1*ueLap+an2*veLap+an3*weLap

   f(0)= f(0) - ( (curlEex- nDotCurlEe*an1)*(1./mu1-1./mu2) )
   f(1)= f(1) - ( (curlEey- nDotCurlEe*an2)*(1./mu1-1./mu2) )
   f(2)= f(2) - ( (curlEez- nDotCurlEe*an3)*(1./mu1-1./mu2) )

   f(3)= f(3) - ( ueLap*(1./epsmu1-1./epsmu2) + nDotLapEe*an1*(cem1-cem2) )
   f(4)= f(4) - ( veLap*(1./epsmu1-1./epsmu2) + nDotLapEe*an2*(cem1-cem2) )
   f(5)= f(5) - ( weLap*(1./epsmu1-1./epsmu2) + nDotLapEe*an3*(cem1-cem2) ) 

 end if

#endMacro

! ----------------------------------------------------------------------------------
!  Macro:
!    Evaluate the interface equations for checking the coefficients
! ----------------------------------------------------------------------------------
#beginMacro evalInterfaceEquations23c()
 evalInterfaceDerivatives3d()
 eval3dJumpOrder2()
#endMacro


! --------------------------------------------------------------------------
! Macro: Assign interface ghost values, DIM=3, ORDER=2, GRID=Curvilinear
! 
! --------------------------------------------------------------------------
#beginMacro assignInterfaceGhost23c()

  ! here are the jump conditions for the ghost points
  !   [ div(E) n + (curl(E)- n.curl(E) n )/mu ] =0                 (3 eqns)
  !   [ Lap(E)/(eps*mu) + (1/mu)*(1-1/eps)*( n.Lap(E) ) n ] = 0    (3 eqns)

  ! These correspond to the 6 conditions:
  !   [ div(E) ] =0 
  !   [ tau. curl(E)/mu ] = 0       (2 tangents)
  !   [ n.Lap(E)/mu ] = 0 
  !   [ tau.Lap(E)/(eps*mu) ] = 0   (2 tangents)

  if( t.le.2*dt )then
    write(*,'("assignInterfaceGhost23c ...")')
  end if

  beginLoopsMask3d()

    ! here is the normal (assumed to be the same on both sides)
    an1=rsxy1(i1,i2,i3,axis1,0)   ! normal (an1,an2)
    an2=rsxy1(i1,i2,i3,axis1,1)
    an3=rsxy1(i1,i2,i3,axis1,2)
    aNorm=max(epsx,sqrt(an1**2+an2**2+an3**2))
    an1=an1/aNorm
    an2=an2/aNorm
    an3=an3/aNorm


    ! --- first evaluate the equations we want to solve with the wrong values at the ghost points:

    cem1=(1.-1./eps1)/mu1
    cem2=(1.-1./eps2)/mu2
    ! evalInterfaceDerivatives3d
    evalInterfaceDerivatives3d()
    eval3dJumpOrder2()

    if( debug.gt.4 )then
     write(debugFile,'(" --> 3d-order2-curv: i1,i2,i3=",3i4," f(start)=",6f8.3)') i1,i2,i3,f(0),f(1),f(2),f(3),f(4),f(5)
     ! '
     write(debugFile,'(" --> u1x,u1y,u1z,v1x,v1y,v1z=",6f8.4)') u1x,u1y,u1z,v1x,v1y,v1z
     write(debugFile,'(" --> u2x,u2y,u2z,v2x,v2y,v2z=",6f8.4)') u2x,u2y,u2z,v2x,v2y,v2z

     write(debugFile,'(" --> vv1r,vv1s,vv1t         =",3e9.2)') vv1r,vv1s,vv1t
     do k3=-1,1
     do k2=-1,1
     write(debugFile,'(" --> v1: =",3f8.4)') u1(i1-1,i2+k2,i3+k3,ey),u1(i1,i2+k2,i3+k3,ey),u1(i1+1,i2+k2,i3+k3,ey)
     end do
     end do
     do k3=-1,1
     do k2=-1,1
     write(debugFile,'(" --> v2: =",3f8.4)') u2(j1-1,j2+k2,j3+k3,ey),u2(j1,j2+k2,j3+k3,ey),u2(j1+1,j2+k2,j3+k3,ey)
     end do
     end do
     ! '
    end if

    ! here is the matrix of coefficients for the unknowns u1(-1),v1(-1),w1(-1),  u2(-1),v2(-1),w2(-1)
    ! Solve:
    !     
    !       A [ U ] = A [ U(old) ] - [ f ]

    c1x = -is*rsxy1(i1,i2,i3,axis1,0)/(2.*dr1(axis1))    ! coeff of u1(-1) from D.x
    c1y = -is*rsxy1(i1,i2,i3,axis1,1)/(2.*dr1(axis1))    ! coeff of u1(-1) from D.y 
    c1z = -is*rsxy1(i1,i2,i3,axis1,2)/(2.*dr1(axis1))    ! coeff of u1(-1) from D.z

    c2x = -js*rsxy2(j1,j2,j3,axis2,0)/(2.*dr2(axis2))
    c2y = -js*rsxy2(j1,j2,j3,axis2,1)/(2.*dr2(axis2))
    c2z = -js*rsxy2(j1,j2,j3,axis2,2)/(2.*dr2(axis2))

    rxx1(0,0,0)=aj1rxx
    rxx1(0,1,1)=aj1ryy
    rxx1(0,2,2)=aj1rzz
    rxx1(1,0,0)=aj1sxx
    rxx1(1,1,1)=aj1syy
    rxx1(1,2,2)=aj1szz
    rxx1(2,0,0)=aj1txx
    rxx1(2,1,1)=aj1tyy
    rxx1(2,2,2)=aj1tzz

    rxx2(0,0,0)=aj2rxx
    rxx2(0,1,1)=aj2ryy
    rxx2(0,2,2)=aj2rzz
    rxx2(1,0,0)=aj2sxx
    rxx2(1,1,1)=aj2syy
    rxx2(1,2,2)=aj2szz
    rxx2(2,0,0)=aj2txx
    rxx2(2,1,1)=aj2tyy
    rxx2(2,2,2)=aj2tzz

    ! clap1 : coeff of u(-1) from lap = u.xx + u.yy + u.zz

    ! clap1=(rsxy1(i1,i2,i3,axis1,0)**2+rsxy1(i1,i2,i3,axis1,1)**2)/(dr1(axis1)**2) \
    !           -is*(rsxy1x22(i1,i2,i3,axis1,0)+rsxy1y22(i1,i2,i3,axis1,1))/(2.*dr1(axis1))
    ! clap2=(rsxy2(j1,j2,j3,axis2,0)**2+rsxy2(j1,j2,j3,axis2,1)**2)/(dr2(axis2)**2) \
    !             -js*(rsxy2x22(j1,j2,j3,axis2,0)+rsxy2y22(j1,j2,j3,axis2,1))/(2.*dr2(axis2)) 
    clap1=(rsxy1(i1,i2,i3,axis1,0)**2+rsxy1(i1,i2,i3,axis1,1)**2+rsxy1(i1,i2,i3,axis1,2)**2)/(dr1(axis1)**2) \
              -is*(rxx1(axis1,0,0)+rxx1(axis1,1,1)+rxx1(axis1,2,2))/(2.*dr1(axis1))
    clap2=(rsxy2(j1,j2,j3,axis2,0)**2+rsxy2(j1,j2,j3,axis2,1)**2+rsxy2(j1,j2,j3,axis2,2)**2)/(dr2(axis2)**2) \
              -js*(rxx2(axis2,0,0)+rxx2(axis2,1,1)+rxx2(axis2,2,2))/(2.*dr2(axis2)) 

    ! cdivE1 =  u.c1x + v.c1y + w.c1z
    ! nDotCurlE1 = (w1y-v1z)*an1 + (u1z-w1x)*an2 + (v1x-u1y)*an3

    ! (u.x+v.y+w.z)*an1 + ( w1y-v1z - nDotCurlE1*an1)/mu1
    a6(0,0) = ( c1x*an1 + (         - (c1z*an2-c1y*an3)*an1 )/mu1 ) ! coeff of u1(-1)
    a6(0,1) = ( c1y*an1 + (    -c1z - (c1x*an3-c1z*an1)*an1 )/mu1 ) ! coeff of v1(-1)
    a6(0,2) = ( c1z*an1 + ( c1y     - (c1y*an1-c1x*an2)*an1 )/mu1 ) ! coeff of w1(-1)

    a6(0,3) =-( c2x*an1 + (         - (c2z*an2-c2y*an3)*an1 )/mu2 ) ! coeff of u2(-1)
    a6(0,4) =-( c2y*an1 + (    -c2z - (c2x*an3-c2z*an1)*an1 )/mu2 ) ! coeff of v2(-1)
    a6(0,5) =-( c2z*an1 + ( c2y     - (c2y*an1-c2x*an2)*an1 )/mu2 ) ! coeff of w2(-1)

    ! Bug fixed - July 2, 2019 wdh : Changed "+" to minus "-" in 4 lines below 

    ! (u.x+v.y+w.z)*an2 + ( u1z-w1x - nDotCurlE1*an2)/mu1
    a6(1,0) = ( c1x*an2 + ( c1z     - (c1z*an2-c1y*an3)*an2 )/mu1 ) ! coeff of u1(-1)
 !  a6(1,1) = ( c1y*an2 + (         + (c1x*an3+c1z*an1)*an2 )/mu1 ) ! coeff of v1(-1)
    a6(1,1) = ( c1y*an2 + (         - (c1x*an3-c1z*an1)*an2 )/mu1 ) ! coeff of v1(-1)
    a6(1,2) = ( c1z*an2 + (    -c1x - (c1y*an1-c1x*an2)*an2 )/mu1 ) ! coeff of w1(-1)

    a6(1,3) =-( c2x*an2 + ( c2z     - (c2z*an2-c2y*an3)*an2 )/mu2 ) ! coeff of u2(-1)
!   a6(1,4) =-( c2y*an2 + (         + (c2x*an3+c2z*an1)*an2 )/mu2 ) ! coeff of v2(-1)
    a6(1,4) =-( c2y*an2 + (         - (c2x*an3-c2z*an1)*an2 )/mu2 ) ! coeff of v2(-1)
    a6(1,5) =-( c2z*an2 + (    -c2x - (c2y*an1-c2x*an2)*an2 )/mu2 ) ! coeff of w2(-1)

    ! (u.x+v.y+w.z)*an3 + ( v1x-u1y - nDotCurlE1*an2)/mu1
    a6(2,0) = ( c1x*an3 + (    -c1y - (c1z*an2-c1y*an3)*an3 )/mu1 ) ! coeff of u1(-1)
!   a6(2,1) = ( c1y*an3 + ( c1x     + (c1x*an3+c1z*an1)*an3 )/mu1 ) ! coeff of v1(-1)
    a6(2,1) = ( c1y*an3 + ( c1x     - (c1x*an3-c1z*an1)*an3 )/mu1 ) ! coeff of v1(-1)
    a6(2,2) = ( c1z*an3 + (         - (c1y*an1-c1x*an2)*an3 )/mu1 ) ! coeff of w1(-1)

    a6(2,3) =-( c2x*an3 + (    -c2y - (c2z*an2-c2y*an3)*an3 )/mu2 ) ! coeff of u2(-1)
!   a6(2,4) =-( c2y*an3 + ( c2x     + (c2x*an3+c2z*an1)*an3 )/mu2 ) ! coeff of v2(-1)
    a6(2,4) =-( c2y*an3 + ( c2x     - (c2x*an3-c2z*an1)*an3 )/mu2 ) ! coeff of v2(-1)
    a6(2,5) =-( c2z*an3 + (         - (c2y*an1-c2x*an2)*an3 )/mu2 ) ! coeff of w2(-1)

    !  u1Lap/(epsmu1) + cem1*( an1*u1Lap + an2*v1Lap + an3*w1Lap )*an1
    a6(3,0) = ( clap1/(epsmu1) + cem1*( an1*clap1                         )*an1 ) ! coeff of u1(-1)
    a6(3,1) = (                  cem1*(             an2*clap1             )*an1 )
    a6(3,2) = (                  cem1*(                         an3*clap1 )*an1 )

    a6(3,3) =-( clap2/(epsmu2) + cem2*( an1*clap2                         )*an1 ) ! coeff of u2(-1)
    a6(3,4) =-(                  cem2*(             an2*clap2             )*an1 )
    a6(3,5) =-(                  cem2*(                         an3*clap2 )*an1 )

    !  v1Lap/(epsmu1) + cem1*( an1*u1Lap + an2*v1Lap + an3*w1Lap )*an2
    a6(4,0) = (                  cem1*( an1*clap1                         )*an2 ) ! coeff of u1(-1)
    a6(4,1) = ( clap1/(epsmu1) + cem1*(             an2*clap1             )*an2 )
    a6(4,2) = (                  cem1*(                         an3*clap1 )*an2 )

    a6(4,3) =-(                  cem2*( an1*clap2                         )*an2 ) ! coeff of u2(-1)
    a6(4,4) =-( clap2/(epsmu2) + cem2*(             an2*clap2             )*an2 )
    a6(4,5) =-(                  cem2*(                         an3*clap2 )*an2 )

    !  w1Lap/(epsmu1) + cem1*( an1*u1Lap + an2*v1Lap + an3*w1Lap )*an3
    a6(5,0) = (                  cem1*( an1*clap1                         )*an3 ) ! coeff of u1(-1)
    a6(5,1) = (                  cem1*(             an2*clap1             )*an3 )
    a6(5,2) = ( clap1/(epsmu1) + cem1*(                         an3*clap1 )*an3 )

    a6(5,3) =-(                  cem2*( an1*clap2                         )*an3 ) ! coeff of u2(-1)
    a6(5,4) =-(                  cem2*(             an2*clap2             )*an3 )
    a6(5,5) =-( clap2/(epsmu2) + cem2*(                         an3*clap2 )*an3 )


    q(0) = u1(i1-is1,i2-is2,i3-is3,ex)
    q(1) = u1(i1-is1,i2-is2,i3-is3,ey)
    q(2) = u1(i1-is1,i2-is2,i3-is3,ez)
    q(3) = u2(j1-js1,j2-js2,j3-js3,ex)
    q(4) = u2(j1-js1,j2-js2,j3-js3,ey)
    q(5) = u2(j1-js1,j2-js2,j3-js3,ez)

    ! --- check matrix coefficients by delta function approach ----
    if( checkCoeff.eq.1 )then
      numberOfEquations=6
      checkCoefficients(i1,i2,i3, j1,j2,j3,numberOfEquations,a6,evalInterfaceEquations23c )

      ! **** OVERWRITE EQUATIONS
      !! evalCoefficients(i1,i2,i3, j1,j2,j3,numberOfEquations,a6,evalInterfaceEquations23c )

    end if

    ! subtract off the contributions from the wrong values at the ghost points:
    do n=0,5
      f(n) = (a6(n,0)*q(0)+a6(n,1)*q(1)+a6(n,2)*q(2)+a6(n,3)*q(3)+a6(n,4)*q(4)+a6(n,5)*q(5)) - f(n)
    end do
    ! write(debugFile,'(" --> 3d:order2-c: f(subtract)=",6f8.3)') f(0),f(1),f(2),f(3),f(4),f(5)
    ! solve A Q = F
    ! factor the matrix
    numberOfEquations=6
    call dgeco( a6(0,0), numberOfEquations, numberOfEquations, ipvt(0),rcond,work(0))
    ! solve
    ! write(debugFile,'(" --> 3d:order2-c: rcond=",e10.2)') rcond
    job=0
    call dgesl( a6(0,0), numberOfEquations, numberOfEquations, ipvt(0), f(0), job)
    ! write(debugFile,'(" --> 3d:order2-c: f(solve)=",6f8.3)') f(0),f(1),f(2),f(3),f(4),f(5)
    ! write(debugFile,'(" --> 3d:order2-c:        q=",6f8.3)') q(0),q(1),q(2),q(3),q(4),q(5)

    u1(i1-is1,i2-is2,i3-is3,ex)=f(0)
    u1(i1-is1,i2-is2,i3-is3,ey)=f(1)
    u1(i1-is1,i2-is2,i3-is3,ez)=f(2)
    u2(j1-js1,j2-js2,j3-js3,ex)=f(3)
    u2(j1-js1,j2-js2,j3-js3,ey)=f(4)
    u2(j1-js1,j2-js2,j3-js3,ez)=f(5)

    if( .false. )then
    u1(i1-is1,i2-is2,i3-is3,ex)=q(0)
    u1(i1-is1,i2-is2,i3-is3,ey)=q(1)
    u1(i1-is1,i2-is2,i3-is3,ez)=q(2)
    u2(j1-js1,j2-js2,j3-js3,ex)=q(3)
    u2(j1-js1,j2-js2,j3-js3,ey)=q(4)
    u2(j1-js1,j2-js2,j3-js3,ez)=q(5)
    end if

    if( debug.gt.3 )then ! re-evaluate
     evalInterfaceDerivatives3d()
     eval3dJumpOrder2()
     write(debugFile,'(" --> 3d-order2-c: i1,i2,i3=",3i4," f(re-eval)=",6e10.2)') i1,i2,i3,f(0),f(1),f(2),f(3),f(4),f(5)
       ! '
    end if

  endLoopsMask3d()

  if( checkCoeff.eq.1 )then
    write(*,'("+++++ I23c: check coeff in interface: max(diff) = ",1pe8.2)') coeffDiff
  end if


#endMacro         
