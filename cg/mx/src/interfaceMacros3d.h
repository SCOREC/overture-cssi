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

  if( t.le.2*dt .and. debug.gt.1 )then
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


! =====================================================================================================
!  Macro:
!     Initialize various variables and loop bounds
!     Called by interface3d.bf and interface3dOrder4.bf 
! =====================================================================================================
#beginMacro initializeInterfaceVariablesMacro(LABEL)
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
            write(*,'("LABEL:WARNING: orderOfAccuracy=",i2," but numParallelGhost=",i3)') orderOfAccuracy,numParallelGhost
            write(*,'("LABEL: Choose numParallelGhost==3 to make answers match to np=1 results")') 
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
            write(*,'("LABEL: bc is inconsistent")')
            stop 178
          end if
        end if
        if( boundaryCondition1(1,0).le.0 )then ! parallel ghost may only have bc<0 on one side
          nn1b=nn1b+numGhost
          if( boundaryCondition2(1,0).gt.0 )then
            write(*,'("LABEL: bc is inconsistent")')
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
            write(*,'("LABEL: bc is inconsistent")')
            stop 180
          end if
        end if
        if( boundaryCondition1(1,1).le.0 )then
          nn2b=nn2b+numGhost
          if( boundaryCondition2(1,1).gt.0 )then
            write(*,'("LABEL: bc is inconsistent")')
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
            write(*,'("LABEL: bc is inconsistent")')
            stop 182
          end if
        end if
        if( boundaryCondition1(1,2).le.0 )then
          nn3b=nn3b+numGhost
          if( boundaryCondition2(1,2).gt.0 )then
            write(*,'("LABEL: bc is inconsistent")')
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
        write(debugFile,'(" ERROR: LABEL: orderOfExtrapolation<3 ")')
        stop 7716
      end if
      if( orderOfAccuracy.eq.4 .and. orderOfExtrapolation.lt.4 )then
        write(debugFile,'(" ERROR: LABEL: orderOfExtrapolation<4 ")')
        stop 7716
      end if


#endMacro 




! ==============================================================================
! Macro: Set index bounds to include extra ghost in tangential directions
!   We do this for stage I in the new 4th order scheme (for parallel) where we
!   assign the first ghost point to second-order
! ==============================================================================
#beginMacro setIndexBoundsExtraGhost()
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
#endMacro

! ==============================================================================
! Macro: reset index bounds 
! ==============================================================================
#beginMacro resetIndexBounds()
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
#endMacro
      
