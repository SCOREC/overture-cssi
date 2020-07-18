! This file automatically generated from forcingOpt.bf with bpp.
! *******************************************************************************
!   Add the forcing to Maxwell's equations such as for moving charges
! *******************************************************************************



! This (cubic) ramp has 1-derivative zero at t=0 and t=tba

! This ramp has 3-derivatives zero at t=0 and t=1
! This is from ramp.maple
! r=-84*t**5+35*t**4-20*t**7+70*t**6
! rt=-420*t**4+140*t**3-140*t**6+420*t**5
! rtt=-1680*t**3+420*t**2-840*t**5+2100*t**4
! rttt=-5040*t**2+840*t-4200*t**4+8400*t**3


! This ramp has 4-derivatives zero at t=0 and t=1
! This is from ramp.maple
! r=126*(t)**5-315*(t)**8+70*(t)**9-420*(t)**6+540*(t)**7
! rt=630*(t)**4-2520*(t)**7+630*(t)**8-2520*(t)**5+3780*(t)**6
! rtt=2520*(t)**3-17640*(t)**6+5040*(t)**7-12600*(t)**4+22680*(t)**5
! rttt=7560*(t)**2-105840*(t)**5+35280*(t)**6-50400*(t)**3+113400*(t)**4





! =====================================================================
! Evaluate the Gaussian pulse 
!
!      rho = a*exp( - [beta* | (xv-xpv) - vv*t |]^p ) )
! =====================================================================




! =====================================================================
! Evaluate the Gaussian pulse and its deriatives
!
!      rho = a*exp( - [beta* | (xv-xpv) - vv*t |]^p ) )
!
!  p should be even and positive: 2,4,6,...
! =====================================================================


! ===========================================================================
!     Forcing for a Moving Gaussian Pulse of Charge rho(x,t)
!
!      E_tt = c^2 [ Delta(E) - grad(rho/eps) - mu J_t ]
!      H_tt = c^2 [ Delta(H) + curl( J ) ]
!
!      rho = a*exp( - [beta* | (xv-xpv) - vv*t |]^p ) )
!      J = rho*vv 
!
!  DIM : 2,3
!  ORDER: 2,4
!  GRIDTYPE : curvilinear, rectangular
! ===========================================================================








      subroutine testGaussianPulse()
! =========================================================================
!  Test the Gaussian pulse and derivatives
! =========================================================================
      implicit none

      real xa(0:2),dx(0:2)
      real amplitude,beta,p,xp0,xp1,xp2,vp0,vp1,vp2
      real x,y,z,h,t,h2
      integer nd,i1,i2,i3
      real rhom2,rhom1,rhop1,rhop2,ux,uxx,uxxx

      real uxm1,uxxm1,uxxxm1,uxp1,uxxp1,uxxxp1
      real uxy,uxxy,uxyy,uxz,uxxz,uxzz,uyz,uyyz,uyzz
      real uxt,uxxt,uxtt, uyt,uyyt,uytt

       real x0,x1,x2,q,r2Min,realMin,betap,a0,ap,app,appp
       real rho,rhop,rhopp,rhoppp
       real r2,r2x,r2y,r2z,r2t,r2f,expr,rhox,rhoy,rhoz,rhot
       real r2xx,r2xy,r2yy,r2xz,r2yz,r2zz,r2xt,r2yt,r2zt,r2tt
       real rhoxx,rhoxy,rhoyy,rhoxz,rhoyz,rhozz
       real r2xxx,r2xxy,r2xyy,r2xxz,r2xzz,r2yyy,r2yyz,r2yzz,r2xyz,
     & r2zzz,r2xxt,r2yyt,r2zzt,r2ttt,r2xtt,r2ytt,r2ztt
       real rhoxxx,rhoxxy,rhoxyy,rhoxxz,rhoxzz,rhoyyy,rhoyyz,rhoyzz,
     & rhoxyz,rhozzz,rhoxxt,rhoyyt,rhozzt
       real rhoxtt,rhoytt,rhoztt,rhottt
       real deltaRhoX,deltaRhoY,deltaRhoZ,deltaRhot
       real deltaJ0t,deltaJ1t,deltaJ2t,j0ttt,j1ttt,j2ttt


      amplitude  =1.
      beta       =2.  ! 1/width
      p          =2. ! 4. ! 2.
      xp0        =.3
      xp1        =.4
      xp2        =.5
      vp0        =.25
      vp1        =.35
      vp2        =.15

      t=.1

      xa(0)=xp0 + vp0*t
      xa(1)=xp1 + vp1*t
      xa(2)=xp2 + vp2*t
      dx(0)=.03
      dx(1)=.04
      dx(2)=.05

      h = (1.e-16)**(1./4.)

      realMin=1.e-30  ! should be REAL_MIN

      betap=beta**p
      q=p/2.
      r2Min=realMin*100.

      do nd=2,3

        ! (x,y,z) must be consistent with (i1,i2,i3)
        i1=1
        i2=1
        i3=1
        x = xa(0)+i1*dx(0)
        y = xa(1)+i2*dx(1)
        z = xa(2)+i3*dx(2)


        if( nd.eq.2 )then

             x0 = xa(0)+i1*dx(0) -xp0
             x1 = xa(1)+i2*dx(1) -xp1
             r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2
             ! radx == r*r_x = .5* (r^2)_x 
             r2x = 2.*(x0-vp0*(t))
             r2y = 2.*(x1-vp1*(t))
               r2t = -2.*( ( x0-vp0*(t) )*vp0 + ( x1-vp1*(t) )*vp1 )
             a0=-betap*r2**q
             expr = exp(a0)  ! exp( -(beta*r)^p )
             rho = amplitude*expr
             r2f=max(r2Min,r2)  ! floor value of r2 to avoid division by zero
             ! rhop = [d(rho)/d(r2)]
             ap=a0*q/r2f  ! d(arg)/d(r2) , arg=-betap*r2**q
             rhop = rho*ap
             rhox = rhop*r2x
             rhoy = rhop*r2y
             rhot = rhop*r2t
               app=ap*(q-1.)/r2f
               rhopp = rhop*ap+rho*app
               r2xx=2.
               r2xy=0.
               r2yy=2.
               r2xt=-2.*vp0
               r2yt=-2.*vp1
                 r2tt=2.*(vp0**2+vp1**2)
               rhoxx = rhopp*r2x*r2x + rhop*r2xx
               rhoxy = rhopp*r2x*r2y + rhop*r2xy
               rhoyy = rhopp*r2y*r2y + rhop*r2yy
               appp=app*(q-2.)/r2f
               rhoppp = rhopp*ap +2.*rhop*app+rho*appp
               r2xxx=0.
               r2xxy=0.
               r2xyy=0.
               r2xxt=0.
               r2yyt=0.
               r2ttt=0.
               rhoxxx = rhoppp*r2x*r2x*r2x + rhopp*( r2x*r2xx + r2x*
     & r2xx + r2x*r2xx ) + rhop*r2xxx
               rhoxxy = rhoppp*r2x*r2x*r2y + rhopp*( r2y*r2xx + r2x*
     & r2xy + r2x*r2xy ) + rhop*r2xxy
               rhoxyy = rhoppp*r2x*r2y*r2y + rhopp*( r2y*r2xy + r2y*
     & r2xy + r2x*r2yy ) + rhop*r2xyy
               rhoyyy = rhoppp*r2y*r2y*r2y + rhopp*( r2y*r2yy + r2y*
     & r2yy + r2y*r2yy ) + rhop*r2yyy
               rhoxxt = rhoppp*r2x*r2x*r2t + rhopp*( r2t*r2xx + r2x*
     & r2xt + r2x*r2xt ) + rhop*r2xxt
               rhoyyt = rhoppp*r2y*r2y*r2t + rhopp*( r2t*r2yy + r2y*
     & r2yt + r2y*r2yt ) + rhop*r2yyt
               rhottt = rhoppp*r2t*r2t*r2t + rhopp*( r2t*r2tt + r2t*
     & r2tt + r2t*r2tt ) + rhop*r2ttt
               rhoxtt = rhoppp*r2x*r2t*r2t + rhopp*( r2t*r2xt + r2t*
     & r2xt + r2x*r2tt ) + rhop*r2xtt
               rhoytt = rhoppp*r2y*r2t*r2t + rhopp*( r2t*r2yt + r2t*
     & r2yt + r2y*r2tt ) + rhop*r2ytt

          ! compare derivatives to those from differences
          write(*,'("2d: rho=",e10.4)') rho

               x0 = x -xp0
               x1 = y -xp1
                 r2 = ( x0-vp0*(t-2.*h) )**2 + ( x1-vp1*(t-2.*h) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom2 = amplitude*expr
               x0 = x -xp0
               x1 = y -xp1
                 r2 = ( x0-vp0*(t-h) )**2 + ( x1-vp1*(t-h) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom1 = amplitude*expr
               x0 = x -xp0
               x1 = y -xp1
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rho = amplitude*expr
               x0 = x -xp0
               x1 = y -xp1
                 r2 = ( x0-vp0*(t+h) )**2 + ( x1-vp1*(t+h) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop1 = amplitude*expr
               x0 = x -xp0
               x1 = y -xp1
                 r2 = ( x0-vp0*(t+2.*h) )**2 + ( x1-vp1*(t+2.*h) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop2 = amplitude*expr
            ux = (rhop1-rhom1)/(2.*h)
            uxx = (rhop1-2.*rho+rhom1)/(h**2)
            uxxx = (-2.*(rhop1-rhom1)+(rhop2-rhom2))/(2.*h**3)
           write(*,'("2 d: rhot  =",e10.3," diff=",e10.3," err=",e8.2)
     & ') rhot,ux,abs(rhot-ux)
           write(*,'("2 d: rhottt=",e10.3," diff=",e10.3," err=",e8.2)
     & ') rhottt,uxxx,abs(rhottt-uxxx)
               x0 = x-2.*h -xp0
               x1 = y -xp1
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom2 = amplitude*expr
               x0 = x-h -xp0
               x1 = y -xp1
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom1 = amplitude*expr
               x0 = x -xp0
               x1 = y -xp1
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rho = amplitude*expr
               x0 = x+h -xp0
               x1 = y -xp1
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop1 = amplitude*expr
               x0 = x+2.*h -xp0
               x1 = y -xp1
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop2 = amplitude*expr
            ux = (rhop1-rhom1)/(2.*h)
            uxx = (rhop1-2.*rho+rhom1)/(h**2)
            uxxx = (-2.*(rhop1-rhom1)+(rhop2-rhom2))/(2.*h**3)
           write(*,'("2 d: rhox  =",e10.3," diff=",e10.3," err=",e8.2)
     & ') rhox,ux,abs(rhox-ux)
           write(*,'("2 d: rhoxx =",e10.3," diff=",e10.3," err=",e8.2)
     & ') rhoxx,uxx,abs(rhoxx-uxx)
           write(*,'("2 d: rhoxxx=",e10.3," diff=",e10.3," err=",e8.2)
     & ') rhoxxx,uxxx,abs(rhoxxx-uxxx)
               x0 = x -xp0
               x1 = y-2.*h -xp1
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom2 = amplitude*expr
               x0 = x -xp0
               x1 = y-h -xp1
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom1 = amplitude*expr
               x0 = x -xp0
               x1 = y -xp1
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rho = amplitude*expr
               x0 = x -xp0
               x1 = y+h -xp1
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop1 = amplitude*expr
               x0 = x -xp0
               x1 = y+2.*h -xp1
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop2 = amplitude*expr
            ux = (rhop1-rhom1)/(2.*h)
            uxx = (rhop1-2.*rho+rhom1)/(h**2)
            uxxx = (-2.*(rhop1-rhom1)+(rhop2-rhom2))/(2.*h**3)
           write(*,'("2 d: rhoy  =",e10.3," diff=",e10.3," err=",e8.2)
     & ') rhoy,ux,abs(rhoy-ux)
           write(*,'("2 d: rhoyy =",e10.3," diff=",e10.3," err=",e8.2)
     & ') rhoyy,uxx,abs(rhoyy-uxx)
           write(*,'("2 d: rhoyyy=",e10.3," diff=",e10.3," err=",e8.2)
     & ') rhoyyy,uxxx,abs(rhoyyy-uxxx)
           h2=h
               x0 = x-2.*h -xp0
               x1 = y -xp1
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom2 = amplitude*expr
               x0 = x-h -xp0
               x1 = y -xp1
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom1 = amplitude*expr
               x0 = x -xp0
               x1 = y -xp1
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rho = amplitude*expr
               x0 = x+h -xp0
               x1 = y -xp1
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop1 = amplitude*expr
               x0 = x+2.*h -xp0
               x1 = y -xp1
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop2 = amplitude*expr
            ux = (rhop1-rhom1)/(2.*h)
            uxx = (rhop1-2.*rho+rhom1)/(h**2)
            uxxx = (-2.*(rhop1-rhom1)+(rhop2-rhom2))/(2.*h**3)
               x0 = x-2.*h -xp0
               x1 = y-h2 -xp1
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom2 = amplitude*expr
               x0 = x-h -xp0
               x1 = y-h2 -xp1
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom1 = amplitude*expr
               x0 = x -xp0
               x1 = y-h2 -xp1
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rho = amplitude*expr
               x0 = x+h -xp0
               x1 = y-h2 -xp1
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop1 = amplitude*expr
               x0 = x+2.*h -xp0
               x1 = y-h2 -xp1
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop2 = amplitude*expr
            uxm1 = (rhop1-rhom1)/(2.*h)
            uxxm1 = (rhop1-2.*rho+rhom1)/(h**2)
            uxxxm1 = (-2.*(rhop1-rhom1)+(rhop2-rhom2))/(2.*h**3)
               x0 = x-2.*h -xp0
               x1 = y+h2 -xp1
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom2 = amplitude*expr
               x0 = x-h -xp0
               x1 = y+h2 -xp1
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom1 = amplitude*expr
               x0 = x -xp0
               x1 = y+h2 -xp1
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rho = amplitude*expr
               x0 = x+h -xp0
               x1 = y+h2 -xp1
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop1 = amplitude*expr
               x0 = x+2.*h -xp0
               x1 = y+h2 -xp1
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop2 = amplitude*expr
            uxp1 = (rhop1-rhom1)/(2.*h)
            uxxp1 = (rhop1-2.*rho+rhom1)/(h**2)
            uxxxp1 = (-2.*(rhop1-rhom1)+(rhop2-rhom2))/(2.*h**3)
           uxy=(uxp1-uxm1)/(2.*h2)
           uxxy=(uxxp1-uxxm1)/(2.*h2)
           uxyy=(uxp1-2.*ux+uxm1)/(h2**2)
          !write(*,'("2 d: ux,uxm1,uxp1=",3e12.5') ux,uxm1,uxp1
          !write(*,'("2 d: uxx,uxxm1,uxxp1=",3e12.5') uxx,uxxm1,uxxp1
          !write(*,'("2 d: uxxx,uxxxm1,uxxxp1=",3e12.5') uxxx,uxxxm1,uxxxp1
           write(*,'("2 d: rhoxy =",e10.3," diff=",e10.3," err=",e8.2)
     & ') rhoxy,uxy,abs(rhoxy-uxy)
           write(*,'("2 d: rhoxxy=",e10.3," diff=",e10.3," err=",e8.2)
     & ') rhoxxy,uxxy,abs(rhoxxy-uxxy)
           write(*,'("2 d: rhoxyy=",e10.3," diff=",e10.3," err=",e8.2)
     & ') rhoxyy,uxyy,abs(rhoxyy-uxyy)
               x0 = x-2.*h -xp0
               x1 = y -xp1
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom2 = amplitude*expr
               x0 = x-h -xp0
               x1 = y -xp1
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom1 = amplitude*expr
               x0 = x -xp0
               x1 = y -xp1
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rho = amplitude*expr
               x0 = x+h -xp0
               x1 = y -xp1
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop1 = amplitude*expr
               x0 = x+2.*h -xp0
               x1 = y -xp1
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop2 = amplitude*expr
            ux = (rhop1-rhom1)/(2.*h)
            uxx = (rhop1-2.*rho+rhom1)/(h**2)
            uxxx = (-2.*(rhop1-rhom1)+(rhop2-rhom2))/(2.*h**3)
               x0 = x-2.*h -xp0
               x1 = y -xp1
                 r2 = ( x0-vp0*(t-h2) )**2 + ( x1-vp1*(t-h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom2 = amplitude*expr
               x0 = x-h -xp0
               x1 = y -xp1
                 r2 = ( x0-vp0*(t-h2) )**2 + ( x1-vp1*(t-h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom1 = amplitude*expr
               x0 = x -xp0
               x1 = y -xp1
                 r2 = ( x0-vp0*(t-h2) )**2 + ( x1-vp1*(t-h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rho = amplitude*expr
               x0 = x+h -xp0
               x1 = y -xp1
                 r2 = ( x0-vp0*(t-h2) )**2 + ( x1-vp1*(t-h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop1 = amplitude*expr
               x0 = x+2.*h -xp0
               x1 = y -xp1
                 r2 = ( x0-vp0*(t-h2) )**2 + ( x1-vp1*(t-h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop2 = amplitude*expr
            uxm1 = (rhop1-rhom1)/(2.*h)
            uxxm1 = (rhop1-2.*rho+rhom1)/(h**2)
            uxxxm1 = (-2.*(rhop1-rhom1)+(rhop2-rhom2))/(2.*h**3)
               x0 = x-2.*h -xp0
               x1 = y -xp1
                 r2 = ( x0-vp0*(t+h2) )**2 + ( x1-vp1*(t+h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom2 = amplitude*expr
               x0 = x-h -xp0
               x1 = y -xp1
                 r2 = ( x0-vp0*(t+h2) )**2 + ( x1-vp1*(t+h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom1 = amplitude*expr
               x0 = x -xp0
               x1 = y -xp1
                 r2 = ( x0-vp0*(t+h2) )**2 + ( x1-vp1*(t+h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rho = amplitude*expr
               x0 = x+h -xp0
               x1 = y -xp1
                 r2 = ( x0-vp0*(t+h2) )**2 + ( x1-vp1*(t+h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop1 = amplitude*expr
               x0 = x+2.*h -xp0
               x1 = y -xp1
                 r2 = ( x0-vp0*(t+h2) )**2 + ( x1-vp1*(t+h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop2 = amplitude*expr
            uxp1 = (rhop1-rhom1)/(2.*h)
            uxxp1 = (rhop1-2.*rho+rhom1)/(h**2)
            uxxxp1 = (-2.*(rhop1-rhom1)+(rhop2-rhom2))/(2.*h**3)
           uxxt=(uxxp1-uxxm1)/(2.*h2)
           uxtt=(uxp1-2.*ux+uxm1)/(h2**2)
           write(*,'("2 d: rhoxxt=",e10.3," diff=",e10.3," err=",e8.2)
     & ') rhoxxt,uxxt,abs(rhoxxt-uxxt)
           write(*,'("2 d: rhoxtt=",e10.3," diff=",e10.3," err=",e8.2)
     & ') rhoxtt,uxtt,abs(rhoxtt-uxtt)
               x0 = x -xp0
               x1 = y-2.*h -xp1
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom2 = amplitude*expr
               x0 = x -xp0
               x1 = y-h -xp1
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom1 = amplitude*expr
               x0 = x -xp0
               x1 = y -xp1
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rho = amplitude*expr
               x0 = x -xp0
               x1 = y+h -xp1
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop1 = amplitude*expr
               x0 = x -xp0
               x1 = y+2.*h -xp1
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop2 = amplitude*expr
            ux = (rhop1-rhom1)/(2.*h)
            uxx = (rhop1-2.*rho+rhom1)/(h**2)
            uxxx = (-2.*(rhop1-rhom1)+(rhop2-rhom2))/(2.*h**3)
               x0 = x -xp0
               x1 = y-2.*h -xp1
                 r2 = ( x0-vp0*(t-h2) )**2 + ( x1-vp1*(t-h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom2 = amplitude*expr
               x0 = x -xp0
               x1 = y-h -xp1
                 r2 = ( x0-vp0*(t-h2) )**2 + ( x1-vp1*(t-h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom1 = amplitude*expr
               x0 = x -xp0
               x1 = y -xp1
                 r2 = ( x0-vp0*(t-h2) )**2 + ( x1-vp1*(t-h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rho = amplitude*expr
               x0 = x -xp0
               x1 = y+h -xp1
                 r2 = ( x0-vp0*(t-h2) )**2 + ( x1-vp1*(t-h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop1 = amplitude*expr
               x0 = x -xp0
               x1 = y+2.*h -xp1
                 r2 = ( x0-vp0*(t-h2) )**2 + ( x1-vp1*(t-h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop2 = amplitude*expr
            uxm1 = (rhop1-rhom1)/(2.*h)
            uxxm1 = (rhop1-2.*rho+rhom1)/(h**2)
            uxxxm1 = (-2.*(rhop1-rhom1)+(rhop2-rhom2))/(2.*h**3)
               x0 = x -xp0
               x1 = y-2.*h -xp1
                 r2 = ( x0-vp0*(t+h2) )**2 + ( x1-vp1*(t+h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom2 = amplitude*expr
               x0 = x -xp0
               x1 = y-h -xp1
                 r2 = ( x0-vp0*(t+h2) )**2 + ( x1-vp1*(t+h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom1 = amplitude*expr
               x0 = x -xp0
               x1 = y -xp1
                 r2 = ( x0-vp0*(t+h2) )**2 + ( x1-vp1*(t+h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rho = amplitude*expr
               x0 = x -xp0
               x1 = y+h -xp1
                 r2 = ( x0-vp0*(t+h2) )**2 + ( x1-vp1*(t+h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop1 = amplitude*expr
               x0 = x -xp0
               x1 = y+2.*h -xp1
                 r2 = ( x0-vp0*(t+h2) )**2 + ( x1-vp1*(t+h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop2 = amplitude*expr
            uxp1 = (rhop1-rhom1)/(2.*h)
            uxxp1 = (rhop1-2.*rho+rhom1)/(h**2)
            uxxxp1 = (-2.*(rhop1-rhom1)+(rhop2-rhom2))/(2.*h**3)
           uxxt=(uxxp1-uxxm1)/(2.*h2)
           uxtt=(uxp1-2.*ux+uxm1)/(h2**2)
           write(*,'("2 d: rhoyyt=",e10.3," diff=",e10.3," err=",e8.2)
     & ') rhoyyt,uxxt,abs(rhoyyt-uxxt)
           write(*,'("2 d: rhoytt=",e10.3," diff=",e10.3," err=",e8.2)
     & ') rhoytt,uxtt,abs(rhoytt-uxtt)

        else

             x0 = xa(0)+i1*dx(0) -xp0
             x1 = xa(1)+i2*dx(1) -xp1
               x2 = xa(2)+i3*dx(2) -xp2
             r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*(t) 
     & )**2
             ! radx == r*r_x = .5* (r^2)_x 
             r2x = 2.*(x0-vp0*(t))
             r2y = 2.*(x1-vp1*(t))
               r2z = 2.*(x2-vp2*(t))
               r2t = -2.*( ( x0-vp0*(t) )*vp0 + ( x1-vp1*(t) )*vp1 + ( 
     & x2-vp2*(t) )*vp2 )
             a0=-betap*r2**q
             expr = exp(a0)  ! exp( -(beta*r)^p )
             rho = amplitude*expr
             r2f=max(r2Min,r2)  ! floor value of r2 to avoid division by zero
             ! rhop = [d(rho)/d(r2)]
             ap=a0*q/r2f  ! d(arg)/d(r2) , arg=-betap*r2**q
             rhop = rho*ap
             rhox = rhop*r2x
             rhoy = rhop*r2y
             rhot = rhop*r2t
              rhoz = rhop*r2z
               app=ap*(q-1.)/r2f
               rhopp = rhop*ap+rho*app
               r2xx=2.
               r2xy=0.
               r2yy=2.
               r2xt=-2.*vp0
               r2yt=-2.*vp1
                 r2xz=0.
                 r2yz=0.
                 r2zz=2.
                 r2zt=-2.*vp2
                 r2tt=2.*(vp0**2+vp1**2+vp2**2)
               rhoxx = rhopp*r2x*r2x + rhop*r2xx
               rhoxy = rhopp*r2x*r2y + rhop*r2xy
               rhoyy = rhopp*r2y*r2y + rhop*r2yy
                 rhoxz = rhopp*r2x*r2z + rhop*r2xz
                 rhoyz = rhopp*r2y*r2z + rhop*r2yz
                 rhozz = rhopp*r2z*r2z + rhop*r2zz
               appp=app*(q-2.)/r2f
               rhoppp = rhopp*ap +2.*rhop*app+rho*appp
               r2xxx=0.
               r2xxy=0.
               r2xyy=0.
               r2xxt=0.
               r2yyt=0.
               r2ttt=0.
               rhoxxx = rhoppp*r2x*r2x*r2x + rhopp*( r2x*r2xx + r2x*
     & r2xx + r2x*r2xx ) + rhop*r2xxx
               rhoxxy = rhoppp*r2x*r2x*r2y + rhopp*( r2y*r2xx + r2x*
     & r2xy + r2x*r2xy ) + rhop*r2xxy
               rhoxyy = rhoppp*r2x*r2y*r2y + rhopp*( r2y*r2xy + r2y*
     & r2xy + r2x*r2yy ) + rhop*r2xyy
               rhoyyy = rhoppp*r2y*r2y*r2y + rhopp*( r2y*r2yy + r2y*
     & r2yy + r2y*r2yy ) + rhop*r2yyy
               rhoxxt = rhoppp*r2x*r2x*r2t + rhopp*( r2t*r2xx + r2x*
     & r2xt + r2x*r2xt ) + rhop*r2xxt
               rhoyyt = rhoppp*r2y*r2y*r2t + rhopp*( r2t*r2yy + r2y*
     & r2yt + r2y*r2yt ) + rhop*r2yyt
               rhottt = rhoppp*r2t*r2t*r2t + rhopp*( r2t*r2tt + r2t*
     & r2tt + r2t*r2tt ) + rhop*r2ttt
               rhoxtt = rhoppp*r2x*r2t*r2t + rhopp*( r2t*r2xt + r2t*
     & r2xt + r2x*r2tt ) + rhop*r2xtt
               rhoytt = rhoppp*r2y*r2t*r2t + rhopp*( r2t*r2yt + r2t*
     & r2yt + r2y*r2tt ) + rhop*r2ytt
                 r2xxz=0.
                 r2xzz=0.
                 r2yyz=0.
                 r2yzz=0.
                 r2zzz=0.
                 r2zzt=0.
                 rhoxxz = rhoppp*r2x*r2x*r2z + rhopp*( r2z*r2xx + r2x*
     & r2xz + r2x*r2xz ) + rhop*r2xxz
                 rhoxzz = rhoppp*r2x*r2z*r2z + rhopp*( r2z*r2xz + r2z*
     & r2xz + r2x*r2zz ) + rhop*r2xzz
                 rhoyzz = rhoppp*r2y*r2z*r2z + rhopp*( r2z*r2yz + r2z*
     & r2yz + r2y*r2zz ) + rhop*r2yzz
                 rhoyyz = rhoppp*r2y*r2y*r2z + rhopp*( r2z*r2yy + r2y*
     & r2yz + r2y*r2yz ) + rhop*r2yyz
                 rhozzz = rhoppp*r2z*r2z*r2z + rhopp*( r2z*r2zz + r2z*
     & r2zz + r2z*r2zz ) + rhop*r2zzz
                 rhozzt = rhoppp*r2z*r2z*r2t + rhopp*( r2t*r2zz + r2z*
     & r2zt + r2z*r2zt ) + rhop*r2zzt
                 rhoztt = rhoppp*r2z*r2t*r2t + rhopp*( r2t*r2zt + r2t*
     & r2zt + r2z*r2tt ) + rhop*r2ztt

          ! compare derivatives to those from differences

               x0 = x -xp0
               x1 = y -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t-2.*h) )**2 + ( x1-vp1*(t-2.*h) )**2 +
     &  ( x2-vp2*(t-2.*h) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom2 = amplitude*expr
               x0 = x -xp0
               x1 = y -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t-h) )**2 + ( x1-vp1*(t-h) )**2 + ( x2-
     & vp2*(t-h) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom1 = amplitude*expr
               x0 = x -xp0
               x1 = y -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rho = amplitude*expr
               x0 = x -xp0
               x1 = y -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t+h) )**2 + ( x1-vp1*(t+h) )**2 + ( x2-
     & vp2*(t+h) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop1 = amplitude*expr
               x0 = x -xp0
               x1 = y -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t+2.*h) )**2 + ( x1-vp1*(t+2.*h) )**2 +
     &  ( x2-vp2*(t+2.*h) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop2 = amplitude*expr
            ux = (rhop1-rhom1)/(2.*h)
            uxx = (rhop1-2.*rho+rhom1)/(h**2)
            uxxx = (-2.*(rhop1-rhom1)+(rhop2-rhom2))/(2.*h**3)
           write(*,'("3 d: rhot  =",e10.3," diff=",e10.3," err=",e8.2)
     & ') rhot,ux,abs(rhot-ux)
           write(*,'("3 d: rhottt=",e10.3," diff=",e10.3," err=",e8.2)
     & ') rhottt,uxxx,abs(rhottt-uxxx)
               x0 = x-2.*h -xp0
               x1 = y -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom2 = amplitude*expr
               x0 = x-h -xp0
               x1 = y -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom1 = amplitude*expr
               x0 = x -xp0
               x1 = y -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rho = amplitude*expr
               x0 = x+h -xp0
               x1 = y -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop1 = amplitude*expr
               x0 = x+2.*h -xp0
               x1 = y -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop2 = amplitude*expr
            ux = (rhop1-rhom1)/(2.*h)
            uxx = (rhop1-2.*rho+rhom1)/(h**2)
            uxxx = (-2.*(rhop1-rhom1)+(rhop2-rhom2))/(2.*h**3)
           write(*,'("3 d: rhox  =",e10.3," diff=",e10.3," err=",e8.2)
     & ') rhox,ux,abs(rhox-ux)
           write(*,'("3 d: rhoxx =",e10.3," diff=",e10.3," err=",e8.2)
     & ') rhoxx,uxx,abs(rhoxx-uxx)
           write(*,'("3 d: rhoxxx=",e10.3," diff=",e10.3," err=",e8.2)
     & ') rhoxxx,uxxx,abs(rhoxxx-uxxx)
               x0 = x -xp0
               x1 = y-2.*h -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom2 = amplitude*expr
               x0 = x -xp0
               x1 = y-h -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom1 = amplitude*expr
               x0 = x -xp0
               x1 = y -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rho = amplitude*expr
               x0 = x -xp0
               x1 = y+h -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop1 = amplitude*expr
               x0 = x -xp0
               x1 = y+2.*h -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop2 = amplitude*expr
            ux = (rhop1-rhom1)/(2.*h)
            uxx = (rhop1-2.*rho+rhom1)/(h**2)
            uxxx = (-2.*(rhop1-rhom1)+(rhop2-rhom2))/(2.*h**3)
           write(*,'("3 d: rhoy  =",e10.3," diff=",e10.3," err=",e8.2)
     & ') rhoy,ux,abs(rhoy-ux)
           write(*,'("3 d: rhoyy =",e10.3," diff=",e10.3," err=",e8.2)
     & ') rhoyy,uxx,abs(rhoyy-uxx)
           write(*,'("3 d: rhoyyy=",e10.3," diff=",e10.3," err=",e8.2)
     & ') rhoyyy,uxxx,abs(rhoyyy-uxxx)
               x0 = x -xp0
               x1 = y -xp1
                 x2 = z-2.*h -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom2 = amplitude*expr
               x0 = x -xp0
               x1 = y -xp1
                 x2 = z-h -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom1 = amplitude*expr
               x0 = x -xp0
               x1 = y -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rho = amplitude*expr
               x0 = x -xp0
               x1 = y -xp1
                 x2 = z+h -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop1 = amplitude*expr
               x0 = x -xp0
               x1 = y -xp1
                 x2 = z+2.*h -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop2 = amplitude*expr
            ux = (rhop1-rhom1)/(2.*h)
            uxx = (rhop1-2.*rho+rhom1)/(h**2)
            uxxx = (-2.*(rhop1-rhom1)+(rhop2-rhom2))/(2.*h**3)
           write(*,'("3 d: rhoz  =",e10.3," diff=",e10.3," err=",e8.2)
     & ') rhoz,ux,abs(rhoz-ux)
           write(*,'("3 d: rhozz =",e10.3," diff=",e10.3," err=",e8.2)
     & ') rhozz,uxx,abs(rhozz-uxx)
           write(*,'("3 d: rhozzz=",e10.3," diff=",e10.3," err=",e8.2)
     & ') rhozzz,uxxx,abs(rhozzz-uxxx)
           h2=h
               x0 = x-2.*h -xp0
               x1 = y -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom2 = amplitude*expr
               x0 = x-h -xp0
               x1 = y -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom1 = amplitude*expr
               x0 = x -xp0
               x1 = y -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rho = amplitude*expr
               x0 = x+h -xp0
               x1 = y -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop1 = amplitude*expr
               x0 = x+2.*h -xp0
               x1 = y -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop2 = amplitude*expr
            ux = (rhop1-rhom1)/(2.*h)
            uxx = (rhop1-2.*rho+rhom1)/(h**2)
            uxxx = (-2.*(rhop1-rhom1)+(rhop2-rhom2))/(2.*h**3)
               x0 = x-2.*h -xp0
               x1 = y-h2 -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom2 = amplitude*expr
               x0 = x-h -xp0
               x1 = y-h2 -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom1 = amplitude*expr
               x0 = x -xp0
               x1 = y-h2 -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rho = amplitude*expr
               x0 = x+h -xp0
               x1 = y-h2 -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop1 = amplitude*expr
               x0 = x+2.*h -xp0
               x1 = y-h2 -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop2 = amplitude*expr
            uxm1 = (rhop1-rhom1)/(2.*h)
            uxxm1 = (rhop1-2.*rho+rhom1)/(h**2)
            uxxxm1 = (-2.*(rhop1-rhom1)+(rhop2-rhom2))/(2.*h**3)
               x0 = x-2.*h -xp0
               x1 = y+h2 -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom2 = amplitude*expr
               x0 = x-h -xp0
               x1 = y+h2 -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom1 = amplitude*expr
               x0 = x -xp0
               x1 = y+h2 -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rho = amplitude*expr
               x0 = x+h -xp0
               x1 = y+h2 -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop1 = amplitude*expr
               x0 = x+2.*h -xp0
               x1 = y+h2 -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop2 = amplitude*expr
            uxp1 = (rhop1-rhom1)/(2.*h)
            uxxp1 = (rhop1-2.*rho+rhom1)/(h**2)
            uxxxp1 = (-2.*(rhop1-rhom1)+(rhop2-rhom2))/(2.*h**3)
           uxy=(uxp1-uxm1)/(2.*h2)
           uxxy=(uxxp1-uxxm1)/(2.*h2)
           uxyy=(uxp1-2.*ux+uxm1)/(h2**2)
          !write(*,'("3 d: ux,uxm1,uxp1=",3e12.5') ux,uxm1,uxp1
          !write(*,'("3 d: uxx,uxxm1,uxxp1=",3e12.5') uxx,uxxm1,uxxp1
          !write(*,'("3 d: uxxx,uxxxm1,uxxxp1=",3e12.5') uxxx,uxxxm1,uxxxp1
           write(*,'("3 d: rhoxy =",e10.3," diff=",e10.3," err=",e8.2)
     & ') rhoxy,uxy,abs(rhoxy-uxy)
           write(*,'("3 d: rhoxxy=",e10.3," diff=",e10.3," err=",e8.2)
     & ') rhoxxy,uxxy,abs(rhoxxy-uxxy)
           write(*,'("3 d: rhoxyy=",e10.3," diff=",e10.3," err=",e8.2)
     & ') rhoxyy,uxyy,abs(rhoxyy-uxyy)
               x0 = x-2.*h -xp0
               x1 = y -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom2 = amplitude*expr
               x0 = x-h -xp0
               x1 = y -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom1 = amplitude*expr
               x0 = x -xp0
               x1 = y -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rho = amplitude*expr
               x0 = x+h -xp0
               x1 = y -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop1 = amplitude*expr
               x0 = x+2.*h -xp0
               x1 = y -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop2 = amplitude*expr
            ux = (rhop1-rhom1)/(2.*h)
            uxx = (rhop1-2.*rho+rhom1)/(h**2)
            uxxx = (-2.*(rhop1-rhom1)+(rhop2-rhom2))/(2.*h**3)
               x0 = x-2.*h -xp0
               x1 = y -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t-h2) )**2 + ( x1-vp1*(t-h2) )**2 + ( 
     & x2-vp2*(t-h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom2 = amplitude*expr
               x0 = x-h -xp0
               x1 = y -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t-h2) )**2 + ( x1-vp1*(t-h2) )**2 + ( 
     & x2-vp2*(t-h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom1 = amplitude*expr
               x0 = x -xp0
               x1 = y -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t-h2) )**2 + ( x1-vp1*(t-h2) )**2 + ( 
     & x2-vp2*(t-h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rho = amplitude*expr
               x0 = x+h -xp0
               x1 = y -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t-h2) )**2 + ( x1-vp1*(t-h2) )**2 + ( 
     & x2-vp2*(t-h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop1 = amplitude*expr
               x0 = x+2.*h -xp0
               x1 = y -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t-h2) )**2 + ( x1-vp1*(t-h2) )**2 + ( 
     & x2-vp2*(t-h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop2 = amplitude*expr
            uxm1 = (rhop1-rhom1)/(2.*h)
            uxxm1 = (rhop1-2.*rho+rhom1)/(h**2)
            uxxxm1 = (-2.*(rhop1-rhom1)+(rhop2-rhom2))/(2.*h**3)
               x0 = x-2.*h -xp0
               x1 = y -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t+h2) )**2 + ( x1-vp1*(t+h2) )**2 + ( 
     & x2-vp2*(t+h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom2 = amplitude*expr
               x0 = x-h -xp0
               x1 = y -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t+h2) )**2 + ( x1-vp1*(t+h2) )**2 + ( 
     & x2-vp2*(t+h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom1 = amplitude*expr
               x0 = x -xp0
               x1 = y -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t+h2) )**2 + ( x1-vp1*(t+h2) )**2 + ( 
     & x2-vp2*(t+h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rho = amplitude*expr
               x0 = x+h -xp0
               x1 = y -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t+h2) )**2 + ( x1-vp1*(t+h2) )**2 + ( 
     & x2-vp2*(t+h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop1 = amplitude*expr
               x0 = x+2.*h -xp0
               x1 = y -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t+h2) )**2 + ( x1-vp1*(t+h2) )**2 + ( 
     & x2-vp2*(t+h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop2 = amplitude*expr
            uxp1 = (rhop1-rhom1)/(2.*h)
            uxxp1 = (rhop1-2.*rho+rhom1)/(h**2)
            uxxxp1 = (-2.*(rhop1-rhom1)+(rhop2-rhom2))/(2.*h**3)
           uxxt=(uxxp1-uxxm1)/(2.*h2)
           uxtt=(uxp1-2.*ux+uxm1)/(h2**2)
           write(*,'("3 d: rhoxxt=",e10.3," diff=",e10.3," err=",e8.2)
     & ') rhoxxt,uxxt,abs(rhoxxt-uxxt)
           write(*,'("3 d: rhoxtt=",e10.3," diff=",e10.3," err=",e8.2)
     & ') rhoxtt,uxtt,abs(rhoxtt-uxtt)
               x0 = x -xp0
               x1 = y-2.*h -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom2 = amplitude*expr
               x0 = x -xp0
               x1 = y-h -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom1 = amplitude*expr
               x0 = x -xp0
               x1 = y -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rho = amplitude*expr
               x0 = x -xp0
               x1 = y+h -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop1 = amplitude*expr
               x0 = x -xp0
               x1 = y+2.*h -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop2 = amplitude*expr
            ux = (rhop1-rhom1)/(2.*h)
            uxx = (rhop1-2.*rho+rhom1)/(h**2)
            uxxx = (-2.*(rhop1-rhom1)+(rhop2-rhom2))/(2.*h**3)
               x0 = x -xp0
               x1 = y-2.*h -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t-h2) )**2 + ( x1-vp1*(t-h2) )**2 + ( 
     & x2-vp2*(t-h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom2 = amplitude*expr
               x0 = x -xp0
               x1 = y-h -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t-h2) )**2 + ( x1-vp1*(t-h2) )**2 + ( 
     & x2-vp2*(t-h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom1 = amplitude*expr
               x0 = x -xp0
               x1 = y -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t-h2) )**2 + ( x1-vp1*(t-h2) )**2 + ( 
     & x2-vp2*(t-h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rho = amplitude*expr
               x0 = x -xp0
               x1 = y+h -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t-h2) )**2 + ( x1-vp1*(t-h2) )**2 + ( 
     & x2-vp2*(t-h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop1 = amplitude*expr
               x0 = x -xp0
               x1 = y+2.*h -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t-h2) )**2 + ( x1-vp1*(t-h2) )**2 + ( 
     & x2-vp2*(t-h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop2 = amplitude*expr
            uxm1 = (rhop1-rhom1)/(2.*h)
            uxxm1 = (rhop1-2.*rho+rhom1)/(h**2)
            uxxxm1 = (-2.*(rhop1-rhom1)+(rhop2-rhom2))/(2.*h**3)
               x0 = x -xp0
               x1 = y-2.*h -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t+h2) )**2 + ( x1-vp1*(t+h2) )**2 + ( 
     & x2-vp2*(t+h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom2 = amplitude*expr
               x0 = x -xp0
               x1 = y-h -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t+h2) )**2 + ( x1-vp1*(t+h2) )**2 + ( 
     & x2-vp2*(t+h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom1 = amplitude*expr
               x0 = x -xp0
               x1 = y -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t+h2) )**2 + ( x1-vp1*(t+h2) )**2 + ( 
     & x2-vp2*(t+h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rho = amplitude*expr
               x0 = x -xp0
               x1 = y+h -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t+h2) )**2 + ( x1-vp1*(t+h2) )**2 + ( 
     & x2-vp2*(t+h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop1 = amplitude*expr
               x0 = x -xp0
               x1 = y+2.*h -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t+h2) )**2 + ( x1-vp1*(t+h2) )**2 + ( 
     & x2-vp2*(t+h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop2 = amplitude*expr
            uxp1 = (rhop1-rhom1)/(2.*h)
            uxxp1 = (rhop1-2.*rho+rhom1)/(h**2)
            uxxxp1 = (-2.*(rhop1-rhom1)+(rhop2-rhom2))/(2.*h**3)
           uxxt=(uxxp1-uxxm1)/(2.*h2)
           uxtt=(uxp1-2.*ux+uxm1)/(h2**2)
           write(*,'("3 d: rhoyyt=",e10.3," diff=",e10.3," err=",e8.2)
     & ') rhoyyt,uxxt,abs(rhoyyt-uxxt)
           write(*,'("3 d: rhoytt=",e10.3," diff=",e10.3," err=",e8.2)
     & ') rhoytt,uxtt,abs(rhoytt-uxtt)
               x0 = x-2.*h -xp0
               x1 = y -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom2 = amplitude*expr
               x0 = x-h -xp0
               x1 = y -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom1 = amplitude*expr
               x0 = x -xp0
               x1 = y -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rho = amplitude*expr
               x0 = x+h -xp0
               x1 = y -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop1 = amplitude*expr
               x0 = x+2.*h -xp0
               x1 = y -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop2 = amplitude*expr
            ux = (rhop1-rhom1)/(2.*h)
            uxx = (rhop1-2.*rho+rhom1)/(h**2)
            uxxx = (-2.*(rhop1-rhom1)+(rhop2-rhom2))/(2.*h**3)
               x0 = x-2.*h -xp0
               x1 = y -xp1
                 x2 = z-h2 -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom2 = amplitude*expr
               x0 = x-h -xp0
               x1 = y -xp1
                 x2 = z-h2 -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom1 = amplitude*expr
               x0 = x -xp0
               x1 = y -xp1
                 x2 = z-h2 -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rho = amplitude*expr
               x0 = x+h -xp0
               x1 = y -xp1
                 x2 = z-h2 -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop1 = amplitude*expr
               x0 = x+2.*h -xp0
               x1 = y -xp1
                 x2 = z-h2 -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop2 = amplitude*expr
            uxm1 = (rhop1-rhom1)/(2.*h)
            uxxm1 = (rhop1-2.*rho+rhom1)/(h**2)
            uxxxm1 = (-2.*(rhop1-rhom1)+(rhop2-rhom2))/(2.*h**3)
               x0 = x-2.*h -xp0
               x1 = y -xp1
                 x2 = z+h2 -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom2 = amplitude*expr
               x0 = x-h -xp0
               x1 = y -xp1
                 x2 = z+h2 -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom1 = amplitude*expr
               x0 = x -xp0
               x1 = y -xp1
                 x2 = z+h2 -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rho = amplitude*expr
               x0 = x+h -xp0
               x1 = y -xp1
                 x2 = z+h2 -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop1 = amplitude*expr
               x0 = x+2.*h -xp0
               x1 = y -xp1
                 x2 = z+h2 -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop2 = amplitude*expr
            uxp1 = (rhop1-rhom1)/(2.*h)
            uxxp1 = (rhop1-2.*rho+rhom1)/(h**2)
            uxxxp1 = (-2.*(rhop1-rhom1)+(rhop2-rhom2))/(2.*h**3)
           uxz=(uxp1-uxm1)/(2.*h2)
           uxxz=(uxxp1-uxxm1)/(2.*h2)
           uxzz=(uxp1-2.*ux+uxm1)/(h2**2)
           write(*,'("3 d: rhoxz =",e10.3," diff=",e10.3," err=",e8.2)
     & ') rhoxz,uxz,abs(rhoxz-uxz)
           write(*,'("3 d: rhoxxz=",e10.3," diff=",e10.3," err=",e8.2)
     & ') rhoxxz,uxxz,abs(rhoxxz-uxxz)
           write(*,'("3 d: rhoxzz=",e10.3," diff=",e10.3," err=",e8.2)
     & ') rhoxzz,uxzz,abs(rhoxzz-uxzz)
               x0 = x -xp0
               x1 = y-2.*h -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom2 = amplitude*expr
               x0 = x -xp0
               x1 = y-h -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom1 = amplitude*expr
               x0 = x -xp0
               x1 = y -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rho = amplitude*expr
               x0 = x -xp0
               x1 = y+h -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop1 = amplitude*expr
               x0 = x -xp0
               x1 = y+2.*h -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop2 = amplitude*expr
            ux = (rhop1-rhom1)/(2.*h)
            uxx = (rhop1-2.*rho+rhom1)/(h**2)
            uxxx = (-2.*(rhop1-rhom1)+(rhop2-rhom2))/(2.*h**3)
               x0 = x -xp0
               x1 = y-2.*h -xp1
                 x2 = z-h2 -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom2 = amplitude*expr
               x0 = x -xp0
               x1 = y-h -xp1
                 x2 = z-h2 -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom1 = amplitude*expr
               x0 = x -xp0
               x1 = y -xp1
                 x2 = z-h2 -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rho = amplitude*expr
               x0 = x -xp0
               x1 = y+h -xp1
                 x2 = z-h2 -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop1 = amplitude*expr
               x0 = x -xp0
               x1 = y+2.*h -xp1
                 x2 = z-h2 -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop2 = amplitude*expr
            uxm1 = (rhop1-rhom1)/(2.*h)
            uxxm1 = (rhop1-2.*rho+rhom1)/(h**2)
            uxxxm1 = (-2.*(rhop1-rhom1)+(rhop2-rhom2))/(2.*h**3)
               x0 = x -xp0
               x1 = y-2.*h -xp1
                 x2 = z+h2 -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom2 = amplitude*expr
               x0 = x -xp0
               x1 = y-h -xp1
                 x2 = z+h2 -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom1 = amplitude*expr
               x0 = x -xp0
               x1 = y -xp1
                 x2 = z+h2 -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rho = amplitude*expr
               x0 = x -xp0
               x1 = y+h -xp1
                 x2 = z+h2 -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop1 = amplitude*expr
               x0 = x -xp0
               x1 = y+2.*h -xp1
                 x2 = z+h2 -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop2 = amplitude*expr
            uxp1 = (rhop1-rhom1)/(2.*h)
            uxxp1 = (rhop1-2.*rho+rhom1)/(h**2)
            uxxxp1 = (-2.*(rhop1-rhom1)+(rhop2-rhom2))/(2.*h**3)
           uyz=(uxp1-uxm1)/(2.*h2)
           uyyz=(uxxp1-uxxm1)/(2.*h2)
           uyzz=(uxp1-2.*ux+uxm1)/(h2**2)
           write(*,'("3 d: rhoyz =",e10.3," diff=",e10.3," err=",e8.2)
     & ') rhoyz,uyz,abs(rhoyz-uyz)
           write(*,'("3 d: rhoyyz=",e10.3," diff=",e10.3," err=",e8.2)
     & ') rhoyyz,uyyz,abs(rhoyyz-uyyz)
           write(*,'("3 d: rhoyzz=",e10.3," diff=",e10.3," err=",e8.2)
     & ') rhoyzz,uyzz,abs(rhoyzz-uyzz)
               x0 = x -xp0
               x1 = y -xp1
                 x2 = z-2.*h -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom2 = amplitude*expr
               x0 = x -xp0
               x1 = y -xp1
                 x2 = z-h -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom1 = amplitude*expr
               x0 = x -xp0
               x1 = y -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rho = amplitude*expr
               x0 = x -xp0
               x1 = y -xp1
                 x2 = z+h -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop1 = amplitude*expr
               x0 = x -xp0
               x1 = y -xp1
                 x2 = z+2.*h -xp2
                 r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*
     & (t) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop2 = amplitude*expr
            ux = (rhop1-rhom1)/(2.*h)
            uxx = (rhop1-2.*rho+rhom1)/(h**2)
            uxxx = (-2.*(rhop1-rhom1)+(rhop2-rhom2))/(2.*h**3)
               x0 = x -xp0
               x1 = y -xp1
                 x2 = z-2.*h -xp2
                 r2 = ( x0-vp0*(t-h2) )**2 + ( x1-vp1*(t-h2) )**2 + ( 
     & x2-vp2*(t-h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom2 = amplitude*expr
               x0 = x -xp0
               x1 = y -xp1
                 x2 = z-h -xp2
                 r2 = ( x0-vp0*(t-h2) )**2 + ( x1-vp1*(t-h2) )**2 + ( 
     & x2-vp2*(t-h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom1 = amplitude*expr
               x0 = x -xp0
               x1 = y -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t-h2) )**2 + ( x1-vp1*(t-h2) )**2 + ( 
     & x2-vp2*(t-h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rho = amplitude*expr
               x0 = x -xp0
               x1 = y -xp1
                 x2 = z+h -xp2
                 r2 = ( x0-vp0*(t-h2) )**2 + ( x1-vp1*(t-h2) )**2 + ( 
     & x2-vp2*(t-h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop1 = amplitude*expr
               x0 = x -xp0
               x1 = y -xp1
                 x2 = z+2.*h -xp2
                 r2 = ( x0-vp0*(t-h2) )**2 + ( x1-vp1*(t-h2) )**2 + ( 
     & x2-vp2*(t-h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop2 = amplitude*expr
            uxm1 = (rhop1-rhom1)/(2.*h)
            uxxm1 = (rhop1-2.*rho+rhom1)/(h**2)
            uxxxm1 = (-2.*(rhop1-rhom1)+(rhop2-rhom2))/(2.*h**3)
               x0 = x -xp0
               x1 = y -xp1
                 x2 = z-2.*h -xp2
                 r2 = ( x0-vp0*(t+h2) )**2 + ( x1-vp1*(t+h2) )**2 + ( 
     & x2-vp2*(t+h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom2 = amplitude*expr
               x0 = x -xp0
               x1 = y -xp1
                 x2 = z-h -xp2
                 r2 = ( x0-vp0*(t+h2) )**2 + ( x1-vp1*(t+h2) )**2 + ( 
     & x2-vp2*(t+h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhom1 = amplitude*expr
               x0 = x -xp0
               x1 = y -xp1
                 x2 = z -xp2
                 r2 = ( x0-vp0*(t+h2) )**2 + ( x1-vp1*(t+h2) )**2 + ( 
     & x2-vp2*(t+h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rho = amplitude*expr
               x0 = x -xp0
               x1 = y -xp1
                 x2 = z+h -xp2
                 r2 = ( x0-vp0*(t+h2) )**2 + ( x1-vp1*(t+h2) )**2 + ( 
     & x2-vp2*(t+h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop1 = amplitude*expr
               x0 = x -xp0
               x1 = y -xp1
                 x2 = z+2.*h -xp2
                 r2 = ( x0-vp0*(t+h2) )**2 + ( x1-vp1*(t+h2) )**2 + ( 
     & x2-vp2*(t+h2) )**2
               expr = exp(-betap*r2**q)  ! exp( -(beta*r)^p )
               rhop2 = amplitude*expr
            uxp1 = (rhop1-rhom1)/(2.*h)
            uxxp1 = (rhop1-2.*rho+rhom1)/(h**2)
            uxxxp1 = (-2.*(rhop1-rhom1)+(rhop2-rhom2))/(2.*h**3)
           uxxt=(uxxp1-uxxm1)/(2.*h2)
           uxtt=(uxp1-2.*ux+uxm1)/(h2**2)
           write(*,'("3 d: rhozzt=",e10.3," diff=",e10.3," err=",e8.2)
     & ') rhozzt,uxxt,abs(rhozzt-uxxt)
           write(*,'("3 d: rhoztt=",e10.3," diff=",e10.3," err=",e8.2)
     & ') rhoztt,uxtt,abs(rhoztt-uxtt)

        end if

      end do

      return
      end


! ===========================================================================
!     Gaussian Source Forcing 
!
!      E_tt = c^2 Delta(E) + f_E 
!      H_tt = c^2 Delta(H) + f_H 
!
!      phi = amp*exp( - beta* | (xv-xv0) |^2 )
!      phiE = sin(2*pi*omega t)* phi 
!      phiH = cos(2*pi*omega t)* phi
!
!   2D:
!       fE1 = -( y-y0 ) phiE  
!       fE2 =  ( x-y0 ) phiE   
!       fH3 =  phiH                       **what should this be ?   
!   3D:
!       fE1 = [ (z-z0) - (y-y0) ] phiE       
!       fE2 = [ (x-x0) - (z-z0) ] phiE       
!       fE3 = [ (y-y0) - (x-x0) ] phiE       
!
!  DIM : 2,3
!  ORDER: 2,4
!  GRIDTYPE : curvilinear, rectangular
! ===========================================================================


      subroutine forcingOptMaxwell( nd, nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,
     & ndf1a,ndf1b,ndf2a,ndf2b,ndf3a,ndf3b,u,f,mask,rsxy, xy,ipar, 
     & rpar, ierr )
! ===================================================================================
!  Optimised Forcing functions for Maxwell's Equations.
!
!   Add the source terms for changes to the forcing RHS f.
!
!             f <- f + (charge source terms)
!
! ===================================================================================

      implicit none

      integer nd, nd1a,nd1b,nd2a,nd2b,nd3a,nd3b, ndf1a,ndf1b,ndf2a,
     & ndf2b,ndf3a,ndf3b

      real u(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
      real f(ndf1a:ndf1b,ndf2a:ndf2b,ndf3a:ndf3b,0:*)
      integer mask(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
      real rsxy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1,0:nd-1)
      real xy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1)

      integer ipar(0:*)
      real rpar(0:*)

!     --- local variables ----

      integer ierr,orderOfExtrapolation
      integer n1a,n1b,n2a,n2b,n3a,n3b
      integer i1,i2,i3

      ! Gaussian source variables: 
      real phi,phiE,phiH, phix,phiy,phiz, phixx,phiyy,phizz, phiEtt, 
     & lapF1,lapF2,lapF3, lapPhiH, lapPhi
      real ampE,ampH,cost,sint,tolGaussian,dtsqBy12 ,omega, pi
      real rampTime,tba,g,rt,rtt,gtt


      ! for macro:
       real x0,x1,x2,q,r2Min,realMin,betap,a0,ap,app,appp
       real rho,rhop,rhopp,rhoppp
       real r2,r2x,r2y,r2z,r2t,r2f,expr,rhox,rhoy,rhoz,rhot
       real r2xx,r2xy,r2yy,r2xz,r2yz,r2zz,r2xt,r2yt,r2zt,r2tt
       real rhoxx,rhoxy,rhoyy,rhoxz,rhoyz,rhozz
       real r2xxx,r2xxy,r2xyy,r2xxz,r2xzz,r2yyy,r2yyz,r2yzz,r2xyz,
     & r2zzz,r2xxt,r2yyt,r2zzt,r2ttt,r2xtt,r2ytt,r2ztt
       real rhoxxx,rhoxxy,rhoxyy,rhoxxz,rhoxzz,rhoyyy,rhoyyz,rhoyzz,
     & rhoxyz,rhozzz,rhoxxt,rhoyyt,rhozzt
       real rhoxtt,rhoytt,rhoztt,rhottt
       real deltaRhoX,deltaRhoY,deltaRhoZ,deltaRhot
       real deltaJ0t,deltaJ1t,deltaJ2t,j0ttt,j1ttt,j2ttt

      integer rectangular,curvilinear
      parameter(rectangular=0,curvilinear=1)

      ! forcing options
      ! forcingOptions -- these should match ForcingEnum in Maxwell.h 
      integer noForcing,magneticSinusoidalPointSource,gaussianSource,
     & twilightZoneForcing, gaussianChargeSource, 
     & userDefinedForcingOption
      integer noBoundaryForcing,planeWaveBoundaryForcing,
     & chirpedPlaneWaveBoundaryForcing
      parameter(noForcing                =0,
     & magneticSinusoidalPointSource =1,gaussianSource                
     & =2,twilightZoneForcing           =3,    gaussianChargeSource   
     &        =4,userDefinedForcingOption      =5 )
      ! boundary forcing options when solved directly for the scattered field:
      parameter( noBoundaryForcing              =0,   
     & planeWaveBoundaryForcing       =1,
     & chirpedPlaneWaveBoundaryForcing=2 )

      integer gridType,orderOfAccuracyInSpace,orderOfAccuracyInTime,
     & useForcing
      integer ex,ey,ez,hx,hy,hz,useWhereMask,grid,debug,forcingOption,
     & icount
      real dx(0:2),dr(0:2),t,ep,dt,c,csq,dtsq,eps,mu,kx,ky,kz,
     & slowStartInterval
      real xa(0:2)
      real amplitude,p,beta,xp0,xp1,xp2,vp0,vp1,vp2
      real tol, tolPulse

       integer method,nfdtd,bamx
       parameter( nfdtd=5, bamx=7 )



      ierr=0

      n1a                   =ipar(0)
      n1b                   =ipar(1)
      n2a                   =ipar(2)
      n2b                   =ipar(3)
      n3a                   =ipar(4)
      n3b                   =ipar(5)
      gridType              =ipar(6)
      orderOfAccuracyInSpace=ipar(7)
      orderOfAccuracyInTime =ipar(8)
      orderOfExtrapolation  =ipar(9)
      useForcing            =ipar(10)
      ex                    =ipar(11)
      ey                    =ipar(12)
      ez                    =ipar(13)
      hx                    =ipar(14)
      hy                    =ipar(15)
      hz                    =ipar(16)
      useWhereMask          =ipar(17)
      grid                  =ipar(18)
      debug                 =ipar(19)
      forcingOption         =ipar(20)
      method                =ipar(21)

      dx(0)                 =rpar(0)
      dx(1)                 =rpar(1)
      dx(2)                 =rpar(2)
      dr(0)                 =rpar(3)
      dr(1)                 =rpar(4)
      dr(2)                 =rpar(5)
      t                     =rpar(6)
      ep                    =rpar(7)
      dt                    =rpar(8)
      c                     =rpar(9)
      eps                   =rpar(10)
      mu                    =rpar(11)
      kx                    =rpar(12)  ! for plane wave forcing
      ky                    =rpar(13)
      kz                    =rpar(14)
      slowStartInterval     =rpar(15)
      xa(0)                 =rpar(16)  ! for rectangular grids
      xa(1)                 =rpar(17)
      xa(2)                 =rpar(18)

      amplitude             =rpar(19)  ! pulse parameters
      beta                  =rpar(20)  ! power
      p                     =rpar(21)  ! pulse parameters
      xp0                   =rpar(22)  ! initial position of the pulse
      xp1                   =rpar(23)
      xp2                   =rpar(24)
      vp0                   =rpar(25)  ! speed of the pulse
      vp1                   =rpar(26)
      vp2                   =rpar(27)
      omega                 =rpar(28)
      rampTime              =rpar(29)

      realMin=1.e-30  ! should be REAL_MIN ******************** fix this ************

      csq=c*c
      dtsq=dt*dt
      pi=atan2(1.,1.)*4.

      betap=beta**p
      q=p/2.
      r2Min=realMin*100.

      tol=1.e-12  ! ignore relative pulse values below this ******************** fix this ************
      ! exp( -(beta*r)^p ) < tol   (don't include amplitude to make relative)
      tolPulse = ( ( log(1./tol)**(1./p) )/beta )**2 ! tolerance for r2=r^2


      if( forcingOption.eq.gaussianChargeSource )then

       if( orderOfAccuracyInTime.eq.2 )then
        ! ***** 2nd order in time *****
        if( gridType.eq.curvilinear .and. nd.eq.2 )then

          !                               (DIM,ORDER,GRIDTYPE)
          ! f = grad(rho) - J_t 
          icount=0
          do i3=n3a,n3b
          do i2=n2a,n2b
          do i1=n1a,n1b
          if( mask(i1,i2,i3).gt.0 )then
              x0 = xy(i1,i2,i3,0)-xp0
              x1 = xy(i1,i2,i3,1)-xp1
              r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2
           if( mask(i1,i2,i3).gt.0 .and. (r2 .lt. tolPulse) )then
          ! if( .true. )then
             ! only evaluate if this is a valid point and pulse value is not too small
             ! f = f + c^2 [- grad(rho/eps) - mu J_t ]
             ! ** evaluate pulse and 1 derivative **
                ! radx == r*r_x = .5* (r^2)_x 
                r2x = 2.*(x0-vp0*(t))
                r2y = 2.*(x1-vp1*(t))
                  r2t = -2.*( ( x0-vp0*(t) )*vp0 + ( x1-vp1*(t) )*vp1 )
                a0=-betap*r2**q
                expr = exp(a0)  ! exp( -(beta*r)^p )
                rho = amplitude*expr
                r2f=max(r2Min,r2)  ! floor value of r2 to avoid division by zero
                ! rhop = [d(rho)/d(r2)]
                ap=a0*q/r2f  ! d(arg)/d(r2) , arg=-betap*r2**q
                rhop = rho*ap
                rhox = rhop*r2x
                rhoy = rhop*r2y
                rhot = rhop*r2t
             f(i1,i2,i3,ex)=f(i1,i2,i3,ex) + csq*( -rhox/eps - mu*vp0*
     & rhot )
             f(i1,i2,i3,ey)=f(i1,i2,i3,ey) + csq*( -rhoy/eps - mu*vp1*
     & rhot )
               f(i1,i2,i3,hz)= f(i1,i2,i3,hz) + csq*( vp1*rhox - vp0*
     & rhoy )
            else
             icount=icount+1 ! this point not evaluated
            end if
           end if
           end do
           end do
           end do
           ! write(*,'("charge source: total points=",i6," points evaluated=",i6," : r2>tolPulse=",e8.2)') ! (n1b-n1a+1)*(n2b-n2a+1)*(n3b-n3a+1),(n1b-n1a+1)*(n2b-n2a+1)*(n3b-n3a+1)-icount,tolPulse
           ! '

        else if( gridType.eq.rectangular .and. nd.eq.2 )then

          !                               (DIM,ORDER,GRIDTYPE)
          ! f = grad(rho) - J_t 
          icount=0
          do i3=n3a,n3b
          do i2=n2a,n2b
          do i1=n1a,n1b
          if( mask(i1,i2,i3).gt.0 )then
              x0 = xa(0)+i1*dx(0) -xp0
              x1 = xa(1)+i2*dx(1) -xp1
              r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2
           if( mask(i1,i2,i3).gt.0 .and. (r2 .lt. tolPulse) )then
          ! if( .true. )then
             ! only evaluate if this is a valid point and pulse value is not too small
             ! f = f + c^2 [- grad(rho/eps) - mu J_t ]
             ! ** evaluate pulse and 1 derivative **
                ! radx == r*r_x = .5* (r^2)_x 
                r2x = 2.*(x0-vp0*(t))
                r2y = 2.*(x1-vp1*(t))
                  r2t = -2.*( ( x0-vp0*(t) )*vp0 + ( x1-vp1*(t) )*vp1 )
                a0=-betap*r2**q
                expr = exp(a0)  ! exp( -(beta*r)^p )
                rho = amplitude*expr
                r2f=max(r2Min,r2)  ! floor value of r2 to avoid division by zero
                ! rhop = [d(rho)/d(r2)]
                ap=a0*q/r2f  ! d(arg)/d(r2) , arg=-betap*r2**q
                rhop = rho*ap
                rhox = rhop*r2x
                rhoy = rhop*r2y
                rhot = rhop*r2t
             f(i1,i2,i3,ex)=f(i1,i2,i3,ex) + csq*( -rhox/eps - mu*vp0*
     & rhot )
             f(i1,i2,i3,ey)=f(i1,i2,i3,ey) + csq*( -rhoy/eps - mu*vp1*
     & rhot )
               f(i1,i2,i3,hz)= f(i1,i2,i3,hz) + csq*( vp1*rhox - vp0*
     & rhoy )
            else
             icount=icount+1 ! this point not evaluated
            end if
           end if
           end do
           end do
           end do
           ! write(*,'("charge source: total points=",i6," points evaluated=",i6," : r2>tolPulse=",e8.2)') ! (n1b-n1a+1)*(n2b-n2a+1)*(n3b-n3a+1),(n1b-n1a+1)*(n2b-n2a+1)*(n3b-n3a+1)-icount,tolPulse
           ! '

        else if( gridType.eq.curvilinear .and. nd.eq.3 )then

          !                               (DIM,ORDER,GRIDTYPE)
          ! f = grad(rho) - J_t 
          icount=0
          do i3=n3a,n3b
          do i2=n2a,n2b
          do i1=n1a,n1b
          if( mask(i1,i2,i3).gt.0 )then
              x0 = xy(i1,i2,i3,0)-xp0
              x1 = xy(i1,i2,i3,1)-xp1
                x2 = xy(i1,i2,i3,2)-xp2
              r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*(t)
     &  )**2
           if( mask(i1,i2,i3).gt.0 .and. (r2 .lt. tolPulse) )then
          ! if( .true. )then
             ! only evaluate if this is a valid point and pulse value is not too small
             ! f = f + c^2 [- grad(rho/eps) - mu J_t ]
             ! ** evaluate pulse and 1 derivative **
                ! radx == r*r_x = .5* (r^2)_x 
                r2x = 2.*(x0-vp0*(t))
                r2y = 2.*(x1-vp1*(t))
                  r2z = 2.*(x2-vp2*(t))
                  r2t = -2.*( ( x0-vp0*(t) )*vp0 + ( x1-vp1*(t) )*vp1 +
     &  ( x2-vp2*(t) )*vp2 )
                a0=-betap*r2**q
                expr = exp(a0)  ! exp( -(beta*r)^p )
                rho = amplitude*expr
                r2f=max(r2Min,r2)  ! floor value of r2 to avoid division by zero
                ! rhop = [d(rho)/d(r2)]
                ap=a0*q/r2f  ! d(arg)/d(r2) , arg=-betap*r2**q
                rhop = rho*ap
                rhox = rhop*r2x
                rhoy = rhop*r2y
                rhot = rhop*r2t
                 rhoz = rhop*r2z
             f(i1,i2,i3,ex)=f(i1,i2,i3,ex) + csq*( -rhox/eps - mu*vp0*
     & rhot )
             f(i1,i2,i3,ey)=f(i1,i2,i3,ey) + csq*( -rhoy/eps - mu*vp1*
     & rhot )
               f(i1,i2,i3,ez)= f(i1,i2,i3,ez) + csq*( -rhoz/eps - mu*
     & vp2*rhot )
            else
             icount=icount+1 ! this point not evaluated
            end if
           end if
           end do
           end do
           end do
           ! write(*,'("charge source: total points=",i6," points evaluated=",i6," : r2>tolPulse=",e8.2)') ! (n1b-n1a+1)*(n2b-n2a+1)*(n3b-n3a+1),(n1b-n1a+1)*(n2b-n2a+1)*(n3b-n3a+1)-icount,tolPulse
           ! '

        else if( gridType.eq.rectangular .and. nd.eq.3 )then

          !                               (DIM,ORDER,GRIDTYPE)
          ! f = grad(rho) - J_t 
          icount=0
          do i3=n3a,n3b
          do i2=n2a,n2b
          do i1=n1a,n1b
          if( mask(i1,i2,i3).gt.0 )then
              x0 = xa(0)+i1*dx(0) -xp0
              x1 = xa(1)+i2*dx(1) -xp1
                x2 = xa(2)+i3*dx(2) -xp2
              r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*(t)
     &  )**2
           if( mask(i1,i2,i3).gt.0 .and. (r2 .lt. tolPulse) )then
          ! if( .true. )then
             ! only evaluate if this is a valid point and pulse value is not too small
             ! f = f + c^2 [- grad(rho/eps) - mu J_t ]
             ! ** evaluate pulse and 1 derivative **
                ! radx == r*r_x = .5* (r^2)_x 
                r2x = 2.*(x0-vp0*(t))
                r2y = 2.*(x1-vp1*(t))
                  r2z = 2.*(x2-vp2*(t))
                  r2t = -2.*( ( x0-vp0*(t) )*vp0 + ( x1-vp1*(t) )*vp1 +
     &  ( x2-vp2*(t) )*vp2 )
                a0=-betap*r2**q
                expr = exp(a0)  ! exp( -(beta*r)^p )
                rho = amplitude*expr
                r2f=max(r2Min,r2)  ! floor value of r2 to avoid division by zero
                ! rhop = [d(rho)/d(r2)]
                ap=a0*q/r2f  ! d(arg)/d(r2) , arg=-betap*r2**q
                rhop = rho*ap
                rhox = rhop*r2x
                rhoy = rhop*r2y
                rhot = rhop*r2t
                 rhoz = rhop*r2z
             f(i1,i2,i3,ex)=f(i1,i2,i3,ex) + csq*( -rhox/eps - mu*vp0*
     & rhot )
             f(i1,i2,i3,ey)=f(i1,i2,i3,ey) + csq*( -rhoy/eps - mu*vp1*
     & rhot )
               f(i1,i2,i3,ez)= f(i1,i2,i3,ez) + csq*( -rhoz/eps - mu*
     & vp2*rhot )
            else
             icount=icount+1 ! this point not evaluated
            end if
           end if
           end do
           end do
           end do
           ! write(*,'("charge source: total points=",i6," points evaluated=",i6," : r2>tolPulse=",e8.2)') ! (n1b-n1a+1)*(n2b-n2a+1)*(n3b-n3a+1),(n1b-n1a+1)*(n2b-n2a+1)*(n3b-n3a+1)-icount,tolPulse
           ! '

        else
          write(*,'("forcingOptMaxwell: case not implemented")')
          stop 7729
        end if

       else if( orderOfAccuracyInTime.eq.4 )then
        ! ***** 4th order in time *****
        if( gridType.eq.curvilinear .and. nd.eq.2 )then
          !                               (DIM,ORDER,GRIDTYPE)
          ! f = grad(rho) - J_t 
          icount=0
          do i3=n3a,n3b
          do i2=n2a,n2b
          do i1=n1a,n1b
          if( mask(i1,i2,i3).gt.0 )then
              x0 = xy(i1,i2,i3,0)-xp0
              x1 = xy(i1,i2,i3,1)-xp1
              r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2
           if( mask(i1,i2,i3).gt.0 .and. (r2 .lt. tolPulse) )then
          ! if( .true. )then
             ! only evaluate if this is a valid point and pulse value is not too small
             ! for 4th order in time we need 
             !    f_tt and Delta(f)
             ! f = f + c^2 [- grad(rho/eps) - mu J_t ]
             !       + c^2*dt^2/12*{ c^2*[ -Delta*grad(rho/eps) - mu*Delta(J_t) ] - grad(rho_tt/eps) - mu*J_ttt } 
             ! ** evaluate pulse and 3 derivatives **
                 ! radx == r*r_x = .5* (r^2)_x 
                 r2x = 2.*(x0-vp0*(t))
                 r2y = 2.*(x1-vp1*(t))
                   r2t = -2.*( ( x0-vp0*(t) )*vp0 + ( x1-vp1*(t) )*vp1 
     & )
                 a0=-betap*r2**q
                 expr = exp(a0)  ! exp( -(beta*r)^p )
                 rho = amplitude*expr
                 r2f=max(r2Min,r2)  ! floor value of r2 to avoid division by zero
                 ! rhop = [d(rho)/d(r2)]
                 ap=a0*q/r2f  ! d(arg)/d(r2) , arg=-betap*r2**q
                 rhop = rho*ap
                 rhox = rhop*r2x
                 rhoy = rhop*r2y
                 rhot = rhop*r2t
                   app=ap*(q-1.)/r2f
                   rhopp = rhop*ap+rho*app
                   r2xx=2.
                   r2xy=0.
                   r2yy=2.
                   r2xt=-2.*vp0
                   r2yt=-2.*vp1
                     r2tt=2.*(vp0**2+vp1**2)
                   rhoxx = rhopp*r2x*r2x + rhop*r2xx
                   rhoxy = rhopp*r2x*r2y + rhop*r2xy
                   rhoyy = rhopp*r2y*r2y + rhop*r2yy
                   appp=app*(q-2.)/r2f
                   rhoppp = rhopp*ap +2.*rhop*app+rho*appp
                   r2xxx=0.
                   r2xxy=0.
                   r2xyy=0.
                   r2xxt=0.
                   r2yyt=0.
                   r2ttt=0.
                   rhoxxx = rhoppp*r2x*r2x*r2x + rhopp*( r2x*r2xx + 
     & r2x*r2xx + r2x*r2xx ) + rhop*r2xxx
                   rhoxxy = rhoppp*r2x*r2x*r2y + rhopp*( r2y*r2xx + 
     & r2x*r2xy + r2x*r2xy ) + rhop*r2xxy
                   rhoxyy = rhoppp*r2x*r2y*r2y + rhopp*( r2y*r2xy + 
     & r2y*r2xy + r2x*r2yy ) + rhop*r2xyy
                   rhoyyy = rhoppp*r2y*r2y*r2y + rhopp*( r2y*r2yy + 
     & r2y*r2yy + r2y*r2yy ) + rhop*r2yyy
                   rhoxxt = rhoppp*r2x*r2x*r2t + rhopp*( r2t*r2xx + 
     & r2x*r2xt + r2x*r2xt ) + rhop*r2xxt
                   rhoyyt = rhoppp*r2y*r2y*r2t + rhopp*( r2t*r2yy + 
     & r2y*r2yt + r2y*r2yt ) + rhop*r2yyt
                   rhottt = rhoppp*r2t*r2t*r2t + rhopp*( r2t*r2tt + 
     & r2t*r2tt + r2t*r2tt ) + rhop*r2ttt
                   rhoxtt = rhoppp*r2x*r2t*r2t + rhopp*( r2t*r2xt + 
     & r2t*r2xt + r2x*r2tt ) + rhop*r2xtt
                   rhoytt = rhoppp*r2y*r2t*r2t + rhopp*( r2t*r2yt + 
     & r2t*r2yt + r2y*r2tt ) + rhop*r2ytt
               deltaRhoX = rhoxxx+rhoxyy
               deltaRhoY = rhoxxy+rhoyyy
               deltaRhot = rhoxxt+rhoyyt
              deltaJ0t = vp0*deltaRhot
              deltaJ1t = vp1*deltaRhot
              j0ttt = vp0*rhottt
              j1ttt = vp1*rhottt
               deltaJ2t = vp2*deltaRhot
               j2ttt = vp2*rhottt
              f(i1,i2,i3,ex)=f(i1,i2,i3,ex) + csq*( -rhox/eps - mu*vp0*
     & rhot ) + csq*dtsq/12.*( -csq*( deltaRhoX/eps + mu*deltaJ0t ) - 
     & rhoxtt/eps -mu*j0ttt )
              f(i1,i2,i3,ey)=f(i1,i2,i3,ey) + csq*( -rhoy/eps - mu*vp1*
     & rhot ) + csq*dtsq/12.*( -csq*( deltaRhoY/eps + mu*deltaJ1t ) - 
     & rhoytt/eps -mu*j1ttt )
                f(i1,i2,i3,hz)= f(i1,i2,i3,hz) + csq*( vp1*rhox - vp0*
     & rhoy ) + csq*dtsq/12.*( vp1*deltaRhox - vp0*deltaRhoy +vp1*
     & rhoxtt - vp0*rhoytt )
            else
             icount=icount+1 ! this point not evaluated
            end if
           end if
           end do
           end do
           end do
           ! write(*,'("charge source: total points=",i6," points evaluated=",i6," : r2>tolPulse=",e8.2)') ! (n1b-n1a+1)*(n2b-n2a+1)*(n3b-n3a+1),(n1b-n1a+1)*(n2b-n2a+1)*(n3b-n3a+1)-icount,tolPulse
           ! '

        else if( gridType.eq.rectangular .and. nd.eq.2 )then

          !                               (DIM,ORDER,GRIDTYPE)
          ! f = grad(rho) - J_t 
          icount=0
          do i3=n3a,n3b
          do i2=n2a,n2b
          do i1=n1a,n1b
          if( mask(i1,i2,i3).gt.0 )then
              x0 = xa(0)+i1*dx(0) -xp0
              x1 = xa(1)+i2*dx(1) -xp1
              r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2
           if( mask(i1,i2,i3).gt.0 .and. (r2 .lt. tolPulse) )then
          ! if( .true. )then
             ! only evaluate if this is a valid point and pulse value is not too small
             ! for 4th order in time we need 
             !    f_tt and Delta(f)
             ! f = f + c^2 [- grad(rho/eps) - mu J_t ]
             !       + c^2*dt^2/12*{ c^2*[ -Delta*grad(rho/eps) - mu*Delta(J_t) ] - grad(rho_tt/eps) - mu*J_ttt } 
             ! ** evaluate pulse and 3 derivatives **
                 ! radx == r*r_x = .5* (r^2)_x 
                 r2x = 2.*(x0-vp0*(t))
                 r2y = 2.*(x1-vp1*(t))
                   r2t = -2.*( ( x0-vp0*(t) )*vp0 + ( x1-vp1*(t) )*vp1 
     & )
                 a0=-betap*r2**q
                 expr = exp(a0)  ! exp( -(beta*r)^p )
                 rho = amplitude*expr
                 r2f=max(r2Min,r2)  ! floor value of r2 to avoid division by zero
                 ! rhop = [d(rho)/d(r2)]
                 ap=a0*q/r2f  ! d(arg)/d(r2) , arg=-betap*r2**q
                 rhop = rho*ap
                 rhox = rhop*r2x
                 rhoy = rhop*r2y
                 rhot = rhop*r2t
                   app=ap*(q-1.)/r2f
                   rhopp = rhop*ap+rho*app
                   r2xx=2.
                   r2xy=0.
                   r2yy=2.
                   r2xt=-2.*vp0
                   r2yt=-2.*vp1
                     r2tt=2.*(vp0**2+vp1**2)
                   rhoxx = rhopp*r2x*r2x + rhop*r2xx
                   rhoxy = rhopp*r2x*r2y + rhop*r2xy
                   rhoyy = rhopp*r2y*r2y + rhop*r2yy
                   appp=app*(q-2.)/r2f
                   rhoppp = rhopp*ap +2.*rhop*app+rho*appp
                   r2xxx=0.
                   r2xxy=0.
                   r2xyy=0.
                   r2xxt=0.
                   r2yyt=0.
                   r2ttt=0.
                   rhoxxx = rhoppp*r2x*r2x*r2x + rhopp*( r2x*r2xx + 
     & r2x*r2xx + r2x*r2xx ) + rhop*r2xxx
                   rhoxxy = rhoppp*r2x*r2x*r2y + rhopp*( r2y*r2xx + 
     & r2x*r2xy + r2x*r2xy ) + rhop*r2xxy
                   rhoxyy = rhoppp*r2x*r2y*r2y + rhopp*( r2y*r2xy + 
     & r2y*r2xy + r2x*r2yy ) + rhop*r2xyy
                   rhoyyy = rhoppp*r2y*r2y*r2y + rhopp*( r2y*r2yy + 
     & r2y*r2yy + r2y*r2yy ) + rhop*r2yyy
                   rhoxxt = rhoppp*r2x*r2x*r2t + rhopp*( r2t*r2xx + 
     & r2x*r2xt + r2x*r2xt ) + rhop*r2xxt
                   rhoyyt = rhoppp*r2y*r2y*r2t + rhopp*( r2t*r2yy + 
     & r2y*r2yt + r2y*r2yt ) + rhop*r2yyt
                   rhottt = rhoppp*r2t*r2t*r2t + rhopp*( r2t*r2tt + 
     & r2t*r2tt + r2t*r2tt ) + rhop*r2ttt
                   rhoxtt = rhoppp*r2x*r2t*r2t + rhopp*( r2t*r2xt + 
     & r2t*r2xt + r2x*r2tt ) + rhop*r2xtt
                   rhoytt = rhoppp*r2y*r2t*r2t + rhopp*( r2t*r2yt + 
     & r2t*r2yt + r2y*r2tt ) + rhop*r2ytt
               deltaRhoX = rhoxxx+rhoxyy
               deltaRhoY = rhoxxy+rhoyyy
               deltaRhot = rhoxxt+rhoyyt
              deltaJ0t = vp0*deltaRhot
              deltaJ1t = vp1*deltaRhot
              j0ttt = vp0*rhottt
              j1ttt = vp1*rhottt
               deltaJ2t = vp2*deltaRhot
               j2ttt = vp2*rhottt
              f(i1,i2,i3,ex)=f(i1,i2,i3,ex) + csq*( -rhox/eps - mu*vp0*
     & rhot ) + csq*dtsq/12.*( -csq*( deltaRhoX/eps + mu*deltaJ0t ) - 
     & rhoxtt/eps -mu*j0ttt )
              f(i1,i2,i3,ey)=f(i1,i2,i3,ey) + csq*( -rhoy/eps - mu*vp1*
     & rhot ) + csq*dtsq/12.*( -csq*( deltaRhoY/eps + mu*deltaJ1t ) - 
     & rhoytt/eps -mu*j1ttt )
                f(i1,i2,i3,hz)= f(i1,i2,i3,hz) + csq*( vp1*rhox - vp0*
     & rhoy ) + csq*dtsq/12.*( vp1*deltaRhox - vp0*deltaRhoy +vp1*
     & rhoxtt - vp0*rhoytt )
            else
             icount=icount+1 ! this point not evaluated
            end if
           end if
           end do
           end do
           end do
           ! write(*,'("charge source: total points=",i6," points evaluated=",i6," : r2>tolPulse=",e8.2)') ! (n1b-n1a+1)*(n2b-n2a+1)*(n3b-n3a+1),(n1b-n1a+1)*(n2b-n2a+1)*(n3b-n3a+1)-icount,tolPulse
           ! '

        else if( gridType.eq.curvilinear .and. nd.eq.3 )then

          !                               (DIM,ORDER,GRIDTYPE)
          ! f = grad(rho) - J_t 
          icount=0
          do i3=n3a,n3b
          do i2=n2a,n2b
          do i1=n1a,n1b
          if( mask(i1,i2,i3).gt.0 )then
              x0 = xy(i1,i2,i3,0)-xp0
              x1 = xy(i1,i2,i3,1)-xp1
                x2 = xy(i1,i2,i3,2)-xp2
              r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*(t)
     &  )**2
           if( mask(i1,i2,i3).gt.0 .and. (r2 .lt. tolPulse) )then
          ! if( .true. )then
             ! only evaluate if this is a valid point and pulse value is not too small
             ! for 4th order in time we need 
             !    f_tt and Delta(f)
             ! f = f + c^2 [- grad(rho/eps) - mu J_t ]
             !       + c^2*dt^2/12*{ c^2*[ -Delta*grad(rho/eps) - mu*Delta(J_t) ] - grad(rho_tt/eps) - mu*J_ttt } 
             ! ** evaluate pulse and 3 derivatives **
                 ! radx == r*r_x = .5* (r^2)_x 
                 r2x = 2.*(x0-vp0*(t))
                 r2y = 2.*(x1-vp1*(t))
                   r2z = 2.*(x2-vp2*(t))
                   r2t = -2.*( ( x0-vp0*(t) )*vp0 + ( x1-vp1*(t) )*vp1 
     & + ( x2-vp2*(t) )*vp2 )
                 a0=-betap*r2**q
                 expr = exp(a0)  ! exp( -(beta*r)^p )
                 rho = amplitude*expr
                 r2f=max(r2Min,r2)  ! floor value of r2 to avoid division by zero
                 ! rhop = [d(rho)/d(r2)]
                 ap=a0*q/r2f  ! d(arg)/d(r2) , arg=-betap*r2**q
                 rhop = rho*ap
                 rhox = rhop*r2x
                 rhoy = rhop*r2y
                 rhot = rhop*r2t
                  rhoz = rhop*r2z
                   app=ap*(q-1.)/r2f
                   rhopp = rhop*ap+rho*app
                   r2xx=2.
                   r2xy=0.
                   r2yy=2.
                   r2xt=-2.*vp0
                   r2yt=-2.*vp1
                     r2xz=0.
                     r2yz=0.
                     r2zz=2.
                     r2zt=-2.*vp2
                     r2tt=2.*(vp0**2+vp1**2+vp2**2)
                   rhoxx = rhopp*r2x*r2x + rhop*r2xx
                   rhoxy = rhopp*r2x*r2y + rhop*r2xy
                   rhoyy = rhopp*r2y*r2y + rhop*r2yy
                     rhoxz = rhopp*r2x*r2z + rhop*r2xz
                     rhoyz = rhopp*r2y*r2z + rhop*r2yz
                     rhozz = rhopp*r2z*r2z + rhop*r2zz
                   appp=app*(q-2.)/r2f
                   rhoppp = rhopp*ap +2.*rhop*app+rho*appp
                   r2xxx=0.
                   r2xxy=0.
                   r2xyy=0.
                   r2xxt=0.
                   r2yyt=0.
                   r2ttt=0.
                   rhoxxx = rhoppp*r2x*r2x*r2x + rhopp*( r2x*r2xx + 
     & r2x*r2xx + r2x*r2xx ) + rhop*r2xxx
                   rhoxxy = rhoppp*r2x*r2x*r2y + rhopp*( r2y*r2xx + 
     & r2x*r2xy + r2x*r2xy ) + rhop*r2xxy
                   rhoxyy = rhoppp*r2x*r2y*r2y + rhopp*( r2y*r2xy + 
     & r2y*r2xy + r2x*r2yy ) + rhop*r2xyy
                   rhoyyy = rhoppp*r2y*r2y*r2y + rhopp*( r2y*r2yy + 
     & r2y*r2yy + r2y*r2yy ) + rhop*r2yyy
                   rhoxxt = rhoppp*r2x*r2x*r2t + rhopp*( r2t*r2xx + 
     & r2x*r2xt + r2x*r2xt ) + rhop*r2xxt
                   rhoyyt = rhoppp*r2y*r2y*r2t + rhopp*( r2t*r2yy + 
     & r2y*r2yt + r2y*r2yt ) + rhop*r2yyt
                   rhottt = rhoppp*r2t*r2t*r2t + rhopp*( r2t*r2tt + 
     & r2t*r2tt + r2t*r2tt ) + rhop*r2ttt
                   rhoxtt = rhoppp*r2x*r2t*r2t + rhopp*( r2t*r2xt + 
     & r2t*r2xt + r2x*r2tt ) + rhop*r2xtt
                   rhoytt = rhoppp*r2y*r2t*r2t + rhopp*( r2t*r2yt + 
     & r2t*r2yt + r2y*r2tt ) + rhop*r2ytt
                     r2xxz=0.
                     r2xzz=0.
                     r2yyz=0.
                     r2yzz=0.
                     r2zzz=0.
                     r2zzt=0.
                     rhoxxz = rhoppp*r2x*r2x*r2z + rhopp*( r2z*r2xx + 
     & r2x*r2xz + r2x*r2xz ) + rhop*r2xxz
                     rhoxzz = rhoppp*r2x*r2z*r2z + rhopp*( r2z*r2xz + 
     & r2z*r2xz + r2x*r2zz ) + rhop*r2xzz
                     rhoyzz = rhoppp*r2y*r2z*r2z + rhopp*( r2z*r2yz + 
     & r2z*r2yz + r2y*r2zz ) + rhop*r2yzz
                     rhoyyz = rhoppp*r2y*r2y*r2z + rhopp*( r2z*r2yy + 
     & r2y*r2yz + r2y*r2yz ) + rhop*r2yyz
                     rhozzz = rhoppp*r2z*r2z*r2z + rhopp*( r2z*r2zz + 
     & r2z*r2zz + r2z*r2zz ) + rhop*r2zzz
                     rhozzt = rhoppp*r2z*r2z*r2t + rhopp*( r2t*r2zz + 
     & r2z*r2zt + r2z*r2zt ) + rhop*r2zzt
                     rhoztt = rhoppp*r2z*r2t*r2t + rhopp*( r2t*r2zt + 
     & r2t*r2zt + r2z*r2tt ) + rhop*r2ztt
               deltaRhoX = rhoxxx+rhoxyy+rhoxzz
               deltaRhoY = rhoxxy+rhoyyy+rhoyzz
               deltaRhoZ = rhoxxz+rhoyyz+rhozzz
               deltaRhot = rhoxxt+rhoyyt+rhozzt
              deltaJ0t = vp0*deltaRhot
              deltaJ1t = vp1*deltaRhot
              j0ttt = vp0*rhottt
              j1ttt = vp1*rhottt
              f(i1,i2,i3,ex)=f(i1,i2,i3,ex) + csq*( -rhox/eps - mu*vp0*
     & rhot ) + csq*dtsq/12.*( -csq*( deltaRhoX/eps + mu*deltaJ0t ) - 
     & rhoxtt/eps -mu*j0ttt )
              f(i1,i2,i3,ey)=f(i1,i2,i3,ey) + csq*( -rhoy/eps - mu*vp1*
     & rhot ) + csq*dtsq/12.*( -csq*( deltaRhoY/eps + mu*deltaJ1t ) - 
     & rhoytt/eps -mu*j1ttt )
               f(i1,i2,i3,ez)=f(i1,i2,i3,ez) + csq*( -rhoz/eps - mu*
     & vp2*rhot ) + csq*dtsq/12.*( -csq*( deltaRhoZ/eps + mu*deltaJ2t 
     & ) - rhoztt/eps -mu*j2ttt )
            else
             icount=icount+1 ! this point not evaluated
            end if
           end if
           end do
           end do
           end do
           ! write(*,'("charge source: total points=",i6," points evaluated=",i6," : r2>tolPulse=",e8.2)') ! (n1b-n1a+1)*(n2b-n2a+1)*(n3b-n3a+1),(n1b-n1a+1)*(n2b-n2a+1)*(n3b-n3a+1)-icount,tolPulse
           ! '

        else if( gridType.eq.rectangular .and. nd.eq.3 )then

          !                               (DIM,ORDER,GRIDTYPE)
          ! f = grad(rho) - J_t 
          icount=0
          do i3=n3a,n3b
          do i2=n2a,n2b
          do i1=n1a,n1b
          if( mask(i1,i2,i3).gt.0 )then
              x0 = xa(0)+i1*dx(0) -xp0
              x1 = xa(1)+i2*dx(1) -xp1
                x2 = xa(2)+i3*dx(2) -xp2
              r2 = ( x0-vp0*(t) )**2 + ( x1-vp1*(t) )**2 + ( x2-vp2*(t)
     &  )**2
           if( mask(i1,i2,i3).gt.0 .and. (r2 .lt. tolPulse) )then
          ! if( .true. )then
             ! only evaluate if this is a valid point and pulse value is not too small
             ! for 4th order in time we need 
             !    f_tt and Delta(f)
             ! f = f + c^2 [- grad(rho/eps) - mu J_t ]
             !       + c^2*dt^2/12*{ c^2*[ -Delta*grad(rho/eps) - mu*Delta(J_t) ] - grad(rho_tt/eps) - mu*J_ttt } 
             ! ** evaluate pulse and 3 derivatives **
                 ! radx == r*r_x = .5* (r^2)_x 
                 r2x = 2.*(x0-vp0*(t))
                 r2y = 2.*(x1-vp1*(t))
                   r2z = 2.*(x2-vp2*(t))
                   r2t = -2.*( ( x0-vp0*(t) )*vp0 + ( x1-vp1*(t) )*vp1 
     & + ( x2-vp2*(t) )*vp2 )
                 a0=-betap*r2**q
                 expr = exp(a0)  ! exp( -(beta*r)^p )
                 rho = amplitude*expr
                 r2f=max(r2Min,r2)  ! floor value of r2 to avoid division by zero
                 ! rhop = [d(rho)/d(r2)]
                 ap=a0*q/r2f  ! d(arg)/d(r2) , arg=-betap*r2**q
                 rhop = rho*ap
                 rhox = rhop*r2x
                 rhoy = rhop*r2y
                 rhot = rhop*r2t
                  rhoz = rhop*r2z
                   app=ap*(q-1.)/r2f
                   rhopp = rhop*ap+rho*app
                   r2xx=2.
                   r2xy=0.
                   r2yy=2.
                   r2xt=-2.*vp0
                   r2yt=-2.*vp1
                     r2xz=0.
                     r2yz=0.
                     r2zz=2.
                     r2zt=-2.*vp2
                     r2tt=2.*(vp0**2+vp1**2+vp2**2)
                   rhoxx = rhopp*r2x*r2x + rhop*r2xx
                   rhoxy = rhopp*r2x*r2y + rhop*r2xy
                   rhoyy = rhopp*r2y*r2y + rhop*r2yy
                     rhoxz = rhopp*r2x*r2z + rhop*r2xz
                     rhoyz = rhopp*r2y*r2z + rhop*r2yz
                     rhozz = rhopp*r2z*r2z + rhop*r2zz
                   appp=app*(q-2.)/r2f
                   rhoppp = rhopp*ap +2.*rhop*app+rho*appp
                   r2xxx=0.
                   r2xxy=0.
                   r2xyy=0.
                   r2xxt=0.
                   r2yyt=0.
                   r2ttt=0.
                   rhoxxx = rhoppp*r2x*r2x*r2x + rhopp*( r2x*r2xx + 
     & r2x*r2xx + r2x*r2xx ) + rhop*r2xxx
                   rhoxxy = rhoppp*r2x*r2x*r2y + rhopp*( r2y*r2xx + 
     & r2x*r2xy + r2x*r2xy ) + rhop*r2xxy
                   rhoxyy = rhoppp*r2x*r2y*r2y + rhopp*( r2y*r2xy + 
     & r2y*r2xy + r2x*r2yy ) + rhop*r2xyy
                   rhoyyy = rhoppp*r2y*r2y*r2y + rhopp*( r2y*r2yy + 
     & r2y*r2yy + r2y*r2yy ) + rhop*r2yyy
                   rhoxxt = rhoppp*r2x*r2x*r2t + rhopp*( r2t*r2xx + 
     & r2x*r2xt + r2x*r2xt ) + rhop*r2xxt
                   rhoyyt = rhoppp*r2y*r2y*r2t + rhopp*( r2t*r2yy + 
     & r2y*r2yt + r2y*r2yt ) + rhop*r2yyt
                   rhottt = rhoppp*r2t*r2t*r2t + rhopp*( r2t*r2tt + 
     & r2t*r2tt + r2t*r2tt ) + rhop*r2ttt
                   rhoxtt = rhoppp*r2x*r2t*r2t + rhopp*( r2t*r2xt + 
     & r2t*r2xt + r2x*r2tt ) + rhop*r2xtt
                   rhoytt = rhoppp*r2y*r2t*r2t + rhopp*( r2t*r2yt + 
     & r2t*r2yt + r2y*r2tt ) + rhop*r2ytt
                     r2xxz=0.
                     r2xzz=0.
                     r2yyz=0.
                     r2yzz=0.
                     r2zzz=0.
                     r2zzt=0.
                     rhoxxz = rhoppp*r2x*r2x*r2z + rhopp*( r2z*r2xx + 
     & r2x*r2xz + r2x*r2xz ) + rhop*r2xxz
                     rhoxzz = rhoppp*r2x*r2z*r2z + rhopp*( r2z*r2xz + 
     & r2z*r2xz + r2x*r2zz ) + rhop*r2xzz
                     rhoyzz = rhoppp*r2y*r2z*r2z + rhopp*( r2z*r2yz + 
     & r2z*r2yz + r2y*r2zz ) + rhop*r2yzz
                     rhoyyz = rhoppp*r2y*r2y*r2z + rhopp*( r2z*r2yy + 
     & r2y*r2yz + r2y*r2yz ) + rhop*r2yyz
                     rhozzz = rhoppp*r2z*r2z*r2z + rhopp*( r2z*r2zz + 
     & r2z*r2zz + r2z*r2zz ) + rhop*r2zzz
                     rhozzt = rhoppp*r2z*r2z*r2t + rhopp*( r2t*r2zz + 
     & r2z*r2zt + r2z*r2zt ) + rhop*r2zzt
                     rhoztt = rhoppp*r2z*r2t*r2t + rhopp*( r2t*r2zt + 
     & r2t*r2zt + r2z*r2tt ) + rhop*r2ztt
               deltaRhoX = rhoxxx+rhoxyy+rhoxzz
               deltaRhoY = rhoxxy+rhoyyy+rhoyzz
               deltaRhoZ = rhoxxz+rhoyyz+rhozzz
               deltaRhot = rhoxxt+rhoyyt+rhozzt
              deltaJ0t = vp0*deltaRhot
              deltaJ1t = vp1*deltaRhot
              j0ttt = vp0*rhottt
              j1ttt = vp1*rhottt
              f(i1,i2,i3,ex)=f(i1,i2,i3,ex) + csq*( -rhox/eps - mu*vp0*
     & rhot ) + csq*dtsq/12.*( -csq*( deltaRhoX/eps + mu*deltaJ0t ) - 
     & rhoxtt/eps -mu*j0ttt )
              f(i1,i2,i3,ey)=f(i1,i2,i3,ey) + csq*( -rhoy/eps - mu*vp1*
     & rhot ) + csq*dtsq/12.*( -csq*( deltaRhoY/eps + mu*deltaJ1t ) - 
     & rhoytt/eps -mu*j1ttt )
               f(i1,i2,i3,ez)=f(i1,i2,i3,ez) + csq*( -rhoz/eps - mu*
     & vp2*rhot ) + csq*dtsq/12.*( -csq*( deltaRhoZ/eps + mu*deltaJ2t 
     & ) - rhoztt/eps -mu*j2ttt )
            else
             icount=icount+1 ! this point not evaluated
            end if
           end if
           end do
           end do
           end do
           ! write(*,'("charge source: total points=",i6," points evaluated=",i6," : r2>tolPulse=",e8.2)') ! (n1b-n1a+1)*(n2b-n2a+1)*(n3b-n3a+1),(n1b-n1a+1)*(n2b-n2a+1)*(n3b-n3a+1)-icount,tolPulse
           ! '

        else
          write(*,'("forcingOptMaxwell: case not implemented")')
          stop 7729
        end if

       else
        write(*,'("forcingOptMaxwell: unknown orderOfAccuracyInTime",
     & i6)') orderOfAccuracyInTime
        stop 5219
       end if


      else if( forcingOption.eq.gaussianSource )then

       ! --- Gaussian Source ------   July 4, 2020

       ! amp = amplitude*beta**2   ! scale by beta^2 to make solution O(1) 
       if( method .eq. nfdtd )then
         ampE = amplitude*(2.*pi*omega)**2*sqrt(2.*beta)   ! scale to make solution O(1)
         ampH = amplitude*(2.*pi*omega)**2                 ! scale to make solution O(1)
       else if( method.eq.bamx )then
         ampE = amplitude*(2.*pi*omega)*sqrt(2.*beta)   ! scale to make solution O(1)
         ampH = amplitude*(2.*pi*omega)                 ! scale to make solution O(1)
       else
          write(*,'("forcingOpt: Unknown method=",i6)') method
          stop 11223
       end if
       ! exp( -beta*r2 ) < tol
       !  beta*r2 < log(1/tol)          
       tolGaussian = log(1./tol)/beta ! tolerance for r2=r^2

       dtsqBy12 = (dt**2)/12.

       sint = sin(2.*pi*omega*t )
       cost = cos(2.*pi*omega*t )

       ! Forcing time dependence is g(t)
       tba = rampTime
       if( t.le.tba )then
         ! slow start time function        
         g = (t)*(t)*(-(t)/3.+.5*tba)*6./(tba*tba*tba)*sint
         rt = (t)*(-(t)+tba)*6./(tba*tba*tba)
         rtt = (-2.*(t)+tba)*6./(tba*tba*tba)
         gtt = rtt*sint + 2.*rt*cost*(2.*pi*omega) - (2.*pi*omega)**2 *
     &  g
       else
         g   = sint
         gtt = -(2.*pi*omega)**2 * sint
       end if

       if( t.le.5*dt )then
         write(*,'(">> forcingOpt:GaussianSource: t=",e10.2," beta=",
     & e10.2," omega=",e10.2)') t,beta,omega
         write(*,'(">> amplitude",e10.2," ampE=",e10.2," ampH=",e10.2)
     & ') amplitude,ampE,ampH
         write(*,'(">> rampTime=",e10.2," g=",e10.2," gtt=",e10.2)') 
     & rampTime,g,gtt
         write(*,'(">> method=",i3," (5=nfdtd, 7=bamx)")') method
       end if

       if( orderOfAccuracyInTime.eq.2 .or. method.eq.bamx )then
        ! ***** NFDTD: 2nd order in time , OR BAMX (RK)*****
        if( gridType.eq.curvilinear .and. nd.eq.2 )then
          icount=0
          do i3=n3a,n3b
          do i2=n2a,n2b
          do i1=n1a,n1b
          if( mask(i1,i2,i3).gt.0 )then
              x0 = xy(i1,i2,i3,0)-xp0
              x1 = xy(i1,i2,i3,1)-xp1
              r2 = ( x0 )**2 + ( x1 )**2
           ! only evaluate if this is a valid point and pulse value is not too small
           if( mask(i1,i2,i3).gt.0 .and. (r2 .lt. tolGaussian) )then
           ! if( .true. )then
            phi = ampE*exp( -beta*( r2 ) )
            phiE = phi*g
               f(i1,i2,i3,ex)=f(i1,i2,i3,ex) -x1*phiE
               f(i1,i2,i3,ey)=f(i1,i2,i3,ey) +x0*phiE
               phiH = phi*g
               f(i1,i2,i3,hz)= f(i1,i2,i3,hz) + (ampH/ampE)*phiH
            else
             icount=icount+1 ! this point not evaluated
            end if
           end if
           end do
           end do
           end do
           if( t.le.5*dt )then
           write(*,'("Gaussian source: total points=",i6," points 
     & evaluated=",i6," : r2>tolGaussian=",e8.2)') (n1b-n1a+1)*(n2b-
     & n2a+1)*(n3b-n3a+1),(n1b-n1a+1)*(n2b-n2a+1)*(n3b-n3a+1)-icount,
     & tolGaussian
           end if

        else if( gridType.eq.rectangular .and. nd.eq.2 )then
          icount=0
          do i3=n3a,n3b
          do i2=n2a,n2b
          do i1=n1a,n1b
          if( mask(i1,i2,i3).gt.0 )then
              x0 = xa(0)+i1*dx(0) -xp0
              x1 = xa(1)+i2*dx(1) -xp1
              r2 = ( x0 )**2 + ( x1 )**2
           ! only evaluate if this is a valid point and pulse value is not too small
           if( mask(i1,i2,i3).gt.0 .and. (r2 .lt. tolGaussian) )then
           ! if( .true. )then
            phi = ampE*exp( -beta*( r2 ) )
            phiE = phi*g
               f(i1,i2,i3,ex)=f(i1,i2,i3,ex) -x1*phiE
               f(i1,i2,i3,ey)=f(i1,i2,i3,ey) +x0*phiE
               phiH = phi*g
               f(i1,i2,i3,hz)= f(i1,i2,i3,hz) + (ampH/ampE)*phiH
            else
             icount=icount+1 ! this point not evaluated
            end if
           end if
           end do
           end do
           end do
           if( t.le.5*dt )then
           write(*,'("Gaussian source: total points=",i6," points 
     & evaluated=",i6," : r2>tolGaussian=",e8.2)') (n1b-n1a+1)*(n2b-
     & n2a+1)*(n3b-n3a+1),(n1b-n1a+1)*(n2b-n2a+1)*(n3b-n3a+1)-icount,
     & tolGaussian
           end if

        else if( gridType.eq.curvilinear .and. nd.eq.3 )then
          icount=0
          do i3=n3a,n3b
          do i2=n2a,n2b
          do i1=n1a,n1b
          if( mask(i1,i2,i3).gt.0 )then
              x0 = xy(i1,i2,i3,0)-xp0
              x1 = xy(i1,i2,i3,1)-xp1
                x2 = xy(i1,i2,i3,2)-xp2
              r2 = ( x0 )**2 + ( x1 )**2 + ( x2 )**2
           ! only evaluate if this is a valid point and pulse value is not too small
           if( mask(i1,i2,i3).gt.0 .and. (r2 .lt. tolGaussian) )then
           ! if( .true. )then
            phi = ampE*exp( -beta*( r2 ) )
            phiE = phi*g
               f(i1,i2,i3,ex)=f(i1,i2,i3,ex) + (x2-x1)*phiE
               f(i1,i2,i3,ey)=f(i1,i2,i3,ey) + (x0-x2)*phiE
               f(i1,i2,i3,ez)=f(i1,i2,i3,ez) + (x1-x0)*phiE
            else
             icount=icount+1 ! this point not evaluated
            end if
           end if
           end do
           end do
           end do
           if( t.le.5*dt )then
           write(*,'("Gaussian source: total points=",i6," points 
     & evaluated=",i6," : r2>tolGaussian=",e8.2)') (n1b-n1a+1)*(n2b-
     & n2a+1)*(n3b-n3a+1),(n1b-n1a+1)*(n2b-n2a+1)*(n3b-n3a+1)-icount,
     & tolGaussian
           end if

        else if( gridType.eq.rectangular .and. nd.eq.3 )then
          icount=0
          do i3=n3a,n3b
          do i2=n2a,n2b
          do i1=n1a,n1b
          if( mask(i1,i2,i3).gt.0 )then
              x0 = xa(0)+i1*dx(0) -xp0
              x1 = xa(1)+i2*dx(1) -xp1
                x2 = xa(2)+i3*dx(2) -xp2
              r2 = ( x0 )**2 + ( x1 )**2 + ( x2 )**2
           ! only evaluate if this is a valid point and pulse value is not too small
           if( mask(i1,i2,i3).gt.0 .and. (r2 .lt. tolGaussian) )then
           ! if( .true. )then
            phi = ampE*exp( -beta*( r2 ) )
            phiE = phi*g
               f(i1,i2,i3,ex)=f(i1,i2,i3,ex) + (x2-x1)*phiE
               f(i1,i2,i3,ey)=f(i1,i2,i3,ey) + (x0-x2)*phiE
               f(i1,i2,i3,ez)=f(i1,i2,i3,ez) + (x1-x0)*phiE
            else
             icount=icount+1 ! this point not evaluated
            end if
           end if
           end do
           end do
           end do
           if( t.le.5*dt )then
           write(*,'("Gaussian source: total points=",i6," points 
     & evaluated=",i6," : r2>tolGaussian=",e8.2)') (n1b-n1a+1)*(n2b-
     & n2a+1)*(n3b-n3a+1),(n1b-n1a+1)*(n2b-n2a+1)*(n3b-n3a+1)-icount,
     & tolGaussian
           end if

        else
          write(*,'("forcingOptMaxwell: case not implemented")')
          stop 7729
        end if

       else if( orderOfAccuracyInTime.eq.4 )then
        ! ***** 4th order in time *****
        if( gridType.eq.curvilinear .and. nd.eq.2 )then
          icount=0
          do i3=n3a,n3b
          do i2=n2a,n2b
          do i1=n1a,n1b
          if( mask(i1,i2,i3).gt.0 )then
              x0 = xy(i1,i2,i3,0)-xp0
              x1 = xy(i1,i2,i3,1)-xp1
              r2 = ( x0 )**2 + ( x1 )**2
           ! only evaluate if this is a valid point and pulse value is not too small
           if( mask(i1,i2,i3).gt.0 .and. (r2 .lt. tolGaussian) )then
           ! if( .true. )then
            phi = ampE*exp( -beta*( r2 ) )
            phiE = phi*g
             ! for 4th order in time we need to add corrections for the modified equation time-stepping
             !     f + (dt^2/12)*( c^2 Delta(f) + f_tt ) 
             phix = -2.*beta*x0*phi
             phiy = -2.*beta*x1*phi
             phixx = -2.*beta*phi + (2.*beta*x0)**2 * phi
             phiyy = -2.*beta*phi + (2.*beta*x1)**2 * phi
             phiEtt = phi*gtt
               lapPhi = phixx+phiyy
               lapF1   = -x1*lapPhi - 2.*phiy  ! Delta( -(y-y0)*phi )
               lapF2   =  x0*lapPhi + 2.*phix  ! Delta(  (x-x0)*phi )
               f(i1,i2,i3,ex)=f(i1,i2,i3,ex) -x1*phiE + dtsqBy12*( csq*
     & lapF1*g - x1*phiEtt )
               f(i1,i2,i3,ey)=f(i1,i2,i3,ey) +x0*phiE + dtsqBy12*( csq*
     & lapF2*g + x0*phiEtt )
               ! phiH = phi*cost
               ! lapPhiH = lapPhi *cost
               ! Do this:
               phiH = phi*g
               lapPhiH = lapPhi *g
               f(i1,i2,i3,hz)= f(i1,i2,i3,hz) + (ampH/ampE)*( phiH + 
     & dtsqBy12*( csq*lapPhiH + phi*gtt ) )
            else
             icount=icount+1 ! this point not evaluated
            end if
           end if
           end do
           end do
           end do
           if( t.le.5*dt )then
           write(*,'("Gaussian source: total points=",i6," points 
     & evaluated=",i6," : r2>tolGaussian=",e8.2)') (n1b-n1a+1)*(n2b-
     & n2a+1)*(n3b-n3a+1),(n1b-n1a+1)*(n2b-n2a+1)*(n3b-n3a+1)-icount,
     & tolGaussian
           end if

        else if( gridType.eq.rectangular .and. nd.eq.2 )then
          icount=0
          do i3=n3a,n3b
          do i2=n2a,n2b
          do i1=n1a,n1b
          if( mask(i1,i2,i3).gt.0 )then
              x0 = xa(0)+i1*dx(0) -xp0
              x1 = xa(1)+i2*dx(1) -xp1
              r2 = ( x0 )**2 + ( x1 )**2
           ! only evaluate if this is a valid point and pulse value is not too small
           if( mask(i1,i2,i3).gt.0 .and. (r2 .lt. tolGaussian) )then
           ! if( .true. )then
            phi = ampE*exp( -beta*( r2 ) )
            phiE = phi*g
             ! for 4th order in time we need to add corrections for the modified equation time-stepping
             !     f + (dt^2/12)*( c^2 Delta(f) + f_tt ) 
             phix = -2.*beta*x0*phi
             phiy = -2.*beta*x1*phi
             phixx = -2.*beta*phi + (2.*beta*x0)**2 * phi
             phiyy = -2.*beta*phi + (2.*beta*x1)**2 * phi
             phiEtt = phi*gtt
               lapPhi = phixx+phiyy
               lapF1   = -x1*lapPhi - 2.*phiy  ! Delta( -(y-y0)*phi )
               lapF2   =  x0*lapPhi + 2.*phix  ! Delta(  (x-x0)*phi )
               f(i1,i2,i3,ex)=f(i1,i2,i3,ex) -x1*phiE + dtsqBy12*( csq*
     & lapF1*g - x1*phiEtt )
               f(i1,i2,i3,ey)=f(i1,i2,i3,ey) +x0*phiE + dtsqBy12*( csq*
     & lapF2*g + x0*phiEtt )
               ! phiH = phi*cost
               ! lapPhiH = lapPhi *cost
               ! Do this:
               phiH = phi*g
               lapPhiH = lapPhi *g
               f(i1,i2,i3,hz)= f(i1,i2,i3,hz) + (ampH/ampE)*( phiH + 
     & dtsqBy12*( csq*lapPhiH + phi*gtt ) )
            else
             icount=icount+1 ! this point not evaluated
            end if
           end if
           end do
           end do
           end do
           if( t.le.5*dt )then
           write(*,'("Gaussian source: total points=",i6," points 
     & evaluated=",i6," : r2>tolGaussian=",e8.2)') (n1b-n1a+1)*(n2b-
     & n2a+1)*(n3b-n3a+1),(n1b-n1a+1)*(n2b-n2a+1)*(n3b-n3a+1)-icount,
     & tolGaussian
           end if

        else if( gridType.eq.curvilinear .and. nd.eq.3 )then
          icount=0
          do i3=n3a,n3b
          do i2=n2a,n2b
          do i1=n1a,n1b
          if( mask(i1,i2,i3).gt.0 )then
              x0 = xy(i1,i2,i3,0)-xp0
              x1 = xy(i1,i2,i3,1)-xp1
                x2 = xy(i1,i2,i3,2)-xp2
              r2 = ( x0 )**2 + ( x1 )**2 + ( x2 )**2
           ! only evaluate if this is a valid point and pulse value is not too small
           if( mask(i1,i2,i3).gt.0 .and. (r2 .lt. tolGaussian) )then
           ! if( .true. )then
            phi = ampE*exp( -beta*( r2 ) )
            phiE = phi*g
             ! for 4th order in time we need to add corrections for the modified equation time-stepping
             !     f + (dt^2/12)*( c^2 Delta(f) + f_tt ) 
             phix = -2.*beta*x0*phi
             phiy = -2.*beta*x1*phi
             phixx = -2.*beta*phi + (2.*beta*x0)**2 * phi
             phiyy = -2.*beta*phi + (2.*beta*x1)**2 * phi
             phiEtt = phi*gtt
               phiz = -2.*beta*x2*phi
               phizz = -2.*beta*phi + (2.*beta*x2)**2 * phi
               lapPhi = phixx+phiyy+phizz
               lapF1   = (x2-x1)*lapPhi + 2.*(phiz-phiy)        ! Delta( [(z-z0)-(y-y0)] *phi )
               lapF2   = (x0-x2)*lapPhi + 2.*(phix-phiz)        ! Delta( [(x-x0)-(z-z0)] *phi )
               lapF3   = (x1-x0)*lapPhi + 2.*(phiy-phix)        ! Delta( [(y-y0)-(x-x0)] *phi )
               f(i1,i2,i3,ex)=f(i1,i2,i3,ex) + (x2-x1)*phiE + dtsqBy12*
     & ( csq*lapF1*g + (x2-x1)*phiEtt )
               f(i1,i2,i3,ey)=f(i1,i2,i3,ey) + (x0-x2)*phiE + dtsqBy12*
     & ( csq*lapF2*g + (x0-x2)*phiEtt )
               f(i1,i2,i3,ez)=f(i1,i2,i3,ez) + (x1-x0)*phiE + dtsqBy12*
     & ( csq*lapF3*g + (x1-x0)*phiEtt )
            else
             icount=icount+1 ! this point not evaluated
            end if
           end if
           end do
           end do
           end do
           if( t.le.5*dt )then
           write(*,'("Gaussian source: total points=",i6," points 
     & evaluated=",i6," : r2>tolGaussian=",e8.2)') (n1b-n1a+1)*(n2b-
     & n2a+1)*(n3b-n3a+1),(n1b-n1a+1)*(n2b-n2a+1)*(n3b-n3a+1)-icount,
     & tolGaussian
           end if

        else if( gridType.eq.rectangular .and. nd.eq.3 )then
          icount=0
          do i3=n3a,n3b
          do i2=n2a,n2b
          do i1=n1a,n1b
          if( mask(i1,i2,i3).gt.0 )then
              x0 = xa(0)+i1*dx(0) -xp0
              x1 = xa(1)+i2*dx(1) -xp1
                x2 = xa(2)+i3*dx(2) -xp2
              r2 = ( x0 )**2 + ( x1 )**2 + ( x2 )**2
           ! only evaluate if this is a valid point and pulse value is not too small
           if( mask(i1,i2,i3).gt.0 .and. (r2 .lt. tolGaussian) )then
           ! if( .true. )then
            phi = ampE*exp( -beta*( r2 ) )
            phiE = phi*g
             ! for 4th order in time we need to add corrections for the modified equation time-stepping
             !     f + (dt^2/12)*( c^2 Delta(f) + f_tt ) 
             phix = -2.*beta*x0*phi
             phiy = -2.*beta*x1*phi
             phixx = -2.*beta*phi + (2.*beta*x0)**2 * phi
             phiyy = -2.*beta*phi + (2.*beta*x1)**2 * phi
             phiEtt = phi*gtt
               phiz = -2.*beta*x2*phi
               phizz = -2.*beta*phi + (2.*beta*x2)**2 * phi
               lapPhi = phixx+phiyy+phizz
               lapF1   = (x2-x1)*lapPhi + 2.*(phiz-phiy)        ! Delta( [(z-z0)-(y-y0)] *phi )
               lapF2   = (x0-x2)*lapPhi + 2.*(phix-phiz)        ! Delta( [(x-x0)-(z-z0)] *phi )
               lapF3   = (x1-x0)*lapPhi + 2.*(phiy-phix)        ! Delta( [(y-y0)-(x-x0)] *phi )
               f(i1,i2,i3,ex)=f(i1,i2,i3,ex) + (x2-x1)*phiE + dtsqBy12*
     & ( csq*lapF1*g + (x2-x1)*phiEtt )
               f(i1,i2,i3,ey)=f(i1,i2,i3,ey) + (x0-x2)*phiE + dtsqBy12*
     & ( csq*lapF2*g + (x0-x2)*phiEtt )
               f(i1,i2,i3,ez)=f(i1,i2,i3,ez) + (x1-x0)*phiE + dtsqBy12*
     & ( csq*lapF3*g + (x1-x0)*phiEtt )
            else
             icount=icount+1 ! this point not evaluated
            end if
           end if
           end do
           end do
           end do
           if( t.le.5*dt )then
           write(*,'("Gaussian source: total points=",i6," points 
     & evaluated=",i6," : r2>tolGaussian=",e8.2)') (n1b-n1a+1)*(n2b-
     & n2a+1)*(n3b-n3a+1),(n1b-n1a+1)*(n2b-n2a+1)*(n3b-n3a+1)-icount,
     & tolGaussian
           end if

        else
          write(*,'("forcingOptMaxwell: case not implemented")')
          stop 7729
        end if

       else
        write(*,'("forcingOptMaxwell: unknown orderOfAccuracyInTime",
     & i6)') orderOfAccuracyInTime
        stop 5219
       end if


      else
        write(*,'("forcingOptMaxwell: unknown forcingOption=",i6)') 
     & forcingOption
        stop 5219
      end if

      return
      end
