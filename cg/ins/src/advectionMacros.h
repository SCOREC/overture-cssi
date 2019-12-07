!  -*- mode: F90 -*-
! Define macros to evaluate an advection term (a.grad) u
!

#beginMacro getBWENOCoeffs(s1,s2,s3,var,vard,drl)

  Apl   = u(i1+s1,i2+s2,i3+s3,var) - 2.*u(i1   ,i2,i3,var)   + u(i1-s1,i2-s2,i3-s3,var)
  Bpl   = u(i1+s1,i2+s2,i3+s3,var) -    u(i1-s1,i2-s2,i3-s3,var)
  
  Apr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) - 2.*u(i1+s1,i2+s2,i3+s3,var) + u(i1,i2,i3,var)
  Bpr   = u(i1+2*s1,i2+2*s2,i3+2*s3,var) -    u(i1   ,i2   ,i3   ,var)
  
  Aml   = u(i1,i2,i3,var) - 2.*u(i1-s1  ,i2-s2  ,i3-s3  ,var) + u(i1-2*s1,i2-2*s2,i3-2*s3,var)
  Bml   = u(i1,i2,i3,var) -    u(i1-2*s1,i2-2*s2,i3-2*s3,var)
  Amr   = Apl
  Bmr   = Bpl

  betapl= 4./3.*Apl**2 + 1./2.*Apl*Bpl + 1./4.*Bpl**2
  betapr= 4./3.*Apr**2 - 1./2.*Apr*Bpr + 1./4.*Bpr**2
  
  betaml= 4./3.*Aml**2 + 1./2.*Aml*Bml + 1./4.*Bml**2
  betamr= 4./3.*Amr**2 - 1./2.*Amr*Bmr + 1./4.*Bmr**2
  
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
        
        aumax = max(abs(au),abs(aup),abs(aum))
        
        expAd = (1./12.)*aumax*( u(i1-2*s1,i2-2*s2,i3-2*s3,var) - 4.*u(i1-s1,i2-s2,i3-s3,var) + 6.*u(i1,i2,i3,var) - 4.*u(i1+s1,i2+s2,i3+s3,var) +  u(i1+2*s1,i2+2*s2,i3+2*s3,var))

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
        
  Fpl = 1./6.*(  -u(i1-s1,i2-s2,i3-s3,var) + 5.*u(i1,i2,i3,var)   + 2.*u(i1+s1,i2+s2,i3+s3,var))
  Fpr = 1./6.*( 2.*u(i1,i2,i3,var)   + 5.*u(i1+s1,i2+s2,i3+s3,var) -   u(i1+2*s1,i2+2*s2,i3+2*s3,var))
  
  Fml = 1./6.*(  -u(i1-2*s1,i2-2*s2,i3-2*s3,var) + 5.*u(i1-s1,i2-s2,i3-s3,var)   + 2.*u(i1,i2,i3,var))
  Fmr = 1./6.*( 2.*u(i1-s1,i2-s2,i3-s3,var) + 5.*u(i1,i2,i3,var) -       u(i1+s1,i2+s2,i3+s3,var))
  
  Fp  = wpl*Fpl + wpr*Fpr
  Fm  = wml*Fml + wmr*Fmr
  
  
     agu(vard,var)= au*(Fp - Fm)/drl + expAd/drl

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
  
#endMacro

! --------------------------------------------------------------------------------
!  Macro: getUpwindAdvection
! 
! --------------------------------------------------------------------------------
#beginMacro getUpwindAdvection(u,uu,i1,i2,i3,SCALAR,DIM,ORDER,GRIDTYPE, agu)

     
 #If #GRIDTYPE eq "rectangular"
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
     agu(uc,uc)= au*(u(i1,i2,i3,uc)-u(i1-1,i2,i3,uc))/(dx(0)) 
     ! u*vx = u*D-x(v)
     agu(uc,vc)= au*(u(i1,i2,i3,vc)-u(i1-1,i2,i3,vc))/(dx(0))
   else
     ! u*ux = u*D+x(u)
     agu(uc,uc)= au*(u(i1+1,i2,i3,uc)-u(i1,i2,i3,uc))/(dx(0))
     ! u*vx = u*D+x(v)
     agu(uc,vc)= au*(u(i1+1,i2,i3,vc)-u(i1,i2,i3,vc))/(dx(0))
   end if
 
   au = uu(i1,i2,i3,vc)
   if( au.gt.0. )then
     ! v*uy = v*D-y(u)
     agu(vc,uc)= au*(u(i1,i2,i3,uc)-u(i1,i2-1,i3,uc))/(dx(1))
     ! v*vy = v*D-y(v)
     agu(vc,vc)= au*(u(i1,i2,i3,vc)-u(i1,i2-1,i3,vc))/(dx(1))
   else
     ! v*uy = v*D+y(u)
     agu(vc,uc)= au*(u(i1,i2+1,i3,uc)-u(i1,i2,i3,uc))/(dx(1))
     ! v*vy = v*D+y(v) 
     agu(vc,vc)= au*(u(i1,i2+1,i3,vc)-u(i1,i2,i3,vc))/(dx(1))
   end if

#If #DIM eq "3"
      ! finish me 

!      stop 777
   au = uu(i1,i2,i3,uc)
   if( au.gt.0. )then
      ! u*ux = u*D-x(u)
      agu(uc,uc)= au*(u(i1,i2,i3,uc)-u(i1-1,i2,i3,uc))/(dx(0)) 
      ! u*vx = u*D-x(v)
      agu(uc,vc)= au*(u(i1,i2,i3,vc)-u(i1-1,i2,i3,vc))/(dx(0))
      ! u*wx = u*D-x(w)
      agu(uc,wc)= au*(u(i1,i2,i3,wc)-u(i1-1,i2,i3,wc))/(dx(0))
   else
      ! u*ux = u*D+x(u)
      agu(uc,uc)= au*(u(i1+1,i2,i3,uc)-u(i1,i2,i3,uc))/(dx(0))
      ! u*vx = u*D+x(v)
      agu(uc,vc)= au*(u(i1+1,i2,i3,vc)-u(i1,i2,i3,vc))/(dx(0))
      ! u*wx = u*D+x(w)
      agu(uc,wc)= au*(u(i1+1,i2,i3,wc)-u(i1,i2,i3,wc))/(dx(0))
   end if
   
   au = uu(i1,i2,i3,vc)
   if( au.gt.0. )then
      ! v*uy = v*D-y(u)
      agu(vc,uc)= au*(u(i1,i2,i3,uc)-u(i1,i2-1,i3,uc))/(dx(1))
      ! v*vy = v*D-y(v)
      agu(vc,vc)= au*(u(i1,i2,i3,vc)-u(i1,i2-1,i3,vc))/(dx(1))
      ! v*wy = v*D-y(w)
      agu(vc,wc)= au*(u(i1,i2,i3,wc)-u(i1,i2-1,i3,wc))/(dx(1))
   else
      ! v*uy = v*D+y(u)
      agu(vc,uc)= au*(u(i1,i2+1,i3,uc)-u(i1,i2,i3,uc))/(dx(1))
      ! v*vy = v*D+y(v) 
      agu(vc,vc)= au*(u(i1,i2+1,i3,vc)-u(i1,i2,i3,vc))/(dx(1))
      ! v*wy = v*D+y(w)
      agu(vc,wc)= au*(u(i1,i2+1,i3,wc)-u(i1,i2,i3,wc))/(dx(1))
   end if
   
   au = uu(i1,i2,i3,wc)
   if( au.gt.0. )then
      ! w*uy = w*D-y(u)
      agu(wc,uc)= au*(u(i1,i2,i3,uc)-u(i1,i2,i3-1,uc))/(dx(2))
      ! w*vy = w*D-y(v)
      agu(wc,vc)= au*(u(i1,i2,i3,vc)-u(i1,i2,i3-1,vc))/(dx(2))
      ! w*wy = w*D-y(w)
      agu(wc,wc)= au*(u(i1,i2,i3,wc)-u(i1,i2,i3-1,wc))/(dx(2))
   else
      ! w*uy = w*D+y(u)
      agu(wc,uc)= au*(u(i1,i2,i3+1,uc)-u(i1,i2,i3,uc))/(dx(2))
      ! w*vy = w*D+y(v)
      agu(wc,vc)= au*(u(i1,i2,i3+1,vc)-u(i1,i2,i3,vc))/(dx(2))
      ! w*wy = w*D+y(w)
      agu(wc,wc)= au*(u(i1,i2,i3+1,wc)-u(i1,i2,i3,wc))/(dx(2))
   end if

#End

  elseif( upwindOrder.eq.2 )then
    ! write(*,'(" finish me, upwindOrder=",i2)') upwindOrder
    ! stop 222

     au = uu(i1,i2,i3,uc)
     if( au.gt.0. )then
        ! u*ux = u*D-x(u)
        agu(uc,uc)= au*(3.*u(i1,i2,i3,uc)-4.*u(i1-1,i2,i3,uc)+u(i1-2,i2,i3,uc))/(2.*dx(0)) 
        ! u*vx = u*D-x(v)
        agu(uc,vc)= au*(3.*u(i1,i2,i3,vc)-4.*u(i1-1,i2,i3,vc)+u(i1-2,i2,i3,vc))/(2.*dx(0))
     else
        ! u*ux = u*D+x(u)
        agu(uc,uc)= au*(-u(i1+2,i2,i3,uc)+4.*u(i1+1,i2,i3,uc)-3.*u(i1,i2,i3,uc))/(2.*dx(0))
        ! u*vx = u*D+x(v)
        agu(uc,vc)= au*(-u(i1+2,i2,i3,vc)+4.*u(i1+1,i2,i3,vc)-3.*u(i1,i2,i3,vc))/(2.*dx(0))
     end if
     
     au = uu(i1,i2,i3,vc)
     if( au.gt.0. )then
        ! v*uy = v*D-y(u)
        agu(vc,uc)= au*(3.*u(i1,i2,i3,uc)-4.*u(i1,i2-1,i3,uc)+u(i1,i2-2,i3,uc))/(2.*dx(1))
        ! v*vy = v*D-y(v)
        agu(vc,vc)= au*(3.*u(i1,i2,i3,vc)-4.*u(i1,i2-1,i3,vc)+u(i1,i2-2,i3,vc))/(2.*dx(1))
     else
        ! v*uy = v*D+y(u)
        agu(vc,uc)= au*(-u(i1,i2+2,i3,uc)+4.*u(i1,i2+1,i3,uc)-3.*u(i1,i2,i3,uc))/(2.*dx(1))
        ! v*vy = v*D+y(v) 
        agu(vc,vc)= au*(-u(i1,i2+2,i3,vc)+4.*u(i1,i2+1,i3,vc)-3.*u(i1,i2,i3,vc))/(2.*dx(1))
     end if
     
    #If #DIM eq "3"
      ! finish me 

      !stop 777
     au = uu(i1,i2,i3,uc)
     if( au.gt.0. )then
        ! u*ux = u*D-x(u)
        agu(uc,uc)= au*(3.*u(i1,i2,i3,uc)-4.*u(i1-1,i2,i3,uc)+u(i1-2,i2,i3,uc))/(2.*dx(0)) 
        ! u*vx = u*D-x(v)
        agu(uc,vc)= au*(3.*u(i1,i2,i3,vc)-4.*u(i1-1,i2,i3,vc)+u(i1-2,i2,i3,vc))/(2.*dx(0))
        ! u*vx = u*D-x(v)
        agu(uc,wc)= au*(3.*u(i1,i2,i3,wc)-4.*u(i1-1,i2,i3,wc)+u(i1-2,i2,i3,wc))/(2.*dx(0))
     else
        ! u*ux = u*D+x(u)
        agu(uc,uc)= au*(-u(i1+2,i2,i3,uc)+4.*u(i1+1,i2,i3,uc)-3.*u(i1,i2,i3,uc))/(2.*dx(0))
        ! u*vx = u*D+x(v)
        agu(uc,vc)= au*(-u(i1+2,i2,i3,vc)+4.*u(i1+1,i2,i3,vc)-3.*u(i1,i2,i3,vc))/(2.*dx(0))
        ! u*vx = u*D-x(v)
        agu(uc,wc)= au*(-u(i1+2,i2,i3,wc)+4.*u(i1+1,i2,i3,wc)-3.*u(i1,i2,i3,wc))/(2.*dx(0))
     end if

     au = uu(i1,i2,i3,vc)
     if( au.gt.0. )then
        ! v*uy = v*D-y(u)
        agu(vc,uc)= au*(3.*u(i1,i2,i3,uc)-4.*u(i1,i2-1,i3,uc)+u(i1,i2-2,i3,uc))/(2.*dx(1))
        ! v*vy = v*D-y(v)
        agu(vc,vc)= au*(3.*u(i1,i2,i3,vc)-4.*u(i1,i2-1,i3,vc)+u(i1,i2-2,i3,vc))/(2.*dx(1))
        ! v*vy = v*D-y(v)
        agu(vc,wc)= au*(3.*u(i1,i2,i3,wc)-4.*u(i1,i2-1,i3,wc)+u(i1,i2-2,i3,wc))/(2.*dx(1))
     else
        ! v*uy = v*D+y(u)
        agu(vc,uc)= au*(-u(i1,i2+2,i3,uc)+4.*u(i1,i2+1,i3,uc)-3.*u(i1,i2,i3,uc))/(2.*dx(1))
        ! v*vy = v*D+y(v) 
        agu(vc,vc)= au*(-u(i1,i2+2,i3,vc)+4.*u(i1,i2+1,i3,vc)-3.*u(i1,i2,i3,vc))/(2.*dx(1))
        ! v*vy = v*D+y(v) 
        agu(vc,wc)= au*(-u(i1,i2+2,i3,wc)+4.*u(i1,i2+1,i3,wc)-3.*u(i1,i2,i3,wc))/(2.*dx(1))
     end if

     au = uu(i1,i2,i3,wc)
     if( au.gt.0. )then
        ! v*uy = v*D-y(u)
        agu(wc,uc)= au*(3.*u(i1,i2,i3,uc)-4.*u(i1,i2,i3-1,uc)+u(i1,i2,i3-2,uc))/(2.*dx(2))
        ! v*vy = v*D-y(v)
        agu(wc,vc)= au*(3.*u(i1,i2,i3,vc)-4.*u(i1,i2,i3-1,vc)+u(i1,i2,i3-2,vc))/(2.*dx(2))
        ! v*vy = v*D-y(v)
        agu(wc,wc)= au*(3.*u(i1,i2,i3,wc)-4.*u(i1,i2,i3-1,wc)+u(i1,i2,i3-2,wc))/(2.*dx(2))
     else
        ! v*uy = v*D+y(u)
        agu(wc,uc)= au*(-u(i1,i2,i3+2,uc)+4.*u(i1,i2,i3+1,uc)-3.*u(i1,i2,i3,uc))/(2.*dx(2))
        ! v*vy = v*D+y(v) 
        agu(wc,vc)= au*(-u(i1,i2,i3+2,vc)+4.*u(i1,i2,i3+1,vc)-3.*u(i1,i2,i3,vc))/(2.*dx(2))
        ! v*vy = v*D+y(v) 
        agu(wc,wc)= au*(-u(i1,i2,i3+2,wc)+4.*u(i1,i2,i3+1,wc)-3.*u(i1,i2,i3,wc))/(2.*dx(2))
     end if
     
    #End

  elseif( upwindOrder.eq.3 )then

     au = uu(i1,i2,i3,uc)
     if( au.gt.0. )then
        ! u*ux = u*D-x(u)
        agu(uc,uc)= au*(2.*u(i1+1,i2,i3,uc)+3.*u(i1,i2,i3,uc)-6.*u(i1-1,i2,i3,uc)+u(i1-2,i2,i3,uc))/(6.*dx(0)) 
        ! u*vx = u*D-x(v)
        agu(uc,vc)= au*(2.*u(i1+1,i2,i3,vc)+3.*u(i1,i2,i3,vc)-6.*u(i1-1,i2,i3,vc)+u(i1-2,i2,i3,vc))/(6.*dx(0))
     else
        ! u*ux = u*D+x(u)
        agu(uc,uc)= au*(-u(i1+2,i2,i3,uc)+6.*u(i1+1,i2,i3,uc)-3.*u(i1,i2,i3,uc)-2.*u(i1-1,i2,i3,uc))/(6.*dx(0))
        ! u*vx = u*D+x(v)
        agu(uc,vc)= au*(-u(i1+2,i2,i3,vc)+6.*u(i1+1,i2,i3,vc)-3.*u(i1,i2,i3,vc)-2.*u(i1-1,i2,i3,vc))/(6.*dx(0))
     end if
     
     au = uu(i1,i2,i3,vc)
     if( au.gt.0. )then
        ! v*uy = v*D-y(u)
        agu(vc,uc)= au*(2.*u(i1,i2+1,i3,uc)+3.*u(i1,i2,i3,uc)-6.*u(i1,i2-1,i3,uc)+u(i1,i2-2,i3,uc))/(6.*dx(1))
        ! v*vy = v*D-y(v)
        agu(vc,vc)= au*(2.*u(i1,i2+1,i3,vc)+3.*u(i1,i2,i3,vc)-6.*u(i1,i2-1,i3,vc)+u(i1,i2-2,i3,vc))/(6.*dx(1))
     else
        ! v*uy = v*D+y(u)
        agu(vc,uc)= au*(-u(i1,i2+2,i3,uc)+6.*u(i1,i2+1,i3,uc)-3.*u(i1,i2,i3,uc)-2.*u(i1,i2-1,i3,uc))/(6.*dx(1))
        ! v*vy = v*D+y(v) 
        agu(vc,vc)= au*(-u(i1,i2+2,i3,vc)+6.*u(i1,i2+1,i3,vc)-3.*u(i1,i2,i3,vc)-2.*u(i1,i2-1,i3,vc))/(6.*dx(1))
     end if
     
    #If #DIM eq "3"

     au = uu(i1,i2,i3,uc)
     ! u*wx
     if( au.gt.0. )then

        agu(uc,wc)= au*(2.*u(i1+1,i2,i3,wc)+3.*u(i1,i2,i3,wc)-6.*u(i1-1,i2,i3,wc)+u(i1-2,i2,i3,wc))/(6.*dx(0)) 

     else

        agu(uc,wc)= au*(-u(i1+2,i2,i3,wc)+6.*u(i1+1,i2,i3,wc)-3.*u(i1,i2,i3,wc)-2.*u(i1-1,i2,i3,wc))/(6.*dx(0))

     end if

     
     au = uu(i1,i2,i3,vc)
     ! v*wy
     if( au.gt.0. )then

        agu(vc,wc)= au*(2.*u(i1,i2+1,i3,wc)+3.*u(i1,i2,i3,wc)-6.*u(i1,i2-1,i3,wc)+u(i1,i2-2,i3,wc))/(6.*dx(1))

     else
        agu(vc,wc)= au*(-u(i1,i2+2,i3,wc)+6.*u(i1,i2+1,i3,wc)-3.*u(i1,i2,i3,wc)-2.*u(i1,i2-1,i3,wc))/(6.*dx(1))

     end if

     au = uu(i1,i2,i3,wc)
     ! w* \f{\p u_i}{ \p z}   i = 1,2,3
     if( au.gt.0. )then

        agu(wc,uc)= au*(2.*u(i1,i2,i3+1,uc)+3.*u(i1,i2,i3,uc)-6.*u(i1,i2,i3-1,uc)+u(i1,i2,i3-2,uc))/(6.*dx(2))

        agu(wc,vc)= au*(2.*u(i1,i2,i3+1,vc)+3.*u(i1,i2,i3,vc)-6.*u(i1,i2,i3-1,vc)+u(i1,i2,i3-2,vc))/(6.*dx(2))

        agu(wc,wc)= au*(2.*u(i1,i2,i3+1,wc)+3.*u(i1,i2,i3,wc)-6.*u(i1,i2,i3-1,wc)+u(i1,i2,i3-2,wc))/(6.*dx(2))

        
     else

        agu(wc,uc)= au*(-u(i1,i2,i3+2,uc)+6.*u(i1,i2,i3+1,uc)-3.*u(i1,i2,i3,uc)-2.*u(i1,i2,i3-1,uc))/(6.*dx(2))

        agu(wc,vc)= au*(-u(i1,i2,i3+2,vc)+6.*u(i1,i2,i3+1,vc)-3.*u(i1,i2,i3,vc)-2.*u(i1,i2,i3-1,vc))/(6.*dx(2))
        
        agu(wc,wc)= au*(-u(i1,i2,i3+2,wc)+6.*u(i1,i2,i3+1,wc)-3.*u(i1,i2,i3,wc)-2.*u(i1,i2,i3-1,wc))/(6.*dx(2))

     end if
     
     
    #End
    
  end if

 #Elif #GRIDTYPE eq "curvilinear"

  ! write(*,'(" we in HERE")' )
  ! stop 7171
  
  
  ! --- CURVILINEAR GRID ---
      if( upwindOrder.eq.1 )then

         au   = rsxy(i1,i2,i3,0,0)*uu(i1,i2,i3,uc)+rsxy(i1,i2,i3,0,1)*uu(i1,i2,i3,vc)
         
         if( au.gt.0. )then
            ! u*ux
            agu(uc,uc)= au*(u(i1,i2,i3,uc)-u(i1-1,i2,i3,uc))/(dr(0))
            ! u*vx
            agu(uc,vc)= au*(u(i1,i2,i3,vc)-u(i1-1,i2,i3,vc))/(dr(0))
         else
            ! u*ux
            agu(uc,uc)= au*(u(i1+1,i2,i3,uc)-u(i1,i2,i3,uc))/(dr(0))
            ! u*vx
            agu(uc,vc)= au*(u(i1+1,i2,i3,vc)-u(i1,i2,i3,vc))/(dr(0))
         end if
         
         au = rsxy(i1,i2,i3,1,0)*uu(i1,i2,i3,uc)+rsxy(i1,i2,i3,1,1)*uu(i1,i2,i3,vc)
         
         if( au.gt.0. )then
            ! v*uy
            agu(vc,uc)= au*(u(i1,i2,i3,uc)-u(i1,i2-1,i3,uc))/(dr(1))
            ! v*vy
            agu(vc,vc)= au*(u(i1,i2,i3,vc)-u(i1,i2-1,i3,vc))/(dr(1))
         else
            ! v*uy
            agu(vc,uc)= au*(u(i1,i2+1,i3,uc)-u(i1,i2,i3,uc))/(dr(1))
            ! v*vy
            agu(vc,vc)= au*(u(i1,i2+1,i3,vc)-u(i1,i2,i3,vc))/(dr(1))
         end if

      elseif( upwindOrder.eq.2 )then

         au = rsxy(i1,i2,i3,0,0)*uu(i1,i2,i3,uc)+rsxy(i1,i2,i3,0,1)*uu(i1,i2,i3,vc)

         if( au.gt.0. )then
          ! u*ux = u*D-x(u)
            agu(uc,uc)= au*(3.*u(i1,i2,i3,uc)-4.*u(i1-1,i2,i3,uc)+u(i1-2,i2,i3,uc))/(2.*dr(0)) 
          ! u*vx = u*D-x(v)
            agu(uc,vc)= au*(3.*u(i1,i2,i3,vc)-4.*u(i1-1,i2,i3,vc)+u(i1-2,i2,i3,vc))/(2.*dr(0))
         else
          ! u*ux = u*D+x(u)
            agu(uc,uc)= au*(-u(i1+2,i2,i3,uc)+4.*u(i1+1,i2,i3,uc)-3.*u(i1,i2,i3,uc))/(2.*dr(0))
          ! u*vx = u*D+x(v)
            agu(uc,vc)= au*(-u(i1+2,i2,i3,vc)+4.*u(i1+1,i2,i3,vc)-3.*u(i1,i2,i3,vc))/(2.*dr(0))
         end if

         au = rsxy(i1,i2,i3,1,0)*uu(i1,i2,i3,uc)+rsxy(i1,i2,i3,1,1)*uu(i1,i2,i3,vc)

         if( au.gt.0. )then
          ! v*uy = v*D-y(u)
            agu(vc,uc)= au*(3.*u(i1,i2,i3,uc)-4.*u(i1,i2-1,i3,uc)+u(i1,i2-2,i3,uc))/(2.*dr(1))
          ! v*vy = v*D-y(v)
            agu(vc,vc)= au*(3.*u(i1,i2,i3,vc)-4.*u(i1,i2-1,i3,vc)+u(i1,i2-2,i3,vc))/(2.*dr(1))
         else
          ! v*uy = v*D+y(u)
            agu(vc,uc)= au*(-u(i1,i2+2,i3,uc)+4.*u(i1,i2+1,i3,uc)-3.*u(i1,i2,i3,uc))/(2.*dr(1))
          ! v*vy = v*D+y(v) 
            agu(vc,vc)= au*(-u(i1,i2+2,i3,vc)+4.*u(i1,i2+1,i3,vc)-3.*u(i1,i2,i3,vc))/(2.*dr(1))
         end if

      elseif( upwindOrder.eq.3 )then

         au = rsxy(i1,i2,i3,0,0)*uu(i1,i2,i3,uc)+rsxy(i1,i2,i3,0,1)*uu(i1,i2,i3,vc)
         
         if( au.gt.0. )then
        ! u*ux = u*D-x(u)
            agu(uc,uc)= au*(2.*u(i1+1,i2,i3,uc)+3.*u(i1,i2,i3,uc)-6.*u(i1-1,i2,i3,uc)+u(i1-2,i2,i3,uc))/(6.*dr(0)) 
        ! u*vx = u*D-x(v)
            agu(uc,vc)= au*(2.*u(i1+1,i2,i3,vc)+3.*u(i1,i2,i3,vc)-6.*u(i1-1,i2,i3,vc)+u(i1-2,i2,i3,vc))/(6.*dr(0))
         else
        ! u*ux = u*D+x(u)
            agu(uc,uc)= au*(-u(i1+2,i2,i3,uc)+6.*u(i1+1,i2,i3,uc)-3.*u(i1,i2,i3,uc)-2.*u(i1-1,i2,i3,uc))/(6.*dr(0))
        ! u*vx = u*D+x(v)
            agu(uc,vc)= au*(-u(i1+2,i2,i3,vc)+6.*u(i1+1,i2,i3,vc)-3.*u(i1,i2,i3,vc)-2.*u(i1-1,i2,i3,vc))/(6.*dr(0))
         end if

         au = rsxy(i1,i2,i3,1,0)*uu(i1,i2,i3,uc)+rsxy(i1,i2,i3,1,1)*uu(i1,i2,i3,vc)

         if( au.gt.0. )then
        ! v*uy = v*D-y(u)
            agu(vc,uc)= au*(2.*u(i1,i2+1,i3,uc)+3.*u(i1,i2,i3,uc)-6.*u(i1,i2-1,i3,uc)+u(i1,i2-2,i3,uc))/(6.*dr(1))
        ! v*vy = v*D-y(v)
            agu(vc,vc)= au*(2.*u(i1,i2+1,i3,vc)+3.*u(i1,i2,i3,vc)-6.*u(i1,i2-1,i3,vc)+u(i1,i2-2,i3,vc))/(6.*dr(1))
         else
        ! v*uy = v*D+y(u)
            agu(vc,uc)= au*(-u(i1,i2+2,i3,uc)+6.*u(i1,i2+1,i3,uc)-3.*u(i1,i2,i3,uc)-2.*u(i1,i2-1,i3,uc))/(6.*dr(1))
        ! v*vy = v*D+y(v) 
            agu(vc,vc)= au*(-u(i1,i2+2,i3,vc)+6.*u(i1,i2+1,i3,vc)-3.*u(i1,i2,i3,vc)-2.*u(i1,i2-1,i3,vc))/(6.*dr(1))
         end if

      end if
!
!      
      #If #DIM eq "3"
!      ! finish me 
!
        stop 777
!      !  au = rsxy(i1,i2,i3,0,0)*u(i1,i2,i3,uc)+rsxy(i1,i2,i3,0,1)*u(i1,i2,i3,vc)+rsxy(i1,i2,i3,0,2)*u(i1,i2,i3,wc)
      #End
      
 #Else
   write(*,'(" getUpwindAdvection: finish me")' )
   stop 7171
 #End

#endMacro 

! --------------------------------------------------------------------------------
!  Macro: getBwenoAdvection
! --------------------------------------------------------------------------------
#beginMacro getBwenoAdvection(u,uu,i1,i2,i3,SCALAR,DIM,ORDER,GRIDTYPE, agu)

   gvU = u(i1,i2,i3,uc) - uu(i1,i2,i3,uc)
   gvV = u(i1,i2,i3,vc) - uu(i1,i2,i3,vc)
   gvW = u(i1,i2,i3,wc) - uu(i1,i2,i3,wc)
   
   
 #If #GRIDTYPE eq "rectangular"

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

      getBWENOCoeffs(s1,s2,s3,var,vard,drl)

      !u*vx
      s1=1
      s2=0
      s3=0
      var = vc
      vard= uc
      drl = dx(0)

      getBWENOCoeffs(s1,s2,s3,var,vard,drl)

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

      getBWENOCoeffs(s1,s2,s3,var,vard,drl)

      !v*vy

      s1=0
      s2=1
      s3=0
      var = vc
      vard= vc
      drl = dx(1)

      getBWENOCoeffs(s1,s2,s3,var,vard,drl)

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
     
    #If #DIM eq "3"
     ! -- fourth order bweno 3D Cartesian--

     if( upwindOrder.eq.4 )then

        !get u
        au = uu(i1,i2,i3,uc)
        aum = u(i1-1,i2,i3,uc) - gvU
        aup = u(i1+1,i2,i3,uc) - gvU

        !u*ux
        s1=1
        s2=0
        s3=0
        var = uc
        vard= uc
        drl = dx(0) 

        getBWENOCoeffs(s1,s2,s3,var,vard,drl)
        
        !u*vx
        s1=1
        s2=0
        s3=0
        var = vc
        vard= uc
        drl = dx(0)

        getBWENOCoeffs(s1,s2,s3,var,vard,drl)

        !u*wx
        s1=1
        s2=0
        s3=0
        var = wc
        vard= uc
        drl = dx(0) 

        getBWENOCoeffs(s1,s2,s3,var,vard,drl)

        !------------------------
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
        
        getBWENOCoeffs(s1,s2,s3,var,vard,drl)

        !v*vy
        s1=0
        s2=1
        s3=0
        var = vc
        vard= vc
        drl = dx(1) 
        
        getBWENOCoeffs(s1,s2,s3,var,vard,drl)

        !v*wy
        s1=0
        s2=1
        s3=0
        var = wc
        vard= vc
        drl = dx(1) 
        
        getBWENOCoeffs(s1,s2,s3,var,vard,drl)

        !------------------------
        !w*uz
        au = uu(i1,i2,i3,wc)
        aum = u(i1,i2,i3-1,wc) - gvW
        aup = u(i1,i2,i3+1,wc) - gvW

        s1=0
        s2=0
        s3=1
        var = uc
        vard= wc
        drl = dx(2)
        
        getBWENOCoeffs(s1,s2,s3,var,vard,drl)

        !w*vz
        s1=0
        s2=0
        s3=1
        var = vc
        vard= wc
        drl = dx(2)
        
        getBWENOCoeffs(s1,s2,s3,var,vard,drl)

        !w*wz
        s1=0
        s2=0
        s3=1
        var = wc
        vard= wc
        drl = dx(2)
        
        getBWENOCoeffs(s1,s2,s3,var,vard,drl)

     end if


    #End

     
  else ! if (upwindOrder.eq.4)
     write(*,'(" getBewnoAdvection: only 4th order is avaliable now")' )
     stop 777
     
  end if
  
 #Elif #GRIDTYPE eq "curvilinear" 
  !bweno curvulinear
  orderOneFix = 1
  augmentPlot = 0

  if( upwindOrder.eq.4 )then

     ! u*ux-------------------------------------------------------------------------------
     ! au = rsxy(i1,i2,i3,0,0)*u(i1,i2,i3,uc)+rsxy(i1,i2,i3,0,1)*u(i1,i2,i3,vc) + rsxy(i1,i2,i3,0,2)*u(i1,i2,i3,wc)
     ! aum   = rsxy(i1-1,i2,i3,0,0)*u(i1-1,i2,i3,uc)+rsxy(i1-1,i2,i3,0,1)*u(i1-1,i2,i3,vc)+ rsxy(i1-1,i2,i3,0,2)*u(i1-1,i2,i3,wc)
     ! aup   = rsxy(i1+1,i2,i3,0,0)*u(i1+1,i2,i3,uc)+rsxy(i1+1,i2,i3,0,1)*u(i1+1,i2,i3,vc)+ rsxy(i1+1,i2,i3,0,2)*u(i1+1,i2,i3,wc)

     au    = rsxy(i1  ,i2,i3,0,0)*uu(i1  ,i2,i3,uc) + rsxy(i1  ,i2,i3,0,1) * uu(i1  ,i2,i3,vc)
     aum   = rsxy(i1-1,i2,i3,0,0)*(u(i1-1,i2,i3,uc) - gvU ) + rsxy(i1-1,i2,i3,0,1) * (u(i1-1,i2,i3,vc)- gvV )
     aup   = rsxy(i1+1,i2,i3,0,0)*(u(i1+1,i2,i3,uc) - gvU ) + rsxy(i1+1,i2,i3,0,1) * (u(i1+1,i2,i3,vc)- gvV )
     
     s1=1
     s2=0
     s3=0
     var = uc
     vard= uc
     
     drl = dr(0)
     
     getBWENOCoeffs(s1,s2,s3,var,vard,drl)

     !u*vx
     s1=1
     s2=0
     s3=0
     var = vc
     vard= uc
     drl = dr(0) 
     
     getBWENOCoeffs(s1,s2,s3,var,vard,drl)

     !v*uy
     au    = rsxy(i1,i2  ,i3,1,0)*uu(i1,i2  ,i3,uc)+rsxy(i1,i2,i3,1,1)*uu(i1,i2,i3,vc)
     aum   = rsxy(i1,i2-1,i3,1,0)*(u(i1,i2-1,i3,uc) - gvU) + rsxy(i1,i2-1,i3,1,1)*(u(i1,i2-1,i3,vc) - gvV)
     aup   = rsxy(i1,i2+1,i3,1,0)*(u(i1,i2+1,i3,uc) - gvU) + rsxy(i1,i2+1,i3,1,1)*(u(i1,i2+1,i3,vc) - gvV)

     s1=0
     s2=1
     s3=0
     var = uc
     vard= vc
     drl = dr(1) 
     
     getBWENOCoeffs(s1,s2,s3,var,vard,drl)
     !v*vy
     s1=0
     s2=1
     s3=0
     var = vc
     vard= vc
     drl = dr(1) 
     
     getBWENOCoeffs(s1,s2,s3,var,vard,drl)

    #If #DIM eq "3"
     ! u*ux-------------------------------------------------------------------------------
     au = rsxy(i1,i2,i3,0,0)*uu(i1,i2,i3,uc)+rsxy(i1,i2,i3,0,1)*uu(i1,i2,i3,vc) + rsxy(i1,i2,i3,0,2)*uu(i1,i2,i3,wc)
     aum   = rsxy(i1-1,i2,i3,0,0)*(u(i1-1,i2,i3,uc)-gvU)+rsxy(i1-1,i2,i3,0,1)*(u(i1-1,i2,i3,vc)-gvV)+ rsxy(i1-1,i2,i3,0,2)*(u(i1-1,i2,i3,wc)-gvW)
     aup   = rsxy(i1+1,i2,i3,0,0)*(u(i1+1,i2,i3,uc)-gvU)+rsxy(i1+1,i2,i3,0,1)*(u(i1+1,i2,i3,vc)-gvV)+ rsxy(i1+1,i2,i3,0,2)*(u(i1+1,i2,i3,wc)-gvW)

     s1=1
     s2=0
     s3=0
     var = uc
     vard= uc
     
     drl = dr(0)
     getBWENOCoeffs(s1,s2,s3,var,vard,drl)

     !u*vx
     s1=1
     s2=0
     s3=0
     var = vc
     vard= uc
     drl = dr(0)
     
     getBWENOCoeffs(s1,s2,s3,var,vard,drl)

     !u*wx
     s1=1
     s2=0
     s3=0
     var = wc
     vard= uc
     drl = dr(0)
     
     getBWENOCoeffs(s1,s2,s3,var,vard,drl)

     
     !v*uy
     au = rsxy(i1,i2,i3,1,0)*uu(i1,i2,i3,uc)+rsxy(i1,i2,i3,1,1)*uu(i1,i2,i3,vc)+rsxy(i1,i2,i3,1,2)*uu(i1,i2,i3,wc)
     aum   = rsxy(i1,i2-1,i3,1,0)*(u(i1,i2-1,i3,uc)-gvU)+rsxy(i1,i2-1,i3,1,1)*(u(i1,i2-1,i3,vc)-gvV)+rsxy(i1,i2-1,i3,1,2)*(u(i1,i2-1,i3,wc)-gvW)
     aup   = rsxy(i1,i2+1,i3,1,0)*(u(i1,i2+1,i3,uc)-gvU)+rsxy(i1,i2+1,i3,1,1)*(u(i1,i2+1,i3,vc)-gvV)+rsxy(i1,i2+1,i3,1,2)*(u(i1,i2+1,i3,wc)-gvW)

     s1=0
     s2=1
     s3=0
     var = uc
     vard= vc
     drl = dr(1)
     
     getBWENOCoeffs(s1,s2,s3,var,vard,drl)
     !v*vy
     s1=0
     s2=1
     s3=0
     var = vc
     vard= vc
     drl = dr(1)
     
     getBWENOCoeffs(s1,s2,s3,var,vard,drl)

     !v*wy
     s1=0
     s2=1
     s3=0
     var = wc
     vard= vc
     drl = dr(1)
     
     getBWENOCoeffs(s1,s2,s3,var,vard,drl)

     !------------------------
     !w*uz
     au = rsxy(i1,i2,i3,2,0)*uu(i1,i2,i3,uc)+rsxy(i1,i2,i3,2,1)*uu(i1,i2,i3,vc)+rsxy(i1,i2,i3,2,2)*uu(i1,i2,i3,wc)
     aum   = rsxy(i1,i2,i3-1,2,0)*(u(i1,i2,i3-1,uc)-gvU)+rsxy(i1,i2,i3-1,2,1)*(u(i1,i2,i3-1,vc)-gvV)+rsxy(i1,i2,i3-1,2,2)*(u(i1,i2,i3-1,wc)-gvW)
     aup   = rsxy(i1,i2,i3+1,2,0)*(u(i1,i2,i3+1,uc)-gvU)+rsxy(i1,i2,i3+1,2,1)*(u(i1,i2,i3+1,vc)-gvV)+rsxy(i1,i2,i3+1,2,2)*(u(i1,i2,i3+1,wc)-gvW)


     s1=0
     s2=0
     s3=1
     var = uc
     vard= wc
     drl = dr(2) 
     
     getBWENOCoeffs(s1,s2,s3,var,vard,drl)
     
     !w*vz
     s1=0
     s2=0
     s3=1
     var = vc
     vard= wc
     drl = dr(2)
     
     getBWENOCoeffs(s1,s2,s3,var,vard,drl)
     
     !w*wz
     s1=0
     s2=0
     s3=1
     var = wc
     vard= wc
     drl = dr(2) 
     
     getBWENOCoeffs(s1,s2,s3,var,vard,drl)

    #End
     
  else ! if (upwindOrder.eq.4)
     write(*,'(" getBewnoAdvection: only 4th order is avaliable now")' )
     stop 777
     
  end if

  
 #End  

#endMacro 


! --------------------------------------------------------------------------------
!  Macro: getAdvection
! 
! u(i1,i2,i3,.) (input) : current solution
! (i1,i2,i3) (input) : get advection terms solution at this point
! advectionOption (input) : 
! SCALAR: NONE
!         PASSIVE - include equations for a passive scalar
!
! DIM,ORDER,GRIDTYPE (input) :
! UPWIND (input) : CENTERED, UPWIND or BWENO
! 
!  agu(m,n) (output) : m=0,1,nd, n=0,1,nd
!     agu(0,0) : u*ux 
!     agu(1,0) : v*uy
!     agu(2,0) : w*uz 
!
!     agu(0,1) : u*vx 
!     agu(1,1) : v*vy
!     agu(2,1) : w*vz 
!
!     agu(0,2) : u*wx 
!     agu(1,2) : v*wy
!     agu(2,2) : w*wz 
! --------------------------------------------------------------------------------
#beginMacro getAdvection(u,uu,i1,i2,i3,SCALAR,DIM,ORDER,GRIDTYPE,UPWIND, agu)


 #If (#UPWIND == "BWENO")
  ! --- Bweno scheme ---
  ! for testing output this message:
  ! if( t.le. 0. )then
  !    write(*,'(" getAdvection BWENO scheme (7)")') 
  ! end if
  getBwenoAdvection(u,uu,i1,i2,i3,SCALAR,DIM,ORDER,GRIDTYPE, agu)
  
 #Elif (#UPWIND == "UPWIND")
  ! --- upwind scheme ---
  ! for testing output this next message:
  ! if( t.le. 0. )then
  !    write(*,'(" getAdvection upwind scheme (7)")') 
  ! end if
  getUpwindAdvection(u,uu,i1,i2,i3,SCALAR,DIM,ORDER,GRIDTYPE, agu)

 #Elif (#UPWIND == "CENTERED") 
 ! -- centered advection ---
 ! write(*,'(" getAdvection -- centered")')
  #If #DIM == "2"
  ! -- 2D --
  agu(uc,uc)=UU(uc)*UX(uc)
  agu(vc,uc)=UU(vc)*UY(uc)
  
  agu(uc,vc)=UU(uc)*UX(vc)
  agu(vc,vc)=UU(vc)*UY(vc)

  #Else
  ! -- 3D -- *check me*
  agu(uc,uc)=UU(uc)*UX(uc)
  agu(vc,uc)=UU(vc)*UY(uc)
  agu(wc,uc)=UU(wc)*UY(uc)
  
  agu(uc,vc)=UU(uc)*UX(vc)
  agu(vc,vc)=UU(vc)*UY(vc)
  agu(wc,vc)=UU(wc)*UY(vc)
  
  agu(uc,wc)=UU(uc)*UX(wc)
  agu(vc,wc)=UU(vc)*UY(wc)
  agu(wc,wc)=UU(wc)*UY(wc)

  #End

 #Else

   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
   stop 999
   
#End
  
! #If (#UPWIND == "CENTERED") 
!  ! -- centered advection ---
!  ! write(*,'(" getAdvection -- centered")')
!  #If #DIM == "2"
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

! #Elif #UPWIND == "UPWIND" 

!   ! --- upwind scheme ---
!   ! for testing output this next message:
!   if( t.le. 0. )then
!     write(*,'(" getAdvection upwind scheme (7)")') 
!   end if
!   getUpwindAdvection(u,i1,i2,i3,SCALAR,DIM,ORDER,GRIDTYPE, agu)

! #Elif #UPWIND == "BWENO" 

!   ! --- Bweno scheme ---
!   ! for testing output this message:
!   if( t.le. 0. )then
!      write(*,'(" getAdvection BWENO scheme (7)")') 
!   end if
  
!   getBwenoAdvection(u,i1,i2,i3,SCALAR,DIM,ORDER,GRIDTYPE, agu)

! #Else

!   write(*,'(" getAdvection:ERROR: unknown advectionOption.")' )
!   stop 999

! #End

#endMacro
