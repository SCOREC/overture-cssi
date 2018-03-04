! -*- mode: f90; -*-

! **************************************************
! Here are macros that define the:
!      dispersive plane wave solution 
! **************************************************



! *************** Here is the 2D dispersive plane wave solution ******************************

! #defineMacro planeWave2Dex0(x,y,t) sint*dpwc(0)
! #defineMacro planeWave2Dey0(x,y,t) sint*dpwc(1)
! #defineMacro planeWave2Dhz0(x,y,t) sint*dpwc() + cost*dpwc()
! 
! ! one time derivative:
! #defineMacro planeWave2Dext0(x,y,t) (-twoPi*cc)*cos(twoPi*(kx*(x)+ky*(y)-cc*(t)))*pwc(0)
! #defineMacro planeWave2Deyt0(x,y,t) (-twoPi*cc)*cos(twoPi*(kx*(x)+ky*(y)-cc*(t)))*pwc(1)
! #defineMacro planeWave2Dhzt0(x,y,t) (-twoPi*cc)*cos(twoPi*(kx*(x)+ky*(y)-cc*(t)))*pwc(5)
! 
! ! two time derivatives:
! #defineMacro planeWave2Dextt0(x,y,t) (-(twoPi*cc)**2*sin(twoPi*(kx*(x)+ky*(y)-cc*(t)))*pwc(0))
! #defineMacro planeWave2Deytt0(x,y,t) (-(twoPi*cc)**2*sin(twoPi*(kx*(x)+ky*(y)-cc*(t)))*pwc(1))
! #defineMacro planeWave2Dhztt0(x,y,t) (-(twoPi*cc)**2*sin(twoPi*(kx*(x)+ky*(y)-cc*(t)))*pwc(5))
! 
! 
! ! Here are the slow start versions
! #defineMacro planeWave2Dex(x,y,t) (ssf*planeWave2Dex0(x,y,t))
! #defineMacro planeWave2Dey(x,y,t) (ssf*planeWave2Dey0(x,y,t))
! #defineMacro planeWave2Dhz(x,y,t) (ssf*planeWave2Dhz0(x,y,t))
! 
! ! one time derivative:
! #defineMacro planeWave2Dext(x,y,t) (ssf*planeWave2Dext0(x,y,t)+ssft*planeWave2Dex0(x,y,t))
! #defineMacro planeWave2Deyt(x,y,t) (ssf*planeWave2Deyt0(x,y,t)+ssft*planeWave2Dey0(x,y,t))
! #defineMacro planeWave2Dhzt(x,y,t) (ssf*planeWave2Dhzt0(x,y,t)+ssft*planeWave2Dhz0(x,y,t))

! --------------------------------------------------------------------
! Macro: Initialize values needed to eval the dispersive plane wave 
! --------------------------------------------------------------------
#beginMacro initializeDispersivePlaneWave()
  ! --- pre-calculations for the dispersive plane wave ---
  ! kk = twoPi*sqrt( kx*kx+ky*ky+kz*kz)
  ! ck2 = (c*kk)**2

  ! write(*,'(" init-dispersive-plane wave: sr=",e10.2," si=",e10.2)') sr,si

  ! si=-si
  ! s^2 E = -(ck)^2 E - (s^2/eps) P --> gives P = -eps*( 1 + (ck)^2/s^2 ) E 
  sNormSq=sr**2+si**2
  ! sNorm4=sNormSq*sNormSq
  ! pc = -eps*( 2.*sr*si*ck2/sNorm4 )    ! check sign 
  ! ps = -eps*( 1. + ck2*(sr*sr-si*si)/sNorm4 )

  ! Hz = (i/s) * (-1) * (kx*Ey - ky*Ex )/mu
  ! *check me*      
  hfactor = -twoPi*( kx*pwc(1) - ky*pwc(0) )/mu  
  ! hr + i*hi = (i/s)*hfactor
  hr = hfactor*si/sNormSq  
  hi = hfactor*sr/sNormSq  
#endMacro

! =================================================================================
!   Macro: add the correction from the polarization vector to the boundary forcing
! ==================================================================================
#beginMacro addDispersiveBoundaryForcing2d()
  if( getDispersiveBoundaryForcing.eq.1 )then
    ! if( t.le.3*dt )then
    !  write(*,'(" addDispersiveBoundaryForcing2d -- add P to boundary forcing")')
    ! end if
    psum(0)=0.
    psum(1)=0.
    do iv=0,numberOfPolarizationVectors-1
      amp=(psir(iv)*costp-psii(iv)*sintp)*sinxi + (psir(iv)*sintp+psii(iv)*costp)*cosxi
      psum(0) = psum(0) + pwc(0)*amp
      psum(1) = psum(1) + pwc(1)*amp
    end do
    ubv(ex) = ubv(ex) + alphaP*psum(0)
    ubv(ey) = ubv(ey) + alphaP*psum(1)

  end if
#endMacro

! =================================================================================
!   Macro: add the correction from the polarization vector to the boundary forcing
! ==================================================================================
#beginMacro addDispersiveBoundaryForcing3d()
  if( getDispersiveBoundaryForcing.eq.1 )then
    ! if( t.le.5*dt )then
    !   write(*,'(" addDispersiveBoundaryForcing3d -- add P to boundary forcing")')
    ! end if
    psum(0)=0.
    psum(1)=0.
    psum(2)=0.
    do iv=0,numberOfPolarizationVectors-1
      amp=(psir(iv)*costp-psii(iv)*sintp)*sinxi + (psir(iv)*sintp+psii(iv)*costp)*cosxi
      psum(0) = psum(0) + pwc(0)*amp
      psum(1) = psum(1) + pwc(1)*amp
      psum(2) = psum(2) + pwc(2)*amp
    end do
    ubv(ex) = ubv(ex) + alphaP*psum(0)
    ubv(ey) = ubv(ey) + alphaP*psum(1)
    ubv(ez) = ubv(ez) + alphaP*psum(2)

  end if
#endMacro

! --------------------------------------------------------------------
! Macro: Evaluate the dispersive plane wave in 2D
! 
!  x,y,t (input) : point to evaluate at 
!  numberOfTimeDerivatives : evaluate this time derivative
!  ubv(.)  (output) : ubc(ex), etc. 
!  pbv(0:2,0...)  (output) : polarization vectors for dispersive models
! --------------------------------------------------------------------
#beginMacro getDispersivePlaneWave2D(x,y,t,numberOfTimeDerivatives,ubv,pbv)

  xi = twoPi*(kx*(x)+ky*(y))
  sinxi = sin(xi)
  cosxi = cos(xi)
  expt = exp(sr*t) 
  cost=cos(si*t)*expt
  sint=sin(si*t)*expt

  if( numberOfTimeDerivatives==0 )then
    if( polarizationOption.eq.0 )then
      ! amp = cosxi*cost-sinxi*sint *wdh* 2018/01/28 
      ! solution is sin( k*x + si*t)*exp(sr*t) *wdh* 2018/01/28
      ! 
      amp = sinxi*cost+cosxi*sint

      ! For testing we first fill in ubv with the non-dispersive answer, and compare here:
      ! write(*,'(" (i1,i2)=(",i3,",",i3,") sr=",e10.2," si=",e10.2," ubv=",2e12.4," ubv(disp)=",2e12.4)') i1,i2,sr,si,ubv(ex),ubv(ey),pwc(0)*amp,pwc(1)*amp 
      ubv(ex) = pwc(0)*amp 
      ubv(ey) = pwc(1)*amp 

      ! amph = Im( (hr+i*hi)*(cosxi + i sinxi)*( cost + i*sint ) )
      !      = Im( (hr+i*hi)*( [cosxi*cost - sinxi*sint] + i[ cosxi*sint + sinxi*cost] )
      !      = [ cosxi*sint + sinxi*cost]*hr + [cosxi*cost - sinxi*sint]*hi
      !      = [ hr*cost - hi*sint ] *sinxi + [ hr*sint+hi*cost ] *cosxi
      ! amph = (hr*cost-hi*sint)*cosxi - (hr*sint+hi*cost)*sinxi *wdh* 2018/01/28 
      amph = (hr*cost-hi*sint)*sinxi + (hr*sint+hi*cost)*cosxi
      ! write(*,'(" (i1,i2)=(",i3,",",i3,") ubv[Hz]=",e12.4," ubv(disp)[Hz]=",e12.4)') i1,i2,ubv(hz),amph

      ubv(hz) = amph

      do iv=0,numberOfPolarizationVectors-1
        ! amp=(psir(iv)*cost-psii(iv)*sint)*cosxi - (psir(iv)*sint+psii(iv)*cost)*sinxi
        amp=(psir(iv)*cost-psii(iv)*sint)*sinxi + (psir(iv)*sint+psii(iv)*cost)*cosxi
        pbv(0,iv) = pwc(0)*amp 
        pbv(1,iv) = pwc(1)*amp 
      end do 

    else
      ! polarization vector: (ex=pxc, ey=pyc) 
      do iv=0,numberOfPolarizationVectors-1
        pxc = ex + iv*nd
        ! amp=(psir(iv)*cost-psii(iv)*sint)*cosxi - (psir(iv)*sint+psii(iv)*cost)*sinxi
        amp=(psir(iv)*cost-psii(iv)*sint)*sinxi + (psir(iv)*sint+psii(iv)*cost)*cosxi
        ubv(pxc  ) = pwc(0)*amp 
        ubv(pxc+1) = pwc(1)*amp 
      end do 

     ! *check me* -- just repeat hz for now 
     !  ubv(hz) = (hc*cosxi+hs*sinxi)*expt
    end if

  else if( numberOfTimeDerivatives==1 )then
    ! Used for boundary forcing for Hz (Neumann, order 4)

    !write(*,'(" GDPW ntd=1 : fix me")')
    !stop 2738

    ! expt = exp(sr*t)
    ! cost=cos(si*t)*expt
    ! sint=sin(si*t)*expt

    costp=-si*sint+sr*cost  ! (d/dt)( cos(st*t)*exp(sr*t) )
    sintp= si*cost+sr*sint  ! (d/dt)( sin(st*t)*exp(sr*t) )
    if( polarizationOption.eq.0 )then
      ! amp = cosxi*costp-sinxi*sintp   *wdh* 2018/01/28
      amp = sinxi*costp+cosxi*sintp
      ! write(*,'(" (i1,i2)=(",i3,",",i3,") ubv.t=",2e12.4," ubv.t(disp)=",2e12.4)') i1,i2,ubv(ex),ubv(ey),pwc(0)*amp,pwc(1)*amp 

      ubv(ex) = pwc(0)*amp 
      ubv(ey) = pwc(1)*amp 

      ! amph = (hr*costp-hi*sintp)*cosxi - (hr*sintp+hi*costp)*sinxi  *wdh* 2018/01/28
      amph = (hr*costp-hi*sintp)*sinxi + (hr*sintp+hi*costp)*cosxi
      ! write(*,'(" (i1,i2)=(",i3,",",i3,") ubv.t[Hz](nd,d)=",2e12.4," diff=",e12.4)') i1,i2,ubv(hz),amph,ubv(hz)-amph
      ubv(hz) = amph

      ! *wdh* Feb 25, 2018 -- return E_t + alphaP*P_t 
      addDispersiveBoundaryForcing2d()

    else
      ! polarization vector: (ex=pxc, ey=pyc) 
      do iv=0,numberOfPolarizationVectors-1
        pxc = ex + iv*nd
        ! amp=(psir(iv)*costp-psii(iv)*sintp)*cosxi - (psir(iv)*sintp+psii(iv)*costp)*sinxi
        amp=(psir(iv)*costp-psii(iv)*sintp)*sinxi + (psir(iv)*sintp+psii(iv)*costp)*cosxi
        ubv(pxc  ) = pwc(0)*amp 
        ubv(pxc+1) = pwc(1)*amp 
      end do
    end if

  else if( numberOfTimeDerivatives==2 )then
    ! Used for boundary forcing for E (Compatibility, order 4)
   
    cet = -si*sint+sr*cost  ! (d/dt)( cos(st*t)*exp(sr*t) )
    set =  si*cost+sr*sint  ! (d/dt)( sin(st*t)*exp(sr*t) )
    costp=-si*set+sr*cet ! d2/dt2( )
    sintp= si*cet+sr*set ! d2/dt2( )

    ! costp=-si*( si*cost+sr*sint)+sr*(-si*sint+sr*cost) ! d2/dt2( cost)
    ! sintp= si*(-si*sint+sr*cost)+sr*( si*cost+sr*sint) ! d2/dt2
   
    if( polarizationOption.eq.0 )then
      ! amp = cosxi*costp-sinxi*sintp   *wdh* 2018/01/28
      amp = sinxi*costp+cosxi*sintp
      ! write(*,'(" (i1,i2)=(",i3,",",i3,") ubv.t=",2e12.4," ubv.t(disp)=",2e12.4)') i1,i2,ubv(ex),ubv(ey),pwc(0)*amp,pwc(1)*amp

      ubv(ex) = pwc(0)*amp
      ubv(ey) = pwc(1)*amp

      ! amph = (hr*costp-hi*sintp)*cosxi - (hr*sintp+hi*costp)*sinxi  *wdh* 2018/01/28
      amph = (hr*costp-hi*sintp)*sinxi + (hr*sintp+hi*costp)*cosxi
      ! write(*,'(" (i1,i2)=(",i3,",",i3,") ubv.t[Hz](nd,d)=",2e12.4," diff=",e12.4)') i1,i2,ubv(hz),amph,ubv(hz)-amph
      ubv(hz) = amph

      ! *wdh* Feb 25, 2018 -- return E_tt + alphaP*P_tt 
      addDispersiveBoundaryForcing2d()

    else
      ! polarization vector: (ex=pxc, ey=pyc)
      do iv=0,numberOfPolarizationVectors-1
        pxc = ex + iv*nd
        ! amp=(psir(iv)*costp-psii(iv)*sintp)*cosxi - (psir(iv)*sintp+psii(iv)*costp)*sinxi
        amp=(psir(iv)*costp-psii(iv)*sintp)*sinxi + (psir(iv)*sintp+psii(iv)*costp)*cosxi
        ubv(pxc  ) = pwc(0)*amp
        ubv(pxc+1) = pwc(1)*amp
      end do
    end if

  else if( numberOfTimeDerivatives==3 )then
    ! Used for boundary forcing for Hz (compatibility, order 4)

    cet =-si*sint+sr*cost   ! (d/dt)( cos(st*t)*exp(sr*t) )
    set = si*cost+sr*sint   ! (d/dt)( sin(st*t)*exp(sr*t) )
    cett=-si*set+sr*cet     ! d2/dt2( cos(st*t)*exp(sr*t) ) = (d/dt)( cet )
    sett= si*cet+sr*set     ! d2/dt2( sin(st*t)*exp(sr*t) )

    costp=-si*sett+sr*cett  ! d3/dt3( )
    sintp= si*cett+sr*sett  ! d3/dt3( )

    ! costp=-si*( si*(-si*sint+sr*cost)+sr*( si*cost+sr*sint))+sr*(-si*( si*cost+sr*sint)+sr*(-si*sint+sr*cost)) ! d3/dt3( cost)
    ! sintp= si*(-si*( si*cost+sr*sint)+sr*(-si*sint+sr*cost))+sr*( si*(-si*sint+sr*cost)+sr*( si*cost+sr*sint)) ! d3/dt3

    if( polarizationOption.eq.0 )then
      ! amp = cosxi*costp-sinxi*sintp   *wdh* 2018/01/28
      amp = sinxi*costp+cosxi*sintp
      ! write(*,'(" (i1,i2)=(",i3,",",i3,") ubv.t=",2e12.4," ubv.t(disp)=",2e12.4)') i1,i2,ubv(ex),ubv(ey),pwc(0)*amp,pwc(1)*amp
      
      ! write(*,'(" You are here. ")')      

      ubv(ex) = pwc(0)*amp
      ubv(ey) = pwc(1)*amp

      ! amph = (hr*costp-hi*sintp)*cosxi - (hr*sintp+hi*costp)*sinxi  *wdh* 2018/01/28
      amph = (hr*costp-hi*sintp)*sinxi + (hr*sintp+hi*costp)*cosxi
      ! write(*,'(" (i1,i2)=(",i3,",",i3,") ubv.t[Hz](nd,d)=",2e12.4," diff=",e12.4)') i1,i2,ubv(hz),amph,ubv(hz)-amph
      ubv(hz) = amph

      ! *wdh* Feb 25, 2018 -- return E_ttt + alphaP*P_ttt 
      addDispersiveBoundaryForcing2d()

    else
      ! polarization vector: (ex=pxc, ey=pyc)
      do iv=0,numberOfPolarizationVectors-1
        pxc = ex + iv*nd
        ! amp=(psir(iv)*costp-psii(iv)*sintp)*cosxi - (psir(iv)*sintp+psii(iv)*costp)*sinxi
        amp=(psir(iv)*costp-psii(iv)*sintp)*sinxi + (psir(iv)*sintp+psii(iv)*costp)*cosxi
        ubv(pxc  ) = pwc(0)*amp
        ubv(pxc+1) = pwc(1)*amp
      end do
    end if

  else if( numberOfTimeDerivatives==4 )then

    cet =-si*sint+sr*cost   ! (d/dt)( cos(st*t)*exp(sr*t) )
    set = si*cost+sr*sint   ! (d/dt)( sin(st*t)*exp(sr*t) )

    cett=-si*set+sr*cet     ! d2/dt2( cos(st*t)*exp(sr*t) ) = (d/dt)( cet )
    sett= si*cet+sr*set     ! d2/dt2( sin(st*t)*exp(sr*t) )

    cettt=-si*sett+sr*cett  ! d3/dt3( )
    settt= si*cett+sr*sett  ! d3/dt3( )

    costp=-si*settt+sr*cettt  ! d4/dt4( )
    sintp= si*cettt+sr*settt  ! d4/dt4( )


    ! costp=-si*( si*(-si*( si*cost+sr*sint)+sr*(-si*sint+sr*cost))+sr*( si*(-si*sint+sr*cost)+sr*( si*cost+sr*sint)))+sr*(-si*( si*(-si*sint+sr*cost)+sr*( si*cost+sr*sint))+sr*(-si*( si*cost+sr*sint)+sr*(-si*sint+sr*cost)))
    ! sintp= si*(-si*( si*(-si*sint+sr*cost)+sr*( si*cost+sr*sint))+sr*(-si*( si*cost+sr*sint)+sr*(-si*sint+sr*cost)))+sr*( si*(-si*( si*cost+sr*sint)+sr*(-si*sint+sr*cost))+sr*( si*(-si*sint+sr*cost)+sr*( si*cost+sr*sint)))

    if( polarizationOption.eq.0 )then
      ! amp = cosxi*costp-sinxi*sintp   *wdh* 2018/01/28
      amp = sinxi*costp+cosxi*sintp
      ! write(*,'(" (i1,i2)=(",i3,",",i3,") ubv.t=",2e12.4," ubv.t(disp)=",2e12.4)') i1,i2,ubv(ex),ubv(ey),pwc(0)*amp,pwc(1)*amp

      ubv(ex) = pwc(0)*amp
      ubv(ey) = pwc(1)*amp

      ! amph = (hr*costp-hi*sintp)*cosxi - (hr*sintp+hi*costp)*sinxi  *wdh* 2018/01/28
      amph = (hr*costp-hi*sintp)*sinxi + (hr*sintp+hi*costp)*cosxi
      ! write(*,'(" (i1,i2)=(",i3,",",i3,") ubv.t[Hz](nd,d)=",2e12.4," diff=",e12.4)') i1,i2,ubv(hz),amph,ubv(hz)-amph
      ubv(hz) = amph

      ! *wdh* Feb 25, 2018 -- return E_tttt + alphaP*P_tttt 
      addDispersiveBoundaryForcing2d()


    else
      ! polarization vector: (ex=pxc, ey=pyc)
      do iv=0,numberOfPolarizationVectors-1
        pxc = ex + iv*nd
        ! amp=(psir(iv)*costp-psii(iv)*sintp)*cosxi - (psir(iv)*sintp+psii(iv)*costp)*sinxi
        amp=(psir(iv)*costp-psii(iv)*sintp)*sinxi + (psir(iv)*sintp+psii(iv)*costp)*cosxi
        ubv(pxc  ) = pwc(0)*amp
        ubv(pxc+1) = pwc(1)*amp
      end do
    end if

  else
    stop 2738
  end if
#endMacro


! --------------------------------------------------------------------
! Evaluate the dispersive plane wave in 3D
! 
!  x,y,z,t (input) : point to evaluate at 
!  numberOfTimeDerivatives : evaluate this time derivative
!  ubv(.)  (output) : ubc(ex), etc. 
!  pbv(0:2,0...)  (output) : polarization vectors for dispersive models
! --------------------------------------------------------------------
#beginMacro getDispersivePlaneWave3D(x,y,z,t,numberOfTimeDerivatives,ubv,pbv)

  xi = twoPi*(kx*(x)+ky*(y)+kz*(z))
  sinxi = sin(xi)
  cosxi = cos(xi)
  expt = exp(sr*t) 
  cost=cos(si*t)*expt
  sint=sin(si*t)*expt

  if( numberOfTimeDerivatives==0 )then
    if( polarizationOption.eq.0 )then
      ! amp = cosxi*cost-sinxi*sint *wdh* 2018/01/28 
      ! solution is sin( k*x + si*t)*exp(sr*t) *wdh* 2018/01/28
      amp = sinxi*cost+cosxi*sint

      ubv(ex) = pwc(0)*amp 
      ubv(ey) = pwc(1)*amp 
      ubv(ez) = pwc(2)*amp 

      do iv=0,numberOfPolarizationVectors-1
        amp=(psir(iv)*cost-psii(iv)*sint)*sinxi + (psir(iv)*sint+psii(iv)*cost)*cosxi
        pbv(0,iv) = pwc(0)*amp 
        pbv(1,iv) = pwc(1)*amp 
        pbv(2,iv) = pwc(2)*amp 
      end do

    else
      ! polarization vector: (ex=pxc, ey=pyc) 
      do iv=0,numberOfPolarizationVectors-1
        pxc = ex + iv*nd
        ! amp=(psir(iv)*cost-psii(iv)*sint)*cosxi - (psir(iv)*sint+psii(iv)*cost)*sinxi
        amp=(psir(iv)*cost-psii(iv)*sint)*sinxi + (psir(iv)*sint+psii(iv)*cost)*cosxi
        ubv(pxc  ) = pwc(0)*amp 
        ubv(pxc+1) = pwc(1)*amp 
        ubv(pxc+2) = pwc(2)*amp 
      end do 

    end if

  else if( numberOfTimeDerivatives==1 )then

    costp=-si*sint+sr*cost  ! d/dt( cost) 
    sintp= si*cost+sr*sint ! d/dt 
    if( polarizationOption.eq.0 )then
      amp = sinxi*costp+cosxi*sintp
      ubv(ex) = pwc(0)*amp 
      ubv(ey) = pwc(1)*amp 
      ubv(ez) = pwc(2)*amp 

      ! *wdh* Feb 25, 2018 -- return E_t + alphaP*P_t 
      addDispersiveBoundaryForcing3d()

    else
      ! polarization vector: (ex=pxc, ey=pyc) 
      do iv=0,numberOfPolarizationVectors-1
        pxc = ex + iv*nd
        ! amp=(psir(iv)*costp-psii(iv)*sintp)*cosxi - (psir(iv)*sintp+psii(iv)*costp)*sinxi
        amp=(psir(iv)*costp-psii(iv)*sintp)*sinxi + (psir(iv)*sintp+psii(iv)*costp)*cosxi
        ubv(pxc  ) = pwc(0)*amp 
        ubv(pxc+1) = pwc(1)*amp 
        ubv(pxc+2) = pwc(2)*amp 
      end do
    end if

  else if( numberOfTimeDerivatives==2 )then
    ! write(*,'(" GDPW ntd=2 : fix me")')
    ! stop 3738

    ! costp=-si*sint+sr*cost  ! d/dt( cost)
    ! sintp= si*cost+sr*sint ! d/dt
   
    costp=-si*( si*cost+sr*sint)+sr*(-si*sint+sr*cost)  ! d2/dt2( cost)
    sintp= si*(-si*sint+sr*cost)+sr*( si*cost+sr*sint) ! d2/dt2

    if( polarizationOption.eq.0 )then
      ! amp = cosxi*costp-sinxi*sintp   *wdh* 2018/01/28
      amp = sinxi*costp+cosxi*sintp
      ubv(ex) = pwc(0)*amp
      ubv(ey) = pwc(1)*amp
      ubv(ez) = pwc(2)*amp

      ! *wdh* Feb 25, 2018 -- return E_tt + alphaP*P_tt 
      addDispersiveBoundaryForcing3d()


    else
      ! polarization vector: (ex=pxc, ey=pyc)
      do iv=0,numberOfPolarizationVectors-1
        pxc = ex + iv*nd
        ! amp=(psir(iv)*costp-psii(iv)*sintp)*cosxi - (psir(iv)*sintp+psii(iv)*costp)*sinxi
        amp=(psir(iv)*costp-psii(iv)*sintp)*sinxi + (psir(iv)*sintp+psii(iv)*costp)*cosxi
        ubv(pxc  ) = pwc(0)*amp
        ubv(pxc+1) = pwc(1)*amp
        ubv(pxc+2) = pwc(2)*amp
      end do
    end if


  else if( numberOfTimeDerivatives==3 )then
    ! write(*,'(" GDPW ntd=3 : fix me")')
    ! stop 3738

    ! costp=-si*sint+sr*cost  ! d/dt( cost)
    ! sintp= si*cost+sr*sint ! d/dt

    ! costp=-si*( si*cost+sr*sint)+sr*(-si*sint+sr*cost)  ! d2/dt2( cost)
    ! sintp= si*(-si*sint+sr*cost)+sr*( si*cost+sr*sint) ! d2/dt2

    costp=-si*( si*(-si*sint+sr*cost)+sr*( si*cost+sr*sint))+sr*(-si*( si*cost+sr*sint)+sr*(-si*sint+sr*cost))  ! d3/dt3( cost)
    sintp= si*(-si*( si*cost+sr*sint)+sr*(-si*sint+sr*cost))+sr*( si*(-si*sint+sr*cost)+sr*( si*cost+sr*sint)) ! d3/dt3

    if( polarizationOption.eq.0 )then
      ! amp = cosxi*costp-sinxi*sintp   *wdh* 2018/01/28
      amp = sinxi*costp+cosxi*sintp
      ubv(ex) = pwc(0)*amp
      ubv(ey) = pwc(1)*amp
      ubv(ez) = pwc(2)*amp

      ! *wdh* Feb 25, 2018 -- return E_ttt + alphaP*P_ttt 
      addDispersiveBoundaryForcing3d()


    else
      ! polarization vector: (ex=pxc, ey=pyc)
      do iv=0,numberOfPolarizationVectors-1
        pxc = ex + iv*nd
        ! amp=(psir(iv)*costp-psii(iv)*sintp)*cosxi - (psir(iv)*sintp+psii(iv)*costp)*sinxi
        amp=(psir(iv)*costp-psii(iv)*sintp)*sinxi + (psir(iv)*sintp+psii(iv)*costp)*cosxi
        ubv(pxc  ) = pwc(0)*amp
        ubv(pxc+1) = pwc(1)*amp
        ubv(pxc+2) = pwc(2)*amp
      end do
    end if


  else if( numberOfTimeDerivatives==4 )then
    write(*,'(" GDPW ntd=4 : fix me")')
    stop 3738
  else
    stop 3738
  end if

#endMacro



