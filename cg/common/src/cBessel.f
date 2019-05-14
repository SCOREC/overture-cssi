!
!  Evaluate the bessel function J_nu(z) with complex argument z=(zr,ri)
!  Return J = (jr,ji)
!
      subroutine cBesselJ( nu, zr,zi, jr,ji )

      double precision nu, zr,zi,jr,ji

      double precision cjr(10), cji(10)
      integer np,nz,ierr,kode

      kode=1 ! do not scale result by exp(-abs(zi))
      np=1 ! we just want J_nu 
      call zbesj( zr,zi,nu,kode,np,cjr,cji,nz,ierr)
      if( nz.ne.0 .or. ierr.ne.0 )then
        write(*,'("WARNING: zbesj: nz,ierr=",2i4)') nz,ierr
        write(*,'("         nu,zr,zi,jr,ji=",5e12.4)') 
     +       nu,zr,zi,cjr(1),cji(1)
      end if

      jr=cjr(1)
      ji=cji(1)

      return 
      end 

!
!  Evaluate the bessel function Y_nu(z) with complex argument z=(zr,ri)
!  Return Y = (yr,yi)
!
      subroutine cBesselY( nu, zr,zi, yr,yi )

      double precision nu, zr,zi,yr,yi

      double precision cyr(10), cyi(10),cwkr(10),cwrki(10)
      integer np,nz,ierr,kode

      kode=1 ! do not scale result by exp(-abs(zi))
      np=1 ! we just want Y_nu 
      call zbesy( zr,zi,nu,kode,np,cyr,cyi,nz,
     &   cwkr,cwrki,ierr)
      if( nz.ne.0 .or. ierr.ne.0 )then
        write(*,'("WARNING: zbesy: nz,ierr=",2i4)') nz,ierr
      end if
      !write(*,'(" zbesy: nz,ierr=",2i4)') nz,ierr
      !write(*,'(" zbesy: Y=",2e12.4)') cyr(1),cyi(1)

      yr=cyr(1)
      yi=cyi(1)

      return 
      end 

!
!  Evaluate the real and imaginary parts of the solution to
!  the oscillating buble.
!    v_r    (r,theta,t) = \tilde{v}_r    (r) exp(i (n theta - omega t))
!    v_theta(r,theta,t) = \tilde{v}_theta(r) exp(i (n theta - omega t))
!    p      (r,theta,t) = \tilde{p}      (r) exp(i (n theta - omega t))
!  
!  Inputs:
!    r:   radius
!    Rb:   radius of boundary
!    n:   number of periods
!    mu:  viscosity
!    lr:  real part of lambda
!    li:  imag part of lambda
!    cr:  real part of integration constant
!    ci:  imag part of integration constant
!    p0r: real part of particular soln constant
!    p0i: imag part of particular soln constant
!  Outputs
!    vrr: real part of \tilde{v}_r
!    vri: imag part of \tilde{v}_r
!    vtr: real part of \tilde{v}_theta
!    vti: imag part of \tilde{v}_theta
!    pr:  real part of \tilde{p}
!    pi:  imag part of \tilde{p}
!

      subroutine evalOscillatingBubble(r,Rb, n, mu, omegar,omegai,
     +     cr,ci,dr,di,vrr,vri,vtr,vti,pr,pi)

      implicit none
      double precision r,Rb,mu,
     +     lr,li,cr,ci,dr,di,omegar,omegai,
     +     vrr,vri,vtr,vti,pr,pi,
     +     jnr,jni,
     +     jnpr,jnpi,jnmr,jnmi,
     +     fnu,zr,zi
      integer n
      complex*16 vr,vt,p,  ! final solution 
     +     l,c,d,nc,omega, ! complex versions of inputs
     +     rl,rc,Rbc,muc,  ! 
     +     jn,jnp,jnm,     ! bessel functions
     +     f,I             ! dummy variable


c     get complex vars
c$$$      l = dcmplx(lr,li)
      c = dcmplx(cr,ci)
      d = dcmplx(dr,di)
      I = dcmplx(0.0d0,1.0d0)
      omega = dcmplx(omegar,omegai)

      rc  = dcmplx(r ,0.0d0) ! complex r
      Rbc = dcmplx(Rb,0.0d0) ! complex Rb
      muc = dcmplx(mu,0.0d0) ! complex mu
      nc  = dcmplx(dble(n),0.0d0) ! complex n

      l = zsqrt(I*omega/muc)
      rl  = l*rc             ! r * lambda
      f = (d * nc) / (muc * (Rbc ** n) * (l ** 2))

c     check if r is close to 0 or not
      if (abs(r) > 1.0d-14) then
c       standard computations when r != 0

c       evaluate bessel funcs
        call cBesselJ( dble(n  ),dreal(rl),dimag(rl),jnr ,jni )
        call cBesselJ( dble(n+1),dreal(rl),dimag(rl),jnpr,jnpi)

        jn  = dcmplx(jnr ,jni )
        jnp = dcmplx(jnpr,jnpi)

c       evaluate vr, vt, and p
c       radial component of velocity
        vr = (c * jn / rc) + (f * (rc ** (n-1)))

c       circumferential component of velocity
        vt = (c * I * jn / rc)
     +       -(c * I * l * jnp / nc) 
     +       +(I * f * (rc ** (n-1)))


      else
c       r is close to 0, so we will evaluate Jn / r 
c       using the limit definition
        
c       evaluate bessel funcs
        call cBesselJ( dble(n-1),dreal(rl),dimag(rl),jnmr,jnmi)
        call cBesselJ( dble(n+1),dreal(rl),dimag(rl),jnpr,jnpi)
        
        jnm = dcmplx(jnmr,jnmi)
        jnp = dcmplx(jnpr,jnpi)

c       evaluate vr, vt, and p
c       radial component of velocity
        vr = (c * l * (jnm + jnp) / (dcmplx(2.0d0,0.0d0)))
     +       +(f * (rc ** (n-1)))

c       circumferential component of velocity
        vt = (c * I * l * (jnm + jnp) / (dcmplx(2.0d0,0.0d0)))
     +       -(c * I * l * jnp / nc) 
     +       +(I * f * (rc ** (n-1)))

c$$$        write(*,'("r=0, jnmr,jnmi=",e12.3,e12.3)') jnmr,jnmi
c$$$        write(*,'("r=0, jnpr,jnpi=",e12.3,e12.3)') jnpr,jnpi
c$$$        write(*,'("r=0, vr,vt=",e12.3,e12.3)') dreal(vr),dreal(vt)
      end if
      

c       pressure
        p = d * (rc ** n) / (Rbc ** n)

c       assign vr, vt, and p
        vrr = dreal(vr)
        vri = dimag(vr)
        vtr = dreal(vt)
        vti = dimag(vt)
        pr  = dreal( p)
        pi  = dimag( p)

      return
      end

!
!  Evaluate the real and imaginary parts of the solution to
!  the cylindrical stream
!    v_r    (r,theta,z,t) = \tilde{v}_r    (r) exp(i (n theta + kz z - omega t))
!    v_theta(r,theta,z,t) = \tilde{v}_theta(r) exp(i (n theta + kz z - omega t))
!    v_z    (r,theta,z,t) = \tilde{v}_z    (r) exp(i (n theta + kz z - omega t))
!    p      (r,theta,z,t) = \tilde{p}      (r) exp(i (n theta + kz z - omega t))
!  
!  Inputs:
!    r:   radius
!    kz:
!    n:   number of periods
!    mu:  viscosity
!    A,B,C: integration constants
!    alpha,lambda: wave numbers
!  Outputs
!    vrr: real part of \tilde{v}_r
!    vri: imag part of \tilde{v}_r
!    vtr: real part of \tilde{v}_theta
!    vti: imag part of \tilde{v}_theta
!    vzr: real part of \tilde{v}_z
!    vzi: imag part of \tilde{v}_z
!    pr:  real part of \tilde{p}
!    pi:  imag part of \tilde{p}
!
      subroutine evalCylindricalStream(
     +     r,kz,n,mu,
     +     Ar,Ai,Br,Bi,Cr,Ci,
     +     alphar,alphai,lambdar,lambdai,
     +     vrr,vri,vtr,vti,vzr,vzi,pr,pi)

      implicit none
      double precision r,kz,n,mu,
     +     Ar,Ai,Br,Bi,Cr,Ci,
     +     alphar,alphai,lambdar,lambdai,
     +     vrr,vri,vtr,vti,vzr,vzi,pr,pi,
     +     jnikr,jnpikr,jnar,jnpar,
     +     jniki,jnpiki,jnai,jnpai
      complex*16 vr,vt,vz,p,          ! final solution
     +     A,B,C,                     ! integration constants
     +     alpha,lambda,              ! wave numbers
     +     I,jnik,jnpik,jna,jnpa,     ! dummy variables
     +     alr,ikr

      A=dcmplx(Ar,Ai)
      B=dcmplx(Br,Bi)
      C=dcmplx(Cr,Ci)
      alpha=dcmplx(alphar,alphai)
      lambda=dcmplx(lambdar,lambdai)
      I=dcmplx(0.0d0,1.0d0)

      ikr=I*kz*r
      alr=alpha*r

c     check if r is close to 0 or not
      if (abs(r) > 1.0d-14) then
c     standard computations when r != 0
        
c     evaluate bessel funcs
        call cBesselJ(dble(n+0),dreal(ikr),dimag(ikr),jnikr,jniki)
        call cBesselJ(dble(n+1),dreal(ikr),dimag(ikr),jnpikr,jnpiki)
        call cBesselJ(dble(n+0),dreal(alr),dimag(alr),jnar,jnai)
        call cBesselJ(dble(n+1),dreal(alr),dimag(alr),jnpar,jnpai)

        jnik =dcmplx(jnikr ,jniki )
        jnpik=dcmplx(jnpikr,jnpiki)
        jna  =dcmplx(jnar  ,jnai  )
        jnpa =dcmplx(jnpar ,jnpai )
        
c     evaluate vr, vt, vz, and p
        p=   +A*jnik
        vz=  +B*jna
     +       +I*A*kz*jnik/(mu*lambda**2)
        vr=  +C*jna/r
     +       +I*B*kz*(-jnpa+n*jna/(alpha*r))/alpha
     +       +I*A*kz*(-jnpik-I*n*jnik/(kz*r))/(mu*lambda**2)
        vt=  +(I/(r*mu*n*(alpha*lambda)**2))*(
     +       +I*B*jna*kz*mu*(lambda*n)**2
     +       -C*jnpa*mu*r*(alpha**3)*(lambda**2)
     +       +C*jna*mu*n*(alpha*lambda)**2
     +       +A*jnik*(alpha*n)**2)
        

      else
c       r is close to 0, so we will evaluate Jn / r 
c       using the limit definition
        p =dcmplx(dble(0),dble(0))
        vz=dcmplx(dble(0),dble(0))
        vr=dcmplx(dble(0),dble(0))
        vt=dcmplx(dble(0),dble(0))
      end if
      
c     assign vr, vt, and p
      vrr = dreal(vr)
      vri = dimag(vr)
      vtr = dreal(vt)
      vti = dimag(vt)
      vzr = dreal(vz)
      vzi = dimag(vz)
      pr  = dreal( p)
      pi  = dimag( p)

      return
      end


      subroutine evalCapillaryFlow(k,y,mu,alphar,alphai,ar,ai,br,bi,
     +     cr,ci,dr,di,uhr,uhi,vhr,vhi,phr,phi)
      implicit none
      double precision k,y,mu,alphar,alphai,ar,ai,br,bi,
     +     cr,ci,dr,di,uhr,uhi,vhr,vhi,phr,phi
      complex*16 uh,vh,ph,
     +     alpha,a,b,c,d,
     +     kc,yc,muc,I

c     get complex vars
      alpha = dcmplx(alphar,alphai)
      a     = dcmplx(ar    ,ai    )
      b     = dcmplx(br    ,bi    )
      c     = dcmplx(cr    ,ci    )
      d     = dcmplx(dr    ,di    )
      I = dcmplx(0.0d0,1.0d0)

      kc  = dcmplx(k ,0.0d0)
      yc  = dcmplx(y ,0.0d0)
      muc = dcmplx(mu,0.0d0)

c     get hat vars
      uh = -c*alpha/(I*kc)*exp( alpha*yc)
     +     +d*alpha/(I*kc)*exp(-alpha*yc)
     +     -kc*a/(I*((k**2)-(alpha**2))*muc)*exp( k*y)
     +     -kc*b/(I*((k**2)-(alpha**2))*muc)*exp(-k*y)

      vh =  c*exp( alpha*yc)
     +     +d*exp(-alpha*yc)
     +     +kc*a/(((k**2)-(alpha**2))*muc)*exp( k*y)
     +     -kc*b/(((k**2)-(alpha**2))*muc)*exp(-k*y)

      ph =  a*exp( kc*yc)
     +     +b*exp(-kc*yc)

c     take real and imaginary parts of solution
      uhr = dreal(uh)
      uhi = dimag(uh)
      vhr = dreal(vh)
      vhi = dimag(vh)
      phr = dreal(ph)
      phi = dimag(ph)

      
      return 
      end

!
!
! Evaluate the solution for shear flow (FSI)
!
! u1  = amp    ( A cos(ks y) + B sin(ks y)) exp(i omega t)
! u1y = amp ks (-A sin(ks y) + B cos(ks y)) exp(i omega t)
!
! Return:
!  ur  = real( A cos(ks y) + B sin(ks y) )
!  ui  = imag( A cos(ks y) + B sin(ks y) )
!  uyr = real( ks (-A sin(ks y) + B cos(ks y)) )
!  uyi = imag( ks (-A sin(ks y) + B cos(ks y)) )
!
      subroutine evalFibShearSolid(ksr,ksi,ar,ai,br,bi,y,ur,ui,uyr,uyi)
      implicit none
      double precision ksr,ksi,ar,ai,br,bi,y,ur,ui,uyr,uyi
      complex*16 ks,a,b,u,uy,I,cosy,siny
      
c     get complex vars
      ks = dcmplx(ksr,ksi)
      a  = dcmplx(ar ,ai )
      b  = dcmplx(br ,bi )
      I  = dcmplx(0.0d0,1.0d0)

c     get solution
      cosy = 0.5d0*(exp(I*ks*y)+exp(-I*ks*y))
      siny = 0.5d0*(exp(I*ks*y)-exp(-I*ks*y))

      u  =    ( a*cos(ks*y)+b*sin(ks*y))
      uy = ks*(-a*sin(ks*y)+b*cos(ks*y))

c     take real and imaginary parts of the solution
      ur  = dreal(u)
      ui  = dimag(u)
      uyr = dreal(uy)
      uyi = dimag(uy)

      return
      end


      subroutine evalFibShearSolidFull(ksr,ksi,ar,ai,br,bi,y,t,
     +     ur,ui,vr,vi,uyr,uyi,omegar,omegai)
      implicit none
      double precision ksr,ksi,ar,ai,br,bi,y,ur,ui,uyr,uyi,
     +     vr,vi,t,omegar,omegai

      complex*16 ks,a,b,u,uy,I,cosy,siny,omega,v
      
c     get complex vars
      ks = dcmplx(ksr,ksi)
      a  = dcmplx(ar ,ai )
      b  = dcmplx(br ,bi )
      I  = dcmplx(0.0d0,1.0d0)
      omega = dcmplx(omegar,omegai)

c     get solution
      cosy = 0.5d0*(exp(I*ks*y)+exp(-I*ks*y))
      siny = 0.5d0*(exp(I*ks*y)-exp(-I*ks*y))

      u  =          ( a*cos(ks*y)+b*sin(ks*y))*exp(I*omega*t)
      v  =  I*omega*( a*cos(ks*y)+b*sin(ks*y))*exp(I*omega*t)
      uy =       ks*(-a*sin(ks*y)+b*cos(ks*y))*exp(I*omega*t)

      

c     take real and imaginary parts of the solution
      ur  = dreal(u)
      ui  = dimag(u)
      vr  = dreal(v)
      vi  = dimag(v)
      uyr = dreal(uy)
      uyi = dimag(uy)

c$$$      write(*,'("ks  = ",2f10.3)')ks
c$$$      write(*,'("a  = ",2f10.3)')a
c$$$      write(*,'("b  = ",2f10.3)')b
c$$$      write(*,'("omega  = ",2f10.3)')omega
c$$$      write(*,'("(y,t)  = ",2f10.3)')y,t
c$$$      write(*,'("v  = ",2f10.3)')vr,vi
c$$$      write(*,'("uy = ",2f10.3)')uyr,uyi

      return
      end


!
! v1  = amp i omega ( C exp(kf y) + D exp(-kf y)) exp(i omega t)
!
! Return:
!  vr = real(i omega ( C exp(kf y) + D exp(-kf y)) exp(i omega t))
!  vi = imag(i omega ( C exp(kf y) + D exp(-kf y)) exp(i omega t))
!
      subroutine evalFibShearFluid(kfr,kfi,omegar,omegai,cr,ci,dr,di,y,
     +     vr,vi)
      implicit none
      double precision kfr,kfi,cr,ci,dr,di,y,vr,vi,omegar,omegai
      complex*16 kf,c,d,v,I,omega
      
c     get complex vars
      kf = dcmplx(kfr,kfi)
      omega = dcmplx(omegar,omegai)
      c  = dcmplx(cr ,ci )
      d  = dcmplx(dr ,di )
      I  = dcmplx(0.0d0,1.0d0)

c     get solution
      v = I*omega*(c*exp(kf*y)+d*exp(-kf*y))
      
c     take real and imaginary parts
      vr = dreal(v)
      vi = dimag(v)


      return
      end


      subroutine evalRadialFibShearFluid(kfr,kfi,omegar,omegai,
     +     ar,ai,br,bi,rho,r,r0,t,vt,p)
      implicit none
      double precision kfr,kfi,ar,ai,br,bi,r,r0,t,omegar,omegai,
     +     rho,vt,p,j1r,j1i,y1r,y1i,j0r,j0i,y0r,y0i,
     +     j1r0,j1i0,y1r0,y1i0,j0r0,j0i0,y0r0,y0i0
      complex*16 kf,rk,rk0,a,b,vc,pc,I,omega,
     +     j1,y1,j0,y0,j10,y10,j00,y00,
     +     term1,term2,term3,term4,term5,term6,term7
      
c     get complex vars
      kf = dcmplx(kfr,kfi)
      omega = dcmplx(omegar,omegai)
      a  = dcmplx(ar ,ai )
      b  = dcmplx(br ,bi )
      I  = dcmplx(0.0d0,1.0d0)
      
      rk  = r*kf
      rk0 = r0*kf

c     evaluate bessel funcs
      call cBesselJ( dble(1),dreal(rk),dimag(rk),j1r,j1i)
      call cBesselY( dble(1),dreal(rk),dimag(rk),y1r,y1i)
      call cBesselJ( dble(0),dreal(rk),dimag(rk),j0r,j0i)
      call cBesselY( dble(0),dreal(rk),dimag(rk),y0r,y0i)
      call cBesselJ( dble(1),dreal(rk0),dimag(rk0),j1r0,j1i0)
      call cBesselY( dble(1),dreal(rk0),dimag(rk0),y1r0,y1i0)
      call cBesselJ( dble(0),dreal(rk0),dimag(rk0),j0r0,j0i0)
      call cBesselY( dble(0),dreal(rk0),dimag(rk0),y0r0,y0i0)

      j1 = dcmplx(j1r,j1i)
      y1 = dcmplx(y1r,y1i)
      j0 = dcmplx(j0r,j0i)
      y0 = dcmplx(y0r,y0i)
      j10 = dcmplx(j1r0,j1i0)
      y10 = dcmplx(y1r0,y1i0)
      j00 = dcmplx(j0r0,j0i0)
      y00 = dcmplx(y0r0,y0i0)

c     calculate terms in pressure
      term1 = dble(.5)*(j00**2 + j10**2 - j0**2 - j1**2)
      term2 = dble(.5)*(j00*y00 + j10*y10 - j0*y0 - j1*y1)
      term3 = dble(.5)*(y00**2 + y10**2 - y0**2 - y1**2)
      term4 = dble(.5)*(j00*conjg(j00) + j10*conjg(j10) 
     +     - j0*conjg(j0) - j1*conjg(j1))
      term5 = dble(.5)*(j00*conjg(y00) + j10*conjg(y10) 
     +     - j0*conjg(y0) - j1*conjg(y1))
      term6 = dble(.5)*(y00*conjg(j00) + y10*conjg(j10) 
     +     - y0*conjg(j0) - y1*conjg(j1))
      term7 = dble(.5)*(y00*conjg(y00) + y10*conjg(y10) 
     +     - y0*conjg(y0) - y1*conjg(y1))


c     get solution
      vc = (a*j1+b*y1)*exp(I*omega*t)
c$$$      pc = dble(1)*rho*((a**2)*term1+dble(2)*a*b*term2+(b**2)*term3 
c$$$     +     +a*conjg(a)*term4+a*conjg(b)*term5
c$$$     +     +b*conjg(a)*term6+b*conjg(b)*term7) / dble(2)
      pc = dble(1)*rho*((a**2)*term1+dble(2)*a*b*term2+(b**2)*term3)

c     take real and imaginary parts
      vt = dreal(vc)
      p  = dreal(pc)
      


      return
      end


      subroutine evalRadialFibShearSolid(ksr,ksi,omegar,omegai,
     +     cr,ci,r,t,ut,vt,at,utr,vtr)
      implicit none
      double precision ksr,ksi,cr,ci,r,t,omegar,omegai,
     +     ut,vt,at,utr,vtr,
     +     j1r,j1i,y1r,y1i,j0r,j0i,y0r,y0i
      complex*16 ks,rk,rk0,c,uc,vc,ac,ucr,vcr,I,omega,
     +     j1,y1,j0,y0

      
c     get complex vars
      ks = dcmplx(ksr,ksi)
      omega = dcmplx(omegar,omegai)
      c  = dcmplx(cr ,ci )
      I  = dcmplx(0.0d0,1.0d0)
      
      rk  = r*ks

c     evaluate bessel funcs
      call cBesselJ( dble(1),dreal(rk),dimag(rk),j1r,j1i)
      call cBesselJ( dble(0),dreal(rk),dimag(rk),j0r,j0i)

      j1 = dcmplx(j1r,j1i)
      j0 = dcmplx(j0r,j0i)

c     get solution
      uc  =         c*j1*exp(I*omega*t)
      vc  = I*omega*uc
      ac  = -(omega**2)*uc

c     calculate the stress, srt, normalized by mubar
      if (abs(r) > 1.0d-14) then
        ucr = c*(ks*j0-2.0d0*j1/r)*exp(I*omega*t)
      else
        ucr = c*(ks*j0-ks)*exp(I*omega*t)
      end if

      vcr = I*omega*ucr

c$$$      write(*,'("exp(i omega t) = ",1f10.3)')abs(exp(I*omega*t))
c$$$      write(*,'("vchat = ",2f10.3)')vc/exp(I*omega*t)

c     take real and imaginary parts
      ut  = dreal(uc)
      vt  = dreal(vc)
      at  = dreal(ac)
      utr = dreal(ucr)
      vtr = dreal(vcr)

      return
      end




      subroutine evalFibCartWaveFluid(omegar,omegai,k,mu,rho,mubar,
     +     lambdabar,Ar,Ai,Br,Bi,epsilon,
     +     x,y,t,H,
     +     v1r,v2r,pr)
      implicit none
      double precision omegar,omegai,k,mu,rho,mubar,
     +     lambdabar,Ar,Ai,Br,Bi,epsilon,
     +     x,y,t,H,v1r,v2r,pr
      complex*16 omega,A,B,A1,B1,A2,B2,r,v1,v2,p,i

!     need 
!     mubar, lambdabar,
!     omega, k, mu, rho, A, B
!     epsilon

c     get complex vars
      omega = dcmplx(omegar,omegai)
      A     = dcmplx(Ar    ,Ai    )
      B     = dcmplx(Br    ,Bi    )
      i  = dcmplx(0.0d0,1.0d0)

c     calculate quantities
      r= i / mu * sqrt(i * mu * (i * k ** 2 * mu + omega * rho))

      B1=-i*B*k**2*(lambdabar+2*mubar)/(rho*omega**2);
      A1=-i*A*k**2*(lambdabar+2*mubar)/(rho*omega**2);
      B2=-A*k**2*(lambdabar+2*mubar)/(rho*omega**2);
      A2=-B*k**2*(lambdabar+2*mubar)/(rho*omega**2);

c     get solution
      v1 = i * epsilon * omega * ((A1 * sinh(k * (y - H))) + 
     +     (B1 * cosh(k * (y - H))) - i * r * B2 / k * sinh(r * 
     +     (y - H)) - (B1 * cosh(r * (y - H)))) * exp(i * (k * x - 
     +     omega * t));
      v2 = i * epsilon * omega * ((A2 * sinh(k * (y - H))) + 
     +     (B2 * cosh(k * (y - H))) + i * k / r * B1 * sinh(r * 
     +     (y - H)) - (B2 * cosh(r * (y - H)))) * exp(i * (k * x - 
     +     omega * t));
      p  = (A * sinh(k * (y - H)) + B * cosh(k * (y - H))) * k * epsilon 
     +     * (lambdabar + 2 * mubar) * exp(i * (k * x - omega * t));

c     take real part
      v1r = dreal(v1)
      v2r = dreal(v2)
      pr  = dreal(p)


      return
      end

      subroutine evalFibCartWaveSolid(omegar,omegai,k,mubar,rhobar,
     +     lambdabar,k1r,k1i,k2r,k2i,epsilon,
     +     x,y,t,Hbar, 
     +     u1barr,u2barr,v1barr,v2barr,s11barr,s12barr,s22barr)
      implicit none
      double precision omegar,omegai,k,mubar,rhobar,
     +     lambdabar,k1r,k1i,k2r,k2i,epsilon,cs,cp,
     +     x,y,t,Hbar,
     +     u1barr,u2barr,v1barr,v2barr,s11barr,s12barr,s22barr
      complex*16 omega,k1,k2,k3,k4,alpha1,alpha2,alpha3,alpha4,
     +     a1,a2,a3,a4,u1bar,u2bar,v1bar,v2bar,s11bar,s12bar,s22bar,i
!     need 
!     mubar, rhobar, lambdabar,
!     omega, k
!     k1,k2

c     get complex vars
      omega = dcmplx(omegar,omegai)
      k1    = dcmplx(k1r   ,k1i   )
      k2    = dcmplx(k2r   ,k2i   )
      i  = dcmplx(0.0d0,1.0d0)

c     calculate quantities
      cs = sqrt(mubar/rhobar)
      cp = sqrt((lambdabar+2*mubar)/rhobar)

      alpha1=sqrt(k**2-omega**2/cs**2)
      alpha2=sqrt(k**2-omega**2/cp**2)
      alpha3=-sqrt(k**2-omega**2/cs**2)
      alpha4=-sqrt(k**2-omega**2/cp**2)

      a1= i*((-2*mubar-lambdabar)*k**2+alpha1**2*mubar+omega**2*rhobar)
     +     /(k*(mubar+lambdabar)*alpha1)
      a2= i*((-2*mubar-lambdabar)*k**2+alpha2**2*mubar+omega**2*rhobar)
     +     /(k*(mubar+lambdabar)*alpha2)
      a3= i*((-2*mubar-lambdabar)*k**2+alpha3**2*mubar+omega**2*rhobar)
     +     /(k*(mubar+lambdabar)*alpha3)
      a4= i*((-2*mubar-lambdabar)*k**2+alpha4**2*mubar+omega**2*rhobar)
     +     /(k*(mubar+lambdabar)*alpha4)

      k3 = (a1*k1+a2*k1+2*a2*k2)/(a1-a2)
      k4 = -(2*a1*k1+a1*k2+a2*k2)/(a1-a2)

c     get solution
      u1bar = epsilon * (k1 * exp(alpha1 * (y + Hbar)) + k2 * exp(alpha2 
     +     * (y + Hbar)) - (a1 * k1 + a2 * k2 - a4 * k1 - a4 * k2) / (a3 
     +     - a4) * exp(-alpha1 * (y + Hbar)) + (a1 * k1 + a2 * k2 - a3 * 
     +     k1 - a3 * k2) / (a3 - a4) * exp(-alpha2 * (y + Hbar))) * 
     +     exp(i * (k * x - omega * t))
      u2bar = epsilon * (a1 * k1 * exp(alpha1 * (y + Hbar)) + a2 * k2 * 
     +     exp(alpha2 * (y + Hbar)) - a3 * (a1 * k1 + a2 * k2 - a4 * k1 
     +     - a4 * k2) / (a3 - a4) * exp(-alpha1 * (y + Hbar)) + a4 * (a1 
     +     * k1 + a2 * k2 - a3 * k1 - a3 * k2) / (a3 - a4) * exp(-alpha2 
     +     * (y + Hbar))) * exp(i * (k * x - omega * t));
      v1bar = -i * epsilon * (k1 * exp(alpha1 * (y + Hbar)) + k2 * 
     +     exp(alpha2 * (y + Hbar)) - (a1 * k1 + a2 * k2 - a4 * k1 - 
     +     a4 * k2) / (a3 - a4) * exp(-alpha1 * (y + Hbar)) + (a1 * k1 + 
     +     a2 * k2 - a3 * k1 - a3 * k2) / (a3 - a4) * exp(-alpha2 * (y + 
     +     Hbar))) * omega * exp(i * (k * x - omega * t))
      v2bar = -i * epsilon * (a1 * k1 * exp(alpha1 * (y + Hbar)) + a2 * 
     +     k2 * exp(alpha2 * (y + Hbar)) - a3 * (a1 * k1 + a2 * k2 - a4
     +     * k1 - a4 * k2) / (a3 - a4) * exp(-alpha1 * (y + Hbar)) + a4 
     +     * (a1 * k1 + a2 * k2 - a3 * k1 - a3 * k2) / (a3 - a4) * 
     +     exp(-alpha2 * (y + Hbar))) * omega * exp(i * (k * x - omega 
     +     * t))
      s11bar = i * (lambdabar + 2 * mubar) * epsilon * (k1 * exp(alpha1 
     +     * (y + Hbar)) + k2 * exp(alpha2 * (y + Hbar)) - (a1 * k1 + a2 
     +     * k2 - a4 * k1 - a4 * k2) / (a3 - a4) * exp(-alpha1 * (y + 
     +     Hbar)) + (a1 * k1 + a2 * k2 - a3 * k1 - a3 * k2) / (a3 - a4) 
     +     * exp(-alpha2 * (y + Hbar))) * k * exp(i * (k * x - omega * t
     +     )) + lambdabar * epsilon * (a1 * k1 * alpha1 * exp(alpha1 * 
     +     (y + Hbar)) + a2 * k2 * alpha2 * exp(alpha2 * (y + Hbar)) + 
     +     a3 * (a1 * k1 + a2 * k2 - a4 * k1 - a4 * k2) / (a3 - a4) * 
     +     alpha1 * exp(-alpha1 * (y + Hbar)) - a4 * (a1 * k1 + a2 * k2
     +     - a3 * k1 - a3 * k2) / (a3 - a4) * alpha2 * exp(-alpha2 * 
     +     (y + Hbar))) * exp(i * (k * x - omega * t))

      s12bar = mubar * (epsilon * (k1 * alpha1 * exp(alpha1 * (y + 
     +     Hbar)) + k2 * alpha2 * exp(alpha2 * (y + Hbar)) + (a1 * k1 
     +     + a2 * k2 - a4 * k1 - a4 * k2) / (a3 - a4) * alpha1 * 
     +     exp(-alpha1 * (y + Hbar)) - (a1 * k1 + a2 * k2 - a3 * k1 - 
     +     a3 * k2) / (a3 - a4) * alpha2 * exp(-alpha2 * (y + Hbar))) * 
     +     exp(i * (k * x - omega * t)) + i * epsilon * (a1 * k1 * 
     +     exp(alpha1 * (y + Hbar)) + a2 * k2 * exp(alpha2 * (y + Hbar))
     +     - a3 * (a1 * k1 + a2 * k2 - a4 * k1 - a4 * k2) / (a3 - a4) * 
     +     exp(-alpha1 * (y + Hbar)) + a4 * (a1 * k1 + a2 * k2 - a3 * k1 
     +     - a3 * k2) / (a3 - a4) * exp(-alpha2 * (y + Hbar))) * k * 
     +     exp(i * (k * x - omega * t)))

      s22bar = (lambdabar + 2 * mubar) * epsilon * (a1 * k1 * alpha1 * 
     +     exp(alpha1 * (y + Hbar)) + a2 * k2 * alpha2 * exp(alpha2 * (y 
     +     + Hbar)) + a3 * (a1 * k1 + a2 * k2 - a4 * k1 - a4 * k2) / 
     +     (a3 - a4) * alpha1 * exp(-alpha1 * (y + Hbar)) - a4 * (a1 * 
     +     k1 + a2 * k2 - a3 * k1 - a3 * k2) / (a3 - a4) * alpha2 * 
     +     exp(-alpha2 * (y + Hbar))) * exp(i * (k * x - omega * t)) + 
     +     i * lambdabar * epsilon * (k1 * exp(alpha1 * (y + Hbar)) + 
     +     k2 * exp(alpha2 * (y + Hbar)) - (a1 * k1 + a2 * k2 - a4 * k1 
     +     - a4 * k2) / (a3 - a4) * exp(-alpha1 * (y + Hbar)) + (a1 * k1 
     +     + a2 * k2 - a3 * k1 - a3 * k2) / (a3 - a4) * exp(-alpha2 * (y 
     +     + Hbar))) * k * exp(i * (k * x - omega * t))


c     take real part
      u1barr = dreal(u1bar)
      u2barr = dreal(u2bar)
      v1barr = dreal(v1bar)
      v2barr = dreal(v2bar)
      s11barr = dreal(s11bar)
      s12barr = dreal(s12bar)
      s22barr = dreal(s22bar)


      return
      end


      subroutine evalFibRadialWaveSolid(omegar,omegai, d1r,d1i,d2r,d2i,
     +     d3r,d3i,d4r,d4i, n, rhoBar,muBar,lambdaBar, r,theta,t,
     +     ur,ut,vr,vt,ar,at,srr,srt,stt,sdrr,sdrt,sdtt)
      implicit none
      double precision omegar,omegai, d1r,d1i,d2r,d2i,
     +     d3r,d3i,d4r,d4i, n, rhoBar,muBar,lambdaBar, r,theta,t,
     +     ur,ut,vr,vt,ar,at,srr,srt,stt,sdrr,sdrt,sdtt,
     +     cp,cs, jpnr,jpni,jpnpr,jpnpi,
     +            ypnr,ypni,ypnpr,ypnpi,
     +            jsnr,jsni,jsnpr,jsnpi,
     +            ysnr,ysni,ysnpr,ysnpi
      complex*16 omega,d1,d2,d3,d4, kp,ks,I,
     +     urc,utc,vrc,vtc,arc,atc,srrc,srtc,sttc,sdrrc,sdrtc,sdttc,
     +     rkp,rks, jpn,jpnp,
     +              ypn,ypnp,
     +              jsn,jsnp,
     +              ysn,ysnp

c     get complex vars
      omega = dcmplx(omegar,omegai)
      d1    = dcmplx(   d1r,   d1i)
      d2    = dcmplx(   d2r,   d2i)
      d3    = dcmplx(   d3r,   d3i)
      d4    = dcmplx(   d4r,   d4i)
      I     = dcmplx( 0.0d0, 1.0d0)

      cp = sqrt((lambdaBar+2.0d0*muBar)/(rhoBar))
      cs = sqrt((muBar)/(rhoBar))
      
      kp = omega/cp
      ks = omega/cs

      rkp = r*kp
      rks = r*ks

      if (r.eq.0.0d0) then
        r = 1.0d-16
      end if

c     get bessel functions
      call cBesselJ(n      ,dreal(rkp),dimag(rkp), jpnr, jpni)
      call cBesselJ(n+1.0d0,dreal(rkp),dimag(rkp),jpnpr,jpnpi)
      call cBesselJ(n      ,dreal(rks),dimag(rks), jsnr, jsni)
      call cBesselJ(n+1.0d0,dreal(rks),dimag(rks),jsnpr,jsnpi)
      call cBesselY(n      ,dreal(rkp),dimag(rkp), ypnr, ypni)
      call cBesselY(n+1.0d0,dreal(rkp),dimag(rkp),ypnpr,ypnpi)
      call cBesselY(n      ,dreal(rks),dimag(rks), ysnr, ysni)
      call cBesselY(n+1.0d0,dreal(rks),dimag(rks),ysnpr,ysnpi)
      
      jpn  = dcmplx( jpnr, jpni)
      jpnp = dcmplx(jpnpr,jpnpi)
      jsn  = dcmplx( jsnr, jsni)
      jsnp = dcmplx(jsnpr,jsnpi)
      ypn  = dcmplx( ypnr, ypni)
      ypnp = dcmplx(ypnpr,ypnpi)
      ysn  = dcmplx( ysnr, ysni)
      ysnp = dcmplx(ysnpr,ysnpi)

c     compute solution
      urc = ( d1*(n*jpn/r-kp*jpnp)
     +      +d2*(n*ypn/r-kp*ypnp)
     +      +I*n/r*(d3*jsn+d4*ysn))*exp(I*(n*theta-omega*t))
      utc = (-d3*(n*jsn/r-ks*jsnp)
     +      -d4*(n*ysn/r-ks*ysnp)
     +      +I*n/r*(d1*jpn+d2*ypn))*exp(I*(n*theta-omega*t))

      vrc = -I*omega*urc
      vtc = -I*omega*utc
      
      arc = -I*omega*vrc
      atc = -I*omega*vtc
      
      srrc = ((2.d0*jpnp*d1*kp*muBar*r-(2.d0*I)*jsnp*d3*ks*muBar*n*r
     +     +2.d0*ypnp*d2*kp*muBar*r-(2.d0*I)*ysnp*d4*ks*muBar*n*r
     +     +2.d0*d1*((-kp**2*r**2+n**2-n)*muBar-.5d0*lambdaBar
     +     *kp**2*r**2)*jpn+(2.d0*((-kp**2*r**2+n**2-n)*muBar-.5d0
     +     *lambdaBar*kp**2*r**2))*d2*ypn+(2.d0*I)*n*muBar*(n-1)
     +     *(d3*jsn+d4*ysn))/r**2)
     +     *exp(I*(n*theta-omega*t))

c$$$      srrc = ((2*jpnp*kp*muBar*r+(2*((-kp**2*r**2+n**2-n)*muBar-
c$$$     +     .5d0*lambdaBar*r**2*kp**2))*jpn)/r**2*d1)
c$$$     +     *exp(I*(n*theta-omega*t))
c$$$      srrc = jpnp*exp(I*(n*theta-omega*t))

      srtc = (muBar*(I*n*(d1*(n*jpn/r-kp*jpnp)+d2*(n*ypn/r-kp*ypnp)
     +     +I*n*(d3*jsn+d4*ysn)/r)/r-d3*(n*(-jsnp+n*jsn/(ks*r))*ks/r
     +     -n*jsn/r**2-ks**2*(jsn-(n+1)*jsnp/(ks*r)))-d4*(n*(-ysnp
     +     +n*ysn/(ks*r))*ks/r-n*ysn/r**2-ks**2*(ysn-(n+1)*ysnp/(ks*r)))
     +     -I*n*(d1*jpn+d2*ypn)/r**2+I*n*(d1*(-jpnp+n*jpn/(kp*r))*kp
     +     +d2*(-ypnp+n*ypn/(kp*r))*kp)/r-(-d3*(n*jsn/r-ks*jsnp)
     +     -d4*(n*ysn/r-ks*ysnp)+I*n*(d1*jpn+d2*ypn)/r)/r))
     +     *exp(I*(n*theta-omega*t))

      sttc = (lambdabar*((d1*(n*jpn/r-kp*jpnp)+d2*(n*ypn/r-kp*ypnp)
     +     +I*n*(d3*jsn+d4*ysn)/r+r*(d1*(n*(-jpnp+n*jpn/(kp*r))*kp/r
     +     -n*jpn/r**2-kp**2*(jpn-(n+1)*jpnp/(kp*r)))+d2*(n*(-ypnp
     +     +n*ypn/(kp*r))*kp/r-n*ypn/r**2-kp**2*(ypn-(n+1)*ypnp/(kp*r)))
     +     -I*n*(d3*jsn+d4*ysn)/r**2+I*n*(d3*(-jsnp+n*jsn/(ks*r))*ks
     +     +d4*(-ysnp+n*ysn/(ks*r))*ks)/r))/r+I*n*(-d3*(n*jsn/r-ks*jsnp)
     +     -d4*(n*ysn/r-ks*ysnp)+I*n*(d1*jpn+d2*ypn)/r)/r)+2*muBar*(I*n*
     +     (-d3*(n*jsn/r-ks*jsnp)-d4*(n*ysn/r-ks*ysnp)+I*n*(d1*jpn
     +     +d2*ypn)/r)+d1*(n*jpn/r-kp*jpnp)+d2*(n*ypn/r-kp*ypnp)
     +     +I*n*(d3*jsn+d4*ysn)/r)/r)
     +     *exp(I*(n*theta-omega*t))

      sdrrc = -I*omega*srrc
      sdrtc = -I*omega*srtc
      sdttc = -I*omega*sttc

c     return real parts
      ur = dreal(urc)
      ut = dreal(utc)
      vr = dreal(vrc)
      vt = dreal(vtc)
      ar = dreal(arc)
      at = dreal(atc)
      srr = dreal(srrc)
      srt = dreal(srtc)
      stt = dreal(sttc)
      sdrr = dreal(srrc)
      sdrt = dreal(sdrtc)
      sdtt = dreal(sdttc)

      return
      end







