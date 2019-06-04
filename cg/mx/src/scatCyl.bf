!===================================================================================
!  Evaluate the bessel functions J_nu, J_{n+1}, ..  with complex argument z=(zr,ri)
!  
! Input:
!    nu : 
!    np : Number of terms
!       J_{nu+k} k=0,1,2,...,np-1 
!  Output
!      J(k) = (jr(k),ji(k)) k=1,2,...,np
!==================================================================================
      subroutine complexBesselJ( nu, np, zr,zi, jr,ji )

      real nu, zr,zi,jr(1:*),ji(1:*)

      ! real cjr(np), cji(np)
      integer np,nz,ierr,kode

      kode=1 ! do not scale result by exp(-abs(zi))
      call zbesj( zr,zi,nu,kode,np,jr,ji,nz,ierr)

      if( nz.ne.0 .or. ierr.ne.0 )then
        write(*,'("WARNING: zbesj: nz,ierr=",2i4)') nz,ierr
        write(*,'("         nu,zr,zi,jr,ji=",5e12.4)') nu,zr,zi,jr(1),ji(1)
      end if

      ! do k=1,np
      !  jr(k)=cjr(k)
      !  ji(k)=cji(k)
      ! end so

      return 
      end 

!===================================================================================
!  Evaluate the bessel functions Y_nu, Y_{n+1}, ..  with complex argument z=(zr,ri)
!  
! Input:
!    nu : 
!    np : Number of terms
!       Y_{nu+k} k=0,1,2,...,np-1 
!  Output
!      Y(k) = (yr(k),yi(k)) k=1,2,...,np
!==================================================================================
      subroutine complexBesselY( nu, np, zr,zi, yr,yi )

      real nu, zr,zi,yr(1:*),yi(1:*)

      real cwkr(np),cwrki(np)
      ! real cyr(10), cyi(10),cwkr(10),cwrki(10)
      integer np,nz,ierr,kode

      kode=1 ! do not scale result by exp(-abs(zi))
      call zbesy( zr,zi,nu,kode,np,yr,yi,nz,cwkr,cwrki,ierr)

      if( nz.ne.0 .or. ierr.ne.0 )then
        write(*,'("WARNING: zbesy: nz,ierr=",2i4)') nz,ierr
      end if
      !write(*,'(" zbesy: nz,ierr=",2i4)') nz,ierr
      !write(*,'(" zbesy: Y=",2e12.4)') cyr(1),cyi(1)

      ! yr=cyr(1)
      ! yi=cyi(1)

      return 
      end 

! ============================================================================================
! Compute the scattered field solution of electromagnetic scattering from a cylinder
!          a : radius of the cylinder
!          k : wavelength of the incident light
!          m : m=c1/c2 - ratio of speed of sounds
!
!  This solution is taken from 
!    Bowman, Senior and Uslemghi, "Electromagnetic And Acoustic Scattering by Simple Shapes"
!
!  The solution is for an incident wave traveling in the positive x-direction (The opposite direction
!  to the above ref):   
!                     Hz =     exp(i(k*x-w*t))
!                     Ey = -Z* exp(i(k*x-w*t))
!
!
!  nd : =2 number of spcae dimensions
!  n1a:n1b, n2a:n2b, n3a:n3b : evaluate solution at these points, i1=n1a..n1b, 
!  nd1a:nd1b,... : dimensions of u, xy
!  xy : grid points 
!
!  ipar(0) = exr : store the Re part of Ex here
!  ipar(1) = eyr : store the Re part of Ey here
!  ipar(2) = hzr : store the Re part of Hz here
!  ipar(3) = exi : store the Im part of Ex here
!  ipar(4) = eyi : store the Im part of Ey here
!  ipar(5) = hzi : store the Im part of Hz here
! 
!  ipar(6) : option : 0=no dielectric, 1=dielectric
!  ipar(7) : inOut : 0=exterior, 1=interior (for dielectric)
!  ipar(9)= staggeredGrid : 0 = return node centered data, 1= return E and H on a "Yee" type staggered grid.
!  
! =============================================================================================
      subroutine scatcyl(nd, n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,\
                         xy,u,ipar,rpar )

      implicit none
      integer nd, n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b

      real u(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
      real xy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1)

      integer ipar(0:*),option
      real rpar(0:*)

!.............local variables
      integer ntermMax
      parameter( ntermMax=50 )
      integer i1,i2,i3,nterm,ncalc,nb,exr,eyr,hzr,exi,eyi,hzi,staggeredGrid,debug
      real k,a,ka,r,theta,kr,x,y,alpha,twoPi
      real jnka(0:ntermMax),ynka(0:ntermMax),jnpka(0:ntermMax),ynpka(0:ntermMax)
      real jn(0:ntermMax),yn(0:ntermMax),jnp(0:ntermMax),ynp(0:ntermMax),an(0:ntermMax)

      real s,sr,si,srr,srt,sir,sit,rx,ry,tx,ty,pm,cnt,cnp1t,sMax,eps,dir
      real sr0,srr0,srt0,si0,sir0,sit0
      integer n,np1
      integer numEdges,edge

!...............end local variables      

      option = ipar(6)
      if( option.eq.1 )then
        ! dielectric case:
        call scatDielectricCyl(nd, n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,\
                               xy,u,ipar,rpar )
        return
      end if

      k = rpar(0)  !
      a = rpar(1)  ! radius

      exr=ipar(0)
      eyr=ipar(1)
      hzr=ipar(2)
      exi=ipar(3)
      eyi=ipar(4)
      hzi=ipar(5)

      staggeredGrid=ipar(9)
      debug= ipar(10)

      ka=k*a
      write(*,'(" scatcyl: k,a=",2f10.6," exr,eyr,hzr=",3i3," staggeredGrid=",i2)') k,a,exr,eyr,hzr,staggeredGrid
      ! ' 
      if( ka.le.0 )then
        stop 11233
      end if

      ! I estimate that the number of terms, N, should satisfy
      !        ??? N * eps**(1/N) > e*k*a/2 * 1/(2*pi)**(1/N)
      ! where eps=size of the final term (series is alternating)
      ! Take N = max( e*k*a, log(1/eps) )

      eps = 1.e-16
      nterm = max( abs(7.*ka), 16. )
      nterm=min(nterm,ntermMax-2)
      nterm = nterm - mod(nterm,2) + 1   ! nterm should be odd
      ! nterm = 25  ! nterm should be odd
      twoPi=atan2(1.,1.)*8.

      ! First evaluate Jn(ka), Yn(ka) n=0,1,...,nterm
      alpha=0. ! fractional part of Bessel order

      nb = nterm+1  ! eval J0, J1, ... J(nb)  -- compute one extra 
      call rjbesl(ka, alpha, nb, jnka, ncalc)
      call rybesl(ka, alpha, nb, ynka, ncalc)

      ! compute the derivatives 
      
      jnpka(0) = -jnka(1)
      ynpka(0) = -ynka(1)
      do n=1,nterm-1
        jnpka(n) = .5*( jnka(n-1)-jnka(n+1) )
        ynpka(n) = .5*( ynka(n-1)-ynka(n+1) )
      end do

      ! precompute some coefficients
      do n=0,nterm-1
        an(n) = jnpka(n)/( jnpka(n)**2 + ynpka(n)**2 )
      end do

      sMax=0. ! keep track of the size of the last term for monitoring convergence

      dir=-1.


  
      numEdges=1
      if( staggeredGrid.eq.1 )then
       numEdges=3
       n1b=min(n1b,nd1b-1)
       n2b=min(n2b,nd2b-1)
      end if

      do edge=1,numEdges  ! for a staggered grid we need to evaluate along edges and the cell-centers 

       i3=n3a
       do i2=n2a,n2b
       do i1=n1a,n1b
         
         x=dir*xy(i1,i2,i3,0)   ! rotate by Pi so incident wave travels in the positive x-direction
         y=dir*xy(i1,i2,i3,1)
         if( staggeredGrid.eq.1 )then
           if( edge.eq.1 )then
             ! Ex lives on this edge
             x=dir*.5*(xy(i1,i2,i3,0)+xy(i1+1,i2,i3,0) )
           else if( edge.eq.2 )then
             ! Ey lives on this edge
             y=dir*.5*(xy(i1,i2,i3,1)+xy(i1,i2+1,i3,1) )
           else
             ! Hz lives at the cell center 
             x=dir*.5*(xy(i1,i2,i3,0)+xy(i1+1,i2,i3,0) )
             y=dir*.5*(xy(i1,i2,i3,1)+xy(i1,i2+1,i3,1) )
           end if
         end if
 
         r=sqrt(x*x+y*y)
         ! r=max(r,.75*a)    !  don't allow r to get too small -- not valid and convergence is poor
         if( r.lt. .75*a )then
           if( r.lt.eps )then
             x=eps  ! avoid atan(0,0)
             y=eps
           end if
           r=.75*a  !  don't allow r to get too small -- not valid and convergence is poor
         end if
 
         theta=atan2(y,x)
         ! if( theta.lt.0. ) then
         !   theta=theta+twoPi
         ! end if
         ! write(*,'(" i=",2i4," x,y,r,theta=",4f10.5)') i1,i2,x,y,r,theta
 
         kr=k*r
         call rjbesl(kr, alpha, nb, jn, ncalc)
         call rybesl(kr, alpha, nb, yn, ncalc)
         ! derivatives:
         jnp(0) = -jn(1)
         ynp(0) = -yn(1)
         do n=1,nterm-1
           jnp(n) = .5*( jn(n-1)-jn(n+1) )
           ynp(n) = .5*( yn(n-1)-yn(n+1) )
         end do
         
         ! compute the hz field: s and it derivatives sr, sTheta
         !   (sr,si) : holds the Re an Im parts of Hz
         !   srr = d(sr)/dr  srt=d(sr)/d(theta)
         sr=0.
         pm=1.  ! +1 or -1
 
         ! n=0: 
         n=0
         sr = .5*(  jn(n)*jnpka(n) +  yn(n)*ynpka(n))*an(n) 
         srr= .5*( jnp(n)*jnpka(n) + ynp(n)*ynpka(n))*an(n)
         srt= 0.
         si = .5*( -jn(n)*ynpka(n)+   yn(n)*jnpka(n))*an(n)
         sir= .5*(-jnp(n)*ynpka(n)+  ynp(n)*jnpka(n))*an(n)
         sit= 0.
         do n=1,nterm-2,2   ! nterm should be odd
 
           cnt=cos(n*theta)
           sr = sr + pm*( -jn(n)*ynpka(n)+  yn(n)*jnpka(n))*an(n)*cnt
           srr= srr+ pm*(-jnp(n)*ynpka(n)+ ynp(n)*jnpka(n))*an(n)*cnt
           srt= srt+ pm*( -jn(n)*ynpka(n)+  yn(n)*jnpka(n))*an(n)*(-n*sin(n*theta))
 
           si = si - pm*(  jn(n)*jnpka(n) +  yn(n)*ynpka(n) )*an(n)*cnt
           sir= sir- pm*( jnp(n)*jnpka(n) + ynp(n)*ynpka(n) )*an(n)*cnt
           sit= sit- pm*(  jn(n)*jnpka(n) +  yn(n)*ynpka(n) )*an(n)*(-n*sin(n*theta))
           
           np1=n+1
           cnp1t=cos(np1*theta)
           sr = sr - pm*(  jn(np1)*jnpka(np1) +  yn(np1)*ynpka(np1) )*an(np1)*cnp1t
           srr= srr- pm*( jnp(np1)*jnpka(np1) + ynp(np1)*ynpka(np1) )*an(np1)*cnp1t
           srt= srt- pm*(  jn(np1)*jnpka(np1) +  yn(np1)*ynpka(np1) )*an(np1)*(-np1*sin(np1*theta))
 
           si = si - pm*( -jn(np1)*ynpka(np1)+  yn(np1)*jnpka(np1))*an(np1)*cnp1t
           sir= sir- pm*(-jnp(np1)*ynpka(np1)+ ynp(np1)*jnpka(np1))*an(np1)*cnp1t
           sit= sit- pm*( -jn(np1)*ynpka(np1)+  yn(np1)*jnpka(np1))*an(np1)*(-np1*sin(np1*theta))
           pm=-pm
 
           if( n.eq.(nterm-4) )then
             sr0=sr
             srr0=srr
             srt0=srt
 
             si0=si
             sir0=sir
             sit0=sit
           end if
         end do
         ! check the size of the last terms
         sMax = max(sMax,max(abs(sr-sr0),max(abs(srr-srr0),abs(srt-srt0))))
         sMax = max(sMax,max(abs(si-si0),max(abs(sir-sir0),abs(sit-sit0))))
 
 
         sr=-2.*sr
         srr=-2.*srr*k   ! note factor k from Dr( Jn(k*r) )
         srt=-2.*srt
 
         si=-2.*si
         sir=-2.*sir*k   ! note factor k from Dr( Jn(k*r) )
         sit=-2.*sit
 
         rx = x/r  ! r.x 
         ry = y/r
         tx=-sin(theta)/r  ! d(theta)/dx
         ty= cos(theta)/r  ! d(theta)/dy
 
         ! Ex.t=(Hz).y =>  -i*k*Ex = (Hz).y -> Ex = i (Hz).y/k  -> Re(Ex) = -Im(Hz.y)/k  Im(Ex) = Re(Hz.y)
         ! k*Ey =-i*(Hz).x
 
         if( staggeredGrid.eq.0 )then
           ! node centered values 
           u(i1,i2,i3,hzr)=sr
           u(i1,i2,i3,exr)=-dir*(ry*sir+ty*sit)/k  
           u(i1,i2,i3,eyr)= dir*(rx*sir+tx*sit)/k  
  
           u(i1,i2,i3,hzi)=si
           u(i1,i2,i3,exi)= dir*(ry*srr+ty*srt)/k  
           u(i1,i2,i3,eyi)=-dir*(rx*srr+tx*srt)/k  
         else 
           if( edge.eq.1 )then
             ! Ex lives on this edge
            u(i1,i2,i3,exr)=-dir*(ry*sir+ty*sit)/k
            u(i1,i2,i3,exi)= dir*(ry*srr+ty*srt)/k  
           else if( edge.eq.2 )then
             ! Ey lives on this edge
            u(i1,i2,i3,eyr)= dir*(rx*sir+tx*sit)/k  
            u(i1,i2,i3,eyi)=-dir*(rx*srr+tx*srt)/k  
           else
             ! Hz lives at the cell center 
             u(i1,i2,i3,hzr)=sr
             u(i1,i2,i3,hzi)=si
           end if
         end if

         ! rotate to the specfied direction of the incident wave (kx,ky) 
 
         !  u(i1,i2,i3,exr)= fexr
         !  u(i1,i2,i3,eyr)= feyr
         !  u(i1,i2,i3,exi)= fexi
         !  u(i1,i2,i3,eyi)= feyi
 
         
       end do  ! end i1
       end do  ! end i2 
 
       write(*,'(" >>>scatcyl: nterm=",i3," largest final term sMax=",e10.2)') nterm,sMax
       ! '
  
      end do ! edge 

      return
      end

! =============================================================================
! Macro to Evaluate the fields scattered by a dielectric cylinder
!  CASE (input): exterior ot interior
! =============================================================================
#beginMacro dielectricScatteringMacro(CASE)
 ! precompute some coefficients
 do n=0,nterm-1
   ! H = H^(2) = J - i Y
   hnc = cmplx(jnka(n),-ynka(n))
   hnpc= cmplx(jnpka(n),-ynpka(n))
   ! Since H becomes large as n gets large, form the ratio Hn/Hn'  (to avoid cancellation)
   hr = hnc/hnpc
   detc = am*jnmka(n)-hr*jnpmka(n)
   detic = 1./(am*jnmka(n)-hr*jnpmka(n))
   #If #CASE eq "exterior"
     an(n)= (jnka(n)*jnpmka(n)-jnpka(n)*jnmka(n)*am)*detic/hnpc
   #Else
     an(n)= ( -2.*ai*am/(pi*ka) )*detic/hnpc
   #End
 !   write(*,'(" n=",i2," h=(",2e10.2,") hp=(",2e10.2,")  hr=(",2e10.2,") detc=",2e10.2," detic=",2e10.2," an=",2e10.2)') n,real(hnc),aimag(hnc),real(hnpc),aimag(hnpc),real(hr),aimag(hr),real(detc),aimag(detc),real(detic),aimag(detic),real(an(n)),aimag(an(n))

   ! aa = 2.*am/( am**n*(am+1./am))  ! asymptotic formula
    ! write(*,'(" n=",i2," h=(",2e10.2,") an=",2e10.2," an/asymp=",2e10.2)') n,real(hnc),aimag(hnc)\
    !  ,real(an(n)),aimag(an(n)),real(an(n)/aa),aimag(an(n)/aa)

  #If #CASE eq "interior"
   !  if( aa.gt.1e8 )then
   !    an(n)=cmplx(aa,0.)
   !  end if
  #End

 end do
  ! #If #CASE eq "exterior"
  !   write(*,'("dieCyl: exterior: an=",50(2e10.2,2x))') (real(an(n)),aimag(an(n)),n=1,nterm-1)
  ! #Else
  !  write(*,'("dieCyl: interior: an=",50(2e10.2,2x))') (real(an(n)),aimag(an(n)),n=1,nterm-1)
  ! #End 


 numEdges=1
 if( staggeredGrid.eq.1 )then
  numEdges=3
  n1b=min(n1b,nd1b-1)
  n2b=min(n2b,nd2b-1)
 end if

 do edge=1,numEdges  ! for a staggered grid we need to evaluate along edges and the cell-centers 


 sMax=0. ! keep track of the size of the last term for monitoring convergence

 i3=n3a
 do i2=n2a,n2b
 do i1=n1a,n1b
   
   x=xy(i1,i2,i3,0)  
   y=xy(i1,i2,i3,1)
   if( staggeredGrid.eq.1 )then
     if( edge.eq.1 )then
       ! Ex lives on this edge
       x=.5*(xy(i1,i2,i3,0)+xy(i1+1,i2,i3,0) )
     else if( edge.eq.2 )then
       ! Ey lives on this edge
       y=.5*(xy(i1,i2,i3,1)+xy(i1,i2+1,i3,1) )
     else
       ! Hz lives at the cell center 
       x=.5*(xy(i1,i2,i3,0)+xy(i1+1,i2,i3,0) )
       y=.5*(xy(i1,i2,i3,1)+xy(i1,i2+1,i3,1) )
     end if
   end if

   r=sqrt(x*x+y*y)
   !  don't allow r to get too small -- not valid and convergence is poor
   if( r.lt. rMin )then
     if( r.lt.eps )then
       x=eps  ! avoid atan(0,0)
       y=eps
     end if
     r=rMin  !  don't allow r to get too small -- not valid and convergence is poor
   end if

   theta=atan2(y,x)
   ! if( theta.lt.0. ) then
   !   theta=theta+twoPi
   ! end if
   ! write(*,'(" i=",2i4," x,y,r,theta=",4f10.5)') i1,i2,x,y,r,theta

   ! eval Jn(kk*r)
   kr=kk*r
   call rjbesl(kr, alpha, nb, jn, ncalc)
   call rybesl(kr, alpha, nb, yn, ncalc)
   ! derivatives:
   jnp(0) = -jn(1)
   ynp(0) = -yn(1)
   do n=1,nterm
     jnp(n) = .5*( jn(n-1)-jn(n+1) )
     ynp(n) = .5*( yn(n-1)-yn(n+1) )
   end do
   
   ! compute the hz field: s and it derivatives sr, sTheta
   !   (sr,si) : holds the Re an Im parts of Hz
   !   srr = d(sr)/dr  srt=d(sr)/d(theta)
   sr=0.

   ! n=0: 
   n=0
   aimn=1.  ! ai**(0)

   expc=cmplx(cos(n*theta),sin(n*theta))

   ! ---- exterior
   #If #CASE eq "exterior"
     hnc=cmplx(jn(n),-yn(n))
     hnpc=cmplx(jnp(n),-ynp(n))
     sc = (aimn*inc*jn(n) + an(n)*hnc)*expc
     scr= (aimn*inc*jnp(n)+ an(n)*hnpc)*expc
   #Else
    ! --- interior
    sc = an(n)*jn(n)*expc
    scr= an(n)*jnp(n)*expc
   #End

   sr=real(sc)
   srr=real(scr)
   srt=0.
   
   si=aimag(sc)
   sir=aimag(scr)
   sit=0.

   aimn=1. ! aimn = (i)^{-n} = (-i)^n
   ain=1.  ! ain  = (i)^{n} 
   do n=1,nterm-1
     if( mod(n,4).eq.0 )then
       aimn=1.
       ain=1.
     else
       aimn=-aimn*ai
       ain =  ain*ai
     end if

     cnt=cos(n*theta)
     snt=sin(n*theta)

     expc=cmplx(cnt,snt)
     exptc=n*cmplx(-snt,cnt)

     expmc=cmplx(cnt,-snt)
     expmtc=n*cmplx(-snt,-cnt)

     #If #CASE eq "exterior"
       hnc=cmplx(jn(n),-yn(n))
       hnpc=cmplx(jnp(n),-ynp(n))
       sc = (inc*jn(n) + an(n)*hnc )*aimn*expc + (inc*jn(n) + an(n)*hnc )*aimn*expmc
       scr= (inc*jnp(n)+ an(n)*hnpc)*aimn*expc + (inc*jnp(n)+ an(n)*hnpc)*aimn*expmc
       sct= (inc*jn(n) + an(n)*hnc )*aimn*exptc+ (inc*jn(n) + an(n)*hnc )*aimn*expmtc
     #Else
       sc = an(n)*jn(n) *(aimn*expc +aimn*expmc)
       scr= an(n)*jnp(n)*(aimn*expc +aimn*expmc)
       sct= an(n)*jn(n) *(aimn*exptc+aimn*expmtc)
     #End
     ! write(*,'("i1,i2=",2i3," n=",i2," sc=",2f6.3)') i1,i2,n,real(sc),aimag(sc)

     sr=sr  +real(sc)
     srr=srr+real(scr)
     srt=srt+real(sct)
     
     si=si  +aimag(sc)
     sir=sir+aimag(scr)
     sit=sit+aimag(sct)

     if( n.eq.(nterm-2) )then
       sr0=sr
       srr0=srr
       srt0=srt

       si0=si
       sir0=sir
       sit0=sit
     end if
   end do
   ! check the size of the last terms
   sMax = max(sMax,max(abs(sr-sr0),max(abs(srr-srr0),abs(srt-srt0))))
   sMax = max(sMax,max(abs(si-si0),max(abs(sir-sir0),abs(sit-sit0))))

   sc=cmplx(sr,si)
   scr=cmplx(srr,sir)
   sct=cmplx(srt,sit)
   ! #If #CASE eq "exterior"
   !   write(*,'("exterior: i1,i2=",2i3," sc=",2f9.5,", scr=",2f9.5,", sct=",2f9.5)') i1,i2,real(sc),aimag(sc),real(scr),aimag(scr),real(sct),aimag(sct)
   ! #Else
   !   write(*,'("interior: i1,i2=",2i3," sc=",2f9.5,", scr=",2f9.5,", sct=",2f9.5)') i1,i2,real(sc),aimag(sc),real(scr),aimag(scr),real(sct),aimag(sct)
   ! #End

   sr = sr 
   srr= srr*kk   ! note factor kk from Dr( Jn(kk*r) )
   srt= srt
        
   si = si 
   sir= sir*kk   ! note factor kk from Dr( Jn(kk*r) )
   sit= sit

   rx = x/r  ! r.x 
   ry = y/r
   tx=-sin(theta)/r  ! d(theta)/dx
   ty= cos(theta)/r  ! d(theta)/dy

   ! Ex.t=(Hz).y =>  -i*k*Ex = (Hz).y -> Ex = i (Hz).y/k  -> Re(Ex) = -Im(Hz.y)/k  Im(Ex) = Re(Hz.y)
   ! k*Ey =-i*(Hz).x

   ! ***** fix this for exterior/interior ******

   #If #CASE eq "exterior"
     kkm=kk
   #Else
     kkm=kk*am
   #End 

   if( staggeredGrid.eq.0 )then
     ! node centered values 
     ! *wdh* 090515 - H should be -H (didn't matter for SOS but does for FOS)
    u(i1,i2,i3,hzr)=-dir*sr                        
    u(i1,i2,i3,hzi)=-si                        

    u(i1,i2,i3,exr)=-dir*(ry*sir+ty*sit)/kkm  
    u(i1,i2,i3,exi)=     (ry*srr+ty*srt)/kkm

    u(i1,i2,i3,eyr)= dir*(rx*sir+tx*sit)/kkm
    u(i1,i2,i3,eyi)=    -(rx*srr+tx*srt)/kkm
 
    ! write(*,'(" exr,exi=",2(1pe14.4)," eyr,eyi=",2(1pe14.4))')  u(i1,i2,i3,exr),u(i1,i2,i3,exi),u(i1,i2,i3,eyr),u(i1,i2,i3,eyi)

    #If #CASE eq "exterior"
     u(i1,i2,i3,hzr)= u(i1,i2,i3,hzr) - dir*(1.-inc)*cos(k*x)
     u(i1,i2,i3,hzi)= u(i1,i2,i3,hzi) +     (1.-inc)*sin(k*x)

     u(i1,i2,i3,eyr)= u(i1,i2,i3,eyr) - dir*(1.-inc)*cos(k*x)
     u(i1,i2,i3,eyi)= u(i1,i2,i3,eyi) +     (1.-inc)*sin(k*x)
    #End

   else 
     if( edge.eq.1 )then
       ! Ex lives on this edge
      u(i1,i2,i3,exr)=-dir*(ry*sir+ty*sit)/kkm  
      u(i1,i2,i3,exi)= (ry*srr+ty*srt)/kkm
     else if( edge.eq.2 )then
       ! Ey lives on this edge
      u(i1,i2,i3,eyr)= dir*(rx*sir+tx*sit)/kkm
      u(i1,i2,i3,eyi)=-(rx*srr+tx*srt)/kkm
      #If #CASE eq "exterior"
       u(i1,i2,i3,eyr)= u(i1,i2,i3,eyr) - dir*(1.-inc)*cos(k*x)
       u(i1,i2,i3,eyi)= u(i1,i2,i3,eyi) + (1.-inc)*sin(k*x)
      #End
     else
       ! Hz lives at the cell center 
      u(i1,i2,i3,hzr)=-dir*sr                        
      u(i1,i2,i3,hzi)=-si                        
      #If #CASE eq "exterior"
       u(i1,i2,i3,hzr)= u(i1,i2,i3,hzr) - dir*(1.-inc)*cos(k*x)
       u(i1,i2,i3,hzi)= u(i1,i2,i3,hzi) +     (1.-inc)*sin(k*x)
      #End

     end if
   end if


 end do
 end do

 write(*,'(" >>>scatDielectricCyl: nterm=",i3," largest final term sMax=",e10.2)') nterm,sMax
 ! ' 
 end do ! edge

#endMacro





! =============================================================================
! Macro to Evaluate the fields scattered by a dielectric cylinder
!          COMPLEX INDEX OF REFRACTION 
!  CASE (input): exterior ot interior
! =============================================================================
#beginMacro dielectricScatteringComplexMacro(CASE)

 ! precompute some coefficients
 do n=0,nterm-1
   ! H = H^(1) = J + i Y
   hnc = cmplx( jnka(n),  ynka(n))
   hnpc= cmplx(jnpka(n), ynpka(n))
   ! Since H becomes large as n gets large, form the ratio Hn/Hn'  (to avoid cancellation)
   hr = hnc/hnpc
   detc =      mc*jcnmka(n)-hr*jcnpmka(n)
   detic = 1./(mc*jcnmka(n)-hr*jcnpmka(n))
   #If #CASE eq "exterior"
     an(n)= (jnka(n)*jcnpmka(n)-jnpka(n)*jcnmka(n)*mc)*detic/hnpc
   #Else
     an(n)= ( 2.*ai*mc/(pi*ka) )*detic/hnpc
   #End
   !   write(*,'(" n=",i2," h=(",2e10.2,") hp=(",2e10.2,")  hr=(",2e10.2,") detc=",2e10.2," detic=",2e10.2," an=",2e10.2)') n,real(hnc),aimag(hnc),real(hnpc),aimag(hnpc),real(hr),aimag(hr),real(detc),aimag(detc),real(detic),aimag(detic),real(an(n)),aimag(an(n))

   ! aa = 2.*am/( am**n*(am+1./am))  ! asymptotic formula
   ! write(*,'(" n=",i2," h=(",2e10.2,") an=",2e10.2," an/asymp=",2e10.2)') n,real(hnc),aimag(hnc)\
   ! ,real(an(n)),aimag(an(n)),real(an(n)/aa),aimag(an(n)/aa)

 end do
 ! #If #CASE eq "exterior"
 !   write(*,'("dieCylComplex: exterior: an=",50(2e10.2,2x))') (real(an(n)),aimag(an(n)),n=1,nterm-1)
 ! #Else
 !   write(*,'("dieCylComplex: interior: an=",50(2e10.2,2x))') (real(an(n)),aimag(an(n)),n=1,nterm-1)
 ! #End 

 numEdges=1
 if( staggeredGrid.eq.1 )then
  numEdges=3
  n1b=min(n1b,nd1b-1)
  n2b=min(n2b,nd2b-1)
 end if

 do edge=1,numEdges  ! for a staggered grid we need to evaluate along edges and the cell-centers 


 sMax=0. ! keep track of the size of the last term for monitoring convergence

 i3=n3a
 do i2=n2a,n2b
 do i1=n1a,n1b
   
   x=xy(i1,i2,i3,0)  
   y=xy(i1,i2,i3,1)
   if( staggeredGrid.eq.1 )then
     if( edge.eq.1 )then
       ! Ex lives on this edge
       x=.5*(xy(i1,i2,i3,0)+xy(i1+1,i2,i3,0) )
     else if( edge.eq.2 )then
       ! Ey lives on this edge
       y=.5*(xy(i1,i2,i3,1)+xy(i1,i2+1,i3,1) )
     else
       ! Hz lives at the cell center 
       x=.5*(xy(i1,i2,i3,0)+xy(i1+1,i2,i3,0) )
       y=.5*(xy(i1,i2,i3,1)+xy(i1,i2+1,i3,1) )
     end if
   end if

   r=sqrt(x*x+y*y)
   !  don't allow r to get too small -- not valid and convergence is poor
   if( r.lt. rMin )then
     if( r.lt.eps )then
       x=eps  ! avoid atan(0,0)
       y=eps
     end if
     r=rMin  !  don't allow r to get too small -- not valid and convergence is poor
   end if

   theta=atan2(y,x)
   ! if( theta.lt.0. ) then
   !   theta=theta+twoPi
   ! end if
   ! write(*,'(" i=",2i4," x,y,r,theta=",4f10.5)') i1,i2,x,y,r,theta

   ! evaluate Bessel functions
   #If #CASE eq "exterior"
     kr=k*r
     call rjbesl(kr, alpha, nb, jn, ncalc)
     call rybesl(kr, alpha, nb, yn, ncalc)
     ! derivatives:
     jnp(0) = -jn(1)
     ynp(0) = -yn(1)
     do n=1,nterm
       jnp(n) = .5*( jn(n-1)-jn(n+1) )
       ynp(n) = .5*( yn(n-1)-yn(n+1) )
     end do
   #Else
     ! interior -- complex bessel
     mckr=mc*k*r
     zr = real(mckr)
     zi = aimag(mckr)
     call complexBesselJ( nu, nb, zr,zi, jnr,jni )
     do n=0,nterm
       jcn(n)=cmplx(jnr(n),jni(n))
     end do
     ! derivatives:
     jcnp(0) = -jcn(1)
     do n=1,nterm
       jcnp(n) = .5*( jcn(n-1)-jcn(n+1) )
     end do

   #End
   
   ! compute the hz field: s and it derivatives sr, sTheta
   !   (sr,si) : holds the Re an Im parts of Hz
   !   srr = d(sr)/dr  srt=d(sr)/d(theta)

   n=0
   aimn=1.  ! ai**(0)
   ain=1.   ! ai**0

   expc=cmplx(cos(n*theta),sin(n*theta))

   ! ---- exterior
   #If #CASE eq "exterior"
     hnc =cmplx( jn(n),  yn(n))
     hnpc=cmplx(jnp(n), ynp(n))
     sc = an(n)*hnc *expc
     scr= an(n)*hnpc*expc

     ! write(*,'("i1,i2=",2i3," jn=",f6.3)') i1,i2,jn(n)

   #Else
     ! --- interior
     sc = an(n)*jcn(n) *expc
     scr= an(n)*jcnp(n)*expc
   #End

   sct=cmplx(0.,0.)

   ! write(*,'("i1,i2=",2i3," n=",i2," sc=",2f6.3)') i1,i2,n,real(sc),aimag(sc)


   ! aimn=1. ! aimn = (i)^{-n} = (-i)^n
   ain=1.  ! ain  = (i)^{n} 

   ! ------------------ SUM MIE SERIES SOLUTION ----------------------------
   do n=1,nterm-1
     if( mod(n,4).eq.0 )then
       ain=1.
       ! aimn=1.
     else
       ain=ain*ai
       ! *wdh* May 15, 2019 try this
       ! ain=-ain*ai


       ! aimn=-aimn*ai
     end if

     cnt=cos(n*theta)
     snt=sin(n*theta)

     ! expc =  cmplx( cnt,snt)
     ! exptc=n*cmplx(-snt,cnt)

     ! expmc =  cmplx( cnt,-snt)
     ! expmtc=n*cmplx(-snt,-cnt)

     #If #CASE eq "exterior"
       hnc =cmplx(jn(n) ,  yn(n))
       hnpc=cmplx(jnp(n), ynp(n))
       ! H: inc=0 
       sc = sc  + an(n)*hnc*ain*2.*cnt
       ! sc = sc  + (inc*jn(n) + an(n)*hnc )*aimn*expc + (inc*jn(n) + an(n)*hnc )*aimn*expmc
       ! H_r (except for a factor of k)
       scr= scr + an(n)*hnpc*ain*2.*cnt
       ! scr= scr + (inc*jnp(n)+ an(n)*hnpc)*aimn*expc + (inc*jnp(n)+ an(n)*hnpc)*aimn*expmc
       ! H_theta
       sct= sct - an(n)*hnc*ain*2.*n*snt
       ! sct= sct + (inc*jn(n) + an(n)*hnc )*aimn*exptc+ (inc*jn(n) + an(n)*hnc )*aimn*expmtc
     #Else
       ! H 
       sc = sc  + an(n)*jcn(n)*ain*2.*cnt
       ! H_r (except for a factor of m*k)
       scr= scr + an(n)*jcnp(n)*ain*2.*cnt
       ! H_theta
       sct= sct - an(n)*jcn(n)*ain*2.*n*snt
     #End
     ! write(*,'("i1,i2=",2i3," n=",i2," sc=",2f6.3)') i1,i2,n,real(sc),aimag(sc)

     if( n.eq.(nterm-2) )then
       ! Save this solution for the convergence test below
       sc0=sc
       scr0=scr
       sct0=sct
     end if
   end do
   ! check the size of the last terms
   sMax = max(sMax,max(cabs(sc-sc0),max(cabs(scr-scr0),cabs(sct-sct0))))


   ! #If #CASE eq "exterior"
   !   write(*,'("exterior: i1,i2=",2i3," sc=",2f9.5,", scr=",2f9.5,", sct=",2f9.5)') i1,i2,real(sc),aimag(sc),real(scr),aimag(scr),real(sct),aimag(sct)
   ! #Else
   !   write(*,'("interior: i1,i2=",2i3," sc=",2f9.5,", scr=",2f9.5,", sct=",2f9.5)') i1,i2,real(sc),aimag(sc),real(scr),aimag(scr),real(sct),aimag(sct)
   ! #End

   ! scr = scr*kkc    ! note factor kk from Dr( Jn(kk*r) )

   ! sr = sr 
   ! srr= srr*kkc   ! note factor kk from Dr( Jn(kk*r) )
   ! srt= srt
        
   ! si = si 
   ! sir= sir*kkc   ! note factor kk from Dr( Jn(kk*r) )
   ! sit= sit


!-    ! Ex.t=(Hz).y =>  -i*k*Ex = (Hz).y -> Ex = i (Hz).y/k  -> Re(Ex) = -Im(Hz.y)/k  Im(Ex) = Re(Hz.y)
!-    ! k*Ey =-i*(Hz).x
!- 
!-    ! To compute E from Hz we use
!-    !   epsHat*Ex.t =  (Hz).y  --> ss*epsHat*Ex = (Hz).y 
!-    !   epsHat*Ey.t = -(Hz).x  --> ss*epsHat*Ey =-(Hz).x 
!- 
!-    ! NOTE; Factor of kkc from space derivatives of Hz was included above
!-    ! NOTE: factor of (-ai) below to make consistent with old code -- fix me 
!-    #If #CASE eq "exterior"
!-      kkmc=kkc
!-      eFactor = 1./( ss*epsHat1 )
!-    #Else
!-      kkmc=kkc*mc
!-      eFactor = 1./( ss*epsHat2 )
!-    #End 
!- 
!-    ! *** TEST ***
!- 
!-    ! 1/efactor = -I*kkmc 
!-    !! eFactor=1./(-ai*kkmc)
!-    eFactorSave=eFactor
!-    ! eFactor=1./kkmc
!- 
!-    ! Ex = eFactor*( H_y )
!-    ! Ey = eFactor*( -H_x)
!-    scr = scr*eFactor
!-    sct = sct*eFactor   
!- 
!-    sr  =  real(sc)
!-    si  = aimag(sc)
!-    srr =  real(scr)
!-    sir = aimag(scr)
!-    srt =  real(sct)
!-    sit = aimag(sct)


   rx = x/r  ! r.x 
   ry = y/r  ! r.y 
   tx=-sin(theta)/r  ! d(theta)/dx
   ty= cos(theta)/r  ! d(theta)/dy

   if( staggeredGrid.eq.0 )then
     ! node centered values 

     ! *wdh* May 15, 2019: CHECK THIS 
     u(i1,i2,i3,hzr)=real(sc)
     u(i1,i2,i3,hzi)=aimag(sc)


     ! eEx = 1/(ss*epsHat) * (Hz)_y = 1/(ss*epsHat) * ( k*r_y*(Hz)_r  + theta_y*(Hz)_theta )
     ! Ey = - 1/(ss*epsHat) * (Hz)_x = - 1/(ss*epsHat) * ( k*r_x*(Hz)_r  + theta_x*(Hz)_theta )
     #If #CASE eq "exterior"
       exField = ( k*ry*scr + ty*sct )/( ss*epsHat1 )
       eyField =-( k*rx*scr + tx*sct )/( ss*epsHat1 )
     #Else
       exField = ( k*mc*ry*scr + ty*sct )/( ss*epsHat2 )
       eyField =-( k*mc*rx*scr + tx*sct )/( ss*epsHat2 )
     #End

     u(i1,i2,i3,exr)=  real(exField)
     u(i1,i2,i3,exi)= aimag(exField)

     u(i1,i2,i3,eyr)=  real(eyField)
     u(i1,i2,i3,eyi)= aimag(eyField)
 
     ! *wdh* May 15, 2019:  try this 
     !u(i1,i2,i3,exr)=-  real(exField)
     !u(i1,i2,i3,exi)= aimag(exField)

     !u(i1,i2,i3,eyr)=-  real(eyField)
     !u(i1,i2,i3,eyi)= aimag(eyField)
 
    ! write(*,'(" exr,exi=",2(1pe14.4)," eyr,eyi=",2(1pe14.4))')  u(i1,i2,i3,exr),u(i1,i2,i3,exi),u(i1,i2,i3,eyr),u(i1,i2,i3,eyi)

    #If #CASE eq "exterior"
     ! Add on the incident field

     ! Hz-incident = e^{i k x}
     coskx=cos(k*x)
     sinkx=sin(k*x)
     u(i1,i2,i3,hzr)= u(i1,i2,i3,hzr) + coskx
     u(i1,i2,i3,hzi)= u(i1,i2,i3,hzi) + sinkx

     !   Ey-incident =   (- i k)/( s*epsHat) * exp(i*k*x )
     eyIncident = (-ai*k)/(ss*epsHat1)*cmplx(coskx,sinkx)

     u(i1,i2,i3,eyr)= u(i1,i2,i3,eyr) +  real(eyIncident)
     u(i1,i2,i3,eyi)= u(i1,i2,i3,eyi) + aimag(eyIncident)

     ! old: 
     ! u(i1,i2,i3,eyr)= u(i1,i2,i3,eyr) - dir*(1.-inc)*cos(k*x)
     ! u(i1,i2,i3,eyi)= u(i1,i2,i3,eyi) +     (1.-inc)*sin(k*x)
    #End

   else 
     if( edge.eq.1 )then
       ! Ex lives on this edge
      u(i1,i2,i3,exr)=-dir*(ry*sir+ty*sit)
      u(i1,i2,i3,exi)=     (ry*srr+ty*srt)
     else if( edge.eq.2 )then
       ! Ey lives on this edge
      u(i1,i2,i3,eyr)= dir*(rx*sir+tx*sit)
      u(i1,i2,i3,eyi)=    -(rx*srr+tx*srt)
      #If #CASE eq "exterior"
       u(i1,i2,i3,eyr)= u(i1,i2,i3,eyr) - dir*(1.-inc)*cos(k*x)
       u(i1,i2,i3,eyi)= u(i1,i2,i3,eyi) + (1.-inc)*sin(k*x)
      #End
     else
       ! Hz lives at the cell center 
      u(i1,i2,i3,hzr)=-dir*sr                        
      u(i1,i2,i3,hzi)=-si                        
      #If #CASE eq "exterior"
       u(i1,i2,i3,hzr)= u(i1,i2,i3,hzr) - dir*(1.-inc)*cos(k*x)
       u(i1,i2,i3,hzi)= u(i1,i2,i3,hzi) +     (1.-inc)*sin(k*x)
      #End

     end if
   end if


 end do
 end do

 write(*,'(" >>>scatDieCylComplex: nterm=",i3," largest final term sMax=",e10.2)') nterm,sMax
 ! ' 
 end do ! edge

#endMacro


! ============================================================================================
! Compute the field of electromagnetic scattering from a *dielectric* cylinder
! Both the field exterior and the field interior to the cylinder can be computed
!          a : radius of the cylinder
!          k : wavelength of the incident light
!
!  This solution is taken from 
!    Balanis, "Advanced Engineering Eletromagnetics", p666 (problem 11.27)
!
!  The solution is for an incident wave traveling in the positive x-direction 
!                     Hz =     exp(i(k*x-w*t))
!                     Ey = -Z* exp(i(k*x-w*t))
!
!
!  nd : =2 number of space dimensions
!  n1a:n1b, n2a:n2b, n3a:n3b : evaluate solution at these points, i1=n1a..n1b, 
!  nd1a:nd1b,... : dimensions of u, xy
!  xy : grid points 
!
!  rpar(0) : k
!  rpar(1) : a 
!  rpar(2) : m 
!
!  ipar(0) = exr : store the Re part of Ex here
!  ipar(1) = eyr : store the Re part of Ey here
!  ipar(2) = hzr : store the Re part of Hz here
!  ipar(3) = exi : store the Im part of Ex here
!  ipar(4) = eyi : store the Im part of Ey here
!  ipar(5) = hzi : store the Im part of Hz here
!
!  ipar(6) : option : 0=no dielectric, 1=dielectric
!  ipar(7) : inOut : 0=exterior, 1=interior
! =============================================================================================
      subroutine scatDielectricCyl(nd, n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,\
                         xy,u,ipar,rpar )

      implicit none
      integer nd, n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b

      real u(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
      real xy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1)

      integer ipar(0:*)
      real rpar(0:*)

      ! Dispersion models
      #Include "dispersionModelsFortranInclude.h"

!.............local variables
      integer ntermMax
      parameter( ntermMax=50 )
      integer i1,i2,i3,nterm,ncalc,nb,exr,eyr,hzr,exi,eyi,hzi,staggeredGrid,debug
      real k,a,ka,r,theta,kr,x,y,alpha,twoPi,pi,coskx,sinkx
      real jnka(0:ntermMax),ynka(0:ntermMax),jnpka(0:ntermMax),ynpka(0:ntermMax)
      real jn(0:ntermMax),yn(0:ntermMax),jnp(0:ntermMax),ynp(0:ntermMax)

      real s,sr,si,srr,srt,sir,sit,rx,ry,tx,ty,pm,cnt,cnp1t,sMax,eps,dir
      real sr0,srr0,srt0,si0,sir0,sit0
      integer n,np1

      integer numEdges,edge
      integer outside,inside
      parameter( outside=0, inside=1)
      integer option,inOut
      integer dispersionModel
      real jnmka(0:ntermMax),jnpmka(0:ntermMax),ynmka(0:ntermMax),ynpmka(0:ntermMax)
      real am,mka,snt,kk,rMin,aa,kkm,inc
      real mr,mi 
      complex*16 ai,hnc,hnpc,detc,detic,aimn,expc,exptc,sc,scr,sct,ain,expmc,expmtc,hr
      complex*16 exField,eyField
      complex*16 an(0:ntermMax)

      ! for complex impedance
      logical useComplexRefractiveIndex
      real ssr,ssi,epsHat1r,epsHat1i,epsHat2r,epsHat2i
      complex*16 mc,mcka,mckr,kkc,kkmc,sc0,scr0,sct0, ss,epsHat1,epsHat2,eFactor,eFactorSave,eScale,eyIncident
      ! for complex Bessel:
      real zr,zi,nu,errr,erri
      real jnr(0:ntermMax),jni(0:ntermMax),ynr(0:ntermMax),yni(0:ntermMax)
      complex*16 jcn(0:ntermMax),jcnp(0:ntermMax),jcnmka(0:ntermMax),jcnpmka(0:ntermMax)

!...............end local variables      

      k        = rpar(0)  !
      a        = rpar(1)  ! radius
      mr       = rpar(2)  ! m = mr + I*mi = index of refraction (complex for dispersive models )
      mi       = rpar(3)  ! 
      ssr      = rpar(4)  ! exp(ss*t) = time-dependence
      ssi      = rpar(5)
      epsHat1r = rpar(6)
      epsHat1i = rpar(7)
      epsHat2r = rpar(8)
      epsHat2i = rpar(9)

      ss = cmplx(ssr,ssi)
      ! epsHat1= eps1*( 1 + chi1(s) )
      epsHat1=cmplx(epsHat1r,epsHat1i)  
      epsHat2=cmplx(epsHat2r,epsHat2i)  

      am = mr      ! for real index of refraction

      exr             = ipar(0)
      eyr             = ipar(1)
      hzr             = ipar(2)
      exi             = ipar(3)
      eyi             = ipar(4)
      hzi             = ipar(5)
      option          = ipar(6)  ! should be 1
      inOut           = ipar(7)  ! 0 = outside, 1=inside

      staggeredGrid   = ipar(9)
      debug           = ipar(10)
      dispersionModel = ipar(11)

      if( dispersionModel.ne.noDispersion )then
        useComplexRefractiveIndex=.true.
      else
        useComplexRefractiveIndex=.false.
        if( mi .ne. 0. )then
          write(*,'(" scatDieCyl:ERROR mi!=0 for non-disersive case")') 
          stop 1111
        end if
      end if

      ai=cmplx(0.,1.)  ! I = sqrt(-1)

      ka=k*a
      mka=am*k*a   ! wave number inside the dielectric is m*k 

      mc=cmplx(mr,mi) ! complex index of refraction
      mcka = mc*k*a


      inc=0.  ! do not compute incident wave by a series

      write(*,'(" scatDieCyl: k,a=",2f10.6,", m=(",f10.6,",",f10.6,") exr,eyr,hzr=",3i3)') k,a,mr,mi,exr,eyr,hzr
      write(*,'(" scatDieCyl: mc*k=(",f10.6,",",f10.6,")")') real(mc*k),aimag(mc*k)
      if( ka.le.0 )then
        stop 12233
      end if

      ! I estimate that the number of terms, N, should satisfy
      !        ??? N * eps**(1/N) > e*k*a/2 * 1/(2*pi)**(1/N)
      ! where eps=size of the final term (series is alternating)
      ! Take N = max( e*k*a, log(1/eps) )

      eps = 1.e-16  ! ********** fix this  -- should be REAL_EPSILON
      nterm = max( abs(7.*ka), 16. )

      nterm=min(nterm,30)  ! do this for now ****

      nterm=min(nterm,ntermMax-2)
      ! ** nterm = nterm - mod(nterm,2) + 1   ! nterm should be odd
      ! nterm = 25  ! nterm should be odd
      twoPi=atan2(1.,1.)*8.
      pi=twoPi*.5

      if( .not.useComplexRefractiveIndex )then
        ! ----------------------------------------
        ! ------- REAL REFRACTIVE INDEX -----------
        ! ----------------------------------------

        ! First evaluate Jn(ka), Yn(ka) n=0,1,...,nterm
        alpha=0. ! fractional part of Bessel order
  
        nb = nterm+1  ! eval J0, J1, ... J(nb)  -- compute one extra 
        call rjbesl(ka, alpha, nb, jnka, ncalc)
        call rybesl(ka, alpha, nb, ynka, ncalc)
        ! also eval Jn(mka), Yn(mka) 
        call rjbesl(mka, alpha, nb, jnmka, ncalc)
        call rybesl(mka, alpha, nb, ynmka, ncalc)
  
        ! compute the derivatives 
        
        jnpka(0) = -jnka(1)
        ynpka(0) = -ynka(1)
        jnpmka(0) = -jnmka(1)
        ynpmka(0) = -ynmka(1)
        do n=1,nterm
          jnpka(n) = .5*( jnka(n-1)-jnka(n+1) )
          ynpka(n) = .5*( ynka(n-1)-ynka(n+1) )
  
          jnpmka(n) = .5*( jnmka(n-1)-jnmka(n+1) )
          ynpmka(n) = .5*( ynmka(n-1)-ynmka(n+1) )
        end do
  
        dir=-1.
        if( inOut.eq.outside )then
          kk=k
          rMin=.75*a
          dielectricScatteringMacro(exterior)
        else
          kk=am*k
          rMin=eps*sqrt(2.)  ! need sqrt(2.) for correct evaluation at origin *wdh* 061008
          dielectricScatteringMacro(interior)
        end if
  
      else

        ! --------------------------------------------
        ! ------- COMPLEX REFRACTIVE INDEX -----------
        ! --------------------------------------------

        ! First evaluate Jn(ka), Yn(ka) n=0,1,...,nterm
        nu=0. ! fractional part of Bessel order
        alpha=nu
        nb = nterm+1  ! eval J0, J1, ... J(nb)  -- compute one extra 
        call rjbesl(ka, nu, nb, jnka, ncalc)
        call rybesl(ka, nu, nb, ynka, ncalc)

        ! also eval complex Jn(mcka) for interior
        zr = real(mcka)
        zi = aimag(mcka)
        call complexBesselJ( nu, nb, zr,zi, jnr,jni )
        do n=0,nterm
          jcnmka(n)=cmplx(jnr(n),jni(n))
        end do
  
       !         if( .false. )then
       !           ! --- test eval of complex bessel routines ---
       !           call rjbesl(mka, alpha, nb, jnmka, ncalc)
       !           call rybesl(mka, alpha, nb, ynmka, ncalc)
       ! 
       ! 
       !           nu=alpha
       !           zr = mka 
       !           zi = 0. 
       !           call complexBesselJ( nu, nb, zr,zi, jnr,jni )
       !           do n=0,nterm
       !             jcn(n)=cmplx(jnr(n),jni(n))
       !           end do
       !   
       !           errr=0.
       !           erri=0.
       !           do n=0,nterm
       !             ! write(*,'(" n=",i2," jnr,jnmka=",2(1pe14.4))') n,jnr(n),jnmka(n)
       !             errr = max(errr,abs((jnr(n)-jnmka(n))/jnr(n)))
       !             erri = max(erri,abs((jni(n)-0.      )))
       !           end do
       !           write(*,'(" scatCyl: Max rel-err in complex Bessel Jn=",2(1pe12.2))') errr,erri
       !   
       !           call complexBesselY( nu, nb, zr,zi, ynr,yni )
       !           do n=0,nterm
       !             ync(n)=cmplx(ynr(n),yni(n))
       !           end do
       !   
       !           errr=0.
       !           erri=0.
       !           do n=0,nterm
       !             ! write(*,'(" n=",i2," ynr,ynmka=",2(1pe14.4))') n,ynr(n),ynmka(n)
       !             write(*,'(" n=",i2," yni      =",1(1pe14.4))') n,yni(n)
       ! 
       !             errr = max(errr,abs((ynr(n)-ynmka(n))/ynr(n)))
       !             erri = max(erri,abs((yni(n)-0.      )/(max(1e-12,abs(ynr(n))))))
       !           end do
       !           write(*,'(" scatCyl: Max rel-err in complex Bessel Yn=",2(1pe12.2))') errr,erri
       !         
       !         end if
  
        ! compute the derivatives 
        
        jnpka(0) = -jnka(1)
        ynpka(0) = -ynka(1)

        jcnpmka(0) = -jcnmka(1)
        do n=1,nterm
          jnpka(n) = .5*( jnka(n-1)-jnka(n+1) )
          ynpka(n) = .5*( ynka(n-1)-ynka(n+1) )
  
          jcnpmka(n) = .5*( jcnmka(n-1)-jcnmka(n+1) )
        end do
  
        dir=-1.
        if( inOut.eq.outside )then
          kkc=k
          rMin=.75*a
          dielectricScatteringComplexMacro(exterior)
        else
          kkc=mc*k
          rMin=eps*sqrt(2.)  ! need sqrt(2.) for correct evaluation at origin *wdh* 061008
          dielectricScatteringComplexMacro(interior)
        end if

      eScale = -ai*k*k*mc*mc/ss

      write(*,'(" ss       =",2(1pe14.4))') real(ss        ),aimag(ss        )
      write(*,'(" epsHat1  =",2(1pe14.4))') real(epsHat1   ),aimag(epsHat1   )
      write(*,'(" epsHat2  =",2(1pe14.4))') real(epsHat2   ),aimag(epsHat2   )
      ! write(*,'(" 1/eFactor=",2(1pe14.4))') real(1./eFactorSave),aimag(1./eFactorSave)
      write(*,'(" kkmc     =",2(1pe14.4))') real(kkmc      ),aimag(kkmc      )
      ! write(*,'(" eScale   =",2(1pe14.4))') real(eScale    ),aimag(eScale    )


      end if ! complex refractive index


      write(*,'(" >>>inOut=",i3," rMin=",e10.2)') inOut,rMin
      flush(6)
      ! stop 5678


      return
      end

