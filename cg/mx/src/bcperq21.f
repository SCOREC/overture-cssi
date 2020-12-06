      ! original calling sequence
      subroutine bcperq21(p,f,ploc,c,len,dt,n,m,md,ns,ord,fold,phi,amc,
     &                    fftsave,bcinit)
      
      implicit none

c
      integer n,m,bcinit,i,j,k,l,ord,ns,md
      complex*16 phi(md,21,m)
      complex*16 fold(0:ord-2,md,21,m),adon,phat,xfact 
      double precision p(n,m),f(n,m),ploc(n),fftsave(ns)
      double precision amc(-1:ord-2),scl,dt,w,len,c


      call  bcperq21d(1,n,p,f,ploc,c,len,dt,n,m,md,ns,
     &                ord,fold,phi,amc,fftsave,bcinit)
      return
      end 


c*wdh* This version takes the bounds on p and f as input
      subroutine bcperq21d(nda,ndb,p,f,ploc,c,len,dt,n,m,md,ns,
     &                     ord,fold,phi,amc,fftsave,bcinit)
c
c  this routine uses an exponential adams-moulton formula to compute -
c  in Fourier variables - 21-pole approximation to the
c  planar kernel
c
c     (d/dt - c beta_j w) phi_j = c alpha_j w^2 phat 
c
c     fhat = sum phi_j
c
c     w = k*scl , k=1, ... 
c      
c  double precision: p(n,m) - m fields to which the operator should be applied
c
c  double precision: f(n,m) - the m results
c
c  double precision: ploc(n) - workspace 
c
c  double precision: c - the wave speed
c
c  double precision: len - the period
c
c  double precision: dt - the time step
c
c  integer: n the number of grid points - most efficient if it has small
c           prime factors, preferably even
c
c  integer: m the number of fields
c
c  integer: md the maximum mode used in the bc md < n/2
c
c  integer: ns>=2*n+15 
c
c  integer: ord - time-stepping order - note that the stability domain for
c                 Adams-Moulton methods gets small if this is too big
c 
c  complex*16: fold(0:ord-2,md,21,m) - stored values for time-stepping
c
c  complex*16: phi(md,21,m) - the auxiliary functions computed here
c
c  double precision: amc(-1:ord-2) - Adams-Moulton coefficients (computed here)
c                              use amcof.f
c
c  double precision: fftsave(ns) - used by fftpack - link to rffti,rfftf,rfftb
c
c  integer bcinit: initialize to zero 
c
      implicit none

c
      integer n,m,bcinit,i,j,k,l,ord,ns,md
      complex*16 alpha(21),beta(21),phi(md,21,m)
      complex*16 fold(0:ord-2,md,21,m),adon,phat,xfact 
c*wdh      double precision p(n,m),f(n,m),ploc(n),fftsave(ns)
      integer nda,ndb
      double precision p(nda:ndb,m),f(nda:ndb,m),ploc(n),fftsave(ns)
      double precision amc(-1:ord-2),scl,dt,w,len,c
c
      alpha(1)=(-.2410467618025768D-06,  -.2431987763837349D-06)
      alpha(2)=(-.2410467618025768D-06,   .2431987763837349D-06) 
      alpha(3)=(-.1617695923999794D-05,  -.1638622585172068D-05)
      alpha(4)=(-.1617695923999794D-05,   .1638622585172068D-05)
      alpha(5)=(-.7723476507531262D-05,  -.7878743138182415D-05)
      alpha(6)=(-.7723476507531262D-05,   .7878743138182415D-05)
      alpha(7)=(-.3400304516975200D-04,  -.3510673092397324D-04)
      alpha(8)=(-.3400304516975200D-04,   .3510673092397324D-04) 
      alpha(9)=(-.1454893381589074D-03,  -.1535469093409158D-03)
      alpha(10)=(-.1454893381589074D-03,   .1535469093409158D-03)
      alpha(11)=(-.6104572904148162D-03,  -.6733883694898616D-03)
      alpha(12)=(-.6104572904148162D-03,   .6733883694898616D-03)
      alpha(13)=(-.2473202929583869D-02,  -.3011442350813045D-02)
      alpha(14)=(-.2473202929583869D-02,   .3011442350813045D-02)
      alpha(15)=(-.8964957513027030D-02,  -.1398751873403249D-01)
      alpha(16)=(-.8964957513027030D-02,   .1398751873403249D-01)
      alpha(17)=(-.1846252520037211D-01,  -.6565858806543060D-01)
      alpha(18)=(-.1846252520037211D-01,   .6565858806543060D-01)
      alpha(19)=( .9181095934161065D-01,  -.2076825633238755D+00)
      alpha(20)=( .9181095934161065D-01,   .2076825633238755D+00)
      alpha(21)=( .3787484004895032D+00,   .0000000000000000D+00)
c
      beta(1)=(-.4998142304334231D-04,   .9999998607359947D+00)
      beta(2)=(-.4998142304334231D-04,  -.9999998607359947D+00)
      beta(3)=(-.2501648855535112D-03,   .9999990907954994D+00)
      beta(4)=(-.2501648855535112D-03,  -.9999990907954994D+00)
      beta(5)=(-.8021925048752190D-03,   .9999958082358295D+00)
      beta(6)=(-.8021925048752190D-03,  -.9999958082358295D+00)
      beta(7)=(-.2263515963206483D-02,   .9999820162287431D+00)
      beta(8)=(-.2263515963206483D-02,  -.9999820162287431D+00)
      beta(9)=(-.6112737916031916D-02,   .9999224860282032D+00)
      beta(10)=(-.6112737916031916D-02,  -.9999224860282032D+00)
      beta(11)=(-.1625071664643320D-01,   .9996497460330479D+00)
      beta(12)=(-.1625071664643320D-01,  -.9996497460330479D+00)
      beta(13)=(-.4295328074381198D-01,   .9982864080633248D+00)
      beta(14)=(-.4295328074381198D-01,  -.9982864080633248D+00)
      beta(15)=(-.1129636068874967D+00,   .9907617913485537D+00)
      beta(16)=(-.1129636068874967D+00,  -.9907617913485537D+00)
      beta(17)=(-.2902222956062986D+00,   .9462036470847180D+00)
      beta(18)=(-.2902222956062986D+00,  -.9462036470847180D+00)
      beta(19)=(-.6548034445533449D+00,   .7077228221122372D+00)
      beta(20)=(-.6548034445533449D+00,  -.7077228221122372D+00)
      beta(21)=(-.9345542777004186D+00,   .0000000000000000D+00)
c
      if (bcinit.eq.0) then
        !  initialize      
        call amcof(amc,ord)
        do i=1,m
          do j=1,21
            do k=1,md
              phi(k,j,i)=0.d0
              do l=0,ord-2
                fold(l,k,j,i)=0.d0
              end do
            end do
          end do
        end do
        CALL rffti(n,fftsave)
        bcinit=1
      end if
c
      scl=2.d0*3.1415926535897932384d0/len
      do i=1,m
c
c loop over fields
c
        do k=1,n
          ploc(k)=p(k,i)
        end do
        CALL rfftf(n,ploc,fftsave)
        ! write(*,*) "bcp[0] ploc=",(ploc(k),k=1,n)
        ! write(*,*) "scl=",scl
        
        do j=1,21
          do k=1,md
            w=scl*DBLE(k)        
            xfact=exp(c*beta(j)*dt*w)
            adon=xfact*phi(k,j,i)
            ! write(*,*) "* adon=",DREAL(adon),DIMAG(adon)          
            do l=0,ord-2
              adon=adon+xfact*amc(l)*fold(l,k,j,i)
              xfact=xfact*exp(c*beta(j)*dt*w)
              ! write(*,*) "l=",l," adon=",DREAL(adon),DIMAG(adon)          
            end do
            phat=DCMPLX(ploc(2*k),ploc(2*k+1)) 
            adon=adon+c*dt*amc(-1)*alpha(j)*w*w*phat
            phi(k,j,i)=adon
 !           write(*,*) "expBeta=",DREAL(exp(c*beta(j)*dt*w)),
 !     &         DIMAG(exp(c*beta(j)*dt*w))
 !      write(*,*) "bcper: phi=",k,j,i,DREAL(phi(k,j,i)),DIMAG(phi(k,j,i))    
            if (ord.gt.2) then
            do l=ord-2,1,-1
              fold(l,k,j,i)=fold(l-1,k,j,i)
            end do
            end if
            fold(0,k,j,i)=c*dt*w*w*alpha(j)*phat
!      write(*,*) "bcper: fold(0,k,j,i)=",k,j,i,fold(0,k,j,i)
          end do
        end do
        do k=1,n
          ploc(k)=0.d0
        end do
        do j=1,21
          do k=1,md
            ploc(2*k)=ploc(2*k)+DREAL(phi(k,j,i))
            ploc(2*k+1)=ploc(2*k+1)+DIMAG(phi(k,j,i))
          end do
        end do 
        ! write(*,*) "bcp[2] ploc=",(ploc(k),k=1,n)
        CALL rfftb(n,ploc,fftsave)
        do k=1,n
          f(k,i)=ploc(k)/DBLE(n)
        end do 
      end do
c
      return
      end

