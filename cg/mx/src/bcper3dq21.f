      subroutine bcper3dq21(nda1,ndb1,nda2,ndb2,p,f,pl,zl,pl1,zl2,c,
     &                       len1,len2,dt,n1,n2,m,md1,md2,ns1,ns2,
     &                       ord,fold,phi,amc,fftsave1,fftsave2,
     &                       bcinit)
c
c wdh: Routine from Tom Hagstrom, June 2020.
c
c  this routine uses an exponential adams-moulton formula to compute -
c  in Fourier variables - 21-pole approximation to the
c  planar kernel
c
c     (d/dt - c beta_j w) phi_j = c alpha_j w^2 phat 
c
c     fhat = sum phi_j
c
c     w = k1*scl1+k2*scl2 , k1=1, ...  k2=1, ... 
c      
c  double precision: p(n1,n2,m) - m fields to which the operator should be applied
c
c  double precision: f(n1,n2,m) - the m results
c
c  double precision: pl(n1,n2),pl1(n1) - workspace
c
c  complex*16 : zl(n1,n2),zl2(n2) - workspace
c
c  double precision: c - the wave speed
c
c  double precision: len1,len2 - the periods
c
c  double precision: dt - the time step
c
c  integer: n1,n2 the number of grid points in directions 1 and 2 -
c           most efficient if it has small prime factors, preferably even
c
c  integer: m the number of fields
c
c  integer: md1,md2 the maximum mode used in the bc md1 < n1/2, md2 < n2/2
c
c  integer: ns1>=2*n1+15, ns2>=4*n2+15 
c
c  integer: ord - time-stepping order - note that the stability domain for
c                 Adams-Moulton methods gets small if this is too big
c 
c  complex*16: fold(0:ord-2,0:md1,-md2:md2,21,m) - stored values for time-stepping
c
c  complex*16: phi(0:md1,-md2:md2,21,m) - the auxiliary functions computed here
c
c  double precision: amc(-1:ord-2) - Adams-Moulton coefficients (computed here)
c                              use amcof.f
c
c  double precision: fftsave1(ns1) fftsave2(ns2) - used by fftpack
c     - link to dfftb.f dfftf.f dffti.f zfftb.f zfftf.f zffti.f cfftb1.f cfftf1.f cffti1.f
c     rfftb1.f rfftf1.f rffti1.f radf2.f radf3.f radf4.f radf5.f radfg.f passb2.f passb3.f
c     passb4.f passb5.f passf2.f passf3.f passf4.f passf5.f passf.f radb2.f radb3.f radb4.f
c     radb5.f radbg.f passb.f
c
c  integer bcinit: initialize to zero 
c
      implicit none

c
      integer n1,n2,m,bcinit,i,j,k1,k1t,k2,k2t,l,ord,ns1,ns2,md1,md2
      complex*16 alpha(21),beta(21),phi(0:md1,-md2:md2,21,m)
      complex*16 fold(0:ord-2,0:md1,-md2:md2,21,m),adon,phat,xfact,ii 
      integer nda1,ndb1,nda2,ndb2 
      double precision p(nda1:ndb1,nda2:ndb2,m),f(nda1:ndb1,nda2:ndb2,m)
      double precision pl(n1,n2),pl1(n1)
      complex*16 zl(n1,n2),zl2(n2) 
      double precision fftsave1(ns1),fftsave2(ns2)
      double precision amc(-1:ord-2),scl1,scl2,dt,w,len1,len2,c
c
      alpha(1)=(-.2410467618025768E-06,  -.2431987763837349E-06)
      alpha(2)=(-.2410467618025768E-06,   .2431987763837349E-06) 
      alpha(3)=(-.1617695923999794E-05,  -.1638622585172068E-05)
      alpha(4)=(-.1617695923999794E-05,   .1638622585172068E-05)
      alpha(5)=(-.7723476507531262E-05,  -.7878743138182415E-05)
      alpha(6)=(-.7723476507531262E-05,   .7878743138182415E-05)
      alpha(7)=(-.3400304516975200E-04,  -.3510673092397324E-04)
      alpha(8)=(-.3400304516975200E-04,   .3510673092397324E-04) 
      alpha(9)=(-.1454893381589074E-03,  -.1535469093409158E-03)
      alpha(10)=(-.1454893381589074E-03,   .1535469093409158E-03)
      alpha(11)=(-.6104572904148162E-03,  -.6733883694898616E-03)
      alpha(12)=(-.6104572904148162E-03,   .6733883694898616E-03)
      alpha(13)=(-.2473202929583869E-02,  -.3011442350813045E-02)
      alpha(14)=(-.2473202929583869E-02,   .3011442350813045E-02)
      alpha(15)=(-.8964957513027030E-02,  -.1398751873403249E-01)
      alpha(16)=(-.8964957513027030E-02,   .1398751873403249E-01)
      alpha(17)=(-.1846252520037211E-01,  -.6565858806543060E-01)
      alpha(18)=(-.1846252520037211E-01,   .6565858806543060E-01)
      alpha(19)=( .9181095934161065E-01,  -.2076825633238755E+00)
      alpha(20)=( .9181095934161065E-01,   .2076825633238755E+00)
      alpha(21)=( .3787484004895032E+00,   .0000000000000000E+00)
c
      beta(1)=(-.4998142304334231E-04,   .9999998607359947E+00)
      beta(2)=(-.4998142304334231E-04,  -.9999998607359947E+00)
      beta(3)=(-.2501648855535112E-03,   .9999990907954994E+00)
      beta(4)=(-.2501648855535112E-03,  -.9999990907954994E+00)
      beta(5)=(-.8021925048752190E-03,   .9999958082358295E+00)
      beta(6)=(-.8021925048752190E-03,  -.9999958082358295E+00)
      beta(7)=(-.2263515963206483E-02,   .9999820162287431E+00)
      beta(8)=(-.2263515963206483E-02,  -.9999820162287431E+00)
      beta(9)=(-.6112737916031916E-02,   .9999224860282032E+00)
      beta(10)=(-.6112737916031916E-02,  -.9999224860282032E+00)
      beta(11)=(-.1625071664643320E-01,   .9996497460330479E+00)
      beta(12)=(-.1625071664643320E-01,  -.9996497460330479E+00)
      beta(13)=(-.4295328074381198E-01,   .9982864080633248E+00)
      beta(14)=(-.4295328074381198E-01,  -.9982864080633248E+00)
      beta(15)=(-.1129636068874967E+00,   .9907617913485537E+00)
      beta(16)=(-.1129636068874967E+00,  -.9907617913485537E+00)
      beta(17)=(-.2902222956062986E+00,   .9462036470847180E+00)
      beta(18)=(-.2902222956062986E+00,  -.9462036470847180E+00)
      beta(19)=(-.6548034445533449E+00,   .7077228221122372E+00)
      beta(20)=(-.6548034445533449E+00,  -.7077228221122372E+00)
      beta(21)=(-.9345542777004186E+00,   .0000000000000000E+00)
c
      if (bcinit.eq.0) then
c  initialize      
        call amcof(amc,ord)
        do i=1,m
          do j=1,21
            do k2=-md2,md2
	      do k1=0,md1
                phi(k1,k2,j,i)=0.d0
                do l=0,ord-2
                  fold(l,k1,k2,j,i)=0.d0
		end do 
              end do
            end do
          end do
        end do
        CALL dffti(n1,fftsave1)
        CALL zffti(n2,fftsave2)
        bcinit=1
      end if

      ! write(*,'("bcper3dq21: c,dt,len1,len2=",4e12.3)') c,dt,len1,len2
      
c
      ii=(0.d0,1.d0)
      scl1=2.d0*3.1415926535897932384d0/len1
      scl2=2.d0*3.1415926535897932384d0/len2
      do i=1,m
c
c loop over fields
c
        do k2=1,n2
	  do k1=1,n1
            ! wdh pl1(k1)=p(nda1+k1-1,nda2+k2-1,i)
            pl1(k1)=p(k1,k2,i)
          end do
          call dfftf(n1,pl1,fftsave1)
	  zl(1,k2)=pl1(1)
          if (mod(n1,2)==0) then
            do k1=2,n1/2
              zl(k1,k2)=pl1(2*k1-2)+ii*pl1(2*k1-1)
              zl((n1/2)+k1-1,k2)=pl1(2*k1-2)-ii*pl1(2*k1-1)
            end do
            zl(n1,k2)=pl1(n1)
          else
            do k1=2,(n1/2)+1
              zl(k1,k2)=pl1(2*k1-2)+ii*pl1(2*k1-1)
              zl((n1/2)+k1,k2)=pl1(2*k1-2)-ii*pl1(2*k1-1)
            end do
          end if
	end do
        if( .false. )then
          write(*,*) "Stage 1: After initial FFTs, zl:"
          do k2=1,n2
            write(*, "(*('('sf7.3xspf7.3x'i)':x))") (zl(k1,k2),k1=1,n1)
          end do
        end if
      
	do k1=1,n1
	  do k2=1,n2
	    zl2(k2)=zl(k1,k2)
	  end do
	  call zfftf(n2,zl2,fftsave2)
          do k2=1,n2
            zl(k1,k2)=zl2(k2)
          end do   
        end do
        if( .false. )then
          write(*,*) "Stage 1: After 2nd FFTs, zl:"
          do k2=1,n2
            write(*, "(*('('sf7.3xspf7.3x'i)':x))") (zl(k1,k2),k1=1,n1)
          end do
        end if

        
        do j=1,21
          do k2=1,n2
            if (k2.le.(md2+1)) then
              k2t=k2-1
            else if ((k2-1-n2).ge.(-md2)) then 
              k2t=k2-1-n2
            else
              goto 100
            end if
	    do k1t=0,md1
              k1=k1t+1
              w=SQRT((scl1*DBLE(k1t))**2+(scl2*DBLE(k2t))**2)         
              xfact=exp(c*beta(j)*dt*w)
              adon=xfact*phi(k1t,k2t,j,i)
              do l=0,ord-2
                adon=adon+xfact*amc(l)*fold(l,k1t,k2t,j,i)
                xfact=xfact*exp(c*beta(j)*dt*w)
              end do
              adon=adon+c*dt*amc(-1)*alpha(j)*w*w*zl(k1,k2)
              phi(k1t,k2t,j,i)=adon
              if (ord.gt.2) then
              do l=ord-2,1,-1
                fold(l,k1t,k2t,j,i)=fold(l-1,k1t,k2t,j,i)
              end do
              end if
              fold(0,k1t,k2t,j,i)=c*dt*w*w*alpha(j)*zl(k1,k2)
            end do
 100        continue 
          end do
        end do
        do k2=1,n2
          do k1=1,n1
            zl(k1,k2)=(0.d0,0.d0)
          end do  
        end do
        do j=1,21
          do k1t=0,md1
            k1=k1t+1  
	    do k2t=0,md2
              k2=k2t+1
              zl(k1,k2)=zl(k1,k2)+phi(k1t,k2t,j,i)
            end do
            do k2t=-md2,-1
              k2=k2t+1+n2 
              zl(k1,k2)=zl(k1,k2)+phi(k1t,k2t,j,i)
            end do 
          end do
        end do
        do k1t=0,md1
          k1=k1t+1
          do k2=1,n2
            zl2(k2)=zl(k1,k2)
          end do
          call zfftb(n2,zl2,fftsave2)
          do k2=1,n2
            zl(k1,k2)=zl2(k2)
          end do
        end do      
        do k2=1,n2
          pl1(1)=dreal(zl(1,k2))
          do k1=2,md1+1
            pl1(2*k1-2)=dreal(zl(k1,k2))
            pl1(2*k1-1)=dimag(zl(k1,k2))
          end do
          do k1=2*md1+2,n1
            pl1(k1)=0.d0
          end do
          call dfftb(n1,pl1,fftsave1)
          do k1=1,n1
            ! *wdh f(nda1+k1-1,nda2+k2-1,i)=pl1(k1)/DBLE(n1*n2)
            f(k1,k2,i)=pl1(k1)/DBLE(n1*n2)
          end do 
        end do 
      end do
c
      return
      end

