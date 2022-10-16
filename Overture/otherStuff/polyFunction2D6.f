! This file automatically generated from polyFunction.bf with bpp.
! polyFun(2D6)
      subroutine polyFunction2D6 (nd,ndra,ndrb,ndsa,ndsb,ndta,ndtb,
     & ndrra,ndrrb,ndrsa,ndrsb,ndrta,ndrtb,ndrca,ndrcb,ndc1,ndc2,ndc3,
     & nra,nrb,nsa,nsb,nta,ntb, nca,ncb, nda, degree, degreeTime, t, 
     & a,c, r,xa,ya,za, dx,dy,dz,dt)
      !==========================================================================
      !   *** Define a polynomial function and it's derivatives ***
      ! nd : number of space dimensions
      ! nra,nrb,nsa,nsb,nta,ntb : return result in this array
      ! c : array of polynomial coefficients
      ! r  : return result in this array
      ! degree: degree of the polynomial
      ! dx,dy,dz: compute this derivative
      !==========================================================================
      integer nca,ncb,dx,dy,dz,dt,degree,degreeTime
      real xa(ndra:ndrb,ndsa:ndsb,ndta:ndtb)
      real ya(ndra:ndrb,ndsa:ndsb,ndta:ndtb)
      real za(ndra:ndrb,ndsa:ndsb,ndta:ndtb)
      ! real r(ndrra:ndrrb,ndrsa:ndrsb,ndrta:ndrtb,nca:ncb)
      real r(ndrra:ndrrb,ndrsa:ndrsb,ndrta:ndrtb,ndrca:ndrcb)
      real c(0:ndc1-1,0:ndc2-1,0:ndc3-1,0:*)
      real a(0:nda-1,0:*)
      real t
      integer n,i1,i2,i3,laplace
      real time
      laplace=0
      if( dx.eq.-2 )then
      laplace=1
      end if
! poly2D6(x1=xa(i1,i2,i3),y1=ya(i1,i2,i3),z1=za(i1,i2,i3),r(i1,i2,i3,n))
      if( dx.eq.0.and.dy.eq.0.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      y1=ya(i1,i2,i3)
      y2=y1*y1
      y3=y2*y1
      y4=y3*y1
      y5=y4*y1
      y6=y5*y1
      x1=xa(i1,i2,i3)
      x1y1=x1*y1
      x1y2=x1y1*y1
      x1y3=x1y2*y1
      x1y4=x1y3*y1
      x1y5=x1y4*y1
      x1y6=x1y5*y1
      x2=x1*x1
      x2y1=x1y1*x1
      x2y2=x1y2*x1
      x2y3=x1y3*x1
      x2y4=x1y4*x1
      x2y5=x1y5*x1
      x2y6=x1y6*x1
      x3=x2*x1
      x3y1=x2y1*x1
      x3y2=x2y2*x1
      x3y3=x2y3*x1
      x3y4=x2y4*x1
      x3y5=x2y5*x1
      x3y6=x2y6*x1
      x4=x3*x1
      x4y1=x3y1*x1
      x4y2=x3y2*x1
      x4y3=x3y3*x1
      x4y4=x3y4*x1
      x4y5=x3y5*x1
      x4y6=x3y6*x1
      x5=x4*x1
      x5y1=x4y1*x1
      x5y2=x4y2*x1
      x5y3=x4y3*x1
      x5y4=x4y4*x1
      x5y5=x4y5*x1
      x5y6=x4y6*x1
      x6=x5*x1
      x6y1=x5y1*x1
      x6y2=x5y2*x1
      x6y3=x5y3*x1
      x6y4=x5y4*x1
      x6y5=x5y5*x1
      x6y6=x5y6*x1
      r(i1,i2,i3,n)=(c(0,0,0,n)+c(0,1,0,n)*y1+c(0,2,0,n)*y2+c(0,3,0,n)*
     & y3+c(0,4,0,n)*y4+c(0,5,0,n)*y5+c(0,6,0,n)*y6+c(1,0,0,n)*x1+c(1,
     & 1,0,n)*x1y1+c(1,2,0,n)*x1y2+c(1,3,0,n)*x1y3+c(1,4,0,n)*x1y4+c(
     & 1,5,0,n)*x1y5+c(1,6,0,n)*x1y6+c(2,0,0,n)*x2+c(2,1,0,n)*x2y1+c(
     & 2,2,0,n)*x2y2+c(2,3,0,n)*x2y3+c(2,4,0,n)*x2y4+c(2,5,0,n)*x2y5+
     & c(2,6,0,n)*x2y6+c(3,0,0,n)*x3+c(3,1,0,n)*x3y1+c(3,2,0,n)*x3y2+
     & c(3,3,0,n)*x3y3+c(3,4,0,n)*x3y4+c(3,5,0,n)*x3y5+c(3,6,0,n)*
     & x3y6+c(4,0,0,n)*x4+c(4,1,0,n)*x4y1+c(4,2,0,n)*x4y2+c(4,3,0,n)*
     & x4y3+c(4,4,0,n)*x4y4+c(4,5,0,n)*x4y5+c(4,6,0,n)*x4y6+c(5,0,0,n)
     & *x5+c(5,1,0,n)*x5y1+c(5,2,0,n)*x5y2+c(5,3,0,n)*x5y3+c(5,4,0,n)*
     & x5y4+c(5,5,0,n)*x5y5+c(5,6,0,n)*x5y6+c(6,0,0,n)*x6+c(6,1,0,n)*
     & x6y1+c(6,2,0,n)*x6y2+c(6,3,0,n)*x6y3+c(6,4,0,n)*x6y4+c(6,5,0,n)
     & *x6y5+c(6,6,0,n)*x6y6)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.eq.0.and.dy.eq.1.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      y1=ya(i1,i2,i3)
      y2=y1*y1
      y3=y2*y1
      y4=y3*y1
      y5=y4*y1
      x1=xa(i1,i2,i3)
      x1y1=x1*y1
      x1y2=x1y1*y1
      x1y3=x1y2*y1
      x1y4=x1y3*y1
      x1y5=x1y4*y1
      x2=x1*x1
      x2y1=x1y1*x1
      x2y2=x1y2*x1
      x2y3=x1y3*x1
      x2y4=x1y4*x1
      x2y5=x1y5*x1
      x3=x2*x1
      x3y1=x2y1*x1
      x3y2=x2y2*x1
      x3y3=x2y3*x1
      x3y4=x2y4*x1
      x3y5=x2y5*x1
      x4=x3*x1
      x4y1=x3y1*x1
      x4y2=x3y2*x1
      x4y3=x3y3*x1
      x4y4=x3y4*x1
      x4y5=x3y5*x1
      x5=x4*x1
      x5y1=x4y1*x1
      x5y2=x4y2*x1
      x5y3=x4y3*x1
      x5y4=x4y4*x1
      x5y5=x4y5*x1
      x6=x5*x1
      x6y1=x5y1*x1
      x6y2=x5y2*x1
      x6y3=x5y3*x1
      x6y4=x5y4*x1
      x6y5=x5y5*x1
      r(i1,i2,i3,n)=(c(0,1,0,n)+c(0,2,0,n)*y1*2.+c(0,3,0,n)*y2*3.+c(0,
     & 4,0,n)*y3*4.+c(0,5,0,n)*y4*5.+c(0,6,0,n)*y5*6.+c(1,1,0,n)*x1+c(
     & 1,2,0,n)*x1y1*2.+c(1,3,0,n)*x1y2*3.+c(1,4,0,n)*x1y3*4.+c(1,5,0,
     & n)*x1y4*5.+c(1,6,0,n)*x1y5*6.+c(2,1,0,n)*x2+c(2,2,0,n)*x2y1*2.+
     & c(2,3,0,n)*x2y2*3.+c(2,4,0,n)*x2y3*4.+c(2,5,0,n)*x2y4*5.+c(2,6,
     & 0,n)*x2y5*6.+c(3,1,0,n)*x3+c(3,2,0,n)*x3y1*2.+c(3,3,0,n)*x3y2*
     & 3.+c(3,4,0,n)*x3y3*4.+c(3,5,0,n)*x3y4*5.+c(3,6,0,n)*x3y5*6.+c(
     & 4,1,0,n)*x4+c(4,2,0,n)*x4y1*2.+c(4,3,0,n)*x4y2*3.+c(4,4,0,n)*
     & x4y3*4.+c(4,5,0,n)*x4y4*5.+c(4,6,0,n)*x4y5*6.+c(5,1,0,n)*x5+c(
     & 5,2,0,n)*x5y1*2.+c(5,3,0,n)*x5y2*3.+c(5,4,0,n)*x5y3*4.+c(5,5,0,
     & n)*x5y4*5.+c(5,6,0,n)*x5y5*6.+c(6,1,0,n)*x6+c(6,2,0,n)*x6y1*2.+
     & c(6,3,0,n)*x6y2*3.+c(6,4,0,n)*x6y3*4.+c(6,5,0,n)*x6y4*5.+c(6,6,
     & 0,n)*x6y5*6.)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.eq.0.and.dy.eq.2.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      y1=ya(i1,i2,i3)
      y2=y1*y1
      y3=y2*y1
      y4=y3*y1
      x1=xa(i1,i2,i3)
      x1y1=x1*y1
      x1y2=x1y1*y1
      x1y3=x1y2*y1
      x1y4=x1y3*y1
      x2=x1*x1
      x2y1=x1y1*x1
      x2y2=x1y2*x1
      x2y3=x1y3*x1
      x2y4=x1y4*x1
      x3=x2*x1
      x3y1=x2y1*x1
      x3y2=x2y2*x1
      x3y3=x2y3*x1
      x3y4=x2y4*x1
      x4=x3*x1
      x4y1=x3y1*x1
      x4y2=x3y2*x1
      x4y3=x3y3*x1
      x4y4=x3y4*x1
      x5=x4*x1
      x5y1=x4y1*x1
      x5y2=x4y2*x1
      x5y3=x4y3*x1
      x5y4=x4y4*x1
      x6=x5*x1
      x6y1=x5y1*x1
      x6y2=x5y2*x1
      x6y3=x5y3*x1
      x6y4=x5y4*x1
      r(i1,i2,i3,n)=(c(0,2,0,n)*2.+c(0,3,0,n)*y1*6.+c(0,4,0,n)*y2*12.+
     & c(0,5,0,n)*y3*20.+c(0,6,0,n)*y4*30.+c(1,2,0,n)*x1*2.+c(1,3,0,n)
     & *x1y1*6.+c(1,4,0,n)*x1y2*12.+c(1,5,0,n)*x1y3*20.+c(1,6,0,n)*
     & x1y4*30.+c(2,2,0,n)*x2*2.+c(2,3,0,n)*x2y1*6.+c(2,4,0,n)*x2y2*
     & 12.+c(2,5,0,n)*x2y3*20.+c(2,6,0,n)*x2y4*30.+c(3,2,0,n)*x3*2.+c(
     & 3,3,0,n)*x3y1*6.+c(3,4,0,n)*x3y2*12.+c(3,5,0,n)*x3y3*20.+c(3,6,
     & 0,n)*x3y4*30.+c(4,2,0,n)*x4*2.+c(4,3,0,n)*x4y1*6.+c(4,4,0,n)*
     & x4y2*12.+c(4,5,0,n)*x4y3*20.+c(4,6,0,n)*x4y4*30.+c(5,2,0,n)*x5*
     & 2.+c(5,3,0,n)*x5y1*6.+c(5,4,0,n)*x5y2*12.+c(5,5,0,n)*x5y3*20.+
     & c(5,6,0,n)*x5y4*30.+c(6,2,0,n)*x6*2.+c(6,3,0,n)*x6y1*6.+c(6,4,
     & 0,n)*x6y2*12.+c(6,5,0,n)*x6y3*20.+c(6,6,0,n)*x6y4*30.)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.eq.0.and.dy.eq.3.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      y1=ya(i1,i2,i3)
      y2=y1*y1
      y3=y2*y1
      x1=xa(i1,i2,i3)
      x1y1=x1*y1
      x1y2=x1y1*y1
      x1y3=x1y2*y1
      x2=x1*x1
      x2y1=x1y1*x1
      x2y2=x1y2*x1
      x2y3=x1y3*x1
      x3=x2*x1
      x3y1=x2y1*x1
      x3y2=x2y2*x1
      x3y3=x2y3*x1
      x4=x3*x1
      x4y1=x3y1*x1
      x4y2=x3y2*x1
      x4y3=x3y3*x1
      x5=x4*x1
      x5y1=x4y1*x1
      x5y2=x4y2*x1
      x5y3=x4y3*x1
      x6=x5*x1
      x6y1=x5y1*x1
      x6y2=x5y2*x1
      x6y3=x5y3*x1
      r(i1,i2,i3,n)=(c(0,3,0,n)*6.+c(0,4,0,n)*y1*24.+c(0,5,0,n)*y2*60.+
     & c(0,6,0,n)*y3*120.+c(1,3,0,n)*x1*6.+c(1,4,0,n)*x1y1*24.+c(1,5,
     & 0,n)*x1y2*60.+c(1,6,0,n)*x1y3*120.+c(2,3,0,n)*x2*6.+c(2,4,0,n)*
     & x2y1*24.+c(2,5,0,n)*x2y2*60.+c(2,6,0,n)*x2y3*120.+c(3,3,0,n)*
     & x3*6.+c(3,4,0,n)*x3y1*24.+c(3,5,0,n)*x3y2*60.+c(3,6,0,n)*x3y3*
     & 120.+c(4,3,0,n)*x4*6.+c(4,4,0,n)*x4y1*24.+c(4,5,0,n)*x4y2*60.+
     & c(4,6,0,n)*x4y3*120.+c(5,3,0,n)*x5*6.+c(5,4,0,n)*x5y1*24.+c(5,
     & 5,0,n)*x5y2*60.+c(5,6,0,n)*x5y3*120.+c(6,3,0,n)*x6*6.+c(6,4,0,
     & n)*x6y1*24.+c(6,5,0,n)*x6y2*60.+c(6,6,0,n)*x6y3*120.)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.eq.0.and.dy.eq.4.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      y1=ya(i1,i2,i3)
      y2=y1*y1
      x1=xa(i1,i2,i3)
      x1y1=x1*y1
      x1y2=x1y1*y1
      x2=x1*x1
      x2y1=x1y1*x1
      x2y2=x1y2*x1
      x3=x2*x1
      x3y1=x2y1*x1
      x3y2=x2y2*x1
      x4=x3*x1
      x4y1=x3y1*x1
      x4y2=x3y2*x1
      x5=x4*x1
      x5y1=x4y1*x1
      x5y2=x4y2*x1
      x6=x5*x1
      x6y1=x5y1*x1
      x6y2=x5y2*x1
      r(i1,i2,i3,n)=(c(0,4,0,n)*24.+c(0,5,0,n)*y1*120.+c(0,6,0,n)*y2*
     & 360.+c(1,4,0,n)*x1*24.+c(1,5,0,n)*x1y1*120.+c(1,6,0,n)*x1y2*
     & 360.+c(2,4,0,n)*x2*24.+c(2,5,0,n)*x2y1*120.+c(2,6,0,n)*x2y2*
     & 360.+c(3,4,0,n)*x3*24.+c(3,5,0,n)*x3y1*120.+c(3,6,0,n)*x3y2*
     & 360.+c(4,4,0,n)*x4*24.+c(4,5,0,n)*x4y1*120.+c(4,6,0,n)*x4y2*
     & 360.+c(5,4,0,n)*x5*24.+c(5,5,0,n)*x5y1*120.+c(5,6,0,n)*x5y2*
     & 360.+c(6,4,0,n)*x6*24.+c(6,5,0,n)*x6y1*120.+c(6,6,0,n)*x6y2*
     & 360.)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.eq.0.and.dy.eq.5.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      y1=ya(i1,i2,i3)
      x1=xa(i1,i2,i3)
      x1y1=x1*y1
      x2=x1*x1
      x2y1=x1y1*x1
      x3=x2*x1
      x3y1=x2y1*x1
      x4=x3*x1
      x4y1=x3y1*x1
      x5=x4*x1
      x5y1=x4y1*x1
      x6=x5*x1
      x6y1=x5y1*x1
      r(i1,i2,i3,n)=(c(0,5,0,n)*120.+c(0,6,0,n)*y1*720.+c(1,5,0,n)*x1*
     & 120.+c(1,6,0,n)*x1y1*720.+c(2,5,0,n)*x2*120.+c(2,6,0,n)*x2y1*
     & 720.+c(3,5,0,n)*x3*120.+c(3,6,0,n)*x3y1*720.+c(4,5,0,n)*x4*
     & 120.+c(4,6,0,n)*x4y1*720.+c(5,5,0,n)*x5*120.+c(5,6,0,n)*x5y1*
     & 720.+c(6,5,0,n)*x6*120.+c(6,6,0,n)*x6y1*720.)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.eq.0.and.dy.eq.6.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      x1=xa(i1,i2,i3)
      x2=x1*x1
      x3=x2*x1
      x4=x3*x1
      x5=x4*x1
      x6=x5*x1
      r(i1,i2,i3,n)=(c(0,6,0,n)*720.+c(1,6,0,n)*x1*720.+c(2,6,0,n)*x2*
     & 720.+c(3,6,0,n)*x3*720.+c(4,6,0,n)*x4*720.+c(5,6,0,n)*x5*720.+
     & c(6,6,0,n)*x6*720.)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.eq.1.and.dy.eq.0.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      y1=ya(i1,i2,i3)
      y2=y1*y1
      y3=y2*y1
      y4=y3*y1
      y5=y4*y1
      y6=y5*y1
      x1=xa(i1,i2,i3)
      x1y1=x1*y1
      x1y2=x1y1*y1
      x1y3=x1y2*y1
      x1y4=x1y3*y1
      x1y5=x1y4*y1
      x1y6=x1y5*y1
      x2=x1*x1
      x2y1=x1y1*x1
      x2y2=x1y2*x1
      x2y3=x1y3*x1
      x2y4=x1y4*x1
      x2y5=x1y5*x1
      x2y6=x1y6*x1
      x3=x2*x1
      x3y1=x2y1*x1
      x3y2=x2y2*x1
      x3y3=x2y3*x1
      x3y4=x2y4*x1
      x3y5=x2y5*x1
      x3y6=x2y6*x1
      x4=x3*x1
      x4y1=x3y1*x1
      x4y2=x3y2*x1
      x4y3=x3y3*x1
      x4y4=x3y4*x1
      x4y5=x3y5*x1
      x4y6=x3y6*x1
      x5=x4*x1
      x5y1=x4y1*x1
      x5y2=x4y2*x1
      x5y3=x4y3*x1
      x5y4=x4y4*x1
      x5y5=x4y5*x1
      x5y6=x4y6*x1
      r(i1,i2,i3,n)=(c(1,0,0,n)+c(1,1,0,n)*y1+c(1,2,0,n)*y2+c(1,3,0,n)*
     & y3+c(1,4,0,n)*y4+c(1,5,0,n)*y5+c(1,6,0,n)*y6+c(2,0,0,n)*x1*2.+
     & c(2,1,0,n)*x1y1*2.+c(2,2,0,n)*x1y2*2.+c(2,3,0,n)*x1y3*2.+c(2,4,
     & 0,n)*x1y4*2.+c(2,5,0,n)*x1y5*2.+c(2,6,0,n)*x1y6*2.+c(3,0,0,n)*
     & x2*3.+c(3,1,0,n)*x2y1*3.+c(3,2,0,n)*x2y2*3.+c(3,3,0,n)*x2y3*3.+
     & c(3,4,0,n)*x2y4*3.+c(3,5,0,n)*x2y5*3.+c(3,6,0,n)*x2y6*3.+c(4,0,
     & 0,n)*x3*4.+c(4,1,0,n)*x3y1*4.+c(4,2,0,n)*x3y2*4.+c(4,3,0,n)*
     & x3y3*4.+c(4,4,0,n)*x3y4*4.+c(4,5,0,n)*x3y5*4.+c(4,6,0,n)*x3y6*
     & 4.+c(5,0,0,n)*x4*5.+c(5,1,0,n)*x4y1*5.+c(5,2,0,n)*x4y2*5.+c(5,
     & 3,0,n)*x4y3*5.+c(5,4,0,n)*x4y4*5.+c(5,5,0,n)*x4y5*5.+c(5,6,0,n)
     & *x4y6*5.+c(6,0,0,n)*x5*6.+c(6,1,0,n)*x5y1*6.+c(6,2,0,n)*x5y2*
     & 6.+c(6,3,0,n)*x5y3*6.+c(6,4,0,n)*x5y4*6.+c(6,5,0,n)*x5y5*6.+c(
     & 6,6,0,n)*x5y6*6.)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.eq.1.and.dy.eq.1.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      y1=ya(i1,i2,i3)
      y2=y1*y1
      y3=y2*y1
      y4=y3*y1
      y5=y4*y1
      x1=xa(i1,i2,i3)
      x1y1=x1*y1
      x1y2=x1y1*y1
      x1y3=x1y2*y1
      x1y4=x1y3*y1
      x1y5=x1y4*y1
      x2=x1*x1
      x2y1=x1y1*x1
      x2y2=x1y2*x1
      x2y3=x1y3*x1
      x2y4=x1y4*x1
      x2y5=x1y5*x1
      x3=x2*x1
      x3y1=x2y1*x1
      x3y2=x2y2*x1
      x3y3=x2y3*x1
      x3y4=x2y4*x1
      x3y5=x2y5*x1
      x4=x3*x1
      x4y1=x3y1*x1
      x4y2=x3y2*x1
      x4y3=x3y3*x1
      x4y4=x3y4*x1
      x4y5=x3y5*x1
      x5=x4*x1
      x5y1=x4y1*x1
      x5y2=x4y2*x1
      x5y3=x4y3*x1
      x5y4=x4y4*x1
      x5y5=x4y5*x1
      r(i1,i2,i3,n)=(c(1,1,0,n)+c(1,2,0,n)*y1*2.+c(1,3,0,n)*y2*3.+c(1,
     & 4,0,n)*y3*4.+c(1,5,0,n)*y4*5.+c(1,6,0,n)*y5*6.+c(2,1,0,n)*x1*
     & 2.+c(2,2,0,n)*x1y1*4.+c(2,3,0,n)*x1y2*6.+c(2,4,0,n)*x1y3*8.+c(
     & 2,5,0,n)*x1y4*10.+c(2,6,0,n)*x1y5*12.+c(3,1,0,n)*x2*3.+c(3,2,0,
     & n)*x2y1*6.+c(3,3,0,n)*x2y2*9.+c(3,4,0,n)*x2y3*12.+c(3,5,0,n)*
     & x2y4*15.+c(3,6,0,n)*x2y5*18.+c(4,1,0,n)*x3*4.+c(4,2,0,n)*x3y1*
     & 8.+c(4,3,0,n)*x3y2*12.+c(4,4,0,n)*x3y3*16.+c(4,5,0,n)*x3y4*20.+
     & c(4,6,0,n)*x3y5*24.+c(5,1,0,n)*x4*5.+c(5,2,0,n)*x4y1*10.+c(5,3,
     & 0,n)*x4y2*15.+c(5,4,0,n)*x4y3*20.+c(5,5,0,n)*x4y4*25.+c(5,6,0,
     & n)*x4y5*30.+c(6,1,0,n)*x5*6.+c(6,2,0,n)*x5y1*12.+c(6,3,0,n)*
     & x5y2*18.+c(6,4,0,n)*x5y3*24.+c(6,5,0,n)*x5y4*30.+c(6,6,0,n)*
     & x5y5*36.)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.eq.1.and.dy.eq.2.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      y1=ya(i1,i2,i3)
      y2=y1*y1
      y3=y2*y1
      y4=y3*y1
      x1=xa(i1,i2,i3)
      x1y1=x1*y1
      x1y2=x1y1*y1
      x1y3=x1y2*y1
      x1y4=x1y3*y1
      x2=x1*x1
      x2y1=x1y1*x1
      x2y2=x1y2*x1
      x2y3=x1y3*x1
      x2y4=x1y4*x1
      x3=x2*x1
      x3y1=x2y1*x1
      x3y2=x2y2*x1
      x3y3=x2y3*x1
      x3y4=x2y4*x1
      x4=x3*x1
      x4y1=x3y1*x1
      x4y2=x3y2*x1
      x4y3=x3y3*x1
      x4y4=x3y4*x1
      x5=x4*x1
      x5y1=x4y1*x1
      x5y2=x4y2*x1
      x5y3=x4y3*x1
      x5y4=x4y4*x1
      r(i1,i2,i3,n)=(c(1,2,0,n)*2.+c(1,3,0,n)*y1*6.+c(1,4,0,n)*y2*12.+
     & c(1,5,0,n)*y3*20.+c(1,6,0,n)*y4*30.+c(2,2,0,n)*x1*4.+c(2,3,0,n)
     & *x1y1*12.+c(2,4,0,n)*x1y2*24.+c(2,5,0,n)*x1y3*40.+c(2,6,0,n)*
     & x1y4*60.+c(3,2,0,n)*x2*6.+c(3,3,0,n)*x2y1*18.+c(3,4,0,n)*x2y2*
     & 36.+c(3,5,0,n)*x2y3*60.+c(3,6,0,n)*x2y4*90.+c(4,2,0,n)*x3*8.+c(
     & 4,3,0,n)*x3y1*24.+c(4,4,0,n)*x3y2*48.+c(4,5,0,n)*x3y3*80.+c(4,
     & 6,0,n)*x3y4*120.+c(5,2,0,n)*x4*10.+c(5,3,0,n)*x4y1*30.+c(5,4,0,
     & n)*x4y2*60.+c(5,5,0,n)*x4y3*100.+c(5,6,0,n)*x4y4*150.+c(6,2,0,
     & n)*x5*12.+c(6,3,0,n)*x5y1*36.+c(6,4,0,n)*x5y2*72.+c(6,5,0,n)*
     & x5y3*120.+c(6,6,0,n)*x5y4*180.)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.eq.1.and.dy.eq.3.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      y1=ya(i1,i2,i3)
      y2=y1*y1
      y3=y2*y1
      x1=xa(i1,i2,i3)
      x1y1=x1*y1
      x1y2=x1y1*y1
      x1y3=x1y2*y1
      x2=x1*x1
      x2y1=x1y1*x1
      x2y2=x1y2*x1
      x2y3=x1y3*x1
      x3=x2*x1
      x3y1=x2y1*x1
      x3y2=x2y2*x1
      x3y3=x2y3*x1
      x4=x3*x1
      x4y1=x3y1*x1
      x4y2=x3y2*x1
      x4y3=x3y3*x1
      x5=x4*x1
      x5y1=x4y1*x1
      x5y2=x4y2*x1
      x5y3=x4y3*x1
      r(i1,i2,i3,n)=(c(1,3,0,n)*6.+c(1,4,0,n)*y1*24.+c(1,5,0,n)*y2*60.+
     & c(1,6,0,n)*y3*120.+c(2,3,0,n)*x1*12.+c(2,4,0,n)*x1y1*48.+c(2,5,
     & 0,n)*x1y2*120.+c(2,6,0,n)*x1y3*240.+c(3,3,0,n)*x2*18.+c(3,4,0,
     & n)*x2y1*72.+c(3,5,0,n)*x2y2*180.+c(3,6,0,n)*x2y3*360.+c(4,3,0,
     & n)*x3*24.+c(4,4,0,n)*x3y1*96.+c(4,5,0,n)*x3y2*240.+c(4,6,0,n)*
     & x3y3*480.+c(5,3,0,n)*x4*30.+c(5,4,0,n)*x4y1*120.+c(5,5,0,n)*
     & x4y2*300.+c(5,6,0,n)*x4y3*600.+c(6,3,0,n)*x5*36.+c(6,4,0,n)*
     & x5y1*144.+c(6,5,0,n)*x5y2*360.+c(6,6,0,n)*x5y3*720.)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.eq.1.and.dy.eq.4.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      y1=ya(i1,i2,i3)
      y2=y1*y1
      x1=xa(i1,i2,i3)
      x1y1=x1*y1
      x1y2=x1y1*y1
      x2=x1*x1
      x2y1=x1y1*x1
      x2y2=x1y2*x1
      x3=x2*x1
      x3y1=x2y1*x1
      x3y2=x2y2*x1
      x4=x3*x1
      x4y1=x3y1*x1
      x4y2=x3y2*x1
      x5=x4*x1
      x5y1=x4y1*x1
      x5y2=x4y2*x1
      r(i1,i2,i3,n)=(c(1,4,0,n)*24.+c(1,5,0,n)*y1*120.+c(1,6,0,n)*y2*
     & 360.+c(2,4,0,n)*x1*48.+c(2,5,0,n)*x1y1*240.+c(2,6,0,n)*x1y2*
     & 720.+c(3,4,0,n)*x2*72.+c(3,5,0,n)*x2y1*360.+c(3,6,0,n)*x2y2*
     & 1080.+c(4,4,0,n)*x3*96.+c(4,5,0,n)*x3y1*480.+c(4,6,0,n)*x3y2*
     & 1440.+c(5,4,0,n)*x4*120.+c(5,5,0,n)*x4y1*600.+c(5,6,0,n)*x4y2*
     & 1800.+c(6,4,0,n)*x5*144.+c(6,5,0,n)*x5y1*720.+c(6,6,0,n)*x5y2*
     & 2160.)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.eq.1.and.dy.eq.5.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      y1=ya(i1,i2,i3)
      x1=xa(i1,i2,i3)
      x1y1=x1*y1
      x2=x1*x1
      x2y1=x1y1*x1
      x3=x2*x1
      x3y1=x2y1*x1
      x4=x3*x1
      x4y1=x3y1*x1
      x5=x4*x1
      x5y1=x4y1*x1
      r(i1,i2,i3,n)=(c(1,5,0,n)*120.+c(1,6,0,n)*y1*720.+c(2,5,0,n)*x1*
     & 240.+c(2,6,0,n)*x1y1*1440.+c(3,5,0,n)*x2*360.+c(3,6,0,n)*x2y1*
     & 2160.+c(4,5,0,n)*x3*480.+c(4,6,0,n)*x3y1*2880.+c(5,5,0,n)*x4*
     & 600.+c(5,6,0,n)*x4y1*3600.+c(6,5,0,n)*x5*720.+c(6,6,0,n)*x5y1*
     & 4320.)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.eq.1.and.dy.eq.6.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      x1=xa(i1,i2,i3)
      x2=x1*x1
      x3=x2*x1
      x4=x3*x1
      x5=x4*x1
      r(i1,i2,i3,n)=(c(1,6,0,n)*720.+c(2,6,0,n)*x1*1440.+c(3,6,0,n)*x2*
     & 2160.+c(4,6,0,n)*x3*2880.+c(5,6,0,n)*x4*3600.+c(6,6,0,n)*x5*
     & 4320.)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.eq.2.and.dy.eq.0.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      y1=ya(i1,i2,i3)
      y2=y1*y1
      y3=y2*y1
      y4=y3*y1
      y5=y4*y1
      y6=y5*y1
      x1=xa(i1,i2,i3)
      x1y1=x1*y1
      x1y2=x1y1*y1
      x1y3=x1y2*y1
      x1y4=x1y3*y1
      x1y5=x1y4*y1
      x1y6=x1y5*y1
      x2=x1*x1
      x2y1=x1y1*x1
      x2y2=x1y2*x1
      x2y3=x1y3*x1
      x2y4=x1y4*x1
      x2y5=x1y5*x1
      x2y6=x1y6*x1
      x3=x2*x1
      x3y1=x2y1*x1
      x3y2=x2y2*x1
      x3y3=x2y3*x1
      x3y4=x2y4*x1
      x3y5=x2y5*x1
      x3y6=x2y6*x1
      x4=x3*x1
      x4y1=x3y1*x1
      x4y2=x3y2*x1
      x4y3=x3y3*x1
      x4y4=x3y4*x1
      x4y5=x3y5*x1
      x4y6=x3y6*x1
      r(i1,i2,i3,n)=(c(2,0,0,n)*2.+c(2,1,0,n)*y1*2.+c(2,2,0,n)*y2*2.+c(
     & 2,3,0,n)*y3*2.+c(2,4,0,n)*y4*2.+c(2,5,0,n)*y5*2.+c(2,6,0,n)*y6*
     & 2.+c(3,0,0,n)*x1*6.+c(3,1,0,n)*x1y1*6.+c(3,2,0,n)*x1y2*6.+c(3,
     & 3,0,n)*x1y3*6.+c(3,4,0,n)*x1y4*6.+c(3,5,0,n)*x1y5*6.+c(3,6,0,n)
     & *x1y6*6.+c(4,0,0,n)*x2*12.+c(4,1,0,n)*x2y1*12.+c(4,2,0,n)*x2y2*
     & 12.+c(4,3,0,n)*x2y3*12.+c(4,4,0,n)*x2y4*12.+c(4,5,0,n)*x2y5*
     & 12.+c(4,6,0,n)*x2y6*12.+c(5,0,0,n)*x3*20.+c(5,1,0,n)*x3y1*20.+
     & c(5,2,0,n)*x3y2*20.+c(5,3,0,n)*x3y3*20.+c(5,4,0,n)*x3y4*20.+c(
     & 5,5,0,n)*x3y5*20.+c(5,6,0,n)*x3y6*20.+c(6,0,0,n)*x4*30.+c(6,1,
     & 0,n)*x4y1*30.+c(6,2,0,n)*x4y2*30.+c(6,3,0,n)*x4y3*30.+c(6,4,0,
     & n)*x4y4*30.+c(6,5,0,n)*x4y5*30.+c(6,6,0,n)*x4y6*30.)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.eq.2.and.dy.eq.1.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      y1=ya(i1,i2,i3)
      y2=y1*y1
      y3=y2*y1
      y4=y3*y1
      y5=y4*y1
      x1=xa(i1,i2,i3)
      x1y1=x1*y1
      x1y2=x1y1*y1
      x1y3=x1y2*y1
      x1y4=x1y3*y1
      x1y5=x1y4*y1
      x2=x1*x1
      x2y1=x1y1*x1
      x2y2=x1y2*x1
      x2y3=x1y3*x1
      x2y4=x1y4*x1
      x2y5=x1y5*x1
      x3=x2*x1
      x3y1=x2y1*x1
      x3y2=x2y2*x1
      x3y3=x2y3*x1
      x3y4=x2y4*x1
      x3y5=x2y5*x1
      x4=x3*x1
      x4y1=x3y1*x1
      x4y2=x3y2*x1
      x4y3=x3y3*x1
      x4y4=x3y4*x1
      x4y5=x3y5*x1
      r(i1,i2,i3,n)=(c(2,1,0,n)*2.+c(2,2,0,n)*y1*4.+c(2,3,0,n)*y2*6.+c(
     & 2,4,0,n)*y3*8.+c(2,5,0,n)*y4*10.+c(2,6,0,n)*y5*12.+c(3,1,0,n)*
     & x1*6.+c(3,2,0,n)*x1y1*12.+c(3,3,0,n)*x1y2*18.+c(3,4,0,n)*x1y3*
     & 24.+c(3,5,0,n)*x1y4*30.+c(3,6,0,n)*x1y5*36.+c(4,1,0,n)*x2*12.+
     & c(4,2,0,n)*x2y1*24.+c(4,3,0,n)*x2y2*36.+c(4,4,0,n)*x2y3*48.+c(
     & 4,5,0,n)*x2y4*60.+c(4,6,0,n)*x2y5*72.+c(5,1,0,n)*x3*20.+c(5,2,
     & 0,n)*x3y1*40.+c(5,3,0,n)*x3y2*60.+c(5,4,0,n)*x3y3*80.+c(5,5,0,
     & n)*x3y4*100.+c(5,6,0,n)*x3y5*120.+c(6,1,0,n)*x4*30.+c(6,2,0,n)*
     & x4y1*60.+c(6,3,0,n)*x4y2*90.+c(6,4,0,n)*x4y3*120.+c(6,5,0,n)*
     & x4y4*150.+c(6,6,0,n)*x4y5*180.)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.eq.2.and.dy.eq.2.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      y1=ya(i1,i2,i3)
      y2=y1*y1
      y3=y2*y1
      y4=y3*y1
      x1=xa(i1,i2,i3)
      x1y1=x1*y1
      x1y2=x1y1*y1
      x1y3=x1y2*y1
      x1y4=x1y3*y1
      x2=x1*x1
      x2y1=x1y1*x1
      x2y2=x1y2*x1
      x2y3=x1y3*x1
      x2y4=x1y4*x1
      x3=x2*x1
      x3y1=x2y1*x1
      x3y2=x2y2*x1
      x3y3=x2y3*x1
      x3y4=x2y4*x1
      x4=x3*x1
      x4y1=x3y1*x1
      x4y2=x3y2*x1
      x4y3=x3y3*x1
      x4y4=x3y4*x1
      r(i1,i2,i3,n)=(c(2,2,0,n)*4.+c(2,3,0,n)*y1*12.+c(2,4,0,n)*y2*24.+
     & c(2,5,0,n)*y3*40.+c(2,6,0,n)*y4*60.+c(3,2,0,n)*x1*12.+c(3,3,0,
     & n)*x1y1*36.+c(3,4,0,n)*x1y2*72.+c(3,5,0,n)*x1y3*120.+c(3,6,0,n)
     & *x1y4*180.+c(4,2,0,n)*x2*24.+c(4,3,0,n)*x2y1*72.+c(4,4,0,n)*
     & x2y2*144.+c(4,5,0,n)*x2y3*240.+c(4,6,0,n)*x2y4*360.+c(5,2,0,n)*
     & x3*40.+c(5,3,0,n)*x3y1*120.+c(5,4,0,n)*x3y2*240.+c(5,5,0,n)*
     & x3y3*400.+c(5,6,0,n)*x3y4*600.+c(6,2,0,n)*x4*60.+c(6,3,0,n)*
     & x4y1*180.+c(6,4,0,n)*x4y2*360.+c(6,5,0,n)*x4y3*600.+c(6,6,0,n)*
     & x4y4*900.)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.eq.2.and.dy.eq.3.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      y1=ya(i1,i2,i3)
      y2=y1*y1
      y3=y2*y1
      x1=xa(i1,i2,i3)
      x1y1=x1*y1
      x1y2=x1y1*y1
      x1y3=x1y2*y1
      x2=x1*x1
      x2y1=x1y1*x1
      x2y2=x1y2*x1
      x2y3=x1y3*x1
      x3=x2*x1
      x3y1=x2y1*x1
      x3y2=x2y2*x1
      x3y3=x2y3*x1
      x4=x3*x1
      x4y1=x3y1*x1
      x4y2=x3y2*x1
      x4y3=x3y3*x1
      r(i1,i2,i3,n)=(c(2,3,0,n)*12.+c(2,4,0,n)*y1*48.+c(2,5,0,n)*y2*
     & 120.+c(2,6,0,n)*y3*240.+c(3,3,0,n)*x1*36.+c(3,4,0,n)*x1y1*144.+
     & c(3,5,0,n)*x1y2*360.+c(3,6,0,n)*x1y3*720.+c(4,3,0,n)*x2*72.+c(
     & 4,4,0,n)*x2y1*288.+c(4,5,0,n)*x2y2*720.+c(4,6,0,n)*x2y3*1440.+
     & c(5,3,0,n)*x3*120.+c(5,4,0,n)*x3y1*480.+c(5,5,0,n)*x3y2*1200.+
     & c(5,6,0,n)*x3y3*2400.+c(6,3,0,n)*x4*180.+c(6,4,0,n)*x4y1*720.+
     & c(6,5,0,n)*x4y2*1800.+c(6,6,0,n)*x4y3*3600.)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.eq.2.and.dy.eq.4.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      y1=ya(i1,i2,i3)
      y2=y1*y1
      x1=xa(i1,i2,i3)
      x1y1=x1*y1
      x1y2=x1y1*y1
      x2=x1*x1
      x2y1=x1y1*x1
      x2y2=x1y2*x1
      x3=x2*x1
      x3y1=x2y1*x1
      x3y2=x2y2*x1
      x4=x3*x1
      x4y1=x3y1*x1
      x4y2=x3y2*x1
      r(i1,i2,i3,n)=(c(2,4,0,n)*48.+c(2,5,0,n)*y1*240.+c(2,6,0,n)*y2*
     & 720.+c(3,4,0,n)*x1*144.+c(3,5,0,n)*x1y1*720.+c(3,6,0,n)*x1y2*
     & 2160.+c(4,4,0,n)*x2*288.+c(4,5,0,n)*x2y1*1440.+c(4,6,0,n)*x2y2*
     & 4320.+c(5,4,0,n)*x3*480.+c(5,5,0,n)*x3y1*2400.+c(5,6,0,n)*x3y2*
     & 7200.+c(6,4,0,n)*x4*720.+c(6,5,0,n)*x4y1*3600.+c(6,6,0,n)*x4y2*
     & 10800.)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.eq.2.and.dy.eq.5.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      y1=ya(i1,i2,i3)
      x1=xa(i1,i2,i3)
      x1y1=x1*y1
      x2=x1*x1
      x2y1=x1y1*x1
      x3=x2*x1
      x3y1=x2y1*x1
      x4=x3*x1
      x4y1=x3y1*x1
      r(i1,i2,i3,n)=(c(2,5,0,n)*240.+c(2,6,0,n)*y1*1440.+c(3,5,0,n)*x1*
     & 720.+c(3,6,0,n)*x1y1*4320.+c(4,5,0,n)*x2*1440.+c(4,6,0,n)*x2y1*
     & 8640.+c(5,5,0,n)*x3*2400.+c(5,6,0,n)*x3y1*14400.+c(6,5,0,n)*x4*
     & 3600.+c(6,6,0,n)*x4y1*21600.)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.eq.2.and.dy.eq.6.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      x1=xa(i1,i2,i3)
      x2=x1*x1
      x3=x2*x1
      x4=x3*x1
      r(i1,i2,i3,n)=(c(2,6,0,n)*1440.+c(3,6,0,n)*x1*4320.+c(4,6,0,n)*
     & x2*8640.+c(5,6,0,n)*x3*14400.+c(6,6,0,n)*x4*21600.)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.eq.3.and.dy.eq.0.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      y1=ya(i1,i2,i3)
      y2=y1*y1
      y3=y2*y1
      y4=y3*y1
      y5=y4*y1
      y6=y5*y1
      x1=xa(i1,i2,i3)
      x1y1=x1*y1
      x1y2=x1y1*y1
      x1y3=x1y2*y1
      x1y4=x1y3*y1
      x1y5=x1y4*y1
      x1y6=x1y5*y1
      x2=x1*x1
      x2y1=x1y1*x1
      x2y2=x1y2*x1
      x2y3=x1y3*x1
      x2y4=x1y4*x1
      x2y5=x1y5*x1
      x2y6=x1y6*x1
      x3=x2*x1
      x3y1=x2y1*x1
      x3y2=x2y2*x1
      x3y3=x2y3*x1
      x3y4=x2y4*x1
      x3y5=x2y5*x1
      x3y6=x2y6*x1
      r(i1,i2,i3,n)=(c(3,0,0,n)*6.+c(3,1,0,n)*y1*6.+c(3,2,0,n)*y2*6.+c(
     & 3,3,0,n)*y3*6.+c(3,4,0,n)*y4*6.+c(3,5,0,n)*y5*6.+c(3,6,0,n)*y6*
     & 6.+c(4,0,0,n)*x1*24.+c(4,1,0,n)*x1y1*24.+c(4,2,0,n)*x1y2*24.+c(
     & 4,3,0,n)*x1y3*24.+c(4,4,0,n)*x1y4*24.+c(4,5,0,n)*x1y5*24.+c(4,
     & 6,0,n)*x1y6*24.+c(5,0,0,n)*x2*60.+c(5,1,0,n)*x2y1*60.+c(5,2,0,
     & n)*x2y2*60.+c(5,3,0,n)*x2y3*60.+c(5,4,0,n)*x2y4*60.+c(5,5,0,n)*
     & x2y5*60.+c(5,6,0,n)*x2y6*60.+c(6,0,0,n)*x3*120.+c(6,1,0,n)*
     & x3y1*120.+c(6,2,0,n)*x3y2*120.+c(6,3,0,n)*x3y3*120.+c(6,4,0,n)*
     & x3y4*120.+c(6,5,0,n)*x3y5*120.+c(6,6,0,n)*x3y6*120.)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.eq.3.and.dy.eq.1.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      y1=ya(i1,i2,i3)
      y2=y1*y1
      y3=y2*y1
      y4=y3*y1
      y5=y4*y1
      x1=xa(i1,i2,i3)
      x1y1=x1*y1
      x1y2=x1y1*y1
      x1y3=x1y2*y1
      x1y4=x1y3*y1
      x1y5=x1y4*y1
      x2=x1*x1
      x2y1=x1y1*x1
      x2y2=x1y2*x1
      x2y3=x1y3*x1
      x2y4=x1y4*x1
      x2y5=x1y5*x1
      x3=x2*x1
      x3y1=x2y1*x1
      x3y2=x2y2*x1
      x3y3=x2y3*x1
      x3y4=x2y4*x1
      x3y5=x2y5*x1
      r(i1,i2,i3,n)=(c(3,1,0,n)*6.+c(3,2,0,n)*y1*12.+c(3,3,0,n)*y2*18.+
     & c(3,4,0,n)*y3*24.+c(3,5,0,n)*y4*30.+c(3,6,0,n)*y5*36.+c(4,1,0,
     & n)*x1*24.+c(4,2,0,n)*x1y1*48.+c(4,3,0,n)*x1y2*72.+c(4,4,0,n)*
     & x1y3*96.+c(4,5,0,n)*x1y4*120.+c(4,6,0,n)*x1y5*144.+c(5,1,0,n)*
     & x2*60.+c(5,2,0,n)*x2y1*120.+c(5,3,0,n)*x2y2*180.+c(5,4,0,n)*
     & x2y3*240.+c(5,5,0,n)*x2y4*300.+c(5,6,0,n)*x2y5*360.+c(6,1,0,n)*
     & x3*120.+c(6,2,0,n)*x3y1*240.+c(6,3,0,n)*x3y2*360.+c(6,4,0,n)*
     & x3y3*480.+c(6,5,0,n)*x3y4*600.+c(6,6,0,n)*x3y5*720.)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.eq.3.and.dy.eq.2.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      y1=ya(i1,i2,i3)
      y2=y1*y1
      y3=y2*y1
      y4=y3*y1
      x1=xa(i1,i2,i3)
      x1y1=x1*y1
      x1y2=x1y1*y1
      x1y3=x1y2*y1
      x1y4=x1y3*y1
      x2=x1*x1
      x2y1=x1y1*x1
      x2y2=x1y2*x1
      x2y3=x1y3*x1
      x2y4=x1y4*x1
      x3=x2*x1
      x3y1=x2y1*x1
      x3y2=x2y2*x1
      x3y3=x2y3*x1
      x3y4=x2y4*x1
      r(i1,i2,i3,n)=(c(3,2,0,n)*12.+c(3,3,0,n)*y1*36.+c(3,4,0,n)*y2*
     & 72.+c(3,5,0,n)*y3*120.+c(3,6,0,n)*y4*180.+c(4,2,0,n)*x1*48.+c(
     & 4,3,0,n)*x1y1*144.+c(4,4,0,n)*x1y2*288.+c(4,5,0,n)*x1y3*480.+c(
     & 4,6,0,n)*x1y4*720.+c(5,2,0,n)*x2*120.+c(5,3,0,n)*x2y1*360.+c(5,
     & 4,0,n)*x2y2*720.+c(5,5,0,n)*x2y3*1200.+c(5,6,0,n)*x2y4*1800.+c(
     & 6,2,0,n)*x3*240.+c(6,3,0,n)*x3y1*720.+c(6,4,0,n)*x3y2*1440.+c(
     & 6,5,0,n)*x3y3*2400.+c(6,6,0,n)*x3y4*3600.)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.eq.3.and.dy.eq.3.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      y1=ya(i1,i2,i3)
      y2=y1*y1
      y3=y2*y1
      x1=xa(i1,i2,i3)
      x1y1=x1*y1
      x1y2=x1y1*y1
      x1y3=x1y2*y1
      x2=x1*x1
      x2y1=x1y1*x1
      x2y2=x1y2*x1
      x2y3=x1y3*x1
      x3=x2*x1
      x3y1=x2y1*x1
      x3y2=x2y2*x1
      x3y3=x2y3*x1
      r(i1,i2,i3,n)=(c(3,3,0,n)*36.+c(3,4,0,n)*y1*144.+c(3,5,0,n)*y2*
     & 360.+c(3,6,0,n)*y3*720.+c(4,3,0,n)*x1*144.+c(4,4,0,n)*x1y1*
     & 576.+c(4,5,0,n)*x1y2*1440.+c(4,6,0,n)*x1y3*2880.+c(5,3,0,n)*x2*
     & 360.+c(5,4,0,n)*x2y1*1440.+c(5,5,0,n)*x2y2*3600.+c(5,6,0,n)*
     & x2y3*7200.+c(6,3,0,n)*x3*720.+c(6,4,0,n)*x3y1*2880.+c(6,5,0,n)*
     & x3y2*7200.+c(6,6,0,n)*x3y3*14400.)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.eq.3.and.dy.eq.4.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      y1=ya(i1,i2,i3)
      y2=y1*y1
      x1=xa(i1,i2,i3)
      x1y1=x1*y1
      x1y2=x1y1*y1
      x2=x1*x1
      x2y1=x1y1*x1
      x2y2=x1y2*x1
      x3=x2*x1
      x3y1=x2y1*x1
      x3y2=x2y2*x1
      r(i1,i2,i3,n)=(c(3,4,0,n)*144.+c(3,5,0,n)*y1*720.+c(3,6,0,n)*y2*
     & 2160.+c(4,4,0,n)*x1*576.+c(4,5,0,n)*x1y1*2880.+c(4,6,0,n)*x1y2*
     & 8640.+c(5,4,0,n)*x2*1440.+c(5,5,0,n)*x2y1*7200.+c(5,6,0,n)*
     & x2y2*21600.+c(6,4,0,n)*x3*2880.+c(6,5,0,n)*x3y1*14400.+c(6,6,0,
     & n)*x3y2*43200.)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.eq.3.and.dy.eq.5.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      y1=ya(i1,i2,i3)
      x1=xa(i1,i2,i3)
      x1y1=x1*y1
      x2=x1*x1
      x2y1=x1y1*x1
      x3=x2*x1
      x3y1=x2y1*x1
      r(i1,i2,i3,n)=(c(3,5,0,n)*720.+c(3,6,0,n)*y1*4320.+c(4,5,0,n)*x1*
     & 2880.+c(4,6,0,n)*x1y1*17280.+c(5,5,0,n)*x2*7200.+c(5,6,0,n)*
     & x2y1*43200.+c(6,5,0,n)*x3*14400.+c(6,6,0,n)*x3y1*86400.)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.eq.3.and.dy.eq.6.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      x1=xa(i1,i2,i3)
      x2=x1*x1
      x3=x2*x1
      r(i1,i2,i3,n)=(c(3,6,0,n)*4320.+c(4,6,0,n)*x1*17280.+c(5,6,0,n)*
     & x2*43200.+c(6,6,0,n)*x3*86400.)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.eq.4.and.dy.eq.0.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      y1=ya(i1,i2,i3)
      y2=y1*y1
      y3=y2*y1
      y4=y3*y1
      y5=y4*y1
      y6=y5*y1
      x1=xa(i1,i2,i3)
      x1y1=x1*y1
      x1y2=x1y1*y1
      x1y3=x1y2*y1
      x1y4=x1y3*y1
      x1y5=x1y4*y1
      x1y6=x1y5*y1
      x2=x1*x1
      x2y1=x1y1*x1
      x2y2=x1y2*x1
      x2y3=x1y3*x1
      x2y4=x1y4*x1
      x2y5=x1y5*x1
      x2y6=x1y6*x1
      r(i1,i2,i3,n)=(c(4,0,0,n)*24.+c(4,1,0,n)*y1*24.+c(4,2,0,n)*y2*
     & 24.+c(4,3,0,n)*y3*24.+c(4,4,0,n)*y4*24.+c(4,5,0,n)*y5*24.+c(4,
     & 6,0,n)*y6*24.+c(5,0,0,n)*x1*120.+c(5,1,0,n)*x1y1*120.+c(5,2,0,
     & n)*x1y2*120.+c(5,3,0,n)*x1y3*120.+c(5,4,0,n)*x1y4*120.+c(5,5,0,
     & n)*x1y5*120.+c(5,6,0,n)*x1y6*120.+c(6,0,0,n)*x2*360.+c(6,1,0,n)
     & *x2y1*360.+c(6,2,0,n)*x2y2*360.+c(6,3,0,n)*x2y3*360.+c(6,4,0,n)
     & *x2y4*360.+c(6,5,0,n)*x2y5*360.+c(6,6,0,n)*x2y6*360.)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.eq.4.and.dy.eq.1.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      y1=ya(i1,i2,i3)
      y2=y1*y1
      y3=y2*y1
      y4=y3*y1
      y5=y4*y1
      x1=xa(i1,i2,i3)
      x1y1=x1*y1
      x1y2=x1y1*y1
      x1y3=x1y2*y1
      x1y4=x1y3*y1
      x1y5=x1y4*y1
      x2=x1*x1
      x2y1=x1y1*x1
      x2y2=x1y2*x1
      x2y3=x1y3*x1
      x2y4=x1y4*x1
      x2y5=x1y5*x1
      r(i1,i2,i3,n)=(c(4,1,0,n)*24.+c(4,2,0,n)*y1*48.+c(4,3,0,n)*y2*
     & 72.+c(4,4,0,n)*y3*96.+c(4,5,0,n)*y4*120.+c(4,6,0,n)*y5*144.+c(
     & 5,1,0,n)*x1*120.+c(5,2,0,n)*x1y1*240.+c(5,3,0,n)*x1y2*360.+c(5,
     & 4,0,n)*x1y3*480.+c(5,5,0,n)*x1y4*600.+c(5,6,0,n)*x1y5*720.+c(6,
     & 1,0,n)*x2*360.+c(6,2,0,n)*x2y1*720.+c(6,3,0,n)*x2y2*1080.+c(6,
     & 4,0,n)*x2y3*1440.+c(6,5,0,n)*x2y4*1800.+c(6,6,0,n)*x2y5*2160.)*
     & time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.eq.4.and.dy.eq.2.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      y1=ya(i1,i2,i3)
      y2=y1*y1
      y3=y2*y1
      y4=y3*y1
      x1=xa(i1,i2,i3)
      x1y1=x1*y1
      x1y2=x1y1*y1
      x1y3=x1y2*y1
      x1y4=x1y3*y1
      x2=x1*x1
      x2y1=x1y1*x1
      x2y2=x1y2*x1
      x2y3=x1y3*x1
      x2y4=x1y4*x1
      r(i1,i2,i3,n)=(c(4,2,0,n)*48.+c(4,3,0,n)*y1*144.+c(4,4,0,n)*y2*
     & 288.+c(4,5,0,n)*y3*480.+c(4,6,0,n)*y4*720.+c(5,2,0,n)*x1*240.+
     & c(5,3,0,n)*x1y1*720.+c(5,4,0,n)*x1y2*1440.+c(5,5,0,n)*x1y3*
     & 2400.+c(5,6,0,n)*x1y4*3600.+c(6,2,0,n)*x2*720.+c(6,3,0,n)*x2y1*
     & 2160.+c(6,4,0,n)*x2y2*4320.+c(6,5,0,n)*x2y3*7200.+c(6,6,0,n)*
     & x2y4*10800.)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.eq.4.and.dy.eq.3.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      y1=ya(i1,i2,i3)
      y2=y1*y1
      y3=y2*y1
      x1=xa(i1,i2,i3)
      x1y1=x1*y1
      x1y2=x1y1*y1
      x1y3=x1y2*y1
      x2=x1*x1
      x2y1=x1y1*x1
      x2y2=x1y2*x1
      x2y3=x1y3*x1
      r(i1,i2,i3,n)=(c(4,3,0,n)*144.+c(4,4,0,n)*y1*576.+c(4,5,0,n)*y2*
     & 1440.+c(4,6,0,n)*y3*2880.+c(5,3,0,n)*x1*720.+c(5,4,0,n)*x1y1*
     & 2880.+c(5,5,0,n)*x1y2*7200.+c(5,6,0,n)*x1y3*14400.+c(6,3,0,n)*
     & x2*2160.+c(6,4,0,n)*x2y1*8640.+c(6,5,0,n)*x2y2*21600.+c(6,6,0,
     & n)*x2y3*43200.)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.eq.4.and.dy.eq.4.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      y1=ya(i1,i2,i3)
      y2=y1*y1
      x1=xa(i1,i2,i3)
      x1y1=x1*y1
      x1y2=x1y1*y1
      x2=x1*x1
      x2y1=x1y1*x1
      x2y2=x1y2*x1
      r(i1,i2,i3,n)=(c(4,4,0,n)*576.+c(4,5,0,n)*y1*2880.+c(4,6,0,n)*y2*
     & 8640.+c(5,4,0,n)*x1*2880.+c(5,5,0,n)*x1y1*14400.+c(5,6,0,n)*
     & x1y2*43200.+c(6,4,0,n)*x2*8640.+c(6,5,0,n)*x2y1*43200.+c(6,6,0,
     & n)*x2y2*129600.)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.eq.4.and.dy.eq.5.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      y1=ya(i1,i2,i3)
      x1=xa(i1,i2,i3)
      x1y1=x1*y1
      x2=x1*x1
      x2y1=x1y1*x1
      r(i1,i2,i3,n)=(c(4,5,0,n)*2880.+c(4,6,0,n)*y1*17280.+c(5,5,0,n)*
     & x1*14400.+c(5,6,0,n)*x1y1*86400.+c(6,5,0,n)*x2*43200.+c(6,6,0,
     & n)*x2y1*259200.)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.eq.4.and.dy.eq.6.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      x1=xa(i1,i2,i3)
      x2=x1*x1
      r(i1,i2,i3,n)=(c(4,6,0,n)*17280.+c(5,6,0,n)*x1*86400.+c(6,6,0,n)*
     & x2*259200.)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.eq.5.and.dy.eq.0.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      y1=ya(i1,i2,i3)
      y2=y1*y1
      y3=y2*y1
      y4=y3*y1
      y5=y4*y1
      y6=y5*y1
      x1=xa(i1,i2,i3)
      x1y1=x1*y1
      x1y2=x1y1*y1
      x1y3=x1y2*y1
      x1y4=x1y3*y1
      x1y5=x1y4*y1
      x1y6=x1y5*y1
      r(i1,i2,i3,n)=(c(5,0,0,n)*120.+c(5,1,0,n)*y1*120.+c(5,2,0,n)*y2*
     & 120.+c(5,3,0,n)*y3*120.+c(5,4,0,n)*y4*120.+c(5,5,0,n)*y5*120.+
     & c(5,6,0,n)*y6*120.+c(6,0,0,n)*x1*720.+c(6,1,0,n)*x1y1*720.+c(6,
     & 2,0,n)*x1y2*720.+c(6,3,0,n)*x1y3*720.+c(6,4,0,n)*x1y4*720.+c(6,
     & 5,0,n)*x1y5*720.+c(6,6,0,n)*x1y6*720.)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.eq.5.and.dy.eq.1.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      y1=ya(i1,i2,i3)
      y2=y1*y1
      y3=y2*y1
      y4=y3*y1
      y5=y4*y1
      x1=xa(i1,i2,i3)
      x1y1=x1*y1
      x1y2=x1y1*y1
      x1y3=x1y2*y1
      x1y4=x1y3*y1
      x1y5=x1y4*y1
      r(i1,i2,i3,n)=(c(5,1,0,n)*120.+c(5,2,0,n)*y1*240.+c(5,3,0,n)*y2*
     & 360.+c(5,4,0,n)*y3*480.+c(5,5,0,n)*y4*600.+c(5,6,0,n)*y5*720.+
     & c(6,1,0,n)*x1*720.+c(6,2,0,n)*x1y1*1440.+c(6,3,0,n)*x1y2*2160.+
     & c(6,4,0,n)*x1y3*2880.+c(6,5,0,n)*x1y4*3600.+c(6,6,0,n)*x1y5*
     & 4320.)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.eq.5.and.dy.eq.2.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      y1=ya(i1,i2,i3)
      y2=y1*y1
      y3=y2*y1
      y4=y3*y1
      x1=xa(i1,i2,i3)
      x1y1=x1*y1
      x1y2=x1y1*y1
      x1y3=x1y2*y1
      x1y4=x1y3*y1
      r(i1,i2,i3,n)=(c(5,2,0,n)*240.+c(5,3,0,n)*y1*720.+c(5,4,0,n)*y2*
     & 1440.+c(5,5,0,n)*y3*2400.+c(5,6,0,n)*y4*3600.+c(6,2,0,n)*x1*
     & 1440.+c(6,3,0,n)*x1y1*4320.+c(6,4,0,n)*x1y2*8640.+c(6,5,0,n)*
     & x1y3*14400.+c(6,6,0,n)*x1y4*21600.)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.eq.5.and.dy.eq.3.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      y1=ya(i1,i2,i3)
      y2=y1*y1
      y3=y2*y1
      x1=xa(i1,i2,i3)
      x1y1=x1*y1
      x1y2=x1y1*y1
      x1y3=x1y2*y1
      r(i1,i2,i3,n)=(c(5,3,0,n)*720.+c(5,4,0,n)*y1*2880.+c(5,5,0,n)*y2*
     & 7200.+c(5,6,0,n)*y3*14400.+c(6,3,0,n)*x1*4320.+c(6,4,0,n)*x1y1*
     & 17280.+c(6,5,0,n)*x1y2*43200.+c(6,6,0,n)*x1y3*86400.)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.eq.5.and.dy.eq.4.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      y1=ya(i1,i2,i3)
      y2=y1*y1
      x1=xa(i1,i2,i3)
      x1y1=x1*y1
      x1y2=x1y1*y1
      r(i1,i2,i3,n)=(c(5,4,0,n)*2880.+c(5,5,0,n)*y1*14400.+c(5,6,0,n)*
     & y2*43200.+c(6,4,0,n)*x1*17280.+c(6,5,0,n)*x1y1*86400.+c(6,6,0,
     & n)*x1y2*259200.)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.eq.5.and.dy.eq.5.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      y1=ya(i1,i2,i3)
      x1=xa(i1,i2,i3)
      x1y1=x1*y1
      r(i1,i2,i3,n)=(c(5,5,0,n)*14400.+c(5,6,0,n)*y1*86400.+c(6,5,0,n)*
     & x1*86400.+c(6,6,0,n)*x1y1*518400.)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.eq.5.and.dy.eq.6.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      x1=xa(i1,i2,i3)
      r(i1,i2,i3,n)=(c(5,6,0,n)*86400.+c(6,6,0,n)*x1*518400.)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.eq.6.and.dy.eq.0.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      y1=ya(i1,i2,i3)
      y2=y1*y1
      y3=y2*y1
      y4=y3*y1
      y5=y4*y1
      y6=y5*y1
      r(i1,i2,i3,n)=(c(6,0,0,n)*720.+c(6,1,0,n)*y1*720.+c(6,2,0,n)*y2*
     & 720.+c(6,3,0,n)*y3*720.+c(6,4,0,n)*y4*720.+c(6,5,0,n)*y5*720.+
     & c(6,6,0,n)*y6*720.)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.eq.6.and.dy.eq.1.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      y1=ya(i1,i2,i3)
      y2=y1*y1
      y3=y2*y1
      y4=y3*y1
      y5=y4*y1
      r(i1,i2,i3,n)=(c(6,1,0,n)*720.+c(6,2,0,n)*y1*1440.+c(6,3,0,n)*y2*
     & 2160.+c(6,4,0,n)*y3*2880.+c(6,5,0,n)*y4*3600.+c(6,6,0,n)*y5*
     & 4320.)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.eq.6.and.dy.eq.2.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      y1=ya(i1,i2,i3)
      y2=y1*y1
      y3=y2*y1
      y4=y3*y1
      r(i1,i2,i3,n)=(c(6,2,0,n)*1440.+c(6,3,0,n)*y1*4320.+c(6,4,0,n)*
     & y2*8640.+c(6,5,0,n)*y3*14400.+c(6,6,0,n)*y4*21600.)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.eq.6.and.dy.eq.3.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      y1=ya(i1,i2,i3)
      y2=y1*y1
      y3=y2*y1
      r(i1,i2,i3,n)=(c(6,3,0,n)*4320.+c(6,4,0,n)*y1*17280.+c(6,5,0,n)*
     & y2*43200.+c(6,6,0,n)*y3*86400.)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.eq.6.and.dy.eq.4.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      y1=ya(i1,i2,i3)
      y2=y1*y1
      r(i1,i2,i3,n)=(c(6,4,0,n)*17280.+c(6,5,0,n)*y1*86400.+c(6,6,0,n)*
     & y2*259200.)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.eq.6.and.dy.eq.5.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      y1=ya(i1,i2,i3)
      r(i1,i2,i3,n)=(c(6,5,0,n)*86400.+c(6,6,0,n)*y1*518400.)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.eq.6.and.dy.eq.6.and.dz.eq.0 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      r(i1,i2,i3,n)=(c(6,6,0,n)*518400.)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( laplace.eq.1 )then
! beginLoops($defineTime())
      do n=nca,ncb
! defineTime()
      if( dt.eq.0 )then
       if( degreeTime.eq.0 )then
        time=a(0,n)
       else if( degreeTime.eq.1 )then
        time=a(0,n)+t*(a(1,n))
       else if( degreeTime.eq.2 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)))
       else if( degreeTime.eq.3 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n))))
       else if( degreeTime.eq.4 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)))))
       else if( degreeTime.eq.5 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)))
     & )))
       else if( degreeTime.eq.6 )then
        time=a(0,n)+t*(a(1,n)+t*(a(2,n)+t*(a(3,n)+t*(a(4,n)+t*(a(5,n)+
     & t*(a(6,n)))))))
       else
        write(*,*) 'ERROR invalid degreeTime'
        stop
       end if
      else if( dt.eq.1 )then
       ! --- first time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=a(1,n)
       else if( degreeTime.eq.2 )then
        time=a(1,n)+t*(2.*a(2,n))
       else if( degreeTime.eq.3 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)))
       else if( degreeTime.eq.4 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n))))
       else if( degreeTime.eq.5 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & ))))
       else if( degreeTime.eq.6 )then
        time=a(1,n)+t*(2.*a(2,n)+t*(3.*a(3,n)+t*(4.*a(4,n)+t*(5.*a(5,n)
     & +t*(6.*a(6,n))))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1834
       end if
      else if( dt.eq.2 )then
       ! --- 2nd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=2.*a(2,n)
       else if( degreeTime.eq.3 )then
        time=2.*a(2,n)+t*(6.*a(3,n))
       else if( degreeTime.eq.4 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)))
       else if( degreeTime.eq.5 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n))))
       else if( degreeTime.eq.6 )then
        time=2.*a(2,n)+t*(6.*a(3,n)+t*(12.*a(4,n)+t*(20.*a(5,n)+t*(30.*
     & a(6,n)))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1835
       end if
      else if( dt.eq.3 )then
       ! --- 3rd time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=6.*a(3,n)
       else if( degreeTime.eq.4 )then
        time=6.*a(3,n)+t*(24.*a(4,n))
       else if( degreeTime.eq.5 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)))
       else if( degreeTime.eq.6 )then
        time=6.*a(3,n)+t*(24.*a(4,n)+t*(60.*a(5,n)+t*(120.*a(6,n))))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1836
       end if
      else if( dt.eq.4 )then
       ! --- 4th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=24.*a(4,n)
       else if( degreeTime.eq.5 )then
        time=24.*a(4,n)+t*(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=24.*a(4,n)+t*(120.*a(5,n)+t*(360.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.5 )then
       ! --- 5th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=(120.*a(5,n))
       else if( degreeTime.eq.6 )then
        time=(120.*a(5,n)+t*(720.*a(6,n)))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else if( dt.eq.6 )then
       ! --- 6th time derivative ---
       if( degreeTime.eq.0 )then
        time=0.
       else if( degreeTime.eq.1 )then
        time=0.
       else if( degreeTime.eq.2 )then
        time=0.
       else if( degreeTime.eq.3 )then
        time=0.
       else if( degreeTime.eq.4 )then
        time=0.
       else if( degreeTime.eq.5 )then
        time=0.
       else if( degreeTime.eq.6 )then
        time=(720.*a(6,n))
       else
        write(*,*) 'PolyFunction:ERROR invalid degreeTime'
        stop 1837
       end if
      else
        write(*,*) 'PolyFunction:ERROR: too many time derivatives=',dt
        stop 1838
      end if
      do i3=nta,ntb
      do i2=nsa,nsb
      do i1=nra,nrb
      y1=ya(i1,i2,i3)
      y2=y1*y1
      y3=y2*y1
      y4=y3*y1
      x1=xa(i1,i2,i3)
      x1y1=x1*y1
      x1y2=x1y1*y1
      x1y3=x1y2*y1
      x1y4=x1y3*y1
      x2=x1*x1
      x2y1=x1y1*x1
      x2y2=x1y2*x1
      x2y3=x1y3*x1
      x2y4=x1y4*x1
      x3=x2*x1
      x3y1=x2y1*x1
      x3y2=x2y2*x1
      x3y3=x2y3*x1
      x3y4=x2y4*x1
      x4=x3*x1
      x4y1=x3y1*x1
      x4y2=x3y2*x1
      x4y3=x3y3*x1
      x4y4=x3y4*x1
      x5=x4*x1
      x5y1=x4y1*x1
      x5y2=x4y2*x1
      x5y3=x4y3*x1
      x5y4=x4y4*x1
      x6=x5*x1
      x6y1=x5y1*x1
      x6y2=x5y2*x1
      x6y3=x5y3*x1
      x6y4=x5y4*x1
      y5=y4*y1
      y6=y5*y1
      x1y5=x1y4*y1
      x1y6=x1y5*y1
      x2y5=x1y5*x1
      x2y6=x1y6*x1
      x3y5=x2y5*x1
      x3y6=x2y6*x1
      x4y5=x3y5*x1
      x4y6=x3y6*x1
      r(i1,i2,i3,n)=(+c(0,2,0,n)*2.+c(0,3,0,n)*y1*6.+c(0,4,0,n)*y2*12.+
     & c(0,5,0,n)*y3*20.+c(0,6,0,n)*y4*30.+c(1,2,0,n)*x1*2.+c(1,3,0,n)
     & *x1y1*6.+c(1,4,0,n)*x1y2*12.+c(1,5,0,n)*x1y3*20.+c(1,6,0,n)*
     & x1y4*30.+c(2,2,0,n)*x2*2.+c(2,3,0,n)*x2y1*6.+c(2,4,0,n)*x2y2*
     & 12.+c(2,5,0,n)*x2y3*20.+c(2,6,0,n)*x2y4*30.+c(3,2,0,n)*x3*2.+c(
     & 3,3,0,n)*x3y1*6.+c(3,4,0,n)*x3y2*12.+c(3,5,0,n)*x3y3*20.+c(3,6,
     & 0,n)*x3y4*30.+c(4,2,0,n)*x4*2.+c(4,3,0,n)*x4y1*6.+c(4,4,0,n)*
     & x4y2*12.+c(4,5,0,n)*x4y3*20.+c(4,6,0,n)*x4y4*30.+c(5,2,0,n)*x5*
     & 2.+c(5,3,0,n)*x5y1*6.+c(5,4,0,n)*x5y2*12.+c(5,5,0,n)*x5y3*20.+
     & c(5,6,0,n)*x5y4*30.+c(6,2,0,n)*x6*2.+c(6,3,0,n)*x6y1*6.+c(6,4,
     & 0,n)*x6y2*12.+c(6,5,0,n)*x6y3*20.+c(6,6,0,n)*x6y4*30.+c(2,0,0,
     & n)*2.+c(2,1,0,n)*y1*2.+c(2,2,0,n)*y2*2.+c(2,3,0,n)*y3*2.+c(2,4,
     & 0,n)*y4*2.+c(2,5,0,n)*y5*2.+c(2,6,0,n)*y6*2.+c(3,0,0,n)*x1*6.+
     & c(3,1,0,n)*x1y1*6.+c(3,2,0,n)*x1y2*6.+c(3,3,0,n)*x1y3*6.+c(3,4,
     & 0,n)*x1y4*6.+c(3,5,0,n)*x1y5*6.+c(3,6,0,n)*x1y6*6.+c(4,0,0,n)*
     & x2*12.+c(4,1,0,n)*x2y1*12.+c(4,2,0,n)*x2y2*12.+c(4,3,0,n)*x2y3*
     & 12.+c(4,4,0,n)*x2y4*12.+c(4,5,0,n)*x2y5*12.+c(4,6,0,n)*x2y6*
     & 12.+c(5,0,0,n)*x3*20.+c(5,1,0,n)*x3y1*20.+c(5,2,0,n)*x3y2*20.+
     & c(5,3,0,n)*x3y3*20.+c(5,4,0,n)*x3y4*20.+c(5,5,0,n)*x3y5*20.+c(
     & 5,6,0,n)*x3y6*20.+c(6,0,0,n)*x4*30.+c(6,1,0,n)*x4y1*30.+c(6,2,
     & 0,n)*x4y2*30.+c(6,3,0,n)*x4y3*30.+c(6,4,0,n)*x4y4*30.+c(6,5,0,
     & n)*x4y5*30.+c(6,6,0,n)*x4y6*30.)*time
! endLoops()
      end do
      end do
      end do
      end do
      else if( dx.gt.6 .or. dy.gt.6 .or. dz.gt.6 )then
! beginLoops(time=0.)
       do n=nca,ncb
       time=0.
       do i3=nta,ntb
       do i2=nsa,nsb
       do i1=nra,nrb
        r(i1,i2,i3,n)=0.
! endLoops()
       end do
       end do
       end do
       end do
      else
       write(*,*) 'polyFunction:ERROR derivative not implemented'
       write(*,*) 'dx,dy,dz=',dx,dy,dz
       stop 6543
      end if
      return
      end
