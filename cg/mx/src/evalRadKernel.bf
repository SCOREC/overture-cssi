      subroutine evalRadKernel( ipar,rpar,amc,alpha,beta,
     &   nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,nd5a,nd5b,
     &   fold2d,fold3d, phi2d,phi3d,
     &   ndzl1a,ndzl1b,ndzl2a,ndzl2b, zl2d, zl3d )
! ========================================================================================
!
! Optimized evaluation of the Radiation Kernel (called by RadiationKernel)
!    See the radbc.pdf notes in CG/DMX/radbc
!
! ========================================================================================

      implicit none
      integer nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,nd5a,nd5b
      integer ndzl1a,ndzl1b,ndzl2a,ndzl2b
      complex*16 zl2d(ndzl1a:ndzl1b)
      complex*16 zl3d(ndzl1a:ndzl1b,ndzl2a:ndzl2b)

      complex*16 fold2d(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
      complex*16 fold3d(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b,nd5a:nd5b)

      complex*16 phi2d(nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
      complex*16 phi3d(nd2a:nd2b,nd3a:nd3b,nd4a:nd4b,nd5a:nd5b)

      ! Coeff of poles 
      complex*16 alpha(0:*), beta(0:*)
      ! Adams-Moulton coeff:
      real amc(-1:*)   

      integer ipar(0:*)
      real rpar(0:*)

      real scl1,scl2,c,dt
      integer nd,numberOfPoles,n1,n2
      integer numberOfModes1,numberOfModes2,orderInTime,mc,currentTimeLevel,numberOfTimeLevels
      integer i1a,i1b,i2a,i2b
      integer k1,k2, k1fft,k2fft

      ! --- local variables  ---
      integer pole,kfft,k,l,n1by2,n2by2,level
      integer lvl(-1:20)
      complex*16 expBeta,xfact,adon,phat
      real w
      logical useCyclicList
      

      ! --- extract parameters 
      nd                 = ipar(0)
      numberOfPoles      = ipar(1)
      orderInTime        = ipar(2)
      n1                 = ipar(3)
      n2                 = ipar(4)
      numberOfModes1     = ipar(5)
      numberOfModes2     = ipar(6)
      i1a                = ipar(7)
      i1b                = ipar(8)
      i2a                = ipar(9)
      i2b                = ipar(10)
      mc                 = ipar(11)
      currentTimeLevel   = ipar(12)
      numberOfTimeLevels = ipar(13)

      scl1            = rpar(0)
      scl2            = rpar(1)
      c               = rpar(2)
      dt              = rpar(3)

      n1by2=n1/2
      n2by2=n2/2

      if( .false. )then
        write(*,'(" evalRadKernel: nd=",i2," mc=",i2," dt=",e10.2)') nd,mc,dt
        write(*,'(" evalRadKernel: numberOfPoles=",i3," i1a,i1b=",2i3)') numberOfPoles,i1a,i1b
        write(*,'(" evalRadKernel: n1,n2=",2i3," numberOfModes=",2i3)') n1,n2,numberOfModes1,numberOfModes2
        write(*,'(" evalRadKernel: ndzl1a:ndzl1b=",2i3)') ndzl1a,ndzl1b
        write(*,'(" evalRadKernel: numberOfTimeLevels,currentTimeLevel=",2i3)') numberOfTimeLevels,currentTimeLevel
      end if

      if( .false. )then
        return
      end if
      
      do l=-1,orderInTime-2
        lvl(l)= mod(currentTimeLevel - l + numberOfTimeLevels, numberOfTimeLevels)
      end do
      
      ! Not clear if useCyclicList is any faster *wdh* Nov 4, 2020
      useCyclicList=.true. ! if true use cyclic list to avoid copies of old data
      
      if( useCyclicList .and. nd.eq.2 )then
        ! ======= TWO DIMENSIONS - CYCLIC LIST ======
        ! -- this version uses cyclic-list "pointers" into fold2d to avoid copying data
        !  -- didn't seem to make mucg difference in timimg  

        ! --- loop over kfft ---
        do pole=0,numberOfPoles-1    
          do kfft=i1a,i1b
            ! k = kfft<n1by2 ? kfft : kfft-n1  // true Fourier mode k
            ! k = kfft - n1*(kfft/n1by2)
            if( kfft .lt. n1by2 )then
               k=kfft
            else
               k=kfft-n1
            end if
             
            ! this next is *wrong*
            ! if( abs(k) > numberOfModes1 )then
            !   continue  ! do not compute this mode
            ! end if   
            if( abs(k) .le.  numberOfModes1 )then
              w = scl1*abs(k)
    
              expBeta = cexp(c*beta(pole)*dt*w)  !  beta's are complex 
              xfact = expBeta
              adon = xfact*phi2d(kfft,pole,mc)     ! add-on 
              do l=0,orderInTime-2
                ! level = mod(currentTimeLevel - l + numberOfTimeLevels, numberOfTimeLevels)
                level=lvl(l)
                adon = adon + xfact*amc(l)*fold2d(level,kfft,pole,mc)
                xfact = xfact * expBeta
              end do 
    
              phat = zl2d(kfft)
              adon = adon + c*dt*amc(-1)*alpha(pole)*w*w*phat
                    
              phi2d(kfft,pole,mc)=adon
  
              ! write(*,*) kfft,phi2d(kfft,pole,mc)
              ! write(*,*) "phat,a,b=",phat,alpha(pole),beta(pole)
                      
              ! do l=orderInTime-2,1,-1
              !   fold2d(l,kfft,pole,mc)=fold2d(l-1,kfft,pole,mc)
              ! end do 
              ! level = mod(currentTimeLevel +1, numberOfTimeLevels) ! next time level 
              level=lvl(-1)
              fold2d(level,kfft,pole,mc)=c*dt*w*w*alpha(pole)*phat 
  
            end if
          end do                   
        end do 
  
        do kfft=i1a,i1b
          zl2d(kfft)=0.
          if( kfft<n1by2 )then
             k=kfft
          else
             k=kfft-n1
          end if
          if( abs(k)<=numberOfModes1 )then
            do pole=0,numberOfPoles-1 
              zl2d(kfft) = zl2d(kfft) + phi2d(kfft,pole,mc)
            end do
          end if
        end do
  
      else if( nd.eq.2 )then
        ! ======= TWO DIMENSIONS -- NO cyclic list  ======
        
        ! --- loop over kfft ---
        do pole=0,numberOfPoles-1    ! "j" loop over number of poles in the approx. of the Kernel
          do kfft=i1a,i1b
            ! k = kfft<n1by2 ? kfft : kfft-n1  // true Fourier mode k
            ! k = kfft - n1*(kfft/n1by2)
            if( kfft .lt. n1by2 )then
               k=kfft
            else
               k=kfft-n1
            end if
             
            ! this next is *wrong*
            ! if( abs(k) > numberOfModes1 )then
            !   continue  ! do not compute this mode
            ! end if   
            if( abs(k) .le.  numberOfModes1 )then
              w = scl1*abs(k)
    
              expBeta = cexp(c*beta(pole)*dt*w)  !  beta's are complex 
              xfact = expBeta
              adon = xfact*phi2d(kfft,pole,mc)     ! add-on 
              do l=0,orderInTime-2
                adon = adon + xfact*amc(l)*fold2d(l,kfft,pole,mc)
                xfact = xfact * expBeta
              end do 
    
              phat = zl2d(kfft)
              adon = adon + c*dt*amc(-1)*alpha(pole)*w*w*phat
                    
              phi2d(kfft,pole,mc)=adon
  
              ! write(*,*) kfft,phi2d(kfft,pole,mc)
              ! write(*,*) "phat,a,b=",phat,alpha(pole),beta(pole)
                      
              do l=orderInTime-2,1,-1
                fold2d(l,kfft,pole,mc)=fold2d(l-1,kfft,pole,mc)
              end do 
              fold2d(0,kfft,pole,mc)=c*dt*w*w*alpha(pole)*phat 
  
            end if
          end do                   
        end do 
  
        do kfft=i1a,i1b
          zl2d(kfft)=0.
          if( kfft<n1by2 )then
             k=kfft
          else
             k=kfft-n1
          end if
          if( abs(k)<=numberOfModes1 )then
            do pole=0,numberOfPoles-1 
              zl2d(kfft) = zl2d(kfft) + phi2d(kfft,pole,mc)
            end do
          end if
        end do
  
      else if( useCyclicList .and. nd.eq.3 ) then
        ! ======= THREE DIMENSIONS USE CYCLIC LIST ======

        ! --- loop over k1fft and k2fft ---
	do pole=0,numberOfPoles-1 ! loop over poles 

	  do  k2fft=i2a,i2b
           ! could use k2 = k2fft -n2*(k2fft/n2by2) 
            if( k2fft .lt. n2by2 )then
               k2=k2fft
            else
               k2=k2fft-n2
            end if
            if( abs(k2) .le. numberOfModes2 )then

              do k1fft=i1a,i1b
                if( k1fft .lt. n1by2 )then
                   k1=k1fft
                else
                   k1=k1fft-n1
                end if
                if( abs(k1) .le. numberOfModes1 )then
  
                  w= sqrt( (scl1*k1)**2 + (scl2*k2)**2 )          ! omega 
    
                  expBeta = exp(c*beta(pole)*dt*w)              ! beta's are complex 
                  xfact = expBeta
                  adon = xfact*phi3d(k1fft,k2fft,pole,mc)          ! add-on
                  do l=0,orderInTime-2
                    level=lvl(l)
                    adon = adon + xfact*amc(l)*fold3d(level,k1fft,k2fft,pole,mc)
                    xfact = xfact * expBeta
                  end do 
    
                  adon = adon + c*dt*amc(-1)*alpha(pole)*w*w*zl3d(k1fft,k2fft)  
                  phi3d(k1fft,k2fft,pole,mc)=adon
    
                  level=lvl(-1)
                  fold3d(level,k1fft,k2fft,pole,mc)=c*dt*w*w*alpha(pole)*zl3d(k1fft,k2fft) 

                end if
                
              end do ! end k1fft

            end if
            
          end do ! end k2fft
        end do ! pole 

        !  ------ Assign zl3d ------
        do k2fft=i2a,i2b
          if( k2fft<n2by2 )then
            k2=k2fft
          else
            k2=k2fft-n2
          end if

          do k1fft=i1a,i1b
            if( k1fft<n1by2 )then
              k1=k1fft
            else
              k1=k1fft-n1
            end if

            zl3d(k1fft,k2fft)=0.

            if( abs(k1)<=numberOfModes1 .and. abs(k2)<=numberOfModes2 )then
              do pole=0,numberOfPoles-1
                zl3d(k1fft,k2fft) = zl3d(k1fft,k2fft) + phi3d(k1fft,k2fft,pole,mc)
              end do
            end if
          end do
        end do


      else if( nd.eq.3 ) then
        ! ======= THREE DIMENSIONS NO CYCLIC LIST ======

        ! --- loop over k1fft and k2fft ---
	do pole=0,numberOfPoles-1 ! loop over poles 

	  do  k2fft=i2a,i2b
           ! could use k2 = k2fft -n2*(k2fft/n2by2) 
            if( k2fft .lt. n2by2 )then
               k2=k2fft
            else
               k2=k2fft-n2
            end if
            if( abs(k2) .le. numberOfModes2 )then

              do k1fft=i1a,i1b
                if( k1fft .lt. n1by2 )then
                   k1=k1fft
                else
                   k1=k1fft-n1
                end if
                if( abs(k1) .le. numberOfModes1 )then
  
                  w= sqrt( (scl1*k1)**2 + (scl2*k2)**2 )          ! omega 
    
                  expBeta = exp(c*beta(pole)*dt*w)              ! beta's are complex 
                  xfact = expBeta
                  adon = xfact*phi3d(k1fft,k2fft,pole,mc)          ! add-on
                  do l=0,orderInTime-2
                    adon = adon + xfact*amc(l)*fold3d(l,k1fft,k2fft,pole,mc)
                    xfact = xfact * expBeta
                  end do 
    
                  adon = adon + c*dt*amc(-1)*alpha(pole)*w*w*zl3d(k1fft,k2fft)  
                  phi3d(k1fft,k2fft,pole,mc)=adon
    
                  do l=orderInTime-2,1,-1
                    fold3d(l,k1fft,k2fft,pole,mc)=fold3d(l-1,k1fft,k2fft,pole,mc)
                  end do 
                  fold3d(0,k1fft,k2fft,pole,mc)=c*dt*w*w*alpha(pole)*zl3d(k1fft,k2fft) 

                end if
                
              end do ! end k1fft

            end if
            
          end do ! end k2fft
        end do ! pole 

        !  ------ Assign zl3d ------
        do k2fft=i2a,i2b
          if( k2fft<n2by2 )then
            k2=k2fft
          else
            k2=k2fft-n2
          end if

          do k1fft=i1a,i1b
            if( k1fft<n1by2 )then
              k1=k1fft
            else
              k1=k1fft-n1
            end if

            zl3d(k1fft,k2fft)=0.

            if( abs(k1)<=numberOfModes1 .and. abs(k2)<=numberOfModes2 )then
              do pole=0,numberOfPoles-1
                zl3d(k1fft,k2fft) = zl3d(k1fft,k2fft) + phi3d(k1fft,k2fft,pole,mc)
              end do
            end if
          end do
        end do

      else ! end 3D
        stop 7777
      end if  

      return
      end 
