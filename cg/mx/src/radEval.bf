! ------------------------------------------------------------------------------------------
! -- evaluate the radiation BC's
! ------------------------------------------------------------------------------------------
#Include "derivMacroDefinitions.h"

#Include "defineParametricDerivMacros.h"

! defineParametricDerivativeMacros(u,dr,dx,DIM,ORDER,COMPONENTS,MAXDERIV)

! 2d, order=4, components=1, derivatives=2:
defineParametricDerivativeMacros(u,dr,dx,2,4,1,2)
! 2d, order=2, components=1, derivatives=4:
defineParametricDerivativeMacros(u,dr,dx,2,2,1,4)

! 3d, order=4, components=1, derivatives=2:
defineParametricDerivativeMacros(u,dr,dx,3,4,1,2)
! 3d, order=2, components=1, derivatives=4:
defineParametricDerivativeMacros(u,dr,dx,3,2,1,4)


defineParametricDerivativeMacros(u1,dr,dx,2,4,1,2)
defineParametricDerivativeMacros(u1,dr,dx,2,2,1,4)

! 2D, order=4, components=2, 1-derivative
defineParametricDerivativeMacros(rsxy,dr,dx,2,4,2,1)
! 2d, odrer=2, components=2, 3-derivatives
defineParametricDerivativeMacros(rsxy,dr,dx,2,2,2,3)



#beginMacro beginLoops(n1a,n1b,n2a,n2b,n3a,n3b,na,nb)
do i3=n3a,n3b
do i2=n2a,n2b
do i1=n1a,n1b
do n=na,nb
  ! write(*,'(" periodic i1,i2,i3,n=",4i4)') i1,i2,i3,n
#endMacro

#beginMacro endLoops()
end do
end do
end do
end do
#endMacro


! ----------------------------------------------------------------------------------------
! Macro: Taylor time-stepping order 2, 3d
!  Note: alpha = coeff of artificial dissipation
! ----------------------------------------------------------------------------------------
#beginMacro taylorbc23d(um1)
 um1 = u(i1,i2,i3,n) + cdt*( ux*is1 - hu3d(i2,i3,n) )  - h*ux \
  +.5*cdt*cdt*utt - cdt*h*( uxx*is1 - hu3dx(i2,i3,n) ) + .5*h*h*uxx\
  +cdt*alpha*( u(im,i2,i3-1,n) + u(im,i2-1,i3,n) -4.*u(im,i2,i3,n) +u(im,i2+1,i3,n) +u(im,i2,i3+1,n) )
#endMacro


! ----------------------------------------------------------------------------------------
! Macro: Taylor time-stepping order 4, 3d
!  Note: alpha = coeff of artificial dissipation
! ----------------------------------------------------------------------------------------
#beginMacro taylorbc43d(um1)
 um1 = u(i1,i2,i3,n) + cdt*ut - h*ux \
  +.5*cdt*cdt*utt -cdt*h*utx + .5*h*h*uxx\
  + cdt*cdt*cdt/6.*uttt -cdt*cdt*h*.5*uttx + cdt*h*h*.5*utxx - h*h*h/6.*uxxx \
  + cdt*cdt*cdt*cdt/24.*utttt -cdt*cdt*cdt*h/6.*utttx + cdt*cdt*h*h/4.*uttxx \
  -cdt*h*h*h/6.*utxxx + h*h*h*h/24.*uxxxx \
  +cdt*(alpha)*( -u(im,i2-2,i3,n)+4.*u(im,i2-1,i3,n)-6.*u(im,i2,i3,n)+4.*u(im,i2+1,i3,n)-u(im,i2+2,i3,n) \
                 -u(im,i2,i3-2,n)+4.*u(im,i2,i3-1,n)-6.*u(im,i2,i3,n)+4.*u(im,i2,i3+1,n)-u(im,i2,i3+2,n)  )
#endMacro

! For 2D
#defineMacro huyy2(i2,n) (hu(i2-1,n)-2.*hu(i2,n)+hu(i2+1,n))/(dx(1)**2)
#defineMacro huxyy2(i2,n) (hux(i2-1,n)-2.*hux(i2,n)+hux(i2+1,n))/(dx(1)**2)

! For 3D
#defineMacro hu3dyy2(i2,i3,n)  ( hu3d(i2-1,i3,n)-2.* hu3d(i2,i3,n)+ hu3d(i2+1,i3,n))/(dx(1)**2)
#defineMacro hu3dxyy2(i2,i3,n) (hu3dx(i2-1,i3,n)-2.*hu3dx(i2,i3,n)+hu3dx(i2+1,i3,n))/(dx(1)**2)

#defineMacro hu3dzz2(i2,i3,n)  ( hu3d(i2,i3-1,n)-2.* hu3d(i2,i3,n)+ hu3d(i2,i3+1,n))/(dx(2)**2)
#defineMacro hu3dxzz2(i2,i3,n) (hu3dx(i2,i3-1,n)-2.*hu3dx(i2,i3,n)+hu3dx(i2,i3+1,n))/(dx(2)**2)



#beginMacro cylTaylorbc2(um1)
 um1 = u(i1,i2,i3,n) + cdt*ut  + hx*ux + hy*uy \
  +.5*cdt*cdt*utt + cdt*hx*utx + cdt*hy*uty + .5*( hx**2*uxx +2.*hx*hy*uxy+hy**2*uyy )\
  +cdt*alpha*( u(i1m-is2,i2m-is1,i3,n)-2.*u(i1m,i2m,i3,n)+u(i1m+is2,i2m+is1,i3,n) )
#endMacro

! May 26, 2020 Fixes for c not equal to 1, change dt to c*dt  
#beginMacro cylTaylorbc4(um1)
 um1 = u(i1,i2,i3,n) + cdt*ut  + hx*ux + hy*uy \
  +.5*cdt*cdt*utt + cdt*hx*utx + cdt*hy*uty + .5*( hx**2*uxx +2.*hx*hy*uxy+hy**2*uyy )\
  +(1./6.)*( cdt**3*uttt + 3.*cdt**2*(hx*uttx+hy*utty) + 3.*cdt*(hx**2*utxx+2.*hx*hy*utxy+hy**2*utyy) \
             +hx**3*uxxx + 3.*hx**2*hy*uxxy + 3.*hx*hy**2*uxyy + hy**3*uyyy ) \
  +(1./24.)*( cdt**4*utttt + 4.*cdt**3*( hx*utttx +hy*uttty) +6.*cdt**2*( hx**2*uttxx+hy**2*uttyy+2.*hx*hy*uttxy )\
              +4.*cdt*( hx**3*utxxx + hy**3*utyyy + 3.*hx**2*hy*utxxy + 3.*hx*hy**2*utxyy )\
              + hx**4*uxxxx + 4.*hx**3*hy*uxxxy + 6.*hx**2*hy**2*uxxyy + 4.*hx*hy**3*uxyyy + hy**4*uyyyy )  \
  +cdt*alpha*( -u(i1m-2*is2,i2m-2*is1,i3,n)+4.*u(i1m-is2,i2m-is1,i3,n)-6.*u(i1m,i2m,i3,n)\
                              +4.*u(i1m+is2,i2m+is1,i3,n)-u(i1m+2*is2,i2m+2*is1,i3,n) )
#endMacro


#beginMacro taylorbc2Old(um1)
 um1 = u(i1,i2,i3,n) + dt*( ux*is1 - hu(i2,n) )  - h*ux \
  +.5*dt*dt*utt - cdt*h*( uxx*is1 - hux(i2,n) ) + .5*h*h*uxx\
  +dt*alpha*( u(im,i2-1,i3,n)-2.*u(im,i2,i3,n)+u(im,i2+1,i3,n) )
#endMacro

! May 26, 2020 Fixes for c not equal to 1, change dt to c*dt  
#beginMacro taylorbc2(um1)
 um1 = u(i1,i2,i3,n) + cdt*( ux*is1 - hu(i2,n) )  - h*ux \
  +.5*cdt*cdt*utt - cdt*h*( uxx*is1 - hux(i2,n) ) + .5*h*h*uxx\
  +cdt*alpha*( u(im,i2-1,i3,n)-2.*u(im,i2,i3,n)+u(im,i2+1,i3,n) )
#endMacro


#beginMacro taylorbc4Old(um1)
 um1 = u(i1,i2,i3,n) + dt*ut - h*ux \
  +.5*dt*dt*utt -dt*h*utx + .5*h*h*uxx\
  + dt*dt*dt/6.*uttt -dt*dt*h*.5*uttx + dt*h*h*.5*utxx - h*h*h/6.*uxxx \
  + dt*dt*dt*dt/24.*utttt -dt*dt*dt*h/6.*utttx + dt*dt*h*h/4.*uttxx \
  -dt*h*h*h/6.*utxxx + h*h*h*h/24.*uxxxx \
  +dt*alpha*( -u(im,i2-2,i3,n)+4.*u(im,i2-1,i3,n)-6.*u(im,i2,i3,n)\
                              +4.*u(im,i2+1,i3,n)-u(im,i2+2,i3,n) )
#endMacro

           ! May 26, 2020 Fixes for c not equal to 1, change dt to c*dt  
#beginMacro taylorbc4(um1)
 um1 = u(i1,i2,i3,n) + cdt*ut - h*ux \
  +.5*cdt*cdt*utt -cdt*h*utx + .5*h*h*uxx\
  + cdt*cdt*cdt/6.*uttt -cdt*cdt*h*.5*uttx + cdt*h*h*.5*utxx - h*h*h/6.*uxxx \
  + cdt*cdt*cdt*cdt/24.*utttt -cdt*cdt*cdt*h/6.*utttx + cdt*cdt*h*h/4.*uttxx \
  -cdt*h*h*h/6.*utxxx + h*h*h*h/24.*uxxxx \
  +cdt*alpha*( -u(im,i2-2,i3,n)+4.*u(im,i2-1,i3,n)-6.*u(im,i2,i3,n)\
                               +4.*u(im,i2+1,i3,n)-u(im,i2+2,i3,n) )
#endMacro


! Not used ?
#beginMacro taylorbc4a(um1)
 um1 = u(i1,i2,i3,n) + dt*ut - h*ux \
  +.5*dt*dt*utt -dt*h*utx + .5*h*h*uxx\
  + dt*dt*dt/6.*uttt -dt*dt*h*.5*uttx + dt*h*h*.5*utxx - h*h*h/6.*uxxx 
#endMacro


! ==================================================================================
!  Macro: update ghost points in 2D using the radiation conditions
! ==================================================================================
#beginMacro updateGhost2d()

  i3=n3a

  is1=0
  is2=0
  is3=0

  if( axis.eq.0 ) then
    is1=1-2*side
    if( side.eq.0 )then 
      nn1b=nn1a
    else
      nn1a=nn1b
    end if
  else
    is2=1-2*side
    if( side.eq.0 )then 
      nn2b=nn2a
    else
      nn2a=nn2b
    end if
  end if

  

  if( kernelType.eq.planar )then
    ! ************* planar interface *******************
    if( axis.ne.0 )then
      write(*,'("radEval:ERROR: not implemented for axis=",i3)') axis
      stop 1163
    end if

   if( orderOfAccuracy.eq.2 )then

    ! write(*,'("radEval: c=",e9.2)') c

    beginLoops(nn1a,nn1b,nn2a,nn2b,nn3a,nn3b,na,nb)
     !     
     ux = ux2(i1,i2,i3,n)
     uxx =uxx2(i1,i2,i3,n)
     uyy =uyy2(i1,i2,i3,n)
     uLap = uxx+uyy
     utt  = uLap  
    

     h=dx(0)*is1
     im=i1-is1
     taylorbc2(um1)

     h=2.*dx(0)*is1
     im=i1-2*is1
     taylorbc2(um2) 

     u1(i1-is1,i2,i3,n)=um1
     u1(i1-2*is1,i2,i3,n)=um2

    endLoops()


   else if( orderOfAccuracy.eq.4 )then


    beginLoops(nn1a,nn1b,nn2a,nn2b,nn3a,nn3b,na,nb)
     !     
     ux = ux4(i1,i2,i3,n)
     uxx =uxx4(i1,i2,i3,n)
     uyy =uyy4(i1,i2,i3,n)
     uLap = uxx+uyy
    
     uxxx = uxxx2(i1,i2,i3,n)
     uxyy = uxyy2(i1,i2,i3,n)
    
     uLapx= uxxx+uxyy

     uxxxx =uxxxx2(i1,i2,i3,n)
     uxxyy =uxxyy2(i1,i2,i3,n)
     uyyyy =uyyyy2(i1,i2,i3,n)

     uLapSq=uxxxx+2.*uxxyy+uyyyy
     uLapxx=uxxxx+uxxyy
    
     huyyi = huyy2(i2,n)
     huxyyi= huxyy2(i2,n)
     huLap = huxx(i2,n)+huyyi

     utx = uxx*is1 - hux(i2,n)
     utxyy = uxxyy*is1 - huxyyi ! utx = uxx - hux  -> (utx)yy = uxxyy - (hux)yy

     ut   = ux*is1 - hu(i2,n)
     utt  = uLap
     uttt = uLapx*is1 -huLap
     uttx = uLapx
     utxx = uxxx*is1 -huxx(i2,n)

     utxxx = uxxxx*is1 - huxxx(i2,n)
     utttt = uLapSq
     utttx = utxxx + utxyy   ! (uxx+uyy).tx
     uttxx = uLapxx

     if( debug.gt.7 )then
      tm=t-dt
      call getExactSolution( xy(i1,i2,i3,0),xy(i1,i2,i3,1),tm,0,0,1,vt)
       write(*,'("** tm=",f6.4," i1,i2=",2i3)') tm,i1,i2
    ! "'
       write(*,'("   ++ ut,vt=",2e10.2," err=",e8.2)') ut,vt,ut-vt

      call getExactSolution( xy(i1,i2,i3,0),xy(i1,i2,i3,1),tm,1,0,0,vx)
      write(*,'("   ++ ux,vx=",2e10.2," err=",e8.2)') ux,vx,ux-vx

      write(*,'("   ++ vt-vx+hu=",e8.2)') vt-vx+hu(i2,n) 


      call getExactSolution( xy(i1,i2,i3,0),xy(i1,i2,i3,1),tm,1,0,1,vtx)
      write(*,'("   ++ utx,vtx=",2e10.2," err=",e8.2)') utx,vtx,utx-vtx

      call getExactSolution( xy(i1,i2,i3,0),xy(i1,i2,i3,1),tm,0,0,3,vttt)
      write(*,'("   ++ uttt,vttt=",2e10.2," err=",e8.2)') uttt,vttt,uttt-vttt

      call getExactSolution( xy(i1,i2,i3,0),xy(i1,i2,i3,1),tm,0,0,4,vtttt)
      write(*,'("   ++ utttt,vtttt=",2e10.2," err=",e8.2)') utttt,vtttt,utttt-vtttt

     end if
    ! ******
    ! uttt=vttt

     h=dx(0)*is1
     im=i1-is1
     taylorbc4(um1)



     h=2.*dx(0)*is1
     im=i1-2*is1
     taylorbc4(um2)

     u1(i1-is1,i2,i3,n)=um1
     u1(i1-2*is1,i2,i3,n)=um2

     if( debug.gt.1 )then
      call getExactSolution( xy(i1-is1,i2,i3,0),xy(i1-is1,i2,i3,1),t,0,0,0,utrue)
      write(*,'(" radEval: t,i2,um1,utrue,err=",e8.2,i3,3e10.3)') t,i2,um1,utrue,um1-utrue
    ! ** u1(i1-1,i2,i3,n)=utrue

      call getExactSolution( xy(i1-2*is1,i2,i3,0),xy(i1-2*is1,i2,i3,1),t,0,0,0,utrue)
      write(*,'("        : t,i2,um2,utrue,err=",e8.2,i3,3e10.3)') t,i2,um2,utrue,um2-utrue
    ! ** u1(i1-2,i2,i3,n)=utrue
     end if

    endLoops()

   else
     stop 8832 ! unknown orderOfAccuracy
   end if

  else if( kernelType.eq.cylindrical )then

     ! **************Cylindrical boundary**********************************
   if( orderOfAccuracy.eq.2 )then
    if( axis.eq.0 )then
      ii=nn2a
    else
      ii=nn1a
    end if
    do i2=nn2a,nn2b
    do i1=nn1a,nn1b
    !     
    !  evalJacobianDerivatives(u,v,DIM,ORDER,MAXDERIV)
    ! NOTE: jacobians need 1 less derivative  ****************** no need to repeat for all components!
    !evalJacobianDerivatives(rsxy,i1,i2,i3,aj4,2,4,1)
     evalJacobianDerivatives(rsxy,i1,i2,i3,aj2,2,2,1)
     do n=na,nb
     
      ! evalParametricDerivativesComponents1(u,i1,i2,i3,m, v,DIM,ORDER,MAXDERIV)
      evalParametricDerivativesComponents1(u,i1,i2,i3,n, uu,2,2,2)

      getDuDx2(uu,aj2,ux)
      getDuDy2(uu,aj2,uy)

      getDuDxx2(uu,aj2,uxx)
      getDuDxy2(uu,aj2,uxy)
      getDuDyy2(uu,aj2,uyy)

      uLap = uxx+uyy
      utt  = uLap
    
      r = sqrt( xy(i1,i2,i3,0)**2 + xy(i1,i2,i3,1)**2 )
      ! outward normal:
      an1= -rsxy(i1,i2,i3,axis,0)*is2
      an2= -rsxy(i1,i2,i3,axis,1)*is2
      aNorm = sqrt( an1**2 + an2**2 )
      an1=an1/aNorm
      an2=an2/aNorm
      uri = an1*ux+an2*uy     ! u.n = u.r 
     
      uxri = an1*uxx+an2*uxy   ! (ux).r = n.grad(ux)
      uyri = an1*uxy+an2*uyy 

      ut = -c*( uri  + u(i1,i2,i3,n)/(2.*r) + hu(ii,n) )
      utx =-c*( uxri +            ux/(2.*r) + hux(ii,n) )
      uty =-c*( uyri +            uy/(2.*r) + huy(ii,n) )

      ! May 26, 2020 Fixes for c not equal to 1, change dt to c*dt  

      if( debug.gt.0 )then
       tm=t-dt
       call ogderiv(ep,1,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),z0,tm,n,ut ) 
       call ogderiv(ep,1,1,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),z0,tm,n,utx ) 
       call ogderiv(ep,1,0,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),z0,tm,n,uty ) 
      end if

      hx = xy(i1-is1,i2-is2,i3,0)- xy(i1,i2,i3,0)
      hy = xy(i1-is1,i2-is2,i3,1)- xy(i1,i2,i3,1)
      i1m=i1-is1
      i2m=i2-is2
      cylTaylorbc2(um1)

      hx = xy(i1-2*is1,i2-2*is2,i3,0)- xy(i1,i2,i3,0)
      hy = xy(i1-2*is1,i2-2*is2,i3,1)- xy(i1,i2,i3,1)
      i1m=i1-2*is1
      i2m=i2-2*is2
      cylTaylorbc2(um2) 

      u1(i1-is1,i2-is2,i3,n)=um1
      u1(i1-2*is1,i2-2*is2,i3,n)=um2

      if( debug.gt.0 )then
       utrue=ogf(ep,xy(i1-is1,i2-is2,i3,0),xy(i1-is1,i2-is2,i3,1),z0,n,t)
       write(*,'(" radEval: t,ii,n,um1,utrue,err=",e8.2,i3,i2,1x,3e10.2)') t,ii,n,um1,utrue,um1-utrue
       u1(i1-is1,i2-is2,i3,n)=utrue

       utrue=ogf(ep,xy(i1-2*is1,i2-2*is2,i3,0),xy(i1-2*is1,i2-2*is2,i3,1),z0,n,t)
       write(*,'("        : t,ii,n,um2,utrue,err=",e8.2,i3,i2,1x,3e10.2)') t,ii,n,um2,utrue,um2-utrue
       u1(i1-2*is1,i2-2*is2,i3,n)=utrue
      end if


      enddo
      ii=ii+1
    enddo
    enddo

   else if( orderOfAccuracy.eq.4 )then
    if( axis.eq.0 )then
      ii=nn2a
    else
      ii=nn1a
    end if
    do i2=nn2a,nn2b
    do i1=nn1a,nn1b
    !     
    !  evalJacobianDerivatives(u,v,DIM,ORDER,MAXDERIV)
    ! NOTE: jacobians need 1 less derivative  ****************** no need to repeat for all components!
     evalJacobianDerivatives(rsxy,i1,i2,i3,aj4,2,4,1)
     evalJacobianDerivatives(rsxy,i1,i2,i3,aj2,2,2,3)
     do n=na,nb
     
      ! *** first eval the first and second derivatives to 4th order
      ! evalParametricDerivativesComponents1(u,i1,i2,i3,m, v,DIM,ORDER,MAXDERIV)
      evalParametricDerivativesComponents1(u,i1,i2,i3,n, uu,2,4,2)

      getDuDx2(uu,aj4,ux)
      getDuDy2(uu,aj4,uy)

      getDuDxx2(uu,aj4,uxx)
      getDuDxy2(uu,aj4,uxy)
      getDuDyy2(uu,aj4,uyy)

      ! Now evaluate the 3rd and 4th derivatives to 2nd order
      evalParametricDerivativesComponents1(u,i1,i2,i3,n, uu,2,2,4)

      getDuDxxx2(uu,aj2,uxxx)
      getDuDxxy2(uu,aj2,uxxy)
      getDuDxyy2(uu,aj2,uxyy)
      getDuDyyy2(uu,aj2,uyyy)

      getDuDxxxx2(uu,aj2,uxxxx)
      getDuDxxxy2(uu,aj2,uxxxy)
      getDuDxxyy2(uu,aj2,uxxyy)
      getDuDxyyy2(uu,aj2,uxyyy)
      getDuDyyyy2(uu,aj2,uyyyy)


      uLap = uxx+uyy
      utt  = uLap
    
      r = sqrt( xy(i1,i2,i3,0)**2 + xy(i1,i2,i3,1)**2 )
      ! outward normal:
      an1= -rsxy(i1,i2,i3,axis,0)*is2
      an2= -rsxy(i1,i2,i3,axis,1)*is2
      aNorm = sqrt( an1**2 + an2**2 )
      an1=an1/aNorm
      an2=an2/aNorm
      uri = an1*ux+an2*uy     ! u.n = u.r 
     
      uxri = an1*uxx+an2*uxy   ! (ux).r = n.grad(ux)
      uyri = an1*uxy+an2*uyy 


      uxxri = an1*uxxx+an2*uxxy   
      uxyri = an1*uxxy+an2*uxyy 
      uyyri = an1*uxyy+an2*uyyy 

      uxxxri = an1*uxxxx+an2*uxxxy   
      uxxyri = an1*uxxxy+an2*uxxyy 
      uxyyri = an1*uxxyy+an2*uxyyy 
      uyyyri = an1*uxyyy+an2*uyyyy 

      ut   =-c*( uri   + u(i1,i2,i3,n)/(2.*r) +   hu(ii,n) )
      utx  =-c*( uxri  +            ux/(2.*r) +  hux(ii,n) )
      uty  =-c*( uyri  +            uy/(2.*r) +  huy(ii,n) )
      utxx =-c*( uxxri +           uxx/(2.*r) + huxx(ii,n) )
      utxy =-c*( uxyri +           uxy/(2.*r) + huxy(ii,n) )
      utyy =-c*( uyyri +           uyy/(2.*r) + huyy(ii,n) )

      utxxx=-c*( uxxxri +         uxxx/(2.*r) + huxxx(ii,n) )
      utxxy=-c*( uxxyri +         uxxy/(2.*r) + huxxy(ii,n) )
      utxyy=-c*( uxyyri +         uxyy/(2.*r) + huxyy(ii,n) )
      utyyy=-c*( uyyyri +         uyyy/(2.*r) + huyyy(ii,n) )

      ! *** note: could also just multiply dt by c ****

      uttt = (utxx+utyy)*c**2

      uttx = (uxxx+uxyy)*c**2
      utty = (uxxy+uyyy)*c**2

      utttt= (uxxxx+2.*uxxyy+uyyyy)*c**4
      uttxx= (uxxxx+uxxyy)*c**2
      uttyy= (uxxyy+uyyyy)*c**2
      uttxy= (uxxxy+uxyyy)*c**2

      utttx= (utxxx+utxyy)*c**2
      uttty= (utxxy+utyyy)*c**2

      ! **** check the 4th-order terms ***
      ! t^4+(4*y+4*x)*t^3+(6*y^2+6*x^2+12*x*y)*t^2+(4*y^3+12*x^2*y+4*x^3+12*x*y^2)*t+x^4+y^4+6*x^2*y^2+4*x*y^3+4*x^3*y


      if( debug.gt.0 )then
       tm=t-dt
       call ogderiv(ep,1,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),z0,tm,n,ut ) 
       call ogderiv(ep,1,1,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),z0,tm,n,utx ) 
       call ogderiv(ep,1,0,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),z0,tm,n,uty ) 

       call ogderiv(ep,3,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),z0,tm,n,uttt )
       call ogderiv(ep,1,2,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),z0,tm,n,utxx) 
       call ogderiv(ep,1,1,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),z0,tm,n,utxy) 
       call ogderiv(ep,1,0,2,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),z0,tm,n,utyy) 

       call ogderiv(ep,3,1,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),z0,tm,n,utttx )
       call ogderiv(ep,3,0,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),z0,tm,n,uttty )

       call ogderiv(ep,1,3,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),z0,tm,n,utxxx) 
       call ogderiv(ep,1,2,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),z0,tm,n,utxxy) 
       call ogderiv(ep,1,1,2,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),z0,tm,n,utxyy) 
       call ogderiv(ep,1,0,3,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),z0,tm,n,utyyy) 


      end if

      hx = xy(i1-is1,i2-is2,i3,0)- xy(i1,i2,i3,0)
      hy = xy(i1-is1,i2-is2,i3,1)- xy(i1,i2,i3,1)
      i1m=i1-is1
      i2m=i2-is2
      cylTaylorbc4(um1)

      hx = xy(i1-2*is1,i2-2*is2,i3,0)- xy(i1,i2,i3,0)
      hy = xy(i1-2*is1,i2-2*is2,i3,1)- xy(i1,i2,i3,1)
      i1m=i1-2*is1
      i2m=i2-2*is2
      cylTaylorbc4(um2) 

      u1(i1-is1,i2-is2,i3,n)=um1
      u1(i1-2*is1,i2-2*is2,i3,n)=um2

      if( debug.gt.0 )then
       utrue=ogf(ep,xy(i1-is1,i2-is2,i3,0),xy(i1-is1,i2-is2,i3,1),z0,n,t)
       write(*,'(" radEval: t,ii,n,um1,utrue,err=",e8.2,i3,i2,1x,3e10.2)') t,ii,n,um1,utrue,um1-utrue
       u1(i1-is1,i2-is2,i3,n)=utrue

       utrue=ogf(ep,xy(i1-2*is1,i2-2*is2,i3,0),xy(i1-2*is1,i2-2*is2,i3,1),z0,n,t)
       write(*,'("        : t,ii,n,um2,utrue,err=",e8.2,i3,i2,1x,3e10.2)') t,ii,n,um2,utrue,um2-utrue
       u1(i1-2*is1,i2-2*is2,i3,n)=utrue
      end if


      enddo
      ii=ii+1
    enddo
    enddo

   else

    write(*,'("Unimplemented order of accuracy")') 
    stop 8263


   end if


  else
    write(*,'("Unknown Kernel type to evaluate")') 
    stop 8264
  end if
#endMacro
! end updateGhost in 2D         


        
! ==================================================================================
!  Macro: update ghost points in 3D using the radiation conditions
! ==================================================================================
#beginMacro updateGhost3d()

  ! alpha=.25
  if( t.le.5*dt .and. debug>1 )then
    write(*,'("radEval:updateGhost for 3D, order=",i2,", alpha=",f5.2)') orderOfAccuracy,alpha
  end if

  is1=0
  is2=0
  is3=0

  if( axis.eq.0 ) then
    is1=1-2*side
    if( side.eq.0 )then 
      nn1b=nn1a
    else
      nn1a=nn1b
    end if
  else if( axis.eq.1 )then
    is2=1-2*side
    if( side.eq.0 )then 
      nn2b=nn2a
    else
      nn2a=nn2b
    end if
  else if( axis.eq.2 )then
    is3=1-2*side
    if( side.eq.0 )then 
      nn3b=nn3a
    else
      nn3a=nn3b
    end if
  else
    stop 6060
  end if

        

  if( kernelType.eq.planar )then
    ! ************* planar interface *******************
    if( axis.ne.0 )then
      write(*,'("radEval:ERROR: not implemented for axis=",i3)') axis
      stop 1163
    end if

    if( orderOfAccuracy.eq.2 )then

      ! write(*,'("radEval: c=",e9.2)') c

      beginLoops(nn1a,nn1b,nn2a,nn2b,nn3a,nn3b,na,nb)
        !     
        ux = ux2(i1,i2,i3,n)
        uxx =uxx2(i1,i2,i3,n)
        uyy =uyy2(i1,i2,i3,n)
        uzz =uzz2(i1,i2,i3,n)
        uLap = uxx+uyy+uzz
        utt  = uLap  

        h=dx(0)*is1
        im=i1-is1
        taylorbc23d(um1)

        h=2.*dx(0)*is1
        im=i1-2*is1
        taylorbc23d(um2) 

        u1(i1-  is1,i2,i3,n)=um1
        u1(i1-2*is1,i2,i3,n)=um2
      endLoops()

    else if( orderOfAccuracy.eq.4 )then

      beginLoops(nn1a,nn1b,nn2a,nn2b,nn3a,nn3b,na,nb)
        !     
        ux = ux4(i1,i2,i3,n)
        uxx =uxx4(i1,i2,i3,n)
        uyy =uyy4(i1,i2,i3,n)
        uzz =uzz4(i1,i2,i3,n)

        ! **TEMP***
        ! uzz=0

        uLap = uxx+uyy+uzz
          
        uxxx = uxxx2(i1,i2,i3,n)
        uxyy = uxyy2(i1,i2,i3,n)
        uxzz = uxzz2(i1,i2,i3,n)
          
        ! **TEMP***
        ! uxzz=0.
        
        uLapx= uxxx+uxyy+uxzz

        uxxxx =uxxxx2(i1,i2,i3,n)
        uxxyy =uxxyy2(i1,i2,i3,n)
        uyyyy =uyyyy2(i1,i2,i3,n)

        uxxzz =uxxzz2(i1,i2,i3,n)
        uyyzz =uyyzz2(i1,i2,i3,n)
        uzzzz =uzzzz2(i1,i2,i3,n)

        ! **TEMP***
        ! uxxzz=0
        ! uyyzz=0
        ! uzzzz=0
        
        uLapSq=uxxxx + 2.*(uxxyy+uxxzz+uyyzz) + uyyyy + uzzzz
        uLapxx=uxxxx+uxxyy+uxxzz
          
        huyyi = hu3dyy2(i2,i3,n)
        huzzi = hu3dzz2(i2,i3,n)

        huxyyi= hu3dxyy2(i2,i3,n)
        huxzzi= hu3dxzz2(i2,i3,n)

        ! **TEMP***
        ! huyyi = 0
        ! huzzi = 0

        ! huxyyi= 0
        ! huxzzi= 0
      

        huLap = hu3dxx(i2,i3,n)+huyyi+huzzi
        utx = uxx*is1 - hu3dx(i2,i3,n)
        utxyy = uxxyy*is1 - huxyyi ! utx = uxx - hux  -> (utx)yy = uxxyy - (hux)yy
        utxzz = uxxzz*is1 - huxzzi ! utx = uxx - hux  -> (utx)zz = uxxzz - (hux)zz

        ut   = ux*is1 - hu3d(i2,i3,n)
        utt  = uLap
        uttt = uLapx*is1 -huLap
        uttx = uLapx
        utxx = uxxx*is1 -hu3dxx(i2,i3,n)

        utxxx = uxxxx*is1 - hu3dxxx(i2,i3,n)
        utttt = uLapSq
        utttx = utxxx + utxyy + utxzz    ! (uxx+uyy+uzz).tx
        uttxx = uLapxx


        ! ghost line 1
        h=dx(0)*is1
        im=i1-is1
        taylorbc43d(um1)

        ! ghost line 2 
  	h=2.*dx(0)*is1
        im=i1-2*is1
        taylorbc43d(um2)

        u1(i1-is1,i2,i3,n)=um1
        u1(i1-2*is1,i2,i3,n)=um2

      endLoops()

    else
      write(*,'("radEval: 3D unknown orderOfAccuracy=",i2)') orderOfACcuracy
      stop 666 
    end if
  else
    write(*,'("radEval:ERROR: kernelType")')
    stop 2222     
  end if   
#endMacro


! =========================================================================================
! Macro: Save the boundary data at time level "m"
! =========================================================================================
#beginMacro saveBoundaryData(AXIS,ORDER)
  beginLoops(nn1a,nn1b,nn2a,nn2b,nn3a,nn3b,na,nb)

    ! here is were we store u.x, u.xx, etc
    kxc   = 1 + numberOfDerivatives*(n)
    #If #ORDER eq "4" 
      kxxc  = 2 + numberOfDerivatives*(n)
      kxxxc = 3 + numberOfDerivatives*(n)
    #End
      
    #If #AXIS eq "0"
      ! -- axis=0 --
      #If #ORDER eq "2" 
       bd(i1,i2,i3,kxc  ,m)=u1x2(i1,i2,i3,n)
      #Else
       bd(i1,i2,i3,kxc  ,m)=u1x4(i1,i2,i3,n)
       bd(i1,i2,i3,kxxc ,m)=u1xx4(i1,i2,i3,n)
       bd(i1,i2,i3,kxxxc,m)=u1xxx2(i1,i2,i3,n)
      #End 
    #Elif #AXIS eq "1"
      ! -- axis=1 --
      #If #ORDER eq "2" 
       bd(i1,i2,i3,kxc  ,m)=u1y2(i1,i2,i3,n)
      #Else
       bd(i1,i2,i3,kxc  ,m)=u1y4(i1,i2,i3,n)
       bd(i1,i2,i3,kxxc ,m)=u1yy4(i1,i2,i3,n)
       bd(i1,i2,i3,kxxxc,m)=u1yyy2(i1,i2,i3,n)
      #End 
    #Else
      ! -- axis=2 --
      #If #ORDER eq "2" 
       bd(i1,i2,i3,kxc  ,m)=u1z2(i1,i2,i3,n)
      #Else
       bd(i1,i2,i3,kxc  ,m)=u1z4(i1,i2,i3,n)
       bd(i1,i2,i3,kxxc ,m)=u1zz4(i1,i2,i3,n)
       bd(i1,i2,i3,kxxxc,m)=u1zzz2(i1,i2,i3,n)
      #End 
    #End

  endLoops()
#endMacro 


      subroutine radEval( nd, nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,\
                          gridIndexRange, u1, u, xy, rsxy, boundaryCondition, \
                          md1a,md1b,md2a,md2b,huv, \
                          ld1a,ld1b,ld2a,ld2b,ld3a,ld3b,huv3d, \
                          sd1a,sd1b,sd2a,sd2b,sd3a,sd3b,sd4a,sd4b,sd5a,sd5b,uSave,\
                          bd1a,bd1b,bd2a,bd2b,bd3a,bd3b,bd4a,bd4b,bd5a,bd5b,bd,\
                          ipar, rpar, ierr )
! ===================================================================================
!  Radition boundary conditions for Maxwell's Equations.
!      
!     Apply the BC of the form  u.t + u.n + H(u) = 0
!
!
!  huv(i,m,n) : Kernel and its derivatives, i=tangential index, m=derivative, n=component
!
! ===================================================================================

      implicit none

      integer nd, \
              nd1a,nd1b,nd2a,nd2b,nd3a,nd3b, \
              n1a,n1b,n2a,n2b,n3a,n3b,na,nb, \
              md1a,md1b,md2a,md2b,md3a,md3b,currentTimeLevel,numberOfTimeLevels,\
              ld1a,ld1b,ld2a,ld2b,ld3a,ld3b,ld4a,ld4b,\
              sd1a,sd1b,sd2a,sd2b,sd3a,sd3b,sd4a,sd4b,sd5a,sd5b,\
              bd1a,bd1b,bd2a,bd2b,bd3a,bd3b,bd4a,bd4b,bd5a,bd5b,\
              ierr

      real u1(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
      real u(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
      real xy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
      real rsxy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1,0:nd-1)

      real uSave(sd1a:sd1b,sd2a:sd2b,sd3a:sd3b,sd4a:sd4b,sd5a:sd5b)

      ! new boundaryData array for parallel vesion
      ! bd(i1,i2,i3,fields,time-levels)
      real bd(bd1a:bd1b,bd2a:bd2b,bd3a:bd3b,bd4a:bd4b,bd5a:bd5b)

      real huv(md1a:md1b,md2a:md2b,0:*)
      real huv3d(ld1a:ld1b,ld2a:ld2b,ld3a:ld3b,0:*)

      integer gridIndexRange,boundaryCondition,useParallelVersion
      integer ipar(0:*)
      real rpar(0:*)
      !     --- local variables ----
      
      real dx(0:2),dr(0:2)
      real t,dt,eps,mu,c,csq,cdt
      integer side,axis,gridType,orderOfAccuracy,grid,kernelType,numberOfComponents,nf,numberOfDerivatives
      integer debug,i1,i2,i3,is1,is2,is3,im,i,ii,j
      integer kxc,kxxc,kxxxc
      !      integer ex,ey,ez,hx,hy,hz
      integer nn1a,nn1b,nn2a,nn2b,nn3a,nn3b, m1a,m1b,m2a,m2b,m3a,m3b,numGhost
      integer mm1,mm2,mm3,mm4
      integer kx,ky,kxx,kxy,kyy,kxxx,kxxy,kxyy,kyyy
      real ux,uxx,uyy,uLap,uxxx,uxxy,uxyy,uyyy,uLapx, uxxxx,uxxxy,uxxyy,uxyyy,uyyyy, uLapSq,uLapxx,uLapyy
      real uzz, uxzz, uzzzz, uxxzz, uyyzz
      real uy,uty,uxy
      real utxyy,ut,utt,utx,uttt,uttx,utxx,utttt, utttx,uttxx,utxxx
      real uxxri,uxyri,uyyri
      real utxy,utyy,utty 
      real uxxxri,uxxyri,uxyyri,uyyyri 
      real utxxy,utyyy, uttxy,uttyy, uttty

      real um1,um2  ,h,alpha
      real utrue
      real v,vt,tm,vtt,vttt,vtttt, vtx, vx
      real hu,hux,huxx,huxxx,huyy,huxy,huLap,huxyy,huxxy,huyyy
      real huyyi,huxyyi
      real huzzi,huxzzi,utxzz

      real hx,hy, huy, z0

      real hu3d,hu3dx,hu3dxx,hu3dxxx
      ! real hu3dy,hu3dxy,hu3dyy,hu3dxyy,hu3dxxy,hu3dyyy

      real ogf
      real ep ! holds the pointer to the TZ function

      integer i1m, i2m, np1, np2
      real r,an1,an2,aNorm,uri,uxri,uyri

      integer rectangular,curvilinear
      parameter(\
        rectangular=0,\
        curvilinear=1)

      integer planar,slab,cylindrical,spherical
      parameter(\
        planar=0,\
        slab=1,\
        cylindrical=2,\
        spherical=3 )

      !  declareTemporaryVariables(DIM,MAXDERIV)
      declareTemporaryVariables(3,4)

      ! declareParametricDerivativeVariables(v,DIM)
      declareParametricDerivativeVariables(uu,3)
      declareJacobianDerivativeVariables(aj2,2)

      declareJacobianDerivativeVariables(aj4,2)
!     declareJacobianDerivativeVariables(aj6,2)

!     --- start statement function ----
      integer kd,m,n
!      declareDifferenceOrder2(u,RX)
!      declareDifferenceOrder4(u,RX)


      hu(i,n)    = huv(i,0,n)
      hux(i,n)   = huv(i,kx,n)
      huxx(i,n)  = huv(i,kxx,n)
      huxxx(i,n) = huv(i,kxxx,n)

      huy(i,n)   = huv(i,ky,n)
      huxy(i,n)  = huv(i,kxy,n)
      huyy(i,n)  = huv(i,kyy,n)
      huxyy(i,n) = huv(i,kxyy,n)
      huxxy(i,n) = huv(i,kxxy,n)
      huyyy(i,n) = huv(i,kyyy,n)

      hu3d(i,j,n)    = huv3d(i,j,0,n)
      hu3dx(i,j,n)   = huv3d(i,j,kx,n)
      hu3dxx(i,j,n)  = huv3d(i,j,kxx,n)
      hu3dxxx(i,j,n) = huv3d(i,j,kxxx,n)

      ! hu3dy(i,j,n)   = huv3d(i,j,ky,n)
      ! hu3dxy(i,j,n)  = huv3d(i,j,kxy,n)
      ! hu3dyy(i,j,n)  = huv3d(i,j,kyy,n)
      ! hu3dxyy(i,j,n) = huv3d(i,j,kxyy,n)
      ! hu3dxxy(i,j,n) = huv3d(i,j,kxxy,n)
      ! hu3dyyy(i,j,n) = huv3d(i,j,kyyy,n)

!      defineDifferenceOrder2Components1(u,RX)
!      defineDifferenceOrder4Components1(u,RX)

!............... end statement functions

      ierr=0

      side                 =ipar(0)
      axis                 =ipar(1)
      grid                 =ipar(2)
      n1a                  =ipar(3)
      n1b                  =ipar(4)
      n2a                  =ipar(5)
      n2b                  =ipar(6)
      n3a                  =ipar(7)
      n3b                  =ipar(8)
      na                   =ipar(9)
      nb                   =ipar(10)
      currentTimeLevel     =ipar(11)
      numberOfTimeLevels   =ipar(12)

      gridType             =ipar(13)
      orderOfAccuracy      =ipar(14)
      debug                =ipar(15)
      kernelType           =ipar(16)
      useParallelVersion   =ipar(17)

      dx(0)                =rpar(0)
      dx(1)                =rpar(1)
      dx(2)                =rpar(2)
      dr(0)                =rpar(3)
      dr(1)                =rpar(4)
      dr(2)                =rpar(5)

      alpha                =rpar(6) ! damping
      ep                   =rpar(7) ! twilight zone pointer

      t                    =rpar(10)
      dt                   =rpar(11)
      eps                  =rpar(12)
      mu                   =rpar(13)
      c                    =rpar(14)
     
      z0=0.

      csq = c**2
      cdt = c*dt
      
      !     numGhost=orderOfAccuracy/2

      ! bounds for loops 
      nn1a=n1a
      nn1b=n1b
      nn2a=n2a
      nn2b=n2b
      nn3a=n3a
      nn3b=n3b

      m = currentTimeLevel
      mm1 = mod(m-1+numberOfTimeLevels,numberOfTimeLevels)
      mm2 = mod(m-2+numberOfTimeLevels,numberOfTimeLevels)
      mm3 = mod(m-3+numberOfTimeLevels,numberOfTimeLevels)
      mm4 = mod(m-4+numberOfTimeLevels,numberOfTimeLevels)


      if( kernelType.eq.planar )then

        ! we may store u.x, u.xx, and u.xxx
        kx=1
        kxx=2
        kxxx=3

      else

        kx=1
        ky=2
        kxx=3
        kxy=4
        kyy=5
        kxxx=6
        kxxy=7
        kxyy=8
        kyyy=9
      end if

      if( nd.eq.2 )then
        ! ------------------ TWO DIMENSIONS -------------------

        updateGhost2d()


      else  
        ! ---------------  THREE DIMENSIONS  --------------------

        updateGhost3d()


      end if


      ! ============ Now save some derivatives ====================

      if( kernelType.eq.cylindrical )then
        ! -- for cylindrical we will need to do a periodic update first ---
        ! update periodic end-pts and ghost-points so we can compute derivatives ux, uy, ... on the boundary
        ! 
        !        |                         |
        !   X  X |                         |  X X
        !   X  X +-------------------------+  X X
        !   X  X n1a                     n1b  X X
        !
        if( axis.eq.1 )then
          ! [n1a,n1b] : indexRange (does not include periodic image) u(n1b+1,i2,i3,n)=u(n1a,i2,i3,n)
          i3=n3a
          np1=n1b-n1a+1
          numGhost=orderOfAccuracy/2
          if( side.eq.0 )then
            m2a=n2a-numGhost
            m2b=n2a+numGhost 
          else
            m2a=n2b-numGhost
            m2b=n2b+numGhost 
          end if
          do n=na,nb
          do i2=m2a,m2b
            do i1=nd1a,n1a-1
              u1(i1,i2,i3,n)=u1(i1+np1,i2,i3,n)
            end do
            do i1=n1b+1,nd1b
              u1(i1,i2,i3,n)=u1(i1-np1,i2,i3,n)
            end do
          end do
          end do
        else
          stop 5532
        end if
      end if


              
      if( useParallelVersion.eq.1 )then

        ! -- save data for new parallel version --- Oct 28, 2020
        if( t.le.2*dt .and. debug.gt.0 )then
          write(*,*) "radEval: useParallelVersion"
       end if
       
        numberOfComponents=nb-na+1

        if( orderOfAccuracy.eq.2 )then
          ! -- 2nd order --
          numberOfDerivatives=2 ! we save u, u.x 

          if( axis.eq.0 )then
            saveBoundaryData(0,2)              
          else if( axis.eq.1 )then
            saveBoundaryData(1,2)              
          else 
            saveBoundaryData(2,2)
          end if 

        else
          ! -- 4th order --
          numberOfDerivatives=4 ! we save u, u.x, u.xx, u.xxx
          if( axis.eq.0 )then
            saveBoundaryData(0,4)              
          else if( axis.eq.1 )then
            saveBoundaryData(1,4)              
          else 
            saveBoundaryData(2,4)
          end if 

        end if    

      else 
        ! --- NOT PARALLEL VERSION ----
      
        i3=n3a
        if( kernelType.eq.planar .and. orderOfAccuracy.eq.2 ) then
  
          ! evaluate one x-derivative to second-order
          if( nd.eq.2 )then
            ! ---- ORDER=2 - TWO DIMENSIONS ---
  
            if( axis.eq.0 )then
              if( side.eq.0 )then
                i1=n1a
              else
                i1=n1b
              end if
              do n=na,nb
                do i2=n2a,n2b
                  uSave(i2,m,1,n,sd5a)=u1x2(i1,i2,i3,n)
                  ! uSave(i2,m,2,n)=u1xx2(i1,i2,i3,n)
                end do
              end do
            else if( axis.eq.1 )then
              if( side.eq.0 )then
                i2=n2a
              else
                i2=n2b
              end if
              do n=na,nb
                do i1=n1a,n1b
                  uSave(i1,m,1,n,sd5a)=u1y2(i1,i2,i3,n)
                  ! uSave(i1,m,2,n)=u1yy2(i1,i2,i3,n)
                end do
              end do
            else
              stop 6639
            end if
  
          else
  
            ! ---- ORDER=2 - THREE DIMENSIONS ---
            if( axis.eq.0 )then
              if( side.eq.0 )then
                i1=n1a
              else
                i1=n1b
              end if
              do n=na,nb
                do i3=n3a,n3b
                do i2=n2a,n2b
                  uSave(i2,i3,m,1,n)=u1x2(i1,i2,i3,n)
                end do
                end do
              end do
            else if( axis.eq.1 )then
              if( side.eq.0 )then
                i2=n2a
              else
                i2=n2b
              end if
              do n=na,nb
                do i3=n3a,n3b
                do i1=n1a,n1b
                  uSave(i1,i3,m,1,n)=u1y2(i1,i2,i3,n)
                end do
                end do
              end do
            else
              stop 6639
            end if
  
          end if
        else if( kernelType.eq.planar .and. orderOfAccuracy.eq.4 ) then
          ! evaluate first 3 derivatives 
  
          if( nd.eq.2 )then
  
            ! ---- ORDER=4 - TWO DIMENSIONS ---
            if( axis.eq.0 )then
              if( side.eq.0 )then
                i1=n1a
              else
                i1=n1b
              end if
              do n=na,nb
                do i2=n2a,n2b
                  uSave(i2,m,1,n,sd5a)=u1x4(i1,i2,i3,n)
                  uSave(i2,m,2,n,sd5a)=u1xx4(i1,i2,i3,n)
                  uSave(i2,m,3,n,sd5a)=u1xxx2(i1,i2,i3,n)
                end do
              end do
            else if( axis.eq.1 )then
              if( side.eq.0 )then
                i2=n2a
              else
                i2=n2b
              end if
              do n=na,nb
                do i1=n1a,n1b
                  uSave(i1,m,1,n,sd5a)=u1y4(i1,i2,i3,n)
                  uSave(i1,m,2,n,sd5a)=u1yy4(i1,i2,i3,n)
                  uSave(i1,m,3,n,sd5a)=u1yyy2(i1,i2,i3,n)
                end do
              end do
            else
              stop 6639
            end if
  
          else
  
            ! ---- ORDER=4 THREE DIMENSIONS ---
            if( axis.eq.0 )then
              if( side.eq.0 )then
                i1=n1a
              else
                i1=n1b
              end if
              do n=na,nb
                do i3=n3a,n3b
                do i2=n2a,n2b
                  uSave(i2,i3,m,1,n)=u1x4(i1,i2,i3,n)
                  uSave(i2,i3,m,2,n)=u1xx4(i1,i2,i3,n)
                  uSave(i2,i3,m,3,n)=u1xxx2(i1,i2,i3,n)
                end do
                end do
              end do
            else if( axis.eq.1 )then
              if( side.eq.0 )then
                i2=n2a
              else
                i2=n2b
              end if
              do n=na,nb
                do i3=n3a,n3b
                do i1=n1a,n1b
                  uSave(i1,i3,m,1,n)=u1y4(i1,i2,i3,n)
                  uSave(i1,i3,m,2,n)=u1yy4(i1,i2,i3,n)
                  uSave(i1,i3,m,3,n)=u1yyy2(i1,i2,i3,n)
                end do
                end do
              end do
            else
              stop 6639
            end if
            
         end if
  
        else if( kernelType.eq.cylindrical .and. orderOfAccuracy.eq.2 ) then
          ! evaluate first 3 derivatives 
  
            if( axis.eq.0 )then
              ii=nn2a
            else
              ii=nn1a
            end if
            do i2=nn2a,nn2b
            do i1=nn1a,nn1b
              evalJacobianDerivatives(rsxy,i1,i2,i3,aj2,2,2,1)
              do n=na,nb
                evalParametricDerivativesComponents1(u1,i1,i2,i3,n, uu,2,2,1)
                getDuDx2(uu,aj2,ux)
                getDuDy2(uu,aj2,uy)
                uSave(ii,m,1,n,sd5a)=ux
                uSave(ii,m,2,n,sd5a)=uy
              end do
              ii=ii+1
            end do
            end do
  
  
        else if( kernelType.eq.cylindrical .and. orderOfAccuracy.eq.4 ) then
            ! evaluate first 3 derivatives 
  
            if( axis.eq.0 )then
              ii=nn2a
            else
              ii=nn1a
            end if
            do i2=nn2a,nn2b
            do i1=nn1a,nn1b
              !  evalJacobianDerivatives(u,v,DIM,ORDER,MAXDERIV)
              evalJacobianDerivatives(rsxy,i1,i2,i3,aj4,2,4,1)
              evalJacobianDerivatives(rsxy,i1,i2,i3,aj2,2,2,3)
              do n=na,nb
                ! evalParametricDerivativesComponents1(u,i1,i2,i3,m, v,DIM,ORDER,MAXDERIV)
                evalParametricDerivativesComponents1(u1,i1,i2,i3,n, uu,2,4,2)
  
                getDuDx2(uu,aj4,ux)
                getDuDy2(uu,aj4,uy)
  
                getDuDxx2(uu,aj4,uxx)
                getDuDxy2(uu,aj4,uxy)
                getDuDyy2(uu,aj4,uyy)
  
                ! Now evaluate the 3rd derivatives to 2nd order
                evalParametricDerivativesComponents1(u1,i1,i2,i3,n, uu,2,2,3)
  
                getDuDxxx2(uu,aj2,uxxx)
                getDuDxxy2(uu,aj2,uxxy)
                getDuDxyy2(uu,aj2,uxyy)
                getDuDyyy2(uu,aj2,uyyy)
  
                evalParametricDerivativesComponents1(u1,i1,i2,i3,n, uu,2,2,1)
  
                uSave(ii,m,1,n,sd5a)=ux
                uSave(ii,m,2,n,sd5a)=uy
                uSave(ii,m,3,n,sd5a)=uxx
                uSave(ii,m,4,n,sd5a)=uxy
                uSave(ii,m,5,n,sd5a)=uyy
                uSave(ii,m,6,n,sd5a)=uxxx
                uSave(ii,m,7,n,sd5a)=uxxy
                uSave(ii,m,8,n,sd5a)=uxyy
                uSave(ii,m,9,n,sd5a)=uyyy
  
  
              end do
              ii=ii+1
            end do
            end do
  
  
        else
         write(*,'("ERROR: unexpected order of accuracy")')
         stop 1132
        end if
      end if ! end not parallel version 

      return
      end
