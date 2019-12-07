!         -*- mode: F90 -*-
! *********************************************************************
! ********** MACROS FOR DISPERSIVE INTERFACE CONDITIONS ***************
!    This file is included into interface3d.bf 
! *********************************************************************


! -------------------------------------------------------------------------
! *NEW* VERSION July 16, 2019 -- use precomputed coefficients - see paper adegdmi.pdf 
!
! Macro: Evaluate DISPERSIVE forcing terms, 2nd-order accuracy 
!   This macro can be used to eval values in either domain 1 or domain 2
!
! Input:
!   fev(n) : forcing on E equation: E_{tt} = c^2 Delta(E) + ... + fev
!   fpv(n,jv) : forcing on equation for P_{n,jv} 
! Output
!   fp(n) : This is the value of P.tt without the term involving L(E) = c^2*Delta(E)
!   beta = 1 - alphaP*Sum_k{ C_k }
! ------------------------------------------------------------------------
#beginMacro getDispersiveForcingOrder2(k1,k2,k3, fp, fpv,fev, p,pn,pm, u,un,um, dispersionModel,numberOfPolarizationVectors,alphaP,beta,a0v,a1v,b0v,b1v)
  do n=0,nd-1
    fp(n)=0.
  end do
  if( dispersionModel.ne.noDispersion )then

   Csum=0.

   do jv=0,numberOfPolarizationVectors-1

     a0=a0v(jv)
     a1=a1v(jv)
     b0=b0v(jv)
     b1=b1v(jv)
     alpha=alphaP
   
     ! Second-order coefficients: 
     ! Ptt = c2PttLE*LE + c2PttE*E + c2PttEm*Em + c2PttP*P + c2PttPm*Pm + c2PttfE*fE + c2PttfP*fP
     #Include interfaceAdeGdmOrder2.h 

     Csum = Csum + c2PttLE

     ! Bk = 1 + .5*dt*( b1v(jv) + alphaP*a1v(jv) )
     ! Ck = (1./Bk)*a1v(jv)*dt*.5
     ! Csum = Csum + Ck 

     do n=0,nd-1
       pc = n + jv*nd 
       ec = ex +n
       ! P at new time t+dt
       ! Pt, Ptt at time t

       pv   =  p(k1,k2,k3,pc)
       pvn  =  pn(k1,k2,k3,pc)
 
       ev    =  u(k1,k2,k3,ec)
       evn   =  un(k1,k2,k3,ec)

       ! Ptt = c2PttLE*LE + c2PttE*E + c2PttEm*Em + c2PttP*P + c2PttPm*Pm + c2PttfE*fE + c2PttfP*fP
       ! Levae off term: c2PttLE*LE(n)
       fp(n) = fp(n) + c2PttE*ev + c2PttEm*evn + c2PttP*pv + c2PttPm*pvn + c2PttfE*fev(n) + c2PttfP*fpv(n,jv)

       ! write(*,'(" k1,k2,k3=",3i3)') k1,k2,k3
       ! write(*,'(" pc=",i3," p,pn,pm=",3e12.2)') pc, p(k1,k2,k3,pc),pn(k1,k2,k3,pc),pm(k1,k2,k3,pc)
       ! write(*,'(" dt=",e12.2," pv,pvt,pvtt, ev,evt,evtt=",6e12.2)') dt,pv,pvt,pvtt, ev,evt,evtt
       ! write(*,'(" jv=",i2," a0,a1,b0,b1=",4e12.2," Bk,Ck=",2e12.2)') jv,a0v(jv),a1v(jv),b0v(jv),b1v(jv),Bk,Ck
       ! write(*,'(" n=",i2," fev(n)=",e12.2," fp(n)=",e12.2," fpv(n,jv)=",e12.2)') n,fev(n),fp(n),fpv(n,jv)
     end do
   end do
   ! we could precompute D
   beta = 1. -alphaP*Csum
  else
   beta = 1.
  end if
#endMacro



! -------------------------------------------------------------------------
! OLD VERSION
! Macro: Evaluate DISPERSIVE forcing terms, 2nd-order accuracy 
!   This macro can be used to eval values in either domain 1 or domain 2
!
! Input:
!   fev(n) : forcing on E equation: E_{tt} = c^2 Delta(E) + ... + fev
!   fpv(n,jv) : forcing on equation for P_{n,jv} 
! Output
!   fp(n) : 
!   beta = 1 - alphaP*Sum_k{ C_k }
! ------------------------------------------------------------------------
#beginMacro getDispersiveForcingOrder2OLD(k1,k2,k3, fp, fpv,fev, p,pn,pm, u,un,um, dispersionModel,numberOfPolarizationVectors,alphaP,beta,a0v,a1v,b0v,b1v)
  do n=0,nd-1
    fp(n)=0.
  end do
  if( dispersionModel.ne.noDispersion )then
   Csum=0.
   do jv=0,numberOfPolarizationVectors-1
     Bk = 1 + .5*dt*( b1v(jv) + alphaP*a1v(jv) )
     Ck = (1./Bk)*a1v(jv)*dt*.5
     Csum = Csum + Ck 

     do n=0,nd-1
       pc = n + jv*nd 
       ec = ex +n
       ! P at new time t+dt
       ! Pt, Ptt at time t
       pv   =  p(k1,k2,k3,pc)
       pvt  = (p(k1,k2,k3,pc)                   -pm(k1,k2,k3,pc))/(2.*dt)
       pvtt = (p(k1,k2,k3,pc)-2.*pn(k1,k2,k3,pc)+pm(k1,k2,k3,pc))/(dt**2)

       ! E at new time t+dt
       ! Et, Ett at time t
       ev    =  u(k1,k2,k3,ec)
       evt   = (u(k1,k2,k3,ec)                    -um(k1,k2,k3,ec))/(2.*dt) 
       evtt  = (u(k1,k2,k3,ec)-2.*un(k1,k2,k3,ec )+um(k1,k2,k3,ec))/(dt**2)

       fp(n) =fp(n) + (1./Bk)*( \
                    - b1v(jv)*( pvt + .5*dt*pvtt )  \
                    - b0v(jv)*pv + a0v(jv)*ev + a1v(jv)*( evt + .5*dt*(evtt+fev(n)) ) \
                   + fpv(n,jv) \
                              )

       ! write(*,'(" k1,k2,k3=",3i3)') k1,k2,k3
       ! write(*,'(" pc=",i3," p,pn,pm=",3e12.2)') pc, p(k1,k2,k3,pc),pn(k1,k2,k3,pc),pm(k1,k2,k3,pc)
       ! write(*,'(" dt=",e12.2," pv,pvt,pvtt, ev,evt,evtt=",6e12.2)') dt,pv,pvt,pvtt, ev,evt,evtt
       ! write(*,'(" jv=",i2," a0,a1,b0,b1=",4e12.2," Bk,Ck=",2e12.2)') jv,a0v(jv),a1v(jv),b0v(jv),b1v(jv),Bk,Ck
       ! write(*,'(" n=",i2," fev(n)=",e12.2," fp(n)=",e12.2," fpv(n,jv)=",e12.2)') n,fev(n),fp(n),fpv(n,jv)
     end do
   end do
   ! we could precompute D
   beta = 1. -alphaP*Csum
  else
   beta = 1.
  end if
#endMacro


!-------------------------------------------------------------------------------------------
! Macro: Eval twilight-zone forcing for GDM equations
! Output
!  fpv(n,jv) : RHS To Pv_{n,jv} equation 
!  fev(n)    : RHS to E_{n} equation
!-------------------------------------------------------------------------------------------
#beginMacro evalTZforcingGDM(xy,i1,i2,i3,dispersionModel,numberOfPolarizationVectors,c,alphaP,a0v,a1v,b0v,b1v,fpv,fpSum,fev,pettSum)
  if( dispersionModel.ne.noDispersion )then
    do n=0,nd-1
      fpSum(n)=0.
      pettSum(n)=0.

      call ogderiv(ep, 0,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,ex+n, es(n)   ) 
      call ogderiv(ep, 1,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,ex+n, est(n)  )
      call ogderiv(ep, 2,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,ex+n, estt(n) )
      call ogderiv(ep, 0,2,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,ex+n, esxx(n) )
      call ogderiv(ep, 0,0,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,ex+n, esyy(n) )

      do jv=0,numberOfPolarizationVectors-1
        ! The TZ component is offset by pxc
        pc = pxc + jv*nd
        call ogderiv(ep, 0,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,pc+n, pe(n)   )
        call ogderiv(ep, 1,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,pc+n, pet(n)  )
        call ogderiv(ep, 2,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,pc+n, pett(n) )
        ! Normal TZ forcing for P_{n,jv} equation: 
        fpv(n,jv) = pett(n) + b1v(jv)*pet(n) + b0v(jv)*pe(n) - a0v(jv)*es(n) - a1v(jv)*est(n)
        ! Keep sum: 
        fpSum(n)  = fpSum(n)  + fpv(n,jv)
        pettSum(n) = pettSum(n) + pett(n) 
      end do 

      ! TZ forcing for E_{n} equation:
      ! E_tt - c^2 Delta E + alphaP*Ptt  = 
      fev(n) = estt(n) - c**2*( esxx(n) + esyy(n) ) + alphaP*pettSum(n)
    end do
  end if
#endMacro

!-------------------------------------------------------------------------------------------
! Macro: Evaluate TZ forcing for dispersive equations in 2D 
!
! Output
!    fpv1(n,jv) : RHS To Pv_{n,jv} equation on domain 1
!    fpv2(n,jv) : RHS To Pv_{n,jv} equation on domain 2
!    fev1(n)    : RHS to E_{n} equation on domain 1
!    fev2(n)    : RHS to E_{n} equation on domain 2
!-------------------------------------------------------------------------------------------
#beginMacro getDispersiveTZForcing(fpv1,fpv2,fev1,fev2)

  if( twilightZone.eq.1 )then
    evalTZforcingGDM(xy1,i1,i2,i3,dispersionModel1,numberOfPolarizationVectors1,c1,alphaP1,a0v1,a1v1,b0v1,b1v1,fpv1,fpSum1,fev1,pettSum1)
    evalTZforcingGDM(xy2,j1,j2,j3,dispersionModel2,numberOfPolarizationVectors2,c2,alphaP2,a0v2,a1v2,b0v2,b1v2,fpv2,fpSum2,fev2,pettSum2)
  end if

#endMacro 

!-------------------------------------------------------------------------------------------
! Macro: Eval twilight-zone forcing for GDM equations THREE-D
! Output
!  fpv(n,jv) : RHS To Pv_{n,jv} equation 
!  fev(n)    : RHS to E_{n} equation
!-------------------------------------------------------------------------------------------
#beginMacro evalTZforcingGDM3d(xy,i1,i2,i3,dispersionModel,numberOfPolarizationVectors,c,alphaP,a0v,a1v,b0v,b1v,fpv,fpSum,fev,pettSum)
  if( dispersionModel.ne.noDispersion )then
    do n=0,nd-1
      fpSum(n)=0.
      pettSum(n)=0.

      call ogderiv(ep, 0,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, es(n)   ) 
      call ogderiv(ep, 1,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, est(n)  )
      call ogderiv(ep, 2,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, estt(n) )
      call ogderiv(ep, 0,2,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esxx(n) )
      call ogderiv(ep, 0,0,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esyy(n) )
      call ogderiv(ep, 0,0,0,2, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, eszz(n) )

      do jv=0,numberOfPolarizationVectors-1
        ! The TZ component is offset by pxc
        pc = pxc + jv*nd
        call ogderiv(ep, 0,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pe(n)   )
        call ogderiv(ep, 1,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pet(n)  )
        call ogderiv(ep, 2,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pett(n) )
        ! Normal TZ forcing for P_{n,jv} equation: 
        fpv(n,jv) = pett(n) + b1v(jv)*pet(n) + b0v(jv)*pe(n) - a0v(jv)*es(n) - a1v(jv)*est(n)
        ! Keep sum: 
        fpSum(n)  = fpSum(n)  + fpv(n,jv)
        pettSum(n) = pettSum(n) + pett(n) 
      end do 

      ! TZ forcing for E_{n} equation:
      ! E_tt - c^2 Delta E + alphaP*Ptt  = 
      fev(n) = estt(n) - c**2*( esxx(n) + esyy(n) + eszz(n) ) + alphaP*pettSum(n)
    end do
  end if
#endMacro

!-------------------------------------------------------------------------------------------
! Macro: Evaluate TZ forcing for dispersive equations in 3D 
!
! Output
!    fpv1(n,jv) : RHS To Pv_{n,jv} equation on domain 1
!    fpv2(n,jv) : RHS To Pv_{n,jv} equation on domain 2
!    fev1(n)    : RHS to E_{n} equation on domain 1
!    fev2(n)    : RHS to E_{n} equation on domain 2
!-------------------------------------------------------------------------------------------
#beginMacro getDispersiveTZForcing3d(fpv1,fpv2,fev1,fev2)

  if( twilightZone.eq.1 )then
    evalTZforcingGDM3d(xy1,i1,i2,i3,dispersionModel1,numberOfPolarizationVectors1,c1,alphaP1,a0v1,a1v1,b0v1,b1v1,fpv1,fpSum1,fev1,pettSum1)
    evalTZforcingGDM3d(xy2,j1,j2,j3,dispersionModel2,numberOfPolarizationVectors2,c2,alphaP2,a0v2,a1v2,b0v2,b1v2,fpv2,fpSum2,fev2,pettSum2)
  end if

#endMacro 

! ---------------------------------------------------------------------------------------
! Macro: Assign DISPERSIVE interface ghost values, DIM=2, ORDER=2, GRID=Rectangular
! 
! Here are the jump conditions (See notes in DMX_ADE)
!   [ u.x + v.y ] = 0
!   [ (1/mu)* tv,.( curl(E) ) ]
!   [ tv.( c^2*Delta(E) -alphaP*P_tt) ] = 0  --> [ tv.( beta*c^2*Delta(E) - alphaP* F) ]=0 
!   [ (1/mu)* nv.( Delta(E) ) ]=0
! 
! -------------------------------------------------------------------------------------------
#beginMacro assignDispersiveInterfaceGhost22r()

 ! ****************************************************
 ! ***********  2D, ORDER=2, RECTANGULAR **************
 ! ****************************************************

INFO("22r-GDM")

! For rectangular, both sides must axis axis1==axis2: 
if( axis1.ne.axis2 )then
  stop 8826
end if

! 
! Solve for the unknowns u1(-1),v1(-1),u2(-1),v2(-1)
!     
!       A [ U ] = A [ U(old) ] - [ f ]
!
!               [ u1(-1) ]
!       [ U ] = [ v1(-1) ]
!               [ u2(-1) ]
!               [ v2(-1) ]
!             

! --- initialize some forcing functions ---
do n=0,nd-1
  fev1(n)=0.
  fev2(n)=0.
  do jv=0,numberOfPolarizationVectors1-1
    fpv1(n,jv)=0.
  end do
  do jv=0,numberOfPolarizationVectors2-1
    fpv2(n,jv)=0.
  end do
end do

! ----------------- START LOOP OVER INTERFACE -------------------------
beginLoopsMask2d()

  ! first evaluate the equations we want to solve with the wrong values at the ghost points:

  evalInterfaceDerivatives2d()
 
  ! Evaluate TZ forcing for dispersive equations in 2D 
  getDispersiveTZForcing(fpv1,fpv2,fev1,fev2)

  ! eval dispersive forcings for domain 1
  getDispersiveForcingOrder2(i1,i2,i3, fp1, fpv1,fev1,p1,p1n,p1m, u1,u1n,u1m, dispersionModel1,numberOfPolarizationVectors1,alphaP1,beta1,a0v1,a1v1,b0v1,b1v1)

  ! eval dispersive forcings for domain 2
  getDispersiveForcingOrder2(j1,j2,j3, fp2, fpv2,fev2,p2,p2n,p2m, u2,u2n,u2m, dispersionModel2,numberOfPolarizationVectors2,alphaP2,beta2,a0v2,a1v2,b0v2,b1v2)


  if( axis1.eq.0 )then
    ! Interface equations for a boundary at x = 0 or x=1

    ! ---- EQUATION 0 -----
    ! NOTE: if mu==mu2 then we do not need TZ forcing for this eqn:
    f(0)=(u1x+v1y) - \
         (u2x+v2y)
    a4(0,0) = -is1/(2.*dx1(axis1))    ! coeff of u1(-1) from [u.x+v.y] 
    a4(0,1) = 0.                      ! coeff of v1(-1) from [u.x+v.y] 
    a4(0,2) =  js1/(2.*dx2(axis2))    ! coeff of u2(-1) from [u.x+v.y] 
    a4(0,3) = 0.                      ! coeff of v2(-1) from [u.x+v.y]
  
    ! ---- EQUATION 1 -----
    ! NOTE: if mu==mu2 then we do not need TZ forcing for this eqn:
    f(1)=(v1x-u1y)/mu1 - \
         (v2x-u2y)/mu2
    a4(1,0) = 0.
    a4(1,1) = -is1/(2.*dx1(axis1))    ! coeff of v1(-1) from [v.x - u.y] 
    a4(1,2) = 0.
    a4(1,3) =  js1/(2.*dx2(axis2))    ! coeff of v2(-1) from [v.x - u.y]
   
    ! ---- EQUATION 2 -----    
    ! NOTE: if mu==mu2 then we do not need TZ forcing for this eqn:
    f(2)=( (u1xx+u1yy)/mu1 ) - \
         ( (u2xx+u2yy)/mu2 )
    a4(2,0) = 1./(dx1(axis1)**2)/mu1   ! coeff of u1(-1) from [(u.xx + u.yy)/mu]
    a4(2,1) = 0. 
    a4(2,2) =-1./(dx2(axis2)**2)/mu2   ! coeff of u2(-1) from [(u.xx + u.yy)/mu]
    a4(2,3) = 0. 
  
    ! ---- EQUATION 3 -----    
    ! The coefficient of Delta(E) in this equation is altered due to Ptt term 
    f(3)=( (v1xx+v1yy)*beta1/epsmu1 -alphaP1*fp1(1) ) - \
         ( (v2xx+v2yy)*beta2/epsmu2 -alphaP2*fp2(1))

    ! TEST 
    if( .false. )then
      f(3)=( (v1xx+v1yy)*beta1/epsmu1 ) - \
           ( (v2xx+v2yy)*beta2/epsmu2 )
    end if
    a4(3,0) = 0.                      
    a4(3,1) = (beta1/epsmu1)/(dx1(axis1)**2) ! coeff of v1(-1) from [beta*c^2*(v.xx+v.yy)]
    a4(3,2) = 0. 
    a4(3,3) =-(beta2/epsmu2)/(dx2(axis2)**2) ! coeff of v2(-1) from [beta*c^2*(v.xx+v.yy)]
  else

    ! Interface equations for a boundary at y = 0 or y=1
    ! Switch u <-> v,  x<-> y in above equations 

    ! ---- EQUATION 0 -----
    f(0)=(v1y+u1x) - \
         (v2y+u2x)
    a4(0,0) = 0.                      ! coeff of u1(-1) from [u.x+v.y] 
    a4(0,1) = -is1/(2.*dx1(axis1))    ! coeff of v1(-1) from [u.x+v.y] 

    a4(0,2) = 0.                      ! coeff of u2(-1) from [u.x+v.y] 
    a4(0,3) = js1/(2.*dx2(axis2))     ! coeff of v2(-1) from [u.x+v.y]
  
    ! ---- EQUATION 1 -----
    f(1)=(u1y-v1x)/mu1 - \
         (u2y-v2x)/mu2
    a4(1,0) = -is1/(2.*dx1(axis1))
    a4(1,1) = 0.
    a4(1,2) =  js1/(2.*dx2(axis2))  
    a4(1,3) = 0.
   
    ! ---- EQUATION 2 -----    
    f(2)=( (v1xx+v1yy)/mu1 ) - \
         ( (v2xx+v2yy)/mu2 )
    a4(2,0) = 0.
    a4(2,1) = 1./(dx1(axis1)**2)/mu1  
    a4(2,2) = 0.
    a4(2,3) =-1./(dx2(axis2)**2)/mu2 
  
    ! ---- EQUATION 3 -----    
    ! The coefficient of Delta(E) in this equation is altered due to Ptt term 
    f(3)=( (u1xx+u1yy)*beta1/epsmu1 -alphaP1*fp1(0) ) - \
         ( (u2xx+u2yy)*beta2/epsmu2 -alphaP2*fp2(0))
    a4(3,0) = (beta1/epsmu1)/(dx1(axis1)**2)
    a4(3,1) = 0.
    a4(3,2) =-(beta2/epsmu2)/(dx2(axis2)**2) 
    a4(3,3) = 0.


  end if


   q(0) = u1(i1-is1,i2-is2,i3,ex)
   q(1) = u1(i1-is1,i2-is2,i3,ey)
   q(2) = u2(j1-js1,j2-js2,j3,ex)
   q(3) = u2(j1-js1,j2-js2,j3,ey)

   if( .false. .or. debug.gt.4 )then 
     write(*,'("BEFORE: --> i1,i2=",2i4," j1,j2=",2i4," f()=",4e10.2)') i1,i2,j1,j2,f(0),f(1),f(2),f(3)
     write(*,'("     beta1,beta2=",2e10.2," fp1=",2e10.2," fp2=",2e10.2)') beta1,beta2,fp1(0),fp1(1),fp2(0),fp2(1)
     write(*,'("     mu1,mu2=",2e10.2," v1y,u1x,v2y,u2x=",4e10.2)') mu1,mu2,v1y,u1x,v2y,u2x
   end if


   ! subtract off the contributions from the wrong values at the ghost points:
   do n=0,3
     f(n) = (a4(n,0)*q(0)+a4(n,1)*q(1)+a4(n,2)*q(2)+a4(n,3)*q(3)) - f(n)
   end do
   ! write(debugFile,'(" --> i1,i2=",2i4," f(subtract)=",4f8.3)') i1,i2,f(0),f(1),f(2),f(3)
   ! solve A Q = F
   ! factor the matrix
   numberOfEquations=4
   call dgeco( a4(0,0), numberOfEquations, numberOfEquations, ipvt(0),rcond,work(0))
   ! solve
   ! write(debugFile,'(" --> i1,i2=",2i4," rcond=",e10.2)') i1,i2,rcond
   job=0
   call dgesl( a4(0,0), numberOfEquations, numberOfEquations, ipvt(0), f(0), job)
   ! write(debugFile,'(" --> i1,i2=",2i4," f(solve)=",4f8.3)') i1,i2,f(0),f(1),f(2),f(3)

   u1(i1-is1,i2-is2,i3,ex)=f(0)
   u1(i1-is1,i2-is2,i3,ey)=f(1)
   u2(j1-js1,j2-js2,j3,ex)=f(2)
   u2(j1-js1,j2-js2,j3,ey)=f(3)

   if( .false. .or. debug.gt.4 )then 
     ! CHECK: re-evaluate the jump conditions
     evalInterfaceDerivatives2d()

     if( axis1.eq.0 )then
        f(0)=(u1x+v1y) - \
             (u2x+v2y)
        f(1)=(v1x-u1y)/mu1 - \
             (v2x-u2y)/mu2
        f(2)=( (u1xx+u1yy)/mu1 ) - \
             ( (u2xx+u2yy)/mu2 )
        f(3)=( (v1xx+v1yy)*beta1/epsmu1 -alphaP1*fp1(1) ) - \
             ( (v2xx+v2yy)*beta2/epsmu2 -alphaP2*fp2(1))
       ! TEST 
        if( .false. )then
          f(3)=( (v1xx+v1yy)*beta1/epsmu1 ) - \
               ( (v2xx+v2yy)*beta2/epsmu2 )
        end if
      else
        f(0)=(v1y+u1x) - \
             (v2y+u2x)
        f(1)=(u1y-v1x)/mu1 - \
             (u2y-v2x)/mu2
        f(2)=( (v1xx+v1yy)/mu1 ) - \
             ( (v2xx+v2yy)/mu2 )    
        f(3)=( (u1xx+u1yy)*beta1/epsmu1 -alphaP1*fp1(0) ) - \
             ( (u2xx+u2yy)*beta2/epsmu2 -alphaP2*fp2(0))
      end if 
      write(*,'("AFTER: --> i1,i2=",2i4," j1,j2=",2i4," f(re-eval)=",4e10.2)') i1,i2,j1,j2,f(0),f(1),f(2),f(3)


      if( twilightZone.eq.1 )then
        ! check errors in the ghost 
          k1=i1-is1
          k2=i2-is2
          k3=i3
          do n=0,nd-1
            call ogderiv(ep, 0,0,0,0, xy1(k1,k2,k3,0),xy1(k1,k2,k3,1),0.,t,ex+n, es(n)   ) 
          end do
          f(0) =  u1(i1-is1,i2-is2,i3,ex) -es(0)
          f(1) =  u1(i1-is1,i2-is2,i3,ey) -es(1)
          k1=j1-js1
          k2=j2-js2
          k3=j3
          do n=0,nd-1
            call ogderiv(ep, 0,0,0,0, xy2(k1,k2,k3,0),xy2(k1,k2,k3,1),0.,t,ex+n, est(n)   ) 
          end do
          f(2) =  u2(j1-js1,j2-js2,j3,ex) -est(0)
          f(3) =  u2(j1-js1,j2-js2,j3,ey) -est(1)
          write(*,'(" ghost err =",4e10.2)') f(0),f(1),f(2),f(3) 
  
      end if
   end if


   ! -------------------------------------------------------
   ! No need to solve for Hz as it is just an ODE
   ! -------------------------------------------------------
 endLoopsMask2d()
 
 ! stop 7777

#endMacro


! --------------------------------------------------------------------------------------------
! Macro:  Evaluate the RHS to the jump conditons: 2D, Order=2, Dispersive
!
! --------------------------------------------------------------------------------------------
#beginMacro eval2dJumpDispersiveOrder2()
 f(0)=(u1x+v1y) - \
      (u2x+v2y)
 f(1)=( an1*u1Lap +an2*v1Lap )/mu1 - \
      ( an1*u2Lap +an2*v2Lap )/mu2 
 f(2)=(v1x-u1y)/mu1 - \
      (v2x-u2y)/mu2
 f(3)=( ( tau1*u1Lap +tau2*v1Lap )*beta1/epsmu1 - alphaP1*(tau1*fp1(0)+tau2*fp1(1)) ) - \
      ( ( tau1*u2Lap +tau2*v2Lap )*beta2/epsmu2 - alphaP2*(tau1*fp2(0)+tau2*fp2(1)) )

 if( twilightZone.eq.1 )then

   ! For now we assume mu1=mu2 and TZ solutions are the same on both sides.

   ! f(3) = [ tv.E.tt] = [ tv.( c^2*Delta(E) - alphaP*P.tt + fev) ] 
   !      = [ tv.( c^2*Delta(E) - alphaP*P.tt] + [ tv.fev ]
  
   ! -- add on the jump in the forcing ---
   f(3) = f(3) + ( tau1*(fev1(0)-fev2(0)) + tau2*(fev1(1)-fev2(1)) )

  !-    ! f(3) = [ tv.( c^2*Delta(E) - alphaP*P.tt ] - [ tv.( c^2*Delta(E^e) - alphaP*(P^e).tt ] =0 
  !-    call ogderiv(ep, 0,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, uexx )
  !-    call ogderiv(ep, 0,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, ueyy )
  !-    call ogderiv(ep, 0,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, vexx )
  !-    call ogderiv(ep, 0,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, veyy )
  !- 
  !-    ueLap = uexx + ueyy
  !-    veLap = vexx + veyy
  !-    f(3) = f(3) - ( tau1*ueLap +tau2*veLap )*(1./epsmu1-1./epsmu2)
  !- 
  !-    f(3) = f(3) +  alphaP1*( tau1*pettSum1(0) + tau2*pettSum1(1) ) - \
  !-                   alphaP2*( tau1*pettSum2(0) + tau2*pettSum2(1) )
  !-    end if 


   ! write(debugFile,'(" u1Lap,ueLap=",2e10.2," v1Lap,veLap=",2e10.2)') u1Lap,ueLap,v1Lap,veLap

   ! write(debugFile,'(" u1Lap,ueLap=",2e10.2," v1Lap,veLap=",2e10.2)') u1Lap,ueLap,v1Lap,veLap

 end if

#endMacro

! ----------------------------------------------------------------------------------
!  Macro:
!    Evaluate the interface equations for checking the coefficients
! ----------------------------------------------------------------------------------
#beginMacro evalDispersiveInterfaceEquations22c()
  evalInterfaceDerivatives2d()
  eval2dJumpDispersiveOrder2()
#endMacro


! --------------------------------------------------------------------
! Macro: Assign  DISPERSIVE interface ghost values, DIM=2, ORDER=2, GRID=Curvilinear
! 
! Here are the jump conditions (See notes in DMX_ADE)
!   [ u.x + v.y ] = 0
!   [ (1/mu)* tv,.( curl(E) ) ]
!   [ tv.( c^2*Delta(E) -alphaP*P_tt) ] = 0  --> [ tv.( beta*c^2*Delta(E) - alphaP* F) ]=0 
!   [ (1/mu)* nv.( Delta(E) ) ]=0
! 
! -------------------------------------------------------------------------------------------
#beginMacro assignDispersiveInterfaceGhost22c()

  ! ****************************************************
  ! ***********  2D, ORDER=2, CURVILINEAR **************
  ! ****************************************************

INFO("22c-GDM")

! --- initialize some forcing functions ---
do n=0,nd-1
  fev1(n)=0.
  fev2(n)=0.
  do jv=0,numberOfPolarizationVectors1-1
    fpv1(n,jv)=0.
  end do
  do jv=0,numberOfPolarizationVectors2-1
    fpv2(n,jv)=0.
  end do
end do

! ----------------- START LOOP OVER INTERFACE -------------------------
beginLoopsMask2d()

  ! here is the normal (assumed to be the same on both sides)
  an1=rsxy1(i1,i2,i3,axis1,0)   ! normal (an1,an2)
  an2=rsxy1(i1,i2,i3,axis1,1)
  aNorm=max(epsx,sqrt(an1**2+an2**2))
  an1=an1/aNorm
  an2=an2/aNorm
  tau1=-an2
  tau2= an1

  ! first evaluate the equations we want to solve with the wrong values at the ghost points:
  evalInterfaceDerivatives2d()
  ! if( .true. .or. debug.gt.4 )then 
  !    write(*,'(" START  (i1,i2)=",2i3," v=Ey(-1:1,-1:1)",9e14.6)') i1,i2, ((u1(i1+k1,i2+k2,i3,ey),k1=-1,1),k2=-1,1)
  !    write(*,'(" START    mu1,mu2=",2e10.2," v1y,u1x,v2y,u2x=",4e10.2)') mu1,mu2,v1y,u1x,v2y,u2x
  !  end if

  ! Evaluate TZ forcing for dispersive equations in 2D 
  getDispersiveTZForcing(fpv1,fpv2,fev1,fev2)

  ! eval dispersive forcings for domain 1
  getDispersiveForcingOrder2(i1,i2,i3, fp1, fpv1,fev1,p1,p1n,p1m, u1,u1n,u1m, dispersionModel1,numberOfPolarizationVectors1,alphaP1,beta1,a0v1,a1v1,b0v1,b1v1)

  ! eval dispersive forcings for domain 2
  getDispersiveForcingOrder2(j1,j2,j3, fp2, fpv2,fev2,p2,p2n,p2m, u2,u2n,u2m, dispersionModel2,numberOfPolarizationVectors2,alphaP2,beta2,a0v2,a1v2,b0v2,b1v2)

  ! Evaulate RHS, f(n),n=0,1,2,3 using current ghost values: 
  eval2dJumpDispersiveOrder2()

  ! write(debugFile,'(" --> order2-curv: i1,i2=",2i4," f(start)=",4f8.3)') i1,i2,f(0),f(1),f(2),f(3)
  ! write(debugFile,'(" --> u1(ghost),u1=",4f8.3)') u1(i1-is1,i2-is2,i3,ex),u1(i1,i2,i3,ex)
  ! write(debugFile,'(" --> u2(ghost),u2=",4f8.3)') u2(j1-js1,j2-js2,j3,ex),u2(j1,j2,j3,ex)
  ! '

  ! here is the matrix of coefficients for the unknowns u1(-1),v1(-1),u2(-1),v2(-1)
  ! Solve:
  !     
  !       A [ U ] = A [ U(old) ] - [ f ]
  ! ---- EQUATION 0 ----- 
  a4(0,0) = -is*rsxy1(i1,i2,i3,axis1,0)/(2.*dr1(axis1))    ! coeff of u1(-1) from [u.x+v.y] 
  a4(0,1) = -is*rsxy1(i1,i2,i3,axis1,1)/(2.*dr1(axis1))    ! coeff of v1(-1) from [u.x+v.y] 
  a4(0,2) =  js*rsxy2(j1,j2,j3,axis2,0)/(2.*dr2(axis2))    ! coeff of u2(-1) from [u.x+v.y] 
  a4(0,3) =  js*rsxy2(j1,j2,j3,axis2,1)/(2.*dr2(axis2))    ! coeff of v2(-1) from [u.x+v.y] 

  ! ---- EQUATION 2 ----- 
  a4(2,0) =  is*rsxy1(i1,i2,i3,axis1,1)/(2.*dr1(axis1))/mu1   ! coeff of u1(-1) from [(v.x - u.y)/mu] 
  a4(2,1) = -is*rsxy1(i1,i2,i3,axis1,0)/(2.*dr1(axis1))/mu1   ! coeff of v1(-1) from [(v.x - u.y)/mu] 

  a4(2,2) = -js*rsxy2(j1,j2,j3,axis2,1)/(2.*dr2(axis2))/mu2   ! coeff of u2(-1) from [(v.x - u.y)/mu] 
  a4(2,3) =  js*rsxy2(j1,j2,j3,axis2,0)/(2.*dr2(axis2))/mu2   ! coeff of v2(-1) from [(v.x - u.y)/mu] 


  ! coeff of u(-1) from lap = u.xx + u.yy
  rxx1(0,0,0)=aj1rxx
  rxx1(1,0,0)=aj1sxx
  rxx1(0,1,1)=aj1ryy
  rxx1(1,1,1)=aj1syy

  rxx2(0,0,0)=aj2rxx
  rxx2(1,0,0)=aj2sxx
  rxx2(0,1,1)=aj2ryy
  rxx2(1,1,1)=aj2syy

  ! clap1=(rsxy1(i1,i2,i3,axis1,0)**2+rsxy1(i1,i2,i3,axis1,1)**2)/(dr1(axis1)**2) \
  !           -is*(rsxy1x22(i1,i2,i3,axis1,0)+rsxy1y22(i1,i2,i3,axis1,1))/(2.*dr1(axis1))
  ! clap2=(rsxy2(j1,j2,j3,axis2,0)**2+rsxy2(j1,j2,j3,axis2,1)**2)/(dr2(axis2)**2) \
  !             -js*(rsxy2x22(j1,j2,j3,axis2,0)+rsxy2y22(j1,j2,j3,axis2,1))/(2.*dr2(axis2)) 
  clap1=(rsxy1(i1,i2,i3,axis1,0)**2+rsxy1(i1,i2,i3,axis1,1)**2)/(dr1(axis1)**2) \
            -is*(rxx1(axis1,0,0)+rxx1(axis1,1,1))/(2.*dr1(axis1))
  clap2=(rsxy2(j1,j2,j3,axis2,0)**2+rsxy2(j1,j2,j3,axis2,1)**2)/(dr2(axis2)**2) \
            -js*(rxx2(axis2,0,0)+rxx2(axis2,1,1))/(2.*dr2(axis2)) 

  ! ---- EQUATION 1 ----- 
  !   [ n.(uv.xx + u.yy)/mu ] = 0
  a4(1,0) = an1*clap1/mu1
  a4(1,1) = an2*clap1/mu1
  a4(1,2) =-an1*clap2/mu2
  a4(1,3) =-an2*clap2/mu2

  ! ---- EQUATION 3 ----- 
  !   [ tau.(uv.xx+uv.yy)*beta/(eps*mu) + ... ] = 0
  a4(3,0) = tau1*clap1*beta1/epsmu1
  a4(3,1) = tau2*clap1*beta1/epsmu1
  a4(3,2) =-tau1*clap2*beta2/epsmu2
  a4(3,3) =-tau2*clap2*beta2/epsmu2
    

   if( .false. .or. debug.gt.4 )then 
     write(*,'("BEFORE: --> i1,i2=",2i4," j1,j2=",2i4," f()=",4e10.2)') i1,i2,j1,j2,f(0),f(1),f(2),f(3)
     write(*,'("     beta1,beta2=",2e10.2," fp1=",2e10.2," fp2=",2e10.2)') beta1,beta2,fp1(0),fp1(1),fp2(0),fp2(1)
     write(*,'("     mu1,mu2=",2e10.2," v1y,u1x,v2y,u2x=",4e10.2)') mu1,mu2,v1y,u1x,v2y,u2x
   end if

  q(0) = u1(i1-is1,i2-is2,i3,ex)
  q(1) = u1(i1-is1,i2-is2,i3,ey)
  q(2) = u2(j1-js1,j2-js2,j3,ex)
  q(3) = u2(j1-js1,j2-js2,j3,ey)

  ! --- check matrix coefficients by delta function approach ----
  if( checkCoeff.eq.1 )then
    numberOfEquations=4
    checkCoefficients(i1,i2,i3, j1,j2,j3,numberOfEquations,a4,evalDispersiveInterfaceEquations22c )
  end if

  ! write(debugFile,'(" --> xy1=",4f8.3)') xy1(i1,i2,i3,0),xy1(i1,i2,i3,1)
  ! write(debugFile,'(" --> rsxy1=",4f8.3)') rsxy1(i1,i2,i3,0,0),rsxy1(i1,i2,i3,1,0),rsxy1(i1,i2,i3,0,1),rsxy1(i1,i2,i3,1,1)
  ! write(debugFile,'(" --> rsxy2=",4f8.3)') rsxy2(j1,j2,j3,0,0),rsxy2(j1,j2,j3,1,0),rsxy2(j1,j2,j3,0,1),rsxy2(j1,j2,j3,1,1)

  ! write(debugFile,'(" --> rxx1=",2f8.3)') rxx1(axis1,0,0),rxx1(axis1,1,1)
  ! write(debugFile,'(" --> rxx2=",2f8.3)') rxx2(axis2,0,0),rxx2(axis1,1,1)

  ! write(debugFile,'(" --> a4(0,.)=",4f8.3)') a4(0,0),a4(0,1),a4(0,2),a4(0,3)
  ! write(debugFile,'(" --> a4(1,.)=",4f8.3)') a4(1,0),a4(1,1),a4(1,2),a4(1,3)
  ! write(debugFile,'(" --> a4(2,.)=",4f8.3)') a4(2,0),a4(2,1),a4(2,2),a4(2,3)
  ! write(debugFile,'(" --> a4(3,.)=",4f8.3)') a4(3,0),a4(3,1),a4(3,2),a4(3,3)
  ! write(debugFile,'(" --> an1,an2=",2f8.3)') an1,an2
  ! write(debugFile,'(" --> clap1,clap2=",2f8.3)') clap1,clap2
  ! subtract off the contributions from the wrong values at the ghost points:
  do n=0,3
    f(n) = (a4(n,0)*q(0)+a4(n,1)*q(1)+a4(n,2)*q(2)+a4(n,3)*q(3)) - f(n)
  end do
  ! write(debugFile,'(" --> order2-curv: i1,i2=",2i4," f(subtract)=",4f8.3)') i1,i2,f(0),f(1),f(2),f(3)
  ! solve A Q = F
  ! factor the matrix
  numberOfEquations=4
  call dgeco( a4(0,0), numberOfEquations, numberOfEquations, ipvt(0),rcond,work(0))
  ! solve
  ! write(debugFile,'(" --> order2-curv: i1,i2=",2i4," rcond=",e10.2)') i1,i2,rcond
  job=0
  call dgesl( a4(0,0), numberOfEquations, numberOfEquations, ipvt(0), f(0), job)
  ! write(debugFile,'(" --> order2-curv: i1,i2=",2i4," f(solve)=",4f8.3)') i1,i2,f(0),f(1),f(2),f(3)

  u1(i1-is1,i2-is2,i3,ex)=f(0)
  u1(i1-is1,i2-is2,i3,ey)=f(1)
  u2(j1-js1,j2-js2,j3,ex)=f(2)
  u2(j1-js1,j2-js2,j3,ey)=f(3)

  if( .false. .or. debug.gt.3 )then ! re-evaluate
    evalInterfaceDerivatives2d()
    eval2dJumpDispersiveOrder2()
    !write(debugFile,'(" --> order2-curv: xy1(ghost)=",2e11.3)') xy1(i1-is1,i2-is2,i3,0),xy1(i1-is1,i2-is2,i3,1)
    !write(debugFile,'(" --> order2-curv: xy2(ghost)=",2e11.3)') xy2(j1-js1,j2-js2,j3,0),xy2(j1-js1,j2-js2,j3,1)

    write(*,'("AFTER: --> i1,i2=",2i4," j1,j2=",2i4," f(re-eval)=",4e10.2)') i1,i2,j1,j2,f(0),f(1),f(2),f(3)

    if( twilightZone.eq.1 )then
      ! check errors in the ghost 
        k1=i1-is1
        k2=i2-is2
        k3=i3
        do n=0,nd-1
          call ogderiv(ep, 0,0,0,0, xy1(k1,k2,k3,0),xy1(k1,k2,k3,1),0.,t,ex+n, es(n)   ) 
        end do
        f(0) =  u1(i1-is1,i2-is2,i3,ex) -es(0)
        f(1) =  u1(i1-is1,i2-is2,i3,ey) -es(1)
        k1=j1-js1
        k2=j2-js2
        k3=j3
        do n=0,nd-1
          call ogderiv(ep, 0,0,0,0, xy2(k1,k2,k3,0),xy2(k1,k2,k3,1),0.,t,ex+n, est(n)   ) 
        end do
        f(2) =  u2(j1-js1,j2-js2,j3,ex) -est(0)
        f(3) =  u2(j1-js1,j2-js2,j3,ey) -est(1)
        write(*,'(" ghost err =",4e10.2)') f(0),f(1),f(2),f(3) 

    end if

  end if

  ! -- Hz has already been filled in by extrapolation ----

endLoopsMask2d()

 if( checkCoeff.eq.1 )then
   write(*,'("+++++ iGDM22c: check coeff in interface: max(diff) = ",1pe8.2)') coeffDiff
 end if


#endMacro


! --------------------------------------------------------------------------------------------
! Macro:  Evaluate the RHS to the jump conditons: 2D, Order=2, Dispersive
!
! Here are the jump conditions (See notes in DMX_ADE)
!   (1) [ div(E) ] = 0
!   (2) [ (1/mu)* nv.( Delta(E) ) ]=0
!   (3) [ (1/mu)* tv.( curl(E) ) ]=0               -->   [ (1/mu)* \nv\times( curl(E) ) ]=0
!   (4) [ tv.(c^2*Delta(E) -alphaP*P_tt) ] = 0    -->   [ \nv X ( c^2*Delta(E) -alphaP*P_tt) ] = 0 
! 
! These 6 equations can be written as 
!   [ div(E) n + (I- n n^T)( curl(E)/mu ) ] =0                                 (3 eqns)
!   [ (1/mu) n n^T Delta(E) + (I-n n^T) ( c^2*Delta(E) -alphaP*P_tt ) ] = 0    (3 eqns)
!
! --------------------------------------------------------------------------------------------
#beginMacro eval3dJumpDispersiveOrder2()

 ! stop 6767

 ! 2D: 
 ! f(0)=(u1x+v1y) - \
 !      (u2x+v2y)
 ! f(1)=( an1*u1Lap +an2*v1Lap )/mu1 - \
 !      ( an1*u2Lap +an2*v2Lap )/mu2 
 ! f(2)=(v1x-u1y)/mu1 - \
 !      (v2x-u2y)/mu2
 ! f(3)=( ( tau1*u1Lap +tau2*v1Lap )*beta1/epsmu1 - alphaP1*(tau1*fp1(0)+tau2*fp1(1)) ) - \
 !      ( ( tau1*u2Lap +tau2*v2Lap )*beta2/epsmu2 - alphaP2*(tau1*fp2(0)+tau2*fp2(1)) )

 divE1 = u1x+v1y+w1z
 curlE1x = w1y-v1z
 curlE1y = u1z-w1x
 curlE1z = v1x-u1y
 nDotCurlE1=an1*curlE1x+an2*curlE1y+an3*curlE1z
 nDotLapE1 = an1*u1Lap + an2*v1Lap + an3*w1Lap

 divE2 = u2x+v2y+w2z
 curlE2x = w2y-v2z
 curlE2y = u2z-w2x
 curlE2z = v2x-u2y
 nDotCurlE2=an1*curlE2x+an2*curlE2y+an3*curlE2z
 nDotLapE2 = an1*u2Lap + an2*v2Lap + an3*w2Lap


 f(0)=( divE1*an1 + (curlE1x- nDotCurlE1*an1)/mu1 ) - ( divE2*an1 + (curlE2x- nDotCurlE2*an1)/mu2 )
 f(1)=( divE1*an2 + (curlE1y- nDotCurlE1*an2)/mu1 ) - ( divE2*an2 + (curlE2y- nDotCurlE2*an2)/mu2 )
 f(2)=( divE1*an3 + (curlE1z- nDotCurlE1*an3)/mu1 ) - ( divE2*an3 + (curlE2z- nDotCurlE2*an3)/mu2 )

 ! cem1 = (1-beta1/eps1)/mu1 
 ! cem2 = (1-beta2/eps2)/mu2 
 ! beta = 1 - alphaP*Sum( C_k )
 !

 nDotFp1 = an1*fp1(0)+an2*fp1(1)+an3*fp1(2)
 nDotFp2 = an1*fp2(0)+an2*fp2(1)+an3*fp2(2)
 f(3)=( u1Lap*beta1/(epsmu1) + cem1*nDotLapE1*an1 - alphaP1*( fp1(0)-an1*nDotFp1) ) - \
      ( u2Lap*beta2/(epsmu2) + cem2*nDotLapE2*an1 - alphaP2*( fp2(0)-an1*nDotFp2) )

 f(4)=( v1Lap*beta1/(epsmu1) + cem1*nDotLapE1*an2 - alphaP1*( fp1(1)-an2*nDotFp1) ) - \
      ( v2Lap*beta2/(epsmu2) + cem2*nDotLapE2*an2 - alphaP2*( fp2(1)-an2*nDotFp2) )

 f(5)=( w1Lap*beta1/(epsmu1) + cem1*nDotLapE1*an3 - alphaP1*( fp1(2)-an3*nDotFp1) ) - \
      ( w2Lap*beta2/(epsmu2) + cem2*nDotLapE2*an3 - alphaP2*( fp2(2)-an3*nDotFp2) )

 ! f(3)=( u1Lap/(epsmu1) + cem1*nDotLapE1*an1 ) - ( u2Lap/(epsmu2) + cem2*nDotLapE2*an1 )
 ! f(4)=( v1Lap/(epsmu1) + cem1*nDotLapE1*an2 ) - ( v2Lap/(epsmu2) + cem2*nDotLapE2*an2 )
 ! f(5)=( w1Lap/(epsmu1) + cem1*nDotLapE1*an3 ) - ( w2Lap/(epsmu2) + cem2*nDotLapE2*an3 )
 if( twilightZone.eq.1 )then

   ! For now we assume mu1=mu2 and TZ solutions are the same on both sides.

   ! TZ forcing goes into equation (4)
   !   (4) [ tv.(c^2*Delta(E) -alphaP*P_tt + fev) ] = 0 
   !
   !   [ (1/mu) n n^T Delta(E) + (I-n n^T) ( c^2*Delta(E) -alphaP*P_tt + fev ) ] = 0    (3 eqns)

  
   ! -- add on the jump in the forcing ---
   f(3) = f(3) + (1.-an1*an1)*(fev1(0)-fev2(0)) + (0.-an1*an2)*(fev1(1)-fev2(1)) + (0.-an1*an3)*(fev1(2)-fev2(2)) 
   f(4) = f(4) + (0.-an2*an1)*(fev1(0)-fev2(0)) + (1.-an2*an2)*(fev1(1)-fev2(1)) + (0.-an2*an3)*(fev1(2)-fev2(2)) 
   f(5) = f(5) + (0.-an3*an1)*(fev1(0)-fev2(0)) + (0.-an3*an2)*(fev1(1)-fev2(1)) + (1.-an3*an3)*(fev1(2)-fev2(2)) 

   ! write(debugFile,'(" u1Lap,ueLap=",2e10.2," v1Lap,veLap=",2e10.2)') u1Lap,ueLap,v1Lap,veLap

 end if


! write(*,'("FINISH ME JUMP 3D GDM..., beta1,beta2=",2(1pe10.2)," nDotFp1,nDotFp2=",2(1pe10.2))') beta1,beta2,nDotFp1,nDotFp2

#endMacro



! ----------------------------------------------------------------------------------
!  Macro:
!    Evaluate the interface equations for checking the coefficients
! ----------------------------------------------------------------------------------
#beginMacro evalDispersiveInterfaceEquations23c()

  evalInterfaceDerivatives3d()
  eval3dJumpDispersiveOrder2()

#endMacro


! --------------------------------------------------------------------------
! Macro: Assign interface ghost values, DIM=3, ORDER=2, GRID=Curvilinear
! 
!                  DISPERSIVE CASE
!
! Here are the jump conditions (See notes in DMX_ADE)
!   (1) [ div(E) ] = 0
!   (2) [ (1/mu)* nv.( Delta(E) ) ]=0
!   (3) [ (1/mu)* tv.( curl(E) ) ]=0               -->   [ (1/mu)* \nv\times( curl(E) ) ]=0
!   (4) [ tv.(c^2*Delta(E) -alphaP*P_tt) ] = 0    -->   [ \nv X ( c^2*Delta(E) -alphaP*P_tt) ] = 0 
! 
! These 6 equations can be written as 
!   [ div(E) n + (I- n n^T)( curl(E)/mu ) ] =0                                 (3 eqns)
!   [ (1/mu) n n^T Delta(E) + (I-n n^T) ( c^2*Delta(E) -alphaP*P_tt ) ] = 0    (3 eqns)
!
! An approximation to P_tt takes the form 
!   P_tt = K Delta(E) + G(E,P)
! --------------------------------------------------------------------------
#beginMacro assignDispersiveInterfaceGhost23c()

  INFO("23c-GDM")

  ! write(*,'(" FINISH ME for non-zero GDM ")')
  ! stop 9876

  ! --- initialize some forcing functions ---
  do n=0,nd-1
    fev1(n)=0.
    fev2(n)=0.
    do jv=0,numberOfPolarizationVectors1-1
      fpv1(n,jv)=0.
    end do
    do jv=0,numberOfPolarizationVectors2-1
      fpv2(n,jv)=0.
    end do
  end do

  beginLoopsMask3d()

    ! here is the normal (assumed to be the same on both sides)
    an1=rsxy1(i1,i2,i3,axis1,0)   ! normal (an1,an2)
    an2=rsxy1(i1,i2,i3,axis1,1)
    an3=rsxy1(i1,i2,i3,axis1,2)
    aNorm=max(epsx,sqrt(an1**2+an2**2+an3**2))
    an1=an1/aNorm
    an2=an2/aNorm
    an3=an3/aNorm


    ! --- first evaluate the equations we want to solve with the wrong values at the ghost points:


    evalInterfaceDerivatives3d()

    ! Evaluate TZ forcing for dispersive equations in 3D 
    getDispersiveTZForcing3d(fpv1,fpv2,fev1,fev2)

    ! eval dispersive forcings for domain 1
    getDispersiveForcingOrder2(i1,i2,i3, fp1, fpv1,fev1,p1,p1n,p1m, u1,u1n,u1m, dispersionModel1,numberOfPolarizationVectors1,alphaP1,beta1,a0v1,a1v1,b0v1,b1v1)

    ! eval dispersive forcings for domain 2
    getDispersiveForcingOrder2(j1,j2,j3, fp2, fpv2,fev2,p2,p2n,p2m, u2,u2n,u2m, dispersionModel2,numberOfPolarizationVectors2,alphaP2,beta2,a0v2,a1v2,b0v2,b1v2)

    cem1=(1.-beta1/eps1)/mu1
    cem2=(1.-beta2/eps2)/mu2

    betac1=beta1/epsmu1
    betac2=beta2/epsmu2


    eval3dJumpDispersiveOrder2()

    if( debug.gt.4 )then
     write(debugFile,'(" --> 3d-order2-curv: i1,i2,i3=",3i4," f(start)=",6f8.3)') i1,i2,i3,f(0),f(1),f(2),f(3),f(4),f(5)
     ! '
     write(debugFile,'(" --> u1x,u1y,u1z,v1x,v1y,v1z=",6f8.4)') u1x,u1y,u1z,v1x,v1y,v1z
     write(debugFile,'(" --> u2x,u2y,u2z,v2x,v2y,v2z=",6f8.4)') u2x,u2y,u2z,v2x,v2y,v2z

     write(debugFile,'(" --> vv1r,vv1s,vv1t         =",3e9.2)') vv1r,vv1s,vv1t
     do k3=-1,1
     do k2=-1,1
     write(debugFile,'(" --> v1: =",3f8.4)') u1(i1-1,i2+k2,i3+k3,ey),u1(i1,i2+k2,i3+k3,ey),u1(i1+1,i2+k2,i3+k3,ey)
     end do
     end do
     do k3=-1,1
     do k2=-1,1
     write(debugFile,'(" --> v2: =",3f8.4)') u2(j1-1,j2+k2,j3+k3,ey),u2(j1,j2+k2,j3+k3,ey),u2(j1+1,j2+k2,j3+k3,ey)
     end do
     end do
     ! '
    end if

    ! here is the matrix of coefficients for the unknowns u1(-1),v1(-1),w1(-1),  u2(-1),v2(-1),w2(-1)
    ! Solve:
    !     
    !       A [ U ] = A [ U(old) ] - [ f ]

    c1x = -is*rsxy1(i1,i2,i3,axis1,0)/(2.*dr1(axis1))    ! coeff of u1(-1) from D.x
    c1y = -is*rsxy1(i1,i2,i3,axis1,1)/(2.*dr1(axis1))    ! coeff of u1(-1) from D.y 
    c1z = -is*rsxy1(i1,i2,i3,axis1,2)/(2.*dr1(axis1))    ! coeff of u1(-1) from D.z

    c2x = -js*rsxy2(j1,j2,j3,axis2,0)/(2.*dr2(axis2))
    c2y = -js*rsxy2(j1,j2,j3,axis2,1)/(2.*dr2(axis2))
    c2z = -js*rsxy2(j1,j2,j3,axis2,2)/(2.*dr2(axis2))

    rxx1(0,0,0)=aj1rxx
    rxx1(0,1,1)=aj1ryy
    rxx1(0,2,2)=aj1rzz
    rxx1(1,0,0)=aj1sxx
    rxx1(1,1,1)=aj1syy
    rxx1(1,2,2)=aj1szz
    rxx1(2,0,0)=aj1txx
    rxx1(2,1,1)=aj1tyy
    rxx1(2,2,2)=aj1tzz

    rxx2(0,0,0)=aj2rxx
    rxx2(0,1,1)=aj2ryy
    rxx2(0,2,2)=aj2rzz
    rxx2(1,0,0)=aj2sxx
    rxx2(1,1,1)=aj2syy
    rxx2(1,2,2)=aj2szz
    rxx2(2,0,0)=aj2txx
    rxx2(2,1,1)=aj2tyy
    rxx2(2,2,2)=aj2tzz

    ! clap1 : coeff of u(-1) from lap = u.xx + u.yy + u.zz

    ! clap1=(rsxy1(i1,i2,i3,axis1,0)**2+rsxy1(i1,i2,i3,axis1,1)**2)/(dr1(axis1)**2) \
    !           -is*(rsxy1x22(i1,i2,i3,axis1,0)+rsxy1y22(i1,i2,i3,axis1,1))/(2.*dr1(axis1))
    ! clap2=(rsxy2(j1,j2,j3,axis2,0)**2+rsxy2(j1,j2,j3,axis2,1)**2)/(dr2(axis2)**2) \
    !             -js*(rsxy2x22(j1,j2,j3,axis2,0)+rsxy2y22(j1,j2,j3,axis2,1))/(2.*dr2(axis2)) 
    clap1=(rsxy1(i1,i2,i3,axis1,0)**2+rsxy1(i1,i2,i3,axis1,1)**2+rsxy1(i1,i2,i3,axis1,2)**2)/(dr1(axis1)**2) \
              -is*(rxx1(axis1,0,0)+rxx1(axis1,1,1)+rxx1(axis1,2,2))/(2.*dr1(axis1))
    clap2=(rsxy2(j1,j2,j3,axis2,0)**2+rsxy2(j1,j2,j3,axis2,1)**2+rsxy2(j1,j2,j3,axis2,2)**2)/(dr2(axis2)**2) \
              -js*(rxx2(axis2,0,0)+rxx2(axis2,1,1)+rxx2(axis2,2,2))/(2.*dr2(axis2)) 

    ! cdivE1 =  u.c1x + v.c1y + w.c1z
    ! nDotCurlE1 = (w1y-v1z)*an1 + (u1z-w1x)*an2 + (v1x-u1y)*an3

    !  f(0)=( divE1*an1 + (curlE1x- nDotCurlE1*an1)/mu1 ) - ( divE2*an1 + (curlE2x- nDotCurlE2*an1)/mu2 )
    !  f(1)=( divE1*an2 + (curlE1y- nDotCurlE1*an2)/mu1 ) - ( divE2*an2 + (curlE2y- nDotCurlE2*an2)/mu2 )
    !  f(2)=( divE1*an3 + (curlE1z- nDotCurlE1*an3)/mu1 ) - ( divE2*an3 + (curlE2z- nDotCurlE2*an3)/mu2 )

    ! ---- EQUATION 0 ----- 
    ! (u.x+v.y+w.z)*an1 + ( w1y-v1z - nDotCurlE1*an1)/mu1
    a6(0,0) = ( c1x*an1 + (         - (c1z*an2-c1y*an3)*an1 )/mu1 ) ! coeff of u1(-1)
    a6(0,1) = ( c1y*an1 + (    -c1z - (c1x*an3-c1z*an1)*an1 )/mu1 ) ! coeff of v1(-1)
    a6(0,2) = ( c1z*an1 + ( c1y     - (c1y*an1-c1x*an2)*an1 )/mu1 ) ! coeff of w1(-1)

    a6(0,3) =-( c2x*an1 + (         - (c2z*an2-c2y*an3)*an1 )/mu2 ) ! coeff of u2(-1)
    a6(0,4) =-( c2y*an1 + (    -c2z - (c2x*an3-c2z*an1)*an1 )/mu2 ) ! coeff of v2(-1)
    a6(0,5) =-( c2z*an1 + ( c2y     - (c2y*an1-c2x*an2)*an1 )/mu2 ) ! coeff of w2(-1)

    ! ---- EQUATION 1 ----- 
    ! (u.x+v.y+w.z)*an2 + ( u1z-w1x - nDotCurlE1*an2)/mu1
    a6(1,0) = ( c1x*an2 + ( c1z     - (c1z*an2-c1y*an3)*an2 )/mu1 ) ! coeff of u1(-1)
    a6(1,1) = ( c1y*an2 + (         - (c1x*an3-c1z*an1)*an2 )/mu1 ) ! coeff of v1(-1)
    a6(1,2) = ( c1z*an2 + (    -c1x - (c1y*an1-c1x*an2)*an2 )/mu1 ) ! coeff of w1(-1)

    a6(1,3) =-( c2x*an2 + ( c2z     - (c2z*an2-c2y*an3)*an2 )/mu2 ) ! coeff of u2(-1)
    a6(1,4) =-( c2y*an2 + (         - (c2x*an3-c2z*an1)*an2 )/mu2 ) ! coeff of v2(-1)
    a6(1,5) =-( c2z*an2 + (    -c2x - (c2y*an1-c2x*an2)*an2 )/mu2 ) ! coeff of w2(-1)

    ! ---- EQUATION 2 ----- 
    ! (u.x+v.y+w.z)*an3 + ( v1x-u1y - nDotCurlE1*an2)/mu1
    a6(2,0) = ( c1x*an3 + (    -c1y - (c1z*an2-c1y*an3)*an3 )/mu1 ) ! coeff of u1(-1)
    a6(2,1) = ( c1y*an3 + ( c1x     - (c1x*an3-c1z*an1)*an3 )/mu1 ) ! coeff of v1(-1)
    a6(2,2) = ( c1z*an3 + (         - (c1y*an1-c1x*an2)*an3 )/mu1 ) ! coeff of w1(-1)

    a6(2,3) =-( c2x*an3 + (    -c2y - (c2z*an2-c2y*an3)*an3 )/mu2 ) ! coeff of u2(-1)
    a6(2,4) =-( c2y*an3 + ( c2x     - (c2x*an3-c2z*an1)*an3 )/mu2 ) ! coeff of v2(-1)
    a6(2,5) =-( c2z*an3 + (         - (c2y*an1-c2x*an2)*an3 )/mu2 ) ! coeff of w2(-1)

    ! betac1=beta1/epsmu1
    ! betac2=beta2/epsmu2

    ! ---- EQUATION 3 ----- 
    !    [ (uv.xx+uv.yy)*beta/(eps*mu) + ... ] = 0
    !  u1Lap*beta1/(epsmu1) + cem1*( an1*u1Lap + an2*v1Lap + an3*w1Lap )*an1
    a6(3,0) = ( clap1*betac1   + cem1*( an1*clap1                         )*an1 ) ! coeff of u1(-1)
    a6(3,1) = (                  cem1*(             an2*clap1             )*an1 ) ! coeff of v1(-1)
    a6(3,2) = (                  cem1*(                         an3*clap1 )*an1 ) ! coeff of w1(-1)

    a6(3,3) =-( clap2*betac2   + cem2*( an1*clap2                         )*an1 ) ! coeff of u2(-1)
    a6(3,4) =-(                  cem2*(             an2*clap2             )*an1 )
    a6(3,5) =-(                  cem2*(                         an3*clap2 )*an1 )

    ! ---- EQUATION 4 ----- 
    !  v1Lap*beta1/(epsmu1) + cem1*( an1*u1Lap + an2*v1Lap + an3*w1Lap )*an2
    a6(4,0) = (                  cem1*( an1*clap1                         )*an2 ) ! coeff of u1(-1)
    a6(4,1) = ( clap1*betac1   + cem1*(             an2*clap1             )*an2 )
    a6(4,2) = (                  cem1*(                         an3*clap1 )*an2 )

    a6(4,3) =-(                  cem2*( an1*clap2                         )*an2 ) ! coeff of u2(-1)
    a6(4,4) =-( clap2*betac2   + cem2*(             an2*clap2             )*an2 )
    a6(4,5) =-(                  cem2*(                         an3*clap2 )*an2 )

    ! ---- EQUATION 5 ----- 
    !  w1Lap*beta1/(epsmu1) + cem1*( an1*u1Lap + an2*v1Lap + an3*w1Lap )*an3
    a6(5,0) = (                  cem1*( an1*clap1                         )*an3 ) ! coeff of u1(-1)
    a6(5,1) = (                  cem1*(             an2*clap1             )*an3 )
    a6(5,2) = ( clap1*betac1   + cem1*(                         an3*clap1 )*an3 )

    ! ---- EQUATION 6 ----- 
    a6(5,3) =-(                  cem2*( an1*clap2                         )*an3 ) ! coeff of u2(-1)
    a6(5,4) =-(                  cem2*(             an2*clap2             )*an3 )
    a6(5,5) =-( clap2*betac2   + cem2*(                         an3*clap2 )*an3 )


    ! --- check matrix coefficients by delta function approach ----
    if( checkCoeff.eq.1 )then
      numberOfEquations=6
      checkCoefficients(i1,i2,i3, j1,j2,j3,numberOfEquations,a6,evalDispersiveInterfaceEquations23c )
    end if

    q(0) = u1(i1-is1,i2-is2,i3-is3,ex)
    q(1) = u1(i1-is1,i2-is2,i3-is3,ey)
    q(2) = u1(i1-is1,i2-is2,i3-is3,ez)
    q(3) = u2(j1-js1,j2-js2,j3-js3,ex)
    q(4) = u2(j1-js1,j2-js2,j3-js3,ey)
    q(5) = u2(j1-js1,j2-js2,j3-js3,ez)

    ! subtract off the contributions from the wrong values at the ghost points:
    do n=0,5
      f(n) = (a6(n,0)*q(0)+a6(n,1)*q(1)+a6(n,2)*q(2)+a6(n,3)*q(3)+a6(n,4)*q(4)+a6(n,5)*q(5)) - f(n)
    end do
!-    if( .true. )then 
!-      do n=0,5
!-        if( isnan(f(n)) )then
!-         write(*,'(" NAN found!")') 
!-         write(*,'(" i1,i2,i3=",3i6)') i1,i2,i3
!-         write(*,'(" j1,j2,j3=",3i6)') j1,j2,j3
!-         write(*,'(" f=",6(e10.2,1x))') (f(m),m=0,5)
!-         write(*,'(" q=",6(e10.2,1x))') (q(m),m=0,5)
!-         write(*,'(" a6=",6(e10.2,1x))') ((a6(n1,n2),n1=0,5),n2=0,5)
!-         write(*,'(" clap1,clap2=",2(e10.2,1x))') clap1,clap2
!-         write(*,'(" betac1,betac2=",2(e10.2,1x))') betac1,betac2
!-         write(*,'(" cem1,cem2=",2(e10.2,1x))') cem1,cem2
!-         write(*,'(" c1x,c1y,c1z=",3(e10.2,1x))') c1x,c1y,c1z
!-         write(*,'(" c2x,c2y,c2z=",3(e10.2,1x))') c2x,c2y,c2z
!-         stop 6666
!-        end if
!-      end do
!-    end if

    ! write(debugFile,'(" --> 3d:order2-c: f(subtract)=",6f8.3)') f(0),f(1),f(2),f(3),f(4),f(5)
    ! solve A Q = F
    ! factor the matrix
    numberOfEquations=6
    call dgeco( a6(0,0), numberOfEquations, numberOfEquations, ipvt(0),rcond,work(0))
    ! solve
    ! write(debugFile,'(" --> 3d:order2-c: rcond=",e10.2)') rcond
    job=0
    call dgesl( a6(0,0), numberOfEquations, numberOfEquations, ipvt(0), f(0), job)
    ! write(debugFile,'(" --> 3d:order2-c: f(solve)=",6f8.3)') f(0),f(1),f(2),f(3),f(4),f(5)
    ! write(debugFile,'(" --> 3d:order2-c:        q=",6f8.3)') q(0),q(1),q(2),q(3),q(4),q(5)

    u1(i1-is1,i2-is2,i3-is3,ex)=f(0)
    u1(i1-is1,i2-is2,i3-is3,ey)=f(1)
    u1(i1-is1,i2-is2,i3-is3,ez)=f(2)
    u2(j1-js1,j2-js2,j3-js3,ex)=f(3)
    u2(j1-js1,j2-js2,j3-js3,ey)=f(4)
    u2(j1-js1,j2-js2,j3-js3,ez)=f(5)

    if( .false. )then
    u1(i1-is1,i2-is2,i3-is3,ex)=q(0)
    u1(i1-is1,i2-is2,i3-is3,ey)=q(1)
    u1(i1-is1,i2-is2,i3-is3,ez)=q(2)
    u2(j1-js1,j2-js2,j3-js3,ex)=q(3)
    u2(j1-js1,j2-js2,j3-js3,ey)=q(4)
    u2(j1-js1,j2-js2,j3-js3,ez)=q(5)
    end if

    if( debug.gt.3 )then ! re-evaluate
     evalInterfaceDerivatives3d()
     eval3dJumpOrder2()
     write(debugFile,'(" --> 3d-order2-c: i1,i2,i3=",3i4," f(re-eval)=",6e10.2)') i1,i2,i3,f(0),f(1),f(2),f(3),f(4),f(5)
       ! '
    end if

  endLoopsMask3d()

 if( checkCoeff.eq.1 )then
   write(*,'("+++++ iGDM23c: check coeff in interface: max(diff) = ",1pe8.2)') coeffDiff
 end if

#endMacro         



! -------------------------------------------------------------------------
! Macro: Evaluate DISPERSIVE forcing terms, FOURTH-ORDER
!   This macro can be usedto eval values in either domain 1 or domain 2
!
! Input:
!   FACE : LEFT or : RIGHT
!   fev(n) : forcing on E equation: E_{tt} = c^2 Delta(E) + ... + fev
!   fpv(n,jv) : forcing on equation for P_{n,jv} 
! Output
!   fp(n) : 
!   c2PttLEsum   : coeff of L(E) in P.tt     (second-order)
!   c4PttLEsum   : coeff of L(E) in P.tt     (fourth-order)
!   c4PttLLEsum  : coeff of L*L(E) in P.tt 
! ------------------------------------------------------------------------
#beginMacro getDispersiveForcingOrder4(FACE,k1,k2,k3, fp, fpv,fev, p,pn,pm, u,un,um, dispersionModel,numberOfPolarizationVectors,alphaP,\
            c2PttEsum,c2PttLEsum,c4PttLEsum,c4PttLLEsum,c2PttttLEsum,c2PttttLLEsum,a0v,a1v,b0v,b1v,LE,LLE,LEm,LfE,LfP,fEt,fEtt,fPt,fPtt,pevtt,pevttx,pevtty,pevtttt,evx,evy,evnx,evny,\
            fevx,fevy,fpvx,fpvy,LEx,LEy,fPttx,fPtty,fLPtt,fPtttt)


 do n=0,nd-1
   fp(n)=0.
   fPttx(n) =0.
   fPtty(n) =0.
   fLPtt(n) =0.
   fPtttt(n)=0.
 end do
 
 c2PttEsum=0.
 c2PttLEsum=0.
 c4PttLEsum=0.
 c4PttLLEsum=0.
 
 c2PttttLEsum=0.
 c2PttttLLEsum=0.
 
 do jv=0,numberOfPolarizationVectors-1
   a0=a0v(jv)
   a1=a1v(jv)
   b0=b0v(jv)
   b1=b1v(jv)
   alpha=alphaP
 
   ! Second-order coefficients: 
   ! Ptt = c2PttLE*LE + c2PttE*E + c2PttEm*Em + c2PttP*P + c2PttPm*Pm + c2PttfE*fE + c2PttfP*fP
   #Include interfaceAdeGdmOrder2.h 
 
   ! Fourth-order coefficients
   ! Ptt = c4PttLE*LE + c4PttE*E + c4PttEm*Em + c4PttP*P + c4PttPm*Pm + c4PttfE*fE + c4PttfP*fP
   !    + c4PttLLE*LLE + c4PttLP*LP + c4PttLEm*LEm + c4PttLPm*LPm+ c4PttLfE*LfE + c4PttLfP*LfP
   !    + c4PttfEt*fEt + c4PttfEtt*fEtt + c4PttfPt*fPt+ c4PttfPtt*fPtt
   ! #Include interfaceAdeGdmOrder4.h 
   ! Nov 4, 2018 *new way* 
   #Include interfaceAdeGdmOrder4New.h 
 
   if( .false. .and. twilightZone.eq.1 )then
     write(*,'(" FACE: alpha,dt,d4,d1=",4e12.4)') alpha,dt,d4,d1
     write(*,'(" a0,a1,b0,b1=",4e12.4)') a0,a1,b0,b1
     write(*,'(" c2EtE,c2PtE,c2PttfP=",3e12.4)') c2EtE,c2PtE,c2PttfP
     write(*,'(" c4PttLE,c4PttE,c4PttEm,c4PttP,c4PttPm,c4PttfE=",6e12.4)') c4PttLE,c4PttE,c4PttEm,c4PttP,c4PttPm,c4PttfE
     write(*,'(" c4PttfP,c4PttLLE,c4PttLP,c4PttLEm,c4PttLPm,c4PttLfE,c4PttLfP=",7e12.4)') c4PttfP,c4PttLLE,c4PttLP,c4PttLEm,c4PttLPm,c4PttLfE,c4PttLfP
     write(*,'(" c4PttfEt,c4PttfEtt,c4PttfPt,c4PttfPtt=",4e12.4)') c4PttfEt,c4PttfEtt,c4PttfPt,c4PttfPtt
   end if
 
 
   ! Coeff of E in P.tt (4th order)
   c2PttEsum = c2PttEsum + c2PttE

   ! Coeff of LE in P.tt (4th order)
   c2PttLEsum  = c2PttLEsum + c2PttLE
 
   ! Coeff of LE in P.tt (4th order)
   c4PttLEsum  = c4PttLEsum + c4PttLE
 
   ! Coeff of LLE in P.tt
   c4PttLLEsum = c4PttLLEsum + c4PttLLE
 
   ! Coeff of LE and LLE in P.tttt
   c2PttttLEsum =c2PttttLEsum +c2PttttLE
   c2PttttLLEsum=c2PttttLLEsum+c2PttttLLE
 
   do n=0,nd-1
     pc = n + jv*nd 
     ec = ex +n
 
     pv   =  p(k1,k2,k3,pc)
     pvn  =  pn(k1,k2,k3,pc)
 
     ev    =  u(k1,k2,k3,ec)
     evn   =  un(k1,k2,k3,ec)
 
     ! Left: u1x,u1y, u1xx, u1yy, u1Lap (ex)
     !       v1x,v1y, v1xx, v1yy, v1Lap (ey) 
 
 
     ! These next derivatives may only be needed to order2, (use order 4 for testing TZ polynomials)
     #perl $ORDER=4;
     ! perl $ORDER=2;
     #If #FACE eq "LEFT"
       opEvalJacobianDerivatives(rsxy1,k1,k2,k3,aj1,1)
       ! uu1 in the next statement defines names of intermediate values
       evalSecondDerivs(rsxy1,aj1,p1,k1,k2,k3,pc,uu1,p1)
       ! write(*,'("FACE: p1x,p1y,p1xx,p1yy,p1Lap=",5(1pe12.4))') p1x,p1y,p1xx,p1yy,p1Lap
 
       LP  = (c1**2)*p1Lap
       evalSecondDerivs(rsxy1,aj1,p1n,k1,k2,k3,pc,uu1,p1n)
       ! write(*,'("FACE: p1nxx,p1nyy,p1nLap=",3e12.4)') p1nxx,p1nyy,p1nLap
       LPm = (c1**2)*p1nLap
 
       pvx  = p1x
       pvy  = p1y
       pvnx = p1nx
       pvny = p1ny
 
     #Elif #FACE eq "RIGHT"
       opEvalJacobianDerivatives(rsxy2,k1,k2,k3,aj2,1)
       ! uu1 in the next statement defines names of intermediate values
       evalSecondDerivs(rsxy2,aj2,p2,k1,k2,k3,pc,uu2,p2)
       LP  = (c2**2)*p2Lap
       evalSecondDerivs(rsxy2,aj2,p2n,k1,k2,k3,pc,uu2,p2n)
       LPm = (c2**2)*p2nLap
 
       pvx  = p2x
       pvy  = p2y
       pvnx = p2nx
       pvny = p2ny
 
     #Else
       write(*,'(" interface3d:ERROR: unknown FACE")')
       stop 7777
     #End
     
 
     ! Accumulate: SUM_m Pm,tt
     ptta = c4PttLE*LE(n) + c4PttE*ev + c4PttEm*evn + c4PttP*pv + c4PttPm*pvn + c4PttfE*fev(n) + c4PttfP*fpv(n,jv) \
        + c4PttLLE*LLE(n) + c4PttLP*LP + c4PttLEm*LEm(n) + c4PttLPm*LPm+ c4PttLfE*LfE(n) + c4PttLfP*LfP(n,jv) \
        + c4PttfEt*fEt(n) + c4PttfEtt*fEtt(n) + c4PttfPt*fPt(n,jv)+ c4PttfPtt*fPtt(n,jv)
 
     ! ---- Compute fp = P.tt 
     fp(n) = fp(n) + ptta
 
     ! ----- Compute fPttx = (P.tt).x , fPtty = (P.tt).y  (second order)
     pttxa = c2PttLE*LEx(n) + c2PttE*evx(n) + c2PttEm*evnx(n) + c2PttP*pvx + c2PttPm*pvnx + c2PttfE*fevx(n) + c2PttfP*fpvx(n,jv)
     pttya = c2PttLE*LEy(n) + c2PttE*evy(n) + c2PttEm*evny(n) + c2PttP*pvy + c2PttPm*pvny + c2PttfE*fevy(n) + c2PttfP*fpvy(n,jv)
 
     fPttx(n) = fPttx(n) + pttxa
     fPtty(n) = fPtty(n) + pttya
 
     ! ----- Compute fLPtt = L(P.tt) (second order)
     Lptta = c2PttLE*LLE(n) + c2PttE*LE(n) + c2PttEm*LEm(n) + c2PttP*LP + c2PttPm*LPm + c2PttfE*LfE(n) + c2PttfP*LfP(n,jv)
 
     fLPtt(n) = fLPtt(n) + Lptta
 
     ! ----- Compute fPtttt = P.tttt
     ptttta= c2PttttLE*LE(n) + c2PttttE*ev + c2PttttEm*evn + c2PttttP*pv + c2PttttPm*pvn + c2PttttfE*fev(n) + c2PttttfP*fpv(n,jv) \
         + c2PttttLLE*LLE(n) + c2PttttLP*LP + c2PttttLEm*LEm(n) + c2PttttLPm*LPm+ c2PttttLfE*LfE(n) + c2PttttLfP*LfP(n,jv) \
         + c2PttttfEt*fEt(n) + c2PttttfEtt*fEtt(n) + c2PttttfPt*fPt(n,jv)+ c2PttttfPtt*fPtt(n,jv)
 
     fPtttt(n) = fPtttt(n) + ptttta
 
     if( .false. .and. twilightZone.eq.1 )then
       write(*,'("")')
       write(*,'("DI4:FACE: k1,k2=",2i3," jv=",i2," n=",i2," ptta,ptte=",2e12.4)') k1,k2,jv,n,ptta,pevtt(n,jv)
       write(*,'("        : pttxa,pttxe=",2(1pe12.4)," pttya,pttye=",2(1pe12.4))') pttxa,pevttx(n,jv),pttya,pevtty(n,jv)
       write(*,'("        : ptttta,ptttte=",2(1pe12.4),/)') ptttta,pevtttt(n,jv)
 
       write(*,'(" c2PttLE,c2PttE,c2PttEm,c2PttP,c2PttPm,c2PttfE,c2PttfP=",7(1pe12.4))') c2PtLE,c2PtE,c2PtEm,c2PtP,c2PtPm,c2PtfE,c2PtfP
       write(*,'(" LEx,evx,evnx,pvx,pvnx,fevx,fpvx=",7(1pe12.4))') LEx(n),evx(n),evnx(n),pvx,pvnx,fevx(n),fpvx(n,jv)
       write(*,'(" LEy,evy,evny,pvy,pvny,fevy,fpvy=",7(1pe12.4))') LEy(n),evy(n),evny(n),pvy,pvny,fevy(n),fpvy(n,jv)
 
       ! write(*,'(" LE,LLE,LEm,LP,LPm=",5e12.4)') LE(n),LLE(n),LEm(n),LP,LPm
       ! write(*,'(" LfE,LfP,fEt,fEtt,fPt,fPtt=",6e12.4)') LfE(n),LfP(n,jv),fEt(n),fEtt(n),fPt(n,jv),fPtt(n,jv)
       ! write(*,'(" ev,evn,pv,pvn,fev,fpv=",6e12.4)')ev,evn,pv,pvn,fev(n),fpv(n,jv)
 
     end if
 
   end do ! end do n 
 
 
 end do

#endMacro


! -------------------------------------------------------------------------------
! Macro: Evaluate the TZ forcings GDM FOURTH-ORDER
! -------------------------------------------------------------------------------
#beginMacro evalTZForcingGDMOrder4(xy,i1,i2,i3,dispersionModel,numberOfPolarizationVectors,c,alphaP,a0v,a1v,b0v,b1v,fpv,fpSum,fev,\
                                   LfE,fEt,fEtt,LfP,fPt,fPtt,pevtt,pevttx,pevtty,pevtttt,fevx,fevy,fpvx,fpvy,pevttSum,pevttxSum,pevttySum,pevttLSum,pevttttSum)

if( dispersionModel.ne.noDispersion )then
  do n=0,nd-1
    fpSum(n)=0.
    pevttSum(n)=0.
    pevttxSum(n)=0.
    pevttySum(n)=0.
    pevttLSum(n)=0.
    pevttttSum(n)=0.

    petttSum=0.

    call ogderiv(ep, 0,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,ex+n, es(n)   ) 
    call ogderiv(ep, 1,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,ex+n, est(n)  )
    call ogderiv(ep, 2,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,ex+n, estt(n) )

    call ogderiv(ep, 0,1,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,ex+n, esx(n) )
    call ogderiv(ep, 0,0,1,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,ex+n, esy(n) )

    call ogderiv(ep, 0,2,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,ex+n, esxx(n) )
    call ogderiv(ep, 0,0,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,ex+n, esyy(n) )

    call ogderiv(ep, 1,1,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,ex+n, estx(n) )
    call ogderiv(ep, 1,0,1,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,ex+n, esty(n) )

    call ogderiv(ep, 2,1,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,ex+n, esttx(n) )
    call ogderiv(ep, 2,0,1,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,ex+n, estty(n) )

    call ogderiv(ep, 0,3,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,ex+n, esxxx(n) )
    call ogderiv(ep, 0,2,1,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,ex+n, esxxy(n) )
    call ogderiv(ep, 0,1,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,ex+n, esxyy(n) )
    call ogderiv(ep, 0,0,3,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,ex+n, esyyy(n) )

    call ogderiv(ep, 3,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,ex+n, esttt(n) )
    call ogderiv(ep, 4,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,ex+n, estttt(n) )

    call ogderiv(ep, 1,2,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,ex+n, estxx(n) )
    call ogderiv(ep, 1,0,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,ex+n, estyy(n) )


    call ogderiv(ep, 2,2,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,ex+n, esttxx(n) )
    call ogderiv(ep, 2,0,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,ex+n, esttyy(n) )

    call ogderiv(ep, 0,4,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,ex+n, esxxxx(n) )
    call ogderiv(ep, 0,2,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,ex+n, esxxyy(n) )
    call ogderiv(ep, 0,0,4,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,ex+n, esyyyy(n) )

    ! L = c^2*Delta
    esL  = c**2*( esxx(n)   + esyy(n) )
    estL = c**2*( estxx(n)  + estyy(n) )
    esttL= c**2*( esttxx(n) + esttyy(n) )

    esLx  = c**2*( esxxx(n)   + esxyy(n) )
    esLy  = c**2*( esxxy(n)   + esyyy(n) )

    ! L^2 : 
    esLL = c**4*( esxxxx(n)  + 2.*esxxyy(n) + esyyyy(n) )

    do jv=0,numberOfPolarizationVectors-1
      ! The TZ component is offset by pxc
      pc = pxc + jv*nd
      call ogderiv(ep, 0,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,pc+n, pe(n)   )
      call ogderiv(ep, 1,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,pc+n, pet(n)  )
      call ogderiv(ep, 2,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,pc+n, pett(n) )
      call ogderiv(ep, 3,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,pc+n, pettt(n) )
      call ogderiv(ep, 4,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,pc+n, petttt(n) )

      call ogderiv(ep, 0,1,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,pc+n, pex(n) )
      call ogderiv(ep, 0,0,1,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,pc+n, pey(n) )

      call ogderiv(ep, 0,2,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,pc+n, pexx(n) )
      call ogderiv(ep, 0,0,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,pc+n, peyy(n) )

      call ogderiv(ep, 1,1,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,pc+n, petx(n) )
      call ogderiv(ep, 1,0,1,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,pc+n, pety(n) )

      call ogderiv(ep, 1,2,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,pc+n, petxx(n) )
      call ogderiv(ep, 1,0,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,pc+n, petyy(n) )

      call ogderiv(ep, 2,1,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,pc+n, pettx(n) )
      call ogderiv(ep, 2,0,1,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,pc+n, petty(n) )

      call ogderiv(ep, 2,2,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,pc+n, pettxx(n) )
      call ogderiv(ep, 2,0,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,pc+n, pettyy(n) )

      peL  = c**2*( pexx(n)   + peyy(n) )
      petL = c**2*( petxx(n)  + petyy(n) )
      pettL= c**2*( pettxx(n) + pettyy(n) )

      ! Normal TZ forcing for P_{n,jv} equation: 
      fpv(n,jv) = pett(n)   + b1v(jv)*pet(n)   + b0v(jv)*pe(n)   - a0v(jv)*es(n)   - a1v(jv)*est(n)
      fPt(n,jv) = pettt(n)  + b1v(jv)*pett(n)  + b0v(jv)*pet(n)  - a0v(jv)*est(n)  - a1v(jv)*estt(n)
      fPtt(n,jv)= petttt(n) + b1v(jv)*pettt(n) + b0v(jv)*pett(n) - a0v(jv)*estt(n) - a1v(jv)*esttt(n)
      LfP(n,jv) = pettL     + b1v(jv)*petL     + b0v(jv)*peL     - a0v(jv)*esL     - a1v(jv)*estL

      fpvx(n,jv)= pettx(n)  + b1v(jv)*petx(n)  + b0v(jv)*pex(n)  - a0v(jv)*esx(n)  - a1v(jv)*estx(n)
      fpvy(n,jv)= petty(n)  + b1v(jv)*pety(n)  + b0v(jv)*pey(n)  - a0v(jv)*esy(n)  - a1v(jv)*esty(n)

      ! write(*,'(" n=",i4," LfP=",e10.4," pettL,petL,peL,esL,estL=",5e12.4)') n,LfP(n,jv),pettL,petL,peL,esL,estL
      ! write(*,'(" pe,pet,pett,pettt,petttt=",5e12.4)') pe(n),pet(n),pett(n),pettt(n),petttt(n)

      ! write(*,'("TZ: n,jv=",2i4," pex,pey,pexx,peyy=",4(1pe12.4))') n,jv,pex(n),pey(n),pexx(n),peyy(n)

      ! Save ptt for checking later
      pevtt(n,jv)=pett(n)
      pevttx(n,jv)=pettx(n)
      pevtty(n,jv)=petty(n)
      pevtttt(n,jv)=petttt(n)
      pevttLSum(n) = pevttLSum(n)  + pettL
      pevttttSum(n)= pevttttSum(n) + petttt(n) 

      ! Keep some sums: 
      fpSum(n)   = fpSum(n)  + fpv(n,jv)
      pevttSum(n)  = pevttSum(n)  + pett(n) 
      pevttxSum(n) = pevttxSum(n) + pettx(n)
      pevttySum(n) = pevttySum(n) + petty(n)

      petttSum  = petttSum  + pettt(n) 


    end do 

    ! TZ forcing for E_{n} equation:
    ! E_tt - c^2 Delta E + alphaP*Ptt  = 
    fev(n) = estt(n)   - esL   + alphaP*pevttSum(n)
    fEt(n) = esttt(n)  - estL  + alphaP*petttSum
    fEtt(n)= estttt(n) - esttL + alphaP*pevttttSum(n)

    fevx(n) = esttx(n) - esLx   + alphaP*pevttxSum(n)
    fevy(n) = estty(n) - esLy   + alphaP*pevttySum(n)
    

    ! write(*,'("--> fEtt=",e10.2," estttt,esttL,pettttSum=",3e10.2)')  fEtt(n),estttt(n),esttL,pettttSum
    LfE(n) = esttL     - esLL  + alphaP*pevttLSum(n)
    
 end do
end if

#endMacro 

!-------------------------------------------------------------------------------------------
! Macro: Evaluate TZ forcing for dispersive equations in 2D 
!
! Output
!    fpv1(n,jv) : RHS To Pv_{n,jv} equation on domain 1
!    fpv2(n,jv) : RHS To Pv_{n,jv} equation on domain 2
!    fev1(n)    : RHS to E_{n} equation on domain 1
!    fev2(n)    : RHS to E_{n} equation on domain 2
!-------------------------------------------------------------------------------------------
#beginMacro getDispersiveTZForcingOrder4(fpv1,fpv2,fev1,fev2)

  if( twilightZone.eq.1 )then
    evalTZForcingGDMOrder4(xy1,i1,i2,i3,dispersionModel1,numberOfPolarizationVectors1,c1,alphaP1,a0v1,a1v1,b0v1,b1v1,fpv1,fpSum1,fev1,\
                          LfE1,fEt1,fEtt1,LfP1,fPt1,fPtt1,pevtt1,pevttx1,pevtty1,pevtttt1,fevx1,fevy1,fpvx1,fpvy1,pevttSum1,pevttxSum1,pevttySum1,pevttLSum1,pevttttSum1)
    evalTZForcingGDMOrder4(xy2,j1,j2,j3,dispersionModel2,numberOfPolarizationVectors2,c2,alphaP2,a0v2,a1v2,b0v2,b1v2,fpv2,fpSum2,fev2,\
                          LfE2,fEt2,fEtt2,LfP2,fPt2,fPtt2,pevtt2,pevttx2,pevtty2,pevtttt2,fevx2,fevy2,fpvx2,fpvy2,pevttSum2,pevttxSum2,pevttySum2,pevttLSum2,pevttttSum2)
  end if

#endMacro 

! --------------------------------------------------------------------------
! Macro: Evaluate the GDM jump conditions in 2D, order=4
! --------------------------------------------------------------------------
#beginMacro eval2dJumpDispersiveOrder4()
 f(0)=(u1x+v1y) - \
      (u2x+v2y)

 ! [ n.Delta( E )/mu ]=0
 f(1)=(an1*u1Lap+an2*v1Lap) - \
      (an1*u2Lap+an2*v2Lap)

 ! [ nv X curl( E) /mu ] = 0 
 f(2)=(v1x-u1y) - \
      (v2x-u2y)

! f(3)=(tau1*u1Lap+tau2*v1Lap)/eps1 - \
!      (tau1*u2Lap+tau2*v2Lap)/eps2

 f(3)=( ( tau1*u1Lap +tau2*v1Lap )/epsmu1 - alphaP1*(tau1*fp1(0)+tau2*fp1(1)) ) - \
      ( ( tau1*u2Lap +tau2*v2Lap )/epsmu2 - alphaP2*(tau1*fp2(0)+tau2*fp2(1)) )

 ! [ Delta( div(E) ) ] = 0 
 f(4)=(u1xxx+u1xyy+v1xxy+v1yyy) - \
      (u2xxx+u2xyy+v2xxy+v2yyy)

 ! Dispersive:
 !   [ (c^2/mu)*{(Delta v).x - (Delta u).y} - (alphaP/mu)*( Py.ttx - Px.tty) ] =0 
 !    fPttx = P.ttx 
 !    fPtty = P.tty 
 f(5)= ( ((v1xxx+v1xyy)-(u1xxy+u1yyy))/epsmu1 - alphaP1*(fPttx1(1) - fPtty1(0)) )/mu1 \
      -( ((v2xxx+v2xyy)-(u2xxy+u2yyy))/epsmu2 - alphaP2*(fPttx2(1) - fPtty2(0)) )/mu2

 if( setDivergenceAtInterfaces.eq.0 )then
  !  [ nv.( c^2*Delta^2(E) - alphaP*Delta(Ptt) )/mu ] = 0 
  ! Note: fLptt = c^2*Delta( Ptt ) 
  f(6)= ( (an1*u1LapSq+an2*v1LapSq)/epsmu1  -alphaP1*( an1*fLPtt1(0)+an2*fLPtt1(1) )*epsmu1  )/mu1 \
       -( (an1*u2LapSq+an2*v2LapSq)/epsmu2  -alphaP2*( an1*fLPtt2(0)+an2*fLPtt2(1) )*epsmu2  )/mu2 
 else
  f(6)=(u1x+v1y)
 end if

 ! [ tv.( c^4*Delta^2(E) - alphaP*c^2*Delta(P.tt) - alphaP*P.tttt) ]=0 
 f(7)=(tau1*u1LapSq+tau2*v1LapSq)/epsmu1**2 - \
      (tau1*u2LapSq+tau2*v2LapSq)/epsmu2**2 \
      -alphaP1*( ( tau1*fLPtt1(0) + tau2*fLPtt1(1) ) + tau1*fPtttt1(0) + tau2*fPtttt1(1) ) \
      +alphaP2*( ( tau1*fLPtt2(0) + tau2*fLPtt2(1) ) + tau1*fPtttt2(0) + tau2*fPtttt2(1) )
     
 if( twilightZone.eq.1 )then
   call ogderiv(ep, 0,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, uexx )
   call ogderiv(ep, 0,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, ueyy )
   call ogderiv(ep, 0,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, vexx )
   call ogderiv(ep, 0,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, veyy )

   call ogderiv(ep, 0,2,1,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, uexxy )
   call ogderiv(ep, 0,0,3,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, ueyyy )

   call ogderiv(ep, 0,3,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, vexxx )
   call ogderiv(ep, 0,1,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, vexyy )

   call ogderiv(ep, 0,4,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, uexxxx )
   call ogderiv(ep, 0,2,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, uexxyy )
   call ogderiv(ep, 0,0,4,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, ueyyyy )

   call ogderiv(ep, 0,4,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, vexxxx )
   call ogderiv(ep, 0,2,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, vexxyy )
   call ogderiv(ep, 0,0,4,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, veyyyy )

   ueLap = uexx + ueyy
   veLap = vexx + veyy
   ueLapSq = uexxxx + 2.*uexxyy + ueyyyy
   veLapSq = vexxxx + 2.*vexxyy + veyyyy

   f(3) = f(3) - ( tau1*ueLap +tau2*veLap )*(1./epsmu1-1./epsmu2) \
               + alphaP1*(tau1*pevttSum1(0)+tau2*pevttSum1(1)) \
               - alphaP2*(tau1*pevttSum2(0)+tau2*pevttSum2(1))

   f(5) = f(5) - ((vexxx+vexyy)-(uexxy+ueyyy))*(1./epsmu1-1./epsmu2) \
               + (alphaP1/mu1)*( pevttxSum1(1) - pevttySum1(0) ) \
               - (alphaP2/mu2)*( pevttxSum2(1) - pevttySum2(0) )


   if( setDivergenceAtInterfaces.eq.0 )then

     f(6) = f(6) - (an1*ueLapSq+an2*veLapSq)*(1./(epsmu1*mu1) - 1./(epsmu2*mu2) ) \
                 + (alphaP1/mu1)*( an1*pevttLSum1(0) + an2*pevttLSum1(1) )*epsmu1 \
                 - (alphaP2/mu2)*( an1*pevttLSum2(0) + an2*pevttLSum2(1) )*epsmu2 

   end if

   f(7) = f(7) - (tau1*ueLapSq+tau2*veLapSq)*(1./epsmu1**2 - 1./epsmu2**2) \
               + alphaP1*( ( tau1*pevttLSum1(0) + tau2*pevttLSum1(1) ) + tau1*pevttttSum1(0) + tau2*pevttttSum1(1) ) \
               - alphaP2*( ( tau1*pevttLSum2(0) + tau2*pevttLSum2(1) ) + tau1*pevttttSum2(0) + tau2*pevttttSum2(1) ) 


 end if
#endMacro

! ========================================================================
! Macro: Getting forcing for GDM
!  Input:
!    ec : E component
!    pc : P component
!  Output:
!    fe : forcing for E ( includes factor of dt^2)
!    fpv(iv) : forcing for polarization vector iv=0,1,2,...( includes factor of dt^2)
! ========================================================================

#beginMacro getGDMForcing(fe,fpv, ec,pc, t, i1,i2,i3,xy,cSq,numberOfPolarizationVectors,alphaP,a0v,a1v,b0v,b1v)

 ! --- compute forcing terms for ADE-GDM equations ----
 if( addForcing.ne.0 )then

   if( forcingOption.eq.twilightZoneForcing )then

     if( nd.eq.2 )then
       call ogderiv(ep, 0,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, ec, ue )
       call ogderiv(ep, 1,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, ec, uet )
       call ogderiv(ep, 2,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, ec, uett )

       call ogderiv(ep, 0,2,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, ec, uexx )
       call ogderiv(ep, 0,0,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, ec, ueyy )

       fe = dtSq*(uett - cSq*(uexx + ueyy))

       ! write(*,'("fe,uett,dtSq=",3(1pe14.4)," c**2*Delta*E=",1pe14.4," exact=",1pe14.2)') fe,uett,dtSq,cSq*(uexx + ueyy),LE(m)
     else
       stop 33387
     end if


     do jv=0,numberOfPolarizationVectors-1
       ! The TZ component is offset by pxc
       pce = pc+jv*nd
       if( nd.eq.2 )then
         call ogderiv(ep, 0,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,pce, p0   )
         call ogderiv(ep, 1,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,pce, p0t  )
         call ogderiv(ep, 2,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,pce, p0tt )
       else
         stop 1111
       end if

       fe = fe + dtSq*alphaP*p0tt
       ! write(*,'(" fe,p0tt=",2e12.4)') fe,p0tt
       fpv(jv) = dtSq*( p0tt + b1v(jv)*p0t + b0v(jv)*p0 - a0v(jv)*ue - a1v(jv)*uet )

     end do

   else

     fe=0.
     do jv=0,numberOfPolarizationVectors-1
       fpv(jv)=0.
     end do

   end if
 end if
#endMacro



! ===========================================================================================
! Macro:     DISPERSIVE: CURVILINEAR, 2D, ORDER 2
!    Evaluate P to second-order
! Input:
!    i1,i2,i3 : evaluate P at this point
!    u,um : u(t-dt), u(t-2*dt)
! Output:
!    evals,pv
! ===========================================================================================
#beginMacro evalPolarization2dOrder2(evals,pvals, t, i1,i2,i3,u,um,p,pm, xy,rsxy,aj,cSq,numberOfPolarizationVectors,alphaP,a0v,a1v,b0v,b1v)

 ! -- first compute some coefficients ---
 beta=0.
 do jv=0,numberOfPolarizationVectors-1
   betav(jv) = 1./( 1.+.5*dt*b1v(jv) )
   beta = beta + .5*dt*a1v(jv)*betav(jv)
   fpv(jv)=0.  ! initialize if not used
 end do

  ! Evaluate u.xx, u.yy, ...
  #perl $ORDER=2;
  opEvalJacobianDerivatives(rsxy,i1,i2,i3,aj,1)
  evalSecondDerivs(rsxy,aj,u,i1,i2,i3,ex,uu1,u1n)
  evalSecondDerivs(rsxy,aj,u,i1,i2,i3,ey,vv1,v1n)
  ! Here are c^2*Delta(E) 
  LE(0) = cSq*u1nLap
  LE(1) = cSq*v1nLap

  ! -- loop over components of the vector --
 do m=0,nd-1
   pc=pxc+m ! FIX ME 
   ec=ex+m

   ! Compute fe, fpv(iv) :
   getGDMForcing(fe,fpv, ec,pc, t, i1,i2,i3,xy,cSq,numberOfPolarizationVectors,alphaP,a0v,a1v,b0v,b1v)


   ev = u(i1,i2,i3,ec)
   evm=um(i1,i2,i3,ec)

   rhsP = 0.
   pSum=0.
   do jv=0,numberOfPolarizationVectors-1
     pvc(jv)= p(i1,i2,i3,m+jv*nd)
     pvm(jv)=pm(i1,i2,i3,m+jv*nd) 

     rhspv(jv) = 2.*pvc(jv)-pvm(jv) + .5*dt*( b1v(jv)*pvm(jv) -a1v(jv)*evm ) + dtSq*( -b0v(jv)*pvc(jv) + a0v(jv)*ev ) + fpv(jv)
     rhsP = rhsP + betav(jv)*rhspv(jv)
     pSum = pSum + 2.*pvc(jv) - pvm(jv)
   end do

   rhsE = 2.*u(i1,i2,i3,ec)-um(i1,i2,i3,ec)+ dtSq*LE(m)+ fe + alphaP*( pSum - rhsP ) 

   evn = rhsE / (1.+ alphaP*beta)
   evals(m)=evn
   ! write(*,'(" ec=",i4," E(t-dt),E(t-2*dt),fe,rhsE=",4(1pe14.4))') ec,u(i1,i2,i3,ec),um(i1,i2,i3,ec),fe,rhsE

   ! un(i1,i2,i3,ec) = evn
   do jv=0,numberOfPolarizationVectors-1
     ! pn(i1,i2,i3,m+jv*nd)  = betav(jv)*( .5*dt*a1v(jv)*evn + rhspv(jv) )
     pvals(m+jv*nd)  = betav(jv)*( .5*dt*a1v(jv)*evn + rhspv(jv) )
  end do


 end do
#endMacro


#beginMacro evalDerivativesForDispersive2dOrder4()
 ! Store c^2*Delta(E) in a vector 
 LE1(0)=(c1**2)*u1Lap
 LE1(1)=(c1**2)*v1Lap
 
 LE2(0)=(c2**2)*u2Lap
 LE2(1)=(c2**2)*v2Lap
 
 ! Store L^2(E) 
 LLE1(0)=(c1**4)*u1LapSq
 LLE1(1)=(c1**4)*v1LapSq
 
 LLE2(0)=(c2**4)*u2LapSq
 LLE2(1)=(c2**4)*v2LapSq
 
 ! Store (LE).x an (LE).y 
 LEx1(0) = (c1**2)*( u1xxx + u1xyy )
 LEx1(1) = (c1**2)*( v1xxx + v1xyy )

 LEy1(0) = (c1**2)*( u1xxy + u1yyy )
 LEy1(1) = (c1**2)*( v1xxy + v1yyy )

 LEx2(0) = (c2**2)*( u2xxx + u2xyy )
 LEx2(1) = (c2**2)*( v2xxx + v2xyy )

 LEy2(0) = (c2**2)*( u2xxy + u2yyy )
 LEy2(1) = (c2**2)*( v2xxy + v2yyy )

 ! We also need derivatives at the old time:
 ! These next derivatives may only be needed to order2, but use order 4 for now so exact for degree 4
 #perl $ORDER=4;
 ! perl $ORDER=2;
 opEvalJacobianDerivatives(rsxy1,i1,i2,i3,aj1,1)
 evalSecondDerivs(rsxy1,aj1,u1n,i1,i2,i3,ex,uu1,u1n)
 evalSecondDerivs(rsxy1,aj1,u1n,i1,i2,i3,ey,vv1,v1n)
 ! Here are c^2*Delta(E) at the old time: 
 LE1m(0) = (c1**2)*u1nLap
 LE1m(1) = (c1**2)*v1nLap

 opEvalJacobianDerivatives(rsxy2,j1,j2,j3,aj2,1)
 evalSecondDerivs(rsxy2,aj2,u2n,j1,j2,j3,ex,uu2,u2n)
 evalSecondDerivs(rsxy2,aj2,u2n,j1,j2,j3,ey,vv2,v2n)
 LE2m(0) = (c2**2)*u2nLap
 LE2m(1) = (c2**2)*v2nLap

 evx1(0) = u1x
 evx1(1) = v1x
 evy1(0) = u1y
 evy1(1) = v1y 

 evnx1(0) = u1nx
 evnx1(1) = v1nx
 evny1(0) = u1ny
 evny1(1) = v1ny 

 evx2(0) = u2x
 evx2(1) = v2x
 evy2(0) = u2y
 evy2(1) = v2y 

 evnx2(0) = u2nx
 evnx2(1) = v2nx
 evny2(0) = u2ny
 evny2(1) = v2ny 
#endMacro


! =============================================================================================
!   Evaluate the jump conditions for the GDM interface equations
! =============================================================================================
#beginMacro evaluateDispersiveInterfaceEquations2dOrder4()

 ! Evaluate TZ forcing for dispersive equations in 2D 
 getDispersiveTZForcingOrder4(fpv1,fpv2,fev1,fev2)

 evalDerivs2dOrder4()

 evalDerivativesForDispersive2dOrder4()

 ! eval dispersive forcings for domain 1
 getDispersiveForcingOrder4(LEFT,i1,i2,i3, fp1, fpv1,fev1,p1,p1n,p1m, u1,u1n,u1m, dispersionModel1,\
    numberOfPolarizationVectors1,alphaP1,c2PttEsum1,c2PttLEsum1,c4PttLEsum1,c4PttLLEsum1,c2PttttLEsum1,c2PttttLLEsum1,\
    a0v1,a1v1,b0v1,b1v1,LE1,LLE1,LE1m,LfE1,LfP1,fEt1,fEtt1,fPt1,fPtt1,pevtt1,pevttx1,pevtty1,pevtttt1,\
    evx1,evy1,evnx1,evny1,fevx1,fevy1,fpvx1,fpvy1,LEx1,LEy1,fPttx1,fPtty1,fLPtt1,fPtttt1) 

 ! eval dispersive forcings for domain 2
 getDispersiveForcingOrder4(RIGHT,j1,j2,j3, fp2, fpv2,fev2,p2,p2n,p2m, u2,u2n,u2m, dispersionModel2,\
    numberOfPolarizationVectors2,alphaP2,c2PttEsum2,c2PttLEsum2,c4PttLEsum2,c4PttLLEsum2,c2PttttLEsum2,c2PttttLLEsum2,\
    a0v2,a1v2,b0v2,b1v2,LE2,LLE2,LE2m,LfE2,LfP2,fEt2,fEtt2,fPt2,fPtt2,pevtt2,pevttx2,pevtty2,pevtttt2,\
    evx2,evy2,evnx2,evny2,fevx2,fevy2,fpvx2,fpvy2,LEx2,LEy2,fPttx2,fPtty2,fLPtt2,fPtttt2) 


 ! first evaluate the equations we want to solve with the wrong values at the ghost points: (assigns f(0:7))
 eval2dJumpDispersiveOrder4()


 if( debug.gt.7 ) write(debugFile,'(" --> 4cth: j1,j2=",2i4," u1xx,u1yy,u2xx,u2yy=",4e10.2)') j1,j2,u1xx,\
  u1yy,u2xx,u2yy
   ! '
  if( debug.gt.3 ) write(debugFile,'(" --> 4cth: i1,i2=",2i4," f(start)=",8e10.2)') i1,i2,f(0),f(1),f(2),f(3),f(4),f(5),f(6),f(7)

#endMacro





! ===========================================================================================
!  Assign P in ghost points for the fourth-order method
!
! *** THIS IS NOT USED CURRENTLY -- was created to fix a bug that was caused by a wrong alphaP
!
! ===========================================================================================
#beginMacro assignPolarizationOnGhostOrder4()

 if( knownSolutionOption.eq.userDefinedKnownSolution )then
   write(*,'(" assignPolarizationOnGhostOrder4: set P in ghost 1 with 2nd-order")') 
 end if

 if( nd.eq.2 )then
  ! --- LOOP over the interface ---
  beginGhostLoopsMask2d()
   if( dispersionModel1.ne.noDispersion )then

     ! Compute P at the first ghost using the 2nd order method
     i1m=i1-is1
     i2m=i2-is2
     i3m=i3
     tm=t-dt
     evalPolarization2dOrder2(evals,pvalsm, tm, i1m,i2m,i3m,u1n,u1m,p1n,p1m, xy1,rsxy1,aj1,cSq1,\
                              numberOfPolarizationVectors1,alphaP1,a0v1,a1v1,b0v1,b1v1)

     evalPolarization2dOrder2(evals,pvals , tm, i1,i2,i3,u1n,u1m,p1n,p1m, xy1,rsxy1,aj1,cSq1,\
                              numberOfPolarizationVectors1,alphaP1,a0v1,a1v1,b0v1,b1v1)

     i1p=i1+is1
     i2p=i2+is2
     i3p=i3
     evalPolarization2dOrder2(evals,pvalsp, tm, i1p,i2p,i3p,u1n,u1m,p1n,p1m, xy1,rsxy1,aj1,cSq1,\
                              numberOfPolarizationVectors1,alphaP1,a0v1,a1v1,b0v1,b1v1)


     if( .false. .and. twilightZone.eq.1 )then
       do jv=0,numberOfPolarizationVectors1-1
         do n=0,nd-1
           pc = n + jv*nd
           ! The TZ component is offset by pxc
           pce = pxc + pc
           call ogderiv(ep, 0,0,0,0, xy1(i1m,i2m,i3m,0),xy1(i1m,i2m,i3m,1),0.,t,pce, pvalse(pc)   )
         end do
       end do

!         call ogderiv(ep, 0,0,0,0, xy1(i1m,i2m,i3m,0),xy1(i1m,i2m,i3m,1),0.,tm,ex, ue   )
! write(*,'(" ghost: Ex(t-dt): u1n,exact=",2(1pe14.4)," error=",1pe12.2)') u1n(i1m,i2m,i3m,ex),ue,u1n(i1m,i2m,i3m,ex)-ue

!         call ogderiv(ep, 0,0,0,0, xy1(i1m,i2m,i3m,0),xy1(i1m,i2m,i3m,1),0.,tm-dt,ex, ue   )
! write(*,'(" ghost: Ex(t-2*dt): u1m,exact=",2(1pe14.4)," error=",1pe12.2)') u1m(i1m,i2m,i3m,ex),ue,u1m(i1m,i2m,i3m,ex)-ue

         call ogderiv(ep, 0,0,0,0, xy1(i1m,i2m,i3m,0),xy1(i1m,i2m,i3m,1),0.,t,ex, ue   )
 write(*,'(" ghost: Ex(t): evn,exact=",2(1pe14.4)," error=",1pe12.2)') evals(0),ue,evals(0)-ue


!         call ogderiv(ep, 0,0,0,0, xy1(i1m,i2m,i3m,0),xy1(i1m,i2m,i3m,1),0.,tm,pxc, p0   )
! write(*,'(" ghost: p1n,exact=",2(1pe14.4)," error=",1pe12.2)') p1n(i1m,i2m,i3m,0),p0,p1n(i1m,i2m,i3m,0)-p0

!         call ogderiv(ep, 0,0,0,0, xy1(i1m,i2m,i3m,0),xy1(i1m,i2m,i3m,1),0.,tm-dt,pxc, p0   )
! write(*,'(" ghost: p1m,exact=",2(1pe14.4)," error=",1pe12.2)') p1m(i1m,i2m,i3m,0),p0,p1m(i1m,i2m,i3m,0)-p0


 write(*,'(" grid1,i1,i2=",3i3," p(-1)=",2(1pe14.4)," exact=",2(1pe14.4)," err=",2(1pe12.2))') \
            grid1,i1,i2,pvals(0),pvals(1),pvalse(0),pvalse(1),pvals(0)-pvalse(0),pvals(1)-pvalse(1)


     end if
     if( .true. .and. knownSolutionOption.eq.userDefinedKnownSolution )then
       ! Evaluate the user defined known solution
       numberOfTimeDerivatives=0
       call evalUserDefinedKnownSolution( t, grid1,i1m,i2m,i3m,evalse,pvalse, numberOfTimeDerivatives )

 !write(*,'(" grid1,i1,i2=",3i3," p(-1)=",2(1pe14.4)," exact=",2(1pe14.4)," err=",2(1pe12.2))') \
 !           grid1,i1,i2,pvals(0),pvals(1),pvalse(0),pvalse(1),pvals(0)-pvalse(0),pvals(1)-pvalse(1)

       ! **** SET EXACT VALUE ON FIRST GHOST ****
       if( .false. )then
        do jv=0,numberOfPolarizationVectors1-1
          do n=0,nd-1
            pc = n + jv*nd
            pvals(pc)=pvalse(pc)
          end do
        end do
       end if
     end if

     do jv=0,numberOfPolarizationVectors1-1
       do n=0,nd-1
         pc = n + jv*nd
         p1(i1m,i2m,i3m,pc) = pvals(pc)
         ! Second-order estimate for "second-derivative": 
         !** dpdm = pvalsm(pc)-2.*pvals(pc)+pvalsp(pc)
         ! Combine 2nd-order values and fourth-order interior+bndry values
         !** p1(i1m     ,i2m     ,i3m,pc) =2.*p1(i1,i2,i3,pc)-p1(i1p     ,i2p     ,i3p,pc) + dpdm
         ! 2nd -ghost : use D+D- with 2h --> leads to dpdm*4.
         ! p1(i1-2*is1,i2-2*is2,i3 ,pc) =2.*p1(i1,i2,i3,pc)-p1(i1+2*is1,i2+2*is2,i3 ,pc) + dpdm*.4
       end do
     end do
   end if
   if( dispersionModel2.ne.noDispersion )then

     ! Compute P at the first ghost using the 2nd order method
     j1m=j1-js1
     j2m=j2-js2
     j3m=j3
     tm=t-dt
     evalPolarization2dOrder2(evals,pvalsm, tm, j1m,j2m,j3m,u2n,u2m,p2n,p2m, xy2,rsxy2,aj2,cSq2,\
                              numberOfPolarizationVectors2,alphaP2,a0v2,a1v2,b0v2,b1v2)

     evalPolarization2dOrder2(evals,pvals , tm, j1 ,j2 ,j3 ,u2n,u2m,p2n,p2m, xy2,rsxy2,aj2,cSq2,\
                              numberOfPolarizationVectors2,alphaP2,a0v2,a1v2,b0v2,b1v2)

     j1p=j1+js1
     j2p=j2+js2
     j3p=j3
     evalPolarization2dOrder2(evals,pvalsp, tm, j1p,j2p,j3p,u2n,u2m,p2n,p2m, xy2,rsxy2,aj2,cSq2,\
                              numberOfPolarizationVectors2,alphaP2,a0v2,a1v2,b0v2,b1v2)

     if( .false. .and. twilightZone.eq.1 )then
       do jv=0,numberOfPolarizationVectors2-1
         do n=0,nd-1
           pc = n + jv*nd
           ! The TZ component is offset by pxc
           pce = pxc + pc
           call ogderiv(ep, 0,0,0,0, xy2(j1m,j2m,j3m,0),xy2(j1m,j2m,j3m,1),0.,t,pce, pvalse(pc)   )
         end do
       end do

 write(*,'(" grid2,j1,j2=",3i3," p(-1)=",2(1pe14.4)," exact=",2(1pe14.4)," err=",2(1pe12.2))') \
            grid2,j1,j2,pvals(0),pvals(1),pvalse(0),pvalse(1),pvals(0)-pvalse(0),pvals(1)-pvalse(1)

     end if
     if( .true. .and. knownSolutionOption.eq.userDefinedKnownSolution )then
       ! Evaluate the user defined known solution
       numberOfTimeDerivatives=0
       call evalUserDefinedKnownSolution( t, grid2,j1m,j2m,j3,evalse,pvalse, numberOfTimeDerivatives )

! write(*,'(" grid2,j1,j2=",3i3," p(-1)=",2(1pe14.4)," exact=",2(1pe14.4)," err=",2(1pe12.2))') \
!            grid2,j1,j2,pvals(0),pvals(1),pvalse(0),pvalse(1),pvals(0)-pvalse(0),pvals(1)-pvalse(1)

       ! **** SET EXACT VALUE ON FIRST GHOST ****
       if( .false. )then
        do jv=0,numberOfPolarizationVectors2-1
          do n=0,nd-1
            pc = n + jv*nd
            pvals(pc)=pvalse(pc)
          end do
        end do
       end if
     end if

     do jv=0,numberOfPolarizationVectors2-1
       do n=0,nd-1
         pc = n + jv*nd
         p2(j1m,j2m,j3m,pc) = pvals(pc)
         ! dpdm = pvalsm(pc)-2.*pvals(pc)+pvalsp(pc)
         ! Combine 2nd-order values and fourth-order interior+bndry values
         ! p2(j1m     ,j2m     ,j3m,pc) =2.*p2(j1,j2,j3,pc)-p2(j1p     ,j2p     ,j3p,pc) + dpdm
         ! p2(j1-2*js1,j2-2*js2,j3 ,pc) =2.*p2(j1,j2,j3,pc)-p2(j1+2*js1,j2+2*js2,j3 ,pc) + dpdm*4.
       end do
     end do
   end if
  endLoopsMask2d()

 else
   write(*,'("cgmx:interface3d:assignPolarizationOnGhostOrder4: finish me 3D")') 
   stop 2975
 end if

 ! stop 9876

#endMacro

! ----------------------------------------------------------------------------------
!  Macro:
!    Evaluate the interface equations for checking the coefficients
! ----------------------------------------------------------------------------------
#beginMacro evalDispersiveInterfaceEquations24c()

 ! Evaluate TZ forcing for dispersive equations in 2D 
 getDispersiveTZForcingOrder4(fpv1,fpv2,fev1,fev2)

 evalDerivs2dOrder4()

 evalDerivativesForDispersive2dOrder4()

 ! eval dispersive forcings for domain 1
 getDispersiveForcingOrder4(LEFT,i1,i2,i3, fp1, fpv1,fev1,p1,p1n,p1m, u1,u1n,u1m, dispersionModel1,\
    numberOfPolarizationVectors1,alphaP1,c2PttEsum1,c2PttLEsum1,c4PttLEsum1,c4PttLLEsum1,c2PttttLEsum1,c2PttttLLEsum1,\
    a0v1,a1v1,b0v1,b1v1,LE1,LLE1,LE1m,LfE1,LfP1,fEt1,fEtt1,fPt1,fPtt1,pevtt1,pevttx1,pevtty1,pevtttt1,\
    evx1,evy1,evnx1,evny1,fevx1,fevy1,fpvx1,fpvy1,LEx1,LEy1,fPttx1,fPtty1,fLPtt1,fPtttt1) 

 ! eval dispersive forcings for domain 2
 getDispersiveForcingOrder4(RIGHT,j1,j2,j3, fp2, fpv2,fev2,p2,p2n,p2m, u2,u2n,u2m, dispersionModel2,\
    numberOfPolarizationVectors2,alphaP2,c2PttEsum2,c2PttLEsum2,c4PttLEsum2,c4PttLLEsum2,c2PttttLEsum2,c2PttttLLEsum2,\
    a0v2,a1v2,b0v2,b1v2,LE2,LLE2,LE2m,LfE2,LfP2,fEt2,fEtt2,fPt2,fPtt2,pevtt2,pevttx2,pevtty2,pevtttt2,\
    evx2,evy2,evnx2,evny2,fevx2,fevy2,fpvx2,fpvy2,LEx2,LEy2,fPttx2,fPtty2,fLPtt2,fPtttt2) 

 ! first evaluate the equations we want to solve with the wrong values at the ghost points: (assigns f(0:7))
 eval2dJumpDispersiveOrder4()

#endMacro


! --------------------------------------------------------------------------
! Macro: Assign interface ghost values, DIM=2, ORDER=4, GRID=Curvilinear
!         DISPERSIVE CASE -- GDM 
! --------------------------------------------------------------------------
#beginMacro assignDispersiveInterfaceGhost24c()

 ! ****************************************************************
 ! ***********  DISPERSIVE, 2D, ORDER=4, CURVILINEAR **************
 ! ****************************************************************

 INFO("24c-GDM")

! --- initialize some forcing functions ---
do n=0,nd-1
  fev1(n)=0.
  LfE1(n)=0.
  fEt1(n)=0.
  fEtt1(n)=0.

  fev2(n)=0.
  LfE2(n)=0.
  fEt2(n)=0.
  fEtt2(n)=0.

  fevx1(n)=0.
  fevy1(n)=0.
  fevx2(n)=0.
  fevy2(n)=0.

  do jv=0,numberOfPolarizationVectors1-1
    fpv1(n,jv)=0.
    LfP1(n,jv)=0.
    fPt1(n,jv)=0.
    fPtt1(n,jv)=0.

    fpvx1(n,jv)=0.
    fpvy1(n,jv)=0.
  end do
  do jv=0,numberOfPolarizationVectors2-1
    fpv2(n,jv)=0.
    LfP2(n,jv)=0.
    fPt2(n,jv)=0.
    fPtt2(n,jv)=0.

    fpvx2(n,jv)=0.
    fpvy2(n,jv)=0.
  end do
end do

 ! write(*,'("p1=",(15(e10.2,1x)))') (((p1(i1,i2,i3,0),i1=nd1a,nd1b),i2=nd2a,nd2b),i3=nd3a,nd3b)

 err=0.
 nn=-1 ! counts points on the interface
 ! =============== start loops ======================
 beginLoopsMask2d() 

   nn=nn+1

   ! here is the normal (assumed to be the same on both sides)
   an1=rsxy1(i1,i2,i3,axis1,0)   ! normal (an1,an2)
   an2=rsxy1(i1,i2,i3,axis1,1)
   aNorm=max(epsx,sqrt(an1**2+an2**2))
   an1=an1/aNorm
   an2=an2/aNorm
   tau1=-an2
   tau2= an1

   ! Evaluate the jump conditions using the wrong values at the ghost points 
   evaluateDispersiveInterfaceEquations2dOrder4()


     ! here is the matrix of coefficients for the unknowns u1(-1),v1(-1),u2(-1),v2(-1)
     ! Solve:
     !     
     !       A [ U ] = A [ U(old) ] - [ f ]
     !      u1r4(i1,i2,i3,kd)=(8.*(u1(i1+1,i2,i3,kd)-u1(i1-1,i2,i3,kd))-(u1(
     !     & i1+2,i2,i3,kd)-u1(i1-2,i2,i3,kd)))*dr114(0)
     !      u1x42(i1,i2,i3,kd)= rsxy1(i1,i2,i3,0,0)*u1r4(i1,i2,i3,kd)+rsxy1(
     !     & i1,i2,i3,1,0)*u1s4(i1,i2,i3,kd)
     !      u1y42(i1,i2,i3,kd)= rsxy1(i1,i2,i3,0,1)*u1r4(i1,i2,i3,kd)+rsxy1(
     !     & i1,i2,i3,1,1)*u1s4(i1,i2,i3,kd)
     !          a4(0,0) = -is*rsxy1(i1,i2,i3,axis1,0)/(2.*dr1(axis1))    ! coeff of u1(-1) from [u.x+v.y] 
     !          a4(0,1) = -is*rsxy1(i1,i2,i3,axis1,1)/(2.*dr1(axis1))    ! coeff of v1(-1) from [u.x+v.y] 
     !
     !          a4(2,0) =  is*rsxy1(i1,i2,i3,axis1,1)/(2.*dr1(axis1))   ! coeff of u1(-1) from [v.x - u.y] 
     !          a4(2,1) = -is*rsxy1(i1,i2,i3,axis1,0)/(2.*dr1(axis1))   ! coeff of v1(-1) from [v.x - u.y] 
     !
     !          a4(0,2) =  js*rsxy2(j1,j2,j3,axis2,0)/(2.*dr2(axis2))    ! coeff of u2(-1) from [u.x+v.y] 
     !          a4(0,3) =  js*rsxy2(j1,j2,j3,axis2,1)/(2.*dr2(axis2))    ! coeff of v2(-1) from [u.x+v.y] 
     !
     !          a4(2,2) = -js*rsxy2(j1,j2,j3,axis2,1)/(2.*dr2(axis2))   ! coeff of u2(-1) from [v.x - u.y] 
     !          a4(2,3) =  js*rsxy2(j1,j2,j3,axis2,0)/(2.*dr2(axis2))   ! coeff of v2(-1) from [v.x - u.y] 


     ! write(debugFile,'(" interface:E: initialized,it=",2i4)') initialized,it
     if( .false. .or. (initialized.eq.0 .and. it.eq.1) )then
       ! form the matrix (and save factor for later use)
       if( nn.eq.0 )then
         write(*,'(" Interface42c: form matrix and factor, it=",i4)') it
       end if

       ! Equation 0: 
       ! 0  [ u.x + v.y ] = 0
       aa8(0,0,0,nn) = -is*8.*rsxy1(i1,i2,i3,axis1,0)*dr114(axis1)     ! coeff of u1(-1) from [u.x+v.y] 
       aa8(0,1,0,nn) = -is*8.*rsxy1(i1,i2,i3,axis1,1)*dr114(axis1)     ! coeff of v1(-1) from [u.x+v.y] 
       aa8(0,4,0,nn) =  is*   rsxy1(i1,i2,i3,axis1,0)*dr114(axis1)     ! u1(-2)
       aa8(0,5,0,nn) =  is*   rsxy1(i1,i2,i3,axis1,1)*dr114(axis1)     ! v1(-2) 

       aa8(0,2,0,nn) =  js*8.*rsxy2(j1,j2,j3,axis2,0)*dr214(axis2)     ! coeff of u2(-1) from [u.x+v.y] 
       aa8(0,3,0,nn) =  js*8.*rsxy2(j1,j2,j3,axis2,1)*dr214(axis2)  
       aa8(0,6,0,nn) = -js*   rsxy2(j1,j2,j3,axis2,0)*dr214(axis2) 
       aa8(0,7,0,nn) = -js*   rsxy2(j1,j2,j3,axis2,1)*dr214(axis2)  

     ! 1  [ u.xx + u.yy ] = 0
     ! this macro comes from deriv.maple
     ! return the coefficient of u(-1) in uxxx+uxyy
     !#defineMacro lapCoeff4a(is,dr,ds) ((-1/3.*rxx*is-1/3.*ryy*is)/dr+(4/3.*rx**2+4/3.*ry**2)/dr**2)
     
     ! return the coefficient of u(-2) in uxxx+uxyy
     !#defineMacro lapCoeff4b(is,dr,ds) ((1/24.*rxx*is+1/24.*ryy*is)/dr+(-1/12.*rx**2-1/12.*ry**2)/dr**2 )

       ! optimize me ** June 27, 2016: 
       setJacobian( aj1, axis1)

       dr0=dr1(axis1)
       ds0=dr1(axis1p1)
       aLap0 = lapCoeff4a(is,dr0,ds0)
       aLap1 = lapCoeff4b(is,dr0,ds0)

       ! if( avoidInterfaceIterations.eq.1 )then
         ! lag cross terms (do not put into the matrix)
       !   ds0=dsBig
       ! end if

       ! dr1a(0:2) = dsBig in tangential directions if avoidInterfaceIterations=1
       ds0 =dr1a(axis1p1)
       aLapSq0 = lapSqCoeff4a(is,dr0,ds0)
       aLapSq1 = lapSqCoeff4b(is,dr0,ds0)

       setJacobian( aj2, axis2)
       dr0=dr2(axis2)
       ds0=dr2(axis2p1)
       bLap0 = lapCoeff4a(js,dr0,ds0)
       bLap1 = lapCoeff4b(js,dr0,ds0)

       ! if( avoidInterfaceIterations.eq.1 )then
         ! lag cross terms (do not put into the matrix)
       !   ds0=dsBig
       ! end if
       ! dr2a(0:2) = dsBig in tangential directions if avoidInterfaceIterations=1
       ds0 = dr2a(axis2p1)
       bLapSq0 = lapSqCoeff4a(js,dr0,ds0)
       bLapSq1 = lapSqCoeff4b(js,dr0,ds0)

      if( debug.gt.8 )then
       aa8(1,0,0,nn) = 16.*dx142(axis1)         ! coeff of u1(-1) from [u.xx + u.yy]
       aa8(1,4,0,nn) =    -dx142(axis1)         ! coeff of u1(-2) from [u.xx + u.yy]
        write(debugFile,'(" 4th: lap4: aLap0: rect=",e12.4," curv=",e12.4)') aLap0,aa8(1,0,0,nn)
        ! '
        write(debugFile,'(" 4th: lap4: aLap1: rect=",e12.4," curv=",e12.4)') aLap1,aa8(1,4,0,nn)
        ! '
      end if

      ! Equation 1:
      aa8(1,0,0,nn) = an1*aLap0       ! coeff of u1(-1) from [n.(u.xx + u.yy)]
      aa8(1,1,0,nn) = an2*aLap0 
      aa8(1,4,0,nn) = an1*aLap1       ! coeff of u1(-2) from [n.(u.xx + u.yy)]
      aa8(1,5,0,nn) = an2*aLap1  
       
      aa8(1,2,0,nn) =-an1*bLap0       ! coeff of u2(-1) from [n.(u.xx + u.yy)]
      aa8(1,3,0,nn) =-an2*bLap0
      aa8(1,6,0,nn) =-an1*bLap1       ! coeff of u2(-2) from [n.(u.xx + u.yy)]
      aa8(1,7,0,nn) =-an2*bLap1

      ! Equation 2: 
      ! 2  [ v.x - u.y ] =0 
      !          a8(2,0) =  is*8.*ry1*dx114(axis1)
      !          a8(2,1) = -is*8.*rx1*dx114(axis1)    ! coeff of v1(-1) from [v.x - u.y] 
      !          a8(2,4) = -is*   ry1*dx114(axis1)
      !          a8(2,5) =  is*   rx1*dx114(axis1)
      !          a8(2,2) = -js*8.*ry2*dx214(axis2)
      !          a8(2,3) =  js*8.*rx2*dx214(axis2)
      !          a8(2,6) =  js*   ry2*dx214(axis2)
      !          a8(2,7) = -js*   rx2*dx214(axis2)

       curl1um1 =  is*8.*rsxy1(i1,i2,i3,axis1,1)*dr114(axis1)   ! coeff of u(-1) from v.x - u.y 
       curl1vm1 = -is*8.*rsxy1(i1,i2,i3,axis1,0)*dr114(axis1)   ! coeff of v(-1) from v.x - u.y 
       curl1um2 = -is*   rsxy1(i1,i2,i3,axis1,1)*dr114(axis1)   ! coeff of u(-2) from v.x - u.y 
       curl1vm2 =  is*   rsxy1(i1,i2,i3,axis1,0)*dr114(axis1)   ! coeff of v(-2) from v.x - u.y

       curl2um1 =  js*8.*rsxy2(j1,j2,j3,axis2,1)*dr214(axis2)   ! coeff of u(-1) from v.x - u.y 
       curl2vm1 = -js*8.*rsxy2(j1,j2,j3,axis2,0)*dr214(axis2)   ! coeff of v(-1) from v.x - u.y 
       curl2um2 = -js*   rsxy2(j1,j2,j3,axis2,1)*dr214(axis2)   ! coeff of u(-2) from v.x - u.y 
       curl2vm2 =  js*   rsxy2(j1,j2,j3,axis2,0)*dr214(axis2)   ! coeff of v(-2) from v.x - u.y

       aa8(2,0,0,nn) =  curl1um1
       aa8(2,1,0,nn) =  curl1vm1
       aa8(2,4,0,nn) =  curl1um2
       aa8(2,5,0,nn) =  curl1vm2

       aa8(2,2,0,nn) = -curl2um1  
       aa8(2,3,0,nn) = -curl2vm1    
       aa8(2,6,0,nn) = -curl2um2  
       aa8(2,7,0,nn) = -curl2vm2 

       ! aa8(2,0,0,nn) =  is*8.*rsxy1(i1,i2,i3,axis1,1)*dr114(axis1)    
       ! aa8(2,1,0,nn) = -is*8.*rsxy1(i1,i2,i3,axis1,0)*dr114(axis1)    
       ! aa8(2,4,0,nn) = -is*   rsxy1(i1,i2,i3,axis1,1)*dr114(axis1)       
       ! aa8(2,5,0,nn) =  is*   rsxy1(i1,i2,i3,axis1,0)*dr114(axis1)       

       ! aa8(2,2,0,nn) = -js*8.*rsxy2(j1,j2,j3,axis2,1)*dr214(axis2)  
       ! aa8(2,3,0,nn) =  js*8.*rsxy2(j1,j2,j3,axis2,0)*dr214(axis2)    
       ! aa8(2,6,0,nn) =  js*   rsxy2(j1,j2,j3,axis2,1)*dr214(axis2)  
       ! aa8(2,7,0,nn) = -js*   rsxy2(j1,j2,j3,axis2,0)*dr214(axis2) 

       ! -------------- Equation 3 -----------------------
       !   [ tau.{ (uv.xx+uv.yy)/eps -alphaP*P.tt } ] = 0
       !    P.tt = c4PttLEsum * L(E) + c4PttLLEsum* L^2(E) + ...

       aa8(3,0,0,nn) = tau1*( aLap0*( 1./epsmu1 -alphaP1*c4PttLEsum1/epsmu1 ) - aLapSq0*alphaP1*c4PttLLEsum1/epsmu1**2 )
       aa8(3,1,0,nn) = tau2*( aLap0*( 1./epsmu1 -alphaP1*c4PttLEsum1/epsmu1 ) - aLapSq0*alphaP1*c4PttLLEsum1/epsmu1**2 )
       aa8(3,4,0,nn) = tau1*( aLap1*( 1./epsmu1 -alphaP1*c4PttLEsum1/epsmu1 ) - aLapSq1*alphaP1*c4PttLLEsum1/epsmu1**2 )
       aa8(3,5,0,nn) = tau2*( aLap1*( 1./epsmu1 -alphaP1*c4PttLEsum1/epsmu1 ) - aLapSq1*alphaP1*c4PttLLEsum1/epsmu1**2 )

       aa8(3,2,0,nn) =-tau1*( bLap0*( 1./epsmu2 -alphaP2*c4PttLEsum2/epsmu2 ) - bLapSq0*alphaP2*c4PttLLEsum2/epsmu2**2 )
       aa8(3,3,0,nn) =-tau2*( bLap0*( 1./epsmu2 -alphaP2*c4PttLEsum2/epsmu2 ) - bLapSq0*alphaP2*c4PttLLEsum2/epsmu2**2 )
       aa8(3,6,0,nn) =-tau1*( bLap1*( 1./epsmu2 -alphaP2*c4PttLEsum2/epsmu2 ) - bLapSq1*alphaP2*c4PttLLEsum2/epsmu2**2 )
       aa8(3,7,0,nn) =-tau2*( bLap1*( 1./epsmu2 -alphaP2*c4PttLEsum2/epsmu2 ) - bLapSq1*alphaP2*c4PttLLEsum2/epsmu2**2 )


      ! -------------- Equation 4 -----------------------
      !    [ (u.xx+u.yy).x + (v.xx+v.yy).y ] = 0

      setJacobian( aj1, axis1)

      ! dr1a(0:2) = dsBig in tangential directions if avoidInterfaceIterations=1
      dr0=dr1a(axis1)
      ds0=dr1a(axis1p1)
      ! if( avoidInterfaceIterations.eq.1 )then
        ! lag cross terms (do not put into the matrix)
      !   ds0=dsBig
      ! end if
      aLapX0 = xLapCoeff4a(is,dr0,ds0)
      aLapX1 = xLapCoeff4b(is,dr0,ds0)

      bLapY0 = yLapCoeff4a(is,dr0,ds0)
      bLapY1 = yLapCoeff4b(is,dr0,ds0)

      setJacobian( aj2, axis2)

      ! dr2a(0:2) = dsBig in tangential directions if avoidInterfaceIterations=1
      dr0=dr2a(axis2)
      ds0=dr2a(axis2p1)
      ! if( avoidInterfaceIterations.eq.1 )then
        ! lag cross terms (do not put into the matrix)
      !   ds0=dsBig
      ! end if
      cLapX0 = xLapCoeff4a(js,dr0,ds0)
      cLapX1 = xLapCoeff4b(js,dr0,ds0)

      dLapY0 = yLapCoeff4a(js,dr0,ds0)
      dLapY1 = yLapCoeff4b(js,dr0,ds0)


      ! 4  [ (u.xx+u.yy).x + (v.xx+v.yy).y ] = 0
      if( debug.gt.8 )then
      aa8(4,0,0,nn)= ( is*rx1*2.*dx122(axis1)*dx112(axis1)+is*rx1*2.*dx122(1)/(2.*dx1(0)))
      aa8(4,1,0,nn)= ( is*ry1*2.*dx122(axis1)*dx112(axis1)+is*ry1*2.*dx122(0)/(2.*dx1(1)))
      aa8(4,4,0,nn)= (-is*rx1   *dx122(axis1)*dx112(axis1) )  
      aa8(4,5,0,nn)= (-is*ry1   *dx122(axis1)*dx112(axis1))
        write(debugFile,'(" 4th: xlap4: aLapX0: rect=",e12.4," curv=",e12.4)') aLapX0,aa8(4,0,0,nn)
        write(debugFile,'(" 4th: xlap4: aLapX1: rect=",e12.4," curv=",e12.4)') aLapX1,aa8(4,4,0,nn)
        write(debugFile,'(" 4th: ylap4: bLapY0: rect=",e12.4," curv=",e12.4)') bLapY0,aa8(4,1,0,nn)
        write(debugFile,'(" 4th: ylap4: bLapY1: rect=",e12.4," curv=",e12.4)') bLapY1,aa8(4,5,0,nn)
        ! '
      end if

      aa8(4,0,0,nn)= aLapX0
      aa8(4,1,0,nn)= bLapY0
      aa8(4,4,0,nn)= aLapX1
      aa8(4,5,0,nn)= bLapY1

      aa8(4,2,0,nn)=-cLapX0
      aa8(4,3,0,nn)=-dLapY0
      aa8(4,6,0,nn)=-cLapX1
      aa8(4,7,0,nn)=-dLapY1

      ! ---------------- Equation 5 (2nd-order) -----------------

      !   [ ( {(Delta v).x - (Delta u).y}/(epsmu) - alphaP*( Py.ttx - Px.tty) )/mu ] =0 
      !
      !     P.tt = c2PttLEsum * L(E)
      eqnCoeff = ( 1./epsmu1 - alphaP1*c2PttLEsum1/epsmu1 )/mu1 
      eqnCoeffb = -alphaP1*c2PttEsum1/mu1 ! added sept 16, 2018 
      aa8(5,0,0,nn)=-bLapY0*eqnCoeff + curl1um1*eqnCoeffb  
      aa8(5,1,0,nn)= aLapX0*eqnCoeff + curl1vm1*eqnCoeffb
      aa8(5,4,0,nn)=-bLapY1*eqnCoeff + curl1um2*eqnCoeffb 
      aa8(5,5,0,nn)= aLapX1*eqnCoeff + curl1vm2*eqnCoeffb

      eqnCoeff = ( 1./epsmu2 - alphaP2*c2PttLEsum2/epsmu2 )/mu2 
      eqnCoeffb = -alphaP2*c2PttEsum2/mu2 ! added sept 16, 2018 
      aa8(5,2,0,nn)=-(-dLapY0*eqnCoeff + curl2um1*eqnCoeffb)
      aa8(5,3,0,nn)=-( cLapX0*eqnCoeff + curl2vm1*eqnCoeffb)
      aa8(5,6,0,nn)=-(-dLapY1*eqnCoeff + curl2um2*eqnCoeffb)
      aa8(5,7,0,nn)=-( cLapX1*eqnCoeff + curl2vm2*eqnCoeffb)


       ! ------- Equation 6 -----
       !  [ nv.( c^2*Delta^2(E) - alphaP*Delta(Ptt) )/mu ] = 0 

       ! 6  [ n.Delta^2 u/eps ] = 0

       if( debug.gt.8 )then
         aa8(6,0,0,nn) = -(4./(dx1(axis1)**4) +4./(dx1(0)**2*dx1(1)**2) )
         aa8(6,4,0,nn) =   1./(dx1(axis1)**4)
         write(debugFile,'(" 4th: lapSq: aLapSq0: rect=",e12.4," curv=",e12.4)') aLapSq0,aa8(6,0,0,nn)
         ! '
         write(debugFile,'(" 4th: lapSq: aLapSq1: rect=",e12.4," curv=",e12.4)') aLapSq1,aa8(6,4,0,nn)
         ! '
       end if

       if( setDivergenceAtInterfaces.eq.0 )then
        ! use Eqn 6 
        ! NOTE: LE = c^2*Delta(E) and LLE = (c^4*Delta^2) E 
        ! Note: the coeff of L(E) in Delta(Ptt) is the coeff of E in Ptt
        ! Note: the coeff of LL(E) in Delta(Ptt) is the coeff of LE in Ptt
        aa8(6,0,0,nn) = an1*( aLapSq0/epsmu1 -alphaP1*( c2PttEsum1*aLap0 + c2PttLEsum1*aLapSq0/epsmu1 ) )/mu1
        aa8(6,1,0,nn) = an2*( aLapSq0/epsmu1 -alphaP1*( c2PttEsum1*aLap0 + c2PttLEsum1*aLapSq0/epsmu1 ) )/mu1
        aa8(6,4,0,nn) = an1*( aLapSq1/epsmu1 -alphaP1*( c2PttEsum1*aLap1 + c2PttLEsum1*aLapSq1/epsmu1 ) )/mu1
        aa8(6,5,0,nn) = an2*( aLapSq1/epsmu1 -alphaP1*( c2PttEsum1*aLap1 + c2PttLEsum1*aLapSq1/epsmu1 ) )/mu1

        aa8(6,2,0,nn) =-an1*( bLapSq0/epsmu2 -alphaP2*( c2PttEsum2*bLap0 + c2PttLEsum2*bLapSq0/epsmu2 ) )/mu2
        aa8(6,3,0,nn) =-an2*( bLapSq0/epsmu2 -alphaP2*( c2PttEsum2*bLap0 + c2PttLEsum2*bLapSq0/epsmu2 ) )/mu2
        aa8(6,6,0,nn) =-an1*( bLapSq1/epsmu2 -alphaP2*( c2PttEsum2*bLap1 + c2PttLEsum2*bLapSq1/epsmu2 ) )/mu2
        aa8(6,7,0,nn) =-an2*( bLapSq1/epsmu2 -alphaP2*( c2PttEsum2*bLap1 + c2PttLEsum2*bLapSq1/epsmu2 ) )/mu2
       end if

       if( setDivergenceAtInterfaces.eq.1 )then
         ! Set div(E)=0 as equation 6
        aa8(6,0,0,nn) = -is*8.*rsxy1(i1,i2,i3,axis1,0)*dr114(axis1)     ! coeff of u1(-1) from [u.x+v.y] 
        aa8(6,1,0,nn) = -is*8.*rsxy1(i1,i2,i3,axis1,1)*dr114(axis1)     ! coeff of v1(-1) from [u.x+v.y] 
        aa8(6,4,0,nn) =  is*   rsxy1(i1,i2,i3,axis1,0)*dr114(axis1)     ! u1(-2)
        aa8(6,5,0,nn) =  is*   rsxy1(i1,i2,i3,axis1,1)*dr114(axis1)     ! v1(-2) 

        aa8(6,2,0,nn) = 0.
        aa8(6,3,0,nn) = 0.
        aa8(6,6,0,nn) = 0.
        aa8(6,7,0,nn) = 0.
       end if 

       ! ------- Equation 7 ------
       ! [ tv.( c^4*Delta^2(E) - alphaP*c^2*Delta(P.tt) - alphaP*P.tttt) ]=0 

       ! 7  [ tau.Delta^2 v/eps^2 ] = 0 
       ! Note: the coeff of L(E) in Delta(Ptt) is the coeff of E in Ptt
       ! Note: the coeff of LL(E) in Delta(Ptt) is the coeff of LE in Ptt
       coeffLap1   =              -alphaP1*(  c2PttEsum1 + c2PttttLEsum1  )/epsmu1
       coeffLapSq1 = 1./epsmu1**2 -alphaP1*( c2PttLEsum1 + c2PttttLLEsum1 )/epsmu1**2

       coeffLap2   =              -alphaP2*(  c2PttEsum2 + c2PttttLEsum2  )/epsmu2
       coeffLapSq2 = 1./epsmu2**2 -alphaP2*( c2PttLEsum2 + c2PttttLLEsum2 )/epsmu2**2

       aa8(7,0,0,nn) = tau1*( coeffLapSq1*aLapSq0 + coeffLap1*aLap0 )
       aa8(7,1,0,nn) = tau2*( coeffLapSq1*aLapSq0 + coeffLap1*aLap0 )
       aa8(7,4,0,nn) = tau1*( coeffLapSq1*aLapSq1 + coeffLap1*aLap1 )
       aa8(7,5,0,nn) = tau2*( coeffLapSq1*aLapSq1 + coeffLap1*aLap1 )

       aa8(7,2,0,nn) =-tau1*( coeffLapSq2*bLapSq0 + coeffLap2*bLap0 )
       aa8(7,3,0,nn) =-tau2*( coeffLapSq2*bLapSq0 + coeffLap2*bLap0 )
       aa8(7,6,0,nn) =-tau1*( coeffLapSq2*bLapSq1 + coeffLap2*bLap1 )
       aa8(7,7,0,nn) =-tau2*( coeffLapSq2*bLapSq1 + coeffLap2*bLap1 )

!       aa8(7,0,0,nn) = tau1*( aLapSq0/epsmu1**2 -alphaP1*( c2PttEsum1*aLap0/epsmu1 +c2PttLEsum1*aLapSq0/epsmu1 + c2PttttLEsum1*aLap0  ))
!       aa8(7,1,0,nn) = tau2*( aLapSq0/epsmu1**2 -alphaP1*( c2PttEsum1*aLap0/epsmu1 +c2PttLEsum1*aLapSq0/epsmu1 + c2PttttLEsum1*aLap0  ))
!       aa8(7,4,0,nn) = tau1*( aLapSq1/epsmu1**2 -alphaP1*( c2PttEsum1*aLap1/epsmu1 +c2PttLEsum1*aLapSq1/epsmu1 + c2PttttLEsum1*aLap1  ))
!       aa8(7,5,0,nn) = tau2*( aLapSq1/epsmu1**2 -alphaP1*( c2PttEsum1*aLap1/epsmu1 +c2PttLEsum1*aLapSq1/epsmu1 + c2PttttLEsum1*aLap1  ))
!
!       aa8(7,2,0,nn) =-tau1*( bLapSq0/epsmu2**2 -alphaP2*( c2PttEsum2*bLap0/epsmu2 +c2PttLEsum2*bLapSq0/epsmu2 + c2PttttLEsum2*bLap0  ))
!       aa8(7,3,0,nn) =-tau2*( bLapSq0/epsmu2**2 -alphaP2*( c2PttEsum2*bLap0/epsmu2 +c2PttLEsum2*bLapSq0/epsmu2 + c2PttttLEsum2*bLap0  ))
!       aa8(7,6,0,nn) =-tau1*( bLapSq1/epsmu2**2 -alphaP2*( c2PttEsum2*bLap1/epsmu2 +c2PttLEsum2*bLapSq1/epsmu2 + c2PttttLEsum2*bLap1  ))
!       aa8(7,7,0,nn) =-tau2*( bLapSq1/epsmu2**2 -alphaP2*( c2PttEsum2*bLap1/epsmu2 +c2PttLEsum2*bLapSq1/epsmu2 + c2PttttLEsum2*bLap1  ))

       ! save a copy of the matrix
       do n2=0,7
       do n1=0,7
         aa8(n1,n2,1,nn)=aa8(n1,n2,0,nn)
       end do
       end do

       ! --- check matrix coefficients by delta function approach ----
       if( checkCoeff.eq.1 )then
        numberOfEquations=8
        do n2=0,7
        do n1=0,7
          a8(n1,n2)=aa8(n1,n2,0,nn)
        end do
        end do                
        checkCoefficients(i1,i2,i3, j1,j2,j3,numberOfEquations,a8,evalDispersiveInterfaceEquations24c )
       end if


       ! solve A Q = F
       ! factor the matrix
       numberOfEquations=8
       call dgeco( aa8(0,0,0,nn), numberOfEquations, numberOfEquations, ipvt8(0,nn),rcond,work(0))

       if( debug.gt.3 ) write(debugFile,'(" --> 4cth: i1,i2=",2i4," rcond=",e10.2)') i1,i2,rcond
       ! '
     end if


     ! Save current solution to compare to new
     q(0) = u1(i1-is1,i2-is2,i3,ex)
     q(1) = u1(i1-is1,i2-is2,i3,ey)
     q(2) = u2(j1-js1,j2-js2,j3,ex)
     q(3) = u2(j1-js1,j2-js2,j3,ey)

     q(4) = u1(i1-2*is1,i2-2*is2,i3,ex)
     q(5) = u1(i1-2*is1,i2-2*is2,i3,ey)
     q(6) = u2(j1-2*js1,j2-2*js2,j3,ex)
     q(7) = u2(j1-2*js1,j2-2*js2,j3,ey)

      if( debug.gt.4 ) write(debugFile,'(" --> 4cth: i1,i2=",2i4," q=",8e10.2)') i1,i2,(q(n),n=0,7)

     ! subtract off the contributions from the initial (wrong) values at the ghost points:
     do n=0,7
       f(n) = (aa8(n,0,1,nn)*q(0)+aa8(n,1,1,nn)*q(1)+aa8(n,2,1,nn)*q(2)+aa8(n,3,1,nn)*q(3)+\
               aa8(n,4,1,nn)*q(4)+aa8(n,5,1,nn)*q(5)+aa8(n,6,1,nn)*q(6)+aa8(n,7,1,nn)*q(7)) - f(n)
     end do

                          ! '

     ! solve A Q = F
     job=0
     numberOfEquations=8
     call dgesl( aa8(0,0,0,nn), numberOfEquations, numberOfEquations, ipvt8(0,nn), f(0), job)


     u1(i1-is1,i2-is2,i3,ex)=(1.-omega)*u1(i1-is1,i2-is2,i3,ex) + omega*f(0)
     u1(i1-is1,i2-is2,i3,ey)=(1.-omega)*u1(i1-is1,i2-is2,i3,ey) + omega*f(1)
     u2(j1-js1,j2-js2,j3,ex)=(1.-omega)*u2(j1-js1,j2-js2,j3,ex) + omega*f(2)
     u2(j1-js1,j2-js2,j3,ey)=(1.-omega)*u2(j1-js1,j2-js2,j3,ey) + omega*f(3)

     u1(i1-2*is1,i2-2*is2,i3,ex)=(1.-omega)*u1(i1-2*is1,i2-2*is2,i3,ex) + omega*f(4)
     u1(i1-2*is1,i2-2*is2,i3,ey)=(1.-omega)*u1(i1-2*is1,i2-2*is2,i3,ey) + omega*f(5)
     u2(j1-2*js1,j2-2*js2,j3,ex)=(1.-omega)*u2(j1-2*js1,j2-2*js2,j3,ex) + omega*f(6)
     u2(j1-2*js1,j2-2*js2,j3,ey)=(1.-omega)*u2(j1-2*js1,j2-2*js2,j3,ey) + omega*f(7)

     ! compute the maximum change in the solution for this iteration
     do n=0,7
       err=max(err,abs(q(n)-f(n)))
     end do

     if( debug.gt.4 )then
      write(debugFile,'(" --> 4cth: i1,i2=",2i4," f(solve)=",8e10.2)') i1,i2,(f(n),n=0,7)
      write(debugFile,'(" --> 4cth: i1,i2=",2i4,"      f-q=",8e10.2,"  err=",e10.2)') i1,i2,(f(n)-q(n),n=0,7),err
     end if
     ! '


    if( .false. .and. knownSolutionOption.eq.userDefinedKnownSolution )then
       ! --- TEST ---
       ! Evaluate the user defined known solution
       numberOfTimeDerivatives=0
       call evalUserDefinedKnownSolution( t, grid1,i1-is1,i2-is2,i3,evals,pvals, numberOfTimeDerivatives )
       write(*,'(" interface:ghost: grid,i1,i2=",3i4," Ex,Ey=",2(1pe12.4)," known=",2(1pe12.4)," err=",2(1pe12.2))') grid1,i1-is1,i2-is2,u1(i1-is1,i2-is2,i3,ex),u1(i1-is1,i2-is2,i3,ey),\
                     evals(0),evals(1),u1(i1-is1,i2-is2,i3,ex)-evals(0),u1(i1-is1,i2-is2,i3,ey)-evals(1)
       if( numberOfPolarizationVectors1.gt.0 )then
         write(*,'(" Ghost: P=",2(1pe12.4)," known=",2(1pe12.4)," err=",2(1pe12.2))') p1(i1-is1,i2-is2,i3,0),p1(i1-is1,i2-is2,i3,1),\
                     pvals(0),pvals(1),p1(i1-is1,i2-is2,i3,0)-pvals(0),p1(i1-is1,i2-is2,i3,ey)-pvals(1)
       end if
    end if


    if( .false. .and. debug.gt.2 )then 

      ! --- check residuals in the jump conditions ----

      ! Evaluate the jump conditions using the new values at the ghost points 
      evaluateDispersiveInterfaceEquations2dOrder4()
 

      write(debugFile,'(" JUMP-residuals: i1,i2=",2i4," f(re-eval)=",8e10.2)') i1,i2,f(0),f(1),f(2),f(3),f(4),f(5),f(6),f(7)

    end if

     ! ******************************************************
     ! solve for Hz
     !  [ w.n/eps ] = 0
     !  [ lap(w)/eps ] = 0
     !  [ lap(w).n/eps**2 ] = 0
     !  [ lapSq(w)/eps**2 ] = 0

     ! first evaluate the equations we want to solve with the wrong values at the ghost points:
     evalMagneticDerivs2dOrder4()
     evalMagneticField2dJumpOrder4()

     if( .false. .or. (initialized.eq.0 .and. it.eq.1) )then
       ! form the matrix for computing Hz (and save factor for later use)

       ! 1: [ w.n/eps ] = 0
       a0 = (an1*rsxy1(i1,i2,i3,axis1,0)+an2*rsxy1(i1,i2,i3,axis1,1))*dr114(axis1)/eps1
       b0 = (an1*rsxy2(j1,j2,j3,axis2,0)+an2*rsxy2(j1,j2,j3,axis2,1))*dr214(axis2)/eps2
       aa4(0,0,0,nn) = -is*8.*a0
       aa4(0,2,0,nn) =  is*   a0
       aa4(0,1,0,nn) =  js*8.*b0
       aa4(0,3,0,nn) = -js*   b0

       ! 2: [ lap(w)/eps ] = 0 
       aa4(1,0,0,nn) = aLap0/eps1
       aa4(1,2,0,nn) = aLap1/eps1
       aa4(1,1,0,nn) =-bLap0/eps2
       aa4(1,3,0,nn) =-bLap1/eps2

       ! 3  [ (an1*(w.xx+w.yy).x + an2.(w.xx+w.yy).y)/eps**2 ] = 0
       aa4(2,0,0,nn)= (an1*aLapX0+an2*bLapY0)/eps1**2
       aa4(2,2,0,nn)= (an1*aLapX1+an2*bLapY1)/eps1**2
       aa4(2,1,0,nn)=-(an1*cLapX0+an2*dLapY0)/eps2**2
       aa4(2,3,0,nn)=-(an1*cLapX1+an2*dLapY1)/eps2**2

       ! 4 [ lapSq(w)/eps**2 ] = 0 
       aa4(3,0,0,nn) = aLapSq0/eps1**2
       aa4(3,2,0,nn) = aLapSq1/eps1**2
       aa4(3,1,0,nn) =-bLapSq0/eps2**2
       aa4(3,3,0,nn) =-bLapSq1/eps2**2

       ! save a copy of the matrix
       do n2=0,3
       do n1=0,3
         aa4(n1,n2,1,nn)=aa4(n1,n2,0,nn)
       end do
       end do

       ! factor the matrix
       numberOfEquations=4
       call dgeco( aa4(0,0,0,nn), numberOfEquations, numberOfEquations, ipvt4(0,nn),rcond,work(0))
     end if

     q(0) = u1(i1-is1,i2-is2,i3,hz)
     q(1) = u2(j1-js1,j2-js2,j3,hz)
     q(2) = u1(i1-2*is1,i2-2*is2,i3,hz)
     q(3) = u2(j1-2*js1,j2-2*js2,j3,hz)

     ! subtract off the contributions from the wrong values at the ghost points:
     do n=0,3
       f(n) = (aa4(n,0,1,nn)*q(0)+aa4(n,1,1,nn)*q(1)+aa4(n,2,1,nn)*q(2)+aa4(n,3,1,nn)*q(3)) - f(n)
     end do
     ! solve
     numberOfEquations=4
     job=0
     call dgesl( aa4(0,0,0,nn), numberOfEquations, numberOfEquations, ipvt4(0,nn), f(0), job)

     u1(i1-  is1,i2-  is2,i3,hz)=(1.-omega)*u1(i1-  is1,i2-  is2,i3,hz) + omega*f(0)
     u2(j1-  js1,j2-  js2,j3,hz)=(1.-omega)*u2(j1-  js1,j2-  js2,j3,hz) + omega*f(1)
     u1(i1-2*is1,i2-2*is2,i3,hz)=(1.-omega)*u1(i1-2*is1,i2-2*is2,i3,hz) + omega*f(2)
     u2(j1-2*js1,j2-2*js2,j3,hz)=(1.-omega)*u2(j1-2*js1,j2-2*js2,j3,hz) + omega*f(3)

    ! compute the maximum change in the solution for this iteration
    if( .false. )then
      ! ************** TURN OFF UNTIL WE FIX Hz *************************
      do n=0,3
        err=max(err,abs(q(n)-f(n)))
      end do
    end if

    if( debug.gt.0 )then ! re-evaluate

     evalMagneticDerivs2dOrder4()
     evalMagneticField2dJumpOrder4()

     if( debug.gt.3 ) write(debugFile,'(" --> 4cth: i1,i2=",2i4," hz-f(re-eval)=",8e10.2)') i1,i2,f(0),f(1),f(2),f(3)
       ! '
    end if


     ! ***********************

     ! u1(i1-is1,i2-is2,i3,hz)=u2(j1+js1,j2+js2,j3,hz) 
     ! u2(j1-js1,j2-js2,j3,hz)=u1(i1+is1,i2+is2,i3,hz)
     ! u1(i1-2*is1,i2-2*is2,i3,hz)=u2(j1+2*js1,j2+2*js2,j3,hz) 
     ! u2(j1-2*js1,j2-2*js2,j3,hz)=u1(i1+2*is1,i2+2*is2,i3,hz)

 endLoopsMask2d()
 ! =============== end loops =======================
      
 if( checkCoeff.eq.1 )then
   write(*,'("+++++ iGDM24c: check coeff in interface: max(diff) = ",1pe8.2)') coeffDiff
 end if

#endMacro
