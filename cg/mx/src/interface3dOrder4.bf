! *******************************************************************************
!   Interface boundary conditions - Fourth Order and 3D
!      See also interface3d.bf 
! *******************************************************************************

! These next include files will define the macros that will define the difference approximations
! The actual macro is called below
!* #Include "defineDiffNewerOrder2f.h"
!* #Include "defineDiffNewerOrder4f.h"

! These next include file will define the macros that will define the difference approximations (in op/src)
! Defines getDuDx2(u,aj,ff), getDuDxx2(u,aj,ff), getDuDx3(u,aj,ff), ...  etc. 
#Include "derivMacroDefinitions.h"

! Define 
!    defineParametricDerivativeMacros(u,dr,dx,DIM,ORDER,COMPONENTS,MAXDERIV)
!       defines -> ur2, us2, ux2, uy2, ...            (2D)
!                  ur3, us3, ut3, ux3, uy3, uz3, ...  (3D)
#Include "defineParametricDerivMacros.h"

! defineParametricDerivativeMacros(u,dr,dx,DIM,ORDER,COMPONENTS,MAXDERIV)
! 2D, order=6, components=1
! defineParametricDerivativeMacros(u,dr,dx,DIM,ORDER,COMPONENTS,MAXDERIV)

!- defineParametricDerivativeMacros(rsxy1,dr1,dx1,3,2,2,4)
!- defineParametricDerivativeMacros(rsxy1,dr1,dx1,3,4,2,2)
!- defineParametricDerivativeMacros(u1,dr1,dx1,3,2,1,4)
!- defineParametricDerivativeMacros(u1,dr1,dx1,3,4,1,2)
!-
!- defineParametricDerivativeMacros(rsxy2,dr2,dx2,3,2,2,4)
!- defineParametricDerivativeMacros(rsxy2,dr2,dx2,3,4,2,2)
!- defineParametricDerivativeMacros(u2,dr2,dx2,3,2,1,4)
!- defineParametricDerivativeMacros(u2,dr2,dx2,3,4,1,2)

 defineParametricDerivativeMacros(rsxy1,dr1,dx1,3,2,2,4)
 defineParametricDerivativeMacros(rsxy1,dr1,dx1,3,4,2,2)
 defineParametricDerivativeMacros(u1,dr1,dx1,3,2,1,4)
 defineParametricDerivativeMacros(u1,dr1,dx1,3,4,1,2)

 ! We need up to 2 derivatives of u1n to order 2 -- do order 4 for TZ
 defineParametricDerivativeMacros(u1n,dr1,dx1,3,4,1,2)
 defineParametricDerivativeMacros(u1n,dr1,dx1,3,2,1,2)

 ! We need up to 2 derivatives of p1 to order 2-- do order 4 for TZ
 defineParametricDerivativeMacros(p1,dr1,dx1,3,4,1,2)
 defineParametricDerivativeMacros(p1n,dr1,dx1,3,4,1,2)
 defineParametricDerivativeMacros(p1,dr1,dx1,3,2,1,2)
 defineParametricDerivativeMacros(p1n,dr1,dx1,3,2,1,2)

 defineParametricDerivativeMacros(rsxy2,dr2,dx2,3,2,2,4)
 defineParametricDerivativeMacros(rsxy2,dr2,dx2,3,4,2,2)
 defineParametricDerivativeMacros(u2,dr2,dx2,3,2,1,4)
 defineParametricDerivativeMacros(u2,dr2,dx2,3,4,1,2)

 ! We need up to 2 derivatives of u1n to order 2-- do order 4 for TZ
 defineParametricDerivativeMacros(u2n,dr2,dx2,3,4,1,2)
 defineParametricDerivativeMacros(u2n,dr2,dx2,3,2,1,2)

 ! We need up to 2 derivatives of p2 to order 2-- do order 4 for TZ
 defineParametricDerivativeMacros(p2,dr1,dx1,3,4,1,2)
 defineParametricDerivativeMacros(p2n,dr1,dx1,3,4,1,2)
 defineParametricDerivativeMacros(p2,dr1,dx1,3,2,1,2)
 defineParametricDerivativeMacros(p2n,dr1,dx1,3,2,1,2)

! ===========================================================================================
! Macro: Output some debug info for the first few time-steps 
! ===========================================================================================
#beginMacro INFO(string)
if( t.le.3.*dt .and. debug.gt.0 )then
  write(*,'("Interface>>>",string)')
end if
#endMacro


! ******************************************************************************************************************
! ************* These are altered version of those from insImp.h ***************************************************
! ******************************************************************************************************************


! ==========================================================================================
!  Evaluate the Jacobian and its derivatives (parametric and spatial). 
!    rsxy   : jacobian matrix name 
!    aj     : prefix for the name of the resulting jacobian variables, 
!             e.g. ajrx, ajsy, ajrxx, ajsxy, ...
!    MAXDER : number of derivatives to evaluate.  
! ==========================================================================================
#beginMacro opEvalJacobianDerivatives(rsxy,i1,i2,i3,aj,MAXDER)

#If $GRIDTYPE eq "curvilinear"
 ! this next call will define the jacobian and its derivatives (parameteric and spatial)
 #peval evalJacobianDerivatives(rsxy,i1,i2,i3,aj,$DIM,$ORDER,MAXDER)

#End

#endMacro 

! ==========================================================================================
!  Evaluate the parametric derivatives of u.
!    u      : evaluate derivatives of this function.
!    uc     : component to evaluate
!    uu     : prefix for the name of the resulting derivatives, e.g. uur, uus, uurr, ...
!    MAXDER : number of derivatives to evaluate.  
! ==========================================================================================
#beginMacro opEvalParametricDerivative(u,i1,i2,i3,uc,uu,MAXDER)
#If $GRIDTYPE eq "curvilinear" 
 #peval evalParametricDerivativesComponents1(u,i1,i2,i3,uc, uu,$DIM,$ORDER,MAXDER)
#Else
 uu=u(i1,i2,i3,uc) ! in the rectangular case just eval the solution
#End
#endMacro


! ==========================================================================================
!  Evaluate a derivative. (assumes parametric derivatives have already been evaluated)
!   DERIV   : name of the derivative. One of 
!                x,y,z,xx,xy,xz,...
!    u      : evaluate derivatives of this function.
!    uc     : component to evaluate
!    uu     : prefix for the name of the resulting derivatives (same name used with opEvalParametricDerivative) 
!    aj     : prefix for the name of the jacobian variables.
!    ud     : derivative is assigned to this variable.
! ==========================================================================================
#beginMacro getOp(DERIV, u,i1,i2,i3,uc,uu,aj,ud )

 #If $GRIDTYPE eq "curvilinear" 
  #peval getDuD ## DERIV ## $DIM(uu,aj,ud)  ! Note: The perl variables are evaluated when the macro is USED. 
 #Else
  #peval ud = u ## DERIV ## $ORDER(i1,i2,i3,uc)
 #End

#endMacro

! ******************************************************************************************************************

! loop over the boundary points
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


! loop over the boundary points
#beginMacro beginLoops2d()
 i3=n3a
 j3=m3a

 j2=m2a
 do i2=n2a,n2b
  j1=m1a
  do i1=n1a,n1b
#endMacro
#beginMacro endLoops2d()
   j1=j1+1
  end do
  j2=j2+1
 end do
#endMacro

! loop over the boundary points with a mask. 
! Assign pts where both mask1 and mask2 are discretization pts.
! If mask1>0 and mask2<0 then we just leave the extrapolated values in u1 and u2 .
#beginMacro beginLoopsMask2d()
 i3=n3a
 j3=m3a

 j2=m2a
 do i2=n2a,n2b
  j1=m1a
  do i1=n1a,n1b
   if( mask1(i1,i2,i3).gt.0 .and. mask2(j1,j2,j3).gt.0 )then
#endMacro
#beginMacro endLoopsMask2d()
   end if
   j1=j1+1
  end do
  j2=j2+1
 end do
#endMacro

! loop over the boundary points that includes ghost points in the tangential direction
#beginMacro beginGhostLoops2d()
 i3=n3a
 j3=m3a
 j2=mm2a
 do i2=nn2a,nn2b
  j1=mm1a
  do i1=nn1a,nn1b
#endMacro

! loop over the boundary points that includes ghost points in the tangential direction.
! Assign pts where both mask1 and mask2 are discretization pts.
! If mask1>0 and mask2<0 then we just leave the extrapolated values in u1 and u2 .
#beginMacro beginGhostLoopsMask2d()
 i3=n3a
 j3=m3a
 j2=mm2a
 do i2=nn2a,nn2b
  j1=mm1a
  do i1=nn1a,nn1b
  if( mask1(i1,i2,i3).gt.0 .and. mask2(j1,j2,j3).gt.0 )then
#endMacro

#beginMacro beginLoops3d()
 j3=m3a
 do i3=n3a,n3b
 j2=m2a
 do i2=n2a,n2b
 j1=m1a
 do i1=n1a,n1b
#endMacro
#beginMacro endLoops3d()
   j1=j1+1
  end do
  j2=j2+1
 end do
  j3=j3+1
 end do
#endMacro

! Assign pts where both mask1 and mask2 are discretization pts.
! If mask1>0 and mask2<0 then we just leave the extrapolated values in u1 and u2 .
#beginMacro beginLoopsMask3d()
 j3=m3a
 do i3=n3a,n3b
 j2=m2a
 do i2=n2a,n2b
 j1=m1a
 do i1=n1a,n1b
 if( mask1(i1,i2,i3).gt.0 .and. mask2(j1,j2,j3).gt.0 )then
#endMacro
#beginMacro endLoopsMask3d()
  end if
   j1=j1+1
  end do
  j2=j2+1
 end do
  j3=j3+1
 end do
#endMacro

#beginMacro beginGhostLoops3d()
 j3=mm3a
 do i3=nn3a,nn3b
 j2=mm2a
 do i2=nn2a,nn2b
  j1=mm1a
  do i1=nn1a,nn1b
#endMacro

#beginMacro beginGhostLoopsMask3d()
 j3=mm3a
 do i3=nn3a,nn3b
 j2=mm2a
 do i2=nn2a,nn2b
  j1=mm1a
  do i1=nn1a,nn1b
  if( mask1(i1,i2,i3).gt.0 .and. mask2(j1,j2,j3).gt.0 )then
#endMacro

#defineMacro extrap2(uu,k1,k2,k3,kc,ks1,ks2,ks3) \
            (2.*uu(k1,k2,k3,kc)-uu(k1+ks1,k2+ks2,k3+ks3,kc))

#defineMacro extrap3(uu,k1,k2,k3,kc,ks1,ks2,ks3) \
            (3.*uu(k1,k2,k3,kc)-3.*uu(k1+ks1,k2+ks2,k3+ks3,kc)\
            +   uu(k1+2*ks1,k2+2*ks2,k3+2*ks3,kc))

#defineMacro extrap4(uu,k1,k2,k3,kc,ks1,ks2,ks3) \
            (4.*uu(k1,k2,k3,kc)-6.*uu(k1+ks1,k2+ks2,k3+ks3,kc)\
            +4.*uu(k1+2*ks1,k2+2*ks2,k3+2*ks3,kc)-uu(k1+3*ks1,k2+3*ks2,k3+3*ks3,kc))

#defineMacro extrap5(uu,k1,k2,k3,kc,ks1,ks2,ks3) \
            (5.*uu(k1,k2,k3,kc)-10.*uu(k1+ks1,k2+ks2,k3+ks3,kc)\
            +10.*uu(k1+2*ks1,k2+2*ks2,k3+2*ks3,kc)-5.*uu(k1+3*ks1,k2+3*ks2,k3+3*ks3,kc)\
            +uu(k1+4*ks1,k2+4*ks2,k3+4*ks3,kc))


! *********************************************************************
! ********** MACROS FOR DISPERSIVE INTERFACE CONDITIONS ***************
! *********************************************************************
#Include "dispersiveInterfaceMacros.h"

! ********************************************************************************
!       INTERFACE MACROS (used in interface3d.bf and interface3dOrder4.bf)
! ********************************************************************************
#Include "interfaceMacros.h"


! ------------------ OLD ------------------------------
! This macro will assign the jump conditions on the boundary
! DIM (input): number of dimensions (2 or 3)
! GRIDTYPE (input) : curvilinear or rectangular
#beginMacro boundaryJumpConditionsOLD(DIM,GRIDTYPE)
 #If #DIM eq "2"
  if( eps1.lt.eps2 )then
    epsRatio=eps1/eps2
    beginGhostLoopsMask2d()
      ! eps2 n.u2 = eps1 n.u1
      !     tau.u2 = tau.u1

      #If #GRIDTYPE eq "curvilinear"
       an1=rsxy1(i1,i2,i3,axis1,0)   ! normal (an1,an2)
       an2=rsxy1(i1,i2,i3,axis1,1)
       aNorm=max(epsx,sqrt(an1**2+an2**2))
       an1=an1/aNorm
       an2=an2/aNorm
      #Elif #GRIDTYPE eq "rectangular"
       an1=an1Cartesian
       an2=an2Cartesian
      #Else
         stop 1111
      #End
      ua=u1(i1,i2,i3,ex)
      ub=u1(i1,i2,i3,ey)
      nDotU = an1*ua+an2*ub
      if( twilightZone.eq.1 )then
       ! adjust for TZ forcing (here we assume the exact solution is the same on both sides)
       call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, ue )
       call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, ve )
       nDotU = nDotU - (an1*ue+an2*ve)
      end if

      ! u2 equals u1 but with normal component = eps1/eps2*(n.u1)
      u2(j1,j2,j3,ex) = ua + (nDotU*epsRatio - nDotU)*an1
      u2(j1,j2,j3,ey) = ub + (nDotU*epsRatio - nDotU)*an2
      u2(j1,j2,j3,hz) = u1(i1,i2,i3,hz)


    endLoopsMask2d()
  else
    epsRatio=eps2/eps1
    beginGhostLoopsMask2d()
      ! eps2 n.u2 = eps1 n.u1
      !     tau.u2 = tau.u1

      #If #GRIDTYPE eq "curvilinear"
       an1=rsxy1(i1,i2,i3,axis1,0)
       an2=rsxy1(i1,i2,i3,axis1,1)
       aNorm=max(epsx,sqrt(an1**2+an2**2))
       an1=an1/aNorm
       an2=an2/aNorm
      #Elif #GRIDTYPE eq "rectangular"
       an1=an1Cartesian
       an2=an2Cartesian
      #Else
        stop 1112
      #End
      ua=u2(j1,j2,j3,ex)
      ub=u2(j1,j2,j3,ey)

      nDotU = an1*ua+an2*ub
      if( twilightZone.eq.1 )then
       ! adjust for TZ forcing (here we assume the exact solution is the same on both sides)
       call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, ue )
       call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, ve )
! write(*,'(" jump: x,y=",2e10.2," ua,ue=",2e10.2," ub,ve=",2e10.2)') xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),ua,ue,ub,ve
       nDotU = nDotU - (an1*ue+an2*ve)
      end if

      u1(i1,i2,i3,ex) = ua + (nDotU*epsRatio - nDotU)*an1
      u1(i1,i2,i3,ey) = ub + (nDotU*epsRatio - nDotU)*an2
      u1(i1,i2,i3,hz) = u2(j1,j2,j3,hz)
    endLoopsMask2d()
  end if
 #Else
  ! *** 3D ***
  if( eps1.lt.eps2 )then
    epsRatio=eps1/eps2
    beginGhostLoopsMask3d()
      ! eps2 n.u2 = eps1 n.u1
      !     tau.u2 = tau.u1

      #If #GRIDTYPE eq "curvilinear"
       an1=rsxy1(i1,i2,i3,axis1,0)   ! normal (an1,an2)
       an2=rsxy1(i1,i2,i3,axis1,1)
       an3=rsxy1(i1,i2,i3,axis1,2)
       aNorm=max(epsx,sqrt(an1**2+an2**2+an3**2))
       an1=an1/aNorm
       an2=an2/aNorm
       an3=an3/aNorm
      #Elif #GRIDTYPE eq "rectangular"
       an1=an1Cartesian
       an2=an2Cartesian
       an3=an3Cartesian
      #Else
         stop 1111
      #End
      ua=u1(i1,i2,i3,ex)
      ub=u1(i1,i2,i3,ey)
      uc=u1(i1,i2,i3,ez)
      nDotU = an1*ua+an2*ub+an3*uc
      if( twilightZone.eq.1 )then
       ! adjust for TZ forcing (here we assume the exact solution is the same on both sides)
       call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ex, ue )
       call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ey, ve )
       call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ez, we )
       nDotU = nDotU - (an1*ue+an2*ve+an3*we)
      end if
      ! u2 equals u1 but with normal component = eps1/eps2*(n.u1)
      u2(j1,j2,j3,ex) = ua + (nDotU*epsRatio - nDotU)*an1
      u2(j1,j2,j3,ey) = ub + (nDotU*epsRatio - nDotU)*an2
      u2(j1,j2,j3,ez) = uc + (nDotU*epsRatio - nDotU)*an3

!   write(*,'(" jump(1): (i1,i2,i3)=",3i3," j1,j2,j3=",3i3)') i1,i2,i3,j1,j2,j3
!   write(*,'(" jump(1): x,y,z=",3e10.2," ua,ue=",2e10.2," ub,ve=",2e10.2," uc,we=",2e10.2)') xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),ua,ue,ub,ve,uc,we

    endLoopsMask3d()
  else
    epsRatio=eps2/eps1
    beginGhostLoopsMask3d()
      ! eps2 n.u2 = eps1 n.u1
      !     tau.u2 = tau.u1

      #If #GRIDTYPE eq "curvilinear"
       an1=rsxy1(i1,i2,i3,axis1,0)   ! normal (an1,an2)
       an2=rsxy1(i1,i2,i3,axis1,1)
       an3=rsxy1(i1,i2,i3,axis1,2)
       aNorm=max(epsx,sqrt(an1**2+an2**2+an3**2))
       an1=an1/aNorm
       an2=an2/aNorm
       an3=an3/aNorm
      #Elif #GRIDTYPE eq "rectangular"
       an1=an1Cartesian
       an2=an2Cartesian
       an3=an3Cartesian
      #Else
         stop 1111
      #End
      ua=u2(j1,j2,j3,ex)
      ub=u2(j1,j2,j3,ey)
      uc=u2(j1,j2,j3,ez)

      nDotU = an1*ua+an2*ub+an3*uc
      if( twilightZone.eq.1 )then
       ! adjust for TZ forcing (here we assume the exact solution is the same on both sides)
       call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ex, ue )
       call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ey, ve )
       call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ez, we )
       nDotU = nDotU - (an1*ue+an2*ve+an3*we)

!   write(*,'(" jump(2): (i1,i2,i3)=",3i3," j1,j2,j3=",3i3)') i1,i2,i3,j1,j2,j3
!   write(*,'(" jump(2): x,y,z=",3e10.2," ua,ue=",2e10.2," ub,ve=",2e10.2," uc,we=",2e10.2)') xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),ua,ue,ub,ve,uc,we

      end if

      u1(i1,i2,i3,ex) = ua + (nDotU*epsRatio - nDotU)*an1
      u1(i1,i2,i3,ey) = ub + (nDotU*epsRatio - nDotU)*an2
      u1(i1,i2,i3,ez) = uc + (nDotU*epsRatio - nDotU)*an3
    endLoopsMask3d()
  end if
 #End
#endMacro


! ********************************************************************************
!     Usage: setJacobianRS( aj1, r, s)
!            setJacobianRS( aj1, s, r)
! ********************************************************************************
#beginMacro setJacobianRS(aj, R, S)
 rx   =aj ## R ## x
 ry   =aj ## R ## y
                    
 rxx  =aj ## R ## xx  
 rxy  =aj ## R ## xy  
 ryy  =aj ## R ## yy  
                    
 rxxx =aj ## R ## xxx 
 rxxy =aj ## R ## xxy 
 rxyy =aj ## R ## xyy 
 ryyy =aj ## R ## yyy 
                    
 rxxxx=aj ## R ## xxxx
 rxxyy=aj ## R ## xxyy
 ryyyy=aj ## R ## yyyy

 sx   =aj ## S ## x   
 sy   =aj ## S ## y   
                    
 sxx  =aj ## S ## xx  
 sxy  =aj ## S ## xy  
 syy  =aj ## S ## yy  
                    
 sxxx =aj ## S ## xxx 
 sxxy =aj ## S ## xxy 
 sxyy =aj ## S ## xyy 
 syyy =aj ## S ## yyy 
                    
 sxxxx=aj ## S ## xxxx
 sxxyy=aj ## S ## xxyy
 syyyy=aj ## S ## yyyy

#endMacro

! ***************************************************************************
! This macro will set the temp variables rx, rxx, ry, ryx, ...
! If axis=0 then
!   rx = ajrx
!   sx = ajsx
!    ...
!  else if axis=1
!    -- permute r <-> s 
!   rx = ajsx
!   sx = ajrx
!    ...
! ***************************************************************************
#beginMacro setJacobian(aj, axis)
if( axis.eq.0 )then
 setJacobianRS( aj, r, s)
else
 setJacobianRS( aj, s, r)
end if

#endMacro

! ===================================================================================
!  Optimized periodic update: (only applied in serial)
!     update the periodic ghost points used by an interface on the grid face (side,axis)
! ===================================================================================
#beginMacro periodicUpdate2d(u,bc,gid,side,axis)
if( parallel.eq.0 )then
 axisp1=mod(axis+1,nd)
 if( bc(0,axisp1).lt.0 )then
  ! direction axisp1 is periodic
  diff(axis)=0
  diff(axisp1)=gid(1,axisp1)-gid(0,axisp1)

  if( side.eq.0 )then
    ! assign 4 ghost points outside lower corner
    np1a=gid(0,0)-2
    np1b=gid(0,0)-1
    np2a=gid(0,1)-2
    np2b=gid(0,1)-1

    beginLoops(np1a,np1b,np2a,np2b,n3a,n3b,ex,hz)
     u(i1,i2,i3,n) = u(i1+diff(0),i2+diff(1),i3,n)
    endLoops()

    ! assign 4 ghost points outside upper corner
    if( axis.eq.0 )then
      np2a=gid(1,axisp1)+1
      np2b=gid(1,axisp1)+2
    else
      np1a=gid(1,axisp1)+1
      np1b=gid(1,axisp1)+2
    end if

    beginLoops(np1a,np1b,np2a,np2b,n3a,n3b,ex,hz)
     u(i1,i2,i3,n) = u(i1-diff(0),i2-diff(1),i3,n)
    endLoops()

  else

    ! assign 4 ghost points outside upper corner
    np1a=gid(1,0)+1
    np1b=gid(1,0)+2
    np2a=gid(1,1)+1
    np2b=gid(1,1)+2

    beginLoops(np1a,np1b,np2a,np2b,n3a,n3b,ex,hz)
     u(i1,i2,i3,n) = u(i1-diff(0),i2-diff(1),i3,n)
    endLoops()

    if( axis.eq.0 )then
      np2a=gid(0,axisp1)-2
      np2b=gid(0,axisp1)-1
    else
      np1a=gid(0,axisp1)-2
      np1b=gid(0,axisp1)-1
    end if

    beginLoops(np1a,np1b,np2a,np2b,n3a,n3b,ex,hz)
     u(i1,i2,i3,n) = u(i1+diff(0),i2+diff(1),i3,n)
    endLoops()
  end if

 endif
end if

#endMacro

! ===================================================================================
!  Optimized periodic update:
!     update the periodic ghost points used by an interface on the grid face (side,axis)
! ===================================================================================
#beginMacro periodicUpdate3d(u,bc,gid,side,axis)
if( parallel.eq.0 )then
 axisp1=mod(axis+1,nd)
 axisp2=mod(axis+2,nd)
 if( bc(0,axisp1).lt.0 .or. bc(0,axisp2).lt.0 )then
  ! this periodic update should be done in the C++ calling routine
  ! write(*,'("periodicUpdate3d: Warning: opt version not implemented")')
  ! stop 
 end if
end if
#endMacro


! ******************************************************************************
!   This next macro is called by other macros to evaluate the first and second derivatives
!   This macro assumes that opEvalJacobianDerivatives has been called
! ******************************************************************************
#beginMacro evalSecondDerivs3d(rsxy1,aj1,u1,i1,i2,i3,ex,uu1,uu)
 opEvalParametricDerivative(u1,i1,i2,i3,ex,uu1,2)    ! computes uu1r, uu1s 
 getOp(x ,u1,i1,i2,i3,ex,uu1,aj1,uu ## x)            ! u1.x
 getOp(y ,u1,i1,i2,i3,ex,uu1,aj1,uu ## y)            ! u1.y
 getOp(z ,u1,i1,i2,i3,ex,uu1,aj1,uu ## z)            ! u1.z
 getOp(xx,u1,i1,i2,i3,ex,uu1,aj1,uu ## xx)
 getOp(yy,u1,i1,i2,i3,ex,uu1,aj1,uu ## yy)
 getOp(zz,u1,i1,i2,i3,ex,uu1,aj1,uu ## zz)
 uu ## Lap = uu ## xx+ uu ## yy+ uu ## zz
#endMacro 


! *********************************************************************************
!   Evaluate derivatives for the 2nd-order 3D interface equations
! *********************************************************************************
#beginMacro evalInterfaceDerivatives3d()
 ! NOTE: the jacobian derivatives can be computed once for all components
 opEvalJacobianDerivatives(rsxy1,i1,i2,i3,aj1,1)
 evalSecondDerivs3d(rsxy1,aj1,u1,i1,i2,i3,ex,uu1,u1)
 evalSecondDerivs3d(rsxy1,aj1,u1,i1,i2,i3,ey,vv1,v1)
 evalSecondDerivs3d(rsxy1,aj1,u1,i1,i2,i3,ez,ww1,w1)

 ! NOTE: the jacobian derivatives can be computed once for all components
 opEvalJacobianDerivatives(rsxy2,j1,j2,j3,aj2,1)
 evalSecondDerivs3d(rsxy2,aj2,u2,j1,j2,j3,ex,uu2,u2)
 evalSecondDerivs3d(rsxy2,aj2,u2,j1,j2,j3,ey,vv2,v2)
 evalSecondDerivs3d(rsxy2,aj2,u2,j1,j2,j3,ez,ww2,w2)
#endMacro


#beginMacro vectorSymmetry(u1,i1,i2,i3,is1,is2,is3)
  ! first set all components to be even 
  u1(i1-is1,i2-is2,i3-is3,ex)=u1(i1+is1,i2+is2,i3+is3,ex)
  u1(i1-is1,i2-is2,i3-is3,ey)=u1(i1+is1,i2+is2,i3+is3,ey)
  u1(i1-is1,i2-is2,i3-is3,ez)=u1(i1+is1,i2+is2,i3+is3,ez)

  ! odd symmetry for the normal component: 
  nDotUm = an1*u1(i1-is1,i2-is2,i3-is3,ex)+an2*u1(i1-is1,i2-is2,i3-is3,ey)+an3*u1(i1-is1,i2-is2,i3-is3,ez)
  nDotU0 = an1*u1(i1    ,i2    ,i3    ,ex)+an2*u1(i1    ,i2    ,i3    ,ey)+an3*u1(i1    ,i2    ,i3    ,ez)
  nDotUp = an1*u1(i1+is1,i2+is2,i3+is3,ex)+an2*u1(i1+is1,i2+is2,i3+is3,ey)+an3*u1(i1+is1,i2+is2,i3+is3,ez)
  nDotU = 2.*nDotU0-nDotUp - nDotUm
  u1(i1-is1,i2-is2,i3-is3,ex)=u1(i1-is1,i2-is2,i3-is3,ex) +  nDotU*an1
  u1(i1-is1,i2-is2,i3-is3,ey)=u1(i1-is1,i2-is2,i3-is3,ey) +  nDotU*an2
  u1(i1-is1,i2-is2,i3-is3,ez)=u1(i1-is1,i2-is2,i3-is3,ez) +  nDotU*an3
#endMacro

#beginMacro eval3dJumpOrder4()
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

 f(3)=( u1Lap/(epsmu1) + cem1*nDotLapE1*an1 ) - ( u2Lap/(epsmu2) + cem2*nDotLapE2*an1 )
 f(4)=( v1Lap/(epsmu1) + cem1*nDotLapE1*an2 ) - ( v2Lap/(epsmu2) + cem2*nDotLapE2*an2 )
 f(5)=( w1Lap/(epsmu1) + cem1*nDotLapE1*an3 ) - ( w2Lap/(epsmu2) + cem2*nDotLapE2*an3 )

 ! For now we extrap the 2nd ghost line 
 f(6)= u1(i1-2*is1,i2-2*is2,i3-2*is3,ex) - extrap5(u1,i1-is1,i2-is2,i3-is3,ex,is1,is2,is3)
 f(7)= u1(i1-2*is1,i2-2*is2,i3-2*is3,ey) - extrap5(u1,i1-is1,i2-is2,i3-is3,ey,is1,is2,is3)
 f(8)= u1(i1-2*is1,i2-2*is2,i3-2*is3,ez) - extrap5(u1,i1-is1,i2-is2,i3-is3,ez,is1,is2,is3)   

 f(9)= u2(j1-2*js1,j2-2*js2,j3-2*js3,ex)-extrap5(u2,j1-js1,j2-js2,j3-js3,ex,js1,js2,js3)
 f(10)=u2(j1-2*js1,j2-2*js2,j3-2*js3,ey)-extrap5(u2,j1-js1,j2-js2,j3-js3,ey,js1,js2,js3)
 f(11)=u2(j1-2*js1,j2-2*js2,j3-2*js3,ez)-extrap5(u2,j1-js1,j2-js2,j3-js3,ez,js1,js2,js3)

 if( twilightZone.eq.1 )then

   call ogderiv(ep, 0,1,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ex, uex  )
   call ogderiv(ep, 0,0,1,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ex, uey  )
   call ogderiv(ep, 0,0,0,1, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ex, uez  )
   call ogderiv(ep, 0,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ex, uexx )
   call ogderiv(ep, 0,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ex, ueyy )
   call ogderiv(ep, 0,0,0,2, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ex, uezz )

   call ogderiv(ep, 0,1,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ey, vex  )
   call ogderiv(ep, 0,0,1,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ey, vey  )
   call ogderiv(ep, 0,0,0,1, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ey, vez  )
   call ogderiv(ep, 0,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ey, vexx )
   call ogderiv(ep, 0,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ey, veyy )
   call ogderiv(ep, 0,0,0,2, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ey, vezz )

   call ogderiv(ep, 0,1,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ez, wex  )
   call ogderiv(ep, 0,0,1,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ez, wey  )
   call ogderiv(ep, 0,0,0,1, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ez, wez  )
   call ogderiv(ep, 0,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ez, wexx )
   call ogderiv(ep, 0,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ez, weyy )
   call ogderiv(ep, 0,0,0,2, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ez, wezz )

   ueLap = uexx+ueyy+uezz
   veLap = vexx+veyy+vezz
   weLap = wexx+weyy+wezz

   curlEex = wey-vez
   curlEey = uez-wex
   curlEez = vex-uey
   nDotCurlEe=an1*curlEex+an2*curlEey+an3*curlEez
   nDotLapEe=an1*ueLap+an2*veLap+an3*weLap

   f(0)= f(0) - ( (curlEex- nDotCurlEe*an1)*(1./mu1-1./mu2) )
   f(1)= f(1) - ( (curlEey- nDotCurlEe*an2)*(1./mu1-1./mu2) )
   f(2)= f(2) - ( (curlEez- nDotCurlEe*an3)*(1./mu1-1./mu2) )

   ! Assign extrap with TZ ??

   f(3)= f(3) - ( ueLap*(1./epsmu1-1./epsmu2) + nDotLapEe*an1*(cem1-cem2) )
   f(4)= f(4) - ( veLap*(1./epsmu1-1./epsmu2) + nDotLapEe*an2*(cem1-cem2) )
   f(5)= f(5) - ( weLap*(1./epsmu1-1./epsmu2) + nDotLapEe*an3*(cem1-cem2) ) 

 end if

#endMacro

! XXXXXXXXXXXXXXXXXXXXXXX from interface3d.bf :  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

! ******************************************************************************
!   This next macro is called by evalDerivs3dOrder4
!   This macro assumes that opEvalJacobianDerivatives has been called
! ******************************************************************************
#beginMacro evalFourthDerivs3d(rsxy1,aj1,u1,i1,i2,i3,ex,uu1,uu)
 opEvalParametricDerivative(u1,i1,i2,i3,ex,uu1,4)    ! computes uu1r, uu1s 
 getOp(xxx ,u1,i1,i2,i3,ex,uu1,aj1,uu ## xxx)       ! u1.xxx
 getOp(xxy ,u1,i1,i2,i3,ex,uu1,aj1,uu ## xxy)       ! u1.xxy
 getOp(xxz ,u1,i1,i2,i3,ex,uu1,aj1,uu ## xxz)
 getOp(xyy ,u1,i1,i2,i3,ex,uu1,aj1,uu ## xyy) 
 getOp(xzz ,u1,i1,i2,i3,ex,uu1,aj1,uu ## xzz) 

 getOp(yyy ,u1,i1,i2,i3,ex,uu1,aj1,uu ## yyy) 
 getOp(yyz ,u1,i1,i2,i3,ex,uu1,aj1,uu ## yyz) 
 getOp(yzz ,u1,i1,i2,i3,ex,uu1,aj1,uu ## yzz) 
 getOp(zzz ,u1,i1,i2,i3,ex,uu1,aj1,uu ## zzz) 

 getOp(xxxx,u1,i1,i2,i3,ex,uu1,aj1,uu ## xxxx) 
 getOp(xxyy,u1,i1,i2,i3,ex,uu1,aj1,uu ## xxyy) 
 getOp(xxzz,u1,i1,i2,i3,ex,uu1,aj1,uu ## xxzz) 
 getOp(yyzz,u1,i1,i2,i3,ex,uu1,aj1,uu ## yyzz) 
 getOp(yyyy,u1,i1,i2,i3,ex,uu1,aj1,uu ## yyyy) 
 getOp(zzzz,u1,i1,i2,i3,ex,uu1,aj1,uu ## zzzz) 

 uu ## LapSq = uu ## xxxx +2.*(uu ## xxyy + uu ## xxzz + uu ## yyzz ) + uu ## yyyy + uu ## zzzz
#endMacro 


! ******************************************************************************
!   Evaluate derivatives for the 4th-order 3D interface equations
! ******************************************************************************
#beginMacro evalDerivs3dOrder4()

#perl $ORDER=2;
 ! These derivatives are computed to 2nd-order accuracy

 ! NOTE: the jacobian derivatives can be computed once for all components
 opEvalJacobianDerivatives(rsxy1,i1,i2,i3,aj1,3)
 evalFourthDerivs3d(rsxy1,aj1,u1,i1,i2,i3,ex,uu1,u1)
 evalFourthDerivs3d(rsxy1,aj1,u1,i1,i2,i3,ey,vv1,v1)
 evalFourthDerivs3d(rsxy1,aj1,u1,i1,i2,i3,ez,ww1,w1)

 ! NOTE: the jacobian derivatives can be computed once for all components
 opEvalJacobianDerivatives(rsxy2,j1,j2,j3,aj2,3)
 evalFourthDerivs3d(rsxy2,aj2,u2,j1,j2,j3,ex,uu2,u2)
 evalFourthDerivs3d(rsxy2,aj2,u2,j1,j2,j3,ey,vv2,v2)
 evalFourthDerivs3d(rsxy2,aj2,u2,j1,j2,j3,ez,ww2,w2)

#perl $ORDER=4;
 ! These derivatives are computed to 4th-order accuracy

 ! NOTE: the jacobian derivatives can be computed once for all components
 opEvalJacobianDerivatives(rsxy1,i1,i2,i3,aj1,1)
 evalSecondDerivs3d(rsxy1,aj1,u1,i1,i2,i3,ex,uu1,u1)
 evalSecondDerivs3d(rsxy1,aj1,u1,i1,i2,i3,ey,vv1,v1)
 evalSecondDerivs3d(rsxy1,aj1,u1,i1,i2,i3,ez,ww1,w1)

 ! NOTE: the jacobian derivatives can be computed once for all components
 opEvalJacobianDerivatives(rsxy2,j1,j2,j3,aj2,1)
 evalSecondDerivs3d(rsxy2,aj2,u2,j1,j2,j3,ex,uu2,u2)
 evalSecondDerivs3d(rsxy2,aj2,u2,j1,j2,j3,ey,vv2,v2)
 evalSecondDerivs3d(rsxy2,aj2,u2,j1,j2,j3,ez,ww2,w2)

#endMacro


! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

! ======================================================================================
!  Evaluate the interface jump conditions when using the full centered scheme
! ======================================================================================
#beginMacro eval3dJumpOrder4New()
 ! ----- START EVAL3D JUMP ORDER 4 ------
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


 f(0)= ( divE1*an1 + (curlE1x- nDotCurlE1*an1)/mu1 ) - ( divE2*an1 + (curlE2x- nDotCurlE2*an1)/mu2 )
 f(1)= ( divE1*an2 + (curlE1y- nDotCurlE1*an2)/mu1 ) - ( divE2*an2 + (curlE2y- nDotCurlE2*an2)/mu2 )
 f(2)= ( divE1*an3 + (curlE1z- nDotCurlE1*an3)/mu1 ) - ( divE2*an3 + (curlE2z- nDotCurlE2*an3)/mu2 )

 f(3)= ( u1Lap/(epsmu1) + cem1*nDotLapE1*an1 ) - ( u2Lap/(epsmu2) + cem2*nDotLapE2*an1 )
 f(4)= ( v1Lap/(epsmu1) + cem1*nDotLapE1*an2 ) - ( v2Lap/(epsmu2) + cem2*nDotLapE2*an2 )
 f(5)= ( w1Lap/(epsmu1) + cem1*nDotLapE1*an3 ) - ( w2Lap/(epsmu2) + cem2*nDotLapE2*an3 )

 ! *new* 
 divLapE1 = u1xxx+v1xxy+w1xxz + u1xyy+v1yyy+w1yyz + u1xzz+v1yzz+w1zzz
 curlLapE1x = w1xxy-v1xxz + w1yyy-v1yyz + w1yzz-v1zzz
 curlLapE1y = u1xxz-w1xxx + u1yyz-w1xyy + u1zzz-w1xzz
 curlLapE1z = v1xxx-u1xxy + v1xyy-u1yyy + v1xzz-u1yzz
 nDotCurlLapE1=an1*curlLapE1x+an2*curlLapE1y+an3*curlLapE1z
 nDotLapSqE1 = an1*u1LapSq + an2*v1LapSq + an3*w1LapSq

 divLapE2 = u2xxx+v2xxy+w2xxz + u2xyy+v2yyy+w2yyz + u2xzz+v2yzz+w2zzz
 curlLapE2x = w2xxy-v2xxz + w2yyy-v2yyz + w2yzz-v2zzz
 curlLapE2y = u2xxz-w2xxx + u2yyz-w2xyy + u2zzz-w2xzz
 curlLapE2z = v2xxx-u2xxy + v2xyy-u2yyy + v2xzz-u2yzz
 nDotCurlLapE2=an1*curlLapE2x+an2*curlLapE2y+an3*curlLapE2z
 nDotLapSqE2 = an1*u2LapSq + an2*v2LapSq + an3*w2LapSq


 ! *wdh* Oct 20, 2018 -- set [div(Lap(E))], not [div(lap(E))/epsmu] 
 ! *new way* FINISH ME -- (1) form coeff's, (2) TZ forcing
 ! f(6)= ( divLapE1*an1 + (curlLapE1x- nDotCurlLapE1*an1)/mu1/(epsmu1) ) - ( divLapE2*an1 + (curlLapE2x- nDotCurlLapE2*an1)/mu2/(epsmu2) )
 ! f(7)= ( divLapE1*an2 + (curlLapE1y- nDotCurlLapE1*an2)/mu1/(epsmu1) ) - ( divLapE2*an2 + (curlLapE2y- nDotCurlLapE2*an2)/mu2/(epsmu2) )
 ! f(8)= ( divLapE1*an3 + (curlLapE1z- nDotCurlLapE1*an3)/mu1/(epsmu1) ) - ( divLapE2*an3 + (curlLapE2z- nDotCurlLapE2*an3)/mu2/(epsmu2) )

 ! *old way*
 f(6)= ( divLapE1*an1 + (curlLapE1x- nDotCurlLapE1*an1)/mu1 )/(epsmu1) - ( divLapE2*an1 + (curlLapE2x- nDotCurlLapE2*an1)/mu2 )/(epsmu2)
 f(7)= ( divLapE1*an2 + (curlLapE1y- nDotCurlLapE1*an2)/mu1 )/(epsmu1) - ( divLapE2*an2 + (curlLapE2y- nDotCurlLapE2*an2)/mu2 )/(epsmu2)
 f(8)= ( divLapE1*an3 + (curlLapE1z- nDotCurlLapE1*an3)/mu1 )/(epsmu1) - ( divLapE2*an3 + (curlLapE2z- nDotCurlLapE2*an3)/mu2 )/(epsmu2)

 f(9)= ( u1LapSq/(epsmu1) + cem1*nDotLapSqE1*an1 )/(epsmu1) - ( u2LapSq/(epsmu2) + cem2*nDotLapSqE2*an1 )/(epsmu2)
 f(10)=( v1LapSq/(epsmu1) + cem1*nDotLapSqE1*an2 )/(epsmu1) - ( v2LapSq/(epsmu2) + cem2*nDotLapSqE2*an2 )/(epsmu2)
 f(11)=( w1LapSq/(epsmu1) + cem1*nDotLapSqE1*an3 )/(epsmu1) - ( w2LapSq/(epsmu2) + cem2*nDotLapSqE2*an3 )/(epsmu2)

 ! ******* TESTING *****
 ! For testing extrap first ghost line
 ! *a396*
 ! f(0)= u1(i1-is1,i2-is2,i3-is3,ex) - extrap5(u1,i1,i2,i3,ex,is1,is2,is3)
 ! f(1)= u1(i1-is1,i2-is2,i3-is3,ey) - extrap5(u1,i1,i2,i3,ey,is1,is2,is3)
 ! f(2)= u1(i1-is1,i2-is2,i3-is3,ez) - extrap5(u1,i1,i2,i3,ez,is1,is2,is3)   

 ! f(0) = u1x - u2x
 ! f(1) = w1x - w2x
 ! f(2) = v1x - v2x

 ! f(3)= u2(j1-js1,j2-js2,j3-js3,ex) - extrap5(u2,j1,j2,j3,ex,js1,js2,js3)
 ! f(4)= u2(j1-js1,j2-js2,j3-js3,ey) - extrap5(u2,j1,j2,j3,ey,js1,js2,js3)
 ! f(5)= u2(j1-js1,j2-js2,j3-js3,ez) - extrap5(u2,j1,j2,j3,ez,js1,js2,js3)

 ! f(3) = u1xx - u2xx 
 ! f(4) = w1xx - w2xx 
 ! f(5) = v1xx - v2xx 

 ! For testing extrap the 2nd ghost line :
 ! *e678*
 ! f(6)= u1(i1-2*is1,i2-2*is2,i3-2*is3,ex) - extrap5(u1,i1-is1,i2-is2,i3-is3,ex,is1,is2,is3)
 ! f(7)= u1(i1-2*is1,i2-2*is2,i3-2*is3,ey) - extrap5(u1,i1-is1,i2-is2,i3-is3,ey,is1,is2,is3)
 ! f(8)= u1(i1-2*is1,i2-2*is2,i3-2*is3,ez) - extrap5(u1,i1-is1,i2-is2,i3-is3,ez,is1,is2,is3)   

 ! f(9)= u2(j1-2*js1,j2-2*js2,j3-2*js3,ex)-extrap5(u2,j1-js1,j2-js2,j3-js3,ex,js1,js2,js3)
 ! f(10)=u2(j1-2*js1,j2-2*js2,j3-2*js3,ey)-extrap5(u2,j1-js1,j2-js2,j3-js3,ey,js1,js2,js3)
 ! f(11)=u2(j1-2*js1,j2-2*js2,j3-2*js3,ez)-extrap5(u2,j1-js1,j2-js2,j3-js3,ez,js1,js2,js3)

 ! Try this Oct 23, 2018 TROUBLE 
 ! f(6) = u1xxx - u2xxx
 ! f(7) = w1xxx - w2xxx 
 ! f(8) = v1xxx - v2xxx 
 ! f(9) = (u1xxxx + 4.*u1xxyy + 0.*u1xxzz)  - (u2xxxx+ 4.*u2xxyy + 0.*u2xxzz)
 ! f(10)= (w1xxxx + 4.*w1xxyy + 0.*w1xxzz)  - (w2xxxx+ 4.*w2xxyy + 0.*w2xxzz)
 ! f(11)= (v1xxxx + 4.*v1xxyy + 0.*v1xxzz)  - (v2xxxx+ 4.*v2xxyy + 0.*v2xxzz)



 ! OK: 
 ! f(6) = u1xxx - u2xxx
 ! f(7) = (w1xxx + w1xyy) - (w2xxx + w2xyy)
 ! f(8) = (v1xxx + v1xyy) - (v2xxx + v2xyy) 
 ! f(9) = (u1xxxx )  - (u2xxxx)
 ! f(10)= (w1xxxx )  - (w2xxxx)
 ! f(11)= (v1xxxx )  - (v2xxxx)

 ! Almost curl(curl(curl)))
 ! f(6) = u1xxx - u2xxx
 ! No: 
 ! f(7) = (w1xxx + 1.*w1xyy + 1.*w1xzz ) - (w2xxx + 1.*w2xyy + 1.*w2xzz)
 ! f(8) = (v1xxx + 1.*v1xyy + 1.*v1xzz ) - (v2xxx + 1.*v2xyy + 1.*v2xzz) 


 ! LapSq() : trouble
 ! f(9) = (u1xxxx + 2.*u1xxyy + 2.*u1xxzz)  - (u2xxxx+ 2.*u2xxyy + 2.*u2xxzz)
 ! f(10)= (w1xxxx + 2.*w1xxyy + 2.*w1xxzz)  - (w2xxxx+ 2.*w2xxyy + 2.*w2xxzz)
 ! f(11)= (v1xxxx + 2.*v1xxyy + 2.*v1xxzz)  - (v2xxxx+ 2.*v2xxyy + 2.*v2xxzz)



 if( twilightZone.eq.1 )then

   x1=xy1(i1,i2,i3,0)
   y1=xy1(i1,i2,i3,1)
   z1=xy1(i1,i2,i3,2)
   call ogderiv(ep, 0,1,0,0, x1,y1,z1,t, ex, uex  )
   call ogderiv(ep, 0,0,1,0, x1,y1,z1,t, ex, uey  )
   call ogderiv(ep, 0,0,0,1, x1,y1,z1,t, ex, uez  )
   call ogderiv(ep, 0,2,0,0, x1,y1,z1,t, ex, uexx )
   call ogderiv(ep, 0,0,2,0, x1,y1,z1,t, ex, ueyy )
   call ogderiv(ep, 0,0,0,2, x1,y1,z1,t, ex, uezz )

   call ogderiv(ep, 0,1,0,0, x1,y1,z1,t, ey, vex  )
   call ogderiv(ep, 0,0,1,0, x1,y1,z1,t, ey, vey  )
   call ogderiv(ep, 0,0,0,1, x1,y1,z1,t, ey, vez  )
   call ogderiv(ep, 0,2,0,0, x1,y1,z1,t, ey, vexx )
   call ogderiv(ep, 0,0,2,0, x1,y1,z1,t, ey, veyy )
   call ogderiv(ep, 0,0,0,2, x1,y1,z1,t, ey, vezz )

   call ogderiv(ep, 0,1,0,0, x1,y1,z1,t, ez, wex  )
   call ogderiv(ep, 0,0,1,0, x1,y1,z1,t, ez, wey  )
   call ogderiv(ep, 0,0,0,1, x1,y1,z1,t, ez, wez  )
   call ogderiv(ep, 0,2,0,0, x1,y1,z1,t, ez, wexx )
   call ogderiv(ep, 0,0,2,0, x1,y1,z1,t, ez, weyy )
   call ogderiv(ep, 0,0,0,2, x1,y1,z1,t, ez, wezz )

   ! finish me: 
   call ogderiv(ep, 0,3,0,0, x1,y1,z1,t, ex, uexxx )
   call ogderiv(ep, 0,0,3,0, x1,y1,z1,t, ex, ueyyy )
   call ogderiv(ep, 0,0,0,3, x1,y1,z1,t, ex, uezzz )
   call ogderiv(ep, 0,2,1,0, x1,y1,z1,t, ex, uexxy )
   call ogderiv(ep, 0,2,0,1, x1,y1,z1,t, ex, uexxz )
   call ogderiv(ep, 0,0,2,1, x1,y1,z1,t, ex, ueyyz )

   call ogderiv(ep, 0,1,2,0, x1,y1,z1,t, ex, uexyy ) ! add
   call ogderiv(ep, 0,1,0,2, x1,y1,z1,t, ex, uexzz )
   call ogderiv(ep, 0,0,1,2, x1,y1,z1,t, ex, ueyzz )


   call ogderiv(ep, 0,4,0,0, x1,y1,z1,t, ex, uexxxx )
   call ogderiv(ep, 0,2,2,0, x1,y1,z1,t, ex, uexxyy )
   call ogderiv(ep, 0,2,0,2, x1,y1,z1,t, ex, uexxzz )
   call ogderiv(ep, 0,0,4,0, x1,y1,z1,t, ex, ueyyyy )
   call ogderiv(ep, 0,0,2,2, x1,y1,z1,t, ex, ueyyzz )
   call ogderiv(ep, 0,0,0,4, x1,y1,z1,t, ex, uezzzz )

   call ogderiv(ep, 0,3,0,0, x1,y1,z1,t, ey, vexxx )
   call ogderiv(ep, 0,0,3,0, x1,y1,z1,t, ey, veyyy )
   call ogderiv(ep, 0,0,0,3, x1,y1,z1,t, ey, vezzz )
   call ogderiv(ep, 0,2,1,0, x1,y1,z1,t, ey, vexxy )
   call ogderiv(ep, 0,2,0,1, x1,y1,z1,t, ey, vexxz )
   call ogderiv(ep, 0,0,2,1, x1,y1,z1,t, ey, veyyz )

   call ogderiv(ep, 0,1,2,0, x1,y1,z1,t, ey, vexyy ) ! add
   call ogderiv(ep, 0,1,0,2, x1,y1,z1,t, ey, vexzz )
   call ogderiv(ep, 0,0,1,2, x1,y1,z1,t, ey, veyzz )

   call ogderiv(ep, 0,4,0,0, x1,y1,z1,t, ey, vexxxx )
   call ogderiv(ep, 0,2,2,0, x1,y1,z1,t, ey, vexxyy )
   call ogderiv(ep, 0,2,0,2, x1,y1,z1,t, ey, vexxzz )
   call ogderiv(ep, 0,0,4,0, x1,y1,z1,t, ey, veyyyy )
   call ogderiv(ep, 0,0,2,2, x1,y1,z1,t, ey, veyyzz )
   call ogderiv(ep, 0,0,0,4, x1,y1,z1,t, ey, vezzzz )

   call ogderiv(ep, 0,3,0,0, x1,y1,z1,t, ez, wexxx )
   call ogderiv(ep, 0,0,3,0, x1,y1,z1,t, ez, weyyy )
   call ogderiv(ep, 0,0,0,3, x1,y1,z1,t, ez, wezzz )
   call ogderiv(ep, 0,2,1,0, x1,y1,z1,t, ez, wexxy )
   call ogderiv(ep, 0,2,0,1, x1,y1,z1,t, ez, wexxz )
   call ogderiv(ep, 0,0,2,1, x1,y1,z1,t, ez, weyyz )

   call ogderiv(ep, 0,1,2,0, x1,y1,z1,t, ez, wexyy ) ! add
   call ogderiv(ep, 0,1,0,2, x1,y1,z1,t, ez, wexzz )
   call ogderiv(ep, 0,0,1,2, x1,y1,z1,t, ez, weyzz )

   call ogderiv(ep, 0,4,0,0, x1,y1,z1,t, ez, wexxxx )
   call ogderiv(ep, 0,2,2,0, x1,y1,z1,t, ez, wexxyy )
   call ogderiv(ep, 0,2,0,2, x1,y1,z1,t, ez, wexxzz )
   call ogderiv(ep, 0,0,4,0, x1,y1,z1,t, ez, weyyyy )
   call ogderiv(ep, 0,0,2,2, x1,y1,z1,t, ez, weyyzz )
   call ogderiv(ep, 0,0,0,4, x1,y1,z1,t, ez, wezzzz )


   ueLap = uexx+ueyy+uezz
   veLap = vexx+veyy+vezz
   weLap = wexx+weyy+wezz

   curlEex = wey-vez
   curlEey = uez-wex
   curlEez = vex-uey
   nDotCurlEe=an1*curlEex+an2*curlEey+an3*curlEez
   nDotLapEe=an1*ueLap+an2*veLap+an3*weLap

   f(0)= f(0) - ( (curlEex- nDotCurlEe*an1)*(1./mu1-1./mu2) )
   f(1)= f(1) - ( (curlEey- nDotCurlEe*an2)*(1./mu1-1./mu2) )
   f(2)= f(2) - ( (curlEez- nDotCurlEe*an3)*(1./mu1-1./mu2) )

   f(3)= f(3) - ( ueLap*(1./epsmu1-1./epsmu2) + nDotLapEe*an1*(cem1-cem2) )
   f(4)= f(4) - ( veLap*(1./epsmu1-1./epsmu2) + nDotLapEe*an2*(cem1-cem2) )
   f(5)= f(5) - ( weLap*(1./epsmu1-1./epsmu2) + nDotLapEe*an3*(cem1-cem2) ) 


   ueLapSq = uexxxx+ueyyyy+uezzzz +2.*( uexxyy+uexxzz+ueyyzz )
   veLapSq = vexxxx+veyyyy+vezzzz +2.*( vexxyy+vexxzz+veyyzz )
   weLapSq = wexxxx+weyyyy+wezzzz +2.*( wexxyy+wexxzz+weyyzz )

   curlLapEex = wexxy-vexxz + weyyy-veyyz + weyzz-vezzz
   curlLapEey = uexxz-wexxx + ueyyz-wexyy + uezzz-wexzz
   curlLapEez = vexxx-uexxy + vexyy-ueyyy + vexzz-ueyzz
   nDotCurlLapEe=an1*curlLapEex+an2*curlLapEey+an3*curlLapEez
   nDotLapSqEe=an1*ueLapSq+an2*veLapSq+an3*weLapSq

 ! *e678*
   f(6) = f(6) - ( (curlLapEex- nDotCurlLapEe*an1)*(1./(mu1*epsmu1)-1./(mu2*epsmu2)) )
   f(7) = f(7) - ( (curlLapEey- nDotCurlLapEe*an2)*(1./(mu1*epsmu1)-1./(mu2*epsmu2)) )
   f(8) = f(8) - ( (curlLapEez- nDotCurlLapEe*an3)*(1./(mu1*epsmu1)-1./(mu2*epsmu2)) )

   f(9) = f(9) - ( ueLapSq*(1./epsmu1**2-1./epsmu2**2) + nDotLapSqEe*an1*(cem1/epsmu1-cem2/epsmu2) )
   f(10)= f(10)- ( veLapSq*(1./epsmu1**2-1./epsmu2**2) + nDotLapSqEe*an2*(cem1/epsmu1-cem2/epsmu2) )
   f(11)= f(11)- ( weLapSq*(1./epsmu1**2-1./epsmu2**2) + nDotLapSqEe*an3*(cem1/epsmu1-cem2/epsmu2) ) 

 end if

#endMacro


! =====================================================================================
! initialization step: assign first two ghost line by extrapolation
! NOTE: assign ghost points outside the ends
! =====================================================================================
#beginMacro extrapGhost()
  if( .true. )then
   beginGhostLoops3d()
     u1(i1-is1,i2-is2,i3-is3,ex)=extrap5(u1,i1,i2,i3,ex,is1,is2,is3)
     u1(i1-is1,i2-is2,i3-is3,ey)=extrap5(u1,i1,i2,i3,ey,is1,is2,is3)
     u1(i1-is1,i2-is2,i3-is3,ez)=extrap5(u1,i1,i2,i3,ez,is1,is2,is3)

     u2(j1-js1,j2-js2,j3-js3,ex)=extrap5(u2,j1,j2,j3,ex,js1,js2,js3)
     u2(j1-js1,j2-js2,j3-js3,ey)=extrap5(u2,j1,j2,j3,ey,js1,js2,js3)
     u2(j1-js1,j2-js2,j3-js3,ez)=extrap5(u2,j1,j2,j3,ez,js1,js2,js3)

     u1(i1-2*is1,i2-2*is2,i3-2*is3,ex)=extrap5(u1,i1-is1,i2-is2,i3-is3,ex,is1,is2,is3)
     u1(i1-2*is1,i2-2*is2,i3-2*is3,ey)=extrap5(u1,i1-is1,i2-is2,i3-is3,ey,is1,is2,is3)
     u1(i1-2*is1,i2-2*is2,i3-2*is3,ez)=extrap5(u1,i1-is1,i2-is2,i3-is3,ez,is1,is2,is3)

     u2(j1-2*js1,j2-2*js2,j3-2*js3,ex)=extrap5(u2,j1-js1,j2-js2,j3-js3,ex,js1,js2,js3)
     u2(j1-2*js1,j2-2*js2,j3-2*js3,ey)=extrap5(u2,j1-js1,j2-js2,j3-js3,ey,js1,js2,js3)
     u2(j1-2*js1,j2-2*js2,j3-2*js3,ez)=extrap5(u2,j1-js1,j2-js2,j3-js3,ez,js1,js2,js3)


     if( dispersive.ne.noDispersion )then
      do jv=0,numberOfPolarizationVectors1-1
        do n=0,nd-1
          pc = n + jv*nd 
          p1(i1  -is1,i2  -is2,i3-  is3,pc)=extrap5(p1,i1    ,i2    ,i3    ,pc,is1,is2,is3)
          p1(i1-2*is1,i2-2*is2,i3-2*is3,pc)=extrap5(p1,i1-is1,i2-is2,i3-is3,pc,is1,is2,is3)

        end do
      end do
      do jv=0,numberOfPolarizationVectors2-1
        do n=0,nd-1
          pc = n + jv*nd 
          p2(j1  -js1,j2  -js2,j3-  js3,pc)=extrap5(p2,j1    ,j2    ,j3    ,pc,js1,js2,js3)
          p2(j1-2*js1,j2-2*js2,j3-2*js3,pc)=extrap5(p2,j1-js1,j2-js2,j3-js3,pc,js1,js2,js3)
        end do
      end do
     end if

   endLoops3d()


  else if( .false. )then
    ! Extrap points outside interp points to lower order:
   beginGhostLoops3d()
    if( mask1(i1,i2,i3).gt.0 )then
     u1(i1-is1,i2-is2,i3-is3,ex)=extrap4(u1,i1,i2,i3,ex,is1,is2,is3)
     u1(i1-is1,i2-is2,i3-is3,ey)=extrap4(u1,i1,i2,i3,ey,is1,is2,is3)
     u1(i1-is1,i2-is2,i3-is3,ez)=extrap4(u1,i1,i2,i3,ez,is1,is2,is3)

     u2(j1-js1,j2-js2,j3-js3,ex)=extrap4(u2,j1,j2,j3,ex,js1,js2,js3)
     u2(j1-js1,j2-js2,j3-js3,ey)=extrap4(u2,j1,j2,j3,ey,js1,js2,js3)
     u2(j1-js1,j2-js2,j3-js3,ez)=extrap4(u2,j1,j2,j3,ez,js1,js2,js3)

     u1(i1-2*is1,i2-2*is2,i3-2*is3,ex)=extrap4(u1,i1-is1,i2-is2,i3-is3,ex,is1,is2,is3)
     u1(i1-2*is1,i2-2*is2,i3-2*is3,ey)=extrap4(u1,i1-is1,i2-is2,i3-is3,ey,is1,is2,is3)
     u1(i1-2*is1,i2-2*is2,i3-2*is3,ez)=extrap4(u1,i1-is1,i2-is2,i3-is3,ez,is1,is2,is3)

     u2(j1-2*js1,j2-2*js2,j3-2*js3,ex)=extrap4(u2,j1-js1,j2-js2,j3-js3,ex,js1,js2,js3)
     u2(j1-2*js1,j2-2*js2,j3-2*js3,ey)=extrap4(u2,j1-js1,j2-js2,j3-js3,ey,js1,js2,js3)
     u2(j1-2*js1,j2-2*js2,j3-2*js3,ez)=extrap4(u2,j1-js1,j2-js2,j3-js3,ez,js1,js2,js3)

     ! u1(i1-2*is1,i2-2*is2,i3-2*is3,ex)=u1(i1+2*is1,i2+2*is2,i3+2*is3,ex)
     ! u1(i1-2*is1,i2-2*is2,i3-2*is3,ey)=u1(i1+2*is1,i2+2*is2,i3+2*is3,ey)
     ! u1(i1-2*is1,i2-2*is2,i3-2*is3,ez)=u1(i1+2*is1,i2+2*is2,i3+2*is3,ez)

     ! u2(j1-2*js1,j2-2*js2,j3-2*js3,ex)=u2(j1+2*js1,j2+2*js2,j3+2*js3,ex)
     ! u2(j1-2*js1,j2-2*js2,j3-2*js3,ey)=u2(j1+2*js1,j2+2*js2,j3+2*js3,ey)
     ! u2(j1-2*js1,j2-2*js2,j3-2*js3,ez)=u2(j1+2*js1,j2+2*js2,j3+2*js3,ez)

    else 
     ! try this for testing stability -- symmetry --
     an1=rsxy1(i1,i2,i3,axis1,0)   ! normal (an1,an2)
     an2=rsxy1(i1,i2,i3,axis1,1)
     an3=rsxy1(i1,i2,i3,axis1,2)
     aNorm=max(epsx,sqrt(an1**2+an2**2+an3**2))
     an1=an1/aNorm
     an2=an2/aNorm
     an3=an3/aNorm


     vectorSymmetry(u1,i1,i2,i3,is1,is2,is3)
     vectorSymmetry(u2,j1,j2,j3,js1,js2,js3)

     vectorSymmetry(u1,i1,i2,i3,2*is1,2*is2,2*is3)
     vectorSymmetry(u2,j1,j2,j3,2*js1,2*js2,2*js3)

    end if
   endLoops3d()

  else if( .true. )then
    ! try this for testing stability -- symmetry --
   beginGhostLoops3d()

     an1=rsxy1(i1,i2,i3,axis1,0)   ! normal (an1,an2)
     an2=rsxy1(i1,i2,i3,axis1,1)
     an3=rsxy1(i1,i2,i3,axis1,2)
     aNorm=max(epsx,sqrt(an1**2+an2**2+an3**2))
     an1=an1/aNorm
     an2=an2/aNorm
     an3=an3/aNorm


     vectorSymmetry(u1,i1,i2,i3,is1,is2,is3)
     vectorSymmetry(u2,j1,j2,j3,js1,js2,js3)

     vectorSymmetry(u1,i1,i2,i3,2*is1,2*is2,2*is3)
     vectorSymmetry(u2,j1,j2,j3,2*js1,2*js2,2*js3)

     ! u1(i1-is1,i2-is2,i3-is3,ex)=u1(i1+is1,i2+is2,i3+is3,ex)
     ! u1(i1-is1,i2-is2,i3-is3,ey)=u1(i1+is1,i2+is2,i3+is3,ey)
     ! u1(i1-is1,i2-is2,i3-is3,ez)=u1(i1+is1,i2+is2,i3+is3,ez)

     ! u2(j1-js1,j2-js2,j3-js3,ex)=u2(j1+js1,j2+js2,j3+js3,ex)
     ! u2(j1-js1,j2-js2,j3-js3,ey)=u2(j1+js1,j2+js2,j3+js3,ey)
     ! u2(j1-js1,j2-js2,j3-js3,ez)=u2(j1+js1,j2+js2,j3+js3,ez)

     ! u1(i1-2*is1,i2-2*is2,i3-2*is3,ex)=u1(i1+2*is1,i2+2*is2,i3+2*is3,ex)
     ! u1(i1-2*is1,i2-2*is2,i3-2*is3,ey)=u1(i1+2*is1,i2+2*is2,i3+2*is3,ey)
     ! u1(i1-2*is1,i2-2*is2,i3-2*is3,ez)=u1(i1+2*is1,i2+2*is2,i3+2*is3,ez)

     ! u2(j1-2*js1,j2-2*js2,j3-2*js3,ex)=u2(j1+2*js1,j2+2*js2,j3+2*js3,ex)
     ! u2(j1-2*js1,j2-2*js2,j3-2*js3,ey)=u2(j1+2*js1,j2+2*js2,j3+2*js3,ey)
     ! u2(j1-2*js1,j2-2*js2,j3-2*js3,ez)=u2(j1+2*js1,j2+2*js2,j3+2*js3,ez)
   endLoops3d()

  else
    stop 3003
  end if

  if( .false. )then
   ! just copy values from ghost points for now -- this will be the true soln if eps1=eps2 and grids match
   beginLoops3d()
    u1(i1-is1,i2-is2,i3-is3,ex)=u2(j1+js1,j2+js2,j3+js3,ex)
    u1(i1-is1,i2-is2,i3-is3,ey)=u2(j1+js1,j2+js2,j3+js3,ey)
    u1(i1-is1,i2-is2,i3-is3,ez)=u2(j1+js1,j2+js2,j3+js3,ez) 
    u2(j1-js1,j2-js2,j3-js3,ex)=u1(i1+is1,i2+is2,i3+is3,ex)
    u2(j1-js1,j2-js2,j3-js3,ey)=u1(i1+is1,i2+is2,i3+is3,ey)
    u2(j1-js1,j2-js2,j3-js3,ez)=u1(i1+is1,i2+is2,i3+is3,ez)
   endLoops3d()
  end if

 if( .false. )then
   beginGhostLoops3d()
    write(debugFile,'(" -->EXTRAP v1(",i2,":",i2,",",i2,",",i2,") =",3f9.4)') i1-1,i1+1,i2,i3,u1(i1-1,i2,i3,ey),u1(i1,i2,i3,ey),u1(i1+1,i2,i3,ey)
    ! '
   endLoops3d()
 end if
#endMacro


! ==========================================================================================
!         Fill in 6 of the interface equations.
! 
!  The equations come in sets of 6 (composed of 2 sets of 3) that are similar in structure so we can fill in
!  6 at a time.
! Input:
!  e0,e1,e2,e3,e4,e5 : equation numbers (0,1,2,3,4,5) or (6,7,8,9,10,11)
! ==========================================================================================
#beginMacro fillEquations3dOrder4(e0,e1,e2,e3,e4,e5)
 ! Equation 0: 
 ! (u.x+v.y+w.z)*an1 + ( w1y-v1z - nDotCurlE1*an1)/mu1
 a12(e0,0) = ( c1m1x*an1 + (           - (c1m1z*an2-c1m1y*an3)*an1 )/mu1 ) ! coeff of u1(-1)
 a12(e0,1) = ( c1m1y*an1 + (    -c1m1z - (c1m1x*an3-c1m1z*an1)*an1 )/mu1 ) ! coeff of v1(-1)
 a12(e0,2) = ( c1m1z*an1 + ( c1m1y     - (c1m1y*an1-c1m1x*an2)*an1 )/mu1 ) ! coeff of w1(-1)

 a12(e0,3) = ( c1m2x*an1 + (           - (c1m2z*an2-c1m2y*an3)*an1 )/mu1 ) ! coeff of u1(-1)
 a12(e0,4) = ( c1m2y*an1 + (    -c1m2z - (c1m2x*an3-c1m2z*an1)*an1 )/mu1 ) ! coeff of v1(-1)
 a12(e0,5) = ( c1m2z*an1 + ( c1m2y     - (c1m2y*an1-c1m2x*an2)*an1 )/mu1 ) ! coeff of w1(-1)

 a12(e0,6) =-( c2m1x*an1 + (           - (c2m1z*an2-c2m1y*an3)*an1 )/mu2 ) ! coeff of u2(-1)
 a12(e0,7) =-( c2m1y*an1 + (    -c2m1z - (c2m1x*an3-c2m1z*an1)*an1 )/mu2 ) ! coeff of v2(-1)
 a12(e0,8) =-( c2m1z*an1 + ( c2m1y     - (c2m1y*an1-c2m1x*an2)*an1 )/mu2 ) ! coeff of w2(-1)

 a12(e0,9) =-( c2m2x*an1 + (           - (c2m2z*an2-c2m2y*an3)*an1 )/mu2 ) ! coeff of u2(-1)
 a12(e0,10)=-( c2m2y*an1 + (    -c2m2z - (c2m2x*an3-c2m2z*an1)*an1 )/mu2 ) ! coeff of v2(-1)
 a12(e0,11)=-( c2m2z*an1 + ( c2m2y     - (c2m2y*an1-c2m2x*an2)*an1 )/mu2 ) ! coeff of w2(-1)

 ! Bug fixed - July 6, 2019 wdh : Changed "+" to minus "-" in 8 lines below 

 ! Equation 1:
 ! (u.x+v.y+w.z)*an2 + ( u1z-w1x - nDotCurlE1*an2)/mu1
 a12(e1,0) = ( c1m1x*an2 + ( c1m1z     - (c1m1z*an2-c1m1y*an3)*an2 )/mu1 ) ! coeff of u1(-1)
!a12(e1,1) = ( c1m1y*an2 + (           - (c1m1x*an3+c1m1z*an1)*an2 )/mu1 ) ! coeff of v1(-1)
 a12(e1,1) = ( c1m1y*an2 + (           - (c1m1x*an3-c1m1z*an1)*an2 )/mu1 ) ! coeff of v1(-1)
 a12(e1,2) = ( c1m1z*an2 + (    -c1m1x - (c1m1y*an1-c1m1x*an2)*an2 )/mu1 ) ! coeff of w1(-1)

 a12(e1,3) = ( c1m2x*an2 + ( c1m2z     - (c1m2z*an2-c1m2y*an3)*an2 )/mu1 ) ! coeff of u1(-1)
!a12(e1,4) = ( c1m2y*an2 + (           - (c1m2x*an3+c1m2z*an1)*an2 )/mu1 ) ! coeff of v1(-1)
 a12(e1,4) = ( c1m2y*an2 + (           - (c1m2x*an3-c1m2z*an1)*an2 )/mu1 ) ! coeff of v1(-1)
 a12(e1,5) = ( c1m2z*an2 + (    -c1m2x - (c1m2y*an1-c1m2x*an2)*an2 )/mu1 ) ! coeff of w1(-1)

 a12(e1,6) =-( c2m1x*an2 + ( c2m1z     - (c2m1z*an2-c2m1y*an3)*an2 )/mu2 ) ! coeff of u2(-1)
!a12(e1,7) =-( c2m1y*an2 + (           - (c2m1x*an3+c2m1z*an1)*an2 )/mu2 ) ! coeff of v2(-1)
 a12(e1,7) =-( c2m1y*an2 + (           - (c2m1x*an3-c2m1z*an1)*an2 )/mu2 ) ! coeff of v2(-1)
 a12(e1,8) =-( c2m1z*an2 + (    -c2m1x - (c2m1y*an1-c2m1x*an2)*an2 )/mu2 ) ! coeff of w2(-1)

 a12(e1,9) =-( c2m2x*an2 + ( c2m2z     - (c2m2z*an2-c2m2y*an3)*an2 )/mu2 ) ! coeff of u2(-1)
!a12(e1,10)=-( c2m2y*an2 + (           - (c2m2x*an3+c2m2z*an1)*an2 )/mu2 ) ! coeff of v2(-1)
 a12(e1,10)=-( c2m2y*an2 + (           - (c2m2x*an3-c2m2z*an1)*an2 )/mu2 ) ! coeff of v2(-1)
 a12(e1,11)=-( c2m2z*an2 + (    -c2m2x - (c2m2y*an1-c2m2x*an2)*an2 )/mu2 ) ! coeff of w2(-1)

 ! Equation 2:
 ! (u.x+v.y+w.z)*an3 + ( v1x-u1y - nDotCurlE1*an3)/mu1
 a12(e2,0) = ( c1m1x*an3 + (    -c1m1y - (c1m1z*an2-c1m1y*an3)*an3 )/mu1 ) ! coeff of u1(-1)
!a12(e2,1) = ( c1m1y*an3 + ( c1m1x     - (c1m1x*an3+c1m1z*an1)*an3 )/mu1 ) ! coeff of v1(-1)
 a12(e2,1) = ( c1m1y*an3 + ( c1m1x     - (c1m1x*an3-c1m1z*an1)*an3 )/mu1 ) ! coeff of v1(-1)
 a12(e2,2) = ( c1m1z*an3 + (           - (c1m1y*an1-c1m1x*an2)*an3 )/mu1 ) ! coeff of w1(-1)

 a12(e2,3) = ( c1m2x*an3 + (    -c1m2y - (c1m2z*an2-c1m2y*an3)*an3 )/mu1 ) ! coeff of u1(-1)
 a12(e2,4) = ( c1m2y*an3 + ( c1m2x     - (c1m2x*an3-c1m2z*an1)*an3 )/mu1 ) ! coeff of v1(-1)
 a12(e2,5) = ( c1m2z*an3 + (           - (c1m2y*an1-c1m2x*an2)*an3 )/mu1 ) ! coeff of w1(-1)

 a12(e2,6) =-( c2m1x*an3 + (    -c2m1y - (c2m1z*an2-c2m1y*an3)*an3 )/mu2 ) ! coeff of u2(-1)
 a12(e2,7) =-( c2m1y*an3 + ( c2m1x     - (c2m1x*an3-c2m1z*an1)*an3 )/mu2 ) ! coeff of v2(-1)
 a12(e2,8) =-( c2m1z*an3 + (           - (c2m1y*an1-c2m1x*an2)*an3 )/mu2 ) ! coeff of w2(-1)

 a12(e2,9) =-( c2m2x*an3 + (    -c2m2y - (c2m2z*an2-c2m2y*an3)*an3 )/mu2 ) ! coeff of u2(-1)
 a12(e2,10)=-( c2m2y*an3 + ( c2m2x     - (c2m2x*an3-c2m2z*an1)*an3 )/mu2 ) ! coeff of v2(-1)
 a12(e2,11)=-( c2m2z*an3 + (           - (c2m2y*an1-c2m2x*an2)*an3 )/mu2 ) ! coeff of w2(-1)

 ! Equation 3:
 !  u1Lap/(epsmu1) + cem1*( an1*u1Lap + an2*v1Lap + an3*w1Lap )*an1
 a12(e3,0) = ( clap1m1/(epsmu1) + cem1*( an1*clap1m1                         )*an1 ) ! coeff of u1(-1)
 a12(e3,1) = (                    cem1*(             an2*clap1m1             )*an1 )
 a12(e3,2) = (                    cem1*(                         an3*clap1m1 )*an1 )

 a12(e3,3) = ( clap1m2/(epsmu1) + cem1*( an1*clap1m2                         )*an1 ) ! coeff of u1(-2)
 a12(e3,4) = (                    cem1*(             an2*clap1m2             )*an1 )
 a12(e3,5) = (                    cem1*(                         an3*clap1m2 )*an1 )

 a12(e3,6) =-( clap2m1/(epsmu2) + cem2*( an1*clap2m1                         )*an1 ) ! coeff of u2(-1)
 a12(e3,7) =-(                    cem2*(             an2*clap2m1             )*an1 )
 a12(e3,8) =-(                    cem2*(                         an3*clap2m1 )*an1 )

 a12(e3,9) =-( clap2m2/(epsmu2) + cem2*( an1*clap2m2                         )*an1 ) ! coeff of u2(-2)
 a12(e3,10)=-(                    cem2*(             an2*clap2m2             )*an1 )
 a12(e3,11)=-(                    cem2*(                         an3*clap2m2 )*an1 )

 ! Equation 4:
 !  v1Lap/(epsmu1) + cem1*( an1*u1Lap + an2*v1Lap + an3*w1Lap )*an2
 a12(e4,0) = (                    cem1*( an1*clap1m1                         )*an2 ) ! coeff of u1(-1)
 a12(e4,1) = ( clap1m1/(epsmu1) + cem1*(             an2*clap1m1             )*an2 )
 a12(e4,2) = (                    cem1*(                         an3*clap1m1 )*an2 )

 a12(e4,3) = (                    cem1*( an1*clap1m2                         )*an2 ) ! coeff of u1(-2)
 a12(e4,4) = ( clap1m2/(epsmu1) + cem1*(             an2*clap1m2             )*an2 )
 a12(e4,5) = (                    cem1*(                         an3*clap1m2 )*an2 )

 a12(e4,6) =-(                    cem2*( an1*clap2m1                         )*an2 ) ! coeff of u2(-1)
 a12(e4,7) =-( clap2m1/(epsmu2) + cem2*(             an2*clap2m1             )*an2 )
 a12(e4,8) =-(                    cem2*(                         an3*clap2m1 )*an2 )

 a12(e4,9) =-(                    cem2*( an1*clap2m2                         )*an2 ) ! coeff of u2(-2)
 a12(e4,10)=-( clap2m2/(epsmu2) + cem2*(             an2*clap2m2             )*an2 )
 a12(e4,11)=-(                    cem2*(                         an3*clap2m2 )*an2 )

 ! Equation 5:
 !  w1Lap/(epsmu1) + cem1*( an1*u1Lap + an2*v1Lap + an3*w1Lap )*an3
 a12(e5,0) = (                    cem1*( an1*clap1m1                         )*an3 ) ! coeff of u1(-1)
 a12(e5,1) = (                    cem1*(             an2*clap1m1             )*an3 )
 a12(e5,2) = ( clap1m1/(epsmu1) + cem1*(                         an3*clap1m1 )*an3 )

 a12(e5,3) = (                    cem1*( an1*clap1m2                         )*an3 ) ! coeff of u1(-2)
 a12(e5,4) = (                    cem1*(             an2*clap1m2             )*an3 )
 a12(e5,5) = ( clap1m2/(epsmu1) + cem1*(                         an3*clap1m2 )*an3 )

 a12(e5,6) =-(                    cem2*( an1*clap2m1                         )*an3 ) ! coeff of u2(-1)
 a12(e5,7) =-(                    cem2*(             an2*clap2m1             )*an3 )
 a12(e5,8) =-( clap2m1/(epsmu2) + cem2*(                         an3*clap2m1 )*an3 )

 a12(e5,9) =-(                    cem2*( an1*clap2m2                         )*an3 ) ! coeff of u2(-2)
 a12(e5,10)=-(                    cem2*(             an2*clap2m2             )*an3 )
 a12(e5,11)=-( clap2m2/(epsmu2) + cem2*(                         an3*clap2m2 )*an3 )
#endMacro


! =================================================================================
! Assign a given value to a set of variables
! =================================================================================
#beginMacro setValues10(value,u1,u2,u3,u4,u5,u6,u7,u8,u9,u10)
 u1=value
 u2=value
 u3=value
 u4=value
 u5=value
 u6=value
 u7=value
 u8=value
 u9=value
 u10=value
#endMacro
#beginMacro setValues9(value,u1,u2,u3,u4,u5,u6,u7,u8,u9)
 u1=value
 u2=value
 u3=value
 u4=value
 u5=value
 u6=value
 u7=value
 u8=value
 u9=value
#endMacro
#beginMacro setValues7(value,u1,u2,u3,u4,u5,u6,u7)
 u1=value
 u2=value
 u3=value
 u4=value
 u5=value
 u6=value
 u7=value
#endMacro
#beginMacro setValues6(value,u1,u2,u3,u4,u5,u6)
 u1=value
 u2=value
 u3=value
 u4=value
 u5=value
 u6=value
#endMacro


! --------------------------------------------------------------------------
! Macro: Assign interface ghost values, DIM=3, ORDER=4, GRID=Curvilinear
!     **OLD VERSION THAT EXTRAPOLATES SECOND GHOST LINE ****
! --------------------------------------------------------------------------
#beginMacro assignInterfaceGhost3dOrder4Old()
  err=0.
  err2=0.
  count=0
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

    cem1=(1.-1./eps1)/mu1
    cem2=(1.-1./eps2)/mu2
    ! evalInterfaceDerivatives3d
    evalInterfaceDerivatives3d()
    eval3dJumpOrder4()

    if( debug.gt.4 )then
     write(debugFile,'(/," --> 3d-order4-curv: i1,i2,i3=",3i4," an1,an2,an3=",3e11.3)') i1,i2,i3,an1,an2,an3
     write(debugFile,'(" --> 3d-order4-curv: i1,i2,i3=",3i4," f(start)=",12e10.2)') i1,i2,i3,f(0),f(1),f(2),f(3),f(4),f(5),f(6),f(7),f(8),f(9),f(10),f(11)
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


    c1m1x = -is*8.*rsxy1(i1,i2,i3,axis1,0)/(12.*dr1(axis1))    ! coeff of u1(-1) from D.x
    c1m1y = -is*8.*rsxy1(i1,i2,i3,axis1,1)/(12.*dr1(axis1))    ! coeff of u1(-1) from D.y 
    c1m1z = -is*8.*rsxy1(i1,i2,i3,axis1,2)/(12.*dr1(axis1))    ! coeff of u1(-1) from D.z

    c1m2x =  is   *rsxy1(i1,i2,i3,axis1,0)/(12.*dr1(axis1))    ! coeff of u1(-2) from D.x
    c1m2y =  is   *rsxy1(i1,i2,i3,axis1,1)/(12.*dr1(axis1))    ! coeff of u1(-2) from D.y 
    c1m2z =  is   *rsxy1(i1,i2,i3,axis1,2)/(12.*dr1(axis1))    ! coeff of u1(-2) from D.z

    c2m1x = -js*8.*rsxy2(j1,j2,j3,axis2,0)/(12.*dr2(axis2))    ! coeff of u2(-1) from D.x
    c2m1y = -js*8.*rsxy2(j1,j2,j3,axis2,1)/(12.*dr2(axis2))
    c2m1z = -js*8.*rsxy2(j1,j2,j3,axis2,2)/(12.*dr2(axis2))

    c2m2x =  js   *rsxy2(j1,j2,j3,axis2,0)/(12.*dr2(axis2))
    c2m2y =  js   *rsxy2(j1,j2,j3,axis2,1)/(12.*dr2(axis2))
    c2m2z =  js   *rsxy2(j1,j2,j3,axis2,2)/(12.*dr2(axis2))

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

    ! clap1m1 : coeff of u(-1) from lap1 = u1.xx + u1.yy + u1.zz
    ! clap1m2 : coeff of u(-2) from lap1 = u1.xx + u1.yy + u1.zz

    ! here are the macros from deriv.maple (file=derivMacros.h)
    ! defineMacro lapCoeff4a(is,dr,ds) ( (-2/3.*rxx*is-2/3.*ryy*is)/dr+(4/3.*rx**2+4/3.*ry**2)/dr**2 )
    ! defineMacro lapCoeff4b(is,dr,ds) ( (1/12.*rxx*is+1/12.*ryy*is)/dr+(-1/12.*rx**2-1/12.*ry**2)/dr**2 )
    
    clap1m1=4./3.*(rsxy1(i1,i2,i3,axis1,0)**2+rsxy1(i1,i2,i3,axis1,1)**2+rsxy1(i1,i2,i3,axis1,2)**2)/(dr1(axis1)**2) \
              -is*2./3.*(rxx1(axis1,0,0)+rxx1(axis1,1,1)+rxx1(axis1,2,2))/(2.*dr1(axis1))
    clap1m2=-1./12.*(rsxy1(i1,i2,i3,axis1,0)**2+rsxy1(i1,i2,i3,axis1,1)**2+rsxy1(i1,i2,i3,axis1,2)**2)/(dr1(axis1)**2) \
              +is*1./12.*(rxx1(axis1,0,0)+rxx1(axis1,1,1)+rxx1(axis1,2,2))/(2.*dr1(axis1)) 

    clap2m1=4/3.*(rsxy2(j1,j2,j3,axis2,0)**2+rsxy2(j1,j2,j3,axis2,1)**2+rsxy2(j1,j2,j3,axis2,2)**2)/(dr2(axis2)**2) \
              -js*2./3.*(rxx2(axis2,0,0)+rxx2(axis2,1,1)+rxx2(axis2,2,2))/(2.*dr2(axis2)) 
    clap2m2=-1./12.*(rsxy2(j1,j2,j3,axis2,0)**2+rsxy2(j1,j2,j3,axis2,1)**2+rsxy2(j1,j2,j3,axis2,2)**2)/(dr2(axis2)**2) \
              +js*1./12.*(rxx2(axis2,0,0)+rxx2(axis2,1,1)+rxx2(axis2,2,2))/(2.*dr2(axis2))

    ! cdivE1 =  u.c1x + v.c1y + w.c1z
    ! nDotCurlE1 = (w1y-v1z)*an1 + (u1z-w1x)*an2 + (v1x-u1y)*an3

    ! 12 Unknowns:
    !   u1(-1), v1(-1), w1(-1), u1(-2), v1(-2), w1(-2), 
    !   u2(-1), v2(-1), w2(-1), u2(-2), v2(-2), w2(-2)  
    ! 12 Equations: 
    !    a12(eqn,unknown) 

    ! Equation 0: 
    ! (u.x+v.y+w.z)*an1 + ( w1y-v1z - nDotCurlE1*an1)/mu1
    a12(0,0) = ( c1m1x*an1 + (           - (c1m1z*an2-c1m1y*an3)*an1 )/mu1 ) ! coeff of u1(-1)
    a12(0,1) = ( c1m1y*an1 + (    -c1m1z - (c1m1x*an3-c1m1z*an1)*an1 )/mu1 ) ! coeff of v1(-1)
    a12(0,2) = ( c1m1z*an1 + ( c1m1y     - (c1m1y*an1-c1m1x*an2)*an1 )/mu1 ) ! coeff of w1(-1)

    a12(0,3) = ( c1m2x*an1 + (           - (c1m2z*an2-c1m2y*an3)*an1 )/mu1 ) ! coeff of u1(-1)
    a12(0,4) = ( c1m2y*an1 + (    -c1m2z - (c1m2x*an3-c1m2z*an1)*an1 )/mu1 ) ! coeff of v1(-1)
    a12(0,5) = ( c1m2z*an1 + ( c1m2y     - (c1m2y*an1-c1m2x*an2)*an1 )/mu1 ) ! coeff of w1(-1)

    a12(0,6) =-( c2m1x*an1 + (           - (c2m1z*an2-c2m1y*an3)*an1 )/mu2 ) ! coeff of u2(-1)
    a12(0,7) =-( c2m1y*an1 + (    -c2m1z - (c2m1x*an3-c2m1z*an1)*an1 )/mu2 ) ! coeff of v2(-1)
    a12(0,8) =-( c2m1z*an1 + ( c2m1y     - (c2m1y*an1-c2m1x*an2)*an1 )/mu2 ) ! coeff of w2(-1)

    a12(0,9) =-( c2m2x*an1 + (           - (c2m2z*an2-c2m2y*an3)*an1 )/mu2 ) ! coeff of u2(-1)
    a12(0,10)=-( c2m2y*an1 + (    -c2m2z - (c2m2x*an3-c2m2z*an1)*an1 )/mu2 ) ! coeff of v2(-1)
    a12(0,11)=-( c2m2z*an1 + ( c2m2y     - (c2m2y*an1-c2m2x*an2)*an1 )/mu2 ) ! coeff of w2(-1)

    ! Equation 1:
    ! (u.x+v.y+w.z)*an2 + ( u1z-w1x - nDotCurlE1*an2)/mu1
    a12(1,0) = ( c1m1x*an2 + ( c1m1z     - (c1m1z*an2-c1m1y*an3)*an2 )/mu1 ) ! coeff of u1(-1)
    a12(1,1) = ( c1m1y*an2 + (           - (c1m1x*an3+c1m1z*an1)*an2 )/mu1 ) ! coeff of v1(-1)
    a12(1,2) = ( c1m1z*an2 + (    -c1m1x - (c1m1y*an1-c1m1x*an2)*an2 )/mu1 ) ! coeff of w1(-1)

    a12(1,3) = ( c1m2x*an2 + ( c1m2z     - (c1m2z*an2-c1m2y*an3)*an2 )/mu1 ) ! coeff of u1(-1)
    a12(1,4) = ( c1m2y*an2 + (           - (c1m2x*an3+c1m2z*an1)*an2 )/mu1 ) ! coeff of v1(-1)
    a12(1,5) = ( c1m2z*an2 + (    -c1m2x - (c1m2y*an1-c1m2x*an2)*an2 )/mu1 ) ! coeff of w1(-1)

    a12(1,6) =-( c2m1x*an2 + ( c2m1z     - (c2m1z*an2-c2m1y*an3)*an2 )/mu2 ) ! coeff of u2(-1)
    a12(1,7) =-( c2m1y*an2 + (           - (c2m1x*an3+c2m1z*an1)*an2 )/mu2 ) ! coeff of v2(-1)
    a12(1,8) =-( c2m1z*an2 + (    -c2m1x - (c2m1y*an1-c2m1x*an2)*an2 )/mu2 ) ! coeff of w2(-1)

    a12(1,9) =-( c2m2x*an2 + ( c2m2z     - (c2m2z*an2-c2m2y*an3)*an2 )/mu2 ) ! coeff of u2(-1)
    a12(1,10)=-( c2m2y*an2 + (           - (c2m2x*an3+c2m2z*an1)*an2 )/mu2 ) ! coeff of v2(-1)
    a12(1,11)=-( c2m2z*an2 + (    -c2m2x - (c2m2y*an1-c2m2x*an2)*an2 )/mu2 ) ! coeff of w2(-1)

    ! Equation 2:
    ! (u.x+v.y+w.z)*an3 + ( v1x-u1y - nDotCurlE1*an2)/mu1
    a12(2,0) = ( c1m1x*an3 + (    -c1m1y - (c1m1z*an2-c1m1y*an3)*an3 )/mu1 ) ! coeff of u1(-1)
    a12(2,1) = ( c1m1y*an3 + ( c1m1x     - (c1m1x*an3+c1m1z*an1)*an3 )/mu1 ) ! coeff of v1(-1)
    a12(2,2) = ( c1m1z*an3 + (           - (c1m1y*an1-c1m1x*an2)*an3 )/mu1 ) ! coeff of w1(-1)

    a12(2,3) = ( c1m2x*an3 + (    -c1m2y - (c1m2z*an2-c1m2y*an3)*an3 )/mu1 ) ! coeff of u1(-1)
    a12(2,4) = ( c1m2y*an3 + ( c1m2x     - (c1m2x*an3+c1m2z*an1)*an3 )/mu1 ) ! coeff of v1(-1)
    a12(2,5) = ( c1m2z*an3 + (           - (c1m2y*an1-c1m2x*an2)*an3 )/mu1 ) ! coeff of w1(-1)

    a12(2,6) =-( c2m1x*an3 + (    -c2m1y - (c2m1z*an2-c2m1y*an3)*an3 )/mu2 ) ! coeff of u2(-1)
    a12(2,7) =-( c2m1y*an3 + ( c2m1x     - (c2m1x*an3+c2m1z*an1)*an3 )/mu2 ) ! coeff of v2(-1)
    a12(2,8) =-( c2m1z*an3 + (           - (c2m1y*an1-c2m1x*an2)*an3 )/mu2 ) ! coeff of w2(-1)

    a12(2,9) =-( c2m2x*an3 + (    -c2m2y - (c2m2z*an2-c2m2y*an3)*an3 )/mu2 ) ! coeff of u2(-1)
    a12(2,10)=-( c2m2y*an3 + ( c2m2x     - (c2m2x*an3+c2m2z*an1)*an3 )/mu2 ) ! coeff of v2(-1)
    a12(2,11)=-( c2m2z*an3 + (           - (c2m2y*an1-c2m2x*an2)*an3 )/mu2 ) ! coeff of w2(-1)

    ! Equation 3:
    !  u1Lap/(epsmu1) + cem1*( an1*u1Lap + an2*v1Lap + an3*w1Lap )*an1
    a12(3,0) = ( clap1m1/(epsmu1) + cem1*( an1*clap1m1                         )*an1 ) ! coeff of u1(-1)
    a12(3,1) = (                    cem1*(             an2*clap1m1             )*an1 )
    a12(3,2) = (                    cem1*(                         an3*clap1m1 )*an1 )

    a12(3,3) = ( clap1m2/(epsmu1) + cem1*( an1*clap1m2                         )*an1 ) ! coeff of u1(-2)
    a12(3,4) = (                    cem1*(             an2*clap1m2             )*an1 )
    a12(3,5) = (                    cem1*(                         an3*clap1m2 )*an1 )

    a12(3,6) =-( clap2m1/(epsmu2) + cem2*( an1*clap2m1                         )*an1 ) ! coeff of u2(-1)
    a12(3,7) =-(                    cem2*(             an2*clap2m1             )*an1 )
    a12(3,8) =-(                    cem2*(                         an3*clap2m1 )*an1 )

    a12(3,9) =-( clap2m2/(epsmu2) + cem2*( an1*clap2m2                         )*an1 ) ! coeff of u2(-2)
    a12(3,10)=-(                    cem2*(             an2*clap2m2             )*an1 )
    a12(3,11)=-(                    cem2*(                         an3*clap2m2 )*an1 )

    ! Equation 4:
    !  v1Lap/(epsmu1) + cem1*( an1*u1Lap + an2*v1Lap + an3*w1Lap )*an2
    a12(4,0) = (                    cem1*( an1*clap1m1                         )*an2 ) ! coeff of u1(-1)
    a12(4,1) = ( clap1m1/(epsmu1) + cem1*(             an2*clap1m1             )*an2 )
    a12(4,2) = (                    cem1*(                         an3*clap1m1 )*an2 )

    a12(4,3) = (                    cem1*( an1*clap1m2                         )*an2 ) ! coeff of u1(-2)
    a12(4,4) = ( clap1m2/(epsmu1) + cem1*(             an2*clap1m2             )*an2 )
    a12(4,5) = (                    cem1*(                         an3*clap1m2 )*an2 )

    a12(4,6) =-(                    cem2*( an1*clap2m1                         )*an2 ) ! coeff of u2(-1)
    a12(4,7) =-( clap2m1/(epsmu2) + cem2*(             an2*clap2m1             )*an2 )
    a12(4,8) =-(                    cem2*(                         an3*clap2m1 )*an2 )

    a12(4,9) =-(                    cem2*( an1*clap2m2                         )*an2 ) ! coeff of u2(-2)
    a12(4,10)=-( clap2m2/(epsmu2) + cem2*(             an2*clap2m2             )*an2 )
    a12(4,11)=-(                    cem2*(                         an3*clap2m2 )*an2 )

    ! Equation 5:
    !  w1Lap/(epsmu1) + cem1*( an1*u1Lap + an2*v1Lap + an3*w1Lap )*an3
    a12(5,0) = (                    cem1*( an1*clap1m1                         )*an3 ) ! coeff of u1(-1)
    a12(5,1) = (                    cem1*(             an2*clap1m1             )*an3 )
    a12(5,2) = ( clap1m1/(epsmu1) + cem1*(                         an3*clap1m1 )*an3 )

    a12(5,3) = (                    cem1*( an1*clap1m2                         )*an3 ) ! coeff of u1(-2)
    a12(5,4) = (                    cem1*(             an2*clap1m2             )*an3 )
    a12(5,5) = ( clap1m2/(epsmu1) + cem1*(                         an3*clap1m2 )*an3 )

    a12(5,6) =-(                    cem2*( an1*clap2m1                         )*an3 ) ! coeff of u2(-1)
    a12(5,7) =-(                    cem2*(             an2*clap2m1             )*an3 )
    a12(5,8) =-( clap2m1/(epsmu2) + cem2*(                         an3*clap2m1 )*an3 )

    a12(5,9) =-(                    cem2*( an1*clap2m2                         )*an3 ) ! coeff of u2(-2)
    a12(5,10)=-(                    cem2*(             an2*clap2m2             )*an3 )
    a12(5,11)=-( clap2m2/(epsmu2) + cem2*(                         an3*clap2m2 )*an3 )

    ! Equation 6..11 :  extrapolate 2nd ghost point 
    cex1=1.
    cex2=-5. ! ** fix me ** orderOfExtrapolation for 2nd ghost point 
    do ii=6,11
      do jj=0,11
        a12(ii,jj)=0.
      end do
    end do
    a12(6,0)  = cex2   ! u1(-1)
    a12(6,3)  = cex1   ! u1(-2)

    a12(7,1)  = cex2   ! v1(-1)
    a12(7,4)  = cex1   ! v1(-2)

    a12(8,2)  = cex2   ! w1(-1)
    a12(8,5)  = cex1   ! w1(-2)

    a12(9,6)  = cex2   ! u2(-1)
    a12(9,9)  = cex1   ! u2(-2)

    a12(10,7) = cex2   ! v2(-1)
    a12(10,10)= cex1   ! v2(-2)

    a12(11,8) = cex2   ! w2(-1)
    a12(11,11)= cex1   ! w2(-2)

    ! fill in the current values for the unknowns: 
    q(0) = u1(i1-is1,i2-is2,i3-is3,ex)
    q(1) = u1(i1-is1,i2-is2,i3-is3,ey)
    q(2) = u1(i1-is1,i2-is2,i3-is3,ez)
    q(3) = u1(i1-2*is1,i2-2*is2,i3-2*is3,ex)
    q(4) = u1(i1-2*is1,i2-2*is2,i3-2*is3,ey)
    q(5) = u1(i1-2*is1,i2-2*is2,i3-2*is3,ez)

    q(6) = u2(j1-js1,j2-js2,j3-js3,ex)
    q(7) = u2(j1-js1,j2-js2,j3-js3,ey)
    q(8) = u2(j1-js1,j2-js2,j3-js3,ez)
    q(9) = u2(j1-2*js1,j2-2*js2,j3-2*js3,ex)
    q(10)= u2(j1-2*js1,j2-2*js2,j3-2*js3,ey)
    q(11)= u2(j1-2*js1,j2-2*js2,j3-2*js3,ez)

    ! subtract off the contributions from the wrong values at the ghost points:
    numberOfEquations=12
    do n=0,numberOfEquations-1
      f(n) = (a12(n,0)*q(0)+a12(n,1)*q(1)+a12(n,2)*q(2)+a12(n,3)*q(3)+a12(n,4)*q(4)+a12(n,5)*q(5)+\
              a12(n,6)*q(6)+a12(n,7)*q(7)+a12(n,8)*q(8)+a12(n,9)*q(9)+a12(n,10)*q(10)+a12(n,11)*q(11) ) - f(n)
    end do
    ! write(debugFile,'(" --> 3d:order2-c: f(subtract)=",6f8.3)') f(0),f(1),f(2),f(3),f(4),f(5)
    ! solve A Q = F
    ! factor the matrix
    call dgeco( a12(0,0), numberOfEquations, numberOfEquations, ipvt(0),rcond,work(0))
    ! solve
     ! write(debugFile,'(" --> 3d:order2-c: rcond=",e10.2)') rcond
    job=0
    call dgesl( a12(0,0), numberOfEquations, numberOfEquations, ipvt(0), f(0), job)
    ! write(debugFile,'(" --> 3d:order2-c: f(solve)=",6f8.3)') f(0),f(1),f(2),f(3),f(4),f(5)
    ! write(debugFile,'(" --> 3d:order2-c:        q=",6f8.3)') q(0),q(1),q(2),q(3),q(4),q(5)

    ! fill in the answer:
    ! fill in the answer:
    if( useJacobiUpdate.eq.0 )then
      u1(i1-is1,i2-is2,i3-is3,ex)      =f(0 )*omega+(1.-omega)*q(0)
      u1(i1-is1,i2-is2,i3-is3,ey)      =f(1 )*omega+(1.-omega)*q(1)
      u1(i1-is1,i2-is2,i3-is3,ez)      =f(2 )*omega+(1.-omega)*q(2)
      u1(i1-2*is1,i2-2*is2,i3-2*is3,ex)=f(3 )*omega+(1.-omega)*q(3)
      u1(i1-2*is1,i2-2*is2,i3-2*is3,ey)=f(4 )*omega+(1.-omega)*q(4)
      u1(i1-2*is1,i2-2*is2,i3-2*is3,ez)=f(5 )*omega+(1.-omega)*q(5)
  
      u2(j1-js1,j2-js2,j3-js3,ex)      =f(6 )*omega+(1.-omega)*q(6)
      u2(j1-js1,j2-js2,j3-js3,ey)      =f(7 )*omega+(1.-omega)*q(7)
      u2(j1-js1,j2-js2,j3-js3,ez)      =f(8 )*omega+(1.-omega)*q(8)
      u2(j1-2*js1,j2-2*js2,j3-2*js3,ex)=f(9 )*omega+(1.-omega)*q(9)
      u2(j1-2*js1,j2-2*js2,j3-2*js3,ey)=f(10)*omega+(1.-omega)*q(10)
      u2(j1-2*js1,j2-2*js2,j3-2*js3,ez)=f(11)*omega+(1.-omega)*q(11)
    else
      ! Jacobi-update -- save results in the work space

      wk1(i1-is1,i2-is2,i3-is3,ex)      =f(0 )*omega+(1.-omega)*q(0)
      wk1(i1-is1,i2-is2,i3-is3,ey)      =f(1 )*omega+(1.-omega)*q(1)
      wk1(i1-is1,i2-is2,i3-is3,ez)      =f(2 )*omega+(1.-omega)*q(2)
      wk1(i1-2*is1,i2-2*is2,i3-2*is3,ex)=f(3 )*omega+(1.-omega)*q(3)
      wk1(i1-2*is1,i2-2*is2,i3-2*is3,ey)=f(4 )*omega+(1.-omega)*q(4)
      wk1(i1-2*is1,i2-2*is2,i3-2*is3,ez)=f(5 )*omega+(1.-omega)*q(5)
  
      wk2(j1-js1,j2-js2,j3-js3,ex)      =f(6 )*omega+(1.-omega)*q(6)
      wk2(j1-js1,j2-js2,j3-js3,ey)      =f(7 )*omega+(1.-omega)*q(7)
      wk2(j1-js1,j2-js2,j3-js3,ez)      =f(8 )*omega+(1.-omega)*q(8)
      wk2(j1-2*js1,j2-2*js2,j3-2*js3,ex)=f(9 )*omega+(1.-omega)*q(9)
      wk2(j1-2*js1,j2-2*js2,j3-2*js3,ey)=f(10)*omega+(1.-omega)*q(10)
      wk2(j1-2*js1,j2-2*js2,j3-2*js3,ez)=f(11)*omega+(1.-omega)*q(11)

    end if
    ! old : Nov 13, 2018
    ! u1(i1-is1,i2-is2,i3-is3,ex)      =f(0)
    ! u1(i1-is1,i2-is2,i3-is3,ey)      =f(1)
    ! u1(i1-is1,i2-is2,i3-is3,ez)      =f(2)
    ! u1(i1-2*is1,i2-2*is2,i3-2*is3,ex)=f(3)
    ! u1(i1-2*is1,i2-2*is2,i3-2*is3,ey)=f(4)
    ! u1(i1-2*is1,i2-2*is2,i3-2*is3,ez)=f(5)

    ! u2(j1-js1,j2-js2,j3-js3,ex)      =f(6)
    ! u2(j1-js1,j2-js2,j3-js3,ey)      =f(7)
    ! u2(j1-js1,j2-js2,j3-js3,ez)      =f(8)
    ! u2(j1-2*js1,j2-2*js2,j3-2*js3,ex)=f(9)
    ! u2(j1-2*js1,j2-2*js2,j3-2*js3,ey)=f(10)
    ! u2(j1-2*js1,j2-2*js2,j3-2*js3,ez)=f(11)

     ! compute the maximum change in the solution for this iteration
    do n=0,11
      err=max(err,abs(q(n)-f(n)))
      err2 = err2 + (q(n)-f(n))**2
      count = count + 1
    end do

    if( debug.gt.3 )then ! re-evaluate
     evalInterfaceDerivatives3d()
     eval3dJumpOrder4()
     write(debugFile,'(" --> 3d-order4-c: i1,i2,i3=",3i4," f(re-eval)=",12e10.2)') i1,i2,i3,f(0),f(1),f(2),f(3),f(4),f(5),f(6),f(7),f(8),f(9),f(10),f(11)
       ! '
    end if

  endLoopsMask3d()
  err2 = sqrt( err2/count )

  if( useJacobiUpdate.ne.0 )then
   ! Jacobi update -- copy work-space to solution arrays
    beginLoopsMask3d()
      u1(i1-is1,i2-is2,i3-is3,ex)      = wk1(i1-is1,i2-is2,i3-is3,ex)      
      u1(i1-is1,i2-is2,i3-is3,ey)      = wk1(i1-is1,i2-is2,i3-is3,ey)      
      u1(i1-is1,i2-is2,i3-is3,ez)      = wk1(i1-is1,i2-is2,i3-is3,ez)      
      u1(i1-2*is1,i2-2*is2,i3-2*is3,ex)= wk1(i1-2*is1,i2-2*is2,i3-2*is3,ex)
      u1(i1-2*is1,i2-2*is2,i3-2*is3,ey)= wk1(i1-2*is1,i2-2*is2,i3-2*is3,ey)
      u1(i1-2*is1,i2-2*is2,i3-2*is3,ez)= wk1(i1-2*is1,i2-2*is2,i3-2*is3,ez)
                                                                           
      u2(j1-js1,j2-js2,j3-js3,ex)      = wk2(j1-js1,j2-js2,j3-js3,ex)      
      u2(j1-js1,j2-js2,j3-js3,ey)      = wk2(j1-js1,j2-js2,j3-js3,ey)      
      u2(j1-js1,j2-js2,j3-js3,ez)      = wk2(j1-js1,j2-js2,j3-js3,ez)      
      u2(j1-2*js1,j2-2*js2,j3-2*js3,ex)= wk2(j1-2*js1,j2-2*js2,j3-2*js3,ex)
      u2(j1-2*js1,j2-2*js2,j3-2*js3,ey)= wk2(j1-2*js1,j2-2*js2,j3-2*js3,ey)
      u2(j1-2*js1,j2-2*js2,j3-2*js3,ez)= wk2(j1-2*js1,j2-2*js2,j3-2*js3,ez)
    endLoopsMask3d()
  end if

#endMacro


 ! here are the macros from deriv.maple (file=derivMacros.h)
 ! defineMacro lapCoeff4a(is,dr,ds) ( (-2/3.*rxx*is-2/3.*ryy*is)/dr+(4/3.*rx**2+4/3.*ry**2)/dr**2 )
 ! defineMacro lapCoeff4b(is,dr,ds) ( (1/12.*rxx*is+1/12.*ryy*is)/dr+(-1/12.*rx**2-1/12.*ry**2)/dr**2 )
 ! -- these are from op/src/derivCoeff.h (generated by op/src/derivCoeff.maple)
 #defineMacro lapSqCoeff3DOrder2a(is,dr,ds,dt,a0,a1,a2,RX,RXX,RXXX,RXXXX) ( -1/2.*(2*RXXXX(a0,0,0,1,1)+RXXXX(a0,0,0,0,0)+2*RXXXX(a0,0,0,2,2)+RXXXX(a0,1,1,1,1)+RXXXX(a0,2,2,2,2)+2*RXXXX(a0,1,1,2,2))*is/dr+(4*RX(a0,2)*RX(a0,0)*RXX(a0,0,2)+2*RX(a0,2)*(2*RX(a0,0)*RXX(a0,0,2)+RX(a0,2)*RXX(a0,0,0))+2*RXX(a0,2,2)*RX(a0,0)**2+4*RX(a0,1)*RXX(a0,0,1)*RX(a0,0)+2*RX(a0,1)*(2*RXX(a0,0,1)*RX(a0,0)+RX(a0,1)*RXX(a0,0,0))+2*RXX(a0,1,1)*RX(a0,0)**2+4*RX(a0,2)*RX(a0,1)*RXX(a0,1,2)+2*RX(a0,2)*(2*RX(a0,1)*RXX(a0,1,2)+RX(a0,2)*RXX(a0,1,1))+2*RXX(a0,2,2)*RX(a0,1)**2+6*RX(a0,1)**2*RXX(a0,1,1)+6*RX(a0,0)**2*RXX(a0,0,0)+6*RX(a0,2)**2*RXX(a0,2,2))*is/dr**3-2*(2*RX(a0,2)*(RX(a0,2)*RX(a1,1)**2+2*RX(a1,2)*RX(a1,1)*RX(a0,1))+2*RX(a1,2)*(RX(a1,2)*RX(a0,1)**2+2*RX(a0,2)*RX(a1,1)*RX(a0,1))+2*RX(a0,2)*(RX(a0,2)*RX(a1,0)**2+2*RX(a1,2)*RX(a1,0)*RX(a0,0))+2*RX(a1,2)*(RX(a1,2)*RX(a0,0)**2+2*RX(a0,2)*RX(a1,0)*RX(a0,0))+6*RX(a0,1)**2*RX(a1,1)**2+2*RX(a1,1)*(RX(a1,1)*RX(a0,0)**2+2*RX(a0,1)*RX(a1,0)*RX(a0,0))+2*RX(a0,1)*(RX(a0,1)*RX(a1,0)**2+2*RX(a1,1)*RX(a1,0)*RX(a0,0))+6*RX(a0,0)**2*RX(a1,0)**2+6*RX(a0,2)**2*RX(a1,2)**2)/dr**2/ds**2-2*(6*RX(a0,2)**2*RX(a2,2)**2+6*RX(a0,0)**2*RX(a2,0)**2+2*RX(a0,2)*(2*RX(a2,2)*RX(a2,0)*RX(a0,0)+RX(a0,2)*RX(a2,0)**2)+2*RX(a2,2)*(RX(a2,2)*RX(a0,0)**2+2*RX(a0,2)*RX(a2,0)*RX(a0,0))+2*RX(a0,2)*(2*RX(a2,2)*RX(a2,1)*RX(a0,1)+RX(a0,2)*RX(a2,1)**2)+2*RX(a2,2)*(RX(a2,2)*RX(a0,1)**2+2*RX(a0,2)*RX(a2,1)*RX(a0,1))+2*RX(a0,1)*(2*RX(a2,1)*RX(a2,0)*RX(a0,0)+RX(a0,1)*RX(a2,0)**2)+2*RX(a2,1)*(RX(a2,1)*RX(a0,0)**2+2*RX(a0,1)*RX(a2,0)*RX(a0,0))+6*RX(a0,1)**2*RX(a2,1)**2)/dr**2/dt**2+(4*RXX(a1,1,1)*RX(a1,0)*RX(a0,0)+2*RX(a1,2)*(2*RX(a0,0)*RXX(a1,0,2)+RX(a0,2)*RXX(a1,0,0)+2*RXX(a0,0,2)*RX(a1,0)+RX(a1,2)*RXX(a0,0,0))+2*RX(a1,2)*(2*RXX(a0,0,2)*RX(a1,0)+2*RX(a0,0)*RXX(a1,0,2))+7*RX(a0,0)*RXX(a1,0,0)*RX(a1,0)+7*RX(a0,1)*RXX(a1,1,1)*RX(a1,1)+4*RX(a0,1)*RXX(a1,0,1)*RX(a1,0)+4*RXX(a1,2,2)*RX(a1,0)*RX(a0,0)+RX(a1,2)*(3*RXX(a0,2,2)*RX(a1,2)+3*RX(a0,2)*RXX(a1,2,2))+RX(a1,2)*(2*RXX(a0,2,2)*RX(a1,2)+2*RX(a0,2)*RXX(a1,2,2))+4*RX(a0,2)*RX(a1,0)*RXX(a1,0,2)+2*RX(a1,2)*(2*RX(a0,1)*RXX(a1,1,2)+RX(a0,2)*RXX(a1,1,1)+2*RXX(a0,1,2)*RX(a1,1)+RX(a1,2)*RXX(a0,1,1))+RXX(a0,0,0)*RX(a1,0)**2+RX(a1,0)*(3*RXX(a1,0,0)*RX(a0,0)+3*RX(a1,0)*RXX(a0,0,0))+RX(a1,0)*(2*RX(a1,0)*RXX(a0,0,0)+2*RXX(a1,0,0)*RX(a0,0))+2*RX(a1,2)*(2*RXX(a0,1,2)*RX(a1,1)+2*RX(a0,1)*RXX(a1,1,2))+7*RX(a0,2)*RXX(a1,2,2)*RX(a1,2)+2*RX(a0,2)*(2*RX(a1,1)*RXX(a1,1,2)+RX(a1,2)*RXX(a1,1,1))+RXX(a0,2,2)*RX(a1,2)**2+2*RXX(a0,2,2)*RX(a1,1)**2+RXX(a0,1,1)*RX(a1,1)**2+RX(a1,1)*(3*RXX(a1,1,1)*RX(a0,1)+3*RX(a1,1)*RXX(a0,1,1))+RX(a1,1)*(2*RX(a1,1)*RXX(a0,1,1)+2*RXX(a1,1,1)*RX(a0,1))+4*RXX(a1,2,2)*RX(a1,1)*RX(a0,1)+4*RX(a0,2)*RX(a1,1)*RXX(a1,1,2)+2*RX(a0,1)*(RX(a1,1)*RXX(a1,0,0)+2*RXX(a1,0,1)*RX(a1,0))+2*RXX(a0,1,1)*RX(a1,0)**2+2*RX(a1,1)*(RX(a0,1)*RXX(a1,0,0)+2*RX(a1,0)*RXX(a0,0,1)+2*RXX(a1,0,1)*RX(a0,0)+RX(a1,1)*RXX(a0,0,0))+2*RX(a1,1)*(2*RX(a1,0)*RXX(a0,0,1)+2*RXX(a1,0,1)*RX(a0,0))+2*RX(a0,2)*(2*RX(a1,0)*RXX(a1,0,2)+RX(a1,2)*RXX(a1,0,0))+2*RXX(a0,2,2)*RX(a1,0)**2)*is/dr/ds**2+(2*RX(a2,2)*(2*RX(a2,0)*RXX(a0,0,2)+2*RXX(a2,0,2)*RX(a0,0))+2*RX(a2,2)*(2*RXX(a2,0,2)*RX(a0,0)+RX(a2,2)*RXX(a0,0,0)+2*RX(a2,0)*RXX(a0,0,2)+RX(a0,2)*RXX(a2,0,0))+4*RX(a0,2)*RXX(a2,0,2)*RX(a2,0)+2*RXX(a0,2,2)*RX(a2,1)**2+2*RX(a2,2)*(2*RX(a2,1)*RXX(a0,1,2)+2*RXX(a2,1,2)*RX(a0,1))+2*RX(a2,2)*(2*RXX(a2,1,2)*RX(a0,1)+RX(a2,2)*RXX(a0,1,1)+2*RX(a2,1)*RXX(a0,1,2)+RX(a0,2)*RXX(a2,1,1))+2*RX(a0,2)*(2*RXX(a2,1,2)*RX(a2,1)+RX(a2,2)*RXX(a2,1,1))+4*RX(a0,2)*RXX(a2,1,2)*RX(a2,1)+4*RXX(a2,2,2)*RX(a2,1)*RX(a0,1)+4*RXX(a2,1,1)*RX(a2,0)*RX(a0,0)+4*RX(a0,1)*RXX(a2,0,1)*RX(a2,0)+2*RX(a0,1)*(2*RXX(a2,0,1)*RX(a2,0)+RX(a2,1)*RXX(a2,0,0))+2*RXX(a0,1,1)*RX(a2,0)**2+2*RX(a2,1)*(2*RX(a2,0)*RXX(a0,0,1)+2*RXX(a2,0,1)*RX(a0,0))+2*RX(a2,1)*(2*RXX(a2,0,1)*RX(a0,0)+RX(a2,1)*RXX(a0,0,0)+2*RX(a2,0)*RXX(a0,0,1)+RX(a0,1)*RXX(a2,0,0))+RXX(a0,1,1)*RX(a2,1)**2+2*RX(a0,2)*(2*RXX(a2,0,2)*RX(a2,0)+RX(a2,2)*RXX(a2,0,0))+RX(a2,1)*(2*RX(a2,1)*RXX(a0,1,1)+2*RXX(a2,1,1)*RX(a0,1))+RX(a2,1)*(3*RXX(a2,1,1)*RX(a0,1)+3*RX(a2,1)*RXX(a0,1,1))+7*RX(a0,0)*RXX(a2,0,0)*RX(a2,0)+7*RX(a0,1)*RXX(a2,1,1)*RX(a2,1)+RXX(a0,2,2)*RX(a2,2)**2+RX(a2,2)*(2*RX(a2,2)*RXX(a0,2,2)+2*RXX(a2,2,2)*RX(a0,2))+RX(a2,2)*(3*RXX(a2,2,2)*RX(a0,2)+3*RX(a2,2)*RXX(a0,2,2))+RXX(a0,0,0)*RX(a2,0)**2+4*RXX(a2,2,2)*RX(a2,0)*RX(a0,0)+RX(a2,0)*(2*RX(a2,0)*RXX(a0,0,0)+2*RXX(a2,0,0)*RX(a0,0))+RX(a2,0)*(3*RXX(a2,0,0)*RX(a0,0)+3*RX(a2,0)*RXX(a0,0,0))+2*RXX(a0,2,2)*RX(a2,0)**2+7*RX(a0,2)*RXX(a2,2,2)*RX(a2,2))*is/dr/dt**2+(4*RX(a0,2)*RXXX(a0,1,1,2)+2*RXX(a0,2,2)*RXX(a0,1,1)+4*RXX(a0,1,2)**2+4*RX(a0,1)*RXXX(a0,1,2,2)+4*RX(a0,2)*RXXX(a0,0,0,2)+2*RXX(a0,2,2)*RXX(a0,0,0)+4*RXX(a0,0,2)**2+4*RX(a0,0)*RXXX(a0,0,2,2)+4*RX(a0,0)*RXXX(a0,0,0,0)+3*RXX(a0,0,0)**2+4*RX(a0,0)*RXXX(a0,0,1,1)+2*RXX(a0,1,1)*RXX(a0,0,0)+4*RX(a0,1)*RXXX(a0,0,0,1)+4*RXX(a0,0,1)**2+4*RX(a0,1)*RXXX(a0,1,1,1)+3*RXX(a0,1,1)**2+3*RXX(a0,2,2)**2+4*RX(a0,2)*RXXX(a0,2,2,2))/dr**2-4*(2*RX(a0,2)**2*RX(a0,1)**2+RX(a0,2)**4+RX(a0,1)**4+2*RX(a0,2)**2*RX(a0,0)**2+2*RX(a0,1)**2*RX(a0,0)**2+RX(a0,0)**4)/dr**4 )
 
 
 #defineMacro lapSqCoeff3DOrder2b(is,dr,ds,dt,a0,a1,a2,RX,RXX,RXXX,RXXXX) ( -1/2.*(4*RX(a0,1)*RXX(a0,0,1)*RX(a0,0)+2*RX(a0,1)*(2*RXX(a0,0,1)*RX(a0,0)+RX(a0,1)*RXX(a0,0,0))+2*RXX(a0,1,1)*RX(a0,0)**2+4*RX(a0,2)*RX(a0,0)*RXX(a0,0,2)+2*RX(a0,2)*(2*RX(a0,0)*RXX(a0,0,2)+RX(a0,2)*RXX(a0,0,0))+2*RXX(a0,2,2)*RX(a0,0)**2+6*RX(a0,1)**2*RXX(a0,1,1)+6*RX(a0,0)**2*RXX(a0,0,0)+4*RX(a0,2)*RX(a0,1)*RXX(a0,1,2)+2*RX(a0,2)*(2*RX(a0,1)*RXX(a0,1,2)+RX(a0,2)*RXX(a0,1,1))+2*RXX(a0,2,2)*RX(a0,1)**2+6*RX(a0,2)**2*RXX(a0,2,2))*is/dr**3+(2*RX(a0,2)**2*RX(a0,1)**2+RX(a0,2)**4+RX(a0,1)**4+2*RX(a0,2)**2*RX(a0,0)**2+2*RX(a0,1)**2*RX(a0,0)**2+RX(a0,0)**4)/dr**4 )
 
 #defineMacro xLapCoeff3DOrder2a(is,dr,ds,dt,a0,a1,a2,RX,RXX,RXXX) ( -1/2.*RXXX(a0,0,2,2)*is/dr-1/2.*RXXX(a0,0,0,0)*is/dr-1/2.*RXXX(a0,0,1,1)*is/dr+3*RX(a0,0)*RXX(a0,0,0)/dr**2+RX(a0,0)**3*is/dr**3+RX(a0,2)**2*RX(a0,0)*is/dr**3+(RX(a0,2)*RX(a1,2)*RX(a1,0)+RX(a1,2)*(RX(a1,2)*RX(a0,0)+RX(a0,2)*RX(a1,0)))*is/dr/ds**2+(RX(a2,2)*(RX(a0,2)*RX(a2,0)+RX(a2,2)*RX(a0,0))+RX(a0,2)*RX(a2,2)*RX(a2,0))*is/dr/dt**2+(RX(a2,1)*(RX(a0,1)*RX(a2,0)+RX(a2,1)*RX(a0,0))+RX(a0,1)*RX(a2,1)*RX(a2,0))*is/dr/dt**2+(RXX(a0,1,1)*RX(a0,0)+2*RX(a0,1)*RXX(a0,0,1))/dr**2+3*RX(a0,0)*RX(a2,0)**2*is/dr/dt**2+RX(a0,1)**2*RX(a0,0)*is/dr**3+3*RX(a0,0)*RX(a1,0)**2*is/dr/ds**2+(RX(a0,1)*RX(a1,1)*RX(a1,0)+RX(a1,1)*(RX(a1,1)*RX(a0,0)+RX(a0,1)*RX(a1,0)))*is/dr/ds**2+(2*RX(a0,2)*RXX(a0,0,2)+RXX(a0,2,2)*RX(a0,0))/dr**2 )
 
 ! #defineMacro xLapCoeff3DOrder2b(is,dr,ds) ( -1/2.*rz**2*rx*is/dr**3-1/2.*ry**2*rx*is/dr**3-1/2.*rx**3*is/dr**3 )
 
 #defineMacro xLapCoeff3DOrder2b(is,dr,ds,dt,a0,a1,a2,RX,RXX,RXXX) ( -1/2.*RX(a0,2)**2*RX(a0,0)*is/dr**3-1/2.*RX(a0,1)**2*RX(a0,0)*is/dr**3-1/2.*RX(a0,0)**3*is/dr**3 )
 
 ! #defineMacro yLapCoeff3DOrder2a(is,dr,ds) ( (2*rxy*rx+ry*rxx)/dr**2+rz**2*ry*is/dr**3+(rz*sz*sy+sz*(sz*ry+rz*sy))*is/dr/ds**2+(2*ty*tx*rx+ry*tx**2)*is/dr/dt**2+(tz*(rz*ty+tz*ry)+rz*tz*ty)*is/dr/dt**2+3*ry*ty**2*is/dr/dt**2+(2*rz*ryz+rzz*ry)/dr**2-1/2.*rxxy*is/dr-1/2.*ryyy*is/dr+3*ry*ryy/dr**2+ry**3*is/dr**3+ry*rx**2*is/dr**3+(ry*sx**2+2*sy*sx*rx)*is/dr/ds**2+3*ry*sy**2*is/dr/ds**2-1/2.*ryzz*is/dr )
 
 #defineMacro yLapCoeff3DOrder2a(is,dr,ds,dt,a0,a1,a2,RX,RXX,RXXX) ( (2*RXX(a0,0,1)*RX(a0,0)+RX(a0,1)*RXX(a0,0,0))/dr**2+RX(a0,2)**2*RX(a0,1)*is/dr**3+(RX(a0,2)*RX(a1,2)*RX(a1,1)+RX(a1,2)*(RX(a1,2)*RX(a0,1)+RX(a0,2)*RX(a1,1)))*is/dr/ds**2+(2*RX(a2,1)*RX(a2,0)*RX(a0,0)+RX(a0,1)*RX(a2,0)**2)*is/dr/dt**2+(RX(a2,2)*(RX(a0,2)*RX(a2,1)+RX(a2,2)*RX(a0,1))+RX(a0,2)*RX(a2,2)*RX(a2,1))*is/dr/dt**2+3*RX(a0,1)*RX(a2,1)**2*is/dr/dt**2+(2*RX(a0,2)*RXX(a0,1,2)+RXX(a0,2,2)*RX(a0,1))/dr**2-1/2.*RXXX(a0,0,0,1)*is/dr-1/2.*RXXX(a0,1,1,1)*is/dr+3*RX(a0,1)*RXX(a0,1,1)/dr**2+RX(a0,1)**3*is/dr**3+RX(a0,1)*RX(a0,0)**2*is/dr**3+(RX(a0,1)*RX(a1,0)**2+2*RX(a1,1)*RX(a1,0)*RX(a0,0))*is/dr/ds**2+3*RX(a0,1)*RX(a1,1)**2*is/dr/ds**2-1/2.*RXXX(a0,1,2,2)*is/dr )
 
 ! #defineMacro yLapCoeff3DOrder2b(is,dr,ds) ( -1/2.*rz**2*ry*is/dr**3-1/2.*ry**3*is/dr**3-1/2.*ry*rx**2*is/dr**3 )
 
 #defineMacro yLapCoeff3DOrder2b(is,dr,ds,dt,a0,a1,a2,RX,RXX,RXXX) ( -1/2.*RX(a0,2)**2*RX(a0,1)*is/dr**3-1/2.*RX(a0,1)**3*is/dr**3-1/2.*RX(a0,1)*RX(a0,0)**2*is/dr**3 )
 
 ! #defineMacro zLapCoeff3DOrder2a(is,dr,ds) ( -1/2.*rxxz*is/dr-1/2.*ryyz*is/dr+3*rz*rzz/dr**2+rz**3*is/dr**3+(2*tz*tx*rx+rz*tx**2)*is/dr/dt**2+rz*rx**2*is/dr**3+3*rz*sz**2*is/dr/ds**2-1/2.*rzzz*is/dr+(rz*sx**2+2*sz*sx*rx)*is/dr/ds**2+3*rz*tz**2*is/dr/dt**2+(rz*sy**2+2*sz*sy*ry)*is/dr/ds**2+(2*tz*ty*ry+rz*ty**2)*is/dr/dt**2+rz*ry**2*is/dr**3+(2*rx*rxz+rz*rxx)/dr**2+(2*ry*ryz+rz*ryy)/dr**2 )
 
 #defineMacro zLapCoeff3DOrder2a(is,dr,ds,dt,a0,a1,a2,RX,RXX,RXXX) ( -1/2.*RXXX(a0,0,0,2)*is/dr-1/2.*RXXX(a0,1,1,2)*is/dr+3*RX(a0,2)*RXX(a0,2,2)/dr**2+RX(a0,2)**3*is/dr**3+(2*RX(a2,2)*RX(a2,0)*RX(a0,0)+RX(a0,2)*RX(a2,0)**2)*is/dr/dt**2+RX(a0,2)*RX(a0,0)**2*is/dr**3+3*RX(a0,2)*RX(a1,2)**2*is/dr/ds**2-1/2.*RXXX(a0,2,2,2)*is/dr+(RX(a0,2)*RX(a1,0)**2+2*RX(a1,2)*RX(a1,0)*RX(a0,0))*is/dr/ds**2+3*RX(a0,2)*RX(a2,2)**2*is/dr/dt**2+(RX(a0,2)*RX(a1,1)**2+2*RX(a1,2)*RX(a1,1)*RX(a0,1))*is/dr/ds**2+(2*RX(a2,2)*RX(a2,1)*RX(a0,1)+RX(a0,2)*RX(a2,1)**2)*is/dr/dt**2+RX(a0,2)*RX(a0,1)**2*is/dr**3+(2*RX(a0,0)*RXX(a0,0,2)+RX(a0,2)*RXX(a0,0,0))/dr**2+(2*RX(a0,1)*RXX(a0,1,2)+RX(a0,2)*RXX(a0,1,1))/dr**2 )
 
 ! #defineMacro zLapCoeff3DOrder2b(is,dr,ds) ( -1/2.*rz**3*is/dr**3-1/2.*rz*ry**2*is/dr**3-1/2.*rz*rx**2*is/dr**3 )
 
 #defineMacro zLapCoeff3DOrder2b(is,dr,ds,dt,a0,a1,a2,RX,RXX,RXXX) ( -1/2.*RX(a0,2)**3*is/dr**3-1/2.*RX(a0,2)*RX(a0,1)**2*is/dr**3-1/2.*RX(a0,2)*RX(a0,0)**2*is/dr**3 )


! ----------------------------------------------------------------------------------
!  Macro:
!    Evaluate the interface equations for checking the coefficients
! ----------------------------------------------------------------------------------
#beginMacro evalInterfaceEquations34c()
  evalDerivs3dOrder4()
  eval3dJumpOrder4New() 
#endMacro


! ------------------------------------------------------------------------------------------------
! ------------------------------------------------------------------------------------------------
! Macro: Assign interface ghost values, DIM=3, ORDER=4, GRID=Curvilinear
! ------------------------------------------------------------------------------------------------
! ------------------------------------------------------------------------------------------------
#beginMacro assignInterfaceGhost3dOrder4()

  err=0.
  err2=0.
  count=0
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

    cem1=(1.-1./eps1)/mu1
    cem2=(1.-1./eps2)/mu2


    ! old evalInterfaceDerivatives3d()
    ! old eval3dJumpOrder4()

    ! *new*
    evalDerivs3dOrder4()
    eval3dJumpOrder4New()

    if( debug.gt.4 )then
     write(debugFile,'(/," --> 3d-order4-curv: i1,i2,i3=",3i4," an1,an2,an3=",3e11.3)') i1,i2,i3,an1,an2,an3
     write(debugFile,'(" --> 3d-order4-curv: i1,i2,i3=",3i4," f(start)=",12e10.2)') i1,i2,i3,(f(n),n=0,11)
     ! '
    end if
    if( debug.gt.8 )then
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


    c1m1x = -is*8.*rsxy1(i1,i2,i3,axis1,0)/(12.*dr1(axis1))    ! coeff of u1(-1) from D.x
    c1m1y = -is*8.*rsxy1(i1,i2,i3,axis1,1)/(12.*dr1(axis1))    ! coeff of u1(-1) from D.y 
    c1m1z = -is*8.*rsxy1(i1,i2,i3,axis1,2)/(12.*dr1(axis1))    ! coeff of u1(-1) from D.z

    c1m2x =  is   *rsxy1(i1,i2,i3,axis1,0)/(12.*dr1(axis1))    ! coeff of u1(-2) from D.x
    c1m2y =  is   *rsxy1(i1,i2,i3,axis1,1)/(12.*dr1(axis1))    ! coeff of u1(-2) from D.y 
    c1m2z =  is   *rsxy1(i1,i2,i3,axis1,2)/(12.*dr1(axis1))    ! coeff of u1(-2) from D.z

    c2m1x = -js*8.*rsxy2(j1,j2,j3,axis2,0)/(12.*dr2(axis2))    ! coeff of u2(-1) from D.x
    c2m1y = -js*8.*rsxy2(j1,j2,j3,axis2,1)/(12.*dr2(axis2))
    c2m1z = -js*8.*rsxy2(j1,j2,j3,axis2,2)/(12.*dr2(axis2))

    c2m2x =  js   *rsxy2(j1,j2,j3,axis2,0)/(12.*dr2(axis2))
    c2m2y =  js   *rsxy2(j1,j2,j3,axis2,1)/(12.*dr2(axis2))
    c2m2z =  js   *rsxy2(j1,j2,j3,axis2,2)/(12.*dr2(axis2))

    ! ----------- These should be the fourth-order accurate jacobian derivatives -------------
    rxx1(0,0,0)=aj1rxx
    rxx1(0,0,1)=aj1rxy
    rxx1(0,0,2)=aj1rxz
    rxx1(0,1,1)=aj1ryy
    rxx1(0,1,2)=aj1ryz
    rxx1(0,2,2)=aj1rzz

    rxx1(1,0,0)=aj1sxx
    rxx1(1,0,1)=aj1sxy
    rxx1(1,0,2)=aj1sxz
    rxx1(1,1,1)=aj1syy
    rxx1(1,1,2)=aj1syz
    rxx1(1,2,2)=aj1szz

    rxx1(2,0,0)=aj1txx
    rxx1(2,0,1)=aj1txy
    rxx1(2,0,2)=aj1txz
    rxx1(2,1,1)=aj1tyy
    rxx1(2,1,2)=aj1tyz
    rxx1(2,2,2)=aj1tzz


    rxx2(0,0,0)=aj2rxx
    rxx2(0,0,1)=aj2rxy
    rxx2(0,0,2)=aj2rxz
    rxx2(0,1,1)=aj2ryy
    rxx2(0,1,2)=aj2ryz
    rxx2(0,2,2)=aj2rzz

    rxx2(1,0,0)=aj2sxx
    rxx2(1,0,1)=aj2sxy
    rxx2(1,0,2)=aj2sxz
    rxx2(1,1,1)=aj2syy
    rxx2(1,1,2)=aj2syz
    rxx2(1,2,2)=aj2szz

    rxx2(2,0,0)=aj2txx
    rxx2(2,0,1)=aj2txy
    rxx2(2,0,2)=aj2txz
    rxx2(2,1,1)=aj2tyy
    rxx2(2,1,2)=aj2tyz
    rxx2(2,2,2)=aj2tzz


    ! clap1m1 : coeff of u(-1) from lap1 = u1.xx + u1.yy + u1.zz
    ! clap1m2 : coeff of u(-2) from lap1 = u1.xx + u1.yy + u1.zz

    
    clap1m1=4./3.*(rsxy1(i1,i2,i3,axis1,0)**2+rsxy1(i1,i2,i3,axis1,1)**2+rsxy1(i1,i2,i3,axis1,2)**2)/(dr1(axis1)**2) \
              -is*2./3.*(rxx1(axis1,0,0)+rxx1(axis1,1,1)+rxx1(axis1,2,2))/(2.*dr1(axis1))
    clap1m2=-1./12.*(rsxy1(i1,i2,i3,axis1,0)**2+rsxy1(i1,i2,i3,axis1,1)**2+rsxy1(i1,i2,i3,axis1,2)**2)/(dr1(axis1)**2) \
              +is*1./12.*(rxx1(axis1,0,0)+rxx1(axis1,1,1)+rxx1(axis1,2,2))/(2.*dr1(axis1)) 

    clap2m1=4/3.*(rsxy2(j1,j2,j3,axis2,0)**2+rsxy2(j1,j2,j3,axis2,1)**2+rsxy2(j1,j2,j3,axis2,2)**2)/(dr2(axis2)**2) \
              -js*2./3.*(rxx2(axis2,0,0)+rxx2(axis2,1,1)+rxx2(axis2,2,2))/(2.*dr2(axis2)) 
    clap2m2=-1./12.*(rsxy2(j1,j2,j3,axis2,0)**2+rsxy2(j1,j2,j3,axis2,1)**2+rsxy2(j1,j2,j3,axis2,2)**2)/(dr2(axis2)**2) \
              +js*1./12.*(rxx2(axis2,0,0)+rxx2(axis2,1,1)+rxx2(axis2,2,2))/(2.*dr2(axis2))

    ! cdivE1 =  u.c1x + v.c1y + w.c1z
    ! nDotCurlE1 = (w1y-v1z)*an1 + (u1z-w1x)*an2 + (v1x-u1y)*an3

    ! 12 Unknowns:
    !   u1(-1), v1(-1), w1(-1), u1(-2), v1(-2), w1(-2), 
    !   u2(-1), v2(-1), w2(-1), u2(-2), v2(-2), w2(-2)  
    ! 12 Equations: 
    !    a12(eqn,unknown) 


    ! fill equations 0,..,5
    fillEquations3dOrder4(0,1,2,3,4,5)


    ! rx1(m,n) = D r_m / D x_n
    ! rxx1(m,n1,n2) = D^2 r /( D x_n1 D X_n2 )
    ! rxxx1(m,n1,n2,n3) 

    ! ----------- These should be the SECOND-order accurate jacobian derivatives -------------

    ! first derivatives
    rx1(0,0)=rsxy1(i1,i2,i3,0,0)
    rx1(0,1)=rsxy1(i1,i2,i3,0,1)
    rx1(0,2)=rsxy1(i1,i2,i3,0,2)
    rx1(1,0)=rsxy1(i1,i2,i3,1,0)
    rx1(1,1)=rsxy1(i1,i2,i3,1,1)
    rx1(1,2)=rsxy1(i1,i2,i3,1,2)
    rx1(2,0)=rsxy1(i1,i2,i3,2,0)
    rx1(2,1)=rsxy1(i1,i2,i3,2,1)
    rx1(2,2)=rsxy1(i1,i2,i3,2,2)

    rx2(0,0)=rsxy2(j1,j2,j3,0,0)
    rx2(0,1)=rsxy2(j1,j2,j3,0,1)
    rx2(0,2)=rsxy2(j1,j2,j3,0,2)
    rx2(1,0)=rsxy2(j1,j2,j3,1,0)
    rx2(1,1)=rsxy2(j1,j2,j3,1,1)
    rx2(1,2)=rsxy2(j1,j2,j3,1,2)
    rx2(2,0)=rsxy2(j1,j2,j3,2,0)
    rx2(2,1)=rsxy2(j1,j2,j3,2,1)
    rx2(2,2)=rsxy2(j1,j2,j3,2,2)


    ! ------ FIX ME : We need SECOND-ORDER ACCURATE rxx1 etc. ****


    ! 3rd derivatives: (only some are needed)
    ! note for last 3 entries - entries must increase or stay the same
    rxxx1(axis1,0,0,0) = aj1rxxx
    rxxx1(axis1,1,1,1) = aj1ryyy
    rxxx1(axis1,2,2,2) = aj1rzzz
    rxxx1(axis1,0,0,1) = aj1rxxy
    rxxx1(axis1,0,0,2) = aj1rxxz
    rxxx1(axis1,0,1,1) = aj1rxyy
    rxxx1(axis1,0,2,2) = aj1rxzz
    rxxx1(axis1,1,1,2) = aj1ryyz
    rxxx1(axis1,1,2,2) = aj1ryzz

    rxxx2(axis2,0,0,0) = aj2rxxx
    rxxx2(axis2,1,1,1) = aj2ryyy
    rxxx2(axis2,2,2,2) = aj2rzzz
    rxxx2(axis2,0,0,1) = aj2rxxy
    rxxx2(axis2,0,0,2) = aj2rxxz
    rxxx2(axis2,0,1,1) = aj2rxyy
    rxxx2(axis2,0,2,2) = aj2rxzz
    rxxx2(axis2,1,1,2) = aj2ryyz
    rxxx2(axis2,1,2,2) = aj2ryzz



    ! Some 4th derivatives are needed by LapSq: 
    rxxxx1(axis1,0,0,0,0) = aj1rxxxx
    rxxxx1(axis1,1,1,1,1) = aj1ryyyy
    rxxxx1(axis1,2,2,2,2) = aj1rzzzz

    rxxxx1(axis1,0,0,1,1) = aj1rxxyy
    rxxxx1(axis1,0,0,2,2) = aj1rxxzz
    rxxxx1(axis1,1,1,2,2) = aj1ryyzz

    rxxxx2(axis2,0,0,0,0) = aj2rxxxx
    rxxxx2(axis2,1,1,1,1) = aj2ryyyy
    rxxxx2(axis2,2,2,2,2) = aj2rzzzz

    rxxxx2(axis2,0,0,1,1) = aj2rxxyy
    rxxxx2(axis2,0,0,2,2) = aj2rxxzz
    rxxxx2(axis2,1,1,2,2) = aj2ryyzz

    ! coeff of u1(-1) from D.x(Lap), D.y(Lap) and D.z(Lap): (divideb by eps*mu)
    ! dr1a, dr2a : used for avoidInterfaceIterations - tangential spacings are dsBig to eliminate mixed derivatives
    c1m1x = xLapCoeff3DOrder2a(is,dr1a(axis1),dr1a(axis1p1),dr1a(axis1p2),axis1,axis1p1,axis1p2,rx1,rxx1,rxxx1)/epsmu1 
    c1m1y = yLapCoeff3DOrder2a(is,dr1a(axis1),dr1a(axis1p1),dr1a(axis1p2),axis1,axis1p1,axis1p2,rx1,rxx1,rxxx1)/epsmu1  
    c1m1z = zLapCoeff3DOrder2a(is,dr1a(axis1),dr1a(axis1p1),dr1a(axis1p2),axis1,axis1p1,axis1p2,rx1,rxx1,rxxx1)/epsmu1  

    ! coeff of u1(-2) from D.x(Lap), D.y(Lap) and D.z(Lap): 
    c1m2x = xLapCoeff3DOrder2b(is,dr1a(axis1),dr1a(axis1p1),dr1a(axis1p2),axis1,axis1p1,axis1p2,rx1,rxx1,rxxx1)/epsmu1  
    c1m2y = yLapCoeff3DOrder2b(is,dr1a(axis1),dr1a(axis1p1),dr1a(axis1p2),axis1,axis1p1,axis1p2,rx1,rxx1,rxxx1)/epsmu1  
    c1m2z = zLapCoeff3DOrder2b(is,dr1a(axis1),dr1a(axis1p1),dr1a(axis1p2),axis1,axis1p1,axis1p2,rx1,rxx1,rxxx1)/epsmu1  


    ! coeff of u2(-1) from D.x(Lap), D.y(Lap) and D.z(Lap): 
    c2m1x = xLapCoeff3DOrder2a(js,dr2a(axis2),dr2a(axis2p1),dr2a(axis2p2),axis2,axis2p1,axis2p2,rx2,rxx2,rxxx2)/epsmu2  
    c2m1y = yLapCoeff3DOrder2a(js,dr2a(axis2),dr2a(axis2p1),dr2a(axis2p2),axis2,axis2p1,axis2p2,rx2,rxx2,rxxx2)/epsmu2 
    c2m1z = zLapCoeff3DOrder2a(js,dr2a(axis2),dr2a(axis2p1),dr2a(axis2p2),axis2,axis2p1,axis2p2,rx2,rxx2,rxxx2)/epsmu2  

    ! coeff of u2(-2) from D.x(Lap), D.y(Lap) and D.z(Lap): 
    c2m2x = xLapCoeff3DOrder2b(js,dr2a(axis2),dr2a(axis2p1),dr2a(axis2p2),axis2,axis2p1,axis2p2,rx2,rxx2,rxxx2)/epsmu2  
    c2m2y = yLapCoeff3DOrder2b(js,dr2a(axis2),dr2a(axis2p1),dr2a(axis2p2),axis2,axis2p1,axis2p2,rx2,rxx2,rxxx2)/epsmu2  
    c2m2z = zLapCoeff3DOrder2b(js,dr2a(axis2),dr2a(axis2p1),dr2a(axis2p2),axis2,axis2p1,axis2p2,rx2,rxx2,rxxx2)/epsmu2  


    ! coeff of u1(-1) and u1(-2) from Lap^2
    clap1m1=lapSqCoeff3DOrder2a(is,dr1a(axis1),dr1a(axis1p1),dr1a(axis1p2),axis1,axis1p1,axis1p2,rx1,rxx1,rxxx1,rxxxx1)/epsmu1
    clap1m2=lapSqCoeff3DOrder2b(is,dr1a(axis1),dr1a(axis1p1),dr1a(axis1p2),axis1,axis1p1,axis1p2,rx1,rxx1,rxxx1,rxxxx1)/epsmu1


    ! coeff of u2(-1) and u2(-2) from Lap^2
    clap2m1=lapSqCoeff3DOrder2a(js,dr2a(axis2),dr2a(axis2p1),dr2a(axis2p2),axis2,axis2p1,axis2p2,rx2,rxx2,rxxx2,rxxxx2)/epsmu2
    clap2m2=lapSqCoeff3DOrder2b(js,dr2a(axis2),dr2a(axis2p1),dr2a(axis2p2),axis2,axis2p1,axis2p2,rx2,rxx2,rxxx2,rxxxx2)/epsmu2

    ! write(debugFile,'(" --> 3d-order4-c: i1,i2,i3=",3i4," c1m1x,c1m1y,c1m1z,clap1m1,clap1m2=",12e10.2)') i1,i2,i3,c1m1x,c1m1y,c1m1z,clap1m1,clap1m2

    ! fill equations 6...11
    fillEquations3dOrder4(6,7,8,9,10,11)

    ! Equation 6..11 :  extrapolate 2nd ghost point 
    cex1=1.
    cex2=-5. ! ** fix me ** orderOfExtrapolation for 2nd ghost point 
    ! *e678*
    ! do ii=6,11
    !   do jj=0,11
    !     a12(ii,jj)=0.
    !   end do
    ! end do
    ! a12(6,0)  = cex2   ! u1(-1)
    ! a12(6,3)  = cex1   ! u1(-2)

    ! a12(7,1)  = cex2   ! v1(-1)
    ! a12(7,4)  = cex1   ! v1(-2)

    ! a12(8,2)  = cex2   ! w1(-1)
    ! a12(8,5)  = cex1   ! w1(-2)

    ! a12(9,6)  = cex2   ! u2(-1)
    ! a12(9,9)  = cex1   ! u2(-2)
    ! a12(10,7) = cex2   ! v2(-1)
    ! a12(10,10)= cex1   ! v2(-2)

    ! a12(11,8) = cex2   ! w2(-1)
    ! a12(11,11)= cex1   ! w2(-2)


    ! --- EVALUATE matrix coefficients by delta function approach ----
    if( .false. .and. checkCoeff.eq.1 )then
      numberOfEquations=12
      evalCoefficients(i1,i2,i3, j1,j2,j3,numberOfEquations,a12,evalInterfaceEquations34c )
    end if
    ! --- check matrix coefficients by delta function approach ----
    if( checkCoeff.eq.1 .and. it.le.1 )then
      numberOfEquations=12
      checkCoefficients(i1,i2,i3, j1,j2,j3,numberOfEquations,a12,evalInterfaceEquations34c )
    end if

    ! fill in the current values for the unknowns: 
    q(0) = u1(i1-is1,i2-is2,i3-is3,ex)
    q(1) = u1(i1-is1,i2-is2,i3-is3,ey)
    q(2) = u1(i1-is1,i2-is2,i3-is3,ez)
    q(3) = u1(i1-2*is1,i2-2*is2,i3-2*is3,ex)
    q(4) = u1(i1-2*is1,i2-2*is2,i3-2*is3,ey)
    q(5) = u1(i1-2*is1,i2-2*is2,i3-2*is3,ez)

    q(6) = u2(j1-js1,j2-js2,j3-js3,ex)
    q(7) = u2(j1-js1,j2-js2,j3-js3,ey)
    q(8) = u2(j1-js1,j2-js2,j3-js3,ez)
    q(9) = u2(j1-2*js1,j2-2*js2,j3-2*js3,ex)
    q(10)= u2(j1-2*js1,j2-2*js2,j3-2*js3,ey)
    q(11)= u2(j1-2*js1,j2-2*js2,j3-2*js3,ez)

    ! subtract off the contributions from the wrong values at the ghost points:
    numberOfEquations=12
    do n=0,numberOfEquations-1
      f(n) = (a12(n,0)*q(0)+a12(n,1)*q(1)+a12(n,2)*q(2)+a12(n,3)*q(3)+a12(n,4)*q(4)+a12(n,5)*q(5)+\
              a12(n,6)*q(6)+a12(n,7)*q(7)+a12(n,8)*q(8)+a12(n,9)*q(9)+a12(n,10)*q(10)+a12(n,11)*q(11) ) - f(n)
    end do
    if( debug.gt.3 )then
      write(debugFile,'(" --> 3d:order4-c: f(subtract)=",12e10.2)') (f(n),n=0,11)
    end if
    ! solve A Q = F
    ! factor the matrix
    call dgeco( a12(0,0), numberOfEquations, numberOfEquations, ipvt(0),rcond,work(0))
    ! solve
     if( debug.gt.3 )then
       write(debugFile,'(" --> 3d:order4-c: rcond=",e10.2)') rcond
     end if
    job=0
    call dgesl( a12(0,0), numberOfEquations, numberOfEquations, ipvt(0), f(0), job)
    if( debug.gt.3 )then
      write(debugFile,'(" --> 3d:order4-c: f(solve)=",12e10.2)') (f(n),n=0,11)
      write(debugFile,'(" --> 3d:order4-c:        q=",12e10.2)') (q(n),n=0,11)
      write(debugFile,'(" --> 3d:order4-c:      f-q=",12e10.2)') (f(n)-q(n),n=0,11)
    end if
    ! fill in the answer:
    if( useJacobiUpdate.eq.0 )then
      u1(i1-is1,i2-is2,i3-is3,ex)      =f(0 )*omega+(1.-omega)*q(0)
      u1(i1-is1,i2-is2,i3-is3,ey)      =f(1 )*omega+(1.-omega)*q(1)
      u1(i1-is1,i2-is2,i3-is3,ez)      =f(2 )*omega+(1.-omega)*q(2)
      u1(i1-2*is1,i2-2*is2,i3-2*is3,ex)=f(3 )*omega+(1.-omega)*q(3)
      u1(i1-2*is1,i2-2*is2,i3-2*is3,ey)=f(4 )*omega+(1.-omega)*q(4)
      u1(i1-2*is1,i2-2*is2,i3-2*is3,ez)=f(5 )*omega+(1.-omega)*q(5)
  
      u2(j1-js1,j2-js2,j3-js3,ex)      =f(6 )*omega+(1.-omega)*q(6)
      u2(j1-js1,j2-js2,j3-js3,ey)      =f(7 )*omega+(1.-omega)*q(7)
      u2(j1-js1,j2-js2,j3-js3,ez)      =f(8 )*omega+(1.-omega)*q(8)
      u2(j1-2*js1,j2-2*js2,j3-2*js3,ex)=f(9 )*omega+(1.-omega)*q(9)
      u2(j1-2*js1,j2-2*js2,j3-2*js3,ey)=f(10)*omega+(1.-omega)*q(10)
      u2(j1-2*js1,j2-2*js2,j3-2*js3,ez)=f(11)*omega+(1.-omega)*q(11)
    else
      ! Jacobi-update -- save results in the work space

      wk1(i1-is1,i2-is2,i3-is3,ex)      =f(0 )*omega+(1.-omega)*q(0)
      wk1(i1-is1,i2-is2,i3-is3,ey)      =f(1 )*omega+(1.-omega)*q(1)
      wk1(i1-is1,i2-is2,i3-is3,ez)      =f(2 )*omega+(1.-omega)*q(2)
      wk1(i1-2*is1,i2-2*is2,i3-2*is3,ex)=f(3 )*omega+(1.-omega)*q(3)
      wk1(i1-2*is1,i2-2*is2,i3-2*is3,ey)=f(4 )*omega+(1.-omega)*q(4)
      wk1(i1-2*is1,i2-2*is2,i3-2*is3,ez)=f(5 )*omega+(1.-omega)*q(5)
  
      wk2(j1-js1,j2-js2,j3-js3,ex)      =f(6 )*omega+(1.-omega)*q(6)
      wk2(j1-js1,j2-js2,j3-js3,ey)      =f(7 )*omega+(1.-omega)*q(7)
      wk2(j1-js1,j2-js2,j3-js3,ez)      =f(8 )*omega+(1.-omega)*q(8)
      wk2(j1-2*js1,j2-2*js2,j3-2*js3,ex)=f(9 )*omega+(1.-omega)*q(9)
      wk2(j1-2*js1,j2-2*js2,j3-2*js3,ey)=f(10)*omega+(1.-omega)*q(10)
      wk2(j1-2*js1,j2-2*js2,j3-2*js3,ez)=f(11)*omega+(1.-omega)*q(11)

    end if


    ! compute the maximum change in the solution for this iteration
    do n=0,11
      err=max(err,abs(q(n)-f(n)))
      err2 = err2 + (q(n)-f(n))**2 
      count=count+1
    end do

    if( debug.gt.3 )then ! re-evaluate

      ! To check residual we need to set answer
     u1(i1-is1,i2-is2,i3-is3,ex)      =f(0 )
     u1(i1-is1,i2-is2,i3-is3,ey)      =f(1 )
     u1(i1-is1,i2-is2,i3-is3,ez)      =f(2 )
     u1(i1-2*is1,i2-2*is2,i3-2*is3,ex)=f(3 )
     u1(i1-2*is1,i2-2*is2,i3-2*is3,ey)=f(4 )
     u1(i1-2*is1,i2-2*is2,i3-2*is3,ez)=f(5 )

     u2(j1-js1,j2-js2,j3-js3,ex)      =f(6 )
     u2(j1-js1,j2-js2,j3-js3,ey)      =f(7 )
     u2(j1-js1,j2-js2,j3-js3,ez)      =f(8 )
     u2(j1-2*js1,j2-2*js2,j3-2*js3,ex)=f(9 )
     u2(j1-2*js1,j2-2*js2,j3-2*js3,ey)=f(10)
     u2(j1-2*js1,j2-2*js2,j3-2*js3,ez)=f(11)

     evalDerivs3dOrder4()
     eval3dJumpOrder4New()

     res=0.
     do n=0,11
       res=max(res,abs(f(n)))
     end do
     if( .false. .and. res.gt.1.e-9 )then
       write(debugFile,'(" --> ERR: 3d-order4-c: it=",i3," i1,i2,i3=",3i4," f(re-eval)=",12e10.2," max-res=",1pe9.2)') it,i1,i2,i3,f(0),f(1),f(2),f(3),f(4),f(5),f(6),f(7),f(8),f(9),f(10),f(11),res
     end if

     ! display differences in iterates
     res=0.
     do n=0,11
       res=max(res,abs(q(n)-f(n)))
     end do
     if( .false. .and. res.gt.1.e-5 )then
       write(debugFile,'(" --> DIFF: 3d-order4-c: it=",i3," i1,i2,i3=",3i4," q-f=",12e10.2," max(q-f)=",1pe9.2)') it,i1,i2,i3,(q(n)-f(n),n=0,11),res
     end if

     if( useJacobiUpdate.ne.0 )then
      u1(i1-is1,i2-is2,i3-is3,ex)      =q(0 )
      u1(i1-is1,i2-is2,i3-is3,ey)      =q(1 )
      u1(i1-is1,i2-is2,i3-is3,ez)      =q(2 )
      u1(i1-2*is1,i2-2*is2,i3-2*is3,ex)=q(3 )
      u1(i1-2*is1,i2-2*is2,i3-2*is3,ey)=q(4 )
      u1(i1-2*is1,i2-2*is2,i3-2*is3,ez)=q(5 )

      u2(j1-js1,j2-js2,j3-js3,ex)      =q(6 )
      u2(j1-js1,j2-js2,j3-js3,ey)      =q(7 )
      u2(j1-js1,j2-js2,j3-js3,ez)      =q(8 )
      u2(j1-2*js1,j2-2*js2,j3-2*js3,ex)=q(9 )
      u2(j1-2*js1,j2-2*js2,j3-2*js3,ey)=q(10)
      u2(j1-2*js1,j2-2*js2,j3-2*js3,ez)=q(11)
     end if

    end if

  endLoopsMask3d()

  err2 = sqrt(err2/count) 

  if( checkCoeff.eq.1 .and. it.le.1 )then
    write(*,'("+++++ I34c: check coeff in interface: max(diff) = ",1pe8.2)') coeffDiff
  end if


  if( useJacobiUpdate.ne.0 )then
   ! Jacobi update -- copy work-space to solution arrays
    beginLoopsMask3d()
      u1(i1-is1,i2-is2,i3-is3,ex)      = wk1(i1-is1,i2-is2,i3-is3,ex)      
      u1(i1-is1,i2-is2,i3-is3,ey)      = wk1(i1-is1,i2-is2,i3-is3,ey)      
      u1(i1-is1,i2-is2,i3-is3,ez)      = wk1(i1-is1,i2-is2,i3-is3,ez)      
      u1(i1-2*is1,i2-2*is2,i3-2*is3,ex)= wk1(i1-2*is1,i2-2*is2,i3-2*is3,ex)
      u1(i1-2*is1,i2-2*is2,i3-2*is3,ey)= wk1(i1-2*is1,i2-2*is2,i3-2*is3,ey)
      u1(i1-2*is1,i2-2*is2,i3-2*is3,ez)= wk1(i1-2*is1,i2-2*is2,i3-2*is3,ez)
                                                                           
      u2(j1-js1,j2-js2,j3-js3,ex)      = wk2(j1-js1,j2-js2,j3-js3,ex)      
      u2(j1-js1,j2-js2,j3-js3,ey)      = wk2(j1-js1,j2-js2,j3-js3,ey)      
      u2(j1-js1,j2-js2,j3-js3,ez)      = wk2(j1-js1,j2-js2,j3-js3,ez)      
      u2(j1-2*js1,j2-2*js2,j3-2*js3,ex)= wk2(j1-2*js1,j2-2*js2,j3-2*js3,ex)
      u2(j1-2*js1,j2-2*js2,j3-2*js3,ey)= wk2(j1-2*js1,j2-2*js2,j3-2*js3,ey)
      u2(j1-2*js1,j2-2*js2,j3-2*js3,ez)= wk2(j1-2*js1,j2-2*js2,j3-2*js3,ez)
    endLoopsMask3d()
  end if

#endMacro

! *********************************************************************
! ********** MACROS FOR DISPERSIVE INTERFACE CONDITIONS ***************
! *********************************************************************
#Include "dispersiveInterfaceMacros.h"

! -------------------------------------------------------------------------------
! Macro: Evaluate the TZ forcings GDM FOURTH-ORDER 3D
! -------------------------------------------------------------------------------
#beginMacro evalTZForcingGDM3dOrder4(xy,i1,i2,i3,dispersionModel,numberOfPolarizationVectors,c,alphaP,a0v,a1v,b0v,b1v,fpv,fpSum,fev,\
                                   LfE,fEt,fEtt,LfP,fPt,fPtt,pevtt,pevttx,pevtty,pevttz,pevtttt,fevx,fevy,fevz,fpvx,fpvy,fpvz,pevttSum,pevttxSum,pevttySum,pevttzSum,pevttLSum,pevttttSum)

if( dispersionModel.ne.noDispersion )then
  do n=0,nd-1

    fpSum(n)=0.
    pevttSum(n)=0.

    pevttxSum(n)=0.
    pevttySum(n)=0.
    pevttzSum(n)=0.

    pevttLSum(n)=0.
    pevttttSum(n)=0.

    petttSum=0.

    call ogderiv(ep, 0,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, es(n)   ) 
    call ogderiv(ep, 1,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, est(n)  )
    call ogderiv(ep, 2,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, estt(n) )

    call ogderiv(ep, 0,1,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esx(n) )
    call ogderiv(ep, 0,0,1,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esy(n) )
    call ogderiv(ep, 0,0,0,1, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esz(n) )

    call ogderiv(ep, 0,2,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esxx(n) )
    call ogderiv(ep, 0,0,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esyy(n) )
    call ogderiv(ep, 0,0,0,2, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, eszz(n) )

    call ogderiv(ep, 1,1,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, estx(n) )
    call ogderiv(ep, 1,0,1,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esty(n) )
    call ogderiv(ep, 1,0,0,1, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, estz(n) )

    call ogderiv(ep, 2,1,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esttx(n) )
    call ogderiv(ep, 2,0,1,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, estty(n) )
    call ogderiv(ep, 2,0,0,1, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esttz(n) )

    call ogderiv(ep, 0,3,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esxxx(n) )
    call ogderiv(ep, 0,2,1,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esxxy(n) )
    call ogderiv(ep, 0,1,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esxyy(n) )
    call ogderiv(ep, 0,0,3,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esyyy(n) )

    call ogderiv(ep, 0,2,0,1, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esxxz(n) )
    call ogderiv(ep, 0,1,0,2, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esxzz(n) )
    call ogderiv(ep, 0,0,2,1, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esyyz(n) )
    call ogderiv(ep, 0,0,1,2, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esyzz(n) )
    call ogderiv(ep, 0,0,0,3, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, eszzz(n) )

    call ogderiv(ep, 3,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esttt(n) )
    call ogderiv(ep, 4,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, estttt(n) )

    call ogderiv(ep, 1,2,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, estxx(n) )
    call ogderiv(ep, 1,0,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, estyy(n) )
    call ogderiv(ep, 1,0,0,2, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, estzz(n) )


    call ogderiv(ep, 2,2,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esttxx(n) )
    call ogderiv(ep, 2,0,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esttyy(n) )
    call ogderiv(ep, 2,0,0,2, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esttzz(n) )

    call ogderiv(ep, 0,4,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esxxxx(n) )
    call ogderiv(ep, 0,2,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esxxyy(n) )
    call ogderiv(ep, 0,0,4,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esyyyy(n) )

    call ogderiv(ep, 0,0,0,4, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, eszzzz(n) )
    call ogderiv(ep, 0,2,0,2, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esxxzz(n) )
    call ogderiv(ep, 0,0,2,2, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esyyzz(n) )

    ! L = c^2*Delta
    esL  = c**2*( esxx(n)   + esyy(n)   + eszz(n)  )
    estL = c**2*( estxx(n)  + estyy(n)  + estzz(n) )
    esttL= c**2*( esttxx(n) + esttyy(n) + esttzz(n) )

    esLx  = c**2*( esxxx(n) + esxyy(n) + esxzz(n))
    esLy  = c**2*( esxxy(n) + esyyy(n) + esyzz(n) )
    esLz  = c**2*( esxxz(n) + esyyz(n) + eszzz(n) )

    ! L^2 : (xx + yy + zz)*( xx + yy + zz )
    esLL = c**4*( esxxxx(n) + esyyyy(n) + eszzzz(n) + 2.*( esxxyy(n) + esxxzz(n) + esyyzz(n)) )

    do jv=0,numberOfPolarizationVectors-1
      ! The TZ component is offset by pxc
      pc = pxc + jv*nd
      call ogderiv(ep, 0,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pe(n)   )
      call ogderiv(ep, 1,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pet(n)  )
      call ogderiv(ep, 2,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pett(n) )
      call ogderiv(ep, 3,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pettt(n) )
      call ogderiv(ep, 4,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, petttt(n) )

      call ogderiv(ep, 0,1,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pex(n) )
      call ogderiv(ep, 0,0,1,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pey(n) )
      call ogderiv(ep, 0,0,0,1, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pez(n) )

      call ogderiv(ep, 0,2,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pexx(n) )
      call ogderiv(ep, 0,0,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, peyy(n) )
      call ogderiv(ep, 0,0,0,2, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pezz(n) )

      call ogderiv(ep, 1,1,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, petx(n) )
      call ogderiv(ep, 1,0,1,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pety(n) )
      call ogderiv(ep, 1,0,0,1, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, petz(n) )

      call ogderiv(ep, 1,2,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, petxx(n) )
      call ogderiv(ep, 1,0,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, petyy(n) )
      call ogderiv(ep, 1,0,0,2, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, petzz(n) )

      call ogderiv(ep, 2,1,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pettx(n) )
      call ogderiv(ep, 2,0,1,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, petty(n) )
      call ogderiv(ep, 2,0,0,1, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pettz(n) )

      call ogderiv(ep, 2,2,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pettxx(n) )
      call ogderiv(ep, 2,0,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pettyy(n) )
      call ogderiv(ep, 2,0,0,2, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pettzz(n) )

      peL  = c**2*( pexx(n)   + peyy(n)   + pezz(n)  )
      petL = c**2*( petxx(n)  + petyy(n)  + petzz(n) )
      pettL= c**2*( pettxx(n) + pettyy(n) + pettzz(n) )

      ! Normal TZ forcing for P_{n,jv} equation: 
      fpv(n,jv) = pett(n)   + b1v(jv)*pet(n)   + b0v(jv)*pe(n)   - a0v(jv)*es(n)   - a1v(jv)*est(n)
      fPt(n,jv) = pettt(n)  + b1v(jv)*pett(n)  + b0v(jv)*pet(n)  - a0v(jv)*est(n)  - a1v(jv)*estt(n)
      fPtt(n,jv)= petttt(n) + b1v(jv)*pettt(n) + b0v(jv)*pett(n) - a0v(jv)*estt(n) - a1v(jv)*esttt(n)
      LfP(n,jv) = pettL     + b1v(jv)*petL     + b0v(jv)*peL     - a0v(jv)*esL     - a1v(jv)*estL

      fpvx(n,jv)= pettx(n)  + b1v(jv)*petx(n)  + b0v(jv)*pex(n)  - a0v(jv)*esx(n)  - a1v(jv)*estx(n)
      fpvy(n,jv)= petty(n)  + b1v(jv)*pety(n)  + b0v(jv)*pey(n)  - a0v(jv)*esy(n)  - a1v(jv)*esty(n)
      fpvz(n,jv)= pettz(n)  + b1v(jv)*petz(n)  + b0v(jv)*pez(n)  - a0v(jv)*esz(n)  - a1v(jv)*estz(n)

      ! write(*,'(" n=",i4," LfP=",e10.4," pettL,petL,peL,esL,estL=",5e12.4)') n,LfP(n,jv),pettL,petL,peL,esL,estL
      ! write(*,'(" pe,pet,pett,pettt,petttt=",5e12.4)') pe(n),pet(n),pett(n),pettt(n),petttt(n)

      ! write(*,'("TZ: n,jv=",2i4," pex,pey,pexx,peyy=",4(1pe12.4))') n,jv,pex(n),pey(n),pexx(n),peyy(n)

      ! Save ptt for checking later
      pevtt(n,jv)=pett(n)
      pevttx(n,jv)=pettx(n)
      pevtty(n,jv)=petty(n)
      pevttz(n,jv)=pettz(n)

      pevtttt(n,jv)=petttt(n)
      pevttLSum(n) = pevttLSum(n)  + pettL
      pevttttSum(n)= pevttttSum(n) + petttt(n) 

      ! Keep some sums: 
      fpSum(n)   = fpSum(n)  + fpv(n,jv)
      pevttSum(n)  = pevttSum(n)  + pett(n) 
      pevttxSum(n) = pevttxSum(n) + pettx(n)
      pevttySum(n) = pevttySum(n) + petty(n)
      pevttzSum(n) = pevttzSum(n) + pettz(n)

      petttSum  = petttSum  + pettt(n) 


    end do 

    ! TZ forcing for E_{n} equation:
    ! E_tt - c^2 Delta E + alphaP*Ptt  = 
    fev(n) = estt(n)   - esL   + alphaP*pevttSum(n)
    fEt(n) = esttt(n)  - estL  + alphaP*petttSum
    fEtt(n)= estttt(n) - esttL + alphaP*pevttttSum(n)

    fevx(n) = esttx(n) - esLx   + alphaP*pevttxSum(n)
    fevy(n) = estty(n) - esLy   + alphaP*pevttySum(n)
    fevz(n) = esttz(n) - esLz   + alphaP*pevttzSum(n)
    

    ! write(*,'("--> fEtt=",e10.2," estttt,esttL,pettttSum=",3e10.2)')  fEtt(n),estttt(n),esttL,pettttSum
    LfE(n) = esttL     - esLL  + alphaP*pevttLSum(n)
    
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
#beginMacro getDispersiveTZForcing3dOrder4(fpv1,fpv2,fev1,fev2)

  if( twilightZone.eq.1 )then
    evalTZForcingGDM3dOrder4(xy1,i1,i2,i3,dispersionModel1,numberOfPolarizationVectors1,c1,alphaP1,a0v1,a1v1,b0v1,b1v1,fpv1,fpSum1,fev1,\
                          LfE1,fEt1,fEtt1,LfP1,fPt1,fPtt1,pevtt1,pevttx1,pevtty1,pevttz1,pevtttt1,fevx1,fevy1,fevz1,fpvx1,fpvy1,fpvz1,pevttSum1,pevttxSum1,pevttySum1,pevttzSum1,pevttLSum1,pevttttSum1)

    evalTZForcingGDM3dOrder4(xy2,j1,j2,j3,dispersionModel2,numberOfPolarizationVectors2,c2,alphaP2,a0v2,a1v2,b0v2,b1v2,fpv2,fpSum2,fev2,\
                          LfE2,fEt2,fEtt2,LfP2,fPt2,fPtt2,pevtt2,pevttx2,pevtty2,pevttz2,pevtttt2,fevx2,fevy2,fevz2,fpvx2,fpvy2,fpvz2,pevttSum2,pevttxSum2,pevttySum2,pevttzSum2,pevttLSum2,pevttttSum2)
  end if

#endMacro 


! --------------------------------------------------------------------------
! Macro: Evaluate the GDM jump conditions in 3D, order=4
! --------------------------------------------------------------------------
#beginMacro eval3dJumpDispersiveOrder4()

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


 ! [ div(E) n + (I-nn^T)( curl(E)/mu ] = 0 
 f(0)= ( divE1*an1 + (curlE1x- nDotCurlE1*an1)/mu1 ) - ( divE2*an1 + (curlE2x- nDotCurlE2*an1)/mu2 )
 f(1)= ( divE1*an2 + (curlE1y- nDotCurlE1*an2)/mu1 ) - ( divE2*an2 + (curlE2y- nDotCurlE2*an2)/mu2 )
 f(2)= ( divE1*an3 + (curlE1z- nDotCurlE1*an3)/mu1 ) - ( divE2*an3 + (curlE2z- nDotCurlE2*an3)/mu2 )

 ! [  n n^T( Delta(E) /mu + (I-n n^T)( c^2 Delta(E) - alphaP*Ptt ] = 0 
 ! [  c^2 Delta(E) + (1/mu -c^2) n n^T Delta(E)+ (I-n n^T)( - alphaP*Ptt ] = 0 
 cem1=(1.-1./eps1)/mu1
 cem2=(1.-1./eps2)/mu2

 nDotFp1 = an1*fp1(0) + an2*fp1(1) + an3*fp1(2)
 nDotFp2 = an1*fp2(0) + an2*fp2(1) + an3*fp2(2)

 f(3)= ( u1Lap/(epsmu1) + cem1*nDotLapE1*an1 - alphaP1*( fp1(0)-an1*nDotFp1) ) - \
       ( u2Lap/(epsmu2) + cem2*nDotLapE2*an1 - alphaP2*( fp2(0)-an1*nDotFp2) )
 f(4)= ( v1Lap/(epsmu1) + cem1*nDotLapE1*an2 - alphaP1*( fp1(1)-an2*nDotFp1) ) - \
       ( v2Lap/(epsmu2) + cem2*nDotLapE2*an2 - alphaP2*( fp2(1)-an2*nDotFp2) )
 f(5)= ( w1Lap/(epsmu1) + cem1*nDotLapE1*an3 - alphaP1*( fp1(2)-an3*nDotFp1) ) - \
       ( w2Lap/(epsmu2) + cem2*nDotLapE2*an3 - alphaP2*( fp2(2)-an3*nDotFp2) )

 divLapE1 = u1xxx+v1xxy+w1xxz + u1xyy+v1yyy+w1yyz + u1xzz+v1yzz+w1zzz
 curlLapE1x = w1xxy-v1xxz + w1yyy-v1yyz + w1yzz-v1zzz
 curlLapE1y = u1xxz-w1xxx + u1yyz-w1xyy + u1zzz-w1xzz
 curlLapE1z = v1xxx-u1xxy + v1xyy-u1yyy + v1xzz-u1yzz
 nDotCurlLapE1=an1*curlLapE1x+an2*curlLapE1y+an3*curlLapE1z
 nDotLapSqE1 = an1*u1LapSq + an2*v1LapSq + an3*w1LapSq

 divLapE2 = u2xxx+v2xxy+w2xxz + u2xyy+v2yyy+w2yyz + u2xzz+v2yzz+w2zzz
 curlLapE2x = w2xxy-v2xxz + w2yyy-v2yyz + w2yzz-v2zzz
 curlLapE2y = u2xxz-w2xxx + u2yyz-w2xyy + u2zzz-w2xzz
 curlLapE2z = v2xxx-u2xxy + v2xyy-u2yyy + v2xzz-u2yzz
 nDotCurlLapE2=an1*curlLapE2x+an2*curlLapE2y+an3*curlLapE2z
 nDotLapSqE2 = an1*u2LapSq + an2*v2LapSq + an3*w2LapSq

 ! Compute curl(fPtt) and n^T curl(fPtt)
 curlfPttx1=fPtty1(2)-fPttz1(1) 
 curlfPtty1=fPttz1(0)-fPttx1(2)
 curlfPttz1=fPttx1(1)-fPtty1(0)
 nDotCurlfPtt1 = an1*curlfPttx1 + an2*curlfPtty1 + an3*curlfPttz1

 curlfPttx2=fPtty2(2)-fPttz2(1) 
 curlfPtty2=fPttz2(0)-fPttx2(2)
 curlfPttz2=fPttx2(1)-fPtty2(0)
 nDotCurlfPtt2 = an1*curlfPttx2 + an2*curlfPtty2 + an3*curlfPttz2

 ! Compute n^T( fLPtt) and n^T( fPtttt )
 nDotfLPtt1 = an1*fLPtt1(0)  + an2*fLPtt1(1)  + an3*fLPtt1(2)
 nDotfPtttt1= an1*fPtttt1(0) + an2*fPtttt1(1) + an3*fPtttt1(2) 

 nDotfLPtt2 = an1*fLPtt2(0)  + an2*fLPtt2(1)  + an3*fLPtt2(2)
 nDotfPtttt2= an1*fPtttt2(0) + an2*fPtttt2(1) + an3*fPtttt2(2) 

 ! [ c^2 div(Delta(E)) n + (1/mu)*( I - n n^T )( curl( c^2 Delta^2(E) -alphaP*P_tt) ] = 0
 !
 f(6)= ( ( divLapE1*an1 + (curlLapE1x- nDotCurlLapE1*an1)/mu1 )/(epsmu1) - (alphaP1/mu1)*( curlfPttx1 - an1*nDotCurlfPtt1 ) ) - \
       ( ( divLapE2*an1 + (curlLapE2x- nDotCurlLapE2*an1)/mu2 )/(epsmu2) - (alphaP2/mu2)*( curlfPttx2 - an1*nDotCurlfPtt2 ) )

 f(7)= ( ( divLapE1*an2 + (curlLapE1y- nDotCurlLapE1*an2)/mu1 )/(epsmu1) - (alphaP1/mu1)*( curlfPtty1 - an2*nDotCurlfPtt1 ) ) - \
       ( ( divLapE2*an2 + (curlLapE2y- nDotCurlLapE2*an2)/mu2 )/(epsmu2) - (alphaP2/mu2)*( curlfPtty2 - an2*nDotCurlfPtt2 ) )

 f(8)= ( ( divLapE1*an3 + (curlLapE1z- nDotCurlLapE1*an3)/mu1 )/(epsmu1) - (alphaP1/mu1)*( curlfPttz1 - an3*nDotCurlfPtt1 ) ) - \
       ( ( divLapE2*an3 + (curlLapE2z- nDotCurlLapE2*an3)/mu2 )/(epsmu2) - (alphaP2/mu2)*( curlfPttz2 - an3*nDotCurlfPtt2 ) )

 ! write(debugFile,'(" f(6)=",1pe10.2," curlfPttx1,nDotCurlfPtt1=",2(1pe10.2))') f(6),curlfPttx1,nDotCurlfPtt1
 ! write(debugFile,'(" f(7)=",1pe10.2," curlfPttx2,nDotCurlfPtt2=",2(1pe10.2))') f(7),curlfPttx2,nDotCurlfPtt2
 ! write(debugFile,'(" fPttx2=",3(1pe10.2))') fPttx2(0),fPttx2(1),fPttx2(2)
 ! write(debugFile,'(" fPtty2=",3(1pe10.2))') fPtty2(0),fPtty2(1),fPtty2(2)
 ! write(debugFile,'(" fPttz2=",3(1pe10.2))') fPttz2(0),fPttz2(1),fPttz2(2)
 ! write(debugFile,'(" divLapE1,curlLapE1x,nDotCurlLapE1=",3(1pe10.2))') divLapE1,curlLapE1x,nDotCurlLapE1
 ! write(debugFile,'(" divLapE2,curlLapE2x,nDotCurlLapE2=",3(1pe10.2))') divLapE2,curlLapE2x,nDotCurlLapE2


 ! [ (1/mu)* n n^T( c^2*Delta^2(E) - alphaP*Delta(P_tt) ) + (I-n n^T)( c^4 Delta^2(E) - c^2 alphaP*Delta(P_tt) - alphaP*P_tttt ]=0
 !  Note: fLptt = c^2*Delta( Ptt ) 
 ! f(9)= ( ( u1LapSq/(epsmu1) + cem1*nDotLapSqE1*an1 )/(epsmu1) - alphaP1*( (fLPtt1(0)-an1*nDotfLPtt1)/epsmu1 + fPtttt1(0)-an1*nDotfPtttt1  ) ) - \
 !       ( ( u2LapSq/(epsmu2) + cem2*nDotLapSqE2*an1 )/(epsmu2) - alphaP2*( (fLPtt2(0)-an1*nDotfLPtt2)/epsmu2 + fPtttt2(0)-an1*nDotfPtttt2  ) )
 !
 ! f(10)=( ( v1LapSq/(epsmu1) + cem1*nDotLapSqE1*an2 )/(epsmu1) - alphaP1*( (fLPtt1(1)-an2*nDotfLPtt1)/epsmu1 + fPtttt1(1)-an2*nDotfPtttt1  ) ) - \
 !       ( ( v2LapSq/(epsmu2) + cem2*nDotLapSqE2*an2 )/(epsmu2) - alphaP2*( (fLPtt2(1)-an2*nDotfLPtt2)/epsmu2 + fPtttt2(1)-an2*nDotfPtttt2  ) )
 !
 ! f(11)=( ( w1LapSq/(epsmu1) + cem1*nDotLapSqE1*an3 )/(epsmu1) - alphaP1*( (fLPtt1(2)-an3*nDotfLPtt1)/epsmu1 + fPtttt1(2)-an3*nDotfPtttt1  ) ) - \
 !       ( ( w2LapSq/(epsmu2) + cem2*nDotLapSqE2*an3 )/(epsmu2) - alphaP2*( (fLPtt2(2)-an3*nDotfLPtt2)/epsmu2 + fPtttt2(2)-an3*nDotfPtttt2  ) )

 ! AGAIN -- Oct 24, 2018:
 ! [ (1/mu)* n n^T( c^2 div(Delta(E) - alphaP* Delta(Ptt) )  + (I-n n^T)( c^4 Delta^2(E) - alphaP*{ c^2 *Delta(P_tt) + P_tttt} ]=0
 !  Note: fLptt = c^2*Delta( Ptt ) 
 nDotEP1 = (nDotLapSqE1/epsmu1 -alphaP1*nDotfLPtt1*epsmu1)/mu1
 nDotEP2 = (nDotLapSqE2/epsmu2 -alphaP2*nDotfLPtt2*epsmu2)/mu2

 f(9)= ( nDotEP1*an1 + (u1LapSq-nDotLapSqE1*an1)/epsmu1**2 -alphaP1*( (fLPtt1(0)-an1*nDotfLPtt1) + fPtttt1(0)-an1*nDotfPtttt1 ) ) - \
       ( nDotEP2*an1 + (u2LapSq-nDotLapSqE2*an1)/epsmu2**2 -alphaP2*( (fLPtt2(0)-an1*nDotfLPtt2) + fPtttt2(0)-an1*nDotfPtttt2 ) )
                                                                                             
 f(10)=( nDotEP1*an2 + (v1LapSq-nDotLapSqE1*an2)/epsmu1**2 -alphaP1*( (fLPtt1(1)-an2*nDotfLPtt1) + fPtttt1(1)-an2*nDotfPtttt1 ) ) - \
       ( nDotEP2*an2 + (v2LapSq-nDotLapSqE2*an2)/epsmu2**2 -alphaP2*( (fLPtt2(1)-an2*nDotfLPtt2) + fPtttt2(1)-an2*nDotfPtttt2 ) )
                                                                                             
 f(11)=( nDotEP1*an3 + (w1LapSq-nDotLapSqE1*an3)/epsmu1**2 -alphaP1*( (fLPtt1(2)-an3*nDotfLPtt1) + fPtttt1(2)-an3*nDotfPtttt1 ) ) - \
       ( nDotEP2*an3 + (w2LapSq-nDotLapSqE2*an3)/epsmu2**2 -alphaP2*( (fLPtt2(2)-an3*nDotfLPtt2) + fPtttt2(2)-an3*nDotfPtttt2 ) )

 ! For testing extrap the 2nd ghost line :
 ! *e678*
 ! f(6)= u1(i1-2*is1,i2-2*is2,i3-2*is3,ex) - extrap5(u1,i1-is1,i2-is2,i3-is3,ex,is1,is2,is3)
 ! f(7)= u1(i1-2*is1,i2-2*is2,i3-2*is3,ey) - extrap5(u1,i1-is1,i2-is2,i3-is3,ey,is1,is2,is3)
 ! f(8)= u1(i1-2*is1,i2-2*is2,i3-2*is3,ez) - extrap5(u1,i1-is1,i2-is2,i3-is3,ez,is1,is2,is3)   

 ! f(9)= u2(j1-2*js1,j2-2*js2,j3-2*js3,ex)-extrap5(u2,j1-js1,j2-js2,j3-js3,ex,js1,js2,js3)
 ! f(10)=u2(j1-2*js1,j2-2*js2,j3-2*js3,ey)-extrap5(u2,j1-js1,j2-js2,j3-js3,ey,js1,js2,js3)
 ! f(11)=u2(j1-2*js1,j2-2*js2,j3-2*js3,ez)-extrap5(u2,j1-js1,j2-js2,j3-js3,ez,js1,js2,js3)

 if( twilightZone.eq.1 )then

   x1=xy1(i1,i2,i3,0)
   y1=xy1(i1,i2,i3,1)
   z1=xy1(i1,i2,i3,2)
   call ogderiv(ep, 0,1,0,0, x1,y1,z1,t, ex, uex  )
   call ogderiv(ep, 0,0,1,0, x1,y1,z1,t, ex, uey  )
   call ogderiv(ep, 0,0,0,1, x1,y1,z1,t, ex, uez  )
   call ogderiv(ep, 0,2,0,0, x1,y1,z1,t, ex, uexx )
   call ogderiv(ep, 0,0,2,0, x1,y1,z1,t, ex, ueyy )
   call ogderiv(ep, 0,0,0,2, x1,y1,z1,t, ex, uezz )

   call ogderiv(ep, 0,1,0,0, x1,y1,z1,t, ey, vex  )
   call ogderiv(ep, 0,0,1,0, x1,y1,z1,t, ey, vey  )
   call ogderiv(ep, 0,0,0,1, x1,y1,z1,t, ey, vez  )
   call ogderiv(ep, 0,2,0,0, x1,y1,z1,t, ey, vexx )
   call ogderiv(ep, 0,0,2,0, x1,y1,z1,t, ey, veyy )
   call ogderiv(ep, 0,0,0,2, x1,y1,z1,t, ey, vezz )

   call ogderiv(ep, 0,1,0,0, x1,y1,z1,t, ez, wex  )
   call ogderiv(ep, 0,0,1,0, x1,y1,z1,t, ez, wey  )
   call ogderiv(ep, 0,0,0,1, x1,y1,z1,t, ez, wez  )
   call ogderiv(ep, 0,2,0,0, x1,y1,z1,t, ez, wexx )
   call ogderiv(ep, 0,0,2,0, x1,y1,z1,t, ez, weyy )
   call ogderiv(ep, 0,0,0,2, x1,y1,z1,t, ez, wezz )

   call ogderiv(ep, 0,3,0,0, x1,y1,z1,t, ex, uexxx )
   call ogderiv(ep, 0,0,3,0, x1,y1,z1,t, ex, ueyyy )
   call ogderiv(ep, 0,0,0,3, x1,y1,z1,t, ex, uezzz )
   call ogderiv(ep, 0,2,1,0, x1,y1,z1,t, ex, uexxy )
   call ogderiv(ep, 0,2,0,1, x1,y1,z1,t, ex, uexxz )
   call ogderiv(ep, 0,0,2,1, x1,y1,z1,t, ex, ueyyz )

   call ogderiv(ep, 0,1,2,0, x1,y1,z1,t, ex, uexyy ) 
   call ogderiv(ep, 0,1,0,2, x1,y1,z1,t, ex, uexzz )
   call ogderiv(ep, 0,0,1,2, x1,y1,z1,t, ex, ueyzz )


   call ogderiv(ep, 0,4,0,0, x1,y1,z1,t, ex, uexxxx )
   call ogderiv(ep, 0,2,2,0, x1,y1,z1,t, ex, uexxyy )
   call ogderiv(ep, 0,2,0,2, x1,y1,z1,t, ex, uexxzz )
   call ogderiv(ep, 0,0,4,0, x1,y1,z1,t, ex, ueyyyy )
   call ogderiv(ep, 0,0,2,2, x1,y1,z1,t, ex, ueyyzz )
   call ogderiv(ep, 0,0,0,4, x1,y1,z1,t, ex, uezzzz )

   call ogderiv(ep, 0,3,0,0, x1,y1,z1,t, ey, vexxx )
   call ogderiv(ep, 0,0,3,0, x1,y1,z1,t, ey, veyyy )
   call ogderiv(ep, 0,0,0,3, x1,y1,z1,t, ey, vezzz )
   call ogderiv(ep, 0,2,1,0, x1,y1,z1,t, ey, vexxy )
   call ogderiv(ep, 0,2,0,1, x1,y1,z1,t, ey, vexxz )
   call ogderiv(ep, 0,0,2,1, x1,y1,z1,t, ey, veyyz )

   call ogderiv(ep, 0,1,2,0, x1,y1,z1,t, ey, vexyy )
   call ogderiv(ep, 0,1,0,2, x1,y1,z1,t, ey, vexzz )
   call ogderiv(ep, 0,0,1,2, x1,y1,z1,t, ey, veyzz )

   call ogderiv(ep, 0,4,0,0, x1,y1,z1,t, ey, vexxxx )
   call ogderiv(ep, 0,2,2,0, x1,y1,z1,t, ey, vexxyy )
   call ogderiv(ep, 0,2,0,2, x1,y1,z1,t, ey, vexxzz )
   call ogderiv(ep, 0,0,4,0, x1,y1,z1,t, ey, veyyyy )
   call ogderiv(ep, 0,0,2,2, x1,y1,z1,t, ey, veyyzz )
   call ogderiv(ep, 0,0,0,4, x1,y1,z1,t, ey, vezzzz )

   call ogderiv(ep, 0,3,0,0, x1,y1,z1,t, ez, wexxx )
   call ogderiv(ep, 0,0,3,0, x1,y1,z1,t, ez, weyyy )
   call ogderiv(ep, 0,0,0,3, x1,y1,z1,t, ez, wezzz )
   call ogderiv(ep, 0,2,1,0, x1,y1,z1,t, ez, wexxy )
   call ogderiv(ep, 0,2,0,1, x1,y1,z1,t, ez, wexxz )
   call ogderiv(ep, 0,0,2,1, x1,y1,z1,t, ez, weyyz )

   call ogderiv(ep, 0,1,2,0, x1,y1,z1,t, ez, wexyy ) 
   call ogderiv(ep, 0,1,0,2, x1,y1,z1,t, ez, wexzz )
   call ogderiv(ep, 0,0,1,2, x1,y1,z1,t, ez, weyzz )

   call ogderiv(ep, 0,4,0,0, x1,y1,z1,t, ez, wexxxx )
   call ogderiv(ep, 0,2,2,0, x1,y1,z1,t, ez, wexxyy )
   call ogderiv(ep, 0,2,0,2, x1,y1,z1,t, ez, wexxzz )
   call ogderiv(ep, 0,0,4,0, x1,y1,z1,t, ez, weyyyy )
   call ogderiv(ep, 0,0,2,2, x1,y1,z1,t, ez, weyyzz )
   call ogderiv(ep, 0,0,0,4, x1,y1,z1,t, ez, wezzzz )


   ueLap = uexx+ueyy+uezz
   veLap = vexx+veyy+vezz
   weLap = wexx+weyy+wezz

   curlEex = wey-vez
   curlEey = uez-wex
   curlEez = vex-uey
   nDotCurlEe=an1*curlEex+an2*curlEey+an3*curlEez
   nDotLapEe=an1*ueLap+an2*veLap+an3*weLap

   f(0)= f(0) - ( (curlEex- nDotCurlEe*an1)*(1./mu1-1./mu2) )
   f(1)= f(1) - ( (curlEey- nDotCurlEe*an2)*(1./mu1-1./mu2) )
   f(2)= f(2) - ( (curlEez- nDotCurlEe*an3)*(1./mu1-1./mu2) )

   nDotPevttSum1 = an1*pevttSum1(0) + an2*pevttSum1(1) + an3*pevttSum1(2)
   nDotPevttSum2 = an1*pevttSum2(0) + an2*pevttSum2(1) + an3*pevttSum2(2)

   f(3)= f(3) - ( ueLap*(1./epsmu1-1./epsmu2) + nDotLapEe*an1*(cem1-cem2) ) \
              + alphaP1*( pevttSum1(0) -an1*nDotPevttSum1 ) \
              - alphaP2*( pevttSum2(0) -an1*nDotPevttSum2 ) 

   f(4)= f(4) - ( veLap*(1./epsmu1-1./epsmu2) + nDotLapEe*an2*(cem1-cem2) ) \
              + alphaP1*( pevttSum1(1) -an2*nDotPevttSum1 ) \
              - alphaP2*( pevttSum2(1) -an2*nDotPevttSum2 )

   f(5)= f(5) - ( weLap*(1./epsmu1-1./epsmu2) + nDotLapEe*an3*(cem1-cem2) )  \
              + alphaP1*( pevttSum1(2) -an3*nDotPevttSum1 ) \
              - alphaP2*( pevttSum2(2) -an3*nDotPevttSum2 )


   ueLapSq = uexxxx+ueyyyy+uezzzz +2.*( uexxyy+uexxzz+ueyyzz )
   veLapSq = vexxxx+veyyyy+vezzzz +2.*( vexxyy+vexxzz+veyyzz )
   weLapSq = wexxxx+weyyyy+wezzzz +2.*( wexxyy+wexxzz+weyyzz )

   curlLapEex = wexxy-vexxz + weyyy-veyyz + weyzz-vezzz
   curlLapEey = uexxz-wexxx + ueyyz-wexyy + uezzz-wexzz
   curlLapEez = vexxx-uexxy + vexyy-ueyyy + vexzz-ueyzz
   nDotCurlLapEe=an1*curlLapEex+an2*curlLapEey+an3*curlLapEez
   nDotLapSqEe=an1*ueLapSq+an2*veLapSq+an3*weLapSq

   curlPevttxSum1 = pevttySum1(2) - pevttzSum1(1)
   curlPevttySum1 = pevttzSum1(0) - pevttxSum1(2)
   curlPevttzSum1 = pevttxSum1(1) - pevttySum1(0)
   nDotCurlPevttSum1 = an1*curlPevttxSum1 + an2*curlPevttySum1 + an3*curlPevttzSum1

   curlPevttxSum2 = pevttySum2(2) - pevttzSum2(1)
   curlPevttySum2 = pevttzSum2(0) - pevttxSum2(2)
   curlPevttzSum2 = pevttxSum2(1) - pevttySum2(0)
   nDotCurlPevttSum2 = an1*curlPevttxSum2 + an2*curlPevttySum2 + an3*curlPevttzSum2

   f(6) = f(6) - ( (curlLapEex- nDotCurlLapEe*an1)*(1./(mu1*epsmu1)-1./(mu2*epsmu2)) )  \
               + (alphaP1/mu1)*( curlPevttxSum1 - an1*nDotCurlPevttSum1 ) \
               - (alphaP2/mu2)*( curlPevttxSum2 - an1*nDotCurlPevttSum2 ) 

   f(7) = f(7) - ( (curlLapEey- nDotCurlLapEe*an2)*(1./(mu1*epsmu1)-1./(mu2*epsmu2)) ) \
               + (alphaP1/mu1)*( curlPevttySum1 - an2*nDotCurlPevttSum1 ) \
               - (alphaP2/mu2)*( curlPevttySum2 - an2*nDotCurlPevttSum2 ) 

   f(8) = f(8) - ( (curlLapEez- nDotCurlLapEe*an3)*(1./(mu1*epsmu1)-1./(mu2*epsmu2)) ) \
               + (alphaP1/mu1)*( curlPevttzSum1 - an3*nDotCurlPevttSum1 ) \
               - (alphaP2/mu2)*( curlPevttzSum2 - an3*nDotCurlPevttSum2 ) 

   nDotPevttLSum1 = an1*pevttLSum1(0)  + an2*pevttLSum1(1)  + an3*pevttLSum1(2)
   nDotPevttttSum1= an1*pevttttSum1(0) + an2*pevttttSum1(1) + an3*pevttttSum1(2)

   nDotPevttLSum2 = an1*pevttLSum2(0)  + an2*pevttLSum2(1)  + an3*pevttLSum2(2)
   nDotPevttttSum2= an1*pevttttSum2(0) + an2*pevttttSum2(1) + an3*pevttttSum2(2)

   f(9) = f(9) - ( ueLapSq*(1./epsmu1**2-1./epsmu2**2) + nDotLapSqEe*an1*(cem1/epsmu1-cem2/epsmu2) ) \
               + alphaP1*(  (pevttLSum1(0)-an1*nDotPevttLSum1)/epsmu1 + (pevttttSum1(0)-an1*nDotPevttttSum1) ) \
               - alphaP2*(  (pevttLSum2(0)-an1*nDotPevttLSum2)/epsmu2 + (pevttttSum2(0)-an1*nDotPevttttSum2) )

   f(10)= f(10)- ( veLapSq*(1./epsmu1**2-1./epsmu2**2) + nDotLapSqEe*an2*(cem1/epsmu1-cem2/epsmu2) ) \
               + alphaP1*(  (pevttLSum1(1)-an2*nDotPevttLSum1)/epsmu1 + (pevttttSum1(1)-an2*nDotPevttttSum1) ) \
               - alphaP2*(  (pevttLSum2(1)-an2*nDotPevttLSum2)/epsmu2 + (pevttttSum2(1)-an2*nDotPevttttSum2) )

   f(11)= f(11)- ( weLapSq*(1./epsmu1**2-1./epsmu2**2) + nDotLapSqEe*an3*(cem1/epsmu1-cem2/epsmu2) )  \
               + alphaP1*(  (pevttLSum1(2)-an3*nDotPevttLSum1)/epsmu1 + (pevttttSum1(2)-an3*nDotPevttttSum1) ) \
               - alphaP2*(  (pevttLSum2(2)-an3*nDotPevttLSum2)/epsmu2 + (pevttttSum2(2)-an3*nDotPevttttSum2) )

 end if

#endMacro
 


! -------------------------------------------------------------------------
! Macro: Evaluate DISPERSIVE forcing terms, FOURTH-ORDER AND 3D
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
#beginMacro getDispersiveForcing3dOrder4(FACE,k1,k2,k3, fp, fpv,fev, p,pn,pm, u,un,um, dispersionModel,numberOfPolarizationVectors,alphaP,\
            c2PttEsum,c2PttLEsum,c4PttLEsum,c4PttLLEsum,c2PttttLEsum,c2PttttLLEsum,a0v,a1v,b0v,b1v,LE,LLE,LEm,LfE,LfP,fEt,fEtt,fPt,fPtt,pevtt,pevttx,pevtty,pevttz,pevtttt,evx,evy,evz,evnx,evny,evnz,\
            fevx,fevy,fevz,fpvx,fpvy,fpvz,LEx,LEy,LEz,fPttx,fPtty,fPttz,fLPtt,fPtttt)

 ! ********************** FINISH ME FOR 3D *******************

 do n=0,nd-1
   fp(n)=0.
   fPttx(n) =0.
   fPtty(n) =0.
   fPttz(n) =0.
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
   #Include interfaceAdeGdmOrder4.h 
 
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
       evalSecondDerivs3d(rsxy1,aj1,p1,k1,k2,k3,pc,uu1,p1)
       ! write(*,'("FACE: p1x,p1y,p1xx,p1yy,p1Lap=",5(1pe12.4))') p1x,p1y,p1xx,p1yy,p1Lap
 
       LP  = (c1**2)*p1Lap
       evalSecondDerivs3d(rsxy1,aj1,p1n,k1,k2,k3,pc,uu1,p1n)
       ! write(*,'("FACE: p1nxx,p1nyy,p1nLap=",3e12.4)') p1nxx,p1nyy,p1nLap
       LPm = (c1**2)*p1nLap
 
       pvx  = p1x
       pvy  = p1y
       pvz  = p1z

       pvnx = p1nx
       pvny = p1ny
       pvnz = p1nz
 
     #Elif #FACE eq "RIGHT"
       opEvalJacobianDerivatives(rsxy2,k1,k2,k3,aj2,1)
       ! uu1 in the next statement defines names of intermediate values
       evalSecondDerivs3d(rsxy2,aj2,p2,k1,k2,k3,pc,uu2,p2)
       LP  = (c2**2)*p2Lap
       evalSecondDerivs3d(rsxy2,aj2,p2n,k1,k2,k3,pc,uu2,p2n)
       LPm = (c2**2)*p2nLap
 
       pvx  = p2x
       pvy  = p2y
       pvz  = p2z

       pvnx = p2nx
       pvny = p2ny
       pvnz = p2nz
 
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
     pttza = c2PttLE*LEz(n) + c2PttE*evz(n) + c2PttEm*evnz(n) + c2PttP*pvz + c2PttPm*pvnz + c2PttfE*fevz(n) + c2PttfP*fpvz(n,jv)
 
  ! write(debugFile,'(" pttxa,fevx(n),fpvx(n,jv)=",3(1pe10.2))') pttxa,fevx(n),fpvx(n,jv)
  ! write(debugFile,'(" LEx(n),evx(n),evnx(n),pvx,pvnx=",5(1pe10.2))') LEx(n),evx(n),evnx(n),pvx,pvnx

     fPttx(n) = fPttx(n) + pttxa
     fPtty(n) = fPtty(n) + pttya
     fPttz(n) = fPttz(n) + pttza
 
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



! ==========================================================================================
!   Evaluate the jump conditions (including compatibility) for the GDM interface equations 
! ==========================================================================================
#beginMacro evaluateDispersiveInterfaceEquations3dOrder4()

 ! Evaluate TZ forcing for dispersive equations in 2=3D 
 getDispersiveTZForcing3dOrder4(fpv1,fpv2,fev1,fev2)

 evalDerivs3dOrder4()

 ! Store c^2*Delta(E) in a vector 
 LE1(0)=(c1**2)*u1Lap
 LE1(1)=(c1**2)*v1Lap
 LE1(2)=(c1**2)*w1Lap
 
 LE2(0)=(c2**2)*u2Lap
 LE2(1)=(c2**2)*v2Lap
 LE2(2)=(c2**2)*w2Lap
 
 ! Store L^2(E) 
 LLE1(0)=(c1**4)*u1LapSq
 LLE1(1)=(c1**4)*v1LapSq
 LLE1(2)=(c1**4)*w1LapSq
 
 LLE2(0)=(c2**4)*u2LapSq
 LLE2(1)=(c2**4)*v2LapSq
 LLE2(2)=(c2**4)*w2LapSq
 
 ! Store (LE).x an (LE).y 
 LEx1(0) = (c1**2)*( u1xxx + u1xyy + u1xzz )
 LEx1(1) = (c1**2)*( v1xxx + v1xyy + v1xzz )
 LEx1(2) = (c1**2)*( w1xxx + w1xyy + w1xzz )

 LEy1(0) = (c1**2)*( u1xxy + u1yyy + u1yzz )
 LEy1(1) = (c1**2)*( v1xxy + v1yyy + v1yzz )
 LEy1(2) = (c1**2)*( w1xxy + w1yyy + w1yzz )

 LEz1(0) = (c1**2)*( u1xxz + u1yyz + u1zzz )
 LEz1(1) = (c1**2)*( v1xxz + v1yyz + v1zzz )
 LEz1(2) = (c1**2)*( w1xxz + w1yyz + w1zzz )

 LEx2(0) = (c2**2)*( u2xxx + u2xyy + u2xzz )
 LEx2(1) = (c2**2)*( v2xxx + v2xyy + v2xzz )
 LEx2(2) = (c2**2)*( w2xxx + w2xyy + w2xzz )

 LEy2(0) = (c2**2)*( u2xxy + u2yyy + u2yzz )
 LEy2(1) = (c2**2)*( v2xxy + v2yyy + v2yzz )
 LEy2(2) = (c2**2)*( w2xxy + w2yyy + w2yzz )

 LEz2(0) = (c2**2)*( u2xxz + u2yyz + u2zzz )
 LEz2(1) = (c2**2)*( v2xxz + v2yyz + v2zzz )
 LEz2(2) = (c2**2)*( w2xxz + w2yyz + w2zzz )

 ! We also need derivatives at the old time:
 ! These next derivatives may only be needed to order2, but use order 4 for now so exact for degree 4
 #perl $ORDER=4;
 ! perl $ORDER=2;
 opEvalJacobianDerivatives(rsxy1,i1,i2,i3,aj1,1)
 evalSecondDerivs3d(rsxy1,aj1,u1n,i1,i2,i3,ex,uu1,u1n)
 evalSecondDerivs3d(rsxy1,aj1,u1n,i1,i2,i3,ey,vv1,v1n)
 evalSecondDerivs3d(rsxy1,aj1,u1n,i1,i2,i3,ez,ww1,w1n)
 ! Here are c^2*Delta(E) at the old time: 
 LE1m(0) = (c1**2)*u1nLap
 LE1m(1) = (c1**2)*v1nLap
 LE1m(2) = (c1**2)*w1nLap

 opEvalJacobianDerivatives(rsxy2,j1,j2,j3,aj2,1)
 evalSecondDerivs3d(rsxy2,aj2,u2n,j1,j2,j3,ex,uu2,u2n)
 evalSecondDerivs3d(rsxy2,aj2,u2n,j1,j2,j3,ey,vv2,v2n)
 evalSecondDerivs3d(rsxy2,aj2,u2n,j1,j2,j3,ez,ww2,w2n)
 LE2m(0) = (c2**2)*u2nLap
 LE2m(1) = (c2**2)*v2nLap
 LE2m(2) = (c2**2)*w2nLap

 evx1(0) = u1x
 evx1(1) = v1x
 evx1(2) = w1x

 evy1(0) = u1y
 evy1(1) = v1y 
 evy1(2) = w1y 

 evz1(0) = u1z
 evz1(1) = v1z 
 evz1(2) = w1z 

 evnx1(0) = u1nx
 evnx1(1) = v1nx
 evnx1(2) = w1nx

 evny1(0) = u1ny
 evny1(1) = v1ny 
 evny1(2) = w1ny 

 evnz1(0) = u1nz
 evnz1(1) = v1nz 
 evnz1(2) = w1nz 

 evx2(0) = u2x
 evx2(1) = v2x
 evx2(2) = w2x

 evy2(0) = u2y
 evy2(1) = v2y 
 evy2(2) = w2y 

 evz2(0) = u2z
 evz2(1) = v2z 
 evz2(2) = w2z 

 evnx2(0) = u2nx
 evnx2(1) = v2nx
 evnx2(2) = w2nx

 evny2(0) = u2ny
 evny2(1) = v2ny 
 evny2(2) = w2ny 

 evnz2(0) = u2nz
 evnz2(1) = v2nz 
 evnz2(2) = w2nz 

 ! eval dispersive forcings for domain 1
 getDispersiveForcing3dOrder4(LEFT,i1,i2,i3, fp1, fpv1,fev1,p1,p1n,p1m, u1,u1n,u1m, dispersionModel1,\
    numberOfPolarizationVectors1,alphaP1,c2PttEsum1,c2PttLEsum1,c4PttLEsum1,c4PttLLEsum1,c2PttttLEsum1,c2PttttLLEsum1,\
    a0v1,a1v1,b0v1,b1v1,LE1,LLE1,LE1m,LfE1,LfP1,fEt1,fEtt1,fPt1,fPtt1,pevtt1,pevttx1,pevtty1,pevttz1,pevtttt1,\
    evx1,evy1,evz1,evnx1,evny1,evnz1,fevx1,fevy1,fevz1,fpvx1,fpvy1,fpvz1,LEx1,LEy1,LEz1,fPttx1,fPtty1,fPttz1,fLPtt1,fPtttt1) 

 ! eval dispersive forcings for domain 2
 getDispersiveForcing3dOrder4(RIGHT,j1,j2,j3, fp2, fpv2,fev2,p2,p2n,p2m, u2,u2n,u2m, dispersionModel2,\
    numberOfPolarizationVectors2,alphaP2,c2PttEsum2,c2PttLEsum2,c4PttLEsum2,c4PttLLEsum2,c2PttttLEsum2,c2PttttLLEsum2,\
    a0v2,a1v2,b0v2,b1v2,LE2,LLE2,LE2m,LfE2,LfP2,fEt2,fEtt2,fPt2,fPtt2,pevtt2,pevttx2,pevtty2,pevttz2,pevtttt2,\
    evx2,evy2,evz2,evnx2,evny2,evnz2,fevx2,fevy2,fevz2,fpvx2,fpvy2,fpvz2,LEx2,LEy2,LEz2,fPttx2,fPtty2,fPttz2,fLPtt2,fPtttt2) 


 eval3dJumpDispersiveOrder4()

 ! eval3dJumpOrder4New()
#endMacro


! ==========================================================================================
!         Fill in FIRST SET of 6 interface equations -- GDM -- .
! 
! Input:
!  e0,e1,e2,e3,e4,e5 : equation numbers (0,1,2,3,4,5) or (6,7,8,9,10,11)
! ==========================================================================================
#beginMacro fillDispersiveEquations3dOrder4a(e0,e1,e2,e3,e4,e5)

 ! Equations 0,1,2:
 ! [ div(E) n  + (1/mu)*(I-n n^T)( curl(E) ] = 0 

 ! Equation 0: 
 ! (u.x+v.y+w.z)*an1 + ( w1y-v1z - nDotCurlE1*an1)/mu1

 ! n^T curl(E) = an1*((Ez)_y - (Ey)_z+ + an2*((Ex)_z - (Ez)_x) + an3*((Ey)_x - (Ex)_y )
 ! (I - n n^T) curl(E) = [ (Ez)_y - (Ey)_z - an1*( n^T curl(E) )]   [ a11*Ex + a12*Ey + a13*Ez ]
 !                       [ (Ex)_z - (Ez)_x - an2*(             )] = [ a21*Ex + a22*Ey + a23*Ez ]
 !                       [ (Ey)_x - (Ex)_y - an3*(             )]   [ a31*Ex + a32*Et + a33*Ez ]
 !       a11 =     -an1*( an2*Dz -an3*Dy )
 !       a12 = -Dz -an1*( an3*Dx -an1*Dz )
 !       a13 =  Dy -an1*( an1*Dy -an2*Dz )

 ! ------- Equation 0 ------
 ! coeffs of u1(-1), v1(-1), w1(-1)
 a12(e0,0) = ( c1m1x*an1 + (           - (c1m1z*an2-c1m1y*an3)*an1 )/mu1 ) ! coeff of u1(-1)
 a12(e0,1) = ( c1m1y*an1 + (    -c1m1z - (c1m1x*an3-c1m1z*an1)*an1 )/mu1 ) ! coeff of v1(-1)
 a12(e0,2) = ( c1m1z*an1 + ( c1m1y     - (c1m1y*an1-c1m1x*an2)*an1 )/mu1 ) ! coeff of w1(-1)

 ! coeffs of u1(-2), v1(-2), w1(-2)
 a12(e0,3) = ( c1m2x*an1 + (           - (c1m2z*an2-c1m2y*an3)*an1 )/mu1 ) ! coeff of u1(-1)
 a12(e0,4) = ( c1m2y*an1 + (    -c1m2z - (c1m2x*an3-c1m2z*an1)*an1 )/mu1 ) ! coeff of v1(-1)
 a12(e0,5) = ( c1m2z*an1 + ( c1m2y     - (c1m2y*an1-c1m2x*an2)*an1 )/mu1 ) ! coeff of w1(-1)

 ! coeffs of u2(-1), v2(-1), w2(-1) 
 a12(e0,6) =-( c2m1x*an1 + (           - (c2m1z*an2-c2m1y*an3)*an1 )/mu2 ) ! coeff of u2(-1)
 a12(e0,7) =-( c2m1y*an1 + (    -c2m1z - (c2m1x*an3-c2m1z*an1)*an1 )/mu2 ) ! coeff of v2(-1)
 a12(e0,8) =-( c2m1z*an1 + ( c2m1y     - (c2m1y*an1-c2m1x*an2)*an1 )/mu2 ) ! coeff of w2(-1)

 ! coeffs of u2(-2), v2(-2), w2(-2)
 a12(e0,9) =-( c2m2x*an1 + (           - (c2m2z*an2-c2m2y*an3)*an1 )/mu2 ) ! coeff of u2(-1)
 a12(e0,10)=-( c2m2y*an1 + (    -c2m2z - (c2m2x*an3-c2m2z*an1)*an1 )/mu2 ) ! coeff of v2(-1)
 a12(e0,11)=-( c2m2z*an1 + ( c2m2y     - (c2m2y*an1-c2m2x*an2)*an1 )/mu2 ) ! coeff of w2(-1)

 ! *bug* fixed wdh July 6, 2019 "+" to minus "-" 8 places 
 ! Equation 1:
 ! (u.x+v.y+w.z)*an2 + ( u1z-w1x - nDotCurlE1*an2)/mu1
 a12(e1,0) = ( c1m1x*an2 + ( c1m1z     - (c1m1z*an2-c1m1y*an3)*an2 )/mu1 ) ! coeff of u1(-1)
 a12(e1,1) = ( c1m1y*an2 + (           - (c1m1x*an3-c1m1z*an1)*an2 )/mu1 ) ! coeff of v1(-1)
 a12(e1,2) = ( c1m1z*an2 + (    -c1m1x - (c1m1y*an1-c1m1x*an2)*an2 )/mu1 ) ! coeff of w1(-1)

 a12(e1,3) = ( c1m2x*an2 + ( c1m2z     - (c1m2z*an2-c1m2y*an3)*an2 )/mu1 ) ! coeff of u1(-1)
 a12(e1,4) = ( c1m2y*an2 + (           - (c1m2x*an3-c1m2z*an1)*an2 )/mu1 ) ! coeff of v1(-1)
 a12(e1,5) = ( c1m2z*an2 + (    -c1m2x - (c1m2y*an1-c1m2x*an2)*an2 )/mu1 ) ! coeff of w1(-1)

 a12(e1,6) =-( c2m1x*an2 + ( c2m1z     - (c2m1z*an2-c2m1y*an3)*an2 )/mu2 ) ! coeff of u2(-1)
 a12(e1,7) =-( c2m1y*an2 + (           - (c2m1x*an3-c2m1z*an1)*an2 )/mu2 ) ! coeff of v2(-1)
 a12(e1,8) =-( c2m1z*an2 + (    -c2m1x - (c2m1y*an1-c2m1x*an2)*an2 )/mu2 ) ! coeff of w2(-1)

 a12(e1,9) =-( c2m2x*an2 + ( c2m2z     - (c2m2z*an2-c2m2y*an3)*an2 )/mu2 ) ! coeff of u2(-1)
 a12(e1,10)=-( c2m2y*an2 + (           - (c2m2x*an3-c2m2z*an1)*an2 )/mu2 ) ! coeff of v2(-1)
 a12(e1,11)=-( c2m2z*an2 + (    -c2m2x - (c2m2y*an1-c2m2x*an2)*an2 )/mu2 ) ! coeff of w2(-1)

 ! Equation 2:
 ! (u.x+v.y+w.z)*an3 + ( v1x-u1y - nDotCurlE1*an3)/mu1
 a12(e2,0) = ( c1m1x*an3 + (    -c1m1y - (c1m1z*an2-c1m1y*an3)*an3 )/mu1 ) ! coeff of u1(-1)
 a12(e2,1) = ( c1m1y*an3 + ( c1m1x     - (c1m1x*an3-c1m1z*an1)*an3 )/mu1 ) ! coeff of v1(-1)
 a12(e2,2) = ( c1m1z*an3 + (           - (c1m1y*an1-c1m1x*an2)*an3 )/mu1 ) ! coeff of w1(-1)

 a12(e2,3) = ( c1m2x*an3 + (    -c1m2y - (c1m2z*an2-c1m2y*an3)*an3 )/mu1 ) ! coeff of u1(-1)
 a12(e2,4) = ( c1m2y*an3 + ( c1m2x     - (c1m2x*an3-c1m2z*an1)*an3 )/mu1 ) ! coeff of v1(-1)
 a12(e2,5) = ( c1m2z*an3 + (           - (c1m2y*an1-c1m2x*an2)*an3 )/mu1 ) ! coeff of w1(-1)

 a12(e2,6) =-( c2m1x*an3 + (    -c2m1y - (c2m1z*an2-c2m1y*an3)*an3 )/mu2 ) ! coeff of u2(-1)
 a12(e2,7) =-( c2m1y*an3 + ( c2m1x     - (c2m1x*an3-c2m1z*an1)*an3 )/mu2 ) ! coeff of v2(-1)
 a12(e2,8) =-( c2m1z*an3 + (           - (c2m1y*an1-c2m1x*an2)*an3 )/mu2 ) ! coeff of w2(-1)

 a12(e2,9) =-( c2m2x*an3 + (    -c2m2y - (c2m2z*an2-c2m2y*an3)*an3 )/mu2 ) ! coeff of u2(-1)
 a12(e2,10)=-( c2m2y*an3 + ( c2m2x     - (c2m2x*an3-c2m2z*an1)*an3 )/mu2 ) ! coeff of v2(-1)
 a12(e2,11)=-( c2m2z*an3 + (           - (c2m2y*an1-c2m2x*an2)*an3 )/mu2 ) ! coeff of w2(-1)

 !  (1/mu1)* u1Lap*ani*anj + (delta_ij -ani*anj)*( u1Lap/epsmu1 - alphaP1*( c4PttLEsum1*u1Lap + c4PttLLEsum1*u1LapSq )= ...

 ! ------------- Equation 3 -------------
 !  u1Lap/(epsmu1) + cem1*( an1*u1Lap + an2*v1Lap + an3*w1Lap )*an1 + (delta_ij -an1*anj)*( - alphaP1*( c4PttLEsum1*u1Lap + c4PttLLEsum1*u1LapSq ) )

 eqn1Coeffm1 = -alphaP1*( (c4PttLEsum1/epsmu1)*clap1m1+(c4PttLLEsum1/epsmu1**2)*clapSq1m1 )
 eqn1Coeffm2 = -alphaP1*( (c4PttLEsum1/epsmu1)*clap1m2+(c4PttLLEsum1/epsmu1**2)*clapSq1m2 )

 eqn2Coeffm1 = -alphaP2*( (c4PttLEsum2/epsmu2)*clap2m1+(c4PttLLEsum2/epsmu2**2)*clapSq2m1 )
 eqn2Coeffm2 = -alphaP2*( (c4PttLEsum2/epsmu2)*clap2m2+(c4PttLLEsum2/epsmu2**2)*clapSq2m2 )

 ! coeffs of u1(-1), v1(-1), w1(-1)
 a12(e3,0) = ( clap1m1/(epsmu1) + cem1*( an1*clap1m1                         )*an1 + (1.-an1*an1)*( eqn1Coeffm1 ) )
 a12(e3,1) = (                    cem1*(             an2*clap1m1             )*an1 + (  -an1*an2)*( eqn1Coeffm1 ) )
 a12(e3,2) = (                    cem1*(                         an3*clap1m1 )*an1 + (  -an1*an3)*( eqn1Coeffm1 ) )
                                                                                                                   
 ! coeffs of u1(-2), v1(-2), w1(-2)                                                                                
 a12(e3,3) = ( clap1m2/(epsmu1) + cem1*( an1*clap1m2                         )*an1 + (1.-an1*an1)*( eqn1Coeffm2 ) )
 a12(e3,4) = (                    cem1*(             an2*clap1m2             )*an1 + (  -an1*an2)*( eqn1Coeffm2 ) )
 a12(e3,5) = (                    cem1*(                         an3*clap1m2 )*an1 + (  -an1*an3)*( eqn1Coeffm2 ) )
                                                                                                                   
 ! coeffs of u2(-1), v2(-1), w2(-1)                                                                                
 a12(e3,6) =-( clap2m1/(epsmu2) + cem2*( an1*clap2m1                         )*an1 + (1.-an1*an1)*( eqn2Coeffm1 ) )
 a12(e3,7) =-(                    cem2*(             an2*clap2m1             )*an1 + (  -an1*an2)*( eqn2Coeffm1 ) )
 a12(e3,8) =-(                    cem2*(                         an3*clap2m1 )*an1 + (  -an1*an3)*( eqn2Coeffm1 ) )
                                                                                                                   
 ! coeffs of u2(-2), v2(-2), w2(-2)                                                                                
 a12(e3,9) =-( clap2m2/(epsmu2) + cem2*( an1*clap2m2                         )*an1 + (1.-an1*an1)*( eqn2Coeffm2 ) ) 
 a12(e3,10)=-(                    cem2*(             an2*clap2m2             )*an1 + (  -an1*an2)*( eqn2Coeffm2 ) )
 a12(e3,11)=-(                    cem2*(                         an3*clap2m2 )*an1 + (  -an1*an3)*( eqn2Coeffm2 ) )

 ! Equations 4,5,6
 !  [ (1/mu)* Lap(E) n n^T + (I-n n^T)* c^2 Lap(E) - alphaP*Ptt ) ] = 0 

 ! -------------- Equation 4 -----------------
 !  v1Lap/(epsmu1) + cem1*( an1*u1Lap + an2*v1Lap + an3*w1Lap )*an2+ (delta_ij -an2*anj)*( - alphaP1*( c4PttLEsum1*u1Lap + c4PttLLEsum1*u1LapSq ) )

 ! coeffs of u1(-1), v1(-1), w1(-1)
 a12(e4,0) = (                    cem1*( an1*clap1m1                         )*an2 + (  -an2*an1)*( eqn1Coeffm1 ) )  
 a12(e4,1) = ( clap1m1/(epsmu1) + cem1*(             an2*clap1m1             )*an2 + (1.-an2*an2)*( eqn1Coeffm1 ) )
 a12(e4,2) = (                    cem1*(                         an3*clap1m1 )*an2 + (  -an2*an3)*( eqn1Coeffm1 ) )
                                                                                                                   
 ! coeffs of u1(-2), v1(-2), w1(-2)                                                                                
 a12(e4,3) = (                    cem1*( an1*clap1m2                         )*an2 + (  -an2*an1)*( eqn1Coeffm2 ) ) 
 a12(e4,4) = ( clap1m2/(epsmu1) + cem1*(             an2*clap1m2             )*an2 + (1.-an2*an2)*( eqn1Coeffm2 ) )
 a12(e4,5) = (                    cem1*(                         an3*clap1m2 )*an2 + (  -an2*an3)*( eqn1Coeffm2 ) )
                                                                                                                   
 ! coeffs of u2(-1), v2(-1), w2(-1)                                                                                
 a12(e4,6) =-(                    cem2*( an1*clap2m1                         )*an2 + (  -an2*an1)*( eqn2Coeffm1 ) ) 
 a12(e4,7) =-( clap2m1/(epsmu2) + cem2*(             an2*clap2m1             )*an2 + (1.-an2*an2)*( eqn2Coeffm1 ) )
 a12(e4,8) =-(                    cem2*(                         an3*clap2m1 )*an2 + (  -an2*an3)*( eqn2Coeffm1 ) )
                                                                                                                   
 ! coeffs of u2(-2), v2(-2), w2(-2)                                                                                
 a12(e4,9) =-(                    cem2*( an1*clap2m2                         )*an2 + (  -an2*an1)*( eqn2Coeffm2 ) ) 
 a12(e4,10)=-( clap2m2/(epsmu2) + cem2*(             an2*clap2m2             )*an2 + (1.-an2*an2)*( eqn2Coeffm2 ) )
 a12(e4,11)=-(                    cem2*(                         an3*clap2m2 )*an2 + (  -an2*an3)*( eqn2Coeffm2 ) )

 ! ----- Equation 5 -----------

 ! coeffs of u1(-1), v1(-1), w1(-1)
 a12(e5,0) = (                    cem1*( an1*clap1m1                         )*an3 + (  -an3*an1)*( eqn1Coeffm1 ) )
 a12(e5,1) = (                    cem1*(             an2*clap1m1             )*an3 + (  -an3*an2)*( eqn1Coeffm1 ) )
 a12(e5,2) = ( clap1m1/(epsmu1) + cem1*(                         an3*clap1m1 )*an3 + (1.-an3*an3)*( eqn1Coeffm1 ) )
                                                                                                                   
 ! coeffs of u1(-2), v1(-2), w1(-2)                                                                                
 a12(e5,3) = (                    cem1*( an1*clap1m2                         )*an3 + (  -an3*an1)*( eqn1Coeffm2 ) )
 a12(e5,4) = (                    cem1*(             an2*clap1m2             )*an3 + (  -an3*an2)*( eqn1Coeffm2 ) )
 a12(e5,5) = ( clap1m2/(epsmu1) + cem1*(                         an3*clap1m2 )*an3 + (1.-an3*an3)*( eqn1Coeffm2 ) )
                                                                                                                   
 ! coeffs of u2(-1), v2(-1), w2(-1)                                                                                
 a12(e5,6) =-(                    cem2*( an1*clap2m1                         )*an3 + (  -an3*an1)*( eqn2Coeffm1 ) )
 a12(e5,7) =-(                    cem2*(             an2*clap2m1             )*an3 + (  -an3*an2)*( eqn2Coeffm1 ) )
 a12(e5,8) =-( clap2m1/(epsmu2) + cem2*(                         an3*clap2m1 )*an3 + (1.-an3*an3)*( eqn2Coeffm1 ) )
                                                                                                                   
 ! coeffs of u2(-2), v2(-2), w2(-2)                                                                                
 a12(e5,9) =-(                    cem2*( an1*clap2m2                         )*an3 + (  -an3*an1)*( eqn2Coeffm2 ) )
 a12(e5,10)=-(                    cem2*(             an2*clap2m2             )*an3 + (  -an3*an2)*( eqn2Coeffm2 ) )
 a12(e5,11)=-( clap2m2/(epsmu2) + cem2*(                         an3*clap2m2 )*an3 + (1.-an3*an3)*( eqn2Coeffm2 ) )
#endMacro



! ==========================================================================================
!         Fill in SECOND SET of 6 interface equations -- GDM -- .
! 
! Input:
!  e0,e1,e2,e3,e4,e5 : equation numbers  (6,7,8,9,10,11)
! ==========================================================================================
#beginMacro fillDispersiveEquations3dOrder4b(e0,e1,e2,e3,e4,e5)

 ! Equations 0,1,2:
 ! [ c^2 Lap(div(E)) n + (1/mu)*(I-n n^T)( c^2 Lap(curl(E)) - alphaP*curl(Ptt) ) ] = 0 
 ! curl(E) = [ (Ez)_y - (Ey)_z,  (Ex)_z - (Ez)_x , (Ey)_x - (Ex)_y ]

 ! n^T curl(E) = an1*((Ez)_y - (Ey)_z+ + an2*((Ex)_z - (Ez)_x) + an3*((Ey)_x - (Ex)_y )
 ! (I - n n^T) curl(E) = [ (Ez)_y - (Ey)_z - an1*( n^T curl(E) )]   [ a11*Ex + a12*Ey + a13*Ez ]
 !                       [ (Ex)_z - (Ez)_x - an2*(             )] = [ a21*Ex + a22*Ey + a23*Ez ]
 !                       [ (Ey)_x - (Ex)_y - an3*(             )]   [ a31*Ex + a32*Et + a33*Ez ]
 !       a11 =     -an1*( an2*Dz -an3*Dy )
 !       a12 = -Dz -an1*( an3*Dx -an1*Dz )
 !       a13 =  Dy -an1*( an1*Dy -an2*Dz )

 eqn1Coeff  = (1./epsmu1 -alphaP1*c2PttLEsum1/epsmu1)/mu1
 eqn1Coeffb = -alphaP1*c2PttEsum1/mu1 

 eqn2Coeff  = (1./epsmu2 -alphaP2*c2PttLEsum2/epsmu2)/mu2 
 eqn2Coeffb = -alphaP2*c2PttEsum2/mu2 

 ! ------- Equation 6 ------
 ! coeffs of u1(-1), v1(-1), w1(-1)
 a12(e0,0) = ( clap1m1x*an1/epsmu1 + (              - (clap1m1z*an2-clap1m1y*an3)*an1 )*eqn1Coeff + (           - (c1m1z*an2-c1m1y*an3)*an1 )*eqn1Coeffb )
 a12(e0,1) = ( clap1m1y*an1/epsmu1 + (    -clap1m1z - (clap1m1x*an3-clap1m1z*an1)*an1 )*eqn1Coeff + (    -c1m1z - (c1m1x*an3-c1m1z*an1)*an1 )*eqn1Coeffb )
 a12(e0,2) = ( clap1m1z*an1/epsmu1 + ( clap1m1y     - (clap1m1y*an1-clap1m1x*an2)*an1 )*eqn1Coeff + ( c1m1y     - (c1m1y*an1-c1m1x*an2)*an1 )*eqn1Coeffb )

 ! coeffs of u1(-2), v1(-2), w1(-2)
 a12(e0,3) = ( clap1m2x*an1/epsmu1 + (              - (clap1m2z*an2-clap1m2y*an3)*an1 )*eqn1Coeff + (           - (c1m2z*an2-c1m2y*an3)*an1 )*eqn1Coeffb ) 
 a12(e0,4) = ( clap1m2y*an1/epsmu1 + (    -clap1m2z - (clap1m2x*an3-clap1m2z*an1)*an1 )*eqn1Coeff + (    -c1m2z - (c1m2x*an3-c1m2z*an1)*an1 )*eqn1Coeffb ) 
 a12(e0,5) = ( clap1m2z*an1/epsmu1 + ( clap1m2y     - (clap1m2y*an1-clap1m2x*an2)*an1 )*eqn1Coeff + ( c1m2y     - (c1m2y*an1-c1m2x*an2)*an1 )*eqn1Coeffb ) 

 ! coeffs of u2(-1), v2(-1), w2(-1) 
 a12(e0,6) =-( clap2m1x*an1/epsmu2 + (              - (clap2m1z*an2-clap2m1y*an3)*an1 )*eqn2Coeff + (           - (c2m1z*an2-c2m1y*an3)*an1 )*eqn2Coeffb )
 a12(e0,7) =-( clap2m1y*an1/epsmu2 + (    -clap2m1z - (clap2m1x*an3-clap2m1z*an1)*an1 )*eqn2Coeff + (    -c2m1z - (c2m1x*an3-c2m1z*an1)*an1 )*eqn2Coeffb )
 a12(e0,8) =-( clap2m1z*an1/epsmu2 + ( clap2m1y     - (clap2m1y*an1-clap2m1x*an2)*an1 )*eqn2Coeff + ( c2m1y     - (c2m1y*an1-c2m1x*an2)*an1 )*eqn2Coeffb )

 ! coeffs of u2(-2), v2(-2), w2(-2)
 a12(e0,9) =-( clap2m2x*an1/epsmu2 + (              - (clap2m2z*an2-clap2m2y*an3)*an1 )*eqn2Coeff + (           - (c2m2z*an2-c2m2y*an3)*an1 )*eqn2Coeffb ) 
 a12(e0,10)=-( clap2m2y*an1/epsmu2 + (    -clap2m2z - (clap2m2x*an3-clap2m2z*an1)*an1 )*eqn2Coeff + (    -c2m2z - (c2m2x*an3-c2m2z*an1)*an1 )*eqn2Coeffb ) 
 a12(e0,11)=-( clap2m2z*an1/epsmu2 + ( clap2m2y     - (clap2m2y*an1-clap2m2x*an2)*an1 )*eqn2Coeff + ( c2m2y     - (c2m2y*an1-c2m2x*an2)*an1 )*eqn2Coeffb ) 

 ! ----- Equation 7 -----
 ! coeffs of u1(-1), v1(-1), w1(-1)
 a12(e1,0) = ( clap1m1x*an2/epsmu1 + ( clap1m1z     - (clap1m1z*an2-clap1m1y*an3)*an2 )*eqn1Coeff + ( c1m1z     - (c1m1z*an2-c1m1y*an3)*an2 )*eqn1Coeffb ) 
 a12(e1,1) = ( clap1m1y*an2/epsmu1 + (              - (clap1m1x*an3-clap1m1z*an1)*an2 )*eqn1Coeff + (           - (c1m1x*an3-c1m1z*an1)*an2 )*eqn1Coeffb ) 
 a12(e1,2) = ( clap1m1z*an2/epsmu1 + (    -clap1m1x - (clap1m1y*an1-clap1m1x*an2)*an2 )*eqn1Coeff + (    -c1m1x - (c1m1y*an1-c1m1x*an2)*an2 )*eqn1Coeffb ) 

 ! coeffs of u1(-2), v1(-2), w1(-2)
 a12(e1,3) = ( clap1m2x*an2/epsmu1 + ( clap1m2z     - (clap1m2z*an2-clap1m2y*an3)*an2 )*eqn1Coeff + ( c1m2z     - (c1m2z*an2-c1m2y*an3)*an2 )*eqn1Coeffb ) 
 a12(e1,4) = ( clap1m2y*an2/epsmu1 + (              - (clap1m2x*an3-clap1m2z*an1)*an2 )*eqn1Coeff + (           - (c1m2x*an3-c1m2z*an1)*an2 )*eqn1Coeffb ) 
 a12(e1,5) = ( clap1m2z*an2/epsmu1 + (    -clap1m2x - (clap1m2y*an1-clap1m2x*an2)*an2 )*eqn1Coeff + (    -c1m2x - (c1m2y*an1-c1m2x*an2)*an2 )*eqn1Coeffb ) 

 ! coeffs of u2(-1), v2(-1), w2(-1) 
 a12(e1,6) =-( clap2m1x*an2/epsmu2 + ( clap2m1z     - (clap2m1z*an2-clap2m1y*an3)*an2 )*eqn2Coeff + ( c2m1z     - (c2m1z*an2-c2m1y*an3)*an2 )*eqn2Coeffb ) 
 a12(e1,7) =-( clap2m1y*an2/epsmu2 + (              - (clap2m1x*an3-clap2m1z*an1)*an2 )*eqn2Coeff + (           - (c2m1x*an3-c2m1z*an1)*an2 )*eqn2Coeffb ) 
 a12(e1,8) =-( clap2m1z*an2/epsmu2 + (    -clap2m1x - (clap2m1y*an1-clap2m1x*an2)*an2 )*eqn2Coeff + (    -c2m1x - (c2m1y*an1-c2m1x*an2)*an2 )*eqn2Coeffb ) 

 ! coeffs of u2(-2), v2(-2), w2(-2)
 a12(e1,9) =-( clap2m2x*an2/epsmu2 + ( clap2m2z     - (clap2m2z*an2-clap2m2y*an3)*an2 )*eqn2Coeff + ( c2m2z     - (c2m2z*an2-c2m2y*an3)*an2 )*eqn2Coeffb ) 
 a12(e1,10)=-( clap2m2y*an2/epsmu2 + (              - (clap2m2x*an3-clap2m2z*an1)*an2 )*eqn2Coeff + (           - (c2m2x*an3-c2m2z*an1)*an2 )*eqn2Coeffb ) 
 a12(e1,11)=-( clap2m2z*an2/epsmu2 + (    -clap2m2x - (clap2m2y*an1-clap2m2x*an2)*an2 )*eqn2Coeff + (    -c2m2x - (c2m2y*an1-c2m2x*an2)*an2 )*eqn2Coeffb ) 

 ! ----- Equation 8 -----

 ! coeffs of u1(-1), v1(-1), w1(-1)
 a12(e2,0) = ( clap1m1x*an3/epsmu1 + (    -clap1m1y - (clap1m1z*an2-clap1m1y*an3)*an3 )*eqn1Coeff + (    -c1m1y - (c1m1z*an2-c1m1y*an3)*an3 )*eqn1Coeffb ) 
 a12(e2,1) = ( clap1m1y*an3/epsmu1 + ( clap1m1x     - (clap1m1x*an3-clap1m1z*an1)*an3 )*eqn1Coeff + ( c1m1x     - (c1m1x*an3-c1m1z*an1)*an3 )*eqn1Coeffb ) 
 a12(e2,2) = ( clap1m1z*an3/epsmu1 + (              - (clap1m1y*an1-clap1m1x*an2)*an3 )*eqn1Coeff + (           - (c1m1y*an1-c1m1x*an2)*an3 )*eqn1Coeffb ) 

 ! coeffs of u1(-2), v1(-2), w1(-2)
 a12(e2,3) = ( clap1m2x*an3/epsmu1 + (    -clap1m2y - (clap1m2z*an2-clap1m2y*an3)*an3 )*eqn1Coeff + (    -c1m2y - (c1m2z*an2-c1m2y*an3)*an3 )*eqn1Coeffb ) 
 a12(e2,4) = ( clap1m2y*an3/epsmu1 + ( clap1m2x     - (clap1m2x*an3-clap1m2z*an1)*an3 )*eqn1Coeff + ( c1m2x     - (c1m2x*an3-c1m2z*an1)*an3 )*eqn1Coeffb ) 
 a12(e2,5) = ( clap1m2z*an3/epsmu1 + (              - (clap1m2y*an1-clap1m2x*an2)*an3 )*eqn1Coeff + (           - (c1m2y*an1-c1m2x*an2)*an3 )*eqn1Coeffb ) 

 ! coeffs of u2(-1), v2(-1), w2(-1) 
 a12(e2,6) =-( clap2m1x*an3/epsmu2 + (    -clap2m1y - (clap2m1z*an2-clap2m1y*an3)*an3 )*eqn2Coeff + (    -c2m1y - (c2m1z*an2-c2m1y*an3)*an3 )*eqn2Coeffb ) 
 a12(e2,7) =-( clap2m1y*an3/epsmu2 + ( clap2m1x     - (clap2m1x*an3-clap2m1z*an1)*an3 )*eqn2Coeff + ( c2m1x     - (c2m1x*an3-c2m1z*an1)*an3 )*eqn2Coeffb ) 
 a12(e2,8) =-( clap2m1z*an3/epsmu2 + (              - (clap2m1y*an1-clap2m1x*an2)*an3 )*eqn2Coeff + (           - (c2m1y*an1-c2m1x*an2)*an3 )*eqn2Coeffb ) 

 ! coeffs of u2(-2), v2(-2), w2(-2)
 a12(e2,9) =-( clap2m2x*an3/epsmu2 + (    -clap2m2y - (clap2m2z*an2-clap2m2y*an3)*an3 )*eqn2Coeff + (    -c2m2y - (c2m2z*an2-c2m2y*an3)*an3 )*eqn2Coeffb ) 
 a12(e2,10)=-( clap2m2y*an3/epsmu2 + ( clap2m2x     - (clap2m2x*an3-clap2m2z*an1)*an3 )*eqn2Coeff + ( c2m2x     - (c2m2x*an3-c2m2z*an1)*an3 )*eqn2Coeffb ) 
 a12(e2,11)=-( clap2m2z*an3/epsmu2 + (              - (clap2m2y*an1-clap2m2x*an2)*an3 )*eqn2Coeff + (           - (c2m2y*an1-c2m2x*an2)*an3 )*eqn2Coeffb ) 

 ! Equations 9,10,11
 !  [ (1/mu)*n^T( c^2*LapSqE - alphaP*Lap(Ptt) )n  + (I-n n^T)*( c^4LapSq(E) - c^2*alphaP*Lap(Ptt) -alphaP*Ptttt ]
 !  Note: LE = c^2*Delta(E)

 ! first term:  (c^2*LapSqE - alphaP*Lap(Ptt) )/mu 
 coeffLap1a   =               -alphaP1*( c2PttEsum1         )/mu1
 coeffLapSq1a = 1./epsmu1/mu1 -alphaP1*( c2PttLEsum1/epsmu1 )/mu1  
 ! second term:  c^4LapSq(E) - alphaP*( c^2*Lap(Ptt) + Ptttt )
 coeffLap1b   =               -alphaP1*( c2PttEsum1  + c2PttttLEsum1  )/epsmu1
 coeffLapSq1b = 1./epsmu1**2  -alphaP1*( c2PttLEsum1 + c2PttttLLEsum1 )/epsmu1**2
 
 eqn1Coeff1a = clap1m1*coeffLap1a + clapSq1m1*coeffLapSq1a
 eqn1Coeff1b = clap1m1*coeffLap1b + cLapSq1m1*coeffLapSq1b 

 eqn1Coeff2a = clap1m2*coeffLap1a + clapSq1m2*coeffLapSq1a
 eqn1Coeff2b = clap1m2*coeffLap1b + cLapSq1m2*coeffLapSq1b 

 ! first term:  c^2*LapSqE - alphaP*Delta(Ptt)
 coeffLap2a   =               -alphaP2*( c2PttEsum2         )/mu2
 coeffLapSq2a = 1./epsmu2/mu2 -alphaP2*( c2PttLEsum2/epsmu2 )/mu2   
 ! second term:  c^4LapSq(E) - alphaP*( c^2*Lap(Ptt) + Ptttt )
 coeffLap2b   =               -alphaP2*( c2PttEsum2  + c2PttttLEsum2  )/epsmu2
 coeffLapSq2b = 1./epsmu2**2  -alphaP2*( c2PttLEsum2 + c2PttttLLEsum2 )/epsmu2**2

 eqn2Coeff1a = clap2m1*coeffLap2a + clapSq2m1*coeffLapSq2a 
 eqn2Coeff1b = clap2m1*coeffLap2b + cLapSq2m1*coeffLapSq2b 

 eqn2Coeff2a = clap2m2*coeffLap2a + clapSq2m2*coeffLapSq2a 
 eqn2Coeff2b = clap2m2*coeffLap2b + cLapSq2m2*coeffLapSq2b

 ! -----Equation 9-----
 ! coeffs of u1(-1), v1(-1), w1(-1)
 a12(e3,0) = ( an1*an1*eqn1Coeff1a + (1.-an1*an1)*eqn1Coeff1b )
 a12(e3,1) = ( an1*an2*eqn1Coeff1a + (  -an1*an2)*eqn1Coeff1b )
 a12(e3,2) = ( an1*an3*eqn1Coeff1a + (  -an1*an3)*eqn1Coeff1b )

 ! coeffs of u1(-2), v1(-2), w1(-2)
 a12(e3,3) = ( an1*an1*eqn1Coeff2a + (1.-an1*an1)*eqn1Coeff2b )
 a12(e3,4) = ( an1*an2*eqn1Coeff2a + (  -an1*an2)*eqn1Coeff2b )
 a12(e3,5) = ( an1*an3*eqn1Coeff2a + (  -an1*an3)*eqn1Coeff2b )

 ! coeffs of u2(-1), v2(-1), w2(-1) 
 a12(e3,6) =-( an1*an1*eqn2Coeff1a + (1.-an1*an1)*eqn2Coeff1b )
 a12(e3,7) =-( an1*an2*eqn2Coeff1a + (  -an1*an2)*eqn2Coeff1b )
 a12(e3,8) =-( an1*an3*eqn2Coeff1a + (  -an1*an3)*eqn2Coeff1b )

 ! coeffs of u2(-2), v2(-2), w2(-2)
 a12(e3,9) =-( an1*an1*eqn2Coeff2a + (1.-an1*an1)*eqn2Coeff2b )
 a12(e3,10)=-( an1*an2*eqn2Coeff2a + (  -an1*an2)*eqn2Coeff2b )
 a12(e3,11)=-( an1*an3*eqn2Coeff2a + (  -an1*an3)*eqn2Coeff2b )

 ! -----Equation 10 -----
 ! coeffs of u1(-1), v1(-1), w1(-1)
 a12(e4,0) = ( an2*an1*eqn1Coeff1a + (  -an2*an1)*eqn1Coeff1b )
 a12(e4,1) = ( an2*an2*eqn1Coeff1a + (1.-an2*an2)*eqn1Coeff1b )
 a12(e4,2) = ( an2*an3*eqn1Coeff1a + (  -an2*an3)*eqn1Coeff1b )

 ! coeffs of u1(-2), v1(-2), w1(-2)
 a12(e4,3) = ( an2*an1*eqn1Coeff2a + (  -an2*an1)*eqn1Coeff2b )
 a12(e4,4) = ( an2*an2*eqn1Coeff2a + (1.-an2*an2)*eqn1Coeff2b )
 a12(e4,5) = ( an2*an3*eqn1Coeff2a + (  -an2*an3)*eqn1Coeff2b )

 ! coeffs of u2(-1), v2(-1), w2(-1) 
 a12(e4,6) =-( an2*an1*eqn2Coeff1a + (  -an2*an1)*eqn2Coeff1b )
 a12(e4,7) =-( an2*an2*eqn2Coeff1a + (1.-an2*an2)*eqn2Coeff1b )
 a12(e4,8) =-( an2*an3*eqn2Coeff1a + (  -an2*an3)*eqn2Coeff1b )

 ! coeffs of u2(-2), v2(-2), w2(-2)
 a12(e4,9) =-( an2*an1*eqn2Coeff2a + (  -an2*an1)*eqn2Coeff2b )
 a12(e4,10)=-( an2*an2*eqn2Coeff2a + (1.-an2*an2)*eqn2Coeff2b )
 a12(e4,11)=-( an2*an3*eqn2Coeff2a + (  -an2*an3)*eqn2Coeff2b )

 ! -----Equation 11 -----
 ! coeffs of u1(-1), v1(-1), w1(-1)
 a12(e5,0) = ( an3*an1*eqn1Coeff1a + (  -an3*an1)*eqn1Coeff1b )
 a12(e5,1) = ( an3*an2*eqn1Coeff1a + (  -an3*an2)*eqn1Coeff1b )
 a12(e5,2) = ( an3*an3*eqn1Coeff1a + (1.-an3*an3)*eqn1Coeff1b )

 ! coeffs of u1(-2), v1(-2), w1(-2)
 a12(e5,3) = ( an3*an1*eqn1Coeff2a + (  -an3*an1)*eqn1Coeff2b )
 a12(e5,4) = ( an3*an2*eqn1Coeff2a + (  -an3*an2)*eqn1Coeff2b )
 a12(e5,5) = ( an3*an3*eqn1Coeff2a + (1.-an3*an3)*eqn1Coeff2b )

 ! coeffs of u2(-1), v2(-1), w2(-1) 
 a12(e5,6) =-( an3*an1*eqn2Coeff1a + (  -an3*an1)*eqn2Coeff1b )
 a12(e5,7) =-( an3*an2*eqn2Coeff1a + (  -an3*an2)*eqn2Coeff1b )
 a12(e5,8) =-( an3*an3*eqn2Coeff1a + (1.-an3*an3)*eqn2Coeff1b )

 ! coeffs of u2(-2), v2(-2), w2(-2)
 a12(e5,9) =-( an3*an1*eqn2Coeff2a + (  -an3*an1)*eqn2Coeff2b )
 a12(e5,10)=-( an3*an2*eqn2Coeff2a + (  -an3*an2)*eqn2Coeff2b )
 a12(e5,11)=-( an3*an3*eqn2Coeff2a + (1.-an3*an3)*eqn2Coeff2b )

#endMacro



! --------------------------------------------------------------------------
! Macro: Assign interface ghost values, DIM=3, ORDER=4, GRID=Curvilinear
!         DISPERSIVE CASE -- GDM 
! --------------------------------------------------------------------------
#beginMacro assignDispersiveInterfaceGhost3dOrder4()
  INFO("34c-GDM")

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
  fevz1(n)=0.

  fevx2(n)=0.
  fevy2(n)=0.
  fevz2(n)=0.

  do jv=0,numberOfPolarizationVectors1-1
    fpv1(n,jv)=0.
    LfP1(n,jv)=0.
    fPt1(n,jv)=0.
    fPtt1(n,jv)=0.

    fpvx1(n,jv)=0.
    fpvy1(n,jv)=0.
    fpvz1(n,jv)=0.
  end do
  do jv=0,numberOfPolarizationVectors2-1
    fpv2(n,jv)=0.
    LfP2(n,jv)=0.
    fPt2(n,jv)=0.
    fPtt2(n,jv)=0.

    fpvx2(n,jv)=0.
    fpvy2(n,jv)=0.
    fpvz2(n,jv)=0.
  end do
end do

  err2=0.
  count=0
  beginLoopsMask3d()

    ! here is the normal (assumed to be the same on both sides)
    an1=rsxy1(i1,i2,i3,axis1,0)   ! normal (an1,an2)
    an2=rsxy1(i1,i2,i3,axis1,1)
    an3=rsxy1(i1,i2,i3,axis1,2)
    aNorm=max(epsx,sqrt(an1**2+an2**2+an3**2))
    an1=an1/aNorm
    an2=an2/aNorm
    an3=an3/aNorm

    cem1=(1.-1./eps1)/mu1
    cem2=(1.-1./eps2)/mu2

    ! ---- Evaluate the jump conditions using the wrong values at the ghost points  ------
    evaluateDispersiveInterfaceEquations3dOrder4()

   !  write(debugFile,'(" f(6)=",1pe10.2," curlfPttx1,nDotCurlfPtt1=",2(1pe10.2))') f(6),curlfPttx1,nDotCurlfPtt1
   !  write(debugFile,'(" f(7)=",1pe10.2," curlfPttx2,nDotCurlfPtt2=",2(1pe10.2))') f(7),curlfPttx2,nDotCurlfPtt2


    if( debug.gt.4 )then
     write(debugFile,'(/," --> 3d-order4-curv: i1,i2,i3=",3i4," an1,an2,an3=",3e11.3)') i1,i2,i3,an1,an2,an3
     write(debugFile,'(" --> 3d-order4-curv: i1,i2,i3=",3i4," f(start)=",12e10.2)') i1,i2,i3,(f(n),n=0,11)
     ! '
    end if
    if( debug.gt.8 )then
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


    c1m1x = -is*8.*rsxy1(i1,i2,i3,axis1,0)/(12.*dr1(axis1))    ! coeff of u1(-1) from D.x
    c1m1y = -is*8.*rsxy1(i1,i2,i3,axis1,1)/(12.*dr1(axis1))    ! coeff of u1(-1) from D.y 
    c1m1z = -is*8.*rsxy1(i1,i2,i3,axis1,2)/(12.*dr1(axis1))    ! coeff of u1(-1) from D.z

    c1m2x =  is   *rsxy1(i1,i2,i3,axis1,0)/(12.*dr1(axis1))    ! coeff of u1(-2) from D.x
    c1m2y =  is   *rsxy1(i1,i2,i3,axis1,1)/(12.*dr1(axis1))    ! coeff of u1(-2) from D.y 
    c1m2z =  is   *rsxy1(i1,i2,i3,axis1,2)/(12.*dr1(axis1))    ! coeff of u1(-2) from D.z

    c2m1x = -js*8.*rsxy2(j1,j2,j3,axis2,0)/(12.*dr2(axis2))    ! coeff of u2(-1) from D.x
    c2m1y = -js*8.*rsxy2(j1,j2,j3,axis2,1)/(12.*dr2(axis2))
    c2m1z = -js*8.*rsxy2(j1,j2,j3,axis2,2)/(12.*dr2(axis2))

    c2m2x =  js   *rsxy2(j1,j2,j3,axis2,0)/(12.*dr2(axis2))
    c2m2y =  js   *rsxy2(j1,j2,j3,axis2,1)/(12.*dr2(axis2))
    c2m2z =  js   *rsxy2(j1,j2,j3,axis2,2)/(12.*dr2(axis2))

    rxx1(0,0,0)=aj1rxx
    rxx1(0,0,1)=aj1rxy
    rxx1(0,0,2)=aj1rxz
    rxx1(0,1,1)=aj1ryy
    rxx1(0,1,2)=aj1ryz
    rxx1(0,2,2)=aj1rzz

    rxx1(1,0,0)=aj1sxx
    rxx1(1,0,1)=aj1sxy
    rxx1(1,0,2)=aj1sxz
    rxx1(1,1,1)=aj1syy
    rxx1(1,1,2)=aj1syz
    rxx1(1,2,2)=aj1szz

    rxx1(2,0,0)=aj1txx
    rxx1(2,0,1)=aj1txy
    rxx1(2,0,2)=aj1txz
    rxx1(2,1,1)=aj1tyy
    rxx1(2,1,2)=aj1tyz
    rxx1(2,2,2)=aj1tzz


    rxx2(0,0,0)=aj2rxx
    rxx2(0,0,1)=aj2rxy
    rxx2(0,0,2)=aj2rxz
    rxx2(0,1,1)=aj2ryy
    rxx2(0,1,2)=aj2ryz
    rxx2(0,2,2)=aj2rzz

    rxx2(1,0,0)=aj2sxx
    rxx2(1,0,1)=aj2sxy
    rxx2(1,0,2)=aj2sxz
    rxx2(1,1,1)=aj2syy
    rxx2(1,1,2)=aj2syz
    rxx2(1,2,2)=aj2szz

    rxx2(2,0,0)=aj2txx
    rxx2(2,0,1)=aj2txy
    rxx2(2,0,2)=aj2txz
    rxx2(2,1,1)=aj2tyy
    rxx2(2,1,2)=aj2tyz
    rxx2(2,2,2)=aj2tzz


    ! clap1m1 : coeff of u(-1) from lap1 = u1.xx + u1.yy + u1.zz
    ! clap1m2 : coeff of u(-2) from lap1 = u1.xx + u1.yy + u1.zz

    
    clap1m1=4./3.*(rsxy1(i1,i2,i3,axis1,0)**2+rsxy1(i1,i2,i3,axis1,1)**2+rsxy1(i1,i2,i3,axis1,2)**2)/(dr1(axis1)**2) \
              -is*2./3.*(rxx1(axis1,0,0)+rxx1(axis1,1,1)+rxx1(axis1,2,2))/(2.*dr1(axis1))
    clap1m2=-1./12.*(rsxy1(i1,i2,i3,axis1,0)**2+rsxy1(i1,i2,i3,axis1,1)**2+rsxy1(i1,i2,i3,axis1,2)**2)/(dr1(axis1)**2) \
              +is*1./12.*(rxx1(axis1,0,0)+rxx1(axis1,1,1)+rxx1(axis1,2,2))/(2.*dr1(axis1)) 

    clap2m1=4/3.*(rsxy2(j1,j2,j3,axis2,0)**2+rsxy2(j1,j2,j3,axis2,1)**2+rsxy2(j1,j2,j3,axis2,2)**2)/(dr2(axis2)**2) \
              -js*2./3.*(rxx2(axis2,0,0)+rxx2(axis2,1,1)+rxx2(axis2,2,2))/(2.*dr2(axis2)) 
    clap2m2=-1./12.*(rsxy2(j1,j2,j3,axis2,0)**2+rsxy2(j1,j2,j3,axis2,1)**2+rsxy2(j1,j2,j3,axis2,2)**2)/(dr2(axis2)**2) \
              +js*1./12.*(rxx2(axis2,0,0)+rxx2(axis2,1,1)+rxx2(axis2,2,2))/(2.*dr2(axis2))


    ! rx1(m,n) = D r_m / D x_n
    ! rxx1(m,n1,n2) = D^2 r /( D x_n1 D X_n2 )
    ! rxxx1(m,n1,n2,n3) 

    ! first derivatives
    rx1(0,0)=rsxy1(i1,i2,i3,0,0)
    rx1(0,1)=rsxy1(i1,i2,i3,0,1)
    rx1(0,2)=rsxy1(i1,i2,i3,0,2)
    rx1(1,0)=rsxy1(i1,i2,i3,1,0)
    rx1(1,1)=rsxy1(i1,i2,i3,1,1)
    rx1(1,2)=rsxy1(i1,i2,i3,1,2)
    rx1(2,0)=rsxy1(i1,i2,i3,2,0)
    rx1(2,1)=rsxy1(i1,i2,i3,2,1)
    rx1(2,2)=rsxy1(i1,i2,i3,2,2)

    rx2(0,0)=rsxy2(j1,j2,j3,0,0)
    rx2(0,1)=rsxy2(j1,j2,j3,0,1)
    rx2(0,2)=rsxy2(j1,j2,j3,0,2)
    rx2(1,0)=rsxy2(j1,j2,j3,1,0)
    rx2(1,1)=rsxy2(j1,j2,j3,1,1)
    rx2(1,2)=rsxy2(j1,j2,j3,1,2)
    rx2(2,0)=rsxy2(j1,j2,j3,2,0)
    rx2(2,1)=rsxy2(j1,j2,j3,2,1)
    rx2(2,2)=rsxy2(j1,j2,j3,2,2)


    ! 3rd derivatives: (only some are needed)
    ! note for last 3 entries - entries must increase or stay the same
    rxxx1(axis1,0,0,0) = aj1rxxx
    rxxx1(axis1,1,1,1) = aj1ryyy
    rxxx1(axis1,2,2,2) = aj1rzzz
    rxxx1(axis1,0,0,1) = aj1rxxy
    rxxx1(axis1,0,0,2) = aj1rxxz
    rxxx1(axis1,0,1,1) = aj1rxyy
    rxxx1(axis1,0,2,2) = aj1rxzz
    rxxx1(axis1,1,1,2) = aj1ryyz
    rxxx1(axis1,1,2,2) = aj1ryzz

    rxxx2(axis2,0,0,0) = aj2rxxx
    rxxx2(axis2,1,1,1) = aj2ryyy
    rxxx2(axis2,2,2,2) = aj2rzzz
    rxxx2(axis2,0,0,1) = aj2rxxy
    rxxx2(axis2,0,0,2) = aj2rxxz
    rxxx2(axis2,0,1,1) = aj2rxyy
    rxxx2(axis2,0,2,2) = aj2rxzz
    rxxx2(axis2,1,1,2) = aj2ryyz
    rxxx2(axis2,1,2,2) = aj2ryzz



    ! Some 4th derivatives are needed by LapSq: 
    rxxxx1(axis1,0,0,0,0) = aj1rxxxx
    rxxxx1(axis1,1,1,1,1) = aj1ryyyy
    rxxxx1(axis1,2,2,2,2) = aj1rzzzz

    rxxxx1(axis1,0,0,1,1) = aj1rxxyy
    rxxxx1(axis1,0,0,2,2) = aj1rxxzz
    rxxxx1(axis1,1,1,2,2) = aj1ryyzz

    rxxxx2(axis2,0,0,0,0) = aj2rxxxx
    rxxxx2(axis2,1,1,1,1) = aj2ryyyy
    rxxxx2(axis2,2,2,2,2) = aj2rzzzz

    rxxxx2(axis2,0,0,1,1) = aj2rxxyy
    rxxxx2(axis2,0,0,2,2) = aj2rxxzz
    rxxxx2(axis2,1,1,2,2) = aj2ryyzz


 
    ! coeff of u1(-1) and u1(-2) from Lap^2
    ! dr1a, dr2a : used for avoidInterfaceIterations - tangential spacings are dsBig to eliminate mixed derivatives
    clapSq1m1=lapSqCoeff3DOrder2a(is,dr1a(axis1),dr1a(axis1p1),dr1a(axis1p2),axis1,axis1p1,axis1p2,rx1,rxx1,rxxx1,rxxxx1)
    clapSq1m2=lapSqCoeff3DOrder2b(is,dr1a(axis1),dr1a(axis1p1),dr1a(axis1p2),axis1,axis1p1,axis1p2,rx1,rxx1,rxxx1,rxxxx1)


    ! coeff of u2(-1) and u2(-2) from Lap^2
    clapSq2m1=lapSqCoeff3DOrder2a(js,dr2a(axis2),dr2a(axis2p1),dr2a(axis2p2),axis2,axis2p1,axis2p2,rx2,rxx2,rxxx2,rxxxx2)
    clapSq2m2=lapSqCoeff3DOrder2b(js,dr2a(axis2),dr2a(axis2p1),dr2a(axis2p2),axis2,axis2p1,axis2p2,rx2,rxx2,rxxx2,rxxxx2)


    ! coeff of u1(-1) from D.x(Lap), D.y(Lap) and D.z(Lap): (divideb by eps*mu)
    clap1m1x = xLapCoeff3DOrder2a(is,dr1a(axis1),dr1a(axis1p1),dr1a(axis1p2),axis1,axis1p1,axis1p2,rx1,rxx1,rxxx1) 
    clap1m1y = yLapCoeff3DOrder2a(is,dr1a(axis1),dr1a(axis1p1),dr1a(axis1p2),axis1,axis1p1,axis1p2,rx1,rxx1,rxxx1)  
    clap1m1z = zLapCoeff3DOrder2a(is,dr1a(axis1),dr1a(axis1p1),dr1a(axis1p2),axis1,axis1p1,axis1p2,rx1,rxx1,rxxx1)  

    ! coeff of u1(-2) from D.x(Lap), D.y(Lap) and D.z(Lap): 
    clap1m2x = xLapCoeff3DOrder2b(is,dr1a(axis1),dr1a(axis1p1),dr1a(axis1p2),axis1,axis1p1,axis1p2,rx1,rxx1,rxxx1)  
    clap1m2y = yLapCoeff3DOrder2b(is,dr1a(axis1),dr1a(axis1p1),dr1a(axis1p2),axis1,axis1p1,axis1p2,rx1,rxx1,rxxx1)  
    clap1m2z = zLapCoeff3DOrder2b(is,dr1a(axis1),dr1a(axis1p1),dr1a(axis1p2),axis1,axis1p1,axis1p2,rx1,rxx1,rxxx1)  


    ! coeff of u2(-1) from D.x(Lap), D.y(Lap) and D.z(Lap): 
    clap2m1x = xLapCoeff3DOrder2a(js,dr2a(axis2),dr2a(axis2p1),dr2a(axis2p2),axis2,axis2p1,axis2p2,rx2,rxx2,rxxx2)  
    clap2m1y = yLapCoeff3DOrder2a(js,dr2a(axis2),dr2a(axis2p1),dr2a(axis2p2),axis2,axis2p1,axis2p2,rx2,rxx2,rxxx2) 
    clap2m1z = zLapCoeff3DOrder2a(js,dr2a(axis2),dr2a(axis2p1),dr2a(axis2p2),axis2,axis2p1,axis2p2,rx2,rxx2,rxxx2)  

    ! coeff of u2(-2) from D.x(Lap), D.y(Lap) and D.z(Lap): 
    clap2m2x = xLapCoeff3DOrder2b(js,dr2a(axis2),dr2a(axis2p1),dr2a(axis2p2),axis2,axis2p1,axis2p2,rx2,rxx2,rxxx2)  
    clap2m2y = yLapCoeff3DOrder2b(js,dr2a(axis2),dr2a(axis2p1),dr2a(axis2p2),axis2,axis2p1,axis2p2,rx2,rxx2,rxxx2)  
    clap2m2z = zLapCoeff3DOrder2b(js,dr2a(axis2),dr2a(axis2p1),dr2a(axis2p2),axis2,axis2p1,axis2p2,rx2,rxx2,rxxx2)  


    ! cdivE1 =  u.c1x + v.c1y + w.c1z
    ! nDotCurlE1 = (w1y-v1z)*an1 + (u1z-w1x)*an2 + (v1x-u1y)*an3

    ! 12 Unknowns:
    !   u1(-1), v1(-1), w1(-1), u1(-2), v1(-2), w1(-2), 
    !   u2(-1), v2(-1), w2(-1), u2(-2), v2(-2), w2(-2)  
    ! 12 Equations: 
    !    a12(eqn,unknown) 


    ! fill equations 0,..,5
    fillDispersiveEquations3dOrder4a(0,1,2,3,4,5)



    ! coeff of u1(-1) and u1(-2) from Lap^2
    ! clapSq1m1=lapSqCoeff3DOrder2a(is,dr1(axis1),dr1(axis1p1),dr1(axis1p2),axis1,axis1p1,axis1p2,rx1,rxx1,rxxx1,rxxxx1)
    ! clapSq1m2=lapSqCoeff3DOrder2b(is,dr1(axis1),dr1(axis1p1),dr1(axis1p2),axis1,axis1p1,axis1p2,rx1,rxx1,rxxx1,rxxxx1)


    ! coeff of u2(-1) and u2(-2) from Lap^2
    ! clapSq2m1=lapSqCoeff3DOrder2a(js,dr2(axis2),dr2(axis2p1),dr2(axis2p2),axis2,axis2p1,axis2p2,rx2,rxx2,rxxx2,rxxxx2)
    ! clapSq2m2=lapSqCoeff3DOrder2b(js,dr2(axis2),dr2(axis2p1),dr2(axis2p2),axis2,axis2p1,axis2p2,rx2,rxx2,rxxx2,rxxxx2)

    ! write(debugFile,'(" --> 3d-order4-c: i1,i2,i3=",3i4," c1m1x,c1m1y,c1m1z,clap1m1,clap1m2=",12e10.2)') i1,i2,i3,c1m1x,c1m1y,c1m1z,clap1m1,clap1m2

    ! fill equations 6...11
    fillDispersiveEquations3dOrder4b(6,7,8,9,10,11)

    ! Equation 6..11 :  extrapolate 2nd ghost point 
    cex1=1.
    cex2=-5. ! ** fix me ** orderOfExtrapolation for 2nd ghost point 
    ! *e678*
    ! do ii=6,11
    !   do jj=0,11
    !     a12(ii,jj)=0.
    !   end do
    ! end do
    ! a12(6,0)  = cex2   ! u1(-1)
    ! a12(6,3)  = cex1   ! u1(-2)

    ! a12(7,1)  = cex2   ! v1(-1)
    ! a12(7,4)  = cex1   ! v1(-2)

    ! a12(8,2)  = cex2   ! w1(-1)
    ! a12(8,5)  = cex1   ! w1(-2)

    ! a12(9,6)  = cex2   ! u2(-1)
    ! a12(9,9)  = cex1   ! u2(-2)
    ! a12(10,7) = cex2   ! v2(-1)
    ! a12(10,10)= cex1   ! v2(-2)

    ! a12(11,8) = cex2   ! w2(-1)
    ! a12(11,11)= cex1   ! w2(-2)


    ! --- check matrix coefficients by delta function approach ----
    if( checkCoeff.eq.1 .and. it.le.1 )then
      numberOfEquations=12
      checkCoefficients(i1,i2,i3, j1,j2,j3,numberOfEquations,a12,evaluateDispersiveInterfaceEquations3dOrder4 )
    end if

    ! fill in the current values for the unknowns: 
    q(0) = u1(i1-is1,i2-is2,i3-is3,ex)
    q(1) = u1(i1-is1,i2-is2,i3-is3,ey)
    q(2) = u1(i1-is1,i2-is2,i3-is3,ez)
    q(3) = u1(i1-2*is1,i2-2*is2,i3-2*is3,ex)
    q(4) = u1(i1-2*is1,i2-2*is2,i3-2*is3,ey)
    q(5) = u1(i1-2*is1,i2-2*is2,i3-2*is3,ez)

    q(6) = u2(j1-js1,j2-js2,j3-js3,ex)
    q(7) = u2(j1-js1,j2-js2,j3-js3,ey)
    q(8) = u2(j1-js1,j2-js2,j3-js3,ez)
    q(9) = u2(j1-2*js1,j2-2*js2,j3-2*js3,ex)
    q(10)= u2(j1-2*js1,j2-2*js2,j3-2*js3,ey)
    q(11)= u2(j1-2*js1,j2-2*js2,j3-2*js3,ez)

    ! subtract off the contributions from the wrong values at the ghost points:
    numberOfEquations=12
    do n=0,numberOfEquations-1
      f(n) = (a12(n,0)*q(0)+a12(n,1)*q(1)+a12(n,2)*q(2)+a12(n,3)*q(3)+a12(n,4)*q(4)+a12(n,5)*q(5)+\
              a12(n,6)*q(6)+a12(n,7)*q(7)+a12(n,8)*q(8)+a12(n,9)*q(9)+a12(n,10)*q(10)+a12(n,11)*q(11) ) - f(n)
    end do
    if( debug.gt.3 )then
      write(debugFile,'(" --> 3d:order4-c: f(subtract)=",12e10.2)') (f(n),n=0,11)
    end if
    if( .false. )then
      do n=0,numberOfEquations-1
        write(debugFile,'("a(i,j)=",12(1pe10.2))') (a12(n,m),m=0,11)
      end do
    end if

    ! solve A Q = F
    ! factor the matrix
    call dgeco( a12(0,0), numberOfEquations, numberOfEquations, ipvt(0),rcond,work(0))


    if( debug.gt.3 )then
      write(debugFile,'(" --> 3d:order4-c: rcond=",e10.2)') rcond
      write(debugFile,'(" eqn1Coeff1a,eqn1Coeff1b,eqn1Coeff2a,eqn1Coeff2b=",4(1pe10.2))') eqn1Coeff1a,eqn1Coeff1b,eqn1Coeff2a,eqn1Coeff2b
      write(debugFile,'(" eqn2Coeff1a,eqn2Coeff1b,eqn2Coeff2a,eqn2Coeff2b=",4(1pe10.2))') eqn2Coeff1a,eqn2Coeff1b,eqn2Coeff2a,eqn2Coeff2b
    end if

    ! solve
    job=0
    call dgesl( a12(0,0), numberOfEquations, numberOfEquations, ipvt(0), f(0), job)
    if( debug.gt.3 )then
      write(debugFile,'(" --> 3d:order4-c: f(solve)=",12e10.2)') (f(n),n=0,11)
      write(debugFile,'(" --> 3d:order4-c:        q=",12e10.2)') (q(n),n=0,11)
      write(debugFile,'(" --> 3d:order4-c:      f-q=",12e10.2)') (f(n)-q(n),n=0,11)
    end if
    ! fill in the answer:
    if( useJacobiUpdate.eq.0 )then
      u1(i1-is1,i2-is2,i3-is3,ex)      =f(0 )*omega+(1.-omega)*q(0)
      u1(i1-is1,i2-is2,i3-is3,ey)      =f(1 )*omega+(1.-omega)*q(1)
      u1(i1-is1,i2-is2,i3-is3,ez)      =f(2 )*omega+(1.-omega)*q(2)
      u1(i1-2*is1,i2-2*is2,i3-2*is3,ex)=f(3 )*omega+(1.-omega)*q(3)
      u1(i1-2*is1,i2-2*is2,i3-2*is3,ey)=f(4 )*omega+(1.-omega)*q(4)
      u1(i1-2*is1,i2-2*is2,i3-2*is3,ez)=f(5 )*omega+(1.-omega)*q(5)
  
      u2(j1-js1,j2-js2,j3-js3,ex)      =f(6 )*omega+(1.-omega)*q(6)
      u2(j1-js1,j2-js2,j3-js3,ey)      =f(7 )*omega+(1.-omega)*q(7)
      u2(j1-js1,j2-js2,j3-js3,ez)      =f(8 )*omega+(1.-omega)*q(8)
      u2(j1-2*js1,j2-2*js2,j3-2*js3,ex)=f(9 )*omega+(1.-omega)*q(9)
      u2(j1-2*js1,j2-2*js2,j3-2*js3,ey)=f(10)*omega+(1.-omega)*q(10)
      u2(j1-2*js1,j2-2*js2,j3-2*js3,ez)=f(11)*omega+(1.-omega)*q(11)
    else
      ! Jacobi-update -- save results in the work space

      wk1(i1-is1,i2-is2,i3-is3,ex)      =f(0 )*omega+(1.-omega)*q(0)
      wk1(i1-is1,i2-is2,i3-is3,ey)      =f(1 )*omega+(1.-omega)*q(1)
      wk1(i1-is1,i2-is2,i3-is3,ez)      =f(2 )*omega+(1.-omega)*q(2)
      wk1(i1-2*is1,i2-2*is2,i3-2*is3,ex)=f(3 )*omega+(1.-omega)*q(3)
      wk1(i1-2*is1,i2-2*is2,i3-2*is3,ey)=f(4 )*omega+(1.-omega)*q(4)
      wk1(i1-2*is1,i2-2*is2,i3-2*is3,ez)=f(5 )*omega+(1.-omega)*q(5)
  
      wk2(j1-js1,j2-js2,j3-js3,ex)      =f(6 )*omega+(1.-omega)*q(6)
      wk2(j1-js1,j2-js2,j3-js3,ey)      =f(7 )*omega+(1.-omega)*q(7)
      wk2(j1-js1,j2-js2,j3-js3,ez)      =f(8 )*omega+(1.-omega)*q(8)
      wk2(j1-2*js1,j2-2*js2,j3-2*js3,ex)=f(9 )*omega+(1.-omega)*q(9)
      wk2(j1-2*js1,j2-2*js2,j3-2*js3,ey)=f(10)*omega+(1.-omega)*q(10)
      wk2(j1-2*js1,j2-2*js2,j3-2*js3,ez)=f(11)*omega+(1.-omega)*q(11)

    end if

    ! compute the maximum change in the solution for this iteration
    do n=0,11
      err=max(err,abs(q(n)-f(n)))
      err2 = err2 + (q(n)-f(n))**2
      count = count + 1
    end do

    if( .true. .or. debug.gt.3 )then ! re-evaluate

     do n=0,11
       errv(n)=abs(q(n)-f(n))
     end do
     ! evalDerivs3dOrder4()
     ! eval3dJumpOrder4New()

     evaluateDispersiveInterfaceEquations3dOrder4()

     res=0.
     do n=0,11
       res=max(res,abs(f(n)))
     end do
     if( .false. .and. res.gt.1.e-9 )then
       write(debugFile,'(" --> ERR: 3d-GDM-order4-c: it=",i3," i1,i2,i3=",3i4," f(re-eval)=",12e10.2)') it,i1,i2,i3,f(0),f(1),f(2),f(3),f(4),f(5),f(6),f(7),f(8),f(9),f(10),f(11)
       write(debugFile,'("     it=",i3," err in ghost=",12e10.2)') it,(errv(n),n=0,11)
     end if
    end if

  endLoopsMask3d()

  err2 = sqrt( err2/count )

  if( checkCoeff.eq.1 .and. it.le.1 )then
    write(*,'("+++++ IGDM34c: check coeff in interface: max(diff) = ",1pe8.2)') coeffDiff
  end if

  if( useJacobiUpdate.ne.0 )then
    ! Jacobi update -- copy work-space to solution arrays
    beginLoopsMask3d()
      u1(i1-is1,i2-is2,i3-is3,ex)      = wk1(i1-is1,i2-is2,i3-is3,ex)      
      u1(i1-is1,i2-is2,i3-is3,ey)      = wk1(i1-is1,i2-is2,i3-is3,ey)      
      u1(i1-is1,i2-is2,i3-is3,ez)      = wk1(i1-is1,i2-is2,i3-is3,ez)      
      u1(i1-2*is1,i2-2*is2,i3-2*is3,ex)= wk1(i1-2*is1,i2-2*is2,i3-2*is3,ex)
      u1(i1-2*is1,i2-2*is2,i3-2*is3,ey)= wk1(i1-2*is1,i2-2*is2,i3-2*is3,ey)
      u1(i1-2*is1,i2-2*is2,i3-2*is3,ez)= wk1(i1-2*is1,i2-2*is2,i3-2*is3,ez)
                                                                           
      u2(j1-js1,j2-js2,j3-js3,ex)      = wk2(j1-js1,j2-js2,j3-js3,ex)      
      u2(j1-js1,j2-js2,j3-js3,ey)      = wk2(j1-js1,j2-js2,j3-js3,ey)      
      u2(j1-js1,j2-js2,j3-js3,ez)      = wk2(j1-js1,j2-js2,j3-js3,ez)      
      u2(j1-2*js1,j2-2*js2,j3-2*js3,ex)= wk2(j1-2*js1,j2-2*js2,j3-2*js3,ex)
      u2(j1-2*js1,j2-2*js2,j3-2*js3,ey)= wk2(j1-2*js1,j2-2*js2,j3-2*js3,ey)
      u2(j1-2*js1,j2-2*js2,j3-2*js3,ez)= wk2(j1-2*js1,j2-2*js2,j3-2*js3,ez)
    endLoopsMask3d()
  end if

#endMacro



!  -- this next include file holds macros for assigning interface in 3D order 2 
! eval3dJumpOrder2()
! evalInterfaceEquations23c()
! assignInterfaceGhost23c()
! initializeInterfaceVariablesMacro(LABEL)
! setIndexBoundsExtraGhost()
! resetIndexBounds()  
#Include "interfaceMacros3d.h"



      subroutine mxInterface3dOrder4( nd, nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,\
                               gridIndexRange1, u1,u1n,u1m, wk1,mask1,rsxy1, xy1, p1,p1n,p1m, boundaryCondition1, \
                               md1a,md1b,md2a,md2b,md3a,md3b,\
                               gridIndexRange2, u2,u2n,u2m, wk2,mask2,rsxy2, xy2, p2,p2n,p2m, boundaryCondition2, \
                               ipar, rpar, \
                               aa2,aa4,aa8, ipvt2,ipvt4,ipvt8, \
                               ierr )
! ===================================================================================
!  Interface boundary conditions for Maxwells Equations in 3D and order 4.
!
!  gridType : 0=rectangular, 1=curvilinear
!
!  u1: solution on the "left" of the interface
!  u2: solution on the "right" of the interface
!
!  aa2,aa4,aa8 : real work space arrays that must be saved from call to call
!  ipvt2,ipvt4,ipvt8: integer work space arrays that must be saved from call to call
! ===================================================================================

      implicit none

      integer nd, \
              nd1a,nd1b,nd2a,nd2b,nd3a,nd3b, \
              md1a,md1b,md2a,md2b,md3a,md3b, \
              n1a,n1b,n2a,n2b,n3a,n3b,  \
              m1a,m1b,m2a,m2b,m3a,m3b,  \
              ierr

      ! ------- arrays for domain 1 -------------
      real u1(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
      real u1n(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
      real u1m(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
      real wk1(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)

      integer mask1(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
      real rsxy1(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1,0:nd-1)
      real xy1(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1)

      ! polarization vectors
      real p1(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
      real p1n(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
      real p1m(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)

      integer gridIndexRange1(0:1,0:2),boundaryCondition1(0:1,0:2)

      ! ------- arrays for domain 2 -------------
      real u2(md1a:md1b,md2a:md2b,md3a:md3b,0:*)
      real u2n(md1a:md1b,md2a:md2b,md3a:md3b,0:*)
      real u2m(md1a:md1b,md2a:md2b,md3a:md3b,0:*)
      real wk2(md1a:md1b,md2a:md2b,md3a:md3b,0:*)

      integer mask2(md1a:md1b,md2a:md2b,md3a:md3b)
      real rsxy2(md1a:md1b,md2a:md2b,md3a:md3b,0:nd-1,0:nd-1)
      real xy2(md1a:md1b,md2a:md2b,md3a:md3b,0:nd-1)

      ! polarization vectors
      real p2(md1a:md1b,md2a:md2b,md3a:md3b,0:*)
      real p2n(md1a:md1b,md2a:md2b,md3a:md3b,0:*)
      real p2m(md1a:md1b,md2a:md2b,md3a:md3b,0:*)

      integer gridIndexRange2(0:1,0:2),boundaryCondition2(0:1,0:2)

      integer ipar(0:*)
      real rpar(0:*)

      ! work space arrays that must be saved from call to call:
      real aa2(0:1,0:1,0:1,0:*),aa4(0:3,0:3,0:1,0:*),aa8(0:7,0:7,0:1,0:*)
      integer ipvt2(0:1,0:*), ipvt4(0:3,0:*), ipvt8(0:7,0:*)

!     --- local variables ----
      
      integer side1,axis1,grid1,side2,axis2,grid2,gridType,orderOfAccuracy,orderOfExtrapolation,useForcing,\
        ex,ey,ez,hx,hy,hz,useWhereMask,debug,solveForE,solveForH,axis1p1,axis1p2,axis2p1,axis2p2,nn,n1,n2,\
        twilightZone
      real dx1(0:2),dr1(0:2),dx2(0:2),dr2(0:2)
      real dr1a(0:2),dr2a(0:2)
      
      real t,ep,dt,eps1,mu1,c1,eps2,mu2,c2,epsmu1,epsmu2,omega
      integer axisp1,axisp2,i1,i2,i3,is1,is2,is3,j1,j2,j3,js1,js2,js3,ks1,ks2,ks3,is,js,it,nit,k1,k2,k3
      integer interfaceOption,interfaceEquationsOption,initialized,forcingOption,numberOfIterations

      integer assignInterfaceValues,assignInterfaceGhostValues,setDivergenceAtInterfaces
      real absoluteErrorTolerance,relativeErrorTolerance

      integer numGhost,giveDiv
      ! grid1
      integer nn1a,nn1b,nn2a,nn2b,nn3a,nn3b
      integer ns1a,ns1b,ns2a,ns2b,ns3a,ns3b  ! save n1a,n1b,...
      integer ne1a,ne1b,ne2a,ne2b,ne3a,ne3b  ! one extra ghost for 4th-order Stage I order2 

      ! grid2
      integer ms1a,ms1b,ms2a,ms2b,ms3a,ms3b  ! save m1a,m1b,...
      integer me1a,me1b,me2a,me2b,me3a,me3b  ! one extra ghost for 4th-order Stage I order2 
      integer mm1a,mm1b,mm2a,mm2b,mm3a,mm3b
      integer m1,m2

      logical extrapolateSecondGhost

      ! real rx1,ry1,rz1,rx2,ry2,rz2
      real x1,y1,z1

      real aLap0,aLap1,bLap0,bLap1,aLapX0,aLapX1,bLapY0,bLapY1,cLapX0,cLapX1,dLapY0,dLapY1,aLapSq0,aLapSq1,bLapSq0,bLapSq1
      real a0,a1,b0,b1,cc0,cc1,d0,d1,dr0,ds0
      real aNormSq,divu

      real epsRatio,an1,an2,an3,aNorm,ua,ub,uc,nDotU
      real epsx

      real tau1,tau2,clap1,clap2,u1Lap,v1Lap,w1Lap,u2Lap,v2Lap,w2Lap,an1Cartesian,an2Cartesian,an3Cartesian
      real u1LapSq,v1LapSq,u2LapSq,v2LapSq,w1LapSq,w2LapSq


      integer np1a,np1b,np2a,np2b,np3a,np3b,diff(0:2)

      real rx,ry,rxx,rxy,ryy,rxxx,rxxy,rxyy,ryyy,rxxxx,rxxyy,ryyyy
      real sx,sy,sxx,sxy,syy,sxxx,sxxy,sxyy,syyy,sxxxx,sxxyy,syyyy

      integer numberOfEquations,job
      real a2(0:1,0:1),a4(0:3,0:3),a6(0:5,0:5),a8(0:7,0:7),a12(0:11,0:11),q(0:11),f(0:11),rcond,work(0:11),errv(0:11)
      integer ipvt(0:11)

      real err,res,errOld,errRatio,ratioAve
      real err2,err2Old,err2Ratio
      real count
      integer debugFile,myid,parallel
      character*20 debugFileName

      ! for new evaluation method:
      real u1x,u1y,u1z,u1xx,u1xy,u1yy,u1xz,u1yz,u1zz
      real u2x,u2y,u2z,u2xx,u2xy,u2yy,u2xz,u2yz,u2zz

      real v1x,v1y,v1z,v1xx,v1xy,v1yy,v1xz,v1yz,v1zz
      real v2x,v2y,v2z,v2xx,v2xy,v2yy,v2xz,v2yz,v2zz

      real w1x,w1y,w1z,w1xx,w1xy,w1yy,w1xz,w1yz,w1zz
      real w2x,w2y,w2z,w2xx,w2xy,w2yy,w2xz,w2yz,w2zz

      real u1nx,u1ny,u1nz,u1nxx,u1nxy,u1nyy,u1nxz,u1nyz,u1nzz,u1nLap
      real v1nx,v1ny,v1nz,v1nxx,v1nxy,v1nyy,v1nxz,v1nyz,v1nzz,v1nLap
      real w1nx,w1ny,w1nz,w1nxx,w1nxy,w1nyy,w1nxz,w1nyz,w1nzz,w1nLap

      real u2nx,u2ny,u2nz,u2nxx,u2nxy,u2nyy,u2nxz,u2nyz,u2nzz,u2nLap
      real v2nx,v2ny,v2nz,v2nxx,v2nxy,v2nyy,v2nxz,v2nyz,v2nzz,v2nLap
      real w2nx,w2ny,w2nz,w2nxx,w2nxy,w2nyy,w2nxz,w2nyz,w2nzz,w2nLap

      real p1x,p1y,p1z,p1xx,p1xy,p1yy,p1xz,p1yz,p1zz,p1Lap
      real p1nx,p1ny,p1nz,p1nxx,p1nxy,p1nyy,p1nxz,p1nyz,p1nzz,p1nLap
      real p2x,p2y,p2z,p2xx,p2xy,p2yy,p2xz,p2yz,p2zz,p2Lap
      real p2nx,p2ny,p2nz,p2nxx,p2nxy,p2nyy,p2nxz,p2nyz,p2nzz,p2nLap

      real u1xxx,u1xxy,u1xyy,u1yyy, u1xxz,u1xzz,u1zzz, u1yyz, u1yzz
      real u2xxx,u2xxy,u2xyy,u2yyy, u2xxz,u2xzz,u2zzz, u2yyz, u2yzz
      real v1xxx,v1xxy,v1xyy,v1yyy, v1xxz,v1xzz,v1zzz, v1yyz, v1yzz
      real v2xxx,v2xxy,v2xyy,v2yyy, v2xxz,v2xzz,v2zzz, v2yyz, v2yzz
      real w1xxx,w1xxy,w1xyy,w1yyy, w1xxz,w1xzz,w1zzz, w1yyz, w1yzz
      real w2xxx,w2xxy,w2xyy,w2yyy, w2xxz,w2xzz,w2zzz, w2yyz, w2yzz

      real u1xxxx,u1xxyy,u1yyyy, u1xxzz,u1zzzz, u1yyzz
      real u2xxxx,u2xxyy,u2yyyy, u2xxzz,u2zzzz, u2yyzz
      real v1xxxx,v1xxyy,v1yyyy, v1xxzz,v1zzzz, v1yyzz
      real v2xxxx,v2xxyy,v2yyyy, v2xxzz,v2zzzz, v2yyzz
      real w1xxxx,w1xxyy,w1yyyy, w1xxzz,w1zzzz, w1yyzz
      real w2xxxx,w2xxyy,w2yyyy, w2xxzz,w2zzzz, w2yyzz

      real rx1(0:2,0:2), rx2(0:2,0:2)
      real rxx1(0:2,0:2,0:2), rxx2(0:2,0:2,0:2)
      real rxxx1(0:2,0:2,0:2,0:2), rxxx2(0:2,0:2,0:2,0:2)
      real rxxxx1(0:2,0:2,0:2,0:2,0:2), rxxxx2(0:2,0:2,0:2,0:2,0:2)

      real dx112(0:2),dx122(0:2),dx212(0:2),dx222(0:2),dx141(0:2),dx142(0:2),dx241(0:2),dx242(0:2)
      real dr114(0:2),dr214(0:2)

      real cem1,divE1,curlE1x,curlE1y,curlE1z,nDotCurlE1,nDotLapE1
      real cem2,divE2,curlE2x,curlE2y,curlE2z,nDotCurlE2,nDotLapE2

      real divLapE1,curlLapE1x,curlLapE1y,curlLapE1z,nDotCurlLapE1,nDotLapSqE1
      real divLapE2,curlLapE2x,curlLapE2y,curlLapE2z,nDotCurlLapE2,nDotLapSqE2

      real curlfPttx1,curlfPtty1,curlfPttz1,nDotCurlfPtt1
      real curlfPttx2,curlfPtty2,curlfPttz2,nDotCurlfPtt2


      integer ii,jj
      real cex1,cex2
      real clap1m1,clap1m2, clap2m1,clap2m2
      real c1m1x,c1m1y,c1m1z, c1m2x,c1m2y,c1m2z
      real c2m1x,c2m1y,c2m1z, c2m2x,c2m2y,c2m2z
      real clap1m1x,clap1m1y,clap1m1z, clap1m2x,clap1m2y,clap1m2z
      real clap2m1x,clap2m1y,clap2m1z, clap2m2x,clap2m2y,clap2m2z
      real clapSq1m1,clapSq1m2,clapSq2m1,clapSq2m2

      real c1x,c1y,c1z
      real c2x,c2y,c2z

      ! these are for the exact solution from TZ flow: 
      real ue,ve,we
      real uex,uey,uez, vex,vey,vez, wex,wey,wez, hex,hey,hez
      real uexx,ueyy,uezz, vexx,veyy,vezz, wexx,weyy,wezz
      real ueLap, veLap, weLap
      real curlEex,curlEey,curlEez,nDotCurlEe,nDotLapEe,nDotCurlLapEe,nDotLapSqEe
      real uexxx,uexxy,uexyy,ueyyy, uexxz,uexzz,uexyz, ueyyz, ueyzz, uezzz
      real vexxx,vexxy,vexyy,veyyy, vexxz,vexzz,vexyz, veyyz, veyzz, vezzz
      real wexxx,wexxy,wexyy,weyyy, wexxz,wexzz,wexyz, weyyz, weyzz, wezzz

      real uexxxx,uexxyy,ueyyyy,uexxzz,ueyyzz,uezzzz,ueLapSq
      real vexxxx,vexxyy,veyyyy,vexxzz,veyyzz,vezzzz,veLapSq
      real wexxxx,wexxyy,weyyyy,wexxzz,weyyzz,wezzzz,weLapSq

      real weLaSqp,curlLapEex,curlLapEey,curlLapEez,nDotLapSq


      real nDotUm,nDotUp,nDotU0

      integer kd1,kd2,kd3,kd4,kd5
      real bigValue

      ! boundary conditions parameters
      #Include "bcDefineFortranInclude.h"
 
      integer rectangular,curvilinear
      parameter(\
        rectangular=0,\
        curvilinear=1)

      integer useImpedanceInterfaceProjection
      real  ex1,ey1,ez1, hz1, ex2,ey2,ez2, hz2, nDotE1, nDotE2, epsNDotEI,  nDotEI, nDotE1I, nDotE2I
      real  exI, eyI, ezI, hzI, g1,g2, nDotEe, eta1,eta2,eta1i,eta2i

      ! Dispersion models
      integer noDispersion,drude,gdm
      parameter( noDispersion=0, drude=1, gdm=2 )

      ! forcing options
      #Include "forcingDefineFortranInclude.h"

      ! Known solution options
      #Include "knownSolutionFortranInclude.h"
      integer knownSolutionOption

      integer useJacobiUpdate
      integer dispersionModel1, dispersionModel2, dispersive, jv, pxc
      integer gdmParOption
      integer maxNumberOfParameters,maxNumberOfPolarizationVectors,maxPolarizationComponents
      parameter( maxNumberOfParameters=4, maxNumberOfPolarizationVectors=20, maxPolarizationComponents=20*3 )
      integer numberOfPolarizationVectors1
      real gdmPar1(0:maxNumberOfParameters-1,0:maxNumberOfPolarizationVectors-1)
      real a0v1,a1v1,b0v1,b1v1

      integer numberOfPolarizationVectors2
      real gdmPar2(0:maxNumberOfParameters-1,0:maxNumberOfPolarizationVectors-1)
      real a0v2,a1v2,b0v2,b1v2

      integer numberOfTimeDerivatives
      real evals(0:2),evalse(0:2),dpdm
      real pvals(0:maxPolarizationComponents-1),pvalse(0:maxPolarizationComponents-1)
      real pvalsm(0:maxPolarizationComponents-1),pvalsp(0:maxPolarizationComponents-1)
     
      real pvc(0:maxNumberOfPolarizationVectors-1),pvm(0:maxNumberOfPolarizationVectors-1),fpv(0:maxNumberOfPolarizationVectors-1)
      real rhspv(0:maxNumberOfPolarizationVectors-1)
      real betav(0:maxNumberOfPolarizationVectors-1)
      real beta,fe,evm,pSum,rhsE,rhsP,tm
      real uet,uett
      integer addForcing

      ! variables for dispersive models
      real alphaP1,alphaP2
      real fp1(0:2), fp2(0:2)
      integer ec,pc,pce
      real pv,pvt,pvtt
      real ev,evt,evtt
      real Bk, Ck
      real Csum, beta1,beta2, c2Sum
      real fev1(0:2), fpSum1(0:2), fpv1(0:2,0:maxNumberOfPolarizationVectors-1)
      real fev2(0:2), fpSum2(0:2), fpv2(0:2,0:maxNumberOfPolarizationVectors-1)
      real es(0:2), est(0:2), estt(0:2), esxx(0:2), esyy(0:2), eszz(0:2)
      real pe(0:2), pet(0:2), pett(0:2), pettSum1(0:2), pettSum2(0:2)
      real p0,p0t,p0tt

      real p1v(0:2), p2v(0:2), D1v(0:2), D2v(0:2)
      real p1ev(0:2), p2ev(0:2)
      real nDotD1,nDotD2, nDotP1, nDotP2, nDotP1e,nDotP2e, nDotDe, nDotDI
      real nDotFp1, nDotFp2
      real betac1,betac2

      real nDotfLPtt1,nDotfPtttt1
      real nDotfLPtt2,nDotfPtttt2

      ! variables for 4th-order GDM interface
      real c2PttLE,c2PttE,c2PttEm,c2PttP,c2PttPm,c2PttfE,c2PttfP
      real c4PttLE,c4PttE,c4PttEm,c4PttP,c4PttPm,c4PttfE,c4PttfP
      real c4PttLLE,c4PttLP,c4PttLEm,c4PttLPm,c4PttLfE,c4PttLfP
      real c4PttfEt,c4PttfEtt,c4PttfPt,c4PttfPtt
      real c2EtLE,c2EtE,c2EtEm,c2EtP,c2EtPm,c2EtfE,c2EtfP
      real c2PtLE,c2PtE,c2PtEm,c2PtP,c2PtPm,c2PtfE,c2PtfP

      real c2PttttLE,c2PttttE,c2PttttEm,c2PttttP,c2PttttPm,c2PttttLLE,c2PttttLEm,c2PttttLP,c2PttttLPm
      real c2PttttfE,c2PttttfP,c2PttttfEt,c2PttttfPt,c2PttttfEtt,c2PttttfPtt,c2PttttLfe,c2PttttLfP

      real c2PttEsum1, c2PttEsum2, c2PttLEsum1, c2PttLEsum2
      real c4PttLEsum1, c4PttLLEsum1, c4PttLEsum2, c4PttLLEsum2
      real c2PttttLEsum1, c2PttttLLEsum1, c2PttttLEsum2, c2PttttLLEsum2

      real nDotPevttLSum1,nDotPevttttSum1
      real nDotPevttLSum2,nDotPevttttSum2


      integer ismooth,nsmooth
      real eqnCoeff

      real evn,pvn
      real alpha,d4,LP,LPm
      real LE(0:2)
      real LE1(0:2),LE1m(0:2),LLE1(0:2), LEx1(0:2),LEy1(0:2),LEz1(0:2)
      real LE2(0:2),LE2m(0:2),LLE2(0:2), LEx2(0:2),LEy2(0:2),LEz2(0:2) 

      real LfE1(0:2),fEt1(0:2),fEtt1(0:2)
      real LfP1(0:2,0:maxNumberOfPolarizationVectors-1),fPt1(0:2,0:maxNumberOfPolarizationVectors-1)
      real fPtt1(0:2,0:maxNumberOfPolarizationVectors-1)

      real LfE2(0:2),fEt2(0:2),fEtt2(0:2)
      real LfP2(0:2,0:maxNumberOfPolarizationVectors-1),fPt2(0:2,0:maxNumberOfPolarizationVectors-1)
      real fPtt2(0:2,0:maxNumberOfPolarizationVectors-1)
      real pevtt1(0:2,0:maxNumberOfPolarizationVectors-1)
      real pevtt2(0:2,0:maxNumberOfPolarizationVectors-1)
      real pevttx1(0:2,0:maxNumberOfPolarizationVectors-1),pevtty1(0:2,0:maxNumberOfPolarizationVectors-1),pevttz1(0:2,0:maxNumberOfPolarizationVectors-1)
      real pevttx2(0:2,0:maxNumberOfPolarizationVectors-1),pevtty2(0:2,0:maxNumberOfPolarizationVectors-1),pevttz2(0:2,0:maxNumberOfPolarizationVectors-1)
      real pevttSum1(0:2),pevttxSum1(0:2),pevttySum1(0:2),pevttzSum1(0:2)
      real pevttSum2(0:2),pevttxSum2(0:2),pevttySum2(0:2),pevttzSum2(0:2)
      real pevtttt1(0:2,0:maxNumberOfPolarizationVectors-1),pevtttt2(0:2,0:maxNumberOfPolarizationVectors-1)
      real pevttLSum1(0:2),pevttLSum2(0:2),pevttttSum1(0:2),pevttttSum2(0:2)

      real petttSum
      real nDotPevttSum1, nDotPevttSum2
      real curlPevttxSum1,curlPevttySum1,curlPevttzSum1,nDotCurlPevttSum1
      real curlPevttxSum2,curlPevttySum2,curlPevttzSum2,nDotCurlPevttSum2


      real esttt(0:2),estxx(0:2),estyy(0:2),estzz(0:2),esttxx(0:2),esttyy(0:2),esttzz(0:2)
      real esxxxx(0:2),esxxyy(0:2),esyyyy(0:2),esxxzz(0:2),esyyzz(0:2),eszzzz(0:2)
      real esttx(0:2),estty(0:2),esttz(0:2),esxxx(0:2),esxyy(0:2),esxzz(0:2),esxxy(0:2),esxxz(0:2),esyyy(0:2),esyyz(0:2),esyzz(0:2),eszzz(0:2)
      real esx(0:2),esy(0:2),esz(0:2), estx(0:2),esty(0:2),estz(0:2)
      real pettt(0:2),pexx(0:2),peyy(0:2),pezz(0:2),petxx(0:2),petyy(0:2),petzz(0:2),pettxx(0:2),pettyy(0:2),pettzz(0:2)
      real pettx(0:2),petty(0:2),pettz(0:2),pex(0:2),pey(0:2),pez(0:2),petx(0:2),pety(0:2),petz(0:2)

      real estttt(0:2),petttt(0:2)
      real esL,estL,esttL,esLL,esLx,esLy,esLz
      real peL,petL,pettL,ptta,ptttta

      real pvx,pvy,pvz, pvnx,pvny,pvnz, pttxa,pttya,pttza, Lptta
      real fPttx1(0:2), fPtty1(0:2), fPttz1(0:2), fPttx2(0:2), fPtty2(0:2), fPttz2(0:2), fLPtt1(0:2), fLPtt2(0:2), fPtttt1(0:2), fPtttt2(0:2)
      real evx1(0:2),evy1(0:2),evz1(0:2), evnx1(0:2),evny1(0:2),evnz1(0:2)
      real evx2(0:2),evy2(0:2),evz2(0:2), evnx2(0:2),evny2(0:2),evnz2(0:2)
      real fevx1(0:2),fevy1(0:2),fevz1(0:2)
      real fevx2(0:2),fevy2(0:2),fevz2(0:2)
      real fpvx1(0:2,0:maxNumberOfPolarizationVectors-1),fpvy1(0:2,0:maxNumberOfPolarizationVectors-1),fpvz1(0:2,0:maxNumberOfPolarizationVectors-1)
      real fpvx2(0:2,0:maxNumberOfPolarizationVectors-1),fpvy2(0:2,0:maxNumberOfPolarizationVectors-1),fpvz2(0:2,0:maxNumberOfPolarizationVectors-1)

      real eqn1Coeffm1,eqn1Coeffm2,eqn2Coeffm1,eqn2Coeffm2
      real eqn1Coeff,eqn2Coeff
      real eqn1Coeff1a,eqn1Coeff1b,eqn1Coeff2a,eqn1Coeff2b
      real eqn2Coeff1a,eqn2Coeff1b,eqn2Coeff2a,eqn2Coeff2b

      real eqn1Coeffb,eqn2Coeffb

      real coeffLap1a,coeffLapSq1a,coeffLap1b,coeffLapSq1b
      real coeffLap2a,coeffLapSq2a,coeffLap2b,coeffLapSq2b
      real nDotEP1,nDotEP2

      integer avoidInterfaceIterations
      real dsBig 

      ! For checking coefficients by delta approach: 
      real coeffDiff
      integer checkCoeff
      integer hw1,hw2,hw3
      real f0(0:11),delta

      integer internalGhostBC,numGhost2,numParallelGhost

      !...................start statement functions

      ! .......statement functions for GDM parameters
      a0v1(jv) = gdmPar1(0,jv)
      a1v1(jv) = gdmPar1(1,jv)
      b0v1(jv) = gdmPar1(2,jv)
      b1v1(jv) = gdmPar1(3,jv)

      a0v2(jv) = gdmPar2(0,jv)
      a1v2(jv) = gdmPar2(1,jv)
      b0v2(jv) = gdmPar2(2,jv)
      b1v2(jv) = gdmPar2(3,jv)
      integer kd,m,n

      declareTemporaryVariables(2,2)
      declareParametricDerivativeVariables(uu1,3)   ! declare temp variables uu, uur, uus, ...
      declareParametricDerivativeVariables(uu2,3) 
      declareParametricDerivativeVariables(vv1,3)   ! declare temp variables uu, uur, uus, ...
      declareParametricDerivativeVariables(vv2,3) 
      declareParametricDerivativeVariables(ww1,3)   ! declare temp variables uu, uur, uus, ...
      declareParametricDerivativeVariables(ww2,3) 
      declareJacobianDerivativeVariables(aj1,3)     ! declareJacobianDerivativeVariables(aj,DIM)
      declareJacobianDerivativeVariables(aj2,3)     ! declareJacobianDerivativeVariables(aj,DIM)

      !................... end statement functions

      ierr=0

      side1                =ipar(0)
      axis1                =ipar(1)
      grid1                =ipar(2)
      n1a                  =ipar(3)
      n1b                  =ipar(4)
      n2a                  =ipar(5)
      n2b                  =ipar(6)
      n3a                  =ipar(7)
      n3b                  =ipar(8)

      side2                =ipar(9)
      axis2                =ipar(10)
      grid2                =ipar(11)
      m1a                  =ipar(12)
      m1b                  =ipar(13)
      m2a                  =ipar(14)
      m2b                  =ipar(15)
      m3a                  =ipar(16)
      m3b                  =ipar(17)

      gridType             =ipar(18)
      orderOfAccuracy      =ipar(19)
      orderOfExtrapolation =ipar(20)  ! maximum allowable order of extrapolation
      useForcing           =ipar(21)
      ex                   =ipar(22)
      ey                   =ipar(23)
      ez                   =ipar(24)
      hx                   =ipar(25)
      hy                   =ipar(26)
      hz                   =ipar(27)
      solveForE            =ipar(28)
      solveForH            =ipar(29)
      useWhereMask         =ipar(30)
      debug                =ipar(31)
      nit                  =ipar(32)
      interfaceOption      =ipar(33)
      initialized          =ipar(34)
      myid                 =ipar(35)
      parallel             =ipar(36)
      forcingOption        =ipar(37)
      interfaceEquationsOption=ipar(38)
      assignInterfaceValues     =ipar(39)
      assignInterfaceGhostValues=ipar(40)
      setDivergenceAtInterfaces =ipar(41)

      useImpedanceInterfaceProjection=ipar(42)
      ! numberOfInterfaceIterationsUsed = ipar(43)  ! returned value 
      ipar(43)=0
      
      dispersionModel1    = ipar(44)
      dispersionModel2    = ipar(45)
      pxc                 = ipar(46)
      knownSolutionOption = ipar(47)
      useJacobiUpdate     = ipar(48)      

      ! nonlinearModel1     = ipar(49)
      ! nonlinearModel2     = ipar(50)

      numParallelGhost    = ipar(51)
      internalGhostBC     = ipar(52)  ! bc value for internal parallel boundaries 

      ! numberOfInterfaceIterationsUsed = ipar(43)  ! returned value 
     
      dx1(0)                =rpar(0)
      dx1(1)                =rpar(1)
      dx1(2)                =rpar(2)
      dr1(0)                =rpar(3)
      dr1(1)                =rpar(4)
      dr1(2)                =rpar(5)

      dx2(0)                =rpar(6)
      dx2(1)                =rpar(7)
      dx2(2)                =rpar(8)
      dr2(0)                =rpar(9)
      dr2(1)                =rpar(10)
      dr2(2)                =rpar(11)

      t                    =rpar(12)
      ep                   =rpar(13) ! pointer for exact solution
      dt                   =rpar(14)
      eps1                 =rpar(15)
      mu1                  =rpar(16)
      c1                   =rpar(17)
      eps2                 =rpar(18)
      mu2                  =rpar(19)
      c2                   =rpar(20)
      omega                =rpar(21)
      ! rpar(22) : averageInterfaceConvergenceRate : return value 
      ! rpar(23) : maxFinalResidual : return value 
     

      eta1=sqrt(mu1/eps1) ! electrical impedance
      eta2=sqrt(mu2/eps2) ! electrical impedance
      eta1i=1./eta1
      eta2i=1./eta2

      epsmu1=eps1*mu1
      epsmu2=eps2*mu2

      twilightZone=useForcing

      avoidInterfaceIterations=1  ! option to avoid interface iterations. *new* July 6 2019 
      dsBig = 1.e20               ! Large value for tangential grid spacing to zero out mixed derivatives
      if( avoidInterfaceIterations.eq.1 )then     
        nit=1
        omega=1.    ! do not under-relax iterations

        axis1p1=mod(axis1+1,3)
        axis1p2=mod(axis1+2,3)
        axis2p1=mod(axis2+1,3)
        axis2p2=mod(axis2+2,3)

        ! dr1a(0:2) = sets grid spacing in tangential directions to dsBig 
        dr1a(axis1  )=dr1(axis1)
        dr1a(axis1p1)=dsBig
        dr1a(axis1p2)=dsBig

        dr2a(axis2  )=dr2(axis2)
        dr2a(axis2p1)=dsBig
        dr2a(axis2p2)=dsBig

      else
        do kd=0,2
          dr1a(kd)=dr1(kd)
          dr2a(kd)=dr2(kd)
        end do 
      end if

      checkCoeff=0 ! 1 ! if non-zero then check coefficients in interface equations using delta approach
      if( t.gt.1.5*dt )then
         checkCoeff=0   ! turn off check coeff for t>0
      end if
      coeffDiff=0.

      ! *wdh* Sept 5, 2016 absoluteErrorTolerance=1.e-10  ! fix me -- need a relative tol
      relativeErrorTolerance = rpar(24)
      absoluteErrorTolerance = rpar(25)

      debugFile=10

      if(  t.le. 1.5*dt .and. debug>0 )then
        write(*,'(" +++++++++cgmx interface3dOrder4 t=",e9.2," dt=",e9.2," nit=",i3," ++++++++")') t,dt,nit
           ! '
        write(*,'("  ... nd=",i2," gridType=",i2," order=",i2," debug=",i3,", ex=",i2)') nd,gridType,orderOfAccuracy,debug,ex
        write(*,'("  ... assignInterface=",i2," assignGhost=",i2)') assignInterfaceValues,assignInterfaceGhostValues
        write(*,'("  ... useJacobiUpdate=",i2," numParallelGhost=",i2)') useJacobiUpdate,numParallelGhost
        write(*,'("  ... useImpedanceInterfaceProjection=",i2," useJacobiUpdate=",i2)') useImpedanceInterfaceProjection,useJacobiUpdate
        write(*,'("  ... avoidInterfaceIterations=",i2)') avoidInterfaceIterations
        write(*,'("  ... interfaceOption=",i2)') interfaceOption
        write(*,'("  ... interface its (4th-order) relativeTol=",e12.3," absoluteTol=",e12.3)') relativeErrorTolerance,absoluteErrorTolerance
        write(*,'("  ... eps1,mu1=",2(1pe10.2)," eps2,mu2=",2(1pe10.2))') eps1,mu1,eps2,mu2
        write(*,'("  ... bc1=",6i6)') ((boundaryCondition1(i1,i2),i1=0,1),i2=0,nd-1)
        write(*,'("  ... bc2=",6i6)') ((boundaryCondition2(i1,i2),i1=0,1),i2=0,nd-1)
        write(*,'("  ... internalGhostBC=",i6)') internalGhostBC

      end if


!*      if( initialized.eq.0 .and. debug.gt.0 )then
!*        ! open debug files
!*        ! open (debugFile,file=filen,status='unknown',form='formatted')
!*        if( myid.lt.10 )then
!*          write(debugFileName,'("mxi",i1,".fdebug")') myid
!*        else
!*          write(debugFileName,'("mxi",i2,".fdebug")') myid
!*        end if
!*        write(*,*) 'interface3d: myid=',myid,' open debug file:',debugFileName
!*        open (debugFile,file=debugFileName,status='unknown',form='formatted')
!*        ! '
!*        ! INQUIRE(FILE=filen, EXIST=filex)
!*      end if
!*
!*      if( t.lt.dt )then
!*        write(debugFile,'(" +++++++++cgmx interface3d t=",e9.2," ++++++++")') t
!*           ! '
!*        write(debugFile,'(" interface3d new: nd=",i2," gridType=",i2)') nd,gridType
!*      end if
!*
!*      if( abs(c1*c1-1./(mu1*eps1)).gt. 1.e-10 )then
!*        write(debugFile,'(" interface3d:ERROR: c1,eps1,mu1=",3e10.2," not consistent")') c1,eps1,mu1
!*           ! '
!*        stop 11
!*      end if
!*      if( abs(c2*c2-1./(mu2*eps2)).gt. 1.e-10 )then
!*        write(debugFile,'(" interface3d:ERROR: c2,eps2,mu2=",3e10.2," not consistent")') c2,eps2,mu2
!*           ! '
!*        stop 11
!*      end if
!*
!*      if( .false. )then
!*        write(debugFile,'(" interface3d: eps1,eps2=",2f10.5," c1,c2=",2f10.5)') eps1,eps2,c1,c2
!*           ! '
!*      end if
!*
!*      if( nit.lt.0 .or. nit.gt.100 )then
!*        write(debugFile,'(" interfaceBC: ERROR: nit=",i9)') nit
!*        nit=max(1,min(100,nit))
!*      end if
!*
!*      if( debug.gt.1 )then
!*        write(debugFile,'("********************************************************************** ")')
!*        write(debugFile,'(" interface3d: **START** t=",e10.2)') t
!*        write(debugFile,'(" interface3d: **START** grid1=",i4," side1,axis1=",2i2," bc=",6i3)') grid1,side1,axis1,\
!*           boundaryCondition1(0,0),boundaryCondition1(1,0),boundaryCondition1(0,1),boundaryCondition1(1,1),boundaryCondition1(0,2),boundaryCondition1(1,2)
!*           ! '
!*        write(debugFile,'(" interface3d: **START** grid2=",i4," side2,axis2=",2i2," bc=",6i3)') grid2,side2,axis2,\
!*           boundaryCondition2(0,0),boundaryCondition2(1,0),boundaryCondition2(0,1),boundaryCondition2(1,1),boundaryCondition2(0,2),boundaryCondition2(1,2)
!*           ! '
!*        write(debugFile,'("n1a,n1b,...=",6i5)') n1a,n1b,n2a,n2b,n3a,n3b
!*        write(debugFile,'("m1a,m1b,...=",6i5)') m1a,m1b,m2a,m2b,m3a,m3b
!*
!*      end if
!*      if( debug.gt.4 )then
!*       write(debugFile,'("start u1=",(3i4,1x,3e11.2))') (((i1,i2,i3,(u1(i1,i2,i3,m),m=0,2),i1=nd1a,nd1b),i2=nd2a,nd2b),i3=nd3a,nd3b)
!*       write(debugFile,'("start u2=",(3i4,1x,3e11.2))') (((i1,i2,i3,(u2(i1,i2,i3,m),m=0,2),i1=md1a,md1b),i2=md2a,md2b),i3=md3a,md3b)
!*      end if
!*     

      epsx=1.e-20  ! fix this 



      dispersive=0 ! >0 implies at least one side is dispersive
      gdmParOption=1 ! scale a0 and a1 parameters by eps
      if( dispersionModel1.ne.noDispersion )then
        dispersive=dispersive+1
  
        ! get the gdm parameters
        !   gdmPar(0:3,jv) = (a0,a1,b0,b1) 
        call getGDMParameters( grid1,alphaP1,gdmPar1,numberOfPolarizationVectors1, maxNumberOfParameters,maxNumberOfPolarizationVectors,gdmParOption )

        if( t.le. 1.5*dt .and. debug.gt.0 )then
          ! ---- Dispersive Maxwell ----
          write(*,'("--interface3d4-- dispersionModel1=",i4," grid1=",i4," pxc=",i4)') dispersionModel1,grid1,pxc
          write(*,'("--interface3d4-- GDM: numberOfPolarizationVectors1=",i4," alphaP1=",e8.2)') numberOfPolarizationVectors1,alphaP1
          do jv=0,numberOfPolarizationVectors1-1
            write(*,'("--interface3d-- GDM: eqn=",i3," a0,a1,b0,b1=",4(1p,e10.2))') jv,a0v1(jv),a1v1(jv),b0v1(jv),b1v1(jv)
          end do 
        end if
      else
        numberOfPolarizationVectors1=0
      end if

      if( dispersionModel2.ne.noDispersion )then
        dispersive=dispersive+1

        ! get the gdm parameters
        !   gdmPar(0:3,jv) = (a0,a1,b0,b1) 
        call getGDMParameters( grid2,alphaP2,gdmPar2,numberOfPolarizationVectors2, maxNumberOfParameters,maxNumberOfPolarizationVectors,gdmParOption )
        if( t.le. 1.5*dt .and. debug.gt.0 )then
          ! ---- Dispersive Maxwell ----
          write(*,'("--interface3d4-- dispersionModel2=",i4," grid2=",i4)') dispersionModel2,grid2
          write(*,'("--interface3d4-- GDM: numberOfPolarizationVectors2=",i4," alphaP2=",e8.2)') numberOfPolarizationVectors2,alphaP2

          do jv=0,numberOfPolarizationVectors2-1
            write(*,'("--interface3d4-- GDM: eqn=",i3," a0,a1,b0,b1=",4(1p,e10.2))') jv,a0v2(jv),a1v2(jv),b0v2(jv),b1v2(jv)
          end do 
        end if

        ! write(*,'(" interface: FINISH ME")') 
        ! stop 1111
      else
        numberOfPolarizationVectors2=0
      end if



      initializeInterfaceVariablesMacro(interface3dOrder4)

!      do kd=0,nd-1
!       dx112(kd) = 1./(2.*dx1(kd))
!       dx122(kd) = 1./(dx1(kd)**2)
!       dx212(kd) = 1./(2.*dx2(kd))
!       dx222(kd) = 1./(dx2(kd)**2)
!
!       dx141(kd) = 1./(12.*dx1(kd))
!       dx142(kd) = 1./(12.*dx1(kd)**2)
!       dx241(kd) = 1./(12.*dx2(kd))
!       dx242(kd) = 1./(12.*dx2(kd)**2)
!
!       dr114(kd) = 1./(12.*dr1(kd))
!       dr214(kd) = 1./(12.*dr2(kd))
!      end do
!
!      numGhost=orderOfAccuracy/2
!      giveDiv=0   ! set to 1 to give div(u) on both sides, rather than setting the jump in div(u)
!
!      ! bounds for loops that include ghost points in the tangential directions:
!      nn1a=n1a
!      nn1b=n1b
!      nn2a=n2a
!      nn2b=n2b
!      nn3a=n3a
!      nn3b=n3b
!
!      mm1a=m1a
!      mm1b=m1b
!      mm2a=m2a
!      mm2b=m2b
!      mm3a=m3a
!      mm3b=m3b
!
!      i3=n3a
!      j3=m3a
!
!      axis1p1=mod(axis1+1,nd)
!      axis1p2=mod(axis1+2,nd)
!      axis2p1=mod(axis2+1,nd)
!      axis2p2=mod(axis2+2,nd)
!
!      is1=0
!      is2=0
!      is3=0
!
!      ! -----
!      ! Set bounds nn1a,nn1b,... used for loops that include ghost points in tangential directions 
!      ! -----
!      if( axis1.ne.0 )then
!        ! include ghost lines in tangential periodic (and parallel) directions (for extrapolating)
!        ! NOTE: ! parallel ghost may only have bc<0 on one side
!        if( boundaryCondition1(0,0).lt.0 .and. boundaryCondition2(0,0).ge.0 )then
!          write(*,'("Interface3d: bc is inconsistent")')
!          stop 178
!        end if
!        if( boundaryCondition1(1,0).lt.0 .and. boundaryCondition2(1,0).ge.0 )then
!          write(*,'("Interface3d: bc is inconsistent")')
!          stop 179
!        end if
!        if( boundaryCondition1(0,0).le.0 )then ! *wdh* 090506 include interp 
!          nn1a=nn1a-numGhost
!        end if
!        if( boundaryCondition1(1,0).le.0 )then ! parallel ghost may only have bc<0 on one side
!          nn1b=nn1b+numGhost
!        end if
!      end if
!      if( axis1.ne.1 )then
!        ! include ghost lines in tangential periodic (and parallel) directions (for extrapolating)
!        if( boundaryCondition1(0,1).lt.0 .and. boundaryCondition2(0,1).ge.0 )then
!          write(*,'("Interface3d: bc is inconsistent")')
!          stop 180
!        end if
!        if( boundaryCondition1(1,1).lt.0 .and. boundaryCondition2(1,1).ge.0 )then
!          write(*,'("Interface3d: bc is inconsistent")')
!          stop 181
!        end if
!        if( boundaryCondition1(0,1).le.0 )then
!          nn2a=nn2a-numGhost
!        end if
!        if( boundaryCondition1(1,1).le.0 )then
!          nn2b=nn2b+numGhost
!        end if
!      end if
!      if( nd.eq.3 .and. axis1.ne.2 )then
!        ! include ghost lines in tangential periodic (and parallel) directions (for extrapolating)
!        if( boundaryCondition1(0,2).lt.0 .and. boundaryCondition2(0,2).ge.0 )then
!          write(*,'("Interface3d: bc is inconsistent")')
!          stop 182
!        end if
!        if( boundaryCondition1(1,2).lt.0 .and. boundaryCondition2(1,2).ge.0 )then
!          write(*,'("Interface3d: bc is inconsistent")')
!          stop 183
!        end if
!        if( boundaryCondition1(0,2).le.0 )then
!          nn3a=nn3a-numGhost
!        end if
!        if( boundaryCondition1(1,2).le.0 )then
!          nn3b=nn3b+numGhost
!        end if
!      end if
!
!      if( axis1.eq.0 ) then
!        is1=1-2*side1
!        an1Cartesian=1. ! normal for a cartesian grid
!        an2Cartesian=0.
!        an3Cartesian=0.
!
!      else if( axis1.eq.1 )then
!        is2=1-2*side1
!        an1Cartesian=0.
!        an2Cartesian=1.
!        an3Cartesian=0.
!
!      else if( axis1.eq.2 )then
!        is3=1-2*side1
!        an1Cartesian=0.
!        an2Cartesian=0.
!        an3Cartesian=1.
!      else
!        stop 5528
!      end if
!
!
!      js1=0
!      js2=0
!      js3=0
!      ! -----
!      ! Set bounds mm1a,mm1b,... used for loops that include ghost points in tangential directions 
!      ! -----
!      if( axis2.ne.0 )then
!        if( boundaryCondition2(0,0).le.0 )then ! *wdh* 090506 include interp 
!          mm1a=mm1a-numGhost
!        end if
!        if( boundaryCondition2(1,0).le.0 )then
!          mm1b=mm1b+numGhost
!        end if
!      end if
!      if( axis2.ne.1 )then
!        if( boundaryCondition2(0,1).le.0 )then
!          mm2a=mm2a-numGhost
!        end if
!        if( boundaryCondition2(1,1).le.0 )then
!          mm2b=mm2b+numGhost
!        end if
!      end if
!      if( nd.eq.3 .and. axis2.ne.2 )then
!        if( boundaryCondition2(0,2).le.0 )then
!          mm3a=mm3a-numGhost
!        end if
!        if( boundaryCondition2(1,2).le.0 )then
!          mm3b=mm3b+numGhost
!        end if
!      end if
!      if( axis2.eq.0 ) then
!        js1=1-2*side2
!      else if( axis2.eq.1 ) then
!        js2=1-2*side2
!      else  if( axis2.eq.2 ) then
!        js3=1-2*side2
!      else
!        stop 3384
!      end if
!
!      is=1-2*side1
!      js=1-2*side2
!
!!$$$      rx1=0.
!!$$$      ry1=0.
!!$$$      rz1=0.
!!$$$      if( axis1.eq.0 )then
!!$$$        rx1=1.
!!$$$      else if( axis1.eq.1 )then
!!$$$        ry1=1.
!!$$$      else 
!!$$$        rz1=1.
!!$$$      endif
!!$$$
!!$$$      rx2=0.
!!$$$      ry2=0.
!!$$$      rz2=0.
!!$$$      if( axis2.eq.0 )then
!!$$$        rx2=1.
!!$$$      else if( axis2.eq.1 )then
!!$$$        ry2=1.
!!$$$      else 
!!$$$        rz2=1.
!!$$$      endif
!
!      if( debug.gt.3 )then
!        write(debugFile,'("nn1a,nn1b,...=",6i5)') nn1a,nn1b,nn2a,nn2b,nn3a,nn3b
!        write(debugFile,'("mm1a,mm1b,...=",6i5)') mm1a,mm1b,mm2a,mm2b,mm3a,mm3b
!
!      end if
!
!      if( orderOfAccuracy.eq.2 .and. orderOfExtrapolation.lt.3 )then
!        write(debugFile,'(" ERROR: interface3d: orderOfExtrapolation<3 ")')
!        stop 7716
!      end if
!      if( orderOfAccuracy.eq.4 .and. orderOfExtrapolation.lt.4 )then
!        write(debugFile,'(" ERROR: interface3d: orderOfExtrapolation<4 ")')
!        stop 7716
!      end if

      ! first time through check that the mask's are consistent
      ! For now we require the masks to both be positive at the same points on the interface
      ! We assign pts where both mask1 and mask2 are discretization pts.
      ! If mask1>0 and mask2<0 then we just leave the extrapolated values in u1 and u2 .  
      if( initialized.eq.0 )then
       if( nd.eq.2 )then
        ! check the consistency of the mask arrays
        beginLoops2d()
          m1 = mask1(i1,i2,i3)
          m2 = mask2(j1,j2,j3)
          if( (m1.gt.0 .and. m2.eq.0) .or. (m1.eq.0 .and. m2.gt.0) )then
            write(debugFile,'(" interface3d:ERROR: mask1 and mask2 do not agree. One is >0 and one =0 ")')
             ! '
            stop 1111
          end if 
        endLoops2d()

       else if( nd.eq.3 )then
        ! check the consistency of the mask arrays
        beginLoops3d()
          m1 = mask1(i1,i2,i3)
          m2 = mask2(j1,j2,j3)
          if( (m1.gt.0 .and. m2.eq.0) .or. (m1.eq.0 .and. m2.gt.0) )then
            write(debugFile,'(" interface3d:ERROR: mask1 and mask2 do not agree. One is >0 and one =0")')
             ! '
            stop 1111
          end if 
        endLoops3d()

       end if
       if( debug.gt.0 )then
         write(debugFile,'("cgmx:interface3d: The mask arrays for grid1=",i3," and grid2=",i3," were found to be consistent")') grid1,grid2
         ! ' 
       end if
      end if


      if( initialized.eq.0 )then
        ! first time through put bogus values in the following arrays to check for un-assigned values
        !! bigValue=1.e99
        bigValue=0.

        setValues9(bigValue, u1x,u1y,u1z,u1xx,u1xy,u1yy,u1xz,u1yz,u1zz)
        setValues9(bigValue, u2x,u2y,u2z,u2xx,u2xy,u2yy,u2xz,u2yz,u2zz)
        setValues9(bigValue, v1x,v1y,v1z,v1xx,v1xy,v1yy,v1xz,v1yz,v1zz)
        setValues9(bigValue, v2x,v2y,v2z,v2xx,v2xy,v2yy,v2xz,v2yz,v2zz)
        setValues9(bigValue, w1x,w1y,w1z,w1xx,w1xy,w1yy,w1xz,w1yz,w1zz)
        setValues9(bigValue, w2x,w2y,w2z,w2xx,w2xy,w2yy,w2xz,w2yz,w2zz)


        setValues9(bigValue, u1xxx,u1xxy,u1xyy,u1yyy, u1xxz,u1xzz,u1zzz, u1yyz, u1yzz)
        setValues9(bigValue, u2xxx,u2xxy,u2xyy,u2yyy, u2xxz,u2xzz,u2zzz, u2yyz, u2yzz)
        setValues9(bigValue, v1xxx,v1xxy,v1xyy,v1yyy, v1xxz,v1xzz,v1zzz, v1yyz, v1yzz)
        setValues9(bigValue, v2xxx,v2xxy,v2xyy,v2yyy, v2xxz,v2xzz,v2zzz, v2yyz, v2yzz)
        setValues9(bigValue, w1xxx,w1xxy,w1xyy,w1yyy, w1xxz,w1xzz,w1zzz, w1yyz, w1yzz)
        setValues9(bigValue, w2xxx,w2xxy,w2xyy,w2yyy, w2xxz,w2xzz,w2zzz, w2yyz, w2yzz)

        setValues6(bigValue, u1xxxx,u1xxyy,u1yyyy, u1xxzz,u1zzzz, u1yyzz)
        setValues6(bigValue, u2xxxx,u2xxyy,u2yyyy, u2xxzz,u2zzzz, u2yyzz)
        setValues6(bigValue, v1xxxx,v1xxyy,v1yyyy, v1xxzz,v1zzzz, v1yyzz)
        setValues6(bigValue, v2xxxx,v2xxyy,v2yyyy, v2xxzz,v2zzzz, v2yyzz)
        setValues6(bigValue, w1xxxx,w1xxyy,w1yyyy, w1xxzz,w1zzzz, w1yyzz)
        setValues6(bigValue, w2xxxx,w2xxyy,w2yyyy, w2xxzz,w2zzzz, w2yyzz)

        setValues6(bigValue, u1LapSq,v1LapSq,u2LapSq,v2LapSq,w1LapSq,w2LapSq)

        setValues9(bigValue, uex,uey,uez, vex,vey,vez, wex,wey,wez)
        setValues9(bigValue, uexx,ueyy,uezz, vexx,veyy,vezz, wexx,weyy,wezz)

        setValues10(bigValue, uexxx,uexxy,uexyy,ueyyy, uexxz,uexzz,uexyz, ueyyz, ueyzz, uezzz)
        setValues10(bigValue, vexxx,vexxy,vexyy,veyyy, vexxz,vexzz,vexyz, veyyz, veyzz, vezzz)
        setValues10(bigValue, wexxx,wexxy,wexyy,weyyy, wexxz,wexzz,wexyz, weyyz, weyzz, wezzz)

        setValues7(bigValue, uexxxx,uexxyy,ueyyyy,uexxzz,ueyyzz,uezzzz,ueLapSq)
        setValues7(bigValue, vexxxx,vexxyy,veyyyy,vexxzz,veyyzz,vezzzz,veLapSq)
        setValues7(bigValue, wexxxx,wexxyy,weyyyy,wexxzz,weyyzz,wezzzz,weLapSq)

        do kd1=0,nd-1
        do kd2=0,nd-1
          rx1(kd1,kd2)=bigValue
          rx2(kd1,kd2)=bigValue
          do kd3=0,nd-1
            rxx1(kd1,kd2,kd3)=bigValue
            rxx2(kd1,kd2,kd3)=bigValue
            do kd4=0,nd-1
              rxxx1(kd1,kd2,kd3,kd4)=bigValue           
              rxxx2(kd1,kd2,kd3,kd4)=bigValue           
              do kd5=0,nd-1
                rxxxx1(kd1,kd2,kd3,kd4,kd5)=bigValue 
                rxxxx2(kd1,kd2,kd3,kd4,kd5)=bigValue 
              end do
            end do
          end do
        end do
        end do
      end if

      extrapolateSecondGhost=interfaceEquationsOption.eq.0  ! apply interface conditions where 2nd ghost line is extrapolated

      if( extrapolateSecondGhost .and. nd.eq.3 .and. orderOfAccuracy.eq.4 .and. gridType.eq.curvilinear )then

        ! *************************************************************************
        ! ***** 3D fourth-order curvilinear case (extrapolate 2nd ghost line) *****
        ! *************************************************************************

        if( t.lt.3*dt .and. debug.gt.2 )then
          write(*,'("interface3d : *OLD* order 4 with extrapolation for 2nd ghost t=",e10.2)') t
        end if

        if( solveForH .ne.0 )then
          stop 3017
        end if

#perl $DIM=3; $GRIDTYPE="curvilinear"; $ORDER=4; 

        if( .false. )then
          beginGhostLoops3d()
           write(debugFile,'(" -->START v1(",i2,":",i2,",",i2,",",i2,") =",3f9.4)') i1-1,i1+1,i2,i3,u1(i1-1,i2,i3,ey),u1(i1,i2,i3,ey),u1(i1+1,i2,i3,ey)
           ! '
          endLoops3d()
        end if

         ! ---- first satisfy the jump conditions on the boundary --------
         !    [ eps n.u ] = 0
         !    [ tau.u ] = 0
         if( assignInterfaceValues.eq.1 )then
           boundaryJumpConditions(3,curvilinear)
         end if

        if( .false. )then
          beginGhostLoops3d()
           write(debugFile,'(" -->JUMP v1(",i2,":",i2,",",i2,",",i2,") =",3f9.4)') i1-1,i1+1,i2,i3,u1(i1-1,i2,i3,ey),u1(i1,i2,i3,ey),u1(i1+1,i2,i3,ey)
           ! '
          endLoops3d()
        end if

        ! ----------------------------------------------
        ! ----- assign ghost using jump conditions -----
        ! ----------------------------------------------
        if( assignInterfaceGhostValues.eq.1 )then
         ! initialization step: assign first two ghost line by extrapolation
         ! NOTE: assign ghost points outside the ends
 
         extrapGhost()
 
 
          ! here are the jump conditions for the ghost points
          !   [ div(E) n + (curl(E)- n.curl(E) n )/mu ] =0                 (3 eqns)
          !   [ Lap(E)/(eps*mu) + (1/mu)*(1-1/eps)*( n.Lap(E) ) n ] = 0    (3 eqns)
 
          ! These correspond to the 6 conditions:
          !   [ div(E) ] =0 
          !   [ tau. curl(E)/mu ] = 0       (2 tangents)
          !   [ n.Lap(E)/mu ] = 0 
          !   [ tau.Lap(E)/(eps*mu) ] = 0   (2 tangents)
 
          errOld=err
          ratioAve=0.
          err2Old=1. 
          do it=1,nit ! *** begin iteration ****
 
            assignInterfaceGhost3dOrder4Old()
 
            if( it.eq.1 )then
              errRatio=1.
              err2Ratio=1.
            else
              errRatio=err/errOld
              ratioAve=ratioAve+errRatio
              err2Ratio=err2/err2Old
            end if 
            errOld=err
            err2Old=err2
 
            if( t.le.5*dt )then
             write(*,'("interface3d[old] : t=",e10.3," (grid1,grid2)=(",i3,",",i3,"), it=",i3,", err[max,l2]=",2(e10.2,1x)," rate[max,l2]=",2(f5.2,1x)," (omega=",f4.2,")")') t,grid1,grid2,it,err,err2,errRatio,err2Ratio,omega
            end if

            numberOfIterations=it
            if( err.lt.absoluteErrorTolerance )then
              exit
            end if

          end do ! end it
 
          ipar(43) = numberOfIterations  ! returned value
          rpar(22) = ratioAve/(numberOfIterations-1)
          rpar(23) = err !  maxFinalResidual : return value 


          ! periodic update
          periodicUpdate3d(u1,boundaryCondition1,gridIndexRange1,side1,axis1)
          periodicUpdate3d(u2,boundaryCondition2,gridIndexRange2,side2,axis2)

        end if ! end assignInterfaceGhostValues
 
       else if( nd.eq.3 .and. orderOfAccuracy.eq.4 .and. gridType.eq.curvilinear )then

        ! *** NEW WAY *** use equations to get values on the 2nd ghost line

        ! ********************************************
        ! ***** 3D fourth-order curvilinear case *****
        ! ********************************************

        if( t.lt.3*dt .and. debug.gt.2 )then
          write(*,'("interface3d : *NEW* order 4 with centered approx for 2nd ghost t=",e10.2)') t
        end if

        if( solveForH .ne.0 )then
          stop 3017
        end if

#perl $DIM=3; $GRIDTYPE="curvilinear"; $ORDER=4; 

        if( .false. )then
          beginGhostLoops3d()
           write(debugFile,'(" -->START v1(",i2,":",i2,",",i2,",",i2,") =",3f9.4)') i1-1,i1+1,i2,i3,u1(i1-1,i2,i3,ey),u1(i1,i2,i3,ey),u1(i1+1,i2,i3,ey)
           ! '
          endLoops3d()
        end if

         ! ---- first satisfy the jump conditions on the boundary --------
         !    [ eps n.u ] = 0
         !    [ tau.u ] = 0
         if( assignInterfaceValues.eq.1 )then
           boundaryJumpConditions(3,curvilinear)
         end if

        ! ----------------------------------------------
        ! ----- assign ghost using jump conditions -----
        ! ----------------------------------------------
        if( assignInterfaceGhostValues.eq.1 )then

         if( .false. )then
           beginGhostLoops3d()
            write(debugFile,'(" -->JUMP v1(",i2,":",i2,",",i2,",",i2,") =",3f9.4)') i1-1,i1+1,i2,i3,u1(i1-1,i2,i3,ey),u1(i1,i2,i3,ey),u1(i1+1,i2,i3,ey)
            ! '
           endLoops3d()
         end if
 
         ! initialization step: assign first two ghost line by extrapolation
         ! NOTE: assign ghost points outside the ends
 
         if( interfaceOption.eq.1 )then
           extrapGhost()
         end if
 
 
         ! perl $DIM=3; $GRIDTYPE="curvilinear"; $ORDER=2; 

         if( avoidInterfaceIterations.eq.1 )then
           ! Note: checkCoeff fails below for 2nd-order coeff's (they are correct)
           ! Check coeff is OK if we set the next line, but then accuracy degrades to 2 !! why?
           !! perl $DIM=3; $GRIDTYPE="curvilinear"; $ORDER=2;

           ! **TEST** assign ghost to 2nd-order first
           if( t .lt. 3.*dt .and. debug.gt.0 )then
             write(*,'(" **** STAGE I: ASSIGN GHOST TO 2ND ORDER ****")')
           end if

           orderOfAccuracy=2 ! temporarily set (for checkCoeff only?)

           ! in parallel we add extra points in the tangential direction on parallel boundaries
           ! (otherwise we would use extrapolated values which is probably ok) 
           setIndexBoundsExtraGhost()

           if( dispersive.eq.0 )then
             assignInterfaceGhost23c()
           else
             assignDispersiveInterfaceGhost23c()
           end if
           resetIndexBounds()         


           orderOfAccuracy=4 ! reset 
           !! perl $DIM=3; $GRIDTYPE="curvilinear"; $ORDER=4; 
         end if
         ! perl $DIM=3; $GRIDTYPE="curvilinear"; $ORDER=4; 

         ! here are the jump conditions for the ghost points
         !   [ div(E) n + (curl(E)- n.curl(E) n )/mu ] =0                 (3 eqns)
         !   [ Lap(E)/(eps*mu) + (1/mu)*(1-1/eps)*( n.Lap(E) ) n ] = 0    (3 eqns)
 
         ! These correspond to the 6 conditions:
         !   [ div(E) ] =0 
         !   [ tau. curl(E)/mu ] = 0       (2 tangents)
         !   [ n.Lap(E)/mu ] = 0 
         !   [ tau.Lap(E)/(eps*mu) ] = 0   (2 tangents)
 
         ! Here are the equations for the 2nd ghost line (approximate to 2nd order):
         !   [ div(Lap(E)) n/(eps*mu) + (curl(Lap(E))- n.curl(Lap(E)) n )/(mu*eps*mu) ] =0   (3 eqns)
         !   [ Lap^2(E)/(eps*mu)^2 + (1/(eps*mu)^2)*(eps-1)*( n.Lap^2(E) ) n ] = 0           (3 eqns)
          
          errOld=1.
          err2Old=1.
          ratioAve=0.
          do it=1,nit ! *** begin iteration ****
            err=0.
            err2=0.

            if( dispersive.eq.0 )then
              assignInterfaceGhost3dOrder4()
            else
              assignDispersiveInterfaceGhost3dOrder4()
            end if            
 
            if( it.eq.1 )then
              errRatio=1.
              err2Ratio=1.
            else
              errRatio=err/errOld
              ratioAve=ratioAve+errRatio
              err2Ratio=err2/err2Old
            end if 
            errOld=err
            err2Old=err2

            if( t.le.5*dt .and. debug.gt.0 )then
             write(*,'("interface3d : t=",e10.3," (grid1,grid2)=(",i3,",",i3,"), it=",i3,", err[max,l2]=",2(e10.2,1x)," rate[max,l2]=",2(f5.2,1x)," (omega=",f4.2,")")') t,grid1,grid2,it,err,err2,errRatio,err2Ratio,omega
            end if

            if( debug.gt.0 )then 
             write(*,'("interface3d : t=",e10.3," (grid1,grid2)=(",i3,",",i3,"), it=",i3,", err[max,l2]=",2(e10.2,1x)," rate[max,l2]=",2(f5.2,1x)," (omega=",f4.2,")")') t,grid1,grid2,it,err,err2,errRatio,err2Ratio,omega

            end if
            if( it.eq.1 )then
               checkCoeff=0 ! turn off after 1 iteration
            end if
            numberOfIterations=it
            if( err.lt.absoluteErrorTolerance )then
              exit
            end if

          end do ! end it
 
          ipar(43) = numberOfIterations  ! returned value
          rpar(22) = ratioAve/(numberOfIterations-1)
          rpar(23) = err !  maxFinalResidual : return value 
 
          ! periodic update
          periodicUpdate3d(u1,boundaryCondition1,gridIndexRange1,side1,axis1)
          periodicUpdate3d(u2,boundaryCondition2,gridIndexRange2,side2,axis2)
        end if ! end assignInterfaceGhostValues

       else
         write(debugFile,'("interface3d: ERROR: unknown options nd,order=",2i3)') nd,orderOfAccuracy
         stop 3214
       end if

      return
      end
