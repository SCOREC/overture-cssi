! This file automatically generated from champResidualOpt.bf90 with bpp.
!
! Optimized routines for evaluating the residual in the CHAMP conditions
!
! These next include files will define the macros that will define the difference approximations
! The actual macro is called below
! Use this next macro to declare the statement functions that are defined below
! To include derivatives of rx use OPTION=RX


! Define statement functions for difference approximations of order 2 
! To include derivatives of rx use OPTION=RX
! To include derivatives of rx use OPTION=RX



! Use this next macro to declare the statement functions that are defined below
! To include derivatives of rx use OPTION=RX


! Define statement functions for difference approximations of order 4 
! To include derivatives of rx use OPTION=RX
! To include derivatives of rx use OPTION=RX


! ======================================================================================
!   Evaluate the TZ exact solution in 2D
! ======================================================================================

! ======================================================================================
!   Evaluate the TZ exact solution in 3D
! ======================================================================================

  

! ---------------------------------------------------------------------------
! Macro : beginLoopsMask
! ---------------------------------------------------------------------------

! ---------------------------------------------------------------------------
! Macro : endLoopsMask
! ---------------------------------------------------------------------------

! ---------------------------------------------------------------------------
! Macro : beginLoops
! ---------------------------------------------------------------------------

! ---------------------------------------------------------------------------
! Macro : endLoopsMask
! ---------------------------------------------------------------------------








! ************************************************************************************************
!  This macro is used for looping over the faces of a grid to assign booundary conditions
!
! extra: extra points to assign
!          Case 1: extra=numberOfGhostPoints -- for assigning extended boundaries
!          Case 2: extra=-1 -- for assigning ghost points but not including extended boundaries
! numberOfGhostPoints : number of ghost points (1 for 2nd order, 2 for fourth-order ...)
!
!
! Output:
!  n1a,n1b,n2a,n2b,n3a,n3b : from gridIndexRange
!  nn1a,nn1b,nn2a,nn2b,nn3a,nn3b : includes "extra" points
! 
! ***********************************************************************************************





! =========================================================================
! Compute the normal on a curvilinear grid.
!
! Assumes is=1-2*side is defined. 
! =========================================================================



! ===========================================================================================
! Macro: Output some debug info for the first few time-steps 
! ===========================================================================================

  
! =========================================================================================
! Macro: Evaluate the residual (or RHS) in the CHAMP interface conditions
!
!   OPTION: RESIDUAL or RHS
!   DIM : dimension : 2 or 3
!   ORDER : 2 or 4   
!   GRIDTYPE : rectangular or curvilinear
! ========================================================================================

! Arguments to champ subroutine

! **********************************************************************************
! Macro EVAL_CHAMP 
! 
!            **** NOT USED ANYMORE*****
!
!  NAME: name of the subroutine
!  DIM : 2 or 3
!  ORDER : 2 ,4, 6 or 8
!  GRIDTYPE : rectangular, curvilinear
! **********************************************************************************


 

! THESE ARE NO LONGER USED
!     ! NOTE: For now 3D versions are just null versions below 

!     buildFile(champResidual2dOrder2r,2,2,rectangular)
!     buildFile(champResidual3dOrder2r,3,2,rectangular)
! !
!     buildFile(champResidual2dOrder2c,2,2,curvilinear)
!     buildFile(champResidual3dOrder2c,3,2,curvilinear)
! !
!*      buildFile(champ2dOrder4r,2,4,rectangular)
!--      buildFile(champ3dOrder4r,3,4,rectangular)
!
!*      buildFile(champ2dOrder4c,2,4,curvilinear)
!--      buildFile(champ3dOrder4c,3,4,curvilinear)
!
!      buildFile(advMx2dOrder6r,2,6,rectangular)
!      buildFile(advMx3dOrder6r,3,6,rectangular)
!
!       ! build these for testing symmetric operators -- BC's not implemented yet
!      buildFile(advMx2dOrder6c,2,6,curvilinear)
!      buildFile(advMx3dOrder6c,3,6,curvilinear)
!
!      buildFile(advMx2dOrder8r,2,8,rectangular)
!      buildFile(advMx3dOrder8r,3,8,rectangular)
!
!       ! build these for testing symmetric operators -- BC's not implemented yet
!      buildFile(advMx2dOrder8c,2,8,curvilinear)
!      buildFile(advMx3dOrder8c,3,8,curvilinear)



subroutine champResidualOpt( nd,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,md1a,md1b,md2a,md2b,md3a,md3b,md4a,md4b,ld1a,ld1b,ld2a,ld2b,ld3a,ld3b,ld4a,ld4b,gridIndexRange,dimRange,isPeriodic,boundaryCondition,mask,xy,rsxy,u,f,coeff,ipar,rpar,ierr )
!======================================================================
!   Evaluate the RESIDUAL in the CHAMP interface conditions.
!     MAIN ROUTINE THAT CALLS THE APPROPRIATE SPECIALIZED VERSION
!======================================================================
  implicit none
  integer nd, nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,ierr
  integer md1a,md1b,md2a,md2b,md3a,md3b,md4a,md4b
  integer ld1a,ld1b,ld2a,ld2b,ld3a,ld3b,ld4a,ld4b


  real u(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
  real f(md1a:md1b,md2a:md2b,md3a:md3b,md4a:md4b)
  real coeff(ld1a:ld1b,ld2a:ld2b,ld3a:ld3b,ld4a:ld4b)

  real xy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1)
  real rsxy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1,0:nd-1)

  integer mask(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
  integer gridIndexRange(0:1,0:2),boundaryCondition(0:1,0:2), dimRange(0:1,0:2), isPeriodic(0:*)

  integer ipar(0:*)
  real rpar(0:*)
  
  !     ---- local variables -----

  integer option,side,axis,grid,tc,debug,myid,orderOfAccuracy,gridType
  integer n1a,n1b,n2a,n2b,n3a,n3b
  integer i1,i2,i3,m1,m2,m3,m
  real res,maxRes,l2Res,count

  integer rectangular,curvilinear
  parameter( rectangular=0, curvilinear=1 )

  write(*,*) 'Inside champResidualOpt... '

   option             = ipar( 0)
   side               = ipar( 1)
   axis               = ipar( 2)
   grid               = ipar( 3)
   gridType           = ipar( 4)
   orderOfAccuracy    = ipar( 5)
   ! twilightZone       = ipar( 6)
   tc                 = ipar( 7)
   debug              = ipar( 8)
   myid               = ipar( 9)

  maxRes   = rpar(13)  ! ... returned
  l2Res    = rpar(14)  ! ... returned

  ! --- get loop bounds -----
  n1a=gridIndexRange(0,0)
  n1b=gridIndexRange(1,0)
  n2a=gridIndexRange(0,1)
  n2b=gridIndexRange(1,1)
  n3a=gridIndexRange(0,2)
  n3b=gridIndexRange(1,2)
  if( axis.eq.0 )then
    n1a=gridIndexRange(side,axis)
    n1b=gridIndexRange(side,axis)
  else if( axis.eq.1 )then
    n2a=gridIndexRange(side,axis)
    n2b=gridIndexRange(side,axis)
  else
    n3a=gridIndexRange(side,axis)
    n3b=gridIndexRange(side,axis)
  end if
   
  ! ---------------- START LOOP OVER FACE -----------------------
  ! if( .true. )then
  ! -- compute the residual using the stencil in the coefficient matrix ---
  maxRes=0.
  l2Res=0.
  count=0.
  if( nd.eq.2 .and. orderOfAccuracy.eq.2 )then
    !  ------------ 2D ORDER=2 ----------------

      do i3=n3a,n3b
      do i2=n2a,n2b
      do i1=n1a,n1b
        if( mask(i1,i2,i3).gt.0 )then

      res = -f(i1,i2,i3,0)
      do m2=-1,1
        do m1=-1,1
          m = (m1+1)+3*(m2+1)
          res = res + coeff(m,i1,i2,i3)*u(i1+m1,i2+m2,i3,tc)
        end do 
      end do

      res = abs(res)
      maxRes = max(maxRes,res)
      l2res = l2res + res**2
      count = count +1 

        end if
      end do
      end do
      end do
  else
    stop 9876
  end if

  l2res = sqrt( l2res )/max(1.0,count)
  rpar(13) = maxRes  ! ... returned
  rpar(14) = l2Res   ! ... returned

  ! else
  !   ! --- OlD WAY ----

  !   if( orderOfAccuracy.eq.2 )then

  !     if( nd.eq.2 .and. gridType.eq.rectangular ) then
  !       call champResidual2dOrder2r( MYARGS() )
  !     else if( nd.eq.2 .and. gridType.eq.curvilinear ) then
  !       call champResidual2dOrder2c( MYARGS() )
  !     else if( nd.eq.3 .and. gridType.eq.rectangular ) then
  !       call champResidual3dOrder2r( MYARGS() )
  !     else if( nd.eq.3 .and. gridType.eq.curvilinear ) then
  !       call champResidual3dOrder2c( MYARGS() )
  !     else
  !       stop 2271
  !     end if

  !   else if( orderOfAccuracy.eq.4 ) then
  !     write(*,'(" champResidualOpt:ERROR: un-implemented order of accuracy =",i6)') orderOfAccuracy
  !     stop 4444

  !   else
  !     write(*,'(" champResidualOpt:ERROR: un-implemented order of accuracy =",i6)') orderOfAccuracy
  !       ! '
  !     stop 11122
  !   end if

  ! end if ! old way 

  return
end








