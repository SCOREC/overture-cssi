! This file automatically generated from advOptIsm.bf90 with bpp.
!
! Advance the equations of solid mechanics
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


! ---------------------------------------------------------------------------
! Macro : beginLoops
! ---------------------------------------------------------------------------

! ---------------------------------------------------------------------------
! Macro : endLoops
! ---------------------------------------------------------------------------

! ---------------------------------------------------------------------------
! Macro : beginLoopsMask
! ---------------------------------------------------------------------------

! ---------------------------------------------------------------------------
! Macro : endLoopsMask
! ---------------------------------------------------------------------------


! ************************************************************************************************
!  This macro is used for looping over the faces of a grid to assign boundary conditions
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
! #beginMacro beginLoopOverSides(extra,numberOfGhostPoints)




! ===================================================================================================
! Macro: advance the solution
!   FORCING : NOFORCING or TWILIGHTZONEFORCING
! ===================================================================================================

! ===================================================================================================
! Macro: advance the solution SOS -- CORRECTOR 
!   FORCING : NOFORCING or TWILIGHTZONEFORCING
! ===================================================================================================


! ===================================================================================================
! Macro: Evaluate the time-derivative of v:
!
!         v.t = -(1/rho) grad(p) + cs * Delta( u) 
!
!   FORCING : NOFORCING or TWILIGHTZONEFORCING
! ===================================================================================================

! ===================================================================================================
! Macro: Method-of-lines update - TWO TERMS 
!
! ===================================================================================================

! ===================================================================================================
! Macro: Method-of-lines update - THREE TERMS 
!
! ===================================================================================================

! ===================================================================================================
! Macro: Method-of-lines update - FOUR TERMS 
!
! ===================================================================================================


! ===========================================================================================
! Macro: compute the coefficients in the sosup dissipation for curvilinear grids
! ===========================================================================================


! ===========================================================================================
! Macro: Output some debug info for the first few time-steps 
! ===========================================================================================


! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Upwind (sosup) dissipation 2nd-order -- used at boundaries of 2nd-order scheme
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Upwind (sosup) dissipation (4th-order difference used with 2nd-order scheme) 
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Upwind (sosup) dissipation (6th-order difference used with 4th-order scheme)
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



! =========================================================================================
! Macro: ADD UPWIND DISSIPATION
!
!   DIM : dimension : 2 or 3
!   ORDER : 2 or 4   
!   GRIDTYPE : rectangular or curvilinear
! ========================================================================================

! =========================================================================================
! Macro: COMPUTE uDot = (un -um )
!
! ========================================================================================

              

! ====== Here are the common subroutine arguments ====

! ==========================================================================================
! Macro: declare input arguements 
! ==========================================================================================


! **********************************************************************************
! NAME: name of the subroutine
! DIM : 2 or 3
! ORDER : 2 ,4, 6 or 8
! GRIDTYPE : rectangular, curvilinear
! **********************************************************************************




!**
!**      buildFile(advIsm22Order6r,2,6,rectangular)
!**      buildFile(advIsm23Order6r,3,6,rectangular)
!**
!**       ! build these for testing symmetric operators -- BC's not implemented yet
!**      buildFile(advIsm22Order6c,2,6,curvilinear)
!**      buildFile(advIsm23Order6c,3,6,curvilinear)
!**
!**      buildFile(advIsm22Order8r,2,8,rectangular)
!**      buildFile(advIsm23Order8r,3,8,rectangular)
!**
!**       ! build these for testing symmetric operators -- BC's not implemented yet
!**      buildFile(advIsm22Order8c,2,8,curvilinear)
!**      buildFile(advIsm23Order8c,3,8,curvilinear)




      subroutine advIsm( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,rsxy,xy,um,u,un,f,vt,vn,v,v1,v2,v3,v4,vt1,vt2,vt3,vt4,gridIndexRange,boundaryCondition,ipar,rpar,ierr )
!======================================================================
!   Advance routines for
!        Incompressible elasticity 
!
! nd : number of space dimensions
!
!======================================================================
  implicit none
   integer nd, n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b
   real um(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
   real u(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
   real un(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
   real f(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
   real vt(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
   real vn(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
   real  v(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
   real v1(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
   real v2(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
   real v3(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
   real v4(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
   real vt1(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
   real vt2(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
   real vt3(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
   real vt4(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
   real rsxy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1,0:nd-1)
   real xy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1)
   integer mask(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
   integer boundaryCondition(0:1,0:2),gridIndexRange(0:1,0:2),ierr
   integer ipar(0:*)
   real rpar(0:*)

  
  ! -- Declare arrays for variable material properties --
  ! include '../declareVarMatProp.h'

  !     ---- local variables -----
  real dt,dtOld
  integer c,i1,i2,i3,n,gridType,orderOfAccuracy,orderOfAccuracyInTime,useConservative
  integer addForcing,orderOfDissipation,option
  integer useWhereMask,solveForE,solveForH,grid
  integer rectangular,curvilinear
  parameter( rectangular=0, curvilinear=1 )
  !...........end   statement functions


  ! write(*,*) 'Inside advIsm...'
  dt    =rpar(0)
  dtOld =rpar(15) ! dt used on the previous time step 

  gridType           =ipar(1)
  orderOfAccuracy    =ipar(2)
  useConservative    =ipar(12)

  ! write(*,'(" advOpt: gridType=",i2," useConservative=",i2)') gridType,useConservative
  if( abs(dt-dtOld).gt.dt*.001 .and. orderOfAccuracy.ne.2 )then
   write(*,'(" advIsm:ERROR: variable dt not implemented yet for this case")')
   write(*,'("            : dt,dtOld,diff=",3e9.3)') dt,dtOld,dt-dtOld
   write(*,'("              orderOfAccuracy=",i4," useConservative=",i4)') orderOfAccuracy,useConservative
   ! '
   stop 9027
  end if


  if( orderOfAccuracy.eq.2 )then

    if( nd.eq.2 .and. gridType.eq.rectangular ) then
      call advIsm2dOrder2r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,rsxy,xy,um,u,un,f,vt,vn,v,v1,v2,v3,v4,vt1,vt2,vt3,vt4,gridIndexRange,boundaryCondition,ipar,rpar,ierr )
    else if( nd.eq.2 .and. gridType.eq.curvilinear ) then
      call advIsm2dOrder2c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,rsxy,xy,um,u,un,f,vt,vn,v,v1,v2,v3,v4,vt1,vt2,vt3,vt4,gridIndexRange,boundaryCondition,ipar,rpar,ierr )
    else if( nd.eq.3 .and. gridType.eq.rectangular ) then
      call advIsm3dOrder2r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,rsxy,xy,um,u,un,f,vt,vn,v,v1,v2,v3,v4,vt1,vt2,vt3,vt4,gridIndexRange,boundaryCondition,ipar,rpar,ierr )
    else if( nd.eq.3 .and. gridType.eq.curvilinear ) then
      call advIsm3dOrder2c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,rsxy,xy,um,u,un,f,vt,vn,v,v1,v2,v3,v4,vt1,vt2,vt3,vt4,gridIndexRange,boundaryCondition,ipar,rpar,ierr )
    else
      write(*,*) 'advOptIsm: unexpected nd,gridType=',nd,gridType
      stop 2271
    end if

  else if( orderOfAccuracy.eq.4 ) then

    if( nd.eq.2 .and. gridType.eq.rectangular ) then
      call advIsm2dOrder4r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,rsxy,xy,um,u,un,f,vt,vn,v,v1,v2,v3,v4,vt1,vt2,vt3,vt4,gridIndexRange,boundaryCondition,ipar,rpar,ierr )
    else if( nd.eq.2 .and. gridType.eq.curvilinear ) then
      call advIsm2dOrder4c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,rsxy,xy,um,u,un,f,vt,vn,v,v1,v2,v3,v4,vt1,vt2,vt3,vt4,gridIndexRange,boundaryCondition,ipar,rpar,ierr )
    else if( nd.eq.3 .and. gridType.eq.rectangular ) then
      call advIsm3dOrder4r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,rsxy,xy,um,u,un,f,vt,vn,v,v1,v2,v3,v4,vt1,vt2,vt3,vt4,gridIndexRange,boundaryCondition,ipar,rpar,ierr )
    else if( nd.eq.3 .and. gridType.eq.curvilinear ) then
      call advIsm3dOrder4c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,rsxy,xy,um,u,un,f,vt,vn,v,v1,v2,v3,v4,vt1,vt2,vt3,vt4,gridIndexRange,boundaryCondition,ipar,rpar,ierr )
    else
      stop 4444
    end if

  else
    write(*,'(" advIsm:ERROR: un-implemented order of accuracy =",i6)') orderOfAccuracy
      ! '
    stop 11222
  end if

  return
  end








