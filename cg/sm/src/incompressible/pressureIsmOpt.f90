! This file automatically generated from pressureIsmOpt.bf90 with bpp.
!
! FILL IN THE PRESSURE RHS AND BOUNDARY CONDITIONS FOR INCOMPRESSIBLE SOLID MECHANICS
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


! ogf2d, ogf3d, ogDeriv2, etc. are foundin forcing.bC


! ntd,nxd,nyd,nzd : number of derivatives to evaluate in t,x,y,z



! ----- define extrapolation formulae ------







! *************************************************************************************
!  Macro : adjust the gridIndexRange  to account for adjacent sides
!    Precedence at corners:
!       displacement over traction
!       dirichlet over traction
!  
!                  |   |   |
!                  |------------
!  bc=traction     |   |   |
!    =displacement |-----------
!                  |   |   |
!                  X------------
!                     bc=dirichletBoundaryCondition
!     X = this point removed from boundary loops for bc=traction, or bc=displacement
! **************************************************************************************

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





! =============================================================
!  Macro: Fill in the RHS to to the press CBC[2]:
!     p_xxxx - p_yyyy = 0 
! =============================================================

! ====== here are the common arguments to the various subroutines =====

! ==============================================================================================
! Macro to declare input arguments to subroutines
! ==============================================================================================

! **********************************************************************************
! NAME: name of the subroutine
! DIM : 2 or 3
! ORDER : 2 ,4, 6 or 8
! GRIDTYPE : rectangular, curvilinear
! **********************************************************************************


 


      ! buildFile(pressureIsm3dOrder4r,3,4,rectangular)
      ! buildFile(pressureIsm3dOrder4c,3,4,curvilinear)
!**
!**      buildFile(pressureIsm22Order6r,2,6,rectangular)
!**      buildFile(pressureIsm23Order6r,3,6,rectangular)
!**
!**       ! build these for testing symmetric operators -- BC's not implemented yet
!**      buildFile(pressureIsm22Order6c,2,6,curvilinear)
!**      buildFile(pressureIsm23Order6c,3,6,curvilinear)
!**
!**      buildFile(pressureIsm22Order8r,2,8,rectangular)
!**      buildFile(pressureIsm23Order8r,3,8,rectangular)
!**
!**       ! build these for testing symmetric operators -- BC's not implemented yet
!**      buildFile(pressureIsm22Order8c,2,8,curvilinear)
!**      buildFile(pressureIsm23Order8c,3,8,curvilinear)






      subroutine pressureIsmOPt( nd,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,rsxy,xy,um,u,un,f,gridIndexRange,dimRange,isPeriodic,boundaryCondition,ipar,rpar,ierr )
!=====================================================================================
!      PRESSURE RHS
!   This function calls the appropriate lower level routine for a 
!  given dimension, gridType and order of accuracy
!
!====================================================================================
      implicit none

       integer nd, nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b
       real um(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
       real u(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
       real un(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
       real f(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
       real rsxy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1,0:nd-1)
       real xy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1)
       integer mask(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
       integer gridIndexRange(0:1,0:2),boundaryCondition(0:1,0:2), dimRange(0:1,0:2), isPeriodic(0:*), ierr
       integer ipar(0:*)
       real rpar(0:*)

      
      ! -- Declare arrays for variable material properties --
      ! include '../declareVarMatProp.h'

      !     ---- local variables -----
      integer option,gridType,orderOfAccuracy,i

      integer rectangular,curvilinear
      parameter( rectangular=0, curvilinear=1 )
      !...........end   statement functions


      ! write(*,*) 'Inside pressureIsmOpt...'
      ! write(*,*) 'rpar=',(rpar(i),i=0,10)
      ! write(*,*) 'ipar=',(ipar(i),i=0,10)
 
      option                        = ipar( 0)
      gridType                      = ipar( 1)
      orderOfAccuracy               = ipar( 2)      
 

      if( orderOfAccuracy.eq.2 )then

        if( nd.eq.2 .and. gridType.eq.rectangular ) then
          call pressureIsm2dOrder2r( nd,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,rsxy,xy,um,u,un,f,gridIndexRange,dimRange,isPeriodic,boundaryCondition,ipar,rpar,ierr )
        else if( nd.eq.2 .and. gridType.eq.curvilinear ) then
          call pressureIsm2dOrder2c(  nd,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,rsxy,xy,um,u,un,f,gridIndexRange,dimRange,isPeriodic,boundaryCondition,ipar,rpar,ierr ) 

        else if( nd.eq.3 .and. gridType.eq.rectangular ) then
          call pressureIsm3dOrder2r( nd,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,rsxy,xy,um,u,un,f,gridIndexRange,dimRange,isPeriodic,boundaryCondition,ipar,rpar,ierr )

        else if( nd.eq.3 .and. gridType.eq.curvilinear ) then
          call pressureIsm3dOrder2c(  nd,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,rsxy,xy,um,u,un,f,gridIndexRange,dimRange,isPeriodic,boundaryCondition,ipar,rpar,ierr ) 

        ! else if( nd.eq.3 .and. gridType.eq.rectangular ) then
        !   call pressureIsm3dOrder2r(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,!                      mask,rx,xy, um,u,un,f, bc, ipar, rpar, ierr )
        ! else if( nd.eq.3 .and. gridType.eq.curvilinear ) then
        !   call pressureIsm3dOrder2c(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,!                      mask,rx,xy, um,u,un,f, bc, ipar, rpar, ierr )

        else

          stop 2271

        end if

      else if( orderOfAccuracy.eq.4 ) then

        if( nd.eq.2 .and. gridType.eq.rectangular ) then
          call pressureIsm2dOrder4r( nd,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,rsxy,xy,um,u,un,f,gridIndexRange,dimRange,isPeriodic,boundaryCondition,ipar,rpar,ierr )
        else if( nd.eq.2 .and. gridType.eq.curvilinear ) then
          call pressureIsm2dOrder4c(  nd,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,rsxy,xy,um,u,un,f,gridIndexRange,dimRange,isPeriodic,boundaryCondition,ipar,rpar,ierr ) 
        ! else if( nd.eq.3 .and. gridType.eq.rectangular ) then
        !   call pressureIsm3dOrder2r(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,!                      mask,rx,xy, um,u,un,f, bc, ipar, rpar, ierr )
        ! else if( nd.eq.3 .and. gridType.eq.curvilinear ) then
        !   call pressureIsm3dOrder2c(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,!                      mask,rx,xy, um,u,un,f, bc, ipar, rpar, ierr )

        else

          stop 4444

        end if        

      else
        write(*,'(" pressureIsmOpt:ERROR: un-implemented order of accuracy =",i6)') orderOfAccuracy
          ! '
        stop 11222
      end if

      return
      end








